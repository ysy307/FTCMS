#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscviewerhdf5.h>

!> Writing what `Output.json` asks for.
!>
!> Fields go through PETSc's own viewers: VTU (`vtu`, `vtk`, `petsc_vtu`) or
!> HDF5 in PETSc's visualisation layout (`hdf5`; `petsc_gen_xdmf.py` writes
!> its XDMF), one file per output time under
!> `Output/Files/`, with a `.pvd` collection for the VTU series. A potential in
!> a `lagrange` space is written as point data; everything else is cell data,
!> evaluated at the cell's reference centroid. Time series at the observation
!> points go to `Output/obsf_*.dat` in the layout the post-processing reads.
submodule(app_ftcms) ftcms_output
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
    use :: stdlib_logger, only:global_logger
    use :: module_core, only:FE_SPACES, COEFFICIENT_EVALUATIONS, COMP_TYPES, TIME_UNITS, FILE_FORMATS, &
                              OUTPUT_TYPES, type_config_output_manager, type_constant_value
    use :: module_input, only:input_translator
    use :: core_interop_petscfe, only:petscfe_evaluate_cell, petscfe_evaluate_cell_gradient, &
                                      petscfe_reference_coordinates
    use :: app_ftcms_darcy, only:PLANE_DIMENSION
    use :: app_ftcms_closure, only:type_closure_response, type_closure_tangents, eval_closure_response, &
                                   eval_closure_tangents, closure_saturated_conductivity
    use :: app_ftcms_kernels, only:AUX_MATERIAL, AUX_INVERSE_MOBILITY
    use :: app_ftcms_kernels_lagrange, only:compute_lagrange_darcy_flux
    use :: app_ftcms_nodal_flux, only:compute_nodal_cell_flux
    implicit none

    !> What one sample of the state at a point carries, in this order.
    integer(int32), parameter :: S_TEMPERATURE = 1
    integer(int32), parameter :: S_PRESSURE = 2
    integer(int32), parameter :: S_WATER = 3
    integer(int32), parameter :: S_ICE = 4
    integer(int32), parameter :: S_VAPOR = 5
    integer(int32), parameter :: S_THERMAL_CONDUCTIVITY = 6
    integer(int32), parameter :: S_HEAT_CAPACITY = 7
    integer(int32), parameter :: S_HYDRAULIC_CONDUCTIVITY = 8
    integer(int32), parameter :: S_FLUX_X = 9
    integer(int32), parameter :: S_FLUX_Z = 10
    integer(int32), parameter :: NUM_SAMPLES = 10

contains

    module subroutine initialize_output(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        type(type_config_output_manager) :: config
        character(len=4096) :: project
        integer(int32) :: length, status, rank, mpi_ierr

        call get_environment_variable("FTCMS_PROJECT_PATH", project, length, status)
        if (status /= 0) error stop "FTCMS_PROJECT_PATH is not set"
        if (project(length:length) /= "/") then
            project(length + 1:length + 1) = "/"
            length = length + 1
        end if
        self%output_directory = project(:length)//"Output/"
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, mpi_ierr)
        if (rank == 0) call execute_command_line("mkdir -p '"//self%output_directory//"Files'", exitstat=status)
        call MPI_Barrier(PETSC_COMM_WORLD, mpi_ierr)

        call input_translator%execute(self%input, OUTPUT_TYPES%FIELD, config)
        if (.not. (config%file_format == FILE_FORMATS%NONE) .and. config%interval_val > 0.0d0) then
            self%field_interval = config%interval_val*config%interval_unit%value
            self%field_time_unit = config%output_unit%value
            self%field_format = merge("hdf5", "vtu ", config%file_format == FILE_FORMATS%HDF5)
            self%field_variables = self%input%output_settings%field_output%variable_names
            call build_output_dm(self)
        end if
        self%next_field_time = self%time

        call input_translator%execute(self%input, OUTPUT_TYPES%HISTORY, config)
        if (.not. (config%file_format == FILE_FORMATS%NONE) .and. config%interval_val > 0.0d0) then
            self%history_interval = config%interval_val*config%interval_unit%value
            self%history_time_unit = config%output_unit%value
            call open_history(self)
        end if
        self%next_history_time = self%time

        associate (standard => self%input%output_settings%standard_output)
            if (standard%print_progress) then
                self%progress_interval = standard%print_progress_interval &
                                         *unit_seconds(standard%print_progress_unit)
            end if
        end associate
        self%next_progress_time = self%time + self%progress_interval
    end subroutine initialize_output

    !> Seconds per unit of a `TIME_UNITS` name.
    function unit_seconds(name) result(seconds)
        implicit none
        character(len=*), intent(in) :: name
        real(real64) :: seconds

        type(type_constant_value) :: unit

        unit = TIME_UNITS%to_object(name)
        seconds = unit%value
    end function unit_seconds

    module subroutine write_output(self, force)
        implicit none
        class(type_ftcms), intent(inout) :: self
        logical, intent(in) :: force

        real(real64), parameter :: TOLERANCE = 1.0d-9

        if (self%field_interval > 0.0d0) then
            if (force .or. self%time >= self%next_field_time - TOLERANCE*self%field_interval) then
                if (self%field_count == 0 .or. abs(self%time - self%last_field_time) > TOLERANCE*self%field_interval) then
                    call write_fields(self)
                    self%last_field_time = self%time
                end if
                do while (self%next_field_time <= self%time + TOLERANCE*self%field_interval)
                    self%next_field_time = self%next_field_time + self%field_interval
                end do
            end if
        end if
        if (self%history_interval > 0.0d0) then
            if (force .or. self%time >= self%next_history_time - TOLERANCE*self%history_interval) then
                if (.not. self%history_written .or. &
                    abs(self%time - self%last_history_time) > TOLERANCE*self%history_interval) then
                    call write_history(self)
                    self%last_history_time = self%time
                    self%history_written = .true.
                end if
                do while (self%next_history_time <= self%time + TOLERANCE*self%history_interval)
                    self%next_history_time = self%next_history_time + self%history_interval
                end do
            end if
        end if
    end subroutine write_output

    !> The output DM: one field per requested variable, P1 or P0.
    subroutine build_output_dm(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscFE :: fe
        PetscBool :: is_simplex
        PetscErrorCode :: ierr
        integer(int32) :: v, components, order, count

        call DMClone(self%dm, self%output_dm, ierr)
        if (ierr == 0) call DMPlexIsSimplex(self%dm, is_simplex, ierr)
        if (ierr /= 0) error stop "PETSc could not clone the mesh for output"
        count = 0
        do v = 1, size(self%field_variables)
            if (.not. is_written(self, self%field_variables(v))) cycle
            components = 1
            if (trim(self%field_variables(v)) == "water_flux") components = PLANE_DIMENSION
            order = 0
            if (is_point_field(self, self%field_variables(v))) order = 1
            call PetscFECreateLagrange(PETSC_COMM_WORLD, PLANE_DIMENSION, &
                                       components, is_simplex, order, &
                                       1, fe, ierr)
            if (ierr == 0) call PetscObjectSetName(PetscObjectCast(fe), trim(self%field_variables(v)), ierr)
            if (ierr == 0) call DMSetField(self%output_dm, count, PETSC_NULL_DMLABEL, &
                                           PetscObjectCast(fe), ierr)
            if (ierr /= 0) error stop "PETSc could not build an output field"
            call PetscFEDestroy(fe, ierr)
            count = count + 1
        end do
        call DMCreateDS(self%output_dm, ierr)
        if (ierr /= 0) error stop "PETSc could not build the output fields"
        self%has_output_dm = .true.
    end subroutine build_output_dm

    !> Whether this run can write `name`; vapour is designed later.
    logical function is_written(self, name)
        implicit none
        class(type_ftcms), intent(in) :: self
        character(len=*), intent(in) :: name

        select case (trim(name))
        case ("temperature", "thermal_conductivity", "volumetric_heat_capacity")
            is_written = self%temperature_field >= 0 .or. self%pressure_field >= 0
        case ("pressure", "water_content", "ice_content", "hydraulic_conductivity", "vapor_content")
            is_written = .true.
        case ("water_flux")
            is_written = self%pressure_field >= 0
        case default
            is_written = .false.
        end select
    end function is_written

    !> A potential in a `lagrange` space is written where it lives, on the points.
    logical function is_point_field(self, name)
        implicit none
        class(type_ftcms), intent(in) :: self
        character(len=*), intent(in) :: name

        select case (trim(name))
        case ("temperature")
            is_point_field = self%temperature_field >= 0
        case ("pressure")
            is_point_field = self%pressure_field >= 0 .and. self%flux_field < 0
        case default
            is_point_field = .false.
        end select
    end function is_point_field

    !> One field file at the accepted state.
    subroutine write_fields(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscViewer :: viewer
        Vec :: local, global
        PetscErrorCode :: ierr
        character(len=:), allocatable :: file_name

        call DMGetLocalVector(self%output_dm, local, ierr)
        if (ierr /= 0) error stop "PETSc could not create an output vector"
        call fill_output(self, local)
        call DMGetGlobalVector(self%output_dm, global, ierr)
        ! The VTU writer prefixes every array with the vector's name: none, so
        ! the arrays carry the Output.json names alone. The HDF5 writer names
        ! its datasets by field and needs a name of its own.
        if (trim(self%field_format) == "hdf5") then
            call PetscObjectSetName(PetscObjectCast(global), "fields", ierr)
        else
            call PetscObjectSetName(PetscObjectCast(global), "", ierr)
        end if
        call DMLocalToGlobalBegin(self%output_dm, local, INSERT_VALUES, global, ierr)
        call DMLocalToGlobalEnd(self%output_dm, local, INSERT_VALUES, global, ierr)
        call DMRestoreLocalVector(self%output_dm, local, ierr)

        if (trim(self%field_format) == "hdf5") then
            file_name = "Out_"//counter(self%field_count)//".h5"
            call PetscViewerHDF5Open(PETSC_COMM_WORLD, self%output_directory//"Files/"//file_name, &
                                     FILE_MODE_WRITE, viewer, ierr)
            ! The visualisation layout: one dataset per field under
            ! /vertex_fields or /cell_fields, which petsc_gen_xdmf.py describes.
            if (ierr == 0) call PetscViewerPushFormat(viewer, PETSC_VIEWER_HDF5_VIZ, ierr)
            if (ierr == 0) call DMView(self%output_dm, viewer, ierr)
            if (ierr == 0) call VecView(global, viewer, ierr)
            if (ierr == 0) call PetscViewerPopFormat(viewer, ierr)
        else
            file_name = "Out_"//counter(self%field_count)//".vtu"
            call PetscViewerVTKOpen(PETSC_COMM_WORLD, self%output_directory//"Files/"//file_name, &
                                    FILE_MODE_WRITE, viewer, ierr)
            if (ierr == 0) call VecView(global, viewer, ierr)
        end if
        if (ierr /= 0) error stop "PETSc could not write a field file"
        call PetscViewerDestroy(viewer, ierr)
        if (ierr == 0) call DMRestoreGlobalVector(self%output_dm, global, ierr)
        if (ierr /= 0) error stop "PETSc could not close a field file"

        if (trim(self%field_format) /= "hdf5") call write_collection(self, file_name)
        self%field_count = self%field_count + 1
    end subroutine write_fields

    pure function counter(value) result(text)
        implicit none
        integer(int32), intent(in) :: value
        character(len=:), allocatable :: text

        character(len=16) :: buffer

        write (buffer, '(i3.3)') value
        text = trim(buffer)
    end function counter

    !> Rewrite `Out.pvd` with every file written so far, in output time units.
    subroutine write_collection(self, file_name)
        implicit none
        class(type_ftcms), intent(inout) :: self
        character(len=*), intent(in) :: file_name

        integer(int32) :: rank, mpi_ierr, unit, status, k

        if (.not. allocated(self%collection_times)) allocate (self%collection_times(0))
        self%collection_times = [self%collection_times, self%time/self%field_time_unit]
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, mpi_ierr)
        if (rank /= 0) return
        open (newunit=unit, file=self%output_directory//"Files/Out.pvd", status="replace", action="write", &
              iostat=status)
        if (status /= 0) error stop "The field collection file could not be written"
        write (unit, '(a)') '<?xml version="1.0"?>'
        write (unit, '(a)') '<VTKFile type="Collection" version="0.1">'
        write (unit, '(a)') '  <Collection>'
        do k = 1, size(self%collection_times)
            write (unit, '(a,es24.16,a,a,a)') '    <DataSet timestep="', self%collection_times(k), &
                '" file="Out_', counter(k - 1), '.vtu"/>'
        end do
        write (unit, '(a)') '  </Collection>'
        write (unit, '(a)') '</VTKFile>'
        close (unit)
        if (len(file_name) < 0) return
    end subroutine write_collection

    !> Fill the output vector: point fields from the solution, cell fields sampled.
    subroutine fill_output(self, local)
        implicit none
        class(type_ftcms), intent(inout) :: self
        Vec, intent(inout) :: local

        PetscSection :: section, solution_section
        PetscScalar, pointer :: raw(:), solution_raw(:)
        PetscInt :: point, p_start, p_end, cell, c_start, c_end, offset, solution_offset, num_dofs, field
        PetscErrorCode :: ierr
        real(real64) :: sample(NUM_SAMPLES), reference(PLANE_DIMENSION)
        integer(int32) :: v, f
        logical :: ok

        call refresh_local_solution(self)
        call DMGetLocalSection(self%output_dm, section, ierr)
        if (ierr == 0) call DMGetLocalSection(self%dm, solution_section, ierr)
        if (ierr == 0) call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(self%dm, 0, c_start, c_end, ierr)
        if (ierr == 0) call VecGetArray(local, raw, ierr)
        if (ierr == 0) call VecGetArrayRead(self%local_solution, solution_raw, ierr)
        if (ierr /= 0) error stop "PETSc could not open the output vectors"
        raw = 0.0d0

        reference = reference_centroid(self)
        f = -1
        do v = 1, size(self%field_variables)
            if (.not. is_written(self, self%field_variables(v))) cycle
            f = f + 1
            field = f
            if (is_point_field(self, self%field_variables(v))) then
                do point = p_start, p_end - 1
                    call PetscSectionGetFieldDof(section, point, field, num_dofs, ierr)
                    if (num_dofs <= 0) cycle
                    call PetscSectionGetFieldOffset(section, point, field, offset, ierr)
                    if (trim(self%field_variables(v)) == "temperature") then
                        call PetscSectionGetFieldOffset(solution_section, point, self%temperature_field, &
                                                        solution_offset, ierr)
                        raw(offset + 1) = self%temperature_scale*solution_raw(solution_offset + 1)
                    else
                        call PetscSectionGetFieldOffset(solution_section, point, self%pressure_field, &
                                                        solution_offset, ierr)
                        raw(offset + 1) = self%pressure_scale*solution_raw(solution_offset + 1)
                    end if
                end do
            end if
        end do
        call VecRestoreArrayRead(self%local_solution, solution_raw, ierr)

        do cell = c_start, c_end - 1
            call sample_state(self, cell, reference, sample, ok)
            if (.not. ok) sample = ieee_nan()
            f = -1
            do v = 1, size(self%field_variables)
                if (.not. is_written(self, self%field_variables(v))) cycle
                f = f + 1
                if (is_point_field(self, self%field_variables(v))) cycle
                call PetscSectionGetFieldOffset(section, cell, f, offset, ierr)
                select case (trim(self%field_variables(v)))
                case ("water_flux")
                    raw(offset + 1) = sample(S_FLUX_X)
                    raw(offset + 2) = sample(S_FLUX_Z)
                case default
                    raw(offset + 1) = sample(sample_index(self%field_variables(v)))
                end select
            end do
        end do
        call VecRestoreArray(local, raw, ierr)
    end subroutine fill_output

    !> Where in a sample an `Output.json` scalar is.
    pure integer(int32) function sample_index(name)
        implicit none
        character(len=*), intent(in) :: name

        select case (trim(name))
        case ("temperature")
            sample_index = S_TEMPERATURE
        case ("pressure")
            sample_index = S_PRESSURE
        case ("water_content")
            sample_index = S_WATER
        case ("ice_content")
            sample_index = S_ICE
        case ("vapor_content")
            sample_index = S_VAPOR
        case ("thermal_conductivity")
            sample_index = S_THERMAL_CONDUCTIVITY
        case ("volumetric_heat_capacity")
            sample_index = S_HEAT_CAPACITY
        case default
            sample_index = S_HYDRAULIC_CONDUCTIVITY
        end select
    end function sample_index

    function ieee_nan() result(value)
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        implicit none
        real(real64) :: value

        value = ieee_value(1.0d0, ieee_quiet_nan)
    end function ieee_nan

    !> The centroid of the reference cell PETSc maps every cell from.
    function reference_centroid(self) result(reference)
        implicit none
        class(type_ftcms), intent(in) :: self
        real(real64) :: reference(PLANE_DIMENSION)

        PetscBool :: is_simplex
        PetscErrorCode :: ierr

        call DMPlexIsSimplex(self%dm, is_simplex, ierr)
        reference = 0.0d0
        if (is_simplex .eqv. PETSC_TRUE) reference = -1.0d0/3.0d0
    end function reference_centroid

    !> The accepted state as a local vector with its boundary values.
    subroutine refresh_local_solution(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscErrorCode :: ierr

        call DMGlobalToLocalBegin(self%dm, self%old_solution, INSERT_VALUES, self%local_solution, ierr)
        if (ierr == 0) call DMGlobalToLocalEnd(self%dm, self%old_solution, INSERT_VALUES, self%local_solution, ierr)
        if (ierr == 0) call DMPlexInsertBoundaryValues(self%dm, PETSC_TRUE, self%local_solution, &
                                                       self%time, PETSC_NULL_VEC, &
                                                       PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
        if (ierr /= 0) error stop "PETSc could not read the accepted state"
    end subroutine refresh_local_solution

    !> Every output quantity at one reference point of one cell.
    subroutine sample_state(self, cell, reference, sample, ok)
        implicit none
        class(type_ftcms), intent(inout) :: self
        PetscInt, intent(in) :: cell
        real(real64), intent(in) :: reference(PLANE_DIMENSION)
        real(real64), intent(inout) :: sample(NUM_SAMPLES)
        logical, intent(inout) :: ok

        type(type_closure_response) :: response
        type(type_closure_tangents) :: tangents
        type(c_ptr) :: dm, local
        real(c_double) :: points(PLANE_DIMENSION, 1), values(PLANE_DIMENSION, 1), divergence(1)
        real(c_double) :: gradients(PLANE_DIMENSION, 1, 1)
        real(real64) :: temperature_gradient(PLANE_DIMENSION), pressure_gradient(PLANE_DIMENSION)
        real(real64) :: inverse_mobility(PLANE_DIMENSION*PLANE_DIMENSION), flux(PLANE_DIMENSION)
        integer(int32) :: slot
        integer(c_int) :: cerr
        PetscErrorCode :: ierr

        sample = 0.0d0
        ok = .false.
        dm = transfer(self%dm%v, dm)
        local = transfer(self%local_solution%v, local)
        points(:, 1) = reference
        cerr = 0
        temperature_gradient = 0.0d0
        pressure_gradient = 0.0d0

        sample(S_TEMPERATURE) = self%fixed_temperature
        if (self%temperature_field >= 0) then
            call petscfe_evaluate_cell_gradient(dm, local, int(self%temperature_field, c_int), int(cell, c_int), &
                                                points, values(1:1, :), gradients, cerr)
            sample(S_TEMPERATURE) = self%temperature_scale*values(1, 1)
            temperature_gradient = self%temperature_scale*gradients(:, 1, 1)
        end if
        sample(S_PRESSURE) = self%fixed_pressure
        if (cerr == 0 .and. self%pressure_field >= 0) then
            call petscfe_evaluate_cell_gradient(dm, local, int(self%pressure_field, c_int), int(cell, c_int), &
                                                points, values(1:1, :), gradients, cerr)
            sample(S_PRESSURE) = self%pressure_scale*values(1, 1)
            pressure_gradient = self%pressure_scale*gradients(:, 1, 1)
        end if
        if (cerr /= 0) return

        call cell_coefficients(self, cell, slot, inverse_mobility)
        call eval_closure_response(slot, sample(S_TEMPERATURE), sample(S_PRESSURE), response, ok)
        if (ok) call eval_closure_tangents(slot, sample(S_TEMPERATURE), sample(S_PRESSURE), tangents, ok)
        if (.not. ok) return
        sample(S_WATER) = response%liquid_content
        sample(S_ICE) = response%ice_content
        sample(S_VAPOR) = response%vapor_content
        sample(S_THERMAL_CONDUCTIVITY) = response%thermal_conductivity
        sample(S_HEAT_CAPACITY) = tangents%denthalpy_dtemperature
        sample(S_HYDRAULIC_CONDUCTIVITY) = response%relative_conductivity*closure_saturated_conductivity(slot)

        flux = 0.0d0
        if (self%flux_field >= 0) then
            call petscfe_evaluate_cell(dm, local, int(self%flux_field, c_int), int(cell, c_int), points, values, &
                                       divergence, cerr)
            flux = self%flux_scale*values(:, 1)
            ok = cerr == 0
        else if (self%nodal_flux) then
            ierr = 0
            call compute_nodal_cell_flux(self%dm, self%local_solution, cell, flux, ok, ierr)
            ok = ok .and. ierr == 0
        else if (self%pressure_field >= 0) then
            call compute_lagrange_darcy_flux(slot, sample(S_TEMPERATURE), sample(S_PRESSURE), temperature_gradient, &
                                             pressure_gradient, inverse_mobility, self%density, self%gravity, flux, ok)
        end if
        sample(S_FLUX_X) = flux(1)
        sample(S_FLUX_Z) = flux(2)
    end subroutine sample_state

    !> The closure slot and saturated inverse mobility of one cell.
    subroutine cell_coefficients(self, cell, slot, inverse_mobility)
        implicit none
        class(type_ftcms), intent(inout) :: self
        PetscInt, intent(in) :: cell
        integer(int32), intent(inout) :: slot
        real(real64), intent(inout) :: inverse_mobility(PLANE_DIMENSION*PLANE_DIMENSION)

        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: offset
        PetscErrorCode :: ierr

        call DMGetLocalSection(self%aux_dm, section, ierr)
        if (ierr == 0) call VecGetArrayRead(self%accepted_auxiliary, raw, ierr)
        if (ierr /= 0) error stop "PETSc could not read the auxiliary vector"
        call PetscSectionGetFieldOffset(section, cell, AUX_MATERIAL - 1, offset, ierr)
        slot = nint(raw(offset + 1), int32)
        call PetscSectionGetFieldOffset(section, cell, AUX_INVERSE_MOBILITY - 1, offset, ierr)
        inverse_mobility = raw(offset + 1:offset + PLANE_DIMENSION**2)
        call VecRestoreArrayRead(self%accepted_auxiliary, raw, ierr)
    end subroutine cell_coefficients

    !> Locate the observation points and open one file per requested variable.
    subroutine open_history(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        integer(int32) :: k, v, rank, mpi_ierr, status, count
        character(len=8) :: code
        character(len=16) :: unit_name

        associate (history => self%input%output_settings%history_output)
            if (trim(history%observation_type) /= "coordinates") then
                error stop "Observation points are given by coordinates; node ids do not survive the mesh distribution"
            end if
            self%num_observations = min(size(history%coordinates), MAX_OBSERVATIONS)
            do k = 1, self%num_observations
                self%observations(k)%coordinate = [history%coordinates(k)%x, history%coordinates(k)%y, &
                                                   history%coordinates(k)%z]
                if (self%computation_type == COMP_TYPES%XZ_2D%ID) then
                    self%observations(k)%plane = self%observations(k)%coordinate([1, 3])
                else
                    self%observations(k)%plane = self%observations(k)%coordinate([1, 2])
                end if
                self%observations(k)%cell = locate_cell(self, self%observations(k)%plane)
            end do

            call MPI_Comm_rank(PETSC_COMM_WORLD, rank, mpi_ierr)
            allocate (self%series(size(history%variable_names)))
            count = 0
            do v = 1, size(history%variable_names)
                if (.not. is_written(self, history%variable_names(v))) cycle
                count = count + 1
                self%series(count)%variable = history%variable_names(v)
                self%series(count)%width = merge(3, 1, trim(history%variable_names(v)) == "water_flux")
                if (rank /= 0) cycle
                call series_file(history%variable_names(v), code, unit_name)
                open (newunit=self%series(count)%unit, file=self%output_directory//"obsf_"//trim(code)//".dat", &
                      status="replace", action="write", iostat=status)
                if (status /= 0) error stop "An observation file could not be opened"
                call write_history_header(self, self%series(count), unit_name, history%output_time_unit)
            end do
            self%series = self%series(:count)
        end associate
    end subroutine open_history

    !> The file code and unit of an observed variable, as the post-processing names them.
    pure subroutine series_file(name, code, unit_name)
        implicit none
        character(len=*), intent(in) :: name
        character(len=8), intent(inout) :: code
        character(len=16), intent(inout) :: unit_name

        select case (trim(name))
        case ("temperature")
            code = "T"; unit_name = "deg C"
        case ("water_content")
            code = "Qw"; unit_name = "m3/m3"
        case ("ice_content")
            code = "Qi"; unit_name = "m3/m3"
        case ("vapor_content")
            code = "Qv"; unit_name = "m3/m3"
        case ("thermal_conductivity")
            code = "TC"; unit_name = "W/m/K"
        case ("volumetric_heat_capacity")
            code = "C"; unit_name = "J/m3/K"
        case ("pressure")
            code = "P"; unit_name = "Pa"
        case ("water_flux")
            code = "Flux"; unit_name = "m/s"
        case default
            code = "K"; unit_name = "m/s"
        end select
    end subroutine series_file

    subroutine write_history_header(self, series, unit_name, time_unit)
        implicit none
        class(type_ftcms), intent(in) :: self
        type(type_history_series), intent(in) :: series
        character(len=*), intent(in) :: unit_name
        character(len=*), intent(in) :: time_unit

        character(len=4096) :: line
        character(len=16) :: label
        integer(int32) :: k

        write (series%unit, '(a)') "# "//trim(series%variable)//" time variation"
        write (series%unit, '(a)') "#"
        write (series%unit, '(a)') "# Observation Coordinate (x,y,z)"
        do k = 1, self%num_observations
            write (series%unit, '(a,1x,i0,a,3(1x,es18.11,a))') "#    Point", k, ": (", &
                self%observations(k)%coordinate(1), ",", self%observations(k)%coordinate(2), ",", &
                self%observations(k)%coordinate(3), ")"
        end do
        write (series%unit, '(a)') "#"
        write (series%unit, '(a)') "# Output Unit: Time ["//trim(time_unit)//"], "//trim(series%variable)// &
            " ["//trim(unit_name)//"]"
        write (series%unit, '(a)') "#"
        line = "Time"
        do k = 1, self%num_observations
            write (label, '(i0)') k
            if (series%width == 3) then
                line = trim(line)//"  Obs"//trim(label)//"_x  Obs"//trim(label)//"_y  Obs"//trim(label)//"_z"
            else
                line = trim(line)//"  Obs"//trim(label)
            end if
        end do
        write (series%unit, '(a)') trim(line)
    end subroutine write_history_header

    !> The local cell holding `plane`, or -1.
    !>
    !> A brute-force test of every local cell against its affine map; the
    !> observation points are few and located once.
    function locate_cell(self, plane) result(found)
        implicit none
        class(type_ftcms), intent(in) :: self
        real(real64), intent(in) :: plane(PLANE_DIMENSION)
        PetscInt :: found

        PetscInt :: cell, c_start, c_end
        PetscBool :: is_simplex
        PetscErrorCode :: ierr
        real(c_double) :: points(PLANE_DIMENSION, 1), reference(PLANE_DIMENSION, 1)
        integer(c_int) :: cerr
        real(real64), parameter :: SLACK = 1.0d-10

        found = -1
        call DMPlexIsSimplex(self%dm, is_simplex, ierr)
        call DMPlexGetHeightStratum(self%dm, 0, c_start, c_end, ierr)
        points(:, 1) = plane
        do cell = c_start, c_end - 1
            call petscfe_reference_coordinates(transfer(self%dm%v, c_ptr_null()), int(cell, c_int), points, &
                                               reference, cerr)
            if (cerr /= 0) cycle
            if (is_simplex .eqv. PETSC_TRUE) then
                if (all(reference(:, 1) >= -1.0d0 - SLACK) .and. &
                    sum(reference(:, 1)) <= SLACK) found = cell
            else
                if (all(abs(reference(:, 1)) <= 1.0d0 + SLACK)) found = cell
            end if
            if (found >= 0) return
        end do
    end function locate_cell

    function c_ptr_null() result(pointer)
        use, intrinsic :: iso_c_binding, only: c_null_ptr
        implicit none
        type(c_ptr) :: pointer

        pointer = c_null_ptr
    end function c_ptr_null

    !> One line per series: each observation sampled on the rank that holds it.
    subroutine write_history(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        real(real64) :: local(NUM_SAMPLES, MAX_OBSERVATIONS), total(NUM_SAMPLES, MAX_OBSERVATIONS)
        real(real64) :: held(MAX_OBSERVATIONS), holders(MAX_OBSERVATIONS)
        real(real64) :: reference(PLANE_DIMENSION, 1), values(3*MAX_OBSERVATIONS)
        real(c_double) :: points(PLANE_DIMENSION, 1)
        integer(int32) :: k, s, rank, mpi_ierr, n
        integer(c_int) :: cerr
        logical :: ok

        n = self%num_observations
        if (n == 0 .or. size(self%series) == 0) return
        call refresh_local_solution(self)
        local = 0.0d0
        held = 0.0d0
        do k = 1, n
            if (self%observations(k)%cell < 0) cycle
            points(:, 1) = self%observations(k)%plane
            call petscfe_reference_coordinates(transfer(self%dm%v, c_ptr_null()), &
                                               int(self%observations(k)%cell, c_int), points, reference, cerr)
            call sample_state(self, self%observations(k)%cell, reference(:, 1), local(:, k), ok)
            if (.not. ok) local(:, k) = ieee_nan()
            held(k) = 1.0d0
        end do
        call MPI_Reduce(local, total, NUM_SAMPLES*MAX_OBSERVATIONS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                        PETSC_COMM_WORLD, mpi_ierr)
        call MPI_Reduce(held, holders, MAX_OBSERVATIONS, MPI_DOUBLE_PRECISION, MPI_SUM, 0, PETSC_COMM_WORLD, mpi_ierr)
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, mpi_ierr)
        if (rank /= 0) return

        do s = 1, size(self%series)
            do k = 1, n
                if (holders(k) > 0.0d0) total(:, k) = total(:, k)/holders(k)
                if (holders(k) <= 0.0d0) total(:, k) = ieee_nan()
                if (self%series(s)%width == 3) then
                    values(3*k - 2:3*k) = physical_vector(self, total(S_FLUX_X:S_FLUX_Z, k))
                else
                    values(k) = total(sample_index(self%series(s)%variable), k)
                end if
            end do
            if (self%series(s)%width == 3) then
                write (self%series(s)%unit, '(es24.15,*(es24.15))') self%time/self%history_time_unit, values(:3*n)
            else
                write (self%series(s)%unit, '(es24.15,*(es24.15))') self%time/self%history_time_unit, values(:n)
            end if
            flush (self%series(s)%unit)
        end do
    end subroutine write_history

    !> A vector of the intrinsic plane in the project's (x, y, z).
    pure function physical_vector(self, plane) result(vector)
        implicit none
        class(type_ftcms), intent(in) :: self
        real(real64), intent(in) :: plane(PLANE_DIMENSION)
        real(real64) :: vector(3)

        vector = 0.0d0
        if (self%computation_type == COMP_TYPES%XZ_2D%ID) then
            vector([1, 3]) = plane
        else
            vector([1, 2]) = plane
        end if
    end function physical_vector

    module subroutine destroy_output(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscErrorCode :: ierr
        integer(int32) :: s

        if (self%has_output_dm) call DMDestroy(self%output_dm, ierr)
        self%has_output_dm = .false.
        if (allocated(self%series)) then
            do s = 1, size(self%series)
                if (self%series(s)%unit >= 0) close (self%series(s)%unit)
            end do
            deallocate (self%series)
        end if
    end subroutine destroy_output

end submodule ftcms_output
