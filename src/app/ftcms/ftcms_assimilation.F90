#include <petsc/finclude/petscdmplex.h>

!> Coupling the atmospheric ensemble (ETKF) to the soil.
!>
!> At the start and after every accepted step, the surface means of \(T\) and
!> \(p_w\) over the boundary `DataAssimilation.json` names are computed as
!> boundary integrals, reduced over all ranks, and handed to the controller.
!> The controller's analysis goes back through the dynamic data providers of
!> the conditions on that boundary, which the boundary kernels evaluate at the
!> next solve.
submodule(app_ftcms) ftcms_assimilation
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
    use :: stdlib_logger, only:global_logger
    use :: governing_atmosphere, only:type_da_config
    use :: core_interop_petscfe, only:petscfe_evaluate_cell, petscfe_reference_coordinates
    use :: app_ftcms_darcy, only:PLANE_DIMENSION
    use :: app_ftcms_boundary_registry, only:boundary_manager
    implicit none

contains

    module subroutine initialize_assimilation(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        type(type_da_config) :: config

        associate (settings => self%input%assimilation)
            self%assimilation_enabled = settings%enabled
            if (.not. settings%enabled) return
            if (.not. (self%is_active(PHYSICS_TYPES%THERMAL%ID) .and. self%is_active(PHYSICS_TYPES%HYDRAULIC%ID))) then
                error stop "Assimilation exchanges heat and water with the surface: both equations must be solved"
            end if
            self%surface_boundary = settings%surface_boundary
            config%ensemble_size = settings%ensemble_size
            config%max_height = settings%max_height
            config%num_nodes = settings%num_nodes
            config%csv_file = settings%csv_file
            config%upper_bc_file = settings%upper_bc_file
            config%interval_seconds = settings%interval_seconds
            config%sigma_T = settings%sigma_T
            config%sigma_q = settings%sigma_q
            config%sigma_U = settings%sigma_U
            config%z0 = settings%z0
            config%albedo = settings%albedo
            config%emissivity = settings%emissivity
            config%Pmin = settings%Pmin
            config%Pmax = settings%Pmax
            config%lambda_soil = settings%lambda_soil
            config%stomatal_resistance = settings%stomatal_resistance
            config%reference_datetime = settings%reference_datetime
            config%latitude = settings%latitude
            config%longitude = settings%longitude
            config%tau_atm = settings%tau_atm
            config%utc_offset_hours = settings%utc_offset_hours
        end associate
        call self%assimilation%initialize(config, bc_entity_thermal=self%surface_boundary, &
                                          bc_entity_hydraulic=self%surface_boundary)
        call global_logger%log_information("Data assimilation on boundary "//trim(int_text(self%surface_boundary)))
    end subroutine initialize_assimilation

    module subroutine assimilate(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        real(real64) :: temperature, pressure

        if (.not. self%assimilation_enabled) return
        call surface_means(self, temperature, pressure)
        call self%assimilation%set_surface_state(temperature, pressure)
        call self%assimilation%execute_assimilation_cycle(self%time, 1.0d0 + self%time/86400.0d0, &
                                                          boundary_manager(PHYSICS_TYPES%THERMAL%ID), &
                                                          boundary_manager(PHYSICS_TYPES%HYDRAULIC%ID))
    end subroutine assimilate

    !> \(\int_\Gamma T\,ds/|\Gamma|\) and the same for \(p_w\), over all ranks.
    !>
    !> Each owned facet contributes its length times the field at its centroid,
    !> evaluated in the facet's one cell.
    subroutine surface_means(self, temperature, pressure)
        implicit none
        class(type_ftcms), intent(inout) :: self
        real(real64), intent(inout) :: temperature, pressure

        IS :: facets
        PetscInt, pointer :: points(:), support(:)
        PetscInt :: k, global_start, global_end
        PetscReal :: area, centroid(3), normal(3)
        PetscErrorCode :: ierr
        real(c_double) :: plane(PLANE_DIMENSION, 1), reference(PLANE_DIMENSION, 1), values(PLANE_DIMENSION, 1)
        real(c_double) :: divergence(1)
        real(real64) :: local(3), total(3)
        integer(c_int) :: cerr
        type(c_ptr) :: dm, solution

        call DMGlobalToLocalBegin(self%dm, self%old_solution, INSERT_VALUES, self%local_solution, ierr)
        call DMGlobalToLocalEnd(self%dm, self%old_solution, INSERT_VALUES, self%local_solution, ierr)
        call DMPlexInsertBoundaryValues(self%dm, PETSC_TRUE, self%local_solution, self%time, &
                                        PETSC_NULL_VEC, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
        dm = transfer(self%dm%v, dm)
        solution = transfer(self%local_solution%v, solution)

        local = 0.0d0
        if (.not. stratum_is_empty(self%dm, self%surface_boundary)) then
            call DMGetStratumIS(self%dm, "Face Sets", self%surface_boundary, facets, ierr)
            if (ierr /= 0) error stop "PETSc could not read the surface boundary"
            call ISGetIndices(facets, points, ierr)
            do k = 1, size(points)
                call DMPlexGetPointGlobal(self%dm, points(k), global_start, global_end, ierr)
                if (global_start < 0) cycle
                call DMPlexGetSupport(self%dm, points(k), support, ierr)
                if (size(support) /= 1) cycle
                call DMPlexComputeCellGeometryFVM(self%dm, points(k), area, centroid, normal, ierr)
                plane(:, 1) = real(centroid(1:PLANE_DIMENSION), c_double)
                call petscfe_reference_coordinates(dm, int(support(1), c_int), plane, reference, cerr)
                call petscfe_evaluate_cell(dm, solution, int(self%temperature_field, c_int), int(support(1), c_int), &
                                           reference, values(1:1, :), divergence, cerr)
                local(1) = local(1) + area*self%temperature_scale*values(1, 1)
                call petscfe_evaluate_cell(dm, solution, int(self%pressure_field, c_int), int(support(1), c_int), &
                                           reference, values(1:1, :), divergence, cerr)
                local(2) = local(2) + area*self%pressure_scale*values(1, 1)
                local(3) = local(3) + area
                call DMPlexRestoreSupport(self%dm, points(k), support, ierr)
            end do
            call ISRestoreIndices(facets, points, ierr)
            call ISDestroy(facets, ierr)
        end if
        call MPI_Allreduce(local, total, 3, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, ierr)
        if (.not. (total(3) > 0.0d0)) error stop "The surface boundary of DataAssimilation.json has no facets"
        temperature = total(1)/total(3)
        pressure = total(2)/total(3)
    end subroutine surface_means

    pure function int_text(value) result(text)
        implicit none
        integer(int32), intent(in) :: value
        character(len=16) :: text

        write (text, '(i0)') value
    end function int_text

end submodule ftcms_assimilation
