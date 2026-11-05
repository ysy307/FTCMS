#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscsnes.h>

!> Building the engine from the input, and taking it down again.
submodule(app_ftcms) ftcms_setup
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_funloc, c_null_funptr, c_null_ptr
    use :: stdlib_logger, only:global_logger, information_level, warning_level
    use :: module_core, only:FE_SPACES, COEFFICIENT_EVALUATIONS, COMP_TYPES, IC_TARGETS, IC_METHODS, &
                              type_config_ic, type_config_bc, type_config_time, type_config_time_ats, &
                              THERMAL_BC_TYPES, HYDRAULIC_BC_TYPES
    use :: core_physical_constants, only:gravity_acceleration, reference_water_density
    use :: module_input, only:input_translator
    use :: core_interop_petscds, only:petscds_set_residual, petscds_set_jacobian, &
                                      petscds_add_boundary_essential, petscds_add_boundary_natural

    use :: app_ftcms_darcy, only:PLANE_DIMENSION, invert_mobility
    use :: app_ftcms_closure, only:install_closure, clear_closure, closure_saturated_conductivity, &
                                   configure_closure_cache
    use :: app_ftcms_kernels
    use :: app_ftcms_kernels_lagrange, only:f0_lagrange_conservation, f1_lagrange_conservation, &
                                            g2_lagrange_pressure_pressure, g3_lagrange_pressure_pressure, &
                                            g2_lagrange_pressure_temperature, g3_lagrange_pressure_temperature
    use :: app_ftcms_nodal_flux, only:install_nodal_flux, clear_nodal_flux
    use :: app_ftcms_kernels_advection, only:f0_temperature_advected, g1_temperature_temperature, &
                                             g0_temperature_flux, g1_temperature_pressure
    use :: app_ftcms_boundary_registry, only:install_boundary_conditions, clear_boundary_registry, &
                                             add_boundary_slot, boundary_slot_context, set_boundary_time, &
                                             boundary_condition_count, boundary_condition_id, &
                                             bd_f0_natural, bd_g0_natural, bd_f0_flux_potential, &
                                             bc_essential_potential, bc_essential_normal_flux, bc_essential_zero, &
                                             bd_f0_seepage, set_slot_facets
    implicit none

    !> The project input, held here so the engine can point at it.
    type(type_input), save, target :: project_input

    !> The quadrature every field and auxiliary field shares: the RT mass
    !> integral needs order two on an affine cell.
    PetscInt, parameter :: QUADRATURE_ORDER = 2

    !> The options prefix of the nonlinear solver and everything inside it.
    character(len=*), parameter :: SOLVER_PREFIX = "ftcms_"

    !> The label the closed default of a mixed field is built on.
    character(len=*), parameter :: CLOSED_LABEL = "ftcms_closed_boundary"

    !> The role a condition plays: a prescribed potential, an inflow, or a switch between them.
    integer(int32), parameter :: CONDITION_POTENTIAL = 1
    integer(int32), parameter :: CONDITION_FLUX = 2
    integer(int32), parameter :: CONDITION_SWITCHING = 3

    ! PetscFPTrapPush/Pop have no Fortran module interface; these are the C
    ! stubs every PETSc Fortran call goes through. Building the trimmed (RT)
    ! dual space touches an inert floating point exception inside PETSc, which
    ! a build that unmasks traps would turn fatal.
    integer(c_int), parameter :: PETSC_FP_TRAP_OFF_LOCAL = 0_c_int
    interface
        subroutine petscfptrappush_(trap, trap_ierr) bind(c, name="petscfptrappush_")
            import :: c_int
            integer(c_int), intent(in) :: trap
            integer(c_int), intent(inout) :: trap_ierr
        end subroutine petscfptrappush_
        subroutine petscfptrappop_(trap_ierr) bind(c, name="petscfptrappop_")
            import :: c_int
            integer(c_int), intent(inout) :: trap_ierr
        end subroutine petscfptrappop_
    end interface

contains

    module subroutine initialize_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout), target :: self

        PetscErrorCode :: ierr
        integer(int32) :: rank, mpi_ierr

        call project_input%initialize()
        self%input => project_input
        call project_input%display()
        ! Information is global and comes from rank 0; a warning is about one rank's own cells.
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, mpi_ierr)
        if (project_input%output_settings%standard_output%print_progress .and. rank == 0) then
            call global_logger%configure(level=information_level, time_stamp=.true., max_width=0)
        else
            call global_logger%configure(level=warning_level, time_stamp=.true., max_width=0)
        end if

        call choose_discretization(self)
        call build_field_table(self)
        call build_dm(self)
        call install_materials(self)
        call choose_scales(self)
        call build_fields(self)
        call register_kernels(self)
        ierr = 0
        call store_constants(self, ierr)
        if (ierr /= 0) error stop "PETSc rejected the run constants"
        call build_auxiliary(self)
        call add_boundaries(self)
        call read_time_controls(self)
        call set_initial_state(self)
        call build_solver(self)
        self%ready = .true.

        call initialize_assimilation(self)
        call initialize_output(self)
        call write_output(self, force=.true.)
    end subroutine initialize_ftcms

    !> Read which equations are solved and how each is discretised.
    !>
    !> Fails by `error stop` on vapour transport, which is designed later and
    !> must not be ignored silently, and on a combination without kernels.
    subroutine choose_discretization(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        integer(int32) :: id

        associate (basic => self%input%basic)
            self%computation_type = basic%simulation_settings%calculate_type
            do id = 1, PHYSICS_TYPES%NUM_ID
                self%is_active(id) = basic%analysis_controls%is_active(id)
            end do
            self%space = basic%discretization_settings%space
            self%degree = basic%discretization_settings%degree
            self%evaluation = basic%discretization_settings%coefficient_evaluation
            if (.not. basic%discretization_settings%is_given) then
                call global_logger%log_information("discretization_settings absent: every equation is lagrange, "// &
                                                   "degree 1, quadrature")
            end if
            if (basic%analysis_controls%enable_vapor_transport) then
                error stop "Vapour transport is not implemented in this engine; set enable_vapor_transport false"
            end if
            ! The segregation sink and its latent heat (ice-lens growth) have
            ! no kernels; a material asking for them must not run without them.
            do id = 1, basic%num_materials
                if (basic%materials(id)%phase_change%equilibrium_model%segregation) then
                    error stop "Segregation (ice-lens growth) is not implemented in this engine; set segregation false"
                end if
            end do
        end associate

        do id = 1, PHYSICS_TYPES%NUM_ID
            if (.not. self%is_active(id)) cycle
            if (.not. is_implemented(self, id)) then
                write (*, '(a)') "Not implemented: "//PHYSICS_TYPES%to_name(id)//" in "// &
                    FE_SPACES%to_name(self%space(id))//" of degree "//to_text(self%degree(id))//" with "// &
                    COEFFICIENT_EVALUATIONS%to_name(self%evaluation(id))//" coefficients"
                error stop "The discretisation_settings name a combination without kernels"
            end if
        end do
        if (.not. (self%is_active(PHYSICS_TYPES%THERMAL%ID) .or. self%is_active(PHYSICS_TYPES%HYDRAULIC%ID))) then
            error stop "Kernels exist for the heat and water equations only"
        end if
    end subroutine choose_discretization

    !> Whether kernels exist for equation `id` in its chosen discretisation.
    logical function is_implemented(self, id)
        implicit none
        class(type_ftcms), intent(in) :: self
        integer(int32), intent(in) :: id

        logical :: lagrange_one, quadrature

        lagrange_one = self%space(id) == FE_SPACES%LAGRANGE%ID .and. self%degree(id) == 1
        quadrature = self%evaluation(id) == COEFFICIENT_EVALUATIONS%QUADRATURE%ID
        is_implemented = .false.
        if (id == PHYSICS_TYPES%THERMAL%ID) then
            is_implemented = lagrange_one .and. quadrature
        else if (id == PHYSICS_TYPES%HYDRAULIC%ID) then
            is_implemented = (self%space(id) == FE_SPACES%MIXED%ID .and. self%degree(id) == 0 .and. quadrature) &
                             .or. (lagrange_one .and. quadrature) &
                             ! The nodal flux takes the nodal temperature with the pressure.
                             .or. (lagrange_one .and. self%evaluation(id) == COEFFICIENT_EVALUATIONS%NODAL%ID &
                                   .and. self%is_active(PHYSICS_TYPES%THERMAL%ID))
        end if
    end function is_implemented

    !> One table entry per field, equations in `PHYSICS_TYPES` order.
    subroutine build_field_table(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        integer(int32) :: id

        self%num_fields = 0
        do id = 1, PHYSICS_TYPES%NUM_ID
            if (.not. self%is_active(id)) cycle
            if (self%space(id) == FE_SPACES%MIXED%ID) call add_field(self, id, ROLE_FLUX)
            call add_field(self, id, ROLE_POTENTIAL)
        end do
        self%nodal_flux = self%is_active(PHYSICS_TYPES%HYDRAULIC%ID) .and. &
                          self%evaluation(PHYSICS_TYPES%HYDRAULIC%ID) == COEFFICIENT_EVALUATIONS%NODAL%ID
    end subroutine build_field_table

    subroutine add_field(self, id, role)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: id
        integer(int32), intent(in) :: role

        character(len=32) :: name

        name = field_name(id, role)
        self%num_fields = self%num_fields + 1
        self%fields(self%num_fields) = type_field_entry(id, role, name, self%num_fields - 1)
        if (id == PHYSICS_TYPES%THERMAL%ID .and. role == ROLE_POTENTIAL) self%temperature_field = self%num_fields - 1
        if (id == PHYSICS_TYPES%HYDRAULIC%ID .and. role == ROLE_POTENTIAL) self%pressure_field = self%num_fields - 1
        if (id == PHYSICS_TYPES%HYDRAULIC%ID .and. role == ROLE_FLUX) self%flux_field = self%num_fields - 1
    end subroutine add_field

    !> The name of a field, which is also its options prefix and output name.
    pure function field_name(id, role) result(name)
        implicit none
        integer(int32), intent(in) :: id
        integer(int32), intent(in) :: role
        character(len=32) :: name

        if (id == PHYSICS_TYPES%THERMAL%ID) then
            name = merge("temperature ", "heat_flux   ", role == ROLE_POTENTIAL)
        else if (id == PHYSICS_TYPES%HYDRAULIC%ID) then
            name = merge("pressure    ", "water_flux  ", role == ROLE_POTENTIAL)
        else
            name = PHYSICS_TYPES%to_name(id)
        end if
    end function field_name

    !> The solve DM: a clone of the input's mesh, which must be intrinsic 2D.
    subroutine build_dm(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscInt :: dimension, coordinate_dimension
        PetscErrorCode :: ierr

        call DMGetDimension(self%input%geometry%mesh%dm, dimension, ierr)
        if (ierr == 0) call DMGetCoordinateDim(self%input%geometry%mesh%dm, coordinate_dimension, ierr)
        if (ierr /= 0) error stop "PETSc could not report the mesh dimensions"
        if (dimension /= PLANE_DIMENSION .or. dimension /= coordinate_dimension) then
            error stop "Kernels exist for two-dimensional analyses only"
        end if
        call DMClone(self%input%geometry%mesh%dm, self%dm, ierr)
        if (ierr /= 0) error stop "PETSc could not clone the mesh"
    end subroutine build_dm

    !> Build the closure from the materials of `Basic.json`.
    !>
    !> One slot per material, in input order; the porosity is the porosity
    !> initial condition, which must be uniform.
    subroutine install_materials(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        type(type_config_ic) :: porosity
        real(real64), allocatable :: porosities(:)
        integer(int32) :: m, count
        PetscInt :: c_start, c_end
        PetscErrorCode :: ierr

        call input_translator%execute(self%input, IC_TARGETS%POROSITY, porosity)
        if (.not. (porosity%ic_kind == IC_METHODS%UNIFORM)) error stop "The porosity initial condition must be uniform"

        count = self%input%basic%num_materials
        allocate (self%material_ids(count), porosities(count))
        do m = 1, count
            self%material_ids(m) = self%input%basic%materials(m)%id
        end do
        porosities = porosity%value
        call install_closure(self%input, self%material_ids, porosities)

        ! Every residual and Jacobian block asks for the closure at the same
        ! quadrature states; the memo answers a repeated (slot, T, p_w) once.
        ! Sized for a few states at each of a few points per local cell.
        call DMPlexGetHeightStratum(self%dm, 0, c_start, c_end, ierr)
        if (ierr /= 0) error stop "PETSc could not count the cells"
        call configure_closure_cache(.true., 32*(c_end - c_start))
    end subroutine install_materials

    !> The units the unknowns and rows are measured in, from the problem itself.
    !>
    !> Length: the largest extent of the mesh. Flux: the largest saturated
    !> conductivity. Pressure: the initial pore pressure, or the hydrostatic
    !> head over the length when that is larger. The heat row is divided by
    !> \(R_T=1/L^2\) so its conduction term is of the order of the conductivity.
    subroutine choose_scales(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        type(type_config_ic) :: initial
        PetscReal :: low(3), high(3)
        PetscErrorCode :: ierr
        integer(int32) :: m

        call DMGetBoundingBox(self%dm, low, high, ierr)
        if (ierr /= 0) error stop "PETSc could not measure the mesh"
        self%length_scale = maxval(high(1:PLANE_DIMENSION) - low(1:PLANE_DIMENSION))
        if (.not. (self%length_scale > 0.0d0)) error stop "The mesh has no extent"

        self%flux_scale = 0.0d0
        do m = 1, size(self%material_ids)
            self%flux_scale = max(self%flux_scale, closure_saturated_conductivity(m))
        end do

        self%density = reference_water_density
        self%gravity = 0.0d0
        if (self%computation_type == COMP_TYPES%XZ_2D%ID) self%gravity(2) = -gravity_acceleration

        self%initial_temperature = 0.0d0
        self%initial_pressure = 0.0d0
        if (self%is_active(PHYSICS_TYPES%THERMAL%ID)) then
            call input_translator%execute(self%input, IC_TARGETS%THERMAL, initial)
            if (.not. (initial%ic_kind == IC_METHODS%UNIFORM)) error stop "The thermal initial condition must be uniform"
            self%initial_temperature = initial%value
        end if
        if (self%is_active(PHYSICS_TYPES%HYDRAULIC%ID)) then
            call input_translator%execute(self%input, IC_TARGETS%HYDRAULIC, initial)
            if (.not. (initial%ic_kind == IC_METHODS%UNIFORM)) error stop "The hydraulic initial condition must be uniform"
            self%initial_pressure = initial%value
        end if
        ! The potential of an equation that is not solved is held at the
        ! initial value Conditions.json states for it.
        self%fixed_temperature = self%initial_temperature
        self%fixed_pressure = self%initial_pressure
        if (.not. self%is_active(PHYSICS_TYPES%THERMAL%ID)) then
            self%fixed_temperature = stated_initial_value(self, IC_TARGETS%THERMAL%ID, "thermal")
        end if
        if (.not. self%is_active(PHYSICS_TYPES%HYDRAULIC%ID)) then
            self%fixed_pressure = stated_initial_value(self, IC_TARGETS%HYDRAULIC%ID, "hydraulic")
        end if

        self%pressure_scale = max(abs(self%initial_pressure), self%density*gravity_acceleration*self%length_scale)
        self%temperature_scale = 1.0d0
        self%thermal_scale = 1.0d0/self%length_scale**2
    end subroutine choose_scales

    !> The uniform initial value `Conditions.json` states for an equation not solved.
    !>
    !> Fails by `error stop` when there is none: the closure needs a value for
    !> every potential, and none is assumed.
    function stated_initial_value(self, target, name) result(value)
        implicit none
        class(type_ftcms), intent(in) :: self
        integer(int32), intent(in) :: target
        character(len=*), intent(in) :: name
        real(real64) :: value

        associate (stated => self%input%conditions%initial_conditions%physics(target))
            if (.not. allocated(stated%type)) then
                write (*, '(a)') "The "//name//" equation is not solved, so initial_conditions."//name// &
                    " must give the value it is held at"
                error stop "An equation that is not solved needs a uniform initial value"
            end if
            if (stated%type /= IC_METHODS%UNIFORM%NAME) error stop "A held potential must be uniform"
            value = stated%value
        end associate
    end function stated_initial_value

    !> Create one PetscFE per table entry and attach it to the DM.
    subroutine build_fields(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscFE :: fe, first
        PetscBool :: is_simplex
        PetscErrorCode :: ierr
        integer(int32) :: f
        integer(c_int) :: cerr

        call DMPlexIsSimplex(self%dm, is_simplex, ierr)
        if (ierr /= 0) error stop "PETSc could not read the cell shape"
        if ((is_simplex .neqv. PETSC_TRUE) .and. self%flux_field >= 0) call verify_rectangular_cells(self%dm)

        do f = 1, self%num_fields
            associate (entry => self%fields(f))
                call petscfptrappush_(PETSC_FP_TRAP_OFF_LOCAL, cerr)
                if (entry%role == ROLE_FLUX) then
                    call set_flux_space_options(trim(entry%name)//"_", is_simplex)
                    call PetscFECreateDefault(PETSC_COMM_WORLD, PLANE_DIMENSION, &
                                              PLANE_DIMENSION, is_simplex, &
                                              trim(entry%name)//"_", QUADRATURE_ORDER, fe, ierr)
                else
                    call PetscFECreateLagrange(PETSC_COMM_WORLD, PLANE_DIMENSION, &
                                               1, is_simplex, &
                                               self%degree(entry%physics_id), &
                                               QUADRATURE_ORDER, fe, ierr)
                end if
                call petscfptrappop_(cerr)
                if (ierr /= 0) error stop "PETSc could not build the basis of a field"
                if (f == 1) then
                    first = fe
                else
                    call PetscFECopyQuadrature(first, fe, ierr)
                    if (ierr /= 0) error stop "PETSc could not share the quadrature between fields"
                end if
                call PetscObjectSetName(PetscObjectCast(fe), trim(entry%name), ierr)
                if (ierr == 0) call DMSetField(self%dm, entry%number, PETSC_NULL_DMLABEL, PetscObjectCast(fe), ierr)
                if (ierr /= 0) error stop "PETSc could not attach a field"
                if (f > 1) call PetscFEDestroy(fe, ierr)
            end associate
        end do
        call PetscFEDestroy(first, ierr)

        call DMCreateDS(self%dm, ierr)
        if (ierr == 0) call DMGetDS(self%dm, self%ds, ierr)
        if (ierr /= 0) error stop "PETSc could not build the discrete system"
    end subroutine build_fields

    !> The options of a lowest-order Raviart-Thomas space under `prefix`.
    !>
    !> PETSc's `snes/tutorials/ex24.c` spaces: the trimmed polynomial space on
    !> a simplex, and the sum of two tensor-product spaces on a quadrilateral.
    !> Each is set only when the run has not named it, so a command line can
    !> still change it.
    subroutine set_flux_space_options(prefix, is_simplex)
        implicit none
        character(len=*), intent(in) :: prefix
        PetscBool, intent(in) :: is_simplex

        if (is_simplex .eqv. PETSC_TRUE) then
            call set_default_option("-"//prefix//"petscspace_type", "ptrimmed")
            call set_default_option("-"//prefix//"petscspace_degree", "0")
            call set_default_option("-"//prefix//"petscspace_components", "2")
            call set_default_option("-"//prefix//"petscspace_ptrimmed_form_degree", "-1")
        else
            call set_default_option("-"//prefix//"petscspace_type", "sum")
            call set_default_option("-"//prefix//"petscspace_degree", "1")
            call set_default_option("-"//prefix//"petscspace_variables", "2")
            call set_default_option("-"//prefix//"petscspace_components", "2")
            call set_default_option("-"//prefix//"petscspace_sum_spaces", "2")
            call set_default_option("-"//prefix//"petscspace_sum_concatenate", "true")
            ! The first component: degree one in x, zero in z; the second the reverse.
            call set_default_option("-"//prefix//"sumcomp_0_petscspace_variables", "2")
            call set_default_option("-"//prefix//"sumcomp_0_petscspace_type", "tensor")
            call set_default_option("-"//prefix//"sumcomp_0_petscspace_tensor_spaces", "2")
            call set_default_option("-"//prefix//"sumcomp_0_petscspace_tensor_uniform", "false")
            call set_default_option("-"//prefix//"sumcomp_0_tensorcomp_0_petscspace_degree", "1")
            call set_default_option("-"//prefix//"sumcomp_0_tensorcomp_1_petscspace_degree", "0")
            call set_default_option("-"//prefix//"sumcomp_1_petscspace_variables", "2")
            call set_default_option("-"//prefix//"sumcomp_1_petscspace_type", "tensor")
            call set_default_option("-"//prefix//"sumcomp_1_petscspace_tensor_spaces", "2")
            call set_default_option("-"//prefix//"sumcomp_1_petscspace_tensor_uniform", "false")
            call set_default_option("-"//prefix//"sumcomp_1_tensorcomp_0_petscspace_degree", "0")
            call set_default_option("-"//prefix//"sumcomp_1_tensorcomp_1_petscspace_degree", "1")
        end if
        call set_default_option("-"//prefix//"petscdualspace_type", "lagrange")
        call set_default_option("-"//prefix//"petscdualspace_order", "1")
        call set_default_option("-"//prefix//"petscdualspace_form_degree", "-1")
        call set_default_option("-"//prefix//"petscdualspace_lagrange_trimmed", "true")
        call set_default_option("-"//prefix//"petscdualspace_lagrange_continuity", "true")
    end subroutine set_flux_space_options

    !> Refuse any cell that is not an axis-aligned rectangle of positive area.
    !>
    !> The facet-flux reconstruction and the RT pushforward assume one affine
    !> map per cell with a diagonal Jacobian; a warped quadrilateral has none.
    subroutine verify_rectangular_cells(dm)
        implicit none
        DM, intent(in) :: dm

        PetscScalar, pointer :: cell_coordinates(:)
        PetscSection :: coordinate_section
        Vec :: coordinates
        DM :: coordinate_dm
        PetscInt :: cell, c_start, c_end, closure_size
        PetscErrorCode :: ierr
        real(real64) :: corners(PLANE_DIMENSION, 4), edge_a(PLANE_DIMENSION), edge_b(PLANE_DIMENSION)
        real(real64) :: tolerance
        integer(int32) :: v

        call DMPlexGetHeightStratum(dm, 0, c_start, c_end, ierr)
        if (ierr == 0) call DMGetCoordinateDM(dm, coordinate_dm, ierr)
        if (ierr == 0) call DMGetCoordinatesLocal(dm, coordinates, ierr)
        if (ierr == 0) call DMGetCoordinateSection(dm, coordinate_section, ierr)
        if (ierr /= 0) error stop "PETSc could not read the cell coordinates"

        do cell = c_start, c_end - 1
            nullify (cell_coordinates)
            closure_size = 0
            call DMPlexVecGetClosure(coordinate_dm, coordinate_section, coordinates, cell, &
                                     closure_size, cell_coordinates, ierr)
            if (ierr /= 0 .or. closure_size /= 4*PLANE_DIMENSION) then
                error stop "A mixed field on quadrilaterals needs cells of four vertices"
            end if
            do v = 1, 4
                corners(:, v) = cell_coordinates((v - 1)*PLANE_DIMENSION + 1:v*PLANE_DIMENSION)
            end do
            call DMPlexVecRestoreClosure(coordinate_dm, coordinate_section, coordinates, cell, &
                                         closure_size, cell_coordinates, ierr)
            edge_a = corners(:, 2) - corners(:, 1)
            edge_b = corners(:, 4) - corners(:, 1)
            ! A geometric tolerance: mesh files carry coordinates rounded near 1e-14.
            tolerance = 1.0d-9*max(norm2(edge_a), norm2(edge_b))
            if (maxval(abs(corners(:, 3) - (corners(:, 1) + edge_a + edge_b))) > tolerance &
                .or. edge_a(1)*edge_b(2) - edge_a(2)*edge_b(1) <= 0.0d0 &
                .or. .not. ((abs(edge_a(2)) <= tolerance .and. abs(edge_b(1)) <= tolerance) .or. &
                            (abs(edge_a(1)) <= tolerance .and. abs(edge_b(2)) <= tolerance))) then
                error stop "A mixed field on quadrilaterals needs axis-aligned rectangles"
            end if
        end do
    end subroutine verify_rectangular_cells

    !> Register the volume kernels the field table calls for.
    subroutine register_kernels(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        type(c_ptr) :: ds
        integer(c_int) :: t, q, p, cerr

        ds = transfer(self%ds%v, ds)
        t = int(self%temperature_field, c_int)
        q = int(self%flux_field, c_int)
        p = int(self%pressure_field, c_int)
        cerr = 0

        ! Heat is conducted, and carried by the liquid wherever water moves.
        if (t >= 0 .and. p < 0) then
            call petscds_set_residual(ds, t, c_funloc(f0_temperature), c_funloc(f1_temperature), cerr)
            call check(cerr)
            call petscds_set_jacobian(ds, t, t, c_funloc(g0_temperature_temperature), c_null_funptr, &
                                      c_funloc(g2_temperature_temperature), c_funloc(g3_temperature_temperature), cerr)
            call check(cerr)
        else if (t >= 0) then
            call petscds_set_residual(ds, t, c_funloc(f0_temperature_advected), c_funloc(f1_temperature), cerr)
            call check(cerr)
            call petscds_set_jacobian(ds, t, t, c_funloc(g0_temperature_temperature), &
                                      c_funloc(g1_temperature_temperature), &
                                      c_funloc(g2_temperature_temperature), c_funloc(g3_temperature_temperature), cerr)
            call check(cerr)
            if (q >= 0) then
                call petscds_set_jacobian(ds, t, p, c_funloc(g0_temperature_pressure), c_null_funptr, &
                                          c_funloc(g2_temperature_pressure), c_null_funptr, cerr)
                call check(cerr)
                call petscds_set_jacobian(ds, t, q, c_funloc(g0_temperature_flux), c_null_funptr, &
                                          c_null_funptr, c_null_funptr, cerr)
            else
                call petscds_set_jacobian(ds, t, p, c_funloc(g0_temperature_pressure), &
                                          c_funloc(g1_temperature_pressure), &
                                          c_funloc(g2_temperature_pressure), c_null_funptr, cerr)
            end if
            call check(cerr)
        end if

        if (q >= 0) then
            call petscds_set_residual(ds, q, c_funloc(f0_flux), c_funloc(f1_flux), cerr)
            call check(cerr)
            call petscds_set_residual(ds, p, c_funloc(f0_conservation), c_null_funptr, cerr)
            call check(cerr)
            call petscds_set_jacobian(ds, q, q, c_funloc(g0_flux_flux), c_null_funptr, c_null_funptr, &
                                      c_null_funptr, cerr)
            call check(cerr)
            call petscds_set_jacobian(ds, q, p, c_funloc(g0_flux_pressure), c_null_funptr, &
                                      c_funloc(g2_flux_pressure), c_null_funptr, cerr)
            call check(cerr)
            call petscds_set_jacobian(ds, p, q, c_null_funptr, c_funloc(g1_pressure_flux), c_null_funptr, &
                                      c_null_funptr, cerr)
            call check(cerr)
            call petscds_set_jacobian(ds, p, p, c_funloc(g0_pressure_pressure), c_null_funptr, c_null_funptr, &
                                      c_null_funptr, cerr)
            call check(cerr)
            if (t >= 0) then
                call petscds_set_jacobian(ds, q, t, c_funloc(g0_flux_temperature), c_null_funptr, &
                                          c_funloc(g2_flux_temperature), c_null_funptr, cerr)
                call check(cerr)
                call petscds_set_jacobian(ds, p, t, c_funloc(g0_pressure_temperature), c_null_funptr, &
                                          c_null_funptr, c_null_funptr, cerr)
                call check(cerr)
            end if
        else if (p >= 0 .and. self%nodal_flux) then
            ! The Darcy term is assembled by the nodal flux, outside PETSc's element loop.
            call petscds_set_residual(ds, p, c_funloc(f0_lagrange_conservation), c_null_funptr, cerr)
            call check(cerr)
            call petscds_set_jacobian(ds, p, p, c_funloc(g0_pressure_pressure), c_null_funptr, c_null_funptr, &
                                      c_null_funptr, cerr)
            call check(cerr)
            call petscds_set_jacobian(ds, p, t, c_funloc(g0_pressure_temperature), c_null_funptr, &
                                      c_null_funptr, c_null_funptr, cerr)
            call check(cerr)
        else if (p >= 0) then
            call petscds_set_residual(ds, p, c_funloc(f0_lagrange_conservation), &
                                      c_funloc(f1_lagrange_conservation), cerr)
            call check(cerr)
            call petscds_set_jacobian(ds, p, p, c_funloc(g0_pressure_pressure), c_null_funptr, &
                                      c_funloc(g2_lagrange_pressure_pressure), &
                                      c_funloc(g3_lagrange_pressure_pressure), cerr)
            call check(cerr)
            if (t >= 0) then
                call petscds_set_jacobian(ds, p, t, c_funloc(g0_pressure_temperature), c_null_funptr, &
                                          c_funloc(g2_lagrange_pressure_temperature), &
                                          c_funloc(g3_lagrange_pressure_temperature), cerr)
                call check(cerr)
            end if
        end if
    contains
        subroutine check(code)
            integer(c_int), intent(in) :: code
            if (code /= 0) error stop "PETSc rejected a pointwise kernel"
        end subroutine check
    end subroutine register_kernels

    module subroutine store_constants(self, ierr)
        implicit none
        class(type_ftcms), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        self%constants(CONST_DT) = self%dt
        self%constants(CONST_FLUX_SCALE) = self%flux_scale
        self%constants(CONST_PRESSURE_SCALE) = self%pressure_scale
        self%constants(CONST_LENGTH_SCALE) = self%length_scale
        self%constants(CONST_DENSITY) = self%density
        self%constants(CONST_GRAVITY_X) = self%gravity(1)
        self%constants(CONST_GRAVITY_Z) = self%gravity(2)
        self%constants(CONST_TEMPERATURE_SCALE) = self%temperature_scale
        self%constants(CONST_THERMAL_SCALE) = self%thermal_scale
        self%constants(CONST_PRESSURE_FIELD) = real(self%pressure_field, real64)
        self%constants(CONST_TEMPERATURE_FIELD) = real(self%temperature_field, real64)
        self%constants(CONST_FLUX_FIELD) = real(self%flux_field, real64)
        self%constants(CONST_FIXED_TEMPERATURE) = self%fixed_temperature
        self%constants(CONST_FIXED_PRESSURE) = self%fixed_pressure
        self%constants(CONST_BDF_CURRENT) = self%bdf_weights(1)
        self%constants(CONST_BDF_ACCEPTED) = self%bdf_weights(2)
        self%constants(CONST_BDF_OLDER) = self%bdf_weights(3)

        call PetscDSSetConstants(self%ds, NUM_CONSTANTS, self%constants, ierr)
    end subroutine store_constants

    !> The cell-wise coefficients and the accepted state, as auxiliary fields.
    subroutine build_auxiliary(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscFE :: first
        PetscObject :: discretisation
        PetscErrorCode :: ierr
        PetscInt :: temperature_order, pressure_order

        call DMGetField(self%dm, 0, PETSC_NULL_DMLABEL, discretisation, ierr)
        if (ierr /= 0) error stop "PETSc could not read the first field"
        first%v = discretisation%v

        temperature_order = 0
        if (self%temperature_field >= 0) temperature_order = self%degree(PHYSICS_TYPES%THERMAL%ID)
        pressure_order = 0
        if (self%pressure_field >= 0 .and. self%flux_field < 0) then
            pressure_order = self%degree(PHYSICS_TYPES%HYDRAULIC%ID)
        end if

        call DMClone(self%dm, self%aux_dm, ierr)
        if (ierr /= 0) error stop "PETSc could not clone the mesh for the auxiliary fields"
        call add_aux_field(self%aux_dm, AUX_INVERSE_MOBILITY, "inverse_mobility", PLANE_DIMENSION**2, &
                           0, first)
        call add_aux_field(self%aux_dm, AUX_MATERIAL, "material", 1, 0, first)
        call add_aux_field(self%aux_dm, AUX_SOURCE, "source", 1, 0, first)
        call add_aux_field(self%aux_dm, AUX_TEMPERATURE_OLD, "temperature_old", 1, temperature_order, first)
        call add_aux_field(self%aux_dm, AUX_PRESSURE_OLD, "pressure_old", 1, pressure_order, first)
        call add_aux_field(self%aux_dm, AUX_BOUNDARY_SLOT, "boundary_slot", 1, 0, first)
        call add_aux_field(self%aux_dm, AUX_TEMPERATURE_OLDER, "temperature_older", 1, temperature_order, first)
        call add_aux_field(self%aux_dm, AUX_PRESSURE_OLDER, "pressure_older", 1, pressure_order, first)

        call DMCreateDS(self%aux_dm, ierr)
        if (ierr == 0) call DMCreateLocalVector(self%aux_dm, self%auxiliary, ierr)
        if (ierr == 0) call VecSet(self%auxiliary, 0.0d0, ierr)
        if (ierr == 0) call DMSetAuxiliaryVec(self%dm, PETSC_NULL_DMLABEL, 0, 0, &
                                              self%auxiliary, ierr)
        if (ierr /= 0) error stop "PETSc could not build the auxiliary vector"

        call store_materials(self)
        ! A potential the run does not solve for stays at its fixed value.
        if (self%temperature_field < 0) then
            call fill_aux_field(self%aux_dm, self%auxiliary, AUX_TEMPERATURE_OLD, self%fixed_temperature)
            call fill_aux_field(self%aux_dm, self%auxiliary, AUX_TEMPERATURE_OLDER, self%fixed_temperature)
        end if
        if (self%pressure_field < 0) then
            call fill_aux_field(self%aux_dm, self%auxiliary, AUX_PRESSURE_OLD, self%fixed_pressure)
            call fill_aux_field(self%aux_dm, self%auxiliary, AUX_PRESSURE_OLDER, self%fixed_pressure)
        end if

        call DMCreateLocalVector(self%aux_dm, self%accepted_auxiliary, ierr)
        if (ierr == 0) call VecCopy(self%auxiliary, self%accepted_auxiliary, ierr)
        if (ierr /= 0) error stop "PETSc could not copy the auxiliary vector"
    end subroutine build_auxiliary

    subroutine add_aux_field(dm, field, name, num_components, order, first)
        implicit none
        DM, intent(inout) :: dm
        integer(int32), intent(in) :: field
        character(len=*), intent(in) :: name
        integer(int32), intent(in) :: num_components
        PetscInt, intent(in) :: order
        PetscFE, intent(in) :: first

        PetscFE :: fe
        PetscBool :: is_simplex
        PetscErrorCode :: ierr
        integer(c_int) :: cerr

        call DMPlexIsSimplex(dm, is_simplex, ierr)
        call petscfptrappush_(PETSC_FP_TRAP_OFF_LOCAL, cerr)
        call PetscFECreateLagrange(PETSC_COMM_WORLD, PLANE_DIMENSION, &
                                   num_components, is_simplex, order, QUADRATURE_ORDER, fe, ierr)
        call petscfptrappop_(cerr)
        if (ierr == 0) call PetscFECopyQuadrature(first, fe, ierr)
        if (ierr == 0) call PetscObjectSetName(PetscObjectCast(fe), name, ierr)
        if (ierr == 0) call DMSetField(dm, field - 1, PETSC_NULL_DMLABEL, PetscObjectCast(fe), ierr)
        if (ierr /= 0) error stop "PETSc could not build an auxiliary field"
        call PetscFEDestroy(fe, ierr)
    end subroutine add_aux_field

    !> Write each cell's closure slot and saturated inverse mobility.
    subroutine store_materials(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: cell, c_start, c_end, offset, material
        PetscErrorCode :: ierr
        real(real64) :: inverse_mobility(PLANE_DIMENSION*PLANE_DIMENSION), conductivity(PLANE_DIMENSION, PLANE_DIMENSION)
        integer(int32) :: slot, k
        logical :: ok

        call DMGetLocalSection(self%aux_dm, section, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(self%aux_dm, 0, c_start, c_end, ierr)
        if (ierr == 0) call VecGetArray(self%auxiliary, raw, ierr)
        if (ierr /= 0) error stop "PETSc could not open the auxiliary vector"

        do cell = c_start, c_end - 1
            call DMGetLabelValue(self%aux_dm, "Cell Sets", cell, material, ierr)
            slot = findloc(self%material_ids, material, dim=1)
            if (slot == 0 .and. size(self%material_ids) == 1) slot = 1
            if (slot == 0) error stop "A mesh cell carries a region no material of Basic.json defines"

            conductivity = 0.0d0
            do k = 1, PLANE_DIMENSION
                conductivity(k, k) = closure_saturated_conductivity(slot)
            end do
            ! A run without gravity still needs the mobility's units.
            call invert_mobility(conductivity, self%density, gravity_acceleration, inverse_mobility, ok)
            if (.not. ok) error stop "An impermeable region is a different model, not a small conductivity"

            call PetscSectionGetFieldOffset(section, cell, AUX_INVERSE_MOBILITY - 1, offset, ierr)
            raw(offset + 1:offset + PLANE_DIMENSION**2) = inverse_mobility
            call PetscSectionGetFieldOffset(section, cell, AUX_MATERIAL - 1, offset, ierr)
            raw(offset + 1) = real(slot, real64)
            call PetscSectionGetFieldOffset(section, cell, AUX_SOURCE - 1, offset, ierr)
            raw(offset + 1) = 0.0d0
        end do
        call VecRestoreArray(self%auxiliary, raw, ierr)
    end subroutine store_materials

    !> Set every degree of freedom of auxiliary field `field` to `value`.
    subroutine fill_aux_field(dm, vec, field, value)
        implicit none
        DM, intent(in) :: dm
        Vec, intent(inout) :: vec
        integer(int32), intent(in) :: field
        real(real64), intent(in) :: value

        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: point, p_start, p_end, offset, components
        PetscErrorCode :: ierr

        call DMGetLocalSection(dm, section, ierr)
        if (ierr == 0) call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr == 0) call VecGetArray(vec, raw, ierr)
        if (ierr /= 0) error stop "PETSc could not open an auxiliary vector"
        do point = p_start, p_end - 1
            call PetscSectionGetFieldDof(section, point, field - 1, components, ierr)
            if (components <= 0) cycle
            call PetscSectionGetFieldOffset(section, point, field - 1, offset, ierr)
            raw(offset + 1:offset + components) = value
        end do
        call VecRestoreArray(vec, raw, ierr)
    end subroutine fill_aux_field

    !> Register every condition of `Conditions.json` in the role its space gives it.
    !>
    !> Fails by `error stop` naming the condition when its role has no
    !> implementation in that space: a state-dependent flux on a mixed field,
    !> and seepage.
    subroutine add_boundaries(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        type(type_config_bc), allocatable :: configs(:)
        integer(int32) :: id, condition, slot, boundary_id
        logical :: mixed

        call clear_boundary_registry()
        do id = 1, PHYSICS_TYPES%NUM_ID
            if (.not. self%is_active(id)) cycle
            call input_translator%execute(self%input, PHYSICS_TYPES%to_object(id), configs)
            if (.not. allocated(configs)) allocate (configs(0))
            call install_boundary_conditions(id, configs)
            call copy_boundary_label(self, id, configs)
            mixed = self%space(id) == FE_SPACES%MIXED%ID
            if (mixed) call close_unlisted_boundary(self, configs)

            do condition = 1, size(configs)
                boundary_id = configs(condition)%boundary_id
                call report_boundary(self, id, boundary_id)
                select case (role_of(configs(condition)%bc_kind%ID))
                case (CONDITION_POTENTIAL)
                    if (mixed) then
                        slot = add_slot(self, id, condition, ROLE_POTENTIAL)
                        call add_natural(self, id, boundary_id, field_of(self, id, ROLE_FLUX), slot, .false.)
                    else
                        slot = add_slot(self, id, condition, ROLE_POTENTIAL)
                        call add_essential_potential(self, id, boundary_id, field_of(self, id, ROLE_POTENTIAL), slot)
                    end if
                case (CONDITION_FLUX)
                    slot = add_slot(self, id, condition, merge(ROLE_FLUX, ROLE_POTENTIAL, mixed))
                    if (mixed) then
                        if (.not. is_constant_flux(configs(condition)%bc_kind%ID)) call refuse(id, boundary_id)
                        call add_essential_flux(self, id, boundary_id, field_of(self, id, ROLE_FLUX), slot)
                    else
                        call add_natural(self, id, boundary_id, field_of(self, id, ROLE_POTENTIAL), slot, .true.)
                    end if
                case (CONDITION_SWITCHING)
                    ! A seepage face of a Lagrange pressure: its flux enters naturally
                    ! and the variable bounds of the solve switch it to a pressure.
                    if (mixed .or. id /= PHYSICS_TYPES%HYDRAULIC%ID) call refuse(id, boundary_id)
                    slot = add_slot(self, id, condition, ROLE_POTENTIAL)
                    call add_seepage(self, id, boundary_id, field_of(self, id, ROLE_POTENTIAL), slot)
                case default
                    call refuse(id, boundary_id)
                end select
            end do
            deallocate (configs)
        end do
    contains
        subroutine refuse(physics_id, boundary)
            integer(int32), intent(in) :: physics_id, boundary
            write (*, '(a,i0)') "Not implemented: this "//PHYSICS_TYPES%to_name(physics_id)// &
                " condition in its space, on boundary ", boundary
            error stop "A boundary condition has no implementation in the space of its equation"
        end subroutine refuse
    end subroutine add_boundaries

    !> Log how many facets of the whole mesh carry a condition's boundary id.
    !>
    !> A condition on an id the mesh does not carry would otherwise do nothing
    !> without saying so.
    subroutine report_boundary(self, id, boundary_id)
        implicit none
        class(type_ftcms), intent(in) :: self
        integer(int32), intent(in) :: id, boundary_id

        PetscInt :: local, total
        PetscErrorCode :: ierr

        call DMGetStratumSize(self%dm, "Face Sets", boundary_id, local, ierr)
        if (ierr /= 0) error stop "PETSc could not count a boundary"
        call MPI_Allreduce(local, total, 1, MPIU_INTEGER, MPI_SUM, PETSC_COMM_WORLD, ierr)
        call global_logger%log_information(trim(PHYSICS_TYPES%to_name(id))//" condition on boundary "//to_text(boundary_id)// &
                                           ": "//to_text(total)//" facets")
        if (total == 0) error stop "A condition of Conditions.json names a boundary the mesh does not carry"
    end subroutine report_boundary

    !> Whether a condition prescribes the potential, an inflow, or switches between them.
    pure integer(int32) function role_of(kind)
        implicit none
        integer(int32), intent(in) :: kind

        if (kind == THERMAL_BC_TYPES%DIRICHLET%ID .or. kind == HYDRAULIC_BC_TYPES%DIRICHLET%ID) then
            role_of = CONDITION_POTENTIAL
        else if (kind == HYDRAULIC_BC_TYPES%SEEPAGE%ID) then
            role_of = CONDITION_SWITCHING
        else
            role_of = CONDITION_FLUX
        end if
    end function role_of

    !> Whether an inflow condition does not depend on the state.
    pure logical function is_constant_flux(kind)
        implicit none
        integer(int32), intent(in) :: kind

        is_constant_flux = (kind == THERMAL_BC_TYPES%FLUX%ID .or. kind == HYDRAULIC_BC_TYPES%FLUX%ID)
    end function is_constant_flux

    !> The zero-based number of the table field of equation `id` and `role`.
    function field_of(self, id, role) result(number)
        implicit none
        class(type_ftcms), intent(in) :: self
        integer(int32), intent(in) :: id, role
        integer(int32) :: number

        integer(int32) :: f

        do f = 1, self%num_fields
            if (self%fields(f)%physics_id == id .and. self%fields(f)%role == role) then
                number = self%fields(f)%number
                return
            end if
        end do
        error stop "The field table has no such field"
    end function field_of

    !> Register a condition's slot with the scales of the field it reads and the row it enters.
    function add_slot(self, id, condition, role) result(slot)
        implicit none
        class(type_ftcms), intent(in) :: self
        integer(int32), intent(in) :: id, condition, role
        integer(int32) :: slot

        real(real64) :: unknown_scale, row_scale

        if (id == PHYSICS_TYPES%THERMAL%ID) then
            unknown_scale = self%temperature_scale
            row_scale = 1.0d0/self%thermal_scale
        else if (role == ROLE_FLUX) then
            unknown_scale = self%flux_scale
            row_scale = 1.0d0
        else
            unknown_scale = self%pressure_scale
            row_scale = self%length_scale/self%flux_scale
        end if
        slot = add_boundary_slot(id, condition, field_of(self, id, role), unknown_scale, row_scale)
    end function add_slot

    subroutine add_essential_potential(self, id, boundary_id, field, slot)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: id, boundary_id, field, slot

        integer(c_int) :: index, cerr

        call petscds_add_boundary_essential(transfer(self%dm%v, c_null_ptr), boundary_name(slot), boundary_label(id), &
                                            int(boundary_id, c_int), int(field, c_int), &
                                            c_funloc(bc_essential_potential), boundary_slot_context(slot), &
                                            index, cerr)
        if (cerr /= 0) error stop "PETSc rejected an essential boundary condition"
    end subroutine add_essential_potential

    subroutine add_essential_flux(self, id, boundary_id, field, slot)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: id, boundary_id, field, slot

        integer(c_int) :: index, cerr

        call record_facet_normals(self, boundary_id, slot)
        call petscds_add_boundary_essential(transfer(self%dm%v, c_null_ptr), boundary_name(slot), boundary_label(id), &
                                            int(boundary_id, c_int), int(field, c_int), &
                                            c_funloc(bc_essential_normal_flux), boundary_slot_context(slot), &
                                            index, cerr)
        if (cerr /= 0) error stop "PETSc rejected an essential flux condition"
    end subroutine add_essential_flux

    !> Hand the registry the end points and outward unit normals of a boundary's local facets.
    !>
    !> The outward side is the one away from the facet's single cell.
    subroutine record_facet_normals(self, boundary_id, slot)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: boundary_id, slot

        IS :: facets
        Vec :: coordinates
        PetscSection :: coordinate_section
        PetscScalar, pointer :: xyz(:)
        PetscInt, pointer :: points(:), cone(:), support(:)
        PetscInt :: k, offset
        PetscReal :: volume, centroid(3), normal(3)
        PetscErrorCode :: ierr
        real(real64), allocatable :: starts(:, :), ends(:, :), normals(:, :)
        real(real64) :: tangent(PLANE_DIMENSION), outward(PLANE_DIMENSION)
        integer(int32) :: count

        allocate (starts(PLANE_DIMENSION, 0), ends(PLANE_DIMENSION, 0), normals(PLANE_DIMENSION, 0))
        if (.not. stratum_is_empty(self%dm, boundary_id)) then
            call DMGetStratumIS(self%dm, "Face Sets", boundary_id, facets, ierr)
            if (ierr /= 0) error stop "PETSc could not read a boundary"
            call DMGetCoordinatesLocal(self%dm, coordinates, ierr)
            call DMGetCoordinateSection(self%dm, coordinate_section, ierr)
            call VecGetArrayRead(coordinates, xyz, ierr)
            call ISGetIndices(facets, points, ierr)
            deallocate (starts, ends, normals)
            allocate (starts(PLANE_DIMENSION, size(points)), ends(PLANE_DIMENSION, size(points)), &
                      normals(PLANE_DIMENSION, size(points)))
            count = 0
            do k = 1, size(points)
                call DMPlexGetSupport(self%dm, points(k), support, ierr)
                if (size(support) /= 1) cycle
                count = count + 1
                call DMPlexGetCone(self%dm, points(k), cone, ierr)
                call PetscSectionGetOffset(coordinate_section, cone(1), offset, ierr)
                starts(:, count) = xyz(offset + 1:offset + PLANE_DIMENSION)
                call PetscSectionGetOffset(coordinate_section, cone(2), offset, ierr)
                ends(:, count) = xyz(offset + 1:offset + PLANE_DIMENSION)
                call DMPlexComputeCellGeometryFVM(self%dm, support(1), volume, centroid, normal, ierr)
                call DMPlexRestoreCone(self%dm, points(k), cone, ierr)
                call DMPlexRestoreSupport(self%dm, points(k), support, ierr)
                tangent = ends(:, count) - starts(:, count)
                outward = [tangent(2), -tangent(1)]/norm2(tangent)
                if (dot_product(outward, starts(:, count) - centroid(1:PLANE_DIMENSION)) < 0.0d0) then
                    outward = -outward
                end if
                normals(:, count) = outward
            end do
            call ISRestoreIndices(facets, points, ierr)
            call ISDestroy(facets, ierr)
            call VecRestoreArrayRead(coordinates, xyz, ierr)
            starts = starts(:, :count)
            ends = ends(:, :count)
            normals = normals(:, :count)
        end if
        call set_slot_facets(slot, starts, ends, normals)
    end subroutine record_facet_normals

    !> A seepage face: its constant potential flux, and a place in the bounds.
    subroutine add_seepage(self, id, boundary_id, field, slot)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: id, boundary_id, field, slot

        integer(c_int) :: index, cerr

        call petscds_add_boundary_natural(transfer(self%dm%v, c_null_ptr), boundary_name(slot), boundary_label(id), &
                                          int(boundary_id, c_int), int(field, c_int), c_funloc(bd_f0_seepage), &
                                          c_null_funptr, index, cerr)
        if (cerr /= 0) error stop "PETSc rejected a seepage face"
        call attach_boundary_auxiliary(self, id, boundary_id, slot)
        self%num_seepage = self%num_seepage + 1
        self%seepage_slots(self%num_seepage) = slot
        self%seepage_ids(self%num_seepage) = boundary_id
    end subroutine add_seepage

    !> A natural condition and the auxiliary vector that tells its kernels the slot.
    subroutine add_natural(self, id, boundary_id, field, slot, potential_row)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: id, boundary_id, field, slot
        !> Whether it enters a potential row (an inflow) or a flux row (a potential).
        logical, intent(in) :: potential_row

        integer(c_int) :: index, cerr

        if (potential_row) then
            call petscds_add_boundary_natural(transfer(self%dm%v, c_null_ptr), boundary_name(slot), boundary_label(id), &
                                              int(boundary_id, c_int), int(field, c_int), &
                                              c_funloc(bd_f0_natural), c_funloc(bd_g0_natural), index, cerr)
        else
            call petscds_add_boundary_natural(transfer(self%dm%v, c_null_ptr), boundary_name(slot), boundary_label(id), &
                                              int(boundary_id, c_int), int(field, c_int), &
                                              c_funloc(bd_f0_flux_potential), c_null_funptr, index, cerr)
        end if
        if (cerr /= 0) error stop "PETSc rejected a natural boundary condition"
        call attach_boundary_auxiliary(self, id, boundary_id, slot)
    end subroutine add_natural

    !> Clone the auxiliary vector for one boundary value and mark its slot.
    subroutine attach_boundary_auxiliary(self, id, boundary_id, slot)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: id, boundary_id, slot

        DMLabel :: label
        PetscErrorCode :: ierr

        self%num_boundary_auxiliary = self%num_boundary_auxiliary + 1
        associate (vec => self%boundary_auxiliary(self%num_boundary_auxiliary))
            call DMCreateLocalVector(self%aux_dm, vec, ierr)
            if (ierr == 0) call VecCopy(self%auxiliary, vec, ierr)
            if (ierr /= 0) error stop "PETSc could not build a boundary auxiliary vector"
            call fill_aux_field(self%aux_dm, vec, AUX_BOUNDARY_SLOT, real(slot, real64))
            call DMGetLabel(self%dm, boundary_label(id), label, ierr)
            if (ierr == 0) call DMSetAuxiliaryVec(self%dm, label, boundary_id, &
                                                  0, vec, ierr)
            if (ierr /= 0) error stop "PETSc could not attach a boundary auxiliary vector"
        end associate
    end subroutine attach_boundary_auxiliary

    !> Close to flow every boundary facet no condition of a mixed equation names.
    !>
    !> In a mixed space the natural condition is a zero potential, not a zero
    !> flux, so the boundary a project leaves unstated would otherwise drain.
    !> This keeps an unstated boundary closed, as it is for a `lagrange` field.
    subroutine close_unlisted_boundary(self, configs)
        implicit none
        class(type_ftcms), intent(inout) :: self
        type(type_config_bc), intent(in) :: configs(:)

        DMLabel :: closed
        IS :: listed
        PetscInt, pointer :: points(:)
        PetscInt :: k
        PetscErrorCode :: ierr
        integer(int32) :: condition
        integer(c_int) :: index, cerr

        call DMCreateLabel(self%dm, CLOSED_LABEL, ierr)
        if (ierr == 0) call DMGetLabel(self%dm, CLOSED_LABEL, closed, ierr)
        if (ierr == 0) call DMPlexMarkBoundaryFaces(self%dm, 1, closed, ierr)
        if (ierr /= 0) error stop "PETSc could not mark the boundary"
        do condition = 1, size(configs)
            if (stratum_is_empty(self%dm, configs(condition)%boundary_id)) cycle
            call DMGetStratumIS(self%dm, "Face Sets", configs(condition)%boundary_id, listed, ierr)
            if (ierr /= 0) error stop "PETSc could not read a boundary"
            call ISGetIndices(listed, points, ierr)
            do k = 1, size(points)
                call DMLabelClearValue(closed, points(k), 1, ierr)
            end do
            call ISRestoreIndices(listed, points, ierr)
            call ISDestroy(listed, ierr)
        end do
        call petscds_add_boundary_essential(transfer(self%dm%v, c_null_ptr), "closed", CLOSED_LABEL, 1_c_int, &
                                            int(self%flux_field, c_int), c_funloc(bc_essential_zero), c_null_ptr, &
                                            index, cerr)
        if (cerr /= 0) error stop "PETSc rejected the closed boundary"
    end subroutine close_unlisted_boundary

    !> The boundary label of one equation: a copy of "Face Sets" for its conditions.
    !>
    !> PETSc keys a boundary's auxiliary vector by (label, value) alone, so two
    !> equations with conditions on the same boundary id need two labels for
    !> their kernels to find their own slots.
    pure function boundary_label(id) result(name)
        implicit none
        integer(int32), intent(in) :: id
        character(len=:), allocatable :: name

        name = trim(PHYSICS_TYPES%to_name(id))//"_boundaries"
    end function boundary_label

    !> Copy the "Face Sets" points of every listed boundary id into the equation's label.
    subroutine copy_boundary_label(self, id, configs)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: id
        type(type_config_bc), intent(in) :: configs(:)

        IS :: points
        PetscInt, pointer :: indices(:)
        PetscInt :: k, value
        PetscErrorCode :: ierr
        integer(int32) :: condition

        call DMCreateLabel(self%dm, boundary_label(id), ierr)
        if (ierr /= 0) error stop "PETSc could not create a boundary label"
        do condition = 1, size(configs)
            value = configs(condition)%boundary_id
            if (stratum_is_empty(self%dm, configs(condition)%boundary_id)) cycle
            call DMGetStratumIS(self%dm, "Face Sets", value, points, ierr)
            if (ierr /= 0) error stop "PETSc could not read a boundary"
            call ISGetIndices(points, indices, ierr)
            do k = 1, size(indices)
                call DMSetLabelValue(self%dm, boundary_label(id), indices(k), value, ierr)
            end do
            call ISRestoreIndices(points, indices, ierr)
            call ISDestroy(points, ierr)
        end do
    end subroutine copy_boundary_label

    pure function boundary_name(slot) result(name)
        implicit none
        integer(int32), intent(in) :: slot
        character(len=24) :: name

        write (name, '(a,i0)') "condition-", slot
    end function boundary_name

    !> The uniform initial state, published as the accepted one.
    subroutine set_initial_state(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: point, p_start, p_end, offset, num_dofs
        PetscErrorCode :: ierr
        Vec :: local

        call DMCreateGlobalVector(self%dm, self%solution, ierr)
        if (ierr == 0) call PetscObjectSetName(PetscObjectCast(self%solution), "solution", ierr)
        if (ierr == 0) call DMCreateGlobalVector(self%dm, self%old_solution, ierr)
        if (ierr == 0) call DMCreateLocalVector(self%dm, self%local_solution, ierr)
        if (ierr == 0) call DMGetLocalVector(self%dm, local, ierr)
        if (ierr == 0) call VecZeroEntries(local, ierr)
        if (ierr == 0) call DMGetLocalSection(self%dm, section, ierr)
        if (ierr == 0) call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr == 0) call VecGetArray(local, raw, ierr)
        if (ierr /= 0) error stop "PETSc could not create the solution"

        do point = p_start, p_end - 1
            if (self%temperature_field >= 0) then
                call PetscSectionGetFieldDof(section, point, self%temperature_field, num_dofs, ierr)
                if (num_dofs > 0) then
                    call PetscSectionGetFieldOffset(section, point, self%temperature_field, offset, ierr)
                    raw(offset + 1:offset + num_dofs) = self%initial_temperature/self%temperature_scale
                end if
            end if
            if (self%pressure_field >= 0) then
                call PetscSectionGetFieldDof(section, point, self%pressure_field, num_dofs, ierr)
                if (num_dofs > 0) then
                    call PetscSectionGetFieldOffset(section, point, self%pressure_field, offset, ierr)
                    raw(offset + 1:offset + num_dofs) = self%initial_pressure/self%pressure_scale
                end if
            end if
        end do
        call VecRestoreArray(local, raw, ierr)

        call DMLocalToGlobalBegin(self%dm, local, INSERT_VALUES, self%solution, ierr)
        if (ierr == 0) call DMLocalToGlobalEnd(self%dm, local, INSERT_VALUES, self%solution, ierr)
        call DMRestoreLocalVector(self%dm, local, ierr)
        if (ierr == 0) call VecCopy(self%solution, self%old_solution, ierr)
        if (ierr /= 0) error stop "PETSc could not set the initial state"

        self%time = self%start_time
        call set_boundary_time(self%time)
        call publish_accepted_state(self, ierr)
        ! The initial state is its own history.
        if (ierr == 0) call publish_accepted_state(self, ierr)
        if (ierr == 0) call VecCopy(self%auxiliary, self%accepted_auxiliary, ierr)
        if (ierr == 0) call VecDuplicate(self%solution, self%older_solution, ierr)
        if (ierr == 0) call VecDuplicate(self%solution, self%oldest_solution, ierr)
        if (ierr == 0) call VecCopy(self%solution, self%older_solution, ierr)
        if (ierr == 0) call VecCopy(self%solution, self%oldest_solution, ierr)
        if (ierr /= 0) error stop "PETSc could not publish the initial state"
        call build_error_weights(self)
    end subroutine set_initial_state

    !> The absolute tolerance of each potential's dofs, in units of the unknown.
    !>
    !> Flux dofs carry a zero mask: the error control measures the potentials,
    !> as `Conditions.json` states its tolerances for temperature and pressure.
    subroutine build_error_weights(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscSection :: section
        PetscScalar, pointer :: tolerance_raw(:), mask_raw(:)
        PetscInt :: point, p_start, p_end, offset, num_dofs
        PetscErrorCode :: ierr
        Vec :: tolerance_local, mask_local

        call VecDuplicate(self%solution, self%error_tolerance, ierr)
        if (ierr == 0) call VecDuplicate(self%solution, self%error_mask, ierr)
        if (ierr == 0) call DMGetLocalVector(self%dm, tolerance_local, ierr)
        if (ierr == 0) call DMGetLocalVector(self%dm, mask_local, ierr)
        if (ierr == 0) call VecSet(tolerance_local, 1.0d0, ierr)
        if (ierr == 0) call VecSet(mask_local, 0.0d0, ierr)
        if (ierr == 0) call DMGetLocalSection(self%dm, section, ierr)
        if (ierr == 0) call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr == 0) call VecGetArray(tolerance_local, tolerance_raw, ierr)
        if (ierr == 0) call VecGetArray(mask_local, mask_raw, ierr)
        if (ierr /= 0) error stop "PETSc could not build the error weights"
        do point = p_start, p_end - 1
            if (self%temperature_field >= 0) then
                call PetscSectionGetFieldDof(section, point, self%temperature_field, num_dofs, ierr)
                if (num_dofs > 0) then
                    call PetscSectionGetFieldOffset(section, point, self%temperature_field, offset, ierr)
                    tolerance_raw(offset + 1:offset + num_dofs) = self%error_control%temperature_tolerance &
                                                                  /self%temperature_scale
                    mask_raw(offset + 1:offset + num_dofs) = 1.0d0
                end if
            end if
            if (self%pressure_field >= 0) then
                call PetscSectionGetFieldDof(section, point, self%pressure_field, num_dofs, ierr)
                if (num_dofs > 0) then
                    call PetscSectionGetFieldOffset(section, point, self%pressure_field, offset, ierr)
                    tolerance_raw(offset + 1:offset + num_dofs) = self%error_control%pressure_tolerance &
                                                                  /self%pressure_scale
                    mask_raw(offset + 1:offset + num_dofs) = 1.0d0
                end if
            end if
        end do
        call VecRestoreArray(tolerance_local, tolerance_raw, ierr)
        call VecRestoreArray(mask_local, mask_raw, ierr)
        call DMLocalToGlobalBegin(self%dm, tolerance_local, INSERT_VALUES, self%error_tolerance, ierr)
        call DMLocalToGlobalEnd(self%dm, tolerance_local, INSERT_VALUES, self%error_tolerance, ierr)
        call DMLocalToGlobalBegin(self%dm, mask_local, INSERT_VALUES, self%error_mask, ierr)
        call DMLocalToGlobalEnd(self%dm, mask_local, INSERT_VALUES, self%error_mask, ierr)
        call DMRestoreLocalVector(self%dm, tolerance_local, ierr)
        call DMRestoreLocalVector(self%dm, mask_local, ierr)
        if (ierr /= 0) error stop "PETSc could not build the error weights"
    end subroutine build_error_weights

    !> The nonlinear solver, configured from `Basic.json` under one prefix.
    !>
    !> JSON options are inserted under `SOLVER_PREFIX` only where the command
    !> line did not name them, so a command line overrides the project. The
    !> defaults are a direct solve with MUMPS and Deuflhard's `nleqerr` line
    !> search, which measures a trial step in the norm of the Newton correction:
    !> across a freezing front the latent rows dwarf their neighbours, and a
    !> residual-norm search shrinks steps that reduce the error.
    subroutine build_solver(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscErrorCode :: ierr
        character(len=16) :: count

        associate (solver => self%input%basic%solver_settings)
            call insert_prefixed_options(trim(solver%linear_solver%petsc_options))
            write (count, '(i0)') solver%nonlinear_solver%max_iterations
            if (solver%nonlinear_solver%max_iterations > 0) call set_default_option("-"//SOLVER_PREFIX//"snes_max_it", &
                                                                                    trim(count))
        end associate
        if (self%num_seepage > 0) then
            ! Seepage faces are bounds on the pressure: a variational inequality.
            call set_default_option("-"//SOLVER_PREFIX//"snes_type", "vinewtonrsls")
        else
            call set_default_option("-"//SOLVER_PREFIX//"snes_linesearch_type", "nleqerr")
        end if
        if (.not. option_is_set("-"//SOLVER_PREFIX//"pc_type")) then
            call set_default_option("-"//SOLVER_PREFIX//"ksp_type", "preonly")
            call set_default_option("-"//SOLVER_PREFIX//"pc_type", "lu")
            call set_default_option("-"//SOLVER_PREFIX//"pc_factor_mat_solver_type", "mumps")
        end if

        call SNESCreate(PETSC_COMM_WORLD, self%snes, ierr)
        if (ierr == 0) call SNESSetOptionsPrefix(self%snes, SOLVER_PREFIX, ierr)
        if (ierr == 0) call SNESSetDM(self%snes, self%dm, ierr)
        if (ierr == 0) call DMPlexSetSNESLocalFEM(self%dm, PETSC_FALSE, PETSC_NULL_INTEGER, ierr)
        if (ierr == 0) call SNESSetFromOptions(self%snes, ierr)
        if (ierr /= 0) error stop "PETSc could not set up the nonlinear solver"
        if (self%num_seepage > 0) then
            call VecDuplicate(self%solution, self%lower_bound, ierr)
            if (ierr == 0) call VecDuplicate(self%solution, self%upper_bound, ierr)
            if (ierr /= 0) error stop "PETSc could not create the variable bounds"
        end if

        if (self%nodal_flux) then
            call install_nodal_flux(self%dm, self%aux_dm, self%auxiliary, self%pressure_scale, &
                                    self%temperature_scale, self%flux_scale, self%length_scale, &
                                    self%density, self%gravity, ierr)
            if (ierr /= 0) error stop "PETSc rejected the nodal flux"
        end if
    end subroutine build_solver

    !> Insert `-key value` pairs under the solver prefix, where not already set.
    subroutine insert_prefixed_options(options)
        implicit none
        character(len=*), intent(in) :: options

        integer(int32) :: start, finish, next
        character(len=:), allocatable :: key, value

        start = 1
        do
            call next_token(options, start, finish)
            if (finish < start) exit
            key = options(start:finish)
            next = finish + 1
            if (key(1:1) /= "-") error stop "A PETSc option of Basic.json does not start with a dash"
            value = ""
            call next_token(options, next, finish)
            if (finish >= next) then
                if (options(next:next) /= "-") then
                    value = options(next:finish)
                    next = finish + 1
                end if
            end if
            call set_default_option("-"//SOLVER_PREFIX//key(2:), value)
            start = next
        end do
    end subroutine insert_prefixed_options

    !> The next blank-separated token of `text` from `start`: text(start:finish).
    pure subroutine next_token(text, start, finish)
        implicit none
        character(len=*), intent(in) :: text
        integer(int32), intent(inout) :: start
        integer(int32), intent(inout) :: finish

        do while (start <= len(text))
            if (text(start:start) /= " ") exit
            start = start + 1
        end do
        finish = start - 1
        do while (finish < len(text))
            if (text(finish + 1:finish + 1) == " ") exit
            finish = finish + 1
        end do
    end subroutine next_token

    logical function option_is_set(key)
        implicit none
        character(len=*), intent(in) :: key

        PetscBool :: found
        PetscErrorCode :: ierr

        call PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, key, found, ierr)
        option_is_set = (ierr == 0 .and. (found .eqv. PETSC_TRUE))
    end function option_is_set

    !> Set one option unless the run already carries it.
    subroutine set_default_option(key, value)
        implicit none
        character(len=*), intent(in) :: key
        !> The value; empty for a flag, which PETSc then reads as given without one.
        character(len=*), intent(in) :: value

        PetscErrorCode :: ierr

        if (option_is_set(key)) return
        if (len(value) == 0) then
            ! A flag such as -snes_test_jacobian takes an optional number, so "true" would be read as a bad one.
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS, key, PETSC_NULL_CHARACTER, ierr)
        else
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS, key, value, ierr)
        end if
        if (ierr /= 0) error stop "PETSc rejected an option"
    end subroutine set_default_option

    !> The time window and the step controller of `Conditions.json`.
    subroutine read_time_controls(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        type(type_config_time) :: time
        type(type_config_time_ats) :: adaptive

        call input_translator%execute(self%input, time)
        call input_translator%execute(self%input, adaptive)
        self%start_time = time%start_time
        self%end_time = time%end_time
        self%target_order = min(max(time%target_bdf_order, 1), 2)

        associate (error => self%error_control)
            error%active = adaptive%active .and. adaptive%use_error_control
            error%relative_tolerance = adaptive%error_rtol
            error%temperature_tolerance = adaptive%error_atol_temperature
            error%pressure_tolerance = adaptive%error_atol_pressure
            error%proportional_gain = adaptive%pi_k_p
            error%integral_gain = adaptive%pi_k_i
            error%safety = adaptive%safety_factor
            error%max_growth = adaptive%max_growth_rate
        end associate

        associate (control => self%step_control)
            control%initial = time%initial_step
            if (adaptive%active) then
                control%minimum = adaptive%dt_min
                control%maximum = adaptive%dt_max
                control%iterations_low = adaptive%iter_min
                control%iterations_high = adaptive%iter_max
                control%grow = adaptive%scale_up
                control%shrink = adaptive%scale_down
                control%retry = adaptive%scale_retry
            else
                control%minimum = time%initial_step
                control%maximum = time%initial_step
                control%grow = 1.0d0
                control%shrink = 1.0d0
            end if
            if (.not. (control%initial > 0.0d0) .or. .not. (control%minimum > 0.0d0) &
                .or. control%maximum < control%minimum) then
                error stop "The time steps of Conditions.json need positive, ordered bounds"
            end if
            if (control%grow < 1.0d0 .or. control%shrink > 1.0d0 .or. .not. (control%shrink > 0.0d0) &
                .or. .not. (control%retry > 0.0d0 .and. control%retry < 1.0d0)) then
                error stop "The step factors of Conditions.json need growth at least one and reductions below it"
            end if
            if (control%iterations_low > control%iterations_high) then
                error stop "The iteration window of Conditions.json is not ordered"
            end if
            self%next_step = min(max(control%initial, control%minimum), control%maximum)
        end associate
        self%dt = self%next_step
    end subroutine read_time_controls

    pure function to_text(value) result(text)
        implicit none
        integer(int32), intent(in) :: value
        character(len=:), allocatable :: text

        character(len=16) :: buffer

        write (buffer, '(i0)') value
        text = trim(buffer)
    end function to_text

    module subroutine destroy_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscErrorCode :: ierr
        integer(int32) :: b

        if (.not. self%ready) return
        if (self%assimilation_enabled) call self%assimilation%destroy()
        call destroy_output(self)
        if (self%nodal_flux) call clear_nodal_flux()
        call SNESDestroy(self%snes, ierr)
        do b = 1, self%num_boundary_auxiliary
            call VecDestroy(self%boundary_auxiliary(b), ierr)
        end do
        call VecDestroy(self%accepted_auxiliary, ierr)
        call VecDestroy(self%auxiliary, ierr)
        call VecDestroy(self%local_solution, ierr)
        call VecDestroy(self%old_solution, ierr)
        call VecDestroy(self%older_solution, ierr)
        call VecDestroy(self%oldest_solution, ierr)
        call VecDestroy(self%error_tolerance, ierr)
        call VecDestroy(self%error_mask, ierr)
        if (self%num_seepage > 0) then
            call VecDestroy(self%lower_bound, ierr)
            call VecDestroy(self%upper_bound, ierr)
        end if
        call VecDestroy(self%solution, ierr)
        call DMDestroy(self%aux_dm, ierr)
        call DMDestroy(self%dm, ierr)
        call clear_boundary_registry()
        call clear_closure()
        call self%input%geometry%mesh%destroy()
        self%ready = .false.
    end subroutine destroy_ftcms

end submodule ftcms_setup
