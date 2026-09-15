#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscsnes.h>

!> Construction of the RT0/P0 discrete system, its coefficients and boundaries.
!>
!> The flux space is built from PETSc's trimmed polynomial spaces rather than
!> from a named Raviart-Thomas type, which this PETSc does not have: a trimmed
!> primal space of form degree minus one with a trimmed, continuous Lagrange
!> dual space of order one is the lowest-order H(div) element on a simplex.
!> The options are set before the element is created because
!> PetscFECreateDefault configures each sub-object from the options database.
!>
!> The material response is carried per cell in auxiliary fields rather than in
!> DS constants, so a mesh of several materials needs one set of kernels. The
!> inverse mobility is inverted once per cell here, from that cell's own
!> conductivity; nothing averages two materials together and nothing inverts an
!> already averaged number.
!>
!> A prescribed-pressure boundary is registered as a natural condition with its
!> own auxiliary vector, because PETSc keys a boundary's auxiliary data by the
!> label value it was registered with. A prescribed-flux boundary is essential:
!> it constrains the RT facet moments, which is the only place a normal flux can
!> be stated in this formulation, and a zero flux is such a constraint too.
submodule(app_mixed_hydraulic) mixed_hydraulic_setup
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_funptr, c_funloc, c_null_funptr, c_null_ptr
    use :: core_interop_petscds, only:petscds_set_residual, petscds_set_jacobian, &
                                      petscds_add_boundary_essential, petscds_add_boundary_natural
    use :: app_mixed_hydraulic_physics, only:invert_mobility, &
                                             MIXED_AUX_FIELDS, &
                                             A_INVERSE_MOBILITY, A_STORAGE, A_SOURCE, &
                                             A_PRESSURE_OLD, A_BC_PRESSURE, &
                                             C_DT, C_FLUX_SCALE, C_PRESSURE_SCALE, C_LENGTH_SCALE, &
                                             C_DENSITY, C_GRAVITY_X, C_GRAVITY_Z
    use :: app_mixed_hydraulic_pointwise, only:f0_flux, f1_flux, f0_conservation, &
                                               g0_flux_flux, g2_flux_pressure, &
                                               g1_pressure_flux, g0_pressure_pressure, &
                                               bd_f0_flux_pressure, bc_essential_flux, &
                                               set_flux_boundary_value, flux_boundary_context, &
                                               reset_flux_boundaries, MAX_FLUX_BOUNDARIES
    implicit none

    ! The RT mass integral is of degree two on an affine cell, which is what a
    ! constant coefficient needs. A variable coefficient raises this
    ! independently of the element degree.
    PetscInt, parameter :: QUADRATURE_ORDER = 2_PETSC_INT_KIND

contains

    module subroutine initialize_mixed_hydraulic(self, dm, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout), target :: self
        DM, intent(in) :: dm
        PetscErrorCode, intent(inout) :: ierr

        PetscFE :: fe_flux, fe_pressure
        PetscObject :: discretisation
        PetscDS :: ds
        PetscInt :: dimension, coordinate_dimension
        PetscBool :: is_simplex
        type(c_ptr) :: dm_handle, ds_handle
        integer(c_int) :: cerr

        if (.not. allocated(self%materials)) error stop "The mixed hydraulic solver has no materials"

        call DMGetDimension(dm, dimension, ierr)
        if (ierr /= 0) error stop "The mesh has no dimension"
        call DMGetCoordinateDim(dm, coordinate_dimension, ierr)
        if (ierr /= 0) error stop "The mesh has no coordinate dimension"
        if (int(dimension, int32) /= MIXED_DIMENSION .or. dimension /= coordinate_dimension) then
            error stop "The mixed hydraulic solver needs an intrinsic two-dimensional mesh: "// &
                "map a planar (x,0,z) mesh to (x,z) before handing it over"
        end if
        call DMPlexIsSimplex(dm, is_simplex, ierr)
        if (ierr /= 0 .or. is_simplex .neqv. PETSC_TRUE) then
            error stop "The mixed hydraulic solver supports simplicial cells only"
        end if

        call DMClone(dm, self%dm, ierr)
        if (ierr /= 0) error stop "PETSc could not clone the mesh for the mixed path"

        call set_mixed_element_options(ierr)
        if (ierr /= 0) error stop "PETSc rejected the RT0/P0 element options"

        call PetscFECreateDefault(PETSC_COMM_WORLD, dimension, 2_PETSC_INT_KIND, PETSC_TRUE, &
                                  "flux_", QUADRATURE_ORDER, fe_flux, ierr)
        if (ierr /= 0) error stop "PETSc could not build the RT0 flux basis"
        call PetscFECreateDefault(PETSC_COMM_WORLD, dimension, 1_PETSC_INT_KIND, PETSC_TRUE, &
                                  "pressure_", QUADRATURE_ORDER, fe_pressure, ierr)
        if (ierr /= 0) error stop "PETSc could not build the P0 pressure basis"

        ! One quadrature for both fields, so the two halves of the saddle-point
        ! system are integrated at the same points.
        call PetscFECopyQuadrature(fe_flux, fe_pressure, ierr)
        if (ierr /= 0) error stop "PETSc could not share the quadrature between the fields"

        discretisation = PetscObjectCast(fe_flux)
        call PetscObjectSetName(discretisation, "flux", ierr)
        if (ierr == 0) call DMSetField(self%dm, 0_PETSC_INT_KIND, PETSC_NULL_DMLABEL, discretisation, ierr)
        if (ierr /= 0) error stop "PETSc could not attach the flux field"

        discretisation = PetscObjectCast(fe_pressure)
        call PetscObjectSetName(discretisation, "pressure", ierr)
        if (ierr == 0) call DMSetField(self%dm, 1_PETSC_INT_KIND, PETSC_NULL_DMLABEL, discretisation, ierr)
        if (ierr /= 0) error stop "PETSc could not attach the pressure field"

        call PetscFEDestroy(fe_flux, ierr)
        call PetscFEDestroy(fe_pressure, ierr)

        call DMCreateDS(self%dm, ierr)
        if (ierr == 0) call DMGetDS(self%dm, ds, ierr)
        if (ierr /= 0) error stop "PETSc could not build the mixed discrete system"
        self%ds = ds

        dm_handle = transfer(self%dm%v, dm_handle)
        ds_handle = transfer(ds%v, ds_handle)

        ! Darcy's law: the resistance term against the test function, the
        ! pressure against its divergence.
        call petscds_set_residual(ds_handle, 0_c_int, c_funloc(f0_flux), c_funloc(f1_flux), cerr)
        if (cerr /= 0) error stop "PETSc rejected the Darcy residual"

        ! Conservation. The divergence is a test-value term because the P0 test
        ! functions have no element gradient at all: written as f1 of this
        ! equation the transport would integrate to exactly zero.
        call petscds_set_residual(ds_handle, 1_c_int, c_funloc(f0_conservation), c_null_funptr, cerr)
        if (cerr /= 0) error stop "PETSc rejected the conservation residual"

        call petscds_set_jacobian(ds_handle, 0_c_int, 0_c_int, &
                                  c_funloc(g0_flux_flux), c_null_funptr, c_null_funptr, c_null_funptr, cerr)
        if (cerr == 0) call petscds_set_jacobian(ds_handle, 0_c_int, 1_c_int, &
                                                 c_null_funptr, c_null_funptr, &
                                                 c_funloc(g2_flux_pressure), c_null_funptr, cerr)
        if (cerr == 0) call petscds_set_jacobian(ds_handle, 1_c_int, 0_c_int, &
                                                 c_null_funptr, c_funloc(g1_pressure_flux), &
                                                 c_null_funptr, c_null_funptr, cerr)
        if (cerr == 0) call petscds_set_jacobian(ds_handle, 1_c_int, 1_c_int, &
                                                 c_funloc(g0_pressure_pressure), c_null_funptr, &
                                                 c_null_funptr, c_null_funptr, cerr)
        if (cerr /= 0) error stop "PETSc rejected the mixed Jacobian"

        call store_constants(self, ierr)
        if (ierr /= 0) error stop "PETSc rejected the mixed run constants"

        call build_auxiliary_fields(self, dimension, ierr)
        if (ierr /= 0) error stop "PETSc could not build the mixed auxiliary fields"

        call add_boundaries(self, dm_handle, ierr)
        if (ierr /= 0) error stop "PETSc rejected a mixed boundary"

        call DMCreateGlobalVector(self%dm, self%solution, ierr)
        if (ierr == 0) call DMCreateGlobalVector(self%dm, self%old_solution, ierr)
        if (ierr == 0) call DMCreateLocalVector(self%dm, self%local_solution, ierr)
        if (ierr /= 0) error stop "PETSc could not create the mixed vectors"
        call VecZeroEntries(self%solution, ierr)
        if (ierr == 0) call VecZeroEntries(self%old_solution, ierr)

        ! The residual and the Jacobian come from the discrete system through
        ! the DM, so SNES needs neither a callback nor a matrix of its own.
        call SNESCreate(PETSC_COMM_WORLD, self%snes, ierr)
        if (ierr == 0) call SNESSetDM(self%snes, self%dm, ierr)
        if (ierr == 0) call DMPlexSetSNESLocalFEM(self%dm, PETSC_FALSE, PETSC_NULL_INTEGER, ierr)
        if (ierr == 0) call SNESSetFromOptions(self%snes, ierr)
        if (ierr /= 0) error stop "PETSc could not set up the mixed nonlinear solver"

        self%ready = .true.
    end subroutine initialize_mixed_hydraulic

    !> The options that turn PetscFECreateDefault into RT0 and P0.
    !>
    !> The primal degree and the dual order are not interchangeable here: the
    !> trimmed space of form degree minus one and primal degree zero carries the
    !> dual space of order one, which is where the one moment per facet comes
    !> from. A vector Lagrange element of degree one would give three degrees of
    !> freedom per vertex instead, and no normal continuity at all.
    subroutine set_mixed_element_options(ierr)
        implicit none
        PetscErrorCode, intent(inout) :: ierr

        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-flux_petscspace_type", "ptrimmed", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-flux_petscspace_degree", "0", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-flux_petscspace_components", "2", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-flux_petscspace_ptrimmed_form_degree", "-1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-flux_petscdualspace_type", "lagrange", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-flux_petscdualspace_order", "1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-flux_petscdualspace_form_degree", "-1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-flux_petscdualspace_lagrange_trimmed", "true", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-flux_petscdualspace_lagrange_continuity", "true", ierr)

        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-pressure_petscspace_type", "poly", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-pressure_petscspace_degree", "0", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-pressure_petscspace_components", "1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-pressure_petscdualspace_type", "lagrange", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-pressure_petscdualspace_order", "0", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-pressure_petscdualspace_form_degree", "0", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-pressure_petscdualspace_lagrange_continuity", "false", ierr)
    end subroutine set_mixed_element_options

    !> Write the run constants into the discrete system.
    subroutine store_constants(self, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        self%constants(C_DT) = real(self%dt, PETSC_REAL_KIND)
        self%constants(C_FLUX_SCALE) = real(self%flux_scale, PETSC_REAL_KIND)
        self%constants(C_PRESSURE_SCALE) = real(self%pressure_scale, PETSC_REAL_KIND)
        self%constants(C_LENGTH_SCALE) = real(self%length_scale, PETSC_REAL_KIND)
        self%constants(C_DENSITY) = real(self%density, PETSC_REAL_KIND)
        self%constants(C_GRAVITY_X) = real(self%gravity(1), PETSC_REAL_KIND)
        self%constants(C_GRAVITY_Z) = real(self%gravity(2), PETSC_REAL_KIND)

        call PetscDSSetConstants(self%ds, int(MIXED_CONSTANTS, PETSC_INT_KIND), self%constants, ierr)
    end subroutine store_constants

    !> The cell-wise coefficients, on a DM of their own.
    subroutine build_auxiliary_fields(self, dimension, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        PetscInt, intent(in) :: dimension
        PetscErrorCode, intent(inout) :: ierr

        call DMClone(self%dm, self%aux_dm, ierr)
        if (ierr /= 0) return

        ! Degree zero everywhere: one material per cell, and the boundary data
        ! of a facet is read from the cell it belongs to.
        call add_auxiliary_field(self%aux_dm, int(A_INVERSE_MOBILITY - 1, PETSC_INT_KIND), "inverse_mobility", &
                                 dimension, int(MIXED_DIMENSION*MIXED_DIMENSION, PETSC_INT_KIND), ierr)
        if (ierr == 0) call add_auxiliary_field(self%aux_dm, int(A_STORAGE - 1, PETSC_INT_KIND), "storage", &
                                                dimension, 1_PETSC_INT_KIND, ierr)
        if (ierr == 0) call add_auxiliary_field(self%aux_dm, int(A_SOURCE - 1, PETSC_INT_KIND), "source", &
                                                dimension, 1_PETSC_INT_KIND, ierr)
        if (ierr == 0) call add_auxiliary_field(self%aux_dm, int(A_PRESSURE_OLD - 1, PETSC_INT_KIND), "pressure_old", &
                                                dimension, 1_PETSC_INT_KIND, ierr)
        if (ierr == 0) call add_auxiliary_field(self%aux_dm, int(A_BC_PRESSURE - 1, PETSC_INT_KIND), "bc_pressure", &
                                                dimension, 1_PETSC_INT_KIND, ierr)
        if (ierr /= 0) return

        call DMCreateDS(self%aux_dm, ierr)
        if (ierr == 0) call DMCreateLocalVector(self%aux_dm, self%auxiliary, ierr)
        if (ierr == 0) call VecSet(self%auxiliary, 0.0_PETSC_REAL_KIND, ierr)
        if (ierr == 0) call DMSetAuxiliaryVec(self%dm, PETSC_NULL_DMLABEL, 0_PETSC_INT_KIND, &
                                              0_PETSC_INT_KIND, self%auxiliary, ierr)
        if (ierr /= 0) return

        call store_material_coefficients(self, ierr)
        if (ierr /= 0) return

        call DMGetLabel(self%dm, "Face Sets", self%face_sets, ierr)
    end subroutine build_auxiliary_fields

    !> One auxiliary field of `num_components` cell constants.
    subroutine add_auxiliary_field(dm, field, name, dimension, num_components, ierr)
        implicit none
        DM, intent(inout) :: dm
        PetscInt, intent(in) :: field
        character(len=*), intent(in) :: name
        PetscInt, intent(in) :: dimension
        PetscInt, intent(in) :: num_components
        PetscErrorCode, intent(inout) :: ierr

        PetscFE :: fe
        PetscObject :: discretisation

        call PetscFECreateLagrange(PETSC_COMM_WORLD, dimension, num_components, PETSC_TRUE, &
                                   0_PETSC_INT_KIND, QUADRATURE_ORDER, fe, ierr)
        if (ierr /= 0) return

        discretisation = PetscObjectCast(fe)
        call PetscObjectSetName(discretisation, name, ierr)
        if (ierr == 0) call DMSetField(dm, field, PETSC_NULL_DMLABEL, discretisation, ierr)
        call PetscFEDestroy(fe, ierr)
    end subroutine add_auxiliary_field

    !> Invert each cell's own conductivity into its auxiliary field.
    !>
    !> The cell's material is its "Cell Sets" value, the same identifier the
    !> input layer uses. A cell without one, or with a conductivity that is not
    !> positive definite, stops the run and says which cell and which material:
    !> an inverse mobility is what the weak form integrates, and there is no
    !> value of it that stands in for an impermeable cell.
    subroutine store_material_coefficients(self, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: cell, c_start, c_end, offset, num_dofs, material
        real(real64) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION), gravity_magnitude
        integer(int32) :: m, found, k
        logical :: ok

        gravity_magnitude = norm2(self%gravity)
        ! A run without gravity still needs the rho g that turns a conductivity
        ! into a pressure mobility; the standard value is the one the
        ! conductivity was measured with.
        if (gravity_magnitude <= 0.0_real64) gravity_magnitude = 9.80665_real64

        call DMGetLocalSection(self%aux_dm, section, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(self%aux_dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return

        call VecGetArray(self%auxiliary, raw, ierr)
        if (ierr /= 0) return

        do cell = c_start, c_end - 1
            call PetscSectionGetFieldDof(section, cell, int(A_INVERSE_MOBILITY - 1, PETSC_INT_KIND), num_dofs, ierr)
            if (ierr /= 0 .or. num_dofs <= 0) cycle

            call DMGetLabelValue(self%aux_dm, "Cell Sets", cell, material, ierr)
            if (ierr /= 0 .or. material < 0) then
                call VecRestoreArray(self%auxiliary, raw, ierr)
                error stop "A cell carries no material: the mesh has a region without a physical group"
            end if

            found = 0
            do m = 1, size(self%materials)
                if (self%materials(m)%id == int(material, int32)) found = m
            end do
            if (found == 0) then
                call VecRestoreArray(self%auxiliary, raw, ierr)
                write (*, '(a,i0)') "No mixed hydraulic material is defined for Cell Sets value ", int(material, int32)
                error stop "Unknown material"
            end if

            call invert_mobility(self%materials(found)%conductivity, self%density, gravity_magnitude, &
                                 inverse_mobility, ok)
            if (.not. ok) then
                call VecRestoreArray(self%auxiliary, raw, ierr)
                write (*, '(a,i0,a,i0,a)') "Material ", self%materials(found)%id, " of cell ", int(cell, int32), &
                    " has a conductivity that is not symmetric positive definite"
                error stop "An impermeable region is a different model, not a small conductivity"
            end if

            call PetscSectionGetFieldOffset(section, cell, int(A_INVERSE_MOBILITY - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr /= 0) cycle
            do k = 1, MIXED_DIMENSION*MIXED_DIMENSION
                raw(offset + k) = real(inverse_mobility(k), PETSC_REAL_KIND)
            end do

            call PetscSectionGetFieldOffset(section, cell, int(A_STORAGE - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr == 0) raw(offset + 1) = real(self%materials(found)%storage, PETSC_REAL_KIND)
            call PetscSectionGetFieldOffset(section, cell, int(A_SOURCE - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr == 0) raw(offset + 1) = real(self%materials(found)%source, PETSC_REAL_KIND)
        end do

        call VecRestoreArray(self%auxiliary, raw, ierr)
    end subroutine store_material_coefficients

    !> Register every boundary the caller asked for.
    subroutine add_boundaries(self, dm_handle, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        type(c_ptr), intent(in) :: dm_handle
        PetscErrorCode, intent(inout) :: ierr

        integer(int32) :: b
        integer(c_int) :: boundary_index, cerr
        character(len=32) :: name

        call reset_flux_boundaries()
        ierr = 0

        do b = 1, self%num_boundaries
            select case (self%boundaries(b)%kind)
            case (MIXED_BC_PRESSURE)
                write (name, '(a,i0)') "pressure-", self%boundaries(b)%boundary_id
                call petscds_add_boundary_natural(dm_handle, trim(name), "Face Sets", &
                                                  int(self%boundaries(b)%boundary_id, c_int), 0_c_int, &
                                                  c_funloc(bd_f0_flux_pressure), c_null_funptr, &
                                                  boundary_index, cerr)
                if (cerr /= 0) then
                    ierr = int(cerr)
                    return
                end if
                call attach_boundary_pressure(self, self%boundaries(b)%boundary_id, &
                                              self%boundaries(b)%pressure, ierr)
                if (ierr /= 0) return

            case (MIXED_BC_FLUX)
                if (b > MAX_FLUX_BOUNDARIES) then
                    error stop "More prescribed-flux boundaries than the kernels can hold"
                end if
                call set_flux_boundary_value(b, self%boundaries(b)%flux, self%flux_scale)
                write (name, '(a,i0)') "flux-", self%boundaries(b)%boundary_id
                call petscds_add_boundary_essential(dm_handle, trim(name), "Face Sets", &
                                                    int(self%boundaries(b)%boundary_id, c_int), 0_c_int, &
                                                    c_funloc(bc_essential_flux), &
                                                    flux_boundary_context(b), &
                                                    boundary_index, cerr)
                if (cerr /= 0) then
                    ierr = int(cerr)
                    return
                end if

            case default
                error stop "A boundary of this kind is not one the mixed hydraulic solver knows"
            end select
        end do
    end subroutine add_boundaries

    !> Give one prescribed-pressure boundary its own auxiliary vector.
    subroutine attach_boundary_pressure(self, boundary_id, pressure, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        integer(int32), intent(in) :: boundary_id
        real(real64), intent(in) :: pressure
        PetscErrorCode, intent(inout) :: ierr

        Vec :: boundary_aux

        if (self%num_boundary_aux >= MAX_MIXED_BOUNDARIES) then
            error stop "More natural boundaries than the mixed hydraulic solver can hold"
        end if

        call DMCreateLocalVector(self%aux_dm, boundary_aux, ierr)
        if (ierr /= 0) return
        ! The coefficients still have to be there: the boundary integral reads
        ! the cell behind the facet.
        call VecCopy(self%auxiliary, boundary_aux, ierr)
        if (ierr == 0) call fill_auxiliary_field(self%aux_dm, boundary_aux, A_BC_PRESSURE, pressure)
        if (ierr == 0) call DMSetAuxiliaryVec(self%dm, self%face_sets, int(boundary_id, PETSC_INT_KIND), &
                                              0_PETSC_INT_KIND, boundary_aux, ierr)
        if (ierr /= 0) return

        self%num_boundary_aux = self%num_boundary_aux + 1
        self%boundary_aux(self%num_boundary_aux) = boundary_aux
    end subroutine attach_boundary_pressure

    !> Set every degree of freedom of one auxiliary field to the same value.
    subroutine fill_auxiliary_field(dm, vec, field, value)
        implicit none
        DM, intent(in) :: dm
        Vec, intent(inout) :: vec
        integer(int32), intent(in) :: field
        real(real64), intent(in) :: value

        PetscErrorCode :: ierr
        PetscSection :: section
        PetscInt :: point, p_start, p_end, offset, num_dofs, k
        PetscScalar, pointer :: raw(:)

        call DMGetLocalSection(dm, section, ierr)
        if (ierr /= 0) return
        call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr /= 0) return

        call VecGetArray(vec, raw, ierr)
        if (ierr /= 0) return
        do point = p_start, p_end - 1
            call PetscSectionGetFieldDof(section, point, int(field - 1, PETSC_INT_KIND), num_dofs, ierr)
            if (ierr /= 0 .or. num_dofs <= 0) cycle
            call PetscSectionGetFieldOffset(section, point, int(field - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr /= 0) cycle
            do k = 1, num_dofs
                raw(offset + k) = real(value, PETSC_REAL_KIND)
            end do
        end do
        call VecRestoreArray(vec, raw, ierr)
    end subroutine fill_auxiliary_field

    module subroutine destroy_mixed_hydraulic(self, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        integer(int32) :: b

        ierr = 0
        if (.not. self%ready) return

        call SNESDestroy(self%snes, ierr)
        call VecDestroy(self%solution, ierr)
        call VecDestroy(self%old_solution, ierr)
        call VecDestroy(self%local_solution, ierr)
        do b = 1, self%num_boundary_aux
            call VecDestroy(self%boundary_aux(b), ierr)
        end do
        call VecDestroy(self%auxiliary, ierr)
        call DMDestroy(self%aux_dm, ierr)
        call DMDestroy(self%dm, ierr)

        self%num_boundary_aux = 0
        self%num_boundaries = 0
        self%ready = .false.
    end subroutine destroy_mixed_hydraulic

end submodule mixed_hydraulic_setup
