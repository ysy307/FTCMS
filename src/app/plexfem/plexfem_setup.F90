#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscsnes.h>

!> Construction of the DM, the discrete system and the boundary conditions.
submodule(app_plexfem) plexfem_setup
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_funptr, c_funloc, c_null_funptr
    use :: module_core, only:PHYSICS_TYPES, THERMAL_BC_TYPES, HYDRAULIC_BC_TYPES, &
                             IC_TARGETS, IC_METHODS, &
                             type_config_bc, type_config_ic
    use :: module_input, only:input_translator
    use :: core_interop_petscds, only:petscds_set_residual, petscds_set_jacobian, &
                                      petscds_add_boundary_essential, petscds_add_boundary_natural
    use :: app_plexfem_physics, only:PLEXFEM_CONSTANTS, C_POROSITY, C_PRESSURE, &
                                     PLEXFEM_AUX_FIELDS, A_MATERIAL_ID, C_MEASURING_MATERIAL, &
                                     A_BC_KIND, A_BC_VALUE_1, A_BC_VALUE_2, &
                                     A_DARCY, A_COUPLING, A_GRAVITY, A_POTENTIAL, &
                                     A_TEMPERATURE_OLD, A_PRESSURE_OLD, &
                                     A_TEMPERATURE_OLDER, A_PRESSURE_OLDER, &
                                     BC_FLUX, BC_ROBIN, &
                                     set_dirichlet_value, dirichlet_context, attach_thermal_model, &
                                     attach_hydraulic_model, hydraulic_is_active
    use :: app_plexfem_pointwise, only:f1_energy, g2_energy, g2_energy_p, g3_energy, f0_material_indicator, &
                                       f0_energy, f0_water, g0_energy, g0_energy_p, g0_water_p, g0_water_t, &
                                       f1_water, g3_water_pp, g3_water_pt, &
                                       g2_water_p, g2_water_t, &
                                       bd_f0_energy, bd_g0_energy, &
                                       bc_essential_temperature, bc_essential_pressure
    implicit none

contains

    module subroutine initialize_plexfem(self)
        implicit none
        class(type_plexfem), intent(inout), target :: self

        type(type_config_bc), allocatable :: config_bcs(:)

        PetscErrorCode :: ierr
        PetscFE :: fe
        PetscDS :: ds
        PetscObject :: discretisation
        PetscBool :: is_simplex
        PetscInt :: dimension, quadrature_order
        DM :: base_dm
        type(c_ptr) :: dm_handle, ds_handle
        integer(c_int) :: cerr
        integer(int32) :: element_order
        integer(int32), allocatable :: active_region_ids(:)

        call self%input%initialize()

        ! One gate for the whole path: PETSc numbers a degree-one basis on the
        ! vertices, which is what every mapping below assumes.
        call self%input%geometry%mesh%get_active_region_info(active_region_ids, &
                                                            target_dim=self%input%basic%simulation_settings%calculate_dimension)
        call self%thermal%initialize(self%input, active_region_ids)
        call attach_thermal_model(self%thermal)

        ! Pressure becomes a field of its own when the project asks for water.
        if (self%input%basic%analysis_controls%is_active(PHYSICS_TYPES%HYDRAULIC%ID)) then
            call self%hydraulic%initialize(self%input, active_region_ids)
            call attach_hydraulic_model(self%hydraulic)
        end if

        call self%input%geometry%mesh%get_element_shape(element_order, is_simplex)
        if (element_order /= 1) then
            write (*, '(a,i0,a)') "The PetscDS path supports first-order meshes only; "// &
                "the mesh coordinates are of degree ", element_order, "."
            error stop "Unsupported mesh order"
        end if

        base_dm = self%input%geometry%mesh%get_dm()
        call DMClone(base_dm, self%dm, ierr)
        if (ierr /= 0) error stop "PETSc could not clone the mesh for the PetscDS path"
        self%dm_ready = .true.

        call DMGetDimension(self%dm, dimension, ierr)
        if (ierr /= 0) error stop "The mesh has no dimension"

        ! The same rule the older path integrates on, so the two never differ
        ! by the number of quadrature points.
        quadrature_order = int(max(self%input%basic%geometry_settings%integration_order, 1), PETSC_INT_KIND)

        call PetscFECreateLagrange(PETSC_COMM_WORLD, dimension, 1_PETSC_INT_KIND, is_simplex, &
                                   1_PETSC_INT_KIND, quadrature_order, fe, ierr)
        if (ierr /= 0) error stop "PETSc could not build the temperature basis"

        discretisation = PetscObjectCast(fe)
        call PetscObjectSetName(discretisation, "temperature", ierr)
        if (ierr == 0) call DMSetField(self%dm, 0_PETSC_INT_KIND, PETSC_NULL_DMLABEL, discretisation, ierr)
        if (ierr /= 0) error stop "PETSc could not attach the temperature field"
        call PetscFEDestroy(fe, ierr)

        if (hydraulic_is_active()) then
            call PetscFECreateLagrange(PETSC_COMM_WORLD, dimension, 1_PETSC_INT_KIND, is_simplex, &
                                       1_PETSC_INT_KIND, quadrature_order, fe, ierr)
            if (ierr /= 0) error stop "PETSc could not build the pressure basis"

            discretisation = PetscObjectCast(fe)
            call PetscObjectSetName(discretisation, "pressure", ierr)
            if (ierr == 0) call DMSetField(self%dm, 1_PETSC_INT_KIND, PETSC_NULL_DMLABEL, discretisation, ierr)
            if (ierr /= 0) error stop "PETSc could not attach the pressure field"
            call PetscFEDestroy(fe, ierr)
        end if

        call DMCreateDS(self%dm, ierr)
        if (ierr == 0) call DMGetDS(self%dm, ds, ierr)
        if (ierr /= 0) error stop "PETSc could not build the discrete system"

        ! PETSc objects reach the C shim as the address their Fortran handle
        ! carries, the same way domain_fe_petsc passes a PetscFE.
        dm_handle = transfer(self%dm%v, dm_handle)
        ds_handle = transfer(ds%v, ds_handle)

        ! Storage and conduction, both integrated by PETSc. The storage term
        ! differences the enthalpy across the history levels, which reach the
        ! kernel as auxiliary fields, so it is an ordinary f0 rather than
        ! something assembled node by node outside the discrete system.
        call petscds_set_residual(ds_handle, 0_c_int, c_funloc(f0_energy), c_funloc(f1_energy), cerr)
        if (cerr /= 0) error stop "PETSc rejected the energy residual"

        ! g0 is the storage term's tangent, consistent rather than lumped: only
        ! the current level depends on the unknown, so it is the BDF coefficient
        ! of that level times the capacity. g1 stays unset: the storage term has
        ! no gradient in it.
        !
        ! g2 = (dlambda/dT) grad T costs thirty percent of the run, because its
        ! central difference needs three more evaluations of the whole
        ! constitutive chain per quadrature point, and on _1Domain it bought
        ! nothing: the iteration count stayed at four to five per step. It is
        ! registered all the same, because _1Domain is thermal-only and its
        ! conductivity barely moves. Where water freezes, lambda travels between
        ! the ice and the liquid value, and this is the only temperature
        ! dependence of the conduction tangent there is.
        call petscds_set_jacobian(ds_handle, 0_c_int, 0_c_int, &
                                  c_funloc(g0_energy), c_null_funptr, &
                                  c_funloc(g2_energy), c_funloc(g3_energy), cerr)
        if (cerr /= 0) error stop "PETSc rejected the energy Jacobian"

        if (hydraulic_is_active()) then
            ! The energy equation's pressure column. Conduction is a function of
            ! the phase composition, and pressure is half of what sets it, so
            ! without this the energy row has no derivative with respect to the
            ! other unknown at all.
            call petscds_set_jacobian(ds_handle, 0_c_int, 1_c_int, &
                                      c_funloc(g0_energy_p), c_null_funptr, &
                                      c_funloc(g2_energy_p), c_null_funptr, cerr)
            if (cerr /= 0) error stop "PETSc rejected the energy-pressure Jacobian"

            ! Storage and the mixed-form flux, and the derivative of each.
            !
            ! The storage tangent is not optional here. Every water boundary is
            ! zero flux, so a tangent made of gradient terms alone sums to zero
            ! over the rows -- the test functions sum to one -- and the matrix
            ! acquires a left null space. That is a singular operator, not a
            ! badly preconditioned one, and no amount of scaling reaches it.
            !
            ! Registering only g3 would also freeze the Darcy conductivity at the
            ! current iterate, making the water block a Picard iteration rather
            ! than a Newton one: the ice impedance moves the conductivity by a
            ! factor of several within one step, and that factor becomes the
            ! contraction rate.
            call petscds_set_residual(ds_handle, 1_c_int, c_funloc(f0_water), c_funloc(f1_water), cerr)
            ! The derivatives of the Darcy coefficients themselves, g2_water_p
            ! and g2_water_t, are written and available but not registered.
            ! Measured, adding them left the water block contracting at the same
            ! two and a half percent per iteration -- they buy nothing here --
            ! while costing most of the Jacobian: each asks for the coefficients
            ! at four perturbed states, and every one of those evaluates the
            ! Mualem relative permeability, whose incomplete beta function is a
            ! continued fraction of a few dozen terms. Sixteen of those per
            ! quadrature point is where the Jacobian was spending four times what
            ! the residual does.
            if (cerr == 0) call petscds_set_jacobian(ds_handle, 1_c_int, 1_c_int, &
                                                     c_funloc(g0_water_p), c_null_funptr, &
                                                     c_null_funptr, c_funloc(g3_water_pp), cerr)
            if (cerr == 0) call petscds_set_jacobian(ds_handle, 1_c_int, 0_c_int, &
                                                     c_funloc(g0_water_t), c_null_funptr, &
                                                     c_null_funptr, c_funloc(g3_water_pt), cerr)
            if (cerr /= 0) error stop "PETSc rejected the water residual"
        end if

        self%ds = ds
        call build_constants(self, self%constants)
        call PetscDSSetConstants(ds, int(PLEXFEM_CONSTANTS, PETSC_INT_KIND), self%constants, ierr)
        if (ierr /= 0) error stop "PETSc rejected the material constants"

        ! The natural conditions read their constants from a boundary's own
        ! auxiliary vector, so the fields have to exist before they are set.
        call build_auxiliary_field(self, is_simplex, dimension, quadrature_order)

        call input_translator%execute(self%input, PHYSICS_TYPES%THERMAL, config_bcs)
        call add_boundaries(self, dm_handle, config_bcs, PHYSICS_TYPES%THERMAL%ID, 0_c_int)
        if (allocated(config_bcs)) deallocate (config_bcs)

        if (hydraulic_is_active()) then
            call input_translator%execute(self%input, PHYSICS_TYPES%HYDRAULIC, config_bcs)
            call add_boundaries(self, dm_handle, config_bcs, PHYSICS_TYPES%HYDRAULIC%ID, 1_c_int)
            if (allocated(config_bcs)) deallocate (config_bcs)
        end if

        call DMCreateGlobalVector(self%dm, self%solution, ierr)
        if (ierr == 0) call DMCreateGlobalVector(self%dm, self%old_solution, ierr)
        if (ierr == 0) call DMCreateGlobalVector(self%dm, self%older_solution, ierr)
        if (ierr /= 0) error stop "PETSc could not create the solution vectors"

        call apply_initial_condition(self)
        ! Both history levels start at the initial state. The first step is
        ! first order and never reads the older one, but a vector the residual
        ! can reach is not left holding whatever memory it was given.
        call VecCopy(self%solution, self%old_solution, ierr)
        if (ierr == 0) call VecCopy(self%solution, self%older_solution, ierr)
        if (ierr /= 0) error stop "PETSc could not seed the previous states"

        ! After the boundaries, because DMAddBoundary is refused once a
        ! residual has been assembled: assembling one creates the local
        ! section. The kernels stand down while this runs, so what comes back
        ! is the volume alone.
        call build_nodal_volumes(self, region_ids_of(self))
        call PetscDSSetConstants(ds, int(PLEXFEM_CONSTANTS, PETSC_INT_KIND), self%constants, ierr)
        if (ierr /= 0) error stop "PETSc rejected the material constants"
    end subroutine initialize_plexfem

    !> The materials this mesh is made of.
    function region_ids_of(self) result(region_ids)
        implicit none
        class(type_plexfem), intent(inout) :: self
        integer(int32), allocatable :: region_ids(:)

        call self%input%geometry%mesh%get_active_region_info( &
            region_ids, target_dim=self%input%basic%simulation_settings%calculate_dimension)
    end function region_ids_of

    !> Nodal volumes A_i = \int psi_i, one vector per material.
    !>
    !> These are the weights of the lumped storage term. They are measured by
    !> integrating an indicator of each material against the test functions,
    !> which is the same quadrature and the same element geometry PETSc uses
    !> for every other term -- there is no second implementation of either.
    !>
    !> One vector per material rather than one in total, because a node on the
    !> boundary between two materials stores enthalpy of both, in the
    !> proportion of the volume each takes up around it.
    subroutine build_nodal_volumes(self, region_ids)
        implicit none
        class(type_plexfem), intent(inout) :: self
        integer(int32), intent(in) :: region_ids(:)

        PetscErrorCode :: ierr
        PetscDS :: ds
        type(c_ptr) :: ds_handle
        integer(c_int) :: cerr
        Vec :: local_x, local_volume
        integer(int32) :: m

        ! The discrete system is fetched again rather than carried down from
        ! the top: registering the boundaries can leave the DM holding a
        ! different one, and a stale handle takes the kernels nowhere.
        call DMGetDS(self%dm, ds, ierr)
        if (ierr /= 0) error stop "The mesh has no discrete system"
        ds_handle = transfer(ds%v, ds_handle)
        self%ds = ds

        self%num_materials = size(region_ids)

        ! The liquid-pressure auxiliary field is one value per node, and the
        ! effective suction a node carries depends on the material it is being
        ! asked about. With one material that is unambiguous; with two meeting at
        ! a node it is not, and a single degree-one field would quietly keep
        ! whichever was written last and drive the water flux with the wrong
        ! potential on one side of the interface. Checked here rather than where
        ! the fields are built, because that runs before the count is known.
        if (hydraulic_is_active() .and. self%num_materials > 1) then
            error stop "The PetscDS water path supports one material until the "// &
                "liquid-pressure potential is carried per material"
        end if
        if (self%num_materials > size(self%material_id)) then
            error stop "More materials than the PetscDS path can hold"
        end if

        ! The indicator stands in for the residual while the weights are made,
        ! and every other equation stands down: what is wanted is the volume
        ! around each node, not a flux evaluated at a state nobody solved for.
        call petscds_set_residual(ds_handle, 0_c_int, c_funloc(f0_material_indicator), &
                                  c_null_funptr, cerr)
        ! Every equation measures, none of them fluxes: an equation left with
        ! no residual at all makes PETSc skip the assembly it is part of.
        if (cerr == 0 .and. hydraulic_is_active()) then
            call petscds_set_residual(ds_handle, 1_c_int, c_funloc(f0_material_indicator), &
                                      c_null_funptr, cerr)
        end if
        if (cerr /= 0) error stop "PETSc rejected the material indicator"

        call DMGetLocalVector(self%dm, local_x, ierr)
        if (ierr == 0) call DMGetLocalVector(self%dm, local_volume, ierr)
        if (ierr /= 0) error stop "PETSc could not take a local vector"
        call VecZeroEntries(local_x, ierr)

        do m = 1, self%num_materials
            self%material_id(m) = region_ids(m)
            self%constants(C_MEASURING_MATERIAL) = real(region_ids(m), PETSC_REAL_KIND)
            call PetscDSSetConstants(ds, int(PLEXFEM_CONSTANTS, PETSC_INT_KIND), self%constants, ierr)
            if (ierr /= 0) error stop "PETSc rejected the material being measured"

            call VecZeroEntries(local_volume, ierr)
            call DMPlexSNESComputeResidualFEM(self%dm, local_x, local_volume, self%dm, ierr)
            if (ierr /= 0) error stop "PETSc could not measure the nodal volumes"

            ! Kept local, because that is the form the storage term is added
            ! in: a node shared by two ranks holds each rank's share, and the
            ! one scatter of the residual is what sums them.
            call DMCreateLocalVector(self%dm, self%nodal_volume(m), ierr)
            if (ierr == 0) call VecCopy(local_volume, self%nodal_volume(m), ierr)
            if (ierr /= 0) error stop "PETSc could not keep the nodal volumes"

        end do

        call DMRestoreLocalVector(self%dm, local_volume, ierr)
        call DMRestoreLocalVector(self%dm, local_x, ierr)

        ! Put the real residuals back -- storage as well as flux. Restoring only
        ! the flux leaves f0 at null and the time derivative silently absent.
        call petscds_set_residual(ds_handle, 0_c_int, c_funloc(f0_energy), c_funloc(f1_energy), cerr)
        if (cerr == 0 .and. hydraulic_is_active()) then
            call petscds_set_residual(ds_handle, 1_c_int, c_funloc(f0_water), c_funloc(f1_water), cerr)
        end if
        if (cerr /= 0) error stop "PETSc rejected the energy residual"
        self%constants(C_MEASURING_MATERIAL) = 0.0_PETSC_REAL_KIND
    end subroutine build_nodal_volumes

    !> The auxiliary field that carries the previous step's temperature.
    !>
    !> It lives on its own clone so the boundary conditions of the solve DM
    !> cannot reach it: PETSc would otherwise overwrite the stored temperature
    !> with the prescribed one wherever the two meet, and the enthalpy
    !> difference would lose the boundary's own history.
    subroutine build_auxiliary_field(self, is_simplex, dimension, quadrature_order)
        implicit none
        class(type_plexfem), intent(inout) :: self
        PetscBool, intent(in) :: is_simplex
        PetscInt, intent(in) :: dimension
        PetscInt, intent(in) :: quadrature_order

        PetscErrorCode :: ierr
        DM :: base_dm

        base_dm = self%input%geometry%mesh%get_dm()
        call DMClone(base_dm, self%aux_dm, ierr)
        if (ierr /= 0) error stop "PETSc could not clone the mesh for the auxiliary fields"

        ! Field order has to match the A_* constants: PETSc hands the kernels
        ! one offset per field, in the order the fields were added.
        ! Degree zero: one material per cell, constant over it.
        call add_auxiliary_field(self%aux_dm, A_MATERIAL_ID - 1, "material_id", &
                                 dimension, is_simplex, 0_PETSC_INT_KIND, quadrature_order)
        ! Constant over a boundary, so degree zero is enough for all three.
        call add_auxiliary_field(self%aux_dm, A_BC_KIND - 1, "bc_kind", &
                                 dimension, is_simplex, 0_PETSC_INT_KIND, quadrature_order)
        call add_auxiliary_field(self%aux_dm, A_BC_VALUE_1 - 1, "bc_value_1", &
                                 dimension, is_simplex, 0_PETSC_INT_KIND, quadrature_order)
        call add_auxiliary_field(self%aux_dm, A_BC_VALUE_2 - 1, "bc_value_2", &
                                 dimension, is_simplex, 0_PETSC_INT_KIND, quadrature_order)
        ! One representative transport coefficient per cell, written afresh from
        ! the current iterate before every assembly.
        call add_auxiliary_field(self%aux_dm, A_DARCY - 1, "darcy_element", &
                                 dimension, is_simplex, 0_PETSC_INT_KIND, quadrature_order)
        call add_auxiliary_field(self%aux_dm, A_COUPLING - 1, "coupling_element", &
                                 dimension, is_simplex, 0_PETSC_INT_KIND, quadrature_order)
        call add_auxiliary_field(self%aux_dm, A_GRAVITY - 1, "gravity_element", &
                                 dimension, is_simplex, 0_PETSC_INT_KIND, quadrature_order)
        ! Degree one: its gradient is what the water flux is driven by, so it has
        ! to be a field the assembly can differentiate, not a cell constant.
        call add_auxiliary_field(self%aux_dm, A_POTENTIAL - 1, "liquid_pressure", &
                                 dimension, is_simplex, 1_PETSC_INT_KIND, quadrature_order)
        ! The history the storage term differences against, degree one so that
        ! PETSc interpolates it to the quadrature points the residual is
        ! integrated on.
        call add_auxiliary_field(self%aux_dm, A_TEMPERATURE_OLD - 1, "temperature_old", &
                                 dimension, is_simplex, 1_PETSC_INT_KIND, quadrature_order)
        call add_auxiliary_field(self%aux_dm, A_PRESSURE_OLD - 1, "pressure_old", &
                                 dimension, is_simplex, 1_PETSC_INT_KIND, quadrature_order)
        call add_auxiliary_field(self%aux_dm, A_TEMPERATURE_OLDER - 1, "temperature_older", &
                                 dimension, is_simplex, 1_PETSC_INT_KIND, quadrature_order)
        call add_auxiliary_field(self%aux_dm, A_PRESSURE_OLDER - 1, "pressure_older", &
                                 dimension, is_simplex, 1_PETSC_INT_KIND, quadrature_order)

        call DMCreateDS(self%aux_dm, ierr)
        if (ierr /= 0) error stop "PETSc could not set up the auxiliary fields"

        ! PETSc reads an auxiliary field from a local vector.
        call DMCreateLocalVector(self%aux_dm, self%auxiliary, ierr)
        if (ierr == 0) call VecSet(self%auxiliary, 0.0_PETSC_REAL_KIND, ierr)
        if (ierr == 0) call DMSetAuxiliaryVec(self%dm, PETSC_NULL_DMLABEL, 0_PETSC_INT_KIND, &
                                              0_PETSC_INT_KIND, self%auxiliary, ierr)
        if (ierr /= 0) error stop "PETSc could not attach the auxiliary fields"

        call store_material_ids(self)

        call DMGetLabel(self%dm, "Face Sets", self%face_sets, ierr)
        if (ierr /= 0) error stop "The mesh carries no boundary label"
    end subroutine build_auxiliary_field

    !> One auxiliary field of the given Lagrange degree.
    !>
    !> It integrates on the same rules as the solution does: PETSc evaluates the
    !> auxiliary fields at the very points it evaluates the solution at.
    subroutine add_auxiliary_field(dm, field, name, dimension, is_simplex, degree, quadrature_order)
        implicit none
        DM, intent(inout) :: dm
        PetscInt, intent(in) :: field
        character(len=*), intent(in) :: name
        PetscInt, intent(in) :: dimension
        PetscBool, intent(in) :: is_simplex
        PetscInt, intent(in) :: degree
        PetscInt, intent(in) :: quadrature_order

        PetscErrorCode :: ierr
        PetscFE :: fe
        PetscObject :: discretisation

        call PetscFECreateLagrange(PETSC_COMM_WORLD, dimension, 1_PETSC_INT_KIND, is_simplex, &
                                   degree, quadrature_order, fe, ierr)
        if (ierr /= 0) error stop "PETSc could not build an auxiliary basis"

        discretisation = PetscObjectCast(fe)
        call PetscObjectSetName(discretisation, name, ierr)
        if (ierr == 0) call DMSetField(dm, field, PETSC_NULL_DMLABEL, discretisation, ierr)
        if (ierr /= 0) error stop "PETSc could not attach an auxiliary field"
        call PetscFEDestroy(fe, ierr)
    end subroutine add_auxiliary_field

    !> Write each cell's material into its auxiliary degree of freedom.
    !>
    !> The material is the "Cell Sets" value PETSc kept from the Gmsh physical
    !> group, which is what the input layer uses as the material id everywhere
    !> else. A cell without one is a mesh the rest of the code would refuse
    !> too, so it is refused here rather than defaulted.
    subroutine store_material_ids(self)
        implicit none
        class(type_plexfem), intent(inout) :: self

        PetscErrorCode :: ierr
        PetscSection :: section
        PetscInt :: cell, c_start, c_end, offset, num_dofs, material
        PetscScalar, pointer :: raw(:)

        call DMGetLocalSection(self%aux_dm, section, ierr)
        if (ierr /= 0) error stop "The auxiliary fields have no layout"

        call DMPlexGetHeightStratum(self%aux_dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) error stop "The mesh has no cells"

        call VecGetArray(self%auxiliary, raw, ierr)
        if (ierr /= 0) error stop "PETSc could not open the auxiliary vector"

        do cell = c_start, c_end - 1
            call PetscSectionGetFieldDof(section, cell, int(A_MATERIAL_ID - 1, PETSC_INT_KIND), num_dofs, ierr)
            if (ierr /= 0 .or. num_dofs <= 0) cycle
            call PetscSectionGetFieldOffset(section, cell, int(A_MATERIAL_ID - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr /= 0) cycle

            call DMGetLabelValue(self%aux_dm, "Cell Sets", cell, material, ierr)
            if (ierr /= 0 .or. material < 0) then
                error stop "A cell carries no material: the mesh has a region without a physical group"
            end if
            raw(offset + 1) = real(material, PETSC_REAL_KIND)
        end do

        call VecRestoreArray(self%auxiliary, raw, ierr)

    end subroutine store_material_ids

    !> The material response of this stage travels as DS constants.
    !>
    !> One material, one porosity and one pressure: the project is thermal-only,
    !> so the pressure never moves from its initial value and the porosity is
    !> whatever the initial condition set. Anything else has to be refused
    !> rather than averaged, because a kernel reading a constant that is not
    !> constant would be quietly wrong everywhere it differs.
    subroutine build_constants(self, constants)
        implicit none
        class(type_plexfem), intent(inout) :: self
        PetscScalar, intent(inout) :: constants(:)

        type(type_config_ic) :: config_porosity
        integer(int32), allocatable :: region_ids(:)
        integer(int32) :: computation_dimension

        ! The material is carried per cell in an auxiliary field, so several
        ! regions are fine; only what is written here has to be one number for
        ! the whole mesh.

        computation_dimension = self%input%basic%simulation_settings%calculate_dimension
        call self%input%geometry%mesh%get_active_region_info(region_ids, target_dim=computation_dimension)

        call input_translator%execute(self%input, IC_TARGETS%POROSITY, config_porosity)
        if (config_porosity%ic_kind%ID /= IC_METHODS%UNIFORM%ID) then
            error stop "The PetscDS path needs a uniform initial porosity"
        end if

        constants(C_POROSITY) = real(config_porosity%value, PETSC_REAL_KIND)

        ! With the hydraulic physics off the older path never fills a pressure
        ! either: governing_base zeroes P_node and only overwrites it from a
        ! state that was set, so its quadrature points see zero. The freezing
        ! curve has to be given the same reference here, or the two paths would
        ! be solving different problems.
        constants(C_PRESSURE) = 0.0_PETSC_REAL_KIND
    end subroutine build_constants

    !> Register every thermal boundary of the project.
    !>
    !> A boundary the project leaves out is a natural boundary with no
    !> residual, which is the zero-flux condition itself, so nothing is
    !> registered for it. A zero flux is treated the same way; only a boundary
    !> that actually carries a flux or a Robin condition needs PETSc to
    !> integrate over it.
    subroutine add_boundaries(self, dm_handle, config_bcs, physics_id, field)
        implicit none
        class(type_plexfem), intent(inout) :: self
        type(c_ptr), intent(in) :: dm_handle
        type(type_config_bc), allocatable, intent(in) :: config_bcs(:)
        integer(int32), intent(in) :: physics_id
        integer(c_int), intent(in) :: field

        integer(int32) :: i
        integer(c_int) :: boundary_index, cerr
        real(real64) :: value_1, value_2
        character(len=32) :: name

        if (.not. allocated(config_bcs)) return

        do i = 1, size(config_bcs)
            if (config_bcs(i)%physics_type%ID /= physics_id) cycle
            if (.not. allocated(config_bcs(i)%values)) then
                error stop "A thermal boundary carries no value"
            end if
            if (config_bcs(i)%num_time_points > 1) then
                error stop "A time-dependent boundary is not supported by the PetscDS path yet"
            end if

            value_1 = config_bcs(i)%values(1, 1)
            value_2 = 0.0d0
            if (size(config_bcs(i)%values, 1) >= 2) value_2 = config_bcs(i)%values(2, 1)

            select case (config_bcs(i)%bc_kind%ID)
            case (THERMAL_BC_TYPES%DIRICHLET%ID, HYDRAULIC_BC_TYPES%DIRICHLET%ID)
                call set_dirichlet_value(field, config_bcs(i)%boundary_id, value_1)
                write (name, '(a,i0,a,i0)') "essential-", field, "-", config_bcs(i)%boundary_id
                call petscds_add_boundary_essential(dm_handle, trim(name), "Face Sets", &
                                                    int(config_bcs(i)%boundary_id, c_int), field, &
                                                    essential_callback(field), &
                                                    dirichlet_context(field, config_bcs(i)%boundary_id), &
                                                    boundary_index, cerr)
                if (cerr /= 0) error stop "PETSc rejected an essential boundary"

            case (THERMAL_BC_TYPES%FLUX%ID, HYDRAULIC_BC_TYPES%FLUX%ID)
                if (value_1 == 0.0d0) cycle
                if (field /= 0_c_int) then
                    error stop "A non-zero water flux boundary is not supported by the PetscDS path yet"
                end if
                call add_natural_boundary(self, dm_handle, config_bcs(i)%boundary_id, &
                                          BC_FLUX, value_1, 0.0d0)

            case (THERMAL_BC_TYPES%ROBIN%ID)
                call add_natural_boundary(self, dm_handle, config_bcs(i)%boundary_id, &
                                          BC_ROBIN, value_1, value_2)

            case default
                error stop "A boundary condition of this kind is not supported by the PetscDS path yet"
            end select
        end do
    end subroutine add_boundaries

    !> The callback that evaluates the prescribed value of this field.
    function essential_callback(field) result(func)
        implicit none
        integer(c_int), intent(in) :: field
        type(c_funptr) :: func

        if (field == 0_c_int) then
            func = c_funloc(bc_essential_temperature)
        else
            func = c_funloc(bc_essential_pressure)
        end if
    end function essential_callback

    !> Register one natural boundary and the constants its kernels read.
    !>
    !> PETSc keys an auxiliary vector by the label value a boundary was
    !> registered with, and hands that vector to the boundary kernels of that
    !> boundary alone. Each one therefore gets a vector of its own, holding its
    !> condition and its constants.
    subroutine add_natural_boundary(self, dm_handle, boundary_id, kind, value_1, value_2)
        implicit none
        class(type_plexfem), intent(inout) :: self
        type(c_ptr), intent(in) :: dm_handle
        integer(int32), intent(in) :: boundary_id
        integer(int32), intent(in) :: kind
        real(real64), intent(in) :: value_1
        real(real64), intent(in) :: value_2

        PetscErrorCode :: ierr
        integer(c_int) :: boundary_index, cerr
        character(len=32) :: name
        Vec :: boundary_aux

        if (self%num_boundary_aux >= size(self%boundary_aux)) then
            error stop "More natural boundaries than the PetscDS path can hold"
        end if

        write (name, '(a,i0)') "flux-", boundary_id
        call petscds_add_boundary_natural(dm_handle, trim(name), "Face Sets", &
                                          int(boundary_id, c_int), 0_c_int, &
                                          c_funloc(bd_f0_energy), c_funloc(bd_g0_energy), &
                                          boundary_index, cerr)
        if (cerr /= 0) error stop "PETSc rejected a thermal natural boundary"

        call DMCreateLocalVector(self%aux_dm, boundary_aux, ierr)
        if (ierr /= 0) error stop "PETSc could not create a boundary auxiliary vector"

        ! The material still has to be there: the boundary integral reads the
        ! adjacent cell, and its kernels ask the same constitutive model.
        call VecCopy(self%auxiliary, boundary_aux, ierr)
        call fill_field(self%aux_dm, boundary_aux, A_BC_KIND, real(kind, real64))
        call fill_field(self%aux_dm, boundary_aux, A_BC_VALUE_1, value_1)
        call fill_field(self%aux_dm, boundary_aux, A_BC_VALUE_2, value_2)

        call DMSetAuxiliaryVec(self%dm, self%face_sets, int(boundary_id, PETSC_INT_KIND), &
                               0_PETSC_INT_KIND, boundary_aux, ierr)
        if (ierr /= 0) error stop "PETSc could not attach a boundary auxiliary vector"

        self%num_boundary_aux = self%num_boundary_aux + 1
        self%boundary_aux(self%num_boundary_aux) = boundary_aux
    end subroutine add_natural_boundary

    !> Set every degree of freedom of one auxiliary field to the same value.
    subroutine fill_field(dm, vec, field, value)
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
        if (ierr /= 0) error stop "The auxiliary fields have no layout"
        call PetscSectionGetChart(section, p_start, p_end, ierr)

        call VecGetArray(vec, raw, ierr)
        if (ierr /= 0) error stop "PETSc could not open an auxiliary vector"
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
    end subroutine fill_field

    !> Fill the solution with the project's initial temperature.
    !> Fill each field with the value the project starts it at.
    !>
    !> Written through the local section rather than with VecSet, because the
    !> two fields start from different values and share the vector.
    subroutine apply_initial_condition(self)
        implicit none
        class(type_plexfem), intent(inout) :: self

        type(type_config_ic) :: config_thermal, config_hydraulic
        PetscErrorCode :: ierr
        Vec :: local
        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: point, p_start, p_end, offset, num_dofs
        real(real64) :: temperature, pressure

        call input_translator%execute(self%input, IC_TARGETS%THERMAL, config_thermal)
        if (config_thermal%ic_kind%ID /= IC_METHODS%UNIFORM%ID) then
            error stop "Only a uniform initial temperature is supported by the PetscDS path yet"
        end if
        temperature = config_thermal%value

        pressure = 0.0d0
        if (hydraulic_is_active()) then
            call input_translator%execute(self%input, IC_TARGETS%HYDRAULIC, config_hydraulic)
            if (config_hydraulic%ic_kind%ID /= IC_METHODS%UNIFORM%ID) then
                error stop "Only a uniform initial pressure is supported by the PetscDS path yet"
            end if
            pressure = config_hydraulic%value
        end if

        call DMGetLocalVector(self%dm, local, ierr)
        if (ierr == 0) call VecZeroEntries(local, ierr)
        if (ierr /= 0) error stop "PETSc could not take a local vector"

        call DMGetLocalSection(self%dm, section, ierr)
        if (ierr /= 0) error stop "The solve layout is missing"
        call PetscSectionGetChart(section, p_start, p_end, ierr)

        call VecGetArray(local, raw, ierr)
        if (ierr /= 0) error stop "PETSc could not open a local vector"
        do point = p_start, p_end - 1
            call PetscSectionGetFieldDof(section, point, 0_PETSC_INT_KIND, num_dofs, ierr)
            if (ierr /= 0 .or. num_dofs <= 0) cycle
            call PetscSectionGetFieldOffset(section, point, 0_PETSC_INT_KIND, offset, ierr)
            if (ierr == 0) raw(offset + 1) = real(temperature, PETSC_REAL_KIND)

            if (.not. hydraulic_is_active()) cycle
            call PetscSectionGetFieldOffset(section, point, 1_PETSC_INT_KIND, offset, ierr)
            if (ierr == 0) raw(offset + 1) = real(pressure, PETSC_REAL_KIND)
        end do
        call VecRestoreArray(local, raw, ierr)

        call VecZeroEntries(self%solution, ierr)
        if (ierr == 0) call DMLocalToGlobalBegin(self%dm, local, INSERT_VALUES, self%solution, ierr)
        if (ierr == 0) call DMLocalToGlobalEnd(self%dm, local, INSERT_VALUES, self%solution, ierr)
        if (ierr /= 0) error stop "PETSc could not set the initial condition"

        call DMRestoreLocalVector(self%dm, local, ierr)
    end subroutine apply_initial_condition

end submodule plexfem_setup
