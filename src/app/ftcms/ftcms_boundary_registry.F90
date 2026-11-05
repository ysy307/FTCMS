#include <petsc/finclude/petscsys.h>

!> The boundary conditions of `Conditions.json`, as PETSc boundary callbacks.
!>
!> Every condition is evaluated through its strategy's `evaluate`, which returns
!> either a prescribed value or an inflow and its derivative. The engine only
!> decides the role each condition plays in the space of its equation:
!>
!> | space      | prescribed potential (Dirichlet) | flux, Robin, atmospheric, ... |
!> | ---------- | -------------------------------- | ----------------------------- |
!> | `lagrange` | essential on the potential       | natural on the potential row  |
!> | `mixed`    | natural on the flux row          | essential on the flux (normal) |
!>
!> A PetscDS pointwise function carries no typed context, so the managers and
!> the slot table are module state, owned here and nowhere else. One engine per
!> rank may hold them.
!>
!> ## Sign and scale
!>
!> `flux_value` is the inflow per unit area, positive into the domain (heat in
!> W/m2, water in m/s). Each row is scaled by its slot's `row_scale` and each
!> unknown by its `unknown_scale`, so a natural condition contributes
!> \(f_0^\partial=-q_\mathrm{in}\,s_\mathrm{row}\) and
!> \(g_0^\partial=-(\partial q_\mathrm{in}/\partial u)\,s_u\,s_\mathrm{row}\).
module app_ftcms_boundary_registry
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_loc, c_f_pointer
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use :: petscsys
    use :: module_core, only:PHYSICS_TYPES, type_config_bc
    use :: module_boundary, only:type_bc_manager, type_bc_result, BC_DATA_SIZE
    use :: app_ftcms_darcy, only:PLANE_DIMENSION
    use :: app_ftcms_closure, only:type_closure_response, eval_closure_response
    use :: app_ftcms_kernels, only:AUX_BOUNDARY_SLOT, CONST_PRESSURE_SCALE, &
                                   material_slot_of, note_kernel_fault, point_state
    implicit none
    private

    public :: install_boundary_conditions
    public :: clear_boundary_registry
    public :: add_boundary_slot
    public :: boundary_slot_context
    public :: set_boundary_time
    public :: update_boundary_data
    public :: boundary_condition_count
    public :: boundary_condition_id
    public :: boundary_manager
    public :: bd_f0_natural, bd_g0_natural
    public :: bd_f0_flux_potential
    public :: bd_f0_seepage
    public :: boundary_data
    public :: bc_essential_potential
    public :: bc_essential_normal_flux
    public :: set_slot_facets
    public :: bc_essential_zero
    public :: MAX_BOUNDARY_SLOTS

    !> How many (equation, boundary) pairs one run may register.
    integer(int32), parameter :: MAX_BOUNDARY_SLOTS = 64

    !> One registered boundary condition.
    type, bind(C) :: type_boundary_slot
        !> This slot's own number, so a callback context leads back to it.
        integer(c_int) :: slot = 0
        !> `PHYSICS_TYPES` id of the equation the condition belongs to.
        integer(c_int) :: physics_id = 0
        !> Index of the condition in that equation's manager.
        integer(c_int) :: condition = 0
        !> Zero-based field number the condition reads its unknown from.
        integer(c_int) :: field = 0
        !> Physical value of the unknown per unit of it.
        real(real64) :: unknown_scale = 1.0d0
        !> What the physical inflow is multiplied by in the residual row.
        real(real64) :: row_scale = 1.0d0
    end type type_boundary_slot

    !> The facets of a slot whose condition needs their outward normals.
    type :: type_slot_facets
        real(real64), allocatable :: starts(:, :)
        real(real64), allocatable :: ends(:, :)
        real(real64), allocatable :: normals(:, :)
    end type type_slot_facets

    type(type_bc_manager), save, target :: managers(PHYSICS_TYPES%NUM_ID)
    type(type_slot_facets), save :: slot_facets(MAX_BOUNDARY_SLOTS)
    type(type_boundary_slot), save, target :: slots(MAX_BOUNDARY_SLOTS)
    integer(int32), save :: num_slots = 0
    !> Time of the state being solved for, s; set before every solve because
    !> PETSc hands SNES callbacks no time of their own.
    real(real64), save :: current_time = 0.0d0

contains

    !> Build the strategy of every condition of one equation.
    subroutine install_boundary_conditions(physics_id, configs)
        implicit none
        integer(int32), intent(in) :: physics_id
        type(type_config_bc), intent(in) :: configs(:)

        call managers(physics_id)%destroy()
        call managers(physics_id)%initialize(configs)
    end subroutine install_boundary_conditions

    !> Forget every condition and slot.
    subroutine clear_boundary_registry()
        implicit none

        integer(int32) :: id

        do id = 1, PHYSICS_TYPES%NUM_ID
            call managers(id)%destroy()
        end do
        do id = 1, MAX_BOUNDARY_SLOTS
            slot_facets(id) = type_slot_facets()
        end do
        num_slots = 0
        current_time = 0.0d0
    end subroutine clear_boundary_registry

    !> How many conditions one equation carries.
    function boundary_condition_count(physics_id) result(count)
        implicit none
        integer(int32), intent(in) :: physics_id
        integer(int32) :: count

        count = 0
        call managers(physics_id)%get_num_boundaries(count)
    end function boundary_condition_count

    !> The "Face Sets" value condition `condition` of one equation is keyed by.
    function boundary_condition_id(physics_id, condition) result(boundary_id)
        implicit none
        integer(int32), intent(in) :: physics_id
        integer(int32), intent(in) :: condition
        integer(int32) :: boundary_id

        boundary_id = managers(physics_id)%get_boundary_id(condition)
    end function boundary_condition_id

    !> The strategies of one equation, for the assimilation controller to update.
    function boundary_manager(physics_id) result(manager)
        implicit none
        integer(int32), intent(in) :: physics_id
        type(type_bc_manager), pointer :: manager

        manager => managers(physics_id)
    end function boundary_manager

    !> Register one condition and return its slot.
    !>
    !> Fails by `error stop` once `MAX_BOUNDARY_SLOTS` slots are taken.
    function add_boundary_slot(physics_id, condition, field, unknown_scale, row_scale) result(slot)
        implicit none
        integer(int32), intent(in) :: physics_id
        integer(int32), intent(in) :: condition
        integer(int32), intent(in) :: field
        real(real64), intent(in) :: unknown_scale
        real(real64), intent(in) :: row_scale
        integer(int32) :: slot

        if (num_slots >= MAX_BOUNDARY_SLOTS) error stop "More boundary conditions than the registry can hold"
        num_slots = num_slots + 1
        slot = num_slots
        slots(slot) = type_boundary_slot(int(slot, c_int), int(physics_id, c_int), int(condition, c_int), int(field, c_int), &
                                         unknown_scale, row_scale)
    end function add_boundary_slot

    !> The context an essential callback receives for slot `slot`.
    function boundary_slot_context(slot) result(ctx)
        implicit none
        integer(int32), intent(in) :: slot
        type(c_ptr) :: ctx

        ctx = c_loc(slots(slot))
    end function boundary_slot_context

    !> Set the time every condition is evaluated at, s.
    subroutine set_boundary_time(time)
        implicit none
        real(real64), intent(in) :: time

        current_time = time
    end subroutine set_boundary_time

    !> Replace the data of the conditions keyed by `boundary_id` of one equation.
    !>
    !> The assimilation controller uses this to hand its analysis to the
    !> dynamic providers.
    subroutine update_boundary_data(physics_id, boundary_id, values)
        implicit none
        integer(int32), intent(in) :: physics_id
        integer(int32), intent(in) :: boundary_id
        real(real64), intent(in) :: values(:)

        call managers(physics_id)%update_bc_data(boundary_id, values)
    end subroutine update_boundary_data

    !> Evaluate slot `slot` at the physical unknown `value`.
    subroutine evaluate_slot(slot, value, result)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64), intent(in) :: value
        type(type_bc_result), intent(inout) :: result

        call managers(slots(slot)%physics_id)%evaluate(slots(slot)%condition, current_time, &
                                                       value, result)
    end subroutine evaluate_slot

    !> The provider data of slot `slot` at the time being solved for.
    subroutine boundary_data(slot, values)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64), intent(inout) :: values(BC_DATA_SIZE)

        call managers(slots(slot)%physics_id)%get_data(slots(slot)%condition, current_time, values)
    end subroutine boundary_data

    !> The slot a boundary kernel is evaluating, from its auxiliary field.
    pure function slot_of(a, a_off) result(slot)
        implicit none
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        integer(int32) :: slot

        slot = nint(a(a_off(AUX_BOUNDARY_SLOT) + 1), int32)
    end function slot_of

    !> Natural condition on a potential row, \(f_0^\partial=-q_\mathrm{in}(u)\,s_\mathrm{row}\).
    subroutine bd_f0_natural( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, normal, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_bd_f0_natural")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscReal, intent(in) :: normal(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f0(*)

        type(type_bc_result) :: result
        integer(int32) :: slot

        slot = slot_of(a, a_off)
        associate (entry => slots(slot))
            call evaluate_slot(slot, entry%unknown_scale*u(u_off(entry%field + 1) + 1), result)
            f0(1) = -result%flux_value*entry%row_scale
        end associate
    end subroutine bd_f0_natural

    !> \(\partial f_0^\partial/\partial u=-(\partial q_\mathrm{in}/\partial u)\,s_u\,s_\mathrm{row}\).
    subroutine bd_g0_natural( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, normal, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_bd_g0_natural")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time
        PetscReal, value :: u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscReal, intent(in) :: normal(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g0(*)

        type(type_bc_result) :: result
        integer(int32) :: slot

        slot = slot_of(a, a_off)
        associate (entry => slots(slot))
            call evaluate_slot(slot, entry%unknown_scale*u(u_off(entry%field + 1) + 1), result)
            g0(1) = -result%flux_derivative*entry%unknown_scale*entry%row_scale
        end associate
    end subroutine bd_g0_natural

    !> The potential flux of a seepage face, \(f_0^\partial=q_\mathrm{pot}\,s_\mathrm{row}\).
    !>
    !> values(1) is the net evaporation \(q_\mathrm{pot}\) (positive out of the
    !> soil). The face is held between values(2) and values(3) by the variable
    !> bounds of the solve, not here: the flux is the same at every state, so
    !> the residual stays smooth and the bound alone switches the face.
    subroutine bd_f0_seepage( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, normal, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_bd_f0_seepage")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscReal, intent(in) :: normal(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f0(*)

        real(real64) :: values(BC_DATA_SIZE)
        integer(int32) :: slot

        slot = slot_of(a, a_off)
        call boundary_data(slot, values)
        f0(1) = values(1)*slots(slot)%row_scale
    end subroutine bd_f0_seepage

    !> A prescribed pore pressure on the flux row of a mixed hydraulic field.
    !>
    !> \(f_0^\partial=p_l(T,p_D)\,\mathbf n/p_\mathrm{scale}\): the boundary term
    !> of \(-(p_l,\nabla\cdot\mathbf v)\) integrated by parts, with the liquid
    !> pressure the closure gives at the boundary temperature and the prescribed
    !> pore pressure. Its tangent with respect to \(T\) is not registered: the
    !> natural-boundary shim registers same-field Jacobians only, so the Newton
    !> iteration sees this term as lagged in \(T\).
    subroutine bd_f0_flux_potential( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, normal, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_bd_f0_flux_potential")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscReal, intent(in) :: normal(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f0(*)

        type(type_bc_result) :: result
        type(type_closure_response) :: response
        real(real64) :: temperature, pressure
        integer(int32) :: slot, material, d
        logical :: ok

        slot = slot_of(a, a_off)
        call evaluate_slot(slot, 0.0d0, result)
        call point_state(u, u_off, constants, temperature, pressure)
        material = material_slot_of(a, a_off)
        ok = .false.
        call eval_closure_response(material, temperature, result%prescribed_value, response, ok)
        if (.not. ok) then
            call note_kernel_fault(material, temperature, result%prescribed_value)
            do d = 1, PLANE_DIMENSION
                f0(d) = ieee_value(1.0d0, ieee_quiet_nan)
            end do
            return
        end if
        do d = 1, PLANE_DIMENSION
            f0(d) = response%liquid_pressure*normal(d) &
                    /constants(CONST_PRESSURE_SCALE)
        end do
    end subroutine bd_f0_flux_potential

    !> Essential value of a potential: the prescribed value in units of the unknown.
    function bc_essential_potential(dim, time, x, num_components, u, ctx) &
        bind(C, name="ftcms_bc_essential_potential") result(ierr)
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        type(type_boundary_slot), pointer :: entry
        type(type_bc_result) :: result
        integer(int32) :: slot

        call c_f_pointer(ctx, entry)
        slot = entry%slot
        call evaluate_slot(slot, 0.0d0, result)
        u(1) = result%prescribed_value/entry%unknown_scale
        ierr = 0
    end function bc_essential_potential

    !> Essential normal flux of an H(div) field, \(\mathbf q=-q_\mathrm{in}\,\mathbf n\).
    !>
    !> Only the normal component reaches the facet degrees of freedom, so the
    !> tangential part of the value is irrelevant. The outward normal at `x`
    !> is that of the slot's facet nearest to `x`, from the table
    !> `set_slot_facets` holds: PETSc's point callbacks carry no normal.
    function bc_essential_normal_flux(dim, time, x, num_components, u, ctx) &
        bind(C, name="ftcms_bc_essential_normal_flux") result(ierr)
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        type(type_boundary_slot), pointer :: entry
        type(type_bc_result) :: result
        real(real64) :: outward(PLANE_DIMENSION)
        integer(int32) :: slot, d

        call c_f_pointer(ctx, entry)
        slot = entry%slot
        call evaluate_slot(slot, 0.0d0, result)
        outward = facet_normal_at(slot, [(x(d), d=1, PLANE_DIMENSION)])
        do d = 1, PLANE_DIMENSION
            u(d) = -result%flux_value*outward(d)/entry%unknown_scale
        end do
        ierr = 0
    end function bc_essential_normal_flux

    !> Keep the facets of one slot: their end points and outward unit normals.
    subroutine set_slot_facets(slot, starts, ends, normals)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64), intent(in) :: starts(:, :), ends(:, :), normals(:, :)

        slot_facets(slot)%starts = starts
        slot_facets(slot)%ends = ends
        slot_facets(slot)%normals = normals
    end subroutine set_slot_facets

    !> The outward normal of the slot's facet closest to `point`.
    pure function facet_normal_at(slot, point) result(normal)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64), intent(in) :: point(PLANE_DIMENSION)
        real(real64) :: normal(PLANE_DIMENSION)

        real(real64) :: edge(PLANE_DIMENSION), offset(PLANE_DIMENSION), along, distance, nearest
        integer(int32) :: k

        normal = 0.0d0
        nearest = huge(1.0d0)
        associate (facets => slot_facets(slot))
            if (.not. allocated(facets%normals)) return
            do k = 1, size(facets%normals, 2)
                edge = facets%ends(:, k) - facets%starts(:, k)
                along = min(max(dot_product(point - facets%starts(:, k), edge)/dot_product(edge, edge), &
                                0.0d0), 1.0d0)
                offset = point - (facets%starts(:, k) + along*edge)
                distance = dot_product(offset, offset)
                if (distance < nearest) then
                    nearest = distance
                    normal = facets%normals(:, k)
                end if
            end do
        end associate
    end function facet_normal_at

    !> Essential zero of every component: the closed default of a mixed field.
    function bc_essential_zero(dim, time, x, num_components, u, ctx) &
        bind(C, name="ftcms_bc_essential_zero") result(ierr)
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        integer(int32) :: c

        do c = 1, num_components
            u(c) = 0.0d0
        end do
        ierr = 0
    end function bc_essential_zero

end module app_ftcms_boundary_registry
