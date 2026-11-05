#include <petsc/finclude/petscsys.h>

!> The heat and water weak forms, as PetscDS pointwise functions.
!>
!> \(T\) in \(H^1\); the water either as the mixed pair \(\mathbf q_l\) in
!> \(H(\mathrm{div})\) and \(p_w\) in \(L^2\), or as \(p_w\) in \(H^1\)
!> (`app_ftcms_kernels_lagrange`):
!>
!> \[
!>   R_T(w)=\bigl(\dot U(T,p_w),w\bigr)+(\lambda\nabla T,\nabla w),
!> \]
!> \[
!>   R_q(\mathbf v)=(\mathbf\Lambda^{-1}(T,p_w)\mathbf q_l-\rho_l\mathbf g,\mathbf v)
!>   -(p_l(T,p_w),\nabla\cdot\mathbf v),
!> \]
!> \[
!>   R_p(\omega)=\bigl(\dot W(T,p_w)+\nabla\cdot\mathbf q_l-s,\omega\bigr),
!> \]
!>
!> with \(\dot X=(c_0X+c_1X^n+c_2X^{n-1})/\Delta t\), the BDF rate whose
!> weights are run constants. Field numbers are run constants too, so the same
!> functions serve every layout the engine's field table builds; an equation
!> that is not solved has the field number -1 and a fixed potential.
!>
!> ## What is coupled, and where
!>
!> Every constitutive value comes from `app_ftcms_closure`, which is
!> the repository's own governing models. This module states no retention
!> curve, no freezing curve, no impedance law, no enthalpy and no mixing rule.
!>
!> Cryosuction enters once, as \(p_l=-\psi_{\mathrm{eff}}(T,p_w)\) in the Darcy
!> row, and nothing adds a second temperature-driven water flux beside it.
!>
!> ## The history
!>
!> \(U^n\) and \(W^n\) are not stored. The accepted temperature and pressure
!> are, in auxiliary fields with the same discretisations as the unknowns, and
!> the kernels re-evaluate the closure at them. The storage rate is then a
!> difference of one function of the state rather than of two projections of
!> it.
!>
!> Heat carried by the liquid is added in `app_ftcms_kernels_advection`.
module app_ftcms_kernels
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_loc, c_f_pointer
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use :: petscsys
    use :: app_ftcms_darcy, only:PLANE_DIMENSION
    use :: app_ftcms_darcy, only:resistance_of
    use :: app_ftcms_closure, only:type_closure_response, type_closure_tangents, &
                                            eval_closure_response, eval_closure_tangents, eval_closure_storage
    implicit none
    private

    public :: f0_temperature, f1_temperature
    public :: f0_flux, f1_flux
    public :: f0_conservation
    public :: g0_temperature_temperature, g2_temperature_temperature
    public :: g3_temperature_temperature
    public :: g0_temperature_pressure, g2_temperature_pressure
    public :: g0_flux_flux, g0_flux_temperature, g2_flux_temperature
    public :: g0_flux_pressure, g2_flux_pressure
    public :: g1_pressure_flux, g0_pressure_pressure, g0_pressure_temperature
    public :: clear_kernel_fault, kernel_fault_recorded, kernel_fault_report
    public :: NUM_AUX_FIELDS
    public :: AUX_INVERSE_MOBILITY, AUX_MATERIAL, AUX_SOURCE, AUX_TEMPERATURE_OLD, AUX_PRESSURE_OLD
    public :: AUX_BOUNDARY_SLOT
    public :: CONST_DT, CONST_FLUX_SCALE, CONST_PRESSURE_SCALE, CONST_LENGTH_SCALE, CONST_DENSITY
    public :: CONST_GRAVITY_X, CONST_GRAVITY_Z, CONST_TEMPERATURE_SCALE, CONST_THERMAL_SCALE, CONST_PRESSURE_FIELD
    public :: CONST_TEMPERATURE_FIELD, CONST_FLUX_FIELD, CONST_FIXED_TEMPERATURE, CONST_FIXED_PRESSURE
    public :: field_number
    public :: CONST_BDF_CURRENT, CONST_BDF_ACCEPTED, CONST_BDF_OLDER
    public :: AUX_TEMPERATURE_OLDER, AUX_PRESSURE_OLDER
    public :: storage_rate_of, STORE_ENTHALPY, STORE_WATER
    public :: NUM_CONSTANTS
    public :: point_state, point_state_old, material_slot_of, inverse_mobility_of
    public :: point_tangents, note_kernel_fault


    ! Auxiliary fields, in the order PETSc offsets them.
    integer(int32), parameter :: AUX_INVERSE_MOBILITY = 1 !< Saturated inverse mobility, row major, P0.
    integer(int32), parameter :: AUX_MATERIAL = 2         !< Slot in the frozen closure table, P0.
    integer(int32), parameter :: AUX_SOURCE = 3           !< Water source, 1/s, P0.
    integer(int32), parameter :: AUX_TEMPERATURE_OLD = 4  !< Accepted temperature, C, P1.
    integer(int32), parameter :: AUX_PRESSURE_OLD = 5     !< Accepted pore pressure, Pa, P0.
    ! Read by the boundary kernels alone. A natural boundary's pointwise
    ! kernels take no context, so what distinguishes one boundary from another
    ! is the auxiliary vector PETSc looks up for that boundary's label value;
    ! it carries the boundary's slot in `app_ftcms_boundary_registry`. The
    ! volume kernels never read it and the default auxiliary vector leaves it
    ! zero.
    integer(int32), parameter :: AUX_BOUNDARY_SLOT = 6   !< Slot of the boundary condition, P0.
    ! The state accepted one step before, for a second-order rate.
    integer(int32), parameter :: AUX_TEMPERATURE_OLDER = 7 !< Temperature before the accepted one, C.
    integer(int32), parameter :: AUX_PRESSURE_OLDER = 8    !< Pore pressure before the accepted one, Pa.
    integer(int32), parameter :: NUM_AUX_FIELDS = 8

    ! Run constants, in the order PetscDSSetConstants receives them.
    integer(int32), parameter :: CONST_DT = 1                !< Time step, s.
    integer(int32), parameter :: CONST_FLUX_SCALE = 2        !< Flux unit of the unknown, m/s.
    integer(int32), parameter :: CONST_PRESSURE_SCALE = 3    !< Pressure unit of the unknown, Pa.
    integer(int32), parameter :: CONST_LENGTH_SCALE = 4      !< Length the conservation row is scaled by, m.
    integer(int32), parameter :: CONST_DENSITY = 5           !< Reference liquid density, kg/m3.
    integer(int32), parameter :: CONST_GRAVITY_X = 6         !< Gravity along intrinsic x, m/s2.
    integer(int32), parameter :: CONST_GRAVITY_Z = 7         !< Gravity along intrinsic z, m/s2.
    integer(int32), parameter :: CONST_TEMPERATURE_SCALE = 8 !< Temperature unit of the unknown, K.
    integer(int32), parameter :: CONST_THERMAL_SCALE = 9     !< The heat row's divisor, W/m3.
    ! Field numbers come from the engine's field table, so every kernel serves
    ! every layout the table builds. A field that the run does not solve for
    ! has the number -1, and its potential is then the fixed value below.
    integer(int32), parameter :: CONST_PRESSURE_FIELD = 10   !< Zero-based field number of the pressure.
    integer(int32), parameter :: CONST_TEMPERATURE_FIELD = 11 !< Zero-based field number of the temperature.
    integer(int32), parameter :: CONST_FLUX_FIELD = 12       !< Zero-based field number of the liquid flux.
    integer(int32), parameter :: CONST_FIXED_TEMPERATURE = 13 !< Temperature of a run without a heat equation, C.
    integer(int32), parameter :: CONST_FIXED_PRESSURE = 14   !< Pore pressure of a run without a water equation, Pa.
    ! The BDF weights of the current, accepted and older storage: the rate is
    ! (c0 X + c1 X^n + c2 X^(n-1))/dt, so backward Euler is (1, -1, 0).
    integer(int32), parameter :: CONST_BDF_CURRENT = 15
    integer(int32), parameter :: CONST_BDF_ACCEPTED = 16
    integer(int32), parameter :: CONST_BDF_OLDER = 17
    integer(int32), parameter :: NUM_CONSTANTS = 17


    !> Which stored quantity `storage_rate_of` differences.
    integer(int32), parameter :: STORE_ENTHALPY = 1
    integer(int32), parameter :: STORE_WATER = 2

    ! What the last unusable closure evaluation was. A pointwise kernel has
    ! nowhere to return an error to, so it records the fault and writes a
    ! nonfinite residual, which is what refuses the step.
    logical, save :: kernel_fault = .false.
    integer(int32), save :: kernel_fault_slot = 0
    real(real64), save :: kernel_fault_temperature = 0.0d0
    real(real64), save :: kernel_fault_pressure = 0.0d0

contains

    !> The cell's saturated inverse mobility, row major, from its auxiliary field.
    !>
    !> This layout's field number, not the isothermal one: the two auxiliary
    !> sections carry different fields in different places, and reading one
    !> through the other's offsets is how a kernel silently gets a zero.
    pure function inverse_mobility_of(a, a_off) result(inverse_mobility)
        implicit none
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        real(real64) :: inverse_mobility(PLANE_DIMENSION*PLANE_DIMENSION)

        integer(int32) :: k

        do k = 1, PLANE_DIMENSION*PLANE_DIMENSION
            inverse_mobility(k) = a(a_off(AUX_INVERSE_MOBILITY) + k)
        end do
    end function inverse_mobility_of

    !> The cell's slot in the frozen closure table, from its auxiliary field.
    pure function material_slot_of(a, a_off) result(slot)
        implicit none
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        integer(int32) :: slot

        slot = nint(a(a_off(AUX_MATERIAL) + 1), int32)
    end function material_slot_of


    !> Forget any closure fault a previous trial recorded.
    subroutine clear_kernel_fault()
        implicit none

        kernel_fault = .false.
        kernel_fault_slot = 0
        kernel_fault_temperature = 0.0d0
        kernel_fault_pressure = 0.0d0
    end subroutine clear_kernel_fault

    !> Whether a kernel could not evaluate the closure since the last clear.
    logical function kernel_fault_recorded()
        implicit none

        kernel_fault_recorded = kernel_fault
    end function kernel_fault_recorded

    !> What the recorded fault was.
    subroutine kernel_fault_report(slot, temperature, pressure)
        implicit none
        !> The material slot of the offending cell, zero when none was recorded.
        integer(int32), intent(inout) :: slot
        !> The state it was evaluated at.
        real(real64), intent(inout) :: temperature, pressure

        slot = kernel_fault_slot
        temperature = kernel_fault_temperature
        pressure = kernel_fault_pressure
    end subroutine kernel_fault_report

    !> Record that a kernel asked the closure for a state it cannot answer.
    subroutine note_kernel_fault(slot, temperature, pressure)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64), intent(in) :: temperature, pressure

        if (kernel_fault) return
        kernel_fault = .true.
        kernel_fault_slot = slot
        kernel_fault_temperature = temperature
        kernel_fault_pressure = pressure
    end subroutine note_kernel_fault

    !> The physical state of this quadrature point, from the unknowns.
    !>
    !> The pressure is read from the field `CONST_PRESSURE_FIELD` names, so the
    !> same accessor serves the mixed layout and the primal one.
    pure subroutine point_state(u, u_off, constants, temperature, pressure)
        implicit none
        PetscScalar, intent(in) :: u(*)
        PetscInt, intent(in) :: u_off(*)
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: temperature, pressure

        integer(int32) :: field

        field = field_number(constants, CONST_TEMPERATURE_FIELD)
        if (field >= 0) then
            temperature = constants(CONST_TEMPERATURE_SCALE)*u(u_off(field + 1) + 1)
        else
            temperature = constants(CONST_FIXED_TEMPERATURE)
        end if
        field = field_number(constants, CONST_PRESSURE_FIELD)
        if (field >= 0) then
            pressure = constants(CONST_PRESSURE_SCALE)*u(u_off(field + 1) + 1)
        else
            pressure = constants(CONST_FIXED_PRESSURE)
        end if
    end subroutine point_state

    !> The zero-based field number a run constant names, -1 when absent.
    pure function field_number(constants, which) result(field)
        implicit none
        PetscScalar, intent(in) :: constants(*)
        !> `CONST_TEMPERATURE_FIELD`, `CONST_FLUX_FIELD` or `CONST_PRESSURE_FIELD`.
        integer(int32), intent(in) :: which
        integer(int32) :: field

        field = nint(constants(which), int32)
    end function field_number

    !> The accepted state of this quadrature point, from the auxiliary fields.
    pure subroutine point_state_old(a, a_off, temperature, pressure)
        implicit none
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        real(real64), intent(inout) :: temperature, pressure

        temperature = a(a_off(AUX_TEMPERATURE_OLD) + 1)
        pressure = a(a_off(AUX_PRESSURE_OLD) + 1)
    end subroutine point_state_old

    !> \((c_0X+c_1X^n+c_2X^{n-1})/\Delta t\) of one stored quantity at this point.
    !>
    !> `current` is \(X\) at the unknown state; the history is evaluated at
    !> the accepted states in the auxiliary fields, so the rate is a difference
    !> of one function of the state. `ok` is false when the closure cannot
    !> answer at a history state.
    subroutine storage_rate_of(a, a_off, constants, slot, current, quantity, rate, ok)
        implicit none
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        PetscScalar, intent(in) :: constants(*)
        integer(int32), intent(in) :: slot
        real(real64), intent(in) :: current
        !> `STORE_ENTHALPY` or `STORE_WATER`.
        integer(int32), intent(in) :: quantity
        real(real64), intent(inout) :: rate
        logical, intent(inout) :: ok

        real(real64) :: temperature, pressure, older_weight, water, enthalpy

        water = 0.0d0
        enthalpy = 0.0d0
        rate = constants(CONST_BDF_CURRENT)*current
        call point_state_old(a, a_off, temperature, pressure)
        call eval_closure_storage(slot, temperature, pressure, water, enthalpy, ok)
        if (.not. ok) return
        rate = rate + constants(CONST_BDF_ACCEPTED)*merge(enthalpy, water, quantity == STORE_ENTHALPY)
        older_weight = constants(CONST_BDF_OLDER)
        if (older_weight /= 0.0d0) then
            temperature = a(a_off(AUX_TEMPERATURE_OLDER) + 1)
            pressure = a(a_off(AUX_PRESSURE_OLDER) + 1)
            call eval_closure_storage(slot, temperature, pressure, water, enthalpy, ok)
            if (.not. ok) return
            rate = rate + older_weight*merge(enthalpy, water, quantity == STORE_ENTHALPY)
        end if
        rate = rate/constants(CONST_DT)
    end subroutine storage_rate_of

    !> Heat storage, \(f_0=(U(T,p_w)-U^n)/(\Delta t\,R_T)\).
    subroutine f0_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_f0_temperature")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f0(*)

        type(type_closure_response) :: now, accepted
        real(real64) :: temperature, pressure, temperature_old, pressure_old
        real(real64) :: dt, thermal_scale, rate
        integer(int32) :: slot
        logical :: ok

        dt = constants(CONST_DT)
        thermal_scale = constants(CONST_THERMAL_SCALE)
        slot = material_slot_of(a, a_off)
        call point_state(u, u_off, constants, temperature, pressure)
        call point_state_old(a, a_off, temperature_old, pressure_old)

        ok = .false.
        call eval_closure_response(slot, temperature, pressure, now, ok)
        if (ok) call storage_rate_of(a, a_off, constants, slot, now%enthalpy, STORE_ENTHALPY, rate, ok)
        if (.not. ok) then
            call note_kernel_fault(slot, temperature, pressure)
            f0(1) = ieee_value(1.0d0, ieee_quiet_nan)
            return
        end if

        f0(1) = rate/thermal_scale
    end subroutine f0_temperature

    !> Conduction, \(f_1=\lambda T_{\mathrm{scale}}\nabla u_T/R_T\).
    subroutine f1_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f1 &
        ) bind(C, name="ftcms_f1_temperature")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f1(*)

        type(type_closure_response) :: now
        real(real64) :: temperature, pressure, weight
        integer(int32) :: d, slot
        logical :: ok

        slot = material_slot_of(a, a_off)
        call point_state(u, u_off, constants, temperature, pressure)

        ok = .false.
        call eval_closure_response(slot, temperature, pressure, now, ok)
        if (.not. ok) then
            call note_kernel_fault(slot, temperature, pressure)
            do d = 1, PLANE_DIMENSION
                f1(d) = ieee_value(1.0d0, ieee_quiet_nan)
            end do
            return
        end if

        weight = now%thermal_conductivity*constants(CONST_TEMPERATURE_SCALE) &
                 /constants(CONST_THERMAL_SCALE)
        do d = 1, PLANE_DIMENSION
            f1(d) = weight*u_x(u_off_x(field_number(constants, CONST_TEMPERATURE_FIELD) + 1) + d)
        end do
    end subroutine f1_temperature

    !> Darcy's law, \(f_0=(\mathbf\Lambda^{-1}(T,p_w)\mathbf q_l-\rho_l\mathbf g)/p_{\mathrm{scale}}\).
    subroutine f0_flux( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_f0_flux")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f0(*)

        type(type_closure_response) :: now
        real(real64) :: inverse_mobility(PLANE_DIMENSION*PLANE_DIMENSION)
        real(real64) :: flux(PLANE_DIMENSION), gravity(PLANE_DIMENSION), resistance(PLANE_DIMENSION)
        real(real64) :: temperature, pressure, flux_scale, pressure_scale, density
        integer(int32) :: i, slot
        logical :: ok

        flux_scale = constants(CONST_FLUX_SCALE)
        pressure_scale = constants(CONST_PRESSURE_SCALE)
        density = constants(CONST_DENSITY)
        gravity(1) = constants(CONST_GRAVITY_X)
        gravity(2) = constants(CONST_GRAVITY_Z)

        slot = material_slot_of(a, a_off)
        call point_state(u, u_off, constants, temperature, pressure)

        ok = .false.
        call eval_closure_response(slot, temperature, pressure, now, ok)
        if (.not. ok) then
            call note_kernel_fault(slot, temperature, pressure)
            do i = 1, PLANE_DIMENSION
                f0(i) = ieee_value(1.0d0, ieee_quiet_nan)
            end do
            return
        end if

        inverse_mobility = inverse_mobility_of(a, a_off)/now%relative_conductivity
        do i = 1, PLANE_DIMENSION
            flux(i) = flux_scale*u(u_off(field_number(constants, CONST_FLUX_FIELD) + 1) + i)
        end do
        resistance = resistance_of(inverse_mobility, flux)

        do i = 1, PLANE_DIMENSION
            f0(i) = (resistance(i) - density*gravity(i))/pressure_scale
        end do
    end subroutine f0_flux

    !> The driver of Darcy's law, \(f_1=-p_l(T,p_w)\mathbf I/p_{\mathrm{scale}}\).
    !>
    !> The liquid pressure, not the solved \(p_w\): that difference is the
    !> cryosuction, and evaluating it here is the only place it enters.
    subroutine f1_flux( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f1 &
        ) bind(C, name="ftcms_f1_flux")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f1(*)

        type(type_closure_response) :: now
        real(real64) :: temperature, pressure, driver
        integer(int32) :: c, d, slot
        logical :: ok

        slot = material_slot_of(a, a_off)
        call point_state(u, u_off, constants, temperature, pressure)

        ok = .false.
        call eval_closure_response(slot, temperature, pressure, now, ok)
        if (.not. ok) then
            call note_kernel_fault(slot, temperature, pressure)
            do c = 1, PLANE_DIMENSION*PLANE_DIMENSION
                f1(c) = ieee_value(1.0d0, ieee_quiet_nan)
            end do
            return
        end if

        driver = now%liquid_pressure/constants(CONST_PRESSURE_SCALE)
        do c = 1, PLANE_DIMENSION
            do d = 1, PLANE_DIMENSION
                if (c == d) then
                    f1((c - 1)*PLANE_DIMENSION + d) = -driver
                else
                    f1((c - 1)*PLANE_DIMENSION + d) = 0.0d0
                end if
            end do
        end do
    end subroutine f1_flux

    !> Water conservation, \(f_0=((W-W^n)/\Delta t+\nabla\cdot\mathbf q_l-s)L/q_{\mathrm{scale}}\).
    subroutine f0_conservation( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_f0_conservation")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f0(*)

        type(type_closure_response) :: now, accepted
        real(real64) :: temperature, pressure, temperature_old, pressure_old
        real(real64) :: flux_scale, length_scale, dt, divergence, source, storage_rate
        integer(int32) :: c, slot
        logical :: ok

        flux_scale = constants(CONST_FLUX_SCALE)
        length_scale = constants(CONST_LENGTH_SCALE)
        dt = constants(CONST_DT)

        ! PETSc stores the gradient component-major: du_c/dx_d at c*dim + d.
        divergence = 0.0d0
        do c = 1, PLANE_DIMENSION
            divergence = divergence + u_x(u_off_x(field_number(constants, CONST_FLUX_FIELD) + 1) &
                                          + (c - 1)*PLANE_DIMENSION + c)
        end do

        source = a(a_off(AUX_SOURCE) + 1)
        slot = material_slot_of(a, a_off)
        call point_state(u, u_off, constants, temperature, pressure)
        call point_state_old(a, a_off, temperature_old, pressure_old)

        storage_rate = 0.0d0
        if (dt > 0.0d0) then
            ok = .false.
            call eval_closure_response(slot, temperature, pressure, now, ok)
            if (ok) call storage_rate_of(a, a_off, constants, slot, now%stored_water, STORE_WATER, storage_rate, ok)
            if (.not. ok) then
                call note_kernel_fault(slot, temperature, pressure)
                f0(1) = ieee_value(1.0d0, ieee_quiet_nan)
                return
            end if
        end if

        f0(1) = (storage_rate + flux_scale*divergence - source)*length_scale/flux_scale
    end subroutine f0_conservation

    !> The tangents of this quadrature point, or a recorded fault.
    subroutine point_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        implicit none
        PetscScalar, intent(in) :: u(*)
        PetscInt, intent(in) :: u_off(*)
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        PetscScalar, intent(in) :: constants(*)
        integer(int32), intent(inout) :: slot
        type(type_closure_tangents), intent(inout) :: tangents
        logical, intent(inout) :: ok

        real(real64) :: temperature, pressure

        slot = material_slot_of(a, a_off)
        call point_state(u, u_off, constants, temperature, pressure)
        ok = .false.
        call eval_closure_tangents(slot, temperature, pressure, tangents, ok)
        if (.not. ok) call note_kernel_fault(slot, temperature, pressure)
    end subroutine point_tangents

    !> \(\partial f_0^{T}/\partial T=(\partial U/\partial T)T_{\mathrm{scale}}/(\Delta t\,R_T)\).
    subroutine g0_temperature_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_g0_TT")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time, u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g0(*)

        type(type_closure_tangents) :: tangents
        integer(int32) :: slot
        logical :: ok

        g0(1) = 0.0d0
        if (constants(CONST_DT) <= 0.0d0) return
        call point_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        g0(1) = tangents%denthalpy_dtemperature*constants(CONST_TEMPERATURE_SCALE) &
                *constants(CONST_BDF_CURRENT) &
                /(constants(CONST_DT)*constants(CONST_THERMAL_SCALE))
    end subroutine g0_temperature_temperature

    !> \(\partial f_1^{T}/\partial T=(\partial\lambda/\partial T)T_{\mathrm{scale}}^2\nabla u_T/R_T\).
    subroutine g2_temperature_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_g2_TT")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time, u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g2(*)

        type(type_closure_tangents) :: tangents
        real(real64) :: weight
        integer(int32) :: d, slot
        logical :: ok

        do d = 1, PLANE_DIMENSION
            g2(d) = 0.0d0
        end do
        call point_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        weight = tangents%dthermal_dtemperature &
                 *constants(CONST_TEMPERATURE_SCALE)**2 &
                 /constants(CONST_THERMAL_SCALE)
        do d = 1, PLANE_DIMENSION
            g2(d) = weight*u_x(u_off_x(field_number(constants, CONST_TEMPERATURE_FIELD) + 1) + d)
        end do
    end subroutine g2_temperature_temperature

    !> \(\partial f_1^{T}/\partial\nabla T=\lambda T_{\mathrm{scale}}\mathbf I/R_T\).
    subroutine g3_temperature_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g3 &
        ) bind(C, name="ftcms_g3_TT")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time, u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g3(*)

        type(type_closure_response) :: now
        real(real64) :: temperature, pressure, weight
        integer(int32) :: c, d, slot
        logical :: ok

        do c = 1, PLANE_DIMENSION*PLANE_DIMENSION
            g3(c) = 0.0d0
        end do

        slot = material_slot_of(a, a_off)
        call point_state(u, u_off, constants, temperature, pressure)
        ok = .false.
        call eval_closure_response(slot, temperature, pressure, now, ok)
        if (.not. ok) then
            call note_kernel_fault(slot, temperature, pressure)
            return
        end if

        weight = now%thermal_conductivity*constants(CONST_TEMPERATURE_SCALE) &
                 /constants(CONST_THERMAL_SCALE)
        do c = 1, PLANE_DIMENSION
            do d = 1, PLANE_DIMENSION
                if (c == d) g3((c - 1)*PLANE_DIMENSION + d) = weight
            end do
        end do
    end subroutine g3_temperature_temperature

    !> \(\partial f_0^{T}/\partial p_w=(\partial U/\partial p_w)p_{\mathrm{scale}}/(\Delta t\,R_T)\).
    subroutine g0_temperature_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_g0_Tp")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time, u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g0(*)

        type(type_closure_tangents) :: tangents
        integer(int32) :: slot
        logical :: ok

        g0(1) = 0.0d0
        if (constants(CONST_DT) <= 0.0d0) return
        call point_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        g0(1) = tangents%denthalpy_dpressure*constants(CONST_PRESSURE_SCALE) &
                *constants(CONST_BDF_CURRENT) &
                /(constants(CONST_DT)*constants(CONST_THERMAL_SCALE))
    end subroutine g0_temperature_pressure

    !> \(\partial f_1^{T}/\partial p_w=(\partial\lambda/\partial p_w)p_{\mathrm{scale}}T_{\mathrm{scale}}\nabla u_T/R_T\).
    subroutine g2_temperature_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_g2_Tp")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time, u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g2(*)

        type(type_closure_tangents) :: tangents
        real(real64) :: weight
        integer(int32) :: d, slot
        logical :: ok

        do d = 1, PLANE_DIMENSION
            g2(d) = 0.0d0
        end do
        call point_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        weight = tangents%dthermal_dpressure &
                 *constants(CONST_PRESSURE_SCALE)*constants(CONST_TEMPERATURE_SCALE) &
                 /constants(CONST_THERMAL_SCALE)
        do d = 1, PLANE_DIMENSION
            g2(d) = weight*u_x(u_off_x(field_number(constants, CONST_TEMPERATURE_FIELD) + 1) + d)
        end do
    end subroutine g2_temperature_pressure

    !> \(\partial f_0^{q}/\partial\mathbf q=\mathbf\Lambda^{-1}(T,p_w)q_{\mathrm{scale}}/p_{\mathrm{scale}}\).
    subroutine g0_flux_flux( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_g0_qq")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time, u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g0(*)

        type(type_closure_response) :: now
        real(real64) :: inverse_mobility(PLANE_DIMENSION*PLANE_DIMENSION)
        real(real64) :: temperature, pressure, weight
        integer(int32) :: k, slot
        logical :: ok

        do k = 1, PLANE_DIMENSION*PLANE_DIMENSION
            g0(k) = 0.0d0
        end do

        slot = material_slot_of(a, a_off)
        call point_state(u, u_off, constants, temperature, pressure)
        ok = .false.
        call eval_closure_response(slot, temperature, pressure, now, ok)
        if (.not. ok) then
            call note_kernel_fault(slot, temperature, pressure)
            return
        end if

        weight = constants(CONST_FLUX_SCALE)/constants(CONST_PRESSURE_SCALE)
        inverse_mobility = inverse_mobility_of(a, a_off)/now%relative_conductivity
        do k = 1, PLANE_DIMENSION*PLANE_DIMENSION
            g0(k) = inverse_mobility(k)*weight
        end do
    end subroutine g0_flux_flux

    !> \(\partial f_0^{q}/\partial T=-(\partial\ln k_r/\partial T)T_{\mathrm{scale}}
    !> \mathbf\Lambda^{-1}\mathbf q/p_{\mathrm{scale}}\).
    subroutine g0_flux_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_g0_qT")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time, u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g0(*)

        type(type_closure_tangents) :: tangents
        type(type_closure_response) :: now
        real(real64) :: inverse_mobility(PLANE_DIMENSION*PLANE_DIMENSION)
        real(real64) :: flux(PLANE_DIMENSION), resistance(PLANE_DIMENSION)
        real(real64) :: temperature, pressure, flux_scale, pressure_scale
        integer(int32) :: i, slot
        logical :: ok

        do i = 1, PLANE_DIMENSION
            g0(i) = 0.0d0
        end do

        flux_scale = constants(CONST_FLUX_SCALE)
        pressure_scale = constants(CONST_PRESSURE_SCALE)
        slot = material_slot_of(a, a_off)
        call point_state(u, u_off, constants, temperature, pressure)

        ok = .false.
        call eval_closure_response(slot, temperature, pressure, now, ok)
        if (ok) call eval_closure_tangents(slot, temperature, pressure, tangents, ok)
        if (.not. ok) then
            call note_kernel_fault(slot, temperature, pressure)
            return
        end if

        inverse_mobility = inverse_mobility_of(a, a_off)/now%relative_conductivity
        do i = 1, PLANE_DIMENSION
            flux(i) = flux_scale*u(u_off(field_number(constants, CONST_FLUX_FIELD) + 1) + i)
        end do
        resistance = resistance_of(inverse_mobility, flux)

        do i = 1, PLANE_DIMENSION
            g0(i) = -tangents%dlog_conductivity_dtemperature &
                    *constants(CONST_TEMPERATURE_SCALE)*resistance(i)/pressure_scale
        end do
    end subroutine g0_flux_temperature

    !> \(\partial f_1^{q}/\partial T=-(\partial p_l/\partial T)T_{\mathrm{scale}}\mathbf I/p_{\mathrm{scale}}\).
    !>
    !> This is the cryosuction block: the temperature gradient's whole effect on
    !> the water flux, and the only one.
    subroutine g2_flux_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_g2_qT")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time, u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g2(*)

        type(type_closure_tangents) :: tangents
        real(real64) :: weight
        integer(int32) :: c, d, slot
        logical :: ok

        do c = 1, PLANE_DIMENSION*PLANE_DIMENSION
            g2(c) = 0.0d0
        end do
        call point_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        weight = tangents%dliquid_dtemperature*constants(CONST_TEMPERATURE_SCALE) &
                 /constants(CONST_PRESSURE_SCALE)
        do c = 1, PLANE_DIMENSION
            do d = 1, PLANE_DIMENSION
                if (c == d) g2((c - 1)*PLANE_DIMENSION + d) = -weight
            end do
        end do
    end subroutine g2_flux_temperature

    !> \(\partial f_0^{q}/\partial p_w=-(\partial\ln k_r/\partial p_w)\mathbf\Lambda^{-1}\mathbf q\).
    subroutine g0_flux_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_g0_qp")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time, u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g0(*)

        type(type_closure_tangents) :: tangents
        type(type_closure_response) :: now
        real(real64) :: inverse_mobility(PLANE_DIMENSION*PLANE_DIMENSION)
        real(real64) :: flux(PLANE_DIMENSION), resistance(PLANE_DIMENSION)
        real(real64) :: temperature, pressure, flux_scale
        integer(int32) :: i, slot
        logical :: ok

        do i = 1, PLANE_DIMENSION
            g0(i) = 0.0d0
        end do

        flux_scale = constants(CONST_FLUX_SCALE)
        slot = material_slot_of(a, a_off)
        call point_state(u, u_off, constants, temperature, pressure)

        ok = .false.
        call eval_closure_response(slot, temperature, pressure, now, ok)
        if (ok) call eval_closure_tangents(slot, temperature, pressure, tangents, ok)
        if (.not. ok) then
            call note_kernel_fault(slot, temperature, pressure)
            return
        end if

        inverse_mobility = inverse_mobility_of(a, a_off)/now%relative_conductivity
        do i = 1, PLANE_DIMENSION
            flux(i) = flux_scale*u(u_off(field_number(constants, CONST_FLUX_FIELD) + 1) + i)
        end do
        resistance = resistance_of(inverse_mobility, flux)

        ! The two scale factors cancel: the row is divided by the pressure scale
        ! and the pressure trial function is measured in it.
        do i = 1, PLANE_DIMENSION
            g0(i) = -tangents%dlog_conductivity_dpressure*resistance(i)
        end do
    end subroutine g0_flux_pressure

    !> \(\partial f_1^{q}/\partial p_w=-(\partial p_l/\partial p_w)\mathbf I\).
    subroutine g2_flux_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_g2_qp")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time, u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g2(*)

        type(type_closure_tangents) :: tangents
        integer(int32) :: c, d, slot
        logical :: ok

        do c = 1, PLANE_DIMENSION*PLANE_DIMENSION
            g2(c) = 0.0d0
        end do
        call point_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        do c = 1, PLANE_DIMENSION
            do d = 1, PLANE_DIMENSION
                if (c == d) g2((c - 1)*PLANE_DIMENSION + d) = -tangents%dliquid_dpressure
            end do
        end do
    end subroutine g2_flux_pressure

    !> \(\partial f_0^{p}/\partial(\nabla\mathbf q)=L\,\mathbf I\).
    subroutine g1_pressure_flux( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g1 &
        ) bind(C, name="ftcms_g1_pq")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time, u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g1(*)

        real(real64) :: length_scale
        integer(int32) :: c, d

        length_scale = constants(CONST_LENGTH_SCALE)
        do c = 1, PLANE_DIMENSION
            do d = 1, PLANE_DIMENSION
                if (c == d) then
                    g1((c - 1)*PLANE_DIMENSION + d) = length_scale
                else
                    g1((c - 1)*PLANE_DIMENSION + d) = 0.0d0
                end if
            end do
        end do
    end subroutine g1_pressure_flux

    !> \(\partial f_0^{p}/\partial p_w=(\partial W/\partial p_w)p_{\mathrm{scale}}L/(\Delta t\,q_{\mathrm{scale}})\).
    subroutine g0_pressure_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_g0_pp")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time, u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g0(*)

        type(type_closure_tangents) :: tangents
        integer(int32) :: slot
        logical :: ok

        g0(1) = 0.0d0
        if (constants(CONST_DT) <= 0.0d0) return
        call point_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        g0(1) = tangents%dwater_dpressure*constants(CONST_PRESSURE_SCALE) &
                *constants(CONST_LENGTH_SCALE) &
                *constants(CONST_BDF_CURRENT) &
                /(constants(CONST_DT)*constants(CONST_FLUX_SCALE))
    end subroutine g0_pressure_pressure

    !> \(\partial f_0^{p}/\partial T=(\partial W/\partial T)T_{\mathrm{scale}}L/(\Delta t\,q_{\mathrm{scale}})\).
    !>
    !> Small but not zero: the total stored water of this closure is a function
    !> of \(p_w\) alone except for the vapour term, which carries the
    !> temperature. Registering the block is what lets a run measure that rather
    !> than assume it.
    subroutine g0_pressure_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_g0_pT")
        implicit none
        PetscInt, value :: dim, nf, nf_aux
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: u(*), u_t(*), u_x(*)
        PetscInt, intent(in) :: a_off(*), a_off_x(*)
        PetscScalar, intent(in) :: a(*), a_t(*), a_x(*)
        PetscReal, value :: time, u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g0(*)

        type(type_closure_tangents) :: tangents
        integer(int32) :: slot
        logical :: ok

        g0(1) = 0.0d0
        if (constants(CONST_DT) <= 0.0d0) return
        call point_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        g0(1) = tangents%dwater_dtemperature*constants(CONST_TEMPERATURE_SCALE) &
                *constants(CONST_LENGTH_SCALE) &
                *constants(CONST_BDF_CURRENT) &
                /(constants(CONST_DT)*constants(CONST_FLUX_SCALE))
    end subroutine g0_pressure_temperature

end module app_ftcms_kernels
