#include <petsc/finclude/petscsys.h>

!> The coupled \([T,\mathbf q_l,p_w]\) weak form, as PetscDS pointwise functions.
!>
!> One backward-Euler step of heat and water on one mesh, with \(T\) in
!> \(H^1\) (P1), \(\mathbf q_l\) in \(H(\mathrm{div})\) (RT0) and \(p_w\) in
!> \(L^2\) (P0):
!>
!> \[
!>   R_T(w)=\Bigl(\frac{U(T,p_w)-U^n}{\Delta t},w\Bigr)+(\lambda\nabla T,\nabla w),
!> \]
!> \[
!>   R_q(\mathbf v)=(\mathbf\Lambda^{-1}(T,p_w)\mathbf q_l-\rho_l\mathbf g,\mathbf v)
!>   -(p_l(T,p_w),\nabla\cdot\mathbf v),
!> \]
!> \[
!>   R_p(\omega)=\Bigl(\frac{W(T,p_w)-W^n}{\Delta t}+\nabla\cdot\mathbf q_l-s,\omega\Bigr).
!> \]
!>
!> ## What is coupled, and where
!>
!> Every constitutive value comes from `app_mixed_hydraulic_frozen`, which is
!> the repository's own governing models. This module states no retention
!> curve, no freezing curve, no impedance law, no enthalpy and no mixing rule.
!>
!> Cryosuction enters once, as \(p_l=-\psi_{\mathrm{eff}}(T,p_w)\) in the Darcy
!> row, and nothing adds a second temperature-driven water flux beside it. The
!> chain rule \(\nabla p_l=(\partial p_l/\partial p_w)\nabla p_w+(\partial
!> p_l/\partial T)\nabla T\) is not written out here: \(p_l\) is evaluated as
!> the function it is and PETSc differentiates the row through the registered
!> \(g_2\) blocks, so the driver cannot end up counted twice.
!>
!> ## The history
!>
!> \(U^n\) and \(W^n\) are not stored. The accepted temperature and pressure
!> are, in auxiliary fields with the same discretisations as the unknowns, and
!> the kernels re-evaluate the closure at them. The storage rate is then a
!> difference of one function of the state rather than of two projections of
!> it, and the quadrature history is exact at every quadrature point instead of
!> being averaged onto a cell.
!>
!> ## What is not here
!>
!> Advective heat transport by the liquid flux. The column below is conduction
!> dominated and the term is deliberately absent rather than approximated; the
!> energy balance reported by a run that uses these kernels is therefore the
!> conduction-plus-latent balance and is documented as such.
module app_mixed_hydraulic_column_pointwise
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_loc, c_f_pointer
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use :: petscsys
    use :: app_mixed_hydraulic_physics, only:MIXED_DIMENSION
    use :: app_mixed_hydraulic_pointwise, only:resistance_of
    use :: app_mixed_hydraulic_frozen, only:type_mixed_frozen_response, type_mixed_frozen_tangents, &
                                            eval_frozen_response, eval_frozen_tangents
    implicit none
    private

    public :: f0_column_temperature, f1_column_temperature
    public :: f0_column_flux, f1_column_flux
    public :: f0_column_conservation
    public :: g0_column_temperature_temperature, g2_column_temperature_temperature
    public :: g3_column_temperature_temperature
    public :: g0_column_temperature_pressure, g2_column_temperature_pressure
    public :: g0_column_flux_flux, g0_column_flux_temperature, g2_column_flux_temperature
    public :: g0_column_flux_pressure, g2_column_flux_pressure
    public :: g1_column_pressure_flux, g0_column_pressure_pressure, g0_column_pressure_temperature
    public :: bc_column_temperature, set_column_temperature_value, column_temperature_context
    public :: reset_column_temperatures
    public :: bd_f0_column_heat, bd_g0_column_heat
    public :: clear_column_fault, column_fault_recorded, column_fault_report
    public :: COLUMN_FIELD_TEMPERATURE, COLUMN_FIELD_FLUX, COLUMN_FIELD_PRESSURE
    public :: COLUMN_AUX_FIELDS
    public :: CA_INVERSE_MOBILITY, CA_MATERIAL, CA_SOURCE, CA_TEMPERATURE_OLD, CA_PRESSURE_OLD
    public :: CA_HEAT_COEFFICIENT, CA_HEAT_AMBIENT
    public :: CA_HEAT_COEFFICIENT_FROZEN, CA_HEAT_SWITCH_TEMPERATURE, CA_HEAT_SWITCH_WIDTH
    public :: CC_DT, CC_FLUX_SCALE, CC_PRESSURE_SCALE, CC_LENGTH_SCALE, CC_DENSITY
    public :: CC_GRAVITY_X, CC_GRAVITY_Z, CC_TEMPERATURE_SCALE, CC_THERMAL_SCALE, CC_PRESSURE_FIELD
    public :: COLUMN_CONSTANTS
    public :: column_state, column_state_old, column_material_slot, column_inverse_mobility
    public :: column_tangents, note_column_fault
    public :: MAX_COLUMN_TEMPERATURES

    !> Field numbers of the coupled layout, PETSc's zero-based numbering.
    !> Temperature comes first so that the hydraulic pair keeps its relative
    !> order; no kernel outside this module may assume the values.
    PetscInt, parameter :: COLUMN_FIELD_TEMPERATURE = 0
    PetscInt, parameter :: COLUMN_FIELD_FLUX = 1
    PetscInt, parameter :: COLUMN_FIELD_PRESSURE = 2

    ! Auxiliary fields, in the order PETSc offsets them.
    integer(int32), parameter :: CA_INVERSE_MOBILITY = 1 !< Saturated inverse mobility, row major, P0.
    integer(int32), parameter :: CA_MATERIAL = 2         !< Slot in the frozen closure table, P0.
    integer(int32), parameter :: CA_SOURCE = 3           !< Water source, 1/s, P0.
    integer(int32), parameter :: CA_TEMPERATURE_OLD = 4  !< Accepted temperature, C, P1.
    integer(int32), parameter :: CA_PRESSURE_OLD = 5     !< Accepted pore pressure, Pa, P0.
    ! The last two are read by the boundary kernels alone. A natural boundary's
    ! pointwise kernels take no context, so what distinguishes one
    ! heat-transfer boundary from another is the auxiliary vector PETSc looks up
    ! for that boundary's label value. The volume kernels never read them and
    ! the default auxiliary vector leaves them zero.
    integer(int32), parameter :: CA_HEAT_COEFFICIENT = 6 !< Heat transfer coefficient, W/(m2 K), P0.
    integer(int32), parameter :: CA_HEAT_AMBIENT = 7     !< Temperature of the medium, C, P0.
    ! The optional switch of the coefficient with the surface temperature: the
    ! coefficient below the switch, the switch temperature and the width the
    ! step is smoothed over. A boundary without a switch carries the same
    ! coefficient in both, which makes the step vanish.
    integer(int32), parameter :: CA_HEAT_COEFFICIENT_FROZEN = 8 !< Coefficient below the switch, W/(m2 K), P0.
    integer(int32), parameter :: CA_HEAT_SWITCH_TEMPERATURE = 9 !< Switch temperature, C, P0.
    integer(int32), parameter :: CA_HEAT_SWITCH_WIDTH = 10      !< Width of the smoothed step, K, P0.
    integer(int32), parameter :: COLUMN_AUX_FIELDS = 10

    ! Run constants, in the order PetscDSSetConstants receives them.
    integer(int32), parameter :: CC_DT = 1                !< Time step, s.
    integer(int32), parameter :: CC_FLUX_SCALE = 2        !< Flux unit of the unknown, m/s.
    integer(int32), parameter :: CC_PRESSURE_SCALE = 3    !< Pressure unit of the unknown, Pa.
    integer(int32), parameter :: CC_LENGTH_SCALE = 4      !< Length the conservation row is scaled by, m.
    integer(int32), parameter :: CC_DENSITY = 5           !< Reference liquid density, kg/m3.
    integer(int32), parameter :: CC_GRAVITY_X = 6         !< Gravity along intrinsic x, m/s2.
    integer(int32), parameter :: CC_GRAVITY_Z = 7         !< Gravity along intrinsic z, m/s2.
    integer(int32), parameter :: CC_TEMPERATURE_SCALE = 8 !< Temperature unit of the unknown, K.
    integer(int32), parameter :: CC_THERMAL_SCALE = 9     !< The heat row's divisor, W/m3.
    ! The pressure's field number differs between the mixed layout [T, q, p]
    ! and the primal one [T, p]; every kernel reads the pressure through it, so
    ! the temperature kernels serve both layouts unchanged.
    integer(int32), parameter :: CC_PRESSURE_FIELD = 10   !< Zero-based field number of the pressure.
    integer(int32), parameter :: COLUMN_CONSTANTS = 10

    !> How many prescribed-temperature boundaries one run may register.
    integer(int32), parameter :: MAX_COLUMN_TEMPERATURES = 8

    !> The prescribed temperature of one boundary, reachable from the essential
    !> callback as its context.
    type, bind(C) :: type_column_temperature_value
        real(c_double) :: temperature
        real(c_double) :: temperature_scale
    end type type_column_temperature_value

    type(type_column_temperature_value), save, target :: column_temperatures(MAX_COLUMN_TEMPERATURES)

    ! What the last unusable closure evaluation was. A pointwise kernel has
    ! nowhere to return an error to, so it records the fault and writes a
    ! nonfinite residual, which is what refuses the step.
    logical, save :: column_fault = .false.
    integer(int32), save :: column_fault_slot = 0
    real(real64), save :: column_fault_temperature = 0.0_real64
    real(real64), save :: column_fault_pressure = 0.0_real64

contains

    !> The cell's saturated inverse mobility, row major, from its auxiliary field.
    !>
    !> This layout's field number, not the isothermal one: the two auxiliary
    !> sections carry different fields in different places, and reading one
    !> through the other's offsets is how a kernel silently gets a zero.
    pure function column_inverse_mobility(a, a_off) result(inverse_mobility)
        implicit none
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        real(real64) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION)

        integer(int32) :: k

        do k = 1, MIXED_DIMENSION*MIXED_DIMENSION
            inverse_mobility(k) = real(a(a_off(CA_INVERSE_MOBILITY) + k), real64)
        end do
    end function column_inverse_mobility

    !> The cell's slot in the frozen closure table, from its auxiliary field.
    pure function column_material_slot(a, a_off) result(slot)
        implicit none
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        integer(int32) :: slot

        slot = nint(real(a(a_off(CA_MATERIAL) + 1), real64), int32)
    end function column_material_slot

    !> Forget every prescribed temperature of a previous run.
    subroutine reset_column_temperatures()
        implicit none

        integer(int32) :: i

        do i = 1, MAX_COLUMN_TEMPERATURES
            column_temperatures(i)%temperature = 0.0_c_double
            column_temperatures(i)%temperature_scale = 1.0_c_double
        end do
    end subroutine reset_column_temperatures

    !> Record the temperature one registered boundary prescribes.
    subroutine set_column_temperature_value(slot, temperature, temperature_scale)
        implicit none
        !> Index in 1..MAX_COLUMN_TEMPERATURES.
        integer(int32), intent(in) :: slot
        !> Prescribed temperature, C.
        real(real64), intent(in) :: temperature
        !> The temperature unit of the unknown, K, positive.
        real(real64), intent(in) :: temperature_scale

        column_temperatures(slot)%temperature = real(temperature, c_double)
        column_temperatures(slot)%temperature_scale = real(temperature_scale, c_double)
    end subroutine set_column_temperature_value

    !> The context PETSc hands back to the essential callback of one boundary.
    function column_temperature_context(slot) result(ctx)
        implicit none
        integer(int32), intent(in) :: slot
        type(c_ptr) :: ctx

        ctx = c_loc(column_temperatures(slot))
    end function column_temperature_context

    !> Newton cooling at a boundary, \(f_0^\partial=h(T)(T-T_\infty)/R_T\).
    !>
    !> The surface loses \(h(T-T_\infty)\) W/m2 to a medium at \(T_\infty\).
    !> The heat row is divided by \(R_T\) and its unknown is
    !> \(T/T_{\mathrm{scale}}\), so the boundary integrand is that flux over
    !> \(R_T\) -- the same outward-flux convention the conduction term carries,
    !> which is what makes a surface colder than the soil take heat out of it.
    !>
    !> \(h\) and \(T_\infty\) come from the auxiliary vector attached to this
    !> boundary's label value, not from the run constants: the constants are one
    !> array for the whole system and two boundaries need two pairs.
    subroutine bd_f0_column_heat( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, normal, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_column_bd_f0_heat")
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

        real(real64) :: coefficient, slope, ambient, temperature

        ambient = real(a(a_off(CA_HEAT_AMBIENT) + 1), real64)
        temperature = real(constants(CC_TEMPERATURE_SCALE), real64)*real(u(u_off(1) + 1), real64)
        call heat_coefficient_at(a, a_off, temperature, coefficient, slope)

        f0(1) = real(coefficient*(temperature - ambient)/real(constants(CC_THERMAL_SCALE), real64), &
                     PETSC_REAL_KIND)
    end subroutine bd_f0_column_heat

    !> The coefficient of one heat-transfer boundary at the surface temperature.
    !>
    !> \[
    !>   h(T)=h_f+(h_w-h_f)\,S(u),\qquad
    !>   u=\min\bigl(1,\max(0,(T-T_s)/w+\tfrac12)\bigr),\qquad S(u)=u^2(3-2u),
    !> \]
    !>
    !> a step from \(h_f\) below the switch temperature \(T_s\) to \(h_w\)
    !> above it, smoothed to C1 over the width \(w\) so that the residual the
    !> Newton iteration differentiates has a derivative. \(h_w\) is
    !> `CA_HEAT_COEFFICIENT`; with \(w\le 0\) or \(h_f=h_w\) the coefficient is
    !> the constant \(h_w\). Returns \(\mathrm dh/\mathrm dT\) with it.
    pure subroutine heat_coefficient_at(a, a_off, temperature, coefficient, slope)
        implicit none
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        !> Surface temperature, C.
        real(real64), intent(in) :: temperature
        !> \(h(T)\), W/(m2 K); overwritten.
        real(real64), intent(inout) :: coefficient
        !> \(\mathrm dh/\mathrm dT\), W/(m2 K2); overwritten.
        real(real64), intent(inout) :: slope

        real(real64) :: warm, frozen, switch, width, fraction

        warm = real(a(a_off(CA_HEAT_COEFFICIENT) + 1), real64)
        frozen = real(a(a_off(CA_HEAT_COEFFICIENT_FROZEN) + 1), real64)
        switch = real(a(a_off(CA_HEAT_SWITCH_TEMPERATURE) + 1), real64)
        width = real(a(a_off(CA_HEAT_SWITCH_WIDTH) + 1), real64)

        coefficient = warm
        slope = 0.0_real64
        if (.not. (width > 0.0_real64) .or. frozen == warm) return

        fraction = (temperature - switch)/width + 0.5_real64
        if (fraction <= 0.0_real64) then
            coefficient = frozen
        else if (fraction < 1.0_real64) then
            coefficient = frozen + (warm - frozen)*fraction*fraction*(3.0_real64 - 2.0_real64*fraction)
            slope = (warm - frozen)*6.0_real64*fraction*(1.0_real64 - fraction)/width
        end if
    end subroutine heat_coefficient_at

    !> \(\partial f_0^\partial/\partial u_T=(h+h'(T)(T-T_\infty))\,T_{\mathrm{scale}}/R_T\).
    subroutine bd_g0_column_heat( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, normal, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_column_bd_g0_heat")
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

        real(real64) :: coefficient, slope, ambient, temperature

        ambient = real(a(a_off(CA_HEAT_AMBIENT) + 1), real64)
        temperature = real(constants(CC_TEMPERATURE_SCALE), real64)*real(u(u_off(1) + 1), real64)
        call heat_coefficient_at(a, a_off, temperature, coefficient, slope)

        ! d/dT [h(T)(T - T_inf)] = h + h'(T)(T - T_inf)
        g0(1) = real((coefficient + slope*(temperature - ambient))*real(constants(CC_TEMPERATURE_SCALE), real64) &
                     /real(constants(CC_THERMAL_SCALE), real64), PETSC_REAL_KIND)
    end subroutine bd_g0_column_heat

    !> Forget any closure fault a previous trial recorded.
    subroutine clear_column_fault()
        implicit none

        column_fault = .false.
        column_fault_slot = 0
        column_fault_temperature = 0.0_real64
        column_fault_pressure = 0.0_real64
    end subroutine clear_column_fault

    !> Whether a kernel could not evaluate the closure since the last clear.
    logical function column_fault_recorded()
        implicit none

        column_fault_recorded = column_fault
    end function column_fault_recorded

    !> What the recorded fault was.
    subroutine column_fault_report(slot, temperature, pressure)
        implicit none
        !> The material slot of the offending cell, zero when none was recorded.
        integer(int32), intent(inout) :: slot
        !> The state it was evaluated at.
        real(real64), intent(inout) :: temperature, pressure

        slot = column_fault_slot
        temperature = column_fault_temperature
        pressure = column_fault_pressure
    end subroutine column_fault_report

    !> Record that a kernel asked the closure for a state it cannot answer.
    subroutine note_column_fault(slot, temperature, pressure)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64), intent(in) :: temperature, pressure

        if (column_fault) return
        column_fault = .true.
        column_fault_slot = slot
        column_fault_temperature = temperature
        column_fault_pressure = pressure
    end subroutine note_column_fault

    !> The physical state of this quadrature point, from the unknowns.
    !>
    !> The pressure is read from the field `CC_PRESSURE_FIELD` names, so the
    !> same accessor serves the mixed layout and the primal one.
    pure subroutine column_state(u, u_off, constants, temperature, pressure)
        implicit none
        PetscScalar, intent(in) :: u(*)
        PetscInt, intent(in) :: u_off(*)
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: temperature, pressure

        integer(int32) :: field

        field = nint(real(constants(CC_PRESSURE_FIELD), real64), int32)
        temperature = real(constants(CC_TEMPERATURE_SCALE), real64)*real(u(u_off(1) + 1), real64)
        pressure = real(constants(CC_PRESSURE_SCALE), real64)*real(u(u_off(field + 1) + 1), real64)
    end subroutine column_state

    !> The accepted state of this quadrature point, from the auxiliary fields.
    pure subroutine column_state_old(a, a_off, temperature, pressure)
        implicit none
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        real(real64), intent(inout) :: temperature, pressure

        temperature = real(a(a_off(CA_TEMPERATURE_OLD) + 1), real64)
        pressure = real(a(a_off(CA_PRESSURE_OLD) + 1), real64)
    end subroutine column_state_old

    !> Heat storage, \(f_0=(U(T,p_w)-U^n)/(\Delta t\,R_T)\).
    subroutine f0_column_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_column_f0_temperature")
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

        type(type_mixed_frozen_response) :: now, accepted
        real(real64) :: temperature, pressure, temperature_old, pressure_old
        real(real64) :: dt, thermal_scale
        integer(int32) :: slot
        logical :: ok

        dt = real(constants(CC_DT), real64)
        thermal_scale = real(constants(CC_THERMAL_SCALE), real64)
        slot = column_material_slot(a, a_off)
        call column_state(u, u_off, constants, temperature, pressure)
        call column_state_old(a, a_off, temperature_old, pressure_old)

        ok = .false.
        call eval_frozen_response(slot, temperature, pressure, now, ok)
        if (ok) call eval_frozen_response(slot, temperature_old, pressure_old, accepted, ok)
        if (.not. ok) then
            call note_column_fault(slot, temperature, pressure)
            f0(1) = real(ieee_value(1.0_real64, ieee_quiet_nan), PETSC_REAL_KIND)
            return
        end if

        f0(1) = real((now%enthalpy - accepted%enthalpy)/(dt*thermal_scale), PETSC_REAL_KIND)
    end subroutine f0_column_temperature

    !> Conduction, \(f_1=\lambda T_{\mathrm{scale}}\nabla u_T/R_T\).
    subroutine f1_column_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f1 &
        ) bind(C, name="ftcms_column_f1_temperature")
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

        type(type_mixed_frozen_response) :: now
        real(real64) :: temperature, pressure, weight
        integer(int32) :: d, slot
        logical :: ok

        slot = column_material_slot(a, a_off)
        call column_state(u, u_off, constants, temperature, pressure)

        ok = .false.
        call eval_frozen_response(slot, temperature, pressure, now, ok)
        if (.not. ok) then
            call note_column_fault(slot, temperature, pressure)
            do d = 1, MIXED_DIMENSION
                f1(d) = real(ieee_value(1.0_real64, ieee_quiet_nan), PETSC_REAL_KIND)
            end do
            return
        end if

        weight = now%thermal_conductivity*real(constants(CC_TEMPERATURE_SCALE), real64) &
                 /real(constants(CC_THERMAL_SCALE), real64)
        do d = 1, MIXED_DIMENSION
            f1(d) = real(weight*real(u_x(u_off_x(1) + d), real64), PETSC_REAL_KIND)
        end do
    end subroutine f1_column_temperature

    !> Darcy's law, \(f_0=(\mathbf\Lambda^{-1}(T,p_w)\mathbf q_l-\rho_l\mathbf g)/p_{\mathrm{scale}}\).
    subroutine f0_column_flux( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_column_f0_flux")
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

        type(type_mixed_frozen_response) :: now
        real(real64) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION)
        real(real64) :: flux(MIXED_DIMENSION), gravity(MIXED_DIMENSION), resistance(MIXED_DIMENSION)
        real(real64) :: temperature, pressure, flux_scale, pressure_scale, density
        integer(int32) :: i, slot
        logical :: ok

        flux_scale = real(constants(CC_FLUX_SCALE), real64)
        pressure_scale = real(constants(CC_PRESSURE_SCALE), real64)
        density = real(constants(CC_DENSITY), real64)
        gravity(1) = real(constants(CC_GRAVITY_X), real64)
        gravity(2) = real(constants(CC_GRAVITY_Z), real64)

        slot = column_material_slot(a, a_off)
        call column_state(u, u_off, constants, temperature, pressure)

        ok = .false.
        call eval_frozen_response(slot, temperature, pressure, now, ok)
        if (.not. ok) then
            call note_column_fault(slot, temperature, pressure)
            do i = 1, MIXED_DIMENSION
                f0(i) = real(ieee_value(1.0_real64, ieee_quiet_nan), PETSC_REAL_KIND)
            end do
            return
        end if

        inverse_mobility = column_inverse_mobility(a, a_off)/now%relative_conductivity
        do i = 1, MIXED_DIMENSION
            flux(i) = flux_scale*real(u(u_off(2) + i), real64)
        end do
        resistance = resistance_of(inverse_mobility, flux)

        do i = 1, MIXED_DIMENSION
            f0(i) = real((resistance(i) - density*gravity(i))/pressure_scale, PETSC_REAL_KIND)
        end do
    end subroutine f0_column_flux

    !> The driver of Darcy's law, \(f_1=-p_l(T,p_w)\mathbf I/p_{\mathrm{scale}}\).
    !>
    !> The liquid pressure, not the solved \(p_w\): that difference is the
    !> cryosuction, and evaluating it here is the only place it enters.
    subroutine f1_column_flux( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f1 &
        ) bind(C, name="ftcms_column_f1_flux")
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

        type(type_mixed_frozen_response) :: now
        real(real64) :: temperature, pressure, driver
        integer(int32) :: c, d, slot
        logical :: ok

        slot = column_material_slot(a, a_off)
        call column_state(u, u_off, constants, temperature, pressure)

        ok = .false.
        call eval_frozen_response(slot, temperature, pressure, now, ok)
        if (.not. ok) then
            call note_column_fault(slot, temperature, pressure)
            do c = 1, MIXED_DIMENSION*MIXED_DIMENSION
                f1(c) = real(ieee_value(1.0_real64, ieee_quiet_nan), PETSC_REAL_KIND)
            end do
            return
        end if

        driver = now%liquid_pressure/real(constants(CC_PRESSURE_SCALE), real64)
        do c = 1, MIXED_DIMENSION
            do d = 1, MIXED_DIMENSION
                if (c == d) then
                    f1((c - 1)*MIXED_DIMENSION + d) = real(-driver, PETSC_REAL_KIND)
                else
                    f1((c - 1)*MIXED_DIMENSION + d) = 0.0_PETSC_REAL_KIND
                end if
            end do
        end do
    end subroutine f1_column_flux

    !> Water conservation, \(f_0=((W-W^n)/\Delta t+\nabla\cdot\mathbf q_l-s)L/q_{\mathrm{scale}}\).
    subroutine f0_column_conservation( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_column_f0_conservation")
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

        type(type_mixed_frozen_response) :: now, accepted
        real(real64) :: temperature, pressure, temperature_old, pressure_old
        real(real64) :: flux_scale, length_scale, dt, divergence, source, storage_rate
        integer(int32) :: c, slot
        logical :: ok

        flux_scale = real(constants(CC_FLUX_SCALE), real64)
        length_scale = real(constants(CC_LENGTH_SCALE), real64)
        dt = real(constants(CC_DT), real64)

        ! PETSc stores the gradient component-major: du_c/dx_d at c*dim + d.
        divergence = 0.0_real64
        do c = 1, MIXED_DIMENSION
            divergence = divergence + real(u_x(u_off_x(2) + (c - 1)*MIXED_DIMENSION + c), real64)
        end do

        source = real(a(a_off(CA_SOURCE) + 1), real64)
        slot = column_material_slot(a, a_off)
        call column_state(u, u_off, constants, temperature, pressure)
        call column_state_old(a, a_off, temperature_old, pressure_old)

        storage_rate = 0.0_real64
        if (dt > 0.0_real64) then
            ok = .false.
            call eval_frozen_response(slot, temperature, pressure, now, ok)
            if (ok) call eval_frozen_response(slot, temperature_old, pressure_old, accepted, ok)
            if (.not. ok) then
                call note_column_fault(slot, temperature, pressure)
                f0(1) = real(ieee_value(1.0_real64, ieee_quiet_nan), PETSC_REAL_KIND)
                return
            end if
            storage_rate = (now%stored_water - accepted%stored_water)/dt
        end if

        f0(1) = real((storage_rate + flux_scale*divergence - source)*length_scale/flux_scale, PETSC_REAL_KIND)
    end subroutine f0_column_conservation

    !> The tangents of this quadrature point, or a recorded fault.
    subroutine column_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        implicit none
        PetscScalar, intent(in) :: u(*)
        PetscInt, intent(in) :: u_off(*)
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        PetscScalar, intent(in) :: constants(*)
        integer(int32), intent(inout) :: slot
        type(type_mixed_frozen_tangents), intent(inout) :: tangents
        logical, intent(inout) :: ok

        real(real64) :: temperature, pressure

        slot = column_material_slot(a, a_off)
        call column_state(u, u_off, constants, temperature, pressure)
        ok = .false.
        call eval_frozen_tangents(slot, temperature, pressure, tangents, ok)
        if (.not. ok) call note_column_fault(slot, temperature, pressure)
    end subroutine column_tangents

    !> \(\partial f_0^{T}/\partial T=(\partial U/\partial T)T_{\mathrm{scale}}/(\Delta t\,R_T)\).
    subroutine g0_column_temperature_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_column_g0_TT")
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

        type(type_mixed_frozen_tangents) :: tangents
        integer(int32) :: slot
        logical :: ok

        g0(1) = 0.0_PETSC_REAL_KIND
        if (real(constants(CC_DT), real64) <= 0.0_real64) return
        call column_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        g0(1) = real(tangents%denthalpy_dtemperature*real(constants(CC_TEMPERATURE_SCALE), real64) &
                     /(real(constants(CC_DT), real64)*real(constants(CC_THERMAL_SCALE), real64)), PETSC_REAL_KIND)
    end subroutine g0_column_temperature_temperature

    !> \(\partial f_1^{T}/\partial T=(\partial\lambda/\partial T)T_{\mathrm{scale}}^2\nabla u_T/R_T\).
    subroutine g2_column_temperature_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_column_g2_TT")
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

        type(type_mixed_frozen_tangents) :: tangents
        real(real64) :: weight
        integer(int32) :: d, slot
        logical :: ok

        do d = 1, MIXED_DIMENSION
            g2(d) = 0.0_PETSC_REAL_KIND
        end do
        call column_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        weight = tangents%dthermal_dtemperature &
                 *real(constants(CC_TEMPERATURE_SCALE), real64)**2 &
                 /real(constants(CC_THERMAL_SCALE), real64)
        do d = 1, MIXED_DIMENSION
            g2(d) = real(weight*real(u_x(u_off_x(1) + d), real64), PETSC_REAL_KIND)
        end do
    end subroutine g2_column_temperature_temperature

    !> \(\partial f_1^{T}/\partial\nabla T=\lambda T_{\mathrm{scale}}\mathbf I/R_T\).
    subroutine g3_column_temperature_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g3 &
        ) bind(C, name="ftcms_column_g3_TT")
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

        type(type_mixed_frozen_response) :: now
        real(real64) :: temperature, pressure, weight
        integer(int32) :: c, d, slot
        logical :: ok

        do c = 1, MIXED_DIMENSION*MIXED_DIMENSION
            g3(c) = 0.0_PETSC_REAL_KIND
        end do

        slot = column_material_slot(a, a_off)
        call column_state(u, u_off, constants, temperature, pressure)
        ok = .false.
        call eval_frozen_response(slot, temperature, pressure, now, ok)
        if (.not. ok) then
            call note_column_fault(slot, temperature, pressure)
            return
        end if

        weight = now%thermal_conductivity*real(constants(CC_TEMPERATURE_SCALE), real64) &
                 /real(constants(CC_THERMAL_SCALE), real64)
        do c = 1, MIXED_DIMENSION
            do d = 1, MIXED_DIMENSION
                if (c == d) g3((c - 1)*MIXED_DIMENSION + d) = real(weight, PETSC_REAL_KIND)
            end do
        end do
    end subroutine g3_column_temperature_temperature

    !> \(\partial f_0^{T}/\partial p_w=(\partial U/\partial p_w)p_{\mathrm{scale}}/(\Delta t\,R_T)\).
    subroutine g0_column_temperature_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_column_g0_Tp")
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

        type(type_mixed_frozen_tangents) :: tangents
        integer(int32) :: slot
        logical :: ok

        g0(1) = 0.0_PETSC_REAL_KIND
        if (real(constants(CC_DT), real64) <= 0.0_real64) return
        call column_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        g0(1) = real(tangents%denthalpy_dpressure*real(constants(CC_PRESSURE_SCALE), real64) &
                     /(real(constants(CC_DT), real64)*real(constants(CC_THERMAL_SCALE), real64)), PETSC_REAL_KIND)
    end subroutine g0_column_temperature_pressure

    !> \(\partial f_1^{T}/\partial p_w=(\partial\lambda/\partial p_w)p_{\mathrm{scale}}T_{\mathrm{scale}}\nabla u_T/R_T\).
    subroutine g2_column_temperature_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_column_g2_Tp")
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

        type(type_mixed_frozen_tangents) :: tangents
        real(real64) :: weight
        integer(int32) :: d, slot
        logical :: ok

        do d = 1, MIXED_DIMENSION
            g2(d) = 0.0_PETSC_REAL_KIND
        end do
        call column_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        weight = tangents%dthermal_dpressure &
                 *real(constants(CC_PRESSURE_SCALE), real64)*real(constants(CC_TEMPERATURE_SCALE), real64) &
                 /real(constants(CC_THERMAL_SCALE), real64)
        do d = 1, MIXED_DIMENSION
            g2(d) = real(weight*real(u_x(u_off_x(1) + d), real64), PETSC_REAL_KIND)
        end do
    end subroutine g2_column_temperature_pressure

    !> \(\partial f_0^{q}/\partial\mathbf q=\mathbf\Lambda^{-1}(T,p_w)q_{\mathrm{scale}}/p_{\mathrm{scale}}\).
    subroutine g0_column_flux_flux( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_column_g0_qq")
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

        type(type_mixed_frozen_response) :: now
        real(real64) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION)
        real(real64) :: temperature, pressure, weight
        integer(int32) :: k, slot
        logical :: ok

        do k = 1, MIXED_DIMENSION*MIXED_DIMENSION
            g0(k) = 0.0_PETSC_REAL_KIND
        end do

        slot = column_material_slot(a, a_off)
        call column_state(u, u_off, constants, temperature, pressure)
        ok = .false.
        call eval_frozen_response(slot, temperature, pressure, now, ok)
        if (.not. ok) then
            call note_column_fault(slot, temperature, pressure)
            return
        end if

        weight = real(constants(CC_FLUX_SCALE), real64)/real(constants(CC_PRESSURE_SCALE), real64)
        inverse_mobility = column_inverse_mobility(a, a_off)/now%relative_conductivity
        do k = 1, MIXED_DIMENSION*MIXED_DIMENSION
            g0(k) = real(inverse_mobility(k)*weight, PETSC_REAL_KIND)
        end do
    end subroutine g0_column_flux_flux

    !> \(\partial f_0^{q}/\partial T=-(\partial\ln k_r/\partial T)T_{\mathrm{scale}}\mathbf\Lambda^{-1}\mathbf q/p_{\mathrm{scale}}\).
    subroutine g0_column_flux_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_column_g0_qT")
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

        type(type_mixed_frozen_tangents) :: tangents
        type(type_mixed_frozen_response) :: now
        real(real64) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION)
        real(real64) :: flux(MIXED_DIMENSION), resistance(MIXED_DIMENSION)
        real(real64) :: temperature, pressure, flux_scale, pressure_scale
        integer(int32) :: i, slot
        logical :: ok

        do i = 1, MIXED_DIMENSION
            g0(i) = 0.0_PETSC_REAL_KIND
        end do

        flux_scale = real(constants(CC_FLUX_SCALE), real64)
        pressure_scale = real(constants(CC_PRESSURE_SCALE), real64)
        slot = column_material_slot(a, a_off)
        call column_state(u, u_off, constants, temperature, pressure)

        ok = .false.
        call eval_frozen_response(slot, temperature, pressure, now, ok)
        if (ok) call eval_frozen_tangents(slot, temperature, pressure, tangents, ok)
        if (.not. ok) then
            call note_column_fault(slot, temperature, pressure)
            return
        end if

        inverse_mobility = column_inverse_mobility(a, a_off)/now%relative_conductivity
        do i = 1, MIXED_DIMENSION
            flux(i) = flux_scale*real(u(u_off(2) + i), real64)
        end do
        resistance = resistance_of(inverse_mobility, flux)

        do i = 1, MIXED_DIMENSION
            g0(i) = real(-tangents%dlog_conductivity_dtemperature &
                         *real(constants(CC_TEMPERATURE_SCALE), real64)*resistance(i)/pressure_scale, PETSC_REAL_KIND)
        end do
    end subroutine g0_column_flux_temperature

    !> \(\partial f_1^{q}/\partial T=-(\partial p_l/\partial T)T_{\mathrm{scale}}\mathbf I/p_{\mathrm{scale}}\).
    !>
    !> This is the cryosuction block: the temperature gradient's whole effect on
    !> the water flux, and the only one.
    subroutine g2_column_flux_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_column_g2_qT")
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

        type(type_mixed_frozen_tangents) :: tangents
        real(real64) :: weight
        integer(int32) :: c, d, slot
        logical :: ok

        do c = 1, MIXED_DIMENSION*MIXED_DIMENSION
            g2(c) = 0.0_PETSC_REAL_KIND
        end do
        call column_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        weight = tangents%dliquid_dtemperature*real(constants(CC_TEMPERATURE_SCALE), real64) &
                 /real(constants(CC_PRESSURE_SCALE), real64)
        do c = 1, MIXED_DIMENSION
            do d = 1, MIXED_DIMENSION
                if (c == d) g2((c - 1)*MIXED_DIMENSION + d) = real(-weight, PETSC_REAL_KIND)
            end do
        end do
    end subroutine g2_column_flux_temperature

    !> \(\partial f_0^{q}/\partial p_w=-(\partial\ln k_r/\partial p_w)\mathbf\Lambda^{-1}\mathbf q\).
    subroutine g0_column_flux_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_column_g0_qp")
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

        type(type_mixed_frozen_tangents) :: tangents
        type(type_mixed_frozen_response) :: now
        real(real64) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION)
        real(real64) :: flux(MIXED_DIMENSION), resistance(MIXED_DIMENSION)
        real(real64) :: temperature, pressure, flux_scale
        integer(int32) :: i, slot
        logical :: ok

        do i = 1, MIXED_DIMENSION
            g0(i) = 0.0_PETSC_REAL_KIND
        end do

        flux_scale = real(constants(CC_FLUX_SCALE), real64)
        slot = column_material_slot(a, a_off)
        call column_state(u, u_off, constants, temperature, pressure)

        ok = .false.
        call eval_frozen_response(slot, temperature, pressure, now, ok)
        if (ok) call eval_frozen_tangents(slot, temperature, pressure, tangents, ok)
        if (.not. ok) then
            call note_column_fault(slot, temperature, pressure)
            return
        end if

        inverse_mobility = column_inverse_mobility(a, a_off)/now%relative_conductivity
        do i = 1, MIXED_DIMENSION
            flux(i) = flux_scale*real(u(u_off(2) + i), real64)
        end do
        resistance = resistance_of(inverse_mobility, flux)

        ! The two scale factors cancel: the row is divided by the pressure scale
        ! and the pressure trial function is measured in it.
        do i = 1, MIXED_DIMENSION
            g0(i) = real(-tangents%dlog_conductivity_dpressure*resistance(i), PETSC_REAL_KIND)
        end do
    end subroutine g0_column_flux_pressure

    !> \(\partial f_1^{q}/\partial p_w=-(\partial p_l/\partial p_w)\mathbf I\).
    subroutine g2_column_flux_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_column_g2_qp")
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

        type(type_mixed_frozen_tangents) :: tangents
        integer(int32) :: c, d, slot
        logical :: ok

        do c = 1, MIXED_DIMENSION*MIXED_DIMENSION
            g2(c) = 0.0_PETSC_REAL_KIND
        end do
        call column_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        do c = 1, MIXED_DIMENSION
            do d = 1, MIXED_DIMENSION
                if (c == d) g2((c - 1)*MIXED_DIMENSION + d) = real(-tangents%dliquid_dpressure, PETSC_REAL_KIND)
            end do
        end do
    end subroutine g2_column_flux_pressure

    !> \(\partial f_0^{p}/\partial(\nabla\mathbf q)=L\,\mathbf I\).
    subroutine g1_column_pressure_flux( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g1 &
        ) bind(C, name="ftcms_column_g1_pq")
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

        length_scale = real(constants(CC_LENGTH_SCALE), real64)
        do c = 1, MIXED_DIMENSION
            do d = 1, MIXED_DIMENSION
                if (c == d) then
                    g1((c - 1)*MIXED_DIMENSION + d) = real(length_scale, PETSC_REAL_KIND)
                else
                    g1((c - 1)*MIXED_DIMENSION + d) = 0.0_PETSC_REAL_KIND
                end if
            end do
        end do
    end subroutine g1_column_pressure_flux

    !> \(\partial f_0^{p}/\partial p_w=(\partial W/\partial p_w)p_{\mathrm{scale}}L/(\Delta t\,q_{\mathrm{scale}})\).
    subroutine g0_column_pressure_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_column_g0_pp")
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

        type(type_mixed_frozen_tangents) :: tangents
        integer(int32) :: slot
        logical :: ok

        g0(1) = 0.0_PETSC_REAL_KIND
        if (real(constants(CC_DT), real64) <= 0.0_real64) return
        call column_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        g0(1) = real(tangents%dwater_dpressure*real(constants(CC_PRESSURE_SCALE), real64) &
                     *real(constants(CC_LENGTH_SCALE), real64) &
                     /(real(constants(CC_DT), real64)*real(constants(CC_FLUX_SCALE), real64)), PETSC_REAL_KIND)
    end subroutine g0_column_pressure_pressure

    !> \(\partial f_0^{p}/\partial T=(\partial W/\partial T)T_{\mathrm{scale}}L/(\Delta t\,q_{\mathrm{scale}})\).
    !>
    !> Small but not zero: the total stored water of this closure is a function
    !> of \(p_w\) alone except for the vapour term, which carries the
    !> temperature. Registering the block is what lets a run measure that rather
    !> than assume it.
    subroutine g0_column_pressure_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_column_g0_pT")
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

        type(type_mixed_frozen_tangents) :: tangents
        integer(int32) :: slot
        logical :: ok

        g0(1) = 0.0_PETSC_REAL_KIND
        if (real(constants(CC_DT), real64) <= 0.0_real64) return
        call column_tangents(u, u_off, a, a_off, constants, slot, tangents, ok)
        if (.not. ok) return

        g0(1) = real(tangents%dwater_dtemperature*real(constants(CC_TEMPERATURE_SCALE), real64) &
                     *real(constants(CC_LENGTH_SCALE), real64) &
                     /(real(constants(CC_DT), real64)*real(constants(CC_FLUX_SCALE), real64)), PETSC_REAL_KIND)
    end subroutine g0_column_pressure_temperature

    !> The temperature a prescribed-temperature boundary constrains.
    function bc_column_temperature(dim, time, x, num_components, u, ctx) &
        bind(C, name="ftcms_column_bc_temperature") result(ierr)
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        type(type_column_temperature_value), pointer :: boundary

        call c_f_pointer(ctx, boundary)
        u(1) = real(boundary%temperature/boundary%temperature_scale, PETSC_REAL_KIND)
        ierr = 0
    end function bc_column_temperature

end module app_mixed_hydraulic_column_pointwise
