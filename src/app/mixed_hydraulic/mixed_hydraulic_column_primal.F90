#include <petsc/finclude/petscsys.h>

!> The primal (continuous Galerkin) water row of the coupled column.
!>
!> The same backward-Euler step as `app_mixed_hydraulic_column_pointwise`, with
!> the pore pressure \(p_w\) in \(H^1\) (P1 on a triangle, Q1 on a rectangle)
!> instead of the RT0 flux and P0 pressure pair:
!>
!> \[
!>   R_p(\omega)=\Bigl(\frac{W(T,p_w)-W^n}{\Delta t}-s,\omega\Bigr)
!>   -(\mathbf q_l,\nabla\omega),\qquad
!>   \mathbf q_l=-k_r(T,p_w)\,\mathbf M\bigl(\nabla p_l(T,p_w)-\rho_l\mathbf g\bigr),
!> \]
!>
!> with \(\mathbf M=\mathbf\Lambda_{\mathrm{sat}}\) the saturated mobility of the
!> cell and \(\nabla p_l=(\partial p_l/\partial p_w)\nabla p_w+(\partial
!> p_l/\partial T)\nabla T\). The closed boundary is the natural condition of
!> this form and needs no constraint.
!>
!> ## What is shared with the mixed form
!>
!> Everything but the hydraulic space. The closure, the temperature row and its
!> kernels, the storage rate, the quadrature, the scaling of the water row by
!> \(L/q_{\mathrm{scale}}\) and the Newton solver are the mixed column's own, so
!> a difference between the two runs is a difference of discretisation.
!>
!> ## What the primal form needs that the mixed one does not
!>
!> The mixed Darcy row tests the liquid pressure against \(\nabla\cdot\mathbf
!> v\), so \(p_l\) enters as a value and only its first derivatives reach the
!> Jacobian. Here \(\nabla p_l\) is in the residual, so the Jacobian needs the
!> derivatives of \(\partial p_l/\partial p_w\) and \(\partial p_l/\partial T\).
!> The closure publishes the first derivatives only; the second ones are
!> central differences of them.
module app_mixed_hydraulic_column_primal
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use :: petscsys
    use :: app_mixed_hydraulic_physics, only:MIXED_DIMENSION
    use :: app_mixed_hydraulic_frozen, only:type_mixed_frozen_response, type_mixed_frozen_tangents, &
                                            eval_frozen_response, eval_frozen_tangents
    use :: app_mixed_hydraulic_column_pointwise, only:column_state, column_state_old, column_material_slot, &
                                                      column_inverse_mobility, note_column_fault, &
                                                      CA_SOURCE, CC_DT, CC_FLUX_SCALE, CC_PRESSURE_SCALE, &
                                                      CC_LENGTH_SCALE, CC_DENSITY, CC_GRAVITY_X, CC_GRAVITY_Z, &
                                                      CC_TEMPERATURE_SCALE, CC_PRESSURE_FIELD
    implicit none
    private

    public :: COLUMN_PRIMAL_FIELD_PRESSURE
    public :: f0_primal_conservation, f1_primal_conservation
    public :: g2_primal_pressure_pressure, g3_primal_pressure_pressure
    public :: g2_primal_pressure_temperature, g3_primal_pressure_temperature
    public :: compute_primal_darcy_flux

    !> Field number of the pressure in the primal layout [T, p], PETSc's
    !> zero-based numbering.
    PetscInt, parameter :: COLUMN_PRIMAL_FIELD_PRESSURE = 1

    ! Offsets of the central differences of the liquid-pressure tangents: the
    ! steps the closure differences its own tangents with.
    real(real64), parameter :: HESSIAN_PRESSURE_STEP = 1.0_real64
    real(real64), parameter :: HESSIAN_TEMPERATURE_STEP = 1.0e-4_real64

    !> Everything one quadrature point of the water row needs.
    type :: type_primal_point
        integer(int32) :: slot = 0
        !> Temperature, C, and pore pressure, Pa.
        real(real64) :: temperature = 0.0_real64
        real(real64) :: pressure = 0.0_real64
        !> Physical gradients of the two unknowns, K/m and Pa/m.
        real(real64) :: temperature_gradient(MIXED_DIMENSION) = 0.0_real64
        real(real64) :: pressure_gradient(MIXED_DIMENSION) = 0.0_real64
        !> Saturated mobility of the cell, m2/(Pa s), row major.
        real(real64) :: mobility(MIXED_DIMENSION, MIXED_DIMENSION) = 0.0_real64
        !> \(\nabla p_l-\rho_l\mathbf g\), Pa/m.
        real(real64) :: driver(MIXED_DIMENSION) = 0.0_real64
        type(type_mixed_frozen_response) :: response
        type(type_mixed_frozen_tangents) :: tangents
    end type type_primal_point

contains

    !> The liquid flux of the primal form at one state, m/s.
    !>
    !> \[
    !>   \mathbf q_l=-k_r(T,p_w)\,\mathbf M\Bigl(\frac{\partial p_l}{\partial p_w}\nabla p_w
    !>   +\frac{\partial p_l}{\partial T}\nabla T-\rho_l\mathbf g\Bigr).
    !> \]
    !>
    !> Assumes the frozen closure is installed and `inverse_mobility` is a cell's
    !> saturated inverse mobility, row major and positive definite. No
    !> theoretical error bound available. \(O(1)\) work and memory. Fails by
    !> returning `ok = .false.` with `flux` untouched when the closure cannot
    !> answer at the state.
    subroutine compute_primal_darcy_flux(slot, temperature, pressure, temperature_gradient, pressure_gradient, &
                                         inverse_mobility, density, gravity, flux, ok)
        implicit none
        !> One-based slot in the frozen closure table.
        integer(int32), intent(in) :: slot
        !> Temperature, C, and pore pressure \(p_w\), Pa.
        real(real64), intent(in) :: temperature, pressure
        !> \(\nabla T\), K/m, and \(\nabla p_w\), Pa/m.
        real(real64), intent(in) :: temperature_gradient(MIXED_DIMENSION)
        real(real64), intent(in) :: pressure_gradient(MIXED_DIMENSION)
        !> Saturated inverse mobility of the cell, Pa s/m2, row major.
        real(real64), intent(in) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION)
        !> Reference liquid density, kg/m3, and gravity, m/s2.
        real(real64), intent(in) :: density
        real(real64), intent(in) :: gravity(MIXED_DIMENSION)
        !> \(\mathbf q_l\), m/s; overwritten when `ok`.
        real(real64), intent(inout) :: flux(MIXED_DIMENSION)
        !> True when the closure answered.
        logical, intent(inout) :: ok

        type(type_mixed_frozen_response) :: response
        type(type_mixed_frozen_tangents) :: tangents
        real(real64) :: mobility(MIXED_DIMENSION, MIXED_DIMENSION), driver(MIXED_DIMENSION)

        ok = .false.
        call eval_frozen_response(slot, temperature, pressure, response, ok)
        if (ok) call eval_frozen_tangents(slot, temperature, pressure, tangents, ok)
        if (.not. ok) return

        mobility = invert_inverse_mobility(inverse_mobility)
        driver = tangents%dliquid_dpressure*pressure_gradient &
                 + tangents%dliquid_dtemperature*temperature_gradient - density*gravity
        flux = -response%relative_conductivity*matmul(mobility, driver)
    end subroutine compute_primal_darcy_flux

    !> \(\mathbf M=(\mathbf\Lambda^{-1})^{-1}\) of a row-major 2x2 tensor.
    pure function invert_inverse_mobility(inverse_mobility) result(mobility)
        implicit none
        real(real64), intent(in) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION)
        real(real64) :: mobility(MIXED_DIMENSION, MIXED_DIMENSION)

        real(real64) :: determinant

        determinant = inverse_mobility(1)*inverse_mobility(4) - inverse_mobility(2)*inverse_mobility(3)
        mobility(1, 1) = inverse_mobility(4)/determinant
        mobility(1, 2) = -inverse_mobility(2)/determinant
        mobility(2, 1) = -inverse_mobility(3)/determinant
        mobility(2, 2) = inverse_mobility(1)/determinant
    end function invert_inverse_mobility

    !> Gather the state, gradients, mobility and closure of one quadrature point.
    subroutine read_primal_point(u, u_off, u_x, u_off_x, a, a_off, constants, point, ok)
        implicit none
        PetscScalar, intent(in) :: u(*), u_x(*)
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        PetscScalar, intent(in) :: constants(*)
        type(type_primal_point), intent(inout) :: point
        logical, intent(inout) :: ok

        real(real64) :: gravity(MIXED_DIMENSION)
        integer(int32) :: d, field

        field = nint(real(constants(CC_PRESSURE_FIELD), real64), int32)
        point%slot = column_material_slot(a, a_off)
        call column_state(u, u_off, constants, point%temperature, point%pressure)
        ! PETSc stores a scalar field's gradient at u_off_x + d.
        do d = 1, MIXED_DIMENSION
            point%temperature_gradient(d) = real(constants(CC_TEMPERATURE_SCALE), real64) &
                                            *real(u_x(u_off_x(1) + d), real64)
            point%pressure_gradient(d) = real(constants(CC_PRESSURE_SCALE), real64) &
                                         *real(u_x(u_off_x(field + 1) + d), real64)
        end do
        point%mobility = invert_inverse_mobility(column_inverse_mobility(a, a_off))

        ok = .false.
        call eval_frozen_response(point%slot, point%temperature, point%pressure, point%response, ok)
        if (ok) call eval_frozen_tangents(point%slot, point%temperature, point%pressure, point%tangents, ok)
        if (.not. ok) then
            call note_column_fault(point%slot, point%temperature, point%pressure)
            return
        end if

        gravity(1) = real(constants(CC_GRAVITY_X), real64)
        gravity(2) = real(constants(CC_GRAVITY_Z), real64)
        point%driver = point%tangents%dliquid_dpressure*point%pressure_gradient &
                       + point%tangents%dliquid_dtemperature*point%temperature_gradient &
                       - real(constants(CC_DENSITY), real64)*gravity
    end subroutine read_primal_point

    !> The row scale \(L/q_{\mathrm{scale}}\) the mixed conservation row carries.
    pure real(real64) function primal_row_scale(constants)
        implicit none
        PetscScalar, intent(in) :: constants(*)

        primal_row_scale = real(constants(CC_LENGTH_SCALE), real64)/real(constants(CC_FLUX_SCALE), real64)
    end function primal_row_scale

    !> Storage and source, \(f_0=((W-W^n)/\Delta t-s)L/q_{\mathrm{scale}}\).
    subroutine f0_primal_conservation( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_column_primal_f0_conservation")
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
        real(real64) :: temperature, pressure, temperature_old, pressure_old, dt, source, storage_rate
        integer(int32) :: slot
        logical :: ok

        dt = real(constants(CC_DT), real64)
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

        f0(1) = real((storage_rate - source)*primal_row_scale(constants), PETSC_REAL_KIND)
    end subroutine f0_primal_conservation

    !> The weak divergence, \(f_1=-\mathbf q_lL/q_{\mathrm{scale}}=k_r\mathbf M(\nabla p_l-\rho_l\mathbf g)L/q_{\mathrm{scale}}\).
    subroutine f1_primal_conservation( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f1 &
        ) bind(C, name="ftcms_column_primal_f1_conservation")
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

        type(type_primal_point) :: point
        real(real64) :: flow(MIXED_DIMENSION)
        integer(int32) :: d
        logical :: ok

        call read_primal_point(u, u_off, u_x, u_off_x, a, a_off, constants, point, ok)
        if (.not. ok) then
            do d = 1, MIXED_DIMENSION
                f1(d) = real(ieee_value(1.0_real64, ieee_quiet_nan), PETSC_REAL_KIND)
            end do
            return
        end if

        flow = point%response%relative_conductivity*matmul(point%mobility, point%driver)*primal_row_scale(constants)
        do d = 1, MIXED_DIMENSION
            f1(d) = real(flow(d), PETSC_REAL_KIND)
        end do
    end subroutine f1_primal_conservation

    !> Central differences of the liquid-pressure tangents in one unknown.
    !>
    !> Returns \(\partial_v(\partial p_l/\partial p_w)\) and
    !> \(\partial_v(\partial p_l/\partial T)\) for \(v=p_w\) (`in_pressure`) or
    !> \(v=T\). A state the closure cannot answer leaves both at zero, which
    !> drops the curvature term from that Jacobian entry and nothing else.
    subroutine liquid_tangent_slopes(point, in_pressure, slope_pressure, slope_temperature)
        implicit none
        type(type_primal_point), intent(in) :: point
        logical, intent(in) :: in_pressure
        !> \(\partial_v(\partial p_l/\partial p_w)\) and \(\partial_v(\partial p_l/\partial T)\).
        real(real64), intent(inout) :: slope_pressure, slope_temperature

        type(type_mixed_frozen_tangents) :: ahead, behind
        real(real64) :: step
        logical :: ok

        slope_pressure = 0.0_real64
        slope_temperature = 0.0_real64
        ok = .false.
        if (in_pressure) then
            step = HESSIAN_PRESSURE_STEP
            call eval_frozen_tangents(point%slot, point%temperature, point%pressure + step, ahead, ok)
            if (ok) call eval_frozen_tangents(point%slot, point%temperature, point%pressure - step, behind, ok)
        else
            step = HESSIAN_TEMPERATURE_STEP
            call eval_frozen_tangents(point%slot, point%temperature + step, point%pressure, ahead, ok)
            if (ok) call eval_frozen_tangents(point%slot, point%temperature - step, point%pressure, behind, ok)
        end if
        if (.not. ok) return

        slope_pressure = (ahead%dliquid_dpressure - behind%dliquid_dpressure)/(2.0_real64*step)
        slope_temperature = (ahead%dliquid_dtemperature - behind%dliquid_dtemperature)/(2.0_real64*step)
    end subroutine liquid_tangent_slopes

    !> \(\partial f_1^p/\partial u_p\): the pressure dependence of \(k_r\) and of the driver.
    !>
    !> \[
    !>   \frac{\partial f_1}{\partial p_w}=k_r\Bigl(\frac{\partial\ln k_r}{\partial p_w}\mathbf M\mathbf d
    !>   +\mathbf M\bigl(\partial_{p_w}p_{l,p}\nabla p_w+\partial_{p_w}p_{l,T}\nabla T\bigr)\Bigr)
    !>   \frac{p_{\mathrm{scale}}L}{q_{\mathrm{scale}}}.
    !> \]
    subroutine g2_primal_pressure_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_column_primal_g2_pp")
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

        call primal_value_block(u, u_off, u_x, u_off_x, a, a_off, constants, .true., g2)
    end subroutine g2_primal_pressure_pressure

    !> \(\partial f_1^p/\partial u_T\): the temperature dependence of \(k_r\) and of the driver.
    subroutine g2_primal_pressure_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_column_primal_g2_pT")
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

        call primal_value_block(u, u_off, u_x, u_off_x, a, a_off, constants, .false., g2)
    end subroutine g2_primal_pressure_temperature

    !> The value derivative of \(f_1^p\) in \(p_w\) or in \(T\), shared by the two \(g_2\) kernels.
    subroutine primal_value_block(u, u_off, u_x, u_off_x, a, a_off, constants, in_pressure, g2)
        implicit none
        PetscScalar, intent(in) :: u(*), u_x(*)
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        PetscScalar, intent(in) :: constants(*)
        logical, intent(in) :: in_pressure
        PetscScalar, intent(inout) :: g2(*)

        type(type_primal_point) :: point
        real(real64) :: slope_pressure, slope_temperature, dlog_conductivity, unit
        real(real64) :: change(MIXED_DIMENSION), block(MIXED_DIMENSION)
        integer(int32) :: d
        logical :: ok

        do d = 1, MIXED_DIMENSION
            g2(d) = 0.0_PETSC_REAL_KIND
        end do
        call read_primal_point(u, u_off, u_x, u_off_x, a, a_off, constants, point, ok)
        if (.not. ok) return

        call liquid_tangent_slopes(point, in_pressure, slope_pressure, slope_temperature)
        if (in_pressure) then
            dlog_conductivity = point%tangents%dlog_conductivity_dpressure
            unit = real(constants(CC_PRESSURE_SCALE), real64)
        else
            dlog_conductivity = point%tangents%dlog_conductivity_dtemperature
            unit = real(constants(CC_TEMPERATURE_SCALE), real64)
        end if

        ! d(k_r M d)/dv = k_r (dln k_r/dv M d + M dd/dv), dd/dv = p_l,pv grad p + p_l,Tv grad T.
        change = slope_pressure*point%pressure_gradient + slope_temperature*point%temperature_gradient
        block = point%response%relative_conductivity &
                *(dlog_conductivity*matmul(point%mobility, point%driver) + matmul(point%mobility, change)) &
                *unit*primal_row_scale(constants)
        do d = 1, MIXED_DIMENSION
            g2(d) = real(block(d), PETSC_REAL_KIND)
        end do
    end subroutine primal_value_block

    !> \(\partial f_1^p/\partial\nabla u_p=k_r(\partial p_l/\partial p_w)\mathbf Mp_{\mathrm{scale}}L/q_{\mathrm{scale}}\).
    subroutine g3_primal_pressure_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g3 &
        ) bind(C, name="ftcms_column_primal_g3_pp")
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

        call primal_gradient_block(u, u_off, u_x, u_off_x, a, a_off, constants, .true., g3)
    end subroutine g3_primal_pressure_pressure

    !> \(\partial f_1^p/\partial\nabla u_T=k_r(\partial p_l/\partial T)\mathbf MT_{\mathrm{scale}}L/q_{\mathrm{scale}}\).
    !>
    !> The cryosuction driver: the liquid pressure falls with the temperature,
    !> so a temperature gradient alone moves water.
    subroutine g3_primal_pressure_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g3 &
        ) bind(C, name="ftcms_column_primal_g3_pT")
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

        call primal_gradient_block(u, u_off, u_x, u_off_x, a, a_off, constants, .false., g3)
    end subroutine g3_primal_pressure_temperature

    !> The gradient derivative of \(f_1^p\) in \(\nabla p_w\) or in \(\nabla T\), shared by the two \(g_3\) kernels.
    subroutine primal_gradient_block(u, u_off, u_x, u_off_x, a, a_off, constants, in_pressure, g3)
        implicit none
        PetscScalar, intent(in) :: u(*), u_x(*)
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        PetscScalar, intent(in) :: constants(*)
        logical, intent(in) :: in_pressure
        PetscScalar, intent(inout) :: g3(*)

        type(type_primal_point) :: point
        real(real64) :: weight
        integer(int32) :: c, d
        logical :: ok

        do c = 1, MIXED_DIMENSION*MIXED_DIMENSION
            g3(c) = 0.0_PETSC_REAL_KIND
        end do
        call read_primal_point(u, u_off, u_x, u_off_x, a, a_off, constants, point, ok)
        if (.not. ok) return

        if (in_pressure) then
            weight = point%tangents%dliquid_dpressure*real(constants(CC_PRESSURE_SCALE), real64)
        else
            weight = point%tangents%dliquid_dtemperature*real(constants(CC_TEMPERATURE_SCALE), real64)
        end if
        weight = weight*point%response%relative_conductivity*primal_row_scale(constants)

        ! PETSc's g3 for one scalar row and one scalar column: df1_c/d(du/dx_d) at (c - 1)*dim + d.
        do c = 1, MIXED_DIMENSION
            do d = 1, MIXED_DIMENSION
                g3((c - 1)*MIXED_DIMENSION + d) = real(weight*point%mobility(c, d), PETSC_REAL_KIND)
            end do
        end do
    end subroutine primal_gradient_block

end module app_mixed_hydraulic_column_primal
