#include <petsc/finclude/petscsys.h>

!> The water row with the pore pressure in a Lagrange space.
!>
!> The same storage rate as `app_ftcms_kernels`, with
!> the pore pressure \(p_w\) in \(H^1\) (P1 on a triangle, Q1 on a rectangle)
!> instead of the RT0 flux and P0 pressure pair:
!>
!> \[
!>   R_p(\omega)=\bigl(\dot W(T,p_w)-s,\omega\bigr)
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
!> \(L/q_{\mathrm{scale}}\) and the Newton solver are shared, so a difference
!> between two runs is a difference of discretisation.
!>
!> ## What the primal form needs that the mixed one does not
!>
!> The mixed Darcy row tests the liquid pressure against \(\nabla\cdot\mathbf
!> v\), so \(p_l\) enters as a value and only its first derivatives reach the
!> Jacobian. Here \(\nabla p_l\) is in the residual, so the Jacobian needs the
!> derivatives of \(\partial p_l/\partial p_w\) and \(\partial p_l/\partial T\).
!> The closure publishes the first derivatives only; the second ones are
!> central differences of them.
module app_ftcms_kernels_lagrange
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use :: petscsys
    use :: app_ftcms_darcy, only:PLANE_DIMENSION
    use :: app_ftcms_closure, only:type_closure_response, type_closure_tangents, &
                                            eval_closure_response, eval_closure_tangents, eval_liquid_slopes
    use :: app_ftcms_kernels, only:point_state, point_state_old, material_slot_of, storage_rate_of, STORE_WATER, &
                                                      inverse_mobility_of, note_kernel_fault, &
                                                      AUX_SOURCE, CONST_DT, CONST_FLUX_SCALE, CONST_PRESSURE_SCALE, &
                                                      CONST_LENGTH_SCALE, CONST_DENSITY, CONST_GRAVITY_X, CONST_GRAVITY_Z, &
                                                      CONST_TEMPERATURE_SCALE, CONST_PRESSURE_FIELD, &
                                                      CONST_TEMPERATURE_FIELD, field_number
    implicit none
    private

    public :: f0_lagrange_conservation, f1_lagrange_conservation
    public :: g2_lagrange_pressure_pressure, g3_lagrange_pressure_pressure
    public :: g2_lagrange_pressure_temperature, g3_lagrange_pressure_temperature
    public :: compute_lagrange_darcy_flux

    ! Offsets of the central differences of the liquid-pressure tangents: the
    ! phase split publishes their values but not their own derivatives.
    real(real64), parameter :: HESSIAN_PRESSURE_STEP = 1.0d0
    real(real64), parameter :: HESSIAN_TEMPERATURE_STEP = 1.0d-4

    !> Everything one quadrature point of the water row needs.
    type :: type_lagrange_point
        integer(int32) :: slot = 0
        !> Temperature, C, and pore pressure, Pa.
        real(real64) :: temperature = 0.0d0
        real(real64) :: pressure = 0.0d0
        !> Physical gradients of the two unknowns, K/m and Pa/m.
        real(real64) :: temperature_gradient(PLANE_DIMENSION) = 0.0d0
        real(real64) :: pressure_gradient(PLANE_DIMENSION) = 0.0d0
        !> Saturated mobility of the cell, m2/(Pa s), row major.
        real(real64) :: mobility(PLANE_DIMENSION, PLANE_DIMENSION) = 0.0d0
        !> \(\nabla p_l-\rho_l\mathbf g\), Pa/m.
        real(real64) :: driver(PLANE_DIMENSION) = 0.0d0
        type(type_closure_response) :: response
        type(type_closure_tangents) :: tangents
    end type type_lagrange_point

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
    subroutine compute_lagrange_darcy_flux(slot, temperature, pressure, temperature_gradient, pressure_gradient, &
                                         inverse_mobility, density, gravity, flux, ok)
        implicit none
        !> One-based slot in the frozen closure table.
        integer(int32), intent(in) :: slot
        !> Temperature, C, and pore pressure \(p_w\), Pa.
        real(real64), intent(in) :: temperature, pressure
        !> \(\nabla T\), K/m, and \(\nabla p_w\), Pa/m.
        real(real64), intent(in) :: temperature_gradient(PLANE_DIMENSION)
        real(real64), intent(in) :: pressure_gradient(PLANE_DIMENSION)
        !> Saturated inverse mobility of the cell, Pa s/m2, row major.
        real(real64), intent(in) :: inverse_mobility(PLANE_DIMENSION*PLANE_DIMENSION)
        !> Reference liquid density, kg/m3, and gravity, m/s2.
        real(real64), intent(in) :: density
        real(real64), intent(in) :: gravity(PLANE_DIMENSION)
        !> \(\mathbf q_l\), m/s; overwritten when `ok`.
        real(real64), intent(inout) :: flux(PLANE_DIMENSION)
        !> True when the closure answered.
        logical, intent(inout) :: ok

        type(type_closure_response) :: response
        real(real64) :: mobility(PLANE_DIMENSION, PLANE_DIMENSION), driver(PLANE_DIMENSION)
        real(real64) :: dliquid_dpressure, dliquid_dtemperature

        ok = .false.
        call eval_closure_response(slot, temperature, pressure, response, ok)
        if (ok) call eval_liquid_slopes(slot, temperature, pressure, dliquid_dpressure, dliquid_dtemperature, ok)
        if (.not. ok) return

        mobility = invert_inverse_mobility(inverse_mobility)
        driver = dliquid_dpressure*pressure_gradient + dliquid_dtemperature*temperature_gradient - density*gravity
        flux = -response%relative_conductivity*matmul(mobility, driver)
    end subroutine compute_lagrange_darcy_flux

    !> \(\mathbf M=(\mathbf\Lambda^{-1})^{-1}\) of a row-major 2x2 tensor.
    pure function invert_inverse_mobility(inverse_mobility) result(mobility)
        implicit none
        real(real64), intent(in) :: inverse_mobility(PLANE_DIMENSION*PLANE_DIMENSION)
        real(real64) :: mobility(PLANE_DIMENSION, PLANE_DIMENSION)

        real(real64) :: determinant

        determinant = inverse_mobility(1)*inverse_mobility(4) - inverse_mobility(2)*inverse_mobility(3)
        mobility(1, 1) = inverse_mobility(4)/determinant
        mobility(1, 2) = -inverse_mobility(2)/determinant
        mobility(2, 1) = -inverse_mobility(3)/determinant
        mobility(2, 2) = inverse_mobility(1)/determinant
    end function invert_inverse_mobility

    !> Gather the state, gradients, mobility and closure of one quadrature point.
    subroutine read_lagrange_point(u, u_off, u_x, u_off_x, a, a_off, constants, point, ok, slopes_only)
        implicit none
        PetscScalar, intent(in) :: u(*), u_x(*)
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        PetscScalar, intent(in) :: constants(*)
        type(type_lagrange_point), intent(inout) :: point
        logical, intent(inout) :: ok
        !> Fill only the liquid-pressure slopes of `point%tangents`, as a
        !> residual needs; the Jacobian blocks ask for every tangent.
        logical, intent(in), optional :: slopes_only

        real(real64) :: gravity(PLANE_DIMENSION)
        integer(int32) :: d, field

        field = field_number(constants, CONST_PRESSURE_FIELD)
        point%slot = material_slot_of(a, a_off)
        call point_state(u, u_off, constants, point%temperature, point%pressure)
        do d = 1, PLANE_DIMENSION
            point%temperature_gradient(d) = 0.0d0
            if (field_number(constants, CONST_TEMPERATURE_FIELD) >= 0) then
                point%temperature_gradient(d) = constants(CONST_TEMPERATURE_SCALE) &
                                                *u_x(u_off_x(field_number(constants, CONST_TEMPERATURE_FIELD) + 1) + d)
            end if
            point%pressure_gradient(d) = constants(CONST_PRESSURE_SCALE) &
                                         *u_x(u_off_x(field + 1) + d)
        end do
        point%mobility = invert_inverse_mobility(inverse_mobility_of(a, a_off))

        ok = .false.
        call eval_closure_response(point%slot, point%temperature, point%pressure, point%response, ok)
        if (ok) then
            if (present(slopes_only)) then
                if (slopes_only) then
                    call eval_liquid_slopes(point%slot, point%temperature, point%pressure, &
                                            point%tangents%dliquid_dpressure, point%tangents%dliquid_dtemperature, ok)
                else
                    call eval_closure_tangents(point%slot, point%temperature, point%pressure, point%tangents, ok)
                end if
            else
                call eval_closure_tangents(point%slot, point%temperature, point%pressure, point%tangents, ok)
            end if
        end if
        if (.not. ok) then
            call note_kernel_fault(point%slot, point%temperature, point%pressure)
            return
        end if

        gravity(1) = constants(CONST_GRAVITY_X)
        gravity(2) = constants(CONST_GRAVITY_Z)
        point%driver = point%tangents%dliquid_dpressure*point%pressure_gradient &
                       + point%tangents%dliquid_dtemperature*point%temperature_gradient &
                       - constants(CONST_DENSITY)*gravity
    end subroutine read_lagrange_point

    !> The row scale \(L/q_{\mathrm{scale}}\) the mixed conservation row carries.
    pure real(real64) function lagrange_row_scale(constants)
        implicit none
        PetscScalar, intent(in) :: constants(*)

        lagrange_row_scale = constants(CONST_LENGTH_SCALE)/constants(CONST_FLUX_SCALE)
    end function lagrange_row_scale

    !> Storage and source, \(f_0=((W-W^n)/\Delta t-s)L/q_{\mathrm{scale}}\).
    subroutine f0_lagrange_conservation( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_lagrange_f0_conservation")
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
        real(real64) :: temperature, pressure, temperature_old, pressure_old, dt, source, storage_rate
        integer(int32) :: slot
        logical :: ok

        dt = constants(CONST_DT)
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

        f0(1) = (storage_rate - source)*lagrange_row_scale(constants)
    end subroutine f0_lagrange_conservation

    !> The weak divergence, \(f_1=-\mathbf q_lL/q_{\mathrm{scale}}=k_r\mathbf M(\nabla p_l-\rho_l\mathbf g)L/q_{\mathrm{scale}}\).
    subroutine f1_lagrange_conservation( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f1 &
        ) bind(C, name="ftcms_lagrange_f1_conservation")
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

        type(type_lagrange_point) :: point
        real(real64) :: flow(PLANE_DIMENSION)
        integer(int32) :: d
        logical :: ok

        call read_lagrange_point(u, u_off, u_x, u_off_x, a, a_off, constants, point, ok, slopes_only=.true.)
        if (.not. ok) then
            do d = 1, PLANE_DIMENSION
                f1(d) = ieee_value(1.0d0, ieee_quiet_nan)
            end do
            return
        end if

        flow = point%response%relative_conductivity*matmul(point%mobility, point%driver)*lagrange_row_scale(constants)
        do d = 1, PLANE_DIMENSION
            f1(d) = flow(d)
        end do
    end subroutine f1_lagrange_conservation

    !> Central differences of the liquid-pressure tangents in one unknown.
    !>
    !> Returns \(\partial_v(\partial p_l/\partial p_w)\) and
    !> \(\partial_v(\partial p_l/\partial T)\) for \(v=p_w\) (`in_pressure`) or
    !> \(v=T\). A state the closure cannot answer leaves both at zero, which
    !> drops the curvature term from that Jacobian entry and nothing else.
    subroutine liquid_tangent_slopes(point, in_pressure, slope_pressure, slope_temperature)
        implicit none
        type(type_lagrange_point), intent(in) :: point
        logical, intent(in) :: in_pressure
        !> \(\partial_v(\partial p_l/\partial p_w)\) and \(\partial_v(\partial p_l/\partial T)\).
        real(real64), intent(inout) :: slope_pressure, slope_temperature

        real(real64) :: step, ahead(2), behind(2)
        logical :: ok

        slope_pressure = 0.0d0
        slope_temperature = 0.0d0
        ok = .false.
        if (in_pressure) then
            step = HESSIAN_PRESSURE_STEP
            call eval_liquid_slopes(point%slot, point%temperature, point%pressure + step, ahead(1), ahead(2), ok)
            if (ok) call eval_liquid_slopes(point%slot, point%temperature, point%pressure - step, &
                                            behind(1), behind(2), ok)
        else
            step = HESSIAN_TEMPERATURE_STEP
            call eval_liquid_slopes(point%slot, point%temperature + step, point%pressure, ahead(1), ahead(2), ok)
            if (ok) call eval_liquid_slopes(point%slot, point%temperature - step, point%pressure, &
                                            behind(1), behind(2), ok)
        end if
        if (.not. ok) return

        slope_pressure = (ahead(1) - behind(1))/(2.0d0*step)
        slope_temperature = (ahead(2) - behind(2))/(2.0d0*step)
    end subroutine liquid_tangent_slopes

    !> \(\partial f_1^p/\partial u_p\): the pressure dependence of \(k_r\) and of the driver.
    !>
    !> \[
    !>   \frac{\partial f_1}{\partial p_w}=k_r\Bigl(\frac{\partial\ln k_r}{\partial p_w}\mathbf M\mathbf d
    !>   +\mathbf M\bigl(\partial_{p_w}p_{l,p}\nabla p_w+\partial_{p_w}p_{l,T}\nabla T\bigr)\Bigr)
    !>   \frac{p_{\mathrm{scale}}L}{q_{\mathrm{scale}}}.
    !> \]
    subroutine g2_lagrange_pressure_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_lagrange_g2_pp")
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

        call lagrange_value_block(u, u_off, u_x, u_off_x, a, a_off, constants, .true., g2)
    end subroutine g2_lagrange_pressure_pressure

    !> \(\partial f_1^p/\partial u_T\): the temperature dependence of \(k_r\) and of the driver.
    subroutine g2_lagrange_pressure_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_lagrange_g2_pT")
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

        call lagrange_value_block(u, u_off, u_x, u_off_x, a, a_off, constants, .false., g2)
    end subroutine g2_lagrange_pressure_temperature

    !> The value derivative of \(f_1^p\) in \(p_w\) or in \(T\), shared by the two \(g_2\) kernels.
    subroutine lagrange_value_block(u, u_off, u_x, u_off_x, a, a_off, constants, in_pressure, g2)
        implicit none
        PetscScalar, intent(in) :: u(*), u_x(*)
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        PetscScalar, intent(in) :: constants(*)
        logical, intent(in) :: in_pressure
        PetscScalar, intent(inout) :: g2(*)

        type(type_lagrange_point) :: point
        real(real64) :: slope_pressure, slope_temperature, dlog_conductivity, unit
        real(real64) :: change(PLANE_DIMENSION), block(PLANE_DIMENSION)
        integer(int32) :: d
        logical :: ok

        do d = 1, PLANE_DIMENSION
            g2(d) = 0.0d0
        end do
        call read_lagrange_point(u, u_off, u_x, u_off_x, a, a_off, constants, point, ok)
        if (.not. ok) return

        call liquid_tangent_slopes(point, in_pressure, slope_pressure, slope_temperature)
        if (in_pressure) then
            dlog_conductivity = point%tangents%dlog_conductivity_dpressure
            unit = constants(CONST_PRESSURE_SCALE)
        else
            dlog_conductivity = point%tangents%dlog_conductivity_dtemperature
            unit = constants(CONST_TEMPERATURE_SCALE)
        end if

        ! d(k_r M d)/dv = k_r (dln k_r/dv M d + M dd/dv), dd/dv = p_l,pv grad p + p_l,Tv grad T.
        change = slope_pressure*point%pressure_gradient + slope_temperature*point%temperature_gradient
        block = point%response%relative_conductivity &
                *(dlog_conductivity*matmul(point%mobility, point%driver) + matmul(point%mobility, change)) &
                *unit*lagrange_row_scale(constants)
        do d = 1, PLANE_DIMENSION
            g2(d) = block(d)
        end do
    end subroutine lagrange_value_block

    !> \(\partial f_1^p/\partial\nabla u_p=k_r(\partial p_l/\partial p_w)\mathbf Mp_{\mathrm{scale}}L/q_{\mathrm{scale}}\).
    subroutine g3_lagrange_pressure_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g3 &
        ) bind(C, name="ftcms_lagrange_g3_pp")
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

        call lagrange_gradient_block(u, u_off, u_x, u_off_x, a, a_off, constants, .true., g3)
    end subroutine g3_lagrange_pressure_pressure

    !> \(\partial f_1^p/\partial\nabla u_T=k_r(\partial p_l/\partial T)\mathbf MT_{\mathrm{scale}}L/q_{\mathrm{scale}}\).
    !>
    !> The cryosuction driver: the liquid pressure falls with the temperature,
    !> so a temperature gradient alone moves water.
    subroutine g3_lagrange_pressure_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g3 &
        ) bind(C, name="ftcms_lagrange_g3_pT")
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

        call lagrange_gradient_block(u, u_off, u_x, u_off_x, a, a_off, constants, .false., g3)
    end subroutine g3_lagrange_pressure_temperature

    !> The gradient derivative of \(f_1^p\) in \(\nabla p_w\) or in \(\nabla T\), shared by the two \(g_3\) kernels.
    subroutine lagrange_gradient_block(u, u_off, u_x, u_off_x, a, a_off, constants, in_pressure, g3)
        implicit none
        PetscScalar, intent(in) :: u(*), u_x(*)
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        PetscScalar, intent(in) :: constants(*)
        logical, intent(in) :: in_pressure
        PetscScalar, intent(inout) :: g3(*)

        type(type_lagrange_point) :: point
        real(real64) :: weight
        integer(int32) :: c, d
        logical :: ok

        do c = 1, PLANE_DIMENSION*PLANE_DIMENSION
            g3(c) = 0.0d0
        end do
        call read_lagrange_point(u, u_off, u_x, u_off_x, a, a_off, constants, point, ok)
        if (.not. ok) return

        if (in_pressure) then
            weight = point%tangents%dliquid_dpressure*constants(CONST_PRESSURE_SCALE)
        else
            weight = point%tangents%dliquid_dtemperature*constants(CONST_TEMPERATURE_SCALE)
        end if
        weight = weight*point%response%relative_conductivity*lagrange_row_scale(constants)

        ! PETSc's g3 for one scalar row and one scalar column: df1_c/d(du/dx_d) at (c - 1)*dim + d.
        do c = 1, PLANE_DIMENSION
            do d = 1, PLANE_DIMENSION
                g3((c - 1)*PLANE_DIMENSION + d) = weight*point%mobility(c, d)
            end do
        end do
    end subroutine lagrange_gradient_block

end module app_ftcms_kernels_lagrange
