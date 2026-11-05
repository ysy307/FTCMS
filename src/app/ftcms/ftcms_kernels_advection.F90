#include <petsc/finclude/petscsys.h>

!> Sensible heat carried by the liquid flux, as PetscDS pointwise functions.
!>
!> \[
!>   R_T(w)\mathrel{+}=\bigl(C_l\,\mathbf q_l\cdot\nabla T,\ w\bigr)/R_T,
!>   \qquad C_l=\rho_w c_w,
!> \]
!>
!> the advective form of the legacy heat row. The flux is the one the water
!> equation's space defines: the RT unknown of a `mixed` field, and Darcy's law
!> from the pressure gradient for a `lagrange` field,
!>
!> \[
!>   \mathbf q_l=-k_r\,\Lambda\Bigl(\frac{\partial p_l}{\partial p_w}\nabla p_w
!>   +\frac{\partial p_l}{\partial T}\nabla T-\rho_l\mathbf g\Bigr),
!> \]
!>
!> evaluated at the quadrature point. A `nodal` water field uses the same
!> pointwise flux here: its element-mean conductivity enters the water row only.
!>
!> The Jacobian carries the dependence on \(\nabla T\), on the RT unknown and on
!> \(\nabla p_w\). The dependence of \(C_l\), \(k_r\) and the liquid-pressure
!> slopes on the state is not differentiated: it is second order in the
!> advection, which is small against conduction in these soils.
module app_ftcms_kernels_advection
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: petscsys
    use :: app_ftcms_darcy, only:PLANE_DIMENSION
    use :: app_ftcms_closure, only:type_closure_response, &
                                   eval_closure_response, eval_liquid_slopes
    use :: app_ftcms_kernels, only:CONST_FLUX_SCALE, CONST_PRESSURE_SCALE, CONST_TEMPERATURE_SCALE, &
                                   CONST_THERMAL_SCALE, CONST_DENSITY, CONST_GRAVITY_X, CONST_GRAVITY_Z, &
                                   CONST_TEMPERATURE_FIELD, CONST_FLUX_FIELD, CONST_PRESSURE_FIELD, &
                                   field_number, point_state, material_slot_of, inverse_mobility_of, &
                                   f0_temperature, note_kernel_fault
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    private

    public :: f0_temperature_advected
    public :: g1_temperature_temperature
    public :: g0_temperature_flux
    public :: g1_temperature_pressure

    !> The liquid flux at one point and its slopes.
    type :: type_advection_point
        !> \(C_l\), J/(m3 K).
        real(real64) :: capacity = 0.0d0
        !> \(\nabla T\), K/m, and \(\mathbf q_l\), m/s.
        real(real64) :: temperature_gradient(PLANE_DIMENSION) = 0.0d0
        real(real64) :: flux(PLANE_DIMENSION) = 0.0d0
        !> \(\partial q_i/\partial(\nabla T)_d\) and \(\partial q_i/\partial(\nabla p_w)_d\).
        real(real64) :: dflux_dtemperature_gradient(PLANE_DIMENSION, PLANE_DIMENSION) = 0.0d0
        real(real64) :: dflux_dpressure_gradient(PLANE_DIMENSION, PLANE_DIMENSION) = 0.0d0
    end type type_advection_point

contains

    !> Heat storage plus \(C_l\,\mathbf q_l\cdot\nabla T/R_T\): the heat row's \(f_0\) when water moves.
    subroutine f0_temperature_advected( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_f0_temperature_advected")
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

        real(real64) :: advection, temperature, pressure
        logical :: ok

        call f0_temperature(dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, a_off, a_off_x, a, a_t, a_x, time, x, &
                            num_constants, constants, f0)
        call advection_residual(u, u_off, u_x, u_off_x, a, a_off, constants, advection, ok)
        if (.not. ok) then
            call point_state(u, u_off, constants, temperature, pressure)
            call note_kernel_fault(material_slot_of(a, a_off), temperature, pressure)
            f0(1) = ieee_value(1.0d0, ieee_quiet_nan)
            return
        end if
        f0(1) = f0(1) + advection
    end subroutine f0_temperature_advected

    !> \(C_l\,\mathbf q_l\cdot\nabla T/R_T\) at one point; zero without a water field.
    !>
    !> `ok` is false when the closure cannot answer at the state.
    subroutine advection_residual(u, u_off, u_x, u_off_x, a, a_off, constants, value, ok)
        implicit none
        PetscScalar, intent(in) :: u(*), u_x(*)
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: value
        logical, intent(inout) :: ok

        type(type_advection_point) :: point

        value = 0.0d0
        call read_advection_point(u, u_off, u_x, u_off_x, a, a_off, constants, point, ok)
        if (.not. ok) return
        value = point%capacity*dot_product(point%flux, point%temperature_gradient) &
                /constants(CONST_THERMAL_SCALE)
    end subroutine advection_residual

    !> The flux, its slopes and \(C_l\) at one point.
    subroutine read_advection_point(u, u_off, u_x, u_off_x, a, a_off, constants, point, ok)
        implicit none
        PetscScalar, intent(in) :: u(*), u_x(*)
        PetscInt, intent(in) :: u_off(*), u_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        PetscScalar, intent(in) :: constants(*)
        type(type_advection_point), intent(inout) :: point
        logical, intent(inout) :: ok

        type(type_closure_response) :: response
        real(real64) :: temperature, pressure, pressure_gradient(PLANE_DIMENSION), gravity(PLANE_DIMENSION)
        real(real64) :: dliquid_dpressure, dliquid_dtemperature
        real(real64) :: mobility(PLANE_DIMENSION, PLANE_DIMENSION), inverse(PLANE_DIMENSION*PLANE_DIMENSION)
        real(real64) :: determinant
        integer(int32) :: slot, t, q, p, d

        point = type_advection_point()
        ok = .true.
        t = field_number(constants, CONST_TEMPERATURE_FIELD)
        q = field_number(constants, CONST_FLUX_FIELD)
        p = field_number(constants, CONST_PRESSURE_FIELD)
        if (t < 0 .or. p < 0) return

        slot = material_slot_of(a, a_off)
        call point_state(u, u_off, constants, temperature, pressure)
        ok = .false.
        call eval_closure_response(slot, temperature, pressure, response, ok)
        if (.not. ok) return
        point%capacity = response%liquid_heat_capacity
        do d = 1, PLANE_DIMENSION
            point%temperature_gradient(d) = constants(CONST_TEMPERATURE_SCALE) &
                                            *u_x(u_off_x(t + 1) + d)
        end do

        if (q >= 0) then
            do d = 1, PLANE_DIMENSION
                point%flux(d) = constants(CONST_FLUX_SCALE)*u(u_off(q + 1) + d)
            end do
            return
        end if

        call eval_liquid_slopes(slot, temperature, pressure, dliquid_dpressure, dliquid_dtemperature, ok)
        if (.not. ok) return
        do d = 1, PLANE_DIMENSION
            pressure_gradient(d) = constants(CONST_PRESSURE_SCALE)*u_x(u_off_x(p + 1) + d)
        end do
        gravity = [constants(CONST_GRAVITY_X), constants(CONST_GRAVITY_Z)]
        inverse = inverse_mobility_of(a, a_off)
        determinant = inverse(1)*inverse(4) - inverse(2)*inverse(3)
        mobility = response%relative_conductivity*reshape([inverse(4), -inverse(3), -inverse(2), inverse(1)], &
                                                          [PLANE_DIMENSION, PLANE_DIMENSION])/determinant
        point%flux = -matmul(mobility, dliquid_dpressure*pressure_gradient &
                             + dliquid_dtemperature*point%temperature_gradient &
                             - constants(CONST_DENSITY)*gravity)
        point%dflux_dpressure_gradient = -dliquid_dpressure*mobility
        point%dflux_dtemperature_gradient = -dliquid_dtemperature*mobility
    end subroutine read_advection_point

    !> \(\partial f_0/\partial(\nabla T)\): \(C_l(\mathbf q_l+(\partial\mathbf q_l/\partial\nabla T)^T\nabla T)\,T_s/R_T\).
    subroutine g1_temperature_temperature( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g1 &
        ) bind(C, name="ftcms_g1_temperature_temperature")
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

        type(type_advection_point) :: point
        real(real64) :: weight
        integer(int32) :: d
        logical :: ok

        g1(1:PLANE_DIMENSION) = 0.0d0
        call read_advection_point(u, u_off, u_x, u_off_x, a, a_off, constants, point, ok)
        if (.not. ok) return
        weight = point%capacity*constants(CONST_TEMPERATURE_SCALE) &
                 /constants(CONST_THERMAL_SCALE)
        do d = 1, PLANE_DIMENSION
            g1(d) = weight*(point%flux(d) + dot_product(point%temperature_gradient, point%dflux_dtemperature_gradient(:, d)))
        end do
    end subroutine g1_temperature_temperature

    !> \(\partial f_0/\partial\mathbf q\) of an RT flux: \(C_l\nabla T\,q_s/R_T\).
    subroutine g0_temperature_flux( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_g0_temperature_flux")
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

        type(type_advection_point) :: point
        integer(int32) :: c
        logical :: ok

        g0(1:PLANE_DIMENSION) = 0.0d0
        call read_advection_point(u, u_off, u_x, u_off_x, a, a_off, constants, point, ok)
        if (.not. ok) return
        do c = 1, PLANE_DIMENSION
            g0(c) = point%capacity*point%temperature_gradient(c)*constants(CONST_FLUX_SCALE) &
                    /constants(CONST_THERMAL_SCALE)
        end do
    end subroutine g0_temperature_flux

    !> \(\partial f_0/\partial(\nabla p_w)\) of a Lagrange pressure:
    !> \(C_l(\partial\mathbf q_l/\partial\nabla p_w)^T\nabla T\,p_s/R_T\).
    subroutine g1_temperature_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g1 &
        ) bind(C, name="ftcms_g1_temperature_pressure")
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

        type(type_advection_point) :: point
        integer(int32) :: d
        logical :: ok

        g1(1:PLANE_DIMENSION) = 0.0d0
        call read_advection_point(u, u_off, u_x, u_off_x, a, a_off, constants, point, ok)
        if (.not. ok) return
        do d = 1, PLANE_DIMENSION
            g1(d) = point%capacity*dot_product(point%temperature_gradient, point%dflux_dpressure_gradient(:, d)) &
                    *constants(CONST_PRESSURE_SCALE) &
                    /constants(CONST_THERMAL_SCALE)
        end do
    end subroutine g1_temperature_pressure

end module app_ftcms_kernels_advection
