#include <petsc/finclude/petscsys.h>

!> The weak form of the energy equation, as PetscDS pointwise functions.
!>
!> PETSc evaluates these at every quadrature point and does the integration,
!> the element loop and the assembly itself. Each is called through a C
!> function pointer, so the argument lists have to match PETSc's own prototypes
!> exactly; they are PetscPointFn, PetscPointJacFn, PetscBdPointFn and
!> PetscBdPointJacFn in petscdstypes.h.
!>
!> The residual solved at each step is
!>
!>   F(T) = A_i (U(T_i) - U(T_old,i))/dt + \int grad v . lambda(T) grad T
!>
!> with U the enthalpy density, which carries the latent heat of the phase
!> change. Differencing the enthalpy rather than writing C(T) dT/dt is what
!> makes this solvable: the apparent capacity C = dU/dT is a near-delta
!> function across the freezing band, so a residual built on it is not
!> differentiable there and Newton cycles between two states instead of
!> converging. U itself is monotone and Lipschitz, so its difference is smooth.
!>
!> Only the conduction and boundary terms are here. The storage term is lumped
!> -- weighted by the nodal volume A_i and evaluated at the node rather than at
!> the quadrature points -- and lumping is a nodal operation, not an
!> integration, so plexfem_run adds it to what PETSc assembles. Lumping it
!> matters: with the consistent form a strongly cooled boundary drives the
!> neighbouring nodes past the initial temperature.
!>
!> Time is discretised by the caller, not by TS, so u_t is never read here and
!> PETSc's u_tShift is always zero.
!>
!> The material of the cell arrives in an auxiliary field, so one set of
!> kernels serves a mesh of several materials.
!>
!> The Jacobian blocks kept here are
!>
!>   g2 = d(f1)/dT        = (dlambda/dT) grad T
!>   g3 = d(f1)/d(grad T) = lambda
!>
!> dlambda/dT has no closed form to reach from here and comes by central
!> difference.
module app_plexfem_pointwise
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_f_pointer
    use :: petscsys
    use :: app_plexfem_physics, only:eval_capacity, eval_conductivity, eval_enthalpy, &
                                     eval_conductivity_tangent, eval_conductivity_pressure_tangent, &
                                     prescribed_value, C_DT, &
                                     eval_water_driver, &
                                     eval_water_flux_tangent, &
                                     A_MATERIAL_ID, A_DARCY, A_COUPLING, A_GRAVITY, A_POTENTIAL, &
                                     A_TEMPERATURE_OLD, A_PRESSURE_OLD, &
                                     A_TEMPERATURE_OLDER, A_PRESSURE_OLDER, &
                                     eval_water_content, bdf_order, bdf_coefficient, &
                                     eval_capacity, eval_water_capacity, &
                                     eval_enthalpy_pressure_tangent, eval_water_temperature_tangent, &
                                     C_MEASURING_MATERIAL, C_PRESSURE, &
                                     A_BC_KIND, A_BC_VALUE_1, A_BC_VALUE_2, &
                                     BC_FLUX, BC_ROBIN
    implicit none
    private

    public :: f0_energy, f0_water, g0_energy, g0_energy_p, g0_water_p, g0_water_t
    public :: f1_energy, g2_energy, g2_energy_p, g3_energy, f0_material_indicator
    public :: f1_water, g3_water_pp, g3_water_pt, g2_water_p, g2_water_t
    public :: bd_f0_energy, bd_g0_energy
    public :: bc_essential_temperature, bc_essential_pressure

contains

    !> The material of the cell being integrated, from its auxiliary field.
    !>
    !> The field is one value per cell, so PETSc hands the same number back at
    !> every quadrature point of it. It is stored as a real, like every other
    !> auxiliary value, and rounded rather than truncated on the way back.
    pure function material_of(a, a_off) result(material_id)
        implicit none
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        integer(int32) :: material_id

        material_id = int(nint(real(a(a_off(A_MATERIAL_ID) + 1), real64)), int32)
    end function material_of

    !> The pore pressure at this quadrature point.
    !>
    !> It is a field of its own when the project solves for water, and a
    !> constant of the run when it does not. PETSc says which by the number of
    !> fields it hands the kernel, so one set of kernels serves both.
    pure function pressure_of(nf, u_off, u, constants) result(pressure)
        implicit none
        PetscInt, value :: nf
        PetscInt, intent(in) :: u_off(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: constants(*)
        real(real64) :: pressure

        if (nf >= 2) then
            pressure = real(u(u_off(2) + 1), real64)
        else
            pressure = real(constants(C_PRESSURE), real64)
        end if
    end function pressure_of

    !> An element-representative coefficient of this cell, from the auxiliary
    !> field the assembly refreshed.
    pure function element_of(a, a_off, field) result(value)
        implicit none
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        integer(int32), intent(in) :: field
        real(real64) :: value

        value = real(a(a_off(field) + 1), real64)
    end function element_of

    !> One where the cell belongs to the material being measured, zero
    !> elsewhere.
    !>
    !> Integrating this against the test functions gives \int psi_i over that
    !> material, which is the lumped nodal volume the storage term is weighted
    !> by. It runs once per material while the run is being set up, so that the
    !> weights come from the same quadrature PETSc integrates everything else
    !> with, rather than from a second implementation of the element geometry.
    subroutine f0_material_indicator( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_plexfem_f0_material_indicator")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f0(*)

        integer(int32) :: target_material

        target_material = int(nint(real(constants(C_MEASURING_MATERIAL), real64)), int32)

        if (material_of(a, a_off) == target_material) then
            f0(1) = 1.0_PETSC_REAL_KIND
        else
            f0(1) = 0.0_PETSC_REAL_KIND
        end if
    end subroutine f0_material_indicator

    !> Storage term of the energy equation, integrated at the quadrature points.
    !>
    !>   f0 = sum_j c_j U(T_j, p_j)
    !>
    !> with j the history levels the BDF coefficients weight, level zero being the
    !> iterate and the rest arriving as auxiliary fields PETSc interpolates.
    !>
    !> Consistent, not lumped -- and that changes the discrete equations, not
    !> only the iteration: what a lumped residual conserves is the nodal sum and
    !> what this conserves is the integral, and on a finite mesh they differ.
    !>
    !> The reason to prefer it is the tangent. The derivative of
    !> int psi_i c_0 U dOmega is c_0 int psi_i C_a psi_j dOmega, while the lumped
    !> form diag(int psi_i C_a dOmega) approximates that only where C_a is smooth
    !> over the element. At a freezing front C_a is dominated by
    !> -L_f rho_i dtheta_i/dT and varies by orders of magnitude between adjacent
    !> nodes, so the lumped tangent drives the contraction factor toward one. The
    !> older path measured exactly that: the thermal residual stalled at eleven
    !> percent of its step-initial value for thirty iterations.
    subroutine f0_energy( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_plexfem_f0_energy")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f0(*)

        real(real64) :: stored, rate
        integer(int32) :: material_id, level

        material_id = material_of(a, a_off)
        rate = 0.0d0
        do level = 0, bdf_order(constants)
            call eval_enthalpy(level_temperature(nf, u_off, u, a, a_off, level), &
                               level_pressure(nf, u_off, u, a, a_off, constants, level), &
                               material_id, constants, stored)
            rate = rate + bdf_coefficient(constants, level) * stored
        end do
        f0(1) = rate
    end subroutine f0_energy

    !> Storage term of the water equation, the same way.
    subroutine f0_water( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_plexfem_f0_water")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f0(*)

        real(real64) :: stored, rate
        integer(int32) :: material_id, level

        material_id = material_of(a, a_off)
        rate = 0.0d0
        do level = 0, bdf_order(constants)
            call eval_water_content(level_temperature(nf, u_off, u, a, a_off, level), &
                                    level_pressure(nf, u_off, u, a, a_off, constants, level), &
                                    material_id, constants, stored)
            rate = rate + bdf_coefficient(constants, level) * stored
        end do
        f0(1) = rate
    end subroutine f0_water

    !> Temperature at this quadrature point, at one history level.
    pure function level_temperature(nf, u_off, u, a, a_off, level) result(temperature)
        implicit none
        PetscInt, value :: nf
        PetscInt, intent(in) :: u_off(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        integer(int32), intent(in) :: level
        real(real64) :: temperature

        select case (level)
        case (0)
            temperature = real(u(u_off(1) + 1), real64)
        case (1)
            temperature = real(a(a_off(A_TEMPERATURE_OLD) + 1), real64)
        case default
            temperature = real(a(a_off(A_TEMPERATURE_OLDER) + 1), real64)
        end select
    end function level_temperature

    !> Pore pressure at this quadrature point, at one history level.
    pure function level_pressure(nf, u_off, u, a, a_off, constants, level) result(pressure)
        implicit none
        PetscInt, value :: nf
        PetscInt, intent(in) :: u_off(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        PetscScalar, intent(in) :: constants(*)
        integer(int32), intent(in) :: level
        real(real64) :: pressure

        if (nf < 2) then
            pressure = real(constants(C_PRESSURE), real64)
            return
        end if

        select case (level)
        case (0)
            pressure = real(u(u_off(2) + 1), real64)
        case (1)
            pressure = real(a(a_off(A_PRESSURE_OLD) + 1), real64)
        case default
            pressure = real(a(a_off(A_PRESSURE_OLDER) + 1), real64)
        end select
    end function level_pressure

    !> d(f0_energy)/dT = c_0 dU/dT, the apparent heat capacity of the step.
    subroutine g0_energy( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_plexfem_g0_energy")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, value :: u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g0(*)

        real(real64) :: tangent

        call eval_capacity(real(u(u_off(1) + 1), real64), &
                         pressure_of(nf, u_off, u, constants), &
                         material_of(a, a_off), constants, tangent)
        ! Only the current level of the storage term depends on the unknown, so
        ! its coefficient is the whole of the derivative.
        g0(1) = bdf_coefficient(constants, 0) * tangent
    end subroutine g0_energy

    !> d(f0_energy)/dp = c_0 dU/dp. Pressure moves water across the freezing curve, so this carries latent heat.
    subroutine g0_energy_p( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_plexfem_g0_energy_p")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, value :: u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g0(*)

        real(real64) :: tangent

        call eval_enthalpy_pressure_tangent(real(u(u_off(1) + 1), real64), &
                         pressure_of(nf, u_off, u, constants), &
                         material_of(a, a_off), constants, tangent)
        ! Only the current level of the storage term depends on the unknown, so
        ! its coefficient is the whole of the derivative.
        g0(1) = bdf_coefficient(constants, 0) * tangent
    end subroutine g0_energy_p

    !> d(f0_water)/dp = c_0 dTheta/dp, the specific moisture capacity of the step.
    subroutine g0_water_p( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_plexfem_g0_water_p")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, value :: u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g0(*)

        real(real64) :: tangent

        call eval_water_capacity(real(u(u_off(1) + 1), real64), &
                         pressure_of(nf, u_off, u, constants), &
                         material_of(a, a_off), constants, tangent)
        ! Only the current level of the storage term depends on the unknown, so
        ! its coefficient is the whole of the derivative.
        g0(1) = bdf_coefficient(constants, 0) * tangent
    end subroutine g0_water_p

    !> d(f0_water)/dT = c_0 dTheta/dT. Zero in the closure active today, but asked of the model rather than assumed.
    subroutine g0_water_t( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_plexfem_g0_water_t")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, value :: u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g0(*)

        real(real64) :: tangent

        call eval_water_temperature_tangent(real(u(u_off(1) + 1), real64), &
                         pressure_of(nf, u_off, u, constants), &
                         material_of(a, a_off), constants, tangent)
        ! Only the current level of the storage term depends on the unknown, so
        ! its coefficient is the whole of the derivative.
        g0(1) = bdf_coefficient(constants, 0) * tangent
    end subroutine g0_water_t

    !> Conduction term of the residual: lambda(T) grad T.
    subroutine f1_energy( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f1 &
        ) bind(C, name="ftcms_plexfem_f1_energy")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f1(*)

        real(real64) :: conductivity(dim, dim)
        PetscInt :: d, e
        call eval_conductivity(real(u(u_off(1) + 1), real64), pressure_of(nf, u_off, u, constants), &
                               material_of(a, a_off), constants, conductivity)


        do d = 1, dim
            f1(d) = 0.0_PETSC_REAL_KIND
            do e = 1, dim
                f1(d) = f1(d) + conductivity(d, e) * u_x(u_off_x(1) + e)
            end do
        end do
    end subroutine f1_energy

    !> d(f1)/dT = (dlambda/dT) grad T.
    subroutine g2_energy( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_plexfem_g2_energy")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, value :: u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g2(*)

        real(real64) :: conductivity(dim, dim), dconductivity_dT(dim, dim)
        PetscInt :: d, e

        call eval_conductivity_tangent(real(u(u_off(1) + 1), real64), pressure_of(nf, u_off, u, constants), &
                                       material_of(a, a_off), constants, conductivity, dconductivity_dT)

        do d = 1, dim
            g2(d) = 0.0_PETSC_REAL_KIND
            do e = 1, dim
                g2(d) = g2(d) + dconductivity_dT(d, e) * u_x(u_off_x(1) + e)
            end do
        end do
    end subroutine g2_energy

    !> d(f1_energy)/dp = (dlambda/dp) grad T.
    !>
    !> The energy equation's pressure column. Conduction is a function of the
    !> phase composition and pressure sets it, so this is not a refinement of
    !> the temperature tangent above but the other half of the same derivative.
    !> Leaving it out costs more than leaving out dlambda/dT, because the
    !> increments a Newton step takes in pressure are six orders larger than the
    !> ones it takes in temperature.
    subroutine g2_energy_p( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_plexfem_g2_energy_p")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, value :: u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g2(*)

        real(real64) :: conductivity(dim, dim), dconductivity_dp(dim, dim)
        PetscInt :: d, e

        call eval_conductivity_pressure_tangent(real(u(u_off(1) + 1), real64), &
                                                pressure_of(nf, u_off, u, constants), &
                                                material_of(a, a_off), constants, &
                                                conductivity, dconductivity_dp)

        do d = 1, dim
            g2(d) = 0.0_PETSC_REAL_KIND
            do e = 1, dim
                g2(d) = g2(d) + dconductivity_dp(d, e) * u_x(u_off_x(1) + e)
            end do
        end do
    end subroutine g2_energy_p

    !> d(f1_water)/dp = (dD/dp) grad p + (dK_lT/dp) grad T.
    subroutine g2_water_p( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_plexfem_g2_water_p")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, value :: u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g2(*)

        real(real64) :: ddiffusion_dT(dim, dim), ddiffusion_dp(dim, dim)
        real(real64) :: dcoupling_dT(dim, dim), dcoupling_dp(dim, dim)
        real(real64) :: dgravity_dT, dgravity_dp
        PetscInt :: d, e

        call eval_water_flux_tangent(real(u(u_off(1) + 1), real64), &
                                     pressure_of(nf, u_off, u, constants), &
                                     material_of(a, a_off), constants, &
                                     element_of(a, a_off, A_DARCY), element_of(a, a_off, A_COUPLING), &
                                     ddiffusion_dT, ddiffusion_dp, dcoupling_dT, dcoupling_dp)
        ! The gravity coefficient is the cell's, held fixed over the assembly
        ! like the conductivities, so it contributes no derivative here.
        dgravity_dT = 0.0d0
        dgravity_dp = 0.0d0

        do d = 1, dim
            g2(d) = 0.0_PETSC_REAL_KIND
            do e = 1, dim
                g2(d) = g2(d) + ddiffusion_dp(d, e) * u_x(u_off_x(2) + e) &
                        + dcoupling_dp(d, e) * u_x(u_off_x(1) + e)
            end do
        end do
        g2(dim) = g2(dim) + dgravity_dp
    end subroutine g2_water_p

    !> d(f1_water)/dT = (dD/dT) grad p + (dK_lT/dT) grad T.
    subroutine g2_water_t( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_plexfem_g2_water_t")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, value :: u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g2(*)

        real(real64) :: ddiffusion_dT(dim, dim), ddiffusion_dp(dim, dim)
        real(real64) :: dcoupling_dT(dim, dim), dcoupling_dp(dim, dim)
        real(real64) :: dgravity_dT, dgravity_dp
        PetscInt :: d, e

        call eval_water_flux_tangent(real(u(u_off(1) + 1), real64), &
                                     pressure_of(nf, u_off, u, constants), &
                                     material_of(a, a_off), constants, &
                                     element_of(a, a_off, A_DARCY), element_of(a, a_off, A_COUPLING), &
                                     ddiffusion_dT, ddiffusion_dp, dcoupling_dT, dcoupling_dp)
        ! The gravity coefficient is the cell's, held fixed over the assembly
        ! like the conductivities, so it contributes no derivative here.
        dgravity_dT = 0.0d0
        dgravity_dp = 0.0d0

        do d = 1, dim
            g2(d) = 0.0_PETSC_REAL_KIND
            do e = 1, dim
                g2(d) = g2(d) + ddiffusion_dT(d, e) * u_x(u_off_x(2) + e) &
                        + dcoupling_dT(d, e) * u_x(u_off_x(1) + e)
            end do
        end do
        g2(dim) = g2(dim) + dgravity_dT
    end subroutine g2_water_t

    !> d(f1)/d(grad T) = lambda(T), stored row by row.
    subroutine g3_energy( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g3 &
        ) bind(C, name="ftcms_plexfem_g3_energy")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, value :: u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g3(*)

        real(real64) :: conductivity(dim, dim)
        PetscInt :: d, e

        call eval_conductivity(real(u(u_off(1) + 1), real64), pressure_of(nf, u_off, u, constants), &
                               material_of(a, a_off), constants, conductivity)

        do d = 1, dim
            do e = 1, dim
                g3((d - 1) * dim + e) = conductivity(d, e)
            end do
        end do
    end subroutine g3_energy

    !> Water flux of the residual: D_HH grad p + D_HT grad T.
    !>
    !> The mixed form of the mass balance. The storage term is lumped and lives
    !> in plexfem_run, the same way the energy equation's does; what is left
    !> here is the Darcy flux and the part of it the temperature drives, which
    !> is what carries water towards a freezing front.
    subroutine f1_water( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f1 &
        ) bind(C, name="ftcms_plexfem_f1_water")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f1(*)

        real(real64) :: darcy_element, coupling_element, gravity
        PetscInt :: d

        ! The liquid pressure the flux follows is not the unknown: it is the
        ! effective suction the phase split publishes, and what drives the flux
        ! is the gradient of its nodal interpolant. That interpolant is an
        ! auxiliary field, so PETSc differentiates it and hands the gradient over
        ! here -- which a pointwise chain rule on grad T and grad p cannot
        ! reproduce where the suction turns over across a freezing front.
        !
        ! The thermo-osmotic coefficient acts on grad T separately, as it does in
        ! the older path: it is the temperature dependence of surface tension and
        ! is independent of the potential.
        darcy_element = element_of(a, a_off, A_DARCY)
        coupling_element = element_of(a, a_off, A_COUPLING)
        gravity = element_of(a, a_off, A_GRAVITY)

        do d = 1, dim
            f1(d) = darcy_element * a_x(a_off_x(A_POTENTIAL) + d) &
                    + coupling_element * u_x(u_off_x(1) + d)
        end do
        ! Gravity acts along the last coordinate: this path reads meshes with a
        ! coordinate for every axis, so the vertical is the last of them whether
        ! the mesh is an x-z plane or a full volume, and the model returns zero
        ! for an analysis whose plane gravity has no component in.
        f1(dim) = f1(dim) + gravity
    end subroutine f1_water

    !> d(f1_water)/d(grad p) = D_HH.
    subroutine g3_water_pp( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g3 &
        ) bind(C, name="ftcms_plexfem_g3_water_pp")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, value :: u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g3(*)

        real(real64) :: darcy(dim, dim), coupling(dim, dim)
        PetscInt :: d, e

        call eval_water_driver(real(u(u_off(1) + 1), real64), real(u(u_off(2) + 1), real64), &
                               material_of(a, a_off), constants, &
                               element_of(a, a_off, A_DARCY), element_of(a, a_off, A_COUPLING), &
                               darcy, coupling)

        do d = 1, dim
            do e = 1, dim
                g3((d - 1) * dim + e) = darcy(d, e)
            end do
        end do
    end subroutine g3_water_pp

    !> d(f1_water)/d(grad T) = D_HT.
    subroutine g3_water_pt( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g3 &
        ) bind(C, name="ftcms_plexfem_g3_water_pt")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, value :: u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g3(*)

        real(real64) :: darcy(dim, dim), coupling(dim, dim)
        PetscInt :: d, e

        call eval_water_driver(real(u(u_off(1) + 1), real64), real(u(u_off(2) + 1), real64), &
                               material_of(a, a_off), constants, &
                               element_of(a, a_off, A_DARCY), element_of(a, a_off, A_COUPLING), &
                               darcy, coupling)

        do d = 1, dim
            do e = 1, dim
                g3((d - 1) * dim + e) = coupling(d, e)
            end do
        end do
    end subroutine g3_water_pt

    !> Boundary term of the residual: the heat leaving through the surface.
    !>
    !> Integrating the conduction term by parts leaves -\oint v lambda grad T . n
    !> in the residual, and with the heat flux q = -lambda grad T that is
    !> +\oint v q.n, the flux out of the domain. The project states a boundary
    !> flux as an influx, the way strategy_flux and strategy_robin do, so the
    !> sign is turned round here rather than in the input layer.
    !>
    !> Which condition this boundary carries, and its constants, arrive in the
    !> auxiliary fields: PETSc looks up a boundary's own auxiliary vector by
    !> the label value it was registered with, so one pair of kernels serves
    !> every natural boundary of the mesh.
    subroutine bd_f0_energy( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, normal, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_plexfem_bd_f0_energy")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscReal, intent(in) :: normal(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f0(*)

        integer(int32) :: kind
        real(real64) :: temperature, value_1, value_2

        ! Nothing while the nodal volumes are being measured: what is wanted
        ! then is the volume alone.
        if (real(constants(C_MEASURING_MATERIAL), real64) /= 0.0d0) then
            f0(1) = 0.0_PETSC_REAL_KIND
            return
        end if

        kind = int(nint(real(a(a_off(A_BC_KIND) + 1), real64)), int32)
        value_1 = real(a(a_off(A_BC_VALUE_1) + 1), real64)
        value_2 = real(a(a_off(A_BC_VALUE_2) + 1), real64)
        temperature = real(u(u_off(1) + 1), real64)

        select case (kind)
        case (BC_FLUX)
            ! The project gives the influx; the residual carries the outflux.
            f0(1) = -value_1
        case (BC_ROBIN)
            ! Newton cooling: h (T - T_ambient) leaves the surface.
            f0(1) = value_1 * (temperature - value_2)
        case default
            f0(1) = 0.0_PETSC_REAL_KIND
        end select
    end subroutine bd_f0_energy

    !> d(bd_f0)/dT. Only a Robin condition depends on the temperature.
    subroutine bd_g0_energy( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, normal, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_plexfem_bd_g0_energy")
        implicit none
        PetscInt, value :: dim
        PetscInt, value :: nf
        PetscInt, value :: nf_aux
        PetscInt, intent(in) :: u_off(*)
        PetscInt, intent(in) :: u_off_x(*)
        PetscScalar, intent(in) :: u(*)
        PetscScalar, intent(in) :: u_t(*)
        PetscScalar, intent(in) :: u_x(*)
        PetscInt, intent(in) :: a_off(*)
        PetscInt, intent(in) :: a_off_x(*)
        PetscScalar, intent(in) :: a(*)
        PetscScalar, intent(in) :: a_t(*)
        PetscScalar, intent(in) :: a_x(*)
        PetscReal, value :: time
        PetscReal, value :: u_t_shift
        PetscReal, intent(in) :: x(*)
        PetscReal, intent(in) :: normal(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: g0(*)

        integer(int32) :: kind

        if (real(constants(C_MEASURING_MATERIAL), real64) /= 0.0d0) then
            g0(1) = 0.0_PETSC_REAL_KIND
            return
        end if

        kind = int(nint(real(a(a_off(A_BC_KIND) + 1), real64)), int32)

        if (kind == BC_ROBIN) then
            g0(1) = real(a(a_off(A_BC_VALUE_1) + 1), real64)
        else
            g0(1) = 0.0_PETSC_REAL_KIND
        end if
    end subroutine bd_g0_energy

    !> Temperature prescribed on an essential boundary.
    !>
    !> PETSc passes back the context the boundary was registered with, which
    !> here is the boundary's own id, so one callback serves every essential
    !> boundary of this field.
    function bc_essential_temperature(dim, time, x, num_components, u, ctx) &
        result(ierr) bind(C, name="ftcms_plexfem_bc_temperature")
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        integer(c_int), pointer :: boundary_id

        call c_f_pointer(ctx, boundary_id)
        u(1) = prescribed_value(0_c_int, boundary_id)
        ierr = 0
    end function bc_essential_temperature

    !> Pore pressure prescribed on an essential boundary.
    function bc_essential_pressure(dim, time, x, num_components, u, ctx) &
        result(ierr) bind(C, name="ftcms_plexfem_bc_pressure")
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        integer(c_int), pointer :: boundary_id

        call c_f_pointer(ctx, boundary_id)
        u(1) = prescribed_value(1_c_int, boundary_id)
        ierr = 0
    end function bc_essential_pressure

end module app_plexfem_pointwise
