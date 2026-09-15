#include <petsc/finclude/petscsys.h>

!> The mixed Darcy and conservation weak form, as PetscDS pointwise functions.
!>
!> PETSc evaluates these at every quadrature point and owns the element loop and
!> the assembly. Each reaches PETSc through a C function pointer, so the
!> argument lists match `PetscPointFn`, `PetscPointJacFn`, `PetscBdPointFn` and
!> `PetscPointExactSolutionFn` in petscdstypes.h exactly.
!>
!> With flux test function \(\mathbf v\), pressure test function \(w\) and
!> \(\chi=p_l\), the residual is
!>
!> \[
!>   R_q(\mathbf v)=(\mathbf\Lambda^{-1}\mathbf q_l-\rho_l\mathbf g,\mathbf v)
!>   -(\chi,\nabla\cdot\mathbf v)+\langle\chi_D,\mathbf v\cdot\mathbf n\rangle_{\Gamma_D},
!>   \qquad
!>   R_p(w)=(D_tW+\nabla\cdot\mathbf q_l-s,w).
!> \]
!>
!> PETSc integrates a test-value term `f0` and a test-gradient term `f1`, so the
!> pressure enters the flux equation as `f1 = -chi I`: its element gradient term
!> is \(-\int\chi\nabla\cdot\mathbf v\). The transport must not instead be
!> written as `f1` of the pressure equation -- the P0 test functions have a
!> vanishing element gradient and the whole term would silently be zero.
!>
!> The divergence read out of `u_x` is the physical one. PETSc has already
!> applied the contravariant Piola map to the RT basis gradients by the time a
!> kernel sees them, so no Jacobian factor appears here.
!>
!> The unknowns are dimensionless: `u` carries flux in units of the flux scale
!> and pressure in units of the pressure scale, both run constants. The flux row
!> is divided by the pressure scale and the conservation row by flux scale over
!> length scale, which is what lets a mobility spanning seven orders of magnitude
!> and a pressure of order \(10^4\) Pa share one matrix.
!>
!> Time is discretised by the caller, so `u_t` is never read and PETSc's
!> `u_t_shift` is always zero.
module app_mixed_hydraulic_pointwise
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_loc, c_f_pointer
    use :: petscsys
    use :: app_mixed_hydraulic_physics, only:MIXED_DIMENSION, &
                                             C_DT, C_FLUX_SCALE, C_PRESSURE_SCALE, C_LENGTH_SCALE, &
                                             C_DENSITY, C_GRAVITY_X, C_GRAVITY_Z, &
                                             A_INVERSE_MOBILITY, A_STORAGE, A_SOURCE, &
                                             A_PRESSURE_OLD, A_BC_PRESSURE
    implicit none
    private

    public :: f0_flux, f1_flux, f0_conservation
    public :: g0_flux_flux, g2_flux_pressure, g1_pressure_flux, g0_pressure_pressure
    public :: bd_f0_flux_pressure
    public :: bc_essential_flux
    public :: set_flux_boundary_value, flux_boundary_context, reset_flux_boundaries
    public :: MAX_FLUX_BOUNDARIES

    !> How many prescribed-flux boundaries one run may register.
    integer(int32), parameter :: MAX_FLUX_BOUNDARIES = 16

    !> The prescribed flux of one boundary, reachable from the essential
    !> callback as its context.
    !>
    !> Interoperable because PETSc hands the pointer back to a `bind(C)`
    !> procedure; the flux is physical and the scale converts it to the
    !> dimensionless unknown the dual functional is applied to.
    type, bind(C) :: type_flux_boundary_value
        real(c_double) :: flux(MIXED_DIMENSION)
        real(c_double) :: flux_scale
    end type type_flux_boundary_value

    type(type_flux_boundary_value), save, target :: flux_boundary_values(MAX_FLUX_BOUNDARIES)

contains

    !> Forget every prescribed flux of a previous run.
    subroutine reset_flux_boundaries()
        implicit none

        integer(int32) :: i

        do i = 1, MAX_FLUX_BOUNDARIES
            flux_boundary_values(i)%flux = 0.0_c_double
            flux_boundary_values(i)%flux_scale = 1.0_c_double
        end do
    end subroutine reset_flux_boundaries

    !> Record the flux vector one registered boundary prescribes.
    subroutine set_flux_boundary_value(slot, flux, flux_scale)
        implicit none
        !> Index in 1..MAX_FLUX_BOUNDARIES.
        integer(int32), intent(in) :: slot
        !> Prescribed flux vector in m/s; only its normal component is used.
        real(real64), intent(in) :: flux(MIXED_DIMENSION)
        !> The flux unit of the unknown, m/s, positive.
        real(real64), intent(in) :: flux_scale

        flux_boundary_values(slot)%flux = real(flux, c_double)
        flux_boundary_values(slot)%flux_scale = real(flux_scale, c_double)
    end subroutine set_flux_boundary_value

    !> The context PETSc hands back to the essential callback of one boundary.
    function flux_boundary_context(slot) result(ctx)
        implicit none
        integer(int32), intent(in) :: slot
        type(c_ptr) :: ctx

        ctx = c_loc(flux_boundary_values(slot))
    end function flux_boundary_context

    !> The cell's inverse mobility, row-major, from its auxiliary field.
    pure function inverse_mobility_of(a, a_off) result(inverse_mobility)
        implicit none
        PetscScalar, intent(in) :: a(*)
        PetscInt, intent(in) :: a_off(*)
        real(real64) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION)

        integer(int32) :: k

        do k = 1, MIXED_DIMENSION*MIXED_DIMENSION
            inverse_mobility(k) = real(a(a_off(A_INVERSE_MOBILITY) + k), real64)
        end do
    end function inverse_mobility_of

    !> Darcy's law at a quadrature point,
    !> \(f_0=(\mathbf\Lambda^{-1}\mathbf q_l-\rho_l\mathbf g)/p_{\mathrm{scale}}\).
    subroutine f0_flux( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_mixed_f0_flux")
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

        real(real64) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION)
        real(real64) :: flux(MIXED_DIMENSION), gravity(MIXED_DIMENSION)
        real(real64) :: flux_scale, pressure_scale, density, resistance
        integer(int32) :: i, j

        flux_scale = real(constants(C_FLUX_SCALE), real64)
        pressure_scale = real(constants(C_PRESSURE_SCALE), real64)
        density = real(constants(C_DENSITY), real64)
        gravity(1) = real(constants(C_GRAVITY_X), real64)
        gravity(2) = real(constants(C_GRAVITY_Z), real64)

        inverse_mobility = inverse_mobility_of(a, a_off)
        do i = 1, MIXED_DIMENSION
            flux(i) = flux_scale*real(u(u_off(1) + i), real64)
        end do

        do i = 1, MIXED_DIMENSION
            resistance = 0.0_real64
            do j = 1, MIXED_DIMENSION
                resistance = resistance + inverse_mobility((i - 1)*MIXED_DIMENSION + j)*flux(j)
            end do
            f0(i) = real((resistance - density*gravity(i))/pressure_scale, PETSC_REAL_KIND)
        end do
    end subroutine f0_flux

    !> The pressure term of Darcy's law, \(f_1=-\chi\mathbf I\).
    !>
    !> Integrated against the flux test gradient, this is
    !> \(-\int\chi\nabla\cdot\mathbf v\), which is where the pressure and the
    !> conservation equation meet.
    subroutine f1_flux( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f1 &
        ) bind(C, name="ftcms_mixed_f1_flux")
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

        real(real64) :: pressure
        integer(int32) :: c, d

        ! Dimensionless already: the row is divided by the pressure scale and
        ! the unknown is measured in it.
        pressure = real(u(u_off(2) + 1), real64)

        do c = 1, MIXED_DIMENSION
            do d = 1, MIXED_DIMENSION
                if (c == d) then
                    f1((c - 1)*MIXED_DIMENSION + d) = real(-pressure, PETSC_REAL_KIND)
                else
                    f1((c - 1)*MIXED_DIMENSION + d) = 0.0_PETSC_REAL_KIND
                end if
            end do
        end do
    end subroutine f1_flux

    !> Water conservation, \(f_0=(D_tW+\nabla\cdot\mathbf q_l-s)L/q_{\mathrm{scale}}\).
    !>
    !> The storage rate differences the stored quantity itself rather than
    !> multiplying a capacity by a pressure increment; with the linear law of
    !> this stage the two agree, and the difference form is what the freezing
    !> closure will need.
    subroutine f0_conservation( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_mixed_f0_conservation")
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

        real(real64) :: flux_scale, pressure_scale, length_scale, dt
        real(real64) :: divergence, storage, pressure_old, source, storage_rate
        integer(int32) :: c

        flux_scale = real(constants(C_FLUX_SCALE), real64)
        pressure_scale = real(constants(C_PRESSURE_SCALE), real64)
        length_scale = real(constants(C_LENGTH_SCALE), real64)
        dt = real(constants(C_DT), real64)

        ! PETSc stores the gradient component-major: du_c/dx_d at c*dim + d.
        divergence = 0.0_real64
        do c = 1, MIXED_DIMENSION
            divergence = divergence + real(u_x(u_off_x(1) + (c - 1)*MIXED_DIMENSION + c), real64)
        end do

        storage = real(a(a_off(A_STORAGE) + 1), real64)
        pressure_old = real(a(a_off(A_PRESSURE_OLD) + 1), real64)
        source = real(a(a_off(A_SOURCE) + 1), real64)

        storage_rate = 0.0_real64
        if (dt > 0.0_real64) then
            storage_rate = storage*(pressure_scale*real(u(u_off(2) + 1), real64) - pressure_old)/dt
        end if

        f0(1) = real((storage_rate + flux_scale*divergence - source)*length_scale/flux_scale, PETSC_REAL_KIND)
    end subroutine f0_conservation

    !> \(\partial f_0^{q}/\partial\mathbf q=\mathbf\Lambda^{-1}q_{\mathrm{scale}}/p_{\mathrm{scale}}\).
    subroutine g0_flux_flux( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_mixed_g0_flux_flux")
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

        real(real64) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION)
        real(real64) :: weight
        integer(int32) :: k

        weight = real(constants(C_FLUX_SCALE), real64)/real(constants(C_PRESSURE_SCALE), real64)
        inverse_mobility = inverse_mobility_of(a, a_off)

        ! g0 is indexed test component major: g0[fc*Nc_trial + gc].
        do k = 1, MIXED_DIMENSION*MIXED_DIMENSION
            g0(k) = real(inverse_mobility(k)*weight, PETSC_REAL_KIND)
        end do
    end subroutine g0_flux_flux

    !> \(\partial f_1^{q}/\partial\chi=-\mathbf I\).
    subroutine g2_flux_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g2 &
        ) bind(C, name="ftcms_mixed_g2_flux_pressure")
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

        integer(int32) :: c, d

        ! One trial component, so g2[(fc*1 + 0)*dim + df] with df the test
        ! gradient direction: the diagonal of -chi I.
        do c = 1, MIXED_DIMENSION
            do d = 1, MIXED_DIMENSION
                if (c == d) then
                    g2((c - 1)*MIXED_DIMENSION + d) = -1.0_PETSC_REAL_KIND
                else
                    g2((c - 1)*MIXED_DIMENSION + d) = 0.0_PETSC_REAL_KIND
                end if
            end do
        end do
    end subroutine g2_flux_pressure

    !> \(\partial f_0^{p}/\partial(\nabla\mathbf q)=L\,\mathbf I\).
    subroutine g1_pressure_flux( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g1 &
        ) bind(C, name="ftcms_mixed_g1_pressure_flux")
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
        PetscScalar, intent(inout) :: g1(*)

        real(real64) :: length_scale
        integer(int32) :: c, d

        length_scale = real(constants(C_LENGTH_SCALE), real64)

        ! One test component, so g1[(0*Nc_trial + gc)*dim + dg]: the trace that
        ! picks the divergence out of the trial gradient.
        do c = 1, MIXED_DIMENSION
            do d = 1, MIXED_DIMENSION
                if (c == d) then
                    g1((c - 1)*MIXED_DIMENSION + d) = real(length_scale, PETSC_REAL_KIND)
                else
                    g1((c - 1)*MIXED_DIMENSION + d) = 0.0_PETSC_REAL_KIND
                end if
            end do
        end do
    end subroutine g1_pressure_flux

    !> \(\partial f_0^{p}/\partial\chi=(\partial W/\partial p)\,p_{\mathrm{scale}}L/(\Delta t\,q_{\mathrm{scale}})\).
    subroutine g0_pressure_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, u_t_shift, x, &
        num_constants, constants, g0 &
        ) bind(C, name="ftcms_mixed_g0_pressure_pressure")
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

        real(real64) :: dt, storage

        dt = real(constants(C_DT), real64)
        storage = real(a(a_off(A_STORAGE) + 1), real64)

        if (dt > 0.0_real64) then
            g0(1) = real(storage*real(constants(C_PRESSURE_SCALE), real64) &
                         *real(constants(C_LENGTH_SCALE), real64) &
                         /(dt*real(constants(C_FLUX_SCALE), real64)), PETSC_REAL_KIND)
        else
            g0(1) = 0.0_PETSC_REAL_KIND
        end if
    end subroutine g0_pressure_pressure

    !> The prescribed-pressure boundary term \(\langle\chi_D,\mathbf v\cdot\mathbf n\rangle\).
    !>
    !> Written as \(f_0=\chi_D\mathbf n/p_{\mathrm{scale}}\), which PETSc
    !> integrates against the flux test function. It carries no unknown, so the
    !> boundary contributes nothing to the Jacobian.
    subroutine bd_f0_flux_pressure( &
        dim, nf, nf_aux, u_off, u_off_x, u, u_t, u_x, &
        a_off, a_off_x, a, a_t, a_x, time, x, n, &
        num_constants, constants, f0 &
        ) bind(C, name="ftcms_mixed_bd_f0_flux_pressure")
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
        PetscReal, intent(in) :: n(*)
        PetscInt, value :: num_constants
        PetscScalar, intent(in) :: constants(*)
        PetscScalar, intent(inout) :: f0(*)

        real(real64) :: boundary_pressure, pressure_scale
        integer(int32) :: i

        pressure_scale = real(constants(C_PRESSURE_SCALE), real64)
        boundary_pressure = real(a(a_off(A_BC_PRESSURE) + 1), real64)

        do i = 1, MIXED_DIMENSION
            f0(i) = real(boundary_pressure*real(n(i), real64)/pressure_scale, PETSC_REAL_KIND)
        end do
    end subroutine bd_f0_flux_pressure

    !> The flux a prescribed-flux boundary constrains, in units of the flux scale.
    !>
    !> PETSc applies the RT facet functional to whatever this returns, so only
    !> the normal component of the given vector reaches the constrained degree
    !> of freedom; the tangential part is discarded by the functional itself.
    function bc_essential_flux(dim, time, x, num_components, u, ctx) &
        bind(C, name="ftcms_mixed_bc_essential_flux") result(ierr)
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        type(type_flux_boundary_value), pointer :: boundary
        integer(int32) :: i

        call c_f_pointer(ctx, boundary)
        do i = 1, MIXED_DIMENSION
            u(i) = real(boundary%flux(i)/boundary%flux_scale, PETSC_REAL_KIND)
        end do
        ierr = 0
    end function bc_essential_flux

end module app_mixed_hydraulic_pointwise
