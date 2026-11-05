!> Darcy's law in the intrinsic plane: the dimension every kernel works in, and
!> the conversion between a conductivity tensor and the inverse mobility the
!> flux row is written with.
module app_ftcms_darcy
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none
    private

    public :: PLANE_DIMENSION
    public :: invert_mobility
    public :: resistance_of

    !> The kernels solve in the two intrinsic coordinates of the analysis plane.
    integer(int32), parameter :: PLANE_DIMENSION = 2

contains

    !> \(\Lambda^{-1}=\rho g\,K^{-1}\), row major, from a saturated conductivity \(K\) in m/s.
    !>
    !> Fails by returning `ok = .false.` on a nonpositive density or gravity, or
    !> a tensor that is not symmetric positive definite.
    subroutine invert_mobility(conductivity, density, gravity_magnitude, inverse_mobility, ok)
        implicit none
        real(real64), intent(in) :: conductivity(PLANE_DIMENSION, PLANE_DIMENSION)
        real(real64), intent(in) :: density
        real(real64), intent(in) :: gravity_magnitude
        real(real64), intent(inout) :: inverse_mobility(PLANE_DIMENSION*PLANE_DIMENSION)
        logical, intent(inout) :: ok

        real(real64) :: determinant, weight, asymmetry, scale

        ok = .false.

        if (density <= 0.0d0 .or. gravity_magnitude <= 0.0d0) return

        scale = maxval(abs(conductivity))
        if (scale <= 0.0d0) return

        asymmetry = abs(conductivity(1, 2) - conductivity(2, 1))
        if (asymmetry > 1.0d-12*scale) return

        determinant = conductivity(1, 1)*conductivity(2, 2) - conductivity(1, 2)*conductivity(2, 1)
        if (conductivity(1, 1) <= 0.0d0 .or. conductivity(2, 2) <= 0.0d0) return
        if (determinant <= 0.0d0) return

        weight = density*gravity_magnitude/determinant
        inverse_mobility(1) = weight*conductivity(2, 2)
        inverse_mobility(2) = -weight*conductivity(1, 2)
        inverse_mobility(3) = -weight*conductivity(2, 1)
        inverse_mobility(4) = weight*conductivity(1, 1)
        ok = .true.
    end subroutine invert_mobility

    !> \(\Lambda^{-1}\mathbf q\), with \(\Lambda^{-1}\) row major.
    pure function resistance_of(inverse_mobility, flux) result(resistance)
        implicit none
        real(real64), intent(in) :: inverse_mobility(PLANE_DIMENSION*PLANE_DIMENSION)
        real(real64), intent(in) :: flux(PLANE_DIMENSION)
        real(real64) :: resistance(PLANE_DIMENSION)

        integer(int32) :: i, j

        do i = 1, PLANE_DIMENSION
            resistance(i) = 0.0d0
            do j = 1, PLANE_DIMENSION
                resistance(i) = resistance(i) + inverse_mobility((i - 1)*PLANE_DIMENSION + j)*flux(j)
            end do
        end do
    end function resistance_of

end module app_ftcms_darcy
