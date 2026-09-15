!> Fortran access to a PetscFE's basis tabulation.
!>
!> PetscFE can tabulate its basis at any reference point, but the generated
!> Fortran binding hands back a descriptor it never builds, so the values
!> cannot be read from Fortran. The C side of this module does the call and
!> copies the result into arrays the caller owns.
module core_interop_petscfe
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
    implicit none
    private

    public :: petscfe_tabulate
    public :: petscfe_field_size
    public :: petscfe_evaluate_cell
    public :: petscfe_reference_coordinates

    interface
        function ftcms_petscfe_tabulate(fe, num_points, points, values, gradients) &
            bind(c, name="ftcms_petscfe_tabulate") result(ierr)
            import :: c_int, c_double, c_ptr
            implicit none
            type(c_ptr), value :: fe
            integer(c_int), value :: num_points
            real(c_double), intent(in) :: points(*)
            real(c_double), intent(inout) :: values(*)
            real(c_double), intent(inout) :: gradients(*)
            integer(c_int) :: ierr
        end function ftcms_petscfe_tabulate

        function ftcms_petscfe_field_size(dm, field, num_basis, num_components) &
            bind(c, name="ftcms_petscfe_field_size") result(ierr)
            import :: c_int, c_ptr
            implicit none
            type(c_ptr), value :: dm
            integer(c_int), value :: field
            integer(c_int), intent(inout) :: num_basis
            integer(c_int), intent(inout) :: num_components
            integer(c_int) :: ierr
        end function ftcms_petscfe_field_size

        function ftcms_petscfe_evaluate_cell(dm, local_x, field, cell, num_points, ref_points, &
                                             values, divergence) &
            bind(c, name="ftcms_petscfe_evaluate_cell") result(ierr)
            import :: c_int, c_double, c_ptr
            implicit none
            type(c_ptr), value :: dm
            type(c_ptr), value :: local_x
            integer(c_int), value :: field
            integer(c_int), value :: cell
            integer(c_int), value :: num_points
            real(c_double), intent(in) :: ref_points(*)
            real(c_double), intent(inout) :: values(*)
            real(c_double), intent(inout) :: divergence(*)
            integer(c_int) :: ierr
        end function ftcms_petscfe_evaluate_cell

        function ftcms_petscfe_reference_coordinates(dm, cell, num_points, points, ref_points) &
            bind(c, name="ftcms_petscfe_reference_coordinates") result(ierr)
            import :: c_int, c_double, c_ptr
            implicit none
            type(c_ptr), value :: dm
            integer(c_int), value :: cell
            integer(c_int), value :: num_points
            real(c_double), intent(in) :: points(*)
            real(c_double), intent(inout) :: ref_points(*)
            integer(c_int) :: ierr
        end function ftcms_petscfe_reference_coordinates
    end interface

contains

    !> Basis values and reference gradients at the given reference points.
    !>
    !> points is (dim, num_points); values is (num_basis, num_points);
    !> gradients is (dim, num_basis, num_points).
    subroutine petscfe_tabulate(fe, points, values, gradients, ierr)
        implicit none
        type(c_ptr), intent(in) :: fe
        ! The C side reads and writes these in place, so a copy made on the way
        ! in would be both wasted work and, for the results, wrong.
        real(c_double), intent(in), contiguous :: points(:, :)
        real(c_double), intent(inout), contiguous :: values(:, :)
        real(c_double), intent(inout), contiguous :: gradients(:, :, :)
        integer(c_int), intent(inout) :: ierr

        ierr = ftcms_petscfe_tabulate(fe, int(size(points, 2), c_int), points, values, gradients)
    end subroutine petscfe_tabulate

    !> Number of basis functions and components of one field of a DM.
    !>
    !> One degree of freedom per basis function, so num_basis is also the
    !> number of entries this field contributes to a cell closure.
    subroutine petscfe_field_size(dm, field, num_basis, num_components, ierr)
        implicit none
        type(c_ptr), intent(in) :: dm
        integer(c_int), intent(in) :: field
        integer(c_int), intent(inout) :: num_basis
        integer(c_int), intent(inout) :: num_components
        integer(c_int), intent(inout) :: ierr

        ierr = ftcms_petscfe_field_size(dm, field, num_basis, num_components)
    end subroutine petscfe_field_size

    !> Physical value and divergence of one field at points of one cell.
    !>
    !> ref_points is (dim, num_points) on the reference cell; values is
    !> (num_components, num_points); divergence is (num_points). The mapping
    !> from the reference cell is PETSc's own pushforward, which is the
    !> contravariant Piola map for an H(div) field.
    subroutine petscfe_evaluate_cell(dm, local_x, field, cell, ref_points, values, divergence, ierr)
        implicit none
        type(c_ptr), intent(in) :: dm
        type(c_ptr), intent(in) :: local_x
        integer(c_int), intent(in) :: field
        integer(c_int), intent(in) :: cell
        real(c_double), intent(in), contiguous :: ref_points(:, :)
        real(c_double), intent(inout), contiguous :: values(:, :)
        real(c_double), intent(inout), contiguous :: divergence(:)
        integer(c_int), intent(inout) :: ierr

        ierr = ftcms_petscfe_evaluate_cell(dm, local_x, field, cell, int(size(ref_points, 2), c_int), &
                                           ref_points, values, divergence)
    end subroutine petscfe_evaluate_cell

    !> Reference coordinates of physical points inside one affine cell.
    !>
    !> points and ref_points are both (dim, num_points).
    subroutine petscfe_reference_coordinates(dm, cell, points, ref_points, ierr)
        implicit none
        type(c_ptr), intent(in) :: dm
        integer(c_int), intent(in) :: cell
        real(c_double), intent(in), contiguous :: points(:, :)
        real(c_double), intent(inout), contiguous :: ref_points(:, :)
        integer(c_int), intent(inout) :: ierr

        ierr = ftcms_petscfe_reference_coordinates(dm, cell, int(size(points, 2), c_int), points, ref_points)
    end subroutine petscfe_reference_coordinates

end module core_interop_petscfe
