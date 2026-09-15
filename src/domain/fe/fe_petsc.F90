#include <petsc/finclude/petscdm.h>

!> A finite element whose basis comes from PETSc.
!>
!> The shape functions are tabulated by PetscFE and the connectivity comes from
!> the DMPlex closure, and PetscFE numbers its basis in that same closure order.
!> There is therefore no node-ordering convention to reconcile: basis i belongs
!> to node i for every element type and order, simplices included.
module domain_fe_petsc
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double
    use :: petscdm
    use :: module_core
    use :: core_interop_petscfe, only: petscfe_tabulate
    use :: domain_base_fe, only: abst_fe
    use :: domain_fe_integration, only: type_gauss_integration_rule
    implicit none
    private

    public :: type_fe_petsc
    public :: construct_fe_petsc

    type, extends(abst_fe) :: type_fe_petsc
        private
        PetscFE :: fe
        type(c_ptr) :: handle
        logical :: ready = .false.
        integer(int32) :: reference_dim = 0
    contains
        procedure, pass(self) :: calc_measure => calc_measure_fe_petsc
        procedure, pass(self) :: is_inside => is_inside_fe_petsc
        procedure, pass(self) :: tabulate => tabulate_fe_petsc
    end type type_fe_petsc

contains

    !> Build the element PETSc would use for this shape, order and dimension.
    function construct_fe_petsc(fe_type, dimension, order, is_simplex, &
                                integration_order) result(element)
        implicit none
        integer(int32), intent(in) :: fe_type, dimension, order
        logical, intent(in) :: is_simplex
        integer(int32), intent(in) :: integration_order
        class(abst_fe), allocatable :: element

        PetscErrorCode :: ierr
        PetscBool :: simplex
        PetscBool :: ready
        PetscInt :: num_basis, reference_dim
        integer(int32) :: num_nodes

        ! The basis is PETSc's, so an element can be built before any mesh has
        ! been read. Whoever asks first initialises PETSc, the same way the
        ! mesh reader does.
        ierr = 0
        call PetscInitialized(ready, ierr)
        if (ready .eqv. PETSC_FALSE) call PetscInitialize(ierr)
        if (ierr /= 0) error stop "PETSc could not be initialised for the element basis"

        allocate (type_fe_petsc :: element)
        select type (element)
        type is (type_fe_petsc)
            simplex = PETSC_FALSE
            if (is_simplex) simplex = PETSC_TRUE

            call PetscFECreateLagrange(PETSC_COMM_SELF, int(dimension, PETSC_INT_KIND), &
                                       1_PETSC_INT_KIND, simplex, int(order, PETSC_INT_KIND), &
                                       int(max(integration_order, order), PETSC_INT_KIND), &
                                       element%fe, ierr)
            if (ierr /= 0) error stop "PETSc could not build the element basis"

            element%handle = transfer(element%fe%v, element%handle)

            ! The basis PETSc built decides how many nodes the element has and
            ! which reference cell it lives on; nothing here counts them again.
            call PetscFEGetDimension(element%fe, num_basis, ierr)
            if (ierr /= 0) error stop "PETSc element basis has no size"
            call PetscFEGetSpatialDimension(element%fe, reference_dim, ierr)
            if (ierr /= 0) error stop "PETSc element basis has no reference cell"
            num_nodes = int(num_basis, int32)
            element%reference_dim = int(reference_dim, int32)
            element%ready = .true.

            call element%initialize(type=fe_type, dimension=dimension, order=order, &
                                    num_nodes=num_nodes, integration_order=integration_order)

            ! The rule the basis is integrated on comes from the same PetscFE,
            ! so points and basis always live on one reference cell.
            block
                type(type_gauss_integration_rule), pointer :: rule
                PetscQuadrature :: quadrature

                nullify (rule)
                call element%get_integration_rule(rule)
                call PetscFEGetQuadrature(element%fe, quadrature, ierr)
                if (ierr /= 0) error stop "PETSc element has no quadrature"
                if (associated(rule)) call rule%adopt_quadrature(quadrature)
            end block
        end select
    end function construct_fe_petsc

    !> Basis values and reference gradients of every node at one reference point.
    subroutine tabulate_fe_petsc(self, r, psi, dpsi_ref)
        implicit none
        class(type_fe_petsc), intent(in) :: self
        type(type_coordinate_dp), intent(in) :: r
        real(real64), intent(inout), optional :: psi(:)
        real(real64), intent(inout), optional :: dpsi_ref(:, :)

        ! Sized to the reference cell rather than to three dimensions: passing
        ! a section of a larger array to the C binding makes the compiler copy
        ! it in and out at every quadrature point.
        real(c_double) :: point(self%reference_dim, 1)
        real(c_double), allocatable :: values(:, :), gradients(:, :, :)
        integer(c_int) :: cerr
        integer(int32) :: i, k, n

        call self%get_num_nodes(n)
        point(1, 1) = r%x
        if (self%reference_dim >= 2) point(2, 1) = r%y
        if (self%reference_dim >= 3) point(3, 1) = r%z

        allocate (values(n, 1))
        allocate (gradients(self%reference_dim, n, 1))
        values = 0.0d0
        gradients = 0.0d0
        cerr = 0
        call petscfe_tabulate(self%handle, point, values, gradients, cerr)
        if (cerr /= 0) error stop "PETSc could not tabulate the element basis"

        if (present(psi)) then
            do i = 1, min(n, size(psi))
                psi(i) = values(i, 1)
            end do
        end if
        if (present(dpsi_ref)) then
            dpsi_ref = 0.0d0
            do i = 1, min(n, size(dpsi_ref, 2))
                do k = 1, min(self%reference_dim, size(dpsi_ref, 1))
                    dpsi_ref(k, i) = gradients(k, i, 1)
                end do
            end do
        end if

        deallocate (values)
        deallocate (gradients)
    end subroutine tabulate_fe_petsc

    subroutine calc_measure_fe_petsc(self, node_coords, measure)
        implicit none
        class(type_fe_petsc), intent(in) :: self
        real(real64), intent(in) :: node_coords(:, :)
        real(real64), intent(inout) :: measure

        integer(int32) :: p, num_gauss
        real(real64) :: det_j
        type(type_coordinate_dp), pointer, contiguous :: gauss(:)
        real(real64), pointer, contiguous :: weight(:)

        measure = 0.0d0
        call self%get_num_gauss(num_gauss)
        nullify (gauss)
        nullify (weight)
        call self%get_gauss(gauss)
        call self%get_weight(weight)
        if (.not. associated(gauss) .or. .not. associated(weight)) return
        do p = 1, num_gauss
            call self%calc_shape_function(gauss(p), node_coords, determinant_jacobian=det_j)
            measure = measure + weight(p) * abs(det_j)
        end do
    end subroutine calc_measure_fe_petsc

    !> Whether a physical point lies in this element, and where in the
    !> reference cell it lands.
    !>
    !> The isoparametric map is inverted by Newton from the reference centre,
    !> using the same basis and reference gradients PetscFE tabulates. For an
    !> affine element one step is exact; a bilinear or trilinear one converges
    !> in a few.
    !>
    !> Containment is then decided by the basis itself: a first-order Lagrange
    !> basis is non-negative exactly on its own reference cell, on a simplex as
    !> the barycentric coordinates and on a tensor cell as the products of
    !> (1 +/- xi)/2. Testing the basis rather than the coordinates keeps this
    !> free of PETSc's reference-cell convention, whose simplex is not the unit
    !> one. For a higher-order basis the criterion is conservative, since such
    !> a basis also turns negative inside its own cell.
    subroutine is_inside_fe_petsc(self, cartesian, normalized, node_coords, is_in)
        implicit none
        class(type_fe_petsc), intent(in) :: self
        type(type_coordinate_dp), intent(in) :: cartesian
        type(type_coordinate_dp), intent(inout) :: normalized
        real(real64), intent(in) :: node_coords(:, :)
        logical, intent(inout) :: is_in

        ! A point on a face belongs to both of its elements, so the tolerance
        ! only has to keep rounding from excluding it from either.
        real(real64), parameter :: BASIS_TOLERANCE = 1.0d-10
        real(real64), parameter :: STEP_TOLERANCE = 1.0d-12
        integer(int32), parameter :: MAX_NEWTON_STEPS = 20

        integer(int32) :: dim, space_dim, n, i, j, k, step
        real(real64) :: xi(3), target(3), residual(3), correction(3)
        real(real64), allocatable :: psi(:), reference_gradient(:, :)
        real(real64) :: jacobian(3, 3)
        logical :: solved, solved_step

        is_in = .false.
        call normalized%reset()

        call self%get_dimension(dim)
        space_dim = size(node_coords, 1)
        call self%get_num_nodes(n)
        if (dim < 1 .or. dim > 3 .or. n < 1) return

        ! A side of a plane problem has no square Jacobian, so the map cannot be
        ! inverted this way and the point is not located on it.
        if (space_dim /= dim) return

        allocate (psi(n))
        allocate (reference_gradient(dim, n))

        target(:) = 0.0d0
        target(1) = cartesian%x
        if (dim >= 2) target(2) = cartesian%y
        if (dim >= 3) target(3) = cartesian%z

        ! The centre of every reference cell PETSc uses, biunit or simplex, is
        ! inside it or close enough for Newton to start from.
        xi(:) = 0.0d0

        solved = .false.
        do step = 1, MAX_NEWTON_STEPS
            call self%tabulate(coordinate_of(xi), psi=psi, dpsi_ref=reference_gradient)

            residual(1:dim) = target(1:dim)
            jacobian(1:dim, 1:dim) = 0.0d0
            do i = 1, n
                do k = 1, dim
                    residual(k) = residual(k) - psi(i) * node_coords(k, i)
                    do j = 1, dim
                        jacobian(k, j) = jacobian(k, j) + node_coords(k, i) * reference_gradient(j, i)
                    end do
                end do
            end do

            call solve_small(jacobian, residual, dim, correction, solved_step)
            if (.not. solved_step) exit

            xi(1:dim) = xi(1:dim) + correction(1:dim)
            if (maxval(abs(correction(1:dim))) < STEP_TOLERANCE) then
                solved = .true.
                exit
            end if
        end do

        if (solved) then
            call self%tabulate(coordinate_of(xi), psi=psi)
            is_in = all(psi(1:n) >= -BASIS_TOLERANCE)
            if (is_in) call normalized%set(xi(1), xi(2), xi(3))
        end if

        deallocate (psi)
        deallocate (reference_gradient)
    end subroutine is_inside_fe_petsc

    !> The reference point as the coordinate type tabulate takes.
    pure function coordinate_of(xi) result(point)
        implicit none
        real(real64), intent(in) :: xi(3)
        type(type_coordinate_dp) :: point

        point%x = xi(1)
        point%y = xi(2)
        point%z = xi(3)
    end function coordinate_of

    !> Solve a 1, 2 or 3 dimensional system by Cramer's rule.
    !>
    !> Written out rather than handed to a library: at this size the
    !> determinant is a handful of products, and the singular case has to be
    !> reported instead of stopping, because an inverted or degenerate element
    !> only means the point is not located here.
    pure subroutine solve_small(a, b, dim, x, ok)
        implicit none
        real(real64), intent(in) :: a(3, 3)
        real(real64), intent(in) :: b(3)
        integer(int32), intent(in) :: dim
        real(real64), intent(inout) :: x(3)
        logical, intent(inout) :: ok

        real(real64) :: det, inv(3, 3)
        integer(int32) :: i, j

        ok = .false.
        x(:) = 0.0d0

        select case (dim)
        case (1)
            det = a(1, 1)
            if (abs(det) <= tiny(1.0d0)) return
            x(1) = b(1) / det
        case (2)
            det = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
            if (abs(det) <= tiny(1.0d0)) return
            x(1) = (b(1) * a(2, 2) - a(1, 2) * b(2)) / det
            x(2) = (a(1, 1) * b(2) - b(1) * a(2, 1)) / det
        case (3)
            det = a(1, 1) * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2)) &
                  - a(1, 2) * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1)) &
                  + a(1, 3) * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
            if (abs(det) <= tiny(1.0d0)) return
            inv(1, 1) = (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2))
            inv(1, 2) = -(a(1, 2) * a(3, 3) - a(1, 3) * a(3, 2))
            inv(1, 3) = (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2))
            inv(2, 1) = -(a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1))
            inv(2, 2) = (a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1))
            inv(2, 3) = -(a(1, 1) * a(2, 3) - a(1, 3) * a(2, 1))
            inv(3, 1) = (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
            inv(3, 2) = -(a(1, 1) * a(3, 2) - a(1, 2) * a(3, 1))
            inv(3, 3) = (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1))
            do i = 1, 3
                do j = 1, 3
                    x(i) = x(i) + inv(i, j) * b(j) / det
                end do
            end do
        case default
            return
        end select

        ok = .true.
    end subroutine solve_small

end module domain_fe_petsc
