submodule(domain_base_fe) domain_base_fe_coefficients
    implicit none
contains

    module subroutine calc_shape_function_abst_fe(self, r, node_coords, psi, dpsi_dx, inverse_jacobian, determinant_jacobian)
        implicit none
        class(abst_fe), intent(in) :: self
        type(type_coordinate_dp), intent(in) :: r
        real(real64), intent(in) :: node_coords(:, :)

        real(real64), intent(inout), optional :: psi(:)
        real(real64), intent(inout), optional :: dpsi_dx(:, :)
        real(real64), intent(inout), optional :: inverse_jacobian(:, :)
        real(real64), intent(inout), optional :: determinant_jacobian

        integer(int32) :: i, j, k, dim
        ! Use fixed-size arrays on the stack, allocated for max dimension and controlled by dim
        real(real64) :: local_J(self%dimension, self%dimension)
        real(real64) :: local_inv_J(self%dimension, self%dimension)
        real(real64) :: local_det_J
        real(real64) :: dpsi_dxi(self%dimension)
        real(real64) :: reference_gradient(self%dimension, self%num_nodes)
        real(real64) :: mapping(3, self%dimension)
        integer(int32) :: space_dim
        logical :: need_jacobian, need_inverse

        dim = self%dimension

        ! 1. Evaluate Shape Functions (Requires no Jacobian)
        if (present(psi)) then
            psi(:) = 0.0d0
            call self%tabulate(r, psi=psi)
        end if

        ! Determine whether Jacobian computation is needed
        need_jacobian = present(determinant_jacobian) .or. &
                        present(inverse_jacobian) .or. &
                        present(dpsi_dx)

        if (.not. need_jacobian) return

        ! --- A. Jacobian, from the same reference gradients the basis gave ---
        ! J(k,l) = sum_i x(k,i) dpsi_i/dxi_l. Built here rather than asked of the
        ! element, so an element only has to know its basis, not how the mesh is
        ! mapped. An element of lower dimension than the space it sits in - a
        ! side of a plane problem - has no square Jacobian, so its measure comes
        ! from the metric sqrt(det(G^T G)) instead, which is the curve length.
        call self%tabulate(r, dpsi_ref=reference_gradient)
        space_dim = size(node_coords, 1)
        mapping(:, :) = 0.0d0
        do i = 1, self%num_nodes
            do k = 1, min(space_dim, 3)
                do j = 1, dim
                    mapping(k, j) = mapping(k, j) + node_coords(k, i) * reference_gradient(j, i)
                end do
            end do
        end do
        local_J(1:dim, 1:dim) = 0.0d0
        if (space_dim == dim) then
            local_J(1:dim, 1:dim) = mapping(1:dim, 1:dim)
        else
            do j = 1, dim
                do k = 1, dim
                    local_J(k, j) = dot_product(mapping(1:space_dim, k), mapping(1:space_dim, j))
                end do
            end do
            ! The metric's determinant is the square of the measure, so take the
            ! root now and keep the same meaning downstream.
            if (dim == 1) local_J(1, 1) = sqrt(max(local_J(1, 1), 0.0d0))
        end if

        ! --- B. Determinant (analytic for dim<=3; avoids LAPACK overhead on tiny J) ---
        if (dim == 1) then
            local_det_J = local_J(1, 1)
        else if (dim == 2) then
            local_det_J = local_J(1, 1) * local_J(2, 2) - local_J(1, 2) * local_J(2, 1)
        else
            local_det_J = local_J(1, 1) * (local_J(2, 2) * local_J(3, 3) - local_J(2, 3) * local_J(3, 2)) &
                          - local_J(1, 2) * (local_J(2, 1) * local_J(3, 3) - local_J(2, 3) * local_J(3, 1)) &
                          + local_J(1, 3) * (local_J(2, 1) * local_J(3, 2) - local_J(2, 2) * local_J(3, 1))
        end if
        if (local_det_J <= 0.0_real64) then
            error stop "Non-positive Jacobian determinant: degenerate or inverted element"
        end if
        if (present(determinant_jacobian)) determinant_jacobian = local_det_J

        ! Determine whether the inverse Jacobian is needed (also required for dpsi_dx)
        need_inverse = present(inverse_jacobian) .or. present(dpsi_dx)

        if (.not. need_inverse) return

        ! --- C. Inverse Jacobian (analytic for dim<=3) ---
        if (dim == 1) then
            local_inv_J(1, 1) = 1.0d0 / local_det_J
        else if (dim == 2) then
            local_inv_J(1, 1) = local_J(2, 2) / local_det_J
            local_inv_J(1, 2) = -local_J(1, 2) / local_det_J
            local_inv_J(2, 1) = -local_J(2, 1) / local_det_J
            local_inv_J(2, 2) = local_J(1, 1) / local_det_J
        else
            local_inv_J(1, 1) = (local_J(2, 2) * local_J(3, 3) - local_J(2, 3) * local_J(3, 2)) / local_det_J
            local_inv_J(1, 2) = -(local_J(1, 2) * local_J(3, 3) - local_J(1, 3) * local_J(3, 2)) / local_det_J
            local_inv_J(1, 3) = (local_J(1, 2) * local_J(2, 3) - local_J(1, 3) * local_J(2, 2)) / local_det_J
            local_inv_J(2, 1) = -(local_J(2, 1) * local_J(3, 3) - local_J(2, 3) * local_J(3, 1)) / local_det_J
            local_inv_J(2, 2) = (local_J(1, 1) * local_J(3, 3) - local_J(1, 3) * local_J(3, 1)) / local_det_J
            local_inv_J(2, 3) = -(local_J(1, 1) * local_J(2, 3) - local_J(1, 3) * local_J(2, 1)) / local_det_J
            local_inv_J(3, 1) = (local_J(2, 1) * local_J(3, 2) - local_J(2, 2) * local_J(3, 1)) / local_det_J
            local_inv_J(3, 2) = -(local_J(1, 1) * local_J(3, 2) - local_J(1, 2) * local_J(3, 1)) / local_det_J
            local_inv_J(3, 3) = (local_J(1, 1) * local_J(2, 2) - local_J(1, 2) * local_J(2, 1)) / local_det_J
        end if

        if (present(inverse_jacobian)) then
            inverse_jacobian(1:dim, 1:dim) = local_inv_J(1:dim, 1:dim)
        end if

        ! --- D. Compute Global Gradients (dpsi_dx) ---
        ! Formula: dpsi_dx = J^{-T} * dpsi_dxi
        if (present(dpsi_dx)) then
            dpsi_dx(:, :) = 0.0d0

            do i = 1, self%num_nodes
                dpsi_dxi(1:dim) = reference_gradient(1:dim, i)

                ! Transform to physical gradients: multiply by J^{-T}
                ! dpsi_dx(k, i) = sum_l ( invJ(l, k) * dpsi_dxi(l) )
                ! Note: invJ(l, k) uses transposed access, hence the (l, k) indexing
                do k = 1, dim
                    do j = 1, dim
                        dpsi_dx(k, i) = dpsi_dx(k, i) + local_inv_J(j, k) * dpsi_dxi(j)
                    end do
                end do
            end do
        end if

    end subroutine calc_shape_function_abst_fe

end submodule domain_base_fe_coefficients
