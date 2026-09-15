submodule(domain_base_fe) domain_base_fe_interpolation
    implicit none
contains

    module subroutine lerp_1d_abst_fe(self, r, values, lerped_value)
        implicit none
        class(abst_fe), intent(in) :: self
        type(type_coordinate_dp), intent(in) :: r
        real(real64), intent(in) :: values(:)
        real(real64), intent(inout) :: lerped_value
        integer(int32) :: i
        ! One tabulation of the whole basis. Asking for one function at a time
        ! costs a round trip per node, and an element whose basis comes from
        ! outside cannot answer that question at all.
        real(real64) :: psi(self%num_nodes)

        psi(:) = 0.0d0
        call self%tabulate(r, psi=psi)

        lerped_value = 0.0d0
        do i = 1, self%num_nodes
            lerped_value = lerped_value + psi(i) * values(i)
        end do
    end subroutine lerp_1d_abst_fe

    module subroutine lerp_2d_abst_fe(self, r, values, lerped_values)
        implicit none
        class(abst_fe), intent(in) :: self
        type(type_coordinate_dp), intent(in) :: r
        real(real64), intent(in) :: values(:, :)
        real(real64), intent(inout) :: lerped_values(:)
        integer(int32) :: i
        real(real64) :: psi(self%num_nodes)

        psi(:) = 0.0d0
        call self%tabulate(r, psi=psi)

        lerped_values(:) = 0.0d0
        do i = 1, self%num_nodes
            lerped_values(:) = lerped_values(:) + psi(i) * values(:, i)
        end do
    end subroutine lerp_2d_abst_fe

    module subroutine lerp_3d_abst_fe(self, r, values, lerped_values)
        implicit none
        class(abst_fe), intent(in) :: self
        type(type_coordinate_dp), intent(in) :: r
        real(real64), intent(in) :: values(:, :, :)
        real(real64), intent(inout) :: lerped_values(:, :)
        integer(int32) :: i
        real(real64) :: psi(self%num_nodes)

        psi(:) = 0.0d0
        call self%tabulate(r, psi=psi)

        lerped_values = 0.0d0
        do i = 1, self%num_nodes
            lerped_values(:, :) = lerped_values(:, :) + psi(i) * values(:, :, i)
        end do
    end subroutine lerp_3d_abst_fe

    !---------------------------------------------------------------------------
    !> Computes the gradient of a field in physical coordinates (nabla u).
    !> Internally calls dpsi_dx (physical gradients of shape functions) and
    !> forms the linear combination.
    !>
    !> @param[in]    self          Element object
    !> @param[in]    r             Local coordinate (xi, eta, zeta)
    !> @param[in]    values        Nodal value array (u_i) [num_nodes]
    !> @param[in]    node_coords   Nodal coordinate array (x, y, z)
    !> @param[in]    plane_axis    Axis specification for 2D (1: XY-plane, 2: XZ-plane). Ignored in 3D.
    !> @param[inout] dlerped_value Computed gradient vector (du/dx, du/dy, du/dz)
    !---------------------------------------------------------------------------
    module subroutine dlerp_abst_fe(self, r, values, node_coords, plane_axis, dlerped_value)
        implicit none
        class(abst_fe), intent(in) :: self
        type(type_coordinate_dp), intent(in) :: r
        real(real64), intent(in) :: values(:)
        real(real64), intent(in) :: node_coords(:, :)
        integer(int32), intent(in) :: plane_axis
        type(type_coordinate_dp), intent(inout) :: dlerped_value

        ! Temporary array for physical gradients of shape functions (dim x num_nodes).
        ! Note: calc_dpsi_dx fills (1:dim, :) when dim=2.
        real(real64) :: shape_grads(self%dimension, self%num_nodes)
        real(real64) :: dlerp_array(self%dimension)
        integer(int32) :: ierr

        ! 1. Initialize
        call dlerped_value%reset()
        shape_grads(:, :) = 0.0d0
        dlerp_array(:) = 0.0d0

        ! 2. Get physical gradients of shape functions (dN_i/dx, dN_i/dy, dN_i/dz)
        call self%calc_shape_function(r, node_coords, dpsi_dx=shape_grads)

        call matvec(shape_grads, values, dlerp_array, ierr)

        ! 3. Compute gradient: grad u = sum( u_i * grad N_i )
        ! Map from shape_grads (always packed into 1:dim) to physical coordinate components.
        if (self%dimension == 3) then
            ! --- 3D (XYZ) ---
            ! shape_grads(1:3, :) contains x, y, z components in order
            dlerped_value%x = dlerp_array(1)
            dlerped_value%y = dlerp_array(2)
            dlerped_value%z = dlerp_array(3)
        else if (self%dimension == 2) then
            ! --- 2D (XY or XZ) ---
            ! shape_grads(1, :) -> 1st component (x)
            ! shape_grads(2, :) -> 2nd component (y or z)

            ! 1st component is always x
            dlerped_value%x = dlerp_array(1)

            if (plane_axis == 2) then
                dlerped_value%y = 0.0d0
                dlerped_value%z = dlerp_array(2)
            else
                dlerped_value%y = dlerp_array(2)
                dlerped_value%z = 0.0d0
            end if
        else if (self%dimension == 1) then
            dlerped_value%x = dlerp_array(1)
            dlerped_value%y = 0.0d0
            dlerped_value%z = 0.0d0
        end if

    end subroutine dlerp_abst_fe
end submodule domain_base_fe_interpolation
