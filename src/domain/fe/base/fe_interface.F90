!>
!> Defines the abstract interface and common functionalities for finite elements.
!> Refactored to perform numerical integration using coefficients evaluated
!> directly at Gauss quadrature points, rather than interpolating nodal coefficients.
!>
module domain_base_fe
    use, intrinsic :: iso_fortran_env, only: int32, real64, output_unit
    use :: core_parallel_mpi
    use :: stdlib_logger
    use :: stdlib_strings
    use :: stdlib_optval, only:optval
    use :: module_core
    use :: module_linalg
    use :: domain_fe_integration

    implicit none
    private

    public :: abst_fe
    public :: holder_fes

    !>
    !> An abstract base type for finite element objects.
    !>
    type, abstract :: abst_fe
        private
        !> The element type identifier.
        integer(int32) :: type
        !> The number of nodes defining the element.
        integer(int32) :: num_nodes
        !> The geometric dimension of the element.
        integer(int32) :: dimension
        !> The gauss integration rule
        type(type_gauss_integration_rule) :: integration_rule
    contains
        !----------------------------------------------------------------------
        ! Public methods with common implementations
        !----------------------------------------------------------------------
        procedure, public, pass(self) :: initialize => initialize_abst_fe
        procedure, public, pass(self) :: destroy => destroy_abst_fe
        procedure, public, pass(self) :: get_type
        procedure, public, pass(self) :: get_num_nodes
        procedure, public, pass(self) :: get_dimension
        procedure, public, pass(self) :: get_order
        procedure, public, pass(self) :: get_num_gauss
        procedure, public, pass(self) :: get_weight
        procedure, public, pass(self) :: get_gauss
        procedure, public, pass(self) :: get_integration_rule
        procedure, public, pass(self) :: display => display_abst_fe

        !> State interpolation methods (used by assembler to get T/P at Gauss points)
        procedure, pass(self), private :: lerp_1d => lerp_1d_abst_fe
        procedure, pass(self), private :: lerp_2d => lerp_2d_abst_fe
        procedure, pass(self), private :: lerp_3d => lerp_3d_abst_fe
        generic, public :: lerp => lerp_1d, lerp_2d, lerp_3d

        procedure, pass(self), public :: dlerp => dlerp_abst_fe

        !>
        !> Calculates shape functions, physical gradients, and Jacobian determinant
        !> at a given local coordinate.
        !>
        !> Basis values and reference gradients of every node at one reference
        !> point. Every discretisation produces its whole basis at once.
        procedure(abst_tabulate), public, pass(self), deferred :: tabulate
        procedure, public, pass(self) :: calc_shape_function => calc_shape_function_abst_fe

        ! !> Integration routines taking Gauss-point values
        procedure, public, pass(self) :: compute_K1 => compute_K1_capacity_abst_fe
        procedure, public, pass(self) :: compute_K1_lumped => compute_K1_lumped_capacity_abst_fe
        procedure, pass(self), private :: compute_K2_diffusion => compute_K2_diffusion_abst_fe
        procedure, pass(self), private :: compute_K2_diffusion_scalar => compute_K2_diffusion_scalar_abst_fe
        generic, public :: compute_K2 => compute_K2_diffusion, compute_K2_diffusion_scalar
        procedure, public, pass(self) :: compute_K3 => compute_K3_mixed_abst_fe

        procedure, public, pass(self) :: compute_R1 => compute_R1_source_abst_fe
        procedure, public, pass(self) :: compute_R1_lumped => compute_R1_lumped_source_abst_fe
        procedure, public, pass(self) :: compute_R2 => compute_R2_flux_abst_fe

        !----------------------------------------------------------------------
        ! Abstract methods to be implemented in derived types
        !----------------------------------------------------------------------
        procedure(abst_calc_measure), public, pass(self), deferred :: calc_measure
        procedure(abst_is_inside), public, pass(self), deferred :: is_inside
    end type

    abstract interface
        subroutine abst_calc_measure(self, node_coords, measure)
            import :: abst_fe, int32, real64
            implicit none
            class(abst_fe), intent(in) :: self
            real(real64), intent(in) :: node_coords(:, :)
            real(real64), intent(inout) :: measure

        end subroutine abst_calc_measure

        !> Basis values and reference gradients of every node at one point.
        subroutine abst_tabulate(self, r, psi, dpsi_ref)
            import :: abst_fe, type_coordinate_dp, real64
            implicit none
            class(abst_fe), intent(in) :: self
            type(type_coordinate_dp), intent(in) :: r
            real(real64), intent(inout), optional :: psi(:)
            real(real64), intent(inout), optional :: dpsi_ref(:, :)
        end subroutine abst_tabulate

        subroutine abst_is_inside(self, cartesian, normalized, node_coords, is_in)
            import abst_fe, type_coordinate_dp, int32, real64
            implicit none
            class(abst_fe), intent(in) :: self
            type(type_coordinate_dp), intent(in) :: cartesian
            type(type_coordinate_dp), intent(inout) :: normalized
            real(real64), intent(in) :: node_coords(:, :)
            logical, intent(inout) :: is_in
        end subroutine abst_is_inside
    end interface

    interface
        module subroutine initialize_abst_fe(self, type, dimension, order, num_nodes, integration_order)
            implicit none
            class(abst_fe), intent(inout) :: self
            integer(int32), intent(in) :: type
            integer(int32), intent(in) :: dimension
            integer(int32), intent(in) :: order
            integer(int32), intent(in) :: num_nodes
            integer(int32), intent(in) :: integration_order

        end subroutine initialize_abst_fe

        module subroutine destroy_abst_fe(self)
            implicit none
            class(abst_fe), intent(inout) :: self

        end subroutine destroy_abst_fe

        module pure elemental subroutine get_type(self, type)
            implicit none
            class(abst_fe), intent(in) :: self
            integer(int32), intent(inout) :: type

        end subroutine get_type

        module pure elemental subroutine get_num_nodes(self, num_nodes)
            implicit none
            class(abst_fe), intent(in) :: self
            integer(int32), intent(inout) :: num_nodes

        end subroutine get_num_nodes

        module pure elemental subroutine get_dimension(self, dimension)
            implicit none
            class(abst_fe), intent(in) :: self
            integer(int32), intent(inout) :: dimension

        end subroutine get_dimension

        module pure elemental subroutine get_order(self, order)
            implicit none
            class(abst_fe), intent(in) :: self
            integer(int32), intent(inout) :: order

        end subroutine get_order

        module pure elemental subroutine get_num_gauss(self, num_gauss)
            implicit none
            class(abst_fe), intent(in) :: self
            integer(int32), intent(inout) :: num_gauss

        end subroutine get_num_gauss

        module subroutine get_weight(self, weight)
            implicit none
            class(abst_fe), intent(in), target :: self
            real(real64), intent(inout), pointer, contiguous, dimension(:) :: weight

        end subroutine get_weight

        module subroutine get_gauss(self, gauss)
            implicit none
            class(abst_fe), intent(in), target :: self
            type(type_coordinate_dp), intent(inout), pointer, contiguous, dimension(:) :: gauss

        end subroutine get_gauss

        module subroutine get_integration_rule(self, integration_rule)
            implicit none
            class(abst_fe), intent(in), target :: self
            type(type_gauss_integration_rule), intent(inout), pointer :: integration_rule

        end subroutine get_integration_rule

        module subroutine display_abst_fe(self, unit_in)
            implicit none
            class(abst_fe), intent(in) :: self
            integer(int32), intent(in), optional :: unit_in

        end subroutine display_abst_fe
    end interface

    interface
        module subroutine lerp_1d_abst_fe(self, r, values, lerped_value)
            implicit none
            class(abst_fe), intent(in) :: self
            type(type_coordinate_dp), intent(in) :: r
            real(real64), intent(in) :: values(:)
            real(real64), intent(inout) :: lerped_value

        end subroutine lerp_1d_abst_fe

        module subroutine lerp_2d_abst_fe(self, r, values, lerped_values)
            implicit none
            class(abst_fe), intent(in) :: self
            type(type_coordinate_dp), intent(in) :: r
            real(real64), intent(in) :: values(:, :)
            real(real64), intent(inout) :: lerped_values(:)

        end subroutine lerp_2d_abst_fe

        module subroutine lerp_3d_abst_fe(self, r, values, lerped_values)
            implicit none
            class(abst_fe), intent(in) :: self
            type(type_coordinate_dp), intent(in) :: r
            real(real64), intent(in) :: values(:, :, :)
            real(real64), intent(inout) :: lerped_values(:, :)

        end subroutine lerp_3d_abst_fe

        module subroutine dlerp_abst_fe(self, r, values, node_coords, plane_axis, dlerped_value)
            implicit none
            class(abst_fe), intent(in) :: self
            type(type_coordinate_dp), intent(in) :: r
            real(real64), intent(in) :: values(:)
            real(real64), intent(in) :: node_coords(:, :)
            integer(int32), intent(in) :: plane_axis
            type(type_coordinate_dp), intent(inout) :: dlerped_value

        end subroutine dlerp_abst_fe
    end interface

    interface
        module subroutine compute_K1_capacity_abst_fe(self, nodes, A_gp, elem_mat, work_psi)
            implicit none
            class(abst_fe), intent(in) :: self
            real(real64), intent(in) :: nodes(:, :)
            real(real64), intent(in) :: A_gp(:)
            real(real64), intent(inout) :: elem_mat(:, :)
            real(real64), intent(inout), optional, target :: work_psi(:)

        end subroutine compute_K1_capacity_abst_fe

        module subroutine compute_K1_lumped_capacity_abst_fe(self, nodes, A_gp, elem_mat, work_psi)
            implicit none
            class(abst_fe), intent(in) :: self
            real(real64), intent(in) :: nodes(:, :)
            real(real64), intent(in) :: A_gp(:)
            real(real64), intent(inout) :: elem_mat(:, :)
            real(real64), intent(inout), optional, target :: work_psi(:)

        end subroutine compute_K1_lumped_capacity_abst_fe

        module subroutine compute_K2_diffusion_abst_fe(self, nodes, M_gp, elem_mat, &
                                                       work_psi, work_dpsi_dx, work_vec)
            implicit none
            class(abst_fe), intent(in) :: self
            real(real64), intent(in) :: nodes(:, :)
            real(real64), intent(in) :: M_gp(:, :, :)
            real(real64), intent(inout) :: elem_mat(:, :)
            real(real64), intent(inout), optional, target :: work_psi(:)
            real(real64), intent(inout), optional, target :: work_dpsi_dx(:, :)
            real(real64), intent(inout), optional, target :: work_vec(:)

        end subroutine compute_K2_diffusion_abst_fe

        module subroutine compute_K2_diffusion_scalar_abst_fe(self, nodes, M_gp, elem_mat, &
                                                              work_psi, work_dpsi_dx)
            implicit none
            class(abst_fe), intent(in) :: self
            real(real64), intent(in) :: nodes(:, :)
            real(real64), intent(in) :: M_gp(:)
            real(real64), intent(inout) :: elem_mat(:, :)
            real(real64), intent(inout), optional, target :: work_psi(:)
            real(real64), intent(inout), optional, target :: work_dpsi_dx(:, :)

        end subroutine compute_K2_diffusion_scalar_abst_fe

        module subroutine compute_K3_mixed_abst_fe(self, nodes, V_gp, elem_mat, work_psi, work_dpsi_dx)
            implicit none
            class(abst_fe), intent(in) :: self
            real(real64), intent(in) :: nodes(:, :)
            real(real64), intent(in) :: V_gp(:, :)
            real(real64), intent(inout) :: elem_mat(:, :)
            real(real64), intent(inout), optional, target :: work_psi(:)
            real(real64), intent(inout), optional, target :: work_dpsi_dx(:, :)

        end subroutine compute_K3_mixed_abst_fe

        module subroutine compute_R1_source_abst_fe(self, nodes, S_gp, elem_vec, work_psi)
            implicit none
            class(abst_fe), intent(in) :: self
            real(real64), intent(in) :: nodes(:, :)
            real(real64), intent(in) :: S_gp(:)
            real(real64), intent(inout) :: elem_vec(:)
            real(real64), intent(inout), optional, target :: work_psi(:)

        end subroutine compute_R1_source_abst_fe

        module subroutine compute_R1_lumped_source_abst_fe(self, nodes, S_node, elem_vec, work_psi)
            implicit none
            class(abst_fe), intent(in) :: self
            real(real64), intent(in) :: nodes(:, :)
            real(real64), intent(in) :: S_node(:)
            real(real64), intent(inout) :: elem_vec(:)
            real(real64), intent(inout), optional, target :: work_psi(:)

        end subroutine compute_R1_lumped_source_abst_fe

        module subroutine compute_R2_flux_abst_fe(self, nodes, F_gp, elem_vec, work_dpsi_dx)
            implicit none
            class(abst_fe), intent(in) :: self
            real(real64), intent(in) :: nodes(:, :)
            real(real64), intent(in) :: F_gp(:, :)
            real(real64), intent(inout) :: elem_vec(:)
            real(real64), intent(inout), optional, target :: work_dpsi_dx(:, :)

        end subroutine compute_R2_flux_abst_fe

    end interface

    interface
        module subroutine calc_shape_function_abst_fe(self, r, node_coords, psi, dpsi_dx, inverse_jacobian, determinant_jacobian)
            implicit none
            class(abst_fe), intent(in) :: self
            type(type_coordinate_dp), intent(in) :: r
            real(real64), intent(in) :: node_coords(:, :)

            real(real64), intent(inout), optional :: psi(:)
            real(real64), intent(inout), optional :: dpsi_dx(:, :)
            real(real64), intent(inout), optional :: inverse_jacobian(:, :)
            real(real64), intent(inout), optional :: determinant_jacobian

        end subroutine calc_shape_function_abst_fe
    end interface

    !> Wrapper for polymorphic FE objects
    type :: holder_fes
        class(abst_fe), allocatable :: fe
    end type holder_fes

end module domain_base_fe
