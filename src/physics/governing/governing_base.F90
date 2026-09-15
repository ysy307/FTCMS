module physics_governing_base
    use, intrinsic :: iso_fortran_env
    use :: module_core
    use :: module_control
    use :: module_domain
    use :: module_constitutive, only:type_constitutive_manager
    implicit none
    private

    public :: type_assemble_workspace
    public :: HYDRAULIC_DRIVER_TOTAL_POTENTIAL

    !> Drive the liquid Darcy flux with p_l = -s_eff (Dall'Amico) instead of
    !> the primary unknown p_w. Every consumer of the driver must scale by
    !> d(p_l)/dX, which is not unity. Lives here because the thermal block
    !> consumes the same flux and must not reference type_hydraulic.
    logical, parameter :: HYDRAULIC_DRIVER_TOTAL_POTENTIAL = .true.

    type :: type_assemble_workspace
        logical, private :: is_initialized = .false.
        integer(int32), private :: fe_type = -1
        integer(int32) :: num_fe_nodes = -1
        integer(int32) :: num_fe_gauss = -1
        integer(int32) :: num_fe_dimension = -1

        class(abst_fe), pointer :: fe

        real(real64), allocatable :: coordinates(:, :)

        ! --- State Management ---
        type(type_state), allocatable :: state(:)
        type(type_state), allocatable :: state_gp(:)

        ! --- Nodal Values (Buffers) ---
        real(real64), allocatable :: T_node(:)
        real(real64), allocatable :: T_gp(:)
        real(real64), allocatable :: P_node(:)
        real(real64), allocatable :: P_gp(:)
        !> Nodal liquid pressure p_l = -s_eff [Pa]. Element local, so its
        !> gradient enters the residual as grad(I_h p_l), the same discrete
        !> object the hydraulic flux applies the FE operator to.
        real(real64), allocatable :: Pl_node(:)
        real(real64), allocatable :: phi_node(:)
        real(real64), allocatable :: phi_gp(:)
        real(real64), allocatable :: Qi_node(:)
        real(real64), allocatable :: Qi_gp(:)
        real(real64), allocatable :: Qw_node(:)
        real(real64), allocatable :: Qw_gp(:)

        ! --- Workspaces ---
        real(real64), allocatable :: work_node(:, :)
        real(real64), allocatable :: work_bdf_buffer(:)

        real(real64), allocatable :: work_psi(:)
        real(real64), allocatable :: work_dpsi_dx(:, :)
        real(real64), allocatable :: work_vec(:)

        real(real64), allocatable :: work_C(:)
        real(real64), allocatable :: work_D(:, :, :)
        real(real64), allocatable :: work_V(:, :)

        real(real64), allocatable :: work_L(:)
        real(real64), allocatable :: work_d_dt(:)
        real(real64), allocatable :: work_matrix(:, :)


        integer(int32) :: material_id
        integer(int32) :: element_id
        integer(int32) :: computation_type

        logical, private :: associated_bdf = .false.
        integer(int32) :: bdf_order = -1
        real(real64), pointer, contiguous, dimension(:) :: bdf_coeffs

    contains
        procedure, public, pass(self) :: initialize => initialize_type_assemble_workspace
        procedure, public, pass(self) :: destroy => destroy_type_assemble_workspace

        procedure, private, pass(self) :: set_basic => set_basic
        procedure, private, pass(self) :: set_bdf_info => set_bdf_info

        procedure, public, pass(self) :: lerp => lerp_states
        procedure, public, pass(self) :: lerp_dqi_dt => lerp_dqi_dt_from_nodes


        procedure, public, pass(self) :: compute_K1 => compute_K1_assemble_workspace
        procedure, public, pass(self) :: compute_K1_lumped => compute_K1_lumped_assemble_workspace
        procedure, private, pass(self) :: compute_K2_scalar_assemble_workspace
        procedure, private, pass(self) :: compute_K2_assemble_workspace
        generic :: compute_K2 => compute_K2_scalar_assemble_workspace, compute_K2_assemble_workspace
        procedure, public, pass(self) :: compute_K3 => compute_K3_assemble_workspace
        procedure, public, pass(self) :: compute_R1 => compute_R1_assemble_workspace
        procedure, public, pass(self) :: compute_R1_lumped => compute_R1_lumped_assemble_workspace
        procedure, public, pass(self) :: compute_R2 => compute_R2_assemble_workspace
    end type type_assemble_workspace

contains

    subroutine initialize_type_assemble_workspace(self, fe, material_id, element_id, computation_type, coordinates, control)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self
        class(abst_fe), intent(in), target :: fe
        integer(int32), intent(in) :: material_id
        integer(int32), intent(in) :: element_id
        integer(int32), intent(in) :: computation_type
        real(real64), intent(in) :: coordinates(:, :)
        type(type_control), intent(in) :: control

        integer(int32) :: fe_type

        if (.not. self%associated_bdf) then
            call self%set_bdf_info(control)
        end if

        if (.not. self%is_initialized) then
            self%fe => fe
            call self%set_basic()

            self%is_initialized = .true.
        else
            call fe%get_type(fe_type)
            if (fe_type /= self%fe_type) then
                call self%destroy()
                self%fe => fe
                if (self%bdf_order == -1) call self%set_bdf_info(control)
                call self%set_basic()
            end if
        end if

        self%material_id = material_id
        self%element_id = element_id
        self%computation_type = computation_type

        ! Reallocate coordinates only when dimensions change
        if (allocated(self%coordinates)) then
            if (size(self%coordinates, 1) /= size(coordinates, 1) .or. &
                size(self%coordinates, 2) /= size(coordinates, 2)) then
                deallocate (self%coordinates)
                allocate (self%coordinates(size(coordinates, 1), size(coordinates, 2)))
            end if
        else
            allocate (self%coordinates(size(coordinates, 1), size(coordinates, 2)))
        end if
        self%coordinates = coordinates

    end subroutine initialize_type_assemble_workspace

    subroutine set_basic(self)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self

        integer(int32) :: i

        call self%fe%get_type(self%fe_type)
        call self%fe%get_num_nodes(self%num_fe_nodes)
        call self%fe%get_num_gauss(self%num_fe_gauss)
        call self%fe%get_dimension(self%num_fe_dimension)

        ! Resize struct arrays only when element type changes
        if (allocated(self%state)) then
            if (size(self%state) /= self%num_fe_nodes) then
                deallocate (self%state)
                allocate (self%state(self%num_fe_nodes))
                do i = 1, self%num_fe_nodes
                    call self%state(i)%reset()
                end do
            else
                do i = 1, self%num_fe_nodes
                    call self%state(i)%reset()
                end do
            end if
        else
            allocate (self%state(self%num_fe_nodes))
            do i = 1, self%num_fe_nodes
                call self%state(i)%reset()
            end do
        end if

        if (allocated(self%state_gp)) then
            if (size(self%state_gp) /= self%num_fe_gauss) then
                deallocate (self%state_gp)
                allocate (self%state_gp(self%num_fe_gauss))
                do i = 1, self%num_fe_gauss
                    call self%state_gp(i)%reset()
                end do
            else
                do i = 1, self%num_fe_gauss
                    call self%state_gp(i)%reset()
                end do
            end if
        else
            allocate (self%state_gp(self%num_fe_gauss))
            do i = 1, self%num_fe_gauss
                call self%state_gp(i)%reset()
            end do
        end if

        ! Allocate workspace arrays
        call allocate_and_init(self%T_node, self%num_fe_nodes)
        call allocate_and_init(self%P_node, self%num_fe_nodes)
        call allocate_and_init(self%Pl_node, self%num_fe_nodes)
        call allocate_and_init(self%phi_node, self%num_fe_nodes)
        call allocate_and_init(self%Qi_node, self%num_fe_nodes)
        call allocate_and_init(self%Qw_node, self%num_fe_nodes)

        call allocate_and_init(self%T_gp, self%num_fe_gauss)
        call allocate_and_init(self%P_gp, self%num_fe_gauss)
        call allocate_and_init(self%phi_gp, self%num_fe_gauss)
        call allocate_and_init(self%Qi_gp, self%num_fe_gauss)
        call allocate_and_init(self%Qw_gp, self%num_fe_gauss)

        call allocate_and_init_2d(self%work_node, self%bdf_order + 1, self%num_fe_nodes)
        call allocate_and_init(self%work_bdf_buffer, self%bdf_order + 1)
        call allocate_and_init(self%work_psi, self%num_fe_nodes)
        call allocate_and_init_2d(self%work_dpsi_dx, self%num_fe_dimension, self%num_fe_nodes)

        call allocate_and_init(self%work_vec, self%num_fe_nodes)
        call allocate_and_init(self%work_C, self%num_fe_gauss)
        call allocate_and_init_3d(self%work_D, self%num_fe_dimension, self%num_fe_dimension, self%num_fe_gauss)
        call allocate_and_init_2d(self%work_V, self%num_fe_dimension, self%num_fe_gauss)

        call allocate_and_init(self%work_L, self%num_fe_gauss)
        call allocate_and_init(self%work_d_dt, self%num_fe_gauss)
        call allocate_and_init_2d(self%work_matrix, self%num_fe_nodes, self%num_fe_nodes)

    contains
        ! Thread-safe allocation helpers
        subroutine allocate_and_init(arr, sz)
            real(real64), allocatable, intent(inout) :: arr(:)
            integer(int32), intent(in) :: sz

            if (allocated(arr)) then
                if (size(arr) == sz) then
                    arr = 0.0d0
                    return
                end if
                deallocate (arr)
            end if
            allocate (arr(sz))
            arr = 0.0d0
        end subroutine allocate_and_init

        subroutine allocate_and_init_2d(arr, d1, d2)
            real(real64), allocatable, intent(inout) :: arr(:, :)
            integer(int32), intent(in) :: d1, d2

            if (allocated(arr)) then
                if (size(arr, 1) == d1 .and. size(arr, 2) == d2) then
                    arr = 0.0d0
                    return
                end if
                deallocate (arr)
            end if
            allocate (arr(d1, d2))
            arr = 0.0d0
        end subroutine allocate_and_init_2d

        subroutine allocate_and_init_3d(arr, d1, d2, d3)
            real(real64), allocatable, intent(inout) :: arr(:, :, :)
            integer(int32), intent(in) :: d1, d2, d3

            if (allocated(arr)) then
                if (size(arr, 1) == d1 .and. size(arr, 2) == d2 .and. size(arr, 3) == d3) then
                    arr = 0.0d0
                    return
                end if
                deallocate (arr)
            end if
            allocate (arr(d1, d2, d3))
            arr = 0.0d0
        end subroutine allocate_and_init_3d

    end subroutine set_basic

    subroutine set_bdf_info(self, control)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self
        type(type_control), intent(in) :: control

        call control%get_bdf_coeffs(self%bdf_order, self%bdf_coeffs)

        self%associated_bdf = .true.
    end subroutine set_bdf_info

    subroutine lerp_states(self)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self

        integer(int32) :: i, j, k
        type(type_coordinate_dp), pointer, contiguous, dimension(:) :: gp
        type(type_coordinate_dp) :: dlerped_value

        real(real64), pointer, contiguous, dimension(:) :: work_history_ptr
        real(real64) :: work_value

        logical :: is_set
        integer(int32) :: history_len
        nullify (gp)
        nullify (work_history_ptr)

        call self%fe%get_gauss(gp)

        ! Initialize all variables to safe defaults
        self%T_node(:) = 273.15d0
        self%P_node(:) = 0.0d0
        self%Pl_node(:) = 0.0d0
        self%phi_node(:) = 0.0d0
        self%Qi_node(:) = 0.0d0
        self%Qw_node(:) = 0.0d0
        self%work_bdf_buffer(:) = 0.0d0

        ! 1. Temperature (THERMAL)
        if (self%associated_bdf .and. allocated(self%work_bdf_buffer)) then
            self%work_node(:, :) = 0.0d0

            do i = 1, self%num_fe_nodes
                nullify (work_history_ptr)
                is_set = .false.
                call self%state(i)%temperature%get(work_value, is_set=is_set)
                if (is_set) self%T_node(i) = work_value

                is_set = .false.
                call self%state(i)%temperature_history%get(work_history_ptr, is_set=is_set)

                if (associated(work_history_ptr) .and. is_set) then
                    history_len = min(size(work_history_ptr), self%bdf_order + 1)
                    self%work_node(1:history_len, i) = work_history_ptr(1:history_len)
                end if
            end do

            do i = 1, self%num_fe_gauss
                call self%fe%lerp(gp(i), self%T_node(1:self%num_fe_nodes), self%T_gp(i))
                call self%state_gp(i)%temperature%set(self%T_gp(i))

                call self%fe%dlerp(gp(i), self%T_node(1:self%num_fe_nodes), self%coordinates, self%computation_type, dlerped_value)
                call self%state_gp(i)%grad_T%set(dlerped_value)
            end do

            do j = 1, self%num_fe_gauss
                self%work_bdf_buffer(:) = 0.0d0
                do k = 1, self%bdf_order + 1
                    call self%fe%lerp(gp(j), self%work_node(k, 1:self%num_fe_nodes), self%work_bdf_buffer(k))
                end do
                call self%state_gp(j)%temperature_history%set(self%work_bdf_buffer)
            end do
        end if

        ! 2. Pressure (HYDRAULIC)
        self%work_node(:, :) = 0.0d0

        do i = 1, self%num_fe_nodes
            nullify (work_history_ptr)
            is_set = .false.
            call self%state(i)%pressure%get(work_value, is_set=is_set)
            if (is_set) self%P_node(i) = work_value

            is_set = .false.
            call self%state(i)%pressure_history%get(work_history_ptr, is_set=is_set)

            if (associated(work_history_ptr) .and. is_set) then
                history_len = min(size(work_history_ptr), self%bdf_order + 1)
                self%work_node(1:history_len, i) = work_history_ptr(1:history_len)
            end if
        end do

        do i = 1, self%num_fe_nodes
            is_set = .false.
            call self%state(i)%effective_suction%get(work_value, is_set=is_set)
            if (is_set) then
                self%Pl_node(i) = -work_value
            else
                self%Pl_node(i) = self%P_node(i)
            end if
        end do

        do i = 1, self%num_fe_gauss
            call self%fe%lerp(gp(i), self%P_node(1:self%num_fe_nodes), self%P_gp(i))
            call self%state_gp(i)%pressure%set(self%P_gp(i))

            call self%fe%dlerp(gp(i), self%P_node(1:self%num_fe_nodes), self%coordinates, self%computation_type, dlerped_value)
            call self%state_gp(i)%grad_P%set(dlerped_value)

            call self%fe%dlerp(gp(i), self%Pl_node(1:self%num_fe_nodes), self%coordinates, self%computation_type, dlerped_value)
            call self%state_gp(i)%grad_Pl%set(dlerped_value)
        end do

        do j = 1, self%num_fe_gauss
            self%work_bdf_buffer(:) = 0.0d0
            do k = 1, self%bdf_order + 1
                call self%fe%lerp(gp(j), self%work_node(k, 1:self%num_fe_nodes), self%work_bdf_buffer(k))
            end do
            call self%state_gp(j)%pressure_history%set(self%work_bdf_buffer)
        end do

        ! 3. Porosity
        self%work_node(:, :) = 0.0d0

        do i = 1, self%num_fe_nodes
            nullify (work_history_ptr)
            is_set = .false.
            call self%state(i)%porosity%get(work_value, is_set=is_set)
            if (is_set) self%phi_node(i) = work_value

            is_set = .false.
            call self%state(i)%porosity_history%get(work_history_ptr, is_set=is_set)

            if (associated(work_history_ptr) .and. is_set) then
                history_len = min(size(work_history_ptr), self%bdf_order + 1)
                self%work_node(1:history_len, i) = work_history_ptr(1:history_len)
            end if
        end do

        do i = 1, self%num_fe_gauss
            call self%fe%lerp(gp(i), self%phi_node(1:self%num_fe_nodes), self%phi_gp(i))
            call self%state_gp(i)%porosity%set(self%phi_gp(i))
        end do

        do j = 1, self%num_fe_gauss
            self%work_bdf_buffer(:) = 0.0d0
            do k = 1, self%bdf_order + 1
                call self%fe%lerp(gp(j), self%work_node(k, 1:self%num_fe_nodes), self%work_bdf_buffer(k))
            end do
            call self%state_gp(j)%porosity_history%set(self%work_bdf_buffer)
        end do

        ! 4. Outer ice state. The current value and its BDF history are
        ! interpolated independently of the monolithic T-p unknown vector.
        self%work_node(:, :) = 0.0d0

        do i = 1, self%num_fe_nodes
            nullify (work_history_ptr)
            is_set = .false.
            call self%state(i)%ice_content%get(work_value, is_set=is_set)
            if (is_set) self%Qi_node(i) = work_value

            is_set = .false.
            call self%state(i)%ice_content_history%get(work_history_ptr, is_set=is_set)
            if (associated(work_history_ptr) .and. is_set) then
                history_len = min(size(work_history_ptr), self%bdf_order + 1)
                self%work_node(1:history_len, i) = work_history_ptr(1:history_len)
            end if
        end do

        do i = 1, self%num_fe_gauss
            call self%fe%lerp(gp(i), self%Qi_node(1:self%num_fe_nodes), self%Qi_gp(i))
            call self%state_gp(i)%ice_content%set(self%Qi_gp(i))
        end do

        ! Liquid water is interpolated too: together with the ice it carries the
        ! total water content the phase split redistributes at fixed mass. Only
        ! the ice used to be interpolated, because the liquid was rebuilt from
        ! the pressure - which is exactly the closure that had to go.
        do i = 1, self%num_fe_nodes
            is_set = .false.
            call self%state(i)%water_content%get(work_value, is_set=is_set)
            if (is_set) self%Qw_node(i) = work_value
        end do
        do i = 1, self%num_fe_gauss
            call self%fe%lerp(gp(i), self%Qw_node(1:self%num_fe_nodes), self%Qw_gp(i))
            call self%state_gp(i)%water_content%set(self%Qw_gp(i))
        end do

        do j = 1, self%num_fe_gauss
            self%work_bdf_buffer(:) = 0.0d0
            do k = 1, self%bdf_order + 1
                call self%fe%lerp(gp(j), self%work_node(k, 1:self%num_fe_nodes), self%work_bdf_buffer(k))
            end do
            call self%state_gp(j)%ice_content_history%set(self%work_bdf_buffer)
        end do

        ! The continuation parameter is a scalar of the solve, not a field, so
        ! it is carried across rather than interpolated.
        is_set = .false.
        call self%state(1)%continuation_lambda%get(work_value, is_set=is_set)
        if (is_set) then
            do i = 1, self%num_fe_gauss
                call self%state_gp(i)%continuation_lambda%set(work_value)
            end do
        end if

        nullify (gp)
        nullify (work_history_ptr)

    end subroutine lerp_states

    ! Interpolate dQi_dT from node states to Gauss point states.
    ! Called after update_physical_properties_bulk to fix the phase-transition
    ! front element: Gauss point T > T_f0 gives dQi_dT=0 (no latent heat),
    ! but node-interpolated values correctly carry latent heat from frozen nodes.
    subroutine lerp_dqi_dt_from_nodes(self)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self

        integer(int32) :: i
        type(type_coordinate_dp), pointer, contiguous, dimension(:) :: gp
        real(real64) :: dqi_dt_gp

        nullify (gp)
        call self%fe%get_gauss(gp)
        do i = 1, self%num_fe_nodes
            call self%state(i)%dQi_dT%get(self%work_vec(i))
        end do

        do i = 1, self%num_fe_gauss
            call self%fe%lerp(gp(i), self%work_vec(1:self%num_fe_nodes), dqi_dt_gp)
            call self%state_gp(i)%dQi_dT%set(dqi_dt_gp)
        end do

        nullify (gp)

    end subroutine lerp_dqi_dt_from_nodes

    subroutine compute_K1_assemble_workspace(self, A_gp, local_matrix)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self
        real(real64), intent(in) :: A_gp(:)
        real(real64), intent(inout) :: local_matrix(:, :)

        local_matrix(:, :) = 0.0d0
        self%work_psi(:) = 0.0d0
        call self%fe%compute_K1(self%coordinates, A_gp, local_matrix, self%work_psi)
    end subroutine compute_K1_assemble_workspace

    subroutine compute_K1_lumped_assemble_workspace(self, A_gp, local_matrix)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self
        real(real64), intent(in) :: A_gp(:)
        real(real64), intent(inout) :: local_matrix(:, :)

        local_matrix(:, :) = 0.0d0
        self%work_psi(:) = 0.0d0
        call self%fe%compute_K1_lumped(self%coordinates, A_gp, local_matrix, self%work_psi)
    end subroutine compute_K1_lumped_assemble_workspace

    subroutine compute_K2_scalar_assemble_workspace(self, D_gp, local_matrix)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self
        real(real64), intent(in) :: D_gp(:)
        real(real64), intent(inout) :: local_matrix(:, :)

        local_matrix(:, :) = 0.0d0
        self%work_psi(:) = 0.0d0
        self%work_dpsi_dx(:, :) = 0.0d0
        call self%fe%compute_K2(self%coordinates, D_gp, local_matrix, self%work_psi, self%work_dpsi_dx)
    end subroutine compute_K2_scalar_assemble_workspace

    subroutine compute_K2_assemble_workspace(self, D_gp, local_matrix)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self
        real(real64), intent(in) :: D_gp(:, :, :)
        real(real64), intent(inout) :: local_matrix(:, :)

        local_matrix(:, :) = 0.0d0
        self%work_psi(:) = 0.0d0
        self%work_dpsi_dx(:, :) = 0.0d0
        self%work_vec(:) = 0.0d0
        call self%fe%compute_K2(self%coordinates, D_gp, local_matrix, self%work_psi, self%work_dpsi_dx, self%work_vec)
    end subroutine compute_K2_assemble_workspace

    subroutine compute_K3_assemble_workspace(self, V_gp, local_matrix)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self
        real(real64), intent(in) :: V_gp(:, :)
        real(real64), intent(inout) :: local_matrix(:, :)

        local_matrix(:, :) = 0.0d0
        self%work_psi(:) = 0.0d0
        self%work_dpsi_dx(:, :) = 0.0d0
        call self%fe%compute_K3(self%coordinates, V_gp, local_matrix, self%work_psi, self%work_dpsi_dx)
    end subroutine compute_K3_assemble_workspace

    subroutine compute_R1_assemble_workspace(self, S_gp, local_vector)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self
        real(real64), intent(in) :: S_gp(:)
        real(real64), intent(inout) :: local_vector(:)

        local_vector(:) = 0.0d0
        self%work_psi(:) = 0.0d0
        call self%fe%compute_R1(self%coordinates, S_gp, local_vector, self%work_psi)
    end subroutine compute_R1_assemble_workspace

    subroutine compute_R1_lumped_assemble_workspace(self, S_node, local_vector)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self
        real(real64), intent(in) :: S_node(:)
        real(real64), intent(inout) :: local_vector(:)

        local_vector(:) = 0.0d0
        self%work_psi(:) = 0.0d0
        call self%fe%compute_R1_lumped(self%coordinates, S_node, local_vector, self%work_psi)
    end subroutine compute_R1_lumped_assemble_workspace

    subroutine compute_R2_assemble_workspace(self, V_gp, local_vector)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self
        real(real64), intent(in) :: V_gp(:, :)
        real(real64), intent(inout) :: local_vector(:)

        local_vector(:) = 0.0d0
        self%work_psi(:) = 0.0d0
        self%work_dpsi_dx(:, :) = 0.0d0
        call self%fe%compute_R2(self%coordinates, V_gp, local_vector, self%work_dpsi_dx)
    end subroutine compute_R2_assemble_workspace

    subroutine destroy_type_assemble_workspace(self)
        implicit none
        class(type_assemble_workspace), intent(inout) :: self

        if (self%is_initialized) then
            self%fe_type = -1
            self%num_fe_nodes = -1
            self%num_fe_gauss = -1
            self%num_fe_dimension = -1

            if (allocated(self%state)) deallocate (self%state)
            if (allocated(self%state_gp)) deallocate (self%state_gp)

            if (allocated(self%T_node)) deallocate (self%T_node)
            if (allocated(self%P_node)) deallocate (self%P_node)
            if (allocated(self%phi_node)) deallocate (self%phi_node)
            if (allocated(self%Qi_node)) deallocate (self%Qi_node)
            if (allocated(self%Qw_node)) deallocate (self%Qw_node)
            if (allocated(self%T_gp)) deallocate (self%T_gp)
            if (allocated(self%P_gp)) deallocate (self%P_gp)
            if (allocated(self%phi_gp)) deallocate (self%phi_gp)
            if (allocated(self%Qi_gp)) deallocate (self%Qi_gp)
            if (allocated(self%Qw_gp)) deallocate (self%Qw_gp)
            if (allocated(self%coordinates)) deallocate (self%coordinates)

            if (allocated(self%work_node)) deallocate (self%work_node)
            if (allocated(self%work_bdf_buffer)) deallocate (self%work_bdf_buffer)
            if (allocated(self%work_psi)) deallocate (self%work_psi)
            if (allocated(self%work_dpsi_dx)) deallocate (self%work_dpsi_dx)
            if (allocated(self%work_vec)) deallocate (self%work_vec)

            if (allocated(self%work_C)) deallocate (self%work_C)
            if (allocated(self%work_D)) deallocate (self%work_D)
            if (allocated(self%work_V)) deallocate (self%work_V)
            if (allocated(self%work_L)) deallocate (self%work_L)
            if (allocated(self%work_d_dt)) deallocate (self%work_d_dt)
            if (allocated(self%work_matrix)) deallocate (self%work_matrix)

            nullify (self%fe)
            self%is_initialized = .false.

            self%material_id = -1
            self%element_id = -1
            self%computation_type = -1

            self%associated_bdf = .false.
            self%bdf_order = -1
            nullify (self%bdf_coeffs)
        end if

    end subroutine destroy_type_assemble_workspace

end module physics_governing_base
