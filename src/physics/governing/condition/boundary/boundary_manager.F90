module physics_governing_boundary_manager
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: module_core
    use :: boundary_data_provider
    use :: boundary_strategy
    use :: boundary_strategy_factory
    implicit none
    private

    public :: type_bc_manager

    type :: holder_bc_strategy
        class(abst_bc), allocatable :: p
    end type holder_bc_strategy

    type :: type_bc_manager
        integer(int32), private :: num_boundaries = 0
        type(holder_bc_strategy), allocatable, private :: strategies(:)
        !> JSON boundary_id for each strategy (same order as strategies(:)).
        integer(int32), allocatable, private :: boundary_ids(:)
    contains
        procedure, public, pass(self) :: initialize => initialize_bc_manager
        procedure, public, pass(self) :: destroy => destroy_bc_manager
        generic :: evaluate => evaluate_bc, evaluate_bcs
        procedure, public, pass(self) :: evaluate_bc => evaluate_bc_manager
        procedure, public, pass(self) :: evaluate_bcs => evaluate_bcs_manager
        procedure, public, pass(self) :: get_num_boundaries => get_num_boundaries_bc_manager
        procedure, public, pass(self) :: get_bc_index => get_bc_index_bc_manager
        procedure, public, pass(self) :: get_boundary_id => get_boundary_id_bc_manager
        procedure, public, pass(self) :: get_data => get_data_bc_manager
        procedure, public, pass(self) :: update_bc_data => update_bc_data_manager
    end type type_bc_manager

contains

    subroutine initialize_bc_manager(self, configs)
        implicit none
        class(type_bc_manager), intent(inout) :: self
        type(type_config_bc), intent(in) :: configs(:)

        integer(int32) :: i

        self%num_boundaries = size(configs)

        if (allocated(self%strategies)) deallocate (self%strategies)
        if (allocated(self%boundary_ids)) deallocate (self%boundary_ids)

        if (self%num_boundaries > 0) then
            allocate (self%strategies(self%num_boundaries))
            allocate (self%boundary_ids(self%num_boundaries))

            do i = 1, self%num_boundaries
                self%strategies(i)%p = create_bc_strategy(configs(i))
                self%boundary_ids(i) = configs(i)%boundary_id
            end do
        end if

    end subroutine initialize_bc_manager

    subroutine destroy_bc_manager(self)
        implicit none
        class(type_bc_manager), intent(inout) :: self
        integer(int32) :: i

        if (allocated(self%strategies)) then
            do i = 1, size(self%strategies)
                if (allocated(self%strategies(i)%p)) then
                    call self%strategies(i)%p%destroy()
                    deallocate (self%strategies(i)%p)
                end if
            end do
            deallocate (self%strategies)
        end if
        if (allocated(self%boundary_ids)) deallocate (self%boundary_ids)

    end subroutine destroy_bc_manager

    subroutine evaluate_bc_manager(self, bc_id, current_time, u_curr, result)
        implicit none
        class(type_bc_manager), intent(inout) :: self
        integer(int32), intent(in) :: bc_id
        real(real64), intent(in) :: current_time
        real(real64), intent(in) :: u_curr
        type(type_bc_result), intent(inout) :: result

        if (bc_id > 0 .and. bc_id <= size(self%strategies)) then
            call self%strategies(bc_id)%p%evaluate(current_time, u_curr, result)
        else
            call result%initialize()
        end if
    end subroutine evaluate_bc_manager

    subroutine evaluate_bcs_manager(self, bc_id, current_time, u_curr, results)
        implicit none
        class(type_bc_manager), intent(inout) :: self
        integer(int32), intent(in) :: bc_id
        real(real64), intent(in) :: current_time
        real(real64), intent(in) :: u_curr(:)
        type(type_bc_result), intent(inout) :: results(:)

        integer(int32) :: i

        if (size(u_curr) /= size(results)) then
            error stop "Size of u_curr and results must be the same in evaluate_bc_manager."
        end if

        if (bc_id > 0 .and. bc_id <= size(self%strategies)) then
            do i = 1, size(u_curr)
                call self%strategies(bc_id)%p%evaluate(current_time, u_curr(i), results(i))
            end do
        else
            do i = 1, size(u_curr)
                call results(i)%initialize()
            end do
        end if
    end subroutine evaluate_bcs_manager

    pure subroutine get_num_boundaries_bc_manager(self, num_boundaries)
        implicit none
        class(type_bc_manager), intent(in) :: self
        integer(int32), intent(inout) :: num_boundaries

        num_boundaries = self%num_boundaries
    end subroutine get_num_boundaries_bc_manager

    subroutine update_bc_data_manager(self, entity_id, new_values)
        implicit none
        class(type_bc_manager), intent(inout) :: self
        integer(int32), intent(in) :: entity_id
        real(real64), intent(in) :: new_values(:)

        integer(int32) :: bc_idx

        call self%get_bc_index(entity_id, bc_idx)
        if (bc_idx > 0 .and. bc_idx <= self%num_boundaries) then
            if (allocated(self%strategies(bc_idx)%p)) then
                call self%strategies(bc_idx)%p%update_provider_data(new_values)
            end if
        end if
    end subroutine update_bc_data_manager

    !> Returns the positional strategy index whose boundary_id matches entity_id.
    !> Returns -1 if no match is found.
    pure subroutine get_bc_index_bc_manager(self, entity_id, bc_index)
        implicit none
        class(type_bc_manager), intent(in) :: self
        integer(int32), intent(in) :: entity_id
        integer(int32), intent(inout) :: bc_index

        integer(int32) :: i

        bc_index = -1
        if (.not. allocated(self%boundary_ids)) return

        do i = 1, self%num_boundaries
            if (self%boundary_ids(i) == entity_id) then
                bc_index = i
                return
            end if
        end do
    end subroutine get_bc_index_bc_manager

    !> The provider data of the condition at index `bc_index`, at `current_time`.
    subroutine get_data_bc_manager(self, bc_index, current_time, values)
        implicit none
        class(type_bc_manager), intent(inout) :: self
        integer(int32), intent(in) :: bc_index
        real(real64), intent(in) :: current_time
        real(real64), intent(inout) :: values(BC_DATA_SIZE)

        call self%strategies(bc_index)%p%get_data(current_time, values)
    end subroutine get_data_bc_manager

    !> The JSON boundary id of the condition at index `bc_index`.
    pure function get_boundary_id_bc_manager(self, bc_index) result(boundary_id)
        implicit none
        class(type_bc_manager), intent(in) :: self
        integer(int32), intent(in) :: bc_index
        integer(int32) :: boundary_id

        boundary_id = self%boundary_ids(bc_index)
    end function get_boundary_id_bc_manager

end module physics_governing_boundary_manager
