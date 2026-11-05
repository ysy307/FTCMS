submodule(boundary_data_provider) data_provider_table
    implicit none
contains

    module subroutine initialize_type_bc_data_table(self, config)
        implicit none
        class(type_bc_data_table), intent(inout) :: self
        type(type_config_bc), intent(in) :: config

        call allocate_array(self%time_points, source=config%time_points)
        call allocate_array(self%table_values, source=config%values)

        self%current_idx = 0

        self%data_kind = BC_DATA_PROVIDERS%TABLE
    end subroutine initialize_type_bc_data_table

    module subroutine destroy_type_bc_data_table(self)
        implicit none
        class(type_bc_data_table), intent(inout) :: self

        call deallocate_array(self%time_points)
        call deallocate_array(self%table_values)
        self%data_kind = type_constant_id("", "", -1)
    end subroutine destroy_type_bc_data_table

    module subroutine update_data_table(self, new_values)
        implicit none
        class(type_bc_data_table), intent(inout) :: self
        real(real64), intent(in) :: new_values(:)

        return
    end subroutine update_data_table

    module pure subroutine get_data_bc_data_table(self, current_time, values)
        implicit none
        class(type_bc_data_table), intent(inout) :: self
        real(real64), intent(in) :: current_time
        real(real64), intent(inout) :: values(BC_DATA_SIZE)

        real(real64) :: coef
        integer(int32) :: idx, num_vars

        values(:) = 0.0d0
        if (.not. allocated(self%time_points) .or. .not. allocated(self%table_values)) return

        call self%calc_time_coefficient(current_time, coef, idx)

        num_vars = min(BC_DATA_SIZE, size(self%table_values, 1))

        if (idx < size(self%time_points)) then
            values(1:num_vars) = self%table_values(1:num_vars, idx) + &
                                 coef * (self%table_values(1:num_vars, idx + 1) - self%table_values(1:num_vars, idx))
        else
            values(1:num_vars) = self%table_values(1:num_vars, size(self%time_points))
        end if
    end subroutine get_data_bc_data_table

    module pure subroutine calc_time_coefficient_bc_data_table(self, current_time, coef, idx)
        implicit none
        class(type_bc_data_table), intent(inout) :: self
        real(real64), intent(in) :: current_time
        real(real64), intent(inout) :: coef
        integer(int32), intent(inout) :: idx

        integer(int32) :: n

        if (.not. allocated(self%time_points)) then
            coef = 0.0d0
            idx = 1
            return
        end if

        n = size(self%time_points)
        if (n < 2 .or. current_time <= self%time_points(1)) then
            coef = 0.0d0
            idx = 1
            self%current_idx = 1
            return
        end if
        if (current_time >= self%time_points(n)) then
            coef = 0.0d0
            idx = n
            self%current_idx = n
            return
        end if

        ! Use cached index for forward marching
        idx = max(1, self%current_idx - 1)

        ! Linear search forward (assuming monotonic time advancement)
        do while (idx < n .and. current_time >= self%time_points(idx + 1))
            idx = idx + 1
        end do

        ! Fallback if time goes backward
        if (current_time < self%time_points(idx)) then
            idx = 1
            do while (idx < n .and. current_time >= self%time_points(idx + 1))
                idx = idx + 1
            end do
        end if

        self%current_idx = idx
        coef = (current_time - self%time_points(idx)) / &
               (self%time_points(idx + 1) - self%time_points(idx))
    end subroutine calc_time_coefficient_bc_data_table

end submodule data_provider_table
