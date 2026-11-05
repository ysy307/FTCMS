submodule(boundary_data_provider) data_provider_constant
    implicit none
contains

    module subroutine initialize_type_bc_data_constant(self, config)
        implicit none
        class(type_bc_data_constant), intent(inout) :: self
        type(type_config_bc), intent(in) :: config

        call allocate_array(self%constant_values, config%num_variables)
        self%constant_values = config%values(1:config%num_variables, 1) ! Assuming the first time
        self%data_size = config%num_variables

        self%data_kind = BC_DATA_PROVIDERS%CONSTANT
    end subroutine initialize_type_bc_data_constant

    module subroutine update_data_constant(self, new_values)
        implicit none
        class(type_bc_data_constant), intent(inout) :: self
        real(real64), intent(in) :: new_values(:)

        integer(int32) :: n
        n = min(size(new_values), size(self%constant_values))
        if (n > 0) self%constant_values(1:n) = new_values(1:n)
    end subroutine update_data_constant

    module subroutine destroy_type_bc_data_constant(self)
        implicit none
        class(type_bc_data_constant), intent(inout) :: self

        call deallocate_array(self%constant_values)

        self%data_kind = type_constant_id("", "", -1)
    end subroutine destroy_type_bc_data_constant

    module pure subroutine get_data_bc_data_constant(self, current_time, values)
        implicit none
        class(type_bc_data_constant), intent(inout) :: self
        real(real64), intent(in) :: current_time
        real(real64), intent(inout) :: values(BC_DATA_SIZE)

        values = 0.0d0
        values(1:self%data_size) = self%constant_values(1:self%data_size)

    end subroutine get_data_bc_data_constant

end submodule data_provider_constant
