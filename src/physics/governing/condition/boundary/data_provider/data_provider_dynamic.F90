submodule(boundary_data_provider) data_provider_dynamic
    implicit none
contains

    module subroutine initialize_type_bc_data_dynamic(self, config)
        implicit none
        class(type_bc_data_dynamic), intent(inout) :: self
        type(type_config_bc), intent(in) :: config

        self%data_kind = BC_DATA_PROVIDERS%DYNAMIC
    end subroutine initialize_type_bc_data_dynamic

    module subroutine destroy_type_bc_data_dynamic(self)
        implicit none
        class(type_bc_data_dynamic), intent(inout) :: self

        call deallocate_array(self%current_buffer)
        self%data_kind = type_constant_id("", "", -1)
    end subroutine destroy_type_bc_data_dynamic

    module subroutine update_data_bc_data_dynamic(self, new_values)
        implicit none
        class(type_bc_data_dynamic), intent(inout) :: self
        real(real64), intent(in) :: new_values(:)

        self%current_buffer = new_values
    end subroutine update_data_bc_data_dynamic

    module pure subroutine get_data_bc_data_dynamic(self, current_time, values)
        implicit none
        class(type_bc_data_dynamic), intent(inout) :: self
        real(real64), intent(in) :: current_time
        real(real64), intent(inout) :: values(BC_DATA_SIZE)

        values = 0.0d0
        values(1:min(BC_DATA_SIZE, size(self%current_buffer))) = &
            self%current_buffer(1:min(BC_DATA_SIZE, size(self%current_buffer)))
    end subroutine get_data_bc_data_dynamic

end submodule data_provider_dynamic
