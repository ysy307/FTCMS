submodule(boundary_strategy) strategy_base
    implicit none
contains

    module pure elemental subroutine initialize_type_bc_result(self)
        implicit none
        class(type_bc_result), intent(inout) :: self

        self%is_dirichlet = .false.
        self%prescribed_value = 0.0d0
        self%flux_value = 0.0d0
        self%flux_derivative = 0.0d0
    end subroutine initialize_type_bc_result

    module subroutine initialize_abst_bc(self, physics_type, config)
        implicit none
        class(abst_bc), intent(inout), target :: self
        type(type_constant_id), intent(in) :: physics_type
        type(type_config_bc), intent(in) :: config

        self%physics_type = physics_type
        self%bc_kind = config%bc_kind

        self%provider = create_bc_data_provider(config)
    end subroutine initialize_abst_bc

    module subroutine destroy_abst_bc(self)
        implicit none
        class(abst_bc), intent(inout) :: self

        if (allocated(self%provider)) then
            call self%provider%destroy()
            deallocate (self%provider)
        end if
    end subroutine destroy_abst_bc

    module subroutine update_provider_data_abst_bc(self, new_values)
        implicit none
        class(abst_bc), intent(inout) :: self
        real(real64), intent(in) :: new_values(:)

        if (allocated(self%provider)) call self%provider%update_data(new_values)
    end subroutine update_provider_data_abst_bc

    !> The provider's data at `current_time`, without the state-dependent part.
    module subroutine get_data_abst_bc(self, current_time, values)
        implicit none
        class(abst_bc), intent(inout) :: self
        real(real64), intent(in) :: current_time
        real(real64), intent(inout) :: values(BC_DATA_SIZE)

        values = 0.0d0
        if (allocated(self%provider)) call self%provider%get_data(current_time, values)
    end subroutine get_data_abst_bc

end submodule strategy_base
