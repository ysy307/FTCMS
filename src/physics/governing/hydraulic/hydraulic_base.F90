submodule(physics_governing_hydraulic) hydraulic_base
    implicit none
contains
    module subroutine initialize_type_hydraulic(self, input, active_region_ids)
        implicit none
        class(type_hydraulic), intent(inout) :: self
        type(type_input), intent(in) :: input
        integer(int32), intent(in) :: active_region_ids(:)

        integer(int32) :: num_materials
        integer(int32) :: num_active_materials
        integer(int32) :: i, j
        integer(int32) :: material_idx

        type(type_config_constitutive), allocatable :: density_info(:)
        type(type_config_wrf), allocatable :: wrf_model_info(:)
        type(type_config_hcf), allocatable :: hcf_model_info(:)
        type(type_config_gcc), allocatable :: gcc_model_info(:)
        real(real64), allocatable :: storage_info(:)

        ! Initialize thermal physics components

        num_materials = input%basic%num_materials
        num_active_materials = size(active_region_ids)

        allocate (density_info(num_active_materials))
        allocate (wrf_model_info(num_active_materials))
        allocate (hcf_model_info(num_active_materials))
        allocate (gcc_model_info(num_active_materials))

        do j = 1, num_active_materials
            material_idx = 0
            do i = 1, num_materials
                if (input%basic%materials(i)%id == active_region_ids(j)) then
                    material_idx = i
                    exit
                end if
            end do

            if (material_idx == 0) then
                if (num_materials == 1) then
                    material_idx = 1
                else
                    write (*, *) 'Error: No material definition matched active region id ', active_region_ids(j)
                    stop 1
                end if
            end if

            call density_info(j)%reset()

            call input_translator%execute(input, material_idx, CONSTITUTIVE_PROPERTIES%DENSITY, &
                                          density_info(j))
            call input_translator%execute(input, material_idx, wrf_model_info(j))
            call input_translator%execute(input, material_idx, hcf_model_info(j))
            call input_translator%execute(input, material_idx, gcc_model_info(j))
        end do

        ! The specific storage of each active region, in the order the
        ! configurations are given, so that construction below is the same for a
        ! project and for a caller that holds the configurations already.
        allocate (storage_info(num_active_materials), source=0.0d0)
        do j = 1, num_active_materials
            material_idx = 0
            do i = 1, num_materials
                if (input%basic%materials(i)%id == active_region_ids(j)) then
                    material_idx = i
                    exit
                end if
            end do
            if (material_idx == 0 .and. num_materials == 1) material_idx = 1
            if (material_idx > 0) then
                storage_info(j) = input%basic%materials(material_idx)%hydraulic_storage%specific_storage
            end if
        end do

        call self%initialize_from_configs(active_region_ids, density_info, wrf_model_info, &
                                          hcf_model_info, gcc_model_info, storage_info, &
                                          input%basic%simulation_settings%calculate_type, &
                                          input%basic%simulation_settings%calculate_dimension, &
                                          input%basic%analysis_controls%enable_vapor_transport)

    end subroutine initialize_type_hydraulic

    module subroutine initialize_type_hydraulic_from_configs(self, active_region_ids, &
                                                             configs_density, configs_wrf, configs_hcf, configs_gcc, &
                                                             specific_storage, computation_type, computation_dimension, &
                                                             enable_vapor_transport)
        implicit none
        class(type_hydraulic), intent(inout) :: self
        integer(int32), intent(in) :: active_region_ids(:)
        type(type_config_constitutive), intent(in) :: configs_density(:)
        type(type_config_wrf), intent(in) :: configs_wrf(:)
        type(type_config_hcf), intent(in) :: configs_hcf(:)
        type(type_config_gcc), intent(in) :: configs_gcc(:)
        real(real64), intent(in) :: specific_storage(:)
        integer(int32), intent(in) :: computation_type
        integer(int32), intent(in) :: computation_dimension
        logical, intent(in) :: enable_vapor_transport

        integer(int32) :: num_active_materials
        integer(int32) :: i

        num_active_materials = size(active_region_ids)
        if (size(configs_density) /= num_active_materials .or. &
            size(configs_wrf) /= num_active_materials .or. &
            size(configs_hcf) /= num_active_materials .or. &
            size(configs_gcc) /= num_active_materials .or. &
            size(specific_storage) /= num_active_materials) then
            error stop "initialize_type_hydraulic_from_configs: one configuration per active region is required"
        end if

        call self%physics%initialize(active_region_ids, &
                                     configs_density=configs_density, &
                                     configs_wrf=configs_wrf, &
                                     configs_hcf=configs_hcf, &
                                     configs_gcc=configs_gcc)

        if (allocated(self%iteration_capacity_bound)) deallocate (self%iteration_capacity_bound)
        if (allocated(self%specific_storage)) deallocate (self%specific_storage)
        if (num_active_materials > 0) then
            allocate (self%iteration_capacity_bound(maxval(active_region_ids)), source=0.0d0)
            allocate (self%specific_storage(maxval(active_region_ids)), source=0.0d0)
            do i = 1, num_active_materials
                call self%physics%calc_lscheme_capacity(active_region_ids(i), &
                                                        self%iteration_capacity_bound(active_region_ids(i)))
                self%specific_storage(active_region_ids(i)) = specific_storage(i)
            end do
        end if

        self%computation_type = computation_type
        self%computation_dimension = computation_dimension
        self%enable_vapor_transport = enable_vapor_transport

    end subroutine initialize_type_hydraulic_from_configs

    module pure function is_vapor_transport_enabled_hydraulic(self) result(enabled)
        implicit none
        class(type_hydraulic), intent(in) :: self
        logical :: enabled

        enabled = self%enable_vapor_transport
    end function is_vapor_transport_enabled_hydraulic

    module subroutine set_disable_capacity_regularization_hydraulic(self, disable)
        implicit none
        class(type_hydraulic), intent(inout) :: self
        logical, intent(in) :: disable

        self%disable_capacity_regularization = disable
    end subroutine set_disable_capacity_regularization_hydraulic

    module pure function get_capacity_regularization_hydraulic(self, material_id, nonlinear_iteration) &
        result(value)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        integer(int32), intent(in) :: nonlinear_iteration
        real(real64) :: value

        value = 0.0d0
        if (allocated(self%iteration_capacity_bound)) then
            if (material_id >= 1 .and. material_id <= size(self%iteration_capacity_bound)) then
                value = self%iteration_capacity_bound(material_id) / real(max(1, nonlinear_iteration), real64)
            end if
        end if
    end function get_capacity_regularization_hydraulic

end submodule hydraulic_base
