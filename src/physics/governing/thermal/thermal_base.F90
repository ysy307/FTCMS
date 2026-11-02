submodule(physics_governing_thermal) thermal_base
    implicit none
contains
    module subroutine initialize_type_thermal(self, input, active_region_ids)
        implicit none
        class(type_thermal), intent(inout) :: self
        type(type_input), intent(in) :: input
        integer(int32), intent(in) :: active_region_ids(:)

        integer(int32) :: num_materials
        integer(int32) :: num_active_materials
        integer(int32) :: i, j
        integer(int32) :: material_idx

        type(type_config_constitutive), allocatable :: density_info(:)
        type(type_config_constitutive), allocatable :: specific_heat_info(:)
        type(type_config_constitutive), allocatable :: heat_capacity_info(:)
        type(type_config_constitutive), allocatable :: thermal_conductivity_info(:)
        type(type_config_wrf), allocatable :: wrf_model_info(:)
        type(type_config_hcf), allocatable :: hcf_model_info(:)
        type(type_config_gcc), allocatable :: gcc_model_info(:)

        ! Initialize thermal physics components

        num_materials = input%basic%num_materials
        num_active_materials = size(active_region_ids)

        allocate (density_info(num_active_materials))
        allocate (specific_heat_info(num_active_materials))
        allocate (heat_capacity_info(num_active_materials))
        allocate (thermal_conductivity_info(num_active_materials))
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
            call specific_heat_info(j)%reset()
            call heat_capacity_info(j)%reset()
            call thermal_conductivity_info(j)%reset()

            call input_translator%execute(input, material_idx, CONSTITUTIVE_PROPERTIES%DENSITY, &
                                          density_info(j))
            call input_translator%execute(input, material_idx, CONSTITUTIVE_PROPERTIES%SPECIFIC_HEAT, &
                                          specific_heat_info(j))
            call input_translator%execute(input, material_idx, CONSTITUTIVE_PROPERTIES%VOLUMETRIC_HEAT_CAPACITY, &
                                          heat_capacity_info(j))
            call input_translator%execute(input, material_idx, CONSTITUTIVE_PROPERTIES%THERMAL_CONDUCTIVITY, &
                                          thermal_conductivity_info(j))

            call input_translator%execute(input, material_idx, wrf_model_info(j))
            call input_translator%execute(input, material_idx, hcf_model_info(j))
            call input_translator%execute(input, material_idx, gcc_model_info(j))
        end do

        call self%initialize_from_configs(active_region_ids, density_info, &
                                          specific_heat_info, heat_capacity_info, thermal_conductivity_info, &
                                          wrf_model_info, hcf_model_info, gcc_model_info, &
                                          input%basic%simulation_settings%calculate_type, &
                                          input%basic%simulation_settings%calculate_dimension, &
                                          input%basic%analysis_controls%enable_vapor_transport)

    end subroutine initialize_type_thermal

    module subroutine initialize_type_thermal_from_configs(self, active_region_ids, &
                                                           configs_density, configs_specific_heat, &
                                                           configs_heat_capacity, configs_thermal_conductivity, &
                                                           configs_wrf, configs_hcf, configs_gcc, &
                                                           computation_type, computation_dimension, &
                                                           enable_vapor_transport)
        implicit none
        class(type_thermal), intent(inout) :: self
        integer(int32), intent(in) :: active_region_ids(:)
        type(type_config_constitutive), intent(in) :: configs_density(:)
        type(type_config_constitutive), intent(in) :: configs_specific_heat(:)
        type(type_config_constitutive), intent(in) :: configs_heat_capacity(:)
        type(type_config_constitutive), intent(in) :: configs_thermal_conductivity(:)
        type(type_config_wrf), intent(in) :: configs_wrf(:)
        type(type_config_hcf), intent(in) :: configs_hcf(:)
        type(type_config_gcc), intent(in) :: configs_gcc(:)
        integer(int32), intent(in) :: computation_type
        integer(int32), intent(in) :: computation_dimension
        logical, intent(in) :: enable_vapor_transport

        integer(int32) :: num_active_materials

        num_active_materials = size(active_region_ids)
        if (size(configs_density) /= num_active_materials .or. &
            size(configs_specific_heat) /= num_active_materials .or. &
            size(configs_heat_capacity) /= num_active_materials .or. &
            size(configs_thermal_conductivity) /= num_active_materials .or. &
            size(configs_wrf) /= num_active_materials .or. &
            size(configs_hcf) /= num_active_materials .or. &
            size(configs_gcc) /= num_active_materials) then
            error stop "initialize_type_thermal_from_configs: one configuration per active region is required"
        end if

        call self%physics%initialize(active_region_ids, configs_density, &
                                     configs_specific_heat, configs_heat_capacity, configs_thermal_conductivity, &
                                     configs_wrf, configs_hcf, configs_gcc)

        self%computation_type = computation_type
        self%computation_dimension = computation_dimension
        self%enable_vapor_transport = enable_vapor_transport

    end subroutine initialize_type_thermal_from_configs

    module subroutine destroy_type_thermal(self)
        implicit none
        class(type_thermal), intent(inout) :: self

        self%computation_dimension = 0
        self%computation_type = 0

    end subroutine destroy_type_thermal
end submodule thermal_base
