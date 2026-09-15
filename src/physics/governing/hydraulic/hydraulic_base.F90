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

        call self%physics%initialize(active_region_ids, &
                                     configs_density=density_info, &
                                     configs_wrf=wrf_model_info, &
                                     configs_hcf=hcf_model_info, &
                                     configs_gcc=gcc_model_info)

        if (allocated(self%iteration_capacity_bound)) deallocate (self%iteration_capacity_bound)
        if (allocated(self%specific_storage)) deallocate (self%specific_storage)
        if (num_active_materials > 0) then
            allocate (self%iteration_capacity_bound(maxval(active_region_ids)), source=0.0d0)
            allocate (self%specific_storage(maxval(active_region_ids)), source=0.0d0)
            do i = 1, num_active_materials
                call self%physics%calc_lscheme_capacity(active_region_ids(i), &
                                                        self%iteration_capacity_bound(active_region_ids(i)))
                material_idx = 0
                do j = 1, num_materials
                    if (input%basic%materials(j)%id == active_region_ids(i)) then
                        material_idx = j
                        exit
                    end if
                end do
                if (material_idx == 0 .and. num_materials == 1) material_idx = 1
                if (material_idx > 0) then
                    self%specific_storage(active_region_ids(i)) = &
                        input%basic%materials(material_idx)%hydraulic_storage%specific_storage
                end if
            end do
        end if

        self%computation_type = input%basic%simulation_settings%calculate_type
        self%computation_dimension = input%basic%simulation_settings%calculate_dimension
        self%enable_vapor_transport = input%basic%analysis_controls%enable_vapor_transport

    end subroutine initialize_type_hydraulic

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
