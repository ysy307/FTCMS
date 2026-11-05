submodule(io_input_basic) input_basic_materials
    implicit none
    !!------------------------------------------------------------------------------------------------------------------------------
    ! JSON key names for materials
    !!------------------------------------------------------------------------------------------------------------------------------
    character(*), parameter :: materials = "materials"
    character(*), parameter :: id = "id"
    character(*), parameter :: name = "name"

    character(*), parameter :: phase = "phase"
    character(*), parameter :: is_frozen = "is_frozen"
    character(*), parameter :: is_dispersed = "is_dispersed"

    character(*), parameter :: value = "value"
    character(*), parameter :: params = "params"

    character(*), parameter :: density = "density"
    character(*), parameter :: specific_heat = "specific_heat"
    character(*), parameter :: thermal_conductivity = "thermal_conductivity"
    character(*), parameter :: volumetric_heat_capacity = "volumetric_heat_capacity"
    character(*), parameter :: dispersivity = "dispersivity"
    character(*), parameter :: phase_change = "phase_change"
    character(*), parameter :: latent_heat = "latent_heat"
    character(*), parameter :: fusion = "fusion"
    character(*), parameter :: freezing_temperature = "freezing_temperature"
    character(*), parameter :: water_characteristic_curve = "water_characteristic_curve"
    character(*), parameter :: model_number = "model_number"
    character(*), parameter :: model = "model"
    character(*), parameter :: theta_s = "theta_s"
    character(*), parameter :: theta_r = "theta_r"
    character(*), parameter :: n1 = "n1"
    character(*), parameter :: n2 = "n2"
    character(*), parameter :: m1 = "m1"
    character(*), parameter :: m2 = "m2"
    character(*), parameter :: alpha1 = "alpha1"
    character(*), parameter :: alpha2 = "alpha2"
    character(*), parameter :: w1 = "w1"
    character(*), parameter :: h_crit = "h_crit"
    character(*), parameter :: equilibrium_model = "equilibrium_model"
    character(*), parameter :: segregation = "segregation"
    character(*), parameter :: segregation_potential_key = "segregation_potential"
    character(*), parameter :: t_fringe_low_key = "T_fringe_low"
    character(*), parameter :: t_fringe_high_key = "T_fringe_high"
    character(*), parameter :: unit = "unit"
    character(*), parameter :: valid_units(3) = [character(len=4) :: "cm", "m", "Pa"]
    character(*), parameter :: hydraulic_conductivity_model = "hydraulic_conductivity_model"
    character(*), parameter :: hydraulic_storage = "hydraulic_storage"
    character(*), parameter :: specific_storage = "specific_storage"
    character(*), parameter :: saturated_conductivity = "saturated_conductivity"
    character(*), parameter :: l = "l"
    character(*), parameter :: impedance_factor = "impedance_factor"
    character(*), parameter :: gain_factor = "gain_factor"
    character(*), parameter :: water_viscosity_model = "water_viscosity_model"
    character(*), parameter :: water_retention_model = "water_retention_model"
contains
    module subroutine read_materials(self, json)
        implicit none
        class(type_input_basic), intent(inout) :: self
        type(json_file), intent(inout) :: json

        logical :: found
        integer(int32) :: i

        call json%info(materials, found=found, n_children=self%num_materials)
        call json%print_error_message(output_unit)
        if (.not. found .or. self%num_materials <= 0) then
            call raise_error(ERROR_CODES%VAR_INVALID, opt=materials)
        end if

        if (allocated(self%materials)) deallocate (self%materials)
        allocate (self%materials(self%num_materials))

        do i = 1, self%num_materials
            call read_materials_basic(self, json, i)

            call read_materials_physics(self, json, i)

            ! if (self%analysis_controls%is_active(get_physics_type(thermal))) then
            !     call read_materials_thermal(self, json, i)
            ! end if
            ! if (self%analysis_controls%is_active(get_physics_type(hydraulic))) then
            !     call read_materials_hydraulic(self, json, i)
            ! end if
            ! if (self%analysis_controls%is_active(get_physics_type(mechanical))) then
            !     ! Mechanical parameters can be added here in the future
            ! end if
        end do

    end subroutine read_materials

    subroutine read_materials_basic(self, json, i_material)
        !> Load the basic material parameters from the JSON file
        implicit none
        class(type_input_basic), intent(inout) :: self
        type(json_file), intent(inout) :: json !! JSON parser
        integer(int32), intent(in) :: i_material !! Material index

        character(256) :: buffer(2)

        buffer(1) = join([materials//"("//to_string(i_material)//")"])
        buffer(2) = id
        call get_json_value(json, join(buffer), self%materials(i_material)%ID, &
                            is_required=.true.)

        buffer(2) = name
        call get_json_value(json, join(buffer), self%materials(i_material)%name, &
                            is_required=.false., default_value="Material_"//to_string(i_material))

        ! if (self%analysis_controls%is_active(get_physics_type(thermal))) then
        buffer(2) = calculate_thermal
        call get_json_value(json, join(buffer), self%materials(i_material)%is_active(PHYSICS_TYPES%THERMAL%ID), &
                            is_required=.false., default_value=.false.)
        ! end if

        ! if (self%analysis_controls%is_active(PHYSICS_TYPES%HYDRAULIC%ID)) then
        buffer(2) = calculate_hydraulic
        call get_json_value(json, join(buffer), self%materials(i_material)%is_active(PHYSICS_TYPES%HYDRAULIC%ID), &
                            is_required=.false., default_value=.false.)
        ! end if

        ! if (self%analysis_controls%is_active(PHYSICS_TYPES%MECHANICAL%ID)) then
        buffer(2) = calculate_mechanical
        call get_json_value(json, join(buffer), self%materials(i_material)%is_active(PHYSICS_TYPES%MECHANICAL%ID), &
                            is_required=.false., default_value=.false.)
        ! end if

        buffer(2) = phase
        call get_json_value(json, join(buffer), self%materials(i_material)%phase, &
                            is_required=.true., valid_range=[1, 4])

    end subroutine read_materials_basic

    subroutine read_materials_physics(self, json, i_material)
        implicit none
        class(type_input_basic), intent(inout) :: self
        type(json_file), intent(inout) :: json !! JSON parser
        integer(int32), intent(in) :: i_material !! Material index

        character(256) :: buffer(5)

        ! Assigned here, not in the declaration: an initialiser would SAVE it across materials.
        buffer(:) = ""
        buffer(1) = materials//"("//to_string(i_material)//")"

        buffer(2) = density
        buffer(3) = value
        call get_json_value(json, join(buffer(1:4)), self%materials(i_material)%density%value, &
                            is_required=.true., valid_range=[0.0d0, huge(0.0d0)], array_size=self%materials(i_material)%phase)

        buffer(2) = specific_heat
        buffer(3) = value
        call get_json_value(json, join(buffer(1:4)), self%materials(i_material)%specific_heat%value, &
                            is_required=.true., valid_range=[0.0d0, huge(0.0d0)], array_size=self%materials(i_material)%phase)

        buffer(2) = volumetric_heat_capacity
        buffer(3) = value
        call get_json_value(json, join(buffer(1:4)), self%materials(i_material)%volumetric_heat_capacity%value, &
                            is_required=.true., valid_range=[0.0d0, huge(0.0d0)], array_size=self%materials(i_material)%phase)

        if (self%materials(i_material)%phase > 3) then
            buffer(3) = params
            call get_json_value(json, join(buffer(1:4)), self%materials(i_material)%volumetric_heat_capacity%params, &
                                is_required=.false.)
        end if

        buffer(2) = thermal_conductivity
        buffer(3) = value
        call get_json_value(json, join(buffer(1:4)), self%materials(i_material)%thermal_conductivity%value, &
                            is_required=.true., valid_range=[0.0d0, huge(0.0d0)], array_size=self%materials(i_material)%phase)

        buffer(3) = is_dispersed
        call get_json_value(json, join(buffer(1:4)), self%materials(i_material)%thermal_conductivity%is_dispersed, &
                            is_required=.false., default_value=.false.)

        if (self%materials(i_material)%thermal_conductivity%is_dispersed) then
            buffer(3) = dispersivity
            call get_json_value(json, join(buffer(1:4)), self%materials(i_material)%thermal_conductivity%dispersivity, &
                                is_required=.true., valid_range=[0.0d0, huge(0.0d0)], array_size=2)
        end if

        if (self%materials(i_material)%phase > 3) then
            buffer(3) = params
            call get_json_value(json, join(buffer(1:4)), self%materials(i_material)%thermal_conductivity%params, &
                                is_required=.false.)
        end if

        if (self%materials(i_material)%phase > 2) then
            buffer(2) = phase_change
            buffer(3) = equilibrium_model
            buffer(4) = segregation
            call get_json_value(json, join(buffer), self%materials(i_material)%phase_change%equilibrium_model%segregation, &
                                is_required=.true., default_value=.false.)

            buffer(3) = segregation_potential_key
            call get_json_value(json, join(buffer(1:3)), self%materials(i_material)%phase_change%segregation_potential, &
                                is_required=.false., default_value=0.0d0, valid_range=[0.0d0, huge(0.0d0)])

            buffer(3) = t_fringe_low_key
            call get_json_value(json, join(buffer(1:3)), self%materials(i_material)%phase_change%T_fringe_low, &
                                is_required=.false., default_value=-1.0d0)

            buffer(3) = t_fringe_high_key
            call get_json_value(json, join(buffer(1:3)), self%materials(i_material)%phase_change%T_fringe_high, &
                                is_required=.false., default_value=0.0d0)

            call read_materials_swcc(self, json, i_material)
            call read_materials_hydraulic_model(self, json, i_material)
            call read_materials_hydraulic_storage(self, json, i_material)
        end if

    end subroutine read_materials_physics

    subroutine read_materials_hydraulic_model(self, json, i_material)
        !> Load the hydraulic parameters from the JSON file
        implicit none
        class(type_input_basic), intent(inout) :: self
        type(json_file), intent(inout) :: json !! JSON parser
        integer(int32), intent(in) :: i_material !! Material index

        character(256) :: buffer(3)

        associate (hydraulic_conductivity => self%materials(i_material)%hydraulic_conductivity_model)
            ! Assigned here, not in the declaration: an initialiser would SAVE it across materials.
            buffer(:) = ""
            buffer(1) = materials//"("//to_string(i_material)//")"
            buffer(2) = hydraulic_conductivity_model

            buffer(3) = model_number
            call get_json_value(json, join(buffer(1:3)), hydraulic_conductivity%model_number, &
                                is_required=.true., valid_range=[1, 7])

            buffer(3) = saturated_conductivity
            call get_json_value(json, join(buffer(1:3)), hydraulic_conductivity%saturated_conductivity, &
                                is_required=.true., valid_range=[0.0d0, huge(0.0d0)])

            buffer(3) = impedance_factor
            call get_json_value(json, join(buffer(1:3)), hydraulic_conductivity%impedance_factor, &
                                is_required=.false., default_value=1.0d0, valid_range=[0.0d0, huge(0.0d0)])

            buffer(3) = water_viscosity_model
            call get_json_value(json, join(buffer(1:3)), hydraulic_conductivity%water_viscosity_model, &
                                is_required=.false., default_value=HCF_VISCOSITY_TYPES%EXPONENTIAL%ID, valid_range=[1, 2])

            buffer(3) = gain_factor
            call get_json_value(json, join(buffer(1:3)), hydraulic_conductivity%gain_factor, &
                                is_required=.false., default_value=1.0d0, valid_range=[0.0d0, huge(0.0d0)])
        end associate
    end subroutine read_materials_hydraulic_model

    subroutine read_materials_hydraulic_storage(self, json, i_material)
        implicit none
        class(type_input_basic), intent(inout) :: self
        type(json_file), intent(inout) :: json
        integer(int32), intent(in) :: i_material

        character(256) :: buffer(3)

        ! Assigned here, not in the declaration: an initialiser would SAVE it across materials.
        buffer(:) = ""
        buffer(1) = materials//"("//to_string(i_material)//")"
        buffer(2) = hydraulic_storage
        buffer(3) = specific_storage
        call get_json_value(json, join(buffer), self%materials(i_material)%hydraulic_storage%specific_storage, &
                            is_required=.false., default_value=0.0d0, valid_range=[0.0d0, huge(0.0d0)])
    end subroutine read_materials_hydraulic_storage

    subroutine read_materials_swcc(self, json, i_material)
        implicit none
        class(type_input_basic), intent(inout) :: self
        type(json_file), intent(inout) :: json !! JSON parser
        integer(int32), intent(in) :: i_material !! Material index

        character(256) :: buffer(3)
        logical :: found
        integer(int32) :: nul_pos

        character(:), allocatable :: tmp_strings

        ! Assigned here, not in the declaration: an initialiser would SAVE it across materials.
        buffer(:) = ""
        buffer(1) = materials//"("//to_string(i_material)//")"
        buffer(2) = water_characteristic_curve

        associate (swcc => self%materials(i_material)%water_characteristic_curve)
            buffer(3) = model
            call get_json_value(json, join(buffer(1:3)), tmp_strings, is_required=.true.)

            ! Normalize parser output to avoid hidden characters (e.g., NUL) causing model lookup failures.
            tmp_strings = strip(tmp_strings)
            nul_pos = index(tmp_strings, achar(0))
            if (nul_pos > 0) then
                if (nul_pos == 1) then
                    tmp_strings = ""
                else
                    tmp_strings = tmp_strings(:nul_pos - 1)
                end if
            end if

            swcc%model_number = SWCC_MODELS%to_id(tmp_strings)
            if (swcc%model_number < 0) then
                write (*, *) 'Error: Invalid SWCC model in input ', trim(tmp_strings)
                stop 1
            end if

            buffer(3) = unit
            call get_json_value(json, join(buffer(1:3)), tmp_strings, is_required=.true., valid_list=valid_units)
            swcc%unit = PHYSICS_UNITS%to_id(tmp_strings)

            buffer(3) = theta_s
            call get_json_value(json, join(buffer(1:3)), swcc%theta_s, is_required=.true., valid_range=[0.0d0, 1.0d0])
            buffer(3) = theta_r
            call get_json_value(json, join(buffer(1:3)), swcc%theta_r, is_required=.true., valid_range=[0.0d0, 1.0d0])

            buffer(3) = alpha1
            call get_json_value(json, join(buffer(1:3)), swcc%alpha1, is_required=.true.)

            buffer(3) = n1
            call get_json_value(json, join(buffer(1:3)), swcc%n1, is_required=.true.)

            buffer(3) = m1
            call get_json_value(json, join(buffer(1:3)), swcc%m1, found=found, is_required=.false., valid_range=[0.0d0, 1.0d0])
            if (.not. found) then
                swcc%m1 = 1.0d0 - 1.0d0 / swcc%n1
            end if

            select case (swcc%model_number)
            case (4)
                buffer(3) = h_crit
                call get_json_value(json, join(buffer(1:3)), swcc%h_crit, &
                                    is_required=.true., valid_range=[-huge(0.0d0), 0.0d0])
            case (5)
                buffer(3) = alpha2
                call get_json_value(json, join(buffer(1:3)), swcc%alpha2, is_required=.true.)
                buffer(3) = n2
                call get_json_value(json, join(buffer(1:3)), swcc%n2, is_required=.true.)

                buffer(3) = m2
                call get_json_value(json, join(buffer(1:3)), swcc%m2, found=found, is_required=.false., valid_range=[0.0d0, 1.0d0])
                if (.not. found) then
                    swcc%m2 = 1.0d0 - 1.0d0 / swcc%n2
                end if

                buffer(3) = w1
                call get_json_value(json, join(buffer(1:3)), swcc%w1, is_required=.true., valid_range=[0.0d0, 1.0d0])
                swcc%w2 = 1.0d0 - swcc%w1
            case (6)
                buffer(3) = alpha2
                call get_json_value(json, join(buffer(1:3)), swcc%alpha2, is_required=.true.)
                buffer(3) = n2
                call get_json_value(json, join(buffer(1:3)), swcc%n2, is_required=.true.)
                buffer(3) = m2
                call get_json_value(json, join(buffer(1:3)), swcc%m2, found=found, is_required=.false., valid_range=[0.0d0, 1.0d0])
                    if (.not. found) then
                        swcc%m2 = 1.0d0 - 1.0d0 / swcc%n2
                    end if
                buffer(3) = w1
                call get_json_value(json, join(buffer(1:3)), swcc%w1, is_required=.true., valid_range=[0.0d0, 1.0d0])
                swcc%w2 = 1.0d0 - swcc%w1
            end select

            buffer(3) = l
            call get_json_value(json, join(buffer(1:3)), swcc%l, is_required=.true., valid_range=[0.0d0, huge(0.0d0)])

        end associate

    end subroutine read_materials_swcc

    ! !---
    ! ! getter
    ! ! ---
    ! module subroutine get_density_info_material_settings(self, density_info)
    !     implicit none
    !     class(type_material_settings), intent(in) :: self
    !     type(type_physics_info), intent(inout) :: density_info

    !     call density_info%reset()

    !     density_info%num_phases = self%phase
    !     if (allocated(self%density%value)) then
    !         density_info%solid = self%density%value(1)
    !         if (self%phase >= 2 .and. size(self%density%value) >= 2) then
    !             density_info%water = self%density%value(2)
    !         end if
    !         if (self%phase >= 3 .and. size(self%density%value) >= 3) then
    !             density_info%ice = self%density%value(3)
    !         end if
    !         if (self%phase >= 4 .and. size(self%density%value) >= 4) then
    !             density_info%vapor = self%density%value(4)
    !         end if
    !     end if
    ! end subroutine get_density_info_material_settings

    ! module subroutine get_specific_heat_info_material_settings(self, specific_heat_info)
    !     implicit none
    !     class(type_material_settings), intent(in) :: self
    !     type(type_physics_info), intent(inout) :: specific_heat_info

    !     call specific_heat_info%reset()

    !     specific_heat_info%num_phases = self%phase
    !     if (allocated(self%specific_heat%value)) then
    !         specific_heat_info%solid = self%specific_heat%value(1)
    !         if (self%phase >= 2 .and. size(self%specific_heat%value) >= 2) then
    !             specific_heat_info%water = self%specific_heat%value(2)
    !         end if
    !         if (self%phase >= 3 .and. size(self%specific_heat%value) >= 3) then
    !             specific_heat_info%ice = self%specific_heat%value(3)
    !         end if
    !         if (self%phase >= 4 .and. size(self%specific_heat%value) >= 4) then
    !             specific_heat_info%vapor = self%specific_heat%value(4)
    !         end if
    !     end if
    ! end subroutine get_specific_heat_info_material_settings

    ! module subroutine get_volumetric_heat_capacity_info_material_settings(self, volumetric_heat_capacity_info)
    !     implicit none
    !     class(type_material_settings), intent(in) :: self
    !     type(type_physics_info), intent(inout) :: volumetric_heat_capacity_info

    !     call volumetric_heat_capacity_info%reset()

    !     volumetric_heat_capacity_info%num_phases = self%phase
    !     if (allocated(self%volumetric_heat_capacity%value)) then
    !         volumetric_heat_capacity_info%solid = self%volumetric_heat_capacity%value(1)
    !         if (self%phase >= 2 .and. size(self%volumetric_heat_capacity%value) >= 2) then
    !             volumetric_heat_capacity_info%water = self%volumetric_heat_capacity%value(2)
    !         end if
    !         if (self%phase >= 3 .and. size(self%volumetric_heat_capacity%value) >= 3) then
    !             volumetric_heat_capacity_info%ice = self%volumetric_heat_capacity%value(3)
    !         end if
    !         if (self%phase >= 4 .and. size(self%volumetric_heat_capacity%value) >= 4) then
    !             volumetric_heat_capacity_info%vapor = self%volumetric_heat_capacity%value(4)
    !             if (allocated(self%volumetric_heat_capacity%params)) then
    !                 call allocate_array(volumetric_heat_capacity_info%params, source=self%volumetric_heat_capacity%params)
    !             end if
    !         end if
    !     end if
    ! end subroutine get_volumetric_heat_capacity_info_material_settings

    ! module subroutine get_thermal_conductivity_info_material_settings(self, thermal_conductivity_info)
    !     implicit none
    !     class(type_material_settings), intent(in) :: self
    !     type(type_physics_info), intent(inout) :: thermal_conductivity_info

    !     call thermal_conductivity_info%reset()

    !     thermal_conductivity_info%num_phases = self%phase
    !     if (allocated(self%thermal_conductivity%value)) then
    !         thermal_conductivity_info%solid = self%thermal_conductivity%value(1)
    !         if (self%phase >= 2 .and. size(self%thermal_conductivity%value) >= 2) then
    !             thermal_conductivity_info%water = self%thermal_conductivity%value(2)
    !         end if
    !         if (self%phase >= 3 .and. size(self%thermal_conductivity%value) >= 3) then
    !             thermal_conductivity_info%ice = self%thermal_conductivity%value(3)
    !         end if
    !         if (self%phase >= 4 .and. size(self%thermal_conductivity%value) >= 4) then
    !             thermal_conductivity_info%vapor = self%thermal_conductivity%value(4)
    !             if (allocated(self%thermal_conductivity%params)) then
    !                 call allocate_array(thermal_conductivity_info%params, source=self%thermal_conductivity%params)
    !             end if
    !         end if
    !     end if

    !     if (self%thermal_conductivity%is_dispersed .and. allocated(self%thermal_conductivity%dispersivity)) then
    !         call allocate_array(thermal_conductivity_info%dispersivity, source=self%thermal_conductivity%dispersivity)
    !     end if
    ! end subroutine get_thermal_conductivity_info_material_settings

    ! module subroutine get_wrf_info_material_settings(self, wrf_info)
    !     implicit none
    !     class(type_material_settings), intent(in) :: self
    !     type(type_wrf_params), intent(inout) :: wrf_info

    !     call wrf_info%reset()

    !     wrf_info%model_number = self%water_characteristic_curve%model_number
    !     wrf_info%unit_id = self%water_characteristic_curve%unit
    !     wrf_info%theta_s = self%water_characteristic_curve%theta_s
    !     wrf_info%theta_r = self%water_characteristic_curve%theta_r
    !     wrf_info%alpha1 = self%water_characteristic_curve%alpha1
    !     wrf_info%n1 = self%water_characteristic_curve%n1
    !     wrf_info%m1 = self%water_characteristic_curve%m1
    !     wrf_info%alpha2 = self%water_characteristic_curve%alpha2
    !     wrf_info%n2 = self%water_characteristic_curve%n2
    !     wrf_info%m2 = self%water_characteristic_curve%m2
    !     wrf_info%w1 = self%water_characteristic_curve%w1
    !     wrf_info%w2 = self%water_characteristic_curve%w2
    !     wrf_info%h_crit = self%water_characteristic_curve%h_crit

    ! end subroutine get_wrf_info_material_settings

    ! module subroutine get_gcc_info_material_settings(self, gcc_info)
    !     implicit none
    !     class(type_material_settings), intent(in) :: self
    !     integer(int32), intent(inout) :: gcc_info

    !     if (self%phase > 2) then
    !         if (self%phase_change%equilibrium_model%segregation) then
    !             gcc_info = GCC_TYPES%SEGREGATION%ID
    !         else
    !             gcc_info = GCC_TYPES%NON_SEGREGATION%ID
    !         end if
    !     else
    !         gcc_info = -1
    !     end if

    ! end subroutine get_gcc_info_material_settings

    ! module subroutine get_hcf_info_material_settings(self, hcf_info)
    !     implicit none
    !     class(type_material_settings), intent(in) :: self
    !     type(type_hcf_params), intent(inout) :: hcf_info

    !     call hcf_info%reset()

    !     hcf_info%model_number = self%hydraulic_conductivity_model%model_number
    !     hcf_info%k_s = self%hydraulic_conductivity_model%saturated_conductivity
    !     hcf_info%omega = self%hydraulic_conductivity_model%impedance_factor
    !     hcf_info%water_viscosity_model = self%hydraulic_conductivity_model%water_viscosity_model
    !     hcf_info%gain_factor = self%hydraulic_conductivity_model%gain_factor

    !     hcf_info%unit_id = self%water_characteristic_curve%unit
    !     hcf_info%hcf_model_number = self%water_characteristic_curve%model_number
    !     hcf_info%theta_r = self%water_characteristic_curve%theta_r
    !     hcf_info%theta_s = self%water_characteristic_curve%theta_s
    !     hcf_info%alpha1 = self%water_characteristic_curve%alpha1
    !     hcf_info%n1 = self%water_characteristic_curve%n1
    !     hcf_info%m1 = self%water_characteristic_curve%m1
    !     hcf_info%h_crit = self%water_characteristic_curve%h_crit
    !     hcf_info%alpha2 = self%water_characteristic_curve%alpha2
    !     hcf_info%n2 = self%water_characteristic_curve%n2
    !     hcf_info%m2 = self%water_characteristic_curve%m2
    !     hcf_info%w1 = self%water_characteristic_curve%w1
    !     hcf_info%w2 = self%water_characteristic_curve%w2
    !     hcf_info%l = self%water_characteristic_curve%l

    ! end subroutine get_hcf_info_material_settings

!---------------------------------------------------------------------------------------------------------------------------
    ! Display Subroutines
    !---------------------------------------------------------------------------------------------------------------------------
    module subroutine display_material_settings(self)
        implicit none
        class(type_material_settings), intent(in) :: self

        ! --- 1. Basic Information ---
        call display_material_basic(self)

        ! --- 2. Thermal Properties ---
        !     The struct is flat, but display is grouped by physics type for readability
        if (self%is_active(PHYSICS_TYPES%THERMAL%ID)) then
            call display_material_thermal(self)
        end if

        ! --- 3. Hydraulic Properties ---
        if (self%is_active(PHYSICS_TYPES%HYDRAULIC%ID)) then
            call display_material_hydraulic(self)
        end if

        ! --- 4. Mechanical Properties ---
        if (self%is_active(PHYSICS_TYPES%MECHANICAL%ID)) then
            write (*, '(a)') "  Mechanical Properties: (Not implemented)"
        end if

    end subroutine display_material_settings

    subroutine display_material_basic(material)
        implicit none
        class(type_material_settings), intent(in) :: material

        write (*, '(a, i0)') "  ID                  : ", material%ID
        if (allocated(material%name)) then
            write (*, '(a, a)') "  Name                : ", trim(material%name)
        end if
        write (*, '(a, i0)') "  Phase Count         : ", material%phase
        ! Remove this line if is_frozen flag is absent from the struct; restore if present
        ! write (*, '(a, g0)') "  Is Frozen           : ", material%is_frozen

        write (*, '(a, L1)') "  Calculate Thermal   : ", material%is_active(PHYSICS_TYPES%THERMAL%ID)
        write (*, '(a, L1)') "  Calculate Hydraulic : ", material%is_active(PHYSICS_TYPES%HYDRAULIC%ID)
        write (*, '(a, L1)') "  Calculate Mechanical: ", material%is_active(PHYSICS_TYPES%MECHANICAL%ID)
    end subroutine display_material_basic

    subroutine display_material_thermal(material)
        implicit none
        class(type_material_settings), intent(in) :: material
        character(:), allocatable :: fmt

        write (*, '(a)') "  --- Thermal Properties ---"

        ! Density
        if (allocated(material%density%value)) then
            fmt = '(a, '//to_string(size(material%density%value))//'(es12.4e2, 2x))'
            write (*, fmt) "    Density             : ", material%density%value
        end if

        ! Specific Heat
        if (allocated(material%specific_heat%value)) then
            fmt = '(a, '//to_string(size(material%specific_heat%value))//'(es12.4e2, 2x))'
            write (*, fmt) "    Specific Heat       : ", material%specific_heat%value
        end if

        ! Volumetric Heat Capacity
        if (allocated(material%volumetric_heat_capacity%value)) then
            fmt = '(a, '//to_string(size(material%volumetric_heat_capacity%value))//'(es12.4e2, 2x))'
            write (*, fmt) "    Vol. Heat Cap.      : ", material%volumetric_heat_capacity%value

            if (allocated(material%volumetric_heat_capacity%params)) then
                fmt = '(a, '//to_string(size(material%volumetric_heat_capacity%params))//'(es12.4e2, 2x))'
                write (*, fmt) "    (Params)            : ", material%volumetric_heat_capacity%params
            end if
        end if

        ! Thermal Conductivity
        if (allocated(material%thermal_conductivity%value)) then
            fmt = '(a, '//to_string(size(material%thermal_conductivity%value))//'(es12.4e2, 2x))'
            write (*, fmt) "    Thermal Conductivity: ", material%thermal_conductivity%value
        end if

        if (allocated(material%thermal_conductivity%params)) then
            fmt = '(a, '//to_string(size(material%thermal_conductivity%params))//'(es12.4e2, 2x))'
            write (*, fmt) "    (Params)            : ", material%thermal_conductivity%params
        end if

        write (*, '(a, L1)') "    Is Dispersed        : ", material%thermal_conductivity%is_dispersed

        if (material%thermal_conductivity%is_dispersed .and. allocated(material%thermal_conductivity%dispersivity)) then
            write (*, '(a, es12.4e2, " / ", es12.4e2)') "    Dispersivity (L/T)  : ", &
                material%thermal_conductivity%dispersivity(1), material%thermal_conductivity%dispersivity(2)
        end if

        ! Phase Change
        write (*, '(a)') "    --- Phase Change ---"
        ! Adjust access path to equilibrium_model to match loading code
        write (*, '(a, L1)') "    Segregation         : ", material%phase_change%equilibrium_model%segregation

    end subroutine display_material_thermal

    subroutine display_material_hydraulic(material)
        implicit none
        class(type_material_settings), intent(in) :: material

        write (*, '(a)') "  --- Hydraulic Properties ---"

        ! Hydraulic Conductivity
        write (*, '(a, i0)') "    HC Model Number     : ", material%hydraulic_conductivity_model%model_number
        write (*, '(a, es12.4e2)') "    Saturated K         : ", material%hydraulic_conductivity_model%saturated_conductivity
        write (*, '(a, es12.4e2)') "    Impedance Factor    : ", material%hydraulic_conductivity_model%impedance_factor
        write (*, '(a, i0)') "    Viscosity Model #   : ", material%hydraulic_conductivity_model%water_viscosity_model
        write (*, '(a, es12.4e2)') "    Gain Factor         : ", material%hydraulic_conductivity_model%gain_factor
        write (*, '(a, es12.4e2)') "    Specific Storage    : ", material%hydraulic_storage%specific_storage

        ! Water Retention Model (WCC)
        call display_wcc(material%water_characteristic_curve)

    end subroutine display_material_hydraulic

    subroutine display_wcc(wcc)
        implicit none
        class(type_materials_water_characteristic_curve), intent(in) :: wcc

        type(type_constant_id) :: constant

        ! Model Name (String)
        constant = SWCC_MODELS%to_object(wcc%model_number)
        write (*, '(a, a)') "    [WCC] Model         : ", strip(constant%NAME)
        constant = PHYSICS_UNITS%to_object(wcc%unit)
        write (*, '(a, a)') "      Unit              : ", strip(constant%NAME)

        write (*, '(a, es12.4e2)') "      theta_s           : ", wcc%theta_s
        write (*, '(a, es12.4e2)') "      theta_r           : ", wcc%theta_r
        write (*, '(a, es12.4e2)') "      alpha1            : ", wcc%alpha1
        write (*, '(a, es12.4e2)') "      n1                : ", wcc%n1
        write (*, '(a, es12.4e2)') "      l                 : ", wcc%l

        ! The following values may be unset; conditionally displaying them could be
        ! considered, but here we simply print them (expecting 0.0 for unset values)
        write (*, '(a, es12.4e2)') "      h_crit            : ", wcc%h_crit
        write (*, '(a, es12.4e2)') "      alpha2            : ", wcc%alpha2
        write (*, '(a, es12.4e2)') "      n2                : ", wcc%n2
        write (*, '(a, es12.4e2)') "      w1                : ", wcc%w1

    end subroutine display_wcc

end submodule input_basic_materials
