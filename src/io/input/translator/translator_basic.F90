submodule(io_input_translator) translator_basic
    implicit none
contains

    module subroutine execute_basic_properties(self, input, material_id, target_property, config)
        implicit none
        class(type_input_translator), intent(in) :: self
        class(type_input), intent(in) :: input
        integer(int32), intent(in) :: material_id
        type(type_constant_id), intent(in) :: target_property
        type(type_config_constitutive), intent(inout) :: config

        if (material_id < 1 .or. material_id > input%basic%num_materials) then
            error stop "Input Error: material_id is out of range."
        end if

        associate (material => input%basic%materials(material_id))
            select case (target_property%ID)
            case (CONSTITUTIVE_PROPERTIES%DENSITY%ID)
                config%material_id = material_id
                config%num_phases = material%phase
                if (allocated(material%density%value)) then
                    config%solid = material%density%value(1)
                    if (material%phase >= 2 .and. size(material%density%value) >= 2) then
                        config%water = material%density%value(2)
                    end if
                    if (material%phase >= 3 .and. size(material%density%value) >= 3) then
                        config%ice = material%density%value(3)
                    end if
                    if (material%phase >= 4 .and. size(material%density%value) >= 4) then
                        config%vapor = material%density%value(4)
                    end if
                end if
            case (CONSTITUTIVE_PROPERTIES%SPECIFIC_HEAT%ID)
                config%material_id = material_id
                config%num_phases = material%phase
                if (allocated(material%specific_heat%value)) then
                    config%solid = material%specific_heat%value(1)
                    if (material%phase >= 2 .and. size(material%specific_heat%value) >= 2) then
                        config%water = material%specific_heat%value(2)
                    end if
                    if (material%phase >= 3 .and. size(material%specific_heat%value) >= 3) then
                        config%ice = material%specific_heat%value(3)
                    end if
                    if (material%phase >= 4 .and. size(material%specific_heat%value) >= 4) then
                        config%vapor = material%specific_heat%value(4)
                    end if
                end if
            case (CONSTITUTIVE_PROPERTIES%VOLUMETRIC_HEAT_CAPACITY%ID)
                config%material_id = material_id
                config%num_phases = material%phase
                if (allocated(material%volumetric_heat_capacity%value)) then
                    config%solid = material%volumetric_heat_capacity%value(1)
                    if (material%phase >= 2 .and. size(material%volumetric_heat_capacity%value) >= 2) then
                        config%water = material%volumetric_heat_capacity%value(2)
                    end if
                    if (material%phase >= 3 .and. size(material%volumetric_heat_capacity%value) >= 3) then
                        config%ice = material%volumetric_heat_capacity%value(3)
                    end if
                    if (material%phase >= 4 .and. size(material%volumetric_heat_capacity%value) >= 4) then
                        config%vapor = material%volumetric_heat_capacity%value(4)
                        if (allocated(material%volumetric_heat_capacity%params)) then
                            call allocate_array(config%params, source=material%volumetric_heat_capacity%params)
                        end if
                    end if
                end if
            case (CONSTITUTIVE_PROPERTIES%THERMAL_CONDUCTIVITY%ID)
                config%material_id = material_id
                config%num_phases = material%phase
                if (allocated(material%thermal_conductivity%value)) then
                    config%solid = material%thermal_conductivity%value(1)
                    if (material%phase >= 2 .and. size(material%thermal_conductivity%value) >= 2) then
                        config%water = material%thermal_conductivity%value(2)
                    end if
                    if (material%phase >= 3 .and. size(material%thermal_conductivity%value) >= 3) then
                        config%ice = material%thermal_conductivity%value(3)
                    end if
                    if (material%phase >= 4 .and. size(material%thermal_conductivity%value) >= 4) then
                        config%vapor = material%thermal_conductivity%value(4)
                        if (allocated(material%thermal_conductivity%params)) then
                            call allocate_array(config%params, source=material%thermal_conductivity%params)
                        end if
                    end if
                end if

                if (material%thermal_conductivity%is_dispersed .and. &
                    allocated(material%thermal_conductivity%dispersivity)) then
                    call allocate_array(config%dispersivity, source=material%thermal_conductivity%dispersivity)
                end if
            case default
                error stop "execute_basic_properties: unsupported target property for type_config_constitutive."
            end select
        end associate

    end subroutine execute_basic_properties

    module subroutine execute_basic_swcc(self, input, material_id, config)
        implicit none
        class(type_input_translator), intent(in) :: self
        class(type_input), intent(in) :: input
        integer(int32), intent(in) :: material_id
        class(type_config_wrf), intent(inout) :: config

        if (material_id < 1 .or. material_id > input%basic%num_materials) then
            error stop "Input Error: material_id is out of range."
        end if

        select type (config)
        type is (type_config_hcf)

            call setup_wrf_basic(input, material_id, config)

            associate (hcf => input%basic%materials(material_id)%hydraulic_conductivity_model, &
                       swcc => input%basic%materials(material_id)%water_characteristic_curve)

                config%model = HCF_MODES%to_object(hcf%model_number)
                config%water_viscosity_model = HCF_VISCOSITY_TYPES%to_object(hcf%water_viscosity_model)

                config%k_sat = hcf%saturated_conductivity
                config%l = swcc%l
                config%omega = hcf%impedance_factor
                config%gain_factor = hcf%gain_factor

            end associate
        type is (type_config_wrf)

            call setup_wrf_basic(input, material_id, config)

        class default
            error stop "execute_basic_swcc: unsupported config type."

        end select

    end subroutine execute_basic_swcc

    subroutine setup_wrf_basic(input, material_id, config)
        implicit none
        class(type_input), intent(in) :: input
        integer(int32), intent(in) :: material_id
        class(type_config_wrf), intent(inout) :: config

        type(type_constant_id) :: L_unit
        ! WRF models operate on pressure head in METERS (fusion converts pore
        ! pressure p [Pa] to head via p/(rho_w g) before every WRF call). The
        ! input WRF parameters are therefore length-based and must be converted
        ! from the declared input length unit to meters: inverse-length params
        ! (van Genuchten alpha) are divided, length params (air-entry head) are
        ! multiplied. scale_length is meters per input unit.
        real(real64) :: scale_length

        associate (material => input%basic%materials(material_id)%water_characteristic_curve)

            call config%reset()

            config%material_id = material_id
            config%swcc_model = SWCC_MODELS%to_object(material%model_number)
            L_unit = PHYSICS_UNITS%to_object(material%unit)

            config%theta_s = material%theta_s
            config%theta_r = material%theta_r
            config%alpha1 = material%alpha1

            select case (config%swcc_model%ID)

            case (SWCC_MODELS%VG%ID)
                config%n1 = material%n1
                config%m1 = material%m1

            case (SWCC_MODELS%MVG%ID)
                config%n1 = material%n1
                config%m1 = material%m1
                config%h_crit = material%h_crit

            case (SWCC_MODELS%DURNER%ID)
                config%n1 = material%n1
                config%m1 = material%m1
                config%alpha2 = material%alpha2
                config%n2 = material%n2
                config%m2 = material%m2
                config%w1 = material%w1
                config%w2 = material%w2

            case (SWCC_MODELS%DVGCH%ID)
                config%n1 = material%n1
                config%m1 = material%m1
                config%n2 = material%n2
                config%m2 = material%m2
                config%w1 = material%w1
                config%w2 = material%w2

            end select

            select case (L_unit%ID)
            case (PHYSICS_UNITS%M%ID)
                scale_length = 1.0d0
            case (PHYSICS_UNITS%CM%ID)
                scale_length = 1.0d-2
            case (PHYSICS_UNITS%PA%ID)
                scale_length = 1.0d0 / (reference_water_density * gravity_acceleration)
            case default
                error stop "Input Error: unsupported SWCC parameter unit."
            end select

            select case (config%swcc_model%ID)

            case (SWCC_MODELS%BC%ID, SWCC_MODELS%KO%ID)
                config%alpha1 = config%alpha1 * scale_length
                config%h_crit = config%h_crit * scale_length
                config%alpha2 = config%alpha2 * scale_length

            case (SWCC_MODELS%VG%ID, SWCC_MODELS%DVGCH%ID)
                config%alpha1 = config%alpha1 / scale_length

            case (SWCC_MODELS%MVG%ID)
                config%alpha1 = config%alpha1 / scale_length
                config%h_crit = config%h_crit * scale_length

            case (SWCC_MODELS%DURNER%ID)
                config%alpha1 = config%alpha1 / scale_length
                config%alpha2 = config%alpha2 / scale_length

            end select

        end associate

    end subroutine setup_wrf_basic

    module subroutine execute_basic_gcc(self, input, material_id, config)
        implicit none
        class(type_input_translator), intent(in) :: self
        class(type_input), intent(in) :: input
        integer(int32), intent(in) :: material_id
        type(type_config_gcc), intent(inout) :: config

        if (material_id < 1 .or. material_id > input%basic%num_materials) then
            error stop "Input Error: material_id is out of range."
        end if

        !==================================================
        ! GCC model processing
        !==================================================

        associate (material => input%basic%materials(material_id)%phase_change%equilibrium_model, &
                   phase => input%basic%materials(material_id)%phase_change)

            call config%reset()

            config%material_id = material_id
            if (material%segregation .or. phase%segregation_potential > 0.0d0) then
                config%gcc_model = GCC_TYPES%SEGREGATION
            else
                config%gcc_model = GCC_TYPES%NON_SEGREGATION
            end if

            config%segregation_potential = phase%segregation_potential
            config%T_fringe_low = phase%T_fringe_low
            config%T_fringe_high = phase%T_fringe_high

        end associate

    end subroutine execute_basic_gcc

end submodule translator_basic
