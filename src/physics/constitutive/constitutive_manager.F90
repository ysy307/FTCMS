module physics_constitutive_manager
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: iapws, only:type_iapws97, type_iapws06
    use :: module_core
    use :: constitutive_constants, only:tangent_pressure, tangent_temperature
    use :: physics_constitutive_base, only:type_iapws_wrapper
    use :: constitutive_materials_manager, only:type_material_manager
    use :: constitutive_models_manager, only:type_models_manager
    implicit none
    private
    public :: type_constitutive_manager

    !-------------------------------------------------------------------------------------------------------------------------------
    ! Main Derived Type with Generic Type-Bound Procedures
    !-------------------------------------------------------------------------------------------------------------------------------
    type :: type_constitutive_manager
        private
        type(type_iapws97) :: water
        type(type_iapws06) :: ice
        type(type_iapws_wrapper) :: iapws_wrapper
        type(type_material_manager), allocatable :: materials(:)
        type(type_models_manager), allocatable :: models(:)
        integer(int32), allocatable :: materials_id_map(:)
        integer(int32) :: num_materials
    contains
        procedure, public :: initialize => initialize_constitutive_manager
        procedure, private :: set_map => set_map_materials

        procedure, public :: calc_density
        procedure, public :: get_density_solid
        procedure, public :: calc_density_water
        procedure, public :: calc_density_ice
        procedure, public :: calc_density_vapor
        procedure, public :: calc_density_vapor_saturation
        procedure, public :: calc_density_water_derivatives
        procedure, public :: calc_density_ice_derivatives
        procedure, public :: calc_density_vapor_derivatives
        procedure, public :: calc_specific_heat
        procedure, public :: get_specific_heat_solid
        procedure, public :: calc_specific_heat_water
        procedure, public :: calc_specific_heat_ice
        procedure, public :: calc_specific_heat_vapor
        procedure, private :: calc_thermal_conductivity_nondispersivity
        procedure, private :: calc_thermal_conductivity_dispersivity
        generic, public :: calc_thermal_conductivity &
            => calc_thermal_conductivity_nondispersivity, calc_thermal_conductivity_dispersivity
        procedure, public :: calc_vol_heat_capacity
        procedure, public :: calc_latent_heat_fusion
        procedure, public :: calc_latent_heat_vaporization
        procedure, public :: calc_saturation_pressure
        procedure, public :: calc_pressure_ice_water_derivative

        procedure, public :: update_water_phases
        procedure, public :: project_ice_content
        procedure, public :: calc_conserved_target
        procedure, public :: solve_local_conserved_equilibrium
        procedure, public :: calc_Kflh
        procedure, public :: calc_KlT
        procedure, public :: calc_cryo_head_dT
        procedure, public :: calc_freezing_level_set
        procedure, public :: calc_Kvh
        procedure, public :: calc_KvT
        procedure, public :: calc_cryo_suction_deriv_T
        procedure, public :: calc_lscheme_capacity
        procedure, public :: calc_suction_weights
        procedure, public :: calc_segregation_sink
        procedure, public :: calc_effective_segregation_sink
        procedure, public :: is_segregation_active
        procedure, public :: has_cryo_transport

    end type type_constitutive_manager

contains

    !-------------------------------------------------------------------------------------------------------------------------------
    ! Initialization
    !-------------------------------------------------------------------------------------------------------------------------------
    subroutine initialize_constitutive_manager(self, unique_material_ids, configs_density, &
                                               configs_specific_heat, configs_heat_capacity, configs_thermal_conductivity, &
                                               configs_wrf, configs_hcf, configs_gcc)
        implicit none
        class(type_constitutive_manager), intent(inout) :: self
        integer(int32), intent(in) :: unique_material_ids(:)
        type(type_config_constitutive), intent(in), optional :: configs_density(:)
        type(type_config_constitutive), intent(in), optional :: configs_specific_heat(:)
        type(type_config_constitutive), intent(in), optional :: configs_heat_capacity(:)
        type(type_config_constitutive), intent(in), optional :: configs_thermal_conductivity(:)
        type(type_config_wrf), intent(in), optional :: configs_wrf(:)
        type(type_config_hcf), intent(in), optional :: configs_hcf(:)
        type(type_config_gcc), intent(in), optional :: configs_gcc(:)

        integer(int32) :: stat
        integer(int32) :: model_idx, current_material_id

        call self%water%initialize()
        call self%ice%initialize()

        call self%set_map(unique_material_ids)
        call self%iapws_wrapper%initialize(self%water, self%ice)
        allocate (self%materials(self%num_materials), stat=stat)
        if (stat /= 0) then
            error stop "Error in properties_manager%initialize: Failed to allocate materials array."
        end if
        allocate (self%models(self%num_materials), stat=stat)
        if (stat /= 0) then
            error stop "Error in properties_manager%initialize: Failed to allocate models array."
        end if

        ! Density
        if (present(configs_density)) then
            do model_idx = 1, self%num_materials
                current_material_id = unique_material_ids(model_idx)
                call self%materials(model_idx)%initialize(material_id=current_material_id, &
                                                          config_den=configs_density(model_idx), &
                                                          water=self%water, ice=self%ice)
            end do
        end if

        ! Specific Heat
        if (present(configs_specific_heat)) then
            do model_idx = 1, self%num_materials
                current_material_id = unique_material_ids(model_idx)
                call self%materials(model_idx)%initialize(material_id=current_material_id, &
                                                          config_sph=configs_specific_heat(model_idx), &
                                                          water=self%water, ice=self%ice)
            end do
        end if

        ! Volumetric Heat Capacity
        if (present(configs_heat_capacity)) then
            do model_idx = 1, self%num_materials
                current_material_id = unique_material_ids(model_idx)
                call self%materials(model_idx)%initialize(material_id=current_material_id, &
                                                          config_vhc=configs_heat_capacity(model_idx), &
                                                          water=self%water, ice=self%ice)
            end do
        end if

        ! Thermal Conductivity
        if (present(configs_thermal_conductivity)) then
            do model_idx = 1, self%num_materials
                current_material_id = unique_material_ids(model_idx)
                call self%materials(model_idx)%initialize(material_id=current_material_id, &
                                                          config_thc=configs_thermal_conductivity(model_idx), &
                                                          water=self%water, ice=self%ice)
            end do
        end if

        ! WRF Model
        if (present(configs_wrf)) then
            do model_idx = 1, self%num_materials
                current_material_id = unique_material_ids(model_idx)
                call self%models(model_idx)%initialize(material_id=current_material_id, &
                                                       config_wrf=configs_wrf(model_idx), &
                                                       water=self%water, ice=self%ice)
            end do
        end if

        ! HCF Model
        if (present(configs_hcf)) then
            do model_idx = 1, self%num_materials
                current_material_id = unique_material_ids(model_idx)
                call self%models(model_idx)%initialize(material_id=current_material_id, &
                                                       config_hcf=configs_hcf(model_idx), &
                                                       water=self%water, ice=self%ice)
            end do
        end if

        ! GCC Model
        if (present(configs_gcc)) then
            do model_idx = 1, self%num_materials
                current_material_id = unique_material_ids(model_idx)
                call self%models(model_idx)%initialize(material_id=current_material_id, &
                                                       config_gcc=configs_gcc(model_idx), &
                                                       water=self%water, ice=self%ice)
            end do
        end if

    end subroutine initialize_constitutive_manager

    subroutine set_map_materials(self, unique_material_ids)
        implicit none
        class(type_constitutive_manager), intent(inout) :: self
        integer(int32), intent(in) :: unique_material_ids(:)

        integer(int32) :: max_region_id
        integer(int32) :: model_idx

        self%num_materials = size(unique_material_ids)
        max_region_id = maxval(unique_material_ids)

        allocate (self%materials_id_map(max_region_id), source=0)

        do model_idx = 1, self%num_materials
            self%materials_id_map(unique_material_ids(model_idx)) = model_idx
        end do

    end subroutine set_map_materials

    !-------------------------------------------------------------------------------------------------------------------------------
    ! Generic Type-Bound Procedures
    !-------------------------------------------------------------------------------------------------------------------------------
    subroutine calc_density(self, material_id, state, density)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: density

        call self%materials(self%materials_id_map(material_id))%calc_density(state, density)

    end subroutine calc_density

    subroutine get_density_solid(self, material_id, density_solid)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        real(real64), intent(inout) :: density_solid

        call self%materials(self%materials_id_map(material_id))%get_density_solid(density_solid)

    end subroutine get_density_solid

    subroutine calc_density_water(self, state, density_water, tangent)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: density_water
        !> \((\partial\rho_w/\partial p_w,\partial\rho_w/\partial T)\); overwritten on exit when present
        real(real64), intent(inout), optional :: tangent(2)

        call self%iapws_wrapper%calc_rho_water(state, density_water)
        if (present(tangent)) then
            call self%iapws_wrapper%calc_drho_water_dP(state, tangent(tangent_pressure))
            call self%iapws_wrapper%calc_drho_water_dT(state, tangent(tangent_temperature))
        end if

    end subroutine calc_density_water

    subroutine calc_density_ice(self, state, density_ice, tangent)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: density_ice
        !> \((\partial\rho_i/\partial p_w,\partial\rho_i/\partial T)\); overwritten on exit when present
        real(real64), intent(inout), optional :: tangent(2)

        call self%iapws_wrapper%calc_rho_ice(state, density_ice)
        if (present(tangent)) then
            call self%iapws_wrapper%calc_drho_ice_dP(state, tangent(tangent_pressure))
            call self%iapws_wrapper%calc_drho_ice_dT(state, tangent(tangent_temperature))
        end if

    end subroutine calc_density_ice

    subroutine calc_density_vapor(self, state, density_vapor)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: density_vapor

        call self%iapws_wrapper%calc_rho_vapor(state, density_vapor)

    end subroutine calc_density_vapor

    subroutine calc_density_vapor_saturation(self, state, density_vapor_saturation)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: density_vapor_saturation

        call self%iapws_wrapper%calc_rho_vapor_saturation(state, density_vapor_saturation)

    end subroutine calc_density_vapor_saturation

    subroutine calc_density_water_derivatives(self, material_id, state, dden_dT, dden_dP)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout), optional :: dden_dT
        real(real64), intent(inout), optional :: dden_dP

        call self%materials(self%materials_id_map(material_id))%calc_density_water_derivatives(state, dden_dT, dden_dP)
    end subroutine calc_density_water_derivatives

    subroutine calc_density_ice_derivatives(self, material_id, state, dden_dT, dden_dP)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout), optional :: dden_dT
        real(real64), intent(inout), optional :: dden_dP

        call self%materials(self%materials_id_map(material_id))%calc_density_ice_derivatives(state, dden_dT, dden_dP)
    end subroutine calc_density_ice_derivatives

    subroutine calc_density_vapor_derivatives(self, material_id, state, dden_dT, dden_dP)

        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout), optional :: dden_dT
        real(real64), intent(inout), optional :: dden_dP

        call self%materials(self%materials_id_map(material_id))%calc_density_vapor_derivatives(state, dden_dT, dden_dP)
    end subroutine calc_density_vapor_derivatives

    subroutine calc_specific_heat(self, material_id, state, specific_heat)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: specific_heat

        call self%materials(self%materials_id_map(material_id))%calc_specific_heat(state, specific_heat)

    end subroutine calc_specific_heat

    subroutine get_specific_heat_solid(self, material_id, cp)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        real(real64), intent(inout) :: cp

        call self%materials(self%materials_id_map(material_id))%get_specific_heat_solid(cp)

    end subroutine get_specific_heat_solid

    subroutine calc_specific_heat_water(self, state, specific_heat_water, tangent)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: specific_heat_water
        !> \((\partial c_w/\partial p_w,\partial c_w/\partial T)\); overwritten on exit when present
        real(real64), intent(inout), optional :: tangent(2)

        call self%iapws_wrapper%calc_cp_water(state, specific_heat_water)
        if (present(tangent)) then
            call self%iapws_wrapper%calc_dcp_water_dP(state, tangent(tangent_pressure))
            call self%iapws_wrapper%calc_dcp_water_dT(state, tangent(tangent_temperature))
        end if

    end subroutine calc_specific_heat_water

    subroutine calc_specific_heat_ice(self, state, specific_heat_ice, tangent)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: specific_heat_ice
        !> \((\partial c_i/\partial p_w,\partial c_i/\partial T)\); overwritten on exit when present
        real(real64), intent(inout), optional :: tangent(2)

        call self%iapws_wrapper%calc_cp_ice(state, specific_heat_ice)
        if (present(tangent)) then
            call self%iapws_wrapper%calc_dcp_ice_dP(state, tangent(tangent_pressure))
            call self%iapws_wrapper%calc_dcp_ice_dT(state, tangent(tangent_temperature))
        end if

    end subroutine calc_specific_heat_ice

    subroutine calc_specific_heat_vapor(self, state, specific_heat_vapor, tangent)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: specific_heat_vapor
        !> \((0,\mathrm{d}c_v/\mathrm{d}T)\): saturated vapour depends on \(T\) alone; overwritten on exit when present
        real(real64), intent(inout), optional :: tangent(2)

        call self%iapws_wrapper%calc_cp_vapor(state, specific_heat_vapor)
        if (present(tangent)) then
            tangent(tangent_pressure) = 0.0d0
            call self%iapws_wrapper%calc_dcp_vapor_dT(state, tangent(tangent_temperature))
        end if

    end subroutine calc_specific_heat_vapor

    subroutine calc_thermal_conductivity_nondispersivity(self, material_id, state, thermal_conductivity, tangent)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: thermal_conductivity
        !> \((\partial\lambda/\partial p_w,\partial\lambda/\partial T)\); overwritten on exit when present
        real(real64), intent(inout), optional :: tangent(2)

        call self%materials(self%materials_id_map(material_id))%calc_thermal_conductivity(state, thermal_conductivity, &
                                                                                          tangent)

    end subroutine calc_thermal_conductivity_nondispersivity

    subroutine calc_thermal_conductivity_dispersivity(self, material_id, state, thermal_conductivity)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        type(type_state_thc), intent(inout) :: thermal_conductivity

        call self%materials(self%materials_id_map(material_id))%calc_thermal_conductivity(state, thermal_conductivity)

    end subroutine calc_thermal_conductivity_dispersivity

    subroutine calc_vol_heat_capacity(self, material_id, state, vol_heat_capacity)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: vol_heat_capacity

        call self%materials(self%materials_id_map(material_id))%calc_vol_heat_capacity(state, vol_heat_capacity)

    end subroutine calc_vol_heat_capacity

    subroutine calc_latent_heat_fusion(self, material_id, state, latent_heat_fusion)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: latent_heat_fusion

        call self%models(self%materials_id_map(material_id))%calc_latent_heat_fusion(state, latent_heat_fusion)
    end subroutine calc_latent_heat_fusion

    subroutine calc_latent_heat_vaporization(self, material_id, state, latent_heat_vaporization, tangent)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: latent_heat_vaporization
        !> \((\partial L_v/\partial p_w,\partial L_v/\partial T)\); overwritten on exit when present
        real(real64), intent(inout), optional :: tangent(2)

        call self%models(self%materials_id_map(material_id))%calc_latent_heat_vaporization(state, latent_heat_vaporization, &
                                                                                           tangent)
    end subroutine calc_latent_heat_vaporization

    subroutine calc_pressure_ice_water_derivative(self, material_id, state, deriv)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: deriv

        call self%models(self%materials_id_map(material_id))%calc_pressure_ice_water_derivative(state, deriv)
    end subroutine calc_pressure_ice_water_derivative

    ! subroutine update_water_phases(self, material_id, state)
    subroutine update_water_phases(self, material_id, state)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(inout) :: state

        call self%models(self%materials_id_map(material_id))%update_water_phases(state)
    end subroutine update_water_phases

    subroutine project_ice_content(self, material_id, state, projected_ice, ice_increment, equilibrium_error, active_bound, &
                                    dice_dT, dice_dP)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: projected_ice
        real(real64), intent(inout) :: ice_increment
        real(real64), intent(inout) :: equilibrium_error
        integer(int32), intent(inout), optional :: active_bound
        real(real64), intent(inout), optional :: dice_dT
        real(real64), intent(inout), optional :: dice_dP

        call self%models(self%materials_id_map(material_id))%project_ice_content( &
            state, projected_ice, ice_increment, equilibrium_error, active_bound, dice_dT, dice_dP)
    end subroutine project_ice_content

    subroutine calc_conserved_target(self, material_id, state, target_total_water, available)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: target_total_water
        logical, intent(inout) :: available

        call self%models(self%materials_id_map(material_id))%calc_conserved_target(state, target_total_water, available)
    end subroutine calc_conserved_target

    subroutine solve_local_conserved_equilibrium(self, material_id, state, target_total_water, new_pressure, new_ice, &
                                                  converged)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(in) :: target_total_water
        real(real64), intent(inout) :: new_pressure
        real(real64), intent(inout) :: new_ice
        logical, intent(inout) :: converged

        call self%models(self%materials_id_map(material_id))%solve_local_conserved_equilibrium( &
            state, target_total_water, new_pressure, new_ice, converged)
    end subroutine solve_local_conserved_equilibrium

    subroutine calc_Kflh(self, material_id, state, Kflh, log_tangent)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: Kflh
        !> \((\partial\ln K/\partial p_w,\partial\ln K/\partial T)\); overwritten on exit when present
        real(real64), intent(inout), optional :: log_tangent(2)

        call self%models(self%materials_id_map(material_id))%calc_Kflh(state, Kflh, log_tangent)
    end subroutine calc_Kflh

    subroutine calc_KlT(self, material_id, state, KlT)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: KlT

        call self%models(self%materials_id_map(material_id))%calc_KlT(state, KlT)
    end subroutine calc_KlT

    subroutine calc_cryo_head_dT(self, material_id, state, dh_dT)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: dh_dT

        call self%models(self%materials_id_map(material_id))%calc_cryo_head_dT(state, dh_dT)
    end subroutine calc_cryo_head_dT

    !> Freezing-interface level set consistent with calc_phase_split's own
    !> switch; see type_fusion%calc_freezing_level_set for the definition.
    subroutine calc_freezing_level_set(self, material_id, state, phi)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: phi

        call self%models(self%materials_id_map(material_id))%calc_freezing_level_set(state, phi)
    end subroutine calc_freezing_level_set

    subroutine calc_Kvh(self, material_id, state, Kvh)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: Kvh

        call self%models(self%materials_id_map(material_id))%calc_Kvh(state, Kvh)
    end subroutine calc_Kvh

    subroutine calc_KvT(self, material_id, state, KvT)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: KvT

        call self%models(self%materials_id_map(material_id))%calc_KvT(state, KvT)
    end subroutine calc_KvT

    subroutine calc_saturation_pressure(self, material_id, state, saturation_pressure, is_saturated, tangent)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: saturation_pressure
        logical, intent(inout) :: is_saturated
        !> \((\partial p_{sat}/\partial p_w,\partial p_{sat}/\partial T)\); overwritten on exit when present
        real(real64), intent(inout), optional :: tangent(2)

        call self%models(self%materials_id_map(material_id))%calc_saturation_pressure( &
            state, saturation_pressure, is_saturated, tangent)
    end subroutine calc_saturation_pressure

    subroutine calc_cryo_suction_deriv_T(self, material_id, state, deriv)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: deriv

        call self%models(self%materials_id_map(material_id))%calc_cryo_suction_deriv_T(state, deriv)
    end subroutine calc_cryo_suction_deriv_T

    subroutine calc_lscheme_capacity(self, material_id, capacity)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        real(real64), intent(inout) :: capacity

        call self%models(self%materials_id_map(material_id))%calc_lscheme_capacity(capacity)
    end subroutine calc_lscheme_capacity

    subroutine calc_suction_weights(self, material_id, state, w_cap, w_cryo)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: w_cap, w_cryo

        call self%models(self%materials_id_map(material_id))%calc_suction_weights(state, w_cap, w_cryo)
    end subroutine calc_suction_weights

    subroutine calc_segregation_sink(self, material_id, state, grad_T_magnitude, S_seg)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(in) :: grad_T_magnitude
        real(real64), intent(inout) :: S_seg

        call self%models(self%materials_id_map(material_id))%calc_segregation_sink(state, grad_T_magnitude, S_seg)
    end subroutine calc_segregation_sink

    subroutine calc_effective_segregation_sink(self, material_id, state, grad_T_magnitude, dt, S_seg)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(in) :: grad_T_magnitude
        real(real64), intent(in) :: dt
        real(real64), intent(inout) :: S_seg

        real(real64) :: S_seg_raw
        real(real64) :: Qw, Qi_pore, porosity_val
        real(real64) :: rho_w, rho_i, density_ratio
        real(real64) :: delta_liquid_raw, delta_liquid_clamped, pore_space_left

        S_seg = 0.0d0
        if (grad_T_magnitude <= 0.0d0 .or. dt <= 0.0d0) return

        S_seg_raw = 0.0d0
        call self%calc_segregation_sink(material_id, state, grad_T_magnitude, S_seg_raw)
        if (S_seg_raw <= 0.0d0) return

        call state%water_content%get(Qw)
        call state%ice_content%get(Qi_pore)
        call state%porosity%get(porosity_val)
        call self%calc_density_water(state, rho_w)
        call self%calc_density_ice(state, rho_i)
        if (rho_w <= tiny(1.0d0)) return
        if (rho_i <= tiny(1.0d0)) return

        density_ratio = rho_w / rho_i

        ! Clamp the liquid-water intake consistently for all equations. S_seg is
        ! a liquid-water-equivalent volumetric removal rate [1/s]; it cannot
        ! consume more liquid than is present, and its ice-volume equivalent
        ! cannot exceed the currently open pore space during this step.
        delta_liquid_raw = S_seg_raw * dt
        delta_liquid_clamped = min(delta_liquid_raw, max(Qw, 0.0d0))
        pore_space_left = max(porosity_val - Qi_pore - max(Qw, 0.0d0), 0.0d0)
        delta_liquid_clamped = min(delta_liquid_clamped, pore_space_left / density_ratio)
        delta_liquid_clamped = max(delta_liquid_clamped, 0.0d0)

        S_seg = delta_liquid_clamped / dt
    end subroutine calc_effective_segregation_sink

    function is_segregation_active(self, material_id) result(active)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        logical :: active

        active = self%models(self%materials_id_map(material_id))%is_segregation_active()
    end function is_segregation_active

    pure function has_cryo_transport(self, material_id) result(active)
        implicit none
        class(type_constitutive_manager), intent(in) :: self
        integer(int32), intent(in) :: material_id
        logical :: active

        active = self%models(self%materials_id_map(material_id))%has_cryo_transport()
    end function has_cryo_transport

end module physics_constitutive_manager
