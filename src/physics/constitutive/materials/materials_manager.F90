module constitutive_materials_manager
    use :: iso_fortran_env, only:int32, real64
    use :: iapws, only:type_iapws97, type_iapws06
    use :: module_core, only:type_config_constitutive, type_state, type_state_thc
    use :: materials_base
    use :: materials_density
    use :: materials_specific_heat
    use :: materials_thermal_conductivity
    use :: materials_heat_capacity
    implicit none
    private

    public :: type_material_manager

    type :: type_material_manager
        private
        type(holder_dens) :: den
        type(holder_sphs) :: sph
        type(holder_vhcs) :: vhc
        type(holder_thcs) :: thc
    contains
        ! ---- Lifecycle ----
        ! initialize, destroy, reset, etc.
        procedure, public :: initialize => initialize_type_material_manager

        ! ---- Mutator ----
        ! set_XXX, increment_XXX, update_XXX, etc.

        ! ---- Algorithm / Operation ----
        ! compute_XXX, check_XXX, solve_XXX, etc.
        procedure, public :: calc_density
        procedure, public :: calc_density_water_derivatives
        procedure, public :: calc_density_ice_derivatives
        procedure, public :: calc_density_vapor_derivatives
        procedure, public :: calc_specific_heat
        generic, public :: calc_thermal_conductivity => calc_thermal_conductivity_nondispersivity, &
            calc_thermal_conductivity_dispersivity
        procedure, private :: calc_thermal_conductivity_nondispersivity
        procedure, private :: calc_thermal_conductivity_dispersivity
        procedure, public :: calc_vol_heat_capacity

        ! ---- Inquiry ----
        ! is_XXX, has_XXX, should_XXX, etc.

        ! ---- Getter ----
        ! get_XXX, etc.
        procedure, public :: get_density_solid
        procedure, public :: get_specific_heat_solid

        ! ---- Meta / Utility ----
        ! display, to_string, etc.

        ! ---- Operator ----
    end type type_material_manager

contains

    subroutine initialize_type_material_manager(self, material_id, config_den, config_sph, config_vhc, config_thc, water, ice)
        implicit none
        class(type_material_manager), intent(inout) :: self
        integer(int32), intent(in) :: material_id
        type(type_config_constitutive), intent(in), optional :: config_den
        type(type_config_constitutive), intent(in), optional :: config_sph
        type(type_config_constitutive), intent(in), optional :: config_vhc
        type(type_config_constitutive), intent(in), optional :: config_thc

        type(type_iapws97), intent(in), optional, target :: water
        type(type_iapws06), intent(in), optional, target :: ice

        logical :: is_present_iapws = .false.

        if (.not. present(water) .or. .not. present(ice)) then
            is_present_iapws = .false.
        else
            is_present_iapws = .true.
        end if

        if (present(config_den)) then
            if (.not. is_present_iapws) then
                error stop "Error in material_manager%initialize: water and ice must be provided when config_den is provided."
            end if
            call self%den%initialize(material_id, config_den, water, ice)
        end if
        if (present(config_sph)) then
            if (.not. is_present_iapws) then
                error stop "Error in material_manager%initialize: water and ice must be provided when config_sph is provided."
            end if
            call self%sph%initialize(material_id, config_sph, water, ice)
        end if
        if (present(config_vhc)) then
            if (.not. is_present_iapws) then
                error stop "Error in material_manager%initialize: water and ice must be provided when config_vhc is provided."
            end if
            call self%vhc%initialize(material_id, config_vhc, water, ice)
        end if
        if (present(config_thc)) then
            if (.not. is_present_iapws) then
                error stop "Error in material_manager%initialize: water and ice must be provided when config_thc is provided."
            end if
            call self%thc%initialize(material_id, config_thc, water, ice)
        end if
    end subroutine initialize_type_material_manager

    subroutine calc_density(self, state, density)
        implicit none
        class(type_material_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: density

        call self%den%p%calc(state, density)

    end subroutine calc_density

    subroutine get_density_solid(self, density_solid)
        implicit none
        class(type_material_manager), intent(in) :: self
        real(real64), intent(inout) :: density_solid

        call self%den%p%get_solid(density_solid)

    end subroutine get_density_solid

    subroutine get_specific_heat_solid(self, cp)
        implicit none
        class(type_material_manager), intent(in) :: self
        real(real64), intent(inout) :: cp

        call self%sph%p%get_solid(cp)

    end subroutine get_specific_heat_solid

    subroutine calc_density_water_derivatives(self, state, dden_dT, dden_dP)
        implicit none
        class(type_material_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout), optional :: dden_dT
        real(real64), intent(inout), optional :: dden_dP

        if (present(dden_dT)) then
            call self%den%p%calc_drho_water_dT(state, dden_dT)
        end if
        if (present(dden_dP)) then
            call self%den%p%calc_drho_water_dP(state, dden_dP)
        end if

    end subroutine calc_density_water_derivatives

    subroutine calc_density_ice_derivatives(self, state, dden_dT, dden_dP)
        implicit none
        class(type_material_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout), optional :: dden_dT
        real(real64), intent(inout), optional :: dden_dP

        if (present(dden_dT)) then
            call self%den%p%calc_drho_ice_dT(state, dden_dT)
        end if
        if (present(dden_dP)) then
            call self%den%p%calc_drho_ice_dP(state, dden_dP)
        end if

    end subroutine calc_density_ice_derivatives

    subroutine calc_density_vapor_derivatives(self, state, dden_dT, dden_dP)
        implicit none
        class(type_material_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout), optional :: dden_dT
        real(real64), intent(inout), optional :: dden_dP

        if (present(dden_dT)) then
            call self%den%p%calc_drho_vapor_dT(state, dden_dT)
        end if
        if (present(dden_dP)) then
            call self%den%p%calc_drho_vapor_dP(state, dden_dP)
        end if

    end subroutine calc_density_vapor_derivatives

    subroutine calc_specific_heat(self, state, specific_heat)
        implicit none
        class(type_material_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: specific_heat

        call self%sph%p%calc(state, specific_heat)

    end subroutine calc_specific_heat

    subroutine calc_thermal_conductivity_nondispersivity(self, state, lambda, tangent)
        implicit none
        class(type_material_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: lambda
        !> \((\partial\lambda/\partial p_w,\partial\lambda/\partial T)\); overwritten on exit when present
        real(real64), intent(inout), optional :: tangent(2)

        call self%thc%p%calc(state, lambda, tangent)

    end subroutine calc_thermal_conductivity_nondispersivity

    subroutine calc_thermal_conductivity_dispersivity(self, state, lambda)
        implicit none
        class(type_material_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        type(type_state_thc), intent(inout) :: lambda

        call self%thc%p%calc(state, lambda)

    end subroutine calc_thermal_conductivity_dispersivity

    subroutine calc_vol_heat_capacity(self, state, vhc)
        implicit none
        class(type_material_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: vhc

        call self%vhc%p%calc(state, vhc)

    end subroutine calc_vol_heat_capacity
end module constitutive_materials_manager
