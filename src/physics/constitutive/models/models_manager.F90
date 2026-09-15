module constitutive_models_manager
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: iapws, only:type_iapws97, type_iapws06
    use :: constitutive_constants, only: SUCTION_BLEND_EPS => suction_blend_epsilon
    use :: module_core
    use :: models_wrf
    use :: models_hcf
    use :: models_phase_change_gcc
    use :: models_phase_change_fusion
    use :: models_phase_change_vaporization
    use :: models_phase_change_manager
    use :: models_segregation_potential, only: type_segregation_potential
    implicit none
    private

    public :: type_models_manager

    type :: type_models_manager
        private
        type(holder_wrfs) :: wrf
        type(holder_hcfs) :: hcf
        type(holder_gccs) :: gcc
        type(type_phase_manager) :: phase_manager
        type(type_segregation_potential) :: segregation
    contains
        procedure, public :: initialize
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
        procedure, public :: calc_latent_heat_fusion
        procedure, public :: calc_latent_heat_vaporization
        procedure, public :: calc_saturation_pressure
        procedure, public :: calc_pressure_ice_water_derivative
        procedure, public :: calc_cryo_suction_deriv_T
        procedure, public :: calc_lscheme_capacity
        procedure, public :: calc_suction_weights
        procedure, public :: calc_segregation_sink
        procedure, public :: is_segregation_active
        procedure, public :: has_cryo_transport
    end type type_models_manager

contains

    subroutine initialize(self, material_id, config_wrf, config_hcf, config_gcc, water, ice)
        implicit none
        class(type_models_manager), intent(inout) :: self
        integer(int32), intent(in) :: material_id
        type(type_config_wrf), intent(in), optional :: config_wrf
        type(type_config_hcf), intent(in), optional :: config_hcf
        type(type_config_gcc), intent(in), optional :: config_gcc
        type(type_iapws97), intent(in), optional, target :: water
        type(type_iapws06), intent(in), optional, target :: ice

        if (present(config_wrf)) then
            call self%wrf%initialize(config_wrf)
        end if
        if (present(config_hcf) .and. present(water)) then
            call self%hcf%initialize(config_hcf, water, ice)
        end if
        if (present(config_gcc) .and. present(water) .and. present(ice)) then
            call self%gcc%initialize(material_id, config_gcc, water, ice)
            if (config_gcc%segregation_potential > 0.0d0) then
                call self%segregation%initialize(config_gcc%segregation_potential, &
                                                 T_fringe_low=config_gcc%T_fringe_low, &
                                                 T_fringe_high=config_gcc%T_fringe_high)
            end if
        end if

        if (allocated(self%wrf%p) .and. allocated(self%gcc%p)) then
            if (self%wrf%p%is_initialized() .and. self%gcc%p%is_initialized()) then
                call self%phase_manager%initialize(self%gcc%p, self%wrf%p, water, ice)
            end if
        end if
    end subroutine initialize

    subroutine update_water_phases(self, state)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(inout) :: state

        call self%phase_manager%update_water_phases(state)
    end subroutine update_water_phases

    subroutine project_ice_content(self, state, projected_ice, ice_increment, equilibrium_error, active_bound, &
                                    dice_dT, dice_dP)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: projected_ice
        real(real64), intent(inout) :: ice_increment
        real(real64), intent(inout) :: equilibrium_error
        integer(int32), intent(inout), optional :: active_bound
        real(real64), intent(inout), optional :: dice_dT
        real(real64), intent(inout), optional :: dice_dP

        call self%phase_manager%project_ice_content(state, projected_ice, ice_increment, equilibrium_error, active_bound, &
                                                      dice_dT, dice_dP)
    end subroutine project_ice_content

    subroutine calc_conserved_target(self, state, target_total_water, available)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: target_total_water
        logical, intent(inout) :: available

        call self%phase_manager%calc_conserved_target(state, target_total_water, available)
    end subroutine calc_conserved_target

    subroutine solve_local_conserved_equilibrium(self, state, target_total_water, new_pressure, new_ice, converged)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(in) :: target_total_water
        real(real64), intent(inout) :: new_pressure
        real(real64), intent(inout) :: new_ice
        logical, intent(inout) :: converged

        call self%phase_manager%solve_local_conserved_equilibrium(state, target_total_water, new_pressure, new_ice, converged)
    end subroutine solve_local_conserved_equilibrium

    subroutine calc_Kflh(self, state, Kflh)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: Kflh

        call self%hcf%p%calc_Kflh(state, Kflh)
    end subroutine calc_Kflh

    !> Thermal liquid conductivity from the surface-tension relation.
    subroutine calc_KlT(self, state, KlT)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: KlT

        call self%hcf%p%calc_KlT(state, KlT)
    end subroutine calc_KlT

    !> dh/dT [m/K] of the generalized head h = -psi_eff/(rho_std g): the
    !> temperature sensitivity of the total-potential head that drives the
    !> liquid flux. Used to assemble the consistent K_HT coupling tangent so
    !> the pressure solve accounts for the cryosuction flux's dependence on
    !> the freezing-front temperature.
    subroutine calc_cryo_head_dT(self, state, dh_dT)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: dh_dT

        call self%phase_manager%calc_cryo_head_dT(state, dh_dT)
    end subroutine calc_cryo_head_dT

    !> Freezing-interface level set consistent with calc_phase_split's own
    !> switch; see type_fusion%calc_freezing_level_set for the definition.
    subroutine calc_freezing_level_set(self, state, phi)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: phi

        call self%phase_manager%calc_freezing_level_set(state, phi)
    end subroutine calc_freezing_level_set

    subroutine calc_Kvh(self, state, Kvh)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: Kvh

        call self%hcf%p%calc_Kvh(state, Kvh)
    end subroutine calc_Kvh

    subroutine calc_KvT(self, state, KvT)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: KvT

        call self%hcf%p%calc_KvT(state, KvT)
    end subroutine calc_KvT

    subroutine calc_latent_heat_fusion(self, state, L_fusion)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: L_fusion

        call self%phase_manager%calc_latent_heat_fusion(state, L_fusion)
    end subroutine calc_latent_heat_fusion

    subroutine calc_latent_heat_vaporization(self, state, L_vaporization)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: L_vaporization

        call self%phase_manager%calc_latent_heat_vaporization(state, L_vaporization)
    end subroutine calc_latent_heat_vaporization

    subroutine calc_saturation_pressure(self, state, saturation_pressure, is_saturated)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: saturation_pressure
        logical, intent(inout) :: is_saturated

        call self%phase_manager%calc_saturation_pressure(state, saturation_pressure, is_saturated)
    end subroutine calc_saturation_pressure

    subroutine calc_pressure_ice_water_derivative(self, state, deriv)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: deriv

        call self%phase_manager%deriv_pressure_ice_water(state, deriv)
    end subroutine calc_pressure_ice_water_derivative

    subroutine calc_cryo_suction_deriv_T(self, state, deriv)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: deriv

        call self%gcc%deriv_temperature(state, deriv)
    end subroutine calc_cryo_suction_deriv_T

    subroutine calc_lscheme_capacity(self, capacity)
        implicit none
        class(type_models_manager), intent(in) :: self
        real(real64), intent(inout) :: capacity

        capacity = 0.0d0
        if (allocated(self%wrf%p)) call self%wrf%p%calc_lscheme_capacity(capacity)
    end subroutine calc_lscheme_capacity

    !> Weights of the generalized suction p_c* = max(P_aw, P_iw) with respect to its
    !> two contributions: w_cap = d(p_c*)/d(P_aw), w_cryo = d(p_c*)/d(P_iw). These are
    !> the fractions with which the capillary (grad p_w) and cryogenic (grad T) terms
    !> enter the Darcy flux driven by grad(p_c*) = grad(mu_w) (the chemical potential).
    !> Uses the same smooth max as the water-retention evaluation. With no cryogenic
    !> model (psi_cryo = 0) it returns w_cap = 1, w_cryo = 0 (ordinary unfrozen flow).
    subroutine calc_suction_weights(self, state, w_cap, w_cryo)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: w_cap, w_cryo

        real(real64) :: pressure, psi_cap, psi_cryo, delta, denom

        call state%pressure%get(pressure)
        psi_cap = max(0.0d0, -pressure)
        psi_cryo = 0.0d0
        call self%gcc%calc(state, psi_cryo)

        delta = psi_cap - psi_cryo
        denom = sqrt(delta*delta + SUCTION_BLEND_EPS*SUCTION_BLEND_EPS)
        w_cap = 0.5d0*(1.0d0 + delta/denom)
        w_cryo = 0.5d0*(1.0d0 - delta/denom)
    end subroutine calc_suction_weights

    subroutine calc_segregation_sink(self, state, grad_T_magnitude, S_seg)
        implicit none
        class(type_models_manager), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(in) :: grad_T_magnitude
        real(real64), intent(inout) :: S_seg

        call self%segregation%calc_sink(state, grad_T_magnitude, S_seg)
    end subroutine calc_segregation_sink

    pure function is_segregation_active(self) result(active)
        implicit none
        class(type_models_manager), intent(in) :: self
        logical :: active

        active = self%segregation%is_active()
    end function is_segregation_active

    ! Returns .true. when a GCC model is active. The hydraulic element uses
    ! it does not add a separate cryogenic Darcy driving term.
    pure function has_cryo_transport(self) result(active)
        implicit none
        class(type_models_manager), intent(in) :: self
        logical :: active

        active = allocated(self%gcc%p)
    end function has_cryo_transport

end module constitutive_models_manager
