!>
!> @brief Physics models for vaporization.
!>
module models_phase_change_vaporization
    use, intrinsic :: iso_fortran_env
    use :: iapws, only:type_iapws97
    use :: module_core, only:type_state
    use :: constitutive_constants, only: &
        g => gravity_acceleration, &
        T_to_K => celsius_to_kelvin, &
        Rg => universal_gas_constant, &
        Mw => molar_mass_water, &
        rho_std => reference_water_density, &
        P_atm => standard_atmospheric_pressure, &
        water_triple_point_temperature, water_critical_point_temperature
    use :: physics_constitutive_base, only:abst_constitutive

    implicit none
    private

    public :: type_evaporation

    type, extends(abst_constitutive) :: type_evaporation
    contains
        procedure, pass(self), public :: initialize => initialize_evaporation_model
        procedure, pass(self), public :: calc_latent_heat_vaporization
        procedure, pass(self), public :: calc_relative_humidity => calc_relative_humidity_evaporation
        procedure, pass(self), public :: deriv_relative_humidity_temperature => deriv_relative_humidity_temperature_evaporation
        procedure, pass(self), public :: deriv_relative_humidity_pressure => deriv_relative_humidity_pressure_evaporation

        !> Calculate vapor content (Liquid equivalent volume fraction)
        procedure, pass(self), public :: calc_vapor_content => calc_vapor_content_vaporization

        !> Calculate derivatives of vapor content w.r.t P and T
        procedure, pass(self), public :: calc_vapor_content_derivatives
        procedure, pass(self), public :: calc_vapor_content_with_derivatives
    end type type_evaporation

contains

    subroutine initialize_evaporation_model(self, water)
        implicit none
        class(type_evaporation), intent(inout) :: self
        type(type_iapws97), intent(in), target, optional :: water

        if (present(water)) self%water => water

        self%initialized = .true.
    end subroutine initialize_evaporation_model

    subroutine calc_latent_heat_vaporization(self, temperature, latent_heat, dlatent_heat_dT)
        implicit none
        class(type_evaporation), intent(in) :: self
        real(real64), intent(in) :: temperature
        real(real64), intent(inout) :: latent_heat
        !> \(\mathrm{d}L_v/\mathrm{d}T\) along the saturation line, J/(kg K); overwritten on exit when present
        real(real64), intent(inout), optional :: dlatent_heat_dT

        real(real64) :: temperature_K

        if (associated(self%water)) then
            call self%shift_temperature_absolute(temperature, temperature_K)
            ! The saturation line of IAPWS-97 runs from the triple to the
            ! critical point; outside it the latent heat is held at the bound,
            ! as every other saturated-vapour property of the closure is.
            if (temperature_K < water_triple_point_temperature .or. &
                temperature_K > water_critical_point_temperature) then
                temperature_K = min(max(temperature_K, water_triple_point_temperature), water_critical_point_temperature)
                call self%calc_latent_heat_saturation(temperature_K, latent_heat)
                if (present(dlatent_heat_dT)) dlatent_heat_dT = 0.0d0
            else
                call self%calc_latent_heat_saturation(temperature_K, latent_heat)
                if (present(dlatent_heat_dT)) call self%calc_dlatent_heat_saturation_dT(temperature_K, dlatent_heat_dT)
            end if
        else
            latent_heat = 2.501d6 - 2369.2d0 * temperature
            if (present(dlatent_heat_dT)) dlatent_heat_dT = -2369.2d0
        end if
    end subroutine calc_latent_heat_vaporization

    !> @brief Calculate Relative Humidity (Kelvin equation)
    subroutine calc_relative_humidity_evaporation(self, state, relative_humidity)
        implicit none
        class(type_evaporation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: relative_humidity

        real(real64) :: temperature, temperature_K
        real(real64) :: pressure

        call state%temperature%get(temperature)
        call state%pressure%get(pressure)

        if (pressure >= 0.0d0) then
            relative_humidity = 1.0d0
            return
        end if

        call self%shift_temperature_absolute(temperature, temperature_K)
        ! Kelvin equation: RH = exp( P * Mw / (rho * R * T) )
        ! Using rho_std for liquid density in the exponent is a standard approximation
        relative_humidity = exp((pressure * Mw) / (rho_std * Rg * temperature_K))

    end subroutine calc_relative_humidity_evaporation

    !> @brief d(RH)/dT
    subroutine deriv_relative_humidity_temperature_evaporation(self, state, deriv_rh_temp)
        implicit none
        class(type_evaporation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: deriv_rh_temp
        real(real64) :: rh, temperature, temperature_K
        real(real64) :: pressure

        call state%temperature%get(temperature)
        call state%pressure%get(pressure)

        if (pressure >= 0.0d0) then
            deriv_rh_temp = 0.0d0
            return
        end if

        call self%shift_temperature_absolute(temperature, temperature_K)
        call self%calc_relative_humidity(state, rh)

        ! d(exp(A/T))/dT = exp(A/T) * (-A/T^2)
        deriv_rh_temp = rh * (-pressure * Mw) / (rho_std * Rg * temperature_K**2)
    end subroutine deriv_relative_humidity_temperature_evaporation

    !> @brief d(RH)/dP
    subroutine deriv_relative_humidity_pressure_evaporation(self, state, deriv_rh_pressure)
        implicit none
        class(type_evaporation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: deriv_rh_pressure
        real(real64) :: rh, temperature, temperature_K
        real(real64) :: pressure

        call state%temperature%get(temperature)
        call state%pressure%get(pressure)

        if (pressure >= 0.0d0) then
            deriv_rh_pressure = 0.0d0
            return
        end if

        call self%shift_temperature_absolute(temperature, temperature_K)
        call self%calc_relative_humidity(state, rh)

        ! d(exp(C*P))/dP = C * exp(C*P)
        deriv_rh_pressure = rh * Mw / (rho_std * Rg * temperature_K)

    end subroutine deriv_relative_humidity_pressure_evaporation

    !> @brief Calculate vapor content (Liquid equivalent).
    !> theta_v = theta_g * (rho_sat * RH / rho_w)
    subroutine calc_vapor_content_vaporization(self, state, vapor_content)
        implicit none
        class(type_evaporation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: vapor_content

        real(real64) :: relative_humidity
        real(real64) :: saturation_density, water_density
        real(real64) :: air_content

        call self%calc_relative_humidity(state, relative_humidity)
        call self%calc_rho_vapor_saturation(state, saturation_density)

        call self%calc_rho_water(state, water_density)
        call state%air_content%get(air_content)

        vapor_content = saturation_density * relative_humidity * air_content / water_density
    end subroutine calc_vapor_content_vaporization

    !>
    !> @brief Compute pressure and temperature derivatives of vapor content (liquid equivalent).
    !>
    !> d(theta_v)/dX = d(theta_g)/dX * (Ratio) + theta_g * d(Ratio)/dX
    !> where Ratio = rho_sat * RH / rho_w
    !>
    subroutine calc_vapor_content_derivatives(self, state, dvapor_dP, dvapor_dT)
        implicit none
        class(type_evaporation), intent(in) :: self
        type(type_state), intent(in) :: state
        !> Output: Gradient of vapor content
        real(real64), intent(inout) :: dvapor_dP, dvapor_dT

        real(real64) :: temperature, pressure, air_content
        real(real64) :: d_air_content_dP, d_air_content_dT

        real(real64) :: temperature_K
        real(real64) :: rh, drh_dP, drh_dT
        real(real64) :: rho_sat, drho_sat_dT
        real(real64) :: rho_w, drho_w_dP, drho_w_dT
        real(real64) :: ratio, dratio_dP, dratio_dT

        ! 1. Prepare state variables
        call state%temperature%get(temperature)
        call state%pressure%get(pressure)
        call state%air_content%get(air_content)
        call state%dQa_dP%get(d_air_content_dP)
        call state%dQa_dT%get(d_air_content_dT)
        call self%shift_temperature_absolute(temperature, temperature_K)

        ! 2. Compute individual terms

        ! Relative Humidity & Derivatives
        call self%calc_relative_humidity(state, rh)
        call self%deriv_relative_humidity_pressure(state, drh_dP)
        call self%deriv_relative_humidity_temperature(state, drh_dT)

        ! Saturation Density (only T dependent)
        call self%calc_rho_vapor_saturation(state, rho_sat)
        call self%calc_drho_vapor_saturation_dT(state, drho_sat_dT)

        ! Water Density (liquid)
        call self%calc_rho_water(state, rho_w)
        call self%calc_drho_water_dT(state, drho_w_dT)
        call self%calc_drho_water_dP(state, drho_w_dP)

        ! 3. Vapor density ratio (Ratio = rho_sat * RH / rho_w) and its derivatives

        ratio = (rho_sat * rh) / rho_w

        ! d(Ratio)/dP = (rho_sat/rho_w)*dRH/dP - (rho_sat*RH/rho_w^2)*drho_w/dP
        dratio_dP = (rho_sat / rho_w) * drh_dP - &
                    (rho_sat * rh / (rho_w**2)) * drho_w_dP

        ! d(Ratio)/dT = (1/rho_w)*(drho_sat/dT*RH + rho_sat*dRH/dT) - (Ratio/rho_w)*drho_w/dT
        dratio_dT = (1.0d0 / rho_w) * (drho_sat_dT * rh + rho_sat * drh_dT) - &
                    (ratio / rho_w) * drho_w_dT

        ! 4. Assemble total derivatives (product rule)
        ! Vapor = Air * Ratio
        ! dVapor/dX = dAir/dX * Ratio + Air * dRatio/dX

        dvapor_dP = d_air_content_dP * ratio + air_content * dratio_dP
        dvapor_dT = d_air_content_dT * ratio + air_content * dratio_dT

    end subroutine calc_vapor_content_derivatives

    !> Compute vapor content and its pressure and temperature derivatives together.
    !>
    !> \[ \theta_v = \theta_g \rho_{v,sat} RH / \rho_w \]
    !> Assumptions: relative humidity and gas-phase state derivatives are current.
    !> Numerical guarantee: algebraically equivalent to the separate value and derivative procedures.
    !> Computational complexity: \(O(1)\) arithmetic and memory.
    !> Failure behavior: propagates the underlying IAPWS property behavior.
    subroutine calc_vapor_content_with_derivatives(self, state, vapor_content, dvapor_dP, dvapor_dT)
        implicit none
        class(type_evaporation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: vapor_content
        real(real64), intent(inout) :: dvapor_dP
        real(real64), intent(inout) :: dvapor_dT

        real(real64) :: temperature, pressure, temperature_K
        real(real64) :: air_content, d_air_content_dP, d_air_content_dT
        real(real64) :: rh, drh_dP, drh_dT
        real(real64) :: rho_sat, drho_sat_dT
        real(real64) :: rho_w, drho_w_dP, drho_w_dT
        real(real64) :: ratio, dratio_dP, dratio_dT

        call state%temperature%get(temperature)
        call state%pressure%get(pressure)
        call state%air_content%get(air_content)
        call state%dQa_dP%get(d_air_content_dP)
        call state%dQa_dT%get(d_air_content_dT)
        call state%relative_humidity%get(rh)
        call self%shift_temperature_absolute(temperature, temperature_K)

        if (pressure >= 0.0d0) then
            drh_dP = 0.0d0
            drh_dT = 0.0d0
        else
            drh_dP = rh * Mw / (rho_std * Rg * temperature_K)
            drh_dT = rh * (-pressure * Mw) / (rho_std * Rg * temperature_K**2)
        end if

        call self%calc_rho_vapor_saturation(state, rho_sat)
        call self%calc_drho_vapor_saturation_dT(state, drho_sat_dT)
        call self%calc_rho_water(state, rho_w)
        call self%calc_drho_water_dT(state, drho_w_dT)
        call self%calc_drho_water_dP(state, drho_w_dP)

        ratio = rho_sat * rh / rho_w
        dratio_dP = (rho_sat / rho_w) * drh_dP - &
                    (rho_sat * rh / rho_w**2) * drho_w_dP
        dratio_dT = (drho_sat_dT * rh + rho_sat * drh_dT) / rho_w - &
                    (ratio / rho_w) * drho_w_dT

        vapor_content = air_content * ratio
        dvapor_dP = d_air_content_dP * ratio + air_content * dratio_dP
        dvapor_dT = d_air_content_dT * ratio + air_content * dratio_dT
    end subroutine calc_vapor_content_with_derivatives

end module models_phase_change_vaporization
