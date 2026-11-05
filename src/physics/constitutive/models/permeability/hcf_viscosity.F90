submodule(models_hcf) hcf_viscosity
    use :: constitutive_laws, only:evaluate_water_viscosity_exponential, evaluate_water_viscosity_exponential_tangent, &
        evaluate_water_viscosity_supercooled, evaluate_water_viscosity_supercooled_tangent
    implicit none
contains
    module subroutine initialize_abst_hcf_viscosity(self, temperature_critical)
        implicit none
        class(abst_hcf_viscosity), intent(inout) :: self
        real(real64), intent(in), optional :: temperature_critical

        real(real64) :: temp_critical

        if (present(temperature_critical)) then
            temp_critical = temperature_critical
        else
            temp_critical = 15.0d0
        end if

        call self%calc_mu(temp_critical, self%mu_zero)
    end subroutine initialize_abst_hcf_viscosity

    module subroutine calc_kr_abst_hcf_viscosity(self, temperature, kr, dkr_dT)
        implicit none
        class(abst_hcf_viscosity), intent(in) :: self
        real(real64), intent(in) :: temperature
        real(real64), intent(inout) :: kr
        real(real64), intent(inout), optional :: dkr_dT

        real(real64) :: mu, dmu_dT

        if (present(dkr_dT)) then
            call self%calc_mu(temperature, mu, dmu_dT)
            kr = self%mu_zero / mu
            dkr_dT = -kr * dmu_dT / mu
        else
            call self%calc_mu(temperature, mu)
            kr = self%mu_zero / mu
        end if

    end subroutine calc_kr_abst_hcf_viscosity

    module subroutine calc_mu_exponential(self, temperature, mu, dmu_dT)
        implicit none
        class(type_hcf_viscosity_exp), intent(in) :: self
        real(real64), intent(in) :: temperature
        real(real64), intent(inout) :: mu
        real(real64), intent(inout), optional :: dmu_dT

        real(real64) :: dmu(1)

        if (present(dmu_dT)) then
            call evaluate_water_viscosity_exponential_tangent(temperature, [1.0d0], mu, dmu)
            dmu_dT = dmu(1)
        else
            call evaluate_water_viscosity_exponential(temperature, mu)
        end if
    end subroutine calc_mu_exponential

    module subroutine calc_mu_exponential_supercooled(self, temperature, mu, dmu_dT)
        implicit none
        class(type_hcf_viscosity_supercool), intent(in) :: self
        real(real64), intent(in) :: temperature
        real(real64), intent(inout) :: mu
        real(real64), intent(inout), optional :: dmu_dT

        real(real64) :: dmu(1)

        if (present(dmu_dT)) then
            call evaluate_water_viscosity_supercooled_tangent(temperature, [1.0d0], mu, dmu)
            dmu_dT = dmu(1)
        else
            call evaluate_water_viscosity_supercooled(temperature, mu)
        end if
    end subroutine calc_mu_exponential_supercooled

end submodule hcf_viscosity
