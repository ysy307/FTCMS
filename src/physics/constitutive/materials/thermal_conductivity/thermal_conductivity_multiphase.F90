submodule(materials_thermal_conductivity) thermal_conductivity_multiphase
    use :: constitutive_constants, only:tangent_pressure, tangent_temperature
    use :: constitutive_laws, only:evaluate_thermal_conductivity_geometric2, &
        evaluate_thermal_conductivity_geometric2_tangent, evaluate_thermal_conductivity_geometric3, &
        evaluate_thermal_conductivity_geometric3_tangent, evaluate_thermal_conductivity_geometric4, &
        evaluate_thermal_conductivity_geometric4_tangent, evaluate_thermal_conductivity_campbell_hansson, &
        evaluate_thermal_conductivity_campbell_hansson_tangent
    implicit none

    !> The solid fraction \(1-\phi\) is fixed, so its tangent is zero.
    real(real64), parameter :: no_tangent(2) = 0.0d0

contains

    !> The \((p_w,T)\) tangents of the liquid and ice contents the split published.
    subroutine get_content_tangents(state, dwater, dice)
        implicit none
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: dwater(2), dice(2)

        call state%dQw_dP%get(dwater(tangent_pressure))
        call state%dQw_dT%get(dwater(tangent_temperature))
        call state%dQi_dP%get(dice(tangent_pressure))
        call state%dQi_dT%get(dice(tangent_temperature))
    end subroutine get_content_tangents

    module subroutine calc_thc_gp_1phase(self, state, lambda, tangent)
        implicit none
        class(type_thc_1phase), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: lambda
        real(real64), intent(inout), optional :: tangent(2)

        lambda = self%material1
        if (present(tangent)) tangent(:) = 0.0d0

    end subroutine calc_thc_gp_1phase

    module subroutine calc_thc_dispersivity_gp_1phase(self, state, lambda)
        implicit none
        class(type_thc_1phase), intent(in) :: self
        type(type_state), intent(in) :: state
        type(type_state_thc), intent(inout) :: lambda

        real(real64) :: lambda_0
        type(type_coordinate_dp), pointer :: water_flux

        call state%get(water_flux=water_flux)

        call self%calc(state, lambda_0)
        call self%calc_lambda_dispersivity(lambda_0, 0.0d0, water_flux%x, water_flux%y, water_flux%z, lambda)

    end subroutine calc_thc_dispersivity_gp_1phase

    module subroutine calc_thc_gp_2phase(self, state, lambda, tangent)
        implicit none
        class(type_thc_2phase), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: lambda
        real(real64), intent(inout), optional :: tangent(2)

        real(real64) :: phi1, phi2
        real(real64) :: dwater(2), dice(2)

        call self%get_phi(state, phi1, phi2)

        if (present(tangent)) then
            call get_content_tangents(state, dwater, dice)
            call evaluate_thermal_conductivity_geometric2_tangent(self%material1, self%material2, &
                                                                  phi1, no_tangent, phi2, dwater, lambda, tangent)
        else
            call evaluate_thermal_conductivity_geometric2(self%material1, self%material2, phi1, phi2, lambda)
        end if

    end subroutine calc_thc_gp_2phase

    module subroutine calc_thc_dispersivity_gp_2phase(self, state, lambda)
        implicit none
        class(type_thc_2phase), intent(in) :: self
        type(type_state), intent(in) :: state
        type(type_state_thc), intent(inout) :: lambda

        real(real64) :: lambda_0
        type(type_coordinate_dp), pointer :: water_flux

        real(real64) :: htc_water, rho_water, cp_water

        call state%get(water_flux=water_flux)

        call self%calc_water_density(state, rho_water)
        call self%calc_water_cp(state, cp_water)
        htc_water = rho_water * cp_water

        call self%calc(state, lambda_0)
        call self%calc_lambda_dispersivity(lambda_0, htc_water, water_flux%x, water_flux%y, water_flux%z, lambda)

    end subroutine calc_thc_dispersivity_gp_2phase

    module subroutine calc_thc_gp_3phase(self, state, lambda, tangent)
        implicit none
        class(type_thc_3phase), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: lambda
        real(real64), intent(inout), optional :: tangent(2)

        real(real64) :: phi1, phi2, phi3
        real(real64) :: dwater(2), dice(2)

        call self%get_phi(state, phi1, phi2, phi3)

        if (present(tangent)) then
            call get_content_tangents(state, dwater, dice)
            call evaluate_thermal_conductivity_geometric3_tangent(self%material1, self%material2, self%material3, &
                                                                  phi1, no_tangent, phi2, dwater, phi3, dice, &
                                                                  lambda, tangent)
        else
            call evaluate_thermal_conductivity_geometric3(self%material1, self%material2, self%material3, &
                                                          phi1, phi2, phi3, lambda)
        end if

    end subroutine calc_thc_gp_3phase

    module subroutine calc_thc_dispersivity_gp_3phase(self, state, lambda)
        implicit none
        class(type_thc_3phase), intent(in) :: self
        type(type_state), intent(in) :: state
        type(type_state_thc), intent(inout) :: lambda

        real(real64) :: lambda_0
        type(type_coordinate_dp), pointer :: water_flux

        real(real64) :: qx, qy, qz
        real(real64) :: htc_water, rho_water, cp_water

        call self%calc_water_density(state, rho_water)
        call self%calc_water_cp(state, cp_water)
        htc_water = rho_water * cp_water

        call state%get(water_flux=water_flux)
        if (associated(water_flux)) then
            qx = water_flux%x
            qy = water_flux%y
            qz = water_flux%z
        else
            qx = 0.0d0
            qy = 0.0d0
            qz = 0.0d0
        end if

        call self%calc(state, lambda_0)
        call self%calc_lambda_dispersivity(lambda_0, htc_water, qx, qy, qz, lambda)

    end subroutine calc_thc_dispersivity_gp_3phase

    module subroutine calc_thc_gp_4phase(self, state, lambda, tangent)
        implicit none
        class(type_thc_4phase), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: lambda
        real(real64), intent(inout), optional :: tangent(2)

        integer(int32), parameter :: CAMPBELL_HANSSON_PARAMS = 6
        real(real64) :: phi1, phi2, phi3, phi4
        real(real64) :: dwater(2), dice(2)
        logical :: has_campbell_hansson

        call self%get_phi(state, phi1, phi2, phi3, phi4)
        if (present(tangent)) call get_content_tangents(state, dwater, dice)

        ! The mixing rule is a property of the material, chosen by whether it
        ! was configured with the Campbell-Hansson coefficients, and not of the
        ! state: that correlation is fitted over the whole range up to a full
        ! pore, so switching to the geometric mean when the gas fraction reaches
        ! zero would make lambda jump at saturation and at a frozen full pore.
        has_campbell_hansson = .false.
        if (allocated(self%params)) has_campbell_hansson = (size(self%params) >= CAMPBELL_HANSSON_PARAMS)

        if (has_campbell_hansson) then
            ! phi2 = liquid water, phi3 = ice. Hansson et al. (2004) count only
            ! the liquid in theta: ice enters through F*theta_i. The gas-filled
            ! fraction phi4 must not be passed: phi2 + phi3 + phi4 is the
            ! porosity, and with it the conductivity stops responding to drying.
            associate (c => self%params)
                if (present(tangent)) then
                    call evaluate_thermal_conductivity_campbell_hansson_tangent(c(1), c(2), c(3), c(4), c(5), c(6), &
                                                                                phi2, dwater, phi3, dice, lambda, tangent)
                else
                    call evaluate_thermal_conductivity_campbell_hansson(c(1), c(2), c(3), c(4), c(5), c(6), &
                                                                        phi2, phi3, lambda)
                end if
            end associate
        else if (present(tangent)) then
            ! phi4 = 1 - phi1 - phi2 - phi3, so the gas takes what the liquid and ice give up.
            call evaluate_thermal_conductivity_geometric4_tangent(self%material1, self%material2, self%material3, &
                                                                  self%material4, phi1, no_tangent, phi2, dwater, &
                                                                  phi3, dice, phi4, -(dwater + dice), lambda, tangent)
        else
            call evaluate_thermal_conductivity_geometric4(self%material1, self%material2, self%material3, &
                                                          self%material4, phi1, phi2, phi3, phi4, lambda)
        end if

    end subroutine calc_thc_gp_4phase

    module subroutine calc_thc_dispersivity_gp_4phase(self, state, lambda)
        implicit none
        class(type_thc_4phase), intent(in) :: self
        type(type_state), intent(in) :: state
        type(type_state_thc), intent(inout) :: lambda

        real(real64) :: lambda_0
        type(type_coordinate_dp), pointer :: water_flux

        real(real64) :: qx, qy, qz
        real(real64) :: htc_water, rho_water, cp_water

        call self%calc_water_density(state, rho_water)
        call self%calc_water_cp(state, cp_water)
        htc_water = rho_water * cp_water
        call state%get(water_flux=water_flux)
        if (associated(water_flux)) then
            qx = water_flux%x
            qy = water_flux%y
            qz = water_flux%z
        else
            qx = 0.0d0
            qy = 0.0d0
            qz = 0.0d0
        end if

        call self%calc(state, lambda_0)

        call self%calc_lambda_dispersivity(lambda_0, htc_water, qx, qy, qz, lambda)

    end subroutine calc_thc_dispersivity_gp_4phase

end submodule thermal_conductivity_multiphase
