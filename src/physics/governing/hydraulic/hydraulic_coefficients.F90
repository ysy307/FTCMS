submodule(physics_governing_hydraulic) hydraulic_coefficients
    use :: constitutive_constants, only:tangent_pressure, tangent_temperature
    use :: constitutive_laws, only:evaluate_total_stored_water, evaluate_total_stored_water_tangent
    implicit none
contains

    !> @brief Calculate Diffusion Term D_HH (Liquid Hydraulic Conductivity Tensor)
    !> @details
    !>   Mixed water-content form uses volumetric Darcy flux:
    !>   q = -K/(rho_w g) * grad P ...
    !>   D_HH = K_flh / (rho_w g)
    !>   The vapour conductivity is a separate tensor: it is driven by the gas
    !>   pressure, not by the liquid potential. See compute_diffusion_term_gas.
    module subroutine compute_diffusion_term_hydraulic(self, material_id, state, D_HH, log_conductivity_tangent)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(inout) :: state
        real(real64), intent(inout) :: D_HH(:, :)
        real(real64), intent(inout), optional :: log_conductivity_tangent(2)

        real(real64) :: K_flh, rho_w
        real(real64) :: coeff_D

        integer(int32) :: i

        call self%physics%calc_Kflh(material_id, state, K_flh, log_conductivity_tangent)
        call self%physics%calc_density_water(state, rho_w)
        coeff_D = K_flh / (rho_w * g)

        D_HH(:, :) = 0.0d0
        do i = 1, self%computation_dimension
            D_HH(i, i) = coeff_D
        end do

    end subroutine compute_diffusion_term_hydraulic

    ! ==========================================================================
    ! Coupling Coefficients
    ! ==========================================================================

    module subroutine update_water_phases_hydraulic(self, material_id, state)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(inout) :: state
        call self%physics%update_water_phases(material_id, state)
    end subroutine update_water_phases_hydraulic


    module subroutine compute_total_storage_hydraulic(self, material_id, state, Theta, tangent)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: Theta
        real(real64), intent(inout), optional :: tangent(2)

        real(real64) :: Qw, Qi, Qv
        real(real64) :: rho_w, rho_i
        real(real64) :: compressive_storage
        real(real64) :: dQw(2), dQi(2), dQv(2), drho_w(2), drho_i(2), dstorage(2)

        call state%water_content%get(Qw)
        call state%ice_content%get(Qi)
        call state%vapor_content%get(Qv)

        ! Theta = Qw + (rho_i/rho_w)*Qi + Qv + S_c  (rho_v ~ rho_w)
        if (present(tangent)) then
            call state%dQw_dP%get(dQw(tangent_pressure))
            call state%dQw_dT%get(dQw(tangent_temperature))
            call state%dQi_dP%get(dQi(tangent_pressure))
            call state%dQi_dT%get(dQi(tangent_temperature))
            call state%dQv_dP%get(dQv(tangent_pressure))
            call state%dQv_dT%get(dQv(tangent_temperature))
            call self%physics%calc_density_water(state, rho_w, drho_w)
            call self%physics%calc_density_ice(state, rho_i, drho_i)
            call self%compute_compressive_storage(material_id, state, compressive_storage, dstorage)
            call evaluate_total_stored_water_tangent(Qw, dQw, Qi, dQi, Qv, dQv, rho_w, drho_w, rho_i, drho_i, &
                                                     compressive_storage, dstorage, Theta, tangent)
        else
            call self%physics%calc_density_water(state, rho_w)
            call self%physics%calc_density_ice(state, rho_i)
            call self%compute_compressive_storage(material_id, state, compressive_storage)
            call evaluate_total_stored_water(Qw, Qi, Qv, rho_w, rho_i, compressive_storage, Theta)
        end if
    end subroutine compute_total_storage_hydraulic

    module subroutine compute_compressive_storage_hydraulic(self, material_id, state, storage, tangent)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: storage
        real(real64), intent(inout), optional :: tangent(2)

        real(real64) :: pressure, saturation_pressure, storage_coefficient
        real(real64) :: dsaturation(2)
        logical :: is_saturated

        storage = 0.0d0
        if (present(tangent)) tangent(:) = 0.0d0
        if (.not. allocated(self%specific_storage)) return
        if (material_id < 1 .or. material_id > size(self%specific_storage)) return
        if (self%specific_storage(material_id) <= 0.0d0) return

        saturation_pressure = 0.0d0
        is_saturated = .false.
        if (present(tangent)) then
            call self%physics%calc_saturation_pressure(material_id, state, saturation_pressure, is_saturated, dsaturation)
        else
            call self%physics%calc_saturation_pressure(material_id, state, saturation_pressure, is_saturated)
        end if
        if (.not. is_saturated) return

        call state%pressure%get(pressure)
        storage_coefficient = self%specific_storage(material_id) / (rho_std * g)
        storage = storage_coefficient * max(0.0d0, pressure - saturation_pressure)
        ! The saturation pressure moves with the state through the ice it holds.
        if (present(tangent) .and. pressure > saturation_pressure) then
            tangent(tangent_pressure) = storage_coefficient * (1.0d0 - dsaturation(tangent_pressure))
            tangent(tangent_temperature) = -storage_coefficient * dsaturation(tangent_temperature)
        end if
    end subroutine compute_compressive_storage_hydraulic

end submodule hydraulic_coefficients
