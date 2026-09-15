submodule(physics_governing_hydraulic) hydraulic_coefficients
    implicit none
contains

    !> @brief Calculate Diffusion Term D_HH (Liquid Hydraulic Conductivity Tensor)
    !> @details
    !>   Mixed water-content form uses volumetric Darcy flux:
    !>   q = -K/(rho_w g) * grad P ...
    !>   D_HH = K_flh / (rho_w g)
    !>   The vapour conductivity is a separate tensor: it is driven by the gas
    !>   pressure, not by the liquid potential. See compute_diffusion_term_gas.
    module subroutine compute_diffusion_term_hydraulic(self, material_id, state, D_HH)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(inout) :: state
        real(real64), intent(inout) :: D_HH(:, :)

        real(real64) :: K_flh, rho_w
        real(real64) :: coeff_D

        integer(int32) :: i

        call self%physics%calc_Kflh(material_id, state, K_flh)
        call self%physics%calc_density_water(state, rho_w)
        coeff_D = K_flh / (rho_w * g)

        D_HH(:, :) = 0.0d0
        do i = 1, self%computation_dimension
            D_HH(i, i) = coeff_D
        end do

    end subroutine compute_diffusion_term_hydraulic

    !> @brief Calculate the vapour conductivity tensor D_HH_gas
    !> @details
    !>   D_HH_gas = K_vP / (rho_w g), zero when vapour transport is disabled.
    module subroutine compute_diffusion_term_gas_hydraulic(self, material_id, state, D_HH_gas)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(inout) :: state
        real(real64), intent(inout) :: D_HH_gas(:, :)

        real(real64) :: K_vP, rho_w
        real(real64) :: coeff_D

        integer(int32) :: i

        call self%calc_K_vP(material_id, state, K_vP)
        call self%physics%calc_density_water(state, rho_w)
        coeff_D = K_vP / (rho_w * g)

        D_HH_gas(:, :) = 0.0d0
        do i = 1, self%computation_dimension
            D_HH_gas(i, i) = coeff_D
        end do

    end subroutine compute_diffusion_term_gas_hydraulic

    !> @brief Calculate Advective (Gravity) Term V_H
    !> @details
    !>   q_grav = V_H
    !>   q_grav = -K * grad z
    !>   V_H = -K * grad z
    module subroutine compute_advective_term_hydraulic(self, material_id, state, V_H)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(inout) :: state
        real(real64), intent(inout) :: V_H(:)

        real(real64) :: K_flh
        real(real64) :: grav_flux_mag

        call self%physics%calc_Kflh(material_id, state, K_flh)

        ! Gravity volumetric flux magnitude in the mixed water-content equation.
        grav_flux_mag = K_flh

        V_H(:) = 0.0d0

        ! Assuming z is the vertical coordinate acting against gravity
        select case (self%computation_type)
        case (COMP_TYPES%XZ_2D%ID)
            V_H(2) = -grav_flux_mag ! z-direction
        case (COMP_TYPES%XYZ_3D%ID)
            V_H(3) = -grav_flux_mag ! z-direction
        case (COMP_TYPES%XY_2D%ID)
            ! No gravity contribution in XY_2D.
            ! Keep the thermo-osmotic term accumulated above.
        end select

    end subroutine compute_advective_term_hydraulic

    ! ==========================================================================
    ! Coupling Coefficients
    ! ==========================================================================

    !> @brief Calculate Coupling Mass Term C_HT = dTheta/dT + vapor storage tangent.
    module subroutine compute_coupling_mass_term_hydraulic(self, material_id, state, C_HT)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: C_HT

        real(real64) :: dTheta_dT, dQv_dT_term

        dTheta_dT = 0.0d0
        dQv_dT_term = 0.0d0

        call state%dTheta_dT%get(dTheta_dT)
        if (self%enable_vapor_transport) then
            call state%dQv_dT%get(dQv_dT_term)
        end if

        ! Modified Picard freezes this physical storage tangent at the current
        ! iterate while solving the coupled temperature-pressure increment.
        ! dTheta_dT is evaluated directly by calc_phase_split (see
        ! phase_systems.F90), not recombined from dQw_dT + (rho_i/rho_w)*dQi_dT
        ! + Qi*d(rho_i/rho_w)/dT: that recombination adds a ratio-derivative
        ! term calc_phase_split's own derivatives do not carry.
        C_HT = dTheta_dT + dQv_dT_term

    end subroutine compute_coupling_mass_term_hydraulic

    !> @brief Calculate Coupling Diffusion Term D_HT
    !> @details
    !>   Total moisture flux driven by temperature gradient:
    !>   \( q^{(T)} = -(K_{lT} + K_{vT}) \nabla T \).
    !>   Clapeyron equilibrium enters phase storage; the solved pore-water pressure
    !>   carries its hydraulic effect through the ordinary Darcy pressure gradient.
    module subroutine compute_coupling_diffusion_term_hydraulic(self, material_id, state, D_HT)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(inout) :: state
        real(real64), intent(inout) :: D_HT(:, :)

        real(real64) :: K_lT, K_vT
        real(real64) :: coeff_D
        integer(int32) :: i

        call self%physics%calc_KlT(material_id, state, K_lT)
        call self%calc_K_vT(material_id, state, K_vT)
        ! Freezing changes storage and therefore the solved pore-water pressure.
        ! Its Clapeyron gradient is not added again as an explicit Darcy force.
        coeff_D = K_lT + K_vT

        D_HT(:, :) = 0.0d0
        do i = 1, self%computation_dimension
            D_HT(i, i) = coeff_D
        end do

    end subroutine compute_coupling_diffusion_term_hydraulic

    ! ==========================================================================
    ! Helper Wrappers (Existing)
    ! ==========================================================================
    module subroutine calc_cryo_head_dT_hydraulic(self, material_id, state, dh_dT)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: dh_dT

        call self%physics%calc_cryo_head_dT(material_id, state, dh_dT)
    end subroutine calc_cryo_head_dT_hydraulic

    module subroutine calc_K_wT_hydraulic(self, target_id, state, K_wT)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: target_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: K_wT
        call self%physics%calc_KlT(target_id, state, K_wT)
    end subroutine calc_K_wT_hydraulic

    module subroutine calc_K_wP_hydraulic(self, target_id, state, K_wP)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: target_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: K_wP
        call self%physics%calc_Kflh(target_id, state, K_wP)
    end subroutine calc_K_wP_hydraulic

    module subroutine calc_K_vT_hydraulic(self, target_id, state, K_vT)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: target_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: K_vT
        if (self%enable_vapor_transport) then
            call self%physics%calc_KvT(target_id, state, K_vT)
        else
            K_vT = 0.0d0
        end if
    end subroutine calc_K_vT_hydraulic

    module subroutine calc_K_vP_hydraulic(self, target_id, state, K_vP)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: target_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: K_vP
        if (self%enable_vapor_transport) then
            call self%physics%calc_Kvh(target_id, state, K_vP)
        else
            K_vP = 0.0d0
        end if
    end subroutine calc_K_vP_hydraulic

    module subroutine update_water_phases_hydraulic(self, material_id, state)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(inout) :: state
        call self%physics%update_water_phases(material_id, state)
    end subroutine update_water_phases_hydraulic

    module subroutine calc_effective_density_hydraulic(self, material_id, state, bdf_coeffs, drho_dt)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(in) :: bdf_coeffs(:)
        real(real64), intent(inout) :: drho_dt

        type(type_state) :: local_state
        real(real64), pointer, dimension(:), contiguous :: temperature_history
        real(real64), pointer, dimension(:), contiguous :: pressure_history
        real(real64), pointer, dimension(:), contiguous :: ice_content_history

        real(real64) :: Qw, Qi, Qv
        real(real64) :: rho_w, rho_i
        real(real64) :: Uj, compressive_storage
        integer(int32) :: j, n

        nullify (temperature_history)
        nullify (pressure_history)
        nullify (ice_content_history)

        call state%temperature_history%get(temperature_history)
        call state%pressure_history%get(pressure_history)
        call state%ice_content_history%get(ice_content_history)

        drho_dt = 0.0d0
        if (.not. associated(temperature_history)) return
        if (.not. associated(pressure_history)) return
        if (.not. associated(ice_content_history)) return

        n = min(size(bdf_coeffs), size(temperature_history), size(pressure_history), size(ice_content_history))
        call local_state%copy(state)

        do j = 1, n
            call local_state%temperature%set(temperature_history(j))
            call local_state%pressure%set(pressure_history(j))
            call local_state%ice_content%set(ice_content_history(j))

            call self%update_water_phases(material_id, local_state)

            call local_state%water_content%get(Qw)
            call local_state%ice_content%get(Qi)
            call local_state%vapor_content%get(Qv)

            call self%physics%calc_density_water(local_state, rho_w)
            call self%physics%calc_density_ice(local_state, rho_i)

            ! Effective Density (Mass per unit volume)
            Uj = rho_w * Qw &
                 + rho_i * Qi &
                 + rho_w * Qv
            call self%compute_compressive_storage(material_id, local_state, compressive_storage)
            Uj = Uj + rho_std * compressive_storage

            drho_dt = drho_dt + bdf_coeffs(j) * Uj
        end do

    end subroutine calc_effective_density_hydraulic

    !> Evaluate the pore-water effective density rho_eff at the supplied state.
    !>
    !> Mathematical definition:
    !> \( \rho_{eff} = \rho_w \theta_w + \rho_{ice} \theta_{ice} + \rho_w \theta_v^{\star} \) [kg/m3]
    !>
    !> This is the conserved storage quantity of the water-mass balance and the plain
    !> counterpart of calc_effective_density (which returns its BDF time-derivative).
    !> Assumptions: the phase contents Qw, Qi, Qv stored in state are already
    !> consistent with (T, p_w) (call update_water_phases beforehand). Used by the
    !> conserved-quantity convergence norm (PDF 6.2.4). Cost: O(1).
    module subroutine calc_effective_density_value_hydraulic(self, material_id, state, rho_eff)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: rho_eff

        real(real64) :: Qw, Qi, Qv
        real(real64) :: rho_w, rho_i
        real(real64) :: compressive_storage

        call state%water_content%get(Qw)
        call state%ice_content%get(Qi)
        call state%vapor_content%get(Qv)

        call self%physics%calc_density_water(state, rho_w)
        call self%physics%calc_density_ice(state, rho_i)

        call self%compute_compressive_storage(material_id, state, compressive_storage)
        rho_eff = rho_w * Qw + rho_i * Qi + rho_w * Qv + rho_std * compressive_storage
    end subroutine calc_effective_density_value_hydraulic

    !> @brief Compute equivalent specific moisture capacity
    !> C_eq = dTheta/dP + vapor storage tangent + compressible-storage capacity.
    module subroutine compute_C_eq_hydraulic(self, material_id, state, C_eq)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: C_eq

        real(real64) :: dTheta_dP, dQv_dP_term
        real(real64) :: compressive_capacity, compressive_storage

        dTheta_dP = 0.0d0
        dQv_dP_term = 0.0d0

        call state%dTheta_dP%get(dTheta_dP)
        if (self%enable_vapor_transport) then
            call state%dQv_dP%get(dQv_dP_term)
        end if

        ! C_eq = dTheta/dP, evaluated directly by calc_phase_split (see
        ! phase_systems.F90) rather than recombined from
        ! dQw_dP + (rho_i/rho_w)*dQi_dP + Qi*d(rho_i/rho_w)/dP: calc_phase_split
        ! treats the density ratio as frozen in its own derivatives, so the
        ! ratio-derivative term is spurious surplus, not part of dTheta/dP.
        ! No post-hoc max(0,.) clamp: dTheta/dP is the analytic tangent of the
        ! (unclamped) values the residual integrates, and clamping it would
        ! break FD-consistency of the assembled Jacobian.
        call self%compute_compressive_storage(material_id, state, compressive_storage, compressive_capacity)
        C_eq = dTheta_dP + dQv_dP_term + compressive_capacity
    end subroutine compute_C_eq_hydraulic

    !> Compute the nonlinear-iteration capacity of the pressure block.
    module subroutine compute_iteration_capacity_hydraulic(self, material_id, state, capacity)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: capacity

        ! Modified Picard requires the mass matrix to be the pressure derivative
        ! of the storage term the residual integrates, i.e. C_eq = dTheta/dp at
        ! the current iterate. This used to return
        !   max(C_eq, iteration_capacity_bound)
        ! where the bound is the material-wide L-scheme constant
        ! L = max_h(dtheta/dh) / (rho_std g). Replacing the tangent by a fixed
        ! upper bound turns the pressure block into an L-scheme, whose
        ! contraction factor is 1 - C_eq/L:
        !   Mizoguchi VG curve (theta_s 0.535, theta_r 0.05, alpha 1.11, n 1.48,
        !   m 0.2) gives L = 8.1e-6 1/Pa, while the column's initial state
        !   (p_w = -5.41e4 Pa, |h| = 5.51 m) has C_eq = 1.4e-6 1/Pa, so
        !   L/C_eq = 5.6 and 1 - C_eq/L = 0.82 - the measured per-iteration ratio
        !   of the hydraulic residual. Where the pore volume binds C_eq collapses
        !   and the ratio reaches 1e4, so the tolerance is unreachable within any
        !   iteration budget. An L-scheme step also reduces only the L-weighted
        !   norm, never ||R||_2, so it is not a descent direction for the line
        !   search either.
        ! The tangent is small but not singular where the pore volume binds: that
        ! state is saturated, so the specific-storage term inside compute_C_eq
        ! carries the remaining capacity.
        call self%compute_C_eq(material_id, state, capacity)
    end subroutine compute_iteration_capacity_hydraulic

    !> @brief Compute BDF approximation of dTheta/dt for Mixed formulation.
    module subroutine compute_transient_term_mixed_hydraulic(self, material_id, state, bdf_coeffs, dTheta_dt)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(in) :: bdf_coeffs(:)
        real(real64), intent(inout) :: dTheta_dt

        type(type_state) :: local_state
        real(real64), pointer, dimension(:), contiguous :: temperature_history
        real(real64), pointer, dimension(:), contiguous :: pressure_history
        real(real64), pointer, dimension(:), contiguous :: ice_content_history

        real(real64) :: Theta_j
        integer(int32) :: j, n

        nullify (temperature_history)
        nullify (pressure_history)
        nullify (ice_content_history)

        call state%temperature_history%get(temperature_history)
        call state%pressure_history%get(pressure_history)
        call state%ice_content_history%get(ice_content_history)

        dTheta_dt = 0.0d0
        if (.not. associated(temperature_history)) return
        if (.not. associated(pressure_history)) return
        if (.not. associated(ice_content_history)) return

        n = min(size(bdf_coeffs), size(temperature_history), size(pressure_history), size(ice_content_history))
        if (n < 1) return

        ! j = 1 is the current state, and must be read from state itself rather
        ! than from history(1): the finite-difference tangent perturbs
        ! state%temperature/%pressure without touching the history arrays, so
        ! taking it from history(1) drops the current storage term from K_HH/K_HT.
        call self%compute_total_storage(material_id, state, Theta_j)
        dTheta_dt = dTheta_dt + bdf_coeffs(1) * Theta_j

        ! j >= 2 are previous-step states: recompute phases at the historical T/P.
        if (n >= 2) then
            call local_state%copy(state)
            ! Known levels carry no continuation; see compute_transient_term_thermal.
            call local_state%continuation_lambda%set(1.0d0)
            do j = 2, n
                call local_state%temperature%set(temperature_history(j))
                call local_state%pressure%set(pressure_history(j))
                call local_state%ice_content%set(ice_content_history(j))

                call self%update_water_phases(material_id, local_state)
                call self%compute_total_storage(material_id, local_state, Theta_j)

                dTheta_dt = dTheta_dt + bdf_coeffs(j) * Theta_j
            end do
        end if

    end subroutine compute_transient_term_mixed_hydraulic

    module subroutine compute_total_storage_hydraulic(self, material_id, state, Theta)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: Theta

        real(real64) :: Qw, Qi, Qv
        real(real64) :: rho_w, rho_i
        real(real64) :: compressive_storage

        call state%water_content%get(Qw)
        call state%ice_content%get(Qi)
        call state%vapor_content%get(Qv)

        call self%physics%calc_density_water(state, rho_w)
        call self%physics%calc_density_ice(state, rho_i)

        ! Theta = Qw + (rho_i/rho_w)*Qi + Qv  (rho_v ~ rho_w)
        Theta = Qw + (rho_i / rho_w) * Qi + Qv
        call self%compute_compressive_storage(material_id, state, compressive_storage)
        Theta = Theta + compressive_storage
    end subroutine compute_total_storage_hydraulic

    module subroutine compute_compressive_storage_hydraulic(self, material_id, state, storage, capacity)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: storage
        real(real64), intent(inout), optional :: capacity

        real(real64) :: pressure, saturation_pressure, storage_coefficient
        logical :: is_saturated

        storage = 0.0d0
        if (present(capacity)) capacity = 0.0d0
        if (.not. allocated(self%specific_storage)) return
        if (material_id < 1 .or. material_id > size(self%specific_storage)) return
        if (self%specific_storage(material_id) <= 0.0d0) return

        saturation_pressure = 0.0d0
        is_saturated = .false.
        call self%physics%calc_saturation_pressure(material_id, state, saturation_pressure, is_saturated)
        if (.not. is_saturated) return

        call state%pressure%get(pressure)
        storage_coefficient = self%specific_storage(material_id) / (rho_std * g)
        storage = storage_coefficient * max(0.0d0, pressure - saturation_pressure)
        if (present(capacity)) capacity = storage_coefficient
    end subroutine compute_compressive_storage_hydraulic

    !> @brief Compute segregation sink from temperature gradient magnitude.
    !>
    !> Returns an **effective** S_seg (volumetric liquid-water removal rate
    !> [1/s]) with available-water and open-pore capacity folded into the rate.
    module subroutine calc_segregation_sink_hydraulic(self, material_id, state, dt, S_seg)
        implicit none
        class(type_hydraulic), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(in) :: dt
        real(real64), intent(inout) :: S_seg

        type(type_coordinate_dp), pointer :: grad_T
        real(real64) :: grad_T_mag

        S_seg = 0.0d0
        nullify (grad_T)
        call state%grad_T%get(grad_T)
        if (.not. associated(grad_T)) return

        grad_T_mag = sqrt(grad_T%x**2 + grad_T%y**2 + grad_T%z**2)
        if (grad_T_mag <= 0.0d0) return

        call self%physics%calc_effective_segregation_sink(material_id, state, grad_T_mag, dt, S_seg)

    end subroutine calc_segregation_sink_hydraulic

end submodule hydraulic_coefficients
