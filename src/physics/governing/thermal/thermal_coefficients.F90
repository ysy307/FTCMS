!> Implementation overview
!>
!> Algorithm:
!> - Calculates thermal coefficients and governing equation terms
!> - Assembles enthalpy, transient, mass, diffusion, advection, and latent terms
submodule(physics_governing_thermal) thermal_coefficients
    implicit none

contains

    !> Calculate enthalpy (internal energy) density per unit volume
    !>
    !> Mathematical definition:
    !> - Equation of state: \( U = f(T, \phi, \dots) \)
    module subroutine calc_enthalpy_density_thermal(self, material_id, state, U)
        implicit none
        !> Thermal governing equation object
        class(type_thermal), intent(in) :: self
        !> Material identifier
        integer(int32), intent(in) :: material_id
        !> Thermodynamic state
        type(type_state), intent(in) :: state
        !> Enthalpy density \( U \) [J/m3]
        !> Overwritten on exit
        real(real64), intent(inout) :: U

        ! Local variables
        real(real64) :: temperature
        real(real64) :: porosity, Qw, Qi, Qv
        real(real64) :: rho_s, rho_w, rho_i
        real(real64) :: c_s, c_w, c_i, c_v
        real(real64) :: Lf, Lv
        logical :: has_rho_w

        call state%temperature%get(temperature)
        call state%porosity%get(porosity)
        call state%water_content%get(Qw)
        call state%ice_content%get(Qi)
        call state%vapor_content%get(Qv)

        U = 0.0d0
        has_rho_w = .false.

        ! Solid phase
        if (porosity > 0.0d0) then
            call self%physics%get_density_solid(material_id, rho_s)
            call self%physics%get_specific_heat_solid(material_id, c_s)
            U = U + rho_s * c_s * (1.0d0 - porosity) * temperature
        end if

        ! Water phase
        if (Qw > 0.0d0) then
            call self%physics%calc_density_water(state, rho_w)
            has_rho_w = .true.
            call self%physics%calc_specific_heat_water(state, c_w)
            U = U + rho_w * c_w * Qw * temperature
        end if

        ! Ice phase (Include Latent Heat of Fusion)
        ! Additional segregation-sink latent heat is handled by an explicit source term.
        if (Qi > 0.0d0) then
            call self%physics%calc_density_ice(state, rho_i)
            call self%physics%calc_specific_heat_ice(state, c_i)
            call self%physics%calc_latent_heat_fusion(material_id, state, Lf)
            U = U + rho_i * c_i * Qi * temperature - Lf * rho_i * Qi
        end if

        ! Vapor phase (Include Latent Heat of Vaporization)
        if (Qv > 0.0d0) then
            if (.not. has_rho_w) then
                ! Reference density for vapor is typically liquid phase
                call self%physics%calc_density_water(state, rho_w)
                has_rho_w = .true.
            end if
            call self%physics%calc_specific_heat_vapor(state, c_v)
            call self%physics%calc_latent_heat_vaporization(material_id, state, Lv)
            U = U + rho_w * c_v * Qv * temperature + Lv * rho_w * Qv
        end if

    end subroutine calc_enthalpy_density_thermal

    ! --------------------------------------------------------------------------
    ! Helper Wrappers
    ! --------------------------------------------------------------------------

    !> Wrapper to calculate water density
    module subroutine calc_density_water_thermal(self, state, rho_water)
        implicit none
        !> Thermal governing equation object
        class(type_thermal), intent(in) :: self
        !> Thermodynamic state
        type(type_state), intent(in) :: state
        !> Water density
        !> Overwritten on exit
        real(real64), intent(inout) :: rho_water

        call self%physics%calc_density_water(state, rho_water)
    end subroutine calc_density_water_thermal

    !> Wrapper to calculate ice density
    module subroutine calc_density_ice_thermal(self, state, rho_ice)
        implicit none
        !> Thermal governing equation object
        class(type_thermal), intent(in) :: self
        !> Thermodynamic state
        type(type_state), intent(in) :: state
        !> Ice density
        !> Overwritten on exit
        real(real64), intent(inout) :: rho_ice

        call self%physics%calc_density_ice(state, rho_ice)
    end subroutine calc_density_ice_thermal

    !> Wrapper to calculate saturation vapor density
    module subroutine calc_density_vapor_saturation_thermal(self, state, rho_vapor_sat)
        implicit none
        !> Thermal governing equation object
        class(type_thermal), intent(in) :: self
        !> Thermodynamic state
        type(type_state), intent(in) :: state
        !> Saturation vapor density
        !> Overwritten on exit
        real(real64), intent(inout) :: rho_vapor_sat

        call self%physics%calc_density_vapor_saturation(state, rho_vapor_sat)
    end subroutine calc_density_vapor_saturation_thermal

    !> Wrapper to update water phases
    module subroutine update_water_phases_thermal(self, material_id, state)
        implicit none
        !> Thermal governing equation object
        class(type_thermal), intent(in) :: self
        !> Material identifier
        integer(int32), intent(in) :: material_id
        !> Thermodynamic state
        !> Overwritten on exit
        type(type_state), intent(inout) :: state

        call self%physics%update_water_phases(material_id, state)
    end subroutine update_water_phases_thermal

    module subroutine project_ice_content_thermal(self, material_id, state, projected_ice, &
                                                   ice_increment, equilibrium_error, active_bound, &
                                                   dice_dT, dice_dP)
        implicit none
        class(type_thermal), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: projected_ice
        real(real64), intent(inout) :: ice_increment
        real(real64), intent(inout) :: equilibrium_error
        integer(int32), intent(inout), optional :: active_bound
        real(real64), intent(inout), optional :: dice_dT
        real(real64), intent(inout), optional :: dice_dP

        call self%physics%project_ice_content(material_id, state, projected_ice, ice_increment, equilibrium_error, &
                                               active_bound, dice_dT, dice_dP)
    end subroutine project_ice_content_thermal

    module subroutine calc_conserved_target_thermal(self, material_id, state, target_total_water, available)
        implicit none
        class(type_thermal), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: target_total_water
        logical, intent(inout) :: available

        call self%physics%calc_conserved_target(material_id, state, target_total_water, available)
    end subroutine calc_conserved_target_thermal

    module subroutine solve_local_conserved_equilibrium_thermal(self, material_id, state, target_total_water, &
                                                                 new_pressure, new_ice, converged)
        implicit none
        class(type_thermal), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(in) :: target_total_water
        real(real64), intent(inout) :: new_pressure
        real(real64), intent(inout) :: new_ice
        logical, intent(inout) :: converged

        call self%physics%solve_local_conserved_equilibrium(material_id, state, target_total_water, new_pressure, &
                                                              new_ice, converged)
    end subroutine solve_local_conserved_equilibrium_thermal

    ! ==========================================================================
    ! Transient / Mass Terms (Refactored for NR & Picard)
    ! ==========================================================================

    !> Calculate time derivative term (Residual): \( \frac{dU}{dt} \)
    !>
    !> Mathematical definition:
    !> - Used for residual calculation in Newton-Raphson method
    !> - \( Res_{transient} = \frac{dU}{dt} = \sum bdf\_coeffs_j U(t_{n+1-j}) \)
    !>
    !> When enthalpy cache is valid, history terms (j >= 2) use pre-computed
    !> nodal U values instead of recomputing phase-change. The cache is
    !> populated once per time step via cache_enthalpy_history.
    module subroutine compute_transient_term_thermal(self, material_id, state, bdf_coeffs, dU_dt)
        implicit none
        !> Thermal governing equation object
        class(type_thermal), intent(in) :: self
        !> Material identifier
        integer(int32), intent(in) :: material_id
        !> Thermodynamic state
        type(type_state), intent(in) :: state
        !> BDF coefficients for time integration
        real(real64), intent(in) :: bdf_coeffs(:)
        !> Time derivative of enthalpy density
        !> Overwritten on exit
        real(real64), intent(inout) :: dU_dt

        type(type_state) :: local_state
        real(real64), pointer, contiguous, dimension(:) :: temperature_history
        real(real64), pointer, contiguous, dimension(:) :: pressure_history
        real(real64), pointer, contiguous, dimension(:) :: porosity_history
        real(real64), pointer, contiguous, dimension(:) :: ice_content_history
        real(real64) :: Uj
        integer(int32) :: j, n_hist

        nullify (ice_content_history)
        call state%get(temperature_history=temperature_history, &
                       pressure_history=pressure_history, &
                       porosity_history=porosity_history)
        call state%ice_content_history%get(ice_content_history)

        dU_dt = 0.0d0
        if (.not. associated(temperature_history)) return
        if (.not. associated(pressure_history)) return
        if (.not. associated(porosity_history)) return
        if (.not. associated(ice_content_history)) return
        n_hist = min(size(bdf_coeffs), size(temperature_history), size(pressure_history), &
                     size(porosity_history), size(ice_content_history))

        ! Copy once: only T/P/porosity vary over the BDF history; the other state
        ! fields (e.g. clay fraction) are constant, so the deep copy is hoisted
        ! out of the history loop.
        ! j = 1 is the current state: phases were already computed during set_state,
        ! so reuse them directly and skip the redundant update_water_phases.
        call self%calc_enthalpy_density(material_id, state, Uj)
        dU_dt = dU_dt + bdf_coeffs(1) * Uj

        ! j >= 2 are previous-step states: recompute phases at the historical T/P.
        if (n_hist >= 2) then
            call local_state%copy(state)
            ! Known levels carry no continuation. The blend is already the
            ! identity at its own reference level, but level n-1 is not that
            ! reference, so U at n-1 would otherwise move with lambda.
            call local_state%continuation_lambda%set(1.0d0)
            do j = 2, n_hist
                call local_state%temperature%set(temperature_history(j))
                call local_state%pressure%set(pressure_history(j))
                call local_state%porosity%set(porosity_history(j))
                call local_state%ice_content%set(ice_content_history(j))
                call self%update_water_phases(material_id, local_state)
                call self%calc_enthalpy_density(material_id, local_state, Uj)

                dU_dt = dU_dt + bdf_coeffs(j) * Uj
            end do
        end if

    end subroutine compute_transient_term_thermal

    !> Calculate the apparent heat capacity \( C_{app} = dU/dT \)
    !>
    !> Mathematical definition:
    !> - The analytical tangent of the enthalpy density, phase by phase, with
    !>   the latent heat of fusion and of vaporization included through
    !>   \( dQ_i/dT \) and \( dQ_v/dT \)
    !>
    !> It reads only the current state, never the time history, so a caller
    !> that discretises the time derivative itself can use it directly.
    !> compute_mass_term_thermal calls it for the tangent scheme.
    module subroutine compute_apparent_capacity_thermal(self, material_id, state, C_app)
        implicit none
        !> Thermal governing equation object
        class(type_thermal), intent(in) :: self
        !> Material identifier
        integer(int32), intent(in) :: material_id
        !> Thermodynamic state
        type(type_state), intent(in) :: state
        !> Apparent heat capacity
        !> Overwritten on exit
        real(real64), intent(inout) :: C_app

        real(real64) :: temperature
        real(real64) :: porosity, Qw, Qi, Qv
        real(real64) :: rho_s, rho_w, rho_i
        real(real64) :: c_s, c_w, c_i, c_v
        real(real64) :: dQw_dT, dQi_dT, dQv_dT, Lf, Lv
        logical :: has_rho_w

        call state%temperature%get(temperature)

        C_app = 0.0d0
        has_rho_w = .false.

        ! --- Analytical tangent (apparent heat capacity) C_app = dU/dT (always) ---
        call state%get(porosity=porosity, water_content=Qw, &
                       ice_content=Qi, vapor_content=Qv)

        ! Solid
        if (porosity > 0.0d0) then
            call self%physics%get_density_solid(material_id, rho_s)
            call self%physics%get_specific_heat_solid(material_id, c_s)
            C_app = C_app + rho_s * c_s * (1.0d0 - porosity)
        end if

        ! Water
        if (Qw > 0.0d0) then
            call self%physics%calc_density_water(state, rho_w)
            has_rho_w = .true.
            call self%physics%calc_specific_heat_water(state, c_w)
            call state%dQw_dT%get(dQw_dT)
            C_app = C_app + rho_w * c_w * (Qw + temperature * dQw_dT)
        end if

        ! Ice (including latent heat derivative).
        ! Additional segregation-sink latent heat is handled by an explicit source.
        call state%dQi_dT%get(dQi_dT)
        if (Qi > 0.0d0 .or. dQi_dT /= 0.0d0) then
            call self%physics%calc_density_ice(state, rho_i)
            call self%physics%calc_specific_heat_ice(state, c_i)
            C_app = C_app + rho_i * c_i * (Qi + temperature * dQi_dT)
            call self%physics%calc_latent_heat_fusion(material_id, state, Lf)
            C_app = C_app - Lf * rho_i * dQi_dT
        end if

        ! Vapor (including latent heat derivative)
        if (Qv > 0.0d0) then
            if (.not. has_rho_w) then
                call self%physics%calc_density_water(state, rho_w)
                has_rho_w = .true.
            end if
            call self%physics%calc_specific_heat_vapor(state, c_v)
            call self%physics%calc_latent_heat_vaporization(material_id, state, Lv)
            call state%dQv_dT%get(dQv_dT)
            C_app = C_app + rho_w * c_v * (Qv + temperature * dQv_dT) + Lv * rho_w * dQv_dT
        end if
    end subroutine compute_apparent_capacity_thermal

    !> Calculate heat capacity coefficient \( C_{TT} = \frac{dU}{dT} \)
    !>
    !> Mathematical definition:
    !> - Used for Jacobian in NR method or mass matrix in Picard method
    !> - Applies Tangent scheme for small \( dT \) or SECANT scheme for large \( dT \)
    module subroutine compute_mass_term_thermal(self, material_id, state, C_TT, scheme_opt)
        implicit none
        !> Thermal governing equation object
        class(type_thermal), intent(in) :: self
        !> Material identifier
        integer(int32), intent(in) :: material_id
        !> Thermodynamic state
        type(type_state), intent(in) :: state
        !> Heat capacity coefficient
        !> Overwritten on exit
        real(real64), intent(inout) :: C_TT
        !> Optional scheme selector (TANGENT or SECANT)
        integer(int32), intent(in), optional :: scheme_opt

        real(real64) :: temperature
        real(real64), pointer, contiguous, dimension(:) :: temperature_history
        real(real64), pointer, contiguous, dimension(:) :: pressure_history
        real(real64), pointer, contiguous, dimension(:) :: porosity_history

        ! Secant variables
        real(real64) :: C_TT_current, C_TT_old, C_TT_secant, dT
        type(type_state) :: temp_state
        integer(int32) :: use_scheme
        ! Minimum |T^m - T^n| for evaluating the chord heat capacity. Below this the
        ! tangent is the correct slope and the chord would be ill-conditioned.
        real(real64), parameter :: DT_SECANT_FLOOR = 1.0d-3

        ! Determine default behavior
        call state%get(temperature=temperature, temperature_history=temperature_history, &
                       pressure_history=pressure_history, porosity_history=porosity_history)

        if (.not. associated(temperature_history)) then
            C_TT = 0.0d0
            return
        end if
        if (.not. associated(pressure_history)) then
            C_TT = 0.0d0
            return
        end if
        if (.not. associated(porosity_history)) then
            C_TT = 0.0d0
            return
        end if

        if (size(temperature_history) >= 2) then
            dT = temperature - temperature_history(2)
        else
            dT = 0.0d0
        end if

        call self%compute_apparent_capacity(material_id, state, C_TT)

        ! Modified Picard uses the physical tangent evaluated at the current
        ! iterate. A pure secant remains available only when explicitly
        ! requested by a caller.
        use_scheme = SCHEME_TANGENT
        if (present(scheme_opt)) use_scheme = scheme_opt

        if (use_scheme == SCHEME_SECANT) then
            if (size(temperature_history) >= 2 .and. abs(dT) > DT_SECANT_FLOOR) then
                call self%calc_enthalpy_density(material_id, state, C_TT_current)
                call temp_state%copy(state)
                call temp_state%temperature%set(temperature_history(2))
                ! The chord's old end is the reference level, where the blend is
                ! the identity for any lambda. Stated explicitly so a pressure-
                ! dependent freezing model cannot make the two ends disagree.
                call temp_state%continuation_lambda%set(1.0d0)
                call self%update_water_phases(material_id, temp_state)
                call self%calc_enthalpy_density(material_id, temp_state, C_TT_old)
                C_TT_secant = (C_TT_current - C_TT_old) / dT
                C_TT = C_TT_secant
            end if
        end if

    end subroutine compute_mass_term_thermal

    !> @brief Calculate Coupling Mass Term C_TH = dU/dP
    !> @details
    !>   U = rho_w*c_w*Qw*T + (rho_i*c_i*T - Lf*rho_i)*Qi + (rho_w*c_v*T + Lv*rho_w)*Qv + solid
    !>   C_TH = (rho_w*c_w*T)*dQw/dP + (rho_i*c_i*T - Lf*rho_i)*dQi/dP
    !>          + (rho_w*c_v*T + Lv*rho_w)*dQv/dP
    module subroutine compute_coupling_mass_term_thermal(self, material_id, state, C_TH)
        implicit none
        class(type_thermal), intent(in) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: C_TH

        real(real64) :: temperature
        real(real64) :: rho_w, rho_i
        real(real64) :: c_w, c_i, c_v
        real(real64) :: Lf, Lv
        real(real64) :: dQw_dP, dQi_dP, dQv_dP

        C_TH = 0.0d0

        call state%temperature%get(temperature)
        call state%dQw_dP%get(dQw_dP)
        call state%dQi_dP%get(dQi_dP)
        call state%dQv_dP%get(dQv_dP)

        if (abs(dQw_dP) + abs(dQi_dP) + abs(dQv_dP) < 1.0d-30) return

        call self%physics%calc_density_water(state, rho_w)
        call self%physics%calc_density_ice(state, rho_i)
        call self%physics%calc_specific_heat_water(state, c_w)
        call self%physics%calc_specific_heat_ice(state, c_i)
        call self%physics%calc_specific_heat_vapor(state, c_v)
        call self%physics%calc_latent_heat_fusion(material_id, state, Lf)
        call self%physics%calc_latent_heat_vaporization(material_id, state, Lv)

        C_TH = rho_w * c_w * temperature * dQw_dP &
               + (rho_i * c_i * temperature - Lf * rho_i) * dQi_dP &
               + (rho_w * c_v * temperature + Lv * rho_w) * dQv_dP

    end subroutine compute_coupling_mass_term_thermal

    ! ==========================================================================
    ! Flux Terms (Diffusion / Advection / Latent Source)
    ! ==========================================================================

    !> Calculate diffusion term matrix
    !>
    !> Mathematical definition:
    !> - Assembles thermal conductivity tensor based on spatial dimensions
    module subroutine compute_diffusion_term_thermal(self, material_id, state, D_TT)
        implicit none
        !> Thermal governing equation object
        class(type_thermal), intent(in) :: self
        !> Material identifier
        integer(int32), intent(in) :: material_id
        !> Thermodynamic state
        !> Overwritten on exit
        type(type_state), intent(inout) :: state
        !> Diffusion matrix components
        !> Overwritten on exit
        real(real64), intent(inout) :: D_TT(:, :)

        type(type_state_thc) :: lambda

        call self%physics%calc_thermal_conductivity(material_id, state, lambda)

        D_TT(:, :) = 0.0d0
        select case (self%computation_type)
        case (COMP_TYPES%XY_2D%ID)
            D_TT(1, 1) = lambda%lambda_xx
            D_TT(1, 2) = lambda%lambda_xy
            D_TT(2, 1) = lambda%lambda_xy
            D_TT(2, 2) = lambda%lambda_yy
        case (COMP_TYPES%XZ_2D%ID)
            D_TT(1, 1) = lambda%lambda_xx
            D_TT(1, 2) = lambda%lambda_zx
            D_TT(2, 1) = lambda%lambda_zx
            D_TT(2, 2) = lambda%lambda_zz
        case (COMP_TYPES%XYZ_3D%ID)
            D_TT(1, 1) = lambda%lambda_xx
            D_TT(1, 2) = lambda%lambda_xy
            D_TT(1, 3) = lambda%lambda_zx
            D_TT(2, 1) = lambda%lambda_xy
            D_TT(2, 2) = lambda%lambda_yy
            D_TT(2, 3) = lambda%lambda_yz
            D_TT(3, 1) = lambda%lambda_zx
            D_TT(3, 2) = lambda%lambda_yz
            D_TT(3, 3) = lambda%lambda_zz
        end select
    end subroutine compute_diffusion_term_thermal

    !> Calculate advective term vector
    !>
    !> Mathematical definition:
    !> - Computes thermal advection driven by water and vapor fluxes
    module subroutine compute_advective_term_thermal(self, material_id, state, V_TT)
        implicit none
        !> Thermal governing equation object
        class(type_thermal), intent(in) :: self
        !> Material identifier
        integer(int32), intent(in) :: material_id
        !> Thermodynamic state
        !> Overwritten on exit
        type(type_state), intent(inout) :: state
        !> Advective velocity vector
        !> Overwritten on exit
        real(real64), intent(inout) :: V_TT(:)

        type(type_coordinate_dp), pointer :: water_flux, vapor_flux
        real(real64) :: rho_w, c_w, c_v
        real(real64) :: flux_norm1

        call state%get(water_flux=water_flux, vapor_flux=vapor_flux)

        V_TT(:) = 0.0d0
        flux_norm1 = abs(water_flux%x) + abs(water_flux%y) + abs(water_flux%z) + &
                     abs(vapor_flux%x) + abs(vapor_flux%y) + abs(vapor_flux%z)
        if (flux_norm1 <= epsilon(1.0d0)) return

        call self%physics%calc_density_water(state, rho_w)
        call self%physics%calc_specific_heat_water(state, c_w)
        call self%physics%calc_specific_heat_vapor(state, c_v)

        select case (self%computation_type)
        case (COMP_TYPES%XY_2D%ID)
            V_TT(1) = rho_w * (c_w * water_flux%x + c_v * vapor_flux%x)
            V_TT(2) = rho_w * (c_w * water_flux%y + c_v * vapor_flux%y)
        case (COMP_TYPES%XZ_2D%ID)
            V_TT(1) = rho_w * (c_w * water_flux%x + c_v * vapor_flux%x)
            V_TT(2) = rho_w * (c_w * water_flux%z + c_v * vapor_flux%z)
        case (COMP_TYPES%XYZ_3D%ID)
            V_TT(1) = rho_w * (c_w * water_flux%x + c_v * vapor_flux%x)
            V_TT(2) = rho_w * (c_w * water_flux%y + c_v * vapor_flux%y)
            V_TT(3) = rho_w * (c_w * water_flux%z + c_v * vapor_flux%z)
        end select
    end subroutine compute_advective_term_thermal

    !> Calculate latent heat source term
    !>
    !> Mathematical definition:
    !> - \( L_{TT} = \rho_w L_v K_{vT} \)
    module subroutine compute_latent_term_thermal(self, material_id, state, L_TT)
        implicit none
        !> Thermal governing equation object
        class(type_thermal), intent(in) :: self
        !> Material identifier
        integer(int32), intent(in) :: material_id
        !> Thermodynamic state
        !> Overwritten on exit
        type(type_state), intent(inout) :: state
        !> Latent heat term
        !> Overwritten on exit
        real(real64), intent(inout) :: L_TT

        real(real64) :: rho_w, L_v, K_vT

        if (.not. self%enable_vapor_transport) then
            L_TT = 0.0d0
            return
        end if

        call self%physics%calc_density_water(state, rho_w)
        call self%physics%calc_latent_heat_vaporization(material_id, state, L_v)
        call self%physics%calc_KvT(material_id, state, K_vT)
        L_TT = rho_w * L_v * K_vT
    end subroutine compute_latent_term_thermal

    !> Calculate history term at integration points
    !>
    !> Mathematical definition:
    !> - Used for linearization in methods like Picard
    !> - \( History = C_{TT} \sum_{k=2}^{order} \alpha_k T_{n+1-k} \)
    !> - Warning: Do not use when constructing residuals using the enthalpy method with NR
    module subroutine compute_history_term_thermal(self, material_id, state, bdf_coeffs, history_term)
        implicit none
        !> Thermal governing equation object
        class(type_thermal), intent(in) :: self
        !> Material identifier
        integer(int32), intent(in) :: material_id
        !> Thermodynamic state
        type(type_state), intent(in) :: state
        !> BDF coefficients for time integration
        real(real64), intent(in) :: bdf_coeffs(:)
        !> History term value
        !> Overwritten on exit
        real(real64), intent(inout) :: history_term

        real(real64), pointer, contiguous, dimension(:) :: temperature_history
        real(real64) :: C_TT
        real(real64) :: T_hist_sum
        integer(int32) :: j, n_steps

        call state%get(temperature_history=temperature_history)

        ! Calculate heat capacity (Picard generally uses Secant at the current iteration temperature)
        call self%compute_mass_term(material_id, state, C_TT)

        T_hist_sum = 0.0d0
        n_steps = size(bdf_coeffs)

        ! Exclude alpha_0 (current) and sum only past terms
        do j = 2, n_steps
            if (j > size(temperature_history)) exit
            T_hist_sum = T_hist_sum + bdf_coeffs(j) * temperature_history(j)
        end do

        history_term = C_TT * T_hist_sum

    end subroutine compute_history_term_thermal

    ! ==========================================================================
    ! Enthalpy Cache Management
    ! ==========================================================================

    !> Pre-compute and cache nodal enthalpy density for BDF history levels
    !>
    !> Called once per time step before the nonlinear iteration begins.
    !> For each node and each history level (k = 1..num_hist), computes
    !> U(T_{n-k}, P_{n-k}, phi_{n-k}) and stores it in enthalpy_cache.
    !> This avoids redundant phase-change recalculation during assembly.
    !>
    !> NOTE (A1 Clapeyron closure): this cache is currently DEAD CODE - no
    !> call site exists anywhere in src/ or test/ (enthalpy_cache_valid is
    !> never set at runtime), so the prognostic-ice substitution applied to
    !> compute_transient_term_thermal / compute_mass_term_thermal is NOT
    !> replicated here. If this cache is ever revived, it must receive the
    !> nodal prognostic ice history (its interface only carries T/P/phi
    !> arrays) and apply the same per-level ice substitution, otherwise it
    !> would silently reintroduce the equilibrium-ice bias at constrained
    !> nodes.
    module subroutine cache_enthalpy_history_thermal(self, num_nodes, num_hist, material_ids, &
                                                     temperature_all, pressure_all, porosity_all)
        implicit none
        class(type_thermal), intent(inout) :: self
        !> Number of mesh nodes
        integer(int32), intent(in) :: num_nodes
        !> Number of history levels to cache
        integer(int32), intent(in) :: num_hist
        !> Material ID for each node (from element adjacency)
        integer(int32), intent(in) :: material_ids(:)
        !> Nodal temperature history: temperature_all(node, hist_level)
        real(real64), intent(in) :: temperature_all(:, :)
        !> Nodal pressure history: pressure_all(node, hist_level)
        real(real64), intent(in) :: pressure_all(:, :)
        !> Nodal porosity history: porosity_all(node, hist_level)
        real(real64), intent(in) :: porosity_all(:, :)

        type(type_state) :: local_state
        real(real64) :: Uj
        integer(int32) :: i, k

        ! Allocate or reallocate cache if dimensions changed
        if (.not. allocated(self%enthalpy_cache) .or. &
            self%enthalpy_cache_num_nodes /= num_nodes .or. &
            self%enthalpy_cache_num_hist /= num_hist) then
            if (allocated(self%enthalpy_cache)) deallocate (self%enthalpy_cache)
            allocate (self%enthalpy_cache(num_nodes, num_hist))
            self%enthalpy_cache_num_nodes = num_nodes
            self%enthalpy_cache_num_hist = num_hist
        end if

        ! Compute enthalpy at each node for each history level
        do k = 1, num_hist
            do i = 1, num_nodes
                call local_state%temperature%set(temperature_all(i, k))
                call local_state%pressure%set(pressure_all(i, k))
                call local_state%porosity%set(porosity_all(i, k))
                call self%update_water_phases(material_ids(i), local_state)
                call self%calc_enthalpy_density(material_ids(i), local_state, Uj)
                self%enthalpy_cache(i, k) = Uj
            end do
        end do

        self%enthalpy_cache_valid = .true.

    end subroutine cache_enthalpy_history_thermal

    !> Invalidate the enthalpy cache
    !>
    !> Called when time step is accepted and variable history is shifted,
    !> since the cached values correspond to the previous step layout.
    module subroutine invalidate_enthalpy_cache_thermal(self)
        implicit none
        class(type_thermal), intent(inout) :: self

        self%enthalpy_cache_valid = .false.

    end subroutine invalidate_enthalpy_cache_thermal

    module subroutine calc_specific_heat_water_thermal(self, state, c_w)
        implicit none
        class(type_thermal), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: c_w
        call self%physics%calc_specific_heat_water(state, c_w)
    end subroutine calc_specific_heat_water_thermal

    module subroutine calc_specific_heat_vapor_thermal(self, state, c_v)
        implicit none
        class(type_thermal), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: c_v
        call self%physics%calc_specific_heat_vapor(state, c_v)
    end subroutine calc_specific_heat_vapor_thermal

end submodule thermal_coefficients
