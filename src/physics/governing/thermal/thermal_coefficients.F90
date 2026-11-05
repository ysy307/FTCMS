!> Implementation overview
!>
!> Algorithm:
!> - Calculates thermal coefficients and governing equation terms
!> - Assembles enthalpy, transient, mass, diffusion, advection, and latent terms
submodule(physics_governing_thermal) thermal_coefficients
    use :: constitutive_constants, only:tangent_pressure, tangent_temperature
    use :: constitutive_laws, only:evaluate_enthalpy_density, evaluate_enthalpy_density_tangent
    implicit none

contains

    !> Calculate enthalpy (internal energy) density per unit volume
    !>
    !> Mathematical definition:
    !> - Equation of state: \( U = f(T, \phi, \dots) \)
    module subroutine calc_enthalpy_density_thermal(self, material_id, state, U, tangent)
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
        !> \((\partial U/\partial p_w,\partial U/\partial T)\) [J/(m3 Pa), J/(m3 K)]
        !> Overwritten on exit when present
        real(real64), intent(inout), optional :: tangent(2)

        ! Local variables
        real(real64) :: temperature
        real(real64) :: porosity, Qw, Qi, Qv
        real(real64) :: rho_s, rho_w, rho_i
        real(real64) :: c_s, c_w, c_i, c_v
        real(real64) :: Lf, Lv, solid_heat_capacity
        real(real64) :: dT(2), dQw(2), dQi(2), dQv(2)
        real(real64) :: drho_w(2), dc_w(2), drho_i(2), dc_i(2), dc_v(2), dLv(2)
        logical :: with_tangent, has_water, has_ice, has_vapor

        call state%temperature%get(temperature)
        call state%porosity%get(porosity)
        call state%water_content%get(Qw)
        call state%ice_content%get(Qi)
        call state%vapor_content%get(Qv)

        with_tangent = present(tangent)
        dT = [0.0d0, 0.0d0]
        dT(tangent_temperature) = 1.0d0
        dQw(:) = 0.0d0
        dQi(:) = 0.0d0
        dQv(:) = 0.0d0
        if (with_tangent) then
            call state%dQw_dP%get(dQw(tangent_pressure))
            call state%dQw_dT%get(dQw(tangent_temperature))
            call state%dQi_dP%get(dQi(tangent_pressure))
            call state%dQi_dT%get(dQi(tangent_temperature))
            call state%dQv_dP%get(dQv(tangent_pressure))
            call state%dQv_dT%get(dQv(tangent_temperature))
        end if

        ! A phase contributes once it is present, or once its content moves:
        ! the term's partial in that content is needed even at zero content.
        has_water = Qw > 0.0d0 .or. any(dQw /= 0.0d0)
        has_ice = Qi > 0.0d0 .or. any(dQi /= 0.0d0)
        has_vapor = Qv > 0.0d0 .or. any(dQv /= 0.0d0)

        ! Properties of an absent phase multiply a zero content and stay zero.
        solid_heat_capacity = 0.0d0
        rho_w = 0.0d0
        c_w = 0.0d0
        rho_i = 0.0d0
        c_i = 0.0d0
        c_v = 0.0d0
        Lf = 0.0d0
        Lv = 0.0d0
        drho_w(:) = 0.0d0
        dc_w(:) = 0.0d0
        drho_i(:) = 0.0d0
        dc_i(:) = 0.0d0
        dc_v(:) = 0.0d0
        dLv(:) = 0.0d0

        ! Solid phase
        if (porosity > 0.0d0) then
            call self%physics%get_density_solid(material_id, rho_s)
            call self%physics%get_specific_heat_solid(material_id, c_s)
            solid_heat_capacity = rho_s * c_s * (1.0d0 - porosity)
        end if

        ! Water phase; vapour is liquid-equivalent, so it carries the liquid density.
        if (has_water .or. has_vapor) then
            if (with_tangent) then
                call self%physics%calc_density_water(state, rho_w, drho_w)
            else
                call self%physics%calc_density_water(state, rho_w)
            end if
        end if
        if (has_water) then
            if (with_tangent) then
                call self%physics%calc_specific_heat_water(state, c_w, dc_w)
            else
                call self%physics%calc_specific_heat_water(state, c_w)
            end if
        end if

        ! Ice phase, with the latent heat of fusion. Additional segregation-sink
        ! latent heat is handled by an explicit source term. The fusion model's
        ! latent heat is a constant, so it carries no tangent.
        if (has_ice) then
            if (with_tangent) then
                call self%physics%calc_density_ice(state, rho_i, drho_i)
                call self%physics%calc_specific_heat_ice(state, c_i, dc_i)
            else
                call self%physics%calc_density_ice(state, rho_i)
                call self%physics%calc_specific_heat_ice(state, c_i)
            end if
            call self%physics%calc_latent_heat_fusion(material_id, state, Lf)
        end if

        ! Vapor phase, with the latent heat of vaporization.
        if (has_vapor) then
            if (with_tangent) then
                call self%physics%calc_specific_heat_vapor(state, c_v, dc_v)
                call self%physics%calc_latent_heat_vaporization(material_id, state, Lv, dLv)
            else
                call self%physics%calc_specific_heat_vapor(state, c_v)
                call self%physics%calc_latent_heat_vaporization(material_id, state, Lv)
            end if
        end if

        if (with_tangent) then
            call evaluate_enthalpy_density_tangent(solid_heat_capacity, Lf, temperature, dT, Qw, dQw, Qi, dQi, Qv, dQv, &
                                                   rho_w, drho_w, c_w, dc_w, rho_i, drho_i, c_i, dc_i, c_v, dc_v, &
                                                   Lv, dLv, U, tangent)
        else
            call evaluate_enthalpy_density(solid_heat_capacity, Lf, temperature, Qw, Qi, Qv, rho_w, c_w, rho_i, c_i, &
                                           c_v, Lv, U)
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

    ! ==========================================================================
    ! Transient / Mass Terms (Refactored for NR & Picard)
    ! ==========================================================================

    ! ==========================================================================
    ! Flux Terms (Diffusion / Advection / Latent Source)
    ! ==========================================================================

    !> Calculate diffusion term matrix
    !>
    !> Mathematical definition:
    !> - Assembles thermal conductivity tensor based on spatial dimensions
    module subroutine compute_diffusion_term_thermal(self, material_id, state, D_TT, isotropic_tangent)
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
        !> \((\partial\lambda_0/\partial p_w,\partial\lambda_0/\partial T)\) of the isotropic conductivity
        !> the tensor is built from; overwritten on exit when present
        real(real64), intent(inout), optional :: isotropic_tangent(2)

        type(type_state_thc) :: lambda
        real(real64) :: lambda_0

        call self%physics%calc_thermal_conductivity(material_id, state, lambda)
        ! Dispersion is driven by the flux, which no (p_w,T) tangent moves.
        if (present(isotropic_tangent)) then
            call self%physics%calc_thermal_conductivity(material_id, state, lambda_0, isotropic_tangent)
        end if

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

    ! ==========================================================================
    ! Enthalpy Cache Management
    ! ==========================================================================

    module subroutine calc_specific_heat_water_thermal(self, state, c_w)
        implicit none
        class(type_thermal), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: c_w
        call self%physics%calc_specific_heat_water(state, c_w)
    end subroutine calc_specific_heat_water_thermal

end submodule thermal_coefficients
