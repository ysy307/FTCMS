module models_phase_change_fusion
    use, intrinsic :: iso_fortran_env
    use :: iapws, only:type_iapws97, type_iapws06
    use :: module_core, only:type_state
    use :: constitutive_constants, only: &
        Tf0 => water_freezing_point_at_standard_atmospheric_pressure, &
        g => gravity_acceleration, rho_std => reference_water_density, &
        Lf0 => latent_heat_fusion_water_0C, &
        TtoK => celsius_to_kelvin
    use :: physics_constitutive_base, only:abst_constitutive
    use :: models_wrf, only:abst_wrf
    use :: models_phase_change_gcc, only:abst_gcc
    implicit none
    private

    public :: type_fusion

    !>
    !> @brief Model for fusion (melting/freezing) physics.
    !>
    type, extends(abst_constitutive) :: type_fusion
        private
        class(abst_wrf), pointer :: wrf
        class(abst_gcc), pointer :: gcc
    contains
        procedure, pass(self), public :: initialize => initialize_type_fusion
        procedure, pass(self), public :: calc_ice_content
        procedure, pass(self), public :: calc_ice_content_derivatives
        procedure, pass(self), public :: calc_water_content
        procedure, pass(self), public :: calc_water_content_derivatives
        procedure, pass(self), public :: calc_cryo_head_dT
        procedure, pass(self), public :: calc_phase_split
        procedure, pass(self), public :: calc_effective_suction
        procedure, pass(self), public :: calc_freezing_level_set
        procedure, pass(self), public :: project_ice_content
        procedure, pass(self), public :: calc_conserved_target
        procedure, pass(self), public :: solve_local_conserved_equilibrium
        procedure, pass(self), public :: calc_saturation_pressure
        procedure, pass(self), public :: deriv_pressure_ice_water

    end type type_fusion

    !> Suction-domain smoothing scale shared by compute_effective_suction's
    !> smooth maximum and compute_smooth_min's smooth minimum, so the storage
    !> split's max-switch and the pore-volume blend's min-switch close at the
    !> same scale.
    !>
    !> epsilon_s is the suction equivalent of 0.01 K through the Clapeyron
    !> slope rho_w L_f/T_m, i.e. below any temperature difference the
    !> discretization resolves; it exists only so both switches are
    !> differentiable for the Newton linearization.
    real(real64), private, parameter :: SUCTION_SMOOTHING = 1.0d-2 * rho_std * Lf0 / (Tf0 + TtoK)

    !> Let a full pore keep freezing and return the displaced water to the mass
    !> balance, instead of pinning the cryogenic suction at the pore-volume root.
    !>
    !> Pinning the suction is what makes d(theta_i)/dT vanish exactly where the
    !> latent heat is largest: the apparent heat capacity collapses from the
    !> latent-dominated ~1e8 to the sensible ~3e6, and the few nodes it happens
    !> at then dominate the WRMS convergence gate. With the constraint carried
    !> by the stored water instead, theta_l keeps falling with temperature,
    !> theta_i = phi - theta_l keeps growing, and the ice's lower density means
    !> the stored water
    !>   Theta_store = theta_l + alpha*theta_i = alpha*phi + (1-alpha)*theta_l
    !> decreases as freezing proceeds - which the transport equation reads as a
    !> source that expels water from the saturated node.
    logical, parameter, public :: PORE_LIMIT_EXPELS_WATER = .false.

    !> Painter-Karra three-phase equilibrium closure (Vadose Zone J. 13, 2014).
    !>
    !> The ice pressure is not fixed at atmospheric. Liquid-ice equilibrium makes
    !> it the dependent quantity
    !>   p_i = alpha*p_w - rho_i*L_f*ln(T/T_0),   alpha = rho_i/rho_w
    !> so the ice-water pressure difference chi = p_i - p_w responds to BOTH
    !> temperature and the pore pressure - the derivative d(chi)/dp_w = alpha-1
    !> that the p_i = p_a assumption sets to zero, which is what drops the
    !> pressure out of the phase split.
    !>
    !> The condensed phases then partition by volume, not by liquid-equivalent
    !> mass:
    !>   theta_c = theta_SWRC(s_m)                condensed volume (liquid+ice)
    !>   theta_l = theta_SWRC(max(s_m, beta*chi)) most restrictive interface
    !>   theta_i = theta_c - theta_l              volume complement
    !>   theta_a = phi - theta_c                  >= 0 automatically
    !> and the density ratio enters the CONSERVED MASS, theta_l + alpha*theta_i,
    !> not the geometric complement. Putting it in the complement (as
    !> theta_i = (rho_w/rho_i)(Theta - theta_l) does) is what lets the condensed
    !> phases exceed the pore volume and forces a limiter, and the limiter is
    !> what pins theta_l to a temperature-independent expression.
    logical, parameter, public :: PHASE_CLOSURE_PAINTER_KARRA = .false.

    !> Ratio of the air-water to the ice-water surface tension. The same pore
    !> geometry drains at s = 2*sigma_aw/r against air and freezes at
    !> chi = 2*sigma_iw/r against ice, so one retention curve serves both once
    !> the ice-water pressure difference is rescaled by this factor.
    !> sigma_aw = 0.0728 N/m, sigma_iw = 0.033 N/m at 0 C.
    real(real64), private, parameter :: SURFACE_TENSION_RATIO = 0.0728d0 / 0.033d0

    !> Width of the C^1 band over which the pore overfill is removed, in water
    !> content. Small against the porosity so the cap is exact on both sides
    !> outside it, wide enough that the tangent transition is resolvable.
    real(real64), private, parameter :: PORE_EXCESS_SMOOTHING = 1.0d-3

contains

    !> Generalized-Clapeyron freezing suction under homotopy continuation.
    !>
    !> \[ s_f^{\lambda} = s_f^{ref} + \lambda\,(s_f(T,p) - s_f^{ref}),
    !>    \qquad s_f^{ref} = s_f(T^n, p^n) \]
    !>
    !> Every consumer of the raw generalized-Clapeyron suction goes through
    !> here, so the split, the transport potential, the level set and the
    !> quadrature rule cannot end up at different values of \( \lambda \).
    !>
    !> Assumptions: the reference is the previous BDF level, so at that level
    !> the blend is the identity and the known history terms stay independent
    !> of \( \lambda \).
    !> Numerical guarantees: at \( \lambda = 1 \), and whenever the state
    !> carries no continuation parameter or no usable history, the raw model
    !> values are returned unchanged.
    !> Computational complexity: O(1) arithmetic; one state copy when
    !> \( \lambda < 1 \), none otherwise.
    !> Failure behavior: none; a missing reference degrades to \( \lambda = 1 \).
    subroutine calc_freezing_suction(self, state, suction_freezing, dfreezing_dP, dfreezing_dT)
        implicit none
        !> Fusion model
        class(type_fusion), intent(in) :: self
        !> Thermodynamic state, carrying the continuation parameter and history
        type(type_state), intent(in) :: state
        !> Freezing suction [Pa]
        !> Overwritten on exit
        real(real64), intent(inout) :: suction_freezing
        !> d(s_f)/dP [-]
        real(real64), intent(inout), optional :: dfreezing_dP
        !> d(s_f)/dT [Pa/K]
        real(real64), intent(inout), optional :: dfreezing_dT

        type(type_state) :: reference_state
        real(real64), pointer, contiguous, dimension(:) :: temperature_history
        real(real64), pointer, contiguous, dimension(:) :: pressure_history
        real(real64) :: lambda, suction_reference
        logical :: lambda_set, history_set

        call self%gcc%calc(state, suction_freezing)
        if (present(dfreezing_dP)) call self%gcc%deriv_pressure(state, dfreezing_dP)
        if (present(dfreezing_dT)) call self%gcc%deriv_temperature(state, dfreezing_dT)

        if (PHASE_CLOSURE_PAINTER_KARRA) then
            ! Map the p_i = p_a suction onto the ice-water pressure difference
            ! rescaled to an equivalent air-water suction, here rather than in
            ! the split: this is the one entry point every consumer of the
            ! freezing suction goes through, so the composition, the level set
            ! and the flux tangent cannot end up on different definitions.
            block
                real(real64) :: alpha_pk, rho_w_pk, rho_i_pk, pressure_pk
                real(real64) :: ds_dP_pk, ds_dT_pk, s_raw_pk
                real(real64) :: drw_dP, drw_dT, dri_dP, dri_dT
                real(real64) :: dalpha_dT

                call self%calc_rho_water(state, rho_w_pk)
                call self%calc_rho_ice(state, rho_i_pk)
                alpha_pk = 1.0d0

                dalpha_dT = 0.0d0
                if (rho_w_pk > tiny(1.0d0)) then
                    alpha_pk = rho_i_pk / rho_w_pk
                    call self%calc_drho_water_dP(state, drw_dP)
                    call self%calc_drho_water_dT(state, drw_dT)
                    call self%calc_drho_ice_dP(state, dri_dP)
                    call self%calc_drho_ice_dT(state, dri_dT)
                    dalpha_dT = dri_dT / rho_w_pk - rho_i_pk * drw_dT / rho_w_pk**2
                end if
                call state%pressure%get(pressure_pk)

                ds_dP_pk = 0.0d0
                ds_dT_pk = 0.0d0
                if (present(dfreezing_dP)) ds_dP_pk = dfreezing_dP
                if (present(dfreezing_dT)) ds_dT_pk = dfreezing_dT
                s_raw_pk = suction_freezing

                ! chi = (alpha-1)p_w + alpha*s_f. alpha is a state function, so
                ! its own derivative multiplies (p_w + s_f); dropping it in T
                ! costs 1.5e-4 relative against the finite difference.
                !
                ! The pressure derivative deliberately omits the matching
                ! dalpha/dp term: the freezing model declares d(s_f)/dp = 0,
                ! discarding the rho_w(p) dependence it would otherwise carry,
                ! and the two omissions cancel. Measured, keeping only
                ! dalpha/dp breaks that cancellation and the pressure
                ! derivatives stop matching their finite differences.
                suction_freezing = SURFACE_TENSION_RATIO * &
                                   ((alpha_pk - 1.0d0) * pressure_pk + alpha_pk * s_raw_pk)
                if (present(dfreezing_dP)) dfreezing_dP = SURFACE_TENSION_RATIO * &
                    ((alpha_pk - 1.0d0) + alpha_pk * ds_dP_pk)
                if (present(dfreezing_dT)) dfreezing_dT = SURFACE_TENSION_RATIO * &
                    (alpha_pk * ds_dT_pk + dalpha_dT * (pressure_pk + s_raw_pk))
            end block
        end if

        lambda = 1.0d0
        lambda_set = .false.
        call state%continuation_lambda%get(lambda, is_set=lambda_set)
        if (.not. lambda_set) return
        if (lambda >= 1.0d0) return

        nullify (temperature_history)
        nullify (pressure_history)
        history_set = .false.
        call state%temperature_history%get(temperature_history, is_set=history_set)
        if (.not. history_set) return
        if (.not. associated(temperature_history)) return
        if (size(temperature_history) < 2) return
        history_set = .false.
        call state%pressure_history%get(pressure_history, is_set=history_set)

        call reference_state%copy(state)
        call reference_state%temperature%set(temperature_history(2))
        if (history_set) then
            if (associated(pressure_history)) then
                if (size(pressure_history) >= 2) call reference_state%pressure%set(pressure_history(2))
            end if
        end if
        call self%gcc%calc(reference_state, suction_reference)

        suction_freezing = suction_reference + lambda * (suction_freezing - suction_reference)
        if (present(dfreezing_dP)) dfreezing_dP = lambda * dfreezing_dP
        if (present(dfreezing_dT)) dfreezing_dT = lambda * dfreezing_dT

        nullify (temperature_history)
        nullify (pressure_history)
    end subroutine calc_freezing_suction

    !> Painter-Karra three-phase equilibrium split; see PHASE_CLOSURE_PAINTER_KARRA.
    !>
    !> \[ \chi = p_i - p_w = (\alpha-1)p_w + \alpha s_f, \qquad
    !>    \theta_l = \theta_{SWRC}(\max(s_m,\beta\chi)), \qquad
    !>    \theta_i = \theta_c - \theta_l \]
    !>
    !> Assumptions: local liquid-ice equilibrium wherever both are present, a
    !> rigid skeleton, and one pore-size distribution shared by the air-water
    !> and ice-water interfaces.
    !> Numerical guarantees: theta_a = phi - theta_c >= 0 holds by construction
    !> for a retention curve bounded by the porosity, so no pore-volume limiter
    !> is needed and no branch pins theta_l.
    !> Computational complexity: O(1) arithmetic, two retention evaluations.
    !> Failure behavior: none; degenerates to the unfrozen split above T_f.
    subroutine phase_split_painter_karra(self, state, suction_matric, suction_freezing, &
                                         dfreezing_dP, dfreezing_dT, &
                                         theta_condensed, dcondensed_dpressure, &
                                         rho_w, rho_i, density_ratio, &
                                         theta_liquid, theta_ice, &
                                         dliquid_dP, dliquid_dT, dice_dP, dice_dT, &
                                         dtotal_dP, dtotal_dT, &
                                         suction_effective_out, dsuction_eff_dP_out, dsuction_eff_dT_out)
        implicit none
        !> Fusion model
        class(type_fusion), intent(in) :: self
        !> Thermodynamic state
        type(type_state), intent(in) :: state
        !> Matric suction s_m = -p_w [Pa]
        real(real64), intent(in) :: suction_matric
        !> Raw generalized-Clapeyron suction at p_i = p_a [Pa]
        real(real64), intent(in) :: suction_freezing
        !> Tangents of the raw freezing suction
        real(real64), intent(in) :: dfreezing_dP, dfreezing_dT
        !> Condensed-phase volume theta_c and its pressure tangent; on exit the
        !> liquid-equivalent conserved mass replaces theta_c
        real(real64), intent(inout) :: theta_condensed
        real(real64), intent(in) :: dcondensed_dpressure
        !> Phase densities and their ratio rho_w/rho_i
        real(real64), intent(in) :: rho_w, rho_i, density_ratio
        !> Liquid and ice volume fractions [-]
        real(real64), intent(inout) :: theta_liquid, theta_ice
        !> Optional tangents, same meaning as calc_phase_split's
        real(real64), intent(inout), optional :: dliquid_dP, dliquid_dT, dice_dP, dice_dT
        real(real64), intent(inout), optional :: dtotal_dP, dtotal_dT
        real(real64), intent(inout), optional :: suction_effective_out
        real(real64), intent(inout), optional :: dsuction_eff_dP_out, dsuction_eff_dT_out

        real(real64) :: alpha
        real(real64) :: suction_effective, sig_a, sig_b
        real(real64) :: dsuction_eff_dP, dsuction_eff_dT
        real(real64) :: dtheta_dhead_effective
        real(real64) :: dliquid_dpressure, dliquid_dtemperature
        real(real64) :: dice_dpressure, dice_dtemperature
        real(real64) :: drho_w_dP, drho_w_dT, drho_i_dP, drho_i_dT
        real(real64) :: dratio_dP, dratio_dT, dalpha_dP, dalpha_dT

        alpha = 1.0d0 / density_ratio

        ! suction_freezing already carries beta*chi: calc_freezing_suction does
        ! the mapping so every consumer shares it.
        call compute_effective_suction(suction_matric, suction_freezing, suction_effective, sig_a, sig_b)
        dsuction_eff_dP = -sig_a + sig_b * dfreezing_dP
        dsuction_eff_dT = sig_b * dfreezing_dT

        call self%wrf%calc(-suction_effective / (rho_std * g), theta_liquid)
        call self%wrf%deriv(-suction_effective / (rho_std * g), dtheta_dhead_effective)

        dliquid_dpressure = dtheta_dhead_effective * (-dsuction_eff_dP) / (rho_std * g)
        dliquid_dtemperature = dtheta_dhead_effective * (-dsuction_eff_dT) / (rho_std * g)

        ! Volume complement: the ice occupies the space the liquid gave up.
        theta_ice = max(0.0d0, theta_condensed - theta_liquid)
        dice_dpressure = dcondensed_dpressure - dliquid_dpressure
        dice_dtemperature = -dliquid_dtemperature
        if (theta_condensed <= theta_liquid) then
            dice_dpressure = 0.0d0
            dice_dtemperature = 0.0d0
        end if

        dratio_dP = 0.0d0
        dratio_dT = 0.0d0
        if (rho_i > tiny(1.0d0)) then
            call self%calc_drho_water_dP(state, drho_w_dP)
            call self%calc_drho_water_dT(state, drho_w_dT)
            call self%calc_drho_ice_dP(state, drho_i_dP)
            call self%calc_drho_ice_dT(state, drho_i_dT)
            dratio_dP = drho_w_dP / rho_i - rho_w * drho_i_dP / rho_i**2
            dratio_dT = drho_w_dT / rho_i - rho_w * drho_i_dT / rho_i**2
        end if
        dalpha_dP = -dratio_dP / density_ratio**2
        dalpha_dT = -dratio_dT / density_ratio**2

        if (present(dliquid_dP)) dliquid_dP = dliquid_dpressure
        if (present(dliquid_dT)) dliquid_dT = dliquid_dtemperature
        if (present(dice_dP)) dice_dP = dice_dpressure
        if (present(dice_dT)) dice_dT = dice_dtemperature
        ! The density ratio belongs to the conserved mass, not the complement.
        if (present(dtotal_dP)) dtotal_dP = dliquid_dpressure + alpha * dice_dpressure + theta_ice * dalpha_dP
        if (present(dtotal_dT)) dtotal_dT = dliquid_dtemperature + alpha * dice_dtemperature + theta_ice * dalpha_dT
        theta_condensed = theta_liquid + alpha * theta_ice

        if (present(suction_effective_out)) suction_effective_out = suction_effective
        if (present(dsuction_eff_dP_out)) dsuction_eff_dP_out = dsuction_eff_dP
        if (present(dsuction_eff_dT_out)) dsuction_eff_dT_out = dsuction_eff_dT
    end subroutine phase_split_painter_karra

    !> One-sided C^1 ramp: exactly zero up to the corner, exactly the identity
    !> past the band, cubic in between.
    !>
    !> \[ h(d,\varepsilon) = \begin{cases} 0 & d \le 0 \\
    !>    2d^2/\varepsilon - d^3/\varepsilon^2 & 0 < d < \varepsilon \\
    !>    d & d \ge \varepsilon \end{cases} \]
    !>
    !> One-sided because the quantity it caps is already inactive below the
    !> corner: a symmetric band would give the unconstrained state a spurious
    !> tangent. Value and slope match at both breakpoints (h(0)=h'(0)=0,
    !> h(eps)=eps, h'(eps)=1), and h is monotone on the band.
    pure elemental function ramp_band(d, band) result(h)
        implicit none
        !> Argument of the ramp
        real(real64), intent(in) :: d
        !> Width of the regularized band, > 0
        real(real64), intent(in) :: band
        !> Regularized max(0, d)
        real(real64) :: h

        if (d <= 0.0d0) then
            h = 0.0d0
        else if (band <= 0.0d0 .or. d >= band) then
            h = d
        else
            h = 2.0d0 * d**2 / band - d**3 / band**2
        end if
    end function ramp_band

    !> Derivative of ramp_band with respect to its first argument.
    pure elemental function ramp_band_deriv(d, band) result(hp)
        implicit none
        !> Argument of the ramp
        real(real64), intent(in) :: d
        !> Width of the regularized band, > 0
        real(real64), intent(in) :: band
        !> dh/dd
        real(real64) :: hp

        if (d <= 0.0d0) then
            hp = 0.0d0
        else if (band <= 0.0d0 .or. d >= band) then
            hp = 1.0d0
        else
            hp = 4.0d0 * d / band - 3.0d0 * d**2 / band**2
        end if
    end function ramp_band_deriv

    !> C^1 compact-support ramp underlying the smooth max/min below.
    !>
    !> \[ h(d) = \begin{cases} 0 & d \le -\varepsilon_s \\
    !>    (d+\varepsilon_s)^2/(4\varepsilon_s) & -\varepsilon_s < d < \varepsilon_s \\
    !>    d & d \ge \varepsilon_s \end{cases},
    !>    \qquad h'(d) = \begin{cases} 0 & d \le -\varepsilon_s \\
    !>    (d+\varepsilon_s)/(2\varepsilon_s) & -\varepsilon_s < d < \varepsilon_s \\
    !>    1 & d \ge \varepsilon_s \end{cases} \]
    !>
    !> h is the C^1 regularization of the ramp max(0,d): value AND slope match
    !> the corner at both breakpoints, h(-eps_s)=0, h'(-eps_s)=0, h(eps_s)=eps_s,
    !> h'(eps_s)=1. Unlike a hyperbolic smoothing (sqrt(d^2+eps^2)-based), h has
    !> COMPACT support: it equals the exact corner max(0,d) once |d| >= eps_s,
    !> not merely in the limit. That compact support is the entire reason for
    !> this replacement - smooth_max/smooth_min built from h below reduce to
    !> the EXACT max/min outside the band, so no smoothing tail can leak ice
    !> into a state arbitrarily far from the freezing interface (the hyperbolic
    !> form's s_eff > max(s_m,s_f) at every finite state, however far, is what
    !> produced ice at every node of the initial condition).
    pure elemental function smooth_ramp(d) result(h)
        implicit none
        real(real64), intent(in) :: d
        real(real64) :: h

        if (d <= -SUCTION_SMOOTHING) then
            h = 0.0d0
        else if (d >= SUCTION_SMOOTHING) then
            h = d
        else
            h = (d + SUCTION_SMOOTHING)**2 / (4.0d0 * SUCTION_SMOOTHING)
        end if
    end function smooth_ramp

    !> Derivative h'(d) of smooth_ramp; see that function's docstring.
    pure elemental function smooth_ramp_deriv(d) result(hp)
        implicit none
        real(real64), intent(in) :: d
        real(real64) :: hp

        if (d <= -SUCTION_SMOOTHING) then
            hp = 0.0d0
        else if (d >= SUCTION_SMOOTHING) then
            hp = 1.0d0
        else
            hp = (d + SUCTION_SMOOTHING) / (2.0d0 * SUCTION_SMOOTHING)
        end if
    end function smooth_ramp_deriv

    !> Effective suction governing the liquid water content, smoothed max.
    !>
    !> \[ s_{eff} = \max_{\varepsilon}(s_m,s_f) = s_f + h(s_m-s_f) \]
    !>
    !> using the compact-support ramp h above (d(smax)/ds_m = h'(s_m-s_f),
    !> d(smax)/ds_f = 1-h'(s_m-s_f)). Because h(d)=0 for d<=-eps_s and h(d)=d
    !> for d>=eps_s, this is EXACTLY s_f when s_f>=s_m+eps_s and EXACTLY s_m
    !> when s_f<=s_m-eps_s: no tail leaks outside the band |s_m-s_f|<eps_s, in
    !> contrast to the hyperbolic form 0.5*(s_m+s_f+sqrt((s_m-s_f)^2+eps_s^2))
    !> this replaced, which satisfies s_eff > max(s_m,s_f) at every finite
    !> state - so ice (theta_i>0 iff s_eff>s_m) was measured nonzero at every
    !> one of 2874 nodes of a fully unfrozen initial condition.
    !>
    !> The matric suction s_m and the freezing-equivalent suction s_f are two
    !> constraints on the SAME liquid chemical potential, not two contributions
    !> to it: capillarity fixes mu_w through the air-water interface, ice-water
    !> equilibrium fixes it through the generalized Clapeyron relation, and at
    !> equilibrium the two agree. The binding one therefore governs, which is a
    !> maximum, not a sum. Adding them counts the same potential twice and makes
    !> the transport potential differ from the storage potential, so the pair is
    !> not conjugate and the equation has no reachable equilibrium: measured,
    !> that drove the pore to complete saturation within 1500 s against an
    !> experiment that redistributes over 50 h.
    pure subroutine compute_effective_suction(psi_cap, psi_cryo, psi_eff, dpsi_eff_dpsi_cap, dpsi_eff_dpsi_cryo)
        implicit none
        real(real64), intent(in) :: psi_cap, psi_cryo
        real(real64), intent(inout) :: psi_eff
        real(real64), intent(inout), optional :: dpsi_eff_dpsi_cap, dpsi_eff_dpsi_cryo

        real(real64) :: difference, hp

        difference = psi_cap - psi_cryo
        psi_eff = psi_cryo + smooth_ramp(difference)

        if (present(dpsi_eff_dpsi_cap) .or. present(dpsi_eff_dpsi_cryo)) then
            hp = smooth_ramp_deriv(difference)
            if (present(dpsi_eff_dpsi_cap)) dpsi_eff_dpsi_cap = hp
            if (present(dpsi_eff_dpsi_cryo)) dpsi_eff_dpsi_cryo = 1.0d0 - hp
        end if
    end subroutine compute_effective_suction

    !> Smooth minimum of two suctions, used to blend the raw cryogenic suction
    !> with the pore-volume-limit root (calc_limited_cryo_suction) without an
    !> if/else activation switch.
    !>
    !> \[ s_{min} = \min_{\varepsilon}(a,b) = a - h(a-b) \]
    !>
    !> using the SAME compact-support ramp h as compute_effective_suction's
    !> smooth maximum (d(smin)/da = 1-h'(a-b), d(smin)/db = h'(a-b)), so the
    !> storage split's two switches share one width AND one compact band. This
    !> is EXACTLY b when a>=b+eps_s and EXACTLY a when a<=b-eps_s: outside the
    !> band |a-b|<eps_s the cut-off carries no residual tail (see
    !> calc_limited_cryo_suction's BLEND_REACH comment, which relies on this).
    !> As with the hyperbolic form this replaced, smin(a,b) <= min(a,b) always
    !> (h(d) >= max(0,d) everywhere): the blend approaches the pore-volume
    !> limit from BELOW, so the composition never violates the pore constraint
    !> on its own account.
    pure subroutine compute_smooth_min(a, b, smin, dsmin_da, dsmin_db)
        implicit none
        real(real64), intent(in) :: a, b
        real(real64), intent(inout) :: smin, dsmin_da, dsmin_db

        real(real64) :: d, hp

        d = a - b
        smin = a - smooth_ramp(d)
        hp = smooth_ramp_deriv(d)

        dsmin_da = 1.0d0 - hp
        dsmin_db = hp
    end subroutine compute_smooth_min


    !> Cryogenic suction the liquid actually feels, limited by the ice pressure.
    !>
    !> The generalized Clapeyron relation with a non-zero ice pressure is
    !>   p_w/rho_w - p_i/rho_i = L_f ln(T/T_0),
    !> so p_w = -psi_cryo + (rho_w/rho_i) p_i: pressure on the ice raises the
    !> liquid pressure and therefore lowers the suction the liquid experiences.
    !>
    !> While the pore has room the ice is unstressed, p_i = 0, and the full
    !> psi_cryo applies. Once theta_l + theta_i reaches the porosity the ice can
    !> no longer expand, p_i rises, and the suction stops growing. Without this
    !> the suction increases without bound as T falls, so there is no state at
    !> which the inflow stops and water is drawn in until the pore fills - which
    !> at fixed porosity is not a state the soil can occupy. Measured before
    !> this limit: theta_l + theta_i reached the porosity exactly at t = 1500 s,
    !> against an experiment that redistributes over 50 h without saturating.
    !>
    !> The limit is found by bisection on the cryogenic part, which is monotone:
    !> f(psi) = theta_l + theta_i - phi increases with psi and is non-positive at
    !> psi = 0, so a root exists whenever the unlimited value would overfill.
    !>
    !> The optional psi_star_out/needs_blend_out report the SAME root a step
    !> earlier, before the constraint is fully active: whenever the root lies
    !> within BLEND_REACH beyond psi_cryo, compute_smooth_min (used by
    !> compute_blended_effective_suction) needs it to keep the blended suction
    !> C^1 as the state approaches the limit, not only once past it. When the
    !> root is farther than that, needs_blend_out is .false. and no caller
    !> needs a blend - the raw suction is returned as psi_star_out unused.
    subroutine calc_limited_cryo_suction(self, state, psi_cap, psi_cryo, psi_cryo_limited, is_limited, &
                                         psi_star_out, needs_blend_out)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(in) :: psi_cap, psi_cryo
        real(real64), intent(inout) :: psi_cryo_limited
        logical, intent(inout) :: is_limited
        real(real64), intent(inout), optional :: psi_star_out
        logical, intent(inout), optional :: needs_blend_out

        integer(int32), parameter :: MAX_BISECTION = 60
        ! Reach beyond the raw suction over which the limit root is searched
        ! when the constraint is not yet active there. With the compact-
        ! support ramp (smooth_ramp), compute_smooth_min(psi_cryo, psi_star)
        ! reduces to psi_cryo EXACTLY once psi_star - psi_cryo >=
        ! SUCTION_SMOOTHING (see compute_smooth_min's docstring): once
        ! pore_excess is non-positive at psi_cryo + SUCTION_SMOOTHING, the true
        ! (unlimited) root already lies far enough that skipping the blend
        ! introduces NO residual discontinuity at all - unlike the hyperbolic
        ! smooth-min this replaced, whose tail never actually reached the
        ! corner at any finite separation, which is why that version needed a
        ! wide reach (W=64) to keep the neglected tail below the line search's
        ! resolution. The factor of 2 here is headroom against floating-point
        ! round-off in pore_excess and the bisection below, not a smoothing-
        ! tail bound: it no longer exists.
        real(real64), parameter :: BLEND_REACH = 2.0d0 * SUCTION_SMOOTHING
        real(real64) :: porosity, rho_w, rho_i, density_ratio
        real(real64) :: theta_total, low, high, mid
        real(real64) :: excess_raw
        integer(int32) :: k

        psi_cryo_limited = psi_cryo
        is_limited = .false.
        if (present(psi_star_out)) psi_star_out = psi_cryo
        if (present(needs_blend_out)) needs_blend_out = .false.
        if (psi_cryo <= 0.0d0) return

        call state%porosity%get(porosity)
        if (porosity <= 0.0d0) return
        call self%calc_rho_water(state, rho_w)
        call self%calc_rho_ice(state, rho_i)
        if (rho_w <= tiny(1.0d0) .or. rho_i <= tiny(1.0d0)) return
        density_ratio = rho_w / rho_i

        call self%wrf%calc(-psi_cap / (rho_std * g), theta_total)

        excess_raw = pore_excess(psi_cryo)
        if (excess_raw > 0.0d0) then
            is_limited = .true.
            low = 0.0d0
            high = psi_cryo
        else if (pore_excess(psi_cryo + BLEND_REACH) > 0.0d0) then
            ! Not active yet, but the root is within blending reach.
            low = psi_cryo
            high = psi_cryo + BLEND_REACH
        else
            ! The limit is far enough away that no blend is needed.
            return
        end if

        do k = 1, MAX_BISECTION
            mid = 0.5d0 * (low + high)
            if (pore_excess(mid) > 0.0d0) then
                high = mid
            else
                low = mid
            end if
        end do

        if (is_limited) psi_cryo_limited = low
        if (present(psi_star_out)) psi_star_out = low
        if (present(needs_blend_out)) needs_blend_out = .true.

    contains

        !> theta_l + theta_i - phi at a trial cryogenic suction, composed by
        !> the SAME smooth max the split uses (compute_effective_suction), not
        !> an additive superposition of psi_cap and psi_trial. This makes the
        !> root of pore_excess solve exactly the constraint calc_phase_split
        !> enforces,
        !>   E(psi) = r*Theta + (1-r)*theta_l(sigma(psi_cap,psi)) - phi,
        !> which is monotone increasing in psi because theta_l is decreasing
        !> in the effective suction and (1-r) < 0.
        function pore_excess(psi_trial) result(excess)
            implicit none
            real(real64), intent(in) :: psi_trial
            real(real64) :: excess
            real(real64) :: psi_eff_trial, theta_liquid

            call compute_effective_suction(psi_cap, psi_trial, psi_eff_trial)
            call self%wrf%calc(-psi_eff_trial / (rho_std * g), theta_liquid)
            excess = theta_liquid + max(0.0d0, theta_total - theta_liquid) * density_ratio - porosity
        end function pore_excess
    end subroutine calc_limited_cryo_suction

    !> Blend the raw cryogenic suction with the pore-volume-limit root through
    !> compute_smooth_min, then compose the result with the matric suction
    !> through compute_effective_suction's smooth max: a single C^1 effective
    !> suction and its (P,T) tangents, with no is_limited branch anywhere.
    !>
    !> \[ \psi_c = \min_\varepsilon(\psi_{raw},\psi^*), \qquad
    !>    s_{eff} = \max_\varepsilon(s_m,\psi_c) \]
    !>
    !> This is the ONE routine calc_phase_split, calc_cryo_head_dT, and
    !> calc_effective_suction all go through, so their composition cannot
    !> diverge from each other.
    !>
    !> Far from the pore limit (calc_limited_cryo_suction reports
    !> needs_blend=.false.), psi_c = psi_raw exactly and d(psi_c)/dX reduces to
    !> the raw GCC tangents dfreezing_dP/dT - reproducing the previous
    !> unlimited formulas exactly. With the compact-support ramp underlying
    !> compute_smooth_min, this is not an approximation: needs_blend=.false.
    !> is only reported once the true pore-volume-limit root lies at least
    !> SUCTION_SMOOTHING beyond psi_raw (see calc_limited_cryo_suction's
    !> BLEND_REACH comment), and at that separation compute_smooth_min(psi_raw,
    !> psi_star) returns psi_raw EXACTLY, with no neglected tail.
    !>
    !> When blending is needed, psi* moves with pressure through the active
    !> pore-volume constraint
    !>   E(s_m,p,psi*) = r*Theta(s_m) + (1-r)*theta_l(sigma(s_m,psi*)) - phi = 0,
    !> (r = rho_w/rho_i). Differentiating in p at fixed psi* (ds_m/dp = -1) and
    !> applying the implicit function theorem gives dpsi_star_dp = -dE_dp/dE_dpsi
    !> below.
    !>
    !> dpsi_star_dT is taken as zero: E's only T dependence is the density
    !> ratio r(T,p), and the dropped term is ~1e-5 of the Clapeyron slope.
    !> Deep in the limit (psi_raw >> psi*), compute_smooth_min's weights
    !> saturate to dsmin_da -> 0, dsmin_db -> 1, so dpsi_c_dT -> 0 and
    !> dpsi_c_dp -> dpsi_star_dp: the previous limited-branch formulas, and
    !> with them d(theta_l+theta_i)/dp -> 0 (see calc_phase_split's docstring).
    subroutine compute_blended_effective_suction(self, state, suction_matric, suction_freezing, &
                                                 dfreezing_dP, dfreezing_dT, &
                                                 suction_effective, dsuction_eff_dP, dsuction_eff_dT)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(in) :: suction_matric, suction_freezing, dfreezing_dP, dfreezing_dT
        real(real64), intent(inout) :: suction_effective
        real(real64), intent(inout), optional :: dsuction_eff_dP, dsuction_eff_dT

        real(real64) :: psi_cryo_limited_unused, psi_star
        logical :: is_limited_unused, needs_blend
        real(real64) :: psi_c, dsmin_da, dsmin_db, dpsi_c_dp, dpsi_c_dT
        real(real64) :: sig_a, sig_b, sig_a_star, sig_b_star, suction_eff_star
        real(real64) :: rho_w, rho_i, density_ratio
        real(real64) :: dtheta_dhead_matric, theta_prime_matric
        real(real64) :: dtheta_dhead_star, theta_prime_star
        real(real64) :: dE_dpsi, dE_dp, dpsi_star_dp

        ! With the pore constraint carried by the stored water, the suction is
        ! never pinned: theta_l must keep responding to temperature for the
        ! latent capacity to survive. calc_phase_split then caps the phases.
        ! Painter-Karra needs no pore-volume limiter: theta_a = phi - theta_c is
        ! non-negative by construction, so pinning the suction would only
        ! reintroduce the tangent collapse the closure exists to remove.
        if (PORE_LIMIT_EXPELS_WATER .or. PHASE_CLOSURE_PAINTER_KARRA) then
            call compute_effective_suction(suction_matric, suction_freezing, suction_effective, sig_a, sig_b)
            if (present(dsuction_eff_dP)) dsuction_eff_dP = -sig_a + sig_b * dfreezing_dP
            if (present(dsuction_eff_dT)) dsuction_eff_dT = sig_b * dfreezing_dT
            return
        end if

        call calc_limited_cryo_suction(self, state, suction_matric, suction_freezing, &
                                       psi_cryo_limited_unused, is_limited_unused, &
                                       psi_star, needs_blend)

        if (.not. needs_blend) then
            ! Far from the pore limit (or no limit reachable at all): reduces
            ! exactly to the raw GCC suction and its raw tangents.
            psi_c = suction_freezing
            dpsi_c_dp = dfreezing_dP
            dpsi_c_dT = dfreezing_dT
        else
            call compute_smooth_min(suction_freezing, psi_star, psi_c, dsmin_da, dsmin_db)

            ! Implicit-function-theorem tangent of the active pore-volume
            ! root E(s_m,p,psi*)=0 (see docstring above), evaluated at fixed T.
            call self%calc_rho_water(state, rho_w)
            call self%calc_rho_ice(state, rho_i)
            density_ratio = 1.0d0
            if (rho_i > tiny(1.0d0)) density_ratio = rho_w / rho_i

            call self%wrf%deriv(-suction_matric / (rho_std * g), dtheta_dhead_matric)
            theta_prime_matric = -dtheta_dhead_matric / (rho_std * g)

            call compute_effective_suction(suction_matric, psi_star, suction_eff_star, sig_a_star, sig_b_star)
            call self%wrf%deriv(-suction_eff_star / (rho_std * g), dtheta_dhead_star)
            theta_prime_star = -dtheta_dhead_star / (rho_std * g)

            dE_dpsi = (1.0d0 - density_ratio) * theta_prime_star * sig_b_star
            dE_dp = -(density_ratio * theta_prime_matric + &
                      (1.0d0 - density_ratio) * theta_prime_star * sig_a_star)

            if (abs(dE_dpsi) > tiny(1.0d0)) then
                dpsi_star_dp = -dE_dp / dE_dpsi
            else
                ! Degenerate only where theta' itself vanishes (a flat WRF
                ! branch); the active constraint then carries no pressure
                ! sensitivity through this path.
                dpsi_star_dp = 0.0d0
            end if

            dpsi_c_dp = dsmin_da * dfreezing_dP + dsmin_db * dpsi_star_dp
            dpsi_c_dT = dsmin_da * dfreezing_dT
        end if

        call compute_effective_suction(suction_matric, psi_c, suction_effective, sig_a, sig_b)
        if (present(dsuction_eff_dP)) dsuction_eff_dP = -sig_a + sig_b * dpsi_c_dp
        if (present(dsuction_eff_dT)) dsuction_eff_dT = sig_b * dpsi_c_dT
    end subroutine compute_blended_effective_suction

    !>
    !> @brief Initialize fusion model.
    !>
    subroutine initialize_type_fusion(self, wrf, gcc, water, ice)
        implicit none
        class(type_fusion), intent(inout) :: self
        class(abst_wrf), intent(in), target :: wrf
        class(abst_gcc), intent(in), target :: gcc
        type(type_iapws97), intent(in), target :: water
        type(type_iapws06), intent(in), target :: ice

        self%wrf => wrf
        self%gcc => gcc
        self%water => water
        self%ice => ice
    end subroutine initialize_type_fusion


    !> Local phase state: the total water follows the pore pressure, the
    !> effective suction decides how much of it stays liquid.
    !>
    !> \[ \Theta = \theta_{SWRC}(s_m),\qquad
    !>    \theta_w = \theta_{SWRC}(s_{eff}),\qquad
    !>    \theta_i = \frac{\rho_w}{\rho_i}\max(0,\Theta-\theta_w) \]
    !>
    !> with \( s_m = p_a - p_w \) and \( s_{eff} = \max(s_m, s_f) \). Freezing
    !> does not change the conserved water at a node, it only moves it between
    !> the phases, so
    !> \( \Theta = \theta_w + (\rho_i/\rho_w)\theta_i \) holds identically and
    !> \( \partial\Theta/\partial T = 0 \): temperature changes the split, not
    !> the sum. The storage tangent \( \partial\Theta/\partial p_w \) therefore
    !> survives into the frozen zone, where a pressure change moves ice rather
    !> than liquid - the pressure equation never loses its diagonal.
    !>
    !> p_w is never assigned from the Clapeyron relation. It stays the unknown
    !> of the mass balance; the relation enters only through s_f, which sets how
    !> much liquid the temperature permits.
    !>
    !> Pore-volume limit: the raw Clapeyron suction s_f is blended with
    !> calc_limited_cryo_suction's limit root through compute_smooth_min, and
    !> s_eff = sigma(s_m, blend) through the SAME smooth max used elsewhere -
    !> all inside compute_blended_effective_suction, with no is_limited
    !> if/else. Far from the limit this reduces exactly to s_eff=sigma(s_m,s_f)
    !> and the raw tangents. Deep in the limit, on the active constraint
    !> r*Theta+(1-r)*theta_l=phi (r=rho_w/rho_i), the blend's weights saturate
    !> and the implicit function theorem gives the branch's tangents (T does
    !> not appear explicitly in the constraint):
    !>   dtheta_l/dp -> r/(r-1) dTheta/dp,     dtheta_l/dT -> 0.
    !> These feed the same complement formulas used off the constraint, so
    !> dtheta_i/dp -> -r/(r-1) dTheta/dp and dtheta_i/dT -> 0 fall out below
    !> without a separate branch: d(theta_l+theta_i)/dp -> 0, i.e. the pore
    !> approaches staying full while a pressure change only moves water
    !> between the phases that are already there (see
    !> compute_blended_effective_suction's docstring for both limits).
    subroutine calc_phase_split(self, state, total_water, theta_liquid, theta_ice, &
                                dliquid_dP, dliquid_dT, dice_dP, dice_dT, dtotal_dP, dtotal_dT, &
                                suction_effective_out, dsuction_eff_dP_out, dsuction_eff_dT_out)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: total_water, theta_liquid, theta_ice
        real(real64), intent(inout), optional :: dliquid_dP, dliquid_dT, dice_dP, dice_dT
        real(real64), intent(inout), optional :: dtotal_dP
        real(real64), intent(inout), optional :: dtotal_dT
        real(real64), intent(inout), optional :: suction_effective_out
        real(real64), intent(inout), optional :: dsuction_eff_dP_out, dsuction_eff_dT_out

        real(real64) :: pressure, suction_matric, suction_freezing, suction_effective
        real(real64) :: dsuction_eff_dP, dsuction_eff_dT
        real(real64) :: dfreezing_dP, dfreezing_dT
        real(real64) :: dtheta_dhead_matric, dtheta_dhead_effective
        real(real64) :: dtotal_dpressure, dliquid_dpressure, dliquid_dtemperature
        real(real64) :: rho_w, rho_i, density_ratio
        real(real64) :: sig_a, sig_b, porosity_limit
        real(real64) :: theta_ice_unconstrained, pore_excess_value, pore_excess_slope
        real(real64) :: total_water_retention
        logical :: pore_expels

        call state%pressure%get(pressure)
        call self%calc_rho_water(state, rho_w)
        call self%calc_rho_ice(state, rho_i)
        density_ratio = 1.0d0
        if (rho_i > tiny(1.0d0)) density_ratio = rho_w / rho_i

        ! Air pressure is the gauge datum, so the matric suction is -p_w. No
        ! clamp at zero: the retention curve is already flat above saturation
        ! and clamping would put a kink where the physics has none.
        suction_matric = -pressure

        call calc_freezing_suction(self, state, suction_freezing, &
                                   dfreezing_dP=dfreezing_dP, dfreezing_dT=dfreezing_dT)

        call self%wrf%calc(-suction_matric / (rho_std * g), total_water)
        call self%wrf%deriv(-suction_matric / (rho_std * g), dtheta_dhead_matric)
        ! d s_m/d p_w = -1, so d h_m/d p_w = 1/(rho g) and the storage tangent is
        ! positive for a monotone retention curve.
        dtotal_dpressure = dtheta_dhead_matric / (rho_std * g)

        if (PHASE_CLOSURE_PAINTER_KARRA) then
            call phase_split_painter_karra(self, state, suction_matric, suction_freezing, &
                                           dfreezing_dP, dfreezing_dT, &
                                           total_water, dtotal_dpressure, &
                                           rho_w, rho_i, density_ratio, &
                                           theta_liquid, theta_ice, &
                                           dliquid_dP, dliquid_dT, dice_dP, dice_dT, &
                                           dtotal_dP, dtotal_dT, &
                                           suction_effective_out, dsuction_eff_dP_out, dsuction_eff_dT_out)
            return
        end if

        ! Blend the raw cryogenic suction with the pore-volume-limit root
        ! (compute_blended_effective_suction), then compose with the matric
        ! suction through the same smooth max used elsewhere. No is_limited
        ! branch: see that routine's docstring for the far-limit and
        ! deep-limit reductions.
        ! One composition for every consumer (see compute_blended_effective_
        ! suction); with PORE_LIMIT_EXPELS_WATER it returns the unpinned form.
        call compute_blended_effective_suction(self, state, suction_matric, suction_freezing, &
                                               dfreezing_dP, dfreezing_dT, &
                                               suction_effective, dsuction_eff_dP, dsuction_eff_dT)

        call self%wrf%calc(-suction_effective / (rho_std * g), theta_liquid)
        call self%wrf%deriv(-suction_effective / (rho_std * g), dtheta_dhead_effective)

        total_water_retention = total_water
        theta_ice_unconstrained = density_ratio * max(0.0d0, total_water - theta_liquid)
        pore_expels = .false.
        pore_excess_value = 0.0d0
        pore_excess_slope = 0.0d0

        if (PORE_LIMIT_EXPELS_WATER) then
            call state%porosity%get(porosity_limit)
            if (porosity_limit > 0.0d0) pore_expels = .true.
        end if

        if (.not. pore_expels) then
            theta_ice = theta_ice_unconstrained
        else
            ! Remove the pore overfill smoothly rather than switching closures at
            ! the crossing: a hard branch leaves theta_i continuous but its
            ! tangent discontinuous, and the element finite-difference tangent
            ! this residual feeds is documented to need a C^1 residual.
            !
            !   e      = theta_l + theta_i_unconstrained - phi
            !   theta_i = theta_i_unconstrained - h(e)
            !
            ! h is the same compact-support C^1 ramp used for the suction
            ! composition, so below the limit the split is EXACTLY the
            ! unconstrained one and above it theta_i is EXACTLY phi - theta_l.
            pore_excess_value = theta_liquid + theta_ice_unconstrained - porosity_limit
            theta_ice = theta_ice_unconstrained - ramp_band(pore_excess_value, PORE_EXCESS_SMOOTHING)
            pore_excess_slope = ramp_band_deriv(pore_excess_value, PORE_EXCESS_SMOOTHING)
            ! Theta_store = theta_l + alpha*theta_i reproduces the retention
            ! total exactly where the cap is inactive, so one expression covers
            ! both regimes.
            total_water = theta_liquid + theta_ice / density_ratio
        end if

        ! Chain rule through the blended effective suction: dtheta_l/dX =
        ! (dtheta_l/dh)*(-d(s_eff)/dX)/(rho_std g), the same convention used
        ! throughout this module (e.g. calc_water_content_derivatives).
        dliquid_dpressure = dtheta_dhead_effective * (-dsuction_eff_dP) / (rho_std * g)
        dliquid_dtemperature = dtheta_dhead_effective * (-dsuction_eff_dT) / (rho_std * g)

        if (pore_expels) then
            ! One chain rule for both regimes. With
            !   theta_i = theta_i_unc - h(e),  e = theta_l + theta_i_unc - phi
            !   Theta   = theta_l + alpha*theta_i,   alpha = 1/density_ratio
            ! h' = 0 below the limit reproduces the unconstrained split exactly,
            ! including its density-ratio terms, and h' = 1 above it gives
            ! theta_i = phi - theta_l.
            block
                real(real64) :: alpha, dalpha_dP, dalpha_dT
                real(real64) :: drho_w_dP, drho_w_dT, drho_i_dP, drho_i_dT
                real(real64) :: dratio_dP, dratio_dT, phase_gap
                real(real64) :: dice_unc_dP, dice_unc_dT, de_dP, de_dT
                real(real64) :: dice_blend_dP, dice_blend_dT

                alpha = 1.0d0 / density_ratio
                dratio_dP = 0.0d0
                dratio_dT = 0.0d0
                phase_gap = max(0.0d0, total_water_retention - theta_liquid)
                if (rho_i > tiny(1.0d0)) then
                    call self%calc_drho_water_dP(state, drho_w_dP)
                    call self%calc_drho_water_dT(state, drho_w_dT)
                    call self%calc_drho_ice_dP(state, drho_i_dP)
                    call self%calc_drho_ice_dT(state, drho_i_dT)
                    dratio_dP = drho_w_dP / rho_i - rho_w * drho_i_dP / rho_i**2
                    dratio_dT = drho_w_dT / rho_i - rho_w * drho_i_dT / rho_i**2
                end if
                dalpha_dP = -dratio_dP / density_ratio**2
                dalpha_dT = -dratio_dT / density_ratio**2

                dice_unc_dP = density_ratio * (dtotal_dpressure - dliquid_dpressure) + phase_gap * dratio_dP
                dice_unc_dT = -density_ratio * dliquid_dtemperature + phase_gap * dratio_dT
                if (total_water_retention <= theta_liquid) then
                    dice_unc_dP = 0.0d0
                    dice_unc_dT = 0.0d0
                end if

                de_dP = dliquid_dpressure + dice_unc_dP
                de_dT = dliquid_dtemperature + dice_unc_dT
                dice_blend_dP = dice_unc_dP - pore_excess_slope * de_dP
                dice_blend_dT = dice_unc_dT - pore_excess_slope * de_dT

                if (present(dliquid_dP)) dliquid_dP = dliquid_dpressure
                if (present(dliquid_dT)) dliquid_dT = dliquid_dtemperature
                if (present(dice_dP)) dice_dP = dice_blend_dP
                if (present(dice_dT)) dice_dT = dice_blend_dT
                if (present(dtotal_dP)) dtotal_dP = dliquid_dpressure + alpha * dice_blend_dP + theta_ice * dalpha_dP
                if (present(dtotal_dT)) dtotal_dT = dliquid_dtemperature + alpha * dice_blend_dT + theta_ice * dalpha_dT
            end block
            if (present(suction_effective_out)) suction_effective_out = suction_effective
            if (present(dsuction_eff_dP_out)) dsuction_eff_dP_out = dsuction_eff_dP
            if (present(dsuction_eff_dT_out)) dsuction_eff_dT_out = dsuction_eff_dT
            return
        end if

        if (present(dtotal_dP)) dtotal_dP = dtotal_dpressure
        ! Theta = theta_SWRC(s_m(p_w)) has no explicit temperature dependence
        ! in this formulation (s_m = -p_w only); the slot exists so a future
        ! T-dependent total-water relation flows to callers without another
        ! interface change.
        if (present(dtotal_dT)) dtotal_dT = 0.0d0
        if (present(dliquid_dP)) dliquid_dP = dliquid_dpressure
        if (present(dliquid_dT)) dliquid_dT = dliquid_dtemperature
        ! Ice takes the complement, so that dQw + (rho_i/rho_w) dQi reproduces
        ! dTheta exactly and the assembled storage tangent is the retention
        ! one. Deep in the pore-volume limit this is also what makes
        ! d(theta_l+theta_i)/dp approach zero (see compute_blended_effective_
        ! suction's docstring): the smooth composition approaches, but unlike
        ! the old hard switch does not reach exactly, that cancellation.
        !
        ! The ratio r = rho_w/rho_i is itself a state function (IAPWS ice
        ! density varies with T and p), so differentiating theta_i = r*(Theta -
        ! theta_l) exactly gives
        !   d(theta_i) = r*(dTheta - dtheta_l) + (Theta - theta_l)*dr.
        ! Dropping the second term is not a small correction everywhere: deep
        ! in the frozen zone dtheta_l/dT is nearly flat while (Theta - theta_l)
        ! is at its largest, and the neglected term then exceeds the retained
        ! one (measured: the state-function dQi/dT check missed its finite
        ! difference by more than 100 percent at T = -1 C). The identity
        ! theta_l + (rho_i/rho_w) theta_i = Theta is preserved by construction,
        ! so the hydraulic storage tangent dTheta/dp is unaffected; the
        ! correction lands where it belongs, in the enthalpy's dH/dT.
        if (present(dice_dP) .or. present(dice_dT)) then
            block
                real(real64) :: drho_w_dP, drho_w_dT, drho_i_dP, drho_i_dT
                real(real64) :: dratio_dP, dratio_dT, phase_difference

                drho_w_dP = 0.0d0
                drho_w_dT = 0.0d0
                drho_i_dP = 0.0d0
                drho_i_dT = 0.0d0
                dratio_dP = 0.0d0
                dratio_dT = 0.0d0
                phase_difference = max(0.0d0, total_water - theta_liquid)

                if (rho_i > tiny(1.0d0)) then
                    call self%calc_drho_water_dP(state, drho_w_dP)
                    call self%calc_drho_water_dT(state, drho_w_dT)
                    call self%calc_drho_ice_dP(state, drho_i_dP)
                    call self%calc_drho_ice_dT(state, drho_i_dT)
                    dratio_dP = drho_w_dP / rho_i - rho_w * drho_i_dP / rho_i**2
                    dratio_dT = drho_w_dT / rho_i - rho_w * drho_i_dT / rho_i**2
                end if

                if (present(dice_dP)) then
                    dice_dP = density_ratio * (dtotal_dpressure - dliquid_dpressure) + &
                              phase_difference * dratio_dP
                end if
                if (present(dice_dT)) then
                    dice_dT = -density_ratio * dliquid_dtemperature + &
                              phase_difference * dratio_dT
                end if
            end block
        end if
        if (total_water <= theta_liquid) then
            if (present(dice_dP)) dice_dP = 0.0d0
            if (present(dice_dT)) dice_dT = 0.0d0
        end if
        if (present(suction_effective_out)) suction_effective_out = suction_effective
        if (present(dsuction_eff_dP_out)) dsuction_eff_dP_out = dsuction_eff_dP
        if (present(dsuction_eff_dT_out)) dsuction_eff_dT_out = dsuction_eff_dT
    end subroutine calc_phase_split

    !---------------------------------------------------------------------------
    ! Ice Calculations
    !---------------------------------------------------------------------------

    !>
    !> @brief Calculate ice content based on thermodynamic state.
    !>
    !> \[ \theta_{tot} = \theta_{WRF}(-\psi_{cap}) + \frac{\rho_i}{\rho_l}\theta_i \]
    !> \[ \theta_{l,new} = \theta_{WRF}(-\psi_{eff}) \]
    !> \[ \theta_i = \left(\theta_{tot} - \theta_{l,new}\right)\frac{\rho_l}{\rho_i} \]
    !> \(\theta_{WRF}(-\psi_{cap})\) approximates the post-hydraulic-solve liquid content
    !> at current pressure without cryogenic suction, avoiding dependence on the stale
    !> nodal water content field (updated only after convergence).
    !> Assumptions: \(T < 0\) for phase change; \(\theta_i\) from state is the previous-iteration value.
    !> Numerical guarantee: No theoretical error bound available.
    !> Computational complexity: \(O(1)\) arithmetic and memory.
    !> Failure behavior: returns zero ice content when \(T \ge 0\) or \(\theta_{l,new} \ge \theta_{tot}\).
    subroutine calc_ice_content(self, state, ice_content)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: ice_content

        real(real64) :: pressure, temperature
        real(real64) :: psi_cap, psi_cryo, psi_eff
        real(real64) :: theta_l_cap, theta_l_new
        real(real64) :: rho_w, rho_i
        real(real64) :: theta_tot, phi

        call state%temperature%get(temperature)
        call state%pressure%get(pressure)

        if (temperature < 0.0d0) then
            call self%calc_rho_water(state, rho_w)
            call self%calc_rho_ice(state, rho_i)

            if (rho_w <= tiny(1.0d0) .or. rho_i <= tiny(1.0d0)) then
                ice_content = 0.0d0
                return
            end if

            psi_cap = max(0.0d0, -pressure)
            call self%wrf%calc(-psi_cap / (rho_std * g), theta_l_cap)

            call calc_freezing_suction(self, state, psi_cryo)
            call compute_effective_suction(psi_cap, psi_cryo, psi_eff)
            call self%wrf%calc(-psi_eff / (rho_std * g), theta_l_new)

            ! Pore-volume bound on the water-equivalent total theta_tot.
            !
            ! The physical constraint is the instantaneous phase-volume balance
            !   theta_l + theta_i <= phi,
            ! which, with theta_i = (theta_tot - theta_l) * rho_w/rho_i, is
            !   theta_tot <= phi * (rho_i/rho_w) + theta_l * (1 - rho_i/rho_w).
            ! The bound therefore relaxes continuously from phi*(rho_i/rho_w)
            ! (pore entirely ice, theta_l = 0) to phi (pore entirely liquid).
            ! The previous form capped theta_tot at phi*(rho_i/rho_w) uncondition-
            ! ally, i.e. it demanded that the water present could still fit after
            ! freezing ALL of it - far too strict whenever liquid remains. Because
            ! that cap made theta_tot a constant, it also drove the storage
            ! capacity dTheta/dP to zero, so any inflow into a capped node turned
            ! into unbounded pressure growth instead of storage (the +MPa pockets).
            theta_tot = theta_l_cap
            call state%porosity%get(phi)
            if (phi > 0.0d0) then
                theta_tot = min(theta_tot, phi * (rho_i / rho_w) + theta_l_new * (1.0d0 - rho_i / rho_w))
            end if

            if (theta_l_new < theta_tot) then
                ice_content = (theta_tot - theta_l_new) * (rho_w / rho_i)
            else
                ice_content = 0.0d0
            end if
        else
            ice_content = 0.0d0
        end if

    end subroutine calc_ice_content

    !>
    !> @brief Calculate derivatives of ice content w.r.t pressure and temperature.
    !>
    !> \[ \frac{\partial \theta_i}{\partial P} = \frac{\rho_l}{\rho_i}
    !>    \left(\frac{\partial \theta_{l,cap}}{\partial P}
    !>          - \frac{\partial \theta_{l,new}}{\partial P}\right) \]
    !> \[ \frac{\partial \theta_i}{\partial T} = -\frac{\rho_l}{\rho_i}
    !>    \frac{\partial \theta_{l,new}}{\partial T} \]
    !> \(\theta_i\) from state is treated as fixed (previous-iteration value).
    !> Numerical guarantee: No theoretical error bound available.
    !> Computational complexity: \(O(1)\) arithmetic and memory.
    !> Failure behavior: returns zero derivatives when \(T \ge 0\) or \(\theta_{l,new} \ge \theta_{tot}\).
    subroutine calc_ice_content_derivatives(self, state, dice_dP, dice_dT)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: dice_dP
        real(real64), intent(inout) :: dice_dT

        real(real64) :: pressure, temperature
        real(real64) :: psi_cap, psi_cryo, psi_eff
        real(real64) :: d_psi_eff_dpsi_cap, d_psi_eff_dpsi_cryo
        real(real64) :: d_psi_cap_dP, d_psi_cryo_dP, d_psi_cryo_dT
        real(real64) :: d_psi_eff_dP, d_psi_eff_dT
        real(real64) :: theta_l_cap, theta_l_new
        real(real64) :: dtheta_dPin_cap, dtheta_dPin_eff
        real(real64) :: d_theta_cap_dP, d_theta_eff_dP, d_theta_eff_dT
        real(real64) :: rho_w, rho_i
        real(real64) :: theta_tot, phi, d_theta_tot_dP, d_theta_tot_dT, theta_bound
        real(real64) :: drho_w_dP, drho_w_dT, drho_i_dP, drho_i_dT
        real(real64) :: density_ratio, dratio_dP, dratio_dT, phase_difference
        logical :: on_volume_bound

        call state%temperature%get(temperature)
        call state%pressure%get(pressure)

        if (temperature < 0.0d0) then
            call self%calc_rho_water(state, rho_w)
            call self%calc_rho_ice(state, rho_i)

            if (rho_w <= tiny(1.0d0) .or. rho_i <= tiny(1.0d0)) then
                dice_dP = 0.0d0
                dice_dT = 0.0d0
                return
            end if

            if (pressure < 0.0d0) then
                psi_cap = -pressure
                d_psi_cap_dP = -1.0d0
            else
                psi_cap = 0.0d0
                d_psi_cap_dP = 0.0d0
            end if

            call self%wrf%calc(-psi_cap / (rho_std * g), theta_l_cap)
            call self%wrf%deriv(-psi_cap / (rho_std * g), dtheta_dPin_cap)
            d_theta_cap_dP = dtheta_dPin_cap * (-d_psi_cap_dP) / (rho_std * g)

            call calc_freezing_suction(self, state, psi_cryo, &
                                       dfreezing_dP=d_psi_cryo_dP, dfreezing_dT=d_psi_cryo_dT)

            call compute_effective_suction(psi_cap, psi_cryo, psi_eff, d_psi_eff_dpsi_cap, d_psi_eff_dpsi_cryo)
            d_psi_eff_dP = d_psi_eff_dpsi_cap * d_psi_cap_dP + d_psi_eff_dpsi_cryo * d_psi_cryo_dP
            d_psi_eff_dT = d_psi_eff_dpsi_cryo * d_psi_cryo_dT

            call self%wrf%calc(-psi_eff / (rho_std * g), theta_l_new)
            call self%wrf%deriv(-psi_eff / (rho_std * g), dtheta_dPin_eff)

            d_theta_eff_dP = dtheta_dPin_eff * (-d_psi_eff_dP) / (rho_std * g)
            d_theta_eff_dT = dtheta_dPin_eff * (-d_psi_eff_dT) / (rho_std * g)

            ! Pore-volume bound, differentiated consistently with calc_ice_content:
            !   theta_tot <= phi*(rho_i/rho_w) + theta_l*(1 - rho_i/rho_w).
            ! On the bound the total follows the liquid, so
            !   d(theta_tot) = (1 - rho_i/rho_w) * d(theta_l),
            ! and the ice derivatives collapse to dQi = -d(theta_l): further ice
            ! can only replace liquid one-for-one in volume. Off the bound the
            ! total is theta(psi_cap), which has no direct T dependence.
            theta_tot = theta_l_cap
            call state%porosity%get(phi)
            d_theta_tot_dP = d_theta_cap_dP
            d_theta_tot_dT = 0.0d0
            on_volume_bound = .false.
            if (phi > 0.0d0) then
                theta_bound = phi * (rho_i / rho_w) + theta_l_new * (1.0d0 - rho_i / rho_w)
                if (theta_l_cap >= theta_bound) then
                    on_volume_bound = .true.
                    theta_tot = theta_bound
                    d_theta_tot_dP = (1.0d0 - rho_i / rho_w) * d_theta_eff_dP
                    d_theta_tot_dT = (1.0d0 - rho_i / rho_w) * d_theta_eff_dT
                end if
            end if

            if (theta_l_new < theta_tot) then
                density_ratio = rho_w / rho_i
                dice_dP = (d_theta_tot_dP - d_theta_eff_dP) * density_ratio
                dice_dT = (d_theta_tot_dT - d_theta_eff_dT) * density_ratio
                if (.not. on_volume_bound) then
                    call self%calc_drho_water_dP(state, drho_w_dP)
                    call self%calc_drho_water_dT(state, drho_w_dT)
                    call self%calc_drho_ice_dP(state, drho_i_dP)
                    call self%calc_drho_ice_dT(state, drho_i_dT)
                    dratio_dP = drho_w_dP / rho_i - rho_w * drho_i_dP / rho_i**2
                    dratio_dT = drho_w_dT / rho_i - rho_w * drho_i_dT / rho_i**2
                    phase_difference = theta_tot - theta_l_new
                    dice_dP = dice_dP + phase_difference * dratio_dP
                    dice_dT = dice_dT + phase_difference * dratio_dT
                end if
            else
                dice_dP = 0.0d0
                dice_dT = 0.0d0
            end if
        else
            dice_dP = 0.0d0
            dice_dT = 0.0d0
        end if

    end subroutine calc_ice_content_derivatives

    !---------------------------------------------------------------------------
    ! Liquid Water Calculations
    !---------------------------------------------------------------------------

    !>
    !> @brief Calculate liquid water content based on thermodynamic state.
    !>
    subroutine calc_water_content(self, state, water_content)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: water_content

        real(real64) :: pressure
        real(real64) :: psi_cap, psi_cryo, psi_eff
        real(real64) :: dpsi_cryo_dP, dpsi_cryo_dT

        call state%pressure%get(pressure)

        if (pressure < 0.0d0) then
            psi_cap = -pressure
        else
            psi_cap = 0.0d0
        end if

        ! Generalized-suction (freezing = drying analogy): the unfrozen liquid
        ! content follows the retention curve evaluated at the generalized
        ! suction psi_eff = max(psi_cap, psi_cryo). Below freezing the
        ! cryogenic suction psi_cryo(T) lowers the unfrozen content, which is
        ! what makes theta_l temperature-dependent and supplies the apparent
        ! heat capacity to the energy equation.
        ! Through compute_blended_effective_suction, not the limiter directly:
        ! it is the one routine that decides whether the pore-volume root pins
        ! the suction, so calling the limiter here left the saturation test
        ! (calc_saturation_pressure -> here) on the old closure while the split
        ! had moved to the new one.
        call calc_freezing_suction(self, state, psi_cryo, &
                                   dfreezing_dP=dpsi_cryo_dP, dfreezing_dT=dpsi_cryo_dT)
        call compute_blended_effective_suction(self, state, psi_cap, psi_cryo, &
                                               dpsi_cryo_dP, dpsi_cryo_dT, psi_eff)
        call self%wrf%calc(-psi_eff / (rho_std * g), water_content)

    end subroutine calc_water_content

    !>
    !> @brief Calculate derivatives of liquid water content w.r.t pressure and temperature.
    !>
    subroutine calc_water_content_derivatives(self, state, dwater_dP, dwater_dT)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: dwater_dP !> d(theta_l)/dP
        real(real64), intent(inout) :: dwater_dT !> d(theta_l)/dT

        real(real64) :: pressure
        real(real64) :: psi_cap, psi_cryo, psi_eff
        real(real64) :: d_psi_cap_dP, d_psi_cryo_dP, d_psi_cryo_dT
        real(real64) :: d_psi_eff_dpsi_cap, d_psi_eff_dpsi_cryo
        real(real64) :: d_psi_eff_dP, d_psi_eff_dT
        real(real64) :: d_theta_liquid_dPress

        call state%pressure%get(pressure)

        ! Capillary suction
        if (pressure < 0.0d0) then
            psi_cap = -pressure
            d_psi_cap_dP = -1.0d0
        else
            psi_cap = 0.0d0
            d_psi_cap_dP = 0.0d0
        end if

        ! Generalized suction psi_eff = max(psi_cap, psi_cryo(T)). The T
        ! dependence flows through psi_cryo, giving a nonzero dwater_dT: this
        ! is the freezing-curve slope that becomes the apparent heat capacity
        ! in the energy equation's C_TT (thermal_coefficients.F90).
        call calc_freezing_suction(self, state, psi_cryo, &
                                   dfreezing_dP=d_psi_cryo_dP, dfreezing_dT=d_psi_cryo_dT)

        call compute_effective_suction(psi_cap, psi_cryo, psi_eff, d_psi_eff_dpsi_cap, d_psi_eff_dpsi_cryo)
        d_psi_eff_dP = d_psi_eff_dpsi_cap * d_psi_cap_dP + d_psi_eff_dpsi_cryo * d_psi_cryo_dP
        d_psi_eff_dT = d_psi_eff_dpsi_cryo * d_psi_cryo_dT

        call self%wrf%deriv(-psi_eff / (rho_std * g), d_theta_liquid_dPress)

        ! Chain rule: dTheta/dX = (dTheta/dh) * (-d_psi_eff/dX) / (rho_std g).
        dwater_dP = d_theta_liquid_dPress * (-d_psi_eff_dP) / (rho_std * g)
        dwater_dT = d_theta_liquid_dPress * (-d_psi_eff_dT) / (rho_std * g)

    end subroutine calc_water_content_derivatives

    !> Temperature derivative of the generalized-suction head, dh/dT [m/K].
    !>
    !> This derivative belongs to the phase-storage relation. It is not the
    !> thermal liquid conductivity in the hydraulic flux.
    !>
    !> Goes through the same compute_blended_effective_suction as
    !> calc_phase_split and calc_effective_suction, so there is no is_limited
    !> branch here either. Far from the pore-volume limit this reduces to the
    !> raw d(psi_eff)/dT below. Deep in the limit compute_smooth_min's weight
    !> on the raw suction saturates to zero, so d(s_eff)/dT -> 0 smoothly: the
    !> feedback that stops the freezing front from pulling in more water once
    !> the pore is full, now continuous instead of a hard drop. This tangent
    !> feeds K_TT/K_HT (ftcms_assemble.F90's transform_to_total_potential).
    subroutine calc_cryo_head_dT(self, state, dh_dT)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: dh_dT

        real(real64) :: pressure, psi_cap, psi_cryo
        real(real64) :: d_psi_cryo_dP, d_psi_cryo_dT
        real(real64) :: suction_effective, dsuction_eff_dP, dsuction_eff_dT

        call state%pressure%get(pressure)
        psi_cap = max(0.0d0, -pressure)

        call calc_freezing_suction(self, state, psi_cryo, &
                                   dfreezing_dP=d_psi_cryo_dP, dfreezing_dT=d_psi_cryo_dT)

        call compute_blended_effective_suction(self, state, psi_cap, psi_cryo, d_psi_cryo_dP, d_psi_cryo_dT, &
                                               suction_effective, dsuction_eff_dP, dsuction_eff_dT)

        dh_dT = -dsuction_eff_dT / (rho_std * g)
    end subroutine calc_cryo_head_dT

    !>
    !> @brief Generalized suction \(\max(\psi_{cap}, \psi_{cryo})\) [Pa].
    !>
    !> The generalized suction combines capillary and cryogenic constraints for
    !> phase storage. psi_cap is clamped at zero here, unclamped in
    !> calc_phase_split, so the assembly reads the published state value.
    !>
    !> Published through the SAME compute_blended_effective_suction as
    !> calc_phase_split and calc_cryo_head_dT, so the reported suction cannot
    !> diverge from what those routines actually used.
    subroutine calc_effective_suction(self, state, psi_eff)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: psi_eff

        real(real64) :: pressure, psi_cap, psi_cryo, d_psi_cryo_dP, d_psi_cryo_dT

        call state%pressure%get(pressure)
        psi_cap = max(0.0d0, -pressure)
        call calc_freezing_suction(self, state, psi_cryo, &
                                   dfreezing_dP=d_psi_cryo_dP, dfreezing_dT=d_psi_cryo_dT)
        call compute_blended_effective_suction(self, state, psi_cap, psi_cryo, d_psi_cryo_dP, d_psi_cryo_dT, psi_eff)
    end subroutine calc_effective_suction

    !> Freezing-interface level set consistent with calc_phase_split's own
    !> ice-existence switch.
    !>
    !> \[ \phi = (s_f - s_m) + \varepsilon_s, \qquad s_m = -p_w, \qquad
    !>    s_f = \Psi_{ice}^{GCC}(T,p_w) \]
    !>
    !> Ice exists, \(\theta_i>0\), iff \(s_{eff}>s_m\) (calc_phase_split), and
    !> with the compact-support smooth max (compute_effective_suction's
    !> smooth_ramp) that holds exactly iff \(s_f>s_m-\varepsilon_s\), i.e.
    !> \(\phi>0\); \(\phi=0\) is exactly ice onset and \(\phi\le0\) is exactly
    !> ice-free. This is the freezing-interface level set
    !> hydraulic_matrix.F90) should cut elements on, in place of a separately
    !> defined critical temperature: computing it here, from the SAME s_m,s_f
    !> calc_phase_split forms, is what keeps the "cut" elements from drifting
    !> away from where the constitutive split actually places the phase
    !> change.
    !>
    !> s_f is the RAW generalized-Clapeyron suction (self%gcc%calc), the same
    !> suction_freezing calc_phase_split forms before any pore-volume
    !> blending - not the blended/limited suction
    !> compute_blended_effective_suction can return under an active
    !> pore-volume constraint. Away from pore saturation
    !> (calc_limited_cryo_suction's needs_blend=.false., the ordinary case)
    !> the two coincide exactly; only near saturation can this level set lag
    !> the blended split by the (small) blend correction.
    !>
    !> Numerical guarantee: exact given exact s_m, s_f (no further
    !> approximation is introduced).
    !> Computational complexity: O(1) arithmetic and memory.
    !> Failure behavior: none; returns a finite phi for any finite T, p_w.
    subroutine calc_freezing_level_set(self, state, phi)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: phi

        real(real64) :: pressure, suction_matric, suction_freezing

        call state%pressure%get(pressure)
        suction_matric = -pressure
        call calc_freezing_suction(self, state, suction_freezing)
        phi = (suction_freezing - suction_matric) + SUCTION_SMOOTHING
    end subroutine calc_freezing_level_set

    !> Project the outer ice state onto local ice-water equilibrium.
    !>
    !> The liquid content is evaluated from the actual pore pressure, while
    !> the equilibrium unfrozen content is evaluated from the Clapeyron
    !> suction. The proposed phase increment is the liquid-equivalent phase
    !> transfer used by the outer water-conserving solve:
    !> \(\Delta\theta_i=(\theta_l(p_w)-\theta_l^{eq}(T))\rho_w/\rho_i\).
    !> Bounds enforce non-negative ice and fixed pore volume. The returned
    !> increment is zero on the admissible side of either active bound.
    !>
    !> The optional \(dice\_dT\), \(dice\_dP\) are the analytic derivatives of
    !> the returned `projected_ice` with respect to the local temperature and
    !> pore pressure, holding `current_ice` fixed. They chain-rule through the
    !> same `wrf%deriv` and `gcc%deriv_temperature`/`deriv_pressure` primitives
    !> already used elsewhere for the monolithic Jacobian (e.g.
    !> `calc_water_content_derivatives`), so no new closed-form relation is
    !> introduced. On an active bound the derivative is that of the bound
    !> itself (zero at the ice-free bound; \(-d(liquid\_equilibrium)/dX\) at
    !> the full-ice bound, since the bound value \(\phi-\theta_l^{eq}(T)\)
    !> still depends on \(T,P\) through the Clapeyron suction).
    subroutine project_ice_content(self, state, projected_ice, ice_increment, equilibrium_error, active_bound, &
                                    dice_dT, dice_dP)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: projected_ice
        real(real64), intent(inout) :: ice_increment
        real(real64), intent(inout) :: equilibrium_error
        integer(int32), intent(inout), optional :: active_bound
        real(real64), intent(inout), optional :: dice_dT
        real(real64), intent(inout), optional :: dice_dP

        real(real64) :: current_ice, liquid_pressure, liquid_equilibrium
        real(real64) :: pressure, psi_cap, psi_cryo, rho_w, rho_i, porosity, upper_bound, unconstrained_ice
        integer(int32) :: bound_state
        real(real64) :: dliquid_dP, dliquid_dT, dtheta_eq_dh, dpsi_cryo_dT, dpsi_cryo_dP
        real(real64) :: dliquid_eq_dT, dliquid_eq_dP
        real(real64) :: drho_w_dT, drho_w_dP, drho_i_dT, drho_i_dP, drho_ratio_dT, drho_ratio_dP
        ! Active-set tolerance for the ice-free/full-ice complementarity
        ! bounds. It matches the outer phase-content discretization tolerance
        ! so a node is not classified simultaneously as bound-converged and
        ! as an interior point requiring pressure equality.
        real(real64), parameter :: BOUND_TOLERANCE = 1.0d-3
        logical :: ice_is_set, at_freezing_onset
        real(real64) :: onset_content_tolerance

        current_ice = 0.0d0
        ice_is_set = .false.
        call state%ice_content%get(current_ice, ice_is_set)
        if (.not. ice_is_set) current_ice = 0.0d0
        call state%porosity%get(porosity)
        call self%calc_water_content(state, liquid_pressure)
        call calc_freezing_suction(self, state, psi_cryo)
        call self%wrf%calc(-max(0.0d0, psi_cryo) / (rho_std * g), liquid_equilibrium)
        call self%calc_rho_water(state, rho_w)
        call self%calc_rho_ice(state, rho_i)

        if (rho_w <= tiny(1.0d0) .or. rho_i <= tiny(1.0d0)) then
            projected_ice = max(0.0d0, current_ice)
            ice_increment = 0.0d0
            equilibrium_error = 0.0d0
            if (present(active_bound)) active_bound = 0
            if (present(dice_dT)) dice_dT = 0.0d0
            if (present(dice_dP)) dice_dP = 0.0d0
            return
        end if

        unconstrained_ice = current_ice + (liquid_pressure - liquid_equilibrium) * rho_w / rho_i
        upper_bound = max(0.0d0, porosity - liquid_equilibrium)
        projected_ice = min(max(unconstrained_ice, 0.0d0), upper_bound)
        ice_increment = projected_ice - current_ice

        onset_content_tolerance = 64.0d0 * epsilon(1.0d0) * &
                                  max(1.0d0, abs(liquid_pressure), abs(liquid_equilibrium))
        at_freezing_onset = current_ice <= BOUND_TOLERANCE .and. &
                            abs(liquid_pressure - liquid_equilibrium) <= onset_content_tolerance
        bound_state = 0
        if (unconstrained_ice <= 0.0d0) then
            bound_state = -1
        else if (unconstrained_ice >= upper_bound) then
            bound_state = 1
        end if
        if (present(active_bound)) active_bound = bound_state

        call state%pressure%get(pressure)
        psi_cap = max(0.0d0, -pressure)
        equilibrium_error = 0.0d0
        if (abs(ice_increment) > BOUND_TOLERANCE .or. &
            (current_ice > BOUND_TOLERANCE .and. current_ice < upper_bound - BOUND_TOLERANCE)) then
            equilibrium_error = abs(max(0.0d0, psi_cryo) - psi_cap)
        end if

        if (present(dice_dT) .or. present(dice_dP)) then
            ! d(liquid_equilibrium)/dX, X in {T,P}: liquid_equilibrium =
            ! wrf(-max(0,psi_cryo)/(rho_std g)). When psi_cryo<=0 the argument
            ! is pinned at h=0 (saturated branch), independent of T,P.
            dliquid_eq_dT = 0.0d0
            dliquid_eq_dP = 0.0d0
            if (psi_cryo > 0.0d0) then
                call self%wrf%deriv(-psi_cryo / (rho_std * g), dtheta_eq_dh)
                ! Re-entering the helper reproduces psi_cryo and returns the
                ! tangents belonging to that same continuation value.
                call calc_freezing_suction(self, state, psi_cryo, &
                                           dfreezing_dP=dpsi_cryo_dP, dfreezing_dT=dpsi_cryo_dT)
                dliquid_eq_dT = -dtheta_eq_dh * dpsi_cryo_dT / (rho_std * g)
                dliquid_eq_dP = -dtheta_eq_dh * dpsi_cryo_dP / (rho_std * g)
            end if

            ! d(liquid_pressure)/dX: liquid_pressure = theta_l(p_w), reuses
            ! the same relation as calc_water_content_derivatives (dT is
            ! identically zero there since psi_cap depends only on p_w).
            call self%calc_water_content_derivatives(state, dliquid_dP, dliquid_dT)

            select case (bound_state)
            case (0)
                ! unconstrained_ice also depends on T,P through rho_w/rho_i:
                ! rho_w is pinned below the freezing point (IAPWS-97 branch
                ! in calc_rho_water) but rho_i is not (IAPWS-06 has no such
                ! floor), so d(rho_w/rho_i)/dX is nonzero and must be
                ! included by the product rule, not just the dtheta term.
                call self%calc_drho_water_dT(state, drho_w_dT)
                call self%calc_drho_water_dP(state, drho_w_dP)
                call self%calc_drho_ice_dT(state, drho_i_dT)
                call self%calc_drho_ice_dP(state, drho_i_dP)
                drho_ratio_dT = drho_w_dT / rho_i - rho_w * drho_i_dT / rho_i**2
                drho_ratio_dP = drho_w_dP / rho_i - rho_w * drho_i_dP / rho_i**2
                if (present(dice_dT)) dice_dT = (dliquid_dT - dliquid_eq_dT) * rho_w / rho_i + &
                    (liquid_pressure - liquid_equilibrium) * drho_ratio_dT
                if (present(dice_dP)) dice_dP = (dliquid_dP - dliquid_eq_dP) * rho_w / rho_i + &
                    (liquid_pressure - liquid_equilibrium) * drho_ratio_dP
            case (-1)
                if (at_freezing_onset) then
                    ! Hansson's crossing reset places the iterate exactly at
                    ! the critical freezing temperature. Use the frozen-side
                    ! (maximum apparent-capacity) tangent there, although the
                    ! admissible ice value is still the lower bound Qi=0.
                    call self%calc_drho_water_dT(state, drho_w_dT)
                    call self%calc_drho_water_dP(state, drho_w_dP)
                    call self%calc_drho_ice_dT(state, drho_i_dT)
                    call self%calc_drho_ice_dP(state, drho_i_dP)
                    drho_ratio_dT = drho_w_dT / rho_i - rho_w * drho_i_dT / rho_i**2
                    drho_ratio_dP = drho_w_dP / rho_i - rho_w * drho_i_dP / rho_i**2
                    if (present(dice_dT)) dice_dT = (dliquid_dT - dliquid_eq_dT) * rho_w / rho_i + &
                        (liquid_pressure - liquid_equilibrium) * drho_ratio_dT
                    if (present(dice_dP)) dice_dP = (dliquid_dP - dliquid_eq_dP) * rho_w / rho_i + &
                        (liquid_pressure - liquid_equilibrium) * drho_ratio_dP
                else
                    if (present(dice_dT)) dice_dT = 0.0d0
                    if (present(dice_dP)) dice_dP = 0.0d0
                end if
            case (1)
                if (present(dice_dT)) dice_dT = -dliquid_eq_dT
                if (present(dice_dP)) dice_dP = -dliquid_eq_dP
            end select
        end if
    end subroutine project_ice_content

    !> Water+ice-equivalent content conserved locally when a node's hydraulic
    !> conductivity has collapsed under freezing impedance (no flux in/out
    !> over the time step): \(\theta_w^n+(\rho_i/\rho_w)\theta_i^n\), read from
    !> the previous accepted step's pressure/ice history (index 2; index 1 is
    !> the current, still-being-solved step). Returns `.false.` if history is
    !> not yet available (e.g. the first step) or the state is degenerate.
    subroutine calc_conserved_target(self, state, target_total_water, available)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: target_total_water
        logical, intent(inout) :: available

        type(type_state) :: local_state
        real(real64), pointer, contiguous, dimension(:) :: pressure_history
        real(real64), pointer, contiguous, dimension(:) :: ice_content_history
        real(real64) :: liquid_content_prev, rho_w, rho_i

        available = .false.
        nullify (pressure_history, ice_content_history)
        call state%get(pressure_history=pressure_history)
        call state%ice_content_history%get(ice_content_history)
        if (.not. associated(pressure_history)) return
        if (.not. associated(ice_content_history)) return
        if (size(pressure_history) < 2 .or. size(ice_content_history) < 2) return

        call self%calc_rho_water(state, rho_w)
        call self%calc_rho_ice(state, rho_i)
        if (rho_w <= tiny(1.0d0) .or. rho_i <= tiny(1.0d0)) return

        call local_state%copy(state)
        call local_state%pressure%set(pressure_history(2))
        call self%calc_water_content(local_state, liquid_content_prev)

        target_total_water = liquid_content_prev + (rho_i / rho_w) * ice_content_history(2)
        available = .true.
    end subroutine calc_conserved_target

    !> Solve the local, flux-free phase equilibrium at fixed temperature.
    !>
    !> When hydraulic conductivity has collapsed under freezing impedance
    !> (see plan spicy-sauteeing-scroll.md, WP3 refinement), a node no longer
    !> exchanges water with its neighbors within a time step. Its converged
    !> outer-loop state is then the solution of two LOCAL equations, not a
    !> spatially coupled one: total water+ice mass is conserved (no flux) and
    !> the pore pressure satisfies the same Clapeyron liquid equilibrium the
    !> ice update is driving toward, \(\theta_l(p_w)=\theta_l^{eq}(T,p_w)\).
    !> This reduces to a SCALAR Newton iteration on \(p_w\) alone (T is fixed,
    !> supplied by the already-converged monolithic solve for this outer
    !> iteration); ice content is then read off directly from the mass
    !> constraint, not iterated. Reuses the same `wrf%deriv`,
    !> `gcc%deriv_pressure` primitives as `project_ice_content`'s `dice_dP`.
    !>
    !> `target_total_water` is the node's water+ice-equivalent content at the
    !> start of the time step (before this step's flux, assumed negligible in
    !> this regime): \(\theta_w^{n}+(\rho_i/\rho_w)\theta_i^{n}\).
    !>
    !> Returns `converged=.false.` (caller should fall back to the existing
    !> outer Picard/Anderson update) if Newton does not reach the pressure
    !> tolerance within the iteration budget, or the state is degenerate.
    subroutine solve_local_conserved_equilibrium(self, state, target_total_water, new_pressure, new_ice, converged)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(in) :: target_total_water
        real(real64), intent(inout) :: new_pressure
        real(real64), intent(inout) :: new_ice
        logical, intent(inout) :: converged

        type(type_state) :: local_state
        real(real64) :: pressure, psi_cryo, liquid_equilibrium, liquid_content
        real(real64) :: dliquid_dP, dliquid_dT_unused, dtheta_eq_dh, dpsi_cryo_dP, dliquid_eq_dP
        real(real64) :: residual, jacobian, delta_pressure, rho_w, rho_i, porosity, upper_bound
        integer(int32) :: iter
        integer(int32), parameter :: MAX_ITER = 30
        real(real64), parameter :: CONTENT_TOL = 1.0d-8
        real(real64), parameter :: MAX_STEP = 5.0d4 ! Pa; safeguards against runaway Newton steps

        converged = .false.
        call local_state%copy(state)
        call state%pressure%get(pressure)
        call self%calc_rho_water(state, rho_w)
        call self%calc_rho_ice(state, rho_i)
        call state%porosity%get(porosity)
        if (rho_w <= tiny(1.0d0) .or. rho_i <= tiny(1.0d0)) return

        do iter = 1, MAX_ITER
            call local_state%pressure%set(pressure)
            call self%calc_water_content(local_state, liquid_content)
            call calc_freezing_suction(self, local_state, psi_cryo)
            call self%wrf%calc(-max(0.0d0, psi_cryo) / (rho_std * g), liquid_equilibrium)

            residual = liquid_content - liquid_equilibrium
            if (abs(residual) <= CONTENT_TOL) then
                converged = .true.
                exit
            end if

            call self%calc_water_content_derivatives(local_state, dliquid_dP, dliquid_dT_unused)
            dliquid_eq_dP = 0.0d0
            if (psi_cryo > 0.0d0) then
                call self%wrf%deriv(-psi_cryo / (rho_std * g), dtheta_eq_dh)
                call calc_freezing_suction(self, local_state, psi_cryo, dfreezing_dP=dpsi_cryo_dP)
                dliquid_eq_dP = -dtheta_eq_dh * dpsi_cryo_dP / (rho_std * g)
            end if
            jacobian = dliquid_dP - dliquid_eq_dP
            if (abs(jacobian) <= tiny(1.0d0)) return

            delta_pressure = -residual / jacobian
            delta_pressure = sign(min(abs(delta_pressure), MAX_STEP), delta_pressure)
            pressure = pressure + delta_pressure
        end do
        if (.not. converged) return

        new_pressure = pressure
        upper_bound = max(0.0d0, porosity - liquid_content)
        new_ice = min(max((target_total_water - liquid_content) * rho_w / rho_i, 0.0d0), upper_bound)
    end subroutine solve_local_conserved_equilibrium

    !> Calculate the pore-water pressure at which the gas-filled pore volume vanishes.
    !>
    !> \[ \theta_l(T,p_b)+\theta_i(T,p_b)=\phi \]
    !> Assumptions: local phase equilibrium, monotone WRF, and fixed porosity.
    !> Numerical guarantee: returns a bracketed bisection solution when a transition exists.
    !> Computational complexity: \(O(N_b)\) arithmetic and \(O(1)\) memory.
    !> Failure behavior: reports an unsaturated state when no finite transition is bracketed.
    subroutine calc_saturation_pressure(self, state, saturation_pressure, is_saturated)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: saturation_pressure
        logical, intent(inout) :: is_saturated

        integer(int32), parameter :: max_bisection_iterations = 64
        real(real64), parameter :: pressure_floor = -1.0d12
        type(type_state) :: local_state
        real(real64) :: pressure, porosity, ice_content
        real(real64) :: pressure_low, pressure_high, pressure_mid
        real(real64) :: volume_current, volume_low, volume_high, volume_mid
        real(real64) :: volume_tolerance
        integer(int32) :: i

        saturation_pressure = 0.0d0
        is_saturated = .false.
        call state%pressure%get(pressure)
        call state%porosity%get(porosity)
        call state%ice_content%get(ice_content)
        if (porosity <= 0.0d0) return

        call local_state%copy(state)
        call evaluate_unconstrained_phase_volume(local_state, pressure, volume_current)
        volume_tolerance = 64.0d0 * epsilon(1.0d0) * max(1.0d0, porosity)
        if (volume_current < porosity - volume_tolerance) return

        pressure_high = 0.0d0
        call evaluate_unconstrained_phase_volume(local_state, pressure_high, volume_high)
        if (volume_high < porosity - volume_tolerance) return

        pressure_low = -rho_std * g
        call evaluate_unconstrained_phase_volume(local_state, pressure_low, volume_low)
        do while (volume_low >= porosity .and. pressure_low > pressure_floor)
            pressure_low = max(10.0d0 * pressure_low, pressure_floor)
            call evaluate_unconstrained_phase_volume(local_state, pressure_low, volume_low)
        end do
        if (volume_low >= porosity) return

        do i = 1, max_bisection_iterations
            pressure_mid = 0.5d0 * (pressure_low + pressure_high)
            call evaluate_unconstrained_phase_volume(local_state, pressure_mid, volume_mid)
            if (volume_mid >= porosity) then
                pressure_high = pressure_mid
            else
                pressure_low = pressure_mid
            end if
        end do

        saturation_pressure = 0.5d0 * (pressure_low + pressure_high)
        is_saturated = pressure >= saturation_pressure - &
                       max(1.0d-8 * max(abs(saturation_pressure), 1.0d0), epsilon(1.0d0))

    contains
        subroutine evaluate_unconstrained_phase_volume(candidate_state, candidate_pressure, phase_volume)
            implicit none
            type(type_state), intent(inout) :: candidate_state
            real(real64), intent(in) :: candidate_pressure
            real(real64), intent(inout) :: phase_volume

            real(real64) :: theta_liquid

            call candidate_state%pressure%set(candidate_pressure)
            call self%calc_water_content(candidate_state, theta_liquid)
            phase_volume = theta_liquid + ice_content
        end subroutine evaluate_unconstrained_phase_volume
    end subroutine calc_saturation_pressure

    !>
    !> @brief Calculate derivative of ice pressure w.r.t water pressure.
    !>
    subroutine deriv_pressure_ice_water(self, state, deriv)
        implicit none
        class(type_fusion), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: deriv

        call self%gcc%deriv_pressure_ice_water(state, deriv)

    end subroutine deriv_pressure_ice_water

end module models_phase_change_fusion
