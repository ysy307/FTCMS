submodule(models_hcf) hcf_vg
    use :: constitutive_laws, only:evaluate_kr_van_genuchten, evaluate_kr_van_genuchten_tangent, &
        evaluate_van_genuchten_mualem_argument, evaluate_kr_van_genuchten_mualem, &
        evaluate_kr_van_genuchten_mualem_tangent
    implicit none
contains

    !> Cache the fixed incomplete-beta normalization for this material.
    module subroutine initialize_hcf_base_vg(self)
        implicit none
        class(type_hcf_base_vg), intent(inout) :: self

        real(real64) :: p, q

        associate (params => self%parent%config)
            if (params%n1 > 1.0d0) then
                p = params%m1 + 1.0d0 / params%n1
                q = 1.0d0 - 1.0d0 / params%n1
            else
                p = 0.0d0
                q = 0.0d0
            end if
            call self%incomplete_beta%initialize(p, q)
        end associate
    end subroutine initialize_hcf_base_vg

    !> @brief Mualem relative permeability for the van Genuchten retention model.
    !>
    !> Mathematical definition (van Genuchten et al., 1991; Hansson et al., 2004,
    !> Eq. [3]) for INDEPENDENT \(m\) and \(n\):
    !> \[ k_r = S_e^{l}\,\bigl[I_{\zeta}(m + 1/n,\; 1 - 1/n)\bigr]^2, \qquad
    !>    \zeta = S_e^{1/m} \]
    !> where \(I_x(a,b)\) is the regularized incomplete beta function.  Only when
    !> \(m = 1 - 1/n\) does this reduce to the closed form
    !> \(k_r = S_e^{l}[1-(1-S_e^{1/m})^m]^2\); using the closed form with an
    !> independent \(m\) (e.g. the Kanagawa sandy loam fit m=0.2, n=1.48)
    !> underestimates \(k_r\) by up to an order of magnitude in the wet-to-medium
    !> suction range and correspondingly suppresses cryosuction redistribution.
    !>
    !> Assumptions: \(n > 1\) (otherwise the Mualem integral is ill-defined and
    !> the closed form is used as a fallback).  Computational complexity: O(1);
    !> the continued fraction converges in a few dozen terms.  Failure behavior:
    !> returns the closed-form value if the beta continued fraction fails to
    !> converge (does not abort). The head derivative is that of the branch the
    !> value came from.
    module subroutine calc_kr_base_vg(self, h, kr, dkr_dh)
        implicit none
        class(type_hcf_base_vg), intent(in) :: self
        real(real64), intent(in) :: h
        real(real64), intent(inout) :: kr
        real(real64), intent(inout), optional :: dkr_dh

        real(real64), parameter :: m_consistent_tolerance = 1.0d-9
        real(real64) :: zeta, beta_value, dkr(1)
        logical :: converged

        associate (params => self%parent%config)
            converged = .false.
            beta_value = 0.0d0
            if (h < 0.0d0 .and. .not. (abs(params%m1 - (1.0d0 - 1.0d0 / params%n1)) <= m_consistent_tolerance .or. &
                                         params%n1 <= 1.0d0)) then
                call evaluate_van_genuchten_mualem_argument(params%alpha1, params%n1, h, zeta)
                call self%incomplete_beta%evaluate(zeta, beta_value, converged)
            end if

            if (.not. (h < 0.0d0)) then
                ! A non-negative head is saturated: k_r = 1 on either branch.
                kr = 1.0d0
                if (present(dkr_dh)) dkr_dh = 0.0d0
            else if (converged) then
                if (present(dkr_dh)) then
                    ! The beta parameters initialize_hcf_base_vg set up.
                    call evaluate_kr_van_genuchten_mualem_tangent(params%alpha1, params%n1, params%m1, params%l, &
                                                                  params%m1 + 1.0d0 / params%n1, 1.0d0 - 1.0d0 / params%n1, &
                                                                  self%incomplete_beta%get_log_normalization(), &
                                                                  beta_value, h, [1.0d0], kr, dkr)
                    dkr_dh = dkr(1)
                else
                    call evaluate_kr_van_genuchten_mualem(params%alpha1, params%n1, params%m1, params%l, beta_value, h, kr)
                end if
            else if (present(dkr_dh)) then
                call evaluate_kr_van_genuchten_tangent(params%alpha1, params%n1, params%m1, params%l, h, [1.0d0], kr, dkr)
                dkr_dh = dkr(1)
            else
                call evaluate_kr_van_genuchten(params%alpha1, params%n1, params%m1, params%l, h, kr)
            end if
        end associate
    end subroutine calc_kr_base_vg

end submodule hcf_vg
