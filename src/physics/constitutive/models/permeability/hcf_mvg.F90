submodule(models_hcf) hcf_mvg
    use :: constitutive_laws, only:evaluate_kr_modified_van_genuchten, evaluate_kr_modified_van_genuchten_tangent
    implicit none
contains
    module subroutine calc_kr_base_mvg(self, h, kr, dkr_dh)
        implicit none
        class(type_hcf_base_mvg), intent(in) :: self
        real(real64), intent(in) :: h
        real(real64), intent(inout) :: kr
        real(real64), intent(inout), optional :: dkr_dh

        real(real64) :: dkr(1)

        associate (params => self%parent%config)
            if (.not. (h < params%h_crit)) then
                ! Above the critical head the modified curve is saturated: k_r = 1.
                kr = 1.0d0
                if (present(dkr_dh)) dkr_dh = 0.0d0
            else if (present(dkr_dh)) then
                call evaluate_kr_modified_van_genuchten_tangent(params%alpha1, params%n1, params%m1, params%l, &
                                                                params%h_crit, h, [1.0d0], kr, dkr)
                dkr_dh = dkr(1)
            else
                call evaluate_kr_modified_van_genuchten(params%alpha1, params%n1, params%m1, params%l, params%h_crit, h, kr)
            end if
        end associate
    end subroutine calc_kr_base_mvg

end submodule hcf_mvg
