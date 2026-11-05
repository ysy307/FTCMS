submodule(models_hcf) hcf_durner
    use :: constitutive_laws, only:evaluate_kr_durner, evaluate_kr_durner_tangent
    implicit none
contains
    module subroutine calc_kr_base_durner(self, h, kr, dkr_dh)
        implicit none
        class(type_hcf_base_durner), intent(in) :: self
        real(real64), intent(in) :: h
        real(real64), intent(inout) :: kr
        real(real64), intent(inout), optional :: dkr_dh

        real(real64) :: dkr(1)

        associate (params => self%parent%config)
            if (.not. (h < 0.0d0)) then
                ! A non-negative head is saturated: k_r = 1.
                kr = 1.0d0
                if (present(dkr_dh)) dkr_dh = 0.0d0
            else if (present(dkr_dh)) then
                call evaluate_kr_durner_tangent(params%alpha1, params%n1, params%m1, params%alpha2, params%n2, params%m2, &
                                                params%w1, params%w2, params%l, h, [1.0d0], kr, dkr)
                dkr_dh = dkr(1)
            else
                call evaluate_kr_durner(params%alpha1, params%n1, params%m1, params%alpha2, params%n2, params%m2, &
                                        params%w1, params%w2, params%l, h, kr)
            end if
        end associate
    end subroutine calc_kr_base_durner

end submodule hcf_durner
