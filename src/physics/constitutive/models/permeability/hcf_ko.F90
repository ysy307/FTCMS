submodule(models_hcf) hcf_ko
    use :: constitutive_laws, only:evaluate_kr_kosugi, evaluate_kr_kosugi_tangent
    implicit none
contains
    module subroutine calc_kr_base_ko(self, h, kr, dkr_dh)
        implicit none
        class(type_hcf_base_ko), intent(in) :: self
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
                call evaluate_kr_kosugi_tangent(params%alpha1, params%n1, params%l, h, [1.0d0], kr, dkr)
                dkr_dh = dkr(1)
            else
                call evaluate_kr_kosugi(params%alpha1, params%n1, params%l, h, kr)
            end if
        end associate
    end subroutine calc_kr_base_ko

end submodule hcf_ko
