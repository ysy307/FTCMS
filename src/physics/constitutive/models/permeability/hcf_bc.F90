submodule(models_hcf) hcf_bc
    use :: constitutive_laws, only:evaluate_kr_brooks_corey, evaluate_kr_brooks_corey_tangent
    implicit none
contains
    module subroutine calc_kr_base_bc(self, h, kr, dkr_dh)
        implicit none
        class(type_hcf_base_bc), intent(in) :: self
        real(real64), intent(in) :: h
        real(real64), intent(inout) :: kr
        real(real64), intent(inout), optional :: dkr_dh

        real(real64) :: dkr(1)

        associate (params => self%parent%config)
            if (.not. (h < params%alpha1)) then
                ! At or above the air-entry head the medium is saturated: k_r = 1.
                kr = 1.0d0
                if (present(dkr_dh)) dkr_dh = 0.0d0
            else if (present(dkr_dh)) then
                call evaluate_kr_brooks_corey_tangent(params%alpha1, params%n1, params%l, h, [1.0d0], kr, dkr)
                dkr_dh = dkr(1)
            else
                call evaluate_kr_brooks_corey(params%alpha1, params%n1, params%l, h, kr)
            end if
        end associate
    end subroutine calc_kr_base_bc

end submodule hcf_bc
