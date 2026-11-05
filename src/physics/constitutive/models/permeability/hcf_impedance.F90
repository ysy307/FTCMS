submodule(models_hcf) hcf_impedance
    use :: constitutive_laws, only:evaluate_ice_impedance, evaluate_ice_impedance_tangent
    implicit none
contains

    module subroutine calc_impedance_exp(self, Qice, kr, dkr_dQ)
        implicit none
        class(type_hcf_impedance_exp), intent(in) :: self
        real(real64), intent(in) :: Qice
        real(real64), intent(inout) :: kr
        real(real64), intent(inout), optional :: dkr_dQ

        real(real64) :: dkr(1)

        associate (params => self%parent%config)
            if (present(dkr_dQ)) then
                call evaluate_ice_impedance_tangent(params%omega, Qice, [1.0d0], kr, dkr)
                dkr_dQ = dkr(1)
            else
                call evaluate_ice_impedance(params%omega, Qice, kr)
            end if
        end associate

    end subroutine calc_impedance_exp

end submodule hcf_impedance
