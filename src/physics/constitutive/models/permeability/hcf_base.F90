submodule(models_hcf) hcf_base
    use, intrinsic :: iso_fortran_env
    use :: constitutive_constants, only:tangent_pressure, tangent_temperature
    use :: constitutive_laws, only:evaluate_ice_impedance_ratio, evaluate_ice_impedance_ratio_tangent
    implicit none
    real(real64), parameter :: gamma_0 = 71.88875d0 ! g/s^2
contains

    !> Suction head [m] controlling the retention-based relative permeability.
    !>
    !> Relative permeability is evaluated at the same generalized suction that
    !> sets the unfrozen liquid content, so the mobile-water saturation it sees
    !> is the saturation actually present. Ice impedance is applied separately
    !> as a multiplicative factor by the configured HCF model. Falls back to the
    !> pore-water head when psi_eff is not set on the state (isothermal callers
    !> that never populate effective_suction).
    module subroutine calc_head_hcf(self, state, head, tangent)
        implicit none
        class(abst_hcf), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: head
        real(real64), intent(inout), optional :: tangent(2)

        real(real64) :: psi_eff, pressure
        logical :: is_set

        psi_eff = 0.0d0
        is_set = .false.
        call state%effective_suction%get(psi_eff, is_set)
        if (is_set) then
            head = -psi_eff / (rho_std * g)
            if (present(tangent)) then
                ! The split that set psi_eff published its tangents with it.
                call state%dSeff_dP%get(tangent(tangent_pressure))
                call state%dSeff_dT%get(tangent(tangent_temperature))
                tangent(:) = -tangent(:) / (rho_std * g)
            end if
        else
            call state%pressure%get(pressure)
            head = pressure / (rho_std * g)
            if (present(tangent)) then
                tangent(tangent_pressure) = 1.0d0 / (rho_std * g)
                tangent(tangent_temperature) = 0.0d0
            end if
        end if
    end subroutine calc_head_hcf

    !> Impedance factor of Hansson et al. (2004):
    !> \( Q = \theta_\mathrm{ice} / (\theta_\mathrm{T} - \theta_\mathrm{r}) \), the ice fraction of the
    !> total (minus residual) water content.  The impedance \(10^{-\Omega Q}\)
    !> with \(\Omega = 7\) was calibrated against the Mizoguchi columns using
    !> THIS definition; passing the raw volumetric ice content instead changes
    !> the meaning (and effective strength) of the calibrated factor.
    module subroutine calc_impedance_ratio_hcf(self, state, ratio, tangent)
        implicit none
        class(abst_hcf), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: ratio
        real(real64), intent(inout), optional :: tangent(2)

        real(real64) :: Qw, Qi, denom
        real(real64) :: dQw(2), dQi(2)

        call state%water_content%get(Qw)
        call state%ice_content%get(Qi)
        if (present(tangent)) tangent(:) = 0.0d0
        denom = Qw + Qi - self%config%theta_r
        if (.not. (denom > tiny(1.0d0) .and. Qi > 0.0d0)) then
            ratio = 0.0d0
            return
        end if

        if (present(tangent)) then
            call state%dQw_dP%get(dQw(tangent_pressure))
            call state%dQw_dT%get(dQw(tangent_temperature))
            call state%dQi_dP%get(dQi(tangent_pressure))
            call state%dQi_dT%get(dQi(tangent_temperature))
            call evaluate_ice_impedance_ratio_tangent(self%config%theta_r, Qw, dQw, Qi, dQi, ratio, tangent)
        else
            call evaluate_ice_impedance_ratio(self%config%theta_r, Qw, Qi, ratio)
        end if
        ! The cap is a branch of its own, constant in both unknowns.
        if (ratio > 1.0d0) then
            ratio = 1.0d0
            if (present(tangent)) tangent(:) = 0.0d0
        end if
    end subroutine calc_impedance_ratio_hcf

    module subroutine apply_base_factor_hcf(self, state, kflh, log_tangent)
        implicit none
        class(abst_hcf), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: kflh
        real(real64), intent(inout), optional :: log_tangent(2)

        real(real64) :: head, dhead(2), kr_base, dkr_dh

        if (present(log_tangent)) then
            call self%calc_head(state, head, dhead)
            call self%base%calc_kr(head, kr_base, dkr_dh)
            ! A vanishing k_r has no logarithm; the conductivity it gives is refused downstream.
            if (kr_base > 0.0d0) log_tangent(:) = log_tangent(:) + (dkr_dh / kr_base) * dhead(:)
        else
            call self%calc_head(state, head)
            call self%base%calc_kr(head, kr_base)
        end if
        kflh = kflh * kr_base
    end subroutine apply_base_factor_hcf

    module subroutine apply_impedance_factor_hcf(self, state, kflh, log_tangent)
        implicit none
        class(abst_hcf), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: kflh
        real(real64), intent(inout), optional :: log_tangent(2)

        real(real64) :: Q_ratio, dQ_ratio(2), kr_impedance, dkr_dQ

        if (present(log_tangent)) then
            call self%calc_impedance_ratio(state, Q_ratio, dQ_ratio)
            call self%impedance%calc_impedance(Q_ratio, kr_impedance, dkr_dQ)
            log_tangent(:) = log_tangent(:) + (dkr_dQ / kr_impedance) * dQ_ratio(:)
        else
            call self%calc_impedance_ratio(state, Q_ratio)
            call self%impedance%calc_impedance(Q_ratio, kr_impedance)
        end if
        kflh = kflh * kr_impedance
    end subroutine apply_impedance_factor_hcf

    module subroutine apply_viscosity_factor_hcf(self, state, kflh, log_tangent)
        implicit none
        class(abst_hcf), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: kflh
        real(real64), intent(inout), optional :: log_tangent(2)

        real(real64) :: temperature, kr_viscosity, dkr_dT

        call state%temperature%get(temperature)
        if (present(log_tangent)) then
            call self%viscosity%calc_viscosity(temperature, kr_viscosity, dkr_dT)
            log_tangent(tangent_temperature) = log_tangent(tangent_temperature) + dkr_dT / kr_viscosity
        else
            call self%viscosity%calc_viscosity(temperature, kr_viscosity)
        end if
        kflh = kflh * kr_viscosity
    end subroutine apply_viscosity_factor_hcf

    module subroutine initialize_holder_hcfs(self, config, water, ice)
        implicit none
        class(holder_hcfs), intent(inout) :: self
        type(type_config_hcf), intent(in) :: config
        type(type_iapws97), intent(in), target :: water
        type(type_iapws06), intent(in), target :: ice

        if (allocated(self%p)) then
            deallocate (self%p)
        end if

        if (config%model == HCF_MODES%BASE) then
            allocate (type_hcf_base :: self%p)
        else if (config%model == HCF_MODES%IMPEDANCE) then
            allocate (type_hcf_impedance :: self%p)
        else if (config%model == HCF_MODES%VISCOSITY) then
            allocate (type_hcf_viscosity :: self%p)
        else if (config%model == HCF_MODES%BASE_IMPEDANCE) then
            allocate (type_hcf_base_impedance :: self%p)
        else if (config%model == HCF_MODES%BASE_VISCOSITY) then
            allocate (type_hcf_base_viscosity :: self%p)
        else if (config%model == HCF_MODES%IMPEDANCE_VISCOSITY) then
            allocate (type_hcf_impedance_viscosity :: self%p)
        else if (config%model == HCF_MODES%BASE_IMPEDANCE_VISCOSITY) then
            allocate (type_hcf_base_impedance_viscosity :: self%p)
        else
            write (*, *) 'Error: Unknown HCF model ', config%model
            stop
        end if

        call self%p%initialize(config, water, ice)

    end subroutine initialize_holder_hcfs

    module subroutine initialize_abst_hcf(self, config, water, ice)
        implicit none
        class(abst_hcf), intent(inout), target :: self
        type(type_config_hcf), intent(in) :: config
        type(type_iapws97), intent(in), target :: water
        type(type_iapws06), intent(in), target :: ice

        call self%config%copy(config)

        if (self%config%model == HCF_MODES%BASE .or. &
            self%config%model == HCF_MODES%BASE_IMPEDANCE .or. &
            self%config%model == HCF_MODES%BASE_VISCOSITY .or. &
            self%config%model == HCF_MODES%BASE_IMPEDANCE_VISCOSITY) then

            if (self%config%swcc_model == SWCC_MODELS%BC) then
                allocate (type_hcf_base_bc :: self%base)
            else if (self%config%swcc_model == SWCC_MODELS%VG) then
                allocate (type_hcf_base_vg :: self%base)
            else if (self%config%swcc_model == SWCC_MODELS%KO) then
                allocate (type_hcf_base_ko :: self%base)
            else if (self%config%swcc_model == SWCC_MODELS%MVG) then
                allocate (type_hcf_base_mvg :: self%base)
            else if (self%config%swcc_model == SWCC_MODELS%DURNER) then
                allocate (type_hcf_base_durner :: self%base)
            else if (self%config%swcc_model == SWCC_MODELS%DVGCH) then
                allocate (type_hcf_base_dvgch :: self%base)
            else
                write (*, *) 'Error: Unknown WRF model ', self%config%swcc_model
                stop
            end if

            self%base%parent => self
            call self%base%initialize()
        end if

        if (self%config%model == HCF_MODES%IMPEDANCE .or. &
            self%config%model == HCF_MODES%BASE_IMPEDANCE .or. &
            self%config%model == HCF_MODES%IMPEDANCE_VISCOSITY .or. &
            self%config%model == HCF_MODES%BASE_IMPEDANCE_VISCOSITY) then
            allocate (type_hcf_impedance_exp :: self%impedance)
            self%impedance%parent => self
        end if

        if (self%config%model == HCF_MODES%VISCOSITY .or. &
            self%config%model == HCF_MODES%BASE_VISCOSITY .or. &
            self%config%model == HCF_MODES%IMPEDANCE_VISCOSITY .or. &
            self%config%model == HCF_MODES%BASE_IMPEDANCE_VISCOSITY) then

            if (self%config%water_viscosity_model == HCF_VISCOSITY_TYPES%EXPONENTIAL) then
                allocate (type_hcf_viscosity_exp :: self%viscosity)
            else if (self%config%water_viscosity_model == HCF_VISCOSITY_TYPES%SUPERCOOLED) then
                allocate (type_hcf_viscosity_supercool :: self%viscosity)
            else
                write (*, *) 'Error: Unknown water viscosity model ', self%config%water_viscosity_model
                stop
            end if
            call self%viscosity%initialize()
            self%viscosity%parent => self
        end if

        self%water => water
        self%ice => ice

        self%vapor%parent => self
        self%vapor%water => water
        self%vapor%ice => ice

        self%initialized = .true.
    end subroutine initialize_abst_hcf

    !> Default initialization for base permeability models without cached coefficients.
    module subroutine initialize_hcf_base(self)
        implicit none
        class(abst_hcf_base), intent(inout) :: self
    end subroutine initialize_hcf_base

    module pure function is_initialized_hcf(self) result(initialized)
        implicit none
        class(abst_hcf), intent(in) :: self
        logical :: initialized

        initialized = self%initialized
    end function is_initialized_hcf

    module subroutine calc_kflh_base(self, state, kflh, log_tangent)
        implicit none
        class(type_hcf_base), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: kflh
        real(real64), intent(inout), optional :: log_tangent(2)

        kflh = self%config%k_sat
        if (present(log_tangent)) log_tangent(:) = 0.0d0
        call self%apply_base_factor(state, kflh, log_tangent)

    end subroutine calc_kflh_base

    module subroutine calc_kflh_impedance(self, state, kflh, log_tangent)
        implicit none
        class(type_hcf_impedance), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: kflh
        real(real64), intent(inout), optional :: log_tangent(2)

        kflh = self%config%k_sat
        if (present(log_tangent)) log_tangent(:) = 0.0d0
        call self%apply_impedance_factor(state, kflh, log_tangent)

    end subroutine calc_kflh_impedance

    module subroutine calc_kflh_viscosity(self, state, kflh, log_tangent)
        implicit none
        class(type_hcf_viscosity), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: kflh
        real(real64), intent(inout), optional :: log_tangent(2)

        kflh = self%config%k_sat
        if (present(log_tangent)) log_tangent(:) = 0.0d0
        call self%apply_viscosity_factor(state, kflh, log_tangent)

    end subroutine calc_kflh_viscosity

    module subroutine calc_kflh_base_impedance(self, state, kflh, log_tangent)
        implicit none
        class(type_hcf_base_impedance), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: kflh
        real(real64), intent(inout), optional :: log_tangent(2)

        kflh = self%config%k_sat
        if (present(log_tangent)) log_tangent(:) = 0.0d0
        call self%apply_base_factor(state, kflh, log_tangent)
        call self%apply_impedance_factor(state, kflh, log_tangent)

    end subroutine calc_kflh_base_impedance

    module subroutine calc_kflh_base_viscosity(self, state, kflh, log_tangent)
        implicit none
        class(type_hcf_base_viscosity), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: kflh
        real(real64), intent(inout), optional :: log_tangent(2)

        kflh = self%config%k_sat
        if (present(log_tangent)) log_tangent(:) = 0.0d0
        call self%apply_base_factor(state, kflh, log_tangent)
        call self%apply_viscosity_factor(state, kflh, log_tangent)

    end subroutine calc_kflh_base_viscosity

    module subroutine calc_kflh_impedance_viscosity(self, state, kflh, log_tangent)
        implicit none
        class(type_hcf_impedance_viscosity), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: kflh
        real(real64), intent(inout), optional :: log_tangent(2)

        kflh = self%config%k_sat
        if (present(log_tangent)) log_tangent(:) = 0.0d0
        call self%apply_impedance_factor(state, kflh, log_tangent)
        call self%apply_viscosity_factor(state, kflh, log_tangent)

    end subroutine calc_kflh_impedance_viscosity

    module subroutine calc_kflh_base_impedance_viscosity(self, state, kflh, log_tangent)
        implicit none
        class(type_hcf_base_impedance_viscosity), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: kflh
        real(real64), intent(inout), optional :: log_tangent(2)

        kflh = self%config%k_sat
        if (present(log_tangent)) log_tangent(:) = 0.0d0
        call self%apply_base_factor(state, kflh, log_tangent)
        call self%apply_impedance_factor(state, kflh, log_tangent)
        call self%apply_viscosity_factor(state, kflh, log_tangent)

    end subroutine calc_kflh_base_impedance_viscosity

    module subroutine calc_klT_hcf(self, state, klT)
        implicit none
        class(abst_hcf), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: klT

        real(real64) :: K_lh, dgamma_dT
        real(real64) :: temperature
        real(real64) :: head

        call state%temperature%get(temperature)

        ! Actual pore-pressure head used by the retention scaling from which
        ! this thermo-osmotic term is derived.
        call self%calc_head(state, head)

        if (allocated(self%base)) then
            call self%calc_Kflh(state, K_lh)
            dgamma_dT = -0.1425d0 - 4.76d-4 * temperature
            klT = K_lh * head * self%config%gain_factor * (dgamma_dT / gamma_0)
        else
            klT = 0.0d0
        end if

    end subroutine calc_klT_hcf

    module subroutine calc_Kvh_hcf(self, state, Kvh)
        implicit none
        class(abst_hcf), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: Kvh

        call self%vapor%calc_Kvh(state, Kvh)

    end subroutine calc_Kvh_hcf

    module subroutine calc_KvT_hcf(self, state, KvT)
        implicit none
        class(abst_hcf), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: KvT

        call self%vapor%calc_KvT(state, KvT)

    end subroutine calc_KvT_hcf
end submodule hcf_base
