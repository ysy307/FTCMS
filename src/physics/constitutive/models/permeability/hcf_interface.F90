module models_hcf
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: iapws, only:type_iapws97, type_iapws06
    use :: module_core
    use :: numerical_special_functions, only: type_regularized_incomplete_beta
    use :: constitutive_constants, only:TtoK => celsius_to_kelvin, Mw => molar_mass_water, &
        Rg => universal_gas_constant, g => gravity_acceleration, rho_std => reference_water_density
    use :: physics_constitutive_base, only:abst_constitutive
    implicit none
    private

    public :: abst_hcf
    ! public :: type_config_hcf
    public :: holder_hcfs
    public :: type_hcf_base
    public :: type_hcf_impedance
    public :: type_hcf_viscosity
    public :: type_hcf_base_impedance
    public :: type_hcf_base_viscosity
    public :: type_hcf_impedance_viscosity
    public :: type_hcf_base_impedance_viscosity

    type :: holder_hcfs
        class(abst_hcf), allocatable :: p
    contains
        procedure, pass(self) :: initialize => initialize_holder_hcfs
    end type holder_hcfs

    interface
        module subroutine initialize_holder_hcfs(self, config, water, ice)
            implicit none
            class(holder_hcfs), intent(inout) :: self
            type(type_config_hcf), intent(in) :: config
            type(type_iapws97), intent(in), target :: water
            type(type_iapws06), intent(in), target :: ice

        end subroutine initialize_holder_hcfs
    end interface

    type, extends(abst_constitutive) :: type_hcf_vapor
        private
        class(abst_hcf), pointer :: parent
    contains
        procedure, pass(self), private :: calc_diffusivity => calc_diffusivity_vapor_in_air
        procedure, pass(self), private :: calc_tortuosity_factor => calc_tortuosity_factor_vapor
        procedure, pass(self), private :: calc_enhancement_factor => calc_enhancement_factor_vapor
        procedure, pass(self), public :: calc_Kvh => calc_Kvh_vapor
        procedure, pass(self), public :: calc_KvT => calc_KvT_vapor
    end type type_hcf_vapor

    interface
        module subroutine calc_diffusivity_vapor_in_air(self, temperature, Da)
            implicit none
            class(type_hcf_vapor), intent(in) :: self
            real(real64), intent(in) :: temperature
            real(real64), intent(inout) :: Da

        end subroutine calc_diffusivity_vapor_in_air

        module subroutine calc_tortuosity_factor_vapor(self, Qa, Qs, tau)
            implicit none
            class(type_hcf_vapor), intent(in) :: self
            real(real64), intent(in) :: Qa
            real(real64), intent(in) :: Qs
            real(real64), intent(inout) :: tau

        end subroutine calc_tortuosity_factor_vapor

        module subroutine calc_enhancement_factor_vapor(self, Qw, Qs, fc, eta)
            implicit none
            class(type_hcf_vapor), intent(in) :: self
            real(real64), intent(in) :: Qw
            real(real64), intent(in) :: Qs
            real(real64), intent(in) :: fc
            real(real64), intent(inout) :: eta

        end subroutine calc_enhancement_factor_vapor

        module subroutine calc_Kvh_vapor(self, state, Kvh)
            implicit none
            class(type_hcf_vapor), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: Kvh

        end subroutine calc_Kvh_vapor

        module subroutine calc_KvT_vapor(self, state, KvT)
            implicit none
            class(type_hcf_vapor), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: KvT

        end subroutine calc_KvT_vapor
    end interface

    type, extends(abst_constitutive), abstract :: abst_hcf
        private
        type(type_config_hcf) :: config
        class(abst_hcf_base), allocatable :: base
        class(abst_hcf_impedance), allocatable :: impedance
        class(abst_hcf_viscosity), allocatable :: viscosity
        type(type_hcf_vapor) :: vapor
    contains
        procedure, pass(self), public :: initialize => initialize_abst_hcf
        procedure(abst_calc_Kflh), pass(self), public, deferred :: calc_Kflh
        procedure, pass(self), public :: calc_KlT => calc_KlT_hcf
        procedure, pass(self), public :: calc_Kvh => calc_Kvh_hcf
        procedure, pass(self), public :: calc_KvT => calc_KvT_hcf
        procedure, pass(self), public :: calc_head => calc_head_hcf
        procedure, pass(self), private :: calc_impedance_ratio => calc_impedance_ratio_hcf
        procedure, pass(self), private :: apply_base_factor => apply_base_factor_hcf
        procedure, pass(self), private :: apply_impedance_factor => apply_impedance_factor_hcf
        procedure, pass(self), private :: apply_viscosity_factor => apply_viscosity_factor_hcf
        procedure, pass(self), public :: is_initialized => is_initialized_hcf
    end type abst_hcf

    interface
        module subroutine initialize_abst_hcf(self, config, water, ice)
            implicit none
            class(abst_hcf), intent(inout), target :: self
            type(type_config_hcf), intent(in) :: config
            type(type_iapws97), intent(in), target :: water
            type(type_iapws06), intent(in), target :: ice

        end subroutine initialize_abst_hcf

        module subroutine calc_KlT_hcf(self, state, KlT)
            implicit none
            class(abst_hcf), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: KlT

        end subroutine calc_KlT_hcf

        module subroutine calc_Kvh_hcf(self, state, Kvh)
            implicit none
            class(abst_hcf), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: Kvh

        end subroutine calc_Kvh_hcf

        module subroutine calc_KvT_hcf(self, state, KvT)
            implicit none
            class(abst_hcf), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: KvT

        end subroutine calc_KvT_hcf

        module subroutine calc_head_hcf(self, state, head, tangent)
            implicit none
            class(abst_hcf), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: head
            !> \((\partial h/\partial p_w,\partial h/\partial T)\), m/Pa and m/K; overwritten on exit when present
            real(real64), intent(inout), optional :: tangent(2)
        end subroutine calc_head_hcf

        module subroutine calc_impedance_ratio_hcf(self, state, ratio, tangent)
            implicit none
            class(abst_hcf), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: ratio
            !> \((\partial Q/\partial p_w,\partial Q/\partial T)\); overwritten on exit when present
            real(real64), intent(inout), optional :: tangent(2)
        end subroutine calc_impedance_ratio_hcf

        !> Multiply \(K\) by the retention factor \(k_r(h)\) and add \(\mathrm{d}\ln k_r\) to `log_tangent`.
        module subroutine apply_base_factor_hcf(self, state, kflh, log_tangent)
            implicit none
            class(abst_hcf), intent(in) :: self
            type(type_state), intent(in) :: state
            !> Conductivity so far, m/s; multiplied on exit
            real(real64), intent(inout) :: kflh
            !> \((\partial\ln K/\partial p_w,\partial\ln K/\partial T)\) so far; incremented on exit when present
            real(real64), intent(inout), optional :: log_tangent(2)
        end subroutine apply_base_factor_hcf

        !> Multiply \(K\) by the ice impedance \(10^{-\Omega Q}\) and add its log-tangent.
        module subroutine apply_impedance_factor_hcf(self, state, kflh, log_tangent)
            implicit none
            class(abst_hcf), intent(in) :: self
            type(type_state), intent(in) :: state
            !> Conductivity so far, m/s; multiplied on exit
            real(real64), intent(inout) :: kflh
            !> \((\partial\ln K/\partial p_w,\partial\ln K/\partial T)\) so far; incremented on exit when present
            real(real64), intent(inout), optional :: log_tangent(2)
        end subroutine apply_impedance_factor_hcf

        !> Multiply \(K\) by the viscosity ratio \(\mu_0/\mu(T)\) and add its log-tangent.
        module subroutine apply_viscosity_factor_hcf(self, state, kflh, log_tangent)
            implicit none
            class(abst_hcf), intent(in) :: self
            type(type_state), intent(in) :: state
            !> Conductivity so far, m/s; multiplied on exit
            real(real64), intent(inout) :: kflh
            !> \((\partial\ln K/\partial p_w,\partial\ln K/\partial T)\) so far; incremented on exit when present
            real(real64), intent(inout), optional :: log_tangent(2)
        end subroutine apply_viscosity_factor_hcf

        module pure function is_initialized_hcf(self) result(initialized)
            implicit none
            class(abst_hcf), intent(in) :: self
            logical :: initialized

        end function is_initialized_hcf
    end interface

    abstract interface
        !> Liquid hydraulic conductivity \(K\) of a state, m/s, and optionally
        !> \(\mathrm{d}\ln K\) in \((p_w,T)\) through the factors the mode applies.
        subroutine abst_calc_Kflh(self, state, Kflh, log_tangent)
            import :: abst_hcf, type_state, real64
            implicit none
            class(abst_hcf), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: Kflh
            !> \((\partial\ln K/\partial p_w,\partial\ln K/\partial T)\); overwritten on exit when present
            real(real64), intent(inout), optional :: log_tangent(2)
        end subroutine abst_calc_Kflh
    end interface

    type, extends(abst_hcf) :: type_hcf_base
    contains
        procedure :: calc_Kflh => calc_Kflh_base
    end type type_hcf_base

    type, extends(abst_hcf) :: type_hcf_impedance
    contains
        procedure :: calc_Kflh => calc_Kflh_impedance
    end type type_hcf_impedance

    type, extends(abst_hcf) :: type_hcf_viscosity
    contains
        procedure :: calc_Kflh => calc_Kflh_viscosity
    end type type_hcf_viscosity

    type, extends(abst_hcf) :: type_hcf_base_impedance
    contains
        procedure :: calc_Kflh => calc_Kflh_base_impedance
    end type type_hcf_base_impedance

    type, extends(abst_hcf) :: type_hcf_base_viscosity
    contains
        procedure :: calc_Kflh => calc_Kflh_base_viscosity
    end type type_hcf_base_viscosity

    type, extends(abst_hcf) :: type_hcf_impedance_viscosity
    contains
        procedure :: calc_Kflh => calc_Kflh_impedance_viscosity
    end type type_hcf_impedance_viscosity

    type, extends(abst_hcf) :: type_hcf_base_impedance_viscosity
    contains
        procedure :: calc_Kflh => calc_Kflh_base_impedance_viscosity
    end type type_hcf_base_impedance_viscosity

    interface
        module subroutine calc_Kflh_base(self, state, Kflh, log_tangent)
            implicit none
            class(type_hcf_base), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: Kflh
            real(real64), intent(inout), optional :: log_tangent(2)
        end subroutine calc_Kflh_base

        module subroutine calc_Kflh_impedance(self, state, Kflh, log_tangent)
            implicit none
            class(type_hcf_impedance), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: Kflh
            real(real64), intent(inout), optional :: log_tangent(2)
        end subroutine calc_Kflh_impedance

        module subroutine calc_Kflh_viscosity(self, state, Kflh, log_tangent)
            implicit none
            class(type_hcf_viscosity), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: Kflh
            real(real64), intent(inout), optional :: log_tangent(2)
        end subroutine calc_Kflh_viscosity

        module subroutine calc_Kflh_base_impedance(self, state, Kflh, log_tangent)
            implicit none
            class(type_hcf_base_impedance), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: Kflh
            real(real64), intent(inout), optional :: log_tangent(2)
        end subroutine calc_Kflh_base_impedance

        module subroutine calc_Kflh_base_viscosity(self, state, Kflh, log_tangent)
            implicit none
            class(type_hcf_base_viscosity), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: Kflh
            real(real64), intent(inout), optional :: log_tangent(2)
        end subroutine calc_Kflh_base_viscosity

        module subroutine calc_Kflh_impedance_viscosity(self, state, Kflh, log_tangent)
            implicit none
            class(type_hcf_impedance_viscosity), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: Kflh
            real(real64), intent(inout), optional :: log_tangent(2)
        end subroutine calc_Kflh_impedance_viscosity

        module subroutine calc_Kflh_base_impedance_viscosity(self, state, Kflh, log_tangent)
            implicit none
            class(type_hcf_base_impedance_viscosity), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: Kflh
            real(real64), intent(inout), optional :: log_tangent(2)
        end subroutine calc_Kflh_base_impedance_viscosity

    end interface

    type, abstract :: abst_hcf_base
        private
        class(abst_hcf), pointer :: parent
    contains
        procedure, pass(self), public :: initialize => initialize_hcf_base
        procedure(abst_calc_base_kr), pass(self), public, deferred :: calc_kr
    end type abst_hcf_base

    abstract interface
        !> Relative permeability \(k_r(h)\) of the retention model, and optionally \(\mathrm{d}k_r/\mathrm{d}h\).
        subroutine abst_calc_base_kr(self, h, kr, dkr_dh)
            import :: abst_hcf_base, real64
            implicit none
            class(abst_hcf_base), intent(in) :: self
            !> Pressure head \(h\), m
            real(real64), intent(in) :: h
            !> \(k_r\); overwritten on exit
            real(real64), intent(inout) :: kr
            !> \(\mathrm{d}k_r/\mathrm{d}h\), 1/m; overwritten on exit when present
            real(real64), intent(inout), optional :: dkr_dh
        end subroutine abst_calc_base_kr
    end interface

    type, extends(abst_hcf_base) :: type_hcf_base_bc
    contains
        procedure :: calc_kr => calc_kr_base_bc
    end type type_hcf_base_bc

    type, extends(abst_hcf_base) :: type_hcf_base_vg
        private
        type(type_regularized_incomplete_beta) :: incomplete_beta
    contains
        procedure :: initialize => initialize_hcf_base_vg
        procedure :: calc_kr => calc_kr_base_vg
    end type type_hcf_base_vg

    type, extends(abst_hcf_base) :: type_hcf_base_ko
    contains
        procedure :: calc_kr => calc_kr_base_ko
    end type type_hcf_base_ko

    type, extends(abst_hcf_base) :: type_hcf_base_mvg
    contains
        procedure :: calc_kr => calc_kr_base_mvg
    end type type_hcf_base_mvg

    type, extends(abst_hcf_base) :: type_hcf_base_durner
    contains
        procedure :: calc_kr => calc_kr_base_durner
    end type type_hcf_base_durner

    type, extends(abst_hcf_base) :: type_hcf_base_dvgch
    contains
        procedure :: calc_kr => calc_kr_base_dvgch
    end type type_hcf_base_dvgch

    interface
        module subroutine initialize_hcf_base(self)
            implicit none
            class(abst_hcf_base), intent(inout) :: self
        end subroutine initialize_hcf_base

        module subroutine initialize_hcf_base_vg(self)
            implicit none
            class(type_hcf_base_vg), intent(inout) :: self
        end subroutine initialize_hcf_base_vg

        module subroutine calc_kr_base_bc(self, h, kr, dkr_dh)
            implicit none
            class(type_hcf_base_bc), intent(in) :: self
            real(real64), intent(in) :: h
            real(real64), intent(inout) :: kr
            real(real64), intent(inout), optional :: dkr_dh
        end subroutine calc_kr_base_bc

        module subroutine calc_kr_base_vg(self, h, kr, dkr_dh)
            implicit none
            class(type_hcf_base_vg), intent(in) :: self
            real(real64), intent(in) :: h
            real(real64), intent(inout) :: kr
            real(real64), intent(inout), optional :: dkr_dh
        end subroutine calc_kr_base_vg

        module subroutine calc_kr_base_ko(self, h, kr, dkr_dh)
            implicit none
            class(type_hcf_base_ko), intent(in) :: self
            real(real64), intent(in) :: h
            real(real64), intent(inout) :: kr
            real(real64), intent(inout), optional :: dkr_dh
        end subroutine calc_kr_base_ko

        module subroutine calc_kr_base_mvg(self, h, kr, dkr_dh)
            implicit none
            class(type_hcf_base_mvg), intent(in) :: self
            real(real64), intent(in) :: h
            real(real64), intent(inout) :: kr
            real(real64), intent(inout), optional :: dkr_dh
        end subroutine calc_kr_base_mvg

        module subroutine calc_kr_base_durner(self, h, kr, dkr_dh)
            implicit none
            class(type_hcf_base_durner), intent(in) :: self
            real(real64), intent(in) :: h
            real(real64), intent(inout) :: kr
            real(real64), intent(inout), optional :: dkr_dh
        end subroutine calc_kr_base_durner

        module subroutine calc_kr_base_dvgch(self, h, kr, dkr_dh)
            implicit none
            class(type_hcf_base_dvgch), intent(in) :: self
            real(real64), intent(in) :: h
            real(real64), intent(inout) :: kr
            real(real64), intent(inout), optional :: dkr_dh
        end subroutine calc_kr_base_dvgch
    end interface

    type, abstract :: abst_hcf_impedance
        private
        class(abst_hcf), pointer :: parent
    contains
        procedure(abst_calc_impedance_kr), pass(self), public, deferred :: calc_impedance
    end type abst_hcf_impedance

    type, extends(abst_hcf_impedance) :: type_hcf_impedance_exp
    contains
        procedure :: calc_impedance => calc_impedance_exp
    end type type_hcf_impedance_exp

    abstract interface
        !> Impedance factor at an ice ratio \(Q\), and optionally its derivative in \(Q\).
        subroutine abst_calc_impedance_kr(self, Qice, kr, dkr_dQ)
            import :: abst_hcf_impedance, real64
            implicit none
            class(abst_hcf_impedance), intent(in) :: self
            !> Ice ratio \(Q\)
            real(real64), intent(in) :: Qice
            !> Impedance factor; overwritten on exit
            real(real64), intent(inout) :: kr
            !> Its derivative in \(Q\); overwritten on exit when present
            real(real64), intent(inout), optional :: dkr_dQ
        end subroutine abst_calc_impedance_kr
    end interface

    interface
        module subroutine calc_impedance_exp(self, Qice, kr, dkr_dQ)
            implicit none
            class(type_hcf_impedance_exp), intent(in) :: self
            real(real64), intent(in) :: Qice
            real(real64), intent(inout) :: kr
            real(real64), intent(inout), optional :: dkr_dQ
        end subroutine calc_impedance_exp
    end interface

    type, abstract :: abst_hcf_viscosity
        private
        real(real64) :: mu_zero
        class(abst_hcf), pointer :: parent
    contains
        procedure, pass(self), public :: initialize => initialize_abst_hcf_viscosity
        procedure(abst_calc_mu), pass(self), public, deferred :: calc_mu
        procedure, pass(self), public :: calc_viscosity => calc_kr_abst_hcf_viscosity
    end type abst_hcf_viscosity

    interface
        module subroutine initialize_abst_hcf_viscosity(self, temperature_critical)
            implicit none
            class(abst_hcf_viscosity), intent(inout) :: self
            real(real64), intent(in), optional :: temperature_critical

        end subroutine initialize_abst_hcf_viscosity
    end interface

    abstract interface
        !> Dynamic viscosity \(\mu(T)\), Pa s, and optionally \(\mathrm{d}\mu/\mathrm{d}T\).
        subroutine abst_calc_mu(self, temperature, mu, dmu_dT)
            import :: abst_hcf_viscosity, real64
            implicit none
            class(abst_hcf_viscosity), intent(in) :: self
            !> Temperature, C
            real(real64), intent(in) :: temperature
            !> \(\mu\); overwritten on exit
            real(real64), intent(inout) :: mu
            !> \(\mathrm{d}\mu/\mathrm{d}T\), Pa s/K; overwritten on exit when present
            real(real64), intent(inout), optional :: dmu_dT
        end subroutine abst_calc_mu
    end interface

    interface
        module subroutine calc_kr_abst_hcf_viscosity(self, temperature, kr, dkr_dT)
            implicit none
            class(abst_hcf_viscosity), intent(in) :: self
            real(real64), intent(in) :: temperature
            real(real64), intent(inout) :: kr
            !> \(\mathrm{d}(\mu_0/\mu)/\mathrm{d}T\), 1/K; overwritten on exit when present
            real(real64), intent(inout), optional :: dkr_dT
        end subroutine calc_kr_abst_hcf_viscosity
    end interface

    type, extends(abst_hcf_viscosity) :: type_hcf_viscosity_exp
    contains
        procedure, pass(self) :: calc_mu => calc_mu_exponential
    end type type_hcf_viscosity_exp

    type, extends(abst_hcf_viscosity) :: type_hcf_viscosity_supercool
    contains
        procedure, pass(self) :: calc_mu => calc_mu_exponential_supercooled
    end type type_hcf_viscosity_supercool

    interface

        module subroutine calc_mu_exponential(self, temperature, mu, dmu_dT)
            implicit none
            class(type_hcf_viscosity_exp), intent(in) :: self
            real(real64), intent(in) :: temperature
            real(real64), intent(inout) :: mu
            real(real64), intent(inout), optional :: dmu_dT
        end subroutine calc_mu_exponential

        module subroutine calc_mu_exponential_supercooled(self, temperature, mu, dmu_dT)
            implicit none
            class(type_hcf_viscosity_supercool), intent(in) :: self
            real(real64), intent(in) :: temperature
            real(real64), intent(inout) :: mu
            real(real64), intent(inout), optional :: dmu_dT
        end subroutine calc_mu_exponential_supercooled

        module subroutine calc_viscosity_exp(self, temperature, kr)
            implicit none
            class(type_hcf_viscosity_exp), intent(in) :: self
            real(real64), intent(in) :: temperature
            real(real64), intent(inout) :: kr

        end subroutine calc_viscosity_exp

        module subroutine calc_viscosity_supercool(self, temperature, kr)
            implicit none
            class(type_hcf_viscosity_supercool), intent(in) :: self
            real(real64), intent(in) :: temperature
            real(real64), intent(inout) :: kr

        end subroutine calc_viscosity_supercool

    end interface

end module models_hcf
