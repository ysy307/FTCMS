!> Base module for constitutive modeling
!>
!> Connects state variables to thermodynamic property models (e.g., IAPWS).
module physics_constitutive_base
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use :: iapws, only:type_iapws97, type_iapws06, type_iapws_property
    use :: module_core, only:type_state
    use :: constitutive_constants, only:TtoK => celsius_to_kelvin, &
        P_atm => standard_atmospheric_pressure, &
        min_vapor_density, &
        reference_water_density, &
        water_triple_point_temperature, &
        water_critical_point_temperature, &
        iapws97_liquid_min_temperature
    implicit none
    private

    public :: abst_constitutive
    public :: type_iapws_wrapper

    !> Abstract base type for constitutive models
    type, abstract :: abst_constitutive
        !> Pointer to IAPWS-97 water model
        type(type_iapws97), pointer :: water
        !> Pointer to IAPWS-06 ice model
        type(type_iapws06), pointer :: ice
        !> Initialization status flag
        logical :: initialized = .false.
    contains
        ! ---- Algorithm / Operation ----
        !> Shift temperature to absolute scale
        procedure, pass(self), public :: shift_temperature_absolute => shift_temperature_absolute_abst_constitutive
        !> Shift pressure to absolute scale
        procedure, pass(self), public :: shift_pressure_absolute => shift_pressure_absolute_abst_constitutive
        !> Calculate liquid water density
        procedure, pass(self), public :: calc_rho_water => calc_rho_water_abst_constitutive
        !> Calculate temperature derivative of water density
        procedure, pass(self), public :: calc_drho_water_dT => calc_drho_water_dT_abst_constitutive
        !> Calculate pressure derivative of water density
        procedure, pass(self), public :: calc_drho_water_dP => calc_drho_water_dP_abst_constitutive
        !> Calculate specific heat capacity of water
        procedure, pass(self), public :: calc_cp_water => calc_cp_water_abst_constitutive
        procedure, pass(self), public :: calc_dcp_water_dT => calc_dcp_water_dT_abst_constitutive
        procedure, pass(self), public :: calc_dcp_water_dP => calc_dcp_water_dP_abst_constitutive
        !> Calculate ice density
        procedure, pass(self), public :: calc_rho_ice => calc_rho_ice_abst_constitutive
        !> Calculate temperature derivative of ice density
        procedure, pass(self), public :: calc_drho_ice_dT => calc_drho_ice_dT_abst_constitutive
        !> Calculate pressure derivative of ice density
        procedure, pass(self), public :: calc_drho_ice_dP => calc_drho_ice_dP_abst_constitutive
        !> Calculate specific heat capacity of ice
        procedure, pass(self), public :: calc_cp_ice => calc_cp_ice_abst_constitutive
        procedure, pass(self), public :: calc_dcp_ice_dT => calc_dcp_ice_dT_abst_constitutive
        procedure, pass(self), public :: calc_dcp_ice_dP => calc_dcp_ice_dP_abst_constitutive
        !> Calculate vapor density based on relative humidity
        procedure, pass(self), public :: calc_rho_vapor => calc_rho_vapor_abst_constitutive
        !> Calculate temperature derivative of vapor density
        procedure, pass(self), public :: calc_drho_vapor_dT => calc_drho_vapor_dT_abst_constitutive
        !> Calculate pressure derivative of vapor density
        procedure, pass(self), public :: calc_drho_vapor_dP => calc_drho_vapor_dP_abst_constitutive
        !> Calculate saturation vapor density
        procedure, pass(self), public :: calc_rho_vapor_saturation => calc_rho_vapor_saturation_abst_constitutive
        !> Calculate temperature derivative of saturation vapor density
        procedure, pass(self), public :: calc_drho_vapor_saturation_dT => calc_drho_vapor_saturation_dT_abst_constitutive
        !> Calculate pressure derivative of saturation vapor density
        procedure, pass(self), public :: calc_drho_vapor_saturation_dP => calc_drho_vapor_saturation_dP_abst_constitutive
        !> Calculate specific heat capacity of vapor
        procedure, pass(self), public :: calc_cp_vapor => calc_cp_vapor_abst_constitutive
        procedure, pass(self), public :: calc_dcp_vapor_dT => calc_dcp_vapor_dT_abst_constitutive
        procedure, pass(self), public :: calc_latent_heat_saturation => calc_latent_heat_saturation_abst_constitutive
        procedure, pass(self), public :: calc_dlatent_heat_saturation_dT => calc_dlatent_heat_saturation_dT_abst_constitutive

        ! ---- Inquiry ----
        !> Check if the model is initialized
        procedure, pass(self), public :: is_initialized => is_initialized_abst_constitutive

        ! ---- Getter ----
        !> Retrieve absolute temperature and pressure from state
        procedure, pass(self), private :: get_thermo_state_TP
        !> Retrieve absolute temperature from state
        procedure, pass(self), private :: get_thermo_state_T
        procedure, pass(self), private :: is_pressure_held

    end type abst_constitutive

    !> Concrete implementation using IAPWS models
    type, extends(abst_constitutive) :: type_iapws_wrapper
    contains
        ! ---- Lifecycle ----
        !> Initialize the IAPWS wrapper
        procedure, public :: initialize => initialize_iapws_wrapper
    end type type_iapws_wrapper

    !> One IAPWS evaluation, remembered with the absolute state it was made at.
    !>
    !> A closure evaluation asks for several properties of one phase at one
    !> state, and IAPWS forms each from the same Gibbs function; the memo
    !> answers the repeats from the first evaluation. The key is exact
    !> equality, so a remembered answer is the one IAPWS gives for that key.
    type :: type_iapws_memo
        !> Absolute temperature of the evaluation, K.
        real(real64) :: T_K = 0.0d0
        !> Absolute pressure of the evaluation, Pa; zero on the saturation line.
        real(real64) :: P_abs = 0.0d0
        !> Whether the memo holds an evaluation.
        logical :: filled = .false.
        !> A remembered scalar.
        real(real64) :: value = 0.0d0
        !> The Gibbs properties of a phase, for the water and ice memos.
        type(type_iapws_property) :: property
    end type type_iapws_memo

    !> The IAPWS-97 saturation quantities at one temperature, from one saturation evaluation.
    !>
    !> Every quantity comes from the same call, so the values are consistent
    !> whichever of them a model asks for first.
    type :: type_saturation_memo
        !> Absolute temperature of the evaluation, K.
        real(real64) :: T_K = 0.0d0
        !> Whether the memo holds an evaluation.
        logical :: filled = .false.
        !> \(\rho''\), kg/m3.
        real(real64) :: rho = 0.0d0
        !> \((\partial\rho''/\partial T)_p\), kg/(m3 K).
        real(real64) :: drho_dT = 0.0d0
        !> \((\partial\rho''/\partial p)_T\), kg/(m3 Pa).
        real(real64) :: drho_dP = 0.0d0
        !> \(\mathrm{d}\rho''/\mathrm{d}T\) along the saturation line, kg/(m3 K).
        real(real64) :: density_dT = 0.0d0
        !> \(c_p''\), J/(kg K).
        real(real64) :: cp = 0.0d0
        !> \(\mathrm{d}c_p''/\mathrm{d}T\) along the saturation line, J/(kg K2).
        real(real64) :: dcp_dT = 0.0d0
        !> \(L_v = h'' - h'\), J/kg.
        real(real64) :: latent_heat = 0.0d0
        !> \(\mathrm{d}L_v/\mathrm{d}T\) along the saturation line, J/(kg K).
        real(real64) :: dlatent_heat_dT = 0.0d0
    end type type_saturation_memo

    ! The last evaluation of each IAPWS quantity. The formulations are fixed,
    ! so the answer depends on the key alone, whichever model asked.
    type(type_iapws_memo) :: water_memo, ice_memo
    type(type_iapws_memo) :: water_dcp_dT_memo, water_dcp_dP_memo, ice_dcp_dT_memo, ice_dcp_dP_memo
    type(type_saturation_memo) :: saturation_memo

contains

    !> Whether `memo` holds the evaluation at \((T_K, P_{abs})\).
    pure logical function holds(memo, T_K, P_abs)
        implicit none
        type(type_iapws_memo), intent(in) :: memo
        real(real64), intent(in) :: T_K
        real(real64), intent(in) :: P_abs

        holds = memo%filled .and. memo%T_K == T_K .and. memo%P_abs == P_abs
    end function holds

    !> Make `water_memo` hold the IAPWS-97 properties of liquid water at \((T_K, P_{abs})\).
    subroutine refresh_water_memo(self, T_K, P_abs)
        implicit none
        class(abst_constitutive), intent(in) :: self
        real(real64), intent(in) :: T_K
        real(real64), intent(in) :: P_abs

        if (holds(water_memo, T_K, P_abs)) return
        call self%water%calc_properties(T_K, P_abs, water_memo%property)
        water_memo = type_iapws_memo(T_K, P_abs, .true., 0.0d0, water_memo%property)
    end subroutine refresh_water_memo

    !> Make `ice_memo` hold the IAPWS-06 properties of ice Ih at \((T_K, P_{abs})\).
    subroutine refresh_ice_memo(self, T_K, P_abs)
        implicit none
        class(abst_constitutive), intent(in) :: self
        real(real64), intent(in) :: T_K
        real(real64), intent(in) :: P_abs

        if (holds(ice_memo, T_K, P_abs)) return
        call self%ice%calc_properties(T_K, P_abs, ice_memo%property)
        ice_memo = type_iapws_memo(T_K, P_abs, .true., 0.0d0, ice_memo%property)
    end subroutine refresh_ice_memo

    !> Make the water heat-capacity derivative memos hold \((\partial c_p/\partial T, \partial c_p/\partial p)\)
    !> at \((T_K, P_{abs})\), both from one IAPWS-97 evaluation.
    subroutine refresh_water_dcp_memos(self, T_K, P_abs)
        implicit none
        class(abst_constitutive), intent(in) :: self
        real(real64), intent(in) :: T_K
        real(real64), intent(in) :: P_abs

        real(real64) :: dcp_dT, dcp_dP

        ! The two memos are only ever filled together.
        if (holds(water_dcp_dT_memo, T_K, P_abs)) return
        call self%water%calc_cp_derivatives(T_K, P_abs, dcp_dT, dcp_dP)
        water_dcp_dT_memo = type_iapws_memo(T_K, P_abs, .true., dcp_dT, water_dcp_dT_memo%property)
        water_dcp_dP_memo = type_iapws_memo(T_K, P_abs, .true., dcp_dP, water_dcp_dP_memo%property)
    end subroutine refresh_water_dcp_memos

    !> Make `saturation_memo` hold the IAPWS-97 saturation quantities at \(T_K\).
    subroutine refresh_saturation_memo(self, T_K)
        implicit none
        class(abst_constitutive), intent(in) :: self
        real(real64), intent(in) :: T_K

        if (saturation_memo%filled .and. saturation_memo%T_K == T_K) return
        call self%water%calc_saturation_properties(T_K, rho_vap=saturation_memo%rho, &
                                                   drho_dT_vap=saturation_memo%drho_dT, &
                                                   drho_dP_vap=saturation_memo%drho_dP, &
                                                   density_dT_vap=saturation_memo%density_dT, &
                                                   cp_vap=saturation_memo%cp, &
                                                   dcp_dT_vap=saturation_memo%dcp_dT, &
                                                   latent_heat=saturation_memo%latent_heat, &
                                                   dlatent_heat_dT=saturation_memo%dlatent_heat_dT)
        saturation_memo%T_K = T_K
        saturation_memo%filled = .true.
    end subroutine refresh_saturation_memo

    pure logical function is_valid_iapws97_liquid_state(T_K, P_abs)
        implicit none
        real(real64), intent(in) :: T_K
        real(real64), intent(in) :: P_abs

        is_valid_iapws97_liquid_state = .false.
        if (.not. ieee_is_finite(T_K) .or. .not. ieee_is_finite(P_abs)) return
        if (T_K < iapws97_liquid_min_temperature) return
        if (T_K > water_critical_point_temperature) return
        if (P_abs < 1.0d3) return
        if (P_abs > 1.0d8) return
        is_valid_iapws97_liquid_state = .true.
    end function is_valid_iapws97_liquid_state

    !> Shift temperature from Celsius to Kelvin
    !>
    !> Mathematical definition:
    !> - \( T = T_C + 273.15 \)
    !>
    !> Assumptions:
    !> - None
    !>
    !> Numerical guarantee:
    !> - Exact conversion
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine shift_temperature_absolute_abst_constitutive(self, temperature_degree, temperature_K)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Temperature in Celsius
        !> Not modified
        real(real64), intent(in) :: temperature_degree
        !> Temperature in Kelvin
        !> Overwritten on exit
        real(real64), intent(inout) :: temperature_K

        temperature_K = temperature_degree + TtoK
    end subroutine shift_temperature_absolute_abst_constitutive

    !> Shift pressure from gauge to absolute
    !>
    !> Mathematical definition:
    !> - \( P = P_{atm} + P_{gauge} \) if \( P_{gauge} \ge 0 \)
    !> - \( P = P_{atm} \) if \( P_{gauge} < 0 \)
    !>
    !> Assumptions:
    !> - Negative gauge pressure implies atmospheric conditions
    !>
    !> Numerical guarantee:
    !> - Exact conversion bounded by \( P_{atm} \)
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine shift_pressure_absolute_abst_constitutive(self, pressure_gauge, pressure_absolute)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Gauge pressure in Pascals
        !> Not modified
        real(real64), intent(in) :: pressure_gauge
        !> Absolute pressure in Pascals
        !> Overwritten on exit
        real(real64), intent(inout) :: pressure_absolute

        if (pressure_gauge < 0.0d0) then
            pressure_absolute = P_atm
        else
            pressure_absolute = P_atm + pressure_gauge
        end if
    end subroutine shift_pressure_absolute_abst_constitutive

    !> Get thermodynamic temperature and pressure
    !>
    !> Mathematical definition:
    !> - Extracts temperature and pressure from state object and converts to absolute
    !>
    !> Assumptions:
    !> - Fallback values used if state variables are not set
    !>
    !> Numerical guarantee:
    !> - No theoretical error bound available
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine get_thermo_state_TP(self, state, T_K, P_abs)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Temperature in Kelvin
        !> Overwritten on exit
        real(real64), intent(inout) :: T_K
        !> Absolute pressure in Pascals
        !> Overwritten on exit
        real(real64), intent(inout) :: P_abs

        real(real64) :: temp_c, press_g
        logical :: is_set

        call state%temperature%get(temp_c, is_set=is_set)
        if (is_set) then
            call self%shift_temperature_absolute(temp_c, T_K)
        else
            T_K = 273.15d0
        end if
        call state%pressure%get(press_g, is_set=is_set)
        if (is_set) then
            call self%shift_pressure_absolute(press_g, P_abs)
        else
            P_abs = P_atm
        end if
    end subroutine get_thermo_state_TP

    !> Get thermodynamic temperature
    !>
    !> Mathematical definition:
    !> - Extracts temperature from state object and converts to absolute
    !>
    !> Assumptions:
    !> - Fallback value used if state variable is not set
    !>
    !> Numerical guarantee:
    !> - No theoretical error bound available
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine get_thermo_state_T(self, state, T_K)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Temperature in Kelvin
        !> Overwritten on exit
        real(real64), intent(inout) :: T_K

        real(real64) :: temp_c
        logical :: is_set

        call state%temperature%get(temp_c, is_set=is_set)
        if (is_set) then
            call self%shift_temperature_absolute(temp_c, T_K)
        else
            T_K = 273.15d0
        end if
    end subroutine get_thermo_state_T

    !> Whether the properties are evaluated at a pressure that does not follow \(p_w\).
    !>
    !> `get_thermo_state_TP` holds the absolute pressure at atmospheric for a
    !> negative gauge pressure and for an unset one, so every pressure
    !> derivative of a property is zero there.
    logical function is_pressure_held(self, state)
        implicit none
        class(abst_constitutive), intent(in) :: self
        type(type_state), intent(in) :: state

        real(real64) :: pressure
        logical :: is_set

        call state%pressure%get(pressure, is_set=is_set)
        is_pressure_held = .true.
        if (.not. is_set) return
        is_pressure_held = pressure < 0.0d0
    end function is_pressure_held

    !> Calculate liquid water density
    !>
    !> Mathematical definition:
    !> - Uses IAPWS-97 formulation
    !>
    !> Assumptions:
    !> - IAPWS-97 also governs the supercooled liquid down to
    !>   iapws97_liquid_min_temperature (243.15 K), covering the unfrozen
    !>   water that persists in pore films below 0 degC; only outside that
    !>   range does the model fall back to a fixed reference value.
    !>
    !> Numerical guarantee:
    !> - Defined by IAPWS-97 standard
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_rho_water_abst_constitutive(self, state, density)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Density
        !> Overwritten on exit
        real(real64), intent(inout) :: density

        real(real64) :: T_K, P_abs

        call self%get_thermo_state_TP(state, T_K, P_abs)
        if (.not. is_valid_iapws97_liquid_state(T_K, P_abs)) then
            density = reference_water_density
            return
        end if
        call refresh_water_memo(self, T_K, P_abs)
        density = water_memo%property%rho
    end subroutine calc_rho_water_abst_constitutive

    !> Calculate temperature derivative of water density
    !>
    !> Mathematical definition:
    !> - \(\frac{\partial \rho}{\partial T}\)
    !>
    !> Assumptions:
    !> - IAPWS-97 also governs the supercooled liquid down to
    !>   iapws97_liquid_min_temperature (243.15 K); only outside that range
    !>   does the model fall back to a fixed (zero) tangent.
    !>
    !> Numerical guarantee:
    !> - Defined by IAPWS-97 standard
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_drho_water_dT_abst_constitutive(self, state, deriv_density)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Density derivative w.r.t temperature
        !> Overwritten on exit
        real(real64), intent(inout) :: deriv_density

        real(real64) :: T_K, P_abs

        call self%get_thermo_state_TP(state, T_K, P_abs)
        if (.not. is_valid_iapws97_liquid_state(T_K, P_abs)) then
            deriv_density = 0.0d0
            return
        end if
        call refresh_water_memo(self, T_K, P_abs)
        deriv_density = water_memo%property%drho_dT
    end subroutine calc_drho_water_dT_abst_constitutive

    !> Calculate pressure derivative of water density
    !>
    !> Mathematical definition:
    !> - \(\frac{\partial \rho}{\partial P}\)
    !>
    !> Assumptions:
    !> - IAPWS-97 also governs the supercooled liquid down to
    !>   iapws97_liquid_min_temperature (243.15 K); only outside that range
    !>   does the model fall back to a fixed (zero) tangent.
    !>
    !> Numerical guarantee:
    !> - Defined by IAPWS-97 standard
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_drho_water_dP_abst_constitutive(self, state, deriv_density)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Density derivative w.r.t pressure
        !> Overwritten on exit
        real(real64), intent(inout) :: deriv_density

        real(real64) :: T_K, P_abs

        call self%get_thermo_state_TP(state, T_K, P_abs)
        if (.not. is_valid_iapws97_liquid_state(T_K, P_abs) .or. self%is_pressure_held(state)) then
            deriv_density = 0.0d0
            return
        end if
        call refresh_water_memo(self, T_K, P_abs)
        deriv_density = water_memo%property%drho_dP
    end subroutine calc_drho_water_dP_abst_constitutive

    !> Calculate specific heat capacity of water
    !>
    !> Mathematical definition:
    !> - Uses IAPWS-97 formulation
    !>
    !> Assumptions:
    !> - IAPWS-97 also governs the supercooled liquid down to
    !>   iapws97_liquid_min_temperature (243.15 K); only outside that range
    !>   does the model fall back to a fixed reference value.
    !>
    !> Numerical guarantee:
    !> - Defined by IAPWS-97 standard
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_cp_water_abst_constitutive(self, state, cp)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Specific heat capacity
        !> Overwritten on exit
        real(real64), intent(inout) :: cp

        real(real64) :: T_K, P_abs

        call self%get_thermo_state_TP(state, T_K, P_abs)
        if (.not. is_valid_iapws97_liquid_state(T_K, P_abs)) then
            cp = 4181.3d0
            return
        end if
        call refresh_water_memo(self, T_K, P_abs)
        cp = water_memo%property%cp
    end subroutine calc_cp_water_abst_constitutive

    !> Temperature derivative of the specific heat of water.
    !>
    !> \[ \left(\frac{\partial c_p}{\partial T}\right)_P \]
    !>
    !> Assumptions: the same IAPWS-97 range as `calc_cp_water`; outside it the
    !> specific heat is a fixed reference value, so the derivative is zero.
    !> Numerical guarantees: defined by the IAPWS-97 Gibbs function.
    !> Computational complexity: \(O(1)\) memory and arithmetic.
    !> Failure behavior: returns without error.
    subroutine calc_dcp_water_dT_abst_constitutive(self, state, deriv_cp)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state; not modified
        type(type_state), intent(in) :: state
        !> \(\partial c_p/\partial T\), J/(kg K2); overwritten on exit
        real(real64), intent(inout) :: deriv_cp

        real(real64) :: T_K, P_abs

        call self%get_thermo_state_TP(state, T_K, P_abs)
        if (.not. is_valid_iapws97_liquid_state(T_K, P_abs)) then
            deriv_cp = 0.0d0
            return
        end if
        call refresh_water_dcp_memos(self, T_K, P_abs)
        deriv_cp = water_dcp_dT_memo%value
    end subroutine calc_dcp_water_dT_abst_constitutive

    !> Pore-pressure derivative of the specific heat of water.
    !>
    !> \[ \left(\frac{\partial c_p}{\partial p_w}\right)_T \]
    !>
    !> Assumptions: as `calc_dcp_water_dT`; the property is evaluated at
    !> atmospheric pressure while the gauge pressure is negative (see
    !> `shift_pressure_absolute`), where the derivative is therefore zero.
    !> Numerical guarantees: defined by the IAPWS-97 Gibbs function.
    !> Computational complexity: \(O(1)\) memory and arithmetic.
    !> Failure behavior: returns without error.
    subroutine calc_dcp_water_dP_abst_constitutive(self, state, deriv_cp)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state; not modified
        type(type_state), intent(in) :: state
        !> \(\partial c_p/\partial p_w\), J/(kg K Pa); overwritten on exit
        real(real64), intent(inout) :: deriv_cp

        real(real64) :: T_K, P_abs

        call self%get_thermo_state_TP(state, T_K, P_abs)
        if (.not. is_valid_iapws97_liquid_state(T_K, P_abs) .or. self%is_pressure_held(state)) then
            deriv_cp = 0.0d0
            return
        end if
        call refresh_water_dcp_memos(self, T_K, P_abs)
        deriv_cp = water_dcp_dP_memo%value
    end subroutine calc_dcp_water_dP_abst_constitutive

    !> Calculate ice density
    !>
    !> Mathematical definition:
    !> - Uses IAPWS-06 formulation
    !>
    !> Assumptions:
    !> - None
    !>
    !> Numerical guarantee:
    !> - Defined by IAPWS-06 standard
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_rho_ice_abst_constitutive(self, state, density)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Density
        !> Overwritten on exit
        real(real64), intent(inout) :: density

        real(real64) :: T_K, P_abs

        call self%get_thermo_state_TP(state, T_K, P_abs)
        call refresh_ice_memo(self, T_K, P_abs)
        density = ice_memo%property%rho
    end subroutine calc_rho_ice_abst_constitutive

    !> Calculate temperature derivative of ice density
    !>
    !> Mathematical definition:
    !> - \(\frac{\partial \rho}{\partial T}\)
    !>
    !> Assumptions:
    !> - None
    !>
    !> Numerical guarantee:
    !> - Defined by IAPWS-06 standard
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_drho_ice_dT_abst_constitutive(self, state, deriv_density)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Density derivative w.r.t temperature
        !> Overwritten on exit
        real(real64), intent(inout) :: deriv_density

        real(real64) :: T_K, P_abs

        call self%get_thermo_state_TP(state, T_K, P_abs)
        call refresh_ice_memo(self, T_K, P_abs)
        deriv_density = ice_memo%property%drho_dT
    end subroutine calc_drho_ice_dT_abst_constitutive

    !> Calculate pressure derivative of ice density
    !>
    !> Mathematical definition:
    !> - \(\frac{\partial \rho}{\partial P}\)
    !>
    !> Assumptions:
    !> - None
    !>
    !> Numerical guarantee:
    !> - Defined by IAPWS-06 standard
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_drho_ice_dP_abst_constitutive(self, state, deriv_density)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Density derivative w.r.t pressure
        !> Overwritten on exit
        real(real64), intent(inout) :: deriv_density

        real(real64) :: T_K, P_abs

        call self%get_thermo_state_TP(state, T_K, P_abs)
        if (self%is_pressure_held(state)) then
            deriv_density = 0.0d0
            return
        end if
        call refresh_ice_memo(self, T_K, P_abs)
        deriv_density = ice_memo%property%drho_dP
    end subroutine calc_drho_ice_dP_abst_constitutive

    !> Calculate specific heat capacity of ice
    !>
    !> Mathematical definition:
    !> - Uses IAPWS-06 formulation
    !>
    !> Assumptions:
    !> - None
    !>
    !> Numerical guarantee:
    !> - Defined by IAPWS-06 standard
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_cp_ice_abst_constitutive(self, state, cp)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Specific heat capacity
        !> Overwritten on exit
        real(real64), intent(inout) :: cp

        real(real64) :: T_K, P_abs

        call self%get_thermo_state_TP(state, T_K, P_abs)
        call refresh_ice_memo(self, T_K, P_abs)
        cp = ice_memo%property%cp
    end subroutine calc_cp_ice_abst_constitutive

    !> Temperature derivative of the specific heat of ice.
    !>
    !> \[ \left(\frac{\partial c_p}{\partial T}\right)_P \]
    !>
    !> Assumptions: none beyond the IAPWS-06 Ih formulation.
    !> Numerical guarantees: defined by the IAPWS-06 Gibbs function.
    !> Computational complexity: \(O(1)\) memory and arithmetic.
    !> Failure behavior: returns without error.
    subroutine calc_dcp_ice_dT_abst_constitutive(self, state, deriv_cp)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state; not modified
        type(type_state), intent(in) :: state
        !> \(\partial c_p/\partial T\), J/(kg K2); overwritten on exit
        real(real64), intent(inout) :: deriv_cp

        real(real64) :: T_K, P_abs

        call self%get_thermo_state_TP(state, T_K, P_abs)
        if (.not. holds(ice_dcp_dT_memo, T_K, P_abs)) then
            call self%ice%calc_dcp_dT(T_K, P_abs, deriv_cp)
            ice_dcp_dT_memo = type_iapws_memo(T_K, P_abs, .true., deriv_cp, ice_dcp_dT_memo%property)
        end if
        deriv_cp = ice_dcp_dT_memo%value
    end subroutine calc_dcp_ice_dT_abst_constitutive

    !> Pore-pressure derivative of the specific heat of ice.
    !>
    !> \[ \left(\frac{\partial c_p}{\partial p_w}\right)_T \]
    !>
    !> Assumptions: the property is evaluated at atmospheric pressure while the
    !> gauge pressure is negative, where the derivative is therefore zero.
    !> Numerical guarantees: defined by the IAPWS-06 Gibbs function.
    !> Computational complexity: \(O(1)\) memory and arithmetic.
    !> Failure behavior: returns without error.
    subroutine calc_dcp_ice_dP_abst_constitutive(self, state, deriv_cp)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state; not modified
        type(type_state), intent(in) :: state
        !> \(\partial c_p/\partial p_w\), J/(kg K Pa); overwritten on exit
        real(real64), intent(inout) :: deriv_cp

        real(real64) :: T_K, P_abs

        call self%get_thermo_state_TP(state, T_K, P_abs)
        if (self%is_pressure_held(state)) then
            deriv_cp = 0.0d0
            return
        end if
        if (.not. holds(ice_dcp_dP_memo, T_K, P_abs)) then
            call self%ice%calc_dcp_dP(T_K, P_abs, deriv_cp)
            ice_dcp_dP_memo = type_iapws_memo(T_K, P_abs, .true., deriv_cp, ice_dcp_dP_memo%property)
        end if
        deriv_cp = ice_dcp_dP_memo%value
    end subroutine calc_dcp_ice_dP_abst_constitutive

    !> Calculate vapor density
    !>
    !> Mathematical definition:
    !> - \(\rho_v = \max(\rho_{sat}(T) \cdot RH, \rho_{min})\)
    !>
    !> Assumptions:
    !> - Density is bounded by minimum vapor density
    !>
    !> Numerical guarantee:
    !> - No theoretical error bound available
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_rho_vapor_abst_constitutive(self, state, density)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Density
        !> Overwritten on exit
        real(real64), intent(inout) :: density

        real(real64) :: T_K, relative_humidity

        call self%get_thermo_state_T(state, T_K)
        T_K = min(max(T_K, water_triple_point_temperature), water_critical_point_temperature)
        call state%relative_humidity%get(relative_humidity)

        call refresh_saturation_memo(self, T_K)
        density = max(saturation_memo%rho * relative_humidity, min_vapor_density)
    end subroutine calc_rho_vapor_abst_constitutive

    !> Calculate temperature derivative of vapor density
    !>
    !> Mathematical definition:
    !> - \(\frac{\partial \rho_v}{\partial T} = \frac{\partial \rho_{sat}}{\partial T} \cdot RH\)
    !>
    !> Assumptions:
    !> - Relative humidity is constant with respect to temperature
    !>
    !> Numerical guarantee:
    !> - No theoretical error bound available
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_drho_vapor_dT_abst_constitutive(self, state, deriv_density)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Density derivative w.r.t temperature
        !> Overwritten on exit
        real(real64), intent(inout) :: deriv_density

        real(real64) :: T_K, relative_humidity, drho_sat_dT

        call self%get_thermo_state_T(state, T_K)
        T_K = min(max(T_K, water_triple_point_temperature), water_critical_point_temperature)
        call state%relative_humidity%get(relative_humidity)

        call refresh_saturation_memo(self, T_K)
        drho_sat_dT = saturation_memo%drho_dT

        deriv_density = drho_sat_dT * relative_humidity
    end subroutine calc_drho_vapor_dT_abst_constitutive

    !> Calculate pressure derivative of vapor density
    !>
    !> Mathematical definition:
    !> - \(\frac{\partial \rho_v}{\partial P} = 0\)
    !>
    !> Assumptions:
    !> - Saturation density depends only on temperature
    !>
    !> Numerical guarantee:
    !> - No theoretical error bound available
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_drho_vapor_dP_abst_constitutive(self, state, deriv_density)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Density derivative w.r.t pressure
        !> Overwritten on exit
        real(real64), intent(inout) :: deriv_density

        deriv_density = 0.0d0
    end subroutine calc_drho_vapor_dP_abst_constitutive

    !> Calculate saturation vapor density
    !>
    !> Mathematical definition:
    !> - Uses IAPWS-97 formulation
    !>
    !> Assumptions:
    !> - None
    !>
    !> Numerical guarantee:
    !> - Defined by IAPWS-97 standard
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_rho_vapor_saturation_abst_constitutive(self, state, density)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Density
        !> Overwritten on exit
        real(real64), intent(inout) :: density

        real(real64) :: T_K

        call self%get_thermo_state_T(state, T_K)
        T_K = min(max(T_K, water_triple_point_temperature), water_critical_point_temperature)
        call refresh_saturation_memo(self, T_K)
        density = saturation_memo%rho
    end subroutine calc_rho_vapor_saturation_abst_constitutive

    !> Calculate temperature derivative of saturation vapor density
    !>
    !> Mathematical definition:
    !> - \(\frac{\mathrm{d}\rho_{sat}}{\mathrm{d}T}\), the total derivative of
    !>   \(\rho''(T,p_s(T))\) along the saturation line, which is the curve
    !>   `calc_rho_vapor_saturation` evaluates on
    !>
    !> Assumptions:
    !> - Outside the triple-to-critical range the value is held at the bound,
    !>   so the derivative there is zero
    !>
    !> Numerical guarantee:
    !> - Defined by IAPWS-97 standard
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_drho_vapor_saturation_dT_abst_constitutive(self, state, deriv_density)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Density derivative w.r.t temperature
        !> Overwritten on exit
        real(real64), intent(inout) :: deriv_density

        real(real64) :: T_K

        call self%get_thermo_state_T(state, T_K)
        if (T_K < water_triple_point_temperature .or. T_K > water_critical_point_temperature) then
            deriv_density = 0.0d0
            return
        end if
        call refresh_saturation_memo(self, T_K)
        deriv_density = saturation_memo%density_dT
    end subroutine calc_drho_vapor_saturation_dT_abst_constitutive

    !> Calculate pressure derivative of saturation vapor density
    !>
    !> Mathematical definition:
    !> - \(\frac{\partial \rho_{sat}}{\partial P}\)
    !>
    !> Assumptions:
    !> - None
    !>
    !> Numerical guarantee:
    !> - Defined by IAPWS-97 standard
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_drho_vapor_saturation_dP_abst_constitutive(self, state, deriv_density)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Density derivative w.r.t pressure
        !> Overwritten on exit
        real(real64), intent(inout) :: deriv_density

        real(real64) :: T_K

        call self%get_thermo_state_T(state, T_K)
        T_K = min(max(T_K, water_triple_point_temperature), water_critical_point_temperature)
        call refresh_saturation_memo(self, T_K)
        deriv_density = saturation_memo%drho_dP
    end subroutine calc_drho_vapor_saturation_dP_abst_constitutive

    !> Calculate specific heat capacity of vapor
    !>
    !> Mathematical definition:
    !> - Uses IAPWS-97 formulation for saturation state
    !>
    !> Assumptions:
    !> - Vapor is at saturation pressure
    !>
    !> Numerical guarantee:
    !> - Defined by IAPWS-97 standard
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine calc_cp_vapor_abst_constitutive(self, state, cp)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state
        !> Not modified
        type(type_state), intent(in) :: state
        !> Specific heat capacity
        !> Overwritten on exit
        real(real64), intent(inout) :: cp

        real(real64) :: T_K

        call self%get_thermo_state_T(state, T_K)
        T_K = min(max(T_K, water_triple_point_temperature), water_critical_point_temperature)
        call refresh_saturation_memo(self, T_K)
        cp = saturation_memo%cp
    end subroutine calc_cp_vapor_abst_constitutive

    !> Temperature derivative of the specific heat of saturated vapour.
    !>
    !> \[ \frac{\mathrm{d}c_p''}{\mathrm{d}T}, \]
    !>
    !> the total derivative along the saturation line, which is the curve
    !> `calc_cp_vapor` evaluates on.
    !>
    !> Assumptions: outside the triple-to-critical range `calc_cp_vapor` holds
    !> the temperature at the bound, so the derivative there is zero.
    !> Numerical guarantees: defined by IAPWS-97 regions 2 and 4.
    !> Computational complexity: \(O(1)\) memory and arithmetic.
    !> Failure behavior: returns without error.
    subroutine calc_dcp_vapor_dT_abst_constitutive(self, state, deriv_cp)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Thermodynamic state; not modified
        type(type_state), intent(in) :: state
        !> \(\mathrm{d}c_p''/\mathrm{d}T\), J/(kg K2); overwritten on exit
        real(real64), intent(inout) :: deriv_cp

        real(real64) :: T_K

        call self%get_thermo_state_T(state, T_K)
        if (T_K < water_triple_point_temperature .or. T_K > water_critical_point_temperature) then
            deriv_cp = 0.0d0
            return
        end if
        call refresh_saturation_memo(self, T_K)
        deriv_cp = saturation_memo%dcp_dT
    end subroutine calc_dcp_vapor_dT_abst_constitutive

    !> Latent heat of vaporization on the IAPWS-97 saturation line.
    !>
    !> \[ L_v(T) = h''(T) - h'(T) \]
    !>
    !> Assumptions: \(T\) lies between the triple and the critical point.
    !> Numerical guarantees: that of IAPWS-97; no theoretical error bound available.
    !> Computational complexity: \(O(1)\) arithmetic and memory.
    !> Failure behavior: IAPWS returns NaN outside the saturation line.
    subroutine calc_latent_heat_saturation_abst_constitutive(self, T_K, latent_heat)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Absolute temperature, K, on the saturation line
        real(real64), intent(in) :: T_K
        !> \(L_v\), J/kg; overwritten on exit
        real(real64), intent(inout) :: latent_heat

        call refresh_saturation_memo(self, T_K)
        latent_heat = saturation_memo%latent_heat
    end subroutine calc_latent_heat_saturation_abst_constitutive

    !> Temperature derivative of the latent heat of vaporization along the saturation line.
    !>
    !> \[ \frac{\mathrm{d}L_v}{\mathrm{d}T} \]
    !>
    !> Assumptions: \(T\) lies between the triple and the critical point.
    !> Numerical guarantees: that of IAPWS-97; no theoretical error bound available.
    !> Computational complexity: \(O(1)\) arithmetic and memory.
    !> Failure behavior: IAPWS returns NaN outside the saturation line.
    subroutine calc_dlatent_heat_saturation_dT_abst_constitutive(self, T_K, deriv_latent_heat)
        implicit none
        !> Constitutive model object
        class(abst_constitutive), intent(in) :: self
        !> Absolute temperature, K, on the saturation line
        real(real64), intent(in) :: T_K
        !> \(\mathrm{d}L_v/\mathrm{d}T\), J/(kg K); overwritten on exit
        real(real64), intent(inout) :: deriv_latent_heat

        call refresh_saturation_memo(self, T_K)
        deriv_latent_heat = saturation_memo%dlatent_heat_dT
    end subroutine calc_dlatent_heat_saturation_dT_abst_constitutive

    !> Initialize the IAPWS wrapper
    !>
    !> Mathematical definition:
    !> - Associates pointers with IAPWS models
    !>
    !> Assumptions:
    !> - Provided models are valid targets
    !>
    !> Numerical guarantee:
    !> - No theoretical error bound available
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    subroutine initialize_iapws_wrapper(self, water_model, ice_model)
        implicit none
        !> IAPWS wrapper object
        class(type_iapws_wrapper), intent(inout) :: self
        !> IAPWS-97 water model
        !> Not modified
        type(type_iapws97), intent(in), target :: water_model
        !> IAPWS-06 ice model
        !> Not modified
        type(type_iapws06), intent(in), target :: ice_model

        self%water => water_model
        self%ice => ice_model
    end subroutine initialize_iapws_wrapper

    !> Check if the constitutive model is initialized
    !>
    !> Mathematical definition:
    !> - Returns the internal initialized flag
    !>
    !> Assumptions:
    !> - None
    !>
    !> Numerical guarantee:
    !> - No theoretical error bound available
    !>
    !> Computational complexity:
    !> - Memory: \(O(1)\)
    !> - Arithmetic: \(O(1)\)
    !>
    !> Failure behavior:
    !> - Returns without error
    function is_initialized_abst_constitutive(self) result(initialized)
        implicit none
        !> Constitutive model object
        !> Not modified
        class(abst_constitutive), intent(in) :: self
        !> Initialization status
        logical :: initialized

        initialized = self%initialized
    end function is_initialized_abst_constitutive

end module physics_constitutive_base
