!> Base module for constitutive modeling
!>
!> Connects state variables to thermodynamic property models (e.g., IAPWS).
module physics_constitutive_base
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use :: iapws, only:type_iapws97, type_iapws06
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
        !> Calculate ice density
        procedure, pass(self), public :: calc_rho_ice => calc_rho_ice_abst_constitutive
        !> Calculate temperature derivative of ice density
        procedure, pass(self), public :: calc_drho_ice_dT => calc_drho_ice_dT_abst_constitutive
        !> Calculate pressure derivative of ice density
        procedure, pass(self), public :: calc_drho_ice_dP => calc_drho_ice_dP_abst_constitutive
        !> Calculate specific heat capacity of ice
        procedure, pass(self), public :: calc_cp_ice => calc_cp_ice_abst_constitutive
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

        ! ---- Inquiry ----
        !> Check if the model is initialized
        procedure, pass(self), public :: is_initialized => is_initialized_abst_constitutive

        ! ---- Getter ----
        !> Retrieve absolute temperature and pressure from state
        procedure, pass(self), private :: get_thermo_state_TP
        !> Retrieve absolute temperature from state
        procedure, pass(self), private :: get_thermo_state_T

    end type abst_constitutive

    !> Concrete implementation using IAPWS models
    type, extends(abst_constitutive) :: type_iapws_wrapper
    contains
        ! ---- Lifecycle ----
        !> Initialize the IAPWS wrapper
        procedure, public :: initialize => initialize_iapws_wrapper
    end type type_iapws_wrapper

contains

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
        call self%water%calc_rho(T_K, P_abs, density)
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
        call self%water%calc_drho_dT(T_K, P_abs, deriv_density)
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
        if (.not. is_valid_iapws97_liquid_state(T_K, P_abs)) then
            deriv_density = 0.0d0
            return
        end if
        call self%water%calc_drho_dP(T_K, P_abs, deriv_density)
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
        call self%water%calc_cp(T_K, P_abs, cp)
    end subroutine calc_cp_water_abst_constitutive

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
        call self%ice%calc_rho(T_K, P_abs, density)
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
        call self%ice%calc_drho_dT(T_K, P_abs, deriv_density)
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
        call self%ice%calc_drho_dP(T_K, P_abs, deriv_density)
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
        call self%ice%calc_cp(T_K, P_abs, cp)
    end subroutine calc_cp_ice_abst_constitutive

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

        call self%water%calc_saturation_density(T_K, density)
        density = max(density * relative_humidity, min_vapor_density)
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

        call self%water%calc_saturation_drho_dT(T_K, drho_sat_dT)

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
        call self%water%calc_saturation_density(T_K, density)
    end subroutine calc_rho_vapor_saturation_abst_constitutive

    !> Calculate temperature derivative of saturation vapor density
    !>
    !> Mathematical definition:
    !> - \(\frac{\partial \rho_{sat}}{\partial T}\)
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
        T_K = min(max(T_K, water_triple_point_temperature), water_critical_point_temperature)
        call self%water%calc_saturation_drho_dT(T_K, deriv_density)
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
        call self%water%calc_saturation_drho_dP(T_K, deriv_density)
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
        call self%water%calc_saturation_cp(T_K, cp)
    end subroutine calc_cp_vapor_abst_constitutive

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
