!>
!> @brief Implementation of GCC models.
!>
!> This submodule implements the methods for calculating suction and its derivatives
!> for both non-segregation and segregation GCC models.
!> Calculations are performed in Pascal [Pa].
!>
submodule(models_phase_change_gcc) gcc_base
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use :: constitutive_laws, only:evaluate_freezing_suction, evaluate_freezing_suction_tangent, &
        evaluate_freezing_suction_segregation, evaluate_freezing_suction_segregation_tangent
    implicit none

    real(real64), parameter :: MAX_REASONABLE_DENSITY = 1.0d6

contains

    !>
    !> @brief Initialize the GCC holder object.
    !>
    module subroutine initialize_holder_gccs(self, material_id, config, water, ice)
        implicit none
        class(holder_gccs), intent(inout) :: self
        integer(int32), intent(in) :: material_id
        class(type_config_gcc), intent(in) :: config
        type(type_iapws97), target, intent(in) :: water
        type(type_iapws06), target, intent(in) :: ice

        select case (config%gcc_model%ID)
        case (GCC_TYPES%NON_SEGREGATION%ID)
            allocate (type_gcc_non_segregation :: self%p)
        case (GCC_TYPES%SEGREGATION%ID)
            allocate (type_gcc_segregation :: self%p)
        end select

        if (allocated(self%p)) then
            call self%p%initialize(material_id, water, ice)
        end if
    end subroutine initialize_holder_gccs

    ! --------------------------------------------------------------------------
    ! Holder Wrapper Methods
    ! --------------------------------------------------------------------------

    module subroutine calc_holder(self, state, suction)
        implicit none
        class(holder_gccs), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: suction
        if (allocated(self%p)) then
            call self%p%calc(state, suction)
        else
            suction = 0.0d0
        end if
    end subroutine calc_holder

    module subroutine deriv_temperature_holder(self, state, deriv)
        implicit none
        class(holder_gccs), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: deriv
        if (allocated(self%p)) then
            call self%p%deriv_temperature(state, deriv)
        else
            deriv = 0.0d0
        end if
    end subroutine deriv_temperature_holder

    module subroutine deriv_pressure_holder(self, state, deriv)
        implicit none
        class(holder_gccs), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: deriv
        if (allocated(self%p)) then
            call self%p%deriv_pressure(state, deriv)
        else
            deriv = 0.0d0
        end if
    end subroutine deriv_pressure_holder

    !>
    !> @brief Initialize the abstract GCC base.
    !>
    module subroutine initialize_abst_gcc(self, material_id, water, ice)
        implicit none
        class(abst_gcc), intent(inout) :: self
        integer(int32), intent(in) :: material_id
        type(type_iapws97), target, intent(in) :: water
        type(type_iapws06), target, intent(in) :: ice

        self%material_id = material_id
        self%water => water
        self%ice => ice

        self%initialized = .true.
    end subroutine initialize_abst_gcc

    ! ==========================================================================
    ! Non-Segregation Implementation
    ! ==========================================================================

    module subroutine calc_gcc_nonseg(self, state, suction)
        implicit none
        class(type_gcc_non_segregation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: suction

        real(real64) :: temperature_K, rho_water

        if (.not. nonseg_state(self, state, temperature_K, rho_water)) then
            suction = 0.0d0
            return
        end if
        ! Keep the signed equilibrium suction above melting. The phase split
        ! selects the binding matric/freezing potential itself; clipping here
        ! makes its smooth maximum create ice in warm, nearly saturated soil.
        call evaluate_freezing_suction(lf, Tf0_K, max(temperature_K, 1.0d0), rho_water, suction)
    end subroutine calc_gcc_nonseg

    !> The absolute temperature and the liquid density the non-segregated
    !> suction is evaluated at, and whether they are usable.
    logical function nonseg_state(self, state, temperature_K, rho_water)
        implicit none
        class(type_gcc_non_segregation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: temperature_K, rho_water

        real(real64) :: temperature

        call state%temperature%get(temperature)
        call self%shift_temperature_absolute(temperature, temperature_K)
        call self%calc_rho_water(state, rho_water)
        nonseg_state = ieee_is_finite(temperature) .and. ieee_is_finite(temperature_K) .and. &
                       ieee_is_finite(rho_water) .and. abs(rho_water) <= MAX_REASONABLE_DENSITY
    end function nonseg_state

    module subroutine deriv_temp_gcc_nonseg(self, state, suction_derivative)
        implicit none
        class(type_gcc_non_segregation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: suction_derivative

        real(real64) :: temperature_K, rho_water, drho_dT, suction, dsuction(1)

        suction_derivative = 0.0d0
        if (.not. nonseg_state(self, state, temperature_K, rho_water)) return
        ! Below 1 K the suction is held at its 1 K value.
        if (temperature_K <= 1.0d0) return
        ! The liquid density is a function of the temperature too.
        call self%calc_drho_water_dT(state, drho_dT)
        call evaluate_freezing_suction_tangent(lf, Tf0_K, temperature_K, [1.0d0], rho_water, [drho_dT], &
                                               suction, dsuction)
        suction_derivative = dsuction(1)
    end subroutine deriv_temp_gcc_nonseg

    module subroutine deriv_pres_gcc_nonseg(self, state, suction_derivative)
        implicit none
        class(type_gcc_non_segregation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: suction_derivative

        real(real64) :: temperature_K, rho_water, drho_dP, suction, dsuction(1)

        suction_derivative = 0.0d0
        if (.not. nonseg_state(self, state, temperature_K, rho_water)) return
        ! The pressure reaches the suction only through the liquid density,
        ! which is held at atmospheric pressure for a negative gauge pressure.
        call self%calc_drho_water_dP(state, drho_dP)
        call evaluate_freezing_suction_tangent(lf, Tf0_K, max(temperature_K, 1.0d0), [0.0d0], rho_water, [drho_dP], &
                                               suction, dsuction)
        suction_derivative = dsuction(1)
    end subroutine deriv_pres_gcc_nonseg

    module subroutine deriv2_temp_gcc_nonseg(self, state, suction_derivative)
        implicit none
        class(type_gcc_non_segregation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: suction_derivative

        real(real64) :: temperature, temperature_K
        real(real64) :: rho_water
        real(real64) :: temperature_safe

        call state%temperature%get(temperature)
        call self%shift_temperature_absolute(temperature, temperature_K)
        call self%calc_rho_water(state, rho_water)

        if (.not. ieee_is_finite(temperature) .or. .not. ieee_is_finite(temperature_K) .or. &
            .not. ieee_is_finite(rho_water) .or. abs(rho_water) > MAX_REASONABLE_DENSITY) then
            suction_derivative = 0.0d0
            return
        end if

        if (.not. ieee_is_finite(temperature_K) .or. temperature_K <= 1.0d0) then
            temperature_safe = 1.0d0
        else
            temperature_safe = temperature_K
        end if
        suction_derivative = lf * rho_water / (temperature_safe * temperature_safe)
    end subroutine deriv2_temp_gcc_nonseg

    module subroutine deriv_pressure_ice_water_nonseg(self, state, deriv)
        implicit none
        class(type_gcc_non_segregation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: deriv

        ! Non-segregation model: dP_ice/dP_w = 0
        deriv = 0.0d0
    end subroutine deriv_pressure_ice_water_nonseg

    ! ==========================================================================
    ! Segregation Implementation
    ! ==========================================================================

    module subroutine calc_gcc_seg(self, state, suction)
        implicit none
        class(type_gcc_segregation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: suction

        real(real64) :: temperature, pressure, temperature_K, rho_water, rho_ice

        suction = 0.0d0
        if (.not. seg_state(self, state, temperature, pressure, temperature_K, rho_water, rho_ice)) return
        if (temperature > Tf0) return
        call evaluate_freezing_suction_segregation(lf, Tf0_K, pressure, max(temperature_K, 1.0d0), rho_water, rho_ice, &
                                                   suction)
    end subroutine calc_gcc_seg

    !> The state the segregated suction is evaluated at, and whether it is usable.
    logical function seg_state(self, state, temperature, pressure, temperature_K, rho_water, rho_ice)
        implicit none
        class(type_gcc_segregation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: temperature, pressure, temperature_K, rho_water, rho_ice

        call state%temperature%get(temperature)
        call state%pressure%get(pressure)
        call self%shift_temperature_absolute(temperature, temperature_K)
        call self%calc_rho_water(state, rho_water)
        call self%calc_rho_ice(state, rho_ice)
        seg_state = ieee_is_finite(temperature) .and. ieee_is_finite(temperature_K) .and. &
                    ieee_is_finite(pressure) .and. ieee_is_finite(rho_water) .and. ieee_is_finite(rho_ice) .and. &
                    abs(rho_water) <= MAX_REASONABLE_DENSITY .and. abs(rho_ice) <= MAX_REASONABLE_DENSITY
    end function seg_state

    module subroutine deriv_temp_gcc_seg(self, state, suction_derivative)
        implicit none
        class(type_gcc_segregation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: suction_derivative

        real(real64) :: temperature, pressure, temperature_K, rho_water, rho_ice
        real(real64) :: drho_water_dT, drho_ice_dT, suction, dsuction(1)

        suction_derivative = 0.0d0
        if (.not. seg_state(self, state, temperature, pressure, temperature_K, rho_water, rho_ice)) return
        if (temperature > Tf0 .or. temperature_K <= 1.0d0) return
        call self%calc_drho_water_dT(state, drho_water_dT)
        call self%calc_drho_ice_dT(state, drho_ice_dT)
        call evaluate_freezing_suction_segregation_tangent(lf, Tf0_K, pressure, [0.0d0], temperature_K, [1.0d0], &
                                                           rho_water, [drho_water_dT], rho_ice, [drho_ice_dT], &
                                                           suction, dsuction)
        suction_derivative = dsuction(1)
    end subroutine deriv_temp_gcc_seg

    module subroutine deriv_pres_gcc_seg(self, state, suction_derivative)
        implicit none
        class(type_gcc_segregation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: suction_derivative

        real(real64) :: temperature, pressure, temperature_K, rho_water, rho_ice
        real(real64) :: drho_water_dP, drho_ice_dP, suction, dsuction(1)

        ! Nonzero only in the frozen state: the suction is fixed at 0 when unfrozen.
        suction_derivative = 0.0d0
        if (.not. seg_state(self, state, temperature, pressure, temperature_K, rho_water, rho_ice)) return
        if (temperature > Tf0) return
        call self%calc_drho_water_dP(state, drho_water_dP)
        call self%calc_drho_ice_dP(state, drho_ice_dP)
        call evaluate_freezing_suction_segregation_tangent(lf, Tf0_K, pressure, [1.0d0], max(temperature_K, 1.0d0), &
                                                           [0.0d0], rho_water, [drho_water_dP], rho_ice, [drho_ice_dP], &
                                                           suction, dsuction)
        suction_derivative = dsuction(1)
    end subroutine deriv_pres_gcc_seg

    module subroutine deriv2_temp_gcc_seg(self, state, suction_derivative)
        implicit none
        class(type_gcc_segregation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: suction_derivative

        real(real64) :: temperature, temperature_K
        real(real64) :: rho_ice
        real(real64) :: temperature_safe

        call state%temperature%get(temperature)
        call self%shift_temperature_absolute(temperature, temperature_K)
        call self%calc_rho_ice(state, rho_ice)

        if (.not. ieee_is_finite(temperature) .or. .not. ieee_is_finite(temperature_K) .or. &
            .not. ieee_is_finite(rho_ice) .or. abs(rho_ice) > MAX_REASONABLE_DENSITY) then
            suction_derivative = 0.0d0
            return
        end if

        if (temperature <= Tf0) then
            if (.not. ieee_is_finite(temperature_K) .or. temperature_K <= 1.0d0) then
                temperature_safe = 1.0d0
            else
                temperature_safe = temperature_K
            end if
            suction_derivative = (lf * rho_ice / (temperature_safe * temperature_safe))
        else
            suction_derivative = 0.0d0
        end if
    end subroutine deriv2_temp_gcc_seg

    module subroutine deriv_pressure_ice_water_seg(self, state, deriv)
        implicit none
        class(type_gcc_segregation), intent(in) :: self
        type(type_state), intent(in) :: state
        real(real64), intent(inout) :: deriv

        real(real64) :: temperature
        real(real64) :: rho_water, rho_ice

        call state%temperature%get(temperature)

        ! Nonzero only in frozen state (zero when unfrozen since suction is fixed at 0)
        if (temperature <= Tf0) then
            call self%calc_rho_water(state, rho_water)
            call self%calc_rho_ice(state, rho_ice)
            ! dP_ice/dP_w = rho_ice / rho_water
            deriv = rho_ice / rho_water
        else
            deriv = 0.0d0
        end if
    end subroutine deriv_pressure_ice_water_seg

end submodule gcc_base
