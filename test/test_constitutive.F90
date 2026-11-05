module test_constitutive_suite
    use, intrinsic :: iso_fortran_env, only: error_unit, int32, int64, output_unit, real64
    use :: iapws, only: type_iapws06, type_iapws97
    use :: module_core, only: GCC_TYPES, HCF_MODES, PHYSICS_UNITS, SWCC_MODELS, &
        type_config_gcc, type_config_hcf, type_config_wrf, type_state
    use :: module_input, only: input_translator, type_input
    use :: models_hcf, only: holder_hcfs
    use :: models_phase_change_vaporization, only: type_evaporation
    use :: models_phase_change_gcc, only: holder_gccs
    use :: models_phase_change_fusion, only: type_fusion, PORE_LIMIT_EXPELS_WATER
    use :: models_phase_change_manager, only: type_phase_manager
    use :: models_wrf, only: holder_wrfs
    use :: numerical_special_functions, only: type_regularized_incomplete_beta
    use :: constitutive_laws
    use :: physics_constitutive_base, only: type_iapws_wrapper
    implicit none
    private

    public :: type_constitutive_test_suite

    type :: type_constitutive_test_suite
        integer(int32), private :: failures = 0
    contains
        procedure, public :: run => run_constitutive_tests
        procedure, private :: test_translator_units
        procedure, private :: test_wrf_derivatives
        procedure, private :: test_pressure_capacity_units
        procedure, private :: test_special_functions
        procedure, private :: test_law_tangents
        procedure, private :: test_hcf_vg_dispatch
        procedure, private :: test_hcf_thermal_conductivity
        procedure, private :: test_freezing_storage_transport_split
        procedure, private :: test_vapor_combined_evaluation
        procedure, private :: test_iapws_memo
        procedure, private :: test_saturation_pressure
        procedure, private :: test_no_spurious_ice_at_warm_state
        procedure, private :: test_blend_derivative_sweep
        procedure, private :: test_liquid_pressure_driver
        procedure, private :: benchmark_hot_paths
        procedure, private :: report_benchmark
        procedure, private :: configure_wrf
        procedure, private :: check_close
        procedure, private :: check_true
    end type type_constitutive_test_suite

contains

    subroutine run_constitutive_tests(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        call self%test_translator_units()
        call self%test_wrf_derivatives()
        call self%test_pressure_capacity_units()
        call self%test_special_functions()
        call self%test_law_tangents()
        call self%test_hcf_vg_dispatch()
        call self%test_hcf_thermal_conductivity()
        call self%test_freezing_storage_transport_split()
        call self%test_vapor_combined_evaluation()
        call self%test_iapws_memo()
        call self%test_saturation_pressure()
        call self%test_no_spurious_ice_at_warm_state()
        call self%test_blend_derivative_sweep()
        call self%test_liquid_pressure_driver()
        call self%benchmark_hot_paths()

        if (self%failures > 0) then
            write (error_unit, '(A,I0)') "FAILED test_constitutive checks: ", self%failures
            error stop 1
        end if
        write (output_unit, '(A)') "All test_constitutive checks passed."
    end subroutine run_constitutive_tests

    subroutine test_translator_units(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        type(type_input) :: input
        type(type_config_wrf) :: config
        real(real64), parameter :: rho_g = 1000.0d0 * 9.80665d0

        input%basic%num_materials = 1
        allocate (input%basic%materials(1))

        input%basic%materials(1)%water_characteristic_curve%model_number = SWCC_MODELS%VG%ID
        input%basic%materials(1)%water_characteristic_curve%theta_s = 0.45d0
        input%basic%materials(1)%water_characteristic_curve%theta_r = 0.05d0
        input%basic%materials(1)%water_characteristic_curve%n1 = 2.0d0
        input%basic%materials(1)%water_characteristic_curve%m1 = 0.5d0
        input%basic%materials(1)%water_characteristic_curve%l = 0.5d0

        input%basic%materials(1)%water_characteristic_curve%unit = PHYSICS_UNITS%M%ID
        input%basic%materials(1)%water_characteristic_curve%alpha1 = 2.0d0
        call input_translator%execute(input, 1_int32, config)
        call self%check_close("VG alpha: m to 1/m", config%alpha1, 2.0d0, 1.0d-14)

        input%basic%materials(1)%water_characteristic_curve%unit = PHYSICS_UNITS%CM%ID
        input%basic%materials(1)%water_characteristic_curve%alpha1 = 0.02d0
        call input_translator%execute(input, 1_int32, config)
        call self%check_close("VG alpha: 1/cm to 1/m", config%alpha1, 2.0d0, 1.0d-14)

        input%basic%materials(1)%water_characteristic_curve%unit = PHYSICS_UNITS%PA%ID
        input%basic%materials(1)%water_characteristic_curve%alpha1 = 2.0d0 / rho_g
        call input_translator%execute(input, 1_int32, config)
        call self%check_close("VG alpha: 1/Pa to 1/m", config%alpha1, 2.0d0, 1.0d-12)

        input%basic%materials(1)%water_characteristic_curve%model_number = SWCC_MODELS%BC%ID
        input%basic%materials(1)%water_characteristic_curve%n1 = 2.0d0

        input%basic%materials(1)%water_characteristic_curve%unit = PHYSICS_UNITS%M%ID
        input%basic%materials(1)%water_characteristic_curve%alpha1 = -0.5d0
        call input_translator%execute(input, 1_int32, config)
        call self%check_close("BC air-entry: m to m", config%alpha1, -0.5d0, 1.0d-14)

        input%basic%materials(1)%water_characteristic_curve%unit = PHYSICS_UNITS%CM%ID
        input%basic%materials(1)%water_characteristic_curve%alpha1 = -50.0d0
        call input_translator%execute(input, 1_int32, config)
        call self%check_close("BC air-entry: cm to m", config%alpha1, -0.5d0, 1.0d-14)

        input%basic%materials(1)%water_characteristic_curve%unit = PHYSICS_UNITS%PA%ID
        input%basic%materials(1)%water_characteristic_curve%alpha1 = -0.5d0 * rho_g
        call input_translator%execute(input, 1_int32, config)
        call self%check_close("BC air-entry: Pa to m", config%alpha1, -0.5d0, 1.0d-12)

        input%basic%materials(1)%water_characteristic_curve%model_number = SWCC_MODELS%KO%ID
        input%basic%materials(1)%water_characteristic_curve%alpha1 = -0.25d0 * rho_g
        call input_translator%execute(input, 1_int32, config)
        call self%check_close("KO median head: Pa to m", config%alpha1, -0.25d0, 1.0d-12)

        input%basic%materials(1)%water_characteristic_curve%model_number = SWCC_MODELS%MVG%ID
        input%basic%materials(1)%water_characteristic_curve%alpha1 = 1.5d0 / rho_g
        input%basic%materials(1)%water_characteristic_curve%h_crit = -0.1d0 * rho_g
        call input_translator%execute(input, 1_int32, config)
        call self%check_close("MVG alpha: 1/Pa to 1/m", config%alpha1, 1.5d0, 1.0d-12)
        call self%check_close("MVG critical head: Pa to m", config%h_crit, -0.1d0, 1.0d-12)

        input%basic%materials(1)%water_characteristic_curve%model_number = SWCC_MODELS%DURNER%ID
        input%basic%materials(1)%water_characteristic_curve%alpha1 = 2.0d0 / rho_g
        input%basic%materials(1)%water_characteristic_curve%alpha2 = 0.4d0 / rho_g
        call input_translator%execute(input, 1_int32, config)
        call self%check_close("Durner alpha1: 1/Pa to 1/m", config%alpha1, 2.0d0, 1.0d-12)
        call self%check_close("Durner alpha2: 1/Pa to 1/m", config%alpha2, 0.4d0, 1.0d-12)

        input%basic%materials(1)%water_characteristic_curve%model_number = SWCC_MODELS%DVGCH%ID
        input%basic%materials(1)%water_characteristic_curve%alpha1 = 0.8d0 / rho_g
        call input_translator%execute(input, 1_int32, config)
        call self%check_close("DVGCH alpha: 1/Pa to 1/m", config%alpha1, 0.8d0, 1.0d-12)
    end subroutine test_translator_units

    subroutine test_wrf_derivatives(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        type(type_config_wrf) :: config
        type(holder_wrfs) :: wrf
        integer(int32) :: model_id
        real(real64) :: h, dh, theta_minus, theta_plus
        real(real64) :: derivative, derivative_fd
        character(len=64) :: label

        h = -1.0d0
        dh = 1.0d-6

        do model_id = SWCC_MODELS%BC%ID, SWCC_MODELS%DVGCH%ID
            call self%configure_wrf(config, model_id)
            call wrf%initialize(config)

            call wrf%p%calc(h - dh, theta_minus)
            call wrf%p%calc(h + dh, theta_plus)
            call wrf%p%deriv(h, derivative)
            derivative_fd = (theta_plus - theta_minus) / (2.0d0 * dh)

            write (label, '(A,I0)') "WRF dtheta/dh finite difference, model ", model_id
            call self%check_close(trim(label), derivative, derivative_fd, 2.0d-6)
            call self%check_true(trim(label)//" is non-negative", derivative >= 0.0d0)
        end do
    end subroutine test_wrf_derivatives

    subroutine test_pressure_capacity_units(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        type(type_config_wrf) :: config
        type(holder_wrfs) :: wrf
        real(real64), parameter :: rho_g = 1000.0d0 * 9.80665d0
        real(real64) :: h, pressure, dp
        real(real64) :: theta_minus, theta_plus
        real(real64) :: dtheta_dh, dtheta_dP, dtheta_dP_fd
        real(real64) :: physical_capacity_saturated, lscheme_capacity

        call self%configure_wrf(config, SWCC_MODELS%VG%ID)
        call wrf%initialize(config)

        h = -1.0d0
        pressure = rho_g * h
        dp = 1.0d-2
        call wrf%p%deriv(h, dtheta_dh)
        dtheta_dP = dtheta_dh / rho_g
        call wrf%p%calc((pressure - dp) / rho_g, theta_minus)
        call wrf%p%calc((pressure + dp) / rho_g, theta_plus)
        dtheta_dP_fd = (theta_plus - theta_minus) / (2.0d0 * dp)
        call self%check_close("dtheta/dP conversion from pressure head", dtheta_dP, dtheta_dP_fd, 2.0d-6)

        call wrf%p%deriv(0.0d0, dtheta_dh)
        physical_capacity_saturated = dtheta_dh / rho_g
        call wrf%p%calc_lscheme_capacity(lscheme_capacity)

        call self%check_close("physical saturated capacity", physical_capacity_saturated, 0.0d0, 1.0d-14)
        call self%check_true("L-scheme capacity is positive at saturation", lscheme_capacity > 0.0d0)
        write (output_unit, '(A,ES24.16)') "physical saturated dtheta/dP [1/Pa] = ", physical_capacity_saturated
        write (output_unit, '(A,ES24.16)') "L-scheme capacity floor [1/Pa]       = ", lscheme_capacity
    end subroutine test_pressure_capacity_units

    !> Every generated tangent against a central difference of its value routine.
    !>
    !> Checks the whole chain sympy.diff -> cse -> printed Fortran, one input
    !> direction at a time: a law with n inputs is called with the n unit
    !> tangents, so its output tangent is its gradient.
    subroutine test_law_tangents(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        real(real64), parameter :: tolerance = 1.0d-6
        type(type_regularized_incomplete_beta) :: beta
        real(real64) :: h, step, value, ahead, behind, zeta, beta_value, a, b
        real(real64) :: tangent(1), gradient(10), x(10), unit(10, 10)
        integer(int32) :: k

        unit(:, :) = 0.0d0
        do k = 1, 10
            unit(k, k) = 1.0d0
        end do

        ! Single-input laws: d/dh or d/dT.
        h = -0.7d0
        step = 1.0d-6
        call evaluate_kr_brooks_corey_tangent(-0.3d0, 2.0d0, 0.5d0, h, [1.0d0], value, tangent)
        call evaluate_kr_brooks_corey(-0.3d0, 2.0d0, 0.5d0, h + step, ahead)
        call evaluate_kr_brooks_corey(-0.3d0, 2.0d0, 0.5d0, h - step, behind)
        call self%check_close("Brooks-Corey kr tangent", tangent(1), (ahead - behind) / (2.0d0 * step), tolerance)

        call evaluate_kr_van_genuchten_tangent(1.11d0, 1.48d0, 1.0d0 - 1.0d0 / 1.48d0, 0.5d0, h, [1.0d0], value, tangent)
        call evaluate_kr_van_genuchten(1.11d0, 1.48d0, 1.0d0 - 1.0d0 / 1.48d0, 0.5d0, h + step, ahead)
        call evaluate_kr_van_genuchten(1.11d0, 1.48d0, 1.0d0 - 1.0d0 / 1.48d0, 0.5d0, h - step, behind)
        call self%check_close("van Genuchten kr tangent", tangent(1), (ahead - behind) / (2.0d0 * step), tolerance)

        ! Independent m: the incomplete beta is evaluated by the caller at each head.
        a = 0.2d0 + 1.0d0 / 1.48d0
        b = 1.0d0 - 1.0d0 / 1.48d0
        call beta%initialize(a, b)
        call evaluate_van_genuchten_mualem_argument(1.11d0, 1.48d0, h, zeta)
        call beta%evaluate(zeta, beta_value)
        call evaluate_kr_van_genuchten_mualem_tangent(1.11d0, 1.48d0, 0.2d0, 0.5d0, a, b, beta%get_log_normalization(), &
                                                      beta_value, h, [1.0d0], value, tangent)
        call evaluate_van_genuchten_mualem_argument(1.11d0, 1.48d0, h + step, zeta)
        call beta%evaluate(zeta, beta_value)
        call evaluate_kr_van_genuchten_mualem(1.11d0, 1.48d0, 0.2d0, 0.5d0, beta_value, h + step, ahead)
        call evaluate_van_genuchten_mualem_argument(1.11d0, 1.48d0, h - step, zeta)
        call beta%evaluate(zeta, beta_value)
        call evaluate_kr_van_genuchten_mualem(1.11d0, 1.48d0, 0.2d0, 0.5d0, beta_value, h - step, behind)
        call self%check_close("van Genuchten-Mualem kr tangent", tangent(1), (ahead - behind) / (2.0d0 * step), tolerance)

        call evaluate_kr_modified_van_genuchten_tangent(1.11d0, 1.48d0, 0.32d0, 0.5d0, -0.02d0, h, [1.0d0], value, tangent)
        call evaluate_kr_modified_van_genuchten(1.11d0, 1.48d0, 0.32d0, 0.5d0, -0.02d0, h + step, ahead)
        call evaluate_kr_modified_van_genuchten(1.11d0, 1.48d0, 0.32d0, 0.5d0, -0.02d0, h - step, behind)
        call self%check_close("modified van Genuchten kr tangent", tangent(1), (ahead - behind) / (2.0d0 * step), &
                              tolerance)

        call evaluate_kr_kosugi_tangent(-0.5d0, 1.2d0, 0.5d0, h, [1.0d0], value, tangent)
        call evaluate_kr_kosugi(-0.5d0, 1.2d0, 0.5d0, h + step, ahead)
        call evaluate_kr_kosugi(-0.5d0, 1.2d0, 0.5d0, h - step, behind)
        call self%check_close("Kosugi kr tangent", tangent(1), (ahead - behind) / (2.0d0 * step), tolerance)

        call evaluate_kr_durner_tangent(1.1d0, 1.5d0, 0.33d0, 5.0d0, 2.0d0, 0.5d0, 0.6d0, 0.4d0, 0.5d0, h, [1.0d0], &
                                        value, tangent)
        call evaluate_kr_durner(1.1d0, 1.5d0, 0.33d0, 5.0d0, 2.0d0, 0.5d0, 0.6d0, 0.4d0, 0.5d0, h + step, ahead)
        call evaluate_kr_durner(1.1d0, 1.5d0, 0.33d0, 5.0d0, 2.0d0, 0.5d0, 0.6d0, 0.4d0, 0.5d0, h - step, behind)
        call self%check_close("Durner kr tangent", tangent(1), (ahead - behind) / (2.0d0 * step), tolerance)

        call evaluate_water_viscosity_exponential_tangent(-2.0d0, [1.0d0], value, tangent)
        call evaluate_water_viscosity_exponential(-2.0d0 + step, ahead)
        call evaluate_water_viscosity_exponential(-2.0d0 - step, behind)
        call self%check_close("exponential viscosity tangent", tangent(1), (ahead - behind) / (2.0d0 * step), tolerance)

        call evaluate_water_viscosity_supercooled_tangent(-2.0d0, [1.0d0], value, tangent)
        call evaluate_water_viscosity_supercooled(-2.0d0 + step, ahead)
        call evaluate_water_viscosity_supercooled(-2.0d0 - step, behind)
        call self%check_close("supercooled viscosity tangent", tangent(1), (ahead - behind) / (2.0d0 * step), tolerance)

        call evaluate_ice_impedance_tangent(7.0d0, 0.3d0, [1.0d0], value, tangent)
        call evaluate_ice_impedance(7.0d0, 0.3d0 + step, ahead)
        call evaluate_ice_impedance(7.0d0, 0.3d0 - step, behind)
        call self%check_close("ice impedance tangent", tangent(1), (ahead - behind) / (2.0d0 * step), tolerance)

        ! Multi-input laws: the gradient, one unit direction per input.
        x(1:2) = [0.2d0, 0.1d0]
        call evaluate_ice_impedance_ratio_tangent(0.05d0, x(1), unit(1:2, 1), x(2), unit(1:2, 2), value, gradient(1:2))
        do k = 1, 2
            call evaluate_ice_impedance_ratio(0.05d0, x(1) + step * unit(1, k), x(2) + step * unit(2, k), ahead)
            call evaluate_ice_impedance_ratio(0.05d0, x(1) - step * unit(1, k), x(2) - step * unit(2, k), behind)
            call self%check_close("ice ratio gradient", gradient(k), (ahead - behind) / (2.0d0 * step), tolerance)
        end do

        x(1:4) = [0.465d0, 0.2d0, 0.15d0, 0.185d0]
        call evaluate_thermal_conductivity_geometric4_tangent(3.78d0, 0.6d0, 2.2d0, 0.0241d0, &
                                                              x(1), unit(1:4, 1), x(2), unit(1:4, 2), &
                                                              x(3), unit(1:4, 3), x(4), unit(1:4, 4), value, gradient(1:4))
        do k = 1, 4
            call evaluate_thermal_conductivity_geometric4(3.78d0, 0.6d0, 2.2d0, 0.0241d0, &
                                                          x(1) + step * unit(1, k), x(2) + step * unit(2, k), &
                                                          x(3) + step * unit(3, k), x(4) + step * unit(4, k), ahead)
            call evaluate_thermal_conductivity_geometric4(3.78d0, 0.6d0, 2.2d0, 0.0241d0, &
                                                          x(1) - step * unit(1, k), x(2) - step * unit(2, k), &
                                                          x(3) - step * unit(3, k), x(4) - step * unit(4, k), behind)
            call self%check_close("geometric conductivity gradient", gradient(k), (ahead - behind) / (2.0d0 * step), &
                                  tolerance)
        end do

        x(1:2) = [0.2d0, 0.15d0]
        call evaluate_thermal_conductivity_campbell_hansson_tangent(0.55d0, 0.8d0, 3.07d0, 0.13d0, 13.05d0, 1.06d0, &
                                                                    x(1), unit(1:2, 1), x(2), unit(1:2, 2), &
                                                                    value, gradient(1:2))
        do k = 1, 2
            call evaluate_thermal_conductivity_campbell_hansson(0.55d0, 0.8d0, 3.07d0, 0.13d0, 13.05d0, 1.06d0, &
                                                                x(1) + step * unit(1, k), x(2) + step * unit(2, k), ahead)
            call evaluate_thermal_conductivity_campbell_hansson(0.55d0, 0.8d0, 3.07d0, 0.13d0, 13.05d0, 1.06d0, &
                                                                x(1) - step * unit(1, k), x(2) - step * unit(2, k), behind)
            call self%check_close("Campbell-Hansson gradient", gradient(k), (ahead - behind) / (2.0d0 * step), tolerance)
        end do

        ! T, theta_w, theta_i, theta_v, rho_w, c_w, rho_i, c_i, c_v, L_v.
        x(:) = [-1.0d0, 0.2d0, 0.15d0, 1.0d-6, 999.8d0, 4217.0d0, 916.7d0, 2093.0d0, 1859.0d0, 2.5d6]
        call evaluate_enthalpy_density_tangent(2.5536d6 * 0.465d0, 334560.0d0, x(1), unit(:, 1), x(2), unit(:, 2), &
                                               x(3), unit(:, 3), x(4), unit(:, 4), x(5), unit(:, 5), x(6), unit(:, 6), &
                                               x(7), unit(:, 7), x(8), unit(:, 8), x(9), unit(:, 9), x(10), unit(:, 10), &
                                               value, gradient)
        do k = 1, 10
            ! U is linear in each input, so a central difference is exact at any
            ! step; a large one keeps it clear of the cancellation in U itself.
            step = 1.0d-1 * max(1.0d0, abs(x(k)))
            associate (p => x + step * unit(:, k), m => x - step * unit(:, k))
                call evaluate_enthalpy_density(2.5536d6 * 0.465d0, 334560.0d0, p(1), p(2), p(3), p(4), p(5), p(6), &
                                               p(7), p(8), p(9), p(10), ahead)
                call evaluate_enthalpy_density(2.5536d6 * 0.465d0, 334560.0d0, m(1), m(2), m(3), m(4), m(5), m(6), &
                                               m(7), m(8), m(9), m(10), behind)
            end associate
            call self%check_close("enthalpy gradient", gradient(k), (ahead - behind) / (2.0d0 * step), tolerance)
        end do

        ! theta_w, theta_i, theta_v, rho_w, rho_i, S_c.
        x(1:6) = [0.2d0, 0.15d0, 1.0d-6, 999.8d0, 916.7d0, 1.0d-7]
        call evaluate_total_stored_water_tangent(x(1), unit(1:6, 1), x(2), unit(1:6, 2), x(3), unit(1:6, 3), &
                                                 x(4), unit(1:6, 4), x(5), unit(1:6, 5), x(6), unit(1:6, 6), &
                                                 value, gradient(1:6))
        do k = 1, 6
            step = 1.0d-6 * max(1.0d0, abs(x(k)))
            associate (p => x(1:6) + step * unit(1:6, k), m => x(1:6) - step * unit(1:6, k))
                call evaluate_total_stored_water(p(1), p(2), p(3), p(4), p(5), p(6), ahead)
                call evaluate_total_stored_water(m(1), m(2), m(3), m(4), m(5), m(6), behind)
            end associate
            call self%check_close("stored water gradient", gradient(k), (ahead - behind) / (2.0d0 * step), tolerance)
        end do
    end subroutine test_law_tangents

    subroutine test_special_functions(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        type(type_regularized_incomplete_beta) :: incomplete_beta
        real(real64) :: value
        logical :: converged

        value = 0.0d0
        converged = .false.
        call incomplete_beta%initialize(2.0d0, 3.0d0)
        call incomplete_beta%evaluate(0.5d0, value, converged)
        call self%check_true("incomplete beta continued fraction converged", converged)
        call self%check_close("regularized incomplete beta I_0.5(2,3)", value, 0.6875d0, 2.0d-14)

        value = 0.0d0
        converged = .false.
        call incomplete_beta%initialize(0.5d0, 0.5d0)
        call incomplete_beta%evaluate(0.5d0, value, converged)
        call self%check_true("incomplete beta symmetric case converged", converged)
        call self%check_close("regularized incomplete beta I_0.5(0.5,0.5)", value, 0.5d0, 2.0d-14)
    end subroutine test_special_functions

    subroutine test_hcf_vg_dispatch(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        type(type_config_hcf) :: config
        type(type_state) :: state
        type(holder_hcfs) :: hcf
        type(type_iapws97), target :: water
        type(type_iapws06), target :: ice
        type(type_regularized_incomplete_beta) :: incomplete_beta
        real(real64), parameter :: rho_g = 1000.0d0 * 9.80665d0
        real(real64) :: effective_saturation, zeta, beta_value
        real(real64) :: conductivity, expected_conductivity
        logical :: converged

        call config%reset()
        config%model = HCF_MODES%BASE
        config%swcc_model = SWCC_MODELS%VG
        config%k_sat = 1.0d0
        config%l = 0.5d0
        config%alpha1 = 2.0d0
        config%n1 = 1.48d0
        config%m1 = 0.2d0
        call hcf%initialize(config, water, ice)

        call state%pressure%set(-rho_g)
        conductivity = 0.0d0
        call hcf%p%calc_Kflh(state, conductivity)

        effective_saturation = (1.0d0 + config%alpha1**config%n1)**(-config%m1)
        zeta = effective_saturation**(1.0d0 / config%m1)
        call incomplete_beta%initialize(config%m1 + 1.0d0 / config%n1, 1.0d0 - 1.0d0 / config%n1)
        beta_value = 0.0d0
        converged = .false.
        call incomplete_beta%evaluate(zeta, beta_value, converged)
        expected_conductivity = config%k_sat * effective_saturation**config%l * beta_value**2

        call self%check_true("HCF VG incomplete beta converged through bound model", converged)
        call self%check_close("HCF VG polymorphic calc_Kflh dispatch", conductivity, expected_conductivity, 2.0d-14)
    end subroutine test_hcf_vg_dispatch

    subroutine test_hcf_thermal_conductivity(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        type(type_config_hcf) :: config
        type(type_state) :: state
        type(holder_hcfs) :: hcf
        type(type_iapws97), target :: water
        type(type_iapws06), target :: ice
        real(real64), parameter :: rho_g = 1000.0d0 * 9.80665d0
        real(real64), parameter :: gamma_0 = 71.88875d0
        real(real64) :: conductivity, thermal_conductivity, expected_thermal_conductivity
        real(real64) :: temperature, head, dgamma_dT

        call config%reset()
        config%model = HCF_MODES%BASE
        config%swcc_model = SWCC_MODELS%VG
        config%k_sat = 3.2d-6
        config%gain_factor = 1.0d0
        config%l = 0.5d0
        config%alpha1 = 1.11d0
        config%n1 = 1.48d0
        config%m1 = 0.2d0
        call water%initialize()
        call ice%initialize()
        call hcf%initialize(config, water, ice)

        temperature = -1.0d0
        head = -0.5d0
        call state%temperature%set(temperature)
        call state%pressure%set(rho_g * head)
        call state%water_content%set(0.3d0)
        call state%ice_content%set(0.0d0)

        conductivity = 0.0d0
        thermal_conductivity = 0.0d0
        call hcf%p%calc_Kflh(state, conductivity)
        call hcf%p%calc_KlT(state, thermal_conductivity)

        dgamma_dT = -0.1425d0 - 4.76d-4 * temperature
        expected_thermal_conductivity = conductivity * head * config%gain_factor * dgamma_dT / gamma_0
        call self%check_close("Hansson thermal liquid conductivity", thermal_conductivity, &
                              expected_thermal_conductivity, 2.0d-14)
        call self%check_true("thermal liquid conductivity is active", thermal_conductivity > 0.0d0)
    end subroutine test_hcf_thermal_conductivity

    subroutine test_freezing_storage_transport_split(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        type(type_config_wrf) :: wrf_config
        type(type_config_gcc) :: gcc_config
        type(holder_wrfs) :: wrf
        type(holder_gccs) :: gcc
        type(type_fusion) :: fusion
        type(type_iapws97), target :: water
        type(type_iapws06), target :: ice
        type(type_state) :: state
        real(real64) :: dwater_dP, dwater_dT, dwater_dT_cryo
        real(real64) :: cryogenic_suction, effective_suction, ice_content
        real(real64) :: theta_total, theta_liquid, theta_ice
        real(real64) :: d_effective_dP, d_effective_dT, d_cryo_dT, freezing_level_set
        integer(int32) :: i

        call self%configure_wrf(wrf_config, SWCC_MODELS%VG%ID)
        wrf_config%theta_s = 0.535d0
        wrf_config%theta_r = 0.05d0
        wrf_config%alpha1 = 1.11d0
        wrf_config%n1 = 1.48d0
        wrf_config%m1 = 0.2d0
        call wrf%initialize(wrf_config)

        call gcc_config%reset()
        gcc_config%gcc_model = GCC_TYPES%NON_SEGREGATION
        call water%initialize()
        call ice%initialize()
        call gcc%initialize(1_int32, gcc_config, water, ice)
        call fusion%initialize(wrf%p, gcc%p, water, ice)

        call state%porosity%set(0.535d0)

        ! The liquid follows whichever chemical-potential lowering binds, so the
        ! generalized suction is max(s_m, s_f), not their sum: a single mu_w is
        ! subject to two constraints and satisfies the stronger one. Their sum
        ! measures how far the state is from equilibrium, which is not a
        ! transport potential. Deep in the frozen zone (here s_f is two orders
        ! above s_m) the max selects s_f exactly.
        call state%temperature%set(-1.0d0)
        call state%pressure%set(-1.0d5)
        call gcc%p%calc(state, cryogenic_suction)
        call fusion%calc_effective_suction(state, effective_suction)
        call self%check_close("generalized suction selects the binding constraint", &
                              effective_suction, max(cryogenic_suction, 1.0d5), 1.0d-12)

        ! On that branch theta_l rides the freezing characteristic and no longer
        ! responds to the pore pressure. The pressure equation does not lose its
        ! diagonal to that: the storage it needs is the TOTAL water dTheta/dp =
        ! dtheta_SWRC/ds_m, which is untouched by freezing, and a pressure
        ! change there moves ice rather than liquid.
        call fusion%calc_water_content_derivatives(state, dwater_dP, dwater_dT_cryo)
        call self%check_close("frozen liquid content is pressure independent", dwater_dP, 0.0d0, 1.0d-14)
        call self%check_true("frozen liquid storage gains a freezing-curve slope", dwater_dT_cryo > 0.0d0)

        ! Above freezing psi_cryo vanishes and the relation reduces to ordinary
        ! unfrozen retention in p_w alone.
        call state%temperature%set(0.5d0)
        call state%pressure%set(-1.0d5)
        call fusion%calc_effective_suction(state, effective_suction)
        call self%check_close("unfrozen suction reduces to the capillary term", &
                              effective_suction, 1.0d5, 1.0d-12)
        call fusion%calc_water_content_derivatives(state, dwater_dP, dwater_dT)
        call self%check_true("unfrozen liquid storage retains its pore-pressure derivative", dwater_dP > 0.0d0)
        call self%check_close("unfrozen liquid storage has no freezing-curve slope", dwater_dT, 0.0d0, 1.0d-14)

        ! A warm, nearly saturated soil holds no ice. The generalized-Clapeyron
        ! suction is SIGNED above melting: ice is stable only where the liquid
        ! pressure is below p_i, and at +6.7 C that would take -8 MPa, so the
        ! max selects the matric branch by a margin of six orders of
        ! magnitude. Clipping s_f to zero there instead lets it tie or beat a
        ! small pore pressure, and the max returns a POSITIVE effective
        ! suction at a positive pore pressure: the retention curve of this soil then reports a liquid
        ! content below theta_s and the complement becomes ice. Measured on the
        ! clipped closure at this state: theta_i = 1.2e-2, in soil seven degrees
        ! above freezing.
        call state%temperature%set(6.7d0)
        call state%pressure%set(2.0d3)
        call gcc%p%calc(state, cryogenic_suction)
        call self%check_true("freezing suction is signed above melting", cryogenic_suction < -1.0d6)
        ! calc_effective_suction clamps the matric term at zero, so above
        ! saturation it publishes zero rather than the negative suction
        ! calc_phase_split works with; either way the freezing branch is not
        ! selected, which is what this measures.
        call fusion%calc_effective_suction(state, effective_suction)
        call self%check_close("warm saturated suction is not the freezing one", effective_suction, 0.0d0, 1.0d-12)
        call fusion%calc_ice_content(state, ice_content)
        call self%check_close("warm saturated soil carries no ice", ice_content, 0.0d0, 0.0d0)

        ! At the capillary/freezing crossing, the compact C1 maximum gives
        ! equal weights to both branches. This is the point where the exact
        ! maximum has a discontinuous tangent.
        call state%temperature%set(-0.1d0)
        do i = 1, 4
            call gcc%p%calc(state, cryogenic_suction)
            call state%pressure%set(-cryogenic_suction)
        end do
        call gcc%p%calc(state, cryogenic_suction)
        call gcc%p%deriv_temperature(state, d_cryo_dT)
        call fusion%calc_phase_split(state, theta_total, theta_liquid, theta_ice, &
                                     suction_effective_out=effective_suction, &
                                     dsuction_eff_dP_out=d_effective_dP, &
                                     dsuction_eff_dT_out=d_effective_dT)
        call self%check_true("C1 maximum rounds the suction value at the crossing", &
                             effective_suction > cryogenic_suction)
        call self%check_close("C1 maximum shares the pressure tangent at the crossing", &
                              d_effective_dP, -0.5d0, 1.0d-12)
        call self%check_close("C1 maximum shares the temperature tangent at the crossing", &
                              d_effective_dT, 0.5d0 * d_cryo_dT, 1.0d-12)
        call fusion%calc_freezing_level_set(state, freezing_level_set)
        call self%check_true("C1 level set places the exact crossing inside the ice band", &
                             freezing_level_set > 0.0d0)
    end subroutine test_freezing_storage_transport_split


    subroutine test_vapor_combined_evaluation(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        type(type_evaporation) :: evaporation
        type(type_iapws97), target :: water
        type(type_state) :: state
        real(real64) :: rh
        real(real64) :: vapor_separate, dP_separate, dT_separate
        real(real64) :: vapor_combined, dP_combined, dT_combined

        call water%initialize()
        call evaporation%initialize(water)
        call state%temperature%set(-1.0d0)
        call state%pressure%set(-10000.0d0)
        call state%air_content%set(0.2d0)
        call state%dQa_dP%set(-1.0d-7)
        call state%dQa_dT%set(1.0d-3)
        call evaporation%calc_relative_humidity(state, rh)
        call state%relative_humidity%set(rh)

        call evaporation%calc_vapor_content(state, vapor_separate)
        call evaporation%calc_vapor_content_derivatives(state, dP_separate, dT_separate)
        call evaporation%calc_vapor_content_with_derivatives(state, vapor_combined, dP_combined, dT_combined)

        call self%check_close("combined vapor content", vapor_combined, vapor_separate, 2.0d-14)
        call self%check_close("combined vapor pressure derivative", dP_combined, dP_separate, 2.0d-14)
        call self%check_close("combined vapor temperature derivative", dT_combined, dT_separate, 2.0d-14)
    end subroutine test_vapor_combined_evaluation

    !> The IAPWS memo answers what IAPWS answers, bit for bit, whatever the order of the states.
    subroutine test_iapws_memo(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        type(type_iapws_wrapper) :: wrapper
        type(type_iapws97), target :: water
        type(type_iapws06), target :: ice
        type(type_state) :: state
        ! A warm saturated state, a frozen suction state, then the first again.
        real(real64), parameter :: temperatures(3) = [5.0d0, -2.0d0, 5.0d0]
        real(real64), parameter :: pressures(3) = [2.0d4, -3.0d4, 2.0d4]
        real(real64), parameter :: T_K(3) = temperatures + 273.15d0
        real(real64), parameter :: P_abs(3) = [101325.0d0 + 2.0d4, 101325.0d0, 101325.0d0 + 2.0d4]
        real(real64) :: memo(9), direct(9), unused(6)
        integer(int32) :: k

        call water%initialize()
        call ice%initialize()
        call wrapper%initialize(water, ice)
        do k = 1, size(temperatures)
            call state%temperature%set(temperatures(k))
            call state%pressure%set(pressures(k))
            call wrapper%calc_rho_water(state, memo(1))
            call wrapper%calc_cp_water(state, memo(2))
            call wrapper%calc_drho_water_dT(state, memo(3))
            call wrapper%calc_rho_ice(state, memo(4))
            call wrapper%calc_cp_ice(state, memo(5))
            call wrapper%calc_rho_vapor_saturation(state, memo(6))
            call wrapper%calc_latent_heat_saturation(max(T_K(k), 273.16d0), memo(7))
            call wrapper%calc_dcp_water_dT(state, memo(8))
            call wrapper%calc_dcp_water_dP(state, memo(9))
            call water%calc_rho(T_K(k), P_abs(k), direct(1))
            call water%calc_cp(T_K(k), P_abs(k), direct(2))
            call water%calc_drho_dT(T_K(k), P_abs(k), direct(3))
            call ice%calc_rho(T_K(k), P_abs(k), direct(4))
            call ice%calc_cp(T_K(k), P_abs(k), direct(5))
            ! The memo fills every saturation quantity from one call with these outputs.
            call water%calc_saturation_properties(max(T_K(k), 273.16d0), rho_vap=direct(6), drho_dT_vap=unused(1), &
                                                  drho_dP_vap=unused(2), density_dT_vap=unused(3), cp_vap=unused(4), &
                                                  dcp_dT_vap=unused(5), latent_heat=direct(7), dlatent_heat_dT=unused(6))
            call water%calc_dcp_dT(T_K(k), P_abs(k), direct(8))
            call water%calc_dcp_dp(T_K(k), P_abs(k), direct(9))
            ! Suction holds the pressure at atmospheric, so the pressure derivative is zero there.
            if (pressures(k) < 0.0d0) direct(9) = 0.0d0
            call self%check_true("IAPWS memo returns the IAPWS values", all(memo == direct))
        end do
    end subroutine test_iapws_memo

    subroutine test_saturation_pressure(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        type(type_config_wrf) :: wrf_config
        type(type_config_gcc) :: gcc_config
        type(holder_wrfs) :: wrf
        type(holder_gccs) :: gcc
        type(type_phase_manager) :: phase_manager
        type(type_fusion) :: fusion_probe
        type(type_iapws97), target :: water
        type(type_iapws06), target :: ice
        type(type_state) :: state
        real(real64) :: saturation_pressure, Qw, Qi, Qa, porosity
        real(real64) :: dQw_dP, dQw_dT, dQi_dP, dQi_dT
        real(real64) :: dQw_dP_fd, dQw_dT_fd, dQi_dP_fd, dQi_dT_fd
        real(real64) :: Qw_fwd, Qw_bwd, Qi_fwd, Qi_bwd
        real(real64) :: t0, p0, fd_step_T, fd_step_P
        real(real64) :: projected_ice, ice_increment, equilibrium_error
        real(real64) :: theta_predicted, rho_w_probe, rho_i_probe
        logical :: is_saturated

        call self%configure_wrf(wrf_config, SWCC_MODELS%VG%ID)
        wrf_config%theta_s = 0.535d0
        wrf_config%theta_r = 0.05d0
        wrf_config%alpha1 = 1.11d0
        wrf_config%n1 = 1.48d0
        wrf_config%m1 = 0.2d0
        call wrf%initialize(wrf_config)

        call gcc_config%reset()
        gcc_config%gcc_model = GCC_TYPES%NON_SEGREGATION
        call water%initialize()
        call ice%initialize()
        call gcc%initialize(1_int32, gcc_config, water, ice)
        call phase_manager%initialize(gcc%p, wrf%p, water, ice)
        call fusion_probe%initialize(wrf%p, gcc%p, water, ice)

        porosity = 0.535d0
        call state%temperature%set(-1.0d0)
        call state%pressure%set(-5.4d4)
        call state%porosity%set(porosity)
        call state%ice_content%set(0.0d0)
        call phase_manager%update_water_phases(state)
        call state%water_content%get(Qw)
        call state%ice_content%get(Qi)
        call state%air_content%get(Qa)
        call state%dQw_dP%get(dQw_dP)
        call state%dQw_dT%get(dQw_dT)
        call state%dQi_dP%get(dQi_dP)
        call state%dQi_dT%get(dQi_dT)
        ! theta_i is a dependent state function of (T,p): a zeroed ice field is
        ! overwritten by the freezing curve, it is not carried forward.
        call self%check_true("cold state function produces ice from T-p alone", Qi > 0.0d0)
        call self%check_close("cold direct phase volumes close the pore volume", Qw + Qi + Qa, porosity, 2.0d-14)
        call state%ice_content%set(0.0d0)
        call phase_manager%update_water_phases(state)
        call state%ice_content%get(Qi_fwd)
        call self%check_close("ice content carries no history of its own", Qi_fwd, Qi, 1.0d-14)

        t0 = -1.0d0
        p0 = -5.4d4
        fd_step_T = 1.0d-5
        fd_step_P = 1.0d2

        call state%temperature%set(t0 + fd_step_T)
        call phase_manager%update_water_phases(state)
        call state%water_content%get(Qw_fwd)
        call state%ice_content%get(Qi_fwd)
        call state%temperature%set(t0 - fd_step_T)
        call phase_manager%update_water_phases(state)
        call state%water_content%get(Qw_bwd)
        call state%ice_content%get(Qi_bwd)
        call state%temperature%set(t0)
        dQw_dT_fd = (Qw_fwd - Qw_bwd) / (2.0d0 * fd_step_T)
        dQi_dT_fd = (Qi_fwd - Qi_bwd) / (2.0d0 * fd_step_T)
        call self%check_close("direct dQw/dT matches finite difference", dQw_dT, dQw_dT_fd, 1.0d-4)
        call self%check_close("state-function dQi/dT matches finite difference", dQi_dT, dQi_dT_fd, 1.0d-4)

        call state%pressure%set(p0 + fd_step_P)
        call phase_manager%update_water_phases(state)
        call state%water_content%get(Qw_fwd)
        call state%ice_content%get(Qi_fwd)
        call state%pressure%set(p0 - fd_step_P)
        call phase_manager%update_water_phases(state)
        call state%water_content%get(Qw_bwd)
        call state%ice_content%get(Qi_bwd)
        call state%pressure%set(p0)
        dQw_dP_fd = (Qw_fwd - Qw_bwd) / (2.0d0 * fd_step_P)
        dQi_dP_fd = (Qi_fwd - Qi_bwd) / (2.0d0 * fd_step_P)
        call self%check_close("direct dQw/dP matches finite difference", dQw_dP, dQw_dP_fd, 1.0d-4)
        call self%check_close("state-function dQi/dP matches finite difference", dQi_dP, dQi_dP_fd, 1.0d-4)

        call state%temperature%set(1.0d0)
        call state%pressure%set(-5.4d4)
        call phase_manager%update_water_phases(state)
        call state%ice_content%get(Qi)
        call self%check_close("warm direct T-p state is ice-free", Qi, 0.0d0, 1.0d-14)

        ! At zero gauge pressure s_m = 0, so Theta = theta_SWRC(s_m) = theta_s =
        ! porosity: the medium is saturated on the total/matric branch, so the
        ! saturation pressure is the pressure itself.
        call state%temperature%set(-1.0d0)
        call state%pressure%set(0.0d0)
        call phase_manager%update_water_phases(state)
        saturation_pressure = 1.0d0
        is_saturated = .false.
        call phase_manager%calc_saturation_pressure(state, saturation_pressure, is_saturated)
        call self%check_true("zero gauge pressure is saturated while frozen", is_saturated)
        call state%water_content%get(Qw)
        call state%ice_content%get(Qi)
        call state%air_content%get(Qa)
        call state%dQw_dP%get(dQw_dP)
        call state%dQi_dP%get(dQi_dP)
        call self%check_true("saturated state has zero gas volume", abs(Qa) < 1.0d-12)
        ! theta_i = (rho_w/rho_i)*(Theta - theta_w) is the complement of theta_w
        ! against Theta, not against the pore volume (see calc_phase_split,
        ! fusion.F90). Without a limit, ice being less dense than water would
        ! push theta_w+theta_i above the nominal pore volume here (the raw
        ! Clapeyron suction is far larger than the vanishing matric suction at
        ! this fully saturated state). calc_limited_cryo_suction is active at
        ! this state precisely because the unlimited split WOULD overfill the
        ! pore: it caps the cryogenic suction so theta_w+theta_i settles AT the
        ! pore volume instead of exceeding it (generalized Clapeyron with a
        ! stressed ice phase, see that routine's docstring), so the check below
        ! is now an equality up to the bisection's residual, not a volume-
        ! expansion overflow. air_content still saturates at zero either way.
        call self%check_true("ice is present in the saturated frozen state", Qi > 0.0d0)
        call self%check_true("saturated phase volumes exceed the pore volume when ice is present", &
                             Qw + Qi >= 0.535d0 - 2.0d-14)
        ! Pore-volume limiter identities (this state activates
        ! calc_limited_cryo_suction, see comment above): the active constraint
        ! r*Theta+(1-r)*theta_l=phi caps theta_w+theta_i AT the pore volume and
        ! makes the p-derivatives of theta_w and theta_i cancel, since neither
        ! phase can grow without the other shrinking once the pore is full.
        ! At psi_cap=0 the compact C1 maximum keeps a bounded residual because
        ! its transition band lies above the limiter's zero-suction floor.
        ! The residual is small relative to the unlimited overfill and vanishes
        ! outside the regularization band.
        call self%check_true("limited saturated state stays within the C1 pore tolerance", &
                             Qw + Qi <= porosity + 1.0d-2)
        call self%check_true("limited-branch volumetric derivatives cancel: d(Qw+Qi)/dP = 0", &
                             abs(dQw_dP + dQi_dP) <= 1.0d-8)
        ! The two closures make different promises here, so the assertion has to
        ! follow whichever one was compiled.
        call fusion_probe%calc_rho_water(state, rho_w_probe)
        call fusion_probe%calc_rho_ice(state, rho_i_probe)
        call wrf%p%calc(0.0d0, theta_predicted)
        if (PORE_LIMIT_EXPELS_WATER) then
            ! Expelling closure: theta_i is defined as phi - theta_l, so the
            ! pore volume is an exact equality, and the stored water
            ! alpha*phi + (1-alpha)*theta_l falls BELOW what the retention curve
            ! asks for. That deficit is the water the mass balance moves out of
            ! the node - closing against theta_SWRC here would mean the
            ! constraint was being met by suppressing the phase change again.
            call self%check_close("saturated frozen state fills the pore exactly", &
                                  Qw + Qi, porosity, 1.0d-12)
            call self%check_true("saturated frozen state stores less than the retention curve demands", &
                                 Qw + (rho_i_probe / rho_w_probe) * Qi < theta_predicted - 1.0d-6)
        else
            ! Suction-pinning closure: the split stays on the retention curve,
            ! and the pore bound is approached within the compact C1 band (see
            ! the note above).
            call self%check_close("saturated frozen state closes the water-conservation identity", &
                                  Qw + (rho_i_probe / rho_w_probe) * Qi, theta_predicted, 1.0d-10)
        end if

        ! Dropping the pressure far enough takes theta(psi_cap) below the bound,
        ! which deactivates it and lets a gas phase appear. The offset has to be
        ! on the retention curve's own pressure scale, not a fixed few Pa.
        call state%pressure%set(-5.4d4)
        call phase_manager%update_water_phases(state)
        call phase_manager%calc_saturation_pressure(state, saturation_pressure, is_saturated)
        call self%check_true("state well below saturation contains gas", .not. is_saturated)
        call state%water_content%get(Qw)
        call state%ice_content%get(Qi)
        call state%air_content%get(Qa)
        call self%check_true("unsaturated state has positive gas volume", Qa > 0.0d0)
        call self%check_close("unsaturated phase volumes close the pore volume", Qw + Qi + Qa, 0.535d0, 2.0d-14)
    end subroutine test_saturation_pressure

    !> Headline regression test for the smoothing defect: at the real case's
    !> warm, fully unfrozen initial condition (T=+6.7 degC, p=-54066 Pa),
    !> s_m exceeds s_f by more than the compact C1 band, so
    !> compute_effective_suction must return s_eff = s_m exactly. A hyperbolic
    !> smooth maximum has infinite tails and reported ice_content ~ 1e-3 by
    !> volume at every one of the real mesh's 2874 nodes at this same state.
    !> Compact support reaches the ice-free branch exactly outside the band,
    !> so ice_content and its P,T tangents must be bitwise zero here.
    subroutine test_no_spurious_ice_at_warm_state(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        type(type_config_wrf) :: wrf_config
        type(type_config_gcc) :: gcc_config
        type(holder_wrfs) :: wrf
        type(holder_gccs) :: gcc
        type(type_phase_manager) :: phase_manager
        type(type_iapws97), target :: water
        type(type_iapws06), target :: ice
        type(type_state) :: state
        real(real64) :: Qi, dQi_dP, dQi_dT

        call self%configure_wrf(wrf_config, SWCC_MODELS%VG%ID)
        wrf_config%theta_s = 0.535d0
        wrf_config%theta_r = 0.05d0
        wrf_config%alpha1 = 1.11d0
        wrf_config%n1 = 1.48d0
        wrf_config%m1 = 0.2d0
        call wrf%initialize(wrf_config)

        call gcc_config%reset()
        gcc_config%gcc_model = GCC_TYPES%NON_SEGREGATION
        call water%initialize()
        call ice%initialize()
        call gcc%initialize(1_int32, gcc_config, water, ice)
        call phase_manager%initialize(gcc%p, wrf%p, water, ice)

        call state%porosity%set(0.535d0)
        call state%ice_content%set(0.0d0)
        call state%temperature%set(6.7d0)
        call state%pressure%set(-54066.0d0)
        call phase_manager%update_water_phases(state)

        call state%ice_content%get(Qi)
        call state%dQi_dP%get(dQi_dP)
        call state%dQi_dT%get(dQi_dT)

        call self%check_close("warm initial condition is exactly ice-free", Qi, 0.0d0, 0.0d0)
        call self%check_close("warm initial condition has zero dQi/dP", dQi_dP, 0.0d0, 0.0d0)
        call self%check_close("warm initial condition has zero dQi/dT", dQi_dT, 0.0d0, 0.0d0)
    end subroutine test_no_spurious_ice_at_warm_state

    !> FD-vs-analytic acceptance sweep across the pore-volume limit
    !> (calc_limited_cryo_suction feeding
    !> compute_pore_limited_effective_suction, fusion.F90): each branch must
    !> report the tangent its own formula promises. This sweep spans p from
    !> deep-unfrozen (-30 kPa) up to 0 at T=-1 C, which crosses from
    !> is_limited=.false. (raw GCC tangents) into is_limited=.true. (the
    !> active pore-volume root's IFT tangent) partway through - the porosity
    !> gap (theta_s=0.52 < porosity=0.535) is chosen so that crossing happens
    !> at a genuine interior root of pore_excess, not pinned at the psi=0
    !> floor (see test_saturation_pressure's theta_s=porosity config, which
    !> deliberately hits that degenerate pinned case instead). The crossing
    !> itself is a kink: a central difference taken ACROSS it compares two
    !> different branches and is not expected to match either one's tangent,
    !> which is why the sweep samples the branches rather than the corner.
    !>
    !> Tolerance: relative error < 1e-4, OR absolute error below a small
    !> floor when the compared derivative itself is small. The floors are not
    !> reverse-fitted to hide a defect: they cover one specific, understood
    !> approximation - dE_dpsi/dE_dp in
    !> compute_pore_limited_effective_suction hold r=rho_w/rho_i fixed, while
    !> both rho_w (IAPWS-97, now evaluated below freezing too) and rho_i
    !> (IAPWS-06) carry a temperature dependence, so the active constraint's
    !> root psi* has a small residual T-sensitivity through r(T) that this
    !> formulation does not carry. It is largest deep in the limited branch,
    !> where dpsi_c_dT is exactly zero, which is where this sweep necessarily
    !> samples it. The floors keep headroom above the measured residual while
    !> a genuinely wrong branch tangent would miss by orders of magnitude
    !> more, which the relative-1e-4 leg still catches.
    !> p_l = -s_eff is the Darcy driver; p_w labels the total water. Guards
    !> the exported d(s_eff)/dX on both branches, and that cooling a closed
    !> point leaves Theta and p_w untouched while p_l follows Clapeyron.
    subroutine test_liquid_pressure_driver(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        type(type_config_wrf) :: wrf_config
        type(type_config_gcc) :: gcc_config
        type(holder_wrfs) :: wrf
        type(holder_gccs) :: gcc
        type(type_fusion) :: fusion
        type(type_iapws97), target :: water
        type(type_iapws06), target :: ice
        type(type_state) :: state
        ! Mizoguchi column initial condition: theta_SWRC(-p_w0) = 0.330.
        real(real64), parameter :: p_w0 = -5.40657412d4
        real(real64), parameter :: delta_p = 10.0d0
        real(real64), parameter :: delta_T = 1.0d-5
        real(real64) :: theta_total, theta_liquid, theta_ice
        real(real64) :: theta_total_warm, theta_liquid_warm, theta_ice_warm
        real(real64) :: s_eff, s_eff_warm, dSeff_dP, dSeff_dT, dtotal_dT
        real(real64) :: s_fwd, s_bwd, dSeff_dP_fd, dSeff_dT_fd
        real(real64) :: cryogenic_suction, rho_w, rho_i

        call self%configure_wrf(wrf_config, SWCC_MODELS%VG%ID)
        wrf_config%theta_s = 0.535d0
        wrf_config%theta_r = 0.05d0
        wrf_config%alpha1 = 1.11d0
        wrf_config%n1 = 1.48d0
        wrf_config%m1 = 0.2d0
        call wrf%initialize(wrf_config)

        call gcc_config%reset()
        gcc_config%gcc_model = GCC_TYPES%NON_SEGREGATION
        call water%initialize()
        call ice%initialize()
        call gcc%initialize(1_int32, gcc_config, water, ice)
        call fusion%initialize(wrf%p, gcc%p, water, ice)

        call state%porosity%set(0.535d0)
        call state%ice_content%set(0.0d0)

        ! --- Unlimited frozen branch: s_eff = s_f(T), so the driver loses its
        !     pressure sensitivity and keeps only the Clapeyron one. ---
        call state%temperature%set(-1.0d0)
        call state%pressure%set(-1.0d5)
        call split_here()
        call finite_difference_here()
        call self%check_close("frozen-branch d(s_eff)/dP matches finite difference", &
                              dSeff_dP, dSeff_dP_fd, 1.0d-6)
        call self%check_close("frozen-branch d(s_eff)/dT matches finite difference", &
                              dSeff_dT, dSeff_dT_fd, 1.0d-4)
        call self%check_close("frozen branch carries no pressure sensitivity", dSeff_dP, 0.0d0, 1.0d-14)
        call self%check_true("frozen branch carries the Clapeyron slope", dSeff_dT < -1.0d5)

        ! --- Pore-volume-limited state: the suction is NOT pinned any more.
        !     The constraint is carried by the stored water instead, so the
        !     Clapeyron slope must survive here - suppressing it is exactly the
        !     apparent-capacity collapse the expelling closure removes. ---
        call state%temperature%set(-0.5d0)
        call state%pressure%set(-1.0d3)
        call split_here()
        call finite_difference_here()
        call self%check_close("pore-limited d(s_eff)/dP matches finite difference", &
                              dSeff_dP, dSeff_dP_fd, 1.0d-4)
        if (PORE_LIMIT_EXPELS_WATER) then
            ! Nothing pins the suction, so the Clapeyron slope survives - that
            ! survival is the whole point of the expelling closure.
            call self%check_close("pore-limited d(s_eff)/dT matches finite difference", &
                                  dSeff_dT, dSeff_dT_fd, 1.0d-4)
            call self%check_true("pore-volume limit keeps the Clapeyron coupling", dSeff_dT < -1.0d5)
        else
            ! The pinned branch drops dpsi_star_dT (~10 Pa/K) and suppresses the
            ! Clapeyron coupling, which is what collapses the apparent capacity.
            call self%check_true("limited-branch d(s_eff)/dT omits only the density-ratio term", &
                                 abs(dSeff_dT - dSeff_dT_fd) <= 5.0d1)
            call self%check_true("pore-volume limit suppresses the Clapeyron coupling", &
                                 abs(dSeff_dT_fd) < 1.0d3)
            call self%check_true("pore-volume limit revives the pressure sensitivity", abs(dSeff_dP) > 1.0d-6)
        end if

        ! --- Closed material point, cooled from the melting point. ---
        call state%temperature%set(0.0d0)
        call state%pressure%set(p_w0)
        call fusion%calc_phase_split(state, theta_total_warm, theta_liquid_warm, theta_ice_warm, &
                                     suction_effective_out=s_eff_warm)
        call self%check_close("melting point is ice-free", theta_ice_warm, 0.0d0, 1.0d-14)
        call self%check_close("melting-point driver reduces to the pore pressure", &
                              -s_eff_warm, p_w0, 1.0d-12)

        call state%temperature%set(-0.1d0)
        call split_here()

        call self%check_close("cooling a closed point leaves the total water invariant", &
                              theta_total, theta_total_warm, 1.0d-12)
        ! The pore cap is inactive at this state under either closure, so the
        ! retention total's zero temperature tangent must survive the expelling
        ! branch's longer cancellation too. Absolute, not relative: check_close
        ! scales its tolerance by max(1e-14, |values|), which against an exact
        ! zero demands 1e-28 and so cannot express "zero to round-off".
        call self%check_true("total water carries no temperature tangent", abs(dtotal_dT) <= 1.0d-14)
        call self%check_true("cooling a closed point produces ice", theta_ice > 0.0d0)
        call fusion%calc_rho_water(state, rho_w)
        call fusion%calc_rho_ice(state, rho_i)
        call self%check_close("closed point closes the water-conservation identity", &
                              theta_liquid + (rho_i / rho_w) * theta_ice, theta_total, 1.0d-12)

        ! The Clapeyron depression lands on p_l, not on p_w.
        call gcc%p%calc(state, cryogenic_suction)
        call self%check_close("Clapeyron suction at -0.1 C", cryogenic_suction, 1.2249d5, 5.0d-3)
        call self%check_close("driver selects the binding constraint", &
                              s_eff, max(cryogenic_suction, -p_w0), 1.0d-12)
        call self%check_true("liquid pressure is far below the pore pressure", &
                             (-s_eff) - p_w0 < -5.0d4)

    contains

        subroutine split_here()
            implicit none

            theta_total = 0.0d0
            theta_liquid = 0.0d0
            theta_ice = 0.0d0
            dtotal_dT = 0.0d0
            dSeff_dP = 0.0d0
            dSeff_dT = 0.0d0
            call fusion%calc_phase_split(state, theta_total, theta_liquid, theta_ice, &
                                         dtotal_dT=dtotal_dT, &
                                         suction_effective_out=s_eff, &
                                         dsuction_eff_dP_out=dSeff_dP, &
                                         dsuction_eff_dT_out=dSeff_dT)
        end subroutine split_here

        !> Central differences of the published effective suction around the
        !> state currently held in `state`, which is restored on exit.
        subroutine finite_difference_here()
            implicit none

            real(real64) :: pressure, temperature

            call state%pressure%get(pressure)
            call state%temperature%get(temperature)

            call state%pressure%set(pressure + delta_p)
            call fusion%calc_effective_suction(state, s_fwd)
            call state%pressure%set(pressure - delta_p)
            call fusion%calc_effective_suction(state, s_bwd)
            call state%pressure%set(pressure)
            dSeff_dP_fd = (s_fwd - s_bwd) / (2.0d0 * delta_p)

            call state%temperature%set(temperature + delta_T)
            call fusion%calc_effective_suction(state, s_fwd)
            call state%temperature%set(temperature - delta_T)
            call fusion%calc_effective_suction(state, s_bwd)
            call state%temperature%set(temperature)
            dSeff_dT_fd = (s_fwd - s_bwd) / (2.0d0 * delta_T)
        end subroutine finite_difference_here
    end subroutine test_liquid_pressure_driver

    subroutine test_blend_derivative_sweep(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        type(type_config_wrf) :: wrf_config
        type(type_config_gcc) :: gcc_config
        type(holder_wrfs) :: wrf
        type(holder_gccs) :: gcc
        type(type_phase_manager) :: phase_manager
        type(type_iapws97), target :: water
        type(type_iapws06), target :: ice
        type(type_state) :: state
        real(real64), parameter :: pressures(7) = &
            [-3.0d4, -2.0d4, -1.2d4, -8.0d3, -4.0d3, -1.0d3, 0.0d0]
        real(real64), parameter :: delta_p = 10.0d0
        real(real64), parameter :: delta_T = 1.0d-5
        real(real64), parameter :: rel_tol = 1.0d-4
        real(real64), parameter :: dP_abs_floor = 1.0d-5
        real(real64), parameter :: dT_abs_floor = 1.0d-3
        real(real64) :: porosity, t0
        real(real64) :: Qw, dQw_dP, dQw_dT
        real(real64) :: Qw_fwd, Qw_bwd, dQw_dP_fd, dQw_dT_fd
        integer(int32) :: i, fail_count
        character(len=96) :: label

        call self%configure_wrf(wrf_config, SWCC_MODELS%VG%ID)
        wrf_config%theta_s = 0.52d0
        wrf_config%theta_r = 0.05d0
        wrf_config%alpha1 = 1.11d0
        wrf_config%n1 = 1.48d0
        wrf_config%m1 = 0.2d0
        call wrf%initialize(wrf_config)

        call gcc_config%reset()
        gcc_config%gcc_model = GCC_TYPES%NON_SEGREGATION
        call water%initialize()
        call ice%initialize()
        call gcc%initialize(1_int32, gcc_config, water, ice)
        call phase_manager%initialize(gcc%p, wrf%p, water, ice)

        porosity = 0.535d0
        t0 = -1.0d0
        call state%porosity%set(porosity)
        call state%ice_content%set(0.0d0)

        fail_count = 0
        do i = 1, size(pressures)
            call state%temperature%set(t0)
            call state%pressure%set(pressures(i))
            call phase_manager%update_water_phases(state)
            call state%water_content%get(Qw)
            call state%dQw_dP%get(dQw_dP)
            call state%dQw_dT%get(dQw_dT)

            call state%pressure%set(pressures(i) + delta_p)
            call phase_manager%update_water_phases(state)
            call state%water_content%get(Qw_fwd)
            call state%pressure%set(pressures(i) - delta_p)
            call phase_manager%update_water_phases(state)
            call state%water_content%get(Qw_bwd)
            call state%pressure%set(pressures(i))
            dQw_dP_fd = (Qw_fwd - Qw_bwd) / (2.0d0 * delta_p)

            call state%temperature%set(t0 + delta_T)
            call phase_manager%update_water_phases(state)
            call state%water_content%get(Qw_fwd)
            call state%temperature%set(t0 - delta_T)
            call phase_manager%update_water_phases(state)
            call state%water_content%get(Qw_bwd)
            call state%temperature%set(t0)
            dQw_dT_fd = (Qw_fwd - Qw_bwd) / (2.0d0 * delta_T)

            write (label, '(A,ES10.3,A)') "blend sweep dQw/dP at p=", pressures(i), " Pa"
            if (.not. within_tolerance(dQw_dP, dQw_dP_fd, dP_abs_floor)) fail_count = fail_count + 1
            call self%check_true(trim(label), within_tolerance(dQw_dP, dQw_dP_fd, dP_abs_floor))

            write (label, '(A,ES10.3,A)') "blend sweep dQw/dT at p=", pressures(i), " Pa"
            if (.not. within_tolerance(dQw_dT, dQw_dT_fd, dT_abs_floor)) fail_count = fail_count + 1
            call self%check_true(trim(label), within_tolerance(dQw_dT, dQw_dT_fd, dT_abs_floor))
        end do

        if (fail_count > 0) write (error_unit, '(A,I0,A)') "blend derivative sweep: ", fail_count, " point(s) failed"

    contains

        !> Relative error < rel_tol, OR absolute error < abs_floor: the second
        !> leg only matters where the compared derivative itself is small
        !> (see the routine's docstring for which residual it covers here).
        function within_tolerance(analytic, fd, abs_floor) result(ok)
            implicit none
            real(real64), intent(in) :: analytic, fd, abs_floor
            logical :: ok
            real(real64) :: diff, scale

            diff = abs(analytic - fd)
            scale = max(abs(analytic), abs(fd), tiny(1.0d0))
            ok = (diff <= rel_tol * scale) .or. (diff <= abs_floor)
        end function within_tolerance
    end subroutine test_blend_derivative_sweep

    subroutine benchmark_hot_paths(self)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self

        integer(int32), parameter :: beta_evaluations = 2000000
        integer(int32), parameter :: hcf_evaluations = 1000000
        integer(int32), parameter :: capacity_gets = 10000000
        integer(int32), parameter :: capacity_scans = 100000
        integer(int32), parameter :: wrf_initializations = 2000
        integer(int32), parameter :: vapor_evaluations = 20000
        real(real64), parameter :: rho_g = 1000.0d0 * 9.80665d0
        type(type_regularized_incomplete_beta) :: incomplete_beta
        type(type_config_wrf) :: wrf_config
        type(type_config_hcf) :: hcf_config
        type(holder_wrfs) :: wrf
        type(holder_hcfs) :: hcf_beta, hcf_closed
        type(type_state) :: state
        type(type_iapws97), target :: water
        type(type_iapws06), target :: ice
        type(type_evaporation) :: evaporation
        integer(int64) :: clock_start, clock_end, clock_rate
        integer(int32) :: i, j
        real(real64) :: value, x, pressure, checksum, elapsed_seconds
        real(real64) :: h_abs_min, h_abs_max, h_abs, fraction
        real(real64) :: derivative, maximum_derivative, x_peak
        real(real64) :: rh, dvapor_dP, dvapor_dT
        logical :: converged

        call system_clock(count_rate=clock_rate)
        call water%initialize()

        call incomplete_beta%initialize(0.2d0 + 1.0d0 / 1.48d0, 1.0d0 - 1.0d0 / 1.48d0)
        checksum = 0.0d0
        call system_clock(clock_start)
        do i = 1, beta_evaluations
            x = 0.01d0 + 0.98d0 * real(mod(i, 997), real64) / 996.0d0
            value = 0.0d0
            converged = .false.
            call incomplete_beta%evaluate(x, value, converged)
            if (converged) checksum = checksum + value
        end do
        call system_clock(clock_end)
        elapsed_seconds = real(clock_end - clock_start, real64) / real(clock_rate, real64)
        call self%report_benchmark("incomplete beta evaluate", elapsed_seconds, beta_evaluations, checksum)

        call self%configure_wrf(wrf_config, SWCC_MODELS%VG%ID)
        call wrf%initialize(wrf_config)
        checksum = 0.0d0
        call system_clock(clock_start)
        do i = 1, capacity_gets
            value = 0.0d0
            call wrf%p%calc_lscheme_capacity(value)
            checksum = checksum + value
        end do
        call system_clock(clock_end)
        elapsed_seconds = real(clock_end - clock_start, real64) / real(clock_rate, real64)
        call self%report_benchmark("cached capacity get", elapsed_seconds, capacity_gets, checksum)

        h_abs_min = min(abs(wrf_config%h_crit), 1.0d0 / wrf_config%alpha1, 1.0d0 / wrf_config%alpha2) * 1.0d-4
        h_abs_max = max(abs(wrf_config%h_crit), 1.0d0 / wrf_config%alpha1, 1.0d0 / wrf_config%alpha2) * 1.0d4
        checksum = 0.0d0
        call system_clock(clock_start)
        do i = 1, capacity_scans
            maximum_derivative = 0.0d0
            do j = 0, 96
                fraction = real(j, real64) / 96.0d0
                h_abs = exp(log(h_abs_min) + fraction * (log(h_abs_max) - log(h_abs_min)))
                derivative = 0.0d0
                call wrf%p%deriv(-h_abs, derivative)
                maximum_derivative = max(maximum_derivative, derivative)
            end do
            x_peak = ((wrf_config%n1 - 1.0d0) / &
                      (wrf_config%m1 * wrf_config%n1 + 1.0d0))**(1.0d0 / wrf_config%n1)
            derivative = 0.0d0
            call wrf%p%deriv(-x_peak / wrf_config%alpha1, derivative)
            maximum_derivative = max(maximum_derivative, derivative)
            x_peak = ((wrf_config%n2 - 1.0d0) / &
                      (wrf_config%m2 * wrf_config%n2 + 1.0d0))**(1.0d0 / wrf_config%n2)
            derivative = 0.0d0
            call wrf%p%deriv(-x_peak / wrf_config%alpha2, derivative)
            maximum_derivative = max(maximum_derivative, derivative)
            derivative = 0.0d0
            call wrf%p%deriv(wrf_config%h_crit * (1.0d0 + 1.0d-8), derivative)
            maximum_derivative = max(maximum_derivative, derivative)
            checksum = checksum + maximum_derivative / rho_g
        end do
        call system_clock(clock_end)
        elapsed_seconds = real(clock_end - clock_start, real64) / real(clock_rate, real64)
        call self%report_benchmark("legacy-equivalent capacity scan", elapsed_seconds, capacity_scans, checksum)

        checksum = 0.0d0
        call system_clock(clock_start)
        do i = 1, wrf_initializations
            call wrf%initialize(wrf_config)
            value = 0.0d0
            call wrf%p%calc_lscheme_capacity(value)
            checksum = checksum + value
        end do
        call system_clock(clock_end)
        elapsed_seconds = real(clock_end - clock_start, real64) / real(clock_rate, real64)
        call self%report_benchmark("WRF initialize and capacity scan", elapsed_seconds, wrf_initializations, checksum)

        call evaporation%initialize(water)
        call state%temperature%set(-1.0d0)
        call state%pressure%set(-10000.0d0)
        call state%air_content%set(0.2d0)
        call state%dQa_dP%set(-1.0d-7)
        call state%dQa_dT%set(1.0d-3)
        call evaporation%calc_relative_humidity(state, rh)
        call state%relative_humidity%set(rh)

        checksum = 0.0d0
        call system_clock(clock_start)
        do i = 1, vapor_evaluations
            value = 0.0d0
            dvapor_dP = 0.0d0
            dvapor_dT = 0.0d0
            call evaporation%calc_vapor_content(state, value)
            call evaporation%calc_vapor_content_derivatives(state, dvapor_dP, dvapor_dT)
            checksum = checksum + value + dvapor_dP + dvapor_dT
        end do
        call system_clock(clock_end)
        elapsed_seconds = real(clock_end - clock_start, real64) / real(clock_rate, real64)
        call self%report_benchmark("separate vapor value and derivatives", elapsed_seconds, vapor_evaluations, checksum)

        checksum = 0.0d0
        call system_clock(clock_start)
        do i = 1, vapor_evaluations
            value = 0.0d0
            dvapor_dP = 0.0d0
            dvapor_dT = 0.0d0
            call evaporation%calc_vapor_content_with_derivatives(state, value, dvapor_dP, dvapor_dT)
            checksum = checksum + value + dvapor_dP + dvapor_dT
        end do
        call system_clock(clock_end)
        elapsed_seconds = real(clock_end - clock_start, real64) / real(clock_rate, real64)
        call self%report_benchmark("combined vapor value and derivatives", elapsed_seconds, vapor_evaluations, checksum)

        call hcf_config%reset()
        hcf_config%model = HCF_MODES%BASE
        hcf_config%swcc_model = SWCC_MODELS%VG
        hcf_config%k_sat = 1.0d0
        hcf_config%l = 0.5d0
        hcf_config%alpha1 = 2.0d0
        hcf_config%n1 = 1.48d0
        hcf_config%m1 = 0.2d0
        call hcf_beta%initialize(hcf_config, water, ice)

        hcf_config%n1 = 2.0d0
        hcf_config%m1 = 0.5d0
        call hcf_closed%initialize(hcf_config, water, ice)

        pressure = -rho_g
        call state%pressure%set(pressure)
        checksum = 0.0d0
        call system_clock(clock_start)
        do i = 1, hcf_evaluations
            value = 0.0d0
            call hcf_closed%p%calc_Kflh(state, value)
            checksum = checksum + value
        end do
        call system_clock(clock_end)
        elapsed_seconds = real(clock_end - clock_start, real64) / real(clock_rate, real64)
        call self%report_benchmark("HCF VG closed form", elapsed_seconds, hcf_evaluations, checksum)

        checksum = 0.0d0
        call system_clock(clock_start)
        do i = 1, hcf_evaluations
            value = 0.0d0
            call hcf_beta%p%calc_Kflh(state, value)
            checksum = checksum + value
        end do
        call system_clock(clock_end)
        elapsed_seconds = real(clock_end - clock_start, real64) / real(clock_rate, real64)
        call self%report_benchmark("HCF VG incomplete beta", elapsed_seconds, hcf_evaluations, checksum)
    end subroutine benchmark_hot_paths

    subroutine report_benchmark(self, name, elapsed_seconds, evaluations, checksum)
        implicit none
        class(type_constitutive_test_suite), intent(in) :: self
        character(len=*), intent(in) :: name
        real(real64), intent(in) :: elapsed_seconds
        integer(int32), intent(in) :: evaluations
        real(real64), intent(in) :: checksum

        real(real64) :: nanoseconds_per_evaluation

        nanoseconds_per_evaluation = elapsed_seconds * 1.0d9 / real(evaluations, real64)
        write (output_unit, '(A,A,A,F10.6,A,F12.3,A,ES16.8)') &
            "BENCH: ", trim(name), " total[s]=", elapsed_seconds, &
            " ns/eval=", nanoseconds_per_evaluation, " checksum=", checksum
    end subroutine report_benchmark

    subroutine configure_wrf(self, config, model_id)
        implicit none
        class(type_constitutive_test_suite), intent(in) :: self
        type(type_config_wrf), intent(inout) :: config
        integer(int32), intent(in) :: model_id

        call config%reset()
        config%material_id = 1
        config%swcc_model = SWCC_MODELS%to_object(model_id)
        config%theta_s = 0.45d0
        config%theta_r = 0.05d0
        config%alpha1 = 2.0d0
        config%n1 = 2.0d0
        config%m1 = 0.5d0
        config%h_crit = -0.05d0
        config%alpha2 = 0.2d0
        config%n2 = 3.0d0
        config%m2 = 2.0d0 / 3.0d0
        config%w1 = 0.4d0
        config%w2 = 0.6d0

        if (model_id == SWCC_MODELS%BC%ID .or. model_id == SWCC_MODELS%KO%ID) then
            config%alpha1 = -0.5d0
            config%n1 = 2.0d0
        end if
    end subroutine configure_wrf

    subroutine check_close(self, name, actual, expected, relative_tolerance)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(real64), intent(in) :: actual
        real(real64), intent(in) :: expected
        real(real64), intent(in) :: relative_tolerance

        real(real64) :: scale

        scale = max(1.0d-14, abs(actual), abs(expected))
        if (abs(actual - expected) <= relative_tolerance * scale) then
            write (output_unit, '(A)') "PASS: "//trim(name)
        else
            self%failures = self%failures + 1
            write (error_unit, '(A)') "FAIL: "//trim(name)
            write (error_unit, '(A,ES24.16)') "  actual   = ", actual
            write (error_unit, '(A,ES24.16)') "  expected = ", expected
        end if
    end subroutine check_close

    subroutine check_true(self, name, condition)
        implicit none
        class(type_constitutive_test_suite), intent(inout) :: self
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition

        if (condition) then
            write (output_unit, '(A)') "PASS: "//trim(name)
        else
            self%failures = self%failures + 1
            write (error_unit, '(A)') "FAIL: "//trim(name)
        end if
    end subroutine check_true

end module test_constitutive_suite

program test_constitutive
    use :: test_constitutive_suite, only: type_constitutive_test_suite
    implicit none

    type(type_constitutive_test_suite) :: suite

    call suite%run()
end program test_constitutive
