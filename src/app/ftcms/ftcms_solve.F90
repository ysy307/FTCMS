submodule(app_ftcms) ftcms_solve
    use :: core_parallel_reduce, only:reduce_sum, reduce_any
    use :: models_phase_change_chemical_potential, only:calc_T_high_celsius
    use :: control_homotopy_manager, only:HOMOTOPY_STAGE_ITERATIONS
    implicit none

    !> Smallest accepted line-search scale that still counts as a healthy
    !> homotopy stage. Below it the increment was too large to be a usable
    !> initial guess for the next one.
    real(real64), parameter :: HOMOTOPY_STAGE_MIN_SCALE = 2.5d-1
    !> A stage that no longer reduces the merit by this factor has nothing left
    !> to contribute to the initial guess, so its iteration budget is dropped.
    real(real64), parameter :: HOMOTOPY_STAGE_PROGRESS = 9.0d-1
    !> Coldest nodal temperature that still leaves the freezing suction inactive
    !> everywhere [C]. Above it the blend is the identity for every lambda, so
    !> the whole ladder is redundant work on one unchanging problem.
    real(real64), parameter :: HOMOTOPY_ONSET_MARGIN = 1.0d-3
    !> One line per time step reporting how far the ladder got.
    logical, parameter :: HOMOTOPY_VERBOSE = .true.

    !> Report which nodes dominate the local error at a failed nonlinear solve.
    !> Diagnostic only: it says whether a stall is a few bad nodes or a
    !> uniformly slow field, which the WRMS gate alone cannot distinguish.
    logical, parameter :: ERROR_NODE_REPORT = .true.
    integer(int32), parameter :: ERROR_NODE_REPORT_LIMIT = 6
    integer(int32), save :: error_node_report_count = 0

    integer(int32), parameter :: SOLVE_STATUS_NOT_RUN = 0
    integer(int32), parameter :: SOLVE_STATUS_CONVERGED = 1
    integer(int32), parameter :: SOLVE_STATUS_LINEAR_FAILURE = 2
    integer(int32), parameter :: SOLVE_STATUS_NONLINEAR_DIVERGED = 3
    integer(int32), parameter :: SOLVE_STATUS_NONLINEAR_LIMIT = 4
    integer(int32), parameter :: SOLVE_STATUS_PHASE_FAILURE = 5
    integer(int32), parameter :: SOLVE_STATUS_LTE_REJECTED = 6

    !> Audit the assembled tangent against the element finite-difference tangent
    !> once, at the first iterate where the line search finds no usable step.
    !> Diagnostic only: it costs O(4 n_node) element residual evaluations per
    !> element and must not be left on for production runs.
    logical, parameter :: FD_JACOBIAN_AUDIT = .false.
    logical, save :: fd_audit_done = .false.

    !> One-shot scan of the residual along the Newton direction at a genuinely
    !> stalled iterate. Diagnostic only.
    logical, parameter :: RESIDUAL_SCAN = .false.
    integer(int32), parameter :: RESIDUAL_SCAN_ITER = 12
    logical, save :: residual_scan_done = .false.
    logical, save :: branch_report_done = .false.

    ! --- Backtracking line search ---
    ! The measured descent region lies below alpha ~ 1e-5; backtracking by
    ! 1/4 reaches 2.4e-7 within MAX_TRIALS. At module scope because both the
    ! nonlinear loop and the homotopy stage loop run the same search.
    integer(int32), parameter :: LINE_SEARCH_MAX_TRIALS = 12
    real(real64), parameter :: LINE_SEARCH_BACKTRACK = 2.5d-1
    real(real64), parameter :: LINE_SEARCH_MIN_SCALE = 1.0d-6
    !> Armijo sufficient-decrease coefficient on the residual norm.
    real(real64), parameter :: LINE_SEARCH_ARMIJO = 1.0d-4

    !> Freezing-band telemetry at the iterations that exhaust the budget.
    integer(int32), parameter :: FREEZE_REPORT_ITER = 10
    integer(int32), parameter :: FREEZE_REPORT_LIMIT = 40
    integer(int32), save :: freeze_report_count = 0
    real(real64), allocatable, save :: freeze_prev_liquid(:), freeze_prev_ice(:)

    !> Predictor used to start each nonlinear iteration.
    !>
    !>   0  previous step (zeroth order)
    !>   1  linear extrapolation of T and p
    !>   2  linear extrapolation clipped at the freezing point
    !>   3  linear extrapolation of the conserved quantities, inverted for T
    !>
    !> Mode 1 is kept only as the measured counter-example: extrapolating T
    !> linearly through the phase boundary drove the second step to p = 4.2e5 Pa
    !> with a pressure increment of -3.5e9 Pa. T(t) has a kink where a node
    !> freezes, so a straight line through it lands far past the transition.
    integer(int32), parameter :: PREDICTOR_MODE = 0

contains
    module subroutine solve_time_step_initial_setup_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        real(real64), pointer, contiguous, dimension(:) :: u

        nullify (u)

        ! Reset iteration control
        ! reset() may set the compute solver to NONE when config is NONE.
        call self%control%reset_iteration()

        ! Anderson(1) history is only meaningful within one nonlinear loop.
        self%aa_has_prev = .false.
        self%aa_gnorm_prev = -1.0d0
        self%last_line_search_failures = 0
        self%last_line_search_trials = 0
        self%last_line_search_scale = 1.0d0
        if (allocated(self%phase_onset_reset)) self%phase_onset_reset = .false.

        ! [Important] Compute solver must always be PICARD or NEWTON if not NONE.
        ! Even for linear config (where iter=1 is forced), Picard discretization
        ! is often the base, but if explicitly NONE, we should respect it.
        if (.not. self%control%is_none()) then
            call self%control%set_nonlinear_solver(NONLINEAR_SOLVER%PICARD)
        end if

        call self%control%increment_total()
        call self%control%reset_acceleration()

        ! Supply the discretization scales the residual floor is derived from.
        ! The mesh measure is fixed, so it is summed once and cached; dt changes
        ! every attempt, so the floor is refreshed here rather than at setup.
        block
            real(real64) :: dt_now, measure
            real(real64) :: atol_temperature_u, atol_pressure_u, rtol_u
            integer(int32) :: i_elem, num_elements
            logical :: first_setup

            first_setup = (self%domain_measure_total <= 0.0d0)
            if (first_setup) then
                call self%domain%get_num_fe(num_elements)
                self%domain_measure_total = 0.0d0
                do i_elem = 1, num_elements
                    measure = 0.0d0
                    call self%domain%calc_measure(i_elem, measure)
                    self%domain_measure_total = self%domain_measure_total + abs(measure)
                end do
                ! Each rank measured only its own cells.
                self%domain_measure_total = reduce_sum(self%domain_measure_total)
                call build_nodal_volume(self, num_elements)
            end if
            call self%control%get_dt(dt_now)
            if (first_setup) then
                ! ATS error-control tolerances (Conditions.json:
                ! adaptive_stepping/error_control) are static run
                ! configuration, so they are read once here and cached inside
                ! convergence_control, reused as the primary-variable
                ! tolerance scale of the local error measure (see
                ! local_error_block in convergence_control.F90).
                call self%control%get_error_control_tolerances(atol_temperature_u, atol_pressure_u, rtol_u)
                call self%control%set_residual_scale(self%domain_measure_total, dt_now, self%nodal_volume, &
                                                     atol_temperature_u=atol_temperature_u, &
                                                     atol_pressure_u=atol_pressure_u, rtol_u=rtol_u)
            else
                call self%control%set_residual_scale(self%domain_measure_total, dt_now, self%nodal_volume)
            end if
        end block

        ! Save previous step values (Previous <- Current)
        call self%porosity%get_previous(u)
        if (associated(u)) then
            call self%porosity%set_current(u)
            nullify (u)
        end if

        if (self%control%is_physics_active(PHYSICS_TYPES%THERMAL)) then
            call self%temperature%get_previous(u)
            if (associated(u)) then
                call self%temperature%set_current(u)
                nullify (u)
            end if
        end if

        if (self%control%is_physics_active(PHYSICS_TYPES%HYDRAULIC)) then
            call self%pressure%get_previous(u)
            if (associated(u)) then
                call self%pressure%set_current(u)
                nullify (u)
            end if
        end if

        ! Restore derived phase fields too. A failed nonlinear attempt may leave
        ! Qw/Qi/Qa/Qv at the rejected iterate; retrying a step must start from the
        ! last accepted thermodynamic state, just like T and p.
        call self%Qw%get_previous(u)
        if (associated(u)) then
            call self%Qw%set_current(u)
            nullify (u)
        end if

        call self%Qi%get_previous(u)
        if (associated(u)) then
            call self%Qi%set_current(u)
            nullify (u)
        end if

        call self%Qa%get_previous(u)
        if (associated(u)) then
            call self%Qa%set_current(u)
            nullify (u)
        end if

        call self%Qv%get_previous(u)
        if (associated(u)) then
            call self%Qv%set_current(u)
            nullify (u)
        end if

    end subroutine solve_time_step_initial_setup_ftcms

    !> Extrapolate the accepted solution history to start the nonlinear
    !> iteration from a predicted state rather than from the previous step.
    !>
    !> The iteration otherwise starts at y_n, which is a zeroth-order predictor:
    !> the whole change over the step has to be produced by the Newton
    !> corrections. At the freezing onset that means one correction has to carry
    !> several hundred nodes across the phase boundary at once, which is exactly
    !> the step that overshoots. A linear predictor
    !>   y^(0) = y_n + (dt_{n+1}/dt_n) (y_n - y_{n-1})
    !> starts with the front already advanced, so the corrections stay small.
    !> After advance() the history holds y_n in slots 1 and 2 and y_{n-1} in
    !> slot 3, which is the convention apply_phase_predictor also uses.
    !>
    !> This changes the starting point only; the converged state is unchanged.
    subroutine apply_solution_predictor(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        !> Never extrapolate further than the step is allowed to grow.
        real(real64), parameter :: PREDICTOR_MAX_RATIO = 2.0d0
        real(real64) :: dt_new, ratio, history(3)
        real(real64), allocatable :: pressure_pred(:)
        real(real64), pointer, contiguous :: current(:)
        real(real64), pointer, contiguous :: bdf_coeffs(:)
        integer(int32) :: node_id, num_nodes, bdf_order

        ! Two genuine accepted steps are required before anything can be
        ! extrapolated. With only one, the oldest history slot still holds the
        ! initial condition, and y_1 - y_0 is the response to switching the
        ! boundary conditions on - a startup transient, not a rate. Predicting
        ! from it doubles that transient: measured, it drove the second step to
        ! p = 4.2e5 Pa. lte_has_second_difference is exactly the flag that says
        ! two steps have been committed.
        if (.not. self%lte_has_second_difference) return
        if (self%last_accepted_dt <= tiny(1.0d0)) return
        call self%control%get_dt(dt_new)
        if (dt_new <= 0.0d0) return
        ratio = min(dt_new / self%last_accepted_dt, PREDICTOR_MAX_RATIO)

        call self%domain%get_num_nodes(num_nodes)
        nullify (current, bdf_coeffs)

        ! Pressure is predicted linearly in every mode: p carries no phase-change
        ! kink of its own, the cryosuction kink lives in psi_cryo(T).
        allocate (pressure_pred(num_nodes), source=0.0d0)
        if (self%is_active_hydraulic()) then
            call self%pressure%get_current(current)
            if (associated(current)) then
                do node_id = 1, min(num_nodes, size(current))
                    call self%pressure%get_history(node_id, history)
                    pressure_pred(node_id) = min(max(history(2) + ratio * (history(2) - history(3)), &
                                                     WALL_PRESS_MIN_PA), WALL_PRESS_MAX_PA)
                    current(node_id) = pressure_pred(node_id)
                end do
            end if
            nullify (current)
        end if

        if (self%is_active_thermal()) then
            select case (PREDICTOR_MODE)
            case (1)
                call predict_temperature_linear(self, num_nodes, ratio)
            case (2)
                call predict_temperature_clipped(self, num_nodes, ratio, pressure_pred)
            case (3)
                call predict_temperature_from_enthalpy(self, num_nodes, dt_new, ratio, pressure_pred)
            end select
        end if

        if (allocated(pressure_pred)) deallocate (pressure_pred)

        ! Everything derived from the primary variables has to follow the
        ! prediction, or the first residual is assembled from a mixed state.
        call self%calc_gradient_temperature()
        call self%calc_gradient_pressure()
        call self%update_nodal_phases()
        call self%control%get_bdf_coeffs(bdf_order, bdf_coeffs)
        if (self%is_active_thermal()) call self%temperature%compute_time_derivative(bdf_coeffs, bdf_order)
        if (self%is_active_hydraulic()) call self%pressure%compute_time_derivative(bdf_coeffs, bdf_order)
    end subroutine apply_solution_predictor

    !> Mode 1: straight line through the last two accepted temperatures.
    subroutine predict_temperature_linear(self, num_nodes, ratio)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: num_nodes
        real(real64), intent(in) :: ratio

        real(real64), pointer, contiguous :: current(:)
        real(real64) :: history(3)
        integer(int32) :: node_id

        nullify (current)
        call self%temperature%get_current(current)
        if (.not. associated(current)) return
        do node_id = 1, min(num_nodes, size(current))
            call self%temperature%get_history(node_id, history)
            current(node_id) = min(max(history(2) + ratio * (history(2) - history(3)), &
                                       WALL_TEMP_MIN_C), WALL_TEMP_MAX_C)
        end do
        nullify (current)
    end subroutine predict_temperature_linear

    !> Mode 2: the same line, stopped at the freezing point.
    !>
    !> A node still above its pressure-dependent freezing temperature is
    !> predicted no further than that temperature. The prediction may bring the
    !> front up to the transition but never through it, which is where the
    !> straight line stops being an extrapolation of anything: the latent heat
    !> holds T near T_crit while enthalpy keeps falling, so the slope on the far
    !> side bears no relation to the slope on this one. Crossing is left to the
    !> Newton correction, which has the apparent heat capacity to price it.
    subroutine predict_temperature_clipped(self, num_nodes, ratio, pressure_pred)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: num_nodes
        real(real64), intent(in) :: ratio
        real(real64), intent(in) :: pressure_pred(:)

        type(type_state) :: state
        real(real64), pointer, contiguous :: current(:)
        real(real64) :: history(3), predicted, rho_water, critical_temperature
        integer(int32) :: node_id, row_start, repr_elem

        nullify (current)
        call self%temperature%get_current(current)
        if (.not. associated(current)) return

        do node_id = 1, min(num_nodes, size(current))
            call self%temperature%get_history(node_id, history)
            predicted = min(max(history(2) + ratio * (history(2) - history(3)), &
                                WALL_TEMP_MIN_C), WALL_TEMP_MAX_C)

            if (predicted < history(2)) then
                row_start = self%node_material_table%ptr(node_id)
                if (row_start < self%node_material_table%ptr(node_id + 1)) then
                    repr_elem = self%node_material_table%repr_element(row_start)
                    call self%set_state(node_id, repr_elem, state, calc_physics=.false., include_fluxes=.false.)
                    call self%thermal%calc_density_water(state, rho_water)
                    call calc_T_high_celsius(pressure_pred(node_id), rho_water, critical_temperature)
                    if (history(2) > critical_temperature .and. predicted < critical_temperature) then
                        predicted = critical_temperature
                    end if
                end if
            end if
            current(node_id) = predicted
        end do
        nullify (current)
    end subroutine predict_temperature_clipped

    !> Mode 3: extrapolate the conserved quantity, then invert it for T.
    !>
    !> Enthalpy has no kink at the transition - that is what makes the enthalpy
    !> formulation well behaved - so a straight line through it is a sound
    !> extrapolation everywhere, including across freezing. Inverting
    !> U(T, p) = H_pred at the predicted pressure then places T wherever the
    !> latent heat says it belongs: a node absorbing latent heat has its
    !> enthalpy advanced while its temperature stays near T_crit, which is the
    !> behaviour a temperature extrapolation cannot produce.
    !>
    !> The enthalpy history is the one the local-error estimator already keeps.
    subroutine predict_temperature_from_enthalpy(self, num_nodes, dt_new, ratio, pressure_pred)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: num_nodes
        real(real64), intent(in) :: dt_new, ratio
        real(real64), intent(in) :: pressure_pred(:)

        integer(int32), parameter :: MAX_INVERSION_ITER = 20
        type(type_state) :: state
        real(real64), pointer, contiguous :: current(:)
        real(real64) :: history(3), target_enthalpy, temperature, enthalpy, capacity
        real(real64) :: increment, tolerance
        integer(int32) :: node_id, row_start, repr_elem, material_id, iter
        logical :: converged

        if (.not. allocated(self%lte_state_prev_thermal)) then
            call predict_temperature_clipped(self, num_nodes, ratio, pressure_pred)
            return
        end if
        if (.not. allocated(self%lte_ydot_prev_thermal)) then
            call predict_temperature_clipped(self, num_nodes, ratio, pressure_pred)
            return
        end if

        nullify (current)
        call self%temperature%get_current(current)
        if (.not. associated(current)) return
        if (size(self%lte_state_prev_thermal) /= size(current)) then
            call predict_temperature_clipped(self, num_nodes, ratio, pressure_pred)
            return
        end if

        do node_id = 1, min(num_nodes, size(current))
            call self%temperature%get_history(node_id, history)
            target_enthalpy = self%lte_state_prev_thermal(node_id) + &
                              dt_new * self%lte_ydot_prev_thermal(node_id)

            row_start = self%node_material_table%ptr(node_id)
            if (row_start >= self%node_material_table%ptr(node_id + 1)) cycle
            repr_elem = self%node_material_table%repr_element(row_start)
            call self%domain%get_material_id(repr_elem, material_id)
            call self%set_state(node_id, repr_elem, state, calc_physics=.false., include_fluxes=.false.)
            call state%pressure%set(pressure_pred(node_id))

            ! Absolute enthalpy tolerance: the same scale the acceptance gate
            ! uses for the thermal conserved quantity.
            tolerance = 1.0d0
            temperature = history(2)
            converged = .false.
            do iter = 1, MAX_INVERSION_ITER
                call state%temperature%set(temperature)
                call self%thermal%update_water_phases(material_id, state)
                call self%thermal%calc_enthalpy_density(material_id, state, enthalpy)
                if (abs(enthalpy - target_enthalpy) <= tolerance) then
                    converged = .true.
                    exit
                end if
                call self%thermal%compute_mass_term(material_id, state, capacity)
                if (abs(capacity) <= tiny(1.0d0)) exit
                increment = -(enthalpy - target_enthalpy) / capacity
                temperature = min(max(temperature + increment, WALL_TEMP_MIN_C), WALL_TEMP_MAX_C)
            end do

            if (converged) then
                current(node_id) = temperature
            else
                current(node_id) = history(2)
            end if
        end do
        nullify (current)
    end subroutine predict_temperature_from_enthalpy

    !> Extrapolate the accepted ice history to provide the next outer-loop
    !> initial guess. The bounded projection and the monolithic T-p solve still
    !> determine the accepted state, so this changes iteration count only.
    subroutine apply_phase_predictor(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        integer(int32) :: node_id
        real(real64) :: dt, ratio, predicted_increment, ice_history(3)
        real(real64), pointer, contiguous, dimension(:) :: current_ice, current_porosity

        if (self%last_accepted_dt <= tiny(1.0d0)) return

        call self%control%get_dt(dt)
        ratio = min(max(dt / self%last_accepted_dt, 0.0d0), 1.0d0)
        nullify (current_ice, current_porosity)
        call self%Qi%get_current(current_ice)
        call self%porosity%get_current(current_porosity)
        do node_id = 1, size(current_ice)
            call self%Qi%get_history(node_id, ice_history)
            predicted_increment = ratio * (ice_history(2) - ice_history(3))
            current_ice(node_id) = min(max(ice_history(2) + predicted_increment, 0.0d0), &
                                       current_porosity(node_id))
        end do
        nullify (current_ice, current_porosity)

        call self%update_nodal_phases()
    end subroutine apply_phase_predictor

    module subroutine solve_time_step_setup_ftcms(self, prescribe_bc)
        implicit none
        class(type_ftcms), intent(inout) :: self
        logical, intent(inout) :: prescribe_bc

        integer(int32) :: iter

        call self%control%increment_nonlinear()
        call self%control%get_nonlinear_iter(iter)

        if (iter == 1) then
            prescribe_bc = .true.
        else
            prescribe_bc = .false.
        end if

        ! No gradient recomputation here: the nodal T/P gradients are already
        ! current at this point. Every T/P mutation is followed by a
        ! projection -- reflect_variables recomputes both gradients after the
        ! update (end of each nonlinear iteration), and the first iteration
        ! recomputes them right after prescribe_dirichlet (before assembly).
        ! Recomputing from unchanged inputs here produced identical values.

    end subroutine solve_time_step_setup_ftcms

    module subroutine solve_time_step_check_convergence_ftcms(self, target_physics)
        implicit none
        class(type_ftcms), intent(inout), target :: self
        type(type_constant_id), intent(in), optional :: target_physics

        real(real64), pointer, contiguous, dimension(:) :: current_value

        real(real64), allocatable :: residual(:)
        real(real64), allocatable :: increment(:)
        real(real64) :: relaxation_factor
        logical, parameter :: diverged = .true.

        logical :: check_thermal, check_hydraulic

        nullify (current_value)

        ! Conserved-quantity mode runs its (post-update) check separately; the
        ! residual/update per-physics check below does not apply to it.
        if (self%control%is_conserved()) return

        check_thermal = self%control%is_physics_active(PHYSICS_TYPES%THERMAL)
        check_hydraulic = self%control%is_physics_active(PHYSICS_TYPES%HYDRAULIC)

        if (present(target_physics)) then
            check_thermal = check_thermal .and. (target_physics%ID == PHYSICS_TYPES%THERMAL%ID)
            check_hydraulic = check_hydraulic .and. (target_physics%ID == PHYSICS_TYPES%HYDRAULIC%ID)

            if (.not. check_thermal .and. self%control%is_physics_active(PHYSICS_TYPES%THERMAL)) then
                call self%control%set_converged(PHYSICS_TYPES%THERMAL, .true.)
                call self%control%set_diverged(PHYSICS_TYPES%THERMAL, .false.)
            end if
            if (.not. check_hydraulic .and. self%control%is_physics_active(PHYSICS_TYPES%HYDRAULIC)) then
                call self%control%set_converged(PHYSICS_TYPES%HYDRAULIC, .true.)
                call self%control%set_diverged(PHYSICS_TYPES%HYDRAULIC, .false.)
            end if
        end if

        ! ----------------------------------------------------------------------
        ! Thermal Convergence Check
        ! ----------------------------------------------------------------------
        if (check_thermal) then
            call self%get_variable_residual(PHYSICS_TYPES%THERMAL, residual)
            call self%get_variable_increment(PHYSICS_TYPES%THERMAL, increment)

            if (.not. allocated(residual) .or. .not. allocated(increment) .or. &
                size(residual) == 0 .or. size(increment) == 0) then
                write (*, *) "Error: Thermal residual/increment is unavailable during convergence check."
                call self%control%set_diverged(PHYSICS_TYPES%THERMAL, diverged)
            else if (has_nan(residual) .or. has_nan(increment)) then
                write (*, *) "Error: NaN detected in thermal variables during convergence check."
                call self%control%set_diverged(PHYSICS_TYPES%THERMAL, diverged)
            else
                call self%control%get_current_relaxation(PHYSICS_TYPES%THERMAL, relaxation_factor)
                increment(:) = relaxation_factor * increment(:)
                call self%control%check_convergence(PHYSICS_TYPES%THERMAL, residual, increment)
            end if

        end if

        ! ----------------------------------------------------------------------
        ! Hydraulic Convergence Check
        ! ----------------------------------------------------------------------
        if (check_hydraulic) then
            call self%get_variable_residual(PHYSICS_TYPES%HYDRAULIC, residual)
            call self%get_variable_increment(PHYSICS_TYPES%HYDRAULIC, increment)

            if (.not. allocated(residual) .or. .not. allocated(increment) .or. &
                size(residual) == 0 .or. size(increment) == 0) then
                write (*, *) "Error: Hydraulic residual/increment is unavailable during convergence check."
                call self%control%set_diverged(PHYSICS_TYPES%HYDRAULIC, diverged)
            else if (has_nan(residual) .or. has_nan(increment)) then
                write (*, *) "Error: NaN detected in hydraulic variables during convergence check."
                call self%control%set_diverged(PHYSICS_TYPES%HYDRAULIC, diverged)
            else
                call self%control%get_current_relaxation(PHYSICS_TYPES%HYDRAULIC, relaxation_factor)
                increment(:) = relaxation_factor * increment(:)
                call self%control%check_convergence(PHYSICS_TYPES%HYDRAULIC, residual, increment)
            end if
        end if

        if (allocated(increment)) call deallocate_array(increment)
        if (allocated(residual)) call deallocate_array(residual)

    end subroutine solve_time_step_check_convergence_ftcms

    !> Conserved-quantity convergence check (PDF 6.2.4), evaluated on the updated
    !> state. Builds the nodal enthalpy/effective-density fields and the per-block
    !> residual vectors, then delegates the coupled decision to the control manager.
    module subroutine solve_time_step_check_convergence_conserved_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        real(real64), allocatable :: enthalpy(:)
        real(real64), allocatable :: density(:)
        real(real64), allocatable :: dH_dT(:)
        real(real64), allocatable :: drho_dp(:)
        real(real64), allocatable :: residual_thermal(:)
        real(real64), allocatable :: residual_hydraulic(:)
        logical :: check_thermal, check_hydraulic
        logical :: hit_wall
        real(real64), pointer, contiguous, dimension(:) :: field
        ! Wall-contact tolerances: the bounded update lands exactly ON the wall
        ! when it binds, so a small absolute buffer is sufficient to detect it.
        real(real64), parameter :: TEMP_WALL_TOL = 1.0d-6
        real(real64), parameter :: PRESS_WALL_TOL = 1.0d-1

        check_thermal = self%is_active_thermal()
        check_hydraulic = self%is_active_hydraulic()

        ! Physical-validity guard: an iterate pinned at the sanity walls by the
        ! bounded update (reflect_variables) is outside the model's validity.
        ! Such a state must never be accepted as converged; declare divergence
        ! so the step fails and the ATS retries with a smaller dt.
        nullify (field)
        if (check_thermal) then
            call self%temperature%get_current(field)
            ! The verdict is reduced before it is acted on: returning here on a
            ! rank-local test would leave the other ranks waiting in the
            ! collective convergence gate below.
            hit_wall = .false.
            if (associated(field)) then
                hit_wall = (minval(field) <= WALL_TEMP_MIN_C + TEMP_WALL_TOL .or. &
                            maxval(field) >= WALL_TEMP_MAX_C - TEMP_WALL_TOL)
            end if
            if (reduce_any(hit_wall)) then
                if (hit_wall) then
                    write (*, '(A,2(ES11.3,1X))') '   [GUARD] temperature pinned at validity wall; T min/max = ', &
                        minval(field), maxval(field)
                end if
                call self%control%set_converged(PHYSICS_TYPES%THERMAL, .false.)
                call self%control%set_diverged(PHYSICS_TYPES%THERMAL, .true.)
                nullify (field)
                return
            end if
            nullify (field)
        end if
        if (check_hydraulic) then
            call self%pressure%get_current(field)
            hit_wall = .false.
            if (associated(field)) then
                hit_wall = (minval(field) <= WALL_PRESS_MIN_PA + PRESS_WALL_TOL .or. &
                            maxval(field) >= WALL_PRESS_MAX_PA - PRESS_WALL_TOL)
            end if
            ! Same reason as the thermal wall above: reduce, then act.
            if (reduce_any(hit_wall)) then
                if (hit_wall) then
                    block
                        integer(int32) :: node_lo, node_hi
                        real(real64), pointer, contiguous, dimension(:) :: T_dbg
                        real(real64) :: T_lo, T_hi
                        node_lo = minloc(field, dim=1)
                        node_hi = maxloc(field, dim=1)
                        T_lo = 0.0d0
                        T_hi = 0.0d0
                        nullify (T_dbg)
                        call self%temperature%get_current(T_dbg)
                        if (associated(T_dbg)) then
                            T_lo = T_dbg(node_lo)
                            T_hi = T_dbg(node_hi)
                            nullify (T_dbg)
                        end if
                        write (*, '(A,2(ES11.3,1X),A,I0,A,ES11.3,A,I0,A,ES11.3)') &
                            '   [GUARD] pressure pinned at validity wall; P min/max = ', &
                            minval(field), maxval(field), ' node_min=', node_lo, ' T=', T_lo, &
                            ' node_max=', node_hi, ' T=', T_hi
                    end block
                end if
                call self%control%set_converged(PHYSICS_TYPES%HYDRAULIC, .false.)
                call self%control%set_diverged(PHYSICS_TYPES%HYDRAULIC, .true.)
                nullify (field)
                return
            end if
            nullify (field)
        end if

        ! Nodal conserved quantities at the updated iterate, together with the
        ! matching nodal storage sensitivities (dH/dT, d rho_eq/dp) the
        ! primary-variable local error measure needs (local_error_block).
        ! Both are as state-dependent as enthalpy/density themselves, so they
        ! are recomputed at the same cadence, here, rather than once at setup.
        call self%compute_nodal_conserved(enthalpy, density, dH_dT, drho_dp)

        ! Per-block residuals from the assembly at the updated iterate.
        if (check_thermal) call self%get_variable_residual(PHYSICS_TYPES%THERMAL, residual_thermal)
        if (check_hydraulic) call self%get_variable_residual(PHYSICS_TYPES%HYDRAULIC, residual_hydraulic)

        ! Refresh the convergence gate's cached sensitivity/primary-variable
        ! snapshot before the gate below reads them. nodal_volume and the ATS
        ! tolerances are left untouched: unsupplied optional arguments keep
        ! the values set once at solve_time_step_initial_setup. A disassociated
        ! pointer (inactive physics) is treated as an absent optional argument.
        block
            real(real64), pointer, contiguous, dimension(:) :: T_current, p_current
            real(real64) :: dt_now, span_thermal, span_hydraulic

            nullify (T_current, p_current)
            if (check_thermal) call self%temperature%get_current(T_current)
            if (check_hydraulic) call self%pressure%get_current(p_current)
            call self%control%get_dt(dt_now)
            ! Reduced here, not inside the gate: the span is a cross-rank
            ! quantity and convergence_control carries no MPI.
            span_thermal = conserved_span(enthalpy)
            span_hydraulic = conserved_span(density)
            call self%control%set_residual_scale(self%domain_measure_total, dt_now, &
                                                 dH_dT=dH_dT, drho_dp=drho_dp, &
                                                 u_thermal=T_current, u_hydraulic=p_current, &
                                                 q_thermal=enthalpy, q_hydraulic=density, &
                                                 q_span_thermal=span_thermal, &
                                                 q_span_hydraulic=span_hydraulic)
            nullify (T_current, p_current)
        end block

        ! Unallocated residual arrays propagate as absent optional arguments.
        call self%control%check_convergence_conserved(enthalpy, density, &
                                                      residual_thermal, residual_hydraulic, &
                                                      check_thermal, check_hydraulic)

        ! Pore-overflow guard: the accepted nodal water+ice content must not
        ! exceed the pore space. Transient iterates within the nonlinear loop
        ! are allowed to overshoot (phase_systems.F90 no longer clamps them),
        ! so the check runs only once the convergence gate above has declared
        ! the iterate acceptable; a violation then vetoes the acceptance via
        ! the same set_diverged/dt-retry mechanism as the wall-pinning guards.
        if (check_hydraulic .and. self%control%is_converged()) then
            block
                real(real64), pointer, contiguous, dimension(:) :: Qw_field, Qi_field, porosity_field
                real(real64) :: excess, worst_excess
                integer(int32) :: i_node, worst_node
                real(real64), parameter :: PORE_OVERFLOW_TOL = 1.0d-8
                ! Diagnostic switch: report the overflow but do not veto the
                ! acceptance. The rigid-pore model has no pressure-relief
                ! physics, so freezing a near-saturated node necessarily
                ! overflows the pore volume by up to (rho_w/rho_i - 1)*Theta;
                ! warn-only runs measure how large that overflow actually
                ! grows before deciding the rejection policy.
                logical, parameter :: PORE_OVERFLOW_WARN_ONLY = .false.

                nullify (Qw_field, Qi_field, porosity_field)
                call self%Qw%get_current(Qw_field)
                call self%Qi%get_current(Qi_field)
                call self%porosity%get_current(porosity_field)
                if (associated(Qw_field) .and. associated(Qi_field) .and. associated(porosity_field)) then
                    worst_excess = PORE_OVERFLOW_TOL
                    worst_node = 0
                    do i_node = 1, size(porosity_field)
                        excess = Qw_field(i_node) + Qi_field(i_node) - porosity_field(i_node)
                        if (excess > worst_excess) then
                            worst_excess = excess
                            worst_node = i_node
                        end if
                    end do
                    if (worst_node > 0) then
                        write (*, '(A,I0,A,ES11.3)') &
                            '   [GUARD] pore-volume overflow at accepted step; node=', worst_node, &
                            ' excess=', worst_excess
                        if (.not. PORE_OVERFLOW_WARN_ONLY) then
                            call self%control%set_converged(PHYSICS_TYPES%HYDRAULIC, .false.)
                            call self%control%set_diverged(PHYSICS_TYPES%HYDRAULIC, .true.)
                        end if
                    end if
                end if
                nullify (Qw_field, Qi_field, porosity_field)
            end block
        end if

        if (allocated(enthalpy)) deallocate (enthalpy)
        if (allocated(density)) deallocate (density)
        if (allocated(dH_dT)) deallocate (dH_dT)
        if (allocated(drho_dp)) deallocate (drho_dp)
        if (allocated(residual_thermal)) deallocate (residual_thermal)
        if (allocated(residual_hydraulic)) deallocate (residual_hydraulic)
    end subroutine solve_time_step_check_convergence_conserved_ftcms

    module subroutine solve_time_step_ftcms(self, is_step_converged)
        implicit none
        class(type_ftcms), intent(inout) :: self
        logical, intent(inout) :: is_step_converged
        logical :: prescribe_bc
        integer(int32) :: iter_nl
        real(real64) :: t_res, t_inc, h_res, h_inc

        ! The primary unknowns remain T and p. Ice is eliminated locally by a
        ! bounded constitutive return map after each coupled Modified Picard
        ! update; no ice DOF is added to the global block structure.
        logical :: do_phase_outer
        integer(int32) :: coupling_iter, node_id
        ! Keep the former outer solver available as a diagnostic fallback.
        ! The production path applies the same local return map at every
        ! Modified Picard iterate, together with its T-p tangent assembled in
        ! phase_systems.F90. Fully converging a fixed-ice T-p problem before
        ! changing ice discards that tangent and causes the measured period-two
        ! front oscillation.
        logical, parameter :: PHASE_USE_OUTER_PROJECTION = .false.
        ! Allow a contracting phase map to finish before ATS reduces dt.
        integer(int32), parameter :: MAX_PHASE_ITER = 240
        ! 5 kPa corresponds to about 4e-3 K through Clapeyron near 0 C,
        ! below the verified spatial temperature discretization error.
        real(real64), parameter :: PHASE_PRESSURE_TOL = 5.0d3
        ! The phase equation is a water-content equality. This tolerance is
        ! consistent with the configured hydraulic conserved-quantity scale.
        real(real64), parameter :: PHASE_CONTENT_TOL = 1.0d-3
        ! A local defect must not be hidden by the volume-weighted norm at the
        ! thin freezing front. Use the same physical content scale locally and
        ! globally; otherwise an ice-free node could accept a finite freezing
        ! target solely because its support is small.
        real(real64), parameter :: PHASE_CONTENT_MAX_TOL = PHASE_CONTENT_TOL
        real(real64), parameter :: PHASE_MIXING = 0.3d0
        real(real64), parameter :: PHASE_AA_GAMMA_MAX = 2.0d0
        real(real64), parameter :: PHASE_AA_STEP_GROWTH_MAX = 4.0d0
        real(real64), parameter :: PHASE_STEP_FLOOR = 1.0d-5
        integer(int32), parameter :: PHASE_STAGNATION_LIMIT = 40
        real(real64) :: phase_increment_max, phase_increment_norm
        real(real64) :: phase_temperature, phase_pressure
        real(real64) :: phase_current_ice, phase_projected_ice
        real(real64) :: phase_equilibrium_error
        real(real64) :: phase_step_limit, phase_step_max, phase_step_factor
        real(real64) :: phase_mixing_current
        real(real64) :: phase_aa_gamma
        real(real64) :: phase_merit, phase_best_merit
        ! --- Backtracking line-search state ---
        integer(int32) :: ls_trial
        real(real64) :: ls_scale, ls_reference_norm
        !> Per-block local error measure at the trial iterate (E_T, E_H; see
        !> local_error_block) and the mean-volume-floor fallback scale used
        !> only when the exact per-node data is not yet available.
        real(real64) :: ls_E_T, ls_E_H
        real(real64) :: ls_ref0_thermal, ls_ref0_hydraulic
        logical :: ls_accepted
        real(real64), allocatable :: initial_residual_thermal(:), initial_residual_hydraulic(:)
        !> Per-node frozen flags, used only by the active-set telemetry.
        logical, allocatable :: frozen_now(:), frozen_previous(:)
        real(real64), allocatable :: phase_increments(:), previous_phase_increments(:)
        real(real64), allocatable :: phase_update(:), previous_phase_update(:)
        integer(int32), allocatable :: phase_active_bounds(:)
        integer(int32) :: phase_max_node, num_phase_nodes
        integer(int32) :: num_active_nodes
        integer(int32) :: phase_stagnation_count
        logical :: linear_failed, phase_is_converged, phase_aa_usable
        logical :: phase_final_correction_applied
        is_step_converged = .false.
        ls_ref0_thermal = 0.0d0
        ls_ref0_hydraulic = 0.0d0
        ls_E_T = 0.0d0
        ls_E_H = 0.0d0

        self%last_phase_iterations = 1
        self%last_inner_iterations = 0
        self%last_nonlinear_work = 0
        self%last_solve_status = SOLVE_STATUS_NOT_RUN
        self%last_phase_metrics_available = .false.
        self%last_phase_increment_max = -1.0d0
        ! Attempt-scoped: a linear failure or an early divergence return can end
        ! the attempt before any iterate is reported, and the history record
        ! would otherwise carry the previous attempt's increments.
        self%last_du_thermal_max = -1.0d0
        self%last_du_hydraulic_max = -1.0d0

        ! Dispatch to true staggered solver (H then T sequential nonlinear loops)
        if (self%control%is_staggered() .and. &
            self%is_active_thermal() .and. self%is_active_hydraulic()) then
            call self%solve_time_step_staggered(is_step_converged)
            return
        end if

        do_phase_outer = self%is_active_thermal() .and. self%is_active_hydraulic() .and. PHASE_USE_OUTER_PROJECTION
        phase_final_correction_applied = .false.

        ! Initialize per-time-step state only once.
        call self%solve_time_step_initial_setup()
        if (PREDICTOR_MODE > 0) call apply_solution_predictor(self)
        if (do_phase_outer) call apply_phase_predictor(self)
        call self%domain%get_num_nodes(num_phase_nodes)
        allocate (phase_update(num_phase_nodes), source=0.0d0)
        allocate (frozen_now(num_phase_nodes), source=.false.)
        allocate (frozen_previous(num_phase_nodes), source=.false.)
        allocate (previous_phase_update(num_phase_nodes), source=0.0d0)
        phase_mixing_current = PHASE_MIXING
        phase_best_merit = huge(1.0d0)
        phase_stagnation_count = 0

        ! Homotopy: march lambda from 0 towards 1 to grow an initial guess for
        ! the accepted solve. No convergence gate is evaluated here - the only
        ! acceptance test is the nonlinear loop below, which runs at lambda = 1.
        call self%control%begin_homotopy()
        if (self%control%is_homotopy_active()) then
            call run_homotopy_stages(self, num_phase_nodes)
        end if
        call self%control%finish_homotopy()
        ! Re-arm the nonlinear controller. reset_iteration is mandatory: the
        ! stages leave a nonzero iteration count and a converged(:) state that
        ! would let the accepted solve exit at its first test, and the count
        ! also reaches the adaptive time step through control%update's max().
        call self%control%reset_iteration()
        call self%control%set_nonlinear_solver(NONLINEAR_SOLVER%PICARD)
        call self%control%reset_acceleration()
        self%aa_has_prev = .false.
        self%aa_gnorm_prev = -1.0d0

        ! Outer phase projection loop
        coupling_loop: do coupling_iter = 1, merge(MAX_PHASE_ITER, 1, do_phase_outer)
            self%last_phase_iterations = coupling_iter
            linear_failed = .false.

            ! For coupling iterations > 1, reset nonlinear controls only.
            if (coupling_iter > 1) then
                call self%control%reset_iteration()
                call self%control%set_nonlinear_solver(NONLINEAR_SOLVER%PICARD)
                call self%control%increment_total()
                call self%control%reset_acceleration()
                ! The local ice projection changes the fixed-point map between
                ! outer iterations. Retain the per-attempt AA telemetry, but
                ! never reuse a T-p secant pair across two different ice states.
                self%aa_has_prev = .false.
                self%aa_gnorm_prev = -1.0d0
            end if

            ! Nonlinear iteration loop
            nonlinear: do
                if (.not. self%control%should_continue()) exit nonlinear

                ! Setup (update iteration counter)
                call self%solve_time_step_setup(prescribe_bc)

                ! Prescribe Dirichlet values before assembly so gradients reflect BCs
                if (prescribe_bc) then
                    call self%prescribe_dirichlet()
                    call self%calc_gradient_temperature()
                    call self%calc_gradient_pressure()
                end if

                ! Assemble matrices and residual
                call self%assemble()

                ! Apply boundary conditions (natural + essential) to the linear system
                call self%apply_bc(prescribed=.false.)

                call self%control%get_nonlinear_iter(iter_nl)

                ! Relative residuals must be normalized by R(x_0), assembled
                ! before any nonlinear update in this fixed-phase subproblem.
                if (self%control%is_conserved() .and. iter_nl == 1) then
                    if (self%is_active_thermal()) then
                        call self%get_variable_residual(PHYSICS_TYPES%THERMAL, initial_residual_thermal)
                    end if
                    if (self%is_active_hydraulic()) then
                        call self%get_variable_residual(PHYSICS_TYPES%HYDRAULIC, initial_residual_hydraulic)
                    end if
                    call self%control%prime_conserved_residual(initial_residual_thermal, &
                                                               initial_residual_hydraulic, &
                                                               self%is_active_thermal(), &
                                                               self%is_active_hydraulic())
                    ! Fallback merit scale for the line search, used only when
                    ! local_error_block's exact per-node measure is not yet
                    ! available (nodal_volume or the storage-sensitivity
                    ! snapshot missing - only possible during the first
                    ! iterations of the very first step, before
                    ! solve_time_step_check_convergence_conserved has run
                    ! once). It must NOT be the step-initial residual norm
                    ! itself: a block that starts the step already satisfied
                    ! has ||R^0|| at round-off - measured, ||R_H^0|| = 2.75e-11
                    ! against ||R_T^0|| = 2.48e-01, ten orders apart - so
                    ! dividing by it makes any later residual in that block
                    ! dominate the merit no matter how small it is physically.
                    ! Scanned along the Newton direction at a stalled iterate,
                    ! the thermal block fell 33 percent at alpha = 0.46 while
                    ! the ratio-normalized hydraulic block outweighed it 2000
                    ! to 1 and pinned the search at alpha = 0.007. The floor
                    ! below is the imbalance a block may carry over the step
                    ! without moving its conserved quantity past its own
                    ! absolute tolerance; it does not collapse when the block
                    ! starts converged.
                    ls_ref0_thermal = 0.0d0
                    ls_ref0_hydraulic = 0.0d0
                    call self%control%get_residual_floors(num_phase_nodes, ls_ref0_thermal, ls_ref0_hydraulic)
                    if (allocated(initial_residual_thermal)) deallocate (initial_residual_thermal)
                    if (allocated(initial_residual_hydraulic)) deallocate (initial_residual_hydraulic)
                end if

                ! Reference merit for the line search, taken BEFORE the linear
                ! solve: solve() calls jacobi_equilibrate_bsr, which scales K and
                ! F in place to fix the 1e13 T/p conditioning disparity. Only du
                ! is unscaled afterwards, so F is left equilibrated and reading it
                ! after the solve compares a scaled residual against a freshly
                ! assembled unscaled one.
                if (self%control%is_conserved()) then
                    call ls_block_norms(self, ls_E_T, ls_E_H, ls_ref0_thermal, ls_ref0_hydraulic)
                    ls_reference_norm = ls_merit(ls_E_T, ls_E_H)
                end if

                ! Linear solve (K * du = F)
                call self%solve()

                ! If linear solver failed, mark as diverged and exit
                if (.not. self%solver%is_success()) then
                    linear_failed = .true.
                    if (self%is_active_thermal()) then
                        call self%control%set_converged( &
                            PHYSICS_TYPES%THERMAL, .false.)
                        call self%control%set_diverged( &
                            PHYSICS_TYPES%THERMAL, .true.)
                    end if
                    if (self%is_active_hydraulic()) then
                        call self%control%set_converged( &
                            PHYSICS_TYPES%HYDRAULIC, .false.)
                        call self%control%set_diverged( &
                            PHYSICS_TYPES%HYDRAULIC, .true.)
                    end if
                    exit nonlinear
                end if

                ! Convergence check; always converged when config is NONE
                call self%solve_time_step_check_convergence()

                if (self%control%is_conserved()) then
                    ! Backtracking line search on the assembled residual norm.
                    !
                    ! The Picard direction is not guaranteed to be a descent
                    ! direction at full step where the freezing curve is steep,
                    ! and a heuristic that damps only after detecting stagnation
                    ! is always one step behind. Choosing the step length so that
                    ! ||R|| actually decreases makes progress monotone by
                    ! construction and, when no admissible step exists, says so
                    ! immediately - that distinguishes a step-length problem from
                    ! a wrong-direction problem instead of hiding both.
                    !
                    ! The system was assembled at x_k, so its residual norm is
                    ! the reference. Each trial restores x_k, applies the scaled
                    ! update, and reassembles: the reassembly was already being
                    ! done once per iteration, so a converging search costs
                    ! nothing extra and only hard iterates pay for retries.
                    if (RESIDUAL_SCAN .and. .not. residual_scan_done .and. iter_nl >= RESIDUAL_SCAN_ITER) then
                        residual_scan_done = .true.
                        call scan_residual_along_direction(self, iter_nl, ls_ref0_thermal, ls_ref0_hydraulic)
                    end if

                    call apply_line_search(self, iter_nl, ls_reference_norm, &
                                           ls_ref0_thermal, ls_ref0_hydraulic, &
                                           ls_accepted, ls_scale, ls_trial, ls_E_T, ls_E_H, &
                                           report_failure=.true.)
                    if (.not. ls_accepted) then
                        if (self%is_active_thermal()) then
                            call self%control%set_converged(PHYSICS_TYPES%THERMAL, .false.)
                            call self%control%set_diverged(PHYSICS_TYPES%THERMAL, .true.)
                        end if
                        if (self%is_active_hydraulic()) then
                            call self%control%set_converged(PHYSICS_TYPES%HYDRAULIC, .false.)
                            call self%control%set_diverged(PHYSICS_TYPES%HYDRAULIC, .true.)
                        end if
                        exit nonlinear
                    end if

                    ! Increment magnitudes, reported on the iteration line below.
                    block
                        real(real64), allocatable :: du_report_T(:), du_report_P(:)

                        call self%get_variable_increment(PHYSICS_TYPES%THERMAL, du_report_T)
                        call self%get_variable_increment(PHYSICS_TYPES%HYDRAULIC, du_report_P)
                        self%last_du_thermal_max = 0.0d0
                        self%last_du_hydraulic_max = 0.0d0
                        if (allocated(du_report_T)) self%last_du_thermal_max = maxval(abs(du_report_T))
                        if (allocated(du_report_P)) self%last_du_hydraulic_max = maxval(abs(du_report_P))
                    end block

                    if (iter_nl >= FREEZE_REPORT_ITER .and. freeze_report_count < FREEZE_REPORT_LIMIT) then
                        freeze_report_count = freeze_report_count + 1
                        call report_freezing_band(self, iter_nl)
                    end if

                    ! Active-set telemetry.
                    !
                    ! The element tangent (all four blocks, audited against the
                    ! finite-difference tangent) and the linear solve (defect
                    ! below 1e-8) are both exact, so a Newton step that still
                    ! needs 14 to 30 iterations cannot be explained by either.
                    ! What is left is that the residual is only piecewise
                    ! smooth: theta_i switches at T = 0, at the pore-volume
                    ! bound, and at the phase-content clamps. A finite-difference
                    ! probe cannot see those - its perturbation is far too small
                    ! to straddle one - but a full Newton step is not. Counting
                    ! how many nodes change frozen state per iteration separates
                    ! a set that settles (smooth Newton) from one that keeps
                    ! flipping (semismooth chattering).
                    block
                        real(real64), pointer, contiguous :: active_set_temperature(:)
                        integer(int32) :: node, num_flipped

                        nullify (active_set_temperature)
                        call self%temperature%get_current(active_set_temperature)
                        if (associated(active_set_temperature) .and. allocated(frozen_now)) then
                            num_flipped = 0
                            do node = 1, min(size(active_set_temperature), num_phase_nodes)
                                frozen_now(node) = active_set_temperature(node) < 0.0d0
                                if (iter_nl > 1) then
                                    if (frozen_now(node) .neqv. frozen_previous(node)) &
                                        num_flipped = num_flipped + 1
                                end if
                            end do
                            if (iter_nl > 1 .and. num_flipped > 0) then
                                write (*, '(A,I0,A,I0,A,I0)') '   [ACTIVE-SET] iter ', iter_nl, &
                                    ': frozen nodes=', count(frozen_now), ', flipped=', num_flipped
                            end if
                            frozen_previous(:) = frozen_now(:)
                        end if
                        nullify (active_set_temperature)
                    end block

                    ! theta_i is a state function of (T,p) evaluated inside the
                    ! reassembly above, so there is no separate phase criterion.
                    call self%solve_time_step_check_convergence_conserved()

                    if (NONLINEAR_VERBOSE) call report_nonlinear_iterate(self, iter_nl)
                else
                    ! Update solution with relaxation (Aitken for legacy Picard,
                    ! damped for Newton).
                    call self%reflect_variables()
                end if

                ! Force exit after one iteration when config is NONE (linear solve)
                if (self%control%is_none()) exit nonlinear

            end do nonlinear

            is_step_converged = self%control%is_converged()
            ! A NONE scheme is one linearised solve per step by definition, so
            ! there is no iteration whose convergence could be tested; asking
            ! the nonlinear gate would reject every step.
            if (self%control%is_none()) is_step_converged = .true.
            call self%control%get_nonlinear_iter(iter_nl)
            self%last_inner_iterations = iter_nl
            self%last_nonlinear_work = self%last_nonlinear_work + max(1_int32, iter_nl)

            if (.not. is_step_converged) then
                if (linear_failed) then
                    self%last_solve_status = SOLVE_STATUS_LINEAR_FAILURE
                else if (self%control%is_diverged()) then
                    self%last_solve_status = SOLVE_STATUS_NONLINEAR_DIVERGED
                else
                    self%last_solve_status = SOLVE_STATUS_NONLINEAR_LIMIT
                end if
                if (ERROR_NODE_REPORT .and. error_node_report_count < ERROR_NODE_REPORT_LIMIT) then
                    error_node_report_count = error_node_report_count + 1
                    block
                        real(real64), allocatable :: rep_residual(:)
                        real(real64), pointer, contiguous :: rep_qw(:), rep_qi(:), rep_phi(:)
                        nullify (rep_qw, rep_qi, rep_phi)
                        call self%Qw%get_current(rep_qw)
                        call self%Qi%get_current(rep_qi)
                        call self%porosity%get_current(rep_phi)
                        if (self%is_active_thermal()) then
                            call self%get_variable_residual(PHYSICS_TYPES%THERMAL, rep_residual)
                            if (allocated(rep_residual)) then
                                if (associated(rep_qw) .and. associated(rep_qi) .and. associated(rep_phi)) then
                                    call self%control%report_local_error_nodes(PHYSICS_TYPES%THERMAL, rep_residual, 'thermal', &
                                                                              theta_w=rep_qw, theta_i=rep_qi, porosity=rep_phi)
                                else
                                    call self%control%report_local_error_nodes(PHYSICS_TYPES%THERMAL, rep_residual, 'thermal')
                                end if
                                deallocate (rep_residual)
                            end if
                        end if
                        if (self%is_active_hydraulic()) then
                            call self%get_variable_residual(PHYSICS_TYPES%HYDRAULIC, rep_residual)
                            if (allocated(rep_residual)) then
                                call self%control%report_local_error_nodes(PHYSICS_TYPES%HYDRAULIC, rep_residual, 'hydraulic')
                                deallocate (rep_residual)
                            end if
                        end if
                        nullify (rep_qw, rep_qi, rep_phi)
                    end block
                end if
                ! One-shot tangent audit at the first attempt the nonlinear
                ! solve gives up on. That state - not the first line-search
                ! failure, which happens before any ice forms - is the one the
                ! linearization has to explain.
                if (FD_JACOBIAN_AUDIT .and. .not. fd_audit_done) then
                    fd_audit_done = .true.
                    call self%report_fd_jacobian('nonlinear failure')
                end if
                t_res = 0.0d0
                t_inc = 0.0d0
                h_res = 0.0d0
                h_inc = 0.0d0
                if (self%is_active_thermal()) then
                    if (.not. linear_failed) then
                        call self%control%get_current_norm(PHYSICS_TYPES%THERMAL, NONLINEAR_NORM_CRITERIA%RESIDUAL, &
                                                           NORM_TYPES%LINF, t_res)
                        call self%control%get_current_norm(PHYSICS_TYPES%THERMAL, NONLINEAR_NORM_CRITERIA%UPDATE, &
                                                           NORM_TYPES%LINF, t_inc)
                    end if
                end if
                if (self%is_active_hydraulic()) then
                    if (.not. linear_failed) then
                        call self%control%get_current_norm(PHYSICS_TYPES%HYDRAULIC, NONLINEAR_NORM_CRITERIA%RESIDUAL, &
                                                           NORM_TYPES%LINF, h_res)
                        call self%control%get_current_norm(PHYSICS_TYPES%HYDRAULIC, NONLINEAR_NORM_CRITERIA%UPDATE, &
                                                           NORM_TYPES%LINF, h_inc)
                    end if
                end if
                if (linear_failed) then
                    write (*, '(A,I0,A,L1,A)') '   [NONLINEAR] failed: iter=', iter_nl, ', diverged=', &
                        self%control%is_diverged(), ', linear solver failure before nonlinear norm update.'
                else
                    write (*, '(A,I0,A,L1,A,4(ES11.3,1X))') '   [NONLINEAR] failed: iter=', iter_nl, ', diverged=', &
                        self%control%is_diverged(), ', T_res/T_inc/H_res/H_inc=', t_res, t_inc, h_res, h_inc
                end if
            end if

            ! If inner solve failed, skip coupling check
            if (.not. is_step_converged) exit coupling_loop

            ! The local phase return was already included in every inner
            ! Modified Picard iterate. Refresh derived fields once and exit;
            ! no fixed-ice outer solve is needed.
            if (.not. do_phase_outer) then
                call self%update_nodal_phases()
                exit coupling_loop
            end if

            call self%project_nodal_ice(.false., phase_update, phase_increment_max, phase_increment_norm, &
                                        phase_max_node, phase_temperature, phase_pressure, &
                                        phase_current_ice, phase_projected_ice, &
                                        phase_equilibrium_error, phase_increments, phase_active_bounds)

            phase_aa_gamma = 0.0d0
            phase_aa_usable = .false.
            if (allocated(previous_phase_increments)) then
                call compute_aa1_coefficient(phase_increments, previous_phase_increments, &
                                             PHASE_AA_GAMMA_MAX, phase_aa_gamma, phase_aa_usable)
            end if
            if (phase_aa_usable) then
                phase_update = phase_mixing_current * phase_increments - phase_aa_gamma * ( &
                    previous_phase_update + phase_mixing_current * &
                    (phase_increments - previous_phase_increments))
            else
                phase_update = phase_mixing_current * phase_increments
            end if

            ! Bound only the accelerated step size. Do not reject a period-two
            ! residual merely because its current branch is larger or because
            ! AA(1) crosses a component-wise Picard direction: the predicted
            ! residual test in compute_aa1_coefficient is the safeguard.
            phase_step_max = maxval(abs(phase_update))
            phase_step_limit = max(PHASE_STEP_FLOOR, &
                                   PHASE_AA_STEP_GROWTH_MAX * phase_mixing_current * phase_increment_max)
            if (phase_step_max > phase_step_limit) phase_update = phase_update * phase_step_limit / phase_step_max

            ! Once a node is within the accepted local content width of an
            ! active ice bound, its semismooth equation is simply Qi=0 or
            ! Qi=Qi_max. Apply that exact projected-Newton step instead of
            ! geometrically damping it with the global AA coefficient. This
            ! removes artificial tail iterations without relaxing any phase
            ! equilibrium or conservation criterion.
            do node_id = 1, num_phase_nodes
                if (phase_active_bounds(node_id) /= 0 .and. &
                    abs(phase_increments(node_id)) <= PHASE_CONTENT_MAX_TOL) then
                    phase_update(node_id) = phase_increments(node_id)
                end if
            end do
            num_active_nodes = count(phase_active_bounds /= 0)
            phase_step_factor = 0.0d0
            if (abs(phase_increments(phase_max_node)) > tiny(1.0d0)) then
                phase_step_factor = phase_update(phase_max_node) / phase_increments(phase_max_node)
            end if
            ! ES for every unbounded field: a diverging iterate can push T or
            ! the step factor past any fixed F width, and an overflowing F
            ! descriptor aborts the run on an output conversion error.
            write (*, '(A,I0,A,ES10.3,A,ES10.3,A,I0,A,ES11.3,A,ES11.3,' // &
                        'A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,I0)') &
                '   [PHASE] outer:', coupling_iter, ' max|dQi|:', phase_increment_max, &
                ' rms|dQi|:', phase_increment_norm, ' node:', phase_max_node, &
                ' T:', phase_temperature, ' p:', phase_pressure, &
                ' Qi:', phase_current_ice, ' target:', phase_projected_ice, &
                ' eq[Pa]:', phase_equilibrium_error, ' omega:', phase_step_factor, &
                ' active_nodes:', num_active_nodes

            ! Both forms of the phase condition are required. In the flat part
            ! of a retention curve a small water-content defect can coexist with
            ! a large chemical-potential (Clapeyron pressure) disequilibrium.
            phase_is_converged = phase_equilibrium_error <= PHASE_PRESSURE_TOL .and. &
                                 phase_increment_norm <= PHASE_CONTENT_TOL .and. &
                                 phase_increment_max <= PHASE_CONTENT_MAX_TOL
            phase_merit = max(phase_equilibrium_error / PHASE_PRESSURE_TOL, &
                              phase_increment_norm / PHASE_CONTENT_TOL, &
                              phase_increment_max / PHASE_CONTENT_MAX_TOL)
            self%last_phase_metrics_available = .true.
            self%last_phase_increment_max = phase_increment_max
            if (phase_merit < 0.9d0 * phase_best_merit) then
                phase_best_merit = phase_merit
                phase_stagnation_count = 0
            else
                phase_stagnation_count = phase_stagnation_count + 1
            end if
            if (phase_is_converged .and. is_step_converged .and. &
                (phase_final_correction_applied .or. phase_increment_max <= tiny(1.0d0))) exit coupling_loop
            if (.not. phase_is_converged) phase_final_correction_applied = .false.

            if (coupling_iter == MAX_PHASE_ITER .or. &
                (phase_stagnation_count >= PHASE_STAGNATION_LIMIT .and. phase_best_merit > 1.0d0)) then
                write (*, '(A,ES11.3,A,ES11.3,A,ES11.3)') &
                    '   [PHASE] failed to reach local equilibrium; max|dQi|=', phase_increment_max, &
                    ', rms|dQi|=', phase_increment_norm, &
                    ', max pressure error [Pa]=', phase_equilibrium_error
                call self%control%set_converged(PHYSICS_TYPES%THERMAL, .false.)
                call self%control%set_converged(PHYSICS_TYPES%HYDRAULIC, .false.)
                call self%control%set_diverged(PHYSICS_TYPES%THERMAL, .true.)
                call self%control%set_diverged(PHYSICS_TYPES%HYDRAULIC, .true.)
                is_step_converged = .false.
                self%last_solve_status = SOLVE_STATUS_PHASE_FAILURE
                exit coupling_loop
            end if

            if (allocated(previous_phase_increments)) deallocate (previous_phase_increments)
            allocate (previous_phase_increments, source=phase_increments)
            phase_final_correction_applied = phase_is_converged
            call self%project_nodal_ice(.true., phase_update, phase_increment_max, phase_increment_norm, &
                                        phase_max_node, phase_temperature, phase_pressure, &
                                        phase_current_ice, phase_projected_ice, &
                                        phase_equilibrium_error, phase_increments, phase_active_bounds)
            previous_phase_update = phase_update
            call self%update_nodal_phases()

        end do coupling_loop

        if (is_step_converged) then
            self%last_solve_status = SOLVE_STATUS_CONVERGED
        else if (self%last_solve_status == SOLVE_STATUS_NOT_RUN) then
            if (self%control%is_diverged()) then
                self%last_solve_status = SOLVE_STATUS_NONLINEAR_DIVERGED
            else
                self%last_solve_status = SOLVE_STATUS_NONLINEAR_LIMIT
            end if
        end if

    end subroutine solve_time_step_ftcms

    !> Local nonlinear-error measure (E_T, E_H) of the thermal and hydraulic
    !> blocks, separately, at the current (trial) iterate.
    !>
    !> Global span max(q) - min(q) of a conserved field.
    !>
    !> Assumptions: q holds one value per local node.
    !> Numerical guarantees: the result is identical on every rank, so the
    !> acceptance gate weights the same way everywhere; a rank-local span
    !> would hand each rank a different criterion and desynchronize them.
    !> Computational complexity: O(n) locally plus one reduction.
    !> Failure behavior: returns zero for an unallocated or empty field.
    function conserved_span(q) result(span)
        implicit none
        !> Conserved quantity per node
        real(real64), allocatable, intent(in) :: q(:)
        !> Span in the quantity's own units
        real(real64) :: span

        real(real64) :: hi, lo
#ifdef _MPI
        real(real64) :: local_pair(2), global_pair(2)
        integer(int32) :: ierr
#endif

        span = 0.0d0
        if (.not. allocated(q)) return
        if (size(q) == 0) return

        hi = maxval(q)
        lo = minval(q)
#ifdef _MPI
        local_pair = [hi, -lo]
        call MPI_Allreduce(local_pair, global_pair, 2, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
        hi = global_pair(1)
        lo = -global_pair(2)
#endif
        span = max(0.0d0, hi - lo)
    end function conserved_span

    !> Delegates to convergence_control:local_error_block, the SAME formula
    !> and SAME cached per-node scales (nodal_volume, dH/dT or d rho_eq/dp,
    !> the primary-variable snapshot, the ATS tolerances) the conserved
    !> acceptance gate's local criterion (locT/locH in the [NL] line)
    !> evaluates, so the line search and the gate cannot disagree about what
    !> progress means. They must stay separate rather than combined into one
    !> raw residual norm: the thermal residual is an energy rate [W] and the
    !> hydraulic one a mass rate [kg/s], and on this problem they differ by
    !> about nine orders of magnitude, so a merit built from their raw sum of
    !> squares would be the thermal norm to machine precision and could not
    !> see the hydraulic block at all - local_error_block's per-block K/Pa
    !> scaling is what makes combining them in ls_merit meaningful.
    !>
    !> Falls back to the block residual ratioed against the mean-volume-floor
    !> scale (floor_thermal/floor_hydraulic, from get_residual_floors) when
    !> local_error_block reports its data as unavailable (sentinel < 0) -
    !> only possible during the first iterations of the very first step,
    !> before solve_time_step_check_convergence_conserved has run once.
    subroutine ls_block_norms(self, e_thermal, e_hydraulic, floor_thermal, floor_hydraulic)
        implicit none
        class(type_ftcms), intent(inout) :: self
        real(real64), intent(inout) :: e_thermal
        real(real64), intent(inout) :: e_hydraulic
        !> Mean-volume-floor fallback scale for each block (see get_residual_floors)
        real(real64), intent(in) :: floor_thermal, floor_hydraulic

        real(real64), allocatable :: r_thermal(:), r_hydraulic(:)
        real(real64) :: norm_thermal, norm_hydraulic
#ifdef _MPI
        real(real64) :: local_errors(2), global_errors(2)
        integer(int32) :: ierr
#endif

        e_thermal = 0.0d0
        e_hydraulic = 0.0d0
        if (self%is_active_thermal()) then
            call self%get_variable_residual(PHYSICS_TYPES%THERMAL, r_thermal)
            if (allocated(r_thermal)) then
                e_thermal = self%control%local_error_block(PHYSICS_TYPES%THERMAL, r_thermal)
                if (e_thermal < 0.0d0) then
                    norm_thermal = sqrt(dot_product(r_thermal, r_thermal))
                    e_thermal = norm_thermal / max(floor_thermal, tiny(1.0d0))
                end if
            end if
        end if
        if (self%is_active_hydraulic()) then
            call self%get_variable_residual(PHYSICS_TYPES%HYDRAULIC, r_hydraulic)
            if (allocated(r_hydraulic)) then
                e_hydraulic = self%control%local_error_block(PHYSICS_TYPES%HYDRAULIC, r_hydraulic)
                if (e_hydraulic < 0.0d0) then
                    norm_hydraulic = sqrt(dot_product(r_hydraulic, r_hydraulic))
                    e_hydraulic = norm_hydraulic / max(floor_hydraulic, tiny(1.0d0))
                end if
            end if
        end if
#ifdef _MPI
        local_errors = [e_thermal, e_hydraulic]
        call MPI_Allreduce(local_errors, global_errors, 2, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
        e_thermal = global_errors(1)
        e_hydraulic = global_errors(2)
#endif
        if (allocated(r_thermal)) deallocate (r_thermal)
        if (allocated(r_hydraulic)) deallocate (r_hydraulic)
    end subroutine ls_block_norms

    !> Dimensionless line-search merit: L2 combination of the two blocks'
    !> local error measures E_T, E_H (ls_block_norms).
    !>
    !> This is deliberately the same quantity the conserved acceptance gate's
    !> local criterion evaluates per block (locT/locH <= 1), so the search and
    !> the gate cannot disagree about what progress means. With the previous
    !> raw-residual-over-R0 merit the search rejected every step that reduced
    !> a hydraulic residual sitting at resH/0 = 3.12 against a tolerance of
    !> 0.1, because the same step raised the nine-orders-larger thermal block
    !> by 1.6-1.9 percent; it then exhausted its backtracking budget and left
    !> alpha = 7.8e-3, which advances nothing. Measured: the iterate was
    !> bit-for-bit unchanged across successive nonlinear iterations while
    !> resT/0 stayed at 0.1145.
    pure function ls_merit(e_thermal, e_hydraulic) result(merit)
        implicit none
        real(real64), intent(in) :: e_thermal, e_hydraulic
        real(real64) :: merit

        merit = sqrt(e_thermal**2 + e_hydraulic**2)
    end function ls_merit

    !> True when at least one node is cold enough for the generalized-Clapeyron
    !> suction to be nonzero, i.e. when lambda can change the physics at all.
    !>
    !> Assumptions: the freezing suction vanishes above the bulk freezing point,
    !> which holds for the non-segregation model; a pressure-dependent model
    !> would need its own test.
    !> Computational complexity: O(n_node) arithmetic, no allocation.
    !> Failure behavior: returns .true. when the field cannot be read, so an
    !> unexpected state costs work rather than silently skipping continuation.
    function homotopy_can_bite(self) result(can_bite)
        implicit none
        !> Solver object
        class(type_ftcms), intent(inout) :: self
        !> True when the ladder can change anything
        logical :: can_bite

        real(real64), pointer, contiguous :: field(:)

        can_bite = .true.
        if (.not. self%is_active_thermal()) return

        nullify (field)
        call self%temperature%get_current(field)
        if (.not. associated(field)) return
        can_bite = minval(field) <= HOMOTOPY_ONSET_MARGIN
        nullify (field)
    end function homotopy_can_bite

    !> March the continuation parameter from 0 towards 1, spending a bounded
    !> number of nonlinear iterations at each value.
    !>
    !> The stages are predictor work, not solves: no convergence gate is
    !> evaluated, and a stage is judged only on whether the line search stayed
    !> healthy and the merit did not grow. A stage that fails is rolled back to
    !> the last banked lambda and retried with half the increment; when the
    !> increment collapses the ladder gives up and leaves the iterate where it
    !> got to, which is still a better starting point than the predictor.
    !>
    !> Rollback is transactional: restoring T and p is not enough, because the
    !> quadrature depth and the assembled system belong to the rejected lambda,
    !> so the state is re-derived and fully re-assembled at the banked value.
    subroutine run_homotopy_stages(self, num_nodes)
        implicit none
        !> Solver object; its iterate and homotopy telemetry are updated
        class(type_ftcms), intent(inout) :: self
        !> Number of domain nodes, for the physical residual floors
        integer(int32), intent(in) :: num_nodes

        real(real64), allocatable :: hc_T(:), hc_P(:)
        real(real64), allocatable :: hc_aa_T(:), hc_aa_P(:), hc_aa_duT(:), hc_aa_duP(:)
        logical, allocatable :: hc_onset(:)
        logical :: hc_has_prev, prescribe_bc, accepted, stage_ok
        real(real64) :: hc_gnorm
        real(real64) :: floor_thermal, floor_hydraulic
        real(real64) :: e_thermal, e_hydraulic
        real(real64) :: reference_norm, entry_merit, exit_merit
        real(real64) :: scale, lambda
        integer(int32) :: stage_iter, trials, stage_work

        ! Nothing is cold enough for the freezing suction to be nonzero, so the
        ! blend is the identity at every lambda and every stage would re-solve
        ! one unchanging problem.
        if (.not. homotopy_can_bite(self)) then
            self%last_phase_iterations = 1
            return
        end if

        floor_thermal = 0.0d0
        floor_hydraulic = 0.0d0
        call self%control%get_residual_floors(num_nodes, floor_thermal, floor_hydraulic)
        stage_work = 0

        stage_loop: do
            if (self%control%is_homotopy_complete()) exit stage_loop
            if (self%control%is_homotopy_exhausted()) exit stage_loop
            lambda = self%control%get_homotopy_lambda()

            call ls_snapshot(self, hc_T, hc_P, hc_onset, &
                             hc_aa_T, hc_aa_P, hc_aa_duT, hc_aa_duP, hc_has_prev, hc_gnorm)

            stage_ok = .true.
            entry_merit = -1.0d0
            exit_merit = -1.0d0
            do stage_iter = 1, HOMOTOPY_STAGE_ITERATIONS
                call self%solve_time_step_setup(prescribe_bc)
                if (prescribe_bc) then
                    call self%prescribe_dirichlet()
                    call self%calc_gradient_temperature()
                    call self%calc_gradient_pressure()
                end if
                call self%assemble()
                call self%apply_bc(prescribed=.false.)

                ! Reference merit before solve(): the linear solve equilibrates
                ! K and F in place and only unscales du.
                call ls_block_norms(self, e_thermal, e_hydraulic, floor_thermal, floor_hydraulic)
                reference_norm = ls_merit(e_thermal, e_hydraulic)
                if (entry_merit < 0.0d0) entry_merit = reference_norm

                call self%solve()
                if (.not. self%solver%is_success()) then
                    stage_ok = .false.
                    exit
                end if

                call apply_line_search(self, stage_iter, reference_norm, &
                                       floor_thermal, floor_hydraulic, &
                                       accepted, scale, trials, e_thermal, e_hydraulic, &
                                       report_failure=.false.)
                stage_work = stage_work + 1
                if (.not. accepted .or. scale < HOMOTOPY_STAGE_MIN_SCALE) then
                    stage_ok = .false.
                    exit
                end if
                exit_merit = ls_merit(e_thermal, e_hydraulic)
                ! Spending the rest of the budget on an iterate that has stopped
                ! moving buys no better initial guess.
                if (exit_merit > HOMOTOPY_STAGE_PROGRESS * reference_norm) exit
            end do

            if (stage_ok .and. entry_merit > 0.0d0 .and. exit_merit > entry_merit) stage_ok = .false.

            if (stage_ok) then
                call self%control%accept_homotopy_stage()
            else
                call ls_restore(self, hc_T, hc_P, hc_onset, &
                                hc_aa_T, hc_aa_P, hc_aa_duT, hc_aa_duP, hc_has_prev, hc_gnorm)
                call self%reflect_variables(step_scale=0.0d0)
                call self%control%reject_homotopy_stage()
                ! depth and the matrix blocks still belong to the rejected one.
                call self%assemble()
                call self%apply_bc(prescribed=.false.)
            end if
        end do stage_loop

        ! Stage work is predictor cost, not step difficulty. It must not reach
        ! the adaptive time step through last_inner_iterations.
        self%last_phase_iterations = max(1, self%control%get_homotopy_stage())
        self%last_nonlinear_work = self%last_nonlinear_work + stage_work
        if (HOMOTOPY_VERBOSE) then
            write (*, '(A,I0,A,ES10.3,A,I0)') '   [HOM] stages=', self%control%get_homotopy_stage(), &
                ', lambda=', self%control%get_homotopy_lambda(), ', iterations=', stage_work
        end if
    end subroutine run_homotopy_stages

    !> Backtracking Armijo search on the assembled residual merit.
    !>
    !> The system was assembled at x_k, so its merit is the reference. Each
    !> trial restores x_k, applies the scaled update and reassembles, so a
    !> converging search costs nothing over the reassembly the iteration was
    !> doing anyway and only hard iterates pay for retries.
    !>
    !> On failure the iterate is put back where the search started: leaving a
    !> rejected trial in place applies a step measured to raise the merit, and
    !> a globalization that keeps a step it rejected is not a globalization.
    !> The caller decides what a failure means - the nonlinear loop marks the
    !> physics diverged, the homotopy stage loop shrinks its increment instead.
    subroutine apply_line_search(self, iter_nl, reference_norm, floor_thermal, floor_hydraulic, &
                                 accepted, scale, trials, e_thermal, e_hydraulic, report_failure)
        implicit none
        !> Solver object; its iterate, residual and line-search telemetry are updated
        class(type_ftcms), intent(inout) :: self
        !> Current nonlinear iterate index, used only for reporting
        integer(int32), intent(in) :: iter_nl
        !> Merit at x_k, taken before the linear solve
        real(real64), intent(in) :: reference_norm
        !> Physical residual floors of the two blocks
        real(real64), intent(in) :: floor_thermal, floor_hydraulic
        !> True when an admissible step was found
        logical, intent(inout) :: accepted
        !> Step scale of the accepted (or last) trial
        real(real64), intent(inout) :: scale
        !> Number of trials used
        integer(int32), intent(inout) :: trials
        !> Per-block local error at the final iterate
        real(real64), intent(inout) :: e_thermal, e_hydraulic
        !> Print the backtracking table when no step is admissible
        logical, intent(in), optional :: report_failure

        real(real64), allocatable :: T_saved(:), P_saved(:)
        real(real64), allocatable :: aa_T(:), aa_P(:), aa_duT(:), aa_duP(:)
        logical, allocatable :: onset_saved(:)
        logical :: aa_has_prev, do_report
        real(real64) :: aa_gnorm, trial_norm
        real(real64) :: trial_alpha(LINE_SEARCH_MAX_TRIALS)
        real(real64) :: trial_ratio(LINE_SEARCH_MAX_TRIALS)
        real(real64) :: trial_E_T(LINE_SEARCH_MAX_TRIALS)
        real(real64) :: trial_E_H(LINE_SEARCH_MAX_TRIALS)
        integer(int32) :: trial, diag

        do_report = .true.
        if (present(report_failure)) do_report = report_failure

        call ls_snapshot(self, T_saved, P_saved, onset_saved, &
                         aa_T, aa_P, aa_duT, aa_duP, aa_has_prev, aa_gnorm)
        scale = 1.0d0
        accepted = .false.
        trial_alpha(:) = -1.0d0
        trial_ratio(:) = -1.0d0
        trial_E_T(:) = -1.0d0
        trial_E_H(:) = -1.0d0

        do trial = 1, LINE_SEARCH_MAX_TRIALS
            if (trial > 1) then
                call ls_restore(self, T_saved, P_saved, onset_saved, &
                                aa_T, aa_P, aa_duT, aa_duP, aa_has_prev, aa_gnorm)
            end if
            call self%reflect_variables(step_scale=scale)
            call self%assemble(residual_only=.true.)
            call self%apply_bc(prescribed=.false.)
            call ls_block_norms(self, e_thermal, e_hydraulic, floor_thermal, floor_hydraulic)
            trial_norm = ls_merit(e_thermal, e_hydraulic)
            if (reference_norm <= 0.0d0) then
                accepted = .true.
                exit
            end if
            trial_alpha(trial) = scale
            trial_ratio(trial) = trial_norm / reference_norm
            trial_E_T(trial) = e_thermal
            trial_E_H(trial) = e_hydraulic
            if (trial_norm <= (1.0d0 - LINE_SEARCH_ARMIJO * scale) * reference_norm) then
                accepted = .true.
                exit
            end if
            if (scale <= LINE_SEARCH_MIN_SCALE) exit
            scale = max(LINE_SEARCH_MIN_SCALE, LINE_SEARCH_BACKTRACK * scale)
        end do

        trials = trial
        self%last_line_search_scale = scale
        self%last_line_search_trials = trial
        if (accepted) return

        self%last_line_search_failures = self%last_line_search_failures + 1
        call ls_restore(self, T_saved, P_saved, onset_saved, &
                        aa_T, aa_P, aa_duT, aa_duP, aa_has_prev, aa_gnorm)
        call self%reflect_variables(step_scale=0.0d0)
        call self%assemble(residual_only=.true.)
        call self%apply_bc(prescribed=.false.)
        if (.not. do_report) return

        if (.not. branch_report_done .and. iter_nl >= 10) then
            branch_report_done = .true.
            call report_active_branches(self)
        end if
        ! Report the whole backtracking sequence, not just its last entry. The
        ! two failure mechanisms are only distinguishable from the trend in
        ! alpha:
        !   ratio - 1 proportional to alpha  -> the direction is genuinely an
        !     ascent direction for ||R||, i.e. the linearization is wrong;
        !   ratio - 1 constant in alpha      -> part of the update does not
        !     scale with alpha, i.e. the search is not controlling the step.
        write (*, '(A,I0,A,ES11.3)') &
            '   [LS] no descent at iter ', iter_nl, ': ||R_k||=', reference_norm
        do diag = 1, LINE_SEARCH_MAX_TRIALS
            if (trial_alpha(diag) < 0.0d0) cycle
            ! ES, not F: a rejected trial can raise the merit by many orders of
            ! magnitude, and an overflowing F field aborts the run on a Fortran
            ! output conversion error.
            write (*, '(A,ES9.2,A,ES12.4,A,ES11.3,A,ES11.3)') &
                '        alpha=', trial_alpha(diag), &
                '  merit/merit_k=', trial_ratio(diag), &
                '  E_T=', trial_E_T(diag), &
                '  E_H=', trial_E_H(diag)
        end do
    end subroutine apply_line_search

    !> Save everything a line-search trial mutates. The AA(1) secant pair is
    !> included because prepare_coupled_aa_step commits it unconditionally: a
    !> rejected trial must not leave its pair in the history, or the next
    !> iteration extrapolates from a step that was never taken. The per-node
    !> phase-onset flags are included for the same reason - they are once-per-step
    !> by design, so a rejected trial firing one would silently change the
    !> remaining trials.
    subroutine ls_snapshot(self, T_saved, P_saved, onset_saved, &
                           aa_T, aa_P, aa_duT, aa_duP, aa_has_prev, aa_gnorm)
        implicit none
        class(type_ftcms), intent(inout) :: self
        real(real64), allocatable, intent(inout) :: T_saved(:), P_saved(:)
        logical, allocatable, intent(inout) :: onset_saved(:)
        real(real64), allocatable, intent(inout) :: aa_T(:), aa_P(:), aa_duT(:), aa_duP(:)
        logical, intent(inout) :: aa_has_prev
        real(real64), intent(inout) :: aa_gnorm

        real(real64), pointer, contiguous :: field(:)

        nullify (field)
        if (self%is_active_thermal()) then
            call self%temperature%get_current(field)
            if (associated(field)) then
                if (allocated(T_saved)) deallocate (T_saved)
                allocate (T_saved, source=field)
            end if
            nullify (field)
        end if
        if (self%is_active_hydraulic()) then
            call self%pressure%get_current(field)
            if (associated(field)) then
                if (allocated(P_saved)) deallocate (P_saved)
                allocate (P_saved, source=field)
            end if
            nullify (field)
        end if
        if (allocated(self%phase_onset_reset)) then
            if (allocated(onset_saved)) deallocate (onset_saved)
            allocate (onset_saved, source=self%phase_onset_reset)
        end if
        call copy_alloc(aa_T, self%aa_T_prev)
        call copy_alloc(aa_P, self%aa_P_prev)
        call copy_alloc(aa_duT, self%aa_duT_prev)
        call copy_alloc(aa_duP, self%aa_duP_prev)
        aa_has_prev = self%aa_has_prev
        aa_gnorm = self%aa_gnorm_prev
    end subroutine ls_snapshot

    !> Restore the state saved by ls_snapshot so the next trial starts from x_k.
    subroutine ls_restore(self, T_saved, P_saved, onset_saved, &
                          aa_T, aa_P, aa_duT, aa_duP, aa_has_prev, aa_gnorm)
        implicit none
        class(type_ftcms), intent(inout) :: self
        real(real64), allocatable, intent(in) :: T_saved(:), P_saved(:)
        logical, allocatable, intent(in) :: onset_saved(:)
        real(real64), allocatable, intent(in) :: aa_T(:), aa_P(:), aa_duT(:), aa_duP(:)
        logical, intent(in) :: aa_has_prev
        real(real64), intent(in) :: aa_gnorm

        real(real64), pointer, contiguous :: field(:)

        nullify (field)
        if (allocated(T_saved) .and. self%is_active_thermal()) then
            call self%temperature%get_current(field)
            if (associated(field)) then
                if (size(field) == size(T_saved)) field(:) = T_saved(:)
            end if
            nullify (field)
        end if
        if (allocated(P_saved) .and. self%is_active_hydraulic()) then
            call self%pressure%get_current(field)
            if (associated(field)) then
                if (size(field) == size(P_saved)) field(:) = P_saved(:)
            end if
            nullify (field)
        end if
        if (allocated(onset_saved) .and. allocated(self%phase_onset_reset)) then
            if (size(self%phase_onset_reset) == size(onset_saved)) &
                self%phase_onset_reset(:) = onset_saved(:)
        end if
        call copy_alloc(self%aa_T_prev, aa_T)
        call copy_alloc(self%aa_P_prev, aa_P)
        call copy_alloc(self%aa_duT_prev, aa_duT)
        call copy_alloc(self%aa_duP_prev, aa_duP)
        self%aa_has_prev = aa_has_prev
        self%aa_gnorm_prev = aa_gnorm
    end subroutine ls_restore

    !> One-shot scan of the block residuals along the Newton direction.
    !>
    !> Eight backtracking points on a combined norm are not enough to say what
    !> shape the merit has. This walks alpha over four decades and reports each
    !> block separately, together with the node each block's maximum sits on and
    !> that node's thermodynamic state, so a smooth-but-ascending direction, a
    !> kink, and a flat plateau can be told apart.
    subroutine scan_residual_along_direction(self, iter_nl, ref_T, ref_H)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: iter_nl
        real(real64), intent(in) :: ref_T, ref_H

        integer(int32), parameter :: NUM_SCAN = 25
        real(real64), allocatable :: T_saved(:), P_saved(:)
        logical, allocatable :: onset_saved(:)
        real(real64), allocatable :: aa_T(:), aa_P(:), aa_duT(:), aa_duP(:)
        logical :: aa_has_prev
        real(real64) :: aa_gnorm
        real(real64), allocatable :: r_thermal(:), r_hydraulic(:), du_T(:), du_P(:)
        real(real64), pointer, contiguous :: field(:)
        real(real64) :: alpha, norm_T, norm_H, linf_T, linf_H
        integer(int32) :: k, node_T, node_H

        call ls_snapshot(self, T_saved, P_saved, onset_saved, &
                         aa_T, aa_P, aa_duT, aa_duP, aa_has_prev, aa_gnorm)

        call self%get_variable_increment(PHYSICS_TYPES%THERMAL, du_T)
        call self%get_variable_increment(PHYSICS_TYPES%HYDRAULIC, du_P)

        write (*, '(A,I0)') '   [SCAN] residual along the Newton direction at iter ', iter_nl
        write (*, '(A,ES13.5,A,ES13.5)') '   [SCAN] step-initial norms: ||R_T^0||=', ref_T, &
            ' ||R_H^0||=', ref_H
        if (allocated(du_T) .and. allocated(du_P)) then
            write (*, '(A,ES11.3,A,ES11.3)') '   [SCAN] max|du_T|=', maxval(abs(du_T)), &
                ' max|du_p|=', maxval(abs(du_P))
        end if
        write (*, '(A)') '   [SCAN]     alpha      ||R_T||/0      ||R_H||/0     maxR_T@node     maxR_H@node'

        do k = 1, NUM_SCAN
            alpha = 10.0d0**(-4.0d0 * real(k - 1, real64) / real(NUM_SCAN - 1, real64))
            call ls_restore(self, T_saved, P_saved, onset_saved, &
                            aa_T, aa_P, aa_duT, aa_duP, aa_has_prev, aa_gnorm)
            call self%reflect_variables(step_scale=alpha)
            call self%assemble()
            call self%apply_bc(prescribed=.false.)

            call self%get_variable_residual(PHYSICS_TYPES%THERMAL, r_thermal)
            call self%get_variable_residual(PHYSICS_TYPES%HYDRAULIC, r_hydraulic)
            norm_T = 0.0d0; norm_H = 0.0d0
            linf_T = 0.0d0; linf_H = 0.0d0
            node_T = 0; node_H = 0
            if (allocated(r_thermal)) then
                norm_T = sqrt(dot_product(r_thermal, r_thermal))
                node_T = maxloc(abs(r_thermal), dim=1)
                linf_T = abs(r_thermal(node_T))
            end if
            if (allocated(r_hydraulic)) then
                norm_H = sqrt(dot_product(r_hydraulic, r_hydraulic))
                node_H = maxloc(abs(r_hydraulic), dim=1)
                linf_H = abs(r_hydraulic(node_H))
            end if
            write (*, '(A,ES10.3,2(2X,ES13.5),2(2X,ES10.3,A,I0))') '   [SCAN] ', alpha, &
                norm_T / max(ref_T, tiny(1.0d0)), norm_H / max(ref_H, tiny(1.0d0)), &
                linf_T, '@', node_T, linf_H, '@', node_H
        end do

        ! Report the state of the nodes that dominate each block at the base
        ! iterate, which is what the direction has to fix.
        call ls_restore(self, T_saved, P_saved, onset_saved, &
                        aa_T, aa_P, aa_duT, aa_duP, aa_has_prev, aa_gnorm)
        call self%reflect_variables(step_scale=0.0d0)
        call self%assemble()
        call self%apply_bc(prescribed=.false.)
        call self%get_variable_residual(PHYSICS_TYPES%THERMAL, r_thermal)
        call self%get_variable_residual(PHYSICS_TYPES%HYDRAULIC, r_hydraulic)
        if (allocated(r_thermal)) node_T = maxloc(abs(r_thermal), dim=1)
        if (allocated(r_hydraulic)) node_H = maxloc(abs(r_hydraulic), dim=1)
        call report_node(self, node_T, 'max R_T', du_T, du_P)
        call report_node(self, node_H, 'max R_H', du_T, du_P)

        call ls_restore(self, T_saved, P_saved, onset_saved, &
                        aa_T, aa_P, aa_duT, aa_duP, aa_has_prev, aa_gnorm)
        call self%reflect_variables(step_scale=0.0d0)
        call self%assemble()
        call self%apply_bc(prescribed=.false.)
    end subroutine scan_residual_along_direction

    !> Verify the analytic constitutive tangent against the state function it
    !> claims to differentiate, at one node, by central differences.
    !>
    !> A large apparent heat capacity is not by itself a fault: with the
    !> Clapeyron slope and a steep retention curve, dH/dT of order 1e8 is
    !> reachable. What matters is whether that tangent is the derivative of the
    !> H the residual actually evaluates. If it is, a correct Jacobian carrying
    !> it will SHRINK the temperature update, not enlarge it - so a single
    !> iteration moving the liquid content by a third of the pore volume points
    !> at a mismatch between the assembled tangent and the state update, not at
    !> the stiffness itself.
    !>
    !> The quotients must approach the analytic value as epsilon falls and then
    !> stop improving when round-off takes over. A ratio that stays away from
    !> one at every epsilon is a wrong derivative; one that drifts with epsilon
    !> without settling is a non-differentiable state function.
    subroutine verify_local_tangent(self, node_id)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: node_id

        integer(int32), parameter :: NUM_EPSILON = 4
        integer(int32), parameter :: NUM_ALPHA = 6
        real(real64), parameter :: ALPHA_SCAN(NUM_ALPHA) = &
            [1.0d0, 1.0d-1, 1.0d-2, 1.0d-3, 1.0d-4, 1.0d-5]
        real(real64), parameter :: EPSILON_TEMPERATURE(NUM_EPSILON) = &
            [1.0d-6, 1.0d-5, 1.0d-4, 1.0d-3]
        type(type_state) :: state
        real(real64) :: temperature, pressure
        real(real64) :: enthalpy_plus, enthalpy_minus, liquid_plus, liquid_minus
        real(real64) :: ice_plus, ice_minus
        real(real64) :: analytic_dH_dT, analytic_dQw_dT, analytic_dQi_dT
        real(real64) :: fd_dH_dT, fd_dQw_dT, fd_dQi_dT
        real(real64) :: analytic_dH_dP, analytic_dQw_dP, step_T, step_P
        real(real64) :: enthalpy_base, liquid_base, predicted_H, predicted_Qw
        real(real64), allocatable :: increment_T(:), increment_P(:)
        integer(int32) :: k, row_start, repr_elem, material_id

        row_start = self%node_material_table%ptr(node_id)
        if (row_start >= self%node_material_table%ptr(node_id + 1)) return
        repr_elem = self%node_material_table%repr_element(row_start)
        call self%domain%get_material_id(repr_elem, material_id)

        call self%set_state(node_id, repr_elem, state, calc_physics=.true., include_fluxes=.false.)
        call state%temperature%get(temperature)
        call state%pressure%get(pressure)
        call state%dQw_dT%get(analytic_dQw_dT)
        call state%dQw_dP%get(analytic_dQw_dP)
        call state%dQi_dT%get(analytic_dQi_dT)
        analytic_dH_dT = 0.0d0
        call self%thermal%compute_mass_term(material_id, state, analytic_dH_dT)

        write (*, '(A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') &
            '   [TANGENT] node ', node_id, ' analytic dH/dT=', analytic_dH_dT, &
            ' dQw/dT=', analytic_dQw_dT, ' dQi/dT=', analytic_dQi_dT
        write (*, '(A)') '   [TANGENT]      eps_T        dH/dT_fd     ratio      dQw/dT_fd    ratio' // &
            '      dQi/dT_fd    ratio'

        ! Linear prediction against the actual change along the Newton
        ! direction. The tangent may be exact and still be the wrong operator
        ! for the step the solver takes; this is the test that separates the
        ! two. r -> 1 as alpha falls means the linearization is valid and the
        ! step is simply too long; r staying away from 1 means the tangent does
        ! not describe the update.
        call self%get_variable_increment(PHYSICS_TYPES%THERMAL, increment_T)
        call self%get_variable_increment(PHYSICS_TYPES%HYDRAULIC, increment_P)
        step_T = 0.0d0
        step_P = 0.0d0
        if (allocated(increment_T)) then
            if (node_id <= size(increment_T)) step_T = increment_T(node_id)
        end if
        if (allocated(increment_P)) then
            if (node_id <= size(increment_P)) step_P = increment_P(node_id)
        end if
        analytic_dH_dP = 0.0d0
        call self%thermal%compute_coupling_mass_term(material_id, state, analytic_dH_dP)
        call self%thermal%calc_enthalpy_density(material_id, state, enthalpy_base)
        call state%water_content%get(liquid_base)

        write (*, '(A,ES11.3,A,ES11.3,A,ES12.4)') '   [TANGENT] step du_T=', step_T, &
            ' du_p=', step_P, ' analytic dH/dp=', analytic_dH_dP
        write (*, '(A)') '   [TANGENT]      alpha        dH_lin       dH_act      r_H' // &
            '        dQw_lin      dQw_act     r_Qw'
        do k = 1, NUM_ALPHA
            call evaluate_step(ALPHA_SCAN(k), enthalpy_plus, liquid_plus, ice_plus)
            predicted_H = (analytic_dH_dT * step_T + analytic_dH_dP * step_P) * ALPHA_SCAN(k)
            predicted_Qw = (analytic_dQw_dT * step_T + analytic_dQw_dP * step_P) * ALPHA_SCAN(k)
            ! ES, not F: safe_ratio spans many decades at a rejected iterate
            ! and an overflowing F field aborts the run from inside a print.
            write (*, '(A,ES11.3,2(2X,ES12.4),1X,ES11.3,2(2X,ES12.4),1X,ES11.3)') '   [TANGENT] ', &
                ALPHA_SCAN(k), predicted_H, enthalpy_plus - enthalpy_base, &
                safe_ratio(enthalpy_plus - enthalpy_base, predicted_H), &
                predicted_Qw, liquid_plus - liquid_base, &
                safe_ratio(liquid_plus - liquid_base, predicted_Qw)
        end do

        do k = 1, NUM_EPSILON
            call evaluate_at(temperature + EPSILON_TEMPERATURE(k), enthalpy_plus, liquid_plus, ice_plus)
            call evaluate_at(temperature - EPSILON_TEMPERATURE(k), enthalpy_minus, liquid_minus, ice_minus)
            fd_dH_dT = (enthalpy_plus - enthalpy_minus) / (2.0d0 * EPSILON_TEMPERATURE(k))
            fd_dQw_dT = (liquid_plus - liquid_minus) / (2.0d0 * EPSILON_TEMPERATURE(k))
            fd_dQi_dT = (ice_plus - ice_minus) / (2.0d0 * EPSILON_TEMPERATURE(k))
            write (*, '(A,ES11.3,3(2X,ES12.4,1X,ES11.3))') '   [TANGENT] ', EPSILON_TEMPERATURE(k), &
                fd_dH_dT, safe_ratio(fd_dH_dT, analytic_dH_dT), &
                fd_dQw_dT, safe_ratio(fd_dQw_dT, analytic_dQw_dT), &
                fd_dQi_dT, safe_ratio(fd_dQi_dT, analytic_dQi_dT)
        end do

    contains

        subroutine evaluate_step(alpha, enthalpy, liquid, ice)
            implicit none
            real(real64), intent(in) :: alpha
            real(real64), intent(inout) :: enthalpy, liquid, ice

            type(type_state) :: trial_state

            call trial_state%copy(state)
            call trial_state%temperature%set(temperature + alpha * step_T)
            call trial_state%pressure%set(pressure + alpha * step_P)
            call self%thermal%update_water_phases(material_id, trial_state)
            call self%thermal%calc_enthalpy_density(material_id, trial_state, enthalpy)
            call trial_state%water_content%get(liquid)
            call trial_state%ice_content%get(ice)
        end subroutine evaluate_step

        subroutine evaluate_at(trial_temperature, enthalpy, liquid, ice)
            implicit none
            real(real64), intent(in) :: trial_temperature
            real(real64), intent(inout) :: enthalpy, liquid, ice

            type(type_state) :: trial_state

            call trial_state%copy(state)
            call trial_state%temperature%set(trial_temperature)
            call self%thermal%update_water_phases(material_id, trial_state)
            call self%thermal%calc_enthalpy_density(material_id, trial_state, enthalpy)
            call trial_state%water_content%get(liquid)
            call trial_state%ice_content%get(ice)
        end subroutine evaluate_at

        pure function safe_ratio(numerator, denominator) result(ratio)
            implicit none
            real(real64), intent(in) :: numerator, denominator
            real(real64) :: ratio

            ratio = 0.0d0
            if (abs(denominator) > tiny(1.0d0)) ratio = numerator / denominator
        end function safe_ratio
    end subroutine verify_local_tangent

    !> One line per nonlinear iterate.
    !>
    !> Carries exactly the quantities the acceptance decision uses: the four
    !> error measures (each passes at <= 1, -1 = not evaluated), the conserved
    !> change measure, the increment magnitudes and the relaxation factor, then
    !> the verdict naming the gates that failed. The R0-normalized residual
    !> ratios are deliberately absent - they steer omega only, and printing
    !> them beside the gates invites reading them as the criterion.
    subroutine report_nonlinear_iterate(self, iter_nl)
        implicit none
        class(type_ftcms), intent(in) :: self
        integer(int32), intent(in) :: iter_nl

        real(real64) :: loc_thermal, loc_hydraulic, bal_thermal, bal_hydraulic, dq_effective
        logical :: local_ok, balance_ok, dq_ok
        character(len=24) :: verdict

        loc_thermal = -1.0d0
        loc_hydraulic = -1.0d0
        bal_thermal = -1.0d0
        bal_hydraulic = -1.0d0
        dq_effective = -1.0d0
        local_ok = .false.
        balance_ok = .false.
        dq_ok = .false.
        call self%control%get_conserved_gates(loc_thermal, loc_hydraulic, &
                                              bal_thermal, bal_hydraulic, dq_effective, &
                                              local_ok, balance_ok, dq_ok)

        if (local_ok .and. balance_ok .and. dq_ok) then
            verdict = 'converged'
        else
            verdict = 'fail:'
            if (.not. local_ok) verdict = trim(verdict)//' loc'
            if (.not. balance_ok) verdict = trim(verdict)//' bal'
            if (.not. dq_ok) verdict = trim(verdict)//' dq'
        end if

        write (*, '(A,I3,9(A,ES10.2),2A)') '   [NL] it', iter_nl, &
            ' locT', loc_thermal, ' locH', loc_hydraulic, &
            ' balT', bal_thermal, ' balH', bal_hydraulic, &
            ' dQ', dq_effective, &
            ' duT', self%last_du_thermal_max, ' dup', self%last_du_hydraulic_max, &
            ' ls', self%last_line_search_scale, &
            ' om', self%control%get_conserved_relaxation(), &
            '  ', trim(verdict)
    end subroutine report_nonlinear_iterate

    !> State and constitutive derivatives across the freezing band.
    !>
    !> Reports, for an iteration that is not converging: how many nodes sit
    !> where the matric and freezing suctions cross (s_m = s_f is where the
    !> effective suction switches branch), the largest phase-content movement
    !> between iterations, and at the node that moves most the suctions and the
    !> three derivatives the linearization depends on. A non-smooth local update
    !> shows up as a large phase movement concentrated on the band with
    !> derivatives that jump between iterations.
    subroutine report_freezing_band(self, iter_nl)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: iter_nl

        real(real64), parameter :: CLAPEYRON_SLOPE = rho_std * Lf0 / (Tf0 + TtoK)
        real(real64), parameter :: BAND_WIDTH_PA = 1.0d4
        type(type_state) :: state
        real(real64), pointer, contiguous :: liquid(:), ice(:)
        real(real64) :: temperature, pressure, suction_matric, suction_freezing, suction_effective
        real(real64) :: change_liquid, change_ice, max_change_liquid, max_change_ice
        real(real64) :: dliquid_dP, dice_dP, dliquid_dT, dice_dT, heat_capacity
        integer(int32) :: node_id, num_nodes, band_count, worst_node
        integer(int32) :: row_start, repr_elem, material_id

        nullify (liquid, ice)
        call self%domain%get_num_nodes(num_nodes)
        call self%Qw%get_current(liquid)
        call self%Qi%get_current(ice)
        if (.not. (associated(liquid) .and. associated(ice))) return

        if (.not. allocated(freeze_prev_liquid)) allocate (freeze_prev_liquid(num_nodes), source=0.0d0)
        if (.not. allocated(freeze_prev_ice)) allocate (freeze_prev_ice(num_nodes), source=0.0d0)

        band_count = 0
        max_change_liquid = 0.0d0
        max_change_ice = 0.0d0
        worst_node = 0
        do node_id = 1, min(num_nodes, size(liquid))
            call self%temperature%get_current(node_id, temperature)
            ! A thermal-only run never allocates the pressure field.
            pressure = 0.0d0
            if (self%is_active_hydraulic()) call self%pressure%get_current(node_id, pressure)
            suction_matric = -pressure
            suction_freezing = 0.0d0
            if (temperature < Tf0) suction_freezing = CLAPEYRON_SLOPE * (Tf0 - temperature)
            if (abs(suction_matric - suction_freezing) < BAND_WIDTH_PA) band_count = band_count + 1

            change_liquid = abs(liquid(node_id) - freeze_prev_liquid(node_id))
            change_ice = abs(ice(node_id) - freeze_prev_ice(node_id))
            max_change_liquid = max(max_change_liquid, change_liquid)
            if (change_ice > max_change_ice) then
                max_change_ice = change_ice
                worst_node = node_id
            end if
        end do

        write (*, '(A,I0,A,I0,A,ES11.3,A,ES11.3)') '   [FREEZE] iter ', iter_nl, &
            '  band|s_m-s_f|<1e4 Pa: ', band_count, '  max|dQw|=', max_change_liquid, &
            '  max|dQi|=', max_change_ice

        if (worst_node >= 1) then
            row_start = self%node_material_table%ptr(worst_node)
            if (row_start < self%node_material_table%ptr(worst_node + 1)) then
                repr_elem = self%node_material_table%repr_element(row_start)
                call self%domain%get_material_id(repr_elem, material_id)
                call self%set_state(worst_node, repr_elem, state, calc_physics=.true., include_fluxes=.false.)
                call self%temperature%get_current(worst_node, temperature)
                ! A thermal-only run never allocates the pressure field.
                pressure = 0.0d0
                if (self%is_active_hydraulic()) call self%pressure%get_current(worst_node, pressure)
                suction_matric = -pressure
                suction_freezing = 0.0d0
                if (temperature < Tf0) suction_freezing = CLAPEYRON_SLOPE * (Tf0 - temperature)
                suction_effective = 0.0d0
                call state%effective_suction%get(suction_effective)
                call state%dQw_dP%get(dliquid_dP)
                call state%dQi_dP%get(dice_dP)
                call state%dQw_dT%get(dliquid_dT)
                call state%dQi_dT%get(dice_dT)
                heat_capacity = 0.0d0
                call self%thermal%compute_mass_term(material_id, state, heat_capacity)
                write (*, '(A,I0,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3)') &
                    '   [FREEZE]   worst dQi node=', worst_node, ' T=', temperature, &
                    ' s_m=', suction_matric, ' s_f=', suction_freezing, ' s_eff=', suction_effective, &
                    ' p=', pressure
                write (*, '(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
                    '   [FREEZE]   dQw/dp=', dliquid_dP, ' dQi/dp=', dice_dP, &
                    ' dQw/dT=', dliquid_dT, ' dQi/dT=', dice_dT, ' dH/dT=', heat_capacity
                if (freeze_report_count <= 2) call verify_local_tangent(self, worst_node)
            end if
        end if

        do node_id = 1, min(num_nodes, size(liquid))
            freeze_prev_liquid(node_id) = liquid(node_id)
            freeze_prev_ice(node_id) = ice(node_id)
        end do
        nullify (liquid, ice)
    end subroutine report_freezing_band

    !> Report which constitutive branch the dominant residual nodes sit on.
    !>
    !> The element tangent matches a finite difference and the linear solve is
    !> exact, yet no step length reduces the merit. For an exact two-sided
    !> Jacobian that cannot happen, so the residual must only be piecewise
    !> differentiable there and the iterate must be sitting on a kink. The
    !> constitutive branches that can produce one are the pore-volume bound on
    !> the total water, the ice-content clamps at 0 and at the porosity, the
    !> saturation of the impedance ratio, and the phase-volume closure. All of
    !> them are decidable from the published nodal phase fields.
    subroutine report_active_branches(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        integer(int32), parameter :: NUM_REPORTED = 5
        real(real64), allocatable :: r_thermal(:), r_hydraulic(:)
        real(real64), pointer, contiguous :: field(:)
        real(real64) :: temperature, pressure, ice, water, porosity
        real(real64) :: theta_total, volume_bound, ratio_denominator, impedance_ratio
        real(real64) :: density_ratio
        integer(int32) :: k, node_id
        integer(int32) :: ranked(NUM_REPORTED)

        density_ratio = 917.0d0 / 1000.0d0
        call self%get_variable_residual(PHYSICS_TYPES%THERMAL, r_thermal)
        call self%get_variable_residual(PHYSICS_TYPES%HYDRAULIC, r_hydraulic)

        write (*, '(A)') '   [BRANCH] nodes carrying the largest residuals'
        write (*, '(A)') '   [BRANCH]  block  node        T          p        Qw       Qi      phi' // &
            '   theta_tot     bound    margin   Q_imp'

        call rank_nodes(r_thermal, ranked)
        do k = 1, NUM_REPORTED
            node_id = ranked(k)
            if (node_id < 1) cycle
            call emit('R_T', node_id)
        end do
        call rank_nodes(r_hydraulic, ranked)
        do k = 1, NUM_REPORTED
            node_id = ranked(k)
            if (node_id < 1) cycle
            call emit('R_H', node_id)
        end do

        if (allocated(r_thermal)) deallocate (r_thermal)
        if (allocated(r_hydraulic)) deallocate (r_hydraulic)

    contains

        subroutine rank_nodes(residual, ranking)
            implicit none
            real(real64), allocatable, intent(in) :: residual(:)
            integer(int32), intent(inout) :: ranking(:)

            real(real64), allocatable :: work(:)
            integer(int32) :: j

            ranking(:) = 0
            if (.not. allocated(residual)) return
            allocate (work, source=abs(residual))
            do j = 1, size(ranking)
                if (size(work) == 0) exit
                ranking(j) = maxloc(work, dim=1)
                work(ranking(j)) = -1.0d0
            end do
        end subroutine rank_nodes

        subroutine emit(label, node)
            implicit none
            character(len=*), intent(in) :: label
            integer(int32), intent(in) :: node

            call self%temperature%get_current(node, temperature)
            call self%pressure%get_current(node, pressure)
            call self%Qw%get_current(node, water)
            call self%Qi%get_current(node, ice)
            call self%porosity%get_current(node, porosity)

            theta_total = water + density_ratio * ice
            volume_bound = porosity * density_ratio + water * (1.0d0 - density_ratio)
            ratio_denominator = water + ice
            impedance_ratio = 0.0d0
            if (ratio_denominator > tiny(1.0d0)) impedance_ratio = ice / ratio_denominator

            write (*, '(A,A,2X,I6,F10.4,ES12.3,5F10.5,F9.4)') '   [BRANCH]  ', label, node, &
                temperature, pressure, water, ice, porosity, theta_total, volume_bound, &
                theta_total - volume_bound, impedance_ratio
        end subroutine emit
    end subroutine report_active_branches

    !> Print the primary and phase state of one node.
    subroutine report_node(self, node_id, label, du_T, du_P)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: node_id
        character(len=*), intent(in) :: label
        real(real64), allocatable, intent(in) :: du_T(:), du_P(:)

        real(real64) :: temperature, pressure, ice, water, porosity
        real(real64) :: step_T, step_P

        if (node_id < 1) return
        call self%temperature%get_current(node_id, temperature)
        ! A thermal-only run never allocates the pressure field.
        pressure = 0.0d0
        if (self%is_active_hydraulic()) call self%pressure%get_current(node_id, pressure)
        call self%Qi%get_current(node_id, ice)
        call self%Qw%get_current(node_id, water)
        call self%porosity%get_current(node_id, porosity)
        step_T = 0.0d0
        step_P = 0.0d0
        if (allocated(du_T)) then
            if (node_id <= size(du_T)) step_T = du_T(node_id)
        end if
        if (allocated(du_P)) then
            if (node_id <= size(du_P)) step_P = du_P(node_id)
        end if
        write (*, '(A,A,A,I0,A,F10.5,A,ES12.4,A,F8.5,A,F8.5,A,F8.5,A,ES11.3,A,ES11.3)') &
            '   [SCAN] ', label, ' node=', node_id, ' T=', temperature, ' p=', pressure, &
            ' Qi=', ice, ' Qw=', water, ' phi=', porosity, ' du_T=', step_T, ' du_p=', step_P
    end subroutine report_node

    !> Copy an allocatable array, propagating the unallocated state.
    subroutine copy_alloc(dst, src)
        implicit none
        real(real64), allocatable, intent(inout) :: dst(:)
        real(real64), allocatable, intent(in) :: src(:)

        if (allocated(dst)) deallocate (dst)
        if (allocated(src)) allocate (dst, source=src)
    end subroutine copy_alloc


    module subroutine solve_time_step_staggered_ftcms(self, is_step_converged)
        implicit none
        class(type_ftcms), intent(inout) :: self
        logical, intent(inout) :: is_step_converged

        logical :: prescribe_bc
        integer(int32) :: iter_nl, coupling_iter, num_nodes, bdf_order
        integer(int32), parameter :: MAX_COUPLING_ITER = 10
        integer(int32), parameter :: MAX_PHASE_NL_ITER = 100
        real(real64), parameter :: COUPLING_TOL = 1.0d-3
        real(real64), parameter :: THERMAL_INCREMENT_GUARD = 1.0d6
        real(real64), parameter :: HYDRAULIC_INCREMENT_GUARD = 1.0d8
        real(real64) :: t_res, t_inc, h_res, h_inc
        real(real64) :: coupling_change_T, coupling_change_P, T_scale, P_scale
        real(real64) :: mean_pressure
        real(real64) :: phase_inc_max
        real(real64), allocatable :: T_old(:), P_old(:)
        real(real64), allocatable :: phase_increment(:)
        real(real64), allocatable :: Qw_save(:), dW_check(:)
        real(real64), allocatable :: hyd_residual_local(:)
        real(real64), pointer, contiguous :: T_cur(:) => null()
        real(real64), pointer, contiguous :: P_cur(:) => null()
        real(real64), pointer, contiguous :: Qw_cur(:) => null()
        logical :: linear_failed
        logical :: excessive_update
        character(len=16) :: phase_label

        is_step_converged = .false.

        call self%domain%get_num_nodes(num_nodes)
        allocate (T_old(num_nodes), P_old(num_nodes))

        call self%solve_time_step_initial_setup()
        self%last_phase_iterations = 1
        self%last_inner_iterations = 0
        self%last_nonlinear_work = 0
        self%last_solve_status = SOLVE_STATUS_NOT_RUN
        self%last_phase_metrics_available = .false.
        self%last_phase_increment_max = -1.0d0
        ! Attempt-scoped: a linear failure or an early divergence return can end
        ! the attempt before any iterate is reported, and the history record
        ! would otherwise carry the previous attempt's increments.
        self%last_du_thermal_max = -1.0d0
        self%last_du_hydraulic_max = -1.0d0

        coupling_loop: do coupling_iter = 1, MAX_COUPLING_ITER
            self%last_phase_iterations = coupling_iter

            if (coupling_iter > 1) then
                call self%temperature%get_current(T_cur)
                call self%pressure%get_current(P_cur)
                if (associated(T_cur)) T_old(:) = T_cur(:)
                if (associated(P_cur)) P_old(:) = P_cur(:)
                nullify (T_cur)
                nullify (P_cur)

                call self%control%reset_iteration()
                call self%control%set_nonlinear_solver(NONLINEAR_SOLVER%PICARD)
                call self%control%increment_total()
                call self%control%reset_acceleration()
            end if

            ! =============================================================
            ! Phase 1: Hydraulic nonlinear loop (T frozen)
            ! =============================================================
            phase_label = '[HYD_NL]'
            linear_failed = .false.

            if (self%is_active_hydraulic()) then
                hydraulic_nl: do while (self%control%should_continue())
                    call self%solve_time_step_setup(prescribe_bc)
                    if (prescribe_bc) then
                        call self%prescribe_dirichlet()
                        call self%calc_gradient_temperature()
                        call self%calc_gradient_pressure()
                    end if
                    call self%assemble()
                    call self%apply_bc(prescribed=.false.)
                    call self%freeze_physics_dofs(PHYSICS_TYPES%THERMAL)
                    call self%control%get_nonlinear_iter(iter_nl)
                    self%current_physics_id = PHYSICS_TYPES%HYDRAULIC%ID
                    call self%solve()
                    call self%zero_frozen_increment(PHYSICS_TYPES%THERMAL)

                    if (.not. self%solver%is_success()) then
                        linear_failed = .true.
                        self%last_solve_status = SOLVE_STATUS_LINEAR_FAILURE
                        call self%control%set_converged(PHYSICS_TYPES%HYDRAULIC, .false.)
                        call self%control%set_diverged(PHYSICS_TYPES%HYDRAULIC, .true.)
                        exit hydraulic_nl
                    end if

                    excessive_update = .false.
                    phase_inc_max = 0.0d0
                    if (allocated(phase_increment)) deallocate (phase_increment)
                    call self%get_variable_increment(PHYSICS_TYPES%HYDRAULIC, phase_increment)
                    if (allocated(phase_increment)) then
                        if (size(phase_increment) > 0) then
                            phase_inc_max = maxval(abs(phase_increment))
                            excessive_update = phase_inc_max > HYDRAULIC_INCREMENT_GUARD
                        end if
                        deallocate (phase_increment)
                    end if
                    if (excessive_update) then
                        write (*, '(A,ES13.5,A,ES13.5,A)') '   [HYD_NL] excessive hydraulic increment detected (> ', &
                            HYDRAULIC_INCREMENT_GUARD, ', max=', phase_inc_max, '). Continue with damped update.'
                    end if

                    ! Save Qw before update to compute dW for convergence
                    call self%Qw%get_current(Qw_cur)
                    if (.not. allocated(Qw_save)) allocate (Qw_save(num_nodes))
                    if (associated(Qw_cur)) then
                        Qw_save(:) = Qw_cur(:)
                    else
                        Qw_save(:) = 0.0d0
                    end if
                    nullify (Qw_cur)

                    ! Apply update (update_nodal_phases recomputes Qw inside)
                    call self%reflect_variables()

                    ! Compute dW = Qw_new - Qw_old
                    call self%Qw%get_current(Qw_cur)
                    if (.not. allocated(dW_check)) allocate (dW_check(num_nodes))
                    if (associated(Qw_cur)) then
                        dW_check(:) = Qw_cur(:) - Qw_save(:)
                    else
                        dW_check(:) = 0.0d0
                    end if
                    nullify (Qw_cur)

                    ! Convergence check: residual + dW (water-content-based update norm)
                    call self%get_variable_residual(PHYSICS_TYPES%HYDRAULIC, hyd_residual_local)
                    if (.not. allocated(hyd_residual_local) .or. size(hyd_residual_local) == 0) then
                        call self%control%set_diverged(PHYSICS_TYPES%HYDRAULIC, .true.)
                    else if (has_nan(hyd_residual_local) .or. has_nan(dW_check)) then
                        call self%control%set_diverged(PHYSICS_TYPES%HYDRAULIC, .true.)
                    else
                        call self%control%check_convergence(PHYSICS_TYPES%HYDRAULIC, hyd_residual_local, dW_check)
                    end if

                    call self%control%get_nonlinear_iter(iter_nl)
                    if ((.not. self%control%is_converged()) .and. iter_nl >= MAX_PHASE_NL_ITER) then
                        self%last_solve_status = SOLVE_STATUS_NONLINEAR_LIMIT
                        write (*, '(A,I0,A)') '   [HYD_NL] reached nonlinear iteration cap (', MAX_PHASE_NL_ITER, &
                            '). Triggering timestep retry.'
                        call self%control%set_converged(PHYSICS_TYPES%HYDRAULIC, .false.)
                        call self%control%set_diverged(PHYSICS_TYPES%HYDRAULIC, .true.)
                        exit hydraulic_nl
                    end if

                    if (self%control%is_none()) exit hydraulic_nl
                end do hydraulic_nl

                call self%control%get_nonlinear_iter(iter_nl)
                self%last_inner_iterations = iter_nl
                self%last_nonlinear_work = self%last_nonlinear_work + max(1_int32, iter_nl)

                if (.not. self%control%is_converged()) then
                    call self%control%get_nonlinear_iter(iter_nl)
                    h_res = 0.0d0
                    h_inc = 0.0d0
                    if (.not. linear_failed) then
                        call self%control%get_current_norm(PHYSICS_TYPES%HYDRAULIC, &
                                                           NONLINEAR_NORM_CRITERIA%RESIDUAL, NORM_TYPES%LINF, h_res)
                        call self%control%get_current_norm(PHYSICS_TYPES%HYDRAULIC, &
                                                           NONLINEAR_NORM_CRITERIA%UPDATE, NORM_TYPES%LINF, h_inc)
                    end if
                    if (linear_failed) then
                        write (*, '(A,A,A,I0,A,L1,A)') '   ', phase_label, &
                            ' failed: iter=', iter_nl, ', diverged=', self%control%is_diverged(), &
                            ', linear solver failure.'
                    else
                        write (*, '(A,A,A,I0,A,L1,A,2(ES11.3,1X))') '   ', phase_label, &
                            ' failed: iter=', iter_nl, ', diverged=', self%control%is_diverged(), &
                            ', H_res/H_inc=', h_res, h_inc
                    end if
                    if (self%last_solve_status == SOLVE_STATUS_NOT_RUN) then
                        if (self%control%is_diverged()) then
                            self%last_solve_status = SOLVE_STATUS_NONLINEAR_DIVERGED
                        else
                            self%last_solve_status = SOLVE_STATUS_NONLINEAR_LIMIT
                        end if
                    end if
                    exit coupling_loop
                end if
            else
                call self%control%set_converged(PHYSICS_TYPES%HYDRAULIC, .true.)
                call self%control%set_diverged(PHYSICS_TYPES%HYDRAULIC, .false.)
            end if

            ! =============================================================
            ! Phase 2: Thermal nonlinear loop (P frozen)
            ! =============================================================
            phase_label = '[THM_NL]'

            call self%control%reset_iteration()
            call self%control%set_nonlinear_solver(NONLINEAR_SOLVER%PICARD)
            call self%control%increment_total()
            call self%control%reset_acceleration()

            linear_failed = .false.

            if (self%is_active_thermal()) then
                thermal_nl: do while (self%control%should_continue())
                    call self%solve_time_step_setup(prescribe_bc)
                    if (prescribe_bc) then
                        call self%prescribe_dirichlet()
                        call self%calc_gradient_temperature()
                        call self%calc_gradient_pressure()
                    end if
                    call self%assemble()
                    call self%apply_bc(prescribed=.false.)
                    call self%freeze_physics_dofs(PHYSICS_TYPES%HYDRAULIC)
                    call self%control%get_nonlinear_iter(iter_nl)
                    self%current_physics_id = PHYSICS_TYPES%THERMAL%ID
                    call self%solve()
                    call self%zero_frozen_increment(PHYSICS_TYPES%HYDRAULIC)

                    if (allocated(self%solver_thermal)) then
                        if (.not. self%solver_thermal%is_success()) then
                            linear_failed = .true.
                            self%last_solve_status = SOLVE_STATUS_LINEAR_FAILURE
                            call self%control%set_converged(PHYSICS_TYPES%THERMAL, .false.)
                            call self%control%set_diverged(PHYSICS_TYPES%THERMAL, .true.)
                            exit thermal_nl
                        end if
                    else if (.not. self%solver%is_success()) then
                        linear_failed = .true.
                        self%last_solve_status = SOLVE_STATUS_LINEAR_FAILURE
                        call self%control%set_converged(PHYSICS_TYPES%THERMAL, .false.)
                        call self%control%set_diverged(PHYSICS_TYPES%THERMAL, .true.)
                        exit thermal_nl
                    end if

                    excessive_update = .false.
                    phase_inc_max = 0.0d0
                    if (allocated(phase_increment)) deallocate (phase_increment)
                    call self%get_variable_increment(PHYSICS_TYPES%THERMAL, phase_increment)
                    if (allocated(phase_increment)) then
                        if (size(phase_increment) > 0) then
                            phase_inc_max = maxval(abs(phase_increment))
                            excessive_update = phase_inc_max > THERMAL_INCREMENT_GUARD
                        end if
                        deallocate (phase_increment)
                    end if
                    if (excessive_update) then
                        write (*, '(A,ES13.5,A,ES13.5,A)') '   [THM_NL] excessive thermal increment detected (> ', &
                            THERMAL_INCREMENT_GUARD, ', max=', phase_inc_max, '). Continue with damped update.'
                    end if

                    call self%solve_time_step_check_convergence(PHYSICS_TYPES%THERMAL)
                    call self%reflect_variables()

                    call self%control%get_nonlinear_iter(iter_nl)
                    if ((.not. self%control%is_converged()) .and. iter_nl >= MAX_PHASE_NL_ITER) then
                        self%last_solve_status = SOLVE_STATUS_NONLINEAR_LIMIT
                        write (*, '(A,I0,A)') '   [THM_NL] reached nonlinear iteration cap (', MAX_PHASE_NL_ITER, &
                            '). Triggering timestep retry.'
                        call self%control%set_converged(PHYSICS_TYPES%THERMAL, .false.)
                        call self%control%set_diverged(PHYSICS_TYPES%THERMAL, .true.)
                        exit thermal_nl
                    end if

                    if (self%control%is_none()) exit thermal_nl
                end do thermal_nl

                call self%control%get_nonlinear_iter(iter_nl)
                self%last_inner_iterations = iter_nl
                self%last_nonlinear_work = self%last_nonlinear_work + max(1_int32, iter_nl)

                is_step_converged = self%control%is_converged()

                if (.not. is_step_converged) then
                    call self%control%get_nonlinear_iter(iter_nl)
                    t_res = 0.0d0
                    t_inc = 0.0d0
                    if (.not. linear_failed) then
                        call self%control%get_current_norm(PHYSICS_TYPES%THERMAL, &
                                                           NONLINEAR_NORM_CRITERIA%RESIDUAL, NORM_TYPES%LINF, t_res)
                        call self%control%get_current_norm(PHYSICS_TYPES%THERMAL, &
                                                           NONLINEAR_NORM_CRITERIA%UPDATE, NORM_TYPES%LINF, t_inc)
                    end if
                    if (linear_failed) then
                        write (*, '(A,A,A,I0,A,L1,A)') '   ', phase_label, &
                            ' failed: iter=', iter_nl, ', diverged=', self%control%is_diverged(), &
                            ', linear solver failure.'
                    else
                        write (*, '(A,A,A,I0,A,L1,A,2(ES11.3,1X))') '   ', phase_label, &
                            ' failed: iter=', iter_nl, ', diverged=', self%control%is_diverged(), &
                            ', T_res/T_inc=', t_res, t_inc
                    end if
                    if (self%last_solve_status == SOLVE_STATUS_NOT_RUN) then
                        if (self%control%is_diverged()) then
                            self%last_solve_status = SOLVE_STATUS_NONLINEAR_DIVERGED
                        else
                            self%last_solve_status = SOLVE_STATUS_NONLINEAR_LIMIT
                        end if
                    end if
                    exit coupling_loop
                end if
            else
                call self%control%set_converged(PHYSICS_TYPES%THERMAL, .true.)
                call self%control%set_diverged(PHYSICS_TYPES%THERMAL, .false.)
                is_step_converged = .true.
            end if

            if (coupling_iter == 1) cycle coupling_loop

            coupling_change_T = 0.0d0
            coupling_change_P = 0.0d0

            call self%temperature%get_current(T_cur)
            call self%pressure%get_current(P_cur)

            if (associated(T_cur)) then
                T_scale = maxval(abs(T_cur)) + 1.0d0
                coupling_change_T = maxval(abs(T_cur - T_old)) / T_scale
            end if
            if (associated(P_cur)) then
                P_scale = maxval(abs(P_cur)) + 1.0d0
                coupling_change_P = maxval(abs(P_cur - P_old)) / P_scale
            end if

            nullify (T_cur)
            nullify (P_cur)

            write (*, '("   [Coupling] Iter:", I2, " dT_rel:", ES10.3, " dP_rel:", ES10.3)') &
                coupling_iter, coupling_change_T, coupling_change_P

            if (coupling_change_T < COUPLING_TOL .and. coupling_change_P < COUPLING_TOL) then
                exit coupling_loop
            end if

        end do coupling_loop

        if (allocated(T_old)) deallocate (T_old)
        if (allocated(P_old)) deallocate (P_old)
        if (allocated(phase_increment)) deallocate (phase_increment)

        if (is_step_converged) then
            self%last_solve_status = SOLVE_STATUS_CONVERGED
        else if (self%last_solve_status == SOLVE_STATUS_NOT_RUN) then
            if (self%control%is_diverged()) then
                self%last_solve_status = SOLVE_STATUS_NONLINEAR_DIVERGED
            else
                self%last_solve_status = SOLVE_STATUS_NONLINEAR_LIMIT
            end if
        end if

    end subroutine solve_time_step_staggered_ftcms

    pure function solve_status_label(status) result(label)
        implicit none
        integer(int32), intent(in) :: status
        character(len=18) :: label

        select case (status)
        case (SOLVE_STATUS_CONVERGED)
            label = "accepted"
        case (SOLVE_STATUS_LINEAR_FAILURE)
            label = "linear_failure"
        case (SOLVE_STATUS_NONLINEAR_DIVERGED)
            label = "nonlinear_diverged"
        case (SOLVE_STATUS_NONLINEAR_LIMIT)
            label = "nonlinear_limit"
        case (SOLVE_STATUS_PHASE_FAILURE)
            label = "phase_failure"
        case (SOLVE_STATUS_LTE_REJECTED)
            label = "lte_rejected"
        case default
            label = "not_run"
        end select
    end function solve_status_label

    !> Write diagnostics after time-control update so dt_next and accepted time
    !> describe the actual state used by the following attempt.
    subroutine write_solver_history_attempt(self, attempt, accepted_step, time_start, time_trial, dt_used, &
                                            accepted, ats_iter, lte_error)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: attempt, accepted_step, ats_iter
        real(real64), intent(in) :: time_start, time_trial, dt_used, lte_error
        logical, intent(in) :: accepted

        character(len=18) :: status
        real(real64) :: time_accepted, dt_next
        real(real64) :: omega_used
        real(real64) :: loc_thermal, loc_hydraulic, bal_thermal, bal_hydraulic, dq_effective
        logical :: local_ok, balance_ok, dq_ok

        if (self%solver_history_unit == -1) return

        call self%control%get_time(time_accepted)
        call self%control%get_dt(dt_next)
        omega_used = self%control%get_conserved_relaxation()
        status = solve_status_label(self%last_solve_status)

        loc_thermal = -1.0d0
        loc_hydraulic = -1.0d0
        bal_thermal = -1.0d0
        bal_hydraulic = -1.0d0
        dq_effective = -1.0d0
        local_ok = .false.
        balance_ok = .false.
        dq_ok = .false.
        call self%control%get_conserved_gates(loc_thermal, loc_hydraulic, &
                                              bal_thermal, bal_hydraulic, dq_effective, &
                                              local_ok, balance_ok, dq_ok)

        write (self%solver_history_unit, &
               '(2(I10,1X),5(ES15.7,1X),I1,1X,A18,1X,2(I10,1X),9(ES13.5,1X))') &
            attempt, accepted_step, time_start, time_trial, time_accepted, dt_used, dt_next, &
            merge(1_int32, 0_int32, accepted), status, &
            self%last_inner_iterations, ats_iter, &
            loc_thermal, loc_hydraulic, bal_thermal, bal_hydraulic, dq_effective, &
            self%last_du_thermal_max, self%last_du_hydraulic_max, omega_used, lte_error
        flush (self%solver_history_unit)
    end subroutine write_solver_history_attempt

    module subroutine run_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        logical :: is_step_converged
        integer(int32) :: consecutive_failures
        integer(int32) :: step_counter
        integer(int32) :: attempt_counter
        integer(int32) :: nl_iter
        integer(int32) :: phase_iter
        integer(int32) :: nonlinear_work
        integer(int32) :: effective_iter
        logical :: failed_at_min_dt
        real(real64) :: time_s, time_start_s, time_trial_s, dt_s, dt_used
        real(real64) :: lte_error
        !> BDF order lte_error was measured at; the retry is resized at it
        integer(int32) :: lte_order
        integer(int32), parameter :: MAX_CONSECUTIVE_FAILURES = 50

        consecutive_failures = 0
        step_counter = 0
        attempt_counter = 0

        ! Loop until end time
        time_loop: do while (.not. self%control%is_end_time())
            call self%control%get_time(time_start_s)
            call self%run_assimilation(time_start_s, 1.0d0 + time_start_s / 86400.0d0)
            call self%solve_time_step(is_step_converged)
            call self%control%get_dt(dt_used)
            time_trial_s = time_start_s + dt_used
            nl_iter = self%last_inner_iterations
            phase_iter = self%last_phase_iterations
            nonlinear_work = self%last_nonlinear_work
            ! Outer phase work is diagnosed separately. LTE controls temporal
            ! accuracy and the final monolithic inner count supplies the existing
            ! nonlinear robustness brake; reducing dt did not reduce outer work.
            effective_iter = nl_iter

            ! Local-truncation-error estimate for error-controlled ATS, evaluated
            ! before the time/variable history is shifted (needs ydot_n and dt_n).
            lte_error = -1.0d0
            lte_order = 1
            if (is_step_converged) lte_error = self%compute_lte_error(order_used=lte_order)
            if (is_step_converged .and. self%lte_error_control_active .and. lte_error > 1.0d0) then
                is_step_converged = .false.
                self%last_solve_status = SOLVE_STATUS_LTE_REJECTED
                write (*, '(A,ES11.3,A,ES11.3)') &
                    '   [LTE] converged step rejected: normalized error=', lte_error, ', dt[s]=', dt_used
            end if

            attempt_counter = attempt_counter + 1

            if (is_step_converged) then
                self%last_accepted_dt = dt_used
                ! Genuine step acceptance: is_step_converged has already
                ! survived the LTE rejection check above, so the conserved-
                ! quantity drift measured at the nonlinear loop's last
                ! iterate (check_conserved_convergence_control, criterion 2)
                ! is the one that actually becomes part of the accepted
                ! history. Fold it into the cumulative budget here, not at
                ! every nonlinear iteration, since an iterate that satisfies
                ! the nonlinear gate can still be discarded by LTE or the
                ! outer phase loop before reaching this point.
                call self%control%commit_conserved_drift()
                call self%commit_lte_history()
                ! Update time and adaptive time stepping
                call self%control%update(is_step_converged, error_estimate=lte_error, &
                                         iteration_count=effective_iter, error_order=lte_order)

                consecutive_failures = 0
                step_counter = step_counter + 1
                call write_solver_history_attempt(self, attempt_counter, step_counter, time_start_s, &
                                                  time_trial_s, dt_used, is_step_converged, effective_iter, &
                                                  lte_error)
                call self%control%get_time(time_s)
                if (step_counter == 1 .or. mod(step_counter, 20) == 0 .or. effective_iter > 8) then
                    write (*, '(A,I0,A,ES13.5,A,I0,A,I0,A,I0)') '   [STEP] converged: n=', step_counter, &
                        ', t[s]=', time_s, ', nonlinear_iter=', nl_iter, &
                        ', outer_iter=', phase_iter, ', nonlinear_work=', nonlinear_work
                end if

                ! Shift variable history on convergence
                call self%shift()

                call self%update_variables()
                call self%output_fields()
                call self%output_history()
            else
                failed_at_min_dt = self%control%is_min_dt()
                ! Update time and adaptive time stepping
                call self%control%update(is_step_converged, error_estimate=lte_error, &
                                         iteration_count=effective_iter, error_order=lte_order)

                ! Retry with smaller dt
                consecutive_failures = consecutive_failures + 1
                call write_solver_history_attempt(self, attempt_counter, step_counter, time_start_s, &
                                                  time_trial_s, dt_used, is_step_converged, effective_iter, &
                                                  lte_error)

                if (failed_at_min_dt) then
                    call self%control%get_dt(dt_s)
                    write (*, '(A,ES13.5,A)') '   [ERROR] Step failed at minimum dt=', dt_s, '. Stopping retry loop.'
                    exit time_loop
                end if

                write (*, '("   [WARNING] Step Failed (",I0,"/",I0,"). Retrying with smaller dt...")') &
                    consecutive_failures, MAX_CONSECUTIVE_FAILURES
                if (consecutive_failures >= MAX_CONSECUTIVE_FAILURES) then
                    write (*, '("   [ERROR] Too many consecutive failures. Stopping.")')
                    exit time_loop
                end if
                cycle time_loop
            end if

        end do time_loop

    end subroutine run_ftcms

    !> Build the nodal control-volume array V_i = sum_e int psi_i dOmega used
    !> by the convergence control's exact per-node local conservation gate
    !> (check_conserved_convergence_control, criterion 1). Computed once per
    !> run - the mesh is static - from the FE row-sum-lumped mass matrix with
    !> a unit coefficient: the same primitive governing_base uses to build a
    !> lumped physical capacity matrix (compute_K1_lumped), just with A_gp=1
    !> so the result is a pure geometric volume.
    !>
    !> No other nodal-volume provider exists in src/domain or src/numerical:
    !> the node_material_table built in ftcms_base.F90 sums whole adjacent-
    !> element measures per node to weight material averages, which is not
    !> integrated against shape functions and double-counts the volume shared
    !> between elements, so it is not a control volume.
    !>
    !> Cost: O(N_elem) FE evaluations, paid once.
    subroutine build_nodal_volume(self, num_elements)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: num_elements

        class(abst_fe), pointer :: fe
        integer(int32), pointer, contiguous, dimension(:) :: connectivity
        real(real64), allocatable :: coordinates(:, :)
        real(real64), allocatable :: unit_gp(:)
        real(real64), allocatable :: elem_volume(:, :)
        integer(int32) :: num_nodes_total, num_nodes_local, num_gauss_local
        integer(int32) :: i_elem, i_node

        nullify (fe)
        nullify (connectivity)

        call self%domain%get_num_nodes(num_nodes_total)
        call allocate_array(self%nodal_volume, num_nodes_total)
        self%nodal_volume = 0.0d0
        if (num_nodes_total <= 0) return

        do i_elem = 1, num_elements
            call self%domain%get_fe(i_elem, fe)
            if (.not. associated(fe)) cycle
            call self%domain%get_fe_connectivity(i_elem, connectivity)
            if (.not. associated(connectivity)) cycle
            call self%domain%get_fe_coordinate(i_elem, coordinates)

            call fe%get_num_nodes(num_nodes_local)
            call fe%get_num_gauss(num_gauss_local)
            call allocate_array(unit_gp, num_gauss_local)
            unit_gp = 1.0d0
            call allocate_array(elem_volume, num_nodes_local, num_nodes_local)

            call fe%compute_K1_lumped(coordinates, unit_gp, elem_volume)
            do i_node = 1, num_nodes_local
                self%nodal_volume(connectivity(i_node)) = self%nodal_volume(connectivity(i_node)) &
                                                          + elem_volume(i_node, i_node)
            end do
        end do

        ! A node on a partition boundary has so far collected only the cells
        ! this rank holds; its control volume is completed here.
        if (associated(self%domain%mesh)) call self%domain%mesh%halo_sum_nodal(self%nodal_volume)

        nullify (fe)
        nullify (connectivity)
        if (allocated(coordinates)) deallocate (coordinates)
        if (allocated(unit_gp)) deallocate (unit_gp)
        if (allocated(elem_volume)) deallocate (elem_volume)
    end subroutine build_nodal_volume

end submodule ftcms_solve
