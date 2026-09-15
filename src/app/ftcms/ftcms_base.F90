submodule(app_ftcms) ftcms_base
    use :: core_parallel_reduce, only:set_node_ownership
    use :: domain_mesh_plex, only:type_dof_layout
    use :: core_types_topology_system_topology, only:type_system_topology
    use :: module_linalg, only:vector_norm2
    implicit none

    !> Highest BDF order the divided-difference local-error estimator implements.
    !> Error control refuses to run above this rather than reusing a lower-order
    !> defect, which would under-estimate the error by a factor O(h).
    integer(int32), parameter :: LTE_MAX_SUPPORTED_ORDER = 2

    !> Pin a node that crosses the freezing point to its critical temperature,
    !> following Hansson et al. (2004).
    !>
    !> Disabling this was tried and is worse, not better. The reset does place
    !> every crossing node exactly on the constitutive kink - 401 nodes in one
    !> Newton step, after which the hydraulic residual jumped 250-fold and the
    !> iterate stopped moving - but removing it does not recover that step and
    !> costs elsewhere: measured over the same interval, nonlinear-limit
    !> failures 1 -> 3, line-search failures 10 -> 41, mean inner iterations
    !> 6.88 -> 8.61, with the same simulated time reached. The kink is a symptom
    !> of the freezing-onset stall, not its cause.
    !> RE-TESTED 2026-08-01 with a true element-FD tangent and the joint step
    !> bound in place. The measurement above was taken when the tangent was
    !> Picard-lagged, where this reset acted as a stabiliser. It cannot serve
    !> that role for a Newton method: overwriting T after the update means the
    !> state the line search evaluates is not u + alpha*du, so the descent
    !> identity d/dalpha ||W R||^2 = -2||W R||^2 < 0 - which needs only an
    !> exact tangent and an exact linear solve, both verified here - no longer
    !> applies. That is exactly the observed failure: every alpha rejected with
    !> the merit approaching 1 from above.
    logical, parameter :: PHASE_ONSET_TEMPERATURE_RESET = .false.

    !> Report the per-block local time-truncation error alongside the max the
    !> controller acts on (compute_lte_error). Cheap: two reals per step.
    logical, parameter :: LTE_BLOCK_REPORT = .true.

    ! Bridge experiment (see plan spicy-sauteeing-scroll.md, WP3 refinement):
    ! where freezing impedance has collapsed hydraulic conductivity, a node no
    ! longer exchanges water with its neighbors within a time step, so its
    ! outer-loop target can be found directly by a local, flux-free Newton
    ! solve (calc_conserved_target + solve_local_conserved_equilibrium in
    ! fusion.F90) instead of many Picard/Anderson iterations toward the same
    ! point.
    !
    ! Gated by two dimensionless, problem-independent numerical-error-control
    ! parameters (not a case-fitted physical threshold - an ice-fraction
    ! proxy was tried first and measurably regressed the A/B test at one
    ! cutoff, then simply never triggered at another picked only to stop that
    ! regression; that case-fitting pattern is exactly what is prohibited):
    !
    ! 1. PHASE_LOCAL_EQ_FOURIER_GATE bounds Fo = (D_HH/C_eq)*dt/h^2, the
    !    Fourier number of the *inner* (fixed-ice) hydraulic transport
    !    operator at the node's representative element - D_HH/C_eq is
    !    verified dimensionally to be m^2/s (compute_diffusion_term's D_HH
    !    alone is m^2 s/kg; C_eq = dTheta/dP is 1/Pa = m s^2/kg; the ratio is
    !    m^2/s), and C_eq's dQi_dP is always zero (phase_systems.F90 - ice is
    !    outer-lagged, zero inner pressure tangent), so this is exactly the
    !    diffusivity the monolithic solve itself already uses at fixed ice.
    !    Fo << 1 means pressure information has not meaningfully diffused
    !    across the element this step.
    ! 2. Fo small only bounds the *internal* transport operator's relaxation,
    !    not actual boundary/source water exchange (gravity flow, prescribed
    !    flux, segregation sink can still move water regardless of Fo) - so
    !    nodes on any hydraulic boundary patch are excluded outright, and
    !    PHASE_LOCAL_EQ_TRANSPORT_EPS bounds the segregation-sink
    !    contribution over dt relative to the pore volume scale.
    !
    ! Both remain genuine, problem-independent safety margins (analogous to a
    ! CFL number or a Newton tolerance), not case-tuned constants.
    logical, parameter :: PHASE_USE_LOCAL_CONSERVED_EQ = .false.
    real(real64), parameter :: PHASE_LOCAL_EQ_FOURIER_GATE = 0.1d0
    real(real64), parameter :: PHASE_LOCAL_EQ_TRANSPORT_EPS = 1.0d-2
contains

    module subroutine initialize_type_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self
        type(type_input), save, target :: input
        type(type_ic_manager) :: ic

        integer(int32) :: max_bdf_order
        integer(int32), allocatable :: active_region_ids(:)
        integer(int32) :: num_nodes
        integer(int32) :: num_elements
        integer(int32) :: i
        integer(int32) :: computation_dimension
        integer(int32) :: num_total_dofs
        integer(int32) :: ierr
        real(real64), pointer, contiguous, dimension(:) :: phase_values

        type(type_config_control_manager) :: config_control_manager
        type(type_config_iteration) :: config_iteration
        type(type_config_time) :: config_time
        type(type_config_time_ats) :: config_time_ats
        type(type_config_output_manager) :: config_output_field
        type(type_config_output_manager) :: config_output_history
        type(type_config_acceleration) :: config_acceleration
        type(type_config_parallel_openmp) :: config_parallel_openmp

        type(type_config_ic) :: configs_ic(IC_TARGETS%NUM_ID)

        type(type_config_nodes) :: config_nodes
        type(type_config_elements) :: config_elements
        type(type_config_multicoloring) :: config_multicoloring
        type(type_config_elements), allocatable :: config_boundary_elements(:)

        type(type_solver_settings) :: solver_info

        type(type_config_output) :: config_output
        type(type_config_observation) :: config_observation
        type(type_config_overall) :: config_overall

        type(type_config_bc), allocatable :: config_bcs(:)

        self%last_phase_iterations = 1
        self%last_inner_iterations = 0
        self%last_nonlinear_work = 0
        self%last_solve_status = 0
        self%last_phase_metrics_available = .false.
        self%last_phase_increment_max = -1.0d0
        self%last_accepted_dt = 0.0d0

        call self%control%initialize()
        call self%control%profiler_record(TIME_RECORDS%START)
        call self%control%profiler_start(PROFILER_TYPES%TOTAL)

        ! call setup_handler()

        call self%control%profiler_start(PROFILER_TYPES%IO)
        call input%initialize()
        call self%control%profiler_stop(PROFILER_TYPES%IO)

        call input_translator%execute(input, config_control_manager)
        call input_translator%execute(input, config_iteration)
        call input_translator%execute(input, config_time)
        call input_translator%execute(input, config_time_ats)
        self%lte_rtol = config_time_ats%error_rtol
        self%lte_atol_thermal = config_time_ats%error_atol_temperature
        self%lte_atol_hydraulic = config_time_ats%error_atol_pressure
        self%lte_error_control_active = config_time_ats%active .and. config_time_ats%use_error_control
        call input_translator%execute(input, OUTPUT_TYPES%FIELD, config_output_field)
        call input_translator%execute(input, OUTPUT_TYPES%HISTORY, config_output_history)
        call input_translator%execute(input, config_acceleration)
        call input_translator%execute(input, config_parallel_openmp)

        call self%control%initialize(config_control_manager, config_iteration, config_time, config_time_ats, &
                                     config_output_field, config_output_history, &
                                     config_acceleration, config_parallel_openmp)

        if (self%is_active_thermal()) then
            call input_translator%execute(input, PHYSICS_TYPES%THERMAL, config_bcs)
            call self%bc(PHYSICS_TYPES%THERMAL%ID)%initialize(config_bcs)
            if (allocated(config_bcs)) deallocate (config_bcs)
        end if

        if (self%is_active_hydraulic()) then
            call input_translator%execute(input, PHYSICS_TYPES%HYDRAULIC, config_bcs)
            call self%bc(PHYSICS_TYPES%HYDRAULIC%ID)%initialize(config_bcs)
            if (allocated(config_bcs)) deallocate (config_bcs)
        end if

        call input_translator%execute(input, IC_TARGETS%POROSITY, configs_ic(IC_TARGETS%POROSITY%ID))
        call input_translator%execute(input, IC_TARGETS%THERMAL, configs_ic(IC_TARGETS%THERMAL%ID))
        call input_translator%execute(input, IC_TARGETS%HYDRAULIC, configs_ic(IC_TARGETS%HYDRAULIC%ID))
        call input_translator%execute(input, IC_TARGETS%MECHANICAL, configs_ic(IC_TARGETS%MECHANICAL%ID))
        call ic%initialize(configs_ic)

        if (input%output_settings%standard_output%print_progress) then
            call global_logger%configure(level=information_level, &
                                         time_stamp=.true., &
                                         max_width=0)
        else
            call global_logger%configure(level=warning_level, &
                                         time_stamp=.true., &
                                         max_width=0)
        end if

        call input_translator%execute(input, config_nodes)
        call input_translator%execute(input, config_elements, config_multicoloring)
        call input_translator%execute(input, config_boundary_elements)
        call self%domain%initialize(config_nodes, config_elements, config_multicoloring, config_boundary_elements, &
                                    input%basic%simulation_settings%calculate_type, config_control_manager%coupling_mode, &
                                    config_control_manager%compute_active)

        ! The DM the input layer read stays the one source of geometry.
        call self%domain%attach_mesh(input%geometry%mesh)

        ! Every global reduction from here on counts a shared node once.
        call set_node_ownership(input%geometry%mesh%node_is_owned)

        ! The mesh (connectivity, materials, adjacency) is static for the run:
        ! build the per-node distinct-material table once, right after the
        ! domain is fully built, for use by update_nodal_phases.
        call self%node_material_table%initialize(self%domain)

        ! Likewise build the static Gauss-point geometry cache used by the
        ! nodal gradient projection (calc_gradient).
        call self%gradient_cache%initialize(self%domain)

        call self%domain%get_num_fe(num_elements)

        call self%domain%get_start_dof_index(PHYSICS_TYPES%THERMAL, self%thermal_start_dof)
        call self%domain%get_start_dof_index(PHYSICS_TYPES%HYDRAULIC, self%hydraulic_start_dof)

        call self%domain%get_total_dofs(num_total_dofs)
        call self%domain%get_num_nodes(num_nodes)
        allocate (self%phase_onset_reset(num_nodes), source=.false.)

        block
            ! Domain-independent carrier injected into the system layer, so the
            ! Jacobian / residual no longer depend on type_domain directly.
            type(type_system_topology) :: topology

            call self%domain%export_topology(topology)
            ! The DM supplies the dof numbering and the matrix; there is no
            ! scatter map to build any more.
            call self%K%initialize(topology, config_control_manager%coupling_mode, input%geometry%mesh)
            call self%F%initialize(topology, config_control_manager%coupling_mode)
            call self%dirichlet_marker%initialize(topology, config_control_manager%coupling_mode)
            call self%du%initialize(topology, config_control_manager%coupling_mode)
        end block

        max_bdf_order = input%basic%solver_settings%bdf_order
        call self%porosity%initialize(num_nodes, max_bdf_order)
        call ic%apply(IC_TARGETS%POROSITY, self%porosity)

        if (self%is_active_thermal()) then
            call self%temperature%initialize(num_nodes, max_bdf_order)
            call ic%apply(IC_TARGETS%THERMAL, self%temperature)
        end if

        if (self%is_active_hydraulic()) then
            call self%pressure%initialize(num_nodes, max_bdf_order)
            call ic%apply(IC_TARGETS%HYDRAULIC, self%pressure)
        end if

        call self%Qw%initialize(num_nodes, max_bdf_order)
        call self%Qi%initialize(num_nodes, max_bdf_order)
        call self%Qa%initialize(num_nodes, max_bdf_order)
        call self%Qv%initialize(num_nodes, max_bdf_order)

        call self%domain%get_computation_dimension(computation_dimension)
        call input%geometry%mesh%get_active_region_info(active_region_ids, target_dim=computation_dimension)

        call self%thermal%initialize(input, active_region_ids)
        call self%hydraulic%initialize(input, active_region_ids)

        ! Apply initial Dirichlet boundary conditions to field variables
        call self%apply_bc()

        ! Initialize solver strictly from input settings.
        associate (linear_solver_settings => input%basic%solver_settings%linear_solver)
            block
                integer(int32) :: solver_size
                if (self%control%is_staggered()) then
                    solver_size = num_nodes
                else
                    solver_size = num_total_dofs
                end if
                call solver_info%set(solver_size, &
                                     linear_solver_settings%tolerance, &
                                     linear_solver_settings%max_iterations, &
                                     petsc_options=trim(linear_solver_settings%petsc_options))
            end block
            call create_solver(self%solver, solver_info, ierr)
            call hand_dof_map_to_solver(self%solver, self%K)

            if (self%is_active_thermal()) then
                if (self%control%is_staggered()) then
                    call solver_info%set(num_nodes, &
                                         linear_solver_settings%tolerance, &
                                         linear_solver_settings%max_iterations, &
                                         relative_tolerance=5.0d-2, &
                                         petsc_options=trim(linear_solver_settings%petsc_options))
                end if
                call create_solver(self%solver_thermal, solver_info, ierr)
                call hand_dof_map_to_solver(self%solver_thermal, self%K, PHYSICS_TYPES%THERMAL)
            end if
        end associate

        ! Populate initial phase variables from initial T/P/porosity before first output.
        call self%update_variables()
        ! Sync previous-step values with IC so the BDF transient term is zero at t=0.
        ! Only active-physics state variables are initialized (see guarded init above);
        ! inactive ones are not retained, matching the is_active_* gating used throughout.
        nullify (phase_values)
        if (self%is_active_thermal()) then
            call self%temperature%get_current(phase_values)
            if (associated(phase_values)) call self%temperature%set_previous(phase_values)
            nullify (phase_values)
        end if
        if (self%is_active_hydraulic()) then
            call self%pressure%get_current(phase_values)
            if (associated(phase_values)) call self%pressure%set_previous(phase_values)
            nullify (phase_values)
        end if
        call self%Qw%get_current(phase_values)
        if (associated(phase_values)) call self%Qw%set_previous(phase_values)
        nullify (phase_values)
        call self%Qi%get_current(phase_values)
        if (associated(phase_values)) call self%Qi%set_previous(phase_values)
        nullify (phase_values)
        call self%Qa%get_current(phase_values)
        if (associated(phase_values)) call self%Qa%set_previous(phase_values)
        nullify (phase_values)
        call self%Qv%get_current(phase_values)
        if (associated(phase_values)) call self%Qv%set_previous(phase_values)
        nullify (phase_values)

        call self%initialize_lte_history()

        call input_translator%execute(input, config_output)
        call input_translator%execute(input, config_observation)
        call input_translator%execute(input, config_overall)

        if (config_observation%point_type == OUTPUT_OBSERVATION_TYPES%COORDINATES) then
            if (allocated(config_observation%observation_geometries)) then
                do i = 1, size(config_observation%observation_geometries)
                    call self%domain%get_config(config_observation%observation_geometries(i))
                end do
            end if
        end if

        call self%output%initialize(config_output, config_observation, config_overall, &
                                    input%geometry%mesh)
        call self%output_fields()
        call self%output_history()

        call initialize_assimilation_ftcms(self, input%input_path)

        call open_solver_history_log(self)

    end subroutine initialize_type_ftcms

    !> Open Output/solver_history.log (rank 0 only): one record per time-step
    !> attempt with the diagnostics that are otherwise invisible from outside
    !> (iterations, acceptance, the gates that decided it, increments, LTE).
    subroutine open_solver_history_log(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        character(:), allocatable :: project_path_env
        character(*), parameter :: PROJECT_ENV = "FTCMS_PROJECT_PATH"
        integer(int32) :: myrank, ierr, ios

        self%solver_history_unit = -1
        call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
        if (myrank /= 0) return

        call get_env_string(PROJECT_ENV, project_path_env)
        call modify_path_format(project_path_env)

        open (newunit=self%solver_history_unit, &
              file=trim(project_path_env)//"Output/solver_history.log", &
              status="replace", action="write", iostat=ios)
        if (ios /= 0) then
            self%solver_history_unit = -1
            return
        end if
        write (self%solver_history_unit, '(A)') &
            "# FTCMS solver history schema=4: one record per time-step attempt"
        write (self%solver_history_unit, '(A)') &
            "# attempt accepted_step time_start_s time_trial_s time_accepted_s dt_used_s dt_next_s" // &
            " accepted status inner_iter ats_iter loc_T loc_H bal_T bal_H dq_eff du_T_max du_p_max" // &
            " omega lte_rel"
        write (self%solver_history_unit, '(A)') &
            "# loc_*/bal_*/dq_eff are the acceptance gates of the last iterate: each passes at <= 1," // &
            " -1 means not evaluated. status is how the attempt ended, not which gate refused it"
        write (self%solver_history_unit, '(A)') &
            "# du_*_max are the last increment magnitudes [K] and [Pa]; lte_rel is -1 when the step" // &
            " never reached the time-error test"
        flush (self%solver_history_unit)
    end subroutine open_solver_history_log

    subroutine initialize_assimilation_ftcms(self, input_path)
        use :: json_module, only: json_file
        implicit none
        class(type_ftcms), intent(inout) :: self
        character(*), intent(in) :: input_path

        type(json_file) :: jf
        type(type_da_config) :: da_cfg
        character(len=512) :: da_file
        character(:), allocatable :: str_val
        logical :: found, file_ok
        integer(int32) :: ierr, myrank
        logical :: da_enabled_flag

        call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

        da_file = trim(input_path)//'DataAssimilation.json'
        inquire (file=trim(da_file), exist=file_ok)
        if (.not. file_ok) return

        call jf%initialize()
        call jf%load(trim(da_file))
        if (jf%failed()) return

        da_enabled_flag = .true.
        call jf%get('enabled', da_enabled_flag, found)
        if (found .and. .not. da_enabled_flag) then
            call jf%destroy()
            return
        end if

        call jf%get('ensemble.size', da_cfg%ensemble_size, found)
        call jf%get('ensemble.max_height', da_cfg%max_height, found)
        call jf%get('ensemble.num_nodes', da_cfg%num_nodes, found)
        call jf%get('observation.csv_file', str_val, found)
        if (found .and. allocated(str_val)) da_cfg%csv_file = str_val
        call jf%get('observation.upper_bc_file', str_val, found)
        if (found .and. allocated(str_val)) da_cfg%upper_bc_file = str_val
        call jf%get('observation.interval_seconds', da_cfg%interval_seconds, found)
        call jf%get('observation.sigma_T', da_cfg%sigma_T, found)
        call jf%get('observation.sigma_q', da_cfg%sigma_q, found)
        call jf%get('observation.sigma_U', da_cfg%sigma_U, found)
        call jf%get('surface.z0', da_cfg%z0, found)
        call jf%get('surface.albedo', da_cfg%albedo, found)
        call jf%get('surface.emissivity', da_cfg%emissivity, found)
        call jf%get('surface.Pmin', da_cfg%Pmin, found)
        call jf%get('surface.Pmax', da_cfg%Pmax, found)
        call jf%get('surface.lambda_soil', da_cfg%lambda_soil, found)
        call jf%get('surface.stomatal_resistance', da_cfg%stomatal_resistance, found)
        call jf%get('reference_datetime', str_val, found)
        if (found .and. allocated(str_val)) da_cfg%reference_datetime = str_val
        call jf%get('solar.latitude', da_cfg%latitude, found)
        call jf%get('solar.longitude', da_cfg%longitude, found)
        call jf%get('solar.tau_atm', da_cfg%tau_atm, found)
        call jf%get('solar.utc_offset_hours', da_cfg%utc_offset_hours, found)
        call jf%destroy()

        call self%assimilation%initialize(da_cfg, &
            bc_entity_thermal=3, bc_entity_hydraulic=3)
        self%assimilation_enabled = .true.

        if (myrank == 0) then
            write (*, '(A)') '[DA] Data assimilation enabled from DataAssimilation.json.'
        end if
    end subroutine initialize_assimilation_ftcms

    module subroutine output_fields_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        integer(int32) :: iter
        real(real64) :: current_time

        real(real64), pointer, contiguous, dimension(:) :: temperature
        real(real64), pointer, contiguous, dimension(:) :: pressure
        real(real64), pointer, contiguous, dimension(:) :: water_content
        real(real64), pointer, contiguous, dimension(:) :: ice_pore
        real(real64), pointer, contiguous, dimension(:) :: vapor_content

        type(type_coordinate_array_dp) :: water_flux_arr
        integer(int32) :: num_nodes_wf, num_fe_wf, i_elem_wf, i_local_wf, node_id_wf
        integer(int32), pointer, contiguous :: conn_wf(:)
        type(type_state) :: state_wf
        type(type_coordinate_dp), pointer :: wf_ptr
        logical :: has_water_flux

        call self%control%profiler_start(PROFILER_TYPES%IO)

        nullify (temperature)
        nullify (pressure)
        nullify (ice_pore)
        nullify (vapor_content)
        nullify (water_content)

        call self%control%get_time(current_time)

        if (self%control%is_output_triggered(OUTPUT_TYPES%FIELD, current_time)) then
            call self%control%get_output_step(OUTPUT_TYPES%FIELD, iter)
            if (self%is_active_thermal()) then
                call self%temperature%get_current(temperature)
            end if
            if (self%is_active_hydraulic()) then
                call self%pressure%get_current(pressure)
            end if
            call self%Qw%get_current(water_content)
            call self%Qi%get_current(ice_pore)
            call self%Qv%get_current(vapor_content)

            has_water_flux = .false.
            if (self%is_active_hydraulic()) then
                nullify (conn_wf)
                nullify (wf_ptr)
                call self%domain%get_num_nodes(num_nodes_wf)
                call self%domain%get_num_fe(num_fe_wf)
                call water_flux_arr%initialize(num_nodes_wf)
                do i_elem_wf = 1, num_fe_wf
                    call self%domain%get_fe_connectivity(i_elem_wf, conn_wf)
                    do i_local_wf = 1, size(conn_wf)
                        node_id_wf = conn_wf(i_local_wf)
                        if (node_id_wf < 1) cycle
                        call self%set_state(node_id_wf, i_elem_wf, state_wf, calc_physics=.true.)
                        nullify (wf_ptr)
                        call state_wf%get(water_flux=wf_ptr)
                        if (associated(wf_ptr)) then
                            water_flux_arr%x(node_id_wf) = wf_ptr%x
                            water_flux_arr%y(node_id_wf) = wf_ptr%y
                            water_flux_arr%z(node_id_wf) = wf_ptr%z
                        end if
                    end do
                    nullify (conn_wf)
                end do
                has_water_flux = .true.
            end if

            if (has_water_flux) then
                call self%output%output_fields(file_counts=iter, &
                                               temperature=temperature, &
                                               water_content=water_content, &
                                               ice_content=ice_pore, &
                                               vapor_content=vapor_content, &
                                               pressure=pressure, &
                                               water_flux=water_flux_arr, &
                                               time=current_time)
            else
                call self%output%output_fields(file_counts=iter, &
                                               temperature=temperature, &
                                               water_content=water_content, &
                                               ice_content=ice_pore, &
                                               vapor_content=vapor_content, &
                                               pressure=pressure, &
                                               time=current_time)
            end if

            call self%control%update_output(OUTPUT_TYPES%FIELD, current_time)

            call water_flux_arr%destroy()
            nullify (temperature)
            nullify (pressure)
            nullify (ice_pore)
            nullify (water_content)
            nullify (vapor_content)
        end if

        call self%control%profiler_stop(PROFILER_TYPES%IO)
    end subroutine output_fields_ftcms

    module subroutine output_history_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        real(real64) :: current_time, current_time_converted
        real(real64), pointer, contiguous, dimension(:) :: temperature
        real(real64), pointer, contiguous, dimension(:) :: pressure
        real(real64), pointer, contiguous, dimension(:) :: water_content
        real(real64), pointer, contiguous, dimension(:) :: ice_pore
        real(real64), pointer, contiguous, dimension(:) :: vapor_content
        type(type_coordinate_array_dp) :: liq_flux_arr, vap_flux_arr
        integer(int32) :: num_nodes_hf, i_obs_fe, i_local_hf, node_id_hf
        integer(int32), pointer, contiguous :: conn_hf(:)
        integer(int32), allocatable :: obs_fe_ids(:)
        type(type_state) :: state_hf
        type(type_coordinate_dp), pointer :: wf_liq_ptr, wf_vap_ptr
        logical :: has_flux

        call self%control%profiler_start(PROFILER_TYPES%IO)

        nullify (temperature)
        nullify (pressure)
        nullify (water_content)
        nullify (ice_pore)
        nullify (vapor_content)

        call self%control%get_time(current_time)

        if (self%control%is_output_triggered(OUTPUT_TYPES%HISTORY, current_time)) then
            if (self%is_active_thermal()) then
                call self%temperature%get_current(temperature)
            end if
            if (self%is_active_hydraulic()) then
                call self%pressure%get_current(pressure)
            end if
            call self%Qw%get_current(water_content)
            call self%Qi%get_current(ice_pore)
            call self%Qv%get_current(vapor_content)

            ! Compute per-node liquid and vapor fluxes for observation elements only
            has_flux = .false.
            if (self%is_active_hydraulic()) then
                call self%output%get_obs_fe_ids(obs_fe_ids)
                if (size(obs_fe_ids) > 0 .and. all(obs_fe_ids > 0)) then
                    nullify (conn_hf)
                    call self%domain%get_num_nodes(num_nodes_hf)
                    call liq_flux_arr%initialize(num_nodes_hf)
                    call vap_flux_arr%initialize(num_nodes_hf)
                    do i_obs_fe = 1, size(obs_fe_ids)
                        call self%domain%get_fe_connectivity(obs_fe_ids(i_obs_fe), conn_hf)
                        do i_local_hf = 1, size(conn_hf)
                            node_id_hf = conn_hf(i_local_hf)
                            if (node_id_hf < 1) cycle
                            nullify (wf_liq_ptr, wf_vap_ptr)
                            call self%set_state(node_id_hf, obs_fe_ids(i_obs_fe), state_hf, calc_physics=.true.)
                            call state_hf%get(water_flux=wf_liq_ptr, vapor_flux=wf_vap_ptr)
                            if (associated(wf_liq_ptr)) then
                                liq_flux_arr%x(node_id_hf) = wf_liq_ptr%x
                                liq_flux_arr%y(node_id_hf) = wf_liq_ptr%y
                                liq_flux_arr%z(node_id_hf) = wf_liq_ptr%z
                            end if
                            if (associated(wf_vap_ptr)) then
                                vap_flux_arr%x(node_id_hf) = wf_vap_ptr%x
                                vap_flux_arr%y(node_id_hf) = wf_vap_ptr%y
                                vap_flux_arr%z(node_id_hf) = wf_vap_ptr%z
                            end if
                        end do
                        nullify (conn_hf)
                    end do
                    has_flux = .true.
                end if
                if (allocated(obs_fe_ids)) deallocate (obs_fe_ids)
            end if

            call self%control%get_output_time(OUTPUT_TYPES%HISTORY, current_time, current_time_converted)
            if (has_flux) then
                call self%output%output_history(time=current_time_converted, &
                                                temperature=temperature, &
                                                water_content=water_content, &
                                                ice_content=ice_pore, &
                                                vapor_content=vapor_content, &
                                                pressure=pressure, &
                                                water_flux=liq_flux_arr, &
                                                vapor_flux=vap_flux_arr)
            else
                call self%output%output_history(time=current_time_converted, &
                                                temperature=temperature, &
                                                water_content=water_content, &
                                                ice_content=ice_pore, &
                                                vapor_content=vapor_content, &
                                                pressure=pressure)
            end if
            call self%control%update_output(OUTPUT_TYPES%HISTORY, current_time)

            if (has_flux) then
                call liq_flux_arr%destroy()
                call vap_flux_arr%destroy()
            end if
            nullify (water_content)
            nullify (ice_pore)
            nullify (vapor_content)
            nullify (temperature)
            nullify (pressure)
        end if

        call self%control%profiler_stop(PROFILER_TYPES%IO)
    end subroutine output_history_ftcms

    module subroutine get_variable_increment_ftcms(self, variable_id, variable)
        implicit none
        class(type_ftcms), intent(inout) :: self
        type(type_constant_id), intent(in) :: variable_id
        real(real64), intent(inout), allocatable :: variable(:)

        real(real64), pointer, dimension(:) :: du

        integer(int32) :: target_dof
        integer(int32) :: num_nodes, num_dofs_per_node
        integer(int32) :: start_idx, end_idx

        nullify (du)
        call deallocate_array(variable)

        if (.not. PHYSICS_TYPES%is_valid(variable_id)) then
            call allocate_array(variable, 0)
            return
        end if

        if (self%control%is_physics_active(variable_id)) then
            call self%domain%get_num_nodes(num_nodes)
            call allocate_array(variable, num_nodes)

            if (self%control%is_staggered()) then
                du => self%du%get_data(variable_id%ID)
                if (.not. associated(du)) then
                    variable(:) = 0.0d0
                else
                    variable(:) = du(1:num_nodes)
                    nullify (du)
                end if
            else
                du => self%du%get_data()
                if (.not. associated(du)) then
                    variable(:) = 0.0d0
                else
                    call self%domain%get_num_dof_per_node(num_dofs_per_node)
                call self%domain%get_start_dof_index(variable_id, target_dof)
                    start_idx = target_dof
                    end_idx = num_dofs_per_node * (num_nodes - 1) + target_dof
                    variable(:) = du(start_idx:end_idx:num_dofs_per_node)
                    nullify (du)
                end if
            end if
        else
            call allocate_array(variable, 0)
        end if

    end subroutine get_variable_increment_ftcms

    module subroutine get_variable_residual_ftcms(self, variable_id, variable)
        implicit none
        class(type_ftcms), intent(inout) :: self
        type(type_constant_id), intent(in) :: variable_id
        real(real64), intent(inout), allocatable :: variable(:)

        real(real64), pointer, dimension(:) :: F

        integer(int32) :: target_dof
        integer(int32) :: num_nodes, num_dofs_per_node
        integer(int32) :: start_idx, end_idx

        nullify (F)
        call deallocate_array(variable)

        if (.not. PHYSICS_TYPES%is_valid(variable_id)) then
            call allocate_array(variable, 0)
            return
        end if

        if (self%control%is_physics_active(variable_id)) then
            call self%domain%get_num_nodes(num_nodes)
            call allocate_array(variable, num_nodes)

            if (self%control%is_staggered()) then
                F => self%F%get_data(variable_id%ID)
                if (.not. associated(F)) then
                    variable(:) = 0.0d0
                else
                    variable(:) = F(1:num_nodes)
                    nullify (F)
                end if
            else
                F => self%F%get_data()
                if (.not. associated(F)) then
                    variable(:) = 0.0d0
                else
                    call self%domain%get_num_dof_per_node(num_dofs_per_node)
                    call self%domain%get_start_dof_index(variable_id, target_dof)
                    start_idx = target_dof
                    end_idx = num_dofs_per_node * (num_nodes - 1) + target_dof
                    variable(:) = F(start_idx:end_idx:num_dofs_per_node)
                    nullify (F)
                end if
            end if

            ! This rank assembled only its own cells, so a node on a partition
            ! boundary holds a partial residual. Every convergence gate reads
            ! this array, so it is completed here rather than at each gate.
            if (associated(self%domain%mesh)) call self%domain%mesh%halo_sum_nodal(variable)
        else
            call allocate_array(variable, 0)
        end if

    end subroutine get_variable_residual_ftcms

    module subroutine set_state_ftcms(self, node_id, element_id, state, calc_physics, include_fluxes)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: node_id
        integer(int32), intent(in) :: element_id
        type(type_state), intent(inout) :: state
        logical, intent(in), optional :: calc_physics
        logical, intent(in), optional :: include_fluxes

        integer(int32) :: material_id
        integer(int32) :: bdf_order
        integer(int32) :: start_dof_thermal, start_dof_hydraulic
        real(real64) :: temperature, pressure, porosity
        type(type_coordinate_dp) :: grad_T, grad_P
        real(real64) :: temperature_history(8), pressure_history(8), porosity_history(8)
        real(real64) :: ice_content_history(8)

        logical :: do_calc, do_fluxes
        logical :: temperature_set, pressure_set

        ! Default: compute physics (for backward compatibility)
        do_calc = .true.
        if (present(calc_physics)) do_calc = calc_physics
        ! Default: compute fluxes (for backward compatibility)
        do_fluxes = .true.
        if (present(include_fluxes)) do_fluxes = include_fluxes
        call state%reset()

        ! Set before any physics is evaluated, and regardless of calc_physics:
        ! states built with calc_physics=.false. are updated in bulk later, and
        ! that update reads the constitutive laws this parameter softens.
        call state%continuation_lambda%set(self%control%get_homotopy_lambda())

        call self%control%get_bdf_coeffs(bdf_order=bdf_order)

        start_dof_thermal = self%thermal_start_dof
        if (start_dof_thermal > 0) then
            call self%temperature%get_current(node_id, temperature)
            call self%temperature%get_current_gradient(node_id, grad_T)
            call self%temperature%get_history(node_id, temperature_history)
        else
            temperature = 0.0d0
            call grad_T%reset()
            temperature_history = 0.0d0
        end if
        call state%set(temperature=temperature, &
                       grad_T=grad_T, &
                       temperature_history=temperature_history(1:bdf_order + 1))

        start_dof_hydraulic = self%hydraulic_start_dof
        ! No hydraulic dof means no liquid-pressure field; the closure then sees
        ! zero matric suction, which is what makes freezing temperature-driven.
        if (start_dof_hydraulic > 0) then
            call self%pressure%get_current(node_id, pressure)
            call self%pressure%get_current_gradient(node_id, grad_P)
            call self%pressure%get_history(node_id, pressure_history)
        else
            pressure = 0.0d0
            call grad_P%reset()
            pressure_history = 0.0d0
        end if
        call state%set(pressure=pressure, &
                       grad_P=grad_P, &
                       pressure_history=pressure_history(1:bdf_order + 1))

        call self%porosity%get_current(node_id, porosity)
        call self%porosity%get_history(node_id, porosity_history)
        call state%set(porosity=porosity, &
                       porosity_history=porosity_history(1:bdf_order + 1))

        block
            real(real64) :: qw_val, qi_val, qa_val, qv_val
            call self%Qw%get_current(node_id, qw_val)
            call self%Qi%get_current(node_id, qi_val)
            call self%Qa%get_current(node_id, qa_val)
            call self%Qv%get_current(node_id, qv_val)
            call state%water_content%set(qw_val)
            call state%ice_content%set(qi_val)
            call state%air_content%set(qa_val)
            call state%vapor_content%set(qv_val)
        end block
        call self%Qi%get_history(node_id, ice_content_history)
        call state%ice_content_history%set(ice_content_history(1:bdf_order + 1))

        call state%temperature%get(temperature, temperature_set)
        call state%pressure%get(pressure, pressure_set)
        if (.not. temperature_set .or. .not. pressure_set) then
            !$omp critical (ftcms_state_diag)
            write (*, '(A,I0,A,I0,A,L1,A,L1,A,I0,A,I0)') 'Error: set_state_ftcms unset primary state. node=', node_id, &
                ', elem=', element_id, ', T_set=', temperature_set, ', P_set=', pressure_set, &
                ', T_dof=', start_dof_thermal, ', H_dof=', start_dof_hydraulic
            !$omp end critical (ftcms_state_diag)
            error stop 'set_state_ftcms: temperature/pressure unset before constitutive update.'
        end if

        ! Run expensive physics computations only when flagged
        if (do_calc) then
            call self%domain%get_material_id(element_id, material_id)
            call self%update_physical_properties(material_id, state, include_fluxes=do_fluxes)
        end if

    end subroutine set_state_ftcms

    module subroutine set_states_from_connectivity_ftcms(self, connectivity, element_id, states, calc_physics)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: connectivity(:)
        integer(int32), intent(in) :: element_id
        type(type_state), intent(inout) :: states(:)
        logical, intent(in), optional :: calc_physics

        integer(int32) :: i, node_id
        integer(int32) :: material_id
        logical :: do_calc

        if (size(states) /= size(connectivity)) then
            error stop 'set_states_from_connectivity_ftcms: size(states) /= size(connectivity)'
        end if

        do_calc = .true.
        if (present(calc_physics)) do_calc = calc_physics

        do i = 1, size(connectivity)
            node_id = connectivity(i)
            if (node_id < 1) then
                write (*, '(A, I0, A, I0, A, I0)') 'invalid node_id=', node_id, ', elem=', element_id, ', local=', i
                error stop 'set_states_from_connectivity_ftcms: node_id out of range'
            end if
            call self%set_state(node_id, element_id, states(i), calc_physics=.false.)
        end do

        if (do_calc) then
            call self%domain%get_material_id(element_id, material_id)
            call self%update_physical_properties_bulk(material_id, states)
        end if
    end subroutine set_states_from_connectivity_ftcms

    ! Update all physical quantities (phase, fluxes) for a given state
    module subroutine update_physical_properties_ftcms(self, material_id, state, include_fluxes)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(inout) :: state
        logical, intent(in), optional :: include_fluxes

        type(type_coordinate_dp) :: water_flux, vapor_flux
        type(type_coordinate_dp), pointer :: grad_T, grad_P
        type(type_coordinate_dp), target :: zero_grad
        logical :: do_fluxes

        ! Default: compute fluxes (for backward compatibility)
        do_fluxes = .true.
        if (present(include_fluxes)) do_fluxes = include_fluxes

        nullify (grad_T)
        nullify (grad_P)
        call zero_grad%set(0.0d0, 0.0d0, 0.0d0)

        ! 1. Phase change calculation (Ice/Water Content)
        call self%thermal%update_water_phases(material_id, state)

        ! Phases-only mode: skip the flux evaluation entirely, leaving the
        ! state's flux fields unset. Only valid for callers that never read
        ! state%water_flux / state%vapor_flux afterwards (e.g. the nodal
        ! phase update, which extracts phase contents only).
        if (.not. do_fluxes) return

        ! 2. Flux calculation (used in advection and diffusion)
        if (self%control%is_target(PHYSICS_TYPES%HYDRAULIC, material_id)) then
            call state%grad_T%get(grad_T)
            call state%grad_P%get(grad_P)
            ! Guard against null pointer when gradient is not set
            if (.not. associated(grad_T)) grad_T => zero_grad
            if (.not. associated(grad_P)) grad_P => zero_grad
            call self%calc_water_flux(material_id, state, grad_T, grad_P, water_flux)
            call self%calc_vapor_flux(material_id, state, grad_T, grad_P, vapor_flux)
        else
            call water_flux%set(0.0d0, 0.0d0, 0.0d0)
            call vapor_flux%set(0.0d0, 0.0d0, 0.0d0)
        end if

        call state%set(water_flux=water_flux, vapor_flux=vapor_flux)
    end subroutine update_physical_properties_ftcms

    module subroutine update_physical_properties_bulk_ftcms(self, material_id, states)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: material_id
        type(type_state), intent(inout) :: states(:)

        integer(int32) :: i

        do i = 1, size(states)
            call self%update_physical_properties(material_id, states(i))
        end do
    end subroutine update_physical_properties_bulk_ftcms

    module subroutine shift_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        call self%control%profiler_start(PROFILER_TYPES%SETUP)

        if (self%is_active_thermal()) then
            call self%temperature%advance()
        end if

        if (self%is_active_hydraulic()) then
            call self%pressure%advance()
        end if

        call self%porosity%advance()
        call self%Qw%advance()
        call self%Qi%advance()
        call self%Qa%advance()
        call self%Qv%advance()

        call self%control%profiler_stop(PROFILER_TYPES%SETUP)
    end subroutine shift_ftcms

    module subroutine reflect_variables_ftcms(self, step_scale)
        implicit none
        class(type_ftcms), intent(inout) :: self
        real(real64), intent(in), optional :: step_scale

        integer(int32) :: iter
        real(real64), pointer, contiguous, dimension(:) :: bdf_coeffs
        integer(int32) :: bdf_order
        real(real64), pointer, contiguous, dimension(:) :: current
        real(real64), allocatable :: du(:)
        real(real64), allocatable :: current_prev(:)
        real(real64), allocatable :: du_eff(:)
        real(real64), allocatable :: temperature_before_phase(:)

        real(real64) :: relaxation_factor
        logical :: is_none
        logical :: is_conserved_mode

        ! --- Anderson(1) mixing of the conserved coupled update ---
        ! x_{k+1} = x_k + w*g_k - gamma*[(x_k - x_{k-1}) + w*(g_k - g_{k-1})],
        ! gamma = <g_k, g_k - g_{k-1}>_W / ||g_k - g_{k-1}||_W^2 (Walker & Ni, 2011),
        ! evaluated jointly over (T, p) with head-equivalent weighting of p so both
        ! blocks are mixed consistently. Pure fixed-point acceleration: no Jacobian,
        ! the modified-Picard structure and the adaptive damping w are unchanged.
        logical :: aa_active
        real(real64) :: aa_gamma
        real(real64) :: aa_alpha_joint
        real(real64), allocatable :: aa_step_T(:), aa_step_P(:)
        ! Anderson(1) replaces the computed direction with its own
        ! extrapolation. With the element tangent verified against a finite
        ! difference and the linear solve verified to 1e-8, -J^-1 R is a descent
        ! direction for ||R||^2 by identity, so a search reporting no descent is
        ! not searching along it. Disabled to test exactly that.
        logical, parameter :: AA_ENABLED = .false.
        real(real64), parameter :: AA_GAMMA_MAX = 2.0d0
        real(real64), parameter :: AA_STEP_GROWTH_MAX = 4.0d0
        real(real64), parameter :: AA_WEIGHT_P = 1.0d0 / 9.81d3  ! pressure-head scaling [m/Pa]
        ! Require the nonsmooth phase return to be one order inside its
        ! 1e-3 content acceptance width before reusing a smooth-map secant.
        real(real64), parameter :: AA_PHASE_RESTART_CONTENT = 1.0d-4

        real(real64) :: max_du, alpha
        real(real64) :: step_bound_lambda
        real(real64), parameter :: PICARD_MAX_DT_STEP = 2.0d1
        real(real64), parameter :: PICARD_MAX_DT_STEP_PHASE = 0.5d0   ! K near T_melt
        real(real64), parameter :: PICARD_PHASE_ZONE = 2.0d0           ! °C half-width
        ! Legacy step caps (used only by the non-conserved convergence modes; the
        ! universal conserved mode below uses adaptive under-relaxation instead).
        real(real64), parameter :: PICARD_MAX_DP_STEP = 5.0d3
        real(real64), parameter :: TEMP_MIN_C = WALL_TEMP_MIN_C
        real(real64), parameter :: TEMP_MAX_C = WALL_TEMP_MAX_C
        real(real64), parameter :: PRESS_MIN_PA = WALL_PRESS_MIN_PA
        real(real64), parameter :: PRESS_MAX_PA = WALL_PRESS_MAX_PA

        call self%control%profiler_start(PROFILER_TYPES%SETUP)

        nullify (current)
        nullify (bdf_coeffs)

        call self%control%get_nonlinear_iter(iter)

        call self%control%get_bdf_coeffs(bdf_order, bdf_coeffs)

        is_none = self%control%is_none()
        is_conserved_mode = self%control%is_conserved()

        ! Anderson(1): compute the joint mixing coefficient from this iterate's
        ! (T, p) increments and the stored previous pair, BEFORE any update.
        aa_gamma = 0.0d0
        aa_active = AA_ENABLED .and. is_conserved_mode .and. (.not. is_none) .and. &
                    (.not. self%control%is_staggered()) .and. &
                    self%is_active_thermal() .and. self%is_active_hydraulic()
        if (aa_active .and. self%last_phase_metrics_available .and. &
            self%last_phase_increment_max > AA_PHASE_RESTART_CONTENT) then
            ! The local active-set return is still moving materially, so the
            ! previous T-p secant pair belongs to a different reduced map.
            ! Restart AA(1) without disabling it once the phase defect has
            ! entered the discretization tolerance.
            self%aa_has_prev = .false.
            self%aa_gnorm_prev = -1.0d0
            aa_active = .false.
        end if
        aa_alpha_joint = 1.0d0
        if (aa_active) then
            relaxation_factor = self%control%get_conserved_relaxation()
            if (present(step_scale)) relaxation_factor = relaxation_factor * step_scale
            call prepare_coupled_aa_step(self, relaxation_factor, aa_gamma, aa_step_T, aa_step_P, &
                                         aa_alpha_joint, aa_active)
        end if

        ! --- Joint variable-wise trust region -------------------------------
        ! One factor for both blocks, computed before either is updated, so the
        ! applied update stays parallel to the direction the line search was
        ! handed (see STEP_BOUND_ENABLED). The validity-wall factors are folded
        ! in here for the same reason. The increment is bounded, never the
        ! resulting value: clamping a value would break the correspondence
        ! between the state and the direction the merit was evaluated along.
        step_bound_lambda = 1.0d0
        if (STEP_BOUND_ENABLED .and. is_conserved_mode) then
            block
                real(real64), allocatable :: du_probe(:)
                real(real64), pointer, contiguous, dimension(:) :: probe_current
                real(real64) :: relax_probe, largest

                ! step_scale is deliberately NOT applied here. The bound must
                ! fix the DIRECTION once, and the line search must then scale
                ! that bounded direction: u_trial = u + alpha*(lambda*du).
                ! Folding step_scale into the probe instead makes
                ! lambda = dp_max/|omega*alpha*du| and the applied update
                ! omega*alpha*lambda*du = dp_max*du/|du| - independent of
                ! alpha, so every trial above the binding threshold evaluates
                ! the SAME state. Measured at the stalled iterate: lambda_bound
                ! = 5e4/4.151e6 = 1.2e-2, so the first seven halvings were
                ! duplicates and the search could never reach the descent
                ! region near an effective 1.2e-5.
                relax_probe = self%control%get_conserved_relaxation()

                nullify (probe_current)
                if (self%is_active_thermal()) then
                    call self%get_variable_increment(PHYSICS_TYPES%THERMAL, du_probe)
                    call self%temperature%get_current(probe_current)
                    if (allocated(du_probe) .and. associated(probe_current)) then
                        if (size(du_probe) > 0) then
                            du_probe(:) = relax_probe * du_probe(:)
                            largest = maxval(abs(du_probe))
                            if (largest > STEP_BOUND_DT) then
                                step_bound_lambda = min(step_bound_lambda, STEP_BOUND_DT / largest)
                            end if
                            step_bound_lambda = min(step_bound_lambda, &
                                                    bounded_step_factor(probe_current, du_probe, &
                                                                        TEMP_MIN_C, TEMP_MAX_C))
                        end if
                    end if
                    nullify (probe_current)
                    if (allocated(du_probe)) deallocate (du_probe)
                end if

                if (self%is_active_hydraulic()) then
                    call self%get_variable_increment(PHYSICS_TYPES%HYDRAULIC, du_probe)
                    call self%pressure%get_current(probe_current)
                    if (allocated(du_probe) .and. associated(probe_current)) then
                        if (size(du_probe) > 0) then
                            du_probe(:) = relax_probe * du_probe(:)
                            largest = maxval(abs(du_probe))
                            if (largest > STEP_BOUND_DP) then
                                step_bound_lambda = min(step_bound_lambda, STEP_BOUND_DP / largest)
                            end if
                            step_bound_lambda = min(step_bound_lambda, &
                                                    bounded_step_factor(probe_current, du_probe, &
                                                                        PRESS_MIN_PA, PRESS_MAX_PA))
                        end if
                    end if
                    nullify (probe_current)
                    if (allocated(du_probe)) deallocate (du_probe)
                end if
            end block
        end if

        if (self%is_active_thermal()) then
            call self%get_variable_increment(PHYSICS_TYPES%THERMAL, du)
            call self%temperature%get_current(current)
            if (associated(current)) then
                call allocate_array(current_prev, size(current))
                current_prev(:) = current(:)
                call allocate_array(temperature_before_phase, size(current))
                temperature_before_phase(:) = current(:)

                if (allocated(du) .and. size(du) > 0) then
                    if (is_conserved_mode) then
                        ! Conserved mode uses one coupled Picard damping factor for
                        ! T and p. The factor is adapted from the conserved-quantity
                        ! contraction rate after each accepted nonlinear iterate.
                        relaxation_factor = self%control%get_conserved_relaxation()
                        if (present(step_scale)) relaxation_factor = relaxation_factor * step_scale
                        call allocate_array(du_eff, size(du))
                        if (aa_active .and. allocated(aa_step_T) .and. size(aa_step_T) == size(du)) then
                            du_eff(:) = aa_step_T(:)
                            alpha = aa_alpha_joint
                        else
                            du_eff(:) = relaxation_factor * du(:)
                            if (STEP_BOUND_ENABLED) then
                                alpha = step_bound_lambda
                            else
                                alpha = min(1.0d0, bounded_step_factor(current, du_eff, TEMP_MIN_C, TEMP_MAX_C))
                            end if
                        end if
                        current(:) = current(:) + alpha * du_eff(:)
                    else
                        max_du = maxval(abs(du))
                        ! Step limiter: tighter limit near T_melt to prevent C-C amplification
                        block
                            real(real64) :: local_max_step
                            if (minval(abs(current)) < PICARD_PHASE_ZONE) then
                                local_max_step = PICARD_MAX_DT_STEP_PHASE
                            else
                                local_max_step = PICARD_MAX_DT_STEP
                            end if
                            if (max_du > local_max_step) then
                                alpha = local_max_step / max_du
                            else
                                alpha = 1.0d0
                            end if
                        end block
                        ! Backtracking line search: halve alpha until T stays within bounds
                        block
                            real(real64) :: alpha_ls
                            real(real64), parameter :: LS_FACTOR = 0.5d0
                            real(real64), parameter :: LS_MIN = 1.0d-4
                            integer(int32) :: ls_iter
                            alpha_ls = alpha
                            do ls_iter = 1, 20
                                if (minval(current + alpha_ls * du) >= TEMP_MIN_C .and. &
                                    maxval(current + alpha_ls * du) <= TEMP_MAX_C) exit
                                alpha_ls = alpha_ls * LS_FACTOR
                                if (alpha_ls < LS_MIN) then
                                    alpha_ls = LS_MIN
                                    exit
                                end if
                            end do
                            alpha = alpha_ls
                        end block
                        if (.not. is_none) then
                            call self%control%compute_relaxation(PHYSICS_TYPES%THERMAL, iter, alpha * du, current)
                            call self%control%get_current_relaxation(PHYSICS_TYPES%THERMAL, relaxation_factor)
                        else
                            relaxation_factor = alpha
                            current(:) = current(:) + alpha * du(:)
                        end if
                    end if
                end if

                if (minval(current(:)) < TEMP_MIN_C .or. maxval(current(:)) > TEMP_MAX_C) then
                    current(:) = current_prev(:)
                    if (allocated(du) .and. size(du) > 0) then
                        alpha = bounded_step_factor(current, du, TEMP_MIN_C, TEMP_MAX_C)
                        current(:) = current(:) + alpha * du(:)
                    end if
                end if
                call self%temperature%set_delta(current(:) - current_prev(:))

                call deallocate_array(current_prev)
            end if

            call self%calc_gradient_temperature()
            call self%temperature%compute_time_derivative(bdf_coeffs, bdf_order)

            call deallocate_array(du)

            if (self%control%is_staggered()) then
                block
                    real(real64), pointer :: du_raw(:) => null()
                    du_raw => self%du%get_data(PHYSICS_TYPES%THERMAL%ID)
                    if (associated(du_raw)) du_raw(:) = 0.0d0
                    nullify (du_raw)
                end block
            end if
        end if

        if (self%is_active_hydraulic()) then
            call self%get_variable_increment(PHYSICS_TYPES%HYDRAULIC, du)
            call self%pressure%get_current(current)
            if (associated(current)) then
                call allocate_array(current_prev, size(current))
                current_prev(:) = current(:)

                if (allocated(du) .and. size(du) > 0) then
                    if (is_conserved_mode) then
                        ! Same coupled Picard damping factor as the thermal branch.
                        relaxation_factor = self%control%get_conserved_relaxation()
                        if (present(step_scale)) relaxation_factor = relaxation_factor * step_scale
                        call allocate_array(du_eff, size(du))
                        if (aa_active .and. allocated(aa_step_P) .and. size(aa_step_P) == size(du)) then
                            du_eff(:) = aa_step_P(:)
                            alpha = aa_alpha_joint
                        else
                            du_eff(:) = relaxation_factor * du(:)
                            if (STEP_BOUND_ENABLED) then
                                alpha = step_bound_lambda
                            else
                                alpha = min(1.0d0, bounded_step_factor(current, du_eff, PRESS_MIN_PA, PRESS_MAX_PA))
                            end if
                        end if
                        current(:) = current(:) + alpha * du_eff(:)
                    else
                        max_du = maxval(abs(du))
                        ! Step limiter: prevent large pressure updates regardless of solver mode
                        if (max_du > PICARD_MAX_DP_STEP) then
                            alpha = PICARD_MAX_DP_STEP / max_du
                        else
                            alpha = 1.0d0
                        end if
                        ! Backtracking line search: halve alpha until P stays within bounds
                        block
                            real(real64) :: alpha_ls
                            real(real64), parameter :: LS_FACTOR = 0.5d0
                            real(real64), parameter :: LS_MIN = 1.0d-4
                            integer(int32) :: ls_iter
                            alpha_ls = alpha
                            do ls_iter = 1, 20
                                if (minval(current + alpha_ls * du) >= PRESS_MIN_PA .and. &
                                    maxval(current + alpha_ls * du) <= PRESS_MAX_PA) exit
                                alpha_ls = alpha_ls * LS_FACTOR
                                if (alpha_ls < LS_MIN) then
                                    alpha_ls = LS_MIN
                                    exit
                                end if
                            end do
                            alpha = alpha_ls
                        end block
                        if (.not. is_none) then
                            call self%control%compute_relaxation(PHYSICS_TYPES%HYDRAULIC, iter, alpha * du, current)
                            call self%control%get_current_relaxation(PHYSICS_TYPES%HYDRAULIC, relaxation_factor)
                        else
                            relaxation_factor = alpha
                            current(:) = current(:) + alpha * du(:)
                        end if
                    end if
                end if

                if (minval(current(:)) < PRESS_MIN_PA .or. maxval(current(:)) > PRESS_MAX_PA) then
                    current(:) = current_prev(:)
                    if (allocated(du) .and. size(du) > 0) then
                        alpha = bounded_step_factor(current, du, PRESS_MIN_PA, PRESS_MAX_PA)
                        current(:) = current(:) + alpha * du(:)
                    end if
                end if
                call self%pressure%set_delta(current(:) - current_prev(:))

                call deallocate_array(current_prev)
            end if

            call self%calc_gradient_pressure()
            call self%pressure%compute_time_derivative(bdf_coeffs, bdf_order)

            call deallocate_array(du)

            if (self%control%is_staggered()) then
                block
                    real(real64), pointer :: du_raw(:) => null()
                    du_raw => self%du%get_data(PHYSICS_TYPES%HYDRAULIC%ID)
                    if (associated(du_raw)) du_raw(:) = 0.0d0
                    nullify (du_raw)
                end block
            end if
        end if

        if (PHASE_ONSET_TEMPERATURE_RESET .and. self%is_active_thermal() .and. self%is_active_hydraulic() .and. &
            allocated(temperature_before_phase)) then
            nullify (current)
            call self%temperature%get_current(current)
            if (associated(current)) then
                call self%apply_phase_change_temperature_correction(temperature_before_phase, current)
                call self%temperature%set_delta(current(:) - temperature_before_phase(:))
                call self%calc_gradient_temperature()
                call self%temperature%compute_time_derivative(bdf_coeffs, bdf_order)
            end if
            call deallocate_array(temperature_before_phase)
        end if

        call self%update_nodal_phases()

        call self%control%profiler_stop(PROFILER_TYPES%SETUP)

    contains
        pure function bounded_step_factor(vec, step, lower, upper) result(factor)
            implicit none
            real(real64), intent(in) :: vec(:)
            real(real64), intent(in) :: step(:)
            real(real64), intent(in) :: lower, upper
            real(real64) :: factor
            real(real64) :: candidate
            integer(int32) :: i, n

            factor = 1.0d0
            n = min(size(vec), size(step))
            do i = 1, n
                if (step(i) > 0.0d0) then
                    candidate = (upper - vec(i)) / step(i)
                    factor = min(factor, candidate)
                else if (step(i) < 0.0d0) then
                    candidate = (lower - vec(i)) / step(i)
                    factor = min(factor, candidate)
                end if
            end do
            factor = max(0.0d0, min(1.0d0, factor))
        end function bounded_step_factor

        !> Form one bounded, joint AA(1) step for temperature and pressure.
        !>
        !> The least-squares safeguard acts on the predicted fixed-point
        !> residual, so a period-two sequence is eligible for acceleration.
        !> A separate step-growth guard catches large secant extrapolations.
        !> One feasibility factor is applied to both fields to preserve the
        !> coupled AA direction at physical bounds.
        subroutine prepare_coupled_aa_step(self, omega, gamma, step_T, step_P, alpha_joint, active)
            implicit none
            class(type_ftcms), intent(inout) :: self
            real(real64), intent(in) :: omega
            real(real64), intent(inout) :: gamma
            real(real64), allocatable, intent(inout) :: step_T(:), step_P(:)
            real(real64), intent(inout) :: alpha_joint
            logical, intent(inout) :: active

            real(real64), allocatable :: g_T(:), g_P(:)
            real(real64), pointer, contiguous :: current_T(:), current_P(:)
            real(real64) :: base_norm, candidate_norm
            logical :: aa_usable

            gamma = 0.0d0
            alpha_joint = 1.0d0
            aa_usable = .false.
            nullify (current_T, current_P)
            call self%get_variable_increment(PHYSICS_TYPES%THERMAL, g_T)
            call self%get_variable_increment(PHYSICS_TYPES%HYDRAULIC, g_P)
            if (.not. (allocated(g_T) .and. allocated(g_P))) then
                active = .false.
                return
            end if
            call self%temperature%get_current(current_T)
            call self%pressure%get_current(current_P)
            if (.not. (associated(current_T) .and. associated(current_P))) then
                active = .false.
                return
            end if
            if (size(current_T) /= size(g_T) .or. size(current_P) /= size(g_P)) then
                active = .false.
                return
            end if

            if (self%aa_has_prev .and. allocated(self%aa_duT_prev) .and. allocated(self%aa_duP_prev)) then
                call compute_coupled_aa1_coefficient(g_T, g_P, self%aa_duT_prev, self%aa_duP_prev, &
                                                     AA_WEIGHT_P, AA_GAMMA_MAX, gamma, aa_usable)
            end if

            allocate (step_T(size(g_T)), step_P(size(g_P)))
            step_T(:) = omega * g_T(:)
            step_P(:) = omega * g_P(:)
            if (aa_usable .and. self%aa_has_prev .and. &
                allocated(self%aa_T_prev) .and. allocated(self%aa_P_prev)) then
                if (size(self%aa_T_prev) == size(current_T) .and. &
                    size(self%aa_P_prev) == size(current_P)) then
                    step_T(:) = step_T(:) - gamma * (current_T(:) - self%aa_T_prev(:) + &
                                                      omega * (g_T(:) - self%aa_duT_prev(:)))
                    step_P(:) = step_P(:) - gamma * (current_P(:) - self%aa_P_prev(:) + &
                                                      omega * (g_P(:) - self%aa_duP_prev(:)))
                else
                    aa_usable = .false.
                end if
            end if

            base_norm = abs(omega) * coupled_weighted_norm(g_T, g_P, AA_WEIGHT_P)
            candidate_norm = coupled_weighted_norm(step_T, step_P, AA_WEIGHT_P)
            if (.not. aa_usable .or. .not. (candidate_norm == candidate_norm) .or. &
                candidate_norm > AA_STEP_GROWTH_MAX * max(base_norm, tiny(1.0d0))) then
                gamma = 0.0d0
                step_T(:) = omega * g_T(:)
                step_P(:) = omega * g_P(:)
            end if

            alpha_joint = min(bounded_step_factor(current_T, step_T, TEMP_MIN_C, TEMP_MAX_C), &
                              bounded_step_factor(current_P, step_P, PRESS_MIN_PA, PRESS_MAX_PA))

            call copy_into(self%aa_T_prev, current_T)
            call copy_into(self%aa_P_prev, current_P)
            call copy_into(self%aa_duT_prev, g_T)
            call copy_into(self%aa_duP_prev, g_P)
            self%aa_gnorm_prev = coupled_weighted_norm(g_T, g_P, AA_WEIGHT_P)
            self%aa_has_prev = .true.
        end subroutine prepare_coupled_aa_step

        !> (Re)allocate dst to the shape of src and copy.
        subroutine copy_into(dst, src)
            implicit none
            real(real64), allocatable, intent(inout) :: dst(:)
            real(real64), intent(in) :: src(:)

            if (allocated(dst)) then
                if (size(dst) /= size(src)) deallocate (dst)
            end if
            if (.not. allocated(dst)) allocate (dst(size(src)))
            dst(:) = src(:)
        end subroutine copy_into

    end subroutine reflect_variables_ftcms

    !> Reset a freezing-onset crossing to its pressure-dependent critical
    !> temperature, following Hansson et al. (2004).
    module subroutine apply_phase_change_temperature_correction_ftcms(self, T_old, T_new)
        implicit none
        class(type_ftcms), intent(inout) :: self
        real(real64), intent(in) :: T_old(:)
        real(real64), intent(inout) :: T_new(:)

        integer(int32) :: node_id, num_nodes, row_start, repr_elem
        type(type_state) :: state
        real(real64) :: pressure, rho_water, critical_temperature

        call self%domain%get_num_nodes(num_nodes)
        if (num_nodes <= 0) return

        do node_id = 1, min(num_nodes, size(T_old), size(T_new))
            ! The pressure-dependent freezing temperature cannot exceed
            ! 0 degC. A node that is still non-negative, or did not cool in
            ! this Picard update, cannot cross the freezing active set.
            if (T_new(node_id) >= 0.0d0 .or. T_new(node_id) >= T_old(node_id)) cycle
            row_start = self%node_material_table%ptr(node_id)
            if (row_start >= self%node_material_table%ptr(node_id + 1)) cycle
            repr_elem = self%node_material_table%repr_element(row_start)
            call self%set_state(node_id, repr_elem, state, calc_physics=.false., include_fluxes=.false.)
            call state%pressure%get(pressure)
            call self%thermal%calc_density_water(state, rho_water)
            call calc_T_high_celsius(pressure, rho_water, critical_temperature)

            if (T_old(node_id) > critical_temperature .and. &
                T_new(node_id) < critical_temperature .and. &
                .not. self%phase_onset_reset(node_id)) then
                T_new(node_id) = critical_temperature
                self%phase_onset_reset(node_id) = .true.
            end if
        end do
    end subroutine apply_phase_change_temperature_correction_ftcms

    !> Implementation strategy: two-pass CSR construction, mirroring
    !> domain_adjacency_node_element's initialize (count-then-fill).
    !>
    !> Pass 1 walks each node's adjacent-element list (ascending element-id
    !> order, from the CSR-backed node->element adjacency) and counts the
    !> number of distinct materials touching that node, using a small
    !> per-node scratch buffer sized to the node's degree. A prefix sum over
    !> the counts yields the CSR row pointers.
    !>
    !> Pass 2 repeats the walk, this time writing each entry's material id,
    !> accumulating measure_sum via domain%calc_measure, and overwriting
    !> repr_element on every repeat visit of a material -- since the walk is
    !> ascending, the value left after the loop is the highest-index element
    !> of that material. Pass 2 uses the (now allocated) material_id row
    !> itself as the scratch buffer for the linear "have we seen this
    !> material yet" search, avoiding a second scratch allocation.
    !>
    !> Parallelization: none (serial, build-once cost; the mesh is static so
    !> this pays for itself many times over across nonlinear iterations).
    !> Memory: O(N_nd + E) as documented at the interface.
    module subroutine initialize_node_material_table(self, domain)
        implicit none
        class(type_node_material_table), intent(inout) :: self
        type(type_domain), intent(in) :: domain

        integer(int32) :: num_nodes, i_node, j, k, n_neighbors, i_elem, material_id, n_distinct
        integer(int32) :: cum_sum, total_entries, row_start
        integer(int32), pointer, contiguous :: element_list(:)
        integer(int32), allocatable :: local_material(:)
        real(real64) :: measure

        call self%destroy()
        nullify (element_list)

        call domain%get_num_nodes(num_nodes)
        if (num_nodes <= 0) return

        self%num_nodes = num_nodes
        allocate (self%ptr(num_nodes + 1))
        self%ptr = 0

        ! --- Pass 1: count distinct materials per node (row sizes) ---
        do i_node = 1, num_nodes
            call domain%element_adjacency%get_list(i_node, element_list)
            n_distinct = 0
            if (associated(element_list)) then
                n_neighbors = size(element_list)
                if (allocated(local_material)) deallocate (local_material)
                allocate (local_material(n_neighbors))
                do j = 1, n_neighbors
                    call domain%get_material_id(element_list(j), material_id)
                    k = 1
                    do while (k <= n_distinct)
                        if (local_material(k) == material_id) exit
                        k = k + 1
                    end do
                    if (k > n_distinct) then
                        n_distinct = n_distinct + 1
                        local_material(n_distinct) = material_id
                    end if
                end do
            end if
            self%ptr(i_node + 1) = n_distinct
            nullify (element_list)
        end do
        if (allocated(local_material)) deallocate (local_material)

        ! --- Prefix sum: distinct-material counts -> CSR row pointers ---
        self%ptr(1) = 1
        cum_sum = 1
        do i_node = 1, num_nodes
            cum_sum = cum_sum + self%ptr(i_node + 1)
            self%ptr(i_node + 1) = cum_sum
        end do
        total_entries = cum_sum - 1

        allocate (self%material_id(max(total_entries, 1)))
        allocate (self%repr_element(max(total_entries, 1)))
        allocate (self%measure_sum(max(total_entries, 1)))
        self%material_id = 0
        self%repr_element = 0
        self%measure_sum = 0.0d0

        ! --- Pass 2: fill entries ---
        do i_node = 1, num_nodes
            call domain%element_adjacency%get_list(i_node, element_list)
            if (.not. associated(element_list)) cycle
            n_neighbors = size(element_list)
            row_start = self%ptr(i_node)
            n_distinct = 0
            do j = 1, n_neighbors
                i_elem = element_list(j)
                call domain%get_material_id(i_elem, material_id)
                call domain%calc_measure(i_elem, measure)

                k = 1
                do while (k <= n_distinct)
                    if (self%material_id(row_start + k - 1) == material_id) exit
                    k = k + 1
                end do
                if (k > n_distinct) then
                    n_distinct = n_distinct + 1
                    self%material_id(row_start + n_distinct - 1) = material_id
                end if
                ! Ascending element order => the last overwrite is the
                ! highest-index adjacent element of this material.
                self%repr_element(row_start + k - 1) = i_elem
                self%measure_sum(row_start + k - 1) = self%measure_sum(row_start + k - 1) + measure
            end do
            nullify (element_list)
        end do

    end subroutine initialize_node_material_table

    module subroutine destroy_node_material_table(self)
        implicit none
        class(type_node_material_table), intent(inout) :: self

        if (allocated(self%ptr)) deallocate (self%ptr)
        if (allocated(self%material_id)) deallocate (self%material_id)
        if (allocated(self%repr_element)) deallocate (self%repr_element)
        if (allocated(self%measure_sum)) deallocate (self%measure_sum)
        self%num_nodes = 0
    end subroutine destroy_node_material_table

    !> Implementation strategy: node-major loop over the precomputed
    !> node_material_table (see type_node_material_table), replacing the
    !> original element-major loop that visited every node once per adjacent
    !> element and let the last (highest-index) visiting element's
    !> set_state() result win.
    !>
    !> - Nodes with no adjacent elements (empty table row) are skipped,
    !>   leaving their Qw/Qi/Qa/Qv untouched -- identical to the original
    !>   loop, which never wrote a node with no incident element.
    !> - Nodes touching a single material (m == 1, the common case) call
    !>   set_state exactly once, at the table's representative element for
    !>   that material (the highest-index adjacent element overall in this
    !>   case). This exactly reproduces the original "last element wins"
    !>   result while visiting the node only once.
    !> - Nodes touching multiple materials (m > 1) call set_state once per
    !>   distinct material and combine the results with a
    !>   material-summed-measure weighted average. The is_set guard is
    !>   applied per component and per material: a component's weighted
    !>   average only includes materials where the constitutive update
    !>   actually produced that component, and the nodal value is written
    !>   only if at least one material produced it -- generalizing the
    !>   original single-element is_set guard to the multi-material case.
    !>
    !> Cost: O(sum of distinct materials per node) set_state calls, versus
    !> O(sum of node degrees) in the original loop (identical when every
    !> node touches a single material, e.g. this project's mesh).
    !>
    !> set_state is called in phases-only mode (include_fluxes=.false.):
    !> only the phase contents Qw/Qi/Qa/Qv are read from the state here, so
    !> the Darcy flux evaluation is skipped. Assembly and output build their
    !> own states with fluxes included and are unaffected.
    !>
    !> Parallelization: the node loop is OpenMP-parallel with one
    !> thread-local scratch state per thread (same scheme as
    !> update_variables_ftcms). Each iteration reads only global nodal
    !> arrays plus the static table and writes only its own node_id entry of
    !> Qw/Qi/Qa/Qv, so the result is independent of the schedule and
    !> bit-identical to the serial loop.
    module subroutine update_nodal_phases_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        integer(int32) :: node_id, num_nodes, row_start, row_end, m, k, repr_elem
        type(type_state), allocatable :: states(:)
        real(real64) :: qw_val, qi_val, qa_val, qv_val
        logical :: qi_set, qw_set, qa_set, qv_set
        real(real64) :: measure
        real(real64) :: sum_qw, sum_qi, sum_qa, sum_qv
        real(real64) :: wsum_qw, wsum_qi, wsum_qa, wsum_qv
        logical :: any_qw, any_qi, any_qa, any_qv
        integer(int32) :: num_threads, tid

        call self%domain%get_num_nodes(num_nodes)
        num_threads = omp_get_max_threads()
        allocate (states(num_threads))

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED(self, num_nodes, states) &
        !$OMP PRIVATE(node_id, row_start, row_end, m, k, repr_elem, &
        !$OMP         qw_val, qi_val, qa_val, qv_val, &
        !$OMP         qi_set, qw_set, qa_set, qv_set, measure, &
        !$OMP         sum_qw, sum_qi, sum_qa, sum_qv, &
        !$OMP         wsum_qw, wsum_qi, wsum_qa, wsum_qv, &
        !$OMP         any_qw, any_qi, any_qa, any_qv, tid)
        tid = omp_get_thread_num() + 1
        !$OMP DO
        do node_id = 1, num_nodes
            row_start = self%node_material_table%ptr(node_id)
            row_end = self%node_material_table%ptr(node_id + 1) - 1
            m = row_end - row_start + 1
            if (m <= 0) cycle ! no adjacent elements: leave Q values untouched

            if (m == 1) then
                repr_elem = self%node_material_table%repr_element(row_start)
                call self%set_state(node_id, repr_elem, states(tid), calc_physics=.true., include_fluxes=.false.)
                qi_set = .false.; qw_set = .false.; qa_set = .false.
                qv_set = .false.
                call states(tid)%ice_content%get(qi_val, qi_set)
                call states(tid)%water_content%get(qw_val, qw_set)
                call states(tid)%air_content%get(qa_val, qa_set)
                call states(tid)%vapor_content%get(qv_val, qv_set)
                if (qi_set) call self%Qi%set_current(node_id, qi_val)
                if (qw_set) call self%Qw%set_current(node_id, qw_val)
                if (qa_set) call self%Qa%set_current(node_id, qa_val)
                if (qv_set) call self%Qv%set_current(node_id, qv_val)
            else
                sum_qw = 0.0d0; sum_qi = 0.0d0; sum_qa = 0.0d0; sum_qv = 0.0d0
                wsum_qw = 0.0d0; wsum_qi = 0.0d0; wsum_qa = 0.0d0; wsum_qv = 0.0d0
                any_qw = .false.; any_qi = .false.; any_qa = .false.; any_qv = .false.

                do k = row_start, row_end
                    repr_elem = self%node_material_table%repr_element(k)
                    measure = self%node_material_table%measure_sum(k)
                    call self%set_state(node_id, repr_elem, states(tid), calc_physics=.true., include_fluxes=.false.)
                    qi_set = .false.; qw_set = .false.; qa_set = .false.
                    qv_set = .false.
                    call states(tid)%ice_content%get(qi_val, qi_set)
                    call states(tid)%water_content%get(qw_val, qw_set)
                    call states(tid)%air_content%get(qa_val, qa_set)
                    call states(tid)%vapor_content%get(qv_val, qv_set)

                    if (qw_set) then
                        sum_qw = sum_qw + qw_val * measure
                        wsum_qw = wsum_qw + measure
                        any_qw = .true.
                    end if
                    if (qi_set) then
                        sum_qi = sum_qi + qi_val * measure
                        wsum_qi = wsum_qi + measure
                        any_qi = .true.
                    end if
                    if (qa_set) then
                        sum_qa = sum_qa + qa_val * measure
                        wsum_qa = wsum_qa + measure
                        any_qa = .true.
                    end if
                    if (qv_set) then
                        sum_qv = sum_qv + qv_val * measure
                        wsum_qv = wsum_qv + measure
                        any_qv = .true.
                    end if
                end do

                if (any_qw .and. wsum_qw > epsilon(1.0d0)) call self%Qw%set_current(node_id, sum_qw / wsum_qw)
                if (any_qi .and. wsum_qi > epsilon(1.0d0)) call self%Qi%set_current(node_id, sum_qi / wsum_qi)
                if (any_qa .and. wsum_qa > epsilon(1.0d0)) call self%Qa%set_current(node_id, sum_qa / wsum_qa)
                if (any_qv .and. wsum_qv > epsilon(1.0d0)) call self%Qv%set_current(node_id, sum_qv / wsum_qv)
            end if
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        if (allocated(states)) deallocate (states)

    end subroutine update_nodal_phases_ftcms

    !> Update the outer nodal ice state by a local Clapeyron projection.
    !> The monolithic unknown vector remains (T, p_w); this projection is a
    !> bounded phase-transfer correction between water-conserving monolithic
    !> solves. A converged outer iteration satisfies both balances.
    module subroutine project_nodal_ice_ftcms(self, apply_update, ice_update, max_increment, increment_norm, &
                                               max_node, max_temperature, max_pressure, &
                                               max_current_ice, max_projected_ice, &
                                               max_equilibrium_error, increments, active_bounds, &
                                               apply_projected_update, projected_update_scale)
        implicit none
        class(type_ftcms), intent(inout) :: self
        logical, intent(in) :: apply_update
        real(real64), intent(inout) :: ice_update(:)
        real(real64), intent(inout) :: max_increment
        real(real64), intent(inout) :: increment_norm
        integer(int32), intent(inout) :: max_node
        real(real64), intent(inout) :: max_temperature
        real(real64), intent(inout) :: max_pressure
        real(real64), intent(inout) :: max_current_ice
        real(real64), intent(inout) :: max_projected_ice
        real(real64), intent(inout) :: max_equilibrium_error
        real(real64), allocatable, intent(inout) :: increments(:)
        integer(int32), allocatable, intent(inout) :: active_bounds(:)
        logical, intent(in), optional :: apply_projected_update
        real(real64), intent(in), optional :: projected_update_scale

        integer(int32) :: node_id, num_nodes, row_start, row_end, k, repr_elem
        integer(int32) :: num_threads, tid, active_bound, node_active_bound
        type(type_state), allocatable :: states(:)
        real(real64), allocatable :: current_ice_values(:), projected_ice_values(:)
        real(real64), allocatable :: equilibrium_errors(:), node_measures(:)
        integer(int32), allocatable :: projected_active_bounds(:)
        real(real64) :: measure, projected_ice, ice_increment, equilibrium_error, porosity, updated_ice
        real(real64) :: node_temperature, projected_scale
        real(real64) :: weighted_ice, weight_sum, current_ice, node_equilibrium_error
        real(real64) :: target_total_water, local_pressure, local_ice
        logical :: target_available, local_converged, use_local_eq, use_projected_update
        logical, allocatable :: is_hydraulic_boundary_node(:)
        integer(int32) :: i_patch, num_patches, i_bc_node, bc_idx, comp_dim, i_coord, j_coord
        type(type_boundary_patch), pointer :: bc_patch
        real(real64) :: dt_local, C_eq_val, hydraulic_diffusivity, element_length, fourier_number, S_seg
        real(real64), allocatable :: D_HH_matrix(:, :), D_HH_gas_matrix(:, :), element_coords(:, :)
        integer(int32) :: local_eq_candidates, local_eq_fires
        real(real64) :: local_eq_fo_min

        call self%domain%get_num_nodes(num_nodes)
        use_projected_update = .false.
        if (present(apply_projected_update)) use_projected_update = apply_projected_update
        projected_scale = 1.0d0
        if (present(projected_update_scale)) projected_scale = min(max(projected_update_scale, 0.0d0), 1.0d0)
        num_threads = omp_get_max_threads()
        allocate (states(num_threads))
        allocate (current_ice_values(num_nodes), projected_ice_values(num_nodes))
        allocate (equilibrium_errors(num_nodes), source=0.0d0)
        allocate (node_measures(num_nodes), source=0.0d0)
        allocate (projected_active_bounds(num_nodes), source=0_int32)
        if (allocated(increments)) deallocate (increments)
        allocate (increments(num_nodes))
        if (allocated(active_bounds)) deallocate (active_bounds)
        allocate (active_bounds(num_nodes))
        max_increment = 0.0d0

        ! Nodes touched by ANY hydraulic boundary patch (Dirichlet or flux) are
        ! excluded from the local-conserved-equilibrium fast path below,
        ! regardless of their Fourier number: a prescribed boundary value or
        ! flux can inject/remove water independent of the internal transport
        ! operator's own relaxation rate (see plan spicy-sauteeing-scroll.md,
        ! WP3 refinement discussion). Built once per call, not per node.
        allocate (is_hydraulic_boundary_node(num_nodes), source=.false.)
        dt_local = 0.0d0
        if (PHASE_USE_LOCAL_CONSERVED_EQ) then
            call self%domain%get_num_bc_patches(num_patches)
            do i_patch = 1, num_patches
                call self%domain%get_bc_patch(i_patch, bc_patch)
                call self%bc(PHYSICS_TYPES%HYDRAULIC%ID)%get_bc_index(bc_patch%entity_id, bc_idx)
                if (bc_idx < 0) cycle
                if (allocated(bc_patch%connectivity%col_ind)) then
                    do i_bc_node = 1, size(bc_patch%connectivity%col_ind)
                        is_hydraulic_boundary_node(bc_patch%connectivity%col_ind(i_bc_node)) = .true.
                    end do
                end if
            end do
            call self%control%get_dt(dt_local)
            call self%domain%get_computation_dimension(comp_dim)
        end if

        local_eq_candidates = 0
        local_eq_fires = 0
        local_eq_fo_min = huge(1.0d0)
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED(self, num_nodes, states, current_ice_values, projected_ice_values, equilibrium_errors, node_measures, &
        !$OMP        projected_active_bounds, is_hydraulic_boundary_node, dt_local, comp_dim) &
        !$OMP PRIVATE(node_id, row_start, row_end, k, repr_elem, tid, measure, porosity, node_temperature, &
        !$OMP         projected_ice, ice_increment, equilibrium_error, weighted_ice, &
        !$OMP         weight_sum, current_ice, node_equilibrium_error, active_bound, node_active_bound, &
        !$OMP         target_total_water, local_pressure, local_ice, target_available, local_converged, use_local_eq, &
        !$OMP         C_eq_val, hydraulic_diffusivity, element_length, fourier_number, S_seg, &
        !$OMP         D_HH_matrix, D_HH_gas_matrix, element_coords, i_coord, j_coord) &
        !$OMP REDUCTION(+:local_eq_candidates, local_eq_fires) REDUCTION(min:local_eq_fo_min)
        tid = omp_get_thread_num() + 1
        !$OMP DO
        do node_id = 1, num_nodes
            call self%Qi%get_current(node_id, current_ice)
            current_ice_values(node_id) = current_ice
            projected_ice_values(node_id) = current_ice
            call self%temperature%get_current(node_id, node_temperature)

            ! Exact warm-side active-set exclusion: the non-segregated
            ! Clapeyron suction is zero for T >= 0 degC, so an ice-free node
            ! has a zero return-map increment without evaluating any material
            ! model. This is only a fast path for a mathematically inactive
            ! branch; all subzero or already-frozen nodes use the full map.
            if (current_ice <= 0.0d0 .and. node_temperature >= 0.0d0) then
                projected_active_bounds(node_id) = -1
                cycle
            end if

            row_start = self%node_material_table%ptr(node_id)
            row_end = self%node_material_table%ptr(node_id + 1) - 1
            if (row_end < row_start) cycle

            call self%porosity%get_current(node_id, porosity)

            weighted_ice = 0.0d0
            weight_sum = 0.0d0
            node_equilibrium_error = 0.0d0
            node_active_bound = 2
            do k = row_start, row_end
                repr_elem = self%node_material_table%repr_element(k)
                measure = self%node_material_table%measure_sum(k)
                call self%set_state(node_id, repr_elem, states(tid), calc_physics=.true., include_fluxes=.false.)

                use_local_eq = .false.
                if (PHASE_USE_LOCAL_CONSERVED_EQ .and. .not. is_hydraulic_boundary_node(node_id)) then
                    ! Fourier number Fo = (D_HH/C_eq) * dt / h^2 at this
                    ! material's representative element: D_HH/C_eq is the
                    ! hydraulic diffusivity [m^2/s] the *inner* (fixed-ice)
                    ! monolithic solve itself uses (C_eq's dQi_dP is always
                    ! zero there - phase_systems.F90 - so this is exactly the
                    ! operator whose relaxation we are judging). h is the
                    ! minimum pairwise node spacing of the element, a safe
                    ! (smallest, not area-averaged) length scale on this
                    ! deliberately anisotropic mesh.
                    if (.not. allocated(D_HH_matrix)) then
                        allocate (D_HH_matrix(comp_dim, comp_dim))
                        allocate (D_HH_gas_matrix(comp_dim, comp_dim))
                    else if (size(D_HH_matrix, 1) /= comp_dim) then
                        deallocate (D_HH_matrix, D_HH_gas_matrix)
                        allocate (D_HH_matrix(comp_dim, comp_dim))
                        allocate (D_HH_gas_matrix(comp_dim, comp_dim))
                    end if
                    call self%hydraulic%compute_diffusion_term(self%node_material_table%material_id(k), &
                                                                states(tid), D_HH_matrix)
                    call self%hydraulic%compute_diffusion_term_gas(self%node_material_table%material_id(k), &
                                                                   states(tid), D_HH_gas_matrix)
                    ! The Fourier gate judges the total hydraulic relaxation.
                    D_HH_matrix(:, :) = D_HH_matrix(:, :) + D_HH_gas_matrix(:, :)
                    call self%hydraulic%compute_C_eq(self%node_material_table%material_id(k), states(tid), C_eq_val)
                    call self%domain%get_fe_coordinate(repr_elem, element_coords)
                    element_length = huge(1.0d0)
                    do i_coord = 1, size(element_coords, 2) - 1
                        do j_coord = i_coord + 1, size(element_coords, 2)
                            element_length = min(element_length, &
                                norm2(element_coords(:, i_coord) - element_coords(:, j_coord)))
                        end do
                    end do
                    fourier_number = huge(1.0d0)
                    if (C_eq_val > tiny(1.0d0) .and. element_length > tiny(1.0d0)) then
                        hydraulic_diffusivity = D_HH_matrix(1, 1) / C_eq_val
                        fourier_number = hydraulic_diffusivity * dt_local / element_length**2
                    end if
                    if (current_ice > 0.05d0) then
                        !$OMP CRITICAL (fo_debug)
                        write (*, '(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
                            '   [FO_DEBUG] node=', node_id, ' Qi=', current_ice, ' D_HH=', D_HH_matrix(1, 1), &
                            ' C_eq=', C_eq_val, ' D=', hydraulic_diffusivity, ' h=', element_length, &
                            ' Fo=', fourier_number
                        !$OMP END CRITICAL (fo_debug)
                    end if
                    local_eq_candidates = local_eq_candidates + 1
                    local_eq_fo_min = min(local_eq_fo_min, fourier_number)

                    S_seg = 0.0d0
                    call self%hydraulic%calc_segregation_sink(self%node_material_table%material_id(k), &
                                                               states(tid), dt_local, S_seg)

                    if (fourier_number < PHASE_LOCAL_EQ_FOURIER_GATE .and. &
                        abs(S_seg * dt_local) < PHASE_LOCAL_EQ_TRANSPORT_EPS * max(porosity, 1.0d-3)) then
                        call self%thermal%calc_conserved_target(self%node_material_table%material_id(k), &
                                                                 states(tid), target_total_water, target_available)
                        if (target_available) then
                            call self%thermal%solve_local_conserved_equilibrium( &
                                self%node_material_table%material_id(k), states(tid), target_total_water, &
                                local_pressure, local_ice, local_converged)
                            if (local_converged) then
                                projected_ice = local_ice
                                ice_increment = projected_ice - current_ice
                                equilibrium_error = 0.0d0
                                active_bound = 0
                                use_local_eq = .true.
                                local_eq_fires = local_eq_fires + 1
                            end if
                        end if
                    end if
                end if

                if (.not. use_local_eq) then
                    call self%thermal%project_ice_content(self%node_material_table%material_id(k), &
                                                          states(tid), projected_ice, ice_increment, &
                                                          equilibrium_error, active_bound)
                end if
                weighted_ice = weighted_ice + measure * projected_ice
                weight_sum = weight_sum + measure
                node_equilibrium_error = max(node_equilibrium_error, equilibrium_error)
                if (node_active_bound == 2) then
                    node_active_bound = active_bound
                else if (node_active_bound /= active_bound) then
                    node_active_bound = 0
                end if
            end do

            if (weight_sum > epsilon(1.0d0)) then
                projected_ice_values(node_id) = weighted_ice / weight_sum
                equilibrium_errors(node_id) = node_equilibrium_error
                node_measures(node_id) = weight_sum
                projected_active_bounds(node_id) = node_active_bound
            end if
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        if (PHASE_USE_LOCAL_CONSERVED_EQ) then
            write (*, '(A,I0,A,I0,A,ES10.3)') '   [LOCAL_EQ_DEBUG] candidates=', local_eq_candidates, &
                ' fires=', local_eq_fires, ' min_Fo=', local_eq_fo_min
        end if

        increments = projected_ice_values - current_ice_values
        active_bounds = projected_active_bounds
        max_node = maxloc(abs(increments), dim=1)
        max_increment = abs(increments(max_node))
        increment_norm = max_increment
        if (sum(node_measures) > tiny(1.0d0)) then
            increment_norm = sqrt(dot_product(node_measures, increments**2) / sum(node_measures))
        end if
        max_equilibrium_error = maxval(equilibrium_errors)
        max_current_ice = current_ice_values(max_node)
        max_projected_ice = projected_ice_values(max_node)
        call self%temperature%get_current(max_node, max_temperature)
        ! A thermal-only run never allocates the pressure field.
        max_pressure = 0.0d0
        if (self%is_active_hydraulic()) call self%pressure%get_current(max_node, max_pressure)

        if (apply_update) then
            do node_id = 1, num_nodes
                call self%porosity%get_current(node_id, porosity)
                if (use_projected_update) then
                    updated_ice = current_ice_values(node_id) + projected_scale * &
                                  (projected_ice_values(node_id) - current_ice_values(node_id))
                else
                    updated_ice = current_ice_values(node_id) + ice_update(node_id)
                end if
                updated_ice = min(max(updated_ice, 0.0d0), porosity)
                ! AA(1) must retain the step that was actually applied after
                ! projection onto the admissible ice-content interval.
                ice_update(node_id) = updated_ice - current_ice_values(node_id)
                call self%Qi%set_current(node_id, updated_ice)
            end do
        end if

        deallocate (states, current_ice_values, projected_ice_values, equilibrium_errors, node_measures, &
                    projected_active_bounds)
    end subroutine project_nodal_ice_ftcms

    !> Evaluate per-node conserved quantities at the current iterate.
    !>
    !> Mathematical definition:
    !> - enthalpy(j)  = volumetric enthalpy density H_j(T_j, p_{w,j}) [J/m3]
    !> - density(j)   = pore-water effective density rho_eff,j         [kg/m3]
    !>
    !> Node-major loop over the precomputed node_material_table (see
    !> type_node_material_table), replacing the original element-major loop
    !> with a processed mask (first visiting element's material wins). The
    !> nodal state is rebuilt with set_state(calc_physics=.true.) so the
    !> phase contents are consistent with the current (T, p_w) before
    !> evaluating H and rho_eff. Fluxes are skipped (include_fluxes=.false.):
    !> the enthalpy/effective-density evaluations read only T, porosity, the
    !> phase contents, and (T, p)-dependent densities/heats -- never the
    !> Darcy fluxes or gradients.
    !>
    !> - Nodes with no adjacent elements keep the pre-zeroed H and rho_eff,
    !>   identical to the original loop.
    !> - m == 1 (single distinct material): one evaluation at the table's
    !>   representative element. The evaluated material is the node's only
    !>   material, so this reproduces the original first-element-wins result
    !>   exactly (the representative element's index does not matter beyond
    !>   its material id).
    !> - m > 1: evaluate once per distinct material and combine H and
    !>   rho_eff by material-summed-measure weighted average -- the same
    !>   multi-material convention as update_nodal_phases (evaluate per
    !>   material first, then average), instead of the original "material of
    !>   the lowest-index adjacent element" pick.
    !>
    !> Inactive physics yields a zero field (contributes nothing to the
    !> weighted norm). Parallelization: OpenMP over nodes with thread-local
    !> scratch states; each iteration writes only its own node_id entries,
    !> so the result is schedule-independent. Cost: O(sum of distinct
    !> materials per node) constitutive evaluations.
    module subroutine compute_nodal_conserved_ftcms(self, enthalpy, density, dH_dT, drho_dp)
        implicit none
        class(type_ftcms), intent(inout) :: self
        real(real64), allocatable, intent(inout) :: enthalpy(:)
        real(real64), allocatable, intent(inout) :: density(:)
        real(real64), allocatable, intent(inout), optional :: dH_dT(:)
        real(real64), allocatable, intent(inout), optional :: drho_dp(:)

        integer(int32) :: node_id, material_id, n_nodes
        integer(int32) :: row_start, row_end, m, k, repr_elem
        type(type_state), allocatable :: states(:)
        logical :: active_thermal, active_hydraulic
        logical :: want_dH_dT, want_drho_dp
        real(real64) :: H_j, rho_j, measure
        real(real64) :: dHdT_j, drhodp_j, rho_w_j
        real(real64) :: sum_H, sum_rho, sum_dHdT, sum_drhodp, wsum
        integer(int32) :: num_threads, tid
        ! Always-allocated scratch for the OpenMP region: OPTIONAL dummy
        ! arguments must not be referenced when absent, including implicitly
        ! through an OMP data-sharing clause, so the sensitivities are
        ! accumulated here and only copied into the caller's (optional)
        ! arrays afterward.
        real(real64), allocatable :: dH_dT_local(:), drho_dp_local(:)

        call self%domain%get_num_nodes(n_nodes)

        if (allocated(enthalpy)) deallocate (enthalpy)
        if (allocated(density)) deallocate (density)
        allocate (enthalpy(max(n_nodes, 1)))
        allocate (density(max(n_nodes, 1)))
        enthalpy = 0.0d0
        density = 0.0d0

        want_dH_dT = present(dH_dT)
        want_drho_dp = present(drho_dp)
        allocate (dH_dT_local(max(n_nodes, 1)), source=0.0d0)
        allocate (drho_dp_local(max(n_nodes, 1)), source=0.0d0)

        if (n_nodes <= 0) then
            if (want_dH_dT) call allocate_array(dH_dT, dH_dT_local)
            if (want_drho_dp) call allocate_array(drho_dp, drho_dp_local)
            return
        end if

        active_thermal = self%is_active_thermal()
        active_hydraulic = self%is_active_hydraulic()

        num_threads = omp_get_max_threads()
        allocate (states(num_threads))

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED(self, n_nodes, states, enthalpy, density, dH_dT_local, drho_dp_local, &
        !$OMP        active_thermal, active_hydraulic, want_dH_dT, want_drho_dp) &
        !$OMP PRIVATE(node_id, material_id, row_start, row_end, m, k, repr_elem, &
        !$OMP         H_j, rho_j, dHdT_j, drhodp_j, rho_w_j, &
        !$OMP         measure, sum_H, sum_rho, sum_dHdT, sum_drhodp, wsum, tid)
        tid = omp_get_thread_num() + 1
        !$OMP DO
        do node_id = 1, n_nodes
            row_start = self%node_material_table%ptr(node_id)
            row_end = self%node_material_table%ptr(node_id + 1) - 1
            m = row_end - row_start + 1
            if (m <= 0) cycle ! no adjacent elements: keep the pre-zeroed values

            if (m == 1) then
                repr_elem = self%node_material_table%repr_element(row_start)
                material_id = self%node_material_table%material_id(row_start)
                call self%set_state(node_id, repr_elem, states(tid), calc_physics=.true., include_fluxes=.false.)
                H_j = 0.0d0
                rho_j = 0.0d0
                if (active_thermal) call self%thermal%calc_enthalpy_density(material_id, states(tid), H_j)
                if (active_hydraulic) then
                    call self%hydraulic%calc_effective_density_value(material_id, states(tid), rho_j)
                end if
                enthalpy(node_id) = H_j
                density(node_id) = rho_j
                if (want_dH_dT .and. active_thermal) then
                    dHdT_j = 0.0d0
                    call self%thermal%compute_mass_term(material_id, states(tid), dHdT_j)
                    dH_dT_local(node_id) = dHdT_j
                end if
                if (want_drho_dp .and. active_hydraulic) then
                    ! d(rho_eq)/dp = rho_w * C_eq, where rho_eq = rho_w*Theta +
                    ! (rho_std-rho_w)*compressive_storage (exact identity, see
                    ! calc_effective_density_value_hydraulic and Theta_j in
                    ! compute_transient_term_mixed_hydraulic) and C_eq =
                    ! d(Theta)/dp. The dropped (rho_std-rho_w)*d(compressive_
                    ! storage)/dp term is at most O(1 kg/m3)*O(1e-7 1/Pa),
                    ! five to six orders below rho_w*C_eq for the retention
                    ! curves this code targets, and is therefore not evaluated.
                    drhodp_j = 0.0d0
                    rho_w_j = 0.0d0
                    call self%hydraulic%compute_C_eq(material_id, states(tid), drhodp_j)
                    call self%thermal%calc_density_water(states(tid), rho_w_j)
                    drho_dp_local(node_id) = rho_w_j * drhodp_j
                end if
            else
                sum_H = 0.0d0
                sum_rho = 0.0d0
                sum_dHdT = 0.0d0
                sum_drhodp = 0.0d0
                wsum = 0.0d0
                do k = row_start, row_end
                    repr_elem = self%node_material_table%repr_element(k)
                    material_id = self%node_material_table%material_id(k)
                    measure = self%node_material_table%measure_sum(k)
                    call self%set_state(node_id, repr_elem, states(tid), calc_physics=.true., include_fluxes=.false.)
                    H_j = 0.0d0
                    rho_j = 0.0d0
                    if (active_thermal) call self%thermal%calc_enthalpy_density(material_id, states(tid), H_j)
                    if (active_hydraulic) then
                        call self%hydraulic%calc_effective_density_value(material_id, states(tid), rho_j)
                    end if
                    sum_H = sum_H + H_j * measure
                    sum_rho = sum_rho + rho_j * measure
                    if (want_dH_dT .and. active_thermal) then
                        dHdT_j = 0.0d0
                        call self%thermal%compute_mass_term(material_id, states(tid), dHdT_j)
                        sum_dHdT = sum_dHdT + dHdT_j * measure
                    end if
                    if (want_drho_dp .and. active_hydraulic) then
                        drhodp_j = 0.0d0
                        rho_w_j = 0.0d0
                        call self%hydraulic%compute_C_eq(material_id, states(tid), drhodp_j)
                        call self%thermal%calc_density_water(states(tid), rho_w_j)
                        sum_drhodp = sum_drhodp + (rho_w_j * drhodp_j) * measure
                    end if
                    wsum = wsum + measure
                end do
                ! Degenerate zero-measure rows leave the pre-zeroed values.
                if (wsum > epsilon(1.0d0)) then
                    enthalpy(node_id) = sum_H / wsum
                    density(node_id) = sum_rho / wsum
                    if (want_dH_dT .and. active_thermal) dH_dT_local(node_id) = sum_dHdT / wsum
                    if (want_drho_dp .and. active_hydraulic) drho_dp_local(node_id) = sum_drhodp / wsum
                end if
            end if
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        if (allocated(states)) deallocate (states)

        if (want_dH_dT) call allocate_array(dH_dT, dH_dT_local)
        if (want_drho_dp) call allocate_array(drho_dp, drho_dp_local)
    end subroutine compute_nodal_conserved_ftcms

    !> See the interface for the mathematical definition. Cost: O(N_dof) per step.
    module function compute_lte_error_ftcms(self, order_used) result(error_rel)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(inout), optional :: order_used
        real(real64) :: error_rel

        real(real64) :: dt_n, e_thermal, e_hydraulic
        integer(int32) :: bdf_order
        !> Nodal enthalpy [J/m3] and effective density [kg/m3]: the quantities
        !> the error is measured on. See the call sites below for why.
        real(real64), allocatable :: enthalpy(:), density(:)

        error_rel = -1.0d0
        if (present(order_used)) order_used = 1
        call self%control%get_dt(dt_n)
        if (dt_n <= 0.0d0) return
        if (.not. self%lte_has_state .or. .not. self%lte_has_derivative) return
        if (dt_n + self%lte_prev_dt <= 0.0d0) return

        ! The estimator must match the integrator's order: reusing the BDF1
        ! defect under BDF2 silently under-estimates the error by O(h), so the
        ! controller would accept steps that violate the tolerance. Orders above
        ! the implemented set are refused rather than approximated.
        call self%control%get_bdf_coeffs(bdf_order=bdf_order)
        if (bdf_order > LTE_MAX_SUPPORTED_ORDER) then
            write (*, '(A,I0,A,I0,A)') '   [LTE] error: no local-error estimator for BDF order ', &
                bdf_order, ' (implemented up to ', LTE_MAX_SUPPORTED_ORDER, ')'
            error stop 'compute_lte_error: unsupported BDF order for error control.'
        end if
        ! BDF2 needs the previous second difference; until it exists (the first
        ! two accepted steps) fall back to the BDF1 estimate, which is a valid
        ! bound because the run is still first order there anyway.
        if (bdf_order >= 2 .and. .not. self%lte_has_second_difference) bdf_order = 1
        if (bdf_order >= 2 .and. dt_n + self%lte_prev_span <= 0.0d0) bdf_order = 1
        if (present(order_used)) order_used = bdf_order

        e_thermal = -1.0d0
        e_hydraulic = -1.0d0

        ! The error is measured on the conserved quantities, not on T and p.
        !
        ! A divided-difference estimator assumes the quantity it differences is
        ! smooth in time. Temperature is not: a node crossing the freezing point
        ! has a kink in T(t), because the latent heat is absorbed at a nearly
        ! constant temperature. Its divided differences then scale like h, not
        ! h^(k+1), so the estimate does not fall when the step is cut. Measured
        ! from one accepted state: dt 111.8 s -> E 1.427, dt 84.2 s -> E 4.427,
        ! dt 36.0 s -> E 1.456, where BDF2 requires E to drop by 30 over that
        ! range. The controller then rejects every step and the run dies by
        ! attrition without the nonlinear solve ever failing.
        !
        ! Enthalpy and water mass are smooth across the transition - that is the
        ! defining property of the enthalpy formulation - so their divided
        ! differences do converge, and they are also what the residual actually
        ! conserves and what the acceptance gate already measures.
        call self%compute_nodal_conserved(enthalpy, density)

        ! atol is left at zero here. An absolute floor in J/m3 or kg/m3 has no
        ! natural value, and the weight already carries the field's own span as
        ! a floor, which is what the temperature atol was standing in for.
        if (self%is_active_thermal() .and. allocated(enthalpy)) then
            call physics_lte(enthalpy, self%lte_state_prev_thermal, &
                             self%lte_ydot_prev_thermal, self%lte_d2_prev_thermal, &
                             0.0d0, bdf_order, e_thermal)
        end if
        if (self%is_active_hydraulic() .and. allocated(density)) then
            call physics_lte(density, self%lte_state_prev_hydraulic, &
                             self%lte_ydot_prev_hydraulic, self%lte_d2_prev_hydraulic, &
                             0.0d0, bdf_order, e_hydraulic)
        end if

        error_rel = max(e_thermal, e_hydraulic)

        if (LTE_BLOCK_REPORT) then
            write (*, '(A,ES11.3,A,ES11.3,A,ES11.3)') '   [LTE] blocks: enthalpy=', e_thermal, &
                ' water_mass=', e_hydraulic, ' dt[s]=', dt_n
        end if

    contains

        subroutine physics_lte(y, state_prev, ydot_prev, d2_prev, atol, order, e_norm)
            implicit none
            !> Conserved quantity at the accepted iterate
            real(real64), intent(in) :: y(:)
            real(real64), allocatable, intent(in) :: state_prev(:)
            real(real64), allocatable, intent(in) :: ydot_prev(:)
            real(real64), allocatable, intent(in) :: d2_prev(:)
            real(real64), intent(in) :: atol
            integer(int32), intent(in) :: order
            real(real64), intent(inout) :: e_norm

            real(real64), allocatable :: ydot(:), d2(:), defect(:), weight(:)
            real(real64) :: span_local(2), span_global(2), y_span
            real(real64) :: accum_local(2), accum_global(2)
#ifdef _MPI
            integer(int32) :: ierr
#endif

            e_norm = -1.0d0
            if (.not. allocated(state_prev) .or. .not. allocated(ydot_prev)) return
            if (size(state_prev) /= size(y) .or. size(ydot_prev) /= size(y)) return
            if (size(y) == 0) return
            if (order >= 2) then
                if (.not. allocated(d2_prev)) return
                if (size(d2_prev) /= size(y)) return
            end if

            ! Reference scale flooring the relative part of the weight. A purely
            ! relative tolerance is meaningless for a variable whose zero is not
            ! a natural zero: Celsius temperature has its origin AT the
            ! phase-change point, so rtol*|T| vanishes exactly at the freezing
            ! front and the weight collapses to atol there while the rest of the
            ! domain is weighted by rtol*|T|. That makes the norm one to two
            ! orders stricter at the front than anywhere else, and a handful of
            ! front nodes then dominate the estimate and make it jump
            ! non-smoothly with dt. Flooring at the field's own span restores a
            ! uniform criterion. Gauge pressure has the same arbitrary origin.
            ! Reduce with MPI_MAX on (max, -min) so the span is the global one.
            span_local(1) = maxval(state_prev)
            span_local(2) = -minval(state_prev)
#ifdef _MPI
            call MPI_Allreduce(span_local, span_global, 2, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
            span_global = span_local
#endif
            y_span = max(0.0d0, span_global(1) + span_global(2))

            allocate (ydot(size(y)), d2(size(y)), defect(size(y)), weight(size(y)))
            ydot(:) = (y(:) - state_prev(:)) / dt_n

            ! Newton divided differences on the accepted times:
            !   y[t_n,t_{n-1}]                 = ydot
            !   y[t_n,t_{n-1},t_{n-2}]         = (ydot - ydot_prev)/(h_n + h_{n-1})
            !   y[t_n,...,t_{n-3}]             = (d2 - d2_prev)/(h_n + span_prev)
            ! and LTE_k ~ h_n**(k+1) * y[t_n,...,t_{n-k-1}].
            d2(:) = (ydot(:) - ydot_prev(:)) / (dt_n + self%lte_prev_dt)
            if (order >= 2) then
                defect(:) = dt_n**3 * (d2(:) - d2_prev(:)) / (dt_n + self%lte_prev_span)
            else
                defect(:) = dt_n**2 * d2(:)
            end if

            weight(:) = atol + self%lte_rtol * max(abs(y(:)), abs(state_prev(:)), y_span)

            ! One weighted RMS over the whole distributed problem: a per-rank
            ! norm would hand each rank a different dt and desynchronize them.
            ! Weighted by nodal volume, not node count, so the estimate does not
            ! grow simply because the surface mesh is graded finer.
            if (allocated(self%nodal_volume)) then
                if (size(self%nodal_volume) == size(y)) then
                    accum_local(1) = sum(self%nodal_volume(:) * &
                                         (defect(:) / max(weight(:), tiny(1.0d0)))**2)
                    accum_local(2) = sum(self%nodal_volume(:))
                else
                    accum_local(1) = sum((defect(:) / max(weight(:), tiny(1.0d0)))**2)
                    accum_local(2) = real(size(y), real64)
                end if
            else
                accum_local(1) = sum((defect(:) / max(weight(:), tiny(1.0d0)))**2)
                accum_local(2) = real(size(y), real64)
            end if
#ifdef _MPI
            call MPI_Allreduce(accum_local, accum_global, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
            accum_global = accum_local
#endif
            if (accum_global(2) <= 0.0d0) return
            e_norm = sqrt(accum_global(1) / accum_global(2))
        end subroutine physics_lte

    end function compute_lte_error_ftcms

    module subroutine initialize_lte_history_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        integer(int32) :: target_order

        self%lte_has_state = .false.
        self%lte_has_derivative = .false.
        self%lte_has_second_difference = .false.
        self%lte_prev_dt = 0.0d0
        self%lte_prev_span = 0.0d0

        ! Refuse an unsupported target order up front rather than mid-run: the
        ! integrator ramps current_bdf_order up to the target, so a target above
        ! what the estimator implements would silently start under-estimating
        ! the error once the ramp reaches it.
        if (self%lte_error_control_active) then
            call self%control%get_target_bdf_order(target_order)
            if (target_order > LTE_MAX_SUPPORTED_ORDER) then
                write (*, '(A,I0,A,I0,A)') ' Error: error-controlled stepping requests BDF order ', &
                    target_order, ' but the local-error estimator implements up to ', &
                    LTE_MAX_SUPPORTED_ORDER, '.'
                error stop 'initialize_lte_history: unsupported bdf_order for error control.'
            end if
        end if

        ! The estimator differences the conserved quantities (see
        ! compute_lte_error), so the history it starts from must hold those.
        block
            real(real64), allocatable :: enthalpy(:), density(:)

            call self%compute_nodal_conserved(enthalpy, density)
            if (self%is_active_thermal()) call copy_current(enthalpy, self%lte_state_prev_thermal)
            if (self%is_active_hydraulic()) call copy_current(density, self%lte_state_prev_hydraulic)
        end block
        self%lte_has_state = allocated(self%lte_state_prev_thermal) .or. &
                             allocated(self%lte_state_prev_hydraulic)

    contains
        subroutine copy_current(current, state_copy)
            implicit none
            real(real64), allocatable, intent(in) :: current(:)
            real(real64), allocatable, intent(inout) :: state_copy(:)

            if (.not. allocated(current)) return
            if (allocated(state_copy)) deallocate (state_copy)
            allocate (state_copy, source=current)
        end subroutine copy_current
    end subroutine initialize_lte_history_ftcms

    module subroutine commit_lte_history_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        real(real64) :: dt_n
        logical :: committed, second_difference_ready

        call self%control%get_dt(dt_n)
        if (dt_n <= 0.0d0) return
        committed = .false.
        ! A second difference can only be formed once ydot_prev already holds
        ! the previous step's slope, i.e. from the second accepted step on.
        second_difference_ready = self%lte_has_derivative .and. dt_n + self%lte_prev_dt > 0.0d0

        block
            real(real64), allocatable :: enthalpy(:), density(:)

            call self%compute_nodal_conserved(enthalpy, density)
            if (self%is_active_thermal() .and. allocated(enthalpy)) then
                call commit_physics(enthalpy, self%lte_state_prev_thermal, &
                                    self%lte_ydot_prev_thermal, self%lte_d2_prev_thermal, committed)
            end if
            if (self%is_active_hydraulic() .and. allocated(density)) then
                call commit_physics(density, self%lte_state_prev_hydraulic, &
                                    self%lte_ydot_prev_hydraulic, self%lte_d2_prev_hydraulic, committed)
            end if
        end block

        if (committed) then
            ! lte_prev_span is the time t_n - t_{n-2} that the second difference
            ! just stored was formed over; after the shift it plays the role of
            ! t_{n-1} - t_{n-3} for the next step's third difference.
            if (second_difference_ready) then
                self%lte_prev_span = dt_n + self%lte_prev_dt
                self%lte_has_second_difference = .true.
            end if
            self%lte_prev_dt = dt_n
            self%lte_has_state = .true.
            self%lte_has_derivative = .true.
        end if

    contains
        subroutine commit_physics(current, state_prev, ydot_prev, d2_prev, did_commit)
            implicit none
            !> Conserved quantity at the newly accepted step
            real(real64), intent(in) :: current(:)
            real(real64), allocatable, intent(inout) :: state_prev(:)
            real(real64), allocatable, intent(inout) :: ydot_prev(:)
            real(real64), allocatable, intent(inout) :: d2_prev(:)
            logical, intent(inout) :: did_commit
            real(real64), allocatable :: ydot_new(:)

            if (size(current) == 0) return
            if (.not. allocated(state_prev)) then
                allocate (state_prev, source=current)
                return
            end if
            if (size(state_prev) /= size(current)) then
                deallocate (state_prev)
                allocate (state_prev, source=current)
                return
            end if
            if (allocated(ydot_prev)) then
                if (size(ydot_prev) /= size(current)) deallocate (ydot_prev)
            end if

            allocate (ydot_new(size(current)))
            ydot_new(:) = (current(:) - state_prev(:)) / dt_n

            ! Store y[t_n,t_{n-1},t_{n-2}] before ydot_prev is overwritten.
            if (allocated(ydot_prev) .and. second_difference_ready) then
                if (allocated(d2_prev)) then
                    if (size(d2_prev) /= size(current)) deallocate (d2_prev)
                end if
                if (.not. allocated(d2_prev)) allocate (d2_prev(size(current)))
                d2_prev(:) = (ydot_new(:) - ydot_prev(:)) / (dt_n + self%lte_prev_dt)
            end if

            if (.not. allocated(ydot_prev)) allocate (ydot_prev(size(current)))
            ydot_prev(:) = ydot_new(:)
            state_prev(:) = current(:)
            did_commit = .true.
        end subroutine commit_physics
    end subroutine commit_lte_history_ftcms

    !> See the interface. Cost: O(N_dof). Combines the energy and water residual
    !> blocks by Euclidean norm; used only as a monotone merit for the line search.
    module function nonlinear_residual_norm_ftcms(self) result(rnorm)
        implicit none
        class(type_ftcms), intent(inout) :: self
        real(real64) :: rnorm

        real(real64), allocatable :: res(:)
        real(real64) :: sumsq

        sumsq = 0.0d0
        if (self%is_active_thermal()) then
            call self%get_variable_residual(PHYSICS_TYPES%THERMAL, res)
            if (allocated(res)) then
                if (size(res) > 0) sumsq = sumsq + vector_norm2(res)**2
            end if
        end if
        if (self%is_active_hydraulic()) then
            call self%get_variable_residual(PHYSICS_TYPES%HYDRAULIC, res)
            if (allocated(res)) then
                if (size(res) > 0) sumsq = sumsq + vector_norm2(res)**2
            end if
        end if
        rnorm = sqrt(sumsq)
    end function nonlinear_residual_norm_ftcms

    module subroutine reset_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        call self%control%reset_iteration()

    end subroutine reset_ftcms

    module function is_active_thermal_ftcms(self) result(is_active)
        implicit none
        class(type_ftcms), intent(in) :: self
        logical :: is_active

        is_active = self%control%is_physics_active(PHYSICS_TYPES%THERMAL)

    end function is_active_thermal_ftcms

    module function is_active_hydraulic_ftcms(self) result(is_active)
        implicit none
        class(type_ftcms), intent(in) :: self
        logical :: is_active

        is_active = self%control%is_physics_active(PHYSICS_TYPES%HYDRAULIC)

    end function is_active_hydraulic_ftcms

    module subroutine destroy_type_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        integer(int32) :: log_io_unit
        ! --- Stop and Record Profiler ---
        call self%control%profiler_stop(PROFILER_TYPES%TOTAL)
        call self%control%profiler_record(TIME_RECORDS%END)

        call self%output%output_system_log()
        call self%output%get_log_io_unit(log_io_unit)
        call self%control%display_profiler(log_io_unit)

        if (self%solver_history_unit /= -1) then
            close (self%solver_history_unit)
            self%solver_history_unit = -1
        end if

    end subroutine destroy_type_ftcms
    !> Give the solver the dof numbering the DM produced, so it can place this
    !> rank's right-hand side into the distributed vectors. The numbering is
    !> the matrix's own, which is what keeps the two consistent.
    subroutine hand_dof_map_to_solver(solver, jacobian, physics_type)
        implicit none
        class(abst_solver), allocatable, intent(inout) :: solver
        type(type_jacobian_matrix), intent(in) :: jacobian
        type(type_constant_id), intent(in), optional :: physics_type

        type(type_dof_layout) :: layout
        integer(int32), allocatable :: global_dof(:)
        integer(int32) :: node, block_offset, index

        if (.not. allocated(solver)) return

        call jacobian%get_dof_layout(layout, physics_type)
        if (.not. allocated(layout%node_dof_base)) return
        if (layout%num_dofs_per_node <= 0) return

        allocate (global_dof(size(layout%node_dof_base) * layout%num_dofs_per_node))
        index = 0
        do node = 1, size(layout%node_dof_base)
            do block_offset = 1, layout%num_dofs_per_node
                index = index + 1
                global_dof(index) = layout%node_dof_base(node) + block_offset - 1
            end do
        end do

        call solver%set_dof_map(global_dof)
        deallocate (global_dof)
    end subroutine hand_dof_map_to_solver

end submodule ftcms_base
