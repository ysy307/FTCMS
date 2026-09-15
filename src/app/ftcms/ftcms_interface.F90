module app_ftcms
    use, intrinsic :: iso_fortran_env
    use :: omp_lib
    use :: core_parallel_mpi

    use :: stdlib_optval, only:optval
    use :: stdlib_logger
    use :: module_core
    use :: module_input, only:type_input, input_translator
    use :: module_output, only:type_output_manager

    use :: module_control, only:type_control
    use :: control_acceleration_aa1_coupled, only:compute_aa1_coefficient, &
        compute_coupled_aa1_coefficient, coupled_weighted_norm
    use :: module_domain
    ! use :: module_boundary, only:
    use :: module_initial, only:type_ic_manager
    use :: module_system, only:type_jacobian_matrix, type_residual_vector
    use :: module_constitutive, only:g => gravity_acceleration, rho_std => reference_water_density, &
        Lf0 => latent_heat_fusion_water_0C, Tf0 => water_freezing_point_at_standard_atmospheric_pressure, &
        TtoK => celsius_to_kelvin
    use :: models_phase_change_chemical_potential, only:calc_T_high_celsius
    use :: module_linalg

    use :: module_governing
    use :: governing_atmosphere, only: type_da_config, type_assimilation_controller
    use :: module_solver
    implicit none

    !> Physical-validity walls shared by the bounded solution update
    !> (reflect_variables) and the acceptance guard of the nonlinear loop.
    !> An iterate pinned at a wall is outside the model's validity and must
    !> never be accepted as a converged step.
    real(real64), parameter :: WALL_TEMP_MIN_C = -80.0d0
    real(real64), parameter :: WALL_TEMP_MAX_C = 80.0d0
    real(real64), parameter :: WALL_PRESS_MIN_PA = -1.0d7
    real(real64), parameter :: WALL_PRESS_MAX_PA = 1.0d7

    !> Overwrite the Gauss-point dQi/dT with the node-interpolated value.
    !>
    !> update_water_phases computes dQi/dT consistently at the Gauss point from
    !> that point's own (T, p); replacing it by an interpolation of the nodal
    !> values makes the apparent heat capacity in K_TT a different number from
    !> the derivative of the enthalpy the residual integrates at that point.
    !> Default is direct Gauss-point evaluation (.false.); the lerp path is
    !> kept only as a regression-comparison switch. Declared here (rather than
    !> in ftcms_assemble.F90, where the production path consumes it) because
    !> the FD-Jacobian audit in ftcms_fd_jacobian.F90 must gate the same
    !> switch, and a parameter private to one submodule is invisible to a
    !> sibling submodule.
    logical, parameter :: LERP_NODAL_DQI_DT = .false.

    !> Replace all four coupled element blocks (K_TT, K_TH, K_HT, K_HH) with
    !> their element finite-difference tangent, instead of only K_TH.
    !>
    !> The residual is now C^1 (smooth pore-limited cryosuction), so the true
    !> element tangent is well-defined everywhere a Newton step is taken. The
    !> remaining hand-written blocks are not that tangent: K_TT, K_HT and K_HH
    !> are Modified-Picard blocks built from transport coefficients (impedance,
    !> thermal conductivity, geometric-mean permeability) and a chord heat
    !> capacity that are lagged at the current iterate rather than
    !> differentiated, so the assembled operator is not dR/du at a developed
    !> freezing front. Measured on such a step, the resulting update was an
    !> ascent direction for the line-search merit (ratio 5194 at alpha=1,
    !> still 8.95 at alpha=0.0078), and the nonlinear solve stalled at dt_min.
    !> Declared here for the same reason as LERP_NODAL_DQI_DT: both
    !> ftcms_assemble.F90 (production wiring) and ftcms_fd_jacobian.F90 (the
    !> new routine) must see the same switch.
    !>
    !> When .true., an element where both physics are active gets all four
    !> blocks from assemble_tangent_fd_ftcms; COUPLING_TH_FINITE_DIFFERENCE is
    !> then moot for such elements (see ftcms_assemble.F90). When .false., the
    !> previous TH-only finite-difference path remains selectable through
    !> COUPLING_TH_FINITE_DIFFERENCE.
    logical, parameter :: FULL_FD_TANGENT = .true.

    !> Bound the coupled Newton update by a variable-wise trust region.
    !>
    !> Measured without it: single updates of 23 K and 13 MPa, and nodal water
    !> contents moving by 0.32 (62 percent of the porosity) in one iteration -
    !> far outside the range where the local linearisation the step came from
    !> is meaningful, even when the direction itself is sound.
    !>
    !> The bound also repairs a structural defect in the unbounded path. The
    !> temperature and pressure blocks each computed their OWN validity-wall
    !> factor, so whenever one block's wall bound bound and the other's did
    !> not, T and p were scaled differently: that does not shorten the Newton
    !> step, it ROTATES it, and the line search is then no longer searching
    !> along the direction it was handed. With the bound enabled a single
    !> factor scales both blocks.
    !> Print one line per nonlinear iterate with the acceptance gates and the
    !> increment magnitudes.
    logical, parameter :: NONLINEAR_VERBOSE = .true.
    logical, parameter :: STEP_BOUND_ENABLED = .true.
    !> Largest temperature change admitted in one nonlinear update [K].
    real(real64), parameter :: STEP_BOUND_DT = 5.0d-1
    !> Largest pore-pressure change admitted in one nonlinear update [Pa].
    !> 50 kPa is about 0.04 K of Clapeyron freezing-point shift, so the two
    !> bounds constrain the phase state by comparable amounts.
    real(real64), parameter :: STEP_BOUND_DP = 5.0d4

    !>
    !> Static per-node table of the distinct materials among a node's
    !> adjacent elements, built once after the domain is initialized (the
    !> mesh and its material assignment are static for the lifetime of a
    !> run). Consumed by update_nodal_phases to evaluate the constitutive
    !> update once per distinct material per node, instead of once per
    !> adjacent element.
    !>
    !> CSR-like layout keyed by node id: node i's entries occupy the range
    !> ptr(i) : ptr(i+1)-1 in material_id / repr_element / measure_sum.
    !> Nodes with no adjacent elements have an empty row (ptr(i) == ptr(i+1)).
    !>
    !> repr_element(k) is the highest-index adjacent element carrying
    !> material_id(k): node->element adjacency is visited in ascending
    !> element-id order, so repeatedly overwriting the representative for a
    !> recurring material naturally leaves the highest-index one recorded.
    !> When a node touches a single material (the common case), this row is
    !> exactly the node's highest-index adjacent element overall, matching
    !> the "last element wins" semantics of the original per-element loop.
    !>
    type :: type_node_material_table
        !> Number of nodes the table was built for (0 if not yet built).
        integer(int32) :: num_nodes = 0
        !> Row pointers (size num_nodes + 1).
        integer(int32), allocatable :: ptr(:)
        !> Distinct material id for each (node, material) entry.
        integer(int32), allocatable :: material_id(:)
        !> Representative (highest-index) adjacent element for each entry.
        integer(int32), allocatable :: repr_element(:)
        !> Summed measure (volume/area) of adjacent elements sharing this material.
        real(real64), allocatable :: measure_sum(:)
    contains
        ! ---- Lifecycle ----
        procedure, pass(self), public :: initialize => initialize_node_material_table
        procedure, pass(self), public :: destroy => destroy_node_material_table
    end type type_node_material_table

    !>
    !> Static per-element Gauss-point geometry cache for the L2 (lumped
    !> mass) nodal gradient projection, built once after the domain is
    !> initialized (the mesh is static for the lifetime of a run). Caches
    !> exactly the quantities the projection re-derived per call: the
    !> shape-function nodal weights psi_k(x_p) * w_p * det(J_p), the global
    !> shape-function derivatives dpsi_k/dx(x_p), and the lumped nodal
    !> volume, so calc_gradient no longer calls calc_shape_function.
    !>
    !> Flattened CSR-like layout: element e's (gauss, node) entries occupy
    !> entry_ptr(e) : entry_ptr(e+1)-1 in shape_weight / dpsi_dx, ordered
    !> Gauss-point-major then local-node (the projection's loop order).
    !>
    type :: type_gradient_geometry_cache
        !> Number of elements the cache was built for (0 if not yet built).
        integer(int32) :: num_elements = 0
        !> Number of global nodes (size of nodal_vol).
        integer(int32) :: num_nodes = 0
        !> Spatial dimension of the computation (rows of dpsi_dx).
        integer(int32) :: dim = 0
        !> Entry offsets per element (size num_elements + 1).
        integer(int32), allocatable :: entry_ptr(:)
        !> Number of Gauss points per element (size num_elements).
        integer(int32), allocatable :: num_gauss(:)
        !> Number of local nodes per element (size num_elements).
        integer(int32), allocatable :: num_nodes_elem(:)
        !> psi_k(x_p) * w_p * det(J_p) per (element, gauss, node) entry.
        real(real64), allocatable :: shape_weight(:)
        !> Global shape-function derivatives dpsi_k/dx_d(x_p), shape (dim, entries).
        real(real64), allocatable :: dpsi_dx(:, :)
        !> Lumped nodal volume sum over adjacent Gauss points of
        !> psi_k w_p det(J_p) (size num_nodes).
        real(real64), allocatable :: nodal_vol(:)
    contains
        ! ---- Lifecycle ----
        procedure, pass(self), public :: initialize => initialize_gradient_geometry_cache
        procedure, pass(self), public :: destroy => destroy_gradient_geometry_cache
    end type type_gradient_geometry_cache

    type :: type_ftcms
        type(type_domain) :: domain

        !> Static per-node distinct-material table, built once from the
        !> domain's node->element adjacency (see type_node_material_table).
        !> Consumed by update_nodal_phases.
        type(type_node_material_table) :: node_material_table

        !> Static Gauss-point geometry cache for the nodal gradient
        !> projection (see type_gradient_geometry_cache). Consumed by
        !> calc_gradient.
        type(type_gradient_geometry_cache) :: gradient_cache

        !> Quadrature depth selected for each element by the latest complete
        !> matrix assembly.  Residual-only and finite-difference assemblies
        !> reuse it so their discrete residual has one fixed integration rule.
        !> True when the depth-d/depth-(d+1) comparison still exceeded the
        !> tolerance at the maximum supported depth.

        type(type_variable) :: porosity
        type(type_variable) :: temperature
        type(type_variable) :: pressure

        !> Nodal d(psi_cryo)/dT [Pa/K]. The hydraulic increment is solved in the
        !> total-potential variable du_g = du_p - c du_T and recovered with this.
        real(real64), allocatable :: cryo_slope(:)

        type(type_variable) :: Qw
        type(type_variable) :: Qi
        type(type_variable) :: Qa
        type(type_variable) :: Qv

        type(type_jacobian_matrix) :: K
        type(type_residual_vector) :: F
        !> 1.0 at every Dirichlet degree of freedom, in the same layout as F.
        type(type_residual_vector) :: dirichlet_marker
        type(type_residual_vector) :: du

        type(type_thermal) :: thermal
        type(type_hydraulic) :: hydraulic

        type(type_bc_manager) :: bc(PHYSICS_TYPES%NUM_ID)

        class(abst_solver), allocatable :: solver
        class(abst_solver), allocatable :: solver_thermal

        integer(int32) :: current_physics_id = 0
        ! Number of outer phase-equilibrium iterations in the latest solve.
        integer(int32) :: last_phase_iterations = 1
        ! Iterations in the final inner nonlinear solve.
        integer(int32) :: last_inner_iterations = 0
        ! Sum of inner nonlinear iterations over all outer phase iterations.
        integer(int32) :: last_nonlinear_work = 0
        integer(int32) :: last_solve_status = 0
        ! Phase-equilibrium metrics feeding the Anderson secant restart test.
        logical :: last_phase_metrics_available = .false.
        real(real64) :: last_phase_increment_max = -1.0d0
        ! Time increment of the latest accepted step. Used only to extrapolate
        ! the initial ice-state guess; it does not enter the governing equations.
        real(real64) :: last_accepted_dt = 0.0d0
        integer(int32) :: thermal_start_dof = 0
        integer(int32) :: hydraulic_start_dof = 0
        logical :: hydraulic_has_dirichlet_bc = .false.
        logical :: thermal_has_dirichlet_bc = .false.

        ! Reference mean pressure captured from the initial condition.
        ! Used to pin the null-mode of all-Neumann hydraulic systems without
        ! distorting the absolute pressure level (WRF relies on absolute P).

        ! DOF column scaling factors for variable non-dimensionalization
        real(real64), allocatable :: col_scale(:)
        real(real64), allocatable :: col_scale_inv(:)

        ! Accepted-step state for transactional BDF1/BDF2 LTE estimation.
        !
        ! The estimator is the Newton divided-difference form
        !   LTE_k ~ h_n**(k+1) * y[t_n, t_{n-1}, ..., t_{n-k-1}],
        ! so order k needs the k-th difference of the previous step as history.
        ! Storing that difference costs one array instead of one more full
        ! solution level per order.
        real(real64), allocatable :: lte_state_prev_thermal(:)
        real(real64), allocatable :: lte_state_prev_hydraulic(:)
        real(real64), allocatable :: lte_ydot_prev_thermal(:)
        real(real64), allocatable :: lte_ydot_prev_hydraulic(:)
        !> Second divided difference y[t_{n-1},t_{n-2},t_{n-3}] of the previous
        !> accepted step, needed only by the BDF2 estimator.
        real(real64), allocatable :: lte_d2_prev_thermal(:)
        real(real64), allocatable :: lte_d2_prev_hydraulic(:)
        !> Total mesh measure, summed once and reused by the residual floor.
        real(real64) :: domain_measure_total = -1.0d0
        !> Nodal control volumes V_i = sum_e int psi_i dOmega [m3], one entry
        !> per domain node (row-sum-lumped mass matrix with unit coefficient),
        !> built once (the mesh is static) and handed to the convergence
        !> control's exact per-node local conservation gate. See
        !> solve_time_step_initial_setup_ftcms.
        real(real64), allocatable :: nodal_volume(:)
        real(real64) :: lte_prev_dt = 0.0d0
        !> t_{n-1} - t_{n-3} spanned by lte_d2_prev.
        real(real64) :: lte_prev_span = 0.0d0
        real(real64) :: lte_rtol = 1.0d-2
        real(real64) :: lte_atol_thermal = 1.0d-3
        real(real64) :: lte_atol_hydraulic = 1.0d1
        logical :: lte_error_control_active = .false.
        logical :: lte_has_state = .false.
        logical :: lte_has_derivative = .false.
        logical :: lte_has_second_difference = .false.

        ! Anderson(1) acceleration state of the conserved coupled Picard loop:
        ! the previous iterate and previous fixed-point increment of (T, p),
        ! used by reflect_variables to form the depth-1 Anderson mixing.
        ! Cleared at the start of every nonlinear loop.
        real(real64), allocatable :: aa_T_prev(:)
        real(real64), allocatable :: aa_P_prev(:)
        real(real64), allocatable :: aa_duT_prev(:)
        real(real64), allocatable :: aa_duP_prev(:)
        logical :: aa_has_prev = .false.
        ! Weighted norm of the previous fixed-point increment ||g_{k-1}||_W.
        ! Retained for nonlinear diagnostics; growth alone does not disable AA(1).
        real(real64) :: aa_gnorm_prev = -1.0d0
        ! Per-attempt active-set marker for Hansson's one-time freezing-onset
        ! temperature reset. It is not a primary unknown or accepted history.
        logical, allocatable :: phase_onset_reset(:)
        !> Line-search telemetry: accepted step scale, trials used and the count
        !> of iterates where no admissible step was found.
        real(real64) :: last_line_search_scale = 1.0d0
        integer(int32) :: last_line_search_trials = 0
        integer(int32) :: last_line_search_failures = 0
        !> Largest per-block increment of the last nonlinear iterate, in primary
        !> variable units [K] and [Pa]. With the tangent and the linear solve
        !> both verified, the size of du is the remaining observable that says
        !> whether the step is bounded. -1 before the first iterate.
        real(real64) :: last_du_thermal_max = -1.0d0
        real(real64) :: last_du_hydraulic_max = -1.0d0

        ! Unit of Output/solver_history.log: one record per time-step attempt,
        ! carrying attempted/accepted time and dt, termination status, inner
        ! work, the acceptance gates, the increments and the LTE ratio.
        ! -1 when the log is not open.
        integer(int32) :: solver_history_unit = -1

        type(type_control) :: control
        type(type_output_manager) :: output

        type(type_assimilation_controller) :: assimilation
        logical :: assimilation_enabled = .false.

    contains
        ! ---- Lifecycle ----
        ! initialize, destroy, reset, etc.
        procedure, public, pass(self) :: initialize => initialize_type_ftcms
        procedure, public, pass(self) :: destroy => destroy_type_ftcms

        ! ---- Mutator ----
        ! set_XXX, increment_XXX, update_XXX, etc.

        ! ---- Algorithm / Operation ----
        ! compute_XXX, check_XXX, solve_XXX, etc.

        ! ---- Inquiry ----
        ! is_XXX, has_XXX, should_XXX, etc.

        ! ---- Getter ----
        ! get_XXX, etc.

        ! ---- Meta / Utility ----
        ! display, to_string, etc.

        ! ---- Operator ----

        procedure, public, pass(self) :: shift => shift_ftcms

        procedure, public, pass(self) :: calc_gradient => calc_gradient_ftcms
        procedure, public, pass(self) :: calc_gradient_temperature => calc_gradient_temperature_ftcms
        procedure, public, pass(self) :: calc_gradient_pressure => calc_gradient_pressure_ftcms

        procedure, public, pass(self) :: calc_water_flux => calc_water_flux_ftcms
        procedure, public, pass(self) :: calc_vapor_flux => calc_vapor_flux_ftcms

        ! --- Boundary Condition Procedures ---
        procedure, public, pass(self) :: apply_bc => apply_bc_ftcms
        procedure, public, pass(self) :: prescribe_dirichlet => prescribe_dirichlet_ftcms
        procedure, private, pass(self) :: freeze_physics_dofs => freeze_physics_dofs_ftcms
        procedure, private, pass(self) :: zero_frozen_increment => zero_frozen_increment_ftcms
        procedure, private, pass(self) :: prescribe_essential_bc_generic
        procedure, private, pass(self) :: apply_natural_bc_generic
        procedure, private, pass(self) :: apply_essential_bc_generic

        procedure, public, pass(self) :: solve => solve_ftcms

        procedure, public, pass(self) :: set_state => set_state_ftcms
        procedure, private, pass(self) :: set_states_from_connectivity => set_states_from_connectivity_ftcms
        procedure, public, pass(self) :: update_physical_properties => update_physical_properties_ftcms
        procedure, private, pass(self) :: update_physical_properties_bulk => update_physical_properties_bulk_ftcms

        procedure, public, pass(self) :: reflect_variables => reflect_variables_ftcms
        procedure, private, pass(self) :: apply_phase_change_temperature_correction => &
            apply_phase_change_temperature_correction_ftcms
        procedure, private, pass(self) :: update_nodal_phases => update_nodal_phases_ftcms
        procedure, private, pass(self) :: project_nodal_ice => project_nodal_ice_ftcms
        procedure, private, pass(self) :: compute_nodal_conserved => compute_nodal_conserved_ftcms
        procedure, public, pass(self) :: compute_lte_error => compute_lte_error_ftcms
        procedure, private, pass(self) :: initialize_lte_history => initialize_lte_history_ftcms
        procedure, private, pass(self) :: commit_lte_history => commit_lte_history_ftcms
        procedure, public, pass(self) :: nonlinear_residual_norm => nonlinear_residual_norm_ftcms
        procedure, public, pass(self) :: update_variables => update_variables_ftcms
        procedure, public, pass(self) :: assemble_local => assemble_local_ftcms
        procedure, public, pass(self) :: assemble => assemble_ftcms
        procedure, public, pass(self) :: report_fd_jacobian => report_fd_jacobian_ftcms
        procedure, private, pass(self) :: assemble_coupling_fd => assemble_coupling_fd_ftcms
        procedure, private, pass(self) :: assemble_tangent_fd => assemble_tangent_fd_ftcms
        procedure, private, pass(self) :: assemble_initialize => assemble_initialize_ftcms
        procedure, private, pass(self) :: assemble_destroy => assemble_destroy_ftcms

        procedure, private, pass(self) :: get_variable_increment => get_variable_increment_ftcms
        procedure, private, pass(self) :: get_variable_residual => get_variable_residual_ftcms

        procedure, public, pass(self) :: reset => reset_ftcms

        !> Solve a single time step including the nonlinear iteration loop
        procedure, public, pass(self) :: solve_time_step => solve_time_step_ftcms
        procedure, private, pass(self) :: solve_time_step_staggered => solve_time_step_staggered_ftcms
        procedure, private, pass(self) :: solve_time_step_initial_setup => solve_time_step_initial_setup_ftcms
        procedure, private, pass(self) :: solve_time_step_setup => solve_time_step_setup_ftcms
        procedure, private, pass(self) :: solve_time_step_check_convergence => solve_time_step_check_convergence_ftcms
        procedure, private, pass(self) :: solve_time_step_check_convergence_conserved => &
            solve_time_step_check_convergence_conserved_ftcms

        procedure, public, pass(self) :: output_fields => output_fields_ftcms
        procedure, public, pass(self) :: output_history => output_history_ftcms

        procedure, public, pass(self) :: is_active_thermal => is_active_thermal_ftcms
        procedure, public, pass(self) :: is_active_hydraulic => is_active_hydraulic_ftcms

        procedure, public, pass(self) :: run => run_ftcms
        procedure, public, pass(self) :: run_assimilation => run_assimilation_ftcms

    end type type_ftcms

    interface
        !> Build the per-node distinct-material table from the (already
        !> initialized) domain's node->element adjacency. See
        !> type_node_material_table for the storage layout and semantics.
        !>
        !> Computational complexity: O(N_nd * d) time and O(N_nd + E) memory,
        !> where d is the maximum node degree and E is the total number of
        !> distinct (node, material) entries (E <= sum of node degrees).
        !> Failure behavior: a domain with zero nodes leaves the table empty
        !> (num_nodes == 0); no error is raised.
        module subroutine initialize_node_material_table(self, domain)
            implicit none
            class(type_node_material_table), intent(inout) :: self
            type(type_domain), intent(in) :: domain
        end subroutine initialize_node_material_table

        !> Release the table's allocatable storage and reset it to empty.
        module subroutine destroy_node_material_table(self)
            implicit none
            class(type_node_material_table), intent(inout) :: self
        end subroutine destroy_node_material_table

        !> Build the Gauss-point geometry cache for the nodal gradient
        !> projection from the (already initialized) domain. See
        !> type_gradient_geometry_cache for the storage layout. The cached
        !> values are computed by the same calc_shape_function evaluations
        !> the projection previously performed per call, in the same order,
        !> so a cache-based projection is bit-identical to the direct one.
        !>
        !> Computational complexity: O(E_gp) time and memory, where E_gp is
        !> the total number of (element, gauss, node) entries.
        !> Failure behavior: a domain with zero elements or nodes leaves the
        !> cache empty (num_elements == 0); no error is raised.
        module subroutine initialize_gradient_geometry_cache(self, domain)
            implicit none
            class(type_gradient_geometry_cache), intent(inout) :: self
            type(type_domain), intent(in) :: domain
        end subroutine initialize_gradient_geometry_cache

        !> Release the cache's allocatable storage and reset it to empty.
        module subroutine destroy_gradient_geometry_cache(self)
            implicit none
            class(type_gradient_geometry_cache), intent(inout) :: self
        end subroutine destroy_gradient_geometry_cache

        module subroutine initialize_type_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self

        end subroutine initialize_type_ftcms

        module subroutine prescribe_essential_bc_generic(self, physics_type, current_time, variable)
            implicit none
            class(type_ftcms), intent(inout), target :: self
            type(type_constant_id), intent(in) :: physics_type
            real(real64), intent(in) :: current_time
            type(type_variable), intent(inout) :: variable
        end subroutine prescribe_essential_bc_generic

        module subroutine apply_natural_bc_generic(self, physics_type, current_time, variable, dof_offset)
            implicit none
            class(type_ftcms), intent(inout), target :: self
            type(type_constant_id), intent(in) :: physics_type
            real(real64), intent(in) :: current_time
            type(type_variable), intent(in) :: variable
            integer(int32), intent(in) :: dof_offset

        end subroutine apply_natural_bc_generic

        module subroutine apply_essential_bc_generic(self, physics_type, current_time, variable, dof_offset)
            implicit none
            class(type_ftcms), intent(inout), target :: self
            type(type_constant_id), intent(in) :: physics_type
            real(real64), intent(in) :: current_time
            type(type_variable), intent(in) :: variable
            integer(int32), intent(in) :: dof_offset

        end subroutine apply_essential_bc_generic

        module subroutine apply_bc_ftcms(self, prescribed)
            implicit none
            class(type_ftcms), intent(inout) :: self
            logical, intent(in), optional :: prescribed

        end subroutine apply_bc_ftcms

        module subroutine prescribe_dirichlet_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self

        end subroutine prescribe_dirichlet_ftcms

        module subroutine freeze_physics_dofs_ftcms(self, physics_type)
            implicit none
            class(type_ftcms), intent(inout) :: self
            type(type_constant_id), intent(in) :: physics_type
        end subroutine freeze_physics_dofs_ftcms

        module subroutine zero_frozen_increment_ftcms(self, frozen_physics)
            import :: type_ftcms, type_constant_id
            implicit none
            class(type_ftcms), intent(inout) :: self
            type(type_constant_id), intent(in) :: frozen_physics
        end subroutine zero_frozen_increment_ftcms

        module subroutine solve_ftcms(self)
            import :: type_ftcms
            implicit none
            class(type_ftcms), intent(inout) :: self
        end subroutine solve_ftcms

        module subroutine set_state_ftcms(self, node_id, element_id, state, calc_physics, include_fluxes)
            implicit none
            class(type_ftcms), intent(inout) :: self
            integer(int32), intent(in) :: node_id
            integer(int32), intent(in) :: element_id
            type(type_state), intent(inout) :: state
            logical, intent(in), optional :: calc_physics
            !> When .false. (and calc_physics applies), skip the Darcy flux
            !> evaluation and leave state%water_flux / state%vapor_flux
            !> unset. Only valid for callers that never read those fields.
            !> Default .true. (full constitutive evaluation).
            logical, intent(in), optional :: include_fluxes
        end subroutine set_state_ftcms

        module subroutine set_states_from_connectivity_ftcms(self, connectivity, element_id, states, calc_physics)
            implicit none
            class(type_ftcms), intent(inout) :: self
            integer(int32), intent(in) :: connectivity(:)
            integer(int32), intent(in) :: element_id
            type(type_state), intent(inout) :: states(:)
            logical, intent(in), optional :: calc_physics
        end subroutine set_states_from_connectivity_ftcms

        module subroutine update_physical_properties_ftcms(self, material_id, state, include_fluxes)
            implicit none
            class(type_ftcms), intent(inout) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state
            !> When .false., perform the phase-change update only and skip
            !> the Darcy flux evaluation, leaving state%water_flux /
            !> state%vapor_flux unset. Default .true. (compute fluxes).
            logical, intent(in), optional :: include_fluxes
        end subroutine update_physical_properties_ftcms

        module subroutine update_physical_properties_bulk_ftcms(self, material_id, states)
            implicit none
            class(type_ftcms), intent(inout) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: states(:)
        end subroutine update_physical_properties_bulk_ftcms

        module subroutine shift_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self

        end subroutine shift_ftcms

        module subroutine update_variables_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self

        end subroutine update_variables_ftcms

        module subroutine reflect_variables_ftcms(self, step_scale)
            implicit none
            class(type_ftcms), intent(inout) :: self
            real(real64), intent(in), optional :: step_scale

        end subroutine reflect_variables_ftcms

        module subroutine update_nodal_phases_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self
        end subroutine update_nodal_phases_ftcms

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
        end subroutine project_nodal_ice_ftcms

        !> Evaluate the per-node conserved quantities (volumetric enthalpy density
        !> and pore-water effective density) at the current iterate, for the
        !> conserved-quantity nonlinear convergence norm. The optional dH_dT/
        !> drho_dp outputs are the matching nodal storage sensitivities (dH/dT,
        !> d rho_eq/dp), used by the primary-variable local error measure (see
        !> control_iteration_convergence:local_error_block); computed only when
        !> the corresponding actual argument is present.
        module subroutine compute_nodal_conserved_ftcms(self, enthalpy, density, dH_dT, drho_dp)
            implicit none
            class(type_ftcms), intent(inout) :: self
            real(real64), allocatable, intent(inout) :: enthalpy(:)
            real(real64), allocatable, intent(inout) :: density(:)
            real(real64), allocatable, intent(inout), optional :: dH_dT(:)
            real(real64), allocatable, intent(inout), optional :: drho_dp(:)
        end subroutine compute_nodal_conserved_ftcms

        !> Normalized local-truncation-error estimate of the tentative BDF1 step.
        !>
        !> Uses the divided-difference (curvature) estimate of the implicit-Euler
        !> LTE: \( \text{LTE} \approx \lVert \dot y_n - \dot y_{n-1}\rVert\,
        !> \Delta t_n^2/(\Delta t_n+\Delta t_{n-1}) \), evaluated with a
        !> component-wise atol/rtol weighted RMS norm.
        !> Returns -1 on the first step (no previous derivative yet) so the caller
        !> skips error control. This function does not mutate accepted history.
        module function compute_lte_error_ftcms(self, order_used) result(error_rel)
            implicit none
            class(type_ftcms), intent(inout) :: self
            !> BDF order the estimate was actually measured at, which is not
            !> always the controller's configured order: the estimator falls
            !> back to 1 until the second difference exists. The step-size
            !> retry must be resized at this order, not the configured one.
            integer(int32), intent(inout), optional :: order_used
            real(real64) :: error_rel
        end function compute_lte_error_ftcms

        module subroutine initialize_lte_history_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self
        end subroutine initialize_lte_history_ftcms

        module subroutine commit_lte_history_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self
        end subroutine commit_lte_history_ftcms

        !> Euclidean norm of the assembled nonlinear conservation residual (energy
        !> and water blocks combined), used as the merit function for the
        !> backtracking line search that globalizes the modified-Picard step.
        module function nonlinear_residual_norm_ftcms(self) result(rnorm)
            implicit none
            class(type_ftcms), intent(inout) :: self
            real(real64) :: rnorm
        end function nonlinear_residual_norm_ftcms

        module subroutine calc_gradient_ftcms(self, values_vec, grad)
            implicit none
            class(type_ftcms), intent(inout) :: self
            real(real64), intent(in) :: values_vec(:)
            type(type_coordinate_array_dp), intent(inout) :: grad

        end subroutine calc_gradient_ftcms

        module subroutine calc_gradient_temperature_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self

        end subroutine calc_gradient_temperature_ftcms

        module subroutine calc_gradient_pressure_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self

        end subroutine calc_gradient_pressure_ftcms

        module subroutine calc_water_flux_ftcms(self, material_id, state, grad_T, grad_P, water_flux)
            implicit none
            class(type_ftcms), intent(inout) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            type(type_coordinate_dp), intent(in) :: grad_T, grad_P
            type(type_coordinate_dp), intent(inout) :: water_flux

        end subroutine calc_water_flux_ftcms

        module subroutine calc_vapor_flux_ftcms(self, material_id, state, grad_T, grad_P, water_flux)
            implicit none
            class(type_ftcms), intent(inout) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            type(type_coordinate_dp), intent(in) :: grad_T, grad_P
            type(type_coordinate_dp), intent(inout) :: water_flux

        end subroutine calc_vapor_flux_ftcms

        module subroutine assemble_local_ftcms(self, workspace, local_K_TT, local_K_TH, &
                                               local_K_HH, local_K_HT, local_F_T, local_F_H)
            implicit none
            class(type_ftcms), intent(inout) :: self
            type(type_assemble_workspace), intent(inout) :: workspace
            type(type_matrix_dense), intent(inout), optional :: local_K_TT, local_K_TH, local_K_HH, local_K_HT
            type(type_vector_dp), intent(inout), optional :: local_F_T, local_F_H

        end subroutine assemble_local_ftcms

        !> Compare every assembled block against the element-level
        !> finite-difference tangent of the same residual.
        !>
        !> Reports, per block and separately for freezing-front and control
        !> elements, the relative Frobenius error and the norm ratio
        !> \( \|J_{fd}\| / \|K\| \). Assumes both physics are active; returns
        !> immediately otherwise. Cost: \(O(N_{elem} \cdot 4 n_{node})\) element
        !> residual evaluations, so this is a diagnostic entry point, not a
        !> per-iteration path.
        module subroutine report_fd_jacobian_ftcms(self, tag, mean_rel_error, element_stride)
            implicit none
            class(type_ftcms), intent(inout) :: self
            !> Short label identifying the calling context in the log
            character(len=*), intent(in) :: tag
            !> Per-block mean relative error in report order
            !> (TT, TH, HT, HH, TH_fd); negative where nothing was sampled
            real(real64), intent(inout), optional :: mean_rel_error(:)
            !> Sample every n-th freezing-front element, >= 1 (default: all)
            integer(int32), intent(in), optional :: element_stride

        end subroutine report_fd_jacobian_ftcms

        !> Overwrite the thermal-pressure coupling block with the element
        !> finite-difference tangent \( \partial R_T/\partial p \).
        !>
        !> Assumes `local_F_T` holds the element thermal residual at the
        !> unperturbed state; it is restored on exit and the workspace is left
        !> on the unperturbed state. Cost: \(n_{node}\) extra thermal residual
        !> evaluations per element. Failure behavior: none - the difference
        !> quotient is always defined.
        module subroutine assemble_coupling_fd_ftcms(self, workspace, local_F_T, K_TH)
            implicit none
            class(type_ftcms), intent(inout) :: self
            type(type_assemble_workspace), intent(inout) :: workspace
            !> Element thermal residual; restored unchanged on exit
            type(type_vector_dp), intent(inout) :: local_F_T
            !> Overwritten with the finite-difference coupling block
            type(type_matrix_dense), intent(inout) :: K_TH

        end subroutine assemble_coupling_fd_ftcms

        !> Overwrite all four coupled element blocks (K_TT, K_TH, K_HT, K_HH)
        !> with the element finite-difference tangent of the same two
        !> residuals \( \partial R_T/\partial(T,p), \partial R_H/\partial(T,p) \).
        !>
        !> Assumes `local_F_T` and `local_F_H` hold the element thermal and
        !> hydraulic residuals at the unperturbed state; both are restored on
        !> exit and the workspace is left on the unperturbed state regardless
        !> of `ok`. Cost: \(2 n_{node}\) perturbations, each reading both
        !> residuals (one thermal and one hydraulic evaluation), against
        !> \(n_{node}\) perturbations of one residual for
        !> assemble_coupling_fd_ftcms. Unlike that routine, each perturbation
        !> here also needs a full nodal constitutive refresh (not only the
        !> Gauss-point one), because the hydraulic residual reads nodal
        !> derived state (effective_suction, transport coefficients) directly;
        !> see the submodule doc comment for why. Failure behavior: if any
        !> perturbed residual evaluation is non-finite, `ok` is returned
        !> .false. and K_TT/K_TH/K_HT/K_HH are left untouched, so the
        !> caller's previously assembled (analytic) blocks survive.
        module subroutine assemble_tangent_fd_ftcms(self, workspace, local_F_T, local_F_H, &
                                                    K_TT, K_TH, K_HT, K_HH, ok)
            implicit none
            class(type_ftcms), intent(inout) :: self
            type(type_assemble_workspace), intent(inout) :: workspace
            !> Element thermal residual; restored unchanged on exit
            type(type_vector_dp), intent(inout) :: local_F_T
            !> Element hydraulic residual; restored unchanged on exit
            type(type_vector_dp), intent(inout) :: local_F_H
            !> Overwritten with the finite-difference d R_T/d T block if ok
            type(type_matrix_dense), intent(inout) :: K_TT
            !> Overwritten with the finite-difference d R_T/d p block if ok
            type(type_matrix_dense), intent(inout) :: K_TH
            !> Overwritten with the finite-difference d R_H/d T block if ok
            type(type_matrix_dense), intent(inout) :: K_HT
            !> Overwritten with the finite-difference d R_H/d p block if ok
            type(type_matrix_dense), intent(inout) :: K_HH
            !> .true. if every perturbation produced a finite residual and the
            !> four blocks were overwritten; .false. if left untouched
            logical, intent(inout) :: ok

        end subroutine assemble_tangent_fd_ftcms

        module subroutine assemble_initialize_ftcms(self, element_id, workspace, local_K_TT, local_K_TH, &
                                                    local_K_HH, local_K_HT, local_F_T, local_F_H, &
                                                    coordinates, raw_coordinates, connectivity)
            implicit none

            class(type_ftcms), intent(inout) :: self
            integer(int32), intent(in) :: element_id
            type(type_assemble_workspace), intent(inout) :: workspace
            type(type_matrix_dense), intent(inout), optional :: local_K_TT, local_K_TH, local_K_HH, local_K_HT
            type(type_vector_dp), intent(inout), optional :: local_F_T, local_F_H
            real(real64), allocatable, intent(inout) :: coordinates(:, :)
            real(real64), allocatable, intent(inout) :: raw_coordinates(:, :)
            integer(int32), pointer, contiguous, intent(inout), optional :: connectivity(:)
        end subroutine assemble_initialize_ftcms

        module subroutine assemble_destroy_ftcms(self, workspace, local_K_TT, local_K_TH, &
                                                 local_K_HH, local_K_HT, local_F_T, local_F_H)
            implicit none
            class(type_ftcms), intent(inout) :: self
            type(type_assemble_workspace), intent(inout) :: workspace
            type(type_matrix_dense), intent(inout), optional :: local_K_TT, local_K_TH, local_K_HH, local_K_HT
            type(type_vector_dp), intent(inout), optional :: local_F_T, local_F_H

        end subroutine assemble_destroy_ftcms

        module subroutine get_variable_increment_ftcms(self, variable_id, variable)
            implicit none
            class(type_ftcms), intent(inout) :: self
            type(type_constant_id), intent(in) :: variable_id
            real(real64), intent(inout), allocatable :: variable(:)

        end subroutine get_variable_increment_ftcms

        module subroutine get_variable_residual_ftcms(self, variable_id, variable)
            implicit none
            class(type_ftcms), intent(inout) :: self
            type(type_constant_id), intent(in) :: variable_id
            real(real64), intent(inout), allocatable :: variable(:)

        end subroutine get_variable_residual_ftcms

        module subroutine reset_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self

        end subroutine reset_ftcms

        module subroutine assemble_ftcms(self, residual_only)
            implicit none
            class(type_ftcms), intent(inout) :: self
            !> Assemble only the residual vector when true.
            logical, intent(in), optional :: residual_only

        end subroutine assemble_ftcms

        module subroutine solve_time_step_initial_setup_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self

        end subroutine solve_time_step_initial_setup_ftcms

        module subroutine solve_time_step_setup_ftcms(self, prescribe_bc)
            implicit none
            class(type_ftcms), intent(inout) :: self
            logical, intent(inout) :: prescribe_bc

        end subroutine solve_time_step_setup_ftcms

        module subroutine solve_time_step_check_convergence_ftcms(self, target_physics)
            implicit none
            class(type_ftcms), intent(inout), target :: self
            type(type_constant_id), intent(in), optional :: target_physics

        end subroutine solve_time_step_check_convergence_ftcms

        !> Conserved-quantity convergence check (PDF 6.2.4), evaluated on the updated
        !> state: computes the nodal enthalpy/effective-density and per-block residual
        !> norms and delegates the coupled decision to the control manager.
        module subroutine solve_time_step_check_convergence_conserved_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self
        end subroutine solve_time_step_check_convergence_conserved_ftcms

        module subroutine solve_time_step_ftcms(self, is_step_converged)
            implicit none
            class(type_ftcms), intent(inout) :: self
            logical, intent(inout) :: is_step_converged

        end subroutine solve_time_step_ftcms

        module subroutine solve_time_step_staggered_ftcms(self, is_step_converged)
            implicit none
            class(type_ftcms), intent(inout) :: self
            logical, intent(inout) :: is_step_converged
        end subroutine solve_time_step_staggered_ftcms

        module subroutine output_fields_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self

        end subroutine output_fields_ftcms

        module subroutine output_history_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self

        end subroutine output_history_ftcms

        module subroutine run_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self

        end subroutine run_ftcms

        module subroutine run_assimilation_ftcms(self, current_time, current_doy)
            implicit none
            class(type_ftcms), intent(inout) :: self
            real(real64), intent(in) :: current_time
            real(real64), intent(in) :: current_doy
        end subroutine run_assimilation_ftcms

        module subroutine apply_phase_change_temperature_correction_ftcms(self, T_old, T_new)
            implicit none
            class(type_ftcms), intent(inout) :: self
            real(real64), intent(in) :: T_old(:)
            real(real64), intent(inout) :: T_new(:)
        end subroutine apply_phase_change_temperature_correction_ftcms

        module subroutine destroy_type_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self

        end subroutine destroy_type_ftcms

        module function is_active_thermal_ftcms(self) result(is_active)
            implicit none
            class(type_ftcms), intent(in) :: self
            logical :: is_active

        end function is_active_thermal_ftcms

        module function is_active_hydraulic_ftcms(self) result(is_active)
            implicit none
            class(type_ftcms), intent(in) :: self
            logical :: is_active

        end function is_active_hydraulic_ftcms

    end interface

end module app_ftcms
