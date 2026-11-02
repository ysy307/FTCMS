#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscsnes.h>

!> A monolithic freezing column: H1 temperature, RT0 liquid flux, P0 pore pressure.
!>
!> `set_hydraulic_space("primal")` puts the pore pressure in H1 instead and the
!> Darcy flux into the weak divergence of the water row
!> (`app_mixed_hydraulic_column_primal`); `"primal_nodal"` does the same with
!> the liquid pressure interpolated from the nodes and a geometric-mean element
!> conductivity, as the repository's CG code builds that term
!> (`app_mixed_hydraulic_column_nodal`). Everything else -- closure, heat row,
!> quadrature, scaling, Newton solver, step control -- is the same object, so
!> the runs differ in the hydraulic discretisation alone.
!>
!> The first coupled increment of the mixed hydraulic path. One backward-Euler
!> step solves heat and water together on one DMPlex, with the frozen
!> constitutive contract of `app_mixed_hydraulic_frozen` supplying every value
!> and every tangent. The hydraulic half is the same mixed pair the isothermal
!> solver verifies -- flux in \(H(\mathrm{div})\) so the normal trace is single
!> valued, pressure in \(L^2\) so the conservation equation holds cell by cell
!> -- and the temperature is added as a continuous nodal field on the same mesh.
!>
!> ## Which cells it runs on
!>
!> Triangles, and axis-aligned rectangles. The rectangle is there so that the
!> same problem can be put to the two element families -- P1/RT0/P0 against
!> Q1 and the tensor-product RT0 -- on meshes that share their vertices, which
!> is a comparison a triangulated mesh alone cannot make. It is a reference
!> geometry and not general quadrilateral support. An affine map is necessary
!> -- a warped cell has none, and the audit and the element loop both assume
!> one -- but axis alignment is asked for on top of it. A parallelogram is
!> affine and is still refused: what the accepted cells buy is a diagonal
!> Jacobian, so the closure's intrinsic x-z frame is each cell's own.
!>
!> ## What this is for
!>
!> A column whose top is cooled below freezing while its sides and both ends
!> are closed to water. If the closure's cryosuction is what the Mizoguchi
!> experiment measures, this is the smallest problem that can show it: with no
!> water crossing any boundary, any redistribution is internal, so an upward
!> flux toward the cold end is evidence about the model rather than about a
!> boundary condition. This solver does not reproduce that experiment, and a
!> run of it is not a validation of one.
!>
!> ## What it does not carry
!>
!> Advective heat transport by the liquid flux, a second time-step order, and
!> field output. Each is a separate increment; none is approximated here. The
!> step length does adapt, but only on the iteration count the solver reports:
!> there is no local error estimate behind it, so a run that meets its horizon
!> has shown that every step converged, not that the step was accurate.
!>
!> ## Ownership
!>
!> The frozen closure is module state, so only one column solver may be live on
!> a rank. It is not the isothermal solver's table and the two do not share
!> one; a run that wants both installs them in turn.
module app_mixed_hydraulic_column
    use, intrinsic :: iso_fortran_env, only: int32, real64, output_unit
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_funloc, c_null_funptr, c_null_ptr
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use :: petscdmplex
    use :: petscsnes
    use :: core_interop_petscds, only:petscds_set_residual, petscds_set_jacobian, &
                                      petscds_add_boundary_essential, petscds_add_boundary_natural
    use :: core_interop_petscfe, only:petscfe_evaluate_cell, petscfe_reference_coordinates, petscfe_tabulate, &
                                      petscfe_field_size
    use :: module_core, only:COMP_TYPES
    use :: app_mixed_hydraulic_physics, only:invert_mobility, MIXED_DIMENSION
    use :: app_mixed_hydraulic_balance, only:MAX_CELL_FACETS, MAX_QUADRATURE_POINTS, &
                                             reconstruct_cell_facet_flux, reconstruct_cell_pressure, &
                                             cell_is_owned, cell_volume, cell_quadrature_rule, &
                                             audit_facet_continuity
    use :: app_mixed_hydraulic_frozen, only:type_mixed_frozen_material, type_mixed_frozen_response, &
                                            install_frozen_laws, clear_frozen_laws, &
                                            eval_frozen_response, frozen_pressure_column
    use :: app_mixed_hydraulic_column_pointwise
    use :: app_mixed_hydraulic_column_primal, only:COLUMN_PRIMAL_FIELD_PRESSURE, &
                                                   f0_primal_conservation, f1_primal_conservation, &
                                                   g2_primal_pressure_pressure, g3_primal_pressure_pressure, &
                                                   g2_primal_pressure_temperature, g3_primal_pressure_temperature, &
                                                   compute_primal_darcy_flux
    use :: app_mixed_hydraulic_column_nodal, only:install_nodal_flux, clear_nodal_flux, compute_nodal_cell_flux
    implicit none
    private

    public :: type_mixed_column
    public :: type_mixed_column_report
    public :: type_mixed_column_step_control
    public :: type_mixed_column_progress
    public :: COLUMN_ERR_INVALID_STEP, COLUMN_ERR_NONFINITE, COLUMN_ERR_CLOSURE
    public :: COLUMN_ERR_STEP_FLOOR, COLUMN_ERR_INVALID_SPAN
    public :: COLUMN_ERR_NO_PROFILE_FILE, COLUMN_ERR_PROFILE_INCOMPLETE
    public :: COLUMN_ERR_NO_HISTORY_FILE
    public :: MAX_COLUMN_BOUNDARIES

    !> Returned when the requested step length is not a positive finite number.
    integer(int32), parameter :: COLUMN_ERR_INVALID_STEP = 2001
    !> Returned when a converged solve came back with a nonfinite entry.
    integer(int32), parameter :: COLUMN_ERR_NONFINITE = 2002
    !> Returned when the frozen closure could not answer at some iterate.
    integer(int32), parameter :: COLUMN_ERR_CLOSURE = 2003
    !> Returned when the step control reached its floor and still had no step.
    integer(int32), parameter :: COLUMN_ERR_STEP_FLOOR = 2004
    !> Returned when an advance was asked for a duration it cannot interpret.
    integer(int32), parameter :: COLUMN_ERR_INVALID_SPAN = 2005
    !> Returned when a profile file could not be opened for writing.
    integer(int32), parameter :: COLUMN_ERR_NO_PROFILE_FILE = 2007
    !> Returned when the gathered profile is not a one-to-one cover of the cells.
    integer(int32), parameter :: COLUMN_ERR_PROFILE_INCOMPLETE = 2008
    !> Returned when a step history file could not be opened for writing.
    integer(int32), parameter :: COLUMN_ERR_NO_HISTORY_FILE = 2009

    !> How many boundaries of each kind one run may register.
    integer(int32), parameter :: MAX_COLUMN_BOUNDARIES = 8

    !> How many numbers a profile row carries beside its time and cell number.
    !>
    !> In the order the header names them: the two centroid coordinates, the
    !> cell volume, the temperature, the two pressures, the three phase
    !> contents, the stored water, the relative conductivity, the two flux
    !> components, the cell's integrated water, and the water the cell's own
    !> flux has carried out through its upward- and downward-facing facets
    !> since the initial state.
    integer(int32), parameter :: COLUMN_PROFILE_VALUES = 16

    !> The quadrature the RT mass integral needs on an affine cell.
    !>
    !> This is the default the temperature element is created with. The element
    !> reads `-temperature_petscfe_default_quadrature_order` over it, so it is
    !> the rule that element ends up with, not this number, that the other
    !> fields and every auxiliary element are made to share.
    PetscInt, parameter :: QUADRATURE_ORDER = 2_PETSC_INT_KIND

    !> The water content one reference pressure change has to move before a
    !> cell's pressure column counts as carrying anything.
    real(real64), parameter :: COLUMN_STORAGE_REFERENCE = 1.0e-3_real64

    !> Wall seconds between the progress lines of a long advance. Longer than
    !> any run in the test suite, so those runs print none of them.
    real(real64), parameter :: ADVANCE_REPORT_SECONDS = 300.0_real64

    !> What one column audit measured, in physical units per unit thickness.
    type :: type_mixed_column_report
        !> Total stored water of the domain, m2 of water-equivalent volume.
        !>
        !> \(\sum_K\int_K W\), integrated with the rule the element loop uses,
        !> so this is the quantity the discretisation conserves and not a
        !> sampled approximation of it.
        real(real64) :: total_water = 0.0_real64
        !> Total ice of the domain, m2, on the same rule.
        real(real64) :: total_ice = 0.0_real64
        !> The same water total from one value per cell at its centroid, m2.
        !>
        !> Kept beside the integrated total because the difference between the
        !> two is the discretisation's own sampling error in \(W\), which a
        !> profile written cell by cell carries and a balance must not.
        real(real64) :: total_water_centroid = 0.0_real64
        !> Largest backward-Euler cell defect, m2.
        !>
        !> \(\max_K|\int_K(W^{n+1}-W^n)+\Delta t\sum_f Q_{Kf}-\Delta t\int_K s|\)
        !> over the last accepted step, on the assembly's own quadrature.
        real(real64) :: max_cell_defect = 0.0_real64
        !> The same defect with \(W\) sampled at the cell centroid instead, m2.
        real(real64) :: max_cell_defect_centroid = 0.0_real64
        !> Sum of the cell defects over all ranks, m2.
        real(real64) :: domain_defect = 0.0_real64
        !> Largest cell storage change of the last accepted step, m2.
        !>
        !> The scale a cell defect is small or large against: a defect is only
        !> evidence of conservation next to the water the step actually moved.
        real(real64) :: max_storage_change = 0.0_real64
        !> Largest \(|\Delta t\sum_f Q_{Kf}|\) of any cell, m2.
        !>
        !> The flux half of the same balance, reported separately because a
        !> defect that is only ever compared with a facet flux cannot tell a
        !> cell that stores water from one that fails to conserve it.
        real(real64) :: max_flux_term = 0.0_real64
        !> Largest \(|\int_K W-|K|W(x_K)|\) of any cell, m2.
        real(real64) :: max_quadrature_gap = 0.0_real64
        !> Largest disagreement between the two sides of an interior facet, m2/s.
        real(real64) :: max_facet_mismatch = 0.0_real64
        !> How many interior facets that maximum was taken over.
        integer(int32) :: interior_facets = 0
        !> Largest ratio of the relative conductivities a single cell's
        !> quadrature points carry, dimensionless and at least one.
        real(real64) :: max_conductivity_spread = 1.0_real64
        !> Largest spread of the liquid pressure inside one cell, Pa.
        real(real64) :: max_liquid_pressure_spread = 0.0_real64
        !> Net outward water flux through the domain boundary, m2/s.
        real(real64) :: boundary_rate = 0.0_real64
        !> Largest reconstructed facet flux, m2/s: the audit's own scale.
        real(real64) :: throughput = 0.0_real64
        !> Largest upward liquid flux, m/s, over the cells that carry ice.
        real(real64) :: upward_flux_frozen = 0.0_real64
        !> Largest downward liquid flux, m/s, over those same cells.
        real(real64) :: downward_flux_frozen = 0.0_real64
        !> The same two over every owned cell.
        real(real64) :: upward_flux = 0.0_real64
        real(real64) :: downward_flux = 0.0_real64
        !> Weakest pressure column of any owned cell, dimensionless.
        real(real64) :: weakest_column = 0.0_real64
        !> Coldest and warmest owned cell, C.
        real(real64) :: coldest = 0.0_real64
        real(real64) :: warmest = 0.0_real64
        !> Lowest liquid pressure of any owned cell, Pa.
        real(real64) :: lowest_liquid_pressure = 0.0_real64
        !> The step the defects belong to, s; zero before any step.
        real(real64) :: step_length = 0.0_real64
        !> How many owned cells carry ice.
        integer(int32) :: frozen_cells = 0
        !> How many owned cells this rank counted.
        integer(int32) :: owned_cells = 0
    end type type_mixed_column_report

    !> How an adaptive run chooses the length of its next step.
    !>
    !> The iteration-count controller of the project's own `Conditions.json`:
    !> a step the solver took easily is lengthened, one it laboured over is
    !> shortened, and one it refused is retried shorter from the state the
    !> refusal left untouched. Nothing here changes what a step converges to --
    !> the tolerances are the solver's -- so a run may be made cheaper or more
    !> careful without changing the answer a step is accepted at.
    type :: type_mixed_column_step_control
        !> The length the first step of a run is attempted at, s.
        real(real64) :: initial = 300.0_real64
        !> Below this the controller gives up rather than retry again, s.
        real(real64) :: minimum = 1.0_real64
        !> No step is lengthened past this, s.
        real(real64) :: maximum = 30000.0_real64
        !> Iteration counts the controller tries to keep a step between.
        integer(int32) :: iterations_low = 5
        integer(int32) :: iterations_high = 16
        !> Factors applied on an easy step, a laboured step and a refusal.
        real(real64) :: grow = 1.25_real64
        real(real64) :: shrink = 0.8_real64
        real(real64) :: retry = 0.5_real64
    end type type_mixed_column_step_control

    !> What one advance did: the ground truth a run reports instead of its plan.
    type :: type_mixed_column_progress
        !> Physical time the accepted steps covered, s.
        real(real64) :: reached = 0.0_real64
        !> The last step that was accepted and the one that would be tried next, s.
        real(real64) :: last_step = 0.0_real64
        real(real64) :: next_step = 0.0_real64
        !> The shortest step any accepted step of this advance was taken at, s.
        real(real64) :: shortest_step = 0.0_real64
        integer(int32) :: accepted_steps = 0
        integer(int32) :: rejected_steps = 0
        !> Newton iterations: the worst accepted step, and the sum over all of them.
        integer(int32) :: max_iterations = 0
        integer(int32) :: total_iterations = 0
    end type type_mixed_column_progress

    !> One registered prescribed-temperature boundary.
    type :: type_column_temperature_boundary
        integer(int32) :: boundary_id = 0
        real(real64) :: temperature = 0.0_real64
        !> Optional piecewise-linear history of the prescribed temperature:
        !> times since the initial state, s, ascending, and the values, C.
        !> Unallocated for a constant temperature.
        real(real64), allocatable :: schedule_times(:)
        real(real64), allocatable :: schedule_values(:)
    end type type_column_temperature_boundary

    !> One registered heat-transfer boundary: Newton cooling against a medium.
    type :: type_column_heat_transfer_boundary
        integer(int32) :: boundary_id = 0
        !> Heat transfer coefficient, W/(m2 K), and the medium's temperature, C.
        real(real64) :: coefficient = 0.0_real64
        real(real64) :: ambient = 0.0_real64
        !> Optional switch of the coefficient with the surface temperature: the
        !> coefficient below `switch_temperature`, C, reached over a smoothed
        !> step `switch_width` K wide. `coefficient` is the one above it. A
        !> negative `frozen_coefficient` means no switch.
        real(real64) :: frozen_coefficient = -1.0_real64
        real(real64) :: switch_temperature = 0.0_real64
        real(real64) :: switch_width = 0.0_real64
        !> Optional piecewise-linear history of the medium's temperature:
        !> times since the initial state, s, ascending, and the values, C.
        !> Unallocated for a constant medium.
        real(real64), allocatable :: schedule_times(:)
        real(real64), allocatable :: schedule_values(:)
    end type type_column_heat_transfer_boundary

    !> The coupled column solver.
    type :: type_mixed_column
        !> The solve DM, a clone of the caller's mesh carrying the three fields.
        DM :: dm
        !> The DM the cell-wise and accepted-state auxiliary data lives on.
        DM :: aux_dm
        PetscDS :: ds
        Vec :: solution
        Vec :: old_solution
        Vec :: local_solution
        Vec :: auxiliary
        Vec :: accepted_auxiliary
        !> The auxiliary state as it was before the last accepted step.
        !>
        !> `accepted_auxiliary` is republished from the new state the moment a
        !> step is accepted, so by the time anything outside the solver looks at
        !> it, it is the state the step ended at and not the one it started
        !> from. An audit differencing a backward-Euler step needs the state it
        !> started from, and this is where that state is kept.
        Vec :: previous_auxiliary
        SNES :: snes
        PetscScalar :: constants(COLUMN_CONSTANTS)
        !> One entry per material of the mesh, in closure-slot order.
        type(type_mixed_frozen_material), allocatable :: materials(:)
        !> The saturated conductivity tensor of each material, m/s.
        real(real64), allocatable :: conductivity(:, :, :)
        !> The registered prescribed-temperature boundaries.
        type(type_column_temperature_boundary) :: temperatures(MAX_COLUMN_BOUNDARIES)
        integer(int32) :: num_temperatures = 0
        !> The registered heat-transfer boundaries and the auxiliary vector each
        !> of them is evaluated through.
        type(type_column_heat_transfer_boundary) :: heat_transfers(MAX_COLUMN_BOUNDARIES)
        Vec :: heat_auxiliary(MAX_COLUMN_BOUNDARIES)
        integer(int32) :: num_heat_transfers = 0
        !> The "Face Sets" values the water boundary is closed on.
        integer(int32) :: closed_ids(MAX_COLUMN_BOUNDARIES) = 0
        integer(int32) :: num_closed = 0
        !> The label the prescribed temperatures are keyed by.
        character(len=32) :: thermal_label = "Thermal Sets"
        real(real64) :: flux_scale = 1.0_real64
        real(real64) :: pressure_scale = 1.0_real64
        real(real64) :: length_scale = 1.0_real64
        real(real64) :: temperature_scale = 1.0_real64
        real(real64) :: thermal_scale = 1.0_real64
        real(real64) :: density = 1000.0_real64
        real(real64) :: gravity(MIXED_DIMENSION) = 0.0_real64
        real(real64) :: dt = 1.0_real64
        !> The length of the last step that was accepted, s.
        !>
        !> `dt` is the length of the last step that was attempted, which after a
        !> refusal is not the one the accepted state was reached by.
        real(real64) :: accepted_dt = 0.0_real64
        !> Time of the accepted state since the initial state, s.
        real(real64) :: time = 0.0_real64
        !> Where the per-step history is written, or empty for none.
        character(len=512) :: history_path = ""
        integer(int32) :: history_unit = -1
        integer(int32) :: history_rows = 0
        !> The adaptive controller and the step it would attempt next, s.
        type(type_mixed_column_step_control) :: step_control
        real(real64) :: next_step = 0.0_real64
        !> Whether the run's constitutive laws may take part in vapour storage.
        logical :: enable_vapor = .false.
        !> Whether the pore pressure is in H1 (the primal form) instead of the
        !> RT0/P0 pair; fixed before `initialize`.
        logical :: primal = .false.
        !> Whether the primal flux takes the nodal liquid pressure and the
        !> geometric-mean conductivity instead of the quadrature-point chain rule.
        logical :: nodal_flux = .false.
        !> Zero-based field number of the pore pressure in the active layout.
        PetscInt :: pressure_field = COLUMN_FIELD_PRESSURE
        !> Water each cell's own flux has carried out through each of its
        !> facets since the initial state, m2 per unit thickness: the time
        !> integral of \(Q_{Kf}\) over the accepted steps, facets in
        !> `DMPlexGetCone` order, cells by local point number.
        real(real64), allocatable :: transported(:, :)
        !> Each cell's water \(\int_K W\) at the accepted state, m2 per unit
        !> thickness, cells by local point number.
        real(real64), allocatable :: cell_water(:)
        !> Mass balance of the last accepted step, with the method's own flux:
        !> the column's water change over its water (`system_balance`), and the
        !> largest \(|\Delta W_K+\Delta t\sum_f Q_{Kf}|/W_K\) of any cell
        !> (`element_balance`). Rank-local on several ranks.
        real(real64) :: system_balance = 0.0_real64
        real(real64) :: element_balance = 0.0_real64
        logical :: laws_installed = .false.
        logical :: ready = .false.
    contains
        ! ---- Lifecycle ----
        procedure, public, pass :: initialize => initialize_mixed_column
        procedure, public, pass :: destroy => destroy_mixed_column
        ! ---- Mutator ----
        procedure, public, pass :: set_materials => set_column_materials
        procedure, public, pass :: set_gravity => set_column_gravity
        procedure, public, pass :: set_scales => set_column_scales
        procedure, public, pass :: add_closed_boundary
        procedure, public, pass :: add_temperature_boundary
        procedure, public, pass :: add_heat_transfer_boundary
        procedure, public, pass :: set_initial_state => set_column_initial_state
        procedure, public, pass :: set_boundary_temperature
        procedure, public, pass :: set_temperature_schedule => set_column_temperature_schedule
        procedure, public, pass :: set_ambient_schedule => set_column_ambient_schedule
        procedure, public, pass :: set_heat_transfer_switch => set_column_heat_transfer_switch
        procedure, public, pass :: set_step_control => set_column_step_control
        procedure, public, pass :: set_step_history => set_column_step_history
        procedure, public, pass :: set_hydraulic_space => set_column_hydraulic_space
        ! ---- Algorithm ----
        procedure, public, pass :: solve_step => solve_column_step
        procedure, public, pass :: advance => advance_column
        procedure, public, pass :: compute_report => compute_column_report
        procedure, public, pass :: write_profile => write_column_profile
        ! ---- Getter ----
        procedure, public, pass :: get_dm => get_column_dm
    end type type_mixed_column

    ! PetscFPTrapPush/Pop have no Fortran .mod interface in this PETSc build,
    ! only the auto-generated C stubs every PETSc Fortran call ultimately
    ! goes through; see
    ! third_party/petsc-src/arch-linux-c-opt/ftn/sys/error/fpf.c for the
    ! exact (by-reference int, by-reference int) signature this mirrors.
    ! Shared by every routine below that constructs a PetscFE/PetscDualSpace:
    ! see initialize_mixed_column's use of it for why.
    integer(c_int), parameter :: PETSC_FP_TRAP_OFF_LOCAL = 0_c_int
    interface
        subroutine petscfptrappush_(trap, trap_ierr) bind(c, name="petscfptrappush_")
            import :: c_int
            integer(c_int), intent(in) :: trap
            integer(c_int), intent(out) :: trap_ierr
        end subroutine petscfptrappush_
        subroutine petscfptrappop_(trap_ierr) bind(c, name="petscfptrappop_")
            import :: c_int
            integer(c_int), intent(out) :: trap_ierr
        end subroutine petscfptrappop_
    end interface

contains

    !> Choose the hydraulic space, before `initialize`.
    !>
    !> `mixed` (the default) is the RT0 flux with the P0 pressure, conservative
    !> cell by cell. `primal` is the continuous Galerkin form: the pressure in
    !> H1 with the same element as the temperature, and no flux unknown.
    !> `primal_nodal` is the same space with the Darcy term built from nodal
    !> liquid pressures and a geometric-mean conductivity.
    !>
    !> \(O(1)\). Fails by `error stop` on any other name, or on a solver whose
    !> fields already exist.
    subroutine set_column_hydraulic_space(self, space)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> `mixed`, `primal` or `primal_nodal`.
        character(len=*), intent(in) :: space

        if (self%ready) error stop "The hydraulic space of a column is chosen before it is initialised"
        self%nodal_flux = .false.
        select case (trim(space))
        case ("mixed")
            self%primal = .false.
            self%pressure_field = COLUMN_FIELD_PRESSURE
        case ("primal", "primal_nodal")
            self%primal = .true.
            self%nodal_flux = (trim(space) == "primal_nodal")
            self%pressure_field = COLUMN_PRIMAL_FIELD_PRESSURE
        case default
            error stop "The hydraulic space of a column is mixed, primal or primal_nodal"
        end select
    end subroutine set_column_hydraulic_space

    !> The materials of the mesh, keyed by their "Cell Sets" value.
    !>
    !> `conductivity(:,:,m)` is material `m`'s saturated hydraulic conductivity
    !> tensor in m/s; the relative conductivity of the state divides it inside
    !> the kernels.
    subroutine set_column_materials(self, materials, conductivity)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        type(type_mixed_frozen_material), intent(in) :: materials(:)
        real(real64), intent(in) :: conductivity(:, :, :)

        if (size(conductivity, 3) /= size(materials)) then
            error stop "The mixed column needs one conductivity tensor per material"
        end if
        if (allocated(self%materials)) deallocate (self%materials)
        if (allocated(self%conductivity)) deallocate (self%conductivity)
        allocate (self%materials(size(materials)))
        allocate (self%conductivity(MIXED_DIMENSION, MIXED_DIMENSION, size(materials)))
        self%materials = materials
        self%conductivity = conductivity
    end subroutine set_column_materials

    !> Gravity in the intrinsic plane, m/s2; downward is negative in z.
    subroutine set_column_gravity(self, gravity, density)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        real(real64), intent(in) :: gravity(MIXED_DIMENSION)
        real(real64), intent(in) :: density

        self%gravity = gravity
        self%density = density
    end subroutine set_column_gravity

    !> The units the dimensionless unknowns are measured in.
    subroutine set_column_scales(self, flux_scale, pressure_scale, length_scale, &
                                 temperature_scale, thermal_scale)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> Flux unit, m/s; pressure unit, Pa; the conservation row's length, m.
        real(real64), intent(in) :: flux_scale, pressure_scale, length_scale
        !> Temperature unit, K.
        real(real64), intent(in) :: temperature_scale
        !> The heat row's divisor \(R_T\), W/m3.
        real(real64), intent(in) :: thermal_scale

        if (flux_scale <= 0.0_real64 .or. pressure_scale <= 0.0_real64 .or. length_scale <= 0.0_real64 &
            .or. temperature_scale <= 0.0_real64 .or. thermal_scale <= 0.0_real64) then
            error stop "The mixed column scales have to be positive"
        end if
        self%flux_scale = flux_scale
        self%pressure_scale = pressure_scale
        self%length_scale = length_scale
        self%temperature_scale = temperature_scale
        self%thermal_scale = thermal_scale
    end subroutine set_column_scales

    !> Close one "Face Sets" boundary to water, before `initialize`.
    subroutine add_closed_boundary(self, boundary_id)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        integer(int32), intent(in) :: boundary_id

        if (self%num_closed >= MAX_COLUMN_BOUNDARIES) then
            error stop "More closed boundaries than the mixed column can hold"
        end if
        self%num_closed = self%num_closed + 1
        self%closed_ids(self%num_closed) = boundary_id
    end subroutine add_closed_boundary

    !> Prescribe the temperature of one boundary of the thermal label.
    subroutine add_temperature_boundary(self, boundary_id, temperature)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        integer(int32), intent(in) :: boundary_id
        !> The prescribed temperature, C.
        real(real64), intent(in) :: temperature

        if (self%num_temperatures >= MAX_COLUMN_BOUNDARIES) then
            error stop "More prescribed temperatures than the mixed column can hold"
        end if
        self%num_temperatures = self%num_temperatures + 1
        self%temperatures(self%num_temperatures)%boundary_id = boundary_id
        self%temperatures(self%num_temperatures)%temperature = temperature
    end subroutine add_temperature_boundary

    !> Exchange heat with a medium across one boundary of the thermal label.
    !>
    !> Newton cooling, \(q=h(T-T_\infty)\) out of the surface. This is the
    !> condition the experiment is run under -- a coolant or a room on the other
    !> side of a wall, not a temperature imposed on the soil itself -- and it
    !> differs from `add_temperature_boundary` in exactly that: the surface
    !> temperature is part of the solution and approaches \(T_\infty\) at a rate
    !> the coefficient sets.
    !>
    !> \[ q\cdot\mathbf n = h\,(T-T_\infty)\quad\text{on the boundary.} \]
    !>
    !> Must be called before `initialize`, which is where the condition is
    !> registered. \(O(1)\). Fails by `error stop` on a negative coefficient or
    !> once `MAX_COLUMN_BOUNDARIES` of them have been registered.
    subroutine add_heat_transfer_boundary(self, boundary_id, coefficient, ambient)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        integer(int32), intent(in) :: boundary_id
        !> Heat transfer coefficient, W/(m2 K); not negative.
        real(real64), intent(in) :: coefficient
        !> The temperature of the medium beyond the boundary, C.
        real(real64), intent(in) :: ambient

        if (self%num_heat_transfers >= MAX_COLUMN_BOUNDARIES) then
            error stop "More heat-transfer boundaries than the mixed column can hold"
        end if
        if (coefficient < 0.0_real64) then
            error stop "A heat-transfer coefficient of the mixed column cannot be negative"
        end if
        self%num_heat_transfers = self%num_heat_transfers + 1
        self%heat_transfers(self%num_heat_transfers)%boundary_id = boundary_id
        self%heat_transfers(self%num_heat_transfers)%coefficient = coefficient
        self%heat_transfers(self%num_heat_transfers)%ambient = ambient
    end subroutine add_heat_transfer_boundary

    !> Change what one registered boundary prescribes, between steps.
    !>
    !> The essential callback reads the value through its context, so a run
    !> cools the top by changing this rather than by rebuilding the solver.
    subroutine set_boundary_temperature(self, boundary_id, temperature)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        integer(int32), intent(in) :: boundary_id
        real(real64), intent(in) :: temperature

        integer(int32) :: b

        do b = 1, self%num_temperatures
            if (self%temperatures(b)%boundary_id /= boundary_id) cycle
            self%temperatures(b)%temperature = temperature
            call set_column_temperature_value(b, temperature, self%temperature_scale)
            return
        end do
        error stop "No prescribed-temperature boundary of the mixed column carries that identifier"
    end subroutine set_boundary_temperature

    !> Let one prescribed-temperature boundary follow a history in time.
    !>
    !> \[
    !>   T_b(t)=T_k+\frac{t-t_k}{t_{k+1}-t_k}\,(T_{k+1}-T_k),\qquad t_k\le t\le t_{k+1},
    !> \]
    !>
    !> held at the first value before \(t_1\) and at the last after \(t_n\), with
    !> \(t\) counted from the initial state. A backward-Euler step to
    !> \(t^{n+1}\) is solved with \(T_b(t^{n+1})\), so the boundary is part of
    !> the implicit step and not lagged. The boundary must have been registered
    !> with `add_temperature_boundary`. \(O(n)\) storage, \(O(n)\) per step.
    !>
    !> Fails by `error stop` on arrays of different or zero length, times that
    !> are not strictly ascending, or an identifier no prescribed-temperature
    !> boundary carries.
    subroutine set_column_temperature_schedule(self, boundary_id, times, values)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        integer(int32), intent(in) :: boundary_id
        !> \(t_k\), s since the initial state; strictly ascending.
        real(real64), intent(in) :: times(:)
        !> \(T_k\), C; one per time.
        real(real64), intent(in) :: values(:)

        integer(int32) :: b

        call check_schedule(times, values)
        do b = 1, self%num_temperatures
            if (self%temperatures(b)%boundary_id /= boundary_id) cycle
            self%temperatures(b)%schedule_times = times
            self%temperatures(b)%schedule_values = values
            self%temperatures(b)%temperature = schedule_value(times, values, self%time)
            if (self%ready) call set_column_temperature_value(b, self%temperatures(b)%temperature, &
                                                              self%temperature_scale)
            return
        end do
        error stop "No prescribed-temperature boundary of the mixed column carries that identifier"
    end subroutine set_column_temperature_schedule

    !> Let the medium beyond one heat-transfer boundary follow a history in time.
    !>
    !> \[
    !>   q\cdot\mathbf n=h\,\bigl(T-T_\infty(t)\bigr),
    !> \]
    !>
    !> with \(T_\infty(t)\) the same piecewise-linear interpolation as
    !> `set_temperature_schedule`, held at its end values outside the table and
    !> evaluated at \(t^{n+1}\) for a step to \(t^{n+1}\). The coefficient
    !> stays the registered one. The boundary must have been registered with
    !> `add_heat_transfer_boundary`. \(O(n)\) storage; \(O(n)\) plus one pass
    !> over the auxiliary vector per step.
    !>
    !> Fails by `error stop` on arrays of different or zero length, times that
    !> are not strictly ascending, or an identifier no heat-transfer boundary
    !> carries.
    subroutine set_column_ambient_schedule(self, boundary_id, times, values)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        integer(int32), intent(in) :: boundary_id
        !> \(t_k\), s since the initial state; strictly ascending.
        real(real64), intent(in) :: times(:)
        !> \(T_{\infty,k}\), C; one per time.
        real(real64), intent(in) :: values(:)

        integer(int32) :: b
        PetscErrorCode :: ierr

        call check_schedule(times, values)
        do b = 1, self%num_heat_transfers
            if (self%heat_transfers(b)%boundary_id /= boundary_id) cycle
            self%heat_transfers(b)%schedule_times = times
            self%heat_transfers(b)%schedule_values = values
            self%heat_transfers(b)%ambient = schedule_value(times, values, self%time)
            if (self%ready) then
                ierr = 0
                call fill_column_aux_field(self%aux_dm, self%heat_auxiliary(b), CA_HEAT_AMBIENT, &
                                           self%heat_transfers(b)%ambient, ierr)
                if (ierr /= 0) error stop "PETSc could not update a heat-transfer boundary's medium"
            end if
            return
        end do
        error stop "No heat-transfer boundary of the mixed column carries that identifier"
    end subroutine set_column_ambient_schedule

    !> Let one heat-transfer coefficient switch with the surface temperature.
    !>
    !> \[
    !>   h(T)=h_f+(h_w-h_f)\,S\Bigl(\frac{T-T_s}{w}+\frac12\Bigr),\qquad
    !>   S(u)=u^2(3-2u)\ \text{on}\ [0,1],
    !> \]
    !>
    !> with \(h_w\) the coefficient the boundary was registered with and \(T\)
    !> the surface temperature: \(h_f\) below \(T_s-w/2\), \(h_w\) above
    !> \(T_s+w/2\), and a C1 step between them, so the Newton iteration never
    !> meets a jump. Must be called before `initialize`. \(O(1)\).
    !>
    !> Fails by `error stop` on a negative coefficient, a width that is not
    !> positive, or an identifier no heat-transfer boundary carries.
    subroutine set_column_heat_transfer_switch(self, boundary_id, frozen_coefficient, &
                                               switch_temperature, switch_width)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        integer(int32), intent(in) :: boundary_id
        !> \(h_f\), W/(m2 K); not negative.
        real(real64), intent(in) :: frozen_coefficient
        !> \(T_s\), C.
        real(real64), intent(in) :: switch_temperature
        !> \(w\), K; positive.
        real(real64), intent(in) :: switch_width

        integer(int32) :: b

        if (frozen_coefficient < 0.0_real64) then
            error stop "A heat-transfer coefficient of the mixed column cannot be negative"
        end if
        if (.not. (switch_width > 0.0_real64)) then
            error stop "The width of a heat-transfer switch has to be positive"
        end if
        do b = 1, self%num_heat_transfers
            if (self%heat_transfers(b)%boundary_id /= boundary_id) cycle
            self%heat_transfers(b)%frozen_coefficient = frozen_coefficient
            self%heat_transfers(b)%switch_temperature = switch_temperature
            self%heat_transfers(b)%switch_width = switch_width
            return
        end do
        error stop "No heat-transfer boundary of the mixed column carries that identifier"
    end subroutine set_column_heat_transfer_switch

    !> Refuse a schedule that is not one value per strictly ascending time.
    subroutine check_schedule(times, values)
        implicit none
        real(real64), intent(in) :: times(:)
        real(real64), intent(in) :: values(:)

        if (size(times) < 1 .or. size(times) /= size(values)) then
            error stop "A boundary schedule needs one value per time and at least one of each"
        end if
        if (size(times) > 1) then
            if (any(times(2:) <= times(:size(times) - 1))) then
                error stop "The times of a boundary schedule have to be strictly ascending"
            end if
        end if
    end subroutine check_schedule

    !> The piecewise-linear value of a schedule at `time`, held at its ends.
    pure function schedule_value(times, values, time) result(value)
        implicit none
        !> Strictly ascending times, s; at least one.
        real(real64), intent(in) :: times(:)
        !> One value per time.
        real(real64), intent(in) :: values(:)
        !> Time since the initial state, s.
        real(real64), intent(in) :: time
        real(real64) :: value

        integer(int32) :: k, n
        real(real64) :: weight

        n = int(size(times), int32)
        if (time <= times(1)) then
            value = values(1)
            return
        end if
        if (time >= times(n)) then
            value = values(n)
            return
        end if
        do k = 1, n - 1
            if (time <= times(k + 1)) exit
        end do
        weight = (time - times(k))/(times(k + 1) - times(k))
        value = values(k) + weight*(values(k + 1) - values(k))
    end function schedule_value

    !> Put every scheduled boundary at its value for the time `time`.
    subroutine apply_boundary_schedules(self, time, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> Time since the initial state, s.
        real(real64), intent(in) :: time
        PetscErrorCode, intent(inout) :: ierr

        integer(int32) :: b

        do b = 1, self%num_temperatures
            if (.not. allocated(self%temperatures(b)%schedule_times)) cycle
            self%temperatures(b)%temperature = schedule_value(self%temperatures(b)%schedule_times, &
                                                              self%temperatures(b)%schedule_values, time)
            call set_column_temperature_value(b, self%temperatures(b)%temperature, self%temperature_scale)
        end do
        do b = 1, self%num_heat_transfers
            if (.not. allocated(self%heat_transfers(b)%schedule_times)) cycle
            self%heat_transfers(b)%ambient = schedule_value(self%heat_transfers(b)%schedule_times, &
                                                            self%heat_transfers(b)%schedule_values, time)
            if (.not. self%ready) cycle
            call fill_column_aux_field(self%aux_dm, self%heat_auxiliary(b), CA_HEAT_AMBIENT, &
                                       self%heat_transfers(b)%ambient, ierr)
            if (ierr /= 0) return
        end do
    end subroutine apply_boundary_schedules

    !> Build the three-field system, its coefficients, boundaries and solver.
    !>
    !> \[
    !>   R_T(w)=\Bigl(\frac{U-U^n}{\Delta t},w\Bigr)+(\lambda\nabla T,\nabla w),\quad
    !>   R_q(\mathbf v)=(\mathbf\Lambda^{-1}\mathbf q_l-\rho_l\mathbf g,\mathbf v)-(p_l,\nabla\cdot\mathbf v),
    !> \]
    !> \[
    !>   R_p(\omega)=\Bigl(\frac{W-W^n}{\Delta t}+\nabla\cdot\mathbf q_l-s,\omega\Bigr).
    !> \]
    !>
    !> Assumes an interpolated DMPlex whose dimension equals its coordinate
    !> dimension, one "Cell Sets" material per cell, and a symmetric positive
    !> definite saturated conductivity for each material. Cells are either all
    !> triangles or all axis-aligned rectangles of positive area; the second is
    !> a reference geometry for comparing the two element families on one
    !> problem, not general quadrilateral support. Materials, gravity, scales
    !> and boundaries must already be set.
    !>
    !> Guarantees a nodal temperature, an RT0 flux with one degree of freedom
    !> per facet and a P0 pressure with one per cell -- P1/RT0/P0 on a triangle,
    !> Q1 and the tensor-product RT0 of PETSc's `snes/tutorials/ex24.c` on a
    !> rectangle. It does not guarantee that the resulting system is
    !> nonsingular. Setup costs \(O(N)\) in cells plus PETSc's section
    !> construction.
    !>
    !> Fails by `error stop` naming the refusing condition: a mesh of the wrong
    !> dimension or element, a quadrilateral that is not an axis-aligned
    !> rectangle, a cell without a material, a conductivity that is not positive
    !> definite, or a closure that cannot evaluate its own reference state.
    subroutine initialize_mixed_column(self, dm, ierr)
        implicit none
        class(type_mixed_column), intent(inout), target :: self
        !> An interpolated DMPlex; the solver clones it and does not own it.
        DM, intent(in) :: dm
        PetscErrorCode, intent(inout) :: ierr

        PetscFE :: fe_temperature, fe_flux, fe_pressure
        PetscObject :: discretisation
        PetscDS :: ds
        PetscInt :: dimension, coordinate_dimension
        PetscBool :: is_simplex
        type(c_ptr) :: dm_handle, ds_handle
        integer(c_int) :: cerr

        if (.not. allocated(self%materials)) error stop "The mixed column solver has no materials"

        call DMGetDimension(dm, dimension, ierr)
        if (ierr /= 0) error stop "The mesh has no dimension"
        call DMGetCoordinateDim(dm, coordinate_dimension, ierr)
        if (ierr /= 0) error stop "The mesh has no coordinate dimension"
        if (int(dimension, int32) /= MIXED_DIMENSION .or. dimension /= coordinate_dimension) then
            error stop "The mixed column solver needs an intrinsic two-dimensional mesh"
        end if
        call DMPlexIsSimplex(dm, is_simplex, ierr)
        if (ierr /= 0) error stop "PETSc could not read the cell shape of the mesh"
        if (is_simplex .neqv. PETSC_TRUE) then
            ! The reconstruction this solver audits itself with inverts one
            ! affine map per cell and pushes the RT basis forward with it, so a
            ! cell whose map is not affine would be reported through a geometry
            ! it does not have. That rules out the warped cells; the gate below
            ! is narrower still and takes only axis-aligned rectangles of
            ! positive area, which is a supported-geometry policy rather than
            ! the affine condition itself.
            call verify_rectangular_cells(dm, ierr)
            if (ierr /= 0) then
                error stop "The mixed column solver supports triangles and axis-aligned rectangles only"
            end if
        end if

        call DMClone(dm, self%dm, ierr)
        if (ierr /= 0) error stop "PETSc could not clone the mesh for the coupled column"

        call set_column_element_options(is_simplex, self%primal, ierr)
        if (ierr /= 0) error stop "PETSc rejected the coupled element options"

        ! Building the trimmed (RT0) dual space's facet subspace touches a
        ! transient, mathematically inert floating point exception inside
        ! PETSc's own PetscDualSpaceLagrangeCreateSimplexNodeMat -- confirmed
        ! by reproducing it in isolation (a bare PetscFECreateDefault call
        ! with these exact options, nothing else) and observing it appear
        ! only once IEEE traps are unmasked (feenableexcept), never
        ! otherwise; PETSc itself runs with traps masked by default, which is
        ! why the library's own test suite (src/snes/tutorials/ex24.c,
        ! suffix 2d_rt0_quad/2d_rt0_tri, whose options this solver's
        ! set_simplex_flux_options/set_tensor_flux_options mirror) never
        ! trips it. This project's Debug build unmasks hardware FP traps,
        ! which turns that inert event fatal before the space construction
        ! can finish. Push the trap off for exactly these three definition
        ! calls -- not any solve -- the same way PETSc's own code crosses
        ! routines it knows do this (PetscFPTrapPush's own doc comment); a
        ! real division by zero anywhere else, including immediately after
        ! Pop, still traps as usual.
        call petscfptrappush_(PETSC_FP_TRAP_OFF_LOCAL, cerr)
        call PetscFECreateDefault(PETSC_COMM_WORLD, dimension, 1_PETSC_INT_KIND, is_simplex, &
                                  "temperature_", QUADRATURE_ORDER, fe_temperature, ierr)
        if (ierr /= 0) error stop "PETSc could not build the nodal temperature basis"
        if (.not. self%primal) then
            call PetscFECreateDefault(PETSC_COMM_WORLD, dimension, 2_PETSC_INT_KIND, is_simplex, &
                                      "colflux_", QUADRATURE_ORDER, fe_flux, ierr)
            if (ierr /= 0) error stop "PETSc could not build the RT0 flux basis"
        end if
        call PetscFECreateDefault(PETSC_COMM_WORLD, dimension, 1_PETSC_INT_KIND, is_simplex, &
                                  "colpressure_", QUADRATURE_ORDER, fe_pressure, ierr)
        if (ierr /= 0) error stop "PETSc could not build the pressure basis"
        call petscfptrappop_(cerr)

        ! One quadrature for all the fields, so every block of the monolithic
        ! system is integrated at the same points.
        if (.not. self%primal) call PetscFECopyQuadrature(fe_temperature, fe_flux, ierr)
        if (ierr == 0) call PetscFECopyQuadrature(fe_temperature, fe_pressure, ierr)
        if (ierr /= 0) error stop "PETSc could not share the quadrature between the coupled fields"

        discretisation = PetscObjectCast(fe_temperature)
        call PetscObjectSetName(discretisation, "temperature", ierr)
        if (ierr == 0) call DMSetField(self%dm, COLUMN_FIELD_TEMPERATURE, PETSC_NULL_DMLABEL, discretisation, ierr)
        if (ierr /= 0) error stop "PETSc could not attach the temperature field"

        if (.not. self%primal) then
            discretisation = PetscObjectCast(fe_flux)
            call PetscObjectSetName(discretisation, "flux", ierr)
            if (ierr == 0) call DMSetField(self%dm, COLUMN_FIELD_FLUX, PETSC_NULL_DMLABEL, discretisation, ierr)
            if (ierr /= 0) error stop "PETSc could not attach the coupled flux field"
        end if

        discretisation = PetscObjectCast(fe_pressure)
        call PetscObjectSetName(discretisation, "pressure", ierr)
        if (ierr == 0) call DMSetField(self%dm, self%pressure_field, PETSC_NULL_DMLABEL, discretisation, ierr)
        if (ierr /= 0) error stop "PETSc could not attach the coupled pressure field"

        ! The temperature element outlives the others here: it carries the
        ! quadrature every auxiliary element below is made to share, and the
        ! options database may have moved that rule off QUADRATURE_ORDER.
        if (.not. self%primal) call PetscFEDestroy(fe_flux, ierr)
        call PetscFEDestroy(fe_pressure, ierr)

        call DMCreateDS(self%dm, ierr)
        if (ierr == 0) call DMGetDS(self%dm, ds, ierr)
        if (ierr /= 0) error stop "PETSc could not build the coupled discrete system"
        self%ds = ds

        dm_handle = transfer(self%dm%v, dm_handle)
        ds_handle = transfer(ds%v, ds_handle)

        if (self%primal) then
            call register_primal_kernels(ds_handle, self%nodal_flux, cerr)
        else
            call register_column_kernels(ds_handle, cerr)
        end if
        if (cerr /= 0) error stop "PETSc rejected a coupled column kernel"

        call store_column_constants(self, ierr)
        if (ierr /= 0) error stop "PETSc rejected the coupled run constants"

        ! The mixed path's intrinsic plane is the x-z one, which is also the
        ! only geometry this solver accepts above.
        call install_frozen_laws(self%materials, COMP_TYPES%XZ_2D%ID, self%enable_vapor)
        self%laws_installed = .true.

        call build_column_auxiliary(self, dimension, is_simplex, fe_temperature, ierr)
        if (ierr /= 0) error stop "PETSc could not build the coupled auxiliary fields"
        call PetscFEDestroy(fe_temperature, ierr)

        call add_column_boundaries(self, dm_handle, ierr)
        if (ierr /= 0) error stop "PETSc rejected a coupled column boundary"

        call DMCreateGlobalVector(self%dm, self%solution, ierr)
        if (ierr == 0) call DMCreateGlobalVector(self%dm, self%old_solution, ierr)
        if (ierr == 0) call DMCreateLocalVector(self%dm, self%local_solution, ierr)
        if (ierr /= 0) error stop "PETSc could not create the coupled vectors"
        call VecZeroEntries(self%solution, ierr)
        if (ierr == 0) call VecZeroEntries(self%old_solution, ierr)

        call SNESCreate(PETSC_COMM_WORLD, self%snes, ierr)
        if (ierr == 0) call SNESSetDM(self%snes, self%dm, ierr)
        if (ierr == 0) call DMPlexSetSNESLocalFEM(self%dm, PETSC_FALSE, PETSC_NULL_INTEGER, ierr)
        if (ierr == 0) call set_column_linesearch_default(ierr)
        if (ierr == 0) call SNESSetFromOptions(self%snes, ierr)
        if (ierr /= 0) error stop "PETSc could not set up the coupled nonlinear solver"

        ! The nodal flux is assembled outside PETSc's element loop, which sees
        ! quadrature points only, and added to what that loop computes.
        if (self%nodal_flux) then
            call install_nodal_flux(self%dm, self%aux_dm, self%auxiliary, self%pressure_scale, &
                                    self%temperature_scale, self%flux_scale, self%length_scale, &
                                    self%density, self%gravity, ierr)
            if (ierr /= 0) error stop "PETSc rejected the nodal column flux"
        end if

        block
            PetscInt :: c_start, c_end

            call DMPlexGetHeightStratum(self%dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
            if (ierr /= 0) error stop "PETSc could not count the column's cells"
            allocate (self%transported(MAX_CELL_FACETS, c_start:c_end - 1))
            self%transported = 0.0_real64
            allocate (self%cell_water(c_start:c_end - 1))
            self%cell_water = 0.0_real64
        end block

        self%ready = .true.
    end subroutine initialize_mixed_column

    !> Refuse any cell that is not an axis-aligned rectangle of positive area.
    !>
    !> Every cell must satisfy, for its closure vertices \(a,b,c,d\) in cyclic
    !> order with \(\mathbf e_1=b-a\) and \(\mathbf e_2=d-a\),
    !>
    !> \[
    !>   c=a+\mathbf e_1+\mathbf e_2,\qquad \mathbf e_1\cdot\mathbf e_2=0,
    !>   \qquad e_{1x}e_{2z}-e_{1z}e_{2x}>0,
    !> \]
    !>
    !> with each edge parallel to a coordinate axis. The first condition makes
    !> the map from the reference square affine, which is what the audit's
    !> `DMPlexComputeCellGeometryAffineFEM` and the element loop's
    !> `isAffine` pushforward both assume; the rest keeps the Jacobian diagonal
    !> and orientation-preserving, so the intrinsic x-z frame of the closure is
    !> the cell's own. A warped or inverted cell is refused here rather than
    !> reported through a geometry it does not have.
    !>
    !> \(O(N)\) in local cells and \(O(1)\) memory. Fails by returning a nonzero
    !> `ierr`, which is the first failing PETSc query or `1` for a cell that is
    !> not such a rectangle.
    subroutine verify_rectangular_cells(dm, ierr)
        implicit none
        DM, intent(in) :: dm
        PetscErrorCode, intent(inout) :: ierr

        PetscScalar, pointer :: cell_coordinates(:)
        PetscSection :: coordinate_section
        Vec :: coordinates
        DM :: coordinate_dm
        PetscInt :: cell, c_start, c_end, cone_size, closure_size
        real(real64) :: corners(MIXED_DIMENSION, 4)
        real(real64) :: edge_a(MIXED_DIMENSION), edge_b(MIXED_DIMENSION)
        real(real64) :: tolerance, area
        integer(int32) :: v, d
        logical :: acceptable

        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr == 0) call DMGetCoordinateDM(dm, coordinate_dm, ierr)
        if (ierr == 0) call DMGetCoordinatesLocal(dm, coordinates, ierr)
        if (ierr == 0) call DMGetCoordinateSection(dm, coordinate_section, ierr)
        if (ierr /= 0) return

        do cell = c_start, c_end - 1
            call DMPlexGetConeSize(dm, cell, cone_size, ierr)
            if (ierr /= 0) return
            if (cone_size /= 4_PETSC_INT_KIND) then
                ierr = 1
                return
            end if

            nullify (cell_coordinates)
            closure_size = 0_PETSC_INT_KIND
            call DMPlexVecGetClosure(coordinate_dm, coordinate_section, coordinates, cell, &
                                     closure_size, cell_coordinates, ierr)
            if (ierr /= 0) return
            acceptable = (closure_size == 4*MIXED_DIMENSION)
            if (acceptable) then
                do v = 1, 4
                    do d = 1, MIXED_DIMENSION
                        corners(d, v) = real(cell_coordinates((v - 1)*MIXED_DIMENSION + d), real64)
                    end do
                end do
            end if
            call DMPlexVecRestoreClosure(coordinate_dm, coordinate_section, coordinates, cell, &
                                         closure_size, cell_coordinates, ierr)
            if (ierr /= 0) return
            if (.not. acceptable) then
                ierr = 1
                return
            end if

            edge_a = corners(:, 2) - corners(:, 1)
            edge_b = corners(:, 4) - corners(:, 1)
            tolerance = 1.0e-12_real64*max(norm2(edge_a), norm2(edge_b))
            area = edge_a(1)*edge_b(2) - edge_a(2)*edge_b(1)

            acceptable = (maxval(abs(corners(:, 3) - (corners(:, 1) + edge_a + edge_b))) <= tolerance)
            acceptable = acceptable .and. (area > 0.0_real64)
            ! Axis aligned: one edge runs along x and the other along z.
            acceptable = acceptable .and. &
                         ((abs(edge_a(2)) <= tolerance .and. abs(edge_b(1)) <= tolerance) .or. &
                          (abs(edge_a(1)) <= tolerance .and. abs(edge_b(2)) <= tolerance))
            if (.not. acceptable) then
                ierr = 1
                return
            end if
        end do
    end subroutine verify_rectangular_cells

    !> Ask for the natural-norm line search unless the run named another.
    !>
    !> PETSc's default `bt` measures a trial step by \(\|F\|_2\). Across a
    !> freezing front that is the wrong measure: the latent term makes the heat
    !> rows of one cell enormous next to the rows of its unfrozen neighbours, so
    !> a step that moves the front reduces the error while raising the residual
    !> norm, and `bt` answers by shrinking the step. Measured on the freezing
    !> column, `bt` cut \(\lambda\) to \(2\times10^{-3}\) and stalled at
    !> `DIVERGED_MAX_IT` after 400 iterations on a step that the same Newton
    !> direction undamped converged in seven. `nleqerr` is Deuflhard's
    !> affine-covariant criterion: it measures the trial step in the norm of the
    !> Newton correction itself, which row scaling cannot distort.
    !>
    !> Setting it here rather than in a caller keeps it a property of this
    !> coupled system; a command line still overrides it.
    subroutine set_column_linesearch_default(ierr)
        implicit none
        PetscErrorCode, intent(inout) :: ierr

        PetscBool :: already_set

        already_set = PETSC_FALSE
        call PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                                 "-snes_linesearch_type", already_set, ierr)
        if (ierr /= 0) return
        if (already_set .eqv. PETSC_TRUE) return
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-snes_linesearch_type", "nleqerr", ierr)
    end subroutine set_column_linesearch_default

    !> The options that turn PetscFECreateDefault into P1/Q1, RT0 and P0.
    !>
    !> The temperature and the pressure differ between the two cell shapes only
    !> in the reference cell `PetscFECreateDefault` is given, so their options
    !> are the same either way: `poly` of degree one is P1 on a triangle and Q1
    !> on a quadrilateral, because the tensor flag follows the reference cell.
    !> The flux does not: there is no trimmed simplex space on a tensor cell, so
    !> the rectangle takes the sum of two tensor spaces instead.
    !>
    !> The primal form asks for the pressure in the temperature's own space,
    !> degree one and continuous.
    subroutine set_column_element_options(is_simplex, primal, ierr)
        implicit none
        !> Whether the mesh's cells are simplices.
        PetscBool, intent(in) :: is_simplex
        !> Whether the pressure is in H1 rather than P0.
        logical, intent(in) :: primal
        PetscErrorCode, intent(inout) :: ierr

        character(len=8) :: degree, continuity

        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-temperature_petscspace_type", "poly", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-temperature_petscspace_degree", "1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-temperature_petscspace_components", "1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-temperature_petscdualspace_type", "lagrange", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-temperature_petscdualspace_order", "1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-temperature_petscdualspace_lagrange_continuity", "true", ierr)

        if (is_simplex .eqv. PETSC_TRUE) then
            call set_simplex_flux_options(ierr)
        else
            call set_tensor_flux_options(ierr)
        end if
        if (ierr /= 0) return

        degree = merge("1    ", "0    ", primal)
        continuity = merge("true ", "false", primal)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colpressure_petscspace_type", "poly", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colpressure_petscspace_degree", trim(degree), ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colpressure_petscspace_components", "1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colpressure_petscdualspace_type", "lagrange", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colpressure_petscdualspace_order", trim(degree), ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colpressure_petscdualspace_form_degree", "0", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colpressure_petscdualspace_lagrange_continuity", trim(continuity), ierr)
    end subroutine set_column_element_options

    !> RT0 on a triangle: the trimmed polynomial space of one-forms.
    !>
    !> The sum-space keys the tensor cell needs are cleared rather than left
    !> behind, because the options database outlives one solver and a stale key
    !> would silently redefine the next column's flux element.
    subroutine set_simplex_flux_options(ierr)
        implicit none
        PetscErrorCode, intent(inout) :: ierr

        call clear_tensor_flux_options(ierr)
        if (ierr /= 0) return

        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscspace_type", "ptrimmed", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscspace_degree", "0", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscspace_components", "2", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_petscspace_ptrimmed_form_degree", "-1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscdualspace_type", "lagrange", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscdualspace_order", "1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscdualspace_form_degree", "-1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_petscdualspace_lagrange_trimmed", "true", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_petscdualspace_lagrange_continuity", "true", ierr)
    end subroutine set_simplex_flux_options

    !> RT0 on a rectangle: \(Q_{1,0}\times Q_{0,1}\), as PETSc's own test builds it.
    !>
    !> \[
    !>   \mathbf{RT}_0(K)=\mathrm{span}\{(1,0),(x,0),(0,1),(0,z)\}
    !> \]
    !>
    !> on the reference square, which is not a trimmed simplex space and cannot
    !> be asked for as one. It is the concatenated sum of two tensor spaces --
    !> \(Q_1(x)\otimes Q_0(z)\) for the first component and \(Q_0(x)\otimes
    !> Q_1(z)\) for the second -- under the same trimmed dual space of one-forms
    !> the triangle uses, so the four functionals are still the four facet normal
    !> moments and the normal trace is still single valued. The keys are those of
    !> `src/snes/tutorials/ex24.c`'s `2d_rt0_quad` test with this field's prefix.
    subroutine set_tensor_flux_options(ierr)
        implicit none
        PetscErrorCode, intent(inout) :: ierr

        call PetscOptionsClearValue(PETSC_NULL_OPTIONS, "-colflux_petscspace_ptrimmed_form_degree", ierr)
        if (ierr /= 0) return

        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscspace_type", "sum", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscspace_degree", "1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscspace_variables", "2", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscspace_components", "2", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscspace_sum_spaces", "2", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_petscspace_sum_concatenate", "true", ierr)

        ! The first component's space: degree one in x, degree zero in z.
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_sumcomp_0_petscspace_variables", "2", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_sumcomp_0_petscspace_type", "tensor", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_sumcomp_0_petscspace_tensor_spaces", "2", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_sumcomp_0_petscspace_tensor_uniform", "false", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_sumcomp_0_tensorcomp_0_petscspace_degree", "1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_sumcomp_0_tensorcomp_1_petscspace_degree", "0", ierr)

        ! The second component's: degree zero in x, degree one in z.
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_sumcomp_1_petscspace_variables", "2", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_sumcomp_1_petscspace_type", "tensor", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_sumcomp_1_petscspace_tensor_spaces", "2", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_sumcomp_1_petscspace_tensor_uniform", "false", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_sumcomp_1_tensorcomp_0_petscspace_degree", "0", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_sumcomp_1_tensorcomp_1_petscspace_degree", "1", ierr)

        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscdualspace_type", "lagrange", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscdualspace_order", "1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-colflux_petscdualspace_form_degree", "-1", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_petscdualspace_lagrange_trimmed", "true", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                                 "-colflux_petscdualspace_lagrange_continuity", "true", ierr)
    end subroutine set_tensor_flux_options

    !> Remove every key only the tensor flux element reads.
    subroutine clear_tensor_flux_options(ierr)
        implicit none
        PetscErrorCode, intent(inout) :: ierr

        character(len=*), parameter :: KEYS(15) = [character(len=52) :: &
            "-colflux_petscspace_variables", &
            "-colflux_petscspace_sum_spaces", &
            "-colflux_petscspace_sum_concatenate", &
            "-colflux_sumcomp_0_petscspace_variables", &
            "-colflux_sumcomp_0_petscspace_type", &
            "-colflux_sumcomp_0_petscspace_tensor_spaces", &
            "-colflux_sumcomp_0_petscspace_tensor_uniform", &
            "-colflux_sumcomp_0_tensorcomp_0_petscspace_degree", &
            "-colflux_sumcomp_0_tensorcomp_1_petscspace_degree", &
            "-colflux_sumcomp_1_petscspace_variables", &
            "-colflux_sumcomp_1_petscspace_type", &
            "-colflux_sumcomp_1_petscspace_tensor_spaces", &
            "-colflux_sumcomp_1_petscspace_tensor_uniform", &
            "-colflux_sumcomp_1_tensorcomp_0_petscspace_degree", &
            "-colflux_sumcomp_1_tensorcomp_1_petscspace_degree"]
        integer(int32) :: k

        do k = 1, size(KEYS)
            call PetscOptionsClearValue(PETSC_NULL_OPTIONS, trim(KEYS(k)), ierr)
            if (ierr /= 0) return
        end do
    end subroutine clear_tensor_flux_options

    !> Register every residual and every active Jacobian block.
    !>
    !> All six cross blocks are registered. An omitted one would leave Newton
    !> with a coupling the residual carries and the matrix does not, which is
    !> exactly the failure a monolithic step is built to avoid.
    subroutine register_column_kernels(ds_handle, cerr)
        implicit none
        type(c_ptr), intent(in) :: ds_handle
        integer(c_int), intent(inout) :: cerr

        integer(c_int), parameter :: T = 0_c_int, Q = 1_c_int, P = 2_c_int

        call petscds_set_residual(ds_handle, T, c_funloc(f0_column_temperature), &
                                  c_funloc(f1_column_temperature), cerr)
        if (cerr == 0) call petscds_set_residual(ds_handle, Q, c_funloc(f0_column_flux), &
                                                 c_funloc(f1_column_flux), cerr)
        if (cerr == 0) call petscds_set_residual(ds_handle, P, c_funloc(f0_column_conservation), &
                                                 c_null_funptr, cerr)
        if (cerr /= 0) return

        call petscds_set_jacobian(ds_handle, T, T, &
                                  c_funloc(g0_column_temperature_temperature), c_null_funptr, &
                                  c_funloc(g2_column_temperature_temperature), &
                                  c_funloc(g3_column_temperature_temperature), cerr)
        if (cerr == 0) call petscds_set_jacobian(ds_handle, T, P, &
                                                 c_funloc(g0_column_temperature_pressure), c_null_funptr, &
                                                 c_funloc(g2_column_temperature_pressure), c_null_funptr, cerr)
        if (cerr == 0) call petscds_set_jacobian(ds_handle, Q, Q, &
                                                 c_funloc(g0_column_flux_flux), c_null_funptr, &
                                                 c_null_funptr, c_null_funptr, cerr)
        if (cerr == 0) call petscds_set_jacobian(ds_handle, Q, T, &
                                                 c_funloc(g0_column_flux_temperature), c_null_funptr, &
                                                 c_funloc(g2_column_flux_temperature), c_null_funptr, cerr)
        if (cerr == 0) call petscds_set_jacobian(ds_handle, Q, P, &
                                                 c_funloc(g0_column_flux_pressure), c_null_funptr, &
                                                 c_funloc(g2_column_flux_pressure), c_null_funptr, cerr)
        if (cerr == 0) call petscds_set_jacobian(ds_handle, P, Q, &
                                                 c_null_funptr, c_funloc(g1_column_pressure_flux), &
                                                 c_null_funptr, c_null_funptr, cerr)
        if (cerr == 0) call petscds_set_jacobian(ds_handle, P, P, &
                                                 c_funloc(g0_column_pressure_pressure), c_null_funptr, &
                                                 c_null_funptr, c_null_funptr, cerr)
        if (cerr == 0) call petscds_set_jacobian(ds_handle, P, T, &
                                                 c_funloc(g0_column_pressure_temperature), c_null_funptr, &
                                                 c_null_funptr, c_null_funptr, cerr)
    end subroutine register_column_kernels

    !> Register the primal layout [T, p]: the heat row as in the mixed form and
    !> the water row of `app_mixed_hydraulic_column_primal`.
    !>
    !> Every block is registered for the same reason as in the mixed layout.
    !> The storage blocks of the water row are the mixed column's own kernels:
    !> they read the pressure through `CC_PRESSURE_FIELD`. With `nodal` the
    !> flux kernels are left out, because `app_mixed_hydraulic_column_nodal`
    !> adds that term from the nodes.
    subroutine register_primal_kernels(ds_handle, nodal, cerr)
        implicit none
        type(c_ptr), intent(in) :: ds_handle
        !> Whether the flux is assembled from the nodes outside this DS.
        logical, intent(in) :: nodal
        integer(c_int), intent(inout) :: cerr

        integer(c_int), parameter :: T = 0_c_int
        integer(c_int), parameter :: P = int(COLUMN_PRIMAL_FIELD_PRESSURE, c_int)

        call petscds_set_residual(ds_handle, T, c_funloc(f0_column_temperature), &
                                  c_funloc(f1_column_temperature), cerr)
        if (cerr /= 0) return
        if (nodal) then
            call petscds_set_residual(ds_handle, P, c_funloc(f0_primal_conservation), c_null_funptr, cerr)
        else
            call petscds_set_residual(ds_handle, P, c_funloc(f0_primal_conservation), &
                                      c_funloc(f1_primal_conservation), cerr)
        end if
        if (cerr /= 0) return

        call petscds_set_jacobian(ds_handle, T, T, &
                                  c_funloc(g0_column_temperature_temperature), c_null_funptr, &
                                  c_funloc(g2_column_temperature_temperature), &
                                  c_funloc(g3_column_temperature_temperature), cerr)
        if (cerr == 0) call petscds_set_jacobian(ds_handle, T, P, &
                                                 c_funloc(g0_column_temperature_pressure), c_null_funptr, &
                                                 c_funloc(g2_column_temperature_pressure), c_null_funptr, cerr)
        if (cerr /= 0) return
        if (nodal) then
            call petscds_set_jacobian(ds_handle, P, P, c_funloc(g0_column_pressure_pressure), c_null_funptr, &
                                      c_null_funptr, c_null_funptr, cerr)
            if (cerr == 0) call petscds_set_jacobian(ds_handle, P, T, c_funloc(g0_column_pressure_temperature), &
                                                     c_null_funptr, c_null_funptr, c_null_funptr, cerr)
        else
            call petscds_set_jacobian(ds_handle, P, P, &
                                      c_funloc(g0_column_pressure_pressure), c_null_funptr, &
                                      c_funloc(g2_primal_pressure_pressure), &
                                      c_funloc(g3_primal_pressure_pressure), cerr)
            if (cerr == 0) call petscds_set_jacobian(ds_handle, P, T, &
                                                     c_funloc(g0_column_pressure_temperature), c_null_funptr, &
                                                     c_funloc(g2_primal_pressure_temperature), &
                                                     c_funloc(g3_primal_pressure_temperature), cerr)
        end if
    end subroutine register_primal_kernels

    !> Write the run constants into the discrete system.
    subroutine store_column_constants(self, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        self%constants(CC_DT) = real(self%dt, PETSC_REAL_KIND)
        self%constants(CC_FLUX_SCALE) = real(self%flux_scale, PETSC_REAL_KIND)
        self%constants(CC_PRESSURE_SCALE) = real(self%pressure_scale, PETSC_REAL_KIND)
        self%constants(CC_LENGTH_SCALE) = real(self%length_scale, PETSC_REAL_KIND)
        self%constants(CC_DENSITY) = real(self%density, PETSC_REAL_KIND)
        self%constants(CC_GRAVITY_X) = real(self%gravity(1), PETSC_REAL_KIND)
        self%constants(CC_GRAVITY_Z) = real(self%gravity(2), PETSC_REAL_KIND)
        self%constants(CC_TEMPERATURE_SCALE) = real(self%temperature_scale, PETSC_REAL_KIND)
        self%constants(CC_THERMAL_SCALE) = real(self%thermal_scale, PETSC_REAL_KIND)
        self%constants(CC_PRESSURE_FIELD) = real(self%pressure_field, PETSC_REAL_KIND)

        call PetscDSSetConstants(self%ds, int(COLUMN_CONSTANTS, PETSC_INT_KIND), self%constants, ierr)
    end subroutine store_column_constants

    !> The cell-wise coefficients and the accepted state, on a DM of their own.
    !>
    !> Every auxiliary element is built on the same reference cell as the
    !> unknowns and takes the temperature element's own quadrature. The element
    !> loop reads an auxiliary field at the main field's quadrature points and
    !> indexes the two tabulations together, so an auxiliary element built on
    !> the wrong reference cell or integrated on a different rule would not
    !> merely be inaccurate: it would carry a different number of degrees of
    !> freedom and a different number of quadrature points than the rows it is
    !> evaluated for, which `PetscFEIntegrateResidual_Basic` refuses outright.
    subroutine build_column_auxiliary(self, dimension, is_simplex, fe_temperature, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        PetscInt, intent(in) :: dimension
        !> Whether the mesh's cells are simplices.
        PetscBool, intent(in) :: is_simplex
        !> The temperature element, whose cell and face quadrature every
        !> auxiliary element copies; the unknowns already share it.
        PetscFE, intent(in) :: fe_temperature
        PetscErrorCode, intent(inout) :: ierr

        call DMClone(self%dm, self%aux_dm, ierr)
        if (ierr /= 0) return

        call add_column_aux_field(self%aux_dm, int(CA_INVERSE_MOBILITY - 1, PETSC_INT_KIND), "inverse_mobility", &
                                  dimension, int(MIXED_DIMENSION*MIXED_DIMENSION, PETSC_INT_KIND), &
                                  0_PETSC_INT_KIND, is_simplex, fe_temperature, ierr)
        if (ierr == 0) call add_column_aux_field(self%aux_dm, int(CA_MATERIAL - 1, PETSC_INT_KIND), "material", &
                                                 dimension, 1_PETSC_INT_KIND, 0_PETSC_INT_KIND, is_simplex, &
                                                 fe_temperature, ierr)
        if (ierr == 0) call add_column_aux_field(self%aux_dm, int(CA_SOURCE - 1, PETSC_INT_KIND), "source", &
                                                 dimension, 1_PETSC_INT_KIND, 0_PETSC_INT_KIND, is_simplex, &
                                                 fe_temperature, ierr)
        ! The accepted temperature carries the same nodal space as the unknown,
        ! so the storage rate differences the closure at the quadrature point
        ! rather than at a cell average of it.
        if (ierr == 0) call add_column_aux_field(self%aux_dm, int(CA_TEMPERATURE_OLD - 1, PETSC_INT_KIND), &
                                                 "temperature_old", dimension, 1_PETSC_INT_KIND, &
                                                 1_PETSC_INT_KIND, is_simplex, fe_temperature, ierr)
        ! The accepted pressure, likewise in the pressure's own space.
        if (ierr == 0) call add_column_aux_field(self%aux_dm, int(CA_PRESSURE_OLD - 1, PETSC_INT_KIND), &
                                                 "pressure_old", dimension, 1_PETSC_INT_KIND, &
                                                 merge(1_PETSC_INT_KIND, 0_PETSC_INT_KIND, self%primal), &
                                                 is_simplex, fe_temperature, ierr)
        ! Read by the heat-transfer boundary kernels only, out of the copy that
        ! boundary's own auxiliary vector holds. Zero in the default vector, so
        ! a boundary that was never registered contributes nothing.
        if (ierr == 0) call add_column_aux_field(self%aux_dm, int(CA_HEAT_COEFFICIENT - 1, PETSC_INT_KIND), &
                                                 "heat_coefficient", dimension, 1_PETSC_INT_KIND, &
                                                 0_PETSC_INT_KIND, is_simplex, fe_temperature, ierr)
        if (ierr == 0) call add_column_aux_field(self%aux_dm, int(CA_HEAT_AMBIENT - 1, PETSC_INT_KIND), &
                                                 "heat_ambient", dimension, 1_PETSC_INT_KIND, &
                                                 0_PETSC_INT_KIND, is_simplex, fe_temperature, ierr)
        if (ierr == 0) call add_column_aux_field(self%aux_dm, int(CA_HEAT_COEFFICIENT_FROZEN - 1, PETSC_INT_KIND), &
                                                 "heat_coefficient_frozen", dimension, 1_PETSC_INT_KIND, &
                                                 0_PETSC_INT_KIND, is_simplex, fe_temperature, ierr)
        if (ierr == 0) call add_column_aux_field(self%aux_dm, int(CA_HEAT_SWITCH_TEMPERATURE - 1, PETSC_INT_KIND), &
                                                 "heat_switch_temperature", dimension, 1_PETSC_INT_KIND, &
                                                 0_PETSC_INT_KIND, is_simplex, fe_temperature, ierr)
        if (ierr == 0) call add_column_aux_field(self%aux_dm, int(CA_HEAT_SWITCH_WIDTH - 1, PETSC_INT_KIND), &
                                                 "heat_switch_width", dimension, 1_PETSC_INT_KIND, &
                                                 0_PETSC_INT_KIND, is_simplex, fe_temperature, ierr)
        if (ierr /= 0) return

        call DMCreateDS(self%aux_dm, ierr)
        if (ierr == 0) call DMCreateLocalVector(self%aux_dm, self%auxiliary, ierr)
        if (ierr == 0) call VecSet(self%auxiliary, 0.0_PETSC_REAL_KIND, ierr)
        if (ierr == 0) call DMSetAuxiliaryVec(self%dm, PETSC_NULL_DMLABEL, 0_PETSC_INT_KIND, &
                                              0_PETSC_INT_KIND, self%auxiliary, ierr)
        if (ierr /= 0) return

        call store_column_materials(self, ierr)
        if (ierr /= 0) return

        call DMCreateLocalVector(self%aux_dm, self%accepted_auxiliary, ierr)
        if (ierr == 0) call VecCopy(self%auxiliary, self%accepted_auxiliary, ierr)
        if (ierr == 0) call DMCreateLocalVector(self%aux_dm, self%previous_auxiliary, ierr)
        if (ierr == 0) call VecCopy(self%auxiliary, self%previous_auxiliary, ierr)
    end subroutine build_column_auxiliary

    !> One auxiliary field of `num_components` Lagrange degrees of freedom.
    subroutine add_column_aux_field(dm, field, name, dimension, num_components, order, is_simplex, &
                                    fe_temperature, ierr)
        implicit none
        DM, intent(inout) :: dm
        PetscInt, intent(in) :: field
        character(len=*), intent(in) :: name
        PetscInt, intent(in) :: dimension
        PetscInt, intent(in) :: num_components
        !> Polynomial order; zero for a cell constant, one for a vertex field.
        PetscInt, intent(in) :: order
        !> Whether the reference cell is a simplex; the unknowns' own.
        PetscBool, intent(in) :: is_simplex
        !> The element whose cell and face quadrature this field takes.
        PetscFE, intent(in) :: fe_temperature
        PetscErrorCode, intent(inout) :: ierr

        PetscFE :: fe
        PetscObject :: discretisation
        PetscErrorCode :: release_ierr
        integer(c_int) :: cerr

        ! `PetscFECreateLagrange` reads no options, so the rule it builds here
        ! is QUADRATURE_ORDER whatever the unknowns were asked for. The copy
        ! that follows replaces it, cell rule and face rule both, before the
        ! field is attached and the auxiliary system is created from it.
        !
        ! An order-0 (cell-constant) Lagrange dual space, built here for every
        ! auxiliary field except temperature_old, hits the same PETSc-internal
        ! transient floating point exception documented at
        ! initialize_mixed_column's PetscFECreateDefault calls -- same
        ! PetscDualSpaceLagrangeCreateSimplexNodeMat, reached this time through
        ! PetscFECreateLagrange rather than PetscFECreateDefault. Same bracket,
        ! same reasoning: only this definition call runs with the trap off.
        call petscfptrappush_(PETSC_FP_TRAP_OFF_LOCAL, cerr)
        call PetscFECreateLagrange(PETSC_COMM_WORLD, dimension, num_components, is_simplex, &
                                   order, QUADRATURE_ORDER, fe, ierr)
        call petscfptrappop_(cerr)
        if (ierr /= 0) return
        call PetscFECopyQuadrature(fe_temperature, fe, ierr)

        if (ierr == 0) then
            discretisation = PetscObjectCast(fe)
            call PetscObjectSetName(discretisation, name, ierr)
            if (ierr == 0) call DMSetField(dm, field, PETSC_NULL_DMLABEL, discretisation, ierr)
        end if
        ! Released through its own code, so that a field that failed above is
        ! still reported as one rather than masked by a successful release.
        call PetscFEDestroy(fe, release_ierr)
        if (ierr == 0) ierr = release_ierr
    end subroutine add_column_aux_field

    !> Invert each cell's own conductivity into its auxiliary field.
    subroutine store_column_materials(self, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: cell, c_start, c_end, offset, num_dofs, material
        real(real64) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION), gravity_magnitude
        integer(int32) :: m, found, k
        logical :: ok

        gravity_magnitude = norm2(self%gravity)
        if (gravity_magnitude <= 0.0_real64) gravity_magnitude = 9.80665_real64

        call DMGetLocalSection(self%aux_dm, section, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(self%aux_dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return

        call VecGetArray(self%auxiliary, raw, ierr)
        if (ierr /= 0) return

        do cell = c_start, c_end - 1
            call PetscSectionGetFieldDof(section, cell, int(CA_INVERSE_MOBILITY - 1, PETSC_INT_KIND), num_dofs, ierr)
            if (ierr /= 0 .or. num_dofs <= 0) cycle

            call DMGetLabelValue(self%aux_dm, "Cell Sets", cell, material, ierr)
            if (ierr /= 0 .or. material < 0) then
                call VecRestoreArray(self%auxiliary, raw, ierr)
                error stop "A cell carries no material: the mesh has a region without a physical group"
            end if

            found = 0
            do m = 1, size(self%materials)
                if (self%materials(m)%id == int(material, int32)) found = m
            end do
            if (found == 0) then
                call VecRestoreArray(self%auxiliary, raw, ierr)
                error stop "No mixed column material is defined for this Cell Sets value"
            end if

            call invert_mobility(self%conductivity(:, :, found), self%density, gravity_magnitude, &
                                 inverse_mobility, ok)
            if (.not. ok) then
                call VecRestoreArray(self%auxiliary, raw, ierr)
                error stop "An impermeable region is a different model, not a small conductivity"
            end if

            call PetscSectionGetFieldOffset(section, cell, int(CA_INVERSE_MOBILITY - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr /= 0) cycle
            do k = 1, MIXED_DIMENSION*MIXED_DIMENSION
                raw(offset + k) = real(inverse_mobility(k), PETSC_REAL_KIND)
            end do

            call PetscSectionGetFieldOffset(section, cell, int(CA_MATERIAL - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr == 0) raw(offset + 1) = real(found, PETSC_REAL_KIND)
            call PetscSectionGetFieldOffset(section, cell, int(CA_SOURCE - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr == 0) raw(offset + 1) = 0.0_PETSC_REAL_KIND
        end do

        call VecRestoreArray(self%auxiliary, raw, ierr)
    end subroutine store_column_materials

    !> Register the closed water boundaries and the prescribed temperatures.
    subroutine add_column_boundaries(self, dm_handle, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        type(c_ptr), intent(in) :: dm_handle
        PetscErrorCode, intent(inout) :: ierr

        integer(int32) :: b
        integer(c_int) :: boundary_index, cerr
        character(len=40) :: name

        call reset_column_temperatures()
        ierr = 0

        ! Every boundary facet has to be constrained. An unconstrained RT facet
        ! is not a no-flow boundary: with no boundary term in the flux row the
        ! weak form's natural condition is a prescribed zero liquid pressure,
        ! which would let the column drain.
        ! The primal form has no flux unknown: a closed boundary is its natural
        ! condition and needs no constraint.
        do b = 1, merge(0, self%num_closed, self%primal)
            write (name, '(a,i0)') "column-closed-", self%closed_ids(b)
            call petscds_add_boundary_essential(dm_handle, trim(name), "Face Sets", &
                                                int(self%closed_ids(b), c_int), &
                                                int(COLUMN_FIELD_FLUX, c_int), &
                                                c_funloc(bc_zero_flux), c_null_ptr, &
                                                boundary_index, cerr)
            if (cerr /= 0) then
                ierr = int(cerr)
                return
            end if
        end do

        do b = 1, self%num_temperatures
            if (b > MAX_COLUMN_TEMPERATURES) error stop "More prescribed temperatures than the kernels can hold"
            call set_column_temperature_value(b, self%temperatures(b)%temperature, self%temperature_scale)
            write (name, '(a,i0)') "column-temperature-", self%temperatures(b)%boundary_id
            call petscds_add_boundary_essential(dm_handle, trim(name), trim(self%thermal_label), &
                                                int(self%temperatures(b)%boundary_id, c_int), &
                                                int(COLUMN_FIELD_TEMPERATURE, c_int), &
                                                c_funloc(bc_column_temperature), &
                                                column_temperature_context(b), &
                                                boundary_index, cerr)
            if (cerr /= 0) then
                ierr = int(cerr)
                return
            end if
        end do

        do b = 1, self%num_heat_transfers
            write (name, '(a,i0)') "column-heat-", self%heat_transfers(b)%boundary_id
            call petscds_add_boundary_natural(dm_handle, trim(name), trim(self%thermal_label), &
                                              int(self%heat_transfers(b)%boundary_id, c_int), &
                                              int(COLUMN_FIELD_TEMPERATURE, c_int), &
                                              c_funloc(bd_f0_column_heat), c_funloc(bd_g0_column_heat), &
                                              boundary_index, cerr)
            if (cerr /= 0) then
                ierr = int(cerr)
                return
            end if
            call attach_heat_transfer_auxiliary(self, b, ierr)
            if (ierr /= 0) return
        end do
    end subroutine add_column_boundaries

    !> Give one heat-transfer boundary the auxiliary vector its kernels read.
    !>
    !> A copy of the run's auxiliary vector with this boundary's coefficient and
    !> ambient temperature written into the two fields the boundary kernels use.
    !> The copy is never refreshed and does not need to be: the boundary kernels
    !> read those two fields and nothing else, so the accepted state the rest of
    !> the vector carries going stale changes no value they return.
    subroutine attach_heat_transfer_auxiliary(self, slot, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        integer(int32), intent(in) :: slot
        PetscErrorCode, intent(inout) :: ierr

        DMLabel :: label

        call DMGetLabel(self%dm, trim(self%thermal_label), label, ierr)
        if (ierr /= 0) return

        call DMCreateLocalVector(self%aux_dm, self%heat_auxiliary(slot), ierr)
        if (ierr /= 0) return
        call VecCopy(self%auxiliary, self%heat_auxiliary(slot), ierr)
        if (ierr == 0) call fill_column_aux_field(self%aux_dm, self%heat_auxiliary(slot), &
                                                  CA_HEAT_COEFFICIENT, &
                                                  self%heat_transfers(slot)%coefficient, ierr)
        if (ierr == 0) call fill_column_aux_field(self%aux_dm, self%heat_auxiliary(slot), &
                                                  CA_HEAT_AMBIENT, &
                                                  self%heat_transfers(slot)%ambient, ierr)
        associate (boundary => self%heat_transfers(slot))
            ! Without a switch the frozen coefficient is the warm one, which the
            ! kernel reads as a constant coefficient.
            if (ierr == 0) call fill_column_aux_field(self%aux_dm, self%heat_auxiliary(slot), &
                                                      CA_HEAT_COEFFICIENT_FROZEN, &
                                                      merge(boundary%frozen_coefficient, boundary%coefficient, &
                                                            boundary%frozen_coefficient >= 0.0_real64), ierr)
            if (ierr == 0) call fill_column_aux_field(self%aux_dm, self%heat_auxiliary(slot), &
                                                      CA_HEAT_SWITCH_TEMPERATURE, boundary%switch_temperature, ierr)
            if (ierr == 0) call fill_column_aux_field(self%aux_dm, self%heat_auxiliary(slot), &
                                                      CA_HEAT_SWITCH_WIDTH, boundary%switch_width, ierr)
        end associate
        if (ierr /= 0) return

        call DMSetAuxiliaryVec(self%dm, label, &
                               int(self%heat_transfers(slot)%boundary_id, PETSC_INT_KIND), &
                               0_PETSC_INT_KIND, self%heat_auxiliary(slot), ierr)
    end subroutine attach_heat_transfer_auxiliary

    !> Set every degree of freedom of one auxiliary field to the same value.
    subroutine fill_column_aux_field(dm, vec, field, value, ierr)
        implicit none
        DM, intent(in) :: dm
        Vec, intent(inout) :: vec
        integer(int32), intent(in) :: field
        real(real64), intent(in) :: value
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: point, p_start, p_end, offset, components, c

        call DMGetLocalSection(dm, section, ierr)
        if (ierr /= 0) return
        call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr /= 0) return
        call VecGetArray(vec, raw, ierr)
        if (ierr /= 0) return

        do point = p_start, p_end - 1
            call PetscSectionGetFieldDof(section, point, int(field - 1, PETSC_INT_KIND), components, ierr)
            if (ierr /= 0 .or. components <= 0) cycle
            call PetscSectionGetFieldOffset(section, point, int(field - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr /= 0) cycle
            do c = 1, components
                raw(offset + c) = real(value, PETSC_REAL_KIND)
            end do
        end do

        call VecRestoreArray(vec, raw, ierr)
    end subroutine fill_column_aux_field

    !> The constraint a closed water boundary imposes: no normal flux.
    function bc_zero_flux(dim, time, x, num_components, u, ctx) &
        bind(C, name="ftcms_column_bc_zero_flux") result(ierr)
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        integer(int32) :: i

        do i = 1, MIXED_DIMENSION
            u(i) = 0.0_PETSC_REAL_KIND
        end do
        ierr = 0
    end function bc_zero_flux

    !> Start the column from one temperature and one pore pressure.
    !>
    !> With `hydrostatic_height` the pore pressure is the hydrostatic profile
    !> through that height,
    !>
    !> \[ p_w(z)=p_w^{\mathrm{ref}}+\rho_l g_z\,(z-z^{\mathrm{ref}}), \]
    !>
    !> which is the state a closed column is at rest in. Without it every cell
    !> starts at the same pressure, which in a column under gravity is a state
    !> that will drain. \(O(N)\) in mesh points.
    subroutine set_column_initial_state(self, temperature, pressure, ierr, hydrostatic_height)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> Initial temperature, C, and pore pressure, Pa.
        real(real64), intent(in) :: temperature, pressure
        PetscErrorCode, intent(inout) :: ierr
        !> The height `pressure` is given at, m.
        real(real64), intent(in), optional :: hydrostatic_height

        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: point, p_start, p_end, offset, num_dofs, k
        PetscInt :: c_start, c_end
        PetscReal :: measure, centroid(3), normal(3)
        real(real64) :: cell_pressure
        Vec :: local

        call VecZeroEntries(self%solution, ierr)
        if (ierr /= 0) return

        ! The initial state is time zero, and its boundary values are those of
        ! time zero.
        self%time = 0.0_real64
        call apply_boundary_schedules(self, self%time, ierr)
        if (ierr /= 0) return

        call DMGetLocalVector(self%dm, local, ierr)
        if (ierr == 0) call VecZeroEntries(local, ierr)
        if (ierr == 0) call DMGetLocalSection(self%dm, section, ierr)
        if (ierr == 0) call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(self%dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return

        call VecGetArray(local, raw, ierr)
        if (ierr /= 0) return
        do point = p_start, p_end - 1
            call PetscSectionGetFieldDof(section, point, COLUMN_FIELD_TEMPERATURE, num_dofs, ierr)
            if (ierr == 0 .and. num_dofs > 0) then
                call PetscSectionGetFieldOffset(section, point, COLUMN_FIELD_TEMPERATURE, offset, ierr)
                if (ierr == 0) then
                    do k = 1, num_dofs
                        raw(offset + k) = real(temperature/self%temperature_scale, PETSC_REAL_KIND)
                    end do
                end if
            end if
            call PetscSectionGetFieldDof(section, point, self%pressure_field, num_dofs, ierr)
            if (ierr == 0 .and. num_dofs > 0) then
                cell_pressure = pressure
                ! The height of the point the degree of freedom sits on: a cell
                ! centroid for P0, a vertex for the primal space.
                if (present(hydrostatic_height) .and. &
                    (self%primal .or. (point >= c_start .and. point < c_end))) then
                    call DMPlexComputeCellGeometryFVM(self%dm, point, measure, centroid, normal, ierr)
                    if (ierr == 0) cell_pressure = pressure + self%density*self%gravity(2) &
                                                   *(real(centroid(2), real64) - hydrostatic_height)
                end if
                call PetscSectionGetFieldOffset(section, point, self%pressure_field, offset, ierr)
                if (ierr == 0) raw(offset + 1) = real(cell_pressure/self%pressure_scale, PETSC_REAL_KIND)
            end if
        end do
        call VecRestoreArray(local, raw, ierr)

        if (ierr == 0) call DMLocalToGlobalBegin(self%dm, local, INSERT_VALUES, self%solution, ierr)
        if (ierr == 0) call DMLocalToGlobalEnd(self%dm, local, INSERT_VALUES, self%solution, ierr)
        call DMRestoreLocalVector(self%dm, local, ierr)
        if (ierr == 0) call VecCopy(self%solution, self%old_solution, ierr)
        if (ierr == 0) call publish_accepted_state(self, ierr)
        if (ierr == 0) call VecCopy(self%auxiliary, self%accepted_auxiliary, ierr)
        ! Nothing has been stepped yet, so the state before the last step is the
        ! state itself and the audit measures a defect of zero.
        if (ierr == 0) call VecCopy(self%auxiliary, self%previous_auxiliary, ierr)
        self%accepted_dt = 0.0_real64
        if (allocated(self%transported)) self%transported = 0.0_real64
        self%system_balance = 0.0_real64
        self%element_balance = 0.0_real64
        if (ierr == 0 .and. allocated(self%cell_water)) call compute_cell_waters(self, self%cell_water, ierr)
    end subroutine set_column_initial_state

    !> Copy the accepted solution into the auxiliary history fields.
    !>
    !> The history is the accepted state itself, not a projection of a quantity
    !> computed from it, so a step re-evaluates the closure at exactly the state
    !> the previous step converged to.
    subroutine publish_accepted_state(self, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: solution_section, aux_section
        PetscScalar, pointer :: solution_raw(:), aux_raw(:)
        PetscInt :: point, p_start, p_end, offset, aux_offset, num_dofs, k
        Vec :: local

        call DMGetLocalVector(self%dm, local, ierr)
        if (ierr /= 0) return
        call DMGlobalToLocalBegin(self%dm, self%old_solution, INSERT_VALUES, local, ierr)
        if (ierr == 0) call DMGlobalToLocalEnd(self%dm, self%old_solution, INSERT_VALUES, local, ierr)
        if (ierr == 0) call DMPlexInsertBoundaryValues(self%dm, PETSC_TRUE, local, 0.0_PETSC_REAL_KIND, &
                                                       PETSC_NULL_VEC, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
        if (ierr == 0) call DMGetLocalSection(self%dm, solution_section, ierr)
        if (ierr == 0) call DMGetLocalSection(self%aux_dm, aux_section, ierr)
        if (ierr == 0) call PetscSectionGetChart(solution_section, p_start, p_end, ierr)
        if (ierr /= 0) then
            call DMRestoreLocalVector(self%dm, local, ierr)
            return
        end if

        call VecGetArrayRead(local, solution_raw, ierr)
        if (ierr == 0) call VecGetArray(self%auxiliary, aux_raw, ierr)
        if (ierr /= 0) then
            call DMRestoreLocalVector(self%dm, local, ierr)
            return
        end if

        do point = p_start, p_end - 1
            call PetscSectionGetFieldDof(solution_section, point, COLUMN_FIELD_TEMPERATURE, num_dofs, ierr)
            if (ierr == 0 .and. num_dofs > 0) then
                call PetscSectionGetFieldOffset(solution_section, point, COLUMN_FIELD_TEMPERATURE, offset, ierr)
                if (ierr == 0) call PetscSectionGetFieldOffset(aux_section, point, &
                                                               int(CA_TEMPERATURE_OLD - 1, PETSC_INT_KIND), &
                                                               aux_offset, ierr)
                if (ierr == 0) then
                    do k = 1, num_dofs
                        aux_raw(aux_offset + k) = real(self%temperature_scale*real(solution_raw(offset + k), real64), &
                                                       PETSC_REAL_KIND)
                    end do
                end if
            end if

            call PetscSectionGetFieldDof(solution_section, point, self%pressure_field, num_dofs, ierr)
            if (ierr == 0 .and. num_dofs > 0) then
                call PetscSectionGetFieldOffset(solution_section, point, self%pressure_field, offset, ierr)
                if (ierr == 0) call PetscSectionGetFieldOffset(aux_section, point, &
                                                               int(CA_PRESSURE_OLD - 1, PETSC_INT_KIND), &
                                                               aux_offset, ierr)
                if (ierr == 0) aux_raw(aux_offset + 1) = real(self%pressure_scale &
                                                              *real(solution_raw(offset + 1), real64), PETSC_REAL_KIND)
            end if
        end do

        call VecRestoreArrayRead(local, solution_raw, ierr)
        call VecRestoreArray(self%auxiliary, aux_raw, ierr)
        call DMRestoreLocalVector(self%dm, local, ierr)
    end subroutine publish_accepted_state

    !> Solve one backward-Euler step of length `dt` seconds.
    !>
    !> \[
    !>   \int_K\frac{W^{n+1}-W^{n}}{\Delta t}+\int_K\nabla\cdot\mathbf q_l^{n+1}=\int_K s,
    !>   \qquad
    !>   \int_\Omega\frac{U^{n+1}-U^{n}}{\Delta t}w+\int_\Omega\lambda\nabla T^{n+1}\cdot\nabla w=0.
    !> \]
    !>
    !> Assumes the accepted history is the state of the previous successful
    !> call. Guarantees that the history advances only when the solve converged
    !> to a finite state; on every other path the trial solution, the accepted
    !> solution and the accepted auxiliary state are restored exactly, so the
    !> same state can be retried at a shorter length. Cost is one nonlinear
    !> solve.
    !>
    !> Fails by returning `accepted = .false.` with `ierr` set: the PETSc error
    !> when the solve failed, zero for a diverged reason, `COLUMN_ERR_INVALID_STEP`
    !> for a step that is not positive and finite, `COLUMN_ERR_NONFINITE` for a
    !> converged but nonfinite state, and `COLUMN_ERR_CLOSURE` when the frozen
    !> closure could not answer at some iterate, which takes precedence.
    subroutine solve_column_step(self, dt, accepted, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> The step length, s; positive.
        real(real64), intent(in) :: dt
        !> True when the step converged and the history advanced.
        logical, intent(inout) :: accepted
        PetscErrorCode, intent(inout) :: ierr

        SNESConvergedReason :: reason
        PetscErrorCode :: restore_ierr
        real(real64), allocatable :: increment(:, :), new_water(:)

        accepted = .false.
        ierr = 0
        restore_ierr = 0
        if (.not. self%ready) then
            ierr = 1
            return
        end if
        if (.not. (dt > 0.0_real64) .or. .not. ieee_is_finite(dt)) then
            ierr = COLUMN_ERR_INVALID_STEP
            return
        end if

        self%dt = dt
        call store_column_constants(self, ierr)
        if (ierr /= 0) then
            call restore_column_state(self, restore_ierr)
            return
        end if

        call clear_column_fault()

        ! The step is implicit in its boundary as well: solve it at the
        ! prescribed temperatures of the time it ends at.
        call apply_boundary_schedules(self, self%time + dt, ierr)
        if (ierr /= 0) then
            call restore_column_state(self, restore_ierr)
            return
        end if

        call SNESSetSolution(self%snes, self%solution, ierr)
        if (ierr == 0) call SNESSolve(self%snes, PETSC_NULL_VEC, self%solution, ierr)
        if (ierr /= 0) then
            call restore_column_state(self, restore_ierr)
            call report_column_fault(ierr)
            return
        end if

        call SNESGetConvergedReason(self%snes, reason, ierr)
        if (ierr /= 0 .or. int(reason%v, int32) <= 0) then
            call restore_column_state(self, restore_ierr)
            call discard_failed_factorization(self, reason, restore_ierr)
            call report_column_fault(ierr)
            return
        end if

        if (.not. column_solution_is_finite(self, ierr)) then
            call restore_column_state(self, restore_ierr)
            if (ierr == 0) ierr = COLUMN_ERR_NONFINITE
            return
        end if

        ! Measured before anything is committed, so that a failure here can
        ! still put the solver back on the accepted state.
        allocate (increment, mold=self%transported)
        call compute_step_transport(self, dt, increment, ierr)
        allocate (new_water, mold=self%cell_water)
        if (ierr == 0) call compute_cell_waters(self, new_water, ierr)
        if (ierr /= 0) then
            call restore_column_state(self, restore_ierr)
            return
        end if

        call VecCopy(self%solution, self%old_solution, ierr)
        ! Before the history is republished, and only on a step that is being
        ! accepted: this is the state the step just taken started from.
        if (ierr == 0) call VecCopy(self%accepted_auxiliary, self%previous_auxiliary, ierr)
        if (ierr == 0) call publish_accepted_state(self, ierr)
        if (ierr == 0) call VecCopy(self%auxiliary, self%accepted_auxiliary, ierr)
        if (ierr /= 0) then
            call restore_column_state(self, restore_ierr)
            return
        end if
        self%transported = self%transported + increment
        call measure_step_balance(new_water - self%cell_water, sum(increment, dim=1), new_water, &
                                  self%system_balance, self%element_balance)
        self%cell_water = new_water
        self%accepted_dt = dt
        self%time = self%time + dt
        accepted = .true.
    end subroutine solve_column_step

    !> The water each cell's own flux carries through each facet in one step, m2.
    !>
    !> \(\Delta t\,Q_{Kf}\) at the converged state of the step, which is the
    !> flux backward Euler integrates over it: the RT0 facet moments for the
    !> mixed form, the cell's own Darcy flux on the facet for the primal forms
    !> (`primal_face_outflows`). \(O(N)\) in owned cells. Fails by returning
    !> the first PETSc error or `COLUMN_ERR_CLOSURE`, with `increment` partly
    !> written.
    subroutine compute_step_transport(self, dt, increment, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> The step length, s.
        real(real64), intent(in) :: dt
        !> \(\Delta t\,Q_{Kf}\), shaped as `transported`; overwritten.
        real(real64), intent(inout) :: increment(:, :)
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: aux_section
        PetscScalar, pointer :: aux_raw(:)
        PetscInt :: cell, c_start, c_end, facets(MAX_CELL_FACETS)
        real(real64) :: moments(MAX_CELL_FACETS)
        integer(int32) :: num_facets, column
        PetscErrorCode :: release_ierr

        increment = 0.0_real64
        call DMGlobalToLocalBegin(self%dm, self%solution, INSERT_VALUES, self%local_solution, ierr)
        if (ierr == 0) call DMGlobalToLocalEnd(self%dm, self%solution, INSERT_VALUES, self%local_solution, ierr)
        if (ierr == 0) call DMPlexInsertBoundaryValues(self%dm, PETSC_TRUE, self%local_solution, &
                                                       0.0_PETSC_REAL_KIND, PETSC_NULL_VEC, &
                                                       PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
        if (ierr == 0) call DMGetLocalSection(self%aux_dm, aux_section, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(self%dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return
        call VecGetArrayRead(self%accepted_auxiliary, aux_raw, ierr)
        if (ierr /= 0) return

        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(self%dm, cell)) cycle
            if (self%primal) then
                call primal_face_outflows(self, self%local_solution, cell, aux_raw, aux_section, facets, moments, &
                                          num_facets, ierr)
            else
                call reconstruct_cell_facet_flux(self%dm, self%local_solution, self%flux_scale, cell, facets, &
                                                 moments, num_facets, ierr, field=COLUMN_FIELD_FLUX)
            end if
            if (ierr /= 0) exit
            column = int(cell - c_start, int32) + 1
            increment(1:num_facets, column) = dt*moments(1:num_facets)
        end do

        call VecRestoreArrayRead(self%accepted_auxiliary, aux_raw, release_ierr)
        if (ierr == 0) ierr = release_ierr
    end subroutine compute_step_transport

    !> Every owned cell's water \(\int_K W\) at the current solution, m2.
    !>
    !> On the assembly's quadrature (`cell_water_content`), so the change of it
    !> over a step is the storage term the water row carries. \(O(N)\) in
    !> owned cells. Fails by returning the first PETSc error or
    !> `COLUMN_ERR_CLOSURE`.
    subroutine compute_cell_waters(self, water, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> \(\int_K W\) per cell, shaped as `cell_water`; overwritten, zero on cells
        !> another rank owns.
        real(real64), intent(inout) :: water(:)
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: aux_section
        PetscScalar, pointer :: aux_raw(:)
        PetscInt :: cell, c_start, c_end, offset
        real(real64) :: pressure, volume
        integer(int32) :: slot
        PetscErrorCode :: release_ierr

        water = 0.0_real64
        call DMGlobalToLocalBegin(self%dm, self%solution, INSERT_VALUES, self%local_solution, ierr)
        if (ierr == 0) call DMGlobalToLocalEnd(self%dm, self%solution, INSERT_VALUES, self%local_solution, ierr)
        if (ierr == 0) call DMPlexInsertBoundaryValues(self%dm, PETSC_TRUE, self%local_solution, &
                                                       0.0_PETSC_REAL_KIND, PETSC_NULL_VEC, &
                                                       PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
        if (ierr == 0) call DMGetLocalSection(self%aux_dm, aux_section, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(self%dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return
        call VecGetArrayRead(self%accepted_auxiliary, aux_raw, ierr)
        if (ierr /= 0) return

        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(self%dm, cell)) cycle
            call PetscSectionGetFieldOffset(aux_section, cell, int(CA_MATERIAL - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr /= 0) exit
            slot = nint(real(aux_raw(offset + 1), real64), int32)
            call cell_volume(self%dm, cell, volume, ierr)
            if (ierr == 0) call reconstruct_cell_pressure(self%dm, self%local_solution, self%pressure_scale, cell, &
                                                          pressure, ierr, field=self%pressure_field)
            if (ierr == 0) call cell_water_content(self, self%local_solution, cell, slot, pressure, volume, &
                                                   water(int(cell - c_start, int32) + 1), ierr)
            if (ierr /= 0) exit
        end do

        call VecRestoreArrayRead(self%accepted_auxiliary, aux_raw, release_ierr)
        if (ierr == 0) ierr = release_ierr
    end subroutine compute_cell_waters

    !> System and element mass balance of one step.
    !>
    !> \[
    !>   e_{\mathrm{sys}}=\frac{|\sum_K\Delta W_K|}{\sum_K W_K},\qquad
    !>   e_{\mathrm{el}}=\max_K\frac{|\Delta W_K+\Delta t\sum_f Q_{Kf}|}{W_K},
    !> \]
    !>
    !> for a closed column, with \(\Delta W_K\) the change of the cell's water
    !> over the step and \(Q_{Kf}\) the method's own facet flux. Both are
    !> relative to water a cell or the column actually holds, so a cell that
    !> barely changed is not divided by its change. \(O(N)\). Cannot fail.
    pure subroutine measure_step_balance(change, outflow, water, system, element)
        implicit none
        !> \(\Delta W_K\), m2, per cell.
        real(real64), intent(in) :: change(:)
        !> \(\Delta t\sum_f Q_{Kf}\), m2, per cell.
        real(real64), intent(in) :: outflow(:)
        !> \(W_K\) at the end of the step, m2, per cell; zero for cells not owned.
        real(real64), intent(in) :: water(:)
        !> \(e_{\mathrm{sys}}\) and \(e_{\mathrm{el}}\); overwritten.
        real(real64), intent(inout) :: system, element

        integer(int32) :: k

        system = 0.0_real64
        element = 0.0_real64
        if (sum(water) > 0.0_real64) system = abs(sum(change))/sum(water)
        do k = 1, size(water)
            if (water(k) > 0.0_real64) element = max(element, abs(change(k) + outflow(k))/water(k))
        end do
    end subroutine measure_step_balance

    !> Throw away a factorization that failed, so the retry builds a new one.
    !>
    !> A direct solver that fails to factorize does not come back to a neutral
    !> state. MUMPS keeps its own error flag, and PETSc's interface to it skips
    !> every later numeric factorization while that flag is set -- returning
    !> success, because the skip is how it avoids refactorizing a matrix it has
    !> already refused. PETSc then clears its own record of the error, finds no
    !> failure to report, and calls the triangular solve on a factorization that
    !> was never completed. That solve is the hard error, and by the time it is
    !> raised the step that caused it has already been rejected and forgotten;
    !> the run dies on the next shorter step instead of taking it.
    !>
    !> Resetting the Krylov solver destroys the factor matrix and with it the
    !> stale flag, which is what makes the shorter step an honest fresh attempt.
    !> It keeps the solver's configured types and options, so the retry solves
    !> the same way the failed attempt was asked to. The symbolic factorization
    !> is rebuilt as well; that is the price of the reset, and it is paid only
    !> on a step the solver already refused.
    !>
    !> Errors are reported through `ierr` and left to the caller: this runs on a
    !> path that is already failing, and it must not replace the reason the step
    !> was rejected with a reason of its own.
    subroutine discard_failed_factorization(self, reason, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> Why the nonlinear solve stopped.
        SNESConvergedReason, intent(in) :: reason
        !> Nonzero if the reset itself failed; never cleared here.
        PetscErrorCode, intent(inout) :: ierr

        KSP :: ksp
        PetscErrorCode :: local_ierr

        ! Only a linear failure can have left a factorization behind. Every
        ! other rejection -- a diverged Newton, a nonfinite residual -- leaves
        ! a usable factor that the next step is entitled to reuse the symbolic
        ! part of.
        if (int(reason%v, int32) /= int(SNES_DIVERGED_LINEAR_SOLVE%v, int32)) return

        call SNESGetKSP(self%snes, ksp, local_ierr)
        if (local_ierr == 0) call KSPReset(ksp, local_ierr)
        if (local_ierr /= 0 .and. ierr == 0) ierr = local_ierr
    end subroutine discard_failed_factorization

    !> Install the adaptive step controller and arm its first step.
    !>
    !> \(O(1)\). Fails by `error stop` on step bounds that are not positive and
    !> ordered, a growth factor below one, a reduction outside \((0,1]\), or an
    !> iteration window whose ends are the wrong way round.
    subroutine set_column_step_control(self, control)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        type(type_mixed_column_step_control), intent(in) :: control

        if (.not. (control%initial > 0.0_real64) .or. .not. (control%minimum > 0.0_real64) &
            .or. .not. (control%maximum >= control%minimum)) then
            error stop "The mixed column step control needs positive ordered step bounds"
        end if
        if (.not. (control%grow >= 1.0_real64) .or. .not. (control%shrink > 0.0_real64) &
            .or. .not. (control%shrink <= 1.0_real64) .or. .not. (control%retry > 0.0_real64) &
            .or. .not. (control%retry < 1.0_real64)) then
            error stop "The mixed column step control needs a growth above one and reductions below it"
        end if
        if (control%iterations_low > control%iterations_high) then
            error stop "The mixed column step control needs an ordered iteration window"
        end if
        self%step_control = control
        self%next_step = min(max(control%initial, control%minimum), control%maximum)
    end subroutine set_column_step_control

    !> Record one row per step attempt of every later advance in `path`.
    !>
    !> A run reports the steps it took only in aggregate -- how many were
    !> accepted, the shortest, the iterations they cost -- which cannot say
    !> whether the controller ever reached its ceiling or what the steps in
    !> between were. The row is written as the attempt ends, so a run killed by
    !> a wall clock cap leaves the history of everything it did get through.
    !>
    !> Written on rank 0 only, one buffered line per attempt, so the cost is
    !> \(O(1)\) per attempt and \(O(1)\) memory: nothing is accumulated in
    !> store. Diagnostic only -- what is recorded is what the controller and the
    !> solver already decided, and no step is chosen differently because a
    !> history is open.
    !>
    !> Fails with `COLUMN_ERR_NO_HISTORY_FILE` when the file cannot be opened,
    !> leaving no history open and the run otherwise untouched.
    subroutine set_column_step_history(self, path, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> Where the rows are written; an existing file is replaced, and an
        !> empty name closes an open history instead of opening one.
        character(len=*), intent(in) :: path
        PetscErrorCode, intent(inout) :: ierr

        integer(int32) :: status, mpi_err
        PetscMPIInt :: rank

        ierr = 0
        call close_column_step_history(self)
        self%history_path = path
        self%history_rows = 0
        if (len_trim(path) == 0) return

        ! Only rank 0 writes, so only rank 0 can fail to open; the failure is
        ! broadcast or the other ranks would go on believing a history exists.
        status = 0
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, mpi_err)
        if (rank == 0) call open_step_history_file(self, status)
        call MPI_Bcast(status, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, mpi_err)
        if (status /= 0) then
            self%history_path = ""
            ierr = status
        end if
    end subroutine set_column_step_history

    !> Open the history on rank 0 and write the header row.
    subroutine open_step_history_file(self, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        integer(int32), intent(inout) :: ierr

        integer(int32) :: io_status

        ierr = 0
        open (newunit=self%history_unit, file=trim(self%history_path), status="replace", &
              action="write", iostat=io_status)
        if (io_status /= 0) then
            self%history_unit = -1
            ierr = COLUMN_ERR_NO_HISTORY_FILE
            return
        end if
        write (self%history_unit, '(a)', iostat=io_status) &
            "attempt,accepted,dt_s,next_step_s,iterations,wall_s,system_balance,element_balance"
        flush (self%history_unit)
    end subroutine open_step_history_file

    !> Close an open step history, if there is one.
    subroutine close_column_step_history(self)
        implicit none
        class(type_mixed_column), intent(inout) :: self

        ! `newunit` hands back a negative number, so -1 and not a sign is what
        ! says whether there is a unit to close.
        if (self%history_unit /= -1) close (self%history_unit)
        self%history_unit = -1
    end subroutine close_column_step_history

    !> Record one attempt, when a history is open.
    !>
    !> The wall clock is recorded with the attempt rather than differenced here
    !> so that the reader can measure the cost of any span of the run, including
    !> one that spans the horizons an advance is cut at.
    subroutine write_column_step_row(self, dt, accepted, iterations)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> The length the attempt was made at, s.
        real(real64), intent(in) :: dt
        !> Whether the solver converged and the accepted state moved.
        logical, intent(in) :: accepted
        !> Newton iterations the attempt cost; zero when it was refused.
        integer(int32), intent(in) :: iterations

        real(real64) :: now
        integer(int32) :: flag
        PetscErrorCode :: local_ierr

        if (self%history_unit == -1) return
        call PetscTime(now, local_ierr)
        if (local_ierr /= 0) now = 0.0_real64
        flag = 0
        if (accepted) flag = 1
        self%history_rows = self%history_rows + 1
        ! The balance belongs to the accepted state; a refused attempt has none.
        write (self%history_unit, '(i0,a,i0,2(a,es16.9),a,i0,a,es20.9,2(a,es16.9))') &
            self%history_rows, ",", flag, ",", dt, ",", self%next_step, ",", iterations, ",", now, &
            ",", merge(self%system_balance, -1.0_real64, accepted), &
            ",", merge(self%element_balance, -1.0_real64, accepted)
        flush (self%history_unit)
    end subroutine write_column_step_row

    !> Carry the accepted state forward by `duration` seconds.
    !>
    !> Steps are chosen by the controller `set_step_control` installed, clipped
    !> so the last one lands exactly on the horizon. A step the solver refuses
    !> costs nothing but time: `solve_step` restores the accepted history before
    !> it returns, so the retry starts from the state the refused step began at,
    !> and the only thing the refusal changes is the length of the next attempt.
    !>
    !> `progress` reports what happened, not what was planned: the time actually
    !> reached, the steps accepted and refused, and the Newton iterations they
    !> cost. A caller that needs a horizon met exactly has to check `reached`.
    !>
    !> No accuracy is guaranteed: the controller reads the iteration count the
    !> solver reports and nothing else, so meeting a horizon says every step
    !> converged and not that the step was short enough to be accurate.
    !>
    !> One nonlinear solve per attempt, so \(O(\text{attempts})\) solves and no
    !> storage beyond the histories `solve_step` already keeps. Fails with
    !> `COLUMN_ERR_INVALID_SPAN` on a duration that is not positive and finite,
    !> and with `COLUMN_ERR_STEP_FLOOR` when a step has been shortened to the
    !> controller's minimum and refused there, leaving the accepted state at the
    !> last step that did converge and `reached` saying where that is.
    subroutine advance_column(self, duration, progress, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> The span to cover, s; positive.
        real(real64), intent(in) :: duration
        !> What the advance actually did.
        type(type_mixed_column_progress), intent(inout) :: progress
        PetscErrorCode, intent(inout) :: ierr

        real(real64) :: remaining, dt
        integer(int32) :: iterations
        logical :: accepted

        progress = type_mixed_column_progress()
        ierr = 0
        if (.not. self%ready) then
            ierr = 1
            return
        end if
        if (.not. (duration > 0.0_real64) .or. .not. ieee_is_finite(duration)) then
            ierr = COLUMN_ERR_INVALID_SPAN
            return
        end if
        if (.not. (self%next_step > 0.0_real64)) then
            call self%set_step_control(self%step_control)
        end if

        remaining = duration
        progress%next_step = self%next_step
        progress%shortest_step = huge(1.0_real64)

        do while (remaining > 0.0_real64)
            ! Never overshoot the horizon, and never let the clipped remainder
            ! become the step the controller carries on with: a horizon that
            ! happens to be short would otherwise halve the run's step for good.
            dt = min(self%next_step, remaining)
            accepted = .false.
            call self%solve_step(dt, accepted, ierr)
            if (accepted) then
                iterations = column_iteration_count(self)
                remaining = remaining - dt
                progress%reached = progress%reached + dt
                progress%last_step = dt
                progress%shortest_step = min(progress%shortest_step, dt)
                progress%accepted_steps = progress%accepted_steps + 1
                progress%total_iterations = progress%total_iterations + iterations
                progress%max_iterations = max(progress%max_iterations, iterations)
                self%next_step = next_column_step(self, self%next_step, iterations)
            else
                iterations = 0
                progress%rejected_steps = progress%rejected_steps + 1
                if (dt <= self%step_control%minimum) then
                    call write_column_step_row(self, dt, .false., iterations)
                    if (ierr == 0) ierr = COLUMN_ERR_STEP_FLOOR
                    progress%next_step = self%next_step
                    if (progress%accepted_steps == 0) progress%shortest_step = 0.0_real64
                    return
                end if
                self%next_step = max(dt*self%step_control%retry, self%step_control%minimum)
                ierr = 0
            end if
            call write_column_step_row(self, dt, accepted, iterations)
            progress%next_step = self%next_step
            call report_advance_progress(progress, duration, self%next_step, .false.)
        end do

        if (progress%accepted_steps == 0) progress%shortest_step = 0.0_real64
    end subroutine advance_column

    !> Say how far a long advance has actually got, at most every few minutes.
    !>
    !> An advance over a 50 h horizon on a 12800 cell mesh runs for hours and,
    !> without this, says nothing at all until it either finishes or is killed:
    !> the simulated time reached is the one number that distinguishes a run
    !> that is slow from a run that has stopped making progress, and the step
    !> the controller has fallen back to is what says which of the two it is.
    !>
    !> Wall clock gated, not step gated, so the cost is bounded however short
    !> the steps have become, and so that nothing is printed at all by a run
    !> shorter than the interval -- which is every run in the test suite, whose
    !> output is therefore unchanged. Reporting only, on rank zero only; it
    !> reads the progress record and changes nothing in it.
    subroutine report_advance_progress(progress, duration, next_step, force)
        implicit none
        !> What the advance has done so far.
        type(type_mixed_column_progress), intent(in) :: progress
        !> The span the advance was asked to cover, s.
        real(real64), intent(in) :: duration
        !> The length the controller would try next, s.
        real(real64), intent(in) :: next_step
        !> Report regardless of how long it is since the last one.
        logical, intent(in) :: force

        real(real64), save :: last_report = -1.0_real64
        real(real64) :: now
        PetscMPIInt :: rank
        PetscErrorCode :: local_ierr

        call PetscTime(now, local_ierr)
        if (local_ierr /= 0) return
        if (last_report < 0.0_real64) last_report = now
        if (.not. force .and. (now - last_report) < ADVANCE_REPORT_SECONDS) return
        last_report = now

        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, local_ierr)
        if (local_ierr /= 0 .or. rank /= 0) return
        write (*, '(a,f9.4,a,f9.4,a,i0,a,i0,a,es11.4,a,i0)') &
            "  progress  reached ", progress%reached/3600.0_real64, &
            " of ", duration/3600.0_real64, &
            " h  accepted ", progress%accepted_steps, &
            "  rejected ", progress%rejected_steps, &
            "  next step ", next_step, &
            "  its total ", progress%total_iterations
        flush (output_unit)
    end subroutine report_advance_progress

    !> The length the controller would try next after a step of `iterations`.
    pure function next_column_step(self, dt, iterations) result(next)
        implicit none
        class(type_mixed_column), intent(in) :: self
        real(real64), intent(in) :: dt
        integer(int32), intent(in) :: iterations
        real(real64) :: next

        next = dt
        if (iterations <= self%step_control%iterations_low) then
            next = dt*self%step_control%grow
        else if (iterations >= self%step_control%iterations_high) then
            next = dt*self%step_control%shrink
        end if
        next = min(max(next, self%step_control%minimum), self%step_control%maximum)
    end function next_column_step

    !> How many Newton iterations the last solve took, or zero if it cannot say.
    function column_iteration_count(self) result(iterations)
        implicit none
        class(type_mixed_column), intent(in) :: self
        integer(int32) :: iterations

        PetscInt :: taken
        PetscErrorCode :: local_ierr

        taken = 0_PETSC_INT_KIND
        call SNESGetIterationNumber(self%snes, taken, local_ierr)
        if (local_ierr /= 0) then
            iterations = 0
            return
        end if
        iterations = int(taken, int32)
    end function column_iteration_count

    !> Turn a recorded closure fault into this step's reason for failing.
    subroutine report_column_fault(ierr)
        implicit none
        PetscErrorCode, intent(inout) :: ierr

        integer(int32) :: slot
        real(real64) :: temperature, pressure
        PetscMPIInt :: rank
        PetscErrorCode :: local_ierr

        if (.not. column_fault_recorded()) return
        slot = 0
        temperature = 0.0_real64
        pressure = 0.0_real64
        call column_fault_report(slot, temperature, pressure)
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, local_ierr)
        if (rank == 0) then
            write (*, '(a,i0,a,es12.5,a,es12.5)') &
                "The frozen closure has no state for material slot ", slot, &
                " at T = ", temperature, " C, p_w = ", pressure
        end if
        ierr = COLUMN_ERR_CLOSURE
    end subroutine report_column_fault

    !> Put the solver back on the last accepted state.
    subroutine restore_column_state(self, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        call VecCopy(self%old_solution, self%solution, ierr)
        if (ierr == 0) call VecCopy(self%accepted_auxiliary, self%auxiliary, ierr)
    end subroutine restore_column_state

    !> Whether every entry of the current solution is a finite number.
    logical function column_solution_is_finite(self, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        PetscReal :: norm

        column_solution_is_finite = .false.
        call VecNorm(self%solution, NORM_INFINITY, norm, ierr)
        if (ierr /= 0) return
        column_solution_is_finite = ieee_is_finite(real(norm, real64))
    end function column_solution_is_finite

    !> The water balance, phase state and Jacobian margin of every owned cell.
    !>
    !> \[
    !>   d_K=\int_K(W^{n+1}-W^{n})+\Delta t\sum_{f\subset\partial K}Q_{Kf}-\Delta t\int_K s.
    !> \]
    !>
    !> The facet fluxes are reconstructed from the degrees of freedom rather
    !> than read out of the assembled residual, and \(W\) is re-evaluated
    !> through the same closure at the cell's own pressure and temperature, so
    !> this is an independent measurement of the step the solver took. The
    !> pressure column of each cell is measured as well: a column that has
    !> collapsed makes the matrix singular, and a run has to be able to say so
    !> rather than infer it from a solver return code. Cost is \(O(N)\) in owned
    !> cells.
    !>
    !> Fails by returning the first PETSc error; every output is initialised on
    !> every path. Assumes the solution vector holds the state the balance is
    !> claimed for and the auxiliary vector its accepted history, which is true
    !> between `solve_step` calls.
    !>
    !> The primal form has no flux unknown. Its facet fluxes are each cell's own
    !> Darcy flux, from that cell's continuous solution, integrated over the
    !> facet (`primal_face_outflows`), and the same defect is measured with
    !> them: how far the Galerkin flux is from balancing a cell's storage. The
    !> facet audit is then the disagreement between the two cells of a facet,
    !> summed on this rank only, so it is exact on one rank and misses facets
    !> on a partition boundary on several.
    subroutine compute_column_report(self, report, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> The measured quantities, reduced over the communicator.
        type(type_mixed_column_report), intent(inout) :: report
        PetscErrorCode, intent(inout) :: ierr

        type(type_mixed_frozen_response) :: now, accepted
        PetscSection :: aux_section
        PetscScalar, pointer :: aux_raw(:)
        Vec :: local
        PetscInt :: cell, c_start, c_end, facets(MAX_CELL_FACETS), offset, support_size
        real(real64) :: moments(MAX_CELL_FACETS)
        real(real64) :: volume, pressure, pressure_old, temperature, temperature_old
        real(real64) :: defect, outflow, margin, sensitivity, capacity
        real(real64) :: water_mean, old_mean, ice_mean, source
        real(real64) :: flux_values(MIXED_DIMENSION)
        real(real64), allocatable :: face_sum(:)
        integer(int32), allocatable :: face_count(:)
        PetscInt :: f_start, f_end
        integer(int32) :: num_facets, f, slot
        logical :: ok

        report = type_mixed_column_report()
        report%step_length = self%accepted_dt
        report%weakest_column = huge(1.0_real64)
        report%coldest = huge(1.0_real64)
        report%warmest = -huge(1.0_real64)
        report%lowest_liquid_pressure = huge(1.0_real64)

        call DMGlobalToLocalBegin(self%dm, self%solution, INSERT_VALUES, self%local_solution, ierr)
        if (ierr == 0) call DMGlobalToLocalEnd(self%dm, self%solution, INSERT_VALUES, self%local_solution, ierr)
        if (ierr == 0) call DMPlexInsertBoundaryValues(self%dm, PETSC_TRUE, self%local_solution, &
                                                       0.0_PETSC_REAL_KIND, PETSC_NULL_VEC, &
                                                       PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
        if (ierr /= 0) return
        local = self%local_solution

        call DMGetLocalSection(self%aux_dm, aux_section, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(self%dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return
        ! The primal facet fluxes of the two cells of a facet, summed per facet.
        if (self%primal) then
            call DMPlexGetHeightStratum(self%dm, 1_PETSC_INT_KIND, f_start, f_end, ierr)
            if (ierr /= 0) return
            allocate (face_sum(f_start:f_end - 1), face_count(f_start:f_end - 1))
            face_sum = 0.0_real64
            face_count = 0
        end if
        call VecGetArrayRead(self%previous_auxiliary, aux_raw, ierr)
        if (ierr /= 0) return

        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(self%dm, cell)) cycle

            call cell_volume(self%dm, cell, volume, ierr)
            if (ierr /= 0) exit
            call reconstruct_cell_pressure(self%dm, local, self%pressure_scale, cell, pressure, ierr, &
                                           field=self%pressure_field)
            if (ierr /= 0) exit
            call cell_temperature(self, local, cell, temperature, ierr)
            if (ierr /= 0) exit
            if (self%primal) then
                call primal_cell_flux(self, local, cell, aux_raw, aux_section, flux_values, ierr)
            else
                call cell_flux(self, local, cell, flux_values, ierr)
            end if
            if (ierr /= 0) exit

            if (self%primal) then
                ! A vertex field: its value at the centroid.
                call reconstruct_cell_pressure(self%aux_dm, self%previous_auxiliary, 1.0_real64, cell, &
                                               pressure_old, ierr, field=int(CA_PRESSURE_OLD - 1, PETSC_INT_KIND))
                if (ierr /= 0) exit
            else
                call PetscSectionGetFieldOffset(aux_section, cell, int(CA_PRESSURE_OLD - 1, PETSC_INT_KIND), &
                                                offset, ierr)
                if (ierr /= 0) exit
                pressure_old = real(aux_raw(offset + 1), real64)
            end if
            call accepted_cell_temperature(self, aux_raw, aux_section, cell, temperature_old, ierr)
            if (ierr /= 0) exit

            call PetscSectionGetFieldOffset(aux_section, cell, int(CA_MATERIAL - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr /= 0) exit
            slot = nint(real(aux_raw(offset + 1), real64), int32)

            ok = .false.
            call eval_frozen_response(slot, temperature, pressure, now, ok)
            if (ok) call eval_frozen_response(slot, temperature_old, pressure_old, accepted, ok)
            if (.not. ok) then
                ierr = COLUMN_ERR_CLOSURE
                exit
            end if

            num_facets = 0
            if (self%primal) then
                call primal_face_outflows(self, local, cell, aux_raw, aux_section, facets, moments, num_facets, ierr)
                if (ierr /= 0) exit
                do f = 1, num_facets
                    face_sum(facets(f)) = face_sum(facets(f)) + moments(f)
                    face_count(facets(f)) = face_count(facets(f)) + 1
                end do
            else
                call reconstruct_cell_facet_flux(self%dm, local, self%flux_scale, cell, facets, moments, &
                                                 num_facets, ierr, field=COLUMN_FIELD_FLUX)
                if (ierr /= 0) exit
            end if

            outflow = 0.0_real64
            do f = 1, num_facets
                outflow = outflow + moments(f)
                report%throughput = max(report%throughput, abs(moments(f)))
                call DMPlexGetSupportSize(self%dm, facets(f), support_size, ierr)
                if (ierr /= 0) exit
                if (support_size == 1_PETSC_INT_KIND) report%boundary_rate = report%boundary_rate + moments(f)
            end do
            if (ierr /= 0) exit

            call PetscSectionGetFieldOffset(aux_section, cell, int(CA_SOURCE - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr /= 0) exit
            source = real(aux_raw(offset + 1), real64)
            call measure_cell_storage(self, local, cell, slot, pressure, pressure_old, &
                                      water_mean, old_mean, ice_mean, report, ierr)
            if (ierr /= 0) exit
            defect = volume*(water_mean - old_mean) + self%accepted_dt*(outflow - volume*source)
            report%max_storage_change = max(report%max_storage_change, abs(volume*(water_mean - old_mean)))
            report%total_water_centroid = report%total_water_centroid + volume*now%stored_water
            report%max_quadrature_gap = max(report%max_quadrature_gap, abs(volume*(water_mean - now%stored_water)))
            report%max_flux_term = max(report%max_flux_term, abs(self%accepted_dt*outflow))
            report%max_cell_defect_centroid = max(report%max_cell_defect_centroid, &
                abs(volume*(now%stored_water - accepted%stored_water) + self%accepted_dt*(outflow - volume*source)))
            report%max_cell_defect = max(report%max_cell_defect, abs(defect))
            report%domain_defect = report%domain_defect + defect
            report%total_water = report%total_water + volume*water_mean
            report%total_ice = report%total_ice + volume*ice_mean
            report%owned_cells = report%owned_cells + 1
            report%coldest = min(report%coldest, temperature)
            report%warmest = max(report%warmest, temperature)
            report%lowest_liquid_pressure = min(report%lowest_liquid_pressure, now%liquid_pressure)

            report%upward_flux = max(report%upward_flux, flux_values(2))
            report%downward_flux = min(report%downward_flux, flux_values(2))

            if (now%ice_content > 0.0_real64) then
                report%frozen_cells = report%frozen_cells + 1
                report%upward_flux_frozen = max(report%upward_flux_frozen, flux_values(2))
                report%downward_flux_frozen = min(report%downward_flux_frozen, flux_values(2))
            end if

            margin = 0.0_real64
            sensitivity = 0.0_real64
            capacity = 0.0_real64
            call frozen_pressure_column(slot, temperature, pressure, self%pressure_scale, &
                                        COLUMN_STORAGE_REFERENCE, margin, sensitivity, capacity, ok)
            if (ok) report%weakest_column = min(report%weakest_column, margin)
        end do

        call VecRestoreArrayRead(self%previous_auxiliary, aux_raw, ierr)
        if (ierr /= 0) return

        if (self%primal) then
            if (any(face_count == 2)) report%max_facet_mismatch = maxval(abs(face_sum), mask=(face_count == 2))
            report%interior_facets = count(face_count == 2)
            deallocate (face_sum, face_count)
        else
            call audit_facet_continuity(self%dm, local, self%flux_scale, report%max_facet_mismatch, &
                                        report%interior_facets, ierr, field=COLUMN_FIELD_FLUX)
            if (ierr /= 0) return
        end if
        call reduce_column_report(report)
    end subroutine compute_column_report

    !> Integrate the two time levels with the assembly rule on an affine cell.
    !> Uses the stored previous accepted state; O(Nq) work and fixed memory.
    !> Returns a PETSc error on reconstruction failure. No continuum error bound is available.
    !> In the primal form both pressures are read at the quadrature points as
    !> well; `pressure` and `old_pressure` are then not used.
    subroutine measure_cell_storage(self, local, cell, slot, pressure, old_pressure, &
                                    water_mean, old_mean, ice_mean, report, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        Vec, intent(in) :: local
        PetscInt, intent(in) :: cell
        integer(int32), intent(in) :: slot
        real(real64), intent(in) :: pressure, old_pressure
        real(real64), intent(inout) :: water_mean, old_mean, ice_mean
        type(type_mixed_column_report), intent(inout) :: report
        PetscErrorCode, intent(inout) :: ierr
        real(c_double) :: points(MIXED_DIMENSION, MAX_QUADRATURE_POINTS), weights(MAX_QUADRATURE_POINTS)
        real(c_double) :: temperatures(1, MAX_QUADRATURE_POINTS), old_temperatures(1, MAX_QUADRATURE_POINTS)
        real(c_double) :: pressures(1, MAX_QUADRATURE_POINTS), old_pressures(1, MAX_QUADRATURE_POINTS)
        real(c_double) :: divergence(MAX_QUADRATURE_POINTS)
        type(c_ptr) :: dm_handle, vector_handle
        integer(c_int) :: cerr
        integer(int32) :: count, order, q
        real(real64) :: kmin, kmax, pmin, pmax
        type(type_mixed_frozen_response) :: now, old
        logical :: ok

        call cell_quadrature_rule(self%dm, COLUMN_FIELD_TEMPERATURE, points, weights, count, order, ierr)
        if (ierr /= 0) return
        dm_handle = transfer(self%dm%v, dm_handle)
        vector_handle = transfer(local%v, vector_handle)
        cerr = 0
        call petscfe_evaluate_cell(dm_handle, vector_handle, int(COLUMN_FIELD_TEMPERATURE, c_int), &
            int(cell, c_int), points(:, :count), temperatures(:, :count), divergence(:count), cerr)
        ierr = int(cerr)
        if (ierr /= 0) return
        dm_handle = transfer(self%aux_dm%v, dm_handle)
        vector_handle = transfer(self%previous_auxiliary%v, vector_handle)
        call petscfe_evaluate_cell(dm_handle, vector_handle, int(CA_TEMPERATURE_OLD - 1, c_int), &
            int(cell, c_int), points(:, :count), old_temperatures(:, :count), divergence(:count), cerr)
        ierr = int(cerr)
        if (ierr /= 0) return
        ! Both pressures in Pa: the cell's own value for P0, the field at each
        ! quadrature point for the primal space.
        pressures(1, :count) = pressure
        old_pressures(1, :count) = old_pressure
        if (self%primal) then
            call petscfe_evaluate_cell(dm_handle, vector_handle, int(CA_PRESSURE_OLD - 1, c_int), &
                int(cell, c_int), points(:, :count), old_pressures(:, :count), divergence(:count), cerr)
            ierr = int(cerr)
            if (ierr /= 0) return
            dm_handle = transfer(self%dm%v, dm_handle)
            vector_handle = transfer(local%v, vector_handle)
            call petscfe_evaluate_cell(dm_handle, vector_handle, int(self%pressure_field, c_int), &
                int(cell, c_int), points(:, :count), pressures(:, :count), divergence(:count), cerr)
            ierr = int(cerr)
            if (ierr /= 0) return
            pressures(1, :count) = self%pressure_scale*pressures(1, :count)
        end if
        water_mean = 0.0_real64
        old_mean = 0.0_real64
        ice_mean = 0.0_real64
        kmin = huge(1.0_real64)
        kmax = 0.0_real64
        pmin = huge(1.0_real64)
        pmax = -huge(1.0_real64)
        do q = 1, count
            ok = .false.
            call eval_frozen_response(slot, self%temperature_scale*temperatures(1,q), pressures(1,q), now, ok)
            if (ok) call eval_frozen_response(slot, old_temperatures(1,q), old_pressures(1,q), old, ok)
            if (.not. ok) then
                ierr = COLUMN_ERR_CLOSURE
                return
            end if
            water_mean = water_mean + weights(q)*now%stored_water
            old_mean = old_mean + weights(q)*old%stored_water
            ice_mean = ice_mean + weights(q)*now%ice_content
            kmin = min(kmin, now%relative_conductivity)
            kmax = max(kmax, now%relative_conductivity)
            pmin = min(pmin, now%liquid_pressure)
            pmax = max(pmax, now%liquid_pressure)
        end do
        report%max_conductivity_spread = max(report%max_conductivity_spread, kmax/kmin)
        report%max_liquid_pressure_spread = max(report%max_liquid_pressure_spread, pmax - pmin)
    end subroutine measure_cell_storage

    !> Write one row per owned cell: where it is, and what it holds.
    !>
    !> The cell is the unit the conservation equation closes on, so the profile
    !> is written in cells and not interpolated to nodes on the way out. A P0
    !> field projected to nodes and read back as if it were the cell value is
    !> not the quantity any balance here is exact in, and averaging one is how a
    !> comparison acquires water the run never moved. What a reader does with
    !> these rows -- depth bands, a moving average -- is its own business, but
    !> it starts from the cells and their volumes.
    !>
    !> Every value is the closure's own at the cell centroid, so the rows carry
    !> no interpolation error of their own; what they inherit is the P0 and RT0
    !> discretisation error of the state they report.
    !>
    !> \(O(N)\) in the domain's cells on rank 0, \(O(N/P)\) on every rank, with
    !> one gather of \(14N/P\) numbers per rank and one row of output per cell.
    !>
    !> Collective, and it writes the whole domain at any rank count. Each rank
    !> measures the cells it owns; the rows are gathered on rank 0 and written
    !> there in the order of the global cell number `DMPlexGetCellNumbering`
    !> gives them. The key is global, so a cell keeps its row wherever it was
    !> computed and a ghost copy is never written twice -- the local point
    !> number names a different cell on every rank and is not what is written.
    !> The file is one rank's share of the domain only if the gather loses
    !> rows, and that cannot pass unnoticed: the global numbers of the owned
    !> cells are the integers \([0,N)\) exactly once, so the ordering pass is
    !> also the cover check.
    !>
    !> Fails with `COLUMN_ERR_PROFILE_INCOMPLETE` when the gathered rows are not
    !> that one-to-one cover, with `COLUMN_ERR_NO_PROFILE_FILE` when the file
    !> cannot be opened or written, and with `COLUMN_ERR_CLOSURE` when a cell's
    !> state cannot be evaluated. Every failure is agreed on by all ranks before
    !> it is returned, so no rank leaves while another waits in a collective.
    subroutine write_column_profile(self, path, time, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> Where to write; an existing file is replaced. Read on rank 0 only.
        character(len=*), intent(in) :: path
        !> The physical time the state belongs to, s. Recorded in every row so a
        !> reader never has to infer it from a file name.
        real(real64), intent(in) :: time
        PetscErrorCode, intent(inout) :: ierr

        integer(int32), allocatable :: local_ids(:), ids(:)
        real(real64), allocatable :: local_values(:), values(:)
        integer(int32) :: num_local, total, status
        integer(int32) :: mpi_err
        PetscMPIInt :: rank

        ierr = 0
        status = 0
        if (.not. self%ready) status = 1
        call agree_on_error(status)
        if (status /= 0) then
            ierr = status
            return
        end if

        call collect_profile_rows(self, local_ids, local_values, num_local, status)
        call agree_on_error(status)
        if (status /= 0) then
            ierr = status
            return
        end if

        call gather_profile_rows(local_ids, local_values, num_local, ids, values, total, status)
        call agree_on_error(status)
        if (status /= 0) then
            ierr = status
            return
        end if

        ! Only rank 0 holds the gathered rows, so only rank 0 can find a gap in
        ! them or a file it cannot write. Both are broadcast rather than
        ! returned locally, or the other ranks would carry on as if the profile
        ! had been written.
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, mpi_err)
        if (rank == 0) call write_profile_file(path, time, ids, values, total, status)
        call MPI_Bcast(status, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, mpi_err)
        ierr = status
    end subroutine write_column_profile

    !> Measure the cells this rank owns into one packed block of rows.
    !>
    !> Ghost cells are skipped, because their owner measures them; what keys
    !> each row is the cell's global number and not its local point number, so
    !> the caller can order the domain without knowing which rank produced what.
    !> The collective calls at the head are made unconditionally so that every
    !> rank makes the same ones whatever any of them finds.
    subroutine collect_profile_rows(self, ids, values, num_rows, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        !> The global cell number of each measured cell, allocated here.
        integer(int32), allocatable, intent(inout) :: ids(:)
        !> `COLUMN_PROFILE_VALUES` numbers per row, row by row, allocated here.
        real(real64), allocatable, intent(inout) :: values(:)
        !> How many cells this rank measured.
        integer(int32), intent(inout) :: num_rows
        integer(int32), intent(inout) :: ierr

        type(type_mixed_frozen_response) :: now
        PetscSection :: aux_section
        PetscScalar, pointer :: aux_raw(:)
        PetscInt, pointer :: cell_numbers(:)
        IS :: numbering
        Vec :: local
        PetscInt :: cell, c_start, c_end, offset, global_cell
        PetscReal :: measure, centroid(3), normal(3)
        real(real64) :: volume, pressure, temperature, flux_values(MIXED_DIMENSION)
        integer(int32) :: slot, base, row
        PetscErrorCode :: status
        logical :: ok

        ierr = 0
        num_rows = 0
        c_start = 0_PETSC_INT_KIND
        c_end = 0_PETSC_INT_KIND

        call DMGlobalToLocalBegin(self%dm, self%solution, INSERT_VALUES, self%local_solution, status)
        if (ierr == 0) ierr = status
        call DMGlobalToLocalEnd(self%dm, self%solution, INSERT_VALUES, self%local_solution, status)
        if (ierr == 0) ierr = status
        call DMPlexInsertBoundaryValues(self%dm, PETSC_TRUE, self%local_solution, &
                                        0.0_PETSC_REAL_KIND, PETSC_NULL_VEC, &
                                        PETSC_NULL_VEC, PETSC_NULL_VEC, status)
        if (ierr == 0) ierr = status
        call DMGetLocalSection(self%aux_dm, aux_section, status)
        if (ierr == 0) ierr = status
        call DMPlexGetHeightStratum(self%dm, 0_PETSC_INT_KIND, c_start, c_end, status)
        if (ierr == 0) ierr = status
        call DMPlexGetCellNumbering(self%dm, numbering, status)
        if (ierr == 0) ierr = status
        if (ierr /= 0) then
            allocate (ids(0), values(0))
            return
        end if
        local = self%local_solution

        do cell = c_start, c_end - 1
            if (cell_is_owned(self%dm, cell)) num_rows = num_rows + 1
        end do
        allocate (ids(num_rows), values(COLUMN_PROFILE_VALUES*num_rows))

        call ISGetIndices(numbering, cell_numbers, ierr)
        if (ierr /= 0) return
        call VecGetArrayRead(self%accepted_auxiliary, aux_raw, ierr)
        if (ierr /= 0) then
            call ISRestoreIndices(numbering, cell_numbers, status)
            return
        end if

        row = 0
        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(self%dm, cell)) cycle

            ! An owned cell always carries a non-negative global number; a
            ! negative one is the encoding of a ghost, and this rank claiming
            ! ownership of one would put a duplicate row in the file.
            global_cell = cell_numbers(cell - c_start + 1)
            if (global_cell < 0_PETSC_INT_KIND .or. &
                global_cell > int(huge(1_int32), PETSC_INT_KIND)) then
                ierr = COLUMN_ERR_PROFILE_INCOMPLETE
                exit
            end if

            call DMPlexComputeCellGeometryFVM(self%dm, cell, measure, centroid, normal, ierr)
            if (ierr /= 0) exit
            volume = real(measure, real64)
            call reconstruct_cell_pressure(self%dm, local, self%pressure_scale, cell, pressure, ierr, &
                                           field=self%pressure_field)
            if (ierr /= 0) exit
            call cell_temperature(self, local, cell, temperature, ierr)
            if (ierr /= 0) exit
            if (self%primal) then
                call primal_cell_flux(self, local, cell, aux_raw, aux_section, flux_values, ierr)
            else
                call cell_flux(self, local, cell, flux_values, ierr)
            end if
            if (ierr /= 0) exit

            call PetscSectionGetFieldOffset(aux_section, cell, int(CA_MATERIAL - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr /= 0) exit
            slot = nint(real(aux_raw(offset + 1), real64), int32)

            ok = .false.
            call eval_frozen_response(slot, temperature, pressure, now, ok)
            if (.not. ok) then
                ierr = COLUMN_ERR_CLOSURE
                exit
            end if

            row = row + 1
            ids(row) = int(global_cell, int32)
            base = COLUMN_PROFILE_VALUES*(row - 1)
            values(base + 1) = real(centroid(1), real64)
            values(base + 2) = real(centroid(2), real64)
            values(base + 3) = volume
            values(base + 4) = temperature
            values(base + 5) = pressure
            values(base + 6) = now%liquid_pressure
            values(base + 7) = now%liquid_content
            values(base + 8) = now%ice_content
            values(base + 9) = now%vapor_content
            values(base + 10) = now%stored_water
            values(base + 11) = now%relative_conductivity
            values(base + 12) = flux_values(1)
            values(base + 13) = flux_values(2)
            call cell_water_content(self, local, cell, slot, pressure, volume, values(base + 14), ierr)
            if (ierr /= 0) exit
            call cell_vertical_transport(self, cell, values(base + 15), values(base + 16), ierr)
            if (ierr /= 0) exit
        end do

        call VecRestoreArrayRead(self%accepted_auxiliary, aux_raw, status)
        if (ierr == 0) ierr = status
        call ISRestoreIndices(numbering, cell_numbers, status)
        if (ierr == 0) ierr = status
        if (ierr == 0 .and. row /= num_rows) ierr = COLUMN_ERR_PROFILE_INCOMPLETE
    end subroutine collect_profile_rows

    !> Bring every rank's rows to rank 0, and say how many there are in all.
    !>
    !> The rows arrive in rank order and in each rank's own cell order, which is
    !> not the order they are written in; `write_profile_file` sorts them by the
    !> global number that travels with them.
    subroutine gather_profile_rows(local_ids, local_values, num_local, ids, values, total, ierr)
        implicit none
        !> The global cell number of each of this rank's rows.
        integer(int32), intent(in) :: local_ids(:)
        !> This rank's rows, `COLUMN_PROFILE_VALUES` numbers each.
        real(real64), intent(in) :: local_values(:)
        !> How many rows this rank contributes.
        integer(int32), intent(in) :: num_local
        !> Every rank's global cell numbers, on rank 0; empty elsewhere.
        integer(int32), allocatable, intent(inout) :: ids(:)
        !> Every rank's rows, on rank 0; empty elsewhere.
        real(real64), allocatable, intent(inout) :: values(:)
        !> How many rows rank 0 received; zero on every other rank.
        integer(int32), intent(inout) :: total
        integer(int32), intent(inout) :: ierr

        integer(int32), allocatable :: counts(:), starts(:), value_counts(:), value_starts(:)
        integer(int32) :: r, ranks, mpi_err, gather_err
        PetscMPIInt :: comm_size, rank

        ierr = 0
        total = 0
        gather_err = 0
        call MPI_Comm_size(PETSC_COMM_WORLD, comm_size, mpi_err)
        if (mpi_err /= 0) gather_err = mpi_err
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, mpi_err)
        if (mpi_err /= 0 .and. gather_err == 0) gather_err = mpi_err
        ranks = int(comm_size, int32)

        allocate (counts(ranks), starts(ranks), value_counts(ranks), value_starts(ranks))
        counts = 0
        call MPI_Gather(num_local, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, mpi_err)
        if (mpi_err /= 0 .and. gather_err == 0) gather_err = mpi_err

        starts(1) = 0
        do r = 2, ranks
            starts(r) = starts(r - 1) + counts(r - 1)
        end do
        value_counts = COLUMN_PROFILE_VALUES*counts
        value_starts = COLUMN_PROFILE_VALUES*starts

        if (rank == 0) then
            total = starts(ranks) + counts(ranks)
            allocate (ids(total), values(COLUMN_PROFILE_VALUES*total))
        else
            allocate (ids(0), values(0))
        end if

        call MPI_Gatherv(local_ids, num_local, MPI_INTEGER, ids, counts, starts, MPI_INTEGER, &
                         0, PETSC_COMM_WORLD, mpi_err)
        if (mpi_err /= 0 .and. gather_err == 0) gather_err = mpi_err
        call MPI_Gatherv(local_values, COLUMN_PROFILE_VALUES*num_local, MPI_DOUBLE_PRECISION, &
                         values, value_counts, value_starts, MPI_DOUBLE_PRECISION, &
                         0, PETSC_COMM_WORLD, mpi_err)
        if (mpi_err /= 0 .and. gather_err == 0) gather_err = mpi_err

        if (gather_err /= 0) ierr = COLUMN_ERR_PROFILE_INCOMPLETE
    end subroutine gather_profile_rows

    !> Write the gathered rows to one file, ordered by global cell number.
    !>
    !> Placing each row at the index its global number names is the sort and the
    !> cover check at once: the owned cells of all the ranks carry the integers
    !> \([0,N)\) exactly once, so a number out of that range, or one that
    !> arrives twice, means a ghost row was collected or a rank's share was
    !> lost. Neither leaves a partial file behind -- the check runs before the
    !> file is opened.
    subroutine write_profile_file(path, time, ids, values, total, ierr)
        implicit none
        !> Where to write; an existing file is replaced.
        character(len=*), intent(in) :: path
        !> The physical time every row belongs to, s.
        real(real64), intent(in) :: time
        !> The global cell number of each gathered row.
        integer(int32), intent(in) :: ids(:)
        !> The gathered rows, `COLUMN_PROFILE_VALUES` numbers each.
        real(real64), intent(in) :: values(:)
        !> How many rows were gathered.
        integer(int32), intent(in) :: total
        integer(int32), intent(inout) :: ierr

        integer(int32), allocatable :: order(:)
        integer(int32) :: k, v, base, unit_number, io_status

        ierr = 0
        allocate (order(total))
        order = 0
        do k = 1, total
            if (ids(k) < 0 .or. ids(k) >= total) then
                ierr = COLUMN_ERR_PROFILE_INCOMPLETE
                return
            end if
            if (order(ids(k) + 1) /= 0) then
                ierr = COLUMN_ERR_PROFILE_INCOMPLETE
                return
            end if
            order(ids(k) + 1) = k
        end do

        open (newunit=unit_number, file=path, status="replace", action="write", iostat=io_status)
        if (io_status /= 0) then
            ierr = COLUMN_ERR_NO_PROFILE_FILE
            return
        end if
        write (unit_number, '(a)', iostat=io_status) "time_s,cell,x_m,z_m,volume_m2,temperature_C,"// &
            "pressure_water_Pa,pressure_liquid_Pa,theta_w,theta_i,theta_v,stored_water,"// &
            "relative_conductivity,flux_x_m_per_s,flux_z_m_per_s,cell_water_m2,transported_up_m2,"// &
            "transported_down_m2"

        do k = 1, total
            if (io_status /= 0) exit
            base = COLUMN_PROFILE_VALUES*(order(k) - 1)
            write (unit_number, '(es22.15,a,i0,16(a,es22.15))', iostat=io_status) &
                time, ",", k - 1, (",", values(base + v), v=1, COLUMN_PROFILE_VALUES)
        end do

        close (unit_number)
        if (io_status /= 0) ierr = COLUMN_ERR_NO_PROFILE_FILE
    end subroutine write_profile_file

    !> Give every rank the same error code, so one rank's fault is not the
    !> rest's deadlock at the next collective.
    subroutine agree_on_error(status)
        implicit none
        !> This rank's code on entry, the whole communicator's largest on exit.
        integer(int32), intent(inout) :: status

        integer(int32) :: reduced, mpi_err

        call MPI_Allreduce(status, reduced, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, mpi_err)
        if (mpi_err == 0) status = reduced
    end subroutine agree_on_error

    !> The temperature at the cell's centroid, C.
    subroutine cell_temperature(self, local, cell, temperature, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        Vec, intent(in) :: local
        PetscInt, intent(in) :: cell
        real(real64), intent(inout) :: temperature
        PetscErrorCode, intent(inout) :: ierr

        call reconstruct_cell_pressure(self%dm, local, self%temperature_scale, cell, temperature, ierr, &
                                       field=COLUMN_FIELD_TEMPERATURE)
    end subroutine cell_temperature

    !> The accepted temperature at the cell's centroid, C.
    !>
    !> The mean of the cell's vertex values, which is what a P1 field is at the
    !> centroid of a simplex and what a Q1 field is at the centre of a
    !> rectangle: the bilinear term is odd about both mid-lines, so it cancels
    !> between opposite corners.
    subroutine accepted_cell_temperature(self, aux_raw, aux_section, cell, temperature, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        PetscScalar, intent(in) :: aux_raw(:)
        PetscSection, intent(in) :: aux_section
        PetscInt, intent(in) :: cell
        real(real64), intent(inout) :: temperature
        PetscErrorCode, intent(inout) :: ierr

        PetscInt, pointer :: closure(:)
        PetscInt :: point, offset, num_dofs, num_closure
        integer(int32) :: k, counted

        temperature = 0.0_real64
        counted = 0
        num_closure = 0_PETSC_INT_KIND
        nullify (closure)
        call DMPlexGetTransitiveClosure(self%aux_dm, cell, PETSC_TRUE, num_closure, closure, ierr)
        if (ierr /= 0) return

        do k = 1, 2*int(num_closure, int32), 2
            point = closure(k)
            call PetscSectionGetFieldDof(aux_section, point, int(CA_TEMPERATURE_OLD - 1, PETSC_INT_KIND), num_dofs, ierr)
            if (ierr /= 0 .or. num_dofs <= 0) cycle
            call PetscSectionGetFieldOffset(aux_section, point, int(CA_TEMPERATURE_OLD - 1, PETSC_INT_KIND), offset, ierr)
            if (ierr /= 0) cycle
            temperature = temperature + real(aux_raw(offset + 1), real64)
            counted = counted + 1
        end do

        call DMPlexRestoreTransitiveClosure(self%aux_dm, cell, PETSC_TRUE, num_closure, closure, ierr)
        if (counted > 0) temperature = temperature/real(counted, real64)
    end subroutine accepted_cell_temperature

    !> The liquid flux at the cell's centroid, m/s.
    subroutine cell_flux(self, local, cell, values, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        Vec, intent(in) :: local
        PetscInt, intent(in) :: cell
        real(real64), intent(inout) :: values(MIXED_DIMENSION)
        PetscErrorCode, intent(inout) :: ierr

        real(real64) :: moments(MAX_CELL_FACETS), volume
        PetscInt :: facets(MAX_CELL_FACETS)
        PetscReal :: measure, centroid(3), normal(3)
        PetscReal :: facet_measure, facet_centroid(3), facet_normal(3)
        integer(int32) :: num_facets, f, d

        values = 0.0_real64
        call reconstruct_cell_facet_flux(self%dm, local, self%flux_scale, cell, facets, moments, &
                                         num_facets, ierr, field=COLUMN_FIELD_FLUX)
        if (ierr /= 0) return
        call DMPlexComputeCellGeometryFVM(self%dm, cell, measure, centroid, normal, ierr)
        if (ierr /= 0) return
        volume = real(measure, real64)
        if (volume <= 0.0_real64) return

        ! The lowest-order divergence theorem reconstruction of a cell-average
        ! vector from its facet moments: sum_f Q_Kf (x_f - x_K) / |K|.
        do f = 1, num_facets
            call DMPlexComputeCellGeometryFVM(self%dm, facets(f), facet_measure, facet_centroid, facet_normal, ierr)
            if (ierr /= 0) return
            do d = 1, MIXED_DIMENSION
                values(d) = values(d) + moments(f)*(real(facet_centroid(d), real64) - real(centroid(d), real64))
            end do
        end do
        values = values/volume
    end subroutine cell_flux

    !> The primal form's liquid flux at the cell's centroid, m/s.
    !>
    !> \[
    !>   \mathbf q_l=-k_r\mathbf M\bigl(p_{l,p}\nabla p_w+p_{l,T}\nabla T-\rho_l\mathbf g\bigr)
    !> \]
    !>
    !> at the centroid, which is what the water row integrates at its
    !> quadrature points. The nodal form reports its own cell-average flux
    !> (`compute_nodal_cell_flux`) instead. There is no flux unknown to read: the two gradients
    !> are central differences across the centroid along x and z, a tenth of
    !> the cell's extent each way. P1 and Q1 are linear along those lines on an
    !> axis-aligned cell, so the difference is the field's own gradient there.
    !> \(O(1)\) per cell. Fails by returning the first PETSc error, or
    !> `COLUMN_ERR_CLOSURE` when the closure cannot answer at the centroid.
    subroutine primal_cell_flux(self, local, cell, aux_raw, aux_section, values, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        Vec, intent(in) :: local
        PetscInt, intent(in) :: cell
        !> Any auxiliary vector of the run, which carries the cell's material
        !> and saturated inverse mobility, and its section.
        PetscScalar, intent(in) :: aux_raw(:)
        PetscSection, intent(in) :: aux_section
        !> \(\mathbf q_l\), m/s; overwritten.
        real(real64), intent(inout) :: values(MIXED_DIMENSION)
        PetscErrorCode, intent(inout) :: ierr

        PetscScalar, pointer :: cell_coordinates(:)
        PetscSection :: coordinate_section
        Vec :: coordinates
        DM :: coordinate_dm
        PetscInt :: offset, closure_size
        PetscReal :: measure, centroid(3), normal(3)
        real(c_double) :: physical(MIXED_DIMENSION, 5), reference(MIXED_DIMENSION, 5)
        real(c_double) :: temperatures(1, 5), pressures(1, 5), divergence(5)
        real(real64) :: lower(MIXED_DIMENSION), upper(MIXED_DIMENSION), offsets(MIXED_DIMENSION)
        real(real64) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION)
        real(real64) :: temperature_gradient(MIXED_DIMENSION), pressure_gradient(MIXED_DIMENSION)
        type(c_ptr) :: dm_handle, vector_handle
        integer(c_int) :: cerr
        integer(int32) :: d, k, slot
        logical :: ok

        values = 0.0_real64
        if (self%nodal_flux) then
            ok = .false.
            call compute_nodal_cell_flux(self%dm, local, cell, values, ok, ierr)
            if (ierr == 0 .and. .not. ok) ierr = COLUMN_ERR_CLOSURE
            return
        end if
        call DMPlexComputeCellGeometryFVM(self%dm, cell, measure, centroid, normal, ierr)
        if (ierr /= 0) return

        ! The cell's extent along each axis, from its vertices.
        call DMGetCoordinateDM(self%dm, coordinate_dm, ierr)
        if (ierr == 0) call DMGetCoordinatesLocal(self%dm, coordinates, ierr)
        if (ierr == 0) call DMGetCoordinateSection(self%dm, coordinate_section, ierr)
        if (ierr /= 0) return
        nullify (cell_coordinates)
        closure_size = 0_PETSC_INT_KIND
        call DMPlexVecGetClosure(coordinate_dm, coordinate_section, coordinates, cell, &
                                 closure_size, cell_coordinates, ierr)
        if (ierr /= 0) return
        lower = huge(1.0_real64)
        upper = -huge(1.0_real64)
        do k = 1, int(closure_size, int32)/MIXED_DIMENSION
            do d = 1, MIXED_DIMENSION
                lower(d) = min(lower(d), real(cell_coordinates((k - 1)*MIXED_DIMENSION + d), real64))
                upper(d) = max(upper(d), real(cell_coordinates((k - 1)*MIXED_DIMENSION + d), real64))
            end do
        end do
        call DMPlexVecRestoreClosure(coordinate_dm, coordinate_section, coordinates, cell, &
                                     closure_size, cell_coordinates, ierr)
        if (ierr /= 0) return
        offsets = 0.1_real64*(upper - lower)

        ! The centroid, then a pair of points across it along x and along z.
        do k = 1, 5
            physical(:, k) = real(centroid(1:MIXED_DIMENSION), c_double)
        end do
        physical(1, 2) = physical(1, 2) + offsets(1)
        physical(1, 3) = physical(1, 3) - offsets(1)
        physical(2, 4) = physical(2, 4) + offsets(2)
        physical(2, 5) = physical(2, 5) - offsets(2)

        dm_handle = transfer(self%dm%v, dm_handle)
        vector_handle = transfer(local%v, vector_handle)
        cerr = 0
        call petscfe_reference_coordinates(dm_handle, int(cell, c_int), physical, reference, cerr)
        if (cerr == 0) call petscfe_evaluate_cell(dm_handle, vector_handle, int(COLUMN_FIELD_TEMPERATURE, c_int), &
                                                  int(cell, c_int), reference, temperatures, divergence, cerr)
        if (cerr == 0) call petscfe_evaluate_cell(dm_handle, vector_handle, int(self%pressure_field, c_int), &
                                                  int(cell, c_int), reference, pressures, divergence, cerr)
        ierr = int(cerr)
        if (ierr /= 0) return

        do d = 1, MIXED_DIMENSION
            temperature_gradient(d) = self%temperature_scale &
                                      *real(temperatures(1, 2*d) - temperatures(1, 2*d + 1), real64) &
                                      /(2.0_real64*offsets(d))
            pressure_gradient(d) = self%pressure_scale &
                                   *real(pressures(1, 2*d) - pressures(1, 2*d + 1), real64) &
                                   /(2.0_real64*offsets(d))
        end do

        call PetscSectionGetFieldOffset(aux_section, cell, int(CA_MATERIAL - 1, PETSC_INT_KIND), offset, ierr)
        if (ierr /= 0) return
        slot = nint(real(aux_raw(offset + 1), real64), int32)
        call PetscSectionGetFieldOffset(aux_section, cell, int(CA_INVERSE_MOBILITY - 1, PETSC_INT_KIND), offset, ierr)
        if (ierr /= 0) return
        do k = 1, MIXED_DIMENSION*MIXED_DIMENSION
            inverse_mobility(k) = real(aux_raw(offset + k), real64)
        end do

        ok = .false.
        call compute_primal_darcy_flux(slot, self%temperature_scale*real(temperatures(1, 1), real64), &
                                       self%pressure_scale*real(pressures(1, 1), real64), &
                                       temperature_gradient, pressure_gradient, inverse_mobility, &
                                       self%density, self%gravity, values, ok)
        if (.not. ok) ierr = COLUMN_ERR_CLOSURE
    end subroutine primal_cell_flux

    !> The water a cell holds, \(\int_K W\), m2 per unit thickness.
    !>
    !> On the quadrature the assembly uses, at the cell's own pressure for the
    !> mixed form and at the pressure field's quadrature-point values for the
    !> primal forms. \(O(N_q)\). Fails by returning the first PETSc error or
    !> `COLUMN_ERR_CLOSURE`.
    subroutine cell_water_content(self, local, cell, slot, pressure, volume, water, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        Vec, intent(in) :: local
        PetscInt, intent(in) :: cell
        integer(int32), intent(in) :: slot
        !> The cell's pressure, Pa (read by the mixed form only), and its area, m2.
        real(real64), intent(in) :: pressure, volume
        !> \(\int_K W\), m2; overwritten.
        real(real64), intent(inout) :: water
        PetscErrorCode, intent(inout) :: ierr

        real(c_double) :: points(MIXED_DIMENSION, MAX_QUADRATURE_POINTS), weights(MAX_QUADRATURE_POINTS)
        real(c_double) :: temperatures(1, MAX_QUADRATURE_POINTS), pressures(1, MAX_QUADRATURE_POINTS)
        real(c_double) :: divergence(MAX_QUADRATURE_POINTS)
        type(type_mixed_frozen_response) :: now
        type(c_ptr) :: dm_handle, vector_handle
        integer(c_int) :: cerr
        integer(int32) :: count, order, q
        real(real64) :: mean
        logical :: ok

        water = 0.0_real64
        call cell_quadrature_rule(self%dm, COLUMN_FIELD_TEMPERATURE, points, weights, count, order, ierr)
        if (ierr /= 0) return
        dm_handle = transfer(self%dm%v, dm_handle)
        vector_handle = transfer(local%v, vector_handle)
        cerr = 0
        call petscfe_evaluate_cell(dm_handle, vector_handle, int(COLUMN_FIELD_TEMPERATURE, c_int), &
            int(cell, c_int), points(:, :count), temperatures(:, :count), divergence(:count), cerr)
        pressures(1, :count) = pressure/self%pressure_scale
        if (cerr == 0 .and. self%primal) then
            call petscfe_evaluate_cell(dm_handle, vector_handle, int(self%pressure_field, c_int), &
                int(cell, c_int), points(:, :count), pressures(:, :count), divergence(:count), cerr)
        end if
        ierr = int(cerr, kind(ierr))
        if (ierr /= 0) return

        mean = 0.0_real64
        do q = 1, count
            ok = .false.
            call eval_frozen_response(slot, self%temperature_scale*temperatures(1, q), &
                                      merge(self%pressure_scale*pressures(1, q), pressure, self%primal), now, ok)
            if (.not. ok) then
                ierr = COLUMN_ERR_CLOSURE
                return
            end if
            mean = mean + weights(q)*now%stored_water
        end do
        water = volume*mean
    end subroutine cell_water_content

    !> Water the cell's own flux has carried out through its upward- and its
    !> downward-facing facets since the initial state, m2 per unit thickness.
    !>
    !> A facet faces up when its outward normal has a positive vertical
    !> component; a vertical facet counts in neither. \(O(1)\) per cell.
    !> Fails by returning the first PETSc error.
    subroutine cell_vertical_transport(self, cell, upward, downward, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        PetscInt, intent(in) :: cell
        !> Out through the upward- and the downward-facing facets, m2; overwritten.
        real(real64), intent(inout) :: upward, downward
        PetscErrorCode, intent(inout) :: ierr

        PetscInt, pointer :: cone(:)
        PetscReal :: cell_measure, cell_centroid(3), cell_normal(3)
        PetscReal :: facet_measure, facet_centroid(3), facet_normal(3)
        real(real64) :: outward_z, orientation
        integer(int32) :: f, d

        upward = 0.0_real64
        downward = 0.0_real64
        call DMPlexComputeCellGeometryFVM(self%dm, cell, cell_measure, cell_centroid, cell_normal, ierr)
        if (ierr /= 0) return
        call DMPlexGetCone(self%dm, cell, cone, ierr)
        if (ierr /= 0) return
        do f = 1, min(int(size(cone), int32), MAX_CELL_FACETS)
            call DMPlexComputeCellGeometryFVM(self%dm, cone(f), facet_measure, facet_centroid, facet_normal, ierr)
            if (ierr /= 0) exit
            orientation = 0.0_real64
            do d = 1, MIXED_DIMENSION
                orientation = orientation + real(facet_normal(d), real64) &
                              *(real(facet_centroid(d), real64) - real(cell_centroid(d), real64))
            end do
            outward_z = real(facet_normal(MIXED_DIMENSION), real64)
            if (orientation < 0.0_real64) outward_z = -outward_z
            if (outward_z > 1.0e-12_real64) upward = upward + self%transported(f, cell)
            if (outward_z < -1.0e-12_real64) downward = downward + self%transported(f, cell)
        end do
        call DMPlexRestoreCone(self%dm, cell, cone, ierr)
    end subroutine cell_vertical_transport

    !> The primal form's outward water flux through each facet of one cell, m2/s.
    !>
    !> \[
    !>   Q_{Kf}=\int_f\mathbf q_l^K\cdot\mathbf n_{Kf}\,\mathrm ds,
    !> \]
    !>
    !> with \(\mathbf q_l^K\) the flux the water row of this cell carries, from
    !> this cell's own nodal values: the chain rule \(-k_r\mathbf M(p_{l,p}\nabla
    !> p_w+p_{l,T}\nabla T-\rho_l\mathbf g)\) at each point for `primal`, and
    !> \(-\bar k_r^K\mathbf M(\nabla\hat p_l-\rho_l\mathbf g)\) for
    !> `primal_nodal`. Two-point Gauss on each facet, the rule the mixed
    !> reconstruction uses; exact for the nodal form on an affine cell, whose
    !> normal flux is linear along a facet. The flux of a continuous Galerkin
    !> solution is discontinuous across facets, so the two cells of one facet
    !> report different values; that difference is what the audit measures.
    !>
    !> Assumes an affine cell of the primal layout. \(O(1)\) per cell. Fails by
    !> returning the first PETSc error, `COLUMN_ERR_CLOSURE` when the closure
    !> cannot answer, or `1` for a cell with more facets than the arrays hold.
    subroutine primal_face_outflows(self, local, cell, aux_raw, aux_section, facets, outflows, num_facets, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        Vec, intent(in) :: local
        PetscInt, intent(in) :: cell
        !> Any auxiliary vector of the run and its section: material and mobility.
        PetscScalar, intent(in) :: aux_raw(:)
        PetscSection, intent(in) :: aux_section
        !> The cell's facets, as DMPlex points.
        PetscInt, intent(inout) :: facets(MAX_CELL_FACETS)
        !> \(Q_{Kf}\) per facet, m2/s, positive out of the cell.
        real(real64), intent(inout) :: outflows(MAX_CELL_FACETS)
        integer(int32), intent(inout) :: num_facets
        PetscErrorCode, intent(inout) :: ierr

        ! Two-point Gauss on the unit interval.
        real(real64), parameter :: NODES(2) = [0.5_real64 - 0.5_real64/sqrt(3.0_real64), &
                                               0.5_real64 + 0.5_real64/sqrt(3.0_real64)]
        real(real64), parameter :: WEIGHTS(2) = [0.5_real64, 0.5_real64]
        integer(int32), parameter :: MAX_BASIS = 4
        PetscInt, pointer :: cone(:)
        PetscScalar, pointer :: closure(:), facet_coordinates(:)
        PetscSection :: section, coordinate_section
        Vec :: coordinates
        DM :: coordinate_dm
        PetscObject :: discretisation
        PetscInt :: closure_size, offset
        PetscReal :: cell_measure, cell_centroid(3), cell_normal(3)
        PetscReal :: facet_measure, facet_centroid(3), facet_normal(3)
        PetscReal :: v0(MIXED_DIMENSION), jac(MIXED_DIMENSION*MIXED_DIMENSION)
        PetscReal :: inverse_jac(MIXED_DIMENSION*MIXED_DIMENSION), determinant
        real(c_double) :: points(MIXED_DIMENSION, 2), reference(MIXED_DIMENSION, 2)
        real(c_double) :: basis(MAX_BASIS, 2), basis_gradients(MIXED_DIMENSION, MAX_BASIS, 2)
        real(real64) :: nodal_temperature(MAX_BASIS), nodal_pressure(MAX_BASIS), nodal_liquid(MAX_BASIS)
        real(real64) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION), mobility(MIXED_DIMENSION, MIXED_DIMENSION)
        real(real64) :: gradients(MIXED_DIMENSION, MAX_BASIS), outward(MIXED_DIMENSION), orientation
        real(real64) :: vertex_a(MIXED_DIMENSION), vertex_b(MIXED_DIMENSION)
        real(real64) :: temperature, pressure, temperature_gradient(MIXED_DIMENSION)
        real(real64) :: pressure_gradient(MIXED_DIMENSION), liquid_gradient(MIXED_DIMENSION), flux(MIXED_DIMENSION)
        real(real64) :: conductivity, log_mean, det_mobility
        type(type_mixed_frozen_response) :: response
        type(c_ptr) :: dm_handle, fe_handle
        integer(c_int) :: cerr, num_basis, num_components
        integer(int32) :: f, q, i, d, e, k, slot, n
        logical :: ok

        outflows = 0.0_real64
        facets = 0_PETSC_INT_KIND
        num_facets = 0

        call DMPlexGetCone(self%dm, cell, cone, ierr)
        if (ierr /= 0) return
        num_facets = int(size(cone), int32)
        if (num_facets > MAX_CELL_FACETS) then
            call DMPlexRestoreCone(self%dm, cell, cone, ierr)
            num_facets = 0
            ierr = 1
            return
        end if
        facets(1:num_facets) = cone(1:num_facets)
        call DMPlexRestoreCone(self%dm, cell, cone, ierr)
        if (ierr /= 0) return

        call DMPlexComputeCellGeometryFVM(self%dm, cell, cell_measure, cell_centroid, cell_normal, ierr)
        if (ierr == 0) call DMPlexComputeCellGeometryAffineFEM(self%dm, cell, v0, jac, inverse_jac, determinant, ierr)
        if (ierr /= 0) return

        ! The cell's material, mobility and nodal state.
        call PetscSectionGetFieldOffset(aux_section, cell, int(CA_MATERIAL - 1, PETSC_INT_KIND), offset, ierr)
        if (ierr /= 0) return
        slot = nint(real(aux_raw(offset + 1), real64), int32)
        call PetscSectionGetFieldOffset(aux_section, cell, int(CA_INVERSE_MOBILITY - 1, PETSC_INT_KIND), offset, ierr)
        if (ierr /= 0) return
        do k = 1, MIXED_DIMENSION*MIXED_DIMENSION
            inverse_mobility(k) = real(aux_raw(offset + k), real64)
        end do
        det_mobility = inverse_mobility(1)*inverse_mobility(4) - inverse_mobility(2)*inverse_mobility(3)
        mobility(1, 1) = inverse_mobility(4)/det_mobility
        mobility(1, 2) = -inverse_mobility(2)/det_mobility
        mobility(2, 1) = -inverse_mobility(3)/det_mobility
        mobility(2, 2) = inverse_mobility(1)/det_mobility

        dm_handle = transfer(self%dm%v, dm_handle)
        cerr = 0
        call petscfe_field_size(dm_handle, 0_c_int, num_basis, num_components, cerr)
        ierr = int(cerr, kind(ierr))
        if (ierr /= 0) return
        n = int(num_basis, int32)
        if (n > MAX_BASIS) then
            ierr = 1
            return
        end if
        call DMGetField(self%dm, 0_PETSC_INT_KIND, PETSC_NULL_DMLABEL, discretisation, ierr)
        if (ierr /= 0) return
        fe_handle = transfer(discretisation%v, fe_handle)

        ! The closure is field-major: the temperature's nodes, then the pressure's.
        call DMGetLocalSection(self%dm, section, ierr)
        if (ierr /= 0) return
        nullify (closure)
        closure_size = 0_PETSC_INT_KIND
        call DMPlexVecGetClosure(self%dm, section, local, cell, closure_size, closure, ierr)
        if (ierr /= 0) return
        do i = 1, n
            nodal_temperature(i) = self%temperature_scale*real(closure(i), real64)
            nodal_pressure(i) = self%pressure_scale*real(closure(n + i), real64)
        end do
        call DMPlexVecRestoreClosure(self%dm, section, local, cell, closure_size, closure, ierr)
        if (ierr /= 0) return

        ! The nodal form's liquid pressure and geometric-mean conductivity.
        conductivity = 0.0_real64
        if (self%nodal_flux) then
            log_mean = 0.0_real64
            do i = 1, n
                ok = .false.
                call eval_frozen_response(slot, nodal_temperature(i), nodal_pressure(i), response, ok)
                if (.not. ok) then
                    ierr = COLUMN_ERR_CLOSURE
                    return
                end if
                nodal_liquid(i) = response%liquid_pressure
                log_mean = log_mean + log(response%relative_conductivity)
            end do
            conductivity = exp(log_mean/real(n, real64))
        end if

        call DMGetCoordinateDM(self%dm, coordinate_dm, ierr)
        if (ierr == 0) call DMGetCoordinatesLocal(self%dm, coordinates, ierr)
        if (ierr == 0) call DMGetCoordinateSection(self%dm, coordinate_section, ierr)
        if (ierr /= 0) return

        do f = 1, num_facets
            call DMPlexComputeCellGeometryFVM(self%dm, facets(f), facet_measure, facet_centroid, facet_normal, ierr)
            if (ierr /= 0) return
            ! Out of this cell: the side of the facet away from the cell's centroid.
            orientation = 0.0_real64
            do d = 1, MIXED_DIMENSION
                orientation = orientation + real(facet_normal(d), real64) &
                              *(real(facet_centroid(d), real64) - real(cell_centroid(d), real64))
            end do
            outward = real(facet_normal(1:MIXED_DIMENSION), real64)
            if (orientation < 0.0_real64) outward = -outward

            nullify (facet_coordinates)
            closure_size = 0_PETSC_INT_KIND
            call DMPlexVecGetClosure(coordinate_dm, coordinate_section, coordinates, facets(f), &
                                     closure_size, facet_coordinates, ierr)
            if (ierr /= 0) return
            do d = 1, MIXED_DIMENSION
                vertex_a(d) = real(facet_coordinates(d), real64)
                vertex_b(d) = real(facet_coordinates(MIXED_DIMENSION + d), real64)
            end do
            call DMPlexVecRestoreClosure(coordinate_dm, coordinate_section, coordinates, facets(f), &
                                         closure_size, facet_coordinates, ierr)
            if (ierr /= 0) return
            do q = 1, 2
                points(:, q) = real(vertex_a + NODES(q)*(vertex_b - vertex_a), c_double)
            end do

            cerr = 0
            call petscfe_reference_coordinates(dm_handle, int(cell, c_int), points, reference, cerr)
            if (cerr == 0) call petscfe_tabulate(fe_handle, reference, basis(:n, :), basis_gradients(:, :n, :), cerr)
            ierr = int(cerr, kind(ierr))
            if (ierr /= 0) return

            do q = 1, 2
                ! Physical basis gradients: d(phi)/dx_d = sum_e d(phi)/d(xi_e) invJ(e, d).
                do i = 1, n
                    do d = 1, MIXED_DIMENSION
                        gradients(d, i) = 0.0_real64
                        do e = 1, MIXED_DIMENSION
                            gradients(d, i) = gradients(d, i) + real(basis_gradients(e, i, q), real64) &
                                              *real(inverse_jac((e - 1)*MIXED_DIMENSION + d), real64)
                        end do
                    end do
                end do
                if (self%nodal_flux) then
                    liquid_gradient = matmul(gradients(:, :n), nodal_liquid(:n))
                    flux = -conductivity*matmul(mobility, liquid_gradient - self%density*self%gravity)
                else
                    temperature = dot_product(real(basis(:n, q), real64), nodal_temperature(:n))
                    pressure = dot_product(real(basis(:n, q), real64), nodal_pressure(:n))
                    temperature_gradient = matmul(gradients(:, :n), nodal_temperature(:n))
                    pressure_gradient = matmul(gradients(:, :n), nodal_pressure(:n))
                    ok = .false.
                    call compute_primal_darcy_flux(slot, temperature, pressure, temperature_gradient, &
                                                   pressure_gradient, inverse_mobility, self%density, &
                                                   self%gravity, flux, ok)
                    if (.not. ok) then
                        ierr = COLUMN_ERR_CLOSURE
                        return
                    end if
                end if
                outflows(f) = outflows(f) + WEIGHTS(q)*real(facet_measure, real64)*dot_product(flux, outward)
            end do
        end do
    end subroutine primal_face_outflows

    !> Reduce one rank's audit into the domain's.
    subroutine reduce_column_report(report)
        implicit none
        type(type_mixed_column_report), intent(inout) :: report

        call reduce_scalar(report%total_water, MPI_SUM)
        call reduce_scalar(report%total_ice, MPI_SUM)
        call reduce_scalar(report%total_water_centroid, MPI_SUM)
        call reduce_scalar(report%max_storage_change, MPI_MAX)
        call reduce_scalar(report%max_flux_term, MPI_MAX)
        call reduce_scalar(report%max_quadrature_gap, MPI_MAX)
        call reduce_scalar(report%max_cell_defect_centroid, MPI_MAX)
        call reduce_scalar(report%max_conductivity_spread, MPI_MAX)
        call reduce_scalar(report%max_liquid_pressure_spread, MPI_MAX)
        call reduce_scalar(report%domain_defect, MPI_SUM)
        call reduce_scalar(report%boundary_rate, MPI_SUM)
        call reduce_scalar(report%max_cell_defect, MPI_MAX)
        call reduce_scalar(report%throughput, MPI_MAX)
        call reduce_scalar(report%upward_flux_frozen, MPI_MAX)
        call reduce_scalar(report%downward_flux_frozen, MPI_MIN)
        call reduce_scalar(report%upward_flux, MPI_MAX)
        call reduce_scalar(report%downward_flux, MPI_MIN)
        call reduce_scalar(report%weakest_column, MPI_MIN)
        call reduce_scalar(report%coldest, MPI_MIN)
        call reduce_scalar(report%warmest, MPI_MAX)
        call reduce_scalar(report%lowest_liquid_pressure, MPI_MIN)
        call reduce_integer(report%frozen_cells)
        call reduce_integer(report%owned_cells)
    end subroutine reduce_column_report

    subroutine reduce_scalar(value, operation)
        implicit none
        real(real64), intent(inout) :: value
        integer, intent(in) :: operation

        real(real64) :: reduced
        integer(int32) :: mpi_err

        call MPI_Allreduce(value, reduced, 1, MPI_DOUBLE_PRECISION, operation, PETSC_COMM_WORLD, mpi_err)
        value = reduced
    end subroutine reduce_scalar

    subroutine reduce_integer(value)
        implicit none
        integer(int32), intent(inout) :: value

        real(real64) :: sent, reduced
        integer(int32) :: mpi_err

        sent = real(value, real64)
        call MPI_Allreduce(sent, reduced, 1, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, mpi_err)
        value = nint(reduced)
    end subroutine reduce_integer

    !> The DM the three fields live on.
    function get_column_dm(self) result(dm)
        implicit none
        class(type_mixed_column), intent(in) :: self
        DM :: dm

        dm = self%dm
    end function get_column_dm

    !> Release every PETSc object the solver made.
    subroutine destroy_mixed_column(self, ierr)
        implicit none
        class(type_mixed_column), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        integer(int32) :: b

        ierr = 0
        if (self%laws_installed) call clear_frozen_laws()
        self%laws_installed = .false.
        if (self%nodal_flux) call clear_nodal_flux()
        call close_column_step_history(self)
        if (.not. self%ready) return

        call SNESDestroy(self%snes, ierr)
        call VecDestroy(self%solution, ierr)
        call VecDestroy(self%old_solution, ierr)
        call VecDestroy(self%local_solution, ierr)
        call VecDestroy(self%accepted_auxiliary, ierr)
        call VecDestroy(self%previous_auxiliary, ierr)
        do b = 1, self%num_heat_transfers
            call VecDestroy(self%heat_auxiliary(b), ierr)
        end do
        call VecDestroy(self%auxiliary, ierr)
        call DMDestroy(self%aux_dm, ierr)
        call DMDestroy(self%dm, ierr)

        self%num_closed = 0
        self%num_temperatures = 0
        self%num_heat_transfers = 0
        self%ready = .false.
    end subroutine destroy_mixed_column

end module app_mixed_hydraulic_column
