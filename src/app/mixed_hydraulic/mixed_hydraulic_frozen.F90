!> The frozen constitutive contract of the mixed hydraulic solver.
!>
!> One map, evaluated cell-locally,
!>
!> \[
!>   (T,p_w,\text{material})\longmapsto
!>   \bigl(W,\;U,\;p_l,\;k_r\bigr),
!> \]
!>
!> together with its tangents in both unknowns. Every value is produced by the
!> repository's own governing models -- `type_hydraulic` and `type_thermal`, the
!> same objects the monolithic path holds -- so the mixed solver evaluates the
!> same freezing closure as the rest of FTCMS rather than a second statement of
!> it. Nothing here restates a retention curve, a freezing curve, an impedance
!> law or a storage formula.
!>
!> ## The two pressures
!>
!> The solved P0 unknown is \(p_w\), the **total-water-equivalent** pore
!> pressure: the closure defines the total stored water as the unfrozen
!> equivalent content at that pressure, and it is \(p_w\) the mass balance
!> transports. It is not the pressure of the liquid phase.
!>
!> The liquid phase sees \(p_l=-\psi_{\mathrm{eff}}\), where
!> \(\psi_{\mathrm{eff}}\) is the effective suction the phase split publishes: a
!> smooth maximum of the matric suction of \(p_w\) and the freezing suction of
!> the generalized Clapeyron relation, limited by the pore volume. Freezing
!> raises that suction, and \(\nabla p_l\) is what draws water toward a front --
!> the cryosuction the Mizoguchi column measures.
!>
!> Cryosuction therefore enters this contract exactly once, through
!> \(\psi_{\mathrm{eff}}\), and it enters in two places that are the same
!> number: the Darcy driver \(p_l\), and the head at which the relative
!> conductivity is evaluated (`calc_head_hcf` reads \(\psi_{\mathrm{eff}}\) off
!> the state when it is set). There is no second cryosuction term anywhere in
!> this module, no thermo-osmotic contribution added on top of it, and no
!> inverse map \(p_l\to p_w\): the unknown stays \(p_w\), and \(p_l\) is a
!> published function of it. A solver that drove the flux with \(\nabla p_w\)
!> and then added a freezing-suction term of its own would count the same
!> physics twice.
!>
!> Writing the liquid driver on the unknowns by the chain rule,
!>
!> \[
!>   \nabla p_l=\frac{\partial p_l}{\partial p_w}\nabla p_w
!>             +\frac{\partial p_l}{\partial T}\nabla T,\qquad
!>   \frac{\partial p_l}{\partial p_w}=-\frac{\partial\psi_{\mathrm{eff}}}{\partial p_w},
!> \]
!>
!> is what the tangents below report, so a consumer composes one driver and not
!> two.
!>
!> ## The stored quantity
!>
!> \(W\) is `compute_total_storage`,
!>
!> \[
!>   W=\theta_w+\frac{\rho_i}{\rho_w}\theta_i+\theta_v+S_c,\qquad
!>   S_c=\frac{S_s}{\rho_{\mathrm{std}}g}\max(0,p_w-p_{\mathrm{sat}}),
!> \]
!>
!> Tangents of W are differenced from this value; the existing Picard capacity
!> is not assumed to be its exact derivative. This is the full frozen quantity: the ice term through its
!> density ratio, the vapour term when the run carries it, and the compressive
!> storage. The isothermal retention-plus-compressive law of
!> `app_mixed_hydraulic_closure` is its unfrozen limit and is not substituted
!> for it here; the two agree where there is no ice, which is what
!> `mixed_freeze_closure` measures.
!>
!> ## Ownership
!>
!> The governing models are module state for the same reason the isothermal
!> tables are: a PetscDS pointwise function is reached through a C function
!> pointer and carries no typed context. One solver per rank owns them, and
!> installing a second while the first holds them is refused.
module app_mixed_hydraulic_frozen
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use :: core_physical_constants, only:gravity_acceleration
    use :: module_core, only:type_state, type_config_constitutive, type_config_wrf, &
                              type_config_hcf, type_config_gcc
    use :: constitutive_constants, only:celsius_to_kelvin, &
        rh_gas_constant => universal_gas_constant, rh_molar_mass_water => molar_mass_water, &
        rh_reference_water_density => reference_water_density
    use :: physics_governing_thermal, only:type_thermal
    use :: physics_governing_hydraulic, only:type_hydraulic
    implicit none
    private

    public :: type_mixed_frozen_material
    public :: type_mixed_frozen_response
    public :: type_mixed_frozen_tangents
    public :: install_frozen_laws
    public :: clear_frozen_laws
    public :: frozen_laws_installed
    public :: publish_frozen_state
    public :: eval_frozen_response
    public :: eval_frozen_tangents
    public :: frozen_pressure_column
    public :: configure_frozen_tangent_cache
    public :: frozen_work_counters
    public :: MAX_FROZEN_LAWS

    !> How many materials one run may carry frozen laws for.
    integer(int32), parameter :: MAX_FROZEN_LAWS = 16

    !> Installation probe pressure, Pa; not a conductivity normalisation state.
    real(real64), parameter :: REFERENCE_PRESSURE = 0.0_real64
    !> Installation probe temperature, C. No unfrozen-state guarantee is assumed.
    real(real64), parameter :: REFERENCE_TEMPERATURE = 20.0_real64

    !> Temperature offset the enthalpy and conductivity tangents are differenced
    !> over, K. The models expose no derivative of either in this direction.
    real(real64), parameter :: TEMPERATURE_STEP = 1.0e-4_real64
    !> Pressure offset of the same differences, Pa.
    real(real64), parameter :: PRESSURE_STEP = 1.0_real64

    !> Relative departure from isotropy a thermal conductivity tensor may carry.
    real(real64), parameter :: ISOTROPY_TOLERANCE = 1.0e-12_real64

    !> One material's frozen constitutive definition.
    !>
    !> The configurations are exactly those the input translator produces for a
    !> project material, so a driver without a project describes a material the
    !> same way the application does.
    type :: type_mixed_frozen_material
        !> The region identifier the governing models key this material by.
        integer(int32) :: id = 0
        !> Porosity of the material, the pore volume the phases share.
        real(real64) :: porosity = 0.0_real64
        !> Specific storage \(S_s\) of the saturated branch, 1/m.
        real(real64) :: specific_storage = 0.0_real64
        !> Phase densities.
        type(type_config_constitutive) :: density
        !> Phase specific heats.
        type(type_config_constitutive) :: specific_heat
        !> Volumetric heat capacities.
        type(type_config_constitutive) :: heat_capacity
        !> Phase thermal conductivities.
        type(type_config_constitutive) :: thermal_conductivity
        !> Retention curve.
        type(type_config_wrf) :: retention
        !> Hydraulic conductivity function, impedance mode included.
        type(type_config_hcf) :: conductivity_function
        !> Freezing curve.
        type(type_config_gcc) :: freezing
    end type type_mixed_frozen_material

    !> What one cell's state is, at one \((T,p_w)\).
    type :: type_mixed_frozen_response
        !> \(W\), total stored water per unit bulk volume, dimensionless.
        real(real64) :: stored_water = 0.0_real64
        !> \(U\), enthalpy density, J/m3.
        real(real64) :: enthalpy = 0.0_real64
        !> \(p_l=-\psi_{\mathrm{eff}}\), the pressure the liquid flux follows, Pa.
        real(real64) :: liquid_pressure = 0.0_real64
        !> \(k_r\), relative conductivity, dimensionless and positive.
        real(real64) :: relative_conductivity = 0.0_real64
        !> \(\lambda\), bulk thermal conductivity, W/(m K).
        !>
        !> The scalar of an isotropic tensor. `install_frozen_laws` refuses a
        !> material whose thermal conductivity tensor is not isotropic rather
        !> than silently reporting one of its entries.
        real(real64) :: thermal_conductivity = 0.0_real64
        !> \(\theta_w\), liquid water content.
        real(real64) :: liquid_content = 0.0_real64
        !> \(\theta_i\), ice content.
        real(real64) :: ice_content = 0.0_real64
        !> \(\theta_v\), vapour content.
        real(real64) :: vapor_content = 0.0_real64
        !> \(\rho_w\), liquid density at this state, kg/m3.
        !>
        !> Reported because \(W\) counts ice through \(\rho_i/\rho_w\): the
        !> conversion is part of the contract and has to be auditable from
        !> outside it.
        real(real64) :: water_density = 0.0_real64
        !> \(\rho_i\), ice density at this state, kg/m3.
        real(real64) :: ice_density = 0.0_real64
        !> \(S_c\), the compressive part of \(W\), dimensionless.
        !>
        !> Read as \(W\) minus the phase sum rather than from a second storage
        !> formula, so it is by construction the term the residual carries.
        real(real64) :: compressive_storage = 0.0_real64
        !> Whether ice is present, \(\theta_i>0\).
        logical :: ice_present = .false.
        !> Whether the compressive branch \(\max(0,p_w-p_{\mathrm{sat}})\) is open.
        logical :: compressive_active = .false.
    end type type_mixed_frozen_response

    !> The tangents of that state in both unknowns.
    type :: type_mixed_frozen_tangents
        !> \(\partial W/\partial p_w\), 1/Pa.
        real(real64) :: dwater_dpressure = 0.0_real64
        !> \(\partial W/\partial T\), 1/K.
        real(real64) :: dwater_dtemperature = 0.0_real64
        !> \(\partial U/\partial T\), the apparent heat capacity, J/(m3 K).
        real(real64) :: denthalpy_dtemperature = 0.0_real64
        !> \(\partial U/\partial p_w\), J/(m3 Pa).
        real(real64) :: denthalpy_dpressure = 0.0_real64
        !> \(\partial p_l/\partial p_w\), dimensionless.
        real(real64) :: dliquid_dpressure = 0.0_real64
        !> \(\partial p_l/\partial T\), Pa/K.
        real(real64) :: dliquid_dtemperature = 0.0_real64
        !> \(\partial\ln k_r/\partial p_w\), 1/Pa.
        real(real64) :: dlog_conductivity_dpressure = 0.0_real64
        !> \(\partial\ln k_r/\partial T\), 1/K.
        real(real64) :: dlog_conductivity_dtemperature = 0.0_real64
        !> \(\partial\lambda/\partial p_w\), W/(m K Pa).
        real(real64) :: dthermal_dpressure = 0.0_real64
        !> \(\partial\lambda/\partial T\), W/(m K K).
        real(real64) :: dthermal_dtemperature = 0.0_real64
        !> Whether the pressure differences were taken one-sided, because an
        !> offset left the branch the point is on.
        logical :: pressure_one_sided = .false.
        !> The same for the temperature differences.
        logical :: temperature_one_sided = .false.
        !> Whether a branch change was found on BOTH sides of the point, so
        !> neither one-sided difference stays on the point's own branch and the
        !> reported tangent is a straddling central difference.
        !>
        !> A consumer that needs a Newton tangent has to treat this as a state
        !> it cannot linearise about, not as an ordinary number.
        logical :: branch_straddled = .false.
    end type type_mixed_frozen_tangents

    ! The governing models of the run, and the per-slot data a kernel reaches
    ! its own material with.
    type(type_hydraulic), save :: hydraulic_model
    type(type_thermal), save :: thermal_model
    integer(int32), save :: material_ids(MAX_FROZEN_LAWS) = 0
    real(real64), save :: material_porosity(MAX_FROZEN_LAWS) = 0.0_real64
    real(real64), save :: saturated_conductivity(MAX_FROZEN_LAWS) = 1.0_real64
    integer(int32), save :: num_frozen_laws = 0

    ! A small exactly-keyed cache of published states. Assembling a Jacobian
    ! visits one point with the base state and the four offsets the tangents
    ! difference, so that is the working set of one point.
    integer(int32), parameter :: STATE_SLOTS = 5
    type(type_state), save :: slot_state(STATE_SLOTS)
    real(real64), save :: slot_temperature(STATE_SLOTS) = 0.0_real64
    real(real64), save :: slot_pressure(STATE_SLOTS) = 0.0_real64
    integer(int32), save :: slot_material(STATE_SLOTS) = 0
    logical, save :: slot_filled(STATE_SLOTS) = .false.
    integer(int32), save :: slot_next = 1

    ! How much constitutive work was asked for and how much was actually done.
    ! Nothing branches on these; they exist so a profile can tell a request
    ! apart from an evaluation.
    integer(int64), save :: tangent_requests = 0_int64
    integer(int64), save :: tangent_evaluations = 0_int64
    integer(int64), save :: response_requests = 0_int64
    integer(int64), save :: response_evaluations = 0_int64
    integer(int64), save :: state_publications = 0_int64

    ! A memo of computed tangents, keyed on the exact (slot, T, p) a kernel asks
    ! with.
    !
    ! Assembling the column Jacobian calls eval_frozen_tangents once per
    ! registered block per quadrature point -- twelve times -- with bit-identical
    ! arguments, and the map is a deterministic function of its key alone, so
    ! every call after the first can only recompute the number the first one
    ! already produced. PETSc integrates one field block over every cell before
    ! it begins the next, so a point is revisited only after a whole sweep of the
    ! rank's cells: the memo is sized for that working set. STATE_SLOTS above
    ! covers the five states inside one tangent evaluation and cannot span
    ! blocks, which is why it does not remove this repetition by itself.
    !
    ! Direct-mapped, so a colliding key overwrites: that costs a recomputation
    ! and never a wrong answer. Only successful evaluations are stored, so a
    ! point the closure refuses is refused again on every visit and its fault is
    ! reported by the caller exactly as it was before.
    !
    ! Like slot_state above, this is module state without a lock, so it inherits
    ! the same single-threaded-assembly assumption the published states already
    ! make. It adds no sharing that was not there.
    integer(int32), parameter :: TANGENT_CACHE_DEFAULT_ENTRIES = 65536
    logical, save :: tangent_cache_enabled = .true.
    integer(int32), save :: tangent_cache_mask = 0
    type(type_mixed_frozen_tangents), allocatable, save :: tangent_cache_value(:)
    real(real64), allocatable, save :: tangent_cache_temperature(:)
    real(real64), allocatable, save :: tangent_cache_pressure(:)
    integer(int32), allocatable, save :: tangent_cache_material(:)
    logical, allocatable, save :: tangent_cache_filled(:)

    ! The same memo for the response, which the residual repeats for the same
    ! reason: its five registered kernels each ask for one at the same key. It
    ! is a separate table because the two are asked for at different keys -- the
    ! tangents also evaluate the four offsets -- and sharing one would evict the
    ! base states the residual reuses.
    type(type_mixed_frozen_response), allocatable, save :: response_cache_value(:)
    real(real64), allocatable, save :: response_cache_temperature(:)
    real(real64), allocatable, save :: response_cache_pressure(:)
    integer(int32), allocatable, save :: response_cache_material(:)
    logical, allocatable, save :: response_cache_filled(:)

contains

    !> Build the frozen closure of every material, in the order given.
    !>
    !> A material's slot is its index in `materials`, which is the number an
    !> auxiliary material field carries; its `id` is what the governing models
    !> key it by. Costs one governing-model construction for the whole table and
    !> one saturated conductivity evaluation per material, and \(O(1)\) memory
    !> each.
    !>
    !> Fails by `error stop` when the table is already held by a live solver,
    !> when there are more materials than slots, when an identifier is not
    !> positive, or when a material has no positive finite saturated
    !> conductivity. Nothing is floored: an impermeable material is a different
    !> model, not a small number.
    subroutine install_frozen_laws(materials, computation_type, enable_vapor_transport)
        implicit none
        !> One entry per material of the mesh.
        type(type_mixed_frozen_material), intent(in) :: materials(:)
        !> The computation type the governing tensors are laid out for; the
        !> mixed solver's intrinsic plane is the x-z one.
        integer(int32), intent(in) :: computation_type
        !> Whether the vapour terms take part in storage, as the project's
        !> analysis controls would say.
        logical, intent(in) :: enable_vapor_transport

        integer(int32) :: m, count
        integer(int32), allocatable :: ids(:)
        real(real64), allocatable :: storages(:)
        type(type_config_constitutive), allocatable :: densities(:), specific_heats(:)
        type(type_config_constitutive), allocatable :: capacities(:), conductivities(:)
        type(type_config_wrf), allocatable :: retentions(:)
        type(type_config_hcf), allocatable :: functions(:)
        type(type_config_gcc), allocatable :: freezings(:)
        type(type_mixed_frozen_response) :: reference
        logical :: ok

        if (num_frozen_laws /= 0) then
            error stop "A second mixed hydraulic solver cannot share the frozen closure of a live one"
        end if
        count = int(size(materials), int32)
        if (count > MAX_FROZEN_LAWS) then
            error stop "More materials than the mixed hydraulic frozen closure can hold"
        end if
        if (count < 1) then
            error stop "The mixed hydraulic frozen closure needs at least one material"
        end if

        ! This table renumbers the slots the memo keys on, and its laws are new.
        call forget_frozen_tangents()

        allocate (ids(count), storages(count))
        allocate (densities(count), specific_heats(count), capacities(count), conductivities(count))
        allocate (retentions(count), functions(count), freezings(count))

        do m = 1, count
            if (materials(m)%id < 1) then
                error stop "A mixed hydraulic frozen material needs a positive region identifier"
            end if
            if (materials(m)%porosity <= 0.0_real64 .or. materials(m)%porosity >= 1.0_real64) then
                error stop "A mixed hydraulic frozen material needs a porosity strictly between zero and one"
            end if
            ids(m) = materials(m)%id
            storages(m) = materials(m)%specific_storage
            densities(m) = materials(m)%density
            specific_heats(m) = materials(m)%specific_heat
            capacities(m) = materials(m)%heat_capacity
            conductivities(m) = materials(m)%thermal_conductivity
            retentions(m) = materials(m)%retention
            functions(m) = materials(m)%conductivity_function
            freezings(m) = materials(m)%freezing
            material_ids(m) = materials(m)%id
            material_porosity(m) = materials(m)%porosity
            saturated_conductivity(m) = 1.0_real64
        end do

        ! The computation dimension of the mixed solver's intrinsic plane. The
        ! governing models fill their tensors over it; this contract reads only
        ! scalars from them, but the models are built for one all the same.
        call thermal_model%initialize_from_configs(ids, densities, specific_heats, capacities, &
                                                   conductivities, retentions, functions, freezings, &
                                                   computation_type, 2_int32, enable_vapor_transport)
        call hydraulic_model%initialize_from_configs(ids, densities, retentions, functions, freezings, &
                                                     storages, computation_type, 2_int32, &
                                                     enable_vapor_transport)
        num_frozen_laws = count

        ! Normalise by the configured saturated conductivity, not a phase-split
        ! state at p=0: the regularised phase split can contain ice there even
        ! above zero Celsius, which would spuriously amplify every mobility.
        do m = 1, count
            saturated_conductivity(m) = materials(m)%conductivity_function%k_sat
            if (.not. ieee_is_finite(saturated_conductivity(m))) then
                error stop "Nonfinite saturated conductivity in frozen closure"
            end if
            if (saturated_conductivity(m) <= 0.0_real64) then
                error stop "Nonpositive saturated conductivity in frozen closure"
            end if
        end do

        ! One evaluation of the reference state, so an installation that cannot
        ! answer at all fails here rather than inside an assembly.
        call eval_frozen_response(1_int32, REFERENCE_TEMPERATURE, REFERENCE_PRESSURE, reference, ok)
        if (.not. ok) then
            num_frozen_laws = 0
            error stop "The mixed hydraulic frozen closure could not evaluate its own reference state"
        end if
    end subroutine install_frozen_laws

    !> Release the closure, so the next solver starts from nothing.
    subroutine clear_frozen_laws()
        implicit none

        if (num_frozen_laws > 0) then
            call hydraulic_model%destroy()
            call thermal_model%destroy()
        end if
        material_ids(:) = 0
        material_porosity(:) = 0.0_real64
        saturated_conductivity(:) = 1.0_real64
        num_frozen_laws = 0
        slot_filled(:) = .false.
        slot_next = 1
        ! The memo keys on a material slot, and the next run renumbers those.
        call forget_frozen_tangents()
    end subroutine clear_frozen_laws

    !> Whether a run currently holds the frozen closure.
    logical function frozen_laws_installed()
        implicit none

        frozen_laws_installed = (num_frozen_laws > 0)
    end function frozen_laws_installed

    !> Publish the phase split of one cell onto a caller-supplied state.
    !>
    !> The state comes back carrying the contents, the effective suction and the
    !> analytic suction and storage tangents the closure computed, which is the
    !> whole of what the freezing model published for \((T,p_w)\). \(O(1)\), one
    !> phase-split evaluation.
    !>
    !> Fails with `ok = .false.` and an untouched state for a slot outside the
    !> table, and for a state outside the domain the freezing model is defined
    !> on: a nonfinite unknown, or a temperature at or below absolute zero.
    subroutine publish_frozen_state(slot, temperature, pressure, state, ok)
        implicit none
        !> One-based index into the frozen closure's material table.
        integer(int32), intent(in) :: slot
        !> Temperature, C.
        real(real64), intent(in) :: temperature
        !> Total-water-equivalent pore pressure \(p_w\), Pa.
        real(real64), intent(in) :: pressure
        !> The state to publish onto; overwritten.
        type(type_state), intent(inout) :: state
        !> True when the slot is one the closure holds.
        logical, intent(inout) :: ok

        ok = .false.
        if (slot < 1 .or. slot > num_frozen_laws) return

        ! A Newton iterate is not a physical state, and this path is reached at
        ! every iterate. The freezing model states its domain as a hard
        ! invariant -- it divides by the absolute temperature to form the
        ! relative humidity and stops the program below it -- so an iterate that
        ! leaves the domain has to be refused here, where refusing is a return
        ! value the step controller acts on, rather than there, where it ends
        ! the run. The state is refused and never moved: nothing is clamped onto
        ! the boundary of the domain, because a clamped state would be reported
        ! as a converged answer for a state the solver never visited.
        if (.not. state_is_admissible(temperature, pressure)) return

        state_publications = state_publications + 1_int64
        call state%temperature%set(temperature)
        call state%pressure%set(pressure)
        call state%porosity%set(material_porosity(slot))
        ! Full freezing physics: this path runs no homotopy continuation.
        call state%continuation_lambda%set(1.0_real64)
        call hydraulic_model%update_water_phases(material_ids(slot), state)
        ok = .true.
    end subroutine publish_frozen_state

    !> The frozen constitutive response of one cell.
    !>
    !> \[
    !>   (T,p_w)\longmapsto(W,U,p_l,k_r),
    !> \]
    !>
    !> with \(W\) dimensionless, \(U\) in J/m3, \(p_l\) in Pa and \(k_r\)
    !> dimensionless and normalised to one at the saturated unfrozen reference.
    !> Exact for the configured models up to their own evaluation. \(O(1)\), one
    !> phase split and a handful of model calls.
    !>
    !> Fails with `ok = .false.` and an untouched result for a slot outside the
    !> table, and when the conductivity model returns a value that is not
    !> positive and finite.
    subroutine eval_frozen_response(slot, temperature, pressure, response, ok)
        implicit none
        !> One-based index into the frozen closure's material table.
        integer(int32), intent(in) :: slot
        !> Temperature, C.
        real(real64), intent(in) :: temperature
        !> Total-water-equivalent pore pressure \(p_w\), Pa.
        real(real64), intent(in) :: pressure
        !> The response; written only when `ok`.
        type(type_mixed_frozen_response), intent(inout) :: response
        !> True when every model supplied a usable value.
        logical, intent(inout) :: ok

        type(type_mixed_frozen_response) :: local
        real(real64) :: suction, value
        integer(int32) :: which

        integer(int32) :: entry

        ok = .false.
        response_requests = response_requests + 1_int64

        entry = response_cache_entry(slot, temperature, pressure)
        if (entry > 0) then
            if (response_cache_filled(entry)) then
                if (response_cache_material(entry) == slot .and. &
                    response_cache_temperature(entry) == temperature .and. &
                    response_cache_pressure(entry) == pressure) then
                    response = response_cache_value(entry)
                    ok = .true.
                    return
                end if
            end if
        end if

        response_evaluations = response_evaluations + 1_int64
        call ensure_state(slot, temperature, pressure, which)
        if (which == 0) return

        call relative_conductivity(slot, temperature, pressure, local%relative_conductivity, ok)
        if (.not. ok) return

        call thermal_conductivity_value(slot, which, local%thermal_conductivity, ok)
        if (.not. ok) return

        call hydraulic_model%compute_total_storage(material_ids(slot), slot_state(which), local%stored_water)
        call thermal_model%calc_enthalpy_density(material_ids(slot), slot_state(which), local%enthalpy)
        call slot_state(which)%effective_suction%get(suction)
        local%liquid_pressure = -suction
        call slot_state(which)%water_content%get(local%liquid_content)
        call slot_state(which)%ice_content%get(local%ice_content)
        call slot_state(which)%vapor_content%get(local%vapor_content)
        call thermal_model%calc_density_water(slot_state(which), local%water_density)
        call thermal_model%calc_density_ice(slot_state(which), local%ice_density)

        ! The compressive term as the difference against the phase sum. Reading
        ! it back out of W rather than re-forming S_s/(rho g) max(0,p-p_sat)
        ! keeps the branch flag tied to the storage the residual differences,
        ! including the saturation pressure the phase model decided on.
        local%compressive_storage = local%stored_water - (local%liquid_content + local%vapor_content)
        if (local%water_density > 0.0_real64) then
            local%compressive_storage = local%compressive_storage &
                                        - (local%ice_density/local%water_density)*local%ice_content
        end if
        local%ice_present = (local%ice_content > 0.0_real64)
        local%compressive_active = (local%compressive_storage > 0.0_real64)

        value = local%stored_water + local%enthalpy + local%liquid_pressure + local%thermal_conductivity
        if (.not. ieee_is_finite(value)) then
            ok = .false.
            return
        end if
        response = local

        ! Stored only on success, so a refused key is refused again every time.
        if (entry > 0) then
            response_cache_value(entry) = local
            response_cache_material(entry) = slot
            response_cache_temperature(entry) = temperature
            response_cache_pressure(entry) = pressure
            response_cache_filled(entry) = .true.
        end if
    end subroutine eval_frozen_response

    !> \(\lambda\) of one material at a published state, W/(m K).
    !>
    !> Read through the governing model's own thermal diffusion tensor, so the
    !> phase mixing rule is the model's. The tensor has to be isotropic: the
    !> coupled column's conduction term carries a scalar, and reporting one
    !> entry of an anisotropic tensor as that scalar would answer a different
    !> problem.
    subroutine thermal_conductivity_value(slot, which, value, ok)
        implicit none
        integer(int32), intent(in) :: slot
        !> The cached state slot, as `ensure_state` returned it.
        integer(int32), intent(in) :: which
        real(real64), intent(inout) :: value
        logical, intent(inout) :: ok

        real(real64) :: tensor(2, 2), scale

        ok = .false.
        tensor(:, :) = 0.0_real64
        call thermal_model%compute_diffusion_term(material_ids(slot), slot_state(which), tensor)

        scale = max(abs(tensor(1, 1)), abs(tensor(2, 2)))
        if (.not. ieee_is_finite(scale)) return
        if (scale <= 0.0_real64) return
        if (abs(tensor(1, 1) - tensor(2, 2)) > ISOTROPY_TOLERANCE*scale) return
        if (max(abs(tensor(1, 2)), abs(tensor(2, 1))) > ISOTROPY_TOLERANCE*scale) return

        value = tensor(1, 1)
        ok = .true.
    end subroutine thermal_conductivity_value

    !> The tangents of that response in both unknowns.
    !>
    !> Liquid-pressure tangents are published by the phase model. Storage,
    !> enthalpy and log-conductivity tangents are differences of the actual
    !> response, including dependencies omitted by Picard capacities.
    !>
    !> \[
    !>   \frac{\partial F}{\partial y}\approx
    !>   \begin{cases}
    !>     \dfrac{F(y+\delta)-F(y-\delta)}{2\delta} & \text{one branch}\\[2ex]
    !>     \dfrac{F(y\pm\delta)-F(y)}{\pm\delta} & \text{branch change on the other side}
    !>   \end{cases}
    !> \]
    !>
    !> The branch of a state is which of the compressive storage and the ice are
    !> present. Where an offset crosses the saturation pressure or the freezing
    !> onset, the central difference is the mean of two different one-sided
    !> slopes and is a tangent of neither; the difference is then taken on the
    !> side that stayed on the point's own branch, and that is reported in
    !> `pressure_one_sided` / `temperature_one_sided`. When both sides leave it,
    !> the central value is returned and `branch_straddled` says so rather than
    !> the caller receiving a number with no derivative behind it.
    !> \(O(1)\), at most five phase splits.
    !>
    !> Fails with `ok = .false.` and untouched results whenever
    !> `eval_frozen_response` does at the point or at either offset.
    subroutine eval_frozen_tangents(slot, temperature, pressure, tangents, ok)
        implicit none
        !> One-based index into the frozen closure's material table.
        integer(int32), intent(in) :: slot
        !> Temperature, C.
        real(real64), intent(in) :: temperature
        !> Total-water-equivalent pore pressure \(p_w\), Pa.
        real(real64), intent(in) :: pressure
        !> The tangents; written only when `ok`.
        type(type_mixed_frozen_tangents), intent(inout) :: tangents
        !> True when every model supplied a usable value.
        logical, intent(inout) :: ok

        type(type_mixed_frozen_tangents) :: local
        type(type_mixed_frozen_response) :: base, ahead, behind
        real(real64) :: dsuction_dp, dsuction_dT
        integer(int32) :: which, entry
        logical :: straddled, ahead_same, behind_same

        ok = .false.
        tangent_requests = tangent_requests + 1_int64

        ! A hit returns the value this same key already produced, so the caller
        ! cannot tell a hit from a miss by the number it receives.
        entry = tangent_cache_entry(slot, temperature, pressure)
        if (entry > 0) then
            if (tangent_cache_filled(entry)) then
                if (tangent_cache_material(entry) == slot .and. &
                    tangent_cache_temperature(entry) == temperature .and. &
                    tangent_cache_pressure(entry) == pressure) then
                    tangents = tangent_cache_value(entry)
                    ok = .true.
                    return
                end if
            end if
        end if

        tangent_evaluations = tangent_evaluations + 1_int64
        call ensure_state(slot, temperature, pressure, which)
        if (which == 0) return

        ! The liquid driver has explicit phase-closure derivatives. Storage and
        ! enthalpy are differentiated below from the actual residual quantities:
        ! the governing Picard capacities omit some dependencies of those values.
        call slot_state(which)%dSeff_dP%get(dsuction_dp)
        call slot_state(which)%dSeff_dT%get(dsuction_dT)
        local%dliquid_dpressure = -dsuction_dp
        local%dliquid_dtemperature = -dsuction_dT

        ! The point's own branch, which the offsets below are compared against.
        call eval_frozen_response(slot, temperature, pressure, base, ok)
        if (.not. ok) return

        ! Differenced, because the models publish no derivative of the enthalpy
        ! in pressure and none of the conductivity in either unknown.
        call eval_frozen_response(slot, temperature, pressure + PRESSURE_STEP, ahead, ok)
        if (.not. ok) return
        call eval_frozen_response(slot, temperature, pressure - PRESSURE_STEP, behind, ok)
        if (.not. ok) return
        ahead_same = same_branch(base, ahead)
        behind_same = same_branch(base, behind)
        call difference_on_branch(base%stored_water, ahead%stored_water, behind%stored_water, &
                                  PRESSURE_STEP, ahead_same, behind_same, &
                                  local%dwater_dpressure, local%pressure_one_sided, straddled)
        call difference_on_branch(base%enthalpy, ahead%enthalpy, behind%enthalpy, &
                                  PRESSURE_STEP, ahead_same, behind_same, &
                                  local%denthalpy_dpressure, local%pressure_one_sided, straddled)
        call difference_on_branch(log(base%relative_conductivity), log(ahead%relative_conductivity), &
                                  log(behind%relative_conductivity), &
                                  PRESSURE_STEP, ahead_same, behind_same, &
                                  local%dlog_conductivity_dpressure, local%pressure_one_sided, straddled)
        call difference_on_branch(base%thermal_conductivity, ahead%thermal_conductivity, behind%thermal_conductivity, &
                                  PRESSURE_STEP, ahead_same, behind_same, &
                                  local%dthermal_dpressure, local%pressure_one_sided, straddled)
        local%branch_straddled = local%branch_straddled .or. straddled

        call eval_frozen_response(slot, temperature + TEMPERATURE_STEP, pressure, ahead, ok)
        if (.not. ok) return
        call eval_frozen_response(slot, temperature - TEMPERATURE_STEP, pressure, behind, ok)
        if (.not. ok) return
        ahead_same = same_branch(base, ahead)
        behind_same = same_branch(base, behind)
        call difference_on_branch(base%stored_water, ahead%stored_water, behind%stored_water, &
                                  TEMPERATURE_STEP, ahead_same, behind_same, &
                                  local%dwater_dtemperature, local%temperature_one_sided, straddled)
        call difference_on_branch(base%enthalpy, ahead%enthalpy, behind%enthalpy, &
                                  TEMPERATURE_STEP, ahead_same, behind_same, &
                                  local%denthalpy_dtemperature, local%temperature_one_sided, straddled)
        call difference_on_branch(log(base%relative_conductivity), log(ahead%relative_conductivity), &
                                  log(behind%relative_conductivity), &
                                  TEMPERATURE_STEP, ahead_same, behind_same, &
                                  local%dlog_conductivity_dtemperature, local%temperature_one_sided, straddled)
        call difference_on_branch(base%thermal_conductivity, ahead%thermal_conductivity, behind%thermal_conductivity, &
                                  TEMPERATURE_STEP, ahead_same, behind_same, &
                                  local%dthermal_dtemperature, local%temperature_one_sided, straddled)
        local%branch_straddled = local%branch_straddled .or. straddled

        tangents = local
        ok = .true.

        ! Stored only now, so that a key the closure refused is never answered
        ! from the memo and keeps reporting its fault.
        if (entry > 0) then
            tangent_cache_value(entry) = local
            tangent_cache_material(entry) = slot
            tangent_cache_temperature(entry) = temperature
            tangent_cache_pressure(entry) = pressure
            tangent_cache_filled(entry) = .true.
        end if
    end subroutine eval_frozen_tangents

    !> Whether two evaluated states are on the same branch of the closure.
    !>
    !> The branch is the pair of switches the closure actually has: whether ice
    !> is present, and whether the compressive storage is open. Both are exact
    !> in the sense that the value, not a tolerance, decides them.
    pure logical function same_branch(left, right)
        implicit none
        type(type_mixed_frozen_response), intent(in) :: left, right

        same_branch = (left%ice_present .eqv. right%ice_present) .and. &
                      (left%compressive_active .eqv. right%compressive_active)
    end function same_branch

    !> One difference quotient, central where the branch allows it.
    subroutine difference_on_branch(base, ahead, behind, step, ahead_same, behind_same, &
                                    derivative, one_sided, straddled)
        implicit none
        !> The quantity at the point and at the two offsets.
        real(real64), intent(in) :: base, ahead, behind
        !> The offset, positive.
        real(real64), intent(in) :: step
        !> Whether each offset stayed on the point's branch.
        logical, intent(in) :: ahead_same, behind_same
        !> The reported derivative.
        real(real64), intent(inout) :: derivative
        !> Set when a one-sided difference was used; never cleared, so one call
        !> per direction accumulates over the quantities of that direction.
        logical, intent(inout) :: one_sided
        !> Whether both offsets left the branch.
        logical, intent(inout) :: straddled

        straddled = .false.
        if (ahead_same .and. behind_same) then
            derivative = (ahead - behind)/(2.0_real64*step)
        else if (ahead_same) then
            derivative = (ahead - base)/step
            one_sided = .true.
        else if (behind_same) then
            derivative = (base - behind)/step
            one_sided = .true.
        else
            derivative = (ahead - behind)/(2.0_real64*step)
            straddled = .true.
        end if
    end subroutine difference_on_branch

    !> How much of the pressure column of one cell's Jacobian survives.
    !>
    !> The unknown \(p_w\) reaches the coupled system through exactly two
    !> quantities: the Darcy driver, through \(\partial p_l/\partial p_w\), and
    !> the storage row, through \(\partial W/\partial p_w\). Where the phase
    !> split is on the freezing branch the first is zero because
    !> \(\psi_{\mathrm{eff}}=\psi_{\mathrm{cryo}}(T)\) does not contain
    !> \(p_w\); where the retention curve is flat and no compressive branch is
    !> open the second is zero as well. When both vanish the cell's pressure
    !> column is zero and the monolithic matrix is singular, whatever the
    !> nonlinear solver reports.
    !>
    !> \[
    !>   m=\max\left(\left|\frac{\partial p_l}{\partial p_w}\right|,\;
    !>   \frac{p_{\mathrm{ref}}}{W_{\mathrm{ref}}}
    !>   \left|\frac{\partial W}{\partial p_w}\right|\right)
    !> \]
    !>
    !> is that column's dimensionless margin: both terms are the response of a
    !> reference-sized quantity to a reference-sized pressure change, so they
    !> are comparable to each other and to one. \(O(1)\), one tangent
    !> evaluation.
    !>
    !> Fails with `ok = .false.` and untouched results wherever
    !> `eval_frozen_tangents` does.
    subroutine frozen_pressure_column(slot, temperature, pressure, pressure_reference, storage_reference, &
                                      margin, liquid_sensitivity, storage_capacity, ok)
        implicit none
        !> One-based index into the frozen closure's material table.
        integer(int32), intent(in) :: slot
        !> Temperature, C.
        real(real64), intent(in) :: temperature
        !> Total-water-equivalent pore pressure \(p_w\), Pa.
        real(real64), intent(in) :: pressure
        !> \(p_{\mathrm{ref}}\), the pressure change the run resolves, Pa; positive.
        real(real64), intent(in) :: pressure_reference
        !> \(W_{\mathrm{ref}}\), the stored water the run resolves; positive.
        real(real64), intent(in) :: storage_reference
        !> The dimensionless margin; written only when `ok`.
        real(real64), intent(inout) :: margin
        !> \(\partial p_l/\partial p_w\); written only when `ok`.
        real(real64), intent(inout) :: liquid_sensitivity
        !> \(\partial W/\partial p_w\), 1/Pa; written only when `ok`.
        real(real64), intent(inout) :: storage_capacity
        !> True when the closure answered at this state.
        logical, intent(inout) :: ok

        type(type_mixed_frozen_tangents) :: tangents

        ok = .false.
        if (pressure_reference <= 0.0_real64 .or. storage_reference <= 0.0_real64) return

        call eval_frozen_tangents(slot, temperature, pressure, tangents, ok)
        if (.not. ok) return

        liquid_sensitivity = tangents%dliquid_dpressure
        storage_capacity = tangents%dwater_dpressure
        margin = max(abs(tangents%dliquid_dpressure), &
                     (pressure_reference/storage_reference)*abs(tangents%dwater_dpressure))
    end subroutine frozen_pressure_column

    !> \(k_r=K_{fl h}(T,p_w)/K_{\mathrm{sat}}\) of one cell.
    subroutine relative_conductivity(slot, temperature, pressure, relative, ok)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        real(real64), intent(inout) :: relative
        logical, intent(inout) :: ok

        real(real64) :: value

        call conductivity_value(slot, temperature, pressure, value, ok)
        if (.not. ok) return

        value = value/saturated_conductivity(slot)
        ok = ieee_is_finite(value) .and. value > 0.0_real64
        if (ok) relative = value
    end subroutine relative_conductivity

    !> \(K_{fl h}\) of one material at a state, m/s.
    !>
    !> The conductivity is read through the governing model's own Darcy tensor,
    !> \(D_{HH}=K_{fl h}/(\rho_w g)\), so the impedance and viscosity factors of
    !> the configured mode are applied by the model rather than reassembled. The
    !> model evaluates its retention argument at the head of the effective
    !> suction the state carries, which is where the cryosuction of this contract
    !> reaches the mobility.
    subroutine conductivity_value(slot, temperature, pressure, value, ok)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        real(real64), intent(inout) :: value
        logical, intent(inout) :: ok

        real(real64) :: tensor(2, 2), density
        integer(int32) :: which

        ok = .false.
        call ensure_state(slot, temperature, pressure, which)
        if (which == 0) return

        tensor(:, :) = 0.0_real64
        call hydraulic_model%compute_diffusion_term(material_ids(slot), slot_state(which), tensor)
        call thermal_model%calc_density_water(slot_state(which), density)
        ! compute_diffusion_term forms K/(rho_w g) with the same two quantities,
        ! so multiplying them back recovers the model's own conductivity without
        ! restating the model.
        value = tensor(1, 1)*density*gravity_acceleration
        ok = ieee_is_finite(value) .and. value > 0.0_real64
    end subroutine conductivity_value

    !> Whether the freezing model is defined at this pair of unknowns.
    !>
    !> The model's own precondition, tested before it is entered rather than
    !> inside it: both unknowns finite, an absolute temperature it can divide
    !> by, and a relative-humidity exponent `update_water_phases`
    !> (phase_systems.F90) can pass to `exp()` without overflowing -- the same
    !> \(|pM_w/(\rho_{std}R T_K)| \le 700\) bound that routine enforces with
    !> `error stop`, restated here so a Newton iterate that would trip it is
    !> refused as an ordinary rejected trial instead of ending the run. Nothing
    !> here is a physical admissibility criterion -- a state this accepts may
    !> still be one no soil is ever in -- it is the arithmetic domain of the
    !> constitutive evaluation and nothing more.
    pure logical function state_is_admissible(temperature, pressure)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure

        real(real64) :: temperature_K, rh_exponent

        ! Fortran does not guarantee short-circuit evaluation, so each
        ! comparison below is reached only after the ones before it pass.
        state_is_admissible = .false.
        if (.not. ieee_is_finite(temperature)) return
        if (.not. ieee_is_finite(pressure)) return
        temperature_K = temperature + celsius_to_kelvin
        if (.not. (temperature_K > tiny(1.0_real64))) return
        rh_exponent = (pressure * rh_molar_mass_water)/(rh_reference_water_density * rh_gas_constant * temperature_K)
        if (.not. ieee_is_finite(rh_exponent)) return
        state_is_admissible = abs(rh_exponent) <= 700.0_real64
    end function state_is_admissible

    !> The slot holding the published state of this point, publishing it only if
    !> no slot holds it already.
    !>
    !> The key is exact equality of the material and both unknowns, so a state
    !> asked for twice is answered from the first answer bit for bit. That is
    !> what lets an accepted water history be re-differenced later and give the
    !> same number.
    subroutine ensure_state(slot, temperature, pressure, which)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        !> Slot of the cached state, zero when the material is not in the table.
        integer(int32), intent(inout) :: which

        logical :: ok
        integer(int32) :: k

        which = 0
        if (slot < 1 .or. slot > num_frozen_laws) return

        do k = 1, STATE_SLOTS
            if (.not. slot_filled(k)) cycle
            if (slot_material(k) /= slot) cycle
            if (slot_temperature(k) /= temperature) cycle
            if (slot_pressure(k) /= pressure) cycle
            which = k
            return
        end do

        which = slot_next
        slot_next = merge(1, slot_next + 1, slot_next >= STATE_SLOTS)

        ok = .false.
        call publish_frozen_state(slot, temperature, pressure, slot_state(which), ok)
        if (.not. ok) then
            slot_filled(which) = .false.
            which = 0
            return
        end if

        slot_material(which) = slot
        slot_temperature(which) = temperature
        slot_pressure(which) = pressure
        slot_filled(which) = .true.
    end subroutine ensure_state

    !> Choose whether tangents are memoised, and on how many entries.
    !>
    !> The memo returns, for a repeated key, the value that key already
    !> produced, so enabling it cannot move a converged answer: the map it
    !> remembers is a deterministic function of \((\text{slot},T,p_w)\) alone.
    !> It is a cost control, not a numerical one.
    !>
    !> \(O(n)\) memory in `entries`, \(O(1)\) per lookup. Entries are rounded up
    !> to a power of two so that the index is a mask rather than a division; a
    !> table too small for the working set costs recomputation and nothing else.
    !>
    !> Cannot fail. A non-positive `entries` asks for the default size.
    subroutine configure_frozen_tangent_cache(enabled, entries)
        implicit none
        !> Whether a repeated key is answered from the memo.
        logical, intent(in) :: enabled
        !> How many keys the memo may hold; rounded up to a power of two. Not
        !> positive asks for the default. Size it to a rank's quadrature points.
        integer(int32), intent(in) :: entries

        integer(int32) :: wanted, size_now

        ! Released first, so that a second call resizes rather than fails.
        if (allocated(tangent_cache_value)) then
            deallocate (tangent_cache_value, tangent_cache_temperature, &
                        tangent_cache_pressure, tangent_cache_material, tangent_cache_filled)
        end if
        if (allocated(response_cache_value)) then
            deallocate (response_cache_value, response_cache_temperature, &
                        response_cache_pressure, response_cache_material, response_cache_filled)
        end if
        tangent_cache_mask = 0

        tangent_cache_enabled = enabled
        if (.not. enabled) return

        wanted = entries
        if (wanted <= 0) wanted = TANGENT_CACHE_DEFAULT_ENTRIES

        size_now = 1
        do while (size_now < wanted .and. size_now < huge(size_now)/2)
            size_now = size_now*2
        end do

        allocate (tangent_cache_value(size_now))
        allocate (tangent_cache_temperature(size_now), source=0.0_real64)
        allocate (tangent_cache_pressure(size_now), source=0.0_real64)
        allocate (tangent_cache_material(size_now), source=0_int32)
        allocate (tangent_cache_filled(size_now), source=.false.)

        allocate (response_cache_value(size_now))
        allocate (response_cache_temperature(size_now), source=0.0_real64)
        allocate (response_cache_pressure(size_now), source=0.0_real64)
        allocate (response_cache_material(size_now), source=0_int32)
        allocate (response_cache_filled(size_now), source=.false.)

        tangent_cache_mask = size_now - 1
    end subroutine configure_frozen_tangent_cache

    !> How much constitutive work was asked for and how much was done.
    !>
    !> Counted since the closure was installed, on the calling rank only.
    !> `requests` counts calls to `eval_frozen_tangents` and `evaluations`
    !> counts those that actually differenced a tangent, so their ratio is the
    !> repetition the memo removes. \(O(1)\). Cannot fail.
    subroutine frozen_work_counters(requests, evaluations, asked, responses, publications)
        implicit none
        !> Calls to `eval_frozen_tangents`.
        integer(int64), intent(inout) :: requests
        !> Those of them that computed a tangent rather than recalling one.
        integer(int64), intent(inout) :: evaluations
        !> Calls to `eval_frozen_response`.
        integer(int64), intent(inout) :: asked
        !> Those of them that computed a response rather than recalling one.
        integer(int64), intent(inout) :: responses
        !> Phase splits actually solved, which `ensure_state` shields.
        integer(int64), intent(inout) :: publications

        requests = tangent_requests
        evaluations = tangent_evaluations
        asked = response_requests
        responses = response_evaluations
        publications = state_publications
    end subroutine frozen_work_counters

    !> Drop every memoised value, keeping the configured size.
    subroutine forget_frozen_tangents()
        implicit none

        if (allocated(tangent_cache_filled)) tangent_cache_filled(:) = .false.
        if (allocated(response_cache_filled)) response_cache_filled(:) = .false.
    end subroutine forget_frozen_tangents

    !> Where a key lives in the memo, or zero when there is no memo.
    !>
    !> Mixing is by exclusive-or and shift only: neighbouring quadrature points
    !> differ in the low bits of an IEEE double, which have to reach the index,
    !> and a product would have to overflow to do the same work.
    function tangent_cache_entry(slot, temperature, pressure) result(entry)
        implicit none
        !> One-based index into the frozen closure's material table.
        integer(int32), intent(in) :: slot
        !> Temperature, C.
        real(real64), intent(in) :: temperature
        !> Total-water-equivalent pore pressure \(p_w\), Pa.
        real(real64), intent(in) :: pressure
        !> The entry to use, or zero when the memo is off.
        integer(int32) :: entry

        integer(int64) :: key

        entry = 0
        if (.not. tangent_cache_enabled) return
        if (.not. allocated(tangent_cache_filled)) return

        entry = mixed_key_entry(slot, temperature, pressure)
    end function tangent_cache_entry

    !> Where a key lives in the response memo, or zero when there is no memo.
    function response_cache_entry(slot, temperature, pressure) result(entry)
        implicit none
        !> One-based index into the frozen closure's material table.
        integer(int32), intent(in) :: slot
        !> Temperature, C.
        real(real64), intent(in) :: temperature
        !> Total-water-equivalent pore pressure \(p_w\), Pa.
        real(real64), intent(in) :: pressure
        !> The entry to use, or zero when the memo is off.
        integer(int32) :: entry

        entry = 0
        if (.not. tangent_cache_enabled) return
        if (.not. allocated(response_cache_filled)) return

        entry = mixed_key_entry(slot, temperature, pressure)
    end function response_cache_entry

    !> The table index one key mixes down to.
    !>
    !> Mixing is by exclusive-or and shift only: neighbouring quadrature points
    !> differ in the low bits of an IEEE double, which have to reach the index,
    !> and a product would have to overflow to do the same work.
    function mixed_key_entry(slot, temperature, pressure) result(entry)
        implicit none
        !> One-based index into the frozen closure's material table.
        integer(int32), intent(in) :: slot
        !> Temperature, C.
        real(real64), intent(in) :: temperature
        !> Total-water-equivalent pore pressure \(p_w\), Pa.
        real(real64), intent(in) :: pressure
        !> The entry to use, one-based.
        integer(int32) :: entry

        integer(int64) :: key

        key = ieor(transfer(temperature, 0_int64), ishft(transfer(pressure, 0_int64), 1))
        key = ieor(key, int(slot, int64))
        key = ieor(key, ishft(key, -33))
        key = ieor(key, ishft(key, -17))
        key = ieor(key, ishft(key, -8))
        entry = int(iand(key, int(tangent_cache_mask, int64)), int32) + 1
    end function mixed_key_entry

end module app_mixed_hydraulic_frozen
