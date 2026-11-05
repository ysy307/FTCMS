!> The frozen constitutive contract every water and heat kernel evaluates.
!>
!> One map, evaluated cell-locally,
!>
!> \[
!>   (T,p_w,\text{material})\longmapsto
!>   \bigl(W,\;U,\;p_l,\;k_r\bigr),
!> \]
!>
!> together with its tangents in both unknowns. Every value is produced by the
!> repository's own governing models -- `type_hydraulic` and `type_thermal`,
!> built from `Basic.json` -- so every discretisation evaluates the same
!> freezing closure rather than a second statement of it. Nothing here restates a retention curve, a freezing curve, an impedance
!> law or a storage formula.
!>
!> ## The two pressures
!>
!> The solved unknown is \(p_w\), the **total-water-equivalent** pore
!> pressure: the closure defines the total stored water as the unfrozen
!> equivalent content at that pressure, and it is \(p_w\) the mass balance
!> transports. It is not the pressure of the liquid phase.
!>
!> The liquid phase sees \(p_l=-\psi_{\mathrm{eff}}\), where
!> \(\psi_{\mathrm{eff}}\) is the effective suction the phase split publishes: a
!> smooth maximum of the matric suction of \(p_w\) and the freezing suction of
!> the generalized Clapeyron relation, limited by the pore volume. Freezing
!> raises that suction, and \(\nabla p_l\) is what draws water toward a front --
!> the cryosuction a freezing soil shows.
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
!> This is the full frozen quantity: the ice term through its density ratio,
!> the vapour term when the run carries it, and the compressive storage.
!>
!> ## Tangents
!>
!> Every tangent is the governing models' own, evaluated at the same published
!> state as the value: the phase split publishes the analytic tangents of the
!> contents and of \(\psi_{\mathrm{eff}}\), and the models compose them through
!> the closed-form laws of `constitutive_laws`, whose value and derivative are
!> generated from one SymPy statement. Nothing is differenced, so a tangent costs
!> one phase split, the one its value needs. Each tangent is the derivative on
!> the branch the state is on (ice present or not, compressive storage open or
!> not), which at a branch change is the one-sided derivative from that side.
!>
!> ## Ownership
!>
!> The governing models are module state for the same reason the isothermal
!> tables are: a PetscDS pointwise function is reached through a C function
!> pointer and carries no typed context. One solver per rank owns them, and
!> installing a second while the first holds them is refused.
module app_ftcms_closure
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use :: core_physical_constants, only:gravity_acceleration
    use :: module_core, only:type_state, type_config_hcf
    use :: module_input, only:type_input, input_translator
    use :: constitutive_constants, only:celsius_to_kelvin, tangent_pressure, tangent_temperature, &
        rh_gas_constant => universal_gas_constant, rh_molar_mass_water => molar_mass_water, &
        rh_reference_water_density => reference_water_density
    use :: physics_governing_thermal, only:type_thermal
    use :: physics_governing_hydraulic, only:type_hydraulic
    implicit none
    private

    public :: type_closure_response
    public :: type_closure_tangents
    public :: install_closure
    public :: closure_saturated_conductivity
    public :: clear_closure
    public :: eval_closure_response
    public :: eval_closure_tangents
    public :: eval_liquid_slopes
    public :: eval_closure_storage
    public :: configure_closure_cache
    public :: closure_work_counters

    !> How many materials one run may carry frozen laws for.
    integer(int32), parameter :: MAX_CLOSURE_MATERIALS = 16

    !> Installation probe pressure, Pa; not a conductivity normalisation state.
    real(real64), parameter :: REFERENCE_PRESSURE = 0.0d0
    !> Installation probe temperature, C. No unfrozen-state guarantee is assumed.
    real(real64), parameter :: REFERENCE_TEMPERATURE = 20.0d0



    !> Relative departure from isotropy a thermal conductivity tensor may carry.
    real(real64), parameter :: ISOTROPY_TOLERANCE = 1.0d-12

    !> What one cell's state is, at one \((T,p_w)\).
    type :: type_closure_response
        !> \(W\), total stored water per unit bulk volume, dimensionless.
        real(real64) :: stored_water = 0.0d0
        !> \(U\), enthalpy density, J/m3.
        real(real64) :: enthalpy = 0.0d0
        !> \(p_l=-\psi_{\mathrm{eff}}\), the pressure the liquid flux follows, Pa.
        real(real64) :: liquid_pressure = 0.0d0
        !> \(k_r\), relative conductivity, dimensionless and positive.
        real(real64) :: relative_conductivity = 0.0d0
        !> \(\lambda\), bulk thermal conductivity, W/(m K).
        !>
        !> The scalar of an isotropic tensor. `install_closure` refuses a
        !> material whose thermal conductivity tensor is not isotropic rather
        !> than silently reporting one of its entries.
        real(real64) :: thermal_conductivity = 0.0d0
        !> \(\theta_w\), liquid water content.
        real(real64) :: liquid_content = 0.0d0
        !> \(\theta_i\), ice content.
        real(real64) :: ice_content = 0.0d0
        !> \(\theta_v\), vapour content.
        real(real64) :: vapor_content = 0.0d0
        !> \(\rho_w\), liquid density at this state, kg/m3.
        !>
        !> Reported because \(W\) counts ice through \(\rho_i/\rho_w\): the
        !> conversion is part of the contract and has to be auditable from
        !> outside it.
        real(real64) :: water_density = 0.0d0
        !> Volumetric heat capacity of the liquid, \(\rho_w c_w\), J/(m3 K): what
        !> a unit of liquid flux carries per kelvin.
        real(real64) :: liquid_heat_capacity = 0.0d0
        !> \(\partial p_l/\partial p_w\) and \(\partial p_l/\partial T\), from the
        !> same phase split: analytic, so a Darcy driver needs no difference.
        real(real64) :: liquid_dpressure = 0.0d0
        real(real64) :: liquid_dtemperature = 0.0d0
        !> \(\rho_i\), ice density at this state, kg/m3.
        real(real64) :: ice_density = 0.0d0
        !> \(S_c\), the compressive part of \(W\), dimensionless.
        !>
        !> Read as \(W\) minus the phase sum rather than from a second storage
        !> formula, so it is by construction the term the residual carries.
        real(real64) :: compressive_storage = 0.0d0
        !> Whether ice is present, \(\theta_i>0\).
        logical :: ice_present = .false.
        !> Whether the compressive branch \(\max(0,p_w-p_{\mathrm{sat}})\) is open.
        logical :: compressive_active = .false.
    end type type_closure_response

    !> The tangents of that state in both unknowns.
    type :: type_closure_tangents
        !> \(\partial W/\partial p_w\), 1/Pa.
        real(real64) :: dwater_dpressure = 0.0d0
        !> \(\partial W/\partial T\), 1/K.
        real(real64) :: dwater_dtemperature = 0.0d0
        !> \(\partial U/\partial T\), the apparent heat capacity, J/(m3 K).
        real(real64) :: denthalpy_dtemperature = 0.0d0
        !> \(\partial U/\partial p_w\), J/(m3 Pa).
        real(real64) :: denthalpy_dpressure = 0.0d0
        !> \(\partial p_l/\partial p_w\), dimensionless.
        real(real64) :: dliquid_dpressure = 0.0d0
        !> \(\partial p_l/\partial T\), Pa/K.
        real(real64) :: dliquid_dtemperature = 0.0d0
        !> \(\partial\ln k_r/\partial p_w\), 1/Pa.
        real(real64) :: dlog_conductivity_dpressure = 0.0d0
        !> \(\partial\ln k_r/\partial T\), 1/K.
        real(real64) :: dlog_conductivity_dtemperature = 0.0d0
        !> \(\partial\lambda/\partial p_w\), W/(m K Pa).
        real(real64) :: dthermal_dpressure = 0.0d0
        !> \(\partial\lambda/\partial T\), W/(m K K).
        real(real64) :: dthermal_dtemperature = 0.0d0
    end type type_closure_tangents

    ! The governing models of the run, and the per-slot data a kernel reaches
    ! its own material with.
    type(type_hydraulic), save :: hydraulic_model
    type(type_thermal), save :: thermal_model
    integer(int32), save :: material_ids(MAX_CLOSURE_MATERIALS) = 0
    real(real64), save :: material_porosity(MAX_CLOSURE_MATERIALS) = 0.0d0
    real(real64), save :: saturated_conductivity(MAX_CLOSURE_MATERIALS) = 1.0d0
    integer(int32), save :: num_closure_materials = 0

    ! A small exactly-keyed cache of published states. One evaluation asks for
    ! the state of its key several times (the response, its conductivity, the
    ! tangents), and a storage kernel visits a point's current state and the
    ! accepted states its BDF rate differences one after another.
    integer(int32), parameter :: STATE_SLOTS = 5
    type(type_state), save :: slot_state(STATE_SLOTS)
    real(real64), save :: slot_temperature(STATE_SLOTS) = 0.0d0
    real(real64), save :: slot_pressure(STATE_SLOTS) = 0.0d0
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
    ! Assembling the Jacobian calls eval_closure_tangents once per
    ! registered block per quadrature point -- twelve times -- with bit-identical
    ! arguments, and the map is a deterministic function of its key alone, so
    ! every call after the first can only recompute the number the first one
    ! already produced. PETSc integrates one field block over every cell before
    ! it begins the next, so a point is revisited only after a whole sweep of the
    ! rank's cells: the memo is sized for that working set. STATE_SLOTS above
    ! holds a few states and cannot span blocks, which is why it does not remove
    ! this repetition by itself.
    !
    ! Four-way set associative: a key may sit in any way of its set, so two
    ! keys that mix to the same set do not evict each other on every visit, as
    ! the accepted states of a step did when the table was direct-mapped. A
    ! full set overwrites one way in turn; that costs a recomputation and
    ! never a wrong answer. Only successful evaluations are stored, so a
    ! point the closure refuses is refused again on every visit and its fault is
    ! reported by the caller exactly as it was before.
    !
    ! Like slot_state above, this is module state without a lock, so it inherits
    ! the same single-threaded-assembly assumption the published states already
    ! make. It adds no sharing that was not there.
    integer(int32), parameter :: TANGENT_CACHE_DEFAULT_ENTRIES = 65536
    integer(int32), parameter :: MEMO_WAYS = 4
    ! The way a full set gives up next, turned on every replacement.
    integer(int32), save :: memo_victim = 0
    logical, save :: tangent_cache_enabled = .true.
    integer(int32), save :: tangent_cache_mask = 0
    type(type_closure_tangents), allocatable, save :: tangent_cache_value(:)
    real(real64), allocatable, save :: tangent_cache_temperature(:)
    real(real64), allocatable, save :: tangent_cache_pressure(:)
    integer(int32), allocatable, save :: tangent_cache_material(:)
    logical, allocatable, save :: tangent_cache_filled(:)

    ! The same memo for the response, which the residual repeats for the same
    ! reason: its five registered kernels each ask for one at the same key. It
    ! is a separate table because the residual is also evaluated at line-search
    ! trials, where no tangent is wanted and none is computed.
    type(type_closure_response), allocatable, save :: response_cache_value(:)
    real(real64), allocatable, save :: response_cache_temperature(:)
    real(real64), allocatable, save :: response_cache_pressure(:)
    integer(int32), allocatable, save :: response_cache_material(:)
    logical, allocatable, save :: response_cache_filled(:)

    ! The stored water and enthalpy of the accepted states a BDF rate
    ! differences. Every residual evaluation of a step revisits the same
    ! accepted states, while its current states change;
    ! in a table of their own they cannot be evicted by those, so they are
    ! computed about once per step instead of once per evaluation.
    real(real64), allocatable, save :: storage_cache_water(:)
    real(real64), allocatable, save :: storage_cache_enthalpy(:)
    real(real64), allocatable, save :: storage_cache_temperature(:)
    real(real64), allocatable, save :: storage_cache_pressure(:)
    integer(int32), allocatable, save :: storage_cache_material(:)
    logical, allocatable, save :: storage_cache_filled(:)

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
    subroutine install_closure(input, ids, porosities)
        implicit none
        !> The project input the governing models are built from.
        type(type_input), intent(in) :: input
        !> The "Cell Sets" value of each material, in slot order.
        integer(int32), intent(in) :: ids(:)
        !> The porosity of each material, in slot order.
        real(real64), intent(in) :: porosities(:)

        integer(int32) :: m, count, material
        type(type_config_hcf) :: conductivity_function
        type(type_closure_response) :: reference
        logical :: ok

        if (num_closure_materials /= 0) then
            error stop "A second solver cannot share the closure of a live one"
        end if
        count = int(size(ids), int32)
        if (count > MAX_CLOSURE_MATERIALS) then
            error stop "More materials than the closure can hold"
        end if
        if (count < 1) then
            error stop "The closure needs at least one material"
        end if
        if (size(porosities) /= count) error stop "The closure needs one porosity per material"

        call forget_closure_tangents()

        do m = 1, count
            if (ids(m) < 1) then
                error stop "A closure material needs a positive region identifier"
            end if
            if (porosities(m) <= 0.0d0 .or. porosities(m) >= 1.0d0) then
                error stop "A closure material needs a porosity strictly between zero and one"
            end if
            material_ids(m) = ids(m)
            material_porosity(m) = porosities(m)
            saturated_conductivity(m) = 1.0d0
        end do

        call thermal_model%initialize(input, ids)
        call hydraulic_model%initialize(input, ids)
        num_closure_materials = count

        do m = 1, count
            material = material_index(input, ids(m))
            call conductivity_function%reset()
            call input_translator%execute(input, material, conductivity_function)
            saturated_conductivity(m) = conductivity_function%k_sat
            if (.not. ieee_is_finite(saturated_conductivity(m))) then
                error stop "Nonfinite saturated conductivity in the closure"
            end if
            if (saturated_conductivity(m) <= 0.0d0) then
                error stop "Nonpositive saturated conductivity in the closure"
            end if
        end do

        call eval_closure_response(1_int32, REFERENCE_TEMPERATURE, REFERENCE_PRESSURE, reference, ok)
        if (.not. ok) then
            num_closure_materials = 0
            error stop "The closure could not evaluate its own reference state"
        end if
    end subroutine install_closure

    !> The index in `input%basic%materials` of the material whose id is `id`.
    !>
    !> A single material serves every region, as the governing models allow.
    !> Fails by `error stop` when several are defined and none matches.
    function material_index(input, id) result(index)
        implicit none
        type(type_input), intent(in) :: input
        integer(int32), intent(in) :: id
        integer(int32) :: index

        do index = 1, input%basic%num_materials
            if (input%basic%materials(index)%id == id) return
        end do
        if (input%basic%num_materials == 1) then
            index = 1
            return
        end if
        error stop "No material of Basic.json carries the region identifier of a mesh cell"
    end function material_index

    !> The saturated hydraulic conductivity of one slot, m/s.
    pure function closure_saturated_conductivity(slot) result(value)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64) :: value

        value = saturated_conductivity(slot)
    end function closure_saturated_conductivity

    !> Release the closure, so the next solver starts from nothing.
    subroutine clear_closure()
        implicit none

        if (num_closure_materials > 0) then
            call hydraulic_model%destroy()
            call thermal_model%destroy()
        end if
        material_ids(:) = 0
        material_porosity(:) = 0.0d0
        saturated_conductivity(:) = 1.0d0
        num_closure_materials = 0
        slot_filled(:) = .false.
        slot_next = 1
        ! The memo keys on a material slot, and the next run renumbers those.
        call forget_closure_tangents()
    end subroutine clear_closure

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
    subroutine publish_closure_state(slot, temperature, pressure, state, ok)
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
        if (slot < 1 .or. slot > num_closure_materials) return

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
        call state%continuation_lambda%set(1.0d0)
        call hydraulic_model%update_water_phases(material_ids(slot), state)
        ok = .true.
    end subroutine publish_closure_state

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
    subroutine eval_closure_response(slot, temperature, pressure, response, ok)
        implicit none
        !> One-based index into the frozen closure's material table.
        integer(int32), intent(in) :: slot
        !> Temperature, C.
        real(real64), intent(in) :: temperature
        !> Total-water-equivalent pore pressure \(p_w\), Pa.
        real(real64), intent(in) :: pressure
        !> The response; written only when `ok`.
        type(type_closure_response), intent(inout) :: response
        !> True when every model supplied a usable value.
        logical, intent(inout) :: ok

        type(type_closure_response) :: local
        integer(int32) :: entry
        logical :: hit

        ok = .false.
        response_requests = response_requests + 1_int64

        entry = 0
        hit = .false.
        if (allocated(response_cache_filled)) then
            call find_memo_entry(response_cache_material, response_cache_temperature, response_cache_pressure, &
                                 response_cache_filled, slot, temperature, pressure, entry, hit)
        end if
        if (hit) then
            response = response_cache_value(entry)
            ok = .true.
            return
        end if

        call compute_closure_response(slot, temperature, pressure, local, ok)
        if (.not. ok) return
        response = local

        ! Stored only on success, so a refused key is refused again every time.
        if (entry > 0) then
            response_cache_value(entry) = local
            response_cache_material(entry) = slot
            response_cache_temperature(entry) = temperature
            response_cache_pressure(entry) = pressure
            response_cache_filled(entry) = .true.
        end if
    end subroutine eval_closure_response

    !> The response of `eval_closure_response`, computed afresh and remembered nowhere.
    !>
    !> For a state no kernel asks for twice -- an accepted state that has a memo
    !> of its own -- so that it evicts nothing from the response memo, and for a
    !> memo miss. With `tangents` present, also the tangents of
    !> `eval_closure_tangents` at the same published state. Fails as
    !> `eval_closure_response` does, and when a tangent is not finite.
    subroutine compute_closure_response(slot, temperature, pressure, response, ok, tangents)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        !> The response; written only when `ok`.
        type(type_closure_response), intent(inout) :: response
        logical, intent(inout) :: ok
        !> The tangents; written only when `ok`.
        type(type_closure_tangents), intent(inout), optional :: tangents

        type(type_closure_response) :: local
        type(type_closure_tangents) :: slopes
        real(real64) :: suction, value, specific_heat
        real(real64) :: dwater(2), denthalpy(2), dlog_conductivity(2), dthermal(2)
        integer(int32) :: which

        ok = .false.
        response_evaluations = response_evaluations + 1_int64
        call ensure_state(slot, temperature, pressure, which)
        if (which == 0) return

        if (present(tangents)) then
            call relative_conductivity(slot, which, local%relative_conductivity, ok, dlog_conductivity)
            if (.not. ok) return
            call thermal_conductivity_value(slot, which, local%thermal_conductivity, ok, dthermal)
            if (.not. ok) return
            call hydraulic_model%compute_total_storage(material_ids(slot), slot_state(which), local%stored_water, dwater)
            call thermal_model%calc_enthalpy_density(material_ids(slot), slot_state(which), local%enthalpy, denthalpy)
        else
            call relative_conductivity(slot, which, local%relative_conductivity, ok)
            if (.not. ok) return
            call thermal_conductivity_value(slot, which, local%thermal_conductivity, ok)
            if (.not. ok) return
            call hydraulic_model%compute_total_storage(material_ids(slot), slot_state(which), local%stored_water)
            call thermal_model%calc_enthalpy_density(material_ids(slot), slot_state(which), local%enthalpy)
        end if
        call slot_state(which)%effective_suction%get(suction)
        local%liquid_pressure = -suction
        call slot_state(which)%dSeff_dP%get(suction)
        local%liquid_dpressure = -suction
        call slot_state(which)%dSeff_dT%get(suction)
        local%liquid_dtemperature = -suction
        call slot_state(which)%water_content%get(local%liquid_content)
        call slot_state(which)%ice_content%get(local%ice_content)
        call slot_state(which)%vapor_content%get(local%vapor_content)
        call thermal_model%calc_density_water(slot_state(which), local%water_density)
        call thermal_model%calc_specific_heat_water(slot_state(which), specific_heat)
        local%liquid_heat_capacity = local%water_density*specific_heat
        call thermal_model%calc_density_ice(slot_state(which), local%ice_density)

        ! The compressive term as the difference against the phase sum. Reading
        ! it back out of W rather than re-forming S_s/(rho g) max(0,p-p_sat)
        ! keeps the branch flag tied to the storage the residual differences,
        ! including the saturation pressure the phase model decided on.
        local%compressive_storage = local%stored_water - (local%liquid_content + local%vapor_content)
        if (local%water_density > 0.0d0) then
            local%compressive_storage = local%compressive_storage &
                                        - (local%ice_density/local%water_density)*local%ice_content
        end if
        local%ice_present = (local%ice_content > 0.0d0)
        local%compressive_active = (local%compressive_storage > 0.0d0)

        value = local%stored_water + local%enthalpy + local%liquid_pressure + local%thermal_conductivity
        if (.not. ieee_is_finite(value)) then
            ok = .false.
            return
        end if

        if (present(tangents)) then
            slopes%dwater_dpressure = dwater(tangent_pressure)
            slopes%dwater_dtemperature = dwater(tangent_temperature)
            slopes%denthalpy_dpressure = denthalpy(tangent_pressure)
            slopes%denthalpy_dtemperature = denthalpy(tangent_temperature)
            slopes%dliquid_dpressure = local%liquid_dpressure
            slopes%dliquid_dtemperature = local%liquid_dtemperature
            slopes%dlog_conductivity_dpressure = dlog_conductivity(tangent_pressure)
            slopes%dlog_conductivity_dtemperature = dlog_conductivity(tangent_temperature)
            slopes%dthermal_dpressure = dthermal(tangent_pressure)
            slopes%dthermal_dtemperature = dthermal(tangent_temperature)
            ! A state the models answer with a value but not with a derivative
            ! is one the Jacobian cannot linearise about.
            value = sum(dwater) + sum(denthalpy) + sum(dlog_conductivity) + sum(dthermal)
            if (.not. ieee_is_finite(value)) then
                ok = .false.
                return
            end if
            tangents = slopes
        end if
        response = local
    end subroutine compute_closure_response

    !> \(W\) and \(U\) of an accepted state, the history a BDF rate differences.
    !>
    !> The same values `eval_closure_response` returns at that state, from a
    !> memo of their own (see `storage_cache_water`). Fails as
    !> `eval_closure_response` does, leaving both results untouched.
    subroutine eval_closure_storage(slot, temperature, pressure, stored_water, enthalpy, ok)
        implicit none
        integer(int32), intent(in) :: slot
        !> Accepted temperature, C.
        real(real64), intent(in) :: temperature
        !> Accepted pore pressure \(p_w\), Pa.
        real(real64), intent(in) :: pressure
        !> \(W\), dimensionless.
        real(real64), intent(inout) :: stored_water
        !> \(U\), J/m3.
        real(real64), intent(inout) :: enthalpy
        logical, intent(inout) :: ok

        type(type_closure_response) :: local
        integer(int32) :: entry
        logical :: hit

        ok = .false.
        response_requests = response_requests + 1_int64

        entry = 0
        hit = .false.
        if (allocated(storage_cache_filled)) then
            call find_memo_entry(storage_cache_material, storage_cache_temperature, storage_cache_pressure, &
                                 storage_cache_filled, slot, temperature, pressure, entry, hit)
        end if
        if (hit) then
            stored_water = storage_cache_water(entry)
            enthalpy = storage_cache_enthalpy(entry)
            ok = .true.
            return
        end if

        call compute_closure_response(slot, temperature, pressure, local, ok)
        if (.not. ok) return
        stored_water = local%stored_water
        enthalpy = local%enthalpy

        if (entry > 0) then
            storage_cache_water(entry) = local%stored_water
            storage_cache_enthalpy(entry) = local%enthalpy
            storage_cache_material(entry) = slot
            storage_cache_temperature(entry) = temperature
            storage_cache_pressure(entry) = pressure
            storage_cache_filled(entry) = .true.
        end if
    end subroutine eval_closure_storage

    !> \(\lambda\) of one material at a published state, W/(m K).
    !>
    !> Read through the governing model's own thermal diffusion tensor, so the
    !> phase mixing rule is the model's. The tensor has to be isotropic: the
    !> coupled conduction term carries a scalar, and reporting one
    !> entry of an anisotropic tensor as that scalar would answer a different
    !> problem.
    subroutine thermal_conductivity_value(slot, which, value, ok, tangent)
        implicit none
        integer(int32), intent(in) :: slot
        !> The cached state slot, as `ensure_state` returned it.
        integer(int32), intent(in) :: which
        real(real64), intent(inout) :: value
        logical, intent(inout) :: ok
        !> \((\partial\lambda/\partial p_w,\partial\lambda/\partial T)\); the tensor is isotropic, so this is
        !> the tangent of its scalar.
        real(real64), intent(inout), optional :: tangent(2)

        real(real64) :: tensor(2, 2), scale

        ok = .false.
        tensor(:, :) = 0.0d0
        call thermal_model%compute_diffusion_term(material_ids(slot), slot_state(which), tensor, tangent)

        scale = max(abs(tensor(1, 1)), abs(tensor(2, 2)))
        if (.not. ieee_is_finite(scale)) return
        if (scale <= 0.0d0) return
        if (abs(tensor(1, 1) - tensor(2, 2)) > ISOTROPY_TOLERANCE*scale) return
        if (max(abs(tensor(1, 2)), abs(tensor(2, 1))) > ISOTROPY_TOLERANCE*scale) return

        value = tensor(1, 1)
        ok = .true.
    end subroutine thermal_conductivity_value

    !> \(\partial p_l/\partial p_w\) and \(\partial p_l/\partial T\) alone.
    !>
    !> The phase split publishes these analytically, so a Darcy driver needs
    !> no differenced tangent: the values are those `eval_closure_tangents`
    !> returns in `dliquid_dpressure` and `dliquid_dtemperature`, at the cost of
    !> one cached response instead of five.
    subroutine eval_liquid_slopes(slot, temperature, pressure, dliquid_dpressure, dliquid_dtemperature, ok)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        real(real64), intent(inout) :: dliquid_dpressure, dliquid_dtemperature
        logical, intent(inout) :: ok

        type(type_closure_response) :: response

        call eval_closure_response(slot, temperature, pressure, response, ok)
        if (.not. ok) return
        dliquid_dpressure = response%liquid_dpressure
        dliquid_dtemperature = response%liquid_dtemperature
    end subroutine eval_liquid_slopes

    !> The tangents of that response in both unknowns.
    !>
    !> \[
    !>   \Bigl(\frac{\partial F}{\partial p_w},\frac{\partial F}{\partial T}\Bigr),\qquad
    !>   F\in\{W,\,U,\,p_l,\,\ln k_r,\,\lambda\},
    !> \]
    !>
    !> the derivatives of the configured models at the published state, through
    !> the analytic tangents of the phase split; see "Tangents" above. Each is
    !> the derivative on the branch the state is on. Exact up to rounding for
    !> the models whose own derivatives are exact. \(O(1)\), one phase split on a
    !> memo miss and none on a hit.
    !>
    !> Fails with `ok = .false.` and untouched results whenever
    !> `eval_closure_response` does at the point, and when a tangent is not
    !> finite.
    subroutine eval_closure_tangents(slot, temperature, pressure, tangents, ok)
        implicit none
        !> One-based index into the frozen closure's material table.
        integer(int32), intent(in) :: slot
        !> Temperature, C.
        real(real64), intent(in) :: temperature
        !> Total-water-equivalent pore pressure \(p_w\), Pa.
        real(real64), intent(in) :: pressure
        !> The tangents; written only when `ok`.
        type(type_closure_tangents), intent(inout) :: tangents
        !> True when every model supplied a usable value.
        logical, intent(inout) :: ok

        type(type_closure_tangents) :: local
        type(type_closure_response) :: base
        integer(int32) :: entry, response_entry
        logical :: hit

        ok = .false.
        tangent_requests = tangent_requests + 1_int64

        ! A hit returns the value this same key already produced, so the caller
        ! cannot tell a hit from a miss by the number it receives.
        entry = 0
        hit = .false.
        if (allocated(tangent_cache_filled)) then
            call find_memo_entry(tangent_cache_material, tangent_cache_temperature, tangent_cache_pressure, &
                                 tangent_cache_filled, slot, temperature, pressure, entry, hit)
        end if
        if (hit) then
            tangents = tangent_cache_value(entry)
            ok = .true.
            return
        end if

        tangent_evaluations = tangent_evaluations + 1_int64
        response_requests = response_requests + 1_int64
        call compute_closure_response(slot, temperature, pressure, base, ok, local)
        if (.not. ok) return
        tangents = local

        ! Stored only now, so that a key the closure refused is never answered
        ! from the memo and keeps reporting its fault. The response came with
        ! the tangents, so it is remembered for the kernels that ask next.
        if (entry > 0) then
            tangent_cache_value(entry) = local
            tangent_cache_material(entry) = slot
            tangent_cache_temperature(entry) = temperature
            tangent_cache_pressure(entry) = pressure
            tangent_cache_filled(entry) = .true.
        end if
        response_entry = 0
        if (allocated(response_cache_filled)) then
            call find_memo_entry(response_cache_material, response_cache_temperature, response_cache_pressure, &
                                 response_cache_filled, slot, temperature, pressure, response_entry, hit)
        end if
        if (response_entry > 0) then
            response_cache_value(response_entry) = base
            response_cache_material(response_entry) = slot
            response_cache_temperature(response_entry) = temperature
            response_cache_pressure(response_entry) = pressure
            response_cache_filled(response_entry) = .true.
        end if
    end subroutine eval_closure_tangents

    !> \(k_r=K_{fl h}(T,p_w)/K_{\mathrm{sat}}\) of one cell.
    subroutine relative_conductivity(slot, which, relative, ok, log_tangent)
        implicit none
        integer(int32), intent(in) :: slot
        !> The cached state slot, as `ensure_state` returned it.
        integer(int32), intent(in) :: which
        real(real64), intent(inout) :: relative
        logical, intent(inout) :: ok
        !> \((\partial\ln k_r/\partial p_w,\partial\ln k_r/\partial T)\); \(K_{\mathrm{sat}}\) is constant.
        real(real64), intent(inout), optional :: log_tangent(2)

        real(real64) :: value

        call conductivity_value(slot, which, value, ok, log_tangent)
        if (.not. ok) return

        value = value/saturated_conductivity(slot)
        ok = ieee_is_finite(value) .and. value > 0.0d0
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
    subroutine conductivity_value(slot, which, value, ok, log_tangent)
        implicit none
        integer(int32), intent(in) :: slot
        !> The cached state slot, as `ensure_state` returned it.
        integer(int32), intent(in) :: which
        real(real64), intent(inout) :: value
        logical, intent(inout) :: ok
        !> \((\partial\ln K/\partial p_w,\partial\ln K/\partial T)\), from the HCF model's own factors.
        real(real64), intent(inout), optional :: log_tangent(2)

        real(real64) :: tensor(2, 2), density

        ok = .false.
        tensor(:, :) = 0.0d0
        call hydraulic_model%compute_diffusion_term(material_ids(slot), slot_state(which), tensor, log_tangent)
        call thermal_model%calc_density_water(slot_state(which), density)
        ! compute_diffusion_term forms K/(rho_w g) with the same two quantities,
        ! so multiplying them back recovers the model's own conductivity without
        ! restating the model.
        value = tensor(1, 1)*density*gravity_acceleration
        ok = ieee_is_finite(value) .and. value > 0.0d0
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
        if (.not. (temperature_K > tiny(1.0d0))) return
        rh_exponent = (pressure * rh_molar_mass_water)/(rh_reference_water_density * rh_gas_constant * temperature_K)
        if (.not. ieee_is_finite(rh_exponent)) return
        state_is_admissible = abs(rh_exponent) <= 700.0d0
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
        if (slot < 1 .or. slot > num_closure_materials) return

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
        call publish_closure_state(slot, temperature, pressure, slot_state(which), ok)
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
    subroutine configure_closure_cache(enabled, entries)
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
        if (allocated(storage_cache_water)) then
            deallocate (storage_cache_water, storage_cache_enthalpy, storage_cache_temperature, &
                        storage_cache_pressure, storage_cache_material, storage_cache_filled)
        end if
        tangent_cache_mask = 0

        tangent_cache_enabled = enabled
        if (.not. enabled) return

        wanted = entries
        if (wanted <= 0) wanted = TANGENT_CACHE_DEFAULT_ENTRIES

        size_now = MEMO_WAYS
        do while (size_now < wanted .and. size_now < huge(size_now)/2)
            size_now = size_now*2
        end do

        allocate (tangent_cache_value(size_now))
        allocate (tangent_cache_temperature(size_now), source=0.0d0)
        allocate (tangent_cache_pressure(size_now), source=0.0d0)
        allocate (tangent_cache_material(size_now), source=0_int32)
        allocate (tangent_cache_filled(size_now), source=.false.)

        allocate (response_cache_value(size_now))
        allocate (response_cache_temperature(size_now), source=0.0d0)
        allocate (response_cache_pressure(size_now), source=0.0d0)
        allocate (response_cache_material(size_now), source=0_int32)
        allocate (response_cache_filled(size_now), source=.false.)

        allocate (storage_cache_water(size_now), source=0.0d0)
        allocate (storage_cache_enthalpy(size_now), source=0.0d0)
        allocate (storage_cache_temperature(size_now), source=0.0d0)
        allocate (storage_cache_pressure(size_now), source=0.0d0)
        allocate (storage_cache_material(size_now), source=0_int32)
        allocate (storage_cache_filled(size_now), source=.false.)

        tangent_cache_mask = size_now - 1
    end subroutine configure_closure_cache

    !> How much constitutive work was asked for and how much was done.
    !>
    !> Counted since the closure was installed, on the calling rank only.
    !> `requests` counts calls to `eval_closure_tangents` and `evaluations`
    !> counts those that actually computed a tangent, so their ratio is the
    !> repetition the memo removes. \(O(1)\). Cannot fail.
    subroutine closure_work_counters(requests, evaluations, asked, responses, publications)
        implicit none
        !> Calls to `eval_closure_tangents`.
        integer(int64), intent(inout) :: requests
        !> Those of them that computed a tangent rather than recalling one.
        integer(int64), intent(inout) :: evaluations
        !> Calls to `eval_closure_response`.
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
    end subroutine closure_work_counters

    !> Drop every memoised value, keeping the configured size.
    subroutine forget_closure_tangents()
        implicit none

        if (allocated(tangent_cache_filled)) tangent_cache_filled(:) = .false.
        if (allocated(response_cache_filled)) response_cache_filled(:) = .false.
        if (allocated(storage_cache_filled)) storage_cache_filled(:) = .false.
    end subroutine forget_closure_tangents

    !> The entry of a memo table that holds a key, or the one a new key is to take.
    !>
    !> The key's set is the `MEMO_WAYS` entries its mixed index falls in. A key
    !> found in the set is a hit. Otherwise the entry is the first empty way,
    !> or, in a full set, the way `memo_victim` names. Zero, and no hit, when
    !> the memo is off. \(O(1)\). Cannot fail.
    subroutine find_memo_entry(materials, temperatures, pressures, filled, slot, temperature, pressure, entry, hit)
        implicit none
        !> Material slot of every entry.
        integer(int32), intent(in) :: materials(:)
        !> Temperature key of every entry, C.
        real(real64), intent(in) :: temperatures(:)
        !> Pressure key of every entry, Pa.
        real(real64), intent(in) :: pressures(:)
        !> Whether each entry holds a value.
        logical, intent(in) :: filled(:)
        !> One-based index into the frozen closure's material table.
        integer(int32), intent(in) :: slot
        !> Temperature, C.
        real(real64), intent(in) :: temperature
        !> Total-water-equivalent pore pressure \(p_w\), Pa.
        real(real64), intent(in) :: pressure
        !> The entry holding the key, or the one to store it in; zero when the memo is off.
        integer(int32), intent(inout) :: entry
        !> Whether `entry` already holds the key.
        logical, intent(inout) :: hit

        integer(int32) :: first, way, candidate

        entry = 0
        hit = .false.
        if (.not. tangent_cache_enabled) return

        first = iand(mixed_key_entry(slot, temperature, pressure) - 1, not(MEMO_WAYS - 1)) + 1
        do way = 0, MEMO_WAYS - 1
            candidate = first + way
            if (.not. filled(candidate)) then
                if (entry == 0) entry = candidate
                cycle
            end if
            if (materials(candidate) == slot .and. temperatures(candidate) == temperature .and. &
                pressures(candidate) == pressure) then
                entry = candidate
                hit = .true.
                return
            end if
        end do
        if (entry > 0) return
        memo_victim = iand(memo_victim + 1, MEMO_WAYS - 1)
        entry = first + memo_victim
    end subroutine find_memo_entry

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

end module app_ftcms_closure
