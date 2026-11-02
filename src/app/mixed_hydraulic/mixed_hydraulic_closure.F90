#include <petsc/finclude/petscsys.h>

!> The cell-local constitutive closure: stored water and relative conductivity.
!>
!> Both are functions of the cell's pressure and its material, and both are
!> defined here once so that the residual, the Jacobian and the balance audit
!> cannot drift apart. A cell reaches its own law by the slot its auxiliary
!> material field carries; nothing interpolates a material and nothing averages
!> two of them.
!>
!> The conservation residual differences a stored quantity,
!>
!> \[
!>   \int_K\frac{W(p^{n+1})-W(p^{n})}{\Delta t},
!> \]
!>
!> and its Jacobian carries \(\partial W/\partial p\) at the current iterate.
!> The two have to come from the same function or the matrix has a diagonal no
!> residual term corresponds to: that is exactly what
!> \(C(p^{n+1})(p^{n+1}-p^{n})/\Delta t\) does, and it is what this module
!> exists to prevent. The residual kernel, the Jacobian kernel and the balance
!> audit all call `eval_storage`; nothing reassembles a storage term of its own.
!>
!> The retention law is the repository's own `holder_wrfs` object, configured
!> from a `type_config_wrf`, so the curve solved here is the curve the
!> application solves. The saturated branch is the compressive storage of
!> `compute_compressive_storage_hydraulic`:
!>
!> \[
!>   W(p)=\theta_{\mathrm{WRF}}\!\left(\frac{p}{\rho_{\mathrm{std}}g}\right)
!>        +\frac{S_s}{\rho_{\mathrm{std}}g}\max(0,p-p_{\mathrm{sat}}),
!> \]
!>
!> with \(h=p/(\rho_{\mathrm{std}}g)\) in m, as `calc_head_hcf` also forms it.
!> Nothing is floored and no curve is replaced: where the retention tangent
!> vanishes the compressive term is the whole capacity, which is the state a
!> saturated cell is actually in.
!>
!> The relative conductivity is the repository's own hydraulic conductivity
!> function, normalised by its own saturated value,
!>
!> \[
!>   k_r(p)=\frac{K_{fl h}(p)}{K_{fl h}(0)},
!> \]
!>
!> so \(k_r(p)=1\) exactly wherever the model is saturated, and the saturated
!> conductivity of the cell is the material's tensor rather than the scalar
!> `k_sat` of the curve. The model exposes no derivative, so the tangent this
!> module reports is a central difference of \(\ln K_{fl h}\) in the head. That
!> is a differentiated evaluation of the law itself, not a second statement of
!> it: re-deriving \(k_r'\) here would duplicate a constitutive model and would
!> be wrong for the incomplete-beta branch of the van Genuchten curve.
!>
!> The table is module state because a PetscDS pointwise function is reached
!> through a C function pointer and carries no typed context. One solver per
!> rank owns it, and PETSc assembles serially within a rank, which is the same
!> argument `app_plexfem_physics` makes for the same reason. A second live
!> solver would share it, so installing one while another holds the table is
!> refused rather than silently answered with the wrong material.
module app_mixed_hydraulic_closure
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use :: petscsys
    use :: iapws, only:type_iapws97, type_iapws06
    use :: core_physical_constants, only:gravity_acceleration, reference_water_density
    use :: module_core, only:type_state, HCF_MODES
    use :: models_wrf, only:holder_wrfs
    use :: models_hcf, only:holder_hcfs
    use :: app_mixed_hydraulic_physics, only:type_mixed_material, &
                                             MIXED_STORAGE_LINEAR, MIXED_STORAGE_RETENTION, &
                                             MIXED_CONDUCTIVITY_FIXED, MIXED_CONDUCTIVITY_HCF
    implicit none
    private

    public :: install_storage_laws
    public :: clear_storage_laws
    public :: storage_laws_installed
    public :: eval_storage
    public :: storage_magnitude_of
    public :: MAX_STORAGE_LAWS
    public :: install_conductivity_laws
    public :: clear_conductivity_laws
    public :: conductivity_laws_installed
    public :: eval_relative_conductivity
    public :: eval_relative_conductivity_tangent
    public :: clear_mobility_fault
    public :: note_mobility_fault
    public :: mobility_fault_recorded
    public :: mobility_fault_report

    !> How many materials one run may carry storage laws for.
    integer(int32), parameter :: MAX_STORAGE_LAWS = 16

    !> \(\rho_{\mathrm{std}}g\), the factor between a pore pressure in Pa and a
    !> suction head in m. The same constants the production closure uses, not
    !> the run's own gravity: the retention parameters were measured with these.
    real(real64), parameter :: SPECIFIC_WEIGHT = reference_water_density*gravity_acceleration

    !> One material's storage law.
    type :: type_mixed_storage_law
        !> MIXED_STORAGE_LINEAR or MIXED_STORAGE_RETENTION.
        integer(int32) :: kind = MIXED_STORAGE_LINEAR
        !> Constant capacity of the linear law, 1/Pa.
        real(real64) :: capacity = 0.0_real64
        !> The retention curve of the nonlinear law.
        type(holder_wrfs) :: retention
        !> \(S_s/(\rho_{\mathrm{std}}g)\) of the compressive branch, 1/Pa.
        real(real64) :: compressive_capacity = 0.0_real64
        !> Pressure the compressive branch starts at, Pa.
        real(real64) :: saturation_pressure = 0.0_real64
        !> A positive number whenever the law can store water, 1/Pa.
        real(real64) :: magnitude = 0.0_real64
    end type type_mixed_storage_law

    !> One material's relative conductivity law.
    type :: type_mixed_conductivity_law
        !> MIXED_CONDUCTIVITY_FIXED or MIXED_CONDUCTIVITY_HCF.
        integer(int32) :: kind = MIXED_CONDUCTIVITY_FIXED
        !> The conductivity function of the nonlinear law.
        type(holder_hcfs) :: curve
        !> \(K_{fl h}(0)\), the saturated value the curve is normalised by, m/s.
        real(real64) :: saturated_value = 1.0_real64
    end type type_mixed_conductivity_law

    type(type_mixed_storage_law), save :: laws(MAX_STORAGE_LAWS)
    integer(int32), save :: num_laws = 0

    type(type_mixed_conductivity_law), save :: conductivity_laws(MAX_STORAGE_LAWS)
    integer(int32), save :: num_conductivity_laws = 0

    !> The water and ice property objects the conductivity models point at.
    !>
    !> They outlive the laws because `abst_hcf` keeps pointers to them, and a
    !> pointer into a temporary would dangle the moment installation returned.
    type(type_iapws97), save, target :: water_properties
    type(type_iapws06), save, target :: ice_properties

    !> Relative step of the head the conductivity tangent is differenced over.
    !>
    !> The cube root of the machine epsilon balances the truncation error of a
    !> central difference against the cancellation of the subtraction.
    real(real64), parameter :: HEAD_STEP_FRACTION = epsilon(1.0_real64)**(1.0_real64/3.0_real64)
    !> The head scale the step falls back to near \(h=0\), m.
    !>
    !> A retention or conductivity curve varies over a head of order one metre,
    !> so a step proportional to \(|h|\) alone would collapse to roundoff there.
    real(real64), parameter :: HEAD_STEP_FLOOR = 1.0_real64

    ! What the last unusable mobility evaluation was, for the diagnostic a
    ! refused step prints. A pointwise kernel has nowhere to return an error to,
    ! so it records the fault and lets the nonfinite residual refuse the step.
    logical, save :: mobility_fault = .false.
    integer(int32), save :: mobility_fault_slot = 0
    real(real64), save :: mobility_fault_pressure = 0.0_real64
    real(real64), save :: mobility_fault_value = 0.0_real64

contains

    !> Build the storage law of every material, in the order given.
    !>
    !> A material's slot is its index in `materials`, which is the number the
    !> auxiliary material field carries. Costs one retention-model construction
    !> per material and \(O(1)\) memory each.
    !>
    !> Fails by `error stop` when the table is already held by a live solver,
    !> when there are more materials than slots, or when a retention law is
    !> asked for without a usable curve.
    subroutine install_storage_laws(materials)
        implicit none
        !> One entry per material of the mesh.
        type(type_mixed_material), intent(in) :: materials(:)

        integer(int32) :: m

        if (num_laws /= 0) then
            error stop "A second mixed hydraulic solver cannot share the storage-law table of a live one"
        end if
        if (size(materials) > MAX_STORAGE_LAWS) then
            error stop "More materials than the mixed hydraulic storage-law table can hold"
        end if

        do m = 1, size(materials)
            laws(m) = type_mixed_storage_law()
            laws(m)%kind = materials(m)%storage_law
            laws(m)%saturation_pressure = materials(m)%saturation_pressure
            laws(m)%compressive_capacity = materials(m)%specific_storage/SPECIFIC_WEIGHT

            select case (materials(m)%storage_law)
            case (MIXED_STORAGE_LINEAR)
                laws(m)%capacity = materials(m)%storage
                ! The linear law has no saturated branch of its own: its
                ! capacity is the whole of it.
                laws(m)%compressive_capacity = 0.0_real64
                laws(m)%magnitude = abs(materials(m)%storage)

            case (MIXED_STORAGE_RETENTION)
                if (materials(m)%specific_storage < 0.0_real64) then
                    error stop "A mixed hydraulic specific storage cannot be negative"
                end if
                call laws(m)%retention%initialize(materials(m)%retention)
                if (.not. allocated(laws(m)%retention%p)) then
                    error stop "A mixed hydraulic retention law was asked for without a curve"
                end if
                call set_retention_magnitude(laws(m))

            case default
                error stop "A storage law of this kind is not one the mixed hydraulic solver knows"
            end select
        end do

        num_laws = int(size(materials), int32)
    end subroutine install_storage_laws

    !> Release the table, so the next solver starts from nothing.
    subroutine clear_storage_laws()
        implicit none

        integer(int32) :: m

        do m = 1, MAX_STORAGE_LAWS
            if (allocated(laws(m)%retention%p)) deallocate (laws(m)%retention%p)
            laws(m) = type_mixed_storage_law()
        end do
        num_laws = 0
    end subroutine clear_storage_laws

    !> Whether a run currently holds the table.
    logical function storage_laws_installed()
        implicit none

        storage_laws_installed = (num_laws > 0)
    end function storage_laws_installed

    !> The stored water of one cell and its pressure tangent.
    !>
    !> \[
    !>   W(p),\qquad \frac{\partial W}{\partial p}(p),
    !> \]
    !>
    !> both per unit bulk volume; \(W\) is dimensionless and the tangent is in
    !> 1/Pa. Assumes `slot` is a material index `install_storage_laws` was given.
    !> Exact for the configured law up to its own evaluation; the tangent is the
    !> analytic derivative of the value returned beside it, except at
    !> \(p=p_{\mathrm{sat}}\), where the compressive branch starts and the
    !> closure's capacity is one-sided. \(O(1)\).
    !>
    !> Fails by returning zeros for a slot outside the table, which is what an
    !> unset material field would ask for.
    subroutine eval_storage(slot, pressure, storage, capacity)
        implicit none
        !> One-based index into the storage-law table.
        integer(int32), intent(in) :: slot
        !> Pore pressure, Pa; the whole real line is valid.
        real(real64), intent(in) :: pressure
        !> \(W\), stored water per unit bulk volume.
        real(real64), intent(inout) :: storage
        !> \(\partial W/\partial p\), 1/Pa.
        real(real64), intent(inout) :: capacity

        real(real64) :: head, water_content, dwater_dhead

        storage = 0.0_real64
        capacity = 0.0_real64
        if (slot < 1 .or. slot > num_laws) return

        select case (laws(slot)%kind)
        case (MIXED_STORAGE_LINEAR)
            storage = laws(slot)%capacity*pressure
            capacity = laws(slot)%capacity

        case (MIXED_STORAGE_RETENTION)
            head = pressure/SPECIFIC_WEIGHT
            water_content = 0.0_real64
            dwater_dhead = 0.0_real64
            call laws(slot)%retention%p%calc(head, water_content)
            call laws(slot)%retention%p%deriv(head, dwater_dhead)
            storage = water_content
            capacity = dwater_dhead/SPECIFIC_WEIGHT

            ! The saturated compressive branch, as the production closure forms
            ! it: zero below the saturation pressure, linear above it. Its
            ! capacity is the whole of the tangent where the retention curve has
            ! flattened, which is what keeps the pressure row of a saturated
            ! cell non-degenerate.
            if (pressure >= laws(slot)%saturation_pressure) then
                storage = storage + laws(slot)%compressive_capacity &
                          *(pressure - laws(slot)%saturation_pressure)
                capacity = capacity + laws(slot)%compressive_capacity
            end if
        end select
    end subroutine eval_storage

    !> How much this material's law can store at all, 1/Pa.
    !>
    !> Zero exactly when the cell is steady. A property of the law rather than
    !> of one state, so a curve whose tangent vanishes at the accepted pressure
    !> is still reported as storing water. \(O(1)\); zero for an unknown slot.
    function storage_magnitude_of(slot) result(magnitude)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64) :: magnitude

        magnitude = 0.0_real64
        if (slot < 1 .or. slot > num_laws) return
        magnitude = laws(slot)%magnitude
    end function storage_magnitude_of

    !> The storage magnitude of a retention law, from the curve's own bound.
    !>
    !> `calc_lscheme_capacity` is the retention model's precomputed
    !> \(\max_h(\mathrm d\theta/\mathrm dh)/(\rho_{\mathrm{std}}g)\) over the
    !> material's characteristic heads, already in 1/Pa, so it is zero only for
    !> a curve that stores nothing anywhere.
    subroutine set_retention_magnitude(law)
        implicit none
        type(type_mixed_storage_law), intent(inout) :: law

        real(real64) :: pressure_capacity

        pressure_capacity = 0.0_real64
        call law%retention%p%calc_lscheme_capacity(pressure_capacity)
        law%magnitude = abs(pressure_capacity) + abs(law%compressive_capacity)
    end subroutine set_retention_magnitude

    !> Build the relative conductivity law of every material, in the order given.
    !>
    !> A material's slot is its index in `materials`, the same number the storage
    !> table is keyed by. Costs one conductivity-model construction and one
    !> saturated evaluation per material, and \(O(1)\) memory each.
    !>
    !> Fails by `error stop` when the table is already held by a live solver,
    !> when there are more materials than slots, when a conductivity function is
    !> asked for in a mode this isothermal path cannot evaluate, or when its
    !> saturated value is not a positive finite number. No value stands in for a
    !> conductivity the model does not supply.
    subroutine install_conductivity_laws(materials)
        implicit none
        !> One entry per material of the mesh.
        type(type_mixed_material), intent(in) :: materials(:)

        integer(int32) :: m
        real(real64) :: saturated

        if (num_conductivity_laws /= 0) then
            error stop "A second mixed hydraulic solver cannot share the conductivity table of a live one"
        end if
        if (size(materials) > MAX_STORAGE_LAWS) then
            error stop "More materials than the mixed hydraulic conductivity table can hold"
        end if

        call water_properties%initialize()
        call ice_properties%initialize()

        do m = 1, size(materials)
            conductivity_laws(m)%kind = materials(m)%conductivity_law
            conductivity_laws(m)%saturated_value = 1.0_real64

            select case (materials(m)%conductivity_law)
            case (MIXED_CONDUCTIVITY_FIXED)
                ! The tensor is the whole law; nothing to build.

            case (MIXED_CONDUCTIVITY_HCF)
                ! Only the retention-driven mode is evaluable here. The
                ! impedance modes read an ice content and the viscosity modes a
                ! temperature, and this path carries neither: refusing says so
                ! instead of handing the model an unset state.
                if (.not. (materials(m)%conductivity_function%model == HCF_MODES%BASE)) then
                    error stop "The isothermal mixed hydraulic path evaluates the base conductivity mode only: "// &
                        "the impedance and viscosity modes need an ice content and a temperature"
                end if
                call conductivity_laws(m)%curve%initialize(materials(m)%conductivity_function, &
                                                           water_properties, ice_properties)
                if (.not. allocated(conductivity_laws(m)%curve%p)) then
                    error stop "A mixed hydraulic conductivity law was asked for without a curve"
                end if

                ! The normaliser, so the saturated conductivity of the cell is
                ! its own tensor whatever scale the curve was configured with.
                saturated = conductivity_value(m, 0.0_real64)
                if (.not. ieee_is_finite(saturated)) then
                    error stop "A mixed hydraulic conductivity function returned a nonfinite saturated value"
                end if
                if (saturated <= 0.0_real64) then
                    error stop "A mixed hydraulic conductivity function has no positive saturated value: "// &
                        "an impermeable material is a different model, not a small number"
                end if
                conductivity_laws(m)%saturated_value = saturated

            case default
                error stop "A conductivity law of this kind is not one the mixed hydraulic solver knows"
            end select
        end do

        num_conductivity_laws = int(size(materials), int32)
    end subroutine install_conductivity_laws

    !> Release the conductivity table, so the next solver starts from nothing.
    subroutine clear_conductivity_laws()
        implicit none

        integer(int32) :: m

        do m = 1, MAX_STORAGE_LAWS
            if (allocated(conductivity_laws(m)%curve%p)) deallocate (conductivity_laws(m)%curve%p)
            conductivity_laws(m)%kind = MIXED_CONDUCTIVITY_FIXED
            conductivity_laws(m)%saturated_value = 1.0_real64
        end do
        num_conductivity_laws = 0
        call clear_mobility_fault()
    end subroutine clear_conductivity_laws

    !> Whether a run currently holds the conductivity table.
    logical function conductivity_laws_installed()
        implicit none

        conductivity_laws_installed = (num_conductivity_laws > 0)
    end function conductivity_laws_installed

    !> The relative conductivity of one cell.
    !>
    !> \[
    !>   k_r(p)=\frac{K_{fl h}(p)}{K_{fl h}(0)},
    !> \]
    !>
    !> dimensionless. Assumes `slot` is a material index `install_conductivity_laws`
    !> was given. Exact for the configured model up to its own evaluation, and
    !> exactly one wherever that model is saturated. \(O(1)\), one model
    !> evaluation.
    !>
    !> Fails with `ok = .false.` and an untouched result when the model returns a
    !> value that is not positive and finite, and for a slot outside the table.
    !> Nothing is floored: a vanishing conductivity has no inverse mobility, and
    !> substituting a small one would answer a different problem silently.
    subroutine eval_relative_conductivity(slot, pressure, relative, ok)
        implicit none
        !> One-based index into the conductivity table.
        integer(int32), intent(in) :: slot
        !> Pore pressure, Pa; the whole real line is valid.
        real(real64), intent(in) :: pressure
        !> \(k_r\), dimensionless and positive; written only when `ok`.
        real(real64), intent(inout) :: relative
        !> True when the model supplied a usable conductivity.
        logical, intent(inout) :: ok

        real(real64) :: value

        ok = .false.
        if (slot < 1 .or. slot > num_conductivity_laws) return

        if (conductivity_laws(slot)%kind == MIXED_CONDUCTIVITY_FIXED) then
            relative = 1.0_real64
            ok = .true.
            return
        end if

        value = conductivity_value(slot, pressure)/conductivity_laws(slot)%saturated_value
        if (.not. ieee_is_finite(value)) return
        if (value <= 0.0_real64) return

        relative = value
        ok = .true.
    end subroutine eval_relative_conductivity

    !> The relative conductivity of one cell and its logarithmic pressure tangent.
    !>
    !> \[
    !>   k_r(p),\qquad \frac{\mathrm d\ln k_r}{\mathrm dp}
    !>   =\frac{\mathrm d\ln K_{fl h}}{\mathrm dp},
    !> \]
    !>
    !> the second in 1/Pa. The logarithm is what the inverse mobility needs,
    !> \(\partial_p\mathbf\Lambda^{-1}=-(\mathrm d\ln k_r/\mathrm dp)\mathbf\Lambda^{-1}\),
    !> and it is also the better conditioned quantity for a curve that spans
    !> decades. Assumes a smooth branch of the model: the tangent is a central
    !> difference of \(\ln K_{fl h}\) over a head step of
    !> \(\epsilon^{1/3}\max(|h|,1\,\mathrm m)\), so it is second-order accurate
    !> away from the saturation head and a smeared one-sided value across it.
    !> \(O(1)\), three model evaluations.
    !>
    !> Fails with `ok = .false.` and untouched results under the same conditions
    !> as `eval_relative_conductivity`, and additionally when either offset
    !> evaluation is not positive and finite.
    subroutine eval_relative_conductivity_tangent(slot, pressure, relative, dlog_dp, ok)
        implicit none
        !> One-based index into the conductivity table.
        integer(int32), intent(in) :: slot
        !> Pore pressure, Pa; the whole real line is valid.
        real(real64), intent(in) :: pressure
        !> \(k_r\), dimensionless and positive; written only when `ok`.
        real(real64), intent(inout) :: relative
        !> \(\mathrm d\ln k_r/\mathrm dp\), 1/Pa; written only when `ok`.
        real(real64), intent(inout) :: dlog_dp
        !> True when the model supplied a usable conductivity and tangent.
        logical, intent(inout) :: ok

        real(real64) :: head, step, ahead, behind, slope

        ok = .false.
        call eval_relative_conductivity(slot, pressure, relative, ok)
        if (.not. ok) return

        if (conductivity_laws(slot)%kind == MIXED_CONDUCTIVITY_FIXED) then
            dlog_dp = 0.0_real64
            return
        end if

        head = pressure/SPECIFIC_WEIGHT
        step = HEAD_STEP_FRACTION*max(abs(head), HEAD_STEP_FLOOR)

        ahead = conductivity_value(slot, (head + step)*SPECIFIC_WEIGHT)
        behind = conductivity_value(slot, (head - step)*SPECIFIC_WEIGHT)
        if (.not. (ieee_is_finite(ahead) .and. ieee_is_finite(behind))) then
            ok = .false.
            return
        end if
        if (ahead <= 0.0_real64 .or. behind <= 0.0_real64) then
            ok = .false.
            return
        end if

        ! d(ln K)/dh over the step, converted to a pressure derivative by the
        ! same specific weight the head was formed with.
        slope = (log(ahead) - log(behind))/(2.0_real64*step)
        dlog_dp = slope/SPECIFIC_WEIGHT
    end subroutine eval_relative_conductivity_tangent

    !> \(K_{fl h}\) of one material at a pore pressure, m/s.
    !>
    !> The model reads the head off the state itself; the effective suction is
    !> deliberately left unset, so `calc_head_hcf` forms \(h=p/(\rho_{\mathrm{std}}g)\),
    !> which is the isothermal unfrozen meaning of this solver's pressure.
    function conductivity_value(slot, pressure) result(value)
        implicit none
        integer(int32), intent(in) :: slot
        real(real64), intent(in) :: pressure
        real(real64) :: value

        type(type_state) :: state

        value = 0.0_real64
        call state%pressure%set(pressure)
        call conductivity_laws(slot)%curve%p%calc_Kflh(state, value)
    end function conductivity_value

    !> Forget any mobility fault a previous trial recorded.
    subroutine clear_mobility_fault()
        implicit none

        mobility_fault = .false.
        mobility_fault_slot = 0
        mobility_fault_pressure = 0.0_real64
        mobility_fault_value = 0.0_real64
    end subroutine clear_mobility_fault

    !> Record that a cell asked for a mobility the model cannot supply.
    !>
    !> Called from a pointwise kernel, which has no error argument and no
    !> communicator: the kernel writes a nonfinite residual so the step is
    !> refused, and this is what lets the refusal name its cause.
    subroutine note_mobility_fault(slot, pressure)
        implicit none
        !> The material slot of the offending cell.
        integer(int32), intent(in) :: slot
        !> The pressure it was evaluated at, Pa.
        real(real64), intent(in) :: pressure

        if (mobility_fault) return
        mobility_fault = .true.
        mobility_fault_slot = slot
        mobility_fault_pressure = pressure

        ! The raw model output is the informative number: it is what the law
        ! answered, before any normalisation this module applies to it.
        mobility_fault_value = 0.0_real64
        if (slot < 1 .or. slot > num_conductivity_laws) return
        if (conductivity_laws(slot)%kind /= MIXED_CONDUCTIVITY_HCF) return
        mobility_fault_value = conductivity_value(slot, pressure)
    end subroutine note_mobility_fault

    !> Whether a kernel refused a mobility since the last clear.
    logical function mobility_fault_recorded()
        implicit none

        mobility_fault_recorded = mobility_fault
    end function mobility_fault_recorded

    !> What the recorded mobility fault was.
    subroutine mobility_fault_report(slot, pressure, conductivity)
        implicit none
        !> The material slot of the offending cell, zero when none was recorded.
        integer(int32), intent(inout) :: slot
        !> The pressure it was evaluated at, Pa.
        real(real64), intent(inout) :: pressure
        !> The conductivity the model answered with, m/s; zero for a fixed law.
        real(real64), intent(inout) :: conductivity

        slot = mobility_fault_slot
        pressure = mobility_fault_pressure
        conductivity = mobility_fault_value
    end subroutine mobility_fault_report

end module app_mixed_hydraulic_closure
