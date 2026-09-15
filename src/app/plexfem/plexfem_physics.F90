#include <petsc/finclude/petscsys.h>

!> State the PetscDS pointwise functions read, and the material response they
!> evaluate through it.
!>
!> A pointwise function is called by PETSc through a C function pointer, so it
!> cannot carry a typed context: PETSc offers only the PetscScalar array of
!> PetscDSSetConstants and a void pointer that Fortran cannot dereference
!> safely. The governing model therefore lives here, as module state owned by
!> the single type_plexfem instance of the run.
!>
!> This is sound because the instance is unique per rank and PETSc assembles
!> the FEM residual serially within a rank. It is not a place to cache anything
!> that changes during an assembly.
module app_plexfem_physics
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_loc
    use :: petscsys
    use :: module_core, only:type_state
    use :: module_governing, only:type_thermal, type_hydraulic
    implicit none
    private

    public :: PLEXFEM_CONSTANTS
    public :: C_POROSITY, C_PRESSURE, C_DT, C_MEASURING_MATERIAL
    public :: C_BDF_ORDER, C_BDF_0, bdf_order, bdf_coefficient
    public :: A_MATERIAL_ID, A_BC_KIND, A_BC_VALUE_1, A_BC_VALUE_2
    public :: A_DARCY, A_COUPLING, A_GRAVITY, A_POTENTIAL
    public :: A_TEMPERATURE_OLD, A_PRESSURE_OLD, A_TEMPERATURE_OLDER, A_PRESSURE_OLDER
    public :: PLEXFEM_AUX_FIELDS
    public :: BC_NONE, BC_FLUX, BC_ROBIN
    public :: attach_thermal_model, detach_thermal_model
    public :: attach_hydraulic_model, hydraulic_is_active
    public :: eval_water_content, eval_water_capacity, eval_water_diffusion, eval_water_coupling
    public :: eval_enthalpy_pressure_tangent, eval_water_temperature_tangent
    public :: eval_phase_contents, eval_water_flux_tangent
    public :: eval_gravity_flux, eval_gravity_flux_tangent
    public :: eval_effective_density, eval_water_driver, eval_liquid_pressure
    public :: eval_capacity, eval_conductivity, eval_enthalpy
    public :: eval_conductivity_tangent, eval_conductivity_pressure_tangent
    public :: prescribed_value, set_dirichlet_value, dirichlet_context

    !> Layout of the PetscDSSetConstants array. The kernels index it by name so
    !> that adding a constant cannot silently shift another one.
    integer(int32), parameter :: C_POROSITY = 1 !< porosity [-]
    integer(int32), parameter :: C_PRESSURE = 2 !< pore pressure [Pa]
    integer(int32), parameter :: C_DT = 3       !< length of the current step [s]
    !> Zero while the run is solving. While the lumped nodal volumes are being
    !> measured it holds the material being measured, and every kernel but the
    !> indicator contributes nothing, so that what comes back is the volume and
    !> not the volume plus a boundary flux.
    integer(int32), parameter :: C_MEASURING_MATERIAL = 4
    !> The BDF coefficients of the step being solved, and the order they belong
    !> to. The storage term is integrated at the quadrature points, so the kernel
    !> that differences the history needs them where it runs.
    integer(int32), parameter :: C_BDF_ORDER = 5
    integer(int32), parameter :: C_BDF_0 = 6
    integer(int32), parameter :: C_BDF_1 = 7
    integer(int32), parameter :: C_BDF_2 = 8
    integer(int32), parameter :: PLEXFEM_CONSTANTS = 8

    !> Layout of the auxiliary fields, in the order they are added to the
    !> auxiliary DM. What varies over the mesh lives here rather than in the
    !> constants, which are one set of numbers for the whole domain.
    integer(int32), parameter :: A_MATERIAL_ID = 1 !< material of the cell
    integer(int32), parameter :: A_BC_KIND = 2     !< which natural condition this boundary is
    integer(int32), parameter :: A_BC_VALUE_1 = 3  !< influx, or the transfer coefficient
    integer(int32), parameter :: A_BC_VALUE_2 = 4  !< ambient value of a Robin condition
    !> Element-representative transport coefficients of the water equation, one
    !> value per cell, refreshed from the current iterate before each assembly.
    !>
    !> The Darcy conductivity spans orders of magnitude inside a single element
    !> at a freezing front -- the ice impedance is 10^(-7Q) and Q crosses from
    !> zero to one across it -- so evaluating it at each quadrature point lets
    !> the warm points of an element decide how much water it passes. The older
    !> path averages the nodal values over the element instead, geometrically,
    !> and records that quadrature-point evaluation makes an element straddling
    !> the front far too permeable. These carry that average.
    integer(int32), parameter :: A_DARCY = 5     !< D_HH, element geometric mean
    integer(int32), parameter :: A_COUPLING = 6  !< D_HT, element signed geometric mean
    integer(int32), parameter :: A_GRAVITY = 7   !< gravity flux, element geometric mean
    !> The liquid pressure the Darcy flux is driven by, p_l = -psi_eff, carried
    !> as a field of its own so that its gradient is the gradient of its nodal
    !> interpolant.
    !>
    !> Degree one, unlike the others. What is wanted is grad I_h[p_l(T_i, p_i)],
    !> and a pointwise function cannot form it from grad T and grad p: the chain
    !> rule gives grad p_l(T_h, p_h) instead, and the two differ wherever p_l is
    !> a strongly nonlinear function of the unknowns -- which at a freezing front
    !> is everywhere, because the effective suction turns over across it. The
    !> older path interpolates the nodal potential, and this is how that reaches
    !> a kernel.
    integer(int32), parameter :: A_POTENTIAL = 8
    !> The last two accepted states, so that the storage term can be integrated
    !> at the quadrature points rather than lumped at the nodes.
    !>
    !> A pointwise function is handed only the current iterate, so the history the
    !> time derivative differences against has to arrive as fields PETSc
    !> interpolates for it. Degree one, like the solution they came from.
    integer(int32), parameter :: A_TEMPERATURE_OLD = 9
    integer(int32), parameter :: A_PRESSURE_OLD = 10
    integer(int32), parameter :: A_TEMPERATURE_OLDER = 11
    integer(int32), parameter :: A_PRESSURE_OLDER = 12
    integer(int32), parameter :: PLEXFEM_AUX_FIELDS = 12

    !> What a boundary's auxiliary values mean. A boundary PETSc integrates but
    !> that carries none of these is left at zero flux, which is what a natural
    !> boundary with no residual already is.
    integer(int32), parameter :: BC_NONE = 0
    integer(int32), parameter :: BC_FLUX = 1
    integer(int32), parameter :: BC_ROBIN = 2

    !> The governing model the kernels evaluate the material response with.
    !>
    !> Only its intent(in) bindings may be called from a kernel: the model also
    !> holds an enthalpy history cache that this path never fills, and calling
    !> anything that reads it would be reading uninitialised state.
    type(type_thermal), pointer, save :: thermal_model => null()

    !> The water model, when the project asks for it. A null pointer is what
    !> tells the kernels that pressure is a constant of the run rather than a
    !> field they solve for.
    type(type_hydraulic), pointer, save :: hydraulic_model => null()

    !> Reusable phase states, kept rather than rebuilt.
    !>
    !> The phase split is a pure function of temperature, material, porosity
    !> and pressure, so a state asked for twice can be answered from the first
    !> answer exactly.
    !>
    !> The number of them is the working set of one quadrature point, not of one
    !> residual evaluation. Assembling the Jacobian there visits the same point
    !> once per registered block, and the blocks whose tangents come by central
    !> difference ask for the state at T +/- h and at p +/- h as well, so the
    !> distinct states in play at a point are the base one and four neighbours.
    !> Two slots hold none of that: every block evicts what the last one just
    !> built, and the Jacobian was measured taking four times as long as the
    !> residual for that reason alone.
    !>
    !> Callers are handed the slot index rather than a copy: type_state carries
    !> allocatable components, and copying one costs more than it saves.
    !> Keeping the states here also spares every kernel call the construction
    !> and destruction of a local one, which is where most of the measured
    !> saving comes from: on _1Domain this takes the step from five seconds to
    !> under three.
    !>
    !> The key is compared for equality, not nearness: these are the same
    !> values passed twice, not two nearby ones.
    integer(int32), parameter :: STATE_SLOTS = 8
    type(type_state), save :: slot_state(STATE_SLOTS)
    real(real64), save :: slot_temperature(STATE_SLOTS) = 0.0d0
    real(real64), save :: slot_porosity(STATE_SLOTS) = 0.0d0
    real(real64), save :: slot_pressure(STATE_SLOTS) = 0.0d0
    integer(int32), save :: slot_material(STATE_SLOTS) = -1
    logical, save :: slot_filled(STATE_SLOTS) = .false.
    integer(int32), save :: slot_next = 1

    !> Step of the central difference that gives dlambda/dT.
    !>
    !> The conductivity's derivative has no closed form to reach from here: it
    !> would need the phase conductivities of the geometric mixing rule, which
    !> the materials keep to themselves. 1.0d-3 K is the same scale below which
    !> thermal_coefficients already refuses to divide by a temperature
    !> difference.
    real(real64), parameter :: TEMPERATURE_STEP = 1.0d-3

    !> Step of the central difference that gives dlambda/dp, for the same reason.
    !>
    !> 1 Pa: far below the pressures a freezing column develops, which run to
    !> 1e4 Pa and beyond, and far above the point where the difference of two
    !> conductivities loses its significant digits.
    real(real64), parameter :: PRESSURE_STEP = 1.0d0

    !> Prescribed temperature of the essential boundary, by "Face Sets" value.
    !> The essential-boundary callback PETSc calls has no constants argument at
    !> all, so this cannot travel with the others.
    integer(int32), parameter :: MAX_BOUNDARY_ID = 64
    integer(int32), parameter :: MAX_FIELDS = 2
    real(real64), save :: dirichlet_value(0:MAX_BOUNDARY_ID, 0:MAX_FIELDS - 1) = 0.0d0

    !> Which boundary a callback was registered for. PETSc hands the essential
    !> callback only the context pointer it was given, so each registration
    !> points at its own element here and one callback per field serves every
    !> boundary of that field.
    integer(c_int), save, target :: dirichlet_owner(0:MAX_BOUNDARY_ID, 0:MAX_FIELDS - 1) = 0

contains

    !> Point the kernels at the model. It has to outlive every assembly, so the
    !> caller owns it and keeps it alive for the run.
    subroutine attach_thermal_model(model)
        implicit none
        type(type_thermal), intent(in), target :: model

        if (associated(thermal_model)) then
            error stop "The PetscDS path already has a thermal model attached"
        end if
        thermal_model => model
        slot_filled(:) = .false.
    end subroutine attach_thermal_model

    subroutine detach_thermal_model()
        implicit none

        nullify (thermal_model)
        nullify (hydraulic_model)
        slot_filled(:) = .false.
    end subroutine detach_thermal_model

    !> Point the kernels at the water model as well, making pressure a field.
    subroutine attach_hydraulic_model(model)
        implicit none
        type(type_hydraulic), intent(in), target :: model

        if (associated(hydraulic_model)) then
            error stop "The PetscDS path already has a water model attached"
        end if
        hydraulic_model => model
        slot_filled(:) = .false.
    end subroutine attach_hydraulic_model

    !> Whether this run solves for pressure as well as temperature.
    logical function hydraulic_is_active()
        implicit none

        hydraulic_is_active = associated(hydraulic_model)
    end function hydraulic_is_active

    !> Total water content Theta at a quadrature point.
    !>
    !> This is what the water balance differences over a step, the way the
    !> energy balance differences the enthalpy: it is monotone in pressure and
    !> stays Lipschitz across the phase change, where its derivative does not.
    !> Ice counts through its density, because what the balance transports is
    !> water mass rather than volume.
    subroutine eval_water_content(temperature, pressure, material_id, constants, water_content)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: water_content

        integer(int32) :: which

        call ensure_state(temperature, pressure, material_id, constants, which)
        ! The water model's own definition of stored water, not a sum of phases
        ! rebuilt here: it is what eval_water_capacity is the derivative of, and
        ! it carries the compressive storage. Where the medium saturates
        ! dTheta_SWRC/dp vanishes and that term is the whole of the capacity, so
        ! a locally reassembled Theta would leave the storage diagonal
        ! corresponding to no residual term at all.
        call hydraulic_model%compute_total_storage(material_id, slot_state(which), water_content)
    end subroutine eval_water_content

    !> Liquid water and ice content at a node, separately.
    !>
    !> The two phases as the retention and freezing curves leave them, not the
    !> stored water eval_water_content differences: what is written out has to
    !> be the physical contents, because the validation forms the total itself
    !> as theta_w + (rho_i/rho_w) theta_i and the compressive storage carried by
    !> the mass balance is not part of that total.
    subroutine eval_phase_contents(temperature, pressure, material_id, constants, liquid, ice)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: liquid
        real(real64), intent(inout) :: ice

        integer(int32) :: which

        call ensure_state(temperature, pressure, material_id, constants, which)
        call slot_state(which)%water_content%get(liquid)
        call slot_state(which)%ice_content%get(ice)
    end subroutine eval_phase_contents

    !> Order the step is being taken at, as the kernels see it.
    pure function bdf_order(constants) result(order)
        implicit none
        PetscScalar, intent(in) :: constants(*)
        integer(int32) :: order

        order = int(nint(real(constants(C_BDF_ORDER), real64)), int32)
    end function bdf_order

    !> BDF coefficient of one history level, index zero being the iterate.
    pure function bdf_coefficient(constants, level) result(coefficient)
        implicit none
        PetscScalar, intent(in) :: constants(*)
        integer(int32), intent(in) :: level
        real(real64) :: coefficient

        coefficient = real(constants(C_BDF_0 + level), real64)
    end function bdf_coefficient

    !> Effective density: the water equation's conserved quantity in the units
    !> its tolerance is stated in.
    !>
    !> rho_w theta_w + rho_i theta_i + rho_w theta_v plus the compressive
    !> storage, in kg/m3. The acceptance criterion weighs the water residual
    !> against a tolerance in those units, so the quantity it is scored against
    !> has to be this one and not the volumetric content.
    subroutine eval_effective_density(temperature, pressure, material_id, constants, density)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: density

        integer(int32) :: which

        call ensure_state(temperature, pressure, material_id, constants, which)
        call hydraulic_model%calc_effective_density_value(material_id, slot_state(which), density)
    end subroutine eval_effective_density

    !> Storage capacity dTheta/dp at a quadrature point.
    subroutine eval_water_capacity(temperature, pressure, material_id, constants, capacity)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: capacity

        integer(int32) :: which

        call ensure_state(temperature, pressure, material_id, constants, which)
        call hydraulic_model%compute_C_eq(material_id, slot_state(which), capacity)
    end subroutine eval_water_capacity

    !> dU/dp: what a change of pore pressure does to the stored energy.
    !>
    !> The off-diagonal of the lumped storage block, not a small correction to
    !> it. Pressure moves water across the freezing curve, so this term carries
    !> latent heat: it is of order -rho_w L_f dTheta/dp, and deep in the
    !> pore-volume limit the phase split amplifies dtheta_w/dp by
    !> 1/(1 - rho_i/rho_w) ~ 12 and turns the sign over.
    subroutine eval_enthalpy_pressure_tangent(temperature, pressure, material_id, constants, tangent)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: tangent

        integer(int32) :: which

        call ensure_state(temperature, pressure, material_id, constants, which)
        call thermal_model%compute_coupling_mass_term(material_id, slot_state(which), tangent)
    end subroutine eval_enthalpy_pressure_tangent

    !> dTheta/dT: what a change of temperature does to the stored water.
    !>
    !> Zero in the closure that is active today, because that one takes the
    !> total water as theta_SWRC(p) with no explicit temperature dependence and
    !> lets ice take the complement. It is asked of the model rather than
    !> assumed, so that enabling a closure which does make Theta depend on T
    !> (PHASE_CLOSURE_PAINTER_KARRA in fusion.F90) reaches the Jacobian without
    !> anyone having to remember this block exists.
    subroutine eval_water_temperature_tangent(temperature, pressure, material_id, constants, tangent)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: tangent

        integer(int32) :: which

        call ensure_state(temperature, pressure, material_id, constants, which)
        call hydraulic_model%compute_coupling_mass_term(material_id, slot_state(which), tangent)
    end subroutine eval_water_temperature_tangent

    !> Darcy conductivity tensor D_HH at a quadrature point.
    subroutine eval_water_diffusion(temperature, pressure, material_id, constants, diffusion)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: diffusion(:, :)

        integer(int32) :: which

        call ensure_state(temperature, pressure, material_id, constants, which)
        call hydraulic_model%compute_diffusion_term(material_id, slot_state(which), diffusion)
        call span_coordinates(diffusion)
    end subroutine eval_water_diffusion

    !> Thermal coupling of the water flux, D_HT, at a quadrature point.
    subroutine eval_water_coupling(temperature, pressure, material_id, constants, coupling)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: coupling(:, :)

        integer(int32) :: which

        call ensure_state(temperature, pressure, material_id, constants, which)
        call hydraulic_model%compute_coupling_diffusion_term(material_id, slot_state(which), coupling)
        call span_coordinates(coupling)
    end subroutine eval_water_coupling

    !> The gravity-driven part of the water flux, as the component that belongs
    !> on the coordinate opposing gravity.
    !>
    !> Darcy's law is q = -K grad(p/(rho g) + z), so beside the pressure
    !> gradient the flux carries a term that is just the conductivity itself,
    !> acting downward. It is the reason an x-z column is meshed in x-z at all:
    !> without it the water has no reason to move down and the only transport is
    !> the suction gradient the freezing front creates.
    !>
    !> The model decides which of its own components is vertical from the
    !> project's computation type, and leaves the rest at zero, so the sum picks
    !> that component out whichever it is and returns zero for a plan-view
    !> analysis where gravity has no component in the mesh.
    subroutine eval_gravity_flux(temperature, pressure, material_id, constants, flux)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: flux

        real(real64) :: advective(3)
        integer(int32) :: which

        call ensure_state(temperature, pressure, material_id, constants, which)
        advective(:) = 0.0d0
        call hydraulic_model%compute_advective_term(material_id, slot_state(which), advective)
        ! The model returns the flux itself; the residual carries it on the
        ! other side, the same way compute_R2 is subtracted in the older path.
        flux = -sum(advective)
    end subroutine eval_gravity_flux

    !> Derivatives of the gravity flux in both unknowns.
    subroutine eval_gravity_flux_tangent(temperature, pressure, material_id, constants, dflux_dT, dflux_dp)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: dflux_dT, dflux_dp

        real(real64) :: plus, minus

        call eval_gravity_flux(temperature + TEMPERATURE_STEP, pressure, material_id, constants, plus)
        call eval_gravity_flux(temperature - TEMPERATURE_STEP, pressure, material_id, constants, minus)
        dflux_dT = (plus - minus) / (2.0d0 * TEMPERATURE_STEP)

        call eval_gravity_flux(temperature, pressure + PRESSURE_STEP, material_id, constants, plus)
        call eval_gravity_flux(temperature, pressure - PRESSURE_STEP, material_id, constants, minus)
        dflux_dp = (plus - minus) / (2.0d0 * PRESSURE_STEP)
    end subroutine eval_gravity_flux_tangent

    !> The liquid pressure at a node, p_l = -psi_eff.
    !>
    !> The potential the Darcy flux follows. The closure carries p_w as the
    !> total-water-equivalent pressure while the liquid sees the effective
    !> suction the phase split published, and freezing raises that suction: its
    !> gradient is what draws water to a front.
    subroutine eval_liquid_pressure(temperature, pressure, material_id, constants, liquid_pressure)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: liquid_pressure

        real(real64) :: suction
        integer(int32) :: which

        call ensure_state(temperature, pressure, material_id, constants, which)
        call slot_state(which)%effective_suction%get(suction)
        liquid_pressure = -suction
    end subroutine eval_liquid_pressure

    !> The coefficients the water flux is actually driven by.
    !>
    !> The liquid Darcy flux follows the gradient of the liquid pressure, and in
    !> a freezing soil that is not the primary unknown: the closure carries p_w
    !> as the total-water-equivalent pressure, while the liquid sees
    !> p_l = -psi_eff, the effective suction the phase split published. Freezing
    !> raises that suction, and its gradient is what draws water to a front --
    !> the whole mechanism the Mizoguchi column measures. Driving the flux with
    !> grad p_w instead leaves only the thermo-osmotic term, which is the
    !> temperature dependence of surface tension and two orders too weak.
    !>
    !> Written on the unknowns by the chain rule,
    !>
    !>   grad p_l = -(dSeff/dp) grad p - (dSeff/dT) grad T
    !>
    !> so what multiplies grad p is D (-dSeff/dp), and what multiplies grad T is
    !> D (-dSeff/dT) plus the thermo-osmotic coefficient, which is independent of
    !> it. Composing them here rather than at each kernel is what keeps the
    !> Jacobian consistent: the tangents below difference these very
    !> coefficients, so a change to the driver reaches them without being
    !> restated.
    subroutine eval_water_driver(temperature, pressure, material_id, constants, &
                                 darcy_element, coupling_element, driver_p, driver_T)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        !> The element's representative conductivities, read from the auxiliary
        !> field the assembly refreshed. The suction derivatives stay pointwise:
        !> what the element average is for is the coefficient, whose variation
        !> inside an element spans orders, and composing the two here keeps the
        !> product the flux is built from in one place.
        real(real64), intent(in) :: darcy_element, coupling_element
        real(real64), intent(inout) :: driver_p(:, :), driver_T(:, :)

        real(real64) :: dsuction_dp, dsuction_dT
        integer(int32) :: which, d

        driver_p(:, :) = 0.0d0
        driver_T(:, :) = 0.0d0
        do d = 1, size(driver_p, 1)
            driver_p(d, d) = darcy_element
            driver_T(d, d) = coupling_element
        end do

        call ensure_state(temperature, pressure, material_id, constants, which)
        call slot_state(which)%dSeff_dP%get(dsuction_dp)
        call slot_state(which)%dSeff_dT%get(dsuction_dT)

        driver_T = driver_T + driver_p * (-dsuction_dT)
        driver_p = driver_p * (-dsuction_dp)
    end subroutine eval_water_driver

    !> The water flux coefficients and their derivatives in both unknowns.
    !>
    !> The Darcy conductivity is not a mild function of the state. The ice
    !> impedance alone is 10^(-Omega Q) with Omega = 7, so
    !> dln(k_r)/dz = -7 ln(10) dQ/dz, and Q crosses from zero to one across the
    !> freezing band: within a single Newton step the conductivity moves by a
    !> factor of several. Holding it fixed makes the water block a Picard
    !> iteration whose contraction is set by exactly that movement, which is the
    !> few-percent-per-iteration decay the water residual showed while these
    !> terms were missing.
    subroutine eval_water_flux_tangent(temperature, pressure, material_id, constants, &
                                       darcy_element, coupling_element, &
                                       ddiffusion_dT, ddiffusion_dp, dcoupling_dT, dcoupling_dp)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        !> Held fixed while the suction derivatives are differenced: the element
        !> average is a property of the cell, not of the quadrature point, so the
        !> tangent this gives is the one of the coefficients frozen over the
        !> assembly -- a Picard term in the conductivity and a Newton term in the
        !> driver, which is the split the older path also makes.
        real(real64), intent(in) :: darcy_element, coupling_element
        real(real64), intent(inout) :: ddiffusion_dT(:, :), ddiffusion_dp(:, :)
        real(real64), intent(inout) :: dcoupling_dT(:, :), dcoupling_dp(:, :)

        real(real64) :: plus(size(ddiffusion_dT, 1), size(ddiffusion_dT, 2))
        real(real64) :: minus(size(ddiffusion_dT, 1), size(ddiffusion_dT, 2))
        real(real64) :: plus_T(size(ddiffusion_dT, 1), size(ddiffusion_dT, 2))
        real(real64) :: minus_T(size(ddiffusion_dT, 1), size(ddiffusion_dT, 2))

        call eval_water_driver(temperature + TEMPERATURE_STEP, pressure, material_id, constants, &
                               darcy_element, coupling_element, plus, plus_T)
        call eval_water_driver(temperature - TEMPERATURE_STEP, pressure, material_id, constants, &
                               darcy_element, coupling_element, minus, minus_T)
        ddiffusion_dT = (plus - minus) / (2.0d0 * TEMPERATURE_STEP)
        dcoupling_dT = (plus_T - minus_T) / (2.0d0 * TEMPERATURE_STEP)

        call eval_water_driver(temperature, pressure + PRESSURE_STEP, material_id, constants, &
                               darcy_element, coupling_element, plus, plus_T)
        call eval_water_driver(temperature, pressure - PRESSURE_STEP, material_id, constants, &
                               darcy_element, coupling_element, minus, minus_T)
        ddiffusion_dp = (plus - minus) / (2.0d0 * PRESSURE_STEP)
        dcoupling_dp = (plus_T - minus_T) / (2.0d0 * PRESSURE_STEP)
    end subroutine eval_water_flux_tangent

    !> Enthalpy density U at a quadrature point, latent heat included.
    !>
    !> This is what the residual differences over a step. U is monotone in
    !> temperature and stays Lipschitz across the phase change even where its
    !> derivative spikes, which is what makes the residual differentiable there
    !> and Newton able to converge on it.
    subroutine eval_enthalpy(temperature, pressure, material_id, constants, enthalpy)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: enthalpy

        integer(int32) :: which

        call ensure_state(temperature, pressure, material_id, constants, which)
        call thermal_model%calc_enthalpy_density(material_id, slot_state(which), enthalpy)
    end subroutine eval_enthalpy

    !> Apparent heat capacity C = dU/dT at a quadrature point.
    subroutine eval_capacity(temperature, pressure, material_id, constants, capacity)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: capacity

        integer(int32) :: which

        call ensure_state(temperature, pressure, material_id, constants, which)
        call thermal_model%compute_apparent_capacity(material_id, slot_state(which), capacity)
    end subroutine eval_capacity

    !> Thermal conductivity tensor at a quadrature point, laid out on the
    !> computation plane the governing model was built for.
    !> Re-express a governing model's tensor in the coordinate frame the
    !> kernels work in.
    !>
    !> PetscDS hands a pointwise function the *coordinate* dimension, and this
    !> path reads meshes with three coordinates so that an x-z column is the
    !> surface it actually is rather than a flat one. The governing models fill
    !> a tensor over the *computation* dimension instead -- two, for a plane --
    !> so on an x-z mesh they set the entries for x and y and leave the entry for
    !> z at zero. z is the direction everything in such a column varies in, so
    !> what was being solved had no vertical conduction and no vertical Darcy
    !> flow at all: heat reached the interior only sideways, and the water could
    !> not redistribute down the column.
    !>
    !> The models return an isotropic tensor -- one coefficient on the diagonal
    !> -- so the frame it is expressed in does not matter and the same
    !> coefficient belongs on every diagonal entry. Anything else cannot be
    !> carried across without knowing which physical axis each entry meant, so it
    !> is refused rather than silently reinterpreted.
    subroutine span_coordinates(tensor)
        implicit none
        real(real64), intent(inout) :: tensor(:, :)

        real(real64) :: coefficient
        integer(int32) :: d, e, n

        n = size(tensor, 1)
        coefficient = tensor(1, 1)
        do d = 1, n
            do e = 1, n
                if (d == e) then
                    ! Every entry the model did set must agree, or the tensor is
                    ! not the isotropic one this may re-express.
                    if (tensor(d, e) /= 0.0d0 .and. tensor(d, e) /= coefficient) then
                        error stop "The PetscDS path cannot embed an anisotropic material tensor"
                    end if
                else if (tensor(d, e) /= 0.0d0) then
                    error stop "The PetscDS path cannot embed a material tensor with off-diagonal terms"
                end if
            end do
        end do

        tensor(:, :) = 0.0d0
        do d = 1, n
            tensor(d, d) = coefficient
        end do
    end subroutine span_coordinates

    subroutine eval_conductivity(temperature, pressure, material_id, constants, conductivity)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: conductivity(:, :)

        integer(int32) :: which

        call ensure_state(temperature, pressure, material_id, constants, which)
        call thermal_model%compute_diffusion_term(material_id, slot_state(which), conductivity)
        call span_coordinates(conductivity)
    end subroutine eval_conductivity

    !> Thermal conductivity tensor and its temperature derivative.
    subroutine eval_conductivity_tangent(temperature, pressure, material_id, constants, conductivity, dconductivity_dT)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: conductivity(:, :)
        real(real64), intent(inout) :: dconductivity_dT(:, :)

        real(real64) :: conductivity_plus(size(conductivity, 1), size(conductivity, 2))
        real(real64) :: conductivity_minus(size(conductivity, 1), size(conductivity, 2))

        call eval_conductivity(temperature, pressure, material_id, constants, conductivity)
        call eval_conductivity(temperature + TEMPERATURE_STEP, pressure, material_id, constants, conductivity_plus)
        call eval_conductivity(temperature - TEMPERATURE_STEP, pressure, material_id, constants, conductivity_minus)

        dconductivity_dT = (conductivity_plus - conductivity_minus) / (2.0d0 * TEMPERATURE_STEP)
    end subroutine eval_conductivity_tangent

    !> Pressure derivative of the thermal conductivity tensor.
    !>
    !> Conduction depends on pore pressure because the composition does: the
    !> geometric mean is taken over the phase fractions, and pressure moves water
    !> between liquid and ice, whose conductivities differ by a factor of four.
    !> Small per unit of pressure, but pressure is measured in Pa while
    !> temperature is measured in K, so over the increments a Newton step takes
    !> this term is the larger of the two.
    subroutine eval_conductivity_pressure_tangent(temperature, pressure, material_id, constants, &
                                                  conductivity, dconductivity_dp)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        real(real64), intent(inout) :: conductivity(:, :)
        real(real64), intent(inout) :: dconductivity_dp(:, :)

        real(real64) :: conductivity_plus(size(conductivity, 1), size(conductivity, 2))
        real(real64) :: conductivity_minus(size(conductivity, 1), size(conductivity, 2))

        call eval_conductivity(temperature, pressure, material_id, constants, conductivity)
        call eval_conductivity(temperature, pressure + PRESSURE_STEP, material_id, constants, conductivity_plus)
        call eval_conductivity(temperature, pressure - PRESSURE_STEP, material_id, constants, conductivity_minus)

        dconductivity_dp = (conductivity_plus - conductivity_minus) / (2.0d0 * PRESSURE_STEP)
    end subroutine eval_conductivity_pressure_tangent

    !> The slot holding the state of this quadrature point, computing it only
    !> if no slot holds it already.
    !>
    !> Temperature is the unknown TS carries; the rest of what the freezing
    !> model needs is constant over this stage of the port and travels in the
    !> DS constants. No time history is set: TS owns the time derivative, so
    !> the kernels only ever ask for quantities of the current state.
    subroutine ensure_state(temperature, pressure, material_id, constants, which)
        implicit none
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(in) :: material_id
        PetscScalar, intent(in) :: constants(*)
        !> Slot the caller should read the state from.
        integer(int32), intent(inout) :: which

        real(real64) :: porosity
        integer(int32) :: k

        if (.not. associated(thermal_model)) then
            error stop "A PetscDS kernel was called before the thermal model was attached"
        end if

        porosity = real(constants(C_POROSITY), real64)

        do k = 1, STATE_SLOTS
            if (.not. slot_filled(k)) cycle
            if (slot_temperature(k) /= temperature) cycle
            if (slot_porosity(k) /= porosity) cycle
            if (slot_pressure(k) /= pressure) cycle
            if (slot_material(k) /= material_id) cycle
            which = k
            return
        end do

        which = slot_next
        slot_next = merge(1, slot_next + 1, slot_next >= STATE_SLOTS)

        call slot_state(which)%temperature%set(temperature)
        call slot_state(which)%pressure%set(pressure)
        call slot_state(which)%porosity%set(porosity)
        ! Full freezing physics: this path runs no homotopy continuation.
        call slot_state(which)%continuation_lambda%set(1.0d0)

        call thermal_model%update_water_phases(material_id, slot_state(which))

        slot_temperature(which) = temperature
        slot_porosity(which) = porosity
        slot_pressure(which) = pressure
        slot_material(which) = material_id
        slot_filled(which) = .true.
    end subroutine ensure_state

    !> Record the value the essential boundary of this field prescribes.
    subroutine set_dirichlet_value(field, boundary_id, value)
        implicit none
        integer(c_int), intent(in) :: field
        integer(int32), intent(in) :: boundary_id
        real(real64), intent(in) :: value

        call check_boundary_id(boundary_id)
        call check_field(field)
        dirichlet_value(boundary_id, field) = value
        dirichlet_owner(boundary_id, field) = int(boundary_id, c_int)
    end subroutine set_dirichlet_value

    !> The value the essential boundary of this field prescribes.
    real(real64) function prescribed_value(field, boundary_id)
        implicit none
        integer(c_int), intent(in) :: field
        integer(c_int), intent(in) :: boundary_id

        prescribed_value = dirichlet_value(boundary_id, field)
    end function prescribed_value

    !> The context to register alongside the essential callback of this
    !> boundary. Reading it back is what tells the callback which value to
    !> prescribe.
    function dirichlet_context(field, boundary_id) result(ctx)
        implicit none
        integer(c_int), intent(in) :: field
        integer(int32), intent(in) :: boundary_id
        type(c_ptr) :: ctx

        call check_boundary_id(boundary_id)
        call check_field(field)
        ctx = c_loc(dirichlet_owner(boundary_id, field))
    end function dirichlet_context

    subroutine check_field(field)
        implicit none
        integer(c_int), intent(in) :: field

        if (field < 0_c_int .or. field >= int(MAX_FIELDS, c_int)) then
            error stop "A field index is outside the range the PetscDS path can hold"
        end if
    end subroutine check_field

    subroutine check_boundary_id(boundary_id)
        implicit none
        integer(int32), intent(in) :: boundary_id

        if (boundary_id < 0 .or. boundary_id > MAX_BOUNDARY_ID) then
            error stop "A boundary id is outside the range the PetscDS path can hold"
        end if
    end subroutine check_boundary_id

end module app_plexfem_physics
