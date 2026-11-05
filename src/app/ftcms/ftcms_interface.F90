#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscsnes.h>

!> The FTCMS engine: one PETSc system built from the project input.
!>
!> The input decides everything the system is. Which governing equations are
!> solved comes from `analysis_controls`, the space each one is discretised in
!> from `discretization_settings`, the materials from `Basic.json`, the
!> boundary and initial conditions and the time window from `Conditions.json`,
!> and what is written from `Output.json`. The engine turns the discretisation
!> into a field table, and the table alone decides the fields of the DM, the
!> pointwise kernels, the role of every boundary condition, the field numbers
!> the kernels read, and what is written.
!>
!> ## Field table
!>
!> Equations enter in `PHYSICS_TYPES` order. A `lagrange` equation contributes
!> one field, its potential; a `mixed` one contributes its flux and then its
!> potential. A combination without kernels is refused by name at setup, never
!> replaced by another one.
!>
!> ## Time integration
!>
!> Backward Euler, with the storage terms differenced at the quadrature points
!> from the accepted state held in auxiliary fields. The step length follows
!> the iteration-count controller of `Conditions.json`: a refused step is
!> retried shorter from the accepted state.
!>
!> ## Ownership
!>
!> The closure, the boundary registry and the nodal flux are module state,
!> because PETSc pointwise functions carry no typed context. One engine may be
!> live per rank.
module app_ftcms
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: petscdmplex
    use :: petscsnes
    use :: module_core, only:PHYSICS_TYPES
    use :: module_input, only:type_input
    use :: app_ftcms_kernels, only:NUM_CONSTANTS
    use :: app_ftcms_boundary_registry, only:MAX_BOUNDARY_SLOTS
    use :: governing_atmosphere, only:type_assimilation_controller
    implicit none
    private

    public :: type_ftcms

    !> What a field of the table carries for its equation.
    integer(int32), parameter, public :: ROLE_POTENTIAL = 1
    integer(int32), parameter, public :: ROLE_FLUX = 2

    !> How many fields one table may hold: a flux and a potential per equation.
    integer(int32), parameter :: MAX_FIELDS = 2*4

    !> How many observation points one run may write.
    integer(int32), parameter, public :: MAX_OBSERVATIONS = 64

    !> One field of the system.
    type :: type_field_entry
        !> `PHYSICS_TYPES` id of the equation the field belongs to.
        integer(int32) :: physics_id = 0
        !> `ROLE_POTENTIAL` or `ROLE_FLUX`.
        integer(int32) :: role = 0
        !> The name PETSc and the output give the field.
        character(len=32) :: name = ""
        !> Zero-based field number in the DM.
        PetscInt :: number = -1
    end type type_field_entry

    !> How the length of the next step is chosen.
    type :: type_step_control
        !> The length the first step is attempted at, the floor and the cap, s.
        real(real64) :: initial = 300.0d0
        real(real64) :: minimum = 1.0d0
        real(real64) :: maximum = 30000.0d0
        !> Iteration counts a step is kept between.
        integer(int32) :: iterations_low = 5
        integer(int32) :: iterations_high = 16
        !> Factors on an easy step, a laboured step and a refusal.
        real(real64) :: grow = 1.25d0
        real(real64) :: shrink = 0.8d0
        real(real64) :: retry = 0.5d0
    end type type_step_control

    !> The local-error control of `Conditions.json`.
    type :: type_error_control
        logical :: active = .false.
        real(real64) :: relative_tolerance = 1.0d-2
        real(real64) :: temperature_tolerance = 1.0d-3
        real(real64) :: pressure_tolerance = 10.0d0
        real(real64) :: proportional_gain = 0.2d0
        real(real64) :: integral_gain = 0.15d0
        real(real64) :: safety = 0.9d0
        real(real64) :: max_growth = 2.0d0
        !> The error of the last accepted step, and whether the last attempt was refused.
        real(real64) :: previous_error = 1.0d0
        logical :: after_rejection = .false.
    end type type_error_control

    !> One observation point and where it lies in the local mesh.
    type :: type_observation
        !> Physical coordinates (x, y, z) as `Output.json` gives them, m.
        real(real64) :: coordinate(3) = 0.0d0
        !> The same point in the two intrinsic coordinates, m.
        real(real64) :: plane(2) = 0.0d0
        !> The local cell holding the point, -1 when another rank holds it.
        PetscInt :: cell = -1
    end type type_observation

    !> One time series of `Output.json` written at the observation points.
    type :: type_history_series
        !> The `Output.json` variable name.
        character(len=32) :: variable = ""
        !> Fortran unit of its file on rank 0, -1 elsewhere.
        integer(int32) :: unit = -1
        !> Values per observation: 1 for a scalar, 3 for a vector.
        integer(int32) :: width = 1
    end type type_history_series

    type :: type_ftcms
        !> The project input; the engine reads it and does not own it.
        type(type_input), pointer :: input => null()

        ! ---- Problem ----
        !> Which governing equations are solved, by `PHYSICS_TYPES` id.
        logical :: is_active(PHYSICS_TYPES%NUM_ID) = .false.
        !> The discretisation of each equation, by `PHYSICS_TYPES` id.
        integer(int32) :: space(PHYSICS_TYPES%NUM_ID) = 0
        integer(int32) :: degree(PHYSICS_TYPES%NUM_ID) = 0
        integer(int32) :: evaluation(PHYSICS_TYPES%NUM_ID) = 0
        !> `COMP_TYPES` id of the analysis.
        integer(int32) :: computation_type = 0
        !> The field table.
        type(type_field_entry) :: fields(MAX_FIELDS)
        integer(int32) :: num_fields = 0
        !> Zero-based field numbers the kernels read, -1 for an absent field.
        PetscInt :: temperature_field = -1
        PetscInt :: flux_field = -1
        PetscInt :: pressure_field = -1
        !> "Cell Sets" value of each closure slot.
        integer(int32), allocatable :: material_ids(:)

        ! ---- Scales and physics ----
        real(real64) :: flux_scale = 1.0d0
        real(real64) :: pressure_scale = 1.0d0
        real(real64) :: length_scale = 1.0d0
        real(real64) :: temperature_scale = 1.0d0
        real(real64) :: thermal_scale = 1.0d0
        real(real64) :: density = 1000.0d0
        real(real64) :: gravity(2) = 0.0d0
        !> The potential of an equation the run does not solve, C and Pa.
        real(real64) :: fixed_temperature = 0.0d0
        real(real64) :: fixed_pressure = 0.0d0
        real(real64) :: initial_temperature = 0.0d0
        real(real64) :: initial_pressure = 0.0d0

        ! ---- PETSc objects ----
        DM :: dm
        DM :: aux_dm
        PetscDS :: ds
        Vec :: solution
        Vec :: old_solution
        Vec :: local_solution
        Vec :: auxiliary
        Vec :: accepted_auxiliary
        Vec :: boundary_auxiliary(MAX_BOUNDARY_SLOTS)
        integer(int32) :: num_boundary_auxiliary = 0
        SNES :: snes
        !> Seepage faces: their slots and "Face Sets" values, and the variable
        !> bounds the solve holds their pressure between.
        integer(int32) :: seepage_slots(MAX_BOUNDARY_SLOTS) = 0
        integer(int32) :: seepage_ids(MAX_BOUNDARY_SLOTS) = 0
        integer(int32) :: num_seepage = 0
        Vec :: lower_bound
        Vec :: upper_bound
        PetscScalar :: constants(NUM_CONSTANTS) = 0.0
        logical :: nodal_flux = .false.

        ! ---- Time ----
        !> Time of the accepted state, s.
        real(real64) :: time = 0.0d0
        real(real64) :: start_time = 0.0d0
        real(real64) :: end_time = 0.0d0
        !> Length of the last attempted step and the next one, s.
        real(real64) :: dt = 0.0d0
        real(real64) :: next_step = 0.0d0
        type(type_step_control) :: step_control
        type(type_error_control) :: error_control
        !> The BDF order asked for in `Basic.json`, and the weights of the step being solved.
        integer(int32) :: target_order = 1
        integer(int32) :: order = 1
        real(real64) :: bdf_weights(3) = [1.0d0, -1.0d0, 0.0d0]
        !> Accepted solutions before the current one, and the steps between them, s:
        !> `old_solution` is y_n, `older_solution` y_(n-1), `oldest_solution` y_(n-2).
        Vec :: older_solution
        Vec :: oldest_solution
        real(real64) :: history_steps(2) = 0.0d0
        integer(int32) :: history_depth = 0
        !> Per-dof absolute tolerance of the error estimate, and 1 on the potentials.
        Vec :: error_tolerance
        Vec :: error_mask
        integer(int32) :: accepted_steps = 0
        integer(int32) :: rejected_steps = 0
        integer(int32) :: total_iterations = 0

        ! ---- Output ----
        character(len=:), allocatable :: output_directory
        !> Field output: interval, s (zero for none), next time, count, format.
        real(real64) :: field_interval = 0.0d0
        real(real64) :: next_field_time = 0.0d0
        integer(int32) :: field_count = 0
        character(len=8) :: field_format = ""
        character(len=32), allocatable :: field_variables(:)
        !> Seconds per output time unit.
        real(real64) :: field_time_unit = 1.0d0
        DM :: output_dm
        logical :: has_output_dm = .false.
        !> Time of the last field file and the output times of all of them.
        real(real64) :: last_field_time = 0.0d0
        real(real64), allocatable :: collection_times(:)
        !> History output at the observation points.
        real(real64) :: history_interval = 0.0d0
        real(real64) :: next_history_time = 0.0d0
        real(real64) :: history_time_unit = 1.0d0
        real(real64) :: last_history_time = 0.0d0
        logical :: history_written = .false.
        type(type_observation) :: observations(MAX_OBSERVATIONS)
        integer(int32) :: num_observations = 0
        type(type_history_series), allocatable :: series(:)
        !> Seconds between progress lines; zero for none.
        real(real64) :: progress_interval = 0.0d0
        real(real64) :: next_progress_time = 0.0d0

        ! ---- Assimilation ----
        type(type_assimilation_controller) :: assimilation
        logical :: assimilation_enabled = .false.
        !> The "Face Sets" value the atmosphere exchanges with the soil across.
        integer(int32) :: surface_boundary = -1

        logical :: ready = .false.
    contains
        ! ---- Lifecycle ----
        procedure, public, pass(self) :: initialize => initialize_ftcms
        procedure, public, pass(self) :: destroy => destroy_ftcms
        ! ---- Algorithm ----
        procedure, public, pass(self) :: run => run_ftcms
    end type type_ftcms

    interface
        ! ---- ftcms_setup ----
        module subroutine initialize_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout), target :: self
        end subroutine initialize_ftcms

        module subroutine destroy_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self
        end subroutine destroy_ftcms

        module subroutine store_constants(self, ierr)
            implicit none
            class(type_ftcms), intent(inout) :: self
            PetscErrorCode, intent(inout) :: ierr
        end subroutine store_constants

        ! ---- ftcms_step ----
        module subroutine run_ftcms(self)
            implicit none
            class(type_ftcms), intent(inout) :: self
        end subroutine run_ftcms

        module subroutine publish_accepted_state(self, ierr)
            implicit none
            class(type_ftcms), intent(inout) :: self
            PetscErrorCode, intent(inout) :: ierr
        end subroutine publish_accepted_state

        ! ---- ftcms_assimilation ----
        module subroutine initialize_assimilation(self)
            implicit none
            class(type_ftcms), intent(inout) :: self
        end subroutine initialize_assimilation

        module subroutine assimilate(self)
            implicit none
            class(type_ftcms), intent(inout) :: self
        end subroutine assimilate

        ! ---- ftcms_output ----
        module subroutine initialize_output(self)
            implicit none
            class(type_ftcms), intent(inout) :: self
        end subroutine initialize_output

        module subroutine write_output(self, force)
            implicit none
            class(type_ftcms), intent(inout) :: self
            !> Write both outputs now, whatever their schedule says.
            logical, intent(in) :: force
        end subroutine write_output

        module subroutine destroy_output(self)
            implicit none
            class(type_ftcms), intent(inout) :: self
        end subroutine destroy_output
    end interface

contains

    !> Whether this rank holds no point of "Face Sets" value `value`.
    !>
    !> Asked before `DMGetStratumIS`: on a rank without such points the IS it
    !> returns is not one the Fortran null object compares equal to.
    function stratum_is_empty(dm, value) result(empty)
        implicit none
        DM, intent(in) :: dm
        integer(int32), intent(in) :: value
        logical :: empty

        PetscInt :: count
        PetscErrorCode :: ierr

        call DMGetStratumSize(dm, "Face Sets", value, count, ierr)
        if (ierr /= 0) error stop "PETSc could not count a boundary"
        empty = count == 0
    end function stratum_is_empty

end module app_ftcms
