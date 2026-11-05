module io_input_conditions
    use, intrinsic :: iso_fortran_env
    use :: core_parallel_mpi
    use :: stdlib_strings, only:to_string, strip, ends_with
    use :: stdlib_logger
    use :: json_module, only:json_file
    use :: module_core
    use :: io_input_base, only:get_json_value, abst_input
    implicit none
    private

    public :: type_conditions

    character(*), parameter :: thermal = "thermal"
    character(*), parameter :: hydraulic = "hydraulic"
    character(*), parameter :: mechanical = "mechanical"
    character(*), parameter :: porosity = "porosity"

    type :: type_time_controls_simulation_period
        character(:), allocatable :: unit
        real(real64) :: start
        real(real64) :: end
    end type type_time_controls_simulation_period

    type :: type_time_controls_time_stepping
        character(:), allocatable :: unit
        real(real64) :: initial_step
        real(real64) :: min_step
        real(real64) :: max_step
    end type type_time_controls_time_stepping

    type :: type_time_controls_ats
        logical :: is_active
        integer(int32) :: iter_min
        integer(int32) :: iter_max
        real(real64) :: scale_up
        real(real64) :: scale_down
        real(real64) :: scale_retry
        logical :: use_error_control = .true.
        real(real64) :: error_relative_tolerance = 1.0d-2
        real(real64) :: error_absolute_tolerance_temperature = 1.0d-3
        real(real64) :: error_absolute_tolerance_pressure = 1.0d1
        real(real64) :: proportional_gain = 0.10d0
        real(real64) :: integral_gain = 0.15d0
        real(real64) :: safety_factor = 0.9d0
        real(real64) :: max_growth_rate = 2.0d0
        real(real64) :: max_temperature_change_per_step = 5.0d0
        real(real64) :: max_relative_change_per_step = 0.3d0
    end type type_time_controls_ats

    type :: type_time_controls
        class(type_conditions), pointer :: parent => null()
        type(type_time_controls_simulation_period) :: simulation_period
        type(type_time_controls_time_stepping) :: time_stepping
        type(type_time_controls_ats) :: adaptive_stepping
        ! real(real64), allocatable :: boundary_time_points(:)
    contains
        procedure, pass(self) :: display => display_time_controls
    end type type_time_controls

    interface
        module subroutine display_time_controls(self)
            implicit none
            class(type_time_controls), intent(in) :: self
        end subroutine display_time_controls
    end interface

    !!------------------------------------------------------------------------------------------------------------------------------
    type :: type_boundary_local_time_dependent
        real(real64) :: time
        character(:), allocatable :: time_iso
        real(real64) :: value
        real(real64), allocatable :: values(:)
    end type type_boundary_local_time_dependent

    type :: type_boundary_local
        logical :: is_active = .false.
        character(:), allocatable :: bc_type
        character(:), allocatable :: bc_value_type
        integer(int32) :: num_time_points
        type(type_boundary_local_time_dependent), allocatable :: values(:)
    contains
        procedure, public, pass(self) :: read => read_conditions_bc_local
    end type type_boundary_local

    interface
        module subroutine read_conditions_bc_local(self, json, buffer_in, end_index, physics_type)
            implicit none
            class(type_boundary_local), intent(inout) :: self
            type(json_file), intent(inout) :: json
            character(*), intent(in) :: buffer_in(:)
            integer(int32), intent(in) :: end_index
            type(type_constant_id), intent(in) :: physics_type

        end subroutine read_conditions_bc_local
    end interface

    type :: type_boundary_conditions
        class(type_conditions), pointer :: parent => null()
        integer(int32) :: id
        type(type_boundary_local) :: physics(PHYSICS_TYPES%NUM_ID)
    contains
        procedure, pass(self) :: display => display_boundary_conditions
    end type type_boundary_conditions

    interface
        module subroutine display_boundary_conditions(self)
            implicit none
            class(type_boundary_conditions), intent(in) :: self
        end subroutine display_boundary_conditions
    end interface
    !!------------------------------------------------------------------------------------------------------------------------------
    type, extends(type_boundary_local) :: type_boundary_local_initial
        integer(int32) :: id
    end type

    type :: type_initial_local
        character(:), allocatable :: type
        real(real64) :: value
        type(type_boundary_local_initial), allocatable :: boundary(:)
        character(:), allocatable :: field_name
    end type type_initial_local

    type :: type_initial_conditions
        class(type_conditions), pointer :: parent => null()
        type(type_initial_local) :: physics(IC_TARGETS%NUM_ID)
    contains
        procedure, pass(self) :: display => display_initial_conditions
    end type type_initial_conditions

    interface
        module subroutine display_initial_conditions(self)
            implicit none
            class(type_initial_conditions), intent(in) :: self
        end subroutine display_initial_conditions
    end interface
    !!------------------------------------------------------------------------------------------------------------------------------
    type :: type_conditions
        class(abst_input), pointer :: parent => null()
        character(:), allocatable :: file_name
        type(type_time_controls) :: time_control
        type(type_boundary_conditions), allocatable :: boundary_conditions(:)
        integer(int32) :: num_boundaries
        type(type_initial_conditions) :: initial_conditions
    contains
        procedure, pass(self), public :: initialize => initialize_type_conditions
        procedure, pass(self), private :: read_time_controls => read_conditions_time_controls
        procedure, pass(self), private :: read_boundary_conditions => read_conditions_bc
        procedure, pass(self), private :: read_initial_conditions => read_conditions_initial_conditions
        procedure, pass(self), public :: display => display_conditions
    end type type_conditions

    interface
        module subroutine initialize_type_conditions(self)
            implicit none
            class(type_conditions), intent(inout), target :: self

        end subroutine initialize_type_conditions

        module subroutine read_conditions_time_controls(self, json)
            implicit none
            class(type_conditions), intent(inout) :: self
            type(json_file), intent(inout) :: json
        end subroutine read_conditions_time_controls

        module subroutine read_conditions_bc(self, json)
            implicit none
            class(type_conditions), intent(inout) :: self
            type(json_file), intent(inout) :: json

        end subroutine read_conditions_bc

        module subroutine read_conditions_initial_conditions(self, json)
            implicit none
            class(type_conditions), intent(inout) :: self
            type(json_file), intent(inout) :: json
        end subroutine read_conditions_initial_conditions

        module subroutine display_conditions(self)
            implicit none
            class(type_conditions), intent(in) :: self
        end subroutine display_conditions

    end interface

contains

end module io_input_conditions
