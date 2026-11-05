module io_input_output_conditions
    use, intrinsic :: iso_fortran_env
    use :: core_parallel_mpi
    use :: stdlib_strings, only:to_string, strip, ends_with
    use :: stdlib_logger
    use :: json_module, only:json_file
    use :: module_core
    use :: io_input_base, only:get_json_value, abst_input
    implicit none
    private

    public :: type_output_settings

    character(*), parameter :: file_format = "file_format"
    character(*), parameter :: unit = "unit"
    character(*), parameter :: value = "value"
    character(*), parameter :: output_interval = "output_interval"
    character(*), parameter :: valid_time_units(5) = [character(len=16) :: "second", "minute", "hour", "day", "year"]
    character(*), parameter :: variables = "variables"
    character(*), parameter :: valid_variables_lists(10) = [character(32) :: &
                                                            ! Keep this list in sync with OUTPUT_VARIABLE_TYPES in include/constants.fypp
                                                            'water_content', 'ice_content', 'vapor_content', &
                                                            'temperature', 'thermal_conductivity', 'volumetric_heat_capacity', &
                                                            'pressure', 'water_flux', 'hydraulic_conductivity', &
                                                            'vapor_flux']

    !!------------------------------------------------------------------------------------------------------------------------------
    type :: type_field_output
        character(:), allocatable :: file_format
        character(:), allocatable :: output_time_unit
        real(real64) :: output_time_offset = 0.0d0
        character(:), allocatable :: output_interval_unit
        real(real64) :: output_interval_step
        character(:), allocatable :: variable_names(:)
    contains
        procedure, pass(self) :: display => display_output_settings_fields
    end type type_field_output

    interface
        module subroutine display_output_settings_fields(self)
            implicit none
            class(type_field_output), intent(in) :: self
        end subroutine display_output_settings_fields
    end interface

    !!------------------------------------------------------------------------------------------------------------------------------
    type :: types_history_output
        character(:), allocatable :: file_format
        character(:), allocatable :: observation_type
        character(:), allocatable :: output_time_unit
        real(real64) :: output_time_offset = 0.0d0
        character(:), allocatable :: output_interval_unit
        real(real64) :: output_interval_step
        character(:), allocatable :: variable_names(:)
        integer(int32) :: num_observations
        type(type_coordinate_dp), allocatable :: coordinates(:)
        integer(int32), allocatable :: node_ids(:)
    contains
        procedure, pass(self) :: display => display_output_settings_history
    end type types_history_output

    interface

        module subroutine display_output_settings_history(self)
            implicit none
            class(types_history_output), intent(in) :: self
        end subroutine display_output_settings_history
    end interface
    !!------------------------------------------------------------------------------------------------------------------------------
    type :: type_standard_output
        logical :: print_progress
        character(:), allocatable :: print_progress_unit
        real(real64) :: print_progress_interval
    contains
        procedure, pass(self) :: display => display_output_settings_standard
    end type type_standard_output
    interface
        module subroutine display_output_settings_standard(self)
            implicit none
            class(type_standard_output), intent(in) :: self
        end subroutine display_output_settings_standard
    end interface
    !!------------------------------------------------------------------------------------------------------------------------------
    type :: type_output_settings
        class(abst_input), pointer :: parent => null()
        character(:), allocatable :: file_name
        type(type_field_output) :: field_output
        type(types_history_output) :: history_output
        type(type_standard_output) :: standard_output
    contains
        procedure, pass(self) :: initialize => initialize_type_output_settings
        procedure, pass(self) :: display => display_output_settings
    end type type_output_settings

    interface
        module subroutine read_output_settings_fields(self, json)
            implicit none
            class(type_output_settings), intent(inout) :: self
            type(json_file), intent(inout) :: json
        end subroutine read_output_settings_fields

        module subroutine read_output_settings_history(self, json)
            implicit none
            class(type_output_settings), intent(inout) :: self
            type(json_file), intent(inout) :: json
        end subroutine read_output_settings_history

        module subroutine read_output_settings_standard(self, json)
            implicit none
            class(type_output_settings), intent(inout) :: self
            type(json_file), intent(inout) :: json
        end subroutine read_output_settings_standard
    end interface

contains

    subroutine initialize_type_output_settings(self)
        implicit none
        class(type_output_settings), intent(inout) :: self
        type(json_file) :: json

        call json%initialize()

        call json%load(filename=self%file_name)
        call json%print_error_message(output_unit)

        call read_output_settings_fields(self, json)
        call read_output_settings_history(self, json)
        call read_output_settings_standard(self, json)

        call json%destroy()
        call json%print_error_message(output_unit)

    end subroutine initialize_type_output_settings

    subroutine display_output_settings(self)
        implicit none
        class(type_output_settings), intent(in) :: self

        integer(int32) :: ierr, myrank

        call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
        if (myrank == 0) then
            write (*, '(A)') "=== Output Settings ==="
            write (*, '(A)') "--- Field Output ---"
            call self%field_output%display()
            write (*, '(A)') "--- History Output ---"
            call self%history_output%display()
            write (*, '(A)') "--- Standard Output ---"
            call self%standard_output%display()
        end if
    end subroutine display_output_settings

end module io_input_output_conditions
