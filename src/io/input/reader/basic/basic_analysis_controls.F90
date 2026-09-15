submodule(io_input_basic) input_basic_analysis_controls
    implicit none
    !!------------------------------------------------------------------------------------------------------------------------------
    ! JSON key names for analysis controls
    !!------------------------------------------------------------------------------------------------------------------------------
    character(*), parameter :: analysis_controls = "analysis_controls"
    character(*), parameter :: calculate_thermal = "calculate_thermal"
    character(*), parameter :: calculate_hydraulic = "calculate_hydraulic"
    character(*), parameter :: calculate_mechanical = "calculate_mechanical"
    character(*), parameter :: enable_vapor_transport = "enable_vapor_transport"
    character(*), parameter :: coupling_mode = "coupling_mode"
    character(*), parameter :: coupling_modes_strings(2) = [character(len=16) :: "staggered", "monolithic"]
    !!------------------------------------------------------------------------------------------------------------------------------
contains
    module subroutine read_analysis_controls(self, json)
        implicit none
        class(type_input_basic), intent(inout) :: self
        type(json_file), intent(inout) :: json

        character(256) :: buffer(2) = [character(256) :: analysis_controls, ""]

        buffer(2) = calculate_thermal
        call get_json_value(json, join(buffer), self%analysis_controls%is_active(PHYSICS_TYPES%THERMAL%ID), &
                            is_required=.true., default_value=.false.)

        buffer(2) = calculate_hydraulic
        call get_json_value(json, join(buffer), self%analysis_controls%is_active(PHYSICS_TYPES%HYDRAULIC%ID), &
                            is_required=.true., default_value=.false.)

        buffer(2) = calculate_mechanical
        call get_json_value(json, join(buffer), self%analysis_controls%is_active(PHYSICS_TYPES%MECHANICAL%ID), &
                            is_required=.true., default_value=.false.)

        buffer(2) = enable_vapor_transport
        call get_json_value(json, join(buffer), self%analysis_controls%enable_vapor_transport, &
                    is_required=.false., default_value=.true.)

        if (.not. any(self%analysis_controls%is_active(:))) then
            call raise_error(ERROR_CODES%VAR_INVALID, opt=analysis_controls)
        end if

        buffer(2) = coupling_mode
        call get_json_value(json, join(buffer), self%analysis_controls%coupling_mode, &
                            is_required=.true., default_value="monolithic", valid_list=coupling_modes_strings)


    end subroutine read_analysis_controls

    module subroutine display_analysis_controls(self)
        implicit none
        class(type_analysis_controls) :: self

        write (*, '(a)') "Calculate Thermal: "//to_string(self%is_active(PHYSICS_TYPES%THERMAL%ID))
        write (*, '(a)') "Calculate Hydraulic: "//to_string(self%is_active(PHYSICS_TYPES%HYDRAULIC%ID))
        write (*, '(a)') "Calculate Mechanical: "//to_string(self%is_active(PHYSICS_TYPES%MECHANICAL%ID))
        write (*, '(a)') "Enable Vapor Transport: "//to_string(self%enable_vapor_transport)
        write (*, '(a)') "Coupling Mode: "//strip(self%coupling_mode)

    end subroutine display_analysis_controls
end submodule input_basic_analysis_controls
