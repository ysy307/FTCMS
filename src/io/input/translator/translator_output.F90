submodule(io_input_translator) translator_output
    implicit none
contains

    module subroutine execute_output_manager(self, input, output_type, config)
        implicit none
        class(type_input_translator), intent(in) :: self
        class(type_input), intent(in) :: input
        type(type_constant_id), intent(in) :: output_type
        type(type_config_output_manager), intent(inout) :: config

        select case (output_type%ID)
        case (OUTPUT_TYPES%FIELD%ID)
            associate (output => input%output_settings%field_output)
                config%interval_val = output%output_interval_step
                config%file_format = FILE_FORMATS%to_object(output%file_format)

                if (config%file_format == FILE_FORMATS%NONE) then
                    config%interval_val = 0.0d0
                    config%interval_unit = TIME_UNITS%to_object("second")
                    config%output_unit = TIME_UNITS%to_object("second")
                    return
                end if

                config%interval_unit = TIME_UNITS%to_object(output%output_interval_unit)
                config%output_unit = TIME_UNITS%to_object(output%output_time_unit)
                config%output_time_offset = output%output_time_offset
            end associate

        case (OUTPUT_TYPES%HISTORY%ID)
            associate (output => input%output_settings%history_output)
                config%interval_val = output%output_interval_step
                config%file_format = FILE_FORMATS%to_object(output%file_format)

                if (config%file_format == FILE_FORMATS%NONE) then
                    config%interval_val = 0.0d0
                    config%interval_unit = TIME_UNITS%to_object("second")
                    config%output_unit = TIME_UNITS%to_object("second")
                    config%output_time_offset = 0.0d0
                    return
                end if

                config%interval_unit = TIME_UNITS%to_object(output%output_interval_unit)
                config%output_unit = TIME_UNITS%to_object(output%output_time_unit)
                config%output_time_offset = output%output_time_offset
            end associate

        case default
            error stop "Only field and history output have output settings"
        end select

    end subroutine execute_output_manager

end submodule translator_output
