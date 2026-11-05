module types_config_control
    use :: types_config_control_output_manager, only: &
        type_config_output_manager
    use :: types_config_control_time, only: &
        type_config_time, &
        type_config_time_ats
    implicit none
    private

    public :: type_config_output_manager

    public :: type_config_time
    public :: type_config_time_ats

end module types_config_control
