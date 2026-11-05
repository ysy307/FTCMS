module core_types
    use :: core_types_config_base, only: &
        abst_config
    use :: types_config_conditions, only: &
        type_config_bc, &
        type_config_ic
    use :: types_config_control, only: &
        type_config_output_manager, &
        type_config_time, &
        type_config_time_ats
    use :: types_config_physics, only: &
        abst_config_physics_model, &
        abst_config_physics_material, &
        type_config_wrf, &
        type_config_hcf, &
        type_config_gcc, &
        type_config_constitutive

    use :: core_types_geometry_coordinate, only: &
        type_coordinate_dp, &
        type_coordinate_int
    use :: core_types_geometry_coordinate_array, only: &
        type_coordinate_array_dp, &
        type_coordinate_array_int

    use :: core_types_physics_state, only: &
        type_state
    use :: core_types_physics_state_thc, only: &
        type_state_thc
    use :: core_types_physics_meteorology, only: &
        type_meteorology, &
        type_meteorology_data
    use :: core_types_physics_solar, only: &
        type_solar_system

    use :: core_types_utils_datetime, only: &
        type_datetime
    implicit none
    public

end module core_types
