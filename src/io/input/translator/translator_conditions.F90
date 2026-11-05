submodule(io_input_translator) translator_conditions
    implicit none
contains

    module subroutine execute_condition_boundary(self, input, target_physics, configs)
        implicit none
        class(type_input_translator), intent(in) :: self
        class(type_input), intent(in) :: input
        type(type_constant_id), intent(in) :: target_physics
        type(type_config_bc), intent(inout), allocatable :: configs(:)

        integer(int32) :: i
        integer(int32) :: index, num_active

        ! --- Pass 1: Count active boundaries ---
        num_active = 0
        do index = 1, input%conditions%num_boundaries
            associate (physics_data => input%conditions%boundary_conditions(index)%physics(target_physics%ID))
                if (physics_data%is_active) num_active = num_active + 1
            end associate
        end do

        ! Allocate compact array for active boundaries only
        if (allocated(configs)) deallocate (configs)
        allocate (configs(num_active))

        if (num_active == 0) return

        ! --- Pass 2: Convert and store only active boundaries ---
        num_active = 0
        do index = 1, input%conditions%num_boundaries
            associate (physics_data => input%conditions%boundary_conditions(index)%physics(target_physics%ID))

                if (.not. physics_data%is_active) cycle

                num_active = num_active + 1

                ! Store JSON id for mapping to mesh entity_id
                configs(num_active)%boundary_id = input%conditions%boundary_conditions(index)%id
                configs(num_active)%physics_type = target_physics
                configs(num_active)%num_time_points = physics_data%num_time_points

                select case (target_physics%ID)
                case (PHYSICS_TYPES%THERMAL%ID)
                    configs(num_active)%bc_kind = THERMAL_BC_TYPES%to_object(physics_data%bc_type)
                case (PHYSICS_TYPES%HYDRAULIC%ID)
                    configs(num_active)%bc_kind = HYDRAULIC_BC_TYPES%to_object(physics_data%bc_type)
                end select

                configs(num_active)%bc_data_kind = BC_DATA_PROVIDERS%to_object(physics_data%bc_value_type)

                select case (configs(num_active)%bc_kind%ID)
                case (THERMAL_BC_TYPES%DIRICHLET%ID, THERMAL_BC_TYPES%FLUX%ID, &
                      HYDRAULIC_BC_TYPES%DIRICHLET%ID, HYDRAULIC_BC_TYPES%FLUX%ID)
                    configs(num_active)%num_variables = 1
                    call allocate_array(configs(num_active)%time_points, &
                                        configs(num_active)%num_time_points)
                    call allocate_array(configs(num_active)%values, &
                                        configs(num_active)%num_variables, &
                                        configs(num_active)%num_time_points)

                    do i = 1, configs(num_active)%num_time_points
                        configs(num_active)%time_points(i) = physics_data%values(i)%time
                        configs(num_active)%values(1, i) = physics_data%values(i)%value
                    end do

                case (THERMAL_BC_TYPES%ROBIN%ID, THERMAL_BC_TYPES%CONVECTIVE%ID, THERMAL_BC_TYPES%RADIATION%ID)
                    ! Two values, or five for a Robin coefficient that switches with the surface value.
                    configs(num_active)%num_variables = 2
                    if (allocated(physics_data%values(1)%values)) then
                        if (size(physics_data%values(1)%values) >= 5) configs(num_active)%num_variables = 5
                    end if
                    call allocate_array(configs(num_active)%time_points, &
                                        configs(num_active)%num_time_points)
                    call allocate_array(configs(num_active)%values, &
                                        configs(num_active)%num_variables, &
                                        configs(num_active)%num_time_points)

                    do i = 1, configs(num_active)%num_time_points
                        configs(num_active)%time_points(i) = physics_data%values(i)%time

                        if (.not. allocated(physics_data%values(i)%values)) then
                            error stop 'Boundary condition requires two values, but none were provided.'
                        end if
                        if (size(physics_data%values(i)%values) < 2) then
                            error stop 'Boundary condition requires at least two values.'
                        end if

                        if (size(physics_data%values(i)%values) < configs(num_active)%num_variables) then
                            error stop 'Every time of a switching Robin condition needs all five values.'
                        end if
                        configs(num_active)%values(:, i) = &
                            physics_data%values(i)%values(1:configs(num_active)%num_variables)
                    end do

                case (HYDRAULIC_BC_TYPES%SEEPAGE%ID)
                    configs(num_active)%num_variables = 3
                    call allocate_array(configs(num_active)%time_points, &
                                        configs(num_active)%num_time_points)
                    call allocate_array(configs(num_active)%values, &
                                        configs(num_active)%num_variables, &
                                        configs(num_active)%num_time_points)

                    do i = 1, configs(num_active)%num_time_points
                        configs(num_active)%time_points(i) = physics_data%values(i)%time

                        if (.not. allocated(physics_data%values(i)%values)) then
                            error stop 'Seepage BC requires three values [q_pot m/s, Pmin Pa, Pmax Pa].'
                        end if
                        if (size(physics_data%values(i)%values) < 3) then
                            error stop 'Seepage BC requires at least three values [q_pot m/s, Pmin Pa, Pmax Pa].'
                        end if

                        configs(num_active)%values(1, i) = physics_data%values(i)%values(1)
                        configs(num_active)%values(2, i) = physics_data%values(i)%values(2)
                        configs(num_active)%values(3, i) = physics_data%values(i)%values(3)
                    end do

                case (THERMAL_BC_TYPES%ATMOSPHERIC%ID, HYDRAULIC_BC_TYPES%ATMOSPHERIC%ID)
                    ! Atmospheric BC: initialized to zeros; DA controller updates each cycle.
                    configs(num_active)%num_variables = 3
                    configs(num_active)%bc_data_kind = BC_DATA_PROVIDERS%CONSTANT
                    configs(num_active)%num_time_points = 1
                    call allocate_array(configs(num_active)%time_points, 1)
                    call allocate_array(configs(num_active)%values, 3, 1)
                    configs(num_active)%time_points(1) = 0.0d0
                    configs(num_active)%values(:, 1) = 0.0d0

                case default
                    configs(num_active)%num_variables = 0
                end select

            end associate
        end do

    end subroutine execute_condition_boundary

    module subroutine execute_condition_initial(self, input, ic_target, config)
        implicit none
        class(type_input_translator), intent(in) :: self
        class(type_input), intent(in) :: input
        type(type_constant_id), intent(in) :: ic_target
        type(type_config_ic), intent(inout) :: config

        if (.not. IC_TARGETS%is_valid(ic_target)) then
            call config%reset()
            return
        end if

        config%physics_type = ic_target
        associate (condition_data => input%conditions%initial_conditions%physics(ic_target%ID))
            if (ic_target == IC_TARGETS%POROSITY) then
                config%active = .true. ! Porosity IC might always be needed
            else
                config%active = input%basic%analysis_controls%is_active(ic_target%ID)
            end if
            if (.not. config%active) return
            
            config%ic_kind = IC_METHODS%to_object(condition_data%type)
            if (config%ic_kind == IC_METHODS%UNIFORM) then
                config%value = condition_data%value
            else if (config%ic_kind == IC_METHODS%FROM_FILE) then
                ! The mesh reader rejects this earlier: PETSc reads the mesh
                ! from the Gmsh file, not the node data written beside it.
                error stop "Initial conditions cannot be read from the mesh file"
            end if
        end associate

    end subroutine execute_condition_initial

    module subroutine execute_condition_time(self, input, config)
        implicit none
        class(type_input_translator), intent(in) :: self
        class(type_input), intent(in) :: input
        type(type_config_time), intent(inout) :: config

        config%target_bdf_order = input%basic%solver_settings%bdf_order

        associate (time_control => input%conditions%time_control)
            config%time_stepping_unit = TIME_UNITS%to_object(time_control%time_stepping%unit)
            config%initial_step = time_control%time_stepping%initial_step * config%time_stepping_unit%value

            config%simulation_period_unit = TIME_UNITS%to_object(time_control%simulation_period%unit)
            config%start_time = time_control%simulation_period%start * config%simulation_period_unit%value
            config%end_time = time_control%simulation_period%end * config%simulation_period_unit%value

        end associate

    end subroutine execute_condition_time

    module subroutine execute_condition_time_ats(self, input, config)
        implicit none
        class(type_input_translator), intent(in) :: self
        class(type_input), intent(in) :: input
        type(type_config_time_ats), intent(inout) :: config

        type(type_constant_value) :: time_unit
        real(real64) :: time_conv_coeff

        associate (time_control => input%conditions%time_control)
            config%active = time_control%adaptive_stepping%is_active

            if (config%active) then
                call config%reset()
                config%active = .true.
                config%iter_min = time_control%adaptive_stepping%iter_min
                config%iter_max = time_control%adaptive_stepping%iter_max
                config%scale_up = time_control%adaptive_stepping%scale_up
                config%scale_down = time_control%adaptive_stepping%scale_down
                config%scale_retry = time_control%adaptive_stepping%scale_retry
                config%use_error_control = time_control%adaptive_stepping%use_error_control
                config%error_rtol = time_control%adaptive_stepping%error_relative_tolerance
                config%error_atol_temperature = &
                    time_control%adaptive_stepping%error_absolute_tolerance_temperature
                config%error_atol_pressure = time_control%adaptive_stepping%error_absolute_tolerance_pressure
                config%pi_k_p = time_control%adaptive_stepping%proportional_gain
                config%pi_k_i = time_control%adaptive_stepping%integral_gain
                config%safety_factor = time_control%adaptive_stepping%safety_factor
                config%max_growth_rate = time_control%adaptive_stepping%max_growth_rate
                config%max_dT_per_step = time_control%adaptive_stepping%max_temperature_change_per_step
                config%max_relative_change_per_step = time_control%adaptive_stepping%max_relative_change_per_step

                time_unit = TIME_UNITS%to_object(time_control%time_stepping%unit)
                time_conv_coeff = time_unit%value
                config%dt_min = time_control%time_stepping%min_step * time_conv_coeff
                config%dt_max = time_control%time_stepping%max_step * time_conv_coeff
            else
                call config%reset()
            end if
        end associate

    end subroutine execute_condition_time_ats
end submodule translator_conditions
