submodule(io_input_basic) input_basic_solver_settings
    implicit none

    !!------------------------------------------------------------------------------------------------------------------------------
    ! JSON key names for solver settings
    !!------------------------------------------------------------------------------------------------------------------------------
    character(*), parameter :: solver_settings = "solver_settings"
    character(*), parameter :: bdf_order = "bdf_order"
    character(*), parameter :: t_nonlinear_solver = "nonlinear_solver"
    character(*), parameter :: method = "method"
    character(*), parameter :: update_frequency = "update_frequency"
    character(*), parameter :: max_iterations = "max_iterations"
    character(*), parameter :: convergence = "convergence"
    character(*), parameter :: use_criteria = "use_criteria"
    character(*), parameter :: norm_type = "norm_type"
    character(*), parameter :: logic_between_criteria = "logic_between_criteria"
    character(*), parameter :: residual = "residual"
    character(*), parameter :: update = "update"
    character(*), parameter :: criteria = "criteria"
    character(*), parameter :: logic = "logic"
    character(*), parameter :: absolute_tolerance = "absolute_tolerance"
    character(*), parameter :: relative_tolerance = "relative_tolerance"
    character(*), parameter :: key_atol_enthalpy = "atol_enthalpy"
    character(*), parameter :: key_atol_density = "atol_density"
    character(*), parameter :: key_rtol = "rtol"
    character(*), parameter :: key_residual_eps = "residual_eps"
    character(*), parameter :: linear_solver = "linear_solver"
    character(*), parameter :: iterative_solver = "iterative_solver"
    character(*), parameter :: tolerance = "tolerance"

    !!------------------------------------------------------------------------------------------------------------------------------
    ! Valid parameter values (for settings not covered by core_constants)
    !!------------------------------------------------------------------------------------------------------------------------------
    ! Linear Solver Methods
    character(len=16), parameter :: LINEAR_METHOD_DIRECT = "direct"
    character(len=16), parameter :: LINEAR_METHOD_ITERATIVE = "iterative"
    character(len=16), parameter :: valid_linear_solver_methods(2) = [LINEAR_METHOD_DIRECT, LINEAR_METHOD_ITERATIVE]

    ! Parallel Schedule Types
    character(len=16), parameter :: SCHEDULE_TYPE_AFFINITY = "affinity"
    character(len=16), parameter :: SCHEDULE_TYPE_AUTO = "auto"
    character(len=16), parameter :: SCHEDULE_TYPE_DYNAMIC = "dynamic"
    character(len=16), parameter :: SCHEDULE_TYPE_GUIDED = "guided"
    character(len=16), parameter :: SCHEDULE_TYPE_RUNTIME = "runtime"
    character(len=16), parameter :: SCHEDULE_TYPE_STATIC = "static"
    character(len=16), parameter :: valid_schedule_types(6) = &
                                    [SCHEDULE_TYPE_AFFINITY, SCHEDULE_TYPE_AUTO, SCHEDULE_TYPE_DYNAMIC, &
                                     SCHEDULE_TYPE_GUIDED, SCHEDULE_TYPE_RUNTIME, SCHEDULE_TYPE_STATIC]

    ! String lists for JSON validation
    character(len=16), parameter :: valid_nonlinear_solver_methods_str(5) = &
                                    [character(len=16) :: "none", "newton", "modified_newton", "picard", "modified_picard"]
    character(len=16), parameter :: valid_norm_types_str(3) = [character(len=4) :: "l1", "l2", "linf"]
    character(len=16), parameter :: valid_criteria_types_str(5) = &
                                    [character(len=16) :: "none", "residual", "update", "both", "conserved"]
    character(len=16), parameter :: valid_logic_types_str(2) = [character(len=16) :: "and", "or"]
    character(len=16), parameter :: valid_local_criteria_types_str(3) = [character(len=16) :: "absolute", "relative", "both"]

contains

    !======================================================================
    ! Read Routines for Solver Settings
    !======================================================================

    module subroutine read_solver_settings(self, json)
        !> Reads all solver settings from the JSON object. This is the main entry point.
        implicit none
        class(type_input_basic), intent(inout) :: self
        type(json_file), intent(inout) :: json
        character(256) :: buffer(2)

        buffer(1) = solver_settings
        buffer(2) = bdf_order
        call get_json_value(json, join(buffer), self%solver_settings%bdf_order, &
                            is_required=.true., default_value=1, valid_range=[1, 6])

        call read_solver_settings_nonlinear(self, json)
        call read_solver_settings_linear(self, json)

    end subroutine read_solver_settings

    subroutine read_solver_settings_nonlinear(self, json)
        !> Reads nonlinear solver settings.
        implicit none
        class(type_input_basic), intent(inout) :: self
        type(json_file), intent(inout) :: json
        character(256) :: buffer(6)
        character(:), allocatable :: temp_string
        real(real64) :: default_atol_enthalpy, default_atol_density
        real(real64) :: default_rtol_conserved, default_residual_eps

        buffer(1) = solver_settings
        buffer(2) = t_nonlinear_solver

        buffer(3) = method
        call get_json_value(json, join(buffer(1:3)), temp_string, &
                            is_required=.false., default_value="picard", &
                            valid_list=valid_nonlinear_solver_methods_str)
        self%solver_settings%nonlinear_solver%method = NONLINEAR_SOLVER%to_id(temp_string)

        buffer(3) = update_frequency
        call get_json_value(json, join(buffer(1:3)), &
                            self%solver_settings%nonlinear_solver%update_frequency, &
                            is_required=.false., default_value=1, valid_range=[1, huge(1)])

        buffer(3) = max_iterations
        call get_json_value(json, join(buffer(1:3)), self%solver_settings%nonlinear_solver%max_iterations, &
                            is_required=.true., default_value=1000, valid_range=[1, huge(1)])

        buffer(3) = convergence
        buffer(4) = use_criteria
        call get_json_value(json, join(buffer(1:4)), temp_string, &
                            is_required=.true., default_value="both", valid_list=valid_criteria_types_str)
        self%solver_settings%nonlinear_solver%convergence%use_criteria = NONLINEAR_NORM_CRITERIA%to_id(temp_string)

        if (self%solver_settings%nonlinear_solver%convergence%use_criteria == NONLINEAR_NORM_CRITERIA%BOTH%ID) then
            buffer(4) = logic_between_criteria
            call get_json_value(json, join(buffer(1:4)), temp_string, &
                                is_required=.true., default_value="and", valid_list=valid_logic_types_str)
            self%solver_settings%nonlinear_solver%convergence%use_logic = NONLINEAR_LOGIC%to_id(temp_string)
        end if

        buffer(4) = norm_type
        call get_json_value(json, join(buffer(1:4)), temp_string, &
                            is_required=.true., default_value="l2", valid_list=valid_norm_types_str)
        self%solver_settings%nonlinear_solver%convergence%norm_type = NORM_TYPES%to_id(temp_string)

        if (any(self%solver_settings%nonlinear_solver%convergence%use_criteria == &
                [NONLINEAR_NORM_CRITERIA%RESIDUAL%ID, NONLINEAR_NORM_CRITERIA%BOTH%ID])) then
            buffer(4) = residual
            call read_solver_settings_nonlinear_convergence( &
                self%solver_settings%nonlinear_solver%convergence%residual, json, buffer, 4)
        end if
        if (any(self%solver_settings%nonlinear_solver%convergence%use_criteria == &
                [NONLINEAR_NORM_CRITERIA%UPDATE%ID, NONLINEAR_NORM_CRITERIA%BOTH%ID])) then
            buffer(4) = update
            call read_solver_settings_nonlinear_convergence( &
                self%solver_settings%nonlinear_solver%convergence%update, json, buffer, 4)
        end if

        ! Conserved-quantity convergence (PDF 6.2.4): weighted-RMS norm of the
        ! nodal enthalpy and effective-density increments plus a per-block
        ! residual ratio.
        if (self%solver_settings%nonlinear_solver%convergence%use_criteria == &
            NONLINEAR_NORM_CRITERIA%CONSERVED%ID) then
            ! Legacy flat keys (backward compatible)
            buffer(4) = key_atol_enthalpy
            call get_json_value(json, join(buffer(1:4)), &
                                self%solver_settings%nonlinear_solver%convergence%atol_enthalpy, &
                                is_required=.false., default_value=1.0d2, valid_range=[0.0d0, huge(0.0d0)])
            buffer(4) = key_atol_density
            call get_json_value(json, join(buffer(1:4)), &
                                self%solver_settings%nonlinear_solver%convergence%atol_density, &
                                is_required=.false., default_value=1.0d-3, valid_range=[0.0d0, huge(0.0d0)])
            buffer(4) = key_rtol
            call get_json_value(json, join(buffer(1:4)), &
                                self%solver_settings%nonlinear_solver%convergence%rtol_conserved, &
                                is_required=.false., default_value=1.0d-3, valid_range=[0.0d0, huge(0.0d0)])
            buffer(4) = key_residual_eps
            call get_json_value(json, join(buffer(1:4)), &
                                self%solver_settings%nonlinear_solver%convergence%residual_eps, &
                                is_required=.false., default_value=1.0d-1, valid_range=[0.0d0, huge(0.0d0)])

            ! Hierarchical keys (preferred)
            ! Keep defaults in distinct variables. Passing the destination as
            ! its own default aliases an INTENT(OUT) argument in json-fortran
            ! and can replace a valid flat-key value with zero.
            default_atol_enthalpy = self%solver_settings%nonlinear_solver%convergence%atol_enthalpy
            default_atol_density = self%solver_settings%nonlinear_solver%convergence%atol_density
            default_rtol_conserved = self%solver_settings%nonlinear_solver%convergence%rtol_conserved
            default_residual_eps = self%solver_settings%nonlinear_solver%convergence%residual_eps
            buffer(4) = "conserved"
            buffer(5) = "thermal"
            buffer(6) = absolute_tolerance
            call get_json_value(json, join(buffer(1:6)), &
                                self%solver_settings%nonlinear_solver%convergence%atol_enthalpy, &
                                is_required=.false., &
                                default_value=default_atol_enthalpy, &
                                valid_range=[0.0d0, huge(0.0d0)])

            buffer(5) = "hydraulic"
            buffer(6) = absolute_tolerance
            call get_json_value(json, join(buffer(1:6)), &
                                self%solver_settings%nonlinear_solver%convergence%atol_density, &
                                is_required=.false., &
                                default_value=default_atol_density, &
                                valid_range=[0.0d0, huge(0.0d0)])

            buffer(5) = ""
            buffer(6) = ""
            buffer(5) = "relative_tolerance"
            call get_json_value(json, join(buffer(1:5)), &
                                self%solver_settings%nonlinear_solver%convergence%rtol_conserved, &
                                is_required=.false., &
                                default_value=default_rtol_conserved, &
                                valid_range=[0.0d0, huge(0.0d0)])

            buffer(5) = "residual_ratio_tolerance"
            call get_json_value(json, join(buffer(1:5)), &
                                self%solver_settings%nonlinear_solver%convergence%residual_eps, &
                                is_required=.false., &
                                default_value=default_residual_eps, &
                                valid_range=[0.0d0, huge(0.0d0)])
        end if

        if (allocated(temp_string)) deallocate (temp_string)
    end subroutine read_solver_settings_nonlinear

    !> Reads a specific convergence criteria block (residual or update).
    subroutine read_solver_settings_nonlinear_convergence(convergence_obj, json, buffer, end_index)
        implicit none
        type(type_convergence_criteria), intent(inout) :: convergence_obj
        type(json_file), intent(inout) :: json
        character(*), intent(in) :: buffer(:)
        integer(int32), intent(in) :: end_index
        character(:), allocatable :: temp_string
        character(len=256), allocatable :: local_buffer(:)

        allocate (local_buffer(size(buffer) + 1))
        local_buffer(1:end_index) = buffer(1:end_index)

        local_buffer(end_index + 1) = criteria
        call get_json_value(json, join(local_buffer), temp_string, &
                            is_required=.true., default_value="both", valid_list=valid_local_criteria_types_str)
        convergence_obj%criteria = NONLINEAR_CRITERIA%to_id(temp_string)

        if (convergence_obj%criteria == NONLINEAR_CRITERIA%BOTH%ID) then
            local_buffer(end_index + 1) = logic
            call get_json_value(json, join(local_buffer), temp_string, &
                                is_required=.true., default_value="and", valid_list=valid_logic_types_str)
            convergence_obj%logic = NONLINEAR_LOGIC%to_id(temp_string)
        end if

        if (any(convergence_obj%criteria == [NONLINEAR_CRITERIA%ABSOLUTE%ID, NONLINEAR_CRITERIA%BOTH%ID])) then
            local_buffer(end_index + 1) = absolute_tolerance
            call get_json_value(json, join(local_buffer), convergence_obj%absolute_tolerance, &
                                is_required=.true., default_value=1.0d-6, valid_range=[0.0d0, huge(0.0d0)])
        end if
        if (any(convergence_obj%criteria == [NONLINEAR_CRITERIA%RELATIVE%ID, NONLINEAR_CRITERIA%BOTH%ID])) then
            local_buffer(end_index + 1) = relative_tolerance
            call get_json_value(json, join(local_buffer), convergence_obj%relative_tolerance, &
                                is_required=.true., default_value=1.0d-6, valid_range=[0.0d0, huge(0.0d0)])
        end if

        if (allocated(temp_string)) deallocate (temp_string)
        if (allocated(local_buffer)) deallocate (local_buffer)
    end subroutine read_solver_settings_nonlinear_convergence

    subroutine read_solver_settings_linear(self, json)
        !> Reads linear solver settings for all active physics.
        implicit none
        class(type_input_basic), intent(inout) :: self
        type(json_file), intent(inout) :: json
        character(256) :: buffer(3)
        character(256) :: legacy_buffer(2)
        character(:), allocatable :: temp_string
        integer(int32) :: legacy_int
        real(real64) :: legacy_real
        logical :: legacy_found

        buffer(1) = solver_settings
        buffer(2) = linear_solver

        buffer(3) = max_iterations
        call get_json_value(json, join(buffer(1:3)), self%solver_settings%linear_solver%max_iterations, &
                    is_required=.false., default_value=10000, valid_range=[1, huge(1)])
        buffer(3) = tolerance
        call get_json_value(json, join(buffer(1:3)), self%solver_settings%linear_solver%tolerance, &
                    is_required=.false., default_value=1.0d-6, valid_range=[0.0d0, huge(0.0d0)])

        buffer(3) = "petsc_options"
        call get_json_value(json, join(buffer(1:3)), temp_string, &
                    is_required=.false., default_value="")
        self%solver_settings%linear_solver%petsc_options = temp_string

        ! `ksp_type` / `pc_type` name a PETSc type directly. They are folded into
        ! the option string, which KSPSetFromOptions applies after the numeric
        ! id has set its default, so a name always wins over an id.
        buffer(3) = "ksp_type"
        call get_json_value(json, join(buffer(1:3)), temp_string, &
                    is_required=.false., default_value="")
        call append_petsc_option(self%solver_settings%linear_solver%petsc_options, "-ksp_type", temp_string)

        buffer(3) = "pc_type"
        call get_json_value(json, join(buffer(1:3)), temp_string, &
                    is_required=.false., default_value="")
        call append_petsc_option(self%solver_settings%linear_solver%petsc_options, "-pc_type", temp_string)

        ! Also accept a top-level "linear_solver" block.
        legacy_buffer(1) = linear_solver

        legacy_buffer(2) = max_iterations
        call get_json_value(json, join(legacy_buffer), legacy_int, found=legacy_found, &
                    is_required=.false., valid_range=[1, huge(1)])
        if (legacy_found) self%solver_settings%linear_solver%max_iterations = legacy_int

        legacy_buffer(2) = tolerance
        call get_json_value(json, join(legacy_buffer), legacy_real, found=legacy_found, &
                    is_required=.false., valid_range=[0.0d0, huge(0.0d0)])
        if (legacy_found) self%solver_settings%linear_solver%tolerance = legacy_real

        legacy_buffer(2) = "petsc_options"
        call get_json_value(json, join(legacy_buffer), temp_string, is_required=.false., &
                    default_value=trim(self%solver_settings%linear_solver%petsc_options))
        self%solver_settings%linear_solver%petsc_options = temp_string

        legacy_buffer(2) = "ksp_type"
        call get_json_value(json, join(legacy_buffer), temp_string, is_required=.false., default_value="")
        call append_petsc_option(self%solver_settings%linear_solver%petsc_options, "-ksp_type", temp_string)

        legacy_buffer(2) = "pc_type"
        call get_json_value(json, join(legacy_buffer), temp_string, is_required=.false., default_value="")
        call append_petsc_option(self%solver_settings%linear_solver%petsc_options, "-pc_type", temp_string)

        if (allocated(temp_string)) deallocate (temp_string)

    end subroutine read_solver_settings_linear

    !> Append `name value` to a PETSc option string, ignoring an empty value.
    subroutine append_petsc_option(options, name, value)
        implicit none
        character(len=*), intent(inout) :: options
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: value

        if (len_trim(value) == 0) return
        options = trim(options)//' '//name//' '//trim(adjustl(value))
    end subroutine append_petsc_option

    !======================================================================
    ! Private Helper Functions for Display Routines
    !======================================================================

    pure function nonlinear_solver_to_string(id) result(str)
        integer(int32), intent(in) :: id
        character(len=16) :: str
        str = trim(NONLINEAR_SOLVER%to_name(id))
        if (len_trim(str) == 0) str = "unknown"
    end function nonlinear_solver_to_string

    pure function norm_criteria_to_string(id) result(str)
        integer(int32), intent(in) :: id
        character(len=16) :: str
        str = trim(NONLINEAR_NORM_CRITERIA%to_name(id))
        if (len_trim(str) == 0) str = "unknown"
    end function norm_criteria_to_string

    pure function local_criteria_to_string(id) result(str)
        integer(int32), intent(in) :: id
        character(len=16) :: str
        str = trim(NONLINEAR_CRITERIA%to_name(id))
        if (len_trim(str) == 0) str = "unknown"
    end function local_criteria_to_string

    pure function logic_to_string(id) result(str)
        integer(int32), intent(in) :: id
        character(len=16) :: str
        str = trim(NONLINEAR_LOGIC%to_name(id))
        if (len_trim(str) == 0) str = "unknown"
    end function logic_to_string

    !======================================================================
    ! Display Routines for Solver Settings
    !======================================================================

    module subroutine display_parameters_solver_settings(self)
        !> Display all solver settings. This is the main entry point.
        implicit none
        class(type_solver_settings), intent(in) :: self

        write (*, '(/a)') "======================================================================"
        write (*, '(a)') "                      Solver Settings"
        write (*, '(a)') "======================================================================"

        write (*, '(a, i0)') "  BDF Order             : ", self%bdf_order

        call display_solver_settings_nonlinear(self%nonlinear_solver)
        call display_solver_settings_linear(self%linear_solver)

        write (*, '(a)') "======================================================================"
    end subroutine display_parameters_solver_settings

    subroutine display_solver_settings_nonlinear(nonlinear)
        !> Displays nonlinear solver settings.
        implicit none
        type(type_nonlinear_solver), intent(in) :: nonlinear

        write (*, '(/a)') "  --- Nonlinear Solver ---"
        write (*, '(a, a)') "    Method                : ", trim(nonlinear_solver_to_string(nonlinear%method))

        if (nonlinear%method == NONLINEAR_SOLVER%MODIFIED_NEWTON%ID) then
            write (*, '(a, i0)') "    Update Frequency      : ", nonlinear%update_frequency
        end if

        if (any(nonlinear%method == [NONLINEAR_SOLVER%NEWTON%ID, &
                                     NONLINEAR_SOLVER%MODIFIED_NEWTON%ID, &
                                     NONLINEAR_SOLVER%PICARD%ID])) then
            write (*, '(a, i0)') "    Max Iterations        : ", nonlinear%max_iterations
            write (*, '(a)') "    --- Convergence ---"
            write (*, '(a, a)') "      Use Criteria          : ", trim(norm_criteria_to_string(nonlinear%convergence%use_criteria))
            if (nonlinear%convergence%use_criteria == NONLINEAR_NORM_CRITERIA%BOTH%ID) then
                write (*, '(a, a)') "      Logic Between         : ", trim(logic_to_string(nonlinear%convergence%use_logic))
            end if

            if (any(nonlinear%convergence%use_criteria == [NONLINEAR_NORM_CRITERIA%RESIDUAL%ID, &
                                                           NONLINEAR_NORM_CRITERIA%BOTH%ID])) then
                call display_solver_convergence_criteria(nonlinear%convergence%residual, "      Residual Criteria")
            end if
            if (any(nonlinear%convergence%use_criteria == [NONLINEAR_NORM_CRITERIA%UPDATE%ID, &
                                                           NONLINEAR_NORM_CRITERIA%BOTH%ID])) then
                call display_solver_convergence_criteria(nonlinear%convergence%update, "      Update Criteria")
            end if
        end if
    end subroutine display_solver_settings_nonlinear

    subroutine display_solver_convergence_criteria(criteria_obj, title)
        !> Reusable helper to display convergence criteria details.
        implicit none
        type(type_convergence_criteria), intent(in) :: criteria_obj
        character(*), intent(in) :: title

        write (*, '(a, a)') trim(title)//'        : ', trim(local_criteria_to_string(criteria_obj%criteria))
        if (criteria_obj%criteria == NONLINEAR_CRITERIA%BOTH%ID) then
            write (*, '(a, a)') "        Logic             : ", trim(logic_to_string(criteria_obj%logic))
        end if
        if (any(criteria_obj%criteria == [NONLINEAR_CRITERIA%ABSOLUTE%ID, NONLINEAR_CRITERIA%BOTH%ID])) then
            write (*, '(a, es12.4e2)') "        Absolute Tol.     : ", criteria_obj%absolute_tolerance
        end if
        if (any(criteria_obj%criteria == [NONLINEAR_CRITERIA%RELATIVE%ID, NONLINEAR_CRITERIA%BOTH%ID])) then
            write (*, '(a, es12.4e2)') "        Relative Tol.     : ", criteria_obj%relative_tolerance
        end if
    end subroutine display_solver_convergence_criteria

    subroutine display_solver_settings_linear(linear)
        !> Displays linear solver settings by checking if each method is defined.
        implicit none
        type(type_linear_solver), intent(in) :: linear

        write (*, '(/a)') "  --- Linear Solver ---"

    end subroutine display_solver_settings_linear

end submodule input_basic_solver_settings
