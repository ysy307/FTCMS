submodule(io_input_basic) input_basic_discretization_settings
    implicit none
    !!------------------------------------------------------------------------------------------------------------------------------
    ! JSON key names for the discretisation of the governing equations
    !!------------------------------------------------------------------------------------------------------------------------------
    character(*), parameter :: discretization_settings = "discretization_settings"
    character(*), parameter :: space = "space"
    character(*), parameter :: degree = "degree"
    character(*), parameter :: coefficient_evaluation = "coefficient_evaluation"
    character(*), parameter :: space_strings(2) = [character(len=16) :: "lagrange", "mixed"]
    character(*), parameter :: evaluation_strings(2) = [character(len=16) :: "quadrature", "nodal"]
    !!------------------------------------------------------------------------------------------------------------------------------
contains
    !> Read `discretization_settings.<physics>.{space, degree, coefficient_evaluation}`.
    !>
    !> One block per active governing equation, keyed by its `PHYSICS_TYPES`
    !> name, with the same three keys for every equation. A missing block or key
    !> keeps the default of `type_discretization_settings`: a degree-one
    !> Lagrange space with the laws evaluated at the quadrature points. Only the
    !> names and the range are checked here; whether the engine implements a
    !> combination is the engine's decision.
    module subroutine read_discretization_settings(self, json)
        implicit none
        class(type_input_basic), intent(inout) :: self
        type(json_file), intent(inout) :: json

        character(256) :: buffer(3) = [character(256) :: discretization_settings, "", ""]
        character(:), allocatable :: name
        integer(int32) :: id
        logical :: found

        associate (settings => self%discretization_settings)
            settings%is_given = .false.
            do id = 1, PHYSICS_TYPES%NUM_ID
                if (.not. self%analysis_controls%is_active(id)) cycle
                buffer(2) = PHYSICS_TYPES%to_name(id)

                buffer(3) = space
                call get_json_value(json, join(buffer), name, found=found, is_required=.false., &
                                    default_value="lagrange", valid_list=space_strings)
                settings%is_given = settings%is_given .or. found
                settings%space(id) = FE_SPACES%to_id(name)

                buffer(3) = degree
                call get_json_value(json, join(buffer), settings%degree(id), found=found, is_required=.false., &
                                    default_value=1_int32, valid_range=[0_int32, 8_int32])
                settings%is_given = settings%is_given .or. found

                buffer(3) = coefficient_evaluation
                call get_json_value(json, join(buffer), name, found=found, is_required=.false., &
                                    default_value="quadrature", valid_list=evaluation_strings)
                settings%is_given = settings%is_given .or. found
                settings%coefficient_evaluation(id) = COEFFICIENT_EVALUATIONS%to_id(name)
            end do
        end associate
    end subroutine read_discretization_settings

    module subroutine display_discretization_settings(self, is_active)
        implicit none
        class(type_discretization_settings), intent(in) :: self
        logical, intent(in) :: is_active(:)

        integer(int32) :: id

        do id = 1, PHYSICS_TYPES%NUM_ID
            if (.not. is_active(id)) cycle
            write (*, '(a)') "Discretization "//PHYSICS_TYPES%to_name(id)//": " &
                //FE_SPACES%to_name(self%space(id))//", degree "//to_string(self%degree(id)) &
                //", "//COEFFICIENT_EVALUATIONS%to_name(self%coefficient_evaluation(id))
        end do
    end subroutine display_discretization_settings
end submodule input_basic_discretization_settings
