submodule(io_input_conditions) input_conditions_initial
    use :: io_input
    implicit none
    !------------------------------------------------------------------------------------------------------------------------------
    ! JSON key names for initial conditions
    !!------------------------------------------------------------------------------------------------------------------------------
    character(*), parameter :: initial_conditions = "initial_conditions"
    character(*), parameter :: type = "type"
    character(*), parameter :: value = "value"
    character(*), parameter :: valid_initial_condition_types(3) = [character(8) :: "uniform", "laplace", "file"]
    character(*), parameter :: field_name = "field_name"
    !!------------------------------------------------------------------------------------------------------------------------------
contains
    module subroutine read_conditions_initial_conditions(self, json)
        implicit none
        class(type_conditions), intent(inout) :: self
        type(json_file), intent(inout) :: json
        integer(int32) :: i
        logical :: required, found
        character(256) :: buffer(3) = [character(256) :: initial_conditions, "", ""]

        select type (p => self%parent)
        type is (type_input)
            ! An equation that is not solved may still state its initial
            ! condition: that is the value its potential is held at.
            do i = 1, IC_TARGETS%NUM_ID

                required = .true.
                select case (i)
                case (IC_TARGETS%THERMAL%ID)
                    required = p%basic%analysis_controls%is_active(i)
                    buffer(2) = thermal
                case (IC_TARGETS%HYDRAULIC%ID)
                    required = p%basic%analysis_controls%is_active(i)
                    buffer(2) = hydraulic
                case (IC_TARGETS%MECHANICAL%ID)
                    if (.not. p%basic%analysis_controls%is_active(i)) cycle
                    buffer(2) = mechanical
                case (IC_TARGETS%POROSITY%ID)
                    buffer(2) = porosity
                end select
                buffer(3) = type

                call get_json_value(json, join(buffer), self%initial_conditions%physics(i)%type, found=found, &
                                    is_required=required, valid_list=valid_initial_condition_types)
                if (.not. found) cycle

                select case (self%initial_conditions%physics(i)%type)
                case (IC_METHODS%UNIFORM%NAME)
                    buffer(3) = value
                    call get_json_value(json, join(buffer), self%initial_conditions%physics(i)%value, &
                                        is_required=.true.)
                case (IC_METHODS%LAPLACE%NAME)
                    ! No additional parameters needed for laplace
                case (IC_METHODS%FROM_FILE%NAME)
                    buffer(3) = field_name
                    call get_json_value(json, join(buffer), self%initial_conditions%physics(i)%field_name, &
                                        is_required=.true.)
                end select
            end do
        end select

    end subroutine read_conditions_initial_conditions

    module subroutine display_initial_conditions(self)
        implicit none
        class(type_initial_conditions), intent(in) :: self

        write (*, '(A)') "  Initial Conditions:"

        select type (p => self%parent%parent)
        type is (type_input)

            if (p%basic%analysis_controls%is_active(PHYSICS_TYPES%THERMAL%ID)) then
                call display_initial_local(self%physics, "Thermal", IC_TARGETS%THERMAL%ID)
            end if
            if (p%basic%analysis_controls%is_active(PHYSICS_TYPES%HYDRAULIC%ID)) then
                call display_initial_local(self%physics, "Hydraulic", IC_TARGETS%HYDRAULIC%ID)
            end if

            call display_initial_local(self%physics, "Porosity", IC_TARGETS%POROSITY%ID)

        end select

    end subroutine display_initial_conditions

    subroutine display_initial_local(fields, title, target_ic_id)
        implicit none
        type(type_initial_local), intent(in) :: fields(:)
        character(*), intent(in) :: title
        integer(int32), intent(in), optional :: target_ic_id

        write (*, '(A, A)') "    ", trim(title), ":"
        write (*, '(A, A)') "      Type: ", trim(fields(target_ic_id)%type)

        select case (fields(target_ic_id)%type)
        case (IC_METHODS%UNIFORM%NAME)
            write (*, '(A, F8.3)') "      Value: ", fields(target_ic_id)%value
        case (IC_METHODS%LAPLACE%NAME)
            write (*, '(A)') "      Value: Laplace equation will be solved"
        case (IC_METHODS%FROM_FILE%NAME)
            write (*, '(A, A)') "      Field Name: ", trim(fields(target_ic_id)%field_name)
        case default
            write (*, '(A)') "      Unknown type"
        end select
    end subroutine display_initial_local

end submodule input_conditions_initial
