module core_types_config_base
    use, intrinsic :: iso_fortran_env
    use :: core_constants, only:type_constant_id, type_constant_value, ERROR_CODES
    use :: core_validation, only:raise_error
    use :: core_types_geometry_coordinate, only:type_coordinate_dp
    use :: core_types_geometry_coordinate_array, only:type_coordinate_array_dp
    use :: core_types_utils_datetime, only:type_datetime

    implicit none
    private

    public :: abst_config

    type, abstract :: abst_config
    contains
        procedure, private, pass(self) :: set_int8
        procedure, private, pass(self) :: set_int8_1d
        procedure, private, pass(self) :: set_int8_2d
        procedure, private, pass(self) :: set_int8_3d
        procedure, private, pass(self) :: set_int32
        procedure, private, pass(self) :: set_int32_1d
        procedure, private, pass(self) :: set_int32_2d
        procedure, private, pass(self) :: set_int32_3d
        procedure, private, pass(self) :: set_real64
        procedure, private, pass(self) :: set_real64_1d
        procedure, private, pass(self) :: set_real64_2d
        procedure, private, pass(self) :: set_real64_3d
        procedure, private, pass(self) :: set_logical
        procedure, private, pass(self) :: set_logical_1d
        procedure, private, pass(self) :: set_logical_2d
        procedure, private, pass(self) :: set_logical_3d
        procedure, private, pass(self) :: set_character
        procedure, private, pass(self) :: set_character_1d
        procedure, private, pass(self) :: set_character_2d
        procedure, private, pass(self) :: set_character_3d

        procedure, private, pass(self) :: set_constant_id
        procedure, private, pass(self) :: set_constant_ids
        procedure, private, pass(self) :: set_constant_value
        procedure, private, pass(self) :: set_constant_values
        procedure, private, pass(self) :: set_coordinate_array_dp
        procedure, private, pass(self) :: set_coordinate_dp
        procedure, private, pass(self) :: set_datetime
        procedure, private, pass(self) :: set_datetimes
        generic :: set => &
            set_int8, set_int8_1d, set_int8_2d, set_int8_3d, &
            set_int32, set_int32_1d, set_int32_2d, set_int32_3d, &
            set_real64, set_real64_1d, set_real64_2d, set_real64_3d, &
            set_logical, set_logical_1d, set_logical_2d, set_logical_3d, &
            set_character, set_character_1d, set_character_2d, set_character_3d, &
            set_constant_id, set_constant_ids, set_constant_value, set_constant_values, &
            set_coordinate_array_dp, set_coordinate_dp, set_datetime, set_datetimes

        procedure(abst_copy_config), public, pass(self), deferred :: copy
        procedure(abst_reset_config), public, pass(self), deferred :: reset
    end type abst_config

    abstract interface
        subroutine abst_copy_config(self, source)
            import :: abst_config
            implicit none
            class(abst_config), intent(inout) :: self
            class(abst_config), intent(in) :: source
        end subroutine abst_copy_config

        subroutine abst_reset_config(self)
            import :: abst_config
            implicit none
            class(abst_config), intent(inout) :: self
        end subroutine abst_reset_config
    end interface

contains

    subroutine set_int8(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        integer(int8), intent(inout) :: member
        integer(int8), intent(in) :: value

        member = value
    end subroutine set_int8

    subroutine set_int8_1d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        integer(int8), allocatable, intent(inout) :: member(:)
        integer(int8), allocatable, intent(in) :: value(:)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_int8_1d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_int8_1d

    subroutine set_int8_2d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        integer(int8), allocatable, intent(inout) :: member(:, :)
        integer(int8), allocatable, intent(in) :: value(:, :)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_int8_2d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_int8_2d

    subroutine set_int8_3d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        integer(int8), allocatable, intent(inout) :: member(:, :, :)
        integer(int8), allocatable, intent(in) :: value(:, :, :)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_int8_3d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_int8_3d

    subroutine set_int32(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        integer(int32), intent(inout) :: member
        integer(int32), intent(in) :: value

        member = value
    end subroutine set_int32

    subroutine set_int32_1d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        integer(int32), allocatable, intent(inout) :: member(:)
        integer(int32), allocatable, intent(in) :: value(:)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_int32_1d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_int32_1d

    subroutine set_int32_2d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        integer(int32), allocatable, intent(inout) :: member(:, :)
        integer(int32), allocatable, intent(in) :: value(:, :)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_int32_2d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_int32_2d

    subroutine set_int32_3d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        integer(int32), allocatable, intent(inout) :: member(:, :, :)
        integer(int32), allocatable, intent(in) :: value(:, :, :)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_int32_3d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_int32_3d

    ! --- real(real64) implementations ---

    subroutine set_real64(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        real(real64), intent(inout) :: member
        real(real64), intent(in) :: value

        member = value
    end subroutine set_real64

    subroutine set_real64_1d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        real(real64), allocatable, intent(inout) :: member(:)
        real(real64), allocatable, intent(in) :: value(:)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_real64_1d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_real64_1d

    subroutine set_real64_2d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        real(real64), allocatable, intent(inout) :: member(:, :)
        real(real64), allocatable, intent(in) :: value(:, :)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_real64_2d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_real64_2d

    subroutine set_real64_3d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        real(real64), allocatable, intent(inout) :: member(:, :, :)
        real(real64), allocatable, intent(in) :: value(:, :, :)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_real64_3d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_real64_3d

    subroutine set_logical(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        logical, intent(inout) :: member
        logical, intent(in) :: value

        member = value
    end subroutine set_logical

    subroutine set_logical_1d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        logical, allocatable, intent(inout) :: member(:)
        logical, allocatable, intent(in) :: value(:)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_logical_1d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_logical_1d

    subroutine set_logical_2d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        logical, allocatable, intent(inout) :: member(:, :)
        logical, allocatable, intent(in) :: value(:, :)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_logical_2d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_logical_2d

    subroutine set_logical_3d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        logical, allocatable, intent(inout) :: member(:, :, :)
        logical, allocatable, intent(in) :: value(:, :, :)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_logical_3d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_logical_3d

    subroutine set_character(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        character(:), allocatable, intent(inout) :: member
        character(:), allocatable, intent(in) :: value

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_character")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_character

    subroutine set_character_1d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        character(:), allocatable, intent(inout) :: member(:)
        character(:), allocatable, intent(in) :: value(:)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_character_1d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_character_1d

    subroutine set_character_2d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        character(:), allocatable, intent(inout) :: member(:, :)
        character(:), allocatable, intent(in) :: value(:, :)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_character_2d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_character_2d

    subroutine set_character_3d(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        character(:), allocatable, intent(inout) :: member(:, :, :)
        character(:), allocatable, intent(in) :: value(:, :, :)

        if (.not. allocated(value)) then
            call raise_error(ERROR_CODES%NOT_ALLOCATED, opt="value", scope="core_types_config_base:set_character_3d")
        end if

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_character_3d

    subroutine set_constant_id(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        type(type_constant_id), intent(inout) :: member
        type(type_constant_id), intent(in) :: value

        member = value
    end subroutine set_constant_id

    subroutine set_constant_ids(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        type(type_constant_id), intent(inout), allocatable :: member(:)
        type(type_constant_id), intent(in) :: value(:)

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_constant_ids

    subroutine set_constant_value(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        type(type_constant_value), intent(inout) :: member
        type(type_constant_value), intent(in) :: value

        member = value
    end subroutine set_constant_value

    subroutine set_constant_values(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        type(type_constant_value), intent(inout), allocatable :: member(:)
        type(type_constant_value), intent(in) :: value(:)

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_constant_values

    subroutine set_coordinate_array_dp(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        type(type_coordinate_array_dp), intent(inout) :: member
        type(type_coordinate_array_dp), intent(in) :: value

        member = value
    end subroutine set_coordinate_array_dp

    subroutine set_coordinate_dp(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        type(type_coordinate_dp), intent(inout) :: member
        type(type_coordinate_dp), intent(in) :: value

        member = value
    end subroutine set_coordinate_dp

    subroutine set_datetime(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        type(type_datetime), intent(inout) :: member
        type(type_datetime), intent(in) :: value

        member = value
    end subroutine set_datetime

    subroutine set_datetimes(self, member, value)
        implicit none
        class(abst_config), intent(in) :: self
        type(type_datetime), intent(inout),allocatable :: member(:)
        type(type_datetime), intent(in) :: value(:)

        if (allocated(member)) deallocate (member)
        allocate (member, source=value)
    end subroutine set_datetimes

end module core_types_config_base
