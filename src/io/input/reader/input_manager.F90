module io_input
    use, intrinsic :: iso_fortran_env, only: int32, real64, output_unit
    use :: core_parallel_mpi
    use :: stdlib_strings, only:to_string, strip, ends_with
    use :: stdlib_logger
    use :: json_module, only:json_file
    use :: module_core
    use :: io_input_base, only:abst_input
    use :: io_input_basic, only:type_input_basic
    use :: io_input_conditions, only:type_conditions
    use :: io_input_output_conditions, only:type_output_settings
    use :: io_input_geometry, only:type_input_geometry
    use :: io_input_assimilation, only:type_assimilation_settings
    ! use :: module_physics
    implicit none
    private

    public :: type_input

    type, extends(abst_input) :: type_input
        character(:), allocatable :: input_path

        type(type_input_basic) :: basic
        type(type_conditions) :: conditions
        type(type_output_settings) :: output_settings
        type(type_input_geometry) :: geometry
        !> The optional DataAssimilation.json.
        type(type_assimilation_settings) :: assimilation
    contains
        procedure, pass(self), public :: initialize => initialize_type_input
        ! procedure, pass(self), public :: get_density_info => get_density_info_input
        ! procedure, pass(self), public :: get_specific_heat_info => get_specific_heat_info_input
        ! procedure, pass(self), public :: get_volumetric_heat_capacity_info => get_volumetric_heat_capacity_info_input
        ! procedure, pass(self), public :: get_thermal_conductivity_info => get_thermal_conductivity_info_input
        ! procedure, pass(self), public :: get_wrf_info => get_wrf_info_input
        ! procedure, pass(self), public :: get_gcc_info => get_gcc_info_input
        ! procedure, pass(self), public :: get_hcf_info => get_hcf_info_input
        procedure, pass(self), public :: display => display_input
    end type type_input

contains

    subroutine initialize_type_input(self)
        implicit none
        class(type_input), intent(inout), target :: self

        character(:), allocatable :: project_path_env

        character(*), parameter :: PROJECT_ENV = "FTCMS_PROJECT_PATH"

        integer(int32) :: status
        integer(int32) :: ierr, myrank
        integer(int32) :: error_flag = 0
        integer(int32) :: input_path_length

        call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

        call get_env_string(PROJECT_ENV, project_path_env, status)
        if (status /= 0) then
            if (myrank == 0) then
                print *, "FATAL ERROR: Environment variable '"//PROJECT_ENV//"' is not set or is empty. Aborting."
            end if
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call modify_path_format(project_path_env)

        if (myrank == 0) then
            self%input_path = strip(project_path_env)//"Input/"
            input_path_length = len_trim(self%input_path)
            if (.not. file_exists(self%input_path//"Basic.json")) error_flag = 2
            if (.not. file_exists(self%input_path//"Conditions.json")) error_flag = 3
            if (.not. file_exists(self%input_path//"Output.json")) error_flag = 4
        end if

        ! 2. Broadcast rank 0 check result (error_flag) to all ranks
        call MPI_Bcast(error_flag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        ! 3. Abort all ranks if an error was detected
        if (error_flag /= 0) then
            if (myrank == 0) then
                print *, "FATAL ERROR: A required input file was not found. Aborting."
                select case (error_flag)
                case (2); print *, "-> Basic.json is missing."
                case (3); print *, "-> Conditions.json is missing."
                case (4); print *, "-> Output.json is missing."
                end select
            end if
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call MPI_Bcast(input_path_length, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (myrank /= 0) then
            allocate (character(len=input_path_length) :: self%input_path)
        end if
        call MPI_Bcast(self%input_path, input_path_length, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        self%basic%file_name = self%input_path//"Basic.json"
        self%conditions%file_name = self%input_path//"Conditions.json"
        self%output_settings%file_name = self%input_path//"Output.json"

        self%basic%parent => self
        self%conditions%parent => self
        self%output_settings%parent => self

        call self%basic%initialize()
        call self%conditions%initialize()
        self%geometry%parent => self
        call self%geometry%initialize()
        call self%output_settings%initialize()
        call self%assimilation%initialize(self%input_path//"DataAssimilation.json")

    end subroutine initialize_type_input

    subroutine display_input(self)
        implicit none
        class(type_input), intent(in) :: self

        integer(int32) :: ierr, myrank

        call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
        if (myrank == 0) then
            write (*, '(A)') "=== Simulation Settings ==="
            call self%basic%display()
            write (*, '(A)') "=== Conditions ==="
            call self%conditions%display()
            write (*, '(A)') "=== Output Settings ==="
            call self%output_settings%display()
        end if

    end subroutine display_input

    ! subroutine get_density_info_input(self, material_id, density_info)
    !     implicit none
    !     class(type_input), intent(in) :: self
    !     integer(int32), intent(in) :: material_id
    !     type(type_physics_info), intent(inout) :: density_info

    !     if (material_id < 1 .or. material_id > self%basic%num_materials) then
    !         error stop "Input Error: material_id is out of range in get_density_info."
    !     end if

    !     call self%basic%materials(material_id)%get_density_info(density_info)

    ! end subroutine get_density_info_input

    ! subroutine get_specific_heat_info_input(self, material_id, specific_heat_info)
    !     implicit none
    !     class(type_input), intent(in) :: self
    !     integer(int32), intent(in) :: material_id
    !     type(type_physics_info), intent(inout) :: specific_heat_info

    !     if (material_id < 1 .or. material_id > self%basic%num_materials) then
    !         error stop "Input Error: material_id is out of range in get_specific_heat_info."
    !     end if

    !     call self%basic%materials(material_id)%get_specific_heat_info(specific_heat_info)

    ! end subroutine get_specific_heat_info_input

    ! subroutine get_volumetric_heat_capacity_info_input(self, material_id, volumetric_heat_capacity_info)
    !     implicit none
    !     class(type_input), intent(in) :: self
    !     integer(int32), intent(in) :: material_id
    !     type(type_physics_info), intent(inout) :: volumetric_heat_capacity_info

    !     if (material_id < 1 .or. material_id > self%basic%num_materials) then
    !         error stop "Input Error: material_id is out of range in get_volumetric_heat_capacity_info."
    !     end if

    !     call self%basic%materials(material_id)%get_volumetric_heat_capacity_info(volumetric_heat_capacity_info)

    ! end subroutine get_volumetric_heat_capacity_info_input

    ! subroutine get_thermal_conductivity_info_input(self, material_id, thermal_conductivity_info)
    !     implicit none
    !     class(type_input), intent(in) :: self
    !     integer(int32), intent(in) :: material_id
    !     type(type_physics_info), intent(inout) :: thermal_conductivity_info

    !     if (material_id < 1 .or. material_id > self%basic%num_materials) then
    !         error stop "Input Error: material_id is out of range in get_thermal_conductivity_info."
    !     end if

    !     call self%basic%materials(material_id)%get_thermal_conductivity_info(thermal_conductivity_info)

    ! end subroutine get_thermal_conductivity_info_input

    ! subroutine get_wrf_info_input(self, material_id, wrf_info)
    !     implicit none
    !     class(type_input), intent(in) :: self
    !     integer(int32), intent(in) :: material_id
    !     type(type_wrf_params), intent(inout) :: wrf_info

    !     if (material_id < 1 .or. material_id > self%basic%num_materials) then
    !         error stop "Input Error: material_id is out of range in get_wrf_info."
    !     end if

    !     call self%basic%materials(material_id)%get_wrf_info(wrf_info)

    ! end subroutine get_wrf_info_input

    ! subroutine get_gcc_info_input(self, material_id, gcc_info)
    !     implicit none
    !     class(type_input), intent(in) :: self
    !     integer(int32), intent(in) :: material_id
    !     integer(int32), intent(inout) :: gcc_info

    !     if (material_id < 1 .or. material_id > self%basic%num_materials) then
    !         error stop "Input Error: material_id is out of range in get_gcc_info."
    !     end if

    !     call self%basic%materials(material_id)%get_gcc_info(gcc_info)

    ! end subroutine get_gcc_info_input

    ! subroutine get_hcf_info_input(self, material_id, hcf_info)
    !     implicit none
    !     class(type_input), intent(in) :: self
    !     integer(int32), intent(in) :: material_id
    !     type(type_hcf_params), intent(inout) :: hcf_info

    !     if (material_id < 1 .or. material_id > self%basic%num_materials) then
    !         error stop "Input Error: material_id is out of range in get_hcf_info."
    !     end if

    !     call self%basic%materials(material_id)%get_hcf_info(hcf_info)

    ! end subroutine get_hcf_info_input

    !>
    !> Checks if a file exists by attempting to open it.
    !> This method is more reliable than using INQUIRE(FILE=...) in some environments.
    !>
    function file_exists(file_path) result(found)
        implicit none
        !> The path of the file to check.
        character(*), intent(in) :: file_path
        !> `.true.` if the file exists and is readable, `.false.` otherwise.
        logical :: found
        integer(int32) :: unit_num, io_status

        open (newunit=unit_num, file=trim(file_path), status="old", action="read", iostat=io_status)

        if (io_status == 0) then
            ! If opened successfully, the file exists. Close it immediately.
            close (unit_num)
            found = .true.
        else
            ! If opened failed, the file does not exist or is not accessible.
            found = .false.
        end if

    end function file_exists

end module io_input
