module io_output_manager
    use, intrinsic :: iso_fortran_env
    use :: omp_lib
    use :: stdlib_strings, only:to_string, strip
    use :: stdlib_io, only:open
    use :: module_core
    use :: module_domain
    use :: domain_mesh_plex, only:type_mesh_plex

    use :: io_output_base, only: &
        setup_directory
    use :: io_output_overall, only: &
        abst_output_overall, &
        type_output_overall_vtk, &
        type_output_overall_vtu, &
        type_output_overall_hdf5
    use :: io_output_observation, only:type_output_observation
    use :: io_output_logging, only:type_output_log

    implicit none
    private

    public :: type_output_manager

    !---------------------------------------------------------------------------
    ! type_output
    !---------------------------------------------------------------------------
    type :: type_output_manager
        logical, private, allocatable :: active(:)

        type(type_output_observation), private :: observation
        class(abst_output_overall), private, allocatable :: overall
        type(type_output_log), private :: log

    contains
        procedure, pass(self), public :: initialize => initialize_type_output_manager
        procedure, pass(self), public :: output_fields => output_fields_output_manager
        procedure, pass(self), public :: output_history => output_history_output_manager
        procedure, pass(self), public :: output_system_log => output_system_log_output_manager
        procedure, pass(self), public :: get_log_io_unit => get_log_io_unit_output_manager
        procedure, pass(self), public :: get_obs_fe_ids => get_obs_fe_ids_output_manager
    end type type_output_manager

contains

    subroutine initialize_type_output_manager(self, config_output, config_observation, config_overall, mesh)
        implicit none
        class(type_output_manager), intent(inout) :: self
        type(type_config_output), intent(in) :: config_output
        type(type_config_observation), intent(in) :: config_observation
        type(type_config_overall), intent(in) :: config_overall
        !> The mesh, for a writer that describes the geometry from the DM.
        type(type_mesh_plex), intent(inout), optional, target :: mesh

        character(:), allocatable :: project_path_env
        character(8) :: output_extentions(3) = [".dat", ".csv", ".log"]
        character(8) :: output_file_extentions(7) = &
            [character(8) :: ".dat", ".csv", ".vtk", ".vtu", ".log", ".h5", ".xmf"]
        character(*), parameter :: PROJECT_ENV = "FTCMS_PROJECT_PATH"

        character(:), allocatable :: dir_output
        character(:), allocatable :: dir_output_field

        call get_env_string(PROJECT_ENV, project_path_env)
        call modify_path_format(project_path_env)
        dir_output = strip(project_path_env)//"Output/"
        call setup_directory(dir_output, output_extentions)
        dir_output_field = strip(project_path_env)//"Output/Files/"
        call setup_directory(dir_output_field, output_file_extentions)

        call allocate_array(self%active, source=config_output%is_output_enabled)

        call self%observation%initialize(dir_output, config_observation)

        select case (config_overall%file_format%ID)
        case (FILE_FORMATS%VTK%ID)
            allocate (type_output_overall_vtk :: self%overall)
        case (FILE_FORMATS%VTU%ID)
            allocate (type_output_overall_vtu :: self%overall)
        case (FILE_FORMATS%HDF5%ID, FILE_FORMATS%PETSC_VTU%ID)
            ! Both are written by PETSc; the writer picks the file from the format.
            allocate (type_output_overall_hdf5 :: self%overall)
        end select
        if (allocated(self%overall)) then
            call self%overall%initialize(dir_output_field, config_overall, mesh)
        end if

        call self%log%initialize(dir_output)
    end subroutine initialize_type_output_manager

    subroutine output_fields_output_manager(self, file_counts, temperature, water_content, ice_content, &
                                            vapor_content, pressure, water_flux, time)
        implicit none
        class(type_output_manager), intent(inout) :: self
        integer(int32), intent(in) :: file_counts
        real(real64), intent(in), optional :: temperature(:)
        real(real64), intent(in), optional :: water_content(:)
        real(real64), intent(in), optional :: ice_content(:)
        real(real64), intent(in), optional :: vapor_content(:)
        real(real64), intent(in), optional :: pressure(:)
        type(type_coordinate_array_dp), intent(in), optional :: water_flux
        real(real64), intent(in), optional :: time

        if (.not. allocated(self%overall)) return

        call self%overall%write_fields(file_counts=file_counts, &
                                       temperature=temperature, &
                                       water_content=water_content, &
                                       ice_content=ice_content, &
                                       vapor_content=vapor_content, &
                                       pressure=pressure, &
                                       water_flux=water_flux, &
                                       time=time)
    end subroutine output_fields_output_manager

    subroutine output_history_output_manager(self, time, temperature, water_content, ice_content, &
                                             vapor_content, pressure, water_flux, vapor_flux)
        implicit none
        class(type_output_manager), intent(inout) :: self
        real(real64), intent(in) :: time
        real(real64), intent(in), optional :: temperature(:)
        real(real64), intent(in), optional :: water_content(:)
        real(real64), intent(in), optional :: ice_content(:)
        real(real64), intent(in), optional :: vapor_content(:)
        real(real64), intent(in), optional :: pressure(:)
        type(type_coordinate_array_dp), intent(in), optional :: water_flux
        type(type_coordinate_array_dp), intent(in), optional :: vapor_flux

        call self%observation%output_history(time=time, &
                                             temperature=temperature, &
                                             water_content=water_content, &
                                             ice_content=ice_content, &
                                             vapor_content=vapor_content, &
                                             pressure=pressure, &
                                             water_flux=water_flux, &
                                             vapor_flux=vapor_flux)
    end subroutine output_history_output_manager

    subroutine get_log_io_unit_output_manager(self, io_unit)
        implicit none
        class(type_output_manager), intent(in) :: self
        integer(int32), intent(inout) :: io_unit

        call self%log%get_io_unit(io_unit)
    end subroutine get_log_io_unit_output_manager

    subroutine output_system_log_output_manager(self)
        implicit none
        class(type_output_manager), intent(inout) :: self

        call self%log%output_system_log()
    end subroutine output_system_log_output_manager

    subroutine get_obs_fe_ids_output_manager(self, fe_ids)
        implicit none
        class(type_output_manager), intent(in) :: self
        integer(int32), allocatable, intent(inout) :: fe_ids(:)

        call self%observation%get_obs_fe_ids(fe_ids)
    end subroutine get_obs_fe_ids_output_manager

end module io_output_manager
