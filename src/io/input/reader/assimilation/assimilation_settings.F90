!> The optional `Input/DataAssimilation.json`: the atmospheric ensemble and its observations.
!>
!> Absent or `"enabled": false` means no assimilation. When enabled, the file
!> names the boundary the ensemble exchanges with the soil (`surface_boundary`,
!> a "Face Sets" value); it is not assumed. The controller reads the
!> observation files from the project's `Input` directory.
module io_input_assimilation
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: json_module, only:json_file
    use :: io_input_base, only:get_json_value
    implicit none
    private

    public :: type_assimilation_settings

    type :: type_assimilation_settings
        logical :: enabled = .false.
        !> The "Face Sets" value of the surface the atmosphere acts on.
        integer(int32) :: surface_boundary = -1
        integer(int32) :: ensemble_size = 30
        real(real64) :: max_height = 200.0d0
        integer(int32) :: num_nodes = 50
        character(len=:), allocatable :: csv_file
        character(len=:), allocatable :: upper_bc_file
        real(real64) :: interval_seconds = 3600.0d0
        real(real64) :: sigma_T = 0.5d0
        real(real64) :: sigma_q = 1.0d-4
        real(real64) :: sigma_U = 0.5d0
        real(real64) :: z0 = 0.01d0
        real(real64) :: albedo = 0.2d0
        real(real64) :: emissivity = 0.95d0
        real(real64) :: Pmin = -1.0d5
        real(real64) :: Pmax = 0.0d0
        real(real64) :: lambda_soil = 1.5d0
        real(real64) :: stomatal_resistance = 50.0d0
        character(len=:), allocatable :: reference_datetime
        real(real64) :: latitude = 35.67d0
        real(real64) :: longitude = 139.48d0
        real(real64) :: tau_atm = 0.6d0
        real(real64) :: utc_offset_hours = 0.0d0
    contains
        procedure, pass(self) :: initialize => initialize_assimilation_settings
    end type type_assimilation_settings

contains

    !> Read `file_name` if it exists; otherwise leave assimilation off.
    !>
    !> Fails by `error stop` when assimilation is enabled without a surface boundary.
    subroutine initialize_assimilation_settings(self, file_name)
        implicit none
        class(type_assimilation_settings), intent(inout) :: self
        character(len=*), intent(in) :: file_name

        type(json_file) :: json
        logical :: exists

        inquire (file=file_name, exist=exists)
        self%enabled = .false.
        if (.not. exists) return

        call json%initialize()
        call json%load(filename=file_name)
        if (json%failed()) error stop "DataAssimilation.json could not be read"

        call get_json_value(json, "enabled", self%enabled, is_required=.false., default_value=.true.)
        if (.not. self%enabled) then
            call json%destroy()
            return
        end if
        call get_json_value(json, "surface_boundary", self%surface_boundary, is_required=.true.)
        call get_json_value(json, "ensemble.size", self%ensemble_size, default_value=self%ensemble_size)
        call get_json_value(json, "ensemble.max_height", self%max_height, default_value=self%max_height)
        call get_json_value(json, "ensemble.num_nodes", self%num_nodes, default_value=self%num_nodes)
        call get_json_value(json, "observation.csv_file", self%csv_file, default_value="obs_data.csv")
        call get_json_value(json, "observation.upper_bc_file", self%upper_bc_file, default_value="upper_bc.csv")
        call get_json_value(json, "observation.interval_seconds", self%interval_seconds, &
                            default_value=self%interval_seconds)
        call get_json_value(json, "observation.sigma_T", self%sigma_T, default_value=self%sigma_T)
        call get_json_value(json, "observation.sigma_q", self%sigma_q, default_value=self%sigma_q)
        call get_json_value(json, "observation.sigma_U", self%sigma_U, default_value=self%sigma_U)
        call get_json_value(json, "surface.z0", self%z0, default_value=self%z0)
        call get_json_value(json, "surface.albedo", self%albedo, default_value=self%albedo)
        call get_json_value(json, "surface.emissivity", self%emissivity, default_value=self%emissivity)
        call get_json_value(json, "surface.Pmin", self%Pmin, default_value=self%Pmin)
        call get_json_value(json, "surface.Pmax", self%Pmax, default_value=self%Pmax)
        call get_json_value(json, "surface.lambda_soil", self%lambda_soil, default_value=self%lambda_soil)
        call get_json_value(json, "surface.stomatal_resistance", self%stomatal_resistance, &
                            default_value=self%stomatal_resistance)
        call get_json_value(json, "reference_datetime", self%reference_datetime, &
                            default_value="2000-01-01T00:00:00")
        call get_json_value(json, "solar.latitude", self%latitude, default_value=self%latitude)
        call get_json_value(json, "solar.longitude", self%longitude, default_value=self%longitude)
        call get_json_value(json, "solar.tau_atm", self%tau_atm, default_value=self%tau_atm)
        call get_json_value(json, "solar.utc_offset_hours", self%utc_offset_hours, &
                            default_value=self%utc_offset_hours)
        call json%destroy()

    end subroutine initialize_assimilation_settings

end module io_input_assimilation
