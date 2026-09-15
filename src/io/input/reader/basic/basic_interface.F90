module io_input_basic
    use, intrinsic :: iso_fortran_env
!$  use :: omp_lib
    use :: core_parallel_mpi
    use :: stdlib_strings, only:to_string, strip, ends_with
    use :: stdlib_logger
    use :: json_module, only:json_file
    use :: module_core
    use :: io_input_base, only:get_json_value, abst_input
    ! use :: module_physics
    implicit none
    private

    public :: type_input_basic

    character(*), parameter :: thermal = "thermal"
    character(*), parameter :: hydraulic = "hydraulic"
    character(*), parameter :: mechanical = "mechanical"
    character(*), parameter :: calculate_thermal = "calculate_thermal"
    character(*), parameter :: calculate_hydraulic = "calculate_hydraulic"
    character(*), parameter :: calculate_mechanical = "calculate_mechanical"
    character(*), parameter :: enable_vapor_transport = "enable_vapor_transport"

    !!------------------------------------------------------------------------------------------------------------------------------
    type :: type_simulation_settings
        character(:), allocatable :: title
        integer(int32) :: calculate_type
        integer(int32) :: calculate_dimension
    contains
        procedure, pass(self) :: display => display_simulation_settings
    end type type_simulation_settings

    interface
        module subroutine display_simulation_settings(self)
            implicit none
            class(type_simulation_settings) :: self
        end subroutine display_simulation_settings
    end interface
    !!------------------------------------------------------------------------------------------------------------------------------
    type :: type_analysis_controls
        logical :: is_active(IC_TARGETS%NUM_ID)
        logical :: enable_vapor_transport = .true.
        character(:), allocatable :: coupling_mode
    contains
        procedure, pass(self) :: display => display_analysis_controls
    end type type_analysis_controls

    interface
        module subroutine display_analysis_controls(self)
            implicit none
            class(type_analysis_controls) :: self
        end subroutine display_analysis_controls
    end interface
    !!------------------------------------------------------------------------------------------------------------------------------
    type :: type_geometry_settings
        character(:), allocatable :: file_name
        integer(int32) :: integration_order
        character(:), allocatable :: integration_type
        real(real64), allocatable :: integration_points(:)
    contains
        procedure, pass(self) :: display => display_geometry_settings
    end type type_geometry_settings

    interface
        module subroutine display_geometry_settings(self)
            implicit none
            class(type_geometry_settings) :: self

        end subroutine display_geometry_settings
    end interface
    !!------------------------------------------------------------------------------------------------------------------------------

    type :: type_materials_water_characteristic_curve
        integer(int32) :: model_number = -1
        integer(int32) :: unit = -1
        real(real64) :: theta_s = 0.0d0
        real(real64) :: theta_r = 0.0d0
        real(real64) :: alpha1 = 0.0d0
        real(real64) :: n1 = 0.0d0
        real(real64) :: m1 = 0.0d0
        real(real64) :: w1 = 0.0d0
        real(real64) :: h_crit = 0.0d0
        real(real64) :: alpha2 = 0.0d0
        real(real64) :: n2 = 0.0d0
        real(real64) :: m2 = 0.0d0
        real(real64) :: w2 = 0.0d0
        real(real64) :: l = 0.0d0
    end type type_materials_water_characteristic_curve

    type :: type_equilibrium_model
        logical :: segregation
    end type type_equilibrium_model

    type :: type_phase_change
        type(type_equilibrium_model) :: equilibrium_model
        real(real64) :: segregation_potential = 0.0d0
        real(real64) :: T_fringe_low = -1.0d0
        real(real64) :: T_fringe_high = 0.0d0
    end type type_phase_change

    type :: type_hydraulic_conductivity_model
        integer(int32) :: model_number
        real(real64) :: saturated_conductivity
        real(real64) :: impedance_factor
        integer(int32) :: water_viscosity_model
        real(real64) :: gain_factor
    end type type_hydraulic_conductivity_model

    type :: type_hydraulic_storage
        real(real64) :: specific_storage = 0.0d0
    end type type_hydraulic_storage

    type :: type_materials_gcc
        logical :: is_segregation
        character(:), allocatable :: unit
    end type type_materials_gcc

    type :: type_materials_density
        real(real64), allocatable :: value(:)
    end type type_materials_density

    type :: type_materials_specific_heat
        real(real64), allocatable :: value(:)
    end type type_materials_specific_heat

    type :: type_materials_volumetric_heat_capacity
        real(real64), allocatable :: value(:)
        real(real64), allocatable :: params(:)
    end type type_materials_volumetric_heat_capacity

    type :: type_materials_thermal_conductivity
        real(real64), allocatable :: value(:)
        logical :: is_dispersed
        real(real64), allocatable :: dispersivity(:)
        real(real64), allocatable :: params(:)
    end type type_materials_thermal_conductivity

    type :: type_material_settings
        integer(int32) :: id
        character(:), allocatable :: name
        logical :: is_active(PHYSICS_TYPES%NUM_ID)
        integer(int32) :: phase
        type(type_materials_density) :: density
        type(type_materials_specific_heat) :: specific_heat
        type(type_materials_volumetric_heat_capacity) :: volumetric_heat_capacity
        type(type_materials_thermal_conductivity) :: thermal_conductivity

        type(type_materials_water_characteristic_curve) :: water_characteristic_curve
        type(type_phase_change) :: phase_change
        type(type_hydraulic_conductivity_model) :: hydraulic_conductivity_model
        type(type_hydraulic_storage) :: hydraulic_storage
    contains
        procedure, pass(self) :: display => display_material_settings
        ! procedure, pass(self) :: get_density_info => get_density_info_material_settings
        ! procedure, pass(self) :: get_specific_heat_info => get_specific_heat_info_material_settings
        ! procedure, pass(self) :: get_volumetric_heat_capacity_info => get_volumetric_heat_capacity_info_material_settings
        ! procedure, pass(self) :: get_thermal_conductivity_info => get_thermal_conductivity_info_material_settings
        ! procedure, pass(self) :: get_wrf_info => get_wrf_info_material_settings
        ! procedure, pass(self) :: get_gcc_info => get_gcc_info_material_settings
        ! procedure, pass(self) :: get_hcf_info => get_hcf_info_material_settings
    end type type_material_settings

    interface
        module subroutine display_material_settings(self)
            implicit none
            class(type_material_settings), intent(in) :: self
        end subroutine display_material_settings

        ! module subroutine get_density_info_material_settings(self, density_info)
        !     implicit none
        !     class(type_material_settings), intent(in) :: self
        !     type(type_physics_info), intent(inout) :: density_info
        ! end subroutine get_density_info_material_settings

        ! module subroutine get_specific_heat_info_material_settings(self, specific_heat_info)
        !     implicit none
        !     class(type_material_settings), intent(in) :: self
        !     type(type_physics_info), intent(inout) :: specific_heat_info

        ! end subroutine get_specific_heat_info_material_settings

        ! module subroutine get_volumetric_heat_capacity_info_material_settings(self, volumetric_heat_capacity_info)
        !     implicit none
        !     class(type_material_settings), intent(in) :: self
        !     type(type_physics_info), intent(inout) :: volumetric_heat_capacity_info

        ! end subroutine get_volumetric_heat_capacity_info_material_settings

        ! module subroutine get_thermal_conductivity_info_material_settings(self, thermal_conductivity_info)
        !     implicit none
        !     class(type_material_settings), intent(in) :: self
        !     type(type_physics_info), intent(inout) :: thermal_conductivity_info

        ! end subroutine get_thermal_conductivity_info_material_settings

        ! module subroutine get_wrf_info_material_settings(self, wrf_info)
        !     implicit none
        !     class(type_material_settings), intent(in) :: self
        !     type(type_wrf_params), intent(inout) :: wrf_info

        ! end subroutine get_wrf_info_material_settings

        ! module subroutine get_gcc_info_material_settings(self, gcc_info)
        !     implicit none
        !     class(type_material_settings), intent(in) :: self
        !     integer(int32), intent(inout) :: gcc_info

        ! end subroutine get_gcc_info_material_settings

        ! module subroutine get_hcf_info_material_settings(self, hcf_info)
        !     implicit none
        !     class(type_material_settings), intent(in) :: self
        !     type(type_hcf_params), intent(inout) :: hcf_info

        ! end subroutine get_hcf_info_material_settings
    end interface
    !!------------------------------------------------------------------------------------------------------------------------------
    type :: type_convergence_criteria
        integer(int32) :: criteria
        integer(int32) :: logic
        real(real64) :: absolute_tolerance
        real(real64) :: relative_tolerance
    end type type_convergence_criteria

    type :: type_convergence
        integer(int32) :: use_criteria
        integer(int32) :: norm_type
        integer(int32) :: use_logic
        type(type_convergence_criteria) :: residual
        type(type_convergence_criteria) :: update
        !> Conserved-quantity convergence tolerances (used when use_criteria == CONSERVED)
        real(real64) :: atol_enthalpy = 1.0d2
        real(real64) :: atol_density = 1.0d-3
        real(real64) :: rtol_conserved = 1.0d-3
        real(real64) :: residual_eps = 1.0d-1
    end type type_convergence

    type :: type_nonlinear_solver
        integer(int32) :: method
        integer(int32) :: update_frequency
        integer(int32) :: max_iterations
        type(type_convergence) :: convergence
    end type type_nonlinear_solver

    type :: type_linear_solver
        integer(int32) :: max_iterations = 10000
        real(real64) :: tolerance = 1.0d-6
        !> Raw PETSc option string, inserted into the option database before
        !> KSPSetFromOptions. Anything -ksp_*/-pc_* accepts can be set here.
        character(len=512) :: petsc_options = ''
    end type type_linear_solver

    type :: type_parallel_threads
        logical :: is_parallel
        integer(int32) :: num_threads
        character(:), allocatable :: schedule
        integer(int32) :: max_active_levels
    end type type_parallel_threads

    type :: type_parallel_settings
        type(type_parallel_threads) :: threads
    end type type_parallel_settings

    type :: type_solver_settings
        integer(int32) :: bdf_order
        character(:), allocatable :: reordering
        character(:), allocatable :: coloring
        type(type_nonlinear_solver) :: nonlinear_solver
        type(type_linear_solver) :: linear_solver
        type(type_parallel_settings) :: parallel_settings
    contains
        procedure, pass(self) :: display => display_parameters_solver_settings
    end type

    interface
        module subroutine display_parameters_solver_settings(self)
            implicit none
            class(type_solver_settings), intent(in) :: self
        end subroutine display_parameters_solver_settings
    end interface
    !!------------------------------------------------------------------------------------------------------------------------------
    type :: type_input_basic
        class(abst_input), pointer :: parent => null()
        character(:), allocatable :: file_name
        type(type_simulation_settings) :: simulation_settings
        type(type_analysis_controls) :: analysis_controls
        type(type_geometry_settings) :: geometry_settings
        integer(int32) :: num_materials
        type(type_material_settings), allocatable :: materials(:)
        type(type_solver_settings) :: solver_settings
    contains
        procedure, pass(self) :: initialize => initialize_type_input_basic
        procedure, pass(self) :: display => display_input_basic
    end type type_input_basic

    interface
        module subroutine read_simulation_settings(self, json)
            implicit none
            class(type_input_basic), intent(inout) :: self
            type(json_file), intent(inout) :: json
        end subroutine read_simulation_settings

        module subroutine read_analysis_controls(self, json)
            implicit none
            class(type_input_basic), intent(inout) :: self
            type(json_file), intent(inout) :: json
        end subroutine read_analysis_controls

        module subroutine read_geometry_settings(self, json)
            implicit none
            class(type_input_basic), intent(inout) :: self
            type(json_file), intent(inout) :: json
        end subroutine read_geometry_settings

        module subroutine read_materials(self, json)
            implicit none
            class(type_input_basic), intent(inout) :: self
            type(json_file), intent(inout) :: json
        end subroutine read_materials

        module subroutine read_solver_settings(self, json)
            implicit none
            class(type_input_basic), intent(inout) :: self
            type(json_file), intent(inout) :: json
        end subroutine read_solver_settings
    end interface

    interface
        module subroutine initialize_type_input_basic(self)
            implicit none
            class(type_input_basic), intent(inout) :: self
        end subroutine initialize_type_input_basic

        module subroutine display_input_basic(self)
            implicit none
            class(type_input_basic), intent(in) :: self
        end subroutine display_input_basic
    end interface

end module io_input_basic
