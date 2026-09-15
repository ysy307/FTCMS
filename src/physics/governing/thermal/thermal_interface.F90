module physics_governing_thermal
    use, intrinsic :: iso_fortran_env
    use :: module_core
    use :: module_control, only:type_control
    use :: module_input, only:type_input, input_translator
    use :: module_linalg
    use :: module_constitutive, only:type_constitutive_manager
    use :: physics_governing_base, only:type_assemble_workspace, HYDRAULIC_DRIVER_TOTAL_POTENTIAL
    implicit none
    private

    public :: type_thermal
    public :: SCHEME_TANGENT

    integer, parameter :: SCHEME_TANGENT = 0 ! Exact derivative (for Newton-Raphson)
    integer, parameter :: SCHEME_SECANT = 1 ! Finite difference (for Picard / phase-change stabilization)

    type :: type_thermal
        private
        integer(int32) :: computation_type
        integer(int32) :: computation_dimension
        logical :: enable_vapor_transport = .true.
        type(type_constitutive_manager) :: physics

        ! Enthalpy cache for BDF history terms
        ! Stores U(T_{n-k}) at each node for each history level (k >= 1)
        ! Avoids recomputing phase-change in compute_transient_term
        real(real64), allocatable :: enthalpy_cache(:, :)
        logical :: enthalpy_cache_valid = .false.
        integer(int32) :: enthalpy_cache_num_nodes = 0
        integer(int32) :: enthalpy_cache_num_hist = 0
    contains
        procedure, pass(self), public :: initialize => initialize_type_thermal
        procedure, pass(self), public :: destroy => destroy_type_thermal

        procedure, pass(self), public :: assemble_local => assemble_local_thermal
        procedure, pass(self), private :: assemble_local_picard => assemble_local_picard_thermal

        procedure, pass(self), public :: compute_mass_term => compute_mass_term_thermal
        procedure, pass(self), public :: compute_apparent_capacity => compute_apparent_capacity_thermal
        procedure, pass(self), public :: compute_diffusion_term => compute_diffusion_term_thermal
        procedure, pass(self), private :: compute_advective_term => compute_advective_term_thermal
        procedure, pass(self), private :: compute_latent_term => compute_latent_term_thermal
        procedure, pass(self), private :: compute_transient_term => compute_transient_term_thermal
        procedure, pass(self), private :: compute_history_term => compute_history_term_thermal
        procedure, pass(self), public :: compute_coupling_mass_term => compute_coupling_mass_term_thermal

        procedure, pass(self), public :: calc_enthalpy_density => calc_enthalpy_density_thermal
        procedure, pass(self), public :: calc_density_water => calc_density_water_thermal
        procedure, pass(self), public :: calc_density_ice => calc_density_ice_thermal
        procedure, pass(self), public :: calc_density_vapor_saturation => calc_density_vapor_saturation_thermal
        procedure, pass(self), public :: calc_specific_heat_water => calc_specific_heat_water_thermal
        procedure, pass(self), public :: calc_specific_heat_vapor => calc_specific_heat_vapor_thermal
        procedure, pass(self), public :: update_water_phases => update_water_phases_thermal
        procedure, pass(self), public :: project_ice_content => project_ice_content_thermal
        procedure, pass(self), public :: calc_conserved_target => calc_conserved_target_thermal
        procedure, pass(self), public :: solve_local_conserved_equilibrium => solve_local_conserved_equilibrium_thermal
        procedure, pass(self), public :: cache_enthalpy_history => cache_enthalpy_history_thermal
        procedure, pass(self), public :: invalidate_enthalpy_cache => invalidate_enthalpy_cache_thermal
    end type type_thermal

    interface
        module subroutine calc_enthalpy_density_thermal(self, material_id, state, U)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: U
        end subroutine calc_enthalpy_density_thermal

        module subroutine calc_density_water_thermal(self, state, rho_water)
            implicit none
            class(type_thermal), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: rho_water

        end subroutine calc_density_water_thermal

        module subroutine calc_density_ice_thermal(self, state, rho_ice)
            implicit none
            class(type_thermal), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: rho_ice

        end subroutine calc_density_ice_thermal

        module subroutine calc_density_vapor_saturation_thermal(self, state, rho_vapor_sat)
            implicit none
            class(type_thermal), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: rho_vapor_sat

        end subroutine calc_density_vapor_saturation_thermal

        module subroutine calc_specific_heat_water_thermal(self, state, c_w)
            implicit none
            class(type_thermal), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: c_w
        end subroutine calc_specific_heat_water_thermal

        module subroutine calc_specific_heat_vapor_thermal(self, state, c_v)
            implicit none
            class(type_thermal), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: c_v
        end subroutine calc_specific_heat_vapor_thermal

        module subroutine update_water_phases_thermal(self, material_id, state)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state

        end subroutine update_water_phases_thermal

        module subroutine project_ice_content_thermal(self, material_id, state, projected_ice, &
                                                       ice_increment, equilibrium_error, active_bound, &
                                                       dice_dT, dice_dP)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: projected_ice
            real(real64), intent(inout) :: ice_increment
            real(real64), intent(inout) :: equilibrium_error
            integer(int32), intent(inout), optional :: active_bound
            real(real64), intent(inout), optional :: dice_dT
            real(real64), intent(inout), optional :: dice_dP
        end subroutine project_ice_content_thermal

        module subroutine calc_conserved_target_thermal(self, material_id, state, target_total_water, available)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: target_total_water
            logical, intent(inout) :: available
        end subroutine calc_conserved_target_thermal

        module subroutine solve_local_conserved_equilibrium_thermal(self, material_id, state, target_total_water, &
                                                                     new_pressure, new_ice, converged)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(in) :: target_total_water
            real(real64), intent(inout) :: new_pressure
            real(real64), intent(inout) :: new_ice
            logical, intent(inout) :: converged
        end subroutine solve_local_conserved_equilibrium_thermal

        module subroutine compute_transient_term_thermal(self, material_id, state, bdf_coeffs, dU_dt)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(in) :: bdf_coeffs(:)
            real(real64), intent(inout) :: dU_dt

        end subroutine compute_transient_term_thermal
        module subroutine compute_mass_term_thermal(self, material_id, state, C_TT, scheme_opt)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: C_TT
            integer(int32), intent(in), optional :: scheme_opt

        end subroutine compute_mass_term_thermal

        module subroutine compute_apparent_capacity_thermal(self, material_id, state, C_app)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: C_app

        end subroutine compute_apparent_capacity_thermal

        module subroutine compute_diffusion_term_thermal(self, material_id, state, D_TT)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state
            real(real64), intent(inout) :: D_TT(:, :)
            type(type_state_thc) :: lambda

        end subroutine compute_diffusion_term_thermal

        module subroutine compute_advective_term_thermal(self, material_id, state, V_TT)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state
            real(real64), intent(inout) :: V_TT(:)

        end subroutine compute_advective_term_thermal

        module subroutine compute_latent_term_thermal(self, material_id, state, L_TT)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state
            real(real64), intent(inout) :: L_TT

        end subroutine compute_latent_term_thermal

        module subroutine compute_history_term_thermal(self, material_id, state, bdf_coeffs, history_term)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(in) :: bdf_coeffs(:)
            real(real64), intent(inout) :: history_term

        end subroutine compute_history_term_thermal

        module subroutine compute_coupling_mass_term_thermal(self, material_id, state, C_TH)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: C_TH
        end subroutine compute_coupling_mass_term_thermal

        module subroutine cache_enthalpy_history_thermal(self, num_nodes, num_hist, material_ids, &
                                                         temperature_all, pressure_all, porosity_all)
            implicit none
            class(type_thermal), intent(inout) :: self
            integer(int32), intent(in) :: num_nodes
            integer(int32), intent(in) :: num_hist
            integer(int32), intent(in) :: material_ids(:)
            real(real64), intent(in) :: temperature_all(:, :)
            real(real64), intent(in) :: pressure_all(:, :)
            real(real64), intent(in) :: porosity_all(:, :)

        end subroutine cache_enthalpy_history_thermal

        module subroutine invalidate_enthalpy_cache_thermal(self)
            implicit none
            class(type_thermal), intent(inout) :: self

        end subroutine invalidate_enthalpy_cache_thermal
    end interface

    interface
        module subroutine initialize_type_thermal(self, input, active_region_ids)
            implicit none
            class(type_thermal), intent(inout) :: self
            type(type_input), intent(in) :: input
            integer(int32), intent(in) :: active_region_ids(:)

        end subroutine initialize_type_thermal

        module subroutine destroy_type_thermal(self)
            implicit none
            class(type_thermal), intent(inout) :: self

        end subroutine destroy_type_thermal

        module subroutine assemble_local_thermal(self, control, workspace, K_TT, K_TH, F_T)
            implicit none
            class(type_thermal), intent(in) :: self
            type(type_control), intent(in) :: control
            type(type_assemble_workspace), intent(inout) :: workspace
            type(type_matrix_dense), intent(inout), optional :: K_TT
            type(type_matrix_dense), intent(inout), optional :: K_TH
            type(type_vector_dp), intent(inout), optional :: F_T

        end subroutine assemble_local_thermal

        module subroutine assemble_local_picard_thermal(self, control, workspace, K_TT, K_TH, F_T)
            implicit none
            class(type_thermal), intent(in) :: self
            type(type_control), intent(in) :: control
            type(type_assemble_workspace), intent(inout) :: workspace
            type(type_matrix_dense), intent(inout), optional :: K_TT
            type(type_matrix_dense), intent(inout), optional :: K_TH
            type(type_vector_dp), intent(inout), optional :: F_T

        end subroutine assemble_local_picard_thermal

    end interface

contains

end module physics_governing_thermal
