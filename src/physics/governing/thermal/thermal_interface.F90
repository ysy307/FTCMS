module physics_governing_thermal
    use, intrinsic :: iso_fortran_env
    use :: module_core
    use :: module_input, only:type_input, input_translator
    use :: module_constitutive, only:type_constitutive_manager
    implicit none
    private

    public :: type_thermal

    type :: type_thermal
        private
        integer(int32) :: computation_type
        integer(int32) :: computation_dimension
        logical :: enable_vapor_transport = .true.
        type(type_constitutive_manager) :: physics

        ! Enthalpy cache for BDF history terms
        ! Stores U(T_{n-k}) at each node for each history level (k >= 1)
        ! Avoids recomputing phase-change in compute_transient_term
    contains
        procedure, pass(self), public :: initialize => initialize_type_thermal
        procedure, pass(self), public :: initialize_from_configs => initialize_type_thermal_from_configs
        procedure, pass(self), public :: destroy => destroy_type_thermal

        procedure, pass(self), public :: compute_diffusion_term => compute_diffusion_term_thermal

        procedure, pass(self), public :: calc_enthalpy_density => calc_enthalpy_density_thermal
        procedure, pass(self), public :: calc_density_water => calc_density_water_thermal
        procedure, pass(self), public :: calc_density_ice => calc_density_ice_thermal
        procedure, pass(self), public :: calc_specific_heat_water => calc_specific_heat_water_thermal
        procedure, pass(self), public :: update_water_phases => update_water_phases_thermal
    end type type_thermal

    interface
        module subroutine calc_enthalpy_density_thermal(self, material_id, state, U, tangent)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: U
            real(real64), intent(inout), optional :: tangent(2)
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

        module subroutine update_water_phases_thermal(self, material_id, state)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state

        end subroutine update_water_phases_thermal

        module subroutine compute_diffusion_term_thermal(self, material_id, state, D_TT, isotropic_tangent)
            implicit none
            class(type_thermal), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state
            real(real64), intent(inout) :: D_TT(:, :)
            real(real64), intent(inout), optional :: isotropic_tangent(2)
            type(type_state_thc) :: lambda

        end subroutine compute_diffusion_term_thermal


        module subroutine calc_specific_heat_water_thermal(self, state, c_w)
            implicit none
            class(type_thermal), intent(in) :: self
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: c_w
        end subroutine calc_specific_heat_water_thermal
    end interface

    interface
        module subroutine initialize_type_thermal(self, input, active_region_ids)
            implicit none
            class(type_thermal), intent(inout) :: self
            type(type_input), intent(in) :: input
            integer(int32), intent(in) :: active_region_ids(:)

        end subroutine initialize_type_thermal

        !> Build the model from constitutive configurations rather than a project.
        !>
        !> The construction `initialize_type_thermal` performs once it has
        !> translated a `type_input`, exposed so that a caller holding the
        !> configurations directly -- a driver with no project directory -- reaches
        !> the same model through the same code. No default is invented for
        !> anything the project would otherwise supply.
        !>
        !> Assumes one configuration entry per active region, in the order of
        !> `active_region_ids`. Costs one constitutive-model construction per
        !> material and \(O(1)\) further memory.
        !>
        !> Fails by `error stop` when the configuration arrays and the region list
        !> have different lengths.
        module subroutine initialize_type_thermal_from_configs(self, active_region_ids, &
                                                               configs_density, configs_specific_heat, &
                                                               configs_heat_capacity, configs_thermal_conductivity, &
                                                               configs_wrf, configs_hcf, configs_gcc, &
                                                               computation_type, computation_dimension, &
                                                               enable_vapor_transport)
            implicit none
            class(type_thermal), intent(inout) :: self
            !> Region identifiers the material laws belong to, one per configuration.
            integer(int32), intent(in) :: active_region_ids(:)
            !> Phase densities of each region.
            type(type_config_constitutive), intent(in) :: configs_density(:)
            !> Phase specific heats of each region.
            type(type_config_constitutive), intent(in) :: configs_specific_heat(:)
            !> Volumetric heat capacities of each region.
            type(type_config_constitutive), intent(in) :: configs_heat_capacity(:)
            !> Phase thermal conductivities of each region.
            type(type_config_constitutive), intent(in) :: configs_thermal_conductivity(:)
            !> Retention curve of each region.
            type(type_config_wrf), intent(in) :: configs_wrf(:)
            !> Hydraulic conductivity function of each region.
            type(type_config_hcf), intent(in) :: configs_hcf(:)
            !> Freezing curve of each region.
            type(type_config_gcc), intent(in) :: configs_gcc(:)
            !> Computation type the coefficient tensors are laid out for.
            integer(int32), intent(in) :: computation_type
            !> Computation dimension of those tensors.
            integer(int32), intent(in) :: computation_dimension
            !> Whether the vapour terms take part in storage and transport.
            logical, intent(in) :: enable_vapor_transport
        end subroutine initialize_type_thermal_from_configs

        module subroutine destroy_type_thermal(self)
            implicit none
            class(type_thermal), intent(inout) :: self

        end subroutine destroy_type_thermal

    end interface

contains

end module physics_governing_thermal
