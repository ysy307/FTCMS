module physics_governing_hydraulic
    use, intrinsic :: iso_fortran_env
    use :: module_core
    use :: module_input, only:type_input, input_translator
    use :: module_constitutive, g => gravity_acceleration, rho_std => reference_water_density
    implicit none
    private

    public :: type_hydraulic

    type :: type_hydraulic
        private
        integer(int32) :: computation_type
        integer(int32) :: computation_dimension
        logical :: enable_vapor_transport = .true.
        real(real64), allocatable :: specific_storage(:)
        type(type_constitutive_manager) :: physics
    contains
        procedure, pass(self), public :: initialize => initialize_type_hydraulic
        procedure, pass(self), public :: initialize_from_configs => initialize_type_hydraulic_from_configs
        procedure, pass(self), public :: destroy => destroy_type_hydraulic

        ! --- Coefficient Computation Procedures ---
        procedure, pass(self), public :: compute_diffusion_term => compute_diffusion_term_hydraulic

        ! --- Helper Procedures ---
        procedure, pass(self), public :: is_vapor_transport_enabled => is_vapor_transport_enabled_hydraulic

        procedure, pass(self), public :: update_water_phases => update_water_phases_hydraulic

        procedure, pass(self), public :: compute_total_storage => compute_total_storage_hydraulic
        procedure, pass(self), private :: compute_compressive_storage => compute_compressive_storage_hydraulic
    end type type_hydraulic

    interface
        module subroutine initialize_type_hydraulic(self, input, active_region_ids)
            implicit none
            class(type_hydraulic), intent(inout) :: self
            type(type_input), intent(in) :: input
            integer(int32), intent(in) :: active_region_ids(:)
        end subroutine initialize_type_hydraulic

        !> Build the model from constitutive configurations rather than a project.
        !>
        !> The construction `initialize_type_hydraulic` performs once it has
        !> translated a `type_input`, exposed so that a caller holding the
        !> configurations directly -- a driver with no project directory -- reaches
        !> the same model through the same code. Every material law, the L-scheme
        !> capacity bound and the specific storage are those of the arguments; no
        !> default is invented for anything the project would otherwise supply.
        !>
        !> Assumes one configuration entry per active region, in the order of
        !> `active_region_ids`, and region identifiers that are positive. Costs one
        !> constitutive-model construction per material and \(O(\max_i \mathrm{id}_i)\)
        !> memory for the per-material arrays.
        !>
        !> Fails by `error stop` when the configuration arrays and the region list
        !> have different lengths.
        module subroutine initialize_type_hydraulic_from_configs(self, active_region_ids, &
                                                                 configs_density, configs_wrf, configs_hcf, configs_gcc, &
                                                                 specific_storage, computation_type, computation_dimension, &
                                                                 enable_vapor_transport)
            implicit none
            class(type_hydraulic), intent(inout) :: self
            !> Region identifiers the material laws belong to, one per configuration.
            integer(int32), intent(in) :: active_region_ids(:)
            !> Phase densities of each region.
            type(type_config_constitutive), intent(in) :: configs_density(:)
            !> Retention curve of each region.
            type(type_config_wrf), intent(in) :: configs_wrf(:)
            !> Hydraulic conductivity function of each region.
            type(type_config_hcf), intent(in) :: configs_hcf(:)
            !> Freezing curve of each region.
            type(type_config_gcc), intent(in) :: configs_gcc(:)
            !> Specific storage \(S_s\) of each region, 1/m; non-negative.
            real(real64), intent(in) :: specific_storage(:)
            !> Computation type the coefficient tensors are laid out for.
            integer(int32), intent(in) :: computation_type
            !> Computation dimension of those tensors.
            integer(int32), intent(in) :: computation_dimension
            !> Whether the vapour terms take part in storage and transport.
            logical, intent(in) :: enable_vapor_transport
        end subroutine initialize_type_hydraulic_from_configs

        module subroutine destroy_type_hydraulic(self)
            implicit none
            class(type_hydraulic), intent(inout) :: self
        end subroutine destroy_type_hydraulic

        ! --- Physics Term Interfaces ---
        module subroutine compute_diffusion_term_hydraulic(self, material_id, state, D_HH, log_conductivity_tangent)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state
            real(real64), intent(inout) :: D_HH(:, :)
            !> \((\partial\ln K_{fl h}/\partial p_w,\partial\ln K_{fl h}/\partial T)\) of the conductivity
            !> \(D_{HH}\) carries, through the published phase split; overwritten on exit when present
            real(real64), intent(inout), optional :: log_conductivity_tangent(2)
        end subroutine compute_diffusion_term_hydraulic

        module pure function is_vapor_transport_enabled_hydraulic(self) result(enabled)
            implicit none
            class(type_hydraulic), intent(in) :: self
            logical :: enabled
        end function is_vapor_transport_enabled_hydraulic

        module subroutine update_water_phases_hydraulic(self, material_id, state)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state
        end subroutine update_water_phases_hydraulic

        !> Total stored water per unit volume at one state,
        !> \[ \Theta = \theta_w + \frac{\rho_i}{\rho_w}\theta_i + \theta_v + S_c \]
        !> where \( S_c \) is the compressive storage. This is the quantity
        !> compute_C_eq is the pressure derivative of, so every level that
        !> differences stored water must difference this one and not a
        !> reassembled sum of phases: dropping \( S_c \) here while keeping its
        !> capacity in compute_C_eq leaves a matrix diagonal that no residual
        !> term corresponds to, which is exactly where the medium saturates and
        !> \( d\theta_{SWRC}/dp \) vanishes.
        module subroutine compute_total_storage_hydraulic(self, material_id, state, Theta, tangent)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: Theta
            !> \((\partial\Theta/\partial p_w,\partial\Theta/\partial T)\), exact through the tangents the
            !> phase split published on `state`; overwritten on exit when present
            real(real64), intent(inout), optional :: tangent(2)
        end subroutine compute_total_storage_hydraulic

        module subroutine compute_compressive_storage_hydraulic(self, material_id, state, storage, tangent)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: storage
            !> \((\partial S_c/\partial p_w,\partial S_c/\partial T)\), including the motion of the saturation
            !> pressure; overwritten on exit when present
            real(real64), intent(inout), optional :: tangent(2)
        end subroutine compute_compressive_storage_hydraulic

    end interface

contains
    module subroutine destroy_type_hydraulic(self)
        implicit none
        class(type_hydraulic), intent(inout) :: self

        if (allocated(self%specific_storage)) deallocate (self%specific_storage)
    end subroutine destroy_type_hydraulic
end module physics_governing_hydraulic
