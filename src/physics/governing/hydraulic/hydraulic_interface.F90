module physics_governing_hydraulic
    use, intrinsic :: iso_fortran_env
    use :: module_core
    use :: module_control, only:type_control
    use :: module_input, only:type_input, input_translator
    use :: module_linalg
    use :: module_constitutive, g => gravity_acceleration, rho_std => reference_water_density
    use :: physics_governing_base, only:type_assemble_workspace, HYDRAULIC_DRIVER_TOTAL_POTENTIAL
    implicit none
    private

    public :: type_hydraulic
    public :: CAPACITY_REGULARIZATION_ENABLED

    !> Add the iteration-decaying L-scheme capacity to K_HH.
    !>
    !> This term is added to the MATRIX only and never to the residual, so
    !> while it is on, K is not the derivative of R by construction. That is
    !> tolerable for a Picard iteration - it is the classical L-scheme bound -
    !> but it invalidates the one identity a line search relies on:
    !> d/dalpha ||W R||^2 = -2||W R||^2 < 0 holds only for K = dR/du. With the
    !> element tangent now assembled by finite differences of the residual and
    !> the linear solve verified to 1e-8, this term is the last remaining
    !> inconsistency between the operator solved and the function measured.
    !>
    !> Removing it has now been measured TWICE and is worse both times, so the
    !> term is load bearing rather than optional. Without it the pressure row
    !> is degenerate where the medium saturates - dtheta_SWRC/ds_m -> 0 leaves
    !> only the compressible storage, ~2e-10 1/Pa - and the Newton direction
    !> for that row is meaningless. The first attempt (before any step bound)
    !> exploded with |du_p| ~ 1e13. The second, with the joint step bound in
    !> place, did not explode: the bound simply shrank the whole update to
    !> lambda ~ 1e-9, the iterate stopped moving (resT/0 pinned at 1.000,
    !> dQeff = 0), and the run died at t = 2.3 s having accepted one step. A
    !> trust region converts the explosion into paralysis; it does not supply
    !> the missing storage.
    !>
    !> The principled replacement is pseudo-transient continuation: the same
    !> matrix-only diagonal, but M/tau with tau driven by the measured residual
    !> reduction, so the added term vanishes as the iteration converges and K
    !> returns to dR/du in the limit. Until that is in place this stays on and
    !> the inconsistency it introduces is a known, bounded one (it decays as
    !> 1/k within a step).
    logical, parameter :: CAPACITY_REGULARIZATION_ENABLED = .false.

    type :: type_hydraulic
        private
        integer(int32) :: computation_type
        integer(int32) :: computation_dimension
        logical :: enable_vapor_transport = .true.
        !> Diagnostic-only switch: skip the pseudo-transient capacity
        !> regularization (see hydraulic_matrix.F90) when assembling K_HH.
        !> Set by the FD-Jacobian audit (report_fd_jacobian_ftcms) so its
        !> analytic comparison target is regularization-free, since that term
        !> is added to the matrix only and the FD residual cannot see it.
        !> Production assembly always runs with this .false.
        logical :: disable_capacity_regularization = .false.
        real(real64), allocatable :: iteration_capacity_bound(:)
        real(real64), allocatable :: specific_storage(:)
        type(type_constitutive_manager) :: physics
    contains
        procedure, pass(self), public :: initialize => initialize_type_hydraulic
        procedure, pass(self), public :: destroy => destroy_type_hydraulic

        ! --- Assembly Procedures ---
        procedure, pass(self), public :: assemble_local => assemble_local_hydraulic
        procedure, pass(self), private :: assemble_local_picard => assemble_local_picard_hydraulic

        ! --- Coefficient Computation Procedures ---
        procedure, pass(self), public :: compute_diffusion_term => compute_diffusion_term_hydraulic
        procedure, pass(self), public :: compute_diffusion_term_gas => compute_diffusion_term_gas_hydraulic
        procedure, pass(self), public :: compute_advective_term => compute_advective_term_hydraulic

        ! --- Coupling Coefficient Procedures ---
        procedure, pass(self), public :: compute_coupling_mass_term => compute_coupling_mass_term_hydraulic
        procedure, pass(self), public :: compute_coupling_diffusion_term => compute_coupling_diffusion_term_hydraulic

        ! --- Helper Procedures ---
        procedure, pass(self), public :: calc_cryo_head_dT => calc_cryo_head_dT_hydraulic
        procedure, pass(self), public :: calc_K_wT => calc_K_wT_hydraulic
        procedure, pass(self), public :: calc_K_wP => calc_K_wP_hydraulic
        procedure, pass(self), public :: calc_K_vT => calc_K_vT_hydraulic
        procedure, pass(self), public :: calc_K_vP => calc_K_vP_hydraulic
        procedure, pass(self), public :: is_vapor_transport_enabled => is_vapor_transport_enabled_hydraulic
        procedure, pass(self), public :: set_disable_capacity_regularization => &
            set_disable_capacity_regularization_hydraulic
        procedure, pass(self), public :: get_capacity_regularization => &
            get_capacity_regularization_hydraulic

        procedure, pass(self), public :: update_water_phases => update_water_phases_hydraulic
        procedure, pass(self), public :: calc_effective_density => calc_effective_density_hydraulic
        procedure, pass(self), public :: calc_effective_density_value => calc_effective_density_value_hydraulic

        procedure, pass(self), public :: compute_C_eq => compute_C_eq_hydraulic
        procedure, pass(self), public :: compute_total_storage => compute_total_storage_hydraulic
        procedure, pass(self), private :: compute_iteration_capacity => compute_iteration_capacity_hydraulic
        procedure, pass(self), private :: compute_transient_term_mixed => compute_transient_term_mixed_hydraulic
        procedure, pass(self), private :: compute_compressive_storage => compute_compressive_storage_hydraulic
        procedure, pass(self), public :: calc_segregation_sink => calc_segregation_sink_hydraulic
        ! ---- Meta / Utility ----
        !> TEMPORARY Step-0 verification diagnostic (see hydraulic_matrix.F90);
    end type type_hydraulic

    interface
        module subroutine initialize_type_hydraulic(self, input, active_region_ids)
            implicit none
            class(type_hydraulic), intent(inout) :: self
            type(type_input), intent(in) :: input
            integer(int32), intent(in) :: active_region_ids(:)
        end subroutine initialize_type_hydraulic

        module subroutine destroy_type_hydraulic(self)
            implicit none
            class(type_hydraulic), intent(inout) :: self
        end subroutine destroy_type_hydraulic

        ! --- Assembly Interfaces ---
        module subroutine assemble_local_hydraulic(self, control, workspace, K_HH, K_HT, F_H)
            implicit none
            class(type_hydraulic), intent(in) :: self
            type(type_control), intent(in) :: control
            type(type_assemble_workspace), intent(inout) :: workspace
            type(type_matrix_dense), intent(inout), optional :: K_HH
            type(type_matrix_dense), intent(inout), optional :: K_HT
            type(type_vector_dp), intent(inout), optional :: F_H
        end subroutine assemble_local_hydraulic

        module subroutine assemble_local_picard_hydraulic(self, control, workspace, K_HH, K_HT, F_H)
            implicit none
            class(type_hydraulic), intent(in) :: self
            type(type_control), intent(in) :: control
            type(type_assemble_workspace), intent(inout) :: workspace
            type(type_matrix_dense), intent(inout), optional :: K_HH
            type(type_matrix_dense), intent(inout), optional :: K_HT
            type(type_vector_dp), intent(inout), optional :: F_H
        end subroutine assemble_local_picard_hydraulic

        ! --- Physics Term Interfaces ---
        module subroutine compute_diffusion_term_hydraulic(self, material_id, state, D_HH)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state
            real(real64), intent(inout) :: D_HH(:, :)
        end subroutine compute_diffusion_term_hydraulic

        module subroutine compute_diffusion_term_gas_hydraulic(self, material_id, state, D_HH_gas)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state
            real(real64), intent(inout) :: D_HH_gas(:, :)
        end subroutine compute_diffusion_term_gas_hydraulic

        module subroutine compute_advective_term_hydraulic(self, material_id, state, V_H)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state
            real(real64), intent(inout) :: V_H(:)
        end subroutine compute_advective_term_hydraulic

        ! --- Coupling Coefficient Interfaces ---
        module subroutine compute_coupling_mass_term_hydraulic(self, material_id, state, C_HT)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: C_HT
        end subroutine compute_coupling_mass_term_hydraulic

        module subroutine compute_coupling_diffusion_term_hydraulic(self, material_id, state, D_HT)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state
            real(real64), intent(inout) :: D_HT(:, :)
        end subroutine compute_coupling_diffusion_term_hydraulic

        ! --- Helper Interfaces ---
        !> Temperature sensitivity of the total-potential head, dh/dT [m/K].
        module subroutine calc_cryo_head_dT_hydraulic(self, material_id, state, dh_dT)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: dh_dT
        end subroutine calc_cryo_head_dT_hydraulic

        module subroutine calc_K_wT_hydraulic(self, target_id, state, K_wT)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: target_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: K_wT
        end subroutine calc_K_wT_hydraulic

        module subroutine calc_K_wP_hydraulic(self, target_id, state, K_wP)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: target_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: K_wP
        end subroutine calc_K_wP_hydraulic

        module subroutine calc_K_vT_hydraulic(self, target_id, state, K_vT)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: target_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: K_vT
        end subroutine calc_K_vT_hydraulic

        module subroutine calc_K_vP_hydraulic(self, target_id, state, K_vP)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: target_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: K_vP
        end subroutine calc_K_vP_hydraulic

        module pure function is_vapor_transport_enabled_hydraulic(self) result(enabled)
            implicit none
            class(type_hydraulic), intent(in) :: self
            logical :: enabled
        end function is_vapor_transport_enabled_hydraulic

        !> Enable/disable the pseudo-transient capacity regularization added to
        !> K_HH by assemble_local_picard_hydraulic. Diagnostic use only (see
        !> the disable_capacity_regularization member doc).
        module subroutine set_disable_capacity_regularization_hydraulic(self, disable)
            implicit none
            class(type_hydraulic), intent(inout) :: self
            logical, intent(in) :: disable
        end subroutine set_disable_capacity_regularization_hydraulic

        !> Value of the pseudo-transient capacity regularization
        !> assemble_local_picard_hydraulic adds to K_HH (hydraulic_matrix.F90),
        !> without needing to assemble K_HH itself.
        !>
        !> \[ C_{reg} = C_{bound}(m) / \max(1, k) \]
        !> where \(k\) is the nonlinear iteration. Returns 0 if the bound array
        !> is not allocated or material_id is out of range, matching
        !> assemble_local_picard_hydraulic's own guard. Reads only the
        !> per-material bound array set once at initialize, so this is safe to
        !> call concurrently from every thread of the parallel element loop.
        module pure function get_capacity_regularization_hydraulic(self, material_id, nonlinear_iteration) &
            result(value)
            implicit none
            class(type_hydraulic), intent(in) :: self
            !> Element material id
            integer(int32), intent(in) :: material_id
            !> Current nonlinear (Newton) iteration, >= 1
            integer(int32), intent(in) :: nonlinear_iteration
            real(real64) :: value
        end function get_capacity_regularization_hydraulic

        module subroutine update_water_phases_hydraulic(self, material_id, state)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(inout) :: state
        end subroutine update_water_phases_hydraulic

        module subroutine calc_effective_density_hydraulic(self, material_id, state, bdf_coeffs, drho_dt)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(in) :: bdf_coeffs(:)
            real(real64), intent(inout) :: drho_dt
        end subroutine calc_effective_density_hydraulic

        !> Evaluate the pore-water effective density at the supplied state.
        !> \( \rho_{eff} = \rho_w \theta_w + \rho_{ice} \theta_{ice} + \rho_w \theta_v^{\star} \) [kg/m3].
        !> Plain (non-time-derivative) conserved storage quantity of the water-mass
        !> balance; phase contents must already be consistent with (T, p_w).
        module subroutine calc_effective_density_value_hydraulic(self, material_id, state, rho_eff)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: rho_eff
        end subroutine calc_effective_density_value_hydraulic

        !> Compute equivalent specific moisture capacity C_eq = dTheta/dP [1/Pa].
        !> \[ C_{eq} = \frac{\partial\Theta}{\partial P} + \frac{\partial\theta_v}{\partial P}
        !>           + C_{compressive} \]
        !> \(\Theta=\theta_{SWRC}(s_m)\) is the retention-curve total water content;
        !> its pressure tangent is read directly from the state (set by
        !> calc_phase_split), not recombined from per-phase derivatives. The vapor
        !> term applies only when vapor transport is enabled, and \(C_{compressive}\)
        !> is the saturated compressible-storage capacity.
        module subroutine compute_C_eq_hydraulic(self, material_id, state, C_eq)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: C_eq
        end subroutine compute_C_eq_hydraulic

        module subroutine compute_iteration_capacity_hydraulic(self, material_id, state, capacity)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: capacity
        end subroutine compute_iteration_capacity_hydraulic

        !> Total stored water per unit volume at one state,
        !> \[ \Theta = \theta_w + \frac{\rho_i}{\rho_w}\theta_i + \theta_v + S_c \]
        !> where \( S_c \) is the compressive storage. This is the quantity
        !> compute_C_eq is the pressure derivative of, so every level that
        !> differences stored water must difference this one and not a
        !> reassembled sum of phases: dropping \( S_c \) here while keeping its
        !> capacity in compute_C_eq leaves a matrix diagonal that no residual
        !> term corresponds to, which is exactly where the medium saturates and
        !> \( d\theta_{SWRC}/dp \) vanishes.
        module subroutine compute_total_storage_hydraulic(self, material_id, state, Theta)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: Theta
        end subroutine compute_total_storage_hydraulic

        !> Compute BDF approximation of dTheta/dt for Mixed formulation.
        !> \[ \frac{d\Theta}{dt} \approx \sum_j \alpha_j \Theta(t_{n+1-j}) \]
        !> where \(\Theta = \theta_l + (\rho_i/\rho_l)\theta_i + \theta_v\).
        module subroutine compute_transient_term_mixed_hydraulic(self, material_id, state, bdf_coeffs, dTheta_dt)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(in) :: bdf_coeffs(:)
            real(real64), intent(inout) :: dTheta_dt
        end subroutine compute_transient_term_mixed_hydraulic

        module subroutine compute_compressive_storage_hydraulic(self, material_id, state, storage, capacity)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(inout) :: storage
            real(real64), intent(inout), optional :: capacity
        end subroutine compute_compressive_storage_hydraulic

        module subroutine calc_segregation_sink_hydraulic(self, material_id, state, dt, S_seg)
            implicit none
            class(type_hydraulic), intent(in) :: self
            integer(int32), intent(in) :: material_id
            type(type_state), intent(in) :: state
            real(real64), intent(in) :: dt
            real(real64), intent(inout) :: S_seg
        end subroutine calc_segregation_sink_hydraulic

    end interface

contains
    module subroutine destroy_type_hydraulic(self)
        implicit none
        class(type_hydraulic), intent(inout) :: self

        if (allocated(self%iteration_capacity_bound)) deallocate (self%iteration_capacity_bound)
        if (allocated(self%specific_storage)) deallocate (self%specific_storage)
    end subroutine destroy_type_hydraulic
end module physics_governing_hydraulic
