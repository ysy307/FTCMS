submodule(physics_governing_hydraulic) hydraulic_matrix
    implicit none


    integer(int32), save :: dbg_elements_considered = 0
    integer(int32), save :: dbg_elements_cut = 0
    real(real64), save :: dbg_phi_min = huge(1.0d0)
    real(real64), save :: dbg_phi_max = -huge(1.0d0)
    real(real64), save :: dbg_abs_phi_min = huge(1.0d0)

    !> Report how far the element conductivity can exceed the series-flow value.
    !>
    !> The element takes the GEOMETRIC mean of its nodal conductivities (see
    !> the comment above D_HH_elem). For flow in series through a nearly
    !> impermeable frozen layer the physical weighting is the HARMONIC mean, so
    !> geo/har bounds how much the discrete operator can over-conduct across
    !> the freezing front. Measured from the implementation rather than from a
    !> re-derivation, because k_r has two branches (closed form only when
    !> m = 1 - 1/n, incomplete beta otherwise) and the wrong one is 26x off in
    !> the frozen range.
    logical, parameter :: KELEM_REPORT = .false.
    !> Weight the element conductivity with the HARMONIC mean of its nodal
    !> values instead of the geometric mean. Flow across a freezing front is in
    !> series through a nearly impermeable layer, for which the harmonic mean
    !> is the physical weighting; measured, the geometric mean the element uses
    !> exceeds it by a factor 4 to 14 at the front. Kept behind a switch
    !> because the geometric mean is the standard internodal weight for an
    !> exponentially varying coefficient away from such an interface.
    logical, parameter :: ELEMENT_CONDUCTIVITY_HARMONIC = .false.
    !> Only elements whose nodal conductivities span at least this factor.
    real(real64), parameter :: KELEM_REPORT_MIN_CONTRAST = 1.0d1
    integer(int32), parameter :: KELEM_REPORT_LIMIT = 40
    integer(int32), save :: kelem_report_count = 0

contains

    !> @brief Assemble local matrix and vector (Picard only)
    module subroutine assemble_local_hydraulic(self, control, workspace, K_HH, K_HT, F_H)
        implicit none
        class(type_hydraulic), intent(in) :: self
        type(type_control), intent(in) :: control
        type(type_assemble_workspace), intent(inout) :: workspace
        type(type_matrix_dense), intent(inout), optional :: K_HH
        type(type_matrix_dense), intent(inout), optional :: K_HT
        type(type_vector_dp), intent(inout), optional :: F_H

        call self%assemble_local_picard(control, workspace, K_HH, K_HT, F_H)

    end subroutine assemble_local_hydraulic

    !> @brief Assemble Picard local components (backward Euler, no BDF history)
    !>
    !> All coefficients are evaluated directly at quadrature points from the
    !> interpolated state (no nodal pre-evaluation / lerp of coefficients).
    !>
    !> The transport coefficients can change sharply across the freezing
    !> interface, where ice first appears (see the cut-detection block below
    !> for the level-set definition).  Elements cut by the
    !> (build_interface_quadrature_points), which REPLACES the standard Gauss
    !> rule for the pressure, temperature, and gravity fluxes in those elements.
    !> Using one quadrature partition for all three parts of the Darcy flux
    !> avoids a discrete imbalance caused solely by sampling the same flux at
    !> different points.
    !>
    !> The storage term (C_eq -> K1, and its residual counterpart d(Theta)/dt)
    !> q_s loop as the fluxes: dTheta/dp collapses across the same front the
    !> flux coefficients do, so sampling it at Gauss points would reintroduce
    !> The segregation sink keeps the standard Gauss rule in all elements.
    !>
    !> In Modified Picard, the transport coefficients are frozen at the current
    !> iterate while both primary increments remain coupled. The K_HT block is
    !> therefore K2(D_HT), matching the temperature-flux term in the residual.
    module subroutine assemble_local_picard_hydraulic(self, control, workspace, K_HH, K_HT, F_H)
        implicit none
        class(type_hydraulic), intent(in) :: self
        type(type_control), intent(in) :: control
        type(type_assemble_workspace), intent(inout) :: workspace
        type(type_matrix_dense), intent(inout), optional :: K_HH
        type(type_matrix_dense), intent(inout), optional :: K_HT
        type(type_vector_dp), intent(inout), optional :: F_H

        integer(int32) :: i, j, d, n_nodes, n_gauss, n_dim, ierr
        real(real64) :: bdf0, dt_local

        real(real64) :: xi_sub, eta_sub, weight_sub
        logical :: is_plus_side_sub
        real(real64) :: phi_nodes(workspace%num_fe_nodes)
        real(real64) :: porosity_nodes(workspace%num_fe_nodes)
        real(real64) :: T_q_sub, P_q_sub, porosity_q_sub, Qi_q_sub
        real(real64) :: D_HH_sub, D_HT_sub, eff_weight_sub, det_J_sub
        real(real64) :: V_sub(workspace%num_fe_dimension)
        real(real64) :: vec_V_sub(workspace%num_fe_nodes)
        real(real64) :: coeff_sub_mat(workspace%num_fe_dimension, workspace%num_fe_dimension)
        real(real64) :: coeff_gas_mat(workspace%num_fe_dimension, workspace%num_fe_dimension)
        real(real64) :: mat_HH_sub(workspace%num_fe_nodes, workspace%num_fe_nodes)
        real(real64) :: mat_HT_sub(workspace%num_fe_nodes, workspace%num_fe_nodes)
        real(real64) :: dpsi_dx_sub(workspace%num_fe_dimension, workspace%num_fe_nodes)
        type(type_coordinate_dp) :: r_sub
        type(type_state) :: state_sub
        !     sub-element front as the flux coefficients above) ---
        real(real64) :: psi_sub(workspace%num_fe_nodes)
        real(real64) :: C_eq_sub, C_HT_storage_sub, d_dt_sub
        real(real64) :: mat_C_sub(workspace%num_fe_nodes, workspace%num_fe_nodes)
        real(real64) :: mat_C_HT_sub(workspace%num_fe_nodes, workspace%num_fe_nodes)
        real(real64) :: vec_ddt_sub(workspace%num_fe_nodes)
        integer(int32) :: k, n_hist, hlen
        real(real64) :: T_hist_node(workspace%bdf_order + 1, workspace%num_fe_nodes)
        real(real64) :: P_hist_node(workspace%bdf_order + 1, workspace%num_fe_nodes)
        real(real64) :: Qi_hist_node(workspace%bdf_order + 1, workspace%num_fe_nodes)
        real(real64) :: T_hist_sub(workspace%bdf_order + 1)
        real(real64) :: P_hist_sub(workspace%bdf_order + 1)
        real(real64) :: Qi_hist_sub(workspace%bdf_order + 1)
        real(real64), pointer, contiguous, dimension(:) :: hist_ptr
        logical :: hist_is_set

        ! --- Standard assembly variables ---
        real(real64) :: local_vec_res(workspace%num_fe_nodes)
        real(real64) :: work_C_HT(workspace%num_fe_gauss)
        real(real64) :: work_sink(workspace%num_fe_gauss)
        real(real64) :: work_D_HT(workspace%num_fe_dimension, workspace%num_fe_dimension, workspace%num_fe_gauss)
        real(real64) :: work_matrix_coupling(workspace%num_fe_nodes, workspace%num_fe_nodes)
        real(real64) :: D_HT_tmp(workspace%num_fe_dimension, workspace%num_fe_dimension)
        ! Nodal transport coefficients and the element values formed from them.
        real(real64) :: D_HH_node(workspace%num_fe_nodes), D_HT_node(workspace%num_fe_nodes)
        real(real64) :: D_HH_gas_node(workspace%num_fe_nodes)
        real(real64) :: V_node_mag(workspace%num_fe_nodes)
        real(real64) :: D_HH_elem, D_HT_elem, V_elem_mag, D_HH_gas_elem
        ! Vapour flux operator: same K2 machinery as the liquid, but driven by
        ! the gas pressure instead of the liquid potential.
        real(real64) :: work_D_gas(workspace%num_fe_dimension, workspace%num_fe_dimension, workspace%num_fe_gauss)
        real(real64) :: work_matrix_gas(workspace%num_fe_nodes, workspace%num_fe_nodes)
        real(real64) :: mat_HH_gas_sub(workspace%num_fe_nodes, workspace%num_fe_nodes)
        real(real64) :: D_HH_gas_sub
        logical :: gas_flux_needed
        real(real64) :: V_dir(workspace%num_fe_dimension), V_node_tmp(workspace%num_fe_dimension)
        ! Nodal Darcy driver and its tangents in the two primary unknowns.
        real(real64) :: P_gen_node(workspace%num_fe_nodes), psi_eff_i
        real(real64) :: dPgen_dP_node(workspace%num_fe_nodes), dSeff_dP_i
        real(real64) :: dPgen_dT_node(workspace%num_fe_nodes), dSeff_dT_i
        integer(int32) :: nonlinear_iteration
        real(real64) :: capacity_regularization
        logical :: thermal_target, coupling_mass_needed, coupling_flux_needed
        logical :: coupling_block_needed

        n_nodes = workspace%num_fe_nodes
        n_gauss = workspace%num_fe_gauss
        n_dim = workspace%num_fe_dimension
        thermal_target = control%is_target(PHYSICS_TYPES%THERMAL, workspace%material_id)
        coupling_mass_needed = present(K_HT) .and. thermal_target
        coupling_flux_needed = (present(F_H) .or. present(K_HT)) .and. thermal_target
        coupling_block_needed = present(K_HT) .and. thermal_target
        gas_flux_needed = self%enable_vapor_transport

        bdf0 = workspace%bdf_coeffs(1)
        dt_local = 0.0d0
        call control%get_dt(dt_local)

        workspace%work_C(:) = 0.0d0
        workspace%work_D(:, :, :) = 0.0d0
        workspace%work_V(:, :) = 0.0d0
        workspace%work_d_dt(:) = 0.0d0
        local_vec_res(:) = 0.0d0
        work_C_HT(:) = 0.0d0
        work_sink(:) = 0.0d0
        work_D_HT(:, :, :) = 0.0d0

        ! ----------------------------------------------------------------
        ! 0. Cut detection: nodal level set phi = (s_f - s_m) + eps_s, the
        !    SAME switch calc_phase_split uses to decide ice existence
        !    (s_m = -p_w, s_f = the raw generalized-Clapeyron suction; see
        !    type_fusion%calc_freezing_level_set for the full derivation and
        !    why phi>0 iff ice_content>0 exactly, given the compact-support
        !    smoothing in fusion.F90). Computing phi from the split's own
        !    s_m,s_f - rather than from a separately defined critical
        !    temperature T_high(p_w) - is what keeps "cut" elements from
        !    drifting away from where the constitutive split actually places
        !    the phase change.
        ! ----------------------------------------------------------------
        ! 1. Gauss loop: storage and sink terms at all Gauss points;
        !    flux coefficients only for uncut elements.
        ! ----------------------------------------------------------------
        ! Pseudo-transient regularization of the pressure block.
        !
        ! At the freezing front the storage capacity dTheta/dp collapses on the
        ! pore-volume bound and the conductivity collapses under the impedance,
        ! so the pressure row loses both of its terms at once and the Jacobian
        ! becomes nearly singular. The Newton step is then unbounded - measured,
        ! max|du_p| = 2.8e5 Pa against an ambient 5.4e4 Pa, with the residual
        ! growing superlinearly in the step length, which is what left the line
        ! search with no admissible alpha.
        !
        ! The time-derivative term cannot regularize this: it enters as
        ! bdf0*K1(C_eq), so when C_eq itself vanishes the row stays singular no
        ! matter how small dt becomes. That is why cutting dt never helped -
        ! measured, dt = 26.9 s and 13.4 s failed where 73.4 s had succeeded.
        !
        ! A capacity is therefore added to the MATRIX only, decaying as 1/k with
        ! the nonlinear iteration so the operator starts strictly diagonally
        ! reinforced and approaches the true tangent as the iteration converges.
        ! The residual keeps the exact storage term, so the converged solution
        ! is untouched; only the path to it is bounded.
        call control%get_nonlinear_iter(nonlinear_iteration)
        capacity_regularization = 0.0d0
        ! Diagnostic-only bypass for the FD-Jacobian audit (see
        ! disable_capacity_regularization member doc): the term below is added
        ! to the MATRIX only, never the residual, so the audit's FD comparison
        ! target must be assembled without it to be comparable.
        if (CAPACITY_REGULARIZATION_ENABLED .and. .not. self%disable_capacity_regularization) then
            if (allocated(self%iteration_capacity_bound)) then
                if (workspace%material_id >= 1 .and. &
                    workspace%material_id <= size(self%iteration_capacity_bound)) then
                    capacity_regularization = self%iteration_capacity_bound(workspace%material_id) / &
                                              real(max(1, nonlinear_iteration), real64)
                end if
            end if
        end if

        do i = 1, n_gauss
            call self%compute_iteration_capacity(workspace%material_id, workspace%state_gp(i), workspace%work_C(i))
            workspace%work_C(i) = workspace%work_C(i) + capacity_regularization
            if (coupling_mass_needed) then
                call self%compute_coupling_mass_term(workspace%material_id, workspace%state_gp(i), work_C_HT(i))
            end if
            call self%compute_transient_term_mixed(workspace%material_id, workspace%state_gp(i), &
                                                   workspace%bdf_coeffs, workspace%work_d_dt(i))


            if (thermal_target) then
                call self%calc_segregation_sink(workspace%material_id, workspace%state_gp(i), dt_local, work_sink(i))
            end if
        end do

        ! ----------------------------------------------------------------
        ! Element transport coefficients from the GEOMETRIC mean of the nodal
        ! values, not from quadrature of the coefficient itself.
        !
        ! K collapses by about eight orders of magnitude over the first kelvin
        ! of undercooling - k_r falls with the cryosuction head and the Hansson
        ! impedance 10^(-Omega Q) falls with the ice fraction - while the
        ! driving potential does the opposite: psi_cryo runs at 1.22e6 Pa/K, so
        ! 0.1 K across a 1 mm element is a gradient of 1.2e8 Pa/m. Sampling K at
        ! Gauss points and integrating gives the element the conductivity of its
        ! warmest sample, which is 1e3 to 1e8 times the cold-side value, so an
        ! element straddling the front conducts as if it were unfrozen. With the
        ! warm-side K that is a Darcy flux of order 10 mm/h through a 200 mm
        ! column, against an experiment that redistributes a few percent of its
        ! water over 50 h - and the resulting imbalance is what left the
        ! hydraulic residual at 0.3 in water content per step, unreachable by
        ! any step length or time step.
        !
        ! The geometric mean is the standard internodal weighting for a
        ! coefficient that varies exponentially: it is exact for exponential K
        ! between two points and cannot exceed the larger nodal value.
        ! ----------------------------------------------------------------
        ! Element transport coefficients.
        ! paths so a cut/uncut flip does not change the discrete operator.
        do i = 1, n_nodes
            coeff_sub_mat(:, :) = 0.0d0
            call self%compute_diffusion_term(workspace%material_id, workspace%state(i), coeff_sub_mat)
            D_HH_node(i) = coeff_sub_mat(1, 1)

            D_HH_gas_node(i) = 0.0d0
            if (gas_flux_needed) then
                coeff_gas_mat(:, :) = 0.0d0
                call self%compute_diffusion_term_gas(workspace%material_id, workspace%state(i), coeff_gas_mat)
                D_HH_gas_node(i) = coeff_gas_mat(1, 1)
            end if

            D_HT_node(i) = 0.0d0
            if (coupling_flux_needed) then
                coeff_sub_mat(:, :) = 0.0d0
                call self%compute_coupling_diffusion_term(workspace%material_id, &
                                                          workspace%state(i), coeff_sub_mat)
                D_HT_node(i) = coeff_sub_mat(1, 1)
            end if

            V_node_tmp(:) = 0.0d0
            call self%compute_advective_term(workspace%material_id, workspace%state(i), V_node_tmp)
            V_node_mag(i) = sqrt(sum(V_node_tmp(:)**2))
            if (i == 1) then
                V_dir(:) = 0.0d0
                if (V_node_mag(i) > 0.0d0) V_dir(:) = V_node_tmp(:) / V_node_mag(i)
            end if
        end do

        if (ELEMENT_CONDUCTIVITY_HARMONIC) then
            D_HH_elem = harmonic_mean(D_HH_node(1:n_nodes))
        else
            D_HH_elem = geometric_mean(D_HH_node(1:n_nodes))
        end if
        if (KELEM_REPORT) call report_element_conductivity(D_HH_node(1:n_nodes), D_HH_elem)
        D_HT_elem = signed_geometric_mean(D_HT_node(1:n_nodes))
        V_elem_mag = geometric_mean(V_node_mag(1:n_nodes))
        D_HH_gas_elem = 0.0d0
        if (gas_flux_needed) D_HH_gas_elem = geometric_mean(D_HH_gas_node(1:n_nodes))

        do i = 1, n_gauss
            do d = 1, n_dim
                workspace%work_D(d, d, i) = D_HH_elem
                work_D_HT(d, d, i) = D_HT_elem
            end do
            workspace%work_V(:, i) = V_elem_mag * V_dir(:)
        end do

        ! p_w is the total-water-equivalent pressure of the closure, so the
        ! liquid pressure is p_l = -s_eff. s_eff is read from the state the
        ! phase split published (unclamped s_m), not from the public accessor.
        do i = 1, n_nodes
            if (HYDRAULIC_DRIVER_TOTAL_POTENTIAL) then
                call workspace%state(i)%effective_suction%get(psi_eff_i)
                call workspace%state(i)%dSeff_dP%get(dSeff_dP_i)
                call workspace%state(i)%dSeff_dT%get(dSeff_dT_i)
                P_gen_node(i) = -psi_eff_i
                dPgen_dP_node(i) = -dSeff_dP_i
                dPgen_dT_node(i) = -dSeff_dT_i
            else
                P_gen_node(i) = workspace%P_node(i)
                dPgen_dP_node(i) = 1.0d0
                dPgen_dT_node(i) = 0.0d0
            end if
        end do

        ! ----------------------------------------------------------------
        ! 2. Mass Matrix K1 (LHS, factor bdf0).
        !    Uncut element: standard Gauss rule (C_eq/C_HT from step 1).
        !    q_s loop as the fluxes below (mat_C_sub/mat_C_HT_sub).
        ! ----------------------------------------------------------------
        call workspace%compute_K1(workspace%work_C, workspace%work_matrix)
        if (present(K_HH)) then
            do j = 1, n_nodes
                do i = 1, n_nodes
                    call K_HH%set(MATRIX_OPS%ADD, i, j, bdf0 * workspace%work_matrix(i, j))
                end do
            end do
        end if

        if (coupling_mass_needed) then
            call workspace%compute_K1(work_C_HT, workspace%work_matrix)
            do j = 1, n_nodes
                do i = 1, n_nodes
                    call K_HT%set(MATRIX_OPS%ADD, i, j, bdf0 * workspace%work_matrix(i, j))
                end do
            end do
        end if

        ! ----------------------------------------------------------------
        ! 4. Flux terms.
        !    Uncut element: standard Gauss rule (coefficients from step 1).
        !    gravity advection (replacement, so nothing is double-counted).
        ! ----------------------------------------------------------------
        if (coupling_flux_needed) then
            call workspace%compute_K2(work_D_HT, work_matrix_coupling)
            if (present(K_HT)) then
                do j = 1, n_nodes
                    do i = 1, n_nodes
                        call K_HT%set(MATRIX_OPS%ADD, i, j, work_matrix_coupling(i, j))
                    end do
                end do
            end if
            if (present(F_H)) then
                workspace%work_vec(:) = 0.0d0
                call matvec(work_matrix_coupling, workspace%T_node, workspace%work_vec, ierr)
                local_vec_res(:) = local_vec_res(:) + workspace%work_vec(:)
            end if
        end if

        call workspace%compute_K2(workspace%work_D, workspace%work_matrix)
        ! Column j of K2(D_HH) scaled by dP_gen/dp_w.
        if (present(K_HH)) then
            do j = 1, n_nodes
                do i = 1, n_nodes
                    call K_HH%set(MATRIX_OPS%ADD, i, j, workspace%work_matrix(i, j) * dPgen_dP_node(j))
                end do
            end do
        end if
        ! K_HT also takes K2(D_HH) scaled by dP_gen/dT, independent of the
        ! thermo-osmotic K2(D_HT) added above.
        if (coupling_block_needed) then
            do j = 1, n_nodes
                do i = 1, n_nodes
                    call K_HT%set(MATRIX_OPS%ADD, i, j, workspace%work_matrix(i, j) * dPgen_dT_node(j))
                end do
            end do
        end if
        if (present(F_H)) then
            do i = 1, n_nodes
                do j = 1, n_nodes
                    local_vec_res(i) = local_vec_res(i) + workspace%work_matrix(i, j) * P_gen_node(j)
                end do
            end do
        end if

        ! Vapour flux is driven by the gas pressure, so its column scale is
        ! unity and it contributes nothing to K_HT.
        if (gas_flux_needed) then
            work_D_gas(:, :, :) = 0.0d0
            do i = 1, n_gauss
                do d = 1, n_dim
                    work_D_gas(d, d, i) = D_HH_gas_elem
                end do
            end do
            call workspace%compute_K2(work_D_gas, work_matrix_gas)
            if (present(K_HH)) then
                do j = 1, n_nodes
                    do i = 1, n_nodes
                        call K_HH%set(MATRIX_OPS%ADD, i, j, work_matrix_gas(i, j))
                    end do
                end do
            end if
            if (present(F_H)) then
                workspace%work_vec(:) = 0.0d0
                call matvec(work_matrix_gas, workspace%P_node, workspace%work_vec, ierr)
                local_vec_res(:) = local_vec_res(:) + workspace%work_vec(:)
            end if
        end if

        ! ----------------------------------------------------------------
        ! 5. Residual Assembly
        ! ----------------------------------------------------------------
        if (present(F_H)) then
            workspace%work_vec(:) = 0.0d0
            call workspace%compute_R1(workspace%work_d_dt, workspace%work_vec)
            local_vec_res(:) = local_vec_res(:) + workspace%work_vec(:)

            workspace%work_vec(:) = 0.0d0
            call workspace%compute_R2(workspace%work_V, workspace%work_vec)
            local_vec_res(:) = local_vec_res(:) - workspace%work_vec(:)

            workspace%work_vec(:) = 0.0d0
            call workspace%compute_R1(work_sink, workspace%work_vec)
            local_vec_res(:) = local_vec_res(:) + workspace%work_vec(:)

            do i = 1, n_nodes
                call F_H%set(VECTOR_OPS%ADD, i, -local_vec_res(i))
            end do
        end if

    end subroutine assemble_local_picard_hydraulic

    !> Geometric mean of non-negative nodal values.
    !>
    !> Returns zero if any node is zero: a single impermeable node blocks the
    !> element, which is the physically correct series behaviour.
    !> Compare the element's geometric mean against the harmonic mean.
    !>
    !> Assumptions: called from inside the assembly's OpenMP region, hence the
    !> critical section; a no-op unless KELEM_REPORT.
    !> Numerical guarantees: none; diagnostic only, the operator is unchanged.
    !> Computational complexity: O(n_nodes).
    !> Failure behavior: silently skipped for a non-positive nodal value.
    subroutine report_element_conductivity(k_node, k_geometric)
        implicit none
        !> Nodal conductivities of one element
        real(real64), intent(in) :: k_node(:)
        !> Geometric mean the element actually uses
        real(real64), intent(in) :: k_geometric

        real(real64) :: k_lo, k_hi, k_harmonic

        if (size(k_node) == 0) return
        k_lo = minval(k_node)
        k_hi = maxval(k_node)
        if (k_lo <= 0.0d0) return
        if (k_hi < KELEM_REPORT_MIN_CONTRAST * k_lo) return

        k_harmonic = real(size(k_node), real64) / sum(1.0d0 / k_node)

        !$OMP CRITICAL (kelem_report_section)
        if (kelem_report_count < KELEM_REPORT_LIMIT) then
            kelem_report_count = kelem_report_count + 1
            write (*, '(A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3)') &
                '   [KELEM] k_min=', k_lo, ' k_max=', k_hi, &
                ' contrast=', k_hi / k_lo, &
                ' geo=', k_geometric, ' har=', k_harmonic, &
                ' geo/har=', k_geometric / k_harmonic
        end if
        !$OMP END CRITICAL (kelem_report_section)
    end subroutine report_element_conductivity

    !> Harmonic mean, the series-flow weight; zero if any value is non-positive.
    pure function harmonic_mean(values) result(mean_value)
        implicit none
        real(real64), intent(in) :: values(:)
        real(real64) :: mean_value

        integer(int32) :: i

        mean_value = 0.0d0
        if (size(values) == 0) return
        do i = 1, size(values)
            if (values(i) <= 0.0d0) return
        end do
        mean_value = real(size(values), real64) / sum(1.0d0 / values)
    end function harmonic_mean

    pure function geometric_mean(values) result(mean_value)
        implicit none
        real(real64), intent(in) :: values(:)
        real(real64) :: mean_value

        integer(int32) :: i

        mean_value = 0.0d0
        if (size(values) == 0) return
        do i = 1, size(values)
            if (values(i) <= 0.0d0) return
        end do
        mean_value = exp(sum(log(values)) / real(size(values), real64))
    end function geometric_mean

    !> Geometric mean of values that share a sign; zero on a sign change.
    !>
    !> The thermo-osmotic coefficient is negative wherever the total-potential
    !> head is negative, so it is averaged on its magnitude and the common sign
    !> is restored. A sign change inside one element has no geometric mean and
    !> is treated as no coupling there rather than as an arbitrary branch.
    pure function signed_geometric_mean(values) result(mean_value)
        implicit none
        real(real64), intent(in) :: values(:)
        real(real64) :: mean_value

        integer(int32) :: i
        logical :: all_positive, all_negative

        mean_value = 0.0d0
        if (size(values) == 0) return
        all_positive = all(values > 0.0d0)
        all_negative = all(values < 0.0d0)
        if (.not. (all_positive .or. all_negative)) return
        mean_value = exp(sum(log(abs(values))) / real(size(values), real64))
        if (all_negative) mean_value = -mean_value
    end function signed_geometric_mean

end submodule hydraulic_matrix
