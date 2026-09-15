!> Element-level finite-difference tangent of the assembled residual.
!>
!> Implementation strategy: rebuild one element's workspace exactly as the
!> production assembly does, then perturb one nodal primary variable at a time
!> and re-evaluate the element residual only. The element residual depends on
!> its own nodes alone - the Gauss-point gradients come from `fe%dlerp` of the
!> element's nodal values, and the nodal projected gradients enter only the
!> nodal fluxes, which the residual never reads - so the difference quotient is
!> the exact element tangent up to truncation, including every constitutive
!> derivative the hand-written blocks omit.
submodule(app_ftcms) ftcms_fd_jacobian
    use, intrinsic :: ieee_arithmetic, only:ieee_is_finite
    implicit none

    !> Relative perturbation. eps^(1/3) balances truncation against round-off
    !> for a central difference in double precision.
    real(real64), parameter :: FD_RELATIVE = 6.0d-6
    !> Difference the production tangent centrally. FD_RELATIVE above is the
    !> step a central difference wants; differencing forward at it leaves an
    !> O(delta) R'' term, and R'' is large at the front. Costs one extra
    !> residual evaluation per perturbed node.
    logical, parameter :: FD_CENTRAL_TANGENT = .true.
    !> Absolute perturbation floors, so a primary variable passing through zero
    !> still receives a meaningful increment: 1 K and 1 kPa.
    real(real64), parameter :: FD_FLOOR_TEMPERATURE = 1.0d0
    real(real64), parameter :: FD_FLOOR_PRESSURE = 1.0d3
    !> An element counts as a freezing-front element when any of its nodes lies
    !> within this band of 0 degC; the remainder form the control group.
    real(real64), parameter :: FD_FRONT_BAND_CELSIUS = 2.0d0
    !> Control elements are subsampled: the front band is what the measurement
    !> is about, and the far field only has to confirm the harness is sane.
    integer(int32), parameter :: FD_CONTROL_STRIDE = 97

contains

    !> Replace the thermal-pressure coupling block by its finite-difference
    !> tangent.
    !>
    !> The hand-written block carries only the storage term
    !> \( bdf_0 K_1(\partial U/\partial p) \). Measured against the element
    !> finite-difference tangent it is 18 times too small in the far field and
    !> 1798 times too small at the freezing front, because the energy residual
    !> also depends on pressure through the advective heat flux: the Darcy
    !> velocity carries \( -K_{wP}\nabla p \), so
    !> \( \partial [K_3(V)T]_i/\partial p_b
    !>    = -\int \psi_i \rho_w (c_w K_{wP} + c_v K_{vP})
    !>      (\nabla\psi_b\cdot\nabla T)\,d\Omega \).
    !> Differencing rather than adding that term by hand also captures the
    !> pressure dependence of \(K_{wP}\) through \(k_r(h)\) and of the densities.
    module subroutine assemble_coupling_fd_ftcms(self, workspace, local_F_T, K_TH)
        implicit none
        class(type_ftcms), intent(inout) :: self
        type(type_assemble_workspace), intent(inout) :: workspace
        !> Element thermal residual on entry; restored unchanged on exit
        type(type_vector_dp), intent(inout) :: local_F_T
        !> Overwritten with d R_T / d p
        type(type_matrix_dense), intent(inout) :: K_TH

        real(real64) :: caller_residual(workspace%num_fe_nodes)
        real(real64) :: base_residual(workspace%num_fe_nodes)
        real(real64), pointer :: values(:, :), residual(:)
        integer(int32) :: b, i, n_nodes
        real(real64) :: base_pressure, delta

        n_nodes = workspace%num_fe_nodes

        ! Keep the residual the caller is about to scatter.
        nullify (residual)
        residual => local_F_T%get_data()
        caller_residual(1:n_nodes) = residual(1:n_nodes)
        nullify (residual)

        base_residual(1:n_nodes) = caller_residual(1:n_nodes)

        nullify (values)
        values => K_TH%get_val()

        do b = 1, n_nodes
            call workspace%state(b)%pressure%get(base_pressure)
            delta = FD_RELATIVE * max(abs(base_pressure), FD_FLOOR_PRESSURE)

            call set_and_refresh(self, workspace, b, PHYSICS_TYPES%HYDRAULIC%ID, base_pressure + delta)

            call local_F_T%zero()
            call self%thermal%assemble_local(control=self%control, workspace=workspace, F_T=local_F_T)

            residual => local_F_T%get_data()
            ! The assembled vector is F = -R, so the residual tangent is -dF/dp.
            do i = 1, n_nodes
                values(i, b) = -(residual(i) - base_residual(i)) / delta
            end do
            nullify (residual)

            call workspace%state(b)%pressure%set(base_pressure)
        end do
        nullify (values)

        ! Put the workspace back on the unperturbed state and hand the caller
        ! back the residual it is about to scatter.
        call self%update_physical_properties_bulk(workspace%material_id, workspace%state)
        call refresh_gauss_states(self, workspace)
        residual => local_F_T%get_data()
        residual(1:n_nodes) = caller_residual(1:n_nodes)
        nullify (residual)
    end subroutine assemble_coupling_fd_ftcms

    !> Replace K_TT, K_TH, K_HT and K_HH together by finite-differencing the
    !> element thermal and hydraulic residuals against every nodal (T, p).
    !>
    !> Extends assemble_coupling_fd_ftcms's argument (see its doc comment) from
    !> the coupling block to all four: once the residual is C^1, the true
    !> element tangent is the difference quotient of the residual regardless of
    !> which block is being read off, and the remaining hand-written blocks are
    !> Modified-Picard blocks (transport coefficients frozen at the current
    !> iterate) rather than that tangent. One perturbation of a nodal variable
    !> changes both residuals at once, so both output columns are read from
    !> each perturbed evaluation: n_nodes temperature perturbations give
    !> K_TT(:,b) and K_HT(:,b) together, and n_nodes pressure perturbations
    !> give K_TH(:,b) and K_HH(:,b) together.
    !>
    !> Each perturbation goes through set_and_refresh (full nodal update, then
    !> Gauss-point update):
    !> hydraulic_matrix.F90's assemble_local_picard_hydraulic reads
    !> workspace%state(i) directly for the diffusion, coupling-diffusion and
    !> advective coefficients and for effective_suction and its tangents, so F_H,
    !> unlike F_T, depends on the perturbed node's own NODAL derived state, not
    !> only on the Gauss-point interpolation of it. Restoring a node's value in
    !> the loop uses a raw %set(); the next perturbation refreshes all nodal
    !> derived state before another residual is evaluated.
    module subroutine assemble_tangent_fd_ftcms(self, workspace, local_F_T, local_F_H, &
                                                K_TT, K_TH, K_HT, K_HH, ok)
        implicit none
        class(type_ftcms), intent(inout) :: self
        type(type_assemble_workspace), intent(inout) :: workspace
        !> Element thermal residual on entry; restored unchanged on exit
        type(type_vector_dp), intent(inout) :: local_F_T
        !> Element hydraulic residual on entry; restored unchanged on exit
        type(type_vector_dp), intent(inout) :: local_F_H
        !> Overwritten with d R_T / d T if ok
        type(type_matrix_dense), intent(inout) :: K_TT
        !> Overwritten with d R_T / d p if ok
        type(type_matrix_dense), intent(inout) :: K_TH
        !> Overwritten with d R_H / d T if ok
        type(type_matrix_dense), intent(inout) :: K_HT
        !> Overwritten with d R_H / d p if ok
        type(type_matrix_dense), intent(inout) :: K_HH
        !> .true. if every column below was finite and the blocks were
        !> overwritten; .false. leaves K_TT/K_TH/K_HT/K_HH untouched
        logical, intent(inout) :: ok

        real(real64) :: caller_residual_T(workspace%num_fe_nodes), caller_residual_H(workspace%num_fe_nodes)
        real(real64) :: base_residual_T(workspace%num_fe_nodes), base_residual_H(workspace%num_fe_nodes)
        ! Staged columns: not written into K_TT/K_TH/K_HT/K_HH until every
        ! column has proven finite, so a failed element keeps its analytic
        ! blocks (the sanity-guarded fallback selected in ftcms_assemble.F90).
        ! stage(:, b, row_physics, col_physics): perturbing physics c fills
        ! column c of every row's block.
        real(real64) :: stage(workspace%num_fe_nodes, workspace%num_fe_nodes, &
                              PHYSICS_TYPES%NUM_ID, PHYSICS_TYPES%NUM_ID)
        ! For a forward difference "minus" holds the base evaluation.
        real(real64) :: plus(workspace%num_fe_nodes, PHYSICS_TYPES%NUM_ID)
        real(real64) :: minus(workspace%num_fe_nodes, PHYSICS_TYPES%NUM_ID)
        real(real64), pointer :: residual_T(:), residual_H(:), values(:, :)
        integer(int32) :: b, i, n_nodes, c_phys, r_phys
        integer(int32), parameter :: ID_THERMAL = PHYSICS_TYPES%THERMAL%ID
        integer(int32), parameter :: ID_HYDRAULIC = PHYSICS_TYPES%HYDRAULIC%ID
        !> Physics perturbed by this routine, in column order
        integer(int32), parameter :: FD_COLUMNS(2) = [ID_THERMAL, ID_HYDRAULIC]
        real(real64) :: base_value, delta, denominator

        n_nodes = workspace%num_fe_nodes

        ! Keep the residuals the caller is about to scatter.
        nullify (residual_T, residual_H)
        residual_T => local_F_T%get_data()
        residual_H => local_F_H%get_data()
        caller_residual_T(1:n_nodes) = residual_T(1:n_nodes)
        caller_residual_H(1:n_nodes) = residual_H(1:n_nodes)
        nullify (residual_T, residual_H)

        base_residual_T(1:n_nodes) = caller_residual_T(1:n_nodes)
        base_residual_H(1:n_nodes) = caller_residual_H(1:n_nodes)

        ok = .true.

        ! Perturbing physics c fills column c of every row block, so one pair of
        ! residual evaluations serves every row.
        do c_phys = 1, size(FD_COLUMNS)
            do b = 1, n_nodes
                call get_state_value(workspace%state(b), FD_COLUMNS(c_phys), base_value)
                delta = FD_RELATIVE * max(abs(base_value), fd_floor(FD_COLUMNS(c_phys)))

                call set_and_refresh(self, workspace, b, FD_COLUMNS(c_phys), base_value + delta)

                call local_F_T%zero()
                call local_F_H%zero()
                call self%thermal%assemble_local(control=self%control, workspace=workspace, F_T=local_F_T)
                call self%hydraulic%assemble_local(control=self%control, workspace=workspace, F_H=local_F_H)

                residual_T => local_F_T%get_data()
                residual_H => local_F_H%get_data()
                plus(1:n_nodes, ID_THERMAL) = residual_T(1:n_nodes)
                plus(1:n_nodes, ID_HYDRAULIC) = residual_H(1:n_nodes)
                nullify (residual_T, residual_H)

                if (FD_CENTRAL_TANGENT) then
                    call set_and_refresh(self, workspace, b, FD_COLUMNS(c_phys), base_value - delta)
                    call local_F_T%zero()
                    call local_F_H%zero()
                    call self%thermal%assemble_local(control=self%control, workspace=workspace, F_T=local_F_T)
                    call self%hydraulic%assemble_local(control=self%control, workspace=workspace, F_H=local_F_H)
                    residual_T => local_F_T%get_data()
                    residual_H => local_F_H%get_data()
                    minus(1:n_nodes, ID_THERMAL) = residual_T(1:n_nodes)
                    minus(1:n_nodes, ID_HYDRAULIC) = residual_H(1:n_nodes)
                    nullify (residual_T, residual_H)
                    denominator = 2.0d0 * delta
                else
                    minus(1:n_nodes, ID_THERMAL) = base_residual_T(1:n_nodes)
                    minus(1:n_nodes, ID_HYDRAULIC) = base_residual_H(1:n_nodes)
                    denominator = delta
                end if

                ! The assembled vector is F = -R, so the residual tangent is -dF/du.
                do r_phys = 1, size(FD_COLUMNS)
                    do i = 1, n_nodes
                        stage(i, b, FD_COLUMNS(r_phys), FD_COLUMNS(c_phys)) = &
                            -(plus(i, FD_COLUMNS(r_phys)) - minus(i, FD_COLUMNS(r_phys))) / denominator
                    end do
                    if (any(.not. ieee_is_finite(stage(1:n_nodes, b, FD_COLUMNS(r_phys), FD_COLUMNS(c_phys))))) ok = .false.
                end do

                call set_state_value(workspace%state(b), FD_COLUMNS(c_phys), base_value)
            end do
        end do

        ! Put the workspace back on the unperturbed state and hand the caller
        ! back both residuals it is about to scatter, on
        ! every exit path (whether or not the staged blocks are accepted).
        !
        ! Every restore above used a raw %set() with no refresh, deferring the
        ! (expensive) full nodal update to whichever set_and_refresh call came
        ! next; after the last one there is none, so one explicit full nodal
        ! refresh is needed here to leave every node's derived state - not only
        ! the Gauss-point one - consistent with the now fully-restored raw
        ! values, matching this routine's "workspace left on the unperturbed
        ! state" contract.
        call self%update_physical_properties_bulk(workspace%material_id, workspace%state)
        call refresh_gauss_states(self, workspace)
        residual_T => local_F_T%get_data()
        residual_H => local_F_H%get_data()
        residual_T(1:n_nodes) = caller_residual_T(1:n_nodes)
        residual_H(1:n_nodes) = caller_residual_H(1:n_nodes)
        nullify (residual_T, residual_H)

        if (.not. ok) return

        nullify (values)
        values => K_TT%get_val()
        values(1:n_nodes, 1:n_nodes) = stage(1:n_nodes, 1:n_nodes, ID_THERMAL, ID_THERMAL)
        nullify (values)
        values => K_TH%get_val()
        values(1:n_nodes, 1:n_nodes) = stage(1:n_nodes, 1:n_nodes, ID_THERMAL, ID_HYDRAULIC)
        nullify (values)
        values => K_HT%get_val()
        values(1:n_nodes, 1:n_nodes) = stage(1:n_nodes, 1:n_nodes, ID_HYDRAULIC, ID_THERMAL)
        nullify (values)
        values => K_HH%get_val()
        values(1:n_nodes, 1:n_nodes) = stage(1:n_nodes, 1:n_nodes, ID_HYDRAULIC, ID_HYDRAULIC)
        nullify (values)
    end subroutine assemble_tangent_fd_ftcms

    !> Absolute perturbation floor for a primary variable passing through zero.
    pure function fd_floor(physics_id) result(floor_value)
        implicit none
        integer(int32), intent(in) :: physics_id
        real(real64) :: floor_value

        select case (physics_id)
        case (PHYSICS_TYPES%THERMAL%ID)
            floor_value = FD_FLOOR_TEMPERATURE
        case (PHYSICS_TYPES%HYDRAULIC%ID)
            floor_value = FD_FLOOR_PRESSURE
        case default
            floor_value = 1.0d0
        end select
    end function fd_floor

    subroutine get_state_value(state, physics_id, value)
        implicit none
        type(type_state), intent(in) :: state
        integer(int32), intent(in) :: physics_id
        real(real64), intent(inout) :: value

        select case (physics_id)
        case (PHYSICS_TYPES%THERMAL%ID)
            call state%temperature%get(value)
        case (PHYSICS_TYPES%HYDRAULIC%ID)
            call state%pressure%get(value)
        case default
            error stop 'get_state_value: no primary variable for this physics id'
        end select
    end subroutine get_state_value

    subroutine set_state_value(state, physics_id, value)
        implicit none
        type(type_state), intent(inout) :: state
        integer(int32), intent(in) :: physics_id
        real(real64), intent(in) :: value

        select case (physics_id)
        case (PHYSICS_TYPES%THERMAL%ID)
            call state%temperature%set(value)
        case (PHYSICS_TYPES%HYDRAULIC%ID)
            call state%pressure%set(value)
        case default
            error stop 'set_state_value: no primary variable for this physics id'
        end select
    end subroutine set_state_value

    !> Re-derive what the element thermal residual reads after a nodal primary
    !> variable changed.
    !>
    !> The thermal residual reads Gauss-point states and the nodal T/P buffers
    !> only; the nodal constitutive update and the dQi_dT projection feed matrix
    !> blocks alone, so both are skipped here.
    subroutine refresh_gauss_states(self, workspace)
        implicit none
        class(type_ftcms), intent(inout) :: self
        type(type_assemble_workspace), intent(inout) :: workspace

        call workspace%lerp()
        call self%update_physical_properties_bulk(workspace%material_id, workspace%state_gp)
    end subroutine refresh_gauss_states

    !> Report how far each assembled block is from the finite-difference tangent.
    module subroutine report_fd_jacobian_ftcms(self, tag, mean_rel_error, element_stride)
        implicit none
        class(type_ftcms), intent(inout) :: self
        !> Short label identifying the calling context in the log
        character(len=*), intent(in) :: tag
        !> Per-block mean relative error, in the block order used by the report;
        !> filled from the freezing-front group, or the control group when no
        !> front element exists. Unwritten entries stay negative.
        real(real64), intent(inout), optional :: mean_rel_error(:)
        !> Sample every n-th freezing-front element (default: all of them)
        integer(int32), intent(in), optional :: element_stride

        type(type_matrix_dense) :: local_K_TT, local_K_TH, local_K_HH, local_K_HT
        type(type_vector_dp) :: local_F_T, local_F_H
        type(type_assemble_workspace) :: workspace
        real(real64), allocatable :: elem_coords(:, :), raw_elem_coords(:, :)
        integer(int32), pointer, contiguous, dimension(:) :: connectivity

        real(real64), allocatable :: fd_TT(:, :), fd_TH(:, :), fd_HH(:, :), fd_HT(:, :)
        real(real64), allocatable :: an_TT(:, :), an_TH(:, :), an_HH(:, :), an_HT(:, :)
        real(real64), allocatable :: res_plus_T(:), res_plus_H(:)
        real(real64), allocatable :: res_minus_T(:), res_minus_H(:)

        integer(int32) :: num_elements, elem_id, n_nodes, b, group, front_stride, reported_group
        logical :: is_front, do_thermal, do_hydraulic
        real(real64) :: t_min, t_max, qi_max

        ! group 1 = freezing front, group 2 = control
        ! Slot 5 is the production coupling block: it re-runs assemble_coupling_fd
        ! and compares it against the full-refresh difference computed here, so
        ! the audit also proves that the lean refresh used in production leaves
        ! the thermal residual unchanged.
        integer(int32), parameter :: NUM_SLOTS = 5
        integer(int32) :: count_group(2)
        real(real64) :: sum_rel(2, NUM_SLOTS), max_rel(2, NUM_SLOTS), sum_ratio(2, NUM_SLOTS)
        integer(int32) :: worst_elem(NUM_SLOTS)
        real(real64) :: worst_rel(NUM_SLOTS), worst_ratio(NUM_SLOTS)
        real(real64) :: worst_T(NUM_SLOTS), worst_Qi(NUM_SLOTS)
        character(len=5), parameter :: BLOCK_NAME(NUM_SLOTS) = &
            ['TT   ', 'TH   ', 'HT   ', 'HH   ', 'TH_fd']
        real(real64), allocatable :: prod_TH(:, :)

        do_thermal = self%is_active_thermal()
        do_hydraulic = self%is_active_hydraulic()
        if (present(mean_rel_error)) mean_rel_error(:) = -1.0d0

        ! The re-assembled "analytic" K_HH block below must be regularization-
        ! free to be comparable against the FD residual tangent: the pseudo-
        ! transient capacity_regularization (hydraulic_matrix.F90) is added to
        ! the matrix only, never the residual, so leaving it enabled would make
        ! the audit compare two different operators. Restored .false. on every
        ! exit path, including this early return.
        call self%hydraulic%set_disable_capacity_regularization(.true.)
        if (.not. (do_thermal .and. do_hydraulic)) then
            call self%hydraulic%set_disable_capacity_regularization(.false.)
            return
        end if

        front_stride = 1
        if (present(element_stride)) front_stride = max(1, element_stride)

        nullify (connectivity)
        count_group(:) = 0
        sum_rel(:, :) = 0.0d0
        max_rel(:, :) = 0.0d0
        sum_ratio(:, :) = 0.0d0
        worst_elem(:) = 0
        worst_rel(:) = -1.0d0
        worst_ratio(:) = 0.0d0
        worst_T(:) = 0.0d0
        worst_Qi(:) = 0.0d0

        call self%domain%get_num_fe(num_elements)

        do elem_id = 1, num_elements
            call self%assemble_initialize(element_id=elem_id, workspace=workspace, &
                                          local_K_TT=local_K_TT, local_K_TH=local_K_TH, &
                                          local_K_HH=local_K_HH, local_K_HT=local_K_HT, &
                                          local_F_T=local_F_T, local_F_H=local_F_H, &
                                          coordinates=elem_coords, raw_coordinates=raw_elem_coords, &
                                          connectivity=connectivity)

            if (.not. self%control%is_target(PHYSICS_TYPES%THERMAL, workspace%material_id)) cycle
            if (.not. self%control%is_target(PHYSICS_TYPES%HYDRAULIC, workspace%material_id)) cycle

            n_nodes = workspace%num_fe_nodes
            call element_state_summary(workspace, n_nodes, t_min, t_max, qi_max)
            is_front = (t_min <= FD_FRONT_BAND_CELSIUS) .and. (t_max >= -FD_FRONT_BAND_CELSIUS)
            if (.not. is_front) then
                if (mod(elem_id, FD_CONTROL_STRIDE) /= 0) cycle
                group = 2
            else
                if (mod(elem_id, front_stride) /= 0) cycle
                group = 1
            end if

            call ensure_block(an_TT, n_nodes)
            call ensure_block(an_TH, n_nodes)
            call ensure_block(an_HH, n_nodes)
            call ensure_block(an_HT, n_nodes)
            call ensure_block(fd_TT, n_nodes)
            call ensure_block(fd_TH, n_nodes)
            call ensure_block(fd_HH, n_nodes)
            call ensure_block(fd_HT, n_nodes)

            ! Assembled (Picard) blocks at the unperturbed element state.
            call zero_blocks(local_K_TT, local_K_TH, local_K_HH, local_K_HT, local_F_T, local_F_H)
            call self%assemble_local(workspace, local_K_TT, local_K_TH, local_K_HH, local_K_HT, &
                                     local_F_T, local_F_H)
            call copy_block(local_K_TT, an_TT, n_nodes)
            call copy_block(local_K_TH, an_TH, n_nodes)
            call copy_block(local_K_HH, an_HH, n_nodes)
            call copy_block(local_K_HT, an_HT, n_nodes)

            ! Production coupling block, built exactly as the assembly does.
            call ensure_block(prod_TH, n_nodes)
            call self%assemble_coupling_fd(workspace, local_F_T, local_K_TH)
            call copy_block(local_K_TH, prod_TH, n_nodes)

            ! Finite-difference columns.
            do b = 1, n_nodes
                call fd_column(self, workspace, local_F_T, local_F_H, b, PHYSICS_TYPES%THERMAL%ID, &
                               res_plus_T, res_plus_H, res_minus_T, res_minus_H, n_nodes, &
                               fd_TT(:, b), fd_HT(:, b))
                call fd_column(self, workspace, local_F_T, local_F_H, b, PHYSICS_TYPES%HYDRAULIC%ID, &
                               res_plus_T, res_plus_H, res_minus_T, res_minus_H, n_nodes, &
                               fd_TH(:, b), fd_HH(:, b))
            end do

            count_group(group) = count_group(group) + 1
            call accumulate(1, an_TT, fd_TT)
            call accumulate(2, an_TH, fd_TH)
            call accumulate(3, an_HT, fd_HT)
            call accumulate(4, an_HH, fd_HH)
            call accumulate(5, prod_TH, fd_TH)
        end do

        call self%assemble_destroy(workspace, local_K_TT, local_K_TH, local_K_HH, local_K_HT, &
                                   local_F_T, local_F_H)
        if (allocated(elem_coords)) deallocate (elem_coords)
        if (allocated(raw_elem_coords)) deallocate (raw_elem_coords)

        call report()

        if (present(mean_rel_error)) then
            reported_group = 0
            if (count_group(1) > 0) then
                reported_group = 1
            else if (count_group(2) > 0) then
                reported_group = 2
            end if
            if (reported_group > 0) then
                do b = 1, min(size(mean_rel_error), NUM_SLOTS)
                    mean_rel_error(b) = sum_rel(reported_group, b) / real(count_group(reported_group), real64)
                end do
            end if
        end if

        call self%hydraulic%set_disable_capacity_regularization(.false.)

    contains

        !> Accumulate the block comparison statistics for one element.
        subroutine accumulate(slot, analytic, finite_difference)
            implicit none
            integer(int32), intent(in) :: slot
            real(real64), intent(in) :: analytic(:, :), finite_difference(:, :)

            real(real64) :: norm_fd, norm_an, rel

            norm_fd = sqrt(sum(finite_difference(1:n_nodes, 1:n_nodes)**2))
            norm_an = sqrt(sum(analytic(1:n_nodes, 1:n_nodes)**2))
            if (norm_fd <= 0.0d0 .and. norm_an <= 0.0d0) return

            rel = sqrt(sum((finite_difference(1:n_nodes, 1:n_nodes) - &
                            analytic(1:n_nodes, 1:n_nodes))**2)) / max(norm_fd, tiny(1.0d0))
            sum_rel(group, slot) = sum_rel(group, slot) + rel
            max_rel(group, slot) = max(max_rel(group, slot), rel)
            sum_ratio(group, slot) = sum_ratio(group, slot) + norm_fd / max(norm_an, tiny(1.0d0))

            if (group == 1 .and. rel > worst_rel(slot)) then
                worst_rel(slot) = rel
                worst_elem(slot) = elem_id
                worst_ratio(slot) = norm_fd / max(norm_an, tiny(1.0d0))
                worst_T(slot) = t_min
                worst_Qi(slot) = qi_max
            end if
        end subroutine accumulate

        !> Print the summary table.
        subroutine report()
            implicit none
            integer(int32) :: g, s
            character(len=7), parameter :: GROUP_NAME(2) = ['front  ', 'control']

            write (*, '(A,A,A,I0,A,I0)') '   [FD-JAC] ', trim(tag), &
                ' elements: front=', count_group(1), ' control=', count_group(2)
            write (*, '(A)') '   [FD-JAC] group   block  mean_rel_err  max_rel_err  mean ||J_fd||/||K||'
            do g = 1, 2
                if (count_group(g) == 0) cycle
                do s = 1, NUM_SLOTS
                    write (*, '(A,A,A,A,3(2X,ES11.3))') '   [FD-JAC] ', GROUP_NAME(g), '   ', BLOCK_NAME(s), &
                        sum_rel(g, s) / real(count_group(g), real64), max_rel(g, s), &
                        sum_ratio(g, s) / real(count_group(g), real64)
                end do
            end do
            do s = 1, NUM_SLOTS
                if (worst_elem(s) == 0) cycle
                write (*, '(A,A,A,I0,A,ES11.3,A,ES11.3,A,F8.3,A,F7.4)') &
                    '   [FD-JAC] worst ', BLOCK_NAME(s), ' elem=', worst_elem(s), &
                    ' rel_err=', worst_rel(s), ' ||J_fd||/||K||=', worst_ratio(s), &
                    ' T_min=', worst_T(s), ' Qi_max=', worst_Qi(s)
            end do
        end subroutine report

    end subroutine report_fd_jacobian_ftcms

    !> One finite-difference column: perturb nodal variable `physics_id` at local
    !> node `b`, re-evaluate the element residual, and return the two block
    !> columns.
    !>
    !> The assembled vector is F = -R, so the tangent of the residual is
    !> -dF/du; the sign is applied here so the result is directly comparable to
    !> the assembled K blocks.
    subroutine fd_column(self, workspace, local_F_T, local_F_H, b, physics_id, &
                         res_plus_T, res_plus_H, res_minus_T, res_minus_H, n_nodes, &
                         column_thermal, column_hydraulic)
        implicit none
        class(type_ftcms), intent(inout) :: self
        type(type_assemble_workspace), intent(inout) :: workspace
        type(type_vector_dp), intent(inout) :: local_F_T, local_F_H
        integer(int32), intent(in) :: b, physics_id, n_nodes
        real(real64), allocatable, intent(inout) :: res_plus_T(:), res_plus_H(:)
        real(real64), allocatable, intent(inout) :: res_minus_T(:), res_minus_H(:)
        real(real64), intent(inout) :: column_thermal(:), column_hydraulic(:)

        real(real64) :: base_value, delta
        integer(int32) :: i

        if (physics_id == PHYSICS_TYPES%THERMAL%ID) then
            call workspace%state(b)%temperature%get(base_value)
            delta = FD_RELATIVE * max(abs(base_value), FD_FLOOR_TEMPERATURE)
        else
            call workspace%state(b)%pressure%get(base_value)
            delta = FD_RELATIVE * max(abs(base_value), FD_FLOOR_PRESSURE)
        end if

        call set_and_refresh(self, workspace, b, physics_id, base_value + delta)
        call element_residual(self, workspace, local_F_T, local_F_H, res_plus_T, res_plus_H, n_nodes)

        call set_and_refresh(self, workspace, b, physics_id, base_value - delta)
        call element_residual(self, workspace, local_F_T, local_F_H, res_minus_T, res_minus_H, n_nodes)

        call set_and_refresh(self, workspace, b, physics_id, base_value)

        do i = 1, n_nodes
            column_thermal(i) = -(res_plus_T(i) - res_minus_T(i)) / (2.0d0 * delta)
            column_hydraulic(i) = -(res_plus_H(i) - res_minus_H(i)) / (2.0d0 * delta)
        end do
    end subroutine fd_column

    !> Set one nodal primary variable and rebuild every derived quantity the
    !> element residual reads, in the same order the production assembly does.
    subroutine set_and_refresh(self, workspace, b, physics_id, value)
        implicit none
        class(type_ftcms), intent(inout) :: self
        type(type_assemble_workspace), intent(inout) :: workspace
        integer(int32), intent(in) :: b, physics_id
        real(real64), intent(in) :: value

        call set_state_value(workspace%state(b), physics_id, value)

        call self%update_physical_properties_bulk(workspace%material_id, workspace%state)
        call workspace%lerp()
        call self%update_physical_properties_bulk(workspace%material_id, workspace%state_gp)
        ! Gate on the same switch the production assembly uses (ftcms_assemble.F90)
        ! so the audit compares the analytic Jacobian against the residual path
        ! actually used in production, whichever LERP_NODAL_DQI_DT selects.
        if (LERP_NODAL_DQI_DT) call workspace%lerp_dqi_dt()
        ! The perturbed node moves the freezing interface, and the element
        ! integration rule follows it. Dropping the cached split makes the
        ! differentiated residual use the same rule the residual itself would
        ! be assembled with at this state; see the invalidate routine's
        ! docstring for why the omitted term is first order, not small.
    end subroutine set_and_refresh

    !> Evaluate the element residual only, leaving the matrix blocks untouched.
    subroutine element_residual(self, workspace, local_F_T, local_F_H, residual_T, residual_H, n_nodes)
        implicit none
        class(type_ftcms), intent(inout) :: self
        type(type_assemble_workspace), intent(inout) :: workspace
        type(type_vector_dp), intent(inout) :: local_F_T, local_F_H
        real(real64), allocatable, intent(inout) :: residual_T(:), residual_H(:)
        integer(int32), intent(in) :: n_nodes

        real(real64), pointer :: data_T(:), data_H(:)

        call local_F_T%zero()
        call local_F_H%zero()
        call self%assemble_local(workspace, local_F_T=local_F_T, local_F_H=local_F_H)

        if (.not. allocated(residual_T)) allocate (residual_T(n_nodes))
        if (.not. allocated(residual_H)) allocate (residual_H(n_nodes))
        if (size(residual_T) /= n_nodes) then
            deallocate (residual_T)
            allocate (residual_T(n_nodes))
        end if
        if (size(residual_H) /= n_nodes) then
            deallocate (residual_H)
            allocate (residual_H(n_nodes))
        end if

        nullify (data_T, data_H)
        data_T => local_F_T%get_data()
        data_H => local_F_H%get_data()
        residual_T(1:n_nodes) = data_T(1:n_nodes)
        residual_H(1:n_nodes) = data_H(1:n_nodes)
        nullify (data_T, data_H)
    end subroutine element_residual

    !> Nodal temperature range and peak ice content, used to classify elements.
    subroutine element_state_summary(workspace, n_nodes, t_min, t_max, qi_max)
        implicit none
        type(type_assemble_workspace), intent(inout) :: workspace
        integer(int32), intent(in) :: n_nodes
        real(real64), intent(inout) :: t_min, t_max, qi_max

        integer(int32) :: i
        real(real64) :: temperature, ice

        t_min = huge(1.0d0)
        t_max = -huge(1.0d0)
        qi_max = 0.0d0
        do i = 1, n_nodes
            call workspace%state(i)%temperature%get(temperature)
            call workspace%state(i)%ice_content%get(ice)
            t_min = min(t_min, temperature)
            t_max = max(t_max, temperature)
            qi_max = max(qi_max, ice)
        end do
    end subroutine element_state_summary

    subroutine ensure_block(block, n_nodes)
        implicit none
        real(real64), allocatable, intent(inout) :: block(:, :)
        integer(int32), intent(in) :: n_nodes

        if (allocated(block)) then
            if (size(block, 1) == n_nodes .and. size(block, 2) == n_nodes) return
            deallocate (block)
        end if
        allocate (block(n_nodes, n_nodes), source=0.0d0)
    end subroutine ensure_block

    subroutine copy_block(matrix, destination, n_nodes)
        implicit none
        type(type_matrix_dense), intent(inout) :: matrix
        real(real64), intent(inout) :: destination(:, :)
        integer(int32), intent(in) :: n_nodes

        real(real64), pointer :: values(:, :)

        nullify (values)
        values => matrix%get_val()
        destination(1:n_nodes, 1:n_nodes) = values(1:n_nodes, 1:n_nodes)
        nullify (values)
    end subroutine copy_block

    subroutine zero_blocks(K_TT, K_TH, K_HH, K_HT, F_T, F_H)
        implicit none
        type(type_matrix_dense), intent(inout) :: K_TT, K_TH, K_HH, K_HT
        type(type_vector_dp), intent(inout) :: F_T, F_H

        call K_TT%zero()
        call K_TH%zero()
        call K_HH%zero()
        call K_HT%zero()
        call F_T%zero()
        call F_H%zero()
    end subroutine zero_blocks

end submodule ftcms_fd_jacobian
