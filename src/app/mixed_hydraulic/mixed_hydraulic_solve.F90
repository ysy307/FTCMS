#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscsnes.h>

!> The backward-Euler step and the balance it is audited with.
!>
!> The step is discretised here rather than by TS because the conservation
!> residual differences a stored quantity between two known times; the accepted
!> history is an auxiliary field, and a rejected step must leave it untouched so
!> that the same state can be retried at a shorter length.
!>
!> The audit rebuilds every flux and divergence from the solved degrees of
!> freedom instead of reading the residual PETSc assembled. Both sides of an
!> internal facet are rebuilt separately, so a wrong orientation shows up as a
!> mismatch rather than cancelling against itself. Only owned cells are summed,
!> which is what leaves the physical boundary and source budget after the
!> internal contributions cancel.
submodule(app_mixed_hydraulic) mixed_hydraulic_solve
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use :: app_mixed_hydraulic_physics, only:A_PRESSURE_OLD, A_SOURCE, A_MATERIAL, C_DT
    use :: app_mixed_hydraulic_balance, only:MAX_CELL_FACETS, reconstruct_cell_facet_flux, &
                                             reconstruct_cell_divergence, reconstruct_cell_pressure, &
                                             cell_is_owned
    use :: app_mixed_hydraulic_closure, only:eval_storage, clear_mobility_fault, &
                                             mobility_fault_recorded, mobility_fault_report
    use :: app_mixed_hydraulic_gauge, only:type_mixed_gauge_report, measure_gauge_state, &
                                           install_pressure_gauge, gauge_data_is_compatible
    implicit none

contains

    module subroutine solve_step(self, time, dt, accepted, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        real(real64), intent(in) :: time
        real(real64), intent(in) :: dt
        logical, intent(inout) :: accepted
        PetscErrorCode, intent(inout) :: ierr

        SNESConvergedReason :: reason
        type(type_mixed_gauge_report) :: gauge_report
        Vec :: boundary_state
        PetscErrorCode :: restore_ierr
        PetscMPIInt :: rank

        accepted = .false.
        ierr = 0
        ! A failure to undo a failed trial must not hide what failed, so the
        ! restore reports separately and the original error is the one returned.
        restore_ierr = 0
        if (.not. self%ready) then
            ierr = 1
            return
        end if
        ! Nothing has been touched yet, so a refused step needs no rollback.
        ! A nonpositive or nonfinite length would divide the storage term by
        ! itself; it is reported rather than accepted with a zero storage rate.
        if (.not. step_is_usable(dt)) then
            ierr = MIXED_ERR_INVALID_STEP
            return
        end if

        call begin_trial(self, dt, ierr)
        if (ierr /= 0) then
            call restore_accepted_state(self, restore_ierr)
            return
        end if

        ! Test compatibility before solving. With a cell that can store water
        ! the pressure constant is physical, so no gauge applies.
        call self%get_local_solution(boundary_state, ierr)
        if (ierr /= 0) then
            call restore_accepted_state(self, restore_ierr)
            return
        end if
        call measure_gauge_state(self%dm, boundary_state, self%aux_dm, self%auxiliary, &
                                 self%flux_scale, gauge_report, ierr)
        if (ierr /= 0) then
            call restore_accepted_state(self, restore_ierr)
            return
        end if
        if (gauge_report%eligible) then
            if (.not. gauge_data_is_compatible(gauge_report)) then
                call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
                if (rank == 0) write (*, '(a)') "Incompatible steady all-flux data: boundary flux does not balance source"
                call restore_accepted_state(self, restore_ierr)
                ierr = MIXED_ERR_INCOMPATIBLE
                return
            end if
            if (.not. self%pressure_gauged) then
                call install_pressure_gauge(self%dm, self%length_scale, gauge_report, ierr)
                if (ierr /= 0) then
                    call restore_accepted_state(self, restore_ierr)
                    return
                end if
                self%pressure_gauged = .true.
            end if
        end if

        ! Whatever a kernel could not evaluate during this trial is this trial's
        ! own, so the record starts empty and is read only on a failed path.
        call clear_mobility_fault()

        call SNESSetSolution(self%snes, self%solution, ierr)
        if (ierr == 0) call SNESSolve(self%snes, PETSC_NULL_VEC, self%solution, ierr)
        if (ierr /= 0) then
            ! A PETSc error leaves the iterate wherever the failed solve
            ! stopped, which is not a state any history describes.
            call restore_accepted_state(self, restore_ierr)
            call report_mobility_fault(ierr)
            return
        end if

        call SNESGetConvergedReason(self%snes, reason, ierr)
        if (ierr /= 0) then
            call restore_accepted_state(self, restore_ierr)
            call report_mobility_fault(ierr)
            return
        end if
        if (int(reason%v, int32) <= 0) then
            call restore_accepted_state(self, ierr)
            call report_mobility_fault(ierr)
            return
        end if

        ! A converged reason is not a finite answer: a solve that overflowed can
        ! satisfy a relative test on a nonfinite norm.
        if (.not. solution_is_finite(self, ierr)) then
            call restore_accepted_state(self, restore_ierr)
            if (ierr == 0) ierr = MIXED_ERR_NONFINITE
            return
        end if

        call VecCopy(self%solution, self%old_solution, ierr)
        if (ierr /= 0) then
            call restore_accepted_state(self, restore_ierr)
            return
        end if
        self%accepted_dt = dt
        accepted = .true.
    end subroutine solve_step

    !> Turn a recorded mobility fault into this step's reason for failing.
    !>
    !> A pointwise kernel cannot refuse a step itself: it writes a nonfinite
    !> residual, which SNES turns into a diverged reason or a PETSc error
    !> depending on the run's settings. Either way the step is already being
    !> rejected by the time this is reached, and the only thing left is to say
    !> which state the model had no conductivity for. The fault takes precedence
    !> over the generic error because it is the cause of it.
    subroutine report_mobility_fault(ierr)
        implicit none
        PetscErrorCode, intent(inout) :: ierr

        integer(int32) :: slot, local_fault, any_fault, mpi_err
        real(real64) :: pressure, value

        ! The fault is per rank but the step is refused on all of them, so the
        ! flag is reduced: every rank has to come back with the same reason.
        local_fault = 0
        if (mobility_fault_recorded()) local_fault = 1
        call MPI_Allreduce(local_fault, any_fault, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, mpi_err)
        if (mpi_err /= 0 .or. any_fault == 0) return

        if (local_fault == 1) then
            slot = 0
            pressure = 0.0_real64
            value = 0.0_real64
            call mobility_fault_report(slot, pressure, value)
            write (*, '(a,i0,a,es13.6,a,es13.6)') &
                "Unusable mobility: material slot ", slot, " at pressure ", pressure, &
                " Pa has conductivity ", value
        end if
        ierr = MIXED_ERR_MOBILITY
    end subroutine report_mobility_fault

    !> Whether this step length can be differenced over at all.
    pure logical function step_is_usable(dt)
        implicit none
        real(real64), intent(in) :: dt

        ! Fortran does not guarantee short-circuit evaluation of .and.; a
        ! comparison with NaN can trap even after a failed finiteness test.
        step_is_usable = .false.
        if (.not. ieee_is_finite(dt)) return
        step_is_usable = (dt > 0.0_real64)
    end function step_is_usable

    !> Put the solver into the state one trial of length `dt` needs.
    !>
    !> Everything this writes is a function of the accepted state, which is why
    !> `restore_accepted_state` can undo it by writing it again.
    subroutine begin_trial(self, dt, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        real(real64), intent(in) :: dt
        PetscErrorCode, intent(inout) :: ierr

        ! Keep the auxiliary data as the accepted state left it, before the
        ! storage history below is written over it.
        call VecCopy(self%auxiliary, self%accepted_auxiliary, ierr)
        if (ierr /= 0) return

        self%dt = dt
        self%constants(C_DT) = real(dt, PETSC_REAL_KIND)
        call PetscDSSetConstants(self%ds, int(MIXED_CONSTANTS, PETSC_INT_KIND), self%constants, ierr)
        if (ierr /= 0) return

        ! The accepted state the storage term differences against. Written from
        ! old_solution, not from the iterate, so a retried step sees the same
        ! history the failed one did.
        call write_accepted_pressure(self, ierr)
    end subroutine begin_trial

    !> Undo a failed trial, exactly.
    !>
    !> The iterate comes back from `old_solution`, which a failed step never
    !> writes, and the auxiliary data from the copy `begin_trial` kept. The
    !> history is restored rather than rewritten from `old_solution`: after an
    !> accepted step the accepted history describes the start of that step, not
    !> the state the step reached, so rewriting it would silently advance the
    !> accepted state that a rejected step is supposed to leave alone -- and the
    !> balance audit of the last accepted step would then report the source term
    !> of a step that never happened. Each boundary's copy is refreshed from the
    !> restored volume data, which is where it came from.
    subroutine restore_accepted_state(self, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        integer(int32) :: b

        call VecCopy(self%old_solution, self%solution, ierr)
        if (ierr == 0) call VecCopy(self%accepted_auxiliary, self%auxiliary, ierr)
        if (ierr /= 0) return

        do b = 1, self%num_boundary_aux
            call copy_history_into(self, self%boundary_aux(b), ierr)
            if (ierr /= 0) return
        end do

        self%dt = self%accepted_dt
        self%constants(C_DT) = real(self%accepted_dt, PETSC_REAL_KIND)
        call PetscDSSetConstants(self%ds, int(MIXED_CONSTANTS, PETSC_INT_KIND), self%constants, ierr)
    end subroutine restore_accepted_state

    !> Whether every entry of the solved state is a finite number.
    logical function solution_is_finite(self, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        PetscReal :: measured

        solution_is_finite = .false.
        call VecNorm(self%solution, NORM_INFINITY, measured, ierr)
        if (ierr /= 0) return
        solution_is_finite = ieee_is_finite(real(measured, real64))
    end function solution_is_finite

    !> Copy the accepted cell pressures into the auxiliary history field.
    subroutine write_accepted_pressure(self, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        Vec :: local_old
        PetscSection :: solve_section, aux_section
        PetscScalar, pointer :: old_raw(:), aux_raw(:)
        PetscInt :: cell, c_start, c_end, solve_offset, aux_offset, num_dofs
        integer(int32) :: b

        call DMGetLocalVector(self%dm, local_old, ierr)
        if (ierr /= 0) return
        call DMGlobalToLocalBegin(self%dm, self%old_solution, INSERT_VALUES, local_old, ierr)
        if (ierr == 0) call DMGlobalToLocalEnd(self%dm, self%old_solution, INSERT_VALUES, local_old, ierr)
        if (ierr == 0) call DMGetLocalSection(self%dm, solve_section, ierr)
        if (ierr == 0) call DMGetLocalSection(self%aux_dm, aux_section, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(self%dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) then
            call DMRestoreLocalVector(self%dm, local_old, ierr)
            return
        end if

        call VecGetArrayRead(local_old, old_raw, ierr)
        if (ierr == 0) call VecGetArray(self%auxiliary, aux_raw, ierr)
        if (ierr == 0) then
            do cell = c_start, c_end - 1
                call PetscSectionGetFieldDof(aux_section, cell, int(A_PRESSURE_OLD - 1, PETSC_INT_KIND), num_dofs, ierr)
                if (ierr /= 0 .or. num_dofs <= 0) cycle
                call PetscSectionGetFieldOffset(aux_section, cell, int(A_PRESSURE_OLD - 1, PETSC_INT_KIND), &
                                                aux_offset, ierr)
                if (ierr /= 0) cycle
                call PetscSectionGetFieldDof(solve_section, cell, 1_PETSC_INT_KIND, num_dofs, ierr)
                if (ierr /= 0 .or. num_dofs <= 0) cycle
                call PetscSectionGetFieldOffset(solve_section, cell, 1_PETSC_INT_KIND, solve_offset, ierr)
                if (ierr /= 0) cycle
                aux_raw(aux_offset + 1) = real(self%pressure_scale, PETSC_REAL_KIND)*old_raw(solve_offset + 1)
            end do
            call VecRestoreArray(self%auxiliary, aux_raw, ierr)
        end if
        call VecRestoreArrayRead(local_old, old_raw, ierr)
        call DMRestoreLocalVector(self%dm, local_old, ierr)

        ! Every boundary's own copy carries the same history, because its
        ! kernels read the cell behind the facet.
        do b = 1, self%num_boundary_aux
            call copy_history_into(self, self%boundary_aux(b), ierr)
            if (ierr /= 0) return
        end do
    end subroutine write_accepted_pressure

    !> Copy the history field of the volume auxiliary vector into a boundary's.
    subroutine copy_history_into(self, target_vec, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        Vec, intent(inout) :: target_vec
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: section
        PetscScalar, pointer :: source_raw(:), target_raw(:)
        PetscInt :: point, p_start, p_end, offset, num_dofs

        call DMGetLocalSection(self%aux_dm, section, ierr)
        if (ierr == 0) call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr /= 0) return

        call VecGetArrayRead(self%auxiliary, source_raw, ierr)
        if (ierr /= 0) return
        call VecGetArray(target_vec, target_raw, ierr)
        if (ierr == 0) then
            do point = p_start, p_end - 1
                call PetscSectionGetFieldDof(section, point, int(A_PRESSURE_OLD - 1, PETSC_INT_KIND), num_dofs, ierr)
                if (ierr /= 0 .or. num_dofs <= 0) cycle
                call PetscSectionGetFieldOffset(section, point, int(A_PRESSURE_OLD - 1, PETSC_INT_KIND), offset, ierr)
                if (ierr /= 0) cycle
                target_raw(offset + 1) = source_raw(offset + 1)
            end do
            call VecRestoreArray(target_vec, target_raw, ierr)
        end if
        call VecRestoreArrayRead(self%auxiliary, source_raw, ierr)
    end subroutine copy_history_into

    module subroutine compute_balance(self, cell_defect, boundary_rate, domain_defect, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        real(real64), intent(inout) :: cell_defect(:)
        real(real64), intent(inout) :: boundary_rate
        real(real64), intent(inout) :: domain_defect
        PetscErrorCode, intent(inout) :: ierr

        type(type_mixed_audit) :: report

        cell_defect = 0.0_real64
        boundary_rate = 0.0_real64
        domain_defect = 0.0_real64

        call audit_cells(self, report, cell_defect, ierr)
        if (ierr /= 0) return
        call reduce_report(report, ierr)
        if (ierr /= 0) return

        boundary_rate = report%boundary_rate
        domain_defect = report%domain_defect
    end subroutine compute_balance

    module subroutine compute_audit(self, report, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        type(type_mixed_audit), intent(inout) :: report
        PetscErrorCode, intent(inout) :: ierr

        real(real64) :: unused(1)

        unused = 0.0_real64
        call audit_cells(self, report, unused, ierr)
        if (ierr /= 0) return
        call reduce_report(report, ierr)
    end subroutine compute_audit

    !> Walk the owned cells once, measuring every quantity of the audit.
    subroutine audit_cells(self, report, cell_defect, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        type(type_mixed_audit), intent(inout) :: report
        !> Filled for as many owned cells as it holds; the rest stay zero.
        real(real64), intent(inout) :: cell_defect(:)
        PetscErrorCode, intent(inout) :: ierr

        Vec :: local
        PetscSection :: aux_section
        PetscScalar, pointer :: aux_raw(:)
        PetscInt :: cell, c_start, c_end, offset, num_dofs
        PetscInt :: facets(MAX_CELL_FACETS), neighbour_facets(MAX_CELL_FACETS)
        PetscInt, pointer :: support(:)
        real(real64) :: moments(MAX_CELL_FACETS), neighbour_moments(MAX_CELL_FACETS)
        real(real64) :: divergence_integral, volume, outward_sum, defect
        real(real64) :: pressure, pressure_old, source
        real(real64) :: stored_now, stored_old, capacity
        integer(int32) :: num_facets, neighbour_facet_count, neighbour_count, f, g, owned_index
        integer(int32) :: slot
        PetscInt :: neighbour

        report%max_facet_mismatch = 0.0_real64
        report%max_divergence_defect = 0.0_real64
        report%max_cell_defect = 0.0_real64
        report%domain_defect = 0.0_real64
        report%boundary_rate = 0.0_real64
        report%source_rate = 0.0_real64
        report%throughput = 0.0_real64
        report%owned_cells = 0
        owned_index = 0

        call self%get_local_solution(local, ierr)
        if (ierr /= 0) return

        call DMGetLocalSection(self%aux_dm, aux_section, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(self%dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return

        call VecGetArrayRead(self%auxiliary, aux_raw, ierr)
        if (ierr /= 0) return

        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(self%dm, cell)) cycle
            report%owned_cells = report%owned_cells + 1

            call reconstruct_cell_facet_flux(self%dm, local, self%flux_scale, cell, facets, moments, num_facets, ierr)
            if (ierr /= 0) exit
            call reconstruct_cell_divergence(self%dm, local, self%flux_scale, cell, divergence_integral, volume, ierr)
            if (ierr /= 0) exit
            call reconstruct_cell_pressure(self%dm, local, self%pressure_scale, cell, pressure, ierr)
            if (ierr /= 0) exit

            outward_sum = sum(moments(1:num_facets))
            report%max_divergence_defect = max(report%max_divergence_defect, abs(outward_sum - divergence_integral))
            do f = 1, num_facets
                report%throughput = max(report%throughput, abs(moments(f)))
            end do

            ! Each facet, from the other side. A facet whose second cell is not
            ! local belongs to a rank that owns the pair, so it is left alone
            ! here; with an overlap of one every owned facet has both sides.
            do f = 1, num_facets
                call DMPlexGetSupport(self%dm, facets(f), support, ierr)
                if (ierr /= 0) exit
                neighbour_count = int(size(support), int32)
                if (neighbour_count == 1) then
                    report%boundary_rate = report%boundary_rate + moments(f)
                    cycle
                end if
                neighbour = support(1)
                if (neighbour == cell) neighbour = support(2)
                if (neighbour < c_start .or. neighbour >= c_end) cycle
                call reconstruct_cell_facet_flux(self%dm, local, self%flux_scale, neighbour, &
                                                 neighbour_facets, neighbour_moments, neighbour_facet_count, ierr)
                if (ierr /= 0) exit
                do g = 1, neighbour_facet_count
                    if (neighbour_facets(g) /= facets(f)) cycle
                    report%max_facet_mismatch = max(report%max_facet_mismatch, &
                                                    abs(moments(f) + neighbour_moments(g)))
                end do
            end do
            if (ierr /= 0) exit

            source = 0.0_real64
            pressure_old = 0.0_real64
            slot = 0
            call PetscSectionGetFieldDof(aux_section, cell, int(A_SOURCE - 1, PETSC_INT_KIND), num_dofs, ierr)
            if (ierr == 0 .and. num_dofs > 0) then
                call PetscSectionGetFieldOffset(aux_section, cell, int(A_SOURCE - 1, PETSC_INT_KIND), offset, ierr)
                if (ierr == 0) source = real(aux_raw(offset + 1), real64)
                call PetscSectionGetFieldOffset(aux_section, cell, int(A_PRESSURE_OLD - 1, PETSC_INT_KIND), offset, ierr)
                if (ierr == 0) pressure_old = real(aux_raw(offset + 1), real64)
                call PetscSectionGetFieldOffset(aux_section, cell, int(A_MATERIAL - 1, PETSC_INT_KIND), offset, ierr)
                if (ierr == 0) slot = nint(real(aux_raw(offset + 1), real64), int32)
            end if

            ! The same stored quantity the residual differences, at the same two
            ! pressures: the accepted one comes from the history field, which a
            ! step advances before it solves and never after, so this still
            ! reads the start of the step once the step has been accepted.
            stored_now = 0.0_real64
            stored_old = 0.0_real64
            capacity = 0.0_real64
            call eval_storage(slot, pressure, stored_now, capacity)
            call eval_storage(slot, pressure_old, stored_old, capacity)

            defect = (stored_now - stored_old)*volume &
                     + self%dt*(outward_sum - source*volume)
            report%source_rate = report%source_rate + source*volume
            report%domain_defect = report%domain_defect + defect
            report%max_cell_defect = max(report%max_cell_defect, abs(defect))

            owned_index = owned_index + 1
            if (owned_index <= size(cell_defect)) cell_defect(owned_index) = defect
        end do

        call VecRestoreArrayRead(self%auxiliary, aux_raw, ierr)
    end subroutine audit_cells

    !> Reduce the maxima and the sums over the communicator.
    subroutine reduce_report(report, ierr)
        implicit none
        type(type_mixed_audit), intent(inout) :: report
        PetscErrorCode, intent(inout) :: ierr

        real(real64) :: maxima_send(4), maxima(4)
        real(real64) :: sums_send(3), sums(3)
        integer(int32) :: mpi_err

        maxima_send = [report%max_facet_mismatch, report%max_divergence_defect, &
                       report%max_cell_defect, report%throughput]
        sums_send = [report%domain_defect, report%boundary_rate, report%source_rate]

        call MPI_Allreduce(maxima_send, maxima, 4, MPI_DOUBLE_PRECISION, MPI_MAX, PETSC_COMM_WORLD, mpi_err)
        if (mpi_err == 0) then
            call MPI_Allreduce(sums_send, sums, 3, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, mpi_err)
        end if
        if (mpi_err /= 0) then
            ierr = 1
            return
        end if

        report%max_facet_mismatch = maxima(1)
        report%max_divergence_defect = maxima(2)
        report%max_cell_defect = maxima(3)
        report%throughput = maxima(4)
        report%domain_defect = sums(1)
        report%boundary_rate = sums(2)
        report%source_rate = sums(3)
    end subroutine reduce_report

end submodule mixed_hydraulic_solve
