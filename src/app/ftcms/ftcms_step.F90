#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscsnes.h>

!> Advancing the engine in time.
submodule(app_ftcms) ftcms_step
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use :: stdlib_logger, only:global_logger
    use :: app_ftcms_kernels, only:AUX_TEMPERATURE_OLD, AUX_PRESSURE_OLD, AUX_TEMPERATURE_OLDER, &
                                   AUX_PRESSURE_OLDER, clear_kernel_fault, &
                                   kernel_fault_recorded, kernel_fault_report
    use :: app_ftcms_boundary_registry, only:set_boundary_time, boundary_data
    use :: module_boundary, only:BC_DATA_SIZE
    use :: app_ftcms_closure, only:closure_work_counters
    implicit none

contains

    !> Advance from the start to the end of the time window.
    !>
    !> Each step is attempted at the controller's length, clipped so that the
    !> window and every output time are met exactly. A step is refused when
    !> SNES does not converge, and, with error control, when its local error
    !> exceeds the tolerance; it is then retried shorter from the accepted
    !> state. A step refused at the floor stops the run with an error.
    module subroutine run_ftcms(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        real(real64) :: dt, error
        integer(int32) :: iterations
        logical :: converged

        if (.not. self%ready) error stop "The engine runs only after it is initialised"

        call assimilate(self)
        do while (self%time < self%end_time)
            dt = min(self%next_step, self%end_time - self%time, next_output_time(self) - self%time)
            dt = max(dt, min(self%step_control%minimum, self%end_time - self%time))
            call choose_order(self, dt)
            call attempt_step(self, dt, converged, iterations)
            error = -1.0d0
            if (converged .and. self%error_control%active) then
                error = local_error(self, dt)
                if (error > 1.0d0) then
                    converged = .false.
                    call refuse_step(self)
                end if
            end if

            if (converged) then
                call accept_step(self, dt)
                call assimilate(self)
                self%accepted_steps = self%accepted_steps + 1
                self%total_iterations = self%total_iterations + iterations
                self%next_step = next_length(self, max(self%next_step, dt), iterations, error)
                call write_output(self, force=.false.)
                call report_progress(self, iterations, .false.)
            else
                self%rejected_steps = self%rejected_steps + 1
                if (dt <= self%step_control%minimum) then
                    call report_progress(self, iterations, .true.)
                    error stop "A step was refused at the minimum step length of Conditions.json"
                end if
                self%next_step = retry_length(self, dt, error)
                self%error_control%after_rejection = .true.
            end if
        end do
        call write_output(self, force=.true.)
        call report_progress(self, 0, .true.)
        call report_closure_work()
    end subroutine run_ftcms

    !> The BDF order and weights of a step of length `dt`.
    !>
    !> Second order once one accepted step lies behind, unless the step grows
    !> by more than \(1+\sqrt2\) over the last one, beyond which variable-step
    !> BDF2 is not zero-stable; that step is taken at first order.
    subroutine choose_order(self, dt)
        implicit none
        class(type_ftcms), intent(inout) :: self
        real(real64), intent(in) :: dt

        real(real64), parameter :: MAX_STEP_RATIO = 1.0d0 + sqrt(2.0d0)
        real(real64) :: ratio

        self%order = 1
        self%bdf_weights = [1.0d0, -1.0d0, 0.0d0]
        if (self%target_order < 2 .or. self%history_depth < 1) return
        ratio = dt/self%history_steps(1)
        if (ratio > MAX_STEP_RATIO) return
        self%order = 2
        self%bdf_weights = [(1.0d0 + 2.0d0*ratio)/(1.0d0 + ratio), -(1.0d0 + ratio), &
                            ratio**2/(1.0d0 + ratio)]
    end subroutine choose_order

    !> The length after an accepted step: PI on the local error, or the iteration count.
    function next_length(self, dt, iterations, error) result(next)
        implicit none
        class(type_ftcms), intent(inout) :: self
        real(real64), intent(in) :: dt
        integer(int32), intent(in) :: iterations
        real(real64), intent(in) :: error
        real(real64) :: next

        real(real64) :: ratio, current

        if (.not. (error > 0.0d0)) then
            next = next_step_length(self%step_control, dt, iterations)
            return
        end if
        associate (control => self%error_control)
            current = max(error, 1.0d-12)
            ratio = control%safety*current**(-(control%integral_gain + control%proportional_gain)) &
                    *max(control%previous_error, 1.0d-12)**control%proportional_gain
            ratio = min(max(ratio, 1.0d0/control%max_growth), control%max_growth)
            ! Right after a refusal the error history describes a step never
            ! taken; growing on it re-enters the step that just failed.
            if (control%after_rejection) ratio = min(ratio, 1.0d0)
            control%after_rejection = .false.
            control%previous_error = current
        end associate
        next = dt*ratio
        if (iterations >= self%step_control%iterations_high) next = min(next, dt*self%step_control%shrink)
        next = min(max(next, self%step_control%minimum), self%step_control%maximum)
    end function next_length

    !> The length a refused step is retried at.
    function retry_length(self, dt, error) result(next)
        implicit none
        class(type_ftcms), intent(in) :: self
        real(real64), intent(in) :: dt
        real(real64), intent(in) :: error
        real(real64) :: next

        real(real64) :: ratio

        if (error > 1.0d0) then
            ratio = self%error_control%safety*error**(-1.0d0/real(self%order + 1, real64))
            next = dt*min(0.8d0, max(0.1d0, ratio))
        else
            next = dt*self%step_control%retry
        end if
        next = max(next, self%step_control%minimum)
    end function retry_length

    !> The weighted RMS of the local truncation error of the step just solved.
    !>
    !> From divided differences of the accepted history: \(h^2\,y[n+1,n,n-1]\)
    !> for backward Euler and \(\tfrac43h^3\,y[n+1,\dots,n-2]\) for BDF2,
    !> weighted by \(a_i+r\max(|y_i|,|y^n_i|)\) over the potentials. Returns -1
    !> when the history is too short to estimate it.
    function local_error(self, dt) result(error)
        implicit none
        class(type_ftcms), intent(inout) :: self
        real(real64), intent(in) :: dt
        real(real64) :: error

        Vec :: new_slope, old_slope, older_slope, curvature, weight
        PetscReal :: norm, count
        PetscErrorCode :: ierr
        real(real64) :: h1, h2

        error = -1.0d0
        if (self%history_depth < self%order) return
        h1 = self%history_steps(1)
        h2 = self%history_steps(2)
        call DMGetGlobalVector(self%dm, new_slope, ierr)
        call DMGetGlobalVector(self%dm, old_slope, ierr)
        call DMGetGlobalVector(self%dm, curvature, ierr)
        call DMGetGlobalVector(self%dm, weight, ierr)

        call divided(self%solution, self%old_solution, dt, new_slope)
        call divided(self%old_solution, self%older_solution, h1, old_slope)
        call divided(new_slope, old_slope, dt + h1, curvature)
        if (self%order == 1) then
            call VecScale(curvature, dt**2, ierr)
        else
            call DMGetGlobalVector(self%dm, older_slope, ierr)
            call divided(self%older_solution, self%oldest_solution, h2, older_slope)
            call divided(old_slope, older_slope, h1 + h2, new_slope)
            call divided(curvature, new_slope, dt + h1 + h2, old_slope)
            call VecCopy(old_slope, curvature, ierr)
            call VecScale(curvature, 4.0d0/3.0d0*dt**3, ierr)
            call DMRestoreGlobalVector(self%dm, older_slope, ierr)
        end if

        ! weight = atol + rtol*max(|y|, |y^n|)
        call VecCopy(self%solution, weight, ierr)
        call VecAbs(weight, ierr)
        call VecCopy(self%old_solution, new_slope, ierr)
        call VecAbs(new_slope, ierr)
        call VecPointwiseMax(weight, weight, new_slope, ierr)
        call VecAYPX(weight, self%error_control%relative_tolerance, self%error_tolerance, ierr)
        call VecPointwiseDivide(curvature, curvature, weight, ierr)
        call VecPointwiseMult(curvature, curvature, self%error_mask, ierr)
        call VecNorm(curvature, NORM_2, norm, ierr)
        call VecSum(self%error_mask, count, ierr)
        if (ierr /= 0) error stop "PETSc could not estimate the local error"
        error = norm/sqrt(max(count, 1.0d0))

        call DMRestoreGlobalVector(self%dm, new_slope, ierr)
        call DMRestoreGlobalVector(self%dm, old_slope, ierr)
        call DMRestoreGlobalVector(self%dm, curvature, ierr)
        call DMRestoreGlobalVector(self%dm, weight, ierr)
    contains
        !> result = (upper - lower)/span
        subroutine divided(upper, lower, span, result)
            Vec, intent(in) :: upper, lower
            real(real64), intent(in) :: span
            Vec, intent(inout) :: result
            PetscErrorCode :: code

            call VecWAXPY(result, -1.0d0, lower, upper, code)
            call VecScale(result, 1.0d0/span, code)
        end subroutine divided
    end function local_error

    !> How much constitutive work the run asked for and how much the memo saved.
    subroutine report_closure_work()
        implicit none

        integer(int64) :: requests, evaluations, asked, responses, publications

        requests = 0
        evaluations = 0
        asked = 0
        responses = 0
        publications = 0
        call closure_work_counters(requests, evaluations, asked, responses, publications)
        call global_logger%log_information("closure: responses "//long_text(responses)//" of "//long_text(asked)// &
                                           " asked, tangents "//long_text(evaluations)//" of "//long_text(requests)// &
                                           ", phase splits "//long_text(publications))
    end subroutine report_closure_work

    pure function long_text(value) result(text)
        implicit none
        integer(int64), intent(in) :: value
        character(len=:), allocatable :: text

        character(len=24) :: buffer

        write (buffer, '(i0)') value
        text = trim(buffer)
    end function long_text

    !> The earliest output time still ahead, or the end of the window.
    pure function next_output_time(self) result(time)
        implicit none
        class(type_ftcms), intent(in) :: self
        real(real64) :: time

        time = self%end_time
        if (self%field_interval > 0.0d0) time = min(time, self%next_field_time)
        if (self%history_interval > 0.0d0) time = min(time, self%next_history_time)
        if (time <= self%time) time = self%end_time
    end function next_output_time

    !> The iteration-count controller: grow an easy step, shrink a laboured one.
    pure function next_step_length(control, dt, iterations) result(next)
        implicit none
        type(type_step_control), intent(in) :: control
        real(real64), intent(in) :: dt
        integer(int32), intent(in) :: iterations
        real(real64) :: next

        next = dt
        if (iterations <= control%iterations_low) then
            next = dt*control%grow
        else if (iterations >= control%iterations_high) then
            next = dt*control%shrink
        end if
        next = min(max(next, control%minimum), control%maximum)
    end function next_step_length

    !> Solve one step of length `dt` from the accepted state, at the chosen order.
    !>
    !> Converged only when SNES converged and the solution is finite; a solve
    !> that did not converge is refused here and leaves the accepted state, the
    !> auxiliary fields and the time as they were.
    subroutine attempt_step(self, dt, converged, iterations)
        implicit none
        class(type_ftcms), intent(inout) :: self
        real(real64), intent(in) :: dt
        logical, intent(inout) :: converged
        integer(int32), intent(inout) :: iterations

        SNESConvergedReason :: reason
        PetscInt :: taken
        PetscReal :: norm
        PetscErrorCode :: ierr

        converged = .false.
        iterations = 0
        if (.not. (dt > 0.0d0) .or. .not. ieee_is_finite(dt)) error stop "A step length must be positive"

        self%dt = dt
        ierr = 0
        call store_constants(self, ierr)
        if (ierr /= 0) error stop "PETSc rejected the run constants"
        call clear_kernel_fault()
        call set_boundary_time(self%time + dt)
        if (self%num_seepage > 0) call update_seepage_bounds(self)

        call SNESSolve(self%snes, PETSC_NULL_VEC, self%solution, ierr)
        call SNESGetIterationNumber(self%snes, taken, ierr)
        iterations = taken
        call SNESGetConvergedReason(self%snes, reason, ierr)
        if (ierr == 0 .and. reason%v > 0) then
            call VecNorm(self%solution, NORM_INFINITY, norm, ierr)
            converged = ierr == 0 .and. ieee_is_finite(norm)
        end if
        if (.not. converged) then
            call report_kernel_fault()
            call discard_failed_factorization(self, reason)
            call refuse_step(self)
        end if
    end subroutine attempt_step

    !> Hold every seepage face's pressure between its data's bounds at the new time.
    !>
    !> values(2) and values(3) of a seepage condition are the drying and the
    !> saturation pressure; every other degree of freedom is unbounded.
    subroutine update_seepage_bounds(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscSection :: section
        IS :: facets
        Vec :: lower_local, upper_local
        PetscInt, pointer :: points(:), closure(:)
        PetscScalar, pointer :: lower_raw(:), upper_raw(:)
        PetscInt :: k, c, closure_size, offset, num_dofs
        PetscErrorCode :: ierr
        real(real64) :: values(BC_DATA_SIZE)
        integer(int32) :: face

        call DMGetLocalVector(self%dm, lower_local, ierr)
        call DMGetLocalVector(self%dm, upper_local, ierr)
        call VecSet(lower_local, PETSC_NINFINITY, ierr)
        call VecSet(upper_local, PETSC_INFINITY, ierr)
        call DMGetLocalSection(self%dm, section, ierr)
        call VecGetArray(lower_local, lower_raw, ierr)
        call VecGetArray(upper_local, upper_raw, ierr)
        if (ierr /= 0) error stop "PETSc could not open the variable bounds"
        do face = 1, self%num_seepage
            call boundary_data(self%seepage_slots(face), values)
            if (stratum_is_empty(self%dm, self%seepage_ids(face))) cycle
            call DMGetStratumIS(self%dm, "Face Sets", self%seepage_ids(face), facets, ierr)
            call ISGetIndices(facets, points, ierr)
            do k = 1, size(points)
                nullify (closure)
                call DMPlexGetTransitiveClosure(self%dm, points(k), PETSC_TRUE, closure_size, closure, ierr)
                do c = 1, closure_size
                    call PetscSectionGetFieldDof(section, closure(2*c - 1), self%pressure_field, num_dofs, ierr)
                    if (num_dofs <= 0) cycle
                    call PetscSectionGetFieldOffset(section, closure(2*c - 1), self%pressure_field, offset, ierr)
                    lower_raw(offset + 1:offset + num_dofs) = values(2)/self%pressure_scale
                    upper_raw(offset + 1:offset + num_dofs) = values(3)/self%pressure_scale
                end do
                call DMPlexRestoreTransitiveClosure(self%dm, points(k), PETSC_TRUE, closure_size, closure, ierr)
            end do
            call ISRestoreIndices(facets, points, ierr)
            call ISDestroy(facets, ierr)
        end do
        call VecRestoreArray(lower_local, lower_raw, ierr)
        call VecRestoreArray(upper_local, upper_raw, ierr)
        call DMLocalToGlobalBegin(self%dm, lower_local, INSERT_VALUES, self%lower_bound, ierr)
        call DMLocalToGlobalEnd(self%dm, lower_local, INSERT_VALUES, self%lower_bound, ierr)
        call DMLocalToGlobalBegin(self%dm, upper_local, INSERT_VALUES, self%upper_bound, ierr)
        call DMLocalToGlobalEnd(self%dm, upper_local, INSERT_VALUES, self%upper_bound, ierr)
        call DMRestoreLocalVector(self%dm, lower_local, ierr)
        call DMRestoreLocalVector(self%dm, upper_local, ierr)
        call SNESVISetVariableBounds(self%snes, self%lower_bound, self%upper_bound, ierr)
        if (ierr /= 0) error stop "PETSc rejected the variable bounds"
    end subroutine update_seepage_bounds

    !> Return to the accepted state.
    subroutine refuse_step(self)
        implicit none
        class(type_ftcms), intent(inout) :: self

        PetscErrorCode :: ierr

        call VecCopy(self%old_solution, self%solution, ierr)
        if (ierr == 0) call VecCopy(self%accepted_auxiliary, self%auxiliary, ierr)
        if (ierr /= 0) error stop "PETSc could not restore the accepted state"
        call set_boundary_time(self%time)
    end subroutine refuse_step

    !> Make the solved state the accepted one and shift the history.
    subroutine accept_step(self, dt)
        implicit none
        class(type_ftcms), intent(inout) :: self
        real(real64), intent(in) :: dt

        PetscErrorCode :: ierr

        call VecCopy(self%older_solution, self%oldest_solution, ierr)
        if (ierr == 0) call VecCopy(self%old_solution, self%older_solution, ierr)
        if (ierr == 0) call VecCopy(self%solution, self%old_solution, ierr)
        self%history_steps = [dt, self%history_steps(1)]
        self%history_depth = min(self%history_depth + 1, 2)
        self%time = self%time + dt
        if (ierr == 0) call publish_accepted_state(self, ierr)
        if (ierr == 0) call VecCopy(self%auxiliary, self%accepted_auxiliary, ierr)
        if (ierr /= 0) error stop "PETSc could not publish the accepted state"
    end subroutine accept_step

    !> A failed direct factorisation would otherwise be reused by the retry.
    subroutine discard_failed_factorization(self, reason)
        implicit none
        class(type_ftcms), intent(inout) :: self
        SNESConvergedReason, intent(in) :: reason

        KSP :: ksp
        PetscErrorCode :: ierr

        if (reason%v /= SNES_DIVERGED_LINEAR_SOLVE%v) return
        call SNESGetKSP(self%snes, ksp, ierr)
        if (ierr == 0) call KSPReset(ksp, ierr)
    end subroutine discard_failed_factorization

    !> Name the state the closure could not answer for, if a kernel recorded one.
    subroutine report_kernel_fault()
        implicit none

        integer(int32) :: slot
        real(real64) :: temperature, pressure

        if (.not. kernel_fault_recorded()) return
        slot = 0
        temperature = 0.0d0
        pressure = 0.0d0
        call kernel_fault_report(slot, temperature, pressure)
        call global_logger%log_warning("The closure has no state for material slot "//to_text(slot)// &
                                       " at T = "//real_text(temperature)//" C, p_w = "//real_text(pressure)//" Pa")
    end subroutine report_kernel_fault

    !> Copy the accepted potentials into the auxiliary history fields.
    !>
    !> The boundary values are inserted first, so an essential degree of
    !> freedom carries the value the step was solved with.
    module subroutine publish_accepted_state(self, ierr)
        implicit none
        class(type_ftcms), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        call DMGlobalToLocalBegin(self%dm, self%old_solution, INSERT_VALUES, self%local_solution, ierr)
        if (ierr == 0) call DMGlobalToLocalEnd(self%dm, self%old_solution, INSERT_VALUES, self%local_solution, ierr)
        if (ierr == 0) call DMPlexInsertBoundaryValues(self%dm, PETSC_TRUE, self%local_solution, &
                                                       self%time, PETSC_NULL_VEC, &
                                                       PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
        if (ierr /= 0) return
        if (self%temperature_field >= 0) call shift_history(self, AUX_TEMPERATURE_OLD, AUX_TEMPERATURE_OLDER, ierr)
        if (ierr == 0 .and. self%pressure_field >= 0) call shift_history(self, AUX_PRESSURE_OLD, AUX_PRESSURE_OLDER, ierr)
        if (ierr /= 0) return
        if (self%temperature_field >= 0) call copy_field(self, self%temperature_field, AUX_TEMPERATURE_OLD, &
                                                         self%temperature_scale, ierr)
        if (ierr == 0 .and. self%pressure_field >= 0) call copy_field(self, self%pressure_field, AUX_PRESSURE_OLD, &
                                                                      self%pressure_scale, ierr)
    end subroutine publish_accepted_state

    !> Auxiliary field `older` takes the values of auxiliary field `old`.
    subroutine shift_history(self, old, older, ierr)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: old, older
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: point, p_start, p_end, old_offset, older_offset, num_dofs

        call DMGetLocalSection(self%aux_dm, section, ierr)
        if (ierr == 0) call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr == 0) call VecGetArray(self%auxiliary, raw, ierr)
        if (ierr /= 0) return
        do point = p_start, p_end - 1
            call PetscSectionGetFieldDof(section, point, old - 1, num_dofs, ierr)
            if (num_dofs <= 0) cycle
            call PetscSectionGetFieldOffset(section, point, old - 1, old_offset, ierr)
            call PetscSectionGetFieldOffset(section, point, older - 1, older_offset, ierr)
            raw(older_offset + 1:older_offset + num_dofs) = raw(old_offset + 1:old_offset + num_dofs)
        end do
        call VecRestoreArray(self%auxiliary, raw, ierr)
    end subroutine shift_history

    !> Auxiliary field `aux` takes the physical values of solution field `field`.
    subroutine copy_field(self, field, aux, scale, ierr)
        implicit none
        class(type_ftcms), intent(inout) :: self
        PetscInt, intent(in) :: field
        integer(int32), intent(in) :: aux
        real(real64), intent(in) :: scale
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: section, aux_section
        PetscScalar, pointer :: solution_raw(:), aux_raw(:)
        PetscInt :: point, p_start, p_end, offset, aux_offset, num_dofs, aux_dofs

        call DMGetLocalSection(self%dm, section, ierr)
        if (ierr == 0) call DMGetLocalSection(self%aux_dm, aux_section, ierr)
        if (ierr == 0) call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr == 0) call VecGetArrayRead(self%local_solution, solution_raw, ierr)
        if (ierr == 0) call VecGetArray(self%auxiliary, aux_raw, ierr)
        if (ierr /= 0) return
        do point = p_start, p_end - 1
            call PetscSectionGetFieldDof(section, point, field, num_dofs, ierr)
            if (num_dofs <= 0) cycle
            call PetscSectionGetFieldDof(aux_section, point, aux - 1, aux_dofs, ierr)
            if (aux_dofs /= num_dofs) error stop "An auxiliary history field does not match its solution field"
            call PetscSectionGetFieldOffset(section, point, field, offset, ierr)
            call PetscSectionGetFieldOffset(aux_section, point, aux - 1, aux_offset, ierr)
            aux_raw(aux_offset + 1:aux_offset + num_dofs) = scale &
                                                            *solution_raw(offset + 1:offset + num_dofs)
        end do
        call VecRestoreArray(self%auxiliary, aux_raw, ierr)
        call VecRestoreArrayRead(self%local_solution, solution_raw, ierr)
    end subroutine copy_field

    !> One progress line at the `Output.json` print interval, and at the end.
    subroutine report_progress(self, iterations, force)
        implicit none
        class(type_ftcms), intent(inout) :: self
        integer(int32), intent(in) :: iterations
        logical, intent(in) :: force

        if (.not. force) then
            if (self%progress_interval <= 0.0d0) return
            if (self%time < self%next_progress_time) return
            do while (self%next_progress_time <= self%time)
                self%next_progress_time = self%next_progress_time + self%progress_interval
            end do
        end if
        call global_logger%log_information("t = "//real_text(self%time)//" s, dt = "//real_text(self%dt)// &
                                           " s, iterations "//to_text(iterations)//", accepted "// &
                                           to_text(self%accepted_steps)//", refused "//to_text(self%rejected_steps)// &
                                           ", next "//real_text(self%next_step)//" s")
    end subroutine report_progress

    pure function to_text(value) result(text)
        implicit none
        integer(int32), intent(in) :: value
        character(len=:), allocatable :: text

        character(len=16) :: buffer

        write (buffer, '(i0)') value
        text = trim(buffer)
    end function to_text

    pure function real_text(value) result(text)
        implicit none
        real(real64), intent(in) :: value
        character(len=:), allocatable :: text

        character(len=24) :: buffer

        write (buffer, '(es12.5)') value
        text = trim(adjustl(buffer))
    end function real_text

end submodule ftcms_step
