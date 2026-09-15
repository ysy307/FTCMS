#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscviewerhdf5.h>

!> The time loop, and the nonlinear solve PETSc does inside it.
!>
!> Only the stepping is here. Each step writes the previous temperature into
!> the auxiliary field and the step length into the DS constants, and then
!> hands the whole nonlinear problem to SNES, which assembles the residual and
!> the Jacobian from the discrete system.
!>
!> The step is discretised here rather than by TS because the residual
!> differences enthalpies between two known times. TS forms the time derivative
!> itself, from stages a residual written this way cannot see.
submodule(app_plexfem) plexfem_run
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: stdlib_strings, only:strip
    use :: module_core, only:type_config_time, type_config_time_ats, type_config_iteration, &
                            type_config_output_manager, OUTPUT_TYPES, &
                            get_env_string, modify_path_format
    use :: module_input, only:input_translator
    use :: control_time_ats, only:type_ats
    use :: app_plexfem_physics, only:detach_thermal_model, eval_enthalpy, eval_capacity, &
                                     eval_water_content, eval_water_capacity, &
                                     eval_enthalpy_pressure_tangent, eval_water_temperature_tangent, &
                                     eval_phase_contents, eval_effective_density, &
                                     eval_water_diffusion, eval_water_coupling, eval_gravity_flux, &
                                     A_MATERIAL_ID, A_DARCY, A_COUPLING, A_GRAVITY, A_POTENTIAL, &
                                     A_TEMPERATURE_OLD, A_PRESSURE_OLD, &
                                     A_TEMPERATURE_OLDER, A_PRESSURE_OLDER, &
                                     eval_liquid_pressure, &
                                     hydraulic_is_active, C_DT, C_PRESSURE, PLEXFEM_CONSTANTS, &
                                     C_BDF_ORDER, C_BDF_0
    implicit none

    !> How far a failed step may be cut before the run is abandoned, and by how
    !> much each time. A step that still will not converge at a thousandth of
    !> the requested length is not going to.
    real(real64), parameter :: STEP_CUT_FACTOR = 0.5d0

    !> Which nodal quantity nodal_quantity should gather.
    integer(int32), parameter :: NODAL_TEMPERATURE = 1
    integer(int32), parameter :: NODAL_PRESSURE = 2
    integer(int32), parameter :: NODAL_WATER_CONTENT = 3
    integer(int32), parameter :: NODAL_ICE_CONTENT = 4
    real(real64), parameter :: MIN_STEP_FRACTION = 1.0d-3

    !> History the local-error estimate differences.
    !>
    !> The error is measured on the conserved quantities -- nodal enthalpy and
    !> nodal stored water -- and not on T and p. A divided-difference estimator
    !> assumes what it differences is smooth in time, and temperature is not: a
    !> node crossing the freezing point has a kink in T(t), because the latent
    !> heat is absorbed at a nearly constant temperature. Its divided
    !> differences then scale like the step rather than its square, so the
    !> estimate does not fall when the step is cut and the controller rejects
    !> every step until the run dies without the nonlinear solve ever failing.
    !> Enthalpy and stored water are smooth across the transition -- that is the
    !> property the enthalpy form is chosen for -- and they are also what the
    !> residual conserves.
    type :: type_lte
        !> Mesh points this rank owns, in the order the histories are indexed.
        !> Only owned points, so the reductions below count each node once.
        PetscInt, allocatable :: point(:)
        !> Nodal volume of each, summed over the materials meeting there.
        real(real64), allocatable :: volume(:)
        !> Nodal enthalpy [J/m3] and nodal effective density [kg/m3].
        real(real64), allocatable :: energy(:), water(:)
        real(real64), allocatable :: energy_rate(:), water_rate(:)
        !> The previous second divided difference, which the second-order defect
        !> differences again.
        real(real64), allocatable :: energy_second(:), water_second(:)
        !> The step just measured, held until it is accepted or discarded.
        real(real64), allocatable :: pending_energy(:), pending_water(:)
        real(real64) :: previous_step = 0.0d0
        !> Span of the two steps before the current one, the divisor of the
        !> third divided difference.
        real(real64) :: previous_span = 0.0d0
        logical :: has_state = .false.
        logical :: has_rate = .false.
        logical :: has_second = .false.
        logical :: has_pending = .false.
    end type type_lte

    !> The instance the SNES callbacks belong to.
    !>
    !> PETSc hands a callback whatever context it was registered with, and a
    !> polymorphic Fortran variable cannot survive that round trip intact. The
    !> run owns a single instance, so the callbacks find it here instead.
    type(type_plexfem), pointer, save :: solving => null()

    !> Where the solver history is written, on rank zero.
    integer(int32), save :: history_unit = -1

    !> Fraction of the time-discretisation tolerance the nonlinear iteration is
    !> allowed to spend. An accepted step carries both errors, so holding the
    !> nonlinear part at a third of what the step controller already accepts
    !> bounds the total at 1.33 times it, inside the controller's own safety
    !> factor. Asking for more buys an accuracy the accepted step cannot use.
    real(real64), parameter :: NONLINEAR_ERROR_FRACTION = 1.0d0 / 3.0d0

    !> Largest first-order cell this path meets: eight for a hexahedron.
    integer(int32), parameter :: MAX_CELL_NODES = 8

    !> Points the transport coefficients are sampled at along a cell, to resolve
    !> the resistance of a cell the freezing front is crossing. Seventeen puts a
    !> sample every hundredth of a kelvin across a cell of the size this mesh
    !> uses, and the conductivity falls three orders in the first tenth of a
    !> kelvin below zero.
    integer(int32), parameter :: SERIES_SAMPLES = 17

    !> How far the ends of a cell must differ before its interior is sampled at
    !> all. Below it the two-sample series equivalent already resolves the
    !> integral to within the difference between the rules.
    real(real64), parameter :: SERIES_CONTRAST = 2.0d0

    !> What counts as a large change of each unknown, and so the unit each is
    !> measured in once the system is scaled.
    !>
    !> Fifty kilopascals is about 0.04 K of Clapeyron freezing-point shift, so
    !> the two constrain the phase state by comparable amounts -- which is what
    !> makes a single bound on the largest scaled component mean the same thing
    !> for both fields.
    real(real64), parameter :: UNKNOWN_SCALE(2) = [5.0d-1, 5.0d4]

contains

    module subroutine run_plexfem(self)
        implicit none
        !> A target because the SNES callbacks reach this instance through a
        !> module pointer rather than through PETSc's context argument.
        class(type_plexfem), intent(inout), target :: self

        type(type_config_time) :: config_time
        type(type_config_time_ats) :: config_ats
        type(type_config_iteration) :: config_iteration
        type(type_ats) :: ats
        type(type_lte) :: lte
        type(type_config_output_manager) :: config_field
        PetscErrorCode :: ierr
        SNESConvergedReason :: reason
        real(real64) :: time, end_time, step, requested_step, min_step, max_step
        real(real64) :: output_interval, next_output
        real(real64) :: error, next_step, previous_step
        integer(int32) :: step_number, output_number, attempt, target_order
        PetscInt :: total_dofs
        character(32) :: status
        PetscInt :: iterations
        logical :: converged
        PetscMPIInt :: rank
        Vec :: residual
        Mat :: jacobian

        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)

        call input_translator%execute(self%input, config_time)
        call input_translator%execute(self%input, config_ats)
        call input_translator%execute(self%input, config_iteration)
        call input_translator%execute(self%input, OUTPUT_TYPES%FIELD, config_field)
        output_interval = config_field%interval_val * config_field%interval_unit%value

        call SNESCreate(PETSC_COMM_WORLD, self%snes, ierr)
        if (ierr /= 0) error stop "PETSc could not create the nonlinear solver"
        self%snes_ready = .true.

        call SNESSetDM(self%snes, self%dm, ierr)
        if (ierr /= 0) error stop "PETSc could not give the nonlinear solver the mesh"

        solving => self
        call DMCreateGlobalVector(self%dm, self%row_scale, ierr)
        if (ierr == 0) call VecSet(self%row_scale, 1.0_PETSC_REAL_KIND, ierr)
        if (ierr /= 0) error stop "PETSc could not create the row scale"
        call build_unknown_scale(self, ierr)
        if (ierr /= 0) error stop "PETSc could not scale the unknowns"
        ! The initial state was built in physical units, before the scale
        ! existed. From here on the solver's vectors are scaled ones.
        call VecPointwiseDivide(self%solution, self%solution, self%unknown_scale, ierr)
        if (ierr == 0) call VecCopy(self%solution, self%old_solution, ierr)
        if (ierr == 0) call VecCopy(self%solution, self%older_solution, ierr)
        if (ierr /= 0) error stop "PETSc could not put the initial state in solver units"
        call DMCreateGlobalVector(self%dm, residual, ierr)
        if (ierr == 0) call DMCreateMatrix(self%dm, jacobian, ierr)
        if (ierr /= 0) error stop "PETSc could not create the residual and the matrix"

        call SNESSetFunction(self%snes, residual, form_function, self%dm, ierr)
        if (ierr == 0) call SNESSetJacobian(self%snes, jacobian, jacobian, form_jacobian, self%dm, ierr)
        if (ierr /= 0) error stop "PETSc rejected the residual and Jacobian"

        ! A step is never accepted without solving it. The scaled residual
        ! carries the step length, so at a short step the very first iterate can
        ! already sit inside the tolerance and SNES would return having done
        ! nothing -- time advances, the state does not, and the run walks off the
        ! range the constitutive models are defined on.
        call SNESSetForceIteration(self%snes, PETSC_TRUE, ierr)
        if (ierr /= 0) error stop "PETSc rejected the forced iteration"

        ! Converged when the project's own criterion is met: the residual is
        ! scaled so that it reads as a multiple of the tolerance at each node, and
        ! its weighted-RMS reaching one is the acceptance the older path applies.
        ! The two-norm PETSc tests is that RMS times the square root of the count.
        call VecGetSize(residual, total_dofs, ierr)
        if (ierr == 0) call SNESSetTolerances(self%snes, sqrt(real(total_dofs, real64)), &
                                              0.0d0, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, &
                                              PETSC_DEFAULT_INTEGER, ierr)
        if (ierr /= 0) error stop "PETSc rejected the convergence tolerances"

        ! How far an update may reach is left to PETSc's trust region, measured
        ! in the units the unknowns are scaled to.
        !
        ! A Newton direction can be sound and its full length still meaningless:
        ! measured without any bound, single updates reached tens of kelvin and
        ! megapascals, far outside the range the linearisation they came from
        ! describes, and the phase state then left the range the constitutive
        ! models are defined on at all.
        !
        ! Both unknowns are scaled by what counts as a large change of each --
        ! see UNKNOWN_SCALE -- so the largest component of a scaled update means
        ! the same thing for either field, which is what makes a region measured
        ! in the infinity norm the right shape. A fixed bound is deliberately NOT
        ! imposed: the radius starts at one scale unit and the trust region grows
        ! or shrinks it on how well the model predicted the last step, which is
        ! the whole reason to use one. Capping it at the starting value instead
        ! made it bind on every iteration and a cold-shock start crawled to its
        ! answer half a kelvin at a time.
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-snes_type", "newtontr", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-snes_tr_norm_type", "infinity", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-snes_tr_delta0", "1.0", ierr)
        if (ierr == 0) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-snes_tr_fallback_type", "newton", ierr)
        if (ierr /= 0) error stop "PETSc rejected the trust-region options"

        ! Last, so every option above can be overridden from the command line.
        call SNESSetFromOptions(self%snes, ierr)
        if (ierr /= 0) error stop "PETSc rejected the nonlinear solver options"

        time = config_time%start_time
        end_time = config_time%end_time
        requested_step = config_time%initial_step

        ! Adaptive stepping. The policy -- the bounds, the factors, the PI gains
        ! and the tolerance -- is the project's, and control_time_ats already
        ! implements it for the older path, so it is used here rather than
        ! restated: what is added below is only the estimate it controls on and
        ! the choice between its branches.
        call ats%initialize(config_ats)
        if (ats%is_active()) then
            min_step = config_ats%dt_min
            max_step = config_ats%dt_max
        else
            min_step = requested_step * MIN_STEP_FRACTION
            max_step = requested_step
        end if
        ! The owned-node list and the nodal volumes behind it are what the row
        ! scale is built from as well, so it is prepared whether or not the step
        ! size is being controlled.
        call lte_prepare(self, lte, ierr)
        if (ierr /= 0) error stop "The nodal conserved quantities could not be prepared"

        step_number = 0
        output_number = 0
        attempt = 0
        status = ""
        previous_step = 0.0d0
        ! The first step has only one state behind it, whatever order the project
        ! asks for; the ramp is what makes the request reachable.
        target_order = 1
        next_output = time
        call open_history(self)

        call open_output(self)
        call write_field(self, time, output_number, next_output, output_interval)

        do while (time < end_time)
            step = min(requested_step, end_time - time)
            ! Land exactly on the next output time. The validation picks a
            ! snapshot by round(time/interval), so a step that oversteps a mark
            ! writes the wrong time under that index, and one that oversteps two
            ! leaves an index without a snapshot at all. Clipping here costs one
            ! short step per output and keeps the file's index the time it says.
            if (output_interval > 0.0d0 .and. next_output > time) then
                step = min(step, next_output - time)
            end if

            ! The states this step differences against are the last two the run
            ! accepted. They are shifted once per step, outside the retry loop,
            ! so a retry differences against the same history the failed attempt
            ! did rather than against the attempt itself.
            call VecCopy(self%old_solution, self%older_solution, ierr)
            if (ierr == 0) call VecCopy(self%solution, self%old_solution, ierr)
            if (ierr /= 0) error stop "PETSc could not keep the previous states"
            call update_history_fields(self, ierr)
            if (ierr /= 0) error stop "PETSc could not publish the previous states"

            converged = .false.
            do while (.not. converged)
                call set_time_step(self, step, previous_step, target_order)
                call build_row_scale(self, lte, step, config_iteration%nonlinear%atol_enthalpy, &
                                     config_iteration%nonlinear%atol_density, &
                                     config_ats%error_rtol, ierr)
                if (ierr /= 0) error stop "PETSc could not scale the residual rows"
                attempt = attempt + 1

                call SNESSolve(self%snes, PETSC_NULL_VEC, self%solution, ierr)
                if (ierr /= 0) error stop "The nonlinear solve failed"
                call SNESGetConvergedReason(self%snes, reason, ierr)
                call SNESGetIterationNumber(self%snes, iterations, ierr)
                if (ierr /= 0) error stop "PETSc could not report the iteration count"

                converged = int(reason%v, int32) > 0
                error = -1.0d0

                ! An accuracy failure is a rejection too. A step Newton solved
                ! is not thereby a step worth keeping: the estimate is what says
                ! whether the answer is one, and a converged step whose error
                ! exceeds the tolerance is resized at the order the estimate was
                ! measured at rather than by the blunter robustness factor.
                if (converged .and. ats%is_active()) then
                    call lte_measure(self%dm, lte, step, self%bdf_order, config_ats%error_rtol, error, ierr)
                    if (ierr /= 0) error stop "The local-error estimate failed"
                    if (ats%use_error_control .and. error > 1.0d0) then
                        converged = .false.
                        status = "lte_rejected"
                        call ats%lte_retry_dt(step, error, next_step, error_order=self%bdf_order)
                    end if
                end if

                if (converged) then
                    status = "accepted"
                else
                    if (status /= "lte_rejected") then
                        status = "nonlinear_rejected"
                        if (ats%is_active()) then
                            call ats%calc_retry_dt(step, next_step)
                        else
                            next_step = step * STEP_CUT_FACTOR
                        end if
                    end if
                end if

                ! An accepted attempt is recorded after the controller has
                ! chosen the next step, since that is one of the columns; a
                ! rejected one already knows its retry length.
                if (.not. converged) then
                    call log_attempt(self, attempt, step_number, time, step, next_step, &
                                     .false., status, int(iterations, int32), error)
                end if

                if (.not. converged) then
                    ! The retry starts from the last accepted state, not from the
                    ! iterate the failed solve left behind: SNES writes its own
                    ! progress into the solution vector, and a diverged iterate is
                    ! a worse initial guess than the state the step begins at.
                    call VecCopy(self%old_solution, self%solution, ierr)
                    if (ierr /= 0) error stop "PETSc could not restart the step"
                    ! A rejected step's history describes lengths that were not
                    ! taken, so the order drops back to one for the retry and
                    ! climbs again on the next clean step.
                    target_order = 1
                    step = next_step
                    ! The PI term may not grow the step until one clean step has
                    ! been taken after this: the error history now describes a
                    ! step that was never accepted.
                    if (ats%is_active()) ats%after_rejection = .true.
                    ! And the next step starts from what worked, not from what
                    ! has just been shown not to.
                    requested_step = step
                    status = ""
                    if (step < min_step) then
                        if (rank == 0) write (*, '(a,es12.5,a)') &
                            "The step could not be made to converge above ", min_step, " s"
                        error stop "The time integration failed"
                    end if
                end if
            end do

            ! What the next step should be. Accuracy leads and the iteration
            ! count is only a brake: a step the nonlinear solve found hard has
            ! its growth capped even when the error allows it, because the
            ! estimate says nothing about whether Newton will converge again.
            ! What was just accepted becomes the history the next step uses, so
            ! from here on a second-order step has the two states it needs.
            previous_step = step
            target_order = min(config_time%target_bdf_order, 2)

            if (ats%is_active()) then
                call lte_commit(lte, step, self%bdf_order)

                if (ats%use_error_control .and. error > 0.0d0) then
                    call ats%pi_controller_dt(step, error, next_step)
                    if (int(iterations, int32) >= ats%iter_max .and. step > min_step) then
                        next_step = min(next_step, step * ats%scale_down)
                    end if
                else
                    call ats%predict_next_dt(step, int(iterations, int32), next_step)
                end if
                requested_step = min(max(next_step, min_step), max_step)
            end if
            call log_attempt(self, attempt, step_number, time, step, requested_step, &
                             .true., "accepted", int(iterations, int32), error)

            time = time + step
            step_number = step_number + 1
            ! Progress is one line per step, not one per rank.
            if (rank == 0) then
                write (*, '(a,i0,a,es12.5,a,es12.5)') "step ", step_number, &
                    "  t = ", time, "  dt = ", step
            end if

            call write_field(self, time, output_number, next_output, output_interval)
        end do

        if (self%output_ready) then
            call PetscViewerHDF5PopTimestepping(self%output_viewer, ierr)
            call PetscViewerDestroy(self%output_viewer, ierr)
            self%output_ready = .false.
        end if

        call VecDestroy(residual, ierr)
        call MatDestroy(jacobian, ierr)
        nullify (solving)
    end subroutine run_plexfem

    !> List the nodes this rank owns and the volume each carries.
    !>
    !> Done once: the mesh does not move, so neither does this. Ownership comes
    !> from the global section, which gives a negative offset at a point another
    !> rank owns, so a shared node is counted by exactly one rank.
    subroutine lte_prepare(self, lte, ierr)
        implicit none
        class(type_plexfem), intent(in) :: self
        type(type_lte), intent(inout) :: lte
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: section, global_section
        PetscScalar, pointer :: volume_raw(:)
        PetscInt :: point, p_start, p_end, offset, global_offset, num_dofs
        integer(int32) :: m, count, k

        call DMGetLocalSection(self%dm, section, ierr)
        if (ierr == 0) call DMGetGlobalSection(self%dm, global_section, ierr)
        if (ierr /= 0) return
        call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr /= 0) return

        count = 0
        do point = p_start, p_end - 1
            call PetscSectionGetFieldDof(section, point, 0_PETSC_INT_KIND, num_dofs, ierr)
            if (ierr /= 0 .or. num_dofs <= 0) cycle
            call PetscSectionGetOffset(global_section, point, global_offset, ierr)
            if (ierr /= 0 .or. global_offset < 0) cycle
            count = count + 1
        end do

        allocate (lte%point(count), lte%volume(count))
        allocate (lte%energy(count), lte%water(count))
        allocate (lte%energy_rate(count), lte%water_rate(count))
        allocate (lte%pending_energy(count), lte%pending_water(count))
        allocate (lte%energy_second(count), lte%water_second(count))
        lte%volume(:) = 0.0d0

        k = 0
        do point = p_start, p_end - 1
            call PetscSectionGetFieldDof(section, point, 0_PETSC_INT_KIND, num_dofs, ierr)
            if (ierr /= 0 .or. num_dofs <= 0) cycle
            call PetscSectionGetOffset(global_section, point, global_offset, ierr)
            if (ierr /= 0 .or. global_offset < 0) cycle
            k = k + 1
            lte%point(k) = point
        end do

        do m = 1, self%num_materials
            call VecGetArrayRead(self%nodal_volume(m), volume_raw, ierr)
            if (ierr /= 0) return
            do k = 1, count
                call PetscSectionGetFieldOffset(section, lte%point(k), 0_PETSC_INT_KIND, offset, ierr)
                if (ierr /= 0) cycle
                lte%volume(k) = lte%volume(k) + real(volume_raw(offset + 1), real64)
            end do
            call VecRestoreArrayRead(self%nodal_volume(m), volume_raw, ierr)
        end do
    end subroutine lte_prepare

    !> The conserved quantities at the accepted solution, node by node.
    !>
    !> Nodal enthalpy [J/m3] and nodal effective density [kg/m3]: the two
    !> quantities the project states its tolerances in, so the same pair serves
    !> the local-error estimate and the row scale.
    !>
    !> Each is the volume-weighted mean over the materials meeting at the node,
    !> so what is differenced is the same density the storage term integrates.
    subroutine lte_sample(dm, lte, energy, water, ierr)
        implicit none
        DM :: dm
        type(type_lte), intent(in) :: lte
        real(real64), intent(inout) :: energy(:), water(:)
        PetscErrorCode, intent(inout) :: ierr

        Vec :: local_x
        PetscSection :: section
        PetscScalar, pointer :: x_raw(:), volume_raw(:)
        PetscInt :: offset, pressure_offset
        real(real64) :: value, temperature, pressure, weight
        integer(int32) :: m, k, material_id
        logical :: water_active

        water_active = hydraulic_is_active()
        energy(:) = 0.0d0
        water(:) = 0.0d0

        call DMGetLocalSection(dm, section, ierr)
        if (ierr /= 0) return
        call DMGetLocalVector(dm, local_x, ierr)
        if (ierr /= 0) return
        call gather(dm, solving%solution, local_x, ierr)
        if (ierr /= 0) return
        call VecGetArrayRead(local_x, x_raw, ierr)
        if (ierr /= 0) return

        do m = 1, solving%num_materials
            material_id = solving%material_id(m)
            call VecGetArrayRead(solving%nodal_volume(m), volume_raw, ierr)
            if (ierr /= 0) return

            do k = 1, size(lte%point)
                call PetscSectionGetFieldOffset(section, lte%point(k), 0_PETSC_INT_KIND, offset, ierr)
                if (ierr /= 0) cycle
                weight = real(volume_raw(offset + 1), real64)
                if (weight == 0.0d0) cycle

                temperature = real(x_raw(offset + 1), real64)
                pressure = real(solving%constants(C_PRESSURE), real64)
                if (water_active) then
                    call PetscSectionGetFieldOffset(section, lte%point(k), 1_PETSC_INT_KIND, pressure_offset, ierr)
                    if (ierr /= 0) cycle
                    pressure = real(x_raw(pressure_offset + 1), real64)
                end if

                call eval_enthalpy(temperature, pressure, material_id, solving%constants, value)
                energy(k) = energy(k) + weight * value
                if (.not. water_active) cycle
                call eval_effective_density(temperature, pressure, material_id, solving%constants, value)
                water(k) = water(k) + weight * value
            end do

            call VecRestoreArrayRead(solving%nodal_volume(m), volume_raw, ierr)
        end do

        call VecRestoreArrayRead(local_x, x_raw, ierr)
        call DMRestoreLocalVector(dm, local_x, ierr)

        do k = 1, size(lte%point)
            if (lte%volume(k) <= 0.0d0) cycle
            energy(k) = energy(k) / lte%volume(k)
            water(k) = water(k) / lte%volume(k)
        end do
    end subroutine lte_sample

    !> Normalised local truncation error of the step just accepted, and the
    !> history rolled forward onto it.
    !>
    !> Newton divided differences on the accepted times give, for the first
    !> order this loop integrates at,
    !>
    !>   ydot   = (y - y_prev)/h
    !>   d2     = (ydot - ydot_prev)/(h + h_prev)
    !>   defect = h^2 d2                       ~ LTE_1
    !>
    !> and the estimate is the volume-weighted RMS of the defect against a
    !> weight of rtol times the field's own magnitude. A value at or below one
    !> means the step met the accuracy target. It is negative until two steps
    !> have been accepted, which is what tells the caller there is no estimate
    !> to control on yet.
    subroutine lte_measure(dm, lte, step, order, rtol, error, ierr)
        implicit none
        DM :: dm
        type(type_lte), intent(inout) :: lte
        real(real64), intent(in) :: step
        !> Order the step was taken at. The estimator has to match it: the
        !> first-order defect measured under a second-order step is smaller than
        !> the true error by a factor of the step, so the controller would accept
        !> steps that violate the tolerance and never learn that it had.
        integer(int32), intent(in) :: order
        real(real64), intent(in) :: rtol
        real(real64), intent(inout) :: error
        PetscErrorCode, intent(inout) :: ierr

        real(real64), allocatable :: energy(:), water(:)
        real(real64) :: error_energy, error_water
        real(real64) :: span_send(4), span(4)
        integer(int32) :: count, mpi_err, measured
        logical :: water_active

        error = -1.0d0
        count = size(lte%point)
        water_active = hydraulic_is_active()

        allocate (energy(count), water(count))
        call lte_sample(dm, lte, energy, water, ierr)
        if (ierr /= 0) return

        if (.not. lte%has_state) then
            lte%pending_energy(:) = energy(:)
            lte%pending_water(:) = water(:)
            lte%has_pending = .true.
            return
        end if

        ! The span of the previous state floors the relative weight. A purely
        ! relative tolerance is meaningless for a quantity whose zero is not a
        ! natural zero, and it is the whole domain's span that makes the
        ! criterion uniform rather than one to two orders stricter at the front
        ! than anywhere else. Reduced with MPI_MAX on (max, -min).
        span_send(1) = maxval(lte%energy)
        span_send(2) = -minval(lte%energy)
        span_send(3) = maxval(lte%water)
        span_send(4) = -minval(lte%water)
        call MPI_Allreduce(span_send, span, 4, MPI_DOUBLE_PRECISION, MPI_MAX, PETSC_COMM_WORLD, mpi_err)

        ! Second order needs the previous second difference as well; until it
        ! exists the first-order defect is used, which is a valid bound because
        ! the run is still first order there anyway.
        measured = order
        if (measured >= 2 .and. .not. lte%has_second) measured = 1
        if (measured >= 2 .and. step + lte%previous_span <= 0.0d0) measured = 1

        error_energy = -1.0d0
        error_water = -1.0d0
        if (lte%has_rate .and. step + lte%previous_step > 0.0d0) then
            call weighted_rms(energy, lte%energy, lte%energy_rate, lte%energy_second, &
                              max(0.0d0, span(1) + span(2)), measured, error_energy)
            if (water_active) then
                call weighted_rms(water, lte%water, lte%water_rate, lte%water_second, &
                                  max(0.0d0, span(3) + span(4)), measured, error_water)
            end if
            error = max(error_energy, error_water)
        end if

        ! The history is NOT rolled forward here. A step whose error exceeds the
        ! tolerance is rejected and retried at a shorter length, and its state
        ! must not become the state the retry differences against.
        lte%pending_energy(:) = energy(:)
        lte%pending_water(:) = water(:)
        lte%has_pending = .true.

    contains

        !> Volume-weighted RMS of the first-order defect over the whole problem.
        !>
        !> One norm for every rank: a per-rank norm would hand each of them a
        !> different step and desynchronise the run. Weighted by nodal volume
        !> rather than node count, so a finer graded region does not inflate it.
        subroutine weighted_rms(y, y_prev, rate_prev, second_prev, y_span, order, norm)
            implicit none
            real(real64), intent(in) :: y(:), y_prev(:), rate_prev(:), second_prev(:)
            real(real64), intent(in) :: y_span
            integer(int32), intent(in) :: order
            real(real64), intent(inout) :: norm

            real(real64) :: rate, second, defect, weight
            real(real64) :: accum_send(2), accum(2)
            integer(int32) :: i

            accum_send(:) = 0.0d0
            do i = 1, size(y)
                ! Newton divided differences on the accepted times:
                !   y[t_n,t_{n-1}]         = rate
                !   y[t_n,..,t_{n-2}]      = (rate - rate_prev)/(h + h_prev)
                !   y[t_n,..,t_{n-3}]      = (second - second_prev)/(h + span_prev)
                ! and LTE_k ~ h**(k+1) times the (k+1)-th of them.
                rate = (y(i) - y_prev(i)) / step
                second = (rate - rate_prev(i)) / (step + lte%previous_step)
                if (order >= 2) then
                    defect = step**3 * (second - second_prev(i)) / (step + lte%previous_span)
                else
                    defect = step * step * second
                end if
                weight = rtol * max(abs(y(i)), abs(y_prev(i)), y_span)
                if (weight <= tiny(1.0d0)) cycle
                accum_send(1) = accum_send(1) + lte%volume(i) * (defect / weight)**2
                accum_send(2) = accum_send(2) + lte%volume(i)
            end do

            call MPI_Allreduce(accum_send, accum, 2, MPI_DOUBLE_PRECISION, MPI_SUM, &
                               PETSC_COMM_WORLD, mpi_err)
            norm = -1.0d0
            if (accum(2) > 0.0d0) norm = sqrt(accum(1) / accum(2))
        end subroutine weighted_rms
    end subroutine lte_measure

    !> Accept the measured step: its state becomes the one the next is
    !> differenced against, and the rate over it becomes the previous rate.
    subroutine lte_commit(lte, step, order)
        implicit none
        type(type_lte), intent(inout) :: lte
        real(real64), intent(in) :: step
        integer(int32), intent(in) :: order

        real(real64) :: rate
        integer(int32) :: k

        if (.not. lte%has_pending) return

        if (lte%has_state) then
            do k = 1, size(lte%point)
                rate = (lte%pending_energy(k) - lte%energy(k)) / step
                if (lte%has_rate) then
                    lte%energy_second(k) = (rate - lte%energy_rate(k)) / (step + lte%previous_step)
                end if
                lte%energy_rate(k) = rate

                rate = (lte%pending_water(k) - lte%water(k)) / step
                if (lte%has_rate) then
                    lte%water_second(k) = (rate - lte%water_rate(k)) / (step + lte%previous_step)
                end if
                lte%water_rate(k) = rate
            end do
            if (lte%has_rate) then
                lte%previous_span = step + lte%previous_step
                lte%has_second = .true.
            end if
            lte%has_rate = .true.
        end if

        lte%energy(:) = lte%pending_energy(:)
        lte%water(:) = lte%pending_water(:)
        lte%previous_step = step
        lte%has_state = .true.
        lte%has_pending = .false.
    end subroutine lte_commit

    !> The residual PETSc's nonlinear solver asks for.
    !>
    !> The conduction and boundary terms are integrated by PETSc from the
    !> discrete system; the storage term is added here, node by node, weighted
    !> by the nodal volume of each material. That is what lumping is: the
    !> enthalpy is differenced at the node rather than at the quadrature
    !> points, so a strongly cooled boundary cannot drive its neighbours past
    !> the temperature they started from.
    subroutine form_function(snes, x, f, dm, ierr)
        implicit none
        SNES :: snes
        Vec :: x
        Vec :: f
        DM :: dm
        PetscErrorCode :: ierr

        Vec :: local_x, local_f

        ierr = 0
        if (.not. associated(solving)) then
            ierr = 1
            return
        end if

        call DMGetLocalVector(dm, local_x, ierr)
        if (ierr == 0) call DMGetLocalVector(dm, local_f, ierr)
        if (ierr /= 0) return

        call gather(dm, x, local_x, ierr)
        ! The element coefficients belong to this iterate, so they are refreshed
        ! before PETSc integrates anything with them.
        if (ierr == 0) call update_element_coefficients(dm, local_x, ierr)
        if (ierr == 0) call VecZeroEntries(local_f, ierr)
        if (ierr == 0) call DMPlexSNESComputeResidualFEM(dm, local_x, local_f, dm, ierr)
        if (ierr /= 0) return

        call VecZeroEntries(f, ierr)
        if (ierr == 0) call DMLocalToGlobalBegin(dm, local_f, ADD_VALUES, f, ierr)
        if (ierr == 0) call DMLocalToGlobalEnd(dm, local_f, ADD_VALUES, f, ierr)
        ! Both equations in the units of their own primary variable's tolerance,
        ! so the norm SNES tests and the merit its line search minimises see the
        ! water equation at all.
        if (ierr == 0) call VecPointwiseMult(f, solving%row_scale, f, ierr)

        call DMRestoreLocalVector(dm, local_f, ierr)
        call DMRestoreLocalVector(dm, local_x, ierr)
    end subroutine form_function

    !> Scatter a global vector into a local one, prescribed values included.
    subroutine gather(dm, global, local, ierr)
        implicit none
        DM :: dm
        Vec :: global
        Vec :: local
        PetscErrorCode, intent(inout) :: ierr

        call DMGlobalToLocalBegin(dm, global, INSERT_VALUES, local, ierr)
        if (ierr == 0) call DMGlobalToLocalEnd(dm, global, INSERT_VALUES, local, ierr)
        ! Back into physical units. The solver's vectors are scaled; everything
        ! downstream of here -- the kernels, the constitutive models, the
        ! prescribed boundary values inserted just below -- is not, and this is
        ! the only place the two meet.
        if (ierr == 0) call VecPointwiseMult(local, solving%unknown_scale_local, local, ierr)
        ! The prescribed values live only in the local vector.
        ! Time zero is safe: a time-dependent essential boundary is refused
        ! when the boundaries are registered.
        if (ierr == 0) call DMPlexInsertBoundaryValues(dm, PETSC_TRUE, local, 0.0_PETSC_REAL_KIND, &
                                                       PETSC_NULL_VEC, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
    end subroutine gather

    !> The Jacobian of the residual above.
    !>
    !> Lumping removes the coupling between neighbouring nodes, not the coupling
    !> between the two fields at one node: the residual differences U(T, p) and
    !> Theta(T, p), so its derivative is a full 2x2 block per node. Writing only
    !> the diagonal of it is what left the energy row without dU/dp, and dU/dp
    !> carries latent heat -- the term is not a correction to the block, it is
    !> comparable to the block. That showed up as ||J - Jfd||_F growing exactly
    !> in proportion to 1/dt.
    subroutine form_jacobian(snes, x, jac, jac_prec, dm, ierr)
        implicit none
        SNES :: snes
        Vec :: x
        Mat :: jac
        Mat :: jac_prec
        DM :: dm
        PetscErrorCode :: ierr

        Vec :: local_x

        ierr = 0
        if (.not. associated(solving)) then
            ierr = 1
            return
        end if

        call DMGetLocalVector(dm, local_x, ierr)
        if (ierr /= 0) return
        call gather(dm, x, local_x, ierr)
        if (ierr == 0) call update_element_coefficients(dm, local_x, ierr)
        if (ierr == 0) call MatZeroEntries(jac_prec, ierr)
        ! Everything is assembled into the preconditioner, and the operator is
        ! only written when SNES has given us the same matrix for both. With
        ! -snes_mf_operator the operator is a matrix-free shell that owns the
        ! true Jacobian action and cannot take values at all, so writing to it
        ! unconditionally is what made that option -- the one diagnostic that
        ! separates a wrong tangent from a hard problem -- unusable here.
        if (ierr == 0) call DMPlexSNESComputeJacobianFEM(dm, local_x, jac_prec, jac_prec, dm, ierr)
        if (ierr /= 0) return

        call DMRestoreLocalVector(dm, local_x, ierr)

        call MatAssemblyBegin(jac_prec, MAT_FINAL_ASSEMBLY, ierr)
        if (ierr == 0) call MatAssemblyEnd(jac_prec, MAT_FINAL_ASSEMBLY, ierr)
        if (ierr /= 0) return

        ! The same scaling the residual carries. The Newton direction is
        ! unchanged by it -- R J d = -R F has the solution J d = -F for any
        ! invertible R -- so this changes what is measured, not where the step
        ! goes: the linear solve's own tolerance, the zero-pivot test, and the
        ! merit function all start seeing both physics.
        ! Rows into tolerances, columns into the units the unknowns are solved
        ! in: what PETSc is handed is R J S, the derivative of the scaled
        ! residual with respect to the scaled unknowns.
        call MatDiagonalScale(jac_prec, solving%row_scale, solving%unknown_scale, ierr)
        if (ierr /= 0) return

        ! The operator is assembled as well when it is a separate object. No
        ! values were put into it, but a matrix-free operator takes the point it
        ! differences about from this call, so skipping it leaves it without a
        ! base rather than merely unassembled.
        if (jac%v /= jac_prec%v) then
            call MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY, ierr)
            if (ierr == 0) call MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY, ierr)
        end if
    end subroutine form_jacobian

    !> Build the row scale of the residual for the step about to be solved.
    !>
    !> This is the project's own acceptance criterion, expressed as a scaling of
    !> the rows rather than as a test applied after the fact:
    !>
    !>   e_i = R_i dt / (V_i w_i),   w_i = f (atol_q + rtol max(|q_i|, span_q))
    !>
    !> with q the block's conserved quantity -- enthalpy for the energy rows,
    !> effective density for the water rows -- and f the fraction of the
    !> time-discretisation tolerance the nonlinear iteration is allowed to spend.
    !> Scaled that way the residual is dimensionless in both blocks and reads as
    !> a multiple of the tolerance, so one norm can see both physics: unscaled
    !> they are an energy rate and a mass rate that differ by some eighteen
    !> orders on this problem, and the norm SNES tests would be the energy rows
    !> to machine precision.
    !>
    !> The weight is deliberately NOT the storage capacity. Dividing by dH/dT or
    !> dTheta/dp scores the imbalance in a currency the node may not be able to
    !> pay in: where the pore volume binds, the split is fixed by the constraint
    !> and dH/dT falls back to the sensible value, so an imbalance the node would
    !> settle through pressure is reported as a large temperature error instead.
    !> The older path records measuring exactly that at its own wall.
    !>
    !> It is built once per step from the state the step starts at and held fixed
    !> for the whole nonlinear solve. A scale that moved with the iterate would
    !> change the function the line search is minimising in the middle of
    !> minimising it.
    subroutine build_row_scale(self, lte, step, atol_energy, atol_water, rtol, ierr)
        implicit none
        class(type_plexfem), intent(inout) :: self
        type(type_lte), intent(in) :: lte
        real(real64), intent(in) :: step
        real(real64), intent(in) :: atol_energy, atol_water, rtol
        PetscErrorCode, intent(inout) :: ierr

        Vec :: local_scale
        PetscSection :: section
        PetscScalar, pointer :: scale_raw(:)
        PetscInt :: offset, pressure_offset
        real(real64), allocatable :: energy(:), water(:)
        real(real64) :: span_send(4), span(4)
        real(real64) :: energy_span, water_span, weight
        integer(int32) :: k, count, mpi_err
        logical :: water_active

        water_active = hydraulic_is_active()
        count = size(lte%point)

        allocate (energy(count), water(count))
        call lte_sample(self%dm, lte, energy, water, ierr)
        if (ierr /= 0) return

        ! The span of each conserved quantity floors its relative weight, so the
        ! criterion is uniform over the domain rather than one to two orders
        ! stricter wherever the quantity happens to pass through zero.
        span_send(1) = maxval(energy)
        span_send(2) = -minval(energy)
        span_send(3) = maxval(water)
        span_send(4) = -minval(water)
        call MPI_Allreduce(span_send, span, 4, MPI_DOUBLE_PRECISION, MPI_MAX, PETSC_COMM_WORLD, mpi_err)
        energy_span = max(0.0d0, span(1) + span(2))
        water_span = max(0.0d0, span(3) + span(4))

        call DMGetLocalSection(self%dm, section, ierr)
        if (ierr /= 0) return
        call DMGetLocalVector(self%dm, local_scale, ierr)
        if (ierr == 0) call VecZeroEntries(local_scale, ierr)
        if (ierr /= 0) return
        call VecGetArray(local_scale, scale_raw, ierr)
        if (ierr /= 0) return

        do k = 1, count
            call PetscSectionGetFieldOffset(section, lte%point(k), 0_PETSC_INT_KIND, offset, ierr)
            if (ierr /= 0) cycle

            weight = NONLINEAR_ERROR_FRACTION * (atol_energy + rtol * max(abs(energy(k)), energy_span))
            scale_raw(offset + 1) = lte%volume(k) * weight / step

            if (.not. water_active) cycle
            call PetscSectionGetFieldOffset(section, lte%point(k), 1_PETSC_INT_KIND, pressure_offset, ierr)
            if (ierr /= 0) cycle
            weight = NONLINEAR_ERROR_FRACTION * (atol_water + rtol * max(abs(water(k)), water_span))
            scale_raw(pressure_offset + 1) = lte%volume(k) * weight / step
        end do

        call VecRestoreArray(local_scale, scale_raw, ierr)

        ! Only the owned nodes were written, so this is an insert rather than the
        ! sum a per-rank quantity would need: lte%volume already holds the whole
        ! node's volume.
        call VecZeroEntries(self%row_scale, ierr)
        if (ierr == 0) call DMLocalToGlobalBegin(self%dm, local_scale, INSERT_VALUES, self%row_scale, ierr)
        if (ierr == 0) call DMLocalToGlobalEnd(self%dm, local_scale, INSERT_VALUES, self%row_scale, ierr)
        if (ierr == 0) call DMRestoreLocalVector(self%dm, local_scale, ierr)
        if (ierr /= 0) return

        ! The residual is multiplied by this, so what is stored is the reciprocal
        ! of the denominator above.
        call VecReciprocal(self%row_scale, ierr)
    end subroutine build_row_scale

    !> Copy the accepted history into the auxiliary fields the storage term reads.
    !>
    !> Once per step, not once per residual: what these hold does not depend on
    !> the iterate. The values are the solver's, so they are scaled back into
    !> physical units on the way, the same way gather does for the iterate.
    subroutine update_history_fields(self, ierr)
        implicit none
        class(type_plexfem), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        Vec :: local_old, local_older
        PetscSection :: section, aux_section
        PetscScalar, pointer :: old_raw(:), older_raw(:), aux_raw(:)
        PetscInt :: point, p_start, p_end, offset, pressure_offset, num_dofs
        PetscInt :: aux_offset
        logical :: water

        water = hydraulic_is_active()

        call DMGetLocalSection(self%dm, section, ierr)
        if (ierr == 0) call DMGetLocalSection(self%aux_dm, aux_section, ierr)
        if (ierr /= 0) return
        call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr /= 0) return

        call DMGetLocalVector(self%dm, local_old, ierr)
        if (ierr == 0) call DMGetLocalVector(self%dm, local_older, ierr)
        if (ierr == 0) call gather(self%dm, self%old_solution, local_old, ierr)
        if (ierr == 0) call gather(self%dm, self%older_solution, local_older, ierr)
        if (ierr /= 0) return

        call VecGetArrayRead(local_old, old_raw, ierr)
        if (ierr == 0) call VecGetArrayRead(local_older, older_raw, ierr)
        if (ierr == 0) call VecGetArray(self%auxiliary, aux_raw, ierr)
        if (ierr /= 0) return

        do point = p_start, p_end - 1
            call PetscSectionGetFieldDof(section, point, 0_PETSC_INT_KIND, num_dofs, ierr)
            if (ierr /= 0 .or. num_dofs <= 0) cycle
            call PetscSectionGetFieldOffset(section, point, 0_PETSC_INT_KIND, offset, ierr)
            if (ierr /= 0) cycle

            call PetscSectionGetFieldOffset(aux_section, point, A_TEMPERATURE_OLD - 1, aux_offset, ierr)
            if (ierr == 0) aux_raw(aux_offset + 1) = old_raw(offset + 1)
            if (ierr == 0) call PetscSectionGetFieldOffset(aux_section, point, &
                                                          A_TEMPERATURE_OLDER - 1, aux_offset, ierr)
            if (ierr == 0) aux_raw(aux_offset + 1) = older_raw(offset + 1)
            if (ierr /= 0) cycle

            if (.not. water) cycle
            call PetscSectionGetFieldOffset(section, point, 1_PETSC_INT_KIND, pressure_offset, ierr)
            if (ierr /= 0) cycle
            call PetscSectionGetFieldOffset(aux_section, point, A_PRESSURE_OLD - 1, aux_offset, ierr)
            if (ierr == 0) aux_raw(aux_offset + 1) = old_raw(pressure_offset + 1)
            if (ierr == 0) call PetscSectionGetFieldOffset(aux_section, point, &
                                                          A_PRESSURE_OLDER - 1, aux_offset, ierr)
            if (ierr == 0) aux_raw(aux_offset + 1) = older_raw(pressure_offset + 1)
        end do

        call VecRestoreArray(self%auxiliary, aux_raw, ierr)
        call VecRestoreArrayRead(local_older, older_raw, ierr)
        call VecRestoreArrayRead(local_old, old_raw, ierr)
        call DMRestoreLocalVector(self%dm, local_older, ierr)
        call DMRestoreLocalVector(self%dm, local_old, ierr)
    end subroutine update_history_fields

    !> Refresh the element-representative transport coefficients from an iterate.
    !>
    !> For each cell the nodal coefficients are collected and reduced to one
    !> value: geometrically for the Darcy conductivity and the gravity flux,
    !> which are positive, and on the magnitude with the common sign restored for
    !> the thermo-osmotic coefficient, which is signed. A sign change inside an
    !> element has no geometric mean and is read as no coupling there rather than
    !> as an arbitrary branch, and a single zero node makes the element
    !> impermeable -- which is the point, because a frozen node does block flow
    !> through the element it belongs to.
    !>
    !> This is what the quadrature points then see, so an element straddling the
    !> freezing front passes what its frozen side allows rather than what its
    !> warm side would.
    !>
    !> Called from the residual and the Jacobian with the iterate each was handed,
    !> not with the stored solution: during a line search the trial point lives in
    !> a vector of its own.
    subroutine update_element_coefficients(dm, local_x, ierr)
        implicit none
        DM :: dm
        Vec :: local_x
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: section, aux_section
        PetscScalar, pointer :: aux_raw(:), x_raw(:)
        PetscInt, pointer :: points(:)
        PetscInt :: cell, c_start, c_end, point_count, p
        PetscInt :: material_offset, darcy_offset, coupling_offset, gravity_offset, potential_offset
        PetscInt :: offset, pressure_offset, num_dofs
        real(real64) :: darcy(SERIES_SAMPLES), coupling(SERIES_SAMPLES), gravity(SERIES_SAMPLES)
        real(real64) :: node_temperature(MAX_CELL_NODES), node_pressure(MAX_CELL_NODES)
        real(real64) :: node_darcy(MAX_CELL_NODES), node_coupling(MAX_CELL_NODES), node_gravity(MAX_CELL_NODES)
        real(real64) :: fraction
        integer(int32) :: warm, cold, samples
        real(real64) :: darcy_element, coupling_element, gravity_element
        real(real64) :: tensor(3, 3), coupling_tensor(3, 3), temperature, pressure, flux, potential
        integer(int32) :: material_id, vertices, v
        logical :: water

        water = hydraulic_is_active()
        if (.not. water) return

        call DMGetLocalSection(dm, section, ierr)
        if (ierr == 0) call DMGetLocalSection(solving%aux_dm, aux_section, ierr)
        if (ierr /= 0) return
        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return

        call VecGetArrayRead(local_x, x_raw, ierr)
        if (ierr == 0) call VecGetArray(solving%auxiliary, aux_raw, ierr)
        if (ierr /= 0) return

        do cell = c_start, c_end - 1
            call PetscSectionGetFieldOffset(aux_section, cell, A_MATERIAL_ID - 1, material_offset, ierr)
            if (ierr /= 0) cycle
            material_id = int(nint(real(aux_raw(material_offset + 1), real64)), int32)

            ! The vertices of this cell, found from its closure rather than from
            ! the order a closure of values would arrive in: the points carry
            ! their own section offsets, so nothing here has to assume how the
            ! two fields of a vertex are interleaved.
            nullify (points)
            call DMPlexGetTransitiveClosure(dm, cell, PETSC_TRUE, point_count, points, ierr)
            if (ierr /= 0) exit

            vertices = 0
            do p = 1, point_count
                call PetscSectionGetFieldDof(section, points(2 * p - 1), 0_PETSC_INT_KIND, num_dofs, ierr)
                if (ierr /= 0 .or. num_dofs <= 0) cycle
                if (vertices >= MAX_CELL_NODES) cycle
                call PetscSectionGetFieldOffset(section, points(2 * p - 1), 0_PETSC_INT_KIND, offset, ierr)
                if (ierr /= 0) cycle
                call PetscSectionGetFieldOffset(section, points(2 * p - 1), 1_PETSC_INT_KIND, pressure_offset, ierr)
                if (ierr /= 0) cycle

                vertices = vertices + 1
                node_temperature(vertices) = real(x_raw(offset + 1), real64)
                node_pressure(vertices) = real(x_raw(pressure_offset + 1), real64)

                ! The nodal liquid pressure, written where the auxiliary field
                ! carries it. One material per element is assumed: a node on a
                ! material interface would have one value per adjoining material
                ! and a single field cannot hold both, which is why more than one
                ! material is refused when the fields are built.
                call PetscSectionGetFieldOffset(aux_section, points(2 * p - 1), &
                                                A_POTENTIAL - 1, potential_offset, ierr)
                if (ierr /= 0) cycle
                call eval_liquid_pressure(node_temperature(vertices), node_pressure(vertices), &
                                          material_id, solving%constants, potential)
                aux_raw(potential_offset + 1) = real(potential, PETSC_REAL_KIND)
            end do

            call DMPlexRestoreTransitiveClosure(dm, cell, PETSC_TRUE, point_count, points, ierr)
            if (ierr /= 0) exit
            if (vertices < 1) cycle

            call PetscSectionGetFieldOffset(aux_section, cell, A_DARCY - 1, darcy_offset, ierr)
            if (ierr == 0) call PetscSectionGetFieldOffset(aux_section, cell, A_COUPLING - 1, coupling_offset, ierr)
            if (ierr == 0) call PetscSectionGetFieldOffset(aux_section, cell, A_GRAVITY - 1, gravity_offset, ierr)
            if (ierr /= 0) cycle

            ! The line the coefficients are sampled along runs between the
            ! vertices of least and greatest conductivity. Choosing it by
            ! conductivity rather than by temperature alone matters because
            ! pressure moves it too, and an isothermal cell would otherwise have
            ! both ends at the same vertex and its pressure spread discarded.
            do v = 1, vertices
                call sample_transport(material_id, node_temperature(v), node_pressure(v), &
                                      node_darcy(v), node_coupling(v), node_gravity(v))
            end do
            warm = 1
            cold = 1
            do v = 2, vertices
                if (node_darcy(v) > node_darcy(warm)) warm = v
                if (node_darcy(v) < node_darcy(cold)) cold = v
            end do

            ! The ends alone where the conductivity barely moves across the cell,
            ! which is almost every cell: a front is a few cells wide at a time.
            ! The interior is sampled only where the ends already differ by more
            ! than a factor, so most cells cost their vertices and no more. The
            ! two-sample form is the harmonic mean of the ends, which is the same
            ! series equivalent evaluated on one interval, so the switch changes
            ! how well the integral is resolved and not which integral it is.
            samples = 2
            darcy(1) = node_darcy(cold); coupling(1) = node_coupling(cold); gravity(1) = node_gravity(cold)
            darcy(2) = node_darcy(warm); coupling(2) = node_coupling(warm); gravity(2) = node_gravity(warm)
            if (darcy(1) > 0.0d0 .and. darcy(2) > 0.0d0) then
                if (max(darcy(1), darcy(2)) > SERIES_CONTRAST * min(darcy(1), darcy(2))) then
                    samples = SERIES_SAMPLES
                    do v = 1, samples
                        fraction = real(v - 1, real64) / real(samples - 1, real64)
                        temperature = (1.0d0 - fraction) * node_temperature(cold) &
                                      + fraction * node_temperature(warm)
                        pressure = (1.0d0 - fraction) * node_pressure(cold) &
                                   + fraction * node_pressure(warm)
                        call sample_transport(material_id, temperature, pressure, &
                                              darcy(v), coupling(v), gravity(v))
                    end do
                end if
            end if

            call series_transport(darcy, coupling, gravity, samples, &
                                  darcy_element, coupling_element, gravity_element)
            aux_raw(darcy_offset + 1) = real(darcy_element, PETSC_REAL_KIND)
            aux_raw(coupling_offset + 1) = real(coupling_element, PETSC_REAL_KIND)
            aux_raw(gravity_offset + 1) = real(gravity_element, PETSC_REAL_KIND)
        end do

        call VecRestoreArrayRead(local_x, x_raw, ierr)
        call VecRestoreArray(solving%auxiliary, aux_raw, ierr)
    end subroutine update_element_coefficients

    !> The transport coefficients at one state.
    subroutine sample_transport(material_id, temperature, pressure, darcy, coupling, gravity)
        implicit none
        integer(int32), intent(in) :: material_id
        real(real64), intent(in) :: temperature, pressure
        real(real64), intent(inout) :: darcy, coupling, gravity

        real(real64) :: tensor(3, 3), coupling_tensor(3, 3)

        tensor = 0.0d0
        coupling_tensor = 0.0d0
        call eval_water_diffusion(temperature, pressure, material_id, solving%constants, tensor)
        call eval_water_coupling(temperature, pressure, material_id, solving%constants, coupling_tensor)
        call eval_gravity_flux(temperature, pressure, material_id, solving%constants, gravity)
        darcy = tensor(1, 1)
        coupling = coupling_tensor(1, 1)
    end subroutine sample_transport

    !> The transport coefficients of a cell whose resistance is in series.
    !>
    !> For one-dimensional steady flow along the sampled line the flux is
    !> constant, so integrating
    !>
    !>   q = -D phi' - C T' - b
    !>
    !> across the cell gives the equivalents
    !>
    !>   D_eff = 1 / <1/D>,   C_eff = D_eff <C/D>,   b_eff = D_eff <b/D>
    !>
    !> and NOT the reciprocal mean of C or of b on their own: what a series
    !> arrangement shares is the flux, so every term is weighted by the local
    !> resistance. The averages are trapezoidal over the samples, which for two
    !> samples is exactly the harmonic mean of the ends.
    !>
    !> Why series rather than a geometric mean: measured on this material, across
    !> a cell straddling the freezing front the geometric mean returns twenty to
    !> thirty times the series value. It averages the logarithm, so a frozen
    !> tenth of a cell counts for as little as its thickness rather than as much
    !> as its resistance. Away from a front the conductivity is flat and the two
    !> agree to the digit, so the difference is confined to the cells the front
    !> is crossing.
    !>
    !> The approximation this makes is directional. The samples run along the
    !> cell's coldest-to-warmest edge, which is where the conductivity varies but
    !> is not in general the front normal, and the single number that comes out
    !> is then used for every direction. A series resistance is the right
    !> equivalent along the flow; across it the right one is an arithmetic mean.
    !> In a column freezing from one end the flow is along the gradient and the
    !> two nearly coincide, which is the case this is written for.
    subroutine series_transport(darcy, coupling, gravity, samples, &
                                darcy_element, coupling_element, gravity_element)
        implicit none
        real(real64), intent(in) :: darcy(:), coupling(:), gravity(:)
        integer(int32), intent(in) :: samples
        real(real64), intent(inout) :: darcy_element, coupling_element, gravity_element

        real(real64) :: resistance, coupling_share, gravity_share, weight
        integer(int32) :: i

        darcy_element = 0.0d0
        coupling_element = 0.0d0
        gravity_element = 0.0d0
        if (samples < 2) return
        if (any(darcy(1:samples) <= 0.0d0)) return

        resistance = 0.0d0
        coupling_share = 0.0d0
        gravity_share = 0.0d0
        do i = 1, samples
            ! Trapezoidal: the ends carry half a interval each.
            weight = merge(0.5d0, 1.0d0, i == 1 .or. i == samples)
            resistance = resistance + weight / darcy(i)
            coupling_share = coupling_share + weight * coupling(i) / darcy(i)
            gravity_share = gravity_share + weight * gravity(i) / darcy(i)
        end do
        if (resistance <= 0.0d0) return

        darcy_element = real(samples - 1, real64) / resistance
        coupling_element = darcy_element * coupling_share / real(samples - 1, real64)
        gravity_element = darcy_element * gravity_share / real(samples - 1, real64)
    end subroutine series_transport

    !> What each degree of freedom is measured in.
    !>
    !> One value per unknown per node, in the order the section lays them out, so
    !> that multiplying a solver vector by it gives physical values and dividing
    !> gives solver ones.
    subroutine build_unknown_scale(self, ierr)
        implicit none
        class(type_plexfem), intent(inout) :: self
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: point, p_start, p_end, offset, num_dofs
        logical :: water

        water = hydraulic_is_active()

        call DMCreateGlobalVector(self%dm, self%unknown_scale, ierr)
        if (ierr == 0) call DMCreateLocalVector(self%dm, self%unknown_scale_local, ierr)
        if (ierr /= 0) return
        ! One, not zero, where nothing writes: a scale of zero would erase the
        ! field it multiplies.
        call VecSet(self%unknown_scale_local, 1.0_PETSC_REAL_KIND, ierr)
        if (ierr /= 0) return

        call DMGetLocalSection(self%dm, section, ierr)
        if (ierr /= 0) return
        call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr == 0) call VecGetArray(self%unknown_scale_local, raw, ierr)
        if (ierr /= 0) return

        do point = p_start, p_end - 1
            call PetscSectionGetFieldDof(section, point, 0_PETSC_INT_KIND, num_dofs, ierr)
            if (ierr /= 0 .or. num_dofs <= 0) cycle
            call PetscSectionGetFieldOffset(section, point, 0_PETSC_INT_KIND, offset, ierr)
            if (ierr /= 0) cycle
            raw(offset + 1) = real(UNKNOWN_SCALE(1), PETSC_REAL_KIND)
            if (.not. water) cycle
            call PetscSectionGetFieldOffset(section, point, 1_PETSC_INT_KIND, offset, ierr)
            if (ierr /= 0) cycle
            raw(offset + 1) = real(UNKNOWN_SCALE(2), PETSC_REAL_KIND)
        end do

        call VecRestoreArray(self%unknown_scale_local, raw, ierr)
        if (ierr /= 0) return

        ! The owned part of the same thing, for scaling the matrix columns.
        call VecZeroEntries(self%unknown_scale, ierr)
        if (ierr == 0) call DMLocalToGlobalBegin(self%dm, self%unknown_scale_local, INSERT_VALUES, &
                                                 self%unknown_scale, ierr)
        if (ierr == 0) call DMLocalToGlobalEnd(self%dm, self%unknown_scale_local, INSERT_VALUES, &
                                               self%unknown_scale, ierr)
    end subroutine build_unknown_scale

    !> Open the solver history and write its schema.
    !>
    !> One record per time-step attempt, in the format the validation reads. The
    !> columns it does not evaluate are written as -1, which the schema line
    !> defines as "not evaluated", so the file stays readable by one parser
    !> whichever path produced it.
    subroutine open_history(self)
        implicit none
        class(type_plexfem), intent(inout) :: self

        character(:), allocatable :: project_path, file_name
        PetscMPIInt :: rank
        PetscErrorCode :: ierr

        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
        if (rank /= 0) return

        call get_env_string("FTCMS_PROJECT_PATH", project_path)
        call modify_path_format(project_path)
        file_name = strip(project_path)//"Output/solver_history.log"

        open (newunit=history_unit, file=file_name, status="replace", action="write")
        write (history_unit, '(a)') "# FTCMS solver history schema=4: one record per time-step attempt"
        write (history_unit, '(a)') "# attempt accepted_step time_start_s time_trial_s time_accepted_s "// &
            "dt_used_s dt_next_s accepted status inner_iter ats_iter loc_T loc_H bal_T bal_H dq_eff "// &
            "du_T_max du_p_max omega lte_rel"
        write (history_unit, '(a)') "# loc_*/bal_*/dq_eff are the acceptance gates of the last iterate: "// &
            "each passes at <= 1, -1 means not evaluated. status is how the attempt ended, "// &
            "not which gate refused it"
        flush (history_unit)
    end subroutine open_history

    !> Record one attempt.
    subroutine log_attempt(self, attempt, step_number, time_start, dt_used, dt_next, &
                           accepted, status, iterations, error)
        implicit none
        class(type_plexfem), intent(in) :: self
        integer(int32), intent(in) :: attempt, step_number
        real(real64), intent(in) :: time_start, dt_used, dt_next
        logical, intent(in) :: accepted
        character(*), intent(in) :: status
        integer(int32), intent(in) :: iterations
        real(real64), intent(in) :: error

        real(real64) :: time_accepted
        PetscMPIInt :: rank
        PetscErrorCode :: ierr

        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
        if (rank /= 0) return

        time_accepted = time_start
        if (accepted) time_accepted = time_start + dt_used

        write (history_unit, '(i0,1x,i0,3(1x,es16.9),2(1x,es13.6),1x,i0,1x,a,2(1x,i0),8(1x,f5.1),1x,es13.6)') &
            attempt, merge(step_number + 1, step_number, accepted), &
            time_start, time_start + dt_used, time_accepted, &
            dt_used, dt_next, merge(1, 0, accepted), trim(status), &
            iterations, iterations, &
            -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
            error
        flush (history_unit)
    end subroutine log_attempt

    !> Open the one file this run writes.
    !>
    !> The mesh goes in first and once; every output time after that is a step
    !> inside the same file, so a run is a single self-contained HDF5 rather
    !> than a directory of snapshots that only mean something together.
    subroutine open_output(self)
        implicit none
        class(type_plexfem), intent(inout) :: self

        PetscErrorCode :: ierr
        character(:), allocatable :: project_path, file_name

        call get_env_string("FTCMS_PROJECT_PATH", project_path)
        call modify_path_format(project_path)
        file_name = strip(project_path)//"Output/Files/PlexFEM.h5"

        call PetscViewerHDF5Open(PETSC_COMM_WORLD, file_name, FILE_MODE_WRITE, self%output_viewer, ierr)
        if (ierr /= 0) error stop "PETSc could not open the output file"
        self%output_ready = .true.

        ! XDMF layout, so the file describes its own geometry and topology and
        ! a reader needs nothing beside it.
        !
        ! The mesh written is the scalar output clone, not the solver's DM: the
        ! HDF5 writer names one dataset per vector, so a two-field vector would
        ! arrive as one array of two components and the quantities would have to
        ! be told apart by position. Written one at a time on a DM carrying a
        ! single value per node, each is a dataset under its own name.
        call self%input%geometry%mesh%get_output_dm(self%output_dm)
        call DMCreateGlobalVector(self%output_dm, self%output_field, ierr)
        if (ierr /= 0) error stop "PETSc could not create the output field"

        call PetscViewerPushFormat(self%output_viewer, PETSC_VIEWER_HDF5_XDMF, ierr)
        if (ierr == 0) call DMView(self%output_dm, self%output_viewer, ierr)
        if (ierr == 0) call PetscViewerHDF5PushTimestepping(self%output_viewer, ierr)
        if (ierr /= 0) error stop "PETSc could not write the mesh to the output file"
    end subroutine open_output

    !> Write one nodal quantity under its own name.
    subroutine write_named_field(self, name, quantity, ierr)
        implicit none
        class(type_plexfem), intent(inout) :: self
        character(*), intent(in) :: name
        !> Which quantity to gather: see nodal_quantity.
        integer(int32), intent(in) :: quantity
        PetscErrorCode, intent(inout) :: ierr

        Vec :: local_field

        call DMGetLocalVector(self%output_dm, local_field, ierr)
        if (ierr /= 0) return
        call VecZeroEntries(local_field, ierr)
        if (ierr == 0) call nodal_quantity(self, quantity, local_field, ierr)
        if (ierr /= 0) return

        call VecZeroEntries(self%output_field, ierr)
        if (ierr == 0) call DMLocalToGlobalBegin(self%output_dm, local_field, INSERT_VALUES, &
                                                 self%output_field, ierr)
        if (ierr == 0) call DMLocalToGlobalEnd(self%output_dm, local_field, INSERT_VALUES, &
                                               self%output_field, ierr)
        if (ierr == 0) call DMRestoreLocalVector(self%output_dm, local_field, ierr)
        if (ierr /= 0) return

        call PetscObjectSetName(PetscObjectCast(self%output_field), name, ierr)
        if (ierr == 0) call VecView(self%output_field, self%output_viewer, ierr)
    end subroutine write_named_field

    !> Gather one nodal quantity of the current solution onto the output DM.
    !>
    !> The two DMs are clones of the same mesh, so a point means the same node
    !> in both; only the offset each section gives it differs. Where materials
    !> meet, the value written is their volume-weighted mean, which is the same
    !> average the storage term integrates.
    subroutine nodal_quantity(self, quantity, local_field, ierr)
        implicit none
        class(type_plexfem), intent(inout) :: self
        integer(int32), intent(in) :: quantity
        Vec :: local_field
        PetscErrorCode, intent(inout) :: ierr

        Vec :: local_x, local_weight
        PetscSection :: section, out_section
        PetscScalar, pointer :: x_raw(:), field_raw(:), weight_raw(:), volume_raw(:)
        PetscInt :: point, p_start, p_end, offset, out_offset, num_dofs, pressure_offset
        real(real64) :: temperature, pressure, weight, liquid, ice, value
        integer(int32) :: m, material_id
        logical :: water

        water = hydraulic_is_active()

        call DMGetLocalSection(self%dm, section, ierr)
        if (ierr == 0) call DMGetLocalSection(self%output_dm, out_section, ierr)
        if (ierr /= 0) return
        call PetscSectionGetChart(section, p_start, p_end, ierr)
        if (ierr /= 0) return

        call DMGetLocalVector(self%dm, local_x, ierr)
        if (ierr == 0) call gather(self%dm, self%solution, local_x, ierr)
        if (ierr /= 0) return
        ! The weights the mean divides by, accumulated alongside it.
        call DMGetLocalVector(self%output_dm, local_weight, ierr)
        if (ierr == 0) call VecZeroEntries(local_weight, ierr)
        if (ierr /= 0) return

        call VecGetArrayRead(local_x, x_raw, ierr)
        if (ierr == 0) call VecGetArray(local_field, field_raw, ierr)
        if (ierr == 0) call VecGetArray(local_weight, weight_raw, ierr)
        if (ierr /= 0) return

        do m = 1, self%num_materials
            material_id = self%material_id(m)
            call VecGetArrayRead(self%nodal_volume(m), volume_raw, ierr)
            if (ierr /= 0) return

            do point = p_start, p_end - 1
                call PetscSectionGetFieldDof(section, point, 0_PETSC_INT_KIND, num_dofs, ierr)
                if (ierr /= 0 .or. num_dofs <= 0) cycle
                call PetscSectionGetFieldOffset(section, point, 0_PETSC_INT_KIND, offset, ierr)
                if (ierr /= 0) cycle
                weight = real(volume_raw(offset + 1), real64)
                if (weight == 0.0d0) cycle
                call PetscSectionGetOffset(out_section, point, out_offset, ierr)
                if (ierr /= 0) cycle

                temperature = real(x_raw(offset + 1), real64)
                pressure = real(self%constants(C_PRESSURE), real64)
                if (water) then
                    call PetscSectionGetFieldOffset(section, point, 1_PETSC_INT_KIND, pressure_offset, ierr)
                    if (ierr /= 0) cycle
                    pressure = real(x_raw(pressure_offset + 1), real64)
                end if

                select case (quantity)
                case (NODAL_TEMPERATURE)
                    value = temperature
                case (NODAL_PRESSURE)
                    value = pressure
                case default
                    call eval_phase_contents(temperature, pressure, material_id, self%constants, liquid, ice)
                    value = merge(liquid, ice, quantity == NODAL_WATER_CONTENT)
                end select

                field_raw(out_offset + 1) = field_raw(out_offset + 1) + weight * value
                weight_raw(out_offset + 1) = weight_raw(out_offset + 1) + weight
            end do

            call VecRestoreArrayRead(self%nodal_volume(m), volume_raw, ierr)
        end do

        do point = 1, size(field_raw)
            if (weight_raw(point) /= 0.0_PETSC_REAL_KIND) then
                field_raw(point) = field_raw(point) / weight_raw(point)
            end if
        end do

        call VecRestoreArray(local_weight, weight_raw, ierr)
        call VecRestoreArray(local_field, field_raw, ierr)
        call VecRestoreArrayRead(local_x, x_raw, ierr)
        call DMRestoreLocalVector(self%output_dm, local_weight, ierr)
        call DMRestoreLocalVector(self%dm, local_x, ierr)
    end subroutine nodal_quantity

    !> Write the temperature, if this step reached the output interval.
    subroutine write_field(self, time, output_number, next_output, output_interval)
        implicit none
        class(type_plexfem), intent(inout) :: self
        real(real64), intent(in) :: time
        integer(int32), intent(inout) :: output_number
        real(real64), intent(inout) :: next_output
        real(real64), intent(in) :: output_interval

        PetscErrorCode :: ierr

        ! A step lands on or past the mark; the next mark is counted from the
        ! interval, not from this time, so rounding cannot make them drift.
        if (output_interval <= 0.0d0) return
        if (time < next_output - 1.0d-9) return
        next_output = next_output + output_interval

        call PetscViewerHDF5SetTimestep(self%output_viewer, int(output_number, PETSC_INT_KIND), ierr)
        if (ierr /= 0) error stop "PETSc could not set the output time step"

        ! The liquid and the ice go out separately rather than as one total: the
        ! total is theta_w + (rho_i/rho_w) theta_i, and forming it here would
        ! decide the density ratio on the reader's behalf.
        call write_named_field(self, "temperature", NODAL_TEMPERATURE, ierr)
        if (ierr == 0 .and. hydraulic_is_active()) then
            call write_named_field(self, "pressure", NODAL_PRESSURE, ierr)
            if (ierr == 0) call write_named_field(self, "WaterContent", NODAL_WATER_CONTENT, ierr)
            if (ierr == 0) call write_named_field(self, "IceContent", NODAL_ICE_CONTENT, ierr)
        end if
        if (ierr /= 0) error stop "PETSc could not write the output fields"

        output_number = output_number + 1
    end subroutine write_field

    !> Tell the kernels how long this step is, and set the BDF coefficients.
    !>
    !> The coefficients are the derivative of the Lagrange polynomial through the
    !> states at t_n, t_{n-1} and, at second order, t_{n-2}, evaluated at t_n:
    !>
    !>   tau_1 = h,  tau_2 = h + h_prev
    !>   c_0 = 1/tau_1 + 1/tau_2
    !>   c_1 = -tau_2/(tau_1 (tau_2 - tau_1))
    !>   c_2 =  tau_1/(tau_2 (tau_2 - tau_1))
    !>
    !> which at equal steps is (3/2, -2, 1/2)/h. The current step's own length
    !> enters tau_1 directly, so a step that was cut takes effect immediately
    !> rather than lagging one step behind the history.
    subroutine set_time_step(self, step, previous_step, order)
        implicit none
        class(type_plexfem), intent(inout) :: self
        real(real64), intent(in) :: step
        real(real64), intent(in) :: previous_step
        integer(int32), intent(in) :: order

        PetscErrorCode :: ierr
        real(real64) :: tau_1, tau_2

        self%bdf_order = 1
        if (order >= 2 .and. previous_step > 0.0d0) self%bdf_order = 2

        self%bdf_coeffs(:) = 0.0d0
        tau_1 = step
        if (self%bdf_order == 1) then
            self%bdf_coeffs(0) = 1.0d0 / tau_1
            self%bdf_coeffs(1) = -1.0d0 / tau_1
        else
            tau_2 = step + previous_step
            self%bdf_coeffs(0) = 1.0d0 / tau_1 + 1.0d0 / tau_2
            self%bdf_coeffs(1) = -tau_2 / (tau_1 * (tau_2 - tau_1))
            self%bdf_coeffs(2) = tau_1 / (tau_2 * (tau_2 - tau_1))
        end if

        self%constants(C_DT) = real(step, PETSC_REAL_KIND)
        self%constants(C_BDF_ORDER) = real(self%bdf_order, PETSC_REAL_KIND)
        self%constants(C_BDF_0) = real(self%bdf_coeffs(0), PETSC_REAL_KIND)
        self%constants(C_BDF_0 + 1) = real(self%bdf_coeffs(1), PETSC_REAL_KIND)
        self%constants(C_BDF_0 + 2) = real(self%bdf_coeffs(2), PETSC_REAL_KIND)
        call PetscDSSetConstants(self%ds, int(PLEXFEM_CONSTANTS, PETSC_INT_KIND), self%constants, ierr)
        if (ierr /= 0) error stop "PETSc rejected the step length"
    end subroutine set_time_step

    module subroutine destroy_plexfem(self)
        implicit none
        class(type_plexfem), intent(inout) :: self

        PetscErrorCode :: ierr
        integer(int32) :: i

        ! Before the model itself goes, so no kernel can be left pointing at it.
        call detach_thermal_model()
        call self%thermal%destroy()
        call self%hydraulic%destroy()

        if (self%snes_ready) then
            call SNESDestroy(self%snes, ierr)
            self%snes_ready = .false.
        end if
        do i = 1, self%num_boundary_aux
            call VecDestroy(self%boundary_aux(i), ierr)
        end do
        self%num_boundary_aux = 0

        if (self%dm_ready) then
            call VecDestroy(self%auxiliary, ierr)
            call VecDestroy(self%old_solution, ierr)
            call VecDestroy(self%older_solution, ierr)
            call VecDestroy(self%row_scale, ierr)
            call VecDestroy(self%unknown_scale, ierr)
            call VecDestroy(self%unknown_scale_local, ierr)
            do i = 1, self%num_materials
                call VecDestroy(self%nodal_volume(i), ierr)
            end do
            self%num_materials = 0
            call VecDestroy(self%solution, ierr)
            call DMDestroy(self%aux_dm, ierr)
            call DMDestroy(self%dm, ierr)
            self%dm_ready = .false.
        end if

        ! The mesh reader is what initialised PETSc, lazily; closing it here is
        ! what releases the objects it still owns and what makes -log_view
        ! report anything at all.
        call self%input%geometry%mesh%destroy()
        call PetscFinalize(ierr)
    end subroutine destroy_plexfem

end submodule plexfem_run
