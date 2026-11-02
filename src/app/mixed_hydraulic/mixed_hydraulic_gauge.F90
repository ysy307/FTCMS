#include <petsc/finclude/petscdmplex.h>

!> The pressure gauge of a steady problem whose whole boundary prescribes flux.
!>
!> When the storage vanishes and every boundary facet constrains its normal
!> flux, the discrete operator has the constant pressure in its kernel and the
!> same vector in its left kernel. The first makes the pressure undetermined up
!> to a constant, the second makes the data solvable only if the prescribed
!> boundary flux and the source balance. Neither is a property of this
!> implementation: both follow from \(\sum_K B_{Ki}=\int_{\partial\Omega}
!> \boldsymbol\phi_i\cdot\mathbf n=0\) for every unconstrained facet function.
!>
!> The gauge fixes the constant by adding one term to the residual and the same
!> term to the Jacobian,
!>
!> \[
!>   F_g(x)=F(x)+\gamma\Big(\frac{1}{V_\Omega}\sum_K V_K\hat p_K\Big)e_{K_0},
!>   \qquad J_g=J+\gamma\,e_{K_0}v^{T},\quad v_K=V_K/V_\Omega .
!> \]
!>
!> Summing every conservation row annihilates \(J\) and leaves \(\gamma\) times
!> the volume-weighted mean pressure, so \(F_g(x)=0\) holds exactly when
!> \(F(x)=0\) and that mean is zero: the added term is zero at the solution and
!> every cell conservation equation still holds. This is an exact rank-one
!> deflation -- the bordered Lagrange system with its multiplier eliminated --
!> not a penalty, a projection or a perturbed diagonal. Its point is that
!> \(J_g\) is nonsingular, which is what a direct factorisation needs; a null
!> space attached to a singular matrix is not.
!>
!> \(\gamma>0\) changes the conditioning and nothing else. It is measured from
!> the assembled row it is added to, so the gauge equation carries the weight of
!> its own physical row.
module app_mixed_hydraulic_gauge
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_funptr, c_funloc, &
                                           c_null_ptr, c_associated
    use :: petscdmplex
    use :: core_interop_petscds, only:petscds_set_snes_local_residual, petscds_set_snes_local_jacobian, &
                                      petscds_snes_local_residual_fem, petscds_snes_local_jacobian_fem
    use :: app_mixed_hydraulic_balance, only:MAX_CELL_FACETS, reconstruct_cell_facet_flux, &
                                             cell_is_owned, cell_volume
    use :: app_mixed_hydraulic_physics, only:A_STORAGE, A_SOURCE, MIXED_FIELD_FLUX, MIXED_FIELD_PRESSURE
    implicit none
    private

    public :: type_mixed_gauge_report
    public :: measure_gauge_state
    public :: install_pressure_gauge
    public :: clear_pressure_gauge
    public :: gauge_compatibility_defect
    public :: gauge_compatibility_scale
    public :: gauge_data_is_compatible
    public :: GAUGE_COMPATIBILITY_TOLERANCE

    !> How far the boundary and source integrals of a steady all-flux problem
    !> may disagree, relative to the larger of their absolute totals.
    !>
    !> The two sides are sums over facets and cells of numbers that are exact in
    !> the data, so the only difference a compatible problem can show is the
    !> rounding of those sums: a few units in the last place times the count.
    !> Anything above this is a statement about the data, not the arithmetic.
    real(real64), parameter :: GAUGE_COMPATIBILITY_TOLERANCE = 1.0e-10_real64

    !> What one measurement of a steady all-flux candidate found, in physical
    !> units per unit out-of-plane thickness.
    type :: type_mixed_gauge_report
        !> True when the constant pressure really is in the kernel: every
        !> physical boundary facet constrains its flux and no cell stores water.
        logical :: eligible = .false.
        !> Prescribed outward flux through the whole boundary, m2/s.
        real(real64) :: prescribed_outflow = 0.0_real64
        !> Sum of the absolute prescribed facet fluxes, m2/s.
        real(real64) :: prescribed_magnitude = 0.0_real64
        !> Integrated source over the domain, m2/s.
        real(real64) :: source_total = 0.0_real64
        !> Sum of the absolute cell sources, m2/s.
        real(real64) :: source_magnitude = 0.0_real64
        !> Sum of the cell storage magnitudes, m2/Pa; zero makes the run steady.
        !>
        !> A magnitude of each cell's storage law rather than of its tangent at
        !> one state, so a nonlinear curve whose capacity happens to vanish at
        !> the accepted pressure still keeps the gauge off.
        real(real64) :: storage_magnitude = 0.0_real64
        !> Area of the domain, m2.
        real(real64) :: total_volume = 0.0_real64
        !> Area of the largest cell, m2.
        real(real64) :: largest_cell_volume = 0.0_real64
        !> How many facets have exactly one adjacent cell over all ranks.
        integer(int32) :: boundary_facets = 0
        !> How many of those do not constrain their normal flux.
        integer(int32) :: unconstrained_boundary_facets = 0
    end type type_mixed_gauge_report

    !> The gauge the registered hooks read.
    !>
    !> One instance, because the hooks reach it by name rather than through a
    !> context pointer. `dm` records which mesh it was installed for, so a
    !> second gauged solver built while the first is alive is refused instead of
    !> silently answering with the wrong geometry.
    type :: type_gauge_state
        logical :: active = .false.
        type(c_ptr) :: dm = c_null_ptr
        !> Global index of the conservation row that carries the constraint.
        PetscInt :: row = -1_PETSC_INT_KIND
        !> Offset of that cell's pressure in a local vector; -1 off its rank.
        PetscInt :: local_offset = -1_PETSC_INT_KIND
        !> \(\gamma\), the weight of the gauge equation.
        real(real64) :: scale = 0.0_real64
        !> \(V_\Omega\), m2.
        real(real64) :: total_volume = 0.0_real64
    end type type_gauge_state

    type(type_gauge_state), save :: gauge

contains

    !> Measure everything that decides whether a run needs a pressure gauge.
    !>
    !> \[
    !>   \sum_{f\subset\partial\Omega}Q_f^{\mathrm{presc}},\qquad
    !>   \sum_K s_KV_K,\qquad \sum_K|c_K|V_K,\qquad V_\Omega=\sum_KV_K .
    !> \]
    !>
    !> A facet is on the physical boundary when exactly one cell is adjacent to
    !> it over all ranks, which is counted from the owned cells and reduced; a
    !> partition interface has two. Assumes `local_x` carries the prescribed
    !> boundary values and nothing else, so the facet moments it reads back are
    !> the ones the solver will impose rather than an assumed formula.
    !>
    !> Guarantees identical results on every rank: every sum is reduced before
    !> it is read, so a decision taken from this report is taken collectively.
    !> Arithmetic is \(O(N)\) in local cells, memory \(O(P)\) in the mesh's
    !> global point count.
    !>
    !> Fails by returning the first PETSc error with `eligible` false.
    subroutine measure_gauge_state(dm, local_x, aux_dm, auxiliary, flux_scale, report, ierr)
        implicit none
        !> The solve DM, carrying the flux and pressure fields.
        DM, intent(in) :: dm
        !> A local vector holding the prescribed boundary values.
        Vec, intent(in) :: local_x
        !> The DM the cell coefficients live on.
        DM, intent(in) :: aux_dm
        !> The cell coefficients.
        Vec, intent(in) :: auxiliary
        !> The flux unit of the unknown, m/s.
        real(real64), intent(in) :: flux_scale
        !> The measured state, reduced over the communicator.
        type(type_mixed_gauge_report), intent(inout) :: report
        PetscErrorCode, intent(inout) :: ierr

        IS :: numbering
        PetscInt, pointer :: numbers(:)
        PetscSection :: section, aux_section
        PetscScalar, pointer :: aux_raw(:)
        PetscInt :: cell, c_start, c_end, p_start, p_end, point
        PetscInt :: facets(MAX_CELL_FACETS)
        PetscInt :: offset, num_dofs, num_constrained
        real(real64) :: moments(MAX_CELL_FACETS)
        real(real64), allocatable :: facet_count(:), facet_count_all(:)
        real(real64) :: sums(8), reduced(8), largest_send, largest
        real(real64) :: volume, source, storage
        integer(int32) :: num_facets, f, index, total_points, mpi_err

        report = type_mixed_gauge_report()

        call DMPlexGetChart(dm, p_start, p_end, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr == 0) call DMGetLocalSection(dm, section, ierr)
        if (ierr == 0) call DMGetLocalSection(aux_dm, aux_section, ierr)
        if (ierr == 0) call DMPlexCreatePointNumbering(dm, numbering, ierr)
        if (ierr /= 0) return
        call ISGetIndices(numbering, numbers, ierr)
        if (ierr /= 0) return

        ! The same indexing the facet-continuity audit uses: a point this rank
        ! does not own carries minus its global number minus one, so both sides
        ! of a shared facet land on one index.
        largest_send = 0.0_real64
        do point = p_start, p_end - 1
            largest_send = max(largest_send, real(global_index_of(numbers(point - p_start + 1)), real64))
        end do
        call MPI_Allreduce(largest_send, largest, 1, MPI_DOUBLE_PRECISION, MPI_MAX, PETSC_COMM_WORLD, mpi_err)
        if (mpi_err /= 0) then
            call ISRestoreIndices(numbering, numbers, ierr)
            call ISDestroy(numbering, ierr)
            ierr = 1
            return
        end if
        total_points = int(largest, int32) + 1

        allocate (facet_count(total_points), facet_count_all(total_points))
        facet_count = 0.0_real64

        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(dm, cell)) cycle
            call reconstruct_cell_facet_flux(dm, local_x, flux_scale, cell, facets, moments, num_facets, ierr)
            if (ierr /= 0) exit
            do f = 1, num_facets
                index = global_index_of(numbers(facets(f) - p_start + 1)) + 1
                facet_count(index) = facet_count(index) + 1.0_real64
            end do
        end do
        if (ierr /= 0) then
            call ISRestoreIndices(numbering, numbers, ierr)
            call ISDestroy(numbering, ierr)
            deallocate (facet_count, facet_count_all)
            return
        end if

        call MPI_Allreduce(facet_count, facet_count_all, total_points, MPI_DOUBLE_PRECISION, MPI_SUM, &
                           PETSC_COMM_WORLD, mpi_err)
        if (mpi_err /= 0) then
            call ISRestoreIndices(numbering, numbers, ierr)
            call ISDestroy(numbering, ierr)
            deallocate (facet_count, facet_count_all)
            ierr = 1
            return
        end if

        sums = 0.0_real64
        largest_send = 0.0_real64

        call VecGetArrayRead(auxiliary, aux_raw, ierr)
        if (ierr == 0) then
            do cell = c_start, c_end - 1
                if (.not. cell_is_owned(dm, cell)) cycle

                call cell_volume(dm, cell, volume, ierr)
                if (ierr /= 0) exit
                sums(1) = sums(1) + volume
                largest_send = max(largest_send, volume)

                source = 0.0_real64
                storage = 0.0_real64
                call PetscSectionGetFieldDof(aux_section, cell, int(A_SOURCE - 1, PETSC_INT_KIND), num_dofs, ierr)
                if (ierr == 0 .and. num_dofs > 0) then
                    call PetscSectionGetFieldOffset(aux_section, cell, int(A_SOURCE - 1, PETSC_INT_KIND), offset, ierr)
                    if (ierr == 0) source = real(aux_raw(offset + 1), real64)
                    call PetscSectionGetFieldOffset(aux_section, cell, int(A_STORAGE - 1, PETSC_INT_KIND), offset, ierr)
                    if (ierr == 0) storage = real(aux_raw(offset + 1), real64)
                end if
                sums(2) = sums(2) + source*volume
                sums(3) = sums(3) + abs(source)*volume
                sums(4) = sums(4) + abs(storage)*volume

                call reconstruct_cell_facet_flux(dm, local_x, flux_scale, cell, facets, moments, num_facets, ierr)
                if (ierr /= 0) exit
                do f = 1, num_facets
                    index = global_index_of(numbers(facets(f) - p_start + 1)) + 1
                    if (nint(facet_count_all(index)) /= 1) cycle
                    sums(5) = sums(5) + moments(f)
                    sums(6) = sums(6) + abs(moments(f))
                    sums(7) = sums(7) + 1.0_real64

                    ! A prescribed normal flux is an essential condition, so the
                    ! facet's degrees of freedom are constrained out of the
                    ! global system. That, and not the label, is what makes the
                    ! constant pressure a kernel vector.
                    call PetscSectionGetFieldDof(section, facets(f), MIXED_FIELD_FLUX, num_dofs, ierr)
                    if (ierr /= 0) exit
                    call PetscSectionGetFieldConstraintDof(section, facets(f), MIXED_FIELD_FLUX, &
                                                           num_constrained, ierr)
                    if (ierr /= 0) exit
                    if (num_constrained < num_dofs) sums(8) = sums(8) + 1.0_real64
                end do
                if (ierr /= 0) exit
            end do
            call VecRestoreArrayRead(auxiliary, aux_raw, ierr)
        end if

        call ISRestoreIndices(numbering, numbers, ierr)
        call ISDestroy(numbering, ierr)
        deallocate (facet_count, facet_count_all)
        if (ierr /= 0) return

        call MPI_Allreduce(sums, reduced, 8, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, mpi_err)
        if (mpi_err == 0) then
            call MPI_Allreduce(largest_send, largest, 1, MPI_DOUBLE_PRECISION, MPI_MAX, PETSC_COMM_WORLD, mpi_err)
        end if
        if (mpi_err /= 0) then
            ierr = 1
            return
        end if

        report%total_volume = reduced(1)
        report%source_total = reduced(2)
        report%source_magnitude = reduced(3)
        report%storage_magnitude = reduced(4)
        report%prescribed_outflow = reduced(5)
        report%prescribed_magnitude = reduced(6)
        report%boundary_facets = nint(reduced(7))
        report%unconstrained_boundary_facets = nint(reduced(8))
        report%largest_cell_volume = largest

        report%eligible = (report%unconstrained_boundary_facets == 0) &
                          .and. (report%boundary_facets > 0) &
                          .and. (report%storage_magnitude <= 0.0_real64) &
                          .and. (report%total_volume > 0.0_real64)
    end subroutine measure_gauge_state

    !> The imbalance between the prescribed boundary flux and the source, m2/s.
    !>
    !> \[
    !>   d=\sum_{f\subset\partial\Omega}Q_f^{\mathrm{presc}}-\sum_K s_KV_K .
    !> \]
    !>
    !> This is the whole solvability condition of the steady all-flux problem:
    !> it is the left kernel vector applied to the data. Exact for the measured
    !> report; \(O(1)\).
    pure function gauge_compatibility_defect(report) result(defect)
        implicit none
        type(type_mixed_gauge_report), intent(in) :: report
        real(real64) :: defect

        defect = report%prescribed_outflow - report%source_total
    end function gauge_compatibility_defect

    !> The size the imbalance is measured against, m2/s.
    pure function gauge_compatibility_scale(report) result(scale)
        implicit none
        type(type_mixed_gauge_report), intent(in) :: report
        real(real64) :: scale

        scale = max(report%prescribed_magnitude, report%source_magnitude)
    end function gauge_compatibility_scale

    !> Whether the data of a steady all-flux run can be solved at all.
    !>
    !> Data with no flux and no source are compatible: both sides are zero.
    !> \(O(1)\); no failure mode of its own.
    pure function gauge_data_is_compatible(report) result(ok)
        implicit none
        type(type_mixed_gauge_report), intent(in) :: report
        logical :: ok

        real(real64) :: scale

        scale = gauge_compatibility_scale(report)
        if (scale <= 0.0_real64) then
            ok = (gauge_compatibility_defect(report) == 0.0_real64)
        else
            ok = abs(gauge_compatibility_defect(report)) <= GAUGE_COMPATIBILITY_TOLERANCE*scale
        end if
    end function gauge_data_is_compatible

    !> Install the volume-weighted mean-pressure gauge on this DM.
    !>
    !> \[
    !>   \frac{1}{V_\Omega}\sum_K V_K\hat p_K=0 ,
    !> \]
    !>
    !> carried by the conservation row of the cell with the lowest global index
    !> and weighted by \(\gamma=\|J^{T}e_{K_0}\|_\infty V_\Omega/\max_KV_K\),
    !> measured from a probe assembly of the ungauged Jacobian.
    !>
    !> Assumes `report` came from `measure_gauge_state` on this DM and is
    !> eligible, and that `DMPlexSetSNESLocalFEM` has already been called: the
    !> hooks replace what it registered. Guarantees the solved state satisfies
    !> both the gauge and every cell conservation equation, to the accuracy of
    !> the linear solve. Setup costs one Jacobian assembly and one matrix; each
    !> later residual costs one reduction and each Jacobian one dense row of
    !> \(N\) entries, which is this gauge's price on a large mesh.
    !>
    !> Fails by `error stop` when a cell does not carry exactly one global
    !> pressure degree of freedom, or when a gauged solver is already installed
    !> on a different mesh.
    subroutine install_pressure_gauge(dm, fallback_scale, report, ierr)
        implicit none
        !> The solve DM.
        DM, intent(in) :: dm
        !> The scale to weight the gauge with if the probe row is empty, m.
        real(real64), intent(in) :: fallback_scale
        !> The eligible measurement of this DM.
        type(type_mixed_gauge_report), intent(in) :: report
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: section
        PetscInt :: cell, c_start, c_end, global_start, global_end, offset, num_dofs
        PetscInt :: lowest_local, lowest
        real(real64) :: send(1), reduced(1), row_scale
        integer(int32) :: mpi_err
        integer(c_int) :: cerr

        if (gauge%active .and. .not. c_associated(gauge%dm, transfer(dm%v, gauge%dm))) then
            error stop "A second gauged mixed hydraulic solver cannot be installed while the first is alive"
        end if

        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr == 0) call DMGetLocalSection(dm, section, ierr)
        if (ierr /= 0) return

        ! The lowest global pressure index in the mesh, so the carrying row is
        ! the same one whatever the partition is.
        lowest_local = huge(lowest_local)
        do cell = c_start, c_end - 1
            call DMPlexGetPointGlobal(dm, cell, global_start, global_end, ierr)
            if (ierr /= 0) return
            if (global_start < 0) cycle
            if (global_end - global_start /= 1_PETSC_INT_KIND) then
                error stop "A cell of the mixed hydraulic mesh does not carry exactly one pressure unknown"
            end if
            lowest_local = min(lowest_local, global_start)
        end do

        send(1) = real(lowest_local, real64)
        call MPI_Allreduce(send, reduced, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PETSC_COMM_WORLD, mpi_err)
        if (mpi_err /= 0) then
            ierr = 1
            return
        end if
        lowest = int(reduced(1), PETSC_INT_KIND)

        gauge%row = lowest
        gauge%local_offset = -1_PETSC_INT_KIND
        do cell = c_start, c_end - 1
            call DMPlexGetPointGlobal(dm, cell, global_start, global_end, ierr)
            if (ierr /= 0) return
            if (global_start /= lowest) cycle
            call PetscSectionGetFieldDof(section, cell, MIXED_FIELD_PRESSURE, num_dofs, ierr)
            if (ierr == 0 .and. num_dofs > 0) then
                call PetscSectionGetFieldOffset(section, cell, MIXED_FIELD_PRESSURE, offset, ierr)
                if (ierr == 0) gauge%local_offset = offset
            end if
            if (ierr /= 0) return
        end do

        ! The weight of the gauge equation, read off the row it is added to:
        ! J^T e_{K0} is that row, so its largest entry is the size the
        ! constraint has to match to be seen by the factorisation.
        call measure_row_scale(dm, lowest, row_scale, ierr)
        if (ierr /= 0) return
        if (row_scale <= 0.0_real64) row_scale = fallback_scale
        if (report%largest_cell_volume <= 0.0_real64) then
            error stop "The mixed hydraulic mesh reports no cell area for the pressure gauge"
        end if

        gauge%scale = row_scale*report%total_volume/report%largest_cell_volume
        gauge%total_volume = report%total_volume
        gauge%dm = transfer(dm%v, gauge%dm)
        gauge%active = .true.

        cerr = 0
        call petscds_set_snes_local_residual(gauge%dm, c_funloc(gauge_residual_local), c_null_ptr, cerr)
        if (cerr == 0) call petscds_set_snes_local_jacobian(gauge%dm, c_funloc(gauge_jacobian_local), &
                                                            c_null_ptr, cerr)
        if (cerr /= 0) then
            gauge%active = .false.
            ierr = int(cerr)
        end if

    end subroutine install_pressure_gauge

    !> Forget the installed gauge, so the next solver starts from nothing.
    subroutine clear_pressure_gauge()
        implicit none

        gauge = type_gauge_state()
    end subroutine clear_pressure_gauge

    !> The largest entry of one conservation row of the ungauged Jacobian.
    !>
    !> Assembled into a matrix of its own and read as \(\|J^{T}e_{K_0}\|_\infty\),
    !> which needs no rank-local matrix query and is the same number on every
    !> rank. The Jacobian of this stage does not depend on the state, so the
    !> probe is taken at zero.
    subroutine measure_row_scale(dm, row, row_scale, ierr)
        implicit none
        DM, intent(in) :: dm
        !> Global index of the row to measure.
        PetscInt, intent(in) :: row
        !> The largest absolute entry of that row.
        real(real64), intent(inout) :: row_scale
        PetscErrorCode, intent(inout) :: ierr

        Mat :: probe
        Vec :: local_state, unit_row, transposed
        PetscReal :: measured
        PetscInt :: rows(1)
        PetscScalar :: values(1)
        integer(c_int) :: cerr
        type(c_ptr) :: dm_handle, state_handle, probe_handle

        row_scale = 0.0_real64

        call DMCreateMatrix(dm, probe, ierr)
        if (ierr /= 0) return
        call DMGetLocalVector(dm, local_state, ierr)
        if (ierr == 0) call VecZeroEntries(local_state, ierr)
        if (ierr /= 0) then
            call MatDestroy(probe, ierr)
            return
        end if

        dm_handle = transfer(dm%v, dm_handle)
        state_handle = transfer(local_state%v, state_handle)
        probe_handle = transfer(probe%v, probe_handle)

        cerr = 0
        call petscds_snes_local_jacobian_fem(dm_handle, state_handle, probe_handle, probe_handle, &
                                             c_null_ptr, cerr)
        call DMRestoreLocalVector(dm, local_state, ierr)
        if (cerr /= 0) then
            call MatDestroy(probe, ierr)
            ierr = int(cerr)
            return
        end if

        call DMCreateGlobalVector(dm, unit_row, ierr)
        if (ierr == 0) call VecDuplicate(unit_row, transposed, ierr)
        if (ierr == 0) call VecZeroEntries(unit_row, ierr)
        if (ierr /= 0) then
            call MatDestroy(probe, ierr)
            return
        end if

        rows(1) = row
        values(1) = 1.0_PETSC_REAL_KIND
        call VecSetValues(unit_row, 1_PETSC_INT_KIND, rows, values, INSERT_VALUES, ierr)
        if (ierr == 0) call VecAssemblyBegin(unit_row, ierr)
        if (ierr == 0) call VecAssemblyEnd(unit_row, ierr)
        if (ierr == 0) call MatMultTranspose(probe, unit_row, transposed, ierr)
        if (ierr == 0) call VecNorm(transposed, NORM_INFINITY, measured, ierr)
        if (ierr == 0) row_scale = real(measured, real64)

        call VecDestroy(unit_row, ierr)
        call VecDestroy(transposed, ierr)
        call MatDestroy(probe, ierr)
    end subroutine measure_row_scale

    !> The local residual with the gauge equation added to it.
    function gauge_residual_local(dm_handle, state_handle, residual_handle, ctx) &
        bind(C, name="ftcms_mixed_gauge_residual") result(ierr)
        implicit none
        type(c_ptr), value :: dm_handle
        type(c_ptr), value :: state_handle
        type(c_ptr), value :: residual_handle
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        DM :: dm
        Vec :: local_state, local_residual
        PetscScalar, pointer :: raw(:)
        integer(c_int) :: cerr
        real(real64) :: weighted

        cerr = 0
        call petscds_snes_local_residual_fem(dm_handle, state_handle, residual_handle, ctx, cerr)
        ierr = int(cerr, kind(ierr))
        if (cerr /= 0) return
        if (.not. gauge%active) return

        dm%v = transfer(dm_handle, dm%v)
        local_state%v = transfer(state_handle, local_state%v)
        local_residual%v = transfer(residual_handle, local_residual%v)
        call check_installed_dm(dm_handle)

        call integrate_cell_pressure(dm, local_state, weighted, ierr)
        if (ierr /= 0) return
        if (gauge%local_offset < 0) return

        call VecGetArray(local_residual, raw, ierr)
        if (ierr /= 0) return
        raw(gauge%local_offset + 1) = raw(gauge%local_offset + 1) &
                                      + real(gauge%scale*weighted/gauge%total_volume, PETSC_REAL_KIND)
        call VecRestoreArray(local_residual, raw, ierr)
    end function gauge_residual_local

    !> The local Jacobian with the gauge row added to it.
    function gauge_jacobian_local(dm_handle, state_handle, jacobian_handle, preconditioner_handle, ctx) &
        bind(C, name="ftcms_mixed_gauge_jacobian") result(ierr)
        implicit none
        type(c_ptr), value :: dm_handle
        type(c_ptr), value :: state_handle
        type(c_ptr), value :: jacobian_handle
        type(c_ptr), value :: preconditioner_handle
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        DM :: dm
        Mat :: jacobian, preconditioner
        integer(c_int) :: cerr

        cerr = 0
        call petscds_snes_local_jacobian_fem(dm_handle, state_handle, jacobian_handle, &
                                             preconditioner_handle, ctx, cerr)
        ierr = int(cerr, kind(ierr))
        if (cerr /= 0) return
        if (.not. gauge%active) return

        dm%v = transfer(dm_handle, dm%v)
        call check_installed_dm(dm_handle)

        preconditioner%v = transfer(preconditioner_handle, preconditioner%v)
        call add_gauge_row(dm, preconditioner, ierr)
        if (ierr /= 0) return

        ! PETSc hands the same matrix twice unless a separate preconditioning
        ! matrix was asked for; the gauge belongs in whichever one is solved.
        if (.not. c_associated(jacobian_handle, preconditioner_handle)) then
            jacobian%v = transfer(jacobian_handle, jacobian%v)
            call add_gauge_row(dm, jacobian, ierr)
        end if
    end function gauge_jacobian_local

    !> Add \(\gamma V_K/V_\Omega\) to the gauge row, one entry per owned cell.
    subroutine add_gauge_row(dm, matrix, ierr)
        implicit none
        DM, intent(in) :: dm
        Mat, intent(inout) :: matrix
        PetscErrorCode, intent(inout) :: ierr

        PetscInt :: cell, c_start, c_end, global_start, global_end
        PetscInt :: rows(1)
        PetscInt, allocatable :: columns(:)
        PetscScalar, allocatable :: values(:)
        real(real64) :: volume
        integer(int32) :: owned

        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return

        allocate (columns(max(int(c_end - c_start, int32), 1)))
        allocate (values(max(int(c_end - c_start, int32), 1)))
        columns = 0_PETSC_INT_KIND
        values = 0.0_PETSC_REAL_KIND

        owned = 0
        do cell = c_start, c_end - 1
            call DMPlexGetPointGlobal(dm, cell, global_start, global_end, ierr)
            if (ierr /= 0) exit
            if (global_start < 0) cycle
            call cell_volume(dm, cell, volume, ierr)
            if (ierr /= 0) exit
            owned = owned + 1
            columns(owned) = global_start
            values(owned) = real(gauge%scale*volume/gauge%total_volume, PETSC_REAL_KIND)
        end do
        if (ierr /= 0) then
            deallocate (columns, values)
            return
        end if

        ! The gauge couples one row to every cell of the mesh, which the
        ! element-loop preallocation cannot know about.
        rows(1) = gauge%row
        call MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
        if (ierr == 0) call MatSetValues(matrix, 1_PETSC_INT_KIND, rows, int(owned, PETSC_INT_KIND), &
                                         columns, values, ADD_VALUES, ierr)
        if (ierr == 0) call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
        if (ierr == 0) call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)

        deallocate (columns, values)
    end subroutine add_gauge_row

    !> \(\sum_K V_K\hat p_K\) over the owned cells of every rank.
    subroutine integrate_cell_pressure(dm, local_state, weighted, ierr)
        implicit none
        DM, intent(in) :: dm
        Vec, intent(in) :: local_state
        !> The volume-weighted sum of the dimensionless cell pressures, m2.
        real(real64), intent(inout) :: weighted
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: cell, c_start, c_end, offset, num_dofs
        real(real64) :: volume, send(1), reduced(1)
        integer(int32) :: mpi_err

        weighted = 0.0_real64
        send(1) = 0.0_real64

        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr == 0) call DMGetLocalSection(dm, section, ierr)
        if (ierr /= 0) return

        call VecGetArrayRead(local_state, raw, ierr)
        if (ierr /= 0) return
        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(dm, cell)) cycle
            call PetscSectionGetFieldDof(section, cell, MIXED_FIELD_PRESSURE, num_dofs, ierr)
            if (ierr /= 0 .or. num_dofs <= 0) cycle
            call PetscSectionGetFieldOffset(section, cell, MIXED_FIELD_PRESSURE, offset, ierr)
            if (ierr /= 0) cycle
            call cell_volume(dm, cell, volume, ierr)
            if (ierr /= 0) exit
            send(1) = send(1) + volume*real(raw(offset + 1), real64)
        end do
        call VecRestoreArrayRead(local_state, raw, ierr)
        if (ierr /= 0) return

        ! Every rank evaluates the residual, so this reduction is reached by all
        ! of them and the gauge term is the same number everywhere.
        call MPI_Allreduce(send, reduced, 1, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, mpi_err)
        if (mpi_err /= 0) then
            ierr = 1
            return
        end if
        weighted = reduced(1)
    end subroutine integrate_cell_pressure

    !> Refuse to evaluate a gauge that belongs to a different mesh.
    subroutine check_installed_dm(dm_handle)
        implicit none
        type(c_ptr), intent(in) :: dm_handle

        if (.not. c_associated(dm_handle, gauge%dm)) then
            error stop "The mixed hydraulic pressure gauge was installed for a different mesh"
        end if
    end subroutine check_installed_dm

    !> PETSc's global number of a point, whether or not this rank owns it.
    pure function global_index_of(number) result(index)
        implicit none
        PetscInt, intent(in) :: number
        integer(int32) :: index

        if (number < 0) then
            index = int(-(number + 1), int32)
        else
            index = int(number, int32)
        end if
    end function global_index_of

end module app_mixed_hydraulic_gauge
