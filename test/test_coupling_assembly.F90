#include <petsc/finclude/petscmat.h>

program test_coupling_assembly
    use, intrinsic :: iso_fortran_env, only:error_unit, int32, output_unit, real64
    use :: petscmat
    use :: module_core, only:PHYSICS_TYPES
    use :: module_ftcms, only:type_ftcms
    implicit none

    type(type_ftcms) :: ftcms
    Mat :: matrix
    integer(int32) :: num_dofs_per_node
    real(real64), pointer, contiguous :: field(:)
    real(real64) :: norm_th, norm_ht
    real(real64) :: thermal_sum_before, thermal_sum_after
    real(real64) :: block_error(5)
    integer(int32) :: ierr, i, num_nodes, n_dim
    real(real64), allocatable :: node_coord(:)
    real(real64) :: z_min, z_max
    !> The finite-difference coupling block is a difference quotient, so it
    !> carries truncation error; three digits is far tighter than the 1.0
    !> relative error the storage-only block showed and far looser than the
    !> 1.8e-6 actually measured.
    real(real64), parameter :: FD_BLOCK_TOLERANCE = 1.0d-3
    !> Sample the freezing-front elements; the audit is O(n_node) element
    !> residual evaluations per element and the test does not need all of them.
    integer(int32), parameter :: AUDIT_STRIDE = 37
#ifdef _MPI
    logical :: mpi_is_initialized, mpi_is_finalized
#endif

#ifdef _MPI
    call MPI_Initialized(mpi_is_initialized, ierr)
    if (.not. mpi_is_initialized) call MPI_Init(ierr)
#endif

    nullify (field)
    call ftcms%initialize()

    ! ------------------------------------------------------------------
    ! Impose a temperature ramp through the freezing point.
    !
    ! The initial condition is uniform, so grad T vanishes and every
    ! temperature-driven coupling term with it - a state in which a coupling
    ! block can be badly wrong and still look right. The ramp restores a
    ! finite gradient and puts part of the mesh below 0 degC, which is where
    ! the constitutive derivatives actually bite.
    ! ------------------------------------------------------------------
    call ftcms%domain%get_num_nodes(num_nodes)
    call ftcms%temperature%get_current(field)
    if (.not. associated(field)) then
        write (error_unit, '(A)') "FAIL: temperature field is not associated"
        error stop 1
    end if
    ! Ramp along the vertical coordinate, not the node index: node numbering
    ! carries no spatial meaning, so indexing by it would impose an oscillating
    ! field with meaningless gradients rather than a freezing front.
    call ftcms%domain%nodes%get_dimension(n_dim)
    allocate (node_coord(n_dim))
    z_min = huge(1.0d0)
    z_max = -huge(1.0d0)
    do i = 1, num_nodes
        call ftcms%domain%nodes%get_coordinate(i, node_coord)
        z_min = min(z_min, node_coord(size(node_coord)))
        z_max = max(z_max, node_coord(size(node_coord)))
    end do
    do i = 1, num_nodes
        call ftcms%domain%nodes%get_coordinate(i, node_coord)
        field(i) = 6.0d0 - 1.2d1 * (node_coord(size(node_coord)) - z_min) / max(z_max - z_min, tiny(1.0d0))
    end do
    nullify (field)

    call ftcms%assemble()

    call ftcms%assemble(residual_only=.true.)
    
    call ftcms%K%assemble()
    matrix = ftcms%K%get_matrix()
    if (matrix%v == PETSC_NULL_MAT%v) then
        write (error_unit, '(A)') "FAIL: monolithic Jacobian matrix does not exist"
        error stop 1
    end if
    call ftcms%K%get_num_dofs_per_node(num_dofs_per_node)

    norm_th = block_extreme(matrix, num_dofs_per_node, PHYSICS_TYPES%THERMAL%ID, &
                            PHYSICS_TYPES%HYDRAULIC%ID, .true.)
    norm_ht = block_extreme(matrix, num_dofs_per_node, PHYSICS_TYPES%HYDRAULIC%ID, &
                            PHYSICS_TYPES%THERMAL%ID, .true.)
    if (norm_th <= tiny(1.0d0)) then
        write (error_unit, '(A,ES12.5)') "FAIL: K_TH is zero, max norm=", norm_th
        error stop 1
    end if
    if (norm_ht <= tiny(1.0d0)) then
        write (error_unit, '(A,ES12.5)') "FAIL: K_HT is zero, max norm=", norm_ht
        error stop 1
    end if
    write (output_unit, '(A,2(1X,ES12.5))') "PASS: coupled Modified-Picard K_TH/K_HT max norms:", norm_th, norm_ht

    ! ------------------------------------------------------------------
    ! A Robin condition must strengthen the diagonal it acts on.
    !
    ! K holds dR/du while F holds -R, so the boundary influx -int(psi*q)
    ! contributes -int(psi * dq/du * psi); with dq/du = -h that is the
    ! standard +h*int(psi*psi). Adding it with the opposite sign weakens the
    ! very diagonal the condition supplies, which is invisible in a max-norm
    ! check but flips the sign of this sum.
    ! ------------------------------------------------------------------
    thermal_sum_before = block_extreme(matrix, num_dofs_per_node, PHYSICS_TYPES%THERMAL%ID, &
                                       PHYSICS_TYPES%THERMAL%ID, .false.)
    call ftcms%apply_bc(prescribed=.false.)
    call ftcms%K%assemble()
    thermal_sum_after = block_extreme(matrix, num_dofs_per_node, PHYSICS_TYPES%THERMAL%ID, &
                                      PHYSICS_TYPES%THERMAL%ID, .false.)
    if (thermal_sum_after <= thermal_sum_before) then
        write (error_unit, '(A,2(1X,ES13.5))') &
            "FAIL: Robin boundary weakens the thermal block; sum before/after=", &
            thermal_sum_before, thermal_sum_after
        error stop 1
    end if
    write (output_unit, '(A,ES13.5)') "PASS: Robin boundary strengthens the thermal block by", &
        thermal_sum_after - thermal_sum_before

    ! ------------------------------------------------------------------
    ! The assembled blocks must match the element finite-difference tangent
    ! of the same residual. Slot 5 is the coupling block the assembly
    ! actually installs.
    ! ------------------------------------------------------------------
    block_error(:) = -1.0d0
    call ftcms%report_fd_jacobian('test_coupling_assembly', mean_rel_error=block_error, &
                                  element_stride=AUDIT_STRIDE)
    if (block_error(5) < 0.0d0) then
        write (error_unit, '(A)') "FAIL: finite-difference audit sampled no element"
        error stop 1
    end if
    if (block_error(5) > FD_BLOCK_TOLERANCE) then
        write (error_unit, '(A,ES12.5)') &
            "FAIL: assembled K_TH disagrees with the finite-difference tangent, mean rel err=", &
            block_error(5)
        error stop 1
    end if
    write (output_unit, '(A,ES12.5)') "PASS: K_TH matches the finite-difference tangent, mean rel err=", &
        block_error(5)

    call ftcms%destroy()

#ifdef _MPI
    call MPI_Initialized(mpi_is_initialized, ierr)
    call MPI_Finalized(mpi_is_finalized, ierr)
    if (mpi_is_initialized .and. (.not. mpi_is_finalized)) call MPI_Finalize(ierr)
#endif

contains

    !> Walk the assembled matrix and reduce one physics block of it: the largest
    !> magnitude when want_max is true, the plain sum otherwise. The block a
    !> matrix entry belongs to follows from its position inside the node block.
    function block_extreme(A, num_dofs_per_node, row_offset, col_offset, want_max) result(value)
        implicit none
        Mat, intent(in) :: A
        integer(int32), intent(in) :: num_dofs_per_node, row_offset, col_offset
        logical, intent(in) :: want_max
        real(real64) :: value

        PetscErrorCode :: ierr_local
        PetscInt :: first_row, last_row, row, num_cols
        PetscInt, pointer :: cols(:)
        PetscScalar, pointer :: vals(:)
        integer(int32) :: k

        value = 0.0d0
        call MatGetOwnershipRange(A, first_row, last_row, ierr_local)
        do row = first_row, last_row - 1
            if (int(mod(row, int(num_dofs_per_node, PETSC_INT_KIND)), int32) + 1 /= row_offset) cycle
            nullify (cols)
            nullify (vals)
            call MatGetRow(A, row, num_cols, cols, vals, ierr_local)
            do k = 1, int(num_cols, int32)
                if (int(mod(cols(k), int(num_dofs_per_node, PETSC_INT_KIND)), int32) + 1 /= col_offset) cycle
                if (want_max) then
                    value = max(value, abs(real(vals(k), real64)))
                else
                    value = value + real(vals(k), real64)
                end if
            end do
            call MatRestoreRow(A, row, num_cols, cols, vals, ierr_local)
        end do
    end function block_extreme

end program test_coupling_assembly
