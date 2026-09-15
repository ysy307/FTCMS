#include <petsc/finclude/petscdmplex.h>

!> For one cell: what does the closure return, what stratum is each point in,
!> and what coordinate does each carry?
program vtu_probe
    use, intrinsic :: iso_fortran_env, only: real64, int32, output_unit
    use :: petscdmplex
    implicit none
    PetscErrorCode :: ierr
    DM :: dm
    PetscSection :: sec
    Vec :: coords
    PetscScalar, pointer :: xy(:)
    PetscInt, pointer :: closure(:)
    PetscInt :: nc, i, dof, off, cdim, vs, ve, es, ee, cs, ce
    character(len=256) :: path
    character(len=8) :: stratum
    integer(int32) :: node

    call PetscInitialize(ierr)
    call get_command_argument(1, path)
    call DMPlexCreateFromFile(PETSC_COMM_WORLD, trim(path), "mesh", PETSC_TRUE, dm, ierr)
    call DMGetCoordinateDim(dm, cdim, ierr)
    call DMGetCoordinateSection(dm, sec, ierr)
    call DMGetCoordinatesLocal(dm, coords, ierr)
    call DMPlexGetDepthStratum(dm, 0_PETSC_INT_KIND, vs, ve, ierr)
    call DMPlexGetDepthStratum(dm, 1_PETSC_INT_KIND, es, ee, ierr)
    call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, cs, ce, ierr)
    write (output_unit,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
        "cells [", cs, ",", ce, ")  vertices [", vs, ",", ve, ")  edges [", es, ",", ee, ")"

    nullify (closure)
    call DMPlexGetTransitiveClosure(dm, cs, PETSC_TRUE, nc, closure, ierr)
    call VecGetArrayRead(coords, xy, ierr)
    write (output_unit,'(A)') "closure of the first cell (point, orientation, kind, node, xy):"
    do i = 1, 2*nc, 2
        call PetscSectionGetDof(sec, closure(i), dof, ierr)
        if (dof <= 0) cycle
        call PetscSectionGetOffset(sec, closure(i), off, ierr)
        node = int(off, int32) / int(cdim, int32) + 1
        if (closure(i) >= vs .and. closure(i) < ve) then
            stratum = "vertex"
        else if (closure(i) >= es .and. closure(i) < ee) then
            stratum = "edge"
        else
            stratum = "cell"
        end if
        write (output_unit,'(A,I7,A,I3,A,A8,A,I7,A,2F11.5)') "  p=", closure(i), " o=", closure(i+1), &
            "  ", stratum, "  node=", node, "  xy=", &
            real(xy(off+1), real64), real(xy(off+2), real64)
    end do
    call VecRestoreArrayRead(coords, xy, ierr)
    call DMPlexRestoreTransitiveClosure(dm, cs, PETSC_TRUE, nc, closure, ierr)

    ! The same cell, gathered the way PETSc gathers a cell's field values.
    block
        PetscScalar, pointer :: cell_xy(:)
        PetscInt :: csize
        integer(int32) :: k
        nullify (cell_xy)
        csize = 0
        call DMPlexVecGetClosure(dm, sec, coords, cs, csize, cell_xy, ierr)
        write (output_unit,'(A,I0,A,I0)') "DMPlexVecGetClosure ierr=", ierr, "  values=", csize
        if (ierr == 0) then
            do k = 1, int(csize, int32) / int(cdim, int32)
                write (output_unit,'(A,I2,A,2F11.5)') "   local node", k, " xy=", &
                    real(cell_xy((k-1)*cdim + 1), real64), real(cell_xy((k-1)*cdim + 2), real64)
            end do
            call DMPlexVecRestoreClosure(dm, sec, coords, cs, csize, cell_xy, ierr)
        end if
    end block
    call PetscFinalize(ierr)
end program vtu_probe
