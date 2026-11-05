program main
    use, intrinsic :: iso_fortran_env, only: int32
    use :: petscsys, only:PetscInitialize, PetscFinalize
    use :: module_ftcms

    implicit none

    integer(int32) :: ierr
    type(type_ftcms) :: ftcms

    ! PETSc starts and stops MPI for the whole run.
    call PetscInitialize(ierr)
    if (ierr /= 0) error stop "PETSc could not be initialised"

    call ftcms%initialize()
    call ftcms%run()
    call ftcms%destroy()

    call PetscFinalize(ierr)

end program main
