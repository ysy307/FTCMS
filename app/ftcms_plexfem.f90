program main
    use, intrinsic :: iso_fortran_env, only: int32
#ifdef _MPI
    use :: core_parallel_mpi
#endif
    use :: app_plexfem, only:type_plexfem

    implicit none

    integer(int32) :: ierr

#ifdef _MPI
    call MPI_Init(ierr)
#endif

    call run_plexfem_solver()

#ifdef _MPI
    call MPI_Finalize(ierr)
#endif

contains

    subroutine run_plexfem_solver()
        implicit none
        type(type_plexfem), target :: plexfem

        call plexfem%initialize()
        call plexfem%run()
        call plexfem%destroy()

    end subroutine run_plexfem_solver

end program main
