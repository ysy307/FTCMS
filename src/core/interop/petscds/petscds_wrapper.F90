!> Fortran access to the PetscDS calls that take a pointwise function.
!>
!> PETSc ships no Fortran stubs for PetscDSSetResidual, PetscDSSetJacobian,
!> DMAddBoundary or the weak form of a boundary, because each of them takes a
!> C function pointer. The C side of this module passes those pointers through
!> unchanged, so the weak form stays in Fortran.
!>
!> The DMSNES local hooks are here for the same reason. They take a function
!> pointer and have no Fortran interface, so a solver that needs to add a term
!> PETSc's element loop cannot see has no other way to reach them. Those four
!> declarations bind straight to the PETSc symbols: none of them needs anything
!> done on the C side, and a pass-through shim would only add a place for the
!> prototypes to disagree.
!>
!> PETSc objects reach C as their opaque handle: a Fortran PetscDS or DM keeps
!> the address in its %v component, so the caller passes transfer(ds%v, handle)
!> the same way domain_fe_petsc does for a PetscFE.
module core_interop_petscds
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_funptr, c_char, c_null_char, c_null_ptr
    implicit none
    private

    public :: petscds_set_residual
    public :: petscds_set_jacobian
    public :: petscds_add_boundary_essential
    public :: petscds_add_boundary_natural
    public :: petscds_project_function
    public :: petscds_set_snes_local_residual
    public :: petscds_set_snes_local_jacobian
    public :: petscds_snes_local_residual_fem
    public :: petscds_snes_local_jacobian_fem

    interface
        function ftcms_ds_set_residual(ds, f, f0, f1) &
            bind(c, name="ftcms_ds_set_residual") result(ierr)
            import :: c_int, c_ptr, c_funptr
            implicit none
            type(c_ptr), value :: ds
            integer(c_int), value :: f
            type(c_funptr), value :: f0
            type(c_funptr), value :: f1
            integer(c_int) :: ierr
        end function ftcms_ds_set_residual

        function ftcms_ds_set_jacobian(ds, f, g, g0, g1, g2, g3) &
            bind(c, name="ftcms_ds_set_jacobian") result(ierr)
            import :: c_int, c_ptr, c_funptr
            implicit none
            type(c_ptr), value :: ds
            integer(c_int), value :: f
            integer(c_int), value :: g
            type(c_funptr), value :: g0
            type(c_funptr), value :: g1
            type(c_funptr), value :: g2
            type(c_funptr), value :: g3
            integer(c_int) :: ierr
        end function ftcms_ds_set_jacobian

        function ftcms_dm_add_boundary_essential(dm, name, label_name, value, field, &
                                                 num_components, components, func, ctx, bd) &
            bind(c, name="ftcms_dm_add_boundary_essential") result(ierr)
            import :: c_int, c_ptr, c_funptr, c_char
            implicit none
            type(c_ptr), value :: dm
            character(kind=c_char), intent(in) :: name(*)
            character(kind=c_char), intent(in) :: label_name(*)
            integer(c_int), value :: value
            integer(c_int), value :: field
            integer(c_int), value :: num_components
            type(c_ptr), value :: components
            type(c_funptr), value :: func
            type(c_ptr), value :: ctx
            integer(c_int), intent(inout) :: bd
            integer(c_int) :: ierr
        end function ftcms_dm_add_boundary_essential

        function ftcms_dm_add_boundary_natural(dm, name, label_name, value, field, f0, g0, bd) &
            bind(c, name="ftcms_dm_add_boundary_natural") result(ierr)
            import :: c_int, c_ptr, c_char, c_funptr
            implicit none
            type(c_ptr), value :: dm
            character(kind=c_char), intent(in) :: name(*)
            character(kind=c_char), intent(in) :: label_name(*)
            integer(c_int), value :: value
            integer(c_int), value :: field
            type(c_funptr), value :: f0
            type(c_funptr), value :: g0
            integer(c_int), intent(inout) :: bd
            integer(c_int) :: ierr
        end function ftcms_dm_add_boundary_natural

        function ftcms_dm_project_function(dm, time, funcs, ctxs, mode, x) &
            bind(c, name="ftcms_dm_project_function") result(ierr)
            import :: c_int, c_double, c_ptr, c_funptr
            implicit none
            type(c_ptr), value :: dm
            real(c_double), value :: time
            type(c_funptr), intent(in) :: funcs(*)
            type(c_ptr), value :: ctxs
            integer(c_int), value :: mode
            type(c_ptr), value :: x
            integer(c_int) :: ierr
        end function ftcms_dm_project_function

        function dm_snes_set_function_local(dm, func, ctx) &
            bind(c, name="DMSNESSetFunctionLocal") result(ierr)
            import :: c_int, c_ptr, c_funptr
            implicit none
            type(c_ptr), value :: dm
            type(c_funptr), value :: func
            type(c_ptr), value :: ctx
            integer(c_int) :: ierr
        end function dm_snes_set_function_local

        function dm_snes_set_jacobian_local(dm, func, ctx) &
            bind(c, name="DMSNESSetJacobianLocal") result(ierr)
            import :: c_int, c_ptr, c_funptr
            implicit none
            type(c_ptr), value :: dm
            type(c_funptr), value :: func
            type(c_ptr), value :: ctx
            integer(c_int) :: ierr
        end function dm_snes_set_jacobian_local

        function dm_plex_snes_residual_fem(dm, local_x, local_f, ctx) &
            bind(c, name="DMPlexSNESComputeResidualFEM") result(ierr)
            import :: c_int, c_ptr
            implicit none
            type(c_ptr), value :: dm
            type(c_ptr), value :: local_x
            type(c_ptr), value :: local_f
            type(c_ptr), value :: ctx
            integer(c_int) :: ierr
        end function dm_plex_snes_residual_fem

        function dm_plex_snes_jacobian_fem(dm, local_x, jacobian, preconditioner, ctx) &
            bind(c, name="DMPlexSNESComputeJacobianFEM") result(ierr)
            import :: c_int, c_ptr
            implicit none
            type(c_ptr), value :: dm
            type(c_ptr), value :: local_x
            type(c_ptr), value :: jacobian
            type(c_ptr), value :: preconditioner
            type(c_ptr), value :: ctx
            integer(c_int) :: ierr
        end function dm_plex_snes_jacobian_fem

    end interface

contains

    !> Residual of field `f`: f0 multiplies the test function, f1 its gradient.
    !> Either may be c_null_funptr.
    subroutine petscds_set_residual(ds, field, f0, f1, ierr)
        implicit none
        type(c_ptr), intent(in) :: ds
        integer(c_int), intent(in) :: field
        type(c_funptr), intent(in) :: f0
        type(c_funptr), intent(in) :: f1
        integer(c_int), intent(inout) :: ierr

        ierr = ftcms_ds_set_residual(ds, field, f0, f1)
    end subroutine petscds_set_residual

    !> Jacobian block (field, field_col) of the residual above:
    !> g0 = d(f0)/du, g1 = d(f0)/d(grad u), g2 = d(f1)/du, g3 = d(f1)/d(grad u).
    subroutine petscds_set_jacobian(ds, field, field_col, g0, g1, g2, g3, ierr)
        implicit none
        type(c_ptr), intent(in) :: ds
        integer(c_int), intent(in) :: field
        integer(c_int), intent(in) :: field_col
        type(c_funptr), intent(in) :: g0
        type(c_funptr), intent(in) :: g1
        type(c_funptr), intent(in) :: g2
        type(c_funptr), intent(in) :: g3
        integer(c_int), intent(inout) :: ierr

        ierr = ftcms_ds_set_jacobian(ds, field, field_col, g0, g1, g2, g3)
    end subroutine petscds_set_jacobian

    !> Constrain `field` on the points of `label_name` that carry `value`.
    !> `func` evaluates the prescribed value at a point; every component is
    !> constrained. `ctx` is handed back to `func` unchanged, which is how one
    !> callback can serve several boundaries.
    subroutine petscds_add_boundary_essential(dm, name, label_name, value, field, func, ctx, bd, ierr)
        implicit none
        type(c_ptr), intent(in) :: dm
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: label_name
        integer(c_int), intent(in) :: value
        integer(c_int), intent(in) :: field
        type(c_funptr), intent(in) :: func
        type(c_ptr), intent(in) :: ctx
        integer(c_int), intent(inout) :: bd
        integer(c_int), intent(inout) :: ierr

        ierr = ftcms_dm_add_boundary_essential(dm, trim(name)//c_null_char, &
                                               trim(label_name)//c_null_char, &
                                               value, field, 0_c_int, c_null_ptr, &
                                               func, ctx, bd)
    end subroutine petscds_add_boundary_essential

    !> Mark the points of `label_name` carrying `value` as a natural boundary
    !> for `field`, and give PETSc the kernels that evaluate it there.
    !>
    !> The kernels belong to the boundary rather than to the discrete system:
    !> PETSc looks them up in that boundary's own weak form, and a boundary
    !> registered without them contributes nothing at all.
    subroutine petscds_add_boundary_natural(dm, name, label_name, value, field, f0, g0, bd, ierr)
        implicit none
        type(c_ptr), intent(in) :: dm
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: label_name
        integer(c_int), intent(in) :: value
        integer(c_int), intent(in) :: field
        type(c_funptr), intent(in) :: f0
        type(c_funptr), intent(in) :: g0
        integer(c_int), intent(inout) :: bd
        integer(c_int), intent(inout) :: ierr

        ierr = ftcms_dm_add_boundary_natural(dm, trim(name)//c_null_char, &
                                             trim(label_name)//c_null_char, &
                                             value, field, f0, g0, bd)
    end subroutine petscds_add_boundary_natural

    !> Interpolate one analytic function per field into `x`.
    !>
    !> Each field receives its own dual functionals, so a prescribed vector
    !> function reaches an H(div) field as its facet normal moments rather than
    !> as point values. `funcs` must hold one entry per field of `dm`; `mode` is
    !> a PETSc InsertMode.
    subroutine petscds_project_function(dm, time, funcs, mode, x, ierr)
        implicit none
        type(c_ptr), intent(in) :: dm
        real(c_double), intent(in) :: time
        type(c_funptr), intent(in) :: funcs(:)
        integer(c_int), intent(in) :: mode
        type(c_ptr), intent(in) :: x
        integer(c_int), intent(inout) :: ierr

        ierr = ftcms_dm_project_function(dm, time, funcs, c_null_ptr, mode, x)
    end subroutine petscds_project_function

    !> Replace the local residual evaluation the DM's SNES uses.
    !>
    !> `func` receives the DM, the local state and the local residual, and is
    !> called inside PETSc's own scatter and boundary-value glue. It must be a
    !> `bind(C)` procedure; registering it after `DMPlexSetSNESLocalFEM` is what
    !> lets a caller add a term to the residual PETSc's element loop assembled.
    subroutine petscds_set_snes_local_residual(dm, func, ctx, ierr)
        implicit none
        type(c_ptr), intent(in) :: dm
        type(c_funptr), intent(in) :: func
        type(c_ptr), intent(in) :: ctx
        integer(c_int), intent(inout) :: ierr

        ierr = dm_snes_set_function_local(dm, func, ctx)
    end subroutine petscds_set_snes_local_residual

    !> Replace the local Jacobian evaluation the DM's SNES uses.
    subroutine petscds_set_snes_local_jacobian(dm, func, ctx, ierr)
        implicit none
        type(c_ptr), intent(in) :: dm
        type(c_funptr), intent(in) :: func
        type(c_ptr), intent(in) :: ctx
        integer(c_int), intent(inout) :: ierr

        ierr = dm_snes_set_jacobian_local(dm, func, ctx)
    end subroutine petscds_set_snes_local_jacobian

    !> PETSc's own local FEM residual, so a replacement hook can still use it.
    subroutine petscds_snes_local_residual_fem(dm, local_x, local_f, ctx, ierr)
        implicit none
        type(c_ptr), intent(in) :: dm
        type(c_ptr), intent(in) :: local_x
        type(c_ptr), intent(in) :: local_f
        type(c_ptr), intent(in) :: ctx
        integer(c_int), intent(inout) :: ierr

        ierr = dm_plex_snes_residual_fem(dm, local_x, local_f, ctx)
    end subroutine petscds_snes_local_residual_fem

    !> PETSc's own local FEM Jacobian. It zeroes and assembles the matrix it
    !> is given, so a caller adding entries afterwards assembles again.
    subroutine petscds_snes_local_jacobian_fem(dm, local_x, jacobian, preconditioner, ctx, ierr)
        implicit none
        type(c_ptr), intent(in) :: dm
        type(c_ptr), intent(in) :: local_x
        type(c_ptr), intent(in) :: jacobian
        type(c_ptr), intent(in) :: preconditioner
        type(c_ptr), intent(in) :: ctx
        integer(c_int), intent(inout) :: ierr

        ierr = dm_plex_snes_jacobian_fem(dm, local_x, jacobian, preconditioner, ctx)
    end subroutine petscds_snes_local_jacobian_fem

end module core_interop_petscds
