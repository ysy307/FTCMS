!>
!> A factory for creating concrete finite element objects based on FE_TYPE IDs.
!>
module domain_fe_factory
    use, intrinsic :: iso_fortran_env
    use :: stdlib_logger, only:global_logger
    use :: stdlib_strings, only:to_string
    use :: module_core, only:FE_TYPE
    use :: domain_base_fe, only:abst_fe
    use :: domain_fe_petsc, only:construct_fe_petsc

    implicit none
    private

    public :: create_fe

contains

    !>
    !> Creates a concrete finite element object based on an FE_TYPE ID.
    !>
    function create_fe(id, integration_order) result(fe)
        implicit none
        !> The FE_TYPE ID for the element to create.
        integer(int32), intent(in) :: id
        !> The integration order for the element.
        integer(int32), intent(in) :: integration_order
        !> The newly allocated finite element object.
        class(abst_fe), allocatable :: fe

        character(len=*), parameter :: func_name = "create_fe"

        ! The basis comes from PETSc, in the same order the DMPlex closure
        ! numbers the nodes, so no element type needs a convention of its own.
        select case (id)
        case (3)      ! line, 2 nodes
            fe = construct_fe_petsc(id, 1, 1, .false., integration_order)
        case (21)     ! quadratic edge, 3 nodes
            fe = construct_fe_petsc(id, 1, 2, .false., integration_order)
        case (68)     ! cubic edge, 4 nodes
            fe = construct_fe_petsc(id, 1, 3, .false., integration_order)
        case (5)      ! triangle, 3 nodes
            fe = construct_fe_petsc(id, 2, 1, .true., integration_order)
        case (22)     ! quadratic triangle, 6 nodes
            fe = construct_fe_petsc(id, 2, 2, .true., integration_order)
        case (69)     ! cubic triangle, 10 nodes
            fe = construct_fe_petsc(id, 2, 3, .true., integration_order)
        case (9)      ! quadrilateral, 4 nodes
            fe = construct_fe_petsc(id, 2, 1, .false., integration_order)
        case (28)     ! biquadratic quadrilateral, 9 nodes
            fe = construct_fe_petsc(id, 2, 2, .false., integration_order)
        case (70)     ! bicubic quadrilateral, 16 nodes
            fe = construct_fe_petsc(id, 2, 3, .false., integration_order)
        case (10)     ! tetrahedron, 4 nodes
            fe = construct_fe_petsc(id, 3, 1, .true., integration_order)
        case (24)     ! quadratic tetrahedron, 10 nodes
            fe = construct_fe_petsc(id, 3, 2, .true., integration_order)
        case (71)     ! cubic tetrahedron, 20 nodes
            fe = construct_fe_petsc(id, 3, 3, .true., integration_order)
        case (12)     ! hexahedron, 8 nodes
            fe = construct_fe_petsc(id, 3, 1, .false., integration_order)
        case (29)     ! triquadratic hexahedron, 27 nodes
            fe = construct_fe_petsc(id, 3, 2, .false., integration_order)
        case (72)     ! cubic hexahedron, 64 nodes
            fe = construct_fe_petsc(id, 3, 3, .false., integration_order)
        case default
            call global_logger%log_error( &
                func_name//": No element for ID = "//to_string(id))
        end select

    end function create_fe

end module domain_fe_factory
