#include <petsc/finclude/petscmat.h>

!> The Jacobian, held as PETSc matrices.
!>
!> PETSc owns the linear algebra: the dof numbering comes from the DM's
!> section, the sparsity and the parallel distribution come from DMCreateMatrix,
!> and a node shared by two ranks is summed by MatAssembly because both ranks
!> add into the same global index. Nothing here keeps a second copy of the
!> matrix, and the solver is handed the Mat itself.
!>
!> The public interface is unchanged: the physics still addresses entries by
!> (physics, node), which this module turns into global dof indices.
module numerical_system_jacobian_matrix
    use, intrinsic :: iso_fortran_env
    use :: petscmat
    use :: module_core
    use :: module_linalg
    use :: domain_mesh_plex, only:type_mesh_plex, type_dof_layout
    implicit none
    private

    public :: type_jacobian_matrix

    !> One linear system: monolithic coupling has a single one, staggered
    !> coupling one per active physics.
    type :: type_jacobian_system
        Mat :: matrix
        type(type_dof_layout) :: layout
        logical :: ready = .false.
        !> Rows the boundary conditions prescribe. They are applied once, on the
        !> assembled matrix, because a row cannot be rewritten while the element
        !> loop is still adding into it.
        integer(int32), allocatable :: prescribed_row(:)
        real(real64), allocatable :: prescribed_diagonal(:)
        integer(int32) :: num_prescribed = 0
    end type type_jacobian_system

    type :: type_jacobian_matrix
        private
        integer(int32) :: num_nodes = 0
        integer(int32) :: num_dofs_per_node = 0
        integer(int32) :: num_dofs_of_physics(PHYSICS_TYPES%NUM_ID) = 0
        integer(int32) :: size = 0

        type(type_constant_id) :: coupling_mode

        type(type_jacobian_system), allocatable :: system(:)
        integer(int32) :: num_system = 0

        !> physics i -> system index (0 if inactive)
        integer(int32) :: physics_to_system(PHYSICS_TYPES%NUM_ID) = 0
        !> physics i -> dof offset inside a monolithic node block (0 if inactive)
        integer(int32) :: physics_to_block(PHYSICS_TYPES%NUM_ID) = 0

        !> Element connectivity, so an element block can be turned into global
        !> row and column indices without going back to the domain.
        integer(int32), allocatable :: element_offset(:)
        integer(int32), allocatable :: element_node(:)
        integer(int32) :: num_elements = 0
    contains
        procedure, public, pass(self) :: initialize => initialize_jacobian_matrix
        procedure, public, pass(self) :: destroy => destroy_jacobian_matrix
        procedure, public, pass(self) :: assemble => assemble_jacobian_matrix

        procedure, public, pass(self) :: get_size => get_size_jacobian_matrix
        procedure, public, pass(self) :: get_num_dofs_per_node => get_num_dofs_per_node
        procedure, public, pass(self) :: get_matrix => get_underlying_matrix
        procedure, public, pass(self) :: get_dof_layout => get_dof_layout_jacobian
        procedure, public, pass(self) :: get_system_index => get_system_index_jacobian

        procedure, private, pass(self) :: resolve_block_target
        procedure, private, pass(self) :: set_value_local => set_value_jacobian_matrix
        procedure, private, pass(self) :: add_value_local => add_value_jacobian_matrix
        procedure, private, pass(self) :: add_local_matrix => add_local_jacobian_matrix
        generic, public :: set => set_value_local
        generic, public :: add => add_value_local, add_local_matrix

        procedure, private, pass(self) :: zero_all => zero_all_jacobian_matrix
        procedure, private, pass(self) :: zero_row => zero_row_jacobian_matrix
        generic, public :: zero => zero_all, zero_row
        procedure, public, pass(self) :: display => display_jacobian_matrix
    end type type_jacobian_matrix

contains

    subroutine initialize_jacobian_matrix(self, topology, coupling_mode, mesh)
        implicit none
        class(type_jacobian_matrix), intent(inout) :: self
        type(type_system_topology), intent(in) :: topology
        type(type_constant_id), intent(in) :: coupling_mode
        type(type_mesh_plex), intent(inout) :: mesh

        integer(int32) :: i, j

        call self%destroy()

        self%coupling_mode = coupling_mode
        self%physics_to_system(:) = 0
        self%physics_to_block(:) = 0
        self%num_dofs_of_physics(:) = 0

        call topology%get_total_dofs(self%size)
        call topology%get_num_nodes(self%num_nodes)

        do i = 1, PHYSICS_TYPES%NUM_ID
            call topology%get_start_dof_index(PHYSICS_TYPES%to_object(i), self%physics_to_block(i))
        end do

        select case (coupling_mode%ID)
        case (COUPLING_MODES%MONOLITHIC%ID)
            call topology%get_num_dof_per_node(self%num_dofs_per_node)
            self%num_system = 1
            allocate (self%system(1))
            call build_system(self%system(1), mesh, self%num_dofs_per_node)

        case (COUPLING_MODES%STAGGERED%ID)
            do i = 1, PHYSICS_TYPES%NUM_ID
                call topology%get_target_dof(PHYSICS_TYPES%to_object(i), self%num_dofs_of_physics(i))
            end do
            self%num_system = count(self%num_dofs_of_physics > 0)
            self%num_dofs_per_node = 0
            allocate (self%system(self%num_system))
            j = 0
            do i = 1, PHYSICS_TYPES%NUM_ID
                if (self%num_dofs_of_physics(i) == 0) cycle
                j = j + 1
                self%physics_to_system(i) = j
                call build_system(self%system(j), mesh, self%num_dofs_of_physics(i))
            end do
        end select

        call cache_element_connectivity(self, topology)
    end subroutine initialize_jacobian_matrix

    !> Ask the DM for the numbering, then for a matrix that matches it.
    subroutine build_system(system, mesh, num_dofs_per_node)
        implicit none
        type(type_jacobian_system), intent(inout) :: system
        type(type_mesh_plex), intent(inout) :: mesh
        integer(int32), intent(in) :: num_dofs_per_node

        call mesh%create_dof_layout(num_dofs_per_node, system%layout)
        call mesh%create_matrix(system%matrix)
        allocate (system%prescribed_row(0))
        allocate (system%prescribed_diagonal(0))
        system%num_prescribed = 0
        system%ready = .true.
    end subroutine build_system

    subroutine cache_element_connectivity(self, topology)
        implicit none
        class(type_jacobian_matrix), intent(inout) :: self
        type(type_system_topology), intent(in), target :: topology

        integer(int32), pointer, contiguous :: connectivity(:)
        integer(int32) :: elem_id, total, i

        call topology%get_num_fe(self%num_elements)
        if (self%num_elements == 0) return

        allocate (self%element_offset(self%num_elements + 1))
        self%element_offset(1) = 1
        do elem_id = 1, self%num_elements
            call topology%get_fe_connectivity(elem_id, connectivity)
            self%element_offset(elem_id + 1) = self%element_offset(elem_id) + size(connectivity)
        end do

        total = self%element_offset(self%num_elements + 1) - 1
        allocate (self%element_node(max(total, 1)))
        do elem_id = 1, self%num_elements
            call topology%get_fe_connectivity(elem_id, connectivity)
            do i = 1, size(connectivity)
                self%element_node(self%element_offset(elem_id) + i - 1) = connectivity(i)
            end do
        end do
    end subroutine cache_element_connectivity

    subroutine destroy_jacobian_matrix(self)
        implicit none
        class(type_jacobian_matrix), intent(inout) :: self

        PetscErrorCode :: ierr
        integer(int32) :: i

        if (allocated(self%system)) then
            do i = 1, size(self%system)
                if (self%system(i)%ready) call MatDestroy(self%system(i)%matrix, ierr)
                self%system(i)%ready = .false.
            end do
            deallocate (self%system)
        end if

        if (allocated(self%element_offset)) deallocate (self%element_offset)
        if (allocated(self%element_node)) deallocate (self%element_node)

        self%size = 0
        self%num_nodes = 0
        self%num_dofs_per_node = 0
        self%num_system = 0
        self%num_elements = 0
        self%physics_to_system(:) = 0
        self%num_dofs_of_physics(:) = 0
    end subroutine destroy_jacobian_matrix

    !> Flush the element loop into the matrix and impose the prescribed rows.
    !>
    !> The rows the boundary conditions fix cannot be rewritten while the
    !> element loop is still adding into them, and on more than one rank the
    !> contributions to a shared row are not even all present until the matrix
    !> is assembled. Both are why this happens here rather than at the call
    !> that recorded the row.
    subroutine assemble_jacobian_matrix(self)
        implicit none
        class(type_jacobian_matrix), intent(inout) :: self

        PetscErrorCode :: ierr
        PetscInt, allocatable :: rows(:)
        PetscInt :: row_index(1)
        PetscScalar :: value(1)
        integer(int32) :: i, k
        logical :: needs_diagonal

        if (.not. allocated(self%system)) return

        do i = 1, size(self%system)
            if (.not. self%system(i)%ready) cycle

            call MatAssemblyBegin(self%system(i)%matrix, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(self%system(i)%matrix, MAT_FINAL_ASSEMBLY, ierr)

            if (self%system(i)%num_prescribed <= 0) cycle

            allocate (rows(self%system(i)%num_prescribed))
            do k = 1, self%system(i)%num_prescribed
                rows(k) = int(self%system(i)%prescribed_row(k), PETSC_INT_KIND)
            end do
            call MatZeroRows(self%system(i)%matrix, int(self%system(i)%num_prescribed, PETSC_INT_KIND), &
                             rows, 1.0_PETSC_REAL_KIND, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
            deallocate (rows)

            ! MatZeroRows already leaves a unit diagonal, so only a prescribed
            ! row that asked for something else needs a second pass.
            needs_diagonal = .false.
            do k = 1, self%system(i)%num_prescribed
                if (abs(self%system(i)%prescribed_diagonal(k) - 1.0d0) > 0.0d0) needs_diagonal = .true.
            end do
            if (.not. needs_diagonal) cycle

            do k = 1, self%system(i)%num_prescribed
                row_index(1) = int(self%system(i)%prescribed_row(k), PETSC_INT_KIND)
                value(1) = self%system(i)%prescribed_diagonal(k)
                call MatSetValues(self%system(i)%matrix, 1_PETSC_INT_KIND, row_index, &
                                  1_PETSC_INT_KIND, row_index, value, INSERT_VALUES, ierr)
            end do
            call MatAssemblyBegin(self%system(i)%matrix, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(self%system(i)%matrix, MAT_FINAL_ASSEMBLY, ierr)
        end do
    end subroutine assemble_jacobian_matrix


    subroutine get_size_jacobian_matrix(self, size)
        implicit none
        class(type_jacobian_matrix), intent(in) :: self
        integer(int32), intent(inout) :: size

        size = self%size
    end subroutine get_size_jacobian_matrix

    subroutine get_num_dofs_per_node(self, num_dofs)
        implicit none
        class(type_jacobian_matrix), intent(in) :: self
        integer(int32), intent(inout) :: num_dofs

        num_dofs = self%num_dofs_per_node
    end subroutine get_num_dofs_per_node

    !> Index of the system a physics is solved in, 0 when it has none.
    integer(int32) function get_system_index_jacobian(self, physics_id_in)
        implicit none
        class(type_jacobian_matrix), intent(in) :: self
        type(type_constant_id), intent(in), optional :: physics_id_in

        integer(int32) :: physics_id

        get_system_index_jacobian = 0
        if (.not. allocated(self%system)) return

        select case (self%coupling_mode%ID)
        case (COUPLING_MODES%MONOLITHIC%ID)
            get_system_index_jacobian = 1
        case (COUPLING_MODES%STAGGERED%ID)
            if (.not. present(physics_id_in)) then
                get_system_index_jacobian = 1
                return
            end if
            physics_id = physics_id_in%ID
            if (physics_id < 1 .or. physics_id > PHYSICS_TYPES%NUM_ID) return
            get_system_index_jacobian = self%physics_to_system(physics_id)
        end select
    end function get_system_index_jacobian

    !> The PETSc matrix itself. The solver operates on this, with no copy.
    function get_underlying_matrix(self, physics_id_in) result(matrix)
        implicit none
        class(type_jacobian_matrix), intent(in) :: self
        type(type_constant_id), intent(in), optional :: physics_id_in
        Mat :: matrix

        integer(int32) :: sys_id

        matrix = PETSC_NULL_MAT
        sys_id = self%get_system_index(physics_id_in)
        if (sys_id <= 0 .or. sys_id > size(self%system)) return
        if (.not. self%system(sys_id)%ready) return
        matrix = self%system(sys_id)%matrix
    end function get_underlying_matrix

    subroutine get_dof_layout_jacobian(self, layout, physics_id_in)
        implicit none
        class(type_jacobian_matrix), intent(in) :: self
        type(type_dof_layout), intent(inout) :: layout
        type(type_constant_id), intent(in), optional :: physics_id_in

        integer(int32) :: sys_id

        sys_id = self%get_system_index(physics_id_in)
        if (sys_id <= 0 .or. sys_id > size(self%system)) return
        layout = self%system(sys_id)%layout
    end subroutine get_dof_layout_jacobian

    !> Which system a (row, column) physics pair belongs to, and the dof offset
    !> of each inside a node block. sys = 0 means the entry is not storable:
    !> an off-diagonal pair under staggered coupling, or an inactive physics.
    subroutine resolve_block_target(self, row_physics_id, col_physics_id, sys, row_blk, col_blk)
        implicit none
        class(type_jacobian_matrix), intent(in) :: self
        integer(int32), intent(in) :: row_physics_id, col_physics_id
        integer(int32), intent(inout) :: sys, row_blk, col_blk

        sys = 0
        row_blk = 0
        col_blk = 0
        if (.not. allocated(self%system)) return

        select case (self%coupling_mode%ID)
        case (COUPLING_MODES%MONOLITHIC%ID)
            if (row_physics_id < 1 .or. row_physics_id > PHYSICS_TYPES%NUM_ID) return
            if (col_physics_id < 1 .or. col_physics_id > PHYSICS_TYPES%NUM_ID) return
            row_blk = self%physics_to_block(row_physics_id)
            col_blk = self%physics_to_block(col_physics_id)
            if (row_blk <= 0 .or. col_blk <= 0) return
            sys = 1
        case (COUPLING_MODES%STAGGERED%ID)
            if (row_physics_id /= col_physics_id) return
            if (row_physics_id < 1 .or. row_physics_id > PHYSICS_TYPES%NUM_ID) return
            sys = self%physics_to_system(row_physics_id)
            if (sys <= 0 .or. sys > size(self%system)) then
                sys = 0
                return
            end if
            row_blk = 1
            col_blk = 1
        end select
    end subroutine resolve_block_target

    !> Global row of one dof of one node.
    integer(int32) function global_dof(system, node, block_offset)
        implicit none
        type(type_jacobian_system), intent(in) :: system
        integer(int32), intent(in) :: node, block_offset

        global_dof = system%layout%node_dof_base(node) + block_offset - 1
    end function global_dof

    !> Record a prescribed row. Applied by assemble, not now.
    subroutine mark_prescribed(system, row, diagonal)
        implicit none
        type(type_jacobian_system), intent(inout) :: system
        integer(int32), intent(in) :: row
        real(real64), intent(in) :: diagonal

        integer(int32), allocatable :: grown_row(:)
        real(real64), allocatable :: grown_diagonal(:)
        integer(int32) :: k

        do k = 1, system%num_prescribed
            if (system%prescribed_row(k) == row) then
                system%prescribed_diagonal(k) = diagonal
                return
            end if
        end do

        if (system%num_prescribed >= size(system%prescribed_row)) then
            allocate (grown_row(max(2 * size(system%prescribed_row), 64)))
            allocate (grown_diagonal(size(grown_row)))
            grown_row = 0
            grown_diagonal = 1.0d0
            if (system%num_prescribed > 0) then
                grown_row(1:system%num_prescribed) = system%prescribed_row(1:system%num_prescribed)
                grown_diagonal(1:system%num_prescribed) = system%prescribed_diagonal(1:system%num_prescribed)
            end if
            call move_alloc(grown_row, system%prescribed_row)
            call move_alloc(grown_diagonal, system%prescribed_diagonal)
        end if

        system%num_prescribed = system%num_prescribed + 1
        system%prescribed_row(system%num_prescribed) = row
        system%prescribed_diagonal(system%num_prescribed) = diagonal
    end subroutine mark_prescribed

    !> Used by the boundary conditions to pin a row to its prescribed value.
    subroutine set_value_jacobian_matrix(self, row_physics_id, col_physics_id, row_node, col_node, value)
        implicit none
        class(type_jacobian_matrix), intent(inout) :: self
        integer(int32), intent(in) :: row_physics_id, col_physics_id, row_node, col_node
        real(real64), intent(in) :: value

        integer(int32) :: sys, row_blk, col_blk

        call self%resolve_block_target(row_physics_id, col_physics_id, sys, row_blk, col_blk)
        if (sys == 0) return

        ! Only the diagonal of a prescribed row is ever set. Anything else would
        ! mix INSERT into an element loop that is adding, which PETSc forbids
        ! between assemblies.
        if (row_node /= col_node .or. row_blk /= col_blk) return

        call mark_prescribed(self%system(sys), global_dof(self%system(sys), row_node, row_blk), value)
    end subroutine set_value_jacobian_matrix

    subroutine add_value_jacobian_matrix(self, row_physics_id, col_physics_id, row_node, col_node, value)
        implicit none
        class(type_jacobian_matrix), intent(inout) :: self
        integer(int32), intent(in) :: row_physics_id, col_physics_id, row_node, col_node
        real(real64), intent(in) :: value

        PetscErrorCode :: ierr
        PetscInt :: row_index(1), col_index(1)
        PetscScalar :: entry(1)
        integer(int32) :: sys, row_blk, col_blk

        call self%resolve_block_target(row_physics_id, col_physics_id, sys, row_blk, col_blk)
        if (sys == 0) return

        row_index(1) = int(global_dof(self%system(sys), row_node, row_blk), PETSC_INT_KIND)
        col_index(1) = int(global_dof(self%system(sys), col_node, col_blk), PETSC_INT_KIND)
        entry(1) = value
        call MatSetValues(self%system(sys)%matrix, 1_PETSC_INT_KIND, row_index, &
                          1_PETSC_INT_KIND, col_index, entry, ADD_VALUES, ierr)
    end subroutine add_value_jacobian_matrix

    !> Add one element block. The global indices come from the dof layout, so a
    !> node shared with another rank lands on the same global row on both and
    !> MatAssembly adds the two contributions together.
    subroutine add_local_jacobian_matrix(self, row_physics_id, col_physics_id, element_id, n_local, local_data)
        implicit none
        class(type_jacobian_matrix), intent(inout) :: self
        integer(int32), intent(in) :: row_physics_id, col_physics_id, element_id, n_local
        type(type_matrix_dense), intent(in) :: local_data

        PetscErrorCode :: ierr
        PetscInt :: rows(n_local), cols(n_local)
        real(real64), pointer, dimension(:, :) :: dense_val
        integer(int32) :: sys, row_blk, col_blk, i, node

        call self%resolve_block_target(row_physics_id, col_physics_id, sys, row_blk, col_blk)
        if (sys == 0) return
        if (element_id < 1 .or. element_id > self%num_elements) return

        dense_val => local_data%get_val()
        if (.not. associated(dense_val)) return

        do i = 1, n_local
            node = self%element_node(self%element_offset(element_id) + i - 1)
            rows(i) = int(global_dof(self%system(sys), node, row_blk), PETSC_INT_KIND)
            cols(i) = int(global_dof(self%system(sys), node, col_blk), PETSC_INT_KIND)
        end do

        call MatSetValues(self%system(sys)%matrix, int(n_local, PETSC_INT_KIND), rows, &
                          int(n_local, PETSC_INT_KIND), cols, dense_val, ADD_VALUES, ierr)
        ! A refused insertion drops the whole block silently, so it must be seen.
        if (ierr /= 0) then
            write (*, '(A,I0,A,I0)') '   [ASSEMBLY] MatSetValues refused element ', element_id, &
                ', ierr = ', ierr
            error stop "Element block rejected by the matrix"
        end if
    end subroutine add_local_jacobian_matrix

    subroutine zero_all_jacobian_matrix(self)
        implicit none
        class(type_jacobian_matrix), intent(inout) :: self

        PetscErrorCode :: ierr
        integer(int32) :: i

        if (.not. allocated(self%system)) return
        do i = 1, size(self%system)
            if (.not. self%system(i)%ready) cycle
            call MatZeroEntries(self%system(i)%matrix, ierr)
            self%system(i)%num_prescribed = 0
        end do
    end subroutine zero_all_jacobian_matrix

    subroutine zero_row_jacobian_matrix(self, row_node, row_physics_id)
        implicit none
        class(type_jacobian_matrix), intent(inout) :: self
        integer(int32), intent(in) :: row_node
        integer(int32), intent(in), optional :: row_physics_id

        integer(int32) :: i, sys, row_blk, col_blk

        if (.not. allocated(self%system)) return

        ! Without a physics id, prescribe the node's row in every system.
        if (.not. present(row_physics_id)) then
            do i = 1, size(self%system)
                if (.not. self%system(i)%ready) cycle
                do sys = 1, self%system(i)%layout%num_dofs_per_node
                    call mark_prescribed(self%system(i), global_dof(self%system(i), row_node, sys), 1.0d0)
                end do
            end do
            return
        end if

        call self%resolve_block_target(row_physics_id, row_physics_id, sys, row_blk, col_blk)
        if (sys == 0) return
        call mark_prescribed(self%system(sys), global_dof(self%system(sys), row_node, row_blk), 1.0d0)
    end subroutine zero_row_jacobian_matrix

    subroutine display_jacobian_matrix(self, unit_in)
        implicit none
        class(type_jacobian_matrix), intent(in) :: self
        integer(int32), intent(in), optional :: unit_in

        PetscErrorCode :: ierr
        integer(int32) :: i

        if (.not. allocated(self%system)) return
        do i = 1, size(self%system)
            if (self%system(i)%ready) call MatView(self%system(i)%matrix, PETSC_VIEWER_STDOUT_WORLD, ierr)
        end do
    end subroutine display_jacobian_matrix

end module numerical_system_jacobian_matrix
