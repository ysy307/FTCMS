#include <petsc/finclude/petscdmplex.h>

!> The mesh, held as a PETSc DMPlex.
!>
!> There is no second copy of the mesh here. Coordinates, connectivity, cell
!> shape and the boundary patches are all answered by querying the DM, so what
!> the code works from is by construction what PETSc read from the file. The
!> only state kept alongside the DM is an index: the dense node and cell
!> numbering FTCMS uses, and the DM point each number refers to.
module domain_mesh_plex
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_ptr
    use petscdmplex
    use :: core_unique, only:unique
    use :: module_core, only:FE_TYPE
    implicit none
    private

    public :: type_mesh_plex

    public :: type_dof_layout

    !> A degree-of-freedom numbering produced by the DM.
    !>
    !> PETSc owns the numbering: the section decides how many dofs sit on each
    !> mesh point and what their global indices are, and the same section is
    !> what preallocates and partitions the matrix. node_dof_base is only the
    !> lookup the element loop needs, cached so it is not a section call per
    !> entry.
    type :: type_dof_layout
        integer(int32) :: num_dofs_per_node = 0
        integer(int32) :: num_local_dofs = 0
        integer(int32) :: num_global_dofs = 0
        !> Global index of the first dof of each local node, counting from zero.
        integer(int32), allocatable :: node_dof_base(:)
    end type type_dof_layout

    type :: type_mesh_plex
        DM :: dm
        logical :: loaded = .false.
        !> Chart of the coordinate section. Every node query needs it.
        PetscInt :: coord_start = 0
        PetscInt :: coord_end = 0
        integer(int32) :: dimension = 0

        !> Dense node numbering, taken from the coordinate section. A mesh
        !> point carries one node per coordinate value on it, which is one for a
        !> vertex but two on the edge of a cubic element and four inside it, so
        !> the count follows the section's degrees of freedom, never the number
        !> of points.
        integer(int32) :: num_nodes = 0
        integer(int32) :: coordinate_dim = 0
        !> One dof per node. PETSc walks the closure against this and hands back
        !> the node indices in the order its own basis is numbered, so the
        !> connectivity and the shape functions agree without either side
        !> knowing the other's convention.
        PetscSection :: node_section
        logical :: node_section_ready = .false.
        !> First node of each point, and how many it carries.
        integer(int32), allocatable :: point_first_node(:)
        integer(int32), allocatable :: point_num_nodes(:)
        integer(int32) :: global_num_nodes = 0
        integer(int32), allocatable :: node_global_id(:)
        logical, allocatable :: node_is_owned(:)

        !> Dense cell numbering: the interior cells first, then the faces that
        !> carry a boundary patch id.
        integer(int32) :: num_cells = 0
        integer(int32), allocatable :: cell_point(:)
        integer(int32), allocatable :: cell_dimension(:)
        integer(int32), allocatable :: cell_entity_id(:)
        integer(int32), allocatable :: cell_fe_type(:)
        integer(int32), allocatable :: cell_num_nodes(:)
        integer(int32), allocatable :: cell_color(:)
        integer(int32) :: num_colors = 0

        !> A clone of the DM carrying one dof per node, for output. It is a
        !> clone so that installing an output layout never disturbs the section
        !> the solver assembles against.
        logical :: output_dm_ready = .false.
        DM :: output_dm

        !> Reusable objects of the halo sum.
        logical :: halo_ready = .false.
        Vec :: halo_global
        Vec :: halo_local
        VecScatter :: halo_scatter
    contains
        procedure, public :: load => load_mesh_plex
        procedure, public :: destroy => destroy_mesh_plex
        procedure, public :: get_dm => get_dm_mesh_plex
        procedure, public :: get_node_coordinates => get_node_coordinates_mesh_plex
        procedure, public :: halo_sum_nodal => halo_sum_nodal_mesh_plex
        procedure, public :: create_dof_layout => create_dof_layout_mesh_plex
        procedure, public :: create_matrix => create_matrix_mesh_plex
        procedure, public :: get_output_dm => get_output_dm_mesh_plex
        procedure, public :: get_element_shape => get_element_shape_mesh_plex
        procedure, private :: coordinate_order => coordinate_order_mesh_plex
        procedure, public :: get_visualisation_shape => get_visualisation_shape_mesh_plex
        procedure, public :: get_vertex_nodes => get_vertex_nodes_mesh_plex
        procedure, public :: get_cell_connectivity => get_cell_connectivity_mesh_plex
        procedure, public :: get_active_region_info => get_active_region_info_mesh_plex
        procedure, private :: build_node_index
        procedure, private :: build_global_node_numbering
        procedure, private :: build_cell_index
        procedure, private :: build_coloring
    end type type_mesh_plex

contains

    !> Read a Gmsh file and distribute it. PETSc chooses the partition and keeps
    !> the physical groups as the "Cell Sets" and "Face Sets" labels.
    subroutine load_mesh_plex(self, file_name)
        implicit none
        class(type_mesh_plex), intent(inout) :: self
        character(*), intent(in) :: file_name

        PetscErrorCode :: ierr
        DM :: distributed, refined
        PetscSection :: section
        PetscMPIInt :: num_procs
        PetscBool :: ready, refine_found
        PetscInt :: raw_dimension, refine_levels, level

        ierr = 0
        call PetscInitialized(ready, ierr)
        if (ready .eqv. PETSC_FALSE) call PetscInitialize(ierr)

        ! Read every coordinate the file carries, not just the first ones.
        !
        ! A two-dimensional analysis of downward flow is meshed in the x-z
        ! plane, because gravity has to act along a coordinate the solve
        ! carries. Left to itself the reader keeps the first two coordinates of
        ! such a file, which are x and a y that never varies, and every cell
        ! comes out flat. Read in three, the mesh is the surface it actually
        ! is, and PETSc measures it with the metric of that surface -- so an
        ! x-y mesh and an x-z mesh both mean what they say.
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-dm_plex_gmsh_spacedim", "3", ierr)
        if (ierr /= 0) error stop "PETSc rejected the mesh reader options"

        ! interpolate = TRUE builds the faces the boundary labels hang off.
        call DMPlexCreateFromFile(PETSC_COMM_WORLD, trim(file_name), "mesh", PETSC_TRUE, self%dm, ierr)
        if (ierr /= 0) error stop "DMPlex could not read the mesh"

        ! A mesh whose coordinates PETSc can only hold one set per cell has no
        ! shared node set to solve on, and is refused rather than bent into one.
        call reject_unsupported_coordinates(self, file_name)

        ! Uniform refinement of the mesh as read, for mesh-dependence studies
        ! that must not touch the mesh file: -ftcms_refine <levels> halves every
        ! cell's length scale that many times. PETSc carries "Cell Sets" and
        ! "Face Sets" onto the children, and the embedding dimension with them,
        ! so materials, boundaries and an x-z plane survive. It runs before the
        ! distribution and before any index is built from the mesh.
        refine_levels = 0
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-ftcms_refine", &
                                refine_levels, refine_found, ierr)
        if (ierr /= 0) error stop "PETSc could not read -ftcms_refine"
        if (refine_levels > 0) then
            call DMPlexSetRefinementUniform(self%dm, PETSC_TRUE, ierr)
            do level = 1, refine_levels
                call DMRefine(self%dm, PETSC_COMM_WORLD, refined, ierr)
                if (ierr /= 0) error stop "DMPlex could not refine the mesh"
                call DMDestroy(self%dm, ierr)
                self%dm = refined
            end do
        end if

        call MPI_Comm_size(PETSC_COMM_WORLD, num_procs, ierr)
        if (num_procs > 1) then
            call DMPlexDistribute(self%dm, 0_PETSC_INT_KIND, PETSC_NULL_SF, distributed, ierr)
            if (ierr == 0 .and. distributed%v /= PETSC_NULL_DM%v) then
                call DMDestroy(self%dm, ierr)
                self%dm = distributed
            end if
        end if

        call DMGetDimension(self%dm, raw_dimension, ierr)
        self%dimension = int(raw_dimension, int32)

        call DMGetCoordinateSection(self%dm, section, ierr)
        call PetscSectionGetChart(section, self%coord_start, self%coord_end, ierr)

        call self%build_node_index()
        call self%build_global_node_numbering()
        call self%build_cell_index()
        call self%build_coloring()

        self%loaded = .true.
    end subroutine load_mesh_plex

    !> Stop on a mesh whose coordinates PETSc cannot hold continuously.
    !>
    !> PETSc reads a high-order simplex mesh with its coordinates held one set
    !> per cell, because the permutation from its own closure ordering to the
    !> Gmsh one is not implemented for simplices
    !> (src/dm/impls/plex/plexgmsh.c: "XXX FIXME Requires
    !> DMPlexSetClosurePermutationLexicographic()"). Tensor cells have that
    !> permutation and are read continuously.
    !>
    !> Forcing continuity through
    !> -dm_plex_gmsh_project_petscdualspace_lagrange_continuity overrides the
    !> guard rather than supplying the permutation: the node count comes out
    !> right but vertices and edges receive coefficients belonging to other
    !> cells, so the elements are silently wrong. A mesh that cannot be read
    !> correctly is refused here instead.
    subroutine reject_unsupported_coordinates(self, file_name)
        implicit none
        class(type_mesh_plex), intent(inout) :: self
        character(*), intent(in) :: file_name

        PetscErrorCode :: ierr
        PetscSection :: section
        PetscInt :: v_start, v_end, c_start, c_end, dof, coord_dim
        integer(int32) :: nodes_on_cell

        ! Every query is checked: a mesh this routine cannot inspect must not be
        ! waved through, since being waved through is what lets a broken one run.
        call DMGetCoordinateSection(self%dm, section, ierr)
        if (ierr /= 0) error stop "PETSc could not report where the mesh keeps its coordinates"
        call DMGetCoordinateDim(self%dm, coord_dim, ierr)
        if (ierr /= 0 .or. coord_dim <= 0) error stop "PETSc could not report the coordinate dimension"
        call DMPlexGetDepthStratum(self%dm, 0_PETSC_INT_KIND, v_start, v_end, ierr)
        if (ierr /= 0) error stop "PETSc could not report the mesh vertices"
        call DMPlexGetHeightStratum(self%dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) error stop "PETSc could not report the mesh cells"
        if (c_end <= c_start .or. v_end <= v_start) return

        ! Coordinates on the vertices mean the set is continuous, which is what
        ! every supported mesh gives.
        call PetscSectionGetDof(section, v_start, dof, ierr)
        if (ierr /= 0) error stop "PETSc could not report the coordinates held on a vertex"
        if (dof > 0) return

        call PetscSectionGetDof(section, c_start, dof, ierr)
        if (ierr /= 0) error stop "PETSc could not report the coordinates held on a cell"
        nodes_on_cell = int(dof, int32) / int(coord_dim, int32)

        write (*, '(A)') "This mesh cannot be read with a continuous set of nodes:"
        write (*, '(A,A)') "   file            : ", trim(file_name)
        write (*, '(A,I0,A)') "   coordinates are held per cell, ", nodes_on_cell, " to a cell"
        write (*, '(A)') "   PETSc reads a high-order simplex mesh this way because the"
        write (*, '(A)') "   ordering it would need for a shared node set is not implemented"
        write (*, '(A)') "   for simplices in this version; tensor cells are read continuously."
        write (*, '(A)') "   Use a linear simplex mesh, or a quadrilateral or hexahedral mesh"
        write (*, '(A)') "   of any order, or a PETSc that carries the simplex ordering."
        error stop "Unsupported mesh: high-order simplex"
    end subroutine reject_unsupported_coordinates

    subroutine destroy_mesh_plex(self)
        implicit none
        class(type_mesh_plex), intent(inout) :: self
        PetscErrorCode :: ierr

        if (.not. self%loaded) return
        if (self%output_dm_ready) then
            call DMDestroy(self%output_dm, ierr)
            self%output_dm_ready = .false.
        end if
        if (self%halo_ready) then
            call VecScatterDestroy(self%halo_scatter, ierr)
            call VecDestroy(self%halo_local, ierr)
            call VecDestroy(self%halo_global, ierr)
            self%halo_ready = .false.
        end if
        call DMDestroy(self%dm, ierr)
        self%loaded = .false.
    end subroutine destroy_mesh_plex

    function get_dm_mesh_plex(self) result(dm)
        implicit none
        class(type_mesh_plex), intent(in) :: self
        DM :: dm
        dm = self%dm
    end function get_dm_mesh_plex

    !> A node is a point that carries coordinates: the vertices, plus the edge
    !> and cell points of a higher-order mesh. Numbering them in DM point order
    !> keeps the numbering reproducible.
    subroutine build_node_index(self)
        implicit none
        class(type_mesh_plex), intent(inout) :: self

        PetscErrorCode :: ierr
        PetscSection :: section
        PetscInt :: point, dof, offset, coord_dim
        integer(int32) :: chart_size, index

        call DMGetCoordinateDim(self%dm, coord_dim, ierr)
        self%coordinate_dim = int(coord_dim, int32)
        call DMGetCoordinateSection(self%dm, section, ierr)

        chart_size = int(self%coord_end - self%coord_start, int32)
        if (allocated(self%point_first_node)) deallocate (self%point_first_node)
        if (allocated(self%point_num_nodes)) deallocate (self%point_num_nodes)
        allocate (self%point_first_node(max(chart_size, 1)))
        allocate (self%point_num_nodes(max(chart_size, 1)))
        self%point_first_node = 0
        self%point_num_nodes = 0

        ! The coordinate values are stored one node after another, so a node's
        ! index is its position in that storage and its coordinates are read
        ! straight out of it.
        self%num_nodes = 0
        do point = self%coord_start, self%coord_end - 1
            call PetscSectionGetDof(section, point, dof, ierr)
            if (dof <= 0) cycle
            call PetscSectionGetOffset(section, point, offset, ierr)
            index = int(point - self%coord_start, int32) + 1
            self%point_first_node(index) = int(offset, int32) / self%coordinate_dim + 1
            self%point_num_nodes(index) = int(dof, int32) / self%coordinate_dim
            self%num_nodes = self%num_nodes + self%point_num_nodes(index)
        end do

        if (self%node_section_ready) then
            call PetscSectionDestroy(self%node_section, ierr)
            self%node_section_ready = .false.
        end if
        call PetscSectionCreate(PETSC_COMM_WORLD, self%node_section, ierr)
        call PetscSectionSetChart(self%node_section, self%coord_start, self%coord_end, ierr)
        do point = self%coord_start, self%coord_end - 1
            index = int(point - self%coord_start, int32) + 1
            if (self%point_num_nodes(index) <= 0) cycle
            call PetscSectionSetDof(self%node_section, point, &
                                    int(self%point_num_nodes(index), PETSC_INT_KIND), ierr)
        end do
        call PetscSectionSetUp(self%node_section, ierr)
        self%node_section_ready = .true.
    end subroutine build_node_index

    !> Contiguous global node ids. A section with one dof per node, pushed
    !> through the point SF, is what makes the numbering agree across ranks:
    !> PETSc hands back a non-negative offset on the owning rank and the
    !> negated owner offset on the ranks that only see the node as a ghost.
    subroutine build_global_node_numbering(self)
        implicit none
        class(type_mesh_plex), intent(inout) :: self

        PetscErrorCode :: ierr
        PetscSF :: point_sf
        PetscSection :: local_section, global_section
        PetscInt :: point, offset
        integer(int32) :: node, num_owned, index, local
        integer(int32) :: send_buffer(1), recv_buffer(1)

        if (allocated(self%node_global_id)) deallocate (self%node_global_id)
        if (allocated(self%node_is_owned)) deallocate (self%node_is_owned)
        allocate (self%node_global_id(max(self%num_nodes, 1)))
        allocate (self%node_is_owned(max(self%num_nodes, 1)))
        self%node_global_id = 0
        self%node_is_owned = .true.

        call DMGetPointSF(self%dm, point_sf, ierr)

        call PetscSectionCreate(PETSC_COMM_WORLD, local_section, ierr)
        call PetscSectionSetChart(local_section, self%coord_start, self%coord_end, ierr)
        do index = 1, size(self%point_num_nodes)
            if (self%point_num_nodes(index) <= 0) cycle
            call PetscSectionSetDof(local_section, &
                                    int(index - 1, PETSC_INT_KIND) + self%coord_start, &
                                    int(self%point_num_nodes(index), PETSC_INT_KIND), ierr)
        end do
        call PetscSectionSetUp(local_section, ierr)

        call PetscSectionCreateGlobalSection(local_section, point_sf, PETSC_FALSE, PETSC_FALSE, &
                                             PETSC_FALSE, global_section, ierr)
        if (ierr /= 0) error stop "PETSc could not build the global node numbering"

        num_owned = 0
        do index = 1, size(self%point_num_nodes)
            if (self%point_num_nodes(index) <= 0) cycle
            point = int(index - 1, PETSC_INT_KIND) + self%coord_start
            call PetscSectionGetOffset(global_section, point, offset, ierr)
            do local = 1, self%point_num_nodes(index)
                node = self%point_first_node(index) + local - 1
                if (offset >= 0) then
                    self%node_global_id(node) = int(offset, int32) + local
                    self%node_is_owned(node) = .true.
                    num_owned = num_owned + 1
                else
                    ! An unowned point stores -(owner offset + 1).
                    self%node_global_id(node) = int(-offset, int32) + local - 1
                    self%node_is_owned(node) = .false.
                end if
            end do
        end do

        call PetscSectionDestroy(global_section, ierr)
        call PetscSectionDestroy(local_section, ierr)

        send_buffer(1) = num_owned
        recv_buffer(1) = num_owned
        call MPI_Allreduce(send_buffer, recv_buffer, 1, MPI_INTEGER, MPI_SUM, PETSC_COMM_WORLD, ierr)
        self%global_num_nodes = recv_buffer(1)
    end subroutine build_global_node_numbering

    !> The cells FTCMS iterates: every interior cell, then every face that
    !> carries a boundary patch id. Both keep their Gmsh physical tag.
    subroutine build_cell_index(self)
        implicit none
        class(type_mesh_plex), intent(inout) :: self

        PetscErrorCode :: ierr
        DMLabel :: face_label
        IS :: value_is, point_is
        PetscInt, pointer :: values(:), points(:)
        PetscInt :: c_start, c_end, num_values, num_points, point
        integer(int32) :: num_faces, index, i, j
        logical :: has_face_label

        call DMPlexGetHeightStratum(self%dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)

        nullify (values)
        nullify (points)
        num_values = 0
        num_faces = 0
        has_face_label = .false.
        call DMGetLabel(self%dm, "Face Sets", face_label, ierr)
        if (ierr == 0 .and. face_label%v /= PETSC_NULL_DMLABEL%v) then
            call DMLabelGetNumValues(face_label, num_values, ierr)
            if (num_values > 0) then
                has_face_label = .true.
                call DMLabelGetValueIS(face_label, value_is, ierr)
                call ISGetIndices(value_is, values, ierr)
                do i = 1, int(num_values, int32)
                    call DMLabelGetStratumSize(face_label, values(i), num_points, ierr)
                    num_faces = num_faces + int(num_points, int32)
                end do
            end if
        end if

        self%num_cells = int(c_end - c_start, int32) + num_faces
        call reallocate(self%cell_point, self%num_cells)
        call reallocate(self%cell_dimension, self%num_cells)
        call reallocate(self%cell_entity_id, self%num_cells)
        call reallocate(self%cell_fe_type, self%num_cells)
        call reallocate(self%cell_num_nodes, self%num_cells)
        call reallocate(self%cell_color, self%num_cells)

        index = 0
        do point = c_start, c_end - 1
            index = index + 1
            call describe_cell(self, point, index, "Cell Sets")
        end do

        if (has_face_label) then
            do i = 1, int(num_values, int32)
                call DMLabelGetStratumIS(face_label, values(i), point_is, ierr)
                if (point_is%v == PETSC_NULL_IS%v) cycle
                call ISGetSize(point_is, num_points, ierr)
                call ISGetIndices(point_is, points, ierr)
                do j = 1, int(num_points, int32)
                    index = index + 1
                    call describe_cell(self, points(j), index, "Face Sets")
                end do
                call ISRestoreIndices(point_is, points, ierr)
                call ISDestroy(point_is, ierr)
            end do
            call ISRestoreIndices(value_is, values, ierr)
            call ISDestroy(value_is, ierr)
        end if

        self%num_cells = index
    end subroutine build_cell_index

    !> Fill one entry of the cell index from the DM.
    subroutine describe_cell(self, point, index, label_name)
        implicit none
        class(type_mesh_plex), intent(inout) :: self
        PetscInt, intent(in) :: point
        integer(int32), intent(in) :: index
        character(*), intent(in) :: label_name

        PetscErrorCode :: ierr
        PetscInt :: entity
        integer(int32), allocatable :: nodes(:)
        integer(int32) :: count

        call get_cell_nodes(self, int(point, int32), nodes, count)

        self%cell_point(index) = int(point, int32)
        self%cell_num_nodes(index) = count
        call classify_cell(self, point, count, self%cell_fe_type(index), self%cell_dimension(index))

        call DMGetLabelValue(self%dm, label_name, point, entity, ierr)
        if (ierr /= 0 .or. entity < 0) then
            error stop "The mesh has a cell with no physical group; every region needs one"
        end if
        self%cell_entity_id(index) = int(entity, int32)

        if (allocated(nodes)) deallocate (nodes)
    end subroutine describe_cell

    !> Map the PETSc cell shape and coordinate degree onto the element type
    !> FTCMS knows. Both come from PETSc: the shape from DMPlexGetCellType and
    !> the degree from the space PETSc keeps the coordinates in, so no node
    !> count is counted here. The list is a whitelist: a shape or a degree that
    !> is not on it is rejected rather than approximated by a neighbour.
    subroutine classify_cell(self, point, num_nodes, element_type, dimension)
        implicit none
        class(type_mesh_plex), intent(in) :: self
        PetscInt, intent(in) :: point
        integer(int32), intent(in) :: num_nodes
        integer(int32), intent(inout) :: element_type, dimension

        PetscErrorCode :: ierr
        DMPolytopeType :: polytope
        integer(int32) :: order

        call DMPlexGetCellType(self%dm, point, polytope, ierr)
        if (ierr /= 0) error stop "PETSc could not report the shape of a cell"
        order = self%coordinate_order()

        element_type = -1
        dimension = -1

        if (polytope%v == DM_POLYTOPE_POINT%v) then
            dimension = 0
            element_type = FE_TYPE%VERTEX%ID
        else if (polytope%v == DM_POLYTOPE_SEGMENT%v) then
            dimension = 1
            if (order == 1) element_type = FE_TYPE%LINE%ID
            if (order == 2) element_type = FE_TYPE%QUADRATIC_EDGE%ID
            if (order == 3) element_type = FE_TYPE%LAGRANGE_CURVE%ID
        else if (polytope%v == DM_POLYTOPE_TRIANGLE%v) then
            dimension = 2
            if (order == 1) element_type = FE_TYPE%TRIANGLE%ID
            if (order == 2) element_type = FE_TYPE%QUADRATIC_TRIANGLE%ID
            if (order == 3) element_type = FE_TYPE%LAGRANGE_TRIANGLE%ID
        else if (polytope%v == DM_POLYTOPE_QUADRILATERAL%v) then
            dimension = 2
            if (order == 1) element_type = FE_TYPE%QUAD%ID
            if (order == 2) element_type = FE_TYPE%BIQUADRATIC_QUAD%ID
            if (order == 3) element_type = FE_TYPE%LAGRANGE_QUADRILATERAL%ID
        else if (polytope%v == DM_POLYTOPE_TETRAHEDRON%v) then
            dimension = 3
            if (order == 1) element_type = FE_TYPE%TETRA%ID
            if (order == 2) element_type = FE_TYPE%QUADRATIC_TETRA%ID
            if (order == 3) element_type = FE_TYPE%LAGRANGE_TETRAHEDRON%ID
        else if (polytope%v == DM_POLYTOPE_HEXAHEDRON%v) then
            dimension = 3
            if (order == 1) element_type = FE_TYPE%HEXAHEDRON%ID
            if (order == 2) element_type = FE_TYPE%TRIQUADRATIC_HEXAHEDRON%ID
            if (order == 3) element_type = FE_TYPE%LAGRANGE_HEXAHEDRON%ID
        end if

        if (element_type < 0) then
            write (*, '(a,i0,a,i0,a,i0)') "Unsupported mesh element: polytope ", polytope%v, &
                " with ", num_nodes, " nodes, at DM point ", point
            error stop "Unsupported mesh element"
        end if
    end subroutine classify_cell

    !> Greedy colouring of the interior cells so threads assembling cells of one
    !> colour never touch the same node.
    subroutine build_coloring(self)
        implicit none
        class(type_mesh_plex), intent(inout) :: self

        integer(int32), allocatable :: node_offset(:), node_cell(:), node_fill(:)
        integer(int32), allocatable :: nodes(:), forbidden(:)
        integer(int32) :: cell, i, j, k, count, node, neighbour, color
        integer(int32) :: total_entries

        self%num_colors = 0
        if (self%num_cells == 0 .or. self%num_nodes == 0) return

        ! Node to cell adjacency in CSR form, over the interior cells only.
        allocate (node_offset(self%num_nodes + 1))
        node_offset = 0
        total_entries = 0
        do cell = 1, self%num_cells
            if (self%cell_dimension(cell) /= self%dimension) cycle
            call get_cell_nodes(self, self%cell_point(cell), nodes, count)
            do i = 1, count
                node = nodes(i)
                node_offset(node + 1) = node_offset(node + 1) + 1
                total_entries = total_entries + 1
            end do
        end do

        node_offset(1) = 1
        do i = 1, self%num_nodes
            node_offset(i + 1) = node_offset(i) + node_offset(i + 1)
        end do

        allocate (node_cell(max(total_entries, 1)))
        allocate (node_fill(self%num_nodes))
        node_fill = 0
        do cell = 1, self%num_cells
            if (self%cell_dimension(cell) /= self%dimension) cycle
            call get_cell_nodes(self, self%cell_point(cell), nodes, count)
            do i = 1, count
                node = nodes(i)
                node_cell(node_offset(node) + node_fill(node)) = cell
                node_fill(node) = node_fill(node) + 1
            end do
        end do

        allocate (forbidden(self%num_cells + 1))
        forbidden = 0

        do cell = 1, self%num_cells
            if (self%cell_dimension(cell) /= self%dimension) cycle
            call get_cell_nodes(self, self%cell_point(cell), nodes, count)
            do i = 1, count
                node = nodes(i)
                do j = node_offset(node), node_offset(node + 1) - 1
                    neighbour = node_cell(j)
                    if (self%cell_color(neighbour) > 0) forbidden(self%cell_color(neighbour)) = cell
                end do
            end do

            color = 0
            do k = 1, self%num_cells
                if (forbidden(k) /= cell) then
                    color = k
                    exit
                end if
            end do
            self%cell_color(cell) = color
            self%num_colors = max(self%num_colors, color)
        end do

        if (allocated(nodes)) deallocate (nodes)
    end subroutine build_coloring

    !> Hand the dof numbering to PETSc.
    !>
    !> A section with num_dofs_per_node dofs on every node point becomes the
    !> DM's layout, so the global numbering, the owner of each dof and the
    !> matrix preallocation all come from the same description. Closure
    !> adjacency is what a finite-element operator couples.
    subroutine create_dof_layout_mesh_plex(self, num_dofs_per_node, layout)
        implicit none
        class(type_mesh_plex), intent(inout) :: self
        integer(int32), intent(in) :: num_dofs_per_node
        type(type_dof_layout), intent(inout) :: layout

        PetscErrorCode :: ierr
        PetscSection :: section, global_section
        PetscInt :: point, offset
        integer(int32) :: node, index, local

        call PetscSectionCreate(PETSC_COMM_WORLD, section, ierr)
        call PetscSectionSetChart(section, self%coord_start, self%coord_end, ierr)
        do index = 1, size(self%point_num_nodes)
            if (self%point_num_nodes(index) <= 0) cycle
            call PetscSectionSetDof(section, int(index - 1, PETSC_INT_KIND) + self%coord_start, &
                                    int(self%point_num_nodes(index) * num_dofs_per_node, PETSC_INT_KIND), ierr)
        end do
        call PetscSectionSetUp(section, ierr)

        call DMSetLocalSection(self%dm, section, ierr)
        call PetscSectionDestroy(section, ierr)
        call DMSetBasicAdjacency(self%dm, PETSC_FALSE, PETSC_TRUE, ierr)

        call DMGetGlobalSection(self%dm, global_section, ierr)
        if (ierr /= 0) error stop "PETSc could not build the global dof numbering"

        layout%num_dofs_per_node = num_dofs_per_node
        if (allocated(layout%node_dof_base)) deallocate (layout%node_dof_base)
        allocate (layout%node_dof_base(max(self%num_nodes, 1)))

        do index = 1, size(self%point_num_nodes)
            if (self%point_num_nodes(index) <= 0) cycle
            point = int(index - 1, PETSC_INT_KIND) + self%coord_start
            call PetscSectionGetOffset(global_section, point, offset, ierr)
            do local = 1, self%point_num_nodes(index)
                node = self%point_first_node(index) + local - 1
                if (offset >= 0) then
                    layout%node_dof_base(node) = int(offset, int32) + (local - 1) * num_dofs_per_node
                else
                    ! A node this rank only sees as a ghost stores
                    ! -(owner offset + 1). Assembling into the owner's index is
                    ! what sums the two ranks' contributions.
                    layout%node_dof_base(node) = int(-offset, int32) - 1 + (local - 1) * num_dofs_per_node
                end if
            end do
        end do

        layout%num_local_dofs = count(self%node_is_owned) * num_dofs_per_node
        layout%num_global_dofs = self%global_num_nodes * num_dofs_per_node
    end subroutine create_dof_layout_mesh_plex

    !> A DM that carries one value per node, for writing fields out.
    !>
    !> It is a clone: setting a section on the solver's own DM would change the
    !> layout the matrix was built against. The clone shares the topology and
    !> the coordinates, so what it describes is the same mesh.
    subroutine get_output_dm_mesh_plex(self, dm)
        implicit none
        class(type_mesh_plex), intent(inout) :: self
        DM, intent(inout) :: dm

        PetscErrorCode :: ierr
        PetscFE :: fe
        PetscObject :: discretisation
        PetscBool :: is_simplex
        PetscInt :: dimension
        integer(int32) :: element_order

        if (.not. self%output_dm_ready) then
            call DMClone(self%dm, self%output_dm, ierr)
            if (ierr /= 0) error stop "PETSc could not clone the mesh for output"

            ! The clone carries a discretisation, not a hand-built section. The
            ! HDF5 writer reaches DMGetDS on its way through
            ! DMPlexInsertBoundaryValues, so a DM without a PetscDS cannot write
            ! a field at all, and the section PETSc derives from the element
            ! places one value on each node this mesh has.
            call self%get_element_shape(element_order, is_simplex)
            call DMGetDimension(self%output_dm, dimension, ierr)

            ! Degree one, whatever the mesh order is. The topology PETSc writes
            ! for visualisation lists cell corners in the vertex numbering, so
            ! the coordinates beside it have to be the vertices too. A
            ! higher-order field here would write every node's coordinate and
            ! leave the corner indices pointing at the wrong points.
            call PetscFECreateLagrange(PETSC_COMM_WORLD, dimension, 1_PETSC_INT_KIND, is_simplex, &
                                       1_PETSC_INT_KIND, PETSC_DETERMINE, fe, ierr)
            if (ierr /= 0) error stop "PETSc could not build the output discretisation"

            discretisation = PetscObjectCast(fe)
            call DMAddField(self%output_dm, PETSC_NULL_DMLABEL, discretisation, ierr)
            if (ierr == 0) call DMCreateDS(self%output_dm, ierr)
            if (ierr /= 0) error stop "PETSc could not set up the output discretisation"
            call PetscFEDestroy(fe, ierr)

            self%output_dm_ready = .true.
        end if

        dm = self%output_dm
    end subroutine get_output_dm_mesh_plex

    !> FTCMS node index of each vertex, in the order the output section lays
    !> its degrees of freedom out. A field written on vertices has to be
    !> gathered from our own numbering into that order.
    subroutine get_vertex_nodes_mesh_plex(self, nodes, count)
        implicit none
        class(type_mesh_plex), intent(in) :: self
        integer(int32), allocatable, intent(inout) :: nodes(:)
        integer(int32), intent(inout) :: count

        PetscErrorCode :: ierr
        PetscInt :: v_start, v_end, point
        integer(int32) :: index

        call DMPlexGetDepthStratum(self%dm, 0_PETSC_INT_KIND, v_start, v_end, ierr)

        count = int(v_end - v_start, int32)
        if (allocated(nodes)) then
            if (size(nodes) < count) deallocate (nodes)
        end if
        if (.not. allocated(nodes)) allocate (nodes(max(count, 1)))

        count = 0
        do point = v_start, v_end - 1
            if (point < self%coord_start .or. point >= self%coord_end) cycle
            index = int(point - self%coord_start, int32) + 1
            if (self%point_num_nodes(index) <= 0) cycle
            count = count + 1
            nodes(count) = self%point_first_node(index)
        end do
    end subroutine get_vertex_nodes_mesh_plex

    !> What the visualisation topology PETSc writes looks like: how many
    !> interior cells it holds, how many corners each has, and the XDMF name of
    !> that shape. PETSc lists the corners of each cell, not its higher-order
    !> nodes, so the corner count is the polytope's own.
    subroutine get_visualisation_shape_mesh_plex(self, num_interior_cells, num_corners, topology_name)
        implicit none
        class(type_mesh_plex), intent(in) :: self
        integer(int32), intent(inout) :: num_interior_cells
        integer(int32), intent(inout) :: num_corners
        character(:), allocatable, intent(inout) :: topology_name

        integer(int32) :: cell, nodes_per_cell

        num_interior_cells = 0
        nodes_per_cell = 0
        do cell = 1, self%num_cells
            if (self%cell_dimension(cell) /= self%dimension) cycle
            num_interior_cells = num_interior_cells + 1
            if (nodes_per_cell == 0) nodes_per_cell = self%cell_num_nodes(cell)
        end do

        num_corners = 0
        topology_name = "Mixed"
        select case (self%dimension)
        case (1)
            num_corners = 2
            topology_name = "Polyline"
        case (2)
            if (nodes_per_cell == 3 .or. nodes_per_cell == 6 .or. nodes_per_cell == 10) then
                num_corners = 3
                topology_name = "Triangle"
            else
                num_corners = 4
                topology_name = "Quadrilateral"
            end if
        case (3)
            if (nodes_per_cell == 4 .or. nodes_per_cell == 10) then
                num_corners = 4
                topology_name = "Tetrahedron"
            else
                num_corners = 8
                topology_name = "Hexahedron"
            end if
        end select
    end subroutine get_visualisation_shape_mesh_plex

    !> Polynomial order and shape of the interior elements, read off the mesh.
    !> The degree and the shape of the mesh, both as PETSc reports them.
    subroutine get_element_shape_mesh_plex(self, element_order, is_simplex)
        implicit none
        class(type_mesh_plex), intent(in) :: self
        integer(int32), intent(inout) :: element_order
        PetscBool, intent(inout) :: is_simplex

        PetscErrorCode :: ierr
        PetscInt :: c_start, c_end
        DMPolytopeType :: polytope

        element_order = self%coordinate_order()
        is_simplex = PETSC_FALSE

        call DMPlexGetHeightStratum(self%dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0 .or. c_end <= c_start) return

        call DMPlexGetCellType(self%dm, c_start, polytope, ierr)
        if (ierr /= 0) return
        if (polytope%v == DM_POLYTOPE_TRIANGLE%v .or. &
            polytope%v == DM_POLYTOPE_TETRAHEDRON%v) is_simplex = PETSC_TRUE
    end subroutine get_element_shape_mesh_plex

    !> Degree of the space PETSc keeps the mesh coordinates in.
    !>
    !> This is the sequence PETSc itself uses in
    !> DMGetCoordinateDegree_Internal(), which it does not expose to Fortran.
    !> A coordinate DM carrying no finite element describes a mesh of degree
    !> one, which is also the answer when any step cannot be taken.
    function coordinate_order_mesh_plex(self) result(order)
        implicit none
        class(type_mesh_plex), intent(in) :: self
        integer(int32) :: order

        PetscErrorCode :: ierr
        DM :: coordinate_dm
        PetscObject :: field
        PetscFE :: fe
        PetscSpace :: space
        PetscInt :: num_fields, degree, max_degree
        character(len=256) :: class_name

        order = 1

        call DMGetCoordinateDM(self%dm, coordinate_dm, ierr)
        if (ierr /= 0) return
        call DMGetNumFields(coordinate_dm, num_fields, ierr)
        if (ierr /= 0 .or. num_fields < 1) return

        call DMGetField(coordinate_dm, 0_PETSC_INT_KIND, PETSC_NULL_DMLABEL, field, ierr)
        if (ierr /= 0) return
        call PetscObjectGetClassName(field, class_name, ierr)
        if (ierr /= 0 .or. trim(class_name) /= "PetscFE") return

        fe%v = field%v
        call PetscFEGetBasisSpace(fe, space, ierr)
        if (ierr /= 0) return
        call PetscSpaceGetDegree(space, degree, max_degree, ierr)
        if (ierr /= 0 .or. degree < 1) return

        order = int(degree, int32)
    end function coordinate_order_mesh_plex

    !> A matrix with the layout, sparsity and parallel distribution the DM
    !> derives from the section installed by create_dof_layout.
    subroutine create_matrix_mesh_plex(self, matrix)
        implicit none
        class(type_mesh_plex), intent(inout) :: self
        Mat, intent(inout) :: matrix

        PetscErrorCode :: ierr

        call DMCreateMatrix(self%dm, matrix, ierr)
        if (ierr /= 0) error stop "PETSc could not create the system matrix"
        ! The element loop writes each node block once; a repeated global index
        ! must add, not overwrite.
        call MatSetOption(matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
        ! An entry outside the preallocated pattern is a real mismatch between
        ! what the DM thinks the operator couples and what the assembly writes,
        ! so it must fail loudly rather than silently malloc.
        call MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE, ierr)
        call MatSetOption(matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)
    end subroutine create_matrix_mesh_plex

    !> Add up the contributions every rank made to the nodes they share.
    !>
    !> A rank assembles only its own cells, so a node on a partition boundary
    !> carries a partial sum on each rank that holds it. Routing the local
    !> values through a distributed Vec with ADD_VALUES leaves the complete sum
    !> on the owner; scattering it back gives every rank the assembled value.
    !> Without this a nodal residual or control volume is short by whatever the
    !> neighbouring rank contributed, and no reduction can recover it.
    subroutine halo_sum_nodal_mesh_plex(self, values)
        implicit none
        class(type_mesh_plex), intent(inout) :: self
        real(real64), intent(inout) :: values(:)

        PetscErrorCode :: ierr
        PetscMPIInt :: num_procs
        PetscInt, allocatable :: global_index(:)
        PetscScalar, pointer :: raw(:)
        IS :: from_is
        integer(int32) :: i

        call MPI_Comm_size(PETSC_COMM_WORLD, num_procs, ierr)
        if (num_procs <= 1) return
        if (size(values) /= self%num_nodes) return

        allocate (global_index(self%num_nodes))
        do i = 1, self%num_nodes
            global_index(i) = int(self%node_global_id(i) - 1, PETSC_INT_KIND)
        end do

        if (.not. self%halo_ready) then
            call VecCreateMPI(PETSC_COMM_WORLD, int(count(self%node_is_owned), PETSC_INT_KIND), &
                              int(self%global_num_nodes, PETSC_INT_KIND), self%halo_global, ierr)
            if (ierr == 0) call VecCreateSeq(PETSC_COMM_SELF, int(self%num_nodes, PETSC_INT_KIND), &
                                             self%halo_local, ierr)
            if (ierr == 0) call ISCreateGeneral(PETSC_COMM_SELF, int(self%num_nodes, PETSC_INT_KIND), &
                                                global_index, PETSC_COPY_VALUES, from_is, ierr)
            if (ierr == 0) call VecScatterCreate(self%halo_global, from_is, self%halo_local, &
                                                 PETSC_NULL_IS, self%halo_scatter, ierr)
            if (ierr == 0) call ISDestroy(from_is, ierr)
            if (ierr /= 0) then
                deallocate (global_index)
                return
            end if
            self%halo_ready = .true.
        end if

        call VecZeroEntries(self%halo_global, ierr)
        if (ierr == 0) call VecSetValues(self%halo_global, int(self%num_nodes, PETSC_INT_KIND), &
                                         global_index, values, ADD_VALUES, ierr)
        if (ierr == 0) call VecAssemblyBegin(self%halo_global, ierr)
        if (ierr == 0) call VecAssemblyEnd(self%halo_global, ierr)

        if (ierr == 0) call VecScatterBegin(self%halo_scatter, self%halo_global, self%halo_local, &
                                            INSERT_VALUES, SCATTER_FORWARD, ierr)
        if (ierr == 0) call VecScatterEnd(self%halo_scatter, self%halo_global, self%halo_local, &
                                          INSERT_VALUES, SCATTER_FORWARD, ierr)

        if (ierr == 0) then
            call VecGetArrayRead(self%halo_local, raw, ierr)
            do i = 1, self%num_nodes
                values(i) = real(raw(i), real64)
            end do
            call VecRestoreArrayRead(self%halo_local, raw, ierr)
        end if

        deallocate (global_index)
    end subroutine halo_sum_nodal_mesh_plex

    !> Coordinates of every node, in the dense node numbering, read through the
    !> coordinate section so a higher-order mesh reads correctly.
    subroutine get_node_coordinates_mesh_plex(self, points)
        implicit none
        class(type_mesh_plex), intent(in) :: self
        real(real64), intent(inout) :: points(:, :)

        PetscErrorCode :: ierr
        Vec :: coordinates
        PetscScalar, pointer :: values(:)
        integer(int32) :: node, k

        points = 0.0d0
        call DMGetCoordinatesLocal(self%dm, coordinates, ierr)

        ! The coordinate storage runs node by node, which is how the nodes were
        ! numbered, so this is a copy rather than a lookup.
        call VecGetArrayRead(coordinates, values, ierr)
        do node = 1, self%num_nodes
            do k = 1, min(self%coordinate_dim, size(points, 1))
                points(k, node) = real(values((node - 1) * self%coordinate_dim + k), real64)
            end do
        end do
        call VecRestoreArrayRead(coordinates, values, ierr)
    end subroutine get_node_coordinates_mesh_plex

    !> Connectivity of one cell in the dense node numbering.
    subroutine get_cell_connectivity_mesh_plex(self, cell, connectivity, count)
        implicit none
        class(type_mesh_plex), intent(in) :: self
        integer(int32), intent(in) :: cell
        integer(int32), allocatable, intent(inout) :: connectivity(:)
        integer(int32), intent(inout) :: count

        integer(int32), allocatable :: nodes(:)
        integer(int32) :: i

        call get_cell_nodes(self, self%cell_point(cell), nodes, count)

        if (allocated(connectivity)) then
            if (size(connectivity) < count) deallocate (connectivity)
        end if
        if (.not. allocated(connectivity)) allocate (connectivity(max(count, 1)))

        do i = 1, count
            connectivity(i) = nodes(i)
        end do

        if (allocated(nodes)) deallocate (nodes)
    end subroutine get_cell_connectivity_mesh_plex

    !> The distinct physical tags of the cells of one dimension. With no
    !> dimension given, the highest dimension present is used.
    subroutine get_active_region_info_mesh_plex(self, unique_ids, target_dim)
        implicit none
        class(type_mesh_plex), intent(in) :: self
        integer(int32), allocatable, intent(inout) :: unique_ids(:)
        integer(int32), intent(in), optional :: target_dim

        integer(int32), allocatable :: collected(:)
        integer(int32) :: wanted, cell, count

        if (allocated(unique_ids)) deallocate (unique_ids)

        if (present(target_dim)) then
            wanted = target_dim
        else
            wanted = 0
            do cell = 1, self%num_cells
                wanted = max(wanted, self%cell_dimension(cell))
            end do
            if (wanted <= 0) then
                allocate (unique_ids(0))
                return
            end if
        end if

        allocate (collected(max(self%num_cells, 1)))
        count = 0
        do cell = 1, self%num_cells
            if (self%cell_dimension(cell) /= wanted) cycle
            count = count + 1
            collected(count) = self%cell_entity_id(cell)
        end do

        if (count > 0) then
            call unique(collected(1:count), unique_ids)
        else
            allocate (unique_ids(0))
        end if
        deallocate (collected)
    end subroutine get_active_region_info_mesh_plex

    !> Nodes of one DM point: the closure points that carry coordinates, taken
    !> vertices first, then edges, then the interior, which is the order an
    !> element numbers its nodes in.
    subroutine get_cell_nodes(self, point, nodes, count)
        implicit none
        class(type_mesh_plex), intent(in) :: self
        integer(int32), intent(in) :: point
        integer(int32), allocatable, intent(inout) :: nodes(:)
        integer(int32), intent(inout) :: count

        PetscErrorCode :: ierr
        PetscInt, pointer :: closure(:)
        PetscInt :: num_closure, i
        integer(int32) :: index, local, capacity

        nullify (closure)
        call DMPlexGetTransitiveClosure(self%dm, int(point, PETSC_INT_KIND), PETSC_TRUE, &
                                        num_closure, closure, ierr)

        ! A closure point can carry several nodes, so the room needed is the
        ! number of values on the points, not the number of points.
        capacity = 0
        do i = 1, 2 * num_closure, 2
            if (closure(i) < self%coord_start .or. closure(i) >= self%coord_end) cycle
            capacity = capacity + self%point_num_nodes(int(closure(i) - self%coord_start, int32) + 1)
        end do
        if (allocated(nodes)) then
            if (size(nodes) < capacity) deallocate (nodes)
        end if
        if (.not. allocated(nodes)) allocate (nodes(max(capacity, 1)))

        ! Straight closure order, which is how PetscFE numbers its basis for the
        ! tensor cells. A point carrying several nodes is read in the direction
        ! the cell sees it, or two of its nodes would swap places.
        count = 0
        do i = 1, 2 * num_closure, 2
            if (closure(i) < self%coord_start .or. closure(i) >= self%coord_end) cycle
            index = int(closure(i) - self%coord_start, int32) + 1
            if (self%point_num_nodes(index) <= 0) cycle
            if (closure(i + 1) >= 0) then
                do local = 1, self%point_num_nodes(index)
                    count = count + 1
                    nodes(count) = self%point_first_node(index) + local - 1
                end do
            else
                do local = self%point_num_nodes(index), 1, -1
                    count = count + 1
                    nodes(count) = self%point_first_node(index) + local - 1
                end do
            end if
        end do

        call DMPlexRestoreTransitiveClosure(self%dm, int(point, PETSC_INT_KIND), PETSC_TRUE, &
                                            num_closure, closure, ierr)
    end subroutine get_cell_nodes

    subroutine reallocate(array, length)
        implicit none
        integer(int32), allocatable, intent(inout) :: array(:)
        integer(int32), intent(in) :: length

        if (allocated(array)) deallocate (array)
        allocate (array(max(length, 1)))
        array = 0
    end subroutine reallocate

end module domain_mesh_plex
