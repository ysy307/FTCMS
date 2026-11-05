#include <petsc/finclude/petscdmplex.h>

!> The mesh, held as a PETSc DMPlex.
!>
!> There is no second copy of the mesh here: the engine works from the DM PETSc
!> read, converted once to the analysis plane and distributed.
module domain_mesh_plex
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use petscdmplex
    use :: module_core, only:COMP_TYPES
    implicit none
    private

    public :: type_mesh_plex

    type :: type_mesh_plex
        DM :: dm
        logical :: loaded = .false.
        integer(int32) :: dimension = 0
        !> Vertices and cells of this rank's part of the mesh.
        integer(int32) :: num_nodes = 0
        integer(int32) :: num_cells = 0
    contains
        procedure, public :: load => load_mesh_plex
        procedure, public :: destroy => destroy_mesh_plex
        procedure, public :: get_dm => get_dm_mesh_plex
    end type type_mesh_plex

contains

    !> Read a Gmsh file and distribute it. PETSc chooses the partition and keeps
    !> the physical groups as the "Cell Sets" and "Face Sets" labels.
    subroutine load_mesh_plex(self, file_name, computation_type)
        implicit none
        class(type_mesh_plex), intent(inout) :: self
        character(*), intent(in) :: file_name
        !> `COMP_TYPES` id of the analysis; a 2D analysis gets its plane's two
        !> coordinates as an intrinsic two-dimensional mesh.
        integer(int32), intent(in) :: computation_type

        PetscErrorCode :: ierr
        DM :: distributed, refined
        PetscMPIInt :: num_procs
        PetscBool :: ready, refine_found, spacedim_given
        PetscInt :: raw_dimension, refine_levels, level, p_start, p_end

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
        call PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dm_plex_gmsh_spacedim", &
                                 spacedim_given, ierr)
        if (.not. (spacedim_given .eqv. PETSC_TRUE)) then
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-dm_plex_gmsh_spacedim", "3", ierr)
        end if
        if (ierr /= 0) error stop "PETSc rejected the mesh reader options"

        ! interpolate = TRUE builds the faces the boundary labels hang off.
        call DMPlexCreateFromFile(PETSC_COMM_WORLD, trim(file_name), "mesh", PETSC_TRUE, self%dm, ierr)
        if (ierr /= 0) error stop "DMPlex could not read the mesh"
        ! The Gmsh reader reads its options unprefixed; leave none behind.
        if (.not. (spacedim_given .eqv. PETSC_TRUE)) then
            call PetscOptionsClearValue(PETSC_NULL_OPTIONS, "-dm_plex_gmsh_spacedim", ierr)
        end if

        ! A mesh whose coordinates PETSc can only hold one set per cell has no
        ! shared node set to solve on, and is refused rather than bent into one.
        call reject_unsupported_coordinates(self, file_name)

        call convert_to_plane(self, computation_type)

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
            call choose_partitioner(self%dm)
            call DMPlexDistribute(self%dm, 0, PETSC_NULL_SF, distributed, ierr)
            if (ierr == 0 .and. distributed%v /= PETSC_NULL_DM%v) then
                call DMDestroy(self%dm, ierr)
                self%dm = distributed
            end if
        end if

        call DMGetDimension(self%dm, raw_dimension, ierr)
        self%dimension = raw_dimension
        call DMPlexGetDepthStratum(self%dm, 0, p_start, p_end, ierr)
        self%num_nodes = p_end - p_start
        call DMPlexGetHeightStratum(self%dm, 0, p_start, p_end, ierr)
        self%num_cells = p_end - p_start

        self%loaded = .true.
    end subroutine load_mesh_plex

    !> ParMETIS unless the run names another partitioner.
    !>
    !> `-petscpartitioner_type` still overrides it, and a PETSc built without
    !> ParMETIS keeps its own default.
    subroutine choose_partitioner(dm)
        implicit none
        DM, intent(inout) :: dm

        PetscPartitioner :: partitioner
        PetscErrorCode :: ierr

        call DMPlexGetPartitioner(dm, partitioner, ierr)
        if (ierr /= 0) return
        call PetscPartitionerSetType(partitioner, PETSCPARTITIONERPARMETIS, ierr)
        call PetscPartitionerSetFromOptions(partitioner, ierr)
    end subroutine choose_partitioner

    !> Rebuild a surface mesh read in three coordinates as an intrinsic 2D mesh.
    !>
    !> A 2D analysis is meshed in its own plane (x-y or x-z) and read with all
    !> three coordinates, so that an x-z mesh keeps its vertical coordinate. The
    !> discretisation needs a mesh whose coordinate dimension is its dimension:
    !> this keeps the plane's two coordinates, orients every cell
    !> counterclockwise in them, and carries "Cell Sets", "Face Sets" and
    !> "Vertex Sets" over. The mesh is still serial here, before distribution,
    !> so the cell list is built on the rank that read the file.
    !>
    !> Leaves a 3D analysis and an already planar mesh untouched. Fails by
    !> `error stop` on vertices off the plane, mixed cell shapes, a cell of no
    !> area, or a failing PETSc call.
    subroutine convert_to_plane(self, computation_type)
        implicit none
        class(type_mesh_plex), intent(inout) :: self
        integer(int32), intent(in) :: computation_type

        character(len=*), parameter :: LABELS(3) = [character(len=11) :: "Cell Sets", "Face Sets", "Vertex Sets"]
        DM :: plane
        Vec :: coordinates
        PetscSection :: coordinate_section
        PetscScalar, pointer :: xyz(:)
        PetscInt, pointer :: cone(:), edge(:), join(:)
        PetscInt, allocatable :: cell_vertices(:)
        PetscReal, allocatable :: points(:)
        PetscInt :: dimension, coordinate_dimension, c_start, c_end, v_start, v_end, f_start, f_end
        PetscInt :: new_v_start, new_v_end, point, offset, cone_size, join_size, value, corners, ends(2)
        PetscBool :: has_label
        PetscErrorCode :: ierr
        integer(int32) :: axes(2), dropped, k, found, listed(4), num_cells, num_vertices, l, labelled
        integer(int32) :: edge_vertices(2, 4)
        logical :: used(4)
        real(real64) :: area, reference, extent

        call DMGetDimension(self%dm, dimension, ierr)
        if (ierr == 0) call DMGetCoordinateDim(self%dm, coordinate_dimension, ierr)
        if (ierr /= 0) error stop "PETSc could not report the mesh dimensions"
        if (dimension /= 2 .or. coordinate_dimension /= 3) return

        select case (computation_type)
        case (COMP_TYPES%XY_2D%ID)
            axes = [1, 2]
            dropped = 3
        case (COMP_TYPES%XZ_2D%ID)
            axes = [1, 3]
            dropped = 2
        case default
            error stop "A surface mesh needs a two-dimensional computation type"
        end select

        call DMPlexGetHeightStratum(self%dm, 0, c_start, c_end, ierr)
        if (ierr == 0) call DMPlexGetDepthStratum(self%dm, 0, v_start, v_end, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(self%dm, 1, f_start, f_end, ierr)
        if (ierr == 0) call DMGetCoordinatesLocal(self%dm, coordinates, ierr)
        if (ierr == 0) call DMGetCoordinateSection(self%dm, coordinate_section, ierr)
        if (ierr /= 0) error stop "PETSc could not report the mesh points"
        num_cells = c_end - c_start
        num_vertices = v_end - v_start

        allocate (points(max(1, 2*num_vertices)))
        reference = 0.0d0
        extent = 0.0d0
        call VecGetArrayRead(coordinates, xyz, ierr)
        if (ierr /= 0) error stop "PETSc could not read the mesh coordinates"
        do point = v_start, v_end - 1
            call PetscSectionGetOffset(coordinate_section, point, offset, ierr)
            if (ierr /= 0) error stop "PETSc could not locate a vertex coordinate"
            k = point - v_start
            points(2*k + 1) = xyz(offset + axes(1))
            points(2*k + 2) = xyz(offset + axes(2))
            if (point == v_start) reference = xyz(offset + dropped)
            extent = max(extent, abs(xyz(offset + axes(1))), abs(xyz(offset + axes(2))))
        end do
        do point = v_start, v_end - 1
            call PetscSectionGetOffset(coordinate_section, point, offset, ierr)
            if (abs(xyz(offset + dropped) - reference) > 1.0d-9*max(extent, 1.0d0)) then
                error stop "The mesh has vertices off the plane of its computation type"
            end if
        end do
        call VecRestoreArrayRead(coordinates, xyz, ierr)

        corners = 0
        allocate (cell_vertices(max(1, 4*num_cells)))
        do point = c_start, c_end - 1
            ! Chain the cell's edges through their shared vertices, so the
            ! vertices come out in cyclic order without relying on how the
            ! closure or the edge orientations are ordered.
            call DMPlexGetConeSize(self%dm, point, cone_size, ierr)
            call DMPlexGetCone(self%dm, point, cone, ierr)
            if (ierr /= 0) error stop "PETSc could not read a cell cone"
            found = cone_size
            if (found /= 3 .and. found /= 4) error stop "The mesh has cells of neither three nor four edges"
            do k = 1, found
                call DMPlexGetCone(self%dm, cone(k), edge, ierr)
                edge_vertices(:, k) = edge(1:2) - v_start
                call DMPlexRestoreCone(self%dm, cone(k), edge, ierr)
            end do
            call DMPlexRestoreCone(self%dm, point, cone, ierr)
            used = .false.
            used(1) = .true.
            listed(1:2) = edge_vertices(:, 1)
            do k = 3, found
                do l = 2, found
                    if (used(l)) cycle
                    if (edge_vertices(1, l) == listed(k - 1)) then
                        listed(k) = edge_vertices(2, l)
                    else if (edge_vertices(2, l) == listed(k - 1)) then
                        listed(k) = edge_vertices(1, l)
                    else
                        cycle
                    end if
                    used(l) = .true.
                    exit
                end do
            end do
            if (corners == 0) corners = found
            if (found /= corners .or. (found /= 3 .and. found /= 4)) then
                error stop "The mesh mixes cell shapes or has cells of neither three nor four vertices"
            end if
            area = 0.0d0
            do k = 1, found
                l = mod(k, found) + 1
                area = area + points(2*listed(k) + 1)*points(2*listed(l) + 2) &
                       - points(2*listed(l) + 1)*points(2*listed(k) + 2)
            end do
            if (area == 0.0d0) error stop "The mesh has a cell of no area in its plane"
            if (area < 0.0d0) listed(2:found) = listed(found:2:-1)
            cell_vertices(corners*(point - c_start) + 1:corners*(point - c_start + 1)) = listed(1:found)
        end do

        ! Only the rank that read the file has cells; the others must still
        ! agree on the cell shape and on which labels exist.
        call MPI_Allreduce(MPI_IN_PLACE, corners, 1, MPIU_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)
        call DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, 2, num_cells, &
                                           num_vertices, corners, PETSC_TRUE, &
                                           cell_vertices, 2, points, plane, ierr)
        if (ierr /= 0) error stop "PETSc could not rebuild the mesh in its plane"
        call DMPlexGetDepthStratum(plane, 0, new_v_start, new_v_end, ierr)
        if (ierr /= 0 .or. new_v_end - new_v_start /= v_end - v_start) then
            error stop "The rebuilt mesh does not have the file's vertices"
        end if

        do l = 1, size(LABELS)
            call DMHasLabel(self%dm, trim(LABELS(l)), has_label, ierr)
            if (ierr /= 0) error stop "PETSc could not query a mesh label"
            labelled = merge(1_int32, 0_int32, has_label .eqv. PETSC_TRUE)
            call MPI_Allreduce(MPI_IN_PLACE, labelled, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)
            if (labelled == 0) cycle
            call DMCreateLabel(plane, trim(LABELS(l)), ierr)
            if (has_label .neqv. PETSC_TRUE) cycle
            select case (l)
            case (1)
                do point = c_start, c_end - 1
                    call DMGetLabelValue(self%dm, trim(LABELS(l)), point, value, ierr)
                    if (value >= 0) call DMSetLabelValue(plane, trim(LABELS(l)), point - c_start, value, ierr)
                end do
            case (2)
                do point = f_start, f_end - 1
                    call DMGetLabelValue(self%dm, trim(LABELS(l)), point, value, ierr)
                    if (value < 0) cycle
                    call DMPlexGetCone(self%dm, point, cone, ierr)
                    ends = cone(1:2) - v_start + new_v_start
                    call DMPlexRestoreCone(self%dm, point, cone, ierr)
                    call DMPlexGetJoin(plane, 2, ends, join_size, join, ierr)
                    if (ierr /= 0 .or. join_size /= 1) error stop "A labelled facet has no edge in the rebuilt mesh"
                    call DMSetLabelValue(plane, trim(LABELS(l)), join(1), value, ierr)
                    call DMPlexRestoreJoin(plane, 2, ends, join_size, join, ierr)
                end do
            case (3)
                do point = v_start, v_end - 1
                    call DMGetLabelValue(self%dm, trim(LABELS(l)), point, value, ierr)
                    if (value >= 0) call DMSetLabelValue(plane, trim(LABELS(l)), point - v_start + new_v_start, &
                                                         value, ierr)
                end do
            end select
            if (ierr /= 0) error stop "PETSc could not carry a label to the rebuilt mesh"
        end do

        call DMDestroy(self%dm, ierr)
        self%dm = plane
    end subroutine convert_to_plane

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
        call DMPlexGetDepthStratum(self%dm, 0, v_start, v_end, ierr)
        if (ierr /= 0) error stop "PETSc could not report the mesh vertices"
        call DMPlexGetHeightStratum(self%dm, 0, c_start, c_end, ierr)
        if (ierr /= 0) error stop "PETSc could not report the mesh cells"
        if (c_end <= c_start .or. v_end <= v_start) return

        ! Coordinates on the vertices mean the set is continuous, which is what
        ! every supported mesh gives.
        call PetscSectionGetDof(section, v_start, dof, ierr)
        if (ierr /= 0) error stop "PETSc could not report the coordinates held on a vertex"
        if (dof > 0) return

        call PetscSectionGetDof(section, c_start, dof, ierr)
        if (ierr /= 0) error stop "PETSc could not report the coordinates held on a cell"
        nodes_on_cell = dof/coord_dim

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
        call DMDestroy(self%dm, ierr)
        self%loaded = .false.
    end subroutine destroy_mesh_plex

    function get_dm_mesh_plex(self) result(dm)
        implicit none
        class(type_mesh_plex), intent(in) :: self
        DM :: dm
        dm = self%dm
    end function get_dm_mesh_plex

end module domain_mesh_plex
