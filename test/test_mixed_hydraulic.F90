#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscsnes.h>

!> The analytic fields the geometry case interpolates into the RT0/P0 spaces.
!>
!> They live in a module rather than inside the test program because PETSc
!> reaches them through a C function pointer, and an internal procedure has no
!> address a C caller can use.
module test_mixed_hydraulic_fields
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_double, c_ptr
    use :: petscsys
    use :: app_mixed_hydraulic_physics, only:MIXED_DIMENSION
    implicit none
    private

    public :: type_linear_field, projected_field
    public :: analytic_flux, analytic_pressure

    !> The field \(\mathbf q=\mathbf b+s\mathbf x\), whose divergence is
    !> \(d\,s\) in \(d\) dimensions and which lies exactly in RT0.
    type, bind(C) :: type_linear_field
        real(c_double) :: base(MIXED_DIMENSION)
        real(c_double) :: slope
        real(c_double) :: pressure
    end type type_linear_field

    type(type_linear_field), save :: projected_field

contains

    !> The flux \(\mathbf q=\mathbf b+s\mathbf x\) at a physical point.
    function analytic_flux(dim, time, x, num_components, u, ctx) &
        bind(C, name="ftcms_test_mixed_analytic_flux") result(ierr)
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        integer(int32) :: i

        do i = 1, MIXED_DIMENSION
            u(i) = real(projected_field%base(i) + projected_field%slope*real(x(i), real64), PETSC_REAL_KIND)
        end do
        ierr = 0
    end function analytic_flux

    !> The constant pressure that goes with it.
    function analytic_pressure(dim, time, x, num_components, u, ctx) &
        bind(C, name="ftcms_test_mixed_analytic_pressure") result(ierr)
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        u(1) = real(projected_field%pressure, PETSC_REAL_KIND)
        ierr = 0
    end function analytic_pressure

end module test_mixed_hydraulic_fields

!> Synthetic cases for the RT0/P0 mixed hydraulic solver.
!>
!> Each case is selected with `-mixed_case <name>` and needs no project
!> directory: the meshes are built here, so the geometry, the materials and the
!> boundaries of every check are visible in one file.
!>
!> The reference solutions are exact rather than refined. For a mesh whose
!> material interface lies on cell boundaries, the lowest-order mixed method
!> reproduces a piecewise-linear pressure's cell averages and a constant flux
!> exactly, because the divergence of an RT0 field is exactly the P0 space the
!> pressure is tested against. That makes a relative flux error of `1e-8` a
!> statement about the solver rather than about the mesh resolution.
program test_mixed_hydraulic
    use, intrinsic :: iso_fortran_env, only: int32, real64, output_unit
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_funptr, c_funloc, c_loc, c_f_pointer
    use :: petscdmplex
    use :: petscsnes
    use :: core_interop_petscds, only:petscds_project_function
    use :: core_interop_petscfe, only:petscfe_field_size
    use :: app_mixed_hydraulic, only:type_mixed_hydraulic, type_mixed_audit, &
                                     MIXED_BC_PRESSURE, MIXED_BC_FLUX
    use :: app_mixed_hydraulic_physics, only:type_mixed_material, MIXED_DIMENSION, &
                                             intrinsic_from_physical, physical_from_intrinsic
    use :: app_mixed_hydraulic_balance, only:MAX_CELL_FACETS, reconstruct_cell_facet_flux, &
                                             reconstruct_cell_divergence, reconstruct_cell_pressure, &
                                             audit_facet_continuity, cell_is_owned, cell_volume
    implicit none

    !> Standard gravity, m/s2. The conductivities of the reference cases were
    !> measured with it, so it is what turns them into pressure mobilities.
    real(real64), parameter :: STANDARD_GRAVITY = 9.80665d0
    !> Reference liquid density, kg/m3.
    real(real64), parameter :: REFERENCE_DENSITY = 1000.0d0

    !> The analytic field a projection case interpolates: q = base + slope*x.
    type, bind(C) :: type_linear_field
        real(c_double) :: base(MIXED_DIMENSION)
        real(c_double) :: slope
        real(c_double) :: pressure
    end type type_linear_field

    type(type_linear_field), save, target :: projected_field

    PetscErrorCode :: ierr
    character(len=64) :: case_name
    integer(int32) :: failures

    call PetscInitialize(ierr)
    if (ierr /= 0) error stop "PETSc could not start"

    failures = 0
    call read_case_name(case_name)

    select case (trim(case_name))
    case ("space")
        call case_space(failures)
    case ("darcy_constant")
        call case_darcy_constant(failures)
    case ("hydrostatic")
        call case_hydrostatic(failures)
    case ("layers")
        call case_layers(failures)
    case ("flux_boundary")
        call case_flux_boundary(failures)
    case ("jacobian")
        call case_jacobian(failures)
    case ("reject_mobility")
        call case_reject_mobility()
    case default
        call report_line("unknown case: "//trim(case_name))
        failures = 1
    end select

    if (failures == 0) then
        call report_line("RESULT pass  case "//trim(case_name))
    else
        call report_line("RESULT fail  case "//trim(case_name))
    end if

    call PetscFinalize(ierr)
    if (failures /= 0) error stop 1
contains

    !> The case this run was asked for, from `-mixed_case <name>`.
    subroutine read_case_name(name)
        implicit none
        character(len=*), intent(inout) :: name

        integer(int32) :: i, count
        character(len=64) :: argument

        name = "space"
        count = command_argument_count()
        do i = 1, count - 1
            call get_command_argument(i, argument)
            if (trim(argument) /= "-mixed_case") cycle
            call get_command_argument(i + 1, argument)
            name = trim(argument)
            return
        end do
    end subroutine read_case_name

    !> A real option read straight from the command line.
    function read_real_option(key, fallback) result(value)
        implicit none
        character(len=*), intent(in) :: key
        real(real64), intent(in) :: fallback
        real(real64) :: value

        integer(int32) :: i, count, status
        character(len=64) :: argument

        value = fallback
        count = command_argument_count()
        do i = 1, count - 1
            call get_command_argument(i, argument)
            if (trim(argument) /= key) cycle
            call get_command_argument(i + 1, argument)
            read (argument, *, iostat=status) value
            if (status /= 0) value = fallback
            return
        end do
    end function read_real_option

    !> Print on rank zero only.
    subroutine report_line(text)
        implicit none
        character(len=*), intent(in) :: text

        PetscMPIInt :: rank
        PetscErrorCode :: local_ierr

        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, local_ierr)
        if (rank == 0) write (output_unit, '(a)') text
        flush (output_unit)
    end subroutine report_line

    !> Record one measurement against its target.
    subroutine check(label, measured, tolerance, failures)
        implicit none
        character(len=*), intent(in) :: label
        real(real64), intent(in) :: measured
        real(real64), intent(in) :: tolerance
        integer(int32), intent(inout) :: failures

        character(len=160) :: line

        if (abs(measured) <= tolerance) then
            write (line, '(a,a,a,es12.5,a,es12.5)') "  ok   ", label, "  measured ", measured, "  target ", tolerance
        else
            write (line, '(a,a,a,es12.5,a,es12.5)') "  FAIL ", label, "  measured ", measured, "  target ", tolerance
            failures = failures + 1
        end if
        call report_line(trim(line))
    end subroutine check

    !> Record one integer measurement against its expected value.
    subroutine check_integer(label, measured, expected, failures)
        implicit none
        character(len=*), intent(in) :: label
        integer(int32), intent(in) :: measured
        integer(int32), intent(in) :: expected
        integer(int32), intent(inout) :: failures

        character(len=160) :: line

        if (measured == expected) then
            write (line, '(a,a,a,i0,a,i0)') "  ok   ", label, "  measured ", measured, "  expected ", expected
        else
            write (line, '(a,a,a,i0,a,i0)') "  FAIL ", label, "  measured ", measured, "  expected ", expected
            failures = failures + 1
        end if
        call report_line(trim(line))
    end subroutine check_integer

    ! ================================================================
    ! Mesh fixtures
    ! ================================================================

    !> Two skew affine triangles sharing one edge, in intrinsic coordinates.
    !>
    !> `reversed` lists the cells in the opposite order with their vertices
    !> rotated, which is the same geometry presented to PETSc differently. Every
    !> physical answer has to be the same on both, or the closure orientations
    !> are not being applied.
    subroutine build_two_triangles(reversed, space_dimension, dm, ierr)
        implicit none
        logical, intent(in) :: reversed
        !> Two for the intrinsic mesh, three for the planar (x,0,z) form.
        PetscInt, intent(in) :: space_dimension
        DM, intent(inout) :: dm
        PetscErrorCode, intent(inout) :: ierr

        PetscInt :: cells(6)
        PetscReal :: coordinates(12)
        real(real64) :: vertices(2, 4)
        integer(int32) :: v

        vertices(:, 1) = [0.0d0, 0.0d0]
        vertices(:, 2) = [1.3d0, 0.2d0]
        vertices(:, 3) = [0.4d0, 1.1d0]
        vertices(:, 4) = [1.7d0, 1.4d0]

        if (reversed) then
            ! The same two triangles, each rotated and listed the other way up.
            cells = [3, 1, 2, 2, 0, 1]
        else
            cells = [0, 1, 2, 1, 3, 2]
        end if

        coordinates = 0.0_PETSC_REAL_KIND
        do v = 1, 4
            if (space_dimension == 2) then
                coordinates((v - 1)*2 + 1) = real(vertices(1, v), PETSC_REAL_KIND)
                coordinates((v - 1)*2 + 2) = real(vertices(2, v), PETSC_REAL_KIND)
            else
                coordinates((v - 1)*3 + 1) = real(vertices(1, v), PETSC_REAL_KIND)
                coordinates((v - 1)*3 + 2) = 0.0_PETSC_REAL_KIND
                coordinates((v - 1)*3 + 3) = real(vertices(2, v), PETSC_REAL_KIND)
            end if
        end do

        call DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, 2_PETSC_INT_KIND, 2_PETSC_INT_KIND, &
                                           4_PETSC_INT_KIND, 3_PETSC_INT_KIND, PETSC_TRUE, cells, &
                                           space_dimension, coordinates, dm, ierr)
    end subroutine build_two_triangles

    !> A rectangle of right triangles, `nx` by `nz` cells.
    !>
    !> Each rectangle is cut along one diagonal, so the material interface at
    !> `x = 0.5` lies on cell boundaries whenever `nx` is even. That is what
    !> makes the layered reference solution exact rather than approximate.
    subroutine build_strip(nx, nz, width, height, dm, ierr)
        implicit none
        integer(int32), intent(in) :: nx
        integer(int32), intent(in) :: nz
        real(real64), intent(in) :: width
        real(real64), intent(in) :: height
        DM, intent(inout) :: dm
        PetscErrorCode, intent(inout) :: ierr

        PetscInt, allocatable :: cells(:)
        PetscReal, allocatable :: coordinates(:)
        integer(int32) :: i, k, corner, cell
        integer(int32) :: num_cells, num_vertices

        num_cells = 2*nx*nz
        num_vertices = (nx + 1)*(nz + 1)
        allocate (cells(3*num_cells), coordinates(2*num_vertices))

        do k = 0, nz
            do i = 0, nx
                corner = k*(nx + 1) + i
                coordinates(2*corner + 1) = real(width*real(i, real64)/real(nx, real64), PETSC_REAL_KIND)
                coordinates(2*corner + 2) = real(height*real(k, real64)/real(nz, real64), PETSC_REAL_KIND)
            end do
        end do

        cell = 0
        do k = 0, nz - 1
            do i = 0, nx - 1
                corner = k*(nx + 1) + i
                ! Counter-clockwise, which is the orientation DMPlex expects.
                cells(3*cell + 1) = corner
                cells(3*cell + 2) = corner + 1
                cells(3*cell + 3) = corner + nx + 2
                cell = cell + 1
                cells(3*cell + 1) = corner
                cells(3*cell + 2) = corner + nx + 2
                cells(3*cell + 3) = corner + nx + 1
                cell = cell + 1
            end do
        end do

        call DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, 2_PETSC_INT_KIND, &
                                           int(num_cells, PETSC_INT_KIND), &
                                           int(num_vertices, PETSC_INT_KIND), 3_PETSC_INT_KIND, &
                                           PETSC_TRUE, cells, 2_PETSC_INT_KIND, coordinates, dm, ierr)
        deallocate (cells, coordinates)
    end subroutine build_strip

    !> Give every cell a material, by the side of `split_x` its centroid is on.
    subroutine label_cells(dm, split_x, id_low, id_high, ierr)
        implicit none
        DM, intent(inout) :: dm
        real(real64), intent(in) :: split_x
        integer(int32), intent(in) :: id_low
        integer(int32), intent(in) :: id_high
        PetscErrorCode, intent(inout) :: ierr

        PetscInt :: cell, c_start, c_end
        PetscReal :: volume, centroid(3), normal(3)

        call DMCreateLabel(dm, "Cell Sets", ierr)
        if (ierr /= 0) return
        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return

        do cell = c_start, c_end - 1
            call DMPlexComputeCellGeometryFVM(dm, cell, volume, centroid, normal, ierr)
            if (ierr /= 0) return
            if (real(centroid(1), real64) < split_x) then
                call DMSetLabelValue(dm, "Cell Sets", cell, int(id_low, PETSC_INT_KIND), ierr)
            else
                call DMSetLabelValue(dm, "Cell Sets", cell, int(id_high, PETSC_INT_KIND), ierr)
            end if
            if (ierr /= 0) return
        end do
    end subroutine label_cells

    !> Mark the boundary facets whose centroid lies in the given box.
    subroutine label_facets(dm, value, x_low, x_high, z_low, z_high, ierr)
        implicit none
        DM, intent(inout) :: dm
        integer(int32), intent(in) :: value
        real(real64), intent(in) :: x_low, x_high, z_low, z_high
        PetscErrorCode, intent(inout) :: ierr

        PetscInt :: facet, f_start, f_end
        PetscInt, pointer :: support(:)
        PetscReal :: measure, centroid(3), normal(3)
        real(real64) :: x, z

        call DMCreateLabel(dm, "Face Sets", ierr)
        if (ierr /= 0) return
        call DMPlexGetHeightStratum(dm, 1_PETSC_INT_KIND, f_start, f_end, ierr)
        if (ierr /= 0) return

        do facet = f_start, f_end - 1
            call DMPlexGetSupport(dm, facet, support, ierr)
            if (ierr /= 0) return
            if (size(support) /= 1) cycle
            call DMPlexComputeCellGeometryFVM(dm, facet, measure, centroid, normal, ierr)
            if (ierr /= 0) return
            x = real(centroid(1), real64)
            z = real(centroid(2), real64)
            if (x < x_low .or. x > x_high .or. z < z_low .or. z > z_high) cycle
            call DMSetLabelValue(dm, "Face Sets", facet, int(value, PETSC_INT_KIND), ierr)
            if (ierr /= 0) return
        end do
    end subroutine label_facets

    !> Hand the mesh to every rank. The labels are already on it, so they travel.
    subroutine distribute(dm, ierr)
        implicit none
        DM, intent(inout) :: dm
        PetscErrorCode, intent(inout) :: ierr

        DM :: distributed
        PetscSF :: migration
        PetscMPIInt :: comm_size

        call MPI_Comm_size(PETSC_COMM_WORLD, comm_size, ierr)
        if (comm_size <= 1) return

        call DMPlexDistribute(dm, 0_PETSC_INT_KIND, migration, distributed, ierr)
        if (ierr /= 0) return
        call PetscSFDestroy(migration, ierr)
        call DMDestroy(dm, ierr)
        dm = distributed
    end subroutine distribute

    ! ================================================================
    ! Analytic fields for projection
    ! ================================================================

    !> The flux q = base + slope*x, whose divergence is `dim*slope`.
    !>
    !> This is exactly the lowest-order Raviart-Thomas space on a simplex, so
    !> its interpolant is itself and any departure measured afterwards is the
    !> mapping's, not the space's.
    function analytic_flux(dim, time, x, num_components, u, ctx) &
        bind(C, name="ftcms_test_mixed_analytic_flux") result(ierr)
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        integer(int32) :: i

        do i = 1, MIXED_DIMENSION
            u(i) = real(projected_field%base(i) + projected_field%slope*real(x(i), real64), PETSC_REAL_KIND)
        end do
        ierr = 0
    end function analytic_flux

    !> The constant pressure of a projection case.
    function analytic_pressure(dim, time, x, num_components, u, ctx) &
        bind(C, name="ftcms_test_mixed_analytic_pressure") result(ierr)
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        u(1) = real(projected_field%pressure, PETSC_REAL_KIND)
        ierr = 0
    end function analytic_pressure

    ! ================================================================
    ! Case: geometry, RT0/P0 construction and boundary flux data
    ! ================================================================

    subroutine case_space(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_material) :: materials(1)
        DM :: dm, dm_reversed, dm_planar
        real(real64) :: straight_flux, reversed_flux
        PetscMPIInt :: comm_size

        call MPI_Comm_size(PETSC_COMM_WORLD, comm_size, ierr)
        call report_line("case space: RT0/P0 construction, Piola geometry and facet data")

        materials(1)%id = 1
        materials(1)%conductivity = reshape([1.0d-5, 0.0d0, 0.0d0, 1.0d-5], [2, 2])
        materials(1)%storage = 0.0d0
        materials(1)%source = 0.0d0

        call prepare_two_triangle_mesh(.false., dm, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL could not build the two-triangle mesh")
            failures = failures + 1
            return
        end if

        call solver%set_materials(materials)
        call solver%set_gravity([0.0d0, 0.0d0], REFERENCE_DENSITY)
        call solver%set_scales(1.0d0, 1.0d0, 1.0d0)
        ! Face set 2 is the sloped lower edge; it carries a flux that is neither
        ! zero nor aligned with an axis, so a sign or a normal that was wrong
        ! would show up in its reconstructed trace.
        call solver%add_boundary(2, MIXED_BC_FLUX, flux=[3.0d-6, -1.0d-6])
        call solver%add_boundary(1, MIXED_BC_PRESSURE, pressure=0.0d0)
        call solver%initialize(dm, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the solver could not be initialised")
            failures = failures + 1
            return
        end if

        call check_element_layout(solver, failures)
        call check_projection(solver, failures, straight_flux)
        call check_prescribed_flux(solver, failures)

        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)

        ! The same geometry with the cells and their vertices permuted.
        call prepare_two_triangle_mesh(.true., dm_reversed, ierr)
        block
            type(type_mixed_hydraulic) :: reordered

            call reordered%set_materials(materials)
            call reordered%set_gravity([0.0d0, 0.0d0], REFERENCE_DENSITY)
            call reordered%set_scales(1.0d0, 1.0d0, 1.0d0)
            call reordered%add_boundary(2, MIXED_BC_FLUX, flux=[3.0d-6, -1.0d-6])
            call reordered%add_boundary(1, MIXED_BC_PRESSURE, pressure=0.0d0)
            call reordered%initialize(dm_reversed, ierr)
            call check_projection(reordered, failures, reversed_flux)
            call reordered%destroy(ierr)
        end block
        call DMDestroy(dm_reversed, ierr)
        call check("reordered boundary flux agrees          ", &
                   reversed_flux - straight_flux, 1.0d-15, failures)

        if (comm_size == 1) then
            ! The planar form of the same mesh, to show that carrying it as
            ! (x,z) is a change of representation and not of geometry.
            call build_two_triangles(.false., 3_PETSC_INT_KIND, dm_planar, ierr)
            call check_planar_mapping(dm_planar, failures)
            call DMDestroy(dm_planar, ierr)
        end if
    end subroutine case_space

    !> Build, label and distribute the two-triangle fixture.
    subroutine prepare_two_triangle_mesh(reversed, dm, ierr)
        implicit none
        logical, intent(in) :: reversed
        DM, intent(inout) :: dm
        PetscErrorCode, intent(inout) :: ierr

        call build_two_triangles(reversed, 2_PETSC_INT_KIND, dm, ierr)
        if (ierr /= 0) return
        call label_cells(dm, 1.0d30, 1, 1, ierr)
        if (ierr == 0) call label_facets(dm, 1, -1.0d0, 3.0d0, -1.0d0, 3.0d0, ierr)
        ! The sloped lower edge runs from (0,0) to (1.3,0.2); its centroid is
        ! the only boundary one below z = 0.2.
        if (ierr == 0) call label_facets(dm, 2, -1.0d0, 3.0d0, -1.0d0, 0.15d0, ierr)
        if (ierr == 0) call distribute(dm, ierr)
    end subroutine prepare_two_triangle_mesh

    !> The element really is RT0 and P0, and nothing else.
    subroutine check_element_layout(solver, failures)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        integer(int32), intent(inout) :: failures

        DM :: dm
        PetscSection :: section
        type(c_ptr) :: dm_handle
        integer(c_int) :: num_basis, num_components, cerr
        PetscInt :: point, v_start, v_end, e_start, e_end, c_start, c_end, dofs
        integer(int32) :: vertex_flux, edge_flux, cell_flux
        integer(int32) :: vertex_pressure, edge_pressure, cell_pressure

        dm = solver%get_dm()
        dm_handle = transfer(dm%v, dm_handle)

        cerr = 0
        num_basis = 0
        num_components = 0
        call petscfe_field_size(dm_handle, 0_c_int, num_basis, num_components, cerr)
        call check_integer("flux basis functions per cell           ", int(num_basis, int32), 3, failures)
        call check_integer("flux components                         ", int(num_components, int32), 2, failures)
        call petscfe_field_size(dm_handle, 1_c_int, num_basis, num_components, cerr)
        call check_integer("pressure basis functions per cell       ", int(num_basis, int32), 1, failures)
        call check_integer("pressure components                     ", int(num_components, int32), 1, failures)

        call DMGetLocalSection(dm, section, ierr)
        call DMPlexGetDepthStratum(dm, 0_PETSC_INT_KIND, v_start, v_end, ierr)
        call DMPlexGetHeightStratum(dm, 1_PETSC_INT_KIND, e_start, e_end, ierr)
        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)

        vertex_flux = 0
        edge_flux = 0
        cell_flux = 0
        vertex_pressure = 0
        edge_pressure = 0
        cell_pressure = 0
        do point = v_start, v_end - 1
            call PetscSectionGetFieldDof(section, point, 0_PETSC_INT_KIND, dofs, ierr)
            vertex_flux = vertex_flux + int(dofs, int32)
            call PetscSectionGetFieldDof(section, point, 1_PETSC_INT_KIND, dofs, ierr)
            vertex_pressure = vertex_pressure + int(dofs, int32)
        end do
        do point = e_start, e_end - 1
            call PetscSectionGetFieldDof(section, point, 0_PETSC_INT_KIND, dofs, ierr)
            edge_flux = max(edge_flux, int(dofs, int32))
            call PetscSectionGetFieldDof(section, point, 1_PETSC_INT_KIND, dofs, ierr)
            edge_pressure = edge_pressure + int(dofs, int32)
        end do
        do point = c_start, c_end - 1
            call PetscSectionGetFieldDof(section, point, 0_PETSC_INT_KIND, dofs, ierr)
            cell_flux = cell_flux + int(dofs, int32)
            call PetscSectionGetFieldDof(section, point, 1_PETSC_INT_KIND, dofs, ierr)
            cell_pressure = max(cell_pressure, int(dofs, int32))
        end do

        call check_integer("flux degrees of freedom on vertices     ", vertex_flux, 0, failures)
        call check_integer("flux degrees of freedom per facet       ", edge_flux, 1, failures)
        call check_integer("flux degrees of freedom in cells        ", cell_flux, 0, failures)
        call check_integer("pressure degrees of freedom on vertices ", vertex_pressure, 0, failures)
        call check_integer("pressure degrees of freedom on facets   ", edge_pressure, 0, failures)
        call check_integer("pressure degrees of freedom per cell    ", cell_pressure, 1, failures)
    end subroutine check_element_layout

    !> Interpolate an exactly representable flux and measure what comes back.
    subroutine check_projection(solver, failures, boundary_flux)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        integer(int32), intent(inout) :: failures
        !> The net outward flux the projected field carries, m2/s.
        real(real64), intent(inout) :: boundary_flux

        DM :: dm
        Vec :: local
        type(c_funptr) :: functions(2)
        type(c_ptr) :: dm_handle, solution_handle
        integer(c_int) :: cerr
        PetscInt :: cell, c_start, c_end
        PetscInt :: facets(MAX_CELL_FACETS)
        PetscInt, pointer :: support(:)
        real(real64) :: moments(MAX_CELL_FACETS)
        real(real64) :: divergence_integral, volume, outward_sum
        real(real64) :: worst_divergence, worst_constant, mismatch, total_area
        real(real64) :: expected_divergence
        integer(int32) :: num_facets, f, interior_facets
        real(real64) :: base(MIXED_DIMENSION), slope

        base = [2.0d-6, -5.0d-7]
        slope = 3.0d-6
        projected_field%base = real(base, c_double)
        projected_field%slope = real(slope, c_double)
        projected_field%pressure = 0.0_c_double

        dm = solver%get_dm()
        dm_handle = transfer(dm%v, dm_handle)
        solution_handle = transfer(solver%solution%v, solution_handle)

        functions(1) = c_funloc(analytic_flux)
        functions(2) = c_funloc(analytic_pressure)
        cerr = 0
        call petscds_project_function(dm_handle, 0.0_c_double, functions, 1_c_int, solution_handle, cerr)
        if (cerr /= 0) then
            call report_line("  FAIL PETSc could not interpolate the analytic flux")
            failures = failures + 1
            boundary_flux = 0.0d0
            return
        end if

        call solver%get_local_solution(local, ierr)
        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)

        worst_divergence = 0.0d0
        worst_constant = 0.0d0
        boundary_flux = 0.0d0
        total_area = 0.0d0
        expected_divergence = real(MIXED_DIMENSION, real64)*slope

        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(dm, cell)) cycle
            call reconstruct_cell_facet_flux(dm, local, solver%flux_scale, cell, facets, moments, num_facets, ierr)
            call reconstruct_cell_divergence(dm, local, solver%flux_scale, cell, divergence_integral, volume, ierr)
            outward_sum = sum(moments(1:num_facets))
            total_area = total_area + volume

            ! The divergence theorem, with the boundary integral and the volume
            ! integral built from the same degrees of freedom by two routes.
            worst_divergence = max(worst_divergence, abs(outward_sum - divergence_integral))
            ! And against the analytic divergence of the field that was put in.
            worst_constant = max(worst_constant, abs(divergence_integral - expected_divergence*volume))

            do f = 1, num_facets
                call DMPlexGetSupport(dm, facets(f), support, ierr)
                if (size(support) == 1) boundary_flux = boundary_flux + moments(f)
            end do
        end do

        call reduce_sum(boundary_flux)
        call reduce_sum(total_area)
        call reduce_max(worst_divergence)
        call reduce_max(worst_constant)

        call check("Piola divergence theorem defect         ", &
                   worst_divergence/(abs(expected_divergence)*total_area), 1.0d-12, failures)
        call check("interpolated divergence against exact   ", &
                   worst_constant/(abs(expected_divergence)*total_area), 1.0d-12, failures)
        call check("net boundary flux against exact         ", &
                   (boundary_flux - expected_divergence*total_area)/(abs(expected_divergence)*total_area), &
                   1.0d-12, failures)

        mismatch = 0.0d0
        interior_facets = 0
        call audit_facet_continuity(dm, local, solver%flux_scale, mismatch, interior_facets, ierr)
        call check_integer("interior facets seen from both sides    ", interior_facets, 1, failures)
        call check("facet normal trace mismatch             ", &
                   mismatch/(abs(expected_divergence)*total_area), 1.0d-12, failures)
    end subroutine check_projection

    !> The prescribed normal flux reaches the constrained facets intact.
    subroutine check_prescribed_flux(solver, failures)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        integer(int32), intent(inout) :: failures

        DM :: dm
        Vec :: local
        DMLabel :: face_sets
        PetscInt :: cell, c_start, c_end, facet_value
        PetscInt :: facets(MAX_CELL_FACETS)
        PetscReal :: measure, centroid(3), normal(3)
        real(real64) :: moments(MAX_CELL_FACETS)
        real(real64) :: prescribed(MIXED_DIMENSION), outward(MIXED_DIMENSION)
        real(real64) :: expected, worst, cell_measure, cell_centroid(3), orientation
        integer(int32) :: num_facets, f, d, marked
        PetscReal :: cell_volume_value, cell_centroid_r(3), cell_normal_r(3)

        prescribed = [3.0d-6, -1.0d-6]
        dm = solver%get_dm()

        ! Nothing has been solved; only the essential data is in the vector.
        call VecZeroEntries(solver%solution, ierr)
        call solver%get_local_solution(local, ierr)

        call DMGetLabel(dm, "Face Sets", face_sets, ierr)
        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)

        worst = 0.0d0
        marked = 0
        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(dm, cell)) cycle
            call reconstruct_cell_facet_flux(dm, local, solver%flux_scale, cell, facets, moments, num_facets, ierr)
            call DMPlexComputeCellGeometryFVM(dm, cell, cell_volume_value, cell_centroid_r, cell_normal_r, ierr)
            cell_measure = real(cell_volume_value, real64)
            cell_centroid = real(cell_centroid_r, real64)
            do f = 1, num_facets
                call DMGetLabelValue(dm, "Face Sets", facets(f), facet_value, ierr)
                if (ierr /= 0 .or. int(facet_value, int32) /= 2) cycle
                marked = marked + 1
                call DMPlexComputeCellGeometryFVM(dm, facets(f), measure, centroid, normal, ierr)
                orientation = 0.0d0
                do d = 1, MIXED_DIMENSION
                    orientation = orientation + real(normal(d), real64)*(real(centroid(d), real64) - cell_centroid(d))
                end do
                outward = real(normal(1:MIXED_DIMENSION), real64)
                if (orientation < 0.0d0) outward = -outward
                expected = real(measure, real64)*dot_product(prescribed, outward)
                worst = max(worst, abs(moments(f) - expected)/abs(expected))
            end do
        end do

        call reduce_max(worst)
        call reduce_count(marked)
        call check_integer("facets carrying the prescribed flux     ", marked, 1, failures)
        call check("prescribed normal flux, relative error  ", worst, 1.0d-12, failures)
    end subroutine check_prescribed_flux

    !> The planar (x,0,z) mesh and its intrinsic (x,z) form are one geometry.
    subroutine check_planar_mapping(dm_planar, failures)
        implicit none
        DM, intent(inout) :: dm_planar
        integer(int32), intent(inout) :: failures

        DM :: dm_intrinsic
        PetscInt :: cell, c_start, c_end, dimension, coordinate_dimension
        real(real64) :: planar_area, intrinsic_area, worst
        real(real64) :: gravity_physical(3), gravity_intrinsic(MIXED_DIMENSION)
        real(real64) :: flux_intrinsic(MIXED_DIMENSION), flux_physical(3)

        call DMGetDimension(dm_planar, dimension, ierr)
        call DMGetCoordinateDim(dm_planar, coordinate_dimension, ierr)
        call check_integer("planar mesh topological dimension       ", int(dimension, int32), 2, failures)
        call check_integer("planar mesh coordinate dimension        ", int(coordinate_dimension, int32), 3, failures)

        call build_two_triangles(.false., 2_PETSC_INT_KIND, dm_intrinsic, ierr)

        worst = 0.0d0
        call DMPlexGetHeightStratum(dm_planar, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        do cell = c_start, c_end - 1
            call cell_volume(dm_planar, cell, planar_area, ierr)
            call cell_volume(dm_intrinsic, cell, intrinsic_area, ierr)
            worst = max(worst, abs(planar_area - intrinsic_area)/abs(planar_area))
        end do
        call DMDestroy(dm_intrinsic, ierr)
        call check("intrinsic against planar cell areas     ", worst, 1.0d-14, failures)

        ! Gravity and a reconstructed flux travel the same map.
        gravity_physical = [0.0d0, 0.0d0, -STANDARD_GRAVITY]
        gravity_intrinsic = intrinsic_from_physical(gravity_physical)
        call check("gravity maps to intrinsic z             ", &
                   gravity_intrinsic(2) + STANDARD_GRAVITY, 1.0d-15, failures)
        call check("gravity has no intrinsic x component    ", gravity_intrinsic(1), 1.0d-15, failures)

        flux_intrinsic = [2.0d-6, -5.0d-7]
        flux_physical = physical_from_intrinsic(flux_intrinsic)
        call check("flux out-of-plane component is zero     ", flux_physical(2), 0.0d0, failures)
        call check("flux vertical component is intrinsic z  ", flux_physical(3) - flux_intrinsic(2), 0.0d0, failures)
    end subroutine check_planar_mapping

    ! ================================================================
    ! Case: uniform Darcy flow through a strip
    ! ================================================================

    subroutine case_darcy_constant(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        real(real64), parameter :: CONDUCTIVITY = 1.0d-5
        real(real64) :: head, exact_flux

        head = 1.0d0
        ! Series resistance of one layer is the layer itself.
        exact_flux = CONDUCTIVITY/1.0d0*head
        call report_line("case darcy_constant: uniform mobility, prescribed head difference")
        call run_layer_case(CONDUCTIVITY, CONDUCTIVITY, head, exact_flux, failures)
    end subroutine case_darcy_constant

    ! ================================================================
    ! Case: two layers in series
    ! ================================================================

    subroutine case_layers(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        real(real64) :: contrast, k_left, k_right, exact_flux, head
        character(len=96) :: line

        contrast = read_real_option("-mixed_contrast", 1.0d7)
        k_left = 1.0d-13
        k_right = k_left/contrast
        head = 1.0d0
        exact_flux = head/(0.5d0/k_left + 0.5d0/k_right)

        write (line, '(a,es9.2,a,es12.5)') "case layers: contrast ", contrast, ", exact flux ", exact_flux
        call report_line(trim(line))
        call run_layer_case(k_left, k_right, head, exact_flux, failures)

        call report_line("  reversed layers")
        call run_layer_case(k_right, k_left, head, exact_flux, failures)

        call report_line("  reversed head")
        call run_layer_case(k_left, k_right, -head, -exact_flux, failures)
    end subroutine case_layers

    !> Solve a two-layer strip and compare with the analytic series solution.
    !>
    !> \[
    !>   q=\frac{\Delta h}{0.5/K_1+0.5/K_2},\qquad
    !>   p(x)=\rho g\Big(\Delta h-\frac{q}{K(x)}x\Big)\ \text{piecewise linear}.
    !> \]
    subroutine run_layer_case(k_left, k_right, head, exact_flux, failures)
        implicit none
        real(real64), intent(in) :: k_left, k_right, head, exact_flux
        integer(int32), intent(inout) :: failures

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_material) :: materials(2)
        type(type_mixed_audit) :: report
        DM :: dm, solve_dm
        Vec :: local
        real(real64), parameter :: WIDTH = 1.0d0, HEIGHT = 0.25d0
        real(real64) :: left_pressure, pressure_scale, flux_scale
        real(real64) :: outflow, computed_flux, worst_pressure
        logical :: accepted

        materials(1)%id = 1
        materials(1)%conductivity = reshape([k_left, 0.0d0, 0.0d0, k_left], [2, 2])
        materials(2)%id = 2
        materials(2)%conductivity = reshape([k_right, 0.0d0, 0.0d0, k_right], [2, 2])

        left_pressure = REFERENCE_DENSITY*STANDARD_GRAVITY*head
        pressure_scale = abs(left_pressure)
        flux_scale = abs(exact_flux)

        call build_strip(8, 2, WIDTH, HEIGHT, dm, ierr)
        call label_cells(dm, 0.5d0, 1, 2, ierr)
        call label_facets(dm, 1, -0.01d0, 0.01d0, -1.0d0, 2.0d0, ierr)
        call label_facets(dm, 2, 0.99d0, 1.01d0, -1.0d0, 2.0d0, ierr)
        call label_facets(dm, 3, -1.0d0, 2.0d0, -0.01d0, 0.01d0, ierr)
        call label_facets(dm, 4, -1.0d0, 2.0d0, HEIGHT - 0.01d0, HEIGHT + 0.01d0, ierr)
        call distribute(dm, ierr)

        call solver%set_materials(materials)
        call solver%set_gravity([0.0d0, 0.0d0], REFERENCE_DENSITY)
        call solver%set_scales(flux_scale, pressure_scale, WIDTH)
        call solver%add_boundary(1, MIXED_BC_PRESSURE, pressure=left_pressure)
        call solver%add_boundary(2, MIXED_BC_PRESSURE, pressure=0.0d0)
        call solver%add_boundary(3, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
        call solver%add_boundary(4, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
        call solver%initialize(dm, ierr)

        accepted = .false.
        call solver%solve_step(0.0d0, 1.0d0, accepted, ierr)
        if (.not. accepted) then
            call report_line("  FAIL the layered solve did not converge")
            failures = failures + 1
            call solver%destroy(ierr)
            call DMDestroy(dm, ierr)
            return
        end if

        solve_dm = solver%get_dm()
        call solver%get_local_solution(local, ierr)

        call measure_boundary_outflow(solver, 2, outflow)
        computed_flux = outflow/HEIGHT
        call check("flux through the outlet, relative error ", &
                   (computed_flux - exact_flux)/abs(exact_flux), 1.0d-8, failures)

        call compare_layer_pressures(solver, k_left, k_right, exact_flux, left_pressure, worst_pressure)
        call check("cell pressure against cell average      ", worst_pressure, 1.0d-8, failures)

        call solver%compute_audit(report, ierr)
        call check("cell defect over throughput             ", &
                   report%max_cell_defect/report%throughput, 1.0d-10, failures)
        call check("domain defect over throughput           ", &
                   report%domain_defect/report%throughput, 1.0d-10, failures)
        call check("facet mismatch over throughput          ", &
                   report%max_facet_mismatch/report%throughput, 1.0d-12, failures)
        call check("divergence theorem over throughput      ", &
                   report%max_divergence_defect/report%throughput, 1.0d-12, failures)
        call check("boundary rate over throughput           ", &
                   report%boundary_rate/report%throughput, 1.0d-10, failures)

        block
            real(real64) :: mismatch
            integer(int32) :: interior_facets

            mismatch = 0.0d0
            interior_facets = 0
            call audit_facet_continuity(solve_dm, local, solver%flux_scale, mismatch, interior_facets, ierr)
            call check("facet continuity across all ranks       ", mismatch/report%throughput, 1.0d-12, failures)
        end block

        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)
    end subroutine run_layer_case

    !> The integrated outward flux through one labelled boundary, m2/s.
    subroutine measure_boundary_outflow(solver, boundary_id, outflow)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        integer(int32), intent(in) :: boundary_id
        real(real64), intent(inout) :: outflow

        DM :: dm
        Vec :: local
        PetscInt :: cell, c_start, c_end, facet_value
        PetscInt :: facets(MAX_CELL_FACETS)
        real(real64) :: moments(MAX_CELL_FACETS)
        integer(int32) :: num_facets, f

        outflow = 0.0d0
        dm = solver%get_dm()
        call solver%get_local_solution(local, ierr)
        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)

        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(dm, cell)) cycle
            call reconstruct_cell_facet_flux(dm, local, solver%flux_scale, cell, facets, moments, num_facets, ierr)
            do f = 1, num_facets
                call DMGetLabelValue(dm, "Face Sets", facets(f), facet_value, ierr)
                if (ierr /= 0) cycle
                if (int(facet_value, int32) /= boundary_id) cycle
                outflow = outflow + moments(f)
            end do
        end do
        call reduce_sum(outflow)
    end subroutine measure_boundary_outflow

    !> Compare the cell pressures with the cell averages of the exact profile.
    subroutine compare_layer_pressures(solver, k_left, k_right, exact_flux, left_pressure, worst)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        real(real64), intent(in) :: k_left, k_right, exact_flux, left_pressure
        real(real64), intent(inout) :: worst

        DM :: dm
        Vec :: local
        PetscInt :: cell, c_start, c_end
        PetscReal :: volume, centroid(3), normal(3)
        real(real64) :: pressure, expected, x

        worst = 0.0d0
        dm = solver%get_dm()
        call solver%get_local_solution(local, ierr)
        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)

        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(dm, cell)) cycle
            call reconstruct_cell_pressure(dm, local, solver%pressure_scale, cell, pressure, ierr)
            call DMPlexComputeCellGeometryFVM(dm, cell, volume, centroid, normal, ierr)
            ! The exact pressure is linear in each layer, so its cell average is
            ! its value at the centroid.
            x = real(centroid(1), real64)
            if (x < 0.5d0) then
                expected = left_pressure - REFERENCE_DENSITY*STANDARD_GRAVITY*exact_flux*x/k_left
            else
                expected = REFERENCE_DENSITY*STANDARD_GRAVITY*exact_flux*(1.0d0 - x)/k_right
            end if
            worst = max(worst, abs(pressure - expected)/abs(left_pressure))
        end do
        call reduce_max(worst)
    end subroutine compare_layer_pressures

    ! ================================================================
    ! Case: hydrostatic equilibrium
    ! ================================================================

    subroutine case_hydrostatic(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_material) :: materials(1)
        type(type_mixed_audit) :: report
        DM :: dm, solve_dm
        Vec :: local
        real(real64), parameter :: WIDTH = 0.4d0, HEIGHT = 1.0d0
        real(real64), parameter :: CONDUCTIVITY = 1.0d-6
        real(real64) :: base_pressure, mobility, flux_scale
        real(real64) :: worst_pressure, worst_flux, pressure, expected
        PetscInt :: cell, c_start, c_end
        PetscReal :: volume, centroid(3), normal(3)
        PetscInt :: facets(MAX_CELL_FACETS)
        real(real64) :: moments(MAX_CELL_FACETS)
        integer(int32) :: num_facets, f
        logical :: accepted

        call report_line("case hydrostatic: gravity balanced by the pressure gradient, no flow")

        materials(1)%id = 1
        materials(1)%conductivity = reshape([CONDUCTIVITY, 0.0d0, 0.0d0, CONDUCTIVITY], [2, 2])

        base_pressure = REFERENCE_DENSITY*STANDARD_GRAVITY*HEIGHT
        mobility = CONDUCTIVITY/(REFERENCE_DENSITY*STANDARD_GRAVITY)
        ! What a flux would be if the gravity term were not balanced at all.
        flux_scale = mobility*REFERENCE_DENSITY*STANDARD_GRAVITY

        call build_strip(4, 8, WIDTH, HEIGHT, dm, ierr)
        call label_cells(dm, 1.0d30, 1, 1, ierr)
        call label_facets(dm, 3, -1.0d0, 2.0d0, -0.01d0, 0.01d0, ierr)
        call label_facets(dm, 1, -0.01d0, 0.01d0, 0.005d0, 2.0d0, ierr)
        call label_facets(dm, 2, WIDTH - 0.01d0, WIDTH + 0.01d0, 0.005d0, 2.0d0, ierr)
        call label_facets(dm, 4, -1.0d0, 2.0d0, HEIGHT - 0.01d0, HEIGHT + 0.01d0, ierr)
        call distribute(dm, ierr)

        call solver%set_materials(materials)
        call solver%set_gravity([0.0d0, -STANDARD_GRAVITY], REFERENCE_DENSITY)
        call solver%set_scales(flux_scale, base_pressure, HEIGHT)
        call solver%add_boundary(3, MIXED_BC_PRESSURE, pressure=base_pressure)
        call solver%add_boundary(1, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
        call solver%add_boundary(2, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
        call solver%add_boundary(4, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
        call solver%initialize(dm, ierr)

        accepted = .false.
        call solver%solve_step(0.0d0, 1.0d0, accepted, ierr)
        if (.not. accepted) then
            call report_line("  FAIL the hydrostatic solve did not converge")
            failures = failures + 1
            call solver%destroy(ierr)
            call DMDestroy(dm, ierr)
            return
        end if

        solve_dm = solver%get_dm()
        call solver%get_local_solution(local, ierr)
        call DMPlexGetHeightStratum(solve_dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)

        worst_pressure = 0.0d0
        worst_flux = 0.0d0
        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(solve_dm, cell)) cycle
            call reconstruct_cell_pressure(solve_dm, local, solver%pressure_scale, cell, pressure, ierr)
            call DMPlexComputeCellGeometryFVM(solve_dm, cell, volume, centroid, normal, ierr)
            expected = base_pressure - REFERENCE_DENSITY*STANDARD_GRAVITY*real(centroid(2), real64)
            worst_pressure = max(worst_pressure, abs(pressure - expected)/base_pressure)

            call reconstruct_cell_facet_flux(solve_dm, local, solver%flux_scale, cell, facets, moments, num_facets, ierr)
            do f = 1, num_facets
                worst_flux = max(worst_flux, abs(moments(f))/(flux_scale*HEIGHT))
            end do
        end do
        call reduce_max(worst_pressure)
        call reduce_max(worst_flux)

        call check("hydrostatic pressure, relative error    ", worst_pressure, 1.0d-10, failures)
        call check("hydrostatic facet flux over its scale   ", worst_flux, 1.0d-10, failures)

        call solver%compute_audit(report, ierr)
        call check("hydrostatic domain defect               ", &
                   report%domain_defect/(flux_scale*HEIGHT), 1.0d-12, failures)

        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)
    end subroutine case_hydrostatic

    ! ================================================================
    ! Case: a prescribed non-zero normal flux drives the flow
    ! ================================================================

    subroutine case_flux_boundary(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_material) :: materials(1)
        type(type_mixed_audit) :: report
        DM :: dm, solve_dm
        Vec :: local
        real(real64), parameter :: WIDTH = 1.0d0, HEIGHT = 0.25d0
        real(real64), parameter :: CONDUCTIVITY = 1.0d-6
        real(real64) :: inflow_flux, mobility, expected_pressure, worst
        real(real64) :: inflow, outflow, pressure, x
        PetscInt :: cell, c_start, c_end
        PetscReal :: volume, centroid(3), normal(3)
        logical :: accepted

        call report_line("case flux_boundary: prescribed inflow, prescribed outlet pressure")

        inflow_flux = 2.0d-7
        mobility = CONDUCTIVITY/(REFERENCE_DENSITY*STANDARD_GRAVITY)
        materials(1)%id = 1
        materials(1)%conductivity = reshape([CONDUCTIVITY, 0.0d0, 0.0d0, CONDUCTIVITY], [2, 2])

        call build_strip(8, 2, WIDTH, HEIGHT, dm, ierr)
        call label_cells(dm, 1.0d30, 1, 1, ierr)
        call label_facets(dm, 1, -0.01d0, 0.01d0, -1.0d0, 2.0d0, ierr)
        call label_facets(dm, 2, WIDTH - 0.01d0, WIDTH + 0.01d0, -1.0d0, 2.0d0, ierr)
        call label_facets(dm, 3, 0.005d0, WIDTH - 0.005d0, -0.01d0, 0.01d0, ierr)
        call label_facets(dm, 4, 0.005d0, WIDTH - 0.005d0, HEIGHT - 0.01d0, HEIGHT + 0.01d0, ierr)
        call distribute(dm, ierr)

        call solver%set_materials(materials)
        call solver%set_gravity([0.0d0, 0.0d0], REFERENCE_DENSITY)
        call solver%set_scales(inflow_flux, inflow_flux*WIDTH/mobility, WIDTH)
        call solver%add_boundary(1, MIXED_BC_FLUX, flux=[inflow_flux, 0.0d0])
        call solver%add_boundary(2, MIXED_BC_PRESSURE, pressure=0.0d0)
        call solver%add_boundary(3, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
        call solver%add_boundary(4, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
        call solver%initialize(dm, ierr)

        accepted = .false.
        call solver%solve_step(0.0d0, 1.0d0, accepted, ierr)
        if (.not. accepted) then
            call report_line("  FAIL the prescribed-flux solve did not converge")
            failures = failures + 1
            call solver%destroy(ierr)
            call DMDestroy(dm, ierr)
            return
        end if

        call measure_boundary_outflow(solver, 1, inflow)
        call measure_boundary_outflow(solver, 2, outflow)
        call check("prescribed inflow reaches the boundary  ", &
                   (inflow + inflow_flux*HEIGHT)/(inflow_flux*HEIGHT), 1.0d-12, failures)
        call check("outflow matches inflow                  ", &
                   (outflow - inflow_flux*HEIGHT)/(inflow_flux*HEIGHT), 1.0d-10, failures)

        solve_dm = solver%get_dm()
        call solver%get_local_solution(local, ierr)
        call DMPlexGetHeightStratum(solve_dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        worst = 0.0d0
        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(solve_dm, cell)) cycle
            call reconstruct_cell_pressure(solve_dm, local, solver%pressure_scale, cell, pressure, ierr)
            call DMPlexComputeCellGeometryFVM(solve_dm, cell, volume, centroid, normal, ierr)
            x = real(centroid(1), real64)
            expected_pressure = inflow_flux*(WIDTH - x)/mobility
            worst = max(worst, abs(pressure - expected_pressure)/(inflow_flux*WIDTH/mobility))
        end do
        call reduce_max(worst)
        call check("pressure against the exact profile      ", worst, 1.0d-10, failures)

        call solver%compute_audit(report, ierr)
        call check("cell defect over throughput             ", &
                   report%max_cell_defect/report%throughput, 1.0d-10, failures)
        call check("boundary rate over throughput           ", &
                   report%boundary_rate/report%throughput, 1.0d-10, failures)

        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)
    end subroutine case_flux_boundary

    ! ================================================================
    ! Case: the assembled Jacobian against a central difference
    ! ================================================================

    subroutine case_jacobian(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_material) :: materials(2)
        DM :: dm
        Mat :: jacobian
        Vec :: state, direction, product, plus, minus, difference
        real(real64), parameter :: WIDTH = 1.0d0, HEIGHT = 0.25d0
        real(real64) :: step, worst, product_norm, difference_norm
        integer(int32) :: k
        PetscReal :: norm_value

        call report_line("case jacobian: assembled action against a central difference")

        materials(1)%id = 1
        materials(1)%conductivity = reshape([1.0d-6, 0.0d0, 0.0d0, 3.0d-6], [2, 2])
        materials(2)%id = 2
        materials(2)%conductivity = reshape([2.0d-7, 5.0d-8, 5.0d-8, 4.0d-7], [2, 2])

        call build_strip(6, 2, WIDTH, HEIGHT, dm, ierr)
        call label_cells(dm, 0.5d0, 1, 2, ierr)
        call label_facets(dm, 1, -0.01d0, 0.01d0, -1.0d0, 2.0d0, ierr)
        call label_facets(dm, 2, WIDTH - 0.01d0, WIDTH + 0.01d0, -1.0d0, 2.0d0, ierr)
        call label_facets(dm, 3, 0.005d0, WIDTH - 0.005d0, -0.01d0, 0.01d0, ierr)
        call label_facets(dm, 4, 0.005d0, WIDTH - 0.005d0, HEIGHT - 0.01d0, HEIGHT + 0.01d0, ierr)
        call distribute(dm, ierr)

        call solver%set_materials(materials)
        call solver%set_gravity([0.0d0, -STANDARD_GRAVITY], REFERENCE_DENSITY)
        call solver%set_scales(1.0d-7, 1.0d4, WIDTH)
        call solver%add_boundary(1, MIXED_BC_PRESSURE, pressure=1.0d4)
        call solver%add_boundary(2, MIXED_BC_PRESSURE, pressure=0.0d0)
        call solver%add_boundary(3, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
        call solver%add_boundary(4, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
        call solver%initialize(dm, ierr)

        call SNESSetUp(solver%snes, ierr)
        call DMCreateMatrix(solver%get_dm(), jacobian, ierr)

        call VecDuplicate(solver%solution, state, ierr)
        call VecDuplicate(solver%solution, direction, ierr)
        call VecDuplicate(solver%solution, product, ierr)
        call VecDuplicate(solver%solution, plus, ierr)
        call VecDuplicate(solver%solution, minus, ierr)
        call VecDuplicate(solver%solution, difference, ierr)

        ! A dimensionless state and a dimensionless direction, because that is
        ! what the residual was scaled for.
        call VecSetRandom(state, PETSC_NULL_RANDOM, ierr)
        call VecSetRandom(direction, PETSC_NULL_RANDOM, ierr)

        call SNESComputeJacobian(solver%snes, state, jacobian, jacobian, ierr)
        call MatMult(jacobian, direction, product, ierr)
        call VecNorm(product, NORM_2, norm_value, ierr)
        product_norm = real(norm_value, real64)

        worst = 1.0d30
        do k = 1, 4
            step = 10.0d0**(-real(k + 2, real64))
            call VecCopy(state, plus, ierr)
            call VecAXPY(plus, real(step, PETSC_REAL_KIND), direction, ierr)
            call VecCopy(state, minus, ierr)
            call VecAXPY(minus, real(-step, PETSC_REAL_KIND), direction, ierr)

            call SNESComputeFunction(solver%snes, plus, difference, ierr)
            call SNESComputeFunction(solver%snes, minus, minus, ierr)
            call VecAXPY(difference, -1.0_PETSC_REAL_KIND, minus, ierr)
            call VecScale(difference, real(1.0d0/(2.0d0*step), PETSC_REAL_KIND), ierr)
            call VecAXPY(difference, -1.0_PETSC_REAL_KIND, product, ierr)
            call VecNorm(difference, NORM_2, norm_value, ierr)
            difference_norm = real(norm_value, real64)
            worst = min(worst, difference_norm/product_norm)
        end do

        call check("Jacobian action against differences     ", worst, 1.0d-6, failures)

        call VecDestroy(state, ierr)
        call VecDestroy(direction, ierr)
        call VecDestroy(product, ierr)
        call VecDestroy(plus, ierr)
        call VecDestroy(minus, ierr)
        call VecDestroy(difference, ierr)
        call MatDestroy(jacobian, ierr)
        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)
    end subroutine case_jacobian

    ! ================================================================
    ! Case: a conductivity that has no inverse is refused
    ! ================================================================

    subroutine case_reject_mobility()
        implicit none

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_material) :: materials(1)
        DM :: dm

        call report_line("case reject_mobility: a zero conductivity has to stop the run")

        materials(1)%id = 1
        materials(1)%conductivity = 0.0d0

        call build_strip(2, 2, 1.0d0, 0.5d0, dm, ierr)
        call label_cells(dm, 1.0d30, 1, 1, ierr)
        call label_facets(dm, 1, -1.0d0, 2.0d0, -1.0d0, 2.0d0, ierr)
        call distribute(dm, ierr)

        call solver%set_materials(materials)
        call solver%set_gravity([0.0d0, 0.0d0], REFERENCE_DENSITY)
        call solver%set_scales(1.0d-6, 1.0d4, 1.0d0)
        call solver%add_boundary(1, MIXED_BC_PRESSURE, pressure=0.0d0)
        ! This must not return.
        call solver%initialize(dm, ierr)

        call report_line("  FAIL a zero conductivity was accepted")
        error stop 2
    end subroutine case_reject_mobility

    ! ================================================================
    ! Reductions
    ! ================================================================

    subroutine reduce_sum(value)
        implicit none
        real(real64), intent(inout) :: value

        real(real64) :: reduced
        integer(int32) :: mpi_err

        call MPI_Allreduce(value, reduced, 1, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, mpi_err)
        value = reduced
    end subroutine reduce_sum

    subroutine reduce_max(value)
        implicit none
        real(real64), intent(inout) :: value

        real(real64) :: reduced
        integer(int32) :: mpi_err

        call MPI_Allreduce(value, reduced, 1, MPI_DOUBLE_PRECISION, MPI_MAX, PETSC_COMM_WORLD, mpi_err)
        value = reduced
    end subroutine reduce_max

    subroutine reduce_count(value)
        implicit none
        integer(int32), intent(inout) :: value

        real(real64) :: sent, reduced
        integer(int32) :: mpi_err

        sent = real(value, real64)
        call MPI_Allreduce(sent, reduced, 1, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, mpi_err)
        value = nint(reduced)
    end subroutine reduce_count

end program test_mixed_hydraulic
