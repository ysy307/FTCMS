#include <petsc/finclude/petscdmplex.h>

!> Reconstruction of facet fluxes and cell divergences, and the water balance
!> they are audited against.
!>
!> Nothing here reads the assembled residual. Each quantity is rebuilt from the
!> solved degrees of freedom, PETSc's basis tabulation and PETSc's own
!> pushforward, so a balance that closes is evidence about the discretisation
!> rather than a restatement of the element loop that produced it.
!>
!> Both sides of an internal facet are reconstructed separately, each from its
!> own cell's closure and with that cell's own orientation flip. An orientation
!> error therefore shows up as a mismatch instead of cancelling.
!>
!> Every integral is per unit out-of-plane thickness: this is a planar domain,
!> so a cell "volume" is an area in m2 and a facet flux is in m2/s.
module app_mixed_hydraulic_balance
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
    use :: petscdmplex
    use :: core_interop_petscfe, only:petscfe_evaluate_cell, petscfe_reference_coordinates, &
                                      petscfe_quadrature
    use :: app_mixed_hydraulic_physics, only:MIXED_DIMENSION, MIXED_FIELD_FLUX, MIXED_FIELD_PRESSURE
    implicit none
    private

    public :: MAX_CELL_FACETS
    public :: MAX_QUADRATURE_POINTS
    public :: cell_quadrature_rule
    public :: reconstruct_cell_facet_flux
    public :: reconstruct_cell_divergence
    public :: reconstruct_cell_pressure
    public :: cell_is_owned
    public :: cell_volume
    public :: audit_facet_continuity

    !> The most facets one supported cell has: a triangle three, a rectangle four.
    integer(int32), parameter :: MAX_CELL_FACETS = 4

    ! Two-point Gauss on the unit interval, exact through degree three. The
    ! normal trace of an RT0 field is constant along a facet, but the rule is
    ! written for the trace of whatever field is handed in.
    integer(int32), parameter :: FACET_POINTS = 2
    real(real64), parameter :: FACET_NODES(FACET_POINTS) = &
        [0.5d0 - 0.5d0/sqrt(3.0d0), 0.5d0 + 0.5d0/sqrt(3.0d0)]
    real(real64), parameter :: FACET_WEIGHTS(FACET_POINTS) = [0.5d0, 0.5d0]

    ! How many facets each supported cell has, which is what tells the two
    ! reference cells apart. PETSc's reference triangle has vertices (-1,-1),
    ! (1,-1) and (-1,1); its reference quadrilateral is the square [-1,1]^2.
    integer(int32), parameter :: TRIANGLE_FACETS = 3
    integer(int32), parameter :: QUADRILATERAL_FACETS = 4

    !> The most quadrature points a rule this module will carry may have.
    integer(int32), parameter :: MAX_QUADRATURE_POINTS = 32

contains

    !> The rule one field of a DM is integrated with, on the reference cell.
    !>
    !> \[
    !>   \int_K f \approx |K|\sum_q w_q f(\xi_q),\qquad \sum_q w_q = 1,
    !> \]
    !>
    !> on an affine cell. The points and the point count are PETSc's own, read
    !> from the `PetscFE` the discrete system holds, so an audit integrating with
    !> this rule integrates what the element loop integrates rather than what the
    !> element was asked for. The weights are rescaled to sum to one so that the
    !> caller may use a measured cell volume and never the reference one.
    !>
    !> Assumes an affine cell. Exact for an integrand of degree `order` or less;
    !> no other error bound applies. \(O(1)\) arithmetic and memory.
    !>
    !> Fails by returning the PETSc error with `num_points = 0`, which includes
    !> a rule with more points than `MAX_QUADRATURE_POINTS`.
    subroutine cell_quadrature_rule(dm, field, points, weights, num_points, order, ierr)
        implicit none
        !> The DM whose discrete system holds the field.
        DM, intent(in) :: dm
        !> The field whose rule is asked for.
        PetscInt, intent(in) :: field
        !> The reference points, (dim, MAX_QUADRATURE_POINTS).
        real(c_double), intent(inout) :: points(MIXED_DIMENSION, MAX_QUADRATURE_POINTS)
        !> The weights, normalised to sum to one.
        real(c_double), intent(inout) :: weights(MAX_QUADRATURE_POINTS)
        !> How many points the rule has.
        integer(int32), intent(inout) :: num_points
        !> The polynomial degree the rule integrates exactly.
        integer(int32), intent(inout) :: order
        PetscErrorCode, intent(inout) :: ierr

        type(c_ptr) :: dm_handle
        integer(c_int) :: counted, rule_order, cerr
        real(real64) :: total

        points = 0.0d0
        weights = 0.0d0
        num_points = 0
        order = -1

        dm_handle = transfer(dm%v, dm_handle)
        counted = 0_c_int
        rule_order = -1_c_int
        cerr = 0
        call petscfe_quadrature(dm_handle, int(field, c_int), points, weights, counted, rule_order, cerr)
        if (cerr /= 0) then
            ierr = int(cerr)
            return
        end if

        num_points = int(counted, int32)
        order = int(rule_order, int32)
        total = sum(real(weights(1:num_points), real64))
        if (.not. (abs(total) > 0.0d0)) then
            num_points = 0
            ierr = 1
            return
        end if
        weights(1:num_points) = real(real(weights(1:num_points), real64)/total, c_double)
    end subroutine cell_quadrature_rule

    !> Whether this rank owns the cell, so that a ghost is counted once overall.
    function cell_is_owned(dm, cell) result(owned)
        implicit none
        !> The solve DM.
        DM, intent(in) :: dm
        !> A cell of the local chart.
        PetscInt, intent(in) :: cell
        logical :: owned

        PetscErrorCode :: ierr
        PetscInt :: global_start, global_end

        call DMPlexGetPointGlobal(dm, cell, global_start, global_end, ierr)
        owned = (ierr == 0) .and. (global_start >= 0)
    end function cell_is_owned

    !> The area of one cell, in m2 per unit thickness.
    subroutine cell_volume(dm, cell, volume, ierr)
        implicit none
        DM, intent(in) :: dm
        PetscInt, intent(in) :: cell
        !> The measured area, m2.
        real(real64), intent(inout) :: volume
        PetscErrorCode, intent(inout) :: ierr

        PetscReal :: measured, centroid(3), normal(3)

        volume = 0.0d0
        centroid = 0.0_PETSC_REAL_KIND
        normal = 0.0_PETSC_REAL_KIND
        call DMPlexComputeCellGeometryFVM(dm, cell, measured, centroid, normal, ierr)
        if (ierr == 0) volume = real(measured, real64)
    end subroutine cell_volume

    !> The outward integrated normal flux through each facet of one cell.
    !>
    !> \[
    !>   Q_{Kf}=\int_f \mathbf q_l\cdot\mathbf n_{Kf}\,\mathrm ds,
    !> \]
    !>
    !> with \(\mathbf n_{Kf}\) the outward unit normal of cell \(K\) on facet
    !> \(f\). Assumes an affine cell and intrinsic coordinates. The rule is
    !> exact for a normal trace of degree three or less; no other error bound
    !> applies. Arithmetic and memory are \(O(1)\) per facet.
    !>
    !> Fails by returning the PETSc error of the first failing query with
    !> `num_facets = 0`; the outputs are initialised on every path.
    subroutine reconstruct_cell_facet_flux(dm, local_x, flux_scale, cell, facets, moments, num_facets, ierr, field)
        implicit none
        !> The solve DM, carrying the flux field.
        DM, intent(in) :: dm
        !> A local vector holding the dimensionless solution.
        Vec, intent(in) :: local_x
        !> The flux unit of the unknown, m/s.
        real(real64), intent(in) :: flux_scale
        !> A cell of the local chart.
        PetscInt, intent(in) :: cell
        !> The cell's facets, in cone order.
        PetscInt, intent(inout) :: facets(MAX_CELL_FACETS)
        !> The outward flux through each facet, m2/s.
        real(real64), intent(inout) :: moments(MAX_CELL_FACETS)
        !> How many facets the cell has.
        integer(int32), intent(inout) :: num_facets
        PetscErrorCode, intent(inout) :: ierr
        !> The field carrying the flux; the isothermal layout's when absent.
        PetscInt, intent(in), optional :: field

        PetscInt :: flux_field
        PetscInt, pointer :: cone(:)
        PetscScalar, pointer :: facet_coordinates(:)
        PetscSection :: coordinate_section
        Vec :: coordinates
        DM :: coordinate_dm
        PetscReal :: facet_measure, facet_centroid(3), facet_normal(3)
        PetscReal :: cell_measure, cell_centroid(3), cell_normal(3)
        PetscInt :: closure_size
        type(c_ptr) :: dm_handle, vec_handle
        real(c_double) :: points(MIXED_DIMENSION, FACET_POINTS)
        real(c_double) :: reference(MIXED_DIMENSION, FACET_POINTS)
        real(c_double) :: values(MIXED_DIMENSION, FACET_POINTS)
        real(c_double) :: divergence(FACET_POINTS)
        real(real64) :: vertex_a(MIXED_DIMENSION), vertex_b(MIXED_DIMENSION)
        real(real64) :: outward(MIXED_DIMENSION), trace, orientation
        integer(int32) :: f, q, d
        integer(c_int) :: cerr

        facets = 0_PETSC_INT_KIND
        moments = 0.0d0
        num_facets = 0
        flux_field = MIXED_FIELD_FLUX
        if (present(field)) flux_field = field

        call DMPlexGetCone(dm, cell, cone, ierr)
        if (ierr /= 0) return
        num_facets = int(size(cone), int32)
        if (num_facets > MAX_CELL_FACETS) then
            num_facets = 0
            ierr = 1
            return
        end if
        facets(1:num_facets) = cone(1:num_facets)

        call DMPlexComputeCellGeometryFVM(dm, cell, cell_measure, cell_centroid, cell_normal, ierr)
        if (ierr /= 0) return

        call DMGetCoordinateDM(dm, coordinate_dm, ierr)
        if (ierr == 0) call DMGetCoordinatesLocal(dm, coordinates, ierr)
        if (ierr == 0) call DMGetCoordinateSection(dm, coordinate_section, ierr)
        if (ierr /= 0) return

        dm_handle = transfer(dm%v, dm_handle)
        vec_handle = transfer(local_x%v, vec_handle)

        do f = 1, num_facets
            call DMPlexComputeCellGeometryFVM(dm, facets(f), facet_measure, facet_centroid, facet_normal, ierr)
            if (ierr /= 0) return

            ! PETSc orients a facet normal by the facet, not by this cell. The
            ! direction from the cell's centroid to the facet's is what decides
            ! which way is out of this cell, and it is the same decision on both
            ! sides of an internal facet, taken independently.
            orientation = 0.0d0
            do d = 1, MIXED_DIMENSION
                orientation = orientation + real(facet_normal(d), real64)* &
                              (real(facet_centroid(d), real64) - real(cell_centroid(d), real64))
            end do
            outward = real(facet_normal(1:MIXED_DIMENSION), real64)
            if (orientation < 0.0d0) outward = -outward

            nullify (facet_coordinates)
            closure_size = 0_PETSC_INT_KIND
            call DMPlexVecGetClosure(coordinate_dm, coordinate_section, coordinates, facets(f), &
                                     closure_size, facet_coordinates, ierr)
            if (ierr /= 0) return
            if (closure_size /= 2*MIXED_DIMENSION) then
                call DMPlexVecRestoreClosure(coordinate_dm, coordinate_section, coordinates, facets(f), &
                                             closure_size, facet_coordinates, ierr)
                ierr = 1
                return
            end if
            do d = 1, MIXED_DIMENSION
                vertex_a(d) = real(facet_coordinates(d), real64)
                vertex_b(d) = real(facet_coordinates(MIXED_DIMENSION + d), real64)
            end do
            call DMPlexVecRestoreClosure(coordinate_dm, coordinate_section, coordinates, facets(f), &
                                         closure_size, facet_coordinates, ierr)
            if (ierr /= 0) return

            do q = 1, FACET_POINTS
                points(:, q) = real(vertex_a + FACET_NODES(q)*(vertex_b - vertex_a), c_double)
            end do

            cerr = 0
            call petscfe_reference_coordinates(dm_handle, int(cell, c_int), points, reference, cerr)
            if (cerr == 0) then
                call petscfe_evaluate_cell(dm_handle, vec_handle, int(flux_field, c_int), &
                                           int(cell, c_int), reference, values, divergence, cerr)
            end if
            if (cerr /= 0) then
                ierr = int(cerr)
                return
            end if

            trace = 0.0d0
            do q = 1, FACET_POINTS
                trace = trace + FACET_WEIGHTS(q)*dot_product(real(values(:, q), real64), outward)
            end do
            moments(f) = flux_scale*real(facet_measure, real64)*trace
        end do
    end subroutine reconstruct_cell_facet_flux

    !> The integrated divergence of the flux over one cell, m2/s.
    !>
    !> \[
    !>   \int_K\nabla\cdot\mathbf q_l\,\mathrm dV.
    !> \]
    !>
    !> The divergence comes from PETSc's pushforward of the reference basis
    !> gradients, which is the same map the assembly uses, evaluated here from
    !> the cell's own degrees of freedom. The points are the flux field's own
    !> rule, read from its `PetscFE`, so the rule belongs to the reference cell
    !> the field was built on rather than to one cell shape this module knows
    !> about. Assumes an affine cell. Exact for a divergence the field's own
    !> quadrature integrates exactly. Arithmetic and memory are \(O(1)\).
    !>
    !> Fails by returning the first PETSc error -- from the cell geometry, the
    !> field's quadrature or the tabulation -- with both outputs at zero.
    subroutine reconstruct_cell_divergence(dm, local_x, flux_scale, cell, integral, volume, ierr, field)
        implicit none
        DM, intent(in) :: dm
        Vec, intent(in) :: local_x
        !> The flux unit of the unknown, m/s.
        real(real64), intent(in) :: flux_scale
        PetscInt, intent(in) :: cell
        !> \(\int_K\nabla\cdot\mathbf q_l\,\mathrm dV\), m2/s.
        real(real64), intent(inout) :: integral
        !> The cell area, m2.
        real(real64), intent(inout) :: volume
        PetscErrorCode, intent(inout) :: ierr
        !> The field carrying the flux; the isothermal layout's when absent.
        PetscInt, intent(in), optional :: field

        PetscInt :: flux_field
        type(c_ptr) :: dm_handle, vec_handle
        real(c_double) :: points(MIXED_DIMENSION, MAX_QUADRATURE_POINTS)
        real(c_double) :: weights(MAX_QUADRATURE_POINTS)
        real(c_double) :: values(MIXED_DIMENSION, MAX_QUADRATURE_POINTS)
        real(c_double) :: divergence(MAX_QUADRATURE_POINTS)
        integer(int32) :: num_points, order, q
        real(real64) :: mean
        integer(c_int) :: cerr

        integral = 0.0d0
        volume = 0.0d0
        flux_field = MIXED_FIELD_FLUX
        if (present(field)) flux_field = field

        call cell_volume(dm, cell, volume, ierr)
        if (ierr /= 0) return
        call cell_quadrature_rule(dm, flux_field, points, weights, num_points, order, ierr)
        if (ierr /= 0) return

        dm_handle = transfer(dm%v, dm_handle)
        vec_handle = transfer(local_x%v, vec_handle)

        cerr = 0
        call petscfe_evaluate_cell(dm_handle, vec_handle, int(flux_field, c_int), &
                                   int(cell, c_int), points(:, :num_points), &
                                   values(:, :num_points), divergence(:num_points), cerr)
        if (cerr /= 0) then
            ierr = int(cerr)
            return
        end if

        ! The weights sum to one, so the measured area carries the scaling and
        ! the reference cell's own volume never enters.
        mean = 0.0d0
        do q = 1, num_points
            mean = mean + real(weights(q), real64)*real(divergence(q), real64)
        end do
        integral = flux_scale*volume*mean
    end subroutine reconstruct_cell_divergence

    !> The cell's pressure, Pa.
    !>
    !> A P0 field is constant over the cell, so any interior point answers. The
    !> point used is the centroid of the reference cell, which is also the
    !> preimage of the physical centroid under an affine map -- so a caller may
    !> pass a field that is not cell-constant and still get the value at the
    !> cell's centre. That is what makes the reference point cell-shape
    !> dependent rather than a constant of this module.
    !>
    !> Fails with a PETSc error for a cell whose reference centroid is unknown,
    !> which is any cell that is neither a triangle nor a quadrilateral.
    subroutine reconstruct_cell_pressure(dm, local_x, pressure_scale, cell, pressure, ierr, field)
        implicit none
        DM, intent(in) :: dm
        Vec, intent(in) :: local_x
        !> The pressure unit of the unknown, Pa.
        real(real64), intent(in) :: pressure_scale
        PetscInt, intent(in) :: cell
        !> The cell pressure, Pa.
        real(real64), intent(inout) :: pressure
        PetscErrorCode, intent(inout) :: ierr
        !> The field carrying the pressure; the isothermal layout's when absent.
        PetscInt, intent(in), optional :: field

        PetscInt :: pressure_field
        type(c_ptr) :: dm_handle, vec_handle
        real(c_double) :: reference(MIXED_DIMENSION, 1)
        real(c_double) :: values(1, 1)
        real(c_double) :: divergence(1)
        integer(c_int) :: cerr

        pressure = 0.0d0
        ierr = 0
        pressure_field = MIXED_FIELD_PRESSURE
        if (present(field)) pressure_field = field

        call cell_reference_centroid(dm, cell, reference(:, 1), ierr)
        if (ierr /= 0) return

        dm_handle = transfer(dm%v, dm_handle)
        vec_handle = transfer(local_x%v, vec_handle)

        cerr = 0
        call petscfe_evaluate_cell(dm_handle, vec_handle, int(pressure_field, c_int), &
                                   int(cell, c_int), reference, values, divergence, cerr)
        if (cerr /= 0) then
            ierr = int(cerr)
            return
        end if
        pressure = pressure_scale*real(values(1, 1), real64)
    end subroutine reconstruct_cell_pressure

    !> The largest disagreement between the two sides of any facet.
    !>
    !> \[
    !>   \max_f\Big|\sum_{K\supset f}Q_{Kf}\Big|,
    !> \]
    !>
    !> where each \(Q_{Kf}\) is reconstructed from cell \(K\)'s own degrees of
    !> freedom with its own outward normal. Facets are matched by PETSc's global
    !> point numbering and the sums are reduced over the communicator, so the two
    !> cells of a facet split across a partition boundary are still compared with
    !> each other -- and compared as two separately computed numbers, not as one
    !> number and its negation.
    !>
    !> Assumes affine cells and that every cell is owned by exactly one rank.
    !> For a conforming H(div) space the exact answer is zero; what is measured
    !> is the reconstruction and arithmetic error. Arithmetic is \(O(N)\) in
    !> local cells, memory \(O(P)\) in the mesh's global point count.
    !>
    !> Fails by returning the first PETSc error with the outputs at zero.
    subroutine audit_facet_continuity(dm, local_x, flux_scale, max_mismatch, interior_facets, ierr, field)
        implicit none
        !> The solve DM.
        DM, intent(in) :: dm
        !> A local vector holding the dimensionless solution.
        Vec, intent(in) :: local_x
        !> The flux unit of the unknown, m/s.
        real(real64), intent(in) :: flux_scale
        !> The largest facet disagreement, m2/s.
        real(real64), intent(inout) :: max_mismatch
        !> How many facets were seen from two cells.
        integer(int32), intent(inout) :: interior_facets
        PetscErrorCode, intent(inout) :: ierr
        !> The field carrying the flux; the isothermal layout's when absent.
        PetscInt, intent(in), optional :: field

        PetscInt :: flux_field

        IS :: numbering
        PetscInt, pointer :: numbers(:)
        PetscInt :: cell, c_start, c_end, p_start, p_end, point
        PetscInt :: facets(MAX_CELL_FACETS)
        real(real64) :: moments(MAX_CELL_FACETS)
        real(real64), allocatable :: facet_sum(:), facet_sum_all(:)
        real(real64), allocatable :: facet_count(:), facet_count_all(:)
        integer(int32) :: num_facets, f, index, total_points, mpi_err
        real(real64) :: highest_send, highest

        max_mismatch = 0.0d0
        interior_facets = 0
        flux_field = MIXED_FIELD_FLUX
        if (present(field)) flux_field = field

        call DMPlexGetChart(dm, p_start, p_end, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr == 0) call DMPlexCreatePointNumbering(dm, numbering, ierr)
        if (ierr /= 0) return
        call ISGetIndices(numbering, numbers, ierr)
        if (ierr /= 0) return

        ! A point this rank does not own carries minus its global number minus
        ! one, so the same facet has the same index on both sides.
        highest_send = 0.0d0
        do point = p_start, p_end - 1
            index = global_index_of(numbers(point - p_start + 1))
            highest_send = max(highest_send, real(index, real64))
        end do
        call MPI_Allreduce(highest_send, highest, 1, MPI_DOUBLE_PRECISION, MPI_MAX, PETSC_COMM_WORLD, mpi_err)
        if (mpi_err /= 0) then
            call ISRestoreIndices(numbering, numbers, ierr)
            call ISDestroy(numbering, ierr)
            ierr = 1
            return
        end if
        total_points = int(highest, int32) + 1

        allocate (facet_sum(total_points), facet_sum_all(total_points))
        allocate (facet_count(total_points), facet_count_all(total_points))
        facet_sum = 0.0d0
        facet_count = 0.0d0

        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(dm, cell)) cycle
            call reconstruct_cell_facet_flux(dm, local_x, flux_scale, cell, facets, moments, num_facets, ierr, &
                                             field=flux_field)
            if (ierr /= 0) exit
            do f = 1, num_facets
                index = global_index_of(numbers(facets(f) - p_start + 1)) + 1
                facet_sum(index) = facet_sum(index) + moments(f)
                facet_count(index) = facet_count(index) + 1.0d0
            end do
        end do

        call ISRestoreIndices(numbering, numbers, ierr)
        call ISDestroy(numbering, ierr)

        call MPI_Allreduce(facet_sum, facet_sum_all, total_points, MPI_DOUBLE_PRECISION, MPI_SUM, &
                           PETSC_COMM_WORLD, mpi_err)
        if (mpi_err == 0) then
            call MPI_Allreduce(facet_count, facet_count_all, total_points, MPI_DOUBLE_PRECISION, MPI_SUM, &
                               PETSC_COMM_WORLD, mpi_err)
        end if
        if (mpi_err /= 0) then
            ierr = 1
            return
        end if

        do index = 1, total_points
            if (nint(facet_count_all(index)) /= 2) cycle
            interior_facets = interior_facets + 1
            max_mismatch = max(max_mismatch, abs(facet_sum_all(index)))
        end do

        deallocate (facet_sum, facet_sum_all, facet_count, facet_count_all)
    end subroutine audit_facet_continuity

    !> The centroid of one cell's reference cell, in reference coordinates.
    !>
    !> Which reference cell that is follows from the facet count, because the
    !> two shapes this module supports have different ones: PETSc's reference
    !> triangle is centred at \((-1/3,-1/3)\) and its reference square at the
    !> origin. Using the triangle's point on a quadrilateral evaluates a field
    !> a sixth of a cell away from the centre in each direction, which is
    !> invisible for a cell-constant field and wrong for every other one.
    subroutine cell_reference_centroid(dm, cell, reference, ierr)
        implicit none
        DM, intent(in) :: dm
        PetscInt, intent(in) :: cell
        !> The reference centroid.
        real(c_double), intent(inout) :: reference(MIXED_DIMENSION)
        PetscErrorCode, intent(inout) :: ierr

        PetscInt :: cone_size

        reference = 0.0_c_double
        call DMPlexGetConeSize(dm, cell, cone_size, ierr)
        if (ierr /= 0) return

        select case (int(cone_size, int32))
        case (TRIANGLE_FACETS)
            reference = real([-1.0d0/3.0d0, -1.0d0/3.0d0], c_double)
        case (QUADRILATERAL_FACETS)
            reference = 0.0_c_double
        case default
            ierr = 1
        end select
    end subroutine cell_reference_centroid

    !> PETSc's global number of a point, whether or not this rank owns it.
    pure function global_index_of(number) result(index)
        implicit none
        PetscInt, intent(in) :: number
        integer(int32) :: index

        if (number < 0) then
            index = int(-(number + 1), int32)
        else
            index = int(number, int32)
        end if
    end function global_index_of

end module app_mixed_hydraulic_balance
