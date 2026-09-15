#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscsnes.h>

!> An isothermal mixed hydraulic solver: RT0 flux, P0 pressure, DMPlex mesh.
!>
!> The flux is the primary unknown and lives in \(H(\mathrm{div})\), so its
!> normal trace is single valued across every facet and the conservation
!> equation holds cell by cell rather than only in a nodal sum. The pressure is
!> discontinuous and constant per cell, which is the space the divergence of the
!> flux lands in; that pairing is what makes the discrete conservation exact.
!>
!> The solver owns a clone of the caller's mesh, its own discrete system and its
!> own nonlinear solver. It does not own the mesh: the caller keeps it, and may
!> carry other fields on it.
!>
!> Everything crossing this interface is in physical units -- flux in m/s,
!> pressure in Pa, balances per unit out-of-plane thickness. The unknowns PETSc
!> stores are dimensionless; the scales that relate the two are set before
!> `initialize` and frozen for the run.
module app_mixed_hydraulic
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr
    use :: petscdmplex
    use :: petscsnes
    use :: app_mixed_hydraulic_physics, only:type_mixed_material, MIXED_CONSTANTS, MIXED_DIMENSION
    implicit none
    private

    public :: type_mixed_hydraulic
    public :: type_mixed_audit
    public :: MIXED_BC_PRESSURE, MIXED_BC_FLUX
    public :: MAX_MIXED_BOUNDARIES

    !> A boundary where the liquid pressure is prescribed. It reaches the weak
    !> form as the natural term \(\langle\chi_D,\mathbf v\cdot\mathbf n\rangle\),
    !> not as a constraint on any degree of freedom.
    integer(int32), parameter :: MIXED_BC_PRESSURE = 1
    !> A boundary where the normal flux is prescribed. It constrains the RT
    !> facet moments there, including when the prescribed flux is zero.
    integer(int32), parameter :: MIXED_BC_FLUX = 2

    !> How many boundaries one run may register.
    integer(int32), parameter :: MAX_MIXED_BOUNDARIES = 16

    !> One registered boundary of the run.
    type :: type_mixed_boundary
        !> The "Face Sets" value the boundary carries.
        integer(int32) :: boundary_id = 0
        !> MIXED_BC_PRESSURE or MIXED_BC_FLUX.
        integer(int32) :: kind = 0
        !> Prescribed liquid pressure, Pa.
        real(real64) :: pressure = 0.0_real64
        !> Prescribed flux vector, m/s; only its normal component is imposed.
        real(real64) :: flux(MIXED_DIMENSION) = 0.0_real64
    end type type_mixed_boundary

    !> What one balance audit measured, in physical units.
    type :: type_mixed_audit
        !> Largest \(|Q_{Kf}+Q_{Lf}|\) over internal facets, m2/s.
        real(real64) :: max_facet_mismatch = 0.0_real64
        !> Largest \(|\sum_f Q_{Kf}-\int_K\nabla\cdot\mathbf q|\), m2/s.
        real(real64) :: max_divergence_defect = 0.0_real64
        !> Largest backward-Euler cell defect, m2.
        real(real64) :: max_cell_defect = 0.0_real64
        !> Sum of the cell defects over owned cells, m2.
        real(real64) :: domain_defect = 0.0_real64
        !> Net outward flux through the domain boundary, m2/s.
        real(real64) :: boundary_rate = 0.0_real64
        !> Integrated source over owned cells, m2/s.
        real(real64) :: source_rate = 0.0_real64
        !> Largest reconstructed facet flux, m2/s: the audit's own scale.
        real(real64) :: throughput = 0.0_real64
        !> How many cells this rank owns.
        integer(int32) :: owned_cells = 0
    end type type_mixed_audit

    !> The mixed hydraulic solver.
    type :: type_mixed_hydraulic
        !> The solve DM, a clone of the caller's mesh carrying flux and pressure.
        DM :: dm
        !> The DM the cell-wise auxiliary data lives on.
        DM :: aux_dm
        !> The discrete system PETSc assembles from.
        PetscDS :: ds
        !> The dimensionless global solution.
        Vec :: solution
        !> The accepted state the storage term differences against.
        Vec :: old_solution
        !> The local form of `solution`, refreshed before every audit.
        Vec :: local_solution
        !> Cell-wise coefficients: inverse mobility, storage, source, history.
        Vec :: auxiliary
        !> One auxiliary vector per prescribed-pressure boundary.
        Vec :: boundary_aux(MAX_MIXED_BOUNDARIES)
        !> The nonlinear solver. It takes its residual, Jacobian and matrix from
        !> the DM, so the solver holds no matrix of its own.
        SNES :: snes
        !> The run constants PETSc hands the kernels.
        PetscScalar :: constants(MIXED_CONSTANTS)
        !> The "Face Sets" label, needed to key a boundary's auxiliary vector.
        DMLabel :: face_sets
        !> One entry per material of the mesh.
        type(type_mixed_material), allocatable :: materials(:)
        !> The registered boundaries.
        type(type_mixed_boundary) :: boundaries(MAX_MIXED_BOUNDARIES)
        integer(int32) :: num_boundaries = 0
        integer(int32) :: num_boundary_aux = 0
        !> The flux unit of the unknown, m/s.
        real(real64) :: flux_scale = 1.0_real64
        !> The pressure unit of the unknown, Pa.
        real(real64) :: pressure_scale = 1.0_real64
        !> The length the conservation row is scaled by, m.
        real(real64) :: length_scale = 1.0_real64
        !> Reference liquid density, kg/m3.
        real(real64) :: density = 1000.0_real64
        !> Gravity in intrinsic coordinates, m/s2.
        real(real64) :: gravity(MIXED_DIMENSION) = 0.0_real64
        !> The step the current residual differences over, s.
        real(real64) :: dt = 1.0_real64
        logical :: ready = .false.
    contains
        ! ---- Lifecycle ----
        procedure, public, pass :: initialize => initialize_mixed_hydraulic
        procedure, public, pass :: destroy => destroy_mixed_hydraulic
        ! ---- Mutator ----
        procedure, public, pass :: set_materials
        procedure, public, pass :: set_gravity
        procedure, public, pass :: set_scales
        procedure, public, pass :: add_boundary
        procedure, public, pass :: set_initial_state
        ! ---- Algorithm ----
        procedure, public, pass :: solve_step
        procedure, public, pass :: compute_balance
        procedure, public, pass :: compute_audit
        ! ---- Getter ----
        procedure, public, pass :: get_dm
        procedure, public, pass :: get_local_solution
        procedure, public, pass :: get_owned_cells
    end type type_mixed_hydraulic

    interface
        !> Build the discrete system, the boundaries and the solver on `dm`.
        !>
        !> \[
        !>   R_q(\mathbf v)=(\mathbf\Lambda^{-1}\mathbf q_l-\rho_l\mathbf g,\mathbf v)
        !>   -(\chi,\nabla\cdot\mathbf v)+\langle\chi_D,\mathbf v\cdot\mathbf n\rangle_{\Gamma_D},
        !>   \qquad R_p(w)=(D_tW+\nabla\cdot\mathbf q_l-s,w).
        !> \]
        !>
        !> Assumes an interpolated simplicial DMPlex whose dimension equals its
        !> coordinate dimension, one "Cell Sets" material per cell, and a
        !> symmetric positive definite conductivity for each material. The
        !> materials, gravity and scales must already be set.
        !>
        !> Guarantees an RT0 flux space with one degree of freedom per facet and
        !> a P0 pressure space with one per cell; it does not guarantee that the
        !> resulting saddle-point system is nonsingular, which depends on the
        !> boundaries the caller registered. Setup costs \(O(N)\) in cells for
        !> the coefficients and whatever PETSc's section construction costs.
        !>
        !> Fails by `error stop` with the refusing condition named: a mesh of the
        !> wrong dimension, a cell without a material, or a conductivity that is
        !> not positive definite. A zero conductivity is refused, never floored.
        module subroutine initialize_mixed_hydraulic(self, dm, ierr)
            implicit none
            class(type_mixed_hydraulic), intent(inout), target :: self
            !> An interpolated DMPlex; the solver clones it and does not own it.
            DM, intent(in) :: dm
            PetscErrorCode, intent(inout) :: ierr
        end subroutine initialize_mixed_hydraulic

        !> Release every PETSc object the solver made.
        module subroutine destroy_mixed_hydraulic(self, ierr)
            implicit none
            class(type_mixed_hydraulic), intent(inout) :: self
            PetscErrorCode, intent(inout) :: ierr
        end subroutine destroy_mixed_hydraulic

        !> Solve one backward-Euler step of length `dt` seconds.
        !>
        !> \[
        !>   \int_K\frac{W^{n+1}-W^{n}}{\Delta t}+\int_K\nabla\cdot\mathbf q_l^{n+1}
        !>   =\int_K s.
        !> \]
        !>
        !> Assumes the accepted history is the state of the previous successful
        !> call. Guarantees that the history is advanced only when SNES converged;
        !> `accepted = .false.` leaves it exactly as it was, so a rejected step
        !> can be retried at a different length. Cost is one nonlinear solve.
        !>
        !> Fails by returning `accepted = .false.` with the PETSc error in `ierr`
        !> when the solve did not converge.
        module subroutine solve_step(self, time, dt, accepted, ierr)
            implicit none
            class(type_mixed_hydraulic), intent(inout) :: self
            !> The time the step ends at, s.
            real(real64), intent(in) :: time
            !> The step length, s; positive.
            real(real64), intent(in) :: dt
            !> True when the step converged and the history advanced.
            logical, intent(inout) :: accepted
            PetscErrorCode, intent(inout) :: ierr
        end subroutine solve_step

        !> The backward-Euler water balance of each owned cell.
        !>
        !> \[
        !>   d_K=\int_K(W^{n+1}-W^{n})+\Delta t\sum_{f\subset\partial K}Q_{Kf}
        !>   -\Delta t\int_K s.
        !> \]
        !>
        !> Assumes the solution vector holds the state the balance is claimed
        !> for. The facet fluxes are reconstructed from the degrees of freedom
        !> rather than read out of the assembled residual, so this is an
        !> independent measurement. Cost is \(O(N)\) in owned cells.
        !>
        !> Fails by returning the first PETSc error; all outputs are initialised
        !> on every path. `cell_defect` beyond the owned-cell count is zeroed.
        module subroutine compute_balance(self, cell_defect, boundary_rate, domain_defect, ierr)
            implicit none
            class(type_mixed_hydraulic), intent(inout) :: self
            !> One defect per owned cell, m2 of water-equivalent volume.
            real(real64), intent(inout) :: cell_defect(:)
            !> Net outward flux through the domain boundary, m2/s.
            real(real64), intent(inout) :: boundary_rate
            !> Sum of the cell defects over all ranks, m2.
            real(real64), intent(inout) :: domain_defect
            PetscErrorCode, intent(inout) :: ierr
        end subroutine compute_balance

        !> Everything `compute_balance` measures, plus the geometry checks.
        !>
        !> Adds the facet continuity mismatch and the divergence-theorem defect,
        !> each reconstructed from both sides independently, and reduces the
        !> maxima and sums over all ranks. Cost is \(O(N)\) in local cells.
        module subroutine compute_audit(self, report, ierr)
            implicit none
            class(type_mixed_hydraulic), intent(inout) :: self
            !> The measured quantities, reduced over the communicator.
            type(type_mixed_audit), intent(inout) :: report
            PetscErrorCode, intent(inout) :: ierr
        end subroutine compute_audit
    end interface

contains

    !> The materials the mesh is made of, keyed by their "Cell Sets" value.
    subroutine set_materials(self, materials)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        type(type_mixed_material), intent(in) :: materials(:)

        if (allocated(self%materials)) deallocate (self%materials)
        allocate (self%materials(size(materials)))
        self%materials = materials
    end subroutine set_materials

    !> Gravity in the intrinsic plane, m/s2; downward is negative in z.
    subroutine set_gravity(self, gravity, density)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        real(real64), intent(in) :: gravity(MIXED_DIMENSION)
        !> Reference liquid density, kg/m3.
        real(real64), intent(in) :: density

        self%gravity = gravity
        self%density = density
    end subroutine set_gravity

    !> The units the dimensionless unknowns are measured in.
    !>
    !> Set them from the throughput and the pressure drop the problem actually
    !> carries. An absolute tolerance of order one in m/s would accept a zero
    !> answer wherever the conductivity is small.
    subroutine set_scales(self, flux_scale, pressure_scale, length_scale)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        !> The flux unit, m/s; positive.
        real(real64), intent(in) :: flux_scale
        !> The pressure unit, Pa; positive.
        real(real64), intent(in) :: pressure_scale
        !> The length the conservation row is scaled by, m; positive.
        real(real64), intent(in) :: length_scale

        if (flux_scale <= 0.0_real64 .or. pressure_scale <= 0.0_real64 .or. length_scale <= 0.0_real64) then
            error stop "The mixed hydraulic scales have to be positive"
        end if
        self%flux_scale = flux_scale
        self%pressure_scale = pressure_scale
        self%length_scale = length_scale
    end subroutine set_scales

    !> Register one boundary of the run, before `initialize`.
    subroutine add_boundary(self, boundary_id, kind, pressure, flux)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        !> The "Face Sets" value the boundary carries.
        integer(int32), intent(in) :: boundary_id
        !> MIXED_BC_PRESSURE or MIXED_BC_FLUX.
        integer(int32), intent(in) :: kind
        !> Prescribed pressure, Pa; read for MIXED_BC_PRESSURE.
        real(real64), intent(in), optional :: pressure
        !> Prescribed flux vector, m/s; read for MIXED_BC_FLUX.
        real(real64), intent(in), optional :: flux(MIXED_DIMENSION)

        if (self%num_boundaries >= MAX_MIXED_BOUNDARIES) then
            error stop "More boundaries than the mixed hydraulic solver can hold"
        end if
        self%num_boundaries = self%num_boundaries + 1
        self%boundaries(self%num_boundaries)%boundary_id = boundary_id
        self%boundaries(self%num_boundaries)%kind = kind
        self%boundaries(self%num_boundaries)%pressure = 0.0_real64
        self%boundaries(self%num_boundaries)%flux = 0.0_real64
        if (present(pressure)) self%boundaries(self%num_boundaries)%pressure = pressure
        if (present(flux)) self%boundaries(self%num_boundaries)%flux = flux
    end subroutine add_boundary

    !> Start every cell from one pressure and a flux of zero.
    subroutine set_initial_state(self, pressure, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        !> The initial liquid pressure, Pa.
        real(real64), intent(in) :: pressure
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: cell, c_start, c_end, offset, num_dofs
        Vec :: local

        call VecZeroEntries(self%solution, ierr)
        if (ierr /= 0) return

        call DMGetLocalVector(self%dm, local, ierr)
        if (ierr == 0) call VecZeroEntries(local, ierr)
        if (ierr == 0) call DMGetLocalSection(self%dm, section, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(self%dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return

        call VecGetArray(local, raw, ierr)
        if (ierr /= 0) return
        do cell = c_start, c_end - 1
            call PetscSectionGetFieldDof(section, cell, 1_PETSC_INT_KIND, num_dofs, ierr)
            if (ierr /= 0 .or. num_dofs <= 0) cycle
            call PetscSectionGetFieldOffset(section, cell, 1_PETSC_INT_KIND, offset, ierr)
            if (ierr == 0) raw(offset + 1) = real(pressure/self%pressure_scale, PETSC_REAL_KIND)
        end do
        call VecRestoreArray(local, raw, ierr)

        if (ierr == 0) call DMLocalToGlobalBegin(self%dm, local, INSERT_VALUES, self%solution, ierr)
        if (ierr == 0) call DMLocalToGlobalEnd(self%dm, local, INSERT_VALUES, self%solution, ierr)
        call DMRestoreLocalVector(self%dm, local, ierr)
        if (ierr == 0) call VecCopy(self%solution, self%old_solution, ierr)
    end subroutine set_initial_state

    !> The DM the flux and pressure live on.
    function get_dm(self) result(dm)
        implicit none
        class(type_mixed_hydraulic), intent(in) :: self
        DM :: dm

        dm = self%dm
    end function get_dm

    !> The local form of the current solution, refreshed from the global one.
    !>
    !> The prescribed normal fluxes are inserted as well. They are constrained
    !> out of the global vector, so a local vector that only came through the
    !> scatter would carry a zero flux on every prescribed boundary and the
    !> balance audit would read an inflow of nothing.
    subroutine get_local_solution(self, local, ierr)
        implicit none
        class(type_mixed_hydraulic), intent(inout) :: self
        !> The refreshed local vector; owned by the solver.
        Vec, intent(inout) :: local
        PetscErrorCode, intent(inout) :: ierr

        call DMGlobalToLocalBegin(self%dm, self%solution, INSERT_VALUES, self%local_solution, ierr)
        if (ierr == 0) call DMGlobalToLocalEnd(self%dm, self%solution, INSERT_VALUES, self%local_solution, ierr)
        if (ierr == 0) call DMPlexInsertBoundaryValues(self%dm, PETSC_TRUE, self%local_solution, &
                                                       0.0_PETSC_REAL_KIND, PETSC_NULL_VEC, &
                                                       PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
        local = self%local_solution
    end subroutine get_local_solution

    !> How many cells this rank owns.
    function get_owned_cells(self) result(owned)
        implicit none
        class(type_mixed_hydraulic), intent(in) :: self
        integer(int32) :: owned

        PetscErrorCode :: ierr
        PetscInt :: cell, c_start, c_end, global_start, global_end

        owned = 0
        call DMPlexGetHeightStratum(self%dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return
        do cell = c_start, c_end - 1
            call DMPlexGetPointGlobal(self%dm, cell, global_start, global_end, ierr)
            if (ierr == 0 .and. global_start >= 0) owned = owned + 1
        end do
    end function get_owned_cells

end module app_mixed_hydraulic
