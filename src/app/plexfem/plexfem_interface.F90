#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscviewerhdf5.h>

!> The PetscDS/TS solution path.
!>
!> Everything the finite element method needs beyond the weak form itself is
!> PETSc's: DMPlex owns the mesh, PetscFE the basis and its quadrature, PetscDS
!> the pointwise residual and Jacobian, and TS the time stepping and the
!> nonlinear solve. What this type does is read the project, build the DM, hand
!> PETSc the pointwise functions and the boundary conditions, and start TS.
!>
!> It exists beside type_ftcms rather than inside it because the two share no
!> machinery: the older path builds its own element loop, matrix, residual
!> vector, nonlinear iteration and time loop, none of which appear here.
module app_plexfem
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: petscdmplex
    use :: petscsnes
    use :: module_input, only:type_input
    use :: module_governing, only:type_thermal, type_hydraulic
    use :: app_plexfem_physics, only:PLEXFEM_CONSTANTS
    implicit none

    !> How many boundaries may carry a flux or a Robin condition at once.
    integer(int32), parameter :: MAX_NATURAL_BOUNDARIES = 32
    !> How many materials a mesh may be made of.
    integer(int32), parameter :: MAX_MATERIALS = 32
    private

    public :: type_plexfem

    type :: type_plexfem
        private
        !> The project, read once. It owns the DMPlex the solve DM clones, so
        !> it has to outlive the DM.
        type(type_input) :: input
        !> The energy equation's material response. The pointwise kernels reach
        !> it through app_plexfem_physics, so it has to live as long as the DM.
        type(type_thermal) :: thermal
        !> The water balance's material response, used when the project solves
        !> for pressure as well as temperature.
        type(type_hydraulic) :: hydraulic
        !> A clone of the mesh DM carrying this path's discretisation. Setting a
        !> field on the mesh's own DM would change the layout every other
        !> consumer sees.
        DM :: dm
        !> A second clone carrying the auxiliary field. PETSc reads the
        !> previous step's temperature from a local vector laid out on it.
        DM :: aux_dm
        Vec :: solution
        !> The last two converged states. The storage term differences the
        !> enthalpy across them: a first-order step uses one, a second-order
        !> step both, weighted by the BDF coefficients below.
        Vec :: old_solution
        Vec :: older_solution
        !> BDF coefficients of the step being solved, index 0 for the current
        !> state. They already carry the division by the step length, so the
        !> storage term multiplies rather than divides by them.
        real(real64) :: bdf_coeffs(0:2) = 0.0d0
        !> Order the step is being taken at. One until two states exist to
        !> difference, and one again after a rejection, because the history a
        !> second-order step needs describes step lengths that were not taken.
        integer(int32) :: bdf_order = 1
        !> Auxiliary values PETSc reads at the quadrature points: the material
        !> of each cell, and the constants of each natural boundary.
        Vec :: auxiliary
        !> Row scale of the residual, one entry per degree of freedom.
        !>
        !> Each equation is divided by the storage it would have for a unit
        !> change of its own primary variable, and then by that variable's
        !> tolerance. What comes back is dimensionless and comparable across the
        !> two physics: how many tolerances of temperature, or of pressure, the
        !> equation at that node is out by.
        !>
        !> Without it the energy rows, which are an energy rate, and the water
        !> rows, which are a mass rate, differ by some eighteen orders on this
        !> problem, so the norm SNES converges on and the merit its line search
        !> minimises are the energy rows to machine precision and the water
        !> equation is not solved at all so much as carried along.
        Vec :: row_scale
        !> What each unknown is measured in, one entry per degree of freedom.
        !>
        !> The nonlinear solver works on y = x / scale and every physical
        !> evaluation is handed x, so the two unknowns reach the solver in
        !> comparable units: a temperature in units of half a kelvin and a
        !> pressure in units of fifty kilopascals, which are the same amount of
        !> phase change. That is what lets one bound on the largest component of
        !> an update bound both fields, and it is what the trust region is set
        !> to. The conversion happens in exactly one place -- gather -- so no
        !> kernel, evaluator or boundary condition sees anything but physical
        !> values.
        Vec :: unknown_scale
        Vec :: unknown_scale_local
        !> Nodal volume \int psi_i of each material, the weight of the lumped
        !> storage term, and the material each one belongs to.
        Vec :: nodal_volume(MAX_MATERIALS)
        integer(int32) :: material_id(MAX_MATERIALS) = -1
        integer(int32) :: num_materials = 0
        !> One auxiliary vector per natural boundary, holding its condition and
        !> its constants. PETSc keys them by the boundary's label value.
        Vec :: boundary_aux(MAX_NATURAL_BOUNDARIES)
        integer(int32) :: num_boundary_aux = 0
        DMLabel :: face_sets
        SNES :: snes
        !> One HDF5 file for the whole run: the mesh is written once and each
        !> output time is a step inside it, so a run is one file rather than a
        !> directory of them.
        PetscViewer :: output_viewer
        logical :: output_ready = .false.
        !> A DM carrying one value per node, and a vector on it. Every quantity
        !> written goes out through this pair, so each arrives in the file as a
        !> dataset under its own name rather than as a component of the solution.
        DM :: output_dm
        Vec :: output_field
        !> The DS constants, kept because the step length is one of them and
        !> changes whenever a step is retried.
        PetscScalar :: constants(PLEXFEM_CONSTANTS)
        PetscDS :: ds
        logical :: dm_ready = .false.
        logical :: snes_ready = .false.
    contains
        procedure, public, pass(self) :: initialize => initialize_plexfem
        procedure, public, pass(self) :: run => run_plexfem
        procedure, public, pass(self) :: destroy => destroy_plexfem
    end type type_plexfem

    interface
        !> Read the project, build the DM and its discrete system, and register
        !> the pointwise functions and boundary conditions.
        !> self is a target because the pointwise kernels reach the governing
        !> model through a module pointer, which has to point at this
        !> instance's own component.
        module subroutine initialize_plexfem(self)
            implicit none
            class(type_plexfem), intent(inout), target :: self
        end subroutine initialize_plexfem

        !> Step from the start time to the end time.
        module subroutine run_plexfem(self)
            implicit none
            class(type_plexfem), intent(inout), target :: self
        end subroutine run_plexfem

        module subroutine destroy_plexfem(self)
            implicit none
            class(type_plexfem), intent(inout) :: self
        end subroutine destroy_plexfem
    end interface

end module app_plexfem
