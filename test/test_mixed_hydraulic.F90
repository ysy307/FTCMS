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
    public :: manufactured_temperature, manufactured_flux, manufactured_pressure
    public :: manufactured_temperature_at, manufactured_flux_at
    public :: MANUFACTURED_T, MANUFACTURED_FLUX_SLOPE, MANUFACTURED_P0

    !> The field \(\mathbf q=\mathbf b+s\mathbf x\), whose divergence is
    !> \(d\,s\) in \(d\) dimensions and which lies exactly in RT0.
    type, bind(C) :: type_linear_field
        real(c_double) :: base(MIXED_DIMENSION)
        real(c_double) :: slope
        real(c_double) :: pressure
    end type type_linear_field

    type(type_linear_field), save :: projected_field

    ! ---- The manufactured state of the rectangle case ----
    !
    ! Each piece is chosen to be exactly representable in the space it is
    ! projected into on an axis-aligned rectangle, and to be outside the
    ! corresponding space on a triangle of the same mesh, so that a value read
    ! back from the discretisation is either exact or wrong rather than
    ! approximate.

    !> \(T=T_1+T_2x+T_3z+T_4xz\), C. The bilinear term is the one a P1 basis
    !> cannot carry and a Q1 basis reproduces.
    real(real64), parameter :: MANUFACTURED_T(4) = [1.0d0, 2.0d0, 3.0d0, 5.0d0]
    !> \(\mathbf q=(0,\,s z)\), m/s: no x component anywhere, and
    !> \(\nabla\cdot\mathbf q=s\), which is nonzero and constant.
    real(real64), parameter :: MANUFACTURED_FLUX_SLOPE = 1.0d0
    !> The cell-constant pressure that goes with it, Pa.
    real(real64), parameter :: MANUFACTURED_P0 = 7.0d0

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

    !> The manufactured temperature at a point, C.
    pure function manufactured_temperature_at(x, z) result(temperature)
        implicit none
        real(real64), intent(in) :: x, z
        real(real64) :: temperature

        temperature = MANUFACTURED_T(1) + MANUFACTURED_T(2)*x + MANUFACTURED_T(3)*z &
                      + MANUFACTURED_T(4)*x*z
    end function manufactured_temperature_at

    !> The manufactured flux at a point, m/s.
    pure function manufactured_flux_at(x, z) result(flux)
        implicit none
        real(real64), intent(in) :: x, z
        real(real64) :: flux(MIXED_DIMENSION)

        flux = [0.0d0, MANUFACTURED_FLUX_SLOPE*z]
    end function manufactured_flux_at

    !> The manufactured temperature, as PETSc's projection calls it.
    function manufactured_temperature(dim, time, x, num_components, u, ctx) &
        bind(C, name="ftcms_test_mixed_manufactured_temperature") result(ierr)
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        u(1) = real(manufactured_temperature_at(real(x(1), real64), real(x(2), real64)), PETSC_REAL_KIND)
        ierr = 0
    end function manufactured_temperature

    !> The manufactured flux, as PETSc's projection calls it.
    function manufactured_flux(dim, time, x, num_components, u, ctx) &
        bind(C, name="ftcms_test_mixed_manufactured_flux") result(ierr)
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        real(real64) :: flux(MIXED_DIMENSION)
        integer(int32) :: i

        flux = manufactured_flux_at(real(x(1), real64), real(x(2), real64))
        do i = 1, MIXED_DIMENSION
            u(i) = real(flux(i), PETSC_REAL_KIND)
        end do
        ierr = 0
    end function manufactured_flux

    !> The manufactured pressure, as PETSc's projection calls it.
    function manufactured_pressure(dim, time, x, num_components, u, ctx) &
        bind(C, name="ftcms_test_mixed_manufactured_pressure") result(ierr)
        implicit none
        PetscInt, value :: dim
        PetscReal, value :: time
        PetscReal, intent(in) :: x(*)
        PetscInt, value :: num_components
        PetscScalar, intent(inout) :: u(*)
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        u(1) = real(MANUFACTURED_P0, PETSC_REAL_KIND)
        ierr = 0
    end function manufactured_pressure

end module test_mixed_hydraulic_fields

!> Synthetic cases for the RT0/P0 mixed hydraulic solver.
!>
!> Each case is selected with `-mixed_case <name>` and needs no project
!> directory: the meshes are built here, so the geometry, the materials and the
!> boundaries of every check are visible in one file.
!>
!> The steady reference solutions are exact. Refinement cases compare numerical
!> profiles and do not treat the finest mesh as exact. For a mesh whose
!> material interface lies on cell boundaries, the lowest-order mixed method
!> reproduces a piecewise-linear pressure's cell averages and a constant flux
!> exactly, because the divergence of an RT0 field is exactly the P0 space the
!> pressure is tested against. That makes a relative flux error of `1e-8` a
!> statement about the solver rather than about the mesh resolution.
program test_mixed_hydraulic
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64, output_unit
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_funptr, c_funloc, c_loc, c_f_pointer
    use, intrinsic :: ieee_arithmetic, only:ieee_value, ieee_quiet_nan
    use :: petscdmplex
    use :: petscsnes
    use :: core_interop_petscds, only:petscds_project_function
    use :: core_interop_petscfe, only:petscfe_field_size, petscfe_evaluate_cell, &
                                      petscfe_reference_coordinates
    use :: iapws, only:type_iapws97, type_iapws06
    use :: module_core, only:SWCC_MODELS, HCF_MODES, HCF_VISCOSITY_TYPES, GCC_TYPES, COMP_TYPES, &
                              type_state, type_config_hcf
    use :: models_hcf, only:holder_hcfs
    use :: app_mixed_hydraulic_closure, only:install_storage_laws, clear_storage_laws, eval_storage, &
                                             install_conductivity_laws, clear_conductivity_laws, &
                                             eval_relative_conductivity
    use :: app_mixed_hydraulic_frozen, only:type_mixed_frozen_material, type_mixed_frozen_response, &
                                            type_mixed_frozen_tangents, install_frozen_laws, &
                                            clear_frozen_laws, publish_frozen_state, &
                                            eval_frozen_response, eval_frozen_tangents, &
                                            frozen_pressure_column, configure_frozen_tangent_cache, &
                                            frozen_work_counters
    use :: app_mixed_hydraulic, only:type_mixed_hydraulic, type_mixed_audit, &
                                     MIXED_BC_PRESSURE, MIXED_BC_FLUX, MIXED_ERR_INCOMPATIBLE, &
                                     MIXED_ERR_INVALID_STEP
    use :: app_mixed_hydraulic_physics, only:type_mixed_material, MIXED_DIMENSION, &
                                             intrinsic_from_physical, physical_from_intrinsic, &
                                             MIXED_STORAGE_RETENTION, &
                                             MIXED_CONDUCTIVITY_FIXED, MIXED_CONDUCTIVITY_HCF
    use :: app_mixed_hydraulic_balance, only:MAX_CELL_FACETS, reconstruct_cell_facet_flux, &
                                             reconstruct_cell_divergence, reconstruct_cell_pressure, &
                                             audit_facet_continuity, cell_is_owned, cell_volume
    use :: app_mixed_hydraulic_column, only:type_mixed_column, type_mixed_column_report, &
                                            type_mixed_column_step_control, type_mixed_column_progress, &
                                            COLUMN_ERR_INVALID_STEP
    use :: app_mixed_hydraulic_column_pointwise, only:COLUMN_FIELD_TEMPERATURE, COLUMN_FIELD_FLUX, &
                                                      COLUMN_FIELD_PRESSURE
    use :: test_mixed_hydraulic_fields, only:projected_field, analytic_flux, analytic_pressure, &
                                             manufactured_temperature, manufactured_flux, &
                                             manufactured_pressure, manufactured_temperature_at, &
                                             manufactured_flux_at, MANUFACTURED_FLUX_SLOPE, &
                                             MANUFACTURED_P0
    implicit none

    !> Standard gravity, m/s2. The conductivities of the reference cases were
    !> measured with it, so it is what turns them into pressure mobilities.
    real(real64), parameter :: STANDARD_GRAVITY = 9.80665d0
    !> Reference liquid density, kg/m3.
    real(real64), parameter :: REFERENCE_DENSITY = 1000.0d0
    !> \(\rho_{\mathrm{std}}g\): a pore pressure in Pa over a suction head in m.
    real(real64), parameter :: SPECIFIC_WEIGHT = REFERENCE_DENSITY*STANDARD_GRAVITY

    ! The retention curve of the nonlinear storage cases: van Genuchten with the
    ! parameters the repository's own notes quote for the Mizoguchi column. The
    ! formula is written out again in `reference_storage` below so that what the
    ! solver is compared against does not come from the code under test.
    real(real64), parameter :: SOIL_THETA_S = 0.535d0
    real(real64), parameter :: SOIL_THETA_R = 0.05d0
    real(real64), parameter :: SOIL_ALPHA = 1.11d0
    real(real64), parameter :: SOIL_N = 1.48d0
    real(real64), parameter :: SOIL_M = 0.2d0
    !> Specific storage of the saturated branch, 1/m.
    real(real64), parameter :: SOIL_SPECIFIC_STORAGE = 1.0d-4
    !> Its capacity in 1/Pa, which is the whole tangent above the air entry.
    real(real64), parameter :: SOIL_COMPRESSIVE_CAPACITY = SOIL_SPECIFIC_STORAGE/SPECIFIC_WEIGHT
    !> Saturated hydraulic conductivity of those cases, m/s.
    real(real64), parameter :: SOIL_CONDUCTIVITY = 1.0d-6
    !> Pore-connectivity exponent of its conductivity curve.
    real(real64), parameter :: SOIL_L = 0.5d0

    ! The frozen-closure case runs the same van Genuchten soil at the rest of the
    ! Mizoguchi column's own material data, so that what the contract is measured
    ! on is the material the validation problem uses.
    !> Region identifier the frozen material is written on.
    integer(int32), parameter :: FROZEN_MATERIAL_ID = 1
    !> Phases the material declares: solid, water, ice and vapour.
    integer(int32), parameter :: FROZEN_PHASES = 4
    !> Saturated conductivity of the column, m/s.
    real(real64), parameter :: FROZEN_CONDUCTIVITY = 3.2d-6
    !> Specific storage of the column, 1/m.
    real(real64), parameter :: FROZEN_SPECIFIC_STORAGE = 2.3d-6
    !> Impedance factor, the physical benchmark parameter of the column.
    real(real64), parameter :: FROZEN_IMPEDANCE_FACTOR = 7.0d0
    !> The column's initial temperature, C.
    real(real64), parameter :: FROZEN_WARM = 6.7d0
    !> A freezing-dominated temperature on a smooth branch, C.
    real(real64), parameter :: FROZEN_MILD_COLD = -0.5d0
    !> The cold-end temperature of the column, C.
    real(real64), parameter :: FROZEN_DEEP_COLD = -6.0d0
    !> Absolute zero in the closure's own temperature unit, C. Written here
    !> rather than imported, so the domain test states the boundary it probes.
    real(real64), parameter :: FROZEN_ABSOLUTE_ZERO = -273.15d0
    !> The column's initial pore pressure, Pa.
    real(real64), parameter :: FROZEN_FIELD_PRESSURE = -44683.6901d0
    !> A nearly saturated pore pressure, where the pore volume limits the ice, Pa.
    real(real64), parameter :: FROZEN_WET_PRESSURE = -1.0d3
    !> A positive pore pressure, where the compressive branch is the capacity, Pa.
    real(real64), parameter :: FROZEN_POSITIVE_PRESSURE = 2.0d3

    !> Height of the coupled column, m. The Mizoguchi cylinder is 20 cm tall.
    real(real64), parameter :: COLUMN_HEIGHT = 0.20d0
    !> Width of the modelled strip, m. The physics is one-dimensional; the
    !> width exists so the mesh is two-dimensional.
    real(real64), parameter :: COLUMN_WIDTH = 0.02d0
    !> Cells across and along the column.
    integer(int32), parameter :: COLUMN_NX = 2
    integer(int32), parameter :: COLUMN_NZ = 16
    !> "Face Sets" value of the closed water boundary: the whole perimeter.
    integer(int32), parameter :: COLUMN_CLOSED_ID = 1
    !> "Thermal Sets" values of the cooled top and the held bottom.
    integer(int32), parameter :: COLUMN_TOP_ID = 10
    integer(int32), parameter :: COLUMN_BOTTOM_ID = 20
    !> The temperature the column starts at and its ends are held at, C.
    real(real64), parameter :: COLUMN_WARM = FROZEN_WARM
    !> The temperature the top is cooled to, C.
    real(real64), parameter :: COLUMN_COLD = -6.0d0
    !> Steps of the warm control, s, and how many of them.
    real(real64), parameter :: COLUMN_CONTROL_STEP = 600.0d0
    integer(int32), parameter :: COLUMN_CONTROL_STEPS = 3
    !> The step the cooling is first attempted at, s, and the span it covers.
    real(real64), parameter :: COLUMN_COOLING_STEP = 900.0d0
    real(real64), parameter :: COLUMN_COOLING_SPAN = 7200.0d0

    ! ---- The Mizoguchi freezing experiment, as project/Mizo-xz-Convex states it ----
    !
    ! Every number below was read out of that project's Input/Conditions.json
    ! and Input/Mizoguchi_unsaturated_v.msh, and the material is the
    ! `set_frozen_soil` above, which is Input/Basic.json's material 1. The
    ! boundary the project calls 5 is the mesh's z = 0.2 edge and the one it
    ! calls 3 is its z = 0 edge, so the coolant is at the top and the leaking
    ! base at the bottom. Hansson et al. (2004) give 28 W/(m2 K) against a
    ! -6 C coolant for the top; the warm base is that paper's heat-leakage
    ! variant and not its adiabatic base case.

    !> The column the experiment was run in: 20 cm of soil.
    real(real64), parameter :: MIZO_HEIGHT = 0.20d0
    !> Width of the modelled strip, m. One material, adiabatic and impermeable
    !> sides and a uniform initial state make the problem one-dimensional, so
    !> this is the thickness the per-unit-thickness reports are taken over and
    !> not a second direction the answer depends on.
    real(real64), parameter :: MIZO_WIDTH = 0.02d0
    !> Cells across and along. The 64 rows are 3.125 mm, which is the depth-wise
    !> edge length of the project's own mesh (its 0.2 m sides carry 66 edges).
    !> One resolution is run, and no claim of mesh convergence is made from it.
    integer(int32), parameter :: MIZO_NX = 2
    integer(int32), parameter :: MIZO_NZ = 64
    !> The uniform initial state: 6.7 C and -44683.6901 Pa of pore water.
    real(real64), parameter :: MIZO_INITIAL_TEMPERATURE = 6.7d0
    real(real64), parameter :: MIZO_INITIAL_PRESSURE = -44683.6901d0
    !> Top, `-mizo_top_bc convective`: a coolant at -6 C behind a surface
    !> resistance of 28 W/(m2 K), Hansson et al. (2004)'s trial-and-error value.
    !> The coefficient is a calibration parameter; set it with
    !> `-mizo_top_coefficient`, and the impedance factor with
    !> `-mizo_impedance_factor`.
    real(real64), parameter :: MIZO_TOP_COEFFICIENT = 28.0d0
    real(real64), parameter :: MIZO_COOLANT = -6.0d0
    !> Top, `-mizo_top_bc plate` (the default): the experiment cooled the soil
    !> through a plate in contact with it, so the soil surface is at the
    !> plate's temperature. The values at 1, 3, 6, 12, 24 and 45 h are the
    !> measured surface temperatures of docs/data/soil_temperature_profiles.csv.
    !> The 0 h value was not measured: it is the intercept of
    !> \(T_\infty+A\,e^{-t/\tau}\) fitted to those six (\(T_\infty=-6.71\) C,
    !> \(A=3.99\) K, \(\tau=2.17\) h, residuals within 0.1 K).
    integer(int32), parameter :: MIZO_PLATE_POINTS = 7
    real(real64), parameter :: MIZO_PLATE_HOURS(MIZO_PLATE_POINTS) = &
                               [0.0d0, 1.0d0, 3.0d0, 6.0d0, 12.0d0, 24.0d0, 45.0d0]
    real(real64), parameter :: MIZO_PLATE_TEMPERATURES(MIZO_PLATE_POINTS) = &
                               [-2.7d0, -4.2d0, -5.7d0, -6.5d0, -6.6d0, -6.7d0, -6.8d0]
    !> Top, `-mizo_top_bc contact`: the same measured plate temperatures as the
    !> medium of a heat-transfer boundary, the coefficient being the thermal
    !> contact conductance between plate and soil, W/(m2 K). Not measured; set
    !> with `-mizo_contact_coefficient`.
    real(real64), parameter :: MIZO_CONTACT_COEFFICIENT = 200.0d0
    !> Base, `-mizo_base_bc measured`: the measured 20 cm temperatures of the
    !> same file at the plate times, with the initial 6.7 C at 0 h. A
    !> diagnostic: with both ends on measurements, what is left to disagree
    !> inside the column is its own thermal physics, not a boundary guess.
    real(real64), parameter :: MIZO_BASE_TEMPERATURES(MIZO_PLATE_POINTS) = &
                               [6.7d0, 7.1d0, 5.9d0, 4.4d0, 2.4d0, 1.4d0, 1.1d0]
    !> Bottom: the room the base leaks to, 6.7 C behind 3 W/(m2 K).
    real(real64), parameter :: MIZO_BASE_COEFFICIENT = 3.0d0
    real(real64), parameter :: MIZO_AMBIENT = 6.7d0
    !> The hours the experiment reports a profile at.
    integer(int32), parameter :: MIZO_HORIZONS = 3
    real(real64), parameter :: MIZO_HOURS(MIZO_HORIZONS) = [12.0d0, 24.0d0, 50.0d0]
    !> How many horizons one run may be asked to stop and report at.
    integer(int32), parameter :: MAX_MIZO_HORIZONS = 12
    !> "Thermal Sets" values of the cooled top and the leaking base.
    integer(int32), parameter :: MIZO_TOP_ID = 5
    integer(int32), parameter :: MIZO_BASE_ID = 3
    !> The step the run starts at and the longest one it may grow to, s. The
    !> first is the project's `initial_step`; the ceiling is well under its
    !> `max_step` of 30000 s because a step that long cannot resolve a front
    !> crossing a 3.125 mm cell.
    real(real64), parameter :: MIZO_INITIAL_STEP = 300.0d0
    real(real64), parameter :: MIZO_MAXIMUM_STEP = 1800.0d0
    real(real64), parameter :: MIZO_MINIMUM_STEP = 1.0d0

    ! ---- The cell profile export ----
    !
    ! The export writes one row per cell of the whole domain whatever the rank
    ! count is, which is a claim about the gather and not about the solve. The
    ! grid is the experiment's, small enough that eight ranks own real cells
    ! each and that a test can afford to solve on it.
    !
    !> How many numbers a profile row carries beside its time and cell number.
    integer(int32), parameter :: PROFILE_COLUMNS = 13
    !> The columns of a profile row that the comparison reports by name.
    integer(int32), parameter :: PROFILE_VOLUME = 3
    integer(int32), parameter :: PROFILE_STORED = 10
    !> Which columns are judged on a common scale.
    !>
    !> Each column is otherwise measured against its own range. The two flux
    !> components are the exception: they are one vector quantity in one unit,
    !> and on a column with no lateral driving force the \(x\) component is
    !> rounding error all the way, so its own range would make a difference of
    !> \(10^{-19}\,\mathrm{m/s}\) look like a relative difference of one.
    integer(int32), parameter :: PROFILE_SCALE_GROUP(PROFILE_COLUMNS) = &
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12]
    !> The default grid of the export check, and how long it is stepped for.
    integer(int32), parameter :: EXPORT_NX = 4
    integer(int32), parameter :: EXPORT_NZ = 16
    real(real64), parameter :: EXPORT_SPAN = 1800.0d0
    !> How closely two rank counts have to agree, as a fraction of each
    !> column's own range. The initial state is the same closure evaluated at
    !> the same uniform state on both, so nothing but the geometry differs; the
    !> solved state additionally carries the distributed factorisation's own
    !> pivot order, which is not the serial one.
    !> Measured across 1, 2, 8 and 16 ranks with the cooling-plate top: the
    !> initial state agrees exactly and the solved state to 2.9e-10. The plate
    !> freezes the top cells within the span, and the latent-heat-dominated
    !> Jacobian amplifies the pivot-order rounding into ~1e-9 K of temperature
    !> and, through the Clapeyron slope, into the liquid pressure; the solved
    !> p_w agrees to 1e-14. (The ice-free convective top agreed to 2.2e-13.)
    !> These leave a factor of a few hundred, far below the O(1) difference a
    !> wrong gather makes.
    real(real64), parameter :: EXPORT_INITIAL_TOLERANCE = 1.0d-12
    real(real64), parameter :: EXPORT_SOLVED_TOLERANCE = 1.0d-7

    ! ---- The rectangle geometry case ----
    !
    ! Four cells, each twenty times wider than it is tall. The elongation is
    ! the point: on a square the reference centroid, the two tensor components
    ! of the flux and the two directions of the Jacobian are all
    ! interchangeable, so a mapping error cancels. Here it does not.
    !
    ! Two cells in each direction rather than one: the manufactured flux is
    ! (0, z), whose trace on a vertical facet is zero, so a single layer would
    ! compare the two sides of its one interior facet at zero and see nothing.
    ! The second layer puts an interior horizontal facet in the mesh, which is
    ! where this field actually transports, and the continuity check then runs
    ! on a nonzero moment.

    integer(int32), parameter :: RECT_NX = 2
    integer(int32), parameter :: RECT_NZ = 2
    real(real64), parameter :: RECT_WIDTH = 1.0d0
    real(real64), parameter :: RECT_HEIGHT = 0.05d0
    !> How many degrees of freedom each field of a rectangle carries.
    integer(int32), parameter :: RECT_TEMPERATURE_DOFS = 4
    integer(int32), parameter :: RECT_FLUX_DOFS = 4
    integer(int32), parameter :: RECT_PRESSURE_DOFS = 1
    !> Interior points one cell's projected flux is sampled at: the four
    !> quarter-points, which lie off both mid-lines of the cell.
    integer(int32), parameter :: RECTANGLE_SAMPLES = 4

    ! The second soil of the material-interface case: a Brooks-Corey curve, so
    ! that one mesh carries two different retention laws and two different
    ! conductivity laws at once. Its saturation head is negative, unlike the van
    ! Genuchten soil's, which is why its saturation pressure is not zero.
    real(real64), parameter :: CLAY_THETA_S = 0.40d0
    real(real64), parameter :: CLAY_THETA_R = 0.03d0
    !> Air-entry head, m; negative by the Brooks-Corey convention.
    real(real64), parameter :: CLAY_ALPHA = -0.5d0
    real(real64), parameter :: CLAY_N = 1.5d0
    real(real64), parameter :: CLAY_L = 0.5d0
    real(real64), parameter :: CLAY_SPECIFIC_STORAGE = 1.0d-4
    real(real64), parameter :: CLAY_COMPRESSIVE_CAPACITY = CLAY_SPECIFIC_STORAGE/SPECIFIC_WEIGHT
    !> The pressure its pore volume fills at, Pa.
    real(real64), parameter :: CLAY_SATURATION_PRESSURE = CLAY_ALPHA*SPECIFIC_WEIGHT
    !> Saturated hydraulic conductivity of that soil, m/s.
    real(real64), parameter :: CLAY_CONDUCTIVITY = 4.0d-7

    ! The closed fixture the nonlinear storage and the rollback cases share.
    integer(int32), parameter :: CLOSED_NX = 4, CLOSED_NZ = 2
    real(real64), parameter :: CLOSED_WIDTH = 1.0d0, CLOSED_HEIGHT = 0.5d0
    !> Area of one cell of that fixture, m2 per unit thickness.
    real(real64), parameter :: CLOSED_CELL_AREA = &
        0.5d0*(CLOSED_WIDTH/real(CLOSED_NX, real64))*(CLOSED_HEIGHT/real(CLOSED_NZ, real64))

    ! The wetting fixture the nonlinear-conductivity cases share. It is short
    ! enough that the front crosses the material interface within the run, which
    ! is what lets a homogeneous control tell the two soils apart.
    integer(int32), parameter :: WET_NX = 8, WET_NZ = 2
    real(real64), parameter :: WET_WIDTH = 0.1d0, WET_HEIGHT = 0.04d0
    !> Where the interface case changes material, m; on a cell boundary.
    real(real64), parameter :: WET_SPLIT = 0.05d0
    !> Area of one of its cells, m2 per unit thickness.
    real(real64), parameter :: WET_CELL_AREA = &
        0.5d0*(WET_WIDTH/real(WET_NX, real64))*(WET_HEIGHT/real(WET_NZ, real64))
    !> The uniform pressure every wetting run starts from, Pa.
    real(real64), parameter :: WET_START_PRESSURE = -4.0d4
    !> The pressure its inlet holds, Pa.
    !>
    !> Below the saturation pressure of both soils, so the whole run stays on
    !> the unsaturated branch and the derivative checks are taken on a smooth
    !> state rather than across the kink where the pore volume fills.
    real(real64), parameter :: WET_INLET_PRESSURE = -8.0d3
    real(real64), parameter :: WET_DT = 600.0d0
    !> The flux unit of those runs, m/s: the order the inlet actually carries.
    real(real64), parameter :: WET_FLUX_SCALE = 1.0d-7

    ! ---- The refinement study ----

    !> How many profile windows the arrays below hold.
    integer(int32), parameter :: MAX_PROFILE_WINDOWS = 64
    !> How many refinement levels one study may hold.
    integer(int32), parameter :: MAX_REFINEMENT_LEVELS = 6

    !> The physical time every refinement level is advanced to, s.
    !>
    !> The same for the spatial and the temporal ladder, and the same as the
    !> interval the earlier `mixed_richards` case measured, so nothing but the
    !> discretisation moves between the runs being compared.
    real(real64), parameter :: REFINE_END_TIME = 6.0d0*WET_DT

    !> The spatial ladder: a uniform refinement in both directions.
    !>
    !> Doubling `nx` and `nz` together splits each rectangle into four and each
    !> coarse triangle into four fine ones, so the meshes are nested, the aspect
    !> ratio is constant, and the material interface at `WET_SPLIT` stays on cell
    !> boundaries at every level.
    !>
    !> The fifth level exists for a longer validation run and is not reached by
    !> the default `SPACE_LEVELS`; `-mixed_refine_space_levels` selects it.
    integer(int32), parameter :: MAX_SPACE_LEVELS = 5
    integer(int32), parameter :: SPACE_LEVELS = 4
    integer(int32), parameter :: SPACE_NX(MAX_SPACE_LEVELS) = [8, 16, 32, 64, 128]
    integer(int32), parameter :: SPACE_NZ(MAX_SPACE_LEVELS) = [2, 4, 8, 16, 32]
    !> The windows the spatial study projects onto: the columns of its coarsest
    !> mesh, so that level loses nothing to the projection and every finer level
    !> is an exact union of cells inside them.
    integer(int32), parameter :: SPACE_WINDOWS = 8
    !> How many steps every spatial level takes. The step is the one thing the
    !> spatial ladder holds fixed, which is what keeps the two studies separate.
    integer(int32), parameter :: SPACE_STEPS = 6

    !> The temporal ladder, all of it on one mesh.
    integer(int32), parameter :: TIME_LEVELS = 4
    integer(int32), parameter :: TIME_NX = 16, TIME_NZ = 4
    integer(int32), parameter :: TIME_WINDOWS = TIME_NX
    integer(int32), parameter :: TIME_STEPS(TIME_LEVELS) = [6, 12, 24, 48]

    !> How much a successive difference has to shrink from one mesh to the next.
    !>
    !> A criterion rather than a measured value. The homogeneous mesh ladder
    !> contracts at worst `0.31` in the measurements recorded with this test, so
    !> this leaves a factor of two while still excluding a sequence that had
    !> stopped converging, which would stagnate near one.
    real(real64), parameter :: REFINE_CONTRACTION_SPACE = 0.60d0
    !> The same for the step ladders, whose worst measured ratio is `0.67`.
    !>
    !> Looser than the mesh criterion because backward Euler is first order and
    !> because a fixed mesh limits how much of the step error is still visible in
    !> a cell average once the step is small; the observed rates are reported.
    real(real64), parameter :: REFINE_CONTRACTION_TIME = 0.80d0
    !> How large the coarsest successive difference has to be, relative to the
    !> problem's own range. Below this the ladder would be comparing a constant
    !> with itself and the contraction above would mean nothing.
    real(real64), parameter :: REFINE_SIGNAL_FLOOR = 1.0d-3
    !> How far the coarsest error has to exceed the reference's own uncertainty
    !> before the finest level is usable as the reference for the ladder.
    real(real64), parameter :: REFINE_REFERENCE_MARGIN = 4.0d0

    !> The mesh and step length one wetting run is discretised with.
    type :: type_wetting_grid
        !> Cells across the strip; even, so the material split is on a boundary.
        integer(int32) :: nx = WET_NX
        !> Cell rows up the strip.
        integer(int32) :: nz = WET_NZ
        !> The target step length, s.
        real(real64) :: dt = WET_DT
        !> When true a refused step ends the interval instead of being halved.
        !>
        !> A refinement ladder in time cannot let the solver answer with a step
        !> it chose itself: the nominal step would then not be the step that ran.
        logical :: fixed_step = .false.
    end type type_wetting_grid

    !> What one advanced wetting interval measured.
    type :: type_wetting_result
        !> Water that entered through the boundary, m2 per unit thickness.
        real(real64) :: influx = 0.0d0
        !> Worst cell defect over the water one cell can hold.
        real(real64) :: worst_inventory = 0.0d0
        !> Worst cell defect over what the step moved through it.
        real(real64) :: worst_throughput = 0.0d0
        !> Worst \(p-\max(p_D)\) over the range of the data, over all steps.
        real(real64) :: overshoot = 0.0d0
        !> Worst \(\min(p_D)-p\) over the range of the data, over all steps.
        real(real64) :: undershoot = 0.0d0
        !> Largest scaled residual norm left by an accepted step.
        real(real64) :: max_residual = 0.0d0
        !> Shortest accepted step, s.
        real(real64) :: shortest_step = 0.0d0
        !> Longest accepted step, s.
        real(real64) :: longest_step = 0.0d0
        !> Wall time of the interval, s.
        real(real64) :: wall_time = 0.0d0
        integer(int32) :: steps = 0
        integer(int32) :: rejections = 0
        !> Total SNES iterations over the accepted steps.
        integer(int32) :: newton = 0
        logical :: ok = .false.
    end type type_wetting_result

    !> One solved state projected onto the study's common windows.
    type :: type_wetting_profile
        integer(int32) :: windows = 0
        !> Area of each window, m2 per unit thickness.
        real(real64) :: volume(MAX_PROFILE_WINDOWS) = 0.0d0
        !> Volume-weighted mean cell pressure of each window, Pa.
        real(real64) :: pressure(MAX_PROFILE_WINDOWS) = 0.0d0
        !> Volume-weighted mean water content of each window, m3/m3.
        real(real64) :: water(MAX_PROFILE_WINDOWS) = 0.0d0
        !> Largest spread of cell pressure inside one window, Pa.
        !>
        !> Exactly what the projection discards. On the fixed-mesh ladder, whose
        !> windows are single cell columns, it is the departure from a
        !> one-dimensional answer; on the mesh ladder, whose windows stay the
        !> coarsest mesh's columns, it also contains the variation in x that the
        !> finer cells resolve inside one window.
        real(real64) :: spread_pressure = 0.0d0
        !> Largest spread of cell water content inside one window, m3/m3.
        real(real64) :: spread_water = 0.0d0
        !> Area of the cells left of the material split, m2.
        real(real64) :: left_volume = 0.0d0
        !> How many cells the whole mesh has.
        integer(int32) :: cells = 0
    end type type_wetting_profile

    !> One refinement level: how it was discretised and what it measured.
    type :: type_refinement_level
        type(type_wetting_grid) :: grid
        type(type_wetting_result) :: run
        type(type_wetting_profile) :: profile
        !> Water the stored quantity gained over the interval, m2.
        real(real64) :: storage_change = 0.0d0
        !> \((\Delta W-\text{influx})/|\text{influx}|\).
        real(real64) :: budget_error = 0.0d0
    end type type_refinement_level

    PetscErrorCode :: ierr
    character(len=64) :: case_name
    !> Where the measured tables are written, from `-mixed_refine_csv`; empty
    !> means no artifact, which is what a CTest run leaves behind.
    character(len=256) :: csv_prefix = ""
    !> Whether the CSV artifact already carries its header rows.
    logical :: csv_started = .false.
    integer(int32) :: failures

    call PetscInitialize(ierr)
    if (ierr /= 0) error stop "PETSc could not start"

    failures = 0
    call read_case_name(case_name)

    ! Everything that solves is verified against a direct factorisation.
    select case (trim(case_name))
    case ("darcy_constant", "hydrostatic", "layers", "flux_boundary", "jacobian", &
          "all_flux", "all_flux_incompatible", "storage_pressure", "storage")
        call configure_direct_solver()
    ! The same factorisation, with a nonlinear failure returned rather than
    ! fatal: these cases measure a rejection, or retry a refused step.
    case ("rollback", "richards", "material_interface", "refinement", "frozen_column", "mizo_column", &
          "profile_export")
        ! The same factorisation, but a failed solve has to come back as a
        ! return value rather than abort the run being measured.
        call configure_recoverable_solver()
    end select

    select case (trim(case_name))
    case ("storage")
        call case_storage(failures)
    case ("refinement")
        call case_refinement(failures)
    case ("rollback")
        call case_rollback(failures)
    case ("richards")
        call case_richards(failures)
    case ("material_interface")
        call case_material_interface(failures)
    case ("storage_pressure")
        call case_storage_pressure(failures)
    case ("all_flux")
        call case_all_flux(failures, .false.)
    case ("all_flux_incompatible")
        call case_all_flux(failures, .true.)
    case ("space")
        call case_space(failures)
    case ("rectangle")
        call case_rectangle(failures)
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
    case ("reject_warped")
        call case_reject_warped()
    case ("freeze_closure")
        call case_freeze_closure(failures)
    case ("frozen_column")
        call case_frozen_column(failures)
    case ("mizo_column")
        call case_mizo_column(failures)
    case ("profile_export")
        call case_profile_export(failures)
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
    !>
    !> Read through the options database rather than off the raw command line,
    !> so that PETSc records the option as consumed. A run that misspells it
    !> then shows up in `-options_left` instead of quietly falling back to the
    !> default case and reporting that one's result.
    subroutine read_case_name(name)
        implicit none
        character(len=*), intent(inout) :: name

        PetscBool :: found
        PetscErrorCode :: local_ierr

        name = "space"
        found = PETSC_FALSE
        call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-mixed_case", &
                                   name, found, local_ierr)
        if (local_ierr /= 0 .or. (found .neqv. PETSC_TRUE)) name = "space"
    end subroutine read_case_name

    !> A real option read from the options database.
    function read_real_option(key, fallback) result(value)
        implicit none
        character(len=*), intent(in) :: key
        real(real64), intent(in) :: fallback
        real(real64) :: value

        PetscReal :: read_value
        PetscBool :: found
        PetscErrorCode :: local_ierr

        value = fallback
        found = PETSC_FALSE
        read_value = real(fallback, PETSC_REAL_KIND)
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, key, read_value, found, local_ierr)
        if (local_ierr == 0 .and. (found .eqv. PETSC_TRUE)) value = real(read_value, real64)
    end function read_real_option

    !> Install the direct solve the steady reference cases are verified with.
    !>
    !> The mixed system is a saddle point, so it has no diagonal for a default
    !> Krylov method to work from; a case that silently fell back to GMRES with
    !> an incomplete factorisation would be measuring the preconditioner rather
    !> than the discretisation. Each setting is only installed when the run did
    !> not already ask for something else, so a command line still overrides it.
    subroutine configure_direct_solver()
        implicit none

        call set_default_option("-ksp_type", "preonly")
        call set_default_option("-pc_type", "lu")
        call set_default_option("-pc_factor_mat_solver_type", "mumps")
        call set_default_option("-ksp_error_if_not_converged", "true")
        call set_default_option("-snes_error_if_not_converged", "true")
    end subroutine configure_direct_solver

    !> The same direct factorisation, with failures reported instead of fatal.
    !>
    !> A rejection test has to survive the rejection it measures, so this leaves
    !> `error_if_not_converged` off and the case switches it on and off around
    !> the one trial that is meant to fail.
    subroutine configure_recoverable_solver()
        implicit none

        call set_default_option("-ksp_type", "preonly")
        call set_default_option("-pc_type", "lu")
        call set_default_option("-pc_factor_mat_solver_type", "mumps")
    end subroutine configure_recoverable_solver

    !> Set one option unless the run already carries it.
    subroutine set_default_option(key, value)
        implicit none
        character(len=*), intent(in) :: key
        character(len=*), intent(in) :: value

        PetscBool :: already_set
        PetscErrorCode :: local_ierr

        already_set = PETSC_FALSE
        call PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, key, already_set, local_ierr)
        if (local_ierr /= 0) return
        if (already_set .eqv. PETSC_TRUE) return
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, key, value, local_ierr)
    end subroutine set_default_option

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

    !> Record one measurement that has to be larger than a floor.
    !>
    !> The discriminating half of a regression: a quantity that has to be small
    !> says the answer is right, and a quantity that has to be large says the
    !> check could have told a wrong answer apart from it.
    subroutine check_exceeds(label, measured, floor, failures)
        implicit none
        character(len=*), intent(in) :: label
        real(real64), intent(in) :: measured
        real(real64), intent(in) :: floor
        integer(int32), intent(inout) :: failures

        character(len=160) :: line

        if (abs(measured) >= floor) then
            write (line, '(a,a,a,es12.5,a,es12.5)') "  ok   ", label, "  measured ", measured, "  floor  ", floor
        else
            write (line, '(a,a,a,es12.5,a,es12.5)') "  FAIL ", label, "  measured ", measured, "  floor  ", floor
            failures = failures + 1
        end if
        call report_line(trim(line))
    end subroutine check_exceeds

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

    !> Record one textual measurement against its expected value.
    subroutine check_text(label, measured, expected, failures)
        implicit none
        character(len=*), intent(in) :: label
        character(len=*), intent(in) :: measured
        character(len=*), intent(in) :: expected
        integer(int32), intent(inout) :: failures

        character(len=240) :: line

        if (trim(measured) == trim(expected)) then
            write (line, '(a,a,a,a,a,a)') "  ok   ", label, "  measured ", trim(measured), &
                "  expected ", trim(expected)
        else
            write (line, '(a,a,a,a,a,a)') "  FAIL ", label, "  measured ", trim(measured), &
                "  expected ", trim(expected)
            failures = failures + 1
        end if
        call report_line(trim(line))
    end subroutine check_text

    !> The linear solver the case actually ran with.
    !>
    !> The tolerances below are only meaningful against a direct factorisation,
    !> so what was used is measured rather than assumed. Reading the factor
    !> matrix back is also what establishes that this build can make a MUMPS
    !> factor, which a configuration macro on its own does not.
    subroutine check_direct_solver(solver, failures)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        integer(int32), intent(inout) :: failures

        KSP :: ksp
        PC :: preconditioner
        Mat :: factor
        character(len=64) :: ksp_type, pc_type, factor_type

        ksp_type = ""
        pc_type = ""
        factor_type = ""

        call SNESGetKSP(solver%snes, ksp, ierr)
        if (ierr == 0) call KSPGetType(ksp, ksp_type, ierr)
        if (ierr == 0) call KSPGetPC(ksp, preconditioner, ierr)
        if (ierr == 0) call PCGetType(preconditioner, pc_type, ierr)
        if (ierr == 0) call PCFactorGetMatrix(preconditioner, factor, ierr)
        if (ierr == 0) call MatGetType(factor, factor_type, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the linear solver of the run could not be inspected")
            failures = failures + 1
            return
        end if

        call check_text("linear solver of the solved step        ", &
                        trim(ksp_type)//"/"//trim(pc_type)//"/"//trim(factor_type), &
                        "preonly/lu/mumps", failures)
    end subroutine check_direct_solver

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

    !> Vertex heights of a strip, geometrically graded towards its top.
    !>
    !> \[ \Delta_j = \Delta_1\,q^{\,j-1},\qquad q = R^{1/(n_z-1)},\qquad
    !>    \sum_{j=1}^{n_z}\Delta_j = H \]
    !> where \(\Delta_j\) is the thickness of the \(j\)-th cell counted down from
    !> the top \(z = H\), so the base cell is \(R\) times as thick as the top one.
    !> Summing the progression gives \(\Delta_j = \Delta_1 + (q-1)\,d_j\) with
    !> \(d_j\) the depth of the cell's upper face: the thickness grows linearly
    !> with the distance from the top. That is the shape of the freezing fringe
    !> under a cooled top: under a millimetre in the first hours, a few
    !> millimetres once the front is deep.
    !>
    !> Assumptions: \(n_z \ge 1\), \(H > 0\), \(R > 0\); \(R < 1\) refines the base.
    !> Numerical guarantees: \(R = 1\) gives \(z_k = Hk/n_z\) exactly, the uniform
    !> grid the strip builders placed before the grading existed; \(z_0 = 0\) and
    !> \(z_{n_z} = H\) exactly for any \(R\).
    !> Computational complexity: O(n_z) arithmetic and memory.
    !> Failure behavior: error stop on a non-positive ratio, and on a ratio so
    !> extreme that round-off leaves a row without positive thickness.
    subroutine calc_strip_heights(nz, height, grading, heights)
        implicit none
        !> Cell count along the strip, \(n_z \ge 1\).
        integer(int32), intent(in) :: nz
        !> Strip height \(H\), m.
        real(real64), intent(in) :: height
        !> Ratio \(R\) of the base cell's thickness to the top cell's, \(R > 0\).
        real(real64), intent(in) :: grading
        !> Vertex heights \(z_0, \dots, z_{n_z}\), m; overwritten, size \(n_z + 1\).
        real(real64), intent(inout) :: heights(0:)

        real(real64), allocatable :: thickness(:)
        real(real64) :: ratio, total
        integer(int32) :: j

        if (.not. (grading > 0.0d0)) error stop "The strip grading ratio must be positive"
        if (grading == 1.0d0 .or. nz == 1) then
            do j = 0, nz
                heights(j) = height*real(j, real64)/real(nz, real64)
            end do
            return
        end if

        ! Relative thicknesses from the top down, scaled to the height at the
        ! end so that both ends are exact whatever the progression rounds to.
        allocate (thickness(nz))
        ratio = grading**(1.0d0/real(nz - 1, real64))
        thickness(1) = 1.0d0
        do j = 2, nz
            thickness(j) = thickness(j - 1)*ratio
        end do
        total = sum(thickness)
        heights(nz) = height
        do j = 1, nz - 1
            heights(nz - j) = heights(nz - j + 1) - height*thickness(j)/total
        end do
        heights(0) = 0.0d0
        deallocate (thickness)

        do j = 1, nz
            if (.not. (heights(j) > heights(j - 1))) error stop "The strip grading leaves a row without thickness"
        end do
    end subroutine calc_strip_heights

    !> A rectangle of right triangles, `nx` by `nz` cells.
    !>
    !> Each rectangle is cut along one diagonal, so the material interface at
    !> `x = 0.5` lies on cell boundaries whenever `nx` is even. That is what
    !> makes the layered reference solution exact rather than approximate.
    subroutine build_strip(nx, nz, width, height, dm, ierr, grading)
        implicit none
        integer(int32), intent(in) :: nx
        integer(int32), intent(in) :: nz
        real(real64), intent(in) :: width
        real(real64), intent(in) :: height
        DM, intent(inout) :: dm
        PetscErrorCode, intent(inout) :: ierr
        !> Ratio of the base row's thickness to the top row's, see
        !> `calc_strip_heights`; absent means uniform rows.
        real(real64), intent(in), optional :: grading

        PetscInt, allocatable :: cells(:)
        PetscReal, allocatable :: coordinates(:)
        real(real64), allocatable :: heights(:)
        integer(int32) :: i, k, corner, cell
        integer(int32) :: num_cells, num_vertices

        num_cells = 2*nx*nz
        num_vertices = (nx + 1)*(nz + 1)
        allocate (cells(3*num_cells), coordinates(2*num_vertices), heights(0:nz))
        if (present(grading)) then
            call calc_strip_heights(nz, height, grading, heights)
        else
            call calc_strip_heights(nz, height, 1.0d0, heights)
        end if

        do k = 0, nz
            do i = 0, nx
                corner = k*(nx + 1) + i
                coordinates(2*corner + 1) = real(width*real(i, real64)/real(nx, real64), PETSC_REAL_KIND)
                coordinates(2*corner + 2) = real(heights(k), PETSC_REAL_KIND)
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
        deallocate (cells, coordinates, heights)
    end subroutine build_strip

    !> The same rectangle, left as `nx` by `nz` rectangles instead of cut.
    !>
    !> The vertices are the ones `build_strip` places, in the same order, so a
    !> quadrilateral run and a triangular run of the same `nx` and `nz` carry
    !> the same vertex grid and therefore the same continuous temperature space.
    !> What differs between them is the cell count -- one per rectangle instead
    !> of two -- and the flux and pressure spaces those cells carry.
    subroutine build_quad_strip(nx, nz, width, height, dm, ierr, grading)
        implicit none
        integer(int32), intent(in) :: nx
        integer(int32), intent(in) :: nz
        real(real64), intent(in) :: width
        real(real64), intent(in) :: height
        DM, intent(inout) :: dm
        PetscErrorCode, intent(inout) :: ierr
        !> Ratio of the base row's thickness to the top row's, see
        !> `calc_strip_heights`; absent means uniform rows.
        real(real64), intent(in), optional :: grading

        PetscInt, allocatable :: cells(:)
        PetscReal, allocatable :: coordinates(:)
        real(real64), allocatable :: heights(:)
        integer(int32) :: i, k, corner, cell
        integer(int32) :: num_cells, num_vertices

        num_cells = nx*nz
        num_vertices = (nx + 1)*(nz + 1)
        allocate (cells(4*num_cells), coordinates(2*num_vertices), heights(0:nz))
        if (present(grading)) then
            call calc_strip_heights(nz, height, grading, heights)
        else
            call calc_strip_heights(nz, height, 1.0d0, heights)
        end if

        do k = 0, nz
            do i = 0, nx
                corner = k*(nx + 1) + i
                coordinates(2*corner + 1) = real(width*real(i, real64)/real(nx, real64), PETSC_REAL_KIND)
                coordinates(2*corner + 2) = real(heights(k), PETSC_REAL_KIND)
            end do
        end do

        cell = 0
        do k = 0, nz - 1
            do i = 0, nx - 1
                corner = k*(nx + 1) + i
                ! Counter-clockwise around the rectangle, which is the cyclic
                ! vertex order DMPlex reads a quadrilateral's geometry from.
                cells(4*cell + 1) = corner
                cells(4*cell + 2) = corner + 1
                cells(4*cell + 3) = corner + nx + 2
                cells(4*cell + 4) = corner + nx + 1
                cell = cell + 1
            end do
        end do

        call DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, 2_PETSC_INT_KIND, &
                                           int(num_cells, PETSC_INT_KIND), &
                                           int(num_vertices, PETSC_INT_KIND), 4_PETSC_INT_KIND, &
                                           PETSC_TRUE, cells, 2_PETSC_INT_KIND, coordinates, dm, ierr)
        deallocate (cells, coordinates, heights)
    end subroutine build_quad_strip

    !> One strip of either cell shape.
    subroutine build_strip_of_type(quadrilateral, nx, nz, width, height, dm, ierr, grading)
        implicit none
        !> True for rectangles, false for the triangles the strip is cut into.
        logical, intent(in) :: quadrilateral
        integer(int32), intent(in) :: nx
        integer(int32), intent(in) :: nz
        real(real64), intent(in) :: width
        real(real64), intent(in) :: height
        DM, intent(inout) :: dm
        PetscErrorCode, intent(inout) :: ierr
        !> Ratio of the base row's thickness to the top row's; absent means uniform.
        real(real64), intent(in), optional :: grading

        if (quadrilateral) then
            call build_quad_strip(nx, nz, width, height, dm, ierr, grading)
        else
            call build_strip(nx, nz, width, height, dm, ierr, grading)
        end if
    end subroutine build_strip_of_type

    !> How many cells a strip of `nx` by `nz` rectangles holds, cut or not.
    pure function strip_cell_count(quadrilateral, nx, nz) result(cells)
        implicit none
        logical, intent(in) :: quadrilateral
        integer(int32), intent(in) :: nx
        integer(int32), intent(in) :: nz
        integer(int32) :: cells

        cells = merge(1, 2, quadrilateral)*nx*nz
    end function strip_cell_count

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

        PetscInt :: facet, f_start, f_end, previous_value
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
            ! Boundary IDs form a partition. PETSc labels can otherwise retain
            ! both the previous pressure ID and the new prescribed-flux ID.
            call DMGetLabelValue(dm, "Face Sets", facet, previous_value, ierr)
            if (ierr /= 0) return
            if (previous_value /= -1_PETSC_INT_KIND) then
                call DMClearLabelValue(dm, "Face Sets", facet, previous_value, ierr)
                if (ierr /= 0) return
            end if
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
        PetscPartitioner :: partitioner
        PetscMPIInt :: comm_size

        call MPI_Comm_size(PETSC_COMM_WORLD, comm_size, ierr)
        if (comm_size <= 1) return

        ! Keep the tiny fixtures balanced, including one triangle per rank
        ! in the two-rank orientation test.
        call DMPlexGetPartitioner(dm, partitioner, ierr)
        if (ierr /= 0) return
        call PetscPartitionerSetType(partitioner, PETSCPARTITIONERSIMPLE, ierr)
        if (ierr /= 0) return
        call DMPlexDistribute(dm, 0_PETSC_INT_KIND, migration, distributed, ierr)
        if (ierr /= 0) return
        call PetscSFDestroy(migration, ierr)
        call DMDestroy(dm, ierr)
        dm = distributed
    end subroutine distribute

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
        PetscMPIInt :: comm_size

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

        call MPI_Comm_size(PETSC_COMM_WORLD, comm_size, ierr)
        call check_integer("owned fixture cells per rank            ", &
                           int(c_end - c_start, int32), 2/int(comm_size, int32), failures)
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
        PetscInt :: cell, c_start, c_end, boundary_label
        PetscInt :: facets(MAX_CELL_FACETS)
        real(real64) :: moments(MAX_CELL_FACETS)
        real(real64) :: divergence_integral, volume, outward_sum
        real(real64) :: worst_divergence, worst_constant, mismatch, total_area
        real(real64) :: expected_divergence
        integer(int32) :: num_facets, f, interior_facets, reference_facet, matched_facets(5)
        real(real64) :: base(MIXED_DIMENSION), slope
        PetscReal :: facet_measure, facet_midpoint(3), facet_normal(3)
        PetscReal :: cell_measure, cell_midpoint(3), cell_normal(3)
        real(real64) :: expected_moment, worst_moment
        ! Independent integrals of (b + s*x) dot n on AB, BC, CA, BD and DC
        ! of the fixture. BC is outward from ABC and inward to BDC.
        real(real64), parameter :: exact_midpoints(2, 5) = reshape( &
            [0.65d0, 0.1d0, 0.85d0, 0.65d0, 0.2d0, 0.55d0, 1.5d0, 0.8d0, 1.05d0, 1.25d0], [2, 5])
        real(real64), parameter :: exact_moments(5) = [1.9d-6, 5.85d-6, -3.7d-6, 8.44d-6, 1.73d-6]
        integer(int32), parameter :: expected_matches(5) = [1, 2, 1, 1, 1]

        ! The constrained edge passes through the origin, so x dot n = 0
        ! there. This affine RT0 field has the prescribed constant normal
        ! trace; inserting essential data therefore preserves the projection.
        base = [3.0d-6, -1.0d-6]
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
        worst_moment = 0.0d0
        matched_facets = 0
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

            call DMPlexComputeCellGeometryFVM(dm, cell, cell_measure, cell_midpoint, cell_normal, ierr)
            do f = 1, num_facets
                call DMPlexComputeCellGeometryFVM(dm, facets(f), facet_measure, facet_midpoint, facet_normal, ierr)
                do reference_facet = 1, 5
                    if (maxval(abs(real(facet_midpoint(1:2), real64) - &
                                   exact_midpoints(:, reference_facet))) > 1.0d-12) cycle
                    expected_moment = exact_moments(reference_facet)
                    if (reference_facet == 2 .and. cell_midpoint(1) > 1.0d0) expected_moment = -expected_moment
                    worst_moment = max(worst_moment, abs(moments(f) - expected_moment)/abs(expected_moment))
                    matched_facets(reference_facet) = matched_facets(reference_facet) + 1
                end do
                ! Local support also has size one on a partition interface.
                ! Physical boundary labels were assigned before distribution.
                call DMGetLabelValue(dm, "Face Sets", facets(f), boundary_label, ierr)
                if (boundary_label > 0) boundary_flux = boundary_flux + moments(f)
            end do
        end do

        call reduce_sum(boundary_flux)
        call reduce_sum(total_area)
        call reduce_max(worst_divergence)
        call reduce_max(worst_constant)
        call reduce_max(worst_moment)
        do reference_facet = 1, 5
            call reduce_count(matched_facets(reference_facet))
        end do

        call check_integer("analytic facet coverage mismatch        ", &
                           sum(abs(matched_facets - expected_matches)), 0, failures)
        call check("all facet moments against exact         ", worst_moment, 1.0d-12, failures)
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
        !> Every cell of the 8 by 2 strip contributes its three facets.
        integer(int32), parameter :: EXPECTED_CELL_FACETS = 3*2*8*2
        real(real64) :: left_pressure, pressure_scale, flux_scale
        real(real64) :: outflow, computed_flux, worst_pressure, worst_facet
        integer(int32) :: facets_compared
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

        call check_direct_solver(solver, failures)

        call measure_boundary_outflow(solver, 2, outflow)
        computed_flux = outflow/HEIGHT
        call check("flux through the outlet, relative error ", &
                   (computed_flux - exact_flux)/abs(exact_flux), 1.0d-8, failures)

        call compare_layer_facet_flux(solver, exact_flux, HEIGHT, worst_facet, facets_compared)
        call check_integer("facets compared with the exact flux      ", &
                           facets_compared, EXPECTED_CELL_FACETS, failures)
        ! The permissive target here is the problem's conditioning, not slack.
        ! Almost all of the head is spent across the low-conductivity layer, so
        ! the high-conductivity one carries a pressure difference smaller than
        ! its own absolute pressure by roughly the contrast. Recovering a flux
        ! from that difference loses about as many digits, measured as 5e-13 at
        ! contrast 1e3, 3e-11 at 1e5 and 2e-9 at 1e7. Reversing the layers puts
        ! the sensitive layer at the zero-pressure outlet and the same run stays
        ! at 2e-15, which is what identifies the pressure datum rather than the
        ! coefficient jump as the limit.
        call check("every facet flux against the exact field", worst_facet, 1.0d-8, failures)

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

    !> Compare every facet's normal flux with the exact uniform Darcy flux.
    !>
    !> \[
    !>   Q_{Kf}=q\,|f|\,(\mathbf e_x\cdot\mathbf n_{Kf}),
    !> \]
    !>
    !> because the exact flux of the series problem is the constant vector
    !> \((q,0)\) in both layers and across the material interface, whatever the
    !> conductivity jump is. Comparing the aggregate through the outlet alone
    !> cannot see a flux that redistributes itself inside the strip, and a facet
    !> whose normal is vertical would pass a relative test for the wrong reason,
    !> so every facet is normalised by the same reference throughput.
    subroutine compare_layer_facet_flux(solver, exact_flux, height, worst, facets_compared, radial_slope)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        !> The exact Darcy flux along x, m/s.
        real(real64), intent(in) :: exact_flux
        !> The height of the strip, m; the reference throughput is `q*height`.
        real(real64), intent(in) :: height
        !> The largest normalised disagreement.
        real(real64), intent(inout) :: worst
        !> How many cell-facet pairs were compared, over all ranks.
        integer(int32), intent(inout) :: facets_compared
        !> Optional isotropic affine part: q = (exact_flux,0) + radial_slope*x.
        real(real64), intent(in), optional :: radial_slope

        DM :: dm
        Vec :: local
        PetscInt :: cell, c_start, c_end
        PetscInt :: facets(MAX_CELL_FACETS)
        real(real64) :: moments(MAX_CELL_FACETS)
        PetscReal :: facet_measure, facet_centroid(3), facet_normal(3)
        PetscReal :: cell_measure, cell_centroid(3), cell_normal(3)
        real(real64) :: outward(MIXED_DIMENSION), expected, orientation, throughput
        integer(int32) :: num_facets, f, d

        worst = 0.0d0
        facets_compared = 0
        throughput = abs(exact_flux)*height

        dm = solver%get_dm()
        call solver%get_local_solution(local, ierr)
        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)

        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(dm, cell)) cycle
            call reconstruct_cell_facet_flux(dm, local, solver%flux_scale, cell, facets, moments, num_facets, ierr)
            if (ierr /= 0) return
            call DMPlexComputeCellGeometryFVM(dm, cell, cell_measure, cell_centroid, cell_normal, ierr)
            do f = 1, num_facets
                call DMPlexComputeCellGeometryFVM(dm, facets(f), facet_measure, facet_centroid, facet_normal, ierr)
                ! The same outward convention the reconstruction uses: away from
                ! this cell's centroid, decided by this cell alone.
                orientation = 0.0d0
                do d = 1, MIXED_DIMENSION
                    orientation = orientation + real(facet_normal(d), real64)* &
                                  (real(facet_centroid(d), real64) - real(cell_centroid(d), real64))
                end do
                outward = real(facet_normal(1:MIXED_DIMENSION), real64)
                if (orientation < 0.0d0) outward = -outward

                expected = exact_flux*real(facet_measure, real64)*outward(1)
                if (present(radial_slope)) expected = expected + radial_slope*real(facet_measure, real64)* &
                    dot_product(real(facet_centroid(1:MIXED_DIMENSION), real64), outward)
                worst = max(worst, abs(moments(f) - expected)/throughput)
                facets_compared = facets_compared + 1
            end do
        end do
        call reduce_max(worst)
        call reduce_count(facets_compared)
    end subroutine compare_layer_facet_flux

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

        call check_direct_solver(solver, failures)

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

        call check_direct_solver(solver, failures)

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

    !> Affine RT0 flux with its quadratic pressure projected to cell averages.
    subroutine case_all_flux(failures, incompatible)
        implicit none
        integer(int32), intent(inout) :: failures
        logical, intent(in) :: incompatible

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_material) :: materials(1)
        type(type_mixed_audit) :: report
        DM :: dm, solve_dm
        Vec :: local
        real(real64), parameter :: width = 1.0d0, height = 0.5d0, conductivity = 1.0d-6
        real(real64), parameter :: base = 2.0d-7
        real(real64) :: slope, mobility, pressure_scale, pressure, expected, mean, area, worst
        PetscInt :: cell, c_start, c_end
        PetscReal :: volume, centroid(3), normal(3)
        integer(int32) :: variant, facets_compared
        logical :: accepted

        call report_line("case all_flux: mean pressure gauge and integrated compatibility")
        do variant = 0, 1
            slope = real(variant, real64)*1.0d-7
            mobility = conductivity/(REFERENCE_DENSITY*STANDARD_GRAVITY)
            pressure_scale = base*width/mobility
            materials(1)%id = 1
            materials(1)%conductivity = reshape([conductivity, 0.0d0, 0.0d0, conductivity], [2, 2])
            materials(1)%source = 2.0d0*slope
            if (incompatible) materials(1)%source = materials(1)%source + base/width
            call build_strip(4, 2, width, height, dm, ierr)
            call label_cells(dm, 1.0d30, 1, 1, ierr)
            call label_facets(dm, 1, -0.01d0, 0.01d0, -1.0d0, 2.0d0, ierr)
            call label_facets(dm, 2, 0.99d0, 1.01d0, -1.0d0, 2.0d0, ierr)
            call label_facets(dm, 3, 0.005d0, 0.995d0, -0.01d0, 0.01d0, ierr)
            call label_facets(dm, 4, 0.005d0, 0.995d0, height - 0.01d0, height + 0.01d0, ierr)
            call distribute(dm, ierr)
            call solver%set_materials(materials)
            call solver%set_gravity([0.0d0, 0.0d0], REFERENCE_DENSITY)
            call solver%set_scales(base, pressure_scale, width)
            call solver%add_boundary(1, MIXED_BC_FLUX, flux=[base, 0.0d0])
            call solver%add_boundary(2, MIXED_BC_FLUX, flux=[base + slope*width, 0.0d0])
            call solver%add_boundary(3, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
            call solver%add_boundary(4, MIXED_BC_FLUX, flux=[0.0d0, slope*height])
            call solver%initialize(dm, ierr)
            call solver%set_initial_state(3.0d0*pressure_scale, ierr)
            accepted = .false.
            call solver%solve_step(0.0d0, 1.0d0, accepted, ierr)
            if (incompatible) then
                call check_integer("incompatible flux/source rejected       ", &
                                   merge(1, 0, ierr == MIXED_ERR_INCOMPATIBLE .and. .not. accepted), &
                                   1, failures)
                call solver%destroy(ierr)
                call DMDestroy(dm, ierr)
                cycle
            end if
            call check_integer("compatible all-flux solve accepted      ", merge(1, 0, accepted .and. ierr == 0), 1, failures)
            if (.not. accepted .or. ierr /= 0) return
            call check_direct_solver(solver, failures)
            solve_dm = solver%get_dm()
            call solver%get_local_solution(local, ierr)
            call DMPlexGetHeightStratum(solve_dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
            mean = 0.0d0
            area = 0.0d0
            worst = 0.0d0
            do cell = c_start, c_end - 1
                if (.not. cell_is_owned(solve_dm, cell)) cycle
                call reconstruct_cell_pressure(solve_dm, local, pressure_scale, cell, pressure, ierr)
                call DMPlexComputeCellGeometryFVM(solve_dm, cell, volume, centroid, normal, ierr)
                ! On each uniform right triangle, Var(x)+Var(z) = (hx^2+hz^2)/18.
                expected = (base*(width/2.0d0 - centroid(1)) + &
                            slope*((width**2 + height**2)/6.0d0 - &
                                   (centroid(1)**2 + centroid(2)**2 + &
                                    ((width/4.0d0)**2 + (height/2.0d0)**2)/18.0d0)/2.0d0))/mobility
                worst = max(worst, abs(pressure - expected)/pressure_scale)
                mean = mean + pressure*volume
                area = area + volume
            end do
            call reduce_max(worst)
            call reduce_sum(mean)
            call reduce_sum(area)
            call check("all-flux pressure cell averages         ", worst, 1.0d-10, failures)
            call check("volume weighted mean pressure           ", mean/(area*pressure_scale), 1.0d-12, failures)
            call compare_layer_facet_flux(solver, base, height, worst, facets_compared, slope)
            call check("all-flux analytic facet moments         ", worst, 1.0d-10, failures)
            call check_integer("all-flux facet coverage                 ", facets_compared, 48, failures)
            call solver%compute_audit(report, ierr)
            call check("all-flux local conservation             ", report%max_cell_defect/(base*height), 1.0d-10, failures)
            call check("all-flux boundary/source balance        ", &
                       (report%boundary_rate - report%source_rate)/(base*height), 1.0d-10, failures)
            call solver%destroy(ierr)
            call DMDestroy(dm, ierr)
        end do
    end subroutine case_all_flux

    !> Positive storage determines the absolute pressure even with a closed boundary.
    subroutine case_storage_pressure(failures)
        implicit none
        integer(int32), intent(inout) :: failures
        type(type_mixed_hydraulic) :: solver
        type(type_mixed_material) :: materials(1)
        type(type_mixed_audit) :: report
        DM :: dm, solve_dm
        Vec :: local
        PetscInt :: cell, first, last
        real(real64) :: pressure, worst
        logical :: accepted

        call report_line("case storage_pressure: absolute pressure set by storage, no gauge")
        materials(1)%id = 1
        materials(1)%conductivity = reshape([1.0d-6, 0.0d0, 0.0d0, 1.0d-6], [2, 2])
        materials(1)%storage = 1.0d-4
        materials(1)%source = 2.0d-4
        call build_strip(4, 2, 1.0d0, 0.5d0, dm, ierr)
        call label_cells(dm, 1.0d30, 1, 1, ierr)
        call label_facets(dm, 1, -1.0d0, 2.0d0, -1.0d0, 2.0d0, ierr)
        call distribute(dm, ierr)
        call solver%set_materials(materials)
        call solver%set_gravity([0.0d0, 0.0d0], REFERENCE_DENSITY)
        call solver%set_scales(1.0d-6, 10.0d0, 1.0d0)
        call solver%add_boundary(1, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
        call solver%initialize(dm, ierr)
        call solver%set_initial_state(3.0d0, ierr)
        accepted = .false.
        call solver%solve_step(0.0d0, 1.0d0, accepted, ierr)
        call check_integer("storage step accepted                   ", merge(1, 0, accepted .and. ierr == 0), 1, failures)
        if (.not. accepted .or. ierr /= 0) return
        call check_integer("storage excludes pressure gauge         ", merge(1, 0, solver%pressure_gauged), 0, failures)
        solve_dm = solver%get_dm()
        call solver%get_local_solution(local, ierr)
        call DMPlexGetHeightStratum(solve_dm, 0_PETSC_INT_KIND, first, last, ierr)
        worst = 0.0d0
        do cell = first, last - 1
            if (.not. cell_is_owned(solve_dm, cell)) cycle
            call reconstruct_cell_pressure(solve_dm, local, solver%pressure_scale, cell, pressure, ierr)
            worst = max(worst, abs(pressure - 5.0d0)/5.0d0)
        end do
        call reduce_max(worst)
        call check("storage sets absolute pressure          ", worst, 1.0d-10, failures)
        call solver%compute_audit(report, ierr)
        call check("storage control local conservation      ", report%max_cell_defect/1.0d-4, 1.0d-10, failures)
        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)
    end subroutine case_storage_pressure

    ! ================================================================
    ! Nonlinear storage: an independent reference closure
    ! ================================================================

    !> The stored water of the reference soil at one pore pressure.
    !>
    !> \[
    !>   W(p)=\theta_{\mathrm{VG}}\!\left(\frac{p}{\rho_{\mathrm{std}}g}\right)
    !>        +\frac{S_s}{\rho_{\mathrm{std}}g}\max(0,p),
    !> \]
    !>
    !> written out here rather than taken from the solver, so that a wrong
    !> closure cannot agree with its own reference. Strictly increasing in
    !> \(p\), which is what makes it invertible below.
    pure function reference_storage(pressure) result(stored)
        implicit none
        real(real64), intent(in) :: pressure
        real(real64) :: stored

        real(real64) :: head

        head = pressure/SPECIFIC_WEIGHT
        if (head < 0.0d0) then
            stored = SOIL_THETA_R + (SOIL_THETA_S - SOIL_THETA_R) &
                     *(1.0d0 + (-SOIL_ALPHA*head)**SOIL_N)**(-SOIL_M)
        else
            stored = SOIL_THETA_S
        end if
        if (pressure > 0.0d0) stored = stored + SOIL_COMPRESSIVE_CAPACITY*pressure
    end function reference_storage

    !> \(\partial W/\partial p\) of the reference soil, 1/Pa.
    pure function reference_capacity(pressure) result(capacity)
        implicit none
        real(real64), intent(in) :: pressure
        real(real64) :: capacity

        real(real64) :: head

        head = pressure/SPECIFIC_WEIGHT
        if (head < 0.0d0) then
            capacity = (SOIL_THETA_S - SOIL_THETA_R)*SOIL_ALPHA**SOIL_N*SOIL_M*SOIL_N &
                       *(-head)**(SOIL_N - 1.0d0) &
                       *(1.0d0 + (-SOIL_ALPHA*head)**SOIL_N)**(-SOIL_M - 1.0d0)
            capacity = capacity/SPECIFIC_WEIGHT
        else
            capacity = 0.0d0
        end if
        if (pressure >= 0.0d0) capacity = capacity + SOIL_COMPRESSIVE_CAPACITY
    end function reference_capacity

    !> The pressure whose reference stored water is `target_storage`.
    !>
    !> Bisection on a strictly increasing function over a bracket wide enough
    !> for every state these cases reach; 300 halvings leave the answer exact in
    !> double precision.
    pure function invert_reference_storage(target_storage) result(pressure)
        implicit none
        real(real64), intent(in) :: target_storage
        real(real64) :: pressure

        real(real64) :: low, high, middle
        integer(int32) :: k

        low = -1.0d9
        high = 1.0d7
        do k = 1, 300
            middle = 0.5d0*(low + high)
            if (reference_storage(middle) < target_storage) then
                low = middle
            else
                high = middle
            end if
        end do
        pressure = 0.5d0*(low + high)
    end function invert_reference_storage

    !> What a capacity-times-increment storage term would have answered.
    !>
    !> \(C(p^{n})(p-p^{n})=\Delta W\) for `at_current = .false.`, and the
    !> converged \(C(p)(p-p^{n})=\Delta W\) for `at_current = .true.`, which is
    !> the form the conservative difference is not allowed to be. Both are
    !> reference arithmetic: the point of computing them is to show that the
    !> pressure check could tell them apart.
    pure function capacity_form_pressure(pressure_old, increment, at_current) result(pressure)
        implicit none
        real(real64), intent(in) :: pressure_old
        real(real64), intent(in) :: increment
        logical, intent(in) :: at_current
        real(real64) :: pressure

        integer(int32) :: k

        pressure = pressure_old + increment/reference_capacity(pressure_old)
        if (.not. at_current) return
        do k = 1, 300
            pressure = pressure_old + increment/reference_capacity(pressure)
        end do
    end function capacity_form_pressure

    ! ================================================================
    ! Nonlinear storage: the closed fixture the two cases share
    ! ================================================================

    !> A closed strip of the reference soil, ready for its first step.
    !>
    !> Every boundary facet prescribes zero normal flux, so nothing can leave:
    !> the source and the storage alone set the pressure, and the answer is a
    !> uniform pressure with a flux of exactly zero. The flux scale is the
    !> throughput the source carries rather than an absolute number, which is
    !> what keeps the conservation row of order one at a conductivity of
    !> \(10^{-6}\) m/s.
    subroutine start_closed_retention_run(solver, dm, source, start_pressure, pressure_scale, ierr)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        DM, intent(inout) :: dm
        !> Volumetric source, 1/s; negative drains the domain.
        real(real64), intent(in) :: source
        !> Uniform initial pore pressure, Pa.
        real(real64), intent(in) :: start_pressure
        !> The pressure unit of the unknown, Pa.
        real(real64), intent(in) :: pressure_scale
        PetscErrorCode, intent(inout) :: ierr

        type(type_mixed_material) :: materials(1)

        ! A nonlinear storage needs a nonlinear solve to be measured at all, so
        ! the step has to be converged far below the quantities checked.
        call set_default_option("-snes_rtol", "1e-12")
        call set_default_option("-snes_atol", "1e-13")

        materials(1)%id = 1
        materials(1)%conductivity = reshape([SOIL_CONDUCTIVITY, 0.0d0, 0.0d0, SOIL_CONDUCTIVITY], [2, 2])
        materials(1)%source = source
        materials(1)%storage_law = MIXED_STORAGE_RETENTION
        call materials(1)%retention%reset()
        materials(1)%retention%swcc_model = SWCC_MODELS%VG
        materials(1)%retention%theta_s = SOIL_THETA_S
        materials(1)%retention%theta_r = SOIL_THETA_R
        materials(1)%retention%alpha1 = SOIL_ALPHA
        materials(1)%retention%n1 = SOIL_N
        materials(1)%retention%m1 = SOIL_M
        materials(1)%specific_storage = SOIL_SPECIFIC_STORAGE
        ! The pressure the pore volume fills at. For this unfrozen isothermal
        ! state with theta_s equal to the porosity, that is where the retention
        ! curve saturates, which every curve here does at a head of zero.
        materials(1)%saturation_pressure = 0.0d0

        call build_strip(CLOSED_NX, CLOSED_NZ, CLOSED_WIDTH, CLOSED_HEIGHT, dm, ierr)
        if (ierr == 0) call label_cells(dm, 1.0d30, 1, 1, ierr)
        if (ierr == 0) call label_facets(dm, 1, -1.0d0, 2.0d0, -1.0d0, 2.0d0, ierr)
        if (ierr == 0) call distribute(dm, ierr)
        if (ierr /= 0) return

        call solver%set_materials(materials)
        call solver%set_gravity([0.0d0, 0.0d0], REFERENCE_DENSITY)
        call solver%set_scales(max(abs(source), tiny(1.0d0))*CLOSED_WIDTH, pressure_scale, CLOSED_WIDTH)
        call solver%add_boundary(1, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
        call solver%initialize(dm, ierr)
        if (ierr == 0) call solver%set_initial_state(start_pressure, ierr)
    end subroutine start_closed_retention_run

    !> The worst relative deviation of the solved cell pressures from one value,
    !> and the lowest pressure any owned cell carries.
    subroutine compare_cell_pressures(solver, expected, pressure_error, lowest)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        !> The pressure every cell should carry, Pa; nonzero.
        real(real64), intent(in) :: expected
        real(real64), intent(inout) :: pressure_error
        real(real64), intent(inout) :: lowest

        DM :: solve_dm
        Vec :: local
        PetscInt :: cell, first, last
        real(real64) :: pressure

        pressure_error = 0.0d0
        lowest = huge(1.0d0)

        solve_dm = solver%get_dm()
        call solver%get_local_solution(local, ierr)
        call DMPlexGetHeightStratum(solve_dm, 0_PETSC_INT_KIND, first, last, ierr)
        do cell = first, last - 1
            if (.not. cell_is_owned(solve_dm, cell)) cycle
            call reconstruct_cell_pressure(solve_dm, local, solver%pressure_scale, cell, pressure, ierr)
            pressure_error = max(pressure_error, abs(pressure - expected)/abs(expected))
            lowest = min(lowest, pressure)
        end do
        call reduce_max(pressure_error)
        lowest = -lowest
        call reduce_max(lowest)
        lowest = -lowest
    end subroutine compare_cell_pressures

    !> The worst deviation of a cell's stored water from its known target.
    !>
    !> The inventory statement itself rather than the pressure one: a closed
    !> cell with a prescribed source has to end the step holding exactly what it
    !> held plus the source's contribution, and the reference closure says how
    !> much pressure that is.
    function worst_cell_inventory_error(solver, target_storage, increment) result(worst)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        !> \(W(p^{n})+s\Delta t\), the stored water the step has to reach.
        real(real64), intent(in) :: target_storage
        !> \(s\Delta t\), the size the deviation is measured against; nonzero.
        real(real64), intent(in) :: increment
        real(real64) :: worst

        DM :: solve_dm
        Vec :: local
        PetscInt :: cell, first, last
        real(real64) :: pressure

        worst = 0.0d0
        solve_dm = solver%get_dm()
        call solver%get_local_solution(local, ierr)
        call DMPlexGetHeightStratum(solve_dm, 0_PETSC_INT_KIND, first, last, ierr)
        do cell = first, last - 1
            if (.not. cell_is_owned(solve_dm, cell)) cycle
            call reconstruct_cell_pressure(solve_dm, local, solver%pressure_scale, cell, pressure, ierr)
            worst = max(worst, abs(reference_storage(pressure) - target_storage)/abs(increment))
        end do
        call reduce_max(worst)
    end function worst_cell_inventory_error

    ! ================================================================
    ! Case: conservative nonlinear storage in a closed domain
    ! ================================================================

    !> A closed cell's backward-Euler inventory, against an independent one.
    subroutine case_storage(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_audit) :: report
        DM :: dm
        real(real64), parameter :: DT = 86400.0d0
        real(real64), parameter :: UNSATURATED_SOURCE = 1.0d-7
        real(real64), parameter :: UNSATURATED_START = -5.41d4
        real(real64), parameter :: SATURATED_SOURCE = 1.0d-10
        real(real64), parameter :: SATURATED_START = 2.0d3
        real(real64), parameter :: DRAINED_SOURCE = -1.0d-10
        real(real64), parameter :: DRAINED_START = 3.0d3
        real(real64) :: increment, expected, previous, worst, lowest, gap
        real(real64) :: total_area
        character(len=12) :: tag
        integer(int32) :: step
        logical :: accepted

        call report_line("case storage: a closed cell's nonlinear backward-Euler inventory")
        total_area = CLOSED_WIDTH*CLOSED_HEIGHT

        ! ---- Unsaturated wetting, where the capacity moves within the step ----
        increment = UNSATURATED_SOURCE*DT
        call start_closed_retention_run(solver, dm, UNSATURATED_SOURCE, UNSATURATED_START, 1.0d4, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the closed retention fixture could not be built")
            failures = failures + 1
            return
        end if
        call check_integer("nonlinear storage excludes the gauge    ", &
                           merge(1, 0, solver%pressure_gauged), 0, failures)

        previous = UNSATURATED_START
        do step = 1, 3
            accepted = .false.
            call solver%solve_step(real(step, real64)*DT, DT, accepted, ierr)
            if (.not. accepted .or. ierr /= 0) then
                call report_line("  FAIL a nonlinear storage step was not accepted")
                failures = failures + 1
                return
            end if
            expected = invert_reference_storage(reference_storage(previous) + increment)
            write (tag, '(a,i0)') "step ", step
            call compare_cell_pressures(solver, expected, worst, lowest)
            call check("unsaturated pressure           "//tag, worst, 1.0d-10, failures)
            call check("unsaturated inventory          "//tag, &
                       worst_cell_inventory_error(solver, reference_storage(previous) + increment, increment), &
                       1.0d-9, failures)

            ! The same data under a capacity-times-increment storage term gives a
            ! visibly different pressure, so the check above is discriminating.
            gap = min(abs(expected - capacity_form_pressure(previous, increment, .false.)), &
                      abs(expected - capacity_form_pressure(previous, increment, .true.)))/abs(expected)
            call check_exceeds("capacity form differs by       "//tag, gap, 1.0d-3, failures)

            previous = expected
        end do

        ! The fixture is otherwise exactly what the steady gauge is for: a
        ! connected domain whose every boundary facet prescribes its flux. Only
        ! the storage magnitude keeps the gauge off, and had it not, the source
        ! would not balance the zero boundary flux and these steps would have
        ! been refused as incompatible rather than solved.
        call check_integer("nonlinear storage keeps the gauge off   ", &
                           merge(1, 0, solver%pressure_gauged), 0, failures)
        call check_direct_solver(solver, failures)
        call solver%compute_audit(report, ierr)
        call check("unsaturated local conservation          ", &
                   report%max_cell_defect/(CLOSED_CELL_AREA*SOIL_THETA_S), 1.0d-10, failures)
        call check("unsaturated domain inventory            ", &
                   report%domain_defect/(increment*total_area), 1.0d-10, failures)
        call compare_jacobian_action(solver, solver%solution, &
                                     "unsaturated storage tangent             ", 1.0d-6, failures)
        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)

        ! ---- Saturated wetting: the compressive branch is the whole capacity ----
        increment = SATURATED_SOURCE*DT
        call start_closed_retention_run(solver, dm, SATURATED_SOURCE, SATURATED_START, 1.0d3, ierr)
        previous = SATURATED_START
        do step = 1, 2
            accepted = .false.
            call solver%solve_step(real(step, real64)*DT, DT, accepted, ierr)
            if (.not. accepted .or. ierr /= 0) then
                call report_line("  FAIL a saturated storage step was not accepted")
                failures = failures + 1
                return
            end if
            ! Above the air entry the retention curve is flat, so the step is
            ! p + s dt rho g / S_s exactly, with no inversion needed.
            expected = previous + increment/SOIL_COMPRESSIVE_CAPACITY
            write (tag, '(a,i0)') "step ", step
            call compare_cell_pressures(solver, expected, worst, lowest)
            call check("saturated pressure             "//tag, worst, 1.0d-10, failures)
            call check_exceeds("saturated p stays positive "//tag, lowest, SATURATED_START, failures)
            call check("saturated inventory            "//tag, &
                       worst_cell_inventory_error(solver, reference_storage(previous) + increment, increment), &
                       1.0d-9, failures)
            previous = expected
        end do
        call solver%compute_audit(report, ierr)
        call check("saturated local conservation            ", &
                   report%max_cell_defect/(CLOSED_CELL_AREA*SOIL_THETA_S), 1.0d-10, failures)
        call compare_jacobian_action(solver, solver%solution, &
                                     "saturated storage tangent               ", 1.0d-6, failures)
        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)

        ! ---- Saturated drainage, still above the air entry ----
        increment = DRAINED_SOURCE*DT
        call start_closed_retention_run(solver, dm, DRAINED_SOURCE, DRAINED_START, 1.0d3, ierr)
        accepted = .false.
        call solver%solve_step(DT, DT, accepted, ierr)
        if (.not. accepted .or. ierr /= 0) then
            call report_line("  FAIL the saturated drainage step was not accepted")
            failures = failures + 1
            return
        end if
        expected = DRAINED_START + increment/SOIL_COMPRESSIVE_CAPACITY
        call compare_cell_pressures(solver, expected, worst, lowest)
        call check("saturated drainage pressure             ", worst, 1.0d-10, failures)
        call check_exceeds("drained pressure stays positive         ", lowest, 1.0d3, failures)
        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)
    end subroutine case_storage

    ! ================================================================
    ! Case: a rejected step restores the accepted state exactly
    ! ================================================================

    !> Copy everything one trial may change.
    subroutine snapshot_accepted_state(solver, solution, old_solution, auxiliary, step)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        Vec, intent(inout) :: solution, old_solution, auxiliary
        real(real64), intent(inout) :: step

        call VecDuplicate(solver%solution, solution, ierr)
        call VecCopy(solver%solution, solution, ierr)
        call VecDuplicate(solver%old_solution, old_solution, ierr)
        call VecCopy(solver%old_solution, old_solution, ierr)
        call VecDuplicate(solver%auxiliary, auxiliary, ierr)
        call VecCopy(solver%auxiliary, auxiliary, ierr)
        step = solver%dt
    end subroutine snapshot_accepted_state

    !> Every snapshot has to come back bit for bit, not merely close.
    !>
    !> The accepted state is what the next trial differences against, so a
    !> rejection that left it a rounding away would leave the retry solving a
    !> slightly different problem, and the inventory of the step would be wrong
    !> by that amount with nothing reporting it.
    subroutine check_state_restored(solver, solution, old_solution, auxiliary, step, label, failures)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        Vec, intent(inout) :: solution, old_solution, auxiliary
        real(real64), intent(in) :: step
        character(len=*), intent(in) :: label
        integer(int32), intent(inout) :: failures

        real(real64) :: worst

        worst = vector_difference(solver%solution, solution)
        worst = max(worst, vector_difference(solver%old_solution, old_solution))
        worst = max(worst, vector_difference(solver%auxiliary, auxiliary))
        call check(label//" state restored   ", worst, 0.0d0, failures)
        call check(label//" step restored    ", solver%dt - step, 0.0d0, failures)
    end subroutine check_state_restored

    !> The largest absolute difference between two vectors of the same layout.
    function vector_difference(left, right) result(worst)
        implicit none
        Vec, intent(in) :: left, right
        real(real64) :: worst

        Vec :: difference
        PetscReal :: measured

        call VecDuplicate(left, difference, ierr)
        call VecCopy(left, difference, ierr)
        call VecAXPY(difference, -1.0_PETSC_REAL_KIND, right, ierr)
        call VecNorm(difference, NORM_INFINITY, measured, ierr)
        worst = real(measured, real64)
        call VecDestroy(difference, ierr)
    end function vector_difference

    !> A failed trial changes nothing, and the retry at half the step works.
    subroutine case_rollback(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_audit) :: report
        DM :: dm
        KSP :: ksp
        SNESConvergedReason :: reason
        Vec :: kept_solution, kept_old, kept_aux
        real(real64), parameter :: DT = 86400.0d0
        real(real64), parameter :: SOURCE = 1.0d-7
        real(real64), parameter :: START = -5.41d4
        real(real64) :: kept_step, accepted_pressure, expected, worst, lowest
        real(real64) :: increment, nonfinite_step
        logical :: accepted

        call report_line("case rollback: a failed trial restores the accepted state exactly")
        increment = SOURCE*DT

        call start_closed_retention_run(solver, dm, SOURCE, START, 1.0d4, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the closed retention fixture could not be built")
            failures = failures + 1
            return
        end if

        ! One accepted step, so there is a nontrivial history to restore.
        accepted = .false.
        call solver%solve_step(DT, DT, accepted, ierr)
        call check_integer("first step accepted                     ", &
                           merge(1, 0, accepted .and. ierr == 0), 1, failures)
        if (.not. accepted .or. ierr /= 0) return
        accepted_pressure = invert_reference_storage(reference_storage(START) + increment)
        call compare_cell_pressures(solver, accepted_pressure, worst, lowest)
        call check("accepted pressure before rejection      ", worst, 1.0d-10, failures)
        call snapshot_accepted_state(solver, kept_solution, kept_old, kept_aux, kept_step)

        call SNESGetKSP(solver%snes, ksp, ierr)
        call SNESSetErrorIfNotConverged(solver%snes, PETSC_FALSE, ierr)
        call KSPSetErrorIfNotConverged(ksp, PETSC_FALSE, ierr)

        ! ---- The diverged-reason path: no iteration is allowed at all ----
        call SNESSetTolerances(solver%snes, PETSC_CURRENT_REAL, PETSC_CURRENT_REAL, PETSC_CURRENT_REAL, &
                               0_PETSC_INT_KIND, PETSC_CURRENT_INTEGER, ierr)
        accepted = .true.
        call solver%solve_step(2.0d0*DT, DT, accepted, ierr)
        call check_integer("diverged reason rejects the step        ", &
                           merge(1, 0, (.not. accepted) .and. ierr == 0), 1, failures)
        call SNESGetConvergedReason(solver%snes, reason, ierr)
        call check_integer("diverged reason was reported            ", &
                           merge(1, 0, int(reason%v, int32) <= 0), 1, failures)
        call check_state_restored(solver, kept_solution, kept_old, kept_aux, kept_step, &
                                  "diverged reason", failures)

        ! What the exactness is for: the balance of the last accepted step is
        ! still the balance of that step. A rollback that rewrote the storage
        ! history from the accepted solution instead would leave this reporting
        ! the source term of a step that never ran.
        call solver%compute_audit(report, ierr)
        call check("accepted balance survives rejection     ", &
                   report%max_cell_defect/(CLOSED_CELL_AREA*SOIL_THETA_S), 1.0d-10, failures)
        call compare_cell_pressures(solver, accepted_pressure, worst, lowest)
        call check("accepted pressure survives rejection    ", worst, 1.0d-10, failures)

        call SNESSetTolerances(solver%snes, PETSC_CURRENT_REAL, PETSC_CURRENT_REAL, PETSC_CURRENT_REAL, &
                               50_PETSC_INT_KIND, PETSC_CURRENT_INTEGER, ierr)

        ! ---- The PETSc error path: a linear solve that refuses to finish ----
        ! A controlled solver setting, not a fault planted in the solver: the
        ! Krylov method is given no iterations and told to treat that as an
        ! error, which is how a real linear failure reaches solve_step.
        call report_line("  the next PETSc traceback is the expected linear failure")
        call KSPSetType(ksp, KSPGMRES, ierr)
        call KSPSetTolerances(ksp, PETSC_CURRENT_REAL, PETSC_CURRENT_REAL, PETSC_CURRENT_REAL, &
                              0_PETSC_INT_KIND, ierr)
        call KSPSetErrorIfNotConverged(ksp, PETSC_TRUE, ierr)
        accepted = .true.
        call solver%solve_step(2.0d0*DT, DT, accepted, ierr)
        call check_integer("PETSc error rejects the step            ", &
                           merge(1, 0, (.not. accepted) .and. ierr /= 0), 1, failures)
        call check_state_restored(solver, kept_solution, kept_old, kept_aux, kept_step, &
                                  "PETSc error    ", failures)

        ! Back to the direct factorisation the reference is measured with.
        call KSPSetErrorIfNotConverged(ksp, PETSC_FALSE, ierr)
        call KSPSetType(ksp, KSPPREONLY, ierr)
        call KSPSetTolerances(ksp, PETSC_CURRENT_REAL, PETSC_CURRENT_REAL, PETSC_CURRENT_REAL, &
                              10000_PETSC_INT_KIND, ierr)

        ! ---- A step length that cannot be differenced over ----
        accepted = .true.
        call solver%solve_step(2.0d0*DT, 0.0d0, accepted, ierr)
        call check_integer("zero step refused                       ", &
                           merge(1, 0, (.not. accepted) .and. ierr == MIXED_ERR_INVALID_STEP), 1, failures)
        call check_state_restored(solver, kept_solution, kept_old, kept_aux, kept_step, &
                                  "zero step      ", failures)

        nonfinite_step = ieee_value(0.0d0, ieee_quiet_nan)
        accepted = .true.
        call solver%solve_step(2.0d0*DT, nonfinite_step, accepted, ierr)
        call check_integer("nonfinite step refused                  ", &
                           merge(1, 0, (.not. accepted) .and. ierr == MIXED_ERR_INVALID_STEP), 1, failures)
        call check_state_restored(solver, kept_solution, kept_old, kept_aux, kept_step, &
                                  "nonfinite step ", failures)

        ! ---- The retry, at half the step ----
        accepted = .false.
        call solver%solve_step(1.5d0*DT, 0.5d0*DT, accepted, ierr)
        call check_integer("retry at half the step accepted         ", &
                           merge(1, 0, accepted .and. ierr == 0), 1, failures)
        if (.not. accepted .or. ierr /= 0) return
        call check_direct_solver(solver, failures)
        expected = invert_reference_storage(reference_storage(accepted_pressure) + 0.5d0*increment)
        call compare_cell_pressures(solver, expected, worst, lowest)
        call check("half-step pressure after the retry      ", worst, 1.0d-10, failures)
        call check("half-step inventory after the retry     ", &
                   worst_cell_inventory_error(solver, reference_storage(accepted_pressure) + 0.5d0*increment, &
                                              0.5d0*increment), 1.0d-9, failures)
        call solver%compute_audit(report, ierr)
        call check("half-step local conservation            ", &
                   report%max_cell_defect/(CLOSED_CELL_AREA*SOIL_THETA_S), 1.0d-10, failures)

        call VecDestroy(kept_solution, ierr)
        call VecDestroy(kept_old, ierr)
        call VecDestroy(kept_aux, ierr)
        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)
    end subroutine case_rollback

    subroutine case_jacobian(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_material) :: materials(2)
        DM :: dm
        Vec :: state
        real(real64), parameter :: WIDTH = 1.0d0, HEIGHT = 0.25d0

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

        ! A dimensionless state, because that is what the residual was scaled for.
        call VecDuplicate(solver%solution, state, ierr)
        call VecSetRandom(state, PETSC_NULL_RANDOM, ierr)
        call compare_jacobian_action(solver, state, &
                                     "Jacobian action against differences     ", 1.0d-6, failures)
        call VecDestroy(state, ierr)

        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)
    end subroutine case_jacobian

    !> The assembled Jacobian action at one state, against central differences.
    !>
    !> The step is swept so that the discretisation error of the difference and
    !> the cancellation error of the subtraction can be told apart: the best of
    !> the four is reported, and a tangent that does not belong to the residual
    !> fails at every step size rather than at one.
    subroutine compare_jacobian_action(solver, state, label, tolerance, failures)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        !> The dimensionless state to linearise about.
        Vec, intent(in) :: state
        character(len=*), intent(in) :: label
        real(real64), intent(in) :: tolerance
        integer(int32), intent(inout) :: failures

        Mat :: jacobian
        Vec :: direction, product, plus, minus
        Vec :: residual_plus, residual_minus, difference
        real(real64) :: step, worst, product_norm, difference_norm
        integer(int32) :: k
        PetscReal :: norm_value

        call SNESSetUp(solver%snes, ierr)
        call DMCreateMatrix(solver%get_dm(), jacobian, ierr)

        call VecDuplicate(state, direction, ierr)
        call VecDuplicate(state, product, ierr)
        call VecDuplicate(state, plus, ierr)
        call VecDuplicate(state, minus, ierr)
        call VecDuplicate(state, residual_plus, ierr)
        call VecDuplicate(state, residual_minus, ierr)
        call VecDuplicate(state, difference, ierr)

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

            ! Each residual needs its own destination: PETSc locks the state it
            ! is given read-only, so evaluating into the perturbed state itself
            ! aborts the evaluation and leaves whatever was in the vector.
            call SNESComputeFunction(solver%snes, plus, residual_plus, ierr)
            call SNESComputeFunction(solver%snes, minus, residual_minus, ierr)
            call VecCopy(residual_plus, difference, ierr)
            call VecAXPY(difference, -1.0_PETSC_REAL_KIND, residual_minus, ierr)
            call VecScale(difference, real(1.0d0/(2.0d0*step), PETSC_REAL_KIND), ierr)
            call VecAXPY(difference, -1.0_PETSC_REAL_KIND, product, ierr)
            call VecNorm(difference, NORM_2, norm_value, ierr)
            difference_norm = real(norm_value, real64)
            worst = min(worst, difference_norm/product_norm)
        end do

        call check(label, worst, tolerance, failures)

        call VecDestroy(direction, ierr)
        call VecDestroy(product, ierr)
        call VecDestroy(plus, ierr)
        call VecDestroy(minus, ierr)
        call VecDestroy(residual_plus, ierr)
        call VecDestroy(residual_minus, ierr)
        call VecDestroy(difference, ierr)
        call MatDestroy(jacobian, ierr)
    end subroutine compare_jacobian_action

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
    ! The two soils of the nonlinear-conductivity cases
    ! ================================================================

    !> The reference van Genuchten soil, retention and conductivity together.
    !>
    !> Both curves carry the same parameters, which is the ordinary way a
    !> Mualem conductivity is paired with its retention curve; they are set
    !> separately here because the solver keeps them as two independent
    !> configurations and a test that shared one object could not tell whether
    !> the solver reads the right one.
    subroutine set_soil_material(material, id, conductivity_law)
        implicit none
        type(type_mixed_material), intent(inout) :: material
        !> The "Cell Sets" value the material is written on.
        integer(int32), intent(in) :: id
        !> MIXED_CONDUCTIVITY_FIXED or MIXED_CONDUCTIVITY_HCF.
        integer(int32), intent(in) :: conductivity_law

        material%id = id
        material%conductivity = reshape([SOIL_CONDUCTIVITY, 0.0d0, 0.0d0, SOIL_CONDUCTIVITY], [2, 2])
        material%source = 0.0d0

        material%storage_law = MIXED_STORAGE_RETENTION
        call material%retention%reset()
        material%retention%swcc_model = SWCC_MODELS%VG
        material%retention%theta_s = SOIL_THETA_S
        material%retention%theta_r = SOIL_THETA_R
        material%retention%alpha1 = SOIL_ALPHA
        material%retention%n1 = SOIL_N
        material%retention%m1 = SOIL_M
        material%specific_storage = SOIL_SPECIFIC_STORAGE
        material%saturation_pressure = 0.0d0

        material%conductivity_law = conductivity_law
        call material%conductivity_function%reset()
        material%conductivity_function%model = HCF_MODES%BASE
        material%conductivity_function%swcc_model = SWCC_MODELS%VG
        material%conductivity_function%theta_s = SOIL_THETA_S
        material%conductivity_function%theta_r = SOIL_THETA_R
        material%conductivity_function%alpha1 = SOIL_ALPHA
        material%conductivity_function%n1 = SOIL_N
        material%conductivity_function%m1 = SOIL_M
        material%conductivity_function%k_sat = SOIL_CONDUCTIVITY
        material%conductivity_function%l = SOIL_L
    end subroutine set_soil_material

    !> The Brooks-Corey soil: a different retention law and a different
    !> conductivity law, with a negative saturation head of its own.
    subroutine set_clay_material(material, id, conductivity_law)
        implicit none
        type(type_mixed_material), intent(inout) :: material
        integer(int32), intent(in) :: id
        integer(int32), intent(in) :: conductivity_law

        material%id = id
        material%conductivity = reshape([CLAY_CONDUCTIVITY, 0.0d0, 0.0d0, CLAY_CONDUCTIVITY], [2, 2])
        material%source = 0.0d0

        material%storage_law = MIXED_STORAGE_RETENTION
        call material%retention%reset()
        material%retention%swcc_model = SWCC_MODELS%BC
        material%retention%theta_s = CLAY_THETA_S
        material%retention%theta_r = CLAY_THETA_R
        material%retention%alpha1 = CLAY_ALPHA
        material%retention%n1 = CLAY_N
        material%specific_storage = CLAY_SPECIFIC_STORAGE
        ! This curve fills its pore volume at its air-entry head, not at zero.
        material%saturation_pressure = CLAY_SATURATION_PRESSURE

        material%conductivity_law = conductivity_law
        call material%conductivity_function%reset()
        material%conductivity_function%model = HCF_MODES%BASE
        material%conductivity_function%swcc_model = SWCC_MODELS%BC
        material%conductivity_function%theta_s = CLAY_THETA_S
        material%conductivity_function%theta_r = CLAY_THETA_R
        material%conductivity_function%alpha1 = CLAY_ALPHA
        material%conductivity_function%n1 = CLAY_N
        material%conductivity_function%k_sat = CLAY_CONDUCTIVITY
        material%conductivity_function%l = CLAY_L
    end subroutine set_clay_material

    !> The Brooks-Corey stored water of this file's second soil.
    !>
    !> \[
    !>   W(p)=\theta_{\mathrm{BC}}\!\left(\frac{p}{\rho_{\mathrm{std}}g}\right)
    !>        +\frac{S_s}{\rho_{\mathrm{std}}g}\max(0,p-p_{\mathrm{sat}}),
    !> \]
    !>
    !> written out here rather than taken from the solver, exactly as the van
    !> Genuchten reference above is.
    pure function reference_storage_clay(pressure) result(stored)
        implicit none
        real(real64), intent(in) :: pressure
        real(real64) :: stored

        real(real64) :: head

        head = pressure/SPECIFIC_WEIGHT
        if (head < CLAY_ALPHA) then
            stored = CLAY_THETA_R + (CLAY_THETA_S - CLAY_THETA_R)*(CLAY_ALPHA/head)**CLAY_N
        else
            stored = CLAY_THETA_S
        end if
        if (pressure >= CLAY_SATURATION_PRESSURE) then
            stored = stored + CLAY_COMPRESSIVE_CAPACITY*(pressure - CLAY_SATURATION_PRESSURE)
        end if
    end function reference_storage_clay

    !> The reference stored water of the owned cells, m2 per unit thickness.
    !>
    !> Each cell is scored with the curve of the material its centroid puts it
    !> in, from this file's own closed forms. A solver that routed the wrong
    !> retention curve to a cell would still conserve its own stored water, but
    !> it would stop agreeing with this sum, which is what makes the budget below
    !> a test of the routing rather than only of the conservation.
    function measure_reference_inventory(solver, split_x) result(total)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        !> Cells left of this use the van Genuchten curve, the rest Brooks-Corey.
        real(real64), intent(in) :: split_x
        real(real64) :: total

        DM :: dm
        Vec :: local
        PetscInt :: cell, c_start, c_end
        PetscReal :: volume, centroid(3), normal(3)
        real(real64) :: pressure

        total = 0.0d0
        dm = solver%get_dm()
        call solver%get_local_solution(local, ierr)
        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)

        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(dm, cell)) cycle
            call reconstruct_cell_pressure(dm, local, solver%pressure_scale, cell, pressure, ierr)
            call DMPlexComputeCellGeometryFVM(dm, cell, volume, centroid, normal, ierr)
            if (real(centroid(1), real64) < split_x) then
                total = total + reference_storage(pressure)*real(volume, real64)
            else
                total = total + reference_storage_clay(pressure)*real(volume, real64)
            end if
        end do
        call reduce_sum(total)
    end function measure_reference_inventory

    !> The lowest and highest pressure any owned cell carries, Pa.
    subroutine measure_pressure_range(solver, lowest, highest)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        real(real64), intent(inout) :: lowest
        real(real64), intent(inout) :: highest

        DM :: dm
        Vec :: local
        PetscInt :: cell, c_start, c_end
        real(real64) :: pressure

        lowest = huge(1.0d0)
        highest = -huge(1.0d0)
        dm = solver%get_dm()
        call solver%get_local_solution(local, ierr)
        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)

        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(dm, cell)) cycle
            call reconstruct_cell_pressure(dm, local, solver%pressure_scale, cell, pressure, ierr)
            lowest = min(lowest, pressure)
            highest = max(highest, pressure)
        end do
        call reduce_max(highest)
        lowest = -lowest
        call reduce_max(lowest)
        lowest = -lowest
    end subroutine measure_pressure_range

    ! ================================================================
    ! The wetting fixture the nonlinear-conductivity cases share
    ! ================================================================

    !> A strip wetted from a prescribed-pressure inlet at \(x=0\).
    !>
    !> Every other boundary facet prescribes zero normal flux, so the only water
    !> that enters comes through the inlet and the boundary budget is a single
    !> number. Gravity is zero: what drives the flow is the suction difference
    !> alone, which is what keeps the pressure of a cell equal to its suction and
    !> the reference curves above applicable without a hydrostatic correction.
    subroutine start_wetting_run(solver, dm, materials, grid, ierr)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        DM, intent(inout) :: dm
        !> One material, or two written either side of `WET_SPLIT`.
        type(type_mixed_material), intent(in) :: materials(:)
        !> The mesh and step this run is discretised with.
        type(type_wetting_grid), intent(in) :: grid
        PetscErrorCode, intent(inout) :: ierr

        real(real64) :: split_x

        ! A nonlinear conductivity needs a converged nonlinear solve to be
        ! measured at all, so the step is driven far below what is checked.
        call set_default_option("-snes_rtol", "1e-12")
        call set_default_option("-snes_atol", "1e-13")

        split_x = WET_SPLIT
        if (size(materials) == 1) split_x = 1.0d30

        call build_strip(grid%nx, grid%nz, WET_WIDTH, WET_HEIGHT, dm, ierr)
        if (ierr == 0) call label_cells(dm, split_x, materials(1)%id, materials(size(materials))%id, ierr)
        ! Closed everywhere, then the inlet labelled over the left edge.
        if (ierr == 0) call label_facets(dm, 2, -1.0d0, 2.0d0, -1.0d0, 2.0d0, ierr)
        if (ierr == 0) call label_facets(dm, 1, -0.001d0, 0.001d0, -1.0d0, 2.0d0, ierr)
        if (ierr == 0) call distribute(dm, ierr)
        if (ierr /= 0) return

        call solver%set_materials(materials)
        call solver%set_gravity([0.0d0, 0.0d0], REFERENCE_DENSITY)
        call solver%set_scales(WET_FLUX_SCALE, 1.0d4, WET_WIDTH)
        call solver%add_boundary(1, MIXED_BC_PRESSURE, pressure=WET_INLET_PRESSURE)
        call solver%add_boundary(2, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
        call solver%initialize(dm, ierr)
        if (ierr == 0) call solver%set_initial_state(WET_START_PRESSURE, ierr)
    end subroutine start_wetting_run

    !> Advance a wetting run to `end_time`, halving the step when one is refused.
    !>
    !> `influx` accumulates \(-\Delta t\sum_f Q_f\) over the domain boundary of
    !> each accepted step, which for backward Euler is exactly the water that
    !> step transported. The cell defects are normalised twice, once by the
    !> inventory a cell can hold and once by what the step actually moved
    !> through it.
    !>
    !> The retry is the solver's own rejection path: a refused trial restores
    !> the accepted state exactly, so the shorter step solves the same problem
    !> over a shorter interval and the accumulated budget is still the budget of
    !> one interval. Every run compared with another is advanced to the same end
    !> time, whatever step lengths it needed to get there. Nothing here changes
    !> the model: a degenerate relative conductivity is met with a shorter step,
    !> not with a modified law.
    !>
    !> A grid that asked for a fixed step gets no retry at all. The accepted
    !> lengths are recorded either way, so a study can say what step actually
    !> ran rather than what it asked for.
    subroutine run_wetting_interval(solver, grid, end_time, measured)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        !> The mesh and target step of the run; only the step is read here.
        type(type_wetting_grid), intent(in) :: grid
        !> The physical time the run has to reach, s.
        real(real64), intent(in) :: end_time
        !> Everything the interval measured; overwritten on every path.
        type(type_wetting_result), intent(inout) :: measured

        type(type_mixed_audit) :: report
        !> What fraction of the target step the interval may fall back to.
        real(real64), parameter :: MIN_STEP_FRACTION = 1.0d-3
        integer(int32), parameter :: MAX_TRIALS = 2000
        Vec :: residual
        PetscReal :: residual_norm
        PetscInt :: iterations
        real(real64) :: time, step, cell_area, started
        real(real64) :: lowest, highest, data_low, data_high, data_span
        integer(int32) :: trial
        logical :: accepted

        measured = type_wetting_result()
        measured%shortest_step = huge(1.0d0)

        cell_area = wetting_cell_area(grid)
        data_low = min(WET_START_PRESSURE, WET_INLET_PRESSURE)
        data_high = max(WET_START_PRESSURE, WET_INLET_PRESSURE)
        data_span = data_high - data_low

        call VecDuplicate(solver%solution, residual, ierr)
        started = MPI_Wtime()

        time = 0.0d0
        step = grid%dt
        do trial = 1, MAX_TRIALS
            if (time >= end_time*(1.0d0 - 1.0d-12)) then
                measured%ok = .true.
                exit
            end if
            step = min(step, end_time - time)

            accepted = .false.
            call solver%solve_step(time + step, step, accepted, ierr)
            if (.not. accepted) then
                measured%rejections = measured%rejections + 1
                ! A fixed-step ladder may not quietly answer with another step.
                if (grid%fixed_step) exit
                step = 0.5d0*step
                if (step < MIN_STEP_FRACTION*grid%dt) exit
                cycle
            end if
            if (ierr /= 0) exit

            call solver%compute_audit(report, ierr)
            if (ierr /= 0) exit
            measured%influx = measured%influx - step*report%boundary_rate
            measured%worst_inventory = max(measured%worst_inventory, &
                                           report%max_cell_defect/(cell_area*SOIL_THETA_S))
            if (report%throughput > 0.0d0) then
                measured%worst_throughput = max(measured%worst_throughput, &
                                                report%max_cell_defect/(step*report%throughput))
            end if

            ! Overshoot is a separate question from the mass budget: a scheme
            ! can conserve exactly and still put a cell outside the range of its
            ! own data, so it is measured at every step rather than at the end.
            call measure_pressure_range(solver, lowest, highest)
            measured%overshoot = max(measured%overshoot, (highest - data_high)/data_span)
            measured%undershoot = max(measured%undershoot, (data_low - lowest)/data_span)

            call SNESGetIterationNumber(solver%snes, iterations, ierr)
            if (ierr == 0) measured%newton = measured%newton + int(iterations, int32)
            ! Into a vector of its own: PETSc locks the state read-only, so a
            ! residual evaluated into its own input would abort the run.
            call SNESComputeFunction(solver%snes, solver%solution, residual, ierr)
            if (ierr == 0) call VecNorm(residual, NORM_2, residual_norm, ierr)
            if (ierr == 0) measured%max_residual = max(measured%max_residual, real(residual_norm, real64))

            time = time + step
            measured%steps = measured%steps + 1
            measured%shortest_step = min(measured%shortest_step, step)
            measured%longest_step = max(measured%longest_step, step)
            ! Grow back towards the target step, so a single hard step does not
            ! decide the resolution of the whole interval.
            step = min(grid%dt, 2.0d0*step)
        end do

        measured%wall_time = MPI_Wtime() - started
        if (measured%steps == 0) measured%shortest_step = 0.0d0
        call VecDestroy(residual, ierr)
    end subroutine run_wetting_interval

    !> The default wetting discretisation: the mesh and step the earlier Task 3
    !> cases were measured on, and the one they still run at.
    function default_wetting_grid() result(grid)
        implicit none
        type(type_wetting_grid) :: grid

        grid%nx = WET_NX
        grid%nz = WET_NZ
        grid%dt = WET_DT
        grid%fixed_step = .false.
    end function default_wetting_grid

    !> The area of one cell of a wetting grid, m2 per unit thickness.
    pure function wetting_cell_area(grid) result(area)
        implicit none
        type(type_wetting_grid), intent(in) :: grid
        real(real64) :: area

        area = 0.5d0*(WET_WIDTH/real(grid%nx, real64))*(WET_HEIGHT/real(grid%nz, real64))
    end function wetting_cell_area

    !> How much water a control run of the same fixture lets in, m2.
    !>
    !> A whole run of its own, because the closure tables are held by one live
    !> solver at a time: the control is built after the run it is compared with
    !> has been destroyed.
    function measure_control_influx(materials, end_time, ok) result(influx)
        implicit none
        type(type_mixed_material), intent(in) :: materials(:)
        real(real64), intent(in) :: end_time
        logical, intent(inout) :: ok
        real(real64) :: influx

        type(type_mixed_hydraulic) :: solver
        type(type_wetting_result) :: measured
        DM :: dm
        character(len=96) :: line

        influx = 0.0d0
        ok = .false.
        call start_wetting_run(solver, dm, materials, default_wetting_grid(), ierr)
        if (ierr /= 0) return
        call run_wetting_interval(solver, default_wetting_grid(), end_time, measured)
        influx = measured%influx
        ok = measured%ok
        write (line, '(a,i0,a,i0)') "    the control took steps ", measured%steps, &
            ", refused ", measured%rejections
        call report_line(trim(line))
        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)
    end function measure_control_influx

    ! ================================================================
    ! The mobility tangent, measured with residual evaluations alone
    ! ================================================================

    !> A global vector that is one on every pressure degree of freedom.
    subroutine build_pressure_mask(solver, mask)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        Vec, intent(inout) :: mask

        DM :: dm
        Vec :: local
        PetscSection :: section
        PetscScalar, pointer :: raw(:)
        PetscInt :: cell, c_start, c_end, offset, num_dofs

        dm = solver%get_dm()
        call VecDuplicate(solver%solution, mask, ierr)
        call VecZeroEntries(mask, ierr)

        call DMGetLocalVector(dm, local, ierr)
        call VecZeroEntries(local, ierr)
        call DMGetLocalSection(dm, section, ierr)
        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)

        call VecGetArray(local, raw, ierr)
        do cell = c_start, c_end - 1
            call PetscSectionGetFieldDof(section, cell, 1_PETSC_INT_KIND, num_dofs, ierr)
            if (ierr /= 0 .or. num_dofs <= 0) cycle
            call PetscSectionGetFieldOffset(section, cell, 1_PETSC_INT_KIND, offset, ierr)
            if (ierr == 0) raw(offset + 1) = 1.0_PETSC_REAL_KIND
        end do
        call VecRestoreArray(local, raw, ierr)

        call DMLocalToGlobalBegin(dm, local, INSERT_VALUES, mask, ierr)
        call DMLocalToGlobalEnd(dm, local, INSERT_VALUES, mask, ierr)
        call DMRestoreLocalVector(dm, local, ierr)
    end subroutine build_pressure_mask

    !> \((F(u+h\,d)-F(u-h\,d))/(2h)\), into `difference`.
    subroutine residual_difference(solver, base, direction, step, difference)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        Vec, intent(in) :: base, direction
        real(real64), intent(in) :: step
        Vec, intent(inout) :: difference

        Vec :: plus, minus, residual_plus, residual_minus

        call VecDuplicate(base, plus, ierr)
        call VecDuplicate(base, minus, ierr)
        call VecDuplicate(base, residual_plus, ierr)
        call VecDuplicate(base, residual_minus, ierr)

        call VecCopy(base, plus, ierr)
        call VecAXPY(plus, real(step, PETSC_REAL_KIND), direction, ierr)
        call VecCopy(base, minus, ierr)
        call VecAXPY(minus, real(-step, PETSC_REAL_KIND), direction, ierr)

        call SNESComputeFunction(solver%snes, plus, residual_plus, ierr)
        call SNESComputeFunction(solver%snes, minus, residual_minus, ierr)
        call VecCopy(residual_plus, difference, ierr)
        call VecAXPY(difference, -1.0_PETSC_REAL_KIND, residual_minus, ierr)
        call VecScale(difference, real(1.0d0/(2.0d0*step), PETSC_REAL_KIND), ierr)

        call VecDestroy(plus, ierr)
        call VecDestroy(minus, ierr)
        call VecDestroy(residual_plus, ierr)
        call VecDestroy(residual_minus, ierr)
    end subroutine residual_difference

    !> What fraction of the pressure sensitivity is the mobility tangent.
    !>
    !> Measured from residual evaluations only, so it is independent of the
    !> assembled Jacobian and can say whether the agreement with it is
    !> discriminating. The transport term \(-(\chi,\nabla\cdot\mathbf v)\) and the
    !> gravity term do not involve the flux and the resistance is linear in it,
    !> so differencing the same pressure-only direction once at the solved state
    !> and once with the flux set to zero leaves exactly
    !> \(\partial_p(\mathbf\Lambda^{-1}\mathbf q)\). A Jacobian without that block
    !> would therefore be wrong by this fraction.
    subroutine measure_mobility_share(solver, state, share)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        !> The dimensionless state to difference about.
        Vec, intent(in) :: state
        real(real64), intent(inout) :: share

        Vec :: mask, direction, flux_free, full, without
        real(real64), parameter :: STEP = 1.0d-5
        PetscReal :: full_norm, mobility_norm

        share = 0.0d0
        call build_pressure_mask(solver, mask)
        call VecDuplicate(state, direction, ierr)
        call VecDuplicate(state, flux_free, ierr)
        call VecDuplicate(state, full, ierr)
        call VecDuplicate(state, without, ierr)

        ! A direction that moves the pressure only, and the same state with its
        ! flux removed. Both come from the same mask, so the comparison below
        ! differs in the flux and in nothing else.
        call VecSetRandom(direction, PETSC_NULL_RANDOM, ierr)
        call VecPointwiseMult(direction, direction, mask, ierr)
        call VecPointwiseMult(flux_free, state, mask, ierr)

        call residual_difference(solver, state, direction, STEP, full)
        call residual_difference(solver, flux_free, direction, STEP, without)

        call VecNorm(full, NORM_2, full_norm, ierr)
        call VecAXPY(full, -1.0_PETSC_REAL_KIND, without, ierr)
        call VecNorm(full, NORM_2, mobility_norm, ierr)
        if (real(full_norm, real64) > 0.0d0) share = real(mobility_norm, real64)/real(full_norm, real64)

        call VecDestroy(mask, ierr)
        call VecDestroy(direction, ierr)
        call VecDestroy(flux_free, ierr)
        call VecDestroy(full, ierr)
        call VecDestroy(without, ierr)
    end subroutine measure_mobility_share

    ! ================================================================
    ! Case: unsaturated flow with a pressure-dependent conductivity
    ! ================================================================

    !> A converged wetting run whose conductivity moves with the pressure.
    subroutine case_richards(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_material) :: materials(1), fixed_materials(1)
        type(type_mixed_audit) :: report
        type(type_wetting_result) :: measured
        DM :: dm
        Vec :: state, kept_solution, kept_old, kept_aux, kept_boundary
        SNESConvergedReason :: reason
        integer(int32), parameter :: STEPS = 6
        real(real64), parameter :: END_TIME = real(STEPS, real64)*WET_DT
        real(real64) :: start_inventory, inventory, influx, control_influx
        real(real64) :: share
        real(real64) :: lowest, highest, kept_step
        logical :: accepted, accepted_all, control_ok

        call report_line("case richards: unsaturated flow with a pressure-dependent conductivity")

        call set_soil_material(materials(1), 1, MIXED_CONDUCTIVITY_HCF)
        call start_wetting_run(solver, dm, materials, default_wetting_grid(), ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the wetting fixture could not be built")
            failures = failures + 1
            return
        end if

        start_inventory = measure_reference_inventory(solver, 1.0d30)
        call run_wetting_interval(solver, default_wetting_grid(), END_TIME, measured)
        influx = measured%influx
        accepted_all = measured%ok
        call check_integer("every wetting step converged            ", &
                           merge(1, 0, accepted_all), 1, failures)
        ! This soil does not need the retry: the interval is the requested one,
        ! step for step, which is what makes the control runs comparable.
        call check_integer("steps the interval needed               ", measured%steps, STEPS, failures)
        call check_integer("steps the interval refused              ", measured%rejections, 0, failures)
        if (.not. accepted_all) then
            call solver%destroy(ierr)
            call DMDestroy(dm, ierr)
            return
        end if
        call check_direct_solver(solver, failures)

        ! ---- The physical budget: what came in is what is stored ----
        inventory = measure_reference_inventory(solver, 1.0d30)
        call check("inventory against the boundary budget   ", &
                   (inventory - start_inventory - influx)/abs(influx), 1.0d-10, failures)
        call check_exceeds("water the run actually transported      ", &
                           influx/(WET_WIDTH*WET_HEIGHT*SOIL_THETA_S), 1.0d-2, failures)
        call check("local conservation over cell inventory  ", measured%worst_inventory, 1.0d-10, failures)
        call check("local conservation over throughput      ", measured%worst_throughput, 1.0d-8, failures)

        ! The front has to be inside the strip: a run that had already filled it
        ! would be a constant state, and a constant state has a constant
        ! conductivity whatever law produced it.
        call measure_pressure_range(solver, lowest, highest)
        call check_exceeds("suction range the front spans, Pa       ", &
                           highest - lowest, 1.0d4, failures)

        ! ---- The tangent of the mobility, and whether the check can see it ----
        call VecDuplicate(solver%solution, state, ierr)
        call VecCopy(solver%solution, state, ierr)
        call measure_mobility_share(solver, state, share)
        call check_exceeds("mobility share of the pressure tangent  ", share, 1.0d-2, failures)
        call compare_jacobian_action(solver, state, &
                                     "Jacobian action against differences     ", 1.0d-6, failures)
        call VecDestroy(state, ierr)

        ! ---- A rejected step leaves the pressure boundary's history alone ----
        call check_integer("one natural boundary carries a history  ", &
                           solver%num_boundary_aux, 1, failures)
        call snapshot_accepted_state(solver, kept_solution, kept_old, kept_aux, kept_step)
        call VecDuplicate(solver%boundary_aux(1), kept_boundary, ierr)
        call VecCopy(solver%boundary_aux(1), kept_boundary, ierr)

        call SNESSetErrorIfNotConverged(solver%snes, PETSC_FALSE, ierr)
        call SNESSetTolerances(solver%snes, PETSC_CURRENT_REAL, PETSC_CURRENT_REAL, PETSC_CURRENT_REAL, &
                               0_PETSC_INT_KIND, PETSC_CURRENT_INTEGER, ierr)
        accepted = .true.
        call solver%solve_step(real(STEPS + 1, real64)*WET_DT, WET_DT, accepted, ierr)
        call check_integer("diverged reason rejects the step        ", &
                           merge(1, 0, (.not. accepted) .and. ierr == 0), 1, failures)
        call SNESGetConvergedReason(solver%snes, reason, ierr)
        call check_integer("diverged reason was reported            ", &
                           merge(1, 0, int(reason%v, int32) <= 0), 1, failures)
        call check_state_restored(solver, kept_solution, kept_old, kept_aux, kept_step, &
                                  "diverged reason", failures)
        ! The one the earlier review found untested: a natural pressure boundary
        ! keeps its own copy of the storage history, and a rejection has to leave
        ! that copy exactly as the accepted step left it too.
        call check("boundary history restored               ", &
                   vector_difference(solver%boundary_aux(1), kept_boundary), 0.0d0, failures)

        call SNESSetTolerances(solver%snes, PETSC_CURRENT_REAL, PETSC_CURRENT_REAL, PETSC_CURRENT_REAL, &
                               50_PETSC_INT_KIND, PETSC_CURRENT_INTEGER, ierr)
        accepted = .false.
        call solver%solve_step(real(STEPS, real64)*WET_DT + 0.5d0*WET_DT, 0.5d0*WET_DT, accepted, ierr)
        call check_integer("retry at half the step accepted         ", &
                           merge(1, 0, accepted .and. ierr == 0), 1, failures)
        call solver%compute_audit(report, ierr)
        call check("half-step local conservation            ", &
                   report%max_cell_defect/(WET_CELL_AREA*SOIL_THETA_S), 1.0d-10, failures)

        call VecDestroy(kept_solution, ierr)
        call VecDestroy(kept_old, ierr)
        call VecDestroy(kept_aux, ierr)
        call VecDestroy(kept_boundary, ierr)
        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)

        ! ---- The control: the same fixture with the conductivity frozen ----
        ! The saturated tensor is the same one, so the only difference is
        ! whether the relative conductivity of the pressure is evaluated at all.
        call set_soil_material(fixed_materials(1), 1, MIXED_CONDUCTIVITY_FIXED)
        control_influx = measure_control_influx(fixed_materials, END_TIME, control_ok)
        call check_integer("the frozen-conductivity control ran     ", &
                           merge(1, 0, control_ok), 1, failures)
        if (control_ok) then
            call check_exceeds("frozen conductivity moves other water   ", &
                               (control_influx - influx)/abs(influx), 1.0d0, failures)
        end if
    end subroutine case_richards

    ! ================================================================
    ! Case: two materials, two retention curves, two conductivity curves
    ! ================================================================

    !> Both laws routed per cell, in the saturated limit and while wetting.
    subroutine case_material_interface(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        call report_line("case material_interface: a retention and a conductivity curve per cell")
        call check_saturated_interface(failures)
        call check_wetting_interface(failures)
    end subroutine case_material_interface

    !> The saturated limit: both curves report one and the answer is analytic.
    !>
    !> \[
    !>   q=\frac{\Delta p/(\rho g)}{L_1/K_1+L_2/K_2},
    !> \]
    !>
    !> with \(K_i\) the saturated tensors. The flux alone is symmetric under
    !> swapping the two conductivities, so the pressure profile is compared as
    !> well: that is what says which side carries which material.
    subroutine check_saturated_interface(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_material) :: materials(2)
        type(type_mixed_audit) :: report
        DM :: dm
        real(real64), parameter :: LEFT_PRESSURE = 2.0d4, RIGHT_PRESSURE = 1.0d4
        !> Every cell of the strip contributes its three facets.
        integer(int32), parameter :: EXPECTED_CELL_FACETS = 3*2*WET_NX*WET_NZ
        real(real64) :: exact_flux, head, outflow, worst_facet, worst_pressure
        integer(int32) :: facets_compared
        logical :: accepted

        call report_line("  saturated limit, against the analytic series resistance")

        call set_soil_material(materials(1), 3, MIXED_CONDUCTIVITY_HCF)
        call set_clay_material(materials(2), 7, MIXED_CONDUCTIVITY_HCF)
        ! Steady: a saturated cell of either soil stores nothing more, so the
        ! backward-Euler difference of its stored water is exactly zero.
        materials(1)%specific_storage = 0.0d0
        materials(2)%specific_storage = 0.0d0

        head = (LEFT_PRESSURE - RIGHT_PRESSURE)/SPECIFIC_WEIGHT
        exact_flux = head/(WET_SPLIT/SOIL_CONDUCTIVITY + (WET_WIDTH - WET_SPLIT)/CLAY_CONDUCTIVITY)

        call build_strip(WET_NX, WET_NZ, WET_WIDTH, WET_HEIGHT, dm, ierr)
        call label_cells(dm, WET_SPLIT, 3, 7, ierr)
        call label_facets(dm, 3, -1.0d0, 2.0d0, -1.0d0, 2.0d0, ierr)
        call label_facets(dm, 1, -0.001d0, 0.001d0, -1.0d0, 2.0d0, ierr)
        call label_facets(dm, 2, WET_WIDTH - 0.001d0, WET_WIDTH + 0.001d0, -1.0d0, 2.0d0, ierr)
        call distribute(dm, ierr)

        call solver%set_materials(materials)
        call solver%set_gravity([0.0d0, 0.0d0], REFERENCE_DENSITY)
        call solver%set_scales(abs(exact_flux), abs(LEFT_PRESSURE), WET_WIDTH)
        call solver%add_boundary(1, MIXED_BC_PRESSURE, pressure=LEFT_PRESSURE)
        call solver%add_boundary(2, MIXED_BC_PRESSURE, pressure=RIGHT_PRESSURE)
        call solver%add_boundary(3, MIXED_BC_FLUX, flux=[0.0d0, 0.0d0])
        call solver%initialize(dm, ierr)
        call solver%set_initial_state(LEFT_PRESSURE, ierr)

        accepted = .false.
        call solver%solve_step(0.0d0, 1.0d0, accepted, ierr)
        call check_integer("the saturated interface solve converged ", &
                           merge(1, 0, accepted .and. ierr == 0), 1, failures)
        if (.not. accepted) then
            call solver%destroy(ierr)
            call DMDestroy(dm, ierr)
            return
        end if
        call check_direct_solver(solver, failures)

        call measure_boundary_outflow(solver, 2, outflow)
        call check("flux through the outlet, relative error ", &
                   (outflow/WET_HEIGHT - exact_flux)/abs(exact_flux), 1.0d-8, failures)

        call compare_layer_facet_flux(solver, exact_flux, WET_HEIGHT, worst_facet, facets_compared)
        call check_integer("facets compared with the exact flux     ", &
                           facets_compared, EXPECTED_CELL_FACETS, failures)
        call check("every facet flux against the exact field", worst_facet, 1.0d-8, failures)

        call compare_interface_pressures(solver, exact_flux, LEFT_PRESSURE, RIGHT_PRESSURE, worst_pressure)
        call check("cell pressure against cell average      ", worst_pressure, 1.0d-8, failures)

        call solver%compute_audit(report, ierr)
        call check("cell defect over throughput             ", &
                   report%max_cell_defect/report%throughput, 1.0d-10, failures)
        call check("facet mismatch over throughput          ", &
                   report%max_facet_mismatch/report%throughput, 1.0d-12, failures)

        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)
    end subroutine check_saturated_interface

    !> Cell pressures against the piecewise-linear profile of the two layers.
    subroutine compare_interface_pressures(solver, exact_flux, left_pressure, right_pressure, worst)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        real(real64), intent(in) :: exact_flux, left_pressure, right_pressure
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
            ! Linear in each layer, so the cell average is the centroid value.
            x = real(centroid(1), real64)
            if (x < WET_SPLIT) then
                expected = left_pressure - SPECIFIC_WEIGHT*exact_flux*x/SOIL_CONDUCTIVITY
            else
                expected = right_pressure + SPECIFIC_WEIGHT*exact_flux*(WET_WIDTH - x)/CLAY_CONDUCTIVITY
            end if
            worst = max(worst, abs(pressure - expected)/abs(left_pressure - right_pressure))
        end do
        call reduce_max(worst)
    end subroutine compare_interface_pressures

    !> Unsaturated wetting across the interface, with both laws per cell.
    subroutine check_wetting_interface(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_material) :: materials(2), control(1)
        type(type_wetting_result) :: measured
        DM :: dm
        ! Long enough that the front is well inside the second soil: while it is
        ! still in the first one the run cannot be told from the control made of
        ! that soil alone, and the comparison below would prove nothing.
        integer(int32), parameter :: STEPS = 12
        real(real64), parameter :: END_TIME = real(STEPS, real64)*WET_DT
        real(real64) :: start_inventory, inventory, influx
        real(real64) :: soil_influx, clay_influx
        real(real64) :: lowest, highest
        character(len=96) :: line
        logical :: accepted_all, control_ok

        call report_line("  unsaturated wetting, with a different law each side")

        call set_soil_material(materials(1), 3, MIXED_CONDUCTIVITY_HCF)
        call set_clay_material(materials(2), 7, MIXED_CONDUCTIVITY_HCF)

        call start_wetting_run(solver, dm, materials, default_wetting_grid(), ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the two-material wetting fixture could not be built")
            failures = failures + 1
            return
        end if

        start_inventory = measure_reference_inventory(solver, WET_SPLIT)
        call run_wetting_interval(solver, default_wetting_grid(), END_TIME, measured)
        influx = measured%influx
        accepted_all = measured%ok
        write (line, '(a,i0,a,i0)') "    the two-material run took steps ", measured%steps, &
            ", refused ", measured%rejections
        call report_line(trim(line))
        call check_integer("every two-material step converged       ", &
                           merge(1, 0, accepted_all), 1, failures)
        if (.not. accepted_all) then
            call solver%destroy(ierr)
            call DMDestroy(dm, ierr)
            return
        end if

        ! Each cell is scored with its own curve, so this only balances if the
        ! solver routed the same curve to the same cell.
        inventory = measure_reference_inventory(solver, WET_SPLIT)
        call check("two-curve inventory against the budget  ", &
                   (inventory - start_inventory - influx)/abs(influx), 1.0d-10, failures)
        call check("local conservation over cell inventory  ", measured%worst_inventory, 1.0d-10, failures)
        call check("local conservation over throughput      ", measured%worst_throughput, 1.0d-8, failures)

        call measure_pressure_range(solver, lowest, highest)
        call check_exceeds("suction range the front spans, Pa       ", &
                           highest - lowest, 1.0d4, failures)

        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)

        ! ---- Controls: the same fixture made of one soil at a time ----
        ! The conductivity of the unsaturated state is what these separate. A
        ! run that answered with one material's curve everywhere would land on
        ! one of these two numbers instead of between them.
        call set_soil_material(control(1), 3, MIXED_CONDUCTIVITY_HCF)
        soil_influx = measure_control_influx(control, END_TIME, control_ok)
        call check_integer("the van Genuchten control ran           ", &
                           merge(1, 0, control_ok), 1, failures)
        call set_clay_material(control(1), 7, MIXED_CONDUCTIVITY_HCF)
        clay_influx = measure_control_influx(control, END_TIME, control_ok)
        call check_integer("the Brooks-Corey control ran            ", &
                           merge(1, 0, control_ok), 1, failures)

        ! Measured at 1.9e-1 and -9.98e-1 over this interval; the floors leave
        ! room for the step pattern to differ between rank counts without
        ! leaving room for the two runs to be the same run.
        call check_exceeds("against the van Genuchten control       ", &
                           (soil_influx - influx)/abs(influx), 5.0d-2, failures)
        call check_exceeds("against the Brooks-Corey control        ", &
                           (clay_influx - influx)/abs(influx), 5.0d-1, failures)
    end subroutine check_wetting_interface

    ! ================================================================
    ! Case: mesh and time refinement of the wetting front
    ! ================================================================

    !> An integer option read from the options database.
    function read_integer_option(key, fallback) result(value)
        implicit none
        character(len=*), intent(in) :: key
        integer(int32), intent(in) :: fallback
        integer(int32) :: value

        PetscInt :: read_value
        PetscBool :: found
        PetscErrorCode :: local_ierr

        value = fallback
        found = PETSC_FALSE
        read_value = int(fallback, PETSC_INT_KIND)
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, key, read_value, found, local_ierr)
        if (local_ierr == 0 .and. (found .eqv. PETSC_TRUE)) value = int(read_value, int32)
    end function read_integer_option

    !> A string option read from the options database.
    subroutine read_text_option(key, value, fallback)
        implicit none
        character(len=*), intent(in) :: key
        character(len=*), intent(inout) :: value
        character(len=*), intent(in) :: fallback

        PetscBool :: found
        PetscErrorCode :: local_ierr

        value = fallback
        found = PETSC_FALSE
        call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, key, value, found, local_ierr)
        if (local_ierr /= 0 .or. (found .neqv. PETSC_TRUE)) value = fallback
    end subroutine read_text_option

    !> Which cell shape a run asked one of the column meshes to be built from.
    !>
    !> The two shapes are named rather than toggled so that a run says which
    !> discretisation it measured, and an unrecognised name stops the run
    !> instead of quietly selecting the default and reporting its numbers.
    function read_cell_type_option(key) result(quadrilateral)
        implicit none
        !> The option the name is read from, e.g. `-mizo_cell_type`.
        character(len=*), intent(in) :: key
        !> True for rectangles, false for the triangles the default cuts them into.
        logical :: quadrilateral

        character(len=32) :: name

        call read_text_option(key, name, "triangle")
        select case (trim(name))
        case ("triangle")
            quadrilateral = .false.
        case ("quadrilateral")
            quadrilateral = .true.
        case default
            error stop "A column cell type is either triangle or quadrilateral"
        end select
    end function read_cell_type_option

    !> The water content the wetting data span, m3/m3.
    !>
    !> The normaliser of every water-content comparison below: the whole change
    !> a cell of this problem can undergo between the initial state and the
    !> inlet. Normalising by \(\theta_s\) instead would divide by water that
    !> never moves and make every difference look small.
    pure function wetting_water_span(two_material) result(span)
        implicit none
        logical, intent(in) :: two_material
        real(real64) :: span

        span = reference_storage(WET_INLET_PRESSURE) - reference_storage(WET_START_PRESSURE)
        if (two_material) then
            span = max(span, reference_storage_clay(WET_INLET_PRESSURE) &
                       - reference_storage_clay(WET_START_PRESSURE))
        end if
    end function wetting_water_span

    !> The solved state projected onto `windows` equal windows across the strip.
    !>
    !> \[
    !>   \bar u_j=\frac{1}{|\Omega_j|}\int_{\Omega_j}u\,\mathrm dV,\qquad
    !>   \Omega_j=\left[\frac{(j-1)W}{M},\frac{jW}{M}\right]\times[0,H].
    !> \]
    !>
    !> Assumes a uniformly refined strip whose `nx` is a multiple of `M`. Each
    !> window is then an exact union of whole cells at every level, so for a P0
    !> field the integral is exact and this is the same \(L^2\) projection on
    !> every mesh: two levels can be compared without interpolating either onto
    !> the other. Arithmetic is \(O(N)\) in owned cells, memory \(O(M)\).
    !>
    !> The projection discards the variation inside one window. That part is
    !> measured separately as the transverse spread rather than left unstated;
    !> for this fixture it is the departure from a one-dimensional answer.
    subroutine measure_wetting_profile(solver, windows, split_x, profile)
        implicit none
        type(type_mixed_hydraulic), intent(inout) :: solver
        !> How many windows to project onto; at most `MAX_PROFILE_WINDOWS`.
        integer(int32), intent(in) :: windows
        !> Cells left of this use the van Genuchten curve, the rest Brooks-Corey.
        real(real64), intent(in) :: split_x
        !> The projected state; overwritten on every path.
        type(type_wetting_profile), intent(inout) :: profile

        DM :: dm
        Vec :: local
        PetscInt :: cell, c_start, c_end
        PetscReal :: volume, centroid(3), normal(3)
        real(real64) :: pressure, water, x, area
        real(real64) :: pressure_sum(MAX_PROFILE_WINDOWS), water_sum(MAX_PROFILE_WINDOWS)
        real(real64) :: volume_sum(MAX_PROFILE_WINDOWS)
        real(real64) :: pressure_high(MAX_PROFILE_WINDOWS), pressure_low(MAX_PROFILE_WINDOWS)
        real(real64) :: water_high(MAX_PROFILE_WINDOWS), water_low(MAX_PROFILE_WINDOWS)
        integer(int32) :: window, counted

        profile%windows = 0
        profile%volume = 0.0d0
        profile%pressure = 0.0d0
        profile%water = 0.0d0
        profile%spread_pressure = 0.0d0
        profile%spread_water = 0.0d0
        profile%left_volume = 0.0d0
        profile%cells = 0
        if (windows < 1 .or. windows > MAX_PROFILE_WINDOWS) return
        profile%windows = windows

        pressure_sum = 0.0d0
        water_sum = 0.0d0
        volume_sum = 0.0d0
        ! The minima travel as their own negation, so one maximum reduction
        ! carries both ends of the spread.
        pressure_high = -huge(1.0d0)
        pressure_low = -huge(1.0d0)
        water_high = -huge(1.0d0)
        water_low = -huge(1.0d0)
        counted = 0

        dm = solver%get_dm()
        call solver%get_local_solution(local, ierr)
        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)

        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(dm, cell)) cycle
            call reconstruct_cell_pressure(dm, local, solver%pressure_scale, cell, pressure, ierr)
            call DMPlexComputeCellGeometryFVM(dm, cell, volume, centroid, normal, ierr)
            area = real(volume, real64)
            x = real(centroid(1), real64)
            if (x < split_x) then
                water = reference_storage(pressure)
                profile%left_volume = profile%left_volume + area
            else
                water = reference_storage_clay(pressure)
            end if

            window = 1 + int(x*real(windows, real64)/WET_WIDTH)
            window = max(1, min(windows, window))
            volume_sum(window) = volume_sum(window) + area
            pressure_sum(window) = pressure_sum(window) + pressure*area
            water_sum(window) = water_sum(window) + water*area
            pressure_high(window) = max(pressure_high(window), pressure)
            pressure_low(window) = max(pressure_low(window), -pressure)
            water_high(window) = max(water_high(window), water)
            water_low(window) = max(water_low(window), -water)
            counted = counted + 1
        end do

        call reduce_sum_array(volume_sum, windows)
        call reduce_sum_array(pressure_sum, windows)
        call reduce_sum_array(water_sum, windows)
        call reduce_max_array(pressure_high, windows)
        call reduce_max_array(pressure_low, windows)
        call reduce_max_array(water_high, windows)
        call reduce_max_array(water_low, windows)
        call reduce_sum(profile%left_volume)
        call reduce_count(counted)
        profile%cells = counted

        do window = 1, windows
            profile%volume(window) = volume_sum(window)
            if (volume_sum(window) <= 0.0d0) cycle
            profile%pressure(window) = pressure_sum(window)/volume_sum(window)
            profile%water(window) = water_sum(window)/volume_sum(window)
            profile%spread_pressure = max(profile%spread_pressure, &
                                              pressure_high(window) + pressure_low(window))
            profile%spread_water = max(profile%spread_water, &
                                           water_high(window) + water_low(window))
        end do
    end subroutine measure_wetting_profile

    !> The volume-weighted RMS and the window maximum of two profiles' gap.
    !>
    !> \[
    !>   \|u-v\|_{\mathrm{rms}}=\sqrt{\frac{\sum_j|\Omega_j|(\bar u_j-\bar v_j)^2}{\sum_j|\Omega_j|}},
    !>   \qquad \|u-v\|_\infty=\max_j|\bar u_j-\bar v_j|.
    !> \]
    !>
    !> The first is the \(L^2\) norm of the difference of the two projections
    !> divided by the square root of the domain measure, retaining the quantity's units.
    subroutine compare_profiles(left, right, pressure_rms, pressure_max, water_rms, water_max)
        implicit none
        type(type_wetting_profile), intent(in) :: left, right
        real(real64), intent(inout) :: pressure_rms, pressure_max
        real(real64), intent(inout) :: water_rms, water_max

        real(real64) :: total, gap
        integer(int32) :: window

        pressure_rms = 0.0d0
        pressure_max = 0.0d0
        water_rms = 0.0d0
        water_max = 0.0d0
        if (left%windows /= right%windows .or. left%windows < 1) return

        total = sum(left%volume(1:left%windows))
        if (total <= 0.0d0) return

        do window = 1, left%windows
            gap = left%pressure(window) - right%pressure(window)
            pressure_rms = pressure_rms + left%volume(window)*gap*gap
            pressure_max = max(pressure_max, abs(gap))
            gap = left%water(window) - right%water(window)
            water_rms = water_rms + left%volume(window)*gap*gap
            water_max = max(water_max, abs(gap))
        end do
        pressure_rms = sqrt(pressure_rms/total)
        water_rms = sqrt(water_rms/total)
    end subroutine compare_profiles

    !> Run one level of a refinement ladder and record everything it measured.
    subroutine run_refinement_level(materials, split_x, grid, windows, level)
        implicit none
        type(type_mixed_material), intent(in) :: materials(:)
        !> Where the material changes, or a number beyond the strip if it does not.
        real(real64), intent(in) :: split_x
        type(type_wetting_grid), intent(in) :: grid
        integer(int32), intent(in) :: windows
        type(type_refinement_level), intent(inout) :: level

        type(type_mixed_hydraulic) :: solver
        DM :: dm
        real(real64) :: start_inventory, inventory

        level%grid = grid
        level%run = type_wetting_result()
        level%storage_change = 0.0d0
        level%budget_error = 0.0d0

        call start_wetting_run(solver, dm, materials, grid, ierr)
        if (ierr /= 0) return

        start_inventory = measure_reference_inventory(solver, split_x)
        call run_wetting_interval(solver, grid, REFINE_END_TIME, level%run)
        if (level%run%ok) then
            inventory = measure_reference_inventory(solver, split_x)
            level%storage_change = inventory - start_inventory
            if (abs(level%run%influx) > 0.0d0) then
                level%budget_error = (level%storage_change - level%run%influx)/abs(level%run%influx)
            end if
            call measure_wetting_profile(solver, windows, split_x, level%profile)
        end if

        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)
    end subroutine run_refinement_level

    !> Print one level's row and check what has to hold at that level alone.
    subroutine check_level(tag, level, expected_steps, failures)
        implicit none
        character(len=*), intent(in) :: tag
        type(type_refinement_level), intent(in) :: level
        !> How many steps of the nominal length reach the end time.
        integer(int32), intent(in) :: expected_steps
        integer(int32), intent(inout) :: failures

        character(len=200) :: line
        real(real64) :: transported, step_spread

        write (line, '(a,i4,a,i4,a,f9.2,a,i4,a,i3,a,i5,a,f8.2,a,es11.4,a,es10.3)') &
            "    nx ", level%grid%nx, "  nz ", level%grid%nz, "  dt ", level%grid%dt, &
            "  steps ", level%run%steps, "  refused ", level%run%rejections, &
            "  newton ", level%run%newton, "  wall/s ", level%run%wall_time, &
            "  influx/m2 ", level%run%influx, "  budget ", level%budget_error
        call report_line(trim(line))
        write (line, '(a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,es10.3)') &
            "      cells/inv ", level%run%worst_inventory, "  cells/thr ", level%run%worst_throughput, &
            "  overshoot ", level%run%overshoot, "  undershoot ", level%run%undershoot, &
            "  residual ", level%run%max_residual, "  spread_p/Pa ", level%profile%spread_pressure
        call report_line(trim(line))

        call check_integer(tag//" reached the end time     ", merge(1, 0, level%run%ok), 1, failures)
        if (.not. level%run%ok) return

        ! The nominal step has to be the step that ran, or the ladder is not a
        ! ladder in the step length at all.
        call check_integer(tag//" took its nominal steps   ", level%run%steps, expected_steps, failures)
        call check_integer(tag//" refused no step          ", level%run%rejections, 0, failures)
        step_spread = 0.0d0
        if (level%run%shortest_step > 0.0d0) then
            step_spread = level%run%longest_step/level%run%shortest_step - 1.0d0
        end if
        call check(tag//" ran one step length     ", step_spread, 1.0d-9, failures)

        ! The physics gates of the earlier cases, at every level of the ladder.
        call check(tag//" inventory against budget ", level%budget_error, 1.0d-10, failures)
        call check(tag//" local conservation/inv   ", level%run%worst_inventory, 1.0d-10, failures)
        call check(tag//" local conservation/thr   ", level%run%worst_throughput, 1.0d-8, failures)
        call check(tag//" accepted residual norm   ", level%run%max_residual, 1.0d-9, failures)

        ! Overshoot is reported and bounded separately from the mass budget:
        ! RT0/P0 carries no discrete maximum principle, so this is a measured
        ! limit on how far outside its own data the answer goes, not a proof
        ! that it cannot happen.
        call check(tag//" pressure overshoot       ", level%run%overshoot, 5.0d-2, failures)
        call check(tag//" pressure undershoot      ", level%run%undershoot, 5.0d-2, failures)

        transported = level%run%influx/(WET_WIDTH*WET_HEIGHT*SOIL_THETA_S)
        call check_exceeds(tag//" water actually moved     ", transported, 1.0d-2, failures)
    end subroutine check_level

    !> Print a ladder's differences and check that it converges.
    !>
    !> Two readings of the same data. The successive differences make no use of
    !> the finest level as an answer; the errors against it do, and are reported
    !> with the reference's own uncertainty beside them. Neither claims the
    !> finest level is exact, and no universal order is required: a wetting front
    !> is not smooth on the scale of these meshes, so the observed rate is
    !> reported and the assertion is only that the sequence keeps contracting.
    subroutine report_refinement_study(title, tag, levels, count, water_span, contraction, &
                                       assert_convergence, failures)
        implicit none
        character(len=*), intent(in) :: title
        !> A short name the assertions are labelled with.
        character(len=*), intent(in) :: tag
        type(type_refinement_level), intent(in) :: levels(:)
        integer(int32), intent(in) :: count
        !> The water content the data span, m3/m3.
        real(real64), intent(in) :: water_span
        !> How much each successive difference has to shrink.
        real(real64), intent(in) :: contraction
        !> False for a ladder whose table is reported without a pass or fail.
        !>
        !> A ladder that is known not to have converged is still worth measuring
        !> and printing; asserting a convergence criterion on it would either
        !> fail every run or have to be weakened until it said nothing.
        logical, intent(in) :: assert_convergence
        integer(int32), intent(inout) :: failures

        real(real64) :: diff_p_rms(MAX_REFINEMENT_LEVELS), diff_p_max(MAX_REFINEMENT_LEVELS)
        real(real64) :: diff_w_rms(MAX_REFINEMENT_LEVELS), diff_w_max(MAX_REFINEMENT_LEVELS)
        real(real64) :: err_p_rms(MAX_REFINEMENT_LEVELS), err_p_max(MAX_REFINEMENT_LEVELS)
        real(real64) :: err_w_rms(MAX_REFINEMENT_LEVELS), err_w_max(MAX_REFINEMENT_LEVELS)
        real(real64) :: pressure_span, rate_p, rate_w, margin
        character(len=200) :: line
        integer(int32) :: l

        pressure_span = abs(WET_INLET_PRESSURE - WET_START_PRESSURE)
        diff_p_rms = 0.0d0
        diff_p_max = 0.0d0
        diff_w_rms = 0.0d0
        diff_w_max = 0.0d0
        err_p_rms = 0.0d0
        err_p_max = 0.0d0
        err_w_rms = 0.0d0
        err_w_max = 0.0d0

        do l = 1, count - 1
            call compare_profiles(levels(l)%profile, levels(l + 1)%profile, &
                                  diff_p_rms(l), diff_p_max(l), diff_w_rms(l), diff_w_max(l))
            diff_p_rms(l) = diff_p_rms(l)/pressure_span
            diff_p_max(l) = diff_p_max(l)/pressure_span
            diff_w_rms(l) = diff_w_rms(l)/water_span
            diff_w_max(l) = diff_w_max(l)/water_span

            call compare_profiles(levels(l)%profile, levels(count)%profile, &
                                  err_p_rms(l), err_p_max(l), err_w_rms(l), err_w_max(l))
            err_p_rms(l) = err_p_rms(l)/pressure_span
            err_p_max(l) = err_p_max(l)/pressure_span
            err_w_rms(l) = err_w_rms(l)/water_span
            err_w_max(l) = err_w_max(l)/water_span
        end do

        call report_line("  "//title//": successive differences of the window profiles")
        call report_line("    l    p_rms      p_max      w_rms      w_max      rate_p  rate_w")
        do l = 1, count - 1
            rate_p = 0.0d0
            rate_w = 0.0d0
            if (l < count - 1) then
                if (diff_p_rms(l + 1) > 0.0d0) rate_p = log(diff_p_rms(l)/diff_p_rms(l + 1))/log(2.0d0)
                if (diff_w_rms(l + 1) > 0.0d0) rate_w = log(diff_w_rms(l)/diff_w_rms(l + 1))/log(2.0d0)
            end if
            write (line, '(a,i2,4(a,es10.3),2(a,f7.3))') "    ", l, "  ", diff_p_rms(l), "  ", diff_p_max(l), &
                "  ", diff_w_rms(l), "  ", diff_w_max(l), "   ", rate_p, " ", rate_w
            call report_line(trim(line))
        end do

        call report_line("  "//title//": against the finest level, which is a reference and not an answer")
        do l = 1, count - 1
            write (line, '(a,i2,4(a,es10.3))') "    ", l, "  ", err_p_rms(l), "  ", err_p_max(l), &
                "  ", err_w_rms(l), "  ", err_w_max(l)
            call report_line(trim(line))
        end do

        margin = 0.0d0
        if (diff_p_rms(count - 1) > 0.0d0) margin = err_p_rms(1)/diff_p_rms(count - 1)

        if (.not. assert_convergence) then
            call report_line("    this ladder is reported, not asserted; see the plan checkpoint for why")
            do l = 2, count - 1
                write (line, '(a,i2,a,f9.4,a,f9.4)') "    ratio at ", l, "  pressure ", &
                    diff_p_rms(l)/max(diff_p_rms(l - 1), tiny(1.0d0)), "  water ", &
                    diff_w_rms(l)/max(diff_w_rms(l - 1), tiny(1.0d0))
                call report_line(trim(line))
            end do
            write (line, '(a,f9.4)') "    reference margin ", margin
            call report_line(trim(line))
            return
        end if

        ! The ladder has to be resolving something. A run that had already
        ! reached a state the mesh or the step could not change would contract
        ! perfectly and prove nothing.
        call check_exceeds(tag//" coarsest p difference    ", diff_p_rms(1), REFINE_SIGNAL_FLOOR, failures)
        call check_exceeds(tag//" coarsest w difference    ", diff_w_rms(1), REFINE_SIGNAL_FLOOR, failures)

        do l = 2, count - 1
            if (diff_p_rms(l - 1) > 0.0d0) then
                call check(tag//" p difference contracts   ", diff_p_rms(l)/diff_p_rms(l - 1), &
                           contraction, failures)
            end if
            if (diff_w_rms(l - 1) > 0.0d0) then
                call check(tag//" w difference contracts   ", diff_w_rms(l)/diff_w_rms(l - 1), &
                           contraction, failures)
            end if
            if (err_p_rms(l - 1) > 0.0d0) then
                call check(tag//" p error contracts        ", err_p_rms(l)/err_p_rms(l - 1), &
                           contraction, failures)
            end if
            if (err_w_rms(l - 1) > 0.0d0) then
                call check(tag//" w error contracts        ", err_w_rms(l)/err_w_rms(l - 1), &
                           contraction, failures)
            end if
        end do

        ! Whether the finest level is a usable reference for the levels it
        ! grades: its own remaining difference has to be well below the error it
        ! is being used to measure. This is an estimate of that uncertainty, not
        ! a claim that the level is converged.
        call check_exceeds(tag//" reference is adequate    ", margin, REFINE_REFERENCE_MARGIN, failures)
    end subroutine report_refinement_study

    !> The water the second material gained, m2 per unit thickness.
    !>
    !> \[
    !>   U=\sum_{j\ge j_0}|\Omega_j|\,(\bar W_j-W(p_0)),
    !> \]
    !>
    !> with \(j_0\) the first window beyond the material split. What the
    !> low-conductivity soil took up is the quantity the two-material mesh ladder
    !> fails to converge in, so it is reported on its own rather than left inside
    !> a norm over the whole strip.
    pure function second_material_uptake(profile, split_x) result(uptake)
        implicit none
        type(type_wetting_profile), intent(in) :: profile
        real(real64), intent(in) :: split_x
        real(real64) :: uptake

        real(real64) :: initial, centre
        integer(int32) :: window

        uptake = 0.0d0
        initial = reference_storage_clay(WET_START_PRESSURE)
        do window = 1, profile%windows
            centre = WET_WIDTH*(real(window, real64) - 0.5d0)/real(profile%windows, real64)
            if (centre < split_x) cycle
            uptake = uptake + profile%volume(window)*(profile%water(window) - initial)
        end do
    end function second_material_uptake

    !> Report how the second material's uptake behaves along a ladder.
    subroutine report_uptake_ladder(levels, count, split_x)
        implicit none
        type(type_refinement_level), intent(in) :: levels(:)
        integer(int32), intent(in) :: count
        real(real64), intent(in) :: split_x

        character(len=160) :: line
        real(real64) :: uptake, previous
        integer(int32) :: l

        call report_line("    uptake of the second material, and its ratio between levels")
        previous = 0.0d0
        do l = 1, count
            uptake = second_material_uptake(levels(l)%profile, split_x)
            if (l == 1) then
                write (line, '(a,i2,a,es12.5)') "      level ", l, "  uptake/m2 ", uptake
            else
                write (line, '(a,i2,a,es12.5,a,f9.3)') "      level ", l, "  uptake/m2 ", uptake, &
                    "  ratio ", uptake/max(abs(previous), tiny(1.0d0))
            end if
            call report_line(trim(line))
            previous = uptake
        end do
    end subroutine report_uptake_ladder

    !> The mesh ladder: the step, the end time, the data and the materials fixed.
    subroutine study_space(failures, two_material)
        implicit none
        integer(int32), intent(inout) :: failures
        !> True for the van Genuchten / Brooks-Corey pair either side of the split.
        logical, intent(in) :: two_material

        type(type_refinement_level) :: levels(MAX_REFINEMENT_LEVELS)
        type(type_mixed_material) :: one_material(1), pair(2)
        type(type_wetting_grid) :: grid
        real(real64) :: split_x
        character(len=24) :: tag
        character(len=64) :: title
        character(len=120) :: line
        integer(int32) :: l, steps, count

        call set_soil_material(one_material(1), 1, MIXED_CONDUCTIVITY_HCF)
        call set_soil_material(pair(1), 3, MIXED_CONDUCTIVITY_HCF)
        call set_clay_material(pair(2), 7, MIXED_CONDUCTIVITY_HCF)
        split_x = 1.0d30
        if (two_material) split_x = WET_SPLIT

        ! Test-only knobs: the registered case runs the defaults, and a longer
        ! validation run extends the ladder without a rebuild.
        steps = max(1, read_integer_option("-mixed_refine_space_steps", SPACE_STEPS))
        count = max(2, min(MAX_SPACE_LEVELS, read_integer_option("-mixed_refine_space_levels", SPACE_LEVELS)))

        if (two_material) then
            tag = "space2"
            title = "two-material mesh ladder"
        else
            tag = "space1"
            title = "homogeneous mesh ladder"
        end if
        call report_line("  "//trim(title)//": same step, end time, data and materials at every level")
        write (line, '(a,i0,a,f8.2,a,i0)') "    levels ", count, ", step/s ", &
            REFINE_END_TIME/real(steps, real64), ", steps per level ", steps
        call report_line(trim(line))

        do l = 1, count
            grid = default_wetting_grid()
            grid%nx = SPACE_NX(l)
            grid%nz = SPACE_NZ(l)
            grid%dt = REFINE_END_TIME/real(steps, real64)
            grid%fixed_step = .true.

            ! The windows have to be unions of whole cells, and the material
            ! interface has to be on a cell boundary, at every level.
            call check_integer(trim(tag)//" windows divide the mesh  ", &
                               mod(grid%nx, SPACE_WINDOWS), 0, failures)
            call check_integer(trim(tag)//" split is on a boundary   ", &
                               mod(grid%nx, 2), 0, failures)

            if (two_material) then
                call run_refinement_level(pair, split_x, grid, SPACE_WINDOWS, levels(l))
            else
                call run_refinement_level(one_material, split_x, grid, SPACE_WINDOWS, levels(l))
            end if
            call check_level(trim(tag), levels(l), steps, failures)
            call check_integer(trim(tag)//" cells of the whole mesh  ", levels(l)%profile%cells, &
                               2*grid%nx*grid%nz, failures)
            if (two_material .and. levels(l)%run%ok) then
                ! Measured, not assumed: the cells of the first material fill
                ! exactly the part of the strip that material occupies, so the
                ! interface is mesh-aligned at this level.
                call check(trim(tag)//" interface is aligned     ", &
                           (levels(l)%profile%left_volume - WET_SPLIT*WET_HEIGHT)/(WET_SPLIT*WET_HEIGHT), &
                           1.0d-14, failures)
            end if
            if (.not. levels(l)%run%ok) return
        end do

        ! The two-material mesh ladder is measured and printed but not asserted:
        ! the wetting front inside the Brooks-Corey soil is a boundary layer of
        ! the order of a millimetre after this interval, which is the size of the
        ! finest cell here, so none of these meshes resolves it and the ladder is
        ! not in a convergent regime. The uptake below is where that shows.
        call report_refinement_study(trim(title), trim(tag), levels, count, &
                                     wetting_water_span(two_material), REFINE_CONTRACTION_SPACE, &
                                     .not. two_material, failures)
        if (two_material) call report_uptake_ladder(levels, count, split_x)
        call write_refinement_csv(trim(tag), levels, count)
    end subroutine study_space

    !> The step ladder: one mesh, the end time, the data and the materials fixed.
    subroutine study_time(failures, two_material)
        implicit none
        integer(int32), intent(inout) :: failures
        logical, intent(in) :: two_material

        type(type_refinement_level) :: levels(MAX_REFINEMENT_LEVELS)
        type(type_mixed_material) :: one_material(1), pair(2)
        type(type_wetting_grid) :: grid
        real(real64) :: split_x
        character(len=24) :: tag
        character(len=64) :: title
        character(len=120) :: line
        integer(int32) :: l, count, nx, nz

        call set_soil_material(one_material(1), 1, MIXED_CONDUCTIVITY_HCF)
        call set_soil_material(pair(1), 3, MIXED_CONDUCTIVITY_HCF)
        call set_clay_material(pair(2), 7, MIXED_CONDUCTIVITY_HCF)
        split_x = 1.0d30
        if (two_material) split_x = WET_SPLIT

        count = max(2, min(TIME_LEVELS, read_integer_option("-mixed_refine_time_levels", TIME_LEVELS)))
        nx = TIME_NX
        nz = TIME_NZ

        if (two_material) then
            tag = "time2"
            title = "two-material step ladder"
        else
            tag = "time1"
            title = "homogeneous step ladder"
        end if
        call report_line("  "//trim(title)//": one mesh, same end time, data and materials")
        write (line, '(a,i0,a,i0,a,i0)') "    levels ", count, ", mesh nx ", nx, " nz ", nz
        call report_line(trim(line))

        do l = 1, count
            grid = default_wetting_grid()
            grid%nx = nx
            grid%nz = nz
            grid%dt = REFINE_END_TIME/real(TIME_STEPS(l), real64)
            grid%fixed_step = .true.

            if (two_material) then
                call run_refinement_level(pair, split_x, grid, TIME_WINDOWS, levels(l))
            else
                call run_refinement_level(one_material, split_x, grid, TIME_WINDOWS, levels(l))
            end if
            call check_level(trim(tag), levels(l), TIME_STEPS(l), failures)
            call check_integer(trim(tag)//" kept the same mesh       ", levels(l)%profile%cells, &
                               2*nx*nz, failures)
            if (.not. levels(l)%run%ok) return
        end do

        call report_refinement_study(trim(title), trim(tag), levels, count, &
                                     wetting_water_span(two_material), REFINE_CONTRACTION_TIME, &
                                     .true., failures)
        if (two_material) call report_uptake_ladder(levels, count, split_x)
        call write_refinement_csv(trim(tag), levels, count)
    end subroutine study_time

    !> What the window projection has to satisfy however the mesh is partitioned.
    !>
    !> The one case registered on two ranks. The projection and its reductions
    !> are the only new sampling logic of this increment, so what is checked is
    !> that they answer the same invariants whatever rank owns which cell: the
    !> whole mesh is counted once, the windows tile the domain exactly, the water
    !> they integrate is the water an independent per-cell sum finds, and the
    !> profile of a front entering from the left still decreases with x.
    subroutine check_profile_sampling(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_hydraulic) :: solver
        type(type_mixed_material) :: pair(2)
        type(type_wetting_profile) :: profile
        type(type_wetting_grid) :: grid
        DM :: dm
        character(len=200) :: line
        real(real64) :: inventory, projected, worst_rise
        integer(int32) :: window
        logical :: accepted

        call report_line("  window projection: the invariants its reductions have to satisfy")

        call set_soil_material(pair(1), 3, MIXED_CONDUCTIVITY_HCF)
        call set_clay_material(pair(2), 7, MIXED_CONDUCTIVITY_HCF)

        grid = default_wetting_grid()
        grid%nx = SPACE_NX(1)
        grid%nz = SPACE_NZ(1)
        grid%dt = REFINE_END_TIME/real(SPACE_STEPS, real64)
        grid%fixed_step = .true.

        call start_wetting_run(solver, dm, pair, grid, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the sampling fixture could not be built")
            failures = failures + 1
            return
        end if

        accepted = .false.
        call solver%solve_step(grid%dt, grid%dt, accepted, ierr)
        call check_integer("sampling step converged                 ", &
                           merge(1, 0, accepted .and. ierr == 0), 1, failures)
        if (.not. accepted) then
            call solver%destroy(ierr)
            call DMDestroy(dm, ierr)
            return
        end if

        call measure_wetting_profile(solver, SPACE_WINDOWS, WET_SPLIT, profile)
        inventory = measure_reference_inventory(solver, WET_SPLIT)

        projected = 0.0d0
        worst_rise = 0.0d0
        do window = 1, profile%windows
            projected = projected + profile%volume(window)*profile%water(window)
            if (window > 1) then
                worst_rise = max(worst_rise, profile%pressure(window) - profile%pressure(window - 1))
            end if
            write (line, '(a,i3,3(a,es22.15))') "    window ", window, "  volume/m2 ", profile%volume(window), &
                "  pressure/Pa ", profile%pressure(window), "  water ", profile%water(window)
            call report_line(trim(line))
        end do

        call check_integer("every cell counted exactly once         ", profile%cells, &
                           2*grid%nx*grid%nz, failures)
        call check("the windows tile the domain             ", &
                   (sum(profile%volume(1:profile%windows)) - WET_WIDTH*WET_HEIGHT)/(WET_WIDTH*WET_HEIGHT), &
                   1.0d-14, failures)
        call check("the first material fills its part       ", &
                   (profile%left_volume - WET_SPLIT*WET_HEIGHT)/(WET_SPLIT*WET_HEIGHT), 1.0d-14, failures)
        ! The projected water and the per-cell sum are two independent walks
        ! over the same owned cells; a ghost counted twice breaks only one.
        call check("projected water against the inventory   ", &
                   (projected - inventory)/abs(inventory), 1.0d-12, failures)
        call check("the profile decreases away from the inlet", &
                   max(0.0d0, worst_rise)/abs(WET_INLET_PRESSURE - WET_START_PRESSURE), 1.0d-12, failures)
        call check_exceeds("the front is inside the strip           ", &
                           (profile%pressure(1) - profile%pressure(profile%windows)) &
                           /abs(WET_INLET_PRESSURE - WET_START_PRESSURE), 1.0d-2, failures)

        call solver%destroy(ierr)
        call DMDestroy(dm, ierr)
    end subroutine check_profile_sampling

    !> Append one ladder's measured table to the CSV artifact, if one was asked for.
    subroutine write_refinement_csv(tag, levels, count)
        implicit none
        character(len=*), intent(in) :: tag
        type(type_refinement_level), intent(in) :: levels(:)
        integer(int32), intent(in) :: count

        PetscMPIInt :: rank
        PetscErrorCode :: local_ierr
        integer(int32) :: unit_number, l, window, open_status
        real(real64) :: centre

        if (len_trim(csv_prefix) == 0) return
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, local_ierr)
        if (rank /= 0) return

        ! An unwritable prefix costs the artifact, not the measurements that
        ! were already made and already printed.
        if (.not. csv_started) then
            open (newunit=unit_number, file=trim(csv_prefix)//"-levels.csv", status="replace", &
                  action="write", iostat=open_status)
            if (open_status /= 0) then
                call report_line("    the CSV artifact could not be opened; the tables above are the record")
                csv_prefix = ""
                return
            end if
            write (unit_number, '(a)') "study,level,nx,nz,cells,dt_s,steps,rejections,newton,"// &
                "shortest_dt_s,longest_dt_s,wall_s,influx_m2,storage_change_m2,budget_rel,"// &
                "defect_over_inventory,defect_over_throughput,overshoot_rel,undershoot_rel,"// &
                "residual,spread_pressure_pa,spread_water"
            close (unit_number)
            open (newunit=unit_number, file=trim(csv_prefix)//"-profiles.csv", status="replace", &
                  action="write", iostat=open_status)
            if (open_status /= 0) then
                call report_line("    the CSV artifact could not be opened; the tables above are the record")
                csv_prefix = ""
                return
            end if
            write (unit_number, '(a)') "study,level,nx,nz,dt_s,window,x_centre_m,window_volume_m2,"// &
                "pressure_pa,water_content"
            close (unit_number)
            csv_started = .true.
        end if

        open (newunit=unit_number, file=trim(csv_prefix)//"-levels.csv", status="old", position="append", &
              action="write")
        do l = 1, count
            write (unit_number, '(a,",",i0,",",i0,",",i0,",",i0,",",es22.15,",",i0,",",i0,",",i0,13(",",es22.15))') &
                tag, l, levels(l)%grid%nx, levels(l)%grid%nz, levels(l)%profile%cells, levels(l)%grid%dt, &
                levels(l)%run%steps, levels(l)%run%rejections, levels(l)%run%newton, &
                levels(l)%run%shortest_step, levels(l)%run%longest_step, levels(l)%run%wall_time, &
                levels(l)%run%influx, levels(l)%storage_change, levels(l)%budget_error, &
                levels(l)%run%worst_inventory, levels(l)%run%worst_throughput, &
                levels(l)%run%overshoot, levels(l)%run%undershoot, levels(l)%run%max_residual, &
                levels(l)%profile%spread_pressure, levels(l)%profile%spread_water
        end do
        close (unit_number)

        open (newunit=unit_number, file=trim(csv_prefix)//"-profiles.csv", status="old", position="append", &
              action="write")
        do l = 1, count
            do window = 1, levels(l)%profile%windows
                centre = WET_WIDTH*(real(window, real64) - 0.5d0)/real(levels(l)%profile%windows, real64)
                write (unit_number, '(a,",",i0,",",i0,",",i0,",",es22.15,",",i0,4(",",es22.15))') &
                    tag, l, levels(l)%grid%nx, levels(l)%grid%nz, levels(l)%grid%dt, window, centre, &
                    levels(l)%profile%volume(window), levels(l)%profile%pressure(window), &
                    levels(l)%profile%water(window)
            end do
        end do
        close (unit_number)
    end subroutine write_refinement_csv

    !> Mesh and step refinement of the wetting front, in separate ladders.
    subroutine case_refinement(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        character(len=32) :: scope

        call report_line("case refinement: mesh and step ladders of the isothermal wetting front")
        call read_text_option("-mixed_refine_scope", scope, "full")
        call read_text_option("-mixed_refine_csv", csv_prefix, "")

        ! The undamped Newton step, because the measured alternative is worse.
        !
        ! PETSc's default `bt` backtracking rejects the full Newton step on this
        ! problem and settles on lambda near 0.02, which reduces the residual by
        ! about half a percent per iteration and runs out of iterations from
        ! `nx = 32` upwards -- so the ladder could not be built with it at all.
        ! The direction is not the cause: applying the true Jacobian matrix-free
        ! with `-snes_mf_operator` reproduces the same iteration counts, and with
        ! the full step every level converges in about five iterations per step
        ! at every mesh of this ladder. The measured numbers and the `l2`
        ! alternative are in the plan's checkpoint. This is the study's own
        ! solver configuration; nothing in the model, the closure or any
        ! tolerance is changed, and a step that does fail still goes through the
        ! solver's own rejection path.
        call set_default_option("-snes_linesearch_type", "basic")

        select case (trim(scope))
        case ("space")
            call study_space(failures, .false.)
        case ("time")
            call study_time(failures, .false.)
        case ("interface")
            call study_space(failures, .true.)
            call study_time(failures, .true.)
        case ("sample")
            call check_profile_sampling(failures)
        case ("full")
            call check_profile_sampling(failures)
            call study_space(failures, .false.)
            call study_time(failures, .false.)
            call study_space(failures, .true.)
            call study_time(failures, .true.)
        case default
            call report_line("  FAIL unknown refinement scope: "//trim(scope))
            failures = failures + 1
        end select
    end subroutine case_refinement

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

    !> Sum the first `count` entries of an array over the communicator.
    subroutine reduce_sum_array(values, count)
        implicit none
        real(real64), intent(inout) :: values(:)
        integer(int32), intent(in) :: count

        real(real64) :: reduced(size(values))
        integer(int32) :: mpi_err

        if (count < 1) return
        call MPI_Allreduce(values, reduced, count, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, mpi_err)
        values(1:count) = reduced(1:count)
    end subroutine reduce_sum_array

    !> The frozen constitutive contract: \((T,p_w,\text{material})\to(W,U,p_l,k_r)\).
    !>
    !> Four states of the Mizoguchi soil are probed -- unfrozen, freezing,
    !> saturated and pore-limited -- and every number is measured against
    !> something that is not the contract's own arithmetic: the already-verified
    !> isothermal closure in the unfrozen limit, a recombination of the reported
    !> phase contents and densities, central differences of the contract itself
    !> away from its kinks, and the repository's conductivity models evaluated
    !> twice on the same published state.
    subroutine case_freeze_closure(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_frozen_material) :: frozen(1)
        type(type_mixed_material) :: isothermal(1)

        call report_line("case freeze_closure: the frozen constitutive contract of the Mizoguchi soil")

        call set_frozen_soil(frozen(1), FROZEN_MATERIAL_ID)
        ! The project the parameters come from runs with vapour transport off.
        call install_frozen_laws(frozen, COMP_TYPES%XZ_2D%ID, .false.)

        ! Every check below states an expected value, so running them against a
        ! memo is a direct test that a recalled value is the computed one:
        ! -mizo_tangent_cache 1 has to leave every one of them passing.
        call configure_frozen_tangent_cache( &
            read_integer_option("-mizo_tangent_cache", 0) /= 0, &
            read_integer_option("-mizo_tangent_cache_entries", 0))

        call check_frozen_states(failures)
        call check_frozen_domain(failures)
        call check_frozen_unfrozen_limit(isothermal, failures)
        call check_frozen_tangents(failures)
        call check_frozen_branch(failures)
        call check_frozen_rank(failures)
        call check_frozen_impedance(frozen(1), failures)
        call check_frozen_history(failures)

        call clear_frozen_laws()
    end subroutine case_freeze_closure

    !> The Mizoguchi soil of `project/Mizo-xz-Convex/Input/Basic.json`.
    !>
    !> Written out here rather than read from that project, so the case needs no
    !> project directory; the values are the ones the archived runs used.
    subroutine set_frozen_soil(material, id)
        implicit none
        type(type_mixed_frozen_material), intent(inout) :: material
        !> The region identifier the governing models key this material by.
        integer(int32), intent(in) :: id

        material%id = id
        material%porosity = SOIL_THETA_S
        material%specific_storage = FROZEN_SPECIFIC_STORAGE

        call material%density%reset()
        material%density%material_id = id
        material%density%num_phases = FROZEN_PHASES
        material%density%solid = 2800.0d0
        material%density%water = 1000.0d0
        material%density%ice = 917.0d0
        material%density%vapor = 0.6d0

        call material%specific_heat%reset()
        material%specific_heat%material_id = id
        material%specific_heat%num_phases = FROZEN_PHASES
        material%specific_heat%solid = 912.0d0
        material%specific_heat%water = 4180.0d0
        material%specific_heat%ice = 2100.0d0
        material%specific_heat%vapor = 2098.0d0

        call material%heat_capacity%reset()
        material%heat_capacity%material_id = id
        material%heat_capacity%num_phases = FROZEN_PHASES
        material%heat_capacity%solid = 2553600.0d0
        material%heat_capacity%water = 4180000.0d0
        material%heat_capacity%ice = 1925700.0d0
        material%heat_capacity%vapor = 1259000.0d0
        material%heat_capacity%params = [1.0d0, 1.0d0, 1.0d0]

        call material%thermal_conductivity%reset()
        material%thermal_conductivity%material_id = id
        material%thermal_conductivity%num_phases = FROZEN_PHASES
        material%thermal_conductivity%solid = 3.78d0
        material%thermal_conductivity%water = 0.6d0
        material%thermal_conductivity%ice = 2.2d0
        material%thermal_conductivity%vapor = 0.0241d0
        material%thermal_conductivity%params = [0.55d0, 0.8d0, 3.07d0, 0.13d0, 13.05d0, 1.06d0]

        call material%retention%reset()
        material%retention%material_id = id
        material%retention%swcc_model = SWCC_MODELS%VG
        material%retention%theta_s = SOIL_THETA_S
        material%retention%theta_r = SOIL_THETA_R
        material%retention%alpha1 = SOIL_ALPHA
        material%retention%n1 = SOIL_N
        material%retention%m1 = SOIL_M

        call material%conductivity_function%reset()
        material%conductivity_function%material_id = id
        ! Model 4 of the project: the retention-driven curve times the ice
        ! impedance, which is the law the impedance factor was calibrated in.
        material%conductivity_function%model = HCF_MODES%BASE_IMPEDANCE
        material%conductivity_function%swcc_model = SWCC_MODELS%VG
        material%conductivity_function%theta_s = SOIL_THETA_S
        material%conductivity_function%theta_r = SOIL_THETA_R
        material%conductivity_function%alpha1 = SOIL_ALPHA
        material%conductivity_function%n1 = SOIL_N
        material%conductivity_function%m1 = SOIL_M
        material%conductivity_function%l = SOIL_L
        material%conductivity_function%k_sat = FROZEN_CONDUCTIVITY
        material%conductivity_function%omega = FROZEN_IMPEDANCE_FACTOR
        material%conductivity_function%water_viscosity_model = HCF_VISCOSITY_TYPES%to_object(1)
        material%conductivity_function%gain_factor = 1.0d0

        call material%freezing%reset()
        material%freezing%material_id = id
        material%freezing%gcc_model = GCC_TYPES%NON_SEGREGATION
    end subroutine set_frozen_soil

    !> The same soil as an isothermal mixed material, for the unfrozen limit.
    !>
    !> The conductivity mode is the base one, which is the only mode the
    !> isothermal closure evaluates: at an unfrozen state the impedance factor of
    !> the frozen law is one, so the two laws are the same function there and the
    !> comparison is between the paths rather than between two models.
    subroutine set_frozen_soil_isothermal(material, id)
        implicit none
        type(type_mixed_material), intent(inout) :: material
        integer(int32), intent(in) :: id

        material%id = id
        material%conductivity = reshape([FROZEN_CONDUCTIVITY, 0.0d0, 0.0d0, FROZEN_CONDUCTIVITY], [2, 2])
        material%source = 0.0d0

        material%storage_law = MIXED_STORAGE_RETENTION
        call material%retention%reset()
        material%retention%swcc_model = SWCC_MODELS%VG
        material%retention%theta_s = SOIL_THETA_S
        material%retention%theta_r = SOIL_THETA_R
        material%retention%alpha1 = SOIL_ALPHA
        material%retention%n1 = SOIL_N
        material%retention%m1 = SOIL_M
        material%specific_storage = FROZEN_SPECIFIC_STORAGE
        material%saturation_pressure = 0.0d0

        material%conductivity_law = MIXED_CONDUCTIVITY_HCF
        call material%conductivity_function%reset()
        material%conductivity_function%model = HCF_MODES%BASE
        material%conductivity_function%swcc_model = SWCC_MODELS%VG
        material%conductivity_function%theta_s = SOIL_THETA_S
        material%conductivity_function%theta_r = SOIL_THETA_R
        material%conductivity_function%alpha1 = SOIL_ALPHA
        material%conductivity_function%n1 = SOIL_N
        material%conductivity_function%m1 = SOIL_M
        material%conductivity_function%k_sat = FROZEN_CONDUCTIVITY
        material%conductivity_function%l = SOIL_L
    end subroutine set_frozen_soil_isothermal

    !> One probe state of the contract, evaluated and reported.
    subroutine probe_frozen(label, temperature, pressure, response, failures)
        implicit none
        character(len=*), intent(in) :: label
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        type(type_mixed_frozen_response), intent(inout) :: response
        integer(int32), intent(inout) :: failures

        character(len=200) :: line
        logical :: ok

        ok = .false.
        call eval_frozen_response(1_int32, temperature, pressure, response, ok)
        if (.not. ok) then
            call report_line("  FAIL "//label//"  the contract could not evaluate this state")
            failures = failures + 1
            return
        end if
        write (line, '(a,a,a,f7.3,a,es12.5,a,es12.5,a,es12.5,a,es12.5,a,es12.5,a,es12.5)') &
            "  state ", label, "  T ", temperature, " C  p_w ", pressure, &
            "  W ", response%stored_water, "  theta_w ", response%liquid_content, &
            "  theta_i ", response%ice_content, "  p_l ", response%liquid_pressure, &
            "  k_r ", response%relative_conductivity
        call report_line(trim(line))
    end subroutine probe_frozen

    !> The four states, and what has to hold at each of them.
    subroutine check_frozen_states(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_frozen_response) :: unfrozen, freezing, saturated, limited
        type(type_mixed_frozen_tangents) :: tangents
        real(real64) :: recombined, compressive
        logical :: ok

        call probe_frozen("unfrozen ", FROZEN_WARM, FROZEN_FIELD_PRESSURE, unfrozen, failures)
        call probe_frozen("freezing ", FROZEN_MILD_COLD, FROZEN_FIELD_PRESSURE, freezing, failures)
        call probe_frozen("saturated", FROZEN_WARM, FROZEN_POSITIVE_PRESSURE, saturated, failures)
        call probe_frozen("limited  ", FROZEN_DEEP_COLD, FROZEN_WET_PRESSURE, limited, failures)

        ! The stored quantity is the one the balance transports: no phase may be
        ! missing from it and none may be counted twice. Recombined from the
        ! contract's own reported contents and densities, which is the density
        ! conversion the ice term carries.
        compressive = FROZEN_SPECIFIC_STORAGE/SPECIFIC_WEIGHT*max(0.0d0, FROZEN_FIELD_PRESSURE)
        recombined = freezing%liquid_content &
                     + freezing%ice_density/freezing%water_density*freezing%ice_content &
                     + freezing%vapor_content + compressive
        call check("W is the phase sum, freezing state  ", freezing%stored_water - recombined, 1.0d-14, failures)

        compressive = FROZEN_SPECIFIC_STORAGE/SPECIFIC_WEIGHT*max(0.0d0, FROZEN_POSITIVE_PRESSURE)
        recombined = saturated%liquid_content &
                     + saturated%ice_density/saturated%water_density*saturated%ice_content &
                     + saturated%vapor_content + compressive
        call check("W is the phase sum, saturated state ", saturated%stored_water - recombined, 1.0d-14, failures)

        ! Freezing has to be doing something at the cold states, or the rest of
        ! this case would pass on an unfrozen soil.
        call check_exceeds("ice at the freezing state          ", freezing%ice_content, 1.0d-3, failures)
        call check_exceeds("ice at the deep cold state         ", limited%ice_content, 1.0d-2, failures)
        call check("no ice at the warm state           ", unfrozen%ice_content, 0.0d0, failures)

        ! Cryosuction, once: the liquid pressure is the published effective
        ! suction, equal to p_w where nothing is frozen and far below it where
        ! the ice is.
        call check("p_l equals p_w unfrozen, relative  ", &
                   (unfrozen%liquid_pressure - FROZEN_FIELD_PRESSURE)/abs(FROZEN_FIELD_PRESSURE), &
                   1.0d-12, failures)
        call check("p_l equals p_w saturated, relative ", &
                   (saturated%liquid_pressure - FROZEN_POSITIVE_PRESSURE)/abs(FROZEN_POSITIVE_PRESSURE), &
                   1.0d-12, failures)
        call check_exceeds("suction gained by freezing, Pa     ", &
                           FROZEN_FIELD_PRESSURE - freezing%liquid_pressure, 1.0d4, failures)

        ! The pore volume binds: the phases cannot ask for more room than there is.
        call check_exceeds("pore volume left at deep cold      ", &
                           SOIL_THETA_S - limited%liquid_content - limited%ice_content + 1.0d0, 1.0d0, failures)

        ! A saturated cell stores through the compressive branch, which is the
        ! whole of its capacity once the retention curve has flattened.
        ok = .false.
        call eval_frozen_tangents(1_int32, FROZEN_WARM, FROZEN_POSITIVE_PRESSURE, tangents, ok)
        if (ok) then
            call check_exceeds("saturated capacity over S_s/rho g  ", &
                               tangents%dwater_dpressure/(FROZEN_SPECIFIC_STORAGE/SPECIFIC_WEIGHT), &
                               0.999d0, failures)
        else
            call report_line("  FAIL saturated tangents could not be evaluated")
            failures = failures + 1
        end if

        ! Every state has a usable mobility, and none exceeds the saturated one.
        call check_exceeds("k_r unfrozen                       ", unfrozen%relative_conductivity, 1.0d-12, failures)
        call check_exceeds("k_r freezing                       ", freezing%relative_conductivity, 1.0d-30, failures)
        call check_exceeds("k_r deep cold                      ", limited%relative_conductivity, 1.0d-30, failures)
        call check("k_r saturated is one               ", saturated%relative_conductivity - 1.0d0, 1.0d-12, failures)
        call check_exceeds("k_r lost to ice, decades           ", &
                           log10(unfrozen%relative_conductivity/freezing%relative_conductivity), 1.0d0, failures)
    end subroutine check_frozen_states

    !> States off the domain are refused, and the refusal returns.
    !>
    !> A Newton iterate is not a physical state. The freezing model divides by
    !> the absolute temperature and ends the program below it, so the contract
    !> in front of it has to answer `ok = .false.` for an iterate that left the
    !> domain -- the answer the residual kernels turn into a NaN and the step
    !> controller turns into a shorter step. If this ever regresses the symptom
    !> is not a failed assertion: the run dies where the iterate was formed, so
    !> the case below reaching its last line is itself the measurement.
    subroutine check_frozen_domain(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_frozen_response) :: response
        type(type_mixed_frozen_tangents) :: tangents
        real(real64) :: nan
        logical :: ok

        nan = ieee_value(1.0d0, ieee_quiet_nan)

        ! Below absolute zero. The iterate that killed the 0.25 mm runs was
        ! -292.9 C, so the excursion this guards against is a real one and not
        ! a contrived extreme.
        ok = .true.
        call eval_frozen_response(1_int32, -292.9d0, FROZEN_FIELD_PRESSURE, response, ok)
        call check_integer("a state below absolute zero is refused  ", merge(1, 0, .not. ok), 1, failures)

        ok = .true.
        call eval_frozen_tangents(1_int32, -292.9d0, FROZEN_FIELD_PRESSURE, tangents, ok)
        call check_integer("tangents below absolute zero refused    ", merge(1, 0, .not. ok), 1, failures)

        ! Exactly at absolute zero, where the division is by zero rather than by
        ! a negative number.
        ok = .true.
        call eval_frozen_response(1_int32, FROZEN_ABSOLUTE_ZERO, FROZEN_FIELD_PRESSURE, response, ok)
        call check_integer("a state at absolute zero is refused     ", merge(1, 0, .not. ok), 1, failures)

        ok = .true.
        call eval_frozen_response(1_int32, nan, FROZEN_FIELD_PRESSURE, response, ok)
        call check_integer("a nonfinite temperature is refused      ", merge(1, 0, .not. ok), 1, failures)

        ok = .true.
        call eval_frozen_response(1_int32, FROZEN_MILD_COLD, nan, response, ok)
        call check_integer("a nonfinite pressure is refused         ", merge(1, 0, .not. ok), 1, failures)

        ! The guard must refuse only what is off the domain. A cold state the
        ! column really reaches still has to be answered, or the refusal would
        ! be shrinking every step instead of the ones that need it.
        ok = .false.
        call eval_frozen_response(1_int32, FROZEN_DEEP_COLD, FROZEN_FIELD_PRESSURE, response, ok)
        call check_integer("a deeply cold physical state answers    ", merge(1, 0, ok), 1, failures)
    end subroutine check_frozen_domain

    !> The unfrozen limit against the already-verified isothermal closure.
    !>
    !> The isothermal path is the one the twenty-five existing cases measure, so
    !> agreeing with it where there is no ice is what says the frozen contract
    !> did not replace the retention or conductivity law on its way through the
    !> governing models.
    subroutine check_frozen_unfrozen_limit(isothermal, failures)
        implicit none
        type(type_mixed_material), intent(inout) :: isothermal(:)
        integer(int32), intent(inout) :: failures

        type(type_mixed_frozen_response) :: frozen_response
        real(real64) :: storage, capacity, relative, difference
        logical :: ok

        call set_frozen_soil_isothermal(isothermal(1), FROZEN_MATERIAL_ID)
        call install_storage_laws(isothermal)
        call install_conductivity_laws(isothermal)

        ok = .false.
        call eval_frozen_response(1_int32, FROZEN_WARM, FROZEN_FIELD_PRESSURE, frozen_response, ok)
        if (.not. ok) then
            call report_line("  FAIL the unfrozen state could not be evaluated")
            failures = failures + 1
            call clear_storage_laws()
            call clear_conductivity_laws()
            return
        end if

        storage = 0.0d0
        capacity = 0.0d0
        call eval_storage(1_int32, FROZEN_FIELD_PRESSURE, storage, capacity)

        ! The isothermal law is retention plus compressive storage; the frozen
        ! one is that plus the vapour content, which the isothermal path does
        ! not carry. Asserting the difference IS the vapour term is stronger
        ! than asserting the two are close.
        difference = frozen_response%stored_water - storage - frozen_response%vapor_content
        call check("W minus the isothermal W and vapour", difference, 1.0d-12, failures)
        call report_value("vapour content at this state       ", frozen_response%vapor_content)

        relative = 0.0d0
        ok = .false.
        call eval_relative_conductivity(1_int32, FROZEN_FIELD_PRESSURE, relative, ok)
        if (ok) then
            call check("k_r against the isothermal closure ", &
                       (frozen_response%relative_conductivity - relative)/relative, 1.0d-12, failures)
        else
            call report_line("  FAIL the isothermal closure has no conductivity at this state")
            failures = failures + 1
        end if

        call clear_storage_laws()
        call clear_conductivity_laws()
    end subroutine check_frozen_unfrozen_limit

    !> The analytic tangents against central differences of the contract itself.
    !>
    !> Taken at the unfrozen state and at a freezing state that is on a smooth
    !> branch: not at the freezing point, not at the saturation pressure and not
    !> at the pore-volume limit, where the closure's own derivative is one-sided.
    subroutine check_frozen_tangents(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        call check_frozen_tangents_at("unfrozen", FROZEN_WARM, FROZEN_FIELD_PRESSURE, failures)
        call check_frozen_tangents_at("freezing", FROZEN_MILD_COLD, FROZEN_FIELD_PRESSURE, failures)
        call check_frozen_tangents_at("saturated", FROZEN_WARM, FROZEN_POSITIVE_PRESSURE, failures)
    end subroutine check_frozen_tangents

    !> One state's tangents, differenced.
    subroutine check_frozen_tangents_at(label, temperature, pressure, failures)
        implicit none
        character(len=*), intent(in) :: label
        real(real64), intent(in) :: temperature
        real(real64), intent(in) :: pressure
        integer(int32), intent(inout) :: failures

        type(type_mixed_frozen_tangents) :: tangents
        type(type_mixed_frozen_response) :: ahead, behind
        real(real64) :: pressure_step, temperature_step
        real(real64) :: differenced, scale
        logical :: ok

        ! Large enough that the difference is not roundoff, small enough that the
        ! curvature of these laws is not what is being measured.
        pressure_step = 2.0d0
        temperature_step = 1.0d-3

        ok = .false.
        call eval_frozen_tangents(1_int32, temperature, pressure, tangents, ok)
        if (.not. ok) then
            call report_line("  FAIL "//label//" tangents could not be evaluated")
            failures = failures + 1
            return
        end if

        call eval_frozen_response(1_int32, temperature, pressure + pressure_step, ahead, ok)
        if (.not. ok) return
        call eval_frozen_response(1_int32, temperature, pressure - pressure_step, behind, ok)
        if (.not. ok) return

        differenced = (ahead%stored_water - behind%stored_water)/(2.0d0*pressure_step)
        scale = max(abs(differenced), tiny(1.0d0))
        call check("dW/dp "//label//", relative        ", &
                   (tangents%dwater_dpressure - differenced)/scale, 2.0d-5, failures)
        call check_exceeds("dW/dp "//label//", magnitude       ", differenced, 1.0d-12, failures)

        differenced = (ahead%enthalpy - behind%enthalpy)/(2.0d0*pressure_step)
        scale = max(abs(differenced), 1.0d-8)
        call check("dU/dp "//label//", relative        ", &
                   (tangents%denthalpy_dpressure - differenced)/scale, 2.0d-4, failures)
        differenced = (log(ahead%relative_conductivity) - log(behind%relative_conductivity))/(2.0d0*pressure_step)
        scale = max(abs(differenced), 1.0d-12)
        call check("dlogK/dp "//label//", relative     ", &
                   (tangents%dlog_conductivity_dpressure - differenced)/scale, 2.0d-4, failures)

        differenced = (ahead%liquid_pressure - behind%liquid_pressure)/(2.0d0*pressure_step)
        scale = max(abs(differenced), tiny(1.0d0))
        call check("dp_l/dp "//label//", relative      ", &
                   (tangents%dliquid_dpressure - differenced)/scale, 2.0d-5, failures)

        call eval_frozen_response(1_int32, temperature + temperature_step, pressure, ahead, ok)
        if (.not. ok) return
        call eval_frozen_response(1_int32, temperature - temperature_step, pressure, behind, ok)
        if (.not. ok) return

        differenced = (ahead%enthalpy - behind%enthalpy)/(2.0d0*temperature_step)
        scale = max(abs(differenced), tiny(1.0d0))
        call check("dU/dT "//label//", relative        ", &
                   (tangents%denthalpy_dtemperature - differenced)/scale, 1.0d-4, failures)

        differenced = (ahead%stored_water - behind%stored_water)/(2.0d0*temperature_step)
        scale = max(abs(differenced), abs(tangents%dwater_dtemperature), tiny(1.0d0))
        call check("dW/dT "//label//", relative        ", &
                   (tangents%dwater_dtemperature - differenced)/scale, 2.0d-4, failures)

        differenced = (log(ahead%relative_conductivity) - log(behind%relative_conductivity))/(2.0d0*temperature_step)
        ! Absolute tolerance 2e-10 /K where the exact derivative vanishes.
        scale = max(abs(differenced), 1.0d-6)
        call check("dlogK/dT "//label//", relative     ", &
                   (tangents%dlog_conductivity_dtemperature - differenced)/scale, 2.0d-4, failures)

        differenced = (ahead%liquid_pressure - behind%liquid_pressure)/(2.0d0*temperature_step)
        scale = max(abs(differenced), abs(tangents%dliquid_dtemperature), tiny(1.0d0))
        call check("dp_l/dT "//label//", relative      ", &
                   (tangents%dliquid_dtemperature - differenced)/scale, 2.0d-4, failures)
    end subroutine check_frozen_tangents_at

    ! ================================================================
    ! Case: the coupled element on axis-aligned rectangles
    ! ================================================================

    !> What a rectangle has to reproduce exactly, and what it may not invent.
    !>
    !> A manufactured state is projected into the three spaces of the coupled
    !> column and read back through the same reconstruction the conservation
    !> audit uses. Each piece of it is exactly representable on this mesh -- the
    !> temperature is bilinear, the flux lies in the tensor-product RT0 space
    !> and has a constant nonzero divergence, the pressure is cell constant --
    !> so every comparison below is against zero rather than against a
    !> discretisation error. A wrong reference cell, a wrong Piola pushforward,
    !> a wrong facet orientation or a simplex element on a tensor cell each
    !> break at least one of them.
    !>
    !> The case never solves. It registers no boundary, because inserting
    !> boundary values would overwrite the projected flux on the perimeter, and
    !> what is under test is the element and the geometry rather than a step.
    subroutine case_rectangle(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_column) :: column
        type(type_mixed_frozen_material) :: frozen(1)
        real(real64) :: conductivity(MIXED_DIMENSION, MIXED_DIMENSION, 1)
        DM :: dm
        Vec :: local
        PetscErrorCode :: ierr
        logical :: projected

        call report_line("case rectangle: Q1/RT0/P0 on an elongated axis-aligned rectangle")

        ierr = 0
        call build_rectangle_mesh(dm, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the rectangle mesh could not be built")
            failures = failures + 1
            return
        end if

        call set_frozen_soil(frozen(1), FROZEN_MATERIAL_ID)
        conductivity = 0.0d0
        conductivity(1, 1, 1) = FROZEN_CONDUCTIVITY
        conductivity(2, 2, 1) = FROZEN_CONDUCTIVITY

        call column%set_materials(frozen, conductivity)
        call column%set_gravity([0.0d0, -STANDARD_GRAVITY], REFERENCE_DENSITY)
        ! Unit scales: what is read back is then in the units it was written in,
        ! so no comparison here passes through a scaling first.
        call column%set_scales(1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)

        call column%initialize(dm, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the rectangle column could not be initialised")
            failures = failures + 1
            call DMDestroy(dm, ierr)
            return
        end if

        call check_rectangle_layout(column, failures)

        projected = .false.
        call project_rectangle_state(column, projected)
        call check_logical("the manufactured state projected ", projected, .true., failures)
        if (projected) then
            call DMGetLocalVector(column%get_dm(), local, ierr)
            if (ierr == 0) call DMGlobalToLocalBegin(column%get_dm(), column%solution, INSERT_VALUES, local, ierr)
            if (ierr == 0) call DMGlobalToLocalEnd(column%get_dm(), column%solution, INSERT_VALUES, local, ierr)
            if (ierr /= 0) then
                call report_line("  FAIL the projected state could not be localised")
                failures = failures + 1
            else
                call check_rectangle_cells(column, local, failures)
                call check_rectangle_continuity(column, local, failures)
            end if
            call DMRestoreLocalVector(column%get_dm(), local, ierr)
        end if

        call finish_frozen_column(column, dm)
    end subroutine case_rectangle

    !> A quadrilateral that is not a rectangle has to stop the run.
    !>
    !> The audit inverts one affine map per cell and the element loop pushes
    !> the basis forward with `isAffine`, and neither is the map a warped cell
    !> actually has. A solver that accepted this mesh would report a state
    !> through a geometry it does not have, so the refusal is the behaviour
    !> under test.
    !>
    !> This covers the warped end of the gate only. The gate itself is
    !> narrower than the affine condition: a parallelogram has an exact affine
    !> map and is refused too, because the supported geometry is the
    !> axis-aligned rectangle.
    subroutine case_reject_warped()
        implicit none

        type(type_mixed_column) :: column
        type(type_mixed_frozen_material) :: frozen(1)
        real(real64) :: conductivity(MIXED_DIMENSION, MIXED_DIMENSION, 1)
        DM :: dm

        call report_line("case reject_warped: a quadrilateral that is not a rectangle has to stop the run")

        call build_warped_quad(dm, ierr)
        call label_cells(dm, 1.0d30, FROZEN_MATERIAL_ID, FROZEN_MATERIAL_ID, ierr)
        call label_facets(dm, COLUMN_CLOSED_ID, -1.0d3, 1.0d3, -1.0d3, 1.0d3, ierr)

        call set_frozen_soil(frozen(1), FROZEN_MATERIAL_ID)
        conductivity = 0.0d0
        conductivity(1, 1, 1) = FROZEN_CONDUCTIVITY
        conductivity(2, 2, 1) = FROZEN_CONDUCTIVITY

        call column%set_materials(frozen, conductivity)
        call column%set_gravity([0.0d0, -STANDARD_GRAVITY], REFERENCE_DENSITY)
        call column%set_scales(1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
        ! This must not return.
        call column%initialize(dm, ierr)

        call report_line("  FAIL a warped quadrilateral was accepted")
        error stop 2
    end subroutine case_reject_warped

    !> One quadrilateral whose fourth vertex is off the parallelogram.
    subroutine build_warped_quad(dm, ierr)
        implicit none
        DM, intent(inout) :: dm
        PetscErrorCode, intent(inout) :: ierr

        PetscInt :: cells(4)
        PetscReal :: coordinates(8)

        cells = [0_PETSC_INT_KIND, 1_PETSC_INT_KIND, 2_PETSC_INT_KIND, 3_PETSC_INT_KIND]
        ! (0,0), (1,0), (1.4,1), (0,1): convex and counter-clockwise, so the
        ! mesh itself is valid, and the map from the reference square is
        ! bilinear rather than affine.
        coordinates = real([0.0d0, 0.0d0, 1.0d0, 0.0d0, 1.4d0, 1.0d0, 0.0d0, 1.0d0], PETSC_REAL_KIND)

        call DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, 2_PETSC_INT_KIND, 1_PETSC_INT_KIND, &
                                           4_PETSC_INT_KIND, 4_PETSC_INT_KIND, PETSC_TRUE, cells, &
                                           2_PETSC_INT_KIND, coordinates, dm, ierr)
    end subroutine build_warped_quad

    !> The two-cell elongated strip, with the column's own labels on it.
    subroutine build_rectangle_mesh(dm, ierr)
        implicit none
        DM, intent(inout) :: dm
        PetscErrorCode, intent(inout) :: ierr

        real(real64) :: reach

        reach = 1.0d3*max(RECT_WIDTH, RECT_HEIGHT)

        call build_quad_strip(RECT_NX, RECT_NZ, RECT_WIDTH, RECT_HEIGHT, dm, ierr)
        if (ierr /= 0) return
        call label_cells(dm, 0.5d0*RECT_WIDTH, FROZEN_MATERIAL_ID, FROZEN_MATERIAL_ID, ierr)
        if (ierr /= 0) return
        call label_facets(dm, COLUMN_CLOSED_ID, -reach, reach, -reach, reach, ierr)
        if (ierr /= 0) return
        call distribute(dm, ierr)
    end subroutine build_rectangle_mesh

    !> The element a rectangle carries: four flux, one pressure, four nodal.
    subroutine check_rectangle_layout(column, failures)
        implicit none
        type(type_mixed_column), intent(inout) :: column
        integer(int32), intent(inout) :: failures

        DM :: dm
        PetscSection :: section
        type(c_ptr) :: dm_handle
        integer(c_int) :: num_basis, num_components, cerr
        PetscInt :: point, v_start, v_end, e_start, e_end, c_start, c_end, dofs
        integer(int32) :: vertex_temperature, vertex_flux, facet_flux, cell_flux
        integer(int32) :: facet_pressure, cell_pressure, cells

        dm = column%get_dm()
        dm_handle = transfer(dm%v, dm_handle)

        cerr = 0
        num_basis = 0
        num_components = 0
        call petscfe_field_size(dm_handle, int(COLUMN_FIELD_TEMPERATURE, c_int), num_basis, num_components, cerr)
        call check_integer("temperature basis per rectangle         ", &
                           int(num_basis, int32), RECT_TEMPERATURE_DOFS, failures)
        call check_integer("temperature components                  ", int(num_components, int32), 1, failures)
        call petscfe_field_size(dm_handle, int(COLUMN_FIELD_FLUX, c_int), num_basis, num_components, cerr)
        call check_integer("flux basis per rectangle                ", &
                           int(num_basis, int32), RECT_FLUX_DOFS, failures)
        call check_integer("flux components                         ", int(num_components, int32), 2, failures)
        call petscfe_field_size(dm_handle, int(COLUMN_FIELD_PRESSURE, c_int), num_basis, num_components, cerr)
        call check_integer("pressure basis per rectangle            ", &
                           int(num_basis, int32), RECT_PRESSURE_DOFS, failures)
        call check_integer("pressure components                     ", int(num_components, int32), 1, failures)
        call check_integer("the element sizes were readable         ", int(cerr, int32), 0, failures)

        call DMGetLocalSection(dm, section, ierr)
        call DMPlexGetDepthStratum(dm, 0_PETSC_INT_KIND, v_start, v_end, ierr)
        call DMPlexGetHeightStratum(dm, 1_PETSC_INT_KIND, e_start, e_end, ierr)
        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)

        vertex_temperature = 0
        vertex_flux = 0
        facet_flux = 0
        facet_pressure = 0
        cell_flux = 0
        cell_pressure = 0
        do point = v_start, v_end - 1
            call PetscSectionGetFieldDof(section, point, COLUMN_FIELD_TEMPERATURE, dofs, ierr)
            vertex_temperature = max(vertex_temperature, int(dofs, int32))
            call PetscSectionGetFieldDof(section, point, COLUMN_FIELD_FLUX, dofs, ierr)
            vertex_flux = vertex_flux + int(dofs, int32)
        end do
        do point = e_start, e_end - 1
            call PetscSectionGetFieldDof(section, point, COLUMN_FIELD_FLUX, dofs, ierr)
            facet_flux = max(facet_flux, int(dofs, int32))
            call PetscSectionGetFieldDof(section, point, COLUMN_FIELD_PRESSURE, dofs, ierr)
            facet_pressure = facet_pressure + int(dofs, int32)
        end do
        cells = 0
        do point = c_start, c_end - 1
            if (.not. cell_is_owned(dm, point)) cycle
            cells = cells + 1
            call PetscSectionGetFieldDof(section, point, COLUMN_FIELD_FLUX, dofs, ierr)
            cell_flux = cell_flux + int(dofs, int32)
            call PetscSectionGetFieldDof(section, point, COLUMN_FIELD_PRESSURE, dofs, ierr)
            cell_pressure = max(cell_pressure, int(dofs, int32))
        end do

        call reduce_count(cells)
        call reduce_count(vertex_flux)
        call reduce_count(cell_flux)
        call reduce_count(facet_pressure)

        call check_integer("rectangles the strip carries            ", cells, RECT_NX*RECT_NZ, failures)
        call check_integer("temperature degrees of freedom per node ", vertex_temperature, 1, failures)
        call check_integer("flux degrees of freedom on vertices     ", vertex_flux, 0, failures)
        call check_integer("flux degrees of freedom per facet       ", facet_flux, 1, failures)
        call check_integer("flux degrees of freedom in cells        ", cell_flux, 0, failures)
        call check_integer("pressure degrees of freedom on facets   ", facet_pressure, 0, failures)
        call check_integer("pressure degrees of freedom per cell    ", cell_pressure, 1, failures)
    end subroutine check_rectangle_layout

    !> Interpolate the manufactured state into the three coupled spaces.
    subroutine project_rectangle_state(column, projected)
        implicit none
        type(type_mixed_column), intent(inout) :: column
        logical, intent(inout) :: projected

        DM :: dm
        type(c_funptr) :: functions(3)
        type(c_ptr) :: dm_handle, solution_handle
        integer(c_int) :: cerr

        dm = column%get_dm()
        dm_handle = transfer(dm%v, dm_handle)
        solution_handle = transfer(column%solution%v, solution_handle)

        functions(1) = c_funloc(manufactured_temperature)
        functions(2) = c_funloc(manufactured_flux)
        functions(3) = c_funloc(manufactured_pressure)

        cerr = 0
        call petscds_project_function(dm_handle, 0.0_c_double, functions, 1_c_int, solution_handle, cerr)
        projected = (cerr == 0)
    end subroutine project_rectangle_state

    !> Read the manufactured state back, cell by cell.
    subroutine check_rectangle_cells(column, local, failures)
        implicit none
        type(type_mixed_column), intent(inout) :: column
        Vec, intent(in) :: local
        integer(int32), intent(inout) :: failures

        DM :: dm
        type(c_ptr) :: dm_handle, vec_handle
        PetscInt :: cell, c_start, c_end
        PetscInt :: facets(MAX_CELL_FACETS)
        PetscReal :: measure, centroid(3), normal(3)
        PetscReal :: facet_measure, facet_centroid(3), facet_normal(3)
        real(real64) :: moments(MAX_CELL_FACETS)
        real(c_double) :: sample(MIXED_DIMENSION, RECTANGLE_SAMPLES)
        real(c_double) :: reference(MIXED_DIMENSION, RECTANGLE_SAMPLES)
        real(c_double) :: values(MIXED_DIMENSION, RECTANGLE_SAMPLES)
        real(c_double) :: divergence(RECTANGLE_SAMPLES)
        real(real64) :: exact(MIXED_DIMENSION), outward(MIXED_DIMENSION)
        real(real64) :: cell_width, cell_height, volume, orientation
        real(real64) :: flux_reference, moment_reference
        real(real64) :: worst_x, worst_z, worst_divergence, worst_area
        real(real64) :: worst_moment, worst_theorem, worst_temperature, worst_pressure
        real(real64) :: pressure, temperature, integral, outflow, exact_moment
        !> The largest normal moment the reconstruction reads on a facet that
        !> is both interior and horizontal. It is what makes the continuity
        !> comparison a statement: on a single layer of cells the only interior
        !> facets are vertical, where this field's trace is identically zero,
        !> and two sides agreeing at zero would agree whatever the element did.
        real(real64) :: transported_moment
        integer(int32) :: num_facets, miscounted_facets, f, q, d, s
        integer(c_int) :: cerr

        dm = column%get_dm()
        dm_handle = transfer(dm%v, dm_handle)
        vec_handle = transfer(local%v, vec_handle)

        cell_width = RECT_WIDTH/real(RECT_NX, real64)
        cell_height = RECT_HEIGHT/real(RECT_NZ, real64)
        ! The scales every comparison below is relative to: the largest flux the
        ! manufactured field reaches, and the largest facet moment it makes.
        flux_reference = MANUFACTURED_FLUX_SLOPE*RECT_HEIGHT
        moment_reference = flux_reference*cell_width

        worst_x = 0.0d0
        worst_z = 0.0d0
        worst_divergence = 0.0d0
        worst_area = 0.0d0
        worst_moment = 0.0d0
        worst_theorem = 0.0d0
        worst_temperature = 0.0d0
        worst_pressure = 0.0d0
        transported_moment = 0.0d0
        miscounted_facets = 0

        call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        do cell = c_start, c_end - 1
            if (.not. cell_is_owned(dm, cell)) cycle

            call DMPlexComputeCellGeometryFVM(dm, cell, measure, centroid, normal, ierr)
            volume = real(measure, real64)
            worst_area = max(worst_area, abs(volume - cell_width*cell_height))

            ! Four points strictly inside the cell, off both mid-lines, so a
            ! field that is only right at the centre does not pass.
            s = 0
            do q = -1, 1, 2
                do d = -1, 1, 2
                    s = s + 1
                    sample(1, s) = real(real(centroid(1), real64) + 0.25d0*real(q, real64)*cell_width, c_double)
                    sample(2, s) = real(real(centroid(2), real64) + 0.25d0*real(d, real64)*cell_height, c_double)
                end do
            end do

            cerr = 0
            call petscfe_reference_coordinates(dm_handle, int(cell, c_int), sample, reference, cerr)
            if (cerr == 0) then
                call petscfe_evaluate_cell(dm_handle, vec_handle, int(COLUMN_FIELD_FLUX, c_int), &
                                           int(cell, c_int), reference, values, divergence, cerr)
            end if
            if (cerr /= 0) then
                call report_line("  FAIL the projected flux could not be evaluated")
                failures = failures + 1
                return
            end if

            do s = 1, RECTANGLE_SAMPLES
                exact = manufactured_flux_at(real(sample(1, s), real64), real(sample(2, s), real64))
                worst_x = max(worst_x, abs(real(values(1, s), real64) - exact(1)))
                worst_z = max(worst_z, abs(real(values(2, s), real64) - exact(2)))
                worst_divergence = max(worst_divergence, &
                                       abs(real(divergence(s), real64) - MANUFACTURED_FLUX_SLOPE))
            end do

            ! Each facet's own moment, against the integral of the exact field.
            ! The normal trace is constant on an axis-aligned facet of this
            ! field, so the centroid value times the length is exact.
            call reconstruct_cell_facet_flux(dm, local, 1.0d0, cell, facets, moments, num_facets, ierr, &
                                             field=COLUMN_FIELD_FLUX)
            miscounted_facets = miscounted_facets + abs(num_facets - 4)
            outflow = 0.0d0
            do f = 1, num_facets
                call DMPlexComputeCellGeometryFVM(dm, facets(f), facet_measure, facet_centroid, facet_normal, ierr)
                orientation = 0.0d0
                do d = 1, MIXED_DIMENSION
                    orientation = orientation + real(facet_normal(d), real64)* &
                                  (real(facet_centroid(d), real64) - real(centroid(d), real64))
                end do
                outward = real(facet_normal(1:MIXED_DIMENSION), real64)
                if (orientation < 0.0d0) outward = -outward
                exact = manufactured_flux_at(real(facet_centroid(1), real64), real(facet_centroid(2), real64))
                exact_moment = dot_product(exact, outward)*real(facet_measure, real64)
                worst_moment = max(worst_moment, abs(moments(f) - exact_moment))
                outflow = outflow + moments(f)
                ! Horizontal, and off both the top and the bottom of the strip.
                if (abs(outward(2)) > 0.5d0 .and. &
                    real(facet_centroid(2), real64) > 0.25d0*cell_height .and. &
                    real(facet_centroid(2), real64) < RECT_HEIGHT - 0.25d0*cell_height) then
                    transported_moment = max(transported_moment, abs(moments(f)))
                end if
            end do

            call reconstruct_cell_divergence(dm, local, 1.0d0, cell, integral, volume, ierr, &
                                             field=COLUMN_FIELD_FLUX)
            worst_divergence = max(worst_divergence, &
                                   abs(integral - MANUFACTURED_FLUX_SLOPE*volume)/volume)
            worst_theorem = max(worst_theorem, abs(outflow - integral))

            ! The centroid value of a field that is not cell constant: this is
            ! the reference point the reconstruction has to get right, and on a
            ! rectangle it is the origin and not the triangle's (-1/3,-1/3).
            call reconstruct_cell_pressure(dm, local, 1.0d0, cell, temperature, ierr, &
                                           field=COLUMN_FIELD_TEMPERATURE)
            worst_temperature = max(worst_temperature, &
                                    abs(temperature - manufactured_temperature_at(real(centroid(1), real64), &
                                                                                  real(centroid(2), real64))))
            call reconstruct_cell_pressure(dm, local, 1.0d0, cell, pressure, ierr, &
                                           field=COLUMN_FIELD_PRESSURE)
            worst_pressure = max(worst_pressure, abs(pressure - MANUFACTURED_P0))
        end do

        call reduce_max(worst_x)
        call reduce_max(worst_z)
        call reduce_max(worst_divergence)
        call reduce_max(worst_area)
        call reduce_max(worst_moment)
        call reduce_max(worst_theorem)
        call reduce_max(worst_temperature)
        call reduce_max(worst_pressure)
        call reduce_max(transported_moment)
        call reduce_count(miscounted_facets)

        call check_integer("cells whose facet count is not four     ", miscounted_facets, 0, failures)
        call check("cell area against the exact rectangle   ", &
                   worst_area/(cell_width*cell_height), 1.0d-14, failures)
        call check("the flux has no x component             ", worst_x/flux_reference, 1.0d-13, failures)
        call check("the flux z component against exact      ", worst_z/flux_reference, 1.0d-13, failures)
        call check("divergence against the exact constant   ", &
                   worst_divergence/MANUFACTURED_FLUX_SLOPE, 1.0d-13, failures)
        call check("every facet moment against exact        ", worst_moment/moment_reference, 1.0d-13, failures)
        call check("Piola divergence theorem defect         ", worst_theorem/moment_reference, 1.0d-13, failures)
        call check("centroid temperature, bilinear term too ", worst_temperature, 1.0d-13, failures)
        call check("cell pressure against exact             ", &
                   worst_pressure/MANUFACTURED_P0, 1.0d-14, failures)
        ! The exact value here is z times the facet length at mid-height, which
        ! is half the reference; the floor only has to exclude zero.
        call check_exceeds("interior horizontal facet moment        ", &
                           transported_moment/moment_reference, 1.0d-1, failures)
    end subroutine check_rectangle_cells

    !> The shared facet carries one normal trace, computed twice.
    subroutine check_rectangle_continuity(column, local, failures)
        implicit none
        type(type_mixed_column), intent(inout) :: column
        Vec, intent(in) :: local
        integer(int32), intent(inout) :: failures

        real(real64) :: mismatch, moment_reference
        integer(int32) :: interior_facets

        mismatch = 0.0d0
        interior_facets = 0
        moment_reference = MANUFACTURED_FLUX_SLOPE*RECT_HEIGHT*RECT_WIDTH/real(RECT_NX, real64)

        call audit_facet_continuity(column%get_dm(), local, 1.0d0, mismatch, interior_facets, ierr, &
                                    field=COLUMN_FIELD_FLUX)
        call check_integer("interior facets seen from both sides    ", &
                           interior_facets, (RECT_NX - 1)*RECT_NZ + RECT_NX*(RECT_NZ - 1), failures)
        call check("facet normal trace mismatch             ", mismatch/moment_reference, 1.0d-13, failures)
    end subroutine check_rectangle_continuity

    ! ================================================================
    ! Case: the coupled [T, q, p_w] freezing column
    ! ================================================================

    !> The first coupled increment: heat and water solved in one Newton step.
    !>
    !> A closed column at field capacity, warm at both ends, then cooled from
    !> the top. Nothing crosses any water boundary, so the total stored water is
    !> a conserved quantity of the exact problem and the cell defects are the
    !> discretisation's own. The warm interval is the control: it establishes
    !> what the column does when no freezing front exists, so the cooled
    !> interval's flux is measured against something rather than asserted.
    !>
    !> This is not a reproduction of the Mizoguchi experiment. It runs for hours
    !> of physical time, not the 12, 24 and 50 of the measurement, it carries no
    !> advective heat transport, and no observed profile is compared against.
    subroutine case_frozen_column(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_column) :: column
        type(type_mixed_column_report) :: start, control, cooled, retried
        type(type_mixed_frozen_material) :: frozen(1)
        real(real64) :: conductivity(MIXED_DIMENSION, MIXED_DIMENSION, 1)
        DM :: dm
        PetscErrorCode :: ierr
        logical :: accepted
        integer(int32) :: step
        real(real64) :: cooling_step, cooling_span
        type(type_mixed_column_progress) :: progress

        call report_line("case frozen_column: one monolithic H1/RT0/P0 step of heat and water")

        ! The registered case runs its own ladder; the options exist so a run can
        ! ask the same column a different step question without a rebuild.
        cooling_step = read_real_option("-mixed_column_cooling_step", COLUMN_COOLING_STEP)
        cooling_span = read_real_option("-mixed_column_cooling_span", COLUMN_COOLING_SPAN)

        ierr = 0
        call build_frozen_column_mesh(dm, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the column mesh could not be built")
            failures = failures + 1
            return
        end if

        call set_frozen_soil(frozen(1), FROZEN_MATERIAL_ID)
        conductivity = 0.0d0
        conductivity(1, 1, 1) = FROZEN_CONDUCTIVITY
        conductivity(2, 2, 1) = FROZEN_CONDUCTIVITY

        call column%set_materials(frozen, conductivity)
        call column%set_gravity([0.0d0, -STANDARD_GRAVITY], REFERENCE_DENSITY)
        ! The flux unit is what an unfrozen cell of this soil carries under a
        ! unit gradient; the pressure unit the suction the column starts at; the
        ! temperature unit one kelvin, and the heat row's divisor the
        ! conduction a kelvin makes over the column's own height.
        call column%set_scales(FROZEN_CONDUCTIVITY, abs(FROZEN_FIELD_PRESSURE), COLUMN_HEIGHT, &
                               1.0d0, 1.0d0/COLUMN_HEIGHT**2)
        call column%add_closed_boundary(COLUMN_CLOSED_ID)
        call column%add_temperature_boundary(COLUMN_TOP_ID, COLUMN_WARM)
        call column%add_temperature_boundary(COLUMN_BOTTOM_ID, COLUMN_WARM)

        call column%initialize(dm, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the coupled column could not be initialised")
            failures = failures + 1
            call DMDestroy(dm, ierr)
            return
        end if

        ! Hydrostatic about the mid-height, which is the state a closed column
        ! is at rest in. A uniform pressure would drain under gravity and the
        ! control interval would be measuring that instead.
        call column%set_initial_state(COLUMN_WARM, FROZEN_FIELD_PRESSURE, ierr, &
                                      hydrostatic_height=0.5d0*COLUMN_HEIGHT)
        if (ierr == 0) call column%compute_report(start, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the initial state could not be published")
            failures = failures + 1
            call finish_frozen_column(column, dm)
            return
        end if
        call report_column("initial ", start)
        call check_integer("cells the column carries          ", start%owned_cells, &
                           strip_cell_count(read_cell_type_option("-column_cell_type"), &
                                            COLUMN_NX, COLUMN_NZ), failures)
        call check_integer("frozen cells at the start         ", start%frozen_cells, 0, failures)

        ! The control: warm at both ends, nothing to drive a flux.
        do step = 1, COLUMN_CONTROL_STEPS
            accepted = .false.
            call column%solve_step(COLUMN_CONTROL_STEP, accepted, ierr)
            if (.not. accepted) then
                call report_line("  FAIL a warm control step was refused")
                failures = failures + 1
                call finish_frozen_column(column, dm)
                return
            end if
        end do
        call column%compute_report(control, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the control report could not be measured")
            failures = failures + 1
            call finish_frozen_column(column, dm)
            return
        end if
        call report_column("control ", control)
        call check_integer("frozen cells after the control    ", control%frozen_cells, 0, failures)
        call check("water the control created, relative", &
                   (control%total_water - start%total_water)/start%total_water, 1.0d-10, failures)
        call check_exceeds("the pressure column stays open    ", control%weakest_column, 1.0d-3, failures)

        ! Cool the top. The bottom stays where it was, which is the
        ! bottom-warm variant of the experiment rather than an adiabatic base.
        !
        ! The cooling is asked for as a span of physical time rather than as a
        ! count of steps: the front the step has to resolve is not the one the
        ! warm control had, so the length that carries it is the controller's to
        ! find. A run that met the horizon is the assertion; the steps it took
        ! to get there are reported, not fixed.
        call column%set_boundary_temperature(COLUMN_TOP_ID, COLUMN_COLD)
        call set_column_cooling_control(column, cooling_step)
        call column%advance(cooling_span, progress, ierr)
        call report_progress("cooling ", progress)
        if (ierr /= 0) then
            call report_line("  FAIL the cooling horizon was not reached")
            failures = failures + 1
            call finish_frozen_column(column, dm)
            return
        end if
        call check("the cooling horizon was met, s    ", progress%reached - cooling_span, &
                   1.0d-9*cooling_span, failures)
        call column%compute_report(cooled, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the cooled report could not be measured")
            failures = failures + 1
            call finish_frozen_column(column, dm)
            return
        end if
        call report_column("cooled  ", cooled)

        call check_frozen_column_physics(start, control, cooled, failures)
        call check_frozen_column_rollback(column, cooled, retried, failures)

        call finish_frozen_column(column, dm)
    end subroutine case_frozen_column

    !> What the cooled column has to show, and what it may not.
    subroutine check_frozen_column_physics(start, control, cooled, failures)
        implicit none
        type(type_mixed_column_report), intent(in) :: start, control, cooled
        integer(int32), intent(inout) :: failures

        ! Freezing happened at all, or nothing below is discriminating.
        call check_exceeds("cells carrying ice after cooling   ", real(cooled%frozen_cells, real64), 1.0d0, failures)
        call check_exceeds("ice the column made, m2            ", cooled%total_ice, 1.0d-6, failures)
        call check_logical("the cold end really is cold      ", cooled%coldest < -1.0d0, .true., failures)

        ! Cryosuction: the liquid pressure of the frozen cells is far below the
        ! pore pressure the column started at, and the warm control never was.
        call check_exceeds("suction the freezing gained, Pa    ", &
                           FROZEN_FIELD_PRESSURE - cooled%lowest_liquid_pressure, 1.0d4, failures)
        call check("the control gained none, relative  ", &
                   (control%lowest_liquid_pressure - start%lowest_liquid_pressure)/abs(FROZEN_FIELD_PRESSURE), &
                   1.0d-3, failures)

        ! Orientation: the ice-bearing cells draw water upward, against gravity,
        ! and by more than they pass downward. The control, with no front, did
        ! not.
        call check_exceeds("upward flux where the ice is, m/s  ", cooled%upward_flux_frozen, 1.0d-16, failures)
        call check_exceeds("it beats the downward flux there   ", &
                           cooled%upward_flux_frozen/max(abs(cooled%downward_flux_frozen), tiny(1.0d0)), 1.0d0, failures)
        call check("the control moved no water, m/s    ", control%upward_flux, 1.0d-12, failures)

        ! Conservation: nothing crosses a closed boundary, so the stored water
        ! is the same as it started and every cell's backward-Euler balance
        ! closes against its own reconstructed facet fluxes.
        call check("water the cooling created, relative", &
                   (cooled%total_water - start%total_water)/start%total_water, 1.0d-9, failures)
        call check("water crossing the boundary, m2/s  ", cooled%boundary_rate, 1.0d-18, failures)
        call check("the domain balance, relative       ", &
                   cooled%domain_defect/start%total_water, 1.0d-9, failures)
        call check_exceeds("the audit had a flux to measure    ", cooled%throughput, 1.0d-20, failures)
        call check_exceeds("the audit retained storage change  ", cooled%max_storage_change, 1.0d-12, failures)
        call check("cell balance relative to transfer  ", &
                   cooled%max_cell_defect/max(cooled%max_storage_change, cooled%max_flux_term, tiny(1.0d0)), &
                   1.0d-9, failures)
        call check("shared facet flux mismatch, m2/s   ", cooled%max_facet_mismatch, 1.0d-18, failures)

        ! The Jacobian the steps above were taken with was not singular.
        call check_exceeds("weakest pressure column, cooled    ", cooled%weakest_column, 1.0d-6, failures)
    end subroutine check_frozen_column_physics

    !> A refused step leaves the accepted history exactly where it was.
    subroutine check_frozen_column_rollback(column, cooled, retried, failures)
        implicit none
        type(type_mixed_column), intent(inout) :: column
        type(type_mixed_column_report), intent(in) :: cooled
        type(type_mixed_column_report), intent(inout) :: retried
        integer(int32), intent(inout) :: failures

        PetscErrorCode :: ierr
        logical :: accepted

        accepted = .true.
        ierr = 0
        call column%solve_step(-1.0d0, accepted, ierr)
        call check_logical("a negative step is refused       ", accepted, .false., failures)
        call check_integer("and refused for that reason      ", int(ierr, int32), COLUMN_ERR_INVALID_STEP, failures)

        ierr = 0
        call column%compute_report(retried, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the state after the refusal could not be measured")
            failures = failures + 1
            return
        end if
        ! Exactly, not nearly: a rejected step that perturbed the history would
        ! make the next accepted step start from a state no report describes.
        call check("the water survived the refusal    ", retried%total_water - cooled%total_water, 0.0d0, failures)
        call check("the ice survived it too           ", retried%total_ice - cooled%total_ice, 0.0d0, failures)
    end subroutine check_frozen_column_rollback

    !> The step controller the cooling runs under.
    !>
    !> The iteration window and the three factors are the project's own
    !> `Conditions.json` values, so the column is stepped the way the
    !> application steps this problem rather than by a rule invented for a test.
    !> The ceiling is the first step: a controller free to grow past it would
    !> make the case measure a different span of the transient every time the
    !> solver's iteration count moved.
    subroutine set_column_cooling_control(column, first_step)
        implicit none
        type(type_mixed_column), intent(inout) :: column
        real(real64), intent(in) :: first_step

        type(type_mixed_column_step_control) :: control

        control%initial = first_step
        control%minimum = 1.0d0
        control%maximum = first_step
        control%iterations_low = 5
        control%iterations_high = 16
        control%grow = 1.25d0
        control%shrink = 0.8d0
        control%retry = 0.5d0
        call column%set_step_control(control)
    end subroutine set_column_cooling_control

    !> What an advance did, printed: the time it reached and what it cost.
    subroutine report_progress(label, progress)
        implicit none
        character(len=*), intent(in) :: label
        type(type_mixed_column_progress), intent(in) :: progress

        character(len=250) :: line

        write (line, '(a,a,a,es12.5,a,i0,a,i0,a,es12.5,a,es12.5,a,i0,a,i0)') &
            "  column ", label, " reached ", progress%reached, " s  accepted ", &
            progress%accepted_steps, "  rejected ", progress%rejected_steps, &
            "  shortest ", progress%shortest_step, "  next ", progress%next_step, &
            "  its max ", progress%max_iterations, "  total ", progress%total_iterations
        call report_line(trim(line))
    end subroutine report_progress

    !> The Mizoguchi freezing experiment on the coupled column.
    !>
    !> The conditions are the experiment's, not a toy: a cooling plate in
    !> contact with the top whose measured temperature the surface follows
    !> (`-mizo_top_bc convective` restores the project's coolant behind a
    !> surface resistance instead), a base leaking to a 6.7 C room, no water
    !> across any boundary, and the soil starting uniform at 6.7 C and
    !> -44683.6901 Pa. What the case asserts is what a run of it can establish
    !> by itself -- that the horizons were reached, that no water was created,
    !> that ice formed and drew water toward the cold end. Agreement with the
    !> measured profiles is not asserted here: the profiles are written out and
    !> compared by `scripts/validation/mizo_column_profile.py`, against a table
    !> whose digitisation this repository cannot verify.
    subroutine case_mizo_column(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_column) :: column
        type(type_mixed_column_report) :: start, reached(MAX_MIZO_HORIZONS)
        type(type_mixed_column_progress) :: progress
        real(real64) :: final_hours, elapsed, span
        real(real64) :: hours(MAX_MIZO_HORIZONS)
        integer(int32) :: horizons
        character(len=256) :: prefix
        character(len=256) :: history
        DM :: dm
        PetscErrorCode :: ierr
        integer(int32) :: h, done, mesh_cells

        call report_line("case mizo_column: the Mizoguchi freezing experiment, monolithic H1/RT0/P0 "// &
                         "(H1/H1 with -mizo_hydraulic_space primal)")

        call read_mizo_horizons(hours, horizons)
        final_hours = read_real_option("-mizo_final_hours", hours(horizons))
        call read_text_option("-mizo_profile_prefix", prefix, "")
        call read_text_option("-mizo_step_history", history, "")

        ierr = 0
        mesh_cells = 0
        call build_mizo_column_mesh(dm, mesh_cells, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the experiment mesh could not be built")
            failures = failures + 1
            return
        end if

        call configure_mizo_column(column)

        call column%initialize(dm, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the experiment column could not be initialised")
            failures = failures + 1
            call DMDestroy(dm, ierr)
            return
        end if

        ! The Jacobian asks for the same (slot, T, p_w) tangent once per
        ! registered block, and that map is a function of the key alone, so a
        ! memo returns the number the first call computed. -mizo_tangent_cache 0
        ! runs the path that differences every one of them again.
        call configure_frozen_tangent_cache( &
            read_integer_option("-mizo_tangent_cache", 1) /= 0, &
            read_integer_option("-mizo_tangent_cache_entries", 0))

        ! Uniform, as Conditions.json states it. The column is not in
        ! hydrostatic balance at that state and will redistribute under gravity;
        ! that redistribution is part of the reference problem, not an artefact
        ! to be initialised away.
        call column%set_initial_state(MIZO_INITIAL_TEMPERATURE, MIZO_INITIAL_PRESSURE, ierr)
        if (ierr == 0) call column%compute_report(start, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the initial state could not be published")
            failures = failures + 1
            call finish_frozen_column(column, dm)
            return
        end if
        call report_column("0 h     ", start)
        call check_integer("cells the column carries          ", start%owned_cells, mesh_cells, failures)
        call check_integer("frozen cells at the start         ", start%frozen_cells, 0, failures)
        call write_mizo_profile(column, prefix, 0.0d0, failures)

        call set_mizo_step_control(column)

        if (len_trim(history) > 0) then
            call column%set_step_history(trim(history), ierr)
            if (ierr /= 0) then
                call report_line("  FAIL the step history could not be opened: "//trim(history))
                failures = failures + 1
                call finish_frozen_column(column, dm)
                return
            end if
        end if

        elapsed = 0.0d0
        done = 0
        do h = 1, horizons
            if (hours(h) > final_hours + 1.0d-9) exit
            span = 3600.0d0*hours(h) - elapsed
            call column%advance(span, progress, ierr)
            elapsed = elapsed + progress%reached
            call report_horizon(hours(h), elapsed, progress)
            if (ierr /= 0) then
                call report_line("  FAIL the run stopped short of a horizon")
                failures = failures + 1
                exit
            end if
            call column%compute_report(reached(h), ierr)
            if (ierr /= 0) then
                call report_line("  FAIL a horizon report could not be measured")
                failures = failures + 1
                exit
            end if
            call report_column(mizo_label(hours(h)), reached(h))
            call write_mizo_profile(column, prefix, elapsed, failures)
            done = h
        end do

        if (done >= 1) call check_mizo_physics(start, reached(1:done), done, failures)

        call report_frozen_work()

        call finish_frozen_column(column, dm)
    end subroutine case_mizo_column

    !> How much constitutive work this rank asked for and how much it did.
    subroutine report_frozen_work()
        implicit none

        integer(int64) :: requests, evaluations, asked, responses, publications
        real(real64) :: tangent_saved, response_saved
        character(len=200) :: line

        requests = 0_int64
        evaluations = 0_int64
        asked = 0_int64
        responses = 0_int64
        publications = 0_int64
        call frozen_work_counters(requests, evaluations, asked, responses, publications)
        if (requests <= 0_int64 .or. asked <= 0_int64) return

        ! What each memo removed, as the share of the calls it answered.
        tangent_saved = real(requests - evaluations, real64)/real(requests, real64)
        response_saved = real(asked - responses, real64)/real(asked, real64)

        write (line, '(a,i0,a,i0,a,f6.2,a,i0,a,i0,a,f6.2,a,i0)') &
            "  work  tangents asked ", requests, " done ", evaluations, &
            " recalled ", 100.0d0*tangent_saved, " %  responses asked ", asked, &
            " done ", responses, " recalled ", 100.0d0*response_saved, &
            " %  phase splits ", publications
        call report_line(trim(line))
    end subroutine report_frozen_work

    !> What a run of the experiment can establish about itself.
    !>
    !> Not agreement with the measurement: these are the properties a coupled
    !> freezing run has to have before its profiles are worth comparing at all.
    subroutine check_mizo_physics(start, reached, done, failures)
        implicit none
        type(type_mixed_column_report), intent(in) :: start
        type(type_mixed_column_report), intent(in) :: reached(:)
        integer(int32), intent(in) :: done
        integer(int32), intent(inout) :: failures

        type(type_mixed_column_report) :: last
        integer(int32) :: h

        last = reached(done)

        ! The column is closed, so every drop it ends with is one it started
        ! with. What the bound allows is the residual each accepted step stops
        ! at, accumulated over the steps a run of this length takes.
        do h = 1, done
            call check("water the run created, relative    ", &
                       (reached(h)%total_water - start%total_water)/start%total_water, &
                       1.0d-7, failures)
        end do
        call check("water crossing the boundary, m2/s  ", last%boundary_rate, 1.0d-18, failures)
        call check("the domain balance, relative       ", &
                   last%domain_defect/start%total_water, 1.0d-7, failures)

        ! Freezing happened, and the cold end is where the experiment puts it.
        call check_exceeds("cells carrying ice                 ", real(last%frozen_cells, real64), 1.0d0, failures)
        call check_exceeds("ice the column made, m2            ", last%total_ice, 1.0d-6, failures)
        call check_logical("the cold end really is cold      ", last%coldest < -1.0d0, .true., failures)
        call check_logical("the warm end never froze         ", last%warmest > 0.0d0, .true., failures)

        ! Cryosuction, and water moving up the column toward it.
        call check_exceeds("suction the freezing gained, Pa    ", &
                           MIZO_INITIAL_PRESSURE - last%lowest_liquid_pressure, 1.0d4, failures)
        call check_exceeds("upward flux where the ice is, m/s  ", last%upward_flux_frozen, 1.0d-16, failures)

        ! The Jacobian every step was taken with was not singular.
        call check_exceeds("weakest pressure column            ", last%weakest_column, 1.0d-6, failures)
    end subroutine check_mizo_physics

    !> The experiment's material, gravity, scales and four boundaries.
    !>
    !> Shared by the experiment itself and by the profile export check, so that
    !> the smaller grid the latter runs on is the same physical problem and not
    !> a copy of it that can drift.
    subroutine configure_mizo_column(column)
        implicit none
        type(type_mixed_column), intent(inout) :: column

        type(type_mixed_frozen_material) :: frozen(1)
        real(real64) :: conductivity(MIXED_DIMENSION, MIXED_DIMENSION, 1)
        character(len=32) :: top_condition, base_condition, hydraulic_space
        real(real64) :: frozen_contact

        ! The hydraulic discretisation: the mixed RT0/P0 pair, or the primal
        ! H1 pressure that compares against it with everything else shared
        ! (`primal`: chain rule at the quadrature points; `primal_nodal`: nodal
        ! liquid pressure and geometric-mean conductivity, as bin/FTCMS).
        call read_text_option("-mizo_hydraulic_space", hydraulic_space, "mixed")
        call column%set_hydraulic_space(trim(hydraulic_space))

        call set_frozen_soil(frozen(1), FROZEN_MATERIAL_ID)
        ! The ice impedance is one of the two quantities the experiment is
        ! calibrated by, so a run may set it; the other tests keep the project's.
        frozen(1)%conductivity_function%omega = read_real_option("-mizo_impedance_factor", &
                                                                 FROZEN_IMPEDANCE_FACTOR)
        conductivity = 0.0d0
        conductivity(1, 1, 1) = FROZEN_CONDUCTIVITY
        conductivity(2, 2, 1) = FROZEN_CONDUCTIVITY

        call column%set_materials(frozen, conductivity)
        call column%set_gravity([0.0d0, -STANDARD_GRAVITY], REFERENCE_DENSITY)
        call column%set_scales(FROZEN_CONDUCTIVITY, abs(MIZO_INITIAL_PRESSURE), MIZO_HEIGHT, &
                               1.0d0, 1.0d0/MIZO_HEIGHT**2)
        call column%add_closed_boundary(COLUMN_CLOSED_ID)
        call read_text_option("-mizo_top_bc", top_condition, "plate")
        select case (trim(top_condition))
        case ("plate")
            call column%add_temperature_boundary(MIZO_TOP_ID, MIZO_PLATE_TEMPERATURES(1))
            call column%set_temperature_schedule(MIZO_TOP_ID, 3600.0d0*MIZO_PLATE_HOURS, &
                                                 MIZO_PLATE_TEMPERATURES)
        case ("contact")
            ! The plate at its measured temperature, and the soil behind the
            ! contact resistance between them: `plate` is the limit of an
            ! infinite contact coefficient.
            call column%add_heat_transfer_boundary(MIZO_TOP_ID, &
                                                   read_real_option("-mizo_contact_coefficient", &
                                                                    MIZO_CONTACT_COEFFICIENT), &
                                                   MIZO_PLATE_TEMPERATURES(1))
            call column%set_ambient_schedule(MIZO_TOP_ID, 3600.0d0*MIZO_PLATE_HOURS, &
                                             MIZO_PLATE_TEMPERATURES)
            ! The contact may change once the soil at the plate has frozen; a
            ! coefficient given for the frozen surface switches it there.
            frozen_contact = read_real_option("-mizo_contact_frozen_coefficient", -1.0d0)
            if (frozen_contact >= 0.0d0) then
                call column%set_heat_transfer_switch(MIZO_TOP_ID, frozen_contact, &
                                                     read_real_option("-mizo_contact_switch_temperature", 0.0d0), &
                                                     read_real_option("-mizo_contact_switch_width", 0.5d0))
            end if
        case ("convective")
            call column%add_heat_transfer_boundary(MIZO_TOP_ID, &
                                                   read_real_option("-mizo_top_coefficient", &
                                                                    MIZO_TOP_COEFFICIENT), &
                                                   MIZO_COOLANT)
        case default
            error stop "The Mizoguchi top condition is plate, contact or convective"
        end select

        call read_text_option("-mizo_base_bc", base_condition, "leak")
        select case (trim(base_condition))
        case ("leak")
            call column%add_heat_transfer_boundary(MIZO_BASE_ID, MIZO_BASE_COEFFICIENT, MIZO_AMBIENT)
        case ("measured")
            call column%add_temperature_boundary(MIZO_BASE_ID, MIZO_BASE_TEMPERATURES(1))
            call column%set_temperature_schedule(MIZO_BASE_ID, 3600.0d0*MIZO_PLATE_HOURS, &
                                                 MIZO_BASE_TEMPERATURES)
        case default
            error stop "The Mizoguchi base condition is either leak or measured"
        end select
    end subroutine configure_mizo_column

    !> The controller the experiment is stepped by.
    !>
    !> The window and the growth factor are the project's own `Conditions.json`
    !> values by default, and `-mizo_iterations_low`, `-mizo_iterations_high`
    !> and `-mizo_grow` are there so that a run can be told how hard to push the
    !> step without the tolerances a step is accepted at moving: the controller
    !> decides only how long the next step is attempted, never what it converges
    !> to. A run given any of them is not comparable with one that was not, for
    !> the same reason two runs with different horizons are not: the step
    !> sequence differs, so the states at a shared time are different numbers.
    !>
    !> `-mizo_min_step` lowers the floor below the project's 1 s. A graded grid
    !> (`-mizo_grading`) puts rows of a few hundredths of a millimetre under the
    !> plate, whose temperature drops to -2.7 C at the first instant; the front
    !> then crosses several such rows within the first second, which one step of
    !> 1 s cannot follow.
    subroutine set_mizo_step_control(column)
        implicit none
        type(type_mixed_column), intent(inout) :: column

        type(type_mixed_column_step_control) :: control

        control%initial = read_real_option("-mizo_initial_step", MIZO_INITIAL_STEP)
        control%minimum = read_real_option("-mizo_min_step", MIZO_MINIMUM_STEP)
        if (.not. (control%minimum > 0.0d0)) error stop "Mizoguchi minimum step must be positive"
        control%maximum = read_real_option("-mizo_max_step", MIZO_MAXIMUM_STEP)
        control%iterations_low = read_integer_option("-mizo_iterations_low", 5)
        control%iterations_high = read_integer_option("-mizo_iterations_high", 16)
        control%grow = read_real_option("-mizo_grow", 1.25d0)
        control%shrink = 0.8d0
        control%retry = 0.5d0
        call column%set_step_control(control)
    end subroutine set_mizo_step_control

    !> Write one horizon's cell profile, when the run asked for profiles.
    subroutine write_mizo_profile(column, prefix, time, failures)
        implicit none
        type(type_mixed_column), intent(inout) :: column
        character(len=*), intent(in) :: prefix
        real(real64), intent(in) :: time
        integer(int32), intent(inout) :: failures

        character(len=320) :: path
        PetscErrorCode :: ierr

        if (len_trim(prefix) == 0) return
        write (path, '(a,a,i0,a)') trim(prefix), "-", nint(time/3600.0d0), "h.csv"
        ierr = 0
        call column%write_profile(trim(path), time, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL a cell profile could not be written: "//trim(path))
            failures = failures + 1
            return
        end if
        call report_line("  wrote "//trim(path))
    end subroutine write_mizo_profile

    !> How far one horizon got and what it cost.
    subroutine report_horizon(hours, elapsed, progress)
        implicit none
        real(real64), intent(in) :: hours, elapsed
        type(type_mixed_column_progress), intent(in) :: progress

        character(len=250) :: line

        write (line, '(a,f6.2,a,f10.4,a,i0,a,i0,a,es12.5,a,i0,a,i0)') &
            "  horizon ", hours, " h  elapsed ", elapsed/3600.0d0, &
            " h  accepted ", progress%accepted_steps, "  rejected ", progress%rejected_steps, &
            "  shortest ", progress%shortest_step, "  its max ", progress%max_iterations, &
            "  total ", progress%total_iterations
        call report_line(trim(line))
    end subroutine report_horizon

    !> The times the run stops at, reports at, and writes a profile at.
    !>
    !> The published 12/24/50 h of the experiment unless `-mizo_profile_hours`
    !> names others, as a comma separated list of hours.
    !>
    !> A horizon is not only an output time: the advance is clipped to land
    !> exactly on it, so two runs given different lists do not take the same
    !> steps and their states at a shared time are not the same number. A mesh
    !> study must therefore give every one of its meshes the same list, and a
    !> run whose list is not the default is not comparable with one whose list
    !> is. That is the whole reason this is an option and not a schedule the
    !> run picks for itself.
    subroutine read_mizo_horizons(hours, horizons)
        implicit none
        !> The horizon times, hours, ascending; written for `1:horizons`.
        real(real64), intent(inout) :: hours(MAX_MIZO_HORIZONS)
        !> How many of `hours` were set.
        integer(int32), intent(inout) :: horizons

        PetscReal :: values(MAX_MIZO_HORIZONS)
        PetscInt :: count
        PetscBool :: found
        PetscErrorCode :: local_ierr
        integer(int32) :: h

        hours = 0.0d0
        hours(1:MIZO_HORIZONS) = MIZO_HOURS
        horizons = MIZO_HORIZONS

        count = int(MAX_MIZO_HORIZONS, PETSC_INT_KIND)
        found = PETSC_FALSE
        values = 0.0_PETSC_REAL_KIND
        call PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                                      "-mizo_profile_hours", values, count, found, local_ierr)
        if (local_ierr /= 0 .or. (found .eqv. PETSC_FALSE) .or. count < 1) return

        ! A list that is not ascending would ask the advance to run backwards.
        do h = 1, int(count, int32)
            if (.not. (values(h) > 0.0_PETSC_REAL_KIND)) return
            if (h > 1) then
                if (.not. (values(h) > values(h - 1))) return
            end if
        end do

        horizons = int(count, int32)
        do h = 1, horizons
            hours(h) = real(values(h), real64)
        end do
    end subroutine read_mizo_horizons

    !> The label one horizon's report is printed under.
    function mizo_label(hours) result(label)
        implicit none
        real(real64), intent(in) :: hours
        character(len=8) :: label

        write (label, '(i0,a)') nint(hours), " h    "
    end function mizo_label

    !> The intrinsic x-z strip the experiment runs on, labelled and distributed.
    !>
    !> `-mizo_cell_type quadrilateral` leaves the rectangles of the grid whole
    !> instead of cutting each into two triangles. At the same `-mizo_nx` and
    !> `-mizo_nz` the two meshes share every vertex, so the temperature space is
    !> the same size on both and what the comparison isolates is the flux and
    !> pressure pair. Nothing else about the problem moves with the flag: the
    !> width, the material, the initial state and all four boundaries are the
    !> same on either shape.
    !>
    !> `-mizo_grading R` makes the base row \(R\) times as thick as the top row,
    !> growing geometrically between them (`calc_strip_heights`); the default 1
    !> is the uniform grid. The freezing fringe under the cooled top is under a
    !> millimetre thick while the front crosses the first centimetre, and a
    !> uniform grid fine enough for it there is far finer than the deeper,
    !> thicker fringe needs.
    !>
    !> `-mizo_mesh_file <path>` reads the mesh from a Gmsh file instead
    !> (`load_mizo_mesh_file`), so that bin/FTCMS and this solver can run on one
    !> mesh file. The file then fixes the grid, and giving any of the grid options
    !> with it is refused rather than silently ignored.
    subroutine build_mizo_column_mesh(dm, cells, ierr)
        implicit none
        DM, intent(inout) :: dm
        !> Cells of the whole mesh, on every rank.
        integer(int32), intent(inout) :: cells
        PetscErrorCode, intent(inout) :: ierr

        character(len=*), parameter :: GRID_KEYS(5) = [character(len=15) :: "-mizo_width", "-mizo_nx", &
                                                       "-mizo_nz", "-mizo_grading", "-mizo_cell_type"]
        character(len=PETSC_MAX_PATH_LEN) :: mesh_file
        real(real64) :: width, grading
        real(real64), allocatable :: heights(:)
        integer(int32) :: nx, nz, k
        PetscBool :: given

        call read_text_option("-mizo_mesh_file", mesh_file, "")
        if (len_trim(mesh_file) > 0) then
            do k = 1, size(GRID_KEYS)
                call PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, trim(GRID_KEYS(k)), given, ierr)
                if (ierr /= 0) return
                if (given .eqv. PETSC_TRUE) then
                    call report_line("  "//trim(GRID_KEYS(k))//" cannot be combined with -mizo_mesh_file")
                    ierr = 1
                    return
                end if
            end do
            call load_mizo_mesh_file(trim(mesh_file), dm, cells, ierr)
            return
        end if

        width = read_real_option("-mizo_width", MIZO_WIDTH)
        if (width <= 0.0d0) error stop "Mizoguchi width must be positive"
        nx = read_integer_option("-mizo_nx", MIZO_NX)
        if (nx <= 0) error stop "Mizoguchi lateral cell count must be positive"
        nz = read_integer_option("-mizo_nz", MIZO_NZ)
        if (nz <= 0) error stop "Mizoguchi vertical cell count must be positive"
        grading = read_real_option("-mizo_grading", 1.0d0)
        if (.not. (grading > 0.0d0)) error stop "Mizoguchi grading ratio must be positive"

        allocate (heights(0:nz))
        call calc_strip_heights(nz, MIZO_HEIGHT, grading, heights)
        call report_value("top row thickness, m               ", heights(nz) - heights(nz - 1))
        call report_value("base row thickness, m              ", heights(1) - heights(0))
        deallocate (heights)

        cells = strip_cell_count(read_cell_type_option("-mizo_cell_type"), nx, nz)
        call build_column_strip(read_cell_type_option("-mizo_cell_type"), nx, nz, width, dm, ierr, grading)
    end subroutine build_mizo_column_mesh

    !> The experiment's strip at a given resolution, labelled and distributed.
    !>
    !> Split out so that a case which needs the same geometry at its own size --
    !> the profile export check, which wants a grid small enough to run at eight
    !> ranks in a test -- gets the experiment's labels and not a copy of them
    !> that could drift away from these.
    subroutine build_column_strip(quadrilateral, nx, nz, width, dm, ierr, grading)
        implicit none
        !> True for rectangles, false for the triangles the strip is cut into.
        logical, intent(in) :: quadrilateral
        integer(int32), intent(in) :: nx
        integer(int32), intent(in) :: nz
        real(real64), intent(in) :: width
        DM, intent(inout) :: dm
        PetscErrorCode, intent(inout) :: ierr
        !> Ratio of the base row's thickness to the top row's; absent means uniform.
        real(real64), intent(in), optional :: grading

        real(real64) :: reach

        reach = 1.0d3*max(width, MIZO_HEIGHT)

        call build_strip_of_type(quadrilateral, nx, nz, width, MIZO_HEIGHT, dm, ierr, grading)
        if (ierr /= 0) return
        call label_cells(dm, 0.5d0*width, FROZEN_MATERIAL_ID, FROZEN_MATERIAL_ID, ierr)
        if (ierr /= 0) return
        ! Every boundary facet is closed to water, which is the experiment's
        ! hydraulic condition on all four edges.
        call label_facets(dm, COLUMN_CLOSED_ID, -reach, reach, -reach, reach, ierr)
        if (ierr /= 0) return
        ! Only the two ends exchange heat; the sides carry no thermal label at
        ! all, which is the weak form's own zero-flux condition.
        call label_facets_named(dm, "Thermal Sets", MIZO_TOP_ID, -reach, reach, &
                                MIZO_HEIGHT - 1.0d-9, reach, ierr)
        if (ierr /= 0) return
        call label_facets_named(dm, "Thermal Sets", MIZO_BASE_ID, -reach, reach, &
                                -reach, 1.0d-9, ierr)
        if (ierr /= 0) return
        call distribute(dm, ierr)
    end subroutine build_column_strip

    !> The experiment's column read from a Gmsh file, labelled and distributed.
    !>
    !> The file is one bin/FTCMS reads: nodes at \((x,0,z)\), the soil as the
    !> physical surface `FROZEN_MATERIAL_ID`, and the edges as physical lines --
    !> `MIZO_TOP_ID` on the cooled top at \(z = 0.2\) m, `MIZO_BASE_ID` on the base
    !> at \(z = 0\) and any other tag on the sides. PETSc reads the groups as
    !> "Cell Sets" and "Face Sets".
    !>
    !> The solver works in the intrinsic \((x,z)\) plane, so the cells are rebuilt
    !> from the file's own vertex lists with \(y\) dropped (`intrinsic_from_physical`),
    !> each listed counter-clockwise, and the labels are carried over facet by
    !> facet: the top and base tags become "Thermal Sets", every boundary facet is
    !> closed to water and the sides carry no thermal label, exactly as on the
    !> built strip. The rebuilt mesh has the file's cells, vertices and facets one
    !> for one; each carried facet is matched by its two end vertices and its
    !> centroid is checked against the file's.
    !>
    !> Cost \(O(N)\) in the cells and facets of the file, on the reading rank.
    !>
    !> A file that is not that column is refused with `ierr = 1`: cells off the
    !> \(y = 0\) plane, a height other than `MIZO_HEIGHT`, mixed or degenerate cells,
    !> a second material, an unlabelled boundary facet, or a top or base group
    !> off its edge.
    subroutine load_mizo_mesh_file(path, dm, cells, ierr)
        implicit none
        !> The Gmsh file.
        character(len=*), intent(in) :: path
        !> The intrinsic mesh, created here.
        DM, intent(inout) :: dm
        !> Cells of the whole mesh, on every rank.
        integer(int32), intent(inout) :: cells
        PetscErrorCode, intent(inout) :: ierr

        real(real64), parameter :: TOLERANCE = 1.0d-9*MIZO_HEIGHT

        DM :: planar
        Vec :: coordinates
        PetscSection :: coordinate_section
        PetscScalar, pointer :: xyz(:)
        PetscInt, pointer :: closure(:), cone(:), support(:), join(:)
        PetscInt, allocatable :: cell_vertices(:)
        PetscReal, allocatable :: plane(:)
        PetscInt :: c_start, c_end, v_start, v_end, f_start, f_end, new_v_start, new_v_end
        PetscInt :: cell, vertex, facet, offset, closure_size, join_size, value, corners
        PetscInt :: dimension, coordinate_dimension, ends(2), listed(4)
        PetscReal :: measure, centroid(3), new_centroid(3), normal(3)
        PetscBool :: has_cells, has_faces
        real(real64) :: point(MIXED_DIMENSION), area, x_low, x_high, z_low, z_high, thinnest, thickest
        real(real64) :: cell_low, cell_high
        integer(int32) :: num_cells, num_vertices, found, k, status, marked(2), mpi_err
        PetscErrorCode :: cleanup_ierr
        character(len=200) :: line

        ! PETSc keeps the first two coordinates of a two-dimensional Gmsh mesh
        ! unless told otherwise, which for an x-z file is x and a constant y.
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-dm_plex_gmsh_spacedim", "3", ierr)
        if (ierr /= 0) return
        call DMPlexCreateFromFile(PETSC_COMM_WORLD, path, "mizo_mesh", PETSC_TRUE, planar, ierr)
        call PetscOptionsClearValue(PETSC_NULL_OPTIONS, "-dm_plex_gmsh_spacedim", cleanup_ierr)
        if (ierr /= 0) then
            call report_line("  the mesh file could not be read: "//path)
            return
        end if

        status = 0
        call DMGetDimension(planar, dimension, ierr)
        if (ierr == 0) call DMGetCoordinateDim(planar, coordinate_dimension, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(planar, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr == 0) call DMPlexGetDepthStratum(planar, 0_PETSC_INT_KIND, v_start, v_end, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(planar, 1_PETSC_INT_KIND, f_start, f_end, ierr)
        if (ierr == 0) call DMHasLabel(planar, "Cell Sets", has_cells, ierr)
        if (ierr == 0) call DMHasLabel(planar, "Face Sets", has_faces, ierr)
        if (ierr /= 0) then
            call DMDestroy(planar, cleanup_ierr)
            return
        end if
        num_cells = int(c_end - c_start, int32)
        num_vertices = int(v_end - v_start, int32)
        if (dimension /= 2 .or. coordinate_dimension /= 3) then
            call report_line("  the mesh file is not a surface mesh in three coordinates")
            status = 1
        else if (num_cells > 0 .and. ((has_cells .neqv. PETSC_TRUE) .or. (has_faces .neqv. PETSC_TRUE))) then
            call report_line("  the mesh file carries no physical groups")
            status = 1
        end if

        ! The vertices in the intrinsic plane, in the file's order.
        allocate (plane(max(1, 2*num_vertices)))
        x_low = huge(1.0d0)
        x_high = -huge(1.0d0)
        z_low = huge(1.0d0)
        z_high = -huge(1.0d0)
        if (status == 0 .and. num_vertices > 0) then
            call DMGetCoordinatesLocal(planar, coordinates, ierr)
            if (ierr == 0) call DMGetCoordinateSection(planar, coordinate_section, ierr)
            if (ierr == 0) call VecGetArrayRead(coordinates, xyz, ierr)
            if (ierr /= 0) then
                call DMDestroy(planar, cleanup_ierr)
                return
            end if
            do vertex = v_start, v_end - 1
                call PetscSectionGetOffset(coordinate_section, vertex, offset, ierr)
                if (ierr /= 0) exit
                if (abs(real(xyz(offset + 2), real64)) > TOLERANCE) status = 1
                point = intrinsic_from_physical(real(xyz(offset + 1:offset + 3), real64))
                plane(2*(vertex - v_start) + 1:2*(vertex - v_start) + 2) = real(point, PETSC_REAL_KIND)
                x_low = min(x_low, point(1))
                x_high = max(x_high, point(1))
                z_low = min(z_low, point(2))
                z_high = max(z_high, point(2))
            end do
            call VecRestoreArrayRead(coordinates, xyz, cleanup_ierr)
            if (ierr == 0) ierr = cleanup_ierr
            if (ierr /= 0) then
                call DMDestroy(planar, cleanup_ierr)
                return
            end if
            if (status /= 0) call report_line("  the mesh file has nodes off the y = 0 plane")
            if (status == 0 .and. (abs(z_low) > TOLERANCE .or. abs(z_high - MIZO_HEIGHT) > TOLERANCE)) then
                call report_line("  the mesh file does not span z = 0 to the column height")
                status = 1
            end if
        end if

        ! Each cell's vertices in the order its closure walks them, which is the
        ! cyclic order round the cell; reversed where that is clockwise in (x, z).
        corners = 0
        allocate (cell_vertices(max(1, 4*num_cells)))
        thinnest = huge(1.0d0)
        thickest = 0.0d0
        do cell = c_start, c_end - 1
            if (status /= 0) exit
            nullify (closure)
            call DMPlexGetTransitiveClosure(planar, cell, PETSC_TRUE, closure_size, closure, ierr)
            if (ierr /= 0) exit
            found = 0
            do k = 1, int(closure_size, int32)
                if (closure(2*k - 1) < v_start .or. closure(2*k - 1) >= v_end) cycle
                found = found + 1
                if (found <= 4) listed(found) = closure(2*k - 1) - v_start
            end do
            call DMPlexRestoreTransitiveClosure(planar, cell, PETSC_TRUE, closure_size, closure, ierr)
            if (ierr /= 0) exit
            if (corners == 0 .and. (found == 3 .or. found == 4)) corners = found
            if (found /= corners) then
                call report_line("  the mesh file mixes cell shapes or has cells of neither three nor four vertices")
                status = 1
                exit
            end if
            area = 0.0d0
            cell_low = huge(1.0d0)
            cell_high = -huge(1.0d0)
            do k = 1, found
                area = area + plane(2*listed(k) + 1)*plane(2*listed(mod(k, found) + 1) + 2) &
                       - plane(2*listed(mod(k, found) + 1) + 1)*plane(2*listed(k) + 2)
                cell_low = min(cell_low, real(plane(2*listed(k) + 2), real64))
                cell_high = max(cell_high, real(plane(2*listed(k) + 2), real64))
            end do
            if (area < 0.0d0) listed(2:found) = listed(found:2:-1)
            if (area == 0.0d0) then
                call report_line("  the mesh file has a cell of no area")
                status = 1
                exit
            end if
            thinnest = min(thinnest, cell_high - cell_low)
            thickest = max(thickest, cell_high - cell_low)
            cell_vertices(corners*(cell - c_start) + 1:corners*(cell - c_start + 1)) = listed(1:found)

            call DMGetLabelValue(planar, "Cell Sets", cell, value, ierr)
            if (ierr /= 0) exit
            if (value /= FROZEN_MATERIAL_ID) then
                call report_line("  the mesh file has a cell outside the soil's physical surface")
                status = 1
                exit
            end if
        end do
        if (ierr /= 0) then
            call DMDestroy(planar, cleanup_ierr)
            return
        end if

        ! Only the reading rank has seen the file; every rank returns together.
        call MPI_Bcast(status, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, mpi_err)
        if (status /= 0) then
            call DMDestroy(planar, cleanup_ierr)
            ierr = 1
            return
        end if

        ! Only the first rank's lists are read; the others hold none.
        call DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, 2_PETSC_INT_KIND, int(num_cells, PETSC_INT_KIND), &
                                           int(num_vertices, PETSC_INT_KIND), corners, PETSC_TRUE, &
                                           cell_vertices, 2_PETSC_INT_KIND, plane, dm, ierr)
        deallocate (cell_vertices, plane)
        if (ierr /= 0) then
            call DMDestroy(planar, cleanup_ierr)
            return
        end if

        call label_cells(dm, 0.5d0*(x_low + x_high), FROZEN_MATERIAL_ID, FROZEN_MATERIAL_ID, ierr)
        if (ierr == 0) call label_facets(dm, COLUMN_CLOSED_ID, -huge(1.0d0), huge(1.0d0), &
                                         -huge(1.0d0), huge(1.0d0), ierr)
        if (ierr == 0) call DMCreateLabel(dm, "Thermal Sets", ierr)
        if (ierr == 0) call DMPlexGetDepthStratum(dm, 0_PETSC_INT_KIND, new_v_start, new_v_end, ierr)
        if (ierr /= 0) then
            call DMDestroy(planar, cleanup_ierr)
            return
        end if
        if (new_v_end - new_v_start /= v_end - v_start) then
            call report_line("  the rebuilt mesh does not have the file's vertices")
            status = 1
        end if

        ! The file's boundary facets, carried to the facet joining the same two
        ! vertices of the rebuilt mesh.
        marked = 0
        do facet = f_start, f_end - 1
            if (status /= 0) exit
            call DMPlexGetSupport(planar, facet, support, ierr)
            if (ierr /= 0) exit
            if (size(support) /= 1) cycle
            call DMGetLabelValue(planar, "Face Sets", facet, value, ierr)
            if (ierr /= 0) exit
            if (value == -1_PETSC_INT_KIND) then
                call report_line("  the mesh file has a boundary facet in no physical group")
                status = 1
                exit
            end if
            if (value /= MIZO_TOP_ID .and. value /= MIZO_BASE_ID) cycle

            call DMPlexGetCone(planar, facet, cone, ierr)
            if (ierr /= 0) exit
            ends = cone(1:2) - v_start + new_v_start
            call DMPlexRestoreCone(planar, facet, cone, ierr)
            if (ierr /= 0) exit
            call DMPlexComputeCellGeometryFVM(planar, facet, measure, centroid, normal, ierr)
            if (ierr /= 0) exit
            call DMPlexGetJoin(dm, 2_PETSC_INT_KIND, ends, join_size, join, ierr)
            if (ierr /= 0) exit
            if (join_size /= 1) then
                status = 1
            else
                call DMPlexComputeCellGeometryFVM(dm, join(1), measure, new_centroid, normal, ierr)
                if (ierr == 0 .and. (abs(new_centroid(1) - centroid(1)) > TOLERANCE .or. &
                                     abs(new_centroid(2) - centroid(3)) > TOLERANCE)) status = 1
                ! The cooled top and the base are the column's two ends.
                if (value == MIZO_TOP_ID .and. abs(centroid(3) - z_high) > TOLERANCE) status = 1
                if (value == MIZO_BASE_ID .and. abs(centroid(3) - z_low) > TOLERANCE) status = 1
                if (ierr == 0 .and. status == 0) call DMSetLabelValue(dm, "Thermal Sets", join(1), value, ierr)
            end if
            call DMPlexRestoreJoin(dm, 2_PETSC_INT_KIND, ends, join_size, join, cleanup_ierr)
            if (ierr /= 0) exit
            if (status /= 0) then
                call report_line("  the mesh file's top or base group could not be carried to its edge")
                exit
            end if
            if (value == MIZO_TOP_ID) marked(1) = marked(1) + 1
            if (value == MIZO_BASE_ID) marked(2) = marked(2) + 1
        end do
        call DMDestroy(planar, cleanup_ierr)
        if (ierr /= 0) return
        if (num_cells > 0 .and. status == 0 .and. (marked(1) == 0 .or. marked(2) == 0)) then
            call report_line("  the mesh file has no top or no base group")
            status = 1
        end if
        call MPI_Bcast(status, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, mpi_err)
        if (status /= 0) then
            ierr = 1
            return
        end if

        call MPI_Bcast(num_cells, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, mpi_err)
        call MPI_Bcast(marked, 2, MPI_INTEGER, 0, PETSC_COMM_WORLD, mpi_err)
        call MPI_Bcast(x_low, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, mpi_err)
        call MPI_Bcast(x_high, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, mpi_err)
        call MPI_Bcast(thinnest, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, mpi_err)
        call MPI_Bcast(thickest, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, mpi_err)
        cells = num_cells
        call report_line("  mesh file "//path)
        write (line, '(a,i0,a,i0,a,i0,a,i0,a)') "  note the file holds ", num_cells, " cells of ", corners, &
            " vertices; ", marked(1), " top and ", marked(2), " base facets"
        call report_line(trim(line))
        call report_value("mesh width, m                      ", x_high - x_low)
        call report_value("thinnest cell height, m            ", thinnest)
        call report_value("thickest cell height, m            ", thickest)

        call distribute(dm, ierr)
    end subroutine load_mizo_mesh_file

    !> What the cell profile export writes, at whatever rank count it is run on.
    !>
    !> The export is a gather, so the question it has to answer is whether the
    !> file is the domain. Three things are checked on the file itself and are
    !> rank-count independent by construction: that it holds one row per cell of
    !> the whole strip, that the cell numbers are \([0,N)\) exactly once -- no
    !> ghost written twice, no rank's share missing -- and that the volumes sum
    !> to the strip's area, which a file missing rows cannot do. The fourth ties
    !> the file to the solver: its volume-weighted water has to equal the
    !> centroid water of the reduced audit, which is computed on a different
    !> path and reduced over all the ranks.
    !>
    !> Then the same run at another rank count is compared against a reference
    !> written by an earlier one, given by `-export_reference`. Both the initial
    !> state and a solved state are compared, because a gather can be right on a
    !> uniform field and wrong on one that varies from rank to rank; the solved
    !> state is asserted to have left the initial one, or the comparison would
    !> be of two copies of the same uniform column.
    subroutine case_profile_export(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_column) :: column
        type(type_mixed_column_report) :: start, solved
        type(type_mixed_column_progress) :: progress
        type(type_mixed_column_step_control) :: control
        character(len=256) :: prefix, reference
        character(len=320) :: initial_path, solved_path
        DM :: dm
        PetscErrorCode :: ierr
        PetscMPIInt :: comm_size
        integer(int32) :: nx, nz, cells
        real(real64) :: width

        call report_line("case profile_export: the cell profile is the domain, at any rank count")

        nx = read_integer_option("-export_nx", EXPORT_NX)
        nz = read_integer_option("-export_nz", EXPORT_NZ)
        width = read_real_option("-mizo_width", MIZO_WIDTH)
        cells = nx*nz
        call read_text_option("-export_prefix", prefix, "export")
        call read_text_option("-export_reference", reference, "")
        write (initial_path, '(a,a)') trim(prefix), "-initial.csv"
        write (solved_path, '(a,a)') trim(prefix), "-solved.csv"

        call MPI_Comm_size(PETSC_COMM_WORLD, comm_size, ierr)
        call report_integer("ranks this run is on              ", int(comm_size, int32))

        ierr = 0
        call build_column_strip(.true., nx, nz, width, dm, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the export strip could not be built")
            failures = failures + 1
            return
        end if

        call configure_mizo_column(column)
        call column%initialize(dm, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the export column could not be initialised")
            failures = failures + 1
            call DMDestroy(dm, ierr)
            return
        end if

        call column%set_initial_state(MIZO_INITIAL_TEMPERATURE, MIZO_INITIAL_PRESSURE, ierr)
        if (ierr == 0) call column%compute_report(start, ierr)
        if (ierr == 0) call column%write_profile(trim(initial_path), 0.0d0, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the initial profile could not be written")
            failures = failures + 1
            call finish_frozen_column(column, dm)
            return
        end if
        call check_integer("cells the strip carries           ", start%owned_cells, cells, failures)
        call check_exported_profile(trim(initial_path), "initial", cells, width*MIZO_HEIGHT, &
                                    start%total_water_centroid, failures)

        control%initial = 75.0d0
        control%minimum = MIZO_MINIMUM_STEP
        control%maximum = 450.0d0
        control%iterations_low = 5
        control%iterations_high = 16
        control%grow = 1.25d0
        control%shrink = 0.8d0
        control%retry = 0.5d0
        call column%set_step_control(control)

        call column%advance(EXPORT_SPAN, progress, ierr)
        if (ierr == 0) call column%compute_report(solved, ierr)
        if (ierr == 0) call column%write_profile(trim(solved_path), progress%reached, ierr)
        if (ierr /= 0) then
            call report_line("  FAIL the solved profile could not be written")
            failures = failures + 1
            call finish_frozen_column(column, dm)
            return
        end if
        call check("the span the run had to cover, s   ", progress%reached - EXPORT_SPAN, 1.0d-9, failures)
        call check_exported_profile(trim(solved_path), "solved ", cells, width*MIZO_HEIGHT, &
                                    solved%total_water_centroid, failures)

        ! The solved state has to be a state, or comparing it across rank
        ! counts would only compare two copies of the uniform initial column.
        call check_exceeds("cooling the solved profile shows, C", &
                           MIZO_INITIAL_TEMPERATURE - solved%coldest, 1.0d0, failures)

        if (len_trim(reference) > 0) then
            call compare_against_reference(trim(initial_path), trim(reference)//"-initial.csv", &
                                           "initial", EXPORT_INITIAL_TOLERANCE, failures)
            call compare_against_reference(trim(solved_path), trim(reference)//"-solved.csv", &
                                           "solved ", EXPORT_SOLVED_TOLERANCE, failures)
        end if

        call finish_frozen_column(column, dm)
        call share_failures(failures)
    end subroutine case_profile_export

    !> Read one written profile back and ask whether it is the whole domain.
    subroutine check_exported_profile(path, label, cells, area, water, failures)
        implicit none
        !> The file just written.
        character(len=*), intent(in) :: path
        !> Which profile this is, for the printed lines.
        character(len=*), intent(in) :: label
        !> How many cells the strip has in all.
        integer(int32), intent(in) :: cells
        !> The strip's area, m2.
        real(real64), intent(in) :: area
        !> The reduced audit's centroid water, m2, for the same state.
        real(real64), intent(in) :: water
        integer(int32), intent(inout) :: failures

        integer(int32), allocatable :: ids(:), seen(:)
        real(real64), allocatable :: values(:, :)
        integer(int32) :: rows, k, duplicates, out_of_range
        real(real64) :: total_area, total_water
        PetscMPIInt :: rank
        PetscErrorCode :: local_ierr
        logical :: ok

        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, local_ierr)
        if (rank /= 0) return

        call read_profile_csv(path, ids, values, rows, ok)
        if (.not. ok) then
            call report_line("  FAIL the "//label//" profile could not be read back: "//path)
            failures = failures + 1
            return
        end if

        call check_integer("rows the "//label//" profile holds  ", rows, cells, failures)

        allocate (seen(max(cells, 1)))
        seen = 0
        duplicates = 0
        out_of_range = 0
        do k = 1, rows
            if (ids(k) < 0 .or. ids(k) >= cells) then
                out_of_range = out_of_range + 1
                cycle
            end if
            if (seen(ids(k) + 1) /= 0) duplicates = duplicates + 1
            seen(ids(k) + 1) = seen(ids(k) + 1) + 1
        end do
        call check_integer("cell numbers out of range         ", out_of_range, 0, failures)
        call check_integer("cells written more than once      ", duplicates, 0, failures)
        call check_integer("cells the "//label//" profile misses ", count(seen == 0), 0, failures)

        total_area = sum(values(PROFILE_VOLUME, 1:rows))
        total_water = sum(values(PROFILE_VOLUME, 1:rows)*values(PROFILE_STORED, 1:rows))
        call check("area the "//label//" profile covers  ", (total_area - area)/area, 1.0d-14, failures)
        call check("its water against the audit's      ", (total_water - water)/water, 1.0d-12, failures)

        deallocate (seen)
    end subroutine check_exported_profile

    !> Compare a profile against one written by a run at another rank count.
    !>
    !> Matched cell by cell on the global cell number, so the comparison does
    !> not depend on the two files listing their rows in the same order. Each
    !> column is judged against its own range, because the columns of a profile
    !> span everything from a coordinate of order \(10^{-1}\) to a relative
    !> conductivity of order \(10^{-14}\) and one absolute tolerance would be
    !> vacuous for some of them and impossible for others.
    subroutine compare_against_reference(path, reference_path, label, tolerance, failures)
        implicit none
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: reference_path
        character(len=*), intent(in) :: label
        !> Allowed difference as a fraction of each column's own range.
        real(real64), intent(in) :: tolerance
        integer(int32), intent(inout) :: failures

        integer(int32), allocatable :: ids(:), reference_ids(:), place(:)
        real(real64), allocatable :: values(:, :), reference_values(:, :)
        real(real64) :: ranges(PROFILE_COLUMNS), scales(PROFILE_COLUMNS)
        integer(int32) :: rows, reference_rows, k, c, worst_column
        real(real64) :: difference, worst
        character(len=32) :: column_label
        PetscMPIInt :: rank
        PetscErrorCode :: local_ierr
        logical :: ok

        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, local_ierr)
        if (rank /= 0) return

        call read_profile_csv(path, ids, values, rows, ok)
        if (ok) call read_profile_csv(reference_path, reference_ids, reference_values, reference_rows, ok)
        if (.not. ok) then
            call report_line("  FAIL the "//label//" comparison could not read "//reference_path)
            failures = failures + 1
            return
        end if
        call check_integer("rows the "//label//" reference has  ", reference_rows, rows, failures)
        if (reference_rows /= rows) return

        ! Where each cell number sits in the reference file, so the two files
        ! may list their rows in any order.
        allocate (place(0:rows - 1))
        place = -1
        do k = 1, reference_rows
            if (reference_ids(k) >= 0 .and. reference_ids(k) < rows) place(reference_ids(k)) = k
        end do

        ! Each scale group takes the largest range of any column in it.
        do c = 1, PROFILE_COLUMNS
            ranges(c) = maxval(abs(reference_values(c, 1:reference_rows)))
        end do
        do c = 1, PROFILE_COLUMNS
            scales(c) = maxval(ranges, mask=(PROFILE_SCALE_GROUP == PROFILE_SCALE_GROUP(c)))
            if (scales(c) <= 0.0d0) scales(c) = 1.0d0
        end do

        worst = 0.0d0
        worst_column = 1
        do c = 1, PROFILE_COLUMNS
            do k = 1, rows
                if (ids(k) < 0 .or. ids(k) >= rows) cycle
                if (place(ids(k)) < 1) cycle
                difference = abs(values(c, k) - reference_values(c, place(ids(k))))/scales(c)
                if (difference > worst) then
                    worst = difference
                    worst_column = c
                end if
            end do
        end do

        write (column_label, '(a,i0,a)') "(column ", worst_column, ")"
        call check("the "//label//" profile against rank 0 ", worst, tolerance, failures)
        call report_line("       worst column of the "//label//" comparison "//trim(column_label))

        deallocate (place)
    end subroutine compare_against_reference

    !> Read a written profile back into its cell numbers and its values.
    !>
    !> List-directed input, which reads the commas of the file as separators.
    !> The file is read twice: once to count its rows and once to take them.
    subroutine read_profile_csv(path, ids, values, rows, ok)
        implicit none
        character(len=*), intent(in) :: path
        !> The cell number of each row, allocated here.
        integer(int32), allocatable, intent(inout) :: ids(:)
        !> `PROFILE_COLUMNS` numbers per row, allocated here.
        real(real64), allocatable, intent(inout) :: values(:, :)
        !> How many rows the file holds.
        integer(int32), intent(inout) :: rows
        logical, intent(inout) :: ok

        character(len=512) :: line
        real(real64) :: time
        integer(int32) :: unit_number, io_status, k

        ok = .false.
        rows = 0
        open (newunit=unit_number, file=path, status="old", action="read", iostat=io_status)
        if (io_status /= 0) return
        read (unit_number, '(a)', iostat=io_status) line
        if (io_status /= 0) then
            close (unit_number)
            return
        end if
        do
            read (unit_number, '(a)', iostat=io_status) line
            if (io_status /= 0) exit
            if (len_trim(line) == 0) cycle
            rows = rows + 1
        end do

        rewind (unit_number)
        allocate (ids(max(rows, 1)), values(PROFILE_COLUMNS, max(rows, 1)))
        read (unit_number, '(a)', iostat=io_status) line
        do k = 1, rows
            read (unit_number, *, iostat=io_status) time, ids(k), values(:, k)
            if (io_status /= 0) then
                close (unit_number)
                return
            end if
        end do
        close (unit_number)
        ok = .true.
    end subroutine read_profile_csv

    !> Give every rank rank 0's verdict, so a check made on one rank's file
    !> stops the whole run and not one process of it.
    subroutine share_failures(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        integer(int32) :: reduced, mpi_err

        call MPI_Allreduce(failures, reduced, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, mpi_err)
        if (mpi_err == 0) failures = reduced
    end subroutine share_failures

    !> Record one integer that is reported and not judged.
    subroutine report_integer(label, value)
        implicit none
        character(len=*), intent(in) :: label
        integer(int32), intent(in) :: value

        character(len=160) :: line

        write (line, '(a,a,a,i0)') "  note ", label, "  ", value
        call report_line(trim(line))
    end subroutine report_integer

    !> One column report, printed.
    subroutine report_column(label, report)
        implicit none
        character(len=*), intent(in) :: label
        type(type_mixed_column_report), intent(in) :: report

        character(len=250) :: line

        write (line, '(a,a,a,f7.3,a,f7.3,a,es12.5,a,es12.5,a,i0)') &
            "  column ", label, " T ", report%coldest, " to ", report%warmest, &
            " C  water ", report%total_water, "  ice ", report%total_ice, &
            "  frozen cells ", report%frozen_cells
        call report_line(trim(line))
        write (line, '(a,a,a,es12.5,a,es12.5,a,es12.5,a,es12.5,a,es12.5)') &
            "  column ", label, " p_l min ", report%lowest_liquid_pressure, &
            "  q_up ", report%upward_flux, "  q_up frozen ", report%upward_flux_frozen, &
            "  q_down frozen ", report%downward_flux_frozen, "  margin ", report%weakest_column
        call report_line(trim(line))
        write (line, '(a,a,a,es12.5,a,es12.5,a,es12.5,a,es12.5)') &
            "  column ", label, " cell defect ", report%max_cell_defect, &
            "  domain defect ", report%domain_defect, "  boundary ", report%boundary_rate, &
            "  throughput ", report%throughput
        call report_line(trim(line))
        write (line, '(a,a,6(a,es12.5))') "  audit ", label, &
            " storage ", report%max_storage_change, " flux ", report%max_flux_term, &
            " centroid gap ", report%max_quadrature_gap, " facet jump ", report%max_facet_mismatch, &
            " K spread ", report%max_conductivity_spread, " p_l spread ", report%max_liquid_pressure_spread
        call report_line(trim(line))
    end subroutine report_column

    !> The intrinsic x-z strip the column runs on, labelled and distributed.
    !>
    !> `-column_cell_type quadrilateral` runs the same column on rectangles, on
    !> the same vertex grid; every other property of the fixture is unchanged.
    subroutine build_frozen_column_mesh(dm, ierr)
        implicit none
        DM, intent(inout) :: dm
        PetscErrorCode, intent(inout) :: ierr

        real(real64) :: reach

        reach = 1.0d3*max(COLUMN_WIDTH, COLUMN_HEIGHT)

        call build_strip_of_type(read_cell_type_option("-column_cell_type"), COLUMN_NX, COLUMN_NZ, &
                                 COLUMN_WIDTH, COLUMN_HEIGHT, dm, ierr)
        if (ierr /= 0) return
        ! One material everywhere: both sides of the split carry the same id.
        call label_cells(dm, 0.5d0*COLUMN_WIDTH, FROZEN_MATERIAL_ID, FROZEN_MATERIAL_ID, ierr)
        if (ierr /= 0) return
        ! Every boundary facet is closed to water. An unlabelled one would be a
        ! prescribed zero liquid pressure, which is an open boundary.
        call label_facets(dm, COLUMN_CLOSED_ID, -reach, reach, -reach, reach, ierr)
        if (ierr /= 0) return
        call label_facets_named(dm, "Thermal Sets", COLUMN_TOP_ID, -reach, reach, &
                                COLUMN_HEIGHT - 1.0d-9, reach, ierr)
        if (ierr /= 0) return
        call label_facets_named(dm, "Thermal Sets", COLUMN_BOTTOM_ID, -reach, reach, &
                                -reach, 1.0d-9, ierr)
        if (ierr /= 0) return
        call distribute(dm, ierr)
    end subroutine build_frozen_column_mesh

    !> Mark the boundary facets of a box into a label of the caller's choosing.
    !>
    !> The water boundary and the thermal boundary partition the perimeter
    !> differently -- every facet is closed to water, only the two ends carry a
    !> prescribed temperature -- so they cannot share one label.
    subroutine label_facets_named(dm, label_name, value, x_low, x_high, z_low, z_high, ierr)
        implicit none
        DM, intent(inout) :: dm
        character(len=*), intent(in) :: label_name
        integer(int32), intent(in) :: value
        real(real64), intent(in) :: x_low, x_high, z_low, z_high
        PetscErrorCode, intent(inout) :: ierr

        PetscInt :: facet, f_start, f_end
        PetscInt, pointer :: support(:)
        PetscReal :: measure, centroid(3), normal(3)
        real(real64) :: x, z

        call DMCreateLabel(dm, label_name, ierr)
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
            call DMSetLabelValue(dm, label_name, facet, int(value, PETSC_INT_KIND), ierr)
            if (ierr /= 0) return
        end do
    end subroutine label_facets_named

    !> Release the column and the mesh it was built on.
    subroutine finish_frozen_column(column, dm)
        implicit none
        type(type_mixed_column), intent(inout) :: column
        DM, intent(inout) :: dm

        PetscErrorCode :: ierr

        ierr = 0
        call column%destroy(ierr)
        call DMDestroy(dm, ierr)
    end subroutine finish_frozen_column

    !> The one-sided difference where the closure changes branch.
    !>
    !> The compressive branch \(S_s/\rho_{\mathrm{std}}g\max(0,p_w-p_{\mathrm{sat}})\)
    !> opens only once the phase model calls the cell saturated, and for this
    !> retention curve that happens within a pascal of \(p_w=0\): the curve
    !> reaches \(\theta_s\) only at zero head. So at \(p_w=0\) the closure's
    !> pressure derivative is genuinely two-valued -- the retention slope from
    !> below, the compressive capacity from above -- and a central difference
    !> across it is a tangent of neither branch. What is measured is that the
    !> contract takes the difference on the side the point itself is on, says
    !> so, and that the two sides really do differ; without the last check the
    !> assertion would pass on a closure with no branch at all.
    subroutine check_frozen_branch(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_frozen_tangents) :: tangents
        type(type_mixed_frozen_response) :: base, ahead, behind
        !> The offset `eval_frozen_tangents` differences over, so that the
        !> one-sided value asserted here is the one it can have formed.
        real(real64), parameter :: STEP = 1.0d0
        real(real64) :: forward, backward, central
        logical :: ok

        ok = .false.
        call eval_frozen_response(1_int32, FROZEN_WARM, 0.0d0, base, ok)
        if (ok) call eval_frozen_response(1_int32, FROZEN_WARM, STEP, ahead, ok)
        if (ok) call eval_frozen_response(1_int32, FROZEN_WARM, -STEP, behind, ok)
        if (ok) call eval_frozen_tangents(1_int32, FROZEN_WARM, 0.0d0, tangents, ok)
        if (.not. ok) then
            call report_line("  FAIL the saturation-pressure branch could not be evaluated")
            failures = failures + 1
            return
        end if

        ! The branch really is a branch, and which side the point is on is what
        ! decides the difference: saturated at and above zero, unsaturated below.
        call check_logical("compressive branch open at p_w=0  ", base%compressive_active, .true., failures)
        call check_logical("compressive branch open above    ", ahead%compressive_active, .true., failures)
        call check_logical("compressive branch shut below    ", behind%compressive_active, .false., failures)
        call check_logical("dW/dp at p_sat taken one-sided   ", tangents%pressure_one_sided, .true., failures)
        call check_logical("no branch straddled at p_sat     ", tangents%branch_straddled, .false., failures)

        forward = (ahead%stored_water - base%stored_water)/STEP
        backward = (base%stored_water - behind%stored_water)/STEP
        central = (ahead%stored_water - behind%stored_water)/(2.0d0*STEP)
        call check("dW/dp is the branch's own slope   ", &
                   (tangents%dwater_dpressure - forward)/max(abs(forward), tiny(1.0d0)), 1.0d-9, failures)
        ! The two sides differ by more than the compressive capacity itself, so
        ! the central difference would not have been a small error.
        call check_exceeds("the two sides of p_sat differ    ", &
                           abs(backward - forward)/max(abs(forward), tiny(1.0d0)), 1.0d1, failures)
        call check_exceeds("central difference error avoided ", &
                           abs(central - forward)/(0.5d0*abs(backward - forward)), 0.99d0, failures)
    end subroutine check_frozen_branch

    !> Where the pressure column of the coupled Jacobian collapses.
    !>
    !> \(p_w\) reaches the monolithic system only through \(p_l\) and \(W\). On
    !> the freezing branch \(\partial p_l/\partial p_w=0\) identically, because
    !> \(\psi_{\mathrm{eff}}=\psi_{\mathrm{cryo}}(T)\) carries no pressure, so
    !> the storage capacity is the whole column there. The sweep covers the
    !> states an unsaturated Mizoguchi column visits, where the retention slope
    !> keeps the column healthy, and then the saturated frozen state, where it
    !> is not small but exactly zero: the retention curve is flat above
    !> saturation and the phase model stops calling an ice-bearing cell
    !> saturated, so the compressive branch is shut too and no quantity in the
    !> closure responds to \(p_w\) at all. That state is recorded as the
    !> supported envelope's boundary, not floored: a capacity invented to make
    !> the matrix invertible would store water no phase holds.
    subroutine check_frozen_rank(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        !> The pressure change and the water content the column is measured
        !> against: a 10 kPa suction step and a milli-unit of water content are
        !> what a Mizoguchi step resolves.
        real(real64), parameter :: PRESSURE_REFERENCE = 1.0d4
        real(real64), parameter :: STORAGE_REFERENCE = 1.0d-3
        !> Below this the column contributes less than a per-mille of a
        !> reference response and the cell is reported degenerate.
        real(real64), parameter :: MARGIN_FLOOR = 1.0d-3

        real(real64) :: temperatures(4), pressures(2)
        real(real64) :: margin, sensitivity, capacity, weakest
        character(len=200) :: line
        integer(int32) :: t, p
        logical :: ok

        temperatures = [FROZEN_WARM, -0.1d0, FROZEN_MILD_COLD, FROZEN_DEEP_COLD]
        ! Unsaturated only: the saturated state is the separate measurement below.
        pressures = [FROZEN_FIELD_PRESSURE, FROZEN_WET_PRESSURE]

        weakest = huge(1.0d0)
        do t = 1, size(temperatures)
            do p = 1, size(pressures)
                call measure_frozen_column(temperatures(t), pressures(p), PRESSURE_REFERENCE, STORAGE_REFERENCE, &
                                           margin, sensitivity, capacity, ok)
                if (.not. ok) then
                    call report_line("  FAIL the pressure column could not be measured")
                    failures = failures + 1
                    cycle
                end if
                weakest = min(weakest, margin)
            end do
        end do
        call check_exceeds("weakest unsaturated column       ", weakest, MARGIN_FLOOR, failures)

        ! Warm and saturated the Darcy driver alone carries the column.
        call measure_frozen_column(FROZEN_WARM, FROZEN_POSITIVE_PRESSURE, PRESSURE_REFERENCE, STORAGE_REFERENCE, &
                                   margin, sensitivity, capacity, ok)
        if (.not. ok) then
            call report_line("  FAIL the saturated warm column could not be measured")
            failures = failures + 1
            return
        end if
        call check_exceeds("saturated warm column survives   ", margin, MARGIN_FLOOR, failures)
        call check("it survives through p_l alone    ", sensitivity - 1.0d0, 1.0d-12, failures)

        ! Frozen and saturated, nothing is left.
        call measure_frozen_column(FROZEN_DEEP_COLD, FROZEN_POSITIVE_PRESSURE, PRESSURE_REFERENCE, STORAGE_REFERENCE, &
                                   margin, sensitivity, capacity, ok)
        if (.not. ok) then
            call report_line("  FAIL the saturated frozen column could not be measured")
            failures = failures + 1
            return
        end if
        call check("saturated frozen column is singular", margin, 1.0d-15, failures)
        call check("no liquid-pressure sensitivity   ", sensitivity, 1.0d-15, failures)
        call check("no storage capacity              ", capacity, 1.0d-20, failures)
    end subroutine check_frozen_rank

    !> One state's pressure column, measured and reported.
    subroutine measure_frozen_column(temperature, pressure, pressure_reference, storage_reference, &
                                     margin, sensitivity, capacity, ok)
        implicit none
        real(real64), intent(in) :: temperature, pressure
        real(real64), intent(in) :: pressure_reference, storage_reference
        real(real64), intent(inout) :: margin, sensitivity, capacity
        logical, intent(inout) :: ok

        character(len=200) :: line

        margin = 0.0d0
        sensitivity = 0.0d0
        capacity = 0.0d0
        ok = .false.
        call frozen_pressure_column(1_int32, temperature, pressure, pressure_reference, storage_reference, &
                                    margin, sensitivity, capacity, ok)
        if (.not. ok) return

        write (line, '(a,f7.3,a,es12.5,a,es12.5,a,es12.5,a,es12.5)') &
            "  column  T ", temperature, " C  p_w ", pressure, &
            "  dp_l/dp ", sensitivity, "  dW/dp ", capacity, "  margin ", margin
        call report_line(trim(line))
    end subroutine measure_frozen_column

    !> Record one logical measurement against what it has to be.
    subroutine check_logical(label, measured, expected, failures)
        implicit none
        character(len=*), intent(in) :: label
        logical, intent(in) :: measured
        logical, intent(in) :: expected
        integer(int32), intent(inout) :: failures

        character(len=160) :: line

        if (measured .eqv. expected) then
            write (line, '(a,a,a,l1,a,l1)') "  ok   ", label, "  measured ", measured, "  expected ", expected
        else
            write (line, '(a,a,a,l1,a,l1)') "  FAIL ", label, "  measured ", measured, "  expected ", expected
            failures = failures + 1
        end if
        call report_line(trim(line))
    end subroutine check_logical

    !> What the ice impedance is actually a function of.
    !>
    !> The factor \(\Omega=7\) was calibrated on the Mizoguchi columns with
    !> \(Q=\theta_i/(\theta_w+\theta_i-\theta_r)\), the ice fraction of the total
    !> water above residual. Passing the raw volumetric ice content instead would
    !> keep the same law and change what it means. Measured here by evaluating the
    !> repository's base and base-impedance models on the SAME published state:
    !> their ratio is the impedance factor the run applies, and it is compared
    !> with \(10^{-\Omega Q}\) formed from the contract's own reported contents.
    subroutine check_frozen_impedance(material, failures)
        implicit none
        type(type_mixed_frozen_material), intent(in) :: material
        integer(int32), intent(inout) :: failures

        type(type_mixed_frozen_response) :: response
        type(type_state) :: state
        type(type_iapws97) :: water_properties
        type(type_iapws06) :: ice_properties
        type(holder_hcfs) :: base_curve, impedance_curve
        type(type_config_hcf) :: base_config
        real(real64) :: base_value, impedance_value, measured
        real(real64) :: normalised_ratio, from_normalised, from_raw
        logical :: ok

        call water_properties%initialize()
        call ice_properties%initialize()

        base_config = material%conductivity_function
        base_config%model = HCF_MODES%BASE
        call base_curve%initialize(base_config, water_properties, ice_properties)
        call impedance_curve%initialize(material%conductivity_function, water_properties, ice_properties)

        ok = .false.
        call eval_frozen_response(1_int32, FROZEN_DEEP_COLD, FROZEN_WET_PRESSURE, response, ok)
        if (ok) call publish_frozen_state(1_int32, FROZEN_DEEP_COLD, FROZEN_WET_PRESSURE, state, ok)
        if (.not. ok) then
            call report_line("  FAIL the deep cold state could not be published")
            failures = failures + 1
            return
        end if

        base_value = 0.0d0
        impedance_value = 0.0d0
        call base_curve%p%calc_Kflh(state, base_value)
        call impedance_curve%p%calc_Kflh(state, impedance_value)
        measured = impedance_value/base_value

        normalised_ratio = response%ice_content &
                           /(response%liquid_content + response%ice_content - SOIL_THETA_R)
        from_normalised = 10.0d0**(-FROZEN_IMPEDANCE_FACTOR*normalised_ratio)
        from_raw = 10.0d0**(-FROZEN_IMPEDANCE_FACTOR*response%ice_content)

        call report_value("normalised ice ratio Q             ", normalised_ratio)
        call check("impedance is 10^(-omega Q), relative", &
                   (measured - from_normalised)/from_normalised, 1.0d-12, failures)
        ! Discriminating: the raw ice content is a different number, so the check
        ! above could have told the two apart.
        call check_exceeds("raw-ice impedance differs, relative", &
                           (from_raw - from_normalised)/from_normalised, 1.0d0, failures)
        call check("omega the run applied              ", &
                   material%conductivity_function%omega - FROZEN_IMPEDANCE_FACTOR, 0.0d0, failures)

        deallocate (base_curve%p)
        deallocate (impedance_curve%p)
    end subroutine check_frozen_impedance

    !> The water history the balance differences has to be exactly reproducible.
    !>
    !> A conservation audit differences \(W(T^n,p^n)\) against \(W(T^{n+1},p^{n+1})\),
    !> and a rejected step restores the first of those. If re-evaluating the same
    !> \((T,p_w)\) later gave a different last bit -- because a cache answered, or
    !> because an intermediate state changed something the closure carries -- the
    !> difference would not close. The requirement is therefore bitwise equality
    !> and not a tolerance.
    subroutine check_frozen_history(failures)
        implicit none
        integer(int32), intent(inout) :: failures

        type(type_mixed_frozen_response) :: first, again, other
        real(real64) :: difference
        logical :: ok

        ok = .false.
        call eval_frozen_response(1_int32, FROZEN_MILD_COLD, FROZEN_FIELD_PRESSURE, first, ok)
        if (.not. ok) then
            call report_line("  FAIL the history state could not be evaluated")
            failures = failures + 1
            return
        end if

        ! Enough other states to cycle every cached slot, at the states the
        ! tangents ask for as well as the four probes.
        call eval_frozen_response(1_int32, FROZEN_WARM, FROZEN_POSITIVE_PRESSURE, other, ok)
        call eval_frozen_response(1_int32, FROZEN_DEEP_COLD, FROZEN_WET_PRESSURE, other, ok)
        call eval_frozen_response(1_int32, FROZEN_WARM, FROZEN_FIELD_PRESSURE, other, ok)
        call eval_frozen_response(1_int32, FROZEN_MILD_COLD, FROZEN_POSITIVE_PRESSURE, other, ok)
        call eval_frozen_response(1_int32, FROZEN_DEEP_COLD, FROZEN_FIELD_PRESSURE, other, ok)
        call eval_frozen_response(1_int32, FROZEN_MILD_COLD, FROZEN_FIELD_PRESSURE, again, ok)
        if (.not. ok) then
            call report_line("  FAIL the history state could not be re-evaluated")
            failures = failures + 1
            return
        end if

        difference = again%stored_water - first%stored_water
        if (difference == 0.0d0) then
            call report_line("  ok   W re-evaluates bit for bit")
        else
            call report_line("  FAIL W does not re-evaluate bit for bit")
            failures = failures + 1
        end if
        call check("p_l re-evaluated, exactly          ", &
                   again%liquid_pressure - first%liquid_pressure, 0.0d0, failures)
        call check("theta_i re-evaluated, exactly      ", &
                   again%ice_content - first%ice_content, 0.0d0, failures)
        call check("k_r re-evaluated, exactly          ", &
                   again%relative_conductivity - first%relative_conductivity, 0.0d0, failures)
    end subroutine check_frozen_history

    !> Report a measured quantity that carries no pass or fail of its own.
    subroutine report_value(label, value)
        implicit none
        character(len=*), intent(in) :: label
        real(real64), intent(in) :: value

        character(len=160) :: line

        write (line, '(a,a,a,es12.5)') "  note ", label, "  measured ", value
        call report_line(trim(line))
    end subroutine report_value

    !> Take the largest of the first `count` entries over the communicator.
    subroutine reduce_max_array(values, count)
        implicit none
        real(real64), intent(inout) :: values(:)
        integer(int32), intent(in) :: count

        real(real64) :: reduced(size(values))
        integer(int32) :: mpi_err

        if (count < 1) return
        call MPI_Allreduce(values, reduced, count, MPI_DOUBLE_PRECISION, MPI_MAX, PETSC_COMM_WORLD, mpi_err)
        values(1:count) = reduced(1:count)
    end subroutine reduce_max_array

end program test_mixed_hydraulic
