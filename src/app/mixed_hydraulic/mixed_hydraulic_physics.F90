#include <petsc/finclude/petscsys.h>

!> Material response and run constants of the mixed hydraulic solver.
!>
!> The solved unknowns are the Darcy flux \(\mathbf q_l\) in m/s and a pressure
!> in Pa, related by
!>
!> \[
!>   \mathbf q_l = -\mathbf\Lambda(\nabla p_l-\rho_l\mathbf g),\qquad
!>   \mathbf\Lambda=\frac{\mathbf K_h}{\rho_{\mathrm{ref}}g},
!> \]
!>
!> with \(\mathbf K_h\) the hydraulic conductivity in m/s. The weak form uses
!> \(\mathbf\Lambda^{-1}\), so the mobility has to be positive definite: a zero
!> conductivity is refused here rather than floored, because a floor would turn
!> an impermeable subdomain into a slightly permeable one without saying so.
!>
!> An unsaturated material splits that conductivity into a saturated tensor and
!> a scalar relative conductivity of the pressure,
!>
!> \[
!>   \mathbf K_h(p)=k_r(p)\,\mathbf K_{h,\mathrm{sat}},\qquad
!>   \mathbf\Lambda^{-1}(p)=\frac{\mathbf\Lambda_{\mathrm{sat}}^{-1}}{k_r(p)},
!> \]
!>
!> so only the tensor is inverted, once per cell, and the state dependence stays
!> a positive scalar. A positive \(k_r\) times a symmetric positive definite
!> tensor is symmetric positive definite, which is what keeps the mobility
!> admissible at every iterate. The curve itself is the repository's own
!> hydraulic conductivity function; see `app_mixed_hydraulic_closure`.
!>
!> The unknowns PETSc stores are dimensionless. A flux scale in m/s and a
!> pressure scale in Pa carry the physical size of the problem, which is what
!> keeps a conductivity of \(10^{-20}\) m/s and a pressure of \(10^{4}\) Pa in
!> the same matrix. Everything leaving this solver -- balances, reconstructed
!> fluxes, reported pressures -- is multiplied back into physical units.
module app_mixed_hydraulic_physics
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: petscsys
    use :: module_core, only:type_config_wrf, type_config_hcf
    implicit none
    private

    public :: type_mixed_material
    public :: invert_mobility
    public :: intrinsic_from_physical
    public :: physical_from_intrinsic
    public :: MIXED_CONSTANTS
    public :: C_DT, C_FLUX_SCALE, C_PRESSURE_SCALE, C_LENGTH_SCALE
    public :: C_DENSITY, C_GRAVITY_X, C_GRAVITY_Z
    public :: MIXED_AUX_FIELDS
    public :: A_INVERSE_MOBILITY, A_STORAGE, A_SOURCE, A_PRESSURE_OLD, A_BC_PRESSURE
    public :: A_MATERIAL
    public :: MIXED_FIELD_FLUX, MIXED_FIELD_PRESSURE
    public :: MIXED_DIMENSION
    public :: MIXED_STORAGE_LINEAR, MIXED_STORAGE_RETENTION
    public :: MIXED_CONDUCTIVITY_FIXED, MIXED_CONDUCTIVITY_HCF

    !> Field numbers of the isothermal proof of concept, PETSc's zero-based
    !> numbering. A coupled run adds temperature ahead of these, so nothing
    !> outside this module may assume the values.
    PetscInt, parameter :: MIXED_FIELD_FLUX = 0
    PetscInt, parameter :: MIXED_FIELD_PRESSURE = 1

    !> The intrinsic geometry this solver works in. A planar mesh given as
    !> \((x,0,z)\) is carried as \((x,z)\); PETSc's H(div) transforms are
    !> defined only where the embedding dimension equals the reference one.
    integer(int32), parameter :: MIXED_DIMENSION = 2

    ! Run constants, in the order PetscDSSetConstants receives them.
    integer(int32), parameter :: C_DT = 1             !< Time step, s.
    integer(int32), parameter :: C_FLUX_SCALE = 2     !< Flux unit of the unknown, m/s.
    integer(int32), parameter :: C_PRESSURE_SCALE = 3 !< Pressure unit of the unknown, Pa.
    integer(int32), parameter :: C_LENGTH_SCALE = 4   !< Length the conservation row is scaled by, m.
    integer(int32), parameter :: C_DENSITY = 5        !< Reference liquid density, kg/m3.
    integer(int32), parameter :: C_GRAVITY_X = 6      !< Gravity along intrinsic x, m/s2.
    integer(int32), parameter :: C_GRAVITY_Z = 7      !< Gravity along intrinsic z, m/s2.
    integer(int32), parameter :: MIXED_CONSTANTS = 7

    ! Auxiliary fields, cell-wise. The order fixes the offsets PETSc hands the
    ! kernels, so it matches the order the fields are added to the auxiliary DM.
    !> \(\mathbf\Lambda_{\mathrm{sat}}^{-1}\), row-major, Pa s/m2.
    !>
    !> The saturated inverse mobility, inverted once per cell from that cell's
    !> own tensor. A kernel divides it by the relative conductivity of the
    !> current pressure, so nothing inverts a state-dependent tensor and nothing
    !> inverts an average of two materials.
    integer(int32), parameter :: A_INVERSE_MOBILITY = 1
    !> Whether the cell's storage law can store water at all, in 1/Pa.
    !>
    !> A magnitude of the law rather than a value of \(\partial W/\partial p\)
    !> at one state: it is written once, at setup, and it is zero exactly when
    !> the cell is steady. The pressure gauge reads it to decide whether the
    !> constant pressure is really in the operator's kernel, and a curve whose
    !> tangent happens to vanish at the accepted state still stores water at the
    !> next one. No kernel reads it; the residual and its tangent both go
    !> through the storage closure.
    integer(int32), parameter :: A_STORAGE = 2
    integer(int32), parameter :: A_SOURCE = 3           !< Volumetric source, 1/s.
    integer(int32), parameter :: A_PRESSURE_OLD = 4     !< Accepted cell pressure, Pa.
    integer(int32), parameter :: A_BC_PRESSURE = 5      !< Prescribed boundary pressure, Pa.
    !> The cell's slot in the storage-law table, one-based, as a whole number.
    !>
    !> A pointwise function receives no context pointer, so this is how a kernel
    !> reaches the law of the cell it is running on.
    integer(int32), parameter :: A_MATERIAL = 6
    integer(int32), parameter :: MIXED_AUX_FIELDS = 6

    !> A constant capacity: \(W=Cp\), so \(\partial W/\partial p=C\).
    integer(int32), parameter :: MIXED_STORAGE_LINEAR = 1
    !> The repository's retention curve plus its saturated compressive storage.
    integer(int32), parameter :: MIXED_STORAGE_RETENTION = 2

    !> \(k_r\equiv1\): the conductivity tensor is the whole law.
    integer(int32), parameter :: MIXED_CONDUCTIVITY_FIXED = 1
    !> The repository's hydraulic conductivity function supplies \(k_r(p)\).
    integer(int32), parameter :: MIXED_CONDUCTIVITY_HCF = 2

    !> One material's hydraulic response, constant over its cells.
    type :: type_mixed_material
        !> The "Cell Sets" value this material is written on.
        integer(int32) :: id = 0
        !> Saturated conductivity \(\mathbf K_{h,\mathrm{sat}}\), m/s, symmetric positive definite.
        real(real64) :: conductivity(MIXED_DIMENSION, MIXED_DIMENSION) = 0.0_real64
        !> Linear storage \(\partial W/\partial p\) in 1/Pa; zero makes the cell steady.
        real(real64) :: storage = 0.0_real64
        !> Volumetric source in 1/s.
        real(real64) :: source = 0.0_real64
        !> MIXED_STORAGE_LINEAR or MIXED_STORAGE_RETENTION.
        integer(int32) :: storage_law = MIXED_STORAGE_LINEAR
        !> Retention parameters, read for MIXED_STORAGE_RETENTION.
        type(type_config_wrf) :: retention
        !> Specific storage \(S_s\) of the saturated branch, 1/m.
        real(real64) :: specific_storage = 0.0_real64
        !> Pressure at which the pore volume fills, Pa.
        !>
        !> The value the closure's own `calc_saturation_pressure` returns for
        !> this material's state. In the isothermal unfrozen case with
        !> \(\theta_s\) equal to the porosity that is zero, because every
        !> retention curve here is saturated exactly at a non-negative head.
        real(real64) :: saturation_pressure = 0.0_real64
        !> MIXED_CONDUCTIVITY_FIXED or MIXED_CONDUCTIVITY_HCF.
        integer(int32) :: conductivity_law = MIXED_CONDUCTIVITY_FIXED
        !> The conductivity function, read for MIXED_CONDUCTIVITY_HCF.
        !>
        !> Its `k_sat` sets the scale of the model's own output only; the
        !> saturated conductivity of the cell is the tensor above, and the
        !> relative conductivity is this model normalised by its own saturated
        !> value. A material may therefore be anisotropic while its curve is not.
        type(type_config_hcf) :: conductivity_function
    end type type_mixed_material

contains

    !> The inverse pressure mobility \(\mathbf\Lambda^{-1}=\rho_{\mathrm{ref}}g\,\mathbf K_h^{-1}\).
    !>
    !> \[
    !>   \mathbf\Lambda^{-1}=\rho_{\mathrm{ref}}g\,\mathbf K_h^{-1},\qquad
    !>   \mathbf K_h\ \text{symmetric positive definite}.
    !> \]
    !>
    !> Assumes a two-dimensional intrinsic tensor. Guarantees an exact inverse of
    !> the given tensor up to the rounding of one two-by-two solve; no iteration
    !> and no conditioning estimate. Arithmetic and memory are \(O(1)\).
    !>
    !> Fails with `ok = .false.` and an untouched result when the tensor is not
    !> symmetric positive definite, which includes a zero conductivity. No floor
    !> is applied: an impermeable region is a different model, not a small number.
    subroutine invert_mobility(conductivity, density, gravity_magnitude, inverse_mobility, ok)
        implicit none
        !> Hydraulic conductivity in m/s.
        real(real64), intent(in) :: conductivity(MIXED_DIMENSION, MIXED_DIMENSION)
        !> Reference liquid density in kg/m3, positive.
        real(real64), intent(in) :: density
        !> Magnitude of gravity in m/s2, positive.
        real(real64), intent(in) :: gravity_magnitude
        !> \(\mathbf\Lambda^{-1}\) in Pa s/m2, row-major, written only when `ok`.
        real(real64), intent(inout) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION)
        !> True when the conductivity was symmetric positive definite.
        logical, intent(inout) :: ok

        real(real64) :: determinant, weight, asymmetry, scale

        ok = .false.

        if (density <= 0.0_real64 .or. gravity_magnitude <= 0.0_real64) return

        scale = maxval(abs(conductivity))
        if (scale <= 0.0_real64) return

        asymmetry = abs(conductivity(1, 2) - conductivity(2, 1))
        if (asymmetry > 1.0e-12_real64*scale) return

        determinant = conductivity(1, 1)*conductivity(2, 2) - conductivity(1, 2)*conductivity(2, 1)
        if (conductivity(1, 1) <= 0.0_real64 .or. conductivity(2, 2) <= 0.0_real64) return
        if (determinant <= 0.0_real64) return

        ! rho g K^{-1}, with the two-by-two inverse written out.
        weight = density*gravity_magnitude/determinant
        inverse_mobility(1) = weight*conductivity(2, 2)
        inverse_mobility(2) = -weight*conductivity(1, 2)
        inverse_mobility(3) = -weight*conductivity(2, 1)
        inverse_mobility(4) = weight*conductivity(1, 1)
        ok = .true.
    end subroutine invert_mobility

    !> The intrinsic \((x,z)\) of a planar point given as \((x,0,z)\).
    pure function intrinsic_from_physical(physical) result(intrinsic_point)
        implicit none
        !> A point of the planar mesh in its three stored coordinates.
        real(real64), intent(in) :: physical(3)
        real(real64) :: intrinsic_point(MIXED_DIMENSION)

        intrinsic_point(1) = physical(1)
        intrinsic_point(2) = physical(3)
    end function intrinsic_from_physical

    !> The \((x,0,z)\) of an intrinsic \((x,z)\) point or vector.
    !>
    !> Used for output and for reading a project's gravity: the out-of-plane
    !> component of a planar field is zero by construction, not by rounding.
    pure function physical_from_intrinsic(intrinsic_point) result(physical)
        implicit none
        !> A point or vector in the intrinsic plane.
        real(real64), intent(in) :: intrinsic_point(MIXED_DIMENSION)
        real(real64) :: physical(3)

        physical(1) = intrinsic_point(1)
        physical(2) = 0.0_real64
        physical(3) = intrinsic_point(2)
    end function physical_from_intrinsic

end module app_mixed_hydraulic_physics
