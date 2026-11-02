#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscsnes.h>

!> The primal water flux with a nodal liquid pressure, as `bin/FTCMS` builds it.
!>
!> The same primal row as `app_mixed_hydraulic_column_primal`, with the two
!> choices the CG code of this repository makes for the Darcy term
!> (`src/physics/governing/hydraulic/hydraulic_matrix.F90`):
!>
!> \[
!>   \mathbf q_l^K=-\bar k_r^K\,\mathbf M\Bigl(\sum_i p_l(T_i,p_{w,i})\nabla\phi_i-\rho_l\mathbf g\Bigr),
!>   \qquad
!>   \bar k_r^K=\Bigl(\prod_{i\in K}k_r(T_i,p_{w,i})\Bigr)^{1/n_K}.
!> \]
!>
!> The liquid pressure is evaluated at the nodes and interpolated, so its
!> gradient is that of a continuous function of the nodal unknowns and not the
!> chain rule \(p_{l,p}\nabla p_w+p_{l,T}\nabla T\) at a quadrature point. The
!> element conductivity is the geometric mean of the nodal values, which cannot
!> exceed the largest of them across a freezing front.
!>
!> ## How it reaches the solver
!>
!> PETSc's pointwise kernels see the fields at quadrature points only, never
!> the nodal values, so this term is assembled here, element by element, and
!> added to the residual and Jacobian PETSc's own element loop computed. The
!> storage term and the heat row stay in that loop. The Jacobian is exact in
!> the nodal unknowns, so the solver is the same Newton iteration as for the
!> mixed form.
!>
!> ## Ownership
!>
!> One column at a time, like the frozen closure it evaluates: the hooks and
!> the reference tabulation are module state.
module app_mixed_hydraulic_column_nodal
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_funloc, c_null_ptr, c_associated
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use :: petscdmplex
    use :: petscsnes
    use :: core_interop_petscds, only:petscds_set_snes_local_residual, petscds_set_snes_local_jacobian, &
                                      petscds_snes_local_residual_fem, petscds_snes_local_jacobian_fem
    use :: core_interop_petscfe, only:petscfe_tabulate, petscfe_quadrature, petscfe_field_size
    use :: app_mixed_hydraulic_physics, only:MIXED_DIMENSION
    use :: app_mixed_hydraulic_frozen, only:type_mixed_frozen_response, type_mixed_frozen_tangents, &
                                            eval_frozen_response, eval_frozen_tangents
    use :: app_mixed_hydraulic_balance, only:MAX_QUADRATURE_POINTS
    use :: app_mixed_hydraulic_column_pointwise, only:CA_MATERIAL, CA_INVERSE_MOBILITY, note_column_fault
    implicit none
    private

    public :: install_nodal_flux
    public :: clear_nodal_flux
    public :: compute_nodal_cell_flux

    !> The largest number of nodes of one cell.
    integer(int32), parameter :: MAX_NODES = 4

    !> What the hooks need beyond the arguments PETSc passes them.
    type :: type_nodal_flux_context
        logical :: active = .false.
        !> The solve DM the hooks were installed on.
        type(c_ptr) :: dm = c_null_ptr
        !> The auxiliary DM and a vector of it carrying each cell's material
        !> slot and saturated inverse mobility.
        type(c_ptr) :: aux_dm = c_null_ptr
        type(c_ptr) :: auxiliary = c_null_ptr
        !> Units of the unknowns, Pa and K; the water row's divisor
        !> \(q_{\mathrm{scale}}/L\) is applied as \(L/q_{\mathrm{scale}}\).
        real(real64) :: pressure_scale = 1.0_real64
        real(real64) :: temperature_scale = 1.0_real64
        real(real64) :: row_scale = 1.0_real64
        !> Reference liquid density, kg/m3, and gravity, m/s2.
        real(real64) :: density = 1000.0_real64
        real(real64) :: gravity(MIXED_DIMENSION) = 0.0_real64
        !> Basis functions per field and quadrature points per cell.
        integer(int32) :: num_basis = 0
        integer(int32) :: num_points = 0
        !> Reference quadrature weights, basis values and reference gradients.
        real(real64) :: weights(MAX_QUADRATURE_POINTS) = 0.0_real64
        real(real64) :: values(MAX_NODES, MAX_QUADRATURE_POINTS) = 0.0_real64
        real(real64) :: gradients(MIXED_DIMENSION, MAX_NODES, MAX_QUADRATURE_POINTS) = 0.0_real64
    end type type_nodal_flux_context

    type(type_nodal_flux_context), save :: nodal

    !> One cell's element quantities.
    type :: type_nodal_cell
        integer(int32) :: slot = 0
        !> \(|\det J|\) of the affine map, and the physical basis gradients at
        !> every quadrature point.
        real(real64) :: jacobian_determinant = 0.0_real64
        real(real64) :: gradients(MIXED_DIMENSION, MAX_NODES, MAX_QUADRATURE_POINTS) = 0.0_real64
        !> Saturated mobility, m2/(Pa s).
        real(real64) :: mobility(MIXED_DIMENSION, MIXED_DIMENSION) = 0.0_real64
        !> Nodal liquid pressure and its derivatives, Pa, 1 and Pa/K.
        real(real64) :: liquid(MAX_NODES) = 0.0_real64
        real(real64) :: liquid_dpressure(MAX_NODES) = 0.0_real64
        real(real64) :: liquid_dtemperature(MAX_NODES) = 0.0_real64
        !> Nodal \(\partial\ln k_r\) in each unknown, 1/Pa and 1/K.
        real(real64) :: log_conductivity_dpressure(MAX_NODES) = 0.0_real64
        real(real64) :: log_conductivity_dtemperature(MAX_NODES) = 0.0_real64
        !> The geometric-mean relative conductivity of the cell.
        real(real64) :: conductivity = 0.0_real64
    end type type_nodal_cell

contains

    !> Install the nodal flux on a primal column whose DS carries the storage only.
    !>
    !> Tabulates the temperature element at its own quadrature once, and
    !> replaces the DM's local SNES residual and Jacobian with PETSc's own plus
    !> the element flux of this module.
    !>
    !> Assumes a primal layout [T, p] with the two fields in the same
    !> continuous degree-one space, affine cells, and an auxiliary vector that
    !> lives as long as the solver. \(O(1)\). Fails by returning a nonzero
    !> `ierr`; by `error stop` on cells with more nodes than this module holds.
    subroutine install_nodal_flux(dm, aux_dm, auxiliary, pressure_scale, temperature_scale, flux_scale, &
                                  length_scale, density, gravity, ierr)
        implicit none
        !> The solve DM, with its DS and SNES already set up.
        DM, intent(in) :: dm
        !> The auxiliary DM and a vector of it.
        DM, intent(in) :: aux_dm
        Vec, intent(in) :: auxiliary
        !> Units of the pressure and temperature unknowns, Pa and K.
        real(real64), intent(in) :: pressure_scale, temperature_scale
        !> The flux unit, m/s, and the length the water row is scaled by, m.
        real(real64), intent(in) :: flux_scale, length_scale
        !> Reference liquid density, kg/m3, and gravity, m/s2.
        real(real64), intent(in) :: density
        real(real64), intent(in) :: gravity(MIXED_DIMENSION)
        PetscErrorCode, intent(inout) :: ierr

        PetscObject :: discretisation
        real(c_double) :: points(MIXED_DIMENSION, MAX_QUADRATURE_POINTS), weights(MAX_QUADRATURE_POINTS)
        real(c_double), allocatable :: values(:, :), gradients(:, :, :)
        type(c_ptr) :: dm_handle, fe_handle
        integer(c_int) :: cerr, count, order, num_basis, num_components

        dm_handle = transfer(dm%v, dm_handle)
        cerr = 0
        call petscfe_field_size(dm_handle, 0_c_int, num_basis, num_components, cerr)
        if (cerr == 0) call petscfe_quadrature(dm_handle, 0_c_int, points, weights, count, order, cerr)
        ierr = int(cerr, kind(ierr))
        if (ierr /= 0) return
        if (num_basis > MAX_NODES) error stop "The nodal column flux holds cells of at most four nodes"

        call DMGetField(dm, 0_PETSC_INT_KIND, PETSC_NULL_DMLABEL, discretisation, ierr)
        if (ierr /= 0) return
        fe_handle = transfer(discretisation%v, fe_handle)
        allocate (values(num_basis, count), gradients(MIXED_DIMENSION, num_basis, count))
        call petscfe_tabulate(fe_handle, points(:, :count), values, gradients, cerr)
        ierr = int(cerr, kind(ierr))
        if (ierr /= 0) return

        nodal%num_basis = int(num_basis, int32)
        nodal%num_points = int(count, int32)
        nodal%weights(:count) = real(weights(:count), real64)
        nodal%values(:num_basis, :count) = real(values, real64)
        nodal%gradients(:, :num_basis, :count) = real(gradients, real64)
        nodal%pressure_scale = pressure_scale
        nodal%temperature_scale = temperature_scale
        nodal%row_scale = length_scale/flux_scale
        nodal%density = density
        nodal%gravity = gravity
        nodal%dm = dm_handle
        nodal%aux_dm = transfer(aux_dm%v, nodal%aux_dm)
        nodal%auxiliary = transfer(auxiliary%v, nodal%auxiliary)

        call petscds_set_snes_local_residual(dm_handle, c_funloc(nodal_residual_local), c_null_ptr, cerr)
        if (cerr == 0) call petscds_set_snes_local_jacobian(dm_handle, c_funloc(nodal_jacobian_local), &
                                                            c_null_ptr, cerr)
        ierr = int(cerr, kind(ierr))
        nodal%active = (ierr == 0)
    end subroutine install_nodal_flux

    !> Forget the installed column. \(O(1)\); cannot fail.
    subroutine clear_nodal_flux()
        implicit none

        nodal = type_nodal_flux_context()
    end subroutine clear_nodal_flux

    !> The local residual with the nodal flux added to it.
    function nodal_residual_local(dm_handle, state_handle, residual_handle, ctx) &
        bind(C, name="ftcms_column_nodal_residual") result(ierr)
        implicit none
        type(c_ptr), value :: dm_handle
        type(c_ptr), value :: state_handle
        type(c_ptr), value :: residual_handle
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        DM :: dm
        Vec :: local_state, local_residual
        PetscSection :: section
        PetscInt :: cell, c_start, c_end
        PetscScalar :: element(2*MAX_NODES)
        real(real64) :: rows(MAX_NODES)
        integer(c_int) :: cerr
        integer(int32) :: n
        logical :: ok

        cerr = 0
        call petscds_snes_local_residual_fem(dm_handle, state_handle, residual_handle, ctx, cerr)
        ierr = int(cerr, kind(ierr))
        if (cerr /= 0) return
        call check_installed_dm(dm_handle)

        dm%v = transfer(dm_handle, dm%v)
        local_state%v = transfer(state_handle, local_state%v)
        local_residual%v = transfer(residual_handle, local_residual%v)
        n = nodal%num_basis

        call DMGetLocalSection(dm, section, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return

        do cell = c_start, c_end - 1
            call element_flux(dm, section, local_state, cell, rows, ok, ierr)
            if (ierr /= 0) return
            ! The closure is field-major: the temperature's nodes, then the pressure's.
            element = 0.0_PETSC_REAL_KIND
            if (ok) then
                element(n + 1:2*n) = real(rows(:n), PETSC_REAL_KIND)
            else
                ! A refused state refuses the step, as a pointwise kernel's does.
                element(n + 1:2*n) = real(ieee_value(1.0_real64, ieee_quiet_nan), PETSC_REAL_KIND)
            end if
            call DMPlexVecSetClosure(dm, section, local_residual, cell, element(:2*n), ADD_VALUES, ierr)
            if (ierr /= 0) return
        end do
    end function nodal_residual_local

    !> The local Jacobian with the nodal flux added to it.
    function nodal_jacobian_local(dm_handle, state_handle, jacobian_handle, preconditioner_handle, ctx) &
        bind(C, name="ftcms_column_nodal_jacobian") result(ierr)
        implicit none
        type(c_ptr), value :: dm_handle
        type(c_ptr), value :: state_handle
        type(c_ptr), value :: jacobian_handle
        type(c_ptr), value :: preconditioner_handle
        type(c_ptr), value :: ctx
        PetscErrorCode :: ierr

        DM :: dm
        Vec :: local_state
        Mat :: jacobian, preconditioner
        integer(c_int) :: cerr

        cerr = 0
        call petscds_snes_local_jacobian_fem(dm_handle, state_handle, jacobian_handle, &
                                             preconditioner_handle, ctx, cerr)
        ierr = int(cerr, kind(ierr))
        if (cerr /= 0) return
        call check_installed_dm(dm_handle)

        dm%v = transfer(dm_handle, dm%v)
        local_state%v = transfer(state_handle, local_state%v)
        preconditioner%v = transfer(preconditioner_handle, preconditioner%v)
        call add_flux_jacobian(dm, local_state, preconditioner, ierr)
        if (ierr /= 0) return

        ! PETSc hands the same matrix twice unless a separate preconditioning
        ! matrix was asked for; the flux belongs in both.
        if (.not. c_associated(jacobian_handle, preconditioner_handle)) then
            jacobian%v = transfer(jacobian_handle, jacobian%v)
            call add_flux_jacobian(dm, local_state, jacobian, ierr)
        end if
    end function nodal_jacobian_local

    !> Add every cell's flux block to an assembled matrix and assemble it again.
    subroutine add_flux_jacobian(dm, local_state, matrix, ierr)
        implicit none
        DM, intent(in) :: dm
        Vec, intent(in) :: local_state
        Mat, intent(inout) :: matrix
        PetscErrorCode, intent(inout) :: ierr

        PetscSection :: section, global_section
        PetscInt :: cell, c_start, c_end
        PetscScalar :: element(4*MAX_NODES*MAX_NODES)
        real(real64) :: block(MAX_NODES, 2*MAX_NODES)
        integer(int32) :: n, r, c
        logical :: ok

        call DMGetLocalSection(dm, section, ierr)
        if (ierr == 0) call DMGetGlobalSection(dm, global_section, ierr)
        if (ierr == 0) call DMPlexGetHeightStratum(dm, 0_PETSC_INT_KIND, c_start, c_end, ierr)
        if (ierr /= 0) return
        n = nodal%num_basis

        do cell = c_start, c_end - 1
            call element_flux_jacobian(dm, section, local_state, cell, block, ok, ierr)
            if (ierr /= 0) return
            if (.not. ok) cycle
            ! Row-major over the closure [T nodes, p nodes]; only the p rows carry the flux.
            element(:4*n*n) = 0.0_PETSC_REAL_KIND
            do r = 1, n
                do c = 1, 2*n
                    element((n + r - 1)*2*n + c) = real(block(r, c), PETSC_REAL_KIND)
                end do
            end do
            call DMPlexMatSetClosure(dm, section, global_section, matrix, cell, element(:4*n*n), ADD_VALUES, ierr)
            if (ierr /= 0) return
        end do
        call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
        if (ierr == 0) call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
    end subroutine add_flux_jacobian

    !> The cell's geometry, mobility and nodal closure, from the local state.
    subroutine read_nodal_cell(dm, section, local_state, cell, element, ok, ierr)
        implicit none
        DM, intent(in) :: dm
        PetscSection, intent(in) :: section
        Vec, intent(in) :: local_state
        PetscInt, intent(in) :: cell
        type(type_nodal_cell), intent(inout) :: element
        !> False when the closure could not answer at one of the nodes.
        logical, intent(inout) :: ok
        PetscErrorCode, intent(inout) :: ierr

        DM :: aux_dm
        Vec :: auxiliary
        PetscSection :: aux_section
        PetscScalar, pointer :: closure(:), aux_raw(:)
        PetscInt :: closure_size, offset
        PetscReal :: v0(MIXED_DIMENSION), jac(MIXED_DIMENSION*MIXED_DIMENSION)
        PetscReal :: inverse_jac(MIXED_DIMENSION*MIXED_DIMENSION), determinant
        type(type_mixed_frozen_response) :: response
        type(type_mixed_frozen_tangents) :: tangents
        real(real64) :: inverse_mobility(MIXED_DIMENSION*MIXED_DIMENSION), log_mean
        real(real64) :: temperature, pressure, det_mobility
        integer(int32) :: i, q, d, e, n, k

        ok = .false.
        n = nodal%num_basis

        call DMPlexComputeCellGeometryAffineFEM(dm, cell, v0, jac, inverse_jac, determinant, ierr)
        if (ierr /= 0) return
        element%jacobian_determinant = abs(real(determinant, real64))
        ! Physical gradient: d(phi)/dx_d = sum_e d(phi)/d(xi_e) invJ(e, d), invJ row major.
        do q = 1, nodal%num_points
            do i = 1, n
                do d = 1, MIXED_DIMENSION
                    element%gradients(d, i, q) = 0.0_real64
                    do e = 1, MIXED_DIMENSION
                        element%gradients(d, i, q) = element%gradients(d, i, q) &
                            + nodal%gradients(e, i, q)*real(inverse_jac((e - 1)*MIXED_DIMENSION + d), real64)
                    end do
                end do
            end do
        end do

        ! Material slot and saturated inverse mobility of the cell.
        aux_dm%v = transfer(nodal%aux_dm, aux_dm%v)
        auxiliary%v = transfer(nodal%auxiliary, auxiliary%v)
        call DMGetLocalSection(aux_dm, aux_section, ierr)
        if (ierr /= 0) return
        call VecGetArrayRead(auxiliary, aux_raw, ierr)
        if (ierr /= 0) return
        call PetscSectionGetFieldOffset(aux_section, cell, int(CA_MATERIAL - 1, PETSC_INT_KIND), offset, ierr)
        if (ierr == 0) element%slot = nint(real(aux_raw(offset + 1), real64), int32)
        if (ierr == 0) call PetscSectionGetFieldOffset(aux_section, cell, int(CA_INVERSE_MOBILITY - 1, PETSC_INT_KIND), &
                                                       offset, ierr)
        if (ierr == 0) then
            do k = 1, MIXED_DIMENSION*MIXED_DIMENSION
                inverse_mobility(k) = real(aux_raw(offset + k), real64)
            end do
        end if
        call VecRestoreArrayRead(auxiliary, aux_raw, ierr)
        if (ierr /= 0) return
        det_mobility = inverse_mobility(1)*inverse_mobility(4) - inverse_mobility(2)*inverse_mobility(3)
        element%mobility(1, 1) = inverse_mobility(4)/det_mobility
        element%mobility(1, 2) = -inverse_mobility(2)/det_mobility
        element%mobility(2, 1) = -inverse_mobility(3)/det_mobility
        element%mobility(2, 2) = inverse_mobility(1)/det_mobility

        ! The nodal unknowns, field-major, and the closure at every node.
        nullify (closure)
        closure_size = 0_PETSC_INT_KIND
        call DMPlexVecGetClosure(dm, section, local_state, cell, closure_size, closure, ierr)
        if (ierr /= 0) return
        log_mean = 0.0_real64
        ok = (int(closure_size, int32) == 2*n)
        do i = 1, n
            if (.not. ok) exit
            temperature = nodal%temperature_scale*real(closure(i), real64)
            pressure = nodal%pressure_scale*real(closure(n + i), real64)
            ok = .false.
            call eval_frozen_response(element%slot, temperature, pressure, response, ok)
            if (ok) call eval_frozen_tangents(element%slot, temperature, pressure, tangents, ok)
            if (.not. ok) then
                call note_column_fault(element%slot, temperature, pressure)
                exit
            end if
            element%liquid(i) = response%liquid_pressure
            element%liquid_dpressure(i) = tangents%dliquid_dpressure
            element%liquid_dtemperature(i) = tangents%dliquid_dtemperature
            element%log_conductivity_dpressure(i) = tangents%dlog_conductivity_dpressure
            element%log_conductivity_dtemperature(i) = tangents%dlog_conductivity_dtemperature
            log_mean = log_mean + log(response%relative_conductivity)
        end do
        call DMPlexVecRestoreClosure(dm, section, local_state, cell, closure_size, closure, ierr)
        if (.not. ok) return
        element%conductivity = exp(log_mean/real(n, real64))
    end subroutine read_nodal_cell

    !> \(\nabla\hat p_l-\rho_l\mathbf g\) at one quadrature point, Pa/m.
    pure function nodal_driver(element, q) result(driver)
        implicit none
        type(type_nodal_cell), intent(in) :: element
        integer(int32), intent(in) :: q
        real(real64) :: driver(MIXED_DIMENSION)

        integer(int32) :: i

        driver = -nodal%density*nodal%gravity
        do i = 1, nodal%num_basis
            driver = driver + element%liquid(i)*element%gradients(:, i, q)
        end do
    end function nodal_driver

    !> The flux rows of one cell: \(\int_K\bar k_r\mathbf M(\nabla\hat p_l-\rho_l\mathbf g)\cdot\nabla\phi_j\,L/q_{\mathrm{scale}}\).
    subroutine element_flux(dm, section, local_state, cell, rows, ok, ierr)
        implicit none
        DM, intent(in) :: dm
        PetscSection, intent(in) :: section
        Vec, intent(in) :: local_state
        PetscInt, intent(in) :: cell
        real(real64), intent(inout) :: rows(MAX_NODES)
        logical, intent(inout) :: ok
        PetscErrorCode, intent(inout) :: ierr

        type(type_nodal_cell) :: element
        real(real64) :: flow(MIXED_DIMENSION), weight
        integer(int32) :: j, q

        rows = 0.0_real64
        call read_nodal_cell(dm, section, local_state, cell, element, ok, ierr)
        if (ierr /= 0 .or. .not. ok) return

        do q = 1, nodal%num_points
            weight = nodal%weights(q)*element%jacobian_determinant*nodal%row_scale
            flow = element%conductivity*matmul(element%mobility, nodal_driver(element, q))
            do j = 1, nodal%num_basis
                rows(j) = rows(j) + weight*dot_product(flow, element%gradients(:, j, q))
            end do
        end do
    end subroutine element_flux

    !> The exact derivative of `element_flux` in the cell's nodal unknowns.
    !>
    !> \[
    !>   \frac{\partial R_j}{\partial u_i}=s_u\int_K\Bigl(\bar k_r\,\mathbf M\,\frac{\partial p_{l,i}}{\partial u_i}
    !>   \nabla\phi_i+\frac{\bar k_r}{n_K}\frac{\partial\ln k_{r,i}}{\partial u_i}\,\mathbf M\mathbf d\Bigr)
    !>   \cdot\nabla\phi_j\,\frac{L}{q_{\mathrm{scale}}},
    !> \]
    !>
    !> with \(u\) the temperature or the pressure of node \(i\) and \(s_u\) its
    !> unit. Columns 1..n are the temperature nodes, n+1..2n the pressure nodes.
    subroutine element_flux_jacobian(dm, section, local_state, cell, block, ok, ierr)
        implicit none
        DM, intent(in) :: dm
        PetscSection, intent(in) :: section
        Vec, intent(in) :: local_state
        PetscInt, intent(in) :: cell
        real(real64), intent(inout) :: block(MAX_NODES, 2*MAX_NODES)
        logical, intent(inout) :: ok
        PetscErrorCode, intent(inout) :: ierr

        type(type_nodal_cell) :: element
        real(real64) :: driven(MIXED_DIMENSION), along(MIXED_DIMENSION), weight, mean_share
        integer(int32) :: i, j, q, n

        block = 0.0_real64
        call read_nodal_cell(dm, section, local_state, cell, element, ok, ierr)
        if (ierr /= 0 .or. .not. ok) return
        n = nodal%num_basis
        mean_share = element%conductivity/real(n, real64)

        do q = 1, nodal%num_points
            weight = nodal%weights(q)*element%jacobian_determinant*nodal%row_scale
            driven = matmul(element%mobility, nodal_driver(element, q))
            do i = 1, n
                along = element%conductivity*matmul(element%mobility, element%gradients(:, i, q))
                do j = 1, n
                    block(j, i) = block(j, i) + weight*nodal%temperature_scale &
                        *dot_product(element%liquid_dtemperature(i)*along &
                                     + mean_share*element%log_conductivity_dtemperature(i)*driven, &
                                     element%gradients(:, j, q))
                    block(j, n + i) = block(j, n + i) + weight*nodal%pressure_scale &
                        *dot_product(element%liquid_dpressure(i)*along &
                                     + mean_share*element%log_conductivity_dpressure(i)*driven, &
                                     element%gradients(:, j, q))
                end do
            end do
        end do
    end subroutine element_flux_jacobian

    !> The cell-average liquid flux of the nodal form, m/s.
    !>
    !> \(\bar{\mathbf q}_l^K=-\bar k_r\mathbf M\,|K|^{-1}\int_K(\nabla\hat p_l-\rho_l\mathbf g)\),
    !> on the quadrature the residual uses; for P1 and Q1 on an axis-aligned
    !> cell this is the flux at the centroid. \(O(1)\) per cell. Fails by
    !> returning a nonzero `ierr`, or `ok = .false.` when the closure cannot
    !> answer at a node.
    subroutine compute_nodal_cell_flux(dm, local_state, cell, flux, ok, ierr)
        implicit none
        !> The solve DM the flux was installed on.
        DM, intent(in) :: dm
        !> A local vector of the solution.
        Vec, intent(in) :: local_state
        PetscInt, intent(in) :: cell
        !> \(\bar{\mathbf q}_l^K\), m/s; overwritten when `ok`.
        real(real64), intent(inout) :: flux(MIXED_DIMENSION)
        logical, intent(inout) :: ok
        PetscErrorCode, intent(inout) :: ierr

        type(type_nodal_cell) :: element
        PetscSection :: section
        real(real64) :: mean(MIXED_DIMENSION), total
        integer(int32) :: q

        call DMGetLocalSection(dm, section, ierr)
        if (ierr /= 0) return
        call read_nodal_cell(dm, section, local_state, cell, element, ok, ierr)
        if (ierr /= 0 .or. .not. ok) return

        mean = 0.0_real64
        total = 0.0_real64
        do q = 1, nodal%num_points
            mean = mean + nodal%weights(q)*nodal_driver(element, q)
            total = total + nodal%weights(q)
        end do
        flux = -element%conductivity*matmul(element%mobility, mean/total)
    end subroutine compute_nodal_cell_flux

    !> Refuse to evaluate the flux of a different mesh.
    subroutine check_installed_dm(dm_handle)
        implicit none
        type(c_ptr), intent(in) :: dm_handle

        if (.not. nodal%active .or. .not. c_associated(dm_handle, nodal%dm)) then
            error stop "The nodal column flux was installed for a different mesh"
        end if
    end subroutine check_installed_dm

end module app_mixed_hydraulic_column_nodal
