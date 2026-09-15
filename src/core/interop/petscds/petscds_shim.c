/*
 * Registration of PetscDS pointwise functions and boundary conditions,
 * reachable from Fortran.
 *
 * PETSc builds no Fortran stubs for the calls that take a pointwise function
 * pointer -- PetscDSSetResidual, PetscDSSetJacobian, DMAddBoundary, the weak
 * form of a boundary, and the DMTS local hooks all exist only as C symbols. This shim passes the pointers straight through, so the weak form
 * itself stays in Fortran where the physics is.
 *
 * The setters are deliberately thin: a Fortran caller hands over a function
 * pointer and nothing else happens on the C side. The one exception is
 * ftcms_ts_use_plex_fem, which wires three PETSc-owned callbacks in one go;
 * routing those through Fortran would add no information and only give the
 * pointers a chance to be mismatched.
 */
#include <petscds.h>
#include <petscdmplex.h>
#include <petscsnes.h>

/* Volume residual of field `f`: f0 is the test-function term, f1 the
 * test-gradient term. Either may be NULL. */
PetscErrorCode ftcms_ds_set_residual(PetscDS ds, PetscInt f, PetscPointFn *f0, PetscPointFn *f1)
{
  PetscFunctionBegin;
  PetscCall(PetscDSSetResidual(ds, f, f0, f1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Volume Jacobian block (f, g). Any of g0..g3 may be NULL:
 *   g0 = d(f0)/du    g1 = d(f0)/d(grad u)
 *   g2 = d(f1)/du    g3 = d(f1)/d(grad u)  */
PetscErrorCode ftcms_ds_set_jacobian(PetscDS ds, PetscInt f, PetscInt g, PetscPointJacFn *g0, PetscPointJacFn *g1, PetscPointJacFn *g2, PetscPointJacFn *g3)
{
  PetscFunctionBegin;
  PetscCall(PetscDSSetJacobian(ds, f, g, g0, g1, g2, g3));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Essential (Dirichlet) condition on the points of `label_name` carrying
 * `value`. The label is looked up here so Fortran never has to hold a DMLabel;
 * "Face Sets" is the name PETSc gives a Gmsh physical group of codimension 1.
 * `comps` may be NULL, which constrains every component. */
PetscErrorCode ftcms_dm_add_boundary_essential(DM dm, const char name[], const char label_name[], PetscInt value, PetscInt field, PetscInt Nc, const PetscInt comps[], void (*func)(void), void *ctx, PetscInt *bd)
{
  DMLabel label;

  PetscFunctionBegin;
  PetscCall(DMGetLabel(dm, label_name, &label));
  PetscCheck(label, PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONG, "The mesh carries no label \"%s\"", label_name);
  PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, name, label, 1, &value, field, Nc, comps, func, NULL, ctx, bd));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Natural (flux, Robin) condition on the points of `label_name` carrying
 * `value`, together with the kernels that evaluate it there.
 *
 * The two go together because PETSc does not consult the DS's own boundary
 * residual for a labelled boundary: it looks the kernels up in the weak form
 * of that boundary, keyed by the label and the value. Registering the boundary
 * without them leaves the surface integral silently at zero. */
PetscErrorCode ftcms_dm_add_boundary_natural(DM dm, const char name[], const char label_name[], PetscInt value, PetscInt field, PetscBdPointFn *f0, PetscBdPointJacFn *g0, PetscInt *bd)
{
  DMLabel       label;
  PetscDS       ds;
  PetscWeakForm wf;

  PetscFunctionBegin;
  PetscCall(DMGetLabel(dm, label_name, &label));
  PetscCheck(label, PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONG, "The mesh carries no label \"%s\"", label_name);
  PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, name, label, 1, &value, field, 0, NULL, NULL, NULL, NULL, bd));

  PetscCall(DMGetDS(dm, &ds));
  PetscCall(PetscDSGetBoundary(ds, *bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
  PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, value, field, 0, 0, f0, 0, NULL));
  if (g0) PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, value, field, field, 0, 0, g0, 0, NULL, 0, NULL, 0, NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Interpolate analytic functions, one per field, into a vector of `dm`.
 *
 * PETSc applies each field's own dual functionals, so an H(div) field receives
 * the facet normal moments of the given vector function rather than its point
 * values. That is the only way to state a flux condition on an RT field, and
 * it is why this cannot be written as a VecSetValues loop on the Fortran side.
 *
 * `funcs` holds one pointer per field of `dm`; `ctxs` may be NULL. */
PetscErrorCode ftcms_dm_project_function(DM dm, PetscReal time, PetscErrorCode (**funcs)(PetscInt, PetscReal, const PetscReal[], PetscInt, PetscScalar *, void *), void **ctxs, InsertMode mode, Vec x)
{
  PetscFunctionBegin;
  PetscCall(DMProjectFunction(dm, time, funcs, ctxs, mode, x));
  PetscFunctionReturn(PETSC_SUCCESS);
}
