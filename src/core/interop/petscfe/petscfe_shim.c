/*
 * Basis tabulation of a PetscFE, reachable from Fortran.
 *
 * PetscFE can tabulate its basis at any set of reference points, but the
 * generated Fortran binding hands back a PetscTabulation whose Fortran
 * descriptor is never built, so reading it from Fortran faults. This shim does
 * the call on the C side and copies the values into arrays the caller owns,
 * which is what the interop layer is for.
 */
#include <petscfe.h>
#include <petscds.h>
#include <petscdmplex.h>
#include <string.h>

/*
 * Tabulate at `num_points` reference points.
 *
 *   points     : num_points * dim, point-major
 *   values     : num_points * num_basis            (may be NULL)
 *   gradients  : num_points * num_basis * dim      (may be NULL)
 *
 * `values` and `gradients` are filled only when non-NULL; ask for gradients and
 * the tabulation is built to first order, otherwise to zeroth.
 */
PetscErrorCode ftcms_petscfe_tabulate(PetscFE fe, PetscInt num_points, const PetscReal points[],
                                      PetscReal values[], PetscReal gradients[])
{
  PetscTabulation tab;
  PetscInt        K = gradients ? 1 : 0;
  PetscInt        Nb, Nc, cdim, count;

  PetscFunctionBegin;
  PetscCall(PetscFECreateTabulation(fe, 1, num_points, points, K, &tab));

  Nb   = tab->Nb;
  Nc   = tab->Nc;
  cdim = tab->cdim;

  if (values) {
    count = num_points * Nb * Nc;
    memcpy(values, tab->T[0], (size_t)count * sizeof(PetscReal));
  }
  if (gradients) {
    count = num_points * Nb * Nc * cdim;
    memcpy(gradients, tab->T[1], (size_t)count * sizeof(PetscReal));
  }

  PetscCall(PetscTabulationDestroy(&tab));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
 * Size of one field's discretisation: `num_basis` basis functions carrying
 * `num_components` components each. One degree of freedom per basis function,
 * so `num_basis` is also this field's share of a cell closure.
 */
PetscErrorCode ftcms_petscfe_field_size(DM dm, PetscInt field, PetscInt *num_basis, PetscInt *num_components)
{
  PetscDS ds;
  PetscFE fe;

  PetscFunctionBegin;
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(PetscDSGetDiscretization(ds, field, (PetscObject *)&fe));
  PetscCall(PetscFEGetDimension(fe, num_basis));
  PetscCall(PetscFEGetNumComponents(fe, num_components));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
 * The physical value and divergence of one field inside one cell.
 *
 * This is a reconstruction, not an assembly: it reads the cell's own degrees of
 * freedom and maps the reference basis to physical space with PETSc's own
 * pushforward, which for an H(div) field is the contravariant Piola map. The
 * conservation audit needs a value that does not come from the same element
 * loop it is checking, and it needs the map to be the one the assembly used,
 * not a second implementation of it.
 *
 *   ref_points : num_points * dim, point-major, on the reference cell
 *   values     : num_points * num_components
 *   divergence : num_points                       (may be NULL)
 *
 * The closure carries PETSc's orientation flips, so a facet degree of freedom
 * shared by two cells is read with each cell's own sign.
 */
PetscErrorCode ftcms_petscfe_evaluate_cell(DM dm, Vec local_x, PetscInt field, PetscInt cell,
                                           PetscInt num_points, const PetscReal ref_points[],
                                           PetscScalar values[], PetscScalar divergence[])
{
  PetscDS         ds;
  PetscFE         fe;
  PetscSection    section;
  PetscTabulation tab;
  PetscFEGeom     fegeom;
  PetscScalar    *closure = NULL;
  PetscScalar     u[3], u_x[9];
  PetscReal       v0[3], J[9], invJ[9], detJ;
  PetscInt        num_fields, Nb, Nc, closure_offset = 0, csize, dim, cdim, f, p, b, c, d;

  PetscFunctionBegin;
  PetscCall(DMGetDimension(dm, &dim));
  PetscCall(DMGetCoordinateDim(dm, &cdim));
  PetscCheck(dim == cdim, PetscObjectComm((PetscObject)dm), PETSC_ERR_SUP,
             "The reconstruction needs intrinsic coordinates: dim %" PetscInt_FMT " but coordinate dim %" PetscInt_FMT, dim, cdim);

  PetscCall(DMGetDS(dm, &ds));
  PetscCall(PetscDSGetNumFields(ds, &num_fields));
  PetscCheck(field >= 0 && field < num_fields, PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_OUTOFRANGE,
             "Field %" PetscInt_FMT " is not one of the %" PetscInt_FMT " fields", field, num_fields);

  /* A closure is field-major, so this field starts past every earlier field. */
  for (f = 0; f < field; ++f) {
    PetscFE  fe_prev;
    PetscInt Nb_prev;

    PetscCall(PetscDSGetDiscretization(ds, f, (PetscObject *)&fe_prev));
    PetscCall(PetscFEGetDimension(fe_prev, &Nb_prev));
    closure_offset += Nb_prev;
  }
  PetscCall(PetscDSGetDiscretization(ds, field, (PetscObject *)&fe));
  PetscCall(PetscFEGetDimension(fe, &Nb));
  PetscCall(PetscFEGetNumComponents(fe, &Nc));
  PetscCheck(Nc <= 3, PetscObjectComm((PetscObject)dm), PETSC_ERR_SUP, "At most three components, not %" PetscInt_FMT, Nc);

  PetscCall(DMGetLocalSection(dm, &section));
  PetscCall(DMPlexVecGetClosure(dm, section, local_x, cell, &csize, &closure));
  PetscCall(PetscFECreateTabulation(fe, 1, num_points, ref_points, 1, &tab));
  PetscCall(DMPlexComputeCellGeometryAffineFEM(dm, cell, v0, J, invJ, &detJ));

  PetscCall(PetscMemzero(&fegeom, sizeof(fegeom)));
  fegeom.mode       = PETSC_FEGEOM_BASIC;
  fegeom.isAffine   = PETSC_TRUE;
  fegeom.dim        = dim;
  fegeom.dimEmbed   = cdim;
  fegeom.numCells   = 1;
  fegeom.numPoints  = 1;
  fegeom.v          = v0;
  fegeom.J          = J;
  fegeom.invJ       = invJ;
  fegeom.detJ       = &detJ;

  for (p = 0; p < num_points; ++p) {
    for (c = 0; c < Nc; ++c) {
      u[c] = 0.0;
      for (d = 0; d < cdim; ++d) u_x[c * cdim + d] = 0.0;
    }
    for (b = 0; b < Nb; ++b) {
      const PetscScalar coefficient = closure[closure_offset + b];

      for (c = 0; c < Nc; ++c) {
        u[c] += coefficient * tab->T[0][(p * Nb + b) * Nc + c];
        for (d = 0; d < cdim; ++d) u_x[c * cdim + d] += coefficient * tab->T[1][((p * Nb + b) * Nc + c) * cdim + d];
      }
    }
    PetscCall(PetscFEPushforward(fe, &fegeom, 1, u));
    PetscCall(PetscFEPushforwardGradient(fe, &fegeom, 1, u_x));

    for (c = 0; c < Nc; ++c) values[p * Nc + c] = u[c];
    if (divergence) {
      PetscScalar div = 0.0;

      for (c = 0; c < Nc && c < cdim; ++c) div += u_x[c * cdim + c];
      divergence[p] = div;
    }
  }

  PetscCall(PetscTabulationDestroy(&tab));
  PetscCall(DMPlexVecRestoreClosure(dm, section, local_x, cell, &csize, &closure));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
 * Reference coordinates of physical points inside one affine cell.
 *
 * The inverse of the affine map the tabulation above is evaluated on, so a
 * caller that knows where a facet quadrature point lies in the mesh can ask
 * for the field there.
 *
 *   points     : num_points * dim, point-major, physical
 *   ref_points : num_points * dim, point-major, on the reference cell
 */
PetscErrorCode ftcms_petscfe_reference_coordinates(DM dm, PetscInt cell, PetscInt num_points,
                                                   const PetscReal points[], PetscReal ref_points[])
{
  PetscReal v0[3], J[9], invJ[9], detJ;
  PetscInt  dim, p, i, j;

  PetscFunctionBegin;
  PetscCall(DMGetDimension(dm, &dim));
  PetscCall(DMPlexComputeCellGeometryAffineFEM(dm, cell, v0, J, invJ, &detJ));
  for (p = 0; p < num_points; ++p) {
    for (i = 0; i < dim; ++i) {
      PetscReal xi = -1.0;

      for (j = 0; j < dim; ++j) xi += invJ[i * dim + j] * (points[p * dim + j] - v0[j]);
      ref_points[p * dim + i] = xi;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
