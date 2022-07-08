#include "space_charge_3d_fd_alias.hpp"

/* --------------------------------------------------------------------- */
/*!
  Initialize global aliases of local vectors
  \param   lctx - local context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
  */
PetscErrorCode init_global_local_aliases(LocalCtx &lctx, GlobalCtx &gctx) {

  PetscInt phi_localsize; /* size of vector on this MPI rank */
  PetscInt rho_localsize; /* size of vector on this MPI rank */

  PetscScalar const
      *d_rho_val; /* read-only pointer to vector's contents on device */
  PetscScalar const
      *d_phi_val; /* read-only pointer to vector's contents on device */

  PetscMemType mtype_phi; /* memory type for phi vector */
  PetscMemType mtype_rho; /* memory type for rho vector */

  PetscFunctionBeginUser;

  /* Get size of local vector */
  PetscCall(VecGetLocalSize(lctx.seqphi, &phi_localsize));
  PetscCall(VecGetLocalSize(lctx.seqrho, &rho_localsize));

  /* Get access to local vector */
  PetscCall(VecGetArrayReadAndMemType(lctx.seqphi, &d_phi_val, &mtype_phi));
  PetscCall(VecGetArrayReadAndMemType(lctx.seqrho, &d_rho_val, &mtype_rho));

  /* consistency check */
  if (gctx.debug) {
    if (mtype_phi != mtype_rho)
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_NOTSAMETYPE,
              "phi and rho vectors don't have the memory type!");
    if (!gctx.VecCreate_type_WithArray)
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_BADPTR,
              "the function pointer in gctx is invalid!");
  }

  PetscCall(gctx.VecCreate_type_WithArray(PETSC_COMM_WORLD, 1, phi_localsize,
                                          PETSC_DECIDE, d_phi_val,
                                          &gctx.phi_global_local));
  PetscCall(gctx.VecCreate_type_WithArray(PETSC_COMM_WORLD, 1, rho_localsize,
                                          PETSC_DECIDE, d_rho_val,
                                          &gctx.rho_global_local));

  PetscCall(PetscObjectSetName((PetscObject)(gctx.phi_global_local),
                               "phi_global_local_on_gctx"));
  PetscCall(PetscObjectSetName((PetscObject)(gctx.rho_global_local),
                               "rho_global_local_on_gctx"));

  /* Restore local vector arrays */
  PetscCall(VecRestoreArrayReadAndMemType(lctx.seqphi, &d_phi_val));
  PetscCall(VecRestoreArrayReadAndMemType(lctx.seqrho, &d_rho_val));

  PetscFunctionReturn(0);
}
