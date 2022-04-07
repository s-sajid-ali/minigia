#include "space_charge_3d_fd_utils.hpp"

/* --------------------------------------------------------------------- */
/*!
  Initialize sequential vectors on each MPI rank
  \param   lctx - local context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
  */

PetscErrorCode init_localvecs(LocalCtx &lctx, GlobalCtx &gctx) {
  PetscFunctionBeginUser;

  PetscCall(VecCreate(PETSC_COMM_SELF, &lctx.seqphi));
  PetscCall(VecSetType(lctx.seqphi, gctx.vectype.c_str()));
  PetscCall(VecSetSizes(lctx.seqphi, PETSC_DECIDE, gctx.nsize));
  PetscCall(VecSetFromOptions(lctx.seqphi));
  PetscCall(PetscObjectSetName((PetscObject)(lctx.seqphi), "seqphi_on_lctx"));

  PetscCall(VecCreate(PETSC_COMM_SELF, &lctx.seqrho));
  PetscCall(VecSetType(lctx.seqrho, gctx.vectype.c_str()));
  PetscCall(VecSetSizes(lctx.seqrho, PETSC_DECIDE, gctx.nsize));
  PetscCall(VecSetFromOptions(lctx.seqrho));
  PetscCall(PetscObjectSetName((PetscObject)(lctx.seqrho), "seqrho_on_lctx"));

  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  finalize by destroying data structures
  \param   lctx - local context
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
*/
PetscErrorCode finalize(LocalCtx &lctx, SubcommCtx &sctx, GlobalCtx &gctx) {
  PetscFunctionBeginUser;

  /* Destroy local vectors */
  PetscCall(VecDestroy(&(lctx.seqphi)));
  PetscCall(VecDestroy(&(lctx.seqrho)));

  PetscFunctionReturn(0);
}
