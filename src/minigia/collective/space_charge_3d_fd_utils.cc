#include "space_charge_3d_fd_impl.hpp"

/* --------------------------------------------------------------------- */
/*!
  Initialize sequential vectors on each MPI rank
  \param   lctx - local context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
  */

PetscErrorCode init_localvecs(LocalCtx &lctx, GlobalCtx &gctx) {

  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  ierr = VecCreate(PETSC_COMM_SELF, &lctx.seqphi);
  CHKERRQ(ierr);
  ierr = VecSetType(lctx.seqphi, gctx.vectype.c_str());
  CHKERRQ(ierr);
  ierr = VecSetSizes(lctx.seqphi, PETSC_DECIDE, gctx.nsize);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(lctx.seqphi);
  CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)(lctx.seqphi), "seqphi_on_lctx");
  CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_SELF, &lctx.seqrho);
  CHKERRQ(ierr);
  ierr = VecSetType(lctx.seqrho, gctx.vectype.c_str());
  CHKERRQ(ierr);
  ierr = VecSetSizes(lctx.seqrho, PETSC_DECIDE, gctx.nsize);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(lctx.seqrho);
  CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)(lctx.seqrho), "seqrho_on_lctx");
  CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
