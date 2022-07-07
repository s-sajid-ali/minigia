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

  lctx.seqphi_view = karray1d_dev("seqphi", gctx.nsize);
  lctx.seqrho_view = karray1d_dev("seqphi", gctx.nsize);

  PetscCall(gctx.VecCreate_type_WithArray(PETSC_COMM_SELF, 1, gctx.nsize,
                                          PETSC_DECIDE, lctx.seqphi_view.data(),
                                          &lctx.seqphi));
  PetscCall(PetscObjectSetName((PetscObject)(lctx.seqphi), "seqphi_on_lctx"));

  PetscCall(gctx.VecCreate_type_WithArray(PETSC_COMM_SELF, 1, gctx.nsize,
                                          PETSC_DECIDE, lctx.seqrho_view.data(),
                                          &lctx.seqrho));
  PetscCall(PetscObjectSetName((PetscObject)(lctx.seqrho), "seqrho_on_lctx"));

  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  Initialize solver subcomms
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
*/
PetscErrorCode init_solversubcomms(SubcommCtx &sctx, GlobalCtx &gctx) {

  PetscFunctionBeginUser;

  PetscCall(PetscSubcommCreate(PETSC_COMM_WORLD, &(sctx.solverpsubcomm)));
  PetscCall(PetscSubcommSetNumber(sctx.solverpsubcomm, gctx.nsubcomms));
  PetscCall(PetscSubcommSetType(sctx.solverpsubcomm, PETSC_SUBCOMM_CONTIGUOUS));

  PetscCall(PetscSubcommGetChild(sctx.solverpsubcomm, &(sctx.solversubcomm)));
  PetscCall(MPIU_Allreduce(&(gctx.global_rank), &(sctx.solversubcommid), 1,
                           MPIU_INT, MPI_MAX, sctx.solversubcomm));
  PetscCallMPI(MPI_Comm_rank(sctx.solversubcomm, &(sctx.solver_rank)));
  PetscCallMPI(MPI_Comm_size(sctx.solversubcomm, &(sctx.solver_size)));

  /* collect the subcommids on each MPI rank & remove duplicates */
  gctx.sids.resize(gctx.global_size);
  PetscCallMPI(MPI_Allgather(&sctx.solversubcommid, 1, MPI_INT,
                             gctx.sids.data(), 1, MPI_INT, PETSC_COMM_WORLD));
  std::sort(gctx.sids.begin(), gctx.sids.end());
  gctx.sids.erase(std::unique(gctx.sids.begin(), gctx.sids.end()),
                  gctx.sids.end());

  if (gctx.debug) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                          "\nsolver-subcomms have been created!\n"));
    PetscCall(PetscSubcommView(sctx.solverpsubcomm, PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n"));
    PetscCall(PetscSynchronizedPrintf(
        PETSC_COMM_WORLD, "size of gct.sids vector on global-rank %d is : %d\n",
        gctx.global_rank, static_cast<int>(gctx.sids.size())));
    PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n"));
  }

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
