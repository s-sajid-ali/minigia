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
  Initialize vectors on each solver-subcommunicator
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
  */
PetscErrorCode init_subcommvecs(SubcommCtx &sctx, GlobalCtx &gctx) {
  PetscInt size;

  PetscFunctionBeginUser;
  PetscCall(VecCreate(sctx.solversubcomm, &sctx.phi_subcomm));
  PetscCall(VecSetType(sctx.phi_subcomm, gctx.vectype.c_str()));
  PetscCall(VecSetSizes(sctx.phi_subcomm, PETSC_DECIDE, gctx.nsize));
  PetscCall(VecSetFromOptions(sctx.phi_subcomm));
  PetscCall(PetscObjectSetName((PetscObject)(sctx.phi_subcomm),
                               "phi_subcomm_on_sctx"));

  PetscCall(VecCreate(sctx.solversubcomm, &sctx.rho_subcomm));
  PetscCall(VecSetType(sctx.rho_subcomm, gctx.vectype.c_str()));
  PetscCall(VecSetSizes(sctx.rho_subcomm, PETSC_DECIDE, gctx.nsize));
  PetscCall(VecSetFromOptions(sctx.rho_subcomm));
  PetscCall(PetscObjectSetName((PetscObject)(sctx.rho_subcomm),
                               "rho_subcomm_on_sctx"));

  if (gctx.debug) {
    PetscCall(VecGetSize(sctx.rho_subcomm, &size));
    PetscCall(
        PetscPrintf(PetscObjectComm((PetscObject)sctx.rho_subcomm),
                    "Hi there from solver-subcomm with id %d, the size of "
                    "subcomm vec here is %d\n",
                    sctx.solversubcommid, size));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n"));
  }

  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  Initialize global (alias of local) to subcomm scatters
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
*/
PetscErrorCode init_global_subcomm_scatters(SubcommCtx &sctx, GlobalCtx &gctx) {

  PetscInt start;
  PetscInt localsize;

  PetscFunctionBeginUser;

  /* localsize should be gctx->nsize, but adding a check doesn't hurt */
  PetscCall(VecGetLocalSize(gctx.rho_global_local, &localsize));
  PetscCall(VecGetOwnershipRange(gctx.rho_global_local, &start, NULL));

  PetscCall(ISCreateStride(PETSC_COMM_SELF, localsize, start, 1,
                           &gctx.ix_scat_glocal_to_subcomms));
  PetscCall(ISCreateStride(PETSC_COMM_SELF, localsize, 0, 1,
                           &gctx.iy_scat_glocal_to_subcomms));

  PetscCall(PetscObjectSetName((PetscObject)(gctx.ix_scat_glocal_to_subcomms),
                               "ix_scat_glocal_to_subcomms"));
  PetscCall(PetscObjectSetName((PetscObject)(gctx.iy_scat_glocal_to_subcomms),
                               "iy_scat_glocal_to_subcomms"));

  /* resize to ensure we only create as many scatters
     as the number of subcomms ! */
  gctx.scat_glocal_to_subcomms.resize(gctx.nsubcomms);

  for (PetscInt i = 0; i < gctx.nsubcomms; i++) {
    PetscCall(VecScatterCreate(
        gctx.rho_global_local, gctx.ix_scat_glocal_to_subcomms,
        gctx.rho_global_subcomm[i], gctx.iy_scat_glocal_to_subcomms,
        &(gctx.scat_glocal_to_subcomms[i])));
  }

  if (gctx.debug) {
    for (PetscInt i = 0; i < gctx.nsubcomms; i++) {
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "global-to-subcomm scatter\n"));
      PetscCall(VecScatterView(gctx.scat_glocal_to_subcomms[i],
                               PETSC_VIEWER_STDOUT_WORLD));
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n"));
  }

  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  Initialize subcomm (alias of local) to local scatters
  \param   lctx - local context
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
*/
PetscErrorCode init_subcomm_local_scatters(LocalCtx &lctx, SubcommCtx &sctx,
                                           GlobalCtx &gctx) {

  PetscInt start;
  PetscInt localsize;

  PetscFunctionBeginUser;

  /* localsize should be gctx->nsize, but adding a check doesn't hurt */
  PetscCall(VecGetLocalSize(sctx.phi_subcomm_local, &localsize));
  PetscCall(VecGetOwnershipRange(sctx.phi_subcomm_local, &start, NULL));

  PetscCall(ISCreateStride(PETSC_COMM_SELF, localsize, start, 1,
                           &sctx.ix_scat_subcomms_to_local));
  PetscCall(ISCreateStride(PETSC_COMM_SELF, localsize, 0, 1,
                           &sctx.iy_scat_subcomms_to_local));

  /* This scatter will always be run as a SCATTER_REVERSE,
     perhaps the naming terminilogy may be updated to prevent
     any confusion in the future.*/
  PetscCall(VecScatterCreate(
      sctx.phi_subcomm_local, sctx.ix_scat_subcomms_to_local, sctx.phi_subcomm,
      sctx.iy_scat_subcomms_to_local, &sctx.scat_subcomm_to_local));

  if (gctx.debug) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "subcomm-to-local scatter\n"));
    PetscCall(VecScatterView(
        sctx.scat_subcomm_to_local,
        PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)sctx.phi_subcomm))));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n"));
    PetscCall(PetscBarrier(NULL));
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

  /* Destroy global aliases of local vectors */
  PetscCall(VecDestroy(&(gctx.phi_global_local)));
  PetscCall(VecDestroy(&(gctx.rho_global_local)));

  /* Destroy global aliases of subcomm vectors */
  for (PetscInt i = 0; i < gctx.nsubcomms; i++) {
    PetscCall(VecDestroy(&(gctx.phi_global_subcomm[i])));
    PetscCall(VecDestroy(&(gctx.rho_global_subcomm[i])));
  }

  /* Destroy subcomm vectors */
  PetscCall(VecDestroy(&(sctx.phi_subcomm)));
  PetscCall(VecDestroy(&(sctx.rho_subcomm)));

  /* Destroy global (alias of local) to subcomm scatters */
  for (PetscInt i = 0; i < gctx.nsubcomms; i++) {
    PetscCall(VecScatterDestroy(&(gctx.scat_glocal_to_subcomms[i])));
  }
  PetscCall(ISDestroy(&gctx.ix_scat_glocal_to_subcomms));
  PetscCall(ISDestroy(&gctx.iy_scat_glocal_to_subcomms));

  /* Destroy subcomm aliases of local vectors */
  PetscCall(VecDestroy(&sctx.phi_subcomm_local));
  PetscCall(VecDestroy(&sctx.rho_subcomm_local));

  /* Destroy subcomm vectors */
  PetscCall(VecDestroy(&(lctx.seqphi)));
  PetscCall(VecDestroy(&(lctx.seqrho)));

  /* Destroy subcomm (alias of local) to local scatters */
  PetscCall(VecScatterDestroy(&sctx.scat_subcomm_to_local));
  PetscCall(ISDestroy(&sctx.ix_scat_subcomms_to_local));
  PetscCall(ISDestroy(&sctx.iy_scat_subcomms_to_local));

  /* Destroy subcomms */
  PetscCall(PetscSubcommDestroy(&(sctx.solverpsubcomm)));
  PetscFunctionReturn(0);
}
