#ifndef SPACE_CHARGE_3D_FD_UTILS_H_
#define SPACE_CHARGE_3D_FD_UTILS_H_

#include "space_charge_3d_fd_impl.hpp"

PetscErrorCode init_solver_subcomms(SubcommCtx &sctx, GlobalCtx &gctx);

PetscErrorCode init_local_vecs(LocalCtx &lctx, GlobalCtx &gctx);
PetscErrorCode init_subcomm_vecs(SubcommCtx &sctx, GlobalCtx &gctx);

PetscErrorCode init_subcomm_mat(SubcommCtx &sctx, GlobalCtx &gctx);

PetscErrorCode compute_mat(SubcommCtx &sctx, GlobalCtx &gctx);

PetscErrorCode init_global_subcomm_scatters(SubcommCtx &sctx, GlobalCtx &gctx);
PetscErrorCode init_subcomm_local_scatters(LocalCtx &lctx, SubcommCtx &sctx,
                                           GlobalCtx &gctx);

PetscErrorCode finalize(LocalCtx &lctx, SubcommCtx &sctx, GlobalCtx &gctx);

#endif
