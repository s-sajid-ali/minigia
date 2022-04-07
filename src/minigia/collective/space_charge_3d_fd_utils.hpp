#include "space_charge_3d_fd_impl.hpp"

PetscErrorCode init_localvecs(LocalCtx &lctx, GlobalCtx &gctx);

PetscErrorCode finalize(LocalCtx &lctx, SubcommCtx &sctx, GlobalCtx &gctx);
