#include "space_charge_3d_fd.hpp"
#include "space_charge_3d_fd_utils.hpp"

// Constructor
Space_charge_3d_fd::Space_charge_3d_fd(Space_charge_3d_fd_options const &ops)
    : Collective_operator("sc_3d_fd", 1.0), options(ops), bunch_sim_id(),
      allocated(false), use_fixed_domain(false) {

  allocate_sc3d_fd(ops);
  allocated = true;
}

// Destructor
Space_charge_3d_fd::~Space_charge_3d_fd() {

  if (allocated) {
    destroy_sc3d_fd();
  }
}

PetscErrorCode
Space_charge_3d_fd::allocate_sc3d_fd(Space_charge_3d_fd_options const &ops) {
  PetscFunctionBeginUser;

  PetscCall(init_localvecs(lctx, gctx));

  PetscFunctionReturn(0);
}

PetscErrorCode Space_charge_3d_fd::destroy_sc3d_fd() {
  PetscFunctionBeginUser;

  PetscCall(finalize(lctx, sctx, gctx));

  PetscFunctionReturn(0);
}
