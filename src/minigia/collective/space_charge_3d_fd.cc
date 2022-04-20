#include "space_charge_3d_fd.hpp"
#include "space_charge_3d_fd_utils.hpp"

namespace {
double get_smallest_non_tiny(double val, double other1, double other2,
                             double tiny) {
  double retval;
  if (val > tiny) {
    retval = val;
  } else {
    if ((other1 > tiny) && (other2 > tiny)) {
      retval = std::min(other1, other2);
    } else {
      retval = std::max(other1, other2);
    }
  }

  return retval;
}

} // namespace

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

// update_domain
void Space_charge_3d_fd::update_domain(Bunch const &bunch) {

  auto mean_stddev = Core_diagnostics::calculate_spatial_mean_std(bunch);

  const double tiny = 1.0e-10;

  const auto ix = Bunch::x;
  const auto iy = Bunch::y;
  const auto iz = Bunch::z;

  if ((mean_stddev[3] < tiny) && (mean_stddev[4] < tiny) &&
      (mean_stddev[5] < tiny)) {
    throw std::runtime_error("space_charge_3d_fd::update_domain: all three "
                             "spatial dimensions have neglible extent");
  }

  std::array<double, 3> size{
      get_smallest_non_tiny(mean_stddev[3], mean_stddev[4], mean_stddev[5],
                            tiny) *
          options.n_sigma,
      get_smallest_non_tiny(mean_stddev[4], mean_stddev[3], mean_stddev[5],
                            tiny) *
          options.n_sigma,
      get_smallest_non_tiny(mean_stddev[5], mean_stddev[3], mean_stddev[4],
                            tiny) *
          options.n_sigma};

  std::array<double, 3> offset{mean_stddev[0], mean_stddev[1], mean_stddev[2]};
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
