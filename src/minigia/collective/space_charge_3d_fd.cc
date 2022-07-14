#include "space_charge_3d_fd.hpp"
#include "deposit.hpp"
#include "space_charge_3d_fd_alias.hpp"
#include "space_charge_3d_fd_utils.hpp"

#include <minigia/utils/simple_timer.hpp>

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
      domain(ops.shape, {1.0, 1.0, 1.0}), use_fixed_domain(false),
      allocated(false) {}

// Destructor
Space_charge_3d_fd::~Space_charge_3d_fd() {

  if (allocated) {
    destroy_sc3d_fd();
  }
}

void Space_charge_3d_fd::get_local_charge_density(Bunch const &bunch) {
  scoped_simple_timer timer("sc3d_fd_local_rho");
  deposit_charge_rectangular_3d_kokkos_scatter_view(
      lctx.seqrho_view, domain, domain.get_grid_shape(), bunch);
}

void Space_charge_3d_fd::apply_impl(Bunch_simulator &sim, double time_step,
                                    Logger &logger) {
  logger << "    Space charge 3d open hockney\n";

  scoped_simple_timer timer("sc3d_total");

  // count number of bunches
  int num_bunches_in_bunch_sim = 0;
  int num_bunches_in_train_0 = 0;
  for (size_t t = 0; t < 2; ++t) {
    for (size_t b = 0; b < sim[t].get_bunch_array_size(); ++b) {
      num_bunches_in_bunch_sim += 1;
      if (t == 0)
        num_bunches_in_train_0 += 1;
    }
  }
  if (num_bunches_in_bunch_sim != num_bunches_in_train_0) {
    throw std::runtime_error(
        "sc3d-fd only works on a single bunch train, work is ongoing \
		    to make it work on multiple bunche trains!");
  }

  // construct the workspace for a new bunch simulator
  if (bunch_sim_id != sim.id()) {
    bunch_sim_id = sim.id();
    // Assumption: When an MPI rank has more than 1 bunch (within the same
    // train), all the bunches on the rank have the same communicator and the
    // same distribution. Reason: when constructing a bunch train, one can have
    // one the two scenarios:
    // [1] number of bunches > number of MPI ranks, each rank always has only 1
    // bunch [2] number of bunches < number of MPI ranks, we require number of
    // bunches to be divisible by number of  MPI ranks and each bunch has a MPI
    // communicator of size 1. Update this bit when enabling bunch sim with two
    // bunch trains!
    allocate_sc3d_fd(sim[0][0]);
    allocated = true;
  }

  // apply to bunches
  for (size_t t = 0; t < 2; ++t) {
    for (size_t b = 0; b < sim[t].get_bunch_array_size(); ++b) {
      apply_bunch(sim[t][b], time_step, logger);
    }
  }
}

void Space_charge_3d_fd::apply_bunch(Bunch &bunch, double time_step,
                                     Logger &logger) {
  // charge density
  get_local_charge_density(bunch); // [C/m^3]
}

// update_domain
void Space_charge_3d_fd::update_domain(Bunch const &bunch) {

  auto spatial_mean_stddev =
      Core_diagnostics::calculate_spatial_mean_std(bunch);
  auto mean_x = spatial_mean_stddev(0);
  auto mean_y = spatial_mean_stddev(1);
  auto mean_z = spatial_mean_stddev(2);
  auto stddev_x = spatial_mean_stddev(3);
  auto stddev_y = spatial_mean_stddev(4);
  auto stddev_z = spatial_mean_stddev(5);

  const double tiny = 1.0e-10;

  if ((stddev_x < tiny) && (stddev_y < tiny) && (stddev_z < tiny)) {
    throw std::runtime_error(
        "Space_charge_3d_open_hockney_eigen::update_domain: "
        "all three spatial dimensions have neglible extent");
  }

  std::array<double, 3> offset{mean_x, mean_y, mean_z};

  std::array<double, 3> size{
      options.n_sigma *
          get_smallest_non_tiny(stddev_x, stddev_y, stddev_z, tiny),
      options.n_sigma *
          get_smallest_non_tiny(stddev_y, stddev_x, stddev_z, tiny),
      options.n_sigma *
          get_smallest_non_tiny(stddev_z, stddev_x, stddev_y, tiny)};

  domain = Rectangular_grid_domain(options.shape, size, offset, false);
}

PetscErrorCode Space_charge_3d_fd::allocate_sc3d_fd(const Bunch &bunch) {
  PetscFunctionBeginUser;

  /* size of seqphi/seqrho vectors/views is size of domain! */
  gctx.nsize = options.shape[0] * options.shape[1] * options.shape[2];

  /* store MPI communicator of bunch in gctx */
  PetscCall(
      PetscCommDuplicate(MPI_Comm(bunch.get_comm()), &gctx.bunch_comm, NULL););
  PetscCallMPI(MPI_Comm_rank(gctx.bunch_comm, &gctx.global_rank));
  PetscCallMPI(MPI_Comm_size(gctx.bunch_comm, &gctx.global_size));

  /* Initialize task subcomms, display task-subcomm details */
  PetscCall(init_solversubcomms(sctx, gctx));

  /* Local rho and phi vectors on each MPI rank */
  PetscCall(init_localvecs(lctx, gctx));

  /* rho and phi vectors on each subcomm */
  PetscCall(init_subcommvecs(sctx, gctx));

  /* create global aliases of local vectors */
  PetscCall(init_global_local_aliases(lctx, gctx));

  /* create global aliases of subcomm vectors */
  PetscCall(init_global_subcomm_aliases(sctx, gctx));

  /* create subcomm aliases of local vectors */
  PetscCall(init_subcomm_local_aliases(lctx, sctx, gctx));

  /* Initialize global (alias of local) to subcomm scatters */
  PetscCall(init_global_subcomm_scatters(sctx, gctx));

  /* Initialize subcomm (alias of local) to local scatters */
  PetscCall(init_subcomm_local_scatters(lctx, sctx, gctx));

  PetscFunctionReturn(0);
}

PetscErrorCode Space_charge_3d_fd::destroy_sc3d_fd() {
  PetscFunctionBeginUser;

  PetscCall(finalize(lctx, sctx, gctx));

  PetscFunctionReturn(0);
}
