#include <catch2/catch_test_macros.hpp>
#include <minigia/collective/space_charge_3d_fd.hpp>

#include "rod_bunch.h"

#include <petsc.h>

TEST_CASE("real_apply_full_lowgamma", "[Rod_bunch]") {
  PetscErrorCode ierr;
  auto logger = Logger(0, LoggerV::DEBUG);
  auto simlogger = Logger(0, LoggerV::INFO_STEP);

  const int gridx = 64;
  const int gridy = 64;
  const int gridz = 64;

  const double time_fraction = 1.0;
  const double step_length = 0.1;

  // bunch
  Rod_bunch_fixture_lowgamma fixture;

  auto &bunch = fixture.bsim.get_bunch();
  auto const &ref = bunch.get_reference_particle();
  auto parts = bunch.get_host_particles();

  const double beta = ref.get_beta();
  const double gamma = ref.get_gamma();
  const double betagamma = beta * gamma;

  const double time_step = step_length / (beta * pconstants::c);
  // const double bunchlen = bunch.get_longitudinal_boundary().second;
  const double bunchlen = 0.1;

  // check out particles before print
  bunch.checkout_particles();

  // print intital coordinates
  logger << "real_apply_full_lowgamma first four particles (x y z):" << '\n';
  for (int k = 0; k < 4; ++k) {
    logger << k << ": " << parts(k, 0) << ", " << parts(k, 2) << ", "
           << parts(k, 4) << '\n';
  }

  logger << "last four particles (x y z):" << '\n';
  for (int k = bunch.get_local_num() - 4; k < bunch.get_local_num(); ++k) {
    logger << k << ": " << parts(k, 0) << ", " << parts(k, 2) << ", "
           << parts(k, 4) << '\n';
  }

  logger << '\n';

  // space charge operator
  auto sc_ops = Space_charge_3d_fd_options(gridx, gridy, gridz);

  sc_ops.comm_group_size = 1;

  auto sc = Space_charge_3d_fd(sc_ops);

  // set domain
  std::array<double, 3> offset = {0, 0, 0};
  std::array<double, 3> size = {parts(0, 0) * 4, parts(0, 0) * 4,
                                bunchlen / beta};

  // apply space charge operator
  sc.apply(fixture.bsim, time_step, simlogger);

  // check out particles
  bunch.checkout_particles();

  // print
  logger << "after sc::apply : bunch.local_particles(0, 0): " << parts(0, 0)
         << '\n';
}
