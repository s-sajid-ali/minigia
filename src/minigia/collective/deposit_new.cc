#include <minigia/bunch/bunch.hpp>
#include <minigia/bunch/core_diagnostics.hpp>
#include <minigia/foundation/physical_constants.hpp>

#include "deposit_new.hpp"

#include <Kokkos_StdAlgorithms.hpp>

namespace KE = Kokkos::Experimental;

void bunch_get_idx(Bunch &bunch, Rectangular_grid_domain &domain,
                   karray1i_dev &idx) {

  std::array<double, 3> h = domain.get_cell_size();
  std::array<double, 3> l = domain.get_left();
  std::array<int, 3> const &dims = domain.get_grid_shape();

  int dx = dims[0];
  int dy = dims[1];
  int dz = dims[2];

  auto parts = bunch.get_local_particles();
  auto masks = bunch.get_local_particle_masks();
  int nparts = bunch.size();
  double ihx = 1 / h[0];
  double ihy = 1 / h[1];
  double ihz = 1 / h[1];
  double lx = l[0];
  double ly = l[1];
  double lz = l[2];

  Kokkos::parallel_for(
      "bunch_sort_idx", nparts, KOKKOS_LAMBDA(const int &i) {
        if (masks(i)) {
          int ix, iy, iz;

          ix = KE::floor((parts(i, 0) - lx) * ihx - 0.5);
          iy = KE::floor((parts(i, 2) - ly) * ihy - 0.5);
          iz = KE::floor((parts(i, 4) - lz) * ihz - 0.5);

          if (ix > 0 && ix < dx - 1 && iy > 0 && iy < dy - 1 && iz > 0 &&
              iz < dz - 1) {
            idx(i) = iz * dx * dy + iy * dx + ix;
          } else {
            idx(i) = dx * dy * dz * 10;
          }
        } else {
          idx(i) = dx * dy * dz * 100;
        }
      });
}

/* The 8 cells that contribute to a deposit location are:
   x/y/z;
   x-1/y/z;
   x/y-1/z;
   x-1/y-1/z;
   ...(repeat the above for z-1)
   this is because we use left bottom index as the cell id for sorting
   */
void deposit_charge_rectangular_3d_sorted_v1(
    karray1d_dev &rho_dev,
    Kokkos::BinSort<karray1i_dev, Kokkos::BinOp1D<karray1i_dev>> &sorter,
    Rectangular_grid_domain &domain, std::array<int, 3> const &dims,
    Bunch const &bunch) {

  typedef Kokkos::TeamPolicy<> team_policy;
  typedef Kokkos::TeamPolicy<>::member_type member_type;

  auto permute_idx = sorter.get_permute_vector();
  double charge_val = 0.0;

  Kokkos::parallel_for(
      "deposit_outer_loop",
      team_policy(dims[0] * dims[1] * dims[2], Kokkos::AUTO),

      KOKKOS_LAMBDA(const member_type &teamMember) {
        for (int cell_idx = 0; cell_idx < 8; cell_idx++) {
          int deposit_loc = teamMember.league_rank();
          int ix, iy, iz;
          int cell_loc;

          iz = KE::trunc(deposit_loc / (dims[0] * dims[1]));
          iy = KE::trunc((deposit_loc - iz * dims[0] * dims[1]) / dims[0]);
          ix = KE::trunc(deposit_loc - iz * dims[0] * dims[1] - iy * dims[0]);

          if (ix > 0 && ix < dims[0] - 1 && iy > 0 && iy < dims[1] - 1 &&
              iz > 0 && iz < dims[2] - 1) {

            if (cell_idx == 0) {
              cell_loc = deposit_loc;
            } else if (cell_idx == 1) {
              cell_loc = deposit_loc - 1;
            } else if (cell_idx == 2) {
              cell_loc = deposit_loc - dims[0];
            } else if (cell_loc == 3) {
              cell_loc = deposit_loc - dims[0] - 1;
            } else if (cell_idx == 4) {
              cell_loc = deposit_loc - dims[0] * dims[1];
            } else if (cell_idx == 5) {
              cell_loc = deposit_loc - 1 - dims[0] * dims[1];
            } else if (cell_idx == 6) {
              cell_loc = deposit_loc - dims[0] - dims[0] * dims[1];
            } else if (cell_loc == 7) {
              cell_loc = deposit_loc - dims[0] - 1 - dims[0] * dims[1];
            }

            auto particles_start_it = KE::find(
                Kokkos::DefaultExecutionSpace(), KE::cbegin(permute_idx),
                KE::cend(permute_idx), deposit_loc);
            auto particles_end_it = KE::find(
                Kokkos::DefaultExecutionSpace(), KE::cbegin(permute_idx),
                KE::cend(permute_idx), deposit_loc + 1);
            size_t num_particles_in_this_cell =
                KE::distance(particles_start_it, particles_end_it);
            size_t first_particle_loc =
                KE::distance(KE::cbegin(permute_idx), particles_start_it);

	    /* TODO */

          } else {
            rho_dev(deposit_loc) = 0.0;
          }
        }
      });
}
