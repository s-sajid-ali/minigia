#include <minigia/bunch/core_diagnostics.hpp>
#include <minigia/foundation/physical_constants.hpp>

#include "deposit.hpp"
#include "utils.hpp"

#include <Kokkos_ScatterView.hpp>

#define VERBOSE 0

namespace deposit_impl {
using scatter_t =
    Kokkos::Experimental::ScatterView<double *, Kokkos::LayoutLeft>;

KOKKOS_INLINE_FUNCTION
bool ingrid(int i, int g) { return i >= 0 && i < g; }

KOKKOS_INLINE_FUNCTION
bool ingrid(int ix, int iy, int iz, int gx, int gy, int gz) {
  // return ix>=0 && ix<gx && iy>=0 && iy<gy && iz>=0 && iz<gz;

  // exclude edges
  return ix > 0 && ix < gx - 1 && iy > 0 && iy < gy - 1 && iz > 0 &&
         iz < gz - 1;
}

struct rho_zeroer {
  karray1d_dev rho;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const { rho(i) = 0.0; }
};

// use scatter view
struct sv_zyx_rho_reducer_non_periodic {
  ConstParticles p;
  ConstParticleMasks masks;

  scatter_t scatter;

  int gx, gy, gz; // original grid size
  int dx, dy, dz; // dimensions of the grid
  double ihx, ihy, ihz;
  double lx, ly, lz;
  double w0;

  sv_zyx_rho_reducer_non_periodic(ConstParticles const &p,
                                  ConstParticleMasks const &masks,
                                  scatter_t const &scatter,
                                  std::array<int, 3> const &g,
                                  std::array<int, 3> const &d,
                                  std::array<double, 3> const &h,
                                  std::array<double, 3> const &l, double w0)
      : p(p), masks(masks), scatter(scatter), gx(g[0]), gy(g[1]), gz(g[2]),
        dx(d[0]), dy(d[1]), dz(d[2]), ihx(1.0 / h[0]), ihy(1.0 / h[1]),
        ihz(1.0 / h[2]), lx(l[0]), ly(l[1]), lz(l[2]), w0(w0) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    if (masks(i)) {
      auto access = scatter.access();

      int ix, iy, iz;
      double ox, oy, oz;

      get_leftmost_indices_offset(p(i, 0), lx, ihx, ix, ox);
      get_leftmost_indices_offset(p(i, 2), ly, ihy, iy, oy);
      get_leftmost_indices_offset(p(i, 4), lz, ihz, iz, oz);

      double aox = 1.0 - ox;
      double aoy = 1.0 - oy;
      double aoz = 1.0 - oz;

      int base = iz * dx * dy;

#if VERBOSE == 1
      std::cout << "particle is at ix : " << ix << ", iy : " << iy
                << ", iz : " << iz << ", with weights aox : " << aox
                << ", aoy : " << aoy << ", aoz : " << aoz << std::endl;
#endif

      if (ingrid(ix, iy, iz, gx, gy, gz)) {
#if VERBOSE == 1
        std::cout << "particle is in at ix, iy, iz! with weights aox : " << aox
                  << ", aoy : " << aoy << ", aoz : " << aoz << std::endl;
#endif
        access(base + iy * dx + ix) += w0 * aox * aoy * aoz;
      }

      if (ingrid(ix + 1, iy, iz, gx, gy, gz)) {
#if VERBOSE == 1
        std::cout << "particle is in at ix+1, iy, iz! with weights ox : " << ox
                  << ", aoy : " << aoy << ", aoz : " << aoz << std::endl;
#endif
        access(base + iy * dx + ix + 1) += w0 * ox * aoy * aoz;
      }

      if (ingrid(ix, iy + 1, iz, gx, gy, gz)) {
#if VERBOSE == 1
        std::cout << "particle is in at ix, iy+1, iz! with weights aox : "
                  << aox << ", oy : " << oy << ", aoz : " << aoz << std::endl;
#endif
        access(base + (iy + 1) * dx + ix) += w0 * aox * oy * aoz;
      }

      if (ingrid(ix + 1, iy + 1, iz, gx, gy, gz)) {
#if VERBOSE == 1
        std::cout << "particle is in at ix+1, iy+1, iz! with weights ox : "
                  << ox << ", oy : " << oy << ", aoz : " << aoz << std::endl;
#endif
        access(base + (iy + 1) * dx + ix + 1) += w0 * ox * oy * aoz;
      }

      base = (iz + 1) * dx * dy;

      if (ingrid(ix, iy, iz + 1, gx, gy, gz)) {
#if VERBOSE == 1
        std::cout << "particle is in at ix, iy, iz+1! with weights aox : "
                  << aox << ", aoy : " << aoy << ", oz : " << oz << std::endl;
#endif
        access(base + iy * dx + ix) += w0 * aox * aoy * oz;
      }

      if (ingrid(ix + 1, iy, iz + 1, gx, gy, gz)) {
#if VERBOSE == 1
        std::cout << "particle is in at ix+1, iy, iz! with weights ox : " << ox
                  << ", aoy : " << aoy << ", oz : " << oz << std::endl;
#endif
        access(base + iy * dx + ix + 1) += w0 * ox * aoy * oz;
      }

      if (ingrid(ix, iy + 1, iz + 1, gx, gy, gz)) {
#if VERBOSE == 1
        std::cout << "particle is in at ix, iy+1, iz+1! with weights aox : "
                  << aox << ", oy : " << oy << ", oz : " << oz << std::endl;
#endif
        access(base + (iy + 1) * dx + ix) += w0 * aox * oy * oz;
      }

      if (ingrid(ix + 1, iy + 1, iz + 1, gx, gy, gz)) {
#if VERBOSE == 1
        std::cout << "particle is in at ix+1, iy+1, iz+1! with weights ox : "
                  << ox << ", oy : " << oy << ", oz : " << oz << std::endl;
#endif
        access(base + (iy + 1) * dx + ix + 1) += w0 * ox * oy * oz;
      }
    }
  }
};

} // namespace deposit_impl

void deposit_charge_rectangular_3d_kokkos_scatter_view(
    karray1d_dev &rho_dev, Rectangular_grid_domain &domain,
    std::array<int, 3> const &dims, Bunch const &bunch) {
  using namespace deposit_impl;

  auto g = domain.get_grid_shape();
  auto h = domain.get_cell_size();
  auto l = domain.get_left();

  auto parts = bunch.get_local_particles();
  auto masks = bunch.get_local_particle_masks();
  int nparts = bunch.size();

  double weight0 = (bunch.get_real_num() / bunch.get_total_num()) *
                   bunch.get_particle_charge() * pconstants::e /
                   (h[0] * h[1] * h[2]);

  if (rho_dev.extent(0) < g[0] * g[1] * g[2])
    throw std::runtime_error("insufficient size for rho in deposit charge");

  // zero first
  rho_zeroer rz{rho_dev};
  Kokkos::parallel_for(rho_dev.extent(0), rz);
  Kokkos::fence();

  // deposit
  scatter_t scatter(rho_dev);
  sv_zyx_rho_reducer_non_periodic rr(parts, masks, scatter, g, dims, h, l,
                                     weight0);

  Kokkos::parallel_for(nparts, rr);
  Kokkos::Experimental::contribute(rho_dev, scatter);

  Kokkos::fence();
}
