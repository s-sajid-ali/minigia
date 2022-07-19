#ifndef SPACE_CHARGE_3D_FD_KERNELS_H_
#define SPACE_CHARGE_3D_FD_KERNELS_H_

#include <Kokkos_Core.hpp>
#include <minigia/utils/multi_array_typedefs.hpp>

struct alg_force_extractor {
  karray1d_dev phi2;
  karray1d_dev enx;
  karray1d_dev eny;
  karray1d_dev enz;

  int gx, gy, gz;
  double ihx, ihy, ihz;
  double igxgy;
  double igx;

  alg_force_extractor(karray1d_dev const &phi2, karray1d_dev const &enx,
                      karray1d_dev const &eny, karray1d_dev const &enz,
                      std::array<int, 3> const &g,
                      std::array<double, 3> const &h)
      : phi2(phi2), enx(enx), eny(eny), enz(enz), gx(g[0]), gy(g[1]), gz(g[2]),
        ihx(0.5 / h[0]), ihy(0.5 / h[1]), ihz(0.5 / h[2]),
        igxgy(1.0 / (gx * gy)), igx(1.0 / gx) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    int iz = i * igxgy;
    int iy = (i - iz * gx * gy) * igx;
    int ix = i - iz * gx * gy - iy * gx;

    int ixl, ixr, iyl, iyr, izl, izr;

    double idx, idy, idz;

    // all boundaries will be skipped (set to 0)
    if (ix == 0 || ix == gx - 1 || iy == 0 || iy == gy - 1 || iz == 0 ||
        iz == gz - 1)
      return;

    ixl = ix - 1;
    ixr = ix + 1;
    iyl = iy - 1;
    iyr = iy + 1;
    izl = iz - 1;
    izr = iz + 1;

    idx = ihx;
    idy = ihy;
    idz = ihz;

    int idx_r = iz * gx * gy + iy * gx + ixr;
    int idx_l = iz * gx * gy + iy * gx + ixl;
    enx(i) = -(phi2(idx_r) - phi2(idx_l)) * idx;

    idx_r = iz * gx * gy + iyr * gx + ix;
    idx_l = iz * gx * gy + iyl * gx + ix;
    eny(i) = -(phi2(idx_r) - phi2(idx_l)) * idy;

    idx_r = izr * gx * gy + iy * gx + ix;
    idx_l = izl * gx * gy + iy * gx + ix;
    enz(i) = -(phi2(idx_r) - phi2(idx_l)) * idz;
  }
};

#endif
