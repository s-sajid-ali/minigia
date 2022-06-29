#ifndef DEPOSIT_NEW_H_
#define DEPOSIT_NEW_H_

#include <minigia/bunch/bunch.hpp>

#include "rectangular_grid_domain.hpp"

#include <Kokkos_Sort.hpp>

/* Fill idx view with locations of particles */
void bunch_get_idx(Bunch &bunch, Rectangular_grid_domain &domain,
                   karray1i_dev &idx);

/* deposit particles on grid using a sorted bunch
   v1: the charge deposited at each grid point is
   calcucated independently */
void deposit_charge_rectangular_3d_sorted_v1(
    karray1d_dev &rho_dev,
    Kokkos::BinSort<karray1i_dev, Kokkos::BinOp1D<karray1i_dev>> &sorter,
    Rectangular_grid_domain &domain, std::array<int, 3> const &dims,
    Bunch const &bunch);

#endif /* DEPOSIT_NEW_H_ */
