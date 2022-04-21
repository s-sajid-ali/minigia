#ifndef DEPOSIT_H_
#define DEPOSIT_H_

#include <minigia/bunch/bunch.hpp>

#include "rectangular_grid_domain.hpp"

void deposit_charge_rectangular_3d_kokkos_scatter_view(
    karray1d_dev &rho_dev, Rectangular_grid_domain &domain,
    std::array<int, 3> const &dims, Bunch const &bunch);

#endif /* DEPOSIT_H_ */
