#include <iostream>

#include <minigia/collective/deposit.hpp>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  Commxx comm;
  Kokkos::ScopeGuard kokkos(argc, argv);

  const double mass = 100.0;
  const double total_energy = 125.0;

  const int total_num = 100;
  const double real_num = 2.0e12;
  auto const id = Bunch::id;

  Four_momentum fm(mass, total_energy);
  Reference_particle ref(pconstants::proton_charge, fm);

  Bunch bunch(ref, 1, 1, Commxx());
  bunch.checkout_particles();
  auto bunch_parts = bunch.get_host_particles();
  bunch_parts.access(0, 0) = 1.55;
  bunch_parts.access(0, 2) = 1.5;
  bunch_parts.access(0, 4) = 1.5;
  bunch.checkin_particles();

  Rectangular_grid_domain domain({6, 6, 6}, {6, 6, 6}, {2, 2, 2}, false);

  auto h = domain.get_cell_size();

  double weight0 = (bunch.get_real_num() / bunch.get_total_num()) *
                   bunch.get_particle_charge() * pconstants::e /
                   (h[0] * h[1] * h[2]);
  std::cout << "\n weight0 is " << weight0 << "\n\n";

  std::cout << "Domain "
            << "\n"
            << "grid shape "
            << " x : " << domain.get_grid_shape()[0]
            << " y : " << domain.get_grid_shape()[1]
            << " z : " << domain.get_grid_shape()[2] << "\n"
            << "left       "
            << " x : " << domain.get_left()[0]
            << " y : " << domain.get_left()[1]
            << " z : " << domain.get_left()[2] << "\n"
            << "cell size  "
            << " x : " << domain.get_cell_size()[0]
            << " y : " << domain.get_cell_size()[1]
            << " z : " << domain.get_cell_size()[2] << "\n";

  const std::array<int, 3> dims{6, 6, 6};
  karray1d_dev rho_dev("rho_dev", 6 * 6 * 6);

  deposit_charge_rectangular_3d_kokkos_scatter_view(rho_dev, domain, dims,
                                                    bunch);

  karray1d_hst rho_dev_hst = Kokkos::create_mirror_view(rho_dev);
  Kokkos::deep_copy(rho_dev_hst, rho_dev);
  Kokkos::fence();

  auto sums = 0.0;

  for (int idx = 0; idx < 6 * 6 * 6; idx++) {

    std::cout << "rho_dev_hst[" << idx << "] : " << rho_dev_hst(idx) / weight0
              << "\n";

    sums += rho_dev_hst(idx) / weight0;
  }

  std::cout << "\nsums of all deposit fractions :" << sums << "\n";

  MPI_Finalize();
  return 0;
}
