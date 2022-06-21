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
  bunch_parts.access(0, 0) = 1.1;
  bunch_parts.access(0, 2) = 1.1;
  bunch_parts.access(0, 4) = 1.1;
  bunch.checkin_particles();

  Rectangular_grid_domain domain({4, 4, 4}, {4, 4, 4}, {1, 1, 1}, false);

  const std::array<int, 3> dims{4, 4, 4};
  karray1d_dev rho_dev("rho_dev", 4 * 4 * 4);

  deposit_charge_rectangular_3d_kokkos_scatter_view(rho_dev, domain, dims,
                                                    bunch);

  karray1d_hst rho_dev_hst = Kokkos::create_mirror_view(rho_dev);
  Kokkos::deep_copy(rho_dev_hst, rho_dev);
  Kokkos::fence();

  for (int idx = 0; idx < 64; idx++) {

    std::cout << "rho_dev_hst[" << idx << "] : " << rho_dev_hst(idx) << "\n";
  }

  MPI_Finalize();
  return 0;
}
