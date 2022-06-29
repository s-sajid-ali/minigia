#include <minigia/bunch/bunch.hpp>
#include <minigia/collective/deposit_new.hpp>
#include <minigia/foundation/physical_constants.hpp>
#include <minigia/utils/hdf5_file.hpp>
#include <minigia/utils/logger.hpp>
#include <minigia/utils/simple_timer.hpp>

#include <Kokkos_Sort.hpp>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  Commxx comm;
  Kokkos::ScopeGuard kokkos(argc, argv);

  const double mass = 100.0;
  const double total_energy = 125.0;
  double real_particles = 10485760;
  double macroparticles = 2940000000000.0;

  Four_momentum four_momentum(pconstants::mp);
  Reference_particle refpart(1, four_momentum);

  Logger screen(0, LoggerV::INFO);
  Bunch bunch;
  bunch.read_file(std::string("particles.h5"));

  Rectangular_grid_domain domain(
      {32, 32, 128},
      {0.045792656585374118, 0.021010777363394857, 0.10253620986616825},
      {1.0473336722046451e-13, 5.7671527053847625e-14, 1.2976955845519457e-7},
      true);

  const std::array<int, 3> dims{32, 32, 128};
  karray1d_dev rho_dev("rho_dev",
                       (dims[0] + 1) * (dims[1] + 1) * (dims[2] + 1));

  karray1i_dev idx("deposit_idx", bunch.size());

  bunch_get_idx(bunch, domain, idx);
  Kokkos::fence();

  {
    Hdf5_file file("idx.h5", Hdf5_file::Flag::truncate, Commxx());
    file.write("idx", idx.data(), idx.size(), true);
  }

  using KeyViewType = karray1i_dev;
  using BinOp = Kokkos::BinOp1D<KeyViewType>;

  int begin = 0;
  int end = idx.size() - 1;
  auto min = 0;
  auto max = std::numeric_limits<int>::max();

  std::cout << "begin : " << begin << " ; "
            << "end   : " << end << "\n";
  std::cout << "min   : " << min << " ; "
            << "max   : " << max << "\n";

  BinOp binner(end - begin, min, max);
  Kokkos::BinSort<KeyViewType, BinOp> Sorter(idx, binner, false);
  Sorter.create_permute_vector(Kokkos::DefaultExecutionSpace());

  auto parts = bunch.get_local_particles();

  {
    auto parts_pos_x = Kokkos::subview(parts, Kokkos::ALL, 0);
    auto parts_pos_y = Kokkos::subview(parts, Kokkos::ALL, 2);
    auto parts_pos_z = Kokkos::subview(parts, Kokkos::ALL, 4);

    Sorter.sort(Kokkos::DefaultExecutionSpace(), parts_pos_x);
    Sorter.sort(Kokkos::DefaultExecutionSpace(), parts_pos_y);
    Sorter.sort(Kokkos::DefaultExecutionSpace(), parts_pos_z);
    Kokkos::fence();
  }

  deposit_charge_rectangular_3d_sorted_v1(rho_dev, Sorter, domain, dims, bunch);
  //  simple_timer_print(screen);

  MPI_Finalize();
  return 0;
}
