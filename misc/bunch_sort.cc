#include <minigia/bunch/bunch.hpp>
#include <minigia/bunch/core_diagnostics.hpp>
#include <minigia/collective/deposit_new.hpp>
#include <minigia/foundation/physical_constants.hpp>
#include <minigia/utils/hdf5_file.hpp>
#include <minigia/utils/logger.hpp>
#include <minigia/utils/simple_timer.hpp>

#include <Kokkos_Sort.hpp>

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

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  Commxx comm;
  Kokkos::ScopeGuard kokkos(argc, argv);

  Logger screen(0, LoggerV::INFO);
  Bunch bunch;
  bunch.read_file(std::string("particles.h5"));

  const std::array<int, 3> dims{32, 32, 128};
  karray1d_dev rho_dev("rho_dev",
                       (dims[0] + 1) * (dims[1] + 1) * (dims[2] + 1));

  auto mean = Core_diagnostics::calculate_mean(bunch);
  auto std = Core_diagnostics::calculate_std(bunch, mean);

  const double tiny = 1.0e-10;

  const auto ix = Bunch::x;
  const auto iy = Bunch::y;
  const auto iz = Bunch::z;
  if ((std[ix] < tiny) && (std[iy] < tiny) && (std[iz] < tiny)) {
    throw std::runtime_error(
        "Space_charge_3d_open_hockney_eigen::update_domain: "
        "all three spatial dimensions have neglible extent");
  }

  std::array<double, 3> offset{mean[ix], mean[iy], mean[iz]};

  int n_sigma = 8;

  std::array<double, 3> size{
      n_sigma * get_smallest_non_tiny(std[0], std[2], std[4], tiny),
      n_sigma * get_smallest_non_tiny(std[2], std[0], std[4], tiny),
      n_sigma * get_smallest_non_tiny(std[4], std[0], std[2], tiny)};

  Rectangular_grid_domain domain =
      Rectangular_grid_domain(dims, size, offset, false);

  karray1i_dev idx("deposit_idx", bunch.capacity());

  bunch_get_idx(bunch, domain, idx);
  Kokkos::fence();

  std::cout << "bunch capacity : " << bunch.capacity() << "\n";
  std::cout << "bunch size : " << bunch.size() << "\n";

  {
    Hdf5_file file("idx.h5", Hdf5_file::Flag::truncate, Commxx());
    file.write("idx", idx.data(), idx.size(), true);
  }

  using KeyViewType = karray1i_dev;
  using BinOp = Kokkos::BinOp1D<KeyViewType>;

  int begin = 0;
  int end = bunch.capacity();
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

  // deposit_charge_rectangular_3d_sorted_v1(rho_dev, Sorter, domain, dims,
  // bunch);
  //  simple_timer_print(screen);

  MPI_Finalize();
  return 0;
}
