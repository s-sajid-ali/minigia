#include <minigia/bunch/core_diagnostics.hpp>
#include <minigia/collective/deposit.hpp>
#include <minigia/utils/logger.hpp>
#include <minigia/utils/simple_timer.hpp>

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
  Kokkos::initialize(argc, argv);

  const double mass = 100.0;
  const double total_energy = 125.0;
  double real_particles = 10485760;
  double macroparticles = 2940000000000.0;

  Four_momentum four_momentum(pconstants::mp);
  Reference_particle refpart(1, four_momentum);

  Logger screen(0, LoggerV::INFO);
  Bunch bunch;
  bunch.read_file(std::string("particles.h5"));
  const std::array<int, 3> dims{32, 32, 128};

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

  karray1d_dev rho_dev("rho2_dev", dims[0] * dims[1] * dims[2]);

  {
    scoped_simple_timer timer("100its of old_method");
    for (int run = 0; run < 100; run++) {
      deposit_charge_rectangular_3d_kokkos_scatter_view(rho_dev, domain, dims,
                                                        bunch);
    }
  }
  simple_timer_print(screen);

  Kokkos::finalize();

  MPI_Finalize();
  return 0;
}
