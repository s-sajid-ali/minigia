#include <minigia/collective/deposit.hpp>
#include <minigia/utils/logger.hpp>
#include <minigia/utils/simple_timer.hpp>

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

  const std::array<int, 3> dims{66, 64, 256};
  karray1d_dev rho_dev("rho2_dev", 1081344);

  {
    scoped_simple_timer timer("100its of old_method");
    for (int run = 0; run < 100; run++) {
      deposit_charge_rectangular_3d_kokkos_scatter_view(rho_dev, domain, dims,
                                                        bunch);
    }
  }
  simple_timer_print(screen);

  MPI_Finalize();
  return 0;
}
