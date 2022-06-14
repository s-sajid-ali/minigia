#include <minigia/bunch/bunch.hpp>
#include <minigia/bunch/core_diagnostics.hpp>
#include <minigia/foundation/physical_constants.hpp>
#include <minigia/utils/logger.hpp>
#include <minigia/utils/simple_timer.hpp>

#include <array>
#include <execution>

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

  {
    scoped_simple_timer timer("100its of old_method");
    for (int run = 0; run < 100; run++) {
      auto mean = Core_diagnostics::calculate_mean(bunch);
    }
  }

  const auto particles = bunch.get_local_particles();
  const auto masks = bunch.get_local_particle_masks();

  const auto total_bunch_particles = bunch.get_total_num();
  const auto local_bunch_capacity = bunch.size();

  auto instances = Kokkos::Experimental::partition_space(
      Kokkos::DefaultExecutionSpace(), 1, 1, 1, 1, 1, 1);
  Kokkos::View<double[6], Kokkos::DefaultHostExecutionSpace::memory_space> mean(
      "cal_mean");
  {
    scoped_simple_timer timer("100its of new_method");

    /* spatial mean and stddev! */
    for (int run = 0; run < 100; run++) {

      std::array<int, 6> arr = {0, 1, 2, 3, 4, 5};

      for (int instance_id = 0; instance_id < 6; instance_id++) {
        int dst_idx = instance_id;
        int src_idx = instance_id;
        Kokkos::parallel_reduce(
            "cal_mean",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(
                instances[instance_id], 0, local_bunch_capacity),
            KOKKOS_LAMBDA(const int &i, double &sum_pos) {
              if (masks(i)) {
                sum_pos += particles(i, src_idx);
              }
            },
            mean[dst_idx]);
      }
      for (int instance_id = 0; instance_id < 6; instance_id++) {
        instances[instance_id].fence();
      }

      int status;
      status = MPI_Allreduce(MPI_IN_PLACE, mean.data(), 6, MPI_DOUBLE, MPI_SUM,
                             bunch.get_comm());
      if (status != MPI_SUCCESS) {
        std::cout << status << "\n";
      }

      for (int j = 0; j < 6; j++) {
        mean[j] = mean[j] / total_bunch_particles;
      }
    }
  }
  simple_timer_print(screen);

  MPI_Finalize();

  return 0;
}
