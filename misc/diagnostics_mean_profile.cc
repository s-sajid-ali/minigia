#include <minigia/bunch/bunch.hpp>
#include <minigia/bunch/core_diagnostics.hpp>
#include <minigia/foundation/physical_constants.hpp>
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

  Kokkos::View<double[6], Kokkos::DefaultHostExecutionSpace> mean("cal_mean");

  {
    scoped_simple_timer timer("100its of new_method");

    auto instances = Kokkos::Experimental::partition_space(
        Kokkos::DefaultExecutionSpace(), 1, 1, 1);

    /* spatial mean and stddev! */
    for (int run = 0; run < 100; run++) {

      for (int instance_id = 0; instance_id < 3; instance_id++) {
        int dst_idx1 = 0 + instance_id;
        int dst_idx2 = 3 + instance_id;
        int src_idx1 = 0 + instance_id;
        int src_idx2 = 3 + instance_id;
        Kokkos::parallel_reduce(
            "position_sum",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(
                instances[instance_id], 0, local_bunch_capacity),
            KOKKOS_LAMBDA(const int &i, double &sum_pos) {
              if (masks(i)) {
                sum_pos += particles(i, src_idx1);
              }
            },
            mean[dst_idx1]);
        Kokkos::parallel_reduce(
            "momenta_sum",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(
                instances[instance_id], 0, local_bunch_capacity),
            KOKKOS_LAMBDA(const int &i, double &sum_momenta) {
              if (masks(i)) {
                sum_momenta += particles(i, src_idx2);
              }
            },
            mean[dst_idx2]);
      }
      for (int instance_id = 0; instance_id < 3; instance_id++) {
        instances[instance_id].fence();
      }

      int status;
      status = MPI_Allreduce(MPI_IN_PLACE, mean.data(), 6, MPI_DOUBLE, MPI_SUM,
                             bunch.get_comm());
      if (status != MPI_SUCCESS) {
        std::cout << status << "\n";
      }

      for (int j = 0; j < 3; j++) {
        mean[j] = mean[j] / total_bunch_particles;
      }
    }
  }
  simple_timer_print(screen);

  MPI_Finalize();

  return 0;
}
