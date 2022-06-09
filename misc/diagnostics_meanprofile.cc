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
      auto stddev = Core_diagnostics::calculate_std(bunch, mean);
    }
  }

  const auto particles = bunch.get_local_particles();
  const auto masks = bunch.get_local_particle_masks();

  const auto total_bunch_particles = bunch.get_total_num();
  const auto local_bunch_capacity = bunch.size();

  Kokkos::View<double[6], Kokkos::DefaultHostExecutionSpace>
      spatial_mean_and_stddev("spatial_mean_stddev");
  Kokkos::View<double[6], Kokkos::DefaultHostExecutionSpace>
      momenta_mean_and_stddev("momenta_mean_stddev");

  {
    scoped_simple_timer timer("100its of new_method");

    auto instances = Kokkos::Experimental::partition_space(
        Kokkos::DefaultExecutionSpace(), 1, 1, 1);

    /* spatial mean and stddev! */
    for (int run = 0; run < 100; run++) {

      for (int instance_id = 0; instance_id < 3; instance_id++) {
        int dst_idx1 = 0 + instance_id;
        int dst_idx2 = 3 + instance_id;
        int src_idx = 2 * instance_id;
        Kokkos::parallel_reduce(
            "position_sum",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(
                instances[instance_id], 0, local_bunch_capacity),
            KOKKOS_LAMBDA(const int &i, double &sum_pos) {
              if (masks(i)) {
                sum_pos += particles(i, src_idx);
              }
            },
            spatial_mean_and_stddev[dst_idx1]);
        Kokkos::parallel_reduce(
            "position_sumsquares",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(
                instances[instance_id], 0, local_bunch_capacity),
            KOKKOS_LAMBDA(const int &i, double &sum_squarepos) {
              if (masks(i)) {
                sum_squarepos += particles(i, src_idx) * particles(i, src_idx);
              }
            },
            spatial_mean_and_stddev[dst_idx2]);
      }
      for (int instance_id = 0; instance_id < 3; instance_id++) {
        instances[instance_id].fence();
      }

      int status;
      status = MPI_Allreduce(MPI_IN_PLACE, spatial_mean_and_stddev.data(), 6,
                             MPI_DOUBLE, MPI_SUM, bunch.get_comm());
      if (status != MPI_SUCCESS) {
        std::cout << status << "\n";
      }

      /* spatial_mean_and_stddev has sum_and_sumsquares on each MPI rank */
      for (int j = 0; j < 3; j++) {
        spatial_mean_and_stddev[j] =
            spatial_mean_and_stddev[j] / total_bunch_particles;
        spatial_mean_and_stddev[j + 3] =
            std::sqrt((spatial_mean_and_stddev[j + 3]) / total_bunch_particles -
                      std::pow(spatial_mean_and_stddev[j], 2));
      }
    }

    /* momenta mean and stddev! */
    for (int run = 0; run < 100; run++) {

      for (int instance_id = 0; instance_id < 3; instance_id++) {
        int dst_idx1 = 0 + instance_id;
        int dst_idx2 = 3 + instance_id;
        int src_idx = 2 * instance_id;
        Kokkos::parallel_reduce(
            "momenta_sum",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(
                instances[instance_id], 0, local_bunch_capacity),
            KOKKOS_LAMBDA(const int &i, double &sum_pos) {
              if (masks(i)) {
                sum_pos += particles(i, src_idx);
              }
            },
            momenta_mean_and_stddev[dst_idx1]);
        Kokkos::parallel_reduce(
            "momenta_sumsquares",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(
                instances[instance_id], 0, local_bunch_capacity),
            KOKKOS_LAMBDA(const int &i, double &sum_squarepos) {
              if (masks(i)) {
                sum_squarepos += particles(i, src_idx) * particles(i, src_idx);
              }
            },
            momenta_mean_and_stddev[dst_idx2]);
      }
      for (int instance_id = 0; instance_id < 3; instance_id++) {
        instances[instance_id].fence();
      }

      int status;
      status = MPI_Allreduce(MPI_IN_PLACE, momenta_mean_and_stddev.data(), 6,
                             MPI_DOUBLE, MPI_SUM, bunch.get_comm());
      if (status != MPI_SUCCESS) {
        std::cout << status << "\n";
      }

      /* momenta_mean_and_stddev has sum_and_sumsquares on each MPI rank */
      for (int j = 0; j < 3; j++) {
        momenta_mean_and_stddev[j] =
            momenta_mean_and_stddev[j] / total_bunch_particles;
        spatial_mean_and_stddev[j + 3] =
            std::sqrt((momenta_mean_and_stddev[j + 3]) / total_bunch_particles -
                      std::pow(momenta_mean_and_stddev[j], 2));
      }
    }
  }
  simple_timer_print(screen);

  MPI_Finalize();

  return 0;
}
