#include <Kokkos_ScatterView.hpp>
#include <chrono>
#include <iostream>
#include <minigia/bunch/bunch.hpp>
#include <minigia/bunch/core_diagnostics.hpp>
#include <minigia/foundation/physical_constants.hpp>
#include <minigia/utils/logger.hpp>
#include <minigia/utils/simple_timer.hpp>
#include <thread>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  Commxx comm;

  {
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
      scoped_simple_timer timer("old_method");
      for (int run = 0; run < 100; run++) {
        auto mean = Core_diagnostics::calculate_mean(bunch);
        auto stddev = Core_diagnostics::calculate_std(bunch, mean);
      }
    }

    const auto particles = bunch.get_local_particles();
    const auto masks = bunch.get_local_particle_masks();

    const auto total_bunch_particles = bunch.get_total_num();
    const auto local_bunch_capacity = bunch.size();

    Kokkos::View<double[6], Kokkos::DefaultHostExecutionSpace> mean_and_stddev(
        "mean_stddev");
    // std::array<double, 6> mean_and_stddev{0, 0, 0, 0, 0, 0};
    {
      scoped_simple_timer timer("new_method");

      auto instances = Kokkos::Experimental::partition_space(
          Kokkos::DefaultExecutionSpace(), 1, 1, 1);
      for (int run = 0; run < 100; run++) {

        for (int instance_id = 0; instance_id < 3; instance_id++) {
          // std::cout << "starting 1st set of instances with id : " <<
          // instance_id
          //                   << " at " << MPI_Wtime() << "\n";
          int idx1 = 0 + instance_id;
          int idx2 = 3 + instance_id;
          int idx3 = 2 * instance_id;
          Kokkos::parallel_reduce(
              "position_sum",
              Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(
                  instances[instance_id], 0, local_bunch_capacity),
              KOKKOS_LAMBDA(const int &i, double &sum_pos) {
                if (masks(i)) {
                  sum_pos += particles(i, idx3);
                }
              },
              mean_and_stddev[idx1]);
          Kokkos::parallel_reduce(
              "position_sumsquares",
              Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(
                  instances[instance_id], 0, local_bunch_capacity),
              KOKKOS_LAMBDA(const int &i, double &sum_squarepos) {
                if (masks(i)) {
                  sum_squarepos += particles(i, idx3) * particles(i, idx3);
                }
              },
              mean_and_stddev[idx2]);
        }
        for (int instance_id = 0; instance_id < 3; instance_id++) {
          // std::cout << "starting 2nd set of instances with id : " <<
          // instance_id
          //           << " at " << MPI_Wtime() << "\n";
          int idx1 = 0 + instance_id;
          int idx2 = 3 + instance_id;
          instances[instance_id].fence();
        }

        int status;
        status = MPI_Allreduce(MPI_IN_PLACE, mean_and_stddev.data(), 6,
                               MPI_DOUBLE, MPI_SUM, bunch.get_comm());
        if (status != MPI_SUCCESS) {
          std::cout << status << "\n";
        }

        /* mean_and_stddev has sum_and_sumsquares on each MPI rank */
        for (int j = 0; j < 3; j++) {
          mean_and_stddev[j] = mean_and_stddev[j] / total_bunch_particles;
          mean_and_stddev[j + 3] =
              std::sqrt((mean_and_stddev[j + 3]) / total_bunch_particles -
                        std::pow(mean_and_stddev[j], 2));
        }
      }
    }
    simple_timer_print(screen);
  }

  Kokkos::finalize();
  MPI_Finalize();

  return 0;
}
