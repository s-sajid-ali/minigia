#include <catch2/catch_session.hpp>

#include <Kokkos_Core.hpp>
#include <mpi.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);

  int result = Catch::Session().run(argc, argv);

  Kokkos::finalize();
  MPI_Finalize();
  return result;
}
