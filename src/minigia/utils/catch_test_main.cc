#include <catch2/catch_session.hpp>

#include <Kokkos_Core.hpp>
#include <mpi.h>
#include <petscsys.h>

static std::string help = "This is a minigia test program!\n\n";

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc, &argv, (char *)0, help.c_str());

  int result = Catch::Session().run(argc, argv);

  ierr = PetscFinalize();
  Kokkos::finalize();
  MPI_Finalize();
  return result;
}
