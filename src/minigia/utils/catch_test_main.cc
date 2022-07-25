#include <catch2/catch_session.hpp>

#include <Kokkos_Core.hpp>
#include <mpi.h>
#include <petscsys.h>

static std::string help = "This is a minigia test program!\n\n";

int main(int argc, char *argv[]) {
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc, &argv, (char *)0, help.c_str());

  auto initargs = Kokkos::InitArguments{}; /* use default constructor */
  Kokkos::initialize(initargs);

  int result = Catch::Session().run();

  Kokkos::finalize();
  ierr = PetscFinalize();
  return result;
}
