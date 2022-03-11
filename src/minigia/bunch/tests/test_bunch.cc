#include <catch2/catch_test_macros.hpp>

#include <minigia/foundation/physical_constants.hpp>

#include "../bunch.hpp"

const double mass = 100.0;
const double total_energy = 125.0;
const int total_num = 100;
const double real_num = 2.0e12;
auto const id = Bunch::id;

TEST_CASE("Bunch", "[Bunch]") {
  Four_momentum fm(mass, total_energy);
  Reference_particle ref(pconstants::proton_charge, fm);
  Bunch bunch(ref, 1024, 1e13, Commxx());

  CHECK(bunch.get_local_num() == 1024);

  bunch.checkout_particles();
  auto parts = bunch.get_host_particles();

  CHECK(parts(3, id) == 3);
  CHECK(parts(307, id) == 307);

  auto idx1 = bunch.search_particle(12);
  CHECK(idx1 == 12);

  auto idx2 = bunch.search_particle(13, 12);
  CHECK(idx2 == 13);

  auto idx3 = bunch.search_particle(15, 15);
  CHECK(idx3 == 15);

  auto p2 = bunch.get_particles_in_range(123, 6);
  CHECK(p2(0, 6) == 123);
  CHECK(p2(1, 6) == 124);
  CHECK(p2(4, 6) == 127);
}
