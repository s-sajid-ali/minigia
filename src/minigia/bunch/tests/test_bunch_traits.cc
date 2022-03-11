#include <catch2/catch_test_macros.hpp>

#include <map>
#include <memory>
#include <type_traits>
#include <vector>

#include <Kokkos_Core.hpp>

struct TinyBunch {
  explicit TinyBunch(int s) : size(s) {}

  TinyBunch(TinyBunch const &) = delete;
  TinyBunch(TinyBunch &&) = default;

  int size;
  Kokkos::View<double *> particles;
  std::map<std::string, std::unique_ptr<int>> diags;
  // std::unique_ptr<int> pi;
};

static_assert(std::is_move_constructible<TinyBunch>::value,
    "TinyBunch isn't move constructible");

TEST_CASE("Bunch traits", "[Bunch]") {
  CHECK(true);

  std::vector<TinyBunch> bunches;
  bunches.emplace_back(3);
}
