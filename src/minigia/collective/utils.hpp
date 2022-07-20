#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

KOKKOS_INLINE_FUNCTION
void get_leftmost_indices_offset(double pos, double left, double inv_cell_size,
                                 int &idx, double &off) {
  double scaled_location = (pos - left) * inv_cell_size - 0.5;
  idx = Kokkos::Experimental::nearbyint(scaled_location);
  off = scaled_location - idx;

#if VERBOSE == 1
  std::cout << "particle with pos : " << pos << ", left : " << left
            << ", inv_cell_size : " << inv_cell_size << ", idx : " << idx
            << ", off : " << off << std::endl;
#endif
}
