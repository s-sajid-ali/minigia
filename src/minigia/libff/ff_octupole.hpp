#ifndef FF_OCTUPOLE_H
#define FF_OCTUPOLE_H

#include <minigia/utils/simple_timer.hpp>

#include "ff_algorithm.hpp"
#include "ff_patterned_propagator.hpp"

namespace FF_octupole {
template <class T>
KOKKOS_INLINE_FUNCTION void kick(T const &x, T &xp, T const &y, T &yp,
                                 T const &, double const *kL) {
  FF_algorithm::thin_octupole_unit(x, xp, y, yp, kL);
}

template <class BunchT>
void apply(Lattice_element_slice const &slice, BunchT &bunch) {
  auto const &element = slice.get_lattice_element();
  double length = slice.get_right() - slice.get_left();

  double k[2] = {element.get_double_attribute("k3", 0.0),
                 element.get_double_attribute("k3s", 0.0)};

  // tilting
  double tilt = element.get_double_attribute("tilt", 0.0);
  if (tilt != 0.0) {
    std::complex<double> ck2(k[0], k[1]);
    ck2 = ck2 * exp(std::complex<double>(0.0, -4.0 * tilt));
    k[0] = ck2.real();
    k[1] = ck2.imag();
  }

  // scaling
  Reference_particle &ref_l = bunch.get_design_reference_particle();
  Reference_particle const &ref_b = bunch.get_reference_particle();

  double brho_l = ref_l.get_momentum() / ref_l.get_charge(); // GV/c
  double brho_b = ref_b.get_momentum() *
                  (1.0 + ref_b.get_state()[Bunch::dpop]) /
                  ref_l.get_charge(); // GV/c

  double scale = brho_l / brho_b;

  k[0] *= scale;
  k[1] *= scale;

  using gsv_t = typename BunchT::gsv_t;
  using pp = FF_patterned_propagator<BunchT, gsv_t, kick<gsv_t>, kick<double>>;

  if (close_to_zero(length)) {
    pp::get_reference_cdt_zero(ref_l, k);

    pp::apply_thin_kick(bunch, ParticleGroup::regular, k);
    pp::apply_thin_kick(bunch, ParticleGroup::spectator, k);
  } else {
    double pref = bunch.get_reference_particle().get_momentum();
    double mass = bunch.get_mass();

    // yoshida steps
    int steps = (int)slice.get_lattice_element().get_double_attribute(
        "yoshida_steps", 4.0);

    double ref_cdt = pp::get_reference_cdt_yoshida(ref_l, length, k, steps);

    // propagate
    pp::apply_yoshida_kick(bunch, ParticleGroup::regular, pref, mass, ref_cdt,
                           length, k, steps);

    pp::apply_yoshida_kick(bunch, ParticleGroup::spectator, pref, mass, ref_cdt,
                           length, k, steps);

    // trajectory
    bunch.get_reference_particle().increment_trajectory(length);
  }
}
} // namespace FF_octupole

#endif // FF_SEXTUPOLE_H
