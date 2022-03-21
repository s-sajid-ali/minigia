#ifndef FF_ELEMENT_H
#define FF_ELEMENT_H

#include <minigia/bunch/bunch.hpp>
#include <minigia/bunch/bunch_particles.hpp>
#include <minigia/foundation/reference_particle.hpp>
#include <minigia/lattice/lattice_element.hpp>
#include <minigia/lattice/lattice_element_slice.hpp>

#include "ff_drift.hpp"
#include "ff_elens.hpp"
#include "ff_foil.hpp"
#include "ff_kicker.hpp"
#include "ff_multipole.hpp"
#include "ff_nllens.hpp"
#include "ff_octupole.hpp"
#include "ff_quadrupole.hpp"
#include "ff_rfcavity.hpp"
#include "ff_sbend.hpp"
#include "ff_sextupole.hpp"
#include "ff_solenoid.hpp"

namespace FF_element {
template <class BUNCH>
void apply(Lattice_element_slice const &slice, BUNCH &b) {
  auto const &elm = slice.get_lattice_element();
  auto t = elm.get_type();

  switch (t) {
  case element_type::drift:
    FF_drift::apply(slice, b);
    break;
  case element_type::sbend:
    FF_sbend::apply(slice, b);
    break;
  case element_type::quadrupole:
    FF_quadrupole::apply(slice, b);
    break;

  case element_type::multipole:
    FF_multipole::apply(slice, b);
    break;
  case element_type::sextupole:
    FF_sextupole::apply(slice, b);
    break;
  case element_type::octupole:
    FF_octupole::apply(slice, b);
    break;

  case element_type::hkicker:
    FF_kicker::apply(slice, b);
    break;
  case element_type::vkicker:
    FF_kicker::apply(slice, b);
    break;
  case element_type::kicker:
    FF_kicker::apply(slice, b);
    break;

  case element_type::solenoid:
    FF_solenoid::apply(slice, b);
    break;
  case element_type::rfcavity:
    FF_rfcavity::apply(slice, b);
    break;
  case element_type::elens:
    FF_elens::apply(slice, b);
    break;
  case element_type::nllens:
    FF_nllens::apply(slice, b);
    break;

  case element_type::monitor:
    FF_drift::apply(slice, b);
    break;
  case element_type::hmonitor:
    FF_drift::apply(slice, b);
    break;
  case element_type::vmonitor:
    FF_drift::apply(slice, b);
    break;
  case element_type::marker:
    FF_drift::apply(slice, b);
    break;
  case element_type::instrument:
    FF_drift::apply(slice, b);
    break;
  case element_type::rcollimator:
    FF_drift::apply(slice, b);
    break;

  case element_type::foil:
    FF_foil::apply(slice, b);
    break;

  default:
    throw std::runtime_error(
        "FF_element::apply() unknown element type = " + elm.get_type_name() +
        ", element name = " + elm.get_name());
  }
}

template <class BUNCH> void apply(Lattice_element const &element, BUNCH &b) {
  apply(Lattice_element_slice(element), b);
}
}; // namespace FF_element

#endif // FF_ELEMENT_H
