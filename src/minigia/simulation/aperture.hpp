#ifndef APERTURE_H_
#define APERTURE_H_

#include <minigia/bunch/bunch.hpp>
#include <minigia/lattice/lattice_element_slice.hpp>

const double default_aperture_radius = 1000.0;

void
apply_circular_aperture(Bunch & bunch, Lattice_element_slices & slices);


#endif /* APERTURE_H_ */
