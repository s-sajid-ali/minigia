#include <cmath>
#include <iostream>

#include <minigia/foundation/reference_particle.hpp>
#include <minigia/utils/multi_array_typedefs.hpp>

#include "bunch.hpp"
#include "fixed_t_z_converter.hpp"

//*********************Fixed_t_z_zeroth*****************************************

void Fixed_t_z_zeroth::from_z_lab_to_t_bunch(Bunch &bunch) {
  double gamma = bunch.get_reference_particle().get_gamma();
  double beta = bunch.get_reference_particle().get_beta();
  double m = bunch.get_mass();
  double p_ref = bunch.get_reference_particle().get_momentum();

  auto particles = bunch.get_local_particles();
  auto s_particles = bunch.get_local_spectator_particles();

  int local_num = bunch.get_local_num();
  int local_s_num = bunch.get_local_spectator_num();

  double gb = -gamma * beta;

  bool exception_flag = false;

  Kokkos::parallel_for(local_num, [particles, gamma, beta, m, p_ref, gb,
                                   &exception_flag](const int i) {
    // z in beam rest frame
    particles(i, Bunch::z) = gb * particles(i, Bunch::cdt);

    // total momentum in accelerator frame
    double p = p_ref + particles(i, Bunch::dpop) * p_ref;

    // E/c in accelerator frame
    double Eoc = std::sqrt(p * p + m * m);

    // p_{x,y,z} in accelerator frame
    double px = particles(i, Bunch::xp) * p_ref;
    double py = particles(i, Bunch::yp) * p_ref;
    double pz2 = p * p - px * px - py * py;

    if (pz2 < 0.0)
      exception_flag = true;
    double pz = std::sqrt(pz2);

    // zp = pz/p_{ref}^{total}
    particles(i, Bunch::zp) = gamma * (pz - beta * Eoc) / p_ref;

    // n.b. in the zeroth approximation, the transformation from
    //      t' = gamma cdt to t' = 0
    //      is a no-op.
  });

  Kokkos::parallel_for(local_s_num, [s_particles, gamma, beta, m, p_ref, gb,
                                     &exception_flag](const int i) {
    // z in beam rest frame
    s_particles(i, Bunch::z) = gb * s_particles(i, Bunch::cdt);

    // total momentum in accelerator frame
    double p = p_ref + s_particles(i, Bunch::dpop) * p_ref;

    // E/c in accelerator frame
    double Eoc = std::sqrt(p * p + m * m);

    // p_{x,y,z} in accelerator frame
    double px = s_particles(i, Bunch::xp) * p_ref;
    double py = s_particles(i, Bunch::yp) * p_ref;
    double pz2 = p * p - px * px - py * py;

    if (pz2 < 0.0)
      exception_flag = true;
    double pz = std::sqrt(pz2);

    // zp = pz/p_{ref}^{total}
    s_particles(i, Bunch::zp) = gamma * (pz - beta * Eoc) / p_ref;

    // n.b. in the zeroth approximation, the transformation from
    //      t' = gamma cdt to t' = 0
    //      is a no-op.
  });

  if (exception_flag) {
    throw std::runtime_error(
        "Fixed_t_z_zeroth::fixed_z_to_fixed_t: particle has negative pz^2");
  }
}

void Fixed_t_z_zeroth::from_t_bunch_to_z_lab(Bunch &bunch) {}

void Fixed_t_z_zeroth::from_z_lab_to_t_lab(Bunch &bunch) {}
}

void Fixed_t_z_zeroth::from_t_lab_to_z_lab(Bunch &bunch) {}

void Fixed_t_z_zeroth::from_t_lab_to_t_bunch(Bunch &bunch) {}

void Fixed_t_z_zeroth::from_t_bunch_to_t_lab(Bunch &bunch) {}
