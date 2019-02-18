#ifndef BUNCH_KOKKOS_H_
#define BUNCH_KOKKOS_H_

#include "commxx.h"
#include "compare_kokkos.h"
#include "reference_particle.h"
#include "restrict_extension.h"
#include <fstream>
#include <iomanip>
#include <iostream>

#if 1
//  J. Beringer et al. (Particle Data Group), PR D86, 010001 (2012) and 2013
//  partial update for the 2014 edition (URL: http://pdg.lbl.gov)
static const double proton_mass = 0.938272046; // Mass of proton [GeV/c^2]
static const double example_gamma = 10.0;
static const int proton_charge = 1; // Charge in units of e

static constexpr int particle_size_workaround = 7;
#endif

typedef Kokkos::View<double*[7], Kokkos::LayoutLeft> Particles;
typedef Particles::HostMirror HostParticles;

struct InitParts
{
    Particles p;
    const int o;
    InitParts(Particles p_, int o_) : p(p_), o(o_) { }

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const
    {
        const int idx = i + o;
        p(i, 0) = 1.0e-6 * idx;
        p(i, 1) = 1.1e-8 * idx;
        p(i, 2) = 1.3e-6 * idx;
        p(i, 3) = 1.4e-8 * idx;
        p(i, 4) = 1.5e-4 * idx;
        p(i, 5) = 1.5e-7 * idx;
        p(i, 6) = idx;
    }
};

class Bunch_kokkos
{
public:

    static const int x = 0;
    static const int xp = 1;
    static const int y = 2;
    static const int yp = 3;
    static const int z = 4;
    static const int zp = 5;
    static const int cdt = 4;
    static const int dpop = 5;
    static const int id = 6;
    static const int particle_size = particle_size_workaround;

private:

    Reference_particle reference_particle;
    long local_num, total_num;
    double real_num;

    Particles local_particles;
    HostParticles host_particles;

    Commxx_sptr comm_sptr;

public:
    Bunch_kokkos(long total_num, double real_num, int mpi_size = 1, int mpi_rank = 0)
        : reference_particle(proton_charge, proton_mass, example_gamma * proton_mass)
        , local_num(total_num / mpi_size)
        , total_num(total_num)
        , real_num(real_num)
        , local_particles("particles", local_num)
        , host_particles(Kokkos::create_mirror_view(local_particles))
        , comm_sptr(new Commxx())
    {
        Kokkos::parallel_for(local_num, InitParts(local_particles, mpi_rank * mpi_size));
        copy_to_host();
    }

    Reference_particle const& get_reference_particle() const
    { return reference_particle; }

    Particles get_local_particles() 
    { return local_particles; }

    HostParticles get_host_particles() 
    { return host_particles; }

    void copy_to_host() 
    { Kokkos::deep_copy(host_particles, local_particles); }

    void copy_to_device() 
    { Kokkos::deep_copy(local_particles, host_particles); }

    double get_mass() const 
    { return reference_particle.get_mass(); }

    long get_local_num() const 
    { return local_num; }

    long get_total_num() const 
    { return total_num; }

    double get_real_num() const 
    { return real_num; }

    Commxx_sptr get_comm_sptr() const 
    { return comm_sptr; }

    void set_arrays(double* RESTRICT& xa, double* RESTRICT& xpa,
                    double* RESTRICT& ya, double* RESTRICT& ypa,
                    double* RESTRICT& cdta, double* RESTRICT& dpopa)
    {
        xa    = &host_particles(0, Bunch_kokkos::x);
        xpa   = &host_particles(0, Bunch_kokkos::xp);
        ya    = &host_particles(0, Bunch_kokkos::y);
        ypa   = &host_particles(0, Bunch_kokkos::yp);
        cdta  = &host_particles(0, Bunch_kokkos::cdt);
        dpopa = &host_particles(0, Bunch_kokkos::dpop);
    }
};

inline bool
check_equal(Bunch_kokkos& b1, Bunch_kokkos& b2, double tolerance)
{
    if (b1.get_local_num() != b2.get_local_num()) {
        std::cerr << "check_equal: bunch 1 has " << b1.get_local_num()
                  << "local particles, ";
        std::cerr << "bunch 2 has " << b2.get_local_num() << "local particles"
                  << std::endl;
        return false;
    }

    return kokkos_check_equal(b1.get_local_particles(), b2.get_local_particles(), tolerance);
}

#endif /* BUNCH_H_ */
