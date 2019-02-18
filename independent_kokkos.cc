#include <chrono>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

//#include <mpi.h>
#if defined(_OPENMP)
#include <omp.h>
#endif

#include <Kokkos_Core.hpp>

#if 0
#define GSV_AVX 1
#include "gsvector_v2.h"
#endif

#include "bunch_kokkos.h"
#include "populate_kokkos.h"
//#include "bunch_data_paths.h"


const int particles_per_rank = 100000000;
const double real_particles = 1.0e12;

const double dummy_length = 2.1;
const double dummy_reference_time = 0.345;

struct libff_drift
{
    libff_drift() {}
    double Length() const { return dummy_length; }
    double getReferenceTime() const { return dummy_reference_time; }
};

using propagator_t = void(*)(Bunch_kokkos &, libff_drift &);

//#pragma omp declare simd
template<typename T>
KOKKOS_INLINE_FUNCTION
void libff_drift_unit
  (T & x, T & y, T & cdt, T const& xp, T const& yp, T const& dpop,
   double length, double reference_momentum, double m, double reference_cdt)
{
    T uni(1.0);
    //T sig((0.0<length) - (length<0.0));

    T vl(length);
    T vm(m);
    T vrm(reference_momentum);
    T vrc(reference_cdt);

    T dp = dpop + uni;
    T inv_npz = uni / sqrt(dp * dp - xp * xp - yp * yp);
    T lxpr = xp * vl * inv_npz;
    T lypr = yp * vl * inv_npz;
    T D2 = lxpr * lxpr + vl * vl + lypr * lypr;
    T p = dp * vrm;
    T E2 = p * p + vm * vm;
    //T beta2 = p*p / E2;
    T ibeta2 = E2 / (p * p);
    x = x + lxpr;
    y = y + lypr;
    cdt = cdt + sqrt(D2 * ibeta2) - vrc;
    //cdt = cdt + sig * sqrt(D2 * ibeta2) - vrc;
}


void
propagate_orig(Bunch_kokkos& bunch, libff_drift& thelibff_drift)
{
    bunch.copy_to_host();

    auto local_num = bunch.get_local_num();
    auto particles = bunch.get_host_particles();

    auto length = thelibff_drift.Length();
    auto reference_momentum = bunch.get_reference_particle().get_momentum();
    auto m = bunch.get_mass();
    auto reference_time = thelibff_drift.getReferenceTime();

    for (int part = 0; part < local_num; ++part) {
        auto dpop(particles(part, Bunch_kokkos::dpop));
        auto xp(particles(part, Bunch_kokkos::xp));
        auto yp(particles(part, Bunch_kokkos::yp));
        auto inv_npz = 1.0 / sqrt((dpop + 1.0) * (dpop + 1.0) - xp * xp - yp * yp);
        auto lxpr = xp * length * inv_npz;
        auto lypr = yp * length * inv_npz;
        auto D = sqrt(lxpr * lxpr + length * length + lypr * lypr);
        auto p = reference_momentum + dpop * reference_momentum;
        auto E = sqrt(p * p + m * m);
        auto beta = p / E;
        auto x(particles(part, Bunch_kokkos::x));
        auto y(particles(part, Bunch_kokkos::y));
        auto cdt(particles(part, Bunch_kokkos::cdt));

        x += lxpr;
        y += lypr;
        cdt += D / beta - reference_time;

        particles(part, Bunch_kokkos::x) = x;
        particles(part, Bunch_kokkos::y) = y;
        particles(part, Bunch_kokkos::cdt) = cdt;
    }

    bunch.copy_to_device();
}

void
propagate_double(Bunch_kokkos& bunch, libff_drift& thelibff_drift)
{
    bunch.copy_to_host();

    auto local_num = bunch.get_local_num();
    const auto length = thelibff_drift.Length();
    const auto reference_momentum = bunch.get_reference_particle().get_momentum();
    const auto m = bunch.get_mass();
    const auto reference_time = thelibff_drift.getReferenceTime();

    double *RESTRICT xa, *RESTRICT xpa, *RESTRICT ya, *RESTRICT ypa,
        *RESTRICT cdta, *RESTRICT dpopa;
    bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

    for (int part = 0; part < local_num; ++part) 
    {
        auto x(xa[part]);
        auto xp(xpa[part]);
        auto y(ya[part]);
        auto yp(ypa[part]);
        auto cdt(cdta[part]);
        auto dpop(dpopa[part]);

        libff_drift_unit(x, y, cdt, xp, yp, dpop, length, reference_momentum, m, reference_time);

        xa[part] = x;
        ya[part] = y;
        cdta[part] = cdt;
    }

    bunch.copy_to_device();
}

struct PropDevice
{
    Particles p;
    double l, ref_p, m, ref_t;

    PropDevice(Particles p_, double l, double ref_p, double m, double ref_t) 
        : p(p_), l(l), ref_p(ref_p), m(m), ref_t(ref_t) { }

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const
    {
        libff_drift_unit(
                p(i, 0 /*Bunch_kokkos::x*/),
                p(i, 2 /*Bunch_kokkos::y*/),
                p(i, 4 /*Bunch_kokkos::cdt*/),
                p(i, 1 /*Bunch_kokkos::xp*/),
                p(i, 3 /*Bunch_kokkos::yp*/),
                p(i, 5 /*Bunch_kokkos::dpop*/),
                l, ref_p, m, ref_t);
    }
};

void
propagate_device(Bunch_kokkos& bunch, libff_drift& thelibff_drift)
{
    auto particles = bunch.get_local_particles();
    auto local_num = bunch.get_local_num();

    const auto length = thelibff_drift.Length();
    const auto reference_momentum = bunch.get_reference_particle().get_momentum();
    const auto m = bunch.get_mass();
    const auto reference_time = thelibff_drift.getReferenceTime();

    Kokkos::parallel_for(
            local_num, 
            PropDevice(particles, length, reference_momentum, m, reference_time) );
}

void
run_check( propagator_t propagator,
           const char * name,
           libff_drift & drift,
           int size, int rank)
{
    const double tolerance = 1.0e-14;
    const int num_test = 104 * size;
    const double real_num = 1.0e12;

    Bunch_kokkos b1(num_test * size, real_num, size, rank);
    Bunch_kokkos b2(num_test * size, real_num, size, rank);

    propagate_orig(b1, drift);
    propagator(b2, drift);

    if (!check_equal(b1, b2, tolerance)) {
        std::cerr << "run_check failed for " << name << std::endl;
    }
}

double
do_timing( propagator_t propagator,
           const char * name,
           Bunch_kokkos & bunch,
           libff_drift & drift,
           double ref_t,
           const int rank )
{
    double t = 0;
    const int num_runs = 4;
    auto best_time = std::numeric_limits<double>::max();
    std::vector<double> times(num_runs);

    for (size_t i = 0; i < num_runs; ++i) 
    {
        //populate_gaussian(bunch);

        const auto start = std::chrono::high_resolution_clock::now();
        propagator(bunch, drift);
        const auto end = std::chrono::high_resolution_clock::now();
        const auto time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start) .count();

        if (time < best_time) best_time = time;
        times.at(i) = time;
        t += time;
    }

    if (rank == 0) {
        std::cout << name << " best time = " << best_time << std::endl;
    }

    if (ref_t > 0.0) {
        if (rank == 0) {
            std::cout << name << " speedup = " << ref_t / best_time
                      << std::endl;
        }
    }

    return best_time;
}



void run()
{
    int error, rank, size;
    error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    error = MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (error) 
    {
        std::cerr << "MPI error" << std::endl;
        return;
    }

    // bunch
    Bunch_kokkos bunch(particles_per_rank, real_particles, 1, 0);

    // populate gaussian
    std::cout << "populating..."; std::flush(std::cout);
    populate_gaussian(bunch);
    std::cout << "done.\n";

    // element
    libff_drift drift;

    // original, always runs on host
    auto ori_t = do_timing(propagate_orig, "orig", bunch, drift, 0.0, rank);

    // serial, runs on host
    run_check(propagate_double, "double", drift, size, rank);
    do_timing(propagate_double, "double", bunch, drift, ori_t, rank);

    // parallel, runs on available devices (GPU or OpenMP)
    run_check(propagate_device, "device", drift, size, rank);
    do_timing(propagate_device, "device", bunch, drift, ori_t, rank);
}



int main(int argc, char ** argv)
{
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    run();

    Kokkos::finalize();
    MPI_Finalize();
 
    return 0;
}



