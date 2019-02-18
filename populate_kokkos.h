#ifndef POPULATE_KOKKOS_H
#define POPULATE_KOKKOS_H

#include "bunch_kokkos.h"
#include <random>

void
fill_random_gaussian(Bunch_kokkos& bunch)
{
    bunch.copy_to_host();

    auto particles = bunch.get_host_particles();
    auto local_num = bunch.get_local_num();

    // primes number 1e6, 2e6, 3e6, 4e6, 5e6
    std::seed_seq seed{ 15485863, 32452843, 49979687, 67867967, 86028121 };
    std::mt19937 generator(seed);
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (size_t  part = 0; part < local_num; ++part) {
        for (size_t index = 0; index < 6; ++index) {
            particles(part, index) = distribution(generator);
        }
    }

    bunch.copy_to_device();
}

#if 0
void
force_unit_covariance(Bunch_kokkos& bunch)
{
    auto particles = bunch.get_local_particles();
    auto local_num = bunch.get_local_num();
    auto particles6(particles.block(0, 0, local_num, 6));

    particles6.rowwise() -= particles6.colwise().mean();
    auto X((particles6.adjoint() * particles6) / particles6.rows());
    Eigen::Matrix<double, 6, 6> H(X.llt().matrixL());
    auto A(H.inverse());
    particles6 *= A.transpose();
}

void
set_covariance(Bunch_kokkos& bunch, Eigen::Matrix<double, 6, 6> const& covariance)
{
    Bunch_kokkos::Particles& particles(bunch.get_local_particles());
    auto local_num = bunch.get_local_num();
    auto particles6(particles.block(0, 0, local_num, 6));
    Eigen::Matrix<double, 6, 6> G(covariance.llt().matrixL());
    particles6 *= G.transpose();
}

// Second moments found at beginning of first space charge step
// in fodo_space_charge
Eigen::Matrix<double, 6, 6>
example_mom2()
{
    Eigen::Matrix<double, 6, 6> mom2;
    mom2(0, 0) = 3.219e-05;
    mom2(0, 1) = 1.141e-06;
    mom2(0, 2) = 0;
    mom2(0, 3) = 0;
    mom2(0, 4) = 0;
    mom2(0, 5) = 0;
    mom2(1, 0) = 1.141e-06;
    mom2(1, 1) = 7.152e-08;
    mom2(1, 2) = 0;
    mom2(1, 3) = 0;
    mom2(1, 4) = 0;
    mom2(1, 5) = 0;
    mom2(2, 0) = 0;
    mom2(2, 1) = 0;
    mom2(2, 2) = 7.058e-06;
    mom2(2, 3) = -3.226e-07;
    mom2(2, 4) = 0;
    mom2(2, 5) = 0;
    mom2(3, 0) = 0;
    mom2(3, 1) = 0;
    mom2(3, 2) = -3.226e-07;
    mom2(3, 3) = 1.564e-07;
    mom2(3, 4) = 0;
    mom2(3, 5) = 0;
    mom2(4, 0) = 0;
    mom2(4, 1) = 0;
    mom2(4, 2) = 0;
    mom2(4, 3) = 0;
    mom2(4, 4) = 0.0001643;
    mom2(4, 5) = -2.507e-09;
    mom2(5, 0) = 0;
    mom2(5, 1) = 0;
    mom2(5, 2) = 0;
    mom2(5, 3) = 0;
    mom2(5, 4) = -2.507e-09;
    mom2(5, 5) = 1e-08;
    return mom2;
}
#endif

void
set_particle_ids(Bunch_kokkos& bunch)
{
    bunch.copy_to_host();

    auto particles = bunch.get_host_particles();
    auto local_num = bunch.get_local_num();

    for (size_t part = 0; part < local_num; ++part) {
        particles(part, 6) = part;
    }

    bunch.copy_to_device();
}

#if 0
void
show_covariance(Bunch_kokkos const& bunch)
{
    auto particles = bunch.get_local_particles();
    auto local_num = bunch.get_local_num();
    auto particles6(particles.block(0, 0, local_num, 6));
    auto X((particles6.adjoint() * particles6) / particles6.rows());
    std::cout << "covariance =\n" << X << std::endl;
}
#endif

void
populate_gaussian(Bunch_kokkos& bunch)
{
    fill_random_gaussian(bunch);
    //force_unit_covariance(bunch);
    //set_covariance(bunch, example_mom2());
    set_particle_ids(bunch);
}

#endif // POPULATE_H
