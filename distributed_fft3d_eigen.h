#ifndef DISTRIBUTED_FFT3D_EIGEN_H_
#define DISTRIBUTED_FFT3D_EIGEN_H_

#include <fftw3.h>
#include <fftw3-mpi.h>
#include <vector>
#include <string>
#include "boost/shared_ptr.hpp"
#include "multi_array_typedefs.h"
#include "commxx.h"
#include "rectangular_grid_eigen.h"

class Distributed_fft3d_eigen
{
private:
    fftw_plan plan, inv_plan;
    double *data;
    fftw_complex *workspace;
    int lower, upper;
    std::vector<int> uppers, lengths;
    int local_size_real, local_size_complex;
    std::array<int, 3> shape;
    bool have_local_data;
    Commxx_sptr comm_sptr;
    void
    calculate_uppers_lengths();
public:
    Distributed_fft3d_eigen(
            std::array<int, 3> const& shape, 
            Commxx_sptr comm_sptr,
            int planner_flags = FFTW_ESTIMATE,
            std::string const& wisdom_filename = "");
    Commxx_sptr
    get_comm_sptr();
    Commxx &
    get_comm();
    int
    get_lower() const;
    int
    get_upper() const;

    std::vector<int > const&
    get_uppers();

    std::vector<int > const&
    get_lengths();

    std::array<int, 3> const&
    get_shape() const;

    std::array<int, 3>
    get_padded_shape_real() const;

    std::array<int, 3 >
    get_padded_shape_complex() const;

    void
    transform(
            Rectangular_grid_eigen<double> & in, 
            Rectangular_grid_eigen<std::complex<double>> & out);

    void
    inv_transform(
            Rectangular_grid_eigen<std::complex<double>> & in,
            Rectangular_grid_eigen<double> & out);

    double
    get_roundtrip_normalization() const;

    ~Distributed_fft3d_eigen();
};

typedef boost::shared_ptr<Distributed_fft3d_eigen> Distributed_fft3d_eigen_sptr; // syndoc:include
#endif /* DISTRIBUTED_FFT3D_H_ */
