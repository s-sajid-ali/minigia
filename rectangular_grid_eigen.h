#ifndef RECTANGULAR_GRID_EIGEN_H_
#define RECTANGULAR_GRID_EIGEN_H_

//#include "rectangular_grid_domain_eigen.h"
#include <unsupported/Eigen/CXX11/Tensor>


#if 0
typedef Eigen::Tensor<double, 3> EArray3d;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> EArray2dc;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> EArray1d;
#endif

template<typename T = double>
class Rectangular_grid_eigen
{
public:

    typedef Eigen::Tensor<T, 3, Eigen::RowMajor> EArray3d;

private:

#if 0
    Rectangular_grid_domain_eigen domain;
#endif

    EArray3d  points_;
    std::array<int, 3> shape_;

#if 0
    EArray2dc grid_points_2dc;
    EArray1d  grid_points_1d;
#endif

    double normalization;

public:

#if 0
    Rectangular_grid_eigen(
            std::array<double, 3> const & physical_size, 
            std::array<double, 3> const & physical_offset,
            std::array<int, 3>    const & grid_shape, 
            bool periodic_z )
        : domain(grid_shape, physical_size, physical_offset, periodic_z)
        , grid_points(grid_shape[0], grid_shape[1], grid_shape[2])
        //, grid_points_2dc(grid_shape[0], grid_shape[1])
        //, grid_points_1d(grid_shape[2])
        , normalization(1.0)
    { }

    Rectangular_grid_eigen(
            Rectangular_grid_domain_eigen const & domain )
        : domain(domain)
        , grid_points(domain.get_grid_shape()[0], domain.get_grid_shape()[1], domain.get_grid_shape()[2])
        //, grid_points_2dc(domain.get_grid_shape()[0], domain.get_grid_shape()[1])
        //, grid_points_1d(domain.get_grid_shape()[2])
        , normalization(1.0)
    { }
#endif

#if 0
    Rectangular_grid_domain_eigen const &
    get_domain() const
    { return domain; }

    Rectangular_grid_domain_eigen &
    get_domain()
    { return domain; }
#endif

    Rectangular_grid_eigen(std::array<int, 3> const & grid_shape, bool zero = true)
        : points_(grid_shape[0], grid_shape[1], grid_shape[2])
        , shape_(grid_shape)
        , normalization(1.0)
    { if (zero) set_zero(); }

    EArray3d const &
    get_grid_points() const
    { return points_; }

    EArray3d &
    get_grid_points()
    { return points_; }

    std::array<int, 3> const &
    shape() const
    { return shape_; }

    typename EArray3d::Scalar const &
    grid(Eigen::Index x, Eigen::Index y, Eigen::Index z) const
    { return points_(x, y, z); }

    typename EArray3d::Scalar &
    grid(Eigen::Index x, Eigen::Index y, Eigen::Index z)
    { return points_(x, y, z); }

    typename EArray3d::Scalar const *
    data(Eigen::Index x = 0, Eigen::Index y = 0, Eigen::Index z = 0) const
    { return points_.data() + x * shape_[1] * shape_[2] + y * shape_[2] + z; }

    void
    set_zero()
    { points_.setZero(); }

#if 0
    EArray2dc const &
    get_grid_points_2dc() const
    { return grid_points_2dc; }

    EArray2dc &
    get_grid_points_2dc()
    { return grid_points_2dc; }

    EArray1d const &
    get_grid_points_1d() const
    { return grid_points_1d; }

    EArray1d &
    get_grid_points_1d()
    { return grid_points_1d; }
#endif

    void
    set_normalization(double val)
    { normalization = val; }

    double
    get_normalization() const
    { return normalization; }

#if 0
    //
    // P.L. addition, Aug 3 2011
    //
    inline double get_interpolated(std::vector<double> const & location) const {
      return get_interpolated_coord(location[0], location[1], location[2]);
    }

    //
    // Tested for the 2D case only, where coordinate X is X and Y is Y. iz = 0, for this case.
    // Other cases needs work most likely. But I am always confused on the convention in
    // index name, Z vs X.  Depends on usage...
    //
    inline double get_interpolated_coord(double x, double y, double z) const 
    {
       // tri-linear interpolation
       int ix, iy, iz;
       double offx, offy, offz;
       domain.get_leftmost_indices_offsets(x, y, z, ix, iy, iz, offx, offy, offz);
       EArray3d const & a = grid_points;

       double val = 0.0;

       if ( (domain.get_grid_shape()[0] > 1) && 
            (domain.get_grid_shape()[1] > 1) && 
            (domain.get_grid_shape()[2] > 1) ) 
       {
           if ( (ix < 0) || (ix >= domain.get_grid_shape()[0] - 1) 
                   || (iy < 0) || (iy >= domain.get_grid_shape()[1] - 1) 
                   || (iz < 0) || (iz >= domain.get_grid_shape()[2] - 1) ) 
           {
              val = 0.0;
           } 
           else 
           {
              val = ( (1.0 - offx) * (1.0 - offy) * (1.0 - offz) * a(ix, iy, iz) 
                      + (1.0 - offx) * (1.0 - offy) * offz * a(ix, iy, iz + 1) 
                      + (1.0 - offx) * offy * (1.0 - offz) * a(ix, iy + 1, iz) 
                      + (1.0 - offx) * offy * offz * a(ix, iy + 1, iz + 1) 
                      + offx * (1.0 - offy) * (1.0 - offz) * a(ix + 1, iy, iz) 
                      + offx * (1.0 - offy) * offz * a(ix + 1, iy, iz + 1) 
                      + offx * offy * (1.0 - offz) * a(ix + 1, iy + 1, iz) 
                      + offx * offy * offz * a(ix + 1, iy + 1, iz + 1));
           }
       } 
       else if (domain.get_grid_shape()[0] == 1) 
       {
	       // 2D,  Y-Z plane
           if ( (iy < 0) || (iy >= domain.get_grid_shape()[1] - 1) 
                   || (iz < 0) || (iz >= domain.get_grid_shape()[2] - 1)) 
           { 
               val = 0.0; 
           } 
           else 
           {
               val = ((1.0 - offz) * (1.0 - offy) * a(ix, iy, iz) 
                       + offy * (1.0 - offz) * a(ix, iy + 1, iz) 
                       + (1.0 - offy) * offz * a(ix, iy, iz + 1) 
                       + offy * offz * a(ix, iy + 1, iz + 1));
           } 
       } 
       else if (domain.get_grid_shape()[1] == 1) 
       {
	       // 2D,  X-Z plane
           if ( (ix < 0) || (ix >= domain.get_grid_shape()[0] - 1) 
                   || (iz < 0) || (iz >= domain.get_grid_shape()[2] - 1) ) 
           { 
               val = 0.0; 
           } 
           else 
           {
               val = ((1.0 - offz) * (1.0 - offx) * a(ix, iy, iz) 
                       + offx * (1.0 - offz) * a(ix+1, iy, iz) 
                       + (1.0 - offx) * offz * a(ix, iy, iz + 1) 
                       + offx * offz * a(ix + 1, iy, iz + 1)); 
           } 
       }  
       else if (domain.get_grid_shape()[2] == 1) 
       {
	       // 2D,  X-Y plane
           if ( (ix < 0) || (ix >= domain.get_grid_shape()[0] - 1) 
                   || (iy < 0) || (iy >= domain.get_grid_shape()[1] - 1) ) 
           { 
               val = 0.0; 
           } 
           else 
           {
               val = ((1.0 - offy) * (1.0 - offx) * a(ix, iy, iz) 
                       + offx * (1.0 - offy) * a(ix+1, iy, iz) 
                       + (1.0 - offx) * offy * a(ix, iy+1, iz) 
                       + offx * offy * a(ix + 1, iy+1, iz)); 
           } 
       } 

       return val; 
    }
#endif
};

//typedef boost::shared_ptr<Rectangular_grid_eigen> Rectangular_grid_eigen_sptr; // syndoc:include

template<class T = double>
using Rectangular_grid_eigen_sptr = boost::shared_ptr<Rectangular_grid_eigen<T>>;

#endif /* RECTANGULAR_GRID_EIGEN_H_ */
