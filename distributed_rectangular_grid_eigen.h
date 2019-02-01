#ifndef DISTRIBUTED_RECTANGULAR_GRID_EIGEN_H_
#define DISTRIBUTED_RECTANGULAR_GRID_EIGEN_H_

#include "rectangular_grid_eigen.h"
#include "commxx.h"
#include <unsupported/Eigen/CXX11/Tensor>

template<typename T = double>
class Distributed_rectangular_grid_eigen
{
private:

    int lower, upper;
    int lower_guard, upper_guard;

    std::array<int, 3> shape_;
    Rectangular_grid_eigen<T> grid_;

    //std::vector<int> uppers, lengths;
    //double normalization;
    Commxx_sptr comm_sptr;

#if 0
    void construct_hockney(int lower, int upper, std::vector<int > const & array_shape);
    void construct_rectangular(int lower, int upper, std::vector<int > const & array_shape);
    void calculate_uppers_lengths();
#endif


public:

    Distributed_rectangular_grid_eigen(
            std::array<int, 3> const & grid_shape,
            int lower, int upper, bool periodic,
            Commxx_sptr comm_sptr,
            bool zero = true)
        : lower( lower )
        , upper( upper )
        , lower_guard( (lower == 0 && !periodic) ? 0 : lower - 1 )
        , upper_guard( (upper == grid_shape[0] && !periodic) ? grid_shape[0] : upper + 1 )
        , shape_( {upper_guard - lower_guard + 1, grid_shape[1], grid_shape[2]} )
        , grid_ ( shape_, zero )
        , comm_sptr( comm_sptr )
    {
    }

    T const &
    grid(Eigen::Index x, Eigen::Index y, Eigen::Index z) const
    { return grid_.grid(x-lower_guard+1, y, z); }

    T &
    grid(Eigen::Index x, Eigen::Index y, Eigen::Index z)
    { return grid_.grid(x-lower_guard+1, y, z); }

    T const *
    data(Eigen::Index x, Eigen::Index y, Eigen::Index z) const
    { return grid_.data(x-lower_guard+1, y, z); }

    T const *
    data() const
    { return grid_.data(lower-lower_guard+1, 0, 0); }

    std::array<int, 3> const &
    shape() const
    { return shape_; }

    int get_lower() const
    { return lower; }

    int get_upper() const
    { return upper; }

    int get_lower_guard() const
    { return lower_guard; }

    int get_upper_guard() const
    { return upper_guard; }

    void set_normalization(double val)
    { grid_.set_normalization(val); }

    double get_normalization() const
    { return grid_.get_normalization(); }

    Commxx const & get_comm() const
    { return *comm_sptr; }

    Commxx_sptr get_comm_sptr() const
    { return comm_sptr; }


#if 0
    Distributed_rectangular_grid_eigen(
            std::array<double, 3> const & physical_size,
            std::array<double, 3> const & physical_offset,
            std::array<int, 3>    const & grid_shape, 
            bool periodic, 
            int lower, int upper, 
            Commxx_sptr comm_sptr, 
            std::string const solver="hockney");

    Distributed_rectangular_grid_eigen(
            Rectangular_grid_domain_eigen const & domain,
            int lower, int upper, 
            Commxx_sptr comm_sptr, 
            std::string const solver="hockney");

    Distributed_rectangular_grid_eigen(
            Rectangular_grid_domain_eigen const & domain,
            int lower, int upper, 
            std::vector<int > const & padded_shape,
            Commxx_sptr comm_sptr);

    Rectangular_grid_domain_eigen const & 
    get_domain() const
    { return domain; }

    Rectangular_grid_domain_eigen & 
    get_domain()
    { return domain; }

    std::vector<int> const & 
    get_uppers()
    { calculate_uppers_lengths(); return uppers; }

    std::vector<int> const & 
    get_lengths()
    { calculate_uppers_lengths(); return lengths; }

    EArray3d const &
    get_grid_points() const
    { return grid_points; }

    EArray3d &
    get_grid_points()
    { return grid_points; }

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
    void
    fill_guards();
#endif
};

#if 0
typedef boost::shared_ptr<Distributed_rectangular_grid_eigen >
        Distributed_rectangular_grid_eigen_sptr;  // syndoc:include
#endif

template<class T = double>
using Distributed_rectangular_grid_eigen_sptr = boost::shared_ptr<Distributed_rectangular_grid_eigen<T>>;

#endif /* DISTRIBUTED_RECTANGULAR_GRID_EIGEN_H_ */
