#ifndef COMPARE_KOKKOS_H
#define COMPARE_KOKKOS_H

#include <Kokkos_Core.hpp>
#include <iostream>

namespace
{
    KOKKOS_INLINE_FUNCTION
    bool kokkos_floating_point_equal(double a, double b, double tolerance)
    {
        if (std::abs(a) < tolerance) {
            return (std::abs(a - b) < tolerance);
        } else {
            return (std::abs((a - b) / a) < tolerance);
        }
    }
}

struct CompareBunch
{
    typedef Kokkos::View<double*[7], Kokkos::LayoutLeft> view_t;

    view_t a;
    view_t b;
    double tolerance;

    CompareBunch(view_t a_, view_t b_, double t) 
        : a(a_), b(b_), tolerance(t) 
    { }

    KOKKOS_INLINE_FUNCTION
    void operator()(const int p, int & ndiff) const
    {
        for (int i=0; i<7; ++i)
        {
            if (!kokkos_floating_point_equal(a(p,i), b(p,i), tolerance))
            {
                ++ndiff;
                break;
            }
        }
    }
};


inline bool
kokkos_check_equal(
        Kokkos::View<double*[7], Kokkos::LayoutLeft> a,
        Kokkos::View<double*[7], Kokkos::LayoutLeft> b,
        double tolerance)
{
    size_t na = a.extent(0);
    size_t nb = b.extent(0);

    if (na != nb)
    {
        std::cerr << "kokkos_check_equal: mismatched dimensions\n";
        return false;
    }


    int diff = 0;
    Kokkos::parallel_reduce(na, CompareBunch(a, b, tolerance), diff);

    return diff == 0;
}

#endif // COMPARE_H
