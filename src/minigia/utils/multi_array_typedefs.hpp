#ifndef MULTI_ARRAY_TYPEDEFS_H_
#define MULTI_ARRAY_TYPEDEFS_H_

#include <Kokkos_Core.hpp>

// column major, non-const arrays
typedef Kokkos::View<double *, Kokkos::LayoutLeft> karray1d_dev;
typedef Kokkos::View<double **, Kokkos::LayoutLeft> karray2d_dev;
typedef Kokkos::View<double ***, Kokkos::LayoutLeft> karray3d_dev;

typedef karray1d_dev::HostMirror karray1d_hst;
typedef karray2d_dev::HostMirror karray2d_hst;
typedef karray3d_dev::HostMirror karray3d_hst;

typedef karray1d_hst karray1d;
typedef karray2d_hst karray2d;
typedef karray3d_hst karray3d;

// column major, const arrays
typedef Kokkos::View<const double *, Kokkos::LayoutLeft> const_karray1d_dev;
typedef Kokkos::View<const double **, Kokkos::LayoutLeft> const_karray2d_dev;
typedef Kokkos::View<const double ***, Kokkos::LayoutLeft> const_karray3d_dev;

typedef const_karray1d_dev::HostMirror const_karray1d_hst;
typedef const_karray2d_dev::HostMirror const_karray2d_hst;
typedef const_karray3d_dev::HostMirror const_karray3d_hst;

typedef const_karray1d_hst const_karray1d;
typedef const_karray2d_hst const_karray2d;
typedef const_karray3d_hst const_karray3d;

// row major, non-const arrays
typedef Kokkos::View<double *, Kokkos::LayoutRight> karray1d_row_dev;
typedef Kokkos::View<double **, Kokkos::LayoutRight> karray2d_row_dev;
typedef Kokkos::View<double ***, Kokkos::LayoutRight> karray3d_row_dev;

typedef karray1d_row_dev::HostMirror karray1d_row_hst;
typedef karray2d_row_dev::HostMirror karray2d_row_hst;
typedef karray3d_row_dev::HostMirror karray3d_row_hst;

typedef karray1d_row_hst karray1d_row;
typedef karray2d_row_hst karray2d_row;
typedef karray3d_row_hst karray3d_row;

// row major, const arrays
typedef Kokkos::View<const double *, Kokkos::LayoutRight>
    const_karray1d_row_dev;
typedef Kokkos::View<const double **, Kokkos::LayoutRight>
    const_karray2d_row_dev;
typedef Kokkos::View<const double ***, Kokkos::LayoutRight>
    const_karray3d_row_dev;

typedef const_karray1d_row_dev::HostMirror const_karray1d_row_hst;
typedef const_karray2d_row_dev::HostMirror const_karray2d_row_hst;
typedef const_karray3d_row_dev::HostMirror const_karray3d_row_hst;

typedef const_karray1d_row_hst const_karray1d_row;
typedef const_karray2d_row_hst const_karray2d_row;
typedef const_karray3d_row_hst const_karray3d_row;

// row major, non-const int arrays
typedef Kokkos::View<int *, Kokkos::LayoutRight> karray1i_row_dev;
typedef Kokkos::View<int **, Kokkos::LayoutRight> karray2i_row_dev;
typedef Kokkos::View<int ***, Kokkos::LayoutRight> karray3i_row_dev;

typedef karray1i_row_dev::HostMirror karray1i_row_hst;
typedef karray2i_row_dev::HostMirror karray2i_row_hst;
typedef karray3i_row_dev::HostMirror karray3i_row_hst;

typedef karray1i_row_hst karray1i_row;
typedef karray2i_row_hst karray2i_row;
typedef karray3i_row_hst karray3i_row;

// row major, non-const complex arrays
typedef Kokkos::View<Kokkos::complex<double> *, Kokkos::LayoutRight>
    karray1dc_row_dev;
typedef Kokkos::View<Kokkos::complex<double> **, Kokkos::LayoutRight>
    karray2dc_row_dev;
typedef Kokkos::View<Kokkos::complex<double> ***, Kokkos::LayoutRight>
    karray3dc_row_dev;

typedef karray1dc_row_dev::HostMirror karray1dc_row_hst;
typedef karray2dc_row_dev::HostMirror karray2dc_row_hst;
typedef karray3dc_row_dev::HostMirror karray3dc_row_hst;

typedef karray1dc_row_hst karray1dc_row;
typedef karray2dc_row_hst karray2dc_row;
typedef karray3dc_row_hst karray3dc_row;

// atomic arrays
typedef Kokkos::View<double *, Kokkos::LayoutLeft,
                     Kokkos::MemoryTraits<Kokkos::Atomic>>
    karray1d_atomic_dev;

#endif /* MULTI_ARRAY_TYPEDEFS_H_ */
