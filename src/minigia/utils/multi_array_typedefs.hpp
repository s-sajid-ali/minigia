#ifndef MULTI_ARRAY_TYPEDEFS_H_
#define MULTI_ARRAY_TYPEDEFS_H_

#include <Kokkos_Core.hpp>

// column major, non-const arrays
typedef Kokkos::View<double *, Kokkos::LayoutLeft,
                     Kokkos::DefaultExecutionSpace>
    karray1d_dev;
typedef Kokkos::View<double **, Kokkos::LayoutLeft,
                     Kokkos::DefaultExecutionSpace>
    karray2d_dev;
typedef Kokkos::View<double ***, Kokkos::LayoutLeft,
                     Kokkos::DefaultExecutionSpace>
    karray3d_dev;

typedef karray1d_dev::HostMirror karray1d_hst;
typedef karray2d_dev::HostMirror karray2d_hst;
typedef karray3d_dev::HostMirror karray3d_hst;

typedef Kokkos::View<double *, Kokkos::LayoutLeft,
                     Kokkos::DefaultHostExecutionSpace>
    karray1d;
typedef Kokkos::View<double **, Kokkos::LayoutLeft,
                     Kokkos::DefaultHostExecutionSpace>
    karray2d;
typedef Kokkos::View<double ***, Kokkos::LayoutLeft,
                     Kokkos::DefaultHostExecutionSpace>
    karray3d;

// column major, const arrays
typedef Kokkos::View<const double *, Kokkos::LayoutLeft,
                     Kokkos::DefaultExecutionSpace>
    const_karray1d_dev;
typedef Kokkos::View<const double **, Kokkos::LayoutLeft,
                     Kokkos::DefaultExecutionSpace>
    const_karray2d_dev;
typedef Kokkos::View<const double ***, Kokkos::LayoutLeft,
                     Kokkos::DefaultExecutionSpace>
    const_karray3d_dev;

typedef const_karray1d_dev::HostMirror const_karray1d_hst;
typedef const_karray2d_dev::HostMirror const_karray2d_hst;
typedef const_karray3d_dev::HostMirror const_karray3d_hst;

typedef Kokkos::View<const double *, Kokkos::LayoutLeft,
                     Kokkos::DefaultHostExecutionSpace>
    const_karray1d;
typedef Kokkos::View<const double **, Kokkos::LayoutLeft,
                     Kokkos::DefaultHostExecutionSpace>
    const_karray2d;
typedef Kokkos::View<const double ***, Kokkos::LayoutLeft,
                     Kokkos::DefaultHostExecutionSpace>
    const_karray3d;

// row major, non-const arrays
typedef Kokkos::View<double *, Kokkos::LayoutRight,
                     Kokkos::DefaultExecutionSpace>
    karray1d_row_dev;
typedef Kokkos::View<double **, Kokkos::LayoutRight,
                     Kokkos::DefaultExecutionSpace>
    karray2d_row_dev;
typedef Kokkos::View<double ***, Kokkos::LayoutRight,
                     Kokkos::DefaultExecutionSpace>
    karray3d_row_dev;

typedef karray1d_row_dev::HostMirror karray1d_row_hst;
typedef karray2d_row_dev::HostMirror karray2d_row_hst;
typedef karray3d_row_dev::HostMirror karray3d_row_hst;

typedef Kokkos::View<double *, Kokkos::LayoutRight,
                     Kokkos::DefaultHostExecutionSpace>
    karray1d_row;
typedef Kokkos::View<double **, Kokkos::LayoutRight,
                     Kokkos::DefaultHostExecutionSpace>
    karray2d_row;
typedef Kokkos::View<double ***, Kokkos::LayoutRight,
                     Kokkos::DefaultHostExecutionSpace>
    karray3d_row;

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

typedef Kokkos::View<const double *, Kokkos::LayoutRight,
                     Kokkos::DefaultHostExecutionSpace>
    const_karray1d_row;
typedef Kokkos::View<const double **, Kokkos::LayoutRight,
                     Kokkos::DefaultHostExecutionSpace>
    const_karray2d_row;
typedef Kokkos::View<const double ***, Kokkos::LayoutRight,
                     Kokkos::DefaultHostExecutionSpace>
    const_karray3d_row;

// non-const int arrays
typedef Kokkos::View<int *, Kokkos::DefaultExecutionSpace::memory_space>
    karray1i_dev;

typedef karray1i_dev::HostMirror karray1i_hst;

// non-const uint32 arrays
typedef Kokkos::View<uint32_t *, Kokkos::DefaultExecutionSpace::memory_space>
    karray1ui32_dev;

typedef karray1ui32_dev::HostMirror karray1ui32_hst;

// non-const long arrays
typedef Kokkos::View<unsigned long *,
                     Kokkos::DefaultExecutionSpace::memory_space>
    karray1ul_dev;

typedef karray1ul_dev::HostMirror karray1ul_hst;

// row major, non-const int arrays
typedef Kokkos::View<int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
    karray1i_row_dev;
typedef Kokkos::View<int **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
    karray2i_row_dev;
typedef Kokkos::View<int ***, Kokkos::LayoutRight,
                     Kokkos::DefaultExecutionSpace>
    karray3i_row_dev;

typedef karray1i_row_dev::HostMirror karray1i_row_hst;
typedef karray2i_row_dev::HostMirror karray2i_row_hst;
typedef karray3i_row_dev::HostMirror karray3i_row_hst;

typedef Kokkos::View<int *, Kokkos::LayoutLeft,
                     Kokkos::DefaultHostExecutionSpace>
    karray1i_row;
typedef Kokkos::View<int **, Kokkos::LayoutLeft,
                     Kokkos::DefaultHostExecutionSpace>
    karray2i_row;
typedef Kokkos::View<int ***, Kokkos::LayoutLeft,
                     Kokkos::DefaultHostExecutionSpace>
    karray3i_row;

// row major, non-const complex arrays
typedef Kokkos::View<Kokkos::complex<double> *, Kokkos::LayoutRight,
                     Kokkos::DefaultExecutionSpace>
    karray1dc_row_dev;
typedef Kokkos::View<Kokkos::complex<double> **, Kokkos::LayoutRight,
                     Kokkos::DefaultExecutionSpace>
    karray2dc_row_dev;
typedef Kokkos::View<Kokkos::complex<double> ***, Kokkos::LayoutRight,
                     Kokkos::DefaultExecutionSpace>
    karray3dc_row_dev;

typedef karray1dc_row_dev::HostMirror karray1dc_row_hst;
typedef karray2dc_row_dev::HostMirror karray2dc_row_hst;
typedef karray3dc_row_dev::HostMirror karray3dc_row_hst;

typedef Kokkos::View<Kokkos::complex<double> *, Kokkos::LayoutLeft,
                     Kokkos::DefaultHostExecutionSpace>
    karray1dc_row;
typedef Kokkos::View<Kokkos::complex<double> **, Kokkos::LayoutLeft,
                     Kokkos::DefaultHostExecutionSpace>
    karray2dc_row;
typedef Kokkos::View<Kokkos::complex<double> ***, Kokkos::LayoutLeft,
                     Kokkos::DefaultHostExecutionSpace>
    karray3dc_row;

// atomic arrays
typedef Kokkos::View<double *, Kokkos::LayoutLeft,
                     Kokkos::MemoryTraits<Kokkos::Atomic>>
    karray1d_atomic_dev;

#endif /* MULTI_ARRAY_TYPEDEFS_H_ */
