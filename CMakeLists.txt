cmake_minimum_required(VERSION 3.19)
project(Minigia LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
enable_testing()

# project specific cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
# create the add_mpi_test macro
include(MPITest)

# MPI
find_package(MPI REQUIRED COMPONENTS C)
message(MPI_VERSION is: ${MPI_VERSION})
message(STATUS "MPI_INCLUDE_PATH : ${MPI_INCLUDE_PATH}")
message(STATUS "MPI_LIBRARIES : ${MPI_LIBRARIES}")
# mpi oversubscribe
if(${MPI_VERSION} MATCHES "OpenRTE.*")
    set(MPIEXEC_PREFLAGS ${MPIEXEC_PREFLAGS} "--oversubscribe")
endif()

# Kokkos
# Follow official kokkos guidelines
# https://github.com/kokkos/kokkos/wiki/Compiling#42-using-general-cmake-build-system
find_package(Kokkos REQUIRED)
set(kokkos_libs "Kokkos::kokkos")

# Ceral
find_package(cereal REQUIRED)

# HDF5
set(HDF5_PREFER_PARALLEL true)
find_package(HDF5)
list(APPEND extra_includes ${HDF5_INCLUDE_DIRS})

if(HDF5_IS_PARALLEL)
    message(STATUS "Using parallel hdf5")
    add_definitions(-DUSE_PARALLEL_HDF5)
else(HDF5_IS_PARALLEL)
    message(STATUS "Using serial hdf5")
endif()

# GSL
find_package(GSL REQUIRED)

# GSVector settings
if(GSV STREQUAL "SSE")
    message(STATUS "using SSE for gsvector")
    add_definitions(-DGSV_SSE)
elseif(GSV STREQUAL "AVX")
    message(STATUS "using AVX for gsvector")
    add_definitions(-DGSV_AVX)
elseif(GSV STREQUAL "AVX512")
    message(STATUS "using AVX512 for gsvector")
    add_definitions(-DGSV_AVX512)
else()
    message(STATUS "using double for gsvector")
endif()


add_subdirectory(src)
add_subdirectory(examples)
