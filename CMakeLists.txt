cmake_minimum_required(VERSION 3.19)
project(Minigia LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
enable_testing()

# generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# and place it in the source directory
# https://stackoverflow.com/a/60910583
add_custom_target(
    copy-compile-commands ALL
    ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_BINARY_DIR}/compile_commands.json
        ${CMAKE_CURRENT_LIST_DIR}
    )

# project specific cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
# create the add_mpi_test macro
include(MPITest)
# create copy_file
include(CopyFile)

# MPI
find_package(MPI REQUIRED COMPONENTS C)
message(STATUS "MPI_VERSION is: " ${MPI_C_VERSION})
# mpi oversubscribe
if(${MPI_VERSION} MATCHES "OpenRTE.*")
    set(MPIEXEC_PREFLAGS ${MPIEXEC_PREFLAGS} "--oversubscribe")
endif()

# Boost, temporarily remove vendored version
#find_package(Boost)

# Kokkos
# Follow official kokkos guidelines
# https://github.com/kokkos/kokkos/wiki/Compiling#42-using-general-cmake-build-system
find_package(Kokkos REQUIRED)
set(kokkos_libs "Kokkos::kokkos")

# Ceral
find_package(cereal REQUIRED)

# Eigen
find_package(Eigen3 REQUIRED CONFIG)

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

# PETSc
# set root of location to find PETSc's pkg-config
set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)
find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)

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


# Catch2-v3
if(USE_EXTERNAL_CATCH2)
    find_package(Catch2 3 REQUIRED CONFIG)
    message(STATUS "Using external catch2-v3")
else()
    Include(FetchContent)
    FetchContent_Declare(
	Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG        v3.0.1
	)
    FetchContent_MakeAvailable(Catch2)
endif()
set(testing_libs "Catch2::Catch2;${kokkos_libs};MPI::MPI_C")
set(test_main "${CMAKE_SOURCE_DIR}/src/minigia/utils/catch_test_main.cc")
enable_testing()

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_subdirectory(src)
add_subdirectory(profiling)
add_subdirectory(examples)
