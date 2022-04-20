set(minigia_parallel_utils_src
    parallel_utils.cc 
    commxx.cc 
    logger.cc
    simple_timer.cc
    base64.cpp
    )

add_library(minigia_parallel_utils SHARED ${minigia_parallel_utils_src})
target_link_libraries(minigia_parallel_utils MPI::MPI_C ${kokkos_libs} cereal::cereal)

target_include_directories(minigia_parallel_utils INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/>  # <prefix>/include/mylib
)

add_library(minigia_test_main STATIC catch_test_main.cc)
target_link_libraries(minigia_test_main PRIVATE Catch2::Catch2 ${kokkos_libs} MPI::MPI_C)

add_library(minigia_serialization SHARED cereal_files.cc)
target_link_libraries(minigia_serialization minigia_parallel_utils)

add_library(minigia_hdf5_utils SHARED hdf5_file.cc hdf5_misc.cc)
target_link_libraries(minigia_hdf5_utils 
    HDF5::HDF5
    minigia_parallel_utils
    ${kokkos_libs}
    )

add_library(minigia_command_line SHARED command_line_arg.cc)
target_link_libraries(minigia_command_line ${kokkos_libs})

add_library(lsexpr SHARED lsexpr.cc)
target_link_libraries(lsexpr ${kokkos_libs})

install(TARGETS
    minigia_parallel_utils
    minigia_hdf5_utils
    minigia_serialization
    minigia_command_line
    lsexpr
    DESTINATION lib)

install(FILES 
    command_line_arg.hpp
    commxx.hpp
    fast_int_floor.hpp
    floating_point.hpp
    gsvector.hpp
    hdf5_misc.hpp
    hdf5_file.hpp
    hdf5_serial_writer.hpp
    kokkos_types.hpp
    multi_array_typedefs.hpp
    minigia_json.hpp
    parallel_utils.hpp
    simple_timer.hpp
    cereal.hpp
    cereal_files.hpp
    digits.hpp
    logger.hpp
    lsexpr.hpp
    restrict_extension.hpp
    string_utils.hpp
    DESTINATION include/minigia/utils)

add_subdirectory(tests)

