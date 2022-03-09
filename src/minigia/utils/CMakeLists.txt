set(minigia_parallel_utils_src
    parallel_utils.cc 
    commxx.cc 
    logger.cc
    simple_timer.cc
    base64.cpp
    )

add_library(minigia_parallel_utils SHARED ${minigia_parallel_utils_src})
target_link_libraries(minigia_parallel_utils MPI::MPI_C ${kokkos_libs})

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
    ${HDF5_LIBRARIES} 
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
    command_line_arg.h
    commxx.h
    commxx_divider.h
    fast_int_floor.h
    floating_point.h
    gsvector.h
    hdf5_misc.h
    hdf5_file.h
    hdf5_serial_writer.h
    hdf5_writer.h
    multi_array_typedefs.h
    parallel_utils.h
    simple_timer.h
    cereal.h
    cereal_files.h
    digits.h
    logger.h
    lsexpr.h
    DESTINATION include/minigia/utils)

add_subdirectory(tests)

