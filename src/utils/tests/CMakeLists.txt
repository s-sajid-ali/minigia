
add_executable(test_command_line_arg test_command_line_arg.cc ${test_main})
target_link_libraries(test_command_line_arg PRIVATE minigia_command_line ${testing_libs})
add_test(test_command_line_arg test_command_line_arg)

add_executable(test_hdf5_file_mpi test_hdf5_file_mpi.cc ${test_main})
target_link_libraries(test_hdf5_file_mpi PRIVATE minigia_hdf5_utils ${testing_libs})

add_mpi_test(test_hdf5_file_mpi 1)
add_mpi_test(test_hdf5_file_mpi 2)
add_mpi_test(test_hdf5_file_mpi 3)
add_mpi_test(test_hdf5_file_mpi 4)

add_executable(test_hdf5_write_mpi test_hdf5_write_mpi.cc ${test_main})
target_link_libraries(test_hdf5_write_mpi PRIVATE minigia_hdf5_utils ${testing_libs})

add_mpi_test(test_hdf5_write_mpi 1)
add_mpi_test(test_hdf5_write_mpi 2)
add_mpi_test(test_hdf5_write_mpi 3)
add_mpi_test(test_hdf5_write_mpi 4)

add_executable(test_hdf5_append_mpi test_hdf5_append_mpi.cc ${test_main})
target_link_libraries(test_hdf5_append_mpi PRIVATE minigia_hdf5_utils ${testing_libs})

add_mpi_test(test_hdf5_append_mpi 1)
add_mpi_test(test_hdf5_append_mpi 2)
add_mpi_test(test_hdf5_append_mpi 3)
add_mpi_test(test_hdf5_append_mpi 4)

add_executable(test_hdf5_read_mpi test_hdf5_read_mpi.cc ${test_main})
target_link_libraries(test_hdf5_read_mpi PRIVATE minigia_hdf5_utils ${testing_libs})

add_mpi_test(test_hdf5_read_mpi 1)
add_mpi_test(test_hdf5_read_mpi 2)
add_mpi_test(test_hdf5_read_mpi 3)
add_mpi_test(test_hdf5_read_mpi 4)


add_executable(test_kokkos test_kokkos.cc ${test_main})
target_link_libraries(test_kokkos PRIVATE ${testing_libs})

add_mpi_test(test_kokkos 1)

add_executable(test_commxx_mpi test_commxx_mpi.cc ${test_main})
target_link_libraries(test_commxx_mpi PRIVATE minigia_parallel_utils ${testing_libs})
add_mpi_test(test_commxx_mpi 1)
add_mpi_test(test_commxx_mpi 2)
add_mpi_test(test_commxx_mpi 3)
add_mpi_test(test_commxx_mpi 4)

add_executable(test_commxx_serdes test_commxx_serdes.cc ${test_main})
target_link_libraries(test_commxx_serdes PRIVATE minigia_parallel_utils ${testing_libs})
add_mpi_test(test_commxx_serdes 1)
add_mpi_test(test_commxx_serdes 2)
add_mpi_test(test_commxx_serdes 3)
add_mpi_test(test_commxx_serdes 4)
