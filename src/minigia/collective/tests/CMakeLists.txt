
add_executable(test_space_charge_3d_fd test_space_charge_3d_fd.cc ${test_main})
target_link_libraries(test_space_charge_3d_fd PRIVATE minigia_collective minigia_simulation minigia_bunch ${testing_libs})
add_mpi_test(test_space_charge_3d_fd 1)
