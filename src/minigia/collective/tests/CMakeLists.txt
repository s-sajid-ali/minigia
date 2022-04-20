
add_executable(test_deposit test_deposit.cc ${test_main})
target_link_libraries(test_deposit PRIVATE minigia_collective ${testing_libs})
add_mpi_test(test_deposit 1)
