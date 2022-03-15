
add_executable(test_bunch_simulator_mpi test_bunch_simulator_mpi.cc ${test_main})
target_link_libraries(test_bunch_simulator_mpi PRIVATE minigia_simulation minigia_bunch ${testing_libs})

add_mpi_test(test_bunch_simulator_mpi 1)
add_mpi_test(test_bunch_simulator_mpi 2)
add_mpi_test(test_bunch_simulator_mpi 3)
add_mpi_test(test_bunch_simulator_mpi 4)
add_mpi_test(test_bunch_simulator_mpi 5)
add_mpi_test(test_bunch_simulator_mpi 6)

add_executable(test_bunch_simulator test_bunch_simulator.cc ${test_main})
target_link_libraries(test_bunch_simulator PRIVATE minigia_simulation minigia_bunch ${testing_libs})
add_mpi_test(test_bunch_simulator 1)



