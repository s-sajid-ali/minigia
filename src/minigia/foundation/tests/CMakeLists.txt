
add_executable(test_four_momentum test_four_momentum.cc ${test_main})
target_link_libraries(test_four_momentum PRIVATE minigia_foundation ${testing_libs})
add_mpi_test(test_four_momentum 1)

add_executable(test_reference_particle test_reference_particle.cc ${test_main})
target_link_libraries(test_reference_particle PRIVATE minigia_foundation ${testing_libs})
add_mpi_test(test_reference_particle 1)
