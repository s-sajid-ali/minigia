
add_executable(test_bunch test_bunch.cc ${test_main})
target_link_libraries(test_bunch PRIVATE minigia_bunch ${testing_libs})
add_mpi_test(test_bunch 1)

add_executable(test_bunch_particles test_bunch_particles.cc ${test_main})
target_link_libraries(test_bunch_particles PRIVATE minigia_bunch ${testing_libs})
add_mpi_test(test_bunch_particles 1)

add_executable(test_bunch_traits test_bunch_traits.cc ${test_main})
target_link_libraries(test_bunch_traits PRIVATE minigia_bunch ${testing_libs})
add_mpi_test(test_bunch_traits 1)
