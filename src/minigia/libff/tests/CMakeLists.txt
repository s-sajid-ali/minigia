add_executable(test_libff_elements test_libff_elements.cc ${test_main})
target_link_libraries(test_libff_elements 
    PRIVATE
    minigia_simulation 
    ${testing_libs})
add_mpi_test(test_libff_elements 1)
copy_file(fodo.madx test_libff_elements)

add_executable(test_madx_multipoles test_madx_multipoles.cc ${test_main})
target_link_libraries(test_madx_multipoles
    PRIVATE
    minigia_simulation 
    minigia_lattice
    minigia_foundation 
    ${testing_libs})
add_mpi_test(test_madx_multipoles 1)

add_executable(test_madx_elements test_madx_elements.cc ${test_main})
target_link_libraries(test_madx_elements
    PRIVATE
    minigia_simulation 
    minigia_lattice
    minigia_foundation 
    ${testing_libs})
add_mpi_test(test_madx_elements 1)

add_executable(test_foils test_foils.cc ${test_main})
target_link_libraries(test_foils
    PRIVATE
    minigia_simulation 
    minigia_lattice
    minigia_test_main
    ${testing_libs})
add_mpi_test(test_foils 1)

add_subdirectory(lattices)
