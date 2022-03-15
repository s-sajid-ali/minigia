
add_executable(test_lattice test_lattice.cc ${test_main})
target_link_libraries(test_lattice PRIVATE minigia_lattice ${testing_libs})
add_mpi_test(test_lattice 1)


add_executable(test_mx_expr test_mx_expr.cc ${test_main})
target_link_libraries(test_mx_expr PRIVATE minigia_lattice ${testing_libs})
add_mpi_test(test_mx_expr 1)


add_executable(test_madx_parser test_madx_parser.cc ${test_main})
target_link_libraries(test_madx_parser PRIVATE minigia_lattice ${testing_libs})
add_mpi_test(test_madx_parser 1)
copy_file(foo.dbx test_madx_parser)

add_executable(test_dynamic_lattice test_dynamic_lattice.cc ${test_main})
target_link_libraries(test_dynamic_lattice PRIVATE minigia_lattice ${testing_libs})
add_mpi_test(test_dynamic_lattice 1)


