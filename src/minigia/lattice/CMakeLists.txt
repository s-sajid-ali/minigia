if(${Kokkos_ENABLE_CUDA})
    set_source_files_properties(mx_parse.cc PROPERTIES COMPILE_FLAGS --host-only)
    set_source_files_properties(mx_tree.cc PROPERTIES COMPILE_FLAGS --host-only)
    set_source_files_properties(madx.cc PROPERTIES COMPILE_FLAGS --host-only)
endif()

add_library(minigia_lattice SHARED
    lattice_element.cc 
    lattice_element_processor.cc 
    lattice.cc
    lattice_element_slice.cc
    dynamic_lattice.cc
    lattice_tree.cc
    madx.cc 
    mx_expr.cc  
    mx_parse.cc  
    mx_tree.cc
    madx_reader.cc 
    )

target_link_libraries(minigia_lattice 
    minigia_foundation 
    minigia_bunch
    )

target_include_directories(minigia_lattice INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/>  # <prefix>/include/mylib
)

install(TARGETS minigia_lattice DESTINATION lib)

install(FILES
    lattice_element.hpp
    lattice_element_slice.hpp
    lattice.hpp
    madx.hpp
    madx_reader.hpp
    mx_expr.hpp
    mx_parse.hpp
    mx_tree.hpp
    DESTINATION include/minigia/lattice)

add_subdirectory(tests)
