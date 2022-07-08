add_library(minigia_collective SHARED 
    space_charge_3d_fd.cc
    space_charge_3d_fd_utils.cc
    space_charge_3d_fd_alias.cc 
    deposit.cc
    deposit_new.cc
    )

target_link_libraries(minigia_collective
    minigia_bunch
    minigia_lattice
    minigia_parallel_utils
    PkgConfig::PETSC
    )

target_include_directories(minigia_foundation INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/>  # <prefix>/include/mylib
  )

install(TARGETS minigia_collective DESTINATION lib)
install(FILES
    deposit.hpp
    deposit_new.hpp
    rectangular_grid_domain.hpp
    rectangular_grid.hpp
    DESTINATION ${INCLUDE_INSTALL_DIR}/minigia/collective)

add_subdirectory(tests)
