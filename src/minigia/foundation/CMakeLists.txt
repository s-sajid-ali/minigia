
add_library(minigia_foundation 
    four_momentum.cc 
    reference_particle.cc
    distribution.cc 
    )

target_link_libraries(minigia_foundation 
    ${GSL_CONFIG_LIBS}
    minigia_parallel_utils 
    lsexpr
    GSL::gsl
    GSL::gslcblas
    ${kokkos_libs}
    )

target_include_directories(minigia_parallel_utils INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/>  # <prefix>/include/mylib
)


install(TARGETS minigia_foundation DESTINATION lib)

install(FILES
    distribution.h
    four_momentum.h
    math_constants.h
    physical_constants.h
    reference_particle.h
    DESTINATION include/minigia/foundation)

add_subdirectory(tests)
