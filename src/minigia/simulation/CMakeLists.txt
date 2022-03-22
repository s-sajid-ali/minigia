add_library(minigia_simulation SHARED
    operator.cc 
    step.cc 
    stepper.cc
    independent_stepper_elements.cc
    split_operator_stepper.cc 
    split_operator_stepper_elements.cc
    propagator.cc
    independent_operation.cc
    lattice_simulator.cc
    lattice_simulator_host.cc
    operation_extractor.cc 
    bunch_simulator.cc
    checkpoint.cc
    checkpoint_json.cc
    )

target_link_libraries(minigia_simulation 
    minigia_foundation 
    minigia_bunch
    minigia_lattice 
    minigia_libff
    )

target_include_directories(minigia_foundation INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/>  # <prefix>/include/mylib
  )

install(TARGETS minigia_simulation 
    DESTINATION lib)

install(FILES
    aperure.hpp
    aperture_operation.hpp
    aperture_operation_extractor.hpp
    bunch_simulator.hpp
    bunch_simulator_impl.hpp
    checkpoint.hpp
    collective_operator_options.hpp
    dummy_collective_operator.hpp
    independent_operation.hpp
    independent_stepper_elements.hpp
    lattice_simulator.hpp
    lattice_simulator_host.hpp
    operation_extractor.hpp
    operator.hpp
    populate_stationary.hpp
    propagate_actions.hpp
    propagator.hpp
    resume.hpp
    split_operator_stepper.hpp
    split_operator_stepper_elements.hpp
    step.hpp
    stepper.hpp
    DESTINATION include/minigia/simulation
    )

add_subdirectory(tests)
