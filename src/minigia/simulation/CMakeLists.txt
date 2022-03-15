if(${Kokkos_ENABLE_CUDA})
    set_source_files_properties(checkpoint_json.cc PROPERTIES COMPILE_FLAGS --host-only)
    set_source_files_properties(lattice_simulator_host.cc PROPERTIES COMPILE_FLAGS --host-only)
endif()

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

install(TARGETS minigia_simulation 
    DESTINATION ${LIB_INSTALL_DIR})

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
    DESTINATION ${INCLUDE_INSTALL_DIR}/minigia/simulation
    )

add_subdirectory(tests)
