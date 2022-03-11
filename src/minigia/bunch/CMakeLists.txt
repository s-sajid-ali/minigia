if(${Kokkos_ENABLE_CUDA})
    set_source_files_properties(populate_host.cc PROPERTIES COMPILE_FLAGS --host-only)
    set_source_files_properties(diagnostics_full2_host.cc PROPERTIES COMPILE_FLAGS --host-only)
endif()

add_library(minigia_bunch 
    bunch_particles.cc
    bunch.cc
    bunch_train.cc
    core_diagnostics.cc 
    diagnostics_worker.cc 
    diagnostics_file.cc 
    diagnostics_full2.cc 
    diagnostics_full2_host.cc 
    diagnostics_particles.cc 
    diagnostics_loss.cc
    diagnostics_bulk_track.cc 
    populate.cc 
    populate_global.cc 
    populate_host.cc 
    period.cc 
    )

target_link_libraries(minigia_bunch 
    minigia_parallel_utils 
    minigia_hdf5_utils 
    minigia_foundation 
    minigia_serialization 
    ${kokkos_libs}
    Eigen3::Eigen
    )

target_link_options(minigia_bunch
    PRIVATE ${LINKER_OPTIONS}
    )

install(TARGETS minigia_bunch DESTINATION ${LIB_INSTALL_DIR})
install(FILES
    bunch.hpp
    bunch_train.hpp
    core_diagnostics.hpp
    diagnostics.hpp
    diagnostics_basic.hpp
    diagnostics_full2.hpp
    diagnostics_track.hpp
    diagnostics_bulk_track.hpp
    diagnostics_bulk_spectator_track.hpp
    diagnostics_phase_space_density.hpp
    diagnostics_particles.hpp
    fixed_t_z_converter.hpp
    populate.hpp
    analysis.hpp
    period.hpp
    DESTINATION ${INCLUDE_INSTALL_DIR}/minigia/bunch)

add_subdirectory(tests)
