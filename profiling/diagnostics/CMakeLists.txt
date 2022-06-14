add_executable(diagnostics_mean_profile diagnostics_mean_profile.cc)
target_link_libraries(diagnostics_mean_profile
    minigia_bunch 
    minigia_lattice
    minigia_simulation
    minigia_command_line
    )

add_executable(diagnostics_mean_stddev_profile diagnostics_mean_stddev_profile.cc)
target_link_libraries(diagnostics_mean_stddev_profile
    minigia_bunch 
    minigia_lattice
    minigia_simulation
    minigia_command_line
    )

