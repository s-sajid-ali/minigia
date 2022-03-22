add_executable(diagnostics_mean diagnostics_mean.cc)
target_link_libraries(diagnostics_mean
    minigia_bunch 
    minigia_lattice
    minigia_simulation
    minigia_command_line
    )
