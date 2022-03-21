add_library(minigia_collective SHARED
    deposit.cc 
    )

target_link_libraries(minigia_collective 
    minigia_bunch
    minigia_parallel_utils
    )

install(TARGETS minigia_collective DESTINATION lib)
install(FILES
    deposit.hpp
    rectangular_grid_domain.hpp
    rectangular_grid.hpp
    DESTINATION ${INCLUDE_INSTALL_DIR}/minigia/collective)
