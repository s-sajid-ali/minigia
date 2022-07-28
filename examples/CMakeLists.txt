add_executable(test_deposit test_deposit.cc)
target_link_libraries(test_deposit
    minigia_collective
    )

add_executable(test_point_charge test_point_charge.cc)
target_link_libraries(test_point_charge
    minigia_collective
    minigia_simulation
    )
