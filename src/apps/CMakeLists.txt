add_executable(wavein main.cc)
target_link_libraries(wavein PRIVATE wavein_core)

add_executable(periodic_wave_tank periodic_wave_tank.cc)
target_link_libraries(periodic_wave_tank PRIVATE wavein_core)

if(BUILD_TESTING)
    add_test(
        NAME periodic_wave_tank_smoke
        COMMAND
            periodic_wave_tank
            -wave_height 1.0
            -wave_period 5.0
            -water_depth 10.0
            -dt 0.1
            -nx 8
            -nz 4
    )
endif()

add_custom_target(
    run-wavein
    COMMAND "$<TARGET_FILE:wavein>"
    DEPENDS wavein
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    USES_TERMINAL
    COMMENT "Running the WaveIn application"
    VERBATIM
)
