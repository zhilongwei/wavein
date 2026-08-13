add_executable(wavein main.cc)
target_link_libraries(wavein PRIVATE wavein_core)

add_library(wave_tank_app_support wave_tank_options.cc wave_tank_setup.cc)
target_include_directories(wave_tank_app_support PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(wave_tank_app_support PUBLIC wavein_core)

add_executable(periodic_wave_tank periodic_wave_tank.cc)
target_link_libraries(periodic_wave_tank PRIVATE wave_tank_app_support)

add_executable(regular_wave_tank regular_wave_tank.cc)
target_link_libraries(regular_wave_tank PRIVATE wave_tank_app_support)

add_executable(irregular_wave_tank irregular_wave_tank.cc)
target_link_libraries(irregular_wave_tank PRIVATE wave_tank_app_support)

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
            -output periodic_wave_tank_smoke.h5
    )

    add_test(
        NAME regular_wave_tank_smoke
        COMMAND
            regular_wave_tank
            -options_file_yaml
            "${PROJECT_SOURCE_DIR}/tests/inputs/regular_wave_tank_smoke.yaml"
    )

    add_test(
        NAME irregular_wave_tank_smoke
        COMMAND
            irregular_wave_tank
            -options_file_yaml
            "${PROJECT_SOURCE_DIR}/tests/inputs/irregular_wave_tank_smoke.yaml"
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
