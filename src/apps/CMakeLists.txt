add_executable(wavein main.cc)
target_link_libraries(wavein PRIVATE wavein_core)

add_library(wave_tank_app_support wave_tank_options.cc wave_tank_setup.cc)
target_include_directories(wave_tank_app_support PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(wave_tank_app_support PUBLIC wavein_core)

add_executable(periodic_wave_tank periodic_wave_tank.cc)
target_link_libraries(periodic_wave_tank PRIVATE wave_tank_app_support)

add_executable(regular_wave_tank regular_wave_tank.cc)
target_link_libraries(regular_wave_tank PRIVATE wave_tank_app_support)

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
            -sim_start_time 0.0
            -sim_end_time 10.5
            -sim_dt 0.015
            -flow_field_output_start_time 9.0
            -flow_field_output_end_time 10.5
            -flow_surface_elevation_output_start_time 9.0
            -flow_surface_elevation_output_end_time 10.5
            -flow_output_interval 0.075
            -wave_height 0.07
            -wave_period 1.50
            -water_depth 0.70
            -xmin -3.1173251893
            -xmax 28.0559267037
            -nx 400
            -nz 20
            -gamma 6.68
            -ramp_up_time 1.5
            -output regular_wave_tank_smoke.h5
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
