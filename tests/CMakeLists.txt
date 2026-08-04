find_package(Catch2 3 CONFIG REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${Catch2_DIR}")
include(Catch)

add_executable(
    wavein_tests
    petsc_info_test.cc
    airy_wave_test.cc
    wave_spectrum_test.cc
    projection_test.cc
)
target_link_libraries(wavein_tests PRIVATE wavein_core Catch2::Catch2)
catch_discover_tests(wavein_tests)
