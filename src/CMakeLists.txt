add_library(
    wavein_core
    petsc_info.cc
    airy_wave.cc
    jonswap.cc
    irregular_waves.cc
)

target_include_directories(
    wavein_core
    PUBLIC
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)
target_compile_features(wavein_core PUBLIC cxx_std_17)
target_compile_options(wavein_core PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(wavein_core PUBLIC PkgConfig::PETSC MPI::MPI_CXX)

add_subdirectory(apps)
