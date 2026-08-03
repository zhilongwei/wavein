cmake_minimum_required(VERSION 3.25)

project(WaveIn VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)

set(PETSC_DIR "$ENV{PETSC_DIR}" CACHE PATH "PETSc source or installation directory")
set(PETSC_ARCH "$ENV{PETSC_ARCH}" CACHE STRING "PETSc in-place build architecture")

if(PETSC_DIR)
    if(PETSC_ARCH)
        set(_petsc_prefix "${PETSC_DIR}/${PETSC_ARCH}")
    else()
        set(_petsc_prefix "${PETSC_DIR}")
    endif()

    list(PREPEND CMAKE_PREFIX_PATH "${_petsc_prefix}")
endif()

find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc petsc)

option(WAVEIN_ENABLE_CLANG_TIDY "Run clang-tidy while compiling" OFF)

if(WAVEIN_ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXECUTABLE NAMES clang-tidy REQUIRED)
    set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
endif()

add_library(wavein_core src/petsc_info.cc)
target_include_directories(
    wavein_core
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_compile_features(wavein_core PUBLIC cxx_std_17)
target_compile_options(wavein_core PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(wavein_core PUBLIC PkgConfig::PETSC)

add_executable(wavein src/apps/main.cc)
target_link_libraries(wavein PRIVATE wavein_core)

add_custom_target(
    run-wavein
    COMMAND "$<TARGET_FILE:wavein>"
    DEPENDS wavein
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    USES_TERMINAL
    COMMENT "Running the WaveIn application"
    VERBATIM
)

if(BUILD_TESTING)
    find_package(Catch2 3 CONFIG REQUIRED)
    list(APPEND CMAKE_MODULE_PATH "${Catch2_DIR}")
    include(Catch)

    add_executable(wavein_tests tests/petsc_info_test.cc)
    target_link_libraries(wavein_tests PRIVATE wavein_core Catch2::Catch2)
    catch_discover_tests(wavein_tests)
endif()

find_program(CLANG_FORMAT_EXECUTABLE NAMES clang-format)
if(CLANG_FORMAT_EXECUTABLE)
    file(
        GLOB_RECURSE
        WAVEIN_FORMAT_SOURCES
        CONFIGURE_DEPENDS
        LIST_DIRECTORIES false
        "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hxx"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hh"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hxx"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cc"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cxx"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.hh"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.hxx"
    )
    list(SORT WAVEIN_FORMAT_SOURCES)

    add_custom_target(
        format
        COMMAND "${CLANG_FORMAT_EXECUTABLE}" -i ${WAVEIN_FORMAT_SOURCES}
        COMMENT "Formatting WaveIn C++ sources"
        VERBATIM
    )
    add_custom_target(
        format-check
        COMMAND "${CLANG_FORMAT_EXECUTABLE}" --dry-run --Werror ${WAVEIN_FORMAT_SOURCES}
        COMMENT "Checking WaveIn C++ formatting"
        VERBATIM
    )
else()
    message(STATUS "clang-format not found; format targets are unavailable")
endif()

find_program(RUN_CLANG_TIDY_EXECUTABLE NAMES run-clang-tidy run-clang-tidy-18)
if(RUN_CLANG_TIDY_EXECUTABLE)
    add_custom_target(
        tidy
        COMMAND "${RUN_CLANG_TIDY_EXECUTABLE}" -quiet -p "${CMAKE_BINARY_DIR}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMENT "Running clang-tidy over the compilation database"
        VERBATIM
    )
else()
    message(STATUS "run-clang-tidy not found; tidy target is unavailable")
endif()
