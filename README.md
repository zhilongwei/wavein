# WaveIn: Wave-Vegetation Interaction

This repository is set up for C++ PETSc development in a VS Code dev container.
The container extends the versioned
[`wavein-petsc` debug image](https://github.com/zhilongwei/wavein-petsc-images) and
adds Clang, clang-format, clang-tidy, clangd, LLDB/GDB, and Catch2.

## Project layout

WaveIn uses `.cc` for C++ source files and `.h` for C++ headers:

```text
include/wavein/       Public headers for users of wavein_core
src/                  Library sources and private headers
src/apps/             Application entry points
tests/                Catch2 tests and test-only helpers
```

Only public API headers belong under `include/wavein`; implementation and
application-only headers should remain under `src`.

## Start in VS Code

Requirements:

- Docker with `linux/amd64` container support
- VS Code with the **Dev Containers** extension

Open this directory in VS Code, press `F1`, and run **Dev Containers: Reopen in
Container**. Select the `dev` configure preset when prompted. CMake Tools configures
that preset when the container opens; use its status controls or the commands below
to build and test. The PETSc image is currently amd64-only; Docker Desktop users on
Arm must enable x86/amd64 emulation.

Inside the container, the usual commands are:

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
./build/dev/src/apps/wavein
```

The shared presets do not hardcode a compiler. Inside the dev container, `CXX=g++`
selects GCC/G++. The `dev` preset exports `build/dev/compile_commands.json` for
clangd and runs clang-tidy as part of each build. Microsoft C/C++ IntelliSense is
disabled so clangd exclusively provides parsing, completion, navigation, and
diagnostics. Format-on-save formats source-control modifications when available and
falls back to the entire file otherwise. The same checks are available as CMake
targets:

```bash
cmake --build --preset dev --target format
cmake --build --preset dev --target format-check
cmake --build --preset dev --target tidy
cmake --build --preset dev --target run-wavein
```

Catch2 tests appear in VS Code's Test Explorer through CMake Tools. Use the
**Debug WaveIn** launch configuration to debug the executable with CodeLLDB.

For an optimized WaveIn application build:

```bash
cmake --preset release
cmake --build --preset release
```

This optimizes WaveIn itself but still links the debug PETSc installation in the
development image. Use the release PETSc image when measuring end-to-end performance.

## Build on an HPC system

Transfer or clone the source tree, but do not copy an existing `build` directory
between machines. Load the compiler, MPI, PETSc, CMake, and Ninja modules provided by
the cluster. The MPI C++ wrapper must belong to the same compiler and MPI stack used
to build PETSc; it is often named `mpicxx`, but the exact module names and wrapper
vary by system.

WaveIn accepts `PETSC_DIR` and `PETSC_ARCH` as environment variables or CMake cache
options. `PETSC_ARCH` is only needed for an in-place PETSc build; leave it empty when
`PETSC_DIR` is an installation prefix. A cluster module may instead expose PETSc's
pkg-config file through `PKG_CONFIG_PATH` or `CMAKE_PREFIX_PATH`.

For a one-off release build in a fresh directory, adapt the module names and MPI
wrapper to the cluster:

```bash
module load gcc mpi petsc cmake ninja

cmake --preset release \
    -B build/hpc-release \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DPETSC_DIR="$PETSC_DIR" \
    -DPETSC_ARCH="$PETSC_ARCH"
cmake --build build/hpc-release --parallel
```

For a cluster used regularly, create an untracked `CMakeUserPresets.json`:

```json
{
  "version": 6,
  "configurePresets": [
    {
      "name": "hpc-release",
      "inherits": "release",
      "binaryDir": "${sourceDir}/build/hpc-release",
      "cacheVariables": {
        "CMAKE_CXX_COMPILER": "mpicxx",
        "PETSC_DIR": "$env{PETSC_DIR}",
        "PETSC_ARCH": "$env{PETSC_ARCH}"
      }
    }
  ],
  "buildPresets": [
    {
      "name": "hpc-release",
      "configurePreset": "hpc-release",
      "jobs": 0
    }
  ]
}
```

Then configure and build with:

```bash
cmake --preset hpc-release
cmake --build --preset hpc-release
```

Use a new binary directory whenever the compiler, MPI implementation, or PETSc
installation changes. CMake caches all three as part of the configured toolchain.

Submit the executable through the cluster scheduler, for example
`srun -n 4 ./build/hpc-release/src/apps/wavein`, rather than running a parallel job on a
login node.
