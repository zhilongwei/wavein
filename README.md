# WaveIn: Wave-Vegetation Interaction

WaveIn is a C++ PETSc project for wave-vegetation simulations. Its VS Code dev
container extends the versioned PETSc debug image from
[`wavein-petsc-images`](https://github.com/zhilongwei/wavein-petsc-images) with
the compiler, debugger, formatting, linting, and Catch2 tools used for development.

## Layout

```text
include/wavein/       Public wavein_core headers
src/                  Library sources and private headers
src/apps/             Applications
tests/                Catch2 tests
scripts/              Python utilities and locked environment
```

## Development

Install Docker with `linux/amd64` support and VS Code's **Dev Containers**
extension, then run **Dev Containers: Reopen in Container**. Inside the container:

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
./build/dev/src/apps/wavein
```

Formatting and static-analysis targets are also available:

```bash
cmake --build --preset dev --target format-check
cmake --build --preset dev --target tidy
```

For an optimized WaveIn build:

```bash
cmake --preset release
cmake --build --preset release
```

The development container still links the debug PETSc image; use the release
PETSc image for end-to-end performance measurements.

## HPC builds

Load matching compiler, MPI, PETSc, CMake, and Ninja modules. Do not copy a
configured `build` directory between systems. For example:

```bash
cmake --preset release \
    -B build/hpc-release \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DPETSC_DIR="$PETSC_DIR" \
    -DPETSC_ARCH="$PETSC_ARCH"
cmake --build build/hpc-release --parallel
```

`PETSC_ARCH` is needed only for an in-place PETSc build. Use a fresh build
directory whenever the compiler, MPI implementation, or PETSc installation
changes, and run parallel applications through the cluster scheduler.

## Continuous integration

Pushes and pull requests test against the PETSc debug image. Image publication,
weekly compatibility checks, and manual runs test both debug and release images.
The Python workflow validates the locked environment on Python 3.10, and
Dependabot maintains Python and GitHub Actions dependencies.
