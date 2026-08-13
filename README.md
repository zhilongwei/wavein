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

## Wave-tank input files

Use `scripts/setup_wave_tank.py` to validate the primary dimensional inputs and write a PETSc YAML
input file. Simulation and recording durations are specified in wave periods; resolution,
forcing-zone lengths, and the buffers around the domain of interest are specified relative to the
wavelength or water depth; and the recording rate is specified in frames per period. The requested
domain length of interest remains dimensional and is rounded upward to a whole number of
wavelengths. The final recording durations are set independently for the complete solution and the
lighter surface-elevation output.

For example:

```bash
uv run --project scripts python scripts/setup_wave_tank.py regular \
    --water-depth 0.7 --wave-period 1.5 --wave-height 0.07 \
    --tank-length 21.0 --cells-per-wavelength 20 --cells-per-depth 10 \
    --inlet-buffer-wavelengths 1 --outlet-buffer-wavelengths 1 \
    --cfl 0.5 --simulation-periods 7 --solution-record-periods 1 \
    --surface-elevation-record-periods 5 \
    --frames-per-period 20 \
    --input-file regular.yaml --output regular_wave_tank.h5

./build/release/src/apps/regular_wave_tank -options_file_yaml regular.yaml
```

Run `setup_wave_tank.py regular -h` or `setup_wave_tank.py irregular -h` for all inputs. The
generated YAML records the derived dimensional quantities as comments and contains every PETSc
option needed by the corresponding application.

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
