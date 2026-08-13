#include "wave_tank_options.h"
#include "wavein/airy_wave.h"
#include "wavein/dmstag_hdf5_writer.h"
#include "wavein/projection.h"
#include "wavein/wave_tank.h"
#include "wavein/wavemaker.h"

#include <petscdmstag.h>
#include <petscviewerhdf5.h>

#include <limits>

const char help[] = "Simulate a periodic wave tank to check the algorithm convergence rate.\n";

namespace
{

struct WaveTankOptions
{
        wavein::app::AiryWaveOptions wave;
        wavein::app::WaveTankGridOptions grid;
        PetscReal dt = 0.0;
        char output[PETSC_MAX_PATH_LEN] = "output.h5";
};

PetscErrorCode read_wave_tank_options(WaveTankOptions *options, PetscBool *should_exit)
{
    PetscFunctionBeginUser;

    PetscBool dt_set = PETSC_FALSE;

    PetscOptionsBegin(PETSC_COMM_WORLD, nullptr, "Periodic wave tank options", nullptr);
    PetscCall(
        PetscOptionsReal("-dt", "Time-step size", nullptr, options->dt, &options->dt, &dt_set));
    PetscCall(PetscOptionsString("-output", "HDF5 output filename", nullptr, options->output,
                                 options->output, sizeof(options->output), nullptr));
    PetscOptionsEnd();

    PetscCall(wavein::app::read_airy_wave_options(PETSC_COMM_WORLD, options->wave));
    PetscCall(wavein::app::read_wave_tank_grid_options(PETSC_COMM_WORLD, options->grid));

    PetscCall(PetscOptionsHasHelp(nullptr, should_exit));
    if (*should_exit)
    {
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscCheck(dt_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "Missing required option -dt");
    PetscCheck(options->dt > 0.0 && !PetscIsInfOrNanReal(options->dt), PETSC_COMM_WORLD,
               PETSC_ERR_USER_INPUT, "Time-step size -dt must be positive");

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode initialize_velocity_mask(DM dm, Vec mask)
{
    PetscFunctionBeginUser;

    Vec mask_local = nullptr;
    PetscCall(DMGetLocalVector(dm, &mask_local));
    PetscCall(VecSet(mask_local, 1.0));

    PetscInt ip;
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));

    PetscReal ***c_arr_mask_local = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagVecGetArray(dm, mask_local, &c_arr_mask_local));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscInt startx, starty, nx, ny, ex, ey;
    PetscCall(DMStagGetCorners(dm, &startx, &starty, nullptr, &nx, &ny, nullptr, nullptr, nullptr,
                               nullptr));

    for (ey = starty; ey != starty + ny; ++ey)
    {
        for (ex = startx; ex != startx + nx; ++ex)
        {
            c_arr_mask_local[ey][ex][ip] = 0.0;
        }
    }

    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagVecRestoreArray(dm, mask_local, &c_arr_mask_local));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMLocalToGlobal(dm, mask_local, INSERT_VALUES, mask));
    PetscCall(DMRestoreLocalVector(dm, &mask_local));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode initialize_velocity_and_surface_elevation(DM dm, const wavein::AiryWave &wave,
                                                         Vec sol, Vec eta, PetscReal t = 0.0)
{
    PetscFunctionBeginUser;

    Vec sol_local = nullptr;
    Vec eta_local = nullptr;
    PetscCall(DMCreateLocalVector(dm, &sol_local));
    PetscCall(DMCreateLocalVector(dm, &eta_local));
    PetscCall(VecZeroEntries(sol_local));
    PetscCall(VecZeroEntries(eta_local));

    PetscInt iu, iv;
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iu));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iv));

    PetscInt iprev, icenter;
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));

    PetscReal **c_arr_x = nullptr;
    PetscReal **c_arr_z = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagGetProductCoordinateArraysRead(dm, &c_arr_x, &c_arr_z, nullptr));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscReal ***c_arr_sol_local = nullptr;
    PetscReal ***c_arr_eta_local = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagVecGetArray(dm, sol_local, &c_arr_sol_local));
    PetscCall(DMStagVecGetArray(dm, eta_local, &c_arr_eta_local));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscInt ndim[2];
    PetscInt startx, starty, nx, ny, n_extra[2], ex, ey;
    PetscCall(DMStagGetGlobalSizes(dm, &ndim[0], &ndim[1], nullptr));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, nullptr, &nx, &ny, nullptr, &n_extra[0],
                               &n_extra[1], nullptr));

    for (ey = starty; ey != starty + ny; ++ey)
    {
        const PetscReal z = c_arr_z[ey][icenter];

        for (ex = startx; ex != startx + nx + n_extra[0]; ++ex)
        {
            const PetscReal x = c_arr_x[ex][iprev];
            c_arr_sol_local[ey][ex][iu] = wave.horizontal_velocity(x, z, t);
        }
    }

    for (ey = starty; ey != starty + ny + n_extra[1]; ++ey)
    {
        const PetscReal z = c_arr_z[ey][iprev];

        for (ex = startx; ex != startx + nx; ++ex)
        {
            const PetscReal x = c_arr_x[ex][icenter];
            c_arr_sol_local[ey][ex][iv] = wave.vertical_velocity(x, z, t);

            if (ey == ndim[1])
            {
                c_arr_eta_local[ey][ex][iv] = wave.surface_elevation(x, t);
            }
        }
    }
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagVecRestoreArray(dm, sol_local, &c_arr_sol_local));
    PetscCall(DMStagVecRestoreArray(dm, eta_local, &c_arr_eta_local));
    PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &c_arr_x, &c_arr_z, nullptr));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscCall(DMLocalToGlobal(dm, sol_local, INSERT_VALUES, sol));
    PetscCall(DMLocalToGlobal(dm, eta_local, INSERT_VALUES, eta));

    PetscCall(VecDestroy(&sol_local));
    PetscCall(VecDestroy(&eta_local));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode run()
{
    PetscFunctionBeginUser;

    WaveTankOptions options;
    PetscBool should_exit = PETSC_FALSE;
    PetscCall(read_wave_tank_options(&options, &should_exit));
    if (should_exit)
    {
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    MPI_Comm comm = PETSC_COMM_WORLD;
    wavein::AiryWave wave(comm, options.wave.h, options.wave.T, options.wave.H);

    const PetscReal L = wave.wavelength();
    DM dm = nullptr;
    const PetscInt dof0 = 0, dof1 = 1, dof2 = 1;
    const PetscInt stencil_width = 1;
    PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE,
                             options.grid.Nx, options.grid.Nz, PETSC_DECIDE, PETSC_DECIDE, dof0,
                             dof1, dof2, DMSTAG_STENCIL_BOX, stencil_width, nullptr, nullptr, &dm));
    PetscCall(DMSetFromOptions(dm));
    PetscCall(DMSetUp(dm));

    PetscInt ndim[2];
    PetscCall(DMStagGetGlobalSizes(dm, &ndim[0], &ndim[1], nullptr));
    const PetscReal dx = L / static_cast<PetscReal>(ndim[0]);
    const PetscReal dz = options.wave.h / static_cast<PetscReal>(ndim[1]);

    const PetscReal xmin = 0.0, xmax = L;
    const PetscReal zmin = -options.wave.h, zmax = 0.0;
    PetscCall(DMStagSetUniformCoordinatesProduct(dm, xmin, xmax, zmin, zmax, 0.0, 0.0));

    const PetscReal nin = 0.1, nout = 0.1, gamma = 1.0;
    wavein::Wavemaker wavemaker(comm, dm, L, xmin, xmax, nin, nout, gamma);

    wavein::Projection projection(comm, dm, dx, dz, wavein::kSeawaterDensity);

    wavein::WaveTank wave_tank(comm, dm, wave, wavemaker, projection, wavein::kSeawaterDensity);

    Vec sol = nullptr, eta = nullptr, initial_sol = nullptr, initial_eta = nullptr;
    PetscCall(DMCreateGlobalVector(dm, &sol));
    PetscCall(DMCreateGlobalVector(dm, &eta));
    PetscCall(VecDuplicate(sol, &initial_sol));
    PetscCall(VecDuplicate(eta, &initial_eta));

    Vec velocity_mask = nullptr;
    PetscCall(DMCreateGlobalVector(dm, &velocity_mask));
    PetscCall(initialize_velocity_mask(dm, velocity_mask));

    PetscReal tt = 0.0; // start time
    PetscCall(initialize_velocity_and_surface_elevation(dm, wave, sol, eta, tt));
    PetscCall(VecCopy(sol, initial_sol));
    PetscCall(VecCopy(eta, initial_eta));

    // Create the viewer
    PetscViewer hdf5viewer = nullptr;
    PetscCall(PetscViewerHDF5Open(comm, options.output, FILE_MODE_WRITE, &hdf5viewer));
    wavein::DMStagHDF5Writer hdf5_writer(comm, dm, hdf5viewer);
    PetscCall(hdf5_writer.write_airy_wave(wave));
    PetscCall(hdf5_writer.write_domain(xmin, xmax, xmin, xmax));

    PetscCall(hdf5_writer.push_group());

    // Write the initial conditions to the HDF5 file
    PetscCall(hdf5_writer.write_solution(sol, tt, "velocities and pressure"));
    PetscCall(hdf5_writer.write_surface_elevation(eta, tt, "surface elevation"));

    const PetscReal duration = options.wave.T;
    const PetscReal requested_num_steps = PetscCeilReal(duration / options.dt);
    PetscCheck(requested_num_steps <= static_cast<PetscReal>(std::numeric_limits<PetscInt>::max()),
               comm, PETSC_ERR_USER_INPUT, "Time-step size -dt is too small");
    const auto num_steps = static_cast<PetscInt>(requested_num_steps);
    const PetscReal dt = duration / static_cast<PetscReal>(num_steps);

    for (PetscInt step = 0; step != num_steps; ++step)
    {
        PetscCall(wave_tank.update(sol, eta, tt, dt, 1.0));
        tt += dt;
    }

    // Write the final conditions to the HDF5 file
    PetscCall(hdf5_writer.write_solution(sol, tt));
    PetscCall(hdf5_writer.write_surface_elevation(eta, tt));

    PetscCall(hdf5_writer.pop_group());

    PetscCall(hdf5_writer.write_simulation(0.0, duration, dt, num_steps));

    // After the simulation, eta and velocities should recover the initial conditions
    Vec velocity_error = nullptr, surface_elevation_error = nullptr;
    PetscCall(VecDuplicate(sol, &velocity_error));
    PetscCall(VecDuplicate(eta, &surface_elevation_error));
    PetscCall(VecWAXPY(velocity_error, -1.0, initial_sol, sol));
    PetscCall(VecPointwiseMult(velocity_error, velocity_error, velocity_mask));
    PetscCall(VecWAXPY(surface_elevation_error, -1.0, initial_eta, eta));

    PetscReal velocity_relative_l2_error = 0.0, surface_elevation_relative_l2_error = 0.0;
    PetscCall(VecNorm(velocity_error, NORM_2, &velocity_relative_l2_error));
    PetscCall(VecNorm(surface_elevation_error, NORM_2, &surface_elevation_relative_l2_error));

    PetscReal initial_velocity_norm = 0.0, initial_surface_elevation_norm = 0.0;
    PetscCall(VecNorm(initial_sol, NORM_2, &initial_velocity_norm));
    PetscCall(VecNorm(initial_eta, NORM_2, &initial_surface_elevation_norm));

    velocity_relative_l2_error /= initial_velocity_norm;
    surface_elevation_relative_l2_error /= initial_surface_elevation_norm;

    PetscCall(PetscPrintf(
        comm,
        "nx=%" PetscInt_FMT " nz=%" PetscInt_FMT " steps=%" PetscInt_FMT
        " dt=%.16e velocity_relative_l2_error=%.16e surface_elevation_relative_l2_error=%.16e\n",
        ndim[0], ndim[1], num_steps, static_cast<double>(dt),
        static_cast<double>(velocity_relative_l2_error),
        static_cast<double>(surface_elevation_relative_l2_error)));

    PetscCall(VecDestroy(&velocity_error));
    PetscCall(VecDestroy(&surface_elevation_error));
    PetscCall(VecDestroy(&initial_sol));
    PetscCall(VecDestroy(&initial_eta));
    PetscCall(VecDestroy(&sol));
    PetscCall(VecDestroy(&eta));
    PetscCall(VecDestroy(&velocity_mask));
    PetscCall(DMDestroy(&dm));
    PetscCall(PetscViewerDestroy(&hdf5viewer));

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace

int main(int argc, char **argv)
{
    const PetscErrorCode initialize_error = PetscInitialize(&argc, &argv, nullptr, help);
    if (initialize_error != PETSC_SUCCESS)
    {
        return static_cast<int>(initialize_error);
    }

    const PetscErrorCode run_error = run();
    const PetscErrorCode finalize_error = PetscFinalize();

    return static_cast<int>(run_error != PETSC_SUCCESS ? run_error : finalize_error);
}

// Example usages (run from the repository root after building the release preset):
//
// Define the executable path:
//   app=./build/release/src/apps/periodic_wave_tank
//
// Show the available PETSc options:
//   "${app}" -help
//
// Run one wave period on a single MPI rank:
//   "${app}" -wave_height 1.0 -wave_period 5.0 -water_depth 10.0 \
//       -dt 0.005 -nx 128 -nz 128 \
//       -output output_single.h5
//
// Run the same case on four MPI ranks:
//   mpiexec -n 4 "${app}" -wave_height 1.0 -wave_period 5.0 -water_depth 10.0 \
//       -dt 0.005 -nx 128 -nz 128 \
//       -output output_mpi.h5
//
// Refine the x resolution while holding nz fixed:
//   for nx in 16 32 64 128 256 512; do
//       "${app}" -wave_height 1.0 -wave_period 5.0 -water_depth 10.0 \
//           -dt 0.005 -nx "${nx}" -nz 512 \
//           -output "output_nx_${nx}_nz_512.h5"
//   done
//
// Refine the z resolution while holding nx fixed:
//   for nz in 16 32 64 128 256 512; do
//       "${app}" -wave_height 1.0 -wave_period 5.0 -water_depth 10.0 \
//           -dt 0.005 -nx 512 -nz "${nz}" \
//           -output "output_nx_512_nz_${nz}.h5"
//   done
