#include "wavein/regular_wave_tank.h"
#include "wavein/airy_wave.h"
#include "wavein/dmstag_hdf5_writer.h"
#include "wavein/projection.h"
#include "wavein/wavemaker.h"

#include <petscdmstag.h>
#include <petscviewerhdf5.h>

const char help[] = "Simulate a regular wave tank using the forcing zone method.\n";

namespace
{

struct WaveTankOptions
{
        // Simulation control parameters
        PetscReal sim_start_time = 0.0;
        PetscReal sim_end_time = 0.0;
        PetscReal sim_dt = 0.0;

        // Output parameters
        char output[PETSC_MAX_PATH_LEN] = "output.h5";
        PetscReal flow_field_output_start_time = 0.0;
        PetscReal flow_field_output_end_time = 0.0;
        PetscReal flow_surface_elevation_output_start_time = 0.0;
        PetscReal flow_surface_elevation_output_end_time = 0.0;
        PetscReal flow_output_interval = 0.0;

        // Wave parameters
        PetscReal H = 0.0;
        PetscReal T = 0.0;
        PetscReal h = 0.0;

        // Simulation domain parameters
        PetscReal xmin = 0.0;
        PetscReal xmax = 0.0;

        // Grid parameters
        PetscInt Nx = 0;
        PetscInt Nz = 0;

        // Wavemaker parameters
        PetscReal nin = 0.0;
        PetscReal nout = 0.0;
        PetscReal gamma = 0.0;
        PetscReal ramp_up_time = 0.0;
};

PetscErrorCode interval_to_steps(PetscReal interval, PetscReal dt, const char *option_name,
                                 PetscInt *steps)
{
    PetscFunctionBeginUser;

    const PetscReal ratio = interval / dt;
    const PetscReal nearest = PetscFloorReal(ratio + 0.5);
    const PetscReal tolerance = 100.0 * PETSC_MACHINE_EPSILON * PetscMax(1.0, PetscAbsReal(ratio));

    PetscCheck(PetscAbsReal(ratio - nearest) <= tolerance, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "%s must be a multiple of -sim_dt", option_name);

    *steps = static_cast<PetscInt>(nearest);

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_wave_tank_options(WaveTankOptions *options, PetscBool *should_exit)
{
    PetscFunctionBeginUser;

    PetscBool sim_start_time_set = PETSC_FALSE;
    PetscBool sim_end_time_set = PETSC_FALSE;
    PetscBool sim_dt_set = PETSC_FALSE;
    PetscBool flow_field_output_start_time_set = PETSC_FALSE;
    PetscBool flow_field_output_end_time_set = PETSC_FALSE;
    PetscBool flow_surface_elevation_output_start_time_set = PETSC_FALSE;
    PetscBool flow_surface_elevation_output_end_time_set = PETSC_FALSE;
    PetscBool flow_output_interval_set = PETSC_FALSE;
    PetscBool H_set = PETSC_FALSE;
    PetscBool T_set = PETSC_FALSE;
    PetscBool h_set = PETSC_FALSE;
    PetscBool xmin_set = PETSC_FALSE;
    PetscBool xmax_set = PETSC_FALSE;
    PetscBool Nx_set = PETSC_FALSE;
    PetscBool Nz_set = PETSC_FALSE;
    PetscBool nin_set = PETSC_FALSE;
    PetscBool nout_set = PETSC_FALSE;
    PetscBool gamma_set = PETSC_FALSE;
    PetscBool ramp_up_time_set = PETSC_FALSE;

    PetscOptionsBegin(PETSC_COMM_WORLD, nullptr, "Regular wave tank options", nullptr);
    PetscCall(PetscOptionsReal("-sim_start_time", "Simulation start time", nullptr,
                               options->sim_start_time, &options->sim_start_time,
                               &sim_start_time_set));
    PetscCall(PetscOptionsReal("-sim_end_time", "Simulation end time", nullptr,
                               options->sim_end_time, &options->sim_end_time, &sim_end_time_set));
    PetscCall(PetscOptionsReal("-sim_dt", "Simulation time-step size", nullptr, options->sim_dt,
                               &options->sim_dt, &sim_dt_set));

    PetscCall(PetscOptionsString("-output", "HDF5 output filename", nullptr, options->output,
                                 options->output, sizeof(options->output), nullptr));
    PetscCall(PetscOptionsReal("-flow_field_output_start_time", "Flow-field output start time",
                               nullptr, options->flow_field_output_start_time,
                               &options->flow_field_output_start_time,
                               &flow_field_output_start_time_set));
    PetscCall(PetscOptionsReal("-flow_field_output_end_time", "Flow-field output end time", nullptr,
                               options->flow_field_output_end_time,
                               &options->flow_field_output_end_time,
                               &flow_field_output_end_time_set));
    PetscCall(PetscOptionsReal("-flow_surface_elevation_output_start_time",
                               "Surface-elevation output start time", nullptr,
                               options->flow_surface_elevation_output_start_time,
                               &options->flow_surface_elevation_output_start_time,
                               &flow_surface_elevation_output_start_time_set));
    PetscCall(PetscOptionsReal("-flow_surface_elevation_output_end_time",
                               "Surface-elevation output end time", nullptr,
                               options->flow_surface_elevation_output_end_time,
                               &options->flow_surface_elevation_output_end_time,
                               &flow_surface_elevation_output_end_time_set));
    PetscCall(PetscOptionsReal("-flow_output_interval", "Flow output interval", nullptr,
                               options->flow_output_interval, &options->flow_output_interval,
                               &flow_output_interval_set));

    PetscCall(
        PetscOptionsReal("-wave_height", "Wave height", nullptr, options->H, &options->H, &H_set));
    PetscCall(
        PetscOptionsReal("-wave_period", "Wave period", nullptr, options->T, &options->T, &T_set));
    PetscCall(
        PetscOptionsReal("-water_depth", "Water depth", nullptr, options->h, &options->h, &h_set));

    PetscCall(PetscOptionsReal("-xmin", "Minimum x-coordinate", nullptr, options->xmin,
                               &options->xmin, &xmin_set));
    PetscCall(PetscOptionsReal("-xmax", "Maximum x-coordinate", nullptr, options->xmax,
                               &options->xmax, &xmax_set));
    PetscCall(PetscOptionsInt("-nx", "Number of cells in the x-direction", nullptr, options->Nx,
                              &options->Nx, &Nx_set));
    PetscCall(PetscOptionsInt("-nz", "Number of cells in the z-direction", nullptr, options->Nz,
                              &options->Nz, &Nz_set));

    PetscCall(PetscOptionsReal("-nin", "Inlet forcing-zone length in wavelengths", nullptr,
                               options->nin, &options->nin, &nin_set));
    PetscCall(PetscOptionsReal("-nout", "Outlet forcing-zone length in wavelengths", nullptr,
                               options->nout, &options->nout, &nout_set));
    PetscCall(PetscOptionsReal("-gamma", "Forcing-zone strength", nullptr, options->gamma,
                               &options->gamma, &gamma_set));
    PetscCall(PetscOptionsReal("-ramp_up_time", "Wavemaker ramp-up time", nullptr,
                               options->ramp_up_time, &options->ramp_up_time, &ramp_up_time_set));
    PetscOptionsEnd();

    PetscCall(PetscOptionsHasHelp(nullptr, should_exit));
    if (*should_exit)
    {
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscCheck(sim_start_time_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Missing required option -sim_start_time");
    PetscCheck(sim_end_time_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Missing required option -sim_end_time");
    PetscCheck(sim_dt_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Missing required option -sim_dt");
    PetscCheck(flow_field_output_start_time_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Missing required option -flow_field_output_start_time");
    PetscCheck(flow_field_output_end_time_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Missing required option -flow_field_output_end_time");
    PetscCheck(flow_surface_elevation_output_start_time_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Missing required option -flow_surface_elevation_output_start_time");
    PetscCheck(flow_surface_elevation_output_end_time_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Missing required option -flow_surface_elevation_output_end_time");
    PetscCheck(flow_output_interval_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Missing required option -flow_output_interval");
    PetscCheck(H_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Missing required option -wave_height");
    PetscCheck(T_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Missing required option -wave_period");
    PetscCheck(h_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Missing required option -water_depth");
    PetscCheck(xmin_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "Missing required option -xmin");
    PetscCheck(xmax_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "Missing required option -xmax");
    PetscCheck(Nx_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "Missing required option -nx");
    PetscCheck(Nz_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "Missing required option -nz");
    PetscCheck(nin_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "Missing required option -nin");
    PetscCheck(nout_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "Missing required option -nout");
    PetscCheck(gamma_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "Missing required option -gamma");
    PetscCheck(ramp_up_time_set, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Missing required option -ramp_up_time");

    PetscCheck(options->sim_start_time >= 0.0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Simulation start time -sim_start_time must be nonnegative");
    PetscCheck(options->sim_end_time > options->sim_start_time, PETSC_COMM_WORLD,
               PETSC_ERR_USER_INPUT,
               "Simulation end time -sim_end_time must be later than -sim_start_time");
    PetscCheck(options->sim_dt > 0.0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Simulation time-step size -sim_dt must be positive");

    PetscCheck(options->flow_field_output_start_time >= options->sim_start_time, PETSC_COMM_WORLD,
               PETSC_ERR_USER_INPUT,
               "Flow-field output must not start before the simulation start time");
    PetscCheck(options->flow_field_output_end_time >= options->flow_field_output_start_time,
               PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Flow-field output end time must not precede its start time");
    PetscCheck(options->flow_field_output_end_time <= options->sim_end_time, PETSC_COMM_WORLD,
               PETSC_ERR_USER_INPUT,
               "Flow-field output must not end after the simulation end time");
    PetscCheck(options->flow_surface_elevation_output_start_time >= options->sim_start_time,
               PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Surface-elevation output must not start before the simulation start time");
    PetscCheck(options->flow_surface_elevation_output_end_time >=
                   options->flow_surface_elevation_output_start_time,
               PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Surface-elevation output end time must not precede its start time");
    PetscCheck(options->flow_surface_elevation_output_end_time <= options->sim_end_time,
               PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Surface-elevation output must not end after the simulation end time");
    PetscCheck(options->flow_output_interval > 0.0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Flow output interval -flow_output_interval must be positive");

    PetscCheck(options->H > 0.0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Wave height -wave_height must be positive");
    PetscCheck(options->T > 0.0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Wave period -wave_period must be positive");
    PetscCheck(options->h > 0.0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Water depth -water_depth must be positive");
    PetscCheck(options->xmax > options->xmin, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Maximum x-coordinate -xmax must be greater than -xmin");
    PetscCheck(options->Nx >= 2, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Cell count -nx must be at least 2");
    PetscCheck(options->Nz >= 2, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Cell count -nz must be at least 2");
    PetscCheck(options->nin > 0.0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Inlet forcing-zone length -nin must be positive");
    PetscCheck(options->nout > 0.0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Outlet forcing-zone length -nout must be positive");
    PetscCheck(options->gamma > 0.0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Forcing-zone strength -gamma must be positive");
    PetscCheck(options->ramp_up_time > 0.0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Wavemaker ramp-up time -ramp_up_time must be positive");

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

    PetscInt num_steps;
    PetscInt output_stride;
    PetscInt flow_field_output_start_step, flow_field_output_end_step;
    PetscInt flow_surface_elevation_output_start_step, flow_surface_elevation_output_end_step;

    PetscCall(interval_to_steps(options.sim_end_time - options.sim_start_time, options.sim_dt,
                                "-sim_dt", &num_steps));
    PetscCall(interval_to_steps(options.flow_output_interval, options.sim_dt,
                                "-flow_output_interval", &output_stride));
    PetscCheck(output_stride >= 1, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "-flow_output_interval must be at least -sim_dt");
    PetscCall(interval_to_steps(options.flow_field_output_start_time - options.sim_start_time,
                                options.sim_dt, "-flow_field_output_start_time",
                                &flow_field_output_start_step));
    PetscCall(interval_to_steps(options.flow_field_output_end_time - options.sim_start_time,
                                options.sim_dt, "-flow_field_output_end_time",
                                &flow_field_output_end_step));
    PetscCall(interval_to_steps(
        options.flow_surface_elevation_output_start_time - options.sim_start_time, options.sim_dt,
        "-flow_surface_elevation_output_start_time", &flow_surface_elevation_output_start_step));
    PetscCall(interval_to_steps(
        options.flow_surface_elevation_output_end_time - options.sim_start_time, options.sim_dt,
        "-flow_surface_elevation_output_end_time", &flow_surface_elevation_output_end_step));

    PetscCheck((flow_field_output_end_step - flow_field_output_start_step) % output_stride == 0,
               PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Flow-field output window must contain an integer number of output intervals");

    PetscCheck(
        (flow_surface_elevation_output_end_step - flow_surface_elevation_output_start_step) %
                output_stride ==
            0,
        PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
        "Surface-elevation output window must contain an integer number of output intervals");

    MPI_Comm comm = PETSC_COMM_WORLD;
    wavein::AiryWave wave(comm, options.h, options.T, options.H);

    DM dm = nullptr;
    const PetscInt dof0 = 0, dof1 = 1, dof2 = 1;
    const PetscInt stencil_width = 1;
    PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, options.Nx,
                             options.Nz, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2,
                             DMSTAG_STENCIL_BOX, stencil_width, nullptr, nullptr, &dm));
    PetscCall(DMSetFromOptions(dm));
    PetscCall(DMSetUp(dm));

    const PetscReal xmin = options.xmin, xmax = options.xmax;
    const PetscReal zmin = -options.h, zmax = 0.0;
    PetscCall(DMStagSetUniformCoordinatesProduct(dm, xmin, xmax, zmin, zmax, 0.0, 0.0));

    PetscInt ndim[2];
    PetscCall(DMStagGetGlobalSizes(dm, &ndim[0], &ndim[1], nullptr));
    const PetscReal dx = (xmax - xmin) / static_cast<PetscReal>(ndim[0]);
    const PetscReal dz = (zmax - zmin) / static_cast<PetscReal>(ndim[1]);

    wavein::Wavemaker wavemaker(comm, dm, wave.wavelength(), xmin, xmax, options.nin, options.nout,
                                options.gamma);
    wavein::Projection projection(comm, dm, dx, dz, wavein::kSeawaterDensity);
    wavein::RegularWaveTank wave_tank(comm, dm, wave, wavemaker, projection,
                                      wavein::kSeawaterDensity);

    Vec sol = nullptr, eta = nullptr;
    PetscCall(DMCreateGlobalVector(dm, &sol));
    PetscCall(DMCreateGlobalVector(dm, &eta));
    PetscCall(VecSet(sol, 0.0));
    PetscCall(VecSet(eta, 0.0));

    // Create the HDF5 viewer and writer
    PetscViewer hdf5viewer = nullptr;
    PetscCall(PetscViewerHDF5Open(comm, options.output, FILE_MODE_WRITE, &hdf5viewer));
    wavein::DMStagHDF5Writer hdf5_writer(comm, dm, hdf5viewer);

    PetscCall(hdf5_writer.push_group());
    PetscCall(hdf5_writer.write(wave));

    for (PetscInt step = 0; step <= num_steps; ++step)
    {
        const PetscReal time = options.sim_start_time + step * options.sim_dt;

        if (step >= flow_field_output_start_step && step <= flow_field_output_end_step &&
            (step - flow_field_output_start_step) % output_stride == 0)
        {
            PetscCall(hdf5_writer.write_solution(sol, time));
        }

        if (step >= flow_surface_elevation_output_start_step &&
            step <= flow_surface_elevation_output_end_step &&
            (step - flow_surface_elevation_output_start_step) % output_stride == 0)
        {
            PetscCall(hdf5_writer.write_surface_elevation(eta, time));
        }

        if (step == num_steps)
        {
            break;
        }

        const PetscReal factor =
            PetscMin(1.0, (time - options.sim_start_time) / options.ramp_up_time);
        PetscCall(wave_tank.update(sol, eta, time, options.sim_dt, factor));

        if (step % output_stride == 0 || step == num_steps)
        {
            PetscCall(PetscPrintf(
                comm, "step=%" PetscInt_FMT "/%" PetscInt_FMT " time=%.3f s progress=%.1f%%\n",
                step, num_steps, static_cast<double>(time),
                100.0 * static_cast<double>(step) / static_cast<double>(num_steps)));
        }
    }

    PetscCall(hdf5_writer.pop_group());

    PetscCall(VecDestroy(&sol));
    PetscCall(VecDestroy(&eta));
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
