#include "wave_tank_options.h"
#include "wave_tank_setup.h"
#include "wavein/airy_wave.h"
#include "wavein/dmstag_hdf5_writer.h"
#include "wavein/projection.h"
#include "wavein/wave_tank.h"
#include "wavein/wavemaker.h"

#include <petscdmstag.h>
#include <petscviewerhdf5.h>

const char help[] =
    "Simulate a regular wave tank using direct inlet generation and source-term outlet "
    "absorption.\n";

namespace
{

struct WaveTankOptions
{
        wavein::app::SimulationOptions simulation;
        wavein::app::WaveTankOutputOptions output;
        wavein::app::AiryWaveOptions wave;
        wavein::app::WaveTankDomainOptions domain;
        wavein::app::WaveTankGridOptions grid;
        wavein::app::WaveTankWavemakerOptions wavemaker;
};

PetscErrorCode read_wave_tank_options(WaveTankOptions *options, PetscBool *should_exit)
{
    PetscFunctionBeginUser;

    PetscCall(wavein::app::read_simulation_options(PETSC_COMM_WORLD, options->simulation));
    PetscCall(wavein::app::read_wave_tank_output_options(PETSC_COMM_WORLD, options->output));
    PetscCall(wavein::app::read_airy_wave_options(PETSC_COMM_WORLD, options->wave));
    PetscCall(wavein::app::read_wave_tank_domain_options(PETSC_COMM_WORLD, options->domain));
    PetscCall(wavein::app::read_wave_tank_grid_options(PETSC_COMM_WORLD, options->grid));
    PetscCall(wavein::app::read_wave_tank_wavemaker_options(PETSC_COMM_WORLD, options->wavemaker));

    PetscCall(PetscOptionsHasHelp(nullptr, should_exit));
    if (*should_exit)
    {
        PetscFunctionReturn(PETSC_SUCCESS);
    }

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
    const wavein::app::WaveTankSchedule schedule =
        wavein::app::create_wave_tank_schedule(options.simulation, options.output);

    const PetscInt num_steps = schedule.num_steps;
    const PetscInt output_stride = schedule.output_stride;
    const PetscInt flow_field_output_start_step = schedule.field_output_start_step;
    const PetscInt flow_field_output_end_step = schedule.field_output_end_step;
    const PetscInt flow_surface_elevation_output_start_step = schedule.elevation_output_start_step;
    const PetscInt flow_surface_elevation_output_end_step = schedule.elevation_output_end_step;

    wavein::AiryWave wave(comm, options.wave.h, options.wave.T, options.wave.H);

    DM dm = nullptr;
    const PetscInt dof0 = 0, dof1 = 1, dof2 = 1;
    const PetscInt stencil_width = 1;
    PetscCall(DMStagCreate2d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, options.grid.Nx,
                             options.grid.Nz, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2,
                             DMSTAG_STENCIL_BOX, stencil_width, nullptr, nullptr, &dm));
    PetscCall(DMSetFromOptions(dm));
    PetscCall(DMSetUp(dm));

    const PetscReal xmin = options.domain.xmin;
    const PetscReal xmax = options.domain.xmax;
    const PetscReal zmin = -options.wave.h;
    const PetscReal zmax = 0.0;
    PetscCall(DMStagSetUniformCoordinatesProduct(dm, xmin, xmax, zmin, zmax, 0.0, 0.0));

    PetscInt ndim[2];
    PetscCall(DMStagGetGlobalSizes(dm, &ndim[0], &ndim[1], nullptr));
    const PetscReal dx = (xmax - xmin) / static_cast<PetscReal>(ndim[0]);
    const PetscReal dz = (zmax - zmin) / static_cast<PetscReal>(ndim[1]);

    const PetscReal forcing_wavelength = wave.wavelength();
    wavein::Wavemaker wavemaker(comm, dm, forcing_wavelength, xmin, xmax, options.wavemaker.nin,
                                options.wavemaker.nout, options.wavemaker.gamma);
    wavein::Projection projection(comm, dm, dx, dz, wavein::kSeawaterDensity);
    wavein::WaveTank wave_tank(comm, dm, wave, wavemaker, projection, wavein::kSeawaterDensity);

    Vec sol = nullptr, eta = nullptr;
    PetscCall(DMCreateGlobalVector(dm, &sol));
    PetscCall(DMCreateGlobalVector(dm, &eta));
    PetscCall(VecZeroEntries(sol));
    PetscCall(VecZeroEntries(eta));

    // Create the HDF5 viewer and writer
    PetscViewer hdf5viewer = nullptr;
    PetscCall(PetscViewerHDF5Open(comm, options.output.output_file, FILE_MODE_WRITE, &hdf5viewer));
    wavein::DMStagHDF5Writer hdf5_writer(comm, dm, hdf5viewer);
    PetscCall(hdf5_writer.write_airy_wave(wave));
    const wavein::ForcingZoneGeometry forcing_zone = wavemaker.forcing_zone_geometry();
    PetscCall(hdf5_writer.write_domain(forcing_zone.xmin, forcing_zone.xmax,
                                       forcing_zone.xmin + forcing_zone.inlet_length,
                                       forcing_zone.xmax - forcing_zone.outlet_length));

    PetscCall(hdf5_writer.write_simulation(options.simulation.start_time,
                                           options.simulation.end_time, options.simulation.dt,
                                           num_steps));

    PetscCall(hdf5_writer.push_group());
    for (PetscInt step = 0; step <= num_steps; ++step)
    {
        const PetscReal time = options.simulation.start_time + step * options.simulation.dt;

        if (step % output_stride == 0 || step == num_steps)
        {
            PetscReal norm;
            PetscCall(VecNorm(eta, NORM_INFINITY, &norm));
            PetscCall(
                PetscPrintf(comm,
                            "step=%" PetscInt_FMT "/%" PetscInt_FMT
                            " time=%.3f s progress=%.1f%% max_eta/a0=%.6e\n",
                            step, num_steps, static_cast<double>(time),
                            100.0 * static_cast<double>(step) / static_cast<double>(num_steps),
                            static_cast<double>(2.0 * norm / wave.wave_height())));
        }

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

        const PetscReal ramp =
            PetscMin(1.0, (time - options.simulation.start_time) / options.wavemaker.ramp_up_time);
        PetscCall(wave_tank.update(sol, eta, time, options.simulation.dt, ramp));
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

// Example usage
// -------------
//
// Prefer scripts/setup_wave_tank.py to validate the physical inputs and write a complete PETSc
// YAML file, then run this app with -options_file_yaml. The explicit options below document the
// derived app interface.
//
// Show the command-line options:
//   app=./build/release/src/apps/regular_wave_tank
//   "${app}" -help
//
// Simulation options:
//   -sim_start_time     Simulation start time in seconds; normally 0.
//   -sim_end_time       Simulation end time in seconds.
//   -sim_dt             Time-step size in seconds.
//
// Output options:
//   -output                                      HDF5 output filename.
//   -flow_field_output_start_time                First output time for u, w, and p.
//   -flow_field_output_end_time                  Last output time for u, w, and p.
//   -flow_surface_elevation_output_start_time    First output time for eta.
//   -flow_surface_elevation_output_end_time      Last output time for eta.
//   -flow_output_interval                        Time between output samples in seconds.
//
// Wave options:
//   -wave_height       Wave height H in metres.
//   -wave_period       Wave period T in seconds.
//   -water_depth       Still-water depth h in metres; the vertical domain is [-h, 0].
//
// Domain and grid options:
//   -xmin              Left computational boundary in metres.
//   -xmax              Right computational boundary in metres.
//   -nx                Total number of cells over the complete computational domain.
//   -nz                Number of cells over the water depth.
//
// Wavemaker options:
//   -nin               Inlet forcing-zone length in wavelengths; default: 1.
//   -nout              Outlet forcing-zone length in wavelengths; default: 2.
//   -gamma             Outlet source-relaxation rate in 1/s.
//   -ramp_up_time      Duration of the linear wavemaker ramp in seconds.
//
// By convention, the domain of interest starts at x=0 and ends at x=Ltank. Let bin and bout be
// the buffer lengths separating it from the inlet and outlet forcing zones, in wavelengths. Then
//   xmin = -(nin+bin)*L
//   xmax = Ltank+(bout+nout)*L
// If Ltank=NL*L and nxl cells per wavelength are required, use
// nx=(NL+nin+bin+bout+nout)*nxl. The inlet forcing zone ends at -bin*L, and the outlet forcing
// zone starts at Ltank+bout*L.
//
// For T=1.5 s, h=0.7 m, and H=0.07 m, the finite-depth Airy wavelength is
// L=3.1173251893 m. A seven-wavelength domain of interest with nin=1, nout=2, bin=bout=1, and
// 20 cells per wavelength therefore uses:
//   Ltank = 7*L = 21.8212763248 m
//   xmin = -2*L = -6.2346503785 m
//   xmax = Ltank+3*L = 10*L = 31.1732518926 m
//   nx = (7+1+1+1+2)*20 = 240
//   nz = 10
//
// Simulation and output times are conveniently defined using the wave period. The target case
// uses a ramp time of 2*T, a simulation duration of 7*T, dt=T/100, and an output interval of
// T/20. The Airy phase celerity is c=L/T, so a wave travels across the seven-wavelength domain of
// interest in 7*T=10.5 s. The solution is recorded over [6*T, 7*T], or [9.0, 10.5] s. The
// lighter surface-elevation output is recorded over [2*T, 7*T], or [3.0, 10.5] s.
//
// Generate the target case. The requested 21 m domain of interest is rounded up to 7*L. The
// simulation lasts 7*T, the solution is recorded for the final 1*T, the surface elevation for
// the final 5*T, and the recording rate is 20 frames/T:
//   uv run --project scripts python scripts/setup_wave_tank.py regular \
//       --water-depth 0.7 --wave-period 1.5 --wave-height 0.07 \
//       --tank-length 21.0 --cells-per-wavelength 20 --cells-per-depth 10 \
//       --inlet-buffer-wavelengths 1 --outlet-buffer-wavelengths 1 \
//       --cfl 0.5 --simulation-periods 7 --solution-record-periods 1 \
//       --surface-elevation-record-periods 5 --frames-per-period 20 \
//       --outlet-relaxation-strength 10.02 \
//       --input-file regular_wave_tank.yaml --output regular_wave_tank.h5
//
// Run the generated case on one MPI rank:
//   "${app}" -options_file_yaml regular_wave_tank.yaml
//
// Run the same case on four MPI ranks:
//   mpiexec -n 4 "${app}" -options_file_yaml regular_wave_tank.yaml \
//       -output regular_wave_tank_mpi.h5
