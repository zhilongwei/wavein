#include "wave_tank_options.h"
#include "wave_tank_setup.h"
#include "wavein/dmstag_hdf5_writer.h"
#include "wavein/irregular_waves.h"
#include "wavein/jonswap.h"
#include "wavein/projection.h"
#include "wavein/wave_tank.h"
#include "wavein/wavemaker.h"

#include <petscdmstag.h>
#include <petscviewerhdf5.h>

const char help[] =
    "Simulate an irregular wave tank using direct inlet generation and source-term outlet "
    "absorption.\n";

namespace
{

struct WaveTankOptions
{
        wavein::app::SimulationOptions simulation;
        wavein::app::WaveTankOutputOptions output;
        wavein::app::IrregularWaveOptions waves;
        wavein::app::WaveTankDomainOptions domain;
        wavein::app::WaveTankGridOptions grid;
        wavein::app::WaveTankWavemakerOptions wavemaker;
};

PetscErrorCode read_wave_tank_options(WaveTankOptions *options, PetscBool *should_exit)
{
    PetscFunctionBeginUser;

    PetscCall(wavein::app::read_simulation_options(PETSC_COMM_WORLD, options->simulation));
    PetscCall(wavein::app::read_wave_tank_output_options(PETSC_COMM_WORLD, options->output));
    PetscCall(wavein::app::read_irregular_wave_options(PETSC_COMM_WORLD, options->waves));
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

    wavein::Jonswap jonswap(options.waves.Hs, options.waves.Tp, options.waves.gamma);
    wavein::IrregularWaves irregular_waves(comm, jonswap, options.waves.h, options.waves.omega_min,
                                           options.waves.omega_max, options.waves.component_count,
                                           options.waves.random_seed);

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
    const PetscReal zmin = -options.waves.h;
    const PetscReal zmax = 0.0;
    PetscCall(DMStagSetUniformCoordinatesProduct(dm, xmin, xmax, zmin, zmax, 0.0, 0.0));

    PetscInt ndim[2];
    PetscCall(DMStagGetGlobalSizes(dm, &ndim[0], &ndim[1], nullptr));
    const PetscReal dx = (xmax - xmin) / static_cast<PetscReal>(ndim[0]);
    const PetscReal dz = (zmax - zmin) / static_cast<PetscReal>(ndim[1]);

    const PetscReal forcing_wavelength = irregular_waves.wavelength_at_peak_period();
    wavein::Wavemaker wavemaker(comm, dm, forcing_wavelength, xmin, xmax, options.wavemaker.nin,
                                options.wavemaker.nout, options.wavemaker.gamma);
    wavein::Projection projection(comm, dm, dx, dz, wavein::kSeawaterDensity);
    wavein::WaveTank wave_tank(comm, dm, irregular_waves, wavemaker, projection,
                               wavein::kSeawaterDensity);

    Vec sol = nullptr, eta = nullptr;
    PetscCall(DMCreateGlobalVector(dm, &sol));
    PetscCall(DMCreateGlobalVector(dm, &eta));
    PetscCall(VecZeroEntries(sol));
    PetscCall(VecZeroEntries(eta));

    // Create the HDF5 viewer and writer
    PetscViewer hdf5viewer = nullptr;
    PetscCall(PetscViewerHDF5Open(comm, options.output.output_file, FILE_MODE_WRITE, &hdf5viewer));
    wavein::DMStagHDF5Writer hdf5_writer(comm, dm, hdf5viewer);
    const wavein::ForcingZoneGeometry forcing_zone = wavemaker.forcing_zone_geometry();
    PetscCall(hdf5_writer.write_irregular_waves(irregular_waves));
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
                            " time=%.3f s progress=%.1f%% max_eta=%.6e\n",
                            step, num_steps, static_cast<double>(time),
                            100.0 * static_cast<double>(step) / static_cast<double>(num_steps),
                            static_cast<double>(norm)));
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
//   app=./build/release/src/apps/irregular_wave_tank
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
// Irregular-wave options:
//   -significant_wave_height     Significant wave height Hs in metres.
//   -peak_period                 Peak period Tp in seconds.
//   -peak_enhancement_factor     JONSWAP peak-enhancement factor; default: 3.3.
//   -water_depth                 Still-water depth h in metres; the vertical domain is [-h, 0].
//   -omega_min                   Lower angular-frequency bound in rad/s.
//   -omega_max                   Upper angular-frequency bound in rad/s.
//   -component_count             Number of uniformly spaced wave components.
//   -random_seed                 Seed for the random component phases; default: 0.
//
// Domain and grid options:
//   -xmin              Left computational boundary in metres.
//   -xmax              Right computational boundary in metres.
//   -nx                Total number of cells over the complete computational domain.
//   -nz                Number of cells over the water depth.
//
// Wavemaker options:
//   -nin               Inlet forcing-zone length in peak wavelengths; default: 1.
//   -nout              Outlet forcing-zone length in peak wavelengths; default: 2.
//   -gamma             Outlet source-relaxation rate in 1/s.
//   -ramp_up_time      Duration of the linear wavemaker ramp in seconds.
//
// The example below uses h=10 m, Tp=15.24 s, Hs=1 m, and 100 JONSWAP components. The
// finite-depth wavelength and phase celerity at the peak period are
//   Lp = 146.5489406424 m
//   cp = Lp/Tp = 9.6160722206 m/s
// The frequency interval uses Tz=0.834*Tp and [omega_min,omega_max]=[1/Tz,20/Tz].
//
// A requested domain length of interest of 700 m is rounded upward to 5*Lp=732.7447032119 m. With
// a one-wavelength buffer on each side, a one-wavelength inlet zone, and a two-wavelength outlet
// zone, the complete domain is 10*Lp. Using 20 cells per peak wavelength and 10 cells per water
// depth gives nx=200, nz=10,
// dx=Lp/20=7.3274470321 m, and dz=h/10=1 m, so dz controls the CFL time step:
//   dt_cfl = 0.5*min(dx,dz)/cp = 0.0519962817 s
// The generator selects 300 steps/Tp, giving dt=Tp/300=0.0508 s. A 240*Tp simulation lasts
// 3657.6 s. The solution is stored for the final 40*Tp, approximately 600 s, while the lighter
// surface-elevation output is stored for the final 200*Tp. A recording rate of 10 frames/Tp gives
// an interval of Tp/10=1.524 s and is exactly 30 time steps.
//
// Generate the target case:
//   uv run --project scripts python scripts/setup_wave_tank.py irregular \
//       --water-depth 10.0 --peak-period 15.24 --significant-wave-height 1.0 \
//       --tank-length 700 --cells-per-wavelength 20 --cells-per-depth 10 --cfl 0.5 \
//       --inlet-buffer-wavelengths 1 --outlet-buffer-wavelengths 1 \
//       --simulation-periods 240 --solution-record-periods 40 \
//       --surface-elevation-record-periods 200 --frames-per-period 10 \
//       --component-count 100 --random-seed 20260805 \
//       --input-file irregular_wave_tank.yaml --output irregular_wave_tank.h5
//
// Run the generated case on one MPI rank:
//   "${app}" -options_file_yaml irregular_wave_tank.yaml
//
// Run the same case on four MPI ranks:
//   mpiexec -n 4 "${app}" -options_file_yaml irregular_wave_tank.yaml \
//       -output irregular_wave_tank_mpi.h5
