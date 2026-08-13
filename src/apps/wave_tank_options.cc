#include "wave_tank_options.h"

namespace wavein::app
{

PetscErrorCode read_simulation_options(MPI_Comm comm, SimulationOptions &options)
{
    PetscFunctionBeginUser;

    PetscOptionsBegin(comm, "sim_", "Simulation options", nullptr);
    PetscCall(PetscOptionsReal("-start_time", "Simulation start time", nullptr, options.start_time,
                               &options.start_time, nullptr));
    PetscCall(PetscOptionsReal("-end_time", "Simulation end time", nullptr, options.end_time,
                               &options.end_time, nullptr));
    PetscCall(PetscOptionsReal("-dt", "Simulation time-step size", nullptr, options.dt, &options.dt,
                               nullptr));
    PetscOptionsEnd();

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_wave_tank_output_options(MPI_Comm comm, WaveTankOutputOptions &options)
{
    PetscFunctionBeginUser;

    PetscOptionsBegin(comm, nullptr, "Wave-tank output options", nullptr);
    PetscCall(PetscOptionsString("-output", "HDF5 output filename", nullptr, options.output_file,
                                 options.output_file, sizeof(options.output_file), nullptr));
    PetscCall(PetscOptionsReal("-flow_field_output_start_time", "Flow-field output start time",
                               nullptr, options.flow_field_output_start_time,
                               &options.flow_field_output_start_time, nullptr));
    PetscCall(PetscOptionsReal("-flow_field_output_end_time", "Flow-field output end time", nullptr,
                               options.flow_field_output_end_time,
                               &options.flow_field_output_end_time, nullptr));
    PetscCall(PetscOptionsReal("-flow_surface_elevation_output_start_time",
                               "Surface-elevation output start time", nullptr,
                               options.flow_surface_elevation_output_start_time,
                               &options.flow_surface_elevation_output_start_time, nullptr));
    PetscCall(PetscOptionsReal("-flow_surface_elevation_output_end_time",
                               "Surface-elevation output end time", nullptr,
                               options.flow_surface_elevation_output_end_time,
                               &options.flow_surface_elevation_output_end_time, nullptr));
    PetscCall(PetscOptionsReal("-flow_output_interval", "Flow output interval", nullptr,
                               options.flow_output_interval, &options.flow_output_interval,
                               nullptr));
    PetscOptionsEnd();

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_airy_wave_options(MPI_Comm comm, AiryWaveOptions &options)
{
    PetscFunctionBeginUser;

    PetscOptionsBegin(comm, nullptr, "Airy-wave options", nullptr);
    PetscCall(
        PetscOptionsReal("-wave_height", "Wave height", nullptr, options.H, &options.H, nullptr));
    PetscCall(
        PetscOptionsReal("-wave_period", "Wave period", nullptr, options.T, &options.T, nullptr));
    PetscCall(
        PetscOptionsReal("-water_depth", "Water depth", nullptr, options.h, &options.h, nullptr));
    PetscOptionsEnd();

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_irregular_wave_options(MPI_Comm comm, IrregularWaveOptions &options)
{
    PetscFunctionBeginUser;

    auto random_seed = static_cast<PetscInt>(options.random_seed);

    PetscOptionsBegin(comm, nullptr, "Irregular-wave options", nullptr);
    PetscCall(PetscOptionsReal("-significant_wave_height", "Significant wave height", nullptr,
                               options.Hs, &options.Hs, nullptr));
    PetscCall(PetscOptionsReal("-peak_period", "Peak wave period", nullptr, options.Tp, &options.Tp,
                               nullptr));
    PetscCall(PetscOptionsReal("-peak_enhancement_factor", "JONSWAP peak enhancement factor",
                               nullptr, options.gamma, &options.gamma, nullptr));
    PetscCall(
        PetscOptionsReal("-water_depth", "Water depth", nullptr, options.h, &options.h, nullptr));
    PetscCall(PetscOptionsReal("-omega_min", "Minimum angular frequency", nullptr,
                               options.omega_min, &options.omega_min, nullptr));
    PetscCall(PetscOptionsReal("-omega_max", "Maximum angular frequency", nullptr,
                               options.omega_max, &options.omega_max, nullptr));
    PetscCall(PetscOptionsInt("-component_count", "Number of wave components", nullptr,
                              options.component_count, &options.component_count, nullptr));
    PetscCall(PetscOptionsInt("-random_seed", "Random phase seed", nullptr, random_seed,
                              &random_seed, nullptr));
    PetscOptionsEnd();

    options.random_seed = static_cast<unsigned long>(random_seed);

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_wave_tank_domain_options(MPI_Comm comm, WaveTankDomainOptions &options)
{
    PetscFunctionBeginUser;

    PetscOptionsBegin(comm, nullptr, "Wave-tank domain options", nullptr);
    PetscCall(PetscOptionsReal("-xmin", "Minimum x-coordinate", nullptr, options.xmin,
                               &options.xmin, nullptr));
    PetscCall(PetscOptionsReal("-xmax", "Maximum x-coordinate", nullptr, options.xmax,
                               &options.xmax, nullptr));
    PetscOptionsEnd();

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_wave_tank_grid_options(MPI_Comm comm, WaveTankGridOptions &options)
{
    PetscFunctionBeginUser;

    PetscOptionsBegin(comm, nullptr, "Wave-tank grid options", nullptr);
    PetscCall(PetscOptionsInt("-nx", "Number of cells in the x-direction", nullptr, options.Nx,
                              &options.Nx, nullptr));
    PetscCall(PetscOptionsInt("-nz", "Number of cells in the z-direction", nullptr, options.Nz,
                              &options.Nz, nullptr));
    PetscOptionsEnd();

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_wave_tank_wavemaker_options(MPI_Comm comm, WaveTankWavemakerOptions &options)
{
    PetscFunctionBeginUser;

    PetscOptionsBegin(comm, nullptr, "Wave-tank wavemaker options", nullptr);
    PetscCall(PetscOptionsReal("-nin", "Inlet forcing-zone length in wavelengths", nullptr,
                               options.nin, &options.nin, nullptr));
    PetscCall(PetscOptionsReal("-nout", "Outlet forcing-zone length in wavelengths", nullptr,
                               options.nout, &options.nout, nullptr));
    PetscCall(PetscOptionsReal("-gamma", "Outlet source-relaxation rate", nullptr, options.gamma,
                               &options.gamma, nullptr));
    PetscCall(PetscOptionsReal("-ramp_up_time", "Wavemaker ramp-up time", nullptr,
                               options.ramp_up_time, &options.ramp_up_time, nullptr));
    PetscOptionsEnd();

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace wavein::app
