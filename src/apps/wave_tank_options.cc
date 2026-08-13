#include "wave_tank_options.h"

namespace wavein::app
{

PetscErrorCode read_simulation_options(MPI_Comm comm, SimulationOptions &options)
{
    PetscFunctionBeginUser;

    PetscBool start_time_set = PETSC_FALSE;
    PetscBool end_time_set = PETSC_FALSE;
    PetscBool dt_set = PETSC_FALSE;

    PetscOptionsBegin(comm, "sim_", "Simulation options", nullptr);
    PetscCall(PetscOptionsReal("-start_time", "Simulation start time", nullptr, options.start_time,
                               &options.start_time, &start_time_set));
    PetscCall(PetscOptionsReal("-end_time", "Simulation end time", nullptr, options.end_time,
                               &options.end_time, &end_time_set));
    PetscCall(PetscOptionsReal("-dt", "Simulation time-step size", nullptr, options.dt, &options.dt,
                               &dt_set));
    PetscOptionsEnd();

    PetscBool help_requested = PETSC_FALSE;
    PetscCall(PetscOptionsHasHelp(nullptr, &help_requested));
    if (help_requested)
    {
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscCheck(start_time_set, comm, PETSC_ERR_USER_INPUT,
               "Missing required option -sim_start_time");
    PetscCheck(end_time_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -sim_end_time");
    PetscCheck(dt_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -sim_dt");

    PetscCheck(options.start_time >= 0.0 && !PetscIsInfOrNanReal(options.start_time), comm,
               PETSC_ERR_USER_INPUT, "Simulation start time -sim_start_time must be nonnegative");
    PetscCheck(options.end_time > options.start_time && !PetscIsInfOrNanReal(options.end_time),
               comm, PETSC_ERR_USER_INPUT,
               "Simulation end time -sim_end_time must be later than -sim_start_time");
    PetscCheck(options.dt > 0.0 && !PetscIsInfOrNanReal(options.dt), comm, PETSC_ERR_USER_INPUT,
               "Simulation time-step size -sim_dt must be positive");

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_wave_tank_output_options(MPI_Comm comm, WaveTankOutputOptions &options)
{
    PetscFunctionBeginUser;

    PetscBool field_start_set = PETSC_FALSE;
    PetscBool field_end_set = PETSC_FALSE;
    PetscBool elevation_start_set = PETSC_FALSE;
    PetscBool elevation_end_set = PETSC_FALSE;
    PetscBool interval_set = PETSC_FALSE;

    PetscOptionsBegin(comm, nullptr, "Wave-tank output options", nullptr);
    PetscCall(PetscOptionsString("-output", "HDF5 output filename", nullptr, options.output_file,
                                 options.output_file, sizeof(options.output_file), nullptr));
    PetscCall(PetscOptionsReal("-flow_field_output_start_time", "Flow-field output start time",
                               nullptr, options.flow_field_output_start_time,
                               &options.flow_field_output_start_time, &field_start_set));
    PetscCall(PetscOptionsReal("-flow_field_output_end_time", "Flow-field output end time", nullptr,
                               options.flow_field_output_end_time,
                               &options.flow_field_output_end_time, &field_end_set));
    PetscCall(PetscOptionsReal(
        "-flow_surface_elevation_output_start_time", "Surface-elevation output start time", nullptr,
        options.flow_surface_elevation_output_start_time,
        &options.flow_surface_elevation_output_start_time, &elevation_start_set));
    PetscCall(PetscOptionsReal(
        "-flow_surface_elevation_output_end_time", "Surface-elevation output end time", nullptr,
        options.flow_surface_elevation_output_end_time,
        &options.flow_surface_elevation_output_end_time, &elevation_end_set));
    PetscCall(PetscOptionsReal("-flow_output_interval", "Flow output interval", nullptr,
                               options.flow_output_interval, &options.flow_output_interval,
                               &interval_set));
    PetscOptionsEnd();

    PetscBool help_requested = PETSC_FALSE;
    PetscCall(PetscOptionsHasHelp(nullptr, &help_requested));
    if (help_requested)
    {
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscCheck(field_start_set, comm, PETSC_ERR_USER_INPUT,
               "Missing required option -flow_field_output_start_time");
    PetscCheck(field_end_set, comm, PETSC_ERR_USER_INPUT,
               "Missing required option -flow_field_output_end_time");
    PetscCheck(elevation_start_set, comm, PETSC_ERR_USER_INPUT,
               "Missing required option -flow_surface_elevation_output_start_time");
    PetscCheck(elevation_end_set, comm, PETSC_ERR_USER_INPUT,
               "Missing required option -flow_surface_elevation_output_end_time");
    PetscCheck(interval_set, comm, PETSC_ERR_USER_INPUT,
               "Missing required option -flow_output_interval");
    PetscCheck(
        options.flow_output_interval > 0.0 && !PetscIsInfOrNanReal(options.flow_output_interval),
        comm, PETSC_ERR_USER_INPUT, "Flow output interval -flow_output_interval must be positive");

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode validate_wave_tank_output_options(MPI_Comm comm, const SimulationOptions &simulation,
                                                 const WaveTankOutputOptions &output)
{
    PetscFunctionBeginUser;

    PetscCheck(output.flow_field_output_start_time >= simulation.start_time &&
                   !PetscIsInfOrNanReal(output.flow_field_output_start_time),
               comm, PETSC_ERR_USER_INPUT,
               "Flow-field output must not start before the simulation start time");
    PetscCheck(output.flow_field_output_end_time >= output.flow_field_output_start_time &&
                   !PetscIsInfOrNanReal(output.flow_field_output_end_time),
               comm, PETSC_ERR_USER_INPUT,
               "Flow-field output end time must not precede its start time");
    PetscCheck(output.flow_field_output_end_time <= simulation.end_time, comm, PETSC_ERR_USER_INPUT,
               "Flow-field output must not end after the simulation end time");
    PetscCheck(output.flow_surface_elevation_output_start_time >= simulation.start_time &&
                   !PetscIsInfOrNanReal(output.flow_surface_elevation_output_start_time),
               comm, PETSC_ERR_USER_INPUT,
               "Surface-elevation output must not start before the simulation start time");
    PetscCheck(output.flow_surface_elevation_output_end_time >=
                       output.flow_surface_elevation_output_start_time &&
                   !PetscIsInfOrNanReal(output.flow_surface_elevation_output_end_time),
               comm, PETSC_ERR_USER_INPUT,
               "Surface-elevation output end time must not precede its start time");
    PetscCheck(output.flow_surface_elevation_output_end_time <= simulation.end_time, comm,
               PETSC_ERR_USER_INPUT,
               "Surface-elevation output must not end after the simulation end time");

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_airy_wave_options(MPI_Comm comm, AiryWaveOptions &options)
{
    PetscFunctionBeginUser;

    PetscBool height_set = PETSC_FALSE;
    PetscBool period_set = PETSC_FALSE;
    PetscBool depth_set = PETSC_FALSE;

    PetscOptionsBegin(comm, nullptr, "Airy-wave options", nullptr);
    PetscCall(PetscOptionsReal("-wave_height", "Wave height", nullptr, options.H, &options.H,
                               &height_set));
    PetscCall(PetscOptionsReal("-wave_period", "Wave period", nullptr, options.T, &options.T,
                               &period_set));
    PetscCall(PetscOptionsReal("-water_depth", "Water depth", nullptr, options.h, &options.h,
                               &depth_set));
    PetscOptionsEnd();

    PetscBool help_requested = PETSC_FALSE;
    PetscCall(PetscOptionsHasHelp(nullptr, &help_requested));
    if (help_requested)
    {
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscCheck(height_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -wave_height");
    PetscCheck(period_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -wave_period");
    PetscCheck(depth_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -water_depth");
    PetscCheck(options.H > 0.0 && !PetscIsInfOrNanReal(options.H), comm, PETSC_ERR_USER_INPUT,
               "Wave height -wave_height must be positive");
    PetscCheck(options.T > 0.0 && !PetscIsInfOrNanReal(options.T), comm, PETSC_ERR_USER_INPUT,
               "Wave period -wave_period must be positive");
    PetscCheck(options.h > 0.0 && !PetscIsInfOrNanReal(options.h), comm, PETSC_ERR_USER_INPUT,
               "Water depth -water_depth must be positive");

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_irregular_wave_options(MPI_Comm comm, IrregularWaveOptions &options)
{
    PetscFunctionBeginUser;

    PetscBool height_set = PETSC_FALSE;
    PetscBool period_set = PETSC_FALSE;
    PetscBool depth_set = PETSC_FALSE;
    PetscBool omega_min_set = PETSC_FALSE;
    PetscBool omega_max_set = PETSC_FALSE;
    PetscBool component_count_set = PETSC_FALSE;
    auto random_seed = static_cast<PetscInt>(options.random_seed);

    PetscOptionsBegin(comm, nullptr, "Irregular-wave options", nullptr);
    PetscCall(PetscOptionsReal("-significant_wave_height", "Significant wave height", nullptr,
                               options.Hs, &options.Hs, &height_set));
    PetscCall(PetscOptionsReal("-peak_period", "Peak wave period", nullptr, options.Tp, &options.Tp,
                               &period_set));
    PetscCall(PetscOptionsReal("-peak_enhancement_factor", "JONSWAP peak enhancement factor",
                               nullptr, options.gamma, &options.gamma, nullptr));
    PetscCall(PetscOptionsReal("-water_depth", "Water depth", nullptr, options.h, &options.h,
                               &depth_set));
    PetscCall(PetscOptionsReal("-omega_min", "Minimum angular frequency", nullptr,
                               options.omega_min, &options.omega_min, &omega_min_set));
    PetscCall(PetscOptionsReal("-omega_max", "Maximum angular frequency", nullptr,
                               options.omega_max, &options.omega_max, &omega_max_set));
    PetscCall(PetscOptionsInt("-component_count", "Number of wave components", nullptr,
                              options.component_count, &options.component_count,
                              &component_count_set));
    PetscCall(PetscOptionsInt("-random_seed", "Random phase seed", nullptr, random_seed,
                              &random_seed, nullptr));
    PetscOptionsEnd();

    PetscBool help_requested = PETSC_FALSE;
    PetscCall(PetscOptionsHasHelp(nullptr, &help_requested));
    if (help_requested)
    {
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscCheck(height_set, comm, PETSC_ERR_USER_INPUT,
               "Missing required option -significant_wave_height");
    PetscCheck(period_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -peak_period");
    PetscCheck(depth_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -water_depth");
    PetscCheck(omega_min_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -omega_min");
    PetscCheck(omega_max_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -omega_max");
    PetscCheck(component_count_set, comm, PETSC_ERR_USER_INPUT,
               "Missing required option -component_count");

    PetscCheck(options.Hs > 0.0 && !PetscIsInfOrNanReal(options.Hs), comm, PETSC_ERR_USER_INPUT,
               "Significant wave height -significant_wave_height must be positive");
    PetscCheck(options.Tp > 0.0 && !PetscIsInfOrNanReal(options.Tp), comm, PETSC_ERR_USER_INPUT,
               "Peak period -peak_period must be positive");
    PetscCheck(options.gamma >= 1.0 && !PetscIsInfOrNanReal(options.gamma), comm,
               PETSC_ERR_USER_INPUT,
               "Peak enhancement factor -peak_enhancement_factor must be at least 1");
    PetscCheck(options.h > 0.0 && !PetscIsInfOrNanReal(options.h), comm, PETSC_ERR_USER_INPUT,
               "Water depth -water_depth must be positive");
    PetscCheck(options.omega_min > 0.0 && !PetscIsInfOrNanReal(options.omega_min), comm,
               PETSC_ERR_USER_INPUT, "Minimum angular frequency -omega_min must be positive");
    PetscCheck(options.omega_max > options.omega_min && !PetscIsInfOrNanReal(options.omega_max),
               comm, PETSC_ERR_USER_INPUT,
               "Maximum angular frequency -omega_max must be greater than -omega_min");
    PetscCheck(options.component_count >= 1, comm, PETSC_ERR_USER_INPUT,
               "Wave component count -component_count must be at least 1");
    PetscCheck(random_seed >= 0, comm, PETSC_ERR_USER_INPUT,
               "Random phase seed -random_seed must be nonnegative");

    options.random_seed = static_cast<unsigned long>(random_seed);

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_wave_tank_domain_options(MPI_Comm comm, WaveTankDomainOptions &options)
{
    PetscFunctionBeginUser;

    PetscBool xmin_set = PETSC_FALSE;
    PetscBool xmax_set = PETSC_FALSE;

    PetscOptionsBegin(comm, nullptr, "Wave-tank domain options", nullptr);
    PetscCall(PetscOptionsReal("-xmin", "Minimum x-coordinate", nullptr, options.xmin,
                               &options.xmin, &xmin_set));
    PetscCall(PetscOptionsReal("-xmax", "Maximum x-coordinate", nullptr, options.xmax,
                               &options.xmax, &xmax_set));
    PetscOptionsEnd();

    PetscBool help_requested = PETSC_FALSE;
    PetscCall(PetscOptionsHasHelp(nullptr, &help_requested));
    if (help_requested)
    {
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscCheck(xmin_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -xmin");
    PetscCheck(xmax_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -xmax");
    PetscCheck(options.xmax > options.xmin && !PetscIsInfOrNanReal(options.xmin) &&
                   !PetscIsInfOrNanReal(options.xmax),
               comm, PETSC_ERR_USER_INPUT, "Maximum x-coordinate -xmax must be greater than -xmin");

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_wave_tank_grid_options(MPI_Comm comm, WaveTankGridOptions &options)
{
    PetscFunctionBeginUser;

    PetscBool nx_set = PETSC_FALSE;
    PetscBool nz_set = PETSC_FALSE;

    PetscOptionsBegin(comm, nullptr, "Wave-tank grid options", nullptr);
    PetscCall(PetscOptionsInt("-nx", "Number of cells in the x-direction", nullptr, options.Nx,
                              &options.Nx, &nx_set));
    PetscCall(PetscOptionsInt("-nz", "Number of cells in the z-direction", nullptr, options.Nz,
                              &options.Nz, &nz_set));
    PetscOptionsEnd();

    PetscBool help_requested = PETSC_FALSE;
    PetscCall(PetscOptionsHasHelp(nullptr, &help_requested));
    if (help_requested)
    {
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscCheck(nx_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -nx");
    PetscCheck(nz_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -nz");
    PetscCheck(options.Nx >= 2, comm, PETSC_ERR_USER_INPUT, "Cell count -nx must be at least 2");
    PetscCheck(options.Nz >= 2, comm, PETSC_ERR_USER_INPUT, "Cell count -nz must be at least 2");

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_wave_tank_wavemaker_options(MPI_Comm comm, WaveTankWavemakerOptions &options)
{
    PetscFunctionBeginUser;

    PetscBool gamma_set = PETSC_FALSE;
    PetscBool ramp_up_time_set = PETSC_FALSE;

    PetscOptionsBegin(comm, nullptr, "Wave-tank wavemaker options", nullptr);
    PetscCall(PetscOptionsReal("-nin", "Inlet forcing-zone length in wavelengths", nullptr,
                               options.nin, &options.nin, nullptr));
    PetscCall(PetscOptionsReal("-nout", "Outlet forcing-zone length in wavelengths", nullptr,
                               options.nout, &options.nout, nullptr));
    PetscCall(PetscOptionsReal("-gamma", "Outlet source-relaxation rate", nullptr, options.gamma,
                               &options.gamma, &gamma_set));
    PetscCall(PetscOptionsReal("-ramp_up_time", "Wavemaker ramp-up time", nullptr,
                               options.ramp_up_time, &options.ramp_up_time, &ramp_up_time_set));
    PetscOptionsEnd();

    PetscBool help_requested = PETSC_FALSE;
    PetscCall(PetscOptionsHasHelp(nullptr, &help_requested));
    if (help_requested)
    {
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscCheck(gamma_set, comm, PETSC_ERR_USER_INPUT, "Missing required option -gamma");
    PetscCheck(ramp_up_time_set, comm, PETSC_ERR_USER_INPUT,
               "Missing required option -ramp_up_time");
    PetscCheck(options.nin > 0.0 && !PetscIsInfOrNanReal(options.nin), comm, PETSC_ERR_USER_INPUT,
               "Inlet forcing-zone length -nin must be positive");
    PetscCheck(options.nout > 0.0 && !PetscIsInfOrNanReal(options.nout), comm, PETSC_ERR_USER_INPUT,
               "Outlet forcing-zone length -nout must be positive");
    PetscCheck(options.gamma > 0.0 && !PetscIsInfOrNanReal(options.gamma), comm,
               PETSC_ERR_USER_INPUT, "Outlet source-relaxation rate -gamma must be positive");
    PetscCheck(options.ramp_up_time > 0.0 && !PetscIsInfOrNanReal(options.ramp_up_time), comm,
               PETSC_ERR_USER_INPUT, "Wavemaker ramp-up time -ramp_up_time must be positive");

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace wavein::app
