#include "wave_tank_setup.h"

namespace wavein::app
{

namespace
{

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
PetscErrorCode interval_to_steps(MPI_Comm comm, PetscReal interval, PetscReal dt,
                                 const char *option_name, PetscInt *steps)
{
    PetscFunctionBeginUser;

    const PetscReal ratio = interval / dt;
    const PetscReal nearest = PetscFloorReal(ratio + 0.5);
    const PetscReal tolerance = 100.0 * PETSC_MACHINE_EPSILON * PetscMax(1.0, PetscAbsReal(ratio));

    PetscCheck(PetscAbsReal(ratio - nearest) <= tolerance, comm, PETSC_ERR_USER_INPUT,
               "%s must be a multiple of -sim_dt", option_name);

    *steps = static_cast<PetscInt>(nearest);

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace

PetscErrorCode create_wave_tank_schedule(MPI_Comm comm, const SimulationOptions &simulation,
                                         const WaveTankOutputOptions &output,
                                         WaveTankSchedule &schedule)
{
    PetscFunctionBeginUser;

    PetscCall(interval_to_steps(comm, simulation.end_time - simulation.start_time, simulation.dt,
                                "-sim_dt", &schedule.num_steps));
    PetscCall(interval_to_steps(comm, output.flow_output_interval, simulation.dt,
                                "-flow_output_interval", &schedule.output_stride));
    PetscCheck(schedule.output_stride >= 1, comm, PETSC_ERR_USER_INPUT,
               "-flow_output_interval must be at least -sim_dt");
    PetscCall(interval_to_steps(comm, output.flow_field_output_start_time - simulation.start_time,
                                simulation.dt, "-flow_field_output_start_time",
                                &schedule.field_output_start_step));
    PetscCall(interval_to_steps(comm, output.flow_field_output_end_time - simulation.start_time,
                                simulation.dt, "-flow_field_output_end_time",
                                &schedule.field_output_end_step));
    PetscCall(interval_to_steps(
        comm, output.flow_surface_elevation_output_start_time - simulation.start_time,
        simulation.dt, "-flow_surface_elevation_output_start_time",
        &schedule.elevation_output_start_step));
    PetscCall(interval_to_steps(
        comm, output.flow_surface_elevation_output_end_time - simulation.start_time, simulation.dt,
        "-flow_surface_elevation_output_end_time", &schedule.elevation_output_end_step));

    PetscCheck((schedule.field_output_end_step - schedule.field_output_start_step) %
                       schedule.output_stride ==
                   0,
               comm, PETSC_ERR_USER_INPUT,
               "Flow-field output window must contain an integer number of output intervals");
    PetscCheck((schedule.elevation_output_end_step - schedule.elevation_output_start_step) %
                       schedule.output_stride ==
                   0,
               comm, PETSC_ERR_USER_INPUT,
               "Surface-elevation output window must contain an integer number of output "
               "intervals");

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace wavein::app
