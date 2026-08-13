#include "wave_tank_setup.h"

namespace wavein::app
{

namespace
{

[[nodiscard]] PetscInt interval_to_steps(PetscReal interval, PetscReal dt)
{
    return static_cast<PetscInt>(PetscFloorReal(interval / dt + 0.5));
}

} // namespace

WaveTankSchedule create_wave_tank_schedule(const SimulationOptions &simulation,
                                           const WaveTankOutputOptions &output)
{
    WaveTankSchedule schedule;

    schedule.num_steps =
        interval_to_steps(simulation.end_time - simulation.start_time, simulation.dt);
    schedule.output_stride = interval_to_steps(output.flow_output_interval, simulation.dt);
    schedule.field_output_start_step = interval_to_steps(
        output.flow_field_output_start_time - simulation.start_time, simulation.dt);
    schedule.field_output_end_step =
        interval_to_steps(output.flow_field_output_end_time - simulation.start_time, simulation.dt);
    schedule.elevation_output_start_step = interval_to_steps(
        output.flow_surface_elevation_output_start_time - simulation.start_time, simulation.dt);
    schedule.elevation_output_end_step = interval_to_steps(
        output.flow_surface_elevation_output_end_time - simulation.start_time, simulation.dt);

    return schedule;
}

} // namespace wavein::app
