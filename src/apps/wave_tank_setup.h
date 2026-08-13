#pragma once

#include "wave_tank_options.h"

namespace wavein::app
{

struct WaveTankSchedule
{
        PetscInt num_steps = 0;
        PetscInt output_stride = 0;
        PetscInt field_output_start_step = 0;
        PetscInt field_output_end_step = 0;
        PetscInt elevation_output_start_step = 0;
        PetscInt elevation_output_end_step = 0;
};

[[nodiscard]] WaveTankSchedule create_wave_tank_schedule(const SimulationOptions &simulation,
                                                         const WaveTankOutputOptions &output);

} // namespace wavein::app
