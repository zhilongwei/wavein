#pragma once

#include <petscsys.h>

namespace wavein
{

inline constexpr PetscReal kGA = 9.80665;               // gravitational acceleration, m/s^2
inline constexpr PetscReal kSeawaterDensity = 1025.0;   // seawater density, kg/m^3
inline constexpr PetscReal kFreshwaterDensity = 1000.0; // freshwater density, kg/m^3

} // namespace wavein