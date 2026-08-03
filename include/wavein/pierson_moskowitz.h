#pragma once

#include "wave_spectrum.h"

namespace wavein
{

class PiersonMoskowitz : public WaveSpectrum
{
    public:
        PiersonMoskowitz() = delete;
        // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
        PiersonMoskowitz(PetscReal Hs, PetscReal Tp) : WaveSpectrum(Hs, Tp){};
        PiersonMoskowitz(const PiersonMoskowitz &) = delete;
        PiersonMoskowitz &operator=(const PiersonMoskowitz &) = delete;
        ~PiersonMoskowitz() override = default;

        [[nodiscard]] PetscReal spectrum(PetscReal omega) const override
        {
            return omega > 0.0
                       ? 5.0 / 16.0 * Hs_ * Hs_ * PetscPowReal((2.0 * PETSC_PI / Tp_) / omega, 4) /
                             omega *
                             PetscExpReal(-1.25 * PetscPowReal((2.0 * PETSC_PI / Tp_) / omega, 4))
                       : 0.0;
        }
};

} // namespace wavein