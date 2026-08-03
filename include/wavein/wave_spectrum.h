#pragma once

#include <petscsys.h>

namespace wavein
{

class WaveSpectrum
{
    public:
        WaveSpectrum() = delete;
        // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
        WaveSpectrum(PetscReal Hs, PetscReal Tp) : Hs_(Hs), Tp_(Tp)
        {
        }
        virtual ~WaveSpectrum() = default;
        [[nodiscard]] virtual PetscReal spectrum(PetscReal omega) const = 0;

        [[nodiscard]] PetscReal significant_wave_height() const
        {
            return Hs_;
        }
        [[nodiscard]] PetscReal peak_period() const
        {
            return Tp_;
        }

    protected:
        PetscReal Hs_; // Significant wave height
        PetscReal Tp_; // Peak period
};

} // namespace wavein