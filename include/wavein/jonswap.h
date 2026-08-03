#pragma once

#include "wave_spectrum.h"

namespace wavein
{

class Jonswap : public WaveSpectrum
{
    public:
        Jonswap() = delete;
        // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
        Jonswap(PetscReal Hs, PetscReal Tp, PetscReal gamma = 3.3, PetscReal sigma_a = 0.07,
                PetscReal sigma_b = 0.09)
            : WaveSpectrum(Hs, Tp), gamma_(gamma), sigma_a_(sigma_a), sigma_b_(sigma_b){};
        Jonswap(const Jonswap &) = delete;
        Jonswap &operator=(const Jonswap &) = delete;
        ~Jonswap() override = default;

        [[nodiscard]] PetscReal spectrum(PetscReal omega) const override;

    private:
        PetscReal gamma_;
        PetscReal sigma_a_;
        PetscReal sigma_b_;
};

} // namespace wavein