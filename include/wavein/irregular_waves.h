#pragma once

#include <petscvec.h>

#include <vector>

#include "wave_spectrum.h"

namespace wavein
{

struct IrregularWaveComponent
{
        PetscReal omega;
        PetscReal wavenumber;
        PetscReal amplitude;
        PetscReal phase;
};

class IrregularWaves
{
    public:
        IrregularWaves() = delete;

        IrregularWaves(MPI_Comm comm, const WaveSpectrum &spectrum, PetscReal h)
            : comm_(comm), spectrum_(spectrum), h_(h)
        {
        }

        IrregularWaves(MPI_Comm comm, const WaveSpectrum &spectrum, PetscReal h,
                       PetscReal omega_min, PetscReal omega_max, PetscInt component_count,
                       unsigned long random_seed);

        IrregularWaves(MPI_Comm, const WaveSpectrum &&, PetscReal) = delete;
        IrregularWaves(MPI_Comm, const WaveSpectrum &&, PetscReal, PetscReal, PetscReal, PetscInt,
                       unsigned long) = delete;

        IrregularWaves(const IrregularWaves &) = delete;
        IrregularWaves &operator=(const IrregularWaves &) = delete;
        ~IrregularWaves() = default;

        PetscErrorCode velocity_spectra(Vec omegas, PetscReal z, Vec horizontal_velocity_spectrum,
                                        Vec vertical_velocity_spectrum) const;

        [[nodiscard]] PetscReal water_depth() const
        {
            return h_;
        }

        [[nodiscard]] const std::vector<IrregularWaveComponent> &components() const noexcept
        {
            return components_;
        }

    private:
        const MPI_Comm comm_;          // MPI communicator
        const WaveSpectrum &spectrum_; // reference to wave spectrum
        const PetscReal h_;            // water depth
        std::vector<IrregularWaveComponent> components_;
};

} // namespace wavein
