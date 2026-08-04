#pragma once

#include <petscvec.h>

#include "wave_spectrum.h"

namespace wavein
{

class IrregularWaves
{
    public:
        IrregularWaves() = delete;

        IrregularWaves(MPI_Comm comm, const WaveSpectrum &spectrum, PetscReal h)
            : comm_(comm), spectrum_(spectrum), h_(h)
        {
        }

        IrregularWaves(MPI_Comm, const WaveSpectrum &&, PetscReal) = delete;

        IrregularWaves(const IrregularWaves &) = delete;
        IrregularWaves &operator=(const IrregularWaves &) = delete;
        ~IrregularWaves() = default;

        PetscErrorCode velocity_spectra(Vec omegas, PetscReal z, Vec horizontal_velocity_spectrum,
                                        Vec vertical_velocity_spectrum) const;

        [[nodiscard]] PetscReal water_depth() const
        {
            return h_;
        }

    private:
        const MPI_Comm comm_;          // MPI communicator
        const WaveSpectrum &spectrum_; // reference to wave spectrum
        const PetscReal h_;            // water depth
};

} // namespace wavein