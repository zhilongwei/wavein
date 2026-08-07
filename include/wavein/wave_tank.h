#pragma once

#include <vector>

#include "wavein/airy_wave.h"
#include "wavein/irregular_waves.h"
#include "wavein/projection.h"
#include "wavein/wavemaker.h"

namespace wavein
{

class WaveTank
{
    public:
        WaveTank() = delete;

        WaveTank(MPI_Comm comm, DM dm, const AiryWave &wave, Wavemaker &wavemaker,
                 Projection &projection, PetscReal rhow = kSeawaterDensity);
        WaveTank(MPI_Comm comm, DM dm, const IrregularWaves &waves, Wavemaker &wavemaker,
                 Projection &projection, PetscReal rhow = kSeawaterDensity);

        WaveTank(const WaveTank &) = delete;
        WaveTank &operator=(const WaveTank &) = delete;
        ~WaveTank() noexcept
        {
            destroy();
        }

        PetscErrorCode update(Vec sol, Vec eta, PetscReal t, PetscReal dt,
                              PetscReal factor) noexcept;
        PetscErrorCode update(Vec sol, Vec eta, Vec source, PetscReal t, PetscReal dt,
                              PetscReal factor) noexcept;

    private:
        WaveTank(MPI_Comm comm, DM dm, std::vector<WaveComponent> components, PetscReal water_depth,
                 Wavemaker &wavemaker, Projection &projection, PetscReal rhow);

        PetscErrorCode precompute_reference_modes(DM dm, PetscReal water_depth) noexcept;

        void destroy() noexcept
        {
            PetscCallAbort(comm_, VecDestroy(&ref_sol_));
            PetscCallAbort(comm_, VecDestroy(&ref_eta_));
            PetscCallAbort(comm_, MatDestroy(&mat_extract_top_w_));
            PetscCallAbort(comm_, VecDestroy(&top_w_));
            PetscCallAbort(comm_, VecDestroy(&top_p_));
            for (Vec &mode : ref_sol_cos_)
            {
                PetscCallAbort(comm_, VecDestroy(&mode));
            }
            for (Vec &mode : ref_sol_sin_)
            {
                PetscCallAbort(comm_, VecDestroy(&mode));
            }
            for (Vec &mode : ref_eta_cos_)
            {
                PetscCallAbort(comm_, VecDestroy(&mode));
            }
            for (Vec &mode : ref_eta_sin_)
            {
                PetscCallAbort(comm_, VecDestroy(&mode));
            }
        }

        MPI_Comm comm_ = MPI_COMM_NULL;
        const std::vector<WaveComponent> components_;
        Wavemaker &wavemaker_;
        Projection &projection_;
        PetscReal rhow_ = 0.0;
        DMBoundaryType xboundary_ = DM_BOUNDARY_NONE;

        std::vector<Vec> ref_sol_cos_;
        std::vector<Vec> ref_sol_sin_;
        std::vector<Vec> ref_eta_cos_;
        std::vector<Vec> ref_eta_sin_;
        std::vector<PetscReal> component_cos_;
        std::vector<PetscReal> component_sin_;

        Vec ref_sol_ = nullptr;
        Vec ref_eta_ = nullptr;

        Mat mat_extract_top_w_ = nullptr;
        Vec top_w_ = nullptr;
        Vec top_p_ = nullptr;
};

} // namespace wavein
