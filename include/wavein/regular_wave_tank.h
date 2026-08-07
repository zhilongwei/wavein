#pragma once

#include "wavein/airy_wave.h"
#include "wavein/projection.h"
#include "wavein/wavemaker.h"

namespace wavein
{

class RegularWaveTank
{
    public:
        RegularWaveTank() = delete;

        RegularWaveTank(MPI_Comm comm, DM dm, const AiryWave &wave, Wavemaker &wavemaker,
                        Projection &projection, PetscReal rhow = kSeawaterDensity);

        RegularWaveTank(const RegularWaveTank &) = delete;
        RegularWaveTank &operator=(const RegularWaveTank &) = delete;
        ~RegularWaveTank() noexcept
        {
            destroy();
        }

        PetscErrorCode update(Vec sol, Vec eta, PetscReal t, PetscReal dt,
                              PetscReal factor) noexcept;
        PetscErrorCode update(Vec sol, Vec eta, Vec source, PetscReal t, PetscReal dt,
                              PetscReal factor) noexcept;

    private:
        void destroy() noexcept
        {
            PetscCallAbort(comm_, VecDestroy(&ref_sol_));
            PetscCallAbort(comm_, VecDestroy(&ref_eta_));
            PetscCallAbort(comm_, MatDestroy(&mat_extract_top_w_));
            PetscCallAbort(comm_, VecDestroy(&top_w_));
            PetscCallAbort(comm_, VecDestroy(&top_p_));
            PetscCallAbort(comm_, VecDestroy(&ref_sol_cos_));
            PetscCallAbort(comm_, VecDestroy(&ref_sol_sin_));
            PetscCallAbort(comm_, VecDestroy(&ref_eta_cos_));
            PetscCallAbort(comm_, VecDestroy(&ref_eta_sin_));
        }
        MPI_Comm comm_ = MPI_COMM_NULL;
        const AiryWave &wave_;
        Wavemaker &wavemaker_;
        Projection &projection_;
        PetscReal rhow_ = 0.0;

        DMBoundaryType xboundary_;

        Vec ref_sol_ = nullptr;
        Vec ref_eta_ = nullptr;

        Mat mat_extract_top_w_ = nullptr;
        Vec top_w_ = nullptr;
        Vec top_p_ = nullptr;

        Vec ref_sol_cos_ = nullptr;
        Vec ref_sol_sin_ = nullptr;
        Vec ref_eta_cos_ = nullptr;
        Vec ref_eta_sin_ = nullptr;
};

} // namespace wavein
