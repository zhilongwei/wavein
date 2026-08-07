#pragma once

#include <vector>

#include "wavein/irregular_waves.h"
#include "wavein/projection.h"
#include "wavein/wavemaker.h"

namespace wavein
{

class IrregularWaveTank
{
    public:
        IrregularWaveTank() = delete;

        IrregularWaveTank(MPI_Comm comm, DM dm, const IrregularWaves &waves, Wavemaker &wavemaker,
                          Projection &projection, PetscReal rhow = kSeawaterDensity);

        IrregularWaveTank(const IrregularWaveTank &) = delete;
        IrregularWaveTank &operator=(const IrregularWaveTank &) = delete;
        ~IrregularWaveTank() noexcept
        {
            destroy();
        }

        PetscErrorCode update(Vec sol, Vec eta, PetscReal t, PetscReal dt,
                              PetscReal factor) noexcept;
        PetscErrorCode update(Vec sol, Vec eta, Vec source, PetscReal t, PetscReal dt,
                              PetscReal factor) noexcept;

        PetscErrorCode reference_fields(Vec sol, Vec eta, PetscReal t) noexcept;

    private:
        PetscErrorCode precompute_reference_modes() noexcept;

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
            PetscCallAbort(comm_, DMDestroy(&dm_));
        }

        MPI_Comm comm_ = MPI_COMM_NULL;
        DM dm_ = nullptr;
        const IrregularWaves &waves_;
        Wavemaker &wavemaker_;
        Projection &projection_;
        PetscReal rhow_ = 0.0;

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
