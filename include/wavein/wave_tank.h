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

        PetscErrorCode precompute_reference_factors(PetscReal water_depth) noexcept;
        PetscErrorCode synthesize_reference_fields(PetscReal t) noexcept;

        void destroy() noexcept
        {
            PetscCallAbort(comm_, VecDestroy(&ref_sol_));
            PetscCallAbort(comm_, VecDestroy(&ref_eta_));
            PetscCallAbort(comm_, MatDestroy(&mat_extract_top_w_));
            PetscCallAbort(comm_, VecDestroy(&top_w_));
            PetscCallAbort(comm_, VecDestroy(&top_p_));
            PetscCallAbort(comm_, DMDestroy(&dm_));
        }

        MPI_Comm comm_ = MPI_COMM_NULL;
        DM dm_ = nullptr;
        const std::vector<WaveComponent> components_;
        Wavemaker &wavemaker_;
        Projection &projection_;
        PetscReal rhow_ = 0.0;
        DMBoundaryType xboundary_ = DM_BOUNDARY_NONE;

        std::vector<PetscInt> inlet_u_x_indices_;
        std::vector<PetscInt> inlet_w_x_indices_;

        // Component-major separable factors. Spatial phases contain exp(-i*k*x), and the
        // transfer arrays contain the corresponding complex Airy velocity factors.
        std::vector<PetscComplex> spatial_phase_u_;
        std::vector<PetscComplex> spatial_phase_w_;
        std::vector<PetscComplex> horizontal_transfer_;
        std::vector<PetscComplex> vertical_transfer_;

        Vec ref_sol_ = nullptr;
        Vec ref_eta_ = nullptr;

        Mat mat_extract_top_w_ = nullptr;
        Vec top_w_ = nullptr;
        Vec top_p_ = nullptr;
};

} // namespace wavein
