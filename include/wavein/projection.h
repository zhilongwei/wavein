#pragma once

#include <petscdmstag.h>
#include <petscksp.h>
#include <petscvec.h>

#include "wavein/constants.h"

namespace wavein
{

class Projection
{
    public:
        Projection() = delete;

        Projection(MPI_Comm comm, DM dm, PetscReal dx, PetscReal dz,
                   PetscReal rhow = kSeawaterDensity);

        Projection(const Projection &) = delete;
        Projection &operator=(const Projection &) = delete;
        ~Projection() noexcept
        {
            destroy();
        }

        PetscErrorCode project(Vec sol, Vec ptop, PetscReal dt) noexcept;

        PetscErrorCode divergence(Vec sol, Vec div) const noexcept;
        PetscErrorCode divergence_norm(Vec sol, PetscReal *norm) noexcept;

    private:
        void destroy() noexcept
        {
            PetscCallAbort(comm_, KSPDestroy(&ksp_));

            PetscCallAbort(comm_, VecDestroy(&rhs_));
            PetscCallAbort(comm_, VecDestroy(&pressure_top_row_element_));
            PetscCallAbort(comm_, VecDestroy(&rhs_pressure_grid_));
            PetscCallAbort(comm_, VecDestroy(&sol_pressure_grid_));
            PetscCallAbort(comm_, VecDestroy(&dvel_));
            PetscCallAbort(comm_, VecDestroy(&velocity_mask_));

            PetscCallAbort(comm_, MatDestroy(&mat_a_));
            PetscCallAbort(comm_, MatDestroy(&mat_div_vel_));
            PetscCallAbort(comm_, MatDestroy(&mat_grad_pressure_));
            PetscCallAbort(comm_, MatDestroy(&mat_d2p_dz2_top_Dirichlet_bc_));
            PetscCallAbort(comm_, MatDestroy(&mat_dp_dz_top_Dirichlet_bc_));
            PetscCallAbort(comm_, MatDestroy(&mat_migrate_top_p_up2element_));

            PetscCallAbort(comm_, DMDestroy(&dm_pressure_));
            PetscCallAbort(comm_, DMDestroy(&dm_));
        }
        MPI_Comm comm_ = MPI_COMM_NULL;
        DM dm_ = nullptr;
        PetscReal dx_ = 0.0;
        PetscReal dz_ = 0.0;
        PetscReal rhow_ = 0.0;

        DM dm_pressure_ = nullptr;
        Mat mat_a_ = nullptr;
        Mat mat_div_vel_ = nullptr;
        Mat mat_grad_pressure_ = nullptr;
        Mat mat_d2p_dz2_top_Dirichlet_bc_ = nullptr;
        Mat mat_dp_dz_top_Dirichlet_bc_ = nullptr;
        Mat mat_migrate_top_p_up2element_ = nullptr;
        Vec velocity_mask_ = nullptr;

        KSP ksp_ = nullptr;

        Vec rhs_ = nullptr;
        Vec pressure_top_row_element_ = nullptr;
        Vec rhs_pressure_grid_ = nullptr;
        Vec sol_pressure_grid_ = nullptr;
        Vec dvel_ = nullptr;
};

} // namespace wavein