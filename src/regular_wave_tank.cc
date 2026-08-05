#include "wavein/regular_wave_tank.h"
#include <petscdmstag.h>

namespace wavein
{

RegularWaveTank::RegularWaveTank(MPI_Comm comm, DM dm, const AiryWave &wave, Wavemaker &wavemaker,
                                 Projection &projection, PetscReal rhow)
    : comm_(comm), wave_(wave), wavemaker_(wavemaker), projection_(projection), rhow_(rhow)
{
    PetscFunctionBeginUser;

    PetscCallAbort(comm_, DMStagGetBoundaryTypes(dm, &xboundary_, nullptr, nullptr));

    // Create vectors for the reference solution and surface elevation
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &ref_sol_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &ref_eta_));

    // Create a matrix to extract the vertical velocity at the top boundary
    PetscCallAbort(comm_, DMCreateMatrix(dm, &mat_extract_top_w_));

    // Create vectors for the vertical velocity and pressure at the top boundary
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &top_w_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &top_p_));

    // Create vectors for the cosine and sine components of the reference solution and surface
    // elevation
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &ref_sol_cos_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &ref_sol_sin_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &ref_eta_cos_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &ref_eta_sin_));

    const auto geometry = wavemaker_.forcing_zone_geometry();

    Vec ref_sol_cos_local = nullptr;
    Vec ref_sol_sin_local = nullptr;
    Vec ref_eta_cos_local = nullptr;
    Vec ref_eta_sin_local = nullptr;

    PetscCallAbort(comm_, DMCreateLocalVector(dm, &ref_sol_cos_local));
    PetscCallAbort(comm_, DMCreateLocalVector(dm, &ref_sol_sin_local));
    PetscCallAbort(comm_, DMCreateLocalVector(dm, &ref_eta_cos_local));
    PetscCallAbort(comm_, DMCreateLocalVector(dm, &ref_eta_sin_local));

    PetscCallAbort(comm_, VecZeroEntries(ref_sol_cos_local));
    PetscCallAbort(comm_, VecZeroEntries(ref_sol_sin_local));
    PetscCallAbort(comm_, VecZeroEntries(ref_eta_cos_local));
    PetscCallAbort(comm_, VecZeroEntries(ref_eta_sin_local));

    PetscInt iu, iv;
    PetscCallAbort(comm_, DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iu));
    PetscCallAbort(comm_, DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iv));

    PetscInt iprev, icenter;
    PetscCallAbort(comm_, DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
    PetscCallAbort(comm_, DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));

    PetscReal **c_arr_x = nullptr;
    PetscReal **c_arr_z = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm_, DMStagGetProductCoordinateArraysRead(dm, &c_arr_x, &c_arr_z, nullptr));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscReal ***c_arr_ref_sol_cos_local = nullptr;
    PetscReal ***c_arr_ref_sol_sin_local = nullptr;
    PetscReal ***c_arr_ref_eta_cos_local = nullptr;
    PetscReal ***c_arr_ref_eta_sin_local = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm_, DMStagVecGetArray(dm, ref_sol_cos_local, &c_arr_ref_sol_cos_local));
    PetscCallAbort(comm_, DMStagVecGetArray(dm, ref_sol_sin_local, &c_arr_ref_sol_sin_local));
    PetscCallAbort(comm_, DMStagVecGetArray(dm, ref_eta_cos_local, &c_arr_ref_eta_cos_local));
    PetscCallAbort(comm_, DMStagVecGetArray(dm, ref_eta_sin_local, &c_arr_ref_eta_sin_local));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscInt ndim[2];
    PetscInt startx, starty, nx, ny, n_extra[2], ex, ey;

    PetscCallAbort(comm_, DMStagGetGlobalSizes(dm, &ndim[0], &ndim[1], nullptr));
    PetscCallAbort(comm_, DMStagGetCorners(dm, &startx, &starty, nullptr, &nx, &ny, nullptr,
                                           &n_extra[0], &n_extra[1], nullptr));

    const PetscReal H = wave_.wave_height();
    const PetscReal k = wave_.wavenumber();

    for (ey = starty; ey != starty + ny; ++ey)
    {
        const PetscReal z = c_arr_z[ey][icenter];
        const PetscReal transfer =
            0.5 * H * PetscRealPartComplex(wave_.horizontal_velocity_transfer(z));

        for (ex = startx; ex != startx + nx + n_extra[0]; ++ex)
        {
            PetscReal x = c_arr_x[ex][iprev];

            if (!geometry.is_in_forcing_zone(x))
            {
                continue;
            }

            c_arr_ref_sol_cos_local[ey][ex][iu] = transfer * PetscCosReal(k * x);
            c_arr_ref_sol_sin_local[ey][ex][iu] = transfer * PetscSinReal(k * x);
        }
    }

    for (ey = starty; ey != starty + ny + n_extra[1]; ++ey)
    {
        const PetscReal z = c_arr_z[ey][iprev];
        const PetscReal transfer =
            0.5 * H * PetscImaginaryPartComplex(wave_.vertical_velocity_transfer(z));

        for (ex = startx; ex != startx + nx; ++ex)
        {
            PetscReal x = c_arr_x[ex][icenter];

            if (!geometry.is_in_forcing_zone(x))
            {
                continue;
            }

            c_arr_ref_sol_cos_local[ey][ex][iv] = transfer * PetscSinReal(k * x);
            c_arr_ref_sol_sin_local[ey][ex][iv] = -transfer * PetscCosReal(k * x);

            if (ey == ndim[1])
            {
                c_arr_ref_eta_cos_local[ey][ex][iv] = 0.5 * H * PetscCosReal(k * x);
                c_arr_ref_eta_sin_local[ey][ex][iv] = 0.5 * H * PetscSinReal(k * x);
            }
        }
    }
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm_, DMStagVecRestoreArray(dm, ref_sol_cos_local, &c_arr_ref_sol_cos_local));
    PetscCallAbort(comm_, DMStagVecRestoreArray(dm, ref_sol_sin_local, &c_arr_ref_sol_sin_local));
    PetscCallAbort(comm_, DMStagVecRestoreArray(dm, ref_eta_cos_local, &c_arr_ref_eta_cos_local));
    PetscCallAbort(comm_, DMStagVecRestoreArray(dm, ref_eta_sin_local, &c_arr_ref_eta_sin_local));
    PetscCallAbort(comm_,
                   DMStagRestoreProductCoordinateArraysRead(dm, &c_arr_x, &c_arr_z, nullptr));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscCallAbort(comm_, DMLocalToGlobal(dm, ref_sol_cos_local, INSERT_VALUES, ref_sol_cos_));
    PetscCallAbort(comm_, DMLocalToGlobal(dm, ref_sol_sin_local, INSERT_VALUES, ref_sol_sin_));
    PetscCallAbort(comm_, DMLocalToGlobal(dm, ref_eta_cos_local, INSERT_VALUES, ref_eta_cos_));
    PetscCallAbort(comm_, DMLocalToGlobal(dm, ref_eta_sin_local, INSERT_VALUES, ref_eta_sin_));

    PetscCallAbort(comm_, VecDestroy(&ref_sol_cos_local));
    PetscCallAbort(comm_, VecDestroy(&ref_sol_sin_local));
    PetscCallAbort(comm_, VecDestroy(&ref_eta_cos_local));
    PetscCallAbort(comm_, VecDestroy(&ref_eta_sin_local));

    // Assemble the matrix to extract the vertical velocity at the top boundary
    DMStagStencil row, col;
    PetscScalar val;

    for (ey = starty; ey != starty + ny; ++ey)
    {
        for (ex = startx; ex != startx + nx; ++ex)
        {
            if (ey == ndim[1] - 1)
            {
                row.i = ex;
                row.j = ey;
                row.loc = DMSTAG_UP;
                row.c = 0;

                col.i = ex;
                col.j = ey;
                col.loc = DMSTAG_UP;
                col.c = 0;

                val = 1.0;

                PetscCallAbort(comm_, DMStagMatSetValuesStencil(dm, mat_extract_top_w_, 1, &row, 1,
                                                                &col, &val, INSERT_VALUES));
            }
        }
    }
    PetscCallAbort(comm_, MatAssemblyBegin(mat_extract_top_w_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm_, MatAssemblyEnd(mat_extract_top_w_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm_,
                   MatSetOption(mat_extract_top_w_, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE));

    PetscFunctionReturnVoid();
}

PetscErrorCode RegularWaveTank::update(Vec sol, Vec eta, Vec source, PetscReal t, PetscReal dt,
                                       PetscReal factor) noexcept
{
    PetscFunctionBeginUser;

    if (xboundary_ == DM_BOUNDARY_NONE)
    {
        PetscCall(reference_fields(ref_sol_, ref_eta_, t));
        PetscCall(wavemaker_.force(ref_sol_, ref_eta_, dt, factor, sol, eta));
    }

    // Update the surface elevation and the presure at the linearized free surface
    PetscCall(MatMult(mat_extract_top_w_, sol, top_w_));
    PetscCall(VecAXPY(eta, dt, top_w_));
    PetscCall(VecZeroEntries(top_p_));
    PetscCall(VecAXPY(top_p_, rhow_ * kGA, eta));

    // Add source term
    if (source != nullptr)
    {
        PetscCall(VecAXPY(sol, 1.0, source));
    }

    // Force the velocity field to be divergence free
    PetscCall(projection_.project(sol, top_p_, dt));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RegularWaveTank::update(Vec sol, Vec eta, PetscReal t, PetscReal dt,
                                       PetscReal factor) noexcept
{
    PetscFunctionBeginUser;

    PetscCall(update(sol, eta, nullptr, t, dt, factor));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RegularWaveTank::reference_fields(Vec sol, Vec eta, PetscReal t) const noexcept
{
    PetscFunctionBeginUser;

    const PetscReal phase = wave_.wave_frequency() * t;
    const PetscReal cos_phase = PetscCosReal(phase);
    const PetscReal sin_phase = PetscSinReal(phase);

    PetscCall(VecCopy(ref_sol_cos_, sol));
    PetscCall(VecScale(sol, cos_phase));
    PetscCall(VecAXPY(sol, sin_phase, ref_sol_sin_));

    PetscCall(VecCopy(ref_eta_cos_, eta));
    PetscCall(VecScale(eta, cos_phase));
    PetscCall(VecAXPY(eta, sin_phase, ref_eta_sin_));

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace wavein
