#include "wavein/wave_tank.h"

#include <petscdmstag.h>

#include <utility>

namespace wavein
{

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
WaveTank::WaveTank(MPI_Comm comm, DM dm, const AiryWave &wave, Wavemaker &wavemaker,
                   Projection &projection, PetscReal rhow)
    : WaveTank(comm, dm, {{wave.wave_frequency(), 0.5 * wave.wave_height(), 0.0}},
               wave.water_depth(), wavemaker, projection, rhow)
{
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
WaveTank::WaveTank(MPI_Comm comm, DM dm, const IrregularWaves &waves, Wavemaker &wavemaker,
                   Projection &projection, PetscReal rhow)
    : WaveTank(comm, dm, waves.components(), waves.water_depth(), wavemaker, projection, rhow)
{
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
WaveTank::WaveTank(MPI_Comm comm, DM dm, std::vector<WaveComponent> components,
                   PetscReal water_depth, Wavemaker &wavemaker, Projection &projection,
                   PetscReal rhow)
    : comm_(comm), components_(std::move(components)), wavemaker_(wavemaker),
      projection_(projection), rhow_(rhow)
{
    PetscFunctionBeginUser;

    PetscCheckAbort(!components_.empty(), comm_, PETSC_ERR_ARG_WRONG,
                    "WaveTank requires at least one wave component");
    PetscCallAbort(comm_, DMStagGetBoundaryTypes(dm, &xboundary_, nullptr, nullptr));
    PetscCheckAbort(xboundary_ == DM_BOUNDARY_NONE || xboundary_ == DM_BOUNDARY_PERIODIC, comm_,
                    PETSC_ERR_ARG_WRONG,
                    "WaveTank requires none or periodic boundary conditions in the x-direction");

    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &ref_sol_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &ref_eta_));
    PetscCallAbort(comm_, DMCreateMatrix(dm, &mat_extract_top_w_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &top_w_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &top_p_));

    PetscCallAbort(comm_, precompute_reference_modes(dm, water_depth));

    // Assemble the matrix to extract the vertical velocity at the top boundary
    PetscInt ndim[2];
    PetscInt startx, starty, nx, ny, ex, ey;
    PetscCallAbort(comm_, DMStagGetGlobalSizes(dm, &ndim[0], &ndim[1], nullptr));
    PetscCallAbort(comm_, DMStagGetCorners(dm, &startx, &starty, nullptr, &nx, &ny, nullptr,
                                           nullptr, nullptr, nullptr));

    DMStagStencil row, col;
    PetscReal val;

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
// NOLINTEND(bugprone-easily-swappable-parameters)

PetscErrorCode WaveTank::precompute_reference_modes(DM dm, PetscReal water_depth) noexcept
{
    PetscFunctionBeginUser;

    const auto component_count = components_.size();
    ref_sol_cos_.assign(component_count, nullptr);
    ref_sol_sin_.assign(component_count, nullptr);
    ref_eta_cos_.assign(component_count, nullptr);
    ref_eta_sin_.assign(component_count, nullptr);
    component_cos_.resize(component_count);
    component_sin_.resize(component_count);

    for (std::size_t component = 0; component != component_count; ++component)
    {
        PetscCall(DMCreateGlobalVector(dm, &ref_sol_cos_[component]));
        PetscCall(DMCreateGlobalVector(dm, &ref_sol_sin_[component]));
        PetscCall(DMCreateGlobalVector(dm, &ref_eta_cos_[component]));
        PetscCall(DMCreateGlobalVector(dm, &ref_eta_sin_[component]));
    }

    Vec ref_sol_cos_local = nullptr;
    Vec ref_sol_sin_local = nullptr;
    Vec ref_eta_cos_local = nullptr;
    Vec ref_eta_sin_local = nullptr;
    PetscCall(DMGetLocalVector(dm, &ref_sol_cos_local));
    PetscCall(DMGetLocalVector(dm, &ref_sol_sin_local));
    PetscCall(DMGetLocalVector(dm, &ref_eta_cos_local));
    PetscCall(DMGetLocalVector(dm, &ref_eta_sin_local));

    PetscInt iu, iv;
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iu));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iv));

    PetscInt iprev = 0, icenter = 0;
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));

    PetscInt ndim[2];
    PetscInt startx, starty, nx, ny, n_extra[2];
    PetscCall(DMStagGetGlobalSizes(dm, &ndim[0], &ndim[1], nullptr));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, nullptr, &nx, &ny, nullptr, &n_extra[0],
                               &n_extra[1], nullptr));

    PetscReal **c_arr_x = nullptr;
    PetscReal **c_arr_z = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagGetProductCoordinateArraysRead(dm, &c_arr_x, &c_arr_z, nullptr));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    const auto geometry = wavemaker_.forcing_zone_geometry();

    for (std::size_t component = 0; component != component_count; ++component)
    {
        PetscCall(VecZeroEntries(ref_sol_cos_local));
        PetscCall(VecZeroEntries(ref_sol_sin_local));
        PetscCall(VecZeroEntries(ref_eta_cos_local));
        PetscCall(VecZeroEntries(ref_eta_sin_local));

        PetscReal ***c_arr_ref_sol_cos_local = nullptr;
        PetscReal ***c_arr_ref_sol_sin_local = nullptr;
        PetscReal ***c_arr_ref_eta_cos_local = nullptr;
        PetscReal ***c_arr_ref_eta_sin_local = nullptr;
        // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
        PetscCall(DMStagVecGetArray(dm, ref_sol_cos_local, &c_arr_ref_sol_cos_local));
        PetscCall(DMStagVecGetArray(dm, ref_sol_sin_local, &c_arr_ref_sol_sin_local));
        PetscCall(DMStagVecGetArray(dm, ref_eta_cos_local, &c_arr_ref_eta_cos_local));
        PetscCall(DMStagVecGetArray(dm, ref_eta_sin_local, &c_arr_ref_eta_sin_local));
        // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

        const PetscReal omega = components_[component].omega;
        const AiryWave wave(comm_, water_depth, 2.0 * PETSC_PI / omega);
        const PetscReal k = wave.wavenumber();

        for (PetscInt ey = starty; ey != starty + ny; ++ey)
        {
            const PetscReal z = c_arr_z[ey][icenter];
            const PetscReal transfer = PetscRealPartComplex(wave.horizontal_velocity_transfer(z));

            for (PetscInt ex = startx; ex != startx + nx + n_extra[0]; ++ex)
            {
                const PetscReal x = c_arr_x[ex][iprev];
                if (!geometry.is_in_inlet_forcing_zone(x))
                {
                    continue;
                }

                c_arr_ref_sol_cos_local[ey][ex][iu] = transfer * PetscCosReal(k * x);
                c_arr_ref_sol_sin_local[ey][ex][iu] = transfer * PetscSinReal(k * x);
            }
        }

        for (PetscInt ey = starty; ey != starty + ny + n_extra[1]; ++ey)
        {
            const PetscReal z = c_arr_z[ey][iprev];
            const PetscReal transfer =
                PetscImaginaryPartComplex(wave.vertical_velocity_transfer(z));

            for (PetscInt ex = startx; ex != startx + nx; ++ex)
            {
                const PetscReal x = c_arr_x[ex][icenter];
                if (!geometry.is_in_inlet_forcing_zone(x))
                {
                    continue;
                }

                const PetscReal cos_kx = PetscCosReal(k * x);
                const PetscReal sin_kx = PetscSinReal(k * x);
                c_arr_ref_sol_cos_local[ey][ex][iv] = transfer * sin_kx;
                c_arr_ref_sol_sin_local[ey][ex][iv] = -transfer * cos_kx;

                if (ey == ndim[1])
                {
                    c_arr_ref_eta_cos_local[ey][ex][iv] = cos_kx;
                    c_arr_ref_eta_sin_local[ey][ex][iv] = sin_kx;
                }
            }
        }

        // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
        PetscCall(DMStagVecRestoreArray(dm, ref_sol_cos_local, &c_arr_ref_sol_cos_local));
        PetscCall(DMStagVecRestoreArray(dm, ref_sol_sin_local, &c_arr_ref_sol_sin_local));
        PetscCall(DMStagVecRestoreArray(dm, ref_eta_cos_local, &c_arr_ref_eta_cos_local));
        PetscCall(DMStagVecRestoreArray(dm, ref_eta_sin_local, &c_arr_ref_eta_sin_local));
        // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

        PetscCall(DMLocalToGlobal(dm, ref_sol_cos_local, INSERT_VALUES, ref_sol_cos_[component]));
        PetscCall(DMLocalToGlobal(dm, ref_sol_sin_local, INSERT_VALUES, ref_sol_sin_[component]));
        PetscCall(DMLocalToGlobal(dm, ref_eta_cos_local, INSERT_VALUES, ref_eta_cos_[component]));
        PetscCall(DMLocalToGlobal(dm, ref_eta_sin_local, INSERT_VALUES, ref_eta_sin_[component]));
    }

    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &c_arr_x, &c_arr_z, nullptr));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscCall(DMRestoreLocalVector(dm, &ref_sol_cos_local));
    PetscCall(DMRestoreLocalVector(dm, &ref_sol_sin_local));
    PetscCall(DMRestoreLocalVector(dm, &ref_eta_cos_local));
    PetscCall(DMRestoreLocalVector(dm, &ref_eta_sin_local));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode WaveTank::update(Vec sol, Vec eta, Vec source, PetscReal t, PetscReal dt,
                                PetscReal factor) noexcept
{
    PetscFunctionBeginUser;

    if (xboundary_ == DM_BOUNDARY_NONE)
    {
        const auto component_count = components_.size();
        for (std::size_t component = 0; component != component_count; ++component)
        {
            const PetscReal phase = components_[component].omega * t + components_[component].phase;
            component_cos_[component] = components_[component].amplitude * PetscCosReal(phase);
            component_sin_[component] = components_[component].amplitude * PetscSinReal(phase);
        }

        PetscInt count = 0;
        PetscCall(PetscIntCast(component_count, &count));
        PetscCall(VecZeroEntries(ref_sol_));
        PetscCall(VecMAXPY(ref_sol_, count, component_cos_.data(), ref_sol_cos_.data()));
        PetscCall(VecMAXPY(ref_sol_, count, component_sin_.data(), ref_sol_sin_.data()));
        PetscCall(VecZeroEntries(ref_eta_));
        PetscCall(VecMAXPY(ref_eta_, count, component_cos_.data(), ref_eta_cos_.data()));
        PetscCall(VecMAXPY(ref_eta_, count, component_sin_.data(), ref_eta_sin_.data()));
        PetscCall(wavemaker_.force(ref_sol_, ref_eta_, dt, factor, sol, eta));
    }

    PetscCall(MatMult(mat_extract_top_w_, sol, top_w_));
    PetscCall(VecAXPY(eta, dt, top_w_));
    PetscCall(VecZeroEntries(top_p_));
    PetscCall(VecAXPY(top_p_, rhow_ * kGA, eta));

    if (source != nullptr)
    {
        PetscCall(VecAXPY(sol, 1.0, source));
    }

    PetscCall(projection_.project(sol, top_p_, dt));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode WaveTank::update(Vec sol, Vec eta, PetscReal t, PetscReal dt,
                                PetscReal factor) noexcept
{
    PetscFunctionBeginUser;

    PetscCall(update(sol, eta, nullptr, t, dt, factor));

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace wavein
