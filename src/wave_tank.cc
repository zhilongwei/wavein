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
    : comm_(comm), dm_(dm), components_(std::move(components)), wavemaker_(wavemaker),
      projection_(projection), rhow_(rhow)
{
    PetscFunctionBeginUser;

    PetscCallAbort(comm_, PetscObjectReference(reinterpret_cast<PetscObject>(dm_)));
    PetscCheckAbort(!components_.empty(), comm_, PETSC_ERR_ARG_WRONG,
                    "WaveTank requires at least one wave component");
    PetscCallAbort(comm_, DMStagGetBoundaryTypes(dm, &xboundary_, nullptr, nullptr));
    PetscCheckAbort(xboundary_ == DM_BOUNDARY_NONE || xboundary_ == DM_BOUNDARY_PERIODIC, comm_,
                    PETSC_ERR_ARG_WRONG,
                    "WaveTank requires none or periodic boundary conditions in the x-direction");

    if (xboundary_ == DM_BOUNDARY_NONE)
    {
        PetscCallAbort(comm_, DMCreateGlobalVector(dm, &ref_sol_));
        PetscCallAbort(comm_, DMCreateGlobalVector(dm, &ref_eta_));
        PetscCallAbort(comm_, precompute_reference_factors(water_depth));
    }

    PetscCallAbort(comm_, DMCreateMatrix(dm, &mat_extract_top_w_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &top_w_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &top_p_));

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

PetscErrorCode WaveTank::precompute_reference_factors(PetscReal water_depth) noexcept
{
    PetscFunctionBeginUser;

    PetscCheck(water_depth > 0.0 && !PetscIsInfOrNanReal(water_depth), comm_,
               PETSC_ERR_ARG_OUTOFRANGE, "WaveTank water depth must be finite and positive");

    const auto component_count = components_.size();

    PetscInt iprev = 0, icenter = 0;
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm_, DMSTAG_LEFT, &iprev));
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm_, DMSTAG_ELEMENT, &icenter));

    PetscInt startx, startz, nx, nz, n_extra[2];
    PetscCall(DMStagGetCorners(dm_, &startx, &startz, nullptr, &nx, &nz, nullptr, &n_extra[0],
                               &n_extra[1], nullptr));

    PetscReal **c_arr_x = nullptr;
    PetscReal **c_arr_z = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagGetProductCoordinateArraysRead(dm_, &c_arr_x, &c_arr_z, nullptr));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    const auto geometry = wavemaker_.forcing_zone_geometry();

    for (PetscInt ex = startx; ex != startx + nx + n_extra[0]; ++ex)
    {
        if (geometry.is_in_inlet_forcing_zone(c_arr_x[ex][iprev]))
        {
            inlet_u_x_indices_.push_back(ex);
        }
    }

    for (PetscInt ex = startx; ex != startx + nx; ++ex)
    {
        if (geometry.is_in_inlet_forcing_zone(c_arr_x[ex][icenter]))
        {
            inlet_w_x_indices_.push_back(ex);
        }
    }

    const std::size_t u_x_count = inlet_u_x_indices_.size();
    const std::size_t w_x_count = inlet_w_x_indices_.size();
    const auto u_z_count = static_cast<std::size_t>(nz);
    const auto w_z_count = static_cast<std::size_t>(nz) + static_cast<std::size_t>(n_extra[1]);

    spatial_phase_u_.resize(component_count * u_x_count);
    spatial_phase_w_.resize(component_count * w_x_count);
    horizontal_transfer_.resize(component_count * u_z_count);
    vertical_transfer_.resize(component_count * w_z_count);

    for (std::size_t component = 0; component != component_count; ++component)
    {
        const PetscReal omega = components_[component].omega;
        PetscCheck(omega > 0.0 && !PetscIsInfOrNanReal(omega), comm_, PETSC_ERR_ARG_OUTOFRANGE,
                   "Wave component angular frequency must be finite and positive");
        PetscCheck(!PetscIsInfOrNanReal(components_[component].amplitude), comm_,
                   PETSC_ERR_ARG_OUTOFRANGE, "Wave component amplitude must be finite");
        PetscCheck(!PetscIsInfOrNanReal(components_[component].phase), comm_,
                   PETSC_ERR_ARG_OUTOFRANGE, "Wave component phase must be finite");

        const AiryWave wave(comm_, water_depth, 2.0 * PETSC_PI / omega);
        const PetscReal k = wave.wavenumber();

        for (std::size_t local_z = 0; local_z != u_z_count; ++local_z)
        {
            const PetscInt ey = startz + static_cast<PetscInt>(local_z);
            const PetscReal z = c_arr_z[ey][icenter];
            horizontal_transfer_[component * u_z_count + local_z] =
                wave.horizontal_velocity_transfer(z);
        }

        for (std::size_t local_z = 0; local_z != w_z_count; ++local_z)
        {
            const PetscInt ey = startz + static_cast<PetscInt>(local_z);
            const PetscReal z = c_arr_z[ey][iprev];
            vertical_transfer_[component * w_z_count + local_z] =
                wave.vertical_velocity_transfer(z);
        }

        for (std::size_t local_x = 0; local_x != u_x_count; ++local_x)
        {
            const PetscInt ex = inlet_u_x_indices_[local_x];
            spatial_phase_u_[component * u_x_count + local_x] =
                PetscExpComplex(-PETSC_i * k * c_arr_x[ex][iprev]);
        }

        for (std::size_t local_x = 0; local_x != w_x_count; ++local_x)
        {
            const PetscInt ex = inlet_w_x_indices_[local_x];
            spatial_phase_w_[component * w_x_count + local_x] =
                PetscExpComplex(-PETSC_i * k * c_arr_x[ex][icenter]);
        }
    }

    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagRestoreProductCoordinateArraysRead(dm_, &c_arr_x, &c_arr_z, nullptr));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode WaveTank::synthesize_reference_fields(PetscReal t) noexcept
{
    PetscFunctionBeginUser;

    PetscInt startz, nz, n_extra_z, ndim_z;
    PetscInt iu, iv;
    PetscCall(DMStagGetCorners(dm_, nullptr, &startz, nullptr, nullptr, &nz, nullptr, nullptr,
                               &n_extra_z, nullptr));
    PetscCall(DMStagGetGlobalSizes(dm_, nullptr, &ndim_z, nullptr));
    PetscCall(DMStagGetLocationSlot(dm_, DMSTAG_LEFT, 0, &iu));
    PetscCall(DMStagGetLocationSlot(dm_, DMSTAG_DOWN, 0, &iv));
    const PetscBool owns_top_boundary = startz + nz == ndim_z ? PETSC_TRUE : PETSC_FALSE;

    Vec ref_sol_local = nullptr;
    Vec ref_eta_local = nullptr;
    PetscCall(DMGetLocalVector(dm_, &ref_sol_local));
    PetscCall(DMGetLocalVector(dm_, &ref_eta_local));
    PetscCall(VecZeroEntries(ref_sol_local));
    PetscCall(VecZeroEntries(ref_eta_local));

    PetscReal ***c_arr_ref_sol_local = nullptr;
    PetscReal ***c_arr_ref_eta_local = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagVecGetArray(dm_, ref_sol_local, &c_arr_ref_sol_local));
    PetscCall(DMStagVecGetArray(dm_, ref_eta_local, &c_arr_ref_eta_local));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    const std::size_t u_x_count = inlet_u_x_indices_.size();
    const std::size_t w_x_count = inlet_w_x_indices_.size();
    const auto u_z_count = static_cast<std::size_t>(nz);
    const auto w_z_count = static_cast<std::size_t>(nz) + static_cast<std::size_t>(n_extra_z);

    for (std::size_t component = 0; component != components_.size(); ++component)
    {
        const WaveComponent &wave_component = components_[component];
        const PetscComplex temporal_phase =
            wave_component.amplitude *
            PetscExpComplex(PETSC_i * (wave_component.omega * t + wave_component.phase));

        for (std::size_t local_x = 0; local_x != u_x_count; ++local_x)
        {
            const PetscInt ex = inlet_u_x_indices_[local_x];
            const PetscComplex complex_phase =
                temporal_phase * spatial_phase_u_[component * u_x_count + local_x];

            for (std::size_t local_z = 0; local_z != u_z_count; ++local_z)
            {
                const PetscInt ez = startz + static_cast<PetscInt>(local_z);
                c_arr_ref_sol_local[ez][ex][iu] += PetscRealPartComplex(
                    horizontal_transfer_[component * u_z_count + local_z] * complex_phase);
            }
        }

        for (std::size_t local_x = 0; local_x != w_x_count; ++local_x)
        {
            const PetscInt ex = inlet_w_x_indices_[local_x];
            const PetscComplex complex_phase =
                temporal_phase * spatial_phase_w_[component * w_x_count + local_x];

            for (std::size_t local_z = 0; local_z != w_z_count; ++local_z)
            {
                const PetscInt ez = startz + static_cast<PetscInt>(local_z);
                c_arr_ref_sol_local[ez][ex][iv] += PetscRealPartComplex(
                    vertical_transfer_[component * w_z_count + local_z] * complex_phase);
            }

            if (owns_top_boundary)
            {
                c_arr_ref_eta_local[ndim_z][ex][iv] += PetscRealPartComplex(complex_phase);
            }
        }
    }

    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagVecRestoreArray(dm_, ref_sol_local, &c_arr_ref_sol_local));
    PetscCall(DMStagVecRestoreArray(dm_, ref_eta_local, &c_arr_ref_eta_local));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscCall(DMLocalToGlobal(dm_, ref_sol_local, INSERT_VALUES, ref_sol_));
    PetscCall(DMLocalToGlobal(dm_, ref_eta_local, INSERT_VALUES, ref_eta_));
    PetscCall(DMRestoreLocalVector(dm_, &ref_sol_local));
    PetscCall(DMRestoreLocalVector(dm_, &ref_eta_local));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode WaveTank::update(Vec sol, Vec eta, Vec source, PetscReal t, PetscReal dt,
                                PetscReal factor) noexcept
{
    PetscFunctionBeginUser;

    if (xboundary_ == DM_BOUNDARY_NONE)
    {
        PetscCall(synthesize_reference_fields(t));
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
