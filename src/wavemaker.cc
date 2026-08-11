#include "wavein/wavemaker.h"

namespace wavein
{

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
Wavemaker::Wavemaker(MPI_Comm comm, DM dm, PetscReal wavelength, PetscReal xmin, PetscReal xmax,
                     PetscReal nin, PetscReal nout, PetscReal gamma)
    : comm_(comm), geometry_{xmin, xmax, nin * wavelength, nout * wavelength}, gamma_(gamma)
{
    PetscFunctionBeginUser;

    PetscCallAbort(comm, DMCreateGlobalVector(dm, &inlet_blender_));
    PetscCallAbort(comm, DMCreateGlobalVector(dm, &outlet_blender_));
    PetscCallAbort(comm, DMCreateGlobalVector(dm, &delta_vel_));
    PetscCallAbort(comm, DMCreateGlobalVector(dm, &delta_eta_));
    PetscCallAbort(comm, VecZeroEntries(delta_vel_));
    PetscCallAbort(comm, VecZeroEntries(delta_eta_));

    Vec inlet_blender_local, outlet_blender_local;
    PetscCallAbort(comm, DMCreateLocalVector(dm, &inlet_blender_local));
    PetscCallAbort(comm, DMCreateLocalVector(dm, &outlet_blender_local));
    PetscCallAbort(comm, VecSet(inlet_blender_local, 0.0));
    PetscCallAbort(comm, VecSet(outlet_blender_local, 0.0));

    PetscInt iu, iv;
    PetscCallAbort(comm, DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iu));
    PetscCallAbort(comm, DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iv));

    PetscReal **c_arr_x = nullptr;
    PetscReal **c_arr_y = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm, DMStagGetProductCoordinateArraysRead(dm, &c_arr_x, &c_arr_y, nullptr));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscInt iprev, icenter;
    PetscCallAbort(comm, DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
    PetscCallAbort(comm, DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));

    PetscReal ***c_arr_inlet_blender_local = nullptr;
    PetscReal ***c_arr_outlet_blender_local = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm, DMStagVecGetArray(dm, inlet_blender_local, &c_arr_inlet_blender_local));
    PetscCallAbort(comm, DMStagVecGetArray(dm, outlet_blender_local, &c_arr_outlet_blender_local));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)
    PetscInt startx, starty, nx, ny, n_extra[2], ex, ey;
    PetscCallAbort(comm, DMStagGetCorners(dm, &startx, &starty, nullptr, &nx, &ny, nullptr,
                                          &n_extra[0], &n_extra[1], nullptr));

    PetscReal x, inlet_xx, outlet_xx;
    for (ey = starty; ey != starty + ny + n_extra[1]; ++ey)
    {
        for (ex = startx; ex != startx + nx + n_extra[0]; ++ex)
        {
            x = c_arr_x[ex][iprev];
            inlet_xx = PetscMin(
                PetscMax((geometry_.xmin + geometry_.inlet_length - x) / geometry_.inlet_length,
                         0.0),
                1.0);
            outlet_xx = PetscMin(
                PetscMax((x - (geometry_.xmax - geometry_.outlet_length)) / geometry_.outlet_length,
                         0.0),
                1.0);
            c_arr_inlet_blender_local[ey][ex][iu] = exp_blender_func(inlet_xx);
            c_arr_outlet_blender_local[ey][ex][iu] = exp_blender_func(outlet_xx);

            x = c_arr_x[ex][icenter];
            inlet_xx = PetscMin(
                PetscMax((geometry_.xmin + geometry_.inlet_length - x) / geometry_.inlet_length,
                         0.0),
                1.0);
            outlet_xx = PetscMin(
                PetscMax((x - (geometry_.xmax - geometry_.outlet_length)) / geometry_.outlet_length,
                         0.0),
                1.0);
            c_arr_inlet_blender_local[ey][ex][iv] = exp_blender_func(inlet_xx);
            c_arr_outlet_blender_local[ey][ex][iv] = exp_blender_func(outlet_xx);
        }
    }
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm, DMStagRestoreProductCoordinateArraysRead(dm, &c_arr_x, &c_arr_y, nullptr));
    PetscCallAbort(comm,
                   DMStagVecRestoreArray(dm, inlet_blender_local, &c_arr_inlet_blender_local));
    PetscCallAbort(comm,
                   DMStagVecRestoreArray(dm, outlet_blender_local, &c_arr_outlet_blender_local));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm, DMLocalToGlobal(dm, inlet_blender_local, INSERT_VALUES, inlet_blender_));
    PetscCallAbort(comm, DMLocalToGlobal(dm, outlet_blender_local, INSERT_VALUES, outlet_blender_));
    PetscCallAbort(comm, VecDestroy(&inlet_blender_local));
    PetscCallAbort(comm, VecDestroy(&outlet_blender_local));

    PetscFunctionReturnVoid();
}
// NOLINTEND(bugprone-easily-swappable-parameters)

Wavemaker::~Wavemaker() noexcept
{
    PetscFunctionBeginUser;

    PetscCallAbort(comm_, VecDestroy(&inlet_blender_));
    PetscCallAbort(comm_, VecDestroy(&outlet_blender_));
    PetscCallAbort(comm_, VecDestroy(&delta_vel_));
    PetscCallAbort(comm_, VecDestroy(&delta_eta_));

    PetscFunctionReturnVoid();
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
PetscErrorCode Wavemaker::force(Vec ref_vel, Vec ref_eta, PetscReal dt, PetscReal ramp, Vec vel,
                                Vec eta) noexcept
{
    PetscFunctionBeginUser;

    // Generate waves by directly blending towards the ramped inlet reference.
    PetscCall(VecAXPBYPCZ(delta_vel_, ramp, -1.0, 0.0, ref_vel, vel));
    PetscCall(VecAXPBYPCZ(delta_eta_, ramp, -1.0, 0.0, ref_eta, eta));

    PetscCall(VecPointwiseMult(delta_vel_, delta_vel_, inlet_blender_));
    PetscCall(VecPointwiseMult(delta_eta_, delta_eta_, inlet_blender_));

    PetscCall(VecAXPY(vel, 1.0, delta_vel_));
    PetscCall(VecAXPY(eta, 1.0, delta_eta_));

    // Absorb outgoing waves with the legacy zero-target source terms in the outlet zone.
    PetscCall(VecPointwiseMult(delta_vel_, vel, outlet_blender_));
    PetscCall(VecPointwiseMult(delta_eta_, eta, outlet_blender_));
    PetscCall(VecAXPY(vel, -gamma_ * dt, delta_vel_));
    PetscCall(VecAXPY(eta, -gamma_ * dt, delta_eta_));

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace wavein
