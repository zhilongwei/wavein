#include "wavein/wavemaker.h"

namespace wavein
{

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
Wavemaker::Wavemaker(MPI_Comm comm, DM dm, PetscReal wavelength, PetscReal xmin, PetscReal xmax,
                     PetscReal nin, PetscReal nout, PetscReal gamma)
    : comm_(comm), geometry_{xmin, xmax, nin * wavelength, nout * wavelength}
{
    PetscFunctionBeginUser;

    PetscCallAbort(comm, DMCreateGlobalVector(dm, &blender_));
    PetscCallAbort(comm, DMCreateGlobalVector(dm, &delta_vel_));
    PetscCallAbort(comm, DMCreateGlobalVector(dm, &delta_eta_));

    Vec blender_local;
    PetscCallAbort(comm, DMCreateLocalVector(dm, &blender_local));
    PetscCallAbort(comm, VecSet(blender_local, 0.0));

    PetscInt iu, iv, ip;
    PetscCallAbort(comm, DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iu));
    PetscCallAbort(comm, DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iv));
    PetscCallAbort(comm, DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));

    PetscReal **c_arr_x = nullptr;
    PetscReal **c_arr_y = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm, DMStagGetProductCoordinateArraysRead(dm, &c_arr_x, &c_arr_y, nullptr));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscInt iprev, icenter;
    PetscCallAbort(comm, DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
    PetscCallAbort(comm, DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));

    PetscReal ***c_arr_blender_local = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm, DMStagVecGetArray(dm, blender_local, &c_arr_blender_local));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)
    PetscInt startx, starty, nx, ny, n_extra[2], ex, ey;
    PetscCallAbort(comm, DMStagGetCorners(dm, &startx, &starty, nullptr, &nx, &ny, nullptr,
                                          &n_extra[0], &n_extra[1], nullptr));

    PetscReal x, xx, inlet_xx, outlet_xx;
    for (ey = starty; ey != starty + ny + n_extra[1]; ++ey)
    {
        for (ex = startx; ex != startx + nx + n_extra[0]; ++ex)
        {
            x = c_arr_x[ex][iprev];
            inlet_xx = (xmin - x) / (nin * wavelength) + 1.0;
            outlet_xx = (x - xmax) / (nout * wavelength) + 1.0;
            xx = PetscMin(PetscMax(PetscMax(inlet_xx, outlet_xx), 0.0), 1.0);
            c_arr_blender_local[ey][ex][iu] = exp_blender_func(xx);

            x = c_arr_x[ex][icenter];
            inlet_xx = (xmin - x) / (nin * wavelength) + 1.0;
            outlet_xx = (x - xmax) / (nout * wavelength) + 1.0;
            xx = PetscMin(PetscMax(PetscMax(inlet_xx, outlet_xx), 0.0), 1.0);
            c_arr_blender_local[ey][ex][iv] = exp_blender_func(xx);
        }
    }
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm, DMStagRestoreProductCoordinateArraysRead(dm, &c_arr_x, &c_arr_y, nullptr));
    PetscCallAbort(comm, DMStagVecRestoreArray(dm, blender_local, &c_arr_blender_local));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm, DMLocalToGlobal(dm, blender_local, INSERT_VALUES, blender_));
    PetscCallAbort(comm, VecDestroy(&blender_local));
    PetscCallAbort(comm, VecScale(blender_, gamma));

    PetscFunctionReturnVoid();
}
// NOLINTEND(bugprone-easily-swappable-parameters)

Wavemaker::~Wavemaker() noexcept
{
    PetscFunctionBeginUser;

    PetscCallAbort(comm_, VecDestroy(&blender_));
    PetscCallAbort(comm_, VecDestroy(&delta_vel_));
    PetscCallAbort(comm_, VecDestroy(&delta_eta_));

    PetscFunctionReturnVoid();
}

PetscErrorCode Wavemaker::force(Vec ref_vel, Vec ref_eta, PetscReal dt, PetscReal factor, Vec vel,
                                Vec eta) noexcept
{
    PetscFunctionBeginUser;

    PetscCall(VecWAXPY(delta_vel_, -1.0, ref_vel, vel));
    PetscCall(VecWAXPY(delta_eta_, -1.0, ref_eta, eta));

    PetscCall(VecPointwiseMult(delta_vel_, delta_vel_, blender_));
    PetscCall(VecPointwiseMult(delta_eta_, delta_eta_, blender_));

    PetscCall(VecAXPY(vel, factor * dt, delta_vel_));
    PetscCall(VecAXPY(eta, factor * dt, delta_eta_));

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace wavein