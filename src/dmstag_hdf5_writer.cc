#include "wavein/dmstag_hdf5_writer.h"

#include <petscviewerhdf5.h>

namespace
{

constexpr const char *kDMStagGroup = "dmstag";
constexpr PetscInt kLocationU = 0;
constexpr PetscInt kLocationW = 1;
constexpr PetscInt kLocationP = 2;

} // namespace

namespace wavein
{

DMStagHDF5Writer::DMStagHDF5Writer(MPI_Comm comm, DM dm, PetscViewer viewer)
    : comm_(comm), viewer_(viewer)
{
    PetscFunctionBeginUser;

    PetscCallAbort(comm_, PetscObjectReference(reinterpret_cast<PetscObject>(viewer_)));

    Vec prototype = nullptr;
    PetscCallAbort(comm_, DMCreateGlobalVector(dm, &prototype));
    PetscCallAbort(comm_, PetscViewerHDF5PushGroup(viewer_, kDMStagGroup));
    PetscCallAbort(comm_, write_metadata(dm, prototype));
    PetscCallAbort(comm_, PetscViewerHDF5PopGroup(viewer_));
    PetscCallAbort(comm_, VecCreateMPI(comm_, PETSC_DECIDE, 1, &sol_time_));
    PetscCallAbort(comm_, PetscObjectSetName(reinterpret_cast<PetscObject>(sol_time_), "sol_time"));
    PetscCallAbort(comm_, create_eta_field(dm, prototype));
    PetscCallAbort(comm_, VecDestroy(&prototype));

    PetscFunctionReturnVoid();
}

PetscErrorCode DMStagHDF5Writer::write_metadata(DM dm, Vec prototype)
{
    PetscFunctionBeginUser;

    Vec x = nullptr;
    Vec z = nullptr;
    Vec location = nullptr;
    Vec logical_i = nullptr;
    Vec logical_j = nullptr;
    PetscCall(VecDuplicate(prototype, &x));
    PetscCall(VecDuplicate(prototype, &z));
    PetscCall(VecDuplicate(prototype, &location));
    PetscCall(VecDuplicate(prototype, &logical_i));
    PetscCall(VecDuplicate(prototype, &logical_j));
    PetscCall(PetscObjectSetName(reinterpret_cast<PetscObject>(x), "x"));
    PetscCall(PetscObjectSetName(reinterpret_cast<PetscObject>(z), "z"));
    PetscCall(PetscObjectSetName(reinterpret_cast<PetscObject>(location), "location"));
    PetscCall(PetscObjectSetName(reinterpret_cast<PetscObject>(logical_i), "i"));
    PetscCall(PetscObjectSetName(reinterpret_cast<PetscObject>(logical_j), "j"));

    Vec x_local = nullptr;
    Vec z_local = nullptr;
    Vec location_local = nullptr;
    Vec logical_i_local = nullptr;
    Vec logical_j_local = nullptr;
    PetscCall(DMCreateLocalVector(dm, &x_local));
    PetscCall(DMCreateLocalVector(dm, &z_local));
    PetscCall(DMCreateLocalVector(dm, &location_local));
    PetscCall(DMCreateLocalVector(dm, &logical_i_local));
    PetscCall(DMCreateLocalVector(dm, &logical_j_local));
    PetscCall(VecZeroEntries(x_local));
    PetscCall(VecZeroEntries(z_local));
    PetscCall(VecZeroEntries(location_local));
    PetscCall(VecZeroEntries(logical_i_local));
    PetscCall(VecZeroEntries(logical_j_local));

    PetscInt iu, iv, ip;
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iu));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iv));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));

    PetscInt iprev, icenter;
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));

    PetscReal **c_arr_coord_x = nullptr;
    PetscReal **c_arr_coord_z = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagGetProductCoordinateArraysRead(dm, &c_arr_coord_x, &c_arr_coord_z, nullptr));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscReal ***c_arr_x_local = nullptr;
    PetscReal ***c_arr_z_local = nullptr;
    PetscReal ***c_arr_location_local = nullptr;
    PetscReal ***c_arr_logical_i_local = nullptr;
    PetscReal ***c_arr_logical_j_local = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagVecGetArray(dm, x_local, &c_arr_x_local));
    PetscCall(DMStagVecGetArray(dm, z_local, &c_arr_z_local));
    PetscCall(DMStagVecGetArray(dm, location_local, &c_arr_location_local));
    PetscCall(DMStagVecGetArray(dm, logical_i_local, &c_arr_logical_i_local));
    PetscCall(DMStagVecGetArray(dm, logical_j_local, &c_arr_logical_j_local));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscInt ndim[2];
    PetscInt startx, starty, nx, ny, n_extra[2], ex, ey;
    PetscCall(DMStagGetGlobalSizes(dm, &ndim[0], &ndim[1], nullptr));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, nullptr, &nx, &ny, nullptr, &n_extra[0],
                               &n_extra[1], nullptr));

    for (ey = starty; ey != starty + ny; ++ey)
    {
        for (ex = startx; ex != startx + nx + n_extra[0]; ++ex)
        {
            c_arr_x_local[ey][ex][iu] = c_arr_coord_x[ex][iprev];
            c_arr_z_local[ey][ex][iu] = c_arr_coord_z[ey][icenter];
            c_arr_location_local[ey][ex][iu] = kLocationU;
            c_arr_logical_i_local[ey][ex][iu] = ex;
            c_arr_logical_j_local[ey][ex][iu] = ey;
        }
    }

    for (ey = starty; ey != starty + ny + n_extra[1]; ++ey)
    {
        for (ex = startx; ex != startx + nx; ++ex)
        {
            c_arr_x_local[ey][ex][iv] = c_arr_coord_x[ex][icenter];
            c_arr_z_local[ey][ex][iv] = c_arr_coord_z[ey][iprev];
            c_arr_location_local[ey][ex][iv] = kLocationW;
            c_arr_logical_i_local[ey][ex][iv] = ex;
            c_arr_logical_j_local[ey][ex][iv] = ey;
        }
    }

    for (ey = starty; ey != starty + ny; ++ey)
    {
        for (ex = startx; ex != startx + nx; ++ex)
        {
            c_arr_x_local[ey][ex][ip] = c_arr_coord_x[ex][icenter];
            c_arr_z_local[ey][ex][ip] = c_arr_coord_z[ey][icenter];
            c_arr_location_local[ey][ex][ip] = kLocationP;
            c_arr_logical_i_local[ey][ex][ip] = ex;
            c_arr_logical_j_local[ey][ex][ip] = ey;
        }
    }

    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCall(DMStagVecRestoreArray(dm, x_local, &c_arr_x_local));
    PetscCall(DMStagVecRestoreArray(dm, z_local, &c_arr_z_local));
    PetscCall(DMStagVecRestoreArray(dm, location_local, &c_arr_location_local));
    PetscCall(DMStagVecRestoreArray(dm, logical_i_local, &c_arr_logical_i_local));
    PetscCall(DMStagVecRestoreArray(dm, logical_j_local, &c_arr_logical_j_local));
    PetscCall(
        DMStagRestoreProductCoordinateArraysRead(dm, &c_arr_coord_x, &c_arr_coord_z, nullptr));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    PetscCall(DMLocalToGlobal(dm, x_local, INSERT_VALUES, x));
    PetscCall(DMLocalToGlobal(dm, z_local, INSERT_VALUES, z));
    PetscCall(DMLocalToGlobal(dm, location_local, INSERT_VALUES, location));
    PetscCall(DMLocalToGlobal(dm, logical_i_local, INSERT_VALUES, logical_i));
    PetscCall(DMLocalToGlobal(dm, logical_j_local, INSERT_VALUES, logical_j));

    PetscCall(VecDestroy(&logical_j_local));
    PetscCall(VecDestroy(&logical_i_local));
    PetscCall(VecDestroy(&location_local));
    PetscCall(VecDestroy(&z_local));
    PetscCall(VecDestroy(&x_local));

    PetscCall(VecView(x, viewer_));
    PetscCall(VecView(z, viewer_));
    PetscCall(VecView(location, viewer_));
    PetscCall(VecView(logical_i, viewer_));
    PetscCall(VecView(logical_j, viewer_));

    PetscCall(PetscViewerHDF5WriteAttribute(viewer_, nullptr, "nx", PETSC_INT, &ndim[0]));
    PetscCall(PetscViewerHDF5WriteAttribute(viewer_, nullptr, "nz", PETSC_INT, &ndim[1]));
    PetscCall(PetscViewerHDF5WriteObjectAttribute(viewer_, reinterpret_cast<PetscObject>(location),
                                                  "u", PETSC_INT, &kLocationU));
    PetscCall(PetscViewerHDF5WriteObjectAttribute(viewer_, reinterpret_cast<PetscObject>(location),
                                                  "w", PETSC_INT, &kLocationW));
    PetscCall(PetscViewerHDF5WriteObjectAttribute(viewer_, reinterpret_cast<PetscObject>(location),
                                                  "p", PETSC_INT, &kLocationP));
    PetscCall(VecDestroy(&logical_j));
    PetscCall(VecDestroy(&logical_i));
    PetscCall(VecDestroy(&location));
    PetscCall(VecDestroy(&z));
    PetscCall(VecDestroy(&x));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DMStagHDF5Writer::create_eta_field(DM dm, Vec prototype)
{
    PetscFunctionBeginUser;

    PetscInt ndim[2];
    PetscCall(DMStagGetGlobalSizes(dm, &ndim[0], &ndim[1], nullptr));
    PetscCall(VecCreateMPI(comm_, PETSC_DECIDE, ndim[0], &eta_.values));
    PetscCall(PetscObjectSetName(reinterpret_cast<PetscObject>(eta_.values), "eta"));

    PetscInt startx, starty, nx, ny;
    PetscCall(DMStagGetCorners(dm, &startx, &starty, nullptr, &nx, &ny, nullptr, nullptr, nullptr,
                               nullptr));
    const PetscBool owns_top_row = starty + ny == ndim[1] ? PETSC_TRUE : PETSC_FALSE;
    const PetscInt count = owns_top_row ? nx : 0;

    DMStagStencil *stencils = nullptr;
    PetscInt *source_indices_local = nullptr;
    PetscInt *source_indices_global = nullptr;
    PetscCall(PetscMalloc3(count, &stencils, count, &source_indices_local, count,
                           &source_indices_global));
    if (owns_top_row)
    {
        const PetscInt ey = ndim[1] - 1;
        for (PetscInt i = 0; i != count; ++i)
        {
            const PetscInt ex = startx + i;
            stencils[i] = {DMSTAG_UP, ex, ey, 0, 0};
        }
    }

    PetscCall(DMStagStencilToIndexLocal(dm, 2, count, stencils, source_indices_local));
    ISLocalToGlobalMapping local_to_global = nullptr;
    PetscCall(DMGetLocalToGlobalMapping(dm, &local_to_global));
    PetscCall(ISLocalToGlobalMappingApply(local_to_global, count, source_indices_local,
                                          source_indices_global));
    for (PetscInt i = 0; i < count; ++i)
    {
        PetscCheck(source_indices_global[i] >= 0, comm_, PETSC_ERR_PLIB,
                   "Top-row eta point has no corresponding global DMSTAG entry");
    }

    IS source_indices = nullptr;
    IS target_indices = nullptr;
    PetscCall(
        ISCreateGeneral(comm_, count, source_indices_global, PETSC_COPY_VALUES, &source_indices));
    PetscCall(ISCreateStride(comm_, count, startx, 1, &target_indices));
    PetscCall(
        VecScatterCreate(prototype, source_indices, eta_.values, target_indices, &eta_.scatter));
    PetscCall(ISDestroy(&target_indices));
    PetscCall(ISDestroy(&source_indices));
    PetscCall(PetscFree3(stencils, source_indices_local, source_indices_global));

    PetscCall(VecCreateMPI(comm_, PETSC_DECIDE, 1, &eta_.time));
    PetscCall(PetscObjectSetName(reinterpret_cast<PetscObject>(eta_.time), "eta_time"));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DMStagHDF5Writer::push_group()
{
    PetscFunctionBeginUser;

    PetscCall(PetscViewerHDF5PushGroup(viewer_, kDMStagGroup));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DMStagHDF5Writer::pop_group()
{
    PetscFunctionBeginUser;

    PetscCall(PetscViewerHDF5PopGroup(viewer_));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DMStagHDF5Writer::write(const AiryWave &wave)
{
    PetscFunctionBeginUser;

    PetscReal value = wave.wave_height();
    PetscCall(PetscViewerHDF5WriteAttribute(viewer_, nullptr, "height", PETSC_REAL, &value));
    value = wave.water_depth();
    PetscCall(PetscViewerHDF5WriteAttribute(viewer_, nullptr, "depth", PETSC_REAL, &value));
    value = wave.wave_period();
    PetscCall(PetscViewerHDF5WriteAttribute(viewer_, nullptr, "period", PETSC_REAL, &value));
    value = wave.wavelength();
    PetscCall(PetscViewerHDF5WriteAttribute(viewer_, nullptr, "wavelength", PETSC_REAL, &value));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DMStagHDF5Writer::write_solution(Vec sol, PetscReal t, const char *info)
{
    PetscFunctionBeginUser;

    PetscCall(PetscObjectSetName(reinterpret_cast<PetscObject>(sol), "sol"));
    PetscCall(VecSet(sol_time_, t));
    PetscCall(PetscViewerHDF5PushTimestepping(viewer_));
    PetscCall(PetscViewerHDF5SetTimestep(viewer_, sol_timestep_));
    PetscCall(VecView(sol, viewer_));
    PetscCall(VecView(sol_time_, viewer_));
    PetscCall(PetscViewerHDF5PopTimestepping(viewer_));
    if (sol_timestep_ == 0)
    {
        PetscCall(PetscViewerHDF5WriteObjectAttribute(viewer_, reinterpret_cast<PetscObject>(sol),
                                                      "layout", PETSC_STRING,
                                                      "dmstag-entry-aligned"));
    }
    if (info != nullptr)
    {
        PetscCall(PetscViewerHDF5WriteObjectAttribute(viewer_, reinterpret_cast<PetscObject>(sol),
                                                      "info", PETSC_STRING, info));
    }
    ++sol_timestep_;

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DMStagHDF5Writer::write_surface_elevation(Vec eta, PetscReal t, const char *info)
{
    PetscFunctionBeginUser;

    PetscCall(VecScatterBegin(eta_.scatter, eta, eta_.values, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(eta_.scatter, eta, eta_.values, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecSet(eta_.time, t));

    PetscCall(PetscViewerHDF5PushTimestepping(viewer_));
    PetscCall(PetscViewerHDF5SetTimestep(viewer_, eta_.timestep));
    PetscCall(VecView(eta_.values, viewer_));
    PetscCall(VecView(eta_.time, viewer_));
    PetscCall(PetscViewerHDF5PopTimestepping(viewer_));

    if (eta_.timestep == 0)
    {
        PetscCall(PetscViewerHDF5WriteObjectAttribute(
            viewer_, reinterpret_cast<PetscObject>(eta_.values), "location", PETSC_STRING, "up"));
    }
    if (info != nullptr)
    {
        PetscCall(PetscViewerHDF5WriteObjectAttribute(
            viewer_, reinterpret_cast<PetscObject>(eta_.values), "info", PETSC_STRING, info));
    }
    ++eta_.timestep;

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DMStagHDF5Writer::destroy()
{
    PetscFunctionBeginUser;

    PetscCall(VecScatterDestroy(&eta_.scatter));
    PetscCall(VecDestroy(&eta_.time));
    PetscCall(VecDestroy(&eta_.values));
    PetscCall(VecDestroy(&sol_time_));
    PetscCall(PetscViewerDestroy(&viewer_));

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace wavein
