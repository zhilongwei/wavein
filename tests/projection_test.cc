#include "wavein/projection.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <petscdmstag.h>
#include <petscsys.h>

namespace
{

constexpr PetscReal xmin = 0.0;
constexpr PetscReal xmax = 1.0;
constexpr PetscReal zmin = 0.0;
constexpr PetscReal zmax = 1.0;
constexpr PetscInt nx = 10;
constexpr PetscInt nz = 15;
constexpr PetscReal dt = 0.01;
constexpr unsigned long random_seed = 20250804UL;
constexpr PetscReal relative_divergence_tolerance = 1.0e-10;

struct ProjectionResult
{
        PetscReal initial_divergence_norm;
        PetscReal projected_divergence_norm;
        PetscReal top_pressure_norm;
};

ProjectionResult project_random_velocity(DMBoundaryType x_boundary)
{
    const MPI_Comm comm = PETSC_COMM_WORLD;

    DM dm = nullptr;
    PetscCallAbort(comm, DMStagCreate2d(comm, x_boundary, DM_BOUNDARY_NONE, nx, nz, PETSC_DECIDE,
                                        PETSC_DECIDE, 0, 1, 1, DMSTAG_STENCIL_BOX, 1, nullptr,
                                        nullptr, &dm));
    PetscCallAbort(comm, DMSetUp(dm));

    const PetscReal dx = (xmax - xmin) / static_cast<PetscReal>(nx);
    const PetscReal dz = (zmax - zmin) / static_cast<PetscReal>(nz);

    PetscReal initial_divergence_norm = 0.0;
    PetscReal projected_divergence_norm = 0.0;
    PetscReal top_pressure_norm = 0.0;

    PetscCallAbort(comm, PetscOptionsSetValue(nullptr, "-pc_type", PCJACOBI));
    PetscCallAbort(comm, PetscOptionsSetValue(nullptr, "-ksp_rtol", "1.0e-12"));

    {
        wavein::Projection projection(comm, dm, dx, dz);

        Vec dvel = nullptr;
        Vec sol = nullptr;
        Vec ptop = nullptr;
        PetscRandom random = nullptr;

        PetscCallAbort(comm, DMCreateGlobalVector(dm, &dvel));
        PetscCallAbort(comm, DMCreateGlobalVector(dm, &sol));
        PetscCallAbort(comm, DMCreateGlobalVector(dm, &ptop));

        PetscCallAbort(comm, PetscRandomCreate(comm, &random));
        PetscCallAbort(comm, PetscRandomSetSeed(random, random_seed));
        PetscCallAbort(comm, PetscRandomSeed(random));
        PetscCallAbort(comm, PetscRandomSetInterval(random, -1.0, 1.0));
        PetscCallAbort(comm, VecSetRandom(sol, random));
        PetscCallAbort(comm, VecSetRandom(ptop, random));

        PetscCallAbort(comm, projection.divergence_norm(sol, &initial_divergence_norm));
        PetscCallAbort(comm, VecNorm(ptop, NORM_2, &top_pressure_norm));

        PetscCallAbort(comm, projection.project(dvel, sol, ptop, dt));
        PetscCallAbort(comm, VecAXPY(sol, 1.0, dvel));
        PetscCallAbort(comm, projection.divergence_norm(sol, &projected_divergence_norm));

        PetscCallAbort(comm, PetscRandomDestroy(&random));
        PetscCallAbort(comm, VecDestroy(&dvel));
        PetscCallAbort(comm, VecDestroy(&sol));
        PetscCallAbort(comm, VecDestroy(&ptop));
    }

    PetscCallAbort(comm, PetscOptionsClearValue(nullptr, "-pc_type"));
    PetscCallAbort(comm, PetscOptionsClearValue(nullptr, "-ksp_rtol"));
    PetscCallAbort(comm, DMDestroy(&dm));

    return {initial_divergence_norm, projected_divergence_norm, top_pressure_norm};
}

} // namespace

TEST_CASE("Projection produces divergence-free velocity with periodic x boundary",
          "[projection][periodic]")
{
    const ProjectionResult result = project_random_velocity(DM_BOUNDARY_PERIODIC);

    REQUIRE(result.initial_divergence_norm > 0.0);
    REQUIRE(result.top_pressure_norm > 0.0);
    REQUIRE(
        result.projected_divergence_norm ==
        Catch::Approx(0.0).margin(relative_divergence_tolerance * result.initial_divergence_norm));
}

TEST_CASE("Projection produces divergence-free velocity with non-periodic x boundary",
          "[projection][none]")
{
    const ProjectionResult result = project_random_velocity(DM_BOUNDARY_NONE);

    REQUIRE(result.initial_divergence_norm > 0.0);
    REQUIRE(result.top_pressure_norm > 0.0);
    REQUIRE(
        result.projected_divergence_norm ==
        Catch::Approx(0.0).margin(relative_divergence_tolerance * result.initial_divergence_norm));
}
