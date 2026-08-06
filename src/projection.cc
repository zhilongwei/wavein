#include "wavein/projection.h"

namespace wavein
{

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
Projection::Projection(MPI_Comm comm, DM dm, PetscReal dx, PetscReal dz, PetscReal rhow)
    : comm_(comm), dm_(dm), dx_(dx), dz_(dz), rhow_(rhow)
{
    PetscFunctionBeginUser;

    DMBoundaryType xboundary, yboundary;
    PetscCallAbort(comm_, DMStagGetBoundaryTypes(dm_, &xboundary, &yboundary, nullptr));
    PetscCheckAbort(xboundary == DM_BOUNDARY_PERIODIC || xboundary == DM_BOUNDARY_NONE, comm_,
                    PETSC_ERR_ARG_WRONG,
                    "DMStag must have none or periodic boundary conditions in the x-direction");
    PetscCheckAbort(yboundary == DM_BOUNDARY_NONE, comm_, PETSC_ERR_ARG_WRONG,
                    "DMStag must have none boundary conditions in the z-direction");

    // Create pressure grid
    PetscCallAbort(comm_, DMStagCreateCompatibleDMStag(dm_, 0, 0, 1, 0, &dm_pressure_));

    // Create the matrices
    PetscCallAbort(comm_, DMCreateMatrix(dm_pressure_, &mat_a_));
    PetscCallAbort(comm_, DMCreateMatrix(dm_, &mat_div_vel_));
    PetscCallAbort(comm_, DMCreateMatrix(dm_, &mat_grad_pressure_));
    PetscCallAbort(comm_, DMCreateMatrix(dm_, &mat_d2p_dz2_top_Dirichlet_bc_));
    PetscCallAbort(comm_, DMCreateMatrix(dm_, &mat_dp_dz_top_Dirichlet_bc_));
    PetscCallAbort(comm_, DMCreateMatrix(dm_, &mat_migrate_top_p_up2element_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm_, &dvel_));

    // Create a mask vector for the velocity field
    PetscCallAbort(comm_, DMCreateGlobalVector(dm_, &velocity_mask_));

    const PetscReal x_scale = 1.0 / dx_ / dx_;
    const PetscReal z_scale = 1.0 / dz_ / dz_;

    PetscInt ndim[2];
    PetscInt startx, starty, nx, ny, n_extra[2], ex, ey;
    PetscInt ncol;
    DMStagStencil row, col[5];
    PetscReal val[5];

    PetscCallAbort(comm_, DMStagGetGlobalSizes(dm_, &ndim[0], &ndim[1], nullptr));
    PetscCallAbort(comm_, DMStagGetCorners(dm_, &startx, &starty, nullptr, &nx, &ny, nullptr,
                                           &n_extra[0], &n_extra[1], nullptr));

    // Assemble the matrix for the Poisson equation
    for (ey = starty; ey != starty + ny; ++ey)
    {
        for (ex = startx; ex != startx + nx; ++ex)
        {
            if (xboundary == DM_BOUNDARY_NONE)
            {
                ncol = 0;

                row.i = ex;
                row.j = ey;
                row.loc = DMSTAG_ELEMENT;
                row.c = 0;

                col[ncol].i = ex;
                col[ncol].j = ey;
                col[ncol].loc = DMSTAG_ELEMENT;
                col[ncol].c = 0;

                val[ncol] = -2.0 * x_scale - 2.0 * z_scale;
                ++ncol;

                if (ex == 0)
                {
                    val[0] += x_scale;
                }
                else
                {
                    col[ncol].i = ex - 1;
                    col[ncol].j = ey;
                    col[ncol].loc = DMSTAG_ELEMENT;
                    col[ncol].c = 0;

                    val[ncol] = x_scale;
                    ++ncol;
                }

                if (ex == ndim[0] - 1)
                {
                    val[0] += x_scale;
                }
                else
                {
                    col[ncol].i = ex + 1;
                    col[ncol].j = ey;
                    col[ncol].loc = DMSTAG_ELEMENT;
                    col[ncol].c = 0;

                    val[ncol] = x_scale;
                    ++ncol;
                }

                if (ey == 0)
                {
                    val[0] += z_scale;
                }
                else
                {
                    col[ncol].i = ex;
                    col[ncol].j = ey - 1;
                    col[ncol].loc = DMSTAG_ELEMENT;
                    col[ncol].c = 0;

                    val[ncol] = z_scale;
                    ++ncol;
                }

                if (ey == ndim[1] - 1)
                {
                    val[0] -= z_scale;
                }
                else
                {
                    col[ncol].i = ex;
                    col[ncol].j = ey + 1;
                    col[ncol].loc = DMSTAG_ELEMENT;
                    col[ncol].c = 0;

                    val[ncol] = z_scale;
                    ++ncol;
                }

                PetscCallAbort(comm_, DMStagMatSetValuesStencil(dm_pressure_, mat_a_, 1, &row, ncol,
                                                                col, val, INSERT_VALUES));
            }

            else
            {
                ncol = 0;

                row.i = ex;
                row.j = ey;
                row.loc = DMSTAG_ELEMENT;
                row.c = 0;

                col[ncol].i = ex;
                col[ncol].j = ey;
                col[ncol].loc = DMSTAG_ELEMENT;
                col[ncol].c = 0;

                val[ncol] = -2.0 * x_scale - 2.0 * z_scale;
                ++ncol;

                col[ncol].i = ex - 1;
                col[ncol].j = ey;
                col[ncol].loc = DMSTAG_ELEMENT;
                col[ncol].c = 0;
                val[ncol] = x_scale;
                ++ncol;

                col[ncol].i = ex + 1;
                col[ncol].j = ey;
                col[ncol].loc = DMSTAG_ELEMENT;
                col[ncol].c = 0;
                val[ncol] = x_scale;
                ++ncol;

                if (ey == 0)
                {
                    val[0] += z_scale;
                }
                else
                {
                    col[ncol].i = ex;
                    col[ncol].j = ey - 1;
                    col[ncol].loc = DMSTAG_ELEMENT;
                    col[ncol].c = 0;

                    val[ncol] = z_scale;
                    ++ncol;
                }

                if (ey == ndim[1] - 1)
                {
                    val[0] -= z_scale;
                }
                else
                {
                    col[ncol].i = ex;
                    col[ncol].j = ey + 1;
                    col[ncol].loc = DMSTAG_ELEMENT;
                    col[ncol].c = 0;

                    val[ncol] = z_scale;
                    ++ncol;
                }

                PetscCallAbort(comm_, DMStagMatSetValuesStencil(dm_pressure_, mat_a_, 1, &row, ncol,
                                                                col, val, INSERT_VALUES));
            }
        }
    }

    PetscCallAbort(comm_, MatAssemblyBegin(mat_a_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm_, MatAssemblyEnd(mat_a_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm_, MatSetOption(mat_a_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));
    PetscCallAbort(comm_, MatSetOption(mat_a_, MAT_SPD, PETSC_TRUE));

    // Assemble the matrix for the divergence operator
    for (ey = starty; ey != starty + ny; ++ey)
    {
        for (ex = startx; ex != startx + nx; ++ex)
        {
            row.i = ex;
            row.j = ey;
            row.loc = DMSTAG_ELEMENT;
            row.c = 0;

            col[0].i = ex;
            col[0].j = ey;
            col[0].loc = DMSTAG_LEFT;
            col[0].c = 0;
            val[0] = -1.0 / dx_;

            col[1].i = ex;
            col[1].j = ey;
            col[1].loc = DMSTAG_RIGHT;
            col[1].c = 0;
            val[1] = 1.0 / dx_;

            col[2].i = ex;
            col[2].j = ey;
            col[2].loc = DMSTAG_DOWN;
            col[2].c = 0;
            val[2] = -1.0 / dz_;

            col[3].i = ex;
            col[3].j = ey;
            col[3].loc = DMSTAG_UP;
            col[3].c = 0;
            val[3] = 1.0 / dz_;

            PetscCallAbort(comm_, DMStagMatSetValuesStencil(dm_, mat_div_vel_, 1, &row, 4, col, val,
                                                            INSERT_VALUES));
        }
    }
    PetscCallAbort(comm_, MatAssemblyBegin(mat_div_vel_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm_, MatAssemblyEnd(mat_div_vel_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm_, MatSetOption(mat_div_vel_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));

    // Assemble the matrix for the gradient operator for interior cells
    if (xboundary == DM_BOUNDARY_NONE)
    {
        for (ey = starty; ey != starty + ny; ++ey)
        {
            for (ex = startx; ex != startx + nx; ++ex)
            {
                if (ex != 0)
                {
                    row.i = ex;
                    row.j = ey;
                    row.loc = DMSTAG_LEFT;
                    row.c = 0;

                    col[0].i = ex - 1;
                    col[0].j = ey;
                    col[0].loc = DMSTAG_ELEMENT;
                    col[0].c = 0;
                    val[0] = -1.0 / dx_;

                    col[1].i = ex;
                    col[1].j = ey;
                    col[1].loc = DMSTAG_ELEMENT;
                    col[1].c = 0;
                    val[1] = 1.0 / dx_;

                    PetscCallAbort(comm_,
                                   DMStagMatSetValuesStencil(dm_, mat_grad_pressure_, 1, &row, 2,
                                                             col, val, INSERT_VALUES));
                }

                if (ey != 0)
                {
                    row.i = ex;
                    row.j = ey;
                    row.loc = DMSTAG_DOWN;
                    row.c = 0;

                    col[0].i = ex;
                    col[0].j = ey - 1;
                    col[0].loc = DMSTAG_ELEMENT;
                    col[0].c = 0;
                    val[0] = -1.0 / dz_;

                    col[1].i = ex;
                    col[1].j = ey;
                    col[1].loc = DMSTAG_ELEMENT;
                    col[1].c = 0;
                    val[1] = 1.0 / dz_;

                    PetscCallAbort(comm_,
                                   DMStagMatSetValuesStencil(dm_, mat_grad_pressure_, 1, &row, 2,
                                                             col, val, INSERT_VALUES));

                    if (ey == ndim[1] - 1)
                    {
                        row.i = ex;
                        row.j = ey;
                        row.loc = DMSTAG_UP;
                        row.c = 0;

                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].loc = DMSTAG_ELEMENT;
                        col[0].c = 0;
                        val[0] = -2.0 / dz_;

                        PetscCallAbort(comm_,
                                       DMStagMatSetValuesStencil(dm_, mat_grad_pressure_, 1, &row,
                                                                 1, col, val, INSERT_VALUES));
                    }
                }
            }
        }
    }
    else // periodic boundary conditions in the x-direction
    {
        for (ey = starty; ey != starty + ny + n_extra[1]; ++ey)
        {
            for (ex = startx; ex != startx + nx + n_extra[0]; ++ex)
            {
                row.i = ex;
                row.j = ey;
                row.loc = DMSTAG_LEFT;
                row.c = 0;

                col[0].i = ex - 1;
                col[0].j = ey;
                col[0].loc = DMSTAG_ELEMENT;
                col[0].c = 0;
                val[0] = -1.0 / dx_;

                col[1].i = ex;
                col[1].j = ey;
                col[1].loc = DMSTAG_ELEMENT;
                col[1].c = 0;
                val[1] = 1.0 / dx_;

                PetscCallAbort(comm_, DMStagMatSetValuesStencil(dm_, mat_grad_pressure_, 1, &row, 2,
                                                                col, val, INSERT_VALUES));

                if (ey != 0)
                {
                    row.i = ex;
                    row.j = ey;
                    row.loc = DMSTAG_DOWN;
                    row.c = 0;

                    col[0].i = ex;
                    col[0].j = ey - 1;
                    col[0].loc = DMSTAG_ELEMENT;
                    col[0].c = 0;
                    val[0] = -1.0 / dz_;

                    col[1].i = ex;
                    col[1].j = ey;
                    col[1].loc = DMSTAG_ELEMENT;
                    col[1].c = 0;
                    val[1] = 1.0 / dz_;

                    PetscCallAbort(comm_,
                                   DMStagMatSetValuesStencil(dm_, mat_grad_pressure_, 1, &row, 2,
                                                             col, val, INSERT_VALUES));

                    if (ey == ndim[1])
                    {
                        row.i = ex;
                        row.j = ey;
                        row.loc = DMSTAG_DOWN;
                        row.c = 0;

                        col[0].i = ex;
                        col[0].j = ey - 1;
                        col[0].loc = DMSTAG_ELEMENT;
                        col[0].c = 0;
                        val[0] = -2.0 / dz_;

                        PetscCallAbort(comm_,
                                       DMStagMatSetValuesStencil(dm_, mat_grad_pressure_, 1, &row,
                                                                 1, col, val, INSERT_VALUES));
                    }
                }
            }
        }
    }
    PetscCallAbort(comm_, MatAssemblyBegin(mat_grad_pressure_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm_, MatAssemblyEnd(mat_grad_pressure_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm_, MatSetOption(mat_grad_pressure_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));

    // Assemble the matrix for the top Dirichlet boundary condition for d2p/dz2
    if (xboundary == DM_BOUNDARY_NONE)
    {
        for (ey = starty; ey != starty + ny; ++ey)
        {
            for (ex = startx; ex != startx + nx; ++ex)
            {
                if (ey == ndim[1] - 1)
                {
                    ncol = 0;

                    row.i = ex;
                    row.j = ey;
                    row.loc = DMSTAG_ELEMENT;
                    row.c = 0;

                    col[ncol].i = ex;
                    col[ncol].j = ey;
                    col[ncol].loc = DMSTAG_ELEMENT;
                    col[ncol].c = 0;
                    val[ncol] = -1.0 / 2.0 / dx_ / dx_ - 2.0 / dz_ / dz_;
                    ++ncol;

                    if (ex == 0)
                    {
                        val[0] += 0.25 / dx_ / dx_;
                    }
                    else
                    {
                        col[ncol].i = ex - 1;
                        col[ncol].j = ey;
                        col[ncol].loc = DMSTAG_ELEMENT;
                        col[ncol].c = 0;
                        val[ncol] = 1.0 / 4.0 / dx_ / dx_;
                        ++ncol;
                    }

                    if (ex == ndim[0] - 1)
                    {
                        val[0] += 0.25 / dx_ / dx_;
                    }
                    else
                    {
                        col[ncol].i = ex + 1;
                        col[ncol].j = ey;
                        col[ncol].loc = DMSTAG_ELEMENT;
                        col[ncol].c = 0;
                        val[ncol] = 1.0 / 4.0 / dx_ / dx_;
                        ++ncol;
                    }

                    PetscCallAbort(comm_,
                                   DMStagMatSetValuesStencil(dm_, mat_d2p_dz2_top_Dirichlet_bc_, 1,
                                                             &row, ncol, col, val, INSERT_VALUES));
                }
            }
        }
    }
    else // periodic boundary conditions in the x-direction
    {
        for (ey = starty; ey != starty + ny + n_extra[1]; ++ey)
        {
            for (ex = startx; ex != startx + nx + n_extra[0]; ++ex)
            {
                if (ey == ndim[1] - 1)
                {
                    row.i = ex;
                    row.j = ey;
                    row.loc = DMSTAG_ELEMENT;
                    row.c = 0;

                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].loc = DMSTAG_ELEMENT;
                    col[0].c = 0;
                    val[0] = -1.0 / 2.0 / dx_ / dx_ - 2.0 / dz_ / dz_;

                    col[1].i = ex + 1;
                    col[1].j = ey;
                    col[1].loc = DMSTAG_ELEMENT;
                    col[1].c = 0;
                    val[1] = 1.0 / 4.0 / dx_ / dx_;

                    col[2].i = ex - 1;
                    col[2].j = ey;
                    col[2].loc = DMSTAG_ELEMENT;
                    col[2].c = 0;
                    val[2] = 1.0 / 4.0 / dx_ / dx_;

                    PetscCallAbort(comm_,
                                   DMStagMatSetValuesStencil(dm_, mat_d2p_dz2_top_Dirichlet_bc_, 1,
                                                             &row, 3, col, val, INSERT_VALUES));
                }
            }
        }
    }
    PetscCallAbort(comm_, MatAssemblyBegin(mat_d2p_dz2_top_Dirichlet_bc_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm_, MatAssemblyEnd(mat_d2p_dz2_top_Dirichlet_bc_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(
        comm_, MatSetOption(mat_d2p_dz2_top_Dirichlet_bc_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));

    // Assemble the matrix for the top Dirichlet boundary condition for dp/dz
    if (xboundary == DM_BOUNDARY_NONE)
    {
        for (ey = starty; ey != starty + ny; ++ey)
        {
            for (ex = startx; ex != startx + nx; ++ex)
            {
                if (ey == ndim[1] - 1)
                {
                    ncol = 0;

                    row.i = ex;
                    row.j = ey;
                    row.loc = DMSTAG_ELEMENT;
                    row.c = 0;

                    col[ncol].i = ex;
                    col[ncol].j = ey;
                    col[ncol].loc = DMSTAG_ELEMENT;
                    col[ncol].c = 0;
                    val[ncol] = 1.0 / 2.0 * dz_ / dx_ / dx_ + 2.0 / dz_;
                    ++ncol;

                    if (ex == 0)
                    {
                        val[0] -= 1.0 / 4.0 * dz_ / dx_ / dx_;
                    }
                    else
                    {
                        col[ncol].i = ex - 1;
                        col[ncol].j = ey;
                        col[ncol].loc = DMSTAG_ELEMENT;
                        col[ncol].c = 0;
                        val[ncol] = -1.0 / 4.0 * dz_ / dx_ / dx_;
                        ++ncol;
                    }

                    if (ex == ndim[0] - 1)
                    {
                        val[0] -= 1.0 / 4.0 * dz_ / dx_ / dx_;
                    }
                    else
                    {
                        col[ncol].i = ex + 1;
                        col[ncol].j = ey;
                        col[ncol].loc = DMSTAG_ELEMENT;
                        col[ncol].c = 0;
                        val[ncol] = -1.0 / 4.0 * dz_ / dx_ / dx_;
                        ++ncol;
                    }
                    PetscCallAbort(comm_,
                                   DMStagMatSetValuesStencil(dm_, mat_dp_dz_top_Dirichlet_bc_, 1,
                                                             &row, ncol, col, val, INSERT_VALUES));
                }
            }
        }
    }
    else // periodic boundary conditions in the x-direction
    {
        for (ey = starty; ey != starty + ny + n_extra[1]; ++ey)
        {
            for (ex = startx; ex != startx + nx + n_extra[0]; ++ex)
            {
                if (ey == ndim[1] - 1)
                {
                    row.i = ex;
                    row.j = ey;
                    row.loc = DMSTAG_ELEMENT;
                    row.c = 0;

                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].loc = DMSTAG_ELEMENT;
                    col[0].c = 0;
                    val[0] = 1.0 / 2.0 * dz_ / dx_ / dx_ + 2.0 / dz_;

                    col[1].i = ex + 1;
                    col[1].j = ey;
                    col[1].loc = DMSTAG_ELEMENT;
                    col[1].c = 0;
                    val[1] = -1.0 / 4.0 * dz_ / dx_ / dx_;

                    col[2].i = ex - 1;
                    col[2].j = ey;
                    col[2].loc = DMSTAG_ELEMENT;
                    col[2].c = 0;
                    val[2] = -1.0 / 4.0 * dz_ / dx_ / dx_;

                    PetscCallAbort(comm_,
                                   DMStagMatSetValuesStencil(dm_, mat_dp_dz_top_Dirichlet_bc_, 1,
                                                             &row, 3, col, val, INSERT_VALUES));
                }
            }
        }
    }
    PetscCallAbort(comm_, MatAssemblyBegin(mat_dp_dz_top_Dirichlet_bc_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm_, MatAssemblyEnd(mat_dp_dz_top_Dirichlet_bc_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm_,
                   MatSetOption(mat_dp_dz_top_Dirichlet_bc_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));

    // Assemble the matrix for migrating the pressure value from DMSTAG_UP to DMSTAG_ELEMENT
    for (ey = starty; ey != starty + ny; ++ey)
    {
        for (ex = startx; ex != startx + nx; ++ex)
        {
            if (ey == ndim[1] - 1)
            {
                row.i = ex;
                row.j = ey;
                row.loc = DMSTAG_ELEMENT;
                row.c = 0;

                col[0].i = ex;
                col[0].j = ey;
                col[0].loc = DMSTAG_UP;
                col[0].c = 0;
                val[0] = 1.0;

                PetscCallAbort(comm_,
                               DMStagMatSetValuesStencil(dm_, mat_migrate_top_p_up2element_, 1,
                                                         &row, 1, col, val, INSERT_VALUES));
            }
        }
    }
    PetscCallAbort(comm_, MatAssemblyBegin(mat_migrate_top_p_up2element_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm_, MatAssemblyEnd(mat_migrate_top_p_up2element_, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(
        comm_, MatSetOption(mat_migrate_top_p_up2element_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));

    PC prec_ = nullptr;
    PetscCallAbort(comm_, KSPCreate(comm_, &ksp_));
    PetscCallAbort(comm_, KSPSetType(ksp_, KSPGMRES));
    PetscCallAbort(comm_, KSPGetPC(ksp_, &prec_));
    PetscCallAbort(comm_, PCSetType(prec_, PCLU));
    PetscCallAbort(comm_, PCFactorSetMatSolverType(prec_, MATSOLVERMUMPS));
    PetscCallAbort(comm_, KSPSetOperators(ksp_, mat_a_, mat_a_));
    PetscCallAbort(comm_, KSPSetFromOptions(ksp_));
    PetscCallAbort(comm_, KSPSetErrorIfNotConverged(ksp_, PETSC_TRUE));

    // Create the auxiliary vectors for the Poisson solver
    PetscCallAbort(comm_, DMCreateGlobalVector(dm_, &rhs_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm_, &pressure_top_row_element_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm_pressure_, &rhs_pressure_grid_));
    PetscCallAbort(comm_, DMCreateGlobalVector(dm_pressure_, &sol_pressure_grid_));

    // Set the velocity mask to zero for all pressure DOFs
    Vec velocity_mask_local = nullptr;
    PetscCallAbort(comm_, DMGetLocalVector(dm_, &velocity_mask_local));
    PetscCallAbort(comm_, VecSet(velocity_mask_local, 1.0));

    PetscInt ip;
    PetscCallAbort(comm, DMStagGetLocationSlot(dm_, DMSTAG_ELEMENT, 0, &ip));

    PetscReal ***c_arr_velocity_mask_local = nullptr;
    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm_, DMStagVecGetArray(dm_, velocity_mask_local, &c_arr_velocity_mask_local));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)

    for (ey = starty; ey != starty + ny; ++ey)
    {
        for (ex = startx; ex != startx + nx; ++ex)
        {
            c_arr_velocity_mask_local[ey][ex][ip] = 0.0;
        }
    }

    // NOLINTBEGIN(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm_,
                   DMStagVecRestoreArray(dm_, velocity_mask_local, &c_arr_velocity_mask_local));
    // NOLINTEND(bugprone-multi-level-implicit-pointer-conversion)
    PetscCallAbort(comm_, DMLocalToGlobal(dm_, velocity_mask_local, INSERT_VALUES, velocity_mask_));
    PetscCallAbort(comm_, DMRestoreLocalVector(dm_, &velocity_mask_local));

    PetscFunctionReturnVoid();
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
PetscErrorCode Projection::project(Vec sol, Vec ptop, PetscReal dt) noexcept
{
    PetscFunctionBeginUser;

    // rhs on collected grid
    PetscCall(MatMult(mat_div_vel_, sol, rhs_));
    PetscCall(VecScale(rhs_, rhow_ / dt));

    // Apply Dirichlet boundary condition on the linearized water surface
    PetscCall(MatMult(mat_migrate_top_p_up2element_, ptop, pressure_top_row_element_));
    PetscCall(MatMultAdd(mat_d2p_dz2_top_Dirichlet_bc_, pressure_top_row_element_, rhs_, rhs_));

    // Migrate rhs to the pressure grid
    PetscCall(DMStagMigrateVec(dm_, rhs_, dm_pressure_, rhs_pressure_grid_));

    // Solve the Poisson equation
    PetscCall(KSPSolve(ksp_, rhs_pressure_grid_, sol_pressure_grid_));

    // Migrate the solution back to the velocity grid
    PetscCall(DMStagMigrateVec(dm_pressure_, sol_pressure_grid_, dm_, rhs_));

    // dvel now stores the pressure gradient
    PetscCall(MatMult(mat_grad_pressure_, rhs_, dvel_));

    // Store solved pressure in sol before it is reused.
    PetscCall(VecPointwiseMult(sol, sol, velocity_mask_));
    PetscCall(VecAXPY(sol, 1.0, rhs_));

    // dp/dz at the top boundary
    PetscCall(MatMult(mat_dp_dz_top_Dirichlet_bc_, pressure_top_row_element_, rhs_));
    PetscCall(MatMultTransposeAdd(mat_migrate_top_p_up2element_, rhs_, dvel_, dvel_));

    // Change of velocity due to pressure gradient
    PetscCall(VecAXPY(sol, -dt / rhow_, dvel_));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode Projection::divergence(Vec sol, Vec div) const noexcept
{
    PetscFunctionBeginUser;

    PetscCall(MatMult(mat_div_vel_, sol, div));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode Projection::divergence_norm(Vec sol, PetscReal *norm) noexcept
{
    PetscFunctionBeginUser;

    PetscCall(MatMult(mat_div_vel_, sol, rhs_));
    PetscCall(VecNorm(rhs_, NORM_2, norm));

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace wavein