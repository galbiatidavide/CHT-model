#ifndef IMPLICIT_POTENTIAL_SOLVER_HPP
#define IMPLICIT_POTENTIAL_SOLVER_HPP

#include "parameters.hpp"

using namespace problem_setting;

class ImplicitPotentialSolver {
public:
// Constructor
ImplicitPotentialSolver()
{
    
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            for(int k = 0; k < nz; k++){
                /*wx[i][j][k] = wx_ref(i + 0.5, j + 0.5, k + 0.5);
                wy[i][j][k] = wy_ref(i + 0.5, j + 0.5, k + 0.5);
                wz[i][j][k] = wz_ref(i + 0.5, j + 0.5, k + 0.5);*/
                continue;
            }
        }
    }

    create_grid();
    DMCreateMatrix(dmGrid_centered, &A);
    DMCreateGlobalVector(dmGrid_centered, &rhs_x);
    DMCreateGlobalVector(dmGrid_centered, &rhs_y);
    DMCreateGlobalVector(dmGrid_centered, &rhs_z);
    DMCreateGlobalVector(dmGrid_centered, &sol_x);
    DMCreateGlobalVector(dmGrid_centered, &sol_y);
    DMCreateGlobalVector(dmGrid_centered, &sol_z);
    DMCreateGlobalVector(dmGrid_centered, &U_x);
    DMCreateGlobalVector(dmGrid_centered, &U_y);
    DMCreateGlobalVector(dmGrid_centered, &U_z);
    DMCreateGlobalVector(dmGrid_centered, &W_x);
    DMCreateGlobalVector(dmGrid_centered, &W_y);
    DMCreateGlobalVector(dmGrid_centered, &W_z);

    assemble_lhs();
    build_petsc_vec(rhs_x, wx);
    build_petsc_vec(rhs_y, wy);
    build_petsc_vec(rhs_z, wz);

}

PetscErrorCode const solve_problem(){

    PetscFunctionBegin;
    VecScale(rhs_x, -1.0);
    VecScale(rhs_y, -1.0);
    VecScale(rhs_z, -1.0);

    solve_laplacian(rhs_x, sol_x);
    solve_laplacian(rhs_y, sol_y);
    solve_laplacian(rhs_z, sol_z);
    velocity_from_potential();

    PetscFunctionReturn(0);
}

PetscErrorCode const build_std_velocity(double (*ux)[ny][nz], double (*uy)[ny][nz], double (*uz)[ny][nz], int rank, int local_nx)
{
    PetscInt iu_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    Vec vecU_x_loc, vecU_y_loc, vecU_z_loc;
    PetscReal ****arrU_x, ****arrU_y, ****arrU_z;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_centered, &N[0], &N[1], &N[2]);

    DMStagGetLocationSlot(dmGrid_centered, ELEMENT, 0, &iu_element);  

    DMCreateLocalVector(dmGrid_centered, &vecU_x_loc);
    DMGlobalToLocalBegin(dmGrid_centered, U_x, INSERT_VALUES, vecU_x_loc);
    DMGlobalToLocalEnd(dmGrid_centered, U_x, INSERT_VALUES, vecU_x_loc);
    DMStagVecGetArrayRead(dmGrid_centered, vecU_x_loc, &arrU_x);

    DMCreateLocalVector(dmGrid_centered, &vecU_y_loc);
    DMGlobalToLocalBegin(dmGrid_centered, U_y, INSERT_VALUES, vecU_y_loc);
    DMGlobalToLocalEnd(dmGrid_centered, U_y, INSERT_VALUES, vecU_y_loc);
    DMStagVecGetArrayRead(dmGrid_centered, vecU_y_loc, &arrU_y);

    DMCreateLocalVector(dmGrid_centered, &vecU_z_loc);
    DMGlobalToLocalBegin(dmGrid_centered, U_z, INSERT_VALUES, vecU_z_loc);
    DMGlobalToLocalEnd(dmGrid_centered, U_z, INSERT_VALUES, vecU_z_loc);
    DMStagVecGetArrayRead(dmGrid_centered, vecU_z_loc, &arrU_z);
        
    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                ux[ex][ey][ez] = arrU_x[ez][ey][ex][iu_element];
                uy[ex][ey][ez] = arrU_y[ez][ey][ex][iu_element];
                uz[ex][ey][ez] = arrU_z[ez][ey][ex][iu_element];
            }
        }
    }

    DMStagVecRestoreArrayRead(dmGrid_centered, vecU_x_loc, &arrU_x);
    DMRestoreLocalVector(dmGrid_centered, &vecU_x_loc);
    DMStagVecRestoreArrayRead(dmGrid_centered, vecU_y_loc, &arrU_y);
    DMRestoreLocalVector(dmGrid_centered, &vecU_y_loc);
    DMStagVecRestoreArrayRead(dmGrid_centered, vecU_z_loc, &arrU_z);
    DMRestoreLocalVector(dmGrid_centered, &vecU_z_loc);

    MPI_Allgather(&ux[rank * local_nx][0][0], local_nx * ny * nz, MPI_DOUBLE, &ux[0][0][0], local_nx * ny * nz, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&uy[rank * local_nx][0][0], local_nx * ny * nz, MPI_DOUBLE, &uy[0][0][0], local_nx * ny * nz, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&uz[rank * local_nx][0][0], local_nx * ny * nz, MPI_DOUBLE, &uz[0][0][0], local_nx * ny * nz, MPI_DOUBLE, MPI_COMM_WORLD);
 
    PetscFunctionReturn(0);
}

~ImplicitPotentialSolver() {
    
    MatDestroy(&A);
    VecDestroy(&rhs_x);
    VecDestroy(&rhs_y);
    VecDestroy(&rhs_z);
    VecDestroy(&sol_x);
    VecDestroy(&sol_y);
    VecDestroy(&sol_z);
    VecDestroy(&U_x);
    VecDestroy(&U_y);
    VecDestroy(&U_z);
    VecDestroy(&W_x);
    VecDestroy(&W_y);
    VecDestroy(&W_z);
    DMDestroy(&dmGrid_centered);
    free(wx);
    free(wy);
    free(wz);
}

private:

DM dmGrid_centered;
Mat A;
Vec rhs_x, rhs_y, rhs_z, sol_x, sol_y, sol_z, U_x, U_y, U_z, W_x, W_y, W_z;

double (*wx)[ny][nz] = (double(*)[ny][nz]) malloc(nx * ny * nz * sizeof(double));
double (*wy)[ny][nz] = (double(*)[ny][nz]) malloc(nx * ny * nz * sizeof(double));
double (*wz)[ny][nz] = (double(*)[ny][nz]) malloc(nx * ny * nz * sizeof(double));

PetscErrorCode const create_grid()
{
    const PetscInt dof0 = 0;
    const PetscInt dof1 = 0;
    const PetscInt dof2 = 0;
    const PetscInt dof3 = 1;
    const PetscInt stencilWidth = 1;

    PetscFunctionBegin;

    DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, nx, ny, nz, PETSC_DECIDE, 1, 1, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, &dmGrid_centered);
    DMSetFromOptions(dmGrid_centered);
    DMSetUp(dmGrid_centered);
    DMStagSetUniformCoordinatesExplicit(dmGrid_centered, 0, nx, 0, ny, 0, nz);

    PetscFunctionReturn(0);
}

// Implicit stiffness matrix (Homogeneous Neumann BC)
PetscErrorCode const assemble_lhs() 
{
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_centered, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0;
    PetscReal const hy = 1.0;
    PetscReal const hz = 1.0;

    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ex == N[0] - 1) {

                    DMStagStencil row, col[6];
                    PetscReal valA[6];
                    PetscInt nEntries;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {                  
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx); 
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);                            
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        }
                    } else if (ez == 0) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez + 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                    } else {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                    }
                    DMStagMatSetValuesStencil(dmGrid_centered, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              
                } 
                else if (ex == 0) {
                    DMStagStencil row, col[6];
                    PetscReal valA[6];
                    PetscInt nEntries;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        }
                    } else if (ez == 0) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex + 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez + 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex + 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                    } else {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex + 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                    }
                    DMStagMatSetValuesStencil(dmGrid_centered, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              
                } else {
                    DMStagStencil row, col[7];
                    PetscReal valA[7];
                    PetscInt nEntries;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = ELEMENT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = ELEMENT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);
                        }
                    } else if (ez == 0) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez - 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                    } else {
                        nEntries = 7;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez - 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        col[6].i = ex;
                        col[6].j = ey;
                        col[6].k = ez + 1;
                        col[6].loc = ELEMENT;
                        col[6].c = 0;
                        valA[6] = 1.0 / (hz * hz);
                    }
                    DMStagMatSetValuesStencil(dmGrid_centered, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                }                
            }
        }
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    PetscFunctionReturn(0); 
}

/*
// Implicit stiffness matrix (Dirichlet boundary conditions)

PetscErrorCode const assemble_lhs() 
{
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_centered, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0;
    PetscReal const hy = 1.0;
    PetscReal const hz = 1.0;

    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ex == N[0] - 1) {

                    DMStagStencil row, col[6];
                    PetscReal valA[6];
                    PetscInt nEntries;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {                  
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -3.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx); 
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -3.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);                            
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -3.0 / (hx * hx) + -3.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -3.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -3.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -3.0 / (hx * hx) + -3.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        }
                    } else if (ez == 0) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -3.0 / (hx * hx) + -2.0 / (hy * hy) - 3.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez + 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -3.0 / (hx * hx) + -2.0 / (hy * hy) - 3.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                    } else {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -3.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                    }
                    DMStagMatSetValuesStencil(dmGrid_centered, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              
                } 
                else if (ex == 0) {
                    DMStagStencil row, col[6];
                    PetscReal valA[6];
                    PetscInt nEntries;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -3.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -3.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -3.0 / (hx * hx) + -3.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -3.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -3.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -3.0 / (hx * hx) + -3.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        }
                    } else if (ez == 0) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -3.0 / (hx * hx) + -2.0 / (hy * hy) - 3.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex + 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez + 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -3.0 / (hx * hx) + -2.0 / (hy * hy) - 3.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex + 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                    } else {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -3.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex + 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                    }
                    DMStagMatSetValuesStencil(dmGrid_centered, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              
                } else {
                    DMStagStencil row, col[7];
                    PetscReal valA[7];
                    PetscInt nEntries;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -3.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = ELEMENT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -3.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = ELEMENT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);
                        }
                    } else if (ez == 0) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 3.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 3.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez - 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                    } else {
                        nEntries = 7;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez - 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        col[6].i = ex;
                        col[6].j = ey;
                        col[6].k = ez + 1;
                        col[6].loc = ELEMENT;
                        col[6].c = 0;
                        valA[6] = 1.0 / (hz * hz);
                    }
                    DMStagMatSetValuesStencil(dmGrid_centered, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                }                
            }
        }
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    PetscFunctionReturn(0); 
}
*/

PetscErrorCode const build_petsc_vec(Vec const & rhs, double (*vec)[ny][nz])
{
    PetscInt iu_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    Vec vecOutLocal;
    PetscReal ****arrOut;  

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_centered, &N[0], &N[1], &N[2]);

    DMStagGetLocationSlot(dmGrid_centered, ELEMENT, 0, &iu_element);  

    DMGetLocalVector(dmGrid_centered, &vecOutLocal);
    DMStagVecGetArray(dmGrid_centered, vecOutLocal, &arrOut);

    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                arrOut[ez][ey][ex][iu_element] = vec[ex][ey][ez];
            }
        }
    }

    DMStagVecRestoreArray(dmGrid_centered, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid_centered, vecOutLocal, INSERT_VALUES, rhs);    
    DMRestoreLocalVector(dmGrid_centered, &vecOutLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const solve_laplacian(Vec const & rhs, Vec & sol)
{
    KSP ksp;
    PC  pc;

    PetscFunctionBegin;    

    AttachNullspace(dmGrid_centered, A);
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, sol);

    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);
    PetscInt iterations;
    KSPGetIterationNumber(ksp, &iterations);
    PetscReal residual_norm;
    KSPGetResidualNorm(ksp, &residual_norm);

    if (reason < 0) {
        PetscPrintf(PETSC_COMM_WORLD, "p-field KSP did not converge. Reason: %s\n", KSPConvergedReasons[reason]);
    } else {
        PetscPrintf(PETSC_COMM_WORLD, 
                    "p-field KSP converged in %d iterations with a final residual norm of %g. Reason: %s\n", 
                    iterations, residual_norm, KSPConvergedReasons[reason]);
    }

    KSPDestroy(&ksp);
    PetscFunctionReturn(0); 
}

PetscErrorCode AttachNullspace(DM dmSol, Mat A)
{
    DM           dmPressure;
    Vec          constantPressure, basis;
    PetscReal    nrm;
    MatNullSpace matNullSpace;

    PetscFunctionBeginUser;
    DMStagCreateCompatibleDMStag(dmSol, 0, 0, 1, 0, &dmPressure);
    DMGetGlobalVector(dmPressure, &constantPressure);
    VecSet(constantPressure, 1.0);
    VecNorm(constantPressure, NORM_2, &nrm);
    VecScale(constantPressure, 1.0 / nrm);
    DMCreateGlobalVector(dmSol, &basis);
    DMStagMigrateVec(dmPressure, constantPressure, dmSol, basis);
    MatNullSpaceCreate(PetscObjectComm((PetscObject)dmSol), PETSC_FALSE, 1, &basis, &matNullSpace);
    VecDestroy(&basis);
    VecDestroy(&constantPressure);
    MatSetNullSpace(A, matNullSpace);
    MatNullSpaceDestroy(&matNullSpace);
    PetscFunctionReturn(0);
}

PetscErrorCode const velocity_from_potential()
{
    PetscInt iu_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    Vec vecA_x, vecA_y, vecA_z, vecU, vecV, vecW;
    PetscReal ****arrA_x, ****arrA_y, ****arrA_z, ****arrU, ****arrV, ****arrW;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_centered, &N[0], &N[1], &N[2]);

    DMStagGetLocationSlot(dmGrid_centered, ELEMENT, 0, &iu_element);

    DMGetLocalVector(dmGrid_centered, &vecU);
    DMStagVecGetArray(dmGrid_centered, vecU, &arrU);
    DMGetLocalVector(dmGrid_centered, &vecV);
    DMStagVecGetArray(dmGrid_centered, vecV, &arrV);
    DMGetLocalVector(dmGrid_centered, &vecW);
    DMStagVecGetArray(dmGrid_centered, vecW, &arrW);

    DMCreateLocalVector(dmGrid_centered, &vecA_x);
    DMGlobalToLocalBegin(dmGrid_centered, sol_x, INSERT_VALUES, vecA_x);
    DMGlobalToLocalEnd(dmGrid_centered, sol_x, INSERT_VALUES, vecA_x);
    DMStagVecGetArrayRead(dmGrid_centered, vecA_x, &arrA_x);

    DMCreateLocalVector(dmGrid_centered, &vecA_y);
    DMGlobalToLocalBegin(dmGrid_centered, sol_y, INSERT_VALUES, vecA_y);
    DMGlobalToLocalEnd(dmGrid_centered, sol_y, INSERT_VALUES, vecA_y);
    DMStagVecGetArrayRead(dmGrid_centered, vecA_y, &arrA_y);

    DMCreateLocalVector(dmGrid_centered, &vecA_z);
    DMGlobalToLocalBegin(dmGrid_centered, sol_z, INSERT_VALUES, vecA_z);
    DMGlobalToLocalEnd(dmGrid_centered, sol_z, INSERT_VALUES, vecA_z);
    DMStagVecGetArrayRead(dmGrid_centered, vecA_z, &arrA_z);

    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if(ex != 0 and ex != N[0] - 1 and ey != 0 and ey != N[1] - 1 and ez != 0 and ez != N[2] - 1) {
                    arrU[ez][ey][ex][iu_element] = (arrA_z[ez][ey + 1][ex][iu_element] - arrA_z[ez][ey - 1][ex][iu_element])/2.0 - (arrA_y[ez + 1][ey][ex][iu_element] - arrA_y[ez - 1][ey][ex][iu_element])/2.0;
                    arrV[ez][ey][ex][iu_element] = (arrA_x[ez + 1][ey][ex][iu_element] - arrA_x[ez - 1][ey][ex][iu_element])/2.0 - (arrA_z[ez][ey][ex + 1][iu_element] - arrA_z[ez][ey][ex - 1][iu_element])/2.0;
                    arrW[ez][ey][ex][iu_element] = (arrA_y[ez][ey][ex + 1][iu_element] - arrA_y[ez][ey][ex - 1][iu_element])/2.0 - (arrA_x[ez][ey + 1][ex][iu_element] - arrA_x[ez][ey - 1][ex][iu_element])/2.0;

                }

                if(ey == 0){
                    if(ez == 0){
                        arrU[ez][ey][ex][iu_element] = (arrA_z[ez][ey + 1][ex][iu_element] - arrA_z[ez][ey][ex][iu_element]) - (arrA_y[ez + 1][ey][ex][iu_element] - arrA_y[ez][ey][ex][iu_element]);
                    }
                    else if(ez == N[2] - 1){
                        arrU[ez][ey][ex][iu_element] = (arrA_z[ez][ey + 1][ex][iu_element] - arrA_z[ez][ey][ex][iu_element]) - (arrA_y[ez][ey][ex][iu_element] - arrA_y[ez - 1][ey][ex][iu_element]);
                    }
                    else {
                        arrU[ez][ey][ex][iu_element] = (arrA_z[ez][ey + 1][ex][iu_element] - arrA_z[ez][ey][ex][iu_element]) - (arrA_y[ez + 1][ey][ex][iu_element] - arrA_y[ez - 1][ey][ex][iu_element])/2;
                    }
                }

                if(ey == N[1] - 1){
                    if(ez == 0){
                        arrU[ez][ey][ex][iu_element] = (arrA_z[ez][ey][ex][iu_element] - arrA_z[ez][ey - 1][ex][iu_element]) - (arrA_y[ez + 1][ey][ex][iu_element] - arrA_y[ez][ey][ex][iu_element]);
                    }
                    else if(ez == N[2] - 1){
                        arrU[ez][ey][ex][iu_element] = (arrA_z[ez][ey][ex][iu_element] - arrA_z[ez][ey - 1][ex][iu_element]) - (arrA_y[ez][ey][ex][iu_element] - arrA_y[ez - 1][ey][ex][iu_element]);
                    }
                    else {
                        arrU[ez][ey][ex][iu_element] = (arrA_z[ez][ey][ex][iu_element] - arrA_z[ez][ey - 1][ex][iu_element]) - (arrA_y[ez + 1][ey][ex][iu_element] - arrA_y[ez - 1][ey][ex][iu_element])/2;
                    }
                }

                if(ez == 0){
                    if(ex == 0){
                        arrV[ez][ey][ex][iu_element] = (arrA_x[ez + 1][ey][ex][iu_element] - arrA_x[ez][ey][ex][iu_element]) - (arrA_z[ez][ey][ex + 1][iu_element] - arrA_z[ez][ey][ex][iu_element]);
                    }
                    else if (ex == N[0] - 1){
                        arrV[ez][ey][ex][iu_element] = (arrA_x[ez + 1][ey][ex][iu_element] - arrA_x[ez][ey][ex][iu_element]) - (arrA_z[ez][ey][ex][iu_element] - arrA_z[ez][ey][ex - 1][iu_element]);
                    }
                    else {
                        arrV[ez][ey][ex][iu_element] = (arrA_x[ez + 1][ey][ex][iu_element] - arrA_x[ez][ey][ex][iu_element]) - (arrA_z[ez][ey][ex + 1][iu_element] - arrA_z[ez][ey][ex - 1][iu_element])/2.0;
                    }
                }

                if(ez == N[2] - 1){
                    if(ex == 0){
                        arrV[ez][ey][ex][iu_element] = (arrA_x[ez][ey][ex][iu_element] - arrA_x[ez - 1][ey][ex][iu_element]) - (arrA_z[ez][ey][ex + 1][iu_element] - arrA_z[ez][ey][ex][iu_element]);
                    }
                    else if (ex == N[0] - 1){
                        arrV[ez][ey][ex][iu_element] = (arrA_x[ez][ey][ex][iu_element] - arrA_x[ez - 1][ey][ex][iu_element]) - (arrA_z[ez][ey][ex][iu_element] - arrA_z[ez][ey][ex - 1][iu_element]);
                    }
                    else {
                        arrV[ez][ey][ex][iu_element] = (arrA_x[ez][ey][ex][iu_element] - arrA_x[ez - 1][ey][ex][iu_element]) - (arrA_z[ez][ey][ex + 1][iu_element] - arrA_z[ez][ey][ex - 1][iu_element])/2.0;
                    }
                }

                if(ex == 0){
                    if(ey == 0){
                        arrW[ez][ey][ex][iu_element] = (arrA_y[ez][ey][ex + 1][iu_element] - arrA_y[ez][ey][ex][iu_element]) - (arrA_x[ez][ey + 1][ex][iu_element] - arrA_x[ez][ey][ex][iu_element]);
                    }
                    else if (ey == N[1] - 1){
                        arrW[ez][ey][ex][iu_element] = (arrA_y[ez][ey][ex + 1][iu_element] - arrA_y[ez][ey][ex][iu_element]) - (arrA_x[ez][ey][ex][iu_element] - arrA_x[ez][ey - 1][ex][iu_element]);
                    }
                    else {
                        arrW[ez][ey][ex][iu_element] = (arrA_y[ez][ey][ex + 1][iu_element] - arrA_y[ez][ey][ex][iu_element]) - (arrA_x[ez][ey + 1][ex][iu_element] - arrA_x[ez][ey - 1][ex][iu_element])/2.0;
                    }
                }    

                if(ex == N[0] - 1){
                    if(ey == 0){
                        arrW[ez][ey][ex][iu_element] = (arrA_y[ez][ey][ex][iu_element] - arrA_y[ez][ey][ex - 1][iu_element]) - (arrA_x[ez][ey + 1][ex][iu_element] - arrA_x[ez][ey][ex][iu_element]);
                    }
                    else if (ey == N[1] - 1){
                        arrW[ez][ey][ex][iu_element] = (arrA_y[ez][ey][ex][iu_element] - arrA_y[ez][ey][ex - 1][iu_element]) - (arrA_x[ez][ey][ex][iu_element] - arrA_x[ez][ey - 1][ex][iu_element]);
                    }
                    else {
                        arrW[ez][ey][ex][iu_element] = (arrA_y[ez][ey][ex][iu_element] - arrA_y[ez][ey][ex - 1][iu_element]) - (arrA_x[ez][ey + 1][ex][iu_element] - arrA_x[ez][ey - 1][ex][iu_element])/2.0;
                    }
                }         
            }
        }
    }

    DMStagVecRestoreArray(dmGrid_centered, vecU, &arrU);
    DMLocalToGlobal(dmGrid_centered, vecU, INSERT_VALUES, U_x);
    DMRestoreLocalVector(dmGrid_centered, &vecU);

    DMStagVecRestoreArray(dmGrid_centered, vecV, &arrV);
    DMLocalToGlobal(dmGrid_centered, vecV, INSERT_VALUES, U_y);
    DMRestoreLocalVector(dmGrid_centered, &vecV);

    DMStagVecRestoreArray(dmGrid_centered, vecW, &arrW);
    DMLocalToGlobal(dmGrid_centered, vecW, INSERT_VALUES, U_z);
    DMRestoreLocalVector(dmGrid_centered, &vecW);

    DMStagVecRestoreArrayRead(dmGrid_centered, vecA_x, &arrA_x);
    DMRestoreLocalVector(dmGrid_centered, &vecA_x);

    DMStagVecRestoreArrayRead(dmGrid_centered, vecA_y, &arrA_y);
    DMRestoreLocalVector(dmGrid_centered, &vecA_y);

    DMStagVecRestoreArrayRead(dmGrid_centered, vecA_z, &arrA_z);
    DMRestoreLocalVector(dmGrid_centered, &vecA_z);

    PetscFunctionReturn(0);
}

};


#endif // IMPLICIT_POTENTIAL_SOLVER_HPP

// Code to compute velocoty from vorticity with the Biot Savart Kernel. Not practical for medium/large simulations

/*double (*wx)[ny][nz] = (double(*)[ny][nz]) malloc(nx * ny * nz * sizeof(double));
double (*wy)[ny][nz] = (double(*)[ny][nz]) malloc(nx * ny * nz * sizeof(double));
double (*wz)[ny][nz] = (double(*)[ny][nz]) malloc(nx * ny * nz * sizeof(double));

for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
        for(int k = 0; k < nz; k++){
            wx[i][j][k] = wx_ref(i, j, k);
            wy[i][j][k] = wy_ref(i, j, k);
            wz[i][j][k] = wz_ref(i, j, k);
        }
    }
}

DM dmGrid_centered;
CreateGrid(&dmGrid_centered, 0, 0, 1, nx, ny, nz, 0.0, nx, 0.0, ny, 0.0, nz);

Mat A;
DMCreateMatrix(dmGrid_centered, &A);
assemble_lhs(dmGrid_centered, A);

Vec rhs_x, rhs_y, rhs_z, sol_x, sol_y, sol_z;
DMCreateGlobalVector(dmGrid_centered, &rhs_x);
DMCreateGlobalVector(dmGrid_centered, &rhs_y);
DMCreateGlobalVector(dmGrid_centered, &rhs_z);
build_petsc_vec(dmGrid_centered, rhs_x, wx);
build_petsc_vec(dmGrid_centered, rhs_y, wy);
build_petsc_vec(dmGrid_centered, rhs_z, wz);
VecScale(rhs_x, -1.0);
VecScale(rhs_y, -1.0);
VecScale(rhs_z, -1.0);

DMCreateGlobalVector(dmGrid_centered, &sol_x);
DMCreateGlobalVector(dmGrid_centered, &sol_y);
DMCreateGlobalVector(dmGrid_centered, &sol_z);
solve_laplacian(dmGrid_centered, A, rhs_x, sol_x);
solve_laplacian(dmGrid_centered, A, rhs_y, sol_y);
solve_laplacian(dmGrid_centered, A, rhs_z, sol_z);

double (*Ax)[ny][nz] = (double(*)[ny][nz]) malloc(nx * ny * nz * sizeof(double));
double (*Ay)[ny][nz] = (double(*)[ny][nz]) malloc(nx * ny * nz * sizeof(double));
double (*Az)[ny][nz] = (double(*)[ny][nz]) malloc(nx * ny * nz * sizeof(double));

build_std_velocity(dmGrid_centered, sol_x, Ax);
build_std_velocity(dmGrid_centered, sol_y, Ay);
build_std_velocity(dmGrid_centered, sol_z, Az);

for (int i = 1; i < nx - 1; ++i) {
    for (int j = 1; j < ny - 1; ++j) {
        for (int k = 1; k < nz - 1; ++k) {
            // Componente ux = ∂Az/∂y - ∂Ay/∂z
            ux[i][j][k] = (Az[i][j+1][k] - Az[i][j-1][k]) / (2.0)
                        - (Ay[i][j][k+1] - Ay[i][j][k-1]) / (2.0);

            // Componente uy = ∂Ax/∂z - ∂Az/∂x
            uy[i][j][k] = (Ax[i][j][k+1] - Ax[i][j][k-1]) / (2.0)
                        - (Az[i+1][j][k] - Az[i-1][j][k]) / (2.0);

            // Componente uz = ∂Ay/∂x - ∂Ax/∂y
            uz[i][j][k] = (Ay[i+1][j][k] - Ay[i-1][j][k]) / (2.0)
                        - (Ax[i][j+1][k] - Ax[i][j-1][k]) / (2.0);
        }
    }
}*/

/*PetscViewer viewer_w;
DM DM_vorticity;
//DMStagCreateCompatibleDMStag(dmGrid_staggered_x, 0, 0, 1, 0, &DM_u);
Vec vortex_x;
DMStagVecSplitToDMDA(dmGrid_centered, sol, ELEMENT, 0, &DM_vorticity, &vortex_x);
PetscObjectSetName((PetscObject)vortex_x, "vortex_x");
char filename_vortex_x[50]; 
sprintf(filename_vortex_x, "%svortex_x%03zu.vtr", base_path, 0);

PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_centered), filename_vortex_x, FILE_MODE_WRITE, &viewer_w);
VecView(vortex_x, viewer_w);
VecDestroy(&vortex_x);
DMDestroy(&DM_vorticity);
PetscViewerDestroy(&viewer_w);*/

// Calcoliamo la velocità usando la formula di Biot-Savart
/*double factor = 1.0 / (4.0 * M_PI);

for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            // Punto target (dove calcoliamo la velocità)

            // Scorriamo su tutti i punti del dominio per sommare i contributi
            for (int ip = 0; ip < nx; ip++) {
                for (int jp = 0; jp < ny; jp++) {
                    for (int kp = 0; kp < nz; kp++) {
                        // Distanza tra punto target e punto sorgente
                        double rxp = i - ip;
                        double ryp = j - jp;
                        double rzp = k - kp;
                        double r2 = rxp * rxp + ryp * ryp + rzp * rzp;
                        double r3 = pow(r2, 1.5);

                        // Evitiamo il contributo quando r2 è zero (punto singolare)
                        if (r2 == 0.0) continue;

                        // Contributo della vorticità
                        double wxp = wx[ip][jp][kp];
                        double wyp = wy[ip][jp][kp];
                        double wzp = wz[ip][jp][kp];

                        // Prodotto vettoriale w x r
                        double vx = wyp * rzp - wzp * ryp;
                        double vy = wzp * rxp - wxp * rzp;
                        double vz = wxp * ryp - wyp * rxp;

                        // Aggiorniamo la velocità
                        ux[i][j][k] += factor * vx / r3;
                        uy[i][j][k] += factor * vy / r3;
                        uz[i][j][k] += factor * vz / r3;
                    }
                }
            }
        }
    }
}*/