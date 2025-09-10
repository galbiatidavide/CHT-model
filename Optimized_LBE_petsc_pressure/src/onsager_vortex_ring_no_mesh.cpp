
#include "implicit_potential_solver.hpp"
#include "helper_functions.hpp"
#include "kernel.hpp"

int main(int argc, char *argv[]) {

    using namespace problem_setting;

    PetscInitialize(&argc, &argv, NULL, NULL);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, nx, ny, nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 27, 1, NULL, NULL, NULL, &grid);
    DMSetFromOptions(grid);
    DMSetUp(grid);
    DMDASetUniformCoordinates(grid, 0, nx, 0, ny, 0, nz);
    DMCreateGlobalVector(grid, &F_eq);
    DMCreateGlobalVector(grid, &F_temp);
    DMCreateGlobalVector(grid, &F);
    DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, nx, ny, nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, NULL, &macro_grid);
    DMSetFromOptions(macro_grid);
    DMSetUp(macro_grid);
    DMDASetUniformCoordinates(macro_grid, 0, nx, 0, ny, 0, nz);
    DMCreateGlobalVector(macro_grid, &Rho);
    DMCreateGlobalVector(macro_grid, &U_y);
    DMCreateGlobalVector(macro_grid, &U_z);
    DMCreateGlobalVector(macro_grid, &U_x);
    DMCreateGlobalVector(macro_grid, &U_mag);
    DMCreateGlobalVector(macro_grid, &U_mag_temp);
    DMCreateGlobalVector(macro_grid, &P);

    //From the vorticity field, compute the velocity field
    /*{           
        ImplicitPotentialSolver laplacian_solver;
        laplacian_solver.solve_problem();
        laplacian_solver.build_std_velocity(ux, uy, uz, rank, local_nx); 
    }*/


    PetscInt xs, ys, zs, xm, ym, zm;
    DMDAGetCorners(grid, &xs, &ys, &zs, &xm, &ym, &zm);
    PetscScalar ***ux_array, ***uy_array, ***uz_array, ***umag_array, ***p_array;
    DMDAVecGetArray(macro_grid, U_x, &ux_array);
    DMDAVecGetArray(macro_grid, U_y, &uy_array);
    DMDAVecGetArray(macro_grid, U_z, &uz_array);
    DMDAVecGetArray(macro_grid, U_mag, &umag_array);
    DMDAVecGetArray(macro_grid, P, &p_array);

    for (PetscInt k = zs; k < zs + zm; ++k) {
        for (PetscInt j = ys; j < ys + ym; ++j) {
            for (PetscInt i = xs; i < xs + xm; ++i) {
                ux_array[k][j][i] = uxRef(i + 0.5, j + 0.5, k + 0.5, 0);
                uy_array[k][j][i] = uyRef(i + 0.5, j + 0.5, k + 0.5, 0);
                uz_array[k][j][i] = uzRef(i + 0.5, j + 0.5, k + 0.5, 0);
                umag_array[k][j][i] = std::sqrt(ux_array[k][j][i] * ux_array[k][j][i] + uy_array[k][j][i] * uy_array[k][j][i] + uz_array[k][j][i] * uz_array[k][j][i]);
                p_array[k][j][i] = -0.5 * umag_array[k][j][i] * umag_array[k][j][i];

            }
        }
    }

    DMDAVecRestoreArray(macro_grid, U_x, &ux_array);
    DMDAVecRestoreArray(macro_grid, U_y, &uy_array);
    DMDAVecRestoreArray(macro_grid, U_z, &uz_array);
    DMDAVecRestoreArray(macro_grid, U_mag, &umag_array);
    DMDAVecRestoreArray(macro_grid, P, &p_array);



    DMDAGetCorners(grid, &xs, &ys, &zs, &xm, &ym, &zm);
    PetscScalar ***rho_array;
    DMDAVecGetArray(macro_grid, Rho, &rho_array);
    DMDAVecGetArray(macro_grid, P, &p_array);
    for (PetscInt k = zs; k < zs + zm; ++k) {
        for (PetscInt j = ys; j < ys + ym; ++j) {
            for (PetscInt i = xs; i < xs + xm; ++i) {
                rho_array[k][j][i] = 1.0;
            }
        }
    }

    DMDAVecRestoreArray(macro_grid, Rho, &rho_array);
    DMDAVecRestoreArray(macro_grid, P, &p_array);

    PetscScalar ****f_eq_array;
    PetscScalar ****f_temp_array;
    DMDAVecGetArrayDOF(grid, F_eq, &f_eq_array);
    DMDAVecGetArrayDOF(grid, F_temp, &f_temp_array);
    DMDAVecGetArray(macro_grid, Rho, &rho_array);
    DMDAVecGetArray(macro_grid, U_x, &ux_array);
    DMDAVecGetArray(macro_grid, U_y, &uy_array);
    DMDAVecGetArray(macro_grid, U_z, &uz_array);
    DMDAGetCorners(grid, &xs, &ys, &zs, &xm, &ym, &zm);
    for (PetscInt k = zs; k < zs + zm; ++k) {
        for (PetscInt j = ys; j < ys + ym; ++j) {
            for (PetscInt i = xs; i < xs + xm; ++i) {
                for (int a = 0; a < q; a++) {
                    f_eq_array[k][j][i][a] = w[a]*rho_array[k][j][i]*(1.0+3.0*(ux_array[k][j][i]*ex[a]+uy_array[k][j][i]*ey[a] + uz_array[k][j][i]*ez[a])+4.5*pow(ux_array[k][j][i]*ex[a]+uy_array[k][j][i]*ey[a] + uz_array[k][j][i]*ez[a],2.0)-1.5*(pow(ux_array[k][j][i],2.0)+pow(uy_array[k][j][i],2.0) + pow(uz_array[k][j][i],2.0)));
                    f_temp_array[k][j][i][a] = f_eq_array[k][j][i][a];

                }
            }
        }
    }
    DMDAVecRestoreArrayDOF(grid, F_eq, &f_eq_array);
    DMDAVecRestoreArrayDOF(grid, F_temp, &f_temp_array);
    DMDAVecRestoreArray(macro_grid, Rho, &rho_array);
    DMDAVecRestoreArray(macro_grid, U_x, &ux_array);
    DMDAVecRestoreArray(macro_grid, U_y, &uy_array);
    DMDAVecRestoreArray(macro_grid, U_z, &uz_array);


    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "Starting simulation" << std::endl;

    
  
    for(int n=0; n<iter; n++){

    //int m = 0;  

    stream_periodic_linkwise();
    compute_macroscopic_quantities();
    compute_equilibrium();
    collide();

        
    compute_convergence(rank, n);

    exodus(n, U_x, U_y, U_z, U_mag, Rho, macro_grid);

    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Simulation completed in " << elapsed_time << " seconds." << std::endl;

    VecDestroy(&U_x);
    VecDestroy(&U_y);
    VecDestroy(&U_z);
    VecDestroy(&U_mag);
    VecDestroy(&Rho);
    VecDestroy(&F_eq);
    VecDestroy(&F_temp);
    VecDestroy(&F);
    VecDestroy(&P);
    VecDestroy(&U_mag_temp);
    
    DMDestroy(&grid);
    DMDestroy(&macro_grid);
    PetscFinalize();
    
    return 0;
}

