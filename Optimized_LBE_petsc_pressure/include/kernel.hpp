
PetscErrorCode stream_periodic_linkwise(){
    int ix, iy, iz;
    PetscFunctionBegin;
    PetscInt xs, ys, zs, xm, ym, zm;
    PetscScalar ****f_temp_array;
    PetscScalar ****f_array;
    PetscScalar ***ux, ***uy, ***uz;
    Vec ux_local, uy_local, uz_local;
    DMCreateLocalVector(macro_grid, &ux_local);
    DMCreateLocalVector(macro_grid, &uy_local);
    DMCreateLocalVector(macro_grid, &uz_local);
    DMGlobalToLocalBegin(macro_grid, U_x, INSERT_VALUES, ux_local);
    DMGlobalToLocalEnd(macro_grid, U_x, INSERT_VALUES, ux_local);
    DMGlobalToLocalBegin(macro_grid, U_y, INSERT_VALUES, uy_local);
    DMGlobalToLocalEnd(macro_grid, U_y, INSERT_VALUES, uy_local);
    DMGlobalToLocalBegin(macro_grid, U_z, INSERT_VALUES, uz_local); 
    DMGlobalToLocalEnd(macro_grid, U_z, INSERT_VALUES, uz_local);
    DMDAVecGetArray(macro_grid, ux_local, &ux);
    DMDAVecGetArray(macro_grid, uy_local, &uy);
    DMDAVecGetArray(macro_grid, uz_local, &uz);
    DMDAVecGetArrayDOF(grid, F, &f_array);
    Vec F_temp_local;
    DMCreateLocalVector(grid, &F_temp_local);
    DMGlobalToLocalBegin(grid, F_temp, INSERT_VALUES, F_temp_local);
    DMGlobalToLocalEnd(grid, F_temp, INSERT_VALUES, F_temp_local);
    DMDAVecGetArrayDOF(grid, F_temp_local, &f_temp_array);
    DMDAGetCorners(grid, &xs, &ys, &zs, &xm, &ym, &zm);

    int b = 0;
    double v_x = 0.0;
    double v_y = 0.0;
    double v_z = 0.0;


    for (PetscInt k = zs; k < zs + zm; ++k) {
        for (PetscInt j = ys; j < ys + ym; ++j) {
            for (PetscInt i = xs; i < xs + xm; ++i) {     
                
                if(i == 0){
                    v_x = ux[k][j][i] + 0.5*(ux[k][j][i] - ux[k][j][i + 1]);
                    v_y = uy[k][j][i] + 0.5*(uy[k][j][i] - uy[k][j][i + 1]);
                    v_z = uz[k][j][i] + 0.5*(uz[k][j][i] - uz[k][j][i + 1]);
                }

                else if(i == nx - 1){
                    v_x = ux[k][j][i] + 0.5*(ux[k][j][i] - ux[k][j][i - 1]);
                    v_y = uy[k][j][i] + 0.5*(uy[k][j][i] - uy[k][j][i - 1]);
                    v_z = uz[k][j][i] + 0.5*(uz[k][j][i] - uz[k][j][i - 1]);

                }

                for (int a=0; a<q; a++) {

                    if (a == 0) {
                        b = 0;
                    } 
                    else if (a == 1) {
                        b = 2;
                    } 
                    else if (a == 2) {
                        b = 1;
                    } 
                    else if (a == 3) {
                        b = 4;
                    } 
                    else if (a == 4) {
                        b = 3;
                    } 
                    else if (a == 5) {
                        b = 6;
                    } 
                    else if (a == 6) {
                        b = 5;
                    } 
                    else if (a == 7) {
                        b = 8;
                    } 
                    else if (a == 8) {
                        b = 7;
                    } 
                    else if (a == 9) {
                        b = 10;
                    } 
                    else if (a == 10) {
                        b = 9;
                    } 
                    else if (a == 11) {
                        b = 12;
                    } 
                    else if (a == 12) {
                        b = 11;
                    } 
                    else if (a == 13) {
                        b = 14;
                    } 
                    else if (a == 14) {
                        b = 13;
                    } 
                    else if (a == 15) {
                        b = 16;
                    } 
                    else if (a == 16) {
                        b = 15;
                    } 
                    else if (a == 17) {
                        b = 18;
                    } 
                    else if (a == 18) {
                        b = 17;
                    } 
                    else if (a == 19) {
                        b = 20;
                    } 
                    else if (a == 20) {
                        b = 19;
                    } 
                    else if (a == 21) {
                        b = 22;
                    } 
                    else if (a == 22) {
                        b = 21;
                    } 
                    else if (a == 23) {
                        b = 24;
                    } 
                    else if (a == 24) {
                        b = 23;
                    } 
                    else if (a == 25) {
                        b = 26;
                    } 
                    else if (a == 26) {
                        b = 25;
                    } 

                    ix = i-ex[a]; 
                    iy = j-ey[a];
                    iz = k-ez[a];

                    if(iy < 0 or iy > ny - 1 or iz < 0 or iz > nz - 1){
                        f_array[k][j][i][a] = f_temp_array[k][j][i][b];
                    }
                    else
                    if(ix < 0){
                        f_array[k][j][i][a] = -f_temp_array[k][j][i][b] + 2*w[b]*(1.0 +2.46e-7)*(1 + (9.0/2.0)*(ex[b]*v_x + ey[b]*v_y + ez[b]*v_z)*(ex[b]*v_x + ey[b]*v_y + ez[b]*v_z) - (3.0/2.0)*(v_x*v_x + v_y*v_y + v_z*v_z));
                    }
                    else if(ix > nx - 1){
                        f_array[k][j][i][a] = -f_temp_array[k][j][i][b] + 2*w[b]*1.0*(1 + (9.0/2.0)*(ex[b]*v_x + ey[b]*v_y + ez[b]*v_z)*(ex[b]*v_x + ey[b]*v_y + ez[b]*v_z) - (3.0/2.0)*(v_x*v_x + v_y*v_y + v_z*v_z));

                    }
                    else
                    {

                    f_array[k][j][i][a] = f_temp_array[iz][iy][ix][a];
                    }
                }                     
            }
        }
    }

    DMDAVecRestoreArrayDOF(grid, F_temp_local, &f_temp_array);
    VecDestroy(&F_temp_local);
    DMDAVecRestoreArray(macro_grid, ux_local, &ux);
    DMDAVecRestoreArray(macro_grid, uy_local, &uy);
    DMDAVecRestoreArray(macro_grid, uz_local, &uz);
    VecDestroy(&ux_local);
    VecDestroy(&uy_local);
    VecDestroy(&uz_local);
    DMDAVecRestoreArrayDOF(grid, F, &f_array);
    PetscFunctionReturn(0);
}

PetscErrorCode compute_macroscopic_quantities(){
    PetscInt xs, ys, zs, xm, ym, zm;
    PetscFunctionBegin;
    PetscScalar ****f_array;
    DMDAVecGetArrayDOF(grid, F, &f_array);
    PetscScalar ***rho_array;
    //Vec local_rho;
    //DMGetLocalVector(macro_grid, &local_rho);
    DMDAVecGetArray(macro_grid, Rho, &rho_array);
    PetscScalar ***ux_array;
    DMDAVecGetArray(macro_grid, U_x, &ux_array);
    PetscScalar ***uy_array;
    DMDAVecGetArray(macro_grid, U_y, &uy_array);
    PetscScalar ***uz_array;
    DMDAVecGetArray(macro_grid, U_z, &uz_array);
    PetscScalar ***umag_array;
    DMDAVecGetArray(macro_grid, U_mag, &umag_array);
    PetscScalar ***umag_temp_array;
    DMDAVecGetArray(macro_grid, U_mag_temp, &umag_temp_array);
    DMDAGetCorners(grid, &xs, &ys, &zs, &xm, &ym, &zm);

    for (PetscInt k = zs; k < zs + zm; ++k) {
        for (PetscInt j = ys; j < ys + ym; ++j) {
            for (PetscInt i = xs; i < xs + xm; ++i) {
                rho_array[k][j][i] = 0;
                ux_array[k][j][i] = 0;
                uy_array[k][j][i] = 0;
                uz_array[k][j][i] = 0;
                
                for (int a=0; a<q; a++) {
                    rho_array[k][j][i] =  rho_array[k][j][i] + f_array[k][j][i][a];
                    ux_array[k][j][i]  =  ux_array[k][j][i]  + ex[a]*f_array[k][j][i][a];
                    uy_array[k][j][i]  =  uy_array[k][j][i]  + ey[a]*f_array[k][j][i][a];
                    uz_array[k][j][i]  =  uz_array[k][j][i]  + ez[a]*f_array[k][j][i][a];
                }
                umag_temp_array[k][j][i] = umag_array[k][j][i];
                double rho_inv = 1.0 / rho_array[k][j][i];
                ux_array[k][j][i] = (ux_array[k][j][i])*rho_inv;
                uy_array[k][j][i] = (uy_array[k][j][i])*rho_inv;
                uz_array[k][j][i] = (uz_array[k][j][i])*rho_inv;
                umag_array[k][j][i] = sqrt(ux_array[k][j][i]*ux_array[k][j][i] + uy_array[k][j][i]*uy_array[k][j][i] + uz_array[k][j][i]*uz_array[k][j][i]);
            }
        }
    }
    DMDAVecRestoreArrayDOF(grid, F, &f_array);
    DMDAVecRestoreArray(macro_grid, Rho, &rho_array);
    //DMRestoreLocalVector(macro_grid, &local_rho);
    DMDAVecRestoreArray(macro_grid, U_y, &uy_array);
    DMDAVecRestoreArray(macro_grid, U_z, &uz_array);
    DMDAVecRestoreArray(macro_grid, U_x, &ux_array);
    DMDAVecRestoreArray(macro_grid, U_mag, &umag_array);
    DMDAVecRestoreArray(macro_grid, U_mag_temp, &umag_temp_array);
    PetscFunctionReturn(0);
    
}

PetscErrorCode compute_equilibrium(){
    PetscInt xs, ys, zs, xm, ym, zm;
    PetscFunctionBegin;
    PetscScalar ****f_eq_array;
    DMDAVecGetArrayDOF(grid, F_eq, &f_eq_array);
    PetscScalar ***rho_array;
    DMDAVecGetArray(macro_grid, Rho, &rho_array);
    PetscScalar ***ux_array;
    DMDAVecGetArray(macro_grid, U_x, &ux_array);
    PetscScalar ***uy_array;
    DMDAVecGetArray(macro_grid, U_y, &uy_array);
    PetscScalar ***uz_array;
    DMDAVecGetArray(macro_grid, U_z, &uz_array);

    DMDAGetCorners(grid, &xs, &ys, &zs, &xm, &ym, &zm);
    for (PetscInt k = zs; k < zs + zm; ++k) {
        for (PetscInt j = ys; j < ys + ym; ++j) {
            for (PetscInt i = xs; i < xs + xm; ++i) {
                for (int a = 0; a < q; a++) {
                    f_eq_array[k][j][i][a] = w[a]*rho_array[k][j][i]*(1.0+3.0*(ux_array[k][j][i]*ex[a]+uy_array[k][j][i]*ey[a] + uz_array[k][j][i]*ez[a])+4.5*pow(ux_array[k][j][i]*ex[a]+uy_array[k][j][i]*ey[a] + uz_array[k][j][i]*ez[a],2.0)-1.5*(pow(ux_array[k][j][i],2.0)+pow(uy_array[k][j][i],2.0) + pow(uz_array[k][j][i],2.0)));
                }
            }
        }
    }
    DMDAVecRestoreArrayDOF(grid, F_eq, &f_eq_array);
    DMDAVecRestoreArray(macro_grid, Rho, &rho_array);
    DMDAVecRestoreArray(macro_grid, U_x, &ux_array);
    DMDAVecRestoreArray(macro_grid, U_y, &uy_array);
    DMDAVecRestoreArray(macro_grid, U_z, &uz_array);
    PetscFunctionReturn(0);
}

PetscErrorCode collide(){

    PetscInt xs, ys, zs, xm, ym, zm;
    PetscFunctionBegin;
    PetscScalar ****f_eq_array;
    PetscScalar ****f_temp_array;
    PetscScalar ****f_array;
    DMDAVecGetArrayDOF(grid, F_eq, &f_eq_array);
    DMDAVecGetArrayDOF(grid, F_temp, &f_temp_array);
    DMDAVecGetArrayDOF(grid, F, &f_array);
    PetscScalar ***rho_array;
    DMDAVecGetArray(macro_grid, Rho, &rho_array);
    PetscScalar ***ux_array;
    DMDAVecGetArray(macro_grid, U_x, &ux_array);
    PetscScalar ***uy_array;
    DMDAVecGetArray(macro_grid, U_y, &uy_array);
    PetscScalar ***uz_array;
    DMDAVecGetArray(macro_grid, U_z, &uz_array);

    DMDAGetCorners(grid, &xs, &ys, &zs, &xm, &ym, &zm);

    for (PetscInt k = zs; k < zs + zm; ++k) {
        for (PetscInt j = ys; j < ys + ym; ++j) {
            for (PetscInt i = xs; i < xs + xm; ++i) {
                double rho_inv = 1.0 / rho_array[k][j][i];
                double fx = 0.0;
                double fy = 0.0;
                double fz = 0.0;
                double vx = ux_array[k][j][i];
                double vy = uy_array[k][j][i];
                double vz = uz_array[k][j][i];

                double second_tensor[3][3] = {0.0}; 

                for (int r = 0; r < q; r++) {
                        double cx = ex[r];
                        double cy = ey[r];
                        double cz = ez[r];

                        double f_neq = f_array[k][j][i][r] - f_eq_array[k][j][i][r];
                        double c2 = cx * cx + cy * cy + cz * cz;
                        double trace_term = (1.0 / 3.0) * f_neq * c2;

                        second_tensor[0][0] += f_neq * cx * cx - trace_term;
                        second_tensor[0][1] += f_neq * cx * cy;
                        second_tensor[0][2] += f_neq * cx * cz;

                        second_tensor[1][0] += f_neq * cy * cx;
                        second_tensor[1][1] += f_neq * cy * cy - trace_term;
                        second_tensor[1][2] += f_neq * cy * cz;

                        second_tensor[2][0] += f_neq * cz * cx;
                        second_tensor[2][1] += f_neq * cz * cy;
                        second_tensor[2][2] += f_neq * cz * cz - trace_term;
                }

                for (int a = 0; a < q; a++) {
                    double d1 = 0.0;
                    double d2 = 0.0;
                    double d3 = 0.0;
                    double cx = ex[a] - vx;
                    double cy = ey[a] - vy;
                    double cz = ez[a] - vz;

                    double c2 = cx * cx + cy * cy + cz * cz;
                    double coeff = 4.5 * f_eq_array[k][j][i][a] * rho_inv;
                    double trace_term = c2 / 3.0;

                    double first_tensor[3][3];

                    first_tensor[0][0] = coeff * (cx * cx - trace_term);
                    first_tensor[0][1] = coeff * (cx * cy);
                    first_tensor[0][2] = coeff * (cx * cz);

                    first_tensor[1][0] = coeff * (cy * cx);
                    first_tensor[1][1] = coeff * (cy * cy - trace_term);
                    first_tensor[1][2] = coeff * (cy * cz);

                    first_tensor[2][0] = coeff * (cz * cx);
                    first_tensor[2][1] = coeff * (cz * cy);
                    first_tensor[2][2] = coeff * (cz * cz - trace_term);



                    second_tensor[0][0] -= 0.5 * (fx * vx + fx * vx);
                    second_tensor[0][1] -= 0.5 * (fx * vy + fy * vx);
                    second_tensor[0][2] -= 0.5 * (fx * vz + fz * vx);

                    second_tensor[1][0] -= 0.5 * (fy * vx + fx * vy);
                    second_tensor[1][1] -= 0.5 * (fy * vy + fy * vy);
                    second_tensor[1][2] -= 0.5 * (fy * vz + fz * vy);

                    second_tensor[2][0] -= 0.5 * (fz * vx + fx * vz);
                    second_tensor[2][1] -= 0.5 * (fz * vy + fy * vz);
                    second_tensor[2][2] -= 0.5 * (fz * vz + fz * vz);

                    double f_1 = 0.0;

                    f_1 += first_tensor[0][0] * second_tensor[0][0];
                    f_1 += first_tensor[0][1] * second_tensor[0][1];
                    f_1 += first_tensor[0][2] * second_tensor[0][2];

                    f_1 += first_tensor[1][0] * second_tensor[1][0];
                    f_1 += first_tensor[1][1] * second_tensor[1][1];
                    f_1 += first_tensor[1][2] * second_tensor[1][2];

                    f_1 += first_tensor[2][0] * second_tensor[2][0];
                    f_1 += first_tensor[2][1] * second_tensor[2][1];
                    f_1 += first_tensor[2][2] * second_tensor[2][2];

                    double omega = (1.0 - 1.0 / tau) * f_1;
                    f_temp_array[k][j][i][a] = f_eq_array[k][j][i][a] + omega + d1 + d2 + d3;
                }
            }
        }
    }

    DMDAVecRestoreArrayDOF(grid, F_eq, &f_eq_array);
    DMDAVecRestoreArrayDOF(grid, F_temp, &f_temp_array);
    DMDAVecRestoreArrayDOF(grid, F, &f_array);
    DMDAVecRestoreArray(macro_grid, Rho, &rho_array);
    DMDAVecRestoreArray(macro_grid, U_x, &ux_array);
    DMDAVecRestoreArray(macro_grid, U_y, &uy_array);
    DMDAVecRestoreArray(macro_grid, U_z, &uz_array);
    PetscFunctionReturn(0);

}

/*PetscErrorCode collide(){

    PetscInt xs, ys, zs, xm, ym, zm;
    PetscFunctionBegin;
    PetscScalar ****f_eq_array;
    PetscScalar ****f_temp_array;
    PetscScalar ****f_array;
    PetscScalar ***ux_array;
    DMDAVecGetArray(macro_grid, U_x, &ux_array);
    DMDAVecGetArrayDOF(grid, F_eq, &f_eq_array);
    DMDAVecGetArrayDOF(grid, F_temp, &f_temp_array);
    DMDAVecGetArrayDOF(grid, F, &f_array);


    DMDAGetCorners(grid, &xs, &ys, &zs, &xm, &ym, &zm);

    for (PetscInt k = zs; k < zs + zm; ++k) {
        for (PetscInt j = ys; j < ys + ym; ++j) {
            for (PetscInt i = xs; i < xs + xm; ++i) {
                double vx = ux_array[k][j][i];

                for (int a = 0; a < q; a++) {
                    double d1 = 0.0;
                    double d2 = 0.0;
                    double d3 = 0.0;
                    double omega = (-1.0 / tau) * (f_array[k][j][i][a] - f_eq_array[k][j][i][a]);
                    f_temp_array[k][j][i][a] = f_array[k][j][i][a] + omega + d1 + d2 + d3;
                }
            }
        }
    }

    DMDAVecRestoreArrayDOF(grid, F_eq, &f_eq_array);
    DMDAVecRestoreArrayDOF(grid, F_temp, &f_temp_array);
    DMDAVecRestoreArrayDOF(grid, F, &f_array);
    DMDAVecRestoreArray(macro_grid, U_x, &ux_array);

    PetscFunctionReturn(0);

}*/

PetscErrorCode compute_convergence(int rank, int n) {
    double temp1_local = 0.0;
    double temp2_local = 0.0;

    PetscInt xs, ys, zs, xm, ym, zm;
    PetscFunctionBegin;
    PetscScalar ***umag_array;
    DMDAVecGetArray(macro_grid, U_mag, &umag_array);
    PetscScalar ***umag_temp_array;
    DMDAVecGetArray(macro_grid, U_mag_temp, &umag_temp_array);
    
    DMDAGetCorners(grid, &xs, &ys, &zs, &xm, &ym, &zm);
    for (PetscInt k = zs; k < zs + zm; ++k) {
        for (PetscInt j = ys; j < ys + ym; ++j) {
            for (PetscInt i = xs; i < xs + xm; ++i) {
                double u_curr = umag_array[k][j][i];
                double u_prev = umag_temp_array[k][j][i];
                double diff = u_curr - u_prev;

                temp1_local += diff * diff;
                temp2_local += u_curr * u_curr;
            }
        }
    }

    DMDAVecRestoreArray(macro_grid, U_mag, &umag_array);
    DMDAVecRestoreArray(macro_grid, U_mag_temp, &umag_temp_array);

    double temp1_global = 0.0, temp2_global = 0.0;
    MPI_Reduce(&temp1_local, &temp1_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&temp2_local, &temp2_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double err_con = std::sqrt(temp1_global / temp2_global);
        std::cout << "[# Iter " << n << "] Convergence error is " << err_con << std::endl;
    }
    PetscFunctionReturn(0);
}








































