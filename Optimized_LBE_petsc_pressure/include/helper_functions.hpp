
PetscErrorCode const migrate_to_petsc(DM const & dmGrid_centered, Vec & rhs, double (*vec)[ny][nz])
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


void write_to_file(const std::string& filename, const std::string& content) {
    std::ofstream file(filename, std::ios::app); // Open file in append mode
    if (file.is_open()) {
        file << content << std::endl;
        file.close();
    } else {
        std::cerr << "Failed to open " << filename << " for writing." << std::endl;
    }
}

PetscErrorCode exodus(size_t n, Vec const & U_x, Vec const & U_y, Vec const & U_z, Vec const & U_mag, Vec const & Rho, DM const & dmGrid_centered)
{    

    PetscFunctionBegin;
    PetscViewer viewer_U_x_vtk;
    PetscObjectSetName((PetscObject)U_x, "U_x_vtk");
    char filename_U_x_vtk[100];
    sprintf(filename_U_x_vtk, "%sU_x%03zu.vtr", base_path, n);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_centered), filename_U_x_vtk, FILE_MODE_WRITE, &viewer_U_x_vtk);
    VecView(U_x, viewer_U_x_vtk);
    PetscViewerDestroy(&viewer_U_x_vtk);

    PetscViewer viewer_U_y_vtk;
    PetscObjectSetName((PetscObject)U_y, "U_y_vtk");
    char filename_U_y_vtk[100];
    sprintf(filename_U_y_vtk, "%sU_y%03zu.vtr", base_path, n);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_centered), filename_U_y_vtk, FILE_MODE_WRITE, &viewer_U_y_vtk);
    VecView(U_y, viewer_U_y_vtk);
    PetscViewerDestroy(&viewer_U_y_vtk);

    PetscViewer viewer_U_z_vtk;
    PetscObjectSetName((PetscObject)U_z, "U_z_vtk");
    char filename_U_z_vtk[100];
    sprintf(filename_U_z_vtk, "%sU_z%03zu.vtr", base_path, n);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_centered), filename_U_z_vtk, FILE_MODE_WRITE, &viewer_U_z_vtk);
    VecView(U_z, viewer_U_z_vtk);
    PetscViewerDestroy(&viewer_U_z_vtk);

    PetscViewer viewer_U_mag_vtk;
    PetscObjectSetName((PetscObject)U_mag, "U_mag_vtk");
    char filename_U_mag_vtk[100];
    sprintf(filename_U_mag_vtk, "%sU_mag%03zu.vtr", base_path, n);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_centered), filename_U_mag_vtk, FILE_MODE_WRITE, &viewer_U_mag_vtk);
    VecView(U_mag, viewer_U_mag_vtk);
    PetscViewerDestroy(&viewer_U_mag_vtk);


    PetscViewer viewer_Rho_vtk;
    PetscObjectSetName((PetscObject)Rho, "Rho_vtk");
    char filename_Rho_vtk[100];
    sprintf(filename_Rho_vtk, "%sRho%03zu.vtr", base_path, n);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_centered), filename_Rho_vtk, FILE_MODE_WRITE, &viewer_Rho_vtk);
    VecView(Rho, viewer_Rho_vtk);
    PetscViewerDestroy(&viewer_Rho_vtk);

    PetscFunctionReturn(0);
}

   
