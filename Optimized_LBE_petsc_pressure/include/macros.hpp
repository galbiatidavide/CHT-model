#include <iostream>
#include <fstream> 
#include <chrono>
#include <vector>
#include <iomanip> 
#include <petscdmstag.h>
#include <petscksp.h>
#include <cstdio>
#include <petscdmda.h> 


#define ELEMENT          DMSTAG_ELEMENT
#define pi 3.14159265358979323846

constexpr int dx = 1;
constexpr int dy = 1;
constexpr int dz = 1;
constexpr int dt = 1;
constexpr double ro = 1.0;

constexpr int q = 27;           
constexpr int ex[] = {0,1,-1,0,0,0,0,1,-1,1,-1,0,0,1,-1,1,-1,0,0,1,-1,1,-1,1,-1,-1,1}; 
constexpr int ey[] = {0,0,0,1,-1,0,0,1,-1,0,0,1,-1,-1,1,0,0,1,-1,1,-1,1,-1,-1,1,1,-1}; 
constexpr int ez[] = {0,0,0,0,0,1,-1,0,0,1,-1,1,-1,0,0,-1,1,-1,1,1,-1,-1,1,1,-1,1,-1}; 
constexpr double w[]  = {8./27.,2./27., 2./27., 2./27., 2./27., 2./27., 2./27., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216.}; 

/**
 * @brief Parameters for the simulation
 * @param dx Grid spacing in x direction
 * @param dy Grid spacing in y direction
 * @param dz Grid spacing in z direction
 * @param iter Number of iterations
 */
constexpr int nx = 100;          
constexpr int ny = 20;         
constexpr int nz = 20;  
constexpr int iter = 3000; // Number of iterations

/**
 * @brief Simulation main quantities.
 */


DM grid;
DM macro_grid;
Vec F_eq;
Vec F_temp;
Vec F;
Vec U_x;
Vec U_y;
Vec U_z;
Vec Rho;
Vec U_mag;
Vec U_mag_temp;
Vec P;
