#include "macros.hpp"

namespace problem_setting
{
constexpr const char *base_path = "results/";


constexpr double u_inlet = 1e-4; // Fixed to avoid compressibility error (Ma<<0.1)
constexpr double Re = 1000.0; // Reynolds number
constexpr double tau = u_inlet * ny / (3*Re) + 0.5; // Relaxation time proportional to Re


/*
* @brief Navier-Stokes solution in x-direction.

*/
constexpr PetscReal uxRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    // Set initial guess for velocity profile
    return u_inlet*(1 - (2*y/double(ny) - 1)*(2*y/(double(ny)) - 1))*(1 - (2*z/(double(nz)) - 1)*(2*z/(double(nz)) - 1));
    //return 0.0;

}

/*
* @brief Navier-Stokes solution in y-direction.

*/
constexpr PetscReal uyRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    return 0.0;
}

/*
* @brief Navier-Stokes solution in z-direction.
*
*/
constexpr PetscReal uzRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    return 0.0;
}



}

