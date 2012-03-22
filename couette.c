/* Sets up a Couette flow background state for the MRI Global code
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "global.h"

ROTATION_STRUCT *couette(PARAMS_STRUCT *params, GRID_STRUCT *grid) {

  ROTATION_STRUCT *rotation = malloc(sizeof(ROTATION_STRUCT));
  
  rotation->omega=malloc((grid->numcells)*sizeof(double));


  /* Calculate a and b coefficients from equation for Couette flow:
   * Omega(r) = a + b/r**2
   */

  rotation->a = (pow(params->r2,2)*params->omega2 - 
		 pow(params->r1,2)*params->omega1)/
    (pow(params->r2,2)-pow(params->r1,2));
  rotation->b = (params->omega2 - rotation->a)*pow(params->r2,2);

  fprintf(stdout, "a = %5.3e, b= %5.3e\n", rotation->a, rotation->b);

  /* Use a and b coefficients to calculate Omega at every point on
   * the grid
   */

  for(int i=0; i <= grid->numcells-1; i++) {
    rotation->omega[i]=rotation->a+rotation->b/pow(grid->r[i],2);
  }

  params->zetabar = 2*rotation->a/sqrt(params->omega2*params->omega1);

  fprintf(stdout, "zetabar = %5.3e\n", params->zetabar);

  return rotation;
}
