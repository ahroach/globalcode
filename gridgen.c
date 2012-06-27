/* Grid generator for the MRI global code
 * Sets up a radial grid that is uniform in x=log(r)
 * With spacing between grid points dx = log(R2/R1)/num
 * First and last grid points are dx/2 away from the boundary
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "global.h"

GRID_STRUCT *gridgen(PARAMS_STRUCT *params) {
  GRID_STRUCT *grid = malloc(sizeof(GRID_STRUCT));
  double x1;
  double x2;
  double dx;

  grid->numcells = params->numcells;
  grid->r = malloc((grid->numcells)*sizeof(double));
  grid->x = malloc((grid->numcells)*sizeof(double));
  grid->r2inv = malloc((grid->numcells)*sizeof(double));
  grid->diffuse = malloc((grid->numcells)*sizeof(double));
  grid->diffuse2 = malloc((grid->numcells)*sizeof(double));
  grid->is = 1; //Index of the first real grid cell
  grid->ie = grid->numcells-2; //Index of the last real grid cells

  x1 = log(params->r1);
  x2 = log(params->r2);
  dx = (x2-x1)/(grid->numcells -2); //Leave 2 cells for ghost zones
  grid->dx = dx;
  double dx2inv = 1.0/pow(dx,2);
  double m2 = pow(params->m, 2);
  double k2 = pow(params->k, 2);
  if(params->VERBOSE) {
    printf("dx = %g\n", dx);
  }
  
  /* Generate the logarithmic grid */
  for(int i=0;i<=grid->numcells-1;i++) {
    grid->r[i]=exp(x1+(i-0.5)*dx);
    grid->x[i]=x1+(i-0.5)*dx;
    grid->r2inv[i]=1.0/pow(grid->r[i],2);
    grid->diffuse[i] = 2*grid->r2inv[i]*dx2inv + grid->r2inv[i] + 
      grid->r2inv[i]*m2 + k2;
    grid->diffuse2[i] = grid->r2inv[i]*dx2inv;
  }

  return grid;
}
