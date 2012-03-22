/* Sets up a curl-free shear layer in the middle of the sample volume
 * as the background flow state
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "global.h"

ROTATION_STRUCT *shearlayer(PARAMS_STRUCT *params, GRID_STRUCT *grid,
			    double shear_width, double shear_radius) {

  ROTATION_STRUCT *rotation = malloc(sizeof(ROTATION_STRUCT));
  rotation->omega=malloc((grid->numcells)*sizeof(double));
  double a,b,c,d;
  
  /* The rotation profile Omega = A + B*tanh(C*(r-D))
   * where A = offset ; found by matching at r_out
   * B = scaling factor : (omega1-omega2)/(tanh(C*(r1-D))-tanh(C*(r2-D)))
   * Weird normalization is to account for the fact that tanh won't go to its
   * full +/- 1 values at the boundaries for large values of the shear layer
   * width.
   * C = 1/width
   * D = shear layer radius
   */
  
  /* Do the inner half */
  
  c = 1.0/shear_width;
  d = shear_radius;
  b = (params->omega1 - params->omega2)/
    (tanh(c*(params->r1-d)) - tanh(c*(params->r2-d)));
  a = params->omega2 - b*tanh(c*(params->r2-d));

  for(int i=0; i<=grid->numcells-1; i++) {
    rotation->omega[i] = a + b*tanh(c*(grid->r[i]-d));
  }

  return rotation;
}
