#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "global.h"


double complex find_sigma(COMPRESSED_MATRIX *matrix,
                          PARAMS_STRUCT *params, GRID_STRUCT *grid,
                          ROTATION_STRUCT *rotation,
                          ARPACK_CONTROL *arpack_params) {
  //Do the iterate solution to determine sigma, making sure not to miss
  //any large eigenvalues.

  double complex finalsigma;
  double complex max_eigenvalue;
  ARPACK_CONTROL *temp_arpack_params;
  RESULTS_STRUCT *results;

  temp_arpack_params = malloc(sizeof(ARPACK_CONTROL));
  temp_arpack_params->nummodes = 10;
  //Need to make initial value of sigma large enough to not miss any growing
  //modes. Outer rotation rate is brought into this for the cases where    
  //inner cylinder is not rotating.         
  temp_arpack_params->sigma = (5.0*fabs(rotation->omega[0]) +
			       5.0*fabs(rotation->omega[grid->ie]))
    - I*params->m*(0.5*rotation->omega[0] + 0.5*rotation->omega[grid->ie]);
  temp_arpack_params->tol = 1e-2;
  temp_arpack_params->maxiters = arpack_params->maxiters;
  strcpy(temp_arpack_params->which, arpack_params->which);
  for (int i = 0; i < 4; i++) {
    results = eigensolve(matrix, params, grid,
                         rotation, temp_arpack_params);
    if (results->nconv < 1) {
      fprintf(stderr, "Error!  No eigenvalues found in iteration.\n");
    } else {
      //Find eigenvalue with largest real part
      max_eigenvalue = results->lambda[0];
      for (int j = 0; j < results->nconv; j++) {
        if (creal(results->lambda[j]) > creal(max_eigenvalue)) {
          max_eigenvalue = results->lambda[j];
        }
      }
      if(params->VERBOSE) {
	fprintf(stdout, "Max eigenvalue of %g + I*%g found\n",
		creal(max_eigenvalue), cimag(max_eigenvalue));
      }

      //Now reset the parameters and run again.
      temp_arpack_params->tol = temp_arpack_params->tol*.5;
      temp_arpack_params->sigma = 0.1*temp_arpack_params->sigma
        + 0.9*max_eigenvalue;
      if(params->VERBOSE) {
	fprintf(stdout, "Choosing new sigma %g + I*%g\n",
		creal(temp_arpack_params->sigma),
		cimag(temp_arpack_params->sigma));
      }
    }
    free(results->lambda);
    free(results->z);
    free(results);
  }

  finalsigma = temp_arpack_params->sigma;
  free(temp_arpack_params);
  return finalsigma;
}
