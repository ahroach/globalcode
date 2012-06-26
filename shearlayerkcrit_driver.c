#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "global.h"

struct FUNCTION_PARAMS
{
  PARAMS_STRUCT *params;
  GRID_STRUCT *grid;
  ROTATION_STRUCT *rotation;
  COMPRESSED_MATRIX *matrix;
  ARPACK_CONTROL *arpack_params;
};

double growthrate(double k, void *params);

void shearlayerkcrit_driver(char *input_file_name)
{
  PARAMS_STRUCT *params;
  GRID_STRUCT *grid;
  ROTATION_STRUCT *rotation;
  COMPRESSED_MATRIX *matrix;
  ARPACK_CONTROL *arpack_params;
  RESULTS_STRUCT *results;
  OUTPUT_CONTROL *output_control;

  double shear_width, shear_radius, E;

  //Parameters needed for the root-finding routine
  int status;
  int iter=0, max_iter=50;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double k_low, k_high;
  double errabs, errrel;
  gsl_function F;
  struct FUNCTION_PARAMS function_params;

  //Get the physical parameters for the computation
  params = malloc(sizeof(PARAMS_STRUCT));
  probgen(input_file_name, params);
  
  //Set up the grid, based on the physical parameters
  grid = gridgen(params);

  //Set up the rotation profile of a shear layer. Derive the width
  //from the Ekman number, E=\nu/\Omega r^2, width = rE^(1/4)
  //Use r = (r2-r1) and Omega = (Omega1-Omega2)/2.

  shear_radius = get_dparam("shear_radius", input_file_name);
  E = params->nu/(0.5*fabs(params->omega1 - params->omega2) * 
		  pow((params->r2-params->r1),2));
  shear_width = (params->r2-params->r1)*pow(E, 0.25);
  printf("Using shear layer width %g cm\n", shear_width);
  rotation = shearlayer(params, grid, shear_width, shear_radius);
  
  //Set up the matrix strcture for the computations.
  matrix = create_matrix(5*grid->numcells);

  //Setup the ARPACK parameters
  arpack_params = setup_arpack(input_file_name);

  //Pull the error params from the input file to decide when
  //we have converged
  errabs = get_dparam("errabs", input_file_name);
  errrel = get_dparam("errrel", input_file_name);
  
  //Put pointers to all of our control structures in function_params
  function_params.params = params;
  function_params.grid = grid;
  function_params.rotation = rotation;
  function_params.matrix = matrix;
  function_params.arpack_params = arpack_params;
  
  //Assign the evaluation function and params structure to
  //the gsl_function
  F.function = &growthrate;
  F.params = &function_params;
  
  //Set up the root solver
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);

  //Set the initial bounds for the search. Use the k specified by
  //the input file for the lower bound on k.
  k_low = 0.001*2*PI/shear_width;
  k_high = 1e5*k_low;
  gsl_root_fsolver_set(s, &F, k_low, k_high);

  //Now iterate!
  do 
    {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      params->k = gsl_root_fsolver_root(s);
      k_low = gsl_root_fsolver_x_lower(s);
      k_high = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(k_low, k_high, errabs, errrel);
      
      if(status == GSL_SUCCESS) {
	printf("Converged with k=%g\n", params->k);
      }
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  //Check to make sure we converged. If not, don't save the results.
  if (status != GSL_SUCCESS) {
    return;
  }

  //Now do a normal run with the chosen k
  //Unfortunately, I don't think I can trust the results from the
  //last call to growthrate(), because I don't know that the last
  //evaluation of the function by the root-finding routine
  //is guaranteed to occur at the optimized k.
  arpack_params->sigma = find_sigma(matrix, params, grid, rotation,
				    arpack_params);
  results = eigensolve(matrix, params, grid, rotation, arpack_params);
  
  //Setup the structures needed to output the data files, and write them.
  output_control = malloc(sizeof(OUTPUT_CONTROL));
  output_control->filenum = 0;
  get_sparam("basefilename", input_file_name, output_control->basefilename);

  wnetcdf(params, grid, rotation, output_control, arpack_params, results);

  free(results->lambda);
  free(results->z);
  free(results->residual);
  free(results);
  free(matrix->A);
  free(matrix->B);
  free(matrix->Bb);
  free(matrix);
  free(params);
  free(grid->r);
  free(grid->x);
  free(grid->r2inv);
  free(grid);
  free(rotation->omega);
  free(rotation);

  return;
}
 

double growthrate(double k, void *fnparams)
{
  //Get a pointer casting the params as a function_params structure
  //so I can address the elements.
  struct FUNCTION_PARAMS *p = (struct FUNCTION_PARAMS *) fnparams;
  
  //Define these things locally so I'm not addressing p all the time.
  //So I should have local pointers to everything.
  
  PARAMS_STRUCT *params = p->params;
  GRID_STRUCT *grid = p->grid;
  ROTATION_STRUCT *rotation = p->rotation;
  COMPRESSED_MATRIX *matrix = p->matrix;
  ARPACK_CONTROL *arpack_params = p->arpack_params;
  RESULTS_STRUCT *results;
  
  double max_gr = 0;

  //Now change k in the params
  params->k = k;
  
  //And I have to recalculate the quantities in the grid,
  //since I have the matrix elements of the diffuse terms in there
  free(grid->r);
  free(grid->x);
  free(grid->r2inv);
  free(grid);
  grid = gridgen(params);
  
  //Now I'm ready to find that eigenvalue!
  arpack_params->sigma = find_sigma(matrix, params, grid, rotation,
				    arpack_params);
  results = eigensolve(matrix, params, grid, rotation, arpack_params);
  
  if (results->nconv < 1) {
    fprintf(stderr, "Error! No eigenvalues found for k=%g.\n", params->k);
    
  } else {
    //Find the eigenvalue with the largest real part
    max_gr = results->lambda[0];
    for (int i = 0; i < results->nconv; i++) {
      if (creal(results->lambda[i]) > creal(max_gr)) {
	max_gr = results->lambda[i];
      }
    }
  }
  
  //Now free up the results structure
  free(results->lambda);
  free(results->z);
  free(results);

  return max_gr;
}
 

  
