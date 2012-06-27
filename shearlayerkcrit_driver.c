#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
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

double mindampingrate(double k, void *params);

void err_handler(const char *reason, const char *file, int line,
		 int gsl_errno);

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
  const gsl_root_fsolver_type *Troot;
  const gsl_min_fminimizer_type *Tmin;
  gsl_root_fsolver *sroot;
  gsl_min_fminimizer *smin;
  double k_low, k_high, k_guess;
  double k_min = NAN;
  double k_max = NAN;
  double k_peak = NAN;
  double gr_peak;
  double errabs, errrel;
  double width_prefactor;
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
  width_prefactor = get_dparam("width_prefactor", input_file_name);
  E = params->nu/(0.5*fabs(params->omega1 - params->omega2) * 
		  pow((params->r2-params->r1),2));
  shear_width = width_prefactor*(params->r2-params->r1)*pow(E, 0.25);
  printf("Using shear layer width %g cm\n", shear_width);
  rotation = shearlayer(params, grid, shear_width, shear_radius);
  
  //Set up the matrix structure for the computations.
  matrix = create_matrix(5*grid->numcells);

  //Setup the ARPACK parameters
  arpack_params = setup_arpack(input_file_name);

  //Setup the output control structure
  output_control = malloc(sizeof(OUTPUT_CONTROL));

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
  F.function = &mindampingrate;
  F.params = &function_params;

  gsl_set_error_handler(&err_handler);

  /* Now we find the peak of the growth rate, by minimizing the
     damping rate. We set what we hope are reasonable numbers
     for the bounds and initial guess.
  */

  k_low = 0.01;
  k_high = 1000;
  k_guess = params->k;
  Tmin = gsl_min_fminimizer_brent;
  smin = gsl_min_fminimizer_alloc(Tmin);
  status = gsl_min_fminimizer_set(smin, &F, k_guess, k_low, k_high);
  //Make sure that we didn't thrown an error on initialization
  if (status == GSL_SUCCESS) {
    //Now iterate!
    iter = 0;
    do 
      {
	iter++;
	status = gsl_min_fminimizer_iterate(smin);
	//Make sure that we didn't thrown an error in the iteration routine
	if (status != GSL_SUCCESS) {
	  break;
	}
	
	params->k = gsl_min_fminimizer_x_minimum(smin);
	k_low = gsl_min_fminimizer_x_lower(smin);
	k_high = gsl_min_fminimizer_x_upper(smin);
	status = gsl_min_test_interval(k_low, k_high, errabs, errrel);
	
	if(status == GSL_SUCCESS) {
	  printf("Converged with k_peak=%g\n", params->k);
	}
      }
    while (status == GSL_CONTINUE && iter < max_iter);
    //Save the peak growth rate for printing later, then free the solver
    gr_peak = -gsl_min_fminimizer_f_minimum(smin);
  }
  gsl_min_fminimizer_free(smin);

  //Check to make sure we converged. If not, don't save the results.
  if (status == GSL_SUCCESS) {
    k_peak = params->k;
    //Now do a normal run with the chosen k
    arpack_params->sigma = find_sigma(matrix, params, grid, rotation,
				      arpack_params);
    results = eigensolve(matrix, params, grid, rotation, arpack_params);
    
    //Setup the structures needed to output the data files, and write them.
    get_sparam("basefilename", input_file_name, output_control->basefilename);
    strcat(output_control->basefilename, "_kpeak");
    wnetcdf(params, grid, rotation, output_control, arpack_params, results);
    
    free(results->lambda);
    free(results->z);
    free(results->residual);
    free(results);
  }


  /* Now do a root finding search for k_min. */

  //Set up the root solver.
  Troot = gsl_root_fsolver_brent;
  sroot = gsl_root_fsolver_alloc(Troot);

  //Set the initial bounds for the search. We're searching for k_min,
  //so search from 0 up to k_peak.
  k_low = 0;
  k_high = k_peak;
  status = gsl_root_fsolver_set(sroot, &F, k_low, k_high);
  //Make sure that we didn't thrown an error on initialization
  if (status == GSL_SUCCESS) {
    //Now iterate!
    iter = 0;
    do 
      {
	iter++;
	printf("iter = %i\n", iter);
	status = gsl_root_fsolver_iterate(sroot);
	//Make sure that we didn't thrown an error in the iteration routine
	if (status != GSL_SUCCESS) {
	  break;
	}
	
	params->k = gsl_root_fsolver_root(sroot);
	k_low = gsl_root_fsolver_x_lower(sroot);
	k_high = gsl_root_fsolver_x_upper(sroot);
	status = gsl_root_test_interval(k_low, k_high, errabs, errrel);
	
	if(status == GSL_SUCCESS) {
	  printf("Converged with k_min=%g\n", params->k);
	}
      }
    while (status == GSL_CONTINUE && iter < max_iter);
  }
  gsl_root_fsolver_free (sroot);

  //Check to make sure we converged. If not, don't save the results.
  if (status == GSL_SUCCESS) {
    k_min = params->k;
    //Now do a normal run with the chosen k
    arpack_params->sigma = find_sigma(matrix, params, grid, rotation,
				      arpack_params);
    results = eigensolve(matrix, params, grid, rotation, arpack_params);
    
    //Set the new file name, and write the output
    get_sparam("basefilename", input_file_name, output_control->basefilename);
    strcat(output_control->basefilename, "_kmin");
    wnetcdf(params, grid, rotation, output_control, arpack_params, results);
    
    free(results->lambda);
    free(results->z);
    free(results->residual);
    free(results);
  }


  /* Now move on to solving for k_max. */
  sroot = gsl_root_fsolver_alloc(Troot);

  //Set the initial bounds for the search. We're searching for k_max,
  //so search from k_peak to a large number
  k_low = k_peak;
  k_high = 10000;
  status = gsl_root_fsolver_set(sroot, &F, k_low, k_high);
  //Make sure that we didn't thrown an error on initialization
  if (status == GSL_SUCCESS) {
    //Now iterate!
    iter = 0;
    do 
      {
	iter++;
	status = gsl_root_fsolver_iterate(sroot);
	//Make sure that we didn't thrown an error in the iteration routine
	if (status != GSL_SUCCESS) {
	  break;
	}
	
	params->k = gsl_root_fsolver_root(sroot);
	k_low = gsl_root_fsolver_x_lower(sroot);
	k_high = gsl_root_fsolver_x_upper(sroot);
	status = gsl_root_test_interval(k_low, k_high, errabs, errrel);
	
	if(status == GSL_SUCCESS) {
	  printf("Converged with k_max=%g\n", params->k);
	}
      }
    while (status == GSL_CONTINUE && iter < max_iter);
  }
  gsl_root_fsolver_free (sroot);
    
  //Check to make sure we converged. If not, don't save the results.
  if (status == GSL_SUCCESS) {
    k_max = params->k;
    //Now do a normal run with the chosen k
    arpack_params->sigma = find_sigma(matrix, params, grid, rotation,
				      arpack_params);
    results = eigensolve(matrix, params, grid, rotation, arpack_params);
    
    //Set the new file name, and write the output
    get_sparam("basefilename", input_file_name, output_control->basefilename);
    strcat(output_control->basefilename, "_kmax");
    wnetcdf(params, grid, rotation, output_control, arpack_params, results);
    
    free(results->lambda);
    free(results->z);
    free(results->residual);
    free(results);
  }

  printf("Found k_min = %g, k_peak = %g, k_max = %g\n", k_min, k_peak, k_max);
  printf("Peak growth rate: %g\n", gr_peak);

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
  free(output_control);

  return;
}
 

double mindampingrate(double k, void *fnparams)
{
  /* This function calculates the minimumdamping rate for given k and
     params. Because all of our routines are set up to calculate
     growth rate, gamma, we find the maximum growth rate, and then
     return the negative of that as the minimum damping rate. This
     form is required since we want to find the peak growth rate, but
     it's easiest to use a function minimizer to find the minimum
     damping rate.
  */

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

  return -max_gr;
}
 

void err_handler(const char *reason, const char *file, int line,
		 int gsl_errno)
{
  fprintf(stderr, "GSL error number %i in file %s at line %i\n", gsl_errno,
	  file, line);
  fprintf(stderr, "Reason: %s\n", reason);
  return;
}

  
