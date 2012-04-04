/* MRI Global Stability Code with Nonaxisymmetric Terms
 * Written by Austin Roach
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <clapack.h>
#include <complex.h>
#include "global.h"

ARPACK_CONTROL *setup_arpack(char *input_file_name);
COMPRESSED_MATRIX *create_matrix(int numelems);

double complex find_sigma(COMPRESSED_MATRIX *matrix,
			  PARAMS_STRUCT *params, GRID_STRUCT *grid,
			  ROTATION_STRUCT *rotation,
			  ARPACK_CONTROL *arpack_params);

int main(int ac, char **av)
{
  char input_file_name[256];
  char profiletype[256];
  int iterate;
  int batch_run;

  PARAMS_STRUCT *params;
  GRID_STRUCT *grid;
  ROTATION_STRUCT *rotation;
  OUTPUT_CONTROL *output_control;
  COMPRESSED_MATRIX *matrix;
  ARPACK_CONTROL *arpack_params;
  RESULTS_STRUCT *results;

  /* Get input file from the command line*/
  
  if (ac < 2) {
    fprintf(stderr,"Please specify an input file name.\n");
    return 1;
  }

  strcpy(input_file_name, av[1]);

  //Get the physical parameters for the computation
  params = malloc(sizeof(PARAMS_STRUCT));  
  probgen(input_file_name, params);

  //Set up the grid, based on the physical parameters
  grid = gridgen(params);

  //Set up the rotation profile, from a number of different options.
  get_sparam("rotation", input_file_name, profiletype);
  if(!(strcmp(profiletype,"couette"))) {
    rotation = couette(params, grid);
  } else if (!(strcmp(profiletype,"dataprofile"))) {
    rotation = dataprofile(params, grid);
  } else if (!(strcmp(profiletype,"shearlayer"))) {
    double shear_width, shear_radius;
    shear_width = get_dparam("shear_width", input_file_name);
    shear_radius = get_dparam("shear_radius", input_file_name);
    rotation = shearlayer(params, grid, shear_width, shear_radius);
  } else {
    fprintf(stderr, "Error: No rotation profile specified.  Choosing Couette.\n");
    rotation = couette(params, grid);
  }

  //Set up the matrix for the computations.
  matrix = create_matrix(5*grid->numcells);

  //Setup the arpack parameters
  arpack_params = setup_arpack(input_file_name);


  //Now figure out if we're running in full mode or not.

  fullmode = get_iparam("full", input_file_name);

  if(fullmode) {
    //First signal that we're not using arpack in the output file
    arpack_params->sigma = NAN;
    arpack_params->tol = NAN;
    arpack_params->maxiters = NAN;
    arpack_params->nummodes = NAN;
    
    results = eigensolve(matrix, params, grid, rotation, arpack_params);
    //Setup the things needed to output data files
    output_control = malloc(sizeof(OUTPUT_CONTROL));
    output_control->filenum = 0;
    get_sparam("basefilename", input_file_name, output_control->basefilename); 
    wnetcdf(params, grid, rotation, output_control, arpack_params, results);
    free(results->lambda);
    free(results->z);
    free(results->residual);
    free(results);
    goto done;
  }
  

  iterate = get_iparam("iterate", input_file_name);
  batch_run = get_iparam("batch", input_file_name);

  if (batch_run == 1) {
    double batch_B0_init;
    double batch_B0_final;
    double B0_stepsize;
    int batch_B0_steps;
    char base_file_name[256];
    char output_file_path[256];
    FILE *outputfile;

    batch_B0_init = get_dparam("batch_B0_init", input_file_name);
    batch_B0_final = get_dparam("batch_B0_final", input_file_name);
    batch_B0_steps = get_iparam("batch_B0_steps", input_file_name);
    B0_stepsize = (batch_B0_final-batch_B0_init)/batch_B0_steps;
    get_sparam("basefilename", input_file_name, base_file_name);
    strcpy(output_file_path, base_file_name);
    outputfile = fopen(output_file_path, "w");
    fprintf(outputfile, "#B in Gauss, followed by real and imag parts of fastest growing modes\n");
    fclose(outputfile);

    for (int i = 0; i <= batch_B0_steps; i++) {
      params->B0 = batch_B0_init + B0_stepsize*i;
      params->va = params->B0/sqrt(4.0*PI*params->rho);
      params->kva = params->k*params->va;
      fprintf(stdout, "Running batch mode with B = %g Gauss\n", params->B0);
      if (iterate == 1) {
	arpack_params->sigma = find_sigma(matrix, params, grid, rotation,
					  arpack_params);
      }

      fprintf(stdout, "Entering main routine\n");
      results = eigensolve(matrix, params, grid, rotation, arpack_params);
      
      outputfile = fopen(output_file_path, "a");
      fprintf(outputfile, "%g ", params->B0);
      for (int j = 0; j < results->nconv; j++) {
	fprintf(outputfile, "%g %g ",
		creal(results->lambda[j]), cimag(results->lambda[j]));
      }
      fprintf(outputfile, "\n");
      fclose(outputfile);
      free(results->lambda);
      free(results->z);
      free(results->residual);
      free(results);
      
    }
      
  }else {
    //Just run this once, and output the eigenvalues and eigenvectors
    if (iterate == 1) {
      arpack_params->sigma = find_sigma(matrix, params, grid, rotation,
					arpack_params);
    }
    
    fprintf(stdout, "Entering main loop\n");

    results = eigensolve(matrix, params, grid, rotation, arpack_params);
    
    //Setup the things needed to output data files
    output_control = malloc(sizeof(OUTPUT_CONTROL));
    output_control->filenum = 0;
    get_sparam("basefilename", input_file_name, output_control->basefilename);
    
    wnetcdf(params, grid, rotation, output_control, arpack_params, results);
    free(results->lambda);
    free(results->z);
    free(results->residual);
    free(results);
  }  


 done:  
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
  return 0;
}


ARPACK_CONTROL *setup_arpack(char *input_file_name) {
  
  // Set up the arpack parameters from the input file
  ARPACK_CONTROL *arpack_params;
  double imagpart;
  
  arpack_params = malloc(sizeof(ARPACK_CONTROL));
			 
  arpack_params->nummodes = get_iparam("nummodes", input_file_name);
  if (!(arpack_params->nummodes > 0)) {
    fprintf(stderr, "Error reading nummodes.  Will return one eigenmode.\n");
    arpack_params->nummodes = 1;
  }

  arpack_params->sigma = get_dparam("sigma_r", input_file_name);
  imagpart = get_dparam("sigma_i", input_file_name);
  arpack_params->sigma = arpack_params->sigma + I*imagpart;
  arpack_params->tol = get_dparam("tol", input_file_name);
  arpack_params->maxiters = get_iparam("maxiters", input_file_name);

  return arpack_params;
}

COMPRESSED_MATRIX *create_matrix(int numelems) {
  COMPRESSED_MATRIX *matrix;
  matrix = malloc(sizeof(COMPRESSED_MATRIX));
  
  if (matrix == NULL) {
    fprintf(stderr, "Unable to allocate memory for matrix structure.\n");
    exit(1);
  } 

  /* Now's where we get into the meat of the computation.  First I need
   * to generate a matrix containing the coefficients of the terms in the
   * finite difference equation.
   * The order of the terms will be beta_r, beta_theta, phi_r, phi_theta,
   * pi at each gridpoint.
   */

  /* The matrix to hold these values will have 5*numcells rows and
   * 5*numcells columns.  These need to be complex numbers, so there
   * will be 2 doubles stored for each value
   */


  matrix->n = numelems; //Number of columns in A
  matrix->m = numelems; //Number of rows in A
  matrix->kl = 7; //Number of subdiagonals within the band of A
  matrix->ku = 7; //Number of superdiagonals within the band of A
  matrix->lda = 2*matrix->kl + matrix->ku + 1; //Leading dimension of array A
  matrix->ldab = matrix->kl + matrix->ku + 1; //Leading dimension of array Bb

  matrix->A = malloc(sizeof(double complex)*matrix->n*matrix->lda);
  if (matrix->A == NULL) {
    fprintf(stderr, "Unable to allocate memory for A.\n");
    exit(1);
  }

  matrix->B = malloc(sizeof(double complex)*matrix->n*matrix->lda);
  if (matrix->B == NULL) {
    fprintf(stderr, "Unable to allocate memory for B.\n");
    exit(1);
  }

  matrix->Bb = malloc(sizeof(double complex)*matrix->n*matrix->ldab);
  if (matrix->Bb == NULL) {
    fprintf(stderr, "Unable to allocate memory for Bb.\n");
    exit(1);
  }

  return matrix;
}


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
  temp_arpack_params->sigma = 10*abs(rotation->omega[0]) +
    2*abs(rotation->omega[grid->ie]);
  temp_arpack_params->tol = 1e-2;
  temp_arpack_params->maxiters = arpack_params->maxiters;
  for (int i = 0; i < 5; i++) {
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
      fprintf(stdout, "Max eigenvalue of %g + I*%g found\n",
	      creal(max_eigenvalue), cimag(max_eigenvalue));
      
      //Now reset the parameters and run again.
      temp_arpack_params->tol = temp_arpack_params->tol*.5;
      temp_arpack_params->sigma = 0.1*temp_arpack_params->sigma
	+ 0.9*max_eigenvalue;
      fprintf(stdout, "Choosing new sigma %g + I*%g\n",
	      creal(temp_arpack_params->sigma),
	      cimag(temp_arpack_params->sigma));
    }
    free(results->lambda);
    free(results->z);
    free(results);
  }
  
  finalsigma = temp_arpack_params->sigma;
  free(temp_arpack_params);
  return finalsigma;
}
