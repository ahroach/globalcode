#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "global.h"


void batch_driver(char *input_file_name)
{
  int iterate;
  char profiletype[256];

  PARAMS_STRUCT *params;
  GRID_STRUCT *grid;
  ROTATION_STRUCT *rotation;
  COMPRESSED_MATRIX *matrix;
  ARPACK_CONTROL *arpack_params;
  RESULTS_STRUCT *results;

  double batch_B0_init;
  double batch_B0_final;
  double B0_stepsize;
  int batch_B0_steps;
  char base_file_name[256];
  char output_file_path[256];
  FILE *outputfile;


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

  //Now set up all of the batch mode parameters

  batch_B0_init = get_dparam("batch_B0_init", input_file_name);
  batch_B0_final = get_dparam("batch_B0_final", input_file_name);
  batch_B0_steps = get_iparam("batch_B0_steps", input_file_name);
  B0_stepsize = (batch_B0_final-batch_B0_init)/batch_B0_steps;
  get_sparam("basefilename", input_file_name, base_file_name);
  strcpy(output_file_path, base_file_name);
  outputfile = fopen(output_file_path, "w");
  fprintf(outputfile, "#B in Gauss, followed by real and imag parts of fastest growing modes\n");
  fclose(outputfile);

  //Cycle through all of the Bs, resetting the relevant derived quantities 
  //before each run.

  for (int i = 0; i <= batch_B0_steps; i++) {
    params->B0 = batch_B0_init + B0_stepsize*i;
    params->va = params->B0/sqrt(4.0*PI*params->rho);
    params->kva = params->k*params->va;
    fprintf(stdout, "Running batch mode with B = %g Gauss\n", params->B0);
    if (arpack_params->iterate == 1) {
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
