#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "global.h"


void arpack_handler(char *input_file_name)
{
  char profiletype[256];

  PARAMS_STRUCT *params;
  GRID_STRUCT *grid;
  ROTATION_STRUCT *rotation;
  OUTPUT_CONTROL *output_control;
  COMPRESSED_MATRIX *matrix;
  ARPACK_CONTROL *arpack_params;
  RESULTS_STRUCT *results;
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

  //Just run this once, and output the eigenvalues and eigenvectors          
  if (arpack_params->iterate == 1) {
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
