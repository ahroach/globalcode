#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "global.h"

ARPACK_CONTROL *setup_arpack(char *input_file_name) {

  // Set up the arpack parameters from the input file

  ARPACK_CONTROL *arpack_params;
  double imagpart;

  arpack_params = malloc(sizeof(ARPACK_CONTROL));

  arpack_params->nummodes = get_iparam("nummodes", input_file_name);
  if (!(arpack_params->nummodes > 0)) {
    fprintf(stderr, "Error reading nummodes.  Will return one eigenmode if in ARPACK mode.\n");
    arpack_params->nummodes = 1;
  }

  arpack_params->sigma = get_dparam("sigma_r", input_file_name);
  imagpart = get_dparam("sigma_i", input_file_name);
  arpack_params->sigma = arpack_params->sigma + I*imagpart;
  arpack_params->tol = get_dparam("tol", input_file_name);
  arpack_params->maxiters = get_iparam("maxiters", input_file_name);

  return arpack_params;
}
