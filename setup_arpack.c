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

  arpack_params->iterate = get_iparam("iterate", input_file_name);

  //We only need to grab sigma from the input file if we're not iterating.
  if(arpack_params->iterate) {
    arpack_params->sigma = 0;
  } else {
    arpack_params->sigma = get_dparam("sigma_r", input_file_name);
    imagpart = get_dparam("sigma_i", input_file_name);
    arpack_params->sigma = arpack_params->sigma + I*imagpart;
  }
  arpack_params->tol = get_dparam("tol", input_file_name);
  arpack_params->maxiters = get_iparam("maxiters", input_file_name);


  //Make sure the mode specified by "which" is valid. If not, just
  //default to "LM"
  get_sparam("which", input_file_name, arpack_params->which);
  if((strcmp(arpack_params->which,"LM")) &&
     (strcmp(arpack_params->which,"SM")) &&
     (strcmp(arpack_params->which,"LR")) &&
     (strcmp(arpack_params->which,"SR")) &&
     (strcmp(arpack_params->which,"LI")) &&
     (strcmp(arpack_params->which,"SI"))) {
    fprintf(stderr, "Defaulting to 'LM' for value of 'which' in call to ZNAUPD.\n");
    strcpy(arpack_params->which, "LM");

  }

  return arpack_params;
}
