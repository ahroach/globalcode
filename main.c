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


int main(int ac, char **av)
{
  char input_file_name[256];
  int batch_run;
  int fullmode;


  /* Get input file from the command line*/
  
  if (ac < 2) {
    fprintf(stderr,"Please specify an input file name.\n");
    return 1;
  }

  strcpy(input_file_name, av[1]);

  //Now figure out if we're running in full mode or batch mode
  //or if we should just fall back to the standard ARPACK mode.

  fullmode = get_iparam("full", input_file_name);
  batch_run = get_iparam("batch", input_file_name);

  if (fullmode != 0) {
    fullmode_handler(input_file_name);
  } else if (batch_run != 0) {
    batchmode_handler(input_file_name);
  } else {
    arpack_handler(input_file_name);
  }

  return 0;
}



