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
  char driver_name[256];
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

  //Find out which driver we're using, and pass off to the appropriate
  //handler function
  get_sparam("driver", input_file_name, driver_name);

  if (!strcmp(driver_name, "full")) {
    full_driver(input_file_name);
  } else if (!strcmp(driver_name, "batch")) {
    batch_driver(input_file_name);
  } else if (!strcmp(driver_name, "arpack")) {
    arpack_driver(input_file_name);
  } else {
    sprintf(stderr, "Error: unrecognized 'driver' in input file.\n");
  }

  return 0;
}



