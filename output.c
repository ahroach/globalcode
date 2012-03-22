/* Dumps output data to a file
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "global.h"

void output(PARAMS_STRUCT *params, GRID_STRUCT *grid,
	    ROTATION_STRUCT *rotation, OUTPUT_CONTROL *output_control,
	    double complex *state, double complex gr)
{
  char output_file_path[256];
  char output_file_number[256];
  FILE *outputfile;
  double complex phiz, betaz;

  /* Construct a string for the filename of the form baseoutput#, i.e.
     output0009.  Then print tstate variables to this file.
  */

  sprintf(output_file_number, "%04i", output_control->filenum);
  strcpy(output_file_path, output_control->basefilename);
  strcat(output_file_path, output_file_number);
  
  outputfile = fopen(output_file_path, "w");
  if (outputfile == NULL) {
    fprintf(stderr, "Output file '%s' failed to open.\n", output_file_path);
    return;
  }

  fprintf(outputfile, "#gr = %g + I*%g, Re = %g, Rm = %g, Pm = %g\n", 
	  creal(gr), cimag(gr), params->Re, params->Rm, params->Pm);
  fprintf(outputfile, "#r, omega, abs(Beta_r), phase(Beta_r), abs(Beta_theta), phase(Beta_theta), abs(Beta_z), phase(Beta_z), abs(phi_r), phase(phi_r), abs(phi_theta), phase(phi_theta), abs(phi_z), phase(phi_z), abs(Pi), phase(Pi)\n");
  for (int i = grid->is-1; i <= grid->ie+1; i++) {    
    if (i == 0 || i == grid->numcells-1) {
      phiz = 0.0;
      betaz = 0.0;
    } else {
      phiz = (I*params->m*state[5*i + 3] + state[5*i + 2]
	      + (state[5*(i+1) + 2] - state[5*(i-1) + 2])/
	      (2.0*grid->dx))/(params->k*grid->r[i]);
      betaz = -(I*params->m*state[5*i + 1] + state[5*i + 0]
	      + (state[5*(i+1) + 0] - state[5*(i-1) + 0])/
	      (2.0*grid->dx))/(params->k*grid->r[i]);
    }

    fprintf(outputfile, "%8.5e %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e\n", 
	    grid->r[i], 
	    rotation->omega[i], cabs(state[5*i]), carg(state[5*i]),
	    cabs(state[5*i+1]), carg(state[5*i+1]), cabs(betaz),
	    carg(betaz), cabs(state[5*i+2]),
	    carg(state[5*i+2]), cabs(state[5*i+3]), carg(state[5*i+3]),
	    cabs(phiz), carg(phiz),
	    cabs(state[5*i+4]), carg(state[5*i+4]));
      }
  output_control->filenum++;
  fclose(outputfile);
}
