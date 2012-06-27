/* Problem generator for the MRI global code
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "global.h"

void probgen(char input_file_name[], PARAMS_STRUCT *params) {
  params->eta = get_dparam("eta", input_file_name);
  params->nu = get_dparam("nu", input_file_name);
  params->rho = get_dparam("rho", input_file_name);
  params->B0 = get_dparam("B0", input_file_name);
  params->r1 = get_dparam("r1", input_file_name);
  params->r2 = get_dparam("r2", input_file_name);
  params->omega1 = get_dparam("omega1", input_file_name);
  params->omega2 = get_dparam("omega2", input_file_name);
  params->height = get_dparam("height", input_file_name);
  params->nmode = get_dparam("nmode", input_file_name);
  params->m = get_dparam("m", input_file_name);
  params->numcells = get_iparam("numcells", input_file_name);
  params->magnetic_bc = get_iparam("magnetic_bc", input_file_name);

  /* Convert Omega in RPM to Omega in rads/sec */
  params->omega1 = params->omega1*2*PI/60;
  params->omega2 = params->omega2*2*PI/60;

  /* Calculate derived quantities*/

  params->va = params->B0/sqrt(4.0*PI*params->rho);
  params->k = 2.0*params->nmode*PI/params->height;
  params->kva = params->k*params->va;
  params->Pm = params->nu/params->eta;
  params->Re = params->omega1*params->r1*(params->r2-params->r1)/
    params->nu;
  params->Rm = params->omega1*params->r1*(params->r2-params->r1)/
    params->eta;
  params->Ha = params->B0*sqrt(params->r1*(params->r2-params->r1))/
    sqrt(4.0*PI*params->rho*params->nu*params->eta);

  params->VERBOSE = get_iparam("verbose", input_file_name);

  /* Print information for user verification */
  fprintf(stdout,"\nMRI Nonaxisymmetric Global Stability Code\n\n");
  fprintf(stdout,"Nu = %g cm^2/sec\n", params->nu);
  fprintf(stdout,"Eta = %g cm^2/sec\n", params->eta);
  fprintf(stdout,"Rho = %g g/cm^3\n", params->rho);
  fprintf(stdout,"B0 = %g Gauss\n", params->B0);
  fprintf(stdout,"r1 = %g cm\n", params->r1);
  fprintf(stdout,"r2 = %g cm\n", params->r2);
  fprintf(stdout,"omega1 = %g rad/sec\n", params->omega1);
  fprintf(stdout,"omega2 = %g rad/sec\n", params->omega2);
  fprintf(stdout,"height = %g cm\n", params->height);
  fprintf(stdout,"nmode = %g\n", params->nmode);
  fprintf(stdout,"m = %g\n", params->m);
  fprintf(stdout,"va = %g\n", params->va);
  fprintf(stdout,"k = %g\n", params->k);
  fprintf(stdout,"Pm = %g\n", params->Pm);
  fprintf(stdout,"Re = %g\n", params->Re);
  fprintf(stdout,"Rm = %g\n", params->Rm);
  fprintf(stdout,"Ha = %g\n", params->Ha);
}
