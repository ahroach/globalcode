/* Sets up a flow background state from azimuthal velocity data
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "global.h"
#define MAXPOINTS 1000

// Number of fit coefficients
#define NCOEFFS 8

//nbreak = ncoeffs + 2 - k = ncoeffs -2 since k = 4
#define NBREAK (NCOEFFS - 2)

ROTATION_STRUCT *dataprofile(PARAMS_STRUCT *params, GRID_STRUCT *grid) {
  
  FILE *datafile;
  int points;
  char datastring[256];
  double r[MAXPOINTS];
  double omega[MAXPOINTS];

  // Allocate memory
  ROTATION_STRUCT *rotation = malloc(sizeof(ROTATION_STRUCT));
  rotation->omega = malloc((grid->numcells)*sizeof(double));

  //Set a and b parameters to zero, since they are not used in this case

  rotation->a = 0;
  rotation->b = 0;

  //Read data from file.  File should be named "profile.dat", and the
  //data should have two columns, the first the radial position in cm,
  //and the second the angular velocity in 1/sec.

  datafile = fopen("profile.dat", "r");
  
  if (datafile == NULL) {
    fprintf(stderr,"Error.  Azimuthal profile file profile.dat failed to open.\n");
  }
  
  points = 0;
  while(fgets(datastring, 256, datafile) != NULL) {
    sscanf(datastring, "%le%le", &r[points], &omega[points]);
    points++;
    if (points >= MAXPOINTS) {
      fprintf(stderr, "Error.  Number of datapoints exceeds MAXPOINTS");
      break;
    }
  }

  // Now fit cubic splines to the data, using gsl routines

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline
    = gsl_spline_alloc (gsl_interp_cspline, points);
  
  gsl_spline_init(spline, r, omega, points);

  
  for (int i = 0; i < grid->numcells; i++) {
    //gsl doesn't like to extrapolate. If we're inside of the furthest measured point
    //just use that innermost point
    if (grid->r[i] < r[0]) {
      rotation->omega[i] = omega[0];
    } else if (grid->r[i] > r[points-1]) {
      //Similarly if we're outside of the furthest measured point    
      rotation->omega[i] = omega[points-1];
    } else {
      //And in the middle interpolate.
      rotation->omega[i] = gsl_spline_eval(spline, grid->r[i], acc);
    }
  }

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  
  return rotation;
}

  
