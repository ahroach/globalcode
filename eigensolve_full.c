/* Step the equations of the global code solving directly for all eigenvalues
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <clapack.h>
#include <complex.h>
#include <gsl/gsl_sf_bessel.h>
#include "global.h"

void wAelbg(int i, int j, COMPRESSED_MATRIX *matrix, double complex value);
void wBelbg(int i, int j, COMPRESSED_MATRIX *matrix, double complex value);


RESULTS_STRUCT *eigensolve_direct(COMPRESSED_MATRIX *matrix, 
				  PARAMS_STRUCT *params, GRID_STRUCT *grid,
				  ROTATION_STRUCT *rotation,
				  ARPACK_CONTROL *arpack_params) {

  //Physical parameters needed to fill in the matrix
  double eta = params->eta;
  double nu = params->nu;
  double dx2inv = 1.0/pow(grid->dx,2);
  double m2 = pow(params->m,2);
  double k2 = pow(params->k,2);
  double kva = params->kva;
  double complex facin, facout;
  double Im, Imp1, Km, Kmp1;
  double scalefactorpi, scalefactorv, scalefactorb, tempscalefactor;

  //Variables needed to run the LAPACK routine
  int info;
  int lwork;
  double complex *lambda;
  double complex *alpha;
  double complex *beta;
  double complex *vr;
  double complex *work;
  double complex *rwork;
  int ldvl;
  int ldvr;

  char jobvl;
  char jobvr;

  //Since we're not using the compressed matrix here, we need to reallocate
  //some more space
  free(matrix->A);
  free(matrix->B);

  matrix->A = malloc(sizeof(double complex)*matrix->n*matrix->n);
  if (matrix->A == NULL) {
    fprintf(stderr, "Unable to allocate memory for A.\n");
    exit(1);
  }

  matrix->B = malloc(sizeof(double complex)*matrix->n*matrix->n);
  if (matrix->A == NULL) {
    fprintf(stderr, "Unable to allocate memory for B.\n");
    exit(1);
  }

  //Zero out all of the elements of A and B.
  for (int i=0; i < matrix->n*matrix->n; i++) {
    matrix->A[i]=0.0;
    matrix->B[i]=0.0;
  }

  // First we'll fill out the matrix for the beta_r terms in the fluid
  for (int i=grid->is; i <= grid->ie; i++) {
    // beta_r_i^(n+1) terms
    wAelbg(5*(i) + 0, 5*(i) + 0, matrix,
	   -eta*grid->diffuse[i] 
	   - I*params->m*rotation->omega[i]);
    

    // beta_r_(i+1)^(n+1) terms
    wAelbg(5*(i) + 0, 5*(i+1) + 0, matrix, eta*grid->diffuse2[i]);

    // beta_r_(i-1)^(n+1) terms
    wAelbg(5*(i) + 0, 5*(i-1) + 0, matrix, eta*grid->diffuse2[i]);
    // beta_theta_i^(n+1) terms
    wAelbg(5*(i) + 0, 5*(i) + 1, matrix, -I*2*eta*params->m
	  *grid->r2inv[i]);

    // phi_r_i^(n+1) terms
    wAelbg(5*(i) + 0, 5*(i) + 2, matrix, kva);
  }

  // Next the beta_theta terms
  for (int i=grid->is; i <= grid->ie; i++) {
    // beta_theta_i^(n+1) terms
    wAelbg(5*(i) + 1, 5*(i) + 1, matrix,
	   -eta*grid->diffuse[i] 
	   - I*params->m*rotation->omega[i]);

    // beta_theta_(i+1)^(n+1) terms
    wAelbg(5*(i) + 1, 5*(i+1) + 1, matrix, eta*grid->diffuse2[i]);
    
    // beta_theta_(i-1)^(n+1) terms
    wAelbg(5*(i) + 1, 5*(i-1) + 1, matrix, eta*grid->diffuse2[i]);

    //beta_r_i^(n+1) terms
    wAelbg(5*(i) + 1, 5*(i) + 0, matrix,
	   0.5*(rotation->omega[i+1] -
		       rotation->omega[i-1])/grid->dx
	   + I*2.0*eta*params->m*grid->r2inv[i]);

    //phi_theta_i^(n+1) terms
    wAelbg(5*(i) + 1, 5*(i) + 3, matrix, kva);
  }

  // Next the phi_r terms
  for (int i=grid->is; i <= grid->ie; i++) {
    // phi_r_i^(n+1) terms
    wAelbg(5*(i) + 2, 5*(i) + 2, matrix,
	   -nu*grid->diffuse[i] - I*params->m*rotation->omega[i]);

    // phi_r_(i+1)^(n+1) terms
    wAelbg(5*(i) + 2, 5*(i+1) + 2, matrix, nu*grid->diffuse2[i]);

    //phi_r_(i-1)^(n+1) terms
    wAelbg(5*(i) + 2, 5*(i-1) + 2, matrix, nu*grid->diffuse2[i]);

    //beta_r_i^(n+1) terms
    wAelbg(5*(i) + 2, 5*(i) + 0, matrix, -kva);
    
    //phi_theta_i^(n+1) terms
    wAelbg(5*(i) + 2, 5*(i) + 3, matrix, 2.0*
	  rotation->omega[i]
	  - I*2.0*nu*params->m*grid->r2inv[i]);

    //pi_(i+1)^(n+1) terms
    wAelbg(5*(i) + 2, 5*(i+1) + 4, matrix, -0.5/
	   (grid->dx*grid->r[i]));

    //pi_(i-1)^(n+1) terms
    wAelbg(5*(i) + 2, 5*(i-1) + 4, matrix, 0.5/
	   (grid->dx*grid->r[i]));
  }

  // Next the phi_theta terms
  for (int i=grid->is; i <= grid->ie; i++) {
    // phi_theta_i^(n+1) terms
    wAelbg(5*(i) + 3, 5*(i) + 3, matrix,
	  -nu*grid->diffuse[i] 
	  - I*params->m*rotation->omega[i]);

    // phi_theta_(i+1)^(n+1) terms
    wAelbg(5*(i) + 3, 5*(i+1) + 3, matrix, nu*grid->diffuse2[i]);

    //phi_theta_(i-1)^(n+1) terms
    wAelbg(5*(i) + 3, 5*(i-1) + 3, matrix, nu*grid->diffuse2[i]);

    //beta_theta_i^(n+1) terms
    wAelbg(5*(i) + 3, 5*(i) + 1, matrix, -kva);

    //phi_r_i^(n+1) terms
    wAelbg(5*(i) + 3, 5*(i) + 2, matrix,
	   -2.0*rotation->omega[i] 
	   - 0.5*(1/grid->dx)*
	   (rotation->omega[i+1]-rotation->omega[i-1])
	   + I*2.0*nu*params->m*grid->r2inv[i]);

    //pi_i^(n+1) terms
    wAelbg(5*(i) + 3, 5*(i) + 4, matrix, -I*params->m/grid->r[i]);
  }

  // Next the pi terms
  for (int i=grid->is; i <= grid->ie; i++) {
    
    scalefactorpi = nu*grid->diffuse[(grid->ie-grid->is)/2]/
      (2.0*dx2inv*grid->r2inv[i] + m2*grid->r2inv[i] + k2);
    //pi_i^(n+1) terms
    wAelbg(5*(i) + 4, 5*(i) + 4, matrix,
           -(2.0*dx2inv*grid->r2inv[i] + m2*grid->r2inv[i] + k2)
	   *scalefactorpi);

    //pi_(i+1)^(n+1) terms
    wAelbg(5*(i) + 4, 5*(i+1) + 4, matrix, dx2inv*grid->r2inv[i]
	   *scalefactorpi);

    //pi_(i-1)^(n+1) terms
    wAelbg(5*(i) + 4, 5*(i-1) + 4, matrix, dx2inv*grid->r2inv[i]
	   *scalefactorpi);

    //phi_r_i^(n+1) terms
    wAelbg(5*(i) + 4, 5*(i) + 2, matrix,
	   (2.0*I*params->m/grid->r[i])*
	   (rotation->omega[i] + (0.5/grid->dx)*
	    (rotation->omega[i+1] - rotation->omega[i-1]))*scalefactorpi);
    
    //phi_theta_i^(n+1) terms
    wAelbg(5*(i) + 4, 5*(i) + 3, matrix,
    	   (-2.0*rotation->omega[i]/grid->r[i])*scalefactorpi);

    //phi_theta_(i+1)^(n+1) terms
    wAelbg(5*(i) + 4, 5*(i+1) + 3, matrix, 
	   (-rotation->omega[i+1]/(grid->r[i]*grid->dx))*scalefactorpi);

    //phi_theta_(i-1)^(n+2) terms
    wAelbg(5*(i) + 4, 5*(i-1) + 3, matrix,
	   (rotation->omega[i-1]/(grid->r[i]*grid->dx))*scalefactorpi);    
  }



  // Now we need to apply boundary conditions. */

  // Scale factors to make coefficients in constraint equations comparable
  // to coefficients in diffusion equation.
  scalefactorv = nu*grid->diffuse[(grid->ie-grid->is)/2];
  scalefactorb = eta*grid->diffuse[(grid->ie-grid->is)/2];

  if (params->magnetic_bc == 0) {
    // Perfectly-conducting boundary condition for beta_r */
    
    //beta_r_0 + beta_r_1 = 0
    wAelbg(5*(grid->is-1) + 0, 5*(grid->is) + 0, matrix, scalefactorb);
    wAelbg(5*(grid->is-1) + 0, 5*(grid->is-1) + 0, matrix, scalefactorb);

    //beta_r_N+1 + beta_r_N = 0
    wAelbg(5*(grid->ie+1), 5*(grid->ie+1), matrix, scalefactorb);
    wAelbg(5*(grid->ie+1), 5*(grid->ie), matrix, scalefactorb);
	   
    // Perfectly-conducting boundary condition for beta_theta
    
    facin = (1+0.5*grid->dx)/(1-0.5*grid->dx);
    facout = (1-0.5*grid->dx)/(1+0.5*grid->dx);

    //beta_theta_0 - facin*beta_theta_1 = 0
    wAelbg(5*(grid->is-1) + 1, 5*(grid->is) + 1, matrix, 
	   -facin*scalefactorb);
    wAelbg(5*(grid->is-1) + 1, 5*(grid->is-1) + 1, matrix, scalefactorb);

    //beta_theta_N+1 - facout*beta_theta_N = 0
    wAelbg(5*(grid->ie+1) + 1, 5*(grid->ie+1) + 1, matrix,
	   scalefactorb);
    wAelbg(5*(grid->ie+1) + 1, 5*(grid->ie) + 1, matrix,
	   -facout*scalefactorb);

  } else if(params->magnetic_bc == 1) {
    /*
    //beta_theta = 0 at inner boundary
    wAelbg(5*(grid->is-1) + 1, 5*(grid->is-1) + 1, matrix, scalefactorb);
    wAelbg(5*(grid->is-1) + 1, 5*(grid->is) + 1, matrix, scalefactorb);

    //beta_theta = 0 at outer boundary
    wAelbg(5*(grid->ie+1) + 1, 5*(grid->ie+1) + 1, matrix, scalefactorb);
    wAelbg(5*(grid->ie+1) + 1, 5*(grid->ie) + 1, matrix, scalefactorb);

    //Condition on beta_r
    Im = gsl_sf_bessel_In(params->m, params->k*params->r1);
    Imp1 = gsl_sf_bessel_In((params->m + 1), params->k*params->r1);

    Km = gsl_sf_bessel_Kn(params->m, params->k*params->r2);
    Kmp1 = gsl_sf_bessel_Kn((params->m + 1), params->k*params->r2);

    facin = 0.5*grid->dx
      *(params->r1*((k2 + m2/pow(params->r1, 2))/params->k)
	*Im/(Imp1 + (params->m/(params->k*params->r1))*Im) - 1.0);

    facout = 0.5*grid->dx
      *(params->r2*((k2 + m2/pow(params->r2, 2))/params->k)
	*Km/(-Kmp1 + (params->m/(params->k*params->r2))*Km) - 1.0);

    tempscalefactor = scalefactorb/facin;

    wAelbg(5*(grid->is-1) + 0, 5*(grid->is-1) + 0, matrix, 
	   tempscalefactor*(facin + 1.0));
    wAelbg(5*(grid->is-1) + 0, 5*(grid->is) + 0, matrix,
	   tempscalefactor*(facin - 1.0));
    
    tempscalefactor = scalefactorb/facout;

    wAelbg(5*(grid->ie+1) + 0, 5*(grid->ie+1) + 0, matrix,
	   tempscalefactor*(facout - 1.0));
    wAelbg(5*(grid->ie+1) + 0, 5*(grid->ie) + 0, matrix,
	   tempscalefactor*(facout + 1.0));
    */

    double constin1, constin2, constout1, constout2;

    //match beta_r to insulating solution
    constin1 = gsl_sf_bessel_In((params->m+1), params->k*grid->r[grid->is-1])
      + (params->m/(params->k*grid->r[grid->is-1]))*
      gsl_sf_bessel_In((params->m), params->k*grid->r[grid->is-1]);
    constin2 = gsl_sf_bessel_In((params->m+1), params->k*grid->r[grid->is])
      + (params->m/(params->k*grid->r[grid->is]))*
      gsl_sf_bessel_In((params->m), params->k*grid->r[grid->is]);
    constout1 = -gsl_sf_bessel_Kn((params->m+1), params->k*grid->r[grid->ie+1])
      + (params->m/(params->k*grid->r[grid->ie+1]))*
      gsl_sf_bessel_Kn((params->m), params->k*grid->r[grid->ie+1]);
    constout2 = -gsl_sf_bessel_Kn((params->m+1), params->k*grid->r[grid->ie])
      + (params->m/(params->k*grid->r[grid->ie]))*
      gsl_sf_bessel_Kn((params->m), params->k*grid->r[grid->ie]);

    wAelbg(5*(grid->is-1) + 0, 5*(grid->is-1) + 0, matrix, 
	   1.0);
    wAelbg(5*(grid->is-1) + 0, 5*(grid->is) + 0, matrix,
	   -constin1/constin2);
    
    wAelbg(5*(grid->ie+1) + 0, 5*(grid->ie+1) + 0, matrix,
	   1.0);
    wAelbg(5*(grid->ie+1) + 0, 5*(grid->ie) + 0, matrix,
	   -constout1/constout2);

    //beta_theta boundary condition
    if(params->m == 0) {
      //beta_theta = 0 at inner boundary
      wAelbg(5*(grid->is-1) + 1, 5*(grid->is-1) + 1, matrix, 1.0);
      wAelbg(5*(grid->is-1) + 1, 5*(grid->is) + 1, matrix, 1.0);
      
      //beta_theta = 0 at outer boundary
      wAelbg(5*(grid->ie+1) + 1, 5*(grid->ie+1) + 1, matrix, 1.0);
      wAelbg(5*(grid->ie+1) + 1, 5*(grid->ie) + 1, matrix, 1.0);
    } else {
      //match to insulating solution
      constin1 = gsl_sf_bessel_In((params->m), params->k*grid->r[grid->is-1])
	*grid->r[grid->is];
      constin2 = gsl_sf_bessel_In((params->m), params->k*grid->r[grid->is])
	*grid->r[grid->is-1];
      constout1 = gsl_sf_bessel_Kn((params->m), params->k*grid->r[grid->ie+1])
	*grid->r[grid->ie];
      constout2 = gsl_sf_bessel_Kn((params->m), params->k*grid->r[grid->ie])
	*grid->r[grid->ie+1];

      wAelbg(5*(grid->is-1) + 1, 5*(grid->is-1) + 1, matrix, 1.0);
      wAelbg(5*(grid->is-1) + 1, 5*(grid->is) + 1, matrix,
	     -constin1/constin2);
      
      wAelbg(5*(grid->ie+1) + 1, 5*(grid->ie+1) + 1, matrix, 1.0);
      wAelbg(5*(grid->ie+1) + 1, 5*(grid->ie) + 1, matrix,
	     -constout1/constout2);
    }

  }

  // No-slip boundary condition for phi_r, phi_theta  

  //phi_theta = 0 at boundary by setting phi_theta_0 + phi_theta_1 = 0
  wAelbg(5*(grid->is-1) + 3, 5*(grid->is-1) + 3, matrix, scalefactorv);
  wAelbg(5*(grid->is-1) + 3, 5*(grid->is) + 3, matrix, scalefactorv);

  //phi_theta_N+1=0, phi_theta_N = 0
  wAelbg(5*(grid->ie+1) + 3, 5*(grid->ie+1) + 3, matrix, scalefactorv);
  wAelbg(5*(grid->ie+1) + 3, 5*(grid->ie) + 3, matrix, scalefactorv);


  //phi_r = 0 = dphi_r/dr at boundary. Set points on both sides of boundary
  //to zero
  wAelbg(5*(grid->is-1) + 2, 5*(grid->is-1) + 2, matrix, scalefactorv);
  //Note that this equation is actually a constraint on dphi_r just inside
  //the boundary.  It is *not* an equation for Pi.
  wAelbg(5*(grid->is-1) + 4, 5*(grid->is) + 2, matrix, scalefactorv);
  
  wAelbg(5*(grid->ie+1) + 2, 5*(grid->ie+1) + 2, matrix, scalefactorv);
  wAelbg(5*(grid->ie+1) + 4, 5*(grid->ie) + 2, matrix, scalefactorv);


  //Now write 1s on the diagonal for the non-constraint equations
  for (int i=grid->is; i <= grid->ie; i++) {  
    wBelbg((5*(i) + 0) , (5*(i) + 0), matrix, 1.0);
    wBelbg((5*(i) + 1) , (5*(i) + 1), matrix, 1.0);
    wBelbg((5*(i) + 2) , (5*(i) + 2), matrix, 1.0);
    wBelbg((5*(i) + 3) , (5*(i) + 3), matrix, 1.0);
  }

  //Now allocate space for all of the LAPACK variables

  results = malloc(sizeof(RESULTS_STRUCT));
  assert(results);
  results->lambda = malloc(sizeof(double complex)*matrix->n);
  assert(results->lambda);
  results->z = malloc(sizeof(double complex)*matrix->n*matrix->n);
  assert(results->z);
  results->residual = malloc(sizeof(double)*matrix->n);
  assert(results->residual);

  //Make sure these things don't end up full of junk.
  for (int i=0; i < matrix->n; i++) {
    results->lambda[i]=0.0;
    results->residual[i]=0.0;
  }


  alpha = malloc(sizeof(double complex)*matrix->n);
  assert(alpha);
  beta = malloc(sizeof(double complex)*matrix->n);
  assert(beta);
  work = malloc(sizeof(double complex)*8*matrix->n);
  assert(work);
  rwork = malloc(sizeof(double)*8*matrix->n);
  assert(rwork);
  
  jobvl = 'N';
  jobvr = 'V';
  ldvl = 1;
  ldvr = matrix->n;
  lwork = -1;

  //Run the routine once to have it return the optimum size of the work array
  zggev_(&jobvl, &jobvr, &matrix->n, matrix->A, &matrix->n, matrix->B,
	 &matrix->n, alpha, beta, NULL, &ldvl, results->z, &ldvr, work,
	 &lwork, rwork, &info);

  lwork = creal(work[0]);
  free(work);
  work = malloc(sizeof(double complex)*lwork);
  assert(work);

  //Now run again to actually solve the eigenproblem.
  zggev_(&jobvl, &jobvr, &matrix->n, matrix->A, &matrix->n, matrix->B,
	 &matrix->n, alpha, beta, NULL, &ldvl, results->z, &ldvr, work,
	 &lwork, rwork, &info);
  
  if (info != 0) {
    printf("Error: ZGGEV returned %i.\n", info);
  }

  //Now store the eigenvalues from the calculation
  for (int i = 0; i < matrix->n; i++) {
    results->lambda[i] = alpha[i]/beta[i];
  }

  //Set the number of converged eigenvalues to matrix->n, since we solve
  //for everything here.

  results->nconv = matrix->n;

  //And lets get out of here!
  free(alpha);
  free(beta);
  free(vr);
  free(work);
  free(rwork);
  return results;
}


void wAelbg(int i, int j, COMPRESSED_MATRIX *matrix, double complex value) {
  //We are writing to the FULL SIZED UNCOMPRESSED MATRIX here.
  //We simply find the position of the elemtn in the *column-major* matrix
  //structure, and write the value.

  int k;
  k  matrix->n*j + i;
  matrix->A[k] = value;
}


void wBelbg(int i, int j, COMPRESSED_MATRIX *matrix, double complex value) {
  //We are writing to the FULL SIZED UNCOMPRESSED MATRIX here.
  //We simply find the position of the elemtn in the *column-major* matrix
  //structure, and write the value.

  int k;
  k  matrix->n*j + i;
  matrix->B[k] = value;
}
