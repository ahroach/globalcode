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

  //Zero out all of the elements of A.
  for (int i=0; i < matrix->n*matrix->lda; i++) {
    matrix->A[i]=0.0;
  }

  // First we'll fill out the matrix for the beta_r terms in the fluid
  for (int i=grid->is; i <= grid->ie; i++) {
    // beta_r_i^(n+1) terms
    wAelem(5*(i) + 0, 5*(i) + 0, matrix,
	   -eta*grid->diffuse[i] 
	   - I*params->m*rotation->omega[i]);
    

    // beta_r_(i+1)^(n+1) terms
    wAelem(5*(i) + 0, 5*(i+1) + 0, matrix, eta*grid->diffuse2[i]);

    // beta_r_(i-1)^(n+1) terms
    wAelem(5*(i) + 0, 5*(i-1) + 0, matrix, eta*grid->diffuse2[i]);
    // beta_theta_i^(n+1) terms
    wAelem(5*(i) + 0, 5*(i) + 1, matrix, -I*2*eta*params->m
	  *grid->r2inv[i]);

    // phi_r_i^(n+1) terms
    wAelem(5*(i) + 0, 5*(i) + 2, matrix, kva);
  }

  // Next the beta_theta terms
  for (int i=grid->is; i <= grid->ie; i++) {
    // beta_theta_i^(n+1) terms
    wAelem(5*(i) + 1, 5*(i) + 1, matrix,
	   -eta*grid->diffuse[i] 
	   - I*params->m*rotation->omega[i]);

    // beta_theta_(i+1)^(n+1) terms
    wAelem(5*(i) + 1, 5*(i+1) + 1, matrix, eta*grid->diffuse2[i]);
    
    // beta_theta_(i-1)^(n+1) terms
    wAelem(5*(i) + 1, 5*(i-1) + 1, matrix, eta*grid->diffuse2[i]);

    //beta_r_i^(n+1) terms
    wAelem(5*(i) + 1, 5*(i) + 0, matrix,
	   0.5*(rotation->omega[i+1] -
		       rotation->omega[i-1])/grid->dx
	   + I*2.0*eta*params->m*grid->r2inv[i]);

    //phi_theta_i^(n+1) terms
    wAelem(5*(i) + 1, 5*(i) + 3, matrix, kva);
  }

  // Next the phi_r terms
  for (int i=grid->is; i <= grid->ie; i++) {
    // phi_r_i^(n+1) terms
    wAelem(5*(i) + 2, 5*(i) + 2, matrix,
	   -nu*grid->diffuse[i] - I*params->m*rotation->omega[i]);

    // phi_r_(i+1)^(n+1) terms
    wAelem(5*(i) + 2, 5*(i+1) + 2, matrix, nu*grid->diffuse2[i]);

    //phi_r_(i-1)^(n+1) terms
    wAelem(5*(i) + 2, 5*(i-1) + 2, matrix, nu*grid->diffuse2[i]);

    //beta_r_i^(n+1) terms
    wAelem(5*(i) + 2, 5*(i) + 0, matrix, -kva);
    
    //phi_theta_i^(n+1) terms
    wAelem(5*(i) + 2, 5*(i) + 3, matrix, 2.0*
	  rotation->omega[i]
	  - I*2.0*nu*params->m*grid->r2inv[i]);

    //pi_(i+1)^(n+1) terms
    wAelem(5*(i) + 2, 5*(i+1) + 4, matrix, -0.5/
	   (grid->dx*grid->r[i]));

    //pi_(i-1)^(n+1) terms
    wAelem(5*(i) + 2, 5*(i-1) + 4, matrix, 0.5/
	   (grid->dx*grid->r[i]));
  }

  // Next the phi_theta terms
  for (int i=grid->is; i <= grid->ie; i++) {
    // phi_theta_i^(n+1) terms
    wAelem(5*(i) + 3, 5*(i) + 3, matrix,
	  -nu*grid->diffuse[i] 
	  - I*params->m*rotation->omega[i]);

    // phi_theta_(i+1)^(n+1) terms
    wAelem(5*(i) + 3, 5*(i+1) + 3, matrix, nu*grid->diffuse2[i]);

    //phi_theta_(i-1)^(n+1) terms
    wAelem(5*(i) + 3, 5*(i-1) + 3, matrix, nu*grid->diffuse2[i]);

    //beta_theta_i^(n+1) terms
    wAelem(5*(i) + 3, 5*(i) + 1, matrix, -kva);

    //phi_r_i^(n+1) terms
    wAelem(5*(i) + 3, 5*(i) + 2, matrix,
	   -2.0*rotation->omega[i] 
	   - 0.5*(1/grid->dx)*
	   (rotation->omega[i+1]-rotation->omega[i-1])
	   + I*2.0*nu*params->m*grid->r2inv[i]);

    //pi_i^(n+1) terms
    wAelem(5*(i) + 3, 5*(i) + 4, matrix, -I*params->m/grid->r[i]);
  }

  // Next the pi terms
  for (int i=grid->is; i <= grid->ie; i++) {
    
    scalefactorpi = nu*grid->diffuse[(grid->ie-grid->is)/2]/
      (2.0*dx2inv*grid->r2inv[i] + m2*grid->r2inv[i] + k2);
    //pi_i^(n+1) terms
    wAelem(5*(i) + 4, 5*(i) + 4, matrix,
           -(2.0*dx2inv*grid->r2inv[i] + m2*grid->r2inv[i] + k2)
	   *scalefactorpi);

    //pi_(i+1)^(n+1) terms
    wAelem(5*(i) + 4, 5*(i+1) + 4, matrix, dx2inv*grid->r2inv[i]
	   *scalefactorpi);

    //pi_(i-1)^(n+1) terms
    wAelem(5*(i) + 4, 5*(i-1) + 4, matrix, dx2inv*grid->r2inv[i]
	   *scalefactorpi);

    //phi_r_i^(n+1) terms
    wAelem(5*(i) + 4, 5*(i) + 2, matrix,
	   (2.0*I*params->m/grid->r[i])*
	   (rotation->omega[i] + (0.5/grid->dx)*
	    (rotation->omega[i+1] - rotation->omega[i-1]))*scalefactorpi);
    
    //phi_theta_i^(n+1) terms
    wAelem(5*(i) + 4, 5*(i) + 3, matrix,
    	   (-2.0*rotation->omega[i]/grid->r[i])*scalefactorpi);

    //phi_theta_(i+1)^(n+1) terms
    wAelem(5*(i) + 4, 5*(i+1) + 3, matrix, 
	   (-rotation->omega[i+1]/(grid->r[i]*grid->dx))*scalefactorpi);

    //phi_theta_(i-1)^(n+2) terms
    wAelem(5*(i) + 4, 5*(i-1) + 3, matrix,
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
    wAelem(5*(grid->is-1) + 0, 5*(grid->is) + 0, matrix, scalefactorb);
    wAelem(5*(grid->is-1) + 0, 5*(grid->is-1) + 0, matrix, scalefactorb);

    //beta_r_N+1 + beta_r_N = 0
    wAelem(5*(grid->ie+1), 5*(grid->ie+1), matrix, scalefactorb);
    wAelem(5*(grid->ie+1), 5*(grid->ie), matrix, scalefactorb);
	   
    // Perfectly-conducting boundary condition for beta_theta
    
    facin = (1+0.5*grid->dx)/(1-0.5*grid->dx);
    facout = (1-0.5*grid->dx)/(1+0.5*grid->dx);

    //beta_theta_0 - facin*beta_theta_1 = 0
    wAelem(5*(grid->is-1) + 1, 5*(grid->is) + 1, matrix, 
	   -facin*scalefactorb);
    wAelem(5*(grid->is-1) + 1, 5*(grid->is-1) + 1, matrix, scalefactorb);

    //beta_theta_N+1 - facout*beta_theta_N = 0
    wAelem(5*(grid->ie+1) + 1, 5*(grid->ie+1) + 1, matrix,
	   scalefactorb);
    wAelem(5*(grid->ie+1) + 1, 5*(grid->ie) + 1, matrix,
	   -facout*scalefactorb);

  } else if(params->magnetic_bc == 1) {
    /*
    //beta_theta = 0 at inner boundary
    wAelem(5*(grid->is-1) + 1, 5*(grid->is-1) + 1, matrix, scalefactorb);
    wAelem(5*(grid->is-1) + 1, 5*(grid->is) + 1, matrix, scalefactorb);

    //beta_theta = 0 at outer boundary
    wAelem(5*(grid->ie+1) + 1, 5*(grid->ie+1) + 1, matrix, scalefactorb);
    wAelem(5*(grid->ie+1) + 1, 5*(grid->ie) + 1, matrix, scalefactorb);

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

    wAelem(5*(grid->is-1) + 0, 5*(grid->is-1) + 0, matrix, 
	   tempscalefactor*(facin + 1.0));
    wAelem(5*(grid->is-1) + 0, 5*(grid->is) + 0, matrix,
	   tempscalefactor*(facin - 1.0));
    
    tempscalefactor = scalefactorb/facout;

    wAelem(5*(grid->ie+1) + 0, 5*(grid->ie+1) + 0, matrix,
	   tempscalefactor*(facout - 1.0));
    wAelem(5*(grid->ie+1) + 0, 5*(grid->ie) + 0, matrix,
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

    wAelem(5*(grid->is-1) + 0, 5*(grid->is-1) + 0, matrix, 
	   1.0);
    wAelem(5*(grid->is-1) + 0, 5*(grid->is) + 0, matrix,
	   -constin1/constin2);
    
    wAelem(5*(grid->ie+1) + 0, 5*(grid->ie+1) + 0, matrix,
	   1.0);
    wAelem(5*(grid->ie+1) + 0, 5*(grid->ie) + 0, matrix,
	   -constout1/constout2);

    //beta_theta boundary condition
    if(params->m == 0) {
      //beta_theta = 0 at inner boundary
      wAelem(5*(grid->is-1) + 1, 5*(grid->is-1) + 1, matrix, 1.0);
      wAelem(5*(grid->is-1) + 1, 5*(grid->is) + 1, matrix, 1.0);
      
      //beta_theta = 0 at outer boundary
      wAelem(5*(grid->ie+1) + 1, 5*(grid->ie+1) + 1, matrix, 1.0);
      wAelem(5*(grid->ie+1) + 1, 5*(grid->ie) + 1, matrix, 1.0);
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

      wAelem(5*(grid->is-1) + 1, 5*(grid->is-1) + 1, matrix, 1.0);
      wAelem(5*(grid->is-1) + 1, 5*(grid->is) + 1, matrix,
	     -constin1/constin2);
      
      wAelem(5*(grid->ie+1) + 1, 5*(grid->ie+1) + 1, matrix, 1.0);
      wAelem(5*(grid->ie+1) + 1, 5*(grid->ie) + 1, matrix,
	     -constout1/constout2);
    }

  }

  // No-slip boundary condition for phi_r, phi_theta  

  //phi_theta = 0 at boundary by setting phi_theta_0 + phi_theta_1 = 0
  wAelem(5*(grid->is-1) + 3, 5*(grid->is-1) + 3, matrix, scalefactorv);
  wAelem(5*(grid->is-1) + 3, 5*(grid->is) + 3, matrix, scalefactorv);

  //phi_theta_N+1=0, phi_theta_N = 0
  wAelem(5*(grid->ie+1) + 3, 5*(grid->ie+1) + 3, matrix, scalefactorv);
  wAelem(5*(grid->ie+1) + 3, 5*(grid->ie) + 3, matrix, scalefactorv);


  //phi_r = 0 = dphi_r/dr at boundary. Set points on both sides of boundary
  //to zero
  wAelem(5*(grid->is-1) + 2, 5*(grid->is-1) + 2, matrix, scalefactorv);
  //Note that this equation is actually a constraint on dphi_r just inside
  //the boundary.  It is *not* an equation for Pi.
  wAelem(5*(grid->is-1) + 4, 5*(grid->is) + 2, matrix, scalefactorv);
  
  wAelem(5*(grid->ie+1) + 2, 5*(grid->ie+1) + 2, matrix, scalefactorv);
  wAelem(5*(grid->ie+1) + 4, 5*(grid->ie) + 2, matrix, scalefactorv);



  //Zero out all of the elements of B and Bb.
  for (int i=0; i < matrix->n*matrix->lda; i++) {
    matrix->B[i]=0.0;
  }


  //Now write 1s on the diagonal for the non-constraint equations
  for (int i=grid->is; i <= grid->ie; i++) {  
    wBelem((5*(i) + 0) , (5*(i) + 0), matrix, 1.0);
    wBelem((5*(i) + 1) , (5*(i) + 1), matrix, 1.0);
    wBelem((5*(i) + 2) , (5*(i) + 2), matrix, 1.0);
    wBelem((5*(i) + 3) , (5*(i) + 3), matrix, 1.0);
  }


  results = malloc(sizeof(RESULTS_STRUCT));
  results->lambda = malloc(sizeof(double complex)
			   *(arpack_params->nummodes+1));
  results->z = malloc(sizeof(double complex)
		      *arpack_params->nummodes*matrix->n);
  results->residual = malloc(sizeof(double)*(arpack_params->nummodes + 1));

  return results;
}
