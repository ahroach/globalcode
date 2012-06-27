/* Step the equations of the global code
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

RESULTS_STRUCT *eigensolve(COMPRESSED_MATRIX *matrix, 
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

  //Things needed for the ARPACK computation
  double complex *OPfac;
  int info;
  int ipiv[matrix->m]; //Integer array of pivot indices, returned by zgbtrf
  int ido;
  char bmat;
  int bmat_len = 1;
  int which_len = 2;
  double complex *resid;
  int ncv;
  double complex *v;
  int iparam[11];
  int ipntr[14];
  int lworkl;
  double complex *workd;
  double complex *workl;
  double *rwork;
  char trans;
  int trans_len;
  double complex alpha;
  double complex beta;
  int incx;
  int incy;
  int nrhs;
  int ierr;
  int rvec;
  char howmny;
  int howmny_len;
  int *select;
  RESULTS_STRUCT *results;
  double complex *workev;

  //Zero out all of the elements of A.
  for (int i=0; i < matrix->n*matrix->lda; i++) {
    matrix->A[i]=0.0;
  }

  // First we'll fill out the matrix for the beta_r terms in the fluid
  for (int i=grid->is; i <= grid->ie; i++) {
    // beta_r_j terms
    wAelem(5*(i) + 0, 5*(i) + 0, matrix,
	   -eta*grid->diffuse[i] 
	   - I*params->m*rotation->omega[i]);
    

    // beta_r_(j+1) terms
    wAelem(5*(i) + 0, 5*(i+1) + 0, matrix, eta*grid->diffuse2[i]);

    // beta_r_(j-1) terms
    wAelem(5*(i) + 0, 5*(i-1) + 0, matrix, eta*grid->diffuse2[i]);
    // beta_theta_j terms
    wAelem(5*(i) + 0, 5*(i) + 1, matrix, -I*2*eta*params->m
	  *grid->r2inv[i]);

    // phi_r_j terms
    wAelem(5*(i) + 0, 5*(i) + 2, matrix, kva);
  }

  // Next the beta_theta terms
  for (int i=grid->is; i <= grid->ie; i++) {
    // beta_theta_j terms
    wAelem(5*(i) + 1, 5*(i) + 1, matrix,
	   -eta*grid->diffuse[i] 
	   - I*params->m*rotation->omega[i]);

    // beta_theta_(j+1) terms
    wAelem(5*(i) + 1, 5*(i+1) + 1, matrix, eta*grid->diffuse2[i]);
    
    // beta_theta_(j-1) terms
    wAelem(5*(i) + 1, 5*(i-1) + 1, matrix, eta*grid->diffuse2[i]);

    //beta_r_j terms
    wAelem(5*(i) + 1, 5*(i) + 0, matrix,
	   0.5*(rotation->omega[i+1] -
		       rotation->omega[i-1])/grid->dx
	   + I*2.0*eta*params->m*grid->r2inv[i]);

    //phi_theta_j terms
    wAelem(5*(i) + 1, 5*(i) + 3, matrix, kva);
  }

  // Next the phi_r terms
  for (int i=grid->is; i <= grid->ie; i++) {
    // phi_r_j terms
    wAelem(5*(i) + 2, 5*(i) + 2, matrix,
	   -nu*grid->diffuse[i] - I*params->m*rotation->omega[i]);

    // phi_r_(j+1) terms
    wAelem(5*(i) + 2, 5*(i+1) + 2, matrix, nu*grid->diffuse2[i]);

    //phi_r_(j-1) terms
    wAelem(5*(i) + 2, 5*(i-1) + 2, matrix, nu*grid->diffuse2[i]);

    //beta_r_j terms
    wAelem(5*(i) + 2, 5*(i) + 0, matrix, -kva);
    
    //phi_theta_j terms
    wAelem(5*(i) + 2, 5*(i) + 3, matrix, 2.0*
	  rotation->omega[i]
	  - I*2.0*nu*params->m*grid->r2inv[i]);

    //pi_(j+1) terms
    wAelem(5*(i) + 2, 5*(i+1) + 4, matrix, -0.5/
	   (grid->dx*grid->r[i]));

    //pi_(j-1) terms
    wAelem(5*(i) + 2, 5*(i-1) + 4, matrix, 0.5/
	   (grid->dx*grid->r[i]));
  }

  // Next the phi_theta terms
  for (int i=grid->is; i <= grid->ie; i++) {
    // phi_theta_j terms
    wAelem(5*(i) + 3, 5*(i) + 3, matrix,
	  -nu*grid->diffuse[i] 
	  - I*params->m*rotation->omega[i]);

    // phi_theta_(j+1) terms
    wAelem(5*(i) + 3, 5*(i+1) + 3, matrix, nu*grid->diffuse2[i]);

    //phi_theta_(j-1) terms
    wAelem(5*(i) + 3, 5*(i-1) + 3, matrix, nu*grid->diffuse2[i]);

    //beta_theta_j terms
    wAelem(5*(i) + 3, 5*(i) + 1, matrix, -kva);

    //phi_r_j terms
    wAelem(5*(i) + 3, 5*(i) + 2, matrix,
	   -2.0*rotation->omega[i] 
	   - 0.5*(1/grid->dx)*
	   (rotation->omega[i+1]-rotation->omega[i-1])
	   + I*2.0*nu*params->m*grid->r2inv[i]);

    //pi_j terms
    wAelem(5*(i) + 3, 5*(i) + 4, matrix, -I*params->m/grid->r[i]);
  }

  // Next the pi terms
  for (int i=grid->is; i <= grid->ie; i++) {
    //pi_j terms
    wAelem(5*(i) + 4, 5*(i) + 4, matrix,
           -(2.0*dx2inv*grid->r2inv[i] + m2*grid->r2inv[i] + k2));

    //pi_(j+1) terms
    wAelem(5*(i) + 4, 5*(i+1) + 4, matrix, dx2inv*grid->r2inv[i]);

    //pi_(j-1) terms
    wAelem(5*(i) + 4, 5*(i-1) + 4, matrix, dx2inv*grid->r2inv[i]);

    //phi_r_j terms
    wAelem(5*(i) + 4, 5*(i) + 2, matrix,
	   (2.0*I*params->m/grid->r[i])*
	   (rotation->omega[i] + (0.5/grid->dx)*
	    (rotation->omega[i+1] - rotation->omega[i-1])));
    
    //phi_theta_j terms
    wAelem(5*(i) + 4, 5*(i) + 3, matrix,
    	   -2.0*rotation->omega[i]/grid->r[i]);

    //phi_theta_(j+1) terms
    wAelem(5*(i) + 4, 5*(i+1) + 3, matrix, 
	   -rotation->omega[i+1]/(grid->r[i]*grid->dx));

    //phi_theta_(j-1) terms
    wAelem(5*(i) + 4, 5*(i-1) + 3, matrix,
	   rotation->omega[i-1]/(grid->r[i]*grid->dx));    
  }



  // Now we need to apply boundary conditions. */

  if (params->magnetic_bc == 0) {
    // Perfectly-conducting boundary condition for beta_r */
    
    //beta_r_0 + beta_r_1 = 0
    wAelem(5*(grid->is-1) + 0, 5*(grid->is) + 0, matrix, 1.0);
    wAelem(5*(grid->is-1) + 0, 5*(grid->is-1) + 0, matrix, 1.0);

    //beta_r_N+1 + beta_r_N = 0
    wAelem(5*(grid->ie+1), 5*(grid->ie+1), matrix, 1.0);
    wAelem(5*(grid->ie+1), 5*(grid->ie), matrix, 1.0);
	   
    // Perfectly-conducting boundary condition for beta_theta
    
    facin = (1+0.5*grid->dx)/(1-0.5*grid->dx);
    facout = (1-0.5*grid->dx)/(1+0.5*grid->dx);

    //beta_theta_0 - facin*beta_theta_1 = 0
    wAelem(5*(grid->is-1) + 1, 5*(grid->is) + 1, matrix, -facin);
    wAelem(5*(grid->is-1) + 1, 5*(grid->is-1) + 1, matrix, 1.0);

    //beta_theta_N+1 - facout*beta_theta_N = 0
    wAelem(5*(grid->ie+1) + 1, 5*(grid->ie+1) + 1, matrix, 1.0);
    wAelem(5*(grid->ie+1) + 1, 5*(grid->ie) + 1, matrix, -facout);

  } else if(params->magnetic_bc == 1) {
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
  wAelem(5*(grid->is-1) + 3, 5*(grid->is-1) + 3, matrix, 1.0);
  wAelem(5*(grid->is-1) + 3, 5*(grid->is) + 3, matrix, 1.0);

  //phi_theta_N+1=0, phi_theta_N = 0
  wAelem(5*(grid->ie+1) + 3, 5*(grid->ie+1) + 3, matrix, 1.0);
  wAelem(5*(grid->ie+1) + 3, 5*(grid->ie) + 3, matrix, 1.0);


  //phi_r = 0 = dphi_r/dr at boundary. Set points on both sides of boundary
  //to zero
  wAelem(5*(grid->is-1) + 2, 5*(grid->is-1) + 2, matrix, 1.0);
  //Note that this equation is actually a constraint on dphi_r just inside
  //the boundary.  It is *not* an equation for Pi.
  wAelem(5*(grid->is-1) + 4, 5*(grid->is) + 2, matrix, 1.0);
  
  wAelem(5*(grid->ie+1) + 2, 5*(grid->ie+1) + 2, matrix, 1.0);
  wAelem(5*(grid->ie+1) + 4, 5*(grid->ie) + 2, matrix, 1.0);



  //Zero out all of the elements of B and Bb.
  for (int i=0; i < matrix->n*matrix->lda; i++) {
    matrix->B[i]=0.0;
  }

  for (int i=0; i < matrix->n*matrix->ldab; i++) {
    matrix->Bb[i]=0.0;
  }


  //Now write 1s on the diagonal for the non-constraint equations
  for (int i=grid->is; i <= grid->ie; i++) {  
    wBelem((5*(i) + 0) , (5*(i) + 0), matrix, 1.0);
    wBelem((5*(i) + 1) , (5*(i) + 1), matrix, 1.0);
    wBelem((5*(i) + 2) , (5*(i) + 2), matrix, 1.0);
    wBelem((5*(i) + 3) , (5*(i) + 3), matrix, 1.0);
  }


  //Set up a copy of OP to be factored by LU decomposition.
  OPfac = malloc(sizeof(double complex)*matrix->n*matrix->lda);  

  //OP = A - sigma*M
  for (int i = 0; i < matrix->n*matrix->lda; i++) {
    OPfac[i] = matrix->A[i] - arpack_params->sigma*matrix->B[i];
  }

  //Compute LU factorization of OP
  zgbtrf_(&matrix->m, &matrix->n, &matrix->kl, &matrix->ku, OPfac,
	  &matrix->lda, ipiv, &info);  
  if (info != 0) {
    fprintf(stderr, "Error: ZGBTRF returned %i.\n", info);
    return 1;
  }

  //Alright, now time for the main loop.
  
  ido = 0;
  bmat = 'G';
  resid = malloc(sizeof(double complex)*matrix->n);
  //This value of ncv based on a suggestion in the znaupd documentation
  ncv = arpack_params->nummodes*2 + 4;
  lworkl = 3*ncv*ncv + 5*ncv + 1;
  v = malloc(sizeof(double complex)*matrix->n*ncv);
  
  //Exact shifts
  iparam[0] = 1;
  //Max number of Arnoldi update iterations
  iparam[2] = arpack_params->maxiters;
  //Blocksize to be used in the recurrence
  iparam[3] = 1;
  //Generalized eigenvalue problem, B symmetric semi-definite
  iparam[6] = 3;

  workd = malloc(sizeof(double complex)*3*matrix->n);
  workl = malloc(sizeof(double complex)*lworkl);
  rwork = malloc(sizeof(double)*ncv);
  
  info = 0;

  while (1) {
    //Note this crazy requirement to pass the length of strings after
    //the string when calling FORTRAN subroutines.
    znaupd_(&ido, &bmat, &matrix->n, arpack_params->which,
	    &arpack_params->nummodes, &arpack_params->tol, resid, &ncv, v,
	    &matrix->n, iparam, ipntr, workd, workl, &lworkl, rwork, &info,
	    &bmat_len, &which_len); 

    if (info != 0) {
      fprintf(stderr, "Error in znaupd.  Info = %i.\n", info);
    }

    if (ido == -1) {
      //Perform y <-- OP*x = inv[A-sigma*B]*B*x
      trans = 'N';
      trans_len = 1;
      alpha = 1.0;
      incx = 1;
      beta = 0;
      incy = 1;
      nrhs = 1;
      zgbmv_(&trans, &matrix->n, &matrix->n, &matrix->kl, &matrix->ku,
	     &alpha, matrix->Bb, &matrix->ldab, &workd[ipntr[0]-1],
	     &incx, &beta, &workd[ipntr[1]-1], &incy, &trans_len);
      zgbtrs_(&trans, &matrix->n, &matrix->kl, &matrix->ku, &nrhs, OPfac,
	      &matrix->lda, ipiv, &workd[ipntr[1]-1], &matrix->n, &ierr,
	      &trans_len);
      if (ierr != 0) {
	fprintf(stderr, "Error %i in zgbtrs\n", ierr);
      }
    } else if (ido == 1) {
      //Perform y <- inv(A-sigma*B)*(B*x)
      //(B*X) has been computed and stored in workd[ipntr[2]-1]
      trans = 'N';
      trans_len = 1;
      nrhs = 1;
      memcpy(&workd[ipntr[1]-1], &workd[ipntr[2]-1],
	     sizeof(double complex)*matrix->n);
      zgbtrs_(&trans, &matrix->n, &matrix->kl, &matrix->ku, &nrhs, OPfac,
	      &matrix->lda, ipiv, &workd[ipntr[1]-1], &matrix->n, &ierr);
      if (ierr != 0) {
	fprintf(stderr, "Error %i in zgbtrs\n", ierr);
      }
    } else if (ido == 2) {
      //Perform y <- B*x
      trans = 'N';
      trans_len = 1;
      alpha = 1.0;
      incx = 1;
      beta = 0;
      incy = 1;
      zgbmv_(&trans, &matrix->n, &matrix->n, &matrix->kl, &matrix->ku,
	     &alpha, matrix->Bb, &matrix->ldab, &workd[ipntr[0]-1],
	     &incx, &beta, &workd[ipntr[1]-1], &incy, &trans_len);
    } else {
      break;
    }
  }
  
  results = malloc(sizeof(RESULTS_STRUCT));
  results->lambda = malloc(sizeof(double complex)
			   *(arpack_params->nummodes+1));
  results->z = malloc(sizeof(double complex)
		      *arpack_params->nummodes*matrix->n);
  results->residual = malloc(sizeof(double)*(arpack_params->nummodes + 1));


  //Look at some information from this run
  results->itersused = iparam[2];
  results->nconv = iparam[4];
  if(params->VERBOSE) {
    fprintf(stdout, "Run completed.\n%i Arnoldi update iterations taken of %i maximum.\n", results->itersused, arpack_params->maxiters);
    fprintf(stdout, "%i converged Ritz values, of %i requested.\n",
	    results->nconv, arpack_params->nummodes);
  }
  

  // Either there is convergence, or there is an error
  if (!(results->nconv > 0)){
    fprintf(stderr, "Error in znaupd.  info = %i.  Aborting.\n", info);
  } else {
    rvec = 1;
    howmny = 'A';
    howmny_len = 1;
    select = malloc(sizeof(int)*ncv);
    workev = malloc(sizeof(double complex)*ncv*2);

    //Again, make note of the funny extra length parameters for the
    //character*n
    zneupd_(&rvec, &howmny, select, results->lambda, results->z, &matrix->n,
	    &arpack_params->sigma, workev, &bmat, &matrix->n,
	    arpack_params->which, &arpack_params->nummodes,
	    &arpack_params->tol, resid, &ncv, v, &matrix->n, iparam, ipntr,
	    workd, workl, &lworkl, rwork, &info, &howmny_len, &bmat_len,
	    &which_len);
    if (info != 0) {
      fprintf(stderr, "Error in zneupd.  Info = %i.\n", info);
    }

    double complex *tempmatrix;
    double complex *residual;
    double residual_score;
    tempmatrix = malloc(sizeof(double complex)*matrix->n*matrix->ldab);  
    residual = malloc(sizeof(double complex)*matrix->n);
    
    
    for (int k = 0; k < results->nconv; k++) {
      //Test the quality of each of these eigenvalues
      //Compute (A-lambda*B)*eigenvalue
      
      //Compute A - lambda*B
      //Keep in mind that A is in the banded form needed by zgbtrf, so
      //there are matrix->kl rows of offset at the top. COLUMN-MAJOR STORAGE!
      for (int i = 0; i < matrix->ldab; i++) {
	for (int j = 0; j < matrix->n; j++) {
	  tempmatrix[j*matrix->ldab + i] = 
	    matrix->A[j*matrix->lda + (i + matrix->kl)] 
	    - results->lambda[k]*matrix->B[j*matrix->ldab + i];
	}
      }
      
      trans='N';
      alpha = 1;
      beta = 0.0;
      incx = 1;
      incy = 1;
      trans_len = 1;
      
      zgbmv_(&trans, &matrix->n, &matrix->n, &matrix->kl, &matrix->ku,
	     &alpha, tempmatrix, &matrix->ldab,
	     &results->z[matrix->n*k],
	     &incx, &beta, residual, &incy, &trans_len);
      
      //Now sum everything, and find the average deviation
      residual_score = 0;
      for (int i = 0; i < matrix->n; i++) {
	residual_score = residual_score + cabs(residual[i]);
      }
      residual_score = residual_score/matrix->n;

      results->residual[k] = residual_score;
    }
    free(tempmatrix);
    free(select);
    free(workev);
    free(residual);
    
  }

  free(OPfac);
  free(resid);
  free(v);
  free(workd);
  free(workl);
  free(rwork);
  return results;
}
