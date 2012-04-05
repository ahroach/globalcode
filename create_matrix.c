#include <stdio.h>
#include <stdlib.h>
#include "global.h"

COMPRESSED_MATRIX *create_matrix(int numelems) {
  COMPRESSED_MATRIX *matrix;
  matrix = malloc(sizeof(COMPRESSED_MATRIX));

  if (matrix == NULL) {
    fprintf(stderr, "Unable to allocate memory for matrix structure.\n");
    exit(1);
  }

  /* I need to generate a matrix containing the coefficients of the
   * terms in the finite difference equation.  The order of the terms
   * will be beta_r, beta_theta, phi_r, phi_theta, pi at each
   * gridpoint.
   */

  /* The matrix to hold these values will have 5*numcells rows and
   * 5*numcells columns.  These need to be complex numbers, so there
   * will be 2 doubles stored for each value
   */


  matrix->n = numelems; //Number of columns in A
  matrix->m = numelems; //Number of rows in A
  matrix->kl = 7; //Number of subdiagonals within the band of A
  matrix->ku = 7; //Number of superdiagonals within the band of A
  matrix->lda = 2*matrix->kl + matrix->ku + 1; //Leading dimension of array A
  matrix->ldab = matrix->kl + matrix->ku + 1; //Leading dimension of array Bb

  matrix->A = malloc(sizeof(double complex)*matrix->n*matrix->lda);
  if (matrix->A == NULL) {
    fprintf(stderr, "Unable to allocate memory for A.\n");
    exit(1);
  }

  matrix->B = malloc(sizeof(double complex)*matrix->n*matrix->lda);
  if (matrix->B == NULL) {
    fprintf(stderr, "Unable to allocate memory for B.\n");
    exit(1);
  }

  matrix->Bb = malloc(sizeof(double complex)*matrix->n*matrix->ldab);
  if (matrix->Bb == NULL) {
    fprintf(stderr, "Unable to allocate memory for Bb.\n");
    exit(1);
  }

  return matrix;
}
