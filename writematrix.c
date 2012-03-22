/* Writes matrix elements of the compressed band-diagonal matrix A.  The
 * functions take the i and j coordinates of the elements in the uncompressed matrix, and
 * write these to the correct positions in the column-major matrix A.  This
 * matrix must be column-major because it is passed to fortran routines, which
 * expect it.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "global.h"

void wAelem(int i, int j, COMPRESSED_MATRIX *matrix, double complex value) {
  //First find position of element in A matrix, given position in uncompressed
  //matrix
  int ai, aj, ak;
  if (j > (i+matrix->ku) || j < (i-matrix->kl)) {
    printf("Error writing to A.  Element i=%i, j=%i is not in band.\n",
           i, j);
    return;
  }
  aj = j;
  ai = matrix->kl + matrix->ku + i - j;

  //Given this information, find position of element in *column-major* 
  //matrix structure A.

  ak = matrix->lda*aj+ai;
  matrix->A[ak] = value;
}

void pAelem(int i, int j, COMPRESSED_MATRIX *matrix) {
  //First find position of element in A matrix, given position in uncompressed
  //matrix
  int ai, aj, ak;
  if (j > (i+matrix->ku) || j < (i-matrix->kl)) {
    printf("Error writing to A.  Element i=%i, j=%i is not in band.\n",
           i, j);
    return;
  }
  aj = j;
  ai = matrix->kl + matrix->ku + i - j;

  //Given this information, find position of element in *column-major* 
  //matrix structure A.

  ak = matrix->lda*aj+ai;
  printf("%g + i*%g\n", creal(matrix->A[ak]), cimag(matrix->A[ak]));
}

void wBelem(int i, int j, COMPRESSED_MATRIX *matrix, double complex value) {
  //First find position of element in B matrix, given position in uncompressed
  //matrix
  int bi, bj, bk, bbi, bbk;
  if (j > (i+matrix->ku) || j < (i-matrix->kl)) {
    printf("Error writing to A.  Element i=%i, j=%i is not in band.\n",
           i, j);
    return;
  }
  bj = j;
  bi = matrix->kl + matrix->ku + i - j;
  bbi = matrix->ku + i - j;

  //Given this information, find position of element in *column-major* 
  //matrix structure B.

  bk = matrix->lda*bj+bi;
  bbk = matrix->ldab*bj + bbi;
  matrix->B[bk] = value;
  matrix->Bb[bbk] = value;
}

