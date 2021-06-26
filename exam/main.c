#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <time.h>
#include "jacobi.h"
#include "utilities.h"


int main(int argc, char const *argv[]) {

  unsigned int seed = time(NULL);                       // Seed to generate random numbers for matrix
  int dim = 3;                                          // Dimension of matrix

  // Alloc matrices
  gsl_matrix* matrix     = gsl_matrix_alloc(dim, dim);  // Matrix (square) to be diagonalized
  gsl_matrix* matrixCopy = gsl_matrix_alloc(dim, dim);
  gsl_matrix* eigVecMat  = gsl_matrix_alloc(dim, dim);  // Matrix which contains eigenvectors
  gsl_matrix* eigValMat  = gsl_matrix_alloc(dim, dim);  // Matrix which contains eigenvalues

  set_random_square_matrix(matrix, &seed, dim);         // Fill matrix with random numbers
  gsl_matrix_memcpy(matrixCopy,matrix);

  print_matrix(dim, matrix, "Matrix A:");                // Print original matrix

  jacobiSVD(matrix, eigVecMat, eigValMat);               // Do two sided Jacobi SVD
  print_matrix(dim, eigVecMat, "Eigenvectors V:");       // Print computed eigenvectors
  print_matrix(dim, eigValMat, "Eigenvalues U:");        // Print computed eigenvalues

  print_matrix(dim, matrix, "Diagonal matrix D:");      // Print original diagonalized matrix
  // Test
  gsl_matrix* tmp = gsl_matrix_alloc(dim, dim);
  gsl_matrix* testIdentityV = gsl_matrix_alloc(dim, dim);
  gsl_matrix* testIdentityU = gsl_matrix_alloc(dim, dim);
  gsl_matrix* testDiagonal = gsl_matrix_alloc(dim, dim);
  gsl_matrix* testMat = gsl_matrix_alloc(dim, dim);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, eigVecMat, eigVecMat, 0.0, testIdentityV);
  print_matrix(dim, testIdentityV, "Identity test, V^T * V = 1:");

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, eigValMat, eigValMat, 0.0, testIdentityU);
  print_matrix(dim, testIdentityU, "Identity test, U^T * U = 1:");

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, matrix, eigVecMat, 0.0, tmp);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigValMat, tmp, 0.0, testMat);
  print_matrix(dim, testMat, "Test of U * D * V^T = A:");

  gsl_matrix_free(tmp);
  gsl_matrix_free(testIdentityV);
  gsl_matrix_free(testIdentityU);
  gsl_matrix_free(testDiagonal);
  gsl_matrix_free(testMat);
  gsl_matrix_free(matrix);
  gsl_matrix_free(matrixCopy);
  gsl_matrix_free(eigValMat);
  gsl_matrix_free(eigVecMat);

  return 0;
}
