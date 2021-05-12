#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include "GSdecomp.h"
#include "backsub.h"


void GS_decomp(gsl_matrix* matrixToQR, gsl_matrix* inputTriangularMatrix){
  int numOfRows = (int) matrixToQR->size1;
  int numOfCols = (int) matrixToQR->size2;

  assert( numOfRows >= numOfCols );

  for (int colId = 0; colId < numOfCols; colId++) {
    gsl_vector* col = gsl_vector_alloc(numOfRows);
    *col = (gsl_matrix_column(matrixToQR, colId)).vector;                        // define cols of matrix as vectors
    double normOfCol = gsl_blas_dnrm2( col );                                   // calc. euclidean R2norm

    gsl_matrix_set(inputTriangularMatrix, colId, colId, normOfCol);              // insert result in matrix
    gsl_vector* orthgnMatCol = gsl_vector_alloc(numOfRows);

    gsl_vector_memcpy(orthgnMatCol, col);
    gsl_vector_scale( orthgnMatCol, 1.0/normOfCol);                                   // scale vector
    gsl_matrix_set_col(matrixToQR, colId, orthgnMatCol);                              // set the new orthogonal col

    for (int nextColId = colId + 1; nextColId < numOfCols; nextColId++) {
      gsl_vector* nextCol = gsl_vector_alloc(numOfRows);
      *nextCol = (gsl_matrix_column(matrixToQR, nextColId)).vector;    // define next col of matrix as vector

      double triangularMatrixElement;
      gsl_blas_ddot(orthgnMatCol, nextCol, &triangularMatrixElement);
      gsl_matrix_set(inputTriangularMatrix, colId, nextColId, triangularMatrixElement); // insert result in matrix

      gsl_vector* ortColScaled = gsl_vector_alloc(numOfRows);
      gsl_vector_memcpy(ortColScaled, orthgnMatCol);                                //  Copy q_i, the i'th (colId) orthogonal column into this new vector
      gsl_vector_scale(ortColScaled, triangularMatrixElement);                             //  Compute q_i*R_ij
      gsl_vector_sub(nextCol, ortColScaled);                                        //  Compute orthogonal complement by subtracting the component along q_i, q_i*R_ij, from a_j,                                                  //  that is, do exactly a_j = a_j - q_i*R_ij

      gsl_matrix_set_col(matrixToQR, nextColId, nextCol);                       // set the new orthogonal col

    }
  }
}


void GS_solve(gsl_matrix* ortMatrix, gsl_matrix* triangularMatrix, gsl_vector* rhsVec, gsl_vector* var){
  gsl_blas_dgemv(CblasTrans, 1.0, ortMatrix, rhsVec, 0.0, var);
  backsub(triangularMatrix, var);
}


void GS_inverse(gsl_matrix* ortMatrix, gsl_matrix* triangularMatrix, gsl_matrix* inverseMatrix){
  int dimension = (int) triangularMatrix->size1;
  gsl_vector* unitVector = gsl_vector_alloc(dimension);
  gsl_matrix* tmpMatrix = gsl_matrix_alloc(dimension, dimension);

  for (int unitVectorId = 0; unitVectorId < dimension; unitVectorId++) {
    gsl_vector_set_basis(unitVector, unitVectorId);
    backsub(triangularMatrix, unitVector);

    gsl_matrix_set_col(tmpMatrix, unitVectorId, unitVector);
  }
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmpMatrix, ortMatrix, 0, inverseMatrix);
  
  gsl_vector_free(unitVector);
  gsl_matrix_free(tmpMatrix);
}
