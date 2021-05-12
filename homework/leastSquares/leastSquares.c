#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "leastSquares.h"
#include "GSdecomp.h"

void print_matrix(int numOfRows, gsl_matrix* matrixToPrint, char* string ){
  printf("\n%s\n", string);
  for (int rowId = 0; rowId < numOfRows; rowId++){
    gsl_vector_view matrixToPrint_row = gsl_matrix_row (matrixToPrint, rowId);
    gsl_vector* vector = &matrixToPrint_row.vector;
  	for(int iter = 0; iter < vector -> size; iter++){
      if ( gsl_vector_get(vector, iter) > 1e-10 ){
  		    printf("%10g\t", gsl_vector_get(vector, iter));
      }
      else { printf("%10g\t", 0.0); }
  	}
    printf("\n");
  }
}

void vector_print(char* string, gsl_vector* vector){
	printf("%s\n", string);
	for(int iter = 0; iter < vector -> size; iter++){
		printf("%10g ", gsl_vector_get(vector, iter));
	}
  printf("\n");
}

void setData(int numOfPts, int numOfFuncs, double (*fitFuncs)(int, double), gsl_matrix* dataMatrix, gsl_vector* dataVector, double* dataX, double* dataY, double* devY){
  for (int row = 0; row < numOfPts; row++) {
    for (int col = 0; col < numOfFuncs; col++) {
      gsl_matrix_set(dataMatrix, row, col, (fitFuncs(col, dataX[row]))/devY[row]);
    }
    gsl_vector_set(dataVector, row, dataY[row]/devY[row]);
  }
}


gsl_matrix* leastSquares(int numOfPts, int numOfFuncs, double (*fitFuncs)(int, double), gsl_matrix* dataMatrix, gsl_vector* dataVector, gsl_vector* coeffVector, double* dataX, double* dataY, double* devY) {

  setData(numOfPts, numOfFuncs, fitFuncs, dataMatrix, dataVector, dataX, dataY, devY);

  gsl_matrix* ortgMatrix          = gsl_matrix_alloc(numOfPts, numOfFuncs);
  gsl_matrix* triangMatrix        = gsl_matrix_alloc(numOfFuncs, numOfFuncs);
  gsl_matrix* triangMatrixInverse = gsl_matrix_alloc(numOfFuncs, numOfFuncs);
  gsl_matrix* covMatrix           = gsl_matrix_alloc(numOfFuncs, numOfFuncs);

  //gsl_vector* tmpRes              = gsl_vector_alloc(numOfFuncs);

  gsl_matrix_memcpy(ortgMatrix, dataMatrix);
  GS_decomp(ortgMatrix, triangMatrix);
  GS_solve(ortgMatrix, triangMatrix, dataVector, coeffVector);

  GS_inverse_triang(triangMatrix, triangMatrixInverse);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, triangMatrixInverse, triangMatrixInverse, 0.0, covMatrix);

  return covMatrix;
}
