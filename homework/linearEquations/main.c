#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "GSdecomp.h"
#include "utilities.h"
#include "diffClock.h"
#include "backsub.h"

#define RND (double)rand()/RAND_MAX

void test_runtime(int numOfReps, int startRep, char* my_outputFilename, char* gsl_outputFilename, unsigned int* seed);


int main (int argc, char* argv[]){
  unsigned int seed   =   time(NULL) ;

  int numOfRows = 5;
	int numOfCols = 4;

  gsl_matrix* ortgnlMatTall    =  gsl_matrix_alloc(numOfRows, numOfCols);
  gsl_matrix* ortgnlMatSquare  =  gsl_matrix_alloc(numOfRows, numOfRows);
  gsl_matrix* testMatTall      =  gsl_matrix_alloc(numOfRows, numOfCols);
  gsl_matrix* testMatSquare    =  gsl_matrix_alloc(numOfRows, numOfRows);
  gsl_matrix* triangMatTall    =  gsl_matrix_alloc(numOfCols, numOfCols);
  gsl_matrix* triangMatSquare  =  gsl_matrix_alloc(numOfRows, numOfRows);
  gsl_matrix* invMatSquare     =  gsl_matrix_alloc(numOfRows, numOfRows);
  gsl_vector* RHSvecTall       =  gsl_vector_alloc(numOfRows           );
  gsl_vector* RHSvecSquare     =  gsl_vector_alloc(numOfRows           );
  gsl_vector* varVec           =  gsl_vector_alloc(numOfRows           );
  gsl_vector* varVecTall       =  gsl_vector_alloc(numOfCols           );
  gsl_matrix* checkOrtgnl      =  gsl_matrix_alloc(numOfCols, numOfCols);
  gsl_matrix* checkRes         =  gsl_matrix_alloc(numOfRows, numOfCols);
  gsl_matrix* checkInv1        =  gsl_matrix_alloc(numOfRows, numOfRows);
  gsl_matrix* checkInv2        =  gsl_matrix_alloc(numOfRows, numOfRows);
  gsl_vector* checkRHS         =  gsl_vector_alloc(numOfRows           );

	set_data_tall(testMatTall, RHSvecTall, &seed, &seed);
  set_data_square(testMatSquare, RHSvecSquare, &seed, &seed);
  gsl_matrix_memcpy(ortgnlMatTall, testMatTall);

  GS_decomp(ortgnlMatTall, triangMatTall);

  print_matrix(numOfCols, triangMatTall, "Triangular matrix:");

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, ortgnlMatTall, ortgnlMatTall, 0, checkOrtgnl);
  print_matrix(numOfCols, checkOrtgnl, "trans(Q) * Q (should be equal to 1):");

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, ortgnlMatTall, triangMatTall, 0, checkRes);
  print_matrix(numOfRows, checkRes, "QR :");
  print_matrix(numOfRows, testMatTall, "Tall test matrix:");

  print_matrix(numOfRows, testMatSquare, "Square test matrix for solver:");
  gsl_matrix_memcpy(ortgnlMatSquare, testMatSquare);
  GS_decomp(ortgnlMatSquare, triangMatSquare);
  GS_solve(ortgnlMatSquare, triangMatSquare, RHSvecSquare, varVec);

  gsl_blas_dgemv(CblasNoTrans, 1, testMatSquare, varVec, 0, checkRHS);
  vector_print( "\nFrom gramSchmidt_solve we get Ax = ", checkRHS);
  vector_print( "\nwhich should be equal to right-hand side  b = ", RHSvecSquare);

  GS_inverse(ortgnlMatSquare, triangMatSquare, invMatSquare);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, invMatSquare, testMatSquare, 0, checkInv1);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, testMatSquare, invMatSquare, 0, checkInv2);
  print_matrix(numOfRows, invMatSquare, "Inverse matrix:");
  print_matrix(numOfRows, checkInv1, "inv(A)*A (should be 1):");
  print_matrix(numOfRows, checkInv2, "A*inv(A) (should be 1):");


  clock_t myStartTime  = clock();
  clock_t myEndTime    = clock();
  clock_t gslStartTime = clock();
  clock_t gslEndTime   = clock();
  double myRunTime     = 0.0;
  double gslRunTime    = 0.0;
  int maxSize          = 100;
  FILE* outputFile     = fopen( argv[1], "w");

  for (int matrixSize = 2; matrixSize < maxSize; matrixSize++) {

    gsl_matrix* myOrtgMatrix  = gsl_matrix_alloc(matrixSize, matrixSize);
    gsl_matrix* myTriMatrix   = gsl_matrix_alloc(matrixSize, matrixSize);
    gsl_vector* tmpVector     = gsl_vector_alloc(matrixSize);
    gsl_matrix* gslOrtgMatrix = gsl_matrix_alloc(matrixSize, matrixSize);
    gsl_matrix* gslTriMatrix  = gsl_matrix_alloc(matrixSize, matrixSize);

    set_data_square(myOrtgMatrix, tmpVector, &seed, &seed);
    set_data_square(myTriMatrix, tmpVector, &seed, &seed);
    set_data_square(gslOrtgMatrix, tmpVector, &seed, &seed);

    gsl_matrix_free(myOrtgMatrix);
    gsl_matrix_free(myTriMatrix);
    gsl_vector_free(tmpVector);
    gsl_matrix_free(gslOrtgMatrix);
    gsl_matrix_free(gslTriMatrix);
  }
  fclose(outputFile);


  gsl_matrix_free (ortgnlMatTall);
  gsl_matrix_free (ortgnlMatSquare);
	gsl_matrix_free (testMatTall);
	gsl_matrix_free (testMatSquare);
	gsl_matrix_free (triangMatTall);
	gsl_matrix_free (triangMatSquare);
  gsl_matrix_free (invMatSquare);
	gsl_vector_free (RHSvecTall);
	gsl_vector_free (RHSvecSquare);
	gsl_vector_free (varVec);
  gsl_matrix_free (checkOrtgnl);
  gsl_matrix_free (checkRes);
	gsl_vector_free (checkRHS);
  gsl_matrix_free (checkInv1);
  gsl_matrix_free (checkInv2);


  int numOfReps = 300;
  int startRep = 200;
  test_runtime(numOfReps, startRep, argv[1], argv[2], &seed);

  return 0;
}
