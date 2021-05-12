#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>

#include "GSdecomp.h"
#include "backsub.h"
#include "leastSquares.h"

double funcs(int order, double x){
	switch(order){
		case  0: return 1  ; break;
		case  1: return x  ; break;
		default: return NAN;
		}
	}

	void calcDevs(int numOfPts, double* dataY, double* devY) {
	  for (int dev = 0; dev < numOfPts; dev++) {
	    devY[dev] = dataY[dev]/20;
	  }
	}

void logTransformation(int numOfPts, double* dataY, double* dataYtrans, double* devY, double* devYtrans) {
	for (int id = 0; id < numOfPts; id++) {
		dataYtrans[id] = log(dataY[id]);
		devYtrans[id] = devY[id]/dataY[id];
	}
}

void expTransformation(gsl_vector* coeffVector, gsl_matrix* covMatrix) {
	gsl_vector_set(coeffVector, 0, exp( gsl_vector_get(coeffVector,0) ));
	gsl_matrix_set(covMatrix, 0, 0, exp( gsl_matrix_get(covMatrix,0,0) ));

}

void fillArrayFromInputData( double* dataX,  double* dataY, char* inputFileName);

void writeCoeffs(char* outputFileName, double scale, double lambda, double dscale, double dlambda){
	FILE* outputFileStream = fopen(outputFileName, "w");
	fprintf(outputFileStream, "a = %g\nlambda = %g\n", scale, lambda);
	fprintf(outputFileStream, "dap = %g\ndlambdap = %g\n", scale+dscale, lambda+dlambda);
	fprintf(outputFileStream, "dam = %g\ndlambdam = %g\n", scale-dscale, lambda-dlambda);

	fclose(outputFileStream);
}

int main(int argc, char* argv[]) {
  int numOfPts = 9;
  int numOfFuncs = 2;

  double* dataX = malloc(numOfPts*sizeof(double));
  double* dataY = malloc(numOfPts*sizeof(double));
  double* devY = malloc(numOfPts*sizeof(double));
	double* dataYtrans = malloc(numOfPts*sizeof(double));
  double* devYtrans = malloc(numOfPts*sizeof(double));

// READ DATA FROM FILE
  if (argc < 2) {
    fprintf(stderr, "Error: no input arguments\n" );
  }
  else {
    char* inputFileName = argv[1];
    fillArrayFromInputData( dataX, dataY, inputFileName);
  }
	calcDevs(numOfPts, dataY, devY);
	logTransformation(numOfPts, dataY, dataYtrans, devY, devYtrans);

	gsl_matrix* dataMatrix   = gsl_matrix_alloc(numOfPts, numOfFuncs);
  gsl_vector* dataVector   = gsl_vector_alloc(numOfPts);
  gsl_vector* coeffVector  = gsl_vector_alloc(numOfFuncs);
	gsl_matrix* covMatrix    = gsl_matrix_alloc(numOfFuncs, numOfFuncs);


	covMatrix = leastSquares(numOfPts, numOfFuncs, &funcs, dataMatrix, dataVector, coeffVector, dataX, dataYtrans, devYtrans);
	expTransformation(coeffVector, covMatrix);

	double scale = gsl_vector_get(coeffVector,0);
	double lambda = gsl_vector_get(coeffVector,1); //- lambda

	char* outputFileName = argv[2];
	writeCoeffs(outputFileName, scale, lambda, sqrt(gsl_matrix_get(covMatrix,0,0)), sqrt(gsl_matrix_get(covMatrix,1,1)) );

	printf("\n-- Found coefficients from fit: ------------------------------------ \n");
  printf("C_0 (a)      = %lg +/- %lg \n", scale, sqrt(gsl_matrix_get(covMatrix,0,0)));
  printf("C_1 (Lambda) = %lg +/- %lg\n\n", lambda, sqrt(gsl_matrix_get(covMatrix,1,1)));
  printf("The half-life of ThX is                        \tt_1/2 = %lg +/- %lg days\n", log(2)/(-lambda), log(2)/(-gsl_matrix_get(covMatrix,1,1)));
  printf("This should be equal to the half-life of 224Ra \tt_1/2 = %lg days\n", 3.66                );
	print_matrix(numOfFuncs, covMatrix, "The covariance matrix is");

  return 0;
}
