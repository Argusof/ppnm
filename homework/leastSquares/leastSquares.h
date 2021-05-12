#ifndef HAVE_LEASTSQUARES_H
#define HAVE_LEASTSQUARES_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void setData(int numOfPts, int numOfFuncs, double (*fitFuncs)(int, double), gsl_matrix* dataMatrix, gsl_vector* dataVector, double* dataX, double* dataY, double* devY);
gsl_matrix* leastSquares(int numOfPts, int numOfFuncs, double (*fitFuncs)(int, double), gsl_matrix* dataMatrix, gsl_vector* dataVector, gsl_vector* coeffVector, double* dataX, double* dataY, double* devY) ;

#endif
