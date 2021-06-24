#ifndef HAVE_MCINTEGRATION_H
#define HAVE_MCINTEGRATION_H


void randomx(int dim, double* lowerBound, double* upperBound, double* pts);
void plainMC(int dim, double* lowerBound, double* upperBound, double f(double* pts), int numOfPts, double* result, double* error );
double corput(int id, int base);
void halton(int id, int dim, double* pts);
void halton2(int id, int dim, double* pts);
void HaltonCorputRandomx(int id, int dim, double* lowerBound, double* upperBound, double* pts);
void HaltonCorputRandomx2(int id, int dim, double* lowerBound, double* upperBound, double* pts);
void HaltonCorputMC(int dim, double* lowerBound, double* upperBound, double f(double* pts), int numOfPts, double* result, double* error );
double stratMC(int dim, double f(double* pts), double* lowerBound, double* upperBound, double absAcc, double relAcc, int numOfRecalls, double meanRecall);

#endif
