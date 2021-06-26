#ifndef HAVE_JACOBI_H
#define HAVE_JACOBI_H

#include <gsl/gsl_matrix.h>


void jacobiMultiply_right(gsl_matrix* matrix, int firstId, int secondId, double angle);
void jacobiMultiply_left(gsl_matrix* matrix, int firstId, int secondId, double angle);
void GivensRotation(gsl_matrix* matrix, int firstId, int secondId, double angle);
void jacobiSVD(gsl_matrix* matrix, gsl_matrix* eigVecMat, gsl_matrix* eigValMat);

#endif
