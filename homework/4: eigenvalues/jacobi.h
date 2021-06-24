#ifndef HAVE_JACOBI_H
#define HAVE_JACOBI_H

#include <gsl/gsl_matrix.h>


void jacobiMultiply_right(gsl_matrix* matrix, int firstId, int secondId, double angle);
void jacobiMultiply_left(gsl_matrix* matrix, int firstId, int secondId, double angle);
void jacobiDiag(gsl_matrix* matrix, gsl_matrix* eigVecMat);
void jacobiMultiply_transleft_right (gsl_matrix* matrix, int firstId, int secondId, double angle);
void jacobiDiag_opt (gsl_matrix* matrix, gsl_matrix* eigVecMat);

#endif
