#ifndef HAVE_GSDECOMP_H
#define HAVE_GSDECOMP_H

#include <gsl/gsl_matrix.h>

void GS_decomp(gsl_matrix* matrixToQR, gsl_matrix* inputTriangularMatrix);
void GS_solve(gsl_matrix* ortMatrix, gsl_matrix* triangularMatrix, gsl_vector* rhsVec, gsl_vector* var);
void GS_inverse(gsl_matrix* ortMatrix, gsl_matrix* triangularMatrix, gsl_matrix* inverseMatrix);
void GS_inverse_triang(gsl_matrix* triangularMatrix, gsl_matrix* inverseMatrix);



#endif
