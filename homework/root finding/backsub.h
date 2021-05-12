#ifndef HAVE_BACKSUB_H
#define HAVE_BACKSUB_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void backsub(gsl_matrix* upTriangMat, gsl_vector* rhsVec);

#endif
