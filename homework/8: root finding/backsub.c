#include "backsub.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void backsub(gsl_matrix* upTriangMat, gsl_vector* rhsVec){
  int numOfRows = (rhsVec -> size);

  for ( int rowId = numOfRows - 1; rowId >= 0; rowId-- ){
    double rhsVal = gsl_vector_get(rhsVec, rowId);

    for ( int varId = rowId + 1; varId < numOfRows; varId++ ){
      rhsVal -= gsl_matrix_get(upTriangMat, rowId, varId) * gsl_vector_get(rhsVec, varId);
    }
    gsl_vector_set(rhsVec, rowId, rhsVal/gsl_matrix_get(upTriangMat, rowId, rowId));
  }
}
