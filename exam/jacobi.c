#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "jacobi.h"

void jacobiMultiply_right(gsl_matrix* matrix, int firstId, int secondId, double angle){
	// J applied from the right

	double c = cos(angle);
  double s = sin(angle);

  for(int row = 0; row < matrix->size1; row++){
		double new_aip = c*gsl_matrix_get(matrix,row,firstId) - s*gsl_matrix_get(matrix,row,secondId);
		double new_aiq = s*gsl_matrix_get(matrix,row,firstId) + c*gsl_matrix_get(matrix,row,secondId);

		gsl_matrix_set(matrix,row,firstId,new_aip);
		gsl_matrix_set(matrix,row,secondId,new_aiq);
		}
}


void jacobiMultiply_left(gsl_matrix* matrix, int firstId, int secondId, double angle){
	// J applied from the right
	double c = cos(angle);
  double s = sin(angle);

	for(int col = 0; col < matrix->size2; col++){
		double new_apj =  c*gsl_matrix_get(matrix,firstId,col) + s*gsl_matrix_get(matrix,secondId,col);
		double new_aqj = -s*gsl_matrix_get(matrix,firstId,col) + c*gsl_matrix_get(matrix,secondId,col);
		gsl_matrix_set(matrix,firstId,col,new_apj);
		gsl_matrix_set(matrix,secondId,col,new_aqj);
	}
}


void jacobiSVD(gsl_matrix* matrix, gsl_matrix* eigVecMat, gsl_matrix* eigValMat){
	int maxIt = 1e5;
	int ItId  = 0;

	gsl_matrix_set_identity(eigVecMat);
	gsl_matrix_set_identity(eigValMat);

	int dims = matrix->size1;
	int changed;

  do{
  	changed = 0;

  	for(int firstId = 0; firstId < dims-1; firstId++){
    	for(int secondId = firstId+1; secondId < dims; secondId++){
    		double apq = gsl_matrix_get(matrix,firstId,secondId);
				double aqp = gsl_matrix_get(matrix,secondId,firstId);
    		double app = gsl_matrix_get(matrix,firstId,firstId);
    		double aqq = gsl_matrix_get(matrix,secondId,secondId);

        double GivensAngle = atan2(apq - aqp, aqq + app);
				double JacobiAngle = 0.5 * atan2(2*apq, aqq-app);

				double c = cos(JacobiAngle);
				double s = sin(JacobiAngle);

        double new_app = c*c*app - 2*s*c*apq + s*s*aqq;
    		double new_aqq = s*s*app + 2*s*c*apq + c*c*aqq;

        if(new_app != app || new_aqq != aqq) // do rotation
    			{
      			changed = 1;
						ItId++;

						// Givens rotation to get symmetric matrix A'. This is the same as applying the
						// Jacobian matrix from the left, however, with the angle defined differently
						jacobiMultiply_left(matrix, firstId, secondId, -GivensAngle);
						// A'←J^T*A'*J
						jacobiMultiply_right(matrix, firstId, secondId, JacobiAngle);
      			jacobiMultiply_left(matrix, firstId, secondId, -JacobiAngle);

						// U-> UGJ
						jacobiMultiply_right(eigValMat, firstId, secondId, GivensAngle);
						jacobiMultiply_right(eigValMat, firstId, secondId, JacobiAngle);

						// V←V*J
      			jacobiMultiply_right(eigVecMat, firstId, secondId, JacobiAngle);

						if (ItId == maxIt){
							break;
						}
    			}

    	}
			if (ItId == maxIt){
				break;
			}
    }
		if (ItId == maxIt){
			break;
		}
  } while(changed != 0);

}
