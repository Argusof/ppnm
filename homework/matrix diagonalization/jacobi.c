#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "jacobi.h"

void jacobiMultiply_right(gsl_matrix* matrix, int firstId, int secondId, double angle){
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
	double c = cos(angle);
  double s = sin(angle);

	for(int col = 0; col < matrix->size2; col++){
		double new_apj =  c*gsl_matrix_get(matrix,firstId,col) + s*gsl_matrix_get(matrix,secondId,col);
		double new_aqj = -s*gsl_matrix_get(matrix,firstId,col) + c*gsl_matrix_get(matrix,secondId,col);
		gsl_matrix_set(matrix,firstId,col,new_apj);
		gsl_matrix_set(matrix,secondId,col,new_aqj);
		}
}


void jacobiDiag(gsl_matrix* matrix, gsl_matrix* eigVecMat){
gsl_matrix_set_identity(eigVecMat);
int dims = matrix->size1;
int changed;

  do{
  	changed = 0;
  	for(int firstId = 0; firstId < dims-1; firstId++){
    	for(int secondId = firstId+1; secondId < dims; secondId++){
    		double apq = gsl_matrix_get(matrix,firstId,secondId);
    		double app = gsl_matrix_get(matrix,firstId,firstId);
    		double aqq = gsl_matrix_get(matrix,secondId,secondId);

        double angle = 0.5 * atan2(2*apq, aqq-app);
				double c = cos(angle);
				double s = sin(angle);

        double new_app = c*c*app - 2*s*c*apq + s*s*aqq;
    		double new_aqq = s*s*app + 2*s*c*apq + c*c*aqq;

        if(new_app != app || new_aqq != aqq) // do rotation
    			{
      			changed = 1;
      			jacobiMultiply_right(matrix, firstId, secondId, angle);
      			jacobiMultiply_left(matrix, firstId, secondId, -angle); // A←J^T*A*J
      			jacobiMultiply_right(eigVecMat, firstId, secondId, angle); // V←V*J
    			}
    	}
    }
  } while(changed != 0);

}
