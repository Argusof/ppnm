#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#define rand (double)rand()/RAND_MAX

void vector_print( char* string, gsl_vector* vector){
	printf("%s\n", string);
	for(int iter = 0; iter < vector->size; iter++){
		printf("%10g ", gsl_vector_get(vector, iter));
	}
}

int main( ) {
  // Dimensions of matrix
	int numOfRows = 3;
	int numOfCols = 3;

  // Allocate all matrices and vectors
	gsl_matrix* matrix_A_copy  =  gsl_matrix_alloc(numOfRows, numOfCols);
	gsl_vector* vector_x       =  gsl_vector_alloc(numOfRows);
	gsl_vector* vector_y       =  gsl_vector_calloc(numOfRows);

	double matrix_A_data[]  =  { 6.13, -2.90,  5.86,
                      			   8.08, -6.31, -3.89,
                      			  -4.36,  1.00,  0.19 };

  double vector_b_data[]  =  { 6.23,  5.37,  2.29 };

  gsl_matrix_view matrix_A  =  gsl_matrix_view_array (matrix_A_data, numOfRows, numOfCols);
  gsl_vector_view vector_b  =  gsl_vector_view_array (vector_b_data, numOfRows);

	gsl_matrix_memcpy(matrix_A_copy, &matrix_A.matrix);
  gsl_linalg_HH_solve(matrix_A_copy, &vector_b.vector, vector_x);
	gsl_blas_dgemv(CblasNoTrans, 1, &matrix_A.matrix, vector_x, 0, vector_y);

	vector_print("Righthand side b = ", &vector_b.vector);
	printf("\n");
	vector_print("Check: A*x should be equal to b...:", vector_y);
	printf("\n");

  gsl_matrix_free(matrix_A_copy);
	gsl_vector_free(vector_x);
  gsl_vector_free(vector_y);

  int order = 4;
  gsl_matrix* hillbertMatrix = gsl_matrix_alloc(order, order);

  for (size_t i = 0; i < order; i++) {
    for (size_t j = 0; j < order; j++) {
      gsl_matrix_set( hillbertMatrix, i, j, 1.0/(i + j + 1) );
    }
  }

  gsl_vector *eigenvals = gsl_vector_alloc (order);
  gsl_matrix *eigenvecs = gsl_matrix_alloc (order, order);

  gsl_eigen_symmv_workspace * workspace = gsl_eigen_symmv_alloc (order);
  gsl_eigen_symmv (hillbertMatrix, eigenvals, eigenvecs, workspace);
  gsl_eigen_symmv_free (workspace);
  gsl_eigen_symmv_sort (eigenvals, eigenvecs, GSL_EIGEN_SORT_ABS_ASC);

  for (int i = 0; i < order; i++)
      {
        double eigenvals_i = gsl_vector_get (eigenvals, i);
        gsl_vector_view eigenvecs_i = gsl_matrix_column (eigenvecs, i);

        printf ("eigenvalue = %g\n", eigenvals_i);
        printf ("eigenvector = \n");
        gsl_vector_fprintf (stdout, &eigenvecs_i.vector, "%g");
      }

	return 0;
}
