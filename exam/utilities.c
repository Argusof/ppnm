#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "utilities.h"

double randomNumber( unsigned int *seed ){
  double maxRand      =   (double)RAND_MAX;           // Maximum random number, cast to double
  double randNum      =   (double)rand_r( seed );     // Generate pseudo-random number from seed, cast to double
  return randNum/maxRand;
}

void vector_print(char* string, gsl_vector* vector){
	printf("%s\n", string);
	for(int iter = 0; iter < vector -> size; iter++){
		printf("%10g ", gsl_vector_get(vector, iter));
	}
  printf("\n");
}

void print_matrix(int numOfRows, gsl_matrix* matrixToPrint, char* string ){
  printf("\n%s\n", string);
  for (int rowId = 0; rowId < numOfRows; rowId++){
    gsl_vector_view matrixToPrint_row = gsl_matrix_row (matrixToPrint, rowId);
    gsl_vector* vector = &matrixToPrint_row.vector;
  	for(int iter = 0; iter < vector -> size; iter++){
      if ( gsl_vector_get(vector, iter) > 1e-10 ){
  		    printf("%10g\t", gsl_vector_get(vector, iter));
      }
      else { printf("%10g\t", 0.0); }
  	}
    printf("\n");
  }
}


void set_random_square_matrix(gsl_matrix* inputMatrix, unsigned int *seed, int dim){
  for(int rowId = 0; rowId < dim; rowId++){
    for(int colId = 0; colId < dim; colId++){
      double matrixElement = randomNumber(seed);
  		gsl_matrix_set(inputMatrix, rowId, colId, matrixElement);
		}
  }
}
