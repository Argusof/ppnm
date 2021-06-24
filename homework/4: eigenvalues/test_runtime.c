#include <time.h>
#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "diffClock.h"
#include "jacobi.h"
#include "utilities.h"


void test_runtime(int numOfReps, int startRep, const char* my_outputFilename, const char* gsl_outputFilename, unsigned int* seed){
  double scale      =  0   ;
  int numOfDims;

  FILE* myOutputFileStream      =  fopen(my_outputFilename,  "w");
  FILE* myOutputFileStream_gsl  =  fopen(gsl_outputFilename, "w");


  for (int rep = startRep; rep < numOfReps; rep++){
    numOfDims = rep;

    gsl_matrix* my_symm    =  gsl_matrix_alloc(numOfDims, numOfDims);
    gsl_matrix* my_eigVec  =  gsl_matrix_alloc(numOfDims, numOfDims);

    set_data_symm(my_symm, seed);

    clock_t my_begin  = clock(); // We define variables to hold the times used for timing the computations.
    clock_t my_end    = clock(); // These are defined in <time.h> and used in diffClock(). It is my own implementation.

    my_begin = clock(); // Begin timing

    // time something....
    jacobiDiag_opt(my_symm, my_eigVec);

    my_end = clock(); // end timing

    if (rep == startRep){
      scale = (double)(diffClock(my_end, my_begin));
    }

    fprintf(myOutputFileStream, "%d\t%g\t%g\n", numOfDims, (double)(diffClock(my_end, my_begin)), pow(((double)numOfDims)/startRep, 3)*scale);

    gsl_matrix_free(my_symm);
    gsl_matrix_free(my_eigVec);
  }

  for (int rep = startRep; rep < numOfReps; rep++){
    numOfDims = rep;

    gsl_matrix* gsl_symm    =  gsl_matrix_alloc(numOfDims, numOfDims);
    gsl_matrix* gsl_eigVec  =  gsl_matrix_alloc(numOfDims, numOfDims);
    gsl_vector* gsl_diag    =  gsl_vector_alloc(numOfDims );


    set_data_symm(gsl_symm, seed);

    clock_t gsl_begin = clock();
    clock_t gsl_end   = clock();

    gsl_begin = clock();
    gsl_linalg_SV_decomp_jacobi(gsl_symm, gsl_eigVec, gsl_diag);
    gsl_end = clock();

    fprintf(myOutputFileStream_gsl, "%d\t%g\n", numOfDims, (double)(diffClock(gsl_end, gsl_begin)));

    gsl_matrix_free(gsl_symm);
    gsl_matrix_free(gsl_eigVec);
    gsl_vector_free(gsl_diag);
  }
  fclose(myOutputFileStream);
  fclose(myOutputFileStream_gsl);
}
