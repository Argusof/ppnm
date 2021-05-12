#include <time.h>
#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "diffClock.h"
#include "GSdecomp.h"
#include "utilities.h"


void test_runtime(int numOfReps, int startRep, char* my_outputFilename, char* gsl_outputFilename, unsigned int* seed){
  double scale      =  0   ;
  int numOfDims;

  FILE* myOutputFileStream      =  fopen(my_outputFilename,  "w");
  FILE* myOutputFileStream_gsl  =  fopen(gsl_outputFilename, "w");


  for (int rep = startRep; rep < numOfReps; rep++){
    numOfDims = rep;

    gsl_matrix* my_ortg    =  gsl_matrix_alloc(numOfDims, numOfDims);
    gsl_matrix* my_triang  =  gsl_matrix_alloc(numOfDims, numOfDims);
    gsl_vector* my_vecTmp  =  gsl_vector_alloc(numOfDims           );

    set_data_tall(my_ortg, my_vecTmp, seed, seed);

    clock_t my_begin  = clock(); // We define variables to hold the times used for timing the computations.
    clock_t my_end    = clock(); // These are defined in <time.h> and used in diffClock(). It is my own implementation.

    my_begin = clock(); // Begin timing

    // time something....
    GS_decomp(my_ortg, my_triang);

    my_end = clock(); // end timing

    if (rep == startRep){
      scale = (double)(diffClock(my_end, my_begin));
    }

    fprintf(myOutputFileStream, "%d\t%g\t%g\n", numOfDims, (double)(diffClock(my_end, my_begin)), pow(((double)numOfDims)/startRep, 3)*scale);

    gsl_matrix_free(my_ortg);
    gsl_matrix_free(my_triang);
    gsl_vector_free(my_vecTmp);
  }

  for (int rep = startRep; rep < numOfReps; rep++){
    numOfDims = rep;

    gsl_matrix* gsl_ortg    =  gsl_matrix_alloc(numOfDims, numOfDims);
    gsl_matrix* gsl_triang  =  gsl_matrix_alloc(numOfDims, numOfDims);
    gsl_vector* my_vecTmp   =  gsl_vector_alloc(numOfDims);

    set_data_tall(gsl_ortg, my_vecTmp , seed, seed);

    clock_t gsl_begin = clock();
    clock_t gsl_end   = clock();

    gsl_begin = clock();
    gsl_linalg_QR_decomp(gsl_ortg, my_vecTmp);
    gsl_end = clock();

    fprintf(myOutputFileStream_gsl, "%d\t%g\n", numOfDims, (double)(diffClock(gsl_end, gsl_begin)));

    gsl_matrix_free(gsl_ortg);
    gsl_matrix_free(gsl_triang);
    gsl_vector_free(my_vecTmp);
  }
  fclose(myOutputFileStream);
  fclose(myOutputFileStream_gsl);
}
