#include <pthread.h>
#include <omp.h>
#include <stdio.h>
#include "piApproxMultithreaded.h"
#include "piApproxOpenMP.h"
#include "diffClock.h"

int main( void ){

  clock_t begin = clock(); // We define variables to hold the times used for timing the computations.
  clock_t end   = clock(); // These are defined in <time.h> and used in diffClock(). It is my own implementation.

  double numOfPts = 1e9;

  begin = clock(); // Begin timing
  piApproxMultithreaded( numOfPts );
  end = clock(); // End timing
  printf("Done! Elapsed time (pThreads) = %g ms \n", (double)(diffClock(end, begin)) );

  begin = clock(); // Begin timing
  piApproxOpenMP( numOfPts );
  end = clock(); // End timing
  printf("Done! Elapsed time (OpenMP) = %g ms \n", (double)(diffClock(end, begin)) );

  /*
  Here we see that using OpenMP is much faster. This is because OpenMP, with the way we implement parallelization
  allows an arbitrary function signature. Hence, we are not bound to use structs and dynamic memory allocations,
  and can thus write much faster (and safer) code.
  */

  return 0;
}
