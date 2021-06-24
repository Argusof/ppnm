#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "rndNum.h"
#include "piApproxOpenMP.h"

unsigned long int placePtsOpenMP( unsigned int seed, double numOfPts, double unitCircleRadius ){

  double xpos   =   0;
  double ypos   =   0;
  unsigned long int returnCount = 0;

  for ( int iteration = 0; iteration < numOfPts/2; iteration++ ){
    xpos = randomNumber(&seed);
    ypos = randomNumber(&seed);

    if ( sqrt( pow(xpos, 2) + pow(ypos, 2) ) <= unitCircleRadius ){
      returnCount++;
    }
  }

  return returnCount;
}


void piApproxOpenMP( int numOfPts ){
  unsigned int masterSeed   =   time(NULL) ;
  unsigned int firstSeed    =   masterSeed ;
  unsigned int secondSeed   =   masterSeed + 1 ;
  double unitCircleRadius   =   1.0        ;

  unsigned long int returnCount_first   =  0;
  unsigned long int returnCount_second  =  0;
  #pragma omp parallel sections
    {
  	#pragma omp section
      {
      returnCount_first = placePtsOpenMP( firstSeed,  numOfPts,  unitCircleRadius );
      }
  	#pragma omp section
      {
      returnCount_second = placePtsOpenMP( secondSeed, numOfPts, unitCircleRadius );
      }
    }

  double ptsInCircle   =   (returnCount_first + returnCount_second);
  double myPiApprox    =   4.0*ptsInCircle/((double)numOfPts);

  printf("My approximation of pi using OpenMP = %g\n"  , myPiApprox                    );
  printf("Deviation is %g percent\n"      , fabs(1.0-myPiApprox/M_PI)*100 );
}
