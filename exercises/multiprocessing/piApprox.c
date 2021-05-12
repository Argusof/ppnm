#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "rndNum.h"
#include "piApprox.h"

void piApprox( int numOfPts ){
  unsigned int seed         =   time(NULL) ;
  int ptsInCircle           =   0          ;
  double unitCircleRadius   =   1.0        ;

  double xpos   =   0;
  double ypos   =   0;

  for ( int iteration = 0; iteration < numOfPts; iteration++ ){
    xpos = randomNumber(&seed);
    ypos = randomNumber(&seed);

    if ( sqrt( pow(xpos, 2) + pow(ypos, 2) ) <= unitCircleRadius ){
      ptsInCircle++;
    }
  }

  double myPiApprox   =   4.0*((double)ptsInCircle)/((double)numOfPts);

  printf("My approximation of pi = %g\n"  , myPiApprox                    );
  printf("Pi is actually = %g\n"          , M_PI                          );
  printf("Deviation is %g percent\n"      , fabs(1.0-myPiApprox/M_PI)*100 );
}
