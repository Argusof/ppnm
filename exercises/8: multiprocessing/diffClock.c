#include <time.h>
#include <math.h>
#include "diffClock.h"

double diffClock(clock_t startTime, clock_t endTime){
  /* Calculates time dirrerence between input arguments
      - startTime : Timestamp of beginning of time interval
      - endTime   : Timestamp of ending of time interval
  */

  double diffTicks  =  fabs( startTime - endTime );
  double diffms     =  (diffTicks * 10) / CLOCKS_PER_SEC;

  return diffms;
}
