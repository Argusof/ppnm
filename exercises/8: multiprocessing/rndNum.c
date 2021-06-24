#include <stdlib.h>
#include "rndNum.h"

double randomNumber( unsigned int *seed ){
  /*_Thread_local*/
  double maxRand      =   (double)RAND_MAX;           // Maximum random number, cast to double
  double randNum      =   (double)rand_r( seed );     // Generate pseudo-random number from seed, cast to double

  return 2 * randNum/maxRand - 1;                     // Recast number between -1 and 1
}
