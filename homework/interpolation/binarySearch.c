#include "binarySearch.h"

int binarySearch ( int numOfPts, double* pts, double evalPt ) {

  //  Define indices of interval bounds
  int leftId   =  0     ;
  int rightId  =  numOfPts - 1 ;

  //  Use a binary search to determine which subinterval evalPt is in
  while ( rightId - leftId > 1 ) {      // While there are two intervals to choose from
    int middleId  =  (leftId + rightId) / 2;

    if ( evalPt > pts[middleId] ) {   // If the point is in the right interval
      leftId   =  middleId;
    }
    else {                              // Else it is in the left interval
      rightId  =  middleId;
    }
  }
  int whichInterval = leftId;

  return whichInterval;
}
