#include <assert.h>
#include "binarySearch.h"
#include "linSpline.h"
#include "math.h"

double linSplineInterp( int numOfPts, double* pts, double* funcVals, double evalPt){
  //  ------------------------------------------------------------------------------
  /*  Do linear spline interpolation on a set of known points (x_i, f(x_i))
      using a binary search method. This asumes the dataset is an ordered array.
      Otherwise use sorting before passing to the function.

      ¤ int       numOfPts  : The number of points to interpolate in.
      ¤ double*   pts       : A pointer to an array of doubles,
                              the known points {x_i}.
      ¤ double*   funcVals  : A pointer to an array of doubles,
                              the corresponding function values {f(x_i)}
      ¤ double    evalPt    : The query point, double values at which to evaluate
                              the interpolant.

      Returns: A double P(evalPt), where P() is the interpolant function             */
  //  ------------------------------------------------------------------------------

  //  We should start by ensuring we have more than one point, and that the
  //  query point, evalPt is within the interval defined by the dataset
  assert( (numOfPts > 1) && (evalPt >= pts[0]) && (evalPt <= pts[numOfPts - 1]) );

  int whichInterval = binarySearch(numOfPts, pts, evalPt );

  double funcValDiff = (funcVals[whichInterval + 1] - funcVals[whichInterval]);
  double ptsDiff = (pts[whichInterval + 1] - pts[whichInterval]);

  assert( ptsDiff > 0 );

  double interpVal  =  funcVals[whichInterval] + funcValDiff / ptsDiff * (evalPt - pts[whichInterval]);

  return interpVal;
}

/*double linSplineInterpIntegral( int numOfPts, double* pts, double* funcVals, double evalPt){

  assert( (numOfPts > 1) && (evalPt >= pts[0]) && (evalPt <= pts[numOfPts - 1]) );

  int whichInterval = binarySearch(numOfPts, pts, evalPt );

  double funcValDiff = (funcVals[whichInterval + 1] - funcVals[whichInterval]);
  double ptsDiff = (pts[whichInterval + 1] - pts[whichInterval]);
  double slope = funcValDiff / ptsDiff;

  assert( ptsDiff > 0 );

  double interpValInt = funcVals[whichInterval] * ( evalPt - pts[whichInterval]) + slope * pow( ( evalPt - pts[whichInterval] ), 2) / 2 ;//+ funcVals[whichInterval];
  return interpValInt;
}
*/

double linSplineDefiniteIntegral( int numOfPts, double* pts, double* funcVals, double evalPt){

  int whichInterval = binarySearch( numOfPts, pts, evalPt ); // Find interval where z is (left)

  double integral = 0;
  for ( int intervalId = 0; intervalId <= whichInterval; intervalId ++) {

    double funcValDiff = (funcVals[intervalId + 1] - funcVals[intervalId]);
    double ptsDiff = (pts[intervalId + 1] - pts[intervalId]);
    double slope = funcValDiff / ptsDiff;

    if ( intervalId < whichInterval){
      integral += funcVals[intervalId] * ( pts[intervalId+1] - pts[intervalId] ) + slope * pow( ( pts[intervalId+1] - pts[intervalId] ), 2) / 2 ;
    }

    else { //we are at the interval where z is
      integral += funcVals[intervalId] * ( evalPt - pts[intervalId] ) + slope * pow( ( evalPt - pts[intervalId] ), 2) / 2 ;
    }
  }

  return integral;
}
