#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "binarySearch.h"
#include "cubicSpline.h"


cubicSpline* cubicSplineAlloc( int numOfPts, double* pts, double* funcVals ){
  //  ------------------------------------------------------------------------------
  /*  cubicSpline constructor to initiallize a cubic spline from a cubicSpline
      struct, by filling the various field values from respective function inputs.

      ¤ int       numOfPts  : The number of points to interpolate in.
      ¤ double*   pts       : A pointer to an array of doubles,
                              the known points {x_i}.
      ¤ double*   funcVals  : A pointer to an array of doubles,
                              the corresponding function values {f(x_i)}

      Returns: An inittialized cubicSpline* struct                                    */
  //  ------------------------------------------------------------------------------

  int numOfEqs           =  numOfPts - 1;
  cubicSpline* spline    =  (cubicSpline*)malloc( sizeof(cubicSpline) );
  spline -> pts          =  (double*)malloc( numOfPts*sizeof(double) );
  spline -> funcVals     =  (double*)malloc( numOfPts*sizeof(double) );
  spline -> firstCoeff   =  (double*)malloc( numOfPts*sizeof(double) );
  spline -> secondCoeff  =  (double*)malloc( numOfEqs*sizeof(double) );
  spline -> thirdCoeff   =  (double*)malloc( numOfEqs*sizeof(double) );
  spline -> numOfPts     =  numOfPts;

  for( int it = 0; it < numOfPts; it++ ){
    spline -> pts[it]       =  pts[it];
    spline -> funcVals[it]  =  funcVals[it];
  }

  double ptsDiff[numOfEqs];
  double slope[numOfEqs];
  for(int it = 0; it < numOfEqs; it++){
    ptsDiff[it] = pts[it + 1] - pts[it];
    assert( ptsDiff[it] > 0 );

    slope[it] = (funcVals[it + 1] - funcVals[it]) / ptsDiff[it];
  }

  double D[numOfPts    ];
  double Q[numOfPts - 1];
  double B[numOfPts    ];

  D[0] = 2            ;
  Q[0] = 1            ;
  B[0] = 3 * slope[0] ;
  for( int it = 0; it < numOfEqs - 1; it++ ) {
    D[it + 1] = 2 * ptsDiff[it] / ptsDiff[it + 1] + 2;
    Q[it + 1] =     ptsDiff[it] / ptsDiff[it + 1];
    B[it + 1] = 3 * (slope[it] + slope[it + 1] * ptsDiff[it] / ptsDiff[it + 1]);
  }
  D[numOfPts - 1] = 2                      ;
  B[numOfPts - 1] = 3 * slope[numOfPts - 2];

  for(int it = 1; it < numOfPts; it++){
    D[it] -= Q[it - 1] / D[it - 1];
    B[it] -= B[it - 1] / D[it - 1];
  }
  spline -> firstCoeff[numOfEqs] = B[numOfPts - 1] / D[numOfPts - 1];
  for( int it = numOfEqs - 1; it >= 0; it--){
    spline -> firstCoeff[it] = (B[it] - Q[it]*(spline -> firstCoeff[it + 1])) / D[it];
  }
  for( int it = 0; it < numOfEqs; it++ ){
    spline -> secondCoeff[it]  =  (-2 * (spline -> firstCoeff[it]) - (spline -> firstCoeff[it + 1]) + 3*slope[it]) / ptsDiff[it];
    spline -> thirdCoeff[it]   =  ((spline -> firstCoeff[it]) + (spline -> firstCoeff[it + 1]) - 2*slope[it]) / ptsDiff[it] /ptsDiff[it];
  }

  return spline;
}

double cubicSplineEval( cubicSpline* spline, double evalPt ){
  //  ------------------------------------------------------------------------------
  /*  Do cubic spline interpolation, using an already initiallized cubicSpline*
      struct. The cubicSpline* struct may be initiallized using cubicSpline_init().
      The cubic interpolant is computed at the evaluation pt. evalPt.

      ¤ cubicSpline*      : A pointer to an initiallized cubicSpline struct,
                            initiallized using cubicSpline_init()
      ¤ double*   evalPt  : Point at which to evaluate interpolant at

      Returns: The function value of the interpolant at evalPt.                      */
  //  ------------------------------------------------------------------------------
  assert( (evalPt >=  (spline -> pts[0])) && (evalPt <= (spline -> pts[spline -> numOfPts -1])) );

  //  Use a binary search to determine which subinterval evalPt is in
  int whichInterval  =  binarySearch(spline -> numOfPts, spline -> pts, evalPt);

  double ptsDiff     =  evalPt  - ( spline -> pts[whichInterval]                      );
  double thirdDiff   =  ptsDiff * ( spline -> thirdCoeff[whichInterval]               );
  double secondDiff  =  ptsDiff * ((spline -> secondCoeff[whichInterval]) + thirdDiff );
  double firstDiff   =  ptsDiff * ((spline -> firstCoeff[whichInterval])  + secondDiff);

  double interpVal = (spline -> funcVals[whichInterval]) + firstDiff ;
  return interpVal;
}

double cubicSplineDefiniteIntegral(cubicSpline *spline, int numOfPts, double* pts, double* funcVals, double evalPt ){
  int whichInterval = binarySearch( numOfPts, pts, evalPt ); // Find interval where z is (left)

  double integral = 0;
  for ( int intervalId = 0; intervalId <= whichInterval; intervalId ++) {

    double funcValDiff = (funcVals[intervalId + 1] - funcVals[intervalId]);
    double ptsDiff = (pts[intervalId + 1] - pts[intervalId]);
    double slope = funcValDiff / ptsDiff;

    if ( intervalId < whichInterval){
      integral += funcVals[intervalId] * ptsDiff + 1.0 / 2.0 * spline->firstCoeff[intervalId] * pow( ptsDiff, 2) + 1.0/3.0 * spline->secondCoeff[intervalId]*pow( ptsDiff, 3) + 1.0/4.0 * spline->thirdCoeff[intervalId]*pow( ptsDiff, 4) ;
    }

    else { //we are at the interval where z is
      integral += funcVals[intervalId] * ( evalPt - pts[intervalId]) + 1.0 / 2.0 * spline->firstCoeff[intervalId] * pow( ( evalPt - pts[intervalId]), 2 ) + 1.0/3.0 * spline->secondCoeff[intervalId] * pow( ( evalPt - pts[intervalId]), 3) + 1.0/4.0 * spline->thirdCoeff[intervalId]*pow( ( evalPt - pts[intervalId]), 4);
    }
  }

  return integral;
}

double cubicSplineDifferential(cubicSpline *spline, int numOfPts, double* pts, double evalPt ){ // evaluates s(evalPt)
   assert( ( evalPt >= (spline->pts[0]) ) && ( evalPt <= spline->pts[spline->numOfPts-1] ));

   int whichInterval = binarySearch ( numOfPts, pts, evalPt );
   double ptsDiff = evalPt - pts[whichInterval] ;
   double differential = (spline->firstCoeff[whichInterval] + 2.0 * ptsDiff * (spline->secondCoeff[whichInterval]) + 3.0 * pow( ptsDiff, 2) * (spline->thirdCoeff[whichInterval]) );

   return differential;
 }

void cubicSpline_free(cubicSpline* spline){
  //  ------------------------------------------------------------------------------
  /*  cubicSpline destructor, to free allocated memory.

      ¤ cubicSpline*       : A pointer to an initiallized cubicSpline struct,
                              initiallized using cubicSpline_init()                  */
  //  ------------------------------------------------------------------------------
  free(spline -> pts);
  free(spline -> funcVals);
  free(spline -> firstCoeff);
  free(spline -> secondCoeff);
  free(spline -> thirdCoeff);
  free(spline);
}
