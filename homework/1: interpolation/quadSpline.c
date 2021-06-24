#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "quadSpline.h"
#include "binarySearch.h"

quadSpline* quadSplineAlloc( int numOfPts, double* pts, double* funcVals ){
  quadSpline* s = (quadSpline*)malloc(sizeof(quadSpline));

  // NumOfEqs is the number of needed equations for the quadratic spline interpolation,
  // or the number of intervals to interpolate in (between a numOfPts amount of points)
  int numOfEqs     =   numOfPts - 1;

  // Fill out field values
  s -> firstCoeff  = (double*)malloc( numOfEqs*sizeof( double ) );
  s -> secondCoeff = (double*)malloc( numOfEqs*sizeof( double ) );
  s -> pts         = (double*)malloc( numOfPts*sizeof( double ) );
  s -> funcVals    = (double*)malloc( numOfPts*sizeof( double ) );
  s -> numOfPts    = numOfPts;

  //int i;

  for ( int i = 0; i < numOfPts; i++ ){
    s -> pts[i]       =  pts[i];
    s -> funcVals[i]  =  funcVals[i];
  }

  double slope[numOfPts-1]; // p
  double ptsDiff[numOfPts-1] ;//VLA from C99 //h

  for( int i = 0; i < numOfPts-1; i++){
    ptsDiff[i]= pts[i+1]-pts[i] ;
    slope[i] = ( funcVals[i+1]-funcVals[i] ) / ptsDiff [i] ;
  }

  s->secondCoeff[0] = 0 ;// recursion  up;

  for(int i = 0; i < numOfPts-2; i++) {
    s->secondCoeff[i+1] = ( slope[i+1] - slope[i] - (s->secondCoeff[i])*ptsDiff[i] ) / ptsDiff[i+1];
     // recursion down :
  }
  s->secondCoeff[numOfPts-2] /= 2;

  for(int i = numOfPts-3; i >= 0; i--){
    s->secondCoeff[i] = ( slope[i+1] - slope[i] - (s->secondCoeff[i+1])*ptsDiff[i+1] ) / ptsDiff[i];
  }

  for(int i = 0; i < numOfPts-1; i++){
    s->firstCoeff[i] = slope[i] - (s->secondCoeff[i])*ptsDiff[i];
  }

  return s ;
}

double quadSplineEval(quadSpline *s , double evalPt ){ // evaluates s(evalPt)
   assert( ( evalPt >= (s->pts[0]) ) && ( evalPt <= s->pts[s->numOfPts-1] ));

   int whichInterval = binarySearch ( s->numOfPts, s->pts, evalPt );
   double ptsDiff = evalPt - ( s->pts[whichInterval] );
   double interpVal = s->funcVals[whichInterval] + ptsDiff*(s->firstCoeff[whichInterval] + ptsDiff*(s->secondCoeff[whichInterval]) );

   return interpVal;
 } // inerpolating  polynomial


 double quadSplineDefiniteIntegral(quadSpline *s, int numOfPts, double* pts, double* funcVals, double evalPt ){
   int whichInterval = binarySearch( numOfPts, pts, evalPt ); // Find interval where z is (left)

   double integral = 0;
   for ( int intervalId = 0; intervalId <= whichInterval; intervalId ++) {

     double funcValDiff = (funcVals[intervalId + 1] - funcVals[intervalId]);
     double ptsDiff = (pts[intervalId + 1] - pts[intervalId]);
     double slope = funcValDiff / ptsDiff;

     if ( intervalId < whichInterval){
       integral += funcVals[intervalId] * ptsDiff + 1.0 / 2.0 * s->firstCoeff[intervalId] * pow( ptsDiff, 2) + 1.0/3.0 * s->secondCoeff[intervalId]*pow( ptsDiff, 3);
     }

     else { //we are at the interval where z is
       integral += funcVals[intervalId] * ( evalPt - pts[intervalId]) + 1.0 / 2.0 * s->firstCoeff[intervalId] * pow( ( evalPt - pts[intervalId]), 2 ) + 1.0/3.0 * s->secondCoeff[intervalId] * pow( ( evalPt - pts[intervalId]), 3);
     }
   }

   return integral;
 }

 double quadSplineDifferential(quadSpline *s , double evalPt ){ // evaluates s(evalPt)
    assert( ( evalPt >= (s->pts[0]) ) && ( evalPt <= s->pts[s->numOfPts-1] ));

    int whichInterval = binarySearch ( s->numOfPts, s->pts, evalPt );
    double ptsDiff = evalPt - ( s->pts[whichInterval] );
    double differential = (s->firstCoeff[whichInterval] + 2.0 * ptsDiff * (s->secondCoeff[whichInterval]) );

    return differential;
  }

 void quadSplinefree(quadSpline *s ){
   // free the allocated memory
   free( s->pts ) ;
   free( s->funcVals ) ;
   free( s->firstCoeff ) ;
   free( s->secondCoeff ) ;
   free( s ) ;
 }
