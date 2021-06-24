#ifndef HAVE_CUBICSPLINE_H
#define HAVE_CUBICSPLINE_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

typedef struct { int      numOfPts    ; /* n   */
                 double*  pts         ; /* x's */
                 double*  funcVals    ; /* y's */
                 double*  firstCoeff  ; /* b_i */
                 double*  secondCoeff ; /* c_i */
                 double*  thirdCoeff  ; /* d_i */ } cubicSpline;

cubicSpline* cubicSplineAlloc( int numOfPts, double* pts, double* funcVals );
double cubicSplineEval( cubicSpline* spline, double evalPt );
double cubicSplineDefiniteIntegral(cubicSpline *spline, int numOfPts, double* pts, double* funcVals, double evalPt );
double cubicSplineDifferential(cubicSpline *spline, int numOfPts, double* pts, double evalPt ); // evaluates s(evalPt)
void cubicSpline_free(cubicSpline* spline);

#endif
