#ifndef HAVE_LINSPINE_H
#define HAVE_LINSPINE_H

double linSplineInterp( int numOfPts, double* pts, double* funcVals, double evalPt);
//double linSplineInterpIntegral( int numOfPts, double* pts, double* funcVals, double evalPt);
double linSplineDefiniteIntegral( int numOfPts, double* pts, double* funcVals, double evalPt);

#endif
