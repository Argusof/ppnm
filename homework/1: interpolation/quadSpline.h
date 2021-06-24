#ifndef HAVE_QUADSPINE_H
#define HAVE_QUADSPINE_H

typedef struct { int      numOfPts    ; /* n   */
                 double*  pts         ; /* x's */
                 double*  funcVals    ; /* y's */
                 double*  firstCoeff  ; /* b_i */
                 double*  secondCoeff ; /* c_i */ } quadSpline;

quadSpline* quadSplineAlloc( int numOfPts, double* pts, double* funcVals );
double quadSplineEval(quadSpline *s , double z );
double quadSplineDefiniteIntegral(quadSpline *s, int numOfPts, double* pts, double* funcVals, double evalPt );
double quadSplineDifferential(quadSpline *s , double evalPt );
void quadSplinefree(quadSpline *s );

#endif
