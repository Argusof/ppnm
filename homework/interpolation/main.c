#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <math.h>
#include "linSpline.h"
#include "binarySearch.h"
#include "gslIntFunc.h"
#include "quadSpline.h"
#include "cubicSpline.h"

void fillArrayFromInputData( double* dataX,  double* dataY, char* inputFileName);

int main( int argc, char* argv[] ) {
  int numOfPts = 20;
  int numOfQueryPts = 1000;
  double evalPt = 4.3;

  double* dataX = malloc(numOfPts*sizeof(double));
  double* dataY = malloc(numOfPts*sizeof(double));

  if (argc < 2) {
    fprintf(stderr, "Error: no input arguments\n" );
  }
  else {
    char* inputFileName = argv[1];
    fillArrayFromInputData( dataX, dataY, inputFileName);
  }

  // ____________________________ LINEAR INTERPOLATION ____________________________ //

  FILE* OutputDataFileLin   = fopen( argv[2], "w");

  double linInterpValTmp    = 0.0;
  double linInterpValIntTmp = 0.0;
  double gslLinInterpTmp    = 0.0;
  double gslLinIntegTmp     = 0.0;
  double resolution         = fabs ( dataX[numOfPts - 1] - dataX[0] ) / numOfQueryPts;

  gsl_interp* linear = gsl_interp_alloc( gsl_interp_linear, numOfPts);
  gsl_interp_init( linear, dataX, dataY, numOfPts);


  for (double evalPt = dataX[0] ; evalPt < dataX[numOfPts]; evalPt += resolution) {
    linInterpValTmp    = linSplineInterp( numOfPts, dataX, dataY, evalPt);
    linInterpValIntTmp = linSplineDefiniteIntegral( numOfPts, dataX, dataY, evalPt);

    gslLinInterpTmp    = gsl_interp_eval( linear, dataX, dataY, evalPt, NULL );
    gslLinIntegTmp     = gsl_interp_eval_integ( linear, dataX, dataY, dataX[0], evalPt ,NULL );

    fprintf(OutputDataFileLin, "%g\t%g\t%g\t%g\t%g\n", evalPt, linInterpValTmp, linInterpValIntTmp, gslLinInterpTmp, gslLinIntegTmp);
  }


  double integral = linSplineDefiniteIntegral( numOfPts, dataX, dataY, evalPt);
  printf("\nIntegrating the function from 0 to %g gives %g\n", evalPt, integral);

  double lowerLimit     = 0;
  double upperLimit     = evalPt;
  double absError       = 1e-6;
  double relError       = 1e-6;
  size_t iterationLimit = 999;

  gsl_function gslFuncCos;
  gslFuncCos.function = &cos;
  gslFuncCos.params   = NULL;

  double gslIntegral = integratedFunction( lowerLimit, upperLimit, &gslFuncCos, absError, relError, iterationLimit );
  printf("Using gsl this yields %g\n\n", gslIntegral);


  // ____________________________ QUADRADIC INTERPOLATION ____________________________ //
  FILE* OutputDataFileQuad   = fopen( argv[3], "w");

  double quadInterpValTmp     = 0.0;
  double quadInterpValIntTmp  = 0.0;
  double quadInterpValDiffTmp = 0.0;
  double gslQuadInterpTmp     = 0.0;
  double gslQuadIntegTmp      = 0.0;
  double gslQuadDiffTmp       = 0.0;
  quadSpline* quad            = quadSplineAlloc( numOfPts, dataX, dataY );


  gsl_interp* gslQuad    = gsl_interp_alloc(gsl_interp_cspline, numOfPts);
  gsl_interp_init(gslQuad, dataX, dataY, numOfPts);

  for (double evalPt = dataX[0] ; evalPt < dataX[numOfPts]; evalPt += resolution) {
    quadInterpValTmp     = quadSplineEval(quad, evalPt);
    quadInterpValIntTmp  = quadSplineDefiniteIntegral(quad, numOfPts, dataX , dataY , evalPt );
    quadInterpValDiffTmp = quadSplineDifferential(quad, evalPt );

    gslQuadInterpTmp     = gsl_interp_eval(gslQuad , dataX, dataY, evalPt, NULL);
    gslQuadIntegTmp      = gsl_interp_eval_integ(gslQuad, dataX, dataY, dataX[0], evalPt, NULL);
    gslQuadDiffTmp       = gsl_interp_eval_deriv(gslQuad , dataX, dataY, evalPt, NULL);

    fprintf(OutputDataFileQuad, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", evalPt, quadInterpValTmp, quadInterpValIntTmp, quadInterpValDiffTmp,gslQuadInterpTmp, gslQuadIntegTmp, gslQuadDiffTmp);
  }

  // ____________________________ CUBIC INTERPOLATION ____________________________ //
  FILE* OutputDataFileCubic   = fopen( argv[4], "w");

  double cubicInterpValTmp     = 0.0;
  double cubicInterpValIntTmp  = 0.0;
  double cubicInterpValDiffTmp = 0.0;
  double gslCubicInterpTmp     = 0.0;
  double gslCubicIntegTmp      = 0.0;
  double gslCubicDiffTmp       = 0.0;
  cubicSpline* cubic           = cubicSplineAlloc( numOfPts, dataX, dataY );

  gsl_interp* gslCubic  = gsl_interp_alloc(gsl_interp_cspline, numOfPts);
  gsl_interp_init(gslCubic, dataX, dataY, numOfPts);

  for (double evalPt = dataX[0] ; evalPt < dataX[numOfPts]; evalPt += resolution) {
    cubicInterpValTmp     = cubicSplineEval(cubic, evalPt);
    cubicInterpValIntTmp  = cubicSplineDefiniteIntegral(cubic, numOfPts, dataX , dataY , evalPt );
    cubicInterpValDiffTmp = cubicSplineDifferential( cubic, numOfPts, dataX, evalPt );

    gslCubicInterpTmp     = gsl_interp_eval(gslCubic , dataX, dataY, evalPt, NULL);
    gslCubicIntegTmp      = gsl_interp_eval_integ(gslCubic, dataX, dataY, dataX[0], evalPt, NULL);
    gslCubicDiffTmp       = gsl_interp_eval_deriv(gslCubic , dataX, dataY, evalPt, NULL);

    fprintf(OutputDataFileCubic, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", evalPt, cubicInterpValTmp, cubicInterpValIntTmp, cubicInterpValDiffTmp, gslCubicInterpTmp, gslCubicIntegTmp, gslCubicDiffTmp);
  }

  // ____________________________ FREE MEMORY ____________________________ //

  gsl_interp_free(linear);
  quadSplinefree(quad);
  cubicSpline_free(cubic);
  free(dataX);
  free(dataY);

  return 0;
}
