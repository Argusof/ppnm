#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <math.h>
#include "lnXsqrtX.h"
#include "erfint.h"
#include "besselint.h"

int main( int argc, char* argv[] ) {
  FILE* myOutputFileStream_error =  fopen("error.txt", "w");                                  // open output file
  if ( argc < 2 ){                                                                            // check whether the file is open
    fprintf( myOutputFileStream_error, "Error! The program did not receive any arguments.");
  }
  else {
    FILE* myOutputFileStream_lnxSqrtx	  =  fopen(argv[1], "w");
    FILE* myOutputFileStream_erfInt     =  fopen(argv[2], "w");
    FILE* myOutputFileStream_besselInt  =  fopen(argv[3], "w");

    // ------------ Exercise A: logxSqrtx ------------ //
    double lowerLimit       =   0       ;
    double upperLimit       =   1       ;
    double absError         =   1e-6    ;
    double relError         =   1e-6    ;
    size_t iterationLimit   =   999     ;

    gsl_function gslFunctionLogxSqrtx;

    gslFunctionLogxSqrtx.function  =  &logxSqrtx;
    gslFunctionLogxSqrtx.params    =  NULL;

    double result  =  integratedFunction( lowerLimit, upperLimit, &gslFunctionLogxSqrtx, absError, relError, iterationLimit );
    fprintf( myOutputFileStream_lnxSqrtx, "A) The integral of ln(x)/sqrt(x) = %g (numerically)\n", result);



    // ------------ Exercise B: Error function ------------ //
    gsl_function gslFunctionErfInt;
    gslFunctionErfInt.function  =  &erfIntegral;
    gslFunctionErfInt.params    =  NULL;

    double xMax  =  2			;
    double xMin  =  -xMax	;

    for( double val = xMin; val <= xMax; val += 1.0/8 ){
      double result  =  integratedErrorFunction( val, &gslFunctionErfInt, absError, relError, iterationLimit );
      fprintf(myOutputFileStream_erfInt, "%10g %10g %10g\n", val, result, erf(val));
    }



    // ------------ Exercise C: Bessel function ------------
    gsl_function gslFunctionBesselInt_n0;
    gsl_function gslFunctionBesselInt_n1;
    gsl_function gslFunctionBesselInt_n2;

    gslFunctionBesselInt_n0.function  =  &besselIntegral;
    gslFunctionBesselInt_n1.function  =  &besselIntegral;
    gslFunctionBesselInt_n2.function  =  &besselIntegral;

    xMax  =  20	;   // lower bound
    xMin  =  0  ;   // upper bound

    for( double val = xMin; val <= xMax; val += 1.0/8 ){
      besselParams params_n0 =  { 0, val };
      besselParams params_n1 =  { 1, val };
      besselParams params_n2 =  { 2, val };
      gslFunctionBesselInt_n0.params    =  &params_n0;
      gslFunctionBesselInt_n1.params    =  &params_n1;
      gslFunctionBesselInt_n2.params    =  &params_n2;

      double result_n0  =  integratedBesselFunction( &gslFunctionBesselInt_n0, absError, relError, iterationLimit );
      double result_n1  =  integratedBesselFunction( &gslFunctionBesselInt_n1, absError, relError, iterationLimit );
      double result_n2  =  integratedBesselFunction( &gslFunctionBesselInt_n2, absError, relError, iterationLimit );

      fprintf(myOutputFileStream_besselInt, "%10g %10g %10g %10g\n", val, result_n0, result_n1, result_n2);
    }

    return 0;
  }
}
