#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MCintegration.h"

int main(int argc, char const *argv[]) {

  printf("_________________________ PARTS A), B) AND C) _________________________\n");


  // Simple func for debugging___________________________________________________
  int dimDebug = 1;
  int numOfPtsDebug = 1e4;
  double* lowerBoundDebug = malloc(sizeof(double)*dimDebug);
  double* upperBoundDebug = malloc(sizeof(double)*dimDebug);
  double resultDebugRand;
  double errorDebugRand;
  double resultDebugHC;
  double errorDebugHC;
  double resultDebugStrat;
  double absAcc     = 0.005;
  double relAcc     = 0;
  int numOfRecalls  = 0;
  double meanRecall = 0;


  resultDebugRand  = 0;
  errorDebugRand   = 0;
  resultDebugHC    = 0;
  errorDebugHC     = 0;
  resultDebugStrat = 0;

  lowerBoundDebug[0] = 0;
  upperBoundDebug[0] = 1;

  double funcDebug( double* x ){
    return sqrt( x[0] );
  }

  plainMC(dimDebug, lowerBoundDebug, upperBoundDebug, funcDebug, numOfPtsDebug, &resultDebugRand, &errorDebugRand);
  HaltonCorputMC(dimDebug, lowerBoundDebug, upperBoundDebug, funcDebug, numOfPtsDebug, &resultDebugHC, &errorDebugHC );
  resultDebugStrat = stratMC(dimDebug, funcDebug, lowerBoundDebug, upperBoundDebug, absAcc, relAcc, numOfRecalls, meanRecall);




  printf("∫_0^1 √(x) dx = 2/3 = %g\n", 2.0/3.0);
  printf("plain MC integration of this yields %g\n", resultDebugRand);
  printf("with error = %g\n\n", errorDebugRand);
  printf("Halton-Corput MC integration of this yields %g\n", resultDebugHC);
  printf("with error = %g\n\n", errorDebugHC);
  printf("Stratified MC integration of this yields %g\n", resultDebugStrat);




  printf("---------------------------------\n\n");


  // Test func from problem_______________________________________________________
  int dimTest = 3;
  int numOfPtsTest = 1e6;
  double* lowerBoundTest = malloc(sizeof(double)*dimTest);
  double* upperBoundTest = malloc(sizeof(double)*dimTest);
  double resultTestRand;
  double errorTestRand;
  double resultTestHC;
  double errorTestHC;
  double resultTestStrat;


  resultTestRand  = 0;
  errorTestRand   = 0;
  resultTestHC    = 0;
  errorTestHC     = 0;
  resultTestStrat = 0;

  lowerBoundTest[0] = 0;
  lowerBoundTest[1] = 0;
  lowerBoundTest[2] = 0;

  upperBoundTest[0] = M_PI;
  upperBoundTest[1] = M_PI;
  upperBoundTest[2] = M_PI;


  double funcTest( double* x ){
    double result = 1/((1 - cos(x[0])*cos(x[1])*cos(x[2]))* M_PI*M_PI*M_PI);
    return result;
  }
  double newTestFunc(double *x){
    double firstVal    =   x[0];
	  double secondVal   =   x[1];
    return cos(firstVal) * cos(secondVal);
  }

  plainMC(dimTest, lowerBoundTest, upperBoundTest, funcTest, numOfPtsTest, &resultTestRand, &errorTestRand );
  HaltonCorputMC(dimTest, lowerBoundTest, upperBoundTest, funcTest, numOfPtsTest, &resultTestHC, &errorTestHC);
  resultTestStrat = stratMC(2, newTestFunc, lowerBoundTest, upperBoundTest, absAcc, relAcc, numOfRecalls, meanRecall);

  printf("∫0π  dx/π ∫0π  dy/π ∫0π  dz/π [1-cos(x)cos(y)cos(z)]^-1 = 1.3932039296856768591842462603255\n");
  printf("plain MC integration of this yields %.25g\n", resultTestRand);
  printf("with error = %g\n\n", errorTestRand);
  printf("Halton-Corput MC integration of this yields %g\n", resultTestHC);
  printf("with error = %g\n\n", errorTestHC);
  printf("Stratified MC integration of this yields %g\n", resultTestStrat);


// Estimation of errors
  double resultRand;
  double errorRand;
  double resultHC;
  double errorHC;

  resultRand = 0;
  errorRand  = 0;
  resultHC   = 0;
  errorHC    = 0;

  FILE* outputStream = fopen("errors.txt","w");
  for (int numOfPts = 10; numOfPts < 1e5; numOfPts+=1000) {

    plainMC(dimTest, lowerBoundTest, upperBoundTest, funcTest, numOfPts, &resultTestRand, &errorRand);
    HaltonCorputMC(dimTest, lowerBoundTest, upperBoundTest, funcTest, numOfPts, &resultTestHC, &errorHC);

    fprintf(outputStream, "%d\t%g\t%g\n",numOfPts,errorRand, errorHC);
  }
  fclose(outputStream);
  return 0;
}
