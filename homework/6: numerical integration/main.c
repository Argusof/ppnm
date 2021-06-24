#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#include "integration.h"

void printTestResults(char* string, double integralVal, double absAcc, double relAcc, double exactVal, double integrationError, int calls){
  printf("\n%s  : %.25g\n",string, integralVal);
  printf("Error goal               : %.25g\n", absAcc + fabs(exactVal) * relAcc);
  printf("Actual error             : %.25g\n", fabs(integralVal - exactVal));
  printf("Computed error estimate  : %.25g\n", integrationError);
  printf("Called function %d times\n\n",calls);
}

int main(int argc, char const *argv[]) {

  // PART A_______________________________________________________________________________________
  double leftEndPt  = 0;
  double rightEndPt = 1;
  double absAcc     = 1e-3;
  double relAcc     = 1e-3;
  int calls         = 0;
  double integrationError = 0;

  double firstTestFunc( double x ){
    calls++;
    return sqrt( x );
  }

  printf("_______________________________________________________________________\n" );
  printf("\nA) Recursive adaptive integrator (Error extimate from part C) included)\n");
  printf("Testing routine on different integrals; \n\n");
  printf("∫_0^1 √(x) dx = 2/3 = %g\n",2.0/3.0);

  double integralVal = integrate(firstTestFunc, leftEndPt, rightEndPt, absAcc, relAcc, &integrationError);
  double exactVal    = 2.0/3.0;
  printTestResults("Numerical integration yields",  integralVal, absAcc, relAcc, exactVal, integrationError, calls);


  calls = 0;
  integrationError = 0;
  double secondTestFunc( double x ){
    calls++;
    return 4 * sqrt(1 - x*x);
  }
  exactVal = M_PI;
  printf("∫_0^1 4√(1-x²) dx = π = %g\n",exactVal);
  integralVal = integrate(secondTestFunc, leftEndPt, rightEndPt, absAcc, relAcc, &integrationError);
  printTestResults("Numerical integration yields",  integralVal, absAcc, relAcc, exactVal, integrationError, calls);


  // PART B_______________________________________________________________________________________
  printf("_______________________________________________________________________\n" );
  printf("\nB) Open quadrature with Clenshaw–Curtis variable transformation \n");
  printf("Tested on different integrals:\n\n");

  calls = 0;
  integrationError = 0;
  double thirdTestFunc( double x ){
    calls++;
    return 1/sqrt(x);
  }
  integralVal = openQuad(thirdTestFunc, leftEndPt, rightEndPt,absAcc, relAcc, &integrationError);
  exactVal    = 2.0;
  printf("\n∫_0^1 1/√(x) dx = %g\n\n", exactVal);
  printTestResults("With Clenshaw Curtis",  integralVal, absAcc, relAcc, exactVal, integrationError, calls);


  calls = 0;
  integrationError = 0;
  integralVal = integrate(thirdTestFunc, leftEndPt, rightEndPt,absAcc, relAcc, &integrationError);
  printTestResults("Without Clenshaw Curtis",  integralVal, absAcc, relAcc, exactVal, integrationError, calls);


  printf("\n--------------------------------------------\n");
  calls = 0;
  integrationError = 0;
  double fourthTestFunc( double x ){
    calls++;
    return log(x)/sqrt(x);
  }
  integralVal = openQuad(fourthTestFunc, leftEndPt, rightEndPt,absAcc, relAcc, &integrationError);
  exactVal    = -4.0;
  printf("\n∫_0^1 ln(x)/√(x) dx = -4 \n\n");
  printTestResults("With Clenshaw Curtis",  integralVal, absAcc, relAcc, exactVal, integrationError, calls);


  calls = 0;
  integrationError = 0;
  integralVal  = integrate(fourthTestFunc, leftEndPt, rightEndPt,absAcc, relAcc, &integrationError);
  printTestResults("Without Clenshaw Curtis",  integralVal, absAcc, relAcc, exactVal, integrationError, calls);


  printf("--------------------------------------------\n");
  calls = 0;
  integrationError = 0;
  integralVal = openQuad(secondTestFunc, leftEndPt, rightEndPt,absAcc, relAcc, &integrationError);
  exactVal     = M_PI;

  printf("\n∫_0^1 4√(1-x²) dx = π = %g\n\n",exactVal);
  printTestResults("With Clenshaw Curtis",  integralVal, absAcc, relAcc, exactVal, integrationError, calls);

  calls =0;
  integrationError = 0;
  integralVal  = integrate(secondTestFunc, leftEndPt, rightEndPt,absAcc, relAcc, &integrationError);
  printTestResults("Without Clenshaw Curtis",  integralVal, absAcc, relAcc, exactVal, integrationError, calls);


  double gsl_test_func( double x, void* params){
    params = NULL;
    return secondTestFunc(x);
  }
  gsl_function gslTestFunc;
  gslTestFunc.function = &gsl_test_func;
  gslTestFunc.params = NULL;
  size_t limit = 999;
  double result;
  double absError;
  size_t numOfEvals;
  calls = 0;
  integrationError = 0;
  gsl_integration_cquad_workspace* workspace = gsl_integration_cquad_workspace_alloc(limit);
  integralVal = gsl_integration_cquad(&gslTestFunc, leftEndPt, rightEndPt, absAcc, relAcc, workspace, &result, &absError, &numOfEvals);
  exactVal     = M_PI;
  printTestResults("With gsl Clenshaw Curtis",  result, absAcc, relAcc, exactVal, absError, (int)numOfEvals);



// PART C_________________________________________________________________________________________
  printf("_______________________________________________________________________\n" );
  printf("\nC) Infinite limits \n");
  printf("Tested on different integrals:\n\n");

  leftEndPt = -INFINITY;
  rightEndPt = INFINITY;
  double sixthTestFunc( double x ){
    calls++;
    return exp(-x*x);
  }
  calls =0;
  integrationError = 0;
  exactVal = sqrt(M_PI);
  printf("\n∫_-inf^inf exp(-x²) dx = √π = %g",exactVal);
  integralVal = integrate(sixthTestFunc, leftEndPt, rightEndPt,absAcc, relAcc, &integrationError);
  printTestResults("Numerical integration yields",  integralVal, absAcc, relAcc, exactVal, integrationError, calls);

  double gsl_test_func6( double x, void* params){
    params = NULL;
    return sixthTestFunc(x);
  }
  gsl_function gslTestFunc6;
  gslTestFunc6.function = &gsl_test_func6;
  gslTestFunc6.params = NULL;
  result = 0;
  absError = 0;
  numOfEvals = 0;
  gsl_integration_workspace* workspace2 = gsl_integration_workspace_alloc(limit);
  integralVal = gsl_integration_qagi(&gslTestFunc6, absAcc, relAcc, limit, workspace2, &result, &absError);
  printTestResults("With gsl this yields",  result, absAcc, relAcc, exactVal, absError, (int)numOfEvals);
  printf("gsl routine cannot return number of calls\n" );
  printf("--------------------------------------------\n");


  leftEndPt = 0;
  rightEndPt = INFINITY;
  double seventhTestFunc( double x ){
    calls++;
    return 1/(1+x*x);
  }
  calls = 0;
  integrationError = 0;
  exactVal = M_PI/2;
  printf("\n∫_0^inf 1/(1+x²) dx = π/2 = %g",exactVal);
  integralVal = integrate(seventhTestFunc, leftEndPt, rightEndPt,absAcc, relAcc, &integrationError);
  printTestResults("Numerical integration yields",  integralVal, absAcc, relAcc, exactVal, integrationError, calls);

  double gsl_test_func7( double x, void* params){
    params = NULL;
    return seventhTestFunc(x);
  }
  gsl_function gslTestFunc7;
  gslTestFunc7.function = &gsl_test_func7;
  gslTestFunc7.params = NULL;
  result = 0;
  absError = 0;
  numOfEvals = 0;
  gsl_integration_workspace* workspace3 = gsl_integration_workspace_alloc(limit);
  integralVal = gsl_integration_qagiu(&gslTestFunc7, leftEndPt, absAcc, relAcc, limit, workspace3, &result, &absError);
  printTestResults("With gsl this yields",  result, absAcc, relAcc, exactVal, absError, (int)numOfEvals);
  printf("gsl routine cannot return number of calls\n");
  printf("--------------------------------------------\n");


  leftEndPt = -INFINITY;
  rightEndPt = 0;
  calls = 0;
  integrationError = 0;
  exactVal = M_PI/2;
  printf("\n∫_-inf^0 1/(1+x²) dx = π/2 = %g",exactVal);
  integralVal = integrate(seventhTestFunc, leftEndPt, rightEndPt,absAcc, relAcc, &integrationError);
  printTestResults("Numerical integration yields",  integralVal, absAcc, relAcc, exactVal, integrationError, calls);

  double gsl_test_func72( double x, void* params){
    params = NULL;
    return seventhTestFunc(x);
  }
  gsl_function gslTestFunc72;
  gslTestFunc72.function = &gsl_test_func72;
  gslTestFunc72.params = NULL;
  result = 0;
  absError = 0;
  numOfEvals = 0;
  gsl_integration_workspace* workspace4 = gsl_integration_workspace_alloc(limit);
  integralVal = gsl_integration_qagil(&gslTestFunc72, rightEndPt, absAcc, relAcc, limit, workspace4, &result, &absError);
  printTestResults("With gsl this yields",  result, absAcc, relAcc, exactVal, absError, (int)numOfEvals);
  printf("gsl routine cannot return number of calls\n");


  return 0;
}
