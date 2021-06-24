#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "rootfinding.h"
#include "rkstep.h"

int writeToFile;

double energy;
double bound_energy;
double maxPt;
double bound_maxPt;

void testFunc(gsl_vector* vals, gsl_vector* funcVals){ //1+(x-a)^2+(y-b)^2
  double x = gsl_vector_get(vals, 0);
  double y = gsl_vector_get(vals, 1);

  double a = 6.0;
  double b = 13.0;

  double gradientX = 2*x*(x-a);
  double gradientY = 2*y*(y-b);

  gsl_vector_set(funcVals, 0, gradientX);
  gsl_vector_set(funcVals, 1, gradientY);
}

void RosenbrockValley(gsl_vector* vals, gsl_vector* funcVals){
  double gradientX = -2*(1-gsl_vector_get(vals,0)) + (-2*gsl_vector_get(vals,0))*2*100*(gsl_vector_get(vals,1)-gsl_vector_get(vals,0)*gsl_vector_get(vals,0));
  double gradientY = 2*100*(gsl_vector_get(vals,1)-gsl_vector_get(vals,0)*gsl_vector_get(vals,0));

  gsl_vector_set(funcVals, 0, gradientX);
  gsl_vector_set(funcVals, 1, gradientY);
}
void schrodingerEq (double var, gsl_vector* funcVal, gsl_vector* funcDeriv){
    double thisFuncVal = gsl_vector_get(funcVal, 0);
    double firstDeriv = gsl_vector_get(funcVal, 1);
    double secondDeriv = (-2)*(1.0/var + energy)*thisFuncVal;

    gsl_vector_set(funcDeriv, 0, firstDeriv);
    gsl_vector_set(funcDeriv, 1, secondDeriv);
}

void schrodingerEq_bound (double var, gsl_vector* funcVal, gsl_vector* funcDeriv){
    double thisFuncVal = gsl_vector_get(funcVal, 0);
    double firstDeriv = gsl_vector_get(funcVal, 1);
    double secondDeriv = (-2)*(1.0/var + energy)*thisFuncVal;

    gsl_vector_set(funcDeriv, 0, firstDeriv);
    gsl_vector_set(funcDeriv, 1, secondDeriv);
}

void wavefunc(gsl_vector* vals, gsl_vector* funcVals){
    int dim  =  2;                                          // The order of the harmonic function
    gsl_vector* funcValRight  =  gsl_vector_alloc(dim);     // Vector to hold final value
    gsl_vector* funcValLeft   =  gsl_vector_calloc(dim);    // Vector to hold initial value

    // The ODE is defined on the interval [ leftEndpt ,  rightEndpt ]
    double  leftEndpt   =   1e-3;                               // ODE diverges at origin, so we choose small values
    double  rightEndpt  =   maxPt;
    double  absAcc      =   1e-3;                               // Absolute accuracy
    double  relAcc      =   1e-3;                               // Relative accuracy
    double  step        =   (rightEndpt - leftEndpt) / 10;      // Initial stepsize
    energy              =   gsl_vector_get(vals, 0);

    gsl_vector_set(funcValLeft, 0, (leftEndpt - leftEndpt*leftEndpt));
    gsl_vector_set(funcValLeft, 1, (1 - 2*leftEndpt));

    rkdriver( schrodingerEq, leftEndpt, funcValLeft, rightEndpt, funcValRight, step, absAcc, relAcc, NULL);
    gsl_vector_set(funcVals, 0, gsl_vector_get(funcValRight,0) );
}

void wavefunc_bound(gsl_vector* vals, gsl_vector* funcVals){
    int dim  =  2;                                          // The order of the harmonic function
    gsl_vector* funcValRight  =  gsl_vector_alloc(dim);     // Vector to hold final value
    gsl_vector* funcValLeft   =  gsl_vector_calloc(dim);    // Vector to hold initial value

    // The ODE is defined on the interval [ leftEndpt ,  rightEndpt ]
    double  leftEndpt   =   1e-3;
    double  rightEndpt  =   bound_maxPt;
    double  absAcc      =   1e-3;                               // Absolute accuracy
    double  relAcc      =   1e-3;                               // Relative accuracy
    double  step        =   (rightEndpt - leftEndpt) / 10;      // Initial stepsize
    bound_energy        =   gsl_vector_get(vals, 0);

    gsl_vector_set(funcValLeft, 0, (leftEndpt - leftEndpt*leftEndpt));
    gsl_vector_set(funcValLeft, 1, (1 - 2.*leftEndpt));

    rkdriver( schrodingerEq_bound, leftEndpt, funcValLeft, rightEndpt, funcValRight, step, absAcc, relAcc, NULL);
    gsl_vector_set(funcVals, 0, gsl_vector_get(funcValRight, 0) - rightEndpt*exp(-sqrt((-2)*bound_energy)*rightEndpt));
}


int main(int argc, char const *argv[]) {

printf("\n---------------------PART A---------------------\n");
printf("Finding roots of 1+(x-a)^2+(y-b)^2\n\n");

  int numOfDims    = 2;
  double tolerance = 1e-5;
  gsl_vector* minimum = gsl_vector_calloc(numOfDims);

  gsl_vector_set(minimum, 0, 4);
  gsl_vector_set(minimum, 1, 10);

  newtonRaphsonMethod(testFunc, minimum, tolerance);

  printf("The initial value is (x, y) = (4, 10)\n");
  printf("The minimum should be at (x, y) = (a, b) = (6, 13)\n");
  printf("The minimum is found at (x, y) = (%g, %g),\n\n\n", gsl_vector_get(minimum,0), gsl_vector_get(minimum,1));

  printf("Finding roots of the Rosenbrock valley function\n\n");


  gsl_vector_set(minimum, 0, 0.5);
  gsl_vector_set(minimum, 1, 0.5);

  newtonRaphsonMethod(RosenbrockValley, minimum, tolerance);

  printf("The initial value is (x, y) = (0.5, 0.5)\n");
  printf("The minimum is found at (x, y) = (%g, %g)\n", gsl_vector_get(minimum,0), gsl_vector_get(minimum,1));



  printf("\n---------------------PART B---------------------\n");
  printf("(Un)bound states of the hydrogen atom\n\n");

  gsl_vector* minimum_hydrogen       = gsl_vector_alloc(1);
  gsl_vector* minimum_hydrogen_bound = gsl_vector_alloc(1);
  gsl_vector_set(minimum_hydrogen, 0,  -3);
  gsl_vector_set(minimum_hydrogen_bound, 0, -1);

  maxPt = 8.0;
  bound_maxPt = 0.5;
  newtonRaphsonMethod(wavefunc, minimum_hydrogen, tolerance);

  newtonRaphsonMethod(wavefunc_bound, minimum_hydrogen_bound, tolerance);

  printf("Energies found by root-finding\n");
  printf("   Unbound: %g\n", energy);


  int dim  =  2;                                          // The order of the harmonic function
  gsl_vector* funcValRight  =  gsl_vector_alloc(dim);     // Vector to hold final value
  gsl_vector* funcValLeft   =  gsl_vector_calloc(dim);    // Vector to hold initial value

  // The ODE is defined on the interval [ leftEndpt ,  rightEndpt ]
  double  leftEndpt   =   1e-3;
  double  rightEndpt  =   maxPt;
  double  absAcc      =   1e-3;                               // Absolute accuracy
  double  relAcc      =   1e-3;                               // Relative accuracy
  double  step        =   (rightEndpt - leftEndpt) / 10;      // Initial stepsize

  gsl_vector_set(funcValLeft, 0, (leftEndpt - leftEndpt*leftEndpt));
  gsl_vector_set(funcValLeft, 1, (1 - 2.*leftEndpt));

  FILE* filePath = fopen("hydrogen.txt", "w"); // Set up filestream to write ODE solution to
  rkdriver(schrodingerEq, leftEndpt, funcValLeft, rightEndpt, funcValRight, step, absAcc, relAcc, filePath);
  fclose(filePath);



  printf("   Bound: %g\n", bound_energy);

  printf("\n---------------------PART C---------------------\n");

  printf("Better boundary conditions for the hydrogen atom\n\n");


  int numOfSteps = 100;
  double exactEnergy = -0.5;

  FILE* convData = fopen("convData.txt", "w");
  for (int i = 1; i < numOfSteps; i++) {
    maxPt = (double)i/20.0;
    bound_maxPt = (double)i/20.0;
    gsl_vector_set(minimum_hydrogen, 0,  -3);
    gsl_vector_set(minimum_hydrogen_bound, 0, -1);

    newtonRaphsonMethod(wavefunc, minimum_hydrogen, tolerance);
    newtonRaphsonMethod(wavefunc_bound, minimum_hydrogen_bound, tolerance);

    fprintf(convData, "%g\t%g\t%g\n", maxPt, fabs(energy-exactEnergy),fabs(bound_energy-exactEnergy));
  }

  printf("See conv.png\n");


  gsl_vector_free(minimum);
  gsl_vector_free(minimum_hydrogen);
  gsl_vector_free(minimum_hydrogen_bound);

  return 0;

}
