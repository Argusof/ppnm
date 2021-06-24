#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "minimization.h"
#include "input2Array.h"

int numOfDataPts;
double* energy;
double* crossSection;
double* error;


double testFunc(gsl_vector* vals){ //1+(x-a)^2+(y-b)^2
  double x = gsl_vector_get(vals, 0);
  double y = gsl_vector_get(vals, 1);
  double a = 6.0;
  double b = 13.0;

  return 1 + pow((x-a),2) + pow((y-b),2);
}

double simplexTestFunc(double* vals){ //1+(x-a)^2+(y-b)^2
  double x = vals[0];
  double y = vals[1];
  double a = 6.0;
  double b = 13.0;

  return 1 + pow((x-a),2) + pow((y-b),2);
}

double RosenbrockValleyFunc(gsl_vector* vals){
  double x = gsl_vector_get(vals, 0);
  double y = gsl_vector_get(vals, 1);

  return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
}

double HimmelblauFunc(gsl_vector* vals){
  double x = gsl_vector_get(vals, 0);
  double y = gsl_vector_get(vals, 1);

  return (x*x+y-11)*(x*x+y-11) + (x+y*y-7)*(x+y*y-7);
}

double BreitWignerFunc(double mass, double width, double scaleFactor, double BWenergy){
  return scaleFactor/( ((BWenergy-mass)*(BWenergy-mass)) + (width*width)/4 );
}

double deviationFunc(gsl_vector* vals){
  double mass        = gsl_vector_get(vals, 0);
  double width       = gsl_vector_get(vals, 1);
  double scaleFactor = gsl_vector_get(vals, 2);
  double sum = 0;

  for (int dataPt = 0; dataPt < numOfDataPts; dataPt++) {
    sum += pow(BreitWignerFunc(mass, width, scaleFactor, energy[dataPt]) - crossSection[dataPt],2) / pow(error[dataPt],2);
  }

  return sum;
}


int main(int argc, char const *argv[]) {

  printf("--------------------- PART A ---------------------\n");
  printf("Testing minimization routine on 1+(x-a)^2+(y-b)^2 \n");
  int dims = 2;
  double tolerance = 1e-5;
  gsl_vector* minimum = gsl_vector_alloc(dims);
  gsl_vector_set(minimum, 0, 4);
  gsl_vector_set(minimum, 1, 10);

  printf("The initial value was (x,y) = (4,10)\n");
  printf("and actual the minimum is at (x,y) = (6,13)\n");

  quasiNewtonMethod(testFunc, minimum, tolerance);
  printf("The found minimum is (%g,%g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
  printf("____________________________________________________\n");


  printf("Testing minimization routine on the Rosenbrock valley function\n");
  gsl_vector_set(minimum, 0, 0);
  gsl_vector_set(minimum, 1, 0);
  printf("The initial value was (x,y) = (0,0)\n");
  //printf("and actual the minimum is at (x,y) = (1,1)\n");
  quasiNewtonMethod(RosenbrockValleyFunc, minimum, tolerance);
  printf("The found minimum is (%g,%g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
  printf("____________________________________________________\n");


  printf("Testing minimization routine on the Himmelblau function\n");
  gsl_vector_set(minimum, 0, 2.8);
  gsl_vector_set(minimum, 1, 1.8);
  printf("The initial value was (x,y) = (1,1)\n");
  printf("and actual the minima are at (x,y) = (3,2), (-2.80,3.13), (-3.78,-3.28), (3.58,-1.84)\n");
  quasiNewtonMethod(HimmelblauFunc, minimum, tolerance);
  printf("The found minimum is (%g,%g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));



  printf("--------------------- PART B ---------------------\n");
  printf("Fitting to the Higgs particle discovery data\n\n");

  numOfDataPts = 30;
  energy       = malloc(numOfDataPts * sizeof(double));
  crossSection = malloc(numOfDataPts * sizeof(double));
  error        = malloc(numOfDataPts * sizeof(double));

  double massMeasured = 125.3;

  input2Array(energy, crossSection, error, "higgsData.txt");

  dims = 3;
  gsl_vector* minimumHiggs = gsl_vector_alloc(dims);
  gsl_vector_set(minimumHiggs, 0, massMeasured+1.0);
  gsl_vector_set(minimumHiggs, 1, 2.5);
  gsl_vector_set(minimumHiggs, 2, 9.0);

  printf("The initial value was (m,Γ,A) = (%g, %g, %g)\n", massMeasured + 1.0, 2.5, 9.0);

  quasiNewtonMethod(deviationFunc, minimumHiggs, tolerance);

  printf("The found minimum is (m,Γ,A) = (%g, %g, %g)\n", gsl_vector_get(minimumHiggs, 0), gsl_vector_get(minimumHiggs, 1), gsl_vector_get(minimumHiggs, 2));

  FILE* outputStream = fopen("higgsFit.txt","w");
  for (int i = 0; i < numOfDataPts; i++) {
    double funcVals = BreitWignerFunc(gsl_vector_get(minimumHiggs, 0), gsl_vector_get(minimumHiggs, 1), gsl_vector_get(minimumHiggs, 2), energy[i]);
    fprintf(outputStream, "%g\t%g\t%g\t%g\n", energy[i], crossSection[i], error[i], funcVals);
  }
  fclose(outputStream);


  printf("--------------------- PART C ---------------------\n");
  printf("Implementing the downhill simplex method");
  printf("\nand testing on 1+(x-a)^2+(y-b)^2\n\n");


  dims                    = 2;
  int numOfPts            = dims + 1;
  double simplex_sizeGoal = tolerance;

  double** simplex = malloc(numOfPts*sizeof(double));
  for (int i = 0; i < numOfPts; i++) {
    simplex[i] = malloc(dims*sizeof(double));
  }

  simplex[0][0] = 4;
  simplex[0][1] = 10;
  simplex[1][0] = 2;
  simplex[1][1] = 2;
  simplex[2][0] = -3;
  simplex[2][1] = 1;
  printf("The initial values are { (x_0, y_0), (x_1, y_1), (x_2, y_2) } = { (%g, %g), (%g, %g), (%g, %g) }\n", simplex[0][0], simplex[0][1], simplex[1][0], simplex[1][1], simplex[2][0], simplex[2][1]);


  int numOfSteps = downhillsimplex(simplexTestFunc, simplex , dims , simplex_sizeGoal);

  printf("Number of steps = %d\n", numOfSteps);
  printf("The actual minimum is at: (x, y) = (%d, %d)\n", 6, 13);
  printf("and the found minimum, is at:  (x, y) = (%g, %g)\n", simplex[0][0], simplex[0][1]);

return 0;

}
