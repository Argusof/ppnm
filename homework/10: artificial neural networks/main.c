#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "minimization.h"
#include "neuralNetwork.h"
#include "input2Array.h"

double nsin( double val ){
    return -sin(val);
}

int main(int argc, char *argv[]) {

  printf("--------------------- PART A ---------------------\n");
  printf("Neural network implemented and tested on cos(x) - see plot.png\n");


  int numOfDataPts = 20;
  int numOfNeurons = 5;
  gsl_vector* XData = gsl_vector_alloc(numOfDataPts);
  gsl_vector* YData = gsl_vector_alloc(numOfDataPts);
  input2Array(numOfDataPts, XData, YData, argv[1]);

  neuralNetwork* network = neuralNetwork_alloc(numOfNeurons, &cos, &nsin, &sin);
  neuralNetwork_train(network, XData, YData);

  double xMin = 0;
  double xMax = 11;
  FILE* outputStream = fopen(argv[2], "w");
  for( double val = xMin; val <= xMax; val += 1.0/8 ){
    fprintf(outputStream, "%10g %10g %10g %10g %10g %10g %10g\n", val, neuralNetwork_response(network, val), cos(val), neuralNetwork_deriv(network, val), -sin(val), neuralNetwork_int(network, 0, val), sin(val));
  }

  neuralNetwork_free(network);

return 0;

}
