#ifndef HAVE_NEURALNETWORK_H
#define HAVE_NEURALNETWORK_H

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>


typedef struct { int numOfNeurons; double(*targetFunc)(double); double (*targetDeriv) (double); double (*targetInt) (double);gsl_vector* params; } neuralNetwork;

neuralNetwork*   neuralNetwork_alloc    (int numOfNeurons, double(*targetFunc_input)(double), double(*targetDeriv_input)(double), double(*targetInt_input)(double));
void             neuralNetwork_free     (neuralNetwork* network);
double           neuralNetwork_response (neuralNetwork* network, double evalPt);
double           neuralNetwork_deriv    (neuralNetwork* network, double evalPt);
double           neuralNetwork_int      (neuralNetwork* network, double rightPt, double leftPt);
void             neuralNetwork_train    (neuralNetwork* network, gsl_vector* inputVals, gsl_vector* labels);

#endif
