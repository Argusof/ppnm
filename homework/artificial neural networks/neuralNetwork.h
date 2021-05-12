#ifndef HAVE_NEURALNETWORK_H
#define HAVE_NEURALNETWORK_H

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>


typedef struct { int numOfNeurons; double(*targetFunc)(double); gsl_vector* params; } neuralNetwork;

neuralNetwork*   neuralNetwork_alloc    (int numOfNeurons, double(*targetFunc)(double));
void             neuralNetwork_free     (neuralNetwork* network);
double           neuralNetwork_response (neuralNetwork* network, double evalPt);
void             neuralNetwork_train    (neuralNetwork* network, gsl_vector* inputVals, gsl_vector* labels);

#endif
