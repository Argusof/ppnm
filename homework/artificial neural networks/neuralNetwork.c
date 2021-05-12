#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <time.h>
#include "minimization.h"
#include "utilities.h"
#include "neuralNetwork.h"


neuralNetwork* neuralNetwork_alloc(int numOfNeurons_input, double(*targetFunc_input)(double)){
  neuralNetwork* network = (neuralNetwork*)malloc(sizeof(neuralNetwork));
  network->params       = gsl_vector_alloc(numOfNeurons_input*3);
  network->targetFunc   = targetFunc_input;
  network->numOfNeurons = numOfNeurons_input;
  return network;
}

void neuralNetwork_free(neuralNetwork* network){
  gsl_vector_free(network->params);
  free(network);
}

double neuralNetwork_response(neuralNetwork* network, double evalPt){
  int numOfNeurons = network->numOfNeurons;
  int numOfParams = 3;
  double response = 0;

  for (int neuronId = 0; neuronId < numOfNeurons; neuronId++) {
    double shift      = gsl_vector_get(network->params, neuronId*numOfParams + 0);
    double scale      = gsl_vector_get(network->params, neuronId*numOfParams + 1);
    double edgeWeight = gsl_vector_get(network->params, neuronId*numOfParams + 2);

    response += (network->targetFunc((evalPt - shift) / scale)) * edgeWeight;
  }
  return response;
}

void neuralNetwork_train(neuralNetwork* network, gsl_vector* inputVals, gsl_vector* labels){
  unsigned int seed = time(NULL);
  int numOfParams = 3;
  int numOfPts = inputVals->size;
  int numOfNeurons = network->numOfNeurons;

  double costFunc(gsl_vector* nextParams){
    int numOfNeurons = network->numOfNeurons;
    gsl_vector* updatedParams = gsl_vector_alloc(numOfNeurons*numOfParams);

    for (int neuronId = 0; neuronId < numOfNeurons; neuronId++) {
      for (int paramId = 0; paramId < numOfParams; paramId++) {
        gsl_vector_set(updatedParams, neuronId*numOfParams + paramId, gsl_vector_get(nextParams, neuronId*numOfParams + paramId));
      }
    }
    network->params = updatedParams;

    double cost = 0;

    for (int pt = 0; pt < numOfPts; pt++){
      double evalPt   = gsl_vector_get(inputVals, pt);
      double labelPt  = gsl_vector_get(labels, pt);
      double response = neuralNetwork_response(network, evalPt);

      cost += (response - labelPt)*(response - labelPt);
    }
    return cost;
  }

  double tolerance = 1e-5;
  gsl_vector* minimizedParams = gsl_vector_alloc(numOfParams*numOfNeurons);

  for (int i = 0; i < numOfNeurons*numOfParams; i++) {
    gsl_vector_set(minimizedParams, i, randomNumber(&seed));
  }


  quasiNewtonMethod(costFunc, minimizedParams, tolerance);

  network->params = minimizedParams;
}
