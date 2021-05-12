#ifndef HAVE_INPUT2ARRAY_H
#define HAVE_INPUT2ARRAY_H

#include <gsl/gsl_vector.h>

void input2Array(int numOfDataPts, gsl_vector* XData, gsl_vector* YData, char* inputFilename );

#endif
