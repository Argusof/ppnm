#ifndef HAVE_ROOTFINDING_H
#define HAVE_ROOTFINDING_H

#include <gsl/gsl_vector.h>

void newtonRaphsonMethod(void func(gsl_vector* vals, gsl_vector* funcVals), gsl_vector* startPt, double tolerance);


#endif
