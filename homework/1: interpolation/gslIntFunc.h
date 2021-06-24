#ifndef HAVE_GSLINTFUNC_H
#define HAVE_GSLINTFUNC_H

#include <gsl/gsl_integration.h>

double integratedFunction( double lowerLimit, double upperLimit, const gsl_function* myGSLFUNC, double epsabs, double epsrel, size_t limit );

#endif
