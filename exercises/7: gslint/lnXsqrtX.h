#ifndef HAVE_LNXSQRTX_H
#define HAVE_LNXSQRTX_H

double logxSqrtx( double x, void* params );
double integratedFunction( double lowerLimit, double upperLimit, const gsl_function* myGSLFUNC, double epsabs, double epsrel, size_t limit );

#endif
