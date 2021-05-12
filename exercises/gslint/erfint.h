#ifndef HAVE_ERFINT_H
#define HAVE_ERFINT_H

double erfIntegral( double x, void* params );
double integratedErrorFunction( double upperLimit, const gsl_function* myGSLFUNC, double epsabs, double epsrel, size_t limit );

#endif
