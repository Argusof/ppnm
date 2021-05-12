#ifndef HAVE_BESSELINT_H
#define HAVE_BESSELINT_H

typedef struct besselParams { int n; double x;} besselParams;

double besselIntegral( double tau, void* params );
double integratedBesselFunction( const gsl_function* myGSLFUNC, double epsabs, double epsrel, size_t limit );

#endif
