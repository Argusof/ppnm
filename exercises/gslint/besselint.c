#include <gsl/gsl_integration.h>
#include <math.h>
#include "besselint.h"

double besselIntegral( double tau, void* params_input ){
    /* Defines a function that will be used to declare a GSL_FUNCTION below

      Inputs
       - tau    : query points
       - params : set of parameters, to obtain the necessary signature of the function
                  when using as a GSL_FUNCTION (NULL in this case)

       Returns function value.
     */
    besselParams* params = (besselParams*) params_input;
    double n = (params->n);
    double x = (params->x);

    double funcVal = (1/M_PI)*cos(n*tau - x*sin(tau));

    return funcVal;
}

double integratedBesselFunction( const gsl_function* myGSLFUNC, double epsabs, double epsrel, size_t limit ){
    // Function to integrate the GSL_FUNCTION on interval [0, pi]

    double lowerLimit = 0;
    double upperLimit = M_PI;

    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc( limit );
    double result;
    double absError;

    gsl_integration_qags(myGSLFUNC, lowerLimit, upperLimit, epsabs, epsrel, limit, workspace, &result, &absError);

    gsl_integration_workspace_free(workspace);
    workspace = NULL;

    return result;
}
