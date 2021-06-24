#include <gsl/gsl_integration.h>

double integratedFunction( double lowerLimit, double upperLimit, const gsl_function* myGSLFUNC, double epsabs, double epsrel, size_t limit ){
    /* Function to integrate the GSL_FUNCTION on (0, 1)
     *
     *  ¤
     */

    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc( limit );
    double result;
    double absError;

    gsl_integration_qags(myGSLFUNC, lowerLimit, upperLimit, epsabs, epsrel, limit, workspace, &result, &absError);

    gsl_integration_workspace_free(workspace);
    workspace = NULL;

    return result;
}
