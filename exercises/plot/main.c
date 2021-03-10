#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include "my_erf.h"
#include "my_gamma.h"

int main(){

	double xmin_erf = -2; // lower bound
  double xmax_erf = 2;  // upper bound

	FILE * output_erf;
	FILE * output_gamma;
  output_erf = fopen ("output_erf.txt", "w");
	output_gamma = fopen ("output_gamma.txt", "w");

	if (output_erf == NULL){
		printf("Could not open erf file\n");
	}
	else {
  	for(double x = xmin_erf; x <= xmax_erf ; x += 1.0/8){
			fprintf(output_erf,"%10g %10g %10g %10g\n",x,erf(x), gsl_sf_erf(x),myerf(x)); // print x and all three errorfuncs
		}
	}


	double xmin_gamma = -4; // lower bound
  double xmax_gamma = 4;  // upper bound

	if (output_gamma == NULL){
		printf("Could not open gamma file\n");
	}
	else {
  	for(double x = xmin_gamma; x <= xmax_gamma ; x += 1.0/20){
			if ( fabs((int)x - x) > 0.01){
				fprintf(output_gamma,"%10g %10g %10g %10g\n",x,tgamma(x), gsl_sf_gamma(x),mygamma(x)); // print x and all three errorfuncs
			}
		}
	}

	fclose(output_erf);
	fclose(output_gamma);

return 0;
}
