#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char** argv) {
  if(argc<2) fprintf(stderr,"%s: there were no arguments\n",argv[0]);
	else {
		for(int i=1;i<argc;i++){
			double x = atof(argv[i]);
			fprintf(stdout,"argument number %i = %g\tsin(%g)=%g\tcos(%g)=%g\n",i,x,x,sin(x),x,cos(x));
		}
	}
  return 0;
}
