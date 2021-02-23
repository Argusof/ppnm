#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char** argv) {
  int input = 1;
  int inputArgs;

  while(input != EOF) {                          // end of input stream
    input = fscanf(stdin, "%i", &inputArgs);     // send input numbers into inputArgs
    if (input) {                                 // input = 1 means succesfull reading 
      fprintf(stdout, "input = %i\tsin(%i)=%g\tcos(%i)=%g\n", inputArgs, inputArgs, sin(inputArgs), inputArgs, cos(inputArgs));
    }
    else {
		    fprintf(stderr, "%s: Could not read input file\n", argv[0]);
    }
  }
  return 0;
}
