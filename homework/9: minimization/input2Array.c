#include <stdio.h>
#include <stdlib.h>
#include "input2Array.h"

void input2Array(double* XData, double* YData, double* ZData, char* inputFilename ){
  int numOfReturnVals = 3;  // return value of fscanf() is the number of fetched variables

  int    input = numOfReturnVals; 			// will receive the return value from fscanf
  double argX;						              // X data from stdin
  double argY;						              // X data from stdin
  double argZ;						              // Z data from stdin

  FILE* myInputFileStream  = fopen(inputFilename, "r");

  int id = 0;
  while( input != EOF) { // While we are not at the end of the input stream
    if ( input == numOfReturnVals){		 // If input has returned success
      input = fscanf(myInputFileStream, "%lg\t%lg\t%lg", &argX, &argY, &argZ);

      XData[id] = argX;
      YData[id] = argY;
      ZData[id] = argZ;

      id++;
    }
    else { // If we are not successfull, say if there was an error with fscanf()
      fprintf(stderr, "Failed to read input.\n"); // Print this to stderr so we know
      exit(-1); 																	// Terminate program
    }
  }

  fclose(myInputFileStream);
}
