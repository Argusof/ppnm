#include <stdio.h>
#include <stdlib.h>

void fillArrayFromInputData( double* dataX,  double* dataY, char* inputFileName) {

  int numOfReturnVals = 3;
  int input = numOfReturnVals;
  double inputArgX;
  double inputArgY;
  double inputArgYdev;

  FILE* inputStream = fopen(inputFileName,"r"); // Open file input.txt for reading

  int id = 0;
  while(input != EOF) {                               // end of input stream
    if (input == numOfReturnVals) {                                 // input = 1 means succesfull reading
      input = fscanf(inputStream, "%lg\t%lg\t%lg", &inputArgX, &inputArgY, &inputArgYdev);  // send input numbers from file into inputArgs

      dataX[id] = inputArgX;
      dataY[id] = inputArgY;

      id++;

    }
    else {
		    fprintf(stderr, "Could not read input file\n");
        exit(-1);
    }
  }
  fclose(inputStream);

}
