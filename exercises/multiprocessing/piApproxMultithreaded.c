#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

#include "rndNum.h"
#include "piApproxMultithreaded.h"

#define NUM_OF_THREADS 2

typedef struct placePtsArgStruct {      //  Begin by defining a struct that may be
    unsigned int  seed              ;   //  passed to the function placePts below.
    double        numOfPts          ;   //  The struct contains fields representing
    int*          ptsInCircle       ;   //  the RNG state seed, number of points to
} placePtsArgStruct;                    //  use in the simulation, and a pointer to
                                        //  keep track of points within the circle.

void* placePts( void* placePtsArgStructInput ){
  /*
     Place points randomly within a square with sides of length 2.
     count how many points fall inside the unit circle enclosed
     by this square.

     ¤ placePtsArgStructInput : Struct of fields;
                                    ¤ unsigned int  seed : Seed to set RNG state
                                    ¤ double numOfPts    : number of points to simulate with
                                    ¤ int* ptsInCircle   : pointer to integer valued variable
                                                           that keeps count of points falling
                                                           within the circle area
  */
  placePtsArgStruct* args   =   ( placePtsArgStruct* )placePtsArgStructInput;

  unsigned int seed         =   args->seed              ;  //  Fill local variables from struct
  double numOfPts           =   args->numOfPts          ;
  int* ptsInCircle          =   args->ptsInCircle       ;

  double unitCircleRadius = 1.0;  //  simply 1.0, the radius of the unit circle enclosed by the square

  double xpos   =   0;  //  Declare variables to temporarily hold RNG values in x- and y- coordinates
  double ypos   =   0;

  for ( int iteration = 0; iteration < numOfPts; iteration++ ){
    xpos = randomNumber(&seed);   //  Generate random points, uniformly distributed on (0, 1) x (0, 1).
    ypos = randomNumber(&seed);

    if ( sqrt( pow(xpos, 2) + pow(ypos, 2) ) <= unitCircleRadius ){   //  If the points are within unit circle
      (*ptsInCircle)++;                                               //  then we update the count at the address
    }
  }

  pthread_exit( (void*) placePtsArgStructInput );   //  We should make sure to exit the thread after it has done its job
  return NULL;
}


void piApproxMultithreaded( double numOfPts ){
  pthread_t threadArray[ NUM_OF_THREADS ];  //  Begin by declaring an array to hold the threads, NUM_OF_THREADS is defined as a preprocessor variable above
  pthread_attr_t attributes;                //  Declare the attribute state.

  placePtsArgStruct* placePtsArgStructArray = malloc( NUM_OF_THREADS*sizeof(placePtsArgStruct) );   //  Declare an array of structures of the form as defined above.

  int   error  ;    //  Error will hold the return value of pthread_ join/create () functions, which is 1/TRUE in case of failure, or 0/FALSE otherwise.
  void* status ;    //  Status is to hold the status state of a thread, such that we may join them correctly below.

  // Initiallize and set thread detached attribute
  pthread_attr_init(&attributes);                                      //  Initiallize attribute state.
  pthread_attr_setdetachstate(&attributes, PTHREAD_CREATE_JOINABLE);   //  It is important that we set the thread state to be JOINABLE and not DETACHED!

  unsigned int masterSeed   =   time(NULL) ;                           //  We initiallize a master seed from the system clock time
  int* ptsInCircleArray     =   malloc( NUM_OF_THREADS*sizeof(int) );  //  Declare memory dynamically for an array of ptsInCircle variables to hold counts in circle for each thread

  for (int threadId = 0; threadId < NUM_OF_THREADS; threadId++){       //  Loops to create and run threads in parallel.
    unsigned int threadSeed = masterSeed + threadId;                   //  Set a unique seed for each thread to use rand_r() in a thread safe manner
    ptsInCircleArray[ threadId ] = 0;                                  //  Initiallize the ptsInCircleArray variable for each thread, we pass a pointer to it in the struct array

    placePtsArgStructArray[ threadId ].seed               =  threadSeed                   ;   //  We fill out the struct array fields for each struct corresponding to a thread
    placePtsArgStructArray[ threadId ].numOfPts           =  numOfPts                     ;
    placePtsArgStructArray[ threadId ].ptsInCircle        =  &ptsInCircleArray[threadId]  ;

    printf("main(): creating thread #%d.\n", threadId + 1);

    error = pthread_create(&threadArray[ threadId ], &attributes, placePts, (void*) &(placePtsArgStructArray[threadId]));   //  Create each thread in turn, let it executer placePts with its own struct as argument
    if (error){                                                                                                             //  We should be careful to check for the success of pthread_create()
      printf("Error: pthread_create() returned failure while creating thread #%d.", threadId + 1 );                         //  Otherwise let us know in case of failure
      exit(-1);                                                                                                             // and terminate program
    }
  }
  printf("main(): threads created successfully. Destroying atrributes and joining.\n");

  pthread_attr_destroy(&attributes);                                                  // Destroy attribute state before joining.
  for (int threadId = 0; threadId < NUM_OF_THREADS; threadId++ ){
    error = pthread_join(threadArray[threadId], &status);                             // Join each thread
    if (error){                                                                       // Check for proper success just as before when we were creating
      printf("Error: pthread_join() returned failure in thread #%d.", threadId + 1);
      exit(-1);
    }
  }
  printf("main(): Joining successful!\n");

  double ptsInCircle = 0;                                           // Declare and initiallize a variable to hold the sum of points in circle, of the different threads
  for (int threadId = 0; threadId < NUM_OF_THREADS; threadId++ ){
    ptsInCircle   +=   ptsInCircleArray[threadId];                  // Compute sum
  }
  ptsInCircle /= NUM_OF_THREADS;                                    // We need to divide by number of threads to get actual value, since we ran numOfPts for each thread
  double myPiApprox    =   4.0*ptsInCircle/((double)numOfPts);

  printf("My approximation of pi using pthreads = %g\n"  , myPiApprox     );  // Now we can finish by writing the results to the terminal
  printf("Deviation is %g percent\n"      , fabs(1.0-myPiApprox/M_PI)*100 );

  free(ptsInCircleArray);   //  Make sure to free the dynamically allocated memory
}
