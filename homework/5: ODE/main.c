#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include "rkstep.h"


int main(int argc, char* argv[]) {
  if (argc < 2 ){
    fprintf(stderr, "No arguments passed");
    exit(1);
  }

  //_________________________HARMONIC DIFFERENTIAL EQUATION_________________________
  int orderHarm     = 2;
  double leftEndpt  = 0.0;
  double rightEndpt = 2 * M_PI;
  double step       = (rightEndpt - leftEndpt)/10;
  double absAcc     = 1e-3;
  double relAcc     = 1e-3;

  gsl_vector* HarmFuncValLeft  = gsl_vector_calloc(orderHarm);
  gsl_vector* HarmFuncValRight = gsl_vector_alloc(orderHarm);
  gsl_vector_set(HarmFuncValLeft, 1, 1);

  FILE* harmOutputStream = fopen(argv[1], "w");
  rkdriver(&harmFunc, leftEndpt, HarmFuncValLeft, rightEndpt, HarmFuncValRight, step, absAcc, relAcc, harmOutputStream);
  fclose(harmOutputStream);

  gsl_vector_free(HarmFuncValLeft);
  gsl_vector_free(HarmFuncValRight);


  //_________________________SIR MODEL_________________________

  int orderSIR = 3;
  leftEndpt    = 0.0;
  rightEndpt   = 100;
  double population   = 5808180;
  double wasInfected  = 237792;
  double dead         = 2441;
  double recovered    = 226630;
  double vaccinated   = 445566;
  double isInfected   = wasInfected - recovered;
  double removed      = dead + recovered + vaccinated;

  gsl_vector* SIRfuncValLeft  = gsl_vector_calloc(orderSIR);
  gsl_vector* SIRfuncValRight = gsl_vector_alloc(orderSIR);
  gsl_vector_set(SIRfuncValLeft, 0, population - isInfected - removed);
  gsl_vector_set(SIRfuncValLeft, 1, isInfected);
  gsl_vector_set(SIRfuncValLeft, 2, removed);

  FILE* SIRoutputStream = fopen(argv[2], "w");
  rkdriver(&SIR, leftEndpt, SIRfuncValLeft, rightEndpt, SIRfuncValRight, step, absAcc, relAcc, SIRoutputStream);
  fclose(SIRoutputStream);

  gsl_vector* SIRfuncValLeft2  = gsl_vector_calloc(orderSIR);
  gsl_vector* SIRfuncValRight2 = gsl_vector_alloc(orderSIR);
  gsl_vector_set(SIRfuncValLeft2, 0, population - isInfected - removed);
  gsl_vector_set(SIRfuncValLeft2, 1, isInfected);
  gsl_vector_set(SIRfuncValLeft2, 2, removed);

  FILE* SIRoutputStream2 = fopen(argv[3], "w");
  rkdriver(&SIR2, leftEndpt, SIRfuncValLeft2, rightEndpt, SIRfuncValRight2, step, absAcc, relAcc, SIRoutputStream2);
  fclose(SIRoutputStream2);

  gsl_vector_free(SIRfuncValLeft);
  gsl_vector_free(SIRfuncValRight);
  gsl_vector_free(SIRfuncValLeft2);
  gsl_vector_free(SIRfuncValRight2);

//_________________________THREE BODY PROBLEM_________________________

  int order3body = 12;
  leftEndpt    = 0.0;
  rightEndpt   = 6.32591398;

  gsl_vector* threebodyFuncValRight  =  gsl_vector_alloc(order3body);
  gsl_vector* threebodyFuncValLeft   =  gsl_vector_calloc(order3body);

    double init_pos_x_1 =  0.97000436;
      double init_pos_y_1 = -0.24308753;
      double init_pos_x_2 = -0.97000436;
      double init_pos_y_2 =  0.24308753;
      double init_pos_x_3 =  0 ;
      double init_pos_y_3 =  0 ;
      double init_vel_x_3 = -0.93240737;
      double init_vel_y_3 = -0.86473146;
      double init_vel_x_1 = -init_vel_x_3/2;
      double init_vel_y_1 = -init_vel_y_3/2;
      double init_vel_x_2 = -init_vel_x_3/2;
      double init_vel_y_2 = -init_vel_y_3/2;

      gsl_vector_set(threebodyFuncValLeft, 0, init_pos_x_1);
      gsl_vector_set(threebodyFuncValLeft, 1, init_pos_y_1);
      gsl_vector_set(threebodyFuncValLeft, 2, init_pos_x_2);
      gsl_vector_set(threebodyFuncValLeft, 3, init_pos_y_2);
      gsl_vector_set(threebodyFuncValLeft, 4, init_pos_x_3);
      gsl_vector_set(threebodyFuncValLeft, 5, init_pos_y_3);
      gsl_vector_set(threebodyFuncValLeft, 6, init_vel_x_1);
      gsl_vector_set(threebodyFuncValLeft, 7, init_vel_y_1);
      gsl_vector_set(threebodyFuncValLeft, 8, init_vel_x_2);
      gsl_vector_set(threebodyFuncValLeft, 9, init_vel_y_2);
      gsl_vector_set(threebodyFuncValLeft, 10, init_vel_x_3);
      gsl_vector_set(threebodyFuncValLeft, 11, init_vel_y_3);

  FILE* outputStream3body = fopen(argv[4], "w");
  rkdriver(&threeBody, leftEndpt, threebodyFuncValLeft, rightEndpt, threebodyFuncValRight, step, absAcc, relAcc, outputStream3body);
  fclose(outputStream3body);

  gsl_vector_free(threebodyFuncValLeft);
  gsl_vector_free(threebodyFuncValRight);

  return 0;
}
