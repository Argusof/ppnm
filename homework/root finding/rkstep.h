#ifndef HAVE_RKSTEP_H
#define HAVE_RKSTEP_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void SchrodingerFunc(double var, gsl_vector* funcVals, gsl_vector* funcDeriv);

void rkstep12(
	void (*func)( double, gsl_vector* , gsl_vector* ), /* the f from dy/dt=f(t,y) */
	double var,                                        /* the current value of the variable */
	gsl_vector* funcVals,                              /* the current value y(t) of the sought function */
	double step,                                       /* the step to be taken */
	gsl_vector* funcStep,                              /* output: y(t+h) */
	gsl_vector* err                                    /* output: error estimate */
);

void rkstep45(
	void (*func)( double, gsl_vector* , gsl_vector* ), /* the f from dy/dt=f(t,y) */
	double var,                                        /* the current value of the variable */
	gsl_vector* funcVals,                              /* the current value y(t) of the sought function */
	double step,                                       /* the step to be taken */
	gsl_vector* funcStep,                              /* output: y(t+h) */
	gsl_vector* err                                    /* output: error estimate */
);

int rkdriver(
    void (*func)(double , gsl_vector* , gsl_vector* ), /* right-hand-side of dy/dt=f(t,y) */
    double leftEndpt,                                  /* the start-point a */
    gsl_vector* funcValLeft,                           /* y(a) */
    double rightEndpt,                                 /* the end-point of the integration */
    gsl_vector* funcValright,                          /* y(b) to be calculated */
    double step,                                       /* initial step-size */
    double absAcc,                                     /* absolute absAccuracy goal */
    double relAcc,                                     /* relative absAccuracy goal */
    FILE* filePath																		 /* path to solution file */
	);



#endif
