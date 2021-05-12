#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include "rkstep.h"

/*
void SchrodingerFunc(double var, gsl_vector* funcVals, gsl_vector* funcDeriv){
	gsl_vector_set(funcDeriv, 0, gsl_vector_get(funcVals,1));
	gsl_vector_set(funcDeriv, 1, -2.0*(energy+1.0/var)*gsl_vector_get(funcVals, 0));
}
*/

void rkstep12(
	void (*func)( double, gsl_vector* , gsl_vector* ), /* the f from dy/dt=f(t,y) */
	double var,                                        /* the current value of the variable */
	gsl_vector* funcVals,                              /* the current value y(t) of the sought function */
	double step,                                       /* the step to be taken */
	gsl_vector* funcStep,                              /* output: y(t+h) */
	gsl_vector* err                                    /* output: error estimate */
) {

	int order = funcVals->size;
	gsl_vector* k0          = gsl_vector_alloc(order);
	gsl_vector* k12         = gsl_vector_alloc(order);
	gsl_vector* tmpFuncVal  = gsl_vector_alloc(order);

	func(var, funcVals, k0);  		// compute k0
	for (int id = 0; id < order; id++) {
		double k0_id = gsl_vector_get(k0, id);
		double funcVals_id = gsl_vector_get(funcVals, id);
		double tmpFVal = funcVals_id + .5 * step * k0_id;
		gsl_vector_set(tmpFuncVal, id, tmpFVal);
	}

	func(var + .5 * step, tmpFuncVal, k12);
	for (int id = 0; id < order; id++) {
		double k12_id = gsl_vector_get(k12, id);
		double funcVals_id = gsl_vector_get(funcVals, id);
		double tmpFuncStep = funcVals_id + step * k12_id;
		gsl_vector_set(funcStep, id, tmpFuncStep);
	}

	for (int id = 0; id < order; id++) {
		double k0_id = gsl_vector_get(k0, id);
		double k12_id = gsl_vector_get(k12, id);
		double tmpErr_id = (k0_id - k12_id) * (step / 2);
		gsl_vector_set(err, id, tmpErr_id);
	}

	gsl_vector_free(k0);
	gsl_vector_free(k12);
}

void rkstep45(
	void (*func)( double, gsl_vector* , gsl_vector* ), /* the f from dy/dt=f(t,y) */
	double var,                                        /* the current value of the variable */
	gsl_vector* funcVals,                              /* the current value y(t) of the sought function */
	double step,                                       /* the step to be taken */
	gsl_vector* funcStep,                              /* output: y(t+h) */
	gsl_vector* err                                    /* output: error estimate */
) {

    int order = funcVals->size;
    gsl_vector* funcDeriv = gsl_vector_alloc(order);
    gsl_vector* k0        = gsl_vector_alloc(order);
    gsl_vector* k1        = gsl_vector_alloc(order);
    gsl_vector* k2        = gsl_vector_alloc(order);
    gsl_vector* k3        = gsl_vector_alloc(order);
    gsl_vector* tangent   = gsl_vector_alloc(order);
    gsl_matrix* identity  = gsl_matrix_alloc(order, order);


    gsl_matrix_set_identity(identity);
    func(var, funcVals, funcDeriv);

    gsl_vector_memcpy(k0, funcDeriv);

    gsl_vector_memcpy(k1,k0);
    gsl_blas_dgemv(CblasNoTrans, 1.0, identity, funcVals, 0.5*step, k1);
    func(var + 0.5*step, k1, k1 );

    gsl_vector_memcpy(k2,k1);
    gsl_blas_dgemv(CblasNoTrans, 1.0, identity, funcVals, 0.5*step, k2);
    func(var + 0.5*step, k2, k2);

    gsl_vector_memcpy(k3,k2);
    gsl_blas_dgemv(CblasNoTrans, 1.0, identity, funcVals, step, k3);
    func(var + step, k3, k3);

     gsl_blas_dgemv(CblasNoTrans, 1.0/6.0, identity, k0, 1.0/3.0, k1);
     gsl_blas_dgemv(CblasNoTrans, 1.0/3.0, identity, k2, 1.0/6.0, k3);
     gsl_blas_dgemv(CblasNoTrans, 1.0, identity, k1, 1.0, k3);
     gsl_vector_memcpy(tangent, k3);

     gsl_vector_free(k0);
     gsl_vector_free(k1);
     gsl_vector_free(k2);
     gsl_vector_free(k3);

     gsl_blas_dgemv(CblasNoTrans, 1.0, identity, funcVals, step, tangent);
     gsl_vector_memcpy(funcStep, tangent);

     gsl_vector_memcpy(err, funcStep);
     gsl_blas_dgemv(CblasNoTrans, 1.0, identity, funcVals, -1.0, err);

     gsl_vector_free(tangent);
     gsl_vector_free(funcDeriv);
     gsl_matrix_free(identity);
}


int rkdriver(
    void (*func)(double , gsl_vector* , gsl_vector* ), /* right-hand-side of dy/dt=f(t,y) */
    double leftEndpt,                                  /* the start-point a */
    gsl_vector* funcValLeft,                           /* y(a) */
    double rightEndpt,                                 /* the end-point of the integration */
    gsl_vector* funcValRight,                          /* y(b) to be calculated */
    double step,                                       /* initial step-size */
    double absAcc,                                     /* absolute absAccuracy goal */
    double relAcc,                                     /* relative absAccuracy goal */
    FILE* filePath																		 /* path to solution file */
	  ) {

    int order = funcValLeft -> size;
    double err;
    double normFuncVal;
    double tol;

    gsl_vector* funcStep    = gsl_vector_alloc(order);   // func val at next step
    gsl_vector* funcErr     = gsl_vector_alloc(order);	 // func val error (dy)
		gsl_vector* thisFuncVal = gsl_vector_alloc(order);	 // func val at current step
		gsl_vector_memcpy(thisFuncVal, funcValLeft); 				 // initialize y (vector) as left endpoint

		double pos = leftEndpt;		// current position
    while(pos < rightEndpt){	// while not at the end of the interval
			if(filePath != NULL){		// write data to file
				fprintf(filePath, "%.5g\t", pos);
				for (int id = 0; id < order; id++) {
					fprintf(filePath, "%.5g\t",gsl_vector_get(thisFuncVal, id));
				}
				fprintf(filePath, "%.5g\n",pos*exp(-pos));
			}

			double trueStep; 				// final stepsize
			double nextStep = step; // start by initial stepsize

    	if(pos + nextStep > rightEndpt){  // if not inside interval
          nextStep = rightEndpt - pos;
      }

			do{
				rkstep12(func, pos, thisFuncVal, nextStep, funcStep, funcErr); // do step

				err 				= gsl_blas_dnrm2(funcErr);     // compute error
				normFuncVal = gsl_blas_dnrm2(funcStep);		 // compure r2 norm

	      tol = (normFuncVal * relAcc + absAcc) * sqrt(nextStep / (rightEndpt - leftEndpt));  // compute tolerance
				trueStep  = nextStep;
				nextStep *= pow(tol/err, .25) * .95;
			} while( err > tol ); 											 // make stepsize smaller until within tolerance

			gsl_vector_memcpy(thisFuncVal, funcStep);
			pos += trueStep;
		}

		gsl_vector_memcpy(funcValRight, funcStep); 		 // last val is the right endpoint

		gsl_vector_free(funcStep);
		gsl_vector_free(funcErr);
		gsl_vector_free(thisFuncVal);
}
