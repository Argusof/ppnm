#include <float.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "rootfinding.h"
#include "GSdecomp.h"

void newtonRaphsonMethod(void func(gsl_vector* vals, gsl_vector* funcVals), gsl_vector* startPt, double tolerance){
    int dim = startPt -> size;
    double stepSize  = sqrt(DBL_EPSILON);
    int count = 0;

    gsl_matrix* jacobianMatrix  = gsl_matrix_alloc(dim, dim);
    gsl_vector* funcVals        = gsl_vector_alloc(dim);
    gsl_vector* funcVals_tmp    = gsl_vector_alloc(dim);
    gsl_vector* funcValsDiff    = gsl_vector_alloc(dim);
    gsl_matrix* triangMat       = gsl_matrix_alloc(dim, dim);
    gsl_vector* solution        = gsl_vector_alloc(dim);
    gsl_vector* solution_scaled = gsl_vector_alloc(dim);
    gsl_vector* nextPt          = gsl_vector_alloc(dim);
    gsl_vector* nextFuncVal     = gsl_vector_alloc(dim);

    func(startPt, funcVals);
    while (gsl_blas_dnrm2(funcVals) > tolerance) {
    count++;
    assert(count < 1e5);

        for (int valId = 0; valId < dim; valId++){
            gsl_vector_set(startPt, valId, gsl_vector_get(startPt, valId) + stepSize);
            func(startPt, funcVals_tmp);

            for (int id = 0; id < dim; id++){
                double funcVals_tmp_id = gsl_vector_get(funcVals_tmp, id);
                double funcVals_id     = gsl_vector_get(funcVals, id);
                double funcValDiff     = funcVals_tmp_id - funcVals_id;
                double matrixElement   = funcValDiff / stepSize;
                gsl_matrix_set(jacobianMatrix, id, valId, matrixElement);
            }

            gsl_vector_set(startPt, valId, gsl_vector_get(startPt, valId) - stepSize);
        }
        gsl_vector_scale(funcVals, -1.0);
        GS_decomp( jacobianMatrix, triangMat );
        GS_solve( jacobianMatrix, triangMat, funcVals, solution);
        gsl_vector_scale(funcVals, -1.0);

        double scale = 2;
        while ((gsl_blas_dnrm2(nextFuncVal) >= (1 - scale / 2) * gsl_blas_dnrm2(funcVals) ) && scale >= 0.02){
            scale /= 2;
            gsl_vector_memcpy(solution_scaled, solution);
            gsl_vector_scale(solution_scaled, scale);
            gsl_vector_memcpy(nextPt, solution_scaled);
            gsl_vector_add(nextPt, startPt);

            func(nextPt, nextFuncVal);

        }
        gsl_vector_memcpy(startPt, nextPt);
        gsl_vector_memcpy(funcVals, nextFuncVal);

        if ( (gsl_blas_dnrm2(solution) < stepSize) ) {
            break;
        }
    }

    gsl_vector_free(funcVals);
    gsl_vector_free(nextFuncVal);
    gsl_vector_free(funcVals_tmp);
    gsl_vector_free(funcValsDiff);
    gsl_vector_free(solution);
    gsl_vector_free(solution_scaled);
    gsl_vector_free(nextPt);
    gsl_matrix_free(triangMat);
    gsl_matrix_free(jacobianMatrix);
}
