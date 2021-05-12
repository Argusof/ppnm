#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <float.h>
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.22045e-16
#endif
#include "minimization.h"




void numeric_gradient(double func(gsl_vector*), gsl_vector* minimum, gsl_vector* gradient){
    double stepSize = 2.22045e-10;

    double funcVal  =   func(minimum);
    int numOfDims   =   minimum -> size;

    for(int dimId = 0; dimId < numOfDims; ++dimId){
        double step;
        double minimum_i   =   gsl_vector_get(minimum, dimId);

        if (fabs(minimum_i) < stepSize) {
            step = stepSize;
        }
        else {
            step = fabs(minimum_i) * stepSize;
        }

        gsl_vector_set(minimum,  dimId,  minimum_i + step                 );
        gsl_vector_set(gradient, dimId,  (func(minimum) - funcVal) / step );
        gsl_vector_set(minimum,  dimId,  minimum_i - step                 );
    }
}

void quasiNewtonMethod(double func(gsl_vector*), gsl_vector* minimum, double tolerance){
    double stepSize = 2.22045e-10;
    int dimensions = minimum -> size;

    int numOfSteps      =   0;
    int numOfResets     =   0;
    int numOfScales     =   0;

    gsl_matrix* inverse_hessianMatrix   =   gsl_matrix_alloc(dimensions, dimensions);
    gsl_matrix_set_identity(inverse_hessianMatrix);

    gsl_matrix* identity   =   gsl_matrix_alloc(dimensions, dimensions);
    gsl_matrix_set_identity(identity);

    gsl_vector* gradientVal             =   gsl_vector_alloc(dimensions);
    gsl_vector* newtonStep              =   gsl_vector_alloc(dimensions);
    gsl_vector* minimum_next            =   gsl_vector_alloc(dimensions);
    gsl_vector* gradientVal_next        =   gsl_vector_alloc(dimensions);
    gsl_vector* solution                =   gsl_vector_alloc(dimensions);
    gsl_vector* solutionChange          =   gsl_vector_alloc(dimensions);
    gsl_vector* broydenVec              =   gsl_vector_alloc(dimensions);

    numeric_gradient(func, minimum, gradientVal);
    double funcVal  =  func(minimum);
    double funcVal_next;

    while(numOfSteps < 1e4){
        numOfSteps++;

        gsl_blas_dgemv(CblasNoTrans, -1, inverse_hessianMatrix, gradientVal, 0, newtonStep);
        if( gsl_blas_dnrm2(newtonStep) < stepSize * gsl_blas_dnrm2(minimum) ) {
            fprintf(stderr,"qnewton: |Dx|<stepSize*|x|\n");
            break;
        }
        if( gsl_blas_dnrm2(gradientVal) < tolerance ) {
            fprintf(stderr,"qnewton: |grad|<acc\n");
            break;
        }

        double scale = 1;

        while(1){
            gsl_vector_memcpy(minimum_next, minimum);
            gsl_vector_add(minimum_next, newtonStep);

            funcVal_next = func(minimum_next);

            double sTransGrad;
            gsl_blas_ddot(newtonStep, gradientVal, &sTransGrad);

            if(funcVal_next < funcVal + 0.01 * sTransGrad){
                numOfScales++;
                break;
            }
            if(scale < stepSize){
                numOfResets++;
                gsl_matrix_set_identity(inverse_hessianMatrix);
                break;
            }
            scale*=0.5;
            gsl_vector_scale(newtonStep, 0.5);
        }

        numeric_gradient(func, minimum_next, gradientVal_next); // Compute gradient in the next step

        gsl_vector_memcpy(solution, gradientVal_next);
        gsl_blas_daxpy(-1, gradientVal, solution); /* y=grad(x+s)-grad(x) */
        gsl_vector_memcpy(solutionChange, newtonStep); /* u=s */
        gsl_blas_dgemv(CblasNoTrans, -1, inverse_hessianMatrix, solution, 1, solutionChange); /* u=s-By */

        gsl_matrix* solChangeSolChangeTrans = gsl_matrix_calloc(dimensions, dimensions); //u*u^T
        gsl_blas_dsyr(CblasUpper, 1.0, solutionChange, solChangeSolChangeTrans);
        double solChangeTransSol; // u^T*y
        gsl_blas_ddot(solutionChange, solution, &solChangeTransSol);
        if(fabs(solChangeTransSol) > 1e-12){ // SR1 update
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0 / solChangeTransSol, solChangeSolChangeTrans, identity, 1.0, inverse_hessianMatrix); // B= B + delta B
        }

        gsl_vector_memcpy(minimum, minimum_next);
        gsl_vector_memcpy(gradientVal, gradientVal_next);
        funcVal = funcVal_next;
    }

    gsl_matrix_free(inverse_hessianMatrix);
    gsl_matrix_free(identity);
    gsl_vector_free(gradientVal);
    gsl_vector_free(newtonStep);
    gsl_vector_free(minimum_next);
    gsl_vector_free(gradientVal_next);
    gsl_vector_free(solution);
    gsl_vector_free(solutionChange);
    gsl_vector_free(broydenVec);

    printf("\nnumber of steps = %i\nnumber of scales = %i \nnumber of resets = %i \nfunction value = %.1e\n\n", numOfSteps, numOfScales, numOfResets, funcVal);

}
    void simplex_reflection(const double* highest_pt , const double* centroid_pt, int dimensions, double* reflected_pt ){
        for(int dimId = 0; dimId < dimensions ; ++dimId) {
            reflected_pt[dimId]    = 2 * centroid_pt[dimId] - highest_pt[dimId];
        }
    }

    void simplex_expansion(const double* highest_pt, const double* centroid_pt , int dimensions , double* expanded_pt ){
        for(int dimId = 0; dimId < dimensions ; ++dimId)  {
            expanded_pt[dimId]     = 3 * centroid_pt[dimId] - 2 * highest_pt[dimId];
        }
    }

    void simplex_contraction(const double* highest_pt, const double* centroid_pt, int dimensions , double* contracted_pt ){
        for(int dimId = 0; dimId < dimensions ; ++dimId) {
            contracted_pt[dimId]   = 0.5 * centroid_pt[dimId] + 0.5 * highest_pt[dimId];
        }
    }

    void simplex_reduction(double** simplex, int dimensions, int low_ptId ){
        for(int ptId = 0; ptId < dimensions + 1; ++ptId) {
            if(ptId != low_ptId ) {
                for (int dimId = 0; dimId < dimensions; ++dimId) {
                    simplex[ptId][dimId] = 0.5 * (simplex[ptId][dimId] + simplex[low_ptId][dimId]);
                }
            }
        }
    }

    double simplex_distance (double* first_pt , double* second_pt, int dimensions){
        double euclidean_distance = 0;
        for(int dimId = 0; dimId < dimensions ; dimId++)  {
            euclidean_distance += pow (second_pt[dimId] - first_pt[dimId], 2);
        }
        return sqrt(euclidean_distance);
    }
    double simplex_size(double** simplex, int dimensions ){
        double tmp_distance = 0;
        for(int ptId = 1; ptId < dimensions + 1; ++ptId){
            double dist = simplex_distance(simplex[0], simplex[ptId], dimensions) ;
            if(dist > tmp_distance){
                tmp_distance = dist;
            }
        }
        return tmp_distance;
    }

    void simplex_update (double** simplex, const double* funcVals, int dimensions , int* high_ptId , int* low_ptId , double* centroid_pt) {
        *high_ptId  =   0;
        *low_ptId   =   0;

        double highest_pt  =   funcVals[0];     // The function value of the point that has the highest function value
        double lowest_pt   =   funcVals[0];     // And the function value of the point that has the lowest function value

        // Find the highest and lowest points
        for (int ptId = 1; ptId < dimensions + 1; ++ptId) {
            double next_pt = funcVals[ptId];
            if (next_pt > highest_pt) {
                highest_pt = next_pt;
                *high_ptId = ptId;
            }
            if (next_pt < lowest_pt) {
                lowest_pt = next_pt;
                *low_ptId = ptId;
            }
        }

        // Find the centroid point
        for (int dimId = 0; dimId < dimensions; ++dimId) {
            double sum = 0;
            for (int ptId = 0; ptId < dimensions + 1; ++ptId) {
                if (ptId != *high_ptId) {
                    sum += simplex[ptId][dimId];
                }
            }
            centroid_pt[dimId] = sum / dimensions;
        }
    }

    void simplex_initiate (double func(double*), double** simplex, double* funcVals, int dimensions, int* high_ptId, int* low_ptId, double* centroid_pt ){
        for(int ptId = 0; ptId < dimensions + 1; ++ptId) {
            funcVals[ptId] = func(simplex[ptId]);
        }
        simplex_update(simplex, funcVals, dimensions, high_ptId, low_ptId, centroid_pt);
    }

    int downhillsimplex(double func(double*), double** simplex , int dimensions , double simplex_sizeGoal ){
        int high_ptId;
        int low_ptId;
        int numOfSteps = 0;

        double centroid[dimensions];
        double funcVals[dimensions + 1];
        double reflectedPt[dimensions];
        double expandedPt[dimensions];

        simplex_initiate(func, simplex, funcVals, dimensions, &high_ptId, &low_ptId, centroid) ;

        while(simplex_size(simplex, dimensions) > simplex_sizeGoal){
            simplex_update(simplex, funcVals, dimensions, &high_ptId, &low_ptId, centroid);

            // Try a reflection
            simplex_reflection(simplex[high_ptId], centroid, dimensions, reflectedPt) ;
            double funcVal_reflected = func(reflectedPt);

            if(funcVal_reflected < funcVals[low_ptId]){
                //  Reflection looks good, try expansion
                simplex_expansion(simplex[high_ptId], centroid, dimensions, expandedPt);
                double funcVal_expanded = func(expandedPt);

                if(funcVal_expanded < funcVal_reflected){
                    // Accept expansion
                    for(int dimId = 0; dimId < dimensions; ++dimId) {
                        simplex[high_ptId][dimId] = expandedPt[dimId];
                    }
                    funcVals[high_ptId] = funcVal_expanded;
                }
                else{
                    // Reject expansion and accept reflection
                    for(int dimId = 0; dimId < dimensions; ++dimId){
                        simplex[high_ptId][dimId] = reflectedPt[dimId];
                    }
                    funcVals[high_ptId] = funcVal_reflected;
                }
            }
            else{
                // Reflection wasn’t good
                if(funcVal_reflected < funcVals[high_ptId]){
                    // Ok, accept reflection
                    for(int dimId = 0; dimId < dimensions; ++dimId) {
                        simplex[high_ptId][dimId] = reflectedPt[dimId];
                    }
                    funcVals[high_ptId] = funcVal_reflected;
                }
                else{
                    // Try  contraction
                    simplex_contraction(simplex[high_ptId], centroid, dimensions, reflectedPt);
                    double funcVal_contracted = func(reflectedPt);

                    if(funcVal_contracted < funcVals[high_ptId]){
                        // Accept contraction
                        for(int dimId = 0; dimId < dimensions; ++dimId){
                            simplex[high_ptId][dimId]=reflectedPt[dimId];
                        }
                        funcVals[high_ptId] = funcVal_contracted;
                    }
                    else{
                        // Do reduction
                        simplex_reduction(simplex, dimensions, low_ptId);
                        simplex_initiate(func, simplex, funcVals, dimensions, &high_ptId, &low_ptId, centroid);
                    }
                }
            }
            numOfSteps++;
        }
        return numOfSteps ;
  }
