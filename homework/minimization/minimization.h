#ifndef HAVE_MINIMIZATION_H
#define HAVE_MINIMIZATION_H

#include <gsl/gsl_vector.h>

void numeric_gradient(double func(gsl_vector*), gsl_vector* minimum, gsl_vector* gradient);
void quasiNewtonMethod(double func(gsl_vector*), gsl_vector* minimum, double tolerance);
void simplex_reflection(const double* highest_pt , const double* centroid_pt, int dimensions, double* reflected_pt );
void simplex_expansion(const double* highest_pt, const double* centroid_pt , int dimensions , double* expanded_pt );
void simplex_contraction(const double* highest_pt, const double* centroid_pt, int dimensions , double* contracted_pt );
void simplex_reduction(double** simplex, int dimensions, int low_ptId );
double simplex_distance (double* first_pt , double* second_pt, int dimensions);
double simplex_size(double** simplex, int dimensions );
void simplex_update (double** simplex, const double* funcVals, int dimensions , int* high_ptId , int* low_ptId , double* centroid_pt);
void simplex_initiate (double func(double*), double** simplex, double* funcVals, int dimensions, int* high_ptId, int* low_ptId, double* centroid_pt );
int downhillsimplex(double func(double*), double** simplex , int dimensions , double simplex_sizeGoal);

#endif
