#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <time.h>
#include <math.h>

#include "jacobi.h"
#include "utilities.h"

void test_runtime(int numOfReps, int startRep, const char* my_outputFilename, const char* gsl_outputFilename, unsigned int* seed);

int main(int argc, char const *argv[]) {
  printf("\nPART A:\n");

//_______________Jacobi diagonalization with cyclic sweeps____________________

  unsigned int seed = time(NULL);
  int dim = 5;

  gsl_matrix* matrix = gsl_matrix_alloc(dim, dim);
  gsl_matrix* eigVecMat = gsl_matrix_alloc(dim, dim);
  gsl_matrix* eigValMat = gsl_matrix_alloc(dim, dim);

  set_data_symm(matrix, &seed);
  gsl_matrix_memcpy(eigValMat,matrix);

  print_matrix(dim, matrix, "Symmetric matrix A:");

  jacobiDiag(eigValMat, eigVecMat);
  print_matrix(dim, eigVecMat, "Eigenvectors:");
  print_matrix(dim, eigValMat, "Eigenvalues:");

  gsl_matrix* tmp = gsl_matrix_alloc(dim, dim);
  gsl_matrix* testIdentity = gsl_matrix_alloc(dim, dim);
  gsl_matrix* testDiagonal = gsl_matrix_alloc(dim, dim);
  gsl_matrix* testMat = gsl_matrix_alloc(dim, dim);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matrix, eigVecMat, 0.0, tmp);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, eigVecMat, tmp, 0.0, testDiagonal);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, eigVecMat, eigVecMat, 0.0, testIdentity);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, eigValMat, eigVecMat, 0.0, tmp);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigVecMat, tmp, 0.0, testMat);

  print_matrix(dim, testDiagonal, "Diagonal test, V^T * A * V = D:");
  print_matrix(dim, testIdentity, "Identity test, V^T * V:");
  print_matrix(dim, testMat, "Test of V * D * V^T = A:");


//________________________Quantum particle in a box_______________________________

  printf("\n\nPART B:\n");

// Build the hamiltonian matrix
  int numOfDivs = 50; // number of divisions
  double step = 1.0/(numOfDivs + 1);
  gsl_matrix* hamiltonian = gsl_matrix_alloc(numOfDivs, numOfDivs);
  for(int i = 0; i < numOfDivs - 1; i++){
    gsl_matrix_set(hamiltonian, i, i, -2);   // set diagonal
    gsl_matrix_set(hamiltonian, i, i+1, 1);  // above/below diagonal
    gsl_matrix_set(hamiltonian, i+1, i, 1);
  }
  gsl_matrix_set(hamiltonian, numOfDivs-1, numOfDivs-1, -2);
  gsl_matrix_scale(hamiltonian, -1 / step / step);

// Diagonalize the matrix using the Jacobi routine
  gsl_matrix* V = gsl_matrix_alloc(numOfDivs, numOfDivs);
  jacobiDiag(hamiltonian, V);

// check that the energies are correct
printf("\ncalculated vs. exact energies\n");
for (int k = 0; k < numOfDivs/3; k++){
  double exact = M_PI * M_PI * (k + 1) * (k + 1);
  double calculated = gsl_matrix_get(hamiltonian, k, k);
  printf("%i %g %g\n", k, calculated, exact);
}

// Plot the lowest energies
FILE* outputStream = fopen(argv[1], "w");


fprintf(outputStream,"0\t");
for(int k = 0; k < 3; k++){
  fprintf(outputStream,"0\t0\t");
}
fprintf(outputStream, "\n");
for(int i = 0; i < numOfDivs; i++){
   fprintf(outputStream,"%.5g\t", (i + 1.0) / (numOfDivs + 1));
   for (int j = 0; j < 3; j++) {
     double wavefunc = sin(((double)j+1.0) * M_PI * ((double)((i + 1.0) / (numOfDivs + 1)))) ;
     if (j == 1){ wavefunc *= -1;}
     fprintf(outputStream,"%.5g\t%g\t", 5.05*gsl_matrix_get(V, i, j), wavefunc);
   }
   fprintf(outputStream, "\n");
}

fprintf(outputStream,"1\t0\t0\t0\t0\t0\t0\t");

fprintf(outputStream, "\n");
fclose(outputStream);

gsl_matrix_free(matrix);
gsl_matrix_free(eigValMat);
gsl_matrix_free(eigVecMat);
gsl_matrix_free(testMat);
gsl_matrix_free(testDiagonal);
gsl_matrix_free(testIdentity);
gsl_matrix_free(hamiltonian);

//______________

int numOfReps = 260;
int startRep = 200;
test_runtime(numOfReps, startRep, argv[2], argv[3], &seed);

  return 0;
}
