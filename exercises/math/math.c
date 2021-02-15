#include <stdio.h>
#include <math.h>
#include <complex.h>

int main(void){
	double g = tgamma(5);
	double b = j0(0.5);
	double complex k = csqrt(-2);
	double complex l = cexp(I*M_PI);
	double complex m = cexp(I);
	double complex n = cpow(I,M_E);
	double complex o = cpow(I,I);
	
	printf("gamma(5) = %g\n", g); 
	printf("J(0.5)= %g\n", b);
	printf("sqrt(-2) = %g + I * %g \n", creal(k), cimag(k));
	printf("exp(i*pi) = %g + I * %g \n", creal(l), cimag(l));
	printf("exp(i) = %g + I * %g \n", creal(m), cimag(m));
	printf("i^e = %g + I * %g \n", creal(n), cimag(n));
	printf("i^i = %g + I * %g \n", creal(o), cimag(o));
	
	float x_float = 1.f/9;
	double x_double = 1./9;
	long double x_long_double = 1.L/9;
	
	printf("float: 1/9 = %.25g\n", x_float);
	printf("double: 1/9 = %.25lg\n", x_double);
	printf("float: 1/9 = %.25Lg\n", x_long_double);
	
return 0;}
