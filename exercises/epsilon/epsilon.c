#include <stdio.h>
#include <limits.h>
#include <float.h>
#include "equal.h"

int main(){

	/////// Exercise 1 ///////
	printf("\n----------- EXERCISE 1 -----------\n");
	// i) Maximum representable integer
	printf("\nMaximum representable integer\n");
	
	int i=1; while(i+1>i) {i++;}
	printf("my max int = %i (while)\n",i);

	int j=1;
	for(j; j<j+1; j++) {;}
	printf("my max int = %i (for)\n",j);
	
	int k=1; 
	do k++;
	while (k<k+1);
	printf("my max int = %i (do while)\n",k);
	
	printf("INT_MAX = %d (limits.h)\n", INT_MAX);
	
	// ii) Minimum representable integer 
	printf("\nMinimum representable integer\n");
	
	int l=1; while(l-1<l) {l--;}
	printf("my min int = %i (while)\n",l);
	
	int m=1;
	for(m; m-1<m; m--) {;}
	printf("my max int = %i (for)\n",m);
	
	int n=1; 
	do n--;
	while (n-1<n);
	printf("my max int = %i (do while)\n",n);
	
	printf("INT_MIN = %d (limits.h)\n", INT_MIN);
	
	// iii) The machine epsilon 
	printf("\nThe machine epsilon\n");
	
	printf("DOUBLE\n");
	double d=1; 
	double prv_d;
	while(1+d!=1){
		prv_d = d;
		d/=2;
		} 
	prv_d*=2;
	printf("e=%.10g (while)\n",prv_d);
	
	for(d=1; 1+d!=1; d/=2){
		prv_d = d; 
		} 
		d*=2;
	printf("e=%.10g (for)\n",prv_d);
	
	d=1;
	do {prv_d = d;
		d/=2;}
		while(1+d!=1);
		prv_d*=2;
	printf("e=%.10g (do)\n",prv_d);
	
	printf("\nFLOAT\n");
	float f=1; 
	float prv_f;
	while(1+f!=1){
		prv_f = f;
		f/=2;
		} 
	prv_f*=2;
	printf("e=%.10g (while)\n",prv_f);
	
	for(f=1; 1+f!=1; f/=2){
		prv_f = f; 
		} 
		f*=2;
	printf("e=%.10g (for)\n",prv_f);
	
	f=1;
	do {prv_f = f;
		f/=2;}
		while(1+f!=1);
		prv_f*=2;
	printf("e=%.10g (do)\n",prv_f);
	
	printf("\nLONG DOUBLE\n");
	long double ld=1.; 
	long double prv_ld;
	while(1+ld!=1){
		prv_ld = ld;
		ld/=2;
		} 
	prv_ld*=2;
	printf("e=%.10Lg (while)\n",prv_ld);
	
	for(l=1; 1+ld!=1; ld/=2){
		prv_ld = ld; 
		} 
		ld*=2;
	printf("e=%.10Lg (for)\n",prv_ld);
	
	ld=1;
	do {prv_ld = ld;
		ld/=2;}
		while(1+ld!=1);
		prv_ld*=2;
	printf("e=%.10Lg (do)\n",prv_ld);
	
	/////// Exercise 2 ///////
	printf("\n----------- EXERCISE 2 -----------\n");
	
	int max=INT_MAX/2;
	
	// i) Calculate sums using floats
	printf("\nSums calculated using float\n");
	
	int fl;
	float sum_up_float = 0;
	for(fl=1; fl<=max; fl++){sum_up_float += 1.0/fl;}
	printf("sum up = %f\n",sum_up_float);

	float sum_down_float = 0;
	for(fl=1; fl<max; fl++){sum_down_float += 1.0/(max-fl);}
	printf("sum down = %f\n",sum_down_float);
	
	// ii) Explain the difference 
	
	// iii) Does this sum converge as a func. of max?
	
	// iv) Calculate sums using doubles
	printf("\nSums calculated using double\n");
	
	double sum_up_double = 0;
	for(i=1; i<=max; i++){sum_up_double += 1.0/i;}
	printf("sum up = %g\n",sum_up_double);

	double sum_down_double = 0;
	for(i=1; i<max; i++){sum_down_double += 1.0/(max-i);}
	printf("sum down = %g\n",sum_down_double);
	
	/////// Exercise 3 ///////
	
	// See function and header files.
	
	return 0;
}
