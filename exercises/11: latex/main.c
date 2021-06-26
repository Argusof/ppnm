#include<stdio.h>
#include<math.h>

double myExp(double);

int main(){

	for( double x = 1./16; x <= 8; x += 1./16 ){
		printf("%g %g %g\n",x,myExp(x),exp(x));
	}

	return 0;
}

