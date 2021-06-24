#include "komplex.h"
#include <stdio.h>
#include <math.h>
#include "funcs.h"
#define TINY 1e-6
#include <complex.h>

int main(){
	int testcount = 0;

	printf("-----Testing of functions from exercise 'Komplex'-----\n");

	printf("Testing komplex_set(), using komplex_print(), by initiallizing a and b as\n");
	komplex a, b;
	komplex_set(&a, 1., 2.);
	komplex_set(&b, 3., 4.);

	komplex_print("a =",a);
	komplex_print("b =",b);

	printf("\nTesting komplex_add() and komplex_equal()\n");
	printf("----------------------------------------------------\n");

	komplex aPlusB = komplex_add(a,b);
	komplex realAPlusB = {4.,6.};

	komplex_print("a + b should be   = ", realAPlusB);
	komplex_print("a + b is actually = ", aPlusB);

	if( komplex_equal(realAPlusB,aPlusB,TINY,TINY) ){
		printf("Test 'add' passed.\n");
		testcount++;
	}
	else{
		printf("Test 'add' failed: debug me, please... \n");
	}

	printf("\nTesting komplex_sub()\n");
	printf("----------------------------------------------------\n");

	komplex aSubB = komplex_sub(a,b);
	komplex realASubB = {-2.,-2.};

	komplex_print("a - b should be   = ", realASubB);
	komplex_print("a - b is actually = ", aSubB);

	if( komplex_equal(realASubB,aSubB,TINY,TINY) ){
		printf("Test 'sub' passed. \n");
		testcount++;
	}
	else{
		printf("Test 'sub' failed: debug me, please... \n");
	}

	printf("\nTesting komplex_mul()\n");
	printf("----------------------------------------------------\n");

	komplex aMulB = komplex_mul(a,b);
	komplex realAMulB = {-5.,10.};

	komplex_print("a * b should be   = ", realAMulB);
	komplex_print("a * b is actually = ", aMulB);

	if( komplex_equal(realAMulB,aMulB,TINY,TINY) ){
		printf("Test 'mul' passed.\n");
		testcount++;
	}
	else{
		printf("Test 'mul' failed: debug me, please... \n");
	}

	printf("\nTesting komplex_div()\n");
	printf("----------------------------------------------------\n");

	komplex adivB = komplex_div(a,b);
	komplex realAdivB = {11./25.,2./25.};

	komplex_print("a / b should be   = ", realAdivB);
	komplex_print("a / b is actually = ", adivB);

	if( komplex_equal(realAdivB,adivB,TINY,TINY) ){
		printf("Test 'div' passed.\n");
		testcount++;
	}
	else{
		printf("Test 'div' failed: debug me, please... \n");
	}

	printf("\nTesting komplex_conj()\n");
	printf("----------------------------------------------------\n");

	komplex aConj = komplex_conjugate(a);
	komplex realAConj = {1., -2.};

	komplex_print("conj(a) should be   = ", realAConj);
	komplex_print("conj(a) is actually = ", aConj);

	if( komplex_equal(realAConj,aConj,TINY,TINY) ){
		printf("Test 'conj' passed.\n");
		testcount++;
	}
	else{
		printf("Test 'conj' failed: debug me, please... \n");
	}


	printf("\nTesting komplex_conj()\n");
	printf("----------------------------------------------------\n");

	double aAbs = komplex_abs(a);
	double realAAbs = sqrt(pow(1,2) + pow(2,2));

	printf("|a| should be   = %lg\n", realAAbs);
	printf("|a| is actually = %lg \n", aAbs);

	if( equal(realAAbs,aAbs,TINY,TINY) ){
		printf("Test 'abs' passed.\n");
		testcount++;
	}
	else{
		printf("Test 'abs' failed: debug me, please... \n");
	}

	printf("\nTesting komplex_exp()\n");
	printf("----------------------------------------------------\n");

	komplex aExp = komplex_exp(a);
	double complex complexA = CMPLX(1, 2);
	double complex realAExp = cexp(complexA);
	komplex aRealExp = {creal(realAExp), cimag(realAExp)};

	komplex_print("exp(a) should be    = ", aRealExp);
	komplex_print("exp(a) is actually  = ", aExp);

	if( komplex_equal(aRealExp,aExp,TINY,TINY) ){
		printf("Test 'exp' passed.\n");
		testcount++;
	}
	else{
		printf("Test 'exp' failed: debug me, please... \n");
	}

	printf("\nTesting komplex_sin()\n");
	printf("----------------------------------------------------\n");

	komplex aSin = komplex_sin(a);
	double complex realASin = csin(complexA);
	komplex aRealSin = {creal(realASin), cimag(realASin)};

	komplex_print("sin(a) should be    = ", aRealSin);
	komplex_print("sin(a) is actually  = ", aSin);

	if( komplex_equal(aRealSin,aSin,TINY,TINY) ){
		printf("Test 'sin' passed.\n");
		testcount++;
	}
	else{
		printf("Test 'sin' failed: debug me, please... \n");
	}

	printf("\nTesting komplex_cos()\n");
	printf("----------------------------------------------------\n");

	komplex aCos = komplex_cos(a);
	double complex realACos = ccos(complexA);
	komplex aRealCos = {creal(realACos), cimag(realACos)};

	komplex_print("cos(a) should be    = ", aRealCos);
	komplex_print("cos(a) is actually  = ", aCos);

	if( komplex_equal(aRealCos,aCos,TINY,TINY) ){
		printf("Test 'cos' passed.\n");
		testcount++;
	}
	else{
		printf("Test 'cos' failed: debug me, please... \n");
	}

	printf("\nTesting komplex_sqrt()\n");
	printf("----------------------------------------------------\n");

	komplex aSqrt = komplex_sqrt(a);
	double complex realASqrt = csqrt(complexA);
	komplex aRealSqrt = {creal(realASqrt), cimag(realASqrt)};

	komplex_print("sqrt(a) should be    = ", aRealSqrt);
	komplex_print("sqrt(a) is actually  = ", aSqrt);

	if( komplex_equal(aRealSqrt,aSqrt,TINY,TINY) ){
		printf("Test 'sqrt' passed.\n");
		testcount++;
	}
	else{
		printf("Test 'sqrt' failed: debug me, please... \n");
	}

	if ( testcount == 10){
		printf("\n-----Done! All functions are working properly!-----\n");
	}
	else {
		printf("\n-----Error, not all functions are working properly!-----\n");
	}

	return 0;
}
