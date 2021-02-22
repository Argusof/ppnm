#include <stdio.h>
#include <math.h>
#include "komplex.h"

void komplex_print(char *s, komplex a) {
  printf("%s(%g,%g)\n",s, a.re, a.im);
}

komplex komplex_new (double x, double y) {
  komplex z = { x, y };
  return z;
}

void komplex_set (komplex* z, double x, double y) {
  (*z).re = x;
  (*z).im = y;
}

komplex komplex_add (komplex a, komplex b) {
  komplex result_add = { a.re + b.re, a.im + b.im };
  return result_add;
}

komplex komplex_sub (komplex a, komplex b) {
  komplex result_sub = { a.re - b.re, a.im - b.im };
  return result_sub;
}

int komplex_equal (komplex a, komplex b, double acc, double eps) {
  if( ((fabs(a.re-b.re) < acc) ||
      (fabs(a.re-b.re)/( fabs(a.re) + fabs(b.re) ) < eps/2 )) &&
      ((fabs(a.im-b.im) < acc) ||
      (fabs(a.im-b.im)/( fabs(a.im) + fabs(b.im) ) < eps/2 ))){
		return 1;
	}
	else {
		return 0;
		}
}

komplex komplex_mul (komplex a, komplex b) {
  komplex result_mul = komplex_new( (a.re * b.re) - (a.im * b.im), (a.re * b.im) + (a.im * b.re) );
  return result_mul;
}

komplex komplex_div (komplex a, komplex b) {
  komplex result_div = komplex_new( (1/(pow(b.re,2)+pow(b.im,2)))*((a.re * b.re) + (a.im * b.im)), (1/(pow(b.re,2)+pow(b.im,2)))*((a.im * b.re) - (a.re * b.im)) );
  return result_div;
}

komplex komplex_conjugate (komplex z) {
  komplex result_conjugate = { z.re, -z.im };
  return result_conjugate;
}

double komplex_abs (komplex z) {
  double result_abs = sqrt( pow(z.re,2) + pow(z.im,2) );
  return result_abs;
}

komplex komplex_exp (komplex z) {
  komplex exp_re = komplex_new( exp(z.re), 0 );
  komplex exp_im = komplex_new( cos(z.im), sin(z.im) );
  komplex result_exp = komplex_mul(exp_re, exp_im);
  return result_exp;
}

komplex komplex_sin (komplex z){
  komplex m_I        = komplex_new(0, 1);
	komplex expArg     = komplex_mul(m_I, z);
	komplex mexpArg    = komplex_mul(komplex_new(-1,0), expArg);
	komplex num 		   = komplex_sub(komplex_exp(expArg), komplex_exp(mexpArg));
	komplex den 		   = komplex_mul(komplex_new(2,0), m_I);
  komplex result_sin = komplex_div(num, den);
  return result_sin;
}

komplex komplex_cos (komplex z){
  komplex m_I        = komplex_new(0, 1);
	komplex expArg     = komplex_mul(m_I, z);
	komplex mexpArg    = komplex_mul(komplex_new(-1,0), expArg);
	komplex num 		   = komplex_add(komplex_exp(expArg), komplex_exp(mexpArg));
	komplex den 		   = komplex_new(2,0);
  komplex result_cos = komplex_div(num, den);
  return result_cos;
}

komplex komplex_sqrt (komplex z){
  komplex result_sqrt = komplex_new(sqrt( (komplex_abs(z)+z.re)/2), (z.im/fabs(z.im))*sqrt((komplex_abs(z)-z.re)/2) );
  return result_sqrt;
}
