#include<math.h>
double myExp(double val){
  if( val < 0    ) {
    double result  =  1.0 / myExp( -val   );
    return result;
  }
  if( val > 1./8 ) {
    double result  =  pow( myExp( val / 2 ), 2 );
    return result;
  }
  double result = 1+val*(1+val/2*(1+val/3*(1+val/4*(1+val/5*(1+val/6*(1+val/7*(1+val/8*(1+val/9*(1+val/10)))))))));
  return result;
}

