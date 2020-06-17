#include <math.h>
#include "polar.h"
double a2, x2;

double 
phiint(double u) {
  return(1/sqrt(a2-(1-2*u)*u*u*x2));
}


double
phix(double MoverR, double x) {

  if (MoverR==0) {
    return(asin(x));
  }
  a2=(1-2*MoverR)*MoverR*MoverR;
  x2=x*x;

  return(x*qromb(phiint,0,MoverR));
  
}
