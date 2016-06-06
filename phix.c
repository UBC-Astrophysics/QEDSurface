#include <math.h>
#include "polar.h"

#pragma omp threadprivate(a2, x2)
double a2, x2;

void
phix_derivs(double u, double s[], double ds[]) {
  ds[1]=1/sqrt(a2-(1-2*u)*u*u*x2);
}

double
phix(double MoverR, double x) {

  double s[2];
  int nok, nbad;
  if (MoverR==0) {
    return(asin(x));
  }
  a2=(1-2*MoverR)*MoverR*MoverR;
  x2=x*x;

  s[1]=0;
  odeint(s,1,0.0,MoverR,EPS,MoverR/10,MoverR/1e3,&nok,&nbad,phix_derivs,bsstep);
  
  return(x*s[1]);
  
}
