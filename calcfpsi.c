#include <math.h>

void
calcfpsi(double u, double *f, double *psi, double *f2, double *psi2) {
  double lghld, onemu,grr, dum1, dum2;

  onemu=1-u;
  grr=1/sqrt(onemu);

  if (u>0.1) {
     lghld=log(onemu);

     /* calculate the dipole auxiliary functions f and psi */
     *f = -3/(dum2=(dum1=u*u)*u)*(lghld+u*(1+0.5*u));
     *psi = 3/(dum1)*(1/onemu+2*lghld/u+1)/grr;
     
     /* calculate the quadrapole auxiliary functions f2 and psi2 */
     *f2 = 10.0/3.0/(dum2*=dum1)*((18*u-24)*lghld+u*(-24.0 + u*(6+u)));
     *psi2 = 10/dum2*grr*(6*(2-u)*onemu*lghld+u*(12 + u*(-12+u)));

  } else {
     *f = 1 + u*(0.75+u*(0.6+u*0.5));
     *psi = 1 + u*(1+u*(0.925+u*0.875));
     
     *f2 = 1 + u*(1.333333333333+u*1.428571429*(1+u));
     *psi2 = 1 + u*(1.5 + u*(1.732142857+u*1.830357143));
  }
  
}
