#include <math.h>
#include <stdio.h>
#include "polar.h"

extern double mass, b;
extern double cosbeta, sinbeta;
/* bigq=2/3 yields a quadrapole with B^2 at its pole and equator    */
/* equal to that of the dipole at its poles and equator (this is in */
/* flat spacetime limit) */
extern double m[4], radius0;
extern double omega_g[4], magomega_g, omega0, omega_soft;

#pragma omp threadprivate(rdotm_start, azimuth_start)
double rdotm_start;
double azimuth_start;

void
outputdata(double lambda, double *s, double *ds) {
  int i;
  printf("%g ",lambda);
  for (i=1;i<=3;i++) {
    printf("%g ",omega_g[i]*magomega_g);
  }
  for (i=S1;i<=LENGTH;i++) {
    printf("%g ",s[i]);
  }
  printf("\n");
}



void
integrate_path(double omega0_p, double mass_p, 
	       double radius0_p, double b_p,
	       double alpha, double beta,
	       double *s, int verbose)
{
#ifdef USE_RK4
  double sout[4];
#endif

  double ds[NVAR+1];
  double lam, lamstep, lamoutput, lamstop;
  extern double rdotm_global, azimuth_global;
  int nok, nbad, stepcnt;

  omega0=omega0_p;  mass=mass_p;  radius0=radius0_p; b=b_p;
  cosbeta=cos(beta); sinbeta=sin(beta);

  m[1]=cos(alpha);
  m[2]=sin(alpha);
  m[3]=0;
   
  s[LENGTH]=0;
  s[RADIUS]=radius0;
  s[OMDL]=0;

#if 0
  /* calculate the initial value of phi */
  MoverR=mass/radius0;
  ghld=1-2*MoverR;
  rinf=radius0/sqrt(ghld);
  x=b/rinf;
  a2=ghld*MoverR*MoverR;
  x2=x*x;

  s[PHI]=-x*qromb(phiint,0,MoverR);
#else
  s[PHI]=-phix(mass/radius0,b/radius0*sqrt(1-2*mass/radius0));
#endif

  /* calculate omega_soft, i.e. the maximum magnitude of Omega */
  omega_soft=5/(radius0*SOFT);

  /* calculate the initial four-velocity of the photon along with Omega */
  lam=0;
  derivs(lam,s,ds);
  rdotm_start=rdotm_global;
  azimuth_start=azimuth_global;
  /* set the polarization to be in an eigenstate */
  lamstep=0;
  for (nok=1;nok<=3;nok++)  {
    s[nok+S1-1]=omega_g[nok]+1e-10;
    lamstep+=s[nok+S1-1]*s[nok+S1-1];
  }

  /* normalize the polarization vector */
  lamstep=sqrt(lamstep);
  for (nok=S1;nok<=S3;nok++) {
    s[nok]/=lamstep;
  }

  lamstop=radius0*100;

  if (verbose) { 
    outputdata(lam,s,ds);
  }

  for (stepcnt=1;s[RADIUS]<lamstop;lam+=lamstep,stepcnt++) {
    if (magomega_g*radius0>0.1) {
      lamstep=0.1/magomega_g;
    } else {
      lamstep=0.1*radius0;
    }
    derivs(lam+lamstep,s,ds);
#ifdef USE_RK4
    rk4(s,ds,NVAR,lam,lamstep,sout,derivs);
    for (nok=1;nok<=NVAR;nok++) {
      s[nok]=sout[nok];
    }
#else
    odeint(s,NVAR,lam,lam+lamstep,EPS,lamstep,lamstep/1e3,&nok,&nbad,derivs,bsstep);
#endif
    
    if ((verbose) && (stepcnt%10==0)) {
      outputdata(lam+lamstep,s,ds);
    }
    
  }
  if (verbose) {
    outputdata(lam,s,ds);
  }
}



