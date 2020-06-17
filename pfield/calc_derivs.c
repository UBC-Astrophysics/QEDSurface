#include <math.h>
#include <stdio.h>
#include "polar.h"

double mass, b;
double cosbeta=0.6, sinbeta=0.8;
/* bigq=2/3 yields a quadrapole with B^2 at its pole and equator    */
/* equal to that of the dipole at its poles and equator (this is in */
/* flat spacetime limit) */
double m[4]={0,1.0,0.0,0.0}, q[4]={0,0,0,1.0}, bigq=0, radius0,
  deltahat[4]={0,0,0,1.0}, delta=0;
double omega_g[4], magomega_g, omega0, omega_soft;
double azimuth_global;
double rdotm_global;

void
printtestdata() {
  int i;
  printf("mass = %g, b = %g, cosbeta = %g, sinbeta = %g\n",
	 mass,b,cosbeta,sinbeta);
  printf("radius0 = %g omega0 = %g omega_soft = %g\n",radius0,omega0,omega_soft);
  printf("m = [ ");
  for (i=1;i<=3;i++) {
    printf("%g ",m[i]);
  }
  printf("]\nbigq = %g q = [ ",bigq);
  for (i=1;i<=3;i++) {
    printf("%g ",q[i]);
  }
  printf("]\n");

}

#if 0
double a2, x2;
double 
phiint(double u) {
  return(1/sqrt(a2-(1-2*u)*u*u*x2));
}
#endif


/* calculate the value of Omega as a function of  */
/* radius and phi of the photon */
double
calc_omega(double *omega, double radius, double phi, double *ds) {
   double magomega;
   double f, psi, lghld, f2, psi2;
   double u, grr;
   /* vectors and norms */
   double r[4], k2, k[4], b2, b[4], phihat[4], kcr2, kcrossr[4], bp2;
   /* dot products */
   double rdotm, rdotq, rdotdelta, rhld, bdotk, bperpdotkcrossr, bperpdotr, 
     kdotr, bdotr;
   /* angles */
   double sin2angkb;
   double cosangkcrossrbperp, sinangkcrossrbperp;
   double num, denom;
   int i;

   /* calculate the distortion of the field due to GR */
   u=2*mass/radius;
   grr=1/sqrt(1-u);

   u=0;
#if 1
   calcfpsi(u,&f,&psi,&f2,&psi2);
#else
   if (u>0.1) {
     lghld=log(onemu);

     /* calculate the dipole auxiliary functions f and psi */
     f = -3/(rdotm=(bdotk=u*u)*u)*(lghld+u*(1+0.5*u));
     psi = 3/(bdotk)*(1/onemu+2*lghld/u+1)/grr;

     /* calculate the quadrapole auxiliary functions f2 and psi2 */
     f2 = 10.0/3.0/(rdotm*=bdotk)*((18*u-24)*lghld+u*(-24.0 + u*(6+u)));
     psi2 = 10/rdotm*grr*(6*(2-u)*onemu*lghld+u*(12 + u*(-12+u)));

   } else {
     f = 1 + u*(0.75+u*(0.6+u*0.5));
     psi = 1 + u*(1+u*(0.925+u*0.875));

     f2 = 1 + u*(1.333333333333+u*1.428571429*(1+u));
     psi2 = 1 + u*(1.5 + u*(1.732142857+u*1.830357143));
   }
#endif   

   /* calculate r hat */
   r[X]=cos(phi);
   r[Y]=sin(phi);

   /* rotate from the trajectory plane */
   r[Z]=r[Y]*sinbeta;
   r[Y]*=cosbeta;

   /* calculate phi hat */
   phihat[X]=-sin(phi);
   phihat[Y]=cos(phi);

   /* rotate from the trajectory plane */
   phihat[Z]=phihat[Y]*sinbeta;
   phihat[Y]*=cosbeta;

   /* calculate k */
   k2=0;
   kdotr=0;
   for (i=X;i<=Z;i++) {
     k[i]=r[i]*ds[RADIUS]*grr+phihat[i]*ds[PHI]*radius;
     k2+=k[i]*k[i];  
     kdotr+=k[i]*r[i];
   }

   /* calculate B */

   /* first calculate the dipole component */
   /* calculate rdotm */
   rdotm=0;
   for (i=X;i<=Z;i++) {
     rdotm+=r[i]*m[i];
   }
   rdotm_global=rdotm;
#ifdef VERBOSE
   printf("rdotm: %g\n",rdotm);
#endif

   u=radius/radius0;
   bp2=1/u/u/u;
   rhld=(2*f+psi)*rdotm;
   for (i=X;i<=Z;i++) {
#if 1
     b[i]=bp2*(rhld*r[i]-psi*m[i]);
#else
     /* force b to be in the sphere just to test */
     b[i]=bp2*(-f*(m[i]-rdotm*r[i]));
#endif
#ifdef VERBOSE
     printf("%g ",b[i]);
#endif
   }

   if (bigq!=0) {
     /* now work on the quadrapole */
     /* calculate rdotq */
     rdotq=0;
     for (i=X;i<=Z;i++) {
       rdotq+=r[i]*q[i];
     }
     
     bp2=3*bigq/u/u/u/u;
     rhld=(1.5*f2+psi2)*rdotq*rdotq-0.5*f2;
     for (i=X;i<=Z;i++) {
       b[i]+=bp2*(rhld*r[i]-rdotq*psi2*q[i]);
#ifdef VERBOSE
       printf("%g ",b[i]);
#endif
     }
   }

   if (delta!=0) {
     /* now work on the offset dipole */
     /* calculate rdotdelta */
     rdotq=0;
     for (i=X;i<=Z;i++) {
       rdotdelta+=r[i]*deltahat[i];
     }
     
     bp2=3*delta/u/u/u/u;
     rhld=(3*f2+psi2)*rdotm*rdotdelta;
     rdotm=rdotq*psi2;
     for (i=X;i<=Z;i++) {
       b[i]+=bp2*(rhld*r[i]-psi2*(rdotdelta*m[i]+rdotm*deltahat[i]));
#ifdef VERBOSE
       printf("%g ",b[i]);
#endif
     }
   }
#if 0
   b[1]=m[1];
   b[2]=m[2];
   b[3]=m[3];
#endif

   /* calculate b2 */
   b2=0;
   for (i=X;i<=Z;i++) {
     b2+=b[i]*b[i];
   }


#ifdef VERBOSE
   printf("%g %g %g\n",b2,bp2,f);
#endif

   /* calculate Bdotk, Bdotr */
   bdotk=0;
   bdotr=0;
   for (i=X;i<=Z;i++) {
     bdotk+=b[i]*k[i];
     bdotr+=b[i]*r[i];
   }

   
   num=bdotk-kdotr*bdotr;
   // is the magnetic field pointing into the surface?
   // if it is, use -B
   
   if (bdotr<0) {
     num*=-1;
   }
   denom=(b2-bdotr*bdotr)*(k2-kdotr*kdotr);
   /* calculate angle between btangential and ktangential; azimuth */
   if ( denom==0 )
     azimuth_global=0.0;
   else {
     num/=sqrt(denom);
     if (num<=-1) {
       azimuth_global=M_PI;
     } else if (num>=1) {
       azimuth_global=0.0;
     } else {
       azimuth_global=acos(num);
     }
   }

   /* calculate sin of angle between k and b */
   sin2angkb=1-bdotk*bdotk/k2/b2;
   
   /* calculate Bperp */
   bp2=0;
   for (i=X;i<=Z;i++) {
     b[i]-=bdotk*k[i]/k2;
     bp2+=b[i]*b[i];
   }

   /* calculate kcrossr */

   kcrossr[X]=k[Y]*r[Z]-k[Z]*r[Y];
   kcrossr[Y]=k[Z]*r[X]-k[X]*r[Z];
   kcrossr[Z]=k[X]*r[Y]-k[Y]*r[X];

   /* calculate kcrossr2 */
   kcr2=0;
   for (i=X;i<=Z;i++) {
     kcr2+=kcrossr[i]*kcrossr[i];
   }
   
   /* calculate bperpdotkcrossr and bperpdotr */
   bperpdotkcrossr=0;
   bperpdotr=0;
   for (i=X;i<=Z;i++) {
     bperpdotkcrossr+=b[i]*kcrossr[i];
     bperpdotr+=b[i]*r[i];
   }
   cosangkcrossrbperp=bperpdotkcrossr/sqrt(kcr2*bp2);

   sinangkcrossrbperp=sqrt(1-cosangkcrossrbperp*cosangkcrossrbperp);
   if (bperpdotr<0) {
     sinangkcrossrbperp*=-1;
   }

   /* the magnitude of omega is proportional to the square of the strength 
      of the field perpendicular to k, plus the photon has a higher frequency
      near to the star than we observe by (1-2*M/r)^{-1/2} */

   magomega=bp2*omega0*grr;
   if (magomega>omega_soft) magomega=omega_soft;
   
   omega[1]=2*cosangkcrossrbperp*cosangkcrossrbperp-1;
   omega[2]=-2*cosangkcrossrbperp*sinangkcrossrbperp;
   omega[3]=0;

   return(magomega);
}

void
derivs(double lambda, double *s, double *ds) {
  double grr, radius, r2, dsdlambda;
  int i;
  
  grr=1/(1-2*mass/(radius=s[RADIUS]));
  r2=radius*radius;

  ds[RADIUS]=sqrt(1-b*b/grr/r2);
  ds[PHI]=b/r2;
  ds[LENGTH]=dsdlambda=sqrt(grr);
  
  magomega_g=calc_omega(omega_g,radius,s[PHI],ds)/dsdlambda;

  ds[S1] = (omega_g[2]*s[S3] - omega_g[3]*s[S2])*magomega_g;
  ds[S2] = (omega_g[3]*s[S1] - omega_g[1]*s[S3])*magomega_g;
  ds[S3] = (omega_g[1]*s[S2] - omega_g[2]*s[S1])*magomega_g;

  ds[OMDL] = magomega_g;

#ifdef VERBOSE
  for (i=S1;i<=S3;i++) {
    printf("%g %g %g\n",s[i],ds[i],omega_g[i]);
  }
#endif

#if 0
  printf("%g %g %g %g %g %g %g\n",r,s[1],s[2],s[3],ds[1],ds[2],ds[3]);
#endif
}

