#include <math.h>
#include <stdio.h>
#include "polar.h"
#include "models.h"
#define STEPFRAC 24

extern double omega_g[4], magomega_g;

int
main(int argc, char *argv[])
{
  data_node parent_node;
  double res[NDATA];
  double args[]={85,0.560472,46.2245, 40.0};
  double atof();
  double s[NVAR+1];
  double B0, nu, k0, omega0, rinf, radius0, mass; 
  double x, b, alpha, beta, step, xmax;
  extern double rdotm_start, azimuth_start;
  double rdotm_newtonian, rdotb2, f, psi, f2, psi2;
  
  int i, nstep;

   if ( argc<5) {
     printf("Format:\n\n pfield _mu_ _nu_ _mass_ _radius_ _alpha_ [models]\n");
     return(-1);
   }

      printf("# ");
   for (i=0;i<argc;i++) {
     printf(" %s",argv[i]);
   }
   printf("\n");


   if (argc>6) {
     initializenode(&parent_node);

     parent_node.description="Parent Node";

     parent_node.parent=NULL;
     parent_node.nchildren=argc-6;
  
     if ((parent_node.children = malloc_data_node(parent_node.nchildren))==NULL) {
       return(-1);
     }
  
     for (i=0;i<parent_node.nchildren;i++) {
       loadtrjfile(argv[i+6],parent_node.children+i);
       loadgnufile(argv[i+6],parent_node.children+i);
       loadintfile(argv[i+6],parent_node.children+i);
       parent_node.children[i].parent=&parent_node;
     }
   } else {
     res[2]=res[1]=1.0;
     res[3]=0.0;
   }

   radius0=atof(argv[4]);

   B0=atof(argv[1])/radius0/radius0/radius0;
   B0/=BCRIT;

   mass=atof(argv[3]);
   alpha=atof(argv[5])*PI/180;

   /* calculate the value of |Omega| for a photon travelling along the */
   /* field near the polar cap */
   nu=atof(argv[2]);
   args[3]=log(nu);
   k0=TWO_PI*nu/c;
   omega0=TWO_OVER_FIFTEEN_ALPHA_OVER_4PI*B0*B0*k0;

   rinf=radius0/sqrt(1-2*mass/radius0);
   if (radius0<3*mass) {
     xmax=3*mass/sqrt(1-2.0/3.0)/rinf;
   } else {
     xmax=1;
   }
   printf("%g %g %g\n",xmax,radius0,mass);
   printf("#   b       beta      s1       s2      s3      mago      o1        o2      o3       mag_lat    theta    phi       X         O\n");


#ifdef PROGRESSIVE
   nstep=2;
   for (x=0.1;x<xmax;x+=0.1,nstep+=2) {
#else
   for (x=1e-5;x<xmax;x+=0.1) {
#endif
     b=x*rinf;
#ifdef DOSPOKES
     for (beta=1e-5;beta<PI+1e-4;beta+=PI/6) {
#else
#ifdef PROGRESSIVE
     step=PI/(2*nstep);
     for (beta=step*0.5;beta<2*PI;beta+=step) {
#else
     step=PI/STEPFRAC/x;
     for (beta=0.5*step;beta<PI/2+1e-4;beta+=step) {
#endif
#endif
       integrate_path(omega0,mass,radius0,b,alpha,beta,s,0);
       printf("%8.0f %8.6f",b,beta);
       for (i=S1;i<=S3;i++) {
	 printf(" %8.5f",s[i]);
       }
       printf(" %8.5g",magomega_g);
       for (i=1;i<=3;i++) {
	 printf(" %8.5f",omega_g[i]);
       }
       
       /* correct the value of rdotm so that the Newtonian calculation
	  of the angle of the field to the normal will yield the GR value */
#if 0
       calcfpsi(2*mass/radius0,&f,&psi,&f2,&psi2);
       f2=f*f; psi2=psi*psi;
       rdotb2=4*f2*rdotm_start*rdotm_start/(rdotm_start*rdotm_start*(4*f2-psi2)+psi2);

       rdotm_newtonian=sqrt( rdotb2 / ( 4 - 3*rdotb2) );
#else
       rdotm_newtonian=rdotm_start;  
#endif

       printf(" %8.4f %8.4f %8.4f",args[0]=acos(rdotm_newtonian)*180.0/PI,asin(x)*180.0/PI,
	      args[2]=azimuth_start*180.0/PI);
       if (argc>6) {
	 if (args[0]>90) args[0]=180-args[0];
	 args[1]=x;
	 evaltree(&parent_node,args,4, res);
#if 0	 
	 for (i=0;i<4;i++) {
	   printf(" %g",args[i]);
	 }
	 for (i=0;i<NDATA;i++) {
	   printf(" %g",res[i]);
	 }
#endif
	 res[2]=exp(res[2]);
	 res[3]=exp(res[3]);
       }
       printf(" %10.4e %10.4e\n",res[2],res[3]);

     }
#ifndef DOSPOKES
#ifndef PROGRESSIVE
     for (beta=PI-0.5*step;beta>PI/2+0.5*step;beta-=step) {
       integrate_path(omega0,mass,radius0,b,alpha,beta,s,0);
       printf("%8.0f %8.6f",b,beta);
       for (i=S1;i<=S3;i++) {
	 printf(" %8.5f",s[i]);
       }
       printf(" %8.6g",magomega_g);
       for (i=1;i<=3;i++) {
	 printf(" %8.5f",omega_g[i]);
       }
       printf(" %8.4f %8.4f %8.4f %8.4f\n",acos(rdotm_start)*180.0/PI,asin(x)*180.0/PI,
	      azimuth_start*180.0/PI,beta*180.0/PI);

     }
#endif
#endif
   }
   return(0);
}
