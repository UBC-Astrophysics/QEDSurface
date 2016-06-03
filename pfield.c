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
  double rdotm_newtonian, rdotb2, f, psi, f2, psi2, qtot;
  
  int i, nstep;

   if ( argc<5) {
     printf("\n\
Format:\n\n\
   pfield _mu_ _nu_ _mass_ _radius_ _alpha_ [models]\n\n\
_mu_     magnetic dipole moment in G cm^3\n\
_nu_     photon frequency in Hz\n\
_mass_   mass in cm, i.e. GM/c^2\n\
_radius_ radius in cm\n\
_alpha_  angle of magnetic moment with line of sight in degrees\n\
[models] optional spectral models to use\n\n");
     
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
   printf("# %g %g %g\n",xmax,radius0,mass);
   #
   printf("\
#  Column 1  - b, impact parameter in cm\n\
#  Column 2  - beta, angle in plane of sky between image element and magnetic moment [radians]\n\
#              beta defines the angle of the photon geodesic plane with respect to the magnetic moment\n\
#  Column 3  - s1 final Stokes Q [relative to photon geodesic plane]\n\
#  Column 4  - s2 final Stokes U\n\
#  Column 5  - s3 final Stokes V\n\
#  Column 6  - mago, initial value of Omega\n\
#  Column 7  - o1 final Omega Q [relative to geodesic plane]\n\
               the final values of Omega give the direction of the final B-field wrt geodesic plane\n\
#  Column 8  - o2 final Omega U\n\
#  Column 9  - o3 final Omega V\n\
#  Column 10 - magnetic colatitude of emission point [degrees]\n\
#  Column 11 - zenith angle [degrees]\n\
#  Column 12 - azimuth angle relative to local B-field [degrees]\n\
#  Column 13 - initial intensity in X mode\n\
#  Column 14 - initial intensity in O mode\n\
#  Column 15 - final intensity in Q [relative to projected magnetic moment]\n");

   printf("#   b       beta      s1       s2      s3      mago      o1        o2      o3   mag_colat    theta    phi       X         O        Q\n");


   nstep=1;
   for (x=0.05;x<xmax;x+=0.1,nstep+=2) {
     b=x*rinf;
     step=PI/(2*nstep);
#pragma omp parallel for schedule(guided)
     /* only do one half of the image, the other half by symmetry has
	s2->-s2, s3->-s3, o2->-o2, o3->-o3, the rest are the same */

     for (beta=step*0.5;beta<PI;beta+=step) {
       integrate_path(omega0,mass,radius0,b,alpha,beta,s,0);
       printf("%8.0f %8.6f",b,beta);
       for (i=S1;i<=S3;i++) {
	 printf(" %8.5f",s[i]);
       }
       printf(" %8.5g",magomega_g);
       qtot=0;
       for (i=1;i<=3;i++) {
         printf(" %8.5f",omega_g[i]);
	 qtot+=omega_g[i]*s[i];
       }
       
       /* correct the value of rdotm so that the Newtonian calculation
	  of the angle of the field to the normal will yield the GR value */
       rdotm_newtonian=rdotm_start;  

       printf(" %8.4f %8.4f %8.4f",args[0]=acos(rdotm_newtonian)*180.0/PI,asin(x)*180.0/PI,
	      args[2]=azimuth_start*180.0/PI);
       if (argc>6) {
	 if (args[0]>90) args[0]=180-args[0];
	 args[1]=x;
	 evaltree(&parent_node,args,4, res);
	 res[2]=exp(res[2]);
	 res[3]=exp(res[3]);
       }
       printf(" %10.4e %10.4e %10.4e\n",res[2],res[3],qtot);

     }
   }
   return(0);
}
