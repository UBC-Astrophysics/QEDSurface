#define EPS 1.0e-3
#define R0 1e6
#define c 3e10
#define TWO_PI 6.283185307
#define PI 3.141592654
#define TWO_OVER_FIFTEEN_ALPHA_OVER_4PI 7.7425e-5
#define BCRIT 4.41405e13
#define SOFT 0.005
#define Planck_h 4.13566e-18 /* keV s */
#define MSUN 1.4770885e5 /* cm */

#define X 1
#define Y 2
#define Z 3

#define S1 1
#define S2 2
#define S3 3
#define PHI 4			/* position of photon */
#define RADIUS 5                /* radius of photon */
#define LENGTH 6		/* length along the trajectory */
#define OMDL 7
#define NVAR 7

void
integrate_path(double omega0_p, double mass_p, 
	       double radius0_p, double b_p, 
	       double alpha, double beta,
	       double *s, int verbose);

void
calcfpsi(double u, double *f, double *psi, double *f2, double *psi2);

double
phix(double MoverR, double x);

void
derivs(double lambda, double *s, double *ds);


#ifdef USE_RK4
void rk4(double y[], double dydx[], int n,
	 double x, double h, double yout[],
	 void (*derivs)(double,double *,double *));
#else
void rkqc(double y[], double dydx[], int nv, double *xx, double htry,
	      double eps, double yscal[], double *hdid, double *hnext,
	      void (*derivs)(double, double [], double []));
void bsstep(double y[], double dydx[], int nv, double *xx, double htry,
	      double eps, double yscal[], double *hdid, double *hnext,
	      void (*derivs)(double, double [], double []));
void odeint(double ystart[], int nvar, double x1, double x2,
	      double eps, double h1, double hmin, int *nok, int *nbad,
	      void (*derivs)(double, double [], double []),
	      void (*rkqs)(double [], double [], int, double *, double, double,
			   double [], double *, double *, void (*)(double, double [], double [])));
#endif

/* double atof(char *); */
