#include <math.h>
#include <stdio.h>
/*
 * If RHL is defined, take eps as a parameter, else use EPS
 */
#ifdef RHL
#  define EPS eps
#else
#  define EPS 1.0e-6
#endif
#define JMAX 20
#define JMAXP JMAX+1
#define K 5

#ifdef RHL
double qromb(func,a,b,eps)
double a,b;
double (*func)();
double eps;
#else
double qromb(func,a,b)
double a,b;
double (*func)();
#endif
{
	double ss,dss,trapzd();
	double s[JMAXP+1],h[JMAXP+1];
	int j;
	void polint(),nrerror();

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=0.25*h[j];
	}
	printf("Too many steps in routine QROMB\n");
	return(-99.0);
}

#undef EPS
#undef JMAX
#undef JMAXP
#undef K

