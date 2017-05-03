#include <math.h>

void
phasefcns_rayleigh(pnum, p, mu)
int pnum;
double *p;
double *mu;
{
  int i;

  for (i=0; i< pnum; i++) {
    p[i]  = (3./4.)*(1. + mu[i]*mu[i]);
  }

}
