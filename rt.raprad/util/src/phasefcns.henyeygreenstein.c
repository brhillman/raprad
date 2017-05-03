#include <math.h>

void phasefcns_henyeygreenstein(int pnum, double g, double *p, double *mu)
{
  int i;
  double pi, dtheta, thetaforward, theta;

  pi = acos(-1.);

  dtheta = pi / (pnum - 1.);

  thetaforward = 0.;

  for (i=1; i<= pnum; i++) {
    theta = thetaforward + (i-1)*dtheta;
    mu[i] = cos(theta);
    p[i]  = (1. - g*g) / exp((3./2.)*log(1. + g*g - 2.*g*mu[i]));
  }

}
