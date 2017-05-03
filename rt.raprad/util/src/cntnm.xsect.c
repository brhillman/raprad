#include <stdio.h>
#include <math.h>

#define GWATERVAPOR 	18.
#define GAIR		28.96
#define A0 		6.022045E+23
#define p0		1013.
#define T0		296.

#define T220		220.
#define T280		280.

double continuum_crosssection(p,T,densh2o,dryairdens,csself,csforeign)
double
  p,
  T,
  densh2o,
  dryairdens,
  *csself,
  *csforeign;
{
  double
    cs,
    csSelfInterpolate,
    csForeignInterpolate,
    pvapor,
    selfscale,
    foreignscale;

  double
    Tdelta,
    Tratio;

  pvapor = (densh2o/dryairdens)*p;

  selfscale    = (pvapor/p0)*(T0/T);
  foreignscale = ((p-pvapor)/p0)*(T0/T);

  Tdelta = T280 - T220;
  Tratio = (T - T220)/Tdelta;

  if ((csself[0] <= 0.) || (csself[1] <= 0.)) {
    csSelfInterpolate = 0.;
  }
  else {
    csSelfInterpolate = exp(Tratio*(log(   csself[1]) - log(   csself[0])) + log(   csself[0]));
  }

  if ((csforeign[0] <= 0.) || (csforeign[1] <= 0.)) {
    csForeignInterpolate = 0.;
  }
  else {
    csForeignInterpolate = exp(Tratio*(log(csforeign[1]) - log(csforeign[0])) + log(csforeign[0]));
  }

  cs = selfscale*csSelfInterpolate + foreignscale*csForeignInterpolate;

  return ((double) cs);

}
