#include <math.h>

#define FDELTA 0.0279     /* Molecular Anisoptropy Factor (unitless) */
#define N0     2.547E+25  /* Reference Number Density at Sea Level (#/m^3) */

double rayleigh_crosssection(lambda, m)
  double
    lambda,
    m;
{
  double
    f,
    pi,
    numerator,
    denominator;

  pi = acos(-1.0);
  lambda *= 1.E-06;

  f = (6. + 3.*FDELTA) / (6. - 7.*FDELTA);

  numerator = 8.*pi*pi*pi*(m*m - 1.)*(m*m - 1.)*f;
  denominator = 3.*lambda*lambda*lambda*lambda*N0*N0;

  return ((double) (numerator/denominator));

}
