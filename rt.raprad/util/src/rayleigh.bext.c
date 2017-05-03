/***************************************************************/
/* Linear interpolate a profile to a specific altitude.        */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PASCALS_TO_ATMOSPHERES  1.013E+05
#define SMALL_NUMBER            1.E-37

double indexrefraction_real_air();
double rayleigh_crosssection();

double
bext_rayleigh_atmosphere(numlayers, airdensity, height, lambda, altitude)
int numlayers;
double *airdensity;
double *height;
double lambda;
double altitude;
{
  double Na, slope, mr, sigma, value;
  int i;

  i = numlayers - 1;

  while ((altitude<(height[i]-SMALL_NUMBER)) && (i >= 0)) { --i; }

  if ((i==(numlayers - 1)) && (altitude>(height[i]+SMALL_NUMBER))) {
    printf("Some altitude was out of range and we\n");
    printf("currently do not allow extrapolations.  Exiting!\n");
    exit(1);
  }

  if (i<0) {
    printf("Some altitude was out of range and we\n");
    printf("currently do not allow extrapolations.  Exiting!\n");
    exit(1);
  }

  /* linear interpolation fit; OK if in the troposphere */


  slope = (airdensity[i+1]-airdensity[i]) / (height[i+1]-height[i]);
  Na    = slope*(altitude - height[i]) + airdensity[i];

  mr  = indexrefraction_real_air(lambda);
  sigma = rayleigh_crosssection(lambda, mr);

  value  = sigma*Na;

  return ((double)value);

}

/***************************************************************/
/***************************************************************/
