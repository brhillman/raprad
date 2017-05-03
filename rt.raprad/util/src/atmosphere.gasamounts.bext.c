/***************************************************************/
/* Linear interpolate a profile to a specific altitude.        */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define SMALL_NUMBER             1.E-37

#define PASCALS_TO_ATMOSPHERES   1.013E+05
#define PASCALS_TO_MILLIBARS     1.E+02

#define A0			 6.022045E+23   /* #/mole  */
#define DENSITY_WATER            1.000   	/* gm/cm^3 */
#define GWATERVAPOR		 18.016 	/* Molecular Weight of Water Vapor (g/mol) */

double mxratio_precmpercm_numberdensity();

double
bext_gasamounts_atmosphere(numlayers,u,p,height,altitude, flag_ucvd, ucvd, km, coverd)
int numlayers;
double *u;
double *p;
double *height;
double altitude;
int flag_ucvd;
double ucvd;
double km;
double coverd;
{
  double value, slope, peff, ueff, ueffcvd, bext;
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

  ueff = u[i];

  if (flag_ucvd) {
    if (ueff < ucvd) {
      ueff = ucvd;
    }
  }

  slope = (p[i+1]-p[i]) / (height[i+1]-height[i]);
  peff  = (slope*(altitude - height[i]) + p[i]) / PASCALS_TO_ATMOSPHERES;

  bext  = km*ueff*exp(coverd*log(peff));

  return ((double)bext);

}

/***************************************************************/
/***************************************************************/
