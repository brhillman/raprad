/***************************************************************/
/* Linear interpolate a profile to a specific altitude.        */

#include <math.h>
#include "../include/Atmosphere.h"
#include <string.h>

#define RGAS                          8.31441 	/* J/mol/K */
#define A0			 6.022045E+23   /* #/mole  */

double linear_interpolate();
double clausiusclapyron_water_variable_L();
double clausiusclapyron_ice_constant_L();

double
compute_cloud_watervapordensity(es, z, atm)
char *es;
double z;
Atmosphere *atm;
{
  double t, p, cvd, svp;

  t = linear_interpolate(atm->temperature, atm->altitude, atm->numlayers, z);
  p = linear_interpolate(atm->pressure   , atm->altitude, atm->numlayers, z);

  if (strcmp(es,"water") == 0) {
    svp = clausiusclapyron_water_variable_L(t);
  }
  else if (strcmp(es,"ice") == 0) {
    svp = clausiusclapyron_ice_constant_L(t);
  }
  else if (strcmp(es,"none") == 0) {
    svp = 0.;
  }
  else {
    svp = 0.;
  }

  cvd = (svp*A0) / (RGAS*t);

  return ((double)cvd);

}

/***************************************************************/
/***************************************************************/
