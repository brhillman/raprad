/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/******************************************************************************/

#define P0	101325.0	/* Pascals */

/******************************************************************************/
/* Routine to compute the total rayleigh transmisson given the rayleigh       */
/* optical depth and the cosine of the incident zenith angle.                 */

double rayleigh_total_transmission(double depth, double mu)
{

  double t, gamma;

  gamma = sqrt(3.0) / 2.0;
  mu = fabs(mu);

  t = exp((double)(-depth/mu)) - 1.0;
  t = (1.0 + t*(0.5 - gamma*mu)) / (1.0 + gamma*depth);

  return t;

}

/******************************************************************************/
/* Routine to compute the rayleigh optical depth of a standard atmosphere     */
/* above a certain pressure (height) given the pressure and wavelength.       */
/* The pressure must be in pascals */

double rayleigh_optical_depth(double p, double lambda)
{

  double tau;

  tau = 0.008569*exp(-4.0*log(lambda));
  tau *= (1 + 0.0113*exp(-2.0*log(lambda)) + 0.00013*exp(-4.0*log(lambda)));
  tau *= (p / P0);

  return tau;

}

/******************************************************************************/
/******************************************************************************/
