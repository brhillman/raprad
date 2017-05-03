/***************************************************************/
/* Linear interpolate a profile to a specific altitude.        */
/***************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/***************************************************************/

#define SMALL_NUMBER  1.E-37

double
linear_interpolate(profile, height, number, altitude)
double *profile;
double *height;
int number;
double altitude;
{
  double value, slope;
  int i;

  i = number - 1;

  while ((altitude<(height[i]-SMALL_NUMBER)) && (i >= 0)) { --i; }

  if ((i==(number - 1)) && (altitude>(height[i]+SMALL_NUMBER))) {
    printf("Some linear interpolation was out of range and we\n");
    printf("currently do not allow extrapolations.  Exiting! 1\n");
    exit(1);
  }

  if (i<0) {
    printf("Some linear interpolation was out of range and we\n");
    printf("currently do not allow extrapolations.  Exiting! 2\n");
    exit(1);
  }

  /* linear interpolation fit; OK if in the troposphere */

  slope = (profile[i+1]-profile[i]) / (height[i+1]-height[i]);
  value = slope*(altitude - height[i]) + profile[i];

  return ((double)value);

}

/***************************************************************/
/* Logarithmically interpolate a profile to a specific altitude.*/
/***************************************************************/

double
log_interpolate(profile, height, number, altitude)
double *profile;
double *height;
int number;
double altitude;
{
  double value, slope;
  int i;

  i = number - 1;

  while ((altitude<(height[i]-SMALL_NUMBER)) && (i >= 0)) { --i; }

  /*  
  printf("%3i\n",i);
  printf("%6.2f\n",altitude);
  printf("%6.2f\n",height[i]);*/

  if ((i==(number - 1)) && (altitude>(height[i]+SMALL_NUMBER))) {
    printf("Some linear interpolation was out of range and we\n");
    printf("currently do not allow extrapolations.  Exiting! 3\n");
    exit(1);
  }

  if (i<0) {
    printf("Some linear interpolation was out of range and we\n");
    printf("currently do not allow extrapolations.  Exiting! 4\n");
    exit(1);
  }

  /* logarithm interpolation fit */

  slope = (log(profile[i+1]+SMALL_NUMBER)-log(profile[i]+SMALL_NUMBER)) / (height[i+1]-height[i]);
  value = slope*(altitude - height[i]) + log(profile[i]+SMALL_NUMBER);
  value = exp(value);

  return ((double)value);

}

/***************************************************************/
/***************************************************************/
