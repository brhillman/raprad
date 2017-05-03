/****************************************************************/
/****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../include/numrec.nrutil.h"
#include "../include/Angles.h"

/*---------------------------------------------------------------*/

extern int setup(FILE *);

/****************************************************************/

void read_geometry_sunzenith(char *geometryfile, Sun *suna)
{
  double pi;
  FILE *finput, *fopen();

  /*-------------------------------------------------------------*/

  pi = acos(-1.);

  /*-------------------------------------------------------------*/

  finput = fopen(geometryfile,"r");

  /*-------------------------------------------------------------*/
  /* Read in the Sun angle information.                          */
  /*-------------------------------------------------------------*/

  if (setup(finput)) { fscanf(finput, "%d", &suna->sunz_number); }

  if (suna->sunz_number != 1) {
    printf("The code is currently configured to use only one\n");
    printf("   solar zenith angle.\n");
    exit(0);
  }

  if (setup(finput)) {

    fscanf(finput, "%f", &suna->sunz_theta);

    /*-----------------------------------------------------------*/
    /* The azimuth of the Sun's radiance is 180 degrees out of   */
    /* phase with the Sun's location azimuth.                    */
    /*-----------------------------------------------------------*/

    suna->suna_theta = 0.;

    /*-----------------------------------------------------------*/

  }

  if (setup(finput)) { fscanf(finput, "%lf", &suna->sun2earth_distance); }

  if (setup(finput)) { fscanf(finput, "%f", &suna->gmt_time); }

  /*-------------------------------------------------------------*/
  /* Convert all of the angles from degrees to radians.          */
  /*-------------------------------------------------------------*/

  suna->sunz_mu = cos((pi/180.0)*suna->sunz_theta);
  suna->suna_mu = 1.;

  /*-------------------------------------------------------------*/

  fclose(finput);

  /*-------------------------------------------------------------*/

}

/****************************************************************/
/****************************************************************/
