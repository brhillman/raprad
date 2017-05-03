/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "../include/Atmosphere.h"

/******************************************************************************/

#define DATA_CHECK_RES 200

#define G	        9.80665 /* Acceleration Due to Gravity */

#define GAIR	       28.964   /* Molecular Weight of Air (g/mol) */
#define GCARBONDIOXIDE 44.00995 /* Molecular Weight of Carbon Dioxide (g/mol) */
#define GWATERVAPOR    18.016   /* Molecular Weight of Water Vapor (g/mol) */
#define GOXYGEN	       31.9988  /* Molecular Weight of Ozone (g/mol) */
#define GOZONE	       47.9982  /* Molecular Weight of Ozone (g/mol) */
#define GNITROUSOXIDE  44.0128  /* Molecular Weight of Nitrous Oxide (g/mol) */
#define GMETHANE       16.0426  /* Molecular Weight of Methane (g/mol) */

#define CO2_MIXING_RATIO  3.55E-04  /* Volume CO2 / Volume Air */
/* #define CO2_MIXING_RATIO  7.10E-04  Volume CO2 / Volume Air */
#define O2_MIXING_RATIO   2.0946E-01 /* Volume  O2 / Volume Air */

/******************************************************************************/

double pVnRT_numberdensity();
double pVnRT_pressure();
double numberdensity_dryair();
double hydrostatic_airmass();

double mxratio_h2o_gmperkg_numberdensity();
double mxratio_o3_gmpergm_numberdensity();
double mxratio_volpervol_numberdensity();
double mxratio_volpervol_gmperkg();

/******************************************************************************/

int setup();
void atmosphere_read_mcclatchey(file_name, atm)
char *file_name;
Atmosphere *atm;

{
  int
    countrequired,
    i,
    layer,
    n;

  double
    dryairdensity,
    vaporpressure,
    dryairpressure;

  FILE
    *fpt,
    *fopen();

  char
    c,
    **collabel;

  /*--------------------------------------------------------------------------*/
  /* Open the atmospheric file.                                               */
  /*--------------------------------------------------------------------------*/

  if( (fpt = fopen(file_name, "r") ) == NULL ) {
    printf ("Cannot find input file: %s\n", file_name);
    exit(0);
  }

  /*--------------------------------------------------------------------------*/
  /* Get the number of gases.                                                 */
  /*--------------------------------------------------------------------------*/

  while((c=getc(fpt)) != '@');

  if (setup(fpt)) { fscanf(fpt, "%d", &atm->numgases); }

  /*--------------------------------------------------------------------------*/
  /* Get the number of vertical levels at which the gases are specified.      */
  /*--------------------------------------------------------------------------*/

  while((c=getc(fpt)) != '@');

  if (setup(fpt)) { fscanf(fpt, "%d", &atm->numlayers); }

  /*--------------------------------------------------------------------------*/
  /* Allocate the space necessary to hold the atmospheric file information.   */
  /*--------------------------------------------------------------------------*/

  atm->altitude    = (double *) malloc(atm->numlayers*sizeof(double));
  atm->pressure    = (double *) malloc(atm->numlayers*sizeof(double));
  atm->temperature = (double *) malloc(atm->numlayers*sizeof(double));
  atm->airdens     = (double *) malloc(atm->numlayers*sizeof(double));
  atm->dryairdens  = (double *) malloc(atm->numlayers*sizeof(double));
  atm->airmass     = (double *) malloc(atm->numlayers*sizeof(double));
  atm->densh2o     = (double *) malloc(atm->numlayers*sizeof(double));
  atm->densco2     = (double *) malloc(atm->numlayers*sizeof(double));
  atm->denso3      = (double *) malloc(atm->numlayers*sizeof(double));
  atm->denso2      = (double *) malloc(atm->numlayers*sizeof(double));
  atm->densn2o     = (double *) malloc(atm->numlayers*sizeof(double));
  atm->densch4     = (double *) malloc(atm->numlayers*sizeof(double));
  atm->mxratioh2o  = (double *) malloc(atm->numlayers*sizeof(double));
  atm->mxratioco2  = (double *) malloc(atm->numlayers*sizeof(double));
  atm->mxratioo3   = (double *) malloc(atm->numlayers*sizeof(double));
  atm->mxratioo2   = (double *) malloc(atm->numlayers*sizeof(double));
  atm->mxration2o  = (double *) malloc(atm->numlayers*sizeof(double));
  atm->mxratioch4  = (double *) malloc(atm->numlayers*sizeof(double));

  /*--------------------------------------------------------------------------*/
  /* Set the allocated space quantities to 0.                                 */
  /*--------------------------------------------------------------------------*/

  for (i=0; i<atm->numlayers; i++) {
    atm->altitude[i]    = 0.0;
    atm->pressure[i]    = 0.0;
    atm->temperature[i] = 0.0;
    atm->airdens[i]     = 0.0;
    atm->dryairdens[i]  = 0.0;
    atm->airmass[i]     = 0.0;
    atm->densh2o[i]     = 0.0;
    atm->densco2[i]     = 0.0;
    atm->denso3[i]      = 0.0;
    atm->denso2[i]      = 0.0;
    atm->densn2o[i]     = 0.0;
    atm->densch4[i]     = 0.0;
    atm->mxratioh2o[i]  = 0.0;
    atm->mxratioco2[i]  = 0.0;
    atm->mxratioo3[i]   = 0.0;
    atm->mxratioo2[i]   = 0.0;
    atm->mxration2o[i]  = 0.0;
    atm->mxratioch4[i]  = 0.0;
  }

  /*--------------------------------------------------------------------------*/
  /* INPUT PARAMETER units are as follows:                                    */
  /*                                                                          */
  /* altitude       km                                                        */
  /*                                                                          */
  /* pressure       mb = hPa                                                  */
  /*                                                                          */
  /* temperature    K                                                         */
  /*                                                                          */
  /* mxratioh2o     Meteorology usage: This is actually the specific          */
  /*                humidity, which is the amount of watervapor in gm divided */
  /*                by the amount of moist air in kg.                         */
  /*                                                                          */
  /* mxratioo3      Mass mixing ratio: The amount of ozone in gm divided by   */
  /*                the amount of moist air in gm.                            */
  /*                                                                          */
  /* mxration2o     Volume mixing ratio in parts per million (ppmv), which is */
  /*		    the number of N2O molecules per million air molecules.    */
  /*                                                                          */
  /* mxratioch4     Volume mixing ratio in parts per million (ppmv), which is */
  /*		    the number of N2O molecules per million air molecules.    */
  /*                                                                          */
  /*                                                                          */
  /* OUTPUT PARAMETER units are as follows:                                   */
  /*                                                                          */
  /* atm->altitude            m                                               */
  /*                                                                          */
  /* atm->pressure            Pa                                              */
  /*                                                                          */
  /* atm->airdens             # m^-3                                          */
  /*                                                                          */
  /* atm->dryairden           # m^-3                                          */
  /*                                                                          */
  /* atm->densh2o             # m^-3                                          */
  /*                                                                          */
  /* atm->denso3              # m^-3                                          */
  /*                                                                          */
  /* atm->densco2             # m^-3                                          */
  /*                                                                          */
  /* atm->densn2o             # m^-3                                          */
  /*                                                                          */
  /* atm->densch4             # m^-3                                          */
  /*                                                                          */
  /* atm->mxratioh2o          g kg^-1  of moist air                           */
  /*                                                                          */
  /* atm->mxratioo3           g kg^-1  of moist air                           */
  /*                                                                          */
  /* atm->mxratioco2          g kg^-1  of moist air                           */
  /*                                                                          */
  /* atm->mxration2o          g kg^-1  of moist air                           */
  /*                                                                          */
  /* atm->mxratioch4          g kg^-1  of moist air                           */
  /*                                                                          */
  /* atm->airmass             kg m^-2  for each layer                         */
  /*                                                                          */
  /*--------------------------------------------------------------------------*/
  /* Read in each column label.  There must be a Z, P, and T label for each   */
  /* atmospheric file; all of the others are optional.                        */
  /*--------------------------------------------------------------------------*/

  while((c=getc(fpt)) != '@');

  collabel = (char **) malloc((atm->numgases+3+1)*sizeof(char *));
  for (n=1; n<=atm->numgases+3; n++) {
    collabel[n] = (char *) malloc(32*sizeof(char));
    fscanf(fpt, "%s", collabel[n]);
  }

  /*--------------------------------------------------------------------------*/
  /* Check to make sure we have valid atmospheric constituent entries.        */
  /*--------------------------------------------------------------------------*/

  countrequired = 0;

  for (n=1; n<=atm->numgases+3; n++) {

    if      (!strcmp(collabel[n],"Z")) {
      ++countrequired;
      continue;
    }
    else if (!strcmp(collabel[n],"P")) {
      ++countrequired;
      continue;
    }
    else if (!strcmp(collabel[n],"T")) {
      ++countrequired;
      continue;
    }
    else if (!strcmp(collabel[n],"H2O")) {
      continue;
    }
    else if (!strcmp(collabel[n],"O3")) {
      continue;
    }
    else if (!strcmp(collabel[n],"N2O")) {
      continue;
    }
    else if (!strcmp(collabel[n],"CH4")) {
      continue;
    }
    else {
      printf("invalid atmospheric file - exiting!\n");
      exit(0);
    }

  }

  if (countrequired != 3) {
    printf("Height (Z), Pressure (P), or Temperature (T) field is missing - exiting!\n");
    exit(0);
  }

  /*--------------------------------------------------------------------------*/
  /* Loop over each atmospheric layer.                                        */
  /*--------------------------------------------------------------------------*/

  while((c=getc(fpt)) != '@');

  for (i=0; i<atm->numlayers; i++) {

    /*------------------------------------------------------------------------*/
    /* Scan in this atmospheric level of data in the appropriate order.       */
    /*------------------------------------------------------------------------*/

    for (n=1; n<=atm->numgases+3; n++) {

      if      (!strcmp(collabel[n],  "Z")) {

        /*--------------------------------------------------------------------*/
        /* Get altitude and convert to meters.                                */
        /*--------------------------------------------------------------------*/

        fscanf(fpt,"%lf", atm->altitude+i);

        atm->altitude[i] *= 1000.;

        /*--------------------------------------------------------------------*/
        /* Done with altitude.                                                */
        /*--------------------------------------------------------------------*/

      }
      else if (!strcmp(collabel[n],  "P")) {

        /*--------------------------------------------------------------------*/
        /* Get pressure and convert to Pascal.                                */
        /*--------------------------------------------------------------------*/

        fscanf(fpt,"%lf", atm->pressure+i);

        atm->pressure[i] *= 100.;

        /*--------------------------------------------------------------------*/
        /* Done with pressure.                                                */
        /*--------------------------------------------------------------------*/

      }
      else if (!strcmp(collabel[n],  "T")) {

        /*--------------------------------------------------------------------*/
        /* Get temperature.                                                   */
        /*--------------------------------------------------------------------*/

        fscanf(fpt,"%lf", atm->temperature+i);

        /*--------------------------------------------------------------------*/
        /* Done with temperature.                                             */
        /*--------------------------------------------------------------------*/

      }
      else if (!strcmp(collabel[n],"H2O")) {

        /*--------------------------------------------------------------------*/
        /* Get water vapor specific humidity.                                 */
        /*--------------------------------------------------------------------*/

        fscanf(fpt,"%lf", atm->mxratioh2o+i);

        /*--------------------------------------------------------------------*/
        /* Done with water vapor specific humidity.                           */
        /*--------------------------------------------------------------------*/

      }
      else if (!strcmp(collabel[n], "O3")) {

        /*--------------------------------------------------------------------*/
        /* Get ozone to wet air mixing ratio.                                 */
        /*--------------------------------------------------------------------*/

        fscanf(fpt,"%lf", atm->mxratioo3+i);

        /*--------------------------------------------------------------------*/
        /* Done with ozone to wet air mixing ratio.                           */
        /*--------------------------------------------------------------------*/

      }
      else if (!strcmp(collabel[n],"N2O")) {

        /*--------------------------------------------------------------------*/
        /* Get nitrous oxide to wet air mixing ratio.                         */
        /*--------------------------------------------------------------------*/

        fscanf(fpt,"%lf", atm->mxration2o+i);

        /*--------------------------------------------------------------------*/
        /* Done with water vapor specific humidity.                           */
        /*--------------------------------------------------------------------*/

      }
      else if (!strcmp(collabel[n],"CH4")) {

        /*--------------------------------------------------------------------*/
        /* Get methane to wet air mixing ratio.                               */
        /*--------------------------------------------------------------------*/

        fscanf(fpt,"%lf", atm->mxratioch4+i);

        /*--------------------------------------------------------------------*/
        /* Done with water vapor specific humidity.                           */
        /*--------------------------------------------------------------------*/

      }

    }

    /*------------------------------------------------------------------------*/
    /* Calculate air density (wet air).                                       */
    /*------------------------------------------------------------------------*/

    atm->airdens[i] = pVnRT_numberdensity(atm->temperature[i],atm->pressure[i]);

    /*------------------------------------------------------------------------*/
    /* Compute water vapor number density.                                    */
    /*------------------------------------------------------------------------*/

    atm->densh2o[i] = mxratio_h2o_gmperkg_numberdensity(atm->mxratioh2o[i], atm->airdens[i]);

    /*------------------------------------------------------------------------*/
    /* Calculate dry air density.                                             */
    /*------------------------------------------------------------------------*/

    vaporpressure      = pVnRT_pressure(atm->temperature[i], atm->densh2o[i]);
    dryairpressure     = atm->pressure[i] - vaporpressure;
    dryairdensity      = pVnRT_numberdensity(atm->temperature[i], dryairpressure);
    atm->dryairdens[i] = dryairdensity;

    /*------------------------------------------------------------------------*/
    /* Compute ozone number density.                                          */
    /*------------------------------------------------------------------------*/

    atm->denso3[i]     = mxratio_o3_gmpergm_numberdensity(atm->mxratioo3[i], atm->airdens[i]);
    atm->mxratioo3[i] *= 1000.;

    /*------------------------------------------------------------------------*/
    /* Compute various nitrous oxide quantities.  The quantity read in is     */
    /* a volume mixing ratio in parts per million.  Because the units are     */
    /* ppmv, we must introduce the factor of 1.0E-06 to convert from number   */
    /* of nitrous oxide molecules per million air molecules to number of      */
    /* nitrous oxide molecules per number of air molecules.                   */
    /*------------------------------------------------------------------------*/

    atm->mxration2o[i] = atm->mxration2o[i] * 1.0E-06;
    atm->densn2o[i]    = mxratio_volpervol_numberdensity(atm->mxration2o[i], atm->airdens[i]);
    atm->mxration2o[i] = mxratio_volpervol_gmperkg(atm->mxration2o[i], GNITROUSOXIDE);

    /*------------------------------------------------------------------------*/
    /* Compute various methane quantities.  The quantity read in is a volume  */
    /* mixing ratio in parts per million.  Because the units are ppmv         */
    /* we must introduce the factor of 1.0E-06 to convert from number of      */
    /* of methane molecules per million air molecules to number of methane    */
    /* molecules per number of air molecules.                                 */
    /*------------------------------------------------------------------------*/

    atm->mxratioch4[i] = atm->mxratioch4[i] * 1.0E-06;
    atm->densch4[i]    = mxratio_volpervol_numberdensity(atm->mxratioch4[i], atm->airdens[i]);
    atm->mxratioch4[i] = mxratio_volpervol_gmperkg(atm->mxratioch4[i], GMETHANE);

    /*------------------------------------------------------------------------*/
    /* Compute various carbon dioxide quantities.                             */
    /*------------------------------------------------------------------------*/

    atm->mxratioco2[i] = CO2_MIXING_RATIO;
    atm->densco2[i]    = mxratio_volpervol_numberdensity(atm->mxratioco2[i], atm->airdens[i]);
    atm->mxratioco2[i] = mxratio_volpervol_gmperkg(atm->mxratioco2[i], GCARBONDIOXIDE);

    /*------------------------------------------------------------------------*/
    /* Compute various oxygen quantities.                                     */
    /*------------------------------------------------------------------------*/

    atm->mxratioo2[i] = O2_MIXING_RATIO;
    atm->denso2[i]    = mxratio_volpervol_numberdensity(atm->mxratioo2[i], atm->airdens[i]);
    atm->mxratioo2[i] = mxratio_volpervol_gmperkg(atm->mxratioo2[i], GOXYGEN);

    /*------------------------------------------------------------------------*/
    /* DONE with this layer so go to the next one.                            */
    /*------------------------------------------------------------------------*/

  }

  /*--------------------------------------------------------------------------*/
  /* DONE with all layers, so now get the airmass in each layer.              */
  /*--------------------------------------------------------------------------*/

  for (i=0; i<atm->numlayers-1; i++) {
    atm->airmass[i] = hydrostatic_airmass(atm->pressure[i]-atm->pressure[i+1]);
  }

  /*--------------------------------------------------------------------------*/
  /* CLOSE the file containing the atmospheric profile information.           */
  /*--------------------------------------------------------------------------*/

  fclose(fpt);

  /*--------------------------------------------------------------------------*/
  /* RETURN to the calling routine.                                           */
  /*--------------------------------------------------------------------------*/

  return;

  /*--------------------------------------------------------------------------*/
  /* DONE with this routine.                                                  */
  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
