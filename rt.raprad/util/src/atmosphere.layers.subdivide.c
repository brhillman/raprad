/******************************************************************************/
/******************************************************************************/

#include <stdio.h>

#include "../include/PhotonSpace.h"
#include "../include/Atmosphere.h"

/*----------------------------------------------------------------------------*/

double linear_interpolate();
double log_interpolate();
double hydrostatic_airmass();

/******************************************************************************/

void atmosphere_layers_subdivide(ps, atm, atmnew)
  PhotonSpace
    *ps;
  Atmosphere
    *atm,
    *atmnew;
{
  int
    i;
  double
    altitude;

  atmnew->numlayers   = ps->gridnumf + 1;
  atmnew->altitude    = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->pressure    = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->temperature = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->airdens     = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->dryairdens  = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->airdensavg  = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->airmass     = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->mxratioh2o  = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->mxratioco2  = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->mxratioo3   = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->mxratioo2   = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->mxration2o  = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->mxratioch4  = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->densh2o     = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->densco2     = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->denso3      = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->denso2      = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->densn2o     = (double *) malloc(atmnew->numlayers*sizeof(double));
  atmnew->densch4     = (double *) malloc(atmnew->numlayers*sizeof(double));

  for (i=0; i<atmnew->numlayers; i++) {

    altitude		   = ps->gridf[i];
    atmnew->altitude[i]    = altitude;
    atmnew->pressure[i]    =    log_interpolate(atm->pressure,atm->altitude,atm->numlayers,altitude);
    atmnew->temperature[i] = linear_interpolate(atm->temperature,atm->altitude,atm->numlayers,altitude);
    atmnew->airdens[i]     =    log_interpolate(atm->airdens,atm->altitude,atm->numlayers,altitude);
    atmnew->dryairdens[i]  =    log_interpolate(atm->dryairdens,atm->altitude,atm->numlayers,altitude);
    atmnew->mxratioh2o[i]  =    log_interpolate(atm->mxratioh2o,atm->altitude,atm->numlayers,altitude);
    atmnew->mxratioco2[i]  =    log_interpolate(atm->mxratioco2,atm->altitude,atm->numlayers,altitude);
    atmnew->mxratioo3[i]   =    log_interpolate(atm->mxratioo3,atm->altitude,atm->numlayers,altitude);
    atmnew->mxratioo2[i]   =    log_interpolate(atm->mxratioo2,atm->altitude,atm->numlayers,altitude);
    atmnew->mxration2o[i]  =    log_interpolate(atm->mxration2o,atm->altitude,atm->numlayers,altitude);
    atmnew->mxratioch4[i]  =    log_interpolate(atm->mxratioch4,atm->altitude,atm->numlayers,altitude);
    atmnew->densh2o[i]     =    log_interpolate(atm->densh2o,atm->altitude,atm->numlayers,altitude);
    atmnew->densco2[i]     =    log_interpolate(atm->densco2,atm->altitude,atm->numlayers,altitude);
    atmnew->denso3[i]      =    log_interpolate(atm->denso3,atm->altitude,atm->numlayers,altitude);
    atmnew->denso2[i]      =    log_interpolate(atm->denso2,atm->altitude,atm->numlayers,altitude);
    atmnew->densn2o[i]     =    log_interpolate(atm->densn2o,atm->altitude,atm->numlayers,altitude);
    atmnew->densch4[i]     =    log_interpolate(atm->densch4,atm->altitude,atm->numlayers,altitude);

  }

  /*---------------------------------------------------------------------*/
  /* Calculate the air mass between atmospheric levels.  The resulting   */
  /* units are kg per meter^3.                                           */
  /*---------------------------------------------------------------------*/

  for (i=0; i<atmnew->numlayers-1; i++) {
    atmnew->airmass[i] = hydrostatic_airmass(atmnew->pressure[i]-atmnew->pressure[i+1]);
  }

  /*---------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
