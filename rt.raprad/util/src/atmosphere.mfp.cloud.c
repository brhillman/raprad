/**********************************************************************/
/**********************************************************************/
#include <stdlib.h>
#include <stdio.h>

#include "../include/PhotonSpace.h"
#include "../include/PhotonPartition.h"
#include "../include/SpectralModel.h"
#include "../include/Atmosphere.h"
#include "../include/Constituents.h"

/*--------------------------------------------------------------------*/

extern double compute_cloud_watervapordensity(char *, double, Atmosphere *);
extern double mxratio_numberdensity_precmpercm(double);

/**********************************************************************/

void
mfp_cloudextinction(
   int i, int j, int k, 
   PhotonSpace *ps, PhotonPartition *pp, SpectralModel *sm, Atmosphere *atm, Constituents *c, 
   double *cvd, double *ucvd, double *bi
)
{
  double dz, z;

  /*--------------------------------------------------------------*/
  /* Layer thickness and mid-layer height.                        */
  /*--------------------------------------------------------------*/

  dz = ps->gridf[k] - ps->gridf[k-1];
  z  = (ps->gridf[k] + ps->gridf[k-1]) / 2.;

  /*--------------------------------------------------------------*/
  /* Caluculate the cloud extinction coefficients.                */
  /*--------------------------------------------------------------*/

  if (0. < c->tau_user[i]) {
    if ((c->base_height == -1.)||(c->top_height == -1.)) {
      printf("vertical boundaries of %s are not set - exiting!\n", c->name);
      exit(0);
    }
    bi[0] = c->tau_user[i]     / (c->top_height - c->base_height);
    bi[1] = c->tau_sca_user[i] / (c->top_height - c->base_height);
    bi[2] = c->tau_abs_user[i] / (c->top_height - c->base_height);
  }
  else if ((0.<c->cdensity)&&(0.<c->sig_ext_user[i])&&(0.<c->sig_sca_user[i])) {
    bi[0] = c->cdensity*c->sig_ext_user[i]*1.E-12;
    bi[1] = c->cdensity*c->sig_sca_user[i]*1.E-12;
    bi[2] = bi[0] - bi[1];
  }
  else {
    printf("Constituent %s optical depth is indeterminate.\n",c->name);
  }

  /*--------------------------------------------------------------*/
  /* Compute the mid-cloud layer water vapor density in kg m^{-3}.*/
  /*--------------------------------------------------------------*/

  *cvd = compute_cloud_watervapordensity(c->es, z, atm);

  /*--------------------------------------------------------------*/
  /* Convert the cloud vapor density into units of precipitable   */
  /* cm.                                                          */
  /*--------------------------------------------------------------*/

  *ucvd = mxratio_numberdensity_precmpercm((*cvd)*1.E-06);

  /*--------------------------------------------------------------*/
  /* DONE with cloud stuff.                                       */
  /*--------------------------------------------------------------*/

}

/**********************************************************************/
/**********************************************************************/
