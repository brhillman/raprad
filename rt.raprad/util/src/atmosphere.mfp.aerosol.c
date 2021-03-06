/**********************************************************************/
/**********************************************************************/

#include "../include/PhotonSpace.h"
#include "../include/PhotonPartition.h"
#include "../include/SpectralModel.h"
#include "../include/Atmosphere.h"
#include "../include/Constituents.h"

/**********************************************************************/

void
mfp_aerosolextinction(i, j, k, ps, pp, sm, atm, c, cvd, ucvd, bi)
  int
    i,
    j,
    k;
  PhotonSpace
    *ps;
  PhotonPartition
    *pp;
  SpectralModel
    *sm;
  Atmosphere
    *atm;
  Constituents
    *c;
  double
    *cvd,
    *ucvd,
    *bi;
{
  double
    dz,
    z;

  /*--------------------------------------------------------------*/
  /* Layer thickness and mid-layer height.                        */
  /*--------------------------------------------------------------*/

  dz = ps->gridf[k] - ps->gridf[k-1];
  z  = (ps->gridf[k] + ps->gridf[k-1]) / 2.;

  /*--------------------------------------------------------------*/
  /* Caluculate the cloud extinction coefficients.                */
  /*--------------------------------------------------------------*/

  if (0. < c->tau_user[i]) {
    if ((c->base_height == -1.)||(c->top_height == -1)) {
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
  /* DONE with aerosol extinction.                                */
  /*--------------------------------------------------------------*/

}

/**********************************************************************/
/**********************************************************************/
