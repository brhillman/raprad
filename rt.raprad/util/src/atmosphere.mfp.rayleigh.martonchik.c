/**********************************************************************/
/**********************************************************************/

#include <stdio.h>
#include <math.h>

#include "../include/PhotonSpace.h"
#include "../include/PhotonPartition.h"
#include "../include/SpectralModel.h"
#include "../include/Atmosphere.h"
#include "../include/Constituents.h"

/**********************************************************************/

void mfp_rayleighscatter_martonchik(
   int i, int j, int k, 
   PhotonSpace *ps, PhotonPartition *pp, SpectralModel *sm, Atmosphere *atm, Constituents *c, 
   double *cvd, double *ucvd, double *bi
)
{
  double taua, taub, x;

  /*--------------------------------------------------------------*/
  /* Compute the optical depth of the atmosphere at this          */
  /* wavelength up to the bottom of the layer from the surface.   */
  /*--------------------------------------------------------------*/

  if ((c->base_height < 0.) || (c->scale_height < 0.)) {
    printf("Constituent %s base height or scale height is invalid.\n",c->name);
  }
  x = (ps->gridf[k-1] - c->base_height) / c->scale_height;

  if (c->tau_user[i] < 0.) {
    printf("Constituent %s optical depth is invalid.\n",c->name);
  }
  taua = c->tau_user[i] * (1. - exp(-x));

  /*--------------------------------------------------------------*/
  /* Compute the optical depth of the atmosphere at this          */
  /* wavelength up to the top of the layer from the surface.      */
  /*--------------------------------------------------------------*/

  x = (ps->gridf[k] - c->base_height) / c->scale_height;

  taub = c->tau_user[i] * (1. - exp(-x));

  /*--------------------------------------------------------------*/
  /* Compute the extinction coeffcient for this layer.            */
  /*--------------------------------------------------------------*/

  bi[0] = (taub - taua) / (ps->gridf[k] - ps->gridf[k-1]);

  /*--------------------------------------------------------------*/
  /* Make the assignments to the remaining coefficients.          */
  /*--------------------------------------------------------------*/

  bi[1] = bi[0];
  bi[2] = 0.;

  /*--------------------------------------------------------------*/
  /* DONE                                                         */
  /*--------------------------------------------------------------*/

}

/**********************************************************************/
/**********************************************************************/
