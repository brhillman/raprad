/**********************************************************************/
/**********************************************************************/

#include "../include/PhotonSpace.h"
#include "../include/PhotonPartition.h"
#include "../include/SpectralModel.h"
#include "../include/Atmosphere.h"
#include "../include/Constituents.h"

/*--------------------------------------------------------------------*/

double rayleigh_optical_depth(double, double);

/**********************************************************************/

void mfp_rayleighscatter(
   int i, int j, int k, 
   PhotonSpace *ps, PhotonPartition *pp, SpectralModel *sm, Atmosphere *atm, Constituents *c, 
   double *cvd, double *ucvd, double *bi
)
{
  double taua, taub;

  /*--------------------------------------------------------------*/
  /* Compute the optical depth of the atmosphere at this          */
  /* wavelength up to the bottom of the layer from the surface.   */
  /*--------------------------------------------------------------*/

  taua = rayleigh_optical_depth(atm->pressure[k-1], pp->wavelength[i]);

  /*--------------------------------------------------------------*/
  /* Compute the optical depth of the atmosphere at this          */
  /* wavelength up to the top of the layer from the surface.      */
  /*--------------------------------------------------------------*/

  taub = rayleigh_optical_depth(atm->pressure[k]  , pp->wavelength[i]);

  /*--------------------------------------------------------------*/
  /* Compute the extinction coeffcient for this layer.            */
  /*--------------------------------------------------------------*/

  bi[0] = (taua - taub) / (ps->gridf[k] - ps->gridf[k-1]);

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
