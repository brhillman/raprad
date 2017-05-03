/**********************************************************************/
/**********************************************************************/

#include <stdio.h>
#include <math.h>

#include "../include/PhotonSpace.h"
#include "../include/PhotonPartition.h"
#include "../include/SpectralModel.h"
#include "../include/Atmosphere.h"
#include "../include/Constituents.h"

/*--------------------------------------------------------------------*/
/* kabs_ozone is extinction (km^{-1}) for ozone assuming Eltemann     */
/* profile and a total column optical depth of unity.  The extinction */
/* is defined as the ozone optical depth per 2 km of altitude.        */
/*--------------------------------------------------------------------*/

static double kabs_ozone[27] =
{ 0.000e-0,
  1.045e-2, 8.602e-3, 6.624e-3, 6.336e-3, 6.691e-3, 
  1.027e-2, 1.822e-2, 2.801e-2, 3.023e-2, 3.579e-2,
  4.801e-2, 5.780e-2, 5.668e-2, 4.779e-2, 3.601e-2,
  2.645e-2, 2.001e-2, 1.418e-2, 1.061e-2, 7.424e-3,
  5.468e-3, 3.489e-3, 2.183e-3, 1.310e-3, 8.180e-4, 
  5.450e-4
};

/**********************************************************************/

void
mfp_ozoneabsorption_martonchik(i, j, k, ps, pp, sm, atm, c, cvd, ucvd, bi)

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
  int
    l,
    ibottom,
    itop;

  double
    dz,
    height,
    tau_ozone;

  /*--------------------------------------------------------------*/
  /* Find the index into kabs_ozone for the bottom and top of the */
  /* current atmospheric layer.                                   */
  /*--------------------------------------------------------------*/

  ibottom = ((int)(ps->gridf[k-1]/2000.)) + 1;
  itop    = ((int)(ps->gridf[k]  /2000.)) + 1;

  /*--------------------------------------------------------------*/
  /* Set height to the bottom layer height.                       */
  /*--------------------------------------------------------------*/

  height = ps->gridf[k-1];

  /*--------------------------------------------------------------*/
  /* Calculate the ozone optical depth from the bottom to the top */
  /* of the layer.                                                */
  /*--------------------------------------------------------------*/

  tau_ozone = 0.;

  for (l=ibottom; l<=itop; l++) {

    /*------------------------------------------------------------*/
    /* Check to make sure that height is less than 52 km, which   */
    /* is the maximum valid height for this parameterization.     */
    /*------------------------------------------------------------*/

    if (52000. <= height) { break; }

    /*------------------------------------------------------------*/
    /* Calculate the appropriate layer thickness.                 */
    /*------------------------------------------------------------*/

    if (l == itop) {
      dz = ps->gridf[k] - height;
    }
    else if (l == ibottom) {
      dz = 2000. - fmod(height,2000.);
    }
    else {
      dz = 2000.;
    }

    /*------------------------------------------------------------*/
    /* Optical depth for dz at the height given by index i.       */
    /*------------------------------------------------------------*/

    tau_ozone += (dz/1000.)*kabs_ozone[l];

    /*------------------------------------------------------------*/
    /* Increment the height by dz.                                */
    /*------------------------------------------------------------*/

    height += dz;

    /*------------------------------------------------------------*/
    /* Go to the next layer.                                      */
    /*------------------------------------------------------------*/

  }

  /*--------------------------------------------------------------*/
  /* Scale the layer optical depth, which is normalized to unity, */
  /* by the total column optical depth.                           */
  /*--------------------------------------------------------------*/

  tau_ozone *= c->tau_user[i];

  /*--------------------------------------------------------------*/
  /* Copy the results into the output array.                      */
  /*--------------------------------------------------------------*/

  if (0. < (height-ps->gridf[k-1])) {
    bi[0] = tau_ozone / (ps->gridf[k]-ps->gridf[k-1]);
  }

  bi[1] = 0.;
  bi[2] = bi[0];

  /*--------------------------------------------------------------*/
  /* DONE                                                         */
  /*--------------------------------------------------------------*/

}

/**********************************************************************/
/**********************************************************************/
