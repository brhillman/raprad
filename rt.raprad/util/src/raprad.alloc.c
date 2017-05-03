/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "../include/numrec.nrutil.h"
#include "../include/PhotonPartition.h"
#include "../include/PhotonSpace.h"
#include "../include/Rt1d.h"

/*----------------------------------------------------------------------------*/

void
raprad_alloc(pp, ps, rt)

PhotonPartition
  *pp;
PhotonSpace
  *ps;
Rt
  *rt;

{

  /*--------------------------------------------------------------------------*/
  /* Arrays to hold the extinction coefficients for each constituent within   */
  /* each layer.                                                              */
  /*--------------------------------------------------------------------------*/

  rt->layers_kext_c = matrix(1, ps->cnumber, 1, ps->gridnumf);
  rt->layers_ksca_c = matrix(1, ps->cnumber, 1, ps->gridnumf);
  rt->layers_kabs_c = matrix(1, ps->cnumber, 1, ps->gridnumf);

  /*--------------------------------------------------------------------------*/
  /* Arrays to hold the optical depths for each constituent within each layer.*/
  /*--------------------------------------------------------------------------*/

  rt->layers_tau_c = matrix(1, ps->cnumber, 1, ps->gridnumf);

  /*--------------------------------------------------------------------------*/
  /* The final optical depth in each layer after all of the constituent       */
  /* optical depths have been combined is given by layers_tau.  The delta     */
  /* function corrected optical depth is given by layers_tau_dfcn_o and the   */
  /* optical depth correction due to phase function smoothing is given by     */
  /* layers_tau_dfcn_f.                                                       */
  /*--------------------------------------------------------------------------*/

  rt->layers_tau = vector(1, ps->gridnumf);

  /*--------------------------------------------------------------------------*/
  /* The single scatter albedos.                                              */ 
  /*--------------------------------------------------------------------------*/

  rt->layers_w0 = vector(1, ps->gridnumf);

  /*--------------------------------------------------------------------------*/
  /* The blackbody radiation emitted within each layer that is so important   */
  /* to infrared radiative transfer calculations.                             */
  /*--------------------------------------------------------------------------*/

  rt->plankbnd  = dvector(0, pp->numintervals);

  rt->planklayd = dvector(0, ps->gridnumf);
  rt->planklayu = dvector(0, ps->gridnumf);
  rt->planklevd = dvector(0, ps->gridnumf);
  rt->planklevu = dvector(0, ps->gridnumf);

  /*--------------------------------------------------------------------------*/
  /* Moments of the Legendre polynomials (normalized), which is defined by    */
  /* the equation 2.3 of Liou et al. (1988) up to 5 terms.                    */
  /*--------------------------------------------------------------------------*/

  rt->layers_legendre_coef = vector(1, 5*ps->gridnumf);

  /*--------------------------------------------------------------------------*/
  /* Second moments of the Legendre polynomials (normalized), which is        */
  /* defined by the equation 2.3 of Liou et al. (1988) which is equal         */
  /* to asymmetry parameter.                                                  */
  /*--------------------------------------------------------------------------*/

  rt->layers_g0 = vector(1, ps->gridnumf);

  /*--------------------------------------------------------------------------*/
  /* DONE with the spectrally dependent memory allocations.                   */
  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
