/**********************************************************************/
/**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/PhotonSpace.h"
#include "../include/PhotonPartition.h"
#include "../include/SpectralModel.h"
#include "../include/Atmosphere.h"
#include "../include/Constituents.h"
#include "../include/Rt1d.h"
#include "../include/Angles.h"
#include "../include/Brdf.h"


double bext_spectralmodelkato(int, double *, double *, double, double);
double bext_gasamounts_atmosphere(int, double *, double *, double, int, double, double, double);
double bext_watervapor_cntnm_atmosphere(int, double *, double *, double *, double *, double *, double, int, double, double *, double *);
double linear_interpolate(double *, double *, int, double);
double log_interpolate(double *, double *, int, double);

void rrtm_driver_setcoef_taumol_(
   int *iband, int *ig, int *lay, float *tbound, double *tz, double *tavel, double *pz, double *pavel,
   double *wkl_1, double *wkl_2, double *wkl_3, double *wkl_4, double *wkl_6, double *wkl_7,
   double *coldry, float *albedo, int *laytrop, int *layswtch, int *laylow,
   double *bi, double *gau_wt, double *plankbnd, double *planklay, double *planklev
);


/**********************************************************************/

void mfp_spectralmodel_mlawer_lw(
   int i, int j, int k, 
   PhotonSpace *ps, PhotonPartition *pp, SpectralModel *sm, Atmosphere *atm, Constituents *c, Rt *rt, 
   double *cvd, double *ucvd, Brdf *d, double *bi, 
   int laytrop, int layswtch, int laylow
)
{
  double p, plev, t, z, dz, dryaircolumn, plkbnd;
  FILE *fopen();

  /*--------------------------------------------------------------*/
  /* Mid-layer height. z and dz are in m.                         */
  /*--------------------------------------------------------------*/

  z = (ps->gridf[k-1] + ps->gridf[k]) / 2.;
  dz  = ps->gridf[k] - ps->gridf[k-1];

  /*--------------------------------------------------------------*/
  /* Mid-layer pressure and level pressure in mb.                 */
  /*--------------------------------------------------------------*/

  p = log_interpolate(atm->pressure, atm->altitude, atm->numlayers, z)/100.;
  plev = atm->pressure[k-1] / 100. ;

  /*--------------------------------------------------------------*/
  /* Mid-layer temperature.                                       */
  /*--------------------------------------------------------------*/

  t = linear_interpolate(atm->temperature, atm->altitude, atm->numlayers, z);

  /*--------------------------------------------------------------*/
  /* column amount of dry air in # cm^-2.                         */
  /*--------------------------------------------------------------*/

  dryaircolumn = log_interpolate(atm->dryairdens, atm->altitude, atm->numlayers, z);
  dryaircolumn = dryaircolumn * dz * 1.E-04;

  /*--------------------------------------------------------------*/
  /* Call Mlawer's subroutine to get the gaseous absorption       */
  /* optical depth, as opposed to getting the coefficients for    */
  /* gaseous absorption.  It seems to me that Mlawer's code needs */
  /* # cm^-2 for absorber concentration.                          */
  /*--------------------------------------------------------------*/

  rrtm_driver_setcoef_taumol_(&i, &j, &k, &d->srftemp[i][1], &atm->temperature[k-1],
    &t, &plev, &p, &atm->uh2o[k-1],  &atm->uco2[k-1], 
    &atm->uo3[k-1], &atm->un2o[k-1], &atm->uch4[k-1], 
    &atm->uo2[k-1], &dryaircolumn, &d->albedo[i][1], &laytrop, &layswtch, &laylow,
    &bi[0], &sm->alphau[i][j][k], &plkbnd, &rt->planklayu[k], &rt->planklevu[k-1]) ;


  /*--------------------------------------------------------------*/
  /* Save the Planck function source radiance for this layer.     */
  /*--------------------------------------------------------------*/

  if (k == 1) {
    rt->plankbnd[i]  = plkbnd;
  }

  /*--------------------------------------------------------------*/
  /* Based on the spectral model gas constituent, go calculate    */
  /* the extinction coefficient.                                  */
  /*--------------------------------------------------------------*/

  bi[0] /= dz ;

  /*--------------------------------------------------------------*/
  /* Make the assignments to the remaining coefficients.          */
  /*--------------------------------------------------------------*/

  bi[1] = 0.;
  bi[2] = bi[0];

  /*------------------------------------------------------------------*/
  /* DONE                                                             */
  /*------------------------------------------------------------------*/

}

/**********************************************************************/
/**********************************************************************/
