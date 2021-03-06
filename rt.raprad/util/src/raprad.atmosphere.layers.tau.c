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

/*--------------------------------------------------------------------*/

#define LARGE_NUMBER    1.E+37

/*--------------------------------------------------------------------*/

void
  mfp_aerosolextinction(),
  mfp_cloudextinction(),
  mfp_rayleighscatter(),
  mfp_spectralmodel_kato(),
  mfp_spectralmodel_mlawer_lw(),
  mfp_spectralmodel_pollack();

double
  compute_cloud_watervapordensity(),
  mxratio_numberdensity_precmpercm(),
  log_interpolate();

/**********************************************************************/

void
atmosphere_layers_tau(i, j, ps, pp, sm, atm, c, rt, d)
int
  i,
  j;
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
Rt
  *rt;
Brdf
  *d;
{

  int
    k,
    kk,
    n;

  double
    bi[3],
    cvd,
    ucvd,
    plog,
    pressure[2],
    z,
    p,
    dz;

  /*------------------------------------------------------------------*/
  /* The following indices are for the Mlawer model only.             */
  /*------------------------------------------------------------------*/

  static int
    layfirst = 1,
    laytrop,
    layswtch,
    laylow;

  /*------------------------------------------------------------------*/
  /* Start from the first layer closest to the ground and work our    */
  /* way to the top layer.  Note that 0 corresponds to the lowest     */
  /* plane in the atmosphere and ps->gridnumf corresponds to the      */
  /* highest plane in the problem.                                    */
  /*------------------------------------------------------------------*/

  for (k=1; k<=ps->gridnumf; k++) {

    /*----------------------------------------------------------------*/
    /* Layer geometrical thickness                                    */
    /*----------------------------------------------------------------*/

    dz  = ps->gridf[k] - ps->gridf[k-1];

    /*----------------------------------------------------------------*/
    /* cvd stands for the "cloud vapor density" and ucvd stands for   */
    /* the cvd in units of precipitable centimeters.  If a layer      */
    /* contains a cloud, then we assume the layer is saturated and we */
    /* these two variables to store the actual amount of water at the */
    /* saturation vapor pressure.                                     */
    /*----------------------------------------------------------------*/

    cvd  = 0.;
    ucvd = 0.;

    /*----------------------------------------------------------------*/
    /* Calculate the extinction, scattering and absorption            */
    /* coefficients for each of the constituents at this wavelength   */
    /* and location.                                                  */
    /*----------------------------------------------------------------*/

    for (n=1; n<=ps->cnumber; n++) {

      /*--------------------------------------------------------------*/
      /* Intitialize the extinction, scattering and absorption        */
      /* coefficients before we actually calculate their value for    */
      /* each constituent.                                            */
      /*--------------------------------------------------------------*/

      bi[0] = 0.; /*    Extinction    */
      bi[1] = 0.; /*     Scatter      */
      bi[2] = 0.; /*    Absorption    */

      /*--------------------------------------------------------------*/
      /* Make sure that the current layer falls within the physical   */
      /* location of the constituent.  In the following IF statements */
      /* mfp stands for mean free path, but it really stands for the  */
      /* extinction, scattering and absorption coefficients.          */
      /*--------------------------------------------------------------*/

      if ((c[n].base_height<=ps->gridf[k-1])&&(ps->gridf[k]<=c[n].top_height)) {

        /*------------------------------------------------------------*/
        /* RAYLEIGH SCATTERING EXTINCTION                             */
        /*------------------------------------------------------------*/

        if (!strcmp(c[n].name,"RayleighScatter")) {

          /*----------------------------------------------------------*/

          mfp_rayleighscatter(i, j, k, ps, pp, sm, atm, &c[n], &cvd, &ucvd, bi);

          /*----------------------------------------------------------*/

        }

        /*------------------------------------------------------------*/
        /* RAYLEIGH SCATTERING MARTONCHIK                             */
        /*------------------------------------------------------------*/

        else if (!strcmp(c[n].name,"RayleighScatterMartonchik")) {

          /*----------------------------------------------------------*/

          mfp_rayleighscatter_martonchik(i, j, k, ps, pp, sm, atm, &c[n], &cvd, &ucvd, bi);

          /*----------------------------------------------------------*/

        }

        /*------------------------------------------------------------*/
        /* GASEOUS ABSORPTION EXTINCTION USING POLLACK MODEL          */
        /*------------------------------------------------------------*/

        else if (!strcmp(c[n].name,"GaseousAbsorptionSpectralModelPollack")) {

          /*----------------------------------------------------------*/

          mfp_spectralmodel_pollack(i, j, k, ps, pp, sm, atm, &c[n], &cvd, &ucvd, bi);

          /*----------------------------------------------------------*/

        }

        /*------------------------------------------------------------*/
        /* GASEOUS ABSORPTION EXTINCTION USING KATO MODEL             */
        /*------------------------------------------------------------*/

        else if (!strcmp(c[n].name,"GaseousAbsorptionSpectralModelKato")) {

          /*----------------------------------------------------------*/

          mfp_spectralmodel_kato(i, j, k, ps, pp, sm, atm, &c[n], &cvd, &ucvd, bi);

          /*----------------------------------------------------------*/

        }

        /*------------------------------------------------------------*/
        /* GASEOUS ABSORPTION EXTINCTION USING MLAWER LONGWAVE MODEL  */
        /*------------------------------------------------------------*/

        else if (!strcmp(c[n].name,"GaseousAbsorptionSpectralModelMlawerLw")) {

          /*----------------------------------------------------------*/
          /* Look for three key heights called LAYTROP, LAYSWTCH, and */
          /* LAYLOW.  For Mlawer's longwave spectral model, the       */
          /* pressure is in mb, so 1.0E-2 converts Pa to mb.          */
          /*----------------------------------------------------------*/

          if (layfirst) {

            laytrop = 0 ;
            layswtch = 0 ;
            laylow = 0 ;

            for (kk=1; kk<=ps->gridnumf; kk++) {

              z = (ps->gridf[kk-1] + ps->gridf[kk]) / 2.;
              p = log_interpolate(atm->pressure, atm->altitude, atm->numlayers, z)/100.;
              plog = log(p);

              if  (plog > 4.56) {
                ++laytrop ;
              }

              if (plog >= 5.76) {
                ++layswtch ;
              }

              if (plog >= 6.62) {
                ++laylow ;
              }

            }

            layfirst = 0;

          }

          /*----------------------------------------------------------*/

          mfp_spectralmodel_mlawer_lw(i, j, k, ps, pp, sm, atm, &c[n],
                     rt, &cvd, &ucvd, d, bi, laytrop, layswtch, laylow);

          /*----------------------------------------------------------*/

        }

        /*------------------------------------------------------------*/
        /* OZONE ABSORPTION MARTONCHIK                                */
        /*------------------------------------------------------------*/

        else if (!strcmp(c[n].name,"OzoneAbsorptionMartonchik")) {

          /*----------------------------------------------------------*/

          mfp_ozoneabsorption_martonchik(i, j, k, ps, pp, sm, atm, &c[n], &cvd, &ucvd, bi);

          /*----------------------------------------------------------*/

        }

        /*------------------------------------------------------------*/
        /* CLOUD SCATTERING AND ABSORPTION EXTINCTION                 */
        /*------------------------------------------------------------*/

        else if (!strcmp(c[n].name,"CloudExtinction")) {

          /*----------------------------------------------------------*/

          mfp_cloudextinction(i, j, k, ps, pp, sm, atm, &c[n], &cvd, &ucvd, bi);

          /*----------------------------------------------------------*/

        }

        /*------------------------------------------------------------*/
        /* AEROSOL SCATTERING AND ABSORPTION EXTINCTION               */
        /*------------------------------------------------------------*/

        else if (!strcmp(c[n].name, "AerosolExtinction")) {

          /*----------------------------------------------------------*/

          mfp_aerosolextinction(i, j, k, ps, pp, sm, atm, &c[n], &cvd, &ucvd, bi);

          /*----------------------------------------------------------*/

        }

        /*------------------------------------------------------------*/
        /* DONE with this constituent at this height.                 */
        /*------------------------------------------------------------*/

      }

      /*--------------------------------------------------------------*/
      /* Save the extinction, scattering and absorption coefficients  */
      /* into their permanent place holders.                          */
      /*--------------------------------------------------------------*/

      rt->layers_kext_c[n][k] = bi[0];   /* Extinction */
      rt->layers_ksca_c[n][k] = bi[1];   /*  Scatter   */
      rt->layers_kabs_c[n][k] = bi[2];   /* Absorption */

      /*--------------------------------------------------------------*/
      /* Calculate the constituent optical depth.                     */
      /*--------------------------------------------------------------*/

      rt->layers_tau_c[n][k] = (bi[0]*dz);

      /*--------------------------------------------------------------*/
      /* DONE with this constituent so go to the next one.            */
      /*--------------------------------------------------------------*/

    }

    /*----------------------------------------------------------------*/
    /* DONE with all constituents at this height so go to the next    */
    /* height.                                                        */
    /*----------------------------------------------------------------*/

  }

  /*------------------------------------------------------------------*/
  /* DONE with all of the layers.                                     */
  /*------------------------------------------------------------------*/

}

/**********************************************************************/
/**********************************************************************/
