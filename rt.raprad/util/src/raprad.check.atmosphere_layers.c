/******************************************************************************/
/******************************************************************************/
#include <stdlib.h>
#include <stdio.h>

#include "../include/SpectralModel.h"
#include "../include/Rt1d.h"
#include "../include/Angles.h"
#include "../include/Brdf.h"


/******************************************************************************/

void check_atmosphere_layers(int i, int j, SpectralModel *sm, Rt *rt, Sun *suna, Brdf *d)
{
  int k;
  FILE *fpta, *fopen();

  /*--------------------------------------------------------------------------*/
  /* Open check.photon.partition                                              */
  /*--------------------------------------------------------------------------*/

  if ((fpta=fopen("../../output/check.atmosphere_layers","a+"))==NULL) {
    printf("check configuration file cannot be opened for some reason - EXITING!\n");
    exit(1);
  }

  /*--------------------------------------------------------------------------*/
  /* Print the solar constant, nlayers, cosine of the zenith angle and the    */
  /* surface albedo.                                                          */
  /*--------------------------------------------------------------------------*/

  fprintf(fpta, "%10.5f", sm->solarinsol[i]);

  fprintf(fpta, "%5d", rt->layers_number);

  fprintf(fpta, "%10.5f", suna->sunz_mu);

  fprintf(fpta, "%10.5f\n", d->albedo[i][1]);

  /*--------------------------------------------------------------------------*/
  /* Print the single-scattering albedo, tau and moments of legendre          */
  /* polynomials.                                                             */
  /*--------------------------------------------------------------------------*/

  for (k=1; k<=rt->layers_number; ++k) {

    fprintf(fpta, "%5d", k);

    /*------------------------------------------------------------------------*/
    /* Print the single-scattering albedo.                                    */
    /*------------------------------------------------------------------------*/

    fprintf(fpta, "%15.5e", rt->layers_w0[k]);

    /*------------------------------------------------------------------------*/
    /* Print the  tau.                                                        */
    /*------------------------------------------------------------------------*/

    fprintf(fpta, " %15.5e", rt->layers_tau[k]);

    /*------------------------------------------------------------------------*/
    /* Print the moment of the legendre polynomials.                          */
    /*------------------------------------------------------------------------*/

    fprintf(fpta, " %15.5f", rt->layers_legendre_coef[1 + (k-1)*5]);
    fprintf(fpta, " %15.5f", rt->layers_legendre_coef[2 + (k-1)*5]);
    fprintf(fpta, " %15.5f", rt->layers_legendre_coef[3 + (k-1)*5]);
    fprintf(fpta, " %15.5f", rt->layers_legendre_coef[4 + (k-1)*5]);
    fprintf(fpta, " %15.5f\n", rt->layers_legendre_coef[5 + (k-1)*5]);

  }

  fclose(fpta);

  /*--------------------------------------------------------------------------*/
  /* DONE checking variables of the atmospheric layers.                       */
  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
