/***********************************************************************/
/***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/PhotonSpace.h"
#include "../include/SpectralModel.h"
#include "../include/Constituents.h"
#include "../include/Rt1d.h"
#include "../include/Atmosphere.h"

/***********************************************************************/

void atmosphere_layers(int i, int j, PhotonSpace *ps, SpectralModel *sm, Constituents *c, Rt *rt, Atmosphere *atm)
{
  int jj,k,kinv,n;
  double w,wght,g,g_legendre;

  /*-------------------------------------------------------------------*/
  /* COMMENT about the current geometry.  As it now stands, the        */
  /* atmospheric layers are numbered from 1 to n, where n stands for   */
  /* the number of layers in the atmosphere.  Furthermore, the layer   */
  /* numbered 1 is closest to the surface and the layers increase in   */
  /* height up to the highest altitude layer numbered n.  For the      */
  /* optical properities developed in this routine, however, we now    */
  /* relabel the layers so that layer 1 is at the top of the           */
  /* atmosphere and layer n is closest to the surface.                 */
  /*                                                                   */
  /*     	Top-To-Bottom             Bottom-To-Top                */
  /*                                                                   */
  /* TOP      	     0----------------------------n                    */
  /*             	 1                    n                        */
  /*          	     1----------------------------n-1                  */     
  /*                     2                   n-1                       */
  /*          	     2----------------------------n-2                  */
  /*                       	   .                                   */
  /*                       	   .                                   */
  /*                       	   .                                   */
  /*          	    n-1----------------------------1                   */
  /*                     n                    1                        */
  /* SURFACE         n-----------------------------0                   */
  /*                                                                   */
  /* Now proceed with the atmosphere setup from TOP TO BOTTOM!         */
  /*-------------------------------------------------------------------*/

  for (k=1; k<=ps->gridnumf; k++) {

    /*-----------------------------------------------------------------*/
    /* Index to go from bottom-to-top to top-to-bottom.                */
    /*-----------------------------------------------------------------*/

    kinv = ps->gridnumf + 1 - k;

    /*-----------------------------------------------------------------*/
    /* Change the planck function from bottom-to-top to top-to-bottom. */
    /* At the end of this routine, rt->planklayd[k] is the kth layer   */
    /* from the top and rt->planklevd[k] is the bottom boundary of     */
    /* the kth layer from the top.                                     */
    /*-----------------------------------------------------------------*/

    rt->planklayd[k]    = rt->planklayu[kinv];
    rt->planklevd[k]    = rt->planklevu[kinv-1];

    /*-----------------------------------------------------------------*/
    /* Compute the optical depth for each constituent in the layer and */
    /* add the results together to get the total optical depth for the */
    /* layer.                                                          */
    /*-----------------------------------------------------------------*/

    w         	      = 0.;
    g                 = 0.;
    rt->layers_tau[k] = 0.;

    for (n=1; n<=ps->cnumber; n++) {
      rt->layers_tau[k] += rt->layers_tau_c[n][kinv];
      w += (1. - 0.) * c[n].w0 * rt->layers_tau_c[n][kinv];
      g += c[n].p_g * c[n].w0 * rt->layers_tau_c[n][kinv];
    }

    /*-----------------------------------------------------------------*/
    /* Store the effective single scatter albedo for the layer.        */
    /*-----------------------------------------------------------------*/

    rt->layers_w0[k] = w / rt->layers_tau[k];


    /*-----------------------------------------------------------------*/
    /* Store the effective asymmetry parameter for the layer.          */
    /*-----------------------------------------------------------------*/

    rt->layers_g0[k] = g / (rt->layers_w0[k]* rt->layers_tau[k]);


    /*-----------------------------------------------------------------*/
    /* Store moments of the legendre polynomials                       */
    /*-----------------------------------------------------------------*/

    wght = 0.0;

    for (jj=1; jj<=5; jj++) {
      rt->layers_legendre_coef[jj + (k-1)*5] = 0.0;
    }

    for (n=1; n<=ps->cnumber; n++) {
      
      if (c[n].p_which == 1) {

        wght = (c[n].w0 * rt->layers_tau_c[n][kinv])
                         / (rt->layers_w0[k] * rt->layers_tau[k]);
    
        for (jj=1; jj<=5; jj++) {
          rt->layers_legendre_coef[jj+(k-1)*5] += wght*c[n].p_legendre_coef[jj];
        }

      }

    }  

  /*-------------------------------------------------------------------*/
  /* Forcing the first term of legendre coefficients id unity.         */
  /*-------------------------------------------------------------------*/

    for (jj=1; jj<=5; jj++) {
      rt->layers_legendre_coef[jj+(k-1)*5] = rt->layers_legendre_coef[jj+(k-1)*5] 
                                           / rt->layers_legendre_coef[1+(k-1)*5];
    }

  /*-------------------------------------------------------------------*/
  /* Checking normalization of legendre coefficients                   */
  /*-------------------------------------------------------------------*/

    g_legendre = rt->layers_legendre_coef[2 + (k-1)*5] / 3.0 ;

  /*-------------------------------------------------------------------*/
  /*                                                                   */
  /*if (g_legendre != rt->layers_g0[k]) {                              */
  /*  printf("Error in Legendre coefficients:\n");                     */
  /*  printf("Exiting!\n");                                            */
  /*  exit(0);                                                         */
  /*}                                                                  */
  /*-------------------------------------------------------------------*/

    /*-----------------------------------------------------------------*/
    /* DONE                                                            */
    /*-----------------------------------------------------------------*/

  }

  /*-------------------------------------------------------------------*/
  /* Additional atmosphere setup from TOP TO BOTTOM for the Mlawer     */
  /* code.                                                             */
  /*-------------------------------------------------------------------*/

  for (n=1; n<=ps->cnumber; n++) {

    if (!strcmp(c[n].name, "GaseousAbsorptionSpectralModelMlawerLw")) {

      /*---------------------------------------------------------------*/
      /* Additional atmosphere setup from TOP TO BOTTOM for the Mlawer */
      /* code.                                                         */
      /*---------------------------------------------------------------*/

      for (k=1; k<=ps->gridnumf; k++) {

        /*-------------------------------------------------------------*/
        /* Index to go from bottom-to-top to top-to-bottom.            */
        /*-------------------------------------------------------------*/

        kinv = ps->gridnumf + 1 - k;

        /*-------------------------------------------------------------*/
        /* Change the coefficients from bottom-to-top to top-to-bottom.*/
        /*-------------------------------------------------------------*/

        sm->alphad[i][j][k] = sm->alphau[i][j][kinv];

        /*-------------------------------------------------------------*/
        /* DONE with this level so go to the next one.                 */
        /*-------------------------------------------------------------*/

      }

      /*---------------------------------------------------------------*/
      /* DONE with all levels so exit the IF.                          */
      /*---------------------------------------------------------------*/

    }

  }

  /*-------------------------------------------------------------------*/
  /* DONE with the subroutine so exit to the calling routine.          */
  /*-------------------------------------------------------------------*/

}

/***********************************************************************/
/***********************************************************************/
