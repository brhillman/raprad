/**************************************************************************/
/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "../include/PhotonSpace.h"
#include "../include/Constituents.h"

/**************************************************************************/

void
read_constituents_phasefcns(i, ps, c)
int
  i;
PhotonSpace
  *ps;
Constituents
  *c;
{
  int
    n;

  /*----------------------------------------------------------------------*/
  /* Loop over all of the constituents and read in their phase function   */
  /* appropriate for this particular spectral interval.                   */
  /*----------------------------------------------------------------------*/

  for (n=1; n<=ps->cnumber; n++) {

    /*--------------------------------------------------------------------*/
    /* Read in the phase function for each constituent.                   */
    /*--------------------------------------------------------------------*/

    switch (c[n].p_which) {

      /* No phase function */

      case -1:
        break;

      /* Phase function in terms of arrays of numbers. */

      case  0:

        if (!strcmp(c[n].phasefcn_format,"ascii")) {
          phasefcns_explicit_read_ascii(c[n].phasefcn_file, i, &c[n]);
        }
        else if (!strcmp(c[n].phasefcn_format,"binary")) {
          phasefcns_explicit_read_binary(c[n].phasefcn_file, i, &c[n]);
        }
        else {
          printf("Constituent %s must either have an ascii or binary format\n", c[n].name);
          printf("   specified for the phase function.\n");
          exit(0);
        }

        if (c[n].sig_ext_user[i] < 0.) {
          c[n].sig_ext_user[i] = c[n].sig_ext;
        }

        if (c[n].sig_sca_user[i] < 0.) {
          c[n].sig_sca_user[i] = c[n].sig_sca;
        }

        if (c[n].w0_user[i] < 0.) {
          c[n].w0_user[i] = c[n].sig_sca_user[i] / c[n].sig_ext_user[i];
        }

        if ((0. <  c[n].tau_user[i]) && (c[n].tau_sca_user[i] < 0.)) {
          c[n].tau_sca_user[i] = c[n].w0_user[i] * c[n].tau_user[i];
        }

        if ((0. <  c[n].tau_user[i]) && (c[n].tau_abs_user[i] < 0.)) {
          c[n].tau_abs_user[i] = c[n].tau_user[i] - c[n].tau_sca_user[i];
        }

        break;

      /* Phase function in terms of legendre polynomials. */

      case  1:

        if (!strcmp(c[n].phasefcn_format,"ascii")) {
          phasefcns_legendre_read_ascii(c[n].phasefcn_file, i, &c[n]);
        }
        else if (!strcmp(c[n].phasefcn_format,"binary")) {
          phasefcns_legendre_read_binary(c[n].phasefcn_file, i, &c[n]);
        }
        else {
          printf("Constituent %s must either have an ascii or binary format\n", c[n].name);
          printf("   specified for the phase function.\n");
          exit(0);
        }

        if (c[n].sig_ext_user[i] < 0.) {
          c[n].sig_ext_user[i] = c[n].sig_ext;
        }

        if (c[n].sig_sca_user[i] < 0.) {
          c[n].sig_sca_user[i] = c[n].sig_sca;
        }

        if (c[n].w0_user[i] < 0.) {
          c[n].w0_user[i] = c[n].sig_sca_user[i] / c[n].sig_ext_user[i];
        }

        if ((0. <  c[n].tau_user[i]) && (c[n].tau_sca_user[i] < 0.)) {
          c[n].tau_sca_user[i] = c[n].w0_user[i] * c[n].tau_user[i];
        }

        if ((0. <  c[n].tau_user[i]) && (c[n].tau_abs_user[i] < 0.)) {
          c[n].tau_abs_user[i] = c[n].tau_user[i] - c[n].tau_sca_user[i];
        }

        break;

      default:
        break;

    }

    /*--------------------------------------------------------------------*/
    /* If the user option values of the asymmetry parameter and single    */
    /* scatter albedo are not -1, then use them to overwrite the values   */
    /* just read in.  Check to make sure that these values are now        */
    /* defined at this point in the program and if they are not then exit.*/
    /*--------------------------------------------------------------------*/

    switch (c[n].p_which) {

      /* No phase function because of no scattering */

      case -1:

        c[n].w0 = 0.;
        break;

      /* Phase function in terms of legendre polynomials */

      case  0:
      case  1:
      case  2:

        if (-0.999 < c[n].p_g_user[i]) {
          c[n].p_g = c[n].p_g_user[i];
        }
        if (c[n].p_g < -0.999) {
          printf("Constituent %s scatters radiation and its asymmetry\n",c[n].name);
          printf("  parameter has not been specified.\n");
          exit(0);
        } 

        if (-0.999 < c[n].w0_user[i]) {
          c[n].w0 = c[n].w0_user[i];
        }
        if (c[n].w0 < -0.999) {
          printf("Constituent %s single scatter albedo has not been specified.\n",c[n].name);
          exit(0);
        }

        break;

      default:
        break;

    }

    /*--------------------------------------------------------------------*/
    /* At this point enough information has been specified to calculate   */
    /* the Henyey-Greenstein function if that is the option in effect.    */
    /*--------------------------------------------------------------------*/

    if (c[n].p_which == 2) {
      if (c[n].p_explicit_number == -1) {
        printf("Constituent %s uses a Henyey-Greenstein phase\n",c[n].name);
        printf("   function with an undetermined number of angles.\n");
        exit(0);
      }
      phasefcns_henyeygreenstein(c[n].p_explicit_number, c[n].p_g, c[n].p, c[n].pmu);
    }

    /*--------------------------------------------------------------------*/
    /* DONE with this constituent so go to the next one.                  */
    /*--------------------------------------------------------------------*/

  }

  /*----------------------------------------------------------------------*/
  /* DONE with all constituents.                                          */
  /*----------------------------------------------------------------------*/

}

/**************************************************************************/
/**************************************************************************/
