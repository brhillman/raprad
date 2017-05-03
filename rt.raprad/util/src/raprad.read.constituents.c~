/****************************************************************/
/****************************************************************/

#include <stdio.h>

#include "../include/numrec.nrutil.h"
#include "../include/Constituents.h"
#include "../include/PhotonPartition.h"
#include "../include/PhotonSpace.h"

/****************************************************************/

Constituents
*read_constituents(constituentsfile, pp, ps)
char
  *constituentsfile;
PhotonPartition
  *pp;
PhotonSpace
  *ps;
{
  int
    ca,
    i,
    j;

  char
    input_type[2];

  Constituents
    *c;

  FILE
    *finput,
    *fopen();

  /*-------------------------------------------------------------*/
  /* Open the constituents input file.                           */
  /*-------------------------------------------------------------*/

  finput = fopen(constituentsfile,"r");

  /*-------------------------------------------------------------*/
  /* Read the number of atmospheric constituents.                */
  /*-------------------------------------------------------------*/

  if (setup(finput)) { fscanf(finput, "%d", &ps->cnumber); }

  /*-------------------------------------------------------------*/

  c = (Constituents *) malloc(ps->cnumber * sizeof(Constituents));
  c -= 1;

  /*-------------------------------------------------------------*/

  for (i=1; i<=ps->cnumber; i++) {

    /*-----------------------------------------------------------*/

    if (setup(finput)) { myfscanf(c[i].name, finput); }
    else {
      printf("User must specify a valid constituent name.\n");
      exit(0);
    }

    c[i].base_height = -1.;
    if (setup(finput)) { fscanf(finput, "%lf", &c[i].base_height); }

    c[i].top_height = -1.;
    if (setup(finput)) { fscanf(finput, "%lf", &c[i].top_height); }

    c[i].scale_height = -1.;
    if (setup(finput)) { fscanf(finput, "%lf", &c[i].scale_height); }

    c[i].cdensity = -1.;
    if (setup(finput)) { fscanf(finput, "%lf", &c[i].cdensity); }

    c[i].tau_abs_user = dvector(1, pp->numintervals);
    c[i].tau_sca_user = dvector(1, pp->numintervals);
    c[i].tau_user     = dvector(1, pp->numintervals);
    c[i].sig_ext_user = dvector(1, pp->numintervals);
    c[i].sig_sca_user = dvector(1, pp->numintervals);
    c[i].sig_abs_user = dvector(1, pp->numintervals);
    c[i].w0_user      = dvector(1, pp->numintervals);
    c[i].p_g_user     = dvector(1, pp->numintervals);

    for (j=1; j<=pp->numintervals; j++) {
      c[i].tau_abs_user[j] = -1.;
      c[i].tau_sca_user[j] = -1.;
      c[i].tau_user[j]     = -1.;
      c[i].sig_ext_user[j] = -1.;
      c[i].sig_sca_user[j] = -1.;
      c[i].sig_abs_user[j] = -1.;
      c[i].w0_user[j]      = -1.;
      c[i].p_g_user[j]     = -1.;
    }

    if (setup(finput)) {
      for (j=1; j<=pp->numintervals; j++) {
        fscanf(finput, "%lf", &c[i].tau_abs_user[j]);
        if (c[i].tau_abs_user[j] == -1) { break; }
      }
    }

    if (setup(finput)) {
      for (j=1; j<=pp->numintervals; j++) {
        fscanf(finput, "%lf", &c[i].tau_sca_user[j]);
        if (c[i].tau_sca_user[j] == -1) { break; }
      }
    }

    if (setup(finput)) {
      for (j=1; j<=pp->numintervals; j++) {
        fscanf(finput, "%lf", &c[i].tau_user[j]);
        if (c[i].tau_user[j] == -1) { break; }
      }
    }

    if (setup(finput)) {
      for (j=1; j<=pp->numintervals; j++) {
        fscanf(finput, "%lf", &c[i].sig_ext_user[j]);
        if (c[i].sig_ext_user[j] == -1) { break; }
      }
    }

    if (setup(finput)) {
      for (j=1; j<=pp->numintervals; j++) {
        fscanf(finput, "%lf", &c[i].sig_sca_user[j]);
        if (c[i].sig_sca_user[j] == -1) { break; }
      }
    }

    if (setup(finput)) {
      for (j=1; j<=pp->numintervals; j++) {
        fscanf(finput, "%lf", &c[i].sig_abs_user[j]);
        if (c[i].sig_abs_user[j] == -1) { break; }
      }
    }

    if (setup(finput)) {
      for (j=1; j<=pp->numintervals; j++) {
        fscanf(finput, "%lf", &c[i].w0_user[j]);
        if (c[i].w0_user[j] == -1) { break; }
      }
    }

    if (setup(finput)) {
      for (j=1; j<=pp->numintervals; j++) {
        fscanf(finput, "%lf", &c[i].p_g_user[j]);
        if (c[i].p_g_user[j] == -1) { break; }
      }
    }

    if (setup(finput)) { myfscanf(c[i].phasefcn_file, finput); }

    if (setup(finput)) {
      myfscanf(c[i].phasefcn_format, finput);
      if ( (strcmp(c[i].phasefcn_format,"ascii")) &&
           (strcmp(c[i].phasefcn_format,"binary"))    ) {
        printf("Invalid choice for the phase function format.\n");
        exit(0);
      }
    }

    c[i].p_which = -1;
    if (setup(finput)) {
      fscanf(finput, "%d", &c[i].p_which);
      if ((c[i].p_which != -1)&&(c[i].p_which != 0)&&
          (c[i].p_which !=  1)&&(c[i].p_which != 2)   ) {
        printf("Invalid choice for phase function.\n");
        exit(0);
      }
    }

    c[i].p_explicit_number = -1;
    if (setup(finput)) { fscanf(finput, "%d", &c[i].p_explicit_number); }

    c[i].p_legendre_number = -1;
    if (setup(finput)) { fscanf(finput, "%d", &c[i].p_legendre_number); }

    if (setup(finput)) {
      myfscanf(c[i].es, finput);
      if ( (strcmp(c[i].es,"water")) && (strcmp(c[i].es,"ice")) &&
           (strcmp(c[i].es,"none"))    ) {
        printf("Invalid choice for saturation over water or ice.\n");
        exit(0);
      }
    }

    /*-----------------------------------------------------------*/
    /* Height information.                                       */
    /*-----------------------------------------------------------*/

    if (c[i].base_height < 0.) {
      c[i].base_height = ps->zm;
    }

    if (c[i].top_height < 0.) {
      c[i].top_height = ps->zp;
    }

    if (c[i].scale_height < 0.) {
      c[i].scale_height = c[i].top_height - c[i].base_height;
    }

    /*-----------------------------------------------------------*/
    /* Based on the input parameters try to sort out as much as  */
    /* possible about the constituent optical depths and single  */
    /* scatter albedos.                                          */
    /*-----------------------------------------------------------*/

    for (j=1; j<=pp->numintervals; j++) {

      if ( (0. <= c[i].tau_abs_user[j]) && (0. <= c[i].tau_sca_user[j]) ) {
        c[i].tau_user[j] = c[i].tau_abs_user[j] + c[i].tau_sca_user[j];
        c[i].w0_user[j] = c[i].tau_sca_user[j] / c[i].tau_user[j];
      }

      else if ( (0. <= c[i].tau_user[j]) && (0. <= c[i].w0_user[j]) ) {
        c[i].tau_sca_user[j] = c[i].w0_user[j] * c[i].tau_user[j];
        c[i].tau_abs_user[j] = c[i].tau_user[j] - c[i].tau_sca_user[j];
      }

      else if ( (0. <= c[i].tau_sca_user[j]) && (0. <= c[i].w0_user[j]) ) {
        c[i].tau_user[j] = c[i].tau_sca_user[j] / c[i].w0_user[j];
        c[i].tau_abs_user[j] = c[i].tau_user[j] - c[i].tau_sca_user[j];
      }

      else if ( (0. <= c[i].tau_abs_user[j]) && (0. <= c[i].w0_user[j]) ) {
        c[i].tau_user[j] = c[i].tau_abs_user[j] / (1. - c[i].w0_user[j]);
        c[i].tau_sca_user[j] = c[i].tau_user[j] - c[i].tau_abs_user[j];
      }

      if (c[i].tau_user[j] < 0.) { break; }

    }

    /*-----------------------------------------------------------*/
    /* DONE with this constituent so go to the next one.         */
    /*-----------------------------------------------------------*/

  }

  /*-------------------------------------------------------------*/
  /* Close the constituents input file.                          */
  /*-------------------------------------------------------------*/

  fclose(finput);

  /*-------------------------------------------------------------*/
  /* Return the constituents structure to the calling routine.   */
  /*-------------------------------------------------------------*/

  return ((Constituents *) c);

  /*-------------------------------------------------------------*/

}

/****************************************************************/
/****************************************************************/
