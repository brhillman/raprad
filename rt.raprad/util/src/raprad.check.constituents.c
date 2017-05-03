/****************************************************************/
/****************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "../include/numrec.nrutil.h"
#include "../include/Constituents.h"
#include "../include/PhotonPartition.h"
#include "../include/PhotonSpace.h"

/****************************************************************/

void
check_constituents(pp, ps, c)
PhotonPartition
  *pp;
PhotonSpace
  *ps;
Constituents
  *c;
{
  int
    i,
    j;

  FILE
    *finput,
    *fopen();

  /*--------------------------------------------------------------------------*/
  /* Open check.constituents                                                  */
  /*--------------------------------------------------------------------------*/

  if((finput=fopen("../../output/check.constituents","w+")) == NULL) {
    printf ("Cannot create file check.constituents - Exiting!\n");
    exit(0);
  }

  /*-------------------------------------------------------------*/
  /* Print the number of atmospheric constituents.               */
  /*-------------------------------------------------------------*/

  fprintf(finput, "[cnumber] %5d", ps->cnumber);

  /*-------------------------------------------------------------*/

  for (i=1; i<=ps->cnumber; i++) {

    /*-----------------------------------------------------------*/

    fprintf(finput, "\n\n");

    fprintf(finput, "[name] %s\n", c[i].name);

    fprintf(finput, "\t [base_height] %12.4e \n", c[i].base_height);

    fprintf(finput, "\t [top_height] %12.4e \n", c[i].top_height);

    fprintf(finput, "\t [scale_height] %12.4e \n", c[i].scale_height);

    fprintf(finput, "\t [cdensity] %12.4e \n", c[i].cdensity);

    fprintf(finput, "\t [tau_abs_user]");
    for (j=1; j<=pp->numintervals; j++) {
      fprintf(finput, " %12.4e", c[i].tau_abs_user[j]);
    }

    fprintf(finput, "\n\t [tau_sca_user]");
    for (j=1; j<=pp->numintervals; j++) {
      fprintf(finput, " %12.4e", c[i].tau_sca_user[j]);
    }

    fprintf(finput, "\n\t [tau_user]");
    for (j=1; j<=pp->numintervals; j++) {
      fprintf(finput, " %12.4e", c[i].tau_user[j]);
    }

    fprintf(finput, "\n\t [w0_user]");
    for (j=1; j<=pp->numintervals; j++) {
      fprintf(finput, " %12.4e", c[i].w0_user[j]);
    }

    fprintf(finput, "\n\t [p_g_user]");
    for (j=1; j<=pp->numintervals; j++) {
      fprintf(finput, " %12.4e", c[i].p_g_user[j]);
    }

    fprintf(finput, "\n\t [phasefcn_file] %s", c[i].phasefcn_file);

    fprintf(finput, "\n\t [phasefcn_format] %s", c[i].phasefcn_format);

    fprintf(finput, "\n\t [p_which] %5d", c[i].p_which);

    fprintf(finput, "\n\t [p_explicit_number] %5d", c[i].p_explicit_number);

    fprintf(finput, "\n\t [p_legendre_number] %5d", c[i].p_legendre_number);

    fprintf(finput, "\n\t [es] %s", c[i].es);

    /*-----------------------------------------------------------*/
    /* DONE with this constituent so go to the next one.         */
    /*-----------------------------------------------------------*/

  }

  /*-------------------------------------------------------------*/
  /* Close the constituents input file.                          */
  /*-------------------------------------------------------------*/

  fclose(finput);

  /*-------------------------------------------------------------*/

}

/****************************************************************/
/****************************************************************/
