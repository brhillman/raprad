/******************************************************************************/
/* Routine to check atmospheric thermodynamic variables */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "../include/Atmosphere.h"

/******************************************************************************/

void check_atmosphere_gasamounts_u(atm)
  Atmosphere
    *atm;
{
  int
    i;
  FILE
    *fpta,
    *fopen();

  /*--------------------------------------------------------------------------*/
  /* Open check.atmosphere.gasamounts.u                                       */
  /*--------------------------------------------------------------------------*/

  if ((fpta=fopen("../../output/check.atmosphere.gasamounts.u","w+"))==NULL) {
    printf("Cannot create check.atmosphere.gasamounts.u for some reason - EXITING!\n");
    exit(1);
  }

  /*--------------------------------------------------------------------------*/
  /* Print the gas amounts in units of cm.                                    */
  /*--------------------------------------------------------------------------*/

  fprintf(fpta, "[index][layer][airdensity(#/m^3)][uh2o(precmpercm)][uco2(cmatmpercm)][uo3(cmatmpercm)][uo2(cmatmpercm)][uh2o(# cm^-2)][uco2(# cm^-2)][uo3(# cm^-2)][uo2(# cm^-2)][un2o(# cm^-2)][uch4(# cm^-2)]\n\n");

  for (i=0; i<atm->numlayers-1; i++) {
    fprintf(fpta, " %5d", i);
    fprintf(fpta, " %12.4e", atm->altitude[i+1]-atm->altitude[i]);
    fprintf(fpta, " %12.4e", atm->airdensavg[i]);
    fprintf(fpta, " %12.4e", atm->uh2oprecm[i]);
    fprintf(fpta, " %12.4e", atm->uco2amagats[i]);
    fprintf(fpta, " %12.4e", atm->uo3amagats[i]);
    fprintf(fpta, " %12.4e", atm->uo2amagats[i]);
    fprintf(fpta, " %12.4e", atm->uh2o[i]);
    fprintf(fpta, " %12.4e", atm->uco2[i]);
    fprintf(fpta, " %12.4e", atm->uo3[i]);
    fprintf(fpta, " %12.4e", atm->uo2[i]);
    fprintf(fpta, " %12.4e", atm->un2o[i]);
    fprintf(fpta, " %12.4e", atm->uch4[i]);

    fprintf(fpta, "\n");
  }

  /*--------------------------------------------------------------------------*/
  /* CLOSE check.atmosphere.gasamounts.u.                                     */
  /*--------------------------------------------------------------------------*/

  fclose(fpta);

  /*--------------------------------------------------------------------------*/
  /* DONE with checking the atmosphere gasamounts u so exit the routine.      */
  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
