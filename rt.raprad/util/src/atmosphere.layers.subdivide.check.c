/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "../include/Atmosphere.h"


/******************************************************************************/

void check_atmosphere_layers_subdivide(atmnew)
  Atmosphere
    *atmnew;
{
  int
    i,
    j;

  FILE
    *fpta,
    *fopen();

  /*--------------------------------------------------------------------------*/
  /* Open check.atmosphere.layers.subdivide                                   */
  /*--------------------------------------------------------------------------*/

  if ((fpta=fopen("../../output/check.atmosphere.layers.subdivide","w+"))==NULL) {
    printf("check configuration file cannot be opened for some reason - EXITING!\n");
    exit(1);
  }

  /*--------------------------------------------------------------------------*/
  /* Open check.atmosphere.layers.subdivide                                   */
  /*--------------------------------------------------------------------------*/

  fprintf(fpta, "[numlayers] %5d\n\n", atmnew->numlayers);
 
  fprintf(fpta, "[altitude][pressure][temperature][airdens][dryairdens][mxratioh2o][mxratioco2] [mxratioo3][mxratioo2][densh2o][densco2][denso3][denso2]\n\n");

  for (i=0; i<atmnew->numlayers; i++) {
    fprintf(fpta,    " %3d", i);
    fprintf(fpta, " %9.2e", atmnew->altitude[i]);
    fprintf(fpta, " %9.2e", atmnew->pressure[i]);
    fprintf(fpta, " %9.2e", atmnew->temperature[i]);
    fprintf(fpta, " %9.2e", atmnew->dryairdens[i]);
    fprintf(fpta, " %9.2e", atmnew->mxratioh2o[i]);
    fprintf(fpta, " %9.2e", atmnew->mxratioco2[i]);
    fprintf(fpta, " %9.2e", atmnew->mxratioo3[i]);
    fprintf(fpta, " %9.2e", atmnew->mxratioo2[i]);
    fprintf(fpta, " %9.2e", atmnew->densh2o[i]);
    fprintf(fpta, " %9.2e", atmnew->densco2[i]);
    fprintf(fpta, " %9.2e", atmnew->denso3[i]);
    fprintf(fpta, " %9.2e", atmnew->altitude[i]);
    fprintf(fpta, " %9.2e", atmnew->denso2[i]);
    if (i<atmnew->numlayers-1) {
      fprintf(fpta, " %9.2e", atmnew->airmass[i]);
    }
    fprintf(fpta, "\n");
  }

  /*---------------------------------------------------------------------*/
  /* CLOSE check.atmosphere.layers.subdivide.                            */
  /*---------------------------------------------------------------------*/

  fclose(fpta);

  /*---------------------------------------------------------------------*/
  /* DONE with the interpolated atmosphere so exit this routine.         */
  /*---------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
