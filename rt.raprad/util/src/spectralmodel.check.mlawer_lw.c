/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../include/numrec.nrutil.h"

#include "../include/PhotonPartition.h"
#include "../include/SpectralModel.h"

/******************************************************************************/

void check_spectralmodel_mlawer_lw(pp, sm)
PhotonPartition
  *pp;
SpectralModel
  *sm;
{
  FILE
    *fpt,
    *fopen();

  int
    i;

  /*--------------------------------------------------------------------------*/
  /* Open check.spectralmodel.mlawer_lw                                       */
  /*--------------------------------------------------------------------------*/

  if((fpt=fopen("../../output/check.spectralmodel.mlawer_lw","w+")) == NULL) {
    printf ("Cannot create file check.spectralmodel.mlawer_lw - Exiting!\n");
    exit(0);
  }

  /*--------------------------------------------------------------------------*/
  /* Print the number of gaseous absorption bands that are used in the        */
  /* simulation and are held in sm->numwavelngths.                            */
  /*--------------------------------------------------------------------------*/

  fprintf(fpt, "[numwavelngths] %5d\n\n", sm->numwavelngths);

  /*--------------------------------------------------------------------------*/
  /* If no band models are being used there is no need to go further so       */
  /* return to the main function.                                             */
  /*--------------------------------------------------------------------------*/

  if (sm->numwavelngths == 0) {
    fclose(fpt);
    return;
  }

  /*--------------------------------------------------------------------------*/
  /* Selectively print those bands specified in pp.                           */
  /*--------------------------------------------------------------------------*/

  for (i=1; i<=sm->numwavelngths; ++i) {

    fprintf(fpt, "%5d", i);

    fprintf(fpt, " %15s", sm->gas[i]);

    fprintf(fpt, " %5d", sm->nummolecules[i]);

    fprintf(fpt, " %5d", sm->cntnm[i]);

    fprintf(fpt, "%12.4e", sm->bandwidth[i][0]);

    fprintf(fpt, "%12.4e", sm->bandwidth[i][1]);

    fprintf(fpt, "%12.4e", sm->wavelngth[i]);

    fprintf(fpt, "%12.4e", sm->solarinsol[i]);

    fprintf(fpt, "%5d", sm->nquadpts[i][1]);

    fprintf(fpt, "%5d", sm->numintrvls[i]);

    fprintf(fpt,"\n\n");

    /*------------------------------------------------------------------------*/

  }

  /*--------------------------------------------------------------------------*/
  /* Close check.spectralmodel.kato                                           */
  /*--------------------------------------------------------------------------*/

  fclose(fpt);

  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
