/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../include/numrec.nrutil.h"

#include "../include/PhotonPartition.h"
#include "../include/SpectralModel.h"

/******************************************************************************/

void check_spectralmodel_pollack(pp, sm)
PhotonPartition
  *pp;
SpectralModel
  *sm;
{
  FILE
    *fpt,
    *fopen();

  int
    i,
    j,
    m;

  /*--------------------------------------------------------------------------*/
  /* Open check.spectralmodel.pollack                                         */
  /*--------------------------------------------------------------------------*/

  if((fpt=fopen("../../output/check.spectralmodel.pollack","w+")) == NULL) {
    printf ("Cannot create file check.spectralmodel.pollack - Exiting!\n");
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

    fprintf(fpt, "%5d", sm->numintrvls[i]);

    for (j=1; j<=sm->nummolecules[i]; ++j) {
      fprintf(fpt, " %12.4e", sm->coverd[i][j]);
    }

    fprintf(fpt,"\n\t");

    for (j=1; j<=sm->numintrvls[i]; ++j) {
      fprintf(fpt," %12.4e", sm->alpha[i][j]);
    }

    fprintf(fpt,"\n\t");

    for (j=1; j<=sm->numintrvls[i]; ++j) {
      if (sm->nummolecules[i] == 1) {
        fprintf(fpt, " %10.4e", sm->abscoef[i][j]);
      }
      else {
        fprintf(fpt, " %10.4e", sm->abscoef[i][sm->nummolecules[i]*(j-1)+1]);
        for (m=2; m<=sm->nummolecules[i]; m++) {
          fprintf(fpt, " %10.4e", sm->abscoef[i][sm->nummolecules[i]*(j-1)+m]);
        }
      }
    }

    fprintf(fpt,"\n\n");

    /*------------------------------------------------------------------------*/

  }

  /*--------------------------------------------------------------------------*/
  /* Close check.spectralmodel.pollack                                        */
  /*--------------------------------------------------------------------------*/

  fclose(fpt);

  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
