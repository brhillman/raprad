/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../include/numrec.nrutil.h"

#include "../include/PhotonPartition.h"
#include "../include/SpectralModel.h"

/******************************************************************************/

int check_spectralmodel_kato(pp, sm)
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
    k,
    m;

  /*--------------------------------------------------------------------------*/
  /* Open check.spectralmodel.pollack                                         */
  /*--------------------------------------------------------------------------*/

  if((fpt=fopen("../../output/check.spectralmodel.kato","w+")) == NULL) {
    printf ("Cannot create file check.spectralmodel.kato - Exiting!\n");
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

    for (j=1; j<=sm->nummolecules[i]; j++) {
      fprintf(fpt, "%5d", sm->nquadpts[i][j]);
    }

    fprintf(fpt, "%5d", sm->numintrvls[i]);

    fprintf(fpt,"\n\t");

    for (j=1; j<=sm->nummolecules[i]; j++) {

      for (k=1; k<=sm->nquadpts[i][j]; k++) {
        fprintf(fpt, "%10.7f", sm->quadwts[i][j][k]);
      }
      fprintf(fpt,"\n\t");

      if (sm->nquadpts[i][j] == 1) {
        fprintf(fpt, "%10.7f\n", sm->abscoef[i][j]);
      }
      else {
        fprintf(fpt, "%s\n", sm->abscoeffile[i][j]);
      }

    }

    fprintf(fpt,"\n\t");

    for (j=1; j<=sm->numintrvls[i]; ++j) {
      fprintf(fpt," %12.4e", sm->alpha[i][j]);
    }

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
