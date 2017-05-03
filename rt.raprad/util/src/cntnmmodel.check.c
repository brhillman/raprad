/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../include/PhotonPartition.h"
#include "../include/SpectralModel.h"

/******************************************************************************/

void check_cntnmmodel(PhotonPartition *pp, SpectralModel *sm)
{
  FILE *fpt, *fopen();
  int i;

  /*--------------------------------------------------------------------------*/
  /* Open check.cntnmmodel                                                    */
  /*--------------------------------------------------------------------------*/

  if((fpt=fopen("../output/check.cntnmmodel", "w+")) == NULL) {
    printf ("Cannot create check.cntnmmodel - Exiting!\n");
    exit(0);
  }

  /*--------------------------------------------------------------------------*/
  /* Print the number of bands used in the current simulation.                */
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

  fprintf(fpt, "\n\t[cntnmscoef0][cntnmscoef1][cntnmfcoef0][cntnmfcoef1]\n\n");

  for (i=1; i<=sm->numwavelngths; ++i) {

    /*------------------------------------------------------------------------*/
    /* If the band is being used, print the continuum band information.       */
    /*------------------------------------------------------------------------*/

    fprintf(fpt, "\n\t");

    fprintf(fpt, " %12.4e", sm->cntnmscoef[i][0]);
    fprintf(fpt, " %12.4e", sm->cntnmscoef[i][1]);
    fprintf(fpt, " %12.4e", sm->cntnmfcoef[i][0]);
    fprintf(fpt, " %12.4e", sm->cntnmfcoef[i][1]);

    /*------------------------------------------------------------------------*/

  }

  /*--------------------------------------------------------------------------*/
  /* Close configfileCntnm                                                    */
  /*--------------------------------------------------------------------------*/

  fclose(fpt);

  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
