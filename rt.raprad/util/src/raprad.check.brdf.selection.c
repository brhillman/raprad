/****************************************************************/
/****************************************************************/

#include <stdio.h>

#include "../include/numrec.nrutil.h"
#include "../include/PhotonPartition.h"
#include "../include/Angles.h"
#include "../include/Brdf.h"

/****************************************************************/

void
check_brdf_selection(pp,b)
PhotonPartition
  *pp;
Brdf
  *b;
{
  int
    i,
    j;

  FILE
    *finput,
    *fopen();

  /*-------------------------------------------------------------*/
  /* Open check.brdf.selection                                   */
  /*-------------------------------------------------------------*/

  if((finput=fopen("../../output/check.brdf.selection","w+")) == NULL) {
    printf ("Cannot create file check.brdf.selection - Exiting!\n");
    exit(0);
  }

  /*-------------------------------------------------------------*/

  fprintf(finput, "[total_number] %5d\n", b->total_number);
  fprintf(finput, "[used_number] %5d\n\n", b->used_number);

  fprintf(finput, "[whichscene] ");
  for (i=1; i<=b->used_number; i++) {
    fprintf(finput, " %5d", b->whichscene[i]);
  }

  fprintf(finput, "\n[channel] ");
  for (i=1; i<=b->used_number; i++) {
    fprintf(finput, " %5d", b->channel[i]);
  }

  fprintf(finput, "\n[symmetry] ");
  for (i=1; i<=b->used_number; i++) {
    fprintf(finput, " %5d", b->symmetry[i]);
  }

  fprintf(finput, "\n[normalize] %5d\n\n", b->normalize);

  for (i=1; i<=pp->numintervals; ++i) {

    fprintf(finput, "%10d", b->index[i]);

    fprintf(finput, " %10d\n", b->nbrdfruns[i]);

    for (j=1; j<=b->nbrdfruns[i]; j++) {
      fprintf(finput, " %10d", b->brdfindex[i][j]);
    }
    fprintf(finput, "\n");

    for (j=1; j<=b->nbrdfruns[i]; j++) {
      fprintf(finput, " %10.3f", b->srftemp[i][j]);
    }
    fprintf(finput, "\n");

    for (j=1; j<=b->nbrdfruns[i]; j++) {
      fprintf(finput, " %10.3f", b->albedo[i][j]);
    }
    fprintf(finput, "\n");

  }

  /*-------------------------------------------------------------*/
  /* CLOSE check.geometry.quadrature                             */
  /*-------------------------------------------------------------*/

  fclose(finput);

  /*-------------------------------------------------------------*/
  /* DONE so exit the routine.                                   */
  /*-------------------------------------------------------------*/

}

/****************************************************************/
/****************************************************************/
