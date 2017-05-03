/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "../include/numrec.nrutil.h"
#include "../include/PhotonSpace.h"
#include "../include/Constituents.h"

/******************************************************************************/

void check_photon_space_final(ps)
  PhotonSpace
    *ps;
{
  int
    i,
    j;

  FILE
    *fpta,
    *fopen();

  /*--------------------------------------------------------------------------*/
  /* Open check.photon.space.final                                            */
  /*--------------------------------------------------------------------------*/

  if ((fpta=fopen("../../output/check.photon.space.final","w+"))==NULL) {
    printf("Cannot create check.photon.space.final for some reason - EXITING!\n");
    exit(1);
  }

  /*--------------------------------------------------------------------------*/
  /* Print the total number of bounding planes in the simulation.             */
  /*--------------------------------------------------------------------------*/

  fprintf(fpta, "[gridnumf] %5d\n\n", ps->gridnumf);

  /*--------------------------------------------------------------------------*/
  /* Print the plane locations and their corresponding rtflag value.          */
  /*--------------------------------------------------------------------------*/

  fprintf(fpta, "[index][gridf][rtflagf]\n");

  for (i=0; i<=ps->gridnumf; i++) {
    fprintf(fpta, "%5d %12.4e %5d\n", i, ps->gridf[i], ps->rtflagf[i]);
  }

  /*--------------------------------------------------------------------------*/
  /* Print the plane indices at which rtflagf is 1.                           */
  /*--------------------------------------------------------------------------*/

  fprintf(fpta, "\n[rt_number] %5d\n", ps->rt_number);

  for (i=1; i<=ps->rt_number; i++) {
    fprintf(fpta, " %5d\n", ps->rtflagfwhich[i]);
  }

  /*--------------------------------------------------------------------------*/
  /* CLOSE check.photon.space.final                                           */
  /*--------------------------------------------------------------------------*/

  fclose(fpta);

  /*--------------------------------------------------------------------------*/
  /* DONE so exit from the routine.                                           */
  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
