/******************************************************************************/
/******************************************************************************/

#include <stdio.h>

#include "../include/numrec.nrutil.h"
#include "../include/PhotonPartition.h"

/******************************************************************************/

void
check_photon_partition(pp)
  PhotonPartition
    *pp;
{
  int
    i,
    j;

  FILE
    *fpta,
    *fopen();

  /*--------------------------------------------------------------------------*/
  /* Open check.photon.partition                                              */
  /*--------------------------------------------------------------------------*/

  if ((fpta=fopen("../../output/check.photon.partition","w+"))==NULL) {
    printf("check configuration file cannot be opened for some reason - EXITING!\n");
    exit(1);
  }

  /*--------------------------------------------------------------------------*/
  /* Print the number of spectral intervals to be used in this run.           */
  /*--------------------------------------------------------------------------*/

  fprintf(fpta, "[numintervals] %d\n\n", pp->numintervals);

  fprintf(fpta, "[index][wavelength][interval][runsperinterval][addscatterflag][bbradflag][numrefl]\n");

  /*--------------------------------------------------------------------------*/
  /* Print the information for each band in the model.                        */
  /*--------------------------------------------------------------------------*/

  for (i=1; i<=pp->numintervals; ++i) {

    /*------------------------------------------------------------------------*/
    /* Print the channel or band interval index which is the same as i.       */
    /*------------------------------------------------------------------------*/

    fprintf(fpta, "%5d", i);

    /*------------------------------------------------------------------------*/
    /* Print in the channel or band wavelength.                               */
    /*------------------------------------------------------------------------*/

    fprintf(fpta, " %12.4e", pp->wavelength[i]);

    /*------------------------------------------------------------------------*/
    /* Print the channel or spectral band index number.                       */
    /*------------------------------------------------------------------------*/

    fprintf(fpta, " %5d", pp->interval[i]);

    /*------------------------------------------------------------------------*/
    /* Print the number of spectral sub-intervals for this band.              */
    /*------------------------------------------------------------------------*/

    fprintf(fpta, " %5d", pp->runsperinterval[i]);

    /*------------------------------------------------------------------------*/
    /* Print the type of scattering calculations to be performed.             */
    /*------------------------------------------------------------------------*/

    fprintf(fpta, " %5d", pp->addscatterflag[i]);

    /*------------------------------------------------------------------------*/
    /* Print whether (1) or not (0) to include layer emission.                */
    /*------------------------------------------------------------------------*/

    fprintf(fpta, " %5d", pp->bbradflag[i]);

    /*------------------------------------------------------------------------*/
    /* Print in the number of reflections between the atmosphere and surface.  */
    /*------------------------------------------------------------------------*/

    fprintf(fpta, " %5d\n", pp->numrefl[i]);

    /*------------------------------------------------------------------------*/
    /* Done printing in a single spectral band.                               */
    /*------------------------------------------------------------------------*/

  }

  /*--------------------------------------------------------------------------*/
  /* Close check.photon.partition                                             */
  /*--------------------------------------------------------------------------*/

  fclose(fpta);

  /*--------------------------------------------------------------------------*/
  /* DONE checking the photon partition so exiting the routine.               */
  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
