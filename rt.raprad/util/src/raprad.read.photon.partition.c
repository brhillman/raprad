/******************************************************************************/
/******************************************************************************/

#include <stdio.h>

#include "../include/numrec.nrutil.h"
#include "../include/PhotonPartition.h"

/******************************************************************************/

void
read_photon_partition(configfilePhotonPartition, pp)
  char
    *configfilePhotonPartition;
  PhotonPartition
    *pp;
{
  int
    i,
    j,
    dummyi,
    smindex;

  FILE
    *fpta,
    *fopen();

  char
    c;

  /*--------------------------------------------------------------------------*/
  /* Open  configfilePhotonPartition                                          */
  /*--------------------------------------------------------------------------*/

  if ((fpta=fopen(configfilePhotonPartition,"r"))==NULL) {
    printf("configuration file cannot be opened for some reason - EXITING!\n");
    exit(1);
  }

  /*--------------------------------------------------------------------------*/
  /* Read in the number of spectral intervals to be used in this run.         */
  /*--------------------------------------------------------------------------*/

  if (setup(fpta)) { fscanf(fpta, "%d", &pp->numintervals); }

  /*--------------------------------------------------------------------------*/
  /* Allocate space to hold the spectral model information.                   */
  /*--------------------------------------------------------------------------*/

  /* Index value of the band as labelled in the spectral model file. */

  pp->wavelength      =  vector(1, pp->numintervals);
  pp->interval        = ivector(1, pp->numintervals);
  pp->smindex         = ivector(1, pp->numintervals);
  pp->runsperinterval = ivector(1, pp->numintervals);
  pp->addscatterflag  = ivector(1, pp->numintervals);
  pp->bbradflag       = ivector(1, pp->numintervals);
  pp->numrefl         = ivector(1, pp->numintervals);

  /*--------------------------------------------------------------------------*/
  /* Position the file pointer to begin reading in the band model information.*/
  /*--------------------------------------------------------------------------*/

  position(fpta);

  /*--------------------------------------------------------------------------*/
  /* Begin reading in the information for each band in the model.             */
  /*--------------------------------------------------------------------------*/

  smindex = 1;

  for (i=1; i<=pp->numintervals; ++i) {

    /*------------------------------------------------------------------------*/
    /* Read in the channel or band interval index which is the same as i.     */
    /*------------------------------------------------------------------------*/

    fscanf(fpta, "%d", &dummyi);

    /*------------------------------------------------------------------------*/
    /* Read in the channel or band wavelength.                                */
    /*------------------------------------------------------------------------*/

    fscanf(fpta, "%f", &pp->wavelength[i]);

    /*------------------------------------------------------------------------*/
    /* Read in the channel or spectral band index number                      */
    /*------------------------------------------------------------------------*/

    fscanf(fpta, "%d", &pp->interval[i]);

    if (pp->interval[i] != -1) {
      pp->smindex[i] = smindex++;
    }
    else {
      pp->smindex[i] = -1;
    }

    /*------------------------------------------------------------------------*/
    /* Read in the number of spectral sub-intervals for this band.            */
    /*------------------------------------------------------------------------*/

    fscanf(fpta, "%d", &pp->runsperinterval[i]);

    /*------------------------------------------------------------------------*/
    /* Read in the type of scattering calculations to be performed.           */
    /*------------------------------------------------------------------------*/

    fscanf(fpta, "%d", &pp->addscatterflag[i]);

    /*------------------------------------------------------------------------*/
    /* Read in whether (1) or not (0) to include layer emission.              */
    /*------------------------------------------------------------------------*/

    fscanf(fpta, "%d", &pp->bbradflag[i]);

    /*------------------------------------------------------------------------*/
    /* Read in the number of reflections between the atmosphere and surface.  */
    /*------------------------------------------------------------------------*/

    fscanf(fpta, "%d", &pp->numrefl[i]);

    /*------------------------------------------------------------------------*/
    /* Done reading in a single spectral band.                                */
    /*------------------------------------------------------------------------*/

  }

  /*--------------------------------------------------------------------------*/
  /* Close configfilePhotonPartition                                          */
  /*--------------------------------------------------------------------------*/

  fclose(fpta);

  /*--------------------------------------------------------------------------*/
  /* DONE reading the photon partition so exiting the routine.                */
  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
