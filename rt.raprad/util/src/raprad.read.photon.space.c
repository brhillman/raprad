/*****************************************************************************/
/*****************************************************************************/

#include <stdio.h>

#include "../include/PhotonSpace.h"

/*****************************************************************************/

void
read_photon_space(configfilename, ps)
char
  *configfilename;
PhotonSpace
  *ps;
{
  int
    index,
    j,
    number;

  FILE
    *fptr,
    *fopen();

  /*--------------------------------------------------------------------------*/
  /* Open configfilename                                                      */
  /*--------------------------------------------------------------------------*/

  if ((fptr=fopen(configfilename,"r"))==NULL) {
    printf("the configuration file cannot be opened for some reason - EXIT!\n");
    exit(1);
  }

  /*--------------------------------------------------------------------------*/
  /* Read in the top and bottom heights (in m) of the modelled atmosphere.    */
  /*--------------------------------------------------------------------------*/

  if (setup(fptr)) { fscanf(fptr, "%lf", &ps->zm); }
  if (setup(fptr)) { fscanf(fptr, "%lf", &ps->zp); }

  /*--------------------------------------------------------------------------*/
  /* Read in the total number of layers in the atmosphere.                    */
  /*--------------------------------------------------------------------------*/

  if (setup(fptr)) { fscanf(fptr, "%d", &ps->gridnum); }

  /*--------------------------------------------------------------------------*/
  /* Allocate memory to hold in the height information for the                */
  /* ps->gridnum+1 planes that bound the ps->gridnum layers.  Notice that     */
  /* this array is not 0 offset; we adapted this convention so that           */
  /* indices 1 to ps->gridnum correspond to the top heights of the            */
  /* ps->gridnum layers.  Also, allocate memory to hold the value of the flag */
  /* as to whether or not to build r and t matrices for the atmosphere both   */
  /* above and below this height.                                             */
  /*--------------------------------------------------------------------------*/

  ps->grid   = (double *) malloc((ps->gridnum+1)*sizeof(double));
  ps->rtflag =    (int *) malloc((ps->gridnum+1)*sizeof(   int));

  /*--------------------------------------------------------------------------*/
  /* Position the pointer to the file to the plane information.               */
  /*--------------------------------------------------------------------------*/

  position(fptr);

  /*--------------------------------------------------------------------------*/
  /* Read in the heights of the ps->gridnum+1 planes that bound the           */
  /* ps->gridnum layers, as well as the flags indicating whether or not to    */
  /* build r and t matrices at the specified planes.                          */
  /*--------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------*/
  /* Initialize the r/t flag counter.                                         */
  /*--------------------------------------------------------------------------*/

  ps->rt_number = 0;

  /*--------------------------------------------------------------------------*/
  /* Read in the photon space data.                                           */
  /*--------------------------------------------------------------------------*/

  for (j=0; j<=ps->gridnum; j++) {

    /*------------------------------------------------------------------------*/
    /* Read in the plane height and r/t matrix flag value.                    */
    /*------------------------------------------------------------------------*/

    fscanf(fptr, "%d",  &number);
    fscanf(fptr, "%lf", &ps->grid[j]);
    fscanf(fptr, "%d", &ps->rtflag[j]);

    /*------------------------------------------------------------------------*/
    /* Count the number of planes that are flagged for r and t matrices.      */
    /*------------------------------------------------------------------------*/

    ps->rt_number += ps->rtflag[j];

    /*------------------------------------------------------------------------*/
      
  }

  /*--------------------------------------------------------------------------*/
  /* DONE reading in the photon space data.                                   */
  /*--------------------------------------------------------------------------*/

}

/*****************************************************************************/
/*****************************************************************************/
