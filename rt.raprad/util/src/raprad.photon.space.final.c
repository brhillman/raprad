/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "../include/numrec.nrutil.h"
#include "../include/PhotonSpace.h"
#include "../include/Constituents.h"

/*----------------------------------------------------------------------------*/

#define DELTA 0.05

/******************************************************************************/
extern void piksrt(int, double []);

void photon_space_final(PhotonSpace *ps, Constituents *c)
{
  int i, j, n, nn, planecount, planes_max, rtonflag;
  double *planesi, *planesf;

  /*--------------------------------------------------------------------------*/
  /* Allocate space to hold all of the heights of all possible constituent    */
  /* boundaries.                                                              */
  /*--------------------------------------------------------------------------*/

  planesi = (double *) malloc((2*ps->cnumber+1)*sizeof(double));
  planesi -= 1;

  /*--------------------------------------------------------------------------*/
  /* Allocate space to hold all of the distinct constituent boundary heights  */
  /* in increasing order of height.                                           */
  /*--------------------------------------------------------------------------*/

  planesf = (double *) malloc((2*ps->cnumber+1)*sizeof(double));
  planesf -= 1;

  /*--------------------------------------------------------------------------*/
  /* Load in all of the constituent boundary heights, including the -1s, into */
  /* planesi.                                                                 */
  /*--------------------------------------------------------------------------*/

  n = 0;
  for (i=1; i<=ps->cnumber; i++) {
    planesi[++n] = c[i].base_height;
    planesi[++n] = c[i].top_height;
  }

  /*--------------------------------------------------------------------------*/
  /* Sort the constituent boundary heights into monotonically increasing      */
  /* order.                                                                   */
  /*--------------------------------------------------------------------------*/

  piksrt(n, planesi);

  /*--------------------------------------------------------------------------*/
  /* If the same height is counted more than once, count it only once.        */
  /*--------------------------------------------------------------------------*/

  nn = 0;
  planesf[++nn] = planesi[1];
  for (i=1; i<n; i++) {
    if (planesi[i]<planesi[i+1]) {
      planesf[++nn] = planesi[i+1];
    }
  }

  /*--------------------------------------------------------------------------*/
  /* Get the total number of possible planes, which is the sum of the user    */
  /* specified planes and the distinct constituent boundaries.                */
  /*--------------------------------------------------------------------------*/

  planes_max = nn + (ps->gridnum+1);

  /*--------------------------------------------------------------------------*/
  /* Allocate space for the final grid.                                       */
  /*--------------------------------------------------------------------------*/

  ps->gridf = (double *) malloc(planes_max*sizeof(double));

  /*--------------------------------------------------------------------------*/
  /* Combine the user requested planes and the distinct constituent boundary  */
  /* planes.                                                                  */
  /*--------------------------------------------------------------------------*/

  planecount = 0;

  for (i=0; i<ps->gridnum; i++) {

    ps->gridf[planecount++] = ps->grid[i];

    for (j=1; j<=nn; j++) {

      if ( ((ps->grid[i]+DELTA)<planesf[j]) &&
           (planesf[j]<(ps->grid[i+1]-DELTA))  ) {
        ps->gridf[planecount++] = planesf[j];
      }

    }

  }

  ps->gridf[planecount++] = ps->grid[ps->gridnum];
  ps->gridnumf = planecount - 1;

  /*--------------------------------------------------------------------------*/
  /* Update the r/t flag value array so that it is the same size as ps->gridf.*/
  /*--------------------------------------------------------------------------*/

  ps->rtflagf = (int *) malloc((ps->gridnumf+1)*sizeof(int));

  for (i=0; i<=ps->gridnumf; i++) {

    rtonflag = 1;
    for (j=0; j<=ps->gridnum; j++) {
      if (ps->gridf[i] == ps->grid[j]) {
        if (ps->rtflag[j]) { ps->rtflagf[i] = 1; }
        else               { ps->rtflagf[i] = 0; }
        rtonflag = 0;
        break;
      }
    }

    if (rtonflag) {
      ps->rtflagf[i] = 1;
      ++(ps->rt_number);
    }

  }

  /*--------------------------------------------------------------------------*/
  /* It is important to know which planes are flagged as well.  We index the  */
  /* planes from the top of the atmosphere down to the surface in keeping     */
  /* with indices in the optical depth arrays, etc...                         */
  /*--------------------------------------------------------------------------*/

  ps->rtflagfwhich = (int *) malloc((ps->rt_number+1)*sizeof(int));

  j = 1;
  for (i=ps->gridnumf; i>=0; i--) {
    if (ps->rtflagf[i]) {
      ps->rtflagfwhich[j++] = i;
    }
  }

  if ((j-1) != ps->rt_number) {
    printf("There is an inconsistency in the planes that are flagged for r\n");
    printf("   and t matrices.\n");
    exit(0);
  }

  /*--------------------------------------------------------------------------*/
/*
  for (i=0; i<=ps->gridnumf; i++) {
    printf("Plane %5d at height %6.2f\n", i, ps->gridf[i]);
  }
*/
  /*--------------------------------------------------------------------------*/
  /* Compute the physical thickness of each layer.                            */
  /*--------------------------------------------------------------------------*/

  ps->dz = dvector(1, ps->gridnumf);

  for (i=1; i<=ps->gridnumf; i++) {
    ps->dz[i] = ps->gridf[i] - ps->gridf[i-1];
  }

  /*--------------------------------------------------------------------------*/
  /* Free the space allocated for planesi.                                    */
  /*--------------------------------------------------------------------------*/

  planesi += 1;
  free((double *)planesi);

  /*--------------------------------------------------------------------------*/
  /* Free the space allocated for planesf.                                    */
  /*--------------------------------------------------------------------------*/

  planesf += 1;
  free((double *)planesf);

  /*--------------------------------------------------------------------------*/
  /* DONE so exit from the routine.                                           */
  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
