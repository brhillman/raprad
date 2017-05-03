/****************************************************************/
/****************************************************************/

#include <stdio.h>

#include "../include/numrec.nrutil.h"
#include "../include/PhotonPartition.h"
#include "../include/Angles.h"
#include "../include/Brdf.h"

/****************************************************************/

int
read_brdf_selection(brdfpropertiesfile, pp, b)
char
  *brdfpropertiesfile;
PhotonPartition
  *pp;
Brdf
  *b;
{
  int
    c,
    i,
    j;

  FILE
    *finput,
    *fopen();

  /*-------------------------------------------------------------*/

  finput = fopen(brdfpropertiesfile,"r");

  /*-------------------------------------------------------------*/

  if (setup(finput)) { fscanf(finput, "%d", &b->total_number); }
  if (setup(finput)) { fscanf(finput, "%d", &b->used_number); }

  /*-------------------------------------------------------------*/

  b->whichscene = ivector(1, b->used_number);
  b->channel    = ivector(1, b->used_number);
  b->symmetry	= ivector(1, b->used_number);

  /*-------------------------------------------------------------*/

  if (setup(finput)) {
    for (i=1; i<=b->used_number; i++) {
      fscanf(finput, "%d", &b->whichscene[i]);
    }
  }

  /*-------------------------------------------------------------*/

  if (setup(finput)) {
    for (i=1; i<=b->used_number; i++) {
      fscanf(finput, "%d", &b->channel[i]);
    }
  }

  /*-------------------------------------------------------------*/

  if (setup(finput)) {
    for (i=1; i<=b->used_number; i++) {
      fscanf(finput, "%d", &b->symmetry[i]);
    }
  }

  /*-------------------------------------------------------------*/

  if (setup(finput)) { fscanf(finput, "%d", &b->normalize); }

  /*-------------------------------------------------------------*/

  b->index     = ivector(1, pp->numintervals);
  b->nbrdfruns = ivector(1, pp->numintervals);

  b->brdfindex = (int **)   malloc((pp->numintervals+1) * sizeof(int *));
  b->srftemp   = (float **) malloc((pp->numintervals+1) * sizeof(float *));
  b->albedo    = (float **) malloc((pp->numintervals+1) * sizeof(float *));

  while ((c=getc(finput)) != ':');

  for (i=1; i<=pp->numintervals; ++i) {

    fscanf(finput, "%d", &b->index[i]);

    fscanf(finput, "%d", &b->nbrdfruns[i]);

    b->brdfindex[i] = (int *) malloc((b->nbrdfruns[i]+1) * sizeof(int));

    for (j=1; j<=b->nbrdfruns[i]; j++) {
      fscanf(finput, "%d", &b->brdfindex[i][j]);
    }

    b->srftemp[i] = (float *) malloc((b->nbrdfruns[i]+1) * sizeof(float));

    for (j=1; j<=b->nbrdfruns[i]; j++) {
      fscanf(finput, "%f", &b->srftemp[i][j]);
    }

    b->albedo[i] = (float *) malloc((b->nbrdfruns[i]+1) * sizeof(float));

    for (j=1; j<=b->nbrdfruns[i]; j++) {
      fscanf(finput, "%f", &b->albedo[i][j]);
    }

  }

  /*-------------------------------------------------------------*/

  fclose(finput);

  /*-------------------------------------------------------------*/

}

/****************************************************************/
/****************************************************************/
