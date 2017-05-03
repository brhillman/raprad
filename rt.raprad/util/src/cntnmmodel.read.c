/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../include/PhotonPartition.h"
#include "../include/SpectralModel.h"

/******************************************************************************/

void cntnmmodel_read(char *filename, PhotonPartition *pp, SpectralModel *sm)
{
  FILE *fpt, *fopen();
  char c;
  int i, j, k, bandnumber, numbands, useband;
  double dummyd;

  /*--------------------------------------------------------------------------*/
  /* If no band models are being used there is no need to go further so       */
  /* return to the main function.                                             */
  /*--------------------------------------------------------------------------*/

  if (sm->numwavelngths == 0) {
    return;
  }

  /*--------------------------------------------------------------------------*/
  /* Opne configfileCntnm                                                     */
  /*--------------------------------------------------------------------------*/

  if((fpt=fopen(filename,"r")) == NULL) {
    printf ("Cannot find input file: %s\n", filename);
    exit(0);
  }

  /*--------------------------------------------------------------------------*/
  /* Allocate space to hold the continuum coefficients for self and foreign   */
  /* broadening.                                                              */
  /*--------------------------------------------------------------------------*/

  sm->cntnmscoef = (double **) malloc(sm->numwavelngths*sizeof(double *));
  sm->cntnmscoef -= 1;
  sm->cntnmfcoef = (double **) malloc(sm->numwavelngths*sizeof(double *));
  sm->cntnmfcoef -= 1;

  /*--------------------------------------------------------------------------*/
  /* Skip down to the number of bands indicator and read it in.  Then skip    */
  /* down to the first record to prepare for extraction of the water vapor    */
  /* continuum coefficients.                                                  */
  /*--------------------------------------------------------------------------*/

  while((c=getc(fpt))!='@');
  while((c=getc(fpt))!='\n');

  fscanf(fpt,"%d", &numbands);

  while((c=getc(fpt))!='@');
  while((c=getc(fpt))!='\n');

  /*--------------------------------------------------------------------------*/
  /* Selectively extract those bands specified in pp.  k is a counter that is */
  /* advanced after each band is tested and hence is an index to the current  */
  /* band in the spectral model that is currently being tested for extraction.*/
  /*--------------------------------------------------------------------------*/

  i = 1;

  for (j=1;j<=numbands;++j) {

    /*------------------------------------------------------------------------*/
    /* Simple test to make sure that you have not hit the end of the file too */
    /* soon.                                                                  */
    /*------------------------------------------------------------------------*/

    if ((c=getc(fpt))=='@') {
      printf("You have specified more wavelengths than are in the file\n");
      printf("Exiting ...\n");
      exit(1);
    }
    else {
      ungetc(c,fpt);
    }

    /*------------------------------------------------------------------------*/
    /* If a band is not being used, skip the rest of the band information.    */
    /*------------------------------------------------------------------------*/

    useband = 0;
    for (k=1; k<=pp->numintervals;k++) {
      if (pp->interval[i] == j) {
        useband = 1;
        break;
      }
    }

    if (!useband) {
      while((c=getc(fpt))!='\n');
      continue;
    }

    /*------------------------------------------------------------------------*/
    /* If the band is being used, read in the continuum band information.     */
    /*------------------------------------------------------------------------*/

    fscanf(fpt, "%d",&bandnumber);

    sm->cntnmscoef[i] = (double *) malloc(2*sizeof(double));
    sm->cntnmfcoef[i] = (double *) malloc(2*sizeof(double));

    fscanf(fpt,"%lf",&dummyd);
    fscanf(fpt,"%lf",&dummyd);
    fscanf(fpt,"%lf",&dummyd);

    fscanf(fpt,"%lf",&sm->cntnmscoef[i][0]);
    fscanf(fpt,"%lf",&sm->cntnmscoef[i][1]);
    fscanf(fpt,"%lf",&sm->cntnmfcoef[i][0]);
    fscanf(fpt,"%lf",&sm->cntnmfcoef[i][1]);

    /*------------------------------------------------------------------------*/
    /* Since we are done loading in this band, align the file to the next one.*/
    /*------------------------------------------------------------------------*/

    while((c=getc(fpt))!='\n');

    /*------------------------------------------------------------------------*/
    /* We have loaded band i so go to band i+1.  If i+1 is greater than       */
    /* sm->numwavelngths, then we have all of the bands so break.               */
    /*------------------------------------------------------------------------*/

    ++i;

    if (sm->numwavelngths < i) {
      break;
    }

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
