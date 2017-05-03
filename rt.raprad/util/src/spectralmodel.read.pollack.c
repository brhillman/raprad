/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../include/numrec.nrutil.h"

#include "../include/PhotonPartition.h"
#include "../include/SpectralModel.h"

/******************************************************************************/

void spectralmodel_read_pollack(char *filename, PhotonPartition *pp, SpectralModel *sm)
{
  FILE *fpt, *fopen();
  char c;
  int i, j, k, m;
  int bandnumber, numbands, numintervals, nummolecules, useband;
   

  /*--------------------------------------------------------------------------*/
  /* Count the number of gaseous absorption bands that are used in the        */
  /* simulation and place the result into sm->numwavelngths.                  */
  /*--------------------------------------------------------------------------*/

  sm->numwavelngths = 0;
  for (i=1; i<=pp->numintervals; i++) {
    if (pp->interval[i] != -1) {
      ++sm->numwavelngths;
    }
  }

  /*--------------------------------------------------------------------------*/
  /* If no band models are being used there is no need to go further so       */
  /* return to the main function.                                             */
  /*--------------------------------------------------------------------------*/

  if (sm->numwavelngths == 0) {
    return;
  }

  /*--------------------------------------------------------------------------*/
  /* Open configfileSpectralModel                                             */
  /*--------------------------------------------------------------------------*/

  if((fpt=fopen(filename,"r")) == NULL) {
    printf ("Cannot find input file: %s\n", filename);
    exit(0);
  }

  /*--------------------------------------------------------------------------*/
  /* Skip past the header stuff to the number number of bands in the model.   */
  /*--------------------------------------------------------------------------*/

  while((c=getc(fpt))!='@');
  while((c=getc(fpt))!='\n');

  /*--------------------------------------------------------------------------*/
  /* Allocate the space needed to hold all of the spectral model information. */
  /*--------------------------------------------------------------------------*/

  /* Name of the constituent or constituents if there are more than one */
  /* active in this spectral region. */

  sm->gas          = (char  **)  malloc(sm->numwavelngths*sizeof(char *));
  sm->gas -= 1;

  /* Number of active constituents in this spectral region */

  sm->nummolecules = ivector(1, sm->numwavelngths);

  /* Use (1) or do not use (0) the continuum in this interval for this run. */

  sm->cntnm = ivector(1, sm->numwavelngths);

  /* Center wavelength of this spectral band */

  sm->wavelngth = dvector(1, sm->numwavelngths);

  /* Beginning and ending wavelengths of the spectral band */

  sm->bandwidth = dmatrix(1, sm->numwavelngths, 1, 2);

  /* Total flux in W m^{-2} across the band */

  sm->solarinsol = dvector(1, sm->numwavelngths);

  /* Number of sub-intervals in the band */

  sm->numintrvls = ivector(1, sm->numwavelngths);

  /* Pressure effect on the absorption coefficient; there is one for each */
  /* constituent comprising the band. */

  sm->coverd     = (double **) malloc(sm->numwavelngths*sizeof(double *));
  sm->coverd -= 1;

  /* Sub-interval weighting coefficient */

  sm->alpha      = (double **) malloc(sm->numwavelngths*sizeof(double *));
  sm->alpha -= 1;

  /* Sub-interval weighting coefficients with multiple entries per subinterval */

  sm->alphad  = (double ***) malloc((sm->numwavelngths+1)*sizeof(double **));
  sm->alphau  = (double ***) malloc((sm->numwavelngths+1)*sizeof(double **));

  /* Sub-interval absorption coefficient */

  sm->abscoef    = (double **) malloc(sm->numwavelngths*sizeof(double *));
  sm->abscoef -= 1;

  /* Sub-interval absorption coefficient data */

  sm->abscoefdata  = (double ***) malloc((sm->numwavelngths+1)*sizeof(double **));

  /*--------------------------------------------------------------------------*/
  /* Align the input file to the first spectral band in the model.            */
  /*--------------------------------------------------------------------------*/

  fscanf(fpt, "%d",&numbands);
  while((c=getc(fpt))!='@');
  while((c=getc(fpt))!='\n');
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
    /* Decrement i as it is a counter for the number of spectral intervals    */
    /* used and the continue statement increments it when it should not be    */
    /* incremented.  Also, increment k since we are moving on to the next     */
    /* band.                                                                  */
    /*------------------------------------------------------------------------*/

    useband = 0;
    for (k=1;k<=pp->numintervals;k++) {
      if (pp->interval[k] == j) {
        useband = 1;
        break;
      }
    }

    if (!useband) {
      while((c=getc(fpt))!='~');
      while((c=getc(fpt))!='\n');
      continue;
    }

    /*------------------------------------------------------------------------*/
    /* If the band is being used, read in the band information.               */
    /*------------------------------------------------------------------------*/

    fscanf(fpt, "%d",&bandnumber);

    sm->gas[i] = (char  *) malloc(32*sizeof(char));
    fscanf(fpt, "%s",sm->gas[i]);

    fscanf(fpt, "%d",&sm->nummolecules[i]);

    fscanf(fpt, "%d",&sm->cntnm[i]);

    fscanf(fpt,"%lf",&sm->bandwidth[i][0]);
    fscanf(fpt,"%lf",&sm->bandwidth[i][1]);

    fscanf(fpt,"%lf",&sm->wavelngth[i]);
    fscanf(fpt,"%lf",&sm->solarinsol[i]);

    fscanf(fpt, "%d",&sm->numintrvls[i]);

    nummolecules = sm->nummolecules[i];
    sm->coverd[i] = (double *) malloc(nummolecules*sizeof(double));
    sm->coverd[i] -= 1;
    for (k=1; k<=nummolecules; ++k) {
      fscanf(fpt,"%lf",&sm->coverd[i][k]);
    }

    numintervals = sm->numintrvls[i];

    sm->alpha[i] = (double *) malloc(numintervals*sizeof(double));
    sm->alpha[i] -= 1;
    for (k=1; k<=numintervals; ++k) {
      fscanf(fpt,"%lf",&sm->alpha[i][k]);
    }

    sm->abscoef[i] = (double*) malloc(nummolecules*numintervals*sizeof(double));
    sm->abscoef[i] -= 1;
    for (k=1; k<=numintervals; ++k) {
      if (nummolecules == 1) {
        fscanf(fpt,"%lf",&sm->abscoef[i][k]);
      }
      else {
        while((c=getc(fpt))!='(');
        fscanf(fpt,"%lf",&sm->abscoef[i][nummolecules*(k-1)+1]);
        for (m=2; m<=nummolecules; m++) {
          while((c=getc(fpt))!=',');
          fscanf(fpt,"%lf",&sm->abscoef[i][nummolecules*(k-1)+m]);
        }
      }
    }

    /*------------------------------------------------------------------------*/
    /* Since we are done loading in this band, align the file to the next one.*/
    /*------------------------------------------------------------------------*/

    while((c=getc(fpt))!='~');
    while((c=getc(fpt))!='\n');

    /*------------------------------------------------------------------------*/
    /* We have loaded band i so go to band i+1.  If i+1 is greater than       */
    /* sm->numwavelngths, then we have all of the bands so break.             */
    /*------------------------------------------------------------------------*/

    ++i;

    if (sm->numwavelngths < i) {
      break;
    }

    /*------------------------------------------------------------------------*/

  }

  /*--------------------------------------------------------------------------*/
  /* Close configfileSpectralModel                                            */
  /*--------------------------------------------------------------------------*/

  fclose(fpt);

  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
