/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "../include/numrec.nrutil.h"
#include "../include/Constituents.h"

/******************************************************************************/

void phasefcns_read_binary_header(char *filename, int pwhich, Constituents *c)
{
  int index, pnum, skip;
  double wavelength, radiusout, radiussig, relrefoutr, relrefoutc;
  FILE *fptr, *fopen();

  /*----------------------------------------------------------------------*/
  /* Open the phase function file.                                        */
  /*----------------------------------------------------------------------*/

  if ((fptr=fopen(filename,"r"))==NULL) {
    printf("could not open the phase functions file - exiting!\n");
    exit(0);
  }

  /*----------------------------------------------------------------------*/
  /* Skip to the first point in the file that contains the number of      */
  /* points in the phase function.                                        */
  /*----------------------------------------------------------------------*/

  fseek(fptr,ftell(fptr)+44,0);

  /*----------------------------------------------------------------------*/
  /* Read in the number of points in the first phase function.  We assume */
  /* that all of the phase functions in this file have the same number of */
  /* points.                                                              */
  /*----------------------------------------------------------------------*/

  fread(&pnum,sizeof(int),1,fptr);

  /*----------------------------------------------------------------------*/
  /* Go back to the 5th byte of the file where the phase function         */
  /* information starts.                                                  */
  /*----------------------------------------------------------------------*/

  fseek(fptr, 0, 0);

  /*----------------------------------------------------------------------*/
  /* Now skip to the desired phase function using the fact that there are */
  /* 92 bytes in the header information for each phase function.          */
  /*----------------------------------------------------------------------*/

  skip = (pwhich-1)*(96 + 2*pnum*sizeof(double));
  fseek(fptr,ftell(fptr)+skip,0);

  /*----------------------------------------------------------------------*/
  /* Read in the current phase function.                                  */
  /*----------------------------------------------------------------------*/

  fread(&index		 	,sizeof(int)   ,1,fptr);
  fread(&wavelength 		,sizeof(double),1,fptr);
  fread(&radiusout  		,sizeof(double),1,fptr);
  fread(&radiussig  		,sizeof(double),1,fptr);
  fread(&relrefoutr 		,sizeof(double),1,fptr);
  fread(&relrefoutc 		,sizeof(double),1,fptr);
  fread(&c->p_explicit_number 	,sizeof(int)   ,1,fptr);
  fread(&c->p_g			,sizeof(double),1,fptr);
  fread(&c->w0			,sizeof(double),1,fptr);
  fread(&c->sig_ext		,sizeof(double),1,fptr);
  fread(&c->sig_sca		,sizeof(double),1,fptr);
  fread(&c->sig_abs		,sizeof(double),1,fptr);
  fread(&c->p_normalization	,sizeof(double),1,fptr);

  /*----------------------------------------------------------------------*/
  /* DONE reading in the phase function so close down the file.           */
  /*----------------------------------------------------------------------*/

  fclose(fptr);

  /*----------------------------------------------------------------------*/
  /* DONE so exit.                                                        */
  /*----------------------------------------------------------------------*/


}

/******************************************************************************/

void phasefcns_legendre_read_binary(char *filename, int pwhich, Constituents *c)
{
  int index, pnum, skip;
  double wavelength, radiusout, radiussig, relrefoutr, relrefoutc;
  FILE *fptr, *fopen();

  /*----------------------------------------------------------------------*/
  /* Open the phase function file.                                        */
  /*----------------------------------------------------------------------*/

  if ((fptr=fopen(filename,"r"))==NULL) {
    printf("could not open the phase functions file - exiting!\n");
    exit(0);
  }

  /*----------------------------------------------------------------------*/
  /* Skip to the first point in the file that contains the number of      */
  /* points in the phase function.                                        */
  /*----------------------------------------------------------------------*/

  fseek(fptr,ftell(fptr)+44,0);

  /*----------------------------------------------------------------------*/
  /* Read in the number of legendre coefficients in the first phase       */
  /* function.  We assume that all of the phase functions in this file    */
  /* have the same number of legendre coefficients.                       */
  /*----------------------------------------------------------------------*/

  fread(&pnum,sizeof(int),1,fptr);

  /*----------------------------------------------------------------------*/
  /* Go back to the 5th byte of the file where the phase function         */
  /* information starts.                                                  */
  /*----------------------------------------------------------------------*/

  fseek(fptr, 0, 0);

  /*----------------------------------------------------------------------*/
  /* Now skip to the desired phase function using the fact that there are */
  /* 92 bytes in the header information for each phase function.          */
  /*----------------------------------------------------------------------*/

  skip = (pwhich-1)*(96 + pnum*sizeof(double));
  fseek(fptr,ftell(fptr)+skip,0);

  /*----------------------------------------------------------------------*/
  /* Read in the current phase function.                                  */
  /*----------------------------------------------------------------------*/

  fread(&index		 	,sizeof(int)   ,1,fptr);
  fread(&wavelength 		,sizeof(double),1,fptr);
  fread(&radiusout  		,sizeof(double),1,fptr);
  fread(&radiussig  		,sizeof(double),1,fptr);
  fread(&relrefoutr 		,sizeof(double),1,fptr);
  fread(&relrefoutc 		,sizeof(double),1,fptr);
  fread(&c->p_legendre_number 	,sizeof(int)   ,1,fptr);
  fread(&c->p_g			,sizeof(double),1,fptr);
  fread(&c->w0			,sizeof(double),1,fptr);
  fread(&c->sig_ext		,sizeof(double),1,fptr);
  fread(&c->sig_sca		,sizeof(double),1,fptr);
  fread(&c->sig_abs		,sizeof(double),1,fptr);
  fread(&c->p_normalization	,sizeof(double),1,fptr);

  /*----------------------------------------------------------------------*/
  /* Read in the legendre polynomial coefficients.                        */
  /*----------------------------------------------------------------------*/

  /* We do this allocation only once. */

  if (pwhich == 1) {
    c->p_legendre_coef = dvector(1, c->p_legendre_number);
  }

  fread(c->p_legendre_coef+1, sizeof(double), c->p_legendre_number, fptr);

  /*----------------------------------------------------------------------*/
  /* DONE reading in the phase function so close down the file.           */
  /*----------------------------------------------------------------------*/

  fclose(fptr);

  /*----------------------------------------------------------------------*/
  /* DONE so exit.                                                        */
  /*----------------------------------------------------------------------*/


}

/******************************************************************************/

void
phasefcns_explicit_read_binary(filename, pwhich, c)
char
  *filename;
int
  pwhich;
Constituents
  *c;
{
  int
    index,
    pnum,
    skip;

  double
    wavelength,
    radiusout,
    radiussig,
    relrefoutr,
    relrefoutc;

  FILE
    *fptr,
    *fopen();

  /*----------------------------------------------------------------------*/
  /* Open the phase function file.                                        */
  /*----------------------------------------------------------------------*/

  if ((fptr=fopen(filename,"r"))==NULL) {
    printf("could not open the phase functions file - exiting!\n");
    exit(0);
  }

  /*----------------------------------------------------------------------*/
  /* Skip to the first point in the file that contains the number of      */
  /* points in the phase function.                                        */
  /*----------------------------------------------------------------------*/

  fseek(fptr,ftell(fptr)+44,0);

  /*----------------------------------------------------------------------*/
  /* Read in the number of points in the first phase function.  We assume */
  /* that all of the phase functions in this file have the same number of */
  /* points.                                                              */
  /*----------------------------------------------------------------------*/

  fread(&pnum,sizeof(int),1,fptr);

  /*----------------------------------------------------------------------*/
  /* Go back to the 0th byte of the file where the phase function         */
  /* information starts.                                                  */
  /*----------------------------------------------------------------------*/

  fseek(fptr, 0, 0);

  /*----------------------------------------------------------------------*/
  /* Now skip to the desired phase function using the fact that there are */
  /* 96 bytes in the header information for each phase function.          */
  /*----------------------------------------------------------------------*/

  skip = (pwhich-1)*(96 + 2*pnum*sizeof(double));
  fseek(fptr,ftell(fptr)+skip,0);

  /*----------------------------------------------------------------------*/
  /* Read in the current phase function.                                  */
  /*----------------------------------------------------------------------*/

  fread(&index                  ,sizeof(int),   1,fptr);
  fread(&wavelength 		,sizeof(double),1,fptr);
  fread(&radiusout  		,sizeof(double),1,fptr);
  fread(&radiussig  		,sizeof(double),1,fptr);
  fread(&relrefoutr 		,sizeof(double),1,fptr);
  fread(&relrefoutc 		,sizeof(double),1,fptr);
  fread(&c->p_explicit_number 	,sizeof(int)   ,1,fptr);
  fread(&c->p_g			,sizeof(double),1,fptr);
  fread(&c->w0			,sizeof(double),1,fptr);
  fread(&c->sig_ext		,sizeof(double),1,fptr);
  fread(&c->sig_sca		,sizeof(double),1,fptr);
  fread(&c->sig_abs		,sizeof(double),1,fptr);
  fread(&c->p_normalization	,sizeof(double),1,fptr);

  /*----------------------------------------------------------------------*/
  /* Read in the phase function.                                          */
  /*----------------------------------------------------------------------*/

  if (pwhich == 1) {
    c->p      = dvector(1, c->p_explicit_number);
    c->pdfcno = dvector(1, c->p_explicit_number);
    c->pdfcnf = dvector(1, c->p_explicit_number);
    c->pmu    = dvector(1, c->p_explicit_number);
  }

  fread(c->pmu+1, sizeof(double), c->p_explicit_number, fptr);
  fread(c->p+1  , sizeof(double), c->p_explicit_number, fptr);

  /*----------------------------------------------------------------------*/
  /* DONE reading in the phase function so close down the file.           */
  /*----------------------------------------------------------------------*/

  fclose(fptr);

  /*----------------------------------------------------------------------*/
  /* DONE so exit.                                                        */
  /*----------------------------------------------------------------------*/


}

/******************************************************************************/
/******************************************************************************/
