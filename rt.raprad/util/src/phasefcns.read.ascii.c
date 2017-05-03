/**************************************************************************/
/**************************************************************************/

#include <stdio.h>

#include "../include/numrec.nrutil.h"

#include "../include/Constituents.h"

/**************************************************************************/

void
phasefcns_read_ascii_header(filename, interval, c)
char
  *filename;
int
  interval;
Constituents
  *c;
{
  int
    ch,
    i;

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
    printf("could not open the phase function file - exiting!\n");
    exit(0);
  }

  /*----------------------------------------------------------------------*/
  /* Skip to the phase function for the current interval.                 */
  /*----------------------------------------------------------------------*/

  for (i=1; i<=interval; i++) {
    while ((ch=getc(fptr)) != '@');
  }

  /*----------------------------------------------------------------------*/
  /* Read in the needed header information.                               */
  /*----------------------------------------------------------------------*/

  fscanf(fptr, "%lf", &wavelength);
  fscanf(fptr, "%lf", &radiusout);
  fscanf(fptr, "%lf", &radiussig);
  fscanf(fptr, "%lf", &relrefoutr);
  fscanf(fptr, "%lf", &relrefoutc);
  fscanf(fptr, "%d",  &c->p_legendre_number);
  fscanf(fptr, "%lf", &c->p_g);
  fscanf(fptr, "%lf", &c->w0);
  fscanf(fptr, "%lf", &c->sig_ext);
  fscanf(fptr, "%lf", &c->sig_sca);
  fscanf(fptr, "%lf", &c->sig_abs);
  fscanf(fptr, "%lf", &c->p_normalization);

  /*----------------------------------------------------------------------*/
  /* Done reading in the phase function.                                  */
  /*----------------------------------------------------------------------*/

}

/**************************************************************************/

void
phasefcns_legendre_read_ascii(filename, interval, c)
char
  *filename;
int
  interval;
Constituents
  *c;
{
  int
    ch,
    i;

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
    printf("could not open the phase function file - exiting!\n");
    exit(0);
  }

  /*----------------------------------------------------------------------*/
  /* Skip to the phase function for the current interval.                 */
  /*----------------------------------------------------------------------*/

  for (i=1; i<=interval; i++) {
    while ((ch=getc(fptr)) != '@');
  }

  /*----------------------------------------------------------------------*/
  /* Read in the needed header information.                               */
  /*----------------------------------------------------------------------*/

  fscanf(fptr, "%lf", &wavelength);
  fscanf(fptr, "%lf", &radiusout);
  fscanf(fptr, "%lf", &radiussig);
  fscanf(fptr, "%lf", &relrefoutr);
  fscanf(fptr, "%lf", &relrefoutc);
  fscanf(fptr, "%d",  &c->p_legendre_number);
  fscanf(fptr, "%lf", &c->p_g);
  fscanf(fptr, "%lf", &c->w0);
  fscanf(fptr, "%lf", &c->sig_ext);
  fscanf(fptr, "%lf", &c->sig_sca);
  fscanf(fptr, "%lf", &c->sig_abs);
  fscanf(fptr, "%lf", &c->p_normalization);

  /*----------------------------------------------------------------------*/
  /* Read in the legendre polynomial coefficients.                        */
  /*----------------------------------------------------------------------*/

  /* We do this allocation only once. */

  if (interval == 1) {
    c->p_legendre_coef = dvector(1, c->p_legendre_number);
  }

  for (i=1; i<=c->p_legendre_number; i++) {
    fscanf(fptr, "%lf", &c->p_legendre_coef[i]);
  }

  /*----------------------------------------------------------------------*/
  /* Done reading in the phase function.                                  */
  /*----------------------------------------------------------------------*/

  fclose(fptr);

  /*----------------------------------------------------------------------*/
  /* Done closing the file so exit the program.                           */
  /*----------------------------------------------------------------------*/

}

/**************************************************************************/

void
phasefcns_explicit_read_ascii(filename, interval, c)
char
  *filename;
int
  interval;
Constituents
  *c;
{
  int
    ch,
    i;

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
    printf("could not open the phase function file - exiting!\n");
    exit(0);
  }

  /*----------------------------------------------------------------------*/
  /* Skip to the phase function for the current interval.                 */
  /*----------------------------------------------------------------------*/

  for (i=1; i<=interval; i++) {
    while ((ch=getc(fptr)) != '@');
  }

  /*----------------------------------------------------------------------*/
  /* Read in the needed header information.                               */
  /*----------------------------------------------------------------------*/

  fscanf(fptr, "%lf", &wavelength);
  fscanf(fptr, "%lf", &radiusout);
  fscanf(fptr, "%lf", &radiussig);
  fscanf(fptr, "%lf", &relrefoutr);
  fscanf(fptr, "%lf", &relrefoutc);
  fscanf(fptr, "%d",  &c->p_explicit_number);
  fscanf(fptr, "%lf", &c->p_g);
  fscanf(fptr, "%lf", &c->w0);
  fscanf(fptr, "%lf", &c->sig_ext);
  fscanf(fptr, "%lf", &c->sig_sca);
  fscanf(fptr, "%lf", &c->sig_abs);
  fscanf(fptr, "%lf", &c->p_normalization);

  /*----------------------------------------------------------------------*/
  /* Read in the phase function.                                          */
  /*----------------------------------------------------------------------*/

  if (interval == 1) {
    c->p      = dvector(1, c->p_explicit_number);
    c->pdfcno = dvector(1, c->p_explicit_number);
    c->pdfcnf = dvector(1, c->p_explicit_number);
    c->pmu    = dvector(1, c->p_explicit_number);
  }

  for (i=1; i<=c->p_explicit_number; i++) {
    fscanf(fptr, "%lf", &c->pmu[i]);
    fscanf(fptr, "%lf", &c->p[i]);
    c->pdfcno[i] = 0.;
    c->pdfcnf[i] = 0.;
  }

  /*----------------------------------------------------------------------*/
  /* Done reading in the phase function.                                  */
  /*----------------------------------------------------------------------*/

  fclose(fptr);

  /*----------------------------------------------------------------------*/

}

/**************************************************************************/
/**************************************************************************/
