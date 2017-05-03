/****************************************************************/
/****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "../include/numrec.nrutil.h"
#include "../include/Angles.h"

/*---------------------------------------------------------------*/


void
check_geometry_sunzenith(suna)
Sun
  *suna;
{
  int
    i,
    j;

  FILE
    *finput,
    *fopen();

  /*-------------------------------------------------------------*/
  /* Open check.geometry.sunview                                 */
  /*-------------------------------------------------------------*/

  if((finput=fopen("../../output/check.geometry.sunview","w+")) == NULL) {
    printf ("Cannot create file check.geometry.sunview - Exiting!\n");
    exit(0);
  }

  /*-------------------------------------------------------------*/
  /* Print the Sun angle information.                            */
  /*-------------------------------------------------------------*/

  fprintf(finput, "[sunz_number]           %5d\n", suna->sunz_number);
  fprintf(finput, "[sunz_theta]          %8.4f\n", suna->sunz_theta);
  fprintf(finput, "[sun2earth_distance]     %e\n", suna->sun2earth_distance);
  fprintf(finput, "[gmt_time]           %10.4f\n", suna->gmt_time);

  /*-------------------------------------------------------------*/

  fclose(finput);

  /*-------------------------------------------------------------*/

}

/****************************************************************/
/****************************************************************/
