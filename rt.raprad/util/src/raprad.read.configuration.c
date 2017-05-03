/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------*/

int setup(FILE *);

void myfscanf(char *, FILE *);

void read_configuration(
   char *configfileName, char *configfilePhotonPartition,
   char *configfileConstituents, char *configfileSpectralModel, char *configfileCntnm,
   char *configfileAtmosphere, char *configfilePhotonSpace, char *configfileSunView,
   char *configfileQuadrature, char *configfileBrdfSelection, char *configfileBrdf, char *outfileResults
)
{
  FILE *fpta;

  /*--------------------------------------------------------------------------*/
  /* Read all of the input file names from configfileName.                    */
  /*--------------------------------------------------------------------------*/

  if ((fpta=fopen(configfileName,"r"))==NULL) {
    printf("configuration file cannot be opened for some reason - EXITING!\n");
    exit(1);
  }
  else {

    if (setup(fpta)) { myfscanf(configfilePhotonPartition, fpta); }
    if (setup(fpta)) { myfscanf(   configfileConstituents, fpta); }
    if (setup(fpta)) { myfscanf(  configfileSpectralModel, fpta); }
    if (setup(fpta)) { myfscanf(          configfileCntnm, fpta); }
    if (setup(fpta)) { myfscanf(     configfileAtmosphere, fpta); }
    if (setup(fpta)) { myfscanf(    configfilePhotonSpace, fpta); }
    if (setup(fpta)) { myfscanf(        configfileSunView, fpta); }
    if (setup(fpta)) { myfscanf(     configfileQuadrature, fpta); }
    if (setup(fpta)) { myfscanf(  configfileBrdfSelection, fpta); }
    if (setup(fpta)) { myfscanf(           configfileBrdf, fpta); }
    if (setup(fpta)) { myfscanf(           outfileResults, fpta); }

  }

  /*-------------------------------------------------------------------------*/

  fclose(fpta);

  /*-------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
