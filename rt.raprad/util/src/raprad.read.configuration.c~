/******************************************************************************/
/******************************************************************************/

#include <stdio.h>

/*----------------------------------------------------------------------------*/

void
read_configuration(configfileName, configfilePhotonPartition,
     configfileConstituents, configfileSpectralModel, configfileCntnm,
     configfileAtmosphere, configfilePhotonSpace, configfileSunView,
     configfileQuadrature, configfileBrdfSelection, configfileBrdf, outfileResults)
  char
    *configfileName,
    *configfilePhotonPartition,
    *configfileConstituents,
    *configfileSpectralModel,
    *configfileCntnm,
    *configfileAtmosphere,
    *configfilePhotonSpace,
    *configfileSunView,
    *configfileQuadrature,
    *configfileBrdfSelection,
    *configfileBrdf,
    *outfileResults;
{
  FILE
    *fpta;

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
