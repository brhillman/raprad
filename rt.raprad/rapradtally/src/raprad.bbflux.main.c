/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../../util/include/PhotonSpace.h"
#include "../../util/include/PhotonPartition.h"
#include "../../util/include/SpectralModel.h"
#include "../../util/include/Atmosphere.h"
#include "../../util/include/Constituents.h"
#include "../../util/include/Angles.h"
#include "../../util/include/Brdf.h"
#include "../../util/include/Rt1d.h"

/******************************************************************************/

extern void atmosphere_gasamounts_u(Atmosphere *atm);
extern void atmosphere_layers(int i, int j, PhotonSpace *ps, SpectralModel *sm, Constituents *c, Rt *rt, Atmosphere *atm);
extern void atmosphere_layers_subdivide(PhotonSpace *ps, Atmosphere *atm, Atmosphere *atmnew);

extern void atmosphere_layers_tau(int i, int j, PhotonSpace *ps, PhotonPartition *pp, SpectralModel *sm, Atmosphere *atm, Constituents *c, Rt *rt, Brdf *d);


extern void atmosphere_read_mcclatchey(char *file_name, Atmosphere *atm);
extern void check_atmosphere_gasamounts_u(Atmosphere *atm);
extern void check_atmosphere_layers(int i, int j, SpectralModel *sm, Rt *rt, Sun *suna, Brdf *d);
extern void check_atmosphere_layers_subdivide(Atmosphere *atmnew);
extern void check_brdf_selection(PhotonPartition *pp, Brdf *b);
extern void check_constituents(PhotonPartition *pp, PhotonSpace *ps, Constituents *c);
extern void check_geometry_sunzenith(Sun *suna);
extern void check_photon_partition(PhotonPartition *pp);
extern void check_photon_space_final(PhotonSpace *ps);

extern void check_spectralmodel_kato(PhotonPartition *pp, SpectralModel *sm);
extern void check_spectralmodel_pollack(PhotonPartition *pp, SpectralModel *sm);
extern void check_cntnmmodel(PhotonPartition *pp, SpectralModel *sm);
extern void cntnmmodel_read(char *filename, PhotonPartition *pp, SpectralModel *sm);
extern void photon_space_final(PhotonSpace *ps, Constituents *c);

extern void read_configuration(
   char *configfileName, char *configfilePhotonPartition,
   char *configfileConstituents, char *configfileSpectralModel, char *configfileCntnm,
   char *configfileAtmosphere, char *configfilePhotonSpace, char *configfileSunView,
   char *configfileQuadrature, char *configfileBrdfSelection, char *configfileBrdf, char *outfileResults
);
extern void read_geometry_sunzenith(char *geometryfile, Sun *suna);
extern void read_photon_partition(char *filename, PhotonPartition *pp);
extern void read_photon_space(char *filename, PhotonSpace *ps);
extern void read_brdf_selection(char *filename, PhotonPartition *pp, Brdf *);
extern void spectralmodel_read_mlawer_lw(char *filename, PhotonPartition *pp, PhotonSpace *ps, SpectralModel *sm);
extern void spectralmodel_read_kato(char *filename, PhotonPartition *pp, SpectralModel *sm);
extern void spectralmodel_read_pollack(char *filename, PhotonPartition *pp, SpectralModel *sm);
extern Constituents *read_constituents(char *constituentsfile, PhotonPartition *pp, PhotonSpace *ps);
extern void bb_flux_(
   double *alt, double *press, float *wave, int *nlayer, int *n_band, int *nsub_band,
   float *u0, float *thz, float *gmt_time, double *sol_const
);


/******************************************************************************/

main(argc, argv)

int
  argc;
char
  **argv;

{
  char
    configfileName[256],		/* Configuration file */
    configfilePhotonPartition[256],	/* Spectral intervals */
    configfilePhotonSpace[256],		/* Atmosphere Layers */
    configfileConstituents[256],        /* Constituents */
    configfileSpectralModel[256],	/* Spectral band model */
    configfileCntnm[256],		/* Continuum model */
    configfileAtmosphere[256],		/* Atmospheric profile data */
    configfileSunView[256],		/* Sun location and view angles */
    configfileQuadrature[256],		/* Rt quadrature angles */
    configfileBrdfSelection[256],	/* BRDF selection menu */
    configfileBrdf[256],		/* BRDF data */
    outfileResults[256];

  int
    i,
    inputcheck,
    j,
    n;

  Atmosphere
    *atm,
    *atmrt1d;

  Constituents
    *c;

  PhotonPartition
    *pp;

  Rt
    *rt;

  PhotonSpace
    *ps;

  SpectralModel
    *sm;

  Angles
    *qa,
    *va;

  Sun
    *suna;

  Brdf
    *d;

  /*--------------------------------------------------------------------------*/
  /* Load the configuration file name and inputcheck flag from the command    */
  /* line.                                                                    */
  /*--------------------------------------------------------------------------*/

  if (argc != 3) {
    printf("You must have 2 command line arguments:\n");
    printf("   1) The configuration file name with some sort of path to it.\n");
    printf("   2) Error check flag: do (1) or do not (0) reflect input to output.\n");
    printf("Exiting!\n");
    exit(0);
  }

  /*--------------------------------------------------------------------------*/
  /* Get the command line arguments.                                          */
  /*--------------------------------------------------------------------------*/

  strcpy(configfileName, argv[1]);

  inputcheck = atoi(argv[2]);

  /*--------------------------------------------------------------------------*/
  /* Read all input file names contained in ../config/rt1d.configuration      */
  /*--------------------------------------------------------------------------*/

  read_configuration(configfileName, configfilePhotonPartition,
     configfileConstituents, configfileSpectralModel, configfileCntnm,
     configfileAtmosphere, configfilePhotonSpace, configfileSunView,
     configfileQuadrature, configfileBrdfSelection, configfileBrdf, 
     outfileResults);

  /*--------------------------------------------------------------------------*/
  /* Structure that contains the users selections of which bands to run in    */
  /* the Spectral Model data base.  It also allows only certain sub-intervals */
  /* within a band to be used if the user so chooses.  The third field about  */
  /* Rayleigh scatter above the cloud layer is a relic from the Monte Carlo   */
  /* code and is never used.                                                  */
  /*--------------------------------------------------------------------------*/

  pp = (PhotonPartition *) malloc(sizeof(PhotonPartition));

  read_photon_partition(configfilePhotonPartition, pp);

  if (inputcheck) { check_photon_partition(pp); }

  /*--------------------------------------------------------------------------*/
  /* Structure that contains information about the number and location of the */
  /* horizontal planes at which transmission and reflection matrices are to   */
  /* be calculated.                                                           */
  /*--------------------------------------------------------------------------*/

  ps = (PhotonSpace *) malloc(sizeof(PhotonSpace));

  read_photon_space(configfilePhotonSpace, ps);

  /*--------------------------------------------------------------------------*/
  /* Read the constituents information from configfileConstituents.  We have  */
  /* allocated space for ps at this point because we need to use its member   */
  /* cnumber which stores the number of constituents in the current run.      */
  /*--------------------------------------------------------------------------*/

  c = read_constituents(configfileConstituents, pp, ps);

  if (inputcheck) { check_constituents(pp, ps, c); }

  /*--------------------------------------------------------------------------*/
  /* Use the constituent locations and the requested layer geometry in ps to  */
  /* build the final set of layers to be used in the radiative transfer code. */
  /* Transmission and reflection matrices are always calculated at planes     */
  /* that bound a constituent layer.                                          */
  /*--------------------------------------------------------------------------*/

  photon_space_final(ps, c);

  if (inputcheck) { check_photon_space_final(ps); }

  /*--------------------------------------------------------------------------*/
  /* Read the spectral domains and the subintervals within the spectral       */
  /* domains.  Currently, the absorption band model constrains the spectral   */
  /* domain; the particle cross-sections, phase functions and so forth are    */
  /* calculated for these intervals.   Furthermore, only those bands flagged  */
  /* in the configfilePhotonPartition data base and stored in the structure   */
  /* pp are extracted from the database.                                      */
  /*--------------------------------------------------------------------------*/

  sm = (SpectralModel *) malloc(sizeof(SpectralModel));

  for (n=1; n<=ps->cnumber; n++) {

    if      (!strcmp(c[n].name, "GaseousAbsorptionSpectralModelKato")) {

      /*----------------------------------------------------------------------*/
      /* Read in the Kato gaseous absorption model.                           */
      /*----------------------------------------------------------------------*/

      spectralmodel_read_kato(configfileSpectralModel, pp, sm);
      if (inputcheck) check_spectralmodel_kato(pp, sm);
      break;

      /*----------------------------------------------------------------------*/
      /* DONE reading in the Kato gaseous absorption model.                   */
      /*----------------------------------------------------------------------*/

    } else if (!strcmp(c[n].name, "GaseousAbsorptionSpectralModelMlawerLw")) {

      /*----------------------------------------------------------------------*/
      /* Read in the Mlawer gaseous absorption model.                         */
      /*----------------------------------------------------------------------*/

      spectralmodel_read_mlawer_lw(configfileSpectralModel, pp, ps, sm);
      if (inputcheck) check_spectralmodel_mlawer_lw(pp, sm);

      /*----------------------------------------------------------------------*/
      /* DONE reading in the Mlawer gaseous absorption model.                 */
      /*----------------------------------------------------------------------*/

    } else if (!strcmp(c[n].name, "GaseousAbsorptionSpectralModelPollack")) {

      /*----------------------------------------------------------------------*/
      /* Read in the Pollack gaseous absorption model.                        */
      /*----------------------------------------------------------------------*/

      spectralmodel_read_pollack(configfileSpectralModel, pp, sm);
      if (inputcheck) check_spectralmodel_pollack(pp, sm);

      /*----------------------------------------------------------------------*/
      /* Extract the water vapor continuum absorption coefficients as         */
      /* specified by the flags in structure sm.  The water vapor continuum   */
      /* model is intimately tied to the gaseous absorption spectral model    */
      /* and must have the same number of bands.                              */
      /*----------------------------------------------------------------------*/

      cntnmmodel_read(configfileCntnm, pp, sm);
      if (inputcheck) check_cntnmmodel(pp, sm);

      break;

      /*----------------------------------------------------------------------*/
      /* DONE reading in the Pollack gaseous absorption model.                */
      /*----------------------------------------------------------------------*/

    }

  }

  /*--------------------------------------------------------------------------*/
  /* Read in the atmospheric profile information.                             */
  /*--------------------------------------------------------------------------*/

  atm  = (Atmosphere *) malloc(sizeof(Atmosphere));

  atmosphere_read_mcclatchey(configfileAtmosphere, atm);

  /*--------------------------------------------------------------------------*/
  /* Based on the planes that emerge from the preceding routine, interpolate  */
  /* all of the atmospheric profile information stored in atm to the plane    */
  /* heights specified in ps->gridf.                                          */
  /*--------------------------------------------------------------------------*/

  atmrt1d = (Atmosphere *) malloc(sizeof(Atmosphere));

  atmosphere_layers_subdivide(ps, atm, atmrt1d);

  if (inputcheck) check_atmosphere_layers_subdivide(atmrt1d);

  /*--------------------------------------------------------------------------*/
  /* Get ILLUMINATION information from configfileSunZenith.                   */
  /*--------------------------------------------------------------------------*/

  suna = (Sun *) malloc(sizeof(Sun));

  read_geometry_sunzenith(configfileSunView, suna);

  if (inputcheck) check_geometry_sunzenith(suna);

  /*--------------------------------------------------------------------------*/
  /* Compute the broadband flux.                                              */
  /*--------------------------------------------------------------------------*/

  bb_flux_(&atmrt1d->altitude[0],   &atmrt1d->pressure[0], &pp->wavelength[1],
           &ps->gridnumf,           &pp->numintervals,
           &pp->runsperinterval[1], &suna->sunz_mu,        &suna->suna_theta,
           &suna->gmt_time,         &sm->solarinsol[1]);

  /*--------------------------------------------------------------------------*/
  /* DONE with the routine so exit.                                           */
  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/


