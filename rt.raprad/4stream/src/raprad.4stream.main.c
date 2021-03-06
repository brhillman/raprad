/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../../util/include/PhotonSpace.h"
#include "../../util/include/PhotonPartition.h"
#include "../../util/include/SpectralModel.h"
#include "../../util/include/Atmosphere.h"
#include "../../util/include/Constituents.h"
#include "../../util/include/Angles.h"
#include "../../util/include/Brdf.h"
#include "../../util/include/Rt1d.h"

/******************************************************************************/

void
  atmosphere_gasamounts_u(),
  atmosphere_layers(),
  atmosphere_layers_subdivide(),
  atmosphere_layers_tau(),
  atmosphere_read_mcclatchey(),
  check_atmosphere_gasamounts_u(),
  check_atmosphere_layers(),
  check_atmosphere_layers_subdivide(),
  check_brdf_selection(),
  check_constituents(),
  check_geometry_sunzenith(),
  check_photon_partition(),
  check_photon_space_final(),
  check_spectralmodel_kato(),
  photon_space_final(),
  rapad_alloc(),
  read_configuration(),
  read_constituents_phasefcns(),
  read_geometry_sunzenith(),
  read_photon_partition(),
  read_photon_space();

Constituents
  *read_constituents();

int
  read_brdf_selection(),
  spectralmodel_read_kato();

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

  Brdf
    *d;

  Sun
    *suna;

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

      spectralmodel_read_kato(configfileSpectralModel, pp, ps, sm);
      if (inputcheck) { check_spectralmodel_kato(pp, sm); }
      break;

      /*----------------------------------------------------------------------*/
      /* DONE reading in the Kato gaseous absorption model.                   */
      /*----------------------------------------------------------------------*/

    }
    else if (!strcmp(c[n].name, "GaseousAbsorptionSpectralModelMlawerLw")) {

      /*----------------------------------------------------------------------*/
      /* Read in the Mlawer gaseous absorption model.                         */
      /*----------------------------------------------------------------------*/

      spectralmodel_read_mlawer_lw(configfileSpectralModel, pp, ps, sm);
      if (inputcheck) { check_spectralmodel_mlawer_lw(pp, sm); }

      /*----------------------------------------------------------------------*/
      /* DONE reading in the Mlawer gaseous absorption model.                 */
      /*----------------------------------------------------------------------*/

    }
    else if (!strcmp(c[n].name, "GaseousAbsorptionSpectralModelPollack")) {

      /*----------------------------------------------------------------------*/
      /* Read in the Pollack gaseous absorption model.                        */
      /*----------------------------------------------------------------------*/

      spectralmodel_read_pollack(configfileSpectralModel, pp, sm);
      if (inputcheck) { check_spectralmodel_pollack(pp, sm); }

      /*----------------------------------------------------------------------*/
      /* Extract the water vapor continuum absorption coefficients as         */
      /* specified by the flags in structure sm.  The water vapor continuum   */
      /* model is intimately tied to the gaseous absorption spectral model    */
      /* and must have the same number of bands.                              */
      /*----------------------------------------------------------------------*/

      cntnmmodel_read(configfileCntnm, pp, sm);
      if (inputcheck) { check_cntnmmodel(pp, sm); }

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

  if (inputcheck) { check_atmosphere_layers_subdivide(atmrt1d); }

  /*--------------------------------------------------------------------------*/
  /* Compute the gas amounts in each layer in units of cm.                    */
  /*--------------------------------------------------------------------------*/

  atmosphere_gasamounts_u(atmrt1d);

  if (inputcheck) { check_atmosphere_gasamounts_u(atmrt1d); }

  /*--------------------------------------------------------------------------*/
  /* Get ILLUMINATION information from configfileSunZenith.                   */
  /*--------------------------------------------------------------------------*/

  suna = (Sun *) malloc(sizeof(Sun));

  read_geometry_sunzenith(configfileSunView, suna);

  if (inputcheck) { check_geometry_sunzenith(suna); }

  /*--------------------------------------------------------------------------*/
  /* Choose the surface BRDFs to be used in the simulation.                   */
  /*--------------------------------------------------------------------------*/

  d = (Brdf *) malloc(sizeof(Brdf));

  read_brdf_selection(configfileBrdfSelection, pp, d);

  if (inputcheck) { check_brdf_selection(d); }

  /*--------------------------------------------------------------------------*/
  /* Allocate space for arrays to hold layer scattering properties.           */
  /*--------------------------------------------------------------------------*/

  rt = (Rt *) malloc(sizeof(Rt));

  rt->rt_number     = ps->rt_number;
  rt->layers_number = ps->gridnumf;

  raprad_alloc(pp, ps, rt);

  /*--------------------------------------------------------------------------*/
  /* Build the atmosphere spectral-interval-by-spectral-interval and run the  */
  /* adding doubling code for each spectral interval .                        */
  /*--------------------------------------------------------------------------*/

  for (i=1; i<=pp->numintervals; ++i) {

    /*------------------------------------------------------------------------*/
    /* Print the interval number to the user as feedback.                     */
    /*------------------------------------------------------------------------*/

    printf("Interval %3d\n", i);

    /*------------------------------------------------------------------------*/
    /* Phase functions for each constituent at this particular wavelength.    */
    /*------------------------------------------------------------------------*/

    read_constituents_phasefcns(i, ps, c);

    /*------------------------------------------------------------------------*/
    /* For each interval or band, go through all of the subintervals.         */
    /*------------------------------------------------------------------------*/

    for (j=1; j<=pp->runsperinterval[i]; ++j) {

      /*----------------------------------------------------------------------*/
      /* Optical depths for each constituent within each layer.               */
      /*----------------------------------------------------------------------*/

      atmosphere_layers_tau(i, j, ps, pp, sm, atmrt1d, c, rt, d);

      /*----------------------------------------------------------------------*/
      /* Setup all of the atmosphere layer information needed in the          */
      /* radiative transfer code.                                             */
      /*----------------------------------------------------------------------*/

      atmosphere_layers(i, j, ps, sm, c, rt, atm);

      if (inputcheck) { check_atmosphere_layers(i, sm, rt, suna, d); }

      /*----------------------------------------------------------------------*/
      /* Go do the adding-doubling radiative transfer for this particular     */
      /* interval or sub-interval.                                            */
      /*----------------------------------------------------------------------*/

      fourstr_(&sm->solarinsol[i], &rt->layers_number, &suna->sunz_mu, &rt->layers_w0[1], 
              &rt->layers_legendre_coef[1], &d->albedo[i][1], &rt->layers_tau[1],
              &sm->alpha[i][j], &suna->sun2earth_distance);

      /*----------------------------------------------------------------------*/
      /* DONE with the radiative transfer.                                    */
      /*----------------------------------------------------------------------*/

    }

    /*------------------------------------------------------------------------*/
    /* DONE with this subinterval so go to the next one.                      */
    /*------------------------------------------------------------------------*/

  }

  /*--------------------------------------------------------------------------*/
  /* DONE with the spectral band calculations.                                */
  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/


