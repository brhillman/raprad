/******************************************************************************/
/* Routine to read atmospheric thermodynamic variables */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "../include/Atmosphere.h"

/*----------------------------------------------------------------------------*/

#define GAIR			28.964 	/* Molecular Weight of Air (g/mol) */
#define GCARBONDIOXIDE         44.00995 /* Molecular Weight of Carbon Dioxide (g/mol) */
#define GWATERVAPOR	        18.016 	/* Molecular Weight of Water Vapor (g/mol) */
#define GOXYGEN		       31.9988  /* Molecular Weight of Oxygen (g/mol) */
#define GOZONE		       47.9982  /* Molecular Weight of Ozone (g/mol) */
#define GNITROUSOXIDE          44.0128  /* Molecular Weight of Nitrous Oxide (g/mol) */
#define GMETHANE	       16.0426  /* Molecular Weight of Methane (g/mol) */
#define A0			6.022045E+23   /* #/mole  */

/*----------------------------------------------------------------------------*/

double mxratio_h2o_gmperkg_numberdensity();
double mxratio_o3_gmperkg_numberdensity();
double mxratio_co2_gmperkg_numberdensity();
double mxratio_o2_gmperkg_numberdensity();

double mxratio_numberdensity_precmpercm();
double mxratio_numberdensity_cmatmpercm();

double mass_gmpervol_numberdensity();

/******************************************************************************/

void atmosphere_gasamounts_u(atm)
  Atmosphere
    *atm;
{
  int
    i;
  double
    airdensity,
    mass,
    massdensity,
    mixingratio,
    numberdensity,
    dz;

  /*--------------------------------------------------------------------------*/
  /* Allocate space for the layer gas amounts in terms of precipitable cm for */
  /* water vapor and amagats for CO2, O3, and O2.                             */
  /*--------------------------------------------------------------------------*/

  atm->uh2oprecm   = (double *) malloc(atm->numlayers*sizeof(double));
  atm->uco2amagats = (double *) malloc(atm->numlayers*sizeof(double));
  atm->uo3amagats  = (double *) malloc(atm->numlayers*sizeof(double));
  atm->uo2amagats  = (double *) malloc(atm->numlayers*sizeof(double));

  /*--------------------------------------------------------------------------*/
  /* Allocate space for the layer gas amounts in terms of # / cm^-2.          */
  /*--------------------------------------------------------------------------*/

  atm->uh2o 	 = (double *) malloc(atm->numlayers*sizeof(double));
  atm->uco2 	 = (double *) malloc(atm->numlayers*sizeof(double));
  atm->uo3  	 = (double *) malloc(atm->numlayers*sizeof(double));
  atm->uo2  	 = (double *) malloc(atm->numlayers*sizeof(double));
  atm->un2o  	 = (double *) malloc(atm->numlayers*sizeof(double));
  atm->uch4  	 = (double *) malloc(atm->numlayers*sizeof(double));

  /*--------------------------------------------------------------------------*/
  /* Loop over the atmospheric layers and compute the gas amounts per layer.  */
  /*--------------------------------------------------------------------------*/

  for (i=0; i<atm->numlayers-1; i++) {

    /*------------------------------------------------------------------------*/
    /* Thickness of the atmospheric layer in meters.                          */
    /*------------------------------------------------------------------------*/

    dz = atm->altitude[i+1]-atm->altitude[i];

    /*------------------------------------------------------------------------*/
    /* Mass density in kg per meter^3.                                        */
    /*------------------------------------------------------------------------*/

    massdensity = atm->airmass[i] / dz;

    /*------------------------------------------------------------------------*/
    /* Mass density in gm per meter^3.                                        */
    /*------------------------------------------------------------------------*/

    massdensity *= 1.E+03;

    /*------------------------------------------------------------------------*/
    /* Number density in Particles per meter^3.                               */
    /*------------------------------------------------------------------------*/

    airdensity = (massdensity/GAIR)*A0;
    atm->airdensavg[i] = airdensity;

    /*------------------------------------------------------------------------*/
    /* Compute the layer gas amounts in amagats and precipitable cm.          */
    /*------------------------------------------------------------------------*/

    mixingratio       = (atm->mxratioh2o[i+1] + atm->mxratioh2o[i]) / 2.;
    numberdensity     = mxratio_h2o_gmperkg_numberdensity(mixingratio, airdensity);
    atm->uh2oprecm[i] = mxratio_numberdensity_precmpercm(numberdensity*1.E-06);

    mixingratio         = (atm->mxratioco2[i+1] + atm->mxratioco2[i]) / 2.;
    numberdensity       = mxratio_co2_gmperkg_numberdensity(mixingratio, airdensity);
    atm->uco2amagats[i] = mxratio_numberdensity_cmatmpercm(numberdensity*1.E-06);

    mixingratio        = (atm->mxratioo3[i+1] + atm->mxratioo3[i]) / 2.;
    numberdensity      = mxratio_o3_gmperkg_numberdensity(mixingratio, airdensity);
    atm->uo3amagats[i] = mxratio_numberdensity_cmatmpercm(numberdensity*1.E-06);

    mixingratio        = (atm->mxratioo2[i+1] + atm->mxratioo2[i]) / 2.;
    numberdensity      = mxratio_o2_gmperkg_numberdensity(mixingratio, airdensity);
    atm->uo2amagats[i] = mxratio_numberdensity_cmatmpercm(numberdensity*1.E-06);

    /*------------------------------------------------------------------------*/
    /* Compute the layer gas amounts in # / cm^2.                             */
    /*------------------------------------------------------------------------*/

    mixingratio = (atm->mxratioh2o[i+1] + atm->mxratioh2o[i]) / 2.;
    mass = atm->airmass[i] * mixingratio;
    numberdensity = mass_gmpervol_numberdensity(mass, GWATERVAPOR);
    atm->uh2o[i] = numberdensity / 1.E04;

    mixingratio = (atm->mxratioco2[i+1] + atm->mxratioco2[i]) / 2.;
    mass = atm->airmass[i] * mixingratio;
    numberdensity = mass_gmpervol_numberdensity(mass, GCARBONDIOXIDE);
    atm->uco2[i] = numberdensity / 1.E04;

    mixingratio = (atm->mxratioo3[i+1] + atm->mxratioo3[i]) / 2.;
    mass = atm->airmass[i] * mixingratio;
    numberdensity = mass_gmpervol_numberdensity(mass, GOZONE);
    atm->uo3[i] = numberdensity / 1.E04;

    mixingratio = (atm->mxratioo2[i+1] + atm->mxratioo2[i]) / 2.;
    mass = atm->airmass[i] * mixingratio;
    numberdensity = mass_gmpervol_numberdensity(mass, GOXYGEN);
    atm->uo2[i] = numberdensity / 1.E04;

    mixingratio = (atm->mxration2o[i+1] + atm->mxration2o[i]) / 2.;
    mass = atm->airmass[i] * mixingratio;
    numberdensity = mass_gmpervol_numberdensity(mass, GNITROUSOXIDE);
    atm->un2o[i] = numberdensity / 1.E04;

    mixingratio = (atm->mxratioch4[i+1] + atm->mxratioch4[i]) / 2.;
    mass = atm->airmass[i] * mixingratio;
    numberdensity = mass_gmpervol_numberdensity(mass, GMETHANE);
    atm->uch4[i] = numberdensity / 1.E04;

    /*------------------------------------------------------------------------*/
    /* DONE with this layer so go to the next one.                            */
    /*------------------------------------------------------------------------*/

  }

  /*--------------------------------------------------------------------------*/
  /* DONE with all layers so exit the routine.                                */
  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
