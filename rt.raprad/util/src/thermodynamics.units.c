/******************************************************************************/
/******************************************************************************/
/* List of Thermodynamic Functions and Their Inputs and Outputs. */
/******************************************************************************/
/*
	double numberdensity_dryair();
		INPUT: Air Density (#/Volume), Water Mixing Ratio (gm/gm)
		OUTPUT: Dry Air Density (#/Volume)

	double mxratio_h2o_gmperkg_numberdensity();
		INPUT: Mixing Ratio (gm/kg), Air Number Density (#/m^3)
		OUTPUT: Number Density (#/m^3)

	double mxratio_o3_gmpergm_numberdensity();
		INPUT: Mixing Ratio (gm/gm), Air Number Density (#/m^3)
		OUTPUT: Number Density (#/m^3)

	double mxratio_o3_gmperkg_numberdensity();
		INPUT: Mixing Ratio (gm/kg), Air Number Density (#/m^3)
		OUTPUT: Number Density (#/m^3)

	double mxratio_co2_gmperkg_numberdensity();
		INPUT: Mixing Ratio (gm/kg), Air Number Density (#/m^3)
		OUTPUT: Number Density (#/m^3)

	double mxratio_o2_gmperkg_numberdensity();
		INPUT: Mixing Ratio (gm/kg), Air Number Density (#/m^3)
		OUTPUT: Number Density (#/m^3)

	double mxratio_volpervol_numberdensity();
		INPUT: Mixing Ratio (volume %), Air Number Density (#/m^3)
		OUTPUT: Number Density (#/m^3)

	double mxratio_volpervol_gmperkg();
		INPUT: Mixing Ratio (volume %), Gram Molecular Weight (gm/mol)
		OUTPUT: Mixing Ratio (gm per kg)

	double mxratio_ppmv_numberdensity();
		INPUT: Mixing Ratio (ppmv), Air Number Density (#/m^3)
		OUTPUT: Number Density (#/m^3)

	double mxratio_numberdensity_cmatmpercm();
		INPUT: Number Density (#/cm^3)
		OUTPUT: Centimeters-Amagats per Centimeter of Gas  (cm-atm/cm)

	double mxratio_numberdensity_precmpercm();
		INPUT: Number Density (#/cm^3)
		OUTPUT:Precipitable-Centimeters per Centimeter of Gas (precm/cm)

	double mxratio_precmpercm_numberdensity();
		INPUT:Precipitable-Centimeters per Centimeter of Gas (precm/cm)
		OUTPUT: Number Density (#/cm^3)

	double pVnRT_numberdensity();
		INPUT: Temperature (K), Pressure (Pa)
		OUTPUT: Number Density (#/m^3)

	double pVnRT_pressure();
		INPUT: Temperature (K), Number Density (#/m^3)
		OUTPUT: Pressure (Pa)

	double hydrostatic_airmass();
		INPUT: Pressure Difference (Pa)
		OUTPUT: Air Mass (kg/m^2)

	double mass_gmpervol_numberdensity();
		INPUT: Mass of substance in column layer (gm/m^2)
		OUTPUT: Number density (#/m^2)

*/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

/******************************************************************************/

#define G                  9.80665        /* Acceleration Due to Gravity */

#define A0                 6.022045E+23   /* #/mole  */
#define RGAS               8.31441 	/* J/mol/K */
#define LOSCMHIDT_VOLUME   2.24E+04   /* cm^3 */

#define GAIR               28.964 	/* Molecular Weight of Air (g/mol) */
#define GCARBONDIOXIDE     44.00995 /* Molecular Weight of Carbon Dioxide (g/mol) */
#define GWATERVAPOR        18.016 	/* Molecular Weight of Water Vapor (g/mol) */
#define GOXYGEN            31.9988  /* Molecular Weight of Ozone (g/mol) */
#define GOZONE             47.9982  /* Molecular Weight of Ozone (g/mol) */
#define GNITROUSOXIDE      44.0128  /* Molecular Weight of Nitrous Oxide (g/mol) */
#define GMETHANE           16.0426  /* Molecular Weight of Methane (g/mol) */

#define DENSITY_WATER      1.000   /* gm/cm^3 */

/******************************************************************************/
/******************************************************************************/

double numberdensity_dryair(double r, double Na)
{
  double alpha, Nw, Nda;

  Nw = (r/GWATERVAPOR)*A0;
  Nda = (1./GAIR)*A0;
  alpha = Na / (Nw+Nda);

  return ((double) (alpha*Nda));
}

/******************************************************************************/

double mxratio_h2o_gmperkg_numberdensity(double r, double Na)
{
  return ((double) ((GAIR/GWATERVAPOR)*Na*(r/1000.)));
}

/******************************************************************************/

double mxratio_o3_gmperkg_numberdensity(double r, double Na)
{
  return ((double) ((GAIR/GOZONE)*Na*(r/1000.)));
}

/******************************************************************************/

double mxratio_o3_gmpergm_numberdensity(double r, double Na)
{
  return ((double) ((GAIR/GOZONE)*Na*r));
}

/******************************************************************************/

double mxratio_co2_gmperkg_numberdensity(double r, double Na)
{
  return ((double) ((GAIR/GCARBONDIOXIDE)*Na*(r/1000.)));
}

/******************************************************************************/

double mxratio_o2_gmperkg_numberdensity(double r, double Na)
{
  return ((double) ((GAIR/GOXYGEN)*Na*(r/1000.)));
}

/******************************************************************************/

double mxratio_volpervol_numberdensity(double r, double Na)
{
  return ((double) (r*Na));
}

/******************************************************************************/

double mxratio_volpervol_gmperkg(double r, double g)
{
  return ((double) ((r*g/GAIR)*1000.));
}

/******************************************************************************/

double mxratio_ppmv_numberdensity(double r, double Na)
{
  return ((double) ((r/1.E+06)*Na));
}

/******************************************************************************/

double mxratio_numberdensity_cmatmpercm(double Na)
{
  return ((double) ((Na*LOSCMHIDT_VOLUME) / A0));
}

/******************************************************************************/

double mxratio_numberdensity_precmpercm(double Na)
{
  return ((double) ((Na*GWATERVAPOR) / (DENSITY_WATER*A0)));
}

/******************************************************************************/

double mxratio_precmpercm_numberdensity(double u)
{
  return ((double) ((u*DENSITY_WATER*A0) / GWATERVAPOR));
}

/******************************************************************************/

double pVnRT_numberdensity(double T, double p)
{
  return ((double) ((p*A0) / (RGAS*T)));
}

/******************************************************************************/

double pVnRT_pressure(double T, double N)
{
  return ((double) (N * (RGAS/A0) * T));
}

/******************************************************************************/

double hydrostatic_airmass(double deltap)
{
  return ((double) (deltap / G));
}

/******************************************************************************/

double mass_gmpervol_numberdensity(double m, double GramWt)
{
  return ((double) ((m*A0) / GramWt));
}

/******************************************************************************/
/******************************************************************************/
