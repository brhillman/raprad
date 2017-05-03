/******************************************************************************/
/******************************************************************************/
/* List of Thermodynamic Functions and Their Inputs and Outputs.
/******************************************************************************/
/*
	double latenheat();
		INPUT: Temperature (K) | OUTPUT: Latent Heat Water-to-Vapor (J/gm)

	double clausiusclapyron_water_integration();
		INPUT: Temperature (K) | OUTPUT: Water Vapor Pressure (Pa)

	double clausiusclapyron_water_bolton();
		INPUT: Temperature (K) | OUTPUT: Water Vapor Pressure (Pa)

	double clausiusclapyron_water_starr();
		INPUT: Temperature (K) | OUTPUT: Water Vapor Pressure (Pa)

	double clausiusclapyron_water_loweficke();
		INPUT: Temperature (K) | OUTPUT: Water Vapor Pressure (Pa)

	double clausiusclapyron_ice_loweficke();
		INPUT: Temperature (K) | OUTPUT: Water Vapor Pressure (Pa)

	double clausiusclapyron_water_constant_L();
		INPUT: Temperature (K) | OUTPUT: Water Vapor Pressure (Pa)

	double clausiusclapyron_ice_constant_L();
		INPUT: Temperature (K) | OUTPUT: Water Vapor Pressure (Pa)

	double clausiusclapyron_water_variable_L();
		INPUT: Temperature (K) | OUTPUT: Water Vapor Pressure (Pa)

*/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

/******************************************************************************/
/* These values come from Rogers and Yau.                                     */

#define RGAS                          8.31441 	/* J/mol/K */
#define RVAPOR                         461.5 	/* J/kg/K  */
#define T0 			       273.15   /* Kelvin  */
#define E0WATER			       611.21   /* Pascals */
#define E0ICE			       611.15   /* Pascals */
#define CL                             4187.0   /* J/kg/K  */
#define CPV                            1870.0   /* J/kg/K  */
#define LWATER                      2.501E+06   /* J/kg */
#define LICE                        2.834E+06   /* J/kg */

#define GWATERVAPOR		       18.016 	/* Molecular Weight of Water Vapor (g/mol) */

#define  PRECISION                        15
#define  MODIFIER                         1000.000000

/******************************************************************************/

double latentheat();

double clausiusclapyron_water_integration();
	double clausiusclapyron_exponent();
	double trapezoid_integration();
	double richardson_extrapolation();

double clausiusclapyron_water_bolton();
double clausiusclapyron_water_starr();
double clausiusclapyron_water_loweficke();
double clausiusclapyron_ice_loweficke();
double clausiusclapyron_water_constant_L();
double clausiusclapyron_ice_constant_L();
double clausiusclapyron_water_variable_L();

/******************************************************************************/
/******************************************************************************/
/* Function to calculate the latent heat as a function of temperature:        */
/* L(T) = L(T0) - (CL - CPV)*(T - T0) = 2500. - (4.187 - 1.870)*(T - T0)      */

double
latentheat(T)
double T;
{
  double lat_heat;

  lat_heat = 2500.0 - 2.317*(T - T0);  /* J/g */

  return lat_heat;

}

/******************************************************************************/
/******************************************************************************/
/* By Romberg Integration and Richardson Extrapolation */

double
clausiusclapyron_water_integration(T)
double T;
{
  int    i;
  double  **rich, saturationvaporpressure, change;

  /*--------------------------------------------------------------------------*/

  rich = (double **) malloc(PRECISION*sizeof(double *));
  for (i=0;i<PRECISION;++i) {
    rich[i] = (double *) malloc(PRECISION*sizeof(double));
  }

  /*--------------------------------------------------------------------------*/

  rich[0][0] = .5*(T-T0)*(clausiusclapyron_exponent(T0)+clausiusclapyron_exponent(T));
  rich[0][0] *= MODIFIER;
  i = 1;
  change = 999.;

  /*--------------------------------------------------------------------------*/

  while (i<PRECISION&&change>.000001) {
    rich[i][0] = trapezoid_integration(clausiusclapyron_exponent, T0, T, i-1);
    rich[i][0] *= MODIFIER;
    richardson_extrapolation(rich, i);
    change = fabs((rich[i][i]-rich[i-1][i-1])/rich[i][i]);
    ++i;
  }

  if (i==PRECISION&&change>.000001) {
    printf("Clausius Clapyron integration did not converge - exiting\n");
    exit(0);
  }

  /*--------------------------------------------------------------------------*/

  saturationvaporpressure = (E0WATER)*exp(rich[i-1][i-1]/MODIFIER); /* Pa */

  /*--------------------------------------------------------------------------*/

  return saturationvaporpressure;

  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/

double
clausiusclapyron_exponent(T)
double T;
{
  double R_watervapor, value;

  R_watervapor = RGAS / GWATERVAPOR; /* J/g/K */

  value = latentheat(T) / (R_watervapor*T*T ); /* 1/K */
  return value;

}

/******************************************************************************/
/* Integration routine using the trapezoidal rule:
       a: lower bound
       b: upper bound
       i: scale factor to determine the number of integration points ( 2^i )
       y: pointer to function containing the integrand
*/

double trapezoid_integration(y,a,b,i)
double (*y)();
double a,b;
int   i;
{

  double trap,sum,x,dx;
  int   j,np;

  np = (int) pow((double) 2,(double) i);
  dx = (b-a)/np;
  x = a+0.5*dx;
  sum = 0.0;

  for (j=0;j<np;++j) {
    sum += (*y)(x);
    x += dx;
  }

  trap = sum*dx;
  return trap;

}

/******************************************************************************/
/* Richardson Extrapolation Scheme 
	i:      vertical position in matrix to place new value
	matrix: the matrix containing the values
*/

double richardson_extrapolation(matrix,i)
double **matrix;
int   i;
{
  int  j,jpwr;
  double constant;

  for(j=1;j<=i;++j) {
    jpwr = 2*j;
    constant = (double) pow((double) 2,(double) jpwr);
    matrix[i][j] = (constant*matrix[i][j-1] - 
      matrix[i-1][j-1])/(constant-1.0);
  }
  return;

}

/******************************************************************************/
/******************************************************************************/
/* Taken from Bolton and as presented in Rogers and Yau.                      */

double
clausiusclapyron_water_bolton(T)
double T;
{
  double A, B, C, saturationvaporpressure;

  A = 611.2;
  B = 17.67;
  C = 243.5;

  T = T - T0;

  if (T < -30. || T > 35. ) {
/*
     printf("Temperature out of range for clausiusclapyron_water_empirical technique.\n");
     printf("  Returning a default value of 999999.0\n");
*/
     saturationvaporpressure = 999999.0;
     return saturationvaporpressure;
  }

  saturationvaporpressure = A*exp((B*T)/(T+C));

  return saturationvaporpressure;

}

/******************************************************************************/
/******************************************************************************/
/* Taken from Starr 1985, Appendix B,  Page 335                               */
/* Routine to compute the saturation vapor pressure at a specific temperture  */
/*   using a 6th order polynomial.                                            */

double
clausiusclapyron_water_starr(T)
double T;
{
  int i;
  double C0,C1,C2,C3,C4,C5,C6;
  double saturationvaporpressure;

  static double coef[3][7] = { 
     {
         0.6165189065e+04,-0.1667439376e+03,
         0.1886839551e+01,-0.1143897897e-01,
         0.3920338183e-04,-0.7204942038e-07,
         0.5550437565e-10
     },
     {
         0.7676749690e+04,-0.2056793389e+03,
         0.2301179803e+01,-0.1377392223e-01,
         0.4656025967e-04,-0.8434798769e-07,
         0.6403236206e-10
     },
     {
         -0.6831653575e+04, 0.1012101269e+03,
         -0.4024185804e+00,-0.1076593467e-02,
          0.1303140958e-04,-0.3714823651e-07,
          0.3635860564e-10
     }
  };

  if (T < 223. || T > 323. ) {
/*
     printf("Temperature out of range for clausiusclapyron_starr technique.\n");
     printf("  Returning a default value of 999999.0\n");
*/
     saturationvaporpressure = 999999.0;
     return saturationvaporpressure;
  }

  if ( T>=223. && T<258. )  { i = 0; }
  if ( T>=258. && T<283. )  { i = 1; }
  if ( T>=283. && T<=323. ) { i = 2; }

  C0 = coef[i][0];
  C1 = coef[i][1];
  C2 = coef[i][2];
  C3 = coef[i][3];
  C4 = coef[i][4];
  C5 = coef[i][5];
  C6 = coef[i][6];

  saturationvaporpressure = (C0+T*(C1+T*(C2+T*(C3+T*(C4+T*(C5+T*C6))))))*100.; /* Pascals */

  return saturationvaporpressure;

}

/******************************************************************************/
/******************************************************************************/
/* From Pruppacher and Klett Appendices (Lowe and Ficke, 1974)                */

double
clausiusclapyron_water_loweficke(T)
double T;
{
  int i;
  double t0,t1,t2,t3,t4,t5,t6;
  double saturationvaporpressure;

  static double C[7] = { 
	6.107799961,
	4.436518521E-01,
	1.428945805E-02,
	2.650648471E-04,
	3.031240396E-06,
	2.034080948E-08,
        6.136820929E-11
  };

  /*--------------------------------------------------------------------------*/
  /* Convert degrees Kelvin to degrees Centigrade.                            */
  /*--------------------------------------------------------------------------*/

  T = T - T0;

  if (T < -50. || T > 50. ) {
/*
     printf("Temperature out of range for clausiusclapyron_water_empirical technique.\n");
     printf("  Returning a default value of 999999.0\n");
*/
     saturationvaporpressure = 999999.0;
     return saturationvaporpressure;
  }

  saturationvaporpressure = (C[0]+T*(C[1]+T*(C[2]+T*(C[3]+T*(C[4]+T*(C[5]+T*C[6]))))))*100.; /* Pascals */

  return saturationvaporpressure;

}

/******************************************************************************/
/******************************************************************************/
/* From Pruppacher and Klett Appendices (Lowe and Ficke, 1974)                */

double
clausiusclapyron_ice_loweficke(T)
double T;
{
  int i;
  double t0,t1,t2,t3,t4,t5,t6;
  double saturationvaporpressure;

  static double C[7] = { 
	6.109177956,
	5.034698970E-01,
	1.886013408E-02,
	4.176223716E-04,
	5.824720280E-06,
	4.838803174E-08,
        1.838826904E-10
  };

  /*--------------------------------------------------------------------------*/
  /* Convert degrees Kelvin to degrees Centigrade.                            */
  /*--------------------------------------------------------------------------*/

  T = T - T0;

  if (T < -50. || T > 50. ) {
/*
     printf("Temperature out of range for clausiusclapyron_ice_empirical technique.\n");
     printf("  Returning a default value of 999999.0\n");
*/
     saturationvaporpressure = 999999.0;
     return saturationvaporpressure;
  }

  saturationvaporpressure = (C[0]+T*(C[1]+T*(C[2]+T*(C[3]+T*(C[4]+T*(C[5]+T*C[6]))))))*100.; /* Pascals */

  return saturationvaporpressure;

}

/******************************************************************************/
/******************************************************************************/
/* Constant Latent Heat of Vaporization.                                      */

double
clausiusclapyron_water_constant_L(T)
double T;
{
  double A, B, saturationvaporpressure;

  A = E0WATER*exp(LWATER/(RVAPOR*T0));
  B = LWATER/RVAPOR;

  saturationvaporpressure = A*exp(-B/T);

  return saturationvaporpressure;

}

/******************************************************************************/
/******************************************************************************/
/* Constant Latent Heat of Sublimation.                                       */

double
clausiusclapyron_ice_constant_L(T)
double T;
{
  double A, B, saturationvaporpressure;

  A = E0ICE*exp(LICE/(RVAPOR*T0));
  B = LICE/RVAPOR;

  saturationvaporpressure = A*exp(-B/T);

  return saturationvaporpressure;

}

/******************************************************************************/
/******************************************************************************/
/* Constant Latent Heat of Vaporization.                                      */

double
clausiusclapyron_water_variable_L(T)
double T;
{
  double A, B, C, saturationvaporpressure;

  A = (CL-CPV) / RVAPOR;
  B = (LWATER+(CL-CPV)*T0) / RVAPOR;
  C = E0ICE*exp(A*log(T0))*exp(B/T0);

  saturationvaporpressure = C*exp(-A*log(T))*exp(-B/T);

  return saturationvaporpressure;

}

/******************************************************************************/
/******************************************************************************/
