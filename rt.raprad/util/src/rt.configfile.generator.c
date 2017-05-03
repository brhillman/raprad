
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>



void rt_configurationfile_generator (num_mapshgt, mapshgt, cloud_type_flag,
       cloud_top_flag, cloud_base_flag, layer_optical_depth, ss_albedo)

short
  *cloud_type_flag;

int
  num_mapshgt,
  *cloud_base_flag,
  *cloud_top_flag;
  
float
  *mapshgt,
  **layer_optical_depth;

double
  ***ss_albedo;
   

{
  char
    *file_name;

  int
    n,
    i,
    num_cloud_layers,
    num_constituents;

  float
    bottom,
    top;

  FILE
    *fptr,
    *fopen();

  file_name = (char *) malloc(256*sizeof(char));
    
  for (n=1; n<=22; n++) {
    if (cloud_type_flag[n]==1) {
      sprintf (file_name, %s%d",
      "../../rt/rt1d.raprad.2str/scenarios/bbflux/cirrus.atm/Rt1d.constituents",n);

      if ((fptr = fopen(file_name, "w")) == NULL) {
        printf ("Cannot open file %d for writing",n);
        exit (0);
      }

      num_cloud_layers = (cloud_top_flag[n] - cloud_base_flag[n] - 1);

      fprintf (fptr, "/*        Constituent Properties           */\n\n");

      num_constituents = num_cloud_layers + 2;
      
      fprintf (fptr, "[constituents_number]   CONSTITUENTS: %d\n",
                                               num_constituents   );
      for (i=1; i<=num_mapshgt; i++) {
        if (i>=cloud_base_flag[n] && i<cloud_top_flag[n]) {

          fprintf (fptr, "   /* CLOUD EXTINCTION CONSTITUENT PROPERTIES */\n");

          fprintf (fptr, "[constituents_name]	CONSTITUENT NAME: CloudExtinction,\n");

      bottom = mapshgt[i]*1000.;

          fprintf (fptr,"[base_height]		 BASE HEIGHT OF THE LAYER (in meters): %f,\n",bottom);

      top = bottom + 250.;

          fprintf (fptr,"[top_height]		  TOP HEIGHT OF THE LAYER (in meters): %f,\n", top);

          fprintf (fptr, "[scale_height]		SCALE HEIGHT OF THE LAYER (in meters): -1.,\n");

          fprintf (fptr, "[cdensity]		NUMBER DENSITY FOR THE CONSTITUENT (#/m^3): -1.,\n");

          fprintf (fptr, "[*tau_abs] 		TOTAL ABSORPTION OPTICAL DEPTH: -1,\n");
  
          fprintf (fptr, "[*tau_sca] 		TOTAL SCATTERING OPTICAL DEPTH: -1.,\n");

          fprintf (fptr, "[*tau] 			TOTAL OPTICAL DEPTH: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d,\n", layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i], layer_optical_depth[n][i] );
 
          fprintf (fptr,"[*sig_ext]		TOTAL EXTINCTION CROSS-SECTION (micron^2): -1.,\n");

          fprintf (fptr, "[*sig_sca]		TOTAL SCATTERING CROSS-SECTION (micron^2): -1.,\n");

          fprintf (fptr, "[*sig_abs]		TOTAL ABSORPTION CROSS-SECTION (micron^2): -1.,\n");

          fprintf (fptr, "[*w0] 			SINGLE SCATTERING ALBEDO:  %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d,\n", ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][1], ss_albedo[n][i][2], ss_albedo[n][i][2], ss_albedo[n][i][2], ss_albedo[n][i][2], ss_albedo[n][i][2], ss_albedo[n][i][2], ss_albedo[n][i][2], ss_albedo[n][i][3], ss_albedo[n][i][3], ss_albedo[n][i][3], ss_albedo[n][i][4], ss_albedo[n][i][4], ss_albedo[n][i][5], ss_albedo[n][i][5], ss_albedo[n][i][6], ss_albedo[n][i][6] ); 

          fprintf (fptr, "[*p_g]			ASYMMETRY PARAMETER: -1.,
          fprintf (fptr, "[phasefcn_file]        	PHASE FUNCTION FILE NAME: ../ancillary/water.951030.1800_1830.dat,

          fprintf (fptr, "[phasefcn_format]       PHASE FUNCTION FORMAT: ascii,\n");

          fprintf (fptr, "[phasefcn_which]	PHASE FUNCTION REPRESENTATION: 1,\n");

          fprintf (fptr, "[p_explicit_number]	NUMBER OF ANGLES IN THE PHASE FUNCTION: -1,\n");

          fprintf (fptr, "[p_legendre_number]	NUMBER OF LEGENDRE POLYNOMIALS IN THE PHASE FUNCTION: -1,\n");

          fprintf (fptr, "[es] 			SATURATION VAPOR PRESSURE OVER WATER OR ICE: ice,\n\n\n");
        }
      }
      fprintf (fptr, "[constituents_name]	CONSTITUENT NAME: RayleighScatter,\n");

      fprintf (fptr, "[base_height]		 BASE HEIGHT OF THE CONSTITUENT (in meters):      0.,\n");

      fprintf (fptr, "[top_height]		  TOP HEIGHT OF THE CONSTITUENT (in meters): 70000.,\n");

      fprintf (fptr, "[scale_height]		SCALE HEIGHT OF THE CONSTITUENT (in meters): -1.,\n");

      fprintf (fptr, "[cdensity]		NUMBER DENSITY FOR THE CONSTITUENT (#/m^3): -1.,\n");

      fprintf (fptr, "[*tau_abs] 		TOTAL ABSORPTION OPTICAL DEPTH: -1.,\n");

      fprintf (fptr, "[*tau_sca] 		TOTAL SCATTERING OPTICAL DEPTH: -1.,\n");

      fprintf (fptr, "[*tau] 			TOTAL OPTICAL DEPTH: -1.,\n");

      fprintf (fptr, "[*sig_ext]		TOTAL EXTINCTION CROSS-SECTION (micron^2): -1.,\n");

      fprintf (fptr, "[*sig_sca]		TOTAL SCATTERING CROSS-SECTION (micron^2): -1.,\n");

      fprintf (fptr, "[*sig_abs]		TOTAL ABSORPTION CROSS-SECTION (micron^2): -1.,\n");

      fprintf (fptr, "[*w0] 			SINGLE SCATTERING ALBEDO: 1.,\n");

      fprintf (fptr, "[*p_g]			ASYMMETRY PARAMETER: -1.,\n");

      fprintf (fptr, "[phasefcn_file]        	PHASE FUNCTION FILE NAME: ../ancillary/rayleigh.phase_function.32,\n");

      fprintf (fptr, "[phasefcn_format]       PHASE FUNCTION FORMAT: ascii,\n");

      fprintf (fptr, "[phasefcn_which]	PHASE FUNCTION REPRESENTATION: 1,\n");

      fprintf (fptr, "[p_explicit_number]	NUMBER OF ANGLES IN THE PHASE FUNCTION: -1.,\n");

      fprintf (fptr, "[p_legendre_number]	NUMBER OF LEGENDRE POLYNOMIALS IN THE PHASE FUNCTION: -1.,\n");

      fprintf (fptr, "[es] 			SATURATION VAPOR PRESSURE OVER WATER OR ICE: none,\n\n\n");


      fprintf (fptr, "[constituents_name]	CONSTITUENT NAME: GaseousAbsorption,\n");

      fprintf (fptr, "[base_height]		 BASE HEIGHT OF THE CONSTITUENT (in meters):      0.,\n");

      fprintf (fptr, "[top_height]		  TOP HEIGHT OF THE CONSTITUENT (in meters): 100000.,\n");

      fprintf (fptr, "[scale_height]		SCALE HEIGHT OF THE CONSTITUENT (in meters): -1.,\n");

      fprintf (fptr, "[cdensity]		NUMBER DENSITY FOR THE CONSTITUENT (#/m^3): -1.,\n");

      fprintf (fptr, "[*tau_abs] 		TOTAL ABSORPTION OPTICAL DEPTH: -1.,\n");

      fprintf (fptr, "[*tau_sca] 		TOTAL SCATTERING OPTICAL DEPTH: -1.,\n");

      fprintf (fptr, "[*tau] 			TOTAL OPTICAL DEPTH: -1.,,\n");

      fprintf (fptr, "[*sig_ext]		TOTAL EXTINCTION CROSS-SECTION (micron^2): -1.,\n");

      fprintf (fptr, "[*sig_sca]		TOTAL SCATTERING CROSS-SECTION (micron^2): -1.,\n");

      fprintf (fptr, "[*sig_abs]		TOTAL ABSORPTION CROSS-SECTION (micron^2): -1.,\n");

      fprintf (fptr, "[*w0] 			SINGLE SCATTERING ALBEDO: -1.,\n");

      fprintf (fptr, "[*p_g]			ASYMMETRY PARAMETER: -1.,,\n");

      fprintf (fptr, "[phasefcn_file]        	PHASE FUNCTION FILE NAME: -1.,\n");

      fprintf (fptr, "[phasefcn_format]       PHASE FUNCTION FORMAT: ascii,\n");

      fprintf (fptr, "[phasefcn_which]	PHASE FUNCTION REPRESENTATION: -1.,\n");

      fprintf (fptr, "[p_explicit_number]	NUMBER OF ANGLES IN THE PHASE FUNCTION: -1.,\n");

      fprintf (fptr, "[p_legendre_number]	NUMBER OF LEGENDRE POLYNOMIALS IN THE PHASE FUNCTION: -1.,\n");

      fprintf (fptr, "[es] 			SATURATION VAPOR PRESSURE OVER WATER OR ICE: none,\n");

    }
    fclose(fptr);
  }
}

     
