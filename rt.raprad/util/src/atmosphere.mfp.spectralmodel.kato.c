/**********************************************************************/
/**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/PhotonSpace.h"
#include "../include/PhotonPartition.h"
#include "../include/SpectralModel.h"
#include "../include/Atmosphere.h"
#include "../include/Constituents.h"

double
    bext_spectralmodelkato(),
    bext_gasamounts_atmosphere(),
    bext_watervapor_cntnm_atmosphere(),
    linear_interpolate(),
    log_interpolate();

/**********************************************************************/

void
mfp_spectralmodel_kato(i, j, k, ps, pp, sm, atm, c, cvd, ucvd, bi)

  int
    i,
    j,
    k;
  PhotonSpace
    *ps;
  PhotonPartition
    *pp;
  SpectralModel
    *sm;
  Atmosphere
    *atm;
  Constituents
    *c;
  double
    *cvd,
    *ucvd,
    *bi;
{
  int
    ii,
    ip,
    it,
    iw,
    m,
    mm,
    npts,
    npts_shift,
    iw_shift,
    operatediv,
    operatemod,
    quadptindex;

  double
    ansyrl,
    ansyru,
    frcnx,
    frcny,
    p,
    t,
    w,
    *xsect[9],
    z;

  char
    *moleculeptr;

  FILE
    *fptr,
    *fopen();

  /*--------------------------------------------------------------*/
  /* Mid-layer height.                                            */
  /*--------------------------------------------------------------*/

  z = (ps->gridf[k-1] + ps->gridf[k]) / 2.;

  /*--------------------------------------------------------------*/
  /* Mid-layer pressure.                                          */
  /*--------------------------------------------------------------*/

  p = log_interpolate(atm->pressure, atm->altitude, atm->numlayers, z)/100.;
  p = log10(p);

  /*--------------------------------------------------------------*/
  /* Mid-layer pressure lower bound index into the lookup table.  */
  /*--------------------------------------------------------------*/

  ip = sm->npres;
  while (1 < ip) {
    if (p < sm->preslog[ip]) { break; }
    else { --ip; }
  }

  /*--------------------------------------------------------------*/
  /* Mid-layer temperature.                                       */
  /*--------------------------------------------------------------*/

  t = linear_interpolate(atm->temperature, atm->altitude, atm->numlayers, z);

  /*--------------------------------------------------------------*/
  /* Mid-layer temperature lower bound index into the lookup      */
  /* table.                                                       */
  /*--------------------------------------------------------------*/

  it = sm->ntemp;
  while (1 < it) {
    if (sm->temp[it] <= t) { break; }
    else { --it; }
  }

  /*--------------------------------------------------------------*/
  /* Mid-layer water vapor concentration.                         */
  /*--------------------------------------------------------------*/

  w = log_interpolate(atm->densh2o, atm->altitude, atm->numlayers, z);
  if (w < *cvd) { w = *cvd; }

  w = log10(w);

  /*--------------------------------------------------------------*/
  /* Mid-layer water vapor concentration lower bound index into   */
  /* the lookup table.                                            */
  /*--------------------------------------------------------------*/

  iw = sm->nwatvap;
  while (1 < iw) {
    if (w < sm->watvaplog[iw]) { break; }
    else { --iw; }
  }

  /*--------------------------------------------------------------*/
  /* Loop through the different molecules for this wavelength and */
  /* subinterval, filling in the missing aborption coefficients.  */
  /*--------------------------------------------------------------*/

  for (m=1; m<=sm->nummolecules[i]; m++) {

    /*------------------------------------------------------------*/
    /* Check to see if there is a file associated with this m.    */
    /*------------------------------------------------------------*/

    if (sm->nquadpts[i][m] == 1) continue;

    /*------------------------------------------------------------*/
    /* Based on the molecule index m, find the index of the       */
    /* quadrature point to be used for this, the j^th, interval.  */
    /*------------------------------------------------------------*/

    operatemod = 1;
    for (mm=m; mm<=sm->nummolecules[i]; mm++) {
      operatemod *= sm->nquadpts[i][mm];
    }

    operatediv = 1;
    for (mm=m+1; mm<=sm->nummolecules[i]; mm++) {
      operatediv *= sm->nquadpts[i][mm];
    }

    quadptindex = ((j-1)%operatemod)/operatediv + 1;

    /*------------------------------------------------------------*/
    /* Open the absorption coefficients file.                     */
    /*------------------------------------------------------------*/

    if (   ((moleculeptr = strstr(sm->abscoeffile[i][m],"H2O")) == NULL) ||
           ((moleculeptr = strstr(sm->abscoeffile[i][m], "09")) != NULL) ||
           ((moleculeptr = strstr(sm->abscoeffile[i][m], "10")) != NULL) ||
           ((moleculeptr = strstr(sm->abscoeffile[i][m], "11")) != NULL) ||
           ((moleculeptr = strstr(sm->abscoeffile[i][m], "13")) != NULL) ||
           ((moleculeptr = strstr(sm->abscoeffile[i][m], "15")) != NULL)     ) {

      /*----------------------------------------------------------*/
      /* Read in the relevant part of the interpolation table.    */
      /*----------------------------------------------------------*/

      npts = (ip-1)*sm->ntemp*sm->nquadpts[i][m] + (it-1)*sm->nquadpts[i][m];

      xsect[1] = &sm->abscoefdata[i][m][npts];
      xsect[1] -= 1;

      npts += sm->nquadpts[i][m];

      xsect[4] = &sm->abscoefdata[i][m][npts];
      xsect[4] -= 1;

      npts += (sm->nquadpts[i][m]+(sm->ntemp-2)*sm->nquadpts[i][m]);

      xsect[2] = &sm->abscoefdata[i][m][npts];
      xsect[2] -= 1;

      npts += sm->nquadpts[i][m];

      xsect[3] = &sm->abscoefdata[i][m][npts];
      xsect[3] -= 1;

      /*----------------------------------------------------------*/
      /* Checking interpolation region.                           */
      /*----------------------------------------------------------*/

      for (ii=1; ii<=4; ii++) {
        if (xsect[ii][1] < 0.0) {
          printf("%3d%3d%3d\n", ii, ip, it);
        }
      }

      /*----------------------------------------------------------*/
      /* Do two-dimensional interpolation on pressure and         */
      /* temperature plane.                                       */
      /*----------------------------------------------------------*/

      frcnx = (p-sm->preslog[ip]) / (sm->preslog[ip+1]-sm->preslog[ip]);
      frcny = (t-sm->temp[it]) / (sm->temp[it+1]-sm->temp[it]);

      sm->abscoef[i][m + (j-1)*sm->nummolecules[i]] =
                 (1.0-frcnx)*(1.0-frcny)*xsect[1][quadptindex]
               + frcnx*(1.0-frcny)*xsect[2][quadptindex]
               + frcnx*frcny*xsect[3][quadptindex]
               + (1.0-frcnx)*frcny*xsect[4][quadptindex];

      /*----------------------------------------------------------*/
      /* DONE with the cross section for non-water moecules.      */
      /*----------------------------------------------------------*/

    }
    else {

      /*----------------------------------------------------------*/
      /* Branch according to whether or not you need to           */
      /* interpolate the water vapor concentrations.              */
      /*----------------------------------------------------------*/

      if (iw <= 10) {

        /*--------------------------------------------------------*/
        /* Read in the relevant part of the interpolation table.  */
        /*--------------------------------------------------------*/

        npts = (ip-1)*sm->ntemp*sm->nwatvap*sm->nquadpts[i][m]
             +           (it-1)*sm->nwatvap*sm->nquadpts[i][m]
             +                       (iw-1)*sm->nquadpts[i][m];

        if (sm->abscoefdata[i][m][npts] <= 0.0) {
          npts_shift = npts ;
          iw_shift = iw ;
          while (sm->abscoefdata[i][m][npts_shift] <= 0.0 && iw_shift <= 11) {
            ++iw_shift ;
            npts_shift  = (ip-1)*sm->ntemp*sm->nwatvap*sm->nquadpts[i][m]
             +                      (it-1)*sm->nwatvap*sm->nquadpts[i][m]
             +                                (iw_shift-1)*sm->nquadpts[i][m];
          }

          xsect[1] = &sm->abscoefdata[i][m][npts_shift];
        }
        else {

          xsect[1] = &sm->abscoefdata[i][m][npts];

        }
       
        xsect[1] -= 1;

        npts += sm->nquadpts[i][m];

        if (sm->abscoefdata[i][m][npts] <= 0.0) {
          npts_shift = npts ;
          iw_shift = iw ;
          while (sm->abscoefdata[i][m][npts_shift] <= 0.0 && iw_shift <= 11) {
            ++iw_shift ;
            npts_shift  = (ip-1)*sm->ntemp*sm->nwatvap*sm->nquadpts[i][m]
             +                      (it-1)*sm->nwatvap*sm->nquadpts[i][m]
             +                                (iw_shift-1)*sm->nquadpts[i][m];
          }

          xsect[5] = &sm->abscoefdata[i][m][npts_shift];
        }
        else {

          xsect[5] = &sm->abscoefdata[i][m][npts];

        }
        xsect[5] -= 1;

        npts += (sm->nquadpts[i][m]+(sm->nwatvap-2)*sm->nquadpts[i][m]);

        if (sm->abscoefdata[i][m][npts] <= 0.0) {
          npts_shift = npts ;
          iw_shift = iw ;
          while (sm->abscoefdata[i][m][npts_shift] <= 0.0 && iw_shift <= 11) {
            ++iw_shift ;
            npts_shift  = (ip-1)*sm->ntemp*sm->nwatvap*sm->nquadpts[i][m]
             +                      (it-1)*sm->nwatvap*sm->nquadpts[i][m]
             +                                (iw_shift-1)*sm->nquadpts[i][m];
          }

          xsect[4] = &sm->abscoefdata[i][m][npts_shift];
        }
        else {

          xsect[4] = &sm->abscoefdata[i][m][npts];
    
        }
        xsect[4] -= 1;

        npts += sm->nquadpts[i][m];

        if (sm->abscoefdata[i][m][npts] <= 0.0) {
          npts_shift = npts ;
          iw_shift = iw ;
          while (sm->abscoefdata[i][m][npts_shift] <= 0.0 && iw_shift <= 11) {
            ++iw_shift ;
            npts_shift  = (ip-1)*sm->ntemp*sm->nwatvap*sm->nquadpts[i][m]
             +                      (it-1)*sm->nwatvap*sm->nquadpts[i][m]
             +                                (iw_shift-1)*sm->nquadpts[i][m];
          }

          xsect[8] = &sm->abscoefdata[i][m][npts_shift];
        }
        else {

          xsect[8] = &sm->abscoefdata[i][m][npts];
        }

        xsect[8] -= 1;

        npts = (ip+1-1)*sm->ntemp*sm->nwatvap*sm->nquadpts[i][m]
             +             (it-1)*sm->nwatvap*sm->nquadpts[i][m]
             +                         (iw-1)*sm->nquadpts[i][m];

        if (sm->abscoefdata[i][m][npts] <= 0.0) {
          npts_shift = npts ;
          iw_shift = iw ;
          while (sm->abscoefdata[i][m][npts_shift] <= 0.0 && iw_shift <= 11) {
            ++iw_shift ;
            npts_shift  = (ip-1)*sm->ntemp*sm->nwatvap*sm->nquadpts[i][m]
             +                      (it-1)*sm->nwatvap*sm->nquadpts[i][m]
             +                                (iw_shift-1)*sm->nquadpts[i][m];
          }

          xsect[2] = &sm->abscoefdata[i][m][npts_shift];
        }
        else {

          xsect[2] = &sm->abscoefdata[i][m][npts];
        }

        xsect[2] -= 1;

        npts += sm->nquadpts[i][m];

        if (sm->abscoefdata[i][m][npts] <= 0.0) {
          npts_shift = npts ;
          iw_shift = iw ;
          while (sm->abscoefdata[i][m][npts_shift] <= 0.0 && iw_shift <= 11) {
            ++iw_shift ;
            npts_shift  = (ip-1)*sm->ntemp*sm->nwatvap*sm->nquadpts[i][m]
             +                      (it-1)*sm->nwatvap*sm->nquadpts[i][m]
             +                                (iw_shift-1)*sm->nquadpts[i][m];
          }

          xsect[6] = &sm->abscoefdata[i][m][npts_shift];
        }
        else {

          xsect[6] = &sm->abscoefdata[i][m][npts];
        }

        xsect[6] -= 1;

        npts += (sm->nquadpts[i][m]+(sm->nwatvap-2)*sm->nquadpts[i][m]);

        if (sm->abscoefdata[i][m][npts] <= 0.0) {
          npts_shift = npts ;
          iw_shift = iw ;
          while (sm->abscoefdata[i][m][npts_shift] <= 0.0 && iw_shift <= 11) {
            ++iw_shift ;
            npts_shift  = (ip-1)*sm->ntemp*sm->nwatvap*sm->nquadpts[i][m]
             +                      (it-1)*sm->nwatvap*sm->nquadpts[i][m]
             +                                (iw_shift-1)*sm->nquadpts[i][m];
          }

          xsect[3] = &sm->abscoefdata[i][m][npts_shift];
        }
        else {

          xsect[3] = &sm->abscoefdata[i][m][npts];
        }
        xsect[3] -= 1;

        npts += sm->nquadpts[i][m];

        if (sm->abscoefdata[i][m][npts] <= 0.0) {
          npts_shift = npts ;
          iw_shift = iw ;
          while (sm->abscoefdata[i][m][npts_shift] <= 0.0 && iw_shift <= 11) {
            ++iw_shift ;
            npts_shift  = (ip-1)*sm->ntemp*sm->nwatvap*sm->nquadpts[i][m]
             +                      (it-1)*sm->nwatvap*sm->nquadpts[i][m]
             +                                (iw_shift-1)*sm->nquadpts[i][m];
          }

          xsect[7] = &sm->abscoefdata[i][m][npts_shift];
        }
        else {

          xsect[7] = &sm->abscoefdata[i][m][npts];
        }
        xsect[7] -= 1;

        /*--------------------------------------------------------*/
        /* Checking interpolation region.                         */
        /*--------------------------------------------------------*/

        for (ii=1; ii<=8; ii++) {
          if (xsect[ii][1] < 0.0) {
            printf("%3d%3d%3d%3d %e\n", ii, ip, it, iw, xsect[ii][1]);
          }
        }

        /*--------------------------------------------------------*/
        /* Do two-dimensional interpolation on pressure and       */
        /* temperature plane.                                     */
        /*--------------------------------------------------------*/

        frcnx = (p-sm->preslog[ip]) / (sm->preslog[ip+1]-sm->preslog[ip]);
        frcny = (t-sm->temp[it]) / (sm->temp[it+1]-sm->temp[it]);

        ansyrl = (1.0-frcnx)*(1.0-frcny)*xsect[1][quadptindex]
               + frcnx*(1.0-frcny)*xsect[2][quadptindex]
               + frcnx*frcny*xsect[3][quadptindex]
               + (1.0-frcnx)*frcny*xsect[4][quadptindex];

        ansyru = (1.0-frcnx)*(1.0-frcny)*xsect[5][quadptindex]
               + frcnx*(1.0-frcny)*xsect[6][quadptindex]
               + frcnx*frcny*xsect[7][quadptindex]
               + (1.0-frcnx)*frcny*xsect[8][quadptindex];

        sm->abscoef[i][m + (j-1)*sm->nummolecules[i]] =
               ((ansyru-ansyrl)*(w-sm->watvaplog[iw]) / 
               (sm->watvaplog[iw+1]-sm->watvaplog[iw])) + ansyrl;

        /*--------------------------------------------------------*/
        /* DONE with the cross section for water moecules.        */
        /*--------------------------------------------------------*/

      }
      else {

        /*--------------------------------------------------------*/
        /* Read in the relevant part of the interpolation table.  */
        /*--------------------------------------------------------*/

        npts = (ip-1)*sm->ntemp*sm->nwatvap*sm->nquadpts[i][m]
             +           (it-1)*sm->nwatvap*sm->nquadpts[i][m]
             +                       (iw-1)*sm->nquadpts[i][m];

        xsect[1] = &sm->abscoefdata[i][m][npts];
        xsect[1] -= 1;

        npts += (sm->nquadpts[i][m]+(sm->nwatvap-1)*sm->nquadpts[i][m]);

        xsect[4] = &sm->abscoefdata[i][m][npts];
        xsect[4] -= 1;

        npts = (ip+1-1)*sm->ntemp*sm->nwatvap*sm->nquadpts[i][m]
             +             (it-1)*sm->nwatvap*sm->nquadpts[i][m]
             +                         (iw-1)*sm->nquadpts[i][m];

        xsect[2] = &sm->abscoefdata[i][m][npts];
        xsect[2] -= 1;

        npts += (sm->nquadpts[i][m]+(sm->nwatvap-1)*sm->nquadpts[i][m]);

        xsect[3] = &sm->abscoefdata[i][m][npts];
        xsect[3] -= 1;

        /*--------------------------------------------------------*/
        /* Checking interpolation region.                         */
        /*--------------------------------------------------------*/

        for (ii=1; ii<=4; ii++) {
          if (xsect[ii][1] < 0.0) {
            printf("%3d%3d%3d%3d %e\n", ii, ip, it, iw, xsect[ii][1]);
          }
        }

        /*--------------------------------------------------------*/
        /* Do two-dimensional interpolation on pressure and       */
        /* temperature plane.                                     */
        /*--------------------------------------------------------*/

        frcnx = (p-sm->preslog[ip]) / (sm->preslog[ip+1]-sm->preslog[ip]);
        frcny = (t-sm->temp[it]) / (sm->temp[it+1]-sm->temp[it]);

        sm->abscoef[i][m + (j-1)*sm->nummolecules[i]] =
                 (1.0-frcnx)*(1.0-frcny)*xsect[1][quadptindex]
               + frcnx*(1.0-frcny)*xsect[2][quadptindex]
               + frcnx*frcny*xsect[3][quadptindex]
               + (1.0-frcnx)*frcny*xsect[4][quadptindex];

        /*--------------------------------------------------------*/
        /* DONE with the cross section for water moecules.        */
        /*--------------------------------------------------------*/

      }

    }

    /*------------------------------------------------------------*/
    /* DONE with this molecule so go to the next one.             */
    /*------------------------------------------------------------*/

  }

  /*--------------------------------------------------------------*/
  /* Based on the spectral model gas constituent, go calculate    */
  /* the extinction coefficient.  The mulitplicate factor of 100. */
  /* is to go from cm^{-1} to m^{-1}.                             */
  /*--------------------------------------------------------------*/

  /* Water Vapor */

  if (!strcmp("WaterVapor",sm->gas[pp->smindex[i]])) {

    bi[0]  = bext_spectralmodelkato(atm->numlayers, atm->densh2o,
                atm->altitude, z, sm->abscoef[pp->smindex[i]][j]);

  }

  /* Carbon Dioxide */

  else if (!strcmp("CarbonDioxide",sm->gas[pp->smindex[i]])) {

    bi[0]  = bext_spectralmodelkato(atm->numlayers, atm->densco2,
                atm->altitude, z, sm->abscoef[pp->smindex[i]][j]);

  }

  /* Ozone */

  else if (!strcmp("Ozone",sm->gas[pp->smindex[i]])) {

    bi[0] = 100.*bext_gasamounts_atmosphere(atm->numlayers, atm->uo3amagats,
                atm->pressure, atm->altitude, z, 0, *ucvd,
                sm->abscoef[pp->smindex[i]][j], 0.);

  }

  /* Oxygen */

  else if (!strcmp("Oxygen",sm->gas[pp->smindex[i]])) {

    bi[0]  = bext_spectralmodelkato(atm->numlayers, atm->denso2,
                atm->altitude, z, sm->abscoef[pp->smindex[i]][j]);

  }

  /* Water Vapor and Carbon Dioxide */

  else if (!strcmp("OverlapH2OCO2",sm->gas[pp->smindex[i]])) {

    bi[0]  = bext_spectralmodelkato(atm->numlayers, atm->densh2o,
                atm->altitude, z, sm->abscoef[pp->smindex[i]][2*j-1]);

    bi[0] += bext_spectralmodelkato(atm->numlayers, atm->densco2,
                atm->altitude, z, sm->abscoef[pp->smindex[i]][2*j]);

  }

  /* Water Vapor and Oxygen and Ozone */

  else if (!strcmp("OverlapH2OO2O3",sm->gas[pp->smindex[i]])) {

    bi[0]  = bext_spectralmodelkato(atm->numlayers, atm->densh2o,
                atm->altitude, z, sm->abscoef[pp->smindex[i]][3*j-2]);

    bi[0] += bext_spectralmodelkato(atm->numlayers, atm->denso2,
                atm->altitude, z, sm->abscoef[pp->smindex[i]][3*j-1]);

    bi[0] += 100.*bext_gasamounts_atmosphere(atm->numlayers, atm->uo3amagats,
                atm->pressure, atm->altitude, z, 0, *ucvd,
                sm->abscoef[pp->smindex[i]][3*j], 0.);

  }

  /* Water Vapor and Ozone */

  else if (!strcmp("OverlapH2OO3",sm->gas[pp->smindex[i]])) {

    bi[0]  = bext_spectralmodelkato(atm->numlayers, atm->densh2o,
                atm->altitude, z, sm->abscoef[pp->smindex[i]][2*j-1]);

    bi[0] += 100.*bext_gasamounts_atmosphere(atm->numlayers, atm->uo3amagats,
                atm->pressure, atm->altitude, z, 0, *ucvd,
                sm->abscoef[pp->smindex[i]][2*j], 0.);
   
    /* printf("%e %e %e from kato after watervapor ozone bi\n", *atm->densh2o, z, bi[0]) ; */
  }

  /* No Gaseous Extinction */

  else if (!strcmp("None",sm->gas[pp->smindex[i]])) {

    bi[0] = 0.;

  }

  /*--------------------------------------------------------------*/
  /* Make the assignments to the remaining coefficients.          */
  /*--------------------------------------------------------------*/

  bi[1] = 0.;
  bi[2] = bi[0];

  /*--------------------------------------------------------------*/
  /* DONE                                                         */
  /*--------------------------------------------------------------*/

}

/**********************************************************************/

double
bext_spectralmodelkato(numlayers,dens,height,altitude,km)
int numlayers;
double *dens;
double *height;
double altitude;
double km;
{
  double
    bext,
    densmidlayer;

  int i;

  /* densmidlayer = linear_interpolate(dens, height, numlayers, altitude); */
  densmidlayer = log_interpolate(dens, height, numlayers, altitude);
  bext = (km/10000.)*densmidlayer;

  return ((double)bext);

}

/**********************************************************************/
/**********************************************************************/
