      subroutine oppr1
     &  (nlayers,ptempg,taulam,slope,ptemp,bandwidth,
     &   plankbnd,gau_wt,planklev)


c     **********************************************************
c     *  purpose             :  calculate planck function and  *
c     *                         and its derivative at ground   *
c     *                         and at all altitudes.          *
c     *  subroutines called  :  none                           *
c     *  input               :  nlayers,taulam,gau_wt,  *
c      *                        plankbnd,planklev              *
c     *  output              :  ptemp, ptempg, slope           *
c     * ********************************************************

c This routine computes gaussian weighted plank function and slope,
c which is defined by Eq. 26 in Toon et al.
c The layer number is 1 at the top of the atmosphere and increases with
c decreasing the height. The plank function is previously computed by 
c Mlawer's code.

      parameter (ivert=200)
      parameter (ilayer=ivert+1, idbl=2*ilayer)

      implicit double precision(a-h, o-z) 


      integer nlayers
      real*8  ptempg
      real    taulam(ilayer)

      dimension slope(ilayer), ptemp(ilayer), bandwidth(2)
      real*8  plankbnd
      dimension gau_wt(*), planklev(*)

cc define the dimensions used by the radiation model but might be
cc specified by an external model
c
cc ivert  = maximum number of layers;
cc ilayer = maximum number of layer boundaries
cc idbl   = twice the maximum number of layer boundaries
cc irad   = maximum number of aerosol radius bins;
c
c      parameter (ivert=200)
c      parameter (irad=20)
c      parameter (ilayer=ivert+1, idbl=2*ilayer)
c      parameter (iradver=irad*ivert)
c      parameter (iradlay=irad*ilayer)
c
cc define the dimensions used only by the radiation model
cc
cc isol   = total number of solar wavelengths
cc iir    = total numberof infrared wavelengths
cc iwave  = isol + iir
c
c      parameter (isol=32, iir=18, iwave=isol+iir)
c
cc itotal = total number of probability intervals;
cc isolp  = number of solar probability intervals;
cc iirp   = number of infrared probability intervals;
c
cc  double precision
c       implicit real*8 (a-h, o-z)
c
c
c      dimension  ptemp2(itotal), pltemp1(itotal)

      deltainvcm = (1.0/bandwidth(1) - 1.0/bandwidth(2)) * 1.0d4

       
c
cc     **************************************
cc     * calculate ptemp and slope          *
cc     **************************************
cc
cc       calculate the wavelength dependent plank function at the ground.
c        itg                 = anint(10.*tgrnd) - nlow
c
c        do 93 i=1,nirp
c         pltemp1(i) = plank(ltemp(i),itg)
c93      continue
c         do 100 l            =   nsolp+1,ntotal
c            ptempg(l)        =   pltemp1(l-nsolp)*weight(l)

c 1.0e4 factor comes in to convert units from W cm^-2 to W m^-2
            ptempg = plankbnd * gau_wt(nlayers) * 1.0e4 * deltainvcm
c
c 100    continue
c
c 1 is the top of the stmosphere and nlayer is the bottom of the atmosphere.
        do 300 j            =   1,nlayers
         if(j.eq.1) then
          kindex = 1
         else
          kindex = j-1
         endif
c
c            it1             = anint(10.*tt(j)) - nlow
c        do 103 i = 1,nirp
c         ptemp2(i) = plank(ltemp(i),it1)
c103     continue

c           kindex makes the top layer isothermal. using kindex, find
c           plank function at bottom of each layer.
c           note: if you force slope=0, then you have isothermal
c           layers with tt(j) corresponding to average temperature
c           of layer and tt(nlayer) should be set to tgrnd.

c            do 200 l        = nsolp+1,ntotal
c               ptemp(j)   = ptemp2*gau_wt
c               slope(j)   = (ptemp(j)-ptemp(kindex))/taulam(j)

c 1.0e4 factor comes in to convert units from W cm^-2 to W m^-2
	       tmpvalue = 1.0e4 * planklev(j) * gau_wt(j) * deltainvcm
	       ptemp(j) = tmpvalue
c              ptemp(j)    = 1.0e4 * planklev(j) * gau_wt(j) * deltainvcm
               slope(j)    = (ptemp(j) - ptemp(kindex))/taulam(j)
c     &           (gau_wt*1.0e4*(planklev(j)-planklev(kindex)))/taulam(j)
               if(taulam(j).le.1.0e-6)   slope(j)=0.0
c 200        continue
 300     continue

      return
      end
