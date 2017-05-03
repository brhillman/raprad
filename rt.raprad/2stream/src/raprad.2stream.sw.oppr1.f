      subroutine oppr1
     &  (tgrand,nlow,plank

c last revison:  25 jan 1995


c     **********************************************************
c     *  purpose             :  calculate planck function and  *
c     *                         and its derivative at ground   *
c     *                         and at all altitudes.          *
c     *  subroutines called  :  none                           *
c     *  input               :  tgrnd, nlow, weight            *
c     *  output              :  ptemp, ptempg, slope           *
c     * ********************************************************

c define the dimensions used by the radiation model but might be
c specified by an external model

c ivert  = maximum number of layers;
c ilayer = maximum number of layer boundaries
c idbl   = twice the maximum number of layer boundaries
c irad   = maximum number of aerosol radius bins;

      parameter (ivert=200)
      parameter (irad=20)
      parameter (ilayer=ivert+1, idbl=2*ilayer)
      parameter (iradver=irad*ivert)
      parameter (iradlay=irad*ilayer)

c define the dimensions used only by the radiation model
c
c isol   = total number of solar wavelengths
c iir    = total numberof infrared wavelengths
c iwave  = isol + iir

      parameter (isol=32, iir=18, iwave=isol+iir)

c itotal = total number of probability intervals;
c isolp  = number of solar probability intervals;
c iirp   = number of infrared probability intervals;

c  double precision

      implicit real*8 (a-h, o-z)


      dimension  ptemp2(itotal), pltemp1(itotal)

c     **************************************
c     * calculate ptemp and slope          *
c     **************************************
c

c       calculate the wavelength dependent plank function at the ground.

        itg                 = anint(10.*tgrnd) - nlow

        do 93 i=1,nirp
         pltemp1(i) = plank(ltemp(i),itg)
93      continue

         do 100 l            =   nsolp+1,ntotal
            ptempg(l)        =   pltemp1(l-nsolp)*weight(l)
 100    continue


        do 300 j            =   1,nlayer
         if(j.eq.1) then
          kindex = 1
         else
          kindex = j-1
         endif
            it1             = anint(10.*tt(j)) - nlow
        do 103 i = 1,nirp
         ptemp2(i) = plank(ltemp(i),it1)
103     continue

c           kindex makes the top layer isothermal. using kindex, find
c           plank function at bottom of each layer.
c           note: if you force slope=0, then you have isothermal
c           layers with tt(j) corresponding to average temperature
c           of layer and tt(nlayer) should be set to tgrnd.

            do 200 l        = nsolp+1,ntotal
               ptemp(l,j)   = ptemp2(l-nsolp)*weight(l)
               slope(l,j)   = (ptemp(l,j)-ptemp(l,kindex))/taul(l,j)
               if(taul(l,j).le.1.0e-6)   slope(l,j)=0.0
 200        continue
 300     continue

      return
      end
