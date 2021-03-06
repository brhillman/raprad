      subroutine fourstr
     &  (solconst,nlayers,amu0,omg0in,pmomts,surfalb,
     &   taulam,gau_wt,dis2sun)

c need to add radcorr

c
c  driver for 4-stream multiple scattering calculation.
c
c All variables are explained in fourstr_vasriables.dat

c Inputs
c dis2sun      The distance to the sun in astronomical units.
c solconst     The direct beam irradiance at the TOA
c nlayers      Total number of layers in the atmosphere
c amu0         Cosine of the solar zenith angle (Radian)
c omg0in(j)    Single scattering albedo for jth layer (j=1,nlayers)
c pmomts(lx(l),j)  (l=0,4) (j=1,nalyers)
c              Coefficients of Legendre polynominals 
c              pmomts(0)=1, and pmomts(1)=3*g
c              pmomts(4) is used for the delta approximation
c surfalb      Surface albedo
c taulam(j)    Optical depth of jth layer (j=1,nlayers)
c iprint       print output flug =0, no outputs

c Output variables

c sol_fluxdn(j)  Downward shortwave irradinace at j level 
c              (W m^-2) (j=1,ntaujs)
c              j level is the top boundary of jth layer
c              the solar constant at TOA is 1.0 so to obtain W m^-2, need to
c              multiply by the solar constant at TOA
c sol_fluxup(j)  Upwnward shortwave irradinace at j level 
c              (W m^-2) (j=1,ntaujs)
c sol_direct(j)  Direct beam irradiance (W m^-2) (j=1,ntaujs)


c ivert  = maximum number of layers;

      parameter (ivert=200)
      parameter (ilayer=ivert+1)                       


      implicit double precision(a-h, o-z)

      dimension fup(ilayer), fdn(ilayer), direct(ilayer)
      dimension sol_fluxup(ilayer), sol_fluxdn(ilayer)
      dimension sol_direct(ilayer)

      real omg0in(*), pmomts(5,*), taulam(*) 
      real amu0, surfalb
      real*8 gau_wt

c test---------------------------------------------------

      if (iopen .eq. 0) then

        open(unit=40,file='out.diagnostic.fourstr.dat',
     &     status='unknown')

        open(unit=41,file=
     &    '../../results/raprad.sw.out',status='unknown')

      endif

      iopen = 1

c     write(*,500) solconst
c     write(*,501) nlayers
c     write(*,500) amu0
c     write(*,502) (omg0in(i), i=1,3)
c     write(*,503) (pmomts(j,1), j=1,5)
c     write(*,503) (pmomts(j,2), j=1,5)
c     write(*,500) surfalb
c     write(*,502) (taulam(i), i=1,6)

  500 format(f15.5)
  501 format(i5)
  502 format(3f15.5)
  503 format(5f15.5)

      write(40,301) solconst, nlayers, amu0, surfalb

      do 60 j = 1, nlayers
        write(40,300) j, omg0in(j),taulam(j),(pmomts(i,j), i=1,5)
   60 continue


c-------------------------------------------------------


c radcorr corrects the solar constant based on the ephemeris derived
c distance to sun 

c      radcorr = dis2sun*dis2sun

c w0max is the maximum allowed asymmetry parameter - values larger 
c than 0.9999 tend to induce model instabilities.

      w0max = 0.9999
      ntaujs = nlayers + 1

      do 10 j = 1,ntaujs
         fup(j) = 0.0
         fdn(j) = 0.0
 10   continue

       do 30 j = 1,nlayers
         if(omg0in(j).gt.w0max) omg0in(j) = w0max
 30    continue

       call solrad
     &  (nlayers,amu0,omg0in,pmomts,surfalb,taulam,
     &     iprint,fdn,fup,direct)


c adjust top of model domain flux to account for middle atmosphere
c extinction according to two-stream results

c         toaflx = (ck1(i,1)*el2(i,1)+ck2(i,1)*em2(i,1)+
c     1               cmb(i,1)+direct(i,1))/amu0

c correction factor of the sun to the earth distance

         dist_cor = 1.0 / (dis2sun)**2

         do 50 j = 1,ntaujs
            sol_fluxup(j) = fup(j)*solconst*gau_wt*dist_cor
            sol_fluxdn(j) = fdn(j)*solconst*gau_wt*dist_cor
            sol_direct(j) = direct(j)*solconst*gau_wt*dist_cor

            write(41,400) j, sol_fluxup(j), sol_fluxdn(j), sol_direct(j)

 50      continue


  300 format(i5, 7e15.5)
  301 format(f10.5, i5, 2f10.5)
  400 format(i5,3e15.5)

      return
      end

 
