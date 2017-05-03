      subroutine rap_twostr_vir
     &  (solconst,nlayers,amu0,w0,g0,suralb,taulam,
     &   gau_wt,dis2sun)


c
c  driver for 2-stream multiple scattering calculation.
c
c All variables are explained in fourstr_vasriables.dat

c Inputs
c dis2sun      The distance to the sun in astronomical units.
c solconst     The direct beam irradiance at the TOA
c nlayers      Total number of layers in the atmosphere
c amu0         Cosine of the solar zenith angle (Radian)
c w0(j)        Single scattering albedo for jth layer (j=1,nlayers)
c g0(j)        (l=0,4) (j=1,nalyers)
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

      real w0(*), g0(*), taulam(*)
      real amu0, surfalb
      real*8 gau_wt

c radcorr corrects the solar constant based on the ephemeris derived
c distance to sun 

c      radcorr = dis2sun*dis2sun

c w0max is the maximum allowed asymmetry parameter - values larger 
c than 0.9999 tend to induce model instabilities.

      w0max = 0.9999

c making the delta approximation

      do 10 j =1, nlayer
        if(j.ne.1) then
          j1=j-1
        else
          j1=j
        endif

        fo        = g0(j)**2
        den       = 1.-w0(j)*fo
        taul(j)   = taul(j) * den
        taulnd(j) = taul(j)
        w0(j)     = (1.-fo)*w0(j)/den
        g0(j)     = g0(j)/(1.+g0(j))
        opd(j)    = 0.0
        opd(j)    = opd(j1) + taul(j)
        opdnd(j)  = opdnd(j1) + taulnd(j)

        if(w0(j).gt.w0max) w0(j) = w0max

 10   continue

c for IR

      if(irs .ne. 0) then

        nlow = 1250
        nhigh = 3250
        ncount = nhigh - nlow

        v         =   1.438e4  /  wave
       
        do 370 j  =   1,ncount
          jj     =   nlow+j
          t1     =   0.1 * float(jj)

          call plnk(v,t1,plank(j))

c 4.546 and 62.5 are in micron

          if (wave.le.4.546) then
            plank(j) = 0.0
          elseif (wave.ge.62.5) then
            plank(j) = t1**4
          endif

 370    continue

      
        do 15 j = 1, nalyer
          call plnk(v,t1,plank(j))
 15     continue

        call oppr1

      endif


      do 20 j = 1,ntaujs
         fup(j) = 0.0
         fdn(j) = 0.0
 20   continue


       call twostr
     & (nlayer,taul,w0,g0,rsfx,b1,b2,el1,el2,em1,em2,af,bf,ef)

       call add_vir
     &  (nlayer,taul,w0,g0,rsfx,opd,opdnd,ak,b1,b2,b3,em1,em2,
     &   el1,el2,af,bf,ef,fnet,sol_fluxup,sol_fluxdn,direct)


c adjust top of model domain flux to account for middle atmosphere
c extinction according to two-stream results

c         toaflx = (ck1(i,1)*el2(i,1)+ck2(i,1)*em2(i,1)+
c     1               cmb(i,1)+direct(i,1))/amu0

c correction factor of the sun to the earth distance

         dist_cor = 1.0 / (dis2sun)**2


c obtain irradiances at the bottom of jth layer

      if (iopen .eq. 0) then

        open(unit=41,file='out.foursrt.dat',status='unknown')

      endif

      iopen = 1

      do 50 j = 1,ntaujs
         sol_fluxup(j) = fup(j)*solconst*gau_wt*dist_cor
         sol_fluxdn(j) = fdn(j)*solconst*gau_wt*dist_cor
         sol_direct(j) = direct(j)*solconst*gau_wt*dist_cor

         write(41,400) j, sol_fluxup(j), sol_fluxdn(j), sol_direct(j)

 50   continue


c      write(40,240) amu0, dis2sun
c      write(40,200)

 200  format('layer z(km)  prs(mb)   sw down    sw up')
 220  format(i4,2x,f6.2,2x,f7.1,2x,f9.3,3x,f9.3)
 240  format('mu_0=',f6.4,5x,'radau=',f6.4)

  400 format(i5,3e15.5)

      return
      end

 
