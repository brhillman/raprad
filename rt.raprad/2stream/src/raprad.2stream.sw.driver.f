      subroutine rap_twostr_sw
     &  (solconst,nlayers,amu0,w0_temp,g0_temp,suralb,
     &   taul_temp,gau_wt,dis2sun)


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
c taul(j)    Optical depth of jth layer (j=1,nlayers)
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

!     implicit double precision(a-h, o-z)
      implicit none

      integer ivert, ilayer, idbl
      parameter (ivert=200)
      parameter (ilayer=ivert+1)  
      parameter (idbl=ilayer*2)                       
                     
!     inputs
      double precision solconst
      integer nlayers
      real amu0
      real w0_temp(*), g0_temp(*), suralb, taul_temp(*)
      double precision gau_wt
      double precision dis2sun

!     local namespace
      double precision fup(ilayer), fdn(ilayer), fdi(ilayer)
      double precision fnet(ilayer)
      double precision sol_fluxup(ilayer), sol_fluxdn(ilayer)
      double precision sol_direct(ilayer)

      double precision b1(ilayer), b2(ilayer), b3(ilayer)
      double precision el1(ilayer), el2(ilayer)
      double precision em1(ilayer), em2(ilayer)
      double precision af(idbl), bf(idbl), ef(idbl)
      double precision ak(ilayer)

      double precision den, denom, dist_cor, fo, g0t

      real w0(ilayer), g0(ilayer)
      real taul(ilayer), g0l(ilayer), w0l(ilayer)
      real taulnd(ivert), opd(ilayer), opdnd(ilayer)

      integer j, j1, ntaujs
      real w0max,w0t
      double precision epsilon
      data epsilon / 1.0e-15 /

c radcorr corrects the solar constant based on the ephemeris derived
c distance to sun 

c      radcorr = dis2sun*dis2sun

c w0max is the maximum allowed asymmetry parameter - values larger 
c than 0.9999 tend to induce model instabilities.

      w0max = 0.9999
      ntaujs = nlayers + 1

      do 1000 j = 1, ntaujs

        if(j .eq. 1) then
          taul(j) = 0.0
          w0l(j)   = 0.0
          g0l(j)   = 0.0
        else
          taul(j) = taul_temp(j-1)
          w0l(j)   = w0_temp(j-1)
          g0l(j)   = g0_temp(j-1)
        end if

 1000 continue


c making the delta approximation

      do 10 j =1, ntaujs

        if(j.ne.1) then
          j1=j-1
        else
          j1=j
        endif

        if(taul(j).lt.epsilon)  taul(j) = epsilon
           w0t = w0l(j)

        if(w0t.gt.1.-epsilon)   w0t=1.-epsilon
           denom = w0l(j) * taul(j)
        if(denom.le.epsilon)  denom=epsilon

        if(denom.gt.epsilon) then
          g0t = g0l(j)
        else
          g0t = 0.0
        endif


        fo        = g0t**2

        den       = 1.-w0t*fo



        taulnd(j) = taul(j)
        taul(j)   = taul(j) * den
        w0(j)     = (1.-fo)*w0t/den



        g0(j)     = g0t/(1.+g0t)
        opd(j)    = 0.0
        opd(j)    = opd(j1) + taul(j)
        opdnd(j)  = opdnd(j1) + taulnd(j)

        if(w0(j).gt.w0max) w0(j) = w0max

 10   continue


      do 20 j = 1,ntaujs
         fup(j) = 0.0
         fdn(j) = 0.0
 20   continue



       if (amu0 .gt. 0.0) then
         call twostr 
     &   (ntaujs,taul,w0,g0,suralb,b1,b2,el1,el2,em1,em2,af,bf,ef, 
     &    ak)

         call add 
     &    (ntaujs,taul,w0,g0,suralb,opd,opdnd,ak,b1,b2,b3,em1,em2, 
     &     el1,el2,af,bf,ef,amu0,fnet,fup,fdn,fdi)

       end if


c correction factor of the sun to the earth distance

         dist_cor = 1.0 / (dis2sun)**2


c obtain irradiances at the bottom of jth layer


        open
     &  (unit=41,file='../../results/raprad.sw.out',status='unknown')

!solconst changes with layer
!gau_wt=1.
!dist_cor=1. 

      do 50 j = 1,ntaujs

         sol_fluxup(j) = fup(j)*solconst*gau_wt*dist_cor
         sol_fluxdn(j) = fdn(j)*solconst*gau_wt*dist_cor
         sol_direct(j) = fdi(j)*solconst*gau_wt*dist_cor*amu0

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

 
