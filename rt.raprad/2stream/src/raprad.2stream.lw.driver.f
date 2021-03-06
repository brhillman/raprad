      subroutine rap_twostr_lw
     &  (solconst,nlayers,amu0,w0_temp,g0_temp,suralb,taul_temp,
     &   gau_wt_temp,dis2sun,plankbnd,planklay_temp,planklev_temp,
     &   bandwidth)


c
c  driver for 2-stream multiple scattering calculation.
c
c All variables are explained in fourstr_vasriables.dat

c Inputs
c solconst     The direct beam irradiance at the TOA
c nlayers      Total number of layers in the atmosphere
c amu0         Cosine of the solar zenith angle (Radian)
c w0(j)        Single scattering albedo for jth layer (j=1,nlayers)
c g0(j)        (l=0,4) (j=1,nalyers)
c              Coefficients of Legendre polynominals 
c              pmomts(0)=1, and pmomts(1)=3*g
c              pmomts(4) is used for the delta approximation
c suralb       Surface albedo in the shortwave region
c              surface emissivity is 1.0 - suralb
c taul(j)      Optical depth of jth layer (j=1,nlayers)
c gau_wt       Gaussian weight to sum up sub-interval weighted by
c              Planck functions
c dis2sun      The distance to the sun in astronomical units.
c plankbnd     Plank function integrated over the width of the band 
c              in the units of W cm^-2/ pi at the surface temperature
c planklay     Plank function integrated over the width of the band 
c              in the units of W cm^-2/ pi at the layer temperature
c planklev     Plank function integrated over the width of the band 
c              in the units of W cm^-2/ pi at the level temperature

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

      parameter (irflag = 1)
      parameter (ivert=200)
      parameter (ilayer=ivert+1, idbl=2*ilayer)
      parameter (ngauss = 3)


      implicit double precision(a-h, o-z)

      dimension fup(ilayer), fdn(ilayer), direct(ilayer)
      dimension fnet(ilayer)
      dimension sol_direct(ilayer)
      dimension planklay_temp(*), planklev_temp(*), bandwidth(2)
      dimension gau_wt_temp(*)

      dimension b1(ilayer), b2(ilayer), b3(ilayer)
      dimension el1(ilayer), el2(ilayer)
      dimension em1(ilayer), em2(ilayer), ee1(ilayer)
      dimension af(idbl), bf(idbl), ef(idbl)
      dimension ak(ilayer), ck1(ilayer), ck2(ilayer)
      dimension slope(ilayer), ptemp(ilayer)
      dimension direc(ilayer), directu(ilayer)
      dimension gangle(ngauss), gweight(ngauss)
      dimension y3(ngauss,ilayer), gami(ilayer)
      dimension fluxird(ilayer), fluxiru(ilayer)

      real w0_temp(*), g0_temp(*), taul_temp(*)
      real w0(ilayer), g0(ilayer)
      real taul(ilayer), g0l(ilayer), w0l(ilayer)
      real taulnd(ivert), opd(ilayer), opdnd(ilayer)
      real amu0, suralb
      real*8 planklay(ilayer), planklev(ilayer)
      real*8 gau_wt(ilayer)


      data epsilon / 1.0e-15 /

c     gauss angles and gauss weights for gaussian integration
c     moments (use first moment values) n=3
c
      data gangle  /  0.2123405382, 0.5905331356,
     1                0.9114120405                       /
      data gweight /  0.0698269799, 0.2292411064,
     1                0.2009319137                       /


c radcorr corrects the solar constant based on the ephemeris derived
c distance to sun 

c      radcorr = dis2sun*dis2sun

c w0max is the maximum allowed asymmetry parameter - values larger 
c than 0.9999 tend to induce model instabilities.

      w0max = 0.9999
      ntaujs = nlayers + 1

c adding one layer on top the atmosphere. Then level 0 becomes level 1
c and layer 1 becomes layer 1.
c
c                        --------
c                           1
c   0---------          1--------
c        1                  2
c   1---------          2--------
c        2                  3
c   2---------    ->    3--------
c        .
c        .
c n-1---------          n--------
c        n                 n+1
c   n---------        n+1--------
c
      do 1000 j = 1, ntaujs
        if(j .eq. 1) then
          taul(j) = 0.0
          w0l(j)   = 0.0
          g0l(j)   = 0.0
          planklev(j) = 0.0
          planklay(j) = 0.0
          gau_wt(j) = 0.0
        else
          taul(j) = taul_temp(j-1)
          w0l(j)   = w0_temp(j-1)
          g0l(j)   = g0_temp(j-1)
          planklev(j) = planklev_temp(j-1)
          planklay(j) = planklay_temp(j-1)
          gau_wt(j) = gau_wt_temp(j-1)
        end if
 1000 continue


c making the delta approximation

      do 1100 j =1, ntaujs

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

        do 1200 i = 1, ngauss
          x         =   taul(j)/gangle(i)
          if(x.gt.1000.)  x=1000.
          if(x.lt.1000.) then
            y3(i,j) = exp(-x)
          else
            y3(i,j) = 0.0
          endif
 1200   continue

 1100 continue

      
      do 2000 j = 1,ntaujs
         fup(j) = 0.0
         fdn(j) = 0.0
         direct(j) = 0.0
         direc(j)  = 0.0
         directu(j) = 0.0
 2000   continue

      if (iopen .eq. 0) then

        open
     &  (unit=41,file='../../results/raprad.lw.out',status='unknown')
c       open(unit=41,file='out.twosrt.dat',status='unknown')

      endif

      iopen = 1


c calculating weighted plank functions and slope defined by Eq. 26 in 
c Toon et al.

      call oppr1
     & (ntaujs,ptempg,taul,slope,ptemp,bandwidth,
     &  plankbnd,gau_wt,planklev)

      call twostr
     & (ntaujs,irflag,taul,w0,g0,suralb,
     &  b1,b2,el1,el2,em1,em2,af,bf,ef,ak,u1i,u1s,gami,ee1)

      call add
     & (ntaujs,taul,w0,g0,suralb,opd,opdnd,ak,b1,b2,b3,em1,em2,
     &  el1,el2,af,bf,ef,amu0,slope,ptempg,ptemp,u1i,u1s,
     &  fnet,fup,fdn,direct,
     &  irflag,ck1,ck2)

      call newflux1
     & (ntaujs,taul,ptempg,ptemp,slope,b3,
     &  gami,ak,gangle,gweight,y3,ee1,
     &  ck1,ck2,suralb,u1i,u1s,
     &  direc,directu)

      do 3000 j = 1,ntaujs

        fluxird(j) = direc(j)
        fluxiru(j) = directu(j)
        sol_direct(j) = 0.0

c        write(*, 400) j, fluxiru(j), fluxird(j), sol_direct(j)
        write(41,400) j, fluxiru(j), fluxird(j), sol_direct(j)

 3000 continue

c      write(40,240) amu0, dis2sun
c      write(40,200)

 200  format('layer z(km)  prs(mb)   sw down    sw up')
 220  format(i4,2x,f6.2,2x,f7.1,2x,f9.3,3x,f9.3)
 240  format('mu_0=',f6.4,5x,'radau=',f6.4)

  400 format(i5,3e15.5)

      return
      end
