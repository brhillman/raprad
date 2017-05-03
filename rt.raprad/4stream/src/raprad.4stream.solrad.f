      subroutine solrad
     &  (nlayers,amu0,omg0in,pmomts,surfalb,taulam,
     &   iprint,fdn,fup,direct)
 
c ldelta fixed (otherwise, successive solrad calls can be
c   in error if omeg0 not reset for each call) (amv, 1993)
c
c plane-parallel version using block tri-diagonal matrix solver
c
c calculates solar fluxes and heating rates using 4s or delta-4s method
c monochromatic wavelength.
c
c  ivert   maximum number of layers allowed
c  ilayer   maximum number of tau levels allowed (= ivert+1)
c  nwrkmax   dimension of input matrix, column vector, and work space
c            array for banded matrix solver; must be  = 4 * ivert

c Inputs
c nlayers      Total number of layers in the atmosphere
c amu0         Cosine of the solar zenith angle (Radian)
c omg0in(j)    Single scattering albedo for jth layer (j=1,nlayers)
c omega(lx(l),j)  (l=0,4) (j=1,nalyers)
c              Coefficients of Legendre polynominals 
c              omega(0)=1, and omega(1)=3*g
c              omega(4) is used for the delta approximation
c surfalb      Surface albedo
c taulam(j)    Optical depth of jth layer (j=1,nlayers)
c iprint       print output flug =0, no outputs

c Outputs
c fdn(j)       Downward shortwave irradinace at j level (Unit less) (j=1,ntaujs)
c              j level is the top boundary of jth layer
c              the solar constant at TOA is 1.0 so to obtain W m^-2, need to
c              multiply by the solar constant at TOA
c fup(j)       Upwnward shortwave irradinace at j level (Unit less) (j=1,ntaujs)
c direct(j)    Direct beam irradiance (Unit less) (j=1,ntaujs)


      implicit double precision( a-h, o-z )
c                                                                     
c define the dimensions used by the radiation model but might be      
c specified by an external model                                      
c                                                                     
c ivert  = maximum number of layers;                                  
c ilayer = maximum number of layer boundaries                         
c idbl   = twice the maximum number of layer boundaries               
c irad   = maximum number of aerosol radius bins;                     
c                                                                     
      parameter (ivert=200)                                            
      parameter (ilayer=ivert+1)                       

c      dimension omg0in(ivert), taulam(ivert)
      dimension fup(ilayer), fdn(ilayer), direct(ilayer)
      dimension omega(5)
c     dimension pmomts(5,ivert)
      dimension a(4,4,ivert), b(4,4,ivert), c(4,4,ivert), d(4,ivert)
c      dimension gauqdt(4), gauwt(4)
      dimension rmu(4), pleg(4,4), efun(2,4,ivert), eval(2,ivert)
      dimension tau(0:ivert), tcor(ivert)
      dimension taund(ivert), omg0(ivert)
      dimension psi(4), capz(4,ivert), smallb(2), deltai(2)
      dimension p(2,ivert), q(2,ivert), expktn(2), expktm(2)
      dimension conp(2), conm(2), ri_qdt(4)

      real omg0in(*), pmomts(5,*), taulam(*) 
      real amu0, surfalb

      integer  nlayers
      integer  leqerr, iprint, layer

c index function for moments and legendre polynomials

      lx( l )  =  l + 1

c define general constants

      pie     =  acos( -1.0 )
      twopie  =  2.0 * pie
      forpie  =  4.0 * pie


c zero out matrix arrays

 
      leqerr = 0
 

c initialize arrays for 4-stream matrix calculation

      do 10 n = 1, ivert
         do 20 i = 1, 4
             d(i,n) = 0.0
             do 25 j= 1,4
               a(i,j,n) = 0.0
               b(i,j,n) = 0.0
               c(i,j,n) = 0.0
25           continue             
20       continue
10    continue


      f0      =  1.0

      amu0i  =  1.0 / amu0
      amu0sq =  amu0 * amu0


c obtaining Gaussian quadratures and weights

c      call gauleg(-1.0d0, 1.0d0, gauqdt, gauwt, 4)

      rmu1 = 0.3399810435848563
      rmu2 = 0.8611363115940526
      rmu1sq  =  rmu1 * rmu1
      rmu2sq  =  rmu2 * rmu2

      a1 = 0.6521451548625464
      a2 = 0.3478548451374476

      swmu  =  a1 * rmu1  +  a2 * rmu2

c energy conservation in the model

      suralb  =  surfalb / swmu

c setting up the legendre polynomial array
 
      call ingaus(rmu1,rmu2,pleg,rmu)

c These are related to Z function defined by equation 27 in Liou (1974)
c Z function is divided by 4*pi so pi appears on the denominetor. 
      
      anum   =  ( rmu1sq - amu0sq ) * ( rmu2sq - amu0sq )
      dnom   =  forpie * rmu1sq * rmu2sq
      hhcoef =  anum / dnom
 

      summ   =  0.0
      summnd =  0.0
      abstst =  0.0
      tau(0) =  0.0
 

      if ( iprint .gt.0 )  write( 6,910 )
 
 
      if ( iprint .gt.0  ) then
         if ( ldelta ) then
            write( 6,920 )
         else
            write( 6,922 )
         endif
      endif

 
c start loop for layers -------------------------------------------------

c This next block is part of the 'Z' term that is not
c dependent on the wavelength.

      do 30 layer = 1, nlayers

         do 40 l = 0, 4
            omega( lx(l) )  =  pmomts( lx(l),layer )
 40      continue

 
c define normalized phase function moments -- array omega
c   note: omg0 is single-scatter albedo, omega(1) = 1.0
c if delta switch is on, correct omgl and tau.
 
c making the delta approximation (Cuzzi et la. 1982)
c following equation (3.23), (3.24), and (3.25) of Liou et al., (1998)

        fmom  =  pmomts( lx(4),layer ) / 9.0
        tcor( layer )  =  1.0 - omg0in( layer ) * fmom
        cofm  =  1.0 - fmom
 
        omg0( layer ) = cofm * omg0in( layer ) / tcor( layer)

c omega(lx(0)) =1.0 is unchanged

        do 50 l = 1, 3
          tlp1  =  2.0 * float( l )  +  1.0
          omega( lx(l) ) = ( omega( lx(l) ) - tlp1* fmom ) / cofm
50      continue
 
c keep no delta approximation in order to compute the direct beam irradiance

        tau( layer ) =  summ + tcor( layer)* taulam( layer )
        taund( layer ) = summnd + taulam( layer )
        abstst = abstst+tcor(layer)*taulam(layer)*(1.0-omg0in(layer))
        summ   = tau( layer )
        summnd = taund( layer)
 
        if ( iprint .gt.0 )
     &    write(6,875) layer, omg0(layer), ( omega(l),l=2,5 )
 

c unnormalize moments by current value of omg0

         do 60 l = 0, 3
            omega( lx(l) )  =  omg0( layer ) * omega( lx(l) )
 60      continue

c obtaining eiganvalue and eiganfunction given by equation 31 and 26,
c respectively, in Liou (1974)

         call calceig(omega,rmu,pleg,a1,a2,eval,efun,layer)
 

c compute  capz function, which is given by equation 27 in Liou (1974)
c X varies as rmu(i)
c capz = Z / (pi*F0)

         v1sq  =  eval( 1,layer ) * eval( 1,layer ) * amu0sq
         v2sq  =  eval( 2,layer ) * eval( 2,layer ) * amu0sq
         prod  =  ( 1.0 - v1sq ) * ( 1.0 - v2sq )


c define co-single-scattering albedo.
c These are also defined in calceig.

      omg00 =  omega( lx(0) )
      omg1  =  omega( lx(1) )
      omg2  =  omega( lx(2) )
      omg3  =  omega( lx(3) )
 
      coomg0  =  1.0 - omg00
      coomg1  =  3.0 - omg1
      coomg2  =  5.0 - omg2



         call psicalc( amu0i, coomg0, coomg1, coomg2, psi )
 
         do 100 i = 1, 4
            sum  =  0.0
            do 70 l = 0, 3
               sum = sum + omega( lx(l) ) * psi( lx(l) )*pleg( lx(l),i )
 70         continue
            dnom     =  prod * ( amu0 + rmu(i) )
            capz( i,layer )  =  hhcoef * sum / dnom
 100     continue


c  If the column absorption optical depth exceeds 50, limit the 4-stream 
c  calculation to layers above that where this limit is reached.

         if(abstst.lt.50) numstop = layer
 
 30   continue
 

c the end of the layer loop ---------------------------------


c This section is related to the bottom boundary
 
      do 110 i = 1, 2
         smallb( i ) = abs( rmu( i) )
110   continue
 

c  compute delta(i)
c The sum of deltai is the part of the bottom boundary condition

      deltai( 1 ) =  suralb * a1 * smallb( 1 )
      deltai( 2 ) =  suralb * a2 * smallb( 2 )
 

c compute zeta ( only for the last layer )
c zeta is the part of the bottom boundary condition.
c The surface is a lumbertion surface.

      zeta = suralb / twopie + deltai( 1 ) * capz( 3,nlayers ) +
     &       deltai( 2 )*capz( 4,nlayers )
 

c fillabc and filld construct the 4-stream matrices from the notes of
c Tom Ackerman (notation of direction is opposit). 
c The matrix has a block tri-diagonal form and a matrix solver has been 
c chosen to exploit this structure.  
c The matrix equation (abc)x = d is solved where abc is a single matrix 
c with a being the lower block diagonal, b the center and c the upper.
 
      call fillabc (eval,efun,tau,deltai,nlayers,a,b,c)


      call filld (ri,nlayers,tau,amu0,capz,zeta,d)

      il = 1
      iu = numstop


c If the absorption optical depth limit is exceeded in the first
c level of the model domain, the matrix cannot be solved using
c this technique and all fluxes are set to zero.

      if(iu.gt.1) then
         call nbtrip(a, b, c, d, il, iu )
      else
         write(*,*) 'top layer is optically thick'
      end if
 
c      leqerr = ier
c      if ( ier .eq. 129 ) write( 6,990 )
 
c
c need to put the solutions into the p and q array
c
      if (iu.gt.1) then
         do 80 layer = 1, iu
            
            p( 1,layer ) = d(1, layer)
            p( 2,layer ) = d(2, layer)
            
            q( 1,layer ) = d(3, layer)
            q( 2,layer ) = d(4, layer)

 80      continue
      else
         iu = 0
      end if

      if(iu.lt.nlayers) then
         do 85 layer = iu+1, nlayers

            p( 1,layer ) = 0.0
            p( 2,layer ) = 0.0
            
            q( 1,layer ) = 0.0
            q( 2,layer ) = 0.0

 85      continue
      end if


c The following block compute upward and downward fluxes at ntau points
c Upwnward irradiance at 2 quadrature points 1 and 2 at the top of 
c the mth lauer is given by
c
c I(m-1)(tau(m-1),rmu(i=1or2) = 
c Sum(j=1,2){0.5[(pj(m-1)-qj(m-1)]exp(-kj(m-1)tau(m-2))gj(-j)
c           +0.5[(pj(m-1)+qj(m-1)]exp(-kj(m-1)tau(m-1))gj(j)}
c           + Z(m-1)(i)exp(-tau(m-1)/mu0)

c Similary, downwnward irradiance at 2 quadrature points 1 and 2 at the top 
c of the mth lauer is given by
c
c I(m-1)(tau(m-1),rmu(1=-1or-2) = 
c Sum(j=1,2){0.5[(pj(m-1)-qj(m-1)]exp(-kj(m-1)tau(m-2))gj(j)
c           +0.5[(pj(m-1)+qj(m-1)]exp(-kj(m-1)tau(m-1))gj(-j)}
c           + Z(m-1)(-i)exp(-tau(m-1)/mu0)
c
c where kj = eval, gj = efun

      summ   = 0.0
      ilyr = 1
      ntaujs = nlayers + 1
 
      do 90 n = 1, ntaujs
 
          if ( n  .eq.1 ) then
             ilyr  =  1
             taun    =  0.0
           else
             ilyr  =  n - 1
             taun    =  tau( n-1 )
           endif
 
           if ( ilyr .eq. 1 ) then
              taumx = tau( ilyr )
           else
              taumx = tau( ilyr ) + tau( ilyr - 1 )
           endif

          if (n.le.(numstop+1)) then
             exptm0  =  exp(  - taun / amu0 )
          else
             exptm0 = 0.0
          end if

          do 120 j = 1, 2
 
             expktn(j)  =  exp( - eval( j,ilyr ) * taun )
             expktm(j)  =  exp( - eval( j,ilyr ) * (taumx - taun) )
 
             conp(j) = 0.5 * ( p(j,ilyr) + q(j,ilyr) )
             conm(j)= 0.5 * ( p(j,ilyr) - q(j,ilyr) )
 
 120     continue
 
         do 130 i = 1, 4
 
            if ( i .le. 2 ) then
               ii  =  i + 2
            else
               ii  =  i - 2
            endif
 
            sum  =  0.0
            do 140 j = 1, 2
               sum = sum +
     &               conp(j) * efun( j,i,ilyr ) * expktn(j) +
     &               conm(j) * efun(j,ii,ilyr ) * expktm(j)
140         continue
 
            ri_qdt(i)  =  sum  +  capz( i,ilyr ) * exptm0
 
130      continue

         amu1  =  a1 * smallb( 1 )
         amu2  =  a2 * smallb( 2 )

c Following two equations assumes the solar constant at TOA is 1.0
 
         fup(n)  =  twopie * ( amu1 * ri_qdt(1) + amu2 * ri_qdt(2) )
         fdn(n)  =  twopie * ( amu1 * ri_qdt(3) + amu2 * ri_qdt(4) )
     &              +  exptm0

c Since all terms are proportional to mu0*pi*F0, we multiplied by mu0
 
         fup(n)  =  amu0 * fup(n)
         fdn(n)  =  amu0 * fdn(n)

c to obtain the direct beam, the delta assumption has to be removed

         if (n .eq. 1) then

           direct(1) = amu0

         else
 
           direct(n) = amu0 * exp(-taund(n-1)/amu0)

         endif
 
90    continue
 
 
      return
 
875   format(' ','  the normalized moments for layer',i5,'  are',
     &             '  omg0 = ',f13.5,'  momts = ',4f13.5)
910   format( 1h0 )
 
920   format(' ','  the delta function is in ')
922   format(' ','  the delta function is not in ')
 
990   format('0',' * * * * * w a r n i n g  * * * * * ',/
     &          10x,' the calculations are incorrect!!!')
 
      end
 
 
