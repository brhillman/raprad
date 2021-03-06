      subroutine add
     &  (nlayer,taul,w0,g0,rsfx,opd,opdnd,ak,b1,b2,b3,em1,em2,
     &   el1,el2,af,bf,ef,u0,fnet,sol_fluxup,sol_fluxdn,direct_nd)

c     ***************************************************************
c     *  purpose             :  defines source terms, form matrix   *
c     *                         for multiple layers and solve tri-  *
c     *                         diagnol equations to obtain mean    *
c     *                         intensity and net flux.             *
c     *  subroutines called  :  none                                *
c     *  input               :  nlayer,taul,w0,g0,rsfx,opd,ak,b1,b2
c     *                         b3,em1,em2,el1,el2,af,bf,ef         *
c     *  output              :  fnet,sol_fluxup,sol_fluxdn          *
c     * *************************************************************

c define the dimensions used by the radiation model but might be
c specified by an external model

c ivert  = maximum number of layers;
c ilayer = maximum number of layer boundaries
c idbl   = twice the maximum number of layer boundaries

      parameter (ivert=200)
      parameter (ilayer=ivert+1, idbl=2*ilayer)

c  double precision

       implicit real*8 (a-h, o-z)

       real taul(*), w0(*), g0(*)
       real opd(ilayer), opdnd(ilayer)
       real rsfx, u0
       dimension ak(ilayer)
       dimension b1(ilayer), b2(ilayer), b3(ilayer)
       dimension el3(ilayer), ee3(ilayer)
       dimension fnet(ilayer), sol_fluxup(ilayer), sol_fluxdn(ilayer)
       dimension direct(ilayer), diffuse(ilayer), tmi(ilayer)
       dimension direct_nd(ilayer)
       dimension cpb(ilayer), cp(ilayer), cmb(ilayer), cm(ilayer)
c       dimension slope(ilayer)
       dimension em1(ilayer), em2(ilayer)
       dimension el1(ilayer), el2(ilayer)
       dimension as(idbl),af(idbl), bf(idbl), df(idbl), ds(idbl)
       dimension ef(idbl), xk(idbl)
       dimension ck1(ilayer), ck2(ilayer)

       data epsilon / 1.0e-15 /


       jdble = nlayer * 2
       jn = jdble - 1
       sol = 1.0
       sq3 = 3.0**(1.0/2.0)

c  flag for solar(isl) and IR(irs)
       isl = 1
       irs = 0

c     this subroutine forms the matrix for the multiple layers and
c     uses a tridiagonal routine to find radiation in the entire
c     atmosphere.
c
c     ******************************
c     *   calculations for solar   *
c     ******************************
      if(isl .ne. 0)  then
        du0                =  1./u0

        do 10 j            =  1,nlayer

          if(j.ne.1) then

            j1=j-1

          else

            j1=j

          endif

          b3(j)     =  0.5*(1.-sq3*g0(j)*u0)
          b4          =  1. - b3(j)
          x2          =  taul(j)*du0

          if(x2.gt.1000.)  x2 = 1000.

          ee3(j)    =  exp(-x2)
          x3          =  opd(j)*du0

          if(x3.gt.1000.)  x3 = 1000.

          el3(j)    =  exp(-x3)*sol

          if(el3(j).ge.1000.)  el3(j)=0.0

          direct(j) = u0*el3(j)
          c1          =  b1(j) - du0
          c2          =  ak(j)*ak(j) - du0*du0

          if(abs(c2).le.epsilon)   c2=epsilon


c equation 23 in Toon et al. (1989)

          cp1         =  w0(j)*(b3(j)*c1+b4*b2(j))/c2
          cpb(j)    =  cp1 * el3(j)

          if(j.ne.1) then
            x4 = el3(j1)
          else
            x4 = sol
          endif

          cp(j)     =  cp1 * x4

c equation 24 in Toon et al. (1989)

          cm1         =  ( cp1*b2(j) + w0(j)*b4 )/c1
          cmb(j)    =  cm1 * el3(j)
          cm(j)     =  cm1 * x4

 10     continue

c       calculate sfcs, the source at the bottom.

        sfcs         =  direct(nlayer) * rsfx
c
       end if
c
c     ******************************
c     * calculations for infrared. *
c     ******************************
c
c      if(irs .ne. 0)  then
c
c        do 30 j           =   1,nlayer
c
c          if(j.eq.1) then
c            kindex = 1
c          else
c            kindex = j-1
c          endif
c
c          b3(j)     = 1.0/(b1(j)+b2(j))
c          cp(j)     = (ptemp(kindex)+slope(j)*b3(j))*u1s
c          cpb(j)    = cp(j) + slope(j)*taul(j)*u1s
c          cm(j)     = (ptemp(kindex)-slope(j)*b3(j))*u1s
c          cmb(j)    = cm(j) + slope(j)*taul(j)*u1s
c          el3(j)    = 0.0
c          direct(j) = 0.0
c          ee3(j)    = 0.0
c
c 30     continue
c
c        sfcs          = emis*ptempg*pi
c
c      end if

      j                =  0

      do 42 jd         =  2,jn,2
        j             =  j + 1

c           here are the even matrix elements
        df(jd) = (cp(j+1) - cpb(j))*em1(j+1) -
     $ (cm(j+1) - cmb(j))*em2(j+1)

c           here are the odd matrix elements except for the top.

        df(jd+1) =  el2(j) * (cp(j+1)-cpb(j)) +
     &                    el1(j) * (cmb(j) - cm(j+1))
 42   continue



c     here are the top and bottom boundary conditions as well as the
c     beginning of the tridiagonal solution definitions. i assume no
c     diffuse radiation is incident at the top.
c
      df(1)     = -cm(1)
      df(jdble) = sfcs+rsfx*cmb(nlayer)-cpb(nlayer)
      ds(jdble) = df(jdble)/bf(jdble)
      as(jdble) = af(jdble)/bf(jdble)
c
c     ********************************************
c     *     we solve the tridiagonal equations   *
c     ********************************************
c
c This block is following eq. 45, 46, and 47 in Toon et al. (1998)

      do 47 j           = 2, jdble
        x               = 1./(bf(jdble+1-j) -
     1                        ef(jdble+1-j)*as(jdble+2-j))
        as(jdble+1-j)   = af(jdble+1-j)*x
        ds(jdble+1-j)   = (df(jdble+1-j) - ef(jdble+1-j)
     2                        *ds(jdble+2-j))*x
  47  continue

      xk(1)    = ds(1)

      do 50 j       = 2, jdble
            xk(j) = ds(j) - as(j)*xk(j-1)
  50  continue

c  ***************************************************************
c     calculate layer coefficients, net flux and mean intensity
c  ***************************************************************
      
      do 60 j = 1, nlayer

        sol_fluxdn(j) = 0.0
        sol_fluxup(j) = 0.0
        direct_nd(j) = 0.0

 60   continue

      do 62 j= 1,nlayer

c Yl; l = odd => Yl = Y1n,   l = even => Yl = Y2n

        ck1(j)   = xk(2*j-1)
        ck2(j)   = xk(2*j)

c equation 48 of Toon et al. (1989)

        fnet(j)  = ck1(j)  *( el1(j) -el2(j))   +
     3                 ck2(j) *( em1(j)-em2(j) ) + cpb(j) -
     4                 cmb(j) - direct(j)

c diffuse component of solar radiation

        diffuse(j) = ck1(j)*el2(j) + ck2(j)*em2(j)
     +                                     + cmb(j)
c
        tmi(j)     =  el3(j) + u1i * ( ck1(j)  *
     5                   (el1(j) + el2(j))   +
     6                    ck2(j) * ( em1(j)+em2(j) ) +
     7                    cpb(j) + cmb(j) )

        sol_fluxup(j) = sol_fluxup(j) + ck1(j)*el1(j)
     &                + ck2(j)*em1(j) + cpb(j)
        sol_fluxdn(j) = sol_fluxdn(j) + ck1(j)*el2(j)
     &                + ck2(j)*em2(j)+cmb(j)+direct(j)

        direct_nd(j) = exp(-opdnd(j)*du0)

   62 continue


 400  format (255f8.1)
 401  format (/, ' fnet for ', i5, 'wavelengths. ')
 402  format (' layer  ', i5, ' of ', i5)
      return
      end
