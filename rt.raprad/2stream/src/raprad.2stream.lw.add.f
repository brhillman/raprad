      subroutine add
     &  (nlayer,taul,w0,g0,rsfx,opd,opdnd,ak,b1,b2,b3,em1,em2,
     &   el1,el2,af,bf,ef,u0,slope,ptempg,ptemp,u1i,u1s,
     &   fnet,sol_fluxup,sol_fluxdn,direct_nd,
     &   irflag,ck1, ck2)

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

       implicit real*8 (a-h, o-z)

       integer ivert,ilayer,idbl
       parameter (ivert=200)
       parameter (ilayer=ivert+1, idbl=2*ilayer)

c  double precision


       real taul(*), w0(*), g0(*)
       real opd(ilayer), opdnd(ilayer)
       real rsfx, u0, emis
       dimension ak(ilayer)
       dimension b1(ilayer), b2(ilayer), b3(ilayer)
       dimension el3(ilayer), ee3(ilayer)
       dimension fnet(ilayer), sol_fluxup(ilayer), sol_fluxdn(ilayer)
       dimension direct(ilayer), diffuse(ilayer), tmi(ilayer)
       dimension direct_nd(ilayer)
       dimension cpb(ilayer), cp(ilayer), cmb(ilayer), cm(ilayer)
       dimension slope(ilayer), ptemp(ilayer)
       dimension em1(ilayer), em2(ilayer)
       dimension el1(ilayer), el2(ilayer)
       dimension as(idbl),af(idbl), bf(idbl), df(idbl), ds(idbl)
       dimension ef(idbl), xk(idbl)
       dimension ck1(ilayer), ck2(ilayer)

       integer jdble,jn,nlayer,irflag,kindex
       real*8 pi
       jdble = nlayer * 2
       jn = jdble - 1
c      sol = 1.0
c      sq3 = 3.0**(1.0/2.0)
       pi = 4.0 * atan(1.0)

c     this subroutine forms the matrix for the multiple layers and
c     uses a tridiagonal routine to find radiation in the entire
c     atmosphere.
c




c
c     ******************************
c     * calculations for infrared. *
c     ******************************
c
      if(irflag .eq. 1)  then
        emis = 1.0 - rsfx

        do 30 j           =   1,nlayer
          if(j.eq.1) then
            kindex = 1
          else
            kindex = j-1
          endif

          b3(j)     = 1.0/(b1(j)+b2(j))
          cp(j)     = (ptemp(kindex)+slope(j)*b3(j))*u1s
          cpb(j)    = cp(j) + slope(j)*taul(j)*u1s
          cm(j)     = (ptemp(kindex)-slope(j)*b3(j))*u1s
          cmb(j)    = cm(j) + slope(j)*taul(j)*u1s
          el3(j)    = 0.0
          direct(j) = 0.0
          ee3(j)    = 0.0
 30     continue

        sfcs          = emis*ptempg*pi
      end if

      j = 0

      do 42 jd =  2,jn,2
        j =  j + 1
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

        tmi(j)     =  el3(j) + u1i * ( ck1(j)  *
     5                   (el1(j) + el2(j))   +
     6                    ck2(j) * ( em1(j)+em2(j) ) +
     7                    cpb(j) + cmb(j) )

        sol_fluxup(j) = sol_fluxup(j) + ck1(j)*el1(j)
     &                + ck2(j)*em1(j) + cpb(j)
        sol_fluxdn(j) = sol_fluxdn(j) + ck1(j)*el2(j)
     &                + ck2(j)*em2(j)+cmb(j)+direct(j)

!        direct_nd(j) = exp(-opdnd(j)*du0)

   62 continue



      return
      end
