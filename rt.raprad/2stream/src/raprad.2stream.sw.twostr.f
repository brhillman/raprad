      subroutine twostr 
     & (nlayer,taul,w0,g0,rsfx,b1,b2,el1,el2,em1,em2,af,bf,ef,ak)

c
c    ******************************************************************
c    *  purpose             :  defines matrix properties and sets up  *
c    *                         matrix coefficients that do not depend *
c    *                         on zenith angle or temperature.        *
c    *  subroutines called  :  none                                   *
c    *  input               :  w0, g0                                 *
c    *  output              :  b1, b2, el1, el2, em1, em2, af, bf, ef *
c    * ****************************************************************

c define the dimensions used by the radiation model but might be
c specified by an external model
c
c ivert  = maximum number of layers;
c ilayer = maximum number of layer boundaries
c idbl   = twice the maximum number of layer boundaries

      parameter (ivert=200)
      parameter (irad=20)
      parameter (ilayer=ivert+1, idbl=2*ilayer)


c  double precision
       implicit real*8 (a-h, o-z)

       real taul(*), w0(*), g0(*)
       real rsfx
       dimension gami(ilayer), ak(ilayer)
       dimension b1(ilayer), b2(ilayer)
       dimension ee1(ilayer)
       dimension el1(ilayer), el2(ilayer), em1(ilayer), em2(ilayer)
       dimension af(idbl), bf(idbl), ef(idbl)


c       if(l.le.nsolp) then
c           u1i = sq3
c       else
c           u1i = 2
c       endif

c for shortwave
       u1i = sqrt(3.0)
       pi = 4.0*atan(1.0)
       tpi = 2.0 * pi
       jdble = 2 * nlayer
       jn = jdble - 1


       u1s  =  tpi/u1i


c      here we define layer properties following general scheme
c      of meador and weavor. then we set up layer properties
c      needed for matrix.

       do 14 j           = 1,nlayer

c      these are for twostream and hemispheric means
c These b1 and b2 are gamma1 and 2 defined in Table 1 in Toon at al. (1989)

         b1(j)    =  0.5*u1i*(2. - w0(j)*(1. + g0(j)))
         b2(j)    =  0.5*u1i*w0(j)*(1. - g0(j))

c equation 21 of Toon et al. (1989)

         ak(j)    =  sqrt(abs(b1(j)**2 - b2(j)**2))

c equation 22 of Toon et al. (1989)

         gami(j)  =  b2(j)/(b1(j) + ak(j))
         x1         =  ak(j)*taul(j)

         if(x1.gt.1000.)    x1=1000.

         ee1(j)   =  exp(-x1)

         if( x1.gt. 1000.) ee1(j)= 0.

c equation 44 of Toon et al., (1989)

         el1(j)   =  1.0 + gami(j) *ee1(j)
         em1(j)   =  1.0 - gami(j) * ee1(j)
         el2(j)   =  gami(j) + ee1(j)
         em2(j)   =  gami(j) - ee1(j)

   14  continue
c
c     we seek to solve ax(l-1)+bx(l)+ex(l+1) = d.
c     l=2n for even l, l=n+1 for odd l. the mean intensity (tmi/4pi)
c     and the net flux (fnet) are related to x's as noted in add.
c     first we set up the coefficients that are independent of solar
c     angle or temparature: a(i),b(i),e(i). d(i) is defined in add.
c
      j                   =  0
      do 18 jd               =  2,jn,2
        j               =  j + 1

c     here are the even matrix elements eq. 42 of Toon et al.

        af(jd)   = em1(j+1)*el1(j)-em2(j+1)*el2(j)
        bf(jd)   = em1(j+1)* em1(j)-em2(j+1)*em2(j)
        ef(jd)   = el1(j+1)*em2(j+1) - el2(j+1)*em1(j+1)

c     here are the odd matrix elements except for the top.
c     eq. 41 in Toon et al.

        af(jd+1) =  em1(j)*el2(j)-el1(j)*em2(j)
        bf(jd+1) =  el1(j+1)*el1(j) - el2(j+1)*el2(j)
        ef(jd+1) =  el2(j)*em2(j+1)-el1(j)*em1(j+1)

   18  continue
c
c     here are the top and bottom boundary conditions as well as the
c     beginning of the tridiagonal solution definitions. I assume
c     no diffuse radiation is incident at upper boundary.
c

      jdble = 2 * nlayer

      af(1)     = 0.0
      bf(1) = el1(1)
      ef(1) = -em1(1)
      af(jdble) = el1(nlayer)-rsfx*el2(nlayer)
      bf(jdble) = em1(nlayer)-rsfx*em2(nlayer)
      ef(jdble) = 0.0


      return
      end
