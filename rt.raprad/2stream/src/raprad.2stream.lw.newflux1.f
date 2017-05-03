
*******************************************************

      SUBROUTINE NEWFLUX1 (nlayer, taul, ptempg, ptemp, slope, 
     &           b3, gami, ak, gangle, gweight, y3, ee1,
     &           ck1, ck2, rsfx,u1i,u1s,
     &           direc, directu)

C Last Revision:  25 JAN 1995

C
C     **************************************************************
C     *  Purpose             :  Calculate upward and downward      *
C     *                         intensities and fluxes using Gauss *
C     *                         Quadrature angles and weights.     *
C     *  Subroutines Called  :  None                               *
C     *  Input               :  PTEMP, SLOPE, Y3, B3, EE1, EE2     *
C     *  Output              *:  DINTENT, UINTENT, DIREC, DIRECTU   *
C     * ************************************************************
C
c
      parameter (ivert=200)
      parameter (ilayer=ivert+1)
      parameter (ngauss = 3)   

      implicit double precision(a-h, o-z)

      real taul(ilayer), emis, rsfx
      dimension y3(ngauss,ilayer), ptemp(ilayer), slope(ilayer)
      dimension ak(ilayer), b3(ilayer), ee1(ilayer)
      dimension ck1(ilayer), ck2(ilayer), gami(ilayer)


      dimension a1(ilayer), a2(ilayer), a3(ilayer), a4(ilayer)
      dimension a7(ilayer)
      dimension y1(ngauss,ilayer), y2(ngauss,ilayer)
      dimension y4(ngauss,ilayer), y5(ilayer), y8(ngauss,ilayer)
      dimension gangle(ngauss), gweight(ngauss)
      dimension direc(ilayer), directu(ilayer), uintent(ngauss,ilayer)
      dimension dintent(ngauss,ilayer)

      tpi = 2.0 * 4.0 * atan(1.0)
      emis = 1.0 - rsfx 

c irs = 0 for no scattering, irs = 1 for including scattering
      irs = 1

      do 200 j           =  1,nlayer
       if(j.eq.1) then
        kindex = 1
       else
        kindex = j-1
       endif
c          do 100  l      =  nsolp+1,ntotal
c            here we do no scattering coefficients
c This part is forming Eq. 27 in Toon et al.
             a3(j)     =  ptemp(kindex)*tpi
             a4(j)     =  tpi*slope(j)
             a7(j)     =  a3(j)
             y5(j)     =  a4(j)*taul(j)
c 100      continue
c
c         here we do scattering
c For this code we always compute scattering.
c This part forms Eq 27.
          if(irs .ne. 0) then
c              do 50 l    =  nsolp+1,ntotal
                x4       =  slope(j)*(tpi*b3(j)-u1s)
                a1(j)  =  u1i - ak(j)
                a2(j)  =  gami(j)*(ak(j)+u1i)
                a3(j)  =  a3(j)+x4
                a7(j)  =  a7(j)-x4
c 50           continue
          endif
  200 continue
c
c     calculations for all gauss points. here we do no scattering
c
      do 400       j         =  1,nlayer
         do 350    i         =  1,ngauss
c            do 300 l         =  nsolp+1,ntotal
               y1(i,j)  =  0.0
               y2(i,j)  =  0.0
               y4(i,j)  =  a7(j) - a4(j)*gangle(i)
               y8(i,j)  =  a3(j)+a4(j)*gangle(i)
c 300        continue
c
c           here we do scattering
            if(irs .ne. 0) then
c              do 325 l    =  nsolp+1,ntotal
               ya     =  a1(j)*(y3(i,j)-ee1(j))/
     &                   (ak(j)*gangle(i)-1.)
               yb     =  a2(j)*(1.- ee1(j)*y3(i,j))/
     &                   (ak(j)*gangle(i)+1.)
               ckp    =  ck1(j)+ck2(j)
               ckm    =  ck1(j) -ck2(j)
               y1(i,j)  =  ckp*yb+ckm*ya
               y2(i,j)  =  ckp*ya+ ckm*yb
c  325          continue
            endif
 350     continue
 400  continue
c
      do 450 j             =  1,nlayer
c         do 425  l         =  nsolp+1,ntotal
            direc(j)     =  0.0
            directu(j)   =  0.0
c 425     continue
 450  continue
c
c     direc is downward flux. directu is upward flux.
c     calculate dintent the downward intensity and direc the downward fl
c
          do 500 i             = 1,ngauss
c             do 475 l          = nsolp+1,ntotal
                dintent(i,1) = (1.-y3(i,1))*y4(i,1) + y1(i,1)
                direc(1)     = direc(1) + dintent(i,1) * gweight(i)
c 475         continue
 500      continue
c
c      dintent is downward intensity * tpi. direc is the downward flux.
c
       do 530        j           = 2,nlayer
           do 520    i           = 1,ngauss
c              do 510 l           = nsolp+1,ntotal
                 dintent(i,j)  = dintent(i,j-1)*y3(i,j)
     &                             +y1(i,j)+y5(j)+
     &                             (1.-y3(i,j))*y4(i,j)

                 direc(j)      = direc(j) + dintent(i,j) * gweight(i)
c 510          continue
 520       continue
 530   continue
c
c     uintent is the upward intensity * tpi. directu is the upward flux.
c     assume that the reflectivity is lambert.
c
       do 570     i               =  1,ngauss
c          do 560  l               =  nsolp+1,ntotal
             uintent(i,nlayer)  =  ptempg*emis
     &                               *tpi+rsfx*direc(nlayer)*2.0
             directu(nlayer)    =  directu(nlayer) 
     &                          +  uintent(i,nlayer)*gweight(i)
c 560      continue
 570   continue
c
      do 650        m              = 2,nlayer
          j                        = nlayer-m+1
          do 640    i              = 1,ngauss
c             do 630 l              = nsolp+1,ntotal
                 uintent(i,j)    = (uintent(i,j+1)-y5(j+1))
     1                               *y3(i,j+1)+y2(i,j+1)+
     2                               (1.-y3(i,j+1))*y8(i,j+1)
                 directu(j)     = directu(j) 
     &                          + uintent(i,j) * gweight(i)
c 630         continue
 640      continue
 650  continue
c
      return
      end
