      subroutine plnk(e,t1,d)
c
c     ******************************************************
c     *  purpose             :  calculate planck function  *
c     *  subroutines called  :  none                       *
c     *  input               :  wave, ncount               *
c     *  output              :  plank                      *
c     * ****************************************************
c
c  this subroutine computes the integral of the planck function between
c  zero and the specified value of lambda.  thus (using xl as lambda)
c  we want to integrate
c  r = integral(xl=0 to xl=xlspec) ( c1*xl**-5* / (exp(c2/xl*t)-1) )*dxl
c  substituting u=c2/(xl*t), the integral becomes
c  r = a constant times integral (uspec to infinity) of
c            ( u**3 / (exp(u) - 1) )*du
c  the approximations shown here are on page 998 of abramowitz and segun
c  under the heading of debye functions.  c2 is the product of planck's
c  constant and the speed of light divided by boltzmann's constant.
c  c2 = 14390 when lambda is in microns.
c  the factor 0.15399 is the reciprocal of six times
c  the sum of (1/n**2) for all n from one to infinity.  it is chosen to
c  normalize the integral to a maximum value of unity.
c  radiation in real units is obtained by multiplying the integral by
c  the stefan-boltzmann constant times t**4.
c
c
c  double precision
       implicit real*8 (a-h, o-z)
c
      dimension am(5)
c
      d            =   0.0
      v1           =   e/t1
c
      if (v1 .le. 1.) then
         d         =  1.0 - 0.15399*v1**3 *
     1                (1./3.-v1/8. + v1**2/60. - v1**4/5040. +
     2                v1**6/272160. - v1**8/13305600         )
      endif
c
      if ( v1 .gt. 1. .and. v1 .le. 50.) then
         do 100 m   =  1,5
            a       =  float(m)*v1
            am(m)   =  0.15399 * exp(-a)/m**4 *
     1                 (((a+3.)*a+6.)*a+6.)
 100     continue
c
         d          =  am(1)+am(2)+am(3)+am(4)+am(5)
      endif
c
      d             =  d*t1**4
c
      return
      end
