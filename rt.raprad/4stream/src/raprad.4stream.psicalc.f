      subroutine psicalc( x, coomg0, coomg1, coomg2, psi )
 

c calculate psi functions for summations that is given equation 29
c in Liou (1974)

 
      implicit double precision( a-h, o-z )

      dimension psi(4)

c index function for moments and legendre polynomials
c
      lx( l )  =  l + 1
c
 
      xsq  =  x * x
      xcu  =  x * xsq
 
      psi( lx(0) )  =  1.0
      psi( lx(1) )  =  - coomg0 / x
      psi( lx(2) )  =    coomg1 * coomg0 / ( 2.0 * xsq )  -  0.5
      psi( lx(3) )  =  - coomg2 * coomg1 * coomg0 / ( 6.0 * xcu )  +
     &                          ( coomg2 + 4.0 * coomg0 ) / ( 6.0 * x )
 
      return
      end
 
