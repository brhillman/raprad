      subroutine calceig(omega,rmu,pleg,a1,a2,eval,efun,layer)

c Input 
c omega,rmu,pleg,a1,a2

c Output
c efun
 

c compute eiginvalues and eigenfunctions

      implicit double precision( a-h, o-z )

c ivert  = maximum number of layers;                                  
                                                                 
      parameter (ivert=200)

                                           
      dimension omega(5), rmu(4), pleg(4,4)
      dimension efun(2,4,ivert), eval(2,ivert)
      dimension tt(2), tp(2), psi(4)
 

c index function for moments and legendre polynomials

      lx( l )  =  l + 1
 
      omg00 =  omega( lx(0) )
      omg1  =  omega( lx(1) )
      omg2  =  omega( lx(2) )
      omg3  =  omega( lx(3) )
 
      coomg0  =  1.0 - omega( lx(0) )
      coomg1  =  3.0 - omega( lx(1) )
      coomg2  =  5.0 - omega( lx(2) )
 
      term  =  coomg2  +  4.0 * coomg0
 

c compute eigenvalues

      do 10 i = 1, 2
 
         rmui   =  rmu(i)

c The follwing two equations are (34) and (35) for Liou (1974)
c tt is t and tp is t prime
c i corresponds the subscripts 1 and 2

         tt(i)  =  omg00  +  omg1 * coomg0 * rmui * pleg( lx(1),i ) -
     &             0.5 * omg2 * pleg( lx(2),i )  -
     &             omg3 * term * rmui * pleg( lx(3),i ) / 6.0
 
         tp(i)  =  0.5 * coomg0 * coomg1 * ( omg2 * pleg( lx(2),i )  +
     &             omg3 * coomg2 * rmui * pleg( lx(3),i ) / 3.0 )
10    continue

      rmu1sq  =  rmu(1) * rmu(1)
      rmu2sq  =  rmu(2) * rmu(2)

      term  =  rmu1sq * rmu2sq

c The following three equations are equation (32) and (33) of Liou 1974
c cc correspondes c in these equations, and bb correspondes b. 

      cc  =  ( 1.0  -  a1 * tt(1)  -  a2 * tt(2) )  /  term
      cc  =  cc  +  a1 * tp(1) / rmu1sq  +  a2 * tp(2) / rmu2sq
      bb  =  ( a1*tt(1) - 1.0 ) / rmu1sq  +  ( a2*tt(2) - 1.0 ) / rmu2sq
 
      term  =  sqrt( bb * bb  -  4.0 * cc  )

c These eigenvalues are given by equation 31 in Liou (1974)
 
      eval( 1,layer )  =  sqrt( - 0.5 * bb  -  0.5 * term )
      eval( 2,layer )  =  sqrt( - 0.5 * bb  +  0.5 * term )
 

c compute eigenfunctions W given by equation 26 in Liou (1974)
c efun correspondes to W in the paper.
c j=1 for k1, j=2 for k2, X varies as rmu(i)

      do 20 j = 1, 2
 
        ekj  =  eval( j,layer )

c psi functions are given by equation 29 in Liou (1974)

        call psicalc( ekj, coomg0, coomg1, coomg2, psi )
 
        do 30 i = 1, 4
 
          rmui  =  rmu(i)
          term  =  1.0  /  ( 1.0  +  rmui * ekj )
          sum   =  0.0
 
          do 50 l = 0, 3
            sum  =  sum  +  omega( lx(l) ) * psi( lx(l) ) *
     &                 pleg( lx(l),i )
 50       continue
 
          efun( j,i,layer )  =  term * sum
 
 30     continue
 
 20   continue
 
 
      return
      end
 
 
