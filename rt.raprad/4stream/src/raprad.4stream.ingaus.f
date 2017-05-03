      subroutine ingaus(rmu1,rmu2,pleg,rmu)
 

c initialize gauss quad values and legendre polynomials
c Equation 28 of Liou (1974)

      implicit double precision( a-h, o-z )

      dimension rmu(4), pleg(4,4)

                            
c index function for moments and legendre polynomials

      lx( l )  =  l + 1
  
      rmu(1)  =  rmu1
      rmu(2)  =  rmu2

      do 10 i = 3, 4
         rmu(i)  =  - rmu( i-2 )
 10   continue
   
c set up legendre polynomial array

      do 20 i = 1, 4
 
         rmui    =  rmu(i)
         rmuisq  =  rmui * rmui
         rmuicu  =  rmui * rmuisq
 
         pleg( lx(0),i )  =  1.0
         pleg( lx(1),i )  =  rmu(i)
         pleg( lx(2),i )  =  0.5 * ( 3.0 * rmuisq - 1.0 )
         pleg( lx(3),i )  =  0.5 * ( 5.0 * rmuicu - 3.0 * rmui )
 
 20   continue
 
      return
      end
