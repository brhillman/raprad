      subroutine fillabc(eval,efun,tau,deltai,nlayers,a,b,c)
 
      implicit double precision( a-h, o-z )
c                                                                     
c define the dimensions used by the radiation model but might be      
c specified by an external model                                      
c                                                                     
c ivert  = maximum number of layers;                                  


c Inputs
c eval(1or2,j) Eigenvalues given by equation (31) in Liou (1974)
c              1 correspondes minus, 2 correapondes plus (j=1,nalyers)
c tau(j)       Total tau from TOA to bottom of jth layer with the delta
c              approximation
c efun(k,i,j)  Eigenfunctions given by equation 26 in Liou (1974)
c              k=1 for k1, k=2 for k2, i varyes 1 to 4 as rmu, (j=1,nlayers)
c deltai       Part of the bottom boundary condition
c nlayers      Total number of layers in the atmosphere

c Outputs
c a(4,4,ivert) a+b+c matrix contains boundary conditions
c b(4,4,ivert) a+b+c matrix contains boundary conditions
c c(4,4,ivert) a+b+c matrix contains boundary conditions

                                                                     
      parameter (ivert=200)                                            

      dimension efun(2,4,ivert), eval(2,ivert)
      dimension tau(0:ivert), deltai(2)
      dimension a(4,4,ivert), b(4,4,ivert), c(4,4,ivert)
      dimension sqgl1(2,2,ivert),sqgl2(2,2,ivert)
      dimension e1(2),e2(2),eta1(2),eta2(2)

	
c index function for moments and legendre polynomials

      lx( l )  =  l + 1
 

c this loop fills all the rows in the
c 'abc' matrix. there is only incindent flux affecting
c the top layer. so the part of the matrix that accounts
c for the layer above it is just zero for the first four rows.
c same idea for  the last layer.
c
c symmetry is relied on heavily.
c
c the eta is added and subtracted separately after the
c layer loop is finshed.
c
 
c functions for the matrix; beta, squigl, and eta.
c Notation is followed by Dr. Ackerman's note, that is downward irradiance
c is positive mu, upward is negative mu
c i corresponds irradiance direction in two quadrature direction
c i of rmu(i); rmu(1)=-rmu(3), rmu(2)=-rmu(4)
c
c The following block fills matrix (2*nlayers * 2*nlayers) of
c  layer # (boundary #)
c
c        1   bc0
c        2   abc0
c        2   0abc0
c        3    0abc0
c        3     0abc0
c                   .....
c        n               0abc0
c        n                0abc
c  bottom boundary          ab'         b' is slightly different from b


c The each a, b, and c matrix is 4*4
c       (k,j)
c
c       (1,1) (1,2) (1,3) (1,4)             (1,1) (1,2)   0     0
c   a = (2,1) (2,2) (2,3) (2,4)         b = (2,1) (2,2)   0     0 
c       (3,1) (3,2) (3,3) (3,4)               0     0   (3,3) (3,4)
c       (4,1) (4,2) (4,3) (4,4)               0     0   (4,3) (4,4)
c
c       (1,1) (1,2) (1,3) (1,4)             
c   c = (2,1) (2,2) (2,3) (2,4)          
c       (3,1) (3,2) (3,3) (3,4)               
c       (4,1) (4,2) (4,3) (4,4)               





      do 10 layer = 1, nlayers
 
         i = 0
 
         do 20 k = 1, 2
 
            kplus = k + 2
            i = i +1
            j = 1
 
            do 30 l = 1, 2

c  define components of the functions

               e1(j) = exp(-eval(j,layer)*tau(layer-1))
               e2(j) = exp(-eval(j,layer)*tau(layer))

c w1 and w2 correspond W of equation 26 in Liou (1974).

               w1 = efun(j,i+2,layer)
               w2 = efun(j,i,layer)
 
               beta1 = w1*e1(j)
               beta2 = w2*e2(j)
               sqgl1(j,i,layer) = w1*e2(j)
               sqgl2(j,i,layer) = w2*e1(j)
c	
               lplus = l+2 
               b( k,l,layer ) = beta1 + beta2
               b( kplus,lplus,layer ) = beta1-beta2

               if ( layer .ne. 1 ) then
 
                a(k,l,layer) = -.5*(sqgl1(j,i,layer-1)
     &                              +sqgl2(j,i,layer-1))
                a( kplus,l,layer ) = a( k,l,layer )
 
                a(k,lplus,layer) = -.5*(sqgl1(j,i,layer-1)-
     &                                   sqgl2(j,i,layer-1))
                a( kplus,lplus,layer ) = a( k,lplus,layer )

                c(k,l,layer-1) = -.5*(sqgl1(j,i,layer)
     &                              +sqgl2(j,i,layer))
                c( kplus,l,layer-1 ) = - c( k,l,layer-1 )
 
                c(k,lplus,layer-1) = .5*(sqgl1(j,i,layer)
     &                                   -sqgl2(j,i,layer)) 
                c( kplus,lplus,layer-1 ) = - c( k,lplus,layer-1 ) 
 
               endif

               j = j + 1
 
 30        continue
 
 20      continue
 
 10   continue
 
c now the eta's are added and subtracted from their
c spots in the matrix.
c bottom boundary condition

      do 80 k = 1, 2

       kplus = k+2
       j = 1
 
       do 90 l = 1, 2

         eta1(j) = .5*( ( deltai( 1 ) * efun( j,3,nlayers ) +
     &               deltai( 2 ) * efun( j,4,nlayers ) )*e2(j))

         eta2(j) = .5*( ( deltai( 1 ) * efun( j,1,nlayers ) +
     &               deltai( 2 ) * efun( j,2,nlayers ) )*e1(j))

        lplus = l+2
        b(k,l,nlayers) = b(k,l,nlayers) - (eta1(j)+eta2(j))
        b(k,lplus,nlayers) = b(k,lplus,nlayers) - (eta1(j)-eta2(j))
 
        b(kplus,l,nlayers) = b(kplus,l,nlayers ) + (eta1(j)+eta2(j))
        b(kplus,lplus,nlayers)=b(kplus,lplus,nlayers)+(eta1(j)-eta2(j))

        j = j + 1
 
90     continue
80    continue
 
 200  format(1p,12(f8.4,1x))
 
      return
      end
 
 
