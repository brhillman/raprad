      subroutine filld (ri,nlayers,tau,amu0,capz,zeta,d)
 
      implicit double precision( a-h, o-z )

c                                                                     
c define the dimensions used by the radiation model but might be      
c specified by an external model                                      
c                                                                     
c ivert  = maximum number of layers;                                  


c Inputs
c ri           Downward diffuse irradiance at the top of the atmosphere
c nlayers      Total number of layers in the atmosphere
c tau(j)       Total tau from TOA to bottom of jth layer with the delta
c              approximation
c amu0         Cosine of the solar zenith angle (Radian)
c capz(i,j)    Z function defined by Liou (1974) 27

c Outputs
c d(4,ivert)   Matrix that contains Z function
                                                                     
      parameter (ivert=200)                                            

      dimension tau(0:ivert), capz(4,ivert), d(4,ivert)
      real amu0


c index function for moments and legendre polynomials

      lx( l )  =  l + 1
 
      ri = 0.0

      do 10 layer = 1, nlayers
 
         expon1 = exp( -tau( layer-1 ) / amu0 )
         expon2 = exp( -tau( layer )   / amu0 )
 
         do 20 j = 1, 2
 
            if ( layer .eq.1 ) then
 
               part1 = ( ri - capz( j+2,layer ) )
               part2 = ( capz( j,layer ) - capz( j,layer+1 ) )
 
            else if ( layer .eq. nlayers ) then
 
               part1 = ( capz( j+2,layer-1 ) - capz( j+2,layer ) )
               part2 = ( capz( j,layer )  - zeta )
 
            else
 
               part1 = ( capz( j+2,layer-1 ) - capz( j+2,layer ) )
               part2 = ( capz( j,layer )     - capz( j,layer+1 ) )
 
            endif
 
            d( j,layer )     = part1 * expon1  -  part2 * expon2
            d( j + 2,layer ) = part1 * expon1  +  part2 * expon2
 
 20      continue
 
 10   continue

 200  format(e11.4,2x,e11.4)
 
 
      return
      end
 
 
