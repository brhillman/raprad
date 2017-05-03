      subroutine gauleg(x1, x2, x, wt, n)

c references Numerical Recipes, W. H. Press, B. P. Flannery, S. A. Teukolsky,
c and W. T. Vettering, Cambridge University Press (1986) p125

c Given the lower and upper limits of integration x1 and x2, and given
c n, this routine returns arrays x and w of length n, containing
c the abscissas and weights of the Gauss-Legendre n-point quadrature
c formular.


c High precision is a good idea for this routine.
      implicit real*8 (a-h, o-z)

      real*8 x1, x2, x(n), wt(n)

c increase if you don't have this floating precision
      parameter (eps = 3.d-14)

c The roots are symmetric in the interval, so we only have to find 
c half of them
      m = (n + 1) / 2

      xm = 0.5d0 * (x2 + x1)
      xl = 0.5d0 * (x2 - x1)

c Loop over the desired roots
      do 12 i = 1, m
        z = cos(3.141592654d0 * (i - 0.25d0) / (n + 0.5d0))

c Starting with the above approximation to the ith root, we enter the 
c main loop of refinment by Newton's method.

    1   continue
    
        p1 = 1.d0
        p2 = 0.d0

c Loop up the recurrence to get the legendre polynomial evaluated at z.
        do 11 j = 1, n
          p3 = p2
          p2 = p1
          p1 = ((2.d0 *j-1.d0)*z*p2-(j - 1.d0)*p3)/j
   11   continue

c p1 is now the desired legendre polynomial. We next compute pp,
c its derivative, by a standard relation involving also p2,
c the polynomial of one lower order.
        pp = n * (z * p1 - p2) / (z * z - 1.d0)
        z1 = z
        z = z1 - p1 / pp

      if (abs(z - z1) .gt. eps) go to 1

c scale the root to the desired interval
      x(i) = xm - xl * z

c and put in its symmetric counterpart
      x(n+1-i) = xm + xl * z

c compute the weight
      wt(i) = 2.d0 * xl / ((1.d0 - z * z) * pp * pp)

c and its symmetric counterpart
      wt(n+1-i) = wt(i)

   12 continue

      return
      end
