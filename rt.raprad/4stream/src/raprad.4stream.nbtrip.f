c...
c...subroutine to solve non-periodic block tridiagonal 
c...system of equations without pivoting strategy
c...with the dimensions of the block matrices being
c...n x n (n is any number greater than 1).
c...
        subroutine nbtrip(a,b,c,d,il,iu)
c...
        implicit double precision (a-h, o-z)
        dimension a(1),b(1),c(1),d(1)
c        integer order,ordsq
c...
c...a = sub diagonal matrix
c...b =     diagonal matrix
c...c = sup diagonal matrix
c...d = right hand side vector
c...il = lower value of index for which matrices are defined
c...iu = upper value of index for which matrices are defined
c...     (solution is sought for btri(a,b,c)*x = d
c...     for indices of x between il and iu (inclusive).
c...     solution written in d vector (original contents
c...     are overwritten)
c...order = order of a,b,c matrices and length of d vector
c...     at each point denoted by index i
c...     (order can be any integer greater than 1).
c...
c...the matrices and vectors are stored in single subscript form
c...
c        order = 4
c        ordsq = 16
c...
c...forward elimination
c...
        i = il
        i0mat = 1+(i-1)*16
        i0vec = 1+(i-1)*4
c

	call ludeco(b(i0mat))

	call lusolv(b(i0mat),d(i0vec),d(i0vec))

   	i0matj = i0mat
        call lusolv(b(i0mat),c(i0matj),c(i0matj))
   	i0matj = i0mat+4
        call lusolv(b(i0mat),c(i0matj),c(i0matj))      
   	i0matj = i0mat+8
        call lusolv(b(i0mat),c(i0matj),c(i0matj))
   	i0matj = i0mat+12
        call lusolv(b(i0mat),c(i0matj),c(i0matj))

  200	continue
	i = i+1
	i0mat = 1+(i-1)*16
	i0vec = 1+(i-1)*4
	i1mat = i0mat-16
	i1vec = i0vec-4
	call smulput(a(i0mat),d(i1vec),d(i0vec))

        i0matj=i0mat
        i1matj=i1mat
        call smulput(a(i0mat),c(i1matj),b(i0matj))
        i0matj=i0mat+4
        i1matj=i1mat+4
        call smulput(a(i0mat),c(i1matj),b(i0matj))
        i0matj=i0mat+8
        i1matj=i1mat+8
        call smulput(a(i0mat),c(i1matj),b(i0matj))
        i0matj=i0mat+12
        i1matj=i1mat+12
        call smulput(a(i0mat),c(i1matj),b(i0matj))

	call ludeco(b(i0mat))
	call lusolv(b(i0mat),d(i0vec),d(i0vec))
	if(i.eq.iu) go to 500

	i0matj = i0mat
	call lusolv(b(i0mat),c(i0matj),c(i0matj))
	i0matj = i0mat+4
	call lusolv(b(i0mat),c(i0matj),c(i0matj))		
	i0matj = i0mat+8
	call lusolv(b(i0mat),c(i0matj),c(i0matj))
	i0matj = i0mat+12
	call lusolv(b(i0mat),c(i0matj),c(i0matj))
c
	goto 200
  500	continue
c...
c...back substitution
c...

	i = iu
  600	continue
	i = i-1
	i0mat = 1+(i-1)*16
	i0vec = 1+(i-1)*4
	i1vec = i0vec+4
	call mulput(c(i0mat),d(i1vec),d(i0vec))
	if (i.gt.il) go to 600
c...
	return 
	end
c...
c...subroutine to calculate l-u decomposition
c...of a given matrix a and store result in a
c...(no pivoting stategy is employed)
c...
	subroutine ludeco(a)
        implicit double precision (a-h, o-z)
c...
	dimension a(4,1)
c	integer order
c...
	do 8 jc=2,4
    8  	a(1,jc) = a(1,jc)/a(1,1)
	jrjc = 1
   10	continue
	jrjc = jrjc+1
	jrjcm1 = jrjc-1
	jrjcp1 = jrjc+1
	do 14 jr=jrjc,4
	  sum = a(jr,jrjc)
	  do 12 jm=1,jrjcm1
   12	  sum = sum-a(jr,jm)*a(jm,jrjc)
   14	a(jr,jrjc) = sum
	if (jrjc.eq.4) return
	do 18 jc=jrjcp1,4
	  sum = a(jrjc,jc)
	  do 16 jm=1,jrjcm1
   16	  sum = sum-a(jrjc,jm)*a(jm,jc)
   18	a(jrjc,jc) = sum/a(jrjc,jrjc)
	go to 10
	end
c...
c...subroutine to multiply a vector b by a matrix a,
c...subtract result from another vector c and store
c...result in c.  thus vector c is overwritten
c...
	subroutine mulput(a,b,c)
        implicit double precision (a-h, o-z)
c...
	dimension a(1),b(1),c(1)
c	integer order
c...
	do 200 jr=1,4
	  sum = 0.0
          sum = a(jr)*b(1)+a(jr+4)*b(2)+a(jr+8)*b(3)+a(jr+12)*b(4)	
c
  200	c(jr) = c(jr)-sum
c...
	return
	end
c...
c...subroutine to multiply a vector b by a matrix a,
c...subtract result from another vector c and store
c...result in c.  thus vector c is overwritten
c...same as mulput but exploits symmetry of 'a' matrix
c...
	subroutine smulput(a,b,c)
        implicit double precision (a-h, o-z)
c...
	dimension a(1),b(1),c(1)
c	integer order
c...
	do 200 jr=1,2
	  sum = 0.0
          sum = a(jr)*b(1)+a(jr+4)*b(2)+a(jr+8)*b(3)+a(jr+12)*b(4)	
c
  	c(jr) = c(jr)-sum
  200   c(jr+2) = c(jr+2)-sum
c...
	return
	end
c...
c...subroutine to solve linear algebraic system of
c...equations a*c=b and store results in vector c.
c...matrix a is input in l-u decomposition form.
c...(no pivoting strategy has been employed to
c...compute the l-u decomposition of the matrix a).
c...
	subroutine lusolv(a,b,c)
c...
        implicit double precision (a-h, o-z)
	dimension a(4,1),b(1),c(1)
c	integer order
c...
c...first l(inv)*b
c...
	c(1) = c(1)/a(1,1)
	do 14 jr=2,4
	  jrm1 = jr-1
	  sum = b(jr)
	  do 12 jm=1,jrm1
   12	  sum = sum-a(jr,jm)*c(jm)
   14	c(jr) = sum/a(jr,jr)
c...
c...next u(inv) of l(inv)*b
c...
	do 18 jrjr=2,4
	  jr = 5-jrjr
	  jrp1 = jr+1
	  sum = c(jr)
	  do 16 jmjm=jrp1,4
	    jm = 4-jmjm+jrp1
   16	  sum = sum-a(jr,jm)*c(jm)
   18	c(jr) = sum
c...
	return
	end
