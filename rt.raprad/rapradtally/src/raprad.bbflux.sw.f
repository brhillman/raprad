      subroutine bb_flux(alt, press, wave, nlayer, n_band, nsub_band,
     &                   u0, thz, gmt_time, sol_const)

c This program computes broad band flux based on the output.
c Currently, no IR is treated.

      parameter(maxsub_band = 1000)
      parameter(maxlayers = 200)
      parameter(maxband = 200)

      real sol_fluxup_sub(maxsub_band,maxlayers)
      real sol_fluxdn_sub(maxsub_band,maxlayers)
      real sol_direct_sub(maxsub_band,maxlayers)
      real sol_fluxup(maxband,maxlayers)
      real sol_fluxdn(maxband,maxlayers)
      real sol_direct(maxband,maxlayers)
      real total_sol_fluxup(maxlayers), total_sol_fluxdn(maxlayers)
      real total_sol_direct(maxlayers), total_sol_diff(maxlayers)
      real sol_absop(maxband)

      real wave(*)

      real*8 sol_const(*)
      real*8 alt(*), press(*) 

      integer nlayer, n_band, nsub_band(*)

      character*100 in_file

      jlayers = nlayer + 1

      ntot_band = 0

      do 1001 i = 1, n_band
        ntot_band = ntot_band + nsub_band(i)
 1001 continue

      in_file(1:14)
     &  = '../../results/'

      in_file(15:100) = 'raprad.sw.out'

      open(11,file=in_file,status='old')
      
      do 1000 n = 1, ntot_band
        do 1100 j = 1, jlayers
          read(11,*) i, sol_fluxup_sub(n,j), sol_fluxdn_sub(n,j),
     &      sol_direct_sub(n,j)
 1100   continue
 1000 continue

      close(unit=11)

c computing band flux at each level

      do 2000 j = 1, jlayers
        do 2100 i = 1, nband
          sol_fluxup(i,j) = 0.0
          sol_fluxdn(i,j) = 0.0
          sol_direct(i,j) = 0.0
 2100   continue

        total_sol_fluxup(j) = 0.0
        total_sol_fluxdn(j) = 0.0
        total_sol_direct(j) = 0.0
        total_sol_diff(j) = 0.0

 2000 continue
          

      do 2300 i = 1, n_band
        if (i .eq. 1) then
          nstart = 1
          nend = nsub_band(i)
        else
          nstart = nend + 1
          nend = nend + nsub_band(i)
        endif
        
        do 2400 j = 1, jlayers

          do 2500 n = nstart, nend
     
            sol_fluxup(i,j) = sol_fluxup(i,j) + sol_fluxup_sub(n,j)
            sol_fluxdn(i,j) = sol_fluxdn(i,j) + sol_fluxdn_sub(n,j)
            sol_direct(i,j) = sol_direct(i,j) + sol_direct_sub(n,j)

 2500     continue
 2400   continue
 2300 continue

      do 2600 j = 1, jlayers

        do 2700 i = 1, n_band

          total_sol_fluxup(j) = total_sol_fluxup(j) + sol_fluxup(i,j)
          total_sol_fluxdn(j) = total_sol_fluxdn(j) + sol_fluxdn(i,j)
          total_sol_direct(j) = total_sol_direct(j) + sol_direct(i,j)

 2700   continue

        total_sol_diff(j) = total_sol_fluxdn(j) - total_sol_direct(j)

 2600 continue

      do 2800 i = 1, n_band
        sol_absop(i) = sol_fluxdn(i,1) - sol_fluxup(i,1)
     &               - (sol_fluxdn(i,jlayers) - 
     &                  sol_fluxup(i,jlayers))
 2800 continue

      solnet = total_sol_fluxdn(1) - total_sol_fluxup(1)
     &               - (total_sol_fluxdn(jlayers) - 
     &                  total_sol_fluxup(jlayers))



c outpus to files
      fill_ir = 0.0
      xirdown = 0.0
      xirup = 0.0
      
      open(unit=12,file='../../results/raprad.sw.flx',status='unknown')
      open(unit=13,file='../../results/raprad.aflx',status='unknown')
      open(unit=14,file='../../results/raprad.bflx',status='unknown')



      write(12,527)
      write(12,530) solnet, xirdown, xirup, u0, thz, gmt_time
      write(12,660)

      do 3000 j = 1, jlayers
        jj = (jlayers + 1) - j

c Changing units to
c        alt  m -> km, press Pa -> mb
   
        alt_km = alt(jj) / 1000.0
        press_mb = press(jj) / 100.0

        write(12,665) j, alt_km, press_mb, total_sol_fluxdn(j),
     &    total_sol_fluxup(j), fill_ir, fill_ir, 
     &    total_sol_diff(j), total_sol_direct(j)
 3000 continue
      close(unit=12)

      do 3100 i = 1, n_band
        write(14,658) wave(i), sol_const(i), 
     &  sol_fluxdn(i,1), sol_fluxup(i,1), sol_fluxdn(i,jlayers),
     &  sol_fluxup(i,jlayers),sol_absop(i), sol_direct(i,jlayers)
 3100 continue
      close(unit=14)

      do 3200 j = 1, jlayers
        do 3300 i = 1, n_band
          write(13,659) wave(i), sol_fluxdn(i,j), sol_fluxup(i,j),
     &      sol_direct(i,j)
 3300 continue
 3200 continue
      close(unit=13)

      stop



  100 format(1x)
  527 FORMAT(' RADOUT:      ',
     +     'SOLNET    XIRDOWN   XIRUP     U0       THETA Z   TIME-GMT')
  530 FORMAT(10X,F8.2,2X,2(F9.3,2X),F8.4,F10.4,3x,f6.1)
  658 format(1x,f7.3,2x,7e16.5)
  659 format(1x,f7.3,2x,3e16.5)
  660 FORMAT(2X,'J',T9,'ALT',T20,'P',T26,'SOL FLUX DN',T40,
     3       'SOL FLUX UP',T55,'IR FLUX UP',T69,'IR FLUX DN',T82,
     4       'SOL DIFFUSE',T97,'SOL DIRECT')
  665 FORMAT(I4,2X,F7.3,2X,F7.1,6(4X,F10.2))
  666 format(2x,'Wave l',7x,'solfx(top)',7x,'Flux dn',
     &       9x,'Flux up',9x,'Absorbed flux by ATM')
  667 format(1x,f7.3,2x,8e16.5)
      end



      
