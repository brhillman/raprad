      program combine_sw_lw_flx

c This program combines SW raprad.flx and LW raprad.flx file to one file.

      open(unit=10,file='../../results/raprad.sw.flx',status='old')
      open(unit=11,file='../../results/raprad.lw.flx',status='old')
      open(unit=12,file='../../results/raprad.flx',status='unknown')

      read(10,100)
      read(10,*) solnet, dummy1, dummy2, u0, thz, gmt_time
      read(10,100)

      read(11,100)
      read(11,*) dummy1, xirdown, xirup, dummy2, dummy3, dummy4
      read(11,100)

      write(12,527)
      write(12,530) solnet, xirdown, xirup, u0, thz, gmt_time
      write(12,660)



 1000 read(10,*,end=2000) j, alt_km, press_mb, total_sol_fluxdn,
     &    total_sol_fluxup, total_ir_fluxup, total_ir_fluxdn, 
     &    total_sol_diff, total_sol_direct


      read(11,*) j_ir, alt_km_ir, press_mb, dummy1,
     &    dummy2, total_ir_fluxup, total_ir_fluxdn, 
     &    dummy3, dymmy4

      if (j .ne. j_ir) then
        write(*,101)
      end if

      if (alt_km .ne. alt_km_ir) then
        write(*,101)
      end if

      write(12,665) j, alt_km, press_mb, total_sol_fluxdn,
     &    total_sol_fluxup, total_ir_fluxup, total_ir_fluxdn, 
     &    total_sol_diff, total_sol_direct

      go to 1000

 2000 close(unit=10)
      close(unit=11)
      close(unit=12)

      stop

  100 format(1x)
  101 format('Layers are different in SW and LW files')
  527 FORMAT(' RADOUT:      ',
     +     'SOLNET    XIRDOWN   XIRUP     U0       THETA Z   TIME-GMT')
  530 FORMAT(10X,F8.2,2X,2(F9.3,2X),F8.4,F10.4,3x,f6.1)
  660 FORMAT(2X,'J',T9,'ALT',T20,'P',T26,'SOL FLUX DN',T40,
     3       'SOL FLUX UP',T55,'IR FLUX UP',T69,'IR FLUX DN',T82,
     4       'SOL DIFFUSE',T97,'SOL DIRECT')
  665 FORMAT(I4,2X,F7.3,2X,F7.1,6(4X,F10.2))

      end


