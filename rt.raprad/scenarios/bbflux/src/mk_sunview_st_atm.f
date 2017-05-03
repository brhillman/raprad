      program mk_sunview_st_atm

c This program is to make input file for rt code;

      real mu0

      integer sunz_number, z_number

      character*5 argnum

      pi = 4.0 * atan(1.0)

      call getarg(1, argnum)
      
      read(argnum, '(i5)') numb
      write(*,111) numb

      sunz_number = 1

      mu0 = real(numb) / 100.0
      sunz_theta = acos(mu0)
      sunz_theta = sunz_theta * 180.0 / pi

      sun2earth_distance = 1.0

      z_number = 2

      a_theta = 0.0
      z_theta = 0.0
      
      open(unit=11,file='rt1d.geometry.sunview.no_header',
     &     status='unknown')

      sunz_number = 1
      z_number = 1

      write(11,120) sunz_number
      write(11,121) sunz_theta
      write(11,122) sun2earth_distance
      write(11,123) gmt_time
      write(11,124) z_number
      write(11,125) z_theta, a_theta 

      close(unit = 11)

  111 format(5x, 'mu0 = ', i5, '/100.0')
  120 format(':', i5, ',')
  121 format(':', f10.5, ',')
  122 format(':', f12.8, ',')
  123 format(':', f10.1, ',')
  124 format(':', i5, ',')
  125 format(':', '(', f10.5, ',', f10.5, ')', ',')

      close(unit=11)


      stop
      end
