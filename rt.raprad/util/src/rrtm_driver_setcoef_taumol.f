      subroutine rrtm_driver_setcoef_taumol
     &  (iband, ig, lay, tbound, tz, tavel, pz, pavel, 
     &   wkl_1, wkl_2, wkl_3, wkl_4, wkl_6, wkl_7, 
     &   coldry, albedo, laytrop, layswtch, laylow,
     &   bi, gau_wt, plankbnd, planklay, planklev)

 
c This subroutine calls two rrtm subroutines, setcoef and taumol.
c Input
c tbound   temperature at the surface
c
c layer is counted from bottom. So the lowest layer is 1 and the top layer is
c n. The bottom boundary of each layer is L-1 so that the lowest boundary 
c starts from 0 and goes upto n.
c 
c tz(L-1) temperature at bottom of layer L-1  -  used by RRTM for 
c          Planck Function Calculation
c
c tavel    layer temperatures (degrees K)
c
c pz(L-1)  pressure at bottom of layer L-1
c
c pavel    layer pressures (mb)

c  wkl_1    column amount of water vapor (molecules/cm**2)
c  wkl_2    column amount of co2         (molecules/cm**2)
c  wkl_3    column amount of ozone       (molecules/cm**2)
c  wkl_4    column amount of n2o         (molecules/cm**2)
c  wkl_5    column amount of carbon monoxide
c           not used                     (molecules/cm**2)
c  wkl_6    column amount of ch4         (molecules/cm**2)
c  wkl_7    column amount of oxygen      (molecules/cm**2)


c
c coldry   column density of dry air (molecules/cm**2)
c
c semiss   surface emissivity, 1.0 - albedo.
c
c jp       the index of the lower (in altitude) of the two appropriate
c          reference pressure levels needed for interpolation
c jt, jt1  the indices of the lower of the two appropriate reference  
c          temperatures needed for interpolation (for pressure    
c          levels JP and JP+1, respectively)
c
c fac00, fac10 fac01, fac11
c          for layer LAY, these are factors that are needed to compute 
c          the interpolation factors that multiply the appropriate 
c          reference k-values. A value of 0 (1) for i,j indicates that 
c          the corresponding factor multiplies reference k-value for the 
c          lower (higher) of the two appropriate temperatures, and 
c          altitudes, respectively.
c
c selffac  scale factor needed to water vapor self-continuum, equals 
c          (water vapor density)/(atmospheric density at 296K and 
c          1013 mb)
c
c forfac    
c
c colh2o, colco2, colo3, coln2o, colch4
c          column amounts of water vapor,carbon dioxide, ozone, 
c          nitrous oxide, methane, respectively (molecules/cm**2)
c          COLH2O(LAY) = 1.E-20 * WKL(1,LAY)
c          COLCO2(LAY) = 1.E-20 * WKL(2,LAY)
c          COLO3(LAY) = 1.E-20 * WKL(3,LAY)
c          COLN2O(LAY) = 1.E-20 * WKL(4,LAY)
c          COLCH4(LAY) = 1.E-20 * WKL(6,LAY)
c          COLO2(LAY) = 1.E-20 * WKL(7,LAY)
c
c
c co2mult  for bands in which carbon dioxide is implemented as a 
c          trace species, this is the factor used to multiply the 
c          band's average CO2 absorption coefficient to get the added
c          contribution to the optical depth relative to 355 ppm.
c
c
c taug     the optical depths per g-value and layer for a band
c
c fracs    fractions needed to compute Planck functions at every layer
c          and g-value
c
c bi       absorption coefficient
c plankbnd   the integrated Planck functions at the surface temperature
c planklay   the integrated Planck functions at the layer temperature
c planklev   the integrated Planck functions at the level temperature



      parameter (nbands = 16)

      real tbound, albedo, semiss
      real*8 coldry, bi, gau_wt 
      real*8 tz, pz, tavel, pavel
      real*8 fac00, fac01, fac10, fac11, selffac, selffrac, forfac
      real*8 wkl_1, wkl_2, wkl_3, wkl_4, wkl_5, wkl_6, wkl_7
      real*8 colh2o, colco2, colo3, coln2o
      real*8 colch4, colo2, co2mult
      real*8 taug, fracs

      real*8 plankbnd, planklay, planklev

      integer iband, ig, lay, jp, jt, jt1
      integer laytrop, layswtch, laylow
      integer ng(nbands), nspa(nbands), nspb(nbands)

      data ng  /16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/
      data nspa /1, 1,10, 9, 9, 1, 9, 1,11, 1, 1, 9, 9, 1, 9, 9/
      data nspb /1, 1, 5, 6, 5, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0/

c      DATA WX /MAXPROD*0.0/
       
c Mlawer's band starts from shortest wavenumber.     
      iband_lw = 16 - iband + 1
      semiss = 1.0 - albedo


      call setcoef 
     &  (iband_lw, ig, lay,
     &   tbound, tz, tavel, pz, pavel,
     &   wkl_1, wkl_2, wkl_3, wkl_4, wkl_5, wkl_6, wkl_7,
     &   semiss, 
     &   jp, jt, jt1, 
     &   fac00, fac01, fac10, fac11, selffac, selffrac, forfac,
     &   indself, colh2o, colco2, colo3, coln2o, colch4, colo2,  
     &   coldry, co2mult,
     &   plankbnd, planklay, planklev)
      


C     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
      if (iband_lw .eq. 1) then
        call taugb1
     &    (ig, lay, jp, jt, jt1, nspa(1), nspb(1), laytrop,
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colh2o, 
     &     taug, fracs)

C     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)
      else if (iband_lw .eq. 2) then
        call taugb2
     &    (ig, lay, jp, jt, jt1, nspa(2), nspb(2), laytrop, 
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colh2o, coldry, 
     &     taug, fracs)

C     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)
      else if (iband_lw .eq. 3) then
        call taugb3
     &    (ig, lay, jp, jt, jt1, nspa(3), nspb(3), laytrop, 
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colh2o, colco2, coln2o,
     &     taug, fracs)

C     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)
      else if (iband_lw .eq. 4) then
        call taugb4
     &    (ig, lay, jp, jt, jt1, nspa(4), nspb(4), laytrop, 
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colh2o, colco2, colo3,
     &     taug, fracs)

C     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
      else if (iband_lw .eq. 5) then

c     CCL4
        wx_1 = 0.0
        call taugb5
     &    (ig, lay, jp, jt, jt1, nspa(5), nspb(5), laytrop, 
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colh2o, colco2, colo3,
     &     wx_1,
     &     taug, fracs)


C     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)
      else if (iband_lw .eq. 6) then

c     CFC11, CFC12
        wx_2 = 0.0
        wx_3 = 0.0
        call taugb6
     &    (ig, lay, jp, jt, jt1, nspa(6), nspb(6), laytrop, 
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colh2o, co2mult,
     &     wx_2, wx_3,
     &     taug, fracs)

C     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)
      else if (iband_lw .eq. 7) then
        call taugb7
     &    (ig, lay, jp, jt, jt1, nspa(7), nspb(7), laytrop, 
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colh2o, colo3,
     &     co2mult, 
     &     taug, fracs)

C     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)
      else if (iband_lw .eq. 8) then

c     CFC12, CFC22
        wx_3 = 0.0
        wx_4 = 0.0
        call taugb8
     &    (ig, lay, jp, jt, jt1, nspa(8), nspb(8), laytrop, layswtch,
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colh2o, 
     &     co2mult, coln2o, colo3, wx_3, wx_4,
     &     taug, fracs)

C     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)
      else if (iband_lw .eq. 9) then
        call taugb9
     &    (ig, lay, jp, jt, jt1, nspa(9), nspb(9), 
     &     laytrop, layswtch, laylow,
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, 
     &     colh2o, colch4, coln2o,
     &     taug, fracs)

C     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)
      else if (iband_lw .eq. 10) then
        call taugb10
     &    (ig, lay, jp, jt, jt1, nspa(10), nspb(10), 
     &     laytrop,
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colh2o,
     &     taug, fracs)

C     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)
      else if (iband_lw .eq. 11) then
        call taugb11
     &    (ig, lay, jp, jt, jt1, nspa(11), nspb(11), 
     &     laytrop,
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colh2o,
     &     taug, fracs)

C     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)
      else if (iband_lw .eq. 12) then
        call taugb12
     &    (ig, lay, jp, jt, jt1, nspa(12), nspb(12), 
     &     laytrop,
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colh2o, colco2,
     &     taug, fracs)

C     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)
      else if (iband_lw .eq. 13) then
        call taugb13
     &    (ig, lay, jp, jt, jt1, nspa(13), nspb(13), 
     &     laytrop,
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colh2o, coln2o,
     &     taug, fracs)

C     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)
      else if (iband_lw .eq. 14) then
        call taugb14
     &    (ig, lay, jp, jt, jt1, nspa(14), nspb(14), 
     &     laytrop,
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colco2,
     &     taug, fracs)

C     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)
      else if (iband_lw .eq. 15) then
        call taugb15
     &    (ig, lay, jp, jt, jt1, nspa(15), nspb(15), 
     &     laytrop,
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, 
     &     colh2o, colco2, coln2o,
     &     taug, fracs)

C     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)
      else if (iband_lw .eq. 16) then
        call taugb16
     &    (ig, lay, jp, jt, jt1, nspa(16), nspb(16), 
     &     laytrop,
     &     fac00, fac01, fac10, fac11,
     &     indself, selffac, selffrac, forfac, colh2o, colch4,
     &     taug, fracs)


      end if

      bi = taug
      gau_wt = fracs

!      if (bi .lt. 0.0) then
!        write(*,*) iband_lw, bi, co2mult, coln2o, colo3, wx_3, wx_4
!      end if

      return
      end
