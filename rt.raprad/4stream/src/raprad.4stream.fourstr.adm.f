      subroutine fourstr
     &  (solconst,nlayers,amu0,omg0in,pmomts,surfalb,
     &   taulam,gau_wt,dis2sun, WVNMLO, WVNMHI)

c need to add radcorr

c
c  driver for 4-stream multiple scattering calculation.
c
c All variables are explained in fourstr_vasriables.dat

c Inputs
c dis2sun      The distance to the sun in astronomical units.
c solconst     The direct beam irradiance at the TOA
c nlayers      Total number of layers in the atmosphere
c amu0         Cosine of the solar zenith angle 
c omg0in(j)    Single scattering albedo for jth layer (j=1,nlayers)
c pmomts(lx(l),j)  (l=0,4) (j=1,nalyers)
c              Coefficients of Legendre polynominals 
c              pmomts(0)=1, and pmomts(1)=3*g
c              pmomts(4) is used for the delta approximation
c surfalb      Surface albedo
c taulam(j)    Optical depth of jth layer (j=1,nlayers)
c iprint       print output flug =0, no outputs

c Output variables

c sol_fluxdn(j)  Downward shortwave irradinace at j level 
c              (W m^-2) (j=1,ntaujs)
c              j level is the top boundary of jth layer
c              the solar constant at TOA is 1.0 so to obtain W m^-2, need to
c              multiply by the solar constant at TOA
c sol_fluxup(j)  Upwnward shortwave irradinace at j level 
c              (W m^-2) (j=1,ntaujs)
c sol_direct(j)  Direct beam irradiance (W m^-2) (j=1,ntaujs)


c ivert  = maximum number of layers;

      parameter (ivert=200)
      parameter (ilayer=ivert+1)  
      parameter(maxumu = 150, maxcmu = 150)
      parameter(maxphi = 100)
                     


      implicit double precision(a-h, o-z)

      dimension fup(ilayer), fdn(ilayer), direct(ilayer)
      dimension sol_fluxup(ilayer), sol_fluxdn(ilayer)
      dimension sol_direct(ilayer)

      real omg0in(*), pmomts(5,*), taulam(*) 
      real amu0, surfalb, WVNMLO, WVNMHI
      real tot_tau, f_dn_black, dir_sfc
      real diff_amu0
      real*8 gau_wt

      real radiance( MAXUMU, MAXPHI )
      real radiance_sfc( MAXUMU, MAXPHI )
      real radiance_trans( MAXUMU, MAXPHI )
      real radiance_toa( MAXUMU, MAXPHI )
      real radiance_toa2( MAXUMU, MAXPHI )
      real dir_radiance_toa( MAXUMU, MAXPHI )
      real radiance_toa_new( MAXUMU, MAXPHI )
      real radiance_toa2_new( MAXUMU, MAXPHI )
      real dir_radiance_toa_new( MAXUMU, MAXPHI )
      real umu(MAXUMU),  phi(MAXPHI), sum_albedo(MAXPHI)
      real BDREF_F, dphi
      real sum_radiance_up(MAXPHI)

      real*8 umu_w(MAXUMU)
      real*8 gauqdt(maxcmu), gau_wt_stream(maxcmu)

c ----------------------------------------------------------------------
c setting up the angle
      pi = 4.0 * atan(1.0)
      numu = 102
      nphi = 61
      numu_half = numu / 2

      call gauleg(0, 1.0d0, gauqdt, gau_wt_stream, numu_half - 1)

      umu(1) = -1.0
      umu(numu) = 1.0
      umu_w(1) = 0.0
      umu_w(numu) = 0.0

      do 1500 i = 2, numu_half
        i_inv = numu_half + 1 - i
        umu(i) = -real(gauqdt(i_inv))
        umu_w(i) = real(gau_wt_stream(i_inv))
 1500 continue

      do 1505 i = numu_half+1, numu-1
        i_for = i - numu_half
        umu(i) = real(gauqdt(i_for))
        umu_w(i) = real(gau_wt_stream(i_for))
 1505 continue

      do 1520 i = 1, nphi
        if (i .eq. 1) then      
          phi(1) = 0.0
        else
          phi(i) = (real(i - 1) * 180.0) / real(nphi - 1)
        end if
 1520 continue


c -------------------------------------------------------------
c computineg surface albedo based on brdf function
      salbedo = 0.0
      do 4300 i = 1, numu / 2
        sum_albedo(i) = 0.0
        do 4400 j = 1, nphi

          dphi = pi * phi(j) / 180.0
          call bdref(WVNMLO, WVNMHI, -umu(i), amu0, dphi, BDREF_F)
          if (j .eq. 1 .or. j .eq. nphi) then
            sum_albedo(i) = sum_albedo(i) 
     &                    + 0.5 * BDREF_F * pi / (nphi - 1)
          else
            sum_albedo(i) = sum_albedo(i) 
     &                    + BDREF_F * pi / (nphi - 1)
          end if
 4400   continue
        salbedo = salbedo - umu(i) * umu_w(i) * sum_albedo(i)
 4300 continue

c      surfalb = 2.0 * salbedo / pi

c test---------------------------------------------------

      if (iopen .eq. 0) then

        open(unit=40,file='out.diagnostic.fourstr.dat',
     &     status='unknown')

        open(unit=41,file=
     &    '../../results/raprad.sw.out',status='unknown')

        open(unit=42,file=
     &    '../../results/disort_radiance/disort.radiance.out',
     &    status='old')

        open(unit=43,
     &    file='../../results/disort.sw.radiance.sfc_contribution.out',
     &    status='unknown')

        open(unit=51,file=
     &    '../../results/disort_radiance/disort.radiance.diff.1.out',
     &    status='old')

        open(unit=52,file=
     &    '../../results/disort_radiance/disort.radiance.diff.2.out',
     &    status='old')

        open(unit=53,file=
     &    '../../results/disort_radiance/disort.radiance.diff.3.out',
     &    status='old')

        open(unit=54,file=
     &    '../../results/disort_radiance/disort.radiance.diff.4.out',
     &    status='old')

        open(unit=55,file=
     &    '../../results/disort_radiance/disort.radiance.diff.5.out',
     &    status='old')

        open(unit=56,file=
     &    '../../results/disort_radiance/disort.radiance.diff.6.out',
     &    status='old')

        open(unit=57,file=
     &    '../../results/disort_radiance/disort.radiance.diff.7.out',
     &    status='old')

        open(unit=58,file=
     &    '../../results/disort_radiance/disort.radiance.diff.8.out',
     &    status='old')

        open(unit=59,file=
     &    '../../results/disort_radiance/disort.radiance.diff.9.out',
     &    status='old')

        open(unit=60,file=
     &    '../../results/disort_radiance/disort.radiance.diff.10.out',
     &    status='old')

      end if

      iopen = 1

c     write(*,500) solconst
c     write(*,501) nlayers
c     write(*,500) amu0
c     write(*,502) (omg0in(i), i=1,3)
c     write(*,503) (pmomts(j,1), j=1,5)
c     write(*,503) (pmomts(j,2), j=1,5)
c     write(*,500) surfalb
c     write(*,502) (taulam(i), i=1,3)

  500 format(f15.5)
  501 format(i5)
  502 format(3f15.5)
  503 format(5f15.5)

      write(40,301) solconst, nlayers, amu0, surfalb

      do 60 j = 1, nlayers
        write(40,300) j, omg0in(j),taulam(j),(pmomts(i,j), i=1,5)
   60 continue


c-------------------------------------------------------


c radcorr corrects the solar constant based on the ephemeris derived
c distance to sun 

c      radcorr = dis2sun*dis2sun

c w0max is the maximum allowed asymmetry parameter - values larger 
c than 0.9999 tend to induce model instabilities.

      w0max = 0.9999
      ntaujs = nlayers + 1

      do 10 j = 1,ntaujs
         fup(j) = 0.0
         fdn(j) = 0.0
 10   continue

       do 30 j = 1,nlayers
         if(omg0in(j).gt.w0max) omg0in(j) = w0max
 30    continue

       call solrad
     &  (nlayers,amu0,omg0in,pmomts,surfalb,taulam,
     &     iprint,fdn,fup,direct)


c adjust top of model domain flux to account for middle atmosphere
c extinction according to two-stream results

c         toaflx = (ck1(i,1)*el2(i,1)+ck2(i,1)*em2(i,1)+
c     1               cmb(i,1)+direct(i,1))/amu0

c correction factor of the sun to the earth distance

         dist_cor = 1.0 / (dis2sun)**2

         do 50 j = 1,ntaujs
            sol_fluxup(j) = fup(j)*solconst*gau_wt*dist_cor
            sol_fluxdn(j) = fdn(j)*solconst*gau_wt*dist_cor
            sol_direct(j) = direct(j)*solconst*gau_wt*dist_cor

            write(41,400) j, sol_fluxup(j), sol_fluxdn(j), sol_direct(j)

 50      continue


c ----------------------------------------------------------------------
c computing toa radiance

      do 2000 iu = 1, numu
        read(42,*) ii, (radiance(iu,j), j = 1, nphi)
 2000 continue

      do 2100 iu = 1, numu / 2
        sum_radiance_up(iu) = 0.0
        do 2200 j = 1, nphi
          radiance_sfc(iu,j) = radiance(iu,j) 
     &                       / solconst * gau_wt * dist_cor
          radiance_trans(iu,j) = radiance(iu,j) 
     &                       / solconst * gau_wt * dist_cor

          if (j .eq. 1 .or. j .eq. nphi) then
            sum_radiance_up(iu) = sum_radiance_up(iu)
     &              + 0.5 * radiance(numu/2+iu,j) * pi / (nphi - 1)
          else
            sum_radiance_up(iu) = sum_radiance_up(iu)
     &                    + radiance(numu/2+iu,j) * pi / (nphi - 1)
          end if

 2200   continue
        f_up_black = f_up_black 
     &         + umu(numu/2+iu) * umu_w(numu/2+iu) * sum_radiance_up(iu)
 2100 continue

      f_up_black = 2.0 * f_up_black 
     &           / solconst * gau_wt * dist_cor

      dir_sfc = direct(ntaujs)

      tot_tau = 0.0
      do 2300 layer = 1, nlayers
        tot_tau = tot_tau + taulam( layer )
 2300 continue

        idiff = 0
        isave = 1

        call radiance_integration
     &       (radiance_sfc, radiance_trans, dir_sfc, amu0, numu, nphi, 
     &        umu, umu_w, phi, tot_tau,
     &        WVNMLO, WVNMHI, idiff, isave, 
     &        radiance_toa, radiance_toa2,
     &        dir_radiance_toa, f_dn_black, salbedo)


      gain =  (fdn(ntaujs) - direct(ntaujs)) / f_dn_black

      write(*,*) 'gain = ', gain

      do 3000 iu = 1, numu / 2
        do 3100 j = 1, nphi
c          radiance_toa_new(iu,j) = gain * solconst * gau_wt * dist_cor
c     &                          * radiance_toa(iu,j)

          radiance_toa_new(iu,j) = 0.0

          radiance_toa2_new(iu,j) = gain * solconst * gau_wt * dist_cor
     &                           * radiance_toa2(iu,j)

          dir_radiance_toa_new(iu,j) = 
     &                      fdn(ntaujs) * solconst * gau_wt * dist_cor
     &                    * dir_radiance_toa(iu,j)
 3100   continue
 3000 continue

c ------------------------------------------------------------------
c computing diffuse-diffuse transmission at all other solar zenith angles

      idiff = 1

      if (isave .eq. 0) then
        do 6000 nfile = 1, 10
          diff_amu0 = real(nfile) / 10.0

          iunit = 50 + nfile

          do 6100 iu = 1, numu
            read(iunit,*) ii, (radiance(iu,j), j = 1, nphi)

            do 6200 j = 1, nphi
              radiance_trans(iu,j) = radiance(iu,j) 
     &                           / solconst * gau_wt * dist_cor
 6200       continue
 6100     continue


          call radiance_integration
     &         (radiance_sfc, radiance_trans, dir_sfc, diff_amu0, 
     &          numu, nphi, umu, umu_w, phi, tot_tau,
     &          WVNMLO, WVNMHI, idiff, isave, 
     &          radiance_toa, radiance_toa2,
     &          dir_radiance_toa, f_dn_black, salbedo)

          do 6300 iu = 1, numu / 2
            do 6400 j = 1, nphi

              radiance_toa_new(iu,j) = radiance_toa_new(iu,j) + 
     &                    gain * solconst * gau_wt * dist_cor
     &                    * radiance_toa(iu,j) / 10.0
 6400       continue
 6300     continue
 6000   continue
      end if

c ------------------------------------------------------------------
c output results

      do 5000 iu = 1, numu / 2
        write(43,410) iu, (radiance_toa_new(iu,j), j = 1, nphi)
 5000 continue
      do 5100 iu = 1, numu / 2
        write(43,410) iu, (radiance_toa2_new(iu,j), j = 1, nphi)
 5100 continue
      do 5200 iu = 1, numu / 2
        write(43,410) iu, (dir_radiance_toa_new(iu,j), j = 1, nphi)
 5200 continue


  300 format(i5, 7e15.5)
  301 format(f10.5, i5, 2f10.5)
  400 format(i5,3e15.5)
  410 format(i5,<nphi>e15.5)


      return
      end

 
