MODULE usrdef_sbc
   !!======================================================================
   !!                     ***  MODULE  usrdef_sbc  ***
   !!
   !!                     ===  GYRE configuration  ===
   !!
   !! User defined :   surface forcing of a user configuration
   !!======================================================================
   !! History :  4.0   ! 2016-03  (S. Flavoni, G. Madec)  user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usrdef_sbc    : user defined surface bounday conditions in GYRE case
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran    !

   !! nch : forcing for NW2
   USE cubic_interp    ! interpolation of piecewise cubic polynomials
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE ice, ONLY       : at_i_b, a_i_b, t_su, jpl, hm_i
   USE icethd_dh       ! for CALL ice_thd_snwblow
   !! end

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usrdef_sbc_oce       ! routine called in sbcmod module
   PUBLIC   usrdef_sbc_ice_tau   ! routine called by icestp.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx   ! routine called by icestp.F90 for ice thermo

   !! * Substitutions

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_sbc.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usrdef_sbc_oce( kt,Kbb )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE usrdef_sbc  ***
      !!              
      !! ** Purpose :   provide at each time-step the GYRE surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   analytical seasonal cycle for GYRE configuration.
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.   
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!
      !! Reference : Hazeleger, W., and S. Drijfhout, JPO, 30, 677-695, 2000.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ,Kbb  ! ocean time step
      !!
      INTEGER  ::   ji, jj                 ! dummy loop indices
      INTEGER  ::   zyear0                 ! initial year 
      INTEGER  ::   zmonth0                ! initial month
      INTEGER  ::   zday0                  ! initial day
      INTEGER  ::   zday_year0             ! initial day since january 1st
      REAL(wp) ::   ztau     , ztau_sais   ! wind intensity and of the seasonal cycle
      REAL(wp) ::   ztime                  ! time in hour
      REAL(wp) ::   ztimemax , ztimemin    ! 21th June, and 21th decem. if date0 = 1st january
      REAL(wp) ::   ztimemax1, ztimemin1   ! 21th June, and 21th decem. if date0 = 1st january
      REAL(wp) ::   ztimemax2, ztimemin2   ! 21th June, and 21th decem. if date0 = 1st january
      REAL(wp) ::   ztaun                  ! intensity
      REAL(wp) ::   zemp_s, zemp_n, zemp_sais, ztstar
      REAL(wp) ::   zcos_sais1, zcos_sais2, ztrp, zconv, t_star
      REAL(wp) ::   zsumemp, zsurf
      REAL(wp) ::   zrhoa  = 1.22         ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3       ! drag coefficient
      REAL(wp) ::   ztx, zty, zmod, zcoef ! temporary variables
      REAL(wp) ::   zyydd                 ! number of days in one year
      !!---------------------------------------------------------------------
      ! nch : forcing for NW2

      REAL(wp) :: ztrel     !relativity time in a year
      REAL(wp) :: zyrel     !relativity y in domain
      integer, parameter :: wind_nsize = 8, sst_nsize =5, emp_nsize = 7 
      real(wp), dimension(wind_nsize) :: ytau,ytau_t   ! latitude of interpolation
      real(wp), dimension(wind_nsize) :: taud   ! wind intensity of interpolation
      real(wp), dimension(wind_nsize) :: taud_t   ! wind intensity of interpolation
      real(wp), dimension(sst_nsize)  :: ysst,ysst_t   ! latitude of interpolation
      real(wp), dimension(sst_nsize)  :: sstd,sstd_t   ! sst of interpolation
      real(wp), dimension(emp_nsize)  :: yemp,yemp_t   ! latitude of interpolation
      real(wp), dimension(emp_nsize)  :: empd,empd_t   ! emp of interpolation
      
      !data ytau /-75.,-65.,-45.,-15.,0.,15.,45.,70./
      !data taud /-0.004,0.,.25,-0.1,-.02,-.1,.1,0./
      !data ysst /-140.,-68.,0.,140./
      !data sstd /-30.,-3.,28.3,-30./
      !data yemp /-75.,-55.,-20.,5.,30.,50.,80./
      ![0,-0.5,1,-1.2,1,-0.5,0]
      !data empd /0.,-0.5,1.1095246481,-1.2,1.,-0.5,0./

     ysst = (/-140.,-78.,-65.,0.,140./) 
      sstd = (/-30.,-3.4,-3.,28.3,-30./)
      yemp = (/-75.,-55.,-20.,5.,30.,50.,80./)
      empd = (/0.,-0.5,1.1095246481,-1.2,1.,-0.5,0./)
       ytau = (/-80.,-70.,-50.,-15.,0.,15.,45.,70./)
      taud = (/-0.017,-0.017,.21,-0.1,-.02,-.1,.1,0./)
      
      !--------------------------
      zyydd = REAL(nyear_len(1),wp)

      ! ---------------------------- !
      !  heat and freshwater fluxes  !
      ! ---------------------------- !
      !same temperature, E-P as in HAZELEGER 2000

      zyear0     =   ndate0 / 10000                             ! initial year
      zmonth0    = ( ndate0 - zyear0 * 10000 ) / 100            ! initial month
      zday0      =   ndate0 - zyear0 * 10000 - zmonth0 * 100    ! initial day betwen 1 and 30
      zday_year0 = ( zmonth0 - 1 ) * 30.+zday0                  ! initial day betwen 1 and 360

      ! current day (in hours) since january the 1st of the current year
      ztime = REAL( kt ) * rdt / (rmmss * rhhmm)   &       !  total incrementation (in hours)
         &      - (nyear  - 1) * rjjhh * zyydd             !  minus years since beginning of experiment (in hours)

      ztimemax1 = ((5.*30.)+21.)* 24.                      ! 21th june     at 24h in hours
      ztimemin1 = ztimemax1 + rjjhh * zyydd / 2            ! 21th december        in hours
      ztimemax2 = ((6.*30.)+21.)* 24.                      ! 21th july     at 24h in hours
      ztimemin2 = ztimemax2 - rjjhh * zyydd / 2            ! 21th january         in hours
      !                                                    ! NB: rjjhh * zyydd / 4 = one seasonal cycle in hours

      ! amplitudes
      zemp_S    = 0.7       ! intensity of COS in the South
      zemp_N    = 0.8       ! intensity of COS in the North
      zemp_sais = 0.1
      zTstar    = 28.3      ! intemsity from 28.3 a -5 deg

      ! 1/2 period between 21th June and 21th December and between 21th July and 21th January
      zcos_sais1 = COS( (ztime - ztimemax1) / (ztimemin1 - ztimemax1) * rpi ) 
      zcos_sais2 = COS( (ztime - ztimemax2) / (ztimemax2 - ztimemin2) * rpi )

      ztrp= - 40.e0        ! retroaction term on heat fluxes (W/m2/K)
      zconv = 3.16e-5      ! convertion factor: 1 m/yr => 3.16e-5 mm/s

      ! nch : forcing for NW2
      ztrel = ztime/(rjjhh * zyydd)
      DO ji=1,sst_nsize
         ysst_t(ji) = ysst(ji)+3.*COS(2.*rpi*(ztrel-0.558))
         sstd_t(ji) = sstd(ji)
         if (ji==2 .or. ji==3) then
            sstd_t(ji) = sstd(ji)-1.6*COS(2.*rpi*(ztrel-0.558))
         endif
      END DO
      DO ji = 1, emp_nsize
         if (ji<4) then
            empd_t(ji) = empd(ji)*(1.-0.1/0.8*COS(2.*rpi*(ztrel-0.475)-rpi))
         else
            empd_t(ji) = empd(ji)*(1.-0.1/0.8*COS(2.*rpi*(ztrel-0.475)))
         endif
      END DO

       DO jj = 1, jpj
         DO ji = 1, jpi
            t_star = itau(ysst_t,sstd_t,gphit(ji,jj))
            ! 23.5 deg : tropics
            qsr (ji,jj) =  230. * COS( 3.1415 * ( gphit(ji,jj) - 23.5 * zcos_sais1 ) / ( 0.9 * 180. ) )
            qsr (ji,jj) = MAX( qsr(ji,jj), 0. )
            qns (ji,jj) = ztrp * ( ts(ji,jj,1,jp_tem,1) - t_star ) - qsr(ji,jj)
            
            ! jj是3的倍数
            ! if (ji==5 .and. mod(jj,5)==0) then
            !   if( mpprank == 1 .and. jpj<60 .and. .FALSE. ) then !mppsize == 2 .and.
            !      open(unit=78,file='a1p_out.txt',position='APPEND',action='write',form='formatted')
            !      write(78,*) 'kt = ', kt
             !     write(78,*) ' nch debug jj,ji      = ',jj,ji
            !      write(78,*) ' nch debug gphit      = ',gphit(ji,jj)
            !      write(78,*) ' nch debug qsr        = ',qsr(ji,jj)
             !     write(78,*) ' nch debug qns        = ',qns(ji,jj)
            !      write(78,*) ' nch debug t_star     = ',t_star
            !      write(78,*) ' nch debug tsb        = ',tsb(ji,jj,1,jp_tem)
            !      close(78)
            !   endif
            !ENDIF

            !IF( gphit(ji,jj) >= 14.845 .AND. 37.2 >= gphit(ji,jj) ) THEN    ! zero at 37.8 deg, max at 24.6 deg
             !  emp  (ji,jj) =   zemp_S * zconv   &
              !    &         * SIN( rpi / 2 * (gphit(ji,jj) - 37.2) / (24.6 - 37.2) )  &
               !   &         * ( 1 - zemp_sais / zemp_S * zcos_sais1)
            !ELSE
             !  emp (ji,jj) =  - zemp_N * zconv   &
              !    &         * SIN( rpi / 2 * (gphit(ji,jj) - 37.2) / (46.8 - 37.2) )  &
               !   &         * ( 1 - zemp_sais / zemp_N * zcos_sais1 )
            !ENDIF
            emp (ji,jj) = itau(yemp,empd_t,gphit(ji,jj)) * zconv
            

         END DO
      END DO
      ! end

      !zsumemp = GLOB_SUM( 'usrdef_sbc', emp  (:,:)   ) 
      !zsurf   = GLOB_SUM( 'usrdef_sbc', tmask(:,:,1) ) 
      !zsumemp = zsumemp / zsurf         ! Default GYRE configuration

      ! freshwater (mass flux) and update of qns with heat content of emp
      !emp (:,:) = emp(:,:) - zsumemp * tmask(:,:,1)        ! freshwater flux (=0 in domain average)
      sfx (:,:) = 0.0_wp                                   ! no salt flux
      qns (:,:) = qns(:,:) - emp(:,:) * sst_m(:,:) * rcp   ! evap and precip are at SST


      ! ---------------------------- !
      !       momentum fluxes        !
      ! ---------------------------- !
      ! same wind as in Wico
      !test date0 : ndate0 = 010203
      zyear0  =   ndate0 / 10000
      zmonth0 = ( ndate0 - zyear0 * 10000 ) / 100
      zday0   =   ndate0 - zyear0 * 10000 - zmonth0 * 100
      !Calculates nday_year, day since january 1st
      zday_year0 = (zmonth0-1)*30.+zday0

      !accumulates days of previous months of this year
      ! day (in hours) since january the 1st
      ztime = REAL( kt ) * rdt / (rmmss * rhhmm)  &  ! incrementation in hour
         &     - (nyear - 1) * rjjhh * zyydd          !  - nber of hours the precedent years
      ztimemax = ((5.*30.)+21.)* 24.               ! 21th june     in hours
      ztimemin = ztimemax + rjjhh * zyydd / 2      ! 21th december in hours
      !                                            ! NB: rjjhh * zyydd / 4 = 1 seasonal cycle in hours
      !!! mean intensity at 0.105 ; srqt(2) because projected with 45deg angle
      !ztau = 0.105 / SQRT( 2. )
      !! seasonal oscillation intensity
      !ztau_sais = 0.015
      !ztaun = ztau - ztau_sais * COS( (ztime - ztimemax) / (ztimemin - ztimemax) * rpi )
      !DO jj = 1, jpj
      !DO ji = 1, jpi
      !! domain from 15deg to 50deg and 1/2 period along 14deg
      !! so 5/4 of half period with seasonal cycle
      !utau(ji,jj) = - ztaun * SIN( rpi * (gphiu(ji,jj) - 15.) / (29.-15.) )
      !vtau(ji,jj) =   ztaun * SIN( rpi * (gphiv(ji,jj) - 15.) / (29.-15.) )
      !END DO
      !END DO

      !! nch : forcing for NW2
      ztrel = ztime/(rjjhh * zyydd)
      DO ji=1,wind_nsize
         ytau_t(ji) = ytau(ji)+3.*COS(2.*rpi*ztrel-0.79-rpi)
         if (ji<5) then
            taud_t(ji) = taud(ji)/144.*(COS(2.*rpi*ztrel-0.79-rpi)+12.)**2.
         else
            taud_t(ji) = taud(ji)/144.*(COS(2.*rpi*ztrel-0.79)+12.)**2.
         endif

      END DO
      DO jj = 1, jpj
         DO ji = 1, jpi
            utau(ji,jj) = itau(ytau_t,taud_t,gphiu(ji,jj))
            vtau(ji,jj) = 0.
            !if (ji==30) then
               !write(numout,*)'                    nch debug utau        = ',utau(ji,jj)
            !ENDIF
         END DO
      END DO

      !!!!!!!!!!!!!!!!!!!!!!!!!

      ! module of wind stress and wind speed at T-point
      zcoef = 1. / ( zrhoa * zcdrag ) 
      DO jj = 2, jpj-1
         DO ji = 2, jpi-1   ! vect. opt.
            ztx = utau(ji-1,jj  ) + utau(ji,jj) 
            zty = vtau(ji  ,jj-1) + vtau(ji,jj) 
            zmod = 0.5 * SQRT( ztx * ztx + zty * zty )
            taum(ji,jj) = zmod
            wndm(ji,jj) = SQRT( zmod * zcoef )
         END DO
      END DO
      !CALL lbc_lnk_multi( 'usrdef_sbc', taum(:,:), 'T', 1. , wndm(:,:), 'T', 1. )

      ! ---------------------------------- !
      !  control print at first time-step  !
      ! ---------------------------------- !
      IF( kt == nit000 .AND. lwp ) THEN 
         WRITE(numout,*)
         WRITE(numout,*)'usrdef_sbc_oce : NW2 analytical surface fluxes for NW2 configuration'               
         WRITE(numout,*)'~~~~~~~~~~~ ' 
         WRITE(numout,*)'           nyear      = ', nyear
         WRITE(numout,*)'           nmonth     = ', nmonth
         WRITE(numout,*)'           nday       = ', nday
         WRITE(numout,*)'           nday_year  = ', nday_year
         WRITE(numout,*)'           ztime      = ', ztime
         WRITE(numout,*)'           ztimemax   = ', ztimemax
         WRITE(numout,*)'           ztimemin   = ', ztimemin
         WRITE(numout,*)'           ztimemax1  = ', ztimemax1
         WRITE(numout,*)'           ztimemin1  = ', ztimemin1
         WRITE(numout,*)'           ztimemax2  = ', ztimemax2
         WRITE(numout,*)'           ztimemin2  = ', ztimemin2
         WRITE(numout,*)'           zyear0     = ', zyear0
         WRITE(numout,*)'           zmonth0    = ', zmonth0
         WRITE(numout,*)'           zday0      = ', zday0
         WRITE(numout,*)'           zday_year0 = ', zday_year0
         WRITE(numout,*)'           zyydd      = ', zyydd
         WRITE(numout,*)'           zemp_S     = ', zemp_S
         WRITE(numout,*)'           zemp_N     = ', zemp_N
         WRITE(numout,*)'           zemp_sais  = ', zemp_sais
         WRITE(numout,*)'           zTstar     = ', zTstar
         WRITE(numout,*)'           zsumemp    = ', zsumemp
         WRITE(numout,*)'           zsurf      = ', zsurf
         WRITE(numout,*)'           ztrp       = ', ztrp
         WRITE(numout,*)'           zconv      = ', zconv
         WRITE(numout,*)'           ndastp     = ', ndastp
         WRITE(numout,*)'           adatrj     = ', adatrj
      ENDIF
      !
   END SUBROUTINE usrdef_sbc_oce


   SUBROUTINE usrdef_sbc_ice_tau( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      IF( kt==nit000 .AND. lwp)   WRITE(numout,*)' usrdef_sbc_ice : NW2 case'
      !
      utau_ice(:,:) = utau(:,:)
      vtau_ice(:,:) = vtau(:,:)
   END SUBROUTINE usrdef_sbc_ice_tau


   SUBROUTINE usrdef_sbc_ice_flx( kt, phs, phi )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_sbc_flx  ***
      !!
      !! ** Purpose : provide surface boundary condition for sea ice (flux)
      !!
      !! ** Action  : It provides the following fields used in sea ice model:
      !!                emp_oce , emp_ice                        = E-P over ocean and sea ice                    [Kg/m2/s]
      !!                sprecip                                  = solid precipitation                           [Kg/m2/s]
      !!                evap_ice                                 = sublimation                                   [Kg/m2/s]
      !!                qsr_tot , qns_tot                        = solar & non solar heat flux (total)           [W/m2]
      !!                qsr_ice , qns_ice                        = solar & non solar heat flux over ice          [W/m2]
      !!                dqns_ice                                 = non solar  heat sensistivity                  [W/m2]
      !!                qemp_oce, qemp_ice, qprec_ice, qevap_ice = sensible heat (associated with evap & precip) [W/m2]
      !!            + some fields that are not used outside this module: 
      !!                qla_ice                                  = latent heat flux over ice                     [W/m2]
      !!                dqla_ice                                 = latent heat sensistivity                      [W/m2]
      !!                tprecip                                  = total  precipitation                          [Kg/m2/s]
      !!                alb_ice                                  = albedo above sea ice
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phs    ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phi    ! ice thickness

      !!
      REAL(wp) ::   zfr1, zfr2, qns_star  ,cldf_ice               ! local variables
      REAL(wp), DIMENSION(jpi,jpj) ::   zsnw   ! snw distribution after wind blowing
      integer ::   ji,jj,jk                      ! dummy loop indices
      !!---------------------------------------------------------------------
      !
      IF( kt==nit000 .AND. lwp)   WRITE(numout,*)' usrdef_sbc_ice : NW2 case'
cldf_ice=0.81
      !
      ! ocean variables (renaming)
      qsr_oce (:,:)   = qsr(:,:)
      emp_oce (:,:)   = emp(:,:)
      qns_oce (:,:)   = qns(:,:)

      ! ice variables
      do ji = 1, size( qsr_ice, 3)
         qsr_ice (:,:,ji) =  qsr(:,:) 
         !qns_ice (:,:,ji) = qns(:,:) 
         
      enddo

      DO jj = 1, jpj
         DO ji = 1, jpi
            DO jk = 1,size( qsr_ice, 3)
            !qns_ice (ji,jj,jk) = qns(ji,jj) -2.5 exp04
            zfr2 = EXP(  - (hm_i(ji,jj)/0.5) ) 
            zfr1 = - 40.e0  * ( ts(ji,jj,1,jp_tem,1) + 1.9 ) - qsr(ji,jj)
            qns_star  = zfr2*(qns(ji,jj)-zfr1)+zfr1
            qns_ice (ji,jj,jk) = qns_star 
            if (qns(ji,jj)>zfr1) qns_ice (ji,jj,jk) = qns(ji,jj) !qns高温
            ! jj是3的倍数
            ! if (ji==17 .and. mod(jj,5)==0) then
            !   write(numout,*)'                    nch debug jj,ji      = ',jj,ji
            !   write(numout,*)'                    nch debug gphit      = ',gphit(ji,jj)
            !   write(numout,*)'                    nch debug qsr ice     = ',qsr_ice(ji,jj,1)
            !   write(numout,*)'                    nch debug qns ice    = ',qns_ice(ji,jj,1)
            !ENDIF
            
         END DO
         END DO
      END DO
      ! end


      emp_ice (:,:) = 0.

      
      qemp_ice (:,:)   = 0._wp  !  - emp_ice(:,:) * SUM(t_su(:,:,:),dim=3)/jpl * rcpi
      sprecip (:,:)   = 0._wp   ! uniform value for snow precip
      evap_ice(:,:,:) = 0._wp   ! uniform value for sublimation
      qevap_ice(:,:,:) =   0._wp
      

      ! ice fields deduced from above
      qemp_oce (:,:)   = - emp_oce(:,:) * sst_m(:,:) * rcp
      qprec_ice(:,:)   =   rhos * ( sst_m(:,:) * rcpi - rLfus ) * tmask(:,:,1) !  in J/m3


      ! total fluxes
      emp_tot (:,:) = emp_oce
      qns_tot (:,:) = (1._wp - at_i_b(:,:)) * qns_oce(:,:) + SUM( a_i_b(:,:,:) * qns_ice(:,:,:), dim=3 )
      qsr_tot (:,:) = (1._wp - at_i_b(:,:)) * qsr_oce(:,:) + SUM( a_i_b(:,:,:) * qsr_ice(:,:,:), dim=3 )

      ! --- shortwave radiation transmitted below the surface (W/m2, see Grenfell Maykut 77) --- !
      zfr1 = ( 0.18 * ( 1.0 - cldf_ice ) + 0.35 * cldf_ice )            ! transmission when hi>10cm
      zfr2 = ( 0.82 * ( 1.0 - cldf_ice ) + 0.65 * cldf_ice )            ! zfr2 such that zfr1 + zfr2 to equal 1
      !
      WHERE    ( phs(:,:,:) <= 0._wp .AND. phi(:,:,:) <  0.1_wp )       ! linear decrease from hi=0 to 10cm  
         qtr_ice_top(:,:,:) = qsr_ice(:,:,:) * ( zfr1 + zfr2 * ( 1._wp - phi(:,:,:) * 10._wp ) )
      ELSEWHERE( phs(:,:,:) <= 0._wp .AND. phi(:,:,:) >= 0.1_wp )       ! constant (zfr1) when hi>10cm
         qtr_ice_top(:,:,:) = qsr_ice(:,:,:) * zfr1
      ELSEWHERE                                                         ! zero when hs>0
         qtr_ice_top(:,:,:) = 0._wp 
      END WHERE
   END SUBROUTINE usrdef_sbc_ice_flx

   !!======================================================================
END MODULE usrdef_sbc
