MODULE usrdef_istate
   !!======================================================================
   !!                   ***  MODULE  usrdef_istate   ***
   !!
   !!                     ===  GYRE configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  4.0 ! 2016-03  (S. Flavoni) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   ! nch : for NW2
   USE dom_oce        ! ocean space and time domain
   USE cubic_interp   ! cubic interpolation
   ! end
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate   ! called in istate.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here GYRE configuration example : (double gyre with rotated domain)
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity   field
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      
      !
      INTEGER :: ji, jj, jk  ! dummy loop indices
      ! ! nch : for NW2
      integer, parameter :: sst_nsize =4, pts_nsize =5
      real(wp), dimension(sst_nsize)  :: ysst   ! latitude of interpolation
      real(wp), dimension(sst_nsize)  :: sstd   ! sst of interpolation
      real(wp), dimension(pts_nsize)  :: ypts   ! latitude of interpolation
      real(wp), dimension(pts_nsize)  :: ptsd   ! emp of interpolation
      real(wp) :: t_star, latw
      !data ysst /-75.-2.803,0.-2.802,70.-2.802/
      sstd = (/-2.9328,-3.4,28.3,-30./)
      ysst = (/-77.803,-70.803,-2.803,137.197/)
      ypts = (/ -75., -25., 0., 25., 70./)
      ptsd = (/ 33.6,36.6,35.6,36.6,34./)
      
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : NW2 analytical definition of initial state todo' 
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   Ocean at rest, with an horizontally uniform T and S profiles'
      !
      pu  (:,:,:) = 0._wp        ! ocean at rest
      pv  (:,:,:) = 0._wp
      pssh(:,:)   = 0._wp 
      !
      DO jk = 1, jpk             ! horizontally uniform T & S profiles
         DO jj = 1, jpj
            DO ji = 1, jpi
               t_star = itau(ysst,sstd,gphit(ji,jj))
               latw = exp(-(((gphit(ji,jj) + 20.))**2) / 550.) + exp(-(((gphit(ji,jj) - 20.))**2) / 550.) -  exp(-(((gphit(ji,jj) + 80.)**1.5)) / 400.)/3.5
              ! t_star = 10.
               pts(ji,jj,jk,jp_tem) =  (  (  16. - 12. * TANH( (pdept(ji,jj,jk) - 400) / 700 ) )   &
                    &           * (-TANH( (500. - pdept(ji,jj,jk)) / 150. ) + 1.) / 2. * latw          &
                    !&           * (-TANH( (500. - pdept(ji,jj,jk)) / 150. ) + 1.) / 2. * (1.-(gphit(ji,jj)/75.)**2)*0.75          &
                    &           + ( 15. * ( 1. - TANH( (pdept(ji,jj,jk)-50.) / 1500.) )            &
                    &           - 1.4 * TANH((pdept(ji,jj,jk)-100.) / 100.)                        &
                    &           + 7.  * (1500. - pdept(ji,jj,jk) ) / 1500.) * t_star   /   23.5                   &
                    &           * (-TANH( (pdept(ji,jj,jk) - 500.) / 150.) + 1.) / 2.  ) * ptmask(ji,jj,jk)

               t_star = 0.2

               pts(ji,jj,jk,jp_sal) = ( 34.8                  &
                    &         + ( 0.55 + 1.25 * (5000. - pdept(ji,jj,jk)) / 5000.                 &
                    &         - 1.62 * TANH( (pdept(ji,jj,jk) - 60.  ) / 650. )                    &
                    &         + 0.2  * TANH( (pdept(ji,jj,jk) - 35.  ) / 100. )                    &
                    &         + 0.2  * TANH( (pdept(ji,jj,jk) - 1000.) / 5000.) ) *t_star /1.84                 &
                    &         * (-TANH( (pdept(ji,jj,jk) - 500.) / 150.) + 1.) / 2  ) * ptmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !   
   END SUBROUTINE usr_def_istate

   !!======================================================================
END MODULE usrdef_istate
