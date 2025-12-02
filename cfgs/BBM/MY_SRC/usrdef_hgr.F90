MODULE usrdef_hgr
   !!======================================================================
   !!                     ***  MODULE usrdef_hgr   ***
   !!
   !!                     ===  GYRE configuration  ===
   !!
   !! User defined :   mesh and Coriolis parameter of a user configuration
   !!======================================================================
   !! History :  4.0 ! 2016-03  (S. Flavoni) 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_hgr   : initialize the horizontal mesh 
   !!----------------------------------------------------------------------
   USE dom_oce  , ONLY: nimpp, njmpp       ! ocean space and time domain
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE usrdef_nam     !
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_hgr   ! called in domhgr.F90
   PUBLIC   usr_def_hgr_scalar
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_hgr.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_hgr( plamt , plamu , plamv  , plamf  ,   &   ! geographic position (required)
      &                    pphit , pphiu , pphiv  , pphif  ,   &   !
      &                    kff   , pff_f , pff_t  ,            &   ! Coriolis parameter  (if domain not on the sphere)
      &                    pe1t  , pe1u  , pe1v   , pe1f   ,   &   ! scale factors       (required)
      &                    pe2t  , pe2u  , pe2v   , pe2f   ,   &   !
      &                    ke1e2u_v      , pe1e2u , pe1e2v     )   ! u- & v-surfaces (if gridsize reduction is used in strait(s))
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE usr_def_hgr  ***
      !!
      !! ** Purpose :   user defined mesh and Coriolis parameter
      !!
      !! ** Method  :   set all intent(out) argument to a proper value
      !!
      !!                Here GYRE configuration : 
      !!          Rectangular mid-latitude domain  
      !!          - with axes rotated by 45 degrees
      !!          - a constant horizontal resolution of 106 km 
      !!          - on a beta-plane
      !!
      !! ** Action  : - define longitude & latitude of t-, u-, v- and f-points (in degrees) 
      !!              - define coriolis parameter at f-point if the domain in not on the sphere (on beta-plane)
      !!              - define i- & j-scale factors at t-, u-, v- and f-points (in meters)
      !!              - define u- & v-surfaces (if gridsize reduction is used in some straits) (in m2)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   plamt, plamu, plamv, plamf   ! longitude outputs                     [degrees]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pphit, pphiu, pphiv, pphif   ! latitude outputs                      [degrees]
      INTEGER                 , INTENT(out) ::   kff                          ! =1 Coriolis parameter computed here, =0 otherwise
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pff_f, pff_t                 ! Coriolis factor at f-point                [1/s]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1t, pe1u, pe1v, pe1f       ! i-scale factors                             [m]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe2t, pe2u, pe2v, pe2f       ! j-scale factors                             [m]
      INTEGER                 , INTENT(out) ::   ke1e2u_v                     ! =1 u- & v-surfaces computed here, =0 otherwise 
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1e2u, pe1e2v               ! u- & v-surfaces (if reduction in strait)   [m2]
      !!external
	!mod par_oce

      !
      !INTEGER, INTENT(in)  ::   jpi, jpj, jpk
      !
      INTEGER  ::   ji, jj, jjc,jjr               ! dummy loop indices
      INTEGER  ::   jpi_nw, jpj_nw, jpjs_nw, AMR_ratio=1
      REAL(wp) ::   zlam1, zlam0, zdisti, zim1 , zjm1 , ze1deg, zf0 ! local scalars
      REAL(wp) ::   zphi1, zphi0, zdistj, zim05, zjm05, ze2deg, zbeta, znorme, zdj      !   -      -  
      !!-------------------------------------------------------------------------------
      !
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : NW2 configuration (geographical mesh on the sphere)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
       !!-------------------------------------------------------------------------------
      !
      !     !== geographical mesh on the sphere with regular (in degree) grid-spacing in longitude  ==! 
      !
      !                       !==  grid point position  ==!
      !
      zlam0 = 0.                           ! position of southwest (left bottom) corner, longitude
      zphi0 = -75.                         ! latitude
      zlam1 = 60.                           ! position of noutheast (right top) corner, longitude
      zphi1 = 70.                         ! latitude
      !  
      !
      ! number of grid points in the i and j directions
      jpi_nw = ceiling(( zlam1 - zlam0 ) * nn_GYRE * AMR_ratio )
      ! equatorial 1/nn_GYRE degree but 1/(2*nn_GYRE) degree at poles
      ! I want it decrease linearly from equator to poles

      ze2deg = 0.5 / REAL( nn_GYRE*AMR_ratio  )   ! gridspacing in degrees at poles
      ! 0 = -75 + [ 1/2n + (1/2n + D) + (1/2n + 2D) + ... + (1/2n + (N-1)D) ]
      ! 1/2n + (N-1)D = 1/n

      ! D = 1/(2n)/(N-1)
      ! 0 = zphi0 + N*[1/n + (N-1)*1/(2n)/(N-1)]/2
      ! 0 = zphi0 + N*[1/n + 1/(2n)]/2
      ! -zphi0 = N*[3/(2n)]/2
      ! N = -zphi0*4n/3
      ! D = 1/(2n)/(-zphi0*4n/3-1)
      ! zphi1 = [ 1/n + (1/n - D) + (1/n - 2D) + ... + (1/n - (N2-1)D) ]
      ! zphi1 = N2*[2/n - (N2-1)D]/2
      !N2 = 2/3 n (4 phi0 n + 3) ((-(4 phi0)/(4 phi0 n + 3) - 9/(4 n (4 phi0 n + 3))) + sqrt((3 phi1)/(n (4 phi0 n + 3)) + ((4 phi0)/(4 phi0 n + 3) + 9/(4 n (4 phi0 n + 3)))^2)) 

      ! for n = 2,
      ! N = 75*4/3*2 = 200
      ! D = 1/796
      ! N2 = 797/2 - (13*sqrt(1121))/2 = 180.9
      jpjs_nw = ceiling(-zphi0*4./3.*nn_GYRE) 
      jpj_nw = jpjs_nw + ceiling (2./3.*nn_GYRE*(4.*zphi0*nn_GYRE+3.)* ((-4.*zphi0)/(4.*zphi0*nn_GYRE+3.)- 9./(4.*nn_GYRE*(4.*zphi0*nn_GYRE+3.)) + sqrt((3.*zphi1)/(nn_GYRE*(4.*zphi0*nn_GYRE+3.)) +  ((4.*zphi0)/(4.*zphi0*nn_GYRE+3.) + 9./(4.*nn_GYRE*(4.*zphi0*nn_GYRE+3.)))**2.))) 

      zdj = 1./2./nn_GYRE/(ceiling(-zphi0*4./3.*nn_GYRE)-1.)

      ze1deg = 1.0 / REAL( nn_GYRE ) / AMR_ratio  ! gridspacing in degrees
      ! for 0.5 degree, jasmin jpjglo = 381, jpiglo = 120
      !  
      DO jj = 1, jpj 
         DO ji = 1, jpi
            zim1 = REAL( ji + nimpp - 2 ) - 1.   ;   zim05 = REAL( ji + nimpp - 2 ) - 1.5 
            zjm1 = REAL( jj + njmpp - 2 ) - 1.   ;   zjm05 = REAL( jj + njmpp - 2 ) - 1.5 

            jjc = floor((zjm1)/AMR_ratio)+1
            jjr = zjm1+1 - (jjc-1)*AMR_ratio
 
            if ((jjc <= jpjs_nw) ) then
               !  -75 + [ 1/2n + (1/2n + D) + (1/2n + 2D) + ... + (1/2n + (N-1)D) ]
               pphif(ji,jj) = zphi0 + (1./nn_GYRE + ((jjc) - 1.)*zdj)*(jjc)/2.

              ! if (ji==2) write(*,*) 'pphif jj',pphif(ji,jj),'jjc=',jjc
               pphif(ji,jj) = pphif(ji,jj) - (0.5/nn_GYRE + ((jjc-1.)*zdj))/AMR_ratio*(AMR_ratio-jjr)
               !if (ji==2) write(*,*) 'pphif jj',pphif(ji,jj),'jjr=',jjr,'zjmi=',zjm1
               pphit(ji,jj) = pphif(ji,jj) - (0.5/nn_GYRE + ((jjc-1.)*zdj))/AMR_ratio/2.
               pphiv(ji,jj) = pphif(ji,jj)
               pphiu(ji,jj) = pphit(ji,jj)
            else 
               ! 1/n + (1/n - D) + (1/n - 2D) + ... + (1/n - (N2-1)D), N2 = zj - jpjs_nw
               !  N2*[2/n - (N2-1)D]/2
               pphif(ji,jj) = (jjc - jpjs_nw)*(2./nn_GYRE - (jjc - 1. - jpjs_nw)*zdj)/2.
              ! if (ji==2) write(*,*) 'pphif jj',pphif(ji,jj),'jjc=',jjc
               pphif(ji,jj) = pphif(ji,jj) - (1./nn_GYRE - (jjc -1. - jpjs_nw)*zdj)/AMR_ratio*(AMR_ratio-jjr)
               ! if (ji==2) write(*,*) 'pphif jj',pphif(ji,jj),'jjr=',jjr,'zjmi=',zjm1
               pphit(ji,jj) = pphif(ji,jj) - (1./nn_GYRE - (jjc -1. - jpjs_nw)*zdj)/AMR_ratio/2.
               pphiv(ji,jj) = pphif(ji,jj)
               pphiu(ji,jj) = pphit(ji,jj)
            end if
            !   
            plamt(ji,jj) = zlam0 + (zim05+1.) * ze1deg 
            !   
            !glamu(i,j) longitude at U-point
            !gphiu(i,j) latitude at U-point
            plamu(ji,jj) = zlam0 + (zim1+1.)  * ze1deg 
            !   
            !glamv(i,j) longitude at V-point
            !gphiv(i,j) latitude at V-point
            plamv(ji,jj) = zlam0 + (zim05+1.) * ze1deg 
            !
            !glamf(i,j) longitude at F-point
            !gphif(i,j) latitude at F-point 
            plamf(ji,jj) = zlam0 + (zim1+1.)  * ze1deg 
         END DO
      END DO
      !
      !                       !== Horizontal scale factors ==! (in meters)
      !                 
      DO jj = 1, jpj 
         pe1t(:,jj) = ra * rad * COS( rad * pphit(1,jj) ) * ze1deg
         pe1u(:,jj) = ra * rad * COS( rad * pphiu(1,jj) ) * ze1deg
         pe1v(:,jj) = ra * rad * COS( rad * pphiv(1,jj) ) * ze1deg
         pe1f(:,jj) = ra * rad * COS( rad * pphif(1,jj) ) * ze1deg
         if (jj==1) then
            pe2t(:,jj) = ra * rad * (pphif(1,jj+1) - pphif(1,jj)) 
            pe2u(:,jj) = ra * rad * (pphif(1,jj+1) - pphif(1,jj)) 
            pe2v(:,jj) = ra * rad * (pphit(1,jj+1) - pphit(1,jj))
            pe2f(:,jj) = ra * rad * (pphit(1,jj+1) - pphit(1,jj))
         ELSE IF (jj==jpj) THEN
            pe2t(:,jj) = ra * rad * (pphif(1,jj) - pphif(1,jj-1)) 
            pe2u(:,jj) = ra * rad * (pphif(1,jj) - pphif(1,jj-1))
            pe2v(:,jj) = ra * rad * (pphit(1,jj) - pphit(1,jj-1))
            pe2f(:,jj) = ra * rad * (pphit(1,jj) - pphit(1,jj-1))
         ELSE
            pe2t(:,jj) = ra * rad * (pphif(1,jj) - pphif(1,jj-1)) 
            pe2u(:,jj) = ra * rad * (pphif(1,jj) - pphif(1,jj-1))
            pe2v(:,jj) = ra * rad * (pphit(1,jj+1) - pphit(1,jj))
            pe2f(:,jj) = ra * rad * (pphit(1,jj+1) - pphit(1,jj))
         end if
      END DO
      !
      !                                         ! NO reduction of grid size in some straits 
      ke1e2u_v = 0                              !    ==>> u_ & v_surfaces will be computed in dom_ghr routine
      pe1e2u( :,:) = 0.0                       !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v( :,:) = 0.0                       !             require an initialization of INTENT(out) arguments
      !
      !
      !                       !==  Coriolis parameter  ==!
      kff = 1                                            !  indicate not to compute ff afterward
      !
      !
      DO jj = 1, jpj 
         DO ji = 1, jpi 
         	pff_f(ji,jj) = 2. * omega * SIN( rad * pphif(ji,jj) )
         	pff_t(ji,jj) = 2. * omega * SIN( rad * pphit(ji,jj) )
         END DO
      END DO


      !
   END SUBROUTINE usr_def_hgr

   SUBROUTINE usr_def_hgr_scalar( kff   , ke1e2u_v )   ! u- & v-surfaces (if gridsize reduction is used in strait(s))
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE usr_def_hgr  ***
      !!
      !! ** Purpose :   user defined mesh and Coriolis parameter
      !!
      !! ** Method  :   set all intent(out) argument to a proper value
      !!
      !!                Here GYRE configuration :
      !!          Rectangular mid-latitude domain 
      !!          - with axes rotated by 45 degrees
      !!          - a constant horizontal resolution of 106 km 
      !!          - on a beta-plane
      !!
      !! ** Action  : - define longitude & latitude of t-, u-, v- and f-points (in degrees) 
      !!              - define coriolis parameter at f-point if the domain in not on the sphere (on beta-plane)
      !!              - define i- & j-scale factors at t-, u-, v- and f-points (in meters)
      !!              - define u- & v-surfaces (if gridsize reduction is used in some straits) (in m2)
      !!----------------------------------------------------------------------
      INTEGER                 , INTENT(out) ::   kff                          ! =1 Coriolis parameter computed here, =0 otherwise
      INTEGER                 , INTENT(out) ::   ke1e2u_v                     ! =1 u- & v-surfaces computed here, =0 otherwise 
      !
      INTEGER  ::   ji, jj               ! dummy loop indices
      INTEGER  ::   jpi_tmp, jpj_tmp
      REAL(wp) ::   zlam1, zlam0, zcos_alpha, zim1 , zjm1 , ze1  , ze1deg, zf0 ! local scalars
      REAL(wp) ::   zphi1, zphi0, zsin_alpha, zim05, zjm05, zbeta, znorme      !   -      -
      !!-------------------------------------------------------------------------------
      !
      !     !==  beta-plane with regular grid-spacing and rotated domain ==!  (GYRE configuration)
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : GYRE configuration (beta-plane with rotated regular grid-spacing)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      !
      !                       !==  grid point position  ==!
      !
      zlam1 = -85._wp                           ! position of gridpoint (i,j) = (1,jpjglo)
      zphi1 =  29._wp
      !
      ze1 = 106000._wp / REAL( nn_GYRE*AMR_ratio , wp )   ! gridspacing in meters
      !
      !zsin_alpha = - SQRT( 2._wp ) * 0.5_wp     ! angle: 45 degrees
      !zcos_alpha =   SQRT( 2._wp ) * 0.5_wp
      !
      zsin_alpha =   0.0_wp     ! angle: 0 degrees
      zcos_alpha =   1.0_wp
      
      ze1deg = ze1 / (ra * rad)
      zlam0 = zlam1 + zcos_alpha * ze1deg * REAL( jpjglo-2 , wp )
      zphi0 = zphi1 + zsin_alpha * ze1deg * REAL( jpjglo-2 , wp )
      !zlam0 = zlam1 + zcos_alpha * ze1deg * REAL( jpjglo-2 , wp ) * AMR_ratio
      !zphi0 = zphi1 + zsin_alpha * ze1deg * REAL( jpjglo-2 , wp ) * AMR_ratio
#if defined key_agrif
      ! ! Upper left longitude and latitude from parent:
      IF (.NOT.Agrif_root()) THEN
         zlam0 = zlam1 + Agrif_irhox() * REAL(Agrif_Parent(jpjglo)-2 , wp) * ze1deg * zcos_alpha  &
                   &   + ( Agrif_Ix()*Agrif_irhox()-(0.5_wp+nbghostcells)) * ze1deg * zcos_alpha  &
                   &   + ( Agrif_Iy()*Agrif_irhoy()-(0.5_wp+nbghostcells)) * ze1deg * zsin_alpha
         zphi0 = zphi1 + Agrif_irhoy() * REAL(Agrif_Parent(jpjglo)-2 , wp) * ze1deg * zsin_alpha  &
                   &   - ( Agrif_Ix()*Agrif_irhox()-nbghostcells )         * ze1deg * zsin_alpha  &
                   &   + ( Agrif_Iy()*Agrif_irhoy()-nbghostcells )         * ze1deg * zcos_alpha
      ENDIF 
#endif
      !   
      IF( ln_bench ) THEN     ! benchmark: forced the resolution to be 106 km 
         ze1 = 106000._wp     ! but keep (lat,lon) at the right nn_GYRE resolution
         CALL ctl_warn( ' GYRE used as Benchmark: e1=e2=106km, no need to adjust rdt, ahm,aht ' )
      ENDIF
      !  
      DO jj = 1, jpj 
         DO ji = 1, jpi
            zim1 = REAL( ji + nimpp - 1 ) - 1.   ;   zim05 = REAL( ji + nimpp - 1 ) - 1.5 
            zjm1 = REAL( jj + njmpp - 1 ) - 1.   ;   zjm05 = REAL( jj + njmpp - 1 ) - 1.5 
            !   
         END DO
      END DO
      !
      !                       !== Horizontal scale factors ==! (in meters)
      !                     
      !
      !                                         ! NO reduction of grid size in some straits 
      ke1e2u_v = 0                              !    ==>> u_ & v_surfaces will be computed in dom_ghr routine
      !
      !
      !                       !==  Coriolis parameter  ==!
      kff = 1                                            !  indicate not to compute ff afterward
      !
      zbeta = 2. * omega * COS( rad * zphi1 ) / ra       ! beta at latitude zphi1
      !SF we overwrite zphi0 (south point in latitude) used just above to define pphif (value of zphi0=15.5190567531966)
      !SF for computation of Coriolis we keep the parameter of Hazeleger, W., and S. S. Drijfhout, JPO 1998.
      zphi0 = 15._wp                                     !  latitude of the most southern grid point  
      zf0   = 2. * omega * SIN( rad * zphi0 )            !  compute f0 1st point south
      !
      !
      IF(lwp) WRITE(numout,*) '                           beta-plane used. beta = ', zbeta, ' 1/(s.m)'
      !
   END SUBROUTINE usr_def_hgr_scalar


   !!======================================================================
END MODULE usrdef_hgr
