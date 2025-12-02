MODULE usrdef_zgr
   !!======================================================================
   !!                       ***  MODULE  usrdef_zgr  ***
   !!
   !!                       ===  GYRE configuration  ===
   !!
   !! User defined : vertical coordinate system of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-06  (G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_zgr   : user defined vertical coordinate system
   !!      zgr_z      : reference 1D z-coordinate 
   !!      zgr_top_bot: ocean top and bottom level indices
   !!      zgr_zco    : 3D verticl coordinate in pure z-coordinate case
   !!---------------------------------------------------------------------
   USE oce            ! ocean variables
   USE dom_oce        ! ocean domain
   USE depth_e3       ! depth <=> e3
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library
   USE usrdef_nam, ONLY: AMR_ratio
   use topo_builder
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_zgr        ! called by domzgr.F90
   PUBLIC   usr_def_zgr_scalar
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_zgr.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS             

   SUBROUTINE usr_def_zgr(  ld_zco  , ld_zps  , ld_sco  , ld_isfcav,    &   ! type of vertical coordinate
      &                    pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d  ,    &   ! 1D reference vertical coordinate
      &                    pdept , pdepw ,                             &   ! 3D t & w-points depth
      &                    pe3t  , pe3u  , pe3v   , pe3f ,             &   ! vertical scale factors
      &                    pe3w  , pe3uw , pe3vw         ,             &   !     -      -      -
      &                    k_top  , k_bot    )                             ! top & bottom ocean level
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE usr_def_zgr  ***
      !!
      !! ** Purpose :   User defined the vertical coordinates
      !!
      !!----------------------------------------------------------------------
     !mod par_oce

	!!self
      LOGICAL                   , INTENT(out) ::   ld_zco, ld_zps, ld_sco      ! vertical coordinate flags
      LOGICAL                   , INTENT(out) ::   ld_isfcav                   ! under iceshelf cavity flag
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d           ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pdept, pdepw                ! grid-point depth        [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3w , pe3uw, pe3vw         ! i-scale factors 
      INTEGER , DIMENSION(:,:)  , INTENT(out) ::   k_top, k_bot                ! first & last ocean level
      !
      INTEGER  ::   inum   ! local logical unit
      REAL(WP) ::   z_zco, z_zps, z_sco, z_cav
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj) ::   z_bot ! bottom depth
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : GYRE configuration (z-coordinate closed flat box ocean without cavities)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      !
      ! type of vertical coordinate
      ! ---------------------------
      ld_zps    = .TRUE.         ! GYRE case:  z-coordinate without ocean cavities
      ld_zco    = .FALSE.
      ld_sco    = .FALSE.
      ld_isfcav = .FALSE.
      !
      !
      ! Build the vertical coordinate system
      ! ------------------------------------
      CALL zgr_z( jpi, jpj, jpk, jpim1, jpjm1, jpkm1, pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! Reference z-coordinate system
      !
      CALL zgr_bot( z_bot ) 
      CALL zgr_msk_top_bot( jpi,jpj, jpkm1, k_top , k_bot, z_bot , pdept_1d  )                 ! masked top and bottom ocean t-level indices
      !
      
      if (ld_zco) then
      !                                                     ! z-coordinate (3D arrays) from the 1D z-coord.
      CALL zgr_zco( jpi, jpj, jpk, pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d,   &   ! in  : 1D reference vertical coordinate
         &          pdept   , pdepw   ,                     &   ! out : 3D t & w-points depth
         &          pe3t    , pe3u    , pe3v   , pe3f   ,   &   !       vertical scale factors
         &          pe3w    , pe3uw   , pe3vw             )     !           -      -      -
      endif

      if (ld_zps) then
      !                                                     ! z-partial step coordinate (3D arrays) 
      CALL zgr_zps( jpi, jpj, jpk,jpim1,jpjm1,jpkm1, pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d,   &   ! in  : 1D reference vertical coordinate
         &          pdept   , pdepw   ,                     &   ! out : 3D t & w-points depth
         &          pe3t    , pe3u    , pe3v   , pe3f   ,   &   !       vertical scale factors
         &          pe3w    , pe3uw   , pe3vw           ,   &
         &          z_bot   , k_bot  )     !           -      -      -
      endif
   END SUBROUTINE usr_def_zgr
   SUBROUTINE zgr_bot( bot )   ! bottom depth
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_bot  ***
      !!
      !! ** Purpose :   compute the bottom depth
      !!
      !! ** Method  :   bottom depth = 4000 m, with boundary 0 liner interpolation to 4000
      !!
      !! ** Action  :   z_bot(:,:)   = bottom depth
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::   bot   ! bottom depth
      REAL(wp),  DIMENSION(jpj,jpi) ::  z_bot
      INTEGER ::   ji, jj   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
    
    real ::  XC2(jpj,jpi), YC2(jpj,jpi)
    real, parameter :: cw = 5., cd = 200., drake = 2500. !cw=5
    real, parameter :: NW2_latS = -75., NW2_latN = 70, NW2_lonW = 0., NW2_lonE = 60.
    real, parameter ::  D0 = 4000.

    bot = -D0

    do ji=1,jpi
      do jj=1,jpj
        XC2(jj,ji) = glamt(ji,jj)
        YC2(jj,ji) = gphit(ji,jj)
        z_bot(jj,ji) = bot(ji,jj)
      enddo
      enddo

        !call add_NS_coast(NW2_lonW, -90., 90., cw, cd, z_bot, XC2, YC2, D0)
    !call add_NS_coast(NW2_lonE, -90., 90., cw, cd, z_bot, XC2, YC2, D0)

    call add_NS_coast(NW2_lonW, -40., 90., cw, cd, z_bot, XC2, YC2, D0)
    call add_NS_coast(NW2_lonE, -40., 90., cw, cd, z_bot, XC2, YC2, D0)
    !call add_NS_coast(NW2_lonW+2, -40., 90., cw, cd, z_bot, XC2, YC2, D0)
   !call add_NS_coast(NW2_lonE-2, -40., 90., cw, cd, z_bot, XC2, YC2, D0)
    call add_NS_coast(NW2_lonW, -90., -60., cw, cd, z_bot, XC2, YC2, D0)
    call add_NS_coast(NW2_lonE, -90., -60., cw, cd, z_bot, XC2, YC2, D0)
      !call add_NS_coast(NW2_lonW+2, -90., -60., cw, cd, z_bot, XC2, YC2, D0)
   !call add_NS_coast(NW2_lonE-2, -90., -60., cw, cd, z_bot, XC2, YC2, D0)
    call add_EW_coast(-360., 360., NW2_latS, cw, cd, z_bot, XC2, YC2, D0)
    call add_EW_coast(-360., 360., NW2_latN, cw, cd, z_bot, XC2, YC2, D0)
    call add_NS_ridge(30., -90., 90., 20., D0/2.,0.,1., z_bot, XC2, YC2, D0)
    call add_circular_ridge(NW2_lonW, -50., 10., 2., drake, 0., z_bot, XC2, YC2, D0)
    
   ! z_bot = -D0

   do ji=1,jpi
      do jj=1,jpj
        bot(ji,jj) = z_bot(jj,ji)
      enddo
      enddo
   END SUBROUTINE zgr_bot

   SUBROUTINE usr_def_zgr_scalar( ld_zco  , ld_zps  , ld_sco  , ld_isfcav,    &   ! type of vertical coordinate
             &                    pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d)                             ! top & bottom ocean level
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE usr_def_zgr  ***
      !!
      !! ** Purpose :   User defined the vertical coordinates
      !!
      !!----------------------------------------------------------------------
	   !!self
      LOGICAL                   , INTENT(out) ::   ld_zco, ld_zps, ld_sco      ! vertical coordinate flags
      LOGICAL                   , INTENT(out) ::   ld_isfcav                   ! under iceshelf cavity flag
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d           ! 1D grid-point depth     [m]
      !
      INTEGER  ::   inum   ! local logical unit
      REAL(WP) ::   z_zco, z_zps, z_sco, z_cav
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : GYRE configuration (z-coordinate closed flat box ocean without cavities)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      !
      ! type of vertical coordinate
      ! ---------------------------
      ld_zco    = .TRUE.         ! GYRE case:  z-coordinate without ocean cavities
      ld_zps    = .FALSE.
      ld_sco    = .FALSE.
      ld_isfcav = .FALSE.
      !
      !
      ! Build the vertical coordinate system
      ! ------------------------------------
      CALL zgr_z_scalar(pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! Reference z-coordinate system
      !
      !           -      -      -
      !
   END SUBROUTINE usr_def_zgr_scalar

   SUBROUTINE zgr_z( jpi, jpj, jpk, jpim1, jpjm1, jpkm1, pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! 1D reference vertical coordinate
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_z  ***
      !!
      !! ** Purpose :   set the 1D depth of model levels and the resulting 
      !!              vertical scale factors.
      !!
      !! ** Method  :   1D z-coordinate system (use in all type of coordinate)
      !!       The depth of model levels is set from dep(k), an analytical function:
      !!                   w-level: depw_1d  = dep(k)
      !!                   t-level: dept_1d  = dep(k+0.5)
      !!       The scale factors are the discrete derivative of the depth:
      !!                   e3w_1d(jk) = dk[ dept_1d ] 
      !!                   e3t_1d(jk) = dk[ depw_1d ]
      !!           with at top and bottom :
      !!                   e3w_1d( 1 ) = 2 * ( dept_1d( 1 ) - depw_1d( 1 ) )
      !!                   e3t_1d(jpk) = 2 * ( dept_1d(jpk) - depw_1d(jpk) )
      !!       The depth are then re-computed from the sum of e3. This ensures 
      !!    that depths are identical when reading domain configuration file. 
      !!    Indeed, only e3. are saved in this file, depth are compute by a call
      !!    to the e3_to_depth subroutine.
      !!
      !!       Here the Madec & Imbard (1996) function is used.
      !!
      !! ** Action  : - pdept_1d, pdepw_1d : depth of T- and W-point (m)
      !!              - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!
      !! Reference : Marti, Madec & Delecluse, 1992, JGR, 97, No8, 12,763-12,766.
      !!             Madec and Imbard, 1996, Clim. Dyn.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d   ! 1D grid-point depth        [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d    ! 1D vertical scale factors  [m]
      !!external
	!mod par_oce
      INTEGER, INTENT(in)  ::   jpi, jpj, jpk
      INTEGER, INTENT(in)  ::   jpim1, jpjm1, jpkm1
      !
      INTEGER  ::   jk       ! dummy loop indices
      REAL(wp) ::   zt, zw   ! local scalars
      REAL(wp) ::   zsur, za0, za1, zkth, zacr   ! Values for the Madec & Imbard (1996) function  
      !!----------------------------------------------------------------------
      !
      ! Set parameters of z(k) function
      ! -------------------------------
      zsur = -2033.194295283385_wp       
      za0  =   155.8325369664153_wp 
      za1  =   146.3615918601890_wp
      zkth =    17.28520372419791_wp   
      zacr =     5.0_wp       
    !  nch 50 levels
    zsur =  -4575.57339036346456850879_wp
    za0  =   160.59150942018106889009_wp
    za1  =   157.34244572826057151360_wp
    zkth =  39.39930405999999862843_wp
    zacr =    15.0000000000000000_wp
      !
      IF(lwp) THEN            ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) '    zgr_z   : Reference vertical z-coordinates '
         WRITE(numout,*) '    ~~~~~~~'
         WRITE(numout,*) '       GYRE case : MI96 function with the following coefficients :'
         WRITE(numout,*) '                 zsur = ', zsur
         WRITE(numout,*) '                 za0  = ', za0
         WRITE(numout,*) '                 za1  = ', za1
         WRITE(numout,*) '                 zkth = ', zkth
         WRITE(numout,*) '                 zacr = ', zacr
      ENDIF

      !
      ! 1D Reference z-coordinate    (using Madec & Imbard 1996 function)
      ! -------------------------
      !
      DO jk = 1, jpk          ! depth at T and W-points
         zw = REAL( jk , wp ) -1.0
         zt = REAL( jk , wp ) - 0.5_wp
         pdepw_1d(jk) = ( zsur + za0 * zw + za1 * zacr *  LOG( COSH( (zw-zkth) / zacr ) )  )
         pdept_1d(jk) = ( zsur + za0 * zt + za1 * zacr *  LOG( COSH( (zt-zkth) / zacr ) )  )
      END DO
      !
      !                       ! e3t and e3w from depth
      CALL depth_to_e3( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d ) 
      !
      !                       ! recompute depths from SUM(e3)  <== needed
      CALL e3_to_depth( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d ) 
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '              Reference 1D z-coordinate depth and scale factors:'
         WRITE(numout, "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')" )
         WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk )
      ENDIF
      !
   END SUBROUTINE zgr_z

   SUBROUTINE zgr_z_scalar(pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! 1D reference vertical coordinate
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_z  ***
      !!
      !! ** Purpose :   set the 1D depth of model levels and the resulting 
      !!              vertical scale factors.
      !!
      !! ** Method  :   1D z-coordinate system (use in all type of coordinate)
      !!       The depth of model levels is set from dep(k), an analytical function:
      !!                   w-level: depw_1d  = dep(k)
      !!                   t-level: dept_1d  = dep(k+0.5)
      !!       The scale factors are the discrete derivative of the depth:
      !!                   e3w_1d(jk) = dk[ dept_1d ] 
      !!                   e3t_1d(jk) = dk[ depw_1d ]
      !!           with at top and bottom :
      !!                   e3w_1d( 1 ) = 2 * ( dept_1d( 1 ) - depw_1d( 1 ) )
      !!                   e3t_1d(jpk) = 2 * ( dept_1d(jpk) - depw_1d(jpk) )
      !!       The depth are then re-computed from the sum of e3. This ensures 
      !!    that depths are identical when reading domain configuration file. 
      !!    Indeed, only e3. are saved in this file, depth are compute by a call
      !!    to the e3_to_depth subroutine.
      !!
      !!       Here the Madec & Imbard (1996) function is used.
      !!
      !! ** Action  : - pdept_1d, pdepw_1d : depth of T- and W-point (m)
      !!              - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!
      !! Reference : Marti, Madec & Delecluse, 1992, JGR, 97, No8, 12,763-12,766.
      !!             Madec and Imbard, 1996, Clim. Dyn.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d   ! 1D grid-point depth        [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d    ! 1D vertical scale factors  [m]
      !
      INTEGER  ::   jk       ! dummy loop indices
      REAL(wp) ::   zt, zw   ! local scalars
      REAL(wp) ::   zsur, za0, za1, zkth, zacr   ! Values for the Madec & Imbard (1996) function  
      !!----------------------------------------------------------------------
      !
      ! Set parameters of z(k) function
      ! -------------------------------
      zsur = -2033.194295283385_wp       
      za0  =   155.8325369664153_wp 
      za1  =   146.3615918601890_wp
      zkth =    17.28520372419791_wp   
      zacr =     5.0_wp   
      zsur =   -1456.379509196609_wp
    za0  =    82.940596004514603_wp
    za1  =    76.029887827318291_wp
    zkth =    25.927805586296863_wp
    zacr =    10.0000000000000000_wp    
      !
      IF(lwp) THEN            ! Parameter print
      ENDIF

      !
      ! 1D Reference z-coordinate    (using Madec & Imbard 1996 function)
      ! -------------------------
      !
      DO jk = 1, jpk          ! depth at T and W-points
         zw = REAL( jk , wp )
         zt = REAL( jk , wp ) + 0.5_wp
         pdepw_1d(jk) = ( zsur + za0 * zw + za1 * zacr *  LOG( COSH( (zw-zkth) / zacr ) )  )
         pdept_1d(jk) = ( zsur + za0 * zt + za1 * zacr *  LOG( COSH( (zt-zkth) / zacr ) )  )
      END DO
      !
      !
      IF(lwp) THEN                        ! control print
      ENDIF
      !
   END SUBROUTINE zgr_z_scalar


   SUBROUTINE zgr_msk_top_bot( jpi,jpj, jpkm1, k_top , k_bot, z_bot , pdept_1d )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE zgr_msk_top_bot  ***
      !!
      !! ** Purpose :   set the masked top and bottom ocean t-levels
      !!
      !! ** Method  :   GYRE case = closed flat box ocean without ocean cavities
      !!                   k_top = 1     except along north, south, east and west boundaries
      !!                   k_bot = jpk-1 except along north, south, east and west boundaries
      !!
      !! ** Action  : - k_top : first wet ocean level index
      !!              - k_bot : last  wet ocean level index
      !!----------------------------------------------------------------------
     !mod par_oce
      INTEGER, INTENT(in) :: jpi, jpj
      INTEGER, INTENT(in) :: jpkm1
	!!elf
      INTEGER , DIMENSION(jpi,jpj), INTENT(out) ::   k_top , k_bot   ! first & last wet ocean level
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D local workspace
!     REAL(wp), DIMENSION(jpiglo,jpjglo) :: z2d_global ! Rocky: 2D global workspace
REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   z_bot   ! bottom depth
      REAL(wp), DIMENSION(:)   , INTENT(in) ::   pdept_1d   ! 1D depth

	INTEGER :: ji,jj,jk
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '    zgr_top_bot : defines the top and bottom wet ocean levels.'
      IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) '       NW2 case : closed flat box ocean without ocean cavities'
      !

      z2d(:,:) = REAL( jpkm1 , wp )          ! flat bottom
 !     !!Rocky: initialize global workspace
 !     z2d_global(:,:) = REAL( jpkm1, wp)
 !     z2d_global(1,:) = 0.
 !     z2d_global(jpiglo,:) = 0.
 !     z2d_global(:,1) = 0.
 !     z2d_global(:,jpjglo) = 0.
 !     !!Rocky: distribute to get z2d in every processor
 !     write(*,*) "usrdef_zgr: ni/jmpp=" ,nimpp, njmpp
 !     z2d(:,:) = z2d_global(nimpp:nimpp+jpi-1,njmpp:njmpp+jpj-1)
     DO jj = 1, jpj
         DO ji = 1, jpi
            DO jk = 1, jpk
               IF( pdept_1d(jk) >= -z_bot(ji,jj) ) EXIT
            END DO
            
            k_bot(ji,jj) = jk-1
         END DO
      END DO
      DO jj = 1, jpj
         DO ji = 1, jpi
         IF (jj+njmpp-1 <= 1 .or. jj+njmpp-1 >= jpjglo) THEN
        ! IF (ji+nimpp-1 <= 1 .or. ji+nimpp-1 >= jpiglo .or. jj+njmpp-1 <= 1 .or. jj+njmpp-1 >= jpjglo) THEN
                 z2d(ji,jj) = 0._wp
         ENDIF
         END DO
      END DO

      !
      !CALL lbc_lnk( 'usrdef_zgr', z2d, 'T', 1. )           ! set surrounding land to zero (here jperio=0 ==>> closed)
      !
      !k_bot(:,:) = NINT( z2d(:,:) )           ! =jpkm1 over the ocean point, =0 elsewhere

      !
      k_top(:,:) = MIN( 1 , k_bot(:,:) )     ! = 1    over the ocean point, =0 elsewhere
      

      !
   END SUBROUTINE zgr_msk_top_bot
   

   SUBROUTINE zgr_zco( jpi, jpj, jpk, pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d,   &   ! in : 1D reference vertical coordinate
      &                pdept   , pdepw   ,                     &   ! out: 3D t & w-points depth
      &                pe3t    , pe3u    , pe3v   , pe3f   ,   &   !      vertical scale factors
      &                pe3w    , pe3uw   , pe3vw             )     !          -      -      -
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_zco  ***
      !!
      !! ** Purpose :   define the reference z-coordinate system
      !!
      !! ** Method  :   set 3D coord. arrays to reference 1D array 
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(in   ) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth       [m]
      REAL(wp), DIMENSION(:)    , INTENT(in   ) ::   pe3t_1d , pe3w_1d           ! 1D vertical scale factors [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! grid-point depth          [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors    [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         !    -       -      -
      !!external
	!mod par_oce
      INTEGER, INTENT(in)  ::   jpi, jpj, jpk
      !
      INTEGER  ::   jk
      !!----------------------------------------------------------------------
      !
      DO jk = 1, jpk
         pdept(:,:,jk) = pdept_1d(jk)
         pdepw(:,:,jk) = pdepw_1d(jk)
         pe3t (:,:,jk) = pe3t_1d (jk)
         pe3u (:,:,jk) = pe3t_1d (jk)
         pe3v (:,:,jk) = pe3t_1d (jk)
         pe3f (:,:,jk) = pe3t_1d (jk)
         pe3w (:,:,jk) = pe3w_1d (jk)
         pe3uw(:,:,jk) = pe3w_1d (jk)
         pe3vw(:,:,jk) = pe3w_1d (jk)
      END DO
      !
   END SUBROUTINE zgr_zco


   SUBROUTINE zgr_zps( jpi, jpj, jpk, jpim1,jpjm1,jpkm1, gdept_1d, gdepw_1d, e3t_1d, e3w_1d,   &   ! in  : 1D reference vertical coordinate
         &          gdept_0   , gdepw_0   ,                     &   ! out : 3D t & w-points depth
         &          e3t_0    , e3u_0    , e3v_0   , e3f_0   ,   &   !       vertical scale factors
         &          e3w_0    , e3uw_0   , e3vw_0           ,   &
         &          bathy, mbathy  )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_zps  ***
      !!                     
      !! ** Purpose :   the depth and vertical scale factor in partial step
      !!              reference z-coordinate case
      !!
      !! ** Method  :   Partial steps : computes the 3D vertical scale factors
      !!      of T-, U-, V-, W-, UW-, VW and F-points that are associated with
      !!      a partial step representation of bottom topography.
      !!
      !!        The reference depth of model levels is defined from an analytical
      !!      function the derivative of which gives the reference vertical
      !!      scale factors.
      !!        From  depth and scale factors reference, we compute there new value
      !!      with partial steps  on 3d arrays ( i, j, k ).
      !!
      !!              w-level: gdepw_0(i,j,k)  = gdep(k)
      !!                       e3w_0(i,j,k) = dk(gdep)(k)     = e3(i,j,k)
      !!              t-level: gdept_0(i,j,k)  = gdep(k+0.5)
      !!                       e3t_0(i,j,k) = dk(gdep)(k+0.5) = e3(i,j,k+0.5)
      !!
      !!        With the help of the bathymetric file ( bathymetry_depth_ORCA_R2.nc),
      !!      we find the mbathy index of the depth at each grid point.
      !!      This leads us to three cases:
      !!
      !!              - bathy = 0 => mbathy = 0
      !!              - 1 < mbathy < jpkm1    
      !!              - bathy > gdepw_0(jpk) => mbathy = jpkm1  
      !!
      !!        Then, for each case, we find the new depth at t- and w- levels
      !!      and the new vertical scale factors at t-, u-, v-, w-, uw-, vw- 
      !!      and f-points.
      !! 
      !!        This routine is given as an example, it must be modified
      !!      following the user s desiderata. nevertheless, the output as
      !!      well as the way to compute the model levels and scale factors
      !!      must be respected in order to insure second order accuracy
      !!      schemes.
      !!
      !!         c a u t i o n : gdept_1d, gdepw_1d and e3._1d are positives
      !!         - - - - - - -   gdept_0, gdepw_0 and e3. are positives
      !!      
      !!  Reference :   Pacanowsky & Gnanadesikan 1997, Mon. Wea. Rev., 126, 3248-3270.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(in   ) ::   gdept_1d, gdepw_1d          ! 1D grid-point depth       [m]
      REAL(wp), DIMENSION(:)    , INTENT(in   ) ::   e3t_1d , e3w_1d           ! 1D vertical scale factors [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   gdept_0, gdepw_0                ! grid-point depth          [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   e3t_0 , e3u_0 , e3v_0 , e3f_0   ! vertical scale factors    [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   e3w_0 , e3uw_0, e3vw_0         !    -       -      -
      REAL(wp), DIMENSION(:,:)    , INTENT(inout) ::   bathy   ! bottom depth
      INTEGER , DIMENSION(:,:)    , INTENT(inout) ::   mbathy   ! last wet ocean level
      INTEGER, INTENT(in)  ::   jpi, jpj, jpk, jpkm1, jpim1, jpjm1
      INTEGER  ::   ji, jj, jk       ! dummy loop indices
      INTEGER  ::   ik, it, ikb, ikt ! temporary integers
      REAL(wp) ::   ze3tp , ze3wp    ! Last ocean level thickness at T- and W-points
      REAL(wp) ::   zdepwp, zdepth   ! Ajusted ocean depth to avoid too small e3t
      REAL(wp) ::   zdiff            ! temporary scalar
      REAL(wp) ::   zmax             ! temporary scalar
      REAL(wp) ::   e3zps_min, e3zps_rat   ! minimum thickness of partial steps and ratio
      LOGICAL ::   ln_dept_mid, ln_e3_dep   ! logical flag for ln_e3_dep.AND.ln_dept_mid
   
      
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  zprt
      bathy(:,:) = -bathy(:,:)

      !!---------------------------------------------------------------------
      !
      e3zps_min = 25._wp
      e3zps_rat = 0.2_wp
      ln_e3_dep   = .true.    ! =T : e3=dk[depth] in discret sens.
   !                       !      ===>>> will become the only possibility in v4.0
   !                       ! =F : e3 analytical derivative of depth function
   !                       !      only there for backward compatibility test with v3.6
   !                       !
      !                       ! if ln_e3_dep = T 
      ln_dept_mid = .false.   ! =T : set T points in the middle of cells 
      !                       ! =F : e3 analytical depth function
      !                       !
      
      ALLOCATE( zprt(jpi,jpj,jpk) )
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '    zgr_zps : z-coordinate with partial steps'
      IF(lwp) WRITE(numout,*) '    ~~~~~~~ '
      IF(lwp) WRITE(numout,*) '              mbathy is recomputed : bathy_level file is NOT used'

      ! compute position of the ice shelf grounding line
      ! set bathy and isfdraft to 0 where grounded
      !IF ( ln_isfcav ) CALL zgr_isf_zspace

      ! bathymetry in level (from bathy_meter)
      ! ===================
      zmax = gdepw_1d(jpk) + e3t_1d(jpk)        ! maximum depth (i.e. the last ocean level thickness <= 2*e3t_1d(jpkm1) )
      bathy(:,:) = MIN( zmax ,  bathy(:,:) )    ! bounded value of bathy (min already set at the end of zgr_bat)
      WHERE( bathy(:,:) == 0._wp )   ;   mbathy(:,:) = 0       ! land  : set mbathy to 0
      ELSE WHERE                     ;   mbathy(:,:) = jpkm1   ! ocean : initialize mbathy to the max ocean level
      END WHERE

      ! Compute mbathy for ocean points (i.e. the number of ocean levels)
      ! find the number of ocean levels such that the last level thickness
      ! is larger than the minimum of e3zps_min and e3zps_rat * e3t_1d (where
      ! e3t_1d is the reference level thickness
      DO jk = jpkm1, 1, -1
         zdepth = gdepw_1d(jk) + MIN( e3zps_min, e3t_1d(jk)*e3zps_rat )
         WHERE( 0._wp < bathy(:,:) .AND. bathy(:,:) <= zdepth )   mbathy(:,:) = jk-1
      END DO

      ! Check compatibility between bathy and iceshelf draft
      ! insure at least 2 wet level on the vertical under an ice shelf
      ! compute misfdep and adjust isf draft if needed
      !IF ( ln_isfcav ) CALL zgr_isf_kspace

      ! Scale factors and depth at T- and W-points
      DO jk = 1, jpk                        ! intitialization to the reference z-coordinate
         gdept_0(:,:,jk) = gdept_1d(jk)
         gdepw_0(:,:,jk) = gdepw_1d(jk)
         e3t_0  (:,:,jk) = e3t_1d  (jk)
         e3w_0  (:,:,jk) = e3w_1d  (jk)
      END DO
      
      ! Scale factors and depth at T- and W-points
      DO jj = 1, jpj
         DO ji = 1, jpi
            ik = mbathy(ji,jj)
            IF( ik > 0 ) THEN               ! ocean point only
               ! max ocean level case
               IF( ik == jpkm1 ) THEN
                  zdepwp = bathy(ji,jj)
                  ze3tp  = bathy(ji,jj) - gdepw_1d(ik)
                  ze3wp = 0.5_wp * e3w_1d(ik) * ( 1._wp + ( ze3tp/e3t_1d(ik) ) )
                  e3t_0(ji,jj,ik  ) = ze3tp
                  e3t_0(ji,jj,ik+1) = ze3tp
                  IF ( ln_e3_dep.AND.ln_dept_mid ) THEN
                     gdept_0(ji,jj,ik) = gdepw_1d(ik) + 0.5_wp * ze3tp
                     e3w_0(ji,jj,ik  ) = gdept_0(ji,jj,ik) - gdept_0(ji,jj,ik-1) 
                  ELSE
                     gdept_0(ji,jj,ik) = gdept_1d(ik-1) + ze3wp
                     e3w_0(ji,jj,ik  ) = ze3wp
                  ENDIF
                  gdept_0(ji,jj,ik+1) = gdept_0(ji,jj,ik) + ze3tp
                  e3w_0(ji,jj,ik+1) = ze3tp
                  gdepw_0(ji,jj,ik+1) = zdepwp
                  !
               ELSE                         ! standard case
                  IF( bathy(ji,jj) <= gdepw_1d(ik+1) ) THEN  ;   gdepw_0(ji,jj,ik+1) = bathy(ji,jj)
                  ELSE                                       ;   gdepw_0(ji,jj,ik+1) = gdepw_1d(ik+1)
                  ENDIF
!gm Bug?  check the gdepw_1d
                  !       ... on ik
                  e3t_0  (ji,jj,ik) = e3t_1d  (ik) * ( gdepw_0 (ji,jj,ik+1) - gdepw_1d(ik) )   & 
                     &                             / ( gdepw_1d(      ik+1) - gdepw_1d(ik) ) 
                  IF ( ln_e3_dep.AND.ln_dept_mid ) THEN
                     gdept_0(ji,jj,ik) = gdepw_1d(ik) + 0.5_wp * e3t_0(ji,jj,ik)
                     e3w_0(ji,jj,ik) = gdept_0(ji,jj,ik) - gdept_0(ji,jj,ik-1)
                  ELSE
                     gdept_0(ji,jj,ik) = gdepw_1d(ik) + ( gdepw_0(ji,jj,ik+1) - gdepw_1d(ik) )   &
                        &                             * ((gdept_1d(     ik  ) - gdepw_1d(ik) )   &
                        &                             / ( gdepw_1d(     ik+1) - gdepw_1d(ik) ))
                     e3w_0(ji,jj,ik) = 0.5_wp * ( gdepw_0(ji,jj,ik+1) + gdepw_1d(ik+1) - 2._wp * gdepw_1d(ik) )   &
                        &                     * ( e3w_1d(ik) / ( gdepw_1d(ik+1) - gdepw_1d(ik) ) )
                  ENDIF
                  !       ... on ik+1
                  e3w_0  (ji,jj,ik+1) = e3t_0  (ji,jj,ik)
                  e3t_0  (ji,jj,ik+1) = e3t_0  (ji,jj,ik)
                  gdept_0(ji,jj,ik+1) = gdept_0(ji,jj,ik) + e3t_0(ji,jj,ik)
               ENDIF
            ENDIF
         END DO
      END DO

      
      !
      !
      ! compute top scale factor if ice shelf
      !IF (ln_isfcav) CALL zps_isf
      !
      ! Scale factors and depth at U-, V-, UW and VW-points
      DO jk = 1, jpk                        ! initialisation to z-scale factors
         e3u_0 (:,:,jk) = e3t_1d(jk)
         e3v_0 (:,:,jk) = e3t_1d(jk)
         e3uw_0(:,:,jk) = e3w_1d(jk)
         e3vw_0(:,:,jk) = e3w_1d(jk)
      END DO

      DO jk = 1,jpk                         ! Computed as the minimum of neighbooring scale factors
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               e3u_0 (ji,jj,jk) = MIN( e3t_0(ji,jj,jk), e3t_0(ji+1,jj,jk) )
               e3v_0 (ji,jj,jk) = MIN( e3t_0(ji,jj,jk), e3t_0(ji,jj+1,jk) )
               e3uw_0(ji,jj,jk) = MIN( e3w_0(ji,jj,jk), e3w_0(ji+1,jj,jk) )
               e3vw_0(ji,jj,jk) = MIN( e3w_0(ji,jj,jk), e3w_0(ji,jj+1,jk) )
            END DO
         END DO
      END DO

      ! update e3uw in case only 2 cells in the water column
      !IF ( ln_isfcav ) CALL zps_isf_e3uv_w
      !
      CALL lbc_lnk('domzgr', e3u_0 , 'U', 1._wp )   ;   CALL lbc_lnk('domzgr', e3uw_0, 'U', 1._wp )   ! lateral boundary conditions
      CALL lbc_lnk('domzgr', e3v_0 , 'V', 1._wp )   ;   CALL lbc_lnk('domzgr', e3vw_0, 'V', 1._wp )
      !
      DO jk = 1, jpk                        ! set to z-scale factor if zero (i.e. along closed boundaries)
         WHERE( e3u_0 (:,:,jk) == 0._wp )   e3u_0 (:,:,jk) = e3t_1d(jk)
         WHERE( e3v_0 (:,:,jk) == 0._wp )   e3v_0 (:,:,jk) = e3t_1d(jk)
         WHERE( e3uw_0(:,:,jk) == 0._wp )   e3uw_0(:,:,jk) = e3w_1d(jk)
         WHERE( e3vw_0(:,:,jk) == 0._wp )   e3vw_0(:,:,jk) = e3w_1d(jk)
      END DO

      ! Scale factor at F-point
      DO jk = 1, jpk                        ! initialisation to z-scale factors
         e3f_0(:,:,jk) = e3t_1d(jk)
      END DO
      DO jk = 1, jpk                        ! Computed as the minimum of neighbooring V-scale factors
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               e3f_0(ji,jj,jk) = MIN( e3v_0(ji,jj,jk), e3v_0(ji+1,jj,jk) )
            END DO
         END DO
      END DO
      CALL lbc_lnk('domzgr', e3f_0, 'F', 1._wp )       ! Lateral boundary conditions
      !
      DO jk = 1, jpk                        ! set to z-scale factor if zero (i.e. along closed boundaries)
         WHERE( e3f_0(:,:,jk) == 0._wp )   e3f_0(:,:,jk) = e3t_1d(jk)
      END DO
!!gm  bug ? :  must be a do loop with mj0,mj1
      ! 

      e3t_0(:,mj0(1),:) = e3t_0(:,mj0(2),:)     ! we duplicate factor scales for jj = 1 and jj = 2
      e3w_0(:,mj0(1),:) = e3w_0(:,mj0(2),:) 
      e3u_0(:,mj0(1),:) = e3u_0(:,mj0(2),:) 
      e3v_0(:,mj0(1),:) = e3v_0(:,mj0(2),:) 
      e3f_0(:,mj0(1),:) = e3f_0(:,mj0(2),:) 

      ! Control of the sign
      IF( MINVAL( e3t_0  (:,:,:) ) <= 0._wp )   CALL ctl_stop( '    zgr_zps :   e r r o r   e3t_0 <= 0' )
      IF( MINVAL( e3w_0  (:,:,:) ) <= 0._wp )   CALL ctl_stop( '    zgr_zps :   e r r o r   e3w_0 <= 0' )
      IF( MINVAL( gdept_0(:,:,:) ) <  0._wp )   CALL ctl_stop( '    zgr_zps :   e r r o r   gdept_0 <  0' )
      IF( MINVAL( gdepw_0(:,:,:) ) <  0._wp )   CALL ctl_stop( '    zgr_zps :   e r r o r   gdepw_0 <  0' )
      !
      ! if in the future gde3w_0 need to be compute, use the function defined in NEMO
      ! for now gde3w_0 computation is removed as not an output of domcfg
     
      DEALLOCATE( zprt )
      bathy(:,:) = -bathy(:,:)

      !
   END SUBROUTINE zgr_zps

   subroutine add_NS_coast(lon, zlat0, zlat1, zdlon, shelf, z_bot, XC2, YC2, D0)
      real, intent(in) :: lon, zlat0, zlat1, zdlon, shelf, XC2(jpj,jpi), YC2(jpj,jpi), D0
      real :: r(jpj,jpi)
      real, intent(inout) :: z_bot(jpj,jpi)
      r = dist_from_line(XC2, lon, YC2, zlat0, zlat1)
      z_bot = max(z_bot, -D0 * coastal_sprofile(r, 0., zdlon, shelf/D0, .125, .125, .5) )
  end subroutine add_NS_coast

  subroutine add_EW_coast(zlon0, zlon1, lat, zdlat, shelf, z_bot, XC2, YC2, D0)
      real, intent(in) :: zlon0, zlon1, lat, zdlat, shelf, XC2(jpj,jpi), YC2(jpj,jpi), D0
      real :: r(jpj,jpi)
      real, intent(inout) :: z_bot(jpj,jpi)
      r = dist_from_line(YC2, lat, XC2, zlon0, zlon1)
      z_bot = max(z_bot, -D0 * coastal_sprofile(r, 0., zdlat, shelf/D0, .125, .125, .5) )
  end subroutine add_EW_coast

  subroutine add_circular_ridge(zlon0, zlat0, radius, dr, dH, clip, z_bot, XC2, YC2, D0)
      real, intent(in) :: zlon0, zlat0, radius, dr, dH, clip, XC2(jpj,jpi), YC2(jpj,jpi), D0
      real :: r(jpj,jpi)
      real, intent(inout) :: z_bot(jpj,jpi)
      r = sqrt( (XC2 - zlon0)**2 + (YC2 - zlat0)**2 )
      r = abs( r - radius)
      z_bot = max(z_bot, min(clip, D0 * ( clipped_cone(r, 0., dr, 1 - dH/D0) - 1 ) ) )
  end subroutine add_circular_ridge

  subroutine add_NS_ridge(lon, zlat0, zlat1, zdlon, dH, clip, p, z_bot, XC2, YC2, D0)
      real, intent(in) :: lon, zlat0, zlat1, zdlon, dH, clip, p, XC2(jpj,jpi), YC2(jpj,jpi), D0
      real :: r(jpj,jpi)
      real, intent(inout) :: z_bot(jpj,jpi)
      r = scurve(dist_from_line(XC2, lon, YC2, zlat0, zlat1), 0., zdlon)**p
      z_bot = max(z_bot, min(clip, (D0 - dH) * ( -r ) - dH))
  end subroutine add_NS_ridge
   !!======================================================================
END MODULE usrdef_zgr
