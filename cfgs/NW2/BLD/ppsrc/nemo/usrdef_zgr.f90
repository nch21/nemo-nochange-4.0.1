/* Copyright (C) 1991-2023 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */


/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */

/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */



/* wchar_t uses Unicode 10.0.0.  Version 10.0 of the Unicode Standard is
   synchronized with ISO/IEC 10646:2017, fifth edition, plus
   the following additions from Amendment 1 to the fifth edition:
   - 56 emoji characters
   - 285 hentaigana
   - 3 additional Zanabazar Square characters */

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
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_zgr        ! called by domzgr.F90
   PUBLIC   usr_def_zgr_scalar
   
   !! * Substitutions
   !!----------------------------------------------------------------------
   !!                   ***  vectopt_loop_substitute  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute the inner loop start/end indices with CPP macro
   !!                allow unrolling of do-loop (useful with vector processors)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: vectopt_loop_substitute.h90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
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
      CALL zgr_z( jpi, jpj, jpk, jpim1, jpjm1, jpkm1, pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! Reference z-coordinate system
      !
      CALL zgr_msk_top_bot( jpi,jpj, jpkm1, k_top , k_bot )                 ! masked top and bottom ocean t-level indices
      !
      !                                                     ! z-coordinate (3D arrays) from the 1D z-coord.
      CALL zgr_zco( jpi, jpj, jpk, pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d,   &   ! in  : 1D reference vertical coordinate
         &          pdept   , pdepw   ,                     &   ! out : 3D t & w-points depth
         &          pe3t    , pe3u    , pe3v   , pe3f   ,   &   !       vertical scale factors
         &          pe3w    , pe3uw   , pe3vw             )     !           -      -      -
      !
   END SUBROUTINE usr_def_zgr

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
     ! wxt 50 levels
    !zsur =   -1456.379509196609_wp
    !za0  =    82.940596004514603_wp
    !za1  =    76.029887827318291_wp
    !zkth =    25.927805586296863_wp
    !zacr =    10.0000000000000000_wp
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
         zw = REAL( jk , wp )
         zt = REAL( jk , wp ) + 0.5_wp
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


   SUBROUTINE zgr_msk_top_bot( jpi,jpj, jpkm1, k_top , k_bot )
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

	INTEGER :: ji,jj
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '    zgr_top_bot : defines the top and bottom wet ocean levels.'
      IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) '       GYRE case : closed flat box ocean without ocean cavities'
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
         IF (ji+nimpp-1 <= 1 .or. ji+nimpp-1 >= jpiglo .or. jj+njmpp-1 <= 1 .or. jj+njmpp-1 >= jpjglo) THEN
                 z2d(ji,jj) = 0._wp
         ENDIF
         END DO
      END DO

      !
      !CALL lbc_lnk( 'usrdef_zgr', z2d, 'T', 1. )           ! set surrounding land to zero (here jperio=0 ==>> closed)
      !
      k_bot(:,:) = NINT( z2d(:,:) )           ! =jpkm1 over the ocean point, =0 elsewhere

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

   !!======================================================================
END MODULE usrdef_zgr
