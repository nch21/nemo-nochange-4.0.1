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
      ze1 = 106000._wp / REAL( nn_GYRE , wp )   ! gridspacing in meters
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
            !glamt(i,j) longitude at T-point
            !gphit(i,j) latitude at T-point  
            plamt(ji,jj) = zlam0 + zim05 * ze1deg * zcos_alpha + zjm05 * ze1deg * zsin_alpha
            pphit(ji,jj) = zphi0 - zim05 * ze1deg * zsin_alpha + zjm05 * ze1deg * zcos_alpha
            !   
            !glamu(i,j) longitude at U-point
            !gphiu(i,j) latitude at U-point
            plamu(ji,jj) = zlam0 + zim1  * ze1deg * zcos_alpha + zjm05 * ze1deg * zsin_alpha
            pphiu(ji,jj) = zphi0 - zim1  * ze1deg * zsin_alpha + zjm05 * ze1deg * zcos_alpha
            !   
            !glamv(i,j) longitude at V-point
            !gphiv(i,j) latitude at V-point
            plamv(ji,jj) = zlam0 + zim05 * ze1deg * zcos_alpha + zjm1  * ze1deg * zsin_alpha
            pphiv(ji,jj) = zphi0 - zim05 * ze1deg * zsin_alpha + zjm1  * ze1deg * zcos_alpha
            !
            !glamf(i,j) longitude at F-point
            !gphif(i,j) latitude at F-point 
            plamf(ji,jj) = zlam0 + zim1  * ze1deg * zcos_alpha + zjm1  * ze1deg * zsin_alpha
            pphif(ji,jj) = zphi0 - zim1  * ze1deg * zsin_alpha + zjm1  * ze1deg * zcos_alpha
         END DO
      END DO
      !
      !                       !== Horizontal scale factors ==! (in meters)
      !                     
      !                                         ! constant grid spacing
      pe1t(:,:) =  ze1     ;      pe2t(:,:) = ze1
      pe1u(:,:) =  ze1     ;      pe2u(:,:) = ze1
      pe1v(:,:) =  ze1     ;      pe2v(:,:) = ze1
      pe1f(:,:) =  ze1     ;      pe2f(:,:) = ze1
      !
      !                                         ! NO reduction of grid size in some straits 
      ke1e2u_v = 0                              !    ==>> u_ & v_surfaces will be computed in dom_ghr routine
      pe1e2u(:,:) = 0._wp                       !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp                       !             require an initialization of INTENT(out) arguments
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
      pff_f(:,:) = ( zf0 + zbeta * ABS( pphif(:,:) - zphi0 ) * rad * ra ) ! f = f0 +beta* y ( y=0 at south)
      pff_t(:,:) = ( zf0 + zbeta * ABS( pphit(:,:) - zphi0 ) * rad * ra ) ! f = f0 +beta* y ( y=0 at south)
      !
      IF(lwp) WRITE(numout,*) '                           beta-plane used. beta = ', zbeta, ' 1/(s.m)'
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
