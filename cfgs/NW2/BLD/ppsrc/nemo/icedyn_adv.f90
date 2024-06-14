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

MODULE icedyn_adv
   !!======================================================================
   !!                       ***  MODULE icedyn_adv   ***
   !!   sea-ice: advection
   !!======================================================================
   !! History :  4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_dyn_adv   : advection of sea ice variables
   !!----------------------------------------------------------------------
   USE phycst         ! physical constant
   USE dom_oce        ! ocean domain
   USE sbc_oce , ONLY : nn_fsbc   ! frequency of sea-ice call
   USE ice            ! sea-ice: variables
   USE icevar         ! sea-ice: operations
   USE icedyn_adv_pra ! sea-ice: advection scheme (Prather)
   USE icedyn_adv_umx ! sea-ice: advection scheme (ultimate-macho)
   USE icectl         ! sea-ice: control prints
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing
   USE prtctl         ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_adv        ! called by icestp
   PUBLIC   ice_dyn_adv_init   ! called by icedyn

   INTEGER ::              nice_adv   ! choice of the type of advection scheme
   !                                        ! associated indices:
   INTEGER, PARAMETER ::   np_advPRA = 1   ! Prather scheme
   INTEGER, PARAMETER ::   np_advUMx = 2   ! Ultimate-Macho scheme
   !
   ! ** namelist (namdyn_adv) **
   INTEGER         ::   nn_UMx       ! order of the UMx advection scheme   
   !
   !! * Substitution
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
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icedyn_adv.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_adv( kt ) 
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_dyn_adv ***
      !!                    
      !! ** purpose : advection of sea ice
      !!
      !! ** method  : One can choose between 
      !!     a) an Ultimate-Macho scheme (with order defined by nn_UMx) => ln_adv_UMx
      !!     b) and a second order Prather scheme => ln_adv_Pra
      !!
      !! ** action :
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! number of iteration
      !!---------------------------------------------------------------------
      !
      ! controls
      IF( ln_timing    )   CALL timing_start('icedyn_adv')                                                             ! timing
      IF( ln_icediachk )   CALL ice_cons_hsm(0, 'icedyn_adv', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_adv: sea-ice advection'
         WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      !---------------!
      !== Advection ==!
      !---------------!
      SELECT CASE( nice_adv )
      !                                !-----------------------!
      CASE( np_advUMx )                ! ULTIMATE-MACHO scheme !
         !                             !-----------------------!
         CALL ice_dyn_adv_umx( nn_UMx, kt, u_ice, v_ice, h_i, h_s, h_ip, &
            &                          ato_i, v_i, v_s, sv_i, oa_i, a_i, a_ip, v_ip, e_s, e_i )
         !                             !-----------------------!
      CASE( np_advPRA )                ! PRATHER scheme        !
         !                             !-----------------------!
         CALL ice_dyn_adv_pra(         kt, u_ice, v_ice, &
            &                          ato_i, v_i, v_s, sv_i, oa_i, a_i, a_ip, v_ip, e_s, e_i )
      END SELECT

      !------------
      ! diagnostics
      !------------
      diag_trp_ei(:,:) = SUM(SUM( e_i (:,:,1:nlay_i,:) - e_i_b (:,:,1:nlay_i,:), dim=4 ), dim=3 ) * r1_rdtice
      diag_trp_es(:,:) = SUM(SUM( e_s (:,:,1:nlay_s,:) - e_s_b (:,:,1:nlay_s,:), dim=4 ), dim=3 ) * r1_rdtice
      diag_trp_sv(:,:) = SUM(     sv_i(:,:,:)          - sv_i_b(:,:,:)                  , dim=3 ) * r1_rdtice
      diag_trp_vi(:,:) = SUM(     v_i (:,:,:)          - v_i_b (:,:,:)                  , dim=3 ) * r1_rdtice
      diag_trp_vs(:,:) = SUM(     v_s (:,:,:)          - v_s_b (:,:,:)                  , dim=3 ) * r1_rdtice
      IF( iom_use('icemtrp') )   CALL iom_put( 'icemtrp' ,  diag_trp_vi * rhoi          )   ! ice mass transport
      IF( iom_use('snwmtrp') )   CALL iom_put( 'snwmtrp' ,  diag_trp_vs * rhos          )   ! snw mass transport
      IF( iom_use('salmtrp') )   CALL iom_put( 'salmtrp' ,  diag_trp_sv * rhoi * 1.e-03 )   ! salt mass transport (kg/m2/s)
      IF( iom_use('dihctrp') )   CALL iom_put( 'dihctrp' , -diag_trp_ei                 )   ! advected ice heat content (W/m2)
      IF( iom_use('dshctrp') )   CALL iom_put( 'dshctrp' , -diag_trp_es                 )   ! advected snw heat content (W/m2)

      ! controls
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icedyn_adv', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icectl    )   CALL ice_prt     (kt, iiceprt, jiceprt,-1, ' - ice dyn & trp - ')                           ! prints
      IF( ln_timing    )   CALL timing_stop ('icedyn_adv')                                                             ! timing
      !
   END SUBROUTINE ice_dyn_adv


   SUBROUTINE ice_dyn_adv_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_adv_init  ***
      !!
      !! ** Purpose :   Physical constants and parameters linked to the ice
      !!                dynamics
      !!
      !! ** Method  :   Read the namdyn_adv namelist and check the ice-dynamic
      !!                parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namdyn_adv
      !!-------------------------------------------------------------------
      INTEGER ::   ios, ioptio   ! Local integer output status for namelist read
      !!
      NAMELIST/namdyn_adv/ ln_adv_Pra, ln_adv_UMx, nn_UMx
      !!-------------------------------------------------------------------
      !
      REWIND( numnam_ice_ref )         ! Namelist namdyn_adv in reference namelist : Ice dynamics
      READ  ( numnam_ice_ref, namdyn_adv, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdyn_adv in reference namelist' )
      REWIND( numnam_ice_cfg )         ! Namelist namdyn_adv in configuration namelist : Ice dynamics
      READ  ( numnam_ice_cfg, namdyn_adv, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdyn_adv in configuration namelist' )
      IF(lwm) WRITE( numoni, namdyn_adv )
      !
      IF(lwp) THEN                     ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_adv_init: ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_adv:'
         WRITE(numout,*) '      type of advection scheme (Prather)             ln_adv_Pra = ', ln_adv_Pra 
         WRITE(numout,*) '      type of advection scheme (Ulimate-Macho)       ln_adv_UMx = ', ln_adv_UMx 
         WRITE(numout,*) '         order of the Ultimate-Macho scheme          nn_UMx     = ', nn_UMx
      ENDIF
      !
      !                             !== set the choice of ice advection ==!
      ioptio = 0 
      IF( ln_adv_Pra ) THEN   ;   ioptio = ioptio + 1   ;   nice_adv = np_advPRA    ;   ENDIF
      IF( ln_adv_UMx ) THEN   ;   ioptio = ioptio + 1   ;   nice_adv = np_advUMx    ;   ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'ice_dyn_adv_init: choose one and only one ice adv. scheme (ln_adv_Pra or ln_adv_UMx)' )
      !
      IF( ln_adv_Pra )   CALL adv_pra_init  !* read or initialize all required files
      !
   END SUBROUTINE ice_dyn_adv_init


   !!======================================================================
END MODULE icedyn_adv

