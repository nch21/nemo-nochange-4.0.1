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

MODULE icedyn_rhg
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg  ***
   !!   Sea-Ice dynamics : master routine for rheology 
   !!======================================================================
   !! history :  4.0  !  2018     (C. Rousset)      Original code
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!    ice_dyn_rhg      : computes ice velocities
   !!    ice_dyn_rhg_init : initialization and namelist read
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain
   USE ice            ! sea-ice: variables
   USE icedyn_rhg_evp ! sea-ice: EVP rheology
   USE icectl         ! sea-ice: control prints
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_rhg        ! called by icestp.F90
   PUBLIC   ice_dyn_rhg_init   ! called by icestp.F90

   INTEGER ::              nice_rhg   ! choice of the type of rheology
   !                                        ! associated indices:
   INTEGER, PARAMETER ::   np_rhgEVP = 1   ! EVP rheology
!! INTEGER, PARAMETER ::   np_rhgEAP = 2   ! EAP rheology

   ! ** namelist (namrhg) **
   LOGICAL ::   ln_rhg_EVP       ! EVP rheology
   !
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
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icedyn_rhg.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_rhg( kt )
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE ice_dyn_rhg  ***
      !!               
      !! ** Purpose :   compute ice velocity
      !!
      !! ** Action  : comupte - ice velocity (u_ice, v_ice)
      !!                      - 3 components of the stress tensor (stress1_i, stress2_i, stress12_i)
      !!                      - shear, divergence and delta (shear_i, divu_i, delta_i)
      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ice time step
      !
      INTEGER  ::   jl   ! dummy loop indices
      !!--------------------------------------------------------------------
      ! controls
      IF( ln_timing    )   CALL timing_start('icedyn_rhg')                                                             ! timing
      IF( ln_icediachk )   CALL ice_cons_hsm(0, 'icedyn_rhg', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icediachk )   CALL ice_cons2D  (0, 'icedyn_rhg',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_dyn_rhg: sea-ice rheology'
         WRITE(numout,*)'~~~~~~~~~~~'
      ENDIF
      !
      !--------------!
      !== Rheology ==!
      !--------------!   
      SELECT CASE( nice_rhg )
      !                                !------------------------!
      CASE( np_rhgEVP )                ! Elasto-Viscous-Plastic !
         !                             !------------------------!
         CALL ice_dyn_rhg_evp( kt, stress1_i, stress2_i, stress12_i, shear_i, divu_i, delta_i )
         !         
      END SELECT
      !
      IF( lrst_ice ) THEN                       !* write EVP fields in the restart file
         IF( ln_rhg_EVP )   CALL rhg_evp_rst( 'WRITE', kt )
      ENDIF
      !
      ! controls
      IF( ln_ctl       )   CALL ice_prt3D   ('icedyn_rhg')                                                             ! prints
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icedyn_rhg', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icediachk )   CALL ice_cons2D  (1, 'icedyn_rhg',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation
      IF( ln_timing    )   CALL timing_stop ('icedyn_rhg')                                                             ! timing
      !
   END SUBROUTINE ice_dyn_rhg


   SUBROUTINE ice_dyn_rhg_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_rhg_init  ***
      !!
      !! ** Purpose : Physical constants and parameters linked to the ice
      !!      dynamics
      !!
      !! ** Method  :  Read the namdyn_rhg namelist and check the ice-dynamic
      !!       parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namdyn_rhg
      !!-------------------------------------------------------------------
      INTEGER ::   ios, ioptio   ! Local integer output status for namelist read
      !!
      NAMELIST/namdyn_rhg/  ln_rhg_EVP, ln_aEVP, rn_creepl, rn_ecc , nn_nevp, rn_relast
      !!-------------------------------------------------------------------
      !
      REWIND( numnam_ice_ref )         ! Namelist namdyn_rhg in reference namelist : Ice dynamics
      READ  ( numnam_ice_ref, namdyn_rhg, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdyn_rhg in reference namelist' )
      REWIND( numnam_ice_cfg )         ! Namelist namdyn_rhg in configuration namelist : Ice dynamics
      READ  ( numnam_ice_cfg, namdyn_rhg, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namdyn_rhg in configuration namelist' )
      IF(lwm) WRITE ( numoni, namdyn_rhg )
      !
      IF(lwp) THEN                     ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_rhg_init: ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist : namdyn_rhg:'
         WRITE(numout,*) '      rheology EVP (icedyn_rhg_evp)                        ln_rhg_EVP = ', ln_rhg_EVP
         WRITE(numout,*) '         use adaptive EVP (aEVP)                           ln_aEVP    = ', ln_aEVP
         WRITE(numout,*) '         creep limit                                       rn_creepl  = ', rn_creepl
         WRITE(numout,*) '         eccentricity of the elliptical yield curve        rn_ecc     = ', rn_ecc
         WRITE(numout,*) '         number of iterations for subcycling               nn_nevp    = ', nn_nevp
         WRITE(numout,*) '         ratio of elastic timescale over ice time step     rn_relast  = ', rn_relast
      ENDIF
      !
      !                             !== set the choice of ice advection ==!
      ioptio = 0 
      IF( ln_rhg_EVP ) THEN   ;   ioptio = ioptio + 1   ;   nice_rhg = np_rhgEVP    ;   ENDIF
!!    IF( ln_rhg_EAP ) THEN   ;   ioptio = ioptio + 1   ;   nice_rhg = np_rhgEAP    ;   ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'ice_dyn_rhg_init: choose one and only one ice rheology' )
      !
      IF( ln_rhg_EVP  )   CALL rhg_evp_rst( 'READ' )  !* read or initialize all required files
      !
   END SUBROUTINE ice_dyn_rhg_init


   !!======================================================================
END MODULE icedyn_rhg
