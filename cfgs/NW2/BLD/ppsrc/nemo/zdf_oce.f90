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

MODULE zdf_oce
   !!======================================================================
   !!              ***  MODULE  zdf_oce  ***
   !! Ocean physics : define vertical mixing variables
   !!=====================================================================
   !! history :  1.0  !  2002-06  (G. Madec)  Original code
   !!            3.2  !  2009-07  (G. Madec)  addition of avm
   !!            4.0  !  2017-05  (G. Madec)  avm and drag coef. defined at t-point
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC  zdf_oce_alloc    ! Called in nemogcm.F90

   !                            !!* namelist namzdf: vertical physics *
   !                             ! Adaptive-implicit vertical advection flag
   LOGICAL , PUBLIC ::   ln_zad_Aimp !: adaptive (Courant number-based) implicit vertical advection
   !                             ! vertical closure scheme flags
   LOGICAL , PUBLIC ::   ln_zdfcst   !: constant coefficients
   LOGICAL , PUBLIC ::   ln_zdfric   !: Richardson depend coefficients
   LOGICAL , PUBLIC ::   ln_zdftke   !: Turbulent Kinetic Energy closure
   LOGICAL , PUBLIC ::   ln_zdfgls   !: Generic Length Scale closure
   LOGICAL , PUBLIC ::   ln_zdfosm   !: OSMOSIS BL closure
   !                             ! convection
   LOGICAL , PUBLIC ::   ln_zdfevd   !: convection: enhanced vertical diffusion flag
   INTEGER , PUBLIC ::      nn_evdm     !: =0/1 flag to apply enhanced avm or not
   REAL(wp), PUBLIC ::      rn_evd      !: vertical eddy coeff. for enhanced vert. diff. (m2/s)
   LOGICAL , PUBLIC ::   ln_zdfnpc   !: convection: non-penetrative convection flag
   INTEGER , PUBLIC ::      nn_npc      !: non penetrative convective scheme call  frequency
   INTEGER , PUBLIC ::      nn_npcp     !: non penetrative convective scheme print frequency
   !                             ! double diffusion
   LOGICAL , PUBLIC ::   ln_zdfddm   !: double diffusive mixing flag
   REAL(wp), PUBLIC ::      rn_avts     !: maximum value of avs for salt fingering
   REAL(wp), PUBLIC ::      rn_hsbfr    !: heat/salt buoyancy flux ratio
   !                             ! gravity wave-induced vertical mixing
   LOGICAL , PUBLIC ::   ln_zdfswm   !: surface  wave-induced mixing flag
   LOGICAL , PUBLIC ::   ln_zdfiwm   !: internal wave-induced mixing flag
   !                             ! coefficients 
   REAL(wp), PUBLIC ::   rn_avm0     !: vertical eddy viscosity (m2/s)
   REAL(wp), PUBLIC ::   rn_avt0     !: vertical eddy diffusivity (m2/s)
   INTEGER , PUBLIC ::   nn_avb      !: constant or profile background on avt (=0/1)
   INTEGER , PUBLIC ::   nn_havtb    !: horizontal shape or not for avtb (=0/1)   !                             ! convection


   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   avm, avt, avs  !: vertical mixing coefficients (w-point) [m2/s]
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   avm_k , avt_k  !: Kz computed by turbulent closure alone [m2/s]
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   en             !: now turbulent kinetic energy          [m2/s2]
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:)     ::   avmb , avtb    !: background profile of avm and avt      [m2/s]
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)   ::   avtb_2d        !: horizontal shape of background Kz profile [-]

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: zdf_oce.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_oce_alloc()
      !!----------------------------------------------------------------------
      !!            *** FUNCTION zdf_oce_alloc ***
      !!----------------------------------------------------------------------
      !
      ALLOCATE( avm (jpi,jpj,jpk) , avm_k(jpi,jpj,jpk) , avs(jpi,jpj,jpk) ,   &
         &      avt (jpi,jpj,jpk) , avt_k(jpi,jpj,jpk) , en (jpi,jpj,jpk) ,   & 
         &      avmb(jpk)         , avtb(jpk)          , avtb_2d(jpi,jpj) , STAT = zdf_oce_alloc )
         !
      IF( zdf_oce_alloc /= 0 )   CALL ctl_stop( 'STOP', 'zdf_oce_alloc: failed to allocate arrays' )
      !
   END FUNCTION zdf_oce_alloc

   !!======================================================================
END MODULE zdf_oce
