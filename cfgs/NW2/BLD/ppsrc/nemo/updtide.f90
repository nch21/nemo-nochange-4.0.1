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

MODULE updtide
   !!======================================================================
   !!                       ***  MODULE  updtide  ***
   !! Initialization of tidal forcing
   !!======================================================================
   !! History :  9.0  !  07  (O. Le Galloudec)  Original code
   !!----------------------------------------------------------------------
   !!   upd_tide       : update tidal potential
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O units
   USE phycst          ! physical constant
   USE sbctide         ! tide potential variable
   USE tideini, ONLY: ln_tide_ramp, rdttideramp

   IMPLICIT NONE
   PUBLIC

   PUBLIC   upd_tide   ! called in dynspg_... modules
  
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: updtide.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE upd_tide( kt, kit, kt_offset )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE upd_tide  ***
      !!
      !! ** Purpose :   provide at each time step the astronomical potential
      !!
      !! ** Method  :   computed from pulsation and amplitude of all tide components
      !!
      !! ** Action  :   pot_astro   actronomical potential
      !!----------------------------------------------------------------------      
      INTEGER, INTENT(in)           ::   kt      ! ocean time-step index
      INTEGER, INTENT(in), OPTIONAL ::   kit     ! external mode sub-time-step index (lk_dynspg_ts=T)
      INTEGER, INTENT(in), OPTIONAL ::   kt_offset ! time offset in number 
                                                     ! of internal steps             (lk_dynspg_ts=F)
                                                     ! of external steps             (lk_dynspg_ts=T)
      !
      INTEGER  ::   ioffset      ! local integer
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zt, zramp    ! local scalar
      REAL(wp), DIMENSION(nb_harmo) ::   zwt 
      !!----------------------------------------------------------------------      
      !
      !                               ! tide pulsation at model time step (or sub-time-step)
      zt = ( kt - kt_tide ) * rdt
      !
      ioffset = 0
      IF( PRESENT( kt_offset ) )   ioffset = kt_offset
      !
      IF( PRESENT( kit ) )   THEN
         zt = zt + ( kit +  ioffset - 1 ) * rdt / REAL( nn_baro, wp )
      ELSE
         zt = zt + ioffset * rdt
      ENDIF
      !
      zwt(:) = omega_tide(:) * zt

      pot_astro(:,:) = 0._wp          ! update tidal potential (sum of all harmonics)
      DO jk = 1, nb_harmo   
         pot_astro(:,:) = pot_astro(:,:) + amp_pot(:,:,jk) * COS( zwt(jk) + phi_pot(:,:,jk) )      
      END DO
      !
      IF( ln_tide_ramp ) THEN         ! linear increase if asked
         zt = ( kt - nit000 ) * rdt
         IF( PRESENT( kit ) )   zt = zt + ( kit + ioffset -1) * rdt / REAL( nn_baro, wp )
         zramp = MIN(  MAX( zt / (rdttideramp*rday) , 0._wp ) , 1._wp  )
         pot_astro(:,:) = zramp * pot_astro(:,:)
      ENDIF
      !
   END SUBROUTINE upd_tide

  !!======================================================================

END MODULE updtide
