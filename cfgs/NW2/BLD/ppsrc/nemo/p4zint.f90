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

MODULE p4zint
   !!======================================================================
   !!                         ***  MODULE p4zint  ***
   !! TOP :   PISCES interpolation and computation of various accessory fields
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
   !!   p4z_int        :  interpolation and computation of various accessory fields
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_int  
   REAL(wp) ::   xksilim = 16.5e-6_wp   ! Half-saturation constant for the Si half-saturation constant computation

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zint.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_int( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_int  ***
      !!
      !! ** Purpose :   interpolation and computation of various accessory fields
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !
      INTEGER  :: ji, jj                 ! dummy loop indices
      REAL(wp) :: zvar                   ! local variable
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_int')
      !
      ! Computation of phyto and zoo metabolic rate
      ! -------------------------------------------
      tgfunc (:,:,:) = EXP( 0.063913 * tsn(:,:,:,jp_tem) )
      tgfunc2(:,:,:) = EXP( 0.07608  * tsn(:,:,:,jp_tem) )

      ! Computation of the silicon dependant half saturation  constant for silica uptake
      ! ---------------------------------------------------
      DO ji = 1, jpi
         DO jj = 1, jpj
            zvar = trb(ji,jj,1,jpsil) * trb(ji,jj,1,jpsil)
            xksimax(ji,jj) = MAX( xksimax(ji,jj), ( 1.+ 7.* zvar / ( xksilim * xksilim + zvar ) ) * 1e-6 )
         END DO
      END DO
      !
      IF( nday_year == nyear_len(1) ) THEN
         xksi   (:,:) = xksimax(:,:)
         xksimax(:,:) = 0._wp
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_int')
      !
   END SUBROUTINE p4z_int

   !!======================================================================
END MODULE p4zint
