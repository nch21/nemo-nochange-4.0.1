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

MODULE sedco3
   !!======================================================================
   !!              ***  MODULE  sedco3  ***
   !!    Sediment : carbonate in sediment pore water
   !!=====================================================================
   !! * Modules used
   USE sed     ! sediment global variable
   USE sedchem
   USE lib_mpp         ! distribued memory computing library


   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC sed_co3     

   !!----------------------------------------------------------------------
   !!   OPA 9.0   !   LODYC-IPSL   (2003)
   !!----------------------------------------------------------------------

   !! $Id: sedco3.F90 10222 2018-10-25 09:42:23Z aumont $
CONTAINS


   SUBROUTINE sed_co3( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_co3  ***
      !!
      !! ** Purpose :  carbonate ion and proton concentration 
      !!               in sediment pore water
      !!
      !! ** Methode :  - solving nonlinear equation for [H+] with given alkalinity
      !!               and total co2 
      !!               - one dimensional newton-raphson algorithm for [H+])
      !!
      !!   History :
      !!        !  98-08 (E. Maier-Reimer, Christoph Heinze )  Original code
      !!        !  04-10 (N. Emprin, M. Gehlen ) coupled with PISCES
      !!        !  06-04 (C. Ethe)  Re-organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(in)  :: kt   ! time step
      !
      !---Local variables
      INTEGER  :: ji, jk           ! dummy loop indices

      REAL(wp), DIMENSION(jpoce,jpksed) :: zhinit, zhi
     !!----------------------------------------------------------------------

      IF( ln_timing )  CALL timing_start('sed_co3')

      IF( kt == nitsed000 ) THEN
         IF (lwp) WRITE(numsed,*) ' sed_co3 : carbonate ion and proton concentration calculation  '
         IF (lwp) WRITE(numsed,*) ' '
      ENDIF

      DO jk = 1, jpksed
         zhinit(:,jk)   = hipor(:,jk) / densSW(:)
      END DO

      !     -------------------------------------------
      !     COMPUTE [CO3--] and [H+] CONCENTRATIONS
      !     -------------------------------------------

      CALL solve_at_general_sed(zhinit, zhi)

      DO jk = 1, jpksed
         DO ji = 1, jpoce
            co3por(ji,jk) = pwcp(ji,jk,jwdic) * ak1s(ji) * ak2s(ji) / (zhi(ji,jk)**2   &
            &               + ak1s(ji) * zhi(ji,jk) + ak1s(ji) * ak2s(ji) + rtrn )
            hipor(ji,jk)  = zhi(ji,jk) * densSW(ji)
         END DO
      END DO

     IF( ln_timing )  CALL timing_stop('sed_co3')

   END SUBROUTINE sed_co3

END MODULE sedco3
