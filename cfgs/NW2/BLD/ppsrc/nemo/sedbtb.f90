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

MODULE sedbtb
   !!======================================================================
   !!              ***  MODULE  sedbtb  ***
   !!    Sediment : bioturbation of the solid components
   !!=====================================================================
   !! * Modules used
   USE sed     ! sediment global variable
   USE sedmat  ! linear system of equations
   USE lib_mpp         ! distribued memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC sed_btb


   !! $Id: sedbtb.F90 10222 2018-10-25 09:42:23Z aumont $
CONTAINS
   
   SUBROUTINE sed_btb( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sed_btb  ***
      !!
      !! ** Purpose :  performs bioturbation of the solid sediment components
      !!
      !! ** Method  :  ``diffusion'' of solid sediment components. 
      !!
      !!   History :
      !!        !  98-08 (E. Maier-Reimer, Christoph Heinze )  Original code
      !!        !  04-10 (N. Emprin, M. Gehlen ) F90
      !!        !  06-04 (C. Ethe)  Re-organization
      !!----------------------------------------------------------------------
      !!* Arguments
      INTEGER, INTENT(in) ::  kt              ! time step

      ! * local variables
      INTEGER :: ji, jk, js
      REAL(wp), DIMENSION(jpoce,jpksedm1,jpsol) ::  zsol  !   solution
      !------------------------------------------------------------------------

      IF( ln_timing )  CALL timing_start('sed_btb')

      IF( kt == nitsed000 ) THEN
         IF (lwp) WRITE(numsed,*) ' sed_btb : Bioturbation  '
         IF (lwp) WRITE(numsed,*) ' '
      ENDIF

      ! Initializations
      !----------------
      zsol(:,:,:) = 0.

      ! right hand side of coefficient matrix
      !--------------------------------------
      DO js = 1, jpsol
         DO jk = 1, jpksedm1
            DO ji = 1, jpoce
               zsol(ji,jk,js) = solcp(ji,jk+1,js)
            ENDDO
         ENDDO
      ENDDO

      CALL sed_mat( jpsol, jpoce, jpksedm1, zsol, dtsed / 2.0 )


      ! store solution of the tridiagonal system
      !------------------------
      DO js = 1, jpsol
         DO jk = 1, jpksedm1
            DO ji = 1, jpoce
               solcp(ji,jk+1,js) = zsol(ji,jk,js)
            ENDDO
         ENDDO
      ENDDO

      IF( ln_timing )  CALL timing_stop('sed_btb')

   END SUBROUTINE sed_btb

END MODULE sedbtb
