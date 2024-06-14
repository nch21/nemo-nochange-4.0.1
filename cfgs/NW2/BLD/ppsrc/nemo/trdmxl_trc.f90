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

MODULE trdmxl_trc
   !!======================================================================
   !!                       ***  MODULE  trdmxl_trc  ***
   !! Ocean diagnostics:  mixed layer passive tracer trends 
   !!======================================================================
   !! History :  9.0  !  06-08  (C. Deltel)  Original code (from trdmxl.F90)
   !!                 !  07-04  (C. Deltel)  Bug fix : add trcrad trends
   !!                 !  07-06  (C. Deltel)  key_gyre : do not call lbc_lnk
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trd_mxl_trc( kt )                                   ! Empty routine
      INTEGER, INTENT( in) ::   kt
      WRITE(*,*) 'trd_mxl_trc: You should not have seen this print! error?', kt
   END SUBROUTINE trd_mxl_trc
   SUBROUTINE trd_mxl_trc_zint( ptrc_trdmxl, ktrd, ctype, kjn )
      INTEGER               , INTENT( in ) ::  ktrd, kjn              ! ocean trend index and passive tracer rank
      CHARACTER(len=2)      , INTENT( in ) ::  ctype                  ! surface/bottom (2D) or interior (3D) physics
      REAL, DIMENSION(:,:,:), INTENT( in ) ::  ptrc_trdmxl            ! passive trc trend
      WRITE(*,*) 'trd_mxl_trc_zint: You should not have seen this print! error?', ptrc_trdmxl(1,1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', ctype
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', ktrd
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kjn
   END SUBROUTINE trd_mxl_trc_zint
   SUBROUTINE trd_mxl_trc_init                                    ! Empty routine
      WRITE(*,*) 'trd_mxl_trc_init: You should not have seen this print! error?'
   END SUBROUTINE trd_mxl_trc_init

   !!======================================================================
END MODULE trdmxl_trc
