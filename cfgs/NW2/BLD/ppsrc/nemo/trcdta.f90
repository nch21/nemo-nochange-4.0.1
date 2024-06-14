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

MODULE trcdta
   !!======================================================================
   !!                     ***  MODULE  trcdta  ***
   !! TOP :  reads passive tracer data 
   !!=====================================================================
   !! History :   1.0  !  2002-04  (O. Aumont)  original code
   !!              -   !  2004-03  (C. Ethe)  module
   !!              -   !  2005-03  (O. Aumont, A. El Moussaoui) F90
   !!            3.4   !  2010-11  (C. Ethe, G. Madec)  use of fldread + dynamical allocation 
   !!            3.5   !  2013-08  (M. Vichi)  generalization for other BGC models
   !!            3.6   !  2015-03  (T. Lovato) revisit code I/O
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                              NO 3D passive tracer data
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_dta( kt, sf_trcdta, ptrcfac, ptrcdta)        ! Empty routine
      WRITE(*,*) 'trc_dta: You should not have seen this print! error?', kt
   END SUBROUTINE trc_dta

   !!======================================================================
END MODULE trcdta
