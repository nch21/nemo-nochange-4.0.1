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

MODULE trcini_age
   !!======================================================================
   !!                         ***  MODULE trcini_age  ***
   !! TOP :   initialisation of the AGE tracer
   !!======================================================================
   !! History :   2.0  !  2007-12  (G. Nurser, G. Madec, C. Ethe ) Original code
   !!----------------------------------------------------------------------
   !! trc_ini_age   : MY_TRC model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc
   USE trc
   USE trcnam_age
   USE trcsms_age

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_age   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcini_age.F90 10070 2018-08-28 14:30:54Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_age
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_age  ***  
      !!
      !! ** Purpose :   initialization for AGE model
      !!
      !!----------------------------------------------------------------------
      INTEGER    ::  jn
      CHARACTER(len = 20)  ::  cltra
      !!----------------------------------------------------------------------
      !
      CALL trc_nam_age
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_age: passive tracer age'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*)

      rryear  = 1._wp / ( nyear_len(1) * rday )    ! recip number of seconds in one year

      !! BUG in s-coordinate this does not work!
      nlb_age = MINLOC( gdepw_1d, mask = gdepw_1d > rn_age_depth, dim = 1 ) ! shallowest W level Below age_depth
                                                                            !  = shallowest T level wholly below age_depth
      nl_age  = nlb_age - 1                                                 ! deepest    W level Above age_depth
                                                                            !  = T level surrounding age_depth

      nla_age = nl_age - 1                                                   ! deepest    T level wholly above age_depth

      frac_kill_age = ( rn_age_depth - gdepw_1d(nl_age) ) / e3t_1d(nl_age)      ! fraction of level nl_age above age_depth
      frac_add_age  = 1._wp -  frac_kill_age                                    ! fraction of level nl_age below age_depth

      
      IF( .NOT. ln_rsttr ) trn(:,:,:,jp_age) = 0.
      !
   END SUBROUTINE trc_ini_age

   !!======================================================================
END MODULE trcini_age
