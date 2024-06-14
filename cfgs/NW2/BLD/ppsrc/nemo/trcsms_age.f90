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

MODULE trcsms_age
   !!======================================================================
   !!                         ***  MODULE trcsms_age  ***
   !! TOP :   Main module of the AGE tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
   !! trc_sms_age       : AGE model main routine
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_age       ! called by trcsms.F90 module

   INTEGER , PUBLIC :: nl_age             ! T level surrounding age_depth
   INTEGER , PUBLIC :: nla_age            ! T level wholly above age_depth
   INTEGER , PUBLIC :: nlb_age            ! T level wholly below age_depth

   REAL(wp), PUBLIC :: rn_age_depth       ! = 10       depth over which age tracer reset to zero
   REAL(wp), PUBLIC :: rn_age_kill_rate   ! = -1./7200  recip of relaxation timescale (s) for  age tracer shallower than age_depth
   
   REAL(wp), PUBLIC :: rryear          !: recip number of seconds in one year
   REAL(wp), PUBLIC :: frac_kill_age   !: fraction of level nl_age above age_depth where it is relaxed towards zero
   REAL(wp), PUBLIC :: frac_add_age    !: fraction of level nl_age below age_depth where it is incremented


   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsms_age.F90 10070 2018-08-28 14:30:54Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_age( kt )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_age  ***
      !!
      !! ** Purpose :   main routine of AGE model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER ::   jn, jk   ! dummy loop index
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_age')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_age:  AGE model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'


      DO jk = 1, nla_age
         tra(:,:,jk,jp_age) = rn_age_kill_rate * trb(:,:,jk,jp_age)
      END DO
      !
      tra(:,:,nl_age,jp_age) = frac_kill_age * rn_age_kill_rate * trb(:,:,nl_age,jp_age)  &
          &                   + frac_add_age  * rryear * tmask(:,:,nl_age)
      !
      DO jk = nlb_age, jpk
         tra(:,:,jk,jp_age) = tmask(:,:,jk) * rryear
      END DO
      !
      IF( l_trdtrc ) CALL trd_trc( tra(:,:,:,jp_age), jn, jptra_sms, kt )   ! save trends
      !
      IF( ln_timing )   CALL timing_stop('trc_sms_age')
      !
   END SUBROUTINE trc_sms_age

   !!======================================================================
END MODULE trcsms_age
