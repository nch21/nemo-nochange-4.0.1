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

MODULE p2zsms
   !!======================================================================
   !!                         ***  MODULE p2zsms  ***
   !! TOP :   Time loop of LOBSTER model
   !!======================================================================
   !! History :   1.0  !            M. Levy
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   p2zsms        :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE oce_trc          !
   USE trc
   USE sms_pisces
   USE p2zbio
   USE p2zopt
   USE p2zsed
   USE p2zexp
   USE trd_oce
   USE trdtrc_oce
   USE trdtrc
   USE trdmxl_trc

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p2z_sms    ! called in p2zsms.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p2zsms.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p2z_sms( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_sms  ***
      !!
      !! ** Purpose :  Managment of the call to Biological sources and sinks 
      !!               routines of LOBSTER bio-model 
      !!
      !! ** Method  : - ???
      !! --------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !
      INTEGER ::   jn   ! dummy loop index
      !! --------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p2z_sms')
      !
      CALL p2z_opt( kt )      ! optical model
      CALL p2z_bio( kt )      ! biological model
      CALL p2z_sed( kt )      ! sedimentation model
      CALL p2z_exp( kt )      ! export
      !
      IF( l_trdtrc ) THEN
         DO jn = jp_pcs0, jp_pcs1
           CALL trd_trc( tra(:,:,:,jn), jn, jptra_sms, kt )   ! save trends
         END DO
      END IF
      !
      IF ( lwm .AND. kt == nittrc000 ) CALL FLUSH    ( numonp )     ! flush output namelist PISCES
      IF( ln_timing )   CALL timing_stop('p2z_sms')
      !
   END SUBROUTINE p2z_sms

   !!======================================================================
END MODULE p2zsms
