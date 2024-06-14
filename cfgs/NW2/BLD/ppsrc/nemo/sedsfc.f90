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

MODULE sedsfc
   !!======================================================================
   !!              ***  MODULE  sedsfc  ***
   !!    Sediment : Data at sediment surface
   !!=====================================================================
   !! * Modules used
   USE sed     ! sediment global variable
   USE sedarr
   USE seddta

   PUBLIC sed_sfc

   !! $Id: sedsfc.F90 10222 2018-10-25 09:42:23Z aumont $
CONTAINS

   SUBROUTINE sed_sfc( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sed_sfc ***
      !!
      !! ** Purpose :  Give data from sediment model to tracer model
      !!
      !!
      !!   History :
      !!        !  06-04 (C. Ethe)  Orginal code
      !!----------------------------------------------------------------------
      !!* Arguments
      INTEGER, INTENT(in) ::  kt              ! time step

      ! * local variables
      INTEGER :: ji, jj, ikt     ! dummy loop indices

      !------------------------------------------------------------------------
      ! reading variables

      IF( ln_timing )  CALL timing_start('sed_sfc')

      CALL unpack_arr ( jpoce, trc_data(1:jpi,1:jpj,1), iarroce(1:jpoce), pwcp(1:jpoce,1,jwalk) )
      CALL unpack_arr ( jpoce, trc_data(1:jpi,1:jpj,2), iarroce(1:jpoce), pwcp(1:jpoce,1,jwdic) )
      CALL unpack_arr ( jpoce, trc_data(1:jpi,1:jpj,3), iarroce(1:jpoce), pwcp(1:jpoce,1,jwno3) )
      CALL unpack_arr ( jpoce, trc_data(1:jpi,1:jpj,4), iarroce(1:jpoce), pwcp(1:jpoce,1,jwpo4) )
      CALL unpack_arr ( jpoce, trc_data(1:jpi,1:jpj,5), iarroce(1:jpoce), pwcp(1:jpoce,1,jwoxy) )
      CALL unpack_arr ( jpoce, trc_data(1:jpi,1:jpj,6), iarroce(1:jpoce), pwcp(1:jpoce,1,jwsil) )
      CALL unpack_arr ( jpoce, trc_data(1:jpi,1:jpj,7), iarroce(1:jpoce), pwcp(1:jpoce,1,jwnh4) )
      CALL unpack_arr ( jpoce, trc_data(1:jpi,1:jpj,8), iarroce(1:jpoce), pwcp(1:jpoce,1,jwfe2) )


      DO jj = 1,jpj
         DO ji = 1, jpi
            ikt = mbkt(ji,jj)
            IF ( tmask(ji,jj,ikt) == 1 ) THEN
               trb(ji,jj,ikt,jptal) = trc_data(ji,jj,1)
               trb(ji,jj,ikt,jpdic) = trc_data(ji,jj,2)
               trb(ji,jj,ikt,jpno3) = trc_data(ji,jj,3) * 7.625
               trb(ji,jj,ikt,jppo4) = trc_data(ji,jj,4) * 122.
               trb(ji,jj,ikt,jpoxy) = trc_data(ji,jj,5)
               trb(ji,jj,ikt,jpsil) = trc_data(ji,jj,6)
               trb(ji,jj,ikt,jpnh4) = trc_data(ji,jj,7) * 7.625
               trb(ji,jj,ikt,jpfer) = trc_data(ji,jj,8)
            ENDIF
         ENDDO
      ENDDO

      IF( ln_timing )  CALL timing_stop('sed_sfc')

   END SUBROUTINE sed_sfc

END MODULE sedsfc
