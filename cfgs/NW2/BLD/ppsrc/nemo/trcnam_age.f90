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

MODULE trcnam_age
   !!======================================================================
   !!                         ***  MODULE trcnam_age  ***
   !! TOP :   initialisation of some run parameters for Age tracer
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) 
   !!----------------------------------------------------------------------
   !! trc_nam_age      : AGE tracer initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE trc             ! Ocean variables
   USE trcsms_age      ! AGE specific variable

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_age   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcnam_age.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_nam_age
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE trc_nam_age  ***
      !!                 
      !! ** Purpose :   Definition some run parameter for AGE model
      !!
      !! ** input   :   Namelist namage
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !!
      NAMELIST/namage/ rn_age_depth, rn_age_kill_rate 
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Sea Age Tracer'
         WRITE(numout,*)
         WRITE(numout,*) 'trc_nam_age : Read namage namelist for Age passive tracer'
         WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      ! Variable setting
      ctrcnm    (jp_age) = 'Age'
      ctrcln    (jp_age) = 'Sea water age since surface contact'
      ctrcun    (jp_age) = 'year'
      ln_trc_ini(jp_age) = .false.
      ln_trc_sbc(jp_age) = .false.
      ln_trc_cbc(jp_age) = .false.
      ln_trc_obc(jp_age) = .false.
      !
      REWIND( numnat_ref )              ! Namelist namagedate in reference namelist : AGE parameters
      READ  ( numnat_ref, namage, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namage in reference namelist' )
      REWIND( numnat_cfg )              ! Namelist namagedate in configuration namelist : AGE parameters
      READ  ( numnat_cfg, namage, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namage in configuration namelist' )
      IF(lwm) WRITE ( numont, namage )
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*) '   Namelist : namage'
         WRITE(numout,*) '      depth over which age tracer reset to zero     rn_age_depth      = ', rn_age_depth 
         WRITE(numout,*) '      recip of relaxation timescale                 rn_age_kill_rate  = ', rn_age_kill_rate, '[s]'
         WRITE(numout,*) '      (for age tracer shallower than age_depth) '
      ENDIF
      !
   END SUBROUTINE trc_nam_age
   
   !!======================================================================
END MODULE trcnam_age
