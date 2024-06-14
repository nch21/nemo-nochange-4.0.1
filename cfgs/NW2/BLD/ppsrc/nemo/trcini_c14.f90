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

MODULE trcini_c14
   !!======================================================================
   !!                         ***  MODULE trcini_c14  ***
   !! TOP :   initialisation of the C14 tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !! History :   3.0  !  2015 (A. Mouchet) C14 Code
   !!----------------------------------------------------------------------
   !! trc_ini_c14   : C14 model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE sms_c14
   USE trcatm_c14
   USE trcnam_c14
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_c14   ! called by trcini.F90 module
   
   !

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcini_c14.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_c14
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_c14  ***  
      !!
      !! ** Purpose :   initialization for C14 model
      !!
      !! ** Method  : 
      !!----------------------------------------------------------------------
      !
      REAL(wp) :: ztrai
      INTEGER  :: jn
      CHARACTER(len = 20)  ::  cltra
      !!----------------------------------------------------------------------
      !
      CALL trc_nam_c14
      !                       ! Allocate c14 arrays
      IF( sms_c14_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_c14: unable to allocate C14 arrays' )

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_c14: initialisation of C14 model'
      !
      IF( .NOT. ln_rsttr )  THEN
         !
         IF(lwp) WRITE(numout,*) '                      ==>    PRESCRIBED initial VALUES'
         IF(lwp) WRITE(numout,*) '                      ==>    Ocean C14/C :', rc14init 
         !
         trn(:,:,:,jp_c14) = rc14init * tmask(:,:,:)
         !
         qtr_c14(:,:) = 0._wp           ! Init of air-sea BC
         !
      ELSE

        IF(lwp) WRITE(numout,*)
        IF(lwp) WRITE(numout,*) ' trc_rst_read_c14 : Read specific variables for c14 model '
        IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
        ! 
        CALL iom_get( numrtr, 'co2sbc', co2sbc ) 
        CALL iom_get( numrtr, jpdom_autoglo, 'c14sbc', c14sbc ) 
        CALL iom_get( numrtr, jpdom_autoglo, 'exch_co2', exch_co2 ) 
        CALL iom_get( numrtr, jpdom_autoglo, 'exch_c14', exch_c14 ) 
        CALL iom_get( numrtr, jpdom_autoglo, 'qtr_c14', qtr_c14 )
        !
      END IF
      !
      IF( ( nn_rsttr == 0 ) .OR. ( .NOT. ln_rsttr ) ) THEN
      !
      !                         !  qint set to zero <=== Initial of transient
      !                         !                   <=== Restart=false
         IF(lwp) WRITE(numout,*) '                    ==>    qint reset to ZERO '
         qint_c14(:,:) = 0._wp
      !
      ELSE
        !
        CALL iom_get( numrtr, jpdom_autoglo, 'qint_c14', qint_c14 ) 
        !
      ENDIF
      !
      CALL trc_atm_c14_ini   ! Init atm values
      !
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      !
   END SUBROUTINE trc_ini_c14

   !!======================================================================
END MODULE trcini_c14
