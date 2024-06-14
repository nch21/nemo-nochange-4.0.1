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

MODULE iscplini
   !!======================================================================
   !!                       ***  MODULE  sbciscpl  ***
   !! Ocean forcing:  ?????
   !!=====================================================================
   !! History :  NEMO  ! 2015-01 P. Mathiot: original 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   iscpl_init     : initialisation routine (namelist)
   !!   iscpl_alloc    : allocation of correction variables
   !!----------------------------------------------------------------------
   USE oce             ! global tra/dyn variable
   USE dom_oce         ! ocean space and time domain
   !
   USE lib_mpp         ! MPP library
   USE lib_fortran     ! MPP library
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   iscpl_init      
   PUBLIC   iscpl_alloc 
   
   !                                 !!* namsbc_iscpl namelist *
   LOGICAL , PUBLIC ::   ln_hsb       !:
   INTEGER , PUBLIC ::   nn_fiscpl    !:
   INTEGER , PUBLIC ::   nn_drown     !:
   
   INTEGER , PUBLIC ::   nstp_iscpl   !:
   REAL(wp), PUBLIC ::   rdt_iscpl    !: 
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hdiv_iscpl   !:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   htsc_iscpl   !:

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: iscplini.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION iscpl_alloc()
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE sbc_iscpl_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( htsc_iscpl(jpi,jpj,jpk,jpts) , hdiv_iscpl(jpi,jpj,jpk) , STAT=iscpl_alloc )
         !
      CALL mpp_sum ( 'iscplini', iscpl_alloc )
      IF( iscpl_alloc > 0 )   CALL ctl_warn('iscpl_alloc: allocation of arrays failed')
   END FUNCTION iscpl_alloc


   SUBROUTINE iscpl_init()
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER ::   ios           ! Local integer output status for namelist read
      NAMELIST/namsbc_iscpl/ nn_fiscpl, ln_hsb, nn_drown
      !!----------------------------------------------------------------------
      !
      nn_fiscpl = 0
      ln_hsb    = .FALSE.
      REWIND( numnam_ref )              ! Namelist namsbc_iscpl in reference namelist : Ice sheet coupling
      READ  ( numnam_ref, namsbc_iscpl, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_iscpl in reference namelist' )
      REWIND( numnam_cfg )              ! Namelist namsbc_iscpl in configuration namelist : Ice Sheet coupling
      READ  ( numnam_cfg, namsbc_iscpl, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_iscpl in configuration namelist' )
      IF(lwm) WRITE ( numond, namsbc_iscpl )
      !
      nstp_iscpl=MIN( nn_fiscpl, nitend-nit000+1 ) ! the coupling period have to be less or egal than the total number of time step
      rdt_iscpl = nstp_iscpl * rn_rdt
      !
      IF (lwp) THEN
         WRITE(numout,*) 'iscpl_rst:'
         WRITE(numout,*) '~~~~~~~~~'
         WRITE(numout,*) ' coupling     flag (ln_iscpl )            = ', ln_iscpl
         WRITE(numout,*) ' conservation flag (ln_hsb   )            = ', ln_hsb
         WRITE(numout,*) ' nb of stp for cons (rn_fiscpl)           = ', nstp_iscpl
         IF (nstp_iscpl .NE. nn_fiscpl) WRITE(numout,*) 'W A R N I N G: nb of stp for cons has been modified &
            &                                           (larger than run length)'
         WRITE(numout,*) ' coupling time step                       = ', rdt_iscpl
         WRITE(numout,*) ' number of call of the extrapolation loop = ', nn_drown
      ENDIF
      !
   END SUBROUTINE iscpl_init

   !!======================================================================
END MODULE iscplini
