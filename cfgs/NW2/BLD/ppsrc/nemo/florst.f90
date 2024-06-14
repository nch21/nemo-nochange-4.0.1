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

MODULE florst
   !!======================================================================
   !!                       ***  MODULE  florst  ***
   !! Ocean floats : write floats restart files
   !!======================================================================
   !!  History :  OPA  !  1999-09  (Y. Drillet)    : Original code
   !!              -   !  2000-06  (J.-M. Molines) : Profiling floats for CLS 
   !!   NEMO      1.0  !  2002-10  (A. Bozec)  F90 : Free form and module
   !!             3.2  !  2010-08  (slaw, cbricaud): netcdf outputs and others 
   !!----------------------------------------------------------------------
   USE flo_oce         ! ocean drifting floats
   USE dom_oce         ! ocean space and time domain
   USE lib_mpp         ! distribued memory computing library
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC flo_rst         ! routine called by floats.F90
   PUBLIC flo_rst_alloc   ! routine called by floats.F90

   INTEGER, ALLOCATABLE, DIMENSION(:) :: iperproc   ! 1D workspace

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: florst.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION flo_rst_alloc()
      !!-------------------------------------------------------------------
      !!                ***  FUNCTION flo_rst_alloc  ***
      !!-------------------------------------------------------------------
      ALLOCATE( iperproc(jpnij), STAT=flo_rst_alloc )
      !  
      CALL mpp_sum ( 'florst', flo_rst_alloc )
      IF( flo_rst_alloc /= 0 )   CALL ctl_stop( 'STOP', 'flo_rst_alloc: failed to allocate arrays.' )
   END FUNCTION flo_rst_alloc


   SUBROUTINE flo_rst( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE flo_rst ***
      !!             
      !! ** Purpose :  
      !!             
      !!      
      !! ** Method  :   The frequency of  ??? is nwritefl
      !!      
      !!----------------------------------------------------------------------
      INTEGER  :: kt                            ! time step
      !
      CHARACTER (len=80)       :: clname             ! restart filename
      INTEGER                  :: ic , jc , jpn ,jfl ! temporary integer
      INTEGER                  :: inum               ! temporary logical unit for restart file
      !!----------------------------------------------------------------------

      IF(  ( MOD(kt,nn_stockfl) == 0 ) .OR. ( kt == nitend )  )THEN      

         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'flo_rst : write in  restart_float file '
            WRITE(numout,*) '~~~~~~~    '
         ENDIF

         ! file is opened and closed every time it is used.

         clname = 'restart.float.'
         ic = 1
         DO jc = 1, 16
            IF( cexper(jc:jc) /= ' ' ) ic = jc
         END DO
         clname = clname(1:14)//cexper(1:ic)
         ic = 1
         DO jc = 1, 48
            IF( clname(jc:jc) /= ' ' ) ic = jc
         END DO

         inum=0
         IF( lwp )THEN
            CALL ctl_opn( inum, clname, 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
            REWIND inum
         ENDIF
         !
         DO jpn = 1, jpnij
            iperproc(jpn) = 0
         END DO
         !
         IF(lwp) THEN
            REWIND(inum)
            WRITE (inum,*) tpifl,tpjfl,tpkfl,nisobfl,ngrpfl
            CLOSE (inum)
         ENDIF
         !
         ! Compute the number of trajectories for each processor
         !
         IF( lk_mpp ) THEN
            DO jfl = 1, jpnfl
               IF( (INT(tpifl(jfl)) >= mig(nldi)) .AND.   &
                  &(INT(tpifl(jfl)) <= mig(nlei)) .AND.   &
                  &(INT(tpjfl(jfl)) >= mjg(nldj)) .AND.   &
                  &(INT(tpjfl(jfl)) <= mjg(nlej)) ) THEN
                  iperproc(narea) = iperproc(narea)+1
               ENDIF
            END DO
            CALL mpp_sum( 'florst', iperproc, jpnij )
            !
            IF(lwp) THEN
               WRITE(numout,*) 'DATE',adatrj
               DO jpn = 1, jpnij
                  IF( iperproc(jpn) /= 0 ) THEN
                     WRITE(numout,*)'PROCESSOR',jpn-1,'compute',iperproc(jpn), 'trajectories.'
                  ENDIF
               END DO
            ENDIF
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE flo_rst

   !!=======================================================================
END MODULE florst
