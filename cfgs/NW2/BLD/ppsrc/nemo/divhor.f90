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

MODULE divhor
   !!==============================================================================
   !!                       ***  MODULE  divhor  ***
   !! Ocean diagnostic variable : now horizontal divergence
   !!==============================================================================
   !! History :  1.0  ! 2002-09  (G. Madec, E. Durand)  Free form, F90
   !!             -   ! 2005-01  (J. Chanut) Unstructured open boundaries
   !!             -   ! 2003-08  (G. Madec)  merged of cur and div, free form, F90
   !!             -   ! 2005-01  (J. Chanut, A. Sellar) unstructured open boundaries
   !!            3.3  ! 2010-09  (D.Storkey and E.O'Dea) bug fixes for BDY module
   !!             -   ! 2010-10  (R. Furner, G. Madec) runoff and cla added directly here
   !!            3.7  ! 2014-01  (G. Madec) suppression of velocity curl from in-core memory
   !!             -   ! 2014-12  (G. Madec) suppression of cross land advection option
   !!             -   ! 2015-10  (G. Madec) add velocity and rnf flag in argument of div_hor
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   div_hor    : Compute the horizontal divergence field
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce, ONLY : ln_rnf, ln_isf ! surface boundary condition: ocean
   USE sbcrnf          ! river runoff 
   USE sbcisf          ! ice shelf
   USE iscplhsb        ! ice sheet / ocean coupling
   USE iscplini        ! ice sheet / ocean coupling
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   div_hor    ! routine called by step.F90 and istate.F90

   !! * Substitutions
   !!----------------------------------------------------------------------
   !!                   ***  vectopt_loop_substitute  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute the inner loop start/end indices with CPP macro
   !!                allow unrolling of do-loop (useful with vector processors)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: vectopt_loop_substitute.h90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: divhor.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE div_hor( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE div_hor  ***
      !!                    
      !! ** Purpose :   compute the horizontal divergence at now time-step
      !!
      !! ** Method  :   the now divergence is computed as :
      !!         hdivn = 1/(e1e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
      !!      and correct with runoff inflow (div_rnf) and cross land flow (div_cla) 
      !!
      !! ** Action  : - update hdivn, the now horizontal divergence
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk    ! dummy loop indices
      REAL(wp) ::   zraur, zdep   ! local scalars
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('div_hor')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'div_hor : horizontal velocity divergence '
         IF(lwp) WRITE(numout,*) '~~~~~~~   '
      ENDIF
      !
      DO jk = 1, jpkm1                                      !==  Horizontal divergence  ==!
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               hdivn(ji,jj,jk) = (  e2u(ji  ,jj) * e3u_n(ji  ,jj,jk) * un(ji  ,jj,jk)      &
                  &               - e2u(ji-1,jj) * e3u_n(ji-1,jj,jk) * un(ji-1,jj,jk)      &
                  &               + e1v(ji,jj  ) * e3v_n(ji,jj  ,jk) * vn(ji,jj  ,jk)      &
                  &               - e1v(ji,jj-1) * e3v_n(ji,jj-1,jk) * vn(ji,jj-1,jk)  )   &
                  &            * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
            END DO  
         END DO  
      END DO
      !
      IF( ln_rnf )   CALL sbc_rnf_div( hdivn )              !==  runoffs    ==!   (update hdivn field)
      !
      IF( ln_isf )   CALL sbc_isf_div( hdivn )      !==  ice shelf  ==!   (update hdivn field)
      !
      IF( ln_iscpl .AND. ln_hsb )   CALL iscpl_div( hdivn ) !==  ice sheet  ==!   (update hdivn field)
      !
      CALL lbc_lnk( 'divhor', hdivn, 'T', 1. )   !   (no sign change)
      !
      IF( ln_timing )   CALL timing_stop('div_hor')
      !
   END SUBROUTINE div_hor
   
   !!======================================================================
END MODULE divhor
