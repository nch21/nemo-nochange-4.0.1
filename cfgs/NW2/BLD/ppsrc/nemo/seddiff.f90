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

MODULE seddiff
   !!======================================================================
   !!              ***  MODULE  seddsr  ***
   !!    Sediment : dissolution and reaction in pore water related 
   !!    related to organic matter
   !!=====================================================================
   !! * Modules used
   USE sed     ! sediment global variable
   USE sed_oce
   USE sedmat  ! linear system of equations
   USE sedini
   USE lib_mpp         ! distribued memory computing library
   USE lib_fortran

   IMPLICIT NONE
   PRIVATE

   PUBLIC sed_diff

   !! * Module variables

   !! $Id: seddsr.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS
   
   SUBROUTINE sed_diff( kt, knt ) 
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_diff  ***
      !! 
      !!  ** Purpose :  computes pore water diffusion
      !!
      !!  ** Methode :  implicit computation of undersaturation
      !!               resulting from diffusive pore water transport.
      !!
      !!  ** Remarks :
      !!              - undersaturation : deviation from saturation concentration
      !!   History :
      !!        !  98-08 (E. Maier-Reimer, Christoph Heinze )  Original code
      !!        !  04-10 (N. Emprin, M. Gehlen ) f90
      !!        !  06-04 (C. Ethe)  Re-organization
      !!        !  19-08 (O. Aumont) Debugging and improvement of the model
      !!----------------------------------------------------------------------
      !! Arguments
      INTEGER, INTENT(in) ::   kt, knt       ! number of iteration
      ! --- local variables
      INTEGER :: ji, jk, js   ! dummy looop indices

      REAL(wp), DIMENSION(jpoce,jpksed) :: zrearat1, zrearat2   ! reaction rate in pore water
      !!
      !!----------------------------------------------------------------------

      IF( ln_timing )  CALL timing_start('sed_diff')
!
      IF( kt == nitsed000 .AND. knt == 1 ) THEN
         IF (lwp) THEN
            WRITE(numsed,*) ' sed_diff : pore-water diffusion '
            WRITE(numsed,*) ' '
         ENDIF
      ENDIF

     ! Initializations
     !----------------------
      zrearat1(:,:)   = 0.
      zrearat2(:,:) = 0.

      !---------------------------
      ! Solves PO4 diffusion 
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwpo4, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwpo4), dtsed2 / 2.0 )

      !---------------------------
      ! Solves NH4 diffusion 
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwnh4, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwnh4), dtsed2 / 2.0 )

      !---------------------------
      ! Solves Fe2+ diffusion 
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwfe2, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwfe2), dtsed2 / 2.0 )

      !---------------------------
      ! Solves H2S diffusion 
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwh2s, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwh2s), dtsed2 / 2.0  )

      !---------------------------
      ! Solves SO4 diffusion 
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwso4, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwso4), dtsed2 / 2.0 )

      !---------------------------
      ! Solves O2 diffusion
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwoxy, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwoxy), dtsed2 / 2.0 )

      !---------------------------
      ! Solves NO3 diffusion
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwno3, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwno3), dtsed2 / 2.0 )

      CALL sed_mat( jwdic, jpoce, jpksed, zrearat1, zrearat2, sedligand(:,:), dtsed2 / 2.0 )

      IF( ln_timing )  CALL timing_stop('sed_diff')
!      
   END SUBROUTINE sed_diff

END MODULE seddiff
