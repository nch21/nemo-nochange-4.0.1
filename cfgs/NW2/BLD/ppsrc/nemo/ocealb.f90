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

MODULE ocealb
   !!======================================================================
   !!                       ***  MODULE  ocealb  ***
   !! Ocean forcing:  bulk thermohaline forcing of the ocean
   !!=====================================================================
   !! History :
   !!   NEMO     4.0  ! 2017-07  (C. Rousset) Split ocean and ice albedos
   !!----------------------------------------------------------------------
   !!   oce_alb    : albedo for ocean (clear and overcast skies)
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   oce_alb   ! routine called by sbccpl
  
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: ocealb.F90 10069 2018-08-28 14:12:24Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE oce_alb( palb_os , palb_cs )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE oce_alb  ***
      !! 
      !! ** Purpose :   Computation of the albedo of the ocean
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   palb_os   !  albedo of ocean under overcast sky
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   palb_cs   !  albedo of ocean under clear sky
      !!
      REAL(wp) ::   zcoef 
      REAL(wp) ::   rmue = 0.40    !  cosine of local solar altitude
      !!----------------------------------------------------------------------
      !
      zcoef = 0.05 / ( 1.1 * rmue**1.4 + 0.15 )   ! Parameterization of Briegled and Ramanathan, 1982
      palb_cs(:,:) = zcoef 
      palb_os(:,:) = 0.06                         ! Parameterization of Kondratyev, 1969 and Payne, 1972
      !
   END SUBROUTINE oce_alb

   !!======================================================================
END MODULE ocealb
