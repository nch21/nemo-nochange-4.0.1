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

MODULE crslbclnk
   !!======================================================================
   !!                       ***  MODULE  crslbclnk  ***
   !!               A temporary solution for lbclnk for coarsened grid.
   !! Ocean        : lateral boundary conditions for grid coarsening
   !!=====================================================================
   !! History :   ! 2012-06  (J. Simeon, G. Madec, C. Ethe, C. Calone)     Original code
   !!----------------------------------------------------------------------
   USE par_kind, ONLY: wp
   USE dom_oce
   USE crs
   !
   USE lbclnk
   USE in_out_manager
   
   INTERFACE crs_lbc_lnk
      MODULE PROCEDURE crs_lbc_lnk_3d, crs_lbc_lnk_2d
   END INTERFACE
   
   PUBLIC crs_lbc_lnk
   
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: crslbclnk.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE crs_lbc_lnk_3d( pt3d1, cd_type1, psgn, kfillmode, pfillval  )
      !!---------------------------------------------------------------------
      !!                  ***  SUBROUTINE crs_lbc_lnk  ***
      !!
      !! ** Purpose :   set lateral boundary conditions for coarsened grid
      !!
      !! ** Method  :   Swap domain indices from full to coarse domain
      !!                before arguments are passed directly to lbc_lnk.
      !!                Upon exiting, switch back to full domain indices.
      !!----------------------------------------------------------------------
      CHARACTER(len=1)                        , INTENT(in   ) ::   cd_type1 ! grid type
      REAL(wp)                                , INTENT(in   ) ::   psgn     ! control of the sign
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(inout) ::   pt3d1    ! 3D array on which the lbc is applied
      INTEGER                     , OPTIONAL  , INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = cst)
      REAL(wp)                    , OPTIONAL  , INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      !
      LOGICAL  ::   ll_grid_crs
      !!----------------------------------------------------------------------
      !
      ll_grid_crs = ( jpi == jpi_crs )
      !
      IF( .NOT.ll_grid_crs )   CALL dom_grid_crs   ! Save the parent grid information  & Switch to coarse grid domain
      !
      CALL lbc_lnk( 'crslbclnk', pt3d1, cd_type1, psgn, kfillmode, pfillval )
      !
      IF( .NOT.ll_grid_crs )   CALL dom_grid_glo   ! Return to parent grid domain
      !
   END SUBROUTINE crs_lbc_lnk_3d
   
   
   SUBROUTINE crs_lbc_lnk_2d(pt2d, cd_type, psgn, kfillmode, pfillval )
      !!---------------------------------------------------------------------
      !!                  ***  SUBROUTINE crs_lbc_lnk  ***
      !!
      !! ** Purpose :   set lateral boundary conditions for coarsened grid
      !!
      !! ** Method  :   Swap domain indices from full to coarse domain
      !!                before arguments are passed directly to lbc_lnk.
      !!                Upon exiting, switch back to full domain indices.
      !!----------------------------------------------------------------------
      CHARACTER(len=1)                    , INTENT(in   ) ::   cd_type  ! grid type
      REAL(wp)                            , INTENT(in   ) ::   psgn     ! control of the sign
      REAL(wp), DIMENSION(jpi_crs,jpj_crs), INTENT(inout) ::   pt2d     ! 2D array on which the lbc is applied
      INTEGER                 , OPTIONAL  , INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = constant)
      REAL(wp)                , OPTIONAL  , INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      !      
      LOGICAL  ::   ll_grid_crs
      !!----------------------------------------------------------------------
      !
      ll_grid_crs = ( jpi == jpi_crs )
      !
      IF( .NOT.ll_grid_crs )   CALL dom_grid_crs   ! Save the parent grid information  & Switch to coarse grid domain
      !
      CALL lbc_lnk( 'crslbclnk', pt2d, cd_type, psgn, kfillmode, pfillval )
      !
      IF( .NOT.ll_grid_crs )   CALL dom_grid_glo   ! Return to parent grid domain
      !
   END SUBROUTINE crs_lbc_lnk_2d

   !!======================================================================
END MODULE crslbclnk
