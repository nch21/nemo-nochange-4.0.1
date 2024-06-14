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

MODULE zdfsh2
   !!======================================================================
   !!                       ***  MODULE  zdfsh2  ***
   !! Ocean physics:  shear production term of TKE 
   !!=====================================================================
   !! History :   -   !  2014-10  (A. Barthelemy, G. Madec)  original code
   !!   NEMO     4.0  !  2017-04  (G. Madec)  remove u-,v-pts avm
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_sh2       : compute mixing the shear production term of TKE
   !!----------------------------------------------------------------------
   USE dom_oce        ! domain: ocean
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_sh2        ! called by zdftke, zdfglf, and zdfric
   
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: zdfsh2.F90 10069 2018-08-28 14:12:24Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE zdf_sh2( pub, pvb, pun, pvn, p_avm, p_sh2  ) 
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zdf_sh2  ***
      !!
      !! ** Purpose :   Compute the shear production term of a TKE equation
      !!
      !! ** Method  : - a stable discretization of this term is linked to the
      !!                time-space discretization of the vertical diffusion
      !!                of the OGCM. NEMO uses C-grid, a leap-frog environment 
      !!                and an implicit computation of vertical mixing term,
      !!                so the shear production at w-point is given by:
      !!                   sh2 = mi[   mi(avm) * dk[ub]/e3ub * dk[un]/e3un   ] 
      !!                       + mj[   mj(avm) * dk[vb]/e3vb * dk[vn]/e3vn   ] 
      !!                NB: wet-point only horizontal averaging of shear
      !!
      !! ** Action  : - p_sh2 shear prod. term at w-point (inner domain only)
      !!                                                   *****
      !! References :   Bruchard, OM 2002
      !! ---------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:) , INTENT(in   ) ::   pub, pvb, pun, pvn   ! before, now horizontal velocities
      REAL(wp), DIMENSION(:,:,:) , INTENT(in   ) ::   p_avm                ! vertical eddy viscosity (w-points)
      REAL(wp), DIMENSION(:,:,:) , INTENT(  out) ::   p_sh2                ! shear production of TKE (w-points)
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop arguments
      REAL(wp), DIMENSION(jpi,jpj) ::   zsh2u, zsh2v   ! 2D workspace
      !!--------------------------------------------------------------------
      !
      DO jk = 2, jpkm1
         DO jj = 1, jpjm1        !* 2 x shear production at uw- and vw-points (energy conserving form)
            DO ji = 1, jpim1
               zsh2u(ji,jj) = ( p_avm(ji+1,jj,jk) + p_avm(ji,jj,jk) ) &
                  &         * (   pun(ji,jj,jk-1) -   pun(ji,jj,jk) ) &
                  &         * (   pub(ji,jj,jk-1) -   pub(ji,jj,jk) ) / ( e3uw_n(ji,jj,jk) * e3uw_b(ji,jj,jk) ) * wumask(ji,jj,jk)
               zsh2v(ji,jj) = ( p_avm(ji,jj+1,jk) + p_avm(ji,jj,jk) ) &
                  &         * (   pvn(ji,jj,jk-1) -   pvn(ji,jj,jk) ) &
                  &         * (   pvb(ji,jj,jk-1) -   pvb(ji,jj,jk) ) / ( e3vw_n(ji,jj,jk) * e3vw_b(ji,jj,jk) ) * wvmask(ji,jj,jk)
            END DO
         END DO
         DO jj = 2, jpjm1        !* shear production at w-point
            DO ji = 2, jpim1           ! coast mask: =2 at the coast ; =1 otherwise (NB: wmask useless as zsh2 are masked)
               p_sh2(ji,jj,jk) = 0.25 * (   ( zsh2u(ji-1,jj) + zsh2u(ji,jj) ) * ( 2. - umask(ji-1,jj,jk) * umask(ji,jj,jk) )   &
                  &                       + ( zsh2v(ji,jj-1) + zsh2v(ji,jj) ) * ( 2. - vmask(ji,jj-1,jk) * vmask(ji,jj,jk) )   )
            END DO
         END DO
      END DO 
      !
   END SUBROUTINE zdf_sh2

   !!======================================================================
END MODULE zdfsh2
