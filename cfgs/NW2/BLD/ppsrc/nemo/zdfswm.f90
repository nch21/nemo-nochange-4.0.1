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

MODULE zdfswm
   !!======================================================================
   !!                       ***  MODULE  zdfswm  ***
   !! vertical physics :   surface wave-induced mixing 
   !!======================================================================
   !! History :  3.6  !  2014-10  (E. Clementi)  Original code
   !!            4.0  !  2017-04  (G. Madec)  debug + simplifications
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_swm      : update Kz due to surface wave-induced mixing 
   !!   zdf_swm_init : initilisation
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean domain variable
   USE zdf_oce        ! vertical physics: mixing coefficients
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbcwave        ! wave module
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)  
   USE lib_mpp        ! distribued memory computing library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC zdf_swm         ! routine called in zdp_phy
   PUBLIC zdf_swm_init    ! routine called in zdf_phy_init

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: zdfswm.F90 10069 2018-08-28 14:12:24Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE zdf_swm( kt, p_avm, p_avt, p_avs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE zdf_swm ***
      !!
      !! ** Purpose :Compute the swm term (qbv) to be added to
      !!             vertical viscosity and diffusivity coeffs.  
      !!
      !! ** Method  :   Compute the swm term Bv (zqb) and added it to
      !!               vertical viscosity and diffusivity coefficients
      !!                   zqb = alpha * A * Us(0) * exp (3 * k * z)
      !!               where alpha is set here to 1
      !!             
      !! ** action  :    avt, avs, avm updated by the surface wave-induced mixing
      !!                               (inner domain only)
      !!               
      !! reference : Qiao et al. GRL, 2004
      !!---------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) ::   kt             ! ocean time step
      REAL(wp), DIMENSION(:,:,:) , INTENT(inout) ::   p_avm          ! momentum Kz (w-points)
      REAL(wp), DIMENSION(:,:,:) , INTENT(inout) ::   p_avt, p_avs   ! tracer   Kz (w-points)
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      REAL(wp)::   zcoef, zqb   ! local scalar
      !!---------------------------------------------------------------------
      !
      zcoef = 1._wp * 0.353553_wp
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zqb = zcoef * hsw(ji,jj) * tsd2d(ji,jj) * EXP( -3. * wnum(ji,jj) * gdepw_n(ji,jj,jk) ) * wmask(ji,jj,jk)
               !
               p_avt(ji,jj,jk) = p_avt(ji,jj,jk) + zqb
               p_avs(ji,jj,jk) = p_avs(ji,jj,jk) + zqb
               p_avm(ji,jj,jk) = p_avm(ji,jj,jk) + zqb
            END DO
         END DO
      END DO
      !
   END SUBROUTINE zdf_swm
   
   
   SUBROUTINE zdf_swm_init
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE zdf_swm_init ***
      !!
      !! ** Purpose :   surface wave-induced mixing initialisation  
      !!
      !! ** Method  :   check the availability of surface wave fields
      !!---------------------------------------------------------------------
      !
      IF(lwp) THEN                  ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_swm_init : surface wave-driven mixing'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      IF(  .NOT.ln_wave .OR.   &
         & .NOT.ln_sdw    )   CALL ctl_stop ( 'zdf_swm_init: ln_zdfswm=T but ln_wave and ln_sdw /= T')
      !
   END SUBROUTINE zdf_swm_init

   !!======================================================================
END MODULE zdfswm
