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

MODULE icethd_pnd 
   !!======================================================================
   !!                     ***  MODULE  icethd_pnd   ***
   !!   sea-ice: Melt ponds on top of sea ice
   !!======================================================================
   !! history :       !  2012     (O. Lecomte)       Adaptation from Flocco and Turner
   !!                 !  2017     (M. Vancoppenolle, O. Lecomte, C. Rousset) Implementation
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_si3' :                                     SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_thd_pnd_init : some initialization and namelist read
   !!   ice_thd_pnd      : main calling routine
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain
   USE ice            ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamics variables
   USE icetab         ! sea-ice: 1D <==> 2D transformation
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_pnd_init    ! routine called by icestp.F90
   PUBLIC   ice_thd_pnd         ! routine called by icestp.F90

   INTEGER ::              nice_pnd    ! choice of the type of pond scheme
   !                                   ! associated indices:
   INTEGER, PARAMETER ::   np_pndNO  = 0   ! No pond scheme
   INTEGER, PARAMETER ::   np_pndCST = 1   ! Constant pond scheme
   INTEGER, PARAMETER ::   np_pndH12 = 2   ! Evolutive pond scheme (Holland et al. 2012)

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
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icethd_pnd.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_pnd
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE ice_thd_pnd   ***
      !!               
      !! ** Purpose :   change melt pond fraction
      !!                
      !! ** Method  :   brut force
      !!-------------------------------------------------------------------
      !
      SELECT CASE ( nice_pnd )
      !
      CASE (np_pndCST)   ;   CALL pnd_CST    !==  Constant melt ponds  ==!
         !
      CASE (np_pndH12)   ;   CALL pnd_H12    !==  Holland et al 2012 melt ponds  ==!
         !
      END SELECT
      !
   END SUBROUTINE ice_thd_pnd 


   SUBROUTINE pnd_CST 
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE pnd_CST  ***
      !!
      !! ** Purpose :   Compute melt pond evolution
      !!
      !! ** Method  :   Melt pond fraction and thickness are prescribed 
      !!                to non-zero values when t_su = 0C
      !!
      !! ** Tunable parameters : pond fraction (rn_apnd), pond depth (rn_hpnd)
      !!                
      !! ** Note   : Coupling with such melt ponds is only radiative
      !!             Advection, ridging, rafting... are bypassed
      !!
      !! ** References : Bush, G.W., and Trump, D.J. (2017)
      !!-------------------------------------------------------------------
      INTEGER  ::   ji        ! loop indices
      !!-------------------------------------------------------------------
      DO ji = 1, npti
         !
         IF( a_i_1d(ji) > 0._wp .AND. t_su_1d(ji) >= rt0 ) THEN
            a_ip_frac_1d(ji) = rn_apnd
            h_ip_1d(ji)      = rn_hpnd    
            a_ip_1d(ji)      = a_ip_frac_1d(ji) * a_i_1d(ji)
         ELSE
            a_ip_frac_1d(ji) = 0._wp
            h_ip_1d(ji)      = 0._wp    
            a_ip_1d(ji)      = 0._wp
         ENDIF
         !
      END DO
      !
   END SUBROUTINE pnd_CST


   SUBROUTINE pnd_H12
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE pnd_H12  ***
      !!
      !! ** Purpose    : Compute melt pond evolution
      !!
      !! ** Method     : Empirical method. A fraction of meltwater is accumulated in ponds 
      !!                 and sent to ocean when surface is freezing
      !!
      !!                 pond growth:      Vp = Vp + dVmelt
      !!                    with dVmelt = R/rhow * ( rhoi*dh_i + rhos*dh_s ) * a_i
      !!                 pond contraction: Vp = Vp * exp(0.01*MAX(Tp-Tsu,0)/Tp)
      !!                    with Tp = -2degC
      !!  
      !! ** Tunable parameters : (no real expertise yet, ideas?)
      !! 
      !! ** Note       : Stolen from CICE for quick test of the melt pond
      !!                 radiation and freshwater interfaces
      !!                 Coupling can be radiative AND freshwater
      !!                 Advection, ridging, rafting are called
      !!
      !! ** References : Holland, M. M. et al (J Clim 2012)
      !!-------------------------------------------------------------------
      REAL(wp), PARAMETER ::   zrmin       = 0.15_wp  ! minimum fraction of available meltwater retained for melt ponding
      REAL(wp), PARAMETER ::   zrmax       = 0.70_wp  ! maximum     -           -         -         -            -
      REAL(wp), PARAMETER ::   zpnd_aspect = 0.8_wp   ! pond aspect ratio
      REAL(wp), PARAMETER ::   zTp         = -2._wp   ! reference temperature
      !
      REAL(wp) ::   zfr_mlt          ! fraction of available meltwater retained for melt ponding
      REAL(wp) ::   zdv_mlt          ! available meltwater for melt ponding
      REAL(wp) ::   z1_Tp            ! inverse reference temperature
      REAL(wp) ::   z1_rhow          ! inverse freshwater density
      REAL(wp) ::   z1_zpnd_aspect   ! inverse pond aspect ratio
      REAL(wp) ::   zfac, zdum
      !
      INTEGER  ::   ji   ! loop indices
      !!-------------------------------------------------------------------
      z1_rhow        = 1._wp / rhow 
      z1_zpnd_aspect = 1._wp / zpnd_aspect
      z1_Tp          = 1._wp / zTp 

      DO ji = 1, npti
         !                                                        !--------------------------------!
         IF( h_i_1d(ji) < rn_himin) THEN                          ! Case ice thickness < rn_himin  !
            !                                                     !--------------------------------!
            !--- Remove ponds on thin ice
            a_ip_1d(ji)      = 0._wp
            a_ip_frac_1d(ji) = 0._wp
            h_ip_1d(ji)      = 0._wp
            !                                                     !--------------------------------!
         ELSE                                                     ! Case ice thickness >= rn_himin !
            !                                                     !--------------------------------!
            v_ip_1d(ji) = h_ip_1d(ji) * a_ip_1d(ji)   ! record pond volume at previous time step
            !
            ! available meltwater for melt ponding [m, >0] and fraction
            zdv_mlt = -( dh_i_sum(ji)*rhoi + dh_s_mlt(ji)*rhos ) * z1_rhow * a_i_1d(ji)
            zfr_mlt = zrmin + ( zrmax - zrmin ) * a_i_1d(ji)  ! from CICE doc
            !zfr_mlt = zrmin + zrmax * a_i_1d(ji)             ! from Holland paper 
            !
            !--- Pond gowth ---!
            ! v_ip should never be negative, otherwise code crashes
            v_ip_1d(ji) = MAX( 0._wp, v_ip_1d(ji) + zfr_mlt * zdv_mlt )
            !
            ! melt pond mass flux (<0)
            IF( zdv_mlt > 0._wp ) THEN
               zfac = zfr_mlt * zdv_mlt * rhow * r1_rdtice
               wfx_pnd_1d(ji) = wfx_pnd_1d(ji) - zfac
               !
               ! adjust ice/snow melting flux to balance melt pond flux (>0)
               zdum = zfac / ( wfx_snw_sum_1d(ji) + wfx_sum_1d(ji) )
               wfx_snw_sum_1d(ji) = wfx_snw_sum_1d(ji) * (1._wp + zdum)
               wfx_sum_1d(ji)     = wfx_sum_1d(ji)     * (1._wp + zdum)
            ENDIF
            !
            !--- Pond contraction (due to refreezing) ---!
            v_ip_1d(ji) = v_ip_1d(ji) * EXP( 0.01_wp * MAX( zTp+rt0 - t_su_1d(ji), 0._wp ) * z1_Tp )
            !
            ! Set new pond area and depth assuming linear relation between h_ip and a_ip_frac
            !    h_ip = zpnd_aspect * a_ip_frac = zpnd_aspect * a_ip/a_i
            a_ip_1d(ji)      = SQRT( v_ip_1d(ji) * z1_zpnd_aspect * a_i_1d(ji) )
            a_ip_frac_1d(ji) = a_ip_1d(ji) / a_i_1d(ji)
            h_ip_1d(ji)      = zpnd_aspect * a_ip_frac_1d(ji)
            !
         ENDIF
      END DO
      !
   END SUBROUTINE pnd_H12


   SUBROUTINE ice_thd_pnd_init 
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_thd_pnd_init   ***
      !!
      !! ** Purpose : Physical constants and parameters linked to melt ponds
      !!              over sea ice
      !!
      !! ** Method  :  Read the namthd_pnd  namelist and check the melt pond  
      !!               parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namthd_pnd  
      !!-------------------------------------------------------------------
      INTEGER  ::   ios, ioptio   ! Local integer
      !!
      NAMELIST/namthd_pnd/  ln_pnd, ln_pnd_H12, ln_pnd_CST, rn_apnd, rn_hpnd, ln_pnd_alb
      !!-------------------------------------------------------------------
      !
      REWIND( numnam_ice_ref )              ! Namelist namthd_pnd  in reference namelist : Melt Ponds  
      READ  ( numnam_ice_ref, namthd_pnd, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namthd_pnd  in reference namelist' )
      REWIND( numnam_ice_cfg )              ! Namelist namthd_pnd  in configuration namelist : Melt Ponds
      READ  ( numnam_ice_cfg, namthd_pnd, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namthd_pnd in configuration namelist' )
      IF(lwm) WRITE ( numoni, namthd_pnd )
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_thd_pnd_init: ice parameters for melt ponds'
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namicethd_pnd:'
         WRITE(numout,*) '      Melt ponds activated or not                                     ln_pnd     = ', ln_pnd
         WRITE(numout,*) '         Evolutive  melt pond fraction and depth (Holland et al 2012) ln_pnd_H12 = ', ln_pnd_H12
         WRITE(numout,*) '         Prescribed melt pond fraction and depth                      ln_pnd_CST = ', ln_pnd_CST
         WRITE(numout,*) '            Prescribed pond fraction                                  rn_apnd    = ', rn_apnd
         WRITE(numout,*) '            Prescribed pond depth                                     rn_hpnd    = ', rn_hpnd
         WRITE(numout,*) '         Melt ponds affect albedo or not                              ln_pnd_alb = ', ln_pnd_alb
      ENDIF
      !
      !                             !== set the choice of ice pond scheme ==!
      ioptio = 0
      IF( .NOT.ln_pnd ) THEN   ;   ioptio = ioptio + 1   ;   nice_pnd = np_pndNO     ;   ENDIF
      IF( ln_pnd_CST  ) THEN   ;   ioptio = ioptio + 1   ;   nice_pnd = np_pndCST    ;   ENDIF
      IF( ln_pnd_H12  ) THEN   ;   ioptio = ioptio + 1   ;   nice_pnd = np_pndH12    ;   ENDIF
      IF( ioptio /= 1 )   &
         & CALL ctl_stop( 'ice_thd_pnd_init: choose either none (ln_pnd=F) or only one pond scheme (ln_pnd_H12 or ln_pnd_CST)' )
      !
      SELECT CASE( nice_pnd )
      CASE( np_pndNO )         
         IF( ln_pnd_alb ) THEN ; ln_pnd_alb = .FALSE. ; CALL ctl_warn( 'ln_pnd_alb=false when no ponds' ) ; ENDIF
      END SELECT
      !
   END SUBROUTINE ice_thd_pnd_init
   

   !!======================================================================
END MODULE icethd_pnd 
