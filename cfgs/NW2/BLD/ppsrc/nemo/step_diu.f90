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

MODULE step_diu
   !!======================================================================
   !!                       ***  MODULE stp_diu  ***
   !! Time-stepping of diurnal cycle models
   !!======================================================================
   !! History :  3.7  ! 2015-11  (J. While)  Original code

   USE diurnal_bulk    ! diurnal SST bulk routines  (diurnal_sst_takaya routine) 
   USE cool_skin       ! diurnal cool skin correction (diurnal_sst_coolskin routine)   
   USE iom
   USE sbc_oce
   USE sbcmod           ! surface boundary condition       (sbc     routine)
   USE diaobs           ! Observation operator
   USE oce
   USE daymod
   USE restart          ! ocean restart                    (rst_wri routine)
   USE timing           ! Timing
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp_diurnal   ! called by nemogcm.F90 or step.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: step_diu.F90 10069 2018-08-28 14:12:24Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

   CONTAINS

   SUBROUTINE stp_diurnal( kstp ) 
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index 
      !!---------------------------------------------------------------------- 
      !!                     ***  ROUTINE stp_diurnal  *** 
      !!                       
      !! ** Purpose : - Time stepping of diurnal SST model only 
      !!   
      !! ** Method  : -1- Update forcings and data   
      !!              -2- Update ocean physics   
      !!              -3- Compute the t and s trends   
      !!              -4- Update t and s   
      !!              -5- Compute the momentum trends 
      !!              -6- Update the horizontal velocity 
      !!              -7- Compute the diagnostics variables (rd,N2, div,cur,w) 
      !!              -8- Outputs and diagnostics 
      !!---------------------------------------------------------------------- 
      INTEGER ::   jk       ! dummy loop indices
      INTEGER ::   indic    ! error indicator if < 0 
      REAL(wp), DIMENSION(jpi,jpj) :: z_fvel_bkginc, z_hflux_bkginc     
      !! --------------------------------------------------------------------- 
      
      IF(ln_diurnal_only) THEN
         indic = 0                                 ! reset to no error condition 
         IF( kstp /= nit000 )   CALL day( kstp )   ! Calendar (day was already called at nit000 in day_init) 
 
         CALL iom_setkt( kstp - nit000 + 1, cxios_context )   ! tell iom we are at time step kstp
         IF( ln_crs ) THEN
            CALL iom_setkt( kstp - nit000 + 1, TRIM(cxios_context)//"_crs" ) ! tell iom we are at time step kstp
         ENDIF
       
            CALL sbc    ( kstp )                      ! Sea Boundary Conditions 
      ENDIF
     
      ! Cool skin
      IF( .NOT.ln_diurnal )   CALL ctl_stop( "stp_diurnal: ln_diurnal not set" )
         
      IF( .NOT. ln_blk    )   CALL ctl_stop( "stp_diurnal: diurnal flux processing only implemented for bulk forcing" ) 

      CALL diurnal_sst_coolskin_step( qns, taum, rhop(:,:,1), rdt)

      CALL iom_put( "sst_wl"   , x_dsst               )    ! warm layer (write out before update below).
      CALL iom_put( "sst_cs"   , x_csdsst             )    ! cool skin

      ! Diurnal warm layer model       
      CALL diurnal_sst_takaya_step( kstp, & 
      &    qsr, qns, taum, rhop(:,:,1), rdt) 

      IF( ln_diurnal_only ) THEN
         IF( ln_diaobs )         CALL dia_obs( kstp )         ! obs-minus-model (assimilation) diagnostics (call after dynamics update)
     
         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
         ! Control and restarts 
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
         IF( kstp == nit000   )   CALL iom_close( numror )     ! close input  ocean restart file 
         IF( lrst_oce         )   CALL rst_write    ( kstp )   ! write output ocean restart file
     
         IF( ln_timing .AND.  kstp == nit000  )   CALL timing_reset 
      ENDIF
       
   END SUBROUTINE stp_diurnal  
   
END MODULE step_diu
