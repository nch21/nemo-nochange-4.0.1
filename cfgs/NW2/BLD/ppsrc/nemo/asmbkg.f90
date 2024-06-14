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

MODULE asmbkg
   !!======================================================================
   !!                       ***  MODULE asmtrj -> asmbkg  ***
   !! Assimilation trajectory interface: Write to file the background state and the model state trajectory
   !!======================================================================
   !! History :       ! 2007-03 (M. Martin)  Met. Office version
   !!                 ! 2007-04 (A. Weaver)  asm_trj_wri, original code
   !!                 ! 2007-03 (K. Mogensen)  Adapt to NEMOVAR and use IOM instead of IOIPSL
   !!                 ! 2007-04 (A. Weaver)  Name change (formally asmbkg.F90). Distinguish
   !!                                        background states in Jb term and at analysis time.
   !!                                        Include state trajectory routine (currently empty)
   !!                 ! 2007-07 (A. Weaver)  Add tke_rst and flt_rst for case nitbkg=0 
   !!                 ! 2009-03 (F. Vigilant)  Add hmlp (zdfmxl) for no tracer nmldp=2 
   !!                 ! 2009-06 (F. Vigilant) asm_trj_wri: special case when kt=nit000-1
   !!                 ! 2009-07 (F. Vigilant) asm_trj_wri: add computation of eiv at restart
   !!                 ! 2010-01 (A. Vidard) split asm_trj_wri into tam_trj_wri and asm_bkg_wri
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   asm_bkg_wri  : Write out the background state
   !!   asm_trj_wri  : Write out the model state trajectory (used with 4D-Var)
   !!----------------------------------------------------------------------
   USE oce                ! Dynamics and active tracers defined in memory
   USE sbc_oce            ! Ocean surface boundary conditions
   USE zdf_oce            ! Vertical mixing variables
   USE zdfddm             ! Double diffusion mixing parameterization
   USE ldftra             ! Lateral diffusion: eddy diffusivity coefficients
   USE ldfslp             ! Lateral diffusion: slopes of neutral surfaces
   USE tradmp             ! Tracer damping
   USE zdftke             ! TKE vertical physics
   USE eosbn2             ! Equation of state (eos_bn2 routine)
   USE zdfmxl             ! Mixed layer depth
   USE dom_oce     , ONLY :   ndastp
   USE in_out_manager     ! I/O manager
   USE iom                ! I/O module
   USE asmpar             ! Parameters for the assmilation interface
   USE zdfmxl             ! mixed layer depth
   USE ice

   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   asm_bkg_wri   !: Write out the background state

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: asmbkg.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE asm_bkg_wri( kt )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE asm_bkg_wri ***
      !!
      !! ** Purpose : Write to file the background state for later use in the
      !!              inner loop of data assimilation or for direct initialization
      !!              in the outer loop.
      !!
      !! ** Method  : Write out the background state for use in the Jb term
      !!              in the cost function and for use with direct initialization
      !!              at analysis time.
      !!-----------------------------------------------------------------------
      INTEGER, INTENT( IN ) :: kt               ! Current time-step
      !
      CHARACTER (LEN=50) :: cl_asmbkg
      CHARACTER (LEN=50) :: cl_asmdin
      LOGICAL :: llok          ! Check if file exists
      INTEGER :: inum          ! File unit number
      REAL(wp) :: zdate        ! Date
      !!-----------------------------------------------------------------------

      !                                !-------------------------------------------
      IF( kt == nitbkg_r ) THEN        ! Write out background at time step nitbkg_r
         !                             !-----------------------------------========
         !
         WRITE(cl_asmbkg, FMT='(A,".nc")' ) TRIM( c_asmbkg )
         cl_asmbkg = TRIM( cl_asmbkg )
         INQUIRE( FILE = cl_asmbkg, EXIST = llok )
         !
         IF( .NOT. llok ) THEN
            IF(lwp) WRITE(numout,*) ' Setting up assimilation background file '// TRIM( c_asmbkg )
            !
            !                                      ! Define the output file        
            CALL iom_open( c_asmbkg, inum, ldwrt = .TRUE. )
            !
            IF( nitbkg_r == nit000 - 1 ) THEN      ! Treat special case when nitbkg = 0
               zdate = REAL( ndastp )
               IF( ln_zdftke ) THEN                   ! read turbulent kinetic energy ( en )
                  IF(lwp) WRITE(numout,*) ' Reading TKE (en) from restart...'
                  CALL tke_rst( nit000, 'READ' )
               ENDIF
            ELSE
               zdate = REAL( ndastp )
            ENDIF
            !
            !                                      ! Write the information
            CALL iom_rstput( kt, nitbkg_r, inum, 'rdastp' , zdate             )
            CALL iom_rstput( kt, nitbkg_r, inum, 'un'     , un                )
            CALL iom_rstput( kt, nitbkg_r, inum, 'vn'     , vn                )
            CALL iom_rstput( kt, nitbkg_r, inum, 'tn'     , tsn(:,:,:,jp_tem) )
            CALL iom_rstput( kt, nitbkg_r, inum, 'sn'     , tsn(:,:,:,jp_sal) )
            CALL iom_rstput( kt, nitbkg_r, inum, 'sshn'   , sshn              )
            IF( ln_zdftke )   CALL iom_rstput( kt, nitbkg_r, inum, 'en'     , en                )
            !
            CALL iom_close( inum )
         ENDIF
         !
      ENDIF

      !                                !-------------------------------------------
      IF( kt == nitdin_r ) THEN        ! Write out background at time step nitdin_r
         !                             !-----------------------------------========
         !
         WRITE(cl_asmdin, FMT='(A,".nc")' ) TRIM( c_asmdin )
         cl_asmdin = TRIM( cl_asmdin )
         INQUIRE( FILE = cl_asmdin, EXIST = llok )
         !
         IF( .NOT. llok ) THEN
            IF(lwp) WRITE(numout,*) ' Setting up assimilation background file '// TRIM( c_asmdin )
            !
            !                                      ! Define the output file        
            CALL iom_open( c_asmdin, inum, ldwrt = .TRUE. )
            !
            IF( nitdin_r == nit000 - 1 ) THEN      ! Treat special case when nitbkg = 0

               zdate = REAL( ndastp )
            ELSE
               zdate = REAL( ndastp )
            ENDIF
            !
            !                                      ! Write the information
            CALL iom_rstput( kt, nitdin_r, inum, 'rdastp' , zdate             )
            CALL iom_rstput( kt, nitdin_r, inum, 'un'     , un                )
            CALL iom_rstput( kt, nitdin_r, inum, 'vn'     , vn                )
            CALL iom_rstput( kt, nitdin_r, inum, 'tn'     , tsn(:,:,:,jp_tem) )
            CALL iom_rstput( kt, nitdin_r, inum, 'sn'     , tsn(:,:,:,jp_sal) )
            CALL iom_rstput( kt, nitdin_r, inum, 'sshn'   , sshn              )
            IF( nn_ice == 2 ) THEN
	            IF( ALLOCATED(at_i) ) THEN
                  CALL iom_rstput( kt, nitdin_r, inum, 'iceconc', at_i(:,:)   )
               ELSE
		            CALL ctl_warn('asm_bkg_wri: Ice concentration not written to background ',   &
		               &          'as ice variable at_i not allocated on this timestep')
	            ENDIF
            ENDIF
            !
            CALL iom_close( inum )
         ENDIF
         !
      ENDIF
      !                    
   END SUBROUTINE asm_bkg_wri

   !!======================================================================
END MODULE asmbkg
