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

MODULE usrdef_fmask
   !!======================================================================
   !!                     ***  MODULE usrdef_fmask   ***
   !!
   !!                      ===  ORCA configuration  ===
   !!                            (2 and 1 degrees)
   !!
   !! User defined : alteration of land/sea f-point mask in some straits
   !!======================================================================
   !! History :  4.0  ! 2016-06  (G. Madec, S. Flavoni)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_fmask  : alteration of f-point land/ocean mask in some straits
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! Massively Parallel Processing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_fmask    ! routine called by dommsk.F90
   PUBLIC   usr_def_fmask_scalar

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
   !! $Id: usrdef_fmask.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_fmask(  cd_cfg, kcfg, pfmsk )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk  ***
      !!
      !! ** Purpose :   User defined alteration of the lateral boundary 
      !!              condition on velocity.
      !!
      !! ** Method  :   Local change of the value of fmask at lateral ocean/land 
      !!              boundary in straits in order to increase the viscous 
      !!              boundary layer and thus reduce the transport through the 
      !!              corresponding straits.
      !!                Here only alterations in ORCA R2 and R1 cases
      !!
      !! ** Action :   fmask : land/ocean mask at f-point with increased value 
      !!                       in some user defined straits
      !!----------------------------------------------------------------------
      CHARACTER(len=*)          , INTENT(in   ) ::   cd_cfg   ! configuration name
      INTEGER                   , INTENT(in   ) ::   kcfg     ! configuration identifier 
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pfmsk    ! Ocean/Land f-point mask including lateral boundary cond.
      !!external
	!mod par_oce

      !
      INTEGER  ::   iif, iil, ii0, ii1, ii   ! local integers
      INTEGER  ::   ijf, ijl, ij0, ij1       !   -       -
      INTEGER  ::   isrow                    ! index for ORCA1 starting row
	INTEGER  ::   ji,jj,jk
      !!----------------------------------------------------------------------
      !
      IF( TRIM( cd_cfg ) == "orca" ) THEN      !==  ORCA Configurations  ==!
         !
         SELECT CASE ( kcfg )
         !
         CASE( 2 )                           !  R2 case
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_fmask : ORCA_R2: increase lateral friction near the following straits:'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
            !
            IF(lwp) WRITE(numout,*) '      Gibraltar '
            ij0 = 101   ;   ij1 = 101           ! Gibraltar strait  : partial slip (pfmsk=0.5)
            ii0 = 139   ;   ii1 = 140   ;   pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  0.5_wp
            ij0 = 102   ;   ij1 = 102
            ii0 = 139   ;   ii1 = 140   ;   pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  0.5_wp
            !
            IF(lwp) WRITE(numout,*) '      Bab el Mandeb '
            ij0 =  87   ;   ij1 =  88           ! Bab el Mandeb : partial slip (pfmsk=1)
            ii0 = 160   ;   ii1 = 160   ;   pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  1._wp
            ij0 =  88   ;   ij1 =  88
            ii0 = 159   ;   ii1 = 159   ;   pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  1._wp
            !
            ! We keep this as an example but it is instable in this case 
            !IF(lwp) WRITE(numout,*) '      Danish straits '
            !         ij0 = 115   ;   ij1 = 115 ! Danish straits  : strong slip (pfmsk > 2)
            !         ii0 = 145   ;   ii1 = 146   ;   pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 4._wp
            !         ij0 = 116   ;   ij1 = 116
            !         ii0 = 145   ;   ii1 = 146   ;   pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 4._wp
            !
         CASE( 1 )                           ! R1 case
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_fmask : ORCA_R1: increase lateral friction near the following straits:'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'   
!!gm    ! This dirty section will be suppressed by simplification process:
!!gm    ! all this will come back in input files
!!gm    ! Currently these hard-wired indices relate to configuration with extend grid (jpjglo=332)
            !
            isrow = 332 - jpjglo
            !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   orca_r1: increase friction near the following straits : '
            IF(lwp) WRITE(numout,*) '      Gibraltar '
            ii0 = 282           ;   ii1 = 283        ! Gibraltar Strait 
            ij0 = 241 - isrow   ;   ij1 = 241 - isrow   ;   pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  
            !
            IF(lwp) WRITE(numout,*) '      Bhosporus '
            ii0 = 314           ;   ii1 = 315        ! Bhosporus Strait 
            ij0 = 248 - isrow   ;   ij1 = 248 - isrow   ;   pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  
            !
            IF(lwp) WRITE(numout,*) '      Makassar (Top) '
            ii0 =  48           ;   ii1 =  48        ! Makassar Strait (Top) 
            ij0 = 189 - isrow   ;   ij1 = 190 - isrow   ;   pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 3._wp  
            !
            IF(lwp) WRITE(numout,*) '      Lombok '
            ii0 =  44           ;   ii1 =  44        ! Lombok Strait 
            ij0 = 164 - isrow   ;   ij1 = 165 - isrow   ;   pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  
            !
            IF(lwp) WRITE(numout,*) '      Ombai '
            ii0 =  53           ;   ii1 =  53        ! Ombai Strait 
            ij0 = 164 - isrow   ;   ij1 = 165 - isrow   ;   pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  
            !
            IF(lwp) WRITE(numout,*) '      Timor Passage '
            ii0 =  56           ;   ii1 =  56        ! Timor Passage 
            ij0 = 164 - isrow   ;   ij1 = 165 - isrow   ;   pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  
            !
            IF(lwp) WRITE(numout,*) '      West Halmahera '
            ii0 =  58           ;   ii1 =  58        ! West Halmahera Strait 
            ij0 = 181 - isrow   ;   ij1 = 182 - isrow   ;   pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 3._wp  
            !
            IF(lwp) WRITE(numout,*) '      East Halmahera '
            ii0 =  55           ;   ii1 =  55        ! East Halmahera Strait 
            ij0 = 181 - isrow   ;   ij1 = 182 - isrow   ;   pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 3._wp  
            !
         CASE DEFAULT
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_fmask : ORCA_R', kcfg,' : NO alteration of fmask in specific straits '
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'   
         END SELECT
      ELSE
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'usr_def_fmask : NO alteration of fmask in specific straits '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      !CALL lbc_lnk( 'usrdef_fmask', pfmsk, 'F', 1._wp )      ! Lateral boundary conditions on fmask
      !
   END SUBROUTINE usr_def_fmask
   
   SUBROUTINE usr_def_fmask_scalar(cd_cfg, kcfg)
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk  ***
      !!
      !! ** Purpose :   User defined alteration of the lateral boundary 
      !!              condition on velocity.
      !!
      !! ** Method  :   Local change of the value of fmask at lateral ocean/land 
      !!              boundary in straits in order to increase the viscous 
      !!              boundary layer and thus reduce the transport through the 
      !!              corresponding straits.
      !!                Here only alterations in ORCA R2 and R1 cases
      !!
      !! ** Action :   fmask : land/ocean mask at f-point with increased value 
      !!                       in some user defined straits
      !!----------------------------------------------------------------------
      CHARACTER(len=*)          , INTENT(in   ) ::   cd_cfg   ! configuration name
      INTEGER                   , INTENT(in   ) ::   kcfg     ! configuration identifier 
      !
      INTEGER  ::   iif, iil, ii0, ii1, ii   ! local integers
      INTEGER  ::   ijf, ijl, ij0, ij1       !   -       -
      INTEGER  ::   isrow                    ! index for ORCA1 starting row
	INTEGER  ::   ji,jj,jk
      !!----------------------------------------------------------------------
      !
      IF( TRIM( cd_cfg ) == "orca" ) THEN      !==  ORCA Configurations  ==!
         !
         SELECT CASE ( kcfg )
         !
         CASE( 2 )                           !  R2 case
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_fmask : ORCA_R2: increase lateral friction near the following straits:'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
            !
            IF(lwp) WRITE(numout,*) '      Gibraltar '
            ij0 = 101   ;   ij1 = 101           ! Gibraltar strait  : partial slip (pfmsk=0.5)
            ii0 = 139   ;   ii1 = 140
            ij0 = 102   ;   ij1 = 102
            ii0 = 139   ;   ii1 = 140
            !
            IF(lwp) WRITE(numout,*) '      Bab el Mandeb '
            ij0 =  87   ;   ij1 =  88           ! Bab el Mandeb : partial slip (pfmsk=1)
            ii0 = 160   ;   ii1 = 160
            ij0 =  88   ;   ij1 =  88
            ii0 = 159   ;   ii1 = 159
            !
            ! We keep this as an example but it is instable in this case 
            !IF(lwp) WRITE(numout,*) '      Danish straits '
            !         ij0 = 115   ;   ij1 = 115 ! Danish straits  : strong slip (pfmsk > 2)
            !         ii0 = 145   ;   ii1 = 146   ;   pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 4._wp
            !         ij0 = 116   ;   ij1 = 116
            !         ii0 = 145   ;   ii1 = 146   ;   pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 4._wp
            !
         CASE( 1 )                           ! R1 case
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_fmask : ORCA_R1: increase lateral friction near the following straits:'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'   
!!gm    ! This dirty section will be suppressed by simplification process:
!!gm    ! all this will come back in input files
!!gm    ! Currently these hard-wired indices relate to configuration with extend grid (jpjglo=332)
            !
            isrow = 332 - jpjglo
            !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   orca_r1: increase friction near the following straits : '
            IF(lwp) WRITE(numout,*) '      Gibraltar '
            ii0 = 282           ;   ii1 = 283        ! Gibraltar Strait 
            ij0 = 241 - isrow   ;   ij1 = 241 - isrow
            !
            IF(lwp) WRITE(numout,*) '      Bhosporus '
            ii0 = 314           ;   ii1 = 315        ! Bhosporus Strait 
            ij0 = 248 - isrow   ;   ij1 = 248 - isrow
            !
            IF(lwp) WRITE(numout,*) '      Makassar (Top) '
            ii0 =  48           ;   ii1 =  48        ! Makassar Strait (Top) 
            ij0 = 189 - isrow   ;   ij1 = 190 - isrow
            !
            IF(lwp) WRITE(numout,*) '      Lombok '
            ii0 =  44           ;   ii1 =  44        ! Lombok Strait 
            ij0 = 164 - isrow   ;   ij1 = 165 - isrow
            !
            IF(lwp) WRITE(numout,*) '      Ombai '
            ii0 =  53           ;   ii1 =  53        ! Ombai Strait 
            ij0 = 164 - isrow   ;   ij1 = 165 - isrow
            !
            IF(lwp) WRITE(numout,*) '      Timor Passage '
            ii0 =  56           ;   ii1 =  56        ! Timor Passage 
            ij0 = 164 - isrow   ;   ij1 = 165 - isrow
            !
            IF(lwp) WRITE(numout,*) '      West Halmahera '
            ii0 =  58           ;   ii1 =  58        ! West Halmahera Strait 
            ij0 = 181 - isrow   ;   ij1 = 182 - isrow
            !
            IF(lwp) WRITE(numout,*) '      East Halmahera '
            ii0 =  55           ;   ii1 =  55        ! East Halmahera Strait 
            ij0 = 181 - isrow   ;   ij1 = 182 - isrow
         CASE DEFAULT 
         END SELECT
      ELSE
      ENDIF
      !
      !CALL lbc_lnk( 'usrdef_fmask', pfmsk, 'F', 1._wp )      ! Lateral boundary conditions on fmask
      !
   END SUBROUTINE usr_def_fmask_scalar

   !!======================================================================
END MODULE usrdef_fmask
