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

MODULE usrdef_nam
   !!======================================================================
   !!                     ***  MODULE usrdef_nam   ***
   !!
   !!                     ===  GYRE configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-03  (S. Flavoni, G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!   usr_def_hgr   : initialize the horizontal mesh 
   !!----------------------------------------------------------------------
   USE dom_oce  , ONLY: nimpp, njmpp       ! ocean space and time domain
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam   ! called in nemogcm.F90 module

   !                              !!* namusr_def namelist *!!
   LOGICAL, PUBLIC ::   ln_bench   ! =T benchmark test with gyre: the gridsize is constant (no need to adjust timestep or viscosity)
   INTEGER, PUBLIC ::   nn_GYRE    ! 1/nn_GYRE = the resolution chosen in degrees and thus defining the horizontal domain size
   integer, PUBLIC :: AMR_ratio
   integer, PUBLIC :: nghost
   integer, PUBLIC :: JASMIN_restart_flag
   integer, PUBLIC :: AMR_time_ratio
   INTEGER, PUBLIC  ::   level_number
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_nam.F90 11536 2019-09-11 13:54:18Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE  usr_def_nam( cd_cfg, kk_cfg, kpi, kpj, kpk, kperio )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!                    
      !! ** Purpose :   read user defined namelist and define the domain size
      !!
      !! ** Method  :   read in namusr_def containing all the user specific namelist parameter
      !!
      !!                Here GYRE configuration
      !!
      !! ** input   : - namusr_def namelist found in namelist_cfg
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(out) ::   cd_cfg          ! configuration name
      INTEGER         , INTENT(out) ::   kk_cfg          ! configuration resolution
      INTEGER         , INTENT(out) ::   kpi, kpj, kpk   ! global domain sizes 
      INTEGER         , INTENT(out) ::   kperio          ! lateral global domain b.c. 
      !
      INTEGER ::   ios   ! Local integer
      !!
      NAMELIST/namusr_def/ nn_GYRE, ln_bench, jpkglo
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_cfg )          ! Namelist namusr_def (exist in namelist_cfg only)
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist' )
      !
      IF(lwm)   WRITE( numond, namusr_def )
      !
      cd_cfg = 'GYRE'               ! name & resolution (not used)
      kk_cfg = nn_GYRE
      !
write(*,*) "nghost=",nghost
      kpi = 30 * nn_GYRE  + 2        ! Global Domain size
      kpj = 20 * nn_GYRE + 2
      !kpi = 30 * nn_GYRE + 2        ! Global Domain size
      !kpj = 20 * nn_GYRE + 2
      kpk = jpkglo
      !                             ! Set the lateral boundary condition of the global domain
      kperio = 0                    ! GYRE configuration : closed domain
      !
      !                             ! control print
      IF(lwp) THEN
         WRITE(numout,*) '   '
         WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'
         WRITE(numout,*) '~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namusr_def : source case'
         WRITE(numout,*) '      GYRE used as Benchmark (=T)                      ln_bench  = ', ln_bench
         WRITE(numout,*) '      inverse resolution & implied domain size         nn_GYRE   = ', nn_GYRE
         WRITE(numout,*) '         jpiglo = 30*nn_GYRE+2                            jpiglo = ', kpi
         WRITE(numout,*) '         jpjglo = 20*nn_GYRE+2                            jpjglo = ', kpj
         WRITE(numout,*) '      number of model levels                              jpkglo = ', kpk
         WRITE(numout,*) '   '
         WRITE(numout,*) '   Lateral b.c. of the global domain set to closed        jperio = ', kperio
      ENDIF
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
