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

MODULE obs_profiles
   !!=====================================================================
   !!                       ***  MODULE  obs_profiles  ***
   !! Observation diagnostics: Storage space for profile observations
   !!                          arrays and additional flags etc.
   !!=====================================================================
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: obs_profiles.F90 10068 2018-08-28 14:09:04Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

   
   !! * Modules used 
   USE obs_profiles_def ! Definition of profile data types and tools

   IMPLICIT NONE
   
   SAVE

   !! * Routine accessibility
   PRIVATE

   PUBLIC nprofsets, nprofvars, nprofextr, profdata, prodatqc
   PUBLIC nvelosets, nvelovars, nveloextr, velodata, veldatqc
   
   !! * Shared Module variables
   INTEGER :: nprofsets                    ! Total number of profile data sets
   INTEGER :: nprofvars                    ! Total number of variables for profiles
   INTEGER :: nprofextr                    ! Extra fields for each variable
   TYPE(obs_prof), POINTER ::  profdata(:) ! Initial profile data
   TYPE(obs_prof), POINTER ::  prodatqc(:) ! Profile data after quality control

   INTEGER :: nvelosets                     ! Total number of velocity profile data sets
   INTEGER :: nvelovars                     ! Total number of variables for profiles
   INTEGER :: nveloextr                     ! Extra fields for each variable
   TYPE(obs_prof), POINTER ::  velodata(:)  ! Initial velocity profile data
   TYPE(obs_prof), POINTER ::  veldatqc(:)  ! Velocity profile data after quality control
END MODULE obs_profiles
