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

MODULE sms_c14
   !!======================================================================
   !!                      ***  MODULE trcsms_c14  ***
   !! TOP :  C14 main module
   !!======================================================================
   !! History     -   ! 1994-05 ( J. Orr ) original code
   !!            1.0  ! 2006-02 ( J.M. Molines )  Free form + modularity
   !!            2.0  ! 2008-12 ( C. Ethe ) reorganisation
   !!            4.0  ! 2011-02 ( A.R. Porter, STFC Daresbury ) Dynamic memory
   !!                 ! 2015    (A. Mouchet) general C14 + update formulas
   !!----------------------------------------------------------------------
   !!   sms_c14 :  compute and add C14 suface forcing to C14 trends
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc


   IMPLICIT NONE
   PUBLIC  


   LOGICAL  :: ln_chemh                             ! Chemical enhancement (yes/no)
   INTEGER  :: kc14typ                              ! C14 tracer type
   REAL(wp) :: tyrc14_beg                           ! year start atmospheric scenario !! See below
   REAL(wp) :: pco2at, rc14at                       ! atm co2, atm 14C ratio (global, reference)
   REAL(wp) :: rc14init                             ! ocean 14C ratio for initialization
   REAL(wp) :: xkwind, xdicsur                      ! wind coeff, ref DIC
   REAL(wp) :: rlam14                               ! C14 decay  rate

   !
   CHARACTER (len=20)                           :: cfileco2, cfilec14  ! Name of atmospheric forcing files
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   c14sbc   ! atmospheric c14 ratio
   REAL(wp)                                     ::   co2sbc   ! atmospheric co2 pressure
  
   REAL(wp),  ALLOCATABLE, SAVE, DIMENSION(:,:) ::   exch_c14   ! exch. vel. for C14/C
   REAL(wp),  ALLOCATABLE, SAVE, DIMENSION(:,:) ::   exch_co2   ! CO2 invasion rate
   REAL(wp),  ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qtr_c14    ! flux at surface
   REAL(wp),  ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qint_c14   ! cumulative flux

   INTEGER , PARAMETER                          ::   nc14zon     = 3  ! number of zones for bomb c14
   !
   INTEGER                                       ::   nrecco2, nrecc14  ! nb record atm co2 & cc14
   REAL(wp)                                      ::   tyrc14_now ! current yr for transient experiment relative to tyrc14_beg
   INTEGER                                       ::   m1_co2, m1_c14  ! index of first co2 and c14 records to consider
   INTEGER                                       ::   m2_co2, m2_c14  ! index of second co2 and c14 records to consider
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   bomb       ! C14 atm data (bomb - 3 zones)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)     ::   atmc14     ! C14 atm data (paleo - 1 zone)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)     ::   tyrc14     ! Time (yr) atmospheric C14 data
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   fareaz     ! Spatial Interpolation Factors
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)     ::   spco2      ! Atmospheric CO2
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)     ::   tyrco2     ! Time (yr) atmospheric CO2 data

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: sms_c14.F90 10071 2018-08-28 14:49:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   INTEGER FUNCTION sms_c14_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_sms_c14_alloc  ***
      !!----------------------------------------------------------------------
      sms_c14_alloc = 0
      ALLOCATE( exch_c14(jpi,jpj)        ,  exch_co2(jpi,jpj)        ,   &
         &      qtr_c14(jpi,jpj)         ,  qint_c14(jpi,jpj)        ,   &
         &      c14sbc(jpi,jpj)          ,  STAT = sms_c14_alloc )
         !
      !
   END FUNCTION sms_c14_alloc

  !!======================================================================
END MODULE sms_c14
