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

MODULE sedarr
   !!======================================================================
   !!                       ***  MODULE sedarr   ***
   !!              transform 1D (2D) array to a 2D (1D) table
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   arr_2d_1d  : 2-D to 1-D
   !!   arr_1d_2d  : 1-D to 2-D
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_sed
   USE dom_oce
   USE sed

   IMPLICIT NONE
   PRIVATE

   INTERFACE pack_arr
      MODULE PROCEDURE pack_arr_2d_1d , pack_arr_3d_2d 
   END INTERFACE

   INTERFACE unpack_arr
      MODULE PROCEDURE unpack_arr_1d_2d , unpack_arr_2d_3d 
   END INTERFACE

   !! * Routine accessibility
   PUBLIC pack_arr
   PUBLIC unpack_arr

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: sedarr.F90 10222 2018-10-25 09:42:23Z aumont $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE pack_arr_2d_1d ( ndim1d, tab1d, tab2d, tab_ind )

      INTEGER, INTENT(in) ::  ndim1d
      REAL(wp), DIMENSION (jpi, jpj), INTENT(in) :: tab2d
      INTEGER, DIMENSION (ndim1d), INTENT (in) :: tab_ind
      REAL(wp), DIMENSION(ndim1d), INTENT (out) ::  tab1d

      INTEGER ::  jn, jid, jjd

      IF( ln_timing )   CALL timing_start('pack_arr_2d_1d')
        
      DO jn = 1, ndim1d
         jid        = MOD( tab_ind(jn) - 1, jpi ) + 1
         jjd        = ( tab_ind(jn) - 1 ) / jpi + 1
         tab1d(jn)  = tab2d(jid, jjd)
      END DO 

      IF( ln_timing )   CALL timing_stop('pack_arr_2d_1d')

   END SUBROUTINE pack_arr_2d_1d

   SUBROUTINE unpack_arr_1d_2d ( ndim1d, tab2d, tab_ind, tab1d )

      INTEGER, INTENT ( in) ::  ndim1d
      INTEGER, DIMENSION (ndim1d) , INTENT (in) ::   tab_ind
      REAL(wp), DIMENSION(ndim1d), INTENT (in) ::   tab1d  
      REAL(wp), DIMENSION (jpi, jpj), INTENT ( out) ::  tab2d
      INTEGER ::  jn, jid, jjd

      IF( ln_timing )   CALL timing_start('unpack_arr_1d_2d')

      DO jn = 1, ndim1d
         jid             = MOD( tab_ind(jn) - 1, jpi) + 1
         jjd             =    ( tab_ind(jn) - 1 ) / jpi  + 1
         tab2d(jid, jjd) = tab1d(jn)
      END DO

      IF( ln_timing )   CALL timing_stop('unpack_arr_1d_2d')

   END SUBROUTINE unpack_arr_1d_2d

   SUBROUTINE pack_arr_3d_2d ( ndim1d, tab2d, tab3d, tab_ind )

      INTEGER, INTENT(in) :: ndim1d      
      REAL(wp), DIMENSION(jpi,jpj,jpksed), INTENT(in) ::   tab3d      
      INTEGER, DIMENSION(ndim1d), INTENT (in) ::   tab_ind      
      REAL(wp), DIMENSION(ndim1d,jpksed), INTENT (out) ::    tab2d 
      INTEGER, DIMENSION(ndim1d) ::  jid, jjd        
      INTEGER ::    jk, jn , ji, jj

      IF( ln_timing )   CALL timing_start('pack_arr_2d_3d')
      
      DO jn = 1, ndim1d
         jid(jn) = MOD( tab_ind(jn) - 1, jpi ) + 1
         jjd(jn) = ( tab_ind(jn) - 1 ) / jpi + 1
      END DO
 
      DO jk = 1, jpksed
         DO jn = 1, ndim1d
            ji = jid(jn)
            jj = jjd(jn)
            tab2d(jn,jk)  = tab3d(ji,jj,jk) 
         ENDDO
      ENDDO

      IF( ln_timing )   CALL timing_stop('pack_arr_2d_3d')
      
   END SUBROUTINE pack_arr_3d_2d


   SUBROUTINE unpack_arr_2d_3d ( ndim1d, tab3d, tab_ind, tab2d )

      INTEGER, INTENT(in) :: ndim1d      
      REAL(wp), DIMENSION(ndim1d,jpksed), INTENT(in) ::   tab2d      
      INTEGER, DIMENSION(ndim1d), INTENT (in) ::   tab_ind      
      REAL(wp), DIMENSION(jpi,jpj,jpksed), INTENT (out) ::    tab3d 
      INTEGER, DIMENSION(ndim1d) ::  jid, jjd        
      INTEGER ::   jk, jn , ji, jj
      !
      IF( ln_timing )   CALL timing_start('unpack_arr_2d_3d')
      !
      DO jn = 1, ndim1d
         jid(jn) = MOD( tab_ind(jn) - 1, jpi ) + 1
         jjd(jn) = ( tab_ind(jn) - 1 ) / jpi + 1
      END DO
 
      DO jk = 1, jpksed
         DO jn = 1, ndim1d
            ji = jid(jn)
            jj = jjd(jn)
            tab3d(ji, jj, jk) = tab2d(jn,jk)
         ENDDO
      ENDDO

      IF( ln_timing )   CALL timing_stop('unpack_arr_2d_3d')

   END SUBROUTINE unpack_arr_2d_3d

END MODULE sedarr
