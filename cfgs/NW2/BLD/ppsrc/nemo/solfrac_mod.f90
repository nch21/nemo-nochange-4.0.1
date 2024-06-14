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

MODULE solfrac_mod
   !!======================================================================
   !!                    ***  MODULE  solfrac  ***
   !!     POSH representation of solar absorption (Gntermann, 2009)
   !!=====================================================================
   !! History :        !  11-10  (J. While)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   solfrac  : function to calculate the solar fraction
   !!----------------------------------------------------------------------
   
   USE par_kind
   IMPLICIT NONE
     
   ! Parameters 
   REAL(wp), PRIVATE, PARAMETER, DIMENSION(9) :: &
   &                                     pp_wgt = (/0.2370, 0.36,  0.1790, &
   &                                                0.087,  0.08,  0.025,  &
   &                                                0.025,  0.007, 0.0004/)
   REAL(wp), PRIVATE, PARAMETER, DIMENSION(9) :: &
   &                                    pp_len = (/34.84,   2.266,   0.0315,  &
   &                                               0.0055,  8.32e-4, 1.26e-4, &
   &                                                3.13e-4, 7.82e-4, 1.44e-5/)
   
   PUBLIC solfrac
                                                   
CONTAINS

   REAL(dp) FUNCTION solfrac(ptop,pbottom)
       !!----------------------------------------------------------------------
      !! *** ROUTINE solfrac ***
      !!
      !! ** Purpose :   Calculate the solar fraction absorbed between two 
      !!                layers
      !!
      !! ** Reference : POSH a model of diurnal warming, Gentemann et al, 
      !!                 JGR, 2009 
      !!----------------------------------------------------------------------
      
      ! Dummy variabes
      REAL(wp), INTENT(IN) :: ptop, pbottom   ! Top and bottom of layer
      
      ! local variables
      INTEGER :: jt
      
      ! Calculate the solar fraction absorbed between the two layers
      solfrac = 0._wp
      DO jt = 1, 9
           solfrac = solfrac + pp_wgt(jt) * ( exp ( -ptop / pp_len(jt) ) &
            &                                 - exp ( -pbottom / pp_len(jt) ) )
      END DO
      
   END FUNCTION
   
END MODULE solfrac_mod
