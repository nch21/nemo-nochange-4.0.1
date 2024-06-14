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

MODULE obs_sort
   !!=====================================================================
   !!                       ***  MODULE obs_sort  ***
   !! Observation diagnostics: Various tools for sorting etc.
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   sort_dp_indx : Get indicies for ascending order for a double prec. array
   !!   index_sort   : Get indicies for ascending order for a double prec. array
   !!---------------------------------------------------------------------
   !! * Modules used
   USE par_kind, ONLY : & ! Precision variables
      & dp
  
   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE index_sort    ! Get indicies for ascending order for a double prec. array
   
   PUBLIC sort_dp_indx   ! Get indicies for ascending order for a double prec. array
  
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: obs_sort.F90 10068 2018-08-28 14:09:04Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE sort_dp_indx( kvals, pvals, kindx )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE sort_dp_indx  ***
      !!          
      !! ** Purpose : Get indicies for ascending order for a double precision array
      !!
      !! ** Method  : Call index_sort routine
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  06-05  (K. Mogensen)  Original code
      !!        !  06-10  (A. Weaver) Cleaning
      !!----------------------------------------------------------------------

      !! * Arguments
      INTEGER, INTENT(IN) :: kvals     ! Number of elements to be sorted
      REAL(KIND=dp), DIMENSION(kvals), INTENT(IN) :: &
         & pvals            ! Array to be sorted
      INTEGER, DIMENSION(kvals), INTENT(OUT) ::  &
         & kindx            ! Indices for ordering of array

      !! * Local declarations

      !-----------------------------------------------------------------------
      ! Call qsort routine
      !-----------------------------------------------------------------------
      IF (kvals>=1) THEN

         CALL index_sort( pvals, kindx, kvals )

      ENDIF

   END SUBROUTINE sort_dp_indx

   SUBROUTINE index_sort( pval, kindx, kvals )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE index_sort  ***
      !!          
      !! ** Purpose : Get indicies for ascending order for a double precision array
      !!
      !! ** Method  : Heapsort
      !!
      !! ** Action  : 
      !!
      !! References : http://en.wikipedia.org/wiki/Heapsort
      !!
      !! History :
      !!        !  06-05  (K. Mogensen)  Original code
      !!        !  06-10  (A. Weaver) Cleaning
      !!----------------------------------------------------------------------

      !! * Arguments
      INTEGER, INTENT(IN) :: kvals         ! Number of values
      REAL(KIND=dp), DIMENSION(kvals), INTENT(IN) :: &
         & pval                            ! Array to be sorted
      INTEGER, DIMENSION(kvals), INTENT(INOUT) :: &
         & kindx                           ! Indicies for ordering

      !! * Local declarations
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jt
      INTEGER :: jn
      INTEGER :: jparent
      INTEGER :: jchild

      DO ji = 1, kvals
         kindx(ji) = ji
      END DO
      
      ji = kvals/2 + 1
      jn = kvals

      main_loop : DO

         IF ( ji > 1 ) THEN
            ji = ji-1
            jt = kindx(ji)
         ELSE
            jt = kindx(jn)
            kindx(jn) = kindx(1)
            jn = jn-1
            IF ( jn <= 1 ) THEN
               kindx(1) = jt
               EXIT main_loop
            ENDIF
         ENDIF

         jparent = ji
         jchild =  2 * ji

         inner_loop : DO

            IF ( jchild > jn ) EXIT inner_loop
            IF ( jchild < jn ) THEN
               IF ( pval(kindx(jchild)) < pval(kindx(jchild+1)) ) THEN
                 jchild = jchild+1
               ENDIF
            ENDIF
            IF  ( pval(jt) < pval(kindx(jchild))) THEN
               kindx(jparent) = kindx(jchild)
               jparent = jchild
               jchild  = jchild*2
            ELSE 
               jchild = jn + 1 
            ENDIF

         END DO inner_loop

         kindx(jparent) = jt

      END DO main_loop
      
   END SUBROUTINE index_sort

END MODULE obs_sort
 
