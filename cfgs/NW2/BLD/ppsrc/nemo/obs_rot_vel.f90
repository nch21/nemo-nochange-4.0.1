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

MODULE obs_rot_vel
   !!======================================================================
   !!                       ***  MODULE obs_rot_vel  ***
   !! Observation diagnostics: Read the velocity profile observations
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   obs_rotvel : Rotate velocity data into N-S,E-W directorions
   !!----------------------------------------------------------------------
   !! * Modules used   
   USE par_kind                 ! Precision variables
   USE par_oce                  ! Ocean parameters
   USE in_out_manager           ! I/O manager
   USE dom_oce                  ! Ocean space and time domain variables
   USE obs_grid                 ! Grid search
   USE obs_utils                ! For error handling
   USE obs_profiles_def         ! Profile definitions
   USE obs_inter_h2d            ! Horizontal interpolation
   USE obs_inter_sup            ! MPP support routines for interpolation
   USE geo2ocean                ! Rotation of vectors
   USE obs_fbm                  ! Feedback definitions

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE

   PUBLIC obs_rotvel            ! Rotate the observations

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: obs_rot_vel.F90 10068 2018-08-28 14:09:04Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE obs_rotvel( profdata, k2dint, pu, pv )
      !!---------------------------------------------------------------------
      !!
      !!                   *** ROUTINE obs_rea_pro_dri ***
      !!
      !! ** Purpose : Rotate velocity data into N-S,E-W directorions
      !!
      !! ** Method  : Interpolation of geo2ocean coefficients on U,V grid
      !!              to observation point followed by a similar computations
      !!              as in geo2ocean.
      !!
      !! ** Action  : Review if there is a better way to do this.
      !!
      !! References : 
      !!
      !! History :  
      !!      ! :  2009-02 (K. Mogensen) : New routine
      !!----------------------------------------------------------------------
      !! * Modules used
      !! * Arguments
      TYPE(obs_prof), INTENT(INOUT) :: profdata    ! Profile data to be read
      INTEGER, INTENT(IN) :: k2dint     ! Horizontal interpolation methed
      REAL(wp), DIMENSION(*) :: &
         & pu, &
         & pv
      !! * Local declarations
      REAL(wp), DIMENSION(2,2,1) :: zweig
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zmasku, &
         & zmaskv, &
         & zcoslu, &
         & zsinlu, &
         & zcoslv, &
         & zsinlv, &
         & zglamu, &
         & zgphiu, &
         & zglamv, &
         & zgphiv
      REAL(wp), DIMENSION(1) :: &
         & zsinu, &
         & zcosu, &
         & zsinv, &
         & zcosv
      REAL(wp) :: zsin
      REAL(wp) :: zcos
      REAL(wp), DIMENSION(1) :: zobsmask
      REAL(wp), DIMENSION(jpi,jpj) :: zsingu,zcosgu,zsingv,zcosgv
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdiu, &
         & igrdju, &
         & igrdiv, &
         & igrdjv
      INTEGER :: ji
      INTEGER :: jk


      !-----------------------------------------------------------------------
      ! Allocate data for message parsing and interpolation
      !-----------------------------------------------------------------------

      ALLOCATE( &
         & igrdiu(2,2,profdata%nprof), &
         & igrdju(2,2,profdata%nprof), &
         & zglamu(2,2,profdata%nprof), &
         & zgphiu(2,2,profdata%nprof), &
         & zmasku(2,2,profdata%nprof), &
         & zcoslu(2,2,profdata%nprof), &
         & zsinlu(2,2,profdata%nprof), &
         & igrdiv(2,2,profdata%nprof), &
         & igrdjv(2,2,profdata%nprof), &
         & zglamv(2,2,profdata%nprof), &
         & zgphiv(2,2,profdata%nprof), &
         & zmaskv(2,2,profdata%nprof), &
         & zcoslv(2,2,profdata%nprof), &
         & zsinlv(2,2,profdata%nprof)  &
         & )

      !-----------------------------------------------------------------------
      ! Receive the angles on the U and V grids.
      !-----------------------------------------------------------------------

      CALL obs_rot( zsingu, zcosgu, zsingv, zcosgv )

      DO ji = 1, profdata%nprof
         igrdiu(1,1,ji) = profdata%mi(ji,1)-1
         igrdju(1,1,ji) = profdata%mj(ji,1)-1
         igrdiu(1,2,ji) = profdata%mi(ji,1)-1
         igrdju(1,2,ji) = profdata%mj(ji,1)
         igrdiu(2,1,ji) = profdata%mi(ji,1)
         igrdju(2,1,ji) = profdata%mj(ji,1)-1
         igrdiu(2,2,ji) = profdata%mi(ji,1)
         igrdju(2,2,ji) = profdata%mj(ji,1)
         igrdiv(1,1,ji) = profdata%mi(ji,2)-1
         igrdjv(1,1,ji) = profdata%mj(ji,2)-1
         igrdiv(1,2,ji) = profdata%mi(ji,2)-1
         igrdjv(1,2,ji) = profdata%mj(ji,2)
         igrdiv(2,1,ji) = profdata%mi(ji,2)
         igrdjv(2,1,ji) = profdata%mj(ji,2)-1
         igrdiv(2,2,ji) = profdata%mi(ji,2)
         igrdjv(2,2,ji) = profdata%mj(ji,2)
      END DO

      CALL obs_int_comm_2d( 2, 2, profdata%nprof, jpi, jpj, igrdiu, igrdju, &
         &                  glamu, zglamu )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, jpi, jpj, igrdiu, igrdju, &
         &                  gphiu, zgphiu )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, jpi, jpj, igrdiu, igrdju, &
         &                  umask(:,:,1), zmasku )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, jpi, jpj, igrdiu, igrdju, &
         &                  zsingu, zsinlu )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, jpi, jpj, igrdiu, igrdju, &
         &                  zcosgu, zcoslu )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, jpi, jpj, igrdiv, igrdjv, &
         &                  glamv, zglamv )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, jpi, jpj, igrdiv, igrdjv, &
         &                  gphiv, zgphiv )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, jpi, jpj, igrdiv, igrdjv, &
         &                  vmask(:,:,1), zmaskv )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, jpi, jpj, igrdiv, igrdjv, &
         &                  zsingv, zsinlv )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, jpi, jpj, igrdiv, igrdjv, &
         &                  zcosgv, zcoslv )

      DO ji = 1, profdata%nprof
            
         CALL obs_int_h2d_init( 1, 1, k2dint, &
            &                   profdata%rlam(ji), profdata%rphi(ji), &
            &                   zglamu(:,:,ji), zgphiu(:,:,ji), &
            &                   zmasku(:,:,ji), zweig, zobsmask )
         
         CALL obs_int_h2d( 1, 1, zweig, zsinlu(:,:,ji),  zsinu )

         CALL obs_int_h2d( 1, 1, zweig, zcoslu(:,:,ji),  zcosu )

         CALL obs_int_h2d_init( 1, 1, k2dint, &
            &                   profdata%rlam(ji), profdata%rphi(ji), &
            &                   zglamv(:,:,ji), zgphiv(:,:,ji), &
            &                   zmaskv(:,:,ji), zweig, zobsmask )
         
         CALL obs_int_h2d( 1, 1, zweig, zsinlv(:,:,ji),  zsinv )

         CALL obs_int_h2d( 1, 1, zweig, zcoslv(:,:,ji),  zcosv )

         ! Assume that the angle at observation point is the 
         ! mean of u and v cosines/sines

         zcos = 0.5_wp * ( zcosu(1) + zcosv(1) )
         zsin = 0.5_wp * ( zsinu(1) + zsinv(1) )
         
         IF ( ( profdata%npvsta(ji,1) /= profdata%npvsta(ji,2) ) .OR. &
            & ( profdata%npvend(ji,1) /= profdata%npvend(ji,2) ) ) THEN
            CALL fatal_error( 'Different number of U and V observations '// &
               'in a profile in obs_rotvel', 190 )
         ENDIF

         DO jk = profdata%npvsta(ji,1), profdata%npvend(ji,1)
            IF ( ( profdata%var(1)%vmod(jk) /= fbrmdi ) .AND. &
               & ( profdata%var(2)%vmod(jk) /= fbrmdi ) ) THEN
               pu(jk) = profdata%var(1)%vmod(jk) * zcos - &
                  &     profdata%var(2)%vmod(jk) * zsin
               pv(jk) = profdata%var(2)%vmod(jk) * zcos + &
                  &     profdata%var(1)%vmod(jk) * zsin
            ELSE
               pu(jk) = fbrmdi
               pv(jk) = fbrmdi
            ENDIF

         END DO

      END DO
      
      DEALLOCATE( &
         & igrdiu, &
         & igrdju, &
         & zglamu, &
         & zgphiu, &
         & zmasku, &
         & zcoslu, &
         & zsinlu, &
         & igrdiv, &
         & igrdjv, &
         & zglamv, &
         & zgphiv, &
         & zmaskv, &
         & zcoslv, &
         & zsinlv  &
         & )

   END SUBROUTINE obs_rotvel

END MODULE obs_rot_vel
