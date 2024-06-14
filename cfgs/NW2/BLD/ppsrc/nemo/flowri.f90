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

MODULE flowri
   !!======================================================================
   !!                       ***  MODULE  flowri  ***
   !!
   !! Ocean floats: write floats trajectory in ascii                    ln_flo_ascii = T
   !!                                    or in netcdf ( IOM or IOSPSL ) ln_flo_ascii = F           
   !!======================================================================
   !!  History :  OPA  !  1999-09  (Y. Drillet)    : Original code
   !!              -   !  2000-06  (J.-M. Molines) : Profiling floats for CLS 
   !!   NEMO      1.0  !  2002-10  (A. Bozec)  F90 : Free form and module
   !!             3.2  !  2010-08  (slaw, cbricaud): netcdf outputs and others 
   !!----------------------------------------------------------------------
   USE flo_oce         ! ocean drifting floats
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE lib_mpp         ! distribued memory computing library
   USE in_out_manager  ! I/O manager
   USE phycst          ! physic constants
   USE dianam          ! build name of file (routine)
   USE ioipsl
   USE iom             ! I/O library

   IMPLICIT NONE
   PRIVATE

   PUBLIC flo_wri         ! routine called by floats.F90
   PUBLIC flo_wri_alloc   ! routine called by floats.F90

   INTEGER :: jfl                            ! number of floats
   CHARACTER (len=80)  :: clname             ! netcdf output filename

   REAL(wp), ALLOCATABLE, DIMENSION(:) ::   zlon , zlat, zdep   ! 2D workspace
   REAL(wp), ALLOCATABLE, DIMENSION(:) ::   ztem , zsal, zrho   ! 2D workspace

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: flowri.F90 11819 2019-10-29 09:25:15Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION flo_wri_alloc()
      !!-------------------------------------------------------------------
      !!                ***  FUNCTION flo_wri_alloc  ***
      !!-------------------------------------------------------------------
      ALLOCATE( ztem(jpnfl) , zsal(jpnfl) , zrho(jpnfl) , &
                zlon(jpnfl) , zlat(jpnfl) , zdep(jpnfl) , STAT=flo_wri_alloc)
      !  
      CALL mpp_sum ( 'flowri', flo_wri_alloc )
      IF( flo_wri_alloc /= 0 )   CALL ctl_stop( 'STOP', 'flo_wri_alloc: failed to allocate arrays.' )
   END FUNCTION flo_wri_alloc

   SUBROUTINE flo_wri( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE flo_wri ***
      !!             
      !! ** Purpose :   Write position of floats in "trajec_float.nc",according
      !!                to ARIANE TOOLS (http://stockage.univ-brest.fr/~grima/Ariane/ )  n
      !!                nomenclature
      !!    
      !!      
      !! ** Method  :   The frequency of  ??? is nwritefl
      !!      
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER  :: kt                               ! time step

      !! * Local declarations
      INTEGER  :: iafl , ibfl , icfl             ! temporary integer
      INTEGER  :: ia1fl, ib1fl, ic1fl            !   "
      INTEGER  :: iafloc,ibfloc,ia1floc,ib1floc  !   "
      INTEGER  :: irec, irecflo

      REAL(wp) :: zafl,zbfl,zcfl                 ! temporary real
      REAL(wp) :: ztime                          !   "

      INTEGER, DIMENSION(2)          :: icount
      INTEGER, DIMENSION(2)          :: istart
      INTEGER, DIMENSION(1)          :: ish
      INTEGER, DIMENSION(2)          :: ish2
      !!----------------------------------------------------------------------
      
      !-----------------------------------------------------
      ! I- Save positions, temperature, salinty and density 
      !-----------------------------------------------------
      zlon(:)=0.0 ; zlat(:)=0.0 ; zdep(:)=0.0 
      ztem(:)=0.0 ; zsal(:)=0.0 ; zrho(:)=0.0 

      DO jfl = 1, jpnfl

         iafl  = INT (tpifl(jfl))            ! I-index of the nearest point before
         ibfl  = INT (tpjfl(jfl))            ! J-index of the nearest point before
         icfl  = INT (tpkfl(jfl))            ! K-index of the nearest point before
         ia1fl = iafl + 1                    ! I-index of the nearest point after
         ib1fl = ibfl + 1                    ! J-index of the nearest point after
         ic1fl = icfl + 1                    ! K-index of the nearest point after
         zafl  = tpifl(jfl) - REAL(iafl,wp)  ! distance  ?????
         zbfl  = tpjfl(jfl) - REAL(ibfl,wp)  ! distance  ?????
         zcfl  = tpkfl(jfl) - REAL(icfl,wp)  ! distance  ?????

         IF( lk_mpp ) THEN
               
            iafloc = mi1( iafl )
            ibfloc = mj1( ibfl )
 
            IF( nldi <= iafloc .AND. iafloc <= nlei .AND. &
              & nldj <= ibfloc .AND. ibfloc <= nlej       ) THEN 

               !the float is inside of current proc's area
               ia1floc = iafloc + 1
               ib1floc = ibfloc + 1
     
               !save position of the float
               zlat(jfl) = (1.-zafl)*(1.-zbfl)*gphit(iafloc ,ibfloc ) + (1.-zafl) * zbfl * gphit(iafloc ,ib1floc)   &
                     +     zafl *(1.-zbfl)*gphit(ia1floc,ibfloc ) +     zafl  * zbfl * gphit(ia1floc,ib1floc)   
               zlon(jfl) = (1.-zafl)*(1.-zbfl)*glamt(iafloc ,ibfloc ) + (1.-zafl) * zbfl * glamt(iafloc ,ib1floc)   &
                     +     zafl *(1.-zbfl)*glamt(ia1floc,ibfloc ) +     zafl  * zbfl * glamt(ia1floc,ib1floc)
               zdep(jfl) = (1.-zcfl)*gdepw_n(iafloc,ibfloc,icfl ) + zcfl * gdepw_n(iafloc,ibfloc,ic1fl)     

               !save temperature, salinity and density at this position
               ztem(jfl) = tsn(iafloc,ibfloc,icfl,jp_tem)
               zsal (jfl) = tsn(iafloc,ibfloc,icfl,jp_sal)
               zrho (jfl) = (rhd(iafloc,ibfloc,icfl)+1)*rau0

            ENDIF

         ELSE  ! mono proc case  

            iafloc  = iafl
            ibfloc  = ibfl
            ia1floc = iafloc + 1
            ib1floc = ibfloc + 1

            !save position of the float               
            zlat(jfl) = (1.-zafl)*(1.-zbfl)*gphit(iafloc ,ibfloc ) + (1.-zafl) * zbfl * gphit(iafloc ,ib1floc)   &
                      +     zafl *(1.-zbfl)*gphit(ia1floc,ibfloc ) +     zafl  * zbfl * gphit(ia1floc,ib1floc)
            zlon(jfl) = (1.-zafl)*(1.-zbfl)*glamt(iafloc ,ibfloc ) + (1.-zafl) * zbfl * glamt(iafloc ,ib1floc)   &
                      +     zafl *(1.-zbfl)*glamt(ia1floc,ibfloc ) +     zafl  * zbfl * glamt(ia1floc,ib1floc)
            zdep(jfl) = (1.-zcfl)*gdepw_n(iafloc,ibfloc,icfl ) + zcfl * gdepw_n(iafloc,ibfloc,ic1fl)

            ztem(jfl) = tsn(iafloc,ibfloc,icfl,jp_tem)
            zsal(jfl) = tsn(iafloc,ibfloc,icfl,jp_sal)
            zrho(jfl) = (rhd(iafloc,ibfloc,icfl)+1)*rau0
          
         ENDIF

      END DO ! loop on float

      !Only proc 0 writes all positions : SUM of positions on all procs
      IF( lk_mpp ) THEN
         CALL mpp_sum( 'flowri', zlon, jpnfl )   ! sums over the global domain
         CALL mpp_sum( 'flowri', zlat, jpnfl )   ! sums over the global domain
         CALL mpp_sum( 'flowri', zdep, jpnfl )   ! sums over the global domain
         CALL mpp_sum( 'flowri', ztem, jpnfl )   ! sums over the global domain
         CALL mpp_sum( 'flowri', zsal, jpnfl )   ! sums over the global domain
         CALL mpp_sum( 'flowri', zrho, jpnfl )   ! sums over the global domain
      ENDIF


      !-------------------------------------!
      ! II- WRITE WRITE WRITE WRITE WRITE   !
      !-------------------------------------!

      !--------------------------!
      ! II-1 Write in ascii file !
      !--------------------------!

      IF( ln_flo_ascii )THEN

         IF( ( kt == nn_it000 .OR. MOD( kt,nn_writefl)== 0 ) .AND. lwp )THEN

            !II-1-a Open ascii file
            !----------------------
            IF( kt == nn_it000 ) THEN
               CALL ctl_opn( numflo, 'trajec_float', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
               irecflo = NINT( (nitend-nn_it000) / FLOAT(nn_writefl) )
               WRITE(numflo,*) cexper, irecflo, jpnfl, nn_writefl
            ENDIF

            !II-1-b Write in ascii file
            !-----------------------------
            WRITE(numflo,*) zlon,zlat,zdep,nisobfl,ngrpfl,ztem,zsal, FLOAT(ndastp)


            !II-1-c Close netcdf file
            !-------------------------
            IF( kt == nitend )   CLOSE( numflo )

         ENDIF

      !-----------------------------------------------------
      ! II-2 Write in netcdf file
      !-----------------------------------------------------

      ELSE

      !II-2-a Write with IOM
      !----------------------

         CALL iom_put( "traj_lon"     , zlon )
         CALL iom_put( "traj_lat"     , zlat )
         CALL iom_put( "traj_dep"     , zdep )
         CALL iom_put( "traj_temp"    , ztem )
         CALL iom_put( "traj_salt"    , zsal  )
         CALL iom_put( "traj_dens"    , zrho )
         CALL iom_put( "traj_group"   , REAL(ngrpfl,wp) )
      ENDIF ! netcdf writing
   
   END SUBROUTINE flo_wri

   !!=======================================================================
END MODULE flowri
