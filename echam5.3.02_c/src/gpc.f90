SUBROUTINE gpc(krow, kglat)

  ! Description:
  !
  ! Grid point computations.
  !
  ! Method:
  !
  ! This subroutine controls parts the computations in grid points, 
  ! that is the physical computations(*phys*) 
  ! and grid point contributions to the semi implicit scheme (*si1*).
  ! 
  !
  ! *gpc* is called from *scan1*.
  ! 
  ! Externals:
  !   *si1*       grid point contributions to the semi implicit scheme.
  !   *physc*     physical computations.
  !
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, January 1982, original source
  ! F. Lunkeit, MI, June 1989, CLSST added
  ! U. Schlese, MPI, July 1989, add seaice computations
  ! U. Schlese, DKRZ, January 1995, initialization of soil temperatures
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schlese, DKRZ, and M. Esch, MPI July 1999, modifications for ECHAM5
  ! I. Kirchner, MPI Sepember 2000, nudging sst update
  ! I. Kirchner, MPI December 2000, time control
  ! S. Legutke, MPI M&D , Jan 2002, coupled case where ocean has no sea ice (FLUXES4)
  ! U. Schlese, M. Esch, MPI, September 2002, mixed layer ocean
  ! 
  ! for more details see file AUTHORS
  !

#ifdef OBSOLETE
  USE mo_control,       ONLY: lnudge, lcouple, lmlo, lso4
#else
  USE mo_control,       ONLY: lnudge, lcouple
#endif
  USE mo_nudging_sst,   ONLY: NudgingSSTnew
!!#D mlocean +
#ifdef MESSY
  USE messy_main_switch, ONLY: USE_MLOCEAN  ! fb_mk_20110209
#endif
!!#D mlocean -

  IMPLICIT NONE

  INTEGER  :: krow, kglat

  !  Executable statements 

!-- 1. Distribute climate values

  IF (.NOT. lcouple) THEN
     IF (lnudge) THEN
        CALL NudgingSSTnew(krow)
     ELSE
#ifdef OBSOLETE
        IF (lmlo) THEN
          CALL ml_flux(krow)
        ELSE
#endif
!!#D mlocean +
#ifndef MESSY
!!#D mlocean -
          CALL clsst(krow)
!!#D mlocean +
#else
          IF (.NOT.USE_MLOCEAN) CALL clsst(krow) ! fb_mk_20110209
#endif
!!#D mlocean -
#ifdef OBSOLETE
        END IF
#endif
     END IF
#ifdef PFLUXES4
  ELSE
     ! Coupled mode when only sst is passed from the ocean (e.g. regional ocean)
     !         seaice is interpolated from climatology
     !         sst is set to ctfreez if sea ice exists in climatology
     CALL clsst(krow)
#endif
  END IF

  CALL clveg(krow)
#ifdef OBSOLETE
  IF(lso4) CALL intaero(krow)
#endif
!-- 2. Parametrisation of diabatic processes

  CALL physc(krow, kglat)

  RETURN
END SUBROUTINE gpc
