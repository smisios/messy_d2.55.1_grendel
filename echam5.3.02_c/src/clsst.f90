SUBROUTINE clsst(krow)

  ! Description:
  !
  ! Passes climate sea-surface-temperatures and sea ice to atmosphere
  !
  ! Method:
  !
  ! This subroutine interpolates the sea-surface temperatures and
  ! sea-ice concentration at each time step and updates *tsw* and *siced*.
  !
  ! *clsst* is called from *gpc*.
  !
  ! Authors: 
  !
  ! U. Schlese, DKRZ, January 1993, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schlese, DKRZ and M. Esch, MPI, Oct. 1999, modifications for ECHAM5
  ! R. Voss, U. Schlese, MPI, Jan 2000, mods for fract. surface coverage
  ! A. Rhodin, MPI, prescribe sst in SCM run
  ! S. Legutke, MPI M&D , Jan 2002, modify coupling when ocean has no ice (FLUXES4)
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! for more details see file AUTHORS
  ! 
  !

  USE mo_kind,          ONLY: dp
  USE mo_memory_g3b,    ONLY: slf, tsw, tsi, seaice, siced, alake
  USE mo_sst,           ONLY: sst, aice
  USE mo_control,       ONLY: lcouple
#ifndef MESSY
  USE mo_physc2,        ONLY: ctfreez
#else
  USE messy_main_constants_mem, ONLY: ctfreez
#endif
  USE mo_constants,     ONLY: tmelt
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_interpo,       ONLY: nmw1, nmw2, wgt1, wgt2
#ifdef OBSOLETE
  USE mo_column,        ONLY: sst_1d ! if > 0 : prescribed SST
#endif
  USE mo_geoloc,        ONLY: philat_2d
#ifndef MESSY
  USE mo_radiation,     ONLY: nmonth
#else
  USE messy_main_timer_bi, ONLY: nmonth
#endif

  IMPLICIT NONE

  INTEGER :: krow

  ! Local variables

  INTEGER :: jl      ! longitude loop index
  INTEGER :: jrow    ! local latitude index
  INTEGER :: nproma  ! number of longitudes on PE
  REAL(dp):: zic
  REAL(dp):: zts
  LOGICAL :: lonorth ! .true. for northern latitude
#ifdef OBSOLETE
  LOGICAL :: lsst_1d
#endif

  !  Executable statements

  jrow    = krow        ! local latitude index

  IF ( jrow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  END IF

!-- 1. Update temperatures and sea ice

!-- 1.1.1 Uncoupled mode

  IF(.NOT.lcouple) THEN

#ifdef OBSOLETE
  lsst_1d = ALLOCATED(sst_1d)
#endif

!DIR$ CONCURRENT
  DO jl=1,nproma
!
!-- 1.1 Annual cycle
!
   IF(nmonth == 0) THEN
     zts=wgt1*sst(jl,jrow,nmw1)+wgt2*sst(jl,jrow,nmw2)
     zic=wgt1*aice(jl,jrow,nmw1)+wgt2*aice(jl,jrow,nmw2)
!
!-- 1.2 Perpetual month
!
   ELSE
     zts=sst(jl,jrow,nmonth)
     zic=aice(jl,jrow,nmonth)
   END IF
!


    IF(alake(jl,jrow).EQ.0._dp) THEN
      IF(slf(jl,jrow).LT.1._dp) THEN

        zic=zic*0.01_dp  !!!  assuming input data is in percent
        seaice(jl,jrow)=MAX(0._dp,MIN(0.99_dp,zic))
        IF (seaice(jl,jrow).LE.0.01_dp) seaice(jl,jrow)=0.0_dp

        IF (seaice(jl,jrow).GT.0._dp) THEN

          tsw(jl,jrow)=ctfreez

          lonorth = philat_2d(jl,jrow).GT.0 ! true in northern hemisphere
          IF (lonorth) THEN
            siced(jl,jrow)=2._dp
          ELSE
            siced(jl,jrow)=1._dp
          END IF

        ELSE                               ! water

          tsw(jl,jrow)=MAX(zts,ctfreez)
          siced(jl,jrow)=0._dp

        END IF
        !
#ifdef OBSOLETE
        ! prescribe sst in SCM
        !
        IF (lsst_1d) THEN
          IF (sst_1d(jl,jrow) > 0._dp) tsw(jl,jrow) = sst_1d(jl,jrow)
        ENDIF
#endif

      ELSE                                 ! land
        seaice(jl,jrow)=0._dp
        siced(jl,jrow)=0._dp
        tsw(jl,jrow)=tmelt  !! dummy setting to some reasonable value
        tsi(jl,jrow)=tmelt
      END IF
    END IF
  END DO

  END IF

#ifdef PFLUXES4
!-- 1.1.2 Coupled mode: with cppoption=PFLUXES4
!         only sst is passed from the ocean.
!         sst is set to ctfreez if sea ice exists in climatology
  DO jl=1,nproma

    zic=wgt1*aice(jl,jrow,nmw1)+wgt2*aice(jl,jrow,nmw2)

    IF(alake(jl,jrow).EQ.0._dp) THEN
      IF(slf(jl,jrow).LT.1._dp) THEN        ! ocean cell

        zic=zic*0.01_dp                     ! assuming input data is in percent
        seaice(jl,jrow)=MAX(0._dp,MIN(1._dp,zic))
        IF (seaice(jl,jrow).LE.0.01_dp) seaice(jl,jrow)=0.0_dp

        IF (seaice(jl,jrow).GT.0._dp) THEN ! sea ice exists

          tsw(jl,jrow)=ctfreez

          lonorth = philat_2d(jl,jrow).GT.0 ! true in northern hemisphere
          IF (lonorth) THEN
            siced(jl,jrow)=2._dp
          ELSE
            siced(jl,jrow)=1._dp
          END IF

        ELSE                             ! no sea ice

          siced(jl,jrow)=0._dp
          tsw(jl,jrow)=MAX(tsw(jl,jrow),ctfreez)

        END IF

      ELSE                               ! no ocean in cell
        seaice(jl,jrow)=0._dp
        siced(jl,jrow)=0._dp
        tsw(jl,jrow)=tmelt  !! dummy setting to some reasonable value
        tsi(jl,jrow)=tmelt
      END IF
    END IF

  END DO
#endif

  RETURN
END SUBROUTINE clsst
