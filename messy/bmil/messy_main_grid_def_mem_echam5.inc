! -*- f90 -*-
 ! - GRID CONTROL
  USE mo_control,       ONLY: ngl, nlev, nlevp1, nlon     &  ! INTEGER
                            , nhgl                        &  ! INTEGER
                            , nsp                         &  ! INTEGER
                            , nmp1                        &  ! INTEGER
                            , nvclev, nn, nm, nk          &  ! INTEGER
                            , vct                         &  ! REAL(:)
                            , lmidatm
  USE mo_hyb,           ONLY: apzero, apsurf              &  ! REAL
                            , NPLVP1, NLMSGL, nlevm1      &  ! INTEGER
                            , NPLVP2, NLMSLP              &  ! INTEGER
                            , nplev, delb, delpr

  IMPLICIT NONE
  PUBLIC
  PRIVATE :: dp

  INTEGER, DIMENSION(3) :: nrow = (/0,0,0/)

  INTEGER :: nllev
  INTEGER :: nllevp1
  INTEGER :: lnsp             ! LS:
  INTEGER :: snsp             ! SP:
  INTEGER :: nlm              ! FAS:
  INTEGER :: jrow             ! GP: CURRENT ROW
  INTEGER :: jglat            ! GP: global continuous latitude index
  INTEGER :: kproma           ! GP: VECTOR LENGTH OF CURRENT ROW
  INTEGER :: nproma          ! GP: VECTOR LENGTH
  INTEGER :: npromz           ! GP: VECTOR LENGTH OF LAST ROW
  INTEGER :: ngpblks          ! GP: NUMBER OF ROWS
  INTEGER :: nglon, nglat     ! GP: NUMBER OF LON AND LAT (1 if lcolumn=T)
