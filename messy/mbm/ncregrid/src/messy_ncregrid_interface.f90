! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_INTERFACE
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
! ******************************************************************

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: NCREGRID_SET_MESSAGEMODE
  PUBLIC :: INTERFACE_GEOHYBGRID

CONTAINS

! ------------------------------------------------------------------
SUBROUTINE NCREGRID_SET_MESSAGEMODE

  USE MESSY_NCREGRID_BASE, ONLY: MSGMODE &
                               , MSGMODE_S, MSGMODE_E, MSGMODE_VL  &
                               , MSGMODE_W, MSGMODE_VM, MSGMODE_I

  IMPLICIT NONE

  ! RESET MESSAGE MODUS
  MSGMODE = MSGMODE_S + MSGMODE_E + MSGMODE_VL  &
          + MSGMODE_W + MSGMODE_VM + MSGMODE_I

END SUBROUTINE NCREGRID_SET_MESSAGEMODE
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE INTERFACE_GEOHYBGRID(g, ok)

  USE MESSY_NCREGRID_GEOHYB,       ONLY: geohybgrid, init_geohybgrid

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid), INTENT(OUT) :: g
  LOGICAL          , INTENT(OUT) :: ok

  CALL INIT_GEOHYBGRID(g)
  ok = .false.

END SUBROUTINE INTERFACE_GEOHYBGRID
! ------------------------------------------------------------------

! ******************************************************************
! ------------------------------------------------------------------
END MODULE MESSY_NCREGRID_INTERFACE
! ------------------------------------------------------------------
! ******************************************************************
