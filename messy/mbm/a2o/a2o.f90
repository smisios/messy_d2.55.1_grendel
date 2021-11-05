! ******************************************************************
! ------------------------------------------------------------------
PROGRAM A2O
! ------------------------------------------------------------------
! Author: Bastian Kern, MPICH, Mainz, July 2009 / September 2010
! Code parts from Patrick Joeckel, MPICH, Mainz, June 2002
! ******************************************************************

  USE mo_f2kcli                    ! command line interface
  USE messy_a2o_gridtrafo          ! regridding interface

  IMPLICIT NONE

  INTRINSIC :: TRIM
 
  ! FOR COMMAND LINE
  CHARACTER(LEN=256) :: EXE          ! program name
  CHARACTER(LEN=80)  :: CMD          ! argument
  INTEGER            :: NARG         ! number of arguments

  NARG = COMMAND_ARGUMENT_COUNT()    ! number of arguments
  CALL GET_COMMAND_ARGUMENT(0,EXE)   ! program name

  IF (NARG > 1) THEN 
     WRITE(*,*) 'Too many arguments !'
     CALL USAGE(TRIM(EXE)) 
     STOP
  END IF

  IF (NARG == 0) THEN 
     CALL USAGE(TRIM(EXE)) 
     STOP
  END IF

  CALL GET_COMMAND_ARGUMENT(1,CMD)  

  ! REGRIDDING ...
  GT_CTRL = GT_PROC
  GT_NML  = NML_NEXT
  DO ! endless DO loop (must be terminated with EXIT)
     CALL GRIDTRAFO_CONTROL(GT_CTRL, GT_NML, GT_STATUS, TRIM(CMD))
     IF (GT_STATUS == GTSTAT_STOP) EXIT ! leave endless DO loop
  END DO
  ! END OF REGRIDDING

CONTAINS

  SUBROUTINE USAGE(EXE)
    CHARACTER (LEN=*) :: EXE
    WRITE(*,*) '-----------------------------------------'
    WRITE(*,*) 'A2O Version ',A2OVERS
    WRITE(*,*) 'Author: Bastian Kern, MPICH'
    WRITE(*,*) 'July 2009 / September 2010'
    WRITE(*,*) '-----------------------------------------'
    WRITE(*,*) 'Usage: '//TRIM(EXE)//' <namelist-file>'
    WRITE(*,*) '-----------------------------------------'
  END SUBROUTINE USAGE

END PROGRAM A2O
! ------------------------------------------------------------------
