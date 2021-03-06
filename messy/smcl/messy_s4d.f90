! **********************************************************************
MODULE messy_s4d
! **********************************************************************

  ! Sample in 4 Dimensions
  ! MODULE FOR HF-OUTPUT OF TRACERS/CHANNEL OBJECTS
  ! ALONG TRACKS
  !
  ! CORE MODULE (MESSy/SMCL)
  !
  ! Author: Patrick Joeckel, MPICH, January 2006
  !

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: ASSOCIATED, TRIM

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 's4d'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '2.0.2'
  ! MAXIMUM LENGTH OF STRING FOR CHANNEL(S)/OBJECT(S)
  INTEGER,         PARAMETER, PUBLIC :: STRLEN = 500

  TYPE T_POSITION
     INTEGER  :: yr, mo, dy, ho, mi, se
     REAL(DP) :: lon, lat, pre
  END TYPE T_POSITION

  PUBLIC :: T_POSITION
  PUBLIC :: read_track

CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE read_track(status, fname, iou, track)

    IMPLICIT NONE

    ! I/O
    INTEGER,                       INTENT(OUT) :: status 
    CHARACTER(LEN=*),              INTENT(IN)  :: fname  ! filename
    INTEGER,                       INTENT(IN)  :: iou    ! I/O unit
    TYPE(T_POSITION), DIMENSION(:), POINTER    :: track  ! track data

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_track'
    LOGICAL :: lex   ! file exists ?
    LOGICAL :: lopn  ! file open ?
    INTEGER :: iout  ! unit ?
    INTEGER :: n     ! number of data lines
    INTEGER :: mstat ! memory status
    INTEGER :: fstat ! file status
    INTEGER :: i     ! counter

    status = 1 ! ERROR

    ! CHECK, IF FILE IS PRESENT
    INQUIRE(file=TRIM(fname), exist=lex, opened=lopn, number=iout)
    IF (.NOT. lex) THEN
       WRITE(*,*) substr,': *** ERROR: File '''//TRIM(fname)//''' not found !'
       RETURN
    END IF
    IF (lopn) THEN
       WRITE(*,*) substr,': *** ERROR: File  '''//TRIM(fname)//&
            &''' already open on unit ',iout
       RETURN ! ERROR
    END IF

    ! OPEN FILE
    OPEN(unit=iou,file=TRIM(fname))
    ! COUNT LINES
    n = 0
    DO
       READ(iou, *, IOSTAT=fstat)
       IF (fstat < 0) EXIT
       n = n + 1
    END DO
    CLOSE(iou)

    ! ALLOCATE DATA
    IF (ASSOCIATED(track)) THEN
       DEALLOCATE(track)
       NULLIFY(track)
    END IF
    ALLOCATE(track(n), STAT=mstat)
    IF (mstat /= 0) THEN
       WRITE(*,*) '*** ERROR: Memory allocation failed !'
       RETURN ! ERROR
    END IF

    ! READ DATA
    OPEN(unit=iou,file=TRIM(fname))
    DO i=1, n

       READ(iou, *, IOSTAT=fstat) &
            track(i)%yr,  track(i)%mo,  track(i)%dy,   &
            track(i)%ho,  track(i)%mi,  track(i)%se,   &
            track(i)%lon, track(i)%lat, track(i)%pre

       IF (fstat /= 0) THEN
          WRITE(*,*) substr,': *** READ ERROR in line ',i,' of file '''&
               &//TRIM(fname)//''''
          RETURN ! ERROR
       END IF

    END DO

    ! END
    CLOSE(iou)
    status = 0

  END SUBROUTINE read_track
! ----------------------------------------------------------------------

! **********************************************************************
END MODULE messy_s4d
! **********************************************************************
