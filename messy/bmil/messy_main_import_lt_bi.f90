! ============================================================================
MODULE messy_main_import_lt_bi
! ============================================================================

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                    , error_bi, info_bi
  USE messy_main_import_lt

  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: import_lt_initialize
  PUBLIC :: import_lt_global_start
  PUBLIC :: import_lt_free_memory

CONTAINS

  ! ----------------------------------------------------------------------
  SUBROUTINE import_lt_initialize

    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,        ONLY: find_next_free_unit
    USE messy_main_timer,        ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

    IMPLICIT NONE
    INTRINSIC :: TRIM, ADJUSTL, LBOUND, UBOUND

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'import_lt_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i, j, n1, n2

    CALL start_message_bi(submodstr, 'INITIALIZE',substr)

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL import_lt_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    IF (p_parallel_io) THEN
       NLT = 1
       DO i=1, NMAXLT

          IF (TRIM(LT(i)%tname) == '') CYCLE
          IF (TRIM(LT(i)%fname) == '') CYCLE
          IF (TRIM(LT(i)%vname) == '') CYCLE

          XLT(NLT)%io%tname = TRIM(ADJUSTL(LT(i)%tname))
          XLT(NLT)%io%fname = TRIM(ADJUSTL(LT(i)%fname))
          XLT(NLT)%io%vname = TRIM(ADJUSTL(LT(i)%vname))


          WRITE(*,*) 'LOOKUP TABLE  : ', TRIM(ADJUSTL(XLT(NLT)%io%tname))
          WRITE(*,*) '   FILE       : ', TRIM(ADJUSTL(XLT(NLT)%io%fname))
          WRITE(*,*) '   VARIABLE   : ', TRIM(ADJUSTL(XLT(NLT)%io%vname))

          ! NEXT TIME SERIES
          NLT = NLT + 1
          WRITE(*,*) '------------------------------------------------------'
       END DO
       NLT = NLT - 1
    END IF
    CALL p_bcast(NLT, p_io)

    ! BROADCAST USER INPUT
    DO i=1, NLT
       CALL p_bcast(XLT(i)%io%tname, p_io)
       CALL p_bcast(XLT(i)%io%fname, p_io)
       CALL p_bcast(XLT(i)%io%vname, p_io)
    END DO

    ! IMPORT LOOKUP TABLES
    DO i=1, NLT

       ! ALLOCATE AND INIT
       ALLOCATE(XLT(i)%dimlen(NRANKMAX))
       XLT(i)%dimlen(:) = 1
       ALLOCATE(XLT(i)%dimname(NRANKMAX))
       XLT(i)%dimname(:) = ''
       ALLOCATE(XLT(i)%dimunit(NRANKMAX))
       XLT(i)%dimunit(:) = ''

       CALL ilt_set_file(status, i &
            , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, .TRUE.)
       IF (status /= 0) &
            CALL error_bi('ERROR IN LOOKUP TABLE FILENAME: '//&
            &TRIM(XLT(i)%io%fname), substr)

       IF (p_parallel_io) CALL ilt_read_lt_netcdf(status, i)
       CALL p_bcast(status, p_io)
       IF (status /= 0) CALL error_bi('ERROR IN LOOKUP TABLE IMPORT', substr)
       IF (p_parallel_io) &
            WRITE(*,*) '------------------------------------------------------'

       CALL p_bcast(XLT(i)%rank, p_io)
       CALL p_bcast(XLT(i)%dimlen, p_io)
       CALL p_bcast(XLT(i)%dimname, p_io)
       CALL p_bcast(XLT(i)%dimunit, p_io)
    END DO

    ! BROADCAST RESULTS
    DO i=1, NLT

       ! ALLOCATE MEMORY ON NON-IO PEs
       IF (.NOT. p_parallel_io) THEN
          ALLOCATE(XLT(i)%dimvar(XLT(i)%rank))
          DO j=1, XLT(i)%rank
             ALLOCATE(XLT(i)%dimvar(j)%ptr(XLT(i)%dimlen(j)))
          END DO
          ALLOCATE( XLT(i)%data( &
               XLT(i)%dimlen(1), &
               XLT(i)%dimlen(2), &
               XLT(i)%dimlen(3), &
               XLT(i)%dimlen(4), &
               XLT(i)%dimlen(5), &
               XLT(i)%dimlen(6)  ) )
       END IF

       ! BROADCAST CONTENTS
       DO j=1, XLT(i)%rank
          CALL p_bcast(XLT(i)%dimvar(j)%ptr(:), p_io)
       END DO
       DO n1 = LBOUND(XLT(i)%data, 5), UBOUND(XLT(i)%data, 5)
          DO n2 = LBOUND(XLT(i)%data, 6), UBOUND(XLT(i)%data, 6)
             CALL p_bcast(XLT(i)%data(:,:,:,:,n1,n2), p_io)
          END DO
       END DO
    END DO

    CALL end_message_bi(submodstr, 'INITIALIZE',substr)

  END SUBROUTINE import_lt_initialize
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE import_lt_global_start

    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_timer,        ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'import_lt_global_start'
    INTEGER :: status
    INTEGER :: i

    ! UPDATE LOOKUP TABLES
    DO i=1, NLT

       CALL ilt_set_file(status, i &
            , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, .FALSE.)
       IF (status /= 0) &
            CALL error_bi('ERROR IN LOOKUP TABLE FILENAME: '//&
            &TRIM(XLT(i)%io%fname), substr)

       IF (XLT(i)%lupdate) THEN
          CALL info_bi('START UPDATING LOOKUP TABLE',substr)
          IF (p_parallel_io) CALL ilt_read_lt_netcdf(status, i)
          CALL p_bcast(status, p_io)
          IF (status /= 0) CALL error_bi('ERROR IN LOOKUP TABLE UPDATE', substr)
          CALL p_bcast(XLT(i)%rank, p_io)
          CALL p_bcast(XLT(i)%dimlen, p_io)
          CALL p_bcast(XLT(i)%dimname, p_io)
          CALL p_bcast(XLT(i)%dimunit, p_io)
          CALL info_bi('END UPDATING LOOKUP TABLE',substr)
       END IF
    END DO

  END SUBROUTINE import_lt_global_start
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE import_lt_free_memory

    IMPLICIT NONE

    INTEGER :: i

    DO i=1, NLT
       CALL ilt_delete_lt(i)
    END DO

  END SUBROUTINE import_lt_free_memory
  ! ----------------------------------------------------------------------

! ============================================================================
END MODULE messy_main_import_lt_bi
! ============================================================================
