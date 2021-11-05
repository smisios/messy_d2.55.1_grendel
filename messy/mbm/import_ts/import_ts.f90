PROGRAM ximport_ts

  USE messy_main_constants_mem,  ONLY: DP
  USE messy_main_timer,          ONLY: JulianMonthLength
  USE messy_main_import_ts

  IMPLICIT NONE
  INTRINSIC :: NULL, TRIM

  CHARACTER(LEN=*), PARAMETER :: modstr = 'import'
  CHARACTER(LEN=*), PARAMETER :: substr = 'ximport_ts'
  INTEGER :: status
  INTEGER :: iou1 = 101
  INTEGER :: ubase = 1000
  CHARACTER(LEN=2) :: str = ''
  INTEGER :: i, j, n
  ! TIME LOOP
  INTEGER :: year_start = 2000
  INTEGER :: year_end   = 2000
  INTEGER :: dpm
  INTEGER :: year, month, day, hr
  REAL(DP), DIMENSION(:), POINTER :: dat => NULL()


  ! --- INITIALISATION PHASE --------------------------------------------------

  CALL  import_ts_read_nml_ctrl(status, iou1)
  IF (status /= 0) STOP 'ERROR IN import_ts_read_nml_ctrl'

  NTS = 1
  DO i=1, NMAXTS

     IF (TRIM(TS(i)%name) == '') CYCLE

     CALL its_copy_io(XTS(NTS)%io, TS(i))

     WRITE(*,*) 'TIME SERIES   : ', TRIM(XTS(NTS)%io%name)
     WRITE(*,*) '   FILE       : ', TRIM(ADJUSTL(XTS(NTS)%io%fname))
     WRITE(*,*) '   VALID RANGE: ', XTS(NTS)%io%vr(:)
     
     DO j=1,2
        SELECT CASE(XTS(NTS)%io%cnt(j))
        CASE(TS_BD_STOP)
           WRITE(*,*) '   BND POLICY : stop'
        CASE(TS_BD_CONT)
           WRITE(*,*) '   BND POLICY : continue'
        CASE DEFAULT
!           CALL error_bi('UNKNOWN BOUNDARY POLICY',substr)
          STOP 'ERRROR: UNKNOWN BOUNDARY POLICY'
        END SELECT
     END DO

     SELECT CASE(XTS(NTS)%io%im)
     CASE(TS_IM_PREV)
        WRITE(*,*) '   SELECTION  : previous'
     CASE(TS_IM_NEXT)
        WRITE(*,*) '   SELECTION  : next'
     CASE(TS_IM_LINT)
        WRITE(*,*) '   SELECTION  : linear interpolation'
     CASE DEFAULT
!        CALL error_bi('UNKNOWN INTERPOLATION METHOD',substr)
        STOP 'ERROR: UNKNOWN INTERPOLATION METHOD'
     END SELECT

     WRITE(*,*) '   READING DATA ...'
     CALL its_read_ts(status, XTS(NTS), .TRUE.)
     IF (status /= 0) &
!          CALL error_bi('its_read_ts reported an error' ,substr)
          STOP 'ERROR: its_read_ts reported an error'
     WRITE(*,*) '   ... DONE!'
     
     ! NEXT TIME SERIES
     NTS = NTS + 1
     WRITE(*,*) '------------------------------------------------------'
  END DO
  NTS = NTS - 1

  DO i=1, NTS
     WRITE(str,'(i2.2)') i
     OPEN(unit=ubase+i,file='output_'//str//'.dat')
  END DO

  ! --- TIME INTEGRATION PHASE ------------------------------------------------

  WRITE(*,*) '======================================================'
  WRITE(*,*) 'TIME INTEGRATION STARTS ...'

  DO year = year_start, year_end
     DO month = 1, 12
        dpm = JulianMonthLength(year, month)
        DO day = 1, dpm
           DO hr=0,23

              DO i=1, NTS
                 CALL its_set_value_ts(status, XTS(i) &
                      , year, month, day, hr, 0, 0)
                 IF (status /= 0) &
!                      CALL error_bi('its_set_value_ts reported an error '&
!                      &'for time series'//TRIM(XTS(i)%io%name),substr)
                      STOP 'ERROR: its_set_value_ts reported an error '
                      
                 ALLOCATE(dat(XTS(i)%np))
                 dat(:) = XTS(i)%obj(:)
                 DO n=1, XTS(i)%np
                    IF (XTS(i)%flg(n) < 1.0_dp) dat(n) = -1.0E+34_dp
                 END DO
                 WRITE(ubase+i,*) year, month, day, hr, dat
                 DEALLOCATE(dat) ; NULLIFY(dat)

              END DO

           END DO
        END DO
     END DO
  END DO

  WRITE(*,*) ' ... DONE!'
  WRITE(*,*) '======================================================'

  ! --- FINILISING PHASE ------------------------------------------------------

  DO i=1, NTS
     WRITE(str,'(i2.2)') i
     CLOSE(unit=ubase+i)
  END DO

  DO i=1, NTS
     CALL its_delete_ts(XTS(i))
  END DO

  !----------------------------------------------------------------------------

END PROGRAM ximport_ts
