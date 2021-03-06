! **********************************************************************
MODULE messy_timepos_e5
! **********************************************************************

! Note for conversion from _e5 to _si:
! ------------------------------------
! The algorithm searches and interpolates in the "global fields", i.e.
! after gather_gp. philon and philat are not available as 1d arrays
! for COSMO, due to the rotated, curvilinear grid. A simple solution
! could be to convert the input lat,lon to rotated rlat,rlon and
! then operate in the entire model domain, i.e. to search in betweem
! rphilat and rphilon.
!

  ! ECHAM5/MESSy
  USE messy_main_timer_event, ONLY: time_event, io_time_event
  USE messy_main_blather_bi,  ONLY: start_message_bi, end_message_bi
  USE messy_main_tools,       ONLY: PTR_1D_ARRAY, PTR_3D_ARRAY
  USE messy_main_channel,     ONLY: STRLEN_OBJECT, STRLEN_CHANNEL
  ! MESSy
  USE messy_timepos

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  TYPE T_OBJ_IO
     CHARACTER(LEN=STRLEN_CHANNEL) :: channel = ''
     CHARACTER(LEN=STRLEN_OBJECT)  :: object  = ''
     CHARACTER(LEN=8)              :: suffix  = ''
  END TYPE T_OBJ_IO

  TYPE T_OBJ
     TYPE(T_OBJ_IO) :: io
     REAL(DP), DIMENSION(:,:,:,:), POINTER :: field => NULL()
  END TYPE T_OBJ

  INTEGER, PARAMETER :: NMAXOBJ = 100
  INTEGER, SAVE      :: NOBJ

  ! CPL-NAMELIST
  TYPE(T_OBJ_IO), DIMENSION(NMAXOBJ), SAVE :: OBJ
  CHARACTER(LEN=200)                , SAVE :: INPUT_PATH = ''
  LOGICAL                           , SAVE :: L_DEBUG = .FALSE.

  ! TIMER
  TYPE(io_time_event), SAVE :: TIMER = io_time_event(1, 'months','last',0)
  TYPE(time_event),    SAVE :: XTIMER 
  LOGICAL,             SAVE :: LTRIG_MONTHLY

  ! WORKSPACE
  TYPE(T_OBJ),     DIMENSION(NMAXOBJ), SAVE     :: XOBJ
  TYPE(T_TIMEPOS), DIMENSION(:),       POINTER  :: LIST
  LOGICAL, SAVE :: l_nofile = .TRUE.    ! file not available
  INTEGER, SAVE :: NMAXPOS = 0          ! max. position in list
  INTEGER, SAVE :: NPOS    = 1          ! current position in list
  
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: press_4d => NULL()

  PUBLIC :: timepos_initialize
  PUBLIC :: timepos_init_coupling
  PUBLIC :: timepos_global_start
  PUBLIC :: timepos_write_output
  PUBLIC :: timepos_free_memory
  !PRIVATE :: timepos_read_nml_cpl
  !PRIVATE :: read_timepos_data_file

CONTAINS

! ======================================================================
! PUBLIC ROUTINES
! ======================================================================

! ----------------------------------------------------------------------
  SUBROUTINE timepos_initialize

    ! TIMEPOS MODULE ROUTINE (ECHAM5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2006

    USE messy_main_timer_bi,   ONLY: timer_event_init
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'timepos_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL timepos_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('',substr)
    END IF
    CALL p_bcast(INPUT_PATH, p_io)
    CALL p_bcast(L_DEBUG, p_io)

    CALL start_message_bi(modstr,'OUTPUT INITIALISATION ',substr)

    IF (p_parallel_io) THEN
       ! GET NUMBER OF OBJECTS
       NOBJ = 1
       ! COPY DATA AND PARSE STR
       DO i=1, NMAXOBJ
          IF (TRIM(OBJ(i)%channel) == '') CYCLE
          IF (TRIM(OBJ(i)%object)  == '') CYCLE
          XOBJ(NOBJ)%io%channel = TRIM(OBJ(i)%channel)
          XOBJ(NOBJ)%io%object  = TRIM(OBJ(i)%object)
          XOBJ(NOBJ)%io%suffix  = TRIM(OBJ(i)%suffix)

          WRITE(*,*) '... ',TRIM(XOBJ(NOBJ)%io%channel),' - ',&
               TRIM(XOBJ(NOBJ)%io%object),' -> ',XOBJ(NOBJ)%io%suffix

          ! NEXT OBJECT
          NOBJ = NOBJ + 1
       END DO
       NOBJ = NOBJ - 1 
    END IF
    CALL p_bcast(NOBJ, p_io)

    ! BROADCAST RESULTS
    DO i=1, NOBJ
       CALL p_bcast(XOBJ(i)%io%channel, p_io)
       CALL p_bcast(XOBJ(i)%io%object,  p_io)
       CALL p_bcast(XOBJ(i)%io%suffix,  p_io)
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NOBJ,' OBJECT(S) REQUESTED !'
    END IF

    CALL timer_event_init(XTIMER, TIMER, 'timepos', 'present')

    CALL end_message_bi(modstr,'OUTPUT INITIALISATION ',substr)

  END SUBROUTINE timepos_initialize
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE timepos_init_coupling

    ! ECHAM5
    USE messy_main_blather_bi,       ONLY: info_bi
    ! ECHAM5/MESSy
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy
    USE messy_main_channel,    ONLY: get_channel_object &
                                   , get_channel_object_info

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'timepos_init_coupling'
    INTEGER :: i
    INTEGER :: reprid
    INTEGER :: status

    CALL start_message_bi(modstr,'COUPLING INITIALISATION ',substr)

    CALL get_channel_object(status &
         , 'ECHAM5', 'press', p4=press_4d )
    CALL channel_halt(substr, status)

    DO i=1, NOBJ

       CALL info_bi('   ... trying '//TRIM(XOBJ(i)%io%channel)//&
            &' - '//TRIM(XOBJ(i)%io%object)//' ... ', substr)
       CALL get_channel_object(status &
                  , TRIM(XOBJ(i)%io%channel), TRIM(XOBJ(i)%io%object) &
                  , p4=XOBJ(i)%field )
       IF (status /= 0) THEN
          CALL info_bi( &
               '            ... channel object not found ... skipping' &
               , substr)
          NULLIFY(XOBJ(i)%field)
          CYCLE
       END IF

       CALL get_channel_object_info(status &
            , TRIM(XOBJ(i)%io%channel), TRIM(XOBJ(i)%io%object) &
            , reprid=reprid)
       CALL channel_halt(substr, status)
       IF (reprid /= GP_3D_MID) THEN
          CALL info_bi( &
               '            ... representation not supported ... skipping' &
               , substr)
          NULLIFY(XOBJ(i)%field)
          CYCLE
       END IF

    END DO

    CALL end_message_bi(modstr,'COUPLING INITIALISATION ',substr)

  END SUBROUTINE timepos_init_coupling
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE timepos_global_start

    USE messy_main_timer_bi,      ONLY: event_state
    USE messy_main_timer,         ONLY: current_date &
                                      , lstart, lresume, YEAR, MONTH

    IMPLICIT NONE
    
    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr = 'timepos_global_start'

    IF (lstart .OR. lresume) THEN
       ! READ DATA FOR CURRENT DATE
       CALL read_timepos_data_file(YEAR, MONTH, .TRUE.)
    END IF

    ! CHECK EVENT
    LTRIG_MONTHLY = event_state(XTIMER, current_date)

  END SUBROUTINE timepos_global_start
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE timepos_write_output

    ! ECHAM5
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast &
                                   , dcg, gather_gp
    USE messy_main_blather_bi, ONLY: error_bi
    ! ECHAM5/MESSy
    USE messy_main_timer,      ONLY: YEAR,MONTH,DAY,HOUR,MINUTE,SECOND &
                                   , YEAR_NEXT,MONTH_NEXT              &
                                   , delta_time
    USE messy_main_grid_def_bi,    ONLY: philon, philat
    USE messy_main_grid_def_mem_bi,ONLY: nlon, nlev, ngl
    ! MESSy
    USE messy_main_timer,      ONLY: julian_day
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    INTRINSIC :: TRIM, ADJUSTL, ALLOCATED, ASSOCIATED, INT, REAL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'timepos_write_output'
    TYPE(t_ipos), DIMENSION(:), ALLOCATABLE :: IPOS
    INTEGER                     :: status
    REAL(dp)                    :: sjdf  ! julian day + fraction
    REAL(dp)                    :: dt    ! model time step (in days)
    REAL(dp)                    :: tjdf  ! julian day + fraction
    CHARACTER(LEN=26)           :: fname ! filename
    CHARACTER(LEN=6)            :: ymstr
    LOGICAL                     :: lex
    INTEGER                     :: iou
    INTEGER                     :: nstart, nstop
    LOGICAL                     :: lcalc
    INTEGER                     :: i,j
    REAL(DP), DIMENSION(:,:,:,:), POINTER  :: zfield => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER  :: press  => NULL()
    REAL(DP)                    :: value
    INTEGER                     :: ip

    IF (l_nofile) RETURN

    ! INIT
    ! JULIAN DAY + FRACTION OF SIMULATION TIME
    sjdf = julian_day(REAL(DAY, DP), MONTH, YEAR) &
         + REAL(HOUR,dp)/24.0_dp        &
         + REAL(MINUTE,dp)/1440.0_dp    &
         + REAL(SECOND,dp)/86400.0_dp
    dt = delta_time/86400.0_dp

    ! CHECK TIME
    nstart = NPOS
    nstart_loop: DO 
       IF (nstart > NMAXPOS) THEN
          ! END OF TIMEPOS LIST REACHED
          lcalc = .FALSE.
          EXIT
       END IF
       ! JULIAN DAY + FRACTION OF TRACK POSITION
       IF (p_parallel_io) THEN
          tjdf = julian_day( &
               REAL(list(nstart)%dy, dp)               &
               , list(nstart)%mo                       &
               , list(nstart)%yr )                     &
               + REAL(list(nstart)%ho,dp)/24.0_dp      &
               + REAL(list(nstart)%mi,dp)/1440.0_dp    &
               + REAL(list(nstart)%se,dp)/86400.0_dp
       END IF
       CALL p_bcast(tjdf, p_io)
       !
       IF (tjdf >= sjdf+dt) THEN
          ! TIMEPOS IS IN FUTURE
          lcalc = .FALSE.
          EXIT
       END IF
       IF (tjdf < sjdf) THEN
          ! TIMEPOS IS IN PAST: TRY NEXT POSITION
          nstart = nstart + 1 ! update counter
          CYCLE
       END IF
       ! NOW: sjdf <= tjdf < sjdf+dt, i.e.,
       ! TIMEPOS nstart IS BETWEEN CURRENT AND NEXT SIMULATION TIME STEP
       lcalc = .TRUE.
       EXIT
    END DO nstart_loop

    IF (lcalc) THEN

    nstop = nstart
    nstop_loop: DO 
       IF (nstop > NMAXPOS) THEN
          nstop = nstop - 1
          EXIT
       END IF
       ! CHECK TIME
       ! JULIAN DAY + FRACTION OF TRACK POSITION
       IF (p_parallel_io) THEN
          tjdf = julian_day( &
               REAL(list(nstop)%dy, dp)               &
               , list(nstop)%mo                       &
               , list(nstop)%yr )                     &
               + REAL(list(nstop)%ho,dp)/24.0_dp      &
               + REAL(list(nstop)%mi,dp)/1440.0_dp    &
               + REAL(list(nstop)%se,dp)/86400.0_dp
       END IF
       CALL p_bcast(tjdf, p_io)
       IF (tjdf >= sjdf+dt) THEN
          ! TIMEPOS IS IN FUTURE
          nstop = nstop - 1
          EXIT
       END IF
       IF (tjdf < sjdf+dt) THEN
          ! TIMEPOS IS IN PAST: TRY NEXT POSITION
          nstop = nstop + 1 ! update counter
       END IF

    END DO nstop_loop
    ! UPDATE POSITION
    NPOS = nstop + 1

    IF (p_parallel_io) WRITE(*,*) substr,': indices ',nstart,' to ',nstop

    ! INTERPOLATION WEIGHTS
    IF (p_parallel_io) THEN
       ALLOCATE(press(nlon, nlev, ngl,1))
    ELSE
       NULLIFY(press)
    END IF
    CALL gather_gp(press, press_4d, dcg)
    !
    IF (p_parallel_io) THEN
       ALLOCATE(IPOS(nstart:nstop))
       DO j=nstart, nstop
          CALL set_hweight(status, philon, philat, list(j), IPOS(j))
          IF (status /= 0) CALL error_bi( &
               ' set_hweight reported an error', substr)
          CALL set_vweight(status, press(:,:,:,1), list(j), IPOS(j))
          IF (status /= 0) CALL error_bi( &
               ' set_vweight reported an error', substr)
       END DO
    END IF
    !
    IF (ASSOCIATED(press)) THEN
       DEALLOCATE(press)
       NULLIFY(press)
    END IF

    IF (l_debug) THEN
       IF (p_parallel_io) THEN
          iou = find_next_free_unit(100,200)
          WRITE(ymstr(1:6),'(i4,i2.2)') YEAR, MONTH
          fname=modstr//'_'//TRIM(ymstr)//'_debug'//'.txt'
          INQUIRE(file=TRIM(fname), exist=lex)
          IF (lex) THEN
             OPEN(unit=iou, file=TRIM(fname) &
                  , status='OLD', position='APPEND')
          ELSE
             OPEN(unit=iou, file=TRIM(fname) &
                  , status='NEW')
             WRITE(iou,'(a45)') 'YYYYMMDD hhmm  lon(deg) lat(deg)  p(Pa)    ID'
          END IF
          DO j=nstart, nstop
             CALL debug_output(iou, list(j), IPOS(j))
          END DO
          CLOSE(iou)
       END IF
    END IF

    ! MEMORY
    IF (p_parallel_io) THEN
       ALLOCATE(zfield(nlon, nlev, ngl, 1))
    ELSE
       NULLIFY(zfield)
    END IF

    object_loop: DO i=1, NOBJ

       IF (.NOT. ASSOCIATED(xobj(i)%field)) CYCLE

       ! GATHER FIELD ON I/O PE
       IF (p_parallel_io) zfield(:,:,:,:) = 0.0_dp
       CALL gather_gp(zfield, xobj(i)%field, dcg)

       io_pe: IF (p_parallel_io) THEN
          ! OUTPUT FILE (ASCII)
          WRITE(ymstr(1:6),'(i4,i2.2)') YEAR, MONTH
          fname=modstr//'_'//&
               &TRIM(ymstr)//TRIM(ADJUSTL(XOBJ(i)%io%suffix))//'.txt'
          INQUIRE(file=TRIM(fname), exist=lex)
          iou = find_next_free_unit(100,200)
          IF (lex) THEN
             OPEN(unit=iou, file=TRIM(fname) &
                  , status='OLD', position='APPEND')
          ELSE
             OPEN(unit=iou, file=TRIM(fname) &
                  , status='NEW')
             WRITE(iou,'(a45)') 'YYYYMMDD hhmm  lon(deg) lat(deg)  p(Pa)    ID'
          END IF

          timepos_loop: DO j=nstart, nstop
          
             ! INTERPOLATE
             CALL interpolate(status, zfield(:,:,:,1), IPOS(j), value)
             IF (status /= 0) CALL error_bi( &
                  ' interpolate reported an error', substr)

             ! OUTPUT TO ASCII FILE
             ip = INT(list(j)%pre)
             WRITE(iou,'('//fstr//',e12.4)') &
                  YEAR, MONTH, DAY, HOUR, MINUTE &
                  , list(j)%lon, list(j)%lat, ip, list(j)%id, value

          END DO timepos_loop

          ! CLOSE OUTPUT ASCII FILE
          CLOSE(iou)
       END IF io_pe

    END DO object_loop

    ! CLEAN MEMORY
    IF (ASSOCIATED(zfield)) THEN
       DEALLOCATE(zfield)
       NULLIFY(zfield)
    END IF   
    !
    IF (ALLOCATED(IPOS)) DEALLOCATE(IPOS)

    END IF

    ! UPDATE DATA
    ! EVENT TRIGGERED, IF current_date IS 'last step of months'
    ! -> USE NEXT DATE FOR INPUT FILE NAMES
    CALL read_timepos_data_file(YEAR_NEXT, MONTH_NEXT, LTRIG_MONTHLY)

  END SUBROUTINE timepos_write_output
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE timepos_free_memory

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    IF (ASSOCIATED(list)) THEN
       DEALLOCATE(list)
       NULLIFY(list)
    END IF
    
  END SUBROUTINE timepos_free_memory
! ----------------------------------------------------------------------

! ======================================================================
! PRIVATE ROUTINES
! ======================================================================

! ----------------------------------------------------------------------
  SUBROUTINE timepos_read_nml_cpl(status, iou)

    ! TIMEPOS MODULE ROUTINE (ECHAM5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2006

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'timepos_read_nml_cpl'

    NAMELIST /CPL/ INPUT_PATH, L_DEBUG, OBJ

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    ! NOTE: already at definition

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE timepos_read_nml_cpl
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE read_timepos_data_file(YEAR, MONTH, lread_now)

    ! ECHAM5
    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,    ONLY: error_bi
    ! ECHAM5/MESSy
    USE messy_main_tools,         ONLY: find_next_free_unit

    IMPLICIT NONE
    
    INTRINSIC :: LEN_TRIM, ADJUSTL, TRIM, SIZE

    ! I/O
    INTEGER, INTENT(IN) :: YEAR
    INTEGER, INTENT(IN) :: MONTH
    LOGICAL, INTENT(IN) :: lread_now

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_timepos_data_file'
    CHARACTER(LEN=250) :: filename = ''
    INTEGER            :: idx
    INTEGER            :: iou
    INTEGER            :: status

    IF (.NOT. lread_now) RETURN

    filename = TRIM(ADJUSTL(INPUT_PATH))
    idx = LEN_TRIM(ADJUSTL(INPUT_PATH))
    IF (filename(idx:idx) == '/') THEN
       idx = idx + 1
    ELSE
       filename = TRIM(filename)//'/'
       idx = idx + 2
    END IF
    WRITE(filename(idx:idx+14),'(i4,i2.2,a8)') YEAR, MONTH, '_all.txt'

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       WRITE(*,*) substr,': READING ',TRIM(filename),' on unit ',iou
       CALL  read_timepos(status, TRIM(filename), iou, list)
       IF (status > 0) CALL error_bi('',substr)
       IF (status == -1) THEN
          l_nofile = .TRUE.  ! no file available
       ELSE
          l_nofile = .FALSE.
       END IF
       ! NOTE: list is only defined on I/O PE
       NMAXPOS = SIZE(list)
       WRITE(*,*) substr,': ',NMAXPOS,' TIME/POSITIONS READ'
    END IF

    CALL p_bcast(l_nofile, p_io)
    CALL p_bcast(NMAXPOS, p_io)
    
    ! RESET POSITION COUNTER
    NPOS = 1

  END SUBROUTINE read_timepos_data_file
! ----------------------------------------------------------------------

! **********************************************************************
END MODULE messy_timepos_e5
! **********************************************************************
