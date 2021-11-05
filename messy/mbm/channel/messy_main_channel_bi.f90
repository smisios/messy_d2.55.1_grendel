! **********************************************************************
MODULE messy_main_channel_bi
! **********************************************************************

  ! -------------------------------------------------------------------
  !  CHANNEL INTERFACE FOR SIMPLIFIED BOX MODEL (BMIL)
  ! -------------------------------------------------------------------

  USE messy_main_channel_error_bi, ONLY: channel_halt
  ! ... base model parameters ...
  USE channel_mem_bml, basemodstr=>modstr, basemodver=>modver

  ! ... set the working precision (kind)
  USE messy_main_constants_mem, ONLY: DP

  ! ... for output messages ...
  USE messy_main_blather,  ONLY: start_message, end_message

  ! ... CHANNELS
  USE messy_main_channel,       ONLY: NMAXCHANNELS, STRLEN_CHANNEL, NCHANNEL &
                                    , modstr, IOMODE_OUT, IOMODE_RST         &
                                    , AMODE_WRITE, AMODE_READ, AMODE_INIT    &
                                    , REPR_UNDEF, DIMID_UNDEF
  USE messy_main_channel_attributes, ONLY: t_attribute_list

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE

  PUBLIC :: IOMODE_OUT, IOMODE_RST

  ! DIMENSION IDs
  ! - TIME
  INTEGER, SAVE, PUBLIC :: DIMID_TIME = DIMID_UNDEF   ! time
  ! - GRIDPOINT
  INTEGER, SAVE, PUBLIC :: DIMID_LEV  = DIMID_UNDEF   ! level
  INTEGER, SAVE, PUBLIC :: DIMID_LON  = DIMID_UNDEF   ! longitude
  INTEGER, SAVE, PUBLIC :: DIMID_LAT  = DIMID_UNDEF   ! latitude
  !
  ! - (ARBITRARY) ARRAY LENGTH
  INTEGER, SAVE, PUBLIC :: DIMID_AL   = DIMID_UNDEF
  !
  ! REPRESENTATION IDs
  ! - GRIDPOINT (longitude x latitude x level)
  INTEGER, SAVE, PUBLIC :: GP_3D_MID        = REPR_UNDEF
  ! mz_ab_20090921+
  ! - GRIDPOINT WITH 2 BOUNDARY BOXES (EACH SIDE) IN LON AND LAT
  INTEGER, SAVE, PUBLIC :: GP_3D_MID_BND    = REPR_UNDEF
  ! mz_ab_20090921-
  ! - SCALAR (zero dimension)
  INTEGER, SAVE, PUBLIC :: SCALAR           = REPR_UNDEF
  ! - ARRAY
  INTEGER, SAVE, PUBLIC :: ARRAY            = REPR_UNDEF

  ! OUTPUT TIME
  INTEGER,          DIMENSION(:), ALLOCATABLE         :: IOUTPUT_FREQ
  LOGICAL,          DIMENSION(:), ALLOCATABLE, PUBLIC :: LOUTPUT_NOW
  LOGICAL, SAVE                               :: LFORCE_NEW_OUTPUT = .FALSE.

  ! TRIGGER NEW FILE
  INTEGER,          DIMENSION(:), ALLOCATABLE         :: INEWFILE_FREQ
  LOGICAL,          DIMENSION(:), ALLOCATABLE, PUBLIC :: LTNF_NOW

  ! SPECIAL RESTART ATTRIBUTES
  TYPE(t_attribute_list), POINTER, SAVE :: restart_att => NULL()

  ! <EXP_NAME> + _NNNN_  = 14 + 6 (+1)
  INTEGER, PARAMETER :: BASENAME_LEN = 21

  ! ====================================================================
  ! FOR CPL-NAMELIST
  !
  TYPE t_channel_timer
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname = '' ! CHANNEL NAME
     INTEGER                       :: nstep = 1  ! DEFAULT: output every step
  END TYPE t_channel_timer
  !
  ! OUTPUT FREQUENCY
  TYPE(t_channel_timer),                          SAVE :: TIMER_DEFAULT
  TYPE(t_channel_timer), DIMENSION(NMAXCHANNELS), SAVE :: TIMER_CHANNEL
  !
  ! NEW FILE-NAME FREQUENCY
  TYPE(t_channel_timer),                          SAVE :: TIMER_TNF_DEFAULT
  TYPE(t_channel_timer), DIMENSION(NMAXCHANNELS), SAVE :: TIMER_TNF_CHANNEL
  ! ====================================================================

  ! exemplary variables for standard basemodel channel
  REAL(DP), POINTER :: yr,mo,dy,hr,mi,se,ms

  ! MAIN ENTRY POINTS (PUBLIC)
  ! ... initialisation phase ....
  PUBLIC   :: main_channel_setup
  !                         !-> channel_init_restart_bi
  !                         !   !-> initialize_restart_attributes
  PUBLIC   :: main_channel_initialize
  !                         |-> main_channel_read_nml_ctrl
  !                         |-> main_channel_read_nml_cpl
  !                         |-> main_channel_initialize_gatts
  !                         |-> main_channel_initialize_dims
  !                         |-> main_channel_initialize_reprs
  !
  PUBLIC   :: main_channel_init_memory
  !
  PUBLIC   :: main_channel_init_coupling
  !                         |-> fixate_channels
  !                         |-> main_channel_init_timer
  !                         |-> write_attribute
  !                         |-> write_dimension
  !                         |-> write_representation
  !                         |-> write_channel

  ! ... time loop ...
  PUBLIC   :: main_channel_write_output
  !                         |-> main_channel_update_vars
  !                         |-> main_channel_update_timer
  PUBLIC   :: main_channel_write_restart
  !                         !-> channel_write_output(IOMODE_RST)
  ! ... finalising phase ....
  PUBLIC   :: main_channel_free_memory
  !                         |-> channel_finish_io
  !                         |-> clean_representations
  !                         |-> clean_dimensions
  !                         |-> clean_channels
  !                         |-> write_channel

  ! ... input / output ...
  PUBLIC   :: channel_write_output
  !                         |-> update_channels(1)
  !                             ( ACCUMULATE 2ndary DATA )
  !                         |-> update_channels(2)
  !                             ( PREPARE FOR OUTPUT )
  !                         |-> update_channels(3)
  !                             ( RESET AFTER OUTPUT )
  !                         !-> initialize_restart_attributes
  PUBLIC   :: channel_init_restart_bi
  !                         !-> initialize_parallel_io
  !                         !-> initialize_restart_attributes
  PUBLIC   :: main_channel_read_restart
  !                         !-> initialize_restart_attributes

  !
  ! SHARED BETWEEN DIFFERENT BASEMODELS
  !PRIVATE :: initialize_restart_attributes  ! INIT RESTART ATTRIBUTES
  !PRIVATE :: main_channel_init_timer        ! INITIALIZE OUTPUT TIMERS
  !PRIVATE :: main_channel_update_timer
  !PRIVATE :: main_channel_read_nml_cpl
  !
  ! BASEMODEL SPECIFIC
  !PRIVATE :: main_channel_initialize_gatts  ! SET GLOBAL ATTRIBUTES
  !PRIVATE :: main_channel_initialize_dims   ! DEFINE DIMENSIONS
  !PRIVATE :: main_channel_initialize_reprs  ! DEFINE REPRESENTATIONS
  !PRIVATE :: main_channel_update_vars       ! update basemodel variables

CONTAINS

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES (MAIN ENTRY POINTS CALLED FROM 'CONTROL')
  ! NOTE: In the simple example, 'CONTROL' and the BMIL are merged
  !       in the file channel_bml.f90
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_setup

    CALL channel_init_restart_bi

  END SUBROUTINE main_channel_setup
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize

    USE messy_main_channel,  ONLY: main_channel_read_nml_ctrl    &
                                 , NMAXCHANNELS, NMAXOBJECTS     &
                                 , NMAXADDCHANNEL, NMAXADDREF    &
                                 , ADD_CHANNEL, ADD_REF, OUT_DEFAULT &
                                 , OUT_PREC &
                                 , OUT_CHANNEL, OUT_OBJECT, EXP_NAME &
                                 , L_FLUSH_IOBUFFER
    USE messy_main_tools,    ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize'
    INTEGER :: status
    INTEGER :: iou
    INTEGER :: i

    CALL start_message(modstr,'INITIALIZE CHANNELS',substr)

    ! READ CTRL-NAMELIST FOR FIXATION ( ... AT END OF INIT_COUPLING)
    ! INITIALIZE CTRL
    iou = find_next_free_unit(100,200)
    CALL main_channel_read_nml_ctrl(status, iou)
    IF (status /= 0) THEN
       WRITE(*,*) 'ERROR IN READING CTRL FROM '//modstr//'.nml'
       STOP
    END IF

    ! INTIALIZE GLOBAL ATTRIBUTES
    CALL main_channel_initialize_gatts

    ! INTIALIZE DIMENSIONS
    CALL main_channel_initialize_dims

    ! INTIALIZE REPRESENTATIONS
    CALL main_channel_initialize_reprs

    ! INITIALIZE CPL
    iou = find_next_free_unit(100,200)
    CALL main_channel_read_nml_cpl(status, iou)
    IF (status /= 0) THEN
       WRITE(*,*) 'ERROR IN READING CPL FROM'//modstr//'.nml'
       STOP
    END IF

    CALL end_message(modstr,'INITIALIZE CHANNELS',substr)

  END SUBROUTINE main_channel_initialize
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_init_memory

    USE messy_main_channel, ONLY: new_channel, new_channel_object &
                                , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_init_memory'
    INTEGER :: status

    CALL start_message(basemodstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

    CALL new_channel(status, basemodstr, reprid=SCALAR)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, basemodstr, 'channel_info' &
         , c = 'standard basemodel channel' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, basemodstr, 'yr', yr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, basemodstr, 'yr' &
         , 'long_name', c = 'year')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, basemodstr, 'mo', mo)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, basemodstr, 'mo' &
         , 'long_name', c = 'month')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, basemodstr, 'dy', dy)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, basemodstr, 'dy' &
         , 'long_name', c = 'day')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, basemodstr, 'hr', hr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, basemodstr, 'hr' &
         , 'long_name', c = 'hour')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, basemodstr, 'mi', mi)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, basemodstr, 'mi' &
         , 'long_name', c = 'minute')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, basemodstr, 'se', se)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, basemodstr, 'se' &
         , 'long_name', c = 'second')
    CALL channel_halt(substr, status)

    ! FORECE RESTART FILE CREATION
    CALL new_channel_object(status, basemodstr, 'ms', ms, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, basemodstr, 'ms' &
         , 'long_name', c = 'millisecond')
    CALL channel_halt(substr, status)

    CALL end_message(basemodstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

  END SUBROUTINE main_channel_init_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_init_coupling(flag)

    USE messy_main_tools,    ONLY: int2str
    USE messy_main_channel,  ONLY: fixate_channels, write_channel &
                                 , write_attribute, set_channel_output &
                                 , modify_attributes
    USE messy_main_channel_dimensions,  ONLY: write_dimension
    USE messy_main_channel_repr,        ONLY: write_representation

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_init_coupling'
    INTEGER                     :: status

    SELECT CASE(flag)

    CASE(1)

       CALL modify_attributes(status)
       CALL channel_halt(substr, status)

    CASE(2)

       CALL start_message(modstr,'FIXATE CHANNELS',substr)

       CALL fixate_channels(status, 1)
       CALL fixate_channels(status, 2)
       CALL channel_halt(substr, status)

       ! INITIALIZE CHANNEL OUTPUT TIMERS (via CPL-NAMELIST)
       CALL main_channel_init_timer

       ! ... diagnostic ouput ...
       CALL write_attribute(status)
       CALL channel_halt(substr, status)

       CALL write_dimension(status)
       CALL channel_halt(substr, status)

       CALL write_representation(status)
       CALL channel_halt(substr, status)

       CALL write_channel(status)
       CALL channel_halt(substr, status)

       CALL end_message(modstr,'FIXATE CHANNELS',substr)

    END SELECT

  END SUBROUTINE main_channel_init_coupling
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_write_output(flag)

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_write_output'

    SELECT CASE(flag)
    CASE(1)
       CALL main_channel_update_vars
       CALL main_channel_update_timer
    CASE(2)
       CALL channel_write_output(IOMODE_OUT)
    CASE DEFAULT
       WRITE(*,*) 'ERROR: unknown value of flag in call to '//substr
       STOP
    END SELECT

  END SUBROUTINE main_channel_write_output
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_write_restart

    IMPLICIT NONE

    CALL channel_write_output(IOMODE_RST)

  END SUBROUTINE main_channel_write_restart
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_free_memory

    USE messy_main_channel_io,          ONLY: channel_finish_io
    USE messy_main_channel_dimensions,  ONLY: clean_dimensions
    USE messy_main_channel_repr,        ONLY: clean_representations
    USE messy_main_channel,             ONLY: clean_channels, write_channel

    IMPLICIT NONE

    INTRINSIC :: ALLOCATED

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_free_memory'
    INTEGER :: status

    CALL start_message(modstr,'FINISH CHANNELS',substr)

    IF (ALLOCATED(LOUTPUT_NOW))   DEALLOCATE(LOUTPUT_NOW)
    IF (ALLOCATED(IOUTPUT_FREQ))  DEALLOCATE(IOUTPUT_FREQ)
    IF (ALLOCATED(LTNF_NOW))      DEALLOCATE(LTNF_NOW)
    IF (ALLOCATED(INEWFILE_FREQ)) DEALLOCATE(INEWFILE_FREQ)

    ! ------------------------------------------
    ! CLOSE ALL OUTPUT FILES
    ! ------------------------------------------
    CALL channel_finish_io(status, IOMODE_OUT, .TRUE., 0)
    CALL channel_halt(substr, status)

    CALL clean_representations(status)
    CALL channel_halt(substr, status)

    CALL clean_dimensions(status)
    CALL channel_halt(substr, status)

    CALL clean_channels(status)
    CALL channel_halt(substr, status)

    ! ... diagostic output ...
    CALL write_channel(status)
    CALL channel_halt(substr, status)

    CALL end_message(modstr,'FINISH CHANNELS',substr)

  END SUBROUTINE main_channel_free_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PUBLIC ENTRY POINTS CALLED DIRECTLY FROM THE BMIL
  ! NOTE: In the simple example, 'CONTROL' and the BMIL are merged
  !       in the file channel.f90
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_write_output(IOMODE)

    USE messy_main_channel_io, ONLY: channel_init_io           &
                                   , channel_write_header      &
                                   , channel_write_time        &
                                   , channel_write_data        &
                                   , channel_finish_io
    USE messy_main_channel,    ONLY: EXP_NAME, update_channels
    USE messy_main_channel_dimensions, ONLY: update_dimension_variable

    IMPLICIT NONE

    INTRINSIC :: ABS, ASSOCIATED, LEN, LEN_TRIM, NULL, REAL

    ! I/O
    INTEGER, INTENT(IN)         :: IOMODE  ! OUTPUT MODE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'channel_write_output'
    INTEGER                     :: status
    CHARACTER(LEN=BASENAME_LEN) :: fnamebase = ''
    INTEGER                     :: i
    REAL(DP)                    :: now
    LOGICAL                     :: lexit
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: ptr  => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: gptr => NULL()
    INTEGER                               :: reprid
    INTEGER, SAVE    :: nrstcount     = 0
    CHARACTER(LEN=4) :: nrstcount_str = ''
    LOGICAL          :: lp
    INTEGER          :: p_io_c

    CALL start_message(modstr,'WRITE OUTPUT',substr)

    nrstcount = cy

    SELECT CASE(IOMODE)
       !
    CASE(IOMODE_OUT)   ! ### ---------------- OUTPUT -------------------- ###
       !
       CALL update_channels(status, 1, time_step_len) ! ACCUMULATE 2ndary
       CALL channel_halt(substr, status)
       !
       CALL update_channels(status, 2, time_step_len) ! PREPARE FOR OUTPUT
       CALL channel_halt(substr, status)
       !
       ! UPDATE FILENAME BASE
       fnamebase = EXP_NAME
       DO i = LEN_TRIM(EXP_NAME)+1, LEN(EXP_NAME)
          WRITE(fnamebase(i:i),'(a1)') '_'
       END DO
       WRITE(fnamebase(LEN(EXP_NAME)+1:),'(a1,i4.4,a1)') '_', NH, '_'
       !
       ! UPDATE TIME INFORMATION
       now = REAL(NH, dp)
       !
       CALL update_dimension_variable(status, 'time', 'time', (/ now /))
       CALL channel_halt(substr, status)
       !
    CASE(IOMODE_RST)   ! ### ---------------- RESTART -------------------- ###
       !
       ! FORCE NEW OUTPUT FILE IN NEXT STEP
       LFORCE_NEW_OUTPUT = .TRUE.
       ! UPDATE FILENAME BASE
       nrstcount = nrstcount + 1
       write(nrstcount_str,'(i4.4)') nrstcount
       ! UPDATE FILENAME BASE
       fnamebase = 'restart_'//nrstcount_str//'_'
       !
       ! SET RESTART ATTRIBUTES
       CALL initialize_restart_attributes(AMODE_WRITE)
       !
    END SELECT

    ! PREPARE OUTPUT / RESTART FILE
    ! NEW FILE, SAVE I/O-UNITS, FILE-IDs etc.
    CALL channel_init_io(status, IOMODE, fnamebase, AMODE_WRITE, 0)
    CALL channel_halt(substr, status)
    !
    ! WRITE HEADER
    SELECT CASE(IOMODE)
    CASE(IOMODE_OUT)
       CALL channel_write_header(status, IOMODE, DIMID_TIME, 0)
    CASE(IOMODE_RST)
       CALL channel_write_header(status, IOMODE, DIMID_TIME, 0 &
            , restart_att)
    END SELECT
    CALL channel_halt(substr, status)
    !
    ! WRITE TIME STEP DATA TO OUTPUT FILE
    CALL channel_write_time(status, IOMODE, DIMID_TIME, 0)
    CALL channel_halt(substr, status)

    ! WRITE DATA
    DO
       ! NEXT POINTER
       CALL channel_write_data(status, lp &
            , IOMODE, lexit, ptr, reprid, p_io_c, 0)
       CALL channel_halt(substr, status)
       IF (lexit) EXIT

       gptr => ptr

       ! OUTPUT
       CALL channel_write_data(status, lp &
            , IOMODE, lexit, gptr, reprid, p_io_c, 0)
       CALL channel_halt(substr, status)
    END DO

    ! FLUSH (OUTPUT) / CLOSE (RESTART) BUFFER
    CALL channel_finish_io(status, IOMODE, (IOMODE==IOMODE_RST), 0)
    CALL channel_halt(substr, status)


    IF (IOMODE == IOMODE_OUT) THEN
       !
       CALL update_channels(status, 3, time_step_len) ! RESET AFTER OUTPUT
       CALL channel_halt(substr, status)
       !
    END IF

    CALL end_message(modstr,'WRITE OUTPUT',substr)

  END SUBROUTINE channel_write_output
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_init_restart_bi

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'messy_channel_init_restart'
    INTEGER :: status

    IF (.NOT. lresume) RETURN

    CALL initialize_restart_attributes(AMODE_INIT)

  END SUBROUTINE channel_init_restart_bi
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_read_restart

    USE messy_main_channel_io, ONLY: channel_init_io   &
                                   , channel_read_data &
                                   , channel_finish_io
    USE messy_main_channel,    ONLY: t_channel_list, t_channel, GCHANNELLIST

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, NULL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'main_channel_read_restart'
    INTEGER                      :: status
    CHARACTER(LEN=BASENAME_LEN)  :: fnamebase = ''
    TYPE(t_channel_list),  POINTER        :: ls
    TYPE(t_channel),       POINTER        :: channel
    LOGICAL                               :: lexit
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: gptr => NULL()
    INTEGER                               :: reprid
    LOGICAL                               :: lok
    LOGICAL                               :: lp

    CALL start_message(modstr,'READ RESTART',substr)

    ! UPDATE FILENAME BASE
    fnamebase = 'restart_'//cystr//'_'

    ! INITIALIZE RESTART ATTRIBUTES
    CALL initialize_restart_attributes(AMODE_READ)

    ! OPEN RESTART FILES AND CHECK HEADER INFORMATION
    ! OPEN FILE FOR READ
    CALL channel_init_io(status &
         , IOMODE_RST, fnamebase, AMODE_READ, 0, restart_att)
    CALL channel_halt(substr, status)

    ! READ DATA
    DO
       ! NEXT POINTER
       CALL channel_read_data(status, .TRUE. &
            , IOMODE_RST, lexit, gptr, reprid, lp)
       CALL channel_halt(substr, status)
       IF (lexit) EXIT

       ! mz_ab_20100120+
       ! DISTRIBUTE
       CALL channel_read_data(status, .TRUE., &
            IOMODE_RST, lexit, gptr, reprid, lp)
       CALL channel_halt(substr, status)
       ! mz_ab_20100120-

       IF (ASSOCIATED(gptr)) THEN
          DEALLOCATE(gptr)
          NULLIFY(gptr)
       END IF
    END DO

    ! CLOSE ALL RESTART FILES
    CALL channel_finish_io(status, IOMODE_RST, .TRUE., 0)
    CALL channel_halt(substr, status)

    CALL end_message(modstr,'READ RESTART',substr)

  END SUBROUTINE main_channel_read_restart
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize_gatts

    USE messy_main_channel,       ONLY: new_attribute, write_attribute &
                                      , AF_RST_CMP, AF_RST_INP, EXP_NAME &
                                      , modstr, modver

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, TRIM, ADJUSTL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize_gatts'
    INTEGER :: status
    CHARACTER (8)   :: ydate
    CHARACTER (10)  :: ytime

    ! example 1: system date and time
    CALL DATE_AND_TIME(ydate, ytime)
    CALL new_attribute(status, 'operating_date_time' &
         , c = ydate(1:8)//' '//ytime(1:6) )
    CALL channel_halt(substr, status)

    ! example 2: experiment name from namelist
    CALL new_attribute(status, 'MESSy_experiment', c=TRIM(EXP_NAME))
    CALL channel_halt(substr, status)

    ! example 3: basemodel version information
    CALL new_attribute(status, 'MESSy_basemodel' &
         , c = basemodstr//' version '//basemodver//&
         &', Max-Planck Institute for Chemistry, Mainz' )
    CALL channel_halt(substr, status)

    ! example 4: basemodel version information
    CALL new_attribute(status, 'MESSy_'//modstr &
         , c = ' version '//basemodver//&
         &', Max-Planck Institute for Chemistry, Mainz' )
    CALL channel_halt(substr, status)

    ! example 4: special flag to be checked from restart file
    CALL new_attribute(status, 'restart_test_flag' &
         , i=31415, iflag=AF_RST_CMP)
    CALL channel_halt(substr, status)

    ! diagnostic output
    CALL write_attribute(status)
    CALL channel_halt(substr, status)

  END SUBROUTINE main_channel_initialize_gatts
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize_dims

    USE messy_main_channel_dimensions,   ONLY: new_dimension            &
                                             , write_dimension          &
                                             , add_dimension_variable   &
                                             , add_dimension_variable_att

    IMPLICIT NONE
    INTRINSIC :: RANDOM_SEED

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize_dims'
    INTEGER                             :: status
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: array
    INTEGER                             :: i
    INTEGER                             :: nlev, nlon, nlat, n

    ! #############
    ! ### TIME ####
    ! #############
    ! ... dimension ...
    CALL new_dimension(status, DIMID_TIME, 'time', 1, .TRUE.)
    CALL channel_halt(substr, status)

    ! ... dimension variable 'time' ...
    CALL add_dimension_variable(status, 'time', 'time', (/ REAL(NH,DP) /))
    CALL channel_halt(substr, status)
    ! ... with attribute 'long_name' ...
    CALL add_dimension_variable_att(status, 'time', 'time', &
         'long_name', c='time')
    CALL channel_halt(substr, status)
    ! ... with attribute 'units' ...
    CALL add_dimension_variable_att(status, 'time', 'time', &
         'units', c='hour since 2009-01-01 00:00:00')
    CALL channel_halt(substr, status)

    ! ###############
    ! ### LEVELS ####
    ! ###############

    nlev = 2  ! number of levels
    !
    ! ... dimension ...
    CALL new_dimension(status, DIMID_LEV, 'lev', nlev)
    CALL channel_halt(substr, status)
    !
    ! ... dimension variable ...
    ALLOCATE(array(nlev))
    DO i=1, nlev
       array(i) = REAL(i, DP)
    END DO
    CALL add_dimension_variable(status, 'lev', 'lev', array)
    CALL channel_halt(substr, status)
    DEALLOCATE(array)

    ! ... with attribute 'long_name' ...
    CALL add_dimension_variable_att(status, 'lev', 'lev', &
         'long_name', c='level index')
    CALL channel_halt(substr, status)

    ! ... with attribute 'units' ...
    CALL add_dimension_variable_att(status, 'lev', 'lev', &
         'units', c='level')
    CALL channel_halt(substr, status)

    ! ###################
    ! ### LONGITUDES ####
    ! ###################

    nlon = 36  ! number of longitudes

    ! ... dimension ...
    CALL new_dimension(status, DIMID_LON, 'lon', nlon)
    CALL channel_halt(substr, status)

    ! ... dimension variable ...
    ALLOCATE(array(nlon))
    DO i=1, nlon
       array(i) = REAL((i-1)*10, DP) + 5.0_dp
    END DO
    CALL add_dimension_variable(status, 'lon', 'lon', array)
    CALL channel_halt(substr, status)
    DEALLOCATE(array)

    ! ... with attribute 'long_name' ...
    CALL add_dimension_variable_att(status, 'lon', 'lon', &
         'long_name', c='longitude')
    CALL channel_halt(substr, status)

    ! ... with attribute 'units' ...
    CALL add_dimension_variable_att(status, 'lon', 'lon', &
         'units', c='degrees_east')
    CALL channel_halt(substr, status)

    ! ##################
    ! ### LATITUDES ####
    ! ##################

    nlat = 18 ! number of latitudes

    ! ... dimension ...
    CALL new_dimension(status, DIMID_LAT, 'lat', nlat)
    CALL channel_halt(substr, status)

    ! ... dimension variable ...
    ALLOCATE(array(nlat))
    DO i=1, nlat
       array(i) = -85.0_dp + REAL((i-1)*10, DP)
    END DO
    CALL add_dimension_variable(status, 'lat', 'lat', array)
    CALL channel_halt(substr, status)
    DEALLOCATE(array)

    ! ... with attribute 'long_name' ...
    CALL add_dimension_variable_att(status, 'lat', 'lat', &
         'long_name', c='latitude')
    CALL channel_halt(substr, status)

    ! ... with attribute 'units' ...
    CALL add_dimension_variable_att(status, 'lat', 'lat', &
         'units', c='degrees_north')
    CALL channel_halt(substr, status)

    ! ###############################
    ! ### ARBITRARY ARRAY LENGTH ####
    ! ###############################

    ! Note: Here the length of the random state vector of the fortran system
    !       is arbitrarily used as array length ...
    CALL RANDOM_SEED(SIZE = n)

    ! ... dimension ...
    CALL new_dimension(status, DIMID_AL, 'n', n)
    CALL channel_halt(substr, status)

    ! #################
    ! DIAGNOSTIC OUTPUT
    ! #################
    CALL write_dimension(status)
    CALL channel_halt(substr, status)

  END SUBROUTINE main_channel_initialize_dims
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize_reprs

    USE messy_main_channel_repr,         ONLY: new_representation        &
                                             , write_representation_dc   &
                                             , set_representation_decomp &
                                             , write_representation      &
                                             , AUTO, IRANK               &
                                             , PIOTYPE_SGL, PIOTYPE_IND  &
                                             , PIOTYPE_COL

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'main_channel_initialize_reprs'
    INTEGER                        :: status
    INTEGER, DIMENSION(:), POINTER :: scdim => NULL()

    ! ##################
    ! ### GRIDPOINT ####
    ! ##################
    CALL new_representation(status, GP_3D_MID, 'GP_3D_MID' &
         , rank = 3, link = 'xxx-', dctype = 0                   &
         , dimension_ids = (/ DIMID_LON, DIMID_LEV, DIMID_LAT /) &
         , ldimlen       = (/ AUTO     , AUTO     , AUTO      /) &
         , output_order  = (/ 1,3,2 /)                           &
         , axis = 'XZY-'                                         &
         )
    CALL channel_halt(substr, status)

! mz_ab_20090921+
    ! ################################
    ! ### GRIDPOINT WITH BOUNDARY ####
    ! ################################
    CALL new_representation(status, GP_3D_MID_BND, 'GP_3D_MID_BND'   &
         , rank = 3, link = 'xxx-', dctype = 0                       &
         , dimension_ids = (/ DIMID_LON, DIMID_LEV, DIMID_LAT /)     &
         , ldimlen       = (/ AUTO     , AUTO     , AUTO      /)     &
         , nbounds       = (/ 2,0,2 /)                               &
         , output_order  = (/ 1,2,3 /)                               &
         , axis = 'XZY-'                                             &
         )
    CALL channel_halt(substr, status)
! mz_ab_20090921-

    ! ###############
    ! ### SCALAR ####
    ! ###############
    ALLOCATE(scdim(0))
    CALL new_representation(status, SCALAR, 'SCALAR'                 &
         , rank = 0, link = '----', dctype = 0                       &
         , dimension_ids = scdim &
         , ldimlen       = scdim &
         , axis = '----'                                             &
         )
    CALL channel_halt(substr, status)
    DEALLOCATE(scdim)
    NULLIFY(scdim)

    ! ##############
    ! ### ARRAY ####
    ! ##############
    CALL new_representation(status, ARRAY, 'ARRAY'                 &
         , rank = 1, link = 'x---', dctype = 0                     &
         , dimension_ids = (/ DIMID_AL /) &
         , ldimlen       = (/ AUTO /)     &
         , axis = 'N---'                                           &
         )
    CALL channel_halt(substr, status)

    ! #################
    ! DIAGNOSTIC OUTPUT
    ! #################
    CALL write_representation(status)
    CALL channel_halt(substr, status)

  END SUBROUTINE main_channel_initialize_reprs
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE initialize_restart_attributes(AMODE)

    USE messy_main_channel_io, ONLY: channel_init_restart
    USE messy_main_channel,    ONLY: new_attribute, get_attribute &
                                   , AF_RST_CMP, AF_RST_INP, AF_RST_NONE
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM

    IMPLICIT NONE
    INTRINSIC :: INT, MOD, SIGN, TRIM

    ! I/O
    INTEGER, INTENT(IN) :: AMODE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='initialize_restart_attributes'
    !
    INTEGER           :: iflag
    INTEGER           :: timestep
    INTEGER           :: status
    LOGICAL           :: lp

    SELECT CASE(AMODE)
       CASE(AMODE_WRITE)
          iflag    = AF_RST_NONE
       CASE(AMODE_READ)
          iflag    = AF_RST_CMP
       CASE(AMODE_INIT)
          iflag    = AF_RST_INP
    END SELECT

    ! - CURRENT TIME STEP
    CALL new_attribute(status, restart_att &
         , 'time', i=NH                    &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)

    ! CONTINUE ONLY DURING INITIALIZATION AFTER RESTART
    IF (AMODE /= AMODE_INIT) RETURN

    CALL channel_init_restart(status, lp, .TRUE. &
         , 'restart_'//cystr//'_'//basemodstr, restart_att)
    CALL channel_halt(substr, status)

    ! - CURRENT TIME STEP
    CALL get_attribute(status, restart_att, 'time' &
         , i=NH)
    CALL channel_halt(substr, status)

  END SUBROUTINE initialize_restart_attributes
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_init_timer

    USE messy_main_channel,      ONLY: get_channel_name
    USE messy_main_tools,        ONLY: match_wild, find_next_free_unit

    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_init_timer'
    INTEGER :: status
    INTEGER :: iou
    INTEGER :: i, js
    CHARACTER(LEN=STRLEN_CHANNEL) :: cname
    LOGICAL :: lexplicit

    ! SPACE FOR ACTIVE OUTPUT EVENTS (ONE PER CHANNEL)
    ALLOCATE(LOUTPUT_NOW(NCHANNEL))
    LOUTPUT_NOW(:) = .FALSE.
    ALLOCATE(IOUTPUT_FREQ(NCHANNEL))
    IOUTPUT_FREQ(:) = 0

    ALLOCATE(LTNF_NOW(NCHANNEL))
    LTNF_NOW(:) = .FALSE.
    ALLOCATE(INEWFILE_FREQ(NCHANNEL))
    INEWFILE_FREQ(:) = 0

    ! PRESET NON-EXPLICIT TO DEFAULT
    channel_loop: DO js=1, NCHANNEL

       CALL get_channel_name(status, js, cname)
       CALL channel_halt(substr, status)

       ! SET EXPLICIT OR DEFAULT
       lexplicit = .FALSE.
       event_loop: DO i=1, NMAXCHANNELS
          IF (TRIM(ADJUSTL(TIMER_CHANNEL(i)%cname)) == '') CYCLE
          IF (match_wild(TRIM(ADJUSTL(TIMER_CHANNEL(i)%cname)),&
               TRIM(cname))) THEN
             lexplicit  = .TRUE.
             EXIT
          END IF
       END DO event_loop

       IF (lexplicit) THEN
          IOUTPUT_FREQ(js) = TIMER_CHANNEL(i)%nstep
       ELSE
          IOUTPUT_FREQ(js) = TIMER_DEFAULT%nstep
       END IF

       ! SET EXPLICIT OR DEFAULT
       lexplicit = .FALSE.
       event_loop_tnf: DO i=1, NMAXCHANNELS
          IF (TRIM(ADJUSTL(TIMER_TNF_CHANNEL(i)%cname)) == '') CYCLE
          IF (match_wild(TRIM(ADJUSTL(TIMER_TNF_CHANNEL(i)%cname)),&
               TRIM(cname))) THEN
             lexplicit  = .TRUE.
             EXIT
          END IF
       END DO event_loop_tnf

       IF (lexplicit) THEN
          INEWFILE_FREQ(js) = TIMER_TNF_CHANNEL(i)%nstep
       ELSE
          INEWFILE_FREQ(js) = TIMER_TNF_DEFAULT%nstep
       END IF

    END DO channel_loop

  END SUBROUTINE main_channel_init_timer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_update_timer

    USE messy_main_channel,     ONLY: trigger_channel_output

    IMPLICIT NONE
    INTRINSIC :: MOD

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_update_timer'
    INTEGER :: status
    INTEGER :: js

    DO js = 1, NCHANNEL
       LOUTPUT_NOW(js) = (MOD(NH, IOUTPUT_FREQ(js)) == 0)
    END DO

    DO js = 1, NCHANNEL
       LTNF_NOW(js) = (MOD(NH, INEWFILE_FREQ(js)) == 0)
    END DO

    CALL trigger_channel_output(status, LOUTPUT_NOW, LTNF_NOW &
         , LFORCE_NEW_OUTPUT)
    CALL channel_halt(substr, status)
    LFORCE_NEW_OUTPUT = .FALSE.

  END SUBROUTINE main_channel_update_timer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_read_nml_cpl(status, iou)

    ! MODULE ROUTINE (SMIL)
    !
    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2004

    USE messy_main_tools,   ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_main_channel, ONLY: modstr

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CPL/ TIMER_DEFAULT, TIMER_CHANNEL &
         , TIMER_TNF_DEFAULT, TIMER_TNF_CHANNEL ! mz_pj_20080905

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_read_nml_cpl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_channel_read_nml_cpl
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_update_vars

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, TRIM, ADJUSTL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_update_vars'
    INTEGER :: status
    CHARACTER (8)   :: ydate
    CHARACTER (10)  :: ytime

    INTEGER :: iyr, imo, idy, ihr, imi, ise, ims

    CALL DATE_AND_TIME(ydate, ytime)

    READ(ydate,'(i4)') iyr
    READ(ydate,'(4x,i2)') imo
    READ(ydate,'(6x,i2)') idy

    READ(ytime,'(i2)') ihr
    READ(ytime,'(2x,i2)') imi
    READ(ytime,'(4x,i2)') ise
    READ(ytime,'(7x,i3)') ims

    yr = REAL(iyr,dp)
    mo = REAL(imo,dp)
    dy = REAL(idy,dp)
    hr = REAL(ihr,dp)
    mi = REAL(imi,dp)
    se = REAL(ise,dp)
    ms = REAL(ims,dp)

  END SUBROUTINE main_channel_update_vars
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_main_channel_bi
! **********************************************************************
