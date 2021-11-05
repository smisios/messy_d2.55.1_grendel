#include "messy_main_ppd_bi.inc"

! **********************************************************************
MODULE messy_s4d_si
! **********************************************************************

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY, PTR_3D_ARRAY
  USE messy_main_channel,       ONLY: STRLEN_OBJECT, STRLEN_CHANNEL
  USE messy_main_timer_event,   ONLY: time_event, io_time_event
  USE messy_main_constants_mem, ONLY: STRLEN_VLONG
  USE messy_s4d

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL
  
  TYPE T_TRACK_IO
!!$  CHARACTER(LEN=8)      :: name  = ''
     CHARACTER(LEN=STRLEN_VLONG) :: name  = ''
     CHARACTER(LEN=200)    :: fname = ''
     INTEGER               :: usw   = -1
     LOGICAL               :: lcol  = .FALSE.
     LOGICAL               :: lall  = .FALSE.
     REAL(DP)              :: rfill = 0.0_dp
     ! mz_sg_20180610+ switch for horizontal interpolation
     LOGICAL               :: linth = .TRUE.
     ! mz_sg_20180610-
     CHARACTER(LEN=STRLEN) :: str   = ''
  END TYPE T_TRACK_IO

  TYPE T_POSITION_DECOMP
     INTEGER,  DIMENSION(4) :: pe    = -1 ! PROCESS ID
     INTEGER,  DIMENSION(4) :: jrow  = -1 ! ROW INDEX
     INTEGER,  DIMENSION(4) :: jp    = -1 ! VECTOR INDEX
     REAL(DP), DIMENSION(4) :: w     = 0.0_dp ! weight for bilin. interpol.
  END TYPE T_POSITION_DECOMP

  TYPE T_TRACK
     !
     ! USER INTERFACE
     TYPE(T_TRACK_IO)                               :: io
     !
     ! op_pj_20190507+
     ! DIMENSION / INTERPOALTION METHOD
     INTEGER                :: d = 4  ! 1: nearest neighbour
     !                                ! 4: bilin. interpol
     ! op_pj_20190507-
     !
     ! POSITION DATA
     ! ... number of positions
     INTEGER                                        :: npos = 0
     ! ... list of positions
     TYPE(T_POSITION),        DIMENSION(:), POINTER :: pos => NULL()
     ! ... in decomposition
     TYPE(T_POSITION_DECOMP), DIMENSION(:), POINTER :: dec => NULL()
     !
     ! OUTPUT OF POSITION
     ! ... exact
     REAL(DP), POINTER :: lon => NULL()
     REAL(DP), POINTER :: lat => NULL()
     REAL(DP), POINTER :: pre => NULL()
     REAL(DP), POINTER :: ps  => NULL()
     !
     ! COUPLING FOR OUTPUT
     ! ... number of objects
     INTEGER                                      :: nobj = 0
     CHARACTER(LEN=STRLEN_CHANNEL), DIMENSION(:), POINTER :: cha  => NULL()
     CHARACTER(LEN=STRLEN_OBJECT),  DIMENSION(:), POINTER :: obj  => NULL()
     LOGICAL,                       DIMENSION(:), POINTER :: lex  => NULL()
     TYPE(PTR_3D_ARRAY),            DIMENSION(:), POINTER :: dat  => NULL()
     INTEGER,                       DIMENSION(:), POINTER :: rid  => NULL()
     TYPE(PTR_1D_ARRAY),            DIMENSION(:), POINTER :: col  => NULL()
     !
     ! WORKSPACE
     ! ... current position (index in list)
     INTEGER                                      :: cpos = 0
     ! ... on track ?
     LOGICAL                                      :: lt = .TRUE.
     ! ... level indices
     INTEGER,                       DIMENSION(:), POINTER :: nk => NULL()
     !
     INTEGER                                      :: domain_idx = 0
     !
  END TYPE T_TRACK

  INTEGER, PARAMETER :: NMAXTRACK = 100
  TYPE(T_TRACK_IO),     DIMENSION(NMAXTRACK), SAVE :: TRACK  ! CPL
  ! mz_ak_20060516+
  !TYPE(T_TRACK),        DIMENSION(NMAXTRACK), SAVE :: XTRACK 
  TYPE(T_TRACK), DIMENSION(:), POINTER, SAVE :: XTRACK => NULL()
  ! mz_ak_20060516-
  INTEGER,                                    SAVE :: NTRACK

  ! TIMER
  TYPE(io_time_event), SAVE :: TIMER_MONTHLY = &
       io_time_event(1, 'months','first',0)
  TYPE(time_event), DIMENSION(:), ALLOCATABLE, SAVE :: XTIMER_MONTHLY
  LOGICAL,          DIMENSION(:), ALLOCATABLE, SAVE :: LTRIG_MONTHLY
  
  TYPE(io_time_event), SAVE :: TIMER_DAILY = &
       io_time_event(1, 'days','first',0)
  TYPE(time_event), DIMENSION(:), ALLOCATABLE, SAVE :: XTIMER_DAILY 
  LOGICAL,          DIMENSION(:), ALLOCATABLE, SAVE :: LTRIG_DAILY
  
  ! NEW REPRESENTATIONS
  INTEGER, SAVE                      :: S4D_GP_1D_COLUMN_MID
  INTEGER, SAVE                      :: S4D_GP_1D_COLUMN_INT

  PUBLIC :: s4d_initialize
  PUBLIC :: s4d_init_memory
  PUBLIC :: s4d_init_coupling
  PUBLIC :: s4d_global_start
  PUBLIC :: s4d_write_output
  PUBLIC :: s4d_free_memory
  !PRIVATE :: s4d_read_nml_cpl
  !PRIVATE :: read_s4d_data_files

CONTAINS

! ======================================================================
! PUBLIC ROUTINES
! ======================================================================

! ----------------------------------------------------------------------
  SUBROUTINE s4d_initialize

    ! S4D MODULE ROUTINE (ECHAM5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2006

    USE messy_main_timer_bi,   ONLY: timer_event_init
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: str2chob, find_next_free_unit &
                                   , domains_from_string
    USE messy_main_channel_bi, ONLY: n_dom
    
    IMPLICIT NONE

    INTRINSIC :: TRIM, ADJUSTL, LEN_TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 's4d_initialize'
    INTEGER                        :: iou    ! I/O unit
    INTEGER                        :: status ! error status
    INTEGER                        :: i, j
    CHARACTER(LEN=32)              :: varname = '  '
    INTEGER, DIMENSION(:), POINTER :: domnum  => NULL()
    INTEGER                        :: nd, num
    INTEGER                        :: jg
    CHARACTER(LEN=*), PARAMETER    :: trigstr = 'next'

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL s4d_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('',substr)
    END IF

    CALL start_message_bi(modstr,'OUTPUT INITIALISATION ',substr)

    IF (p_parallel_io) THEN
       ! GET NUMBER OF TRACKS
       NTRACK = 1
       DO i=1, NMAXTRACK
          IF (TRIM(TRACK(i)%name) == '') CYCLE
          IF (TRIM(TRACK(i)%fname) == '') CYCLE
          IF (TRIM(TRACK(i)%str) == '') CYCLE
          CALL domains_from_string(status,TRACK(i)%name, n_dom, num)
          IF (status /= 0) CALL error_bi('error in namelist parsing', substr)
          NTRACK = NTRACK + num
       END DO
       NTRACK = NTRACK - 1
    END IF
    CALL p_bcast(NTRACK, p_io)
    ALLOCATE(XTRACK(NTRACK))
    
    IF (p_parallel_io) THEN
       ! GET NUMBER OF LOCATIONS
       NTRACK = 1
       ! COPY DATA AND PARSE STR
       DO i=1, NMAXTRACK
          IF (TRIM(TRACK(i)%name) == '') CYCLE
          IF (TRIM(TRACK(i)%fname) == '') CYCLE
          IF (TRIM(TRACK(i)%str) == '') CYCLE
          CALL domains_from_string(status,TRACK(i)%name, n_dom, num &
                ,varname, dnums=domnum)
          IF (status /= 0) CALL error_bi('error in namelist parsing', substr)
          IF (LEN_TRIM(ADJUSTL(varname)) > 9) &
               CALL error_bi('location name too long (max. 5(+4) characters)' &
               , substr)

          domain_loop: DO nd=1, SIZE(domnum)
             XTRACK(NTRACK)%io%name    = TRIM(ADJUSTL(varname))
             XTRACK(NTRACK)%domain_idx = domnum(nd)
             XTRACK(NTRACK)%io%fname   = TRIM(TRACK(i)%fname)
             XTRACK(NTRACK)%io%str     = TRIM(TRACK(i)%str)

             WRITE(*,*) 'TRACK        : ',TRIM(XTRACK(NTRACK)%io%name)
             WRITE(*,*) 'DOMAIN       : ',XTRACK(NTRACK)%domain_idx
             WRITE(*,*) '   FILE      : ',TRIM(XTRACK(NTRACK)%io%fname)&
                  &//'<DATE>.pos'

             XTRACK(NTRACK)%io%usw = TRACK(i)%usw
             SELECT CASE(XTRACK(NTRACK)%io%usw)
             CASE(-1)
                WRITE(*,*) '   SWITCHED OFF ... skipping'
                CYCLE
             CASE(0)
                WRITE(*,*) '   UPDATE    : daily'
             CASE(1)
                WRITE(*,*) '   UPDATE    : monthly'
             CASE DEFAULT
                WRITE(*,*) '   UNKNOWN UPDATE SWITCH ... skipping'
                CYCLE
             END SELECT

             XTRACK(NTRACK)%io%lcol  = TRACK(i)%lcol
             IF (XTRACK(NTRACK)%io%lcol) THEN
                WRITE(*,*) '   COLUMN    : YES'
             ELSE
                WRITE(*,*) '   COLUMN    : NO'
             END IF

             XTRACK(NTRACK)%io%lall  = TRACK(i)%lall
             XTRACK(NTRACK)%io%rfill  = TRACK(i)%rfill
             IF (XTRACK(NTRACK)%io%lall) THEN
                WRITE(*,*) '   FILL GAPS : YES (',XTRACK(NTRACK)%io%rfill,')'
             ELSE
                WRITE(*,*) '   FILL GAPS : NO'
             END IF

             ! mz_sg_20180610+: switch for horizontal interpolation
             XTRACK(NTRACK)%io%linth = TRACK(i)%linth
             IF (XTRACK(NTRACK)%io%linth) THEN
                WRITE(*,*) '   HORIZONTAL INTERPOLATION : YES'
             ELSE
                WRITE(*,*) '   HORIZONTAL INTERPOLATION : NO'
             END IF
             ! mz_sg_20180610-

             ! SET PRELIMINARY %nobj
             CALL str2chob(status, XTRACK(NTRACK)%io%str, XTRACK(NTRACK)%nobj &
                  , XTRACK(NTRACK)%cha, XTRACK(NTRACK)%obj)
             IF (status /= 0) THEN
                WRITE(*,*) '   ... ERROR IN STRING ... skipping'
                CYCLE
             END IF

             IF (XTRACK(NTRACK)%nobj == 0) THEN
                WRITE(*,*) '   ... EMPTY OUTPUT LIST ... skipping'
                CYCLE
             ELSE
                WRITE(*,*) '   REQUESTS  : ',XTRACK(NTRACK)%nobj
             END IF
             !
             WRITE(*,'(1x,a16,1x,a32)') 'CHANNEL','OBJECT(S)'
             WRITE(*,'(1x,a16,1x,a32)') '-------','---------'
             DO j=1, XTRACK(NTRACK)%nobj 
                WRITE(*,'(1x,a16,1x,a32)') TRIM(XTRACK(NTRACK)%cha(j)) &
                     , TRIM(XTRACK(NTRACK)%obj(j))
             END DO

             ! NEXT TRACK
             NTRACK = NTRACK + 1
             WRITE(*,*) '------------------------------------------------------'
          END DO domain_loop
          DEALLOCATE(domnum); NULLIFY(domnum)
       END DO
       NTRACK = NTRACK - 1
       IF (NTRACK > SIZE(XTRACK)) &
            CALL error_bi('error in namelist parsing', substr)
    END IF

    ! BROADCAST ALL RESULTS
    DO i=1, NTRACK
       ! I/O USER INTERFACE
       CALL p_bcast(XTRACK(i)%io%name, p_io)
       CALL p_bcast(XTRACK(i)%io%fname, p_io)
       CALL p_bcast(XTRACK(i)%io%usw, p_io)
       CALL p_bcast(XTRACK(i)%io%str, p_io)
       CALL p_bcast(XTRACK(i)%io%lcol, p_io)
       CALL p_bcast(XTRACK(i)%io%lall, p_io)
       CALL p_bcast(XTRACK(i)%io%rfill, p_io)
       CALL p_bcast(XTRACK(i)%io%linth, p_io) ! mz_sg_20180610
       !
       CALL p_bcast(XTRACK(i)%domain_idx, p_io)
       !
       ! CHANNEL/OBJECT NAMES
       CALL p_bcast(XTRACK(i)%nobj, p_io)
       IF (.NOT. p_parallel_io) THEN
          ALLOCATE(XTRACK(i)%cha(XTRACK(i)%nobj))
          ALLOCATE(XTRACK(i)%obj(XTRACK(i)%nobj))
       END IF
       DO j=1, XTRACK(i)%nobj
          CALL p_bcast(XTRACK(i)%cha(j), p_io)
          CALL p_bcast(XTRACK(i)%obj(j), p_io)
       END DO
       !
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NTRACK,' TRACK(S) INITIALIZED !'
    END IF

    ! INITALIZE EVENTS
    ALLOCATE(XTIMER_DAILY(n_dom))
    ALLOCATE(XTIMER_MONTHLY(n_dom))
    DO jg=1, n_dom
       CALL timer_event_init(XTIMER_DAILY(jg),   TIMER_DAILY &
            , 's4d_daily',   trigstr, jg)
       CALL timer_event_init(XTIMER_MONTHLY(jg), TIMER_MONTHLY &
            , 's4d_monthly', trigstr, jg)
    END DO
    ALLOCATE(LTRIG_DAILY(n_dom))
    ALLOCATE(LTRIG_MONTHLY(n_dom))
    
    CALL end_message_bi(modstr,'OUTPUT INITIALISATION ',substr)

  END SUBROUTINE s4d_initialize
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE s4d_init_memory

    ! S4D MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2006

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: DIMID_LEV, DIMID_ILEV, DC_BC
    USE messy_main_grid_def_mem_bi,  ONLY: nlev
    USE messy_main_channel_repr,     ONLY: new_representation, AUTO &
                                         , set_representation_decomp &
                                         , IRANK, PIOTYPE_COL

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 's4d_init_memory'
    INTEGER                      :: status
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    CALL start_message_bi(modstr, 'NEW REPRESENTATIONS', substr)

    nseg = 1
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:) = 1
    meml(:,:) = 1
    memu(:,:) = 1

    ! -------------------

    ! NEW REPRESENTATIONS
    CALL new_representation(status, S4D_GP_1D_COLUMN_MID, 'S4D_GP_1D_COL_MID' &
         , rank = 1, link = 'x---', dctype = DC_BC                   &
         , dimension_ids = (/ DIMID_LEV /) &
         , ldimlen       = (/ AUTO  /)     &
         , axis = 'Z---'                   &
         )
    CALL channel_halt(substr, status)

    start(:,1) = 1
    cnt(:,1)   = nlev
    meml(:,1)  = 1
    memu(:,1)  = nlev
    
    CALL set_representation_decomp(status, S4D_GP_1D_COLUMN_MID &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -------------------

    CALL new_representation(status, S4D_GP_1D_COLUMN_INT, 'S4D_GP_1D_COL_INT' &
         , rank = 1, link = 'x---', dctype = DC_BC                       &
         , dimension_ids = (/ DIMID_ILEV /) &
         , ldimlen       = (/ AUTO  /)     &
         , axis = 'Z---'                   &
         )
    CALL channel_halt(substr, status)

    start(:,1) = 1
    cnt(:,1)   = nlev+1
    meml(:,1)  = 1
    memu(:,1)  = nlev+1
    
    CALL set_representation_decomp(status, S4D_GP_1D_COLUMN_INT &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -------------------

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    CALL end_message_bi(modstr, 'NEW REPRESENTATIONS', substr)

  END SUBROUTINE s4d_init_memory
! ----------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE s4d_init_coupling

    ! S4D MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! define specific channel(s) and allocate memory for
    ! global fields
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2006

    USE messy_main_blather_bi,       ONLY: error_bi, warning_bi, info_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: SCALAR, GP_3D_MID, GP_3D_INT  &
                                         , GP_3D_1LEV, GP_2D_HORIZONTAL
    USE messy_main_grid_def_mem_bi,  ONLY: nlev
    USE messy_main_channel,          ONLY: new_channel, new_channel_object   &
                                         , new_attribute, get_channel_object &
                                         , get_attribute, get_channel_info &
                                         , get_channel_object_info
    USE messy_main_constants_mem,    ONLY: STRLEN_ULONG
    USE messy_main_tools,            ONLY: match_wild
    USE messy_main_channel_mem,      ONLY: dom_curid
    USE messy_main_channel_repr,     ONLY: get_representation_id

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 's4d_init_coupling'
    INTEGER                      :: status
    INTEGER                      :: reprid
    INTEGER                      :: reprid_new
    INTEGER                      :: i, jr, n, nptr, jo
    CHARACTER(LEN=STRLEN_ULONG)  :: charatt
    CHARACTER(LEN=STRLEN_OBJECT), DIMENSION(:), POINTER :: ONAMES => NULL()

    CHARACTER(LEN=*), DIMENSION(4), PARAMETER :: olist = &
         (/ 'tlon  ', 'tlat  ', 'tpress', 'tps   ' /)

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    track_loop: DO i=1, NTRACK

       IF (XTRACK(i)%domain_idx /= dom_curid) CYCLE

       ! OPEN ONE CHANNEL PER TRACK
       CALL new_channel(status, modstr//'_'//TRIM(XTRACK(i)%io%name))
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//TRIM(XTRACK(i)%io%name) &
            , 's4d_track', c=TRIM(XTRACK(i)%io%fname)//'<DATE>.pos' )
       CALL channel_halt(substr, status)

       ! ADD OBJECTS FOR EXACT POSITION
       ! ... longitude
       CALL new_channel_object(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name) &
            , 'tlon' &
            , p0=XTRACK(i)%lon &
            , reprid=SCALAR )
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name) &
            , 'tlon'  &
            , 'long_name', c='track longitude')
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name) &
            , 'tlon'  &
            , 'units', c='degrees_east')
       CALL channel_halt(substr, status)

       ! ... latitude
       CALL new_channel_object(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name) &
            , 'tlat' &
            , p0=XTRACK(i)%lat &
            , reprid=SCALAR )
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name) &
            , 'tlat'  &
            , 'long_name', c='track latitude')
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name) &
            , 'tlat'  &
            , 'units', c='degrees_north')
       CALL channel_halt(substr, status)

       ! ... pressure
       CALL new_channel_object(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name) &
            , 'tpress' &
            , p0=XTRACK(i)%pre &
            , reprid=SCALAR )
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name) &
            , 'tpress'  &
            , 'long_name', c='track pressure')
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name) &
            , 'tpress'  &
            , 'units', c='hPa')
       CALL channel_halt(substr, status)

       ! ... surface pressure
       CALL new_channel_object(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name) &
            , 'tps' &
            , p0=XTRACK(i)%ps &
            , reprid=SCALAR )
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name) &
            , 'tps'  &
            , 'long_name', c='track surface pressure')
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name) &
            , 'tps'  &
            , 'units', c='Pa')
       CALL channel_halt(substr, status)

       DO jo=1, SIZE(olist)
          CALL new_attribute(status, modstr//'_'//TRIM(XTRACK(i)%io%name) &
               , TRIM(olist(jo))                                          &
               , 'missing_value', r=XTRACK(i)%io%rfill                    &
               , loverwrite=.TRUE.)
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, modstr//'_'//TRIM(XTRACK(i)%io%name) &
               , TRIM(olist(jo))                                          &
               , '_FillValue', r=XTRACK(i)%io%rfill                       &
               , loverwrite=.TRUE.)
          CALL channel_halt(substr, status)
       END DO

       ! CALCULATE REQUIRED NUMBER OF POINTERS
       nptr = 0
       DO jr=1,  XTRACK(i)%nobj  ! LOOP OVER REQUESTS
          CALL get_channel_info(status, TRIM(XTRACK(i)%cha(jr)) &
               , ONAMES = ONAMES)
          IF (status /= 3003) THEN ! CHANNEL (NAME) DOES NOT EXIST
             CALL channel_halt(substr, status)
          ELSE
             CALL warning_bi( &
                  ' ... channel '''//TRIM(XTRACK(i)%cha(jr))//''' does not exist ... skipping', substr)
             CYCLE
          END IF
          DO jo = 1, SIZE(ONAMES)
             IF ( match_wild(TRIM(XTRACK(i)%obj(jr)), TRIM(ONAMES(jo))) ) THEN
                nptr = nptr + 1
             END IF
          END DO
       END DO
       !
       IF (ASSOCIATED(ONAMES)) THEN
          DEALLOCATE(ONAMES)
          NULLIFY(ONAMES)
       END IF

       ! op_pj_20190507+
       ! SET DIMENSION / INTERPOLATION METHOD
       IF (XTRACK(i)%io%linth) THEN
          XTRACK(i)%d = 4 ! bilinear interpolation
       ELSE
          XTRACK(i)%d = 1 ! nearest neighbour
       END IF
       ! op_pj_20190507-

       ! ALLOCATE SPACE FOR POINTER ARRAYS (DATA, REPR-ID, COLUMN OBJECT)
       ALLOCATE(XTRACK(i)%dat(nptr))   ! POINTER TO CHANNEL OBJECT
       ALLOCATE(XTRACK(i)%lex(nptr))   ! channel / object exsists
       XTRACK(i)%lex(:) = .FALSE.
       ALLOCATE(XTRACK(i)%rid(nptr))   ! REPRESENTATION OF CHANNEL OBJECT
       ALLOCATE(XTRACK(i)%col(nptr))   ! POINTER TO COLUMN CHANNEL OBJECT
       ALLOCATE(XTRACK(i)%nk(nptr))    ! NUMBER OF DATA LEVELS
       XTRACK(i)%nk(:) = 0

       ! RESET NUMBER OF REQUESTS
       n = XTRACK(i)%nobj    ! SAVE NUMBER OF REQUESTS FOR LOOP BELOW       
       XTRACK(i)%nobj = nptr ! NUMBER OF MATCHING OBJECTS
       ! (SO FAR XTRACK(i)%nobj CONTAINED NUMBER OF REQUESTS)

       ! ADD REQUESTED CHANNEL OBJECTS
       nptr = 0
       request_loop: DO jr=1, n  ! LOOP OVER REQUESTS

          ! GET ALL POTENTIAL OBJECTS
          CALL get_channel_info(status, TRIM(XTRACK(i)%cha(jr)) &
               , ONAMES = ONAMES)
          IF (status == 3003) THEN ! CHANNEL (NAME) DOES NOT EXIST
             CYCLE
          ELSE
             CALL channel_halt(substr, status)
          END IF

          object_loop: DO jo = 1, SIZE(ONAMES)

             ! CHECK IF POTENTIAL OBJECT IS MATCHING
             IF ( match_wild(TRIM(XTRACK(i)%obj(jr)), TRIM(ONAMES(jo))) ) THEN
                nptr = nptr + 1      ! NEXT POINTER INDEX
             ELSE
                CYCLE
             END IF
             
             CALL info_bi( &
                  '          '//TRIM(XTRACK(i)%cha(jr))//' - '//&
                  &TRIM(ONAMES(jo)), substr)

             ! INIT
             NULLIFY(XTRACK(i)%col(nptr)%ptr)
             NULLIFY(XTRACK(i)%dat(nptr)%ptr)

             CALL get_channel_object(status &
                  , TRIM(XTRACK(i)%cha(jr)), TRIM(ONAMES(jo)) &
                  , p3=XTRACK(i)%dat(nptr)%ptr )
             IF (status /= 0) THEN
                CALL warning_bi( &
                     ' ... channel object '''//TRIM(ONAMES(jo)) &
                     &//''' does not exist ... skipping', substr)
                NULLIFY(XTRACK(i)%col(nptr)%ptr)
                NULLIFY(XTRACK(i)%dat(nptr)%ptr)
                CYCLE
             END IF

             CALL get_channel_object_info(status &
                  , TRIM(XTRACK(i)%cha(jr)), TRIM(ONAMES(jo)) &
                  , reprid=reprid)
             CALL channel_halt(substr, status)
             !
             XTRACK(i)%rid(nptr) = reprid
             !
             IF (reprid == GP_3D_MID) THEN
                XTRACK(i)%nk(nptr) = nlev
                IF (XTRACK(i)%io%lcol) THEN
                   CALL get_representation_id(status, 'S4D_GP_1D_COL_MID' &
                        , reprid=reprid_new)
                   CALL channel_halt(substr, status)
                ELSE
                   reprid_new = SCALAR
                END IF
             ELSE IF (reprid == GP_3D_INT) THEN
                XTRACK(i)%nk(nptr) = nlev+1
                IF (XTRACK(i)%io%lcol) THEN
                   CALL get_representation_id(status, 'S4D_GP_1D_COL_INT' &
                        , reprid=reprid_new)
                   CALL channel_halt(substr, status)
                ELSE
                   reprid_new = SCALAR
                END IF
             ELSE IF (reprid == GP_3D_1LEV) THEN
                XTRACK(i)%nk(nptr) = 1
                reprid_new = SCALAR
             ELSE IF (reprid == GP_2D_HORIZONTAL) THEN
                XTRACK(i)%nk(nptr) = 1
                reprid_new = SCALAR
             ELSE
                CALL warning_bi( &
                     '            ... representation not supported ... skipping', substr)
                NULLIFY(XTRACK(i)%col(nptr)%ptr)
                NULLIFY(XTRACK(i)%dat(nptr)%ptr)
                CYCLE
             END IF

             CALL new_channel_object(status &
                  , modstr//'_'//TRIM(XTRACK(i)%io%name) &
                  , TRIM((XTRACK(i)%cha(jr)))//'_'//TRIM(ONAMES(jo))  &
                  , p1=XTRACK(i)%col(nptr)%ptr &
                  , reprid=reprid_new )
             CALL channel_halt(substr, status)

             ! COPY ATTRIBUTES
             CALL get_attribute(status &
                  , TRIM(XTRACK(i)%cha(jr)), TRIM(ONAMES(jo)) &
                  , 'long_name', c=charatt)
             IF (status == 0) THEN
                CALL new_attribute(status &
                     , modstr//'_'//TRIM(XTRACK(i)%io%name) &
                     , TRIM((XTRACK(i)%cha(jr)))//'_'//TRIM(ONAMES(jo))  &
                     , 'long_name', c=TRIM(charatt))
                CALL channel_halt(substr, status)
             END IF
             CALL get_attribute(status &
                  , TRIM(XTRACK(i)%cha(jr)), TRIM(ONAMES(jo)) &
                  , 'units', c=charatt)
             IF (status == 0) THEN
                CALL new_attribute(status &
                     , modstr//'_'//TRIM(XTRACK(i)%io%name) &
                     , TRIM((XTRACK(i)%cha(jr)))//'_'//TRIM(ONAMES(jo))  &
                     , 'units', c=TRIM(charatt))
                CALL channel_halt(substr, status)
             END IF

             ! op_pj_20180509+
             CALL new_attribute(status, modstr//'_'//TRIM(XTRACK(i)%io%name) &
                  , TRIM((XTRACK(i)%cha(jr)))//'_'//TRIM(ONAMES(jo))         &
                  , 'missing_value', r=XTRACK(i)%io%rfill                    &
                  , loverwrite=.TRUE.)
             CALL channel_halt(substr, status)
             !
             CALL new_attribute(status, modstr//'_'//TRIM(XTRACK(i)%io%name) &
                  , TRIM((XTRACK(i)%cha(jr)))//'_'//TRIM(ONAMES(jo))         &
                  , '_FillValue', r=XTRACK(i)%io%rfill                       &
                  , loverwrite=.TRUE.)
             CALL channel_halt(substr, status)
             ! op_pj_20180509-

             XTRACK(i)%lex(nptr) = .TRUE.

          END DO object_loop

       END DO request_loop

       IF (nptr /= XTRACK(i)%nobj) THEN
          CALL error_bi('something went wrong with pointer counting', substr)
       END IF

       IF (ASSOCIATED(ONAMES)) DEALLOCATE(ONAMES)

    END DO track_loop

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE s4d_init_coupling
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE s4d_global_start

    USE messy_main_timer_bi,      ONLY: event_state
    USE messy_main_timer,         ONLY: YEAR, MONTH, DAY, lresume, lstart &
                                      , HOUR, MINUTE, SECOND
    USE messy_main_blather_bi,    ONLY: info_bi
    USE messy_main_channel_mem,   ONLY: dom_curid
    USE messy_main_timer,         ONLY: trigger_date => next_date

    IMPLICIT NONE
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 's4d_global_start'
    CHARACTER(LEN=25)           :: dtstr

    IF (lstart .OR. lresume) THEN
       ! READ DATA FOR CURRENT DATE
       CALL read_s4d_data_files(YEAR, MONTH, DAY, .TRUE., .TRUE.)
    END IF

    ! CHECK EVENTS: daily, monthly
    LTRIG_DAILY(dom_curid)   = event_state(XTIMER_DAILY(dom_curid) &
         , trigger_date)
    LTRIG_MONTHLY(dom_curid) = event_state(XTIMER_MONTHLY(dom_curid) &
         , trigger_date)

    IF (LTRIG_DAILY(dom_curid)) THEN
       write(dtstr,'(i2,1a,1x,i4.4,2(a1,i2.2),1x,3(a1,i2.2),a1)') &
            dom_curid, ':' &
            , YEAR, '-', MONTH, '-', DAY &
            , '(', HOUR, ':', MINUTE, ':', SECOND,')'
       CALL info_bi('TRIGGER DAILY   UPDATE, DOMAIN '//dtstr//' ', substr)
    ENDIF
    IF (LTRIG_MONTHLY(dom_curid)) THEN
       write(dtstr,'(i2,1a,1x,i4.4,2(a1,i2.2),1x,3(a1,i2.2),a1)') &
            dom_curid, ':' &
            , YEAR, '-', MONTH, '-', DAY &
            , '(', HOUR, ':', MINUTE, ':', SECOND,')'
       CALL info_bi('TRIGGER MONTHLY UPDATE, DOMAIN '//dtstr//' ', substr)
    ENDIF

  END SUBROUTINE s4d_global_start
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE s4d_write_output

    USE messy_main_mpi_bi,            ONLY: p_pe, p_bcast
    USE messy_main_timer,             ONLY: YEAR,MONTH,DAY,HOUR,MINUTE,SECOND &
                                          , YEAR_NEXT,MONTH_NEXT,DAY_NEXT     &
                                          , delta_time
    USE messy_main_data_bi,           ONLY: press_3d, pressi_3d, aps
    USE messy_main_channel_bi,        ONLY: GP_3D_MID, GP_3D_INT, GP_3D_1LEV &
                                          , GP_2D_HORIZONTAL
    USE messy_main_channel_error_bi,  ONLY: channel_halt
    USE messy_main_timer,             ONLY: julian_day
    USE messy_main_tools,             ONLY: nn_index
    USE messy_main_channel,           ONLY: set_channel_output
    USE messy_main_channel_mem,       ONLY: dom_curid

    IMPLICIT NONE

    INTRINSIC :: ABS, REAL, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 's4d_write_output'
    INTEGER                     :: i, j       ! loop counter
    REAL(dp)                    :: sjdf, tjdf ! julian day + fraction
    REAL(dp)                    :: dt         ! model time step (in days)
    REAL(dp)                    :: dtmin      ! time distance (in days)
    INTEGER                     :: cpos, npos, ccpos ! position counters
    LOGICAL                     :: lcalc      ! calculate this time step ?
    INTEGER                     :: jn         ! neighbour counter (1...4)
    REAL(dp)                    :: zval       ! scalar variable on one PE
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: zcol, pres ! column and pressure
    !                                                   ! on one PE
    REAL(dp), DIMENSION(:),   ALLOCATABLE :: zsum       ! result (interpolated)
    INTEGER                     :: jp, jrow, p_pe_track ! indices on this PE
    INTEGER                     :: nk       ! number of levels
    INTEGER                     :: k1, k2   ! level indices
    REAL(DP)                    :: wp       ! weight for pressure interpol.
    LOGICAL                     :: lout     ! output now ?
    INTEGER                     :: status
    REAL(DP)                    :: zpre    ! pressure in Pa

    ! INIT
    ! JULIAN DAY + FRACTION OF SIMULATION TIME
    sjdf = julian_day(REAL(DAY, DP), MONTH, YEAR) &
         + REAL(HOUR,dp)/24.0_dp        &
         + REAL(MINUTE,dp)/1440.0_dp    &
         + REAL(SECOND,dp)/86400.0_dp
    dt = delta_time/86400.0_dp
    
    track_loop: DO i=1, NTRACK

       IF (XTRACK(i)%domain_idx /= dom_curid) CYCLE

       ! STILL ON TRACK ?
       IF (.NOT. XTRACK(i)%lt) THEN
          CALL set_channel_output(status &
               , modstr//'_'//TRIM(XTRACK(i)%io%name), .FALSE.)
          CALL channel_halt(substr, status)
          CYCLE
       END IF

       ! PERFORM CALCULATION IN THIS TIME STEP
       lcalc = .TRUE.

       ! SEARCH TRACK TIME STEP IN THE MIDDLE OF SIMULATION TIME STEP
       npos = XTRACK(i)%npos
       cpos = XTRACK(i)%cpos  ! = 0 at start and after restart
       time_position: DO

          ! INCREMENT POSITION COUNTER
          cpos = cpos + 1

          ! END OF TRACK REACHED
          IF (cpos > npos) THEN
             XTRACK(i)%cpos = npos + 1   ! STORE COUNTER
             XTRACK(i)%lt = .FALSE.  ! skip earlier in next time step
             lcalc = .FALSE.         ! no calculation below
             ! op_pj_2017032: The next line is just to avoid an "out of array
             !                bounds" error (with run-time checks) below:
             !        ( lcalc = lcalc .AND. (XTRACK(i)%dec(cpos)%pe(jn) >= 0) )
             cpos = cpos - 1 ! op_pj_20170322
             EXIT
          END IF

          ! CALCULATE TIME
          ! JULIAN DAY + FRACTION OF TRACK POSITION
          tjdf = julian_day( &
               REAL(XTRACK(i)%pos(cpos)%dy, dp) &
               , XTRACK(i)%pos(cpos)%mo         &
               , XTRACK(i)%pos(cpos)%yr )                        &
               + REAL(XTRACK(i)%pos(cpos)%ho,dp)/24.0_dp        &
               + REAL(XTRACK(i)%pos(cpos)%mi,dp)/1440.0_dp      &
               + REAL(XTRACK(i)%pos(cpos)%se,dp)/86400.0_dp

          ! FUTURE TRACK POSITION
          IF (tjdf >= sjdf+dt) THEN
             ! do not change counter: XTRACK(i)%cpos
             lcalc = .FALSE.
             EXIT ! exit time_position loop
          END IF

          ! PAST TRACK POSITION: TRY NEXT POSITION
          IF (tjdf < sjdf) THEN
             XTRACK(i)%cpos = cpos ! update counter
             CYCLE ! next time_position
          END IF

          ! NOW: sjdf <= tjdf <= sjdf+dt, i.e.,
          ! BETWEEN CURRENT AND NEXT SIMULATION TIME STEP
          lcalc = .TRUE.

          ! CURRENT TRACK POSITION: SEARCH APPROX. MID POSITION
          !                         IN CURRENT TIME INTERVAL
          dtmin = ABS(tjdf - (sjdf+dt))
          ccpos = cpos
          DO
             ccpos = ccpos + 1
             IF (ccpos > npos) EXIT ! END OF TRACK REACHED
             !
             tjdf = julian_day( &
                  REAL(XTRACK(i)%pos(ccpos)%dy, dp) &
                  , XTRACK(i)%pos(ccpos)%mo         &
                  , XTRACK(i)%pos(ccpos)%yr )                        &
                  + REAL(XTRACK(i)%pos(ccpos)%ho,dp)/24.0_dp         &
                  + REAL(XTRACK(i)%pos(ccpos)%mi,dp)/1440.0_dp       &
                  + REAL(XTRACK(i)%pos(ccpos)%se,dp)/86400.0_dp
             IF (ABS(tjdf - (sjdf+dt)) < dtmin) THEN
                dtmin = ABS(tjdf - (sjdf+dt))
                cpos = ccpos
             ELSE
                EXIT
             END IF
          END DO

          ! update counter
          XTRACK(i)%cpos = cpos

          ! exit time position loop
          EXIT

       END DO time_position

       ! op_pj_20110524+
       ! just in time, but probably out of region ...
       ! if this should work correctly, locate_in_decomp(_4) must
       ! return negative PE, if point is outside region
       DO jn=1, XTRACK(i)%d
          lcalc = lcalc .AND. (XTRACK(i)%dec(cpos)%pe(jn) >= 0)
       END DO
       ! op_pj_20110524-

       ! TRIGGER/SUPPRESS CHANNEL OUTPUT
       lout = lcalc .OR. &
            (XTRACK(i)%io%lall .AND. &
            ((XTRACK(i)%cpos >0) .AND. (XTRACK(i)%cpos <= XTRACK(i)%npos)) )
       CALL set_channel_output(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name), lout)
       CALL channel_halt(substr, status)

       IF (.NOT. lcalc) THEN
          IF (XTRACK(i)%io%lall) THEN
             ! FILL VALUE
             XTRACK(i)%lat = XTRACK(i)%io%rfill
             XTRACK(i)%lon = XTRACK(i)%io%rfill
             XTRACK(i)%pre = XTRACK(i)%io%rfill
             XTRACK(i)%ps  = XTRACK(i)%io%rfill
             !
             obj_loop: DO j=1, XTRACK(i)%nobj
                IF (.NOT. XTRACK(i)%lex(j)) CYCLE
                XTRACK(i)%col(j)%ptr(:) = XTRACK(i)%io%rfill
             END DO obj_loop
          END IF
          CYCLE  ! NEXT TRACK
       END IF

       ! FROM NOW ON ONLY, IF TRACK TIME STEP IS IN CURRENT TIME INTERVAL

       ! COPY EXCACT POSITION TO CHANNEL OBJECT
       XTRACK(i)%lat = XTRACK(i)%pos(cpos)%lat
       XTRACK(i)%lon = XTRACK(i)%pos(cpos)%lon
       XTRACK(i)%pre = XTRACK(i)%pos(cpos)%pre
       !
       ! CALCULATE SURFACE PRESSURE
       XTRACK(i)%ps = 0.0_dp ! INIT
       DO jn=1, XTRACK(i)%d
          ! TRACK NEIGHBOUR IS ON THIS PE
          p_pe_track = XTRACK(i)%dec(cpos)%pe(jn)
          jp         = XTRACK(i)%dec(cpos)%jp(jn)
          jrow       = XTRACK(i)%dec(cpos)%jrow(jn)
          zval       = 0.0_dp
          ! GET VALUE
          IF (p_pe == p_pe_track) zval = aps(jp,jrow)
          ! BROADCAST TO ALL OTHER PEs
          CALL p_bcast(zval, p_pe_track)
          ! WEIGHTED SUM ON ALL PEs
          XTRACK(i)%ps = XTRACK(i)%ps + zval * XTRACK(i)%dec(cpos)%w(jn)
       END DO

       object_loop: DO j=1, XTRACK(i)%nobj

          IF (.NOT. XTRACK(i)%lex(j)) CYCLE

          ! NUMBER OF LEVELS IN DATA COLUMN
          nk = XTRACK(i)%nk(j)
          ALLOCATE(zcol(4,nk))  ! DATA IN COLUMNS
          zcol(:,:)  = 0.0_dp ! INIT
          ALLOCATE(pres(4,nk))  ! PRESSURE / SURFACE PRESSURE COLUMNS
          pres(:,:)  = 0.0_dp ! INIT

          ! COMPLETE COLUMN
          DO jn=1, XTRACK(i)%d
             ! TRACK NEIGHBOUR IS ON THIS PE
             p_pe_track = XTRACK(i)%dec(cpos)%pe(jn)
             jp         = XTRACK(i)%dec(cpos)%jp(jn)
             jrow       = XTRACK(i)%dec(cpos)%jrow(jn)
             ! GET VALUE ...
             IF (p_pe == p_pe_track) THEN
                IF (XTRACK(i)%rid(j) == GP_2D_HORIZONTAL) THEN
                   zcol(jn,:) = XTRACK(i)%dat(j)%ptr(jp,jrow,1)
                ELSE
                   zcol(jn,:) = XTRACK(i)%dat(j)%ptr(_RI_XYZ__(jp,jrow,:))
                END IF
                ! ... AND PRESSURE
                IF (XTRACK(i)%rid(j) == GP_3D_MID) THEN
                   pres(jn,:) = press_3d(_RI_XYZ__(jp,jrow,:))
                ELSE IF (XTRACK(i)%rid(j) == GP_3D_INT) THEN
                   pres(jn,:) = pressi_3d(_RI_XYZ__(jp,jrow,:))
                ELSE IF (XTRACK(i)%rid(j) == GP_3D_1LEV) THEN
                   pres(jn,:) = aps(jp,jrow)
                ELSE IF (XTRACK(i)%rid(j) == GP_2D_HORIZONTAL) THEN
                   pres(jn,:) = aps(jp,jrow)
                ELSE
                   ! ERROR
                END IF
             END IF
             ! BROADCAST TO ALL OTHER PEs
             CALL p_bcast(zcol(jn,:), p_pe_track)
             CALL p_bcast(pres(jn,:), p_pe_track)
          END DO

          IF (XTRACK(i)%io%lcol .OR. (nk==1)) THEN
             ! WEIGHTED SUM OF 1 OR 4 COLUMNS ON ALL PEs
             ALLOCATE(zsum(nk))       ! COLUMN
             zsum(:) = 0.0_dp
             DO jn=1, XTRACK(i)%d
                zsum(:) = zsum(:) &
                     + zcol(jn,:) * XTRACK(i)%dec(cpos)%w(jn)
             END DO
          ELSE
             ! LEVEL INTERPOLATION
             ALLOCATE(zsum(1))
             zsum(:) = 0.0_dp
             DO jn=1, XTRACK(i)%d
                zpre = XTRACK(i)%pos(cpos)%pre * 100.0_DP ! Pa
                !
                IF (zpre < pres(jn,  1)) THEN
                   zpre = pres(jn,  1)
                   k1 = 1
                   k2 = 1
                   wp = 0.5_dp
                ELSE IF (zpre > pres(jn,  nk)) THEN
                   zpre = pres(jn,  nk)
                   k1 = nk
                   k2 = nk
                   wp = 0.5_dp
                ELSE
                   CALL nn_index(pres(jn,:), zpre, k1, k2)
                   wp = ABS( (zpre - pres(jn,k1)) / &
                        ( pres(jn, k2) - pres(jn, k1) ) )
                END IF
                !
                zsum(:) = zsum(:) &
                     + ( (1._dp - wp) * zcol(jn,k1)   &
                     +   wp           * zcol(jn,k2) ) &
                     * XTRACK(i)%dec(cpos)%w(jn)
             END DO
          END IF

          ! TRANSFER COLUMN / INTERPOLATED VALUE
          XTRACK(i)%col(j)%ptr(:) = zsum(:)

          ! CLEAN UP
          DEALLOCATE(zcol)
          DEALLOCATE(pres)
          DEALLOCATE(zsum)

       END DO object_loop

    END DO track_loop

    ! NOW CHECK IF FILES MUST BE UPDATED
    IF (.NOT. (LTRIG_DAILY(dom_curid) .OR. LTRIG_MONTHLY(dom_curid))) RETURN

    ! UPDATE DATA
    ! EVENT TRIGGERED, IF current_date IS 'last step of day'
    ! or 'last step of month', respectively
    ! -> USE NEXT DATE FOR INPUT FILE NAMES
    CALL read_s4d_data_files(YEAR_NEXT, MONTH_NEXT, DAY_NEXT &
         , LTRIG_DAILY(dom_curid), LTRIG_MONTHLY(dom_curid))

  END SUBROUTINE s4d_write_output
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE s4d_free_memory

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! LOCAL
    INTEGER :: i

    track_loop: DO i=1, NTRACK
       !
       IF (ASSOCIATED(XTRACK(i)%pos)) DEALLOCATE(XTRACK(i)%pos)
       NULLIFY(XTRACK(i)%pos)
       IF (ASSOCIATED(XTRACK(i)%dec)) DEALLOCATE(XTRACK(i)%dec)
       NULLIFY(XTRACK(i)%dec)
       !
       IF (ASSOCIATED(XTRACK(i)%cha)) DEALLOCATE(XTRACK(i)%cha)
       NULLIFY(XTRACK(i)%cha)
       IF (ASSOCIATED(XTRACK(i)%obj)) DEALLOCATE(XTRACK(i)%obj)
       NULLIFY(XTRACK(i)%obj)
       !
       IF (ASSOCIATED(XTRACK(i)%dat)) DEALLOCATE(XTRACK(i)%dat)
       NULLIFY(XTRACK(i)%dat)
       IF (ASSOCIATED(XTRACK(i)%rid)) DEALLOCATE(XTRACK(i)%rid)
       NULLIFY(XTRACK(i)%rid)
       IF (ASSOCIATED(XTRACK(i)%col)) DEALLOCATE(XTRACK(i)%col)
       NULLIFY(XTRACK(i)%col)
       IF (ASSOCIATED(XTRACK(i)%lex)) DEALLOCATE(XTRACK(i)%lex)
       NULLIFY(XTRACK(i)%lex)
       IF (ASSOCIATED(XTRACK(i)%nk)) DEALLOCATE(XTRACK(i)%nk)
       NULLIFY(XTRACK(i)%nk)
       !
    END DO track_loop

    IF (ASSOCIATED(XTRACK)) DEALLOCATE(XTRACK)
    NULLIFY(XTRACK)

    DEALLOCATE(LTRIG_DAILY)
    DEALLOCATE(LTRIG_MONTHLY)
    DEALLOCATE(XTIMER_DAILY)
    DEALLOCATE(XTIMER_MONTHLY)
    
  END SUBROUTINE s4d_free_memory
! ---------------------------------------------------------------------

! ======================================================================
! PRIVATE ROUTINES
! ======================================================================

! ----------------------------------------------------------------------
  SUBROUTINE s4d_read_nml_cpl(status, iou)

    ! S4D MODULE ROUTINE (ECHAM5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2006

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 's4d_read_nml_cpl'

    NAMELIST /CPL/ TRACK

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

  END SUBROUTINE s4d_read_nml_cpl
! ----------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE read_s4d_data_files(YEAR, MONTH, DAY, ldaily, lmonthly)

    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_bcast, p_io
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_transform_bi,     ONLY: locate_in_decomp
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: set_channel_newfile
    USE messy_main_tools,            ONLY: find_next_free_unit
    USE messy_main_channel_mem,      ONLY: dom_curid
    
    IMPLICIT NONE
    
    INTRINSIC :: TRIM, ADJUSTL, ASSOCIATED, SIZE

    ! I/O
    INTEGER, INTENT(IN)  :: YEAR, MONTH, DAY
    LOGICAL, INTENT(IN)  :: ldaily, lmonthly

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_s4d_data_files'
    ! 200 + YYYYMMDD + '.pos'
    CHARACTER(LEN=214)          :: fname
    CHARACTER(LEN=8)            :: dstr
    CHARACTER(LEN=6)            :: mstr
    LOGICAL                     :: lex
    INTEGER                     :: iou
    LOGICAL                     :: lskip
    INTEGER                     :: i,j
    INTEGER                     :: status

    track_loop: DO i=1, NTRACK

       IF (XTRACK(i)%domain_idx /= dom_curid) CYCLE

       ! GENERATE FILENAME
       SELECT CASE(XTRACK(i)%io%usw)
       CASE(-1)
          CYCLE
       CASE(0)
          IF (.NOT. ldaily) CYCLE
          WRITE(dstr,'(i4,i2.2,i2.2)') YEAR, MONTH, DAY
          fname = TRIM(ADJUSTL(XTRACK(i)%io%fname))//dstr//'.pos'
       CASE(1)
          IF (.NOT. lmonthly) CYCLE
          WRITE(mstr,'(i4,i2.2)') YEAR, MONTH
          fname = TRIM(ADJUSTL(XTRACK(i)%io%fname))//mstr//'.pos'
       END SELECT

       ! op_pj_20110616+
       ! RESET
       IF (ASSOCIATED(XTRACK(i)%pos)) THEN
          DEALLOCATE(XTRACK(i)%pos)
          NULLIFY(XTRACK(i)%pos)
       END IF

       IF (ASSOCIATED(XTRACK(i)%dec)) THEN
          DEALLOCATE(XTRACK(i)%dec)
          NULLIFY(XTRACK(i)%dec)
       END IF

       XTRACK(i)%cpos = 0      
       XTRACK(i)%npos = 0      
       XTRACK(i)%lt   = .FALSE.
       ! op_pj_20110616-

       ! READ NEW TRACK DATA
       lskip = .FALSE.
       IF (p_parallel_io) THEN
          WRITE(*,*) substr,': (D=',dom_curid &
               ,') CHECKING FOR FILE         (' &
               , i,') '''&
               &//TRIM(fname)//''' ...'
          INQUIRE(file=TRIM(fname), exist=lex)
          IF (lex) THEN
             ! FILE EXISTS: READ DATE
!!$             IF (ASSOCIATED(XTRACK(i)%pos)) THEN
!!$                DEALLOCATE(XTRACK(i)%pos)
!!$                NULLIFY(XTRACK(i)%pos)
!!$             END IF
             WRITE(*,*) substr,': (D=',dom_curid &
                  ,') IMPORTING TRACK DATA FROM ('&
                  , i,') '''&
                  &//TRIM(fname)//''' ...'
             iou = find_next_free_unit(100,200)
             CALL read_track(status, TRIM(fname), iou, XTRACK(i)%pos)
             IF (status /= 0) CALL error_bi('',substr)
             XTRACK(i)%npos = SIZE(XTRACK(i)%pos)
             WRITE(*,*) '   ... ',XTRACK(i)%npos,' positions read'      
          ELSE
             ! FILE DOES NOT EXIST:
             WRITE(*,*) substr,': FILE '''//TRIM(fname)//&
                  &''' does not exist ... skipping'
             lskip = .TRUE.
          END IF
       END IF

       ! BROADCAST RESULTS
       CALL p_bcast(lskip, p_io)
       IF (lskip) CYCLE

       ! BROADCAST TRACK DATA
       CALL p_bcast(XTRACK(i)%npos, p_io)
       IF (.NOT. p_parallel_io) THEN
!!$          IF (ASSOCIATED(XTRACK(i)%pos)) THEN
!!$             DEALLOCATE(XTRACK(i)%pos)
!!$             NULLIFY(XTRACK(i)%pos)
!!$          END IF
          ALLOCATE(XTRACK(i)%pos(XTRACK(i)%npos))
       END IF
       DO j=1, XTRACK(i)%npos
          CALL p_bcast(XTRACK(i)%pos(j)%yr,  p_io)
          CALL p_bcast(XTRACK(i)%pos(j)%mo,  p_io)
          CALL p_bcast(XTRACK(i)%pos(j)%dy,  p_io)
          CALL p_bcast(XTRACK(i)%pos(j)%ho,  p_io)
          CALL p_bcast(XTRACK(i)%pos(j)%mi,  p_io)
          CALL p_bcast(XTRACK(i)%pos(j)%se,  p_io)
          CALL p_bcast(XTRACK(i)%pos(j)%lon, p_io)
          CALL p_bcast(XTRACK(i)%pos(j)%lat, p_io)
          CALL p_bcast(XTRACK(i)%pos(j)%pre, p_io)
       END DO

       ! WORKSPACE
!!$       IF (ASSOCIATED(XTRACK(i)%dec)) THEN
!!$          DEALLOCATE(XTRACK(i)%dec)
!!$          NULLIFY(XTRACK(i)%dec)
!!$       END IF
       ALLOCATE(XTRACK(i)%dec(XTRACK(i)%npos))

       position_loop: DO j=1, XTRACK(i)%npos

          IF ( XTRACK(i)%io%linth ) THEN
             ! horizontal interpolation
             CALL locate_in_decomp(status                           &
                  , XTRACK(i)%pos(j)%lon,   XTRACK(i)%pos(j)%lat    &
                  , XTRACK(i)%dec(j)%pe(:), XTRACK(i)%dec(j)%jp(:)  &
                  , XTRACK(i)%dec(j)%jrow(:)                        &
                  , XTRACK(i)%dec(j)%w(:))
          ELSE
             ! nearest neighbour
             CALL locate_in_decomp(status                           &
                  , XTRACK(i)%pos(j)%lon,   XTRACK(i)%pos(j)%lat    &
                  , XTRACK(i)%dec(j)%pe(1), XTRACK(i)%dec(j)%jp(1)  &
                  , XTRACK(i)%dec(j)%jrow(1) )
             XTRACK(i)%dec(j)%pe(2:4)   = -1
             XTRACK(i)%dec(j)%jp(2:4)   = -1
             XTRACK(i)%dec(j)%jrow(2:4) = -1
             XTRACK(i)%dec(j)%w(1)      = 1.0_dp
             XTRACK(i)%dec(j)%w(2:4)    = 0.0_dp
          ENDIF

       END DO position_loop

       ! RESET
       XTRACK(i)%cpos = 0        ! current position
       XTRACK(i)%lt   = .TRUE.   ! still on track ...

       ! FORCE NEW OUTPUT FILE AFTER UPDATE
       CALL set_channel_newfile(status &
            , modstr//'_'//TRIM(XTRACK(i)%io%name), .TRUE.) 
       CALL channel_halt(substr, status)

    END DO track_loop

  END SUBROUTINE read_s4d_data_files
! ---------------------------------------------------------------------

! **********************************************************************
END MODULE messy_s4d_si
! **********************************************************************
