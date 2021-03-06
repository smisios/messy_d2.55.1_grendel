#include "messy_main_ppd_bi.inc"

! **********************************************************************
MODULE messy_sorbit_si
! **********************************************************************

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY
  USE messy_main_channel,       ONLY: STRLEN_OBJECT, STRLEN_CHANNEL
  USE messy_main_constants_mem, ONLY: FLAGGED_BAD, STRLEN_VLONG
  USE messy_sorbit

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  TYPE T_ORBIT_IO
     CHARACTER(LEN=STRLEN_VLONG) :: name = ''   ! satellite (orbit) name
     LOGICAL                  :: latc = .FALSE. ! latitude correction
     REAL(DP)                 :: incl = 98.2_dp ! orbit: inclination [deg]
     LOGICAL                  :: lod  = .FALSE. ! TRUE:  ascending  (+1)
     !                                            FALSE: descending (-1)
     INTEGER                  :: LTH = 0     ! local time (hour)
     INTEGER                  :: LTM = 0     ! local time (minute)
     LOGICAL                  :: limit = .FALSE. ! limit dt to LT box-distance
     CHARACTER(LEN=STRLEN)    :: str = ''    ! list of cha/obj
  END TYPE T_ORBIT_IO

  TYPE T_ORBIT
     TYPE(T_ORBIT_IO)              :: io
     !
     REAL(DP)                      :: direction ! 1: ascending, -1: descending
     ! second of day (local time)
     INTEGER,                       DIMENSION(:,:), POINTER :: ltsofd => NULL()
     REAL(DP),                      DIMENSION(:,:), POINTER :: lsthr  => NULL()
     INTEGER                       :: deltat
     !
     INTEGER                       :: nobj        ! number of objects
     CHARACTER(LEN=STRLEN_CHANNEL), DIMENSION(:), POINTER :: cha  => NULL()
     CHARACTER(LEN=STRLEN_OBJECT),  DIMENSION(:), POINTER :: obj  => NULL()
     LOGICAL,                       DIMENSION(:), POINTER :: lex  => NULL()
     ! POINTER TO DATA
     TYPE(PTR_3D_ARRAY),            DIMENSION(:), POINTER :: dat  => NULL()
     ! REPRESENTATION ID
     INTEGER,                       DIMENSION(:), POINTER :: rid  => NULL()
     ! POINTER FOR SAMPLED DATA
     TYPE(PTR_3D_ARRAY),            DIMENSION(:), POINTER :: smp  => NULL()
     !
     INTEGER                                              :: domain_idx
  END TYPE T_ORBIT

  ! WORKSPACE
  INTEGER, PARAMETER :: spd = 86400 ! seconds per day
  LOGICAL            :: lout_now = .FALSE.

  ! NAMELIST
  LOGICAL            :: lout_auto = .TRUE.
  REAL(DP)           :: r_init = FLAGGED_BAD

  !relaxes the local time condition
  REAL(DP)           :: deltat_scale = 1._dp

  INTEGER, PARAMETER :: NMAXORB = 50
  TYPE(T_ORBIT_IO),     DIMENSION(NMAXORB),    SAVE :: ORB  ! CPL
  TYPE(T_ORBIT),        DIMENSION(:), POINTER, SAVE :: XORB => NULL()
  INTEGER,                                     SAVE :: NORB

  PUBLIC :: sorbit_initialize
  PUBLIC :: sorbit_init_coupling
  PUBLIC :: sorbit_global_start
  PUBLIC :: sorbit_write_output
  PUBLIC :: sorbit_free_memory
  !PRIVATE :: sorbit_read_nml_cpl
  
CONTAINS

! ======================================================================
! PUBLIC ROUTINES
! ======================================================================

! ----------------------------------------------------------------------
  SUBROUTINE sorbit_initialize

    ! SORBIT MODULE ROUTINE (ECHAM5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2007

    USE messy_main_channel_bi, ONLY: n_dom
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: str2chob, find_next_free_unit &
                                   , domains_from_string

    IMPLICIT NONE

    INTRINSIC :: TRIM, ADJUSTL, LEN_TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'sorbit_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i, j
    CHARACTER(LEN=32)                        :: varname = ''
    INTEGER,           DIMENSION(:), POINTER :: domnum  => NULL()
    INTEGER                                  :: nd, num

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL sorbit_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('', substr)
    END IF

    CALL p_bcast(lout_auto, p_io)
    CALL p_bcast(r_init, p_io)    
    CALL p_bcast(deltat_scale, p_io)

    CALL start_message_bi(modstr,'OUTPUT INITIALISATION ',substr)

    IF (p_parallel_io) THEN
       NORB = 1
       DO i=1, NMAXORB
          IF (TRIM(ORB(i)%name) == '') CYCLE
          CALL domains_from_string(status,ORB(i)%name,n_dom,num)
          NORB = NORB + num
       END DO
       NORB = NORB - 1

       ALLOCATE(XORB(NORB))

       ! GET NUMBER OF LOCATIONS
       NORB = 1
       ! COPY DATA AND PARSE STR
       orbit_loop: DO i=1, NMAXORB
          IF (TRIM(ORB(i)%name) == '') CYCLE
          CALL domains_from_string(status,ORB(i)%name,n_dom,num &
                ,varname,dnums=domnum )
          IF (status /= 0) CALL error_bi('error in namelist parsing', substr)
          IF (LEN_TRIM(ADJUSTL(varname)) > 12) &
               CALL error_bi('orbit name too long (max. 8(+4) characters)' &
               , substr)
          domain_loop: DO nd = 1, SIZE(domnum)
             XORB(NORB)%io%name    = TRIM(ADJUSTL(varname))
             XORB(NORB)%domain_idx = domnum(nd)
             XORB(NORB)%io%lth   = ORB(i)%lth
             XORB(NORB)%io%ltm   = ORB(i)%ltm
             XORB(NORB)%io%limit = ORB(i)%limit
             XORB(NORB)%io%latc  = ORB(i)%latc
             XORB(NORB)%io%incl  = ORB(i)%incl
             XORB(NORB)%io%lod   = ORB(i)%lod
             XORB(NORB)%io%str   = TRIM(ORB(i)%str)

             WRITE(*,*) 'ORBIT NAME     : ',TRIM(XORB(NORB)%io%name)
             WRITE(*,*) '  EQ-LT HH:MM  : ',ORB(i)%lth,':',ORB(i)%ltm
             IF (XORB(NORB)%io%limit) THEN
                WRITE(*,*) '  LIMIT DELTA-t: YES'
             ELSE
                WRITE(*,*) '  LIMIT DELTA-t: NO'
             END IF
             IF (XORB(NORB)%io%latc) THEN
                WRITE(*,*) '  LAT. CORRECT.: YES'
                WRITE(*,*) '  INCLINATION  : ',XORB(NORB)%io%incl,' deg'
                IF (XORB(NORB)%io%lod) THEN
                   WRITE(*,*) '  DIRECTION    : ASCENDING  (+1)'
                ELSE
                   WRITE(*,*) '  DIRECTION    : DESCENDING (-1)'
                END IF
             ELSE
                WRITE(*,*) '  LAT. CORRECT.: NO'
             ENDIF
             
             ! SET PRELIMINARY %nobj
             CALL str2chob(status, XORB(NORB)%io%str, XORB(NORB)%nobj &
                  , XORB(NORB)%cha, XORB(NORB)%obj)
             IF (status /= 0) THEN
                WRITE(*,*) '   ... ERROR IN STRING ... skipping'
                CYCLE
             END IF

             IF (XORB(NORB)%nobj == 0) THEN
                WRITE(*,*) '   ... EMPTY OUTPUT LIST ... skipping'
                CYCLE
             ELSE
                WRITE(*,*) '  REQUESTS     : ',XORB(NORB)%nobj
             END IF
             !
             WRITE(*,'(1x,a16,1x,a32)') 'CHANNEL','OBJECT(S)'
             WRITE(*,'(1x,a16,1x,a32)') '-------','---------'
             DO j=1, XORB(NORB)%nobj 
                WRITE(*,'(1x,a16,1x,a32)') TRIM(XORB(NORB)%cha(j)) &
                     , TRIM(XORB(NORB)%obj(j))
             END DO

             ! NEXT ORB
             NORB = NORB + 1
             WRITE(*,*) '------------------------------------------------------'
          END DO domain_loop
          DEALLOCATE(domnum) ; NULLIFY(domnum)
       END DO orbit_loop
       NORB = NORB - 1 
       IF (NORB > SIZE(XORB)) CALL error_bi('error parsing namelist',substr)
    END IF
    CALL p_bcast(NORB, p_io)
    IF (.NOT. p_parallel_io) ALLOCATE(XORB(NORB))

    ! BROADCAST ALL RESULTS
    DO i=1, NORB
       ! I/O USER INTERFACE
       CALL p_bcast(XORB(i)%io%name, p_io)
       CALL p_bcast(XORB(i)%io%lth,  p_io)
       CALL p_bcast(XORB(i)%io%ltm,  p_io)
       CALL p_bcast(XORB(i)%io%limit, p_io)
       CALL p_bcast(XORB(i)%io%latc,  p_io)
       CALL p_bcast(XORB(i)%io%incl,  p_io)
       CALL p_bcast(XORB(i)%io%lod,   p_io)
       CALL p_bcast(XORB(i)%io%str,  p_io)
       !
       ! CHANNEL/OBJECT NAMES
       CALL p_bcast(XORB(i)%nobj, p_io)
       IF (.NOT. p_parallel_io) THEN
          ALLOCATE(XORB(i)%cha(XORB(i)%nobj))
          ALLOCATE(XORB(i)%obj(XORB(i)%nobj))
       END IF
       DO j=1, XORB(i)%nobj
          CALL p_bcast(XORB(i)%cha(j), p_io)
          CALL p_bcast(XORB(i)%obj(j), p_io)
       END DO
       !
       ! NOTE: THE INITIALISATION OF THE LOCAL TIME FIELD IS EXECUTED IN
       !       sorbit_init_coupling, SINCE HERE philat_2d IS NOT YET AVAILABLE
       !
       CALL p_bcast(XORB(i)%domain_idx, p_io)
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NORB,' ORBIT(S) INITIALIZED !'
    END IF

    CALL end_message_bi(modstr,'OUTPUT INITIALISATION ',substr)

  END SUBROUTINE sorbit_initialize
! ----------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE sorbit_init_coupling

    ! SORBIT MODULE ROUTINE
    !
    ! define specific channel(s) and allocate memory for
    ! global fields
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2007

    USE messy_main_blather_bi,       ONLY: info_bi, error_bi, warning_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, GP_3D_INT, GP_3D_1LEV &
                                   , GP_2D_HORIZONTAL
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, ngpblks
    USE messy_main_grid_def_bi,      ONLY: philat_2d
#ifdef ECHAM5
    USE messy_main_grid_def_mem_bi,          ONLY: nlon
#endif
    USE messy_main_timer,            ONLY: lstart, delta_time

    USE messy_main_channel,       ONLY: new_channel, new_channel_object   &
                                      , new_attribute, get_channel_object &
                                      , get_attribute, get_channel_info   &
                                      , get_channel_object_info
    USE messy_main_channel_mem,   ONLY: dom_curid
    USE messy_main_constants_mem, ONLY: STRLEN_ULONG
    USE messy_main_channel_repr,  ONLY: REPR_UNDEF
    USE messy_main_tools,         ONLY: match_wild

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE, TRIM, NINT, REAL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'sorbit_init_coupling'
    INTEGER                      :: status
    INTEGER                      :: reprid
    INTEGER                      :: reprid_new
    INTEGER                      :: i, jr, n, nptr, jo
    LOGICAL                      :: lstatic
    CHARACTER(LEN=STRLEN_ULONG)  :: charatt
    CHARACTER(LEN=STRLEN_OBJECT), DIMENSION(:), POINTER :: ONAMES => NULL()

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    orbit_loop: DO i=1, NORB
       IF (dom_curid /= XORB(i)%domain_idx) CYCLE

       CALL info_bi('        '//TRIM(XORB(i)%io%name)//': ', substr)

       ! OPEN ONE CHANNEL PER ORB
       CALL new_channel(status, modstr//'_'//TRIM(XORB(i)%io%name))
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//TRIM(XORB(i)%io%name) &
            , 'sorbit_lteq_hour', i=XORB(i)%io%lth )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//TRIM(XORB(i)%io%name) &
            , 'sorbit_lteq_minute', i=XORB(i)%io%ltm )
       CALL channel_halt(substr, status)

       ! DEFAULT (XORB(i)%io%limit = .FALSE.)
       XORB(i)%deltat = NINT(delta_time*deltat_scale) / 2
       !
       IF (XORB(i)%io%limit) THEN
#ifdef ECHAM5
          ! This limit must be  (spd/360)*dlon, 
          ! where dlon is the width of the box (in degrees);
          ! i.e., for ECHAM5 dlon=360/nlon:
          XORB(i)%deltat = MIN(NINT(delta_time*deltat_scale), spd/nlon) / 2
#endif
#ifdef COSMO
          ! For COSMO this does not make sense; therefore we
          ! use the time step length as limit in both cases!
          CALL warning_bi('LIMITATION OF DT HAS NO EFFECT IN COSMO ('&
               &//TRIM(XORB(i)%io%name)//')', substr)
          XORB(i)%deltat = NINT(delta_time*deltat_scale) / 2
#endif
       END IF
       !
       CALL new_attribute(status, modstr//'_'//TRIM(XORB(i)%io%name) &
            , 'sorbit_dt', i=XORB(i)%deltat )
       CALL channel_halt(substr, status)

       IF (XORB(i)%io%latc) THEN
          CALL new_attribute(status, modstr//'_'//TRIM(XORB(i)%io%name) &
               , 'sorbit_inclination', r=XORB(i)%io%incl )
          IF (XORB(i)%io%lod) THEN
             CALL new_attribute(status, modstr//'_'//TRIM(XORB(i)%io%name) &
                  , 'sorbit_direction', c='ascending' )
          ELSE
             CALL new_attribute(status, modstr//'_'//TRIM(XORB(i)%io%name) &
                  , 'sorbit_direction', c='descending' )
          END IF
       ENDIF

       ! CALCULATE REQUIRED NUMBER OF POINTERS
       nptr = 0
       DO jr=1,  XORB(i)%nobj  ! LOOP OVER REQUESTS
          CALL get_channel_info(status, TRIM(XORB(i)%cha(jr)) &
               , ONAMES = ONAMES)
          IF (status /= 3003) THEN ! CHANNEL (NAME) DOES NOT EXIST
             CALL channel_halt(substr, status)
          ELSE
             CALL warning_bi(' ... channel '''&
                  &//TRIM(XORB(i)%cha(jr))//''' does not exist ... skipping' &
                  , substr)
             CYCLE
          END IF
          DO jo = 1, SIZE(ONAMES)
             IF ( match_wild(TRIM(XORB(i)%obj(jr)), TRIM(ONAMES(jo))) ) THEN
                nptr = nptr + 1
             END IF
          END DO
       END DO
       !
       IF (ASSOCIATED(ONAMES)) THEN
          DEALLOCATE(ONAMES)
          NULLIFY(ONAMES)
       END IF

       ! ALLOCATE SPACE FOR POINTER ARRAYS (DATA, REPR-ID, COLUMN OBJECT)
       ALLOCATE(XORB(i)%dat(nptr))   ! POINTER TO CHANNEL OBJECT
       ALLOCATE(XORB(i)%lex(nptr))   ! channel / object exsists
       XORB(i)%lex(:) = .FALSE.
       ALLOCATE(XORB(i)%rid(nptr))   ! REPRESENTATION OF CHANNEL OBJECT
       XORB(i)%rid(:) = REPR_UNDEF
       ALLOCATE(XORB(i)%smp(nptr))   ! POINTER TO SAMPLED CHANNEL OBJECT

       ! RESET NUMBER OF REQUESTS
       n = XORB(i)%nobj    ! SAVE NUMBER OF REQUESTS FOR LOOP BELOW       
       XORB(i)%nobj = nptr ! NUMBER OF MATCHING OBJECTS
       ! (SO FAR XORB(i)%nobj CONTAINED NUMBER OF REQUESTS)

       ! ADD REQUESTED CHANNEL OBJECTS
       nptr = 0
       request_loop: DO jr=1, n  ! LOOP OVER REQUESTS

          ! GET ALL POTENTIAL OBJECTS
          CALL get_channel_info(status, TRIM(XORB(i)%cha(jr)) &
               , ONAMES = ONAMES)
          IF (status == 3003) THEN ! CHANNEL (NAME) DOES NOT EXIST
             CYCLE
          ELSE
             CALL channel_halt(substr, status)
          END IF

          object_loop: DO jo = 1, SIZE(ONAMES)

             ! CHECK IF POTENTIAL OBJECT IS MATCHING
             IF ( match_wild(TRIM(XORB(i)%obj(jr)), TRIM(ONAMES(jo))) ) THEN
                nptr = nptr + 1      ! NEXT POINTER INDEX
             ELSE
                CYCLE
             END IF
             
             CALL info_bi('          '//TRIM(XORB(i)%cha(jr))//' - '//&
                  &TRIM(ONAMES(jo)) &
                  , substr)

             ! INIT
             NULLIFY(XORB(i)%dat(nptr)%ptr)
             NULLIFY(XORB(i)%smp(nptr)%ptr)

             CALL get_channel_object(status &
                  , TRIM(XORB(i)%cha(jr)), TRIM(ONAMES(jo)) &
                  , p3=XORB(i)%dat(nptr)%ptr )
             IF (status /= 0) THEN
                CALL warning_bi( &
                     ' ... channel object '''//TRIM(ONAMES(jo)) &
                     &//''' does not exist ... skipping', substr)
                NULLIFY(XORB(i)%dat(nptr)%ptr)
                NULLIFY(XORB(i)%smp(nptr)%ptr)
                CYCLE
             END IF

             CALL get_channel_object_info(status &
                  , TRIM(XORB(i)%cha(jr)), TRIM(ONAMES(jo)) &
                  , reprid=reprid, lstatic=lstatic)
             CALL channel_halt(substr, status)
             !
             XORB(i)%rid(nptr) = reprid
             !
             IF ( (reprid == GP_3D_MID)       .OR. &
                  (reprid == GP_3D_INT)       .OR. &
                  (reprid == GP_3D_1LEV)      .OR. &
                  (reprid == GP_2D_HORIZONTAL) ) THEN
                reprid_new = reprid
             ELSE
                CALL warning_bi( &
                  '            ... representation not supported ... skipping' &
                  , substr)
                NULLIFY(XORB(i)%dat(nptr)%ptr)
                NULLIFY(XORB(i)%smp(nptr)%ptr)
                CYCLE
             END IF

             CALL new_channel_object(status &
                  , modstr//'_'//TRIM(XORB(i)%io%name) &
                  , TRIM((XORB(i)%cha(jr)))//'_'//TRIM(ONAMES(jo))  &
                  , p3=XORB(i)%smp(nptr)%ptr &
                  , reprid=reprid_new, lrestreq=.TRUE., lstatic=lstatic)
             CALL channel_halt(substr, status)

             ! COPY ATTRIBUTES
             CALL get_attribute(status &
                  , TRIM(XORB(i)%cha(jr)), TRIM(ONAMES(jo)) &
                  , 'long_name', c=charatt)
             IF (status == 0) THEN
                CALL new_attribute(status &
                     , modstr//'_'//TRIM(XORB(i)%io%name) &
                     , TRIM((XORB(i)%cha(jr)))//'_'//TRIM(ONAMES(jo))  &
                     , 'long_name', c=TRIM(charatt))
                CALL channel_halt(substr, status)
             END IF
             CALL get_attribute(status &
                  , TRIM(XORB(i)%cha(jr)), TRIM(ONAMES(jo)) &
                  , 'units', c=charatt)
             IF (status == 0) THEN
                CALL new_attribute(status &
                     , modstr//'_'//TRIM(XORB(i)%io%name) &
                     , TRIM((XORB(i)%cha(jr)))//'_'//TRIM(ONAMES(jo))  &
                     , 'units', c=TRIM(charatt))
                CALL channel_halt(substr, status)
             END IF

             XORB(i)%lex(nptr) = .TRUE.

             IF (lstart) THEN
                XORB(i)%smp(nptr)%ptr(:,:,:) = r_init
             END IF

          END DO object_loop

       END DO request_loop

       IF (nptr /= XORB(i)%nobj) THEN
          CALL error_bi('something went wrong with pointer counting', substr)
       END IF

       IF (ASSOCIATED(ONAMES)) DEALLOCATE(ONAMES)

       ! NOTE: THE FOLLOWING OPERATION MUST BE EXECUTED HERE, SINCE IN
       !       messy_initialize philat_2d IS NOT YET AVAILABLE
       IF (XORB(i)%io%lod) THEN
          XORB(i)%direction = 1.0_dp
       ELSE
          XORB(i)%direction = -1.0_dp
       END IF
       !
       ALLOCATE(XORB(i)%ltsofd(nproma, ngpblks))
       !
       IF (XORB(i)%io%latc) THEN
          !
          CALL new_channel_object(status &
               , modstr//'_'//TRIM(XORB(i)%io%name) &
               , 'lsthr', p2=XORB(i)%lsthr &
               , reprid=GP_2D_HORIZONTAL, lrestreq=.FALSE.)
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status &
               , modstr//'_'//TRIM(XORB(i)%io%name), 'lsthr' &
               , 'long_name', c='local solar time')
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status &
               , modstr//'_'//TRIM(XORB(i)%io%name), 'lsthr' &
               , 'units', c='hour of day')
          CALL channel_halt(substr, status)
          !
          ! calculate latitude dependent local solar time
          ! - time at equator given in hh:mi, converted to hour of day
          CALL polst(XORB(i)%lsthr, XORB(i)%io%incl, XORB(i)%direction  &
               , REAL(XORB(i)%io%lth,dp)+REAL(XORB(i)%io%ltm, DP)/60_dp &
               , philat_2d)
          ! convert from hour of day to second of day
          XORB(i)%ltsofd(:,:) = NINT(XORB(i)%lsthr(:,:)*3600._dp)
       ELSE
          XORB(i)%ltsofd(:,:) = XORB(i)%io%lth * 3600 + XORB(i)%io%ltm * 60
       END IF

    END DO orbit_loop

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE sorbit_init_coupling
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE sorbit_global_start

    USE messy_main_timer,            ONLY: HOUR, MINUTE, SECOND, delta_time
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: set_channel_output
    USE messy_main_channel_mem,      ONLY: dom_curid
    USE messy_main_blather_bi,       ONLY: info_bi

    IMPLICIT NONE

    INTRINSIC :: NINT

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'sorbit_global_start'
    INTEGER :: utsofd  ! universal time (second of day)
    INTEGER :: i, jo
    INTEGER :: status

    ! universal time (seconds of day)
    utsofd = HOUR*3600 + MINUTE*60 + SECOND

    ! ALLOW/FORCE OUTPUT AT LAST TIME STEP OF DAY ONLY
    lout_now =  ( (utsofd + NINT(delta_time)) >= spd)

    ! TRIGGER/SUPPRESS CHANNEL OUTPUT
    IF (lout_auto) THEN
       IF (lout_now) THEN
          CALL info_bi('TRIGGER DATA OUTPUT ...', substr)
       !ELSE
       !   CALL info_bi('SUPPRESS DATA OUTPUT ...', substr)
       END IF
       !
       DO i=1, NORB
          IF (dom_curid /= XORB(i)%domain_idx) CYCLE
          CALL set_channel_output(status &
               , modstr//'_'//TRIM(XORB(i)%io%name), lout_now)
          CALL channel_halt(substr, status)
       END DO
    END IF

    ! RESET ONLY AT FIRST TIME STEP OF DAY
    IF ( utsofd >= NINT(delta_time) ) RETURN
    !
    CALL info_bi('RESET DATA FIELD ...', substr)
    !
    orbit_loop: DO i = 1, NORB
       IF (dom_curid /= XORB(i)%domain_idx) CYCLE
       object_loop: DO jo=1, XORB(i)%nobj
          IF (XORB(i)%lex(jo)) THEN
             XORB(i)%smp(jo)%ptr(:,:,:) = r_init          
          END IF
       END DO object_loop
    END DO orbit_loop

  END SUBROUTINE sorbit_global_start
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE sorbit_write_output

    USE messy_main_timer,       ONLY: HOUR, MINUTE, SECOND
    USE messy_main_channel_mem, ONLY: dom_curid
    USE messy_main_grid_def_mem_bi, ONLY: nproma, npromz, ngpblks
    USE messy_main_grid_def_bi,     ONLY: philon_2d
    USE messy_main_channel_bi,  ONLY: GP_2D_HORIZONTAL

    IMPLICIT NONE

    INTRINSIC :: NINT, REAL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'sorbit_write_output'
    INTEGER                     :: utsofd  ! second of day (universal time)
    INTEGER                     :: i
    INTEGER                     :: jo
    INTEGER                     :: jrow, kproma, jp
    INTEGER                     :: ltsofd     ! second of day (local time) 
    INTEGER                     :: nd
    INTEGER                     :: ll, ul  ! lower and upper limit LT [SofD]
    LOGICAL                     :: lodd, lcl1, lcl2, lcu1, lcu2

    ! universal time (seconds of day)
    utsofd = HOUR*3600 + MINUTE*60 + SECOND

    ! CALCULATE LOCAL TIME AND FILL FLAG
    row: DO jrow = 1, ngpblks
#ifndef CESM1
       IF (jrow == ngpblks) THEN
          kproma = npromz
       ELSE
          kproma = nproma
       END IF
#else
        kproma = npromz(jrow)
#endif
       orbit_loop: DO i=1, NORB
          IF (dom_curid /= XORB(i)%domain_idx) CYCLE

          vector: DO jp=1, kproma

             ! local solar time (second of day) at center of longitude interval
             ltsofd = utsofd + &
                  NINT((philon_2d(jp,jrow)/360._dp) * REAL(spd,dp))
             nd = ltsofd / spd ! = 0, 1
             ltsofd = ltsofd - nd * spd

             ! LT [second of day] IS MODULO ...
             ! ... lower limit
             ll = ltadd(ltsofd, -1*XORB(i)%deltat) ! [0,spd]
             ! ... upper limit
             ul = ltadd(ltsofd,    XORB(i)%deltat) ! [0,spd]
             ! 
             lodd = (ll > ul)  ! modulo
             !
             lcl1 = XORB(i)%ltsofd(jp,jrow) > ll
             lcu1 = XORB(i)%ltsofd(jp,jrow) < ul
             lcl2 = (XORB(i)%ltsofd(jp,jrow) + spd) > ll
             lcu2 = (XORB(i)%ltsofd(jp,jrow) - spd) < ul
             
             IF (.NOT. ( &
                  ( lcl1 .AND. lcu1 )            .OR. &
                  ( lcl2 .AND. lcu1 .AND. lodd ) .OR. &
                  ( lcl1 .AND. lcu2 .AND. lodd ) &
                  ) ) CYCLE

             object_loop: DO jo=1, XORB(i)%nobj
                
                IF (.NOT. XORB(i)%lex(jo)) CYCLE
                
                IF (XORB(i)%rid(jo) == GP_2D_HORIZONTAL) THEN
                   XORB(i)%smp(jo)%ptr(jp, jrow, 1) = &
                        XORB(i)%dat(jo)%ptr(jp, jrow, 1)
                ELSE
                   XORB(i)%smp(jo)%ptr(_RI_XYZ__(jp,jrow,:)) = &
                        XORB(i)%dat(jo)%ptr(_RI_XYZ__(jp,jrow,:))
                END IF
                
             END DO object_loop
             
          END DO vector

       END DO orbit_loop
       
    END DO row

  END SUBROUTINE sorbit_write_output
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE sorbit_free_memory

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! LOCAL
    INTEGER :: i
    
    orb_loop: DO i=1, NORB
       !
       IF (ASSOCIATED(XORB(i)%ltsofd)) DEALLOCATE(XORB(i)%ltsofd)
       NULLIFY(XORB(i)%ltsofd)
       IF (ASSOCIATED(XORB(i)%cha)) DEALLOCATE(XORB(i)%cha)
       NULLIFY(XORB(i)%cha)
       IF (ASSOCIATED(XORB(i)%obj)) DEALLOCATE(XORB(i)%obj)
       NULLIFY(XORB(i)%obj)
       !
       IF (ASSOCIATED(XORB(i)%dat)) DEALLOCATE(XORB(i)%dat)
       NULLIFY(XORB(i)%dat)
       IF (ASSOCIATED(XORB(i)%rid)) DEALLOCATE(XORB(i)%rid)
       NULLIFY(XORB(i)%rid)
       IF (ASSOCIATED(XORB(i)%smp)) DEALLOCATE(XORB(i)%smp)
       NULLIFY(XORB(i)%smp)
       IF (ASSOCIATED(XORB(i)%lex)) DEALLOCATE(XORB(i)%lex)
       NULLIFY(XORB(i)%lex)
       !
    END DO orb_loop

    IF (ASSOCIATED(XORB)) DEALLOCATE(XORB)
    NULLIFY(XORB)
    
  END SUBROUTINE sorbit_free_memory
! ---------------------------------------------------------------------

! ======================================================================
! PRIVATE ROUTINES
! ======================================================================

! ----------------------------------------------------------------------
  SUBROUTINE sorbit_read_nml_cpl(status, iou)

    ! SORBIT MODULE ROUTINE (ECHAM5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling'
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2007

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'sorbit_read_nml_cpl'

    NAMELIST /CPL/ lout_auto, r_init, ORB, deltat_scale

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

  END SUBROUTINE sorbit_read_nml_cpl
! ----------------------------------------------------------------------

! **********************************************************************
END MODULE messy_sorbit_si
! **********************************************************************

