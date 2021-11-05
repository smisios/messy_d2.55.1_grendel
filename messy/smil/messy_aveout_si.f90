! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL AVEOUT 
!
! Author : Astrid Kerkweg, University Bonn, 2018
!
! References: 
!
! **********************************************************************

! **********************************************************************
MODULE messy_aveout_si
! **********************************************************************

  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      error_bi, warning_bi

  ! SMCL
  USE messy_main_timer_event,   ONLY: time_event, io_time_event
  USE messy_main_channel,       ONLY: t_chaobj_cpl
  USE messy_aveout

  IMPLICIT NONE
  INTRINSIC :: NULL
  SAVE
  PRIVATE

  TYPE T_AVEOUT_IO
     TYPE(t_chaobj_cpl)        :: Fld
     TYPE(io_time_event)       :: &
        AVE_IOEVENT = io_time_event(6,'hours','first',0)
  END TYPE T_AVEOUT_IO

  INTEGER,           PARAMETER          :: MAX_OUT = 2000
  TYPE(T_AVEOUT_IO), DIMENSION(MAX_OUT) :: OUTPUT

  TYPE T_AVEOUT
     INTEGER                               :: domain_idx  = 0
     TYPE(t_chaobj_cpl)                    :: Fld
     TYPE(io_time_event)                   :: &
          AVE_IOEVENT = io_time_event(6,'hours','first',0)
     TYPE(time_event)                      :: AVE_EVENT
     REAL(dp), DIMENSION(:,:,:,:), POINTER :: ptr    => NULL()
     REAL(dp), DIMENSION(:,:,:,:), POINTER :: ptrout => NULL()
  END TYPE T_AVEOUT

  TYPE(T_AVEOUT), DIMENSION(:), ALLOCATABLE :: OUT
  INTEGER                                   :: num_out = 0


  ! PUBLIC SUBROUTINES (called from messy_main_control_si.f90)
  ! NOTE: in case you activate further entry points, make sure to call them
  !       in messy_main_control_si.f90
  PUBLIC :: aveout_initialize    ! initialize submodel
  PUBLIC :: aveout_init_coupling ! set pointers for coupling to BM and other SMs
  PUBLIC :: aveout_global_end    ! entry point in time loop (all vectors)

  ! PRIVATE SUBROUTINES
  !PRIVATE :: aveout_read_nml_cpl

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE aveout_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_channel_bi, ONLY: n_dom
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_timer_bi,   ONLY: p_bcast_event, timer_event_init
    USE messy_main_tools,      ONLY: find_next_free_unit, str &
                                   , domains_from_string

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'aveout_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    INTEGER                     :: ix
    CHARACTER(LEN=32)                        :: varname = '  '
    INTEGER,           DIMENSION(:), POINTER :: domnum  => NULL()
    INTEGER                                  :: nd, num

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

    ! READ CPL namelist
    IF (p_parallel_io) THEN                   ! read only on I/O-PE

       iou = find_next_free_unit(100,200)     ! find next free I/O unit
       CALL aveout_read_nml_cpl(status, iou)  ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)

       NUM_OUT = 0
       DO ix = 1,MAX_OUT
         IF (ADJUSTL(TRIM(OUTPUT(ix)%Fld%cha)) == '') CYCLE
         IF (ADJUSTL(TRIM(OUTPUT(ix)%Fld%obj)) == '') CYCLE
         CALL domains_from_string(status,OUTPUT(ix)%Fld%cha,n_dom,num)
         NUM_OUT = NUM_OUT + num
       END DO

       ALLOCATE(OUT(NUM_OUT))
    
       NUM_OUT = 0
       DO ix = 1, MAX_OUT
          IF (ADJUSTL(TRIM(OUTPUT(ix)%Fld%cha)) == ' ') CYCLE
          IF (ADJUSTL(TRIM(OUTPUT(ix)%Fld%obj)) == '')  CYCLE
          CALL domains_from_string(status,OUTPUT(ix)%Fld%cha,n_dom,num &
               ,varname,dnums=domnum )

          DO nd = 1, SIZE(domnum)
             NUM_OUT = NUM_OUT + 1
             OUT(NUM_OUT)%Fld%cha     = ADJUSTL(TRIM(varname))
             OUT(NUM_OUT)%domain_idx  = domnum(nd)
             OUT(NUM_OUT)%Fld%obj     = ADJUSTL(TRIM(OUTPUT(ix)%Fld%obj))
             OUT(NUM_OUT)%AVE_IOEVENT = OUTPUT(ix)%AVE_IOEVENT
             write (*,*) 'OUTPUT AVERAGE FIELD ' , ix, nd &
                  , TRIM( OUT(NUM_OUT)%Fld%cha), TRIM( OUT(NUM_OUT)%Fld%obj) &
                  , OUT(NUM_OUT)%AVE_IOEVENT
          END DO
          DEALLOCATE(domnum); NULLIFY(domnum)
       END DO
       IF  (NUM_OUT > SIZE(OUT)) CALL error_bi('error parsing namelist', substr)
    END IF
    CALL p_bcast(NUM_OUT, p_io)

    IF (.NOT. p_parallel_io) THEN
       ALLOCATE(OUT(NUM_OUT))
    END IF
    ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs
    DO ix = 1, NUM_OUT
       CALL p_bcast(OUT(ix)%Fld%cha, p_io)
       CALL p_bcast(OUT(ix)%Fld%obj, p_io)
       CALL p_bcast(OUT(ix)%domain_idx, p_io)
       CALL p_bcast_event(OUT(ix)%AVE_IOEVENT,  p_io)

       CALL timer_event_init(OUT(ix)%AVE_EVENT, OUT(ix)%AVE_IOEVENT,&
            'AVEOUT_EVENT'//str(ix), 'present')
    END DO

    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE aveout_initialize
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE aveout_init_coupling

    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the 
    ! basemodel and to other submodels.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_channel_error_bi,  ONLY: channel_halt
    USE messy_main_channel,     ONLY: get_channel_object      &
                                    , get_channel_object_info &
                                    , new_channel             &
                                    , new_channel_object
    USE messy_main_channel_mem, ONLY: dom_curid
    USE messy_main_tools,       ONLY: str
    USE messy_main_timer_event, ONLY: event_count, event_unit, event_offset
    
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'aveout_init_coupling'
    INTEGER          :: ix
    INTEGER          :: rep_id
    INTEGER          :: status
    CHARACTER(LEN=7) :: ustr
    CHARACTER(LEN=4) :: cstr
    CHARACTER(LEN=5) :: ostr

    CALL new_channel(status, modstr, lrestreq=.TRUE.)

    DO ix = 1, NUM_OUT
       IF (OUT(ix)%domain_idx /= dom_curid ) CYCLE
       NULLIFY(OUT(ix)%ptr)
       NULLIFY(OUT(ix)%ptrout)
       CALL get_channel_object (status, TRIM(OUT(ix)%Fld%Cha) &
            , TRIM(OUT(ix)%Fld%Obj), p4=OUT(ix)%ptr)

       IF (status /= 0) THEN
          CALL warning_bi('Field '//TRIM(OUT(ix)%Fld%Cha)//' '//TRIM(OUT(ix)%Fld%Obj)//' not available .. skipping', substr)
       ELSE
          
          CALL get_channel_object_info (status, TRIM(OUT(ix)%Fld%Cha) &
               , TRIM(OUT(ix)%Fld%Obj), reprid=rep_id)
          CALL channel_halt(substr, status)

          ! Define object name based on output interval
          cstr = TRIM(ADJUSTL(str(event_count(OUT(ix)%AVE_EVENT))))
          ustr = event_unit(OUT(ix)%AVE_EVENT)

          CALL new_channel_object (status, modstr &
            , TRIM(OUT(ix)%Fld%Obj)//'_'//TRIM(cstr)//ustr(1:2)&
            , p4=OUT(ix)%ptrout, reprid=rep_id)
          IF (status == 3102) THEN ! CHANNEL OBJECT EXISTS ALREADY
             ostr=str(event_offset(OUT(ix)%AVE_EVENT))
            CALL new_channel_object (status, modstr &
            , TRIM(OUT(ix)%Fld%Obj)//'_'//TRIM(cstr)//ustr(1:2)//'off'//ostr&
            , p4=OUT(ix)%ptrout, reprid=rep_id)
          END IF
          CALL channel_halt(substr, status)
       END IF

    END DO

  END SUBROUTINE aveout_init_coupling
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE aveout_global_end

    ! MESSY
    USE messy_main_timer_bi,      ONLY: event_state
    USE messy_main_timer,         ONLY: current_date, lstart
    USE messy_main_channel_mem,   ONLY: dom_curid
    
    IMPLICIT NONE
 
    INTEGER :: ix 
    LOGICAL :: levent
 
    DO ix = 1, NUM_OUT
       IF (OUT(ix)%domain_idx /= dom_curid ) CYCLE
       IF (.NOT. ASSOCIATED(OUT(ix)%ptr))    CYCLE
       levent = event_state(OUT(ix)%AVE_EVENT, current_date) .OR. lstart 
       IF (levent) THEN
          OUT(ix)%ptrout = OUT(ix)%ptr
       END IF
    END DO

    ! ----------------------------------------------------------------------

  END SUBROUTINE aveout_global_end
  ! ====================================================================


  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE aveout_read_nml_cpl(status, iou)
   
    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    
    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ OUTPUT

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='aveout_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! ### ADD HERE DIAGNOSTIC OUTPUT FOR LOG-FILE
    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE aveout_read_nml_cpl
  ! ====================================================================

! **********************************************************************
END MODULE messy_aveout_si
! **********************************************************************
