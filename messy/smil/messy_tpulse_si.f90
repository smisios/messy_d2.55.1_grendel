! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL TPULSE 
!
! Author : Patrick Joeckel, DLR-IPA, October  2009
!
! References: see messy_tpulse.f90
!
! **********************************************************************

! **********************************************************************
MODULE messy_tpulse_si
! **********************************************************************

  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      error_bi, info_bi, warning_bi
  USE messy_main_timer_event,   ONLY: io_time_event, TRIG_FIRST, time_event

  ! SMCL
  USE messy_tpulse

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE

  ! GLOBAL PARAMETERS
  INTEGER, DIMENSION(:), ALLOCATABLE :: idx
  REAL(DP), POINTER :: counterptr => NULL()
  !INTEGER :: inct
  TYPE(io_time_event), PUBLIC :: trigreset = &
                                  io_time_event (1,'months',TRIG_FIRST,0)
  TYPE(time_event),    PUBLIC :: ev_trigreset
  LOGICAL                     :: l_trigreset

  ! CPL-NAMELIST PARAMETERS
  INTEGER :: npulse = 0

  ! PUBLIC SUBROUTINES (called from messy_main_control_echam5.f90)
  ! NOTE: in case you activate further entry points, make sure to call them
  !       in messy_main_control_echam5.f90
  PUBLIC :: tpulse_initialize    ! initialize submodel
  PUBLIC :: tpulse_new_tracer    ! define new tracers
  PUBLIC :: tpulse_init_memory   ! request memory and channel objects
 !!$ PUBLIC :: tpulse_init_coupling ! set pointers for coupling to BM and other SMs
!!$  PUBLIC :: tpulse_init_tracer   ! initialize tracers
  PUBLIC :: tpulse_global_start  ! entry point in time loop (all vectors)
!!$  PUBLIC :: tpulse_physc         ! entry point in time loop (current vector)
  PUBLIC :: tpulse_free_memory   ! free allocated memory

  ! PRIVATE SUBROUTINES
  !PRIVATE :: tpulse_read_nml_cpl

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE tpulse_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,     ONLY: find_next_free_unit
    USE messy_main_timer_bi,  ONLY: p_bcast_event, timer_event_init

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tpulse_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit
    ! check event with present date
    CHARACTER(LEN=*), PARAMETER :: EV_TLEV_PRES = 'present'

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

!!$    ! READ CTRL namelist
!!$    IF (p_parallel_io) THEN                  ! read only on I/O-PE
!!$       iou = find_next_free_unit(100,200)    ! find free I/O unit
!!$       CALL tpulse_read_nml_ctrl(status, iou)  ! read CTRL-namelist
!!$       ! terminate if error
!!$       IF (status /= 0) CALL error_bi('Error in reading CTRL namelist',substr)
!!$    END IF
!!$    ! BROADCAST CTRL namleist entries from I/O-PE to ALL OTHER PEs
!!$    CALL p_bcast(, p_io)

    ! READ CPL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find next free I/O unit
       CALL tpulse_read_nml_cpl(status, iou)  ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
    END IF
    ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(npulse, p_io)
    CALL p_bcast_event(trigreset,p_io)

    ! ### PERFORM INITIAL SETUP (CALL RESPECTIVE SMCL ROUTINE(S)) HERE
    ! initialize reset event
    CALL timer_event_init(ev_trigreset, trigreset, &
         'reset tracer', EV_TLEV_PRES)

    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE tpulse_initialize
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE tpulse_new_tracer

    ! ------------------------------------------------------------------
    ! This subroutine is used to define new tracers. See
    ! http://www.atmos-chem-phys.net/8/1677   (including supplement !)
    ! for full documentation.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: new_tracer, set_tracer &
                                      , R_molarmass ! ,ON, OFF
    USE messy_main_tools,         ONLY: int2str

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tpulse_new_tracer'
    INTEGER                     :: status
    INTEGER                     :: i
    CHARACTER(LEN=2)            :: istr

    CALL start_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output

    ALLOCATE(idx(npulse))

    ! ### define new tracers here
    DO i=1, npulse
       CALL int2str(istr, i)
       CALL new_tracer(status, GPTRSTR, 'PULSE'//istr, modstr, idx=idx(i))
!tracer manual describes option for subname
       CALL tracer_halt(substr, status)   ! terminate if error
       CALL set_tracer(status, GPTRSTR, idx(i) &
            , R_molarmass, r=146.07_dp) ! as SF6
       CALL tracer_halt(substr, status)   ! terminate if error
    END DO

    CALL end_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output

  END SUBROUTINE tpulse_new_tracer
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE tpulse_init_memory

    ! ------------------------------------------------------------------
    ! This subroutine is used to request memory for the submodel.
    ! The preferable method is to use "channel objects".
    ! Allocate your own memory, only if absolutely required.
    ! Info on CHANNELs and CHANNEL OBJECTS:
    ! http://www.geosci-model-dev.net/3/717/2010/
    ! supplement!!!
    ! ------------------------------------------------------------------

    ! BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: SCALAR
    USE messy_main_channel,          ONLY: new_channel, new_channel_object, &
                                            new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tpulse_init_memory'
    INTEGER                     :: status

    CALL start_message_bi(modstr,'MEMORY ALLOCATION',substr)  ! log-output
    ! ### ALLOCATE OWN MEMORY HERE, BUT ONLY IF ABSOLUTELY REQUIRED!
    CALL end_message_bi(modstr,'MEMORY ALLOCATION',substr)  ! log-output

    ! CHANNEL AND CHANNEL OBJECTS
    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

    ! new channel
    CALL new_channel(status, modstr, reprid=SCALAR)
    CALL channel_halt(substr, status)

    ! object with attributes
    CALL new_channel_object(status, modstr, 'counter', &
         p0=counterptr, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'counter', &
         'long_name', c='trigger counter')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'counter', &
         'units', c=' ')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object '//modstr// &
         'counter'//' was created')

    ! ### ADD MORE CHANNEL OBJECTS HERE

    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output


  END SUBROUTINE tpulse_init_memory
  ! ====================================================================

!!$  ! ====================================================================
!!$  SUBROUTINE tpulse_init_coupling
!!$
!!$    ! ------------------------------------------------------------------
!!$    ! This soubroutine is used to set pointers
!!$    ! (channel objects and/or tracers) for coupling to the 
!!$    ! basemodel and to other submodels.
!!$    ! ------------------------------------------------------------------
!!$
!!$    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
!!$    USE messy_main_channel_bi, ONLY: channel_halt
!!$    USE messy_main_channel,    ONLY: get_channel_object
!!$    !
!!$!!!$    USE messy_main_tracer_bi, ONLY: tracer_halt
!!$!!!$    USE messy_main_tracer,    ONLY: get_tracer
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! LOCAL
!!$    CHARACTER(LEN=*), PARAMETER :: substr = 'tpulse_init_coupling'
!!$    INTEGER                     :: status
!!$#if defined(E5202A) || defined(E5301)
!!$    CHARACTER(LEN=*), PARAMETER :: C_TTE    = 'tte_scb'
!!$    CHARACTER(LEN=*), PARAMETER :: C_ALPS   = 'alps_scb'
!!$    CHARACTER(LEN=*), PARAMETER :: C_ALPSTE = 'alpste_scb'
!!$#endif
!!$#if defined(E5302) 
!!$    CHARACTER(LEN=*), PARAMETER :: C_TTE    = 'tte'
!!$    CHARACTER(LEN=*), PARAMETER :: C_ALPS   = 'alps'
!!$    CHARACTER(LEN=*), PARAMETER :: C_ALPSTE = 'alpste'
!!$#endif
!!$
!!$    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output
!!$
!!$    ! ### set pointers to channel objects here
!!$
!!$    CALL get_channel_object(status, 'scnbuf', C_TTE, p3=tte_ptr)
!!$    CALL channel_halt(substr, status)
!!$
!!$
!!$    CALL get_channel_object(status, 'scnbuf', C_ALPS, p2=alps_ptr)
!!$    CALL channel_halt(substr, status)
!!$
!!$
!!$    CALL get_channel_object(status, 'scnbuf', C_ALPSTE, p2=alpste_ptr)
!!$    CALL channel_halt(substr, status)
!!$
!!$    ! ### set pointers to tracers here
!!$!!$    CALL get_tracer(status, ...)
!!$!!$    CALL tracer_halt(substr, status)    ! terminate on error
!!$
!!$    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output
!!$
!!$  END SUBROUTINE tpulse_init_coupling
!!$  ! ====================================================================

!!$  ! ====================================================================
!!$  SUBROUTINE tpulse_init_tracer
!!$
!!$    ! ------------------------------------------------------------------
!!$    ! This subroutine is used to initialise tracers (via NCREGRID)
!!$    ! according to the REGRID-namelists in tpulse_t.nml.
!!$    ! For a full documention of NCREGRID see
!!$    ! http://www.atmos-chem-phys.net/6/3557  (including supplement !)
!!$    ! ------------------------------------------------------------------
!!$
!!$    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
!!$    USE messy_main_tracer_bi, ONLY: tracer_init
!!$
!!$    IMPLICIT NONE
!!$    
!!$    CALL tracer_init(modstr) ! initialise tracers via NCREGRID
!!$
!!$  END SUBROUTINE tpulse_init_tracer
!!$  ! ====================================================================

  ! ====================================================================
  SUBROUTINE tpulse_global_start

    USE messy_main_timer,           ONLY: lstart, time_step_len
    USE messy_main_timer_bi,        ONLY: event_state
    USE messy_main_timer,           ONLY: current_date
    USE messy_main_tracer_mem_bi,   ONLY: xt, xtte, xtm1, xtf

    IMPLICIT NONE
    INTRINSIC :: NINT, MODULO, MOD

    CHARACTER(LEN=*), PARAMETER :: substr = 'tpulse_global_start'
    INTEGER :: status
    INTEGER, TARGET :: icnt
    INTEGER ::idt, iselect
    
    !icnt=0
    !iselect=0

    l_trigreset = event_state(ev_trigreset, current_date) .OR. lstart

    IF (l_trigreset) THEN 
       icnt = NINT(counterptr)!actually, NINT
       icnt = icnt + 1
!qqq
       iselect=MOD(icnt , npulse)+1 !!

       idt = idx(iselect)

       xt(:,:,idt,:) = 0.0_dp
       xtte(:,:,idt,:) = 0.0_dp
       xtm1(:,:,idt,:) = 0.0_dp
       xtf(:,:,idt,:) = 0.0_dp
       
       !icnt = icnt + 1
       counterptr = REAL(icnt, dp)
    ENDIF

    ! ----------------------------------------------------------------------

  END SUBROUTINE tpulse_global_start
  ! ====================================================================

!!$  ! ====================================================================
!!$  SUBROUTINE tpulse_physc
!!$
!!$    ! ------------------------------------------------------------------
!!$    ! This subroutine is called within the time loop.
!!$    ! It constitutes the main entry point for additional processes 
!!$    ! or diagnostics.
!!$    ! Here, only the current vector of the grid-point-fields is
!!$    ! accessible.
!!$    ! ------------------------------------------------------------------
!!$
!!$    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
!!$    USE messy_main_timer,         ONLY: lstart, time_step_len
!!$    USE messy_main_data_bi,       ONLY: kproma, jrow, nlev &
!!$                                      , philon_2d, philat_2d
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! LOCAL
!!$    CHARACTER(LEN=*), PARAMETER :: substr = 'tpulse_physc'
!!$    INTEGER                     :: status
!!$    REAL(DP), DIMENSION(:), POINTER  :: my_tte => NULL() ! perturbation tendency
!!$
!!$    IF (i_do_it /= I_DO_IT_IN_PHYSC) RETURN
!!$    IF (.NOT. lstart) RETURN
!!$
!!$    ! PART I: temperature perturbation
!!$
!!$    ! temporary memory
!!$    ALLOCATE(my_tte(kproma))   ! allocate with actual vector length
!!$                               ! (depends on jrow)
!!$
!!$    ! call SMCL routine to calculate perturbation tendency
!!$    CALL perturb(my_tte(:)           &     ! OUT: perturbation [K]
!!$         , philon_2d(1:kproma,jrow)  &     ! IN: geolocation
!!$         , philat_2d(1:kproma,jrow)  )     ! IN: geolocation
!!$    
!!$    ! add perturbation tendency to overall tendency
!!$    tte_ptr(1:kproma,nlev,jrow) = tte_ptr(1:kproma,nlev,jrow) &
!!$         + my_tte(:)/time_step_len  ! -> tendency [K/s]
!!$
!!$    DEALLOCATE(my_tte)
!!$    NULLIFY(my_tte)
!!$
!!$  END SUBROUTINE tpulse_physc
!!$  ! ====================================================================

  ! ====================================================================
  SUBROUTINE tpulse_free_memory

    ! ------------------------------------------------------------------
    ! This subroutine is used to deallocate the memory, which has
    ! been "manually" allocated in tpulse_init_memory.
    ! Note: channel object memory must not be deallocated! This is
    !       performed centrally.
    ! ------------------------------------------------------------------

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tpulse_free_memory'
    INTEGER                     :: status

    DEALLOCATE(idx)

  END SUBROUTINE tpulse_free_memory
  ! ====================================================================

  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE tpulse_read_nml_cpl(status, iou)
   
    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    
    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ npulse, trigreset

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='tpulse_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! ### ADD HERE DIAGNOSTIC OUTPUT FOR LOG-FILE
!qqq

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE tpulse_read_nml_cpl
  ! ====================================================================

! **********************************************************************
END MODULE messy_tpulse_si
! **********************************************************************
