! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL COSMOTOY 
!
! Author : Mariano Mertens, DLR IPA, 2019 (largely copied from BUFLY by Patrick Joeckel)
!
! ONLY IMPLEMETED FOR COSMO
! 
! References: see messy_cosmotoy.f90
!
! **********************************************************************

! **********************************************************************
MODULE messy_cosmotoy_si
! **********************************************************************

#ifdef COSMO


  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      error_bi, info_bi, warning_bi

  ! SMCL
  USE messy_cosmotoy

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE

  ! GLOBAL PARAMETERS
  INTEGER, PARAMETER :: I_DO_IT_NOT       = 0
  

  ! CPL-NAMELIST PARAMETERS
  INTEGER :: i_do_it = I_DO_IT_NOT

  ! POINTER TO 2 meter temperature
  REAL(dp), DIMENSION(:,:), POINTER :: t2m_ptr => NULL()
  
  !Pointer to modified 2 meter temperature
  REAL(dp), DIMENSION(:,:), POINTER :: t2mmod_ptr => NULL()

  ! PUBLIC SUBROUTINES (called from messy_main_control_e5.f90)
  ! NOTE: in case you activate further entry points, make sure to call them
  !       in messy_main_control_e5.f90
  PUBLIC :: cosmotoy_initialize    ! initialize submodel
  PUBLIC :: cosmotoy_init_memory   ! request memory
!!$  PUBLIC :: cosmotoy_new_tracer    ! define new tracers
  PUBLIC :: cosmotoy_init_coupling ! set pointers for coupling to BM and other SMs
!!$  PUBLIC :: cosmotoy_init_tracer   ! initialize tracers
!!$  PUBLIC :: cosmotoy_global_start  ! entry point in time loop (all vectors)
  PUBLIC :: cosmotoy_physc         ! entry point in time loop (current vector)
!!$  PUBLIC :: cosmotoy_free_memory   ! free allocated memory

  ! PRIVATE SUBROUTINES
  !PRIVATE :: cosmotoy_read_nml_cpl

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE cosmotoy_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,     ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cosmotoy_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

    ! READ CTRL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find free I/O unit
       CALL cosmotoy_read_nml_ctrl(status, iou)  ! read CTRL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CTRL namelist',substr)
    END IF
    ! BROADCAST CTRL namleist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(scalingfac, p_io)


    ! READ CPL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find next free I/O unit
       CALL cosmotoy_read_nml_cpl(status, iou)  ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
    END IF


    ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(i_do_it, p_io)

    ! ### PERFORM INITIAL SETUP (CALL RESPECTIVE SMCL ROUTINE(S)) HERE

    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE cosmotoy_initialize
  ! ====================================================================

!!$  ! ====================================================================
!!$  SUBROUTINE cosmotoy_new_tracer
!!$
!!$    ! ------------------------------------------------------------------
!!$    ! This subroutine is used to define new tracers. See
!!$    ! http://www.atmos-chem-phys.net/8/1677   (including supplement !)
!!$    ! for full documentation.
!!$    ! ------------------------------------------------------------------
!!$
!!$    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
!!$    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
!!$    USE messy_main_tracer_bi,     ONLY: tracer_halt
!!$    ! MESSy
!!$    USE messy_main_tracer,        ONLY: new_tracer, set_tracer &
!!$                                      , R_molarmass ! ,ON, OFF
!!$    USE messy_main_constants_mem, ONLY: MO
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! LOCAL
!!$    CHARACTER(LEN=*), PARAMETER :: substr = 'cosmotoy_new_tracer'
!!$    INTEGER                     :: status
!!$    INTEGER                     :: idt_myO3   ! tracer index (identifier)
!!$
!!$    CALL start_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output
!!$
!!$    ! ### define new tracers here
!!$    CALL new_tracer(status, GPTRSTR, 'myO3', modstr, idx=idt_myO3)
!!$    CALL tracer_halt(substr, status)   ! terminate if error
!!$    CALL set_tracer(status, GPTRSTR, idt_myO3 &
!!$         , R_molarmass, r=3.0_dp*MO)
!!$    CALL tracer_halt(substr, status)   ! terminate if error
!!$
!!$    CALL end_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output
!!$
!!$  END SUBROUTINE cosmotoy_new_tracer
!!$  ! ====================================================================

!!$  ! ====================================================================
  SUBROUTINE cosmotoy_init_memory

    ! ------------------------------------------------------------------
    ! This subroutine is used to request memory for the submodel.
    ! The preferable method is to use "channel objects".
    ! Allocate your own memory, only if absolutely required.
    ! ------------------------------------------------------------------

    ! BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL
    USE messy_main_channel,          ONLY: new_channel, new_channel_object, &
                                           new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cosmotoy_init_memory'
    INTEGER                     :: status

    CALL start_message_bi(modstr,'MEMORY ALLOCATION',substr)  ! log-output
    ! ### ALLOCATE OWN MEMORY HERE, BUT ONLY IF ABSOLUTELY REQUIRED!
    CALL end_message_bi(modstr,'MEMORY ALLOCATION',substr)  ! log-output

    ! CHANNEL AND CHANNEL OBJECTS
    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

    ! new channel
    CALL new_channel(status, modstr//'_gp', reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)

    ! object with attributes
    CALL new_channel_object(status, modstr//'_gp', 'T_2M_MOD', &
         p2=t2mmod_ptr)

    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_gp', 'T_2M_MOD', &
         'long_name', c='modified temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_gp', 'T_2M_MOD', &
         'units', c='K ')
!!$    CALL channel_halt(substr, status)
!!$    CALL info_bi('channel/object '//modstr//'_gp/'// &
!!$         'objname'//' was created')
!!$
!!$    ! ### ADD MORE CHANNEL OBJECTS HERE
!!$
    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output


  END SUBROUTINE cosmotoy_init_memory
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE cosmotoy_init_coupling

    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the 
    ! basemodel and to other submodels.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: get_channel_object
    !
!!$    USE messy_main_tracer_bi, ONLY: tracer_halt
!!$    USE messy_main_tracer,    ONLY: get_tracer

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cosmotoy_init_coupling'
    INTEGER                     :: status

    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output


    ! ### set pointers to channel objects here
    ! In the example we take the 2-meter temperature
    CALL get_channel_object(status, 'COSMO_ORI', 'T_2M', p2=t2m_ptr)
    CALL channel_halt(substr, status)




    ! ### set pointers to tracers here
!!$    CALL get_tracer(status, ...)
!!$    CALL tracer_halt(substr, status)    ! terminate on error

    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output

  END SUBROUTINE cosmotoy_init_coupling
  ! ====================================================================

!!$  ! ====================================================================
!!$  SUBROUTINE cosmotoy_init_tracer
!!$
!!$    ! ------------------------------------------------------------------
!!$    ! This subroutine is used to initialise tracers (via NCREGRID)
!!$    ! according to the REGRID-namelists in cosmotoy_t.nml.
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
!!$  END SUBROUTINE cosmotoy_init_tracer
!!$  ! ====================================================================

  ! ====================================================================
!!  SUBROUTINE cosmotoy_global_start


    ! nothin to do here

!!  END SUBROUTINE cosmotoy_global_start
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE cosmotoy_physc

    ! ------------------------------------------------------------------
    ! This subroutine is called within the time loop.
    ! It constitutes the main entry point for additional processes 
    ! or diagnostics.
    ! Here, only the current vector of the grid-point-fields is
    ! accessible.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_timer,           ONLY: lstart, time_step_len
    USE messy_main_grid_def_mem_bi, ONLY: kproma, jrow, nlev

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cosmotoy_physc'
    INTEGER                     :: status

    
    t2mmod_ptr(1:kproma,jrow)=0.0_dp

! only for COSMO t2m_ptr is defined
    IF (i_do_it .eq. 1) then
       !t2mmod is modified 2 meter temperature
       t2mmod_ptr(1:kproma,jrow)=t2m_ptr(1:kproma,jrow)*scalingfac
    ELSE
       t2mmod_ptr(1:kproma,jrow)=t2m_ptr(1:kproma,jrow)      
    END IF

  END SUBROUTINE cosmotoy_physc
  ! ====================================================================

!!$  ! ====================================================================
!!$  SUBROUTINE cosmotoy_free_memory
!!$
!!$    ! ------------------------------------------------------------------
!!$    ! This subroutine is used to deallocate the memory, which has
!!$    ! been "manually" allocated in cosmotoy_init_memory.
!!$    ! Note: channel object memory must not be deallocated! This is
!!$    !       performed centrally.
!!$    ! ------------------------------------------------------------------
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! LOCAL
!!$    CHARACTER(LEN=*), PARAMETER :: substr = 'cosmotoy_free_memory'
!!$    INTEGER                     :: status
!!$
!!$  END SUBROUTINE cosmotoy_free_memory
!!$  ! ====================================================================

  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE cosmotoy_read_nml_cpl(status, iou)
   
    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    
    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ i_do_it

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='cosmotoy_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! ### ADD HERE DIAGNOSTIC OUTPUT FOR LOG-FILE
    SELECT CASE(i_do_it)
       CASE(0)
          WRITE(*,*) 'No calculation of modified T_2M'
       CASE(1)
          WRITE(*,*) 'Calculation of modified T_2M'
       CASE DEFAULT
          WRITE(*,*) '*** ERROR: UNRECOGNIZED SELECTION OF i_do_it (must be'//&
               &' 0, or  1)!'
          RETURN  ! with status = 1 => error
    END SELECT

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE cosmotoy_read_nml_cpl
  ! ====================================================================

#endif

! **********************************************************************
END MODULE messy_cosmotoy_si
! **********************************************************************
