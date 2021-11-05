! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL BUFLY 
!
! Author : Patrick Joeckel, DLR-IPA, October  2009
!
! References: see messy_bufly.f90
!
! **********************************************************************

! **********************************************************************
MODULE messy_bufly_si
! **********************************************************************

  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      error_bi, info_bi, warning_bi

  ! SMCL
  USE messy_bufly

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE

  ! GLOBAL PARAMETERS
  INTEGER, PARAMETER :: I_DO_IT_NOT       = 0
  INTEGER, PARAMETER :: I_DO_IT_IN_PHYSC  = 1
  INTEGER, PARAMETER :: I_DO_IT_IN_GLOBAL = 2

  ! CPL-NAMELIST PARAMETERS
  INTEGER :: i_do_it = I_DO_IT_NOT

#ifdef ECHAM5
  ! POINTER TO TEMPERATURE TENDENCY
  REAL(dp), DIMENSION(:,:,:), POINTER :: tte_ptr => NULL()
  ! POINTERS TO LOG-SURFACE PRESSURE AND LOG-SURFACE PRESSURE TENDENCY
  REAL(dp), DIMENSION(:,:),   POINTER :: alps_ptr   => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: alpste_ptr => NULL()
#endif

  ! PUBLIC SUBROUTINES (called from messy_main_control_e5.f90)
  ! NOTE: in case you activate further entry points, make sure to call them
  !       in messy_main_control_e5.f90
  PUBLIC :: bufly_initialize    ! initialize submodel
!!$  PUBLIC :: bufly_init_memory   ! request memory
!!$  PUBLIC :: bufly_new_tracer    ! define new tracers
  PUBLIC :: bufly_init_coupling ! set pointers for coupling to BM and other SMs
!!$  PUBLIC :: bufly_init_tracer   ! initialize tracers
  PUBLIC :: bufly_global_start  ! entry point in time loop (all vectors)
  PUBLIC :: bufly_physc         ! entry point in time loop (current vector)
!!$  PUBLIC :: bufly_free_memory   ! free allocated memory

  ! PRIVATE SUBROUTINES
  !PRIVATE :: bufly_read_nml_cpl

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE bufly_initialize

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
    CHARACTER(LEN=*), PARAMETER :: substr = 'bufly_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

    ! READ CTRL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find free I/O unit
       CALL bufly_read_nml_ctrl(status, iou)  ! read CTRL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CTRL namelist',substr)
    END IF
    ! BROADCAST CTRL namleist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(r_lat, p_io)
    CALL p_bcast(r_lon, p_io)
    CALL p_bcast(r_dt,  p_io)
    CALL p_bcast(r_dps, p_io)

    ! READ CPL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find next free I/O unit
       CALL bufly_read_nml_cpl(status, iou)  ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
    END IF
    ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(i_do_it, p_io)

    ! ### PERFORM INITIAL SETUP (CALL RESPECTIVE SMCL ROUTINE(S)) HERE

    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE bufly_initialize
  ! ====================================================================

!!$  ! ====================================================================
!!$  SUBROUTINE bufly_new_tracer
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
!!$    CHARACTER(LEN=*), PARAMETER :: substr = 'bufly_new_tracer'
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
!!$  END SUBROUTINE bufly_new_tracer
!!$  ! ====================================================================

!!$  ! ====================================================================
!!$  SUBROUTINE bufly_init_memory
!!$
!!$    ! ------------------------------------------------------------------
!!$    ! This subroutine is used to request memory for the submodel.
!!$    ! The preferable method is to use "channel objects".
!!$    ! Allocate your own memory, only if absolutely required.
!!$    ! ------------------------------------------------------------------
!!$
!!$    ! BMIL
!!$    USE messy_main_channel_bi,   ONLY: channel_halt, GP_3D_MID
!!$    USE messy_main_channel,      ONLY: new_channel, new_channel_object, &
!!$                                       new_attribute
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! LOCAL
!!$    CHARACTER(LEN=*), PARAMETER :: substr = 'bufly_init_memory'
!!$    INTEGER                     :: status
!!$
!!$    CALL start_message_bi(modstr,'MEMORY ALLOCATION',substr)  ! log-output
!!$    ! ### ALLOCATE OWN MEMORY HERE, BUT ONLY IF ABSOLUTELY REQUIRED!
!!$    CALL end_message_bi(modstr,'MEMORY ALLOCATION',substr)  ! log-output
!!$
!!$    ! CHANNEL AND CHANNEL OBJECTS
!!$    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output
!!$
!!$    ! new channel
!!$    CALL new_channel(status, modstr//'_gp', reprid=GP_3D_MID)
!!$    CALL channel_halt(substr, status)
!!$
!!$    ! object with attributes
!!$    CALL new_channel_object(status, modstr//'_gp', 'objname', &
!!$         p3=objptr)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr//'_gp', 'objname', &
!!$         'long_name', c='explanation')
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr//'_gp', 'objname', &
!!$         'units', c=' ')
!!$    CALL channel_halt(substr, status)
!!$    CALL info_bi('channel/object '//modstr//'_gp/'// &
!!$         'objname'//' was created')
!!$
!!$    ! ### ADD MORE CHANNEL OBJECTS HERE
!!$
!!$    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output
!!$
!!$
!!$  END SUBROUTINE bufly_init_memory
!!$  ! ====================================================================

  ! ====================================================================
  SUBROUTINE bufly_init_coupling

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
    CHARACTER(LEN=*), PARAMETER :: substr = 'bufly_init_coupling'
    INTEGER                     :: status
#ifdef ECHAM5
    CHARACTER(LEN=*), PARAMETER :: C_TTE    = 'tte'
    CHARACTER(LEN=*), PARAMETER :: C_ALPS   = 'alps'
    CHARACTER(LEN=*), PARAMETER :: C_ALPSTE = 'alpste'
#endif

    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output

    ! ### set pointers to channel objects here

#ifdef ECHAM5
    CALL get_channel_object(status, 'scnbuf', C_TTE, p3=tte_ptr)
    CALL channel_halt(substr, status)


    CALL get_channel_object(status, 'scnbuf', C_ALPS, p2=alps_ptr)
    CALL channel_halt(substr, status)


    CALL get_channel_object(status, 'scnbuf', C_ALPSTE, p2=alpste_ptr)
    CALL channel_halt(substr, status)
#endif

    ! ### set pointers to tracers here
!!$    CALL get_tracer(status, ...)
!!$    CALL tracer_halt(substr, status)    ! terminate on error

    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output

  END SUBROUTINE bufly_init_coupling
  ! ====================================================================

!!$  ! ====================================================================
!!$  SUBROUTINE bufly_init_tracer
!!$
!!$    ! ------------------------------------------------------------------
!!$    ! This subroutine is used to initialise tracers (via NCREGRID)
!!$    ! according to the REGRID-namelists in bufly_t.nml.
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
!!$  END SUBROUTINE bufly_init_tracer
!!$  ! ====================================================================

#ifdef ECHAM5
  ! ====================================================================
  SUBROUTINE bufly_global_start

    USE messy_main_timer,           ONLY: lstart, time_step_len
    USE messy_main_grid_def_mem_bi, ONLY: nproma, npromz, ngpblks, nlev
    USE messy_main_grid_def_bi,     ONLY: philon_2d, philat_2d
    USE messy_main_transform_bi,    ONLY: locate_in_decomp
    USE messy_main_mpi_bi,          ONLY: p_pe

    IMPLICIT NONE
    INTRINSIC :: EXP, LOG 

    CHARACTER(LEN=*), PARAMETER :: substr = 'bufly_global_start'
    INTEGER :: status
    ! for temperature perturbation
    INTEGER                          :: jjrow     ! vector row loop index
    INTEGER                          :: zproma    ! vector length
    REAL(DP), DIMENSION(:), POINTER  :: my_t  => NULL() ! temp. perturbation
    ! for pressure perturbation
    INTEGER,  DIMENSION(4) :: pe      ! PROCESS ID
    INTEGER,  DIMENSION(4) :: jp      ! COLUMN
    INTEGER,  DIMENSION(4) :: jrow    ! ROW
    REAL(DP), DIMENSION(4) :: w       ! WEIGHT
    REAL(DP) :: lnaps1, lnaps2  ! ln(pressure [Pa])
    REAL(DP) :: aps             ! pressure [Pa]
    REAL(DP) :: dlnaps          ! delta of ln(pressure [Pa])

    IF (.NOT. lstart) RETURN

    ! ----------------------------------------------------------------------
    ! PART I: surface pressure perturbation (always in global_start)
    ! ----------------------------------------------------------------------

    CALL locate_in_decomp(status, r_lon(1), r_lat(1), pe, jp, jrow, w)
    IF (status /= 0) CALL error_bi('Error in locate_in_decomp.', substr)

    ! modify surface pressure tendency at lon,lat position
    IF (p_pe == pe(1)) THEN

       lnaps1 = alps_ptr(jp(1),jrow(1)) + &
            alpste_ptr(jp(1),jrow(1))*time_step_len

       aps = exp(lnaps1)

       aps = aps + r_dps

       lnaps2 = log(aps)

       dlnaps = (lnaps2 - lnaps1)

       alpste_ptr(jp(1),jrow(1)) = alpste_ptr(jp(1),jrow(1)) + &
            dlnaps/time_step_len

       WRITE(*,*) substr,': dlnaps=', dlnaps, &
            ' at (jp,jrow) = (', jp(1),',', jrow(1),')'

    END IF

    IF (p_pe == pe(2)) THEN

       lnaps1 = alps_ptr(jp(2),jrow(2)) + &
            alpste_ptr(jp(2),jrow(2))*time_step_len

       aps = exp(lnaps1)

       aps = aps - r_dps

       lnaps2 = log(aps)

       dlnaps = (lnaps2 - lnaps1)

       alpste_ptr(jp(2),jrow(2)) = alpste_ptr(jp(2),jrow(2)) + &
            dlnaps/time_step_len

       WRITE(*,*) substr,': dlnaps=', dlnaps, &
            ' at (jp,jrow) = (', jp(2),',', jrow(2),')'

    END IF

    ! ----------------------------------------------------------------------
    ! PART II: temperature perturbation (either in global_start or in physc)
    ! ----------------------------------------------------------------------
    IF (i_do_it /= I_DO_IT_IN_GLOBAL) RETURN

    ! temporary memory
    ALLOCATE(my_t(nproma))   ! allocate in full vector length

    ! loop over all vector rows
    DO jjrow=1, ngpblks

       ! set vector length
       IF (jjrow == ngpblks) THEN
          zproma = npromz   ! potentially shorter vector in last row
       ELSE
          zproma = nproma   ! 'normal' vector length
       ENDIF

       ! call SMCL routine to calculate perturbation tendency
       CALL perturb(my_t(1:zproma)       &     ! OUT: perturbation [K]
            , philon_2d(1:zproma,jjrow)  &     ! IN: geolocation
            , philat_2d(1:zproma,jjrow)  )     ! IN: geolocation

       ! add perturbation tendency to overall tendency
       tte_ptr(1:zproma,nlev,jjrow) = tte_ptr(1:zproma,nlev,jjrow) + &
            my_t(1:zproma)/time_step_len  ! -> tendency [K/s]
    END DO

    ! free temporary memory
    DEALLOCATE(my_t)
    NULLIFY(my_t)

    ! ----------------------------------------------------------------------

  END SUBROUTINE bufly_global_start
#endif
  ! ====================================================================

#ifdef ICON
  ! ====================================================================
  SUBROUTINE bufly_global_start

    USE messy_main_timer,           ONLY: lstart, time_step_len
    USE messy_main_grid_def_mem_bi, ONLY: nproma, npromz, ngpblks, nlev
    USE messy_main_transform_bi,    ONLY: locate_in_decomp
    USE messy_main_mpi_bi,          ONLY: p_pe
    USE messy_main_channel_mem,     ONLY: dom_current
    USE messy_main_constants_mem,   ONLY: RTD_icon
    ! ICON
    USE messy_main_bmluse_bi,       ONLY: p_patch

    IMPLICIT NONE
    INTRINSIC :: EXP, LOG 

    CHARACTER(LEN=*), PARAMETER :: substr = 'bufly_global_start'
    INTEGER :: status
    INTEGER :: ii
    INTEGER,  DIMENSION(4) :: pe      ! PROCESS ID
    INTEGER,  DIMENSION(4) :: jp      ! COLUMN
    INTEGER,  DIMENSION(4) :: jrow    ! ROW
    REAL(DP), DIMENSION(4) :: w       ! WEIGHT

    IF (.NOT. lstart) RETURN

    ! ----------------------------------------------------------------------
    ! PART I: surface pressure perturbation (always in global_start)
    ! ----------------------------------------------------------------------

    CALL locate_in_decomp(status, r_lon(1), r_lat(1), pe, jp, jrow, w)
    IF (status /= 0) CALL error_bi('Error in locate_in_decomp.', substr)

    IF (p_pe == pe(1)) THEN
       WRITE(*,'(A47,1X,I6)') "MESSy DEBUG: locate_in_decomp   dom_current =", dom_current
       WRITE(*,'(A47,2(1X,F6.2))') "MESSy DEBUG: locate_in_decomp   r_lon(1)/rlat(1) =" &
          & , r_lon(1), r_lat(1)
       WRITE(*,'(A36,4(1X,I6))') "MESSy DEBUG: locate_in_decomp   pe =", (pe(ii), ii=1,4)
       WRITE(*,'(A36,4(1X,I6))') "MESSy DEBUG: locate_in_decomp   jp =", (jp(ii), ii=1,4)
       WRITE(*,'(A36,4(1X,I6))') "MESSy DEBUG: locate_in_decomp jrow =", (jrow(ii), ii=1,4)
       WRITE(*,'(A36,4(1X,F6.4))') "MESSy DEBUG: locate_in_decomp    w =", (w(ii), ii=1,4)
       WRITE(*,'(A36,4(1X,F6.4))') "MESSy DEBUG: locate_in_decomp  lon =" &
          & , (p_patch(dom_current)%cells%center(jp(ii),jrow(ii))%lon, ii=1,4)
       WRITE(*,'(A36,4(1X,F6.4))') "MESSy DEBUG: locate_in_decomp  lat =" &
          & , (p_patch(dom_current)%cells%center(jp(ii),jrow(ii))%lat, ii=1,4)
       WRITE(*,'(A36,4(1X,F6.2))') "MESSy DEBUG: locate_in_decomp  lon =" &
          & , (p_patch(dom_current)%cells%center(jp(ii),jrow(ii))%lon*RTD_icon, ii=1,4)
       WRITE(*,'(A36,4(1X,F6.2))') "MESSy DEBUG: locate_in_decomp  lat =" &
          & , (p_patch(dom_current)%cells%center(jp(ii),jrow(ii))%lat*RTD_icon, ii=1,4)
    ENDIF

    ! ----------------------------------------------------------------------

  END SUBROUTINE bufly_global_start
  ! ====================================================================
#endif

#if defined(CESM1) || defined(COSMO) || defined(MESSYDWARF)
  ! ====================================================================
  SUBROUTINE bufly_global_start

  END SUBROUTINE bufly_global_start
  ! ====================================================================
#endif

  ! ====================================================================
  SUBROUTINE bufly_physc

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
#ifdef ECHAM5
    USE messy_main_grid_def_bi,     ONLY: philon_2d, philat_2d
#endif

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'bufly_physc'
    INTEGER                     :: status
#ifdef ECHAM5
    REAL(DP), DIMENSION(:), POINTER  :: my_tte => NULL() ! perturbation tendency

    IF (i_do_it /= I_DO_IT_IN_PHYSC) RETURN
    IF (.NOT. lstart) RETURN

    ! PART I: temperature perturbation

    ! temporary memory
    ALLOCATE(my_tte(kproma))   ! allocate with actual vector length
                               ! (depends on jrow)

    ! call SMCL routine to calculate perturbation tendency
    CALL perturb(my_tte(:)           &     ! OUT: perturbation [K]
         , philon_2d(1:kproma,jrow)  &     ! IN: geolocation
         , philat_2d(1:kproma,jrow)  )     ! IN: geolocation
    
    ! add perturbation tendency to overall tendency
    tte_ptr(1:kproma,nlev,jrow) = tte_ptr(1:kproma,nlev,jrow) &
         + my_tte(:)/time_step_len  ! -> tendency [K/s]

    DEALLOCATE(my_tte)
    NULLIFY(my_tte)
#endif

  END SUBROUTINE bufly_physc
  ! ====================================================================

!!$  ! ====================================================================
!!$  SUBROUTINE bufly_free_memory
!!$
!!$    ! ------------------------------------------------------------------
!!$    ! This subroutine is used to deallocate the memory, which has
!!$    ! been "manually" allocated in bufly_init_memory.
!!$    ! Note: channel object memory must not be deallocated! This is
!!$    !       performed centrally.
!!$    ! ------------------------------------------------------------------
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! LOCAL
!!$    CHARACTER(LEN=*), PARAMETER :: substr = 'bufly_free_memory'
!!$    INTEGER                     :: status
!!$
!!$  END SUBROUTINE bufly_free_memory
!!$  ! ====================================================================

  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE bufly_read_nml_cpl(status, iou)
   
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
    CHARACTER(LEN=*), PARAMETER :: substr='bufly_read_nml_cpl'
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
       CASE(I_DO_IT_NOT)
          WRITE(*,*) 'NO INITIAL TEMPERATURE PERTURBATION!'
       CASE(I_DO_IT_IN_PHYSC)
          WRITE(*,*) 'INITIAL TEMPERATURE PERTURBATION CALCULATED '//&
               &'IN BUFLY_PHYSC!'
       CASE(I_DO_IT_IN_GLOBAL)
          WRITE(*,*) 'INITIAL TEMPERATURE PERTURBATION CALCULATED '//&
               &'IN BUFLY_GLOBAL_START!'
       CASE DEFAULT
          WRITE(*,*) '*** ERROR: UNRECOGNIZED SELECTION OF i_do_it (must be'//&
               &' 0, 1, or 2)!'
          RETURN  ! with status = 1 => error
    END SELECT

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE bufly_read_nml_cpl
  ! ====================================================================

! **********************************************************************
END MODULE messy_bufly_si
! **********************************************************************
