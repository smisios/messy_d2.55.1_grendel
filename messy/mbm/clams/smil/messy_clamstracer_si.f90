!**********************************************************************
MODULE messy_clamstracer_si
!**********************************************************************
!  Submodel interface for clamstracer
!**********************************************************************

  USE messy_clamstracer
  USE messy_clamstracer_global, ONLY: tracer_type

  USE messy_clams_global,       ONLY: PREC

  USE messy_main_timer_event,   ONLY: time_event, io_time_event

  IMPLICIT NONE
  PRIVATE
  SAVE ! op_pj_20180613+
       ! variable 'io_tracerevent' at (1) with a component initialization must have the SAVE attribute 
       ! ... module variables should anyway have the SAVE attribute ...
       ! op_pj_20180613-

  INTRINSIC :: NULL

  TYPE(time_event),    PUBLIC :: tracerevent
  TYPE(io_time_event), PUBLIC :: io_tracerevent

 ! MODULE VARIABLES
  REAL(PREC),         DIMENSION(:), POINTER :: LAT   => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LON   => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LEV   => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: THETA_TROP1 => NULL()

  REAL(PREC),         DIMENSION(:), POINTER :: SGM   => NULL()

  type(tracer_type),  DIMENSION(:), POINTER :: TRACERS => NULL()

  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: clamstracer_initialize
  PUBLIC :: clamstracer_init_memory
  PUBLIC :: clamstracer_init_coupling
  PUBLIC :: clamstracer_global_start
  PUBLIC :: clamstracer_global_end
  PUBLIC :: clamstracer_free_memory
 

!--------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------

  SUBROUTINE clamstracer_initialize

    ! BMIL
    USE messy_main_blather_bi,     ONLY: error_bi
    USE messy_main_mpi_bi,         ONLY: p_pe
    USE messy_main_timer_bi,       ONLY: timer_event_init

    ! SMCL
    USE messy_main_tools,          ONLY: find_next_free_unit
    USE messy_main_timer,          ONLY: delta_time

    USE messy_clams_global,        ONLY: maxspec
    USE messy_clamstracer_global,  ONLY: timestep_tracer
    USE messy_clams_tools_utils,   ONLY: uppercase

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamstracer_initialize'

    INTEGER :: iou, status

    IF (p_pe==0) WRITE(*,*) uppercase(substr)


    allocate (TRACERS(maxspec))

    ! Read namelist
    iou = find_next_free_unit(100,200)
    CALL clamstracer_read_nml(status, iou)
    IF (status /= 0) CALL error_bi('Error in clamstracer_read_nml ',substr)

    if (mod(timestep_tracer*3600,int(delta_time)) /= 0) then
       call error_bi ("Wrong tracer timestep !!!",substr)
    endif

    ! Read list of tracers:
    iou = find_next_free_unit(100,200)
    CALL clamstracer_read_tracerlist(status, iou, TRACERS)
    IF (status /= 0) CALL error_bi('Error in clamstracer_read_tracerlist ',substr)

    ! Define TRACER event:
    io_tracerevent%counter = timestep_tracer   
    io_tracerevent%unit = 'hours'
    io_tracerevent%adjustment = 'exact'
    io_tracerevent%offset = -delta_time
    CALL timer_event_init (tracerevent, io_tracerevent, 'TRACER_Event', 'present')


  END SUBROUTINE clamstracer_initialize
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamstracer_init_memory

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: REPR_LG_CLAMS

    ! SMCL
    USE messy_main_channel,       ONLY: new_channel, new_attribute, &
                                        new_channel_object, &
                                        get_channel_object
    USE messy_clams_global,       ONLY: ldiagout, mdi, username
    USE messy_clamstracer_global, ONLY: ntracer, tracertype
    USE messy_clams_tools_utils,  ONLY: uppercase

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamstracer_init_memory'
    CHARACTER(40) :: cr_date
    CHARACTER(8)  :: ydate
    CHARACTER(10) :: ytime
    INTEGER       :: itracer, status

    IF (p_pe==0) WRITE(*,*) uppercase(substr)

    CALL DATE_AND_TIME(ydate, ytime)
    WRITE(cr_date,'(A4,5(A,A2))')   &
         ydate(1:4),'-',ydate(5:6),'-',ydate(7:8),' ',  &
         ytime(1:2),':',ytime(3:4),':',ytime(5:6) 

    !-----------------------------------------------------------------
    ! Define channel TRACER
    !-----------------------------------------------------------------

    CALL new_channel  (status, modstr, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)
    CALL new_attribute(status, modstr, modstr//'_version', c = modver)
    CALL channel_halt (substr, status)

    ! Define channel objects for all tracers

    DO itracer = 1, ntracer
       if (p_pe==0 .and. ldiagout) &
            write (*,*) 'create channel for ',trim(TRACERS(itracer)%name)
       CALL get_channel_object( status, 'clams', trim(TRACERS(itracer)%name), &
            p1=TRACERS(itracer)%values )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'creator_of_parameter', c = trim(username))
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'param_creation_time',  c = TRIM(cr_date))
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'param_modification_time', c = TRIM(cr_date))
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'long_name', c= trim(TRACERS(itracer)%name))
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'units', c= '')
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'missing_value', r=mdi )
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'flag', c='NONE')
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'pulse_date',c=trim(TRACERS(itracer)%pulse_date))
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'pulse_time',r=TRACERS(itracer)%pulse_time)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'pulse_length',i=TRACERS(itracer)%pulse_length)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'pulse_reset_period',i=TRACERS(itracer)%pulse_reset_period)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'pulse_time_chk',r=TRACERS(itracer)%pulse_time_chk)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'overwrite',i=TRACERS(itracer)%overwrite)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'latmin',r=TRACERS(itracer)%latmin)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'latmax',r=TRACERS(itracer)%latmax)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'lonmin',r=TRACERS(itracer)%lonmin)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'lonmax',r=TRACERS(itracer)%lonmax)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'levmin',r=TRACERS(itracer)%levmin)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'levmax',r=TRACERS(itracer)%levmax)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'seg_no',i=TRACERS(itracer)%seg_no)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'const',r=TRACERS(itracer)%const)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'lin',r=TRACERS(itracer)%lin)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'quad',r=TRACERS(itracer)%quad)
       CALL new_attribute(status, 'clams', trim(TRACERS(itracer)%name), &
            'tau',r=TRACERS(itracer)%tau)
       CALL channel_halt(substr, status)
    ENDDO

    
    CALL new_channel_object(status, modstr, "SGM", &
         p1=SGM, reprid=REPR_LG_CLAMS, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)


 END SUBROUTINE clamstracer_init_memory
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamstracer_init_coupling

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_channel_error_bi, ONLY: channel_halt

    ! SMCL
    USE messy_main_channel,        ONLY: get_channel_object, get_channel_object_info
    USE messy_clams_global,        ONLY: mdi
    USE messy_clams_tools_utils,   ONLY: uppercase

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamstracer_init_coupling'
    INTEGER :: status

    IF (p_pe==0) WRITE(*,*) uppercase(substr)

    ! get current airparcel positions:
    CALL get_channel_object(status, 'clams', 'LAT', p1=LAT)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON', p1=LON)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LEV', p1=LEV)
    CALL channel_halt(substr, status)

    ! get theta at tropopause
    CALL get_channel_object_info (status, 'clams', 'THETA_TROP1')
    if (status == 0) then
       CALL get_channel_object (status, 'clams', 'THETA_TROP1', p1=THETA_TROP1)
       CALL channel_halt (substr, status)
    else
       if (p_pe==0) WRITE(*,*) 'WARNING: THETA_TROP1 set to missing values !!!'
       allocate (THETA_TROP1(size(LAT)))
       THETA_TROP1 = mdi
    endif

  END SUBROUTINE clamstracer_init_coupling
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamstracer_global_start

    USE netcdf

    ! BMIL
    USE messy_main_mpi_bi,         ONLY: p_pe
    USE messy_main_blather_bi,     ONLY: error_bi

    !SMCL
    USE messy_main_timer,          ONLY: lstart, lresume
    USE messy_clams_global,        ONLY: ldiagout, initfile, buffersize, &
                                         dnparts, idx
    USE messy_clamstracer_global,  ONLY: ntracer, tracertype, &
                                         file_lsm, file_gph, lsm, gph
    USE messy_clamstracer_regional_masks, ONLY: read_era_data
    USE messy_clams_tools_utils,   ONLY: uppercase

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamstracer_global_start'
    INTEGER :: itracer
    INTEGER :: status, ncid, varid

    IF (p_pe==0 .and. ldiagout) WRITE(*,*) uppercase(substr)


    !*******************************************************
    ! Read INIT-File only in the first cycle:
    !*******************************************************
    if (lstart) then

       IF (p_pe==0 .and. ldiagout) write(*,*) uppercase(substr), ': in lstart'
 
       ! Open init file
       status = nf90_open (initfile, nf90_nowrite, ncid, buffersize)
       if (status /= nf90_noerr) call error_bi(substr,"Cannot open file "//trim(initfile))

       !--------------------------------------------------------------
       !  Read tracers (if present)
       !--------------------------------------------------------------
       do itracer = 1, ntracer

          status = nf90_inq_varid (ncid, TRACERS(itracer)%name, varid)

          if (status == 0) then

             status = nf90_get_var (ncid, varid, TRACERS(itracer)%values, &
                                    start=(/idx(p_pe,1)/), count=(/dnparts/))
             IF (status /= 0) CALL error_bi("Cannot read variable "//TRACERS(itracer)%name,substr)
          
             ! read 'pulse_time' attribute
             status = nf90_get_att (ncid, varid, 'pulse_time', TRACERS(itracer)%pulse_time_chk)

          else

             if (p_pe==0) write (*,*) 'Variable ',trim(TRACERS(itracer)%name), &
                                       ' not found, initialize with 0. !'
             TRACERS(itracer)%values = 0.

          endif

       enddo

       !--------------------------------------------------------------
       !  Read SGM 
       !--------------------------------------------------------------
       if (tracertype=='map') then  

          status = nf90_inq_varid (ncid, 'SGM', varid)
          IF (status /= 0) CALL error_bi("Cannot find variable SGM",substr)
          status = nf90_get_var (ncid, varid, SGM, start=(/idx(p_pe,1)/), count=(/dnparts/))
          IF (status /= 0) CALL error_bi("Cannot read variable SGM",substr)

       else

          SGM = 0.
          
       endif

       ! Close init file
       status = nf90_close (ncid)

    endif

    !***********************************************************
    ! Read land-sea-mask and orographie for tracertype "region"
    !***********************************************************
    if (lstart .or. lresume) then

       if (tracertype=='region') then
       
          call read_era_data (file_lsm, lsm, "LSM")
          call read_era_data (file_gph, gph, "Z")

       endif

    endif


  END SUBROUTINE clamstracer_global_start
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamstracer_global_end

    USE messy_main_mpi_bi,         ONLY: p_pe

    USE messy_clams_global,        ONLY: DP, ltracerevent, ldiagout
    USE messy_clams_tools_utils,   ONLY: uppercase

    IMPLICIT NONE


    CHARACTER(LEN=*), PARAMETER :: substr = 'clamstracer_global_end'
    INTEGER :: status

    IF (p_pe==0 .and. ldiagout) WRITE(*,*) uppercase(substr)

    IF (ltracerevent) THEN

       IF (p_pe==0) THEN
          WRITE (*,*)
          WRITE(*,*) 'Active tracerevent'
          WRITE (*,*)
       ENDIF
 
       !call set_artificial_tracers (status, LAT, LON, LEV, TRACERS, SGM)
       call set_artificial_tracers (status, LAT, LON, LEV, THETA_TROP1, TRACERS, SGM)

    ENDIF


  END SUBROUTINE clamstracer_global_end
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamstracer_free_memory

    USE messy_main_mpi_bi,         ONLY: p_pe
    USE messy_clams_tools_utils,   ONLY: uppercase

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamstracer_free_memory'

    IF (p_pe==0) WRITE(*,*) uppercase(substr)

    if (associated(TRACERS)) deallocate (TRACERS)


  END SUBROUTINE clamstracer_free_memory

!**********************************************************************
END MODULE messy_clamstracer_si
!**********************************************************************

