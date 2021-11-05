!**********************************************************************
MODULE messy_clamsdeepconv_si
!**********************************************************************
! Submodel interface for clamsdeepconv  
!**********************************************************************

  USE messy_clamsdeepconv
  USE messy_clamsdeepconv_global!, ONLY: 
  
  USE messy_clams_global,    ONLY: PREC

  USE messy_main_timer_event,   ONLY: time_event, io_time_event

  IMPLICIT NONE
  PRIVATE
  
  TYPE(time_event),    PUBLIC, SAVE :: deepconvevent
  TYPE(io_time_event), PUBLIC, SAVE :: io_deepconvevent

  ! MODULE VARIABLES
  REAL(PREC), DIMENSION(:), POINTER :: LAT        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: ZETA       => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: THETA      => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: TEMP       => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: PRESS      => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: BVF_WET    => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: H2O        => NULL()

  REAL(PREC), DIMENSION(:), POINTER :: DC_local   => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: DC         => NULL()
  

  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: clamsdeepconv_setup
  PUBLIC :: clamsdeepconv_initialize
  PUBLIC :: clamsdeepconv_init_memory
  PUBLIC :: clamsdeepconv_init_coupling
  PUBLIC :: clamsdeepconv_global_end
  PUBLIC :: clamsdeepconv_free_memory
  
!--------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE clamsdeepconv_setup

    ! BMIL
    USE messy_main_blather_bi,      ONLY: error_bi

    ! SMCL
    USE messy_clams_global,         ONLY: initfile, buffersize, init_vertcoorname, &
                                          rank, ldiagout
    USE messy_clamsdeepconv_global, ONLY: zeta_min
    USE messy_clams_tools_utils,    ONLY: uppercase

    USE netcdf
    
    IMPLICIT NONE
   
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsdeepconv_setup'
    integer :: status

    character(80) :: helpstr
    integer       :: ncid, varid
    real          :: rarr1(1)

    ! open init file
    status = nf90_open (initfile,nf90_nowrite,ncid, buffersize)
    IF (status /= 0) CALL error_bi('Error on open file '//trim(initfile),substr)

    ! read ZETA_DELTA(1) -> zeta_min
    helpstr = TRIM(uppercase(init_vertcoorname))//'_DELTA'
    status = nf90_inq_varid (ncid,TRIM(helpstr),varid)
    IF (status /= 0) CALL error_bi('Cannot find variable '//trim(helpstr),substr)
    status = nf90_get_var (ncid,varid,rarr1,start=(/1/),count=(/1/))
    IF (status /= 0) CALL error_bi('Cannot read variable '//trim(helpstr),substr)
    zeta_min = rarr1(1)
    if (rank==0 .and. ldiagout) write (*,*) 'clamsdeepconv_setup: zeta_min=',zeta_min
    
    ! close init file
    status = nf90_close(ncid)
    
  END SUBROUTINE clamsdeepconv_setup
  
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsdeepconv_initialize

    ! BMIL
    USE messy_main_blather_bi,     ONLY: error_bi
    USE messy_main_mpi_bi,         ONLY: p_pe
    USE messy_main_timer_bi,       ONLY: timer_event_init
    
    ! SMCL
    USE messy_main_tools,           ONLY: find_next_free_unit
    USE messy_main_timer,           ONLY: delta_time
    USE messy_clamsdeepconv_global, ONLY: timestep_deepconv
    USE messy_clams_global,         ONLY: met_freq
    USE messy_clams_tools_utils,    ONLY: uppercase
    
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsdeepconv_initialize'

    INTEGER :: iou, status

    IF (p_pe==0) WRITE(*,*) uppercase(substr)

    ! Read namelist
    iou = find_next_free_unit(100,200)
    CALL clamsdeepconv_read_nml(status, iou)
    IF (status /= 0) CALL error_bi('Error in clamsdeepconv_read_nml ',substr)

    if (mod(timestep_deepconv*3600,int(delta_time)) /= 0) then
       call error_bi ("Wrong deepconv timestep !!!",substr)
    endif

    ! Define DEEPCONV event:
    if (timestep_deepconv*3600 < delta_time) then
       if (p_pe==0) then
          write (*,*)
          write (*,*) 'timestep_deepconv < delta_time'
          write (*,*) 'timestep_deepconv:', timestep_deepconv
          write (*,*) 'delta_time:', delta_time
          write (*,*) 'timestep_deepconv is set to:', delta_time
          write (*,*)
       endif
       timestep_deepconv = NINT(delta_time)
    elseif (mod(timestep_deepconv*3600,int(delta_time)) /= 0) then
       call error_bi ("Wrong deep convection timestep: timestep_deepconv must be multiple of delta_time !!!",substr)
    elseif (timestep_deepconv > met_freq) then
       call error_bi ("Wrong deep convection timestep: "// &
            "timestep_deepconv must not be greater than the interval between met. files (met_freq) !!!", &
            substr)
    endif
    
    io_deepconvevent%counter = timestep_deepconv   
    io_deepconvevent%unit = 'hours'
    io_deepconvevent%adjustment = 'exact'
    io_deepconvevent%offset = 0.
    CALL timer_event_init (deepconvevent, io_deepconvevent, 'DEEPCONV_Event', 'present')
    
  END SUBROUTINE clamsdeepconv_initialize
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsdeepconv_init_memory

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_channel_bi,       ONLY: REPR_LG_CLAMS
    USE messy_main_channel_error_bi, ONLY: channel_halt

    ! SMCL
    USE messy_main_channel,       ONLY: new_channel,new_channel_object, &
                                        new_attribute
    USE messy_clams_tools_utils,  ONLY: uppercase
    USE messy_clams_global,       ONLY: mdi

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsdeepconv_init_memory'

    integer :: status
    
    IF (p_pe==0) WRITE(*,*) uppercase(substr)

    !-----------------------------------------------------------------
    ! Define channel CLAMSDEEPCONV
    !-----------------------------------------------------------------

    CALL new_channel  (status, modstr, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)
    CALL new_attribute(status, modstr, modstr//'_version', c = modver)
    CALL channel_halt (substr, status)

    !-----------------------------------------------------------------
    ! Define channel object DC_local
    !-----------------------------------------------------------------
    
    CALL new_channel_object(status, modstr, 'DC_local', p1=DC_local, reprid=REPR_LG_CLAMS)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'DC_local', 'missing_value', r = mdi)
    CALL channel_halt(substr, status)


  END SUBROUTINE clamsdeepconv_init_memory
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsdeepconv_init_coupling

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_blather_bi,       ONLY: error_bi

    ! SMCL
    USE messy_main_channel,        ONLY: get_channel_object
    USE messy_clams_global,        ONLY: ldiagout, nspec, specnames
    USE messy_clams_tools_utils,   ONLY: uppercase, str_found

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsdeepconv_init_coupling'
    INTEGER :: status

    IF (p_pe==0 .and. ldiagout) WRITE(*,*) uppercase(substr)

    ! Couple positions (LAT, LON, LEV/ZETA)
    CALL get_channel_object(status, 'clams', 'LAT', p1=LAT)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON', p1=LON)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LEV', p1=ZETA)
    CALL channel_halt(substr, status)

    ! Couple parameters (THETA, TEMP, PRESS, BVF_WET)
    CALL get_channel_object(status, 'clams', 'THETA', p1=THETA)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'TEMP', p1=TEMP)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'PRESS', p1=PRESS)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'BVF_WET', p1=BVF_WET)
    CALL channel_halt(substr, status)

  
    ! Couple H2O 
    ! => CHEM must be switched on: checked in clams_setup  
    if (.not. str_found (nspec, specnames, 'H2O')) then
       call error_bi ('H2O for DEEPCONV missing !!!', substr)
    endif
    CALL get_channel_object(status, 'clams', 'H2O', p1=H2O)
    CALL channel_halt(substr, status)     

    ! Couple DC 
    dc_mix = .false.
    if (str_found (nspec, specnames, 'DC')) then
       CALL get_channel_object(status, 'clams', 'DC', p1=DC)
       CALL channel_halt(substr, status)
       dc_mix = .true.
    endif
   
    
  END SUBROUTINE clamsdeepconv_init_coupling
    
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsdeepconv_global_end

    ! BMIL
    USE messy_main_mpi_bi,         ONLY: p_pe
    USE messy_main_blather_bi,     ONLY: error_bi

    ! SMCL
    USE messy_main_timer,          ONLY: lstart, lfirst_cycle
    USE messy_clams_global,        ONLY: ldeepconvevent, ldiagout, rank, mdi, &
                                         buffersize, initfile, idx, dnparts
    USE messy_clams_tools_utils,   ONLY: uppercase

    USE netcdf

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsdeepconv_global_end'

    REAL(PREC), ALLOCATABLE, DIMENSION(:) :: helparr
    INTEGER :: status, ncid, varid

    IF (p_pe==0 .and. ldiagout) WRITE(*,*) uppercase(substr)

    if (lstart) then

       ! Open init file
       status = nf90_open (initfile, nf90_nowrite, ncid, buffersize)
       if (status /= nf90_noerr) call error_bi("Cannot open file "//trim(initfile),substr)

       ! read DC
       allocate (helparr(dnparts))
       helparr = mdi
       DC_local = 0.
       status = nf90_inq_varid (ncid,'DC',varid)
       if (status /= nf90_noerr) CALL error_bi("Cannot find variable DC",substr)
       status = nf90_get_var (ncid,varid,helparr, start=(/idx(rank,1)/), count=(/dnparts/))
       if (status /= nf90_noerr) CALL error_bi("Cannot read variable DC",substr)
       DC_local(1:dnparts) = helparr
       if (dc_mix) DC = DC_local
       deallocate (helparr)
       
       ! close init file
       status = nf90_close (ncid)

    endif
    
    IF (ldeepconvevent .or. lfirst_cycle) THEN

       IF (p_pe==0) THEN
          WRITE (*,*)
          WRITE(*,*) 'Active deepconvevent'
          WRITE (*,*)
       ENDIF

       ! DC was possibly shuffled in MIX
       if (dc_mix) DC_local = DC
       
       call deep_conv (status, LAT, LON, ZETA, THETA, &
                       TEMP, PRESS, BVF_WET, H2O, DC_local)

       if (dc_mix) DC = DC_local
      
    ENDIF


  END SUBROUTINE clamsdeepconv_global_end
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsdeepconv_free_memory

    USE messy_main_mpi_bi,         ONLY: p_pe
    USE messy_clams_tools_utils,   ONLY: uppercase

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsdeepconv_free_memory'

    IF (p_pe==0) WRITE(*,*) uppercase(substr)


    
  END SUBROUTINE clamsdeepconv_free_memory

!**********************************************************************
END MODULE messy_clamsdeepconv_si
!**********************************************************************
