!**********************************************************************
MODULE messy_clamssedi_si

#if defined(ECHAM5) || defined(MBM_CLAMS)
!**********************************************************************
!  Submodel interface for sedi -- calculation of sedimentation
! Authors:
! Thomas Breuer,     IEK-7, Mar 2015
! Nicole Thomas,     IEK-7, Jan 2015
!
! SUBROUTINE clamssedi_setup
! SUBROUTINE clamssedi_initialize
! SUBROUTINE clamssedi_init_memory
! SUBROUTINE clamssedi_init_coupling
! SUBROUTINE clamssedi_global_start
! SUBROUTINE clamssedi_global_end
! SUBROUTINE clamssedi_free_memory
! SUBROUTINE clamssedi_get_airparcels
! SUBROUTINE clamssedi_prepare_airparcels
! SUBROUTINE clamssedi_get_particles
! SUBROUTINE clamssedi_update_particles
! SUBROUTINE clamssedi_update_airparcels
! SUBROUTINE clamssedi_allocate
! SUBROUTINE clamssedi_deallocate
! SUBROUTINE clamssedi_error_handler (status, substr)
!     
!**********************************************************************

  USE messy_clamssedi
  USE messy_clamssedi_global

  USE messy_main_blather_bi,  ONLY: start_message_bi, end_message_bi

  USE messy_main_timer_event, ONLY: time_event, io_time_event!, event_is_active 

  IMPLICIT NONE
  PRIVATE

! op_pj_20160621+ compiler error workaround for gfortran 4.8, 4.7, ... 6.2, 7.1
!!$INTRINSIC :: NULL
#if defined(__GFORTRAN__) || defined(__G95__)
#define GF_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
! Test for GF
#if (GF_VERSION > 71300) && (GF_VERSION < 80000)
  INTRINSIC :: NULL
#endif
#else
#if !defined(LF) && !defined(NAGFOR)
  INTRINSIC :: NULL
#endif
#endif
! op_pj_20160621-

! op_pj_20160606+
!!$  TYPE(time_event),    PUBLIC :: sedievent
!!$  TYPE(io_time_event), PUBLIC :: io_sedievent
  TYPE(time_event),    PUBLIC, SAVE :: sedievent
  TYPE(io_time_event), PUBLIC, SAVE :: io_sedievent
! op_pj_20160606-

!!!!!  
!!$  logical,             PUBLIC :: lsedioutevent
! op_pj_20160606+
!!$  TYPE(time_event),    PUBLIC :: sedioutevent
!!$  TYPE(io_time_event), PUBLIC :: io_sedioutevent
!!!!!
!!$  TYPE(time_event),    PUBLIC, SAVE :: sedioutevent
!!$  TYPE(io_time_event), PUBLIC, SAVE :: io_sedioutevent
! op_pj_20160606-

  ! MODULE VARIABLES
  REAL(PREC),         DIMENSION(:), POINTER :: LAT         => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LON         => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LEV         => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: HNO3        => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: H2O         => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: TEMP        => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: PRESS       => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: NATbin      => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: ICEbin      => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: aer_H2SO4   => NULL()

  PUBLIC :: clamssedi_setup
  PUBLIC :: clamssedi_initialize
  PUBLIC :: clamssedi_init_memory
  PUBLIC :: clamssedi_init_coupling
  PUBLIC :: clamssedi_global_start  
  PUBLIC :: clamssedi_global_end
  PUBLIC :: clamssedi_free_memory

!--------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE clamssedi_setup

    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit

    implicit none 
    
    INTEGER :: iou, status
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamssedi_setup'

!!!!! nparticle_max muss hier bereits eingelesen werden !
!!!!! vor clamssedi_initialize werden bereits die Channels initialisiert !

    ! Read namelist variables:
    iou = find_next_free_unit(100,200)

    ! Read namelist and set default values:
    CALL clamssedi_read_nml(status, iou)
    IF (status /= 0) call error_bi  ('Error in subroutine clamssedi_read_nml!', substr)

  END SUBROUTINE clamssedi_setup
!--------------------------------------------------------------------
  SUBROUTINE clamssedi_initialize

    USE messy_main_mpi_bi,         ONLY: p_pe
    USE messy_main_timer_bi,       ONLY: timer_event_init
    USE messy_main_blather_bi,     ONLY: error_bi

    USE messy_main_timer,          ONLY: delta_time
    USE messy_clams_global,        ONLY: met_freq
    USE messy_clams_tools_utils,   ONLY: uppercase
    USE messy_clamssedi_data_io,   ONLY: get_init_grid, get_nairparcels_init, &
                                         read_nucleation_table

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamssedi_initialize'
    INTEGER :: status, iou

    IF (p_pe==0) WRITE(*,*) uppercase(substr)

    ! Define SEDI event:
!!$    io_sedievent%counter = timestep_sedi / 3600. 
!!$    io_sedievent%unit = 'hours'
!!$    io_sedievent%adjustment = 'exact'
!!$    io_sedievent%offset = 0.
!!$    CALL timer_event_init (sedievent, io_sedievent, 'SEDI_Event', 'present')
    if (timestep_sedi < delta_time) then
       if (p_pe==0) then
          write (*,*)
          write (*,*) 'timestep_sedi < delta_time'
          write (*,*) 'timestep_sedi:', timestep_sedi
          write (*,*) 'delta_time:', delta_time
          write (*,*) 'timestep_sedi is set to:', delta_time
          write (*,*)
       endif
       ! ub_ak_20180705+
       !timestep_sedi = delta_time
       timestep_sedi = NINT(delta_time)
       ! ub_ak_20180705-
    elseif (mod(timestep_sedi,int(delta_time)) /= 0) then
       call error_bi ("Wrong sedimentation timestep: timestep_sedi must be multiple of delta_time !!!",substr)
    elseif (timestep_sedi > met_freq*3600) then
       call error_bi ("Wrong sedimentation timestep: "// &
            "timestep_sedi must not be greater than the interval between met. files (met_freq) !!!", &
            substr)
    endif
    if (timestep_sedi > delta_time) then
       io_sedievent%counter = timestep_sedi    
       io_sedievent%unit = 'seconds'
       io_sedievent%adjustment = 'exact'
       io_sedievent%offset = 0.
       CALL timer_event_init (sedievent, io_sedievent, 'SEDI_Event', 'present')
    endif


    ! get grid information (nlevs, lev_grid, lev_delta, lev_window etc)
    call get_init_grid (status)
    IF (status /= 0) call error_bi  ('Error in subroutine get_init_grid !',substr)

    ! get number of airparcels per level in initfile
    call get_nairparcels_init (status)
    IF (status/=0) call error_bi ("Error in SEDI !!!", substr)

    ! read table for ice and nat nucleation
    call read_nucleation_table

!!!!!
!!$    ! set output event
!!$    if (loutput_sedi) then
!!$       if (mod(timestep_sediout*3600,int(delta_time)) /= 0) then
!!$          call error_bi ("SEDI output timestep must be multiple of delta_time !!!",substr)
!!$       endif
!!$       if (timestep_sediout*3600. /= timestep_sedi) then
!!$          io_sedioutevent%counter = timestep_sediout
!!$          io_sedioutevent%unit = 'hours'
!!$          io_sedioutevent%adjustment = 'exact'
!!$          io_sedioutevent%offset = -timestep_sedi
!!$          CALL timer_event_init (sedioutevent, io_sedioutevent, 'SEDIOUT_Event', 'present')
!!$       endif
!!$    endif     

  END SUBROUTINE clamssedi_initialize
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamssedi_init_memory
    
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: REPR_SEDI_PARTICLE, SCALAR
    USE messy_clamssedi,             ONLY: modstr
    USE messy_main_channel,          ONLY: new_channel, new_attribute, &
                                           new_channel_object
   ! USE messy_clamssedi_data_io,  ONLY: get_nairparcels_init

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamssedi_init_memory'
    INTEGER       :: status
    CHARACTER(30) :: username


    if (p_pe==0) write(*,*) substr

    username = 'j.-u.grooss'
    CALL new_channel(status, modstr, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, modstr//'_version', c = modver)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'nparticles', p0=nparticles_co, &
                            reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)
    CALL new_channel_object(status, modstr, 'part_id_max', p0=part_id_max_co, &
                            reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)

    CALL new_channel_object(status, modstr, 'lat', p1=particles%lat, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'lon', p1=particles%lon, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'lev', p1=particles%lev, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'pressure', p1=particles%pressure, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'temperature', p1=particles%temperature, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'radius', p1=particles%radius, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'density', p1=particles%density, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'sedimentation', p1=particles%sedimentation, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'tsv', p1=particles%tsv, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'hno3', p1=particles%hno3, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'h2o', p1=particles%h2o, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'sice', p1=particles%sice, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'snat', p1=particles%snat, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'icebin', p1=particles%icebin, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'natbin', p1=particles%natbin, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'airparcel_density_change', p1=particles%airparcel_density_change, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'class', p1=particles%class, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'particle_id', p1=particles%particle_id, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'snatmax', p1=particles%snatmax, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'sicemax', p1=particles%sicemax, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'tmin', p1=particles%tmin, reprid=REPR_SEDI_PARTICLE)
    CALL channel_halt(substr, status)

    
    allocate (ntriang_lev(nlevs))              ! number of triangles per level 
    
  END SUBROUTINE clamssedi_init_memory
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamssedi_init_coupling

    ! BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_mpi_bi,           ONLY: p_pe

    ! SMCL
    USE messy_main_channel,       ONLY: get_channel_object
    USE messy_clams_global,       ONLY: nspec, specnames
    USE messy_clams_tools_utils,  ONLY: uppercase, str_found

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamssedi_init_coupling'
    INTEGER  :: status

    IF (p_pe==0) WRITE(*,*) uppercase(substr)

    ! Get arrays from channel CLAMS (LAT, LON, LEV, TEMP, PRESS)
    CALL get_channel_object(status, 'clams', 'LAT', p1=LAT)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON', p1=LON)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LEV', p1=LEV)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'TEMP', p1=TEMP)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'PRESS', p1=PRESS)
    CALL channel_halt(substr, status)

    ! Check if the species, which will be coupled, exist
    if (.not. str_found (nspec, specnames, 'HONO2') .or. &
        .not. str_found (nspec, specnames, 'H2O') .or. &
        .not. str_found (nspec, specnames, 'NATbin') .or. &
        .not. str_found (nspec, specnames, 'ICEbin') .or. &
        .not. str_found (nspec, specnames, 'aer_H2SO4') ) then
       call error_bi &
            ('One or more of the species used for SEDI missing: '// &
             'HONO2, H2O, NATbin, ICEbin, aer_H2SO4 !!!', substr)
    endif

    ! Get species from channel clams (HNO3, H2O)
    CALL get_channel_object(status, 'clams', 'HONO2', p1=HNO3)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'H2O', p1=H2O)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'NATbin', p1=NATbin)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'ICEbin', p1=ICEbin)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'aer_H2SO4', p1=aer_H2SO4)
    CALL channel_halt(substr, status)
    

  END SUBROUTINE clamssedi_init_coupling
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamssedi_global_start

    USE messy_main_mpi_bi,         ONLY: p_pe
    USE messy_main_timer,          ONLY: lfirst_cycle, current_date, lresume
    USE messy_main_timer_bi,    ONLY: event_state
    
    USE messy_clams_global,        ONLY: ldiagout
    USE messy_clams_tools_utils,   ONLY: uppercase

    IMPLICIT NONE
    
    INTEGER :: status

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamssedi_global_start'

!!!!!
!!$    if (.not. loutput_sedi) then
!!$       lsedioutevent = .FALSE.
!!$    elseif (timestep_sediout*3600. == timestep_sedi) then
!!$       lsedioutevent = .TRUE.
!!$    else
!!$       lsedioutevent = event_state(sedioutevent, current_date)
!!$    endif
    
    ! restart: set nparticles
    if (lresume) then 
       ! ub_ak_20180705+
       !nparticles = nparticles_co
       !part_id_max = part_id_max_co
       nparticles = NINT(nparticles_co)
       part_id_max = NINT(part_id_max_co)
       ! ub_ak_20180705-
    endif


!!$    IF (p_pe==0 .and. ldiagout) write(*,*) uppercase(substr)
!!$
!!$    IF (lfirst_cycle) THEN
!!$
!!$       IF (p_pe==0) write(*,*) uppercase(substr), ': in lfirst_cycle'
!!$
!!$
!!$    ENDIF

  END SUBROUTINE clamssedi_global_start
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamssedi_global_end

    ! BMIL
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_mpi_bi,           ONLY: p_barrier, p_all_comm, p_sum
    USE messy_main_channel_error_bi, ONLY: channel_halt

    ! SMCL
    USE messy_main_channel,         ONLY: set_channel_output
    USE messy_main_timer,           ONLY: lfirst_cycle
    USE messy_clams_global,         ONLY: rank, lsedievent, pi, mdi, eps, lcoupled
    USE messy_clamssedi_create_pos, ONLY: create_position
    USE messy_clamssedi_triang,     ONLY: triangulation
    use messy_clamssedi_data_io,    only: read_clams_particles
    USE messy_main_timer,           ONLY: HOUR
    
    IMPLICIT NONE
    
    integer :: status = 0
    integer :: ierr

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamssedi_global_end'

    integer       :: startclock, endclock, counts_per_sec
    real(DP)      :: seconds

    CALL system_clock (count=startclock)


    if (rank==0) write (*,*) 'in clamssedi_global_start: lsedievent=', lsedievent

    if (lsedievent .or. lfirst_cycle) then

       IF (rank==0) THEN
          WRITE (*,*)
          WRITE(*,*) 'Active sedievent'
          WRITE (*,*)
       ENDIF
       
       ! if desired: read initial particles from file before the first time step 
       if (lfirst_cycle) then
          if (part_init_file == '') then 
             if (rank == 0) write(*,*) 'Don''t read initial particles from file!'
          else
             if (rank == 0) write(*,*) 'read initial particles from: ', part_init_file
             call read_clams_particles(status)
          endif
          IF (status /= 0) call error_bi  ('Error in subroutine read_clams_particles!',substr)
       endif

       ! copy airparcel positions (LAT,LON,LEV), HNO3 and H2O 
       ! from all ranks to global arrays
       write (*,*) 'call clamssedi_get_airparcels', rank
       call clamssedi_get_airparcels 
  
! !!$!!!!! TESTAUSGABE:
! !!$write (30+rank,*) airparcels%lat
! !!$CALL p_barrier(p_all_comm) 
! !!$call error_bi  ('stop nach clamssedi_get_airparcels',substr) 

       if (rank==0) write (*,*) 'call create_position'
       ! create new particles (%lat/%lon are in degree)
       call create_position (status)
       IF (status /= 0) call error_bi  ('Error in subroutine create_position!',substr)

       !!! each particles%lat and %lon have to be in degree!!!
       write (*,*) 'call clamssedi_get_particles', rank
       ! gather and broadcast all particles and convert lat and lon from degree to radian
       call clamssedi_get_particles
       
       write (*,*) 'call triangulation', rank
       ! create triangles
       call triangulation(status)
       IF (status /= 0) call error_bi  ('Error in subroutine get_triangles!',substr)
       
       write (*,*) 'call clamssedi_allocate',rank
       call clamssedi_allocate
       
       write (*,*) 'call sedi_prepare', rank
       call sedi_prepare (status)
       if (status /= 0) call error_bi  ('Error in subroutine sedi_prepare!',substr)
       
       ! all_reduce: nparticles_nat/ice of each triangles
       triangles%nparticles_ice = p_sum(triangles%nparticles_ice)
       triangles%nparticles_nat = p_sum(triangles%nparticles_nat)
       
       write (*,*) 'call clamssedi', rank
       call clamssedi (status, lcoupled)
       if (status /= 0) call error_bi  ('Error in subroutine clamssedi!',substr)


       ! update particles on all ranks 
       call clamssedi_update_particles


       ! sum up changes on airparcels (HNO3, H20) from all ranks and 
       ! write changes to channel clams
       call clamssedi_update_airparcels


       write (*,*) 'call clamssedi_deallocate', rank
       call clamssedi_deallocate

       IF (rank==0) THEN
          WRITE (*,*) 'Normal termination of sedi'
          
          CALL system_clock (COUNT_RATE=counts_per_sec)
          CALL system_clock (count=endclock)
          IF (endclock > startclock) THEN
             seconds = float(endclock-startclock) / float(counts_per_sec)
          ELSE
             seconds = 0.0
          ENDIF
          
          WRITE (*,*)
          WRITE (*,'(A,F10.2,A)') 'This job has taken ',seconds,' seconds to execute.'
       ENDIF
       
#ifdef MBM_CLAMS
!!$       IF (loutput_sedi .and. lsedioutevent .and. (nparticles /= 0 .or. HOUR == 11)) then
          nparticles_co = nparticles
          part_id_max_co = part_id_max
!!!!! => Output in BMIL clams_main.f90           
          ! CALL set_channel_output(status, 'clams', .FALSE.)
          ! CALL channel_halt(substr, status)
          !CALL set_channel_output(status, 'clamssedi', .TRUE.)
          !CALL channel_halt(substr, status)
          ! CALL messy_write_output
          ! if (lclamsoutevent) then
          !    CALL set_channel_output(status, 'clams', .TRUE.)
          !    CALL channel_halt(substr, status)
          ! endif
!!$       ENDIF
#endif

    endif ! lsedievent

  END SUBROUTINE clamssedi_global_end
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamssedi_free_memory

    USE messy_clams_global,     ONLY: rank

    IMPLICIT NONE
    
    deallocate (nairparcels50_lev_init)
    deallocate (ntriang_lev)

    deallocate(lev_grid)
    deallocate(lev_window)
    if (associated(pos_r_grid)) deallocate(pos_r_grid)
    deallocate(pos_lev_grid)
    deallocate(pos_lev_delta)

    deallocate(snat_table)
    deallocate(xnnat_table)
    deallocate(sice_table)
    deallocate(xnice_table)
    
  END SUBROUTINE clamssedi_free_memory
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamssedi_get_airparcels

    ! BMIL
    USE messy_main_mpi_bi,      ONLY: p_recv, p_send, p_bcast, p_sum

    ! SMCL
    USE messy_clams_global,     ONLY: ntasks, rank, nparts, dnparts, &
                                      pi, mdi, eps
    USE messy_clamssedi_global, ONLY: airparcels
    
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamssedi_get_airparcels'

    integer,    dimension(:), allocatable :: nparts_pe
 
    integer                               :: irank, pos, i
    integer :: startpos, iirank

    ! get number of airparcels on all tasks
    allocate (nparts_pe(0:ntasks-1))
    nparts_pe = 0
    nparts_pe(rank) = dnparts
    do irank = 0, ntasks-1
       call p_bcast (nparts_pe(irank),irank)
    enddo
    
    if (rank==0) write (*,*) 'nparts_pe:',nparts_pe, nparts

    ! check if sum of airparcels is correct
    if (p_sum(dnparts) /= nparts) call clamssedi_error_handler (101, substr)

    ! allocate global arrays for airparcels
    allocate (airparcels%lat(nparts))
    allocate (airparcels%lon(nparts))
    allocate (airparcels%lev(nparts))
    allocate (airparcels%hno3(nparts))
    allocate (airparcels%h2o(nparts))
    allocate (airparcels%diff_hno3(nparts))
    allocate (airparcels%diff_h2o(nparts))
    allocate (airparcels%natbin(nparts))
    allocate (airparcels%icebin(nparts))
    allocate (airparcels%natbin_diff(nparts))
    allocate (airparcels%icebin_diff(nparts))
    allocate (airparcels%coor(3,nparts))
    allocate (airparcels%ntriang(nparts))
    allocate (airparcels%triang_ids(max_nb,nparts))
    allocate (airparcels%ilev(nparts))
    allocate (airparcels%temp(nparts))
    allocate (airparcels%press(nparts))
    allocate (airparcels%aer_h2so4(nparts))

    ! initialize components of airparcels
    airparcels%diff_hno3  = 0.
    airparcels%diff_h2o   = 0.
    airparcels%natbin_diff= 0.
    airparcels%icebin_diff= 0.
    airparcels%ntriang    = 0
    airparcels%ilev       = 0 
    airparcels%triang_ids = 0

    ! copy airparcel positions (LAT,LON.LEV), HNO3 and H2O 
    ! from all ranks to global arrays on rank 0
    if (rank==0) then

       ! write own data to global fields
       airparcels%lat   (1:dnparts) = LAT   (1:dnparts)
       airparcels%lon   (1:dnparts) = LON   (1:dnparts)
       airparcels%lev   (1:dnparts) = LEV   (1:dnparts)
       airparcels%hno3  (1:dnparts) = HNO3  (1:dnparts)
       airparcels%h2o   (1:dnparts) = H2O   (1:dnparts)
       airparcels%temp  (1:dnparts) = TEMP  (1:dnparts)
       airparcels%press (1:dnparts) = PRESS (1:dnparts)
       airparcels%natbin(1:dnparts) = NATbin(1:dnparts)
       airparcels%icebin(1:dnparts) = ICEbin(1:dnparts)
       airparcels%aer_h2so4(1:dnparts) = aer_H2SO4(1:dnparts)

       pos = dnparts 

       ! get data from all other ranks and add to global field
       do irank = 1, ntasks-1

          !write (*,*) 'startpos,endpos:',pos+1,pos+nparts_pe(irank)
          call p_recv (airparcels%lat   (pos+1:pos+nparts_pe(irank)), irank, 100*irank+1)
          call p_recv (airparcels%lon   (pos+1:pos+nparts_pe(irank)), irank, 100*irank+2)
          call p_recv (airparcels%lev   (pos+1:pos+nparts_pe(irank)), irank, 100*irank+3)
          call p_recv (airparcels%hno3  (pos+1:pos+nparts_pe(irank)), irank, 100*irank+4)
          call p_recv (airparcels%h2o   (pos+1:pos+nparts_pe(irank)), irank, 100*irank+5)
          call p_recv (airparcels%temp  (pos+1:pos+nparts_pe(irank)), irank, 100*irank+6)
          call p_recv (airparcels%press (pos+1:pos+nparts_pe(irank)), irank, 100*irank+7)
          call p_recv (airparcels%natbin(pos+1:pos+nparts_pe(irank)), irank, 100*irank+8)
          call p_recv (airparcels%icebin(pos+1:pos+nparts_pe(irank)), irank, 100*irank+9)
          call p_recv (airparcels%aer_h2so4(pos+1:pos+nparts_pe(irank)), irank, 100*irank+10)

          pos = pos + nparts_pe(irank) 

       enddo

       ! convert lat and lon to radian
       where (abs((airparcels%lon-mdi)/mdi)>eps) 
          airparcels%lat = airparcels%lat * pi/180.
          airparcels%lon = airparcels%lon * pi/180.
       end where
       
    else

       ! send local field to rank 0
       call p_send (LAT   (1:dnparts),0,100*rank+1)
       call p_send (LON   (1:dnparts),0,100*rank+2)
       call p_send (LEV   (1:dnparts),0,100*rank+3)
       call p_send (HNO3  (1:dnparts),0,100*rank+4)
       call p_send (H2O   (1:dnparts),0,100*rank+5)
       call p_send (TEMP  (1:dnparts),0,100*rank+6)
       call p_send (PRESS (1:dnparts),0,100*rank+7)
       call p_send (NATbin(1:dnparts),0,100*rank+8)
       call p_send (ICEbin(1:dnparts),0,100*rank+9)
       call p_send (aer_H2SO4(1:dnparts),0,100*rank+10)

    endif
    
    ! broadcast global arrays to all ranks
    call p_bcast (airparcels%lat, 0)
    call p_bcast (airparcels%lon, 0)
    call p_bcast (airparcels%lev, 0)
    call p_bcast (airparcels%hno3, 0)
    call p_bcast (airparcels%h2o, 0)
    call p_bcast (airparcels%temp, 0)
    call p_bcast (airparcels%press, 0)
    call p_bcast (airparcels%natbin, 0)
    call p_bcast (airparcels%icebin, 0)
    call p_bcast (airparcels%aer_h2so4, 0)
  

    startpos = 1
    do iirank = 1, rank
       startpos = startpos + nparts_pe(iirank-1)
    enddo

    ! clean up
    deallocate (nparts_pe)

    ! set level index, convert coordinates, determine hno3 in gas phase
    call clamssedi_prepare_airparcels


    
    do i = 1, dnparts
       if (HNO3(i) < 0.) then
          write(*,*) "THOM1: ", rank, i, 'HNO3', LEV(i), HNO3(i), startpos, startpos+i-1, airparcels%ntriang(startpos+i-1), airparcels%lev(startpos+i-1), airparcels%hno3(startpos+i-1)
       endif
       if (H2O(i) < 0.) then
          write(*,*) "THOM1: ", rank, i, 'H2O', LEV(i), H2O(i), startpos, startpos+i-1, airparcels%ntriang(startpos+i-1), airparcels%lev(startpos+i-1), airparcels%h2o(startpos+i-1)
       endif
    enddo

    
  END SUBROUTINE clamssedi_get_airparcels
!--------------------------------------------------------------------

  subroutine clamssedi_prepare_airparcels

    use messy_clamssedi_global,     only: nhemi, airparcels, nlevs, &
                                          nairparcels_lev, nairparcels50_lev, &
                                          lev_window, max_ntriang_lev

    use messy_clams_global,         only: nparts, pi, prec
    use messy_clamssedi_hetero_shi, only: calc_parthno3
    
    implicit none
    
    integer :: ilev, ipart 
    integer :: max_nairparcels_lev

    real(kind=prec) :: parthno3
    
    ! counter for airparcels per level
    allocate(nairparcels_lev(nlevs))
    allocate(nairparcels50_lev(nlevs))
    nairparcels_lev = 0
    nairparcels50_lev = 0

    ! initialize level index
    airparcels%ilev(:) = -1

    do ipart = 1, nparts       ! airparcel loop
       do ilev = 1, nlevs      ! level loop
          if (lev_window(1,ilev) <= airparcels%lev(ipart) .and. &
               airparcels%lev(ipart) < lev_window(2,ilev)) then

             nairparcels_lev(ilev) = nairparcels_lev(ilev) + 1

             if ( nhemi .and. (airparcels%lat(ipart) > 50.*(pi/180.)) ) then
                nairparcels50_lev(ilev) = nairparcels50_lev(ilev) + 1
             elseif ( airparcels%lat(ipart) < -50.*(pi/180.) ) then
                nairparcels50_lev(ilev) = nairparcels50_lev(ilev) + 1
             endif

             !! set level index
             airparcels%ilev(ipart) = ilev

             !! Transfer from spheric coordinates to cartesian coordinates on a unit sphere (r=1)
             airparcels%coor(1,ipart) = cos(airparcels%lon(ipart)) * sin(pi/2. - airparcels%lat(ipart))
             airparcels%coor(2,ipart) = sin(airparcels%lon(ipart)) * sin(pi/2. - airparcels%lat(ipart))
             airparcels%coor(3,ipart) = cos(pi/2. - airparcels%lat(ipart))

             ! determine hno3 in gas phase
             call calc_parthno3(airparcels%press(ipart), airparcels%temp(ipart), &
                                airparcels%hno3(ipart), airparcels%h2o(ipart), &
                                airparcels%aer_h2so4(ipart), parthno3)
             airparcels%hno3(ipart) = airparcels%hno3(ipart) * parthno3

             EXIT     ! process next airparcel
          endif
       enddo
    enddo
    
    max_nairparcels_lev = MAXVAL(nairparcels_lev)
    max_ntriang_lev = max_nairparcels_lev * 3

    ! only used for subroutine calc_parthno3
    deallocate(airparcels%temp)
    deallocate(airparcels%press)
    deallocate(airparcels%aer_h2so4)
    
  end subroutine clamssedi_prepare_airparcels

!--------------------------------------------------------------------
  SUBROUTINE clamssedi_get_particles

    ! BMIL
    USE messy_main_mpi_bi,      ONLY: p_recv, p_send, p_bcast, p_sum

    ! SMCL
    USE messy_clams_global,     ONLY: rank, ntasks, mdi, eps, pi
    USE messy_clamssedi_global, ONLY: particles, &
                                      nparticles, nparticles_added, nparticle_max, &
                                      hno3_background, h2o_background, &
                                      densnat_val, radius_default

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamssedi_get_particles'

    integer,    dimension(:), allocatable :: nparts_pe
    integer                               :: nparticles_old, irank, pos

    nparticles_old = nparticles
    if (rank==0) write (*,*) 'nparticles_old=',nparticles_old
    
    ! get number of particles on all tasks
    allocate (nparts_pe(0:ntasks-1))
    nparts_pe = 0
    nparts_pe(rank) = nparticles_added
    do irank = 0, ntasks-1
       call p_bcast (nparts_pe(irank),irank)
       if (rank==0) write (*,*) 'irank, nparticles_added(irank):',irank, nparts_pe(irank)
    enddo

    ! check sum of parcticles  
!!!!! Alternative zu p_sum:: addiere alle Elemente von nparts_pe
    if (p_sum(nparticles_added)+nparticles > nparticle_max) call clamssedi_error_handler (201, substr)

    ! copy particle positions (LAT,LON.LEV), temperature and pressure 
    ! from all ranks to global arrays on rank 0
    if (rank==0) then

              
       ! added particles on rank 0:
       nparticles = nparticles + nparticles_added

       ! get data from all other ranks and add to global field
       do irank = 1, ntasks-1

          call p_recv (particles%lat (nparticles+1:nparticles+nparts_pe(irank)), irank, 100*irank+1)
          call p_recv (particles%lon (nparticles+1:nparticles+nparts_pe(irank)), irank, 100*irank+2)
          call p_recv (particles%lev (nparticles+1:nparticles+nparts_pe(irank)), irank, 100*irank+3)
          call p_recv (particles%pressure (nparticles+1:nparticles+nparts_pe(irank)), irank, 100*irank+4)
          call p_recv (particles%temperature (nparticles+1:nparticles+nparts_pe(irank)), irank, 100*irank+5)

          nparticles = nparticles + nparts_pe(irank)
 
       enddo

       ! lat,lon: from degree to radians
!! Evtl einschraenken auf den Bereich (1:nparticles)?!
       where (abs((particles%lon-mdi)/mdi)>eps) 
          particles%lon = particles%lon*(pi/180.)
          particles%lat = particles%lat*(pi/180.)       
       end where
       
    else

       ! send local field to rank 0
       call p_send (particles%lat (nparticles+1:nparticles+nparticles_added),0,100*rank+1)
       call p_send (particles%lon (nparticles+1:nparticles+nparticles_added),0,100*rank+2)
       call p_send (particles%lev (nparticles+1:nparticles+nparticles_added),0,100*rank+3)
       call p_send (particles%pressure (nparticles+1:nparticles+nparticles_added),0,100*rank+4)
       call p_send (particles%temperature (nparticles+1:nparticles+nparticles_added),0,100*rank+5)

    endif

    ! broadcast nparticles
    call p_bcast (nparticles, 0)

    if (rank==0) write (*,*) 'nparticles=',nparticles


    ! broadcast global arrays to all ranks
    call p_bcast (particles%lat,0)
    call p_bcast (particles%lon,0)
    call p_bcast (particles%lev,0)
    call p_bcast (particles%pressure,0)
    call p_bcast (particles%temperature,0)

    ! set default values for new particles
    particles%radius(nparticles_old+1:nparticles) = radius_default
    particles%density(nparticles_old+1:nparticles) = densnat_val
    particles%sedimentation(nparticles_old+1:nparticles) = 0.
    particles%tsv(nparticles_old+1:nparticles) = 0.
    particles%hno3(nparticles_old+1:nparticles) = hno3_background
    particles%h2o(nparticles_old+1:nparticles) = h2o_background
    particles%sice(nparticles_old+1:nparticles) = 0.
    particles%snat(nparticles_old+1:nparticles) = 0.
    particles%icebin(nparticles_old+1:nparticles) = 0.
    particles%natbin(nparticles_old+1:nparticles) = 0.
    particles%airparcel_density_change(nparticles_old+1:nparticles) = 1.
    particles%class(nparticles_old+1:nparticles) = 0.
    particles%particle_id(nparticles_old+1:nparticles) = -1.
    particles%snatmax(nparticles_old+1:nparticles) = 0.
    particles%sicemax(nparticles_old+1:nparticles) = 0.
    particles%tmin(nparticles_old+1:nparticles) = 0.
    ! clean up
    deallocate (nparts_pe)
    
  END SUBROUTINE clamssedi_get_particles

!--------------------------------------------------------------------
  
!--------------------------------------------------------------------
  SUBROUTINE clamssedi_update_particles

    ! BMIL
    USE messy_main_mpi_bi,       ONLY: p_recv, p_send, p_bcast

    ! SMCL
    USE messy_clams_global,      ONLY: rank, ntasks, mdi, eps
    USE messy_clamssedi_global,  ONLY: particles, nparticles, &
                                       part_ind_start, part_ind_end
    USE messy_clams_tools_utils, ONLY: pack_values

    implicit none

    logical, dimension(:), pointer :: existing(:)
    integer                        :: ipart_ind_start, ipart_ind_end, &
                                      irank, i   
    
    if (rank==0) then
       
       ! rank 0: receive particles from all ranks 
       do irank =  1, ntasks-1

          ipart_ind_start = ((nparticles/ntasks)* irank )  + 1
          ipart_ind_end = (nparticles/ntasks)*(irank+1) + ((irank+1)/ntasks)*modulo(nparticles,ntasks)
          
          call p_recv (particles%lat(ipart_ind_start:ipart_ind_end), irank, 100*irank+1)
          call p_recv (particles%lon(ipart_ind_start:ipart_ind_end), irank, 100*irank+2)
          call p_recv (particles%lev(ipart_ind_start:ipart_ind_end), irank, 100*irank+3)
          call p_recv (particles%pressure(ipart_ind_start:ipart_ind_end), irank, 100*irank+4)
          call p_recv (particles%temperature(ipart_ind_start:ipart_ind_end), irank, 100*irank+5)
          call p_recv (particles%radius(ipart_ind_start:ipart_ind_end), irank, 100*irank+6)
          call p_recv (particles%density(ipart_ind_start:ipart_ind_end), irank, 100*irank+7)
          call p_recv (particles%sedimentation(ipart_ind_start:ipart_ind_end), irank, 100*irank+8)
          call p_recv (particles%tsv(ipart_ind_start:ipart_ind_end), irank, 100*irank+9)
          call p_recv (particles%hno3(ipart_ind_start:ipart_ind_end), irank, 100*irank+10)
          call p_recv (particles%h2o(ipart_ind_start:ipart_ind_end), irank, 100*irank+11)
          call p_recv (particles%sice(ipart_ind_start:ipart_ind_end), irank, 100*irank+12)
          call p_recv (particles%snat(ipart_ind_start:ipart_ind_end), irank, 100*irank+13)
          call p_recv (particles%icebin(ipart_ind_start:ipart_ind_end), irank, 100*irank+14)
          call p_recv (particles%natbin(ipart_ind_start:ipart_ind_end), irank, 100*irank+15)
          call p_recv (particles%airparcel_density_change(ipart_ind_start:ipart_ind_end), irank, 100*irank+16)
          call p_recv (particles%class(ipart_ind_start:ipart_ind_end), irank, 100*irank+18)
          call p_recv (particles%particle_id(ipart_ind_start:ipart_ind_end), irank, 100*irank+19)
          call p_recv (particles%snatmax(ipart_ind_start:ipart_ind_end), irank, 100*irank+20)
          call p_recv (particles%sicemax(ipart_ind_start:ipart_ind_end), irank, 100*irank+21)
          call p_recv (particles%tmin(ipart_ind_start:ipart_ind_end), irank, 100*irank+22)
       enddo

       ! remove all particles with radius<=0 or density<=0

       allocate (existing(size(particles%lat)))
       existing = .false.
       where (particles%radius>0. .and. particles%density>0.) 
          existing = .true.
       end where

       ! remove all particles with lat, lon or lev = missing_value

       where (abs((particles%lat-mdi)/mdi)<eps)
          existing = .false.
       end where
       where (abs((particles%lon-mdi)/mdi)<eps)
          existing = .false.
       end where
       where (abs((particles%lev-mdi)/mdi)<eps)
          existing = .false.
       end where

       call pack_values (particles%lat, existing, 0._dp)
       call pack_values (particles%lon, existing, 0._dp)
       call pack_values (particles%lev, existing, 0._dp)
       call pack_values (particles%pressure, existing, 0._dp)
       call pack_values (particles%temperature, existing, 0._dp)
       call pack_values (particles%radius, existing, 0._dp)
       call pack_values (particles%density, existing, 0._dp)
       call pack_values (particles%sedimentation, existing, 0._dp)
       call pack_values (particles%tsv, existing, 0._dp)
       call pack_values (particles%hno3, existing, 0._dp)
       call pack_values (particles%h2o, existing, 0._dp)
       call pack_values (particles%sice, existing, 0._dp)
       call pack_values (particles%snat, existing, 0._dp)
       call pack_values (particles%icebin, existing, 0._dp)
       call pack_values (particles%natbin, existing, 0._dp)
       call pack_values (particles%airparcel_density_change, existing, 0._dp)
       call pack_values (particles%class, existing, 0._dp)
       call pack_values (particles%particle_id, existing, 0._dp)
       call pack_values (particles%snatmax, existing, 0._dp)
       call pack_values (particles%sicemax, existing, 0._dp)
       call pack_values (particles%tmin, existing, 0._dp)

       nparticles = count(existing)
       
       deallocate (existing)
       
       ! set particle id for new particles 
       do i = 1, nparticles
          if (particles%particle_id(i) <= 0) then
             part_id_max = part_id_max + 1
             particles%particle_id(i) = part_id_max
          endif
       enddo
       write (*,*) 'part_id_max=',part_id_max


    else ! other ranks: send particles to rank 0
       
       call p_send (particles%lat(part_ind_start:part_ind_end),0,100*rank+1)
       call p_send (particles%lon(part_ind_start:part_ind_end),0,100*rank+2)
       call p_send (particles%lev(part_ind_start:part_ind_end),0,100*rank+3)
       call p_send (particles%pressure(part_ind_start:part_ind_end),0,100*rank+4)
       call p_send (particles%temperature(part_ind_start:part_ind_end),0,100*rank+5)
       call p_send (particles%radius(part_ind_start:part_ind_end),0,100*rank+6)
       call p_send (particles%density(part_ind_start:part_ind_end),0,100*rank+7)
       call p_send (particles%sedimentation(part_ind_start:part_ind_end),0,100*rank+8)
       call p_send (particles%tsv(part_ind_start:part_ind_end),0,100*rank+9)
       call p_send (particles%hno3(part_ind_start:part_ind_end),0,100*rank+10)
       call p_send (particles%h2o(part_ind_start:part_ind_end),0,100*rank+11)
       call p_send (particles%sice(part_ind_start:part_ind_end),0,100*rank+12)
       call p_send (particles%snat(part_ind_start:part_ind_end),0,100*rank+13)
       call p_send (particles%icebin(part_ind_start:part_ind_end),0,100*rank+14)
       call p_send (particles%natbin(part_ind_start:part_ind_end),0,100*rank+15)
       call p_send (particles%airparcel_density_change(part_ind_start:part_ind_end),0,100*rank+16)
       call p_send (particles%class(part_ind_start:part_ind_end),0,100*rank+18)
       call p_send (particles%particle_id(part_ind_start:part_ind_end),0,100*rank+19)
       call p_send (particles%snatmax(part_ind_start:part_ind_end), 0, 100*rank+20)
       call p_send (particles%sicemax(part_ind_start:part_ind_end), 0, 100*rank+21)
       call p_send (particles%tmin(part_ind_start:part_ind_end), 0, 100*rank+22)
    endif


    ! broadcast global array particle to all ranks
    call p_bcast (particles%lat,0)
    call p_bcast (particles%lon,0)
    call p_bcast (particles%lev,0)
    call p_bcast (particles%pressure,0)
    call p_bcast (particles%temperature,0)
    call p_bcast (particles%radius,0)
    call p_bcast (particles%density,0)
    call p_bcast (particles%sedimentation,0)
    call p_bcast (particles%tsv,0)
    call p_bcast (particles%hno3,0)
    call p_bcast (particles%h2o,0)
    call p_bcast (particles%sice,0)
    call p_bcast (particles%snat,0)
    call p_bcast (particles%icebin,0)
    call p_bcast (particles%natbin,0)
    call p_bcast (particles%airparcel_density_change,0)
    call p_bcast (particles%class,0)
    call p_bcast (particles%particle_id,0)
    call p_bcast (particles%snatmax,0)
    call p_bcast (particles%sicemax,0)
    call p_bcast (particles%tmin,0)

    ! broadcast number of particles to all ranks
    call p_bcast (nparticles,0)
    call p_bcast (part_id_max,0)
    
  END SUBROUTINE clamssedi_update_particles

!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamssedi_update_airparcels

    ! BMIL
    USE messy_main_mpi_bi,      ONLY: p_bcast, p_sum

    ! SMCL
    USE messy_clams_global,     ONLY: nparts, dnparts, rank, ntasks
    
    implicit none

    real(prec), dimension(:), pointer     :: diff_hno3, diff_h2o
    real(prec), dimension(:), pointer     :: natbin_diff, icebin_diff
    integer,    dimension(:), allocatable :: nparts_pe
    integer                               :: startpos, irank, i


    ! Sum changes in HNO3 and H2O 
    allocate (diff_hno3(size(airparcels%diff_hno3)))  ! nparts
    allocate (diff_h2o (size(airparcels%diff_h2o)))
    diff_hno3 = p_sum(airparcels%diff_hno3)
    diff_h2o  = p_sum(airparcels%diff_h2o)

    ! Sum changes in NATbin and ICEbin 
    allocate (natbin_diff(size(airparcels%natbin_diff)))  ! nparts
    allocate (icebin_diff(size(airparcels%icebin_diff)))
    natbin_diff = p_sum(airparcels%natbin_diff)
    icebin_diff = p_sum(airparcels%icebin_diff)

    do i = 1, nparts
       if (natbin_diff(i) < 0. .or. icebin_diff(i) < 0.) then
          airparcels%natbin(i) = 0.
          airparcels%icebin(i) = 0.
          if (rank ==0) write(*,*) 'HALLO set airparcel nat/icebin ', i, 'to 0.'
       else
          airparcels%natbin(i) = min(airparcels%natbin(i) + natbin_diff(i), 1000.)
          airparcels%icebin(i) = min(airparcels%icebin(i) + icebin_diff(i), 1000.)
          ! ju_jug_160530 set minimum indices to zero
          airparcels%natbin(i) = max(airparcels%natbin(i) + natbin_diff(i), 0.)
          airparcels%icebin(i) = max(airparcels%icebin(i) + icebin_diff(i), 0.)
       endif

       if (airparcels%icebin(i) < 0.) then
          write(*,* ) 'warning in clamssedi_update_airparcels i,airparcels%icebin(i), diff:',i,airparcels%icebin(i),airparcels%icebin_diff(i)
       endif
    enddo

    ! get number of airparcels on all tasks
    allocate (nparts_pe(0:ntasks-1))
    nparts_pe = 0
    nparts_pe(rank) = dnparts
    do irank = 0, ntasks-1
       call p_bcast (nparts_pe(irank),irank)
    enddo
    !if (rank==0) write (*,*) 'nparts_pe:',nparts_pe

    ! Write changes on HNO3 and H2O to channel variables
    startpos = 1
    do irank = 1, rank
       startpos = startpos + nparts_pe(irank-1)
    enddo
    !write (*,*) 'in clamssedi_update_airparcels: rank, start, end',rank,startpos,startpos+dnparts-1
    
    HNO3(1:dnparts) = HNO3(1:dnparts) - diff_hno3(startpos:startpos+dnparts-1)
    H2O (1:dnparts) = H2O (1:dnparts) - diff_h2o (startpos:startpos+dnparts-1)
    
    NATbin (1:dnparts) = airparcels%natbin(startpos:startpos+dnparts-1)
    ICEbin (1:dnparts) = airparcels%icebin(startpos:startpos+dnparts-1)

    do i = 1, dnparts
       if (HNO3(i) < 0.) then
          HNO3(i) = hno3_minval
          cnt_neg_hno3 = cnt_neg_hno3 + 1
          write(*,*) 'WARNING: set HNO3 to hno3_minval!!! rank=', rank, ' cnt=', cnt_neg_hno3
       endif
       if (H2O(i) < 0.) then
          H2O(i) = h2o_minval
          cnt_neg_h2o = cnt_neg_h2o + 1
          write(*,*) 'WARNING: set H2O to h2o_minval!!! rank=', rank, ' cnt=', cnt_neg_h2o
       endif
    enddo

    deallocate(diff_hno3)
    deallocate(diff_h2o)
    deallocate(natbin_diff)
    deallocate(icebin_diff)
    deallocate(nparts_pe)    
    
  End Subroutine clamssedi_update_airparcels

!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamssedi_allocate

    USE messy_clams_global,     only: levelno, ntasks, rank
    use messy_clamssedi_global, only: part_ind_start, part_ind_end, &
                                      dnparticles, nparticles
    
    implicit none
    
    part_ind_start= ((nparticles/ntasks)* rank )  + 1
    part_ind_end  = (nparticles/ntasks)*(rank+1) + ((rank+1)/ntasks)*modulo(nparticles,ntasks)

    dnparticles = part_ind_end - part_ind_start + 1
    
    ! set in get_triangles_for_particle
    allocate(particles%natbin_diff(dnparticles))
    particles%natbin_diff = 0.
    allocate(particles%icebin_diff(dnparticles))
    particles%icebin_diff = 0.
    allocate(particles%tr_ind_up(dnparticles))
    allocate(particles%tr_ind_down(dnparticles))
    allocate(particles%lev_up(dnparticles))
    allocate(particles%lev_down(dnparticles))

    ! Used for trajectory calculation
    allocate (levelno(nparticles))
    levelno = 1


  END SUBROUTINE clamssedi_allocate

!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamssedi_deallocate

    USE messy_clams_global,     ONLY: rank, levelno
    USE messy_clamssedi_global, ONLY: airparcels

    IMPLICIT NONE

    deallocate (airparcels%lat)
    deallocate (airparcels%lon)
    deallocate (airparcels%lev)
    deallocate (airparcels%hno3)
    deallocate (airparcels%h2o)
    deallocate (airparcels%diff_hno3)
    deallocate (airparcels%diff_h2o)
    deallocate (airparcels%natbin)
    deallocate (airparcels%icebin)
    deallocate (airparcels%natbin_diff)
    deallocate (airparcels%icebin_diff)
    deallocate (airparcels%coor)
    deallocate (airparcels%ntriang)
    deallocate (airparcels%triang_ids)
    deallocate (airparcels%ilev)

    deallocate (nairparcels_lev)
    deallocate (nairparcels50_lev)
    deallocate (triangles%airparcel_indices)
    deallocate (triangles%nparticles_nat)
    deallocate (triangles%nparticles_ice)
    
    deallocate(particles%tr_ind_up)
    deallocate(particles%tr_ind_down)
    deallocate(particles%lev_up)
    deallocate(particles%lev_down)
    deallocate(particles%icebin_diff)
    deallocate(particles%natbin_diff)

    ! used for trajectory calculation
    deallocate (levelno)

  END SUBROUTINE clamssedi_deallocate
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamssedi_error_handler (status, substr)

    USE messy_main_blather_bi, ONLY: error_bi

    IMPLICIT NONE

    INTEGER      :: status
    CHARACTER(*) :: substr

    CHARACTER(80) :: errorstr

    SELECT CASE (status)
    CASE (101)
       errorstr = 'Error in clamssedi_get_airparcels: sum(dnparts) /= nparts !!'
    CASE (201)
       errorstr = 'Error in clamssedi_get_particles: sum(nparticles_addes) > nparticles_max !!'
    END SELECT

    CALL error_bi(errorstr,substr)

  END SUBROUTINE clamssedi_error_handler

!--------------------------------------------------------------------
#endif
END MODULE messy_clamssedi_si
