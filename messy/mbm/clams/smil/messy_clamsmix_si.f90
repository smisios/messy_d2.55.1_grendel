!**********************************************************************
MODULE messy_clamsmix_si

#if defined(ECHAM5) || defined(MBM_CLAMS)
!**********************************************************************
!  Submodel interface for clamsmix 
!**********************************************************************
! SMCL
  USE messy_clams_global,     ONLY: prec, dp, species_type, specnamelen, &
                                    paramtype
  USE messy_clamsmix

  USE messy_main_timer_event, ONLY: time_event, io_time_event, event_is_active 

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

! op_pj_20160606+
!!$  TYPE(time_event),    PUBLIC :: mixevent
!!$  TYPE(io_time_event), PUBLIC :: io_mixevent
  TYPE(time_event),    PUBLIC, SAVE :: mixevent
  TYPE(io_time_event), PUBLIC, SAVE :: io_mixevent
! op_pj_20160606-

 ! MODULE VARIABLES
  REAL(PREC),         DIMENSION(:), POINTER :: LAT         => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LON         => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: ZETA        => NULL()
  REAL(DP),                         POINTER :: JULTIME
  REAL(PREC),         DIMENSION(:), POINTER :: LAT_OLD => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LON_OLD => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LEV_OLD => NULL()
  REAL(DP),                         POINTER :: JULTIME_OLD
  REAL(PREC),         DIMENSION(:), POINTER :: LAT_OLD_MIX => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LON_OLD_MIX => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: STATE       => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: STATE_VERT  => NULL()
  TYPE(species_type), DIMENSION(:), POINTER :: MIXSPECARR
  REAL(PREC),         DIMENSION(:), POINTER :: THETA       => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: BVF      => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: THETA_OLD_MIX  => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: BVF_OLD_MIX => NULL()

  ! Shuffle arrays:
  REAL(PREC),         DIMENSION(:), POINTER :: ZETA_SHUFFLED        => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LAT_SHUFFLED         => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LON_SHUFFLED         => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LAT_OLD_MIX_SHUFFLED => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LON_OLD_MIX_SHUFFLED => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: STATE_SHUFFLED       => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: STATE_VERT_SHUFFLED  => NULL()
  TYPE(species_type), DIMENSION(:), POINTER :: MIXSPECARR_SHUFFLED 
  REAL(PREC),         DIMENSION(:), POINTER :: THETA_SHUFFLED          => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: BVF_SHUFFLED         => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: THETA_OLD_MIX_SHUFFLED  => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: BVF_OLD_MIX_SHUFFLED => NULL()

  REAL(PREC),         DIMENSION(:), POINTER :: LEV_GRID   => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LEV_DELTA  => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: R_GRID     => NULL()
 

  TYPE(paramtype),    DIMENSION(:), POINTER :: PARAM     ! from CLAMS
  TYPE(paramtype),    DIMENSION(:), POINTER :: PARAM_OLD ! from CLAMS

  
  CHARACTER(LEN=specnamelen), dimension(:), allocatable :: helpstrarr
 
  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: clamsmix_setup
  PUBLIC :: clamsmix_initialize
  PUBLIC :: clamsmix_init_memory
  PUBLIC :: clamsmix_init_coupling
  PUBLIC :: clamsmix_global_start
  PUBLIC :: clamsmix_global_end
  PUBLIC :: clamsmix_free_memory
  
  PUBLIC :: clamsmix_reshuffle   ! used in BMIX
 
!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
  SUBROUTINE clamsmix_setup

    ! BMIL
    USE messy_main_blather_bi,     ONLY: error_bi

    ! SMCL
    USE messy_main_tools,          ONLY: find_next_free_unit
    USE messy_clams_global,        ONLY: rank, initfile, init_vertcoorname, &
                                         ldiagout
    USE messy_clamsmix_global,     ONLY: adapt_par, timestep_mix, &
                                         no_steps, lexp,  &
                                         fac_limit_outside, fac_limit_inside, &
                                         fac_limit_lev_down, fac_limit_lev_up, &
                                         fac_eliminate, fac_bvf_min, &
                                         delta_lev, r_dev
    USE messy_clamsmix_lib_io,     ONLY: read_init_mix_config

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsmix_setup'
    integer :: status, iou

!!!!! adapt_par%nlevs wird bereits in main_channel_initialize_dims genutzt !!!

    ! prepare adapt_par for vertical grid info...
    nullify(adapt_par%lev_grid,adapt_par%lev_delta)
    nullify(adapt_par%r_grid)

    ! read global attributes and grid information to adapt_par
    CALL read_init_mix_config (status, initfile, &
            adapt_par%lev_min, adapt_par%lev_max, &
            adapt_par%nlevs, &
            adapt_par%lat_down, adapt_par%lat_up, &
            adapt_par%lat_min, adapt_par%lat_max, &
            adapt_par%r_mean_c, adapt_par%r_mean_h, &
            adapt_par%lev_grid,adapt_par%lev_delta, &
            adapt_par%r_grid,init_vertcoorname)
    IF (status /= 0) CALL error_bi('Error in sub. read_global_mix_atts!',substr)

    ! Read namelist:
    iou = find_next_free_unit(100,200)
    CALL clamsmix_read_nml(status, iou)
    IF (status /= 0) CALL error_bi('Error in clamsmix_read_nml ',substr)

    ! Save namelist parameters on adapt_par
    adapt_par%lexp = lexp
    adapt_par%fac_limit_outside  = fac_limit_outside
    adapt_par%fac_limit_inside   = fac_limit_inside
    adapt_par%fac_limit_lev_down = fac_limit_lev_down
    adapt_par%fac_limit_lev_up   = fac_limit_lev_up
    adapt_par%fac_eliminate      = fac_eliminate
    adapt_par%fac_bvf_min        = fac_bvf_min
    adapt_par%delta_lev          = delta_lev
    adapt_par%no_steps           = no_steps
    adapt_par%r_dev              = r_dev
    
    adapt_par%timestep = timestep_mix * 3600
     
    ! Control output
    if (rank == 0  .and. ldiagout) then
       write(*,*) 'ADAPT_PAR:'
       write(*,*) 'no_steps', no_steps
       write(*,*) 'nlevs', adapt_par%nlevs
       write(*,*) 'lat_up', adapt_par%lat_up
       write(*,*) 'lat_down', adapt_par%lat_down
       write(*,*) 'lat_min', adapt_par%lat_min
       write(*,*) 'lat_max', adapt_par%lat_max
       write(*,*) 'lev_min', adapt_par%lev_min
       write(*,*) 'lev_max', adapt_par%lev_max
       write(*,*) 'lexp', adapt_par%lexp
       write(*,*) 'time_step', adapt_par%timestep
       write(*,*) 'fac_limit_outside', adapt_par%fac_limit_outside
       write(*,*) 'fac_limit_inside', adapt_par%fac_limit_inside
       write(*,*) 'fac_limit_lev_down', adapt_par%fac_limit_lev_down
       write(*,*) 'fac_limit_lev_up', adapt_par%fac_limit_lev_up
       write(*,*) 'fac_eliminate', adapt_par%fac_eliminate
       write(*,*) 'fac_bvf_min', adapt_par%fac_bvf_min
       write(*,*) 'adapt_par%r_mean_c', adapt_par%r_mean_c
       write(*,*) 'adapt_par%r_mean_h', adapt_par%r_mean_h
       write(*,*) 'delta_lev', adapt_par%delta_lev
       write(*,*) 'r_dev', adapt_par%r_dev
       write(*,*) 'lev_grid', adapt_par%lev_grid
       write(*,*) 'lev_delta', adapt_par%lev_delta
       
    end if

  END SUBROUTINE clamsmix_setup

!-----------------------------------------------------------------------

  SUBROUTINE clamsmix_initialize

    ! BMIL
    USE messy_main_blather_bi,     ONLY: error_bi
    USE messy_main_timer_bi,       ONLY: timer_event_init
    USE messy_main_mpi_bi,         ONLY: p_pe

    ! SMCL
    USE messy_main_timer,          ONLY: delta_time, current_date, &
                                         YEAR_START, MONTH_START, DAY_START, &
                                         HOUR_START
    USE messy_clams_global,        ONLY: initfile, nspec, nchemspec, maxspec, &
                                         SPECARR, species_type, specnames, &
                                         nparams, paramnames, &
                                         nparams_old, paramnames_old, &
                                         buffersize
    USE messy_clamsmix_global,     ONLY: timestep_mix, nmixspec, mixspec, &
                                         switch_mixing, vert_mix_param
    USE messy_clamsmix,            ONLY: clamsmix_read_nml
    USE messy_clams_tools_utils,   ONLY: uppercase, str_found
    USE messy_clams_tools_ncutils, ONLY: nc_check_error

    USE netcdf

    IMPLICIT NONE
    
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsmix_initialize'

    TYPE(species_type), DIMENSION(:), POINTER :: helpspecarr
    
    INTEGER            :: status, nhelpspec
    integer            :: i, imix
    integer            :: ncid, varid, rcode
    logical            :: reading, found

    INTRINSIC :: TRIM, ADJUSTL

    IF (p_pe==0) THEN
       WRITE(*,*)
       WRITE(*,*) uppercase(substr)
    ENDIF

!!!!!    
!    ALLOCATE (mixspec(maxspec))
    
    if (mod(timestep_mix*3600,int(delta_time)) /= 0) then
       call error_bi ("Wrong mixing timestep !!!",substr)
    endif

    ! Define MIX event:
    io_mixevent%counter = timestep_mix   
    io_mixevent%unit = 'hours'
    io_mixevent%adjustment = 'exact'
    io_mixevent%offset = -delta_time
    CALL timer_event_init (mixevent, io_mixevent, 'MIX_Event', 'present')


    ! If vertical mixing is switched on:
    !       THETA and BVF_WET (or BVF) must be specified as parameters
    if (switch_mixing==2) then
       if (.not. str_found (nparams, paramnames, 'THETA')) then
          call error_bi ('If vertical mixing is switched on, '// &
               'THETA must be specified as parameters in clams.nml !!!', &
               substr)
       endif
       if (uppercase(vert_mix_param)=='WET') then
          if (.not. str_found (nparams, paramnames, 'BVF_WET')) then
             call error_bi ('If vertical mixing (wet) is switched on, '// &
                  'BVF_WET must be specified as parameters in clams.nml !!!', &
                  substr)
          endif
       elseif (.not. str_found (nparams, paramnames, 'BVF')) then
          call error_bi ('If vertical mixing (dry) is switched on, '// &
               'BVF must be specified as parameters in clams.nml !!!', &
               substr)
       endif
    endif


    ! Add species to SPECARR
    rcode = nf90_open (initfile, nf90_nowrite, ncid, buffersize)
    do imix = 1, nmixspec
       found = .false.
       i = 1
       do while (.not. found .and. i<=nspec)
          if (mixspec(imix)==SPECARR(i)%name) found = .true. 
          i = i + 1
       enddo
       if (.not. found) then
          nspec = nspec + 1
          if (nspec > maxspec) &
               call error_bi (substr,"Too many species: increase MAXSPEC !!!")
          specnames(nspec)        = mixspec(imix)
          SPECARR(nspec)%name     = mixspec(imix)
          rcode = nf90_inq_varid (ncid,trim(SPECARR(nspec)%name),varid)
          if (rcode == 0) then
             rcode = nf90_get_att (ncid,varid,"units",SPECARR(nspec)%units)
             if (rcode/=0) SPECARR(nspec)%units = " "
             rcode = nf90_get_att (ncid,varid,"long_name",SPECARR(nspec)%longname)
             if (rcode/=0) SPECARR(nspec)%longname = SPECARR(nspec)%name
          else
             SPECARR(nspec)%units = " "
             SPECARR(nspec)%longname = SPECARR(nspec)%name
          endif
          SPECARR(nspec)%ctype    = " "
      endif
    enddo
    rcode = nf90_close (ncid)

    if (p_pe==0) then
       write (*,*) 'NSPEC=', NSPEC
       write (*,*) 'NCHEMSPEC=', NCHEMSPEC
       do i = 1, nspec
          write (*,*) specarr(i)%name,'  ', specarr(i)%ctype,'  ', trim(specarr(i)%units)
       enddo
    endif
    ! op_sb_20191021 (moved to init_memory)
!!$    ALLOCATE (MIXSPECARR(nspec))
!!$    ALLOCATE (MIXSPECARR_SHUFFLED(nspec))


    ! aus clams_si: allocate param array:
    IF (nparams > 0) THEN
       ALLOCATE (PARAM(nparams))
       DO i = 1, nparams
          PARAM(i)%name = paramnames(i)
      ENDDO
    ENDIF

    IF (nparams_old > 0) THEN
       ALLOCATE (PARAM_OLD(nparams_old))
       DO i = 1, nparams_old
          PARAM_OLD(i)%name = trim(paramnames_old(i))//'_OLD'
      ENDDO
    ENDIF


    
  END SUBROUTINE clamsmix_initialize

!-----------------------------------------------------------------------

  SUBROUTINE clamsmix_init_memory

    ! BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: REPR_CLAMS_SHUFFLE
    USE messy_main_mpi_bi,           ONLY: p_pe

    ! SMCL
    USE messy_main_channel,     ONLY: new_channel, new_attribute, &
                                      new_channel_object
    USE messy_clams_global,     ONLY: mdi, username, nspec, SPECARR
    USE messy_clams_tools_utils,ONLY: uppercase
    USE messy_clamsmix_global,  ONLY: switch_mixing
    USE messy_clamsmix,         ONLY: modstr, modver

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, TRIM

    ! LOCAL 
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsmix_init_memory'
    INTEGER       :: status
    CHARACTER(40) :: cr_date
    CHARACTER(8)  :: ydate
    CHARACTER(10) :: ytime
    INTEGER       :: ispec

    IF (p_pe==0) WRITE(*,*) uppercase(substr)

    CALL DATE_AND_TIME(ydate, ytime)
    WRITE(cr_date,'(A4,5(A,A2))')   &
         ydate(1:4),'-',ydate(5:6),'-',ydate(7:8),' ',  &
         ytime(1:2),':',ytime(3:4),':',ytime(5:6) 

    ! op_sb_20191021+ (nspec updated in clams_init_memory with number of MESSy tracer)
    ALLOCATE (MIXSPECARR(nspec))
    ALLOCATE (MIXSPECARR_SHUFFLED(nspec))
    ! op_sb_20191021-

     
    ! Define channel MIX
    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, modstr//'_version', c = modver)
    CALL channel_halt(substr, status)

    ! SHUFFLED channel objects
    CALL new_channel_object(status, modstr, 'ZETA_SHUFFLED', p1=ZETA_SHUFFLED, &
                            reprid=REPR_CLAMS_SHUFFLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'LAT_SHUFFLED', p1=LAT_SHUFFLED, &
                            reprid=REPR_CLAMS_SHUFFLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'LON_SHUFFLED', p1=LON_SHUFFLED, &
                            reprid=REPR_CLAMS_SHUFFLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'LAT_OLD_MIX_SHUFFLED', &
                            p1=LAT_OLD_MIX_SHUFFLED, reprid=REPR_CLAMS_SHUFFLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'LON_OLD_MIX_SHUFFLED', &
                            p1=LON_OLD_MIX_SHUFFLED, reprid=REPR_CLAMS_SHUFFLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'STATE_SHUFFLED', &
                            p1=STATE_SHUFFLED, reprid=REPR_CLAMS_SHUFFLE)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'STATE_VERT_SHUFFLED', &
                            p1=STATE_VERT_SHUFFLED, reprid=REPR_CLAMS_SHUFFLE)
    CALL channel_halt(substr, status)
    DO ispec = 1, nspec
       CALL new_channel_object(status, modstr, trim(SPECARR(ispec)%name)//'_SHUFFLED', &
            p1=MIXSPECARR_SHUFFLED(ispec)%values, reprid=REPR_CLAMS_SHUFFLE)
       CALL channel_halt(substr, status)
    ENDDO

    IF (switch_mixing==2) THEN
       CALL new_channel_object(status, modstr, 'THETA_SHUFFLED', &
                    p1=THETA_SHUFFLED, reprid=REPR_CLAMS_SHUFFLE)
       CALL channel_halt(substr, status)
       CALL new_channel_object(status, modstr, 'THETA_OLD_MIX_SHUFFLED', &
                    p1=THETA_OLD_MIX_SHUFFLED, reprid=REPR_CLAMS_SHUFFLE)
       CALL channel_halt(substr, status)
       CALL new_channel_object(status, modstr, 'BVF_SHUFFLED', &
                    p1=BVF_SHUFFLED, reprid=REPR_CLAMS_SHUFFLE)
       CALL channel_halt(substr, status)
       CALL new_channel_object(status, modstr, 'BVF_OLD_MIX_SHUFFLED', &
                    p1=BVF_OLD_MIX_SHUFFLED, reprid=REPR_CLAMS_SHUFFLE)
       CALL channel_halt(substr, status)
    ENDIF
   

  END SUBROUTINE clamsmix_init_memory

!-----------------------------------------------------------------------

  SUBROUTINE clamsmix_init_coupling

    ! BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_mpi_bi,           ONLY: p_pe

    ! SMCL
    USE messy_main_channel,       ONLY: get_channel_object

    USE messy_clams_global,       ONLY: nspec, SPECARR, &
                                        nparams, nparams_old
    USE messy_clamsmix_global,    ONLY: switch_mixing, adapt_par, &
                                        vert_mix_param
    USE messy_clams_tools_utils,  ONLY: uppercase

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsmix_init_coupling'
    INTEGER  :: status
    INTEGER  :: i

    IF (p_pe==0) WRITE(*,*) uppercase(substr)

    ! Get arrays from CLAMS submodel:

    ! get current airparcel positions:
    CALL get_channel_object(status, 'clams', 'LAT', p1=LAT)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON', p1=LON)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LEV', p1=ZETA)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'JULTIME', p0=JULTIME)
    CALL channel_halt(substr, status)

    ! get airparcel positions at last CHEM call
    CALL get_channel_object(status, 'clams', 'LAT_OLD', p1=LAT_OLD)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON_OLD', p1=LON_OLD)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LEV_OLD', p1=LEV_OLD)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'JULTIME_OLD', p0=JULTIME_OLD)
    CALL channel_halt(substr, status)

    ! get airparcel positions at last MIX call
    CALL get_channel_object(status, 'clams', 'LAT_OLD_MIX', p1=LAT_OLD_MIX)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON_OLD_MIX', p1=LON_OLD_MIX)
    CALL channel_halt(substr, status)

    ! Get species (from CLAMS):
    DO i = 1, nspec
       CALL get_channel_object(status, 'clams', &
            trim(SPECARR(i)%name), p1=MIXSPECARR(i)%values)
       CALL channel_halt(substr, status)
       MIXSPECARR(i)%name     = SPECARR(i)%name
       MIXSPECARR(i)%longname = SPECARR(i)%longname
       MIXSPECARR(i)%units    = SPECARR(i)%units
       MIXSPECARR(i)%ctype    = SPECARR(i)%ctype
    ENDDO

    ! STATE, STATE_VERT
    CALL get_channel_object(status, 'clams', 'STATE', p1=STATE)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'STATE_VERT', p1=STATE_VERT)
    CALL channel_halt(substr, status)

    ! THETA, BVF_WET
!!!!! Wenn switch_mixing==2: BVF und THETA mit Channelobjecten koppeln
!!!!! sonst: Felder fuer BVF und THETA anlegen und mit 0. initialisieren
    IF (switch_mixing==2) THEN
       CALL get_channel_object(status, 'clams', 'THETA', p1=THETA)
       CALL channel_halt(substr, status)
       if (uppercase(vert_mix_param) == 'WET') then
          CALL get_channel_object(status, 'clams', 'BVF_WET', p1=BVF)
       else
          CALL get_channel_object(status, 'clams', 'BVF', p1=BVF)
       endif
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'clams', 'THETA_OLD_MIX', p1=THETA_OLD_MIX)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'clams', 'BVF_OLD_MIX', p1=BVF_OLD_MIX)
       CALL channel_halt(substr, status)
    ENDIF

    ! LEV_GRID, LEV_DELTA and R_GRID
    CALL get_channel_object(status, 'clams', 'LEV_GRID', p1=LEV_GRID)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LEV_DELTA', p1=LEV_DELTA)
    CALL channel_halt(substr, status)
    IF (ASSOCIATED(adapt_par%r_grid)) THEN
       CALL get_channel_object(status, 'clams', 'R_GRID', p1=R_GRID)
       CALL channel_halt(substr, status)
    ENDIF

    ! Couple Parameters (clams.nml)
    DO i = 1, nparams
       CALL get_channel_object(status, 'clams',  &
                   trim(PARAM(i)%name), p1=PARAM(i)%values)
       CALL channel_halt(substr, status)
    END DO

    ! TEMP_OLD, PRESS_OLD
    DO i = 1, nparams_old
       !write (*,*) 'couple ', trim(PARAM_OLD(i)%name)
       CALL get_channel_object(status, 'clams',  &
                   trim(PARAM_OLD(i)%name), p1=PARAM_OLD(i)%values)
       CALL channel_halt(substr, status)
    END DO

  END SUBROUTINE clamsmix_init_coupling
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clamsmix_global_start

    USE messy_main_blather_bi,     ONLY: error_bi
    USE messy_main_timer,          ONLY: lfirst_cycle, lresume
    USE messy_clams_global,        ONLY: rank, grid_switch_co, ldiagout
    USE messy_clamsmix_global,     ONLY: grid_switch, adapt_par

    USE netcdf

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsmix_global_start'
    INTEGER :: status

    IF (rank==0 .and. ldiagout) WRITE(*,*) substr

    IF (lresume) THEN
       grid_switch = grid_switch_co
    ELSE
       grid_switch_co = grid_switch
    ENDIF

    IF (lfirst_cycle) THEN

       IF (rank==0) write(*,*) substr, ':  in lfirst_cycle'

       ! set channel objects LEV_GRID, LEV_DELTA and R_GRID
       LEV_GRID = adapt_par%lev_grid
       LEV_DELTA = adapt_par%lev_delta
       IF (ASSOCIATED(adapt_par%r_grid))  R_GRID = adapt_par%r_grid

    ENDIF

  END SUBROUTINE clamsmix_global_start

!-----------------------------------------------------------------------

  SUBROUTINE clamsmix_global_end

    ! BMIL
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_mpi_bi,           ONLY: p_pe, p_nprocs, p_sum, p_bcast

    ! SMCL: MESSy
    USE messy_main_channel,    ONLY: set_channel_output
    USE messy_main_timer,      ONLY: current_date, &
                                     YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, &
                                     delta_time
    
    USE messy_main_switch,     ONLY: USE_CLAMSBMIX

    ! SMCL: CLaMS
    USE messy_clams_global,    ONLY: nparts, nparts_max, dnparts, dnparts_max,  &
                                     nspec,  asad_gfirst, &
                                     mdi, lmixevent, lbmixevent, dnparts_max_shuffle, &
                                     dnparts_co, grid_switch_co, ldiagout, specnames, &
                                     lcoupled, &
                                     nparams, paramnames, &
                                     nparams_old, paramnames_old, &
                                     fut_year, fut_month, fut_day, fut_sec, &
                                     irdsec, irdday, dates30, timetype
    USE messy_clamsmix_global, ONLY: grid_switch, l_nlevs, adapt_par, &
                                     levelrange, vert_mix_param, &
                                     lev_min_act, lev_max_act, lev_delta_act,&
                                     l_min_act, l_max_act, l_delta_act, switch_mixing
    USE messy_clamsmix,        ONLY: mix,  &
                                     clamsmix_set_config, clamsmix_set_bounds, &
                                     clamsmix_nparts_per_level, clamsmix_set_level
    USE messy_clams_tools_interpolreg,  ONLY: interpolate_param
    USE messy_clams_tools_dateconv,     ONLY: incdat_interface
    USE messy_clams_tools_utils,        ONLY: str_pos, uppercase
! op_pj_20170110+
!!$ USE messy_clamschem_global,ONLY: asad_gfirst
!    USE messy_clams_global,    ONLY: asad_gfirst
! op_pj_20170110-

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsmix_global_end'
    INTEGER :: status, irank, i, ilev, ispec
    integer :: p_result, p_error
    integer :: my_nparts

    integer,    dimension(:),    pointer :: indarr        ! (1:dnparts)
    integer,    dimension(:,:,:),pointer :: irange        ! (2,nlev,0:ntasks-1)
    integer,    dimension(:,:),  pointer :: dnparts_level ! (nlev,0:ntasks-1)
    integer,    dimension(:),    pointer :: nparts_level  ! (nlev)

    REAL(PREC),         DIMENSION(:), POINTER :: ZETA_SHUFFLED_TEMP        => NULL()
    REAL(PREC),         DIMENSION(:), POINTER :: LAT_SHUFFLED_TEMP         => NULL()
    REAL(PREC),         DIMENSION(:), POINTER :: LON_SHUFFLED_TEMP         => NULL()
    REAL(PREC),         DIMENSION(:), POINTER :: LAT_OLD_MIX_SHUFFLED_TEMP => NULL()
    REAL(PREC),         DIMENSION(:), POINTER :: LON_OLD_MIX_SHUFFLED_TEMP => NULL()
    TYPE(species_type), DIMENSION(:), POINTER :: MIXSPECARR_SHUFFLED_TEMP  => NULL()

    REAL(PREC),         DIMENSION(:), POINTER :: THETA_SHUFFLED_TEMP          => NULL()
    REAL(PREC),         DIMENSION(:), POINTER :: THETA_OLD_MIX_SHUFFLED_TEMP  => NULL()
    REAL(PREC),         DIMENSION(:), POINTER :: BVF_SHUFFLED_TEMP         => NULL()
    REAL(PREC),         DIMENSION(:), POINTER :: BVF_OLD_MIX_SHUFFLED_TEMP => NULL()

    logical, dimension(:), pointer :: mask
    
    REAL(PREC)     :: pi

    TYPE(timetype) :: itime, ipasttime

    
    pi=4.*atan(1.)                                                        

!!!!! ???
!    IF (lbmixevent) switch_mixing = .FALSE.
!    IF (lmixevent)  switch_mixing = .TRUE.

    IF (lmixevent .OR. lbmixevent) THEN

       IF (p_pe==0) THEN
          WRITE (*,*)
          WRITE(*,*) 'Active mixevent'
          WRITE (*,*)
       ENDIF

       IF (p_pe==0 .and. ldiagout) write (*,*) 'nparts_max vor mix:',nparts_max 
       !write (*,*) 'vor mix: p_pe, dnparts=',p_pe,dnparts

       ! Set nlevs, l_nlevs, lev_min, lev_max
       IF (p_pe==0 .and. ldiagout) write (*,*) 'call clamsmix_set_config'
       call clamsmix_set_config (status)
       IF (status /= 0) CALL error_bi('Error in clamsmix_set_config !!!',substr)

       ! Set boundaries (lev_min_act,lev_max_act,l_min_act,l_max_act)
       ! and adapt_par arrays for all levels
       IF (p_pe==0 .and. ldiagout) write (*,*) 'call clamsmix_set_bounds'
       call clamsmix_set_bounds 

       if (dnparts > 0) then
          ! Index sort of ZETA (-> indarr) 
          ! Get start- and endposition of each zetalevel in sorted array (-> irange)
          IF (p_pe==0 .and. ldiagout) write (*,*) 'call clamsmix_sort'
          call clamsmix_sort (status, ZETA(1:dnparts), indarr, irange)
          IF (status /= 0) CALL error_bi('Error in clamsmix_sort !!!',substr)
       else! dnparts=0
          allocate (irange(2,l_nlevs,0:p_nprocs-1))
          irange(:,:,p_pe) = -1

          ! Dummy values for indarr
          allocate(indarr(l_nlevs))
          indarr = -1
       endif
       
       ! Broadcast irange
       do i = 0, p_nprocs-1
          call p_bcast (irange(:,:,i),i)
       enddo
           
       ! Get number of particles per level
       ! nparts_level (1:l_nlevs)  : total number of particles on each level (for all ranks)
       ! dnparts_level(1:l_nlevs,i): number of particles on each level on rank i
       IF (p_pe==0 .and. ldiagout) write (*,*) 'call clamsmix_nparts_per_level'
       call clamsmix_nparts_per_level (irange, dnparts_level, nparts_level)


       ! Set start- and endlevel for each task
       IF (p_pe==0 .and. ldiagout) write (*,*) 'call clamsmix_set_level'
       call clamsmix_set_level (levelrange, nparts_level)
       !write (*,*) 'rank, startlevel, endlevel=',p_pe,levelrange(1,p_pe),levelrange(2,p_pe)

!!!!! status = 201: Felder nicht gross genug!!! 

       ! Shuffle all positions and species:
       ! Each rank gets all points between levelrange(1,irank) and levelrange(2,irank)
       if (p_pe==0 .and. ldiagout) write (*,*) 'call clamsmix_shuffle '
       call clamsmix_shuffle (status, my_nparts, indarr, irange, levelrange, &
            dnparts_level, nparts_level, ZETA, ZETA, ZETA_SHUFFLED_TEMP)
       IF (status /= 0) call clamsmix_error_handler (status, substr)

       call clamsmix_shuffle (status, my_nparts, indarr, irange, levelrange, &
            dnparts_level, nparts_level, ZETA, LAT, LAT_SHUFFLED_TEMP)
       IF (status /= 0) call clamsmix_error_handler (status, substr)
 
       call clamsmix_shuffle (status, my_nparts, indarr, irange, levelrange, &
            dnparts_level, nparts_level, ZETA, LON, LON_SHUFFLED_TEMP)
       IF (status /= 0) call clamsmix_error_handler (status, substr)

       call clamsmix_shuffle (status, my_nparts, indarr, irange, levelrange, &
            dnparts_level, nparts_level, ZETA, LAT_OLD_MIX, LAT_OLD_MIX_SHUFFLED_TEMP)
       IF (status /= 0) call clamsmix_error_handler (status, substr)

       call clamsmix_shuffle (status, my_nparts, indarr, irange, levelrange, &
            dnparts_level, nparts_level, ZETA, LON_OLD_MIX, LON_OLD_MIX_SHUFFLED_TEMP)
       IF (status /= 0) call clamsmix_error_handler (status, substr)

       allocate (MIXSPECARR_SHUFFLED_TEMP(size(MIXSPECARR)))
       do i = 1, nspec
          MIXSPECARR_SHUFFLED_TEMP(i)%name = MIXSPECARR(i)%name
          call clamsmix_shuffle (status, my_nparts, indarr, irange, levelrange, &
               dnparts_level, nparts_level, ZETA, MIXSPECARR(i)%values, MIXSPECARR_SHUFFLED_TEMP(i)%values)
          IF (status /= 0) call clamsmix_error_handler (status, substr)
       enddo

       ! If vertical mixing is switched on:
       !    Set SHUFFLE arrays for THETA, BVF, THETA_OLD, BVF_OLD
       IF (switch_mixing==2) THEN
          call clamsmix_shuffle (status, my_nparts, indarr, irange, levelrange, &
               dnparts_level, nparts_level, ZETA, THETA, THETA_SHUFFLED_TEMP)
          IF (status /= 0) call clamsmix_error_handler (status, substr)
          call clamsmix_shuffle (status, my_nparts, indarr, irange, levelrange, &
               dnparts_level, nparts_level, ZETA, THETA_OLD_MIX, THETA_OLD_MIX_SHUFFLED_TEMP)
          IF (status /= 0) call clamsmix_error_handler (status, substr)
          call clamsmix_shuffle (status, my_nparts, indarr, irange, levelrange, &
               dnparts_level, nparts_level, ZETA, BVF, BVF_SHUFFLED_TEMP)
          IF (status /= 0) call clamsmix_error_handler (status, substr)
          call clamsmix_shuffle (status, my_nparts, indarr, irange, levelrange, &
               dnparts_level, nparts_level, ZETA, BVF_OLD_MIX, BVF_OLD_MIX_SHUFFLED_TEMP)
          IF (status /= 0) call clamsmix_error_handler (status, substr)

       ! without vertical mixing:
       !    create dummy arrays for THETA, BVF, THETA_OLD, BVF_OLD 
       !    and initialize these arrays (lenght=1) with 0.
       ELSE
          allocate (THETA_SHUFFLED_TEMP(1))
          allocate (THETA_OLD_MIX_SHUFFLED_TEMP(1))
          allocate (BVF_SHUFFLED_TEMP(1))
          allocate (BVF_OLD_MIX_SHUFFLED_TEMP(1))
          THETA_SHUFFLED_TEMP = 0.
          THETA_OLD_MIX_SHUFFLED_TEMP = 0.
          BVF_SHUFFLED_TEMP = 0.
          BVF_OLD_MIX_SHUFFLED_TEMP = 0.
       ENDIF


       if (p_pe==0 .and. ldiagout) write (*,*) 'After clamsmix_shuffle '

       ! Set STATE to 0
       STATE      = 0
       STATE_VERT = 0
       STATE_SHUFFLED      = 0
       STATE_VERT_SHUFFLED = 0

!!!!!
       ! ispec = str_pos(nspec, specnames, 'DC')
       ! write (*,*) 'str_pos(DC)=', ispec
       ! allocate (mask(dnparts_max_shuffle))
       ! mask = .false.
       ! where (mixspecarr_shuffled_temp(ispec)%values==1)
       !    mask = .true.
       ! end where
       ! write (*,*) 'count(DC=1):',count(mask)
       ! deallocate (mask)
       
       if (my_nparts > 0) then
          CALL mix (status, my_nparts,  &
               LAT_SHUFFLED_TEMP, LAT_OLD_MIX_SHUFFLED_TEMP, LAT_SHUFFLED, &
               LON_SHUFFLED_TEMP, LON_OLD_MIX_SHUFFLED_TEMP, LON_SHUFFLED, &
               ZETA_SHUFFLED_TEMP, ZETA_SHUFFLED, &
               STATE_SHUFFLED, STATE_VERT_SHUFFLED, &
               MIXSPECARR_SHUFFLED_TEMP, MIXSPECARR_SHUFFLED, &
               THETA_SHUFFLED_TEMP, THETA_OLD_MIX_SHUFFLED_TEMP, &
               BVF_SHUFFLED_TEMP, BVF_OLD_MIX_SHUFFLED_TEMP)
          IF (status /= 0) CALL error_bi('Error in MIX !!!',substr)
       else
          LAT_SHUFFLED  = LAT_SHUFFLED_TEMP
          LON_SHUFFLED  = LON_SHUFFLED_TEMP
          ZETA_SHUFFLED = ZETA_SHUFFLED_TEMP
          do i = 1, nspec
             MIXSPECARR_SHUFFLED(i)%values = MIXSPECARR_SHUFFLED_TEMP(i)%values
          enddo
          if (switch_mixing==2) then
             THETA_SHUFFLED = THETA_SHUFFLED_TEMP
             BVF_SHUFFLED = BVF_SHUFFLED_TEMP
          endif
       endif

       dnparts = my_nparts
       dnparts_co(1) = dnparts
       nparts = p_sum(dnparts)
       if (p_pe==0) write (*,*) 'nach mix: rank, dnparts, nparts =',p_pe,dnparts,nparts

!!!!!  ????
       if (dnparts > dnparts_max_shuffle) &
            CALL error_bi &
            ('Number of particles after mixing greater than dnparts_max_shuffle !!! ', &
            substr)

       IF (.NOT. lbmixevent) THEN
          
          !write (*,*) 'call clamsmix_reshuffle', p_pe
          call clamsmix_reshuffle

          dnparts_co(1) = dnparts
          nparts = p_sum(dnparts)
          write (*,*) 'nach reshuffle: rank, dnparts, nparts =',p_pe,dnparts,nparts
          
          !--------------------------------------------------------------------
          ! Interpolate TEMP, PRESS (and other traj parameters)
          !--------------------------------------------------------------------
          
          ! current time
          itime%year    = YEAR
          itime%month   = MONTH
          itime%day     = DAY
          itime%sec     = HOUR*3600+MINUTE*60+SECOND+delta_time
          
          ! Calculate time of previous windfile
          ipasttime%sec = fut_sec
          ipasttime%day = fut_day
          ipasttime%month = fut_month
          ipasttime%year = fut_year
          
          CALL incdat_interface( ipasttime%sec,ipasttime%day,ipasttime%month,ipasttime%year, &
               -irdsec,-irdday,0,0, dates30)     
          
          DO i = 1, nparams
             IF (paramnames(i)=='H2O' ) CYCLE
             IF (paramnames(i)=='IWC' ) CYCLE
             IF (paramnames(i)=='CLWC') CYCLE
             PARAM(i)%values = mdi
             CALL interpolate_param(LON, LAT, ZETA, itime, ipasttime, &
                  i, paramnames(i), PARAM(i)%values, lcoupled)
          END DO
          
          !--------------------------------------------------------------------
          ! Set _OLD arrays 
          !--------------------------------------------------------------------

          ! if MIX was executed in current timestep:
          IF (lmixevent) THEN
             LON_OLD_MIX = LON
             LAT_OLD_MIX = LAT
             ! for vertical mixing: THETA_OLD_MIX and BVF_OLD_MIX
             if (switch_mixing==2) then
                THETA_OLD_MIX  = PARAM(str_pos(nparams,paramnames,'THETA'))%values
                if (uppercase(vert_mix_param)=='WET') then
                   BVF_OLD_MIX = PARAM(str_pos(nparams,paramnames,'BVF_WET'))%values
                else
                   BVF_OLD_MIX = PARAM(str_pos(nparams,paramnames,'BVF'))%values
                endif
             endif
          END IF
          
          LAT_OLD = LAT
          LON_OLD = LON
          LEV_OLD = ZETA
          JULTIME_OLD = JULTIME
          
          ! TEMP_OLD and PRESS_OLD
          DO i = 1, nparams_old
             PARAM_OLD(i)%values = PARAM(str_pos(nparams,paramnames,paramnames_old(i)))%values
             !write (*,*) 'iparam_old, iparam:', i, str_pos(nparams,paramnames,paramnames_old(i))
          END DO

       ENDIF

       IF (lmixevent) THEN

          asad_gfirst = .TRUE.

!!!!! ???
          LON_OLD_MIX_SHUFFLED = LON_SHUFFLED
          LAT_OLD_MIX_SHUFFLED = LAT_SHUFFLED


          IF (p_pe==0) write(*,*) 'Old grid_switch: ', grid_switch 
          IF (grid_switch == 0) THEN
             grid_switch = 1
          ELSE
             grid_switch = 0
          END IF
          grid_switch_co = grid_switch
          IF (p_pe==0) write(*,*) 'New grid_switch: ', grid_switch 

!!!!!
!!$          IF (loutput_mix) THEN
!!$             CALL set_channel_output(status, 'clams', .FALSE.)
!!$             CALL channel_halt(substr, status)
!!$             CALL set_channel_output(status, 'clamsmix', .TRUE.)
!!$             CALL channel_halt(substr, status)
!!$             CALL messy_write_output
!!$             if (lclamsoutevent) then
!!$                CALL set_channel_output(status, 'clams', .TRUE.)
!!$                CALL channel_halt(substr, status)
!!$             endif
!!$          ENDIF

       END IF

       ! Deallocate arrays:
       DO i= 1, nspec
          DEALLOCATE(MIXSPECARR_SHUFFLED_TEMP(i)%values)
       END DO
       DEALLOCATE(MIXSPECARR_SHUFFLED_TEMP)
       DEALLOCATE(LAT_SHUFFLED_TEMP)
       DEALLOCATE(LON_SHUFFLED_TEMP)
       DEALLOCATE(ZETA_SHUFFLED_TEMP)
       DEALLOCATE(LAT_OLD_MIX_SHUFFLED_TEMP)
       DEALLOCATE(LON_OLD_MIX_SHUFFLED_TEMP)
       DEALLOCATE(THETA_SHUFFLED_TEMP)
       DEALLOCATE(BVF_SHUFFLED_TEMP)
       DEALLOCATE(THETA_OLD_MIX_SHUFFLED_TEMP)
       DEALLOCATE(BVF_OLD_MIX_SHUFFLED_TEMP)

       DEALLOCATE(indarr)
       DEALLOCATE(irange)
       DEALLOCATE(dnparts_level, nparts_level)
       DEALLOCATE(lev_min_act, lev_max_act, lev_delta_act)
       
       DEALLOCATE(adapt_par%fac_min)
       DEALLOCATE(adapt_par%fac_max)
       DEALLOCATE(adapt_par%r_min_c)
       DEALLOCATE(adapt_par%r_max_c)
       DEALLOCATE(adapt_par%r_min_h)
       DEALLOCATE(adapt_par%r_max_h)
       DEALLOCATE(adapt_par%r_lim_c_outside)
       DEALLOCATE(adapt_par%r_lim_h_outside)
       DEALLOCATE(adapt_par%r_lim_c_inside)
       DEALLOCATE(adapt_par%r_lim_h_inside)

       ! Wird in BMIX genutzt !
       IF (.NOT. USE_CLAMSBMIX) THEN
          DEALLOCATE(levelrange)
          DEALLOCATE(l_min_act, l_max_act, l_delta_act)
       ENDIF

    ENDIF

  END SUBROUTINE clamsmix_global_end

!-----------------------------------------------------------------------

  SUBROUTINE clamsmix_free_memory

    USE messy_clamsmix_global, ONLY: adapt_par

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsmix_free_memory'

    !WRITE(*,*) substr

    IF (ASSOCIATED(adapt_par%lev_grid)) &
         DEALLOCATE(adapt_par%lev_grid,adapt_par%lev_delta)
    IF (ASSOCIATED(adapt_par%r_grid)) &
         DEALLOCATE(adapt_par%r_grid)

    DEALLOCATE(MIXSPECARR)
    DEALLOCATE(MIXSPECARR_SHUFFLED)

  END SUBROUTINE clamsmix_free_memory

!-----------------------------------------------------------------------

  !****************************************************************************
  ! 
  !****************************************************************************
  SUBROUTINE clamsmix_error_handler (status, substr)

    USE messy_main_mpi_bi,     ONLY: p_pe
    USE messy_main_blather_bi, ONLY: error_bi

    IMPLICIT NONE

    INTEGER      :: status
    CHARACTER(*) :: substr

    CHARACTER(80) :: errorstr

    SELECT CASE (status)
    CASE (201)
       errorstr = 'Error in clamsmix_shuffle: size of field too small => increase rres_shuffle !!'
    END SELECT

    CALL error_bi(errorstr,substr)


  END SUBROUTINe clamsmix_error_handler

  !****************************************************************************
  ! 
  !****************************************************************************
  SUBROUTINE clamsmix_sort (status, level, indarr, irange)

    USE messy_clams_global,      ONLY: prec, ntasks, rank, mdi, eps, dnparts
    USE messy_clamsmix_global,   ONLY: l_min_act, l_max_act, l_nlevs
    USE messy_clams_tools_utils, ONLY: bubble_sort_index, quick_sort_index
    USE messy_main_mpi_bi,       ONLY: p_bcast, p_pe

    IMPLICIT NONE

    integer,                     intent(out) :: status
    ! level    : THETA/ZETA-Values 
    ! levelind : Index of THETA/ZETA-Values in THETA/ZETA-Grid
    ! indarr   : Sorted index of levelind
    ! irange(1,i,k) :: first point on level i in levelind(indarr) on rank k
    ! irange(2,i,k) :: last point on level i in levelind(indarr) on rank k
    real(prec), dimension(:),    intent(in)  :: level    ! (1:dnparts)
    integer,    dimension(:),    pointer     :: levelind ! (1:dnparts)
    integer,    dimension(:),    pointer     :: indarr   ! (1:dnparts)
    integer,    dimension(:,:,:),pointer     :: irange   ! (2,nlev,0:ntasks-1)
 
    integer :: npoints

    integer :: ilev, i, start


    status = 0

    npoints = size(level)

    ! allocate arrays
    allocate (indarr(npoints))
    allocate (levelind(npoints))
    allocate (irange(2,l_nlevs,0:ntasks-1))

    ! initialize arrays
    indarr   = -1
    levelind = -1
    irange   = -1

    ! for all theta/zeta-values: get level index
    do ilev = 1, l_nlevs
       where (l_min_act(ilev)<=level .and. level<l_max_act(ilev)) 
          levelind = ilev
       end where
    enddo
!!!!! ACHTUNG folgendes gilt auch fuer missing values
    where (level<l_min_act(1))
       levelind = 1
    endwhere
    where (level>=l_max_act(l_nlevs))
       levelind = l_nlevs
    endwhere
!!!!! welchen levelind fuer mdi setzen ???
    where (abs((level-mdi)/mdi)<eps)
       levelind = 9999
    endwhere

    ! get index of sorted array 
!     call bubble_sort_index (levelind,indarr,size(levelind))
    call quick_sort_index (levelind,indarr)

    ! get first and last point on level ilev in levelind(indarr)

    if (levelind(indarr(1))/=9999) then
       start = 1
       ilev = levelind(indarr(1))
       do i = 2, npoints
          if (levelind(indarr(i))>ilev) then
             irange(1,ilev,rank) = start
             irange(2,ilev,rank) = i-1
             ilev = levelind(indarr(i))
             start = i
             if (ilev==9999) exit
             if (i==npoints) then
                irange(1,ilev,rank) = i
                irange(2,ilev,rank) = i
             endif
          elseif (i==npoints) then
             irange(1,ilev,rank) = start
             irange(2,ilev,rank) = npoints
          endif
       enddo
    endif       
     
    ! clean up
    deallocate (levelind)


!!$    ! Broadcast irange
!!$    do i = 0, ntasks-1
!!$       call p_bcast (irange(:,:,i),i)
!!$    enddo

!!!!! status bei Fehlern setzen !?!


  END SUBROUTINE clamsmix_sort


  !****************************************************************************
  ! Shuffle all positions and parameters:
  ! Each rank gets all points between levelrange(1,irank) and levelrange(2,irank)
  !****************************************************************************
  SUBROUTINE clamsmix_shuffle (status, my_nparts, indarr, irange, levelrange, &
                        dnparts_level, nparts_level, level, infield, outfield)

    USE messy_clams_global,    ONLY: rank, ntasks, prec, mdi, dnparts_max_shuffle
    USE messy_clamsmix_global, ONLY: l_nlevs, l_min_act, l_max_act, &
                                     lev_min_act, lev_max_act, nintervals
    USE messy_main_mpi_bi,     ONLY: p_recv, p_send

    IMPLICIT NONE

    integer,                    intent(out) :: status
    integer,                    intent(out) :: my_nparts
    integer,    dimension(:),    pointer    :: indarr        ! (1:dnparts)
    integer,    dimension(:,:,:),pointer    :: irange        ! (2,nlev,ntasks)
    integer,    dimension(:,:),  pointer    :: levelrange    ! (2,0:ntasks-1)
    integer,    dimension(:,:),  pointer    :: dnparts_level ! (nlev,0:ntasks-1)
    integer,    dimension(:),    pointer    :: nparts_level  ! (nlev)
    REAL(PREC), DIMENSION(:),    intent(in) :: level
    REAL(PREC), DIMENSION(:),    intent(in) :: infield
    real(prec), dimension(:),    pointer    :: outfield

    real(prec), dimension(:),   allocatable :: sendarr
    logical,    dimension(:),   allocatable :: mask

    INTEGER :: recvrank, sendrank, anz,  &
               startpos, startind, endind, ind_irange

    integer, save :: tag = 100.

    status = 0

    if (levelrange(1,rank)==-1) then  
       my_nparts = 0
    else
       my_nparts = SUM (nparts_level(levelrange(1,rank):levelrange(2,rank)))
    end if
    if (my_nparts > dnparts_max_shuffle) then
       status = 201
       return
    endif

    startpos = 1
   
    allocate (outfield(dnparts_max_shuffle))
    
    outfield = mdi

    tag = tag + ntasks

    do recvrank = 0, ntasks-1

       ! rank recvrank empfaengt von allen anderen, alle anderen senden

       if (rank==recvrank) then   ! receive

          !!!!! eigene Daten zuerst auf outfield schreiben

          if (my_nparts>0 .AND. levelrange(1,rank)/=-1) then
             
             ind_irange = levelrange(1,rank)
             do while (irange(1,ind_irange,rank)==-1 .and. ind_irange<levelrange(2,rank))
                ind_irange = ind_irange+1
             end do
             startind = irange(1,ind_irange,rank)
             
             ind_irange = levelrange(2,rank)
             do while (irange(2,ind_irange,rank)==-1 .and. ind_irange>levelrange(1,rank))
                ind_irange = ind_irange-1
             enddo
             endind = irange(2,ind_irange,rank)
             
! op_sb_20200204
!!$         if (endind>=startind) then
            if (endind>=startind .and. startind/=-1) then

                anz = SUM(dnparts_level(levelrange(1,rank):levelrange(2,rank),rank))
                if (startpos+anz-1>size(outfield)) then
                   status = 201
                   return
                endif
                outfield(startpos:startpos+anz-1) = infield(indarr(startind:endind))
                startpos = startpos + anz
             endif
             
          endif
          
          do sendrank = 0, ntasks-1
             
             if (sendrank/=recvrank) then
!!!!! hinten an outfield anhaengen
                if (levelrange(1,rank)/=-1) then 
                   
                   anz = SUM(dnparts_level(levelrange(1,rank):levelrange(2,rank),sendrank))
                   ! write (*,*) 'vor receive: recvrank, sendrank, anz=',rank,sendrank,anz
                   
                   if (anz > 0) then
                      
                      if (startpos+anz-1>size(outfield)) then
                         status = 201
                         return
                      endif
                      
                      !write (*,*) 'receive von rank: recvrank, sendrank: ',rank, sendrank, anz
                      call p_recv (outfield(startpos:startpos+anz-1),sendrank,tag+sendrank)  
                      
                      startpos = startpos + anz
                      
                   endif
                   
                endif
             endif
          enddo

       else   ! send
          
          if (levelrange(1,recvrank)/=-1) then

             anz = SUM(dnparts_level(levelrange(1,recvrank):levelrange(2,recvrank),rank))
             
             if (anz > 0) then
                
                ind_irange = levelrange(1,recvrank)
                do while (irange(1,ind_irange,rank)==-1 .and. ind_irange<levelrange(2,recvrank))
                   ind_irange = ind_irange+1
                end do
                startind = irange(1,ind_irange,rank)
                
                ind_irange = levelrange(2,recvrank)
                do while (irange(2,ind_irange,rank)==-1 .and. ind_irange>levelrange(1,recvrank))
                   ind_irange = ind_irange-1
                enddo
                endind = irange(2,ind_irange,rank)
                
                if (endind>=startind) then
                   
                   allocate (sendarr(anz))
                   sendarr = infield(indarr(startind:endind)) 
                   
                   ! write (*,*) 'send von rank: recvrank, sendrank: ',recvrank, rank, anz
                   call p_send (sendarr,recvrank,tag+rank)
                   
                   deallocate (sendarr)
                   
                endif
                
             endif
             
          endif
          
       endif

    enddo

!!!!! Grenzbereiche

    if (nintervals>1) then

       tag = tag + ntasks
       
       do recvrank = 0, ntasks-1
          
          if (rank == recvrank) then
             
             ! eigene Werte aus Grenzbereichen an outfield anhaengen
             
             ! lev_min_act(levelrange(1,rank))<=val<l_min_act(levelrange(1,rank))
             ! l_max_act(levelrange(2,rank))<val<=lev_max_act(levelrange(2,rank))
             
             allocate (mask(size(level)))
             mask = .false.
             where ((lev_min_act(levelrange(1,rank))<=level .and. &
                  level<l_min_act(levelrange(1,rank))) .or. &
                  (l_max_act(levelrange(2,rank))<level .and. &
                  level<=lev_max_act(levelrange(2,rank))))
                mask = .true.
             end where
             if (count(mask) > 0) then
                if (startpos+count(mask)-1>size(outfield)) then
                   status = 201
                   return
                endif
                outfield(startpos:startpos+count(mask)-1) = pack(infield,mask)
                startpos = startpos+count(mask)
             endif
             deallocate (mask)
             
             
             ! Werte von allen anderen CPUs aus den Grenzbereichen empfangen
             do sendrank = 0, ntasks-1
                if (sendrank /= recvrank) then
                   
                   ! empfange count
                   call p_recv (anz,sendrank,tag+sendrank)
                   
                   ! wenn count>0: empfange sendarr
                   if (anz>0) then
                      
                      if (startpos+anz-1>size(outfield)) then
                         status = 201
                         return
                      endif
                      
                      call p_recv (outfield(startpos:startpos+anz-1),sendrank,tag+sendrank)  
                       
                      startpos = startpos + anz
                      
                   endif
                   
                endif
             enddo
             
          else
             
             ! Sende Werte aus Grenzbereichen an recvrank
             
             allocate (mask(size(level)))
             mask = .false.
             where ((lev_min_act(levelrange(1,recvrank))<=level .and. &
                  level<l_min_act(levelrange(1,recvrank))) .or. &
                  (l_max_act(levelrange(2,recvrank))<level .and. &
                  level<=lev_max_act(levelrange(2,recvrank))))
                mask = .true.
             end where
             
             ! sende count(mask) !?!
             call p_send (count(mask),recvrank,tag+rank)
             
             if (count(mask) > 0) then
                allocate (sendarr(count(mask)))
                sendarr = pack(infield,mask)
                
                ! sende sendarr
                call p_send (sendarr,recvrank,tag+rank)
                
                deallocate (sendarr)
             endif
             deallocate (mask)
             
          endif

       enddo

    endif

    my_nparts = startpos - 1

  END SUBROUTINE clamsmix_shuffle

  !****************************************************************************
  ! 
  !****************************************************************************
  SUBROUTINE clamsmix_reshuffle
!DEC$ NOOPTIMIZE

    USE messy_main_mpi_bi,     ONLY: p_sum, p_bcast, p_recv, p_send, p_barrier
    USE messy_main_blather_bi, ONLY: error_bi

    USE messy_clams_global,    ONLY: rank, ntasks, prec, mdi, nspec, &
                                     nparts, nparts_max, dnparts, dnparts_max
    USE messy_clamsmix_global, ONLY: l_nlevs, switch_mixing

    IMPLICIT NONE

    real(prec), parameter :: maxdev = 0.2  !!! Wie gross darf die Abweichung sein?
    integer,    parameter :: minpoints = 1000 !!! Bei welcher Anz. lohnt der Aufwand?

    real(prec), dimension(:), allocatable :: sendarr
    integer,    dimension(:), allocatable :: nparts_pe

    integer :: dcounter, counter, npoints
    integer :: irank, sendrank, recvrank, i, ispec

    integer :: tag = 100.

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsmix_reshuffle'

    ! Anzahl Punkte auf aktuellem rank: dnparts
    ! Anzahl Punkte insgesamt:          nparts
    ! Anzahl Punkte nach Gleichverteilung pro rank: ~ nparts/ntasks (+ntasks-1)

    !write (*,*) rank, 'in clamsmix_reshuffle'


!!!!! ACHTUNG: kein Reshuffle fuer THETA und BVF noetig, da diese nur 
!!!!!          als Input fuer MIX genutzt und nicht veraendert werden !!!

    ! Sicherstellen, dass Anzahl Punkte nicht zu gross geworden ist:
    if (nparts > nparts_max) then
       CALL error_bi &
            ('nparts > nparts_max => increase rres in clams.nml !!!',substr)
    endif
!!$    write (*,*)'nparts, ntasks, nparts/ntasks, dnparts_max=', &
!!$         nparts, ntasks, nparts/ntasks, dnparts_max
!!$    write (*,*) 'nparts - ((nparts/ntasks)*(ntasks-1))',nparts - ((nparts/ntasks)*(ntasks-1))

    if (nparts/ntasks > dnparts_max) then
       CALL error_bi &
            ('nparts/ntasks > dnparts_max => increase rres in clams.nml !!!',substr)
    elseif (nparts - ((nparts/ntasks)*(ntasks-1)) > dnparts_max) then
       CALL error_bi &
            ('dnparts > dnparts_max => increase rres in clams.nml !!!',substr)
    endif

    ! Ueberpruefe Abweichung von der Idealverteilung fuer alle PEs 
    ! Sicherstellen, dass dnparts <= dnparts_max !
    dcounter = 0
    !write (*,*) 'rank, nparts, minpoints, maxdev', rank, nparts, minpoints, maxdev
    if (nparts > ntasks*minpoints) then  
       !write (*,*) rank, 'nparts/ntasks,dnparts,dnparts_max,Abweichung:', &
       !     real(nparts/ntasks), dnparts, dnparts_max, &
       !     abs(real(nparts)/ntasks-dnparts)/real(nparts/ntasks)
       if (abs(real(nparts)/ntasks-dnparts)/real(nparts/ntasks) >= maxdev) &
            dcounter = 1
       if (dnparts > dnparts_max) dcounter = 1
    endif
    counter = p_sum(dcounter)
    if (rank==0) write (*,*) 'counter:', counter

    ! Mindestens eine PE hat viel zu wenige oder viel zu viele Punkte:
    if (counter > 0) then

       ! Alle PEs kennen die Anz. Punkte der anderen
       allocate (nparts_pe(0:ntasks-1))
       nparts_pe = 0
       nparts_pe(rank) = dnparts
       do irank = 0, ntasks-1
          call p_bcast (nparts_pe(irank),irank)
       enddo
       !write (50+rank,*) nparts_pe 
       if (rank==0) write (*,*) 'nparts_pe vor reshuffle:', nparts_pe


       tag = 100.

       ! Die Anzahl Punkte fuer PE irank soll wieder nparts/ntasks sein:
       do irank = 0, ntasks-2

          ! zuviele Punkte auf aktueller PE: verschiebe Punkte auf andere PEs
          if (nparts_pe(irank) > nparts/ntasks) then
             
!!!!! Wenn folgende Ausgabe fehlt, kann es bei Compilation mit -O2 zu einem Absturz kommen ?!?        
!             write (*,*) rank, 'Zuviele Punkte auf rank',irank

             recvrank = irank + 1
             
             do while (nparts_pe(irank)>nparts/ntasks .and. recvrank<ntasks)
                
                ! Wenn Anz. Punkte auf recvrank zu klein: Punkte dorthin schieben
                if (nparts_pe(recvrank)<nparts/ntasks) then
                   npoints = min(nparts/ntasks-nparts_pe(recvrank),nparts_pe(irank)-nparts/ntasks)
                   !if (rank==0) write (*,*) 'irank, recvrank, npoints:', irank, recvrank, npoints

                   if (rank==irank) then
                      ! send
                      !write (*,*) 'send von rank: recvrank, sendrank: ',recvrank, rank, npoints
                      allocate (sendarr(npoints))
                      sendarr = LAT_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,recvrank,tag+1)
                      sendarr = LON_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,recvrank,tag+2)
                      sendarr = ZETA_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,recvrank,tag+3)
                      sendarr = LAT_OLD_MIX_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,recvrank,tag+4)
                      sendarr = LON_OLD_MIX_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,recvrank,tag+5)
                      sendarr = STATE_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,recvrank,tag+6)
                      sendarr = STATE_VERT_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,recvrank,tag+7)
                      do ispec = 1, nspec
                         sendarr = MIXSPECARR_SHUFFLED(ispec)%values(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                         call p_send (sendarr,recvrank,tag+10+ispec)
                      enddo
                      deallocate (sendarr)

                   elseif (rank==recvrank) then
                      ! receive
                      !write (*,*) 'receive von rank: recvrank, sendrank: ',rank, irank, npoints
                      call p_recv (LAT_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   irank,tag+1)  
                      call p_recv (LON_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   irank,tag+2)  
                      call p_recv (ZETA_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   irank,tag+3)  
                      call p_recv (LAT_OLD_MIX_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   irank,tag+4)  
                      call p_recv (LON_OLD_MIX_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   irank,tag+5)  
                      call p_recv (STATE_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   irank,tag+6)  
                      call p_recv (STATE_VERT_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   irank,tag+7)  
                      do ispec = 1, nspec
                         call p_recv (MIXSPECARR_SHUFFLED(ispec)%values &
                                        (nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                      irank,tag+10+ispec)  
                      enddo
                   endif
                   
                   tag = tag + 10 + nspec

                   !call p_barrier
                   nparts_pe(irank) = nparts_pe(irank) - npoints
                   nparts_pe(recvrank) = nparts_pe(recvrank) + npoints
                   !if (rank==0) write (*,*) 'nparts_pe:', nparts_pe

                endif
                
                recvrank = recvrank + 1

             enddo
              
            
          ! zuwenige Punkte auf aktueller PE: hole Punkte von anderen PEs
          elseif (nparts_pe(irank) < nparts/ntasks) then
                
!!!!! Wenn folgende Ausgabe fehlt, kann es bei Compilation mit -O2 zu einem Absturz kommen ?!?        
!             write (*,*) 'Zuwenige Punkte auf rank',irank

             sendrank = irank + 1
             
             do while (nparts_pe(irank)<nparts/ntasks .and. sendrank<ntasks)
                
                ! Wenn Anz. Punkte auf sendrank zu gross: Punkte von sendrank holen
                if (nparts_pe(sendrank)>nparts/ntasks) then
                   npoints = min(nparts_pe(sendrank)-nparts/ntasks,nparts/ntasks-nparts_pe(irank))
                   !if (rank==0) write (*,*) 'irank, sendrank, npoints:', irank, sendrank, npoints

                   if (rank==irank) then
                      ! receive
                      !write (*,*) 'receive: recvrank, sendrank: ',rank, sendrank, npoints
                      call p_recv (LAT_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   sendrank,tag+1)  
                      call p_recv (LON_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   sendrank,tag+2)  
                      call p_recv (ZETA_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   sendrank,tag+3)  
                      call p_recv (LAT_OLD_MIX_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   sendrank,tag+4)  
                      call p_recv (LON_OLD_MIX_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   sendrank,tag+5)  
                      call p_recv (STATE_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   sendrank,tag+6)  
                       call p_recv (STATE_VERT_SHUFFLED(nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                   sendrank,tag+7)  
                     do ispec = 1, nspec
                         call p_recv (MIXSPECARR_SHUFFLED(ispec)%values &
                                        (nparts_pe(rank)+1:nparts_pe(rank)+npoints), &
                                      sendrank,tag+10+ispec)  
                      enddo

                   elseif (rank==sendrank) then
                      ! send
                      !write (*,*) 'send: recvrank, sendrank: ',irank, sendrank, npoints
                      allocate (sendarr(npoints))
                      sendarr = LAT_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,irank,tag+1)
                      sendarr = LON_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,irank,tag+2)
                      sendarr = ZETA_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,irank,tag+3)
                      sendarr = LAT_OLD_MIX_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,irank,tag+4)
                      sendarr = LON_OLD_MIX_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,irank,tag+5)
                      sendarr = STATE_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,irank,tag+6)
                      sendarr = STATE_VERT_SHUFFLED(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                      call p_send (sendarr,irank,tag+7)
                      do ispec = 1, nspec
                         sendarr = MIXSPECARR_SHUFFLED(ispec)%values(nparts_pe(rank)-npoints+1:nparts_pe(rank)) 
                         call p_send (sendarr,irank,tag+10+ispec)
                      enddo
                      deallocate (sendarr)

                   endif
 
                   tag = tag + 10 + nspec

                   !call p_barrier 
                   nparts_pe(irank) = nparts_pe(irank) + npoints
                   nparts_pe(sendrank) = nparts_pe(sendrank) - npoints
                   !write (*,*) rank, 'nparts_pe:', nparts_pe

                endif
                    
                sendrank = sendrank + 1

             enddo
                          
          endif
  
       enddo

       if (rank==0) write (*,*) 'nparts_pe nach reshuffle:', nparts_pe
       
       dnparts = nparts_pe(rank)

       deallocate (nparts_pe)

    endif

    write (*,*) 'vor barrier',rank
    call p_barrier 
    
    ! Write from SHUFFLED to non-SHUFFLED arrays(max.nparts/ntasks)
    IF (dnparts > 0) THEN
       LAT(1:dnparts)         = LAT_SHUFFLED(1:dnparts)
       LON(1:dnparts)         = LON_SHUFFLED(1:dnparts)
       ZETA(1:dnparts)        = ZETA_SHUFFLED(1:dnparts)
       LAT_OLD_MIX(1:dnparts) = LAT_OLD_MIX_SHUFFLED(1:dnparts)
       LON_OLD_MIX(1:dnparts) = LON_OLD_MIX_SHUFFLED(1:dnparts)
       STATE(1:dnparts)       = STATE_SHUFFLED(1:dnparts)
       STATE_VERT(1:dnparts)  = STATE_VERT_SHUFFLED(1:dnparts)
       DO ispec = 1, nspec
          MIXSPECARR(ispec)%values(1:dnparts) = MIXSPECARR_SHUFFLED(ispec)%values(1:dnparts)
       END DO
    END IF
    
    LAT(dnparts+1:dnparts_max)         = mdi
    LON(dnparts+1:dnparts_max)         = mdi
    ZETA(dnparts+1:dnparts_max)        = mdi
    LAT_OLD_MIX(dnparts+1:dnparts_max) = mdi
    LON_OLD_MIX(dnparts+1:dnparts_max) = mdi
    STATE(dnparts+1:dnparts_max)       = 0
    STATE_VERT(dnparts+1:dnparts_max)  = 0
    do ispec = 1, nspec
       MIXSPECARR(ispec)%values(dnparts+1:dnparts_max) = mdi
    enddo
   

  END SUBROUTINE clamsmix_reshuffle
!-----------------------------------------------------------------------

#endif
END MODULE messy_clamsmix_si
