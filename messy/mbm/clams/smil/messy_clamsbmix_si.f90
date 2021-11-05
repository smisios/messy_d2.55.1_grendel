!**********************************************************************
MODULE messy_clamsbmix_si

#if defined(ECHAM5) || defined(MBM_CLAMS)
!**********************************************************************
!  Submodel interface for clamsbmix 
!**********************************************************************

  USE messy_clamsbmix

  USE messy_main_timer_event, ONLY: time_event, io_time_event, event_is_active 

  USE messy_clams_global,     ONLY: prec, dp, species_type, species_type_3d, &
                                    paramtype

  IMPLICIT NONE
  PRIVATE

! op_pj_20160606+ 
!!$  TYPE(time_event),    PUBLIC   :: bmixevent
!!$  TYPE(io_time_event), PUBLIC:: io_bmixevent
  TYPE(time_event),    PUBLIC, SAVE :: bmixevent
  TYPE(io_time_event), PUBLIC, SAVE :: io_bmixevent
! op_pj_20160606-
  TYPE(species_type_3d), DIMENSION(:), POINTER , PUBLIC, SAVE :: E5CHEMSPECARR_tte => NULL() ! ECHAMTRACER Tendency

  PUBLIC :: clamsbmix_initialize
  PUBLIC :: clamsbmix_init_memory
  PUBLIC :: clamsbmix_init_coupling
  PUBLIC :: clamsbmix_global_end
  PUBLIC :: clamsbmix_free_memory

  ! MODULE VARIABLES

  REAL(PREC), DIMENSION(:), POINTER :: LAT     => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON     => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LEV     => NULL()
  REAL(DP),   POINTER               :: JULTIME
  REAL(PREC), DIMENSION(:), POINTER :: LAT_OLD => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON_OLD => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LEV_OLD => NULL()
  REAL(DP),   POINTER               :: JULTIME_OLD
  REAL(PREC), DIMENSION(:), POINTER :: LAT_OLD_MIX => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON_OLD_MIX => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LAT_LOW  => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON_LOW  => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LEV_LOW  => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LAT_UPPER => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON_UPPER => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LEV_UPPER => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: STATE     => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: THETA_OLD_MIX  => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: BVF_OLD_MIX => NULL()
  TYPE(paramtype),    DIMENSION(:), POINTER :: PARAM     ! from CLAMS
  TYPE(paramtype),    DIMENSION(:), POINTER :: PARAM_OLD ! from CLAMS
  TYPE(species_type), DIMENSION(:), POINTER :: MIXSPECARR      ! from MIX
  TYPE(species_type), DIMENSION(:), POINTER :: SPECARR_LOW 
  TYPE(species_type), DIMENSION(:), POINTER :: SPECARR_UPPER
  ! Shuffle arrays:
  REAL(PREC),         DIMENSION(:), POINTER :: LEV_SHUFFLED        => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LAT_SHUFFLED         => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: LON_SHUFFLED         => NULL()
!!!!! ???
!  REAL(PREC),         DIMENSION(:), POINTER :: LAT_OLD_MIX_SHUFFLED => NULL()
!  REAL(PREC),         DIMENSION(:), POINTER :: LON_OLD_MIX_SHUFFLED => NULL()
  REAL(PREC),         DIMENSION(:), POINTER :: STATE_SHUFFLED       => NULL()
  TYPE(species_type), DIMENSION(:), POINTER :: MIXSPECARR_SHUFFLED 

#ifdef ECHAM5
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_LAT3D_D    => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_LON3D_D    => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_ZETA_D     => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_PRESS_D    => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_SH_D       => NULL()
  REAL(DP), DIMENSION(:),     POINTER :: E5_LAT_1D     => NULL()
  REAL(DP), DIMENSION(:),     POINTER :: E5_LON_1D     => NULL()
  REAL(DP), DIMENSION(:),     POINTER :: E5_ZETA_1D    => NULL()
  REAL(DP), DIMENSION(:),     POINTER :: E5_PRESS_1D   => NULL()
  REAL(DP), DIMENSION(:),     POINTER :: E5_TRACER_1D  => NULL()
#endif
  ! op_pj_20170110+
  LOGICAL, SAVE :: L_CLAMSCHEME5 = .FALSE.
  ! op_pj_20170110-

!--------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsbmix_initialize

    ! BMIL
    USE messy_main_blather_bi,   ONLY: error_bi
    USE messy_main_timer_bi,     ONLY: timer_event_init

    ! SMCL
    USE messy_main_tools,        ONLY: find_next_free_unit
    USE messy_main_timer,        ONLY: delta_time
!!$ USE messy_main_switch,       ONLY: USE_CLAMSCHEME5 ! op_pj_20170110

    USE messy_clams_global,      ONLY: initfile, nspec, species_type, SPECARR, &
                                       rank, nchemspec, H2O_index, H2O_100_index, &
                                       IWC_index, IWC_100_index, CLWC_index, &
                                       maxspec, nparams, paramnames, &
                                       nparams_old, paramnames_old
    USE messy_clamsbmix_global,  ONLY: timestep_bmix

!!#D clamschem +
    USE messy_clamschem_defs_mod, ONLY: chch_defs 
!!#D clamschem +
    USE messy_clams_tools_utils, ONLY: uppercase
    USE netcdf

    USE messy_main_mpi_bi, ONLY: p_pe

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsbmix_initialize'
    INTEGER :: status, iou, i, n, nodd, ichem, ios
    INTEGER :: nhelpspec
    CHARACTER(2)  :: ctype
    CHARACTER(20) :: specname

    INTEGER        :: rcode, ncid, varid
    CHARACTER(150) :: ncdf_file                      
    INTEGER,DIMENSION(NF90_MAX_VAR_DIMS) :: dim_array

    IF (p_pe==0) THEN
       WRITE(*,*)
       write(*,*) uppercase(substr)
    ENDIF

    ! Read namelist variables:
    iou = find_next_free_unit(100,200)

    ! Read namelist and set default values:
    !   (he list of boundfiles are read too)
    CALL clamsbmix_read_nml(status, iou)
    IF (status /= 0) CALL error_bi('Error in clamsbmix_read_nml ',substr)

    if (mod(timestep_bmix*3600,int(delta_time)) /= 0) then
       call error_bi ("Wrong bmix timestep !!!",substr)
    endif

    ! Define BMIX event:
    io_bmixevent%counter = timestep_bmix 
    io_bmixevent%unit = 'hours'
    io_bmixevent%adjustment = 'exact'
    io_bmixevent%offset = -delta_time
    CALL timer_event_init (bmixevent, io_bmixevent, 'BMIX_Event', 'present')

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

    DO i=1,nspec
       IF (TRIM(specarr(i)%name)=='H2O') H2O_index = i
       IF (TRIM(specarr(i)%name)=='H2O_100') H2O_100_index = i
       IF (TRIM(specarr(i)%name)=='IWC') IWC_index = i
       IF (TRIM(specarr(i)%name)=='IWC_100') IWC_100_index = i
       IF (TRIM(specarr(i)%name)=='CLWC') CLWC_index = i
    END DO
! op_sb_20191021 (moved to init_coupling)
!!$    ALLOCATE (MIXSPECARR(nspec))
!!$    ALLOCATE (MIXSPECARR_SHUFFLED(nspec))

! op_pj_20170110+
!!$    IF (USE_CLAMSCHEME5) THEN
!!$       ALLOCATE (E5CHEMSPECARR_tte(nchemspec))
!!$       DO i = 1, nchemspec
!!$          E5CHEMSPECARR_tte(i)%name = chch_defs(i)%speci
!!$          E5CHEMSPECARR_tte(i)%ctype = chch_defs(i)%ctype
!!$          E5CHEMSPECARR_tte(i)%longname = chch_defs(i)%speci
!!$          !E5CHEMSPECARR_tte(i)%units = 'm^3/m^3'
!!$       ENDDO
!!$       !!$ CLOSE(iou)
!!$    ENDIF
!!$    IF (USE_CLAMSCHEME5) THEN
!!#D clamschem +
    IF (L_CLAMSCHEME5) THEN
       ALLOCATE (E5CHEMSPECARR_tte(nchemspec))
       DO i = 1, nchemspec
          E5CHEMSPECARR_tte(i)%name = chch_defs(i)%speci
          E5CHEMSPECARR_tte(i)%ctype = chch_defs(i)%ctype
          E5CHEMSPECARR_tte(i)%longname = chch_defs(i)%speci
          !E5CHEMSPECARR_tte(i)%units = 'm^3/m^3'
       ENDDO
       CLOSE(iou)
    ENDIF
!!#D clamschem -
! op_pj_20170110-

  END SUBROUTINE clamsbmix_initialize
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsbmix_init_memory

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    USE messy_main_channel,          ONLY: new_channel, new_attribute &
                                         , new_channel_object
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_clams_tools_utils,     ONLY: uppercase

    IMPLICIT NONE

    ! LOCAL 
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsbmix_init_memory'
    INTEGER :: status

    IF (p_pe==0) write(*,*) uppercase(substr)

    ! Define channel BMIX
    CALL new_channel  (status, modstr)
    CALL channel_halt (substr, status)
    CALL new_attribute(status, modstr, modstr//'_version', c = modver)
    CALL channel_halt (substr, status)

  END SUBROUTINE clamsbmix_init_memory
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsbmix_init_coupling

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: get_channel_object, get_channel_info
    USE messy_main_switch,           ONLY: USE_CLAMSCIRRUS, &
                                           USE_CLAMSMIX

    USE messy_clams_global,          ONLY: nspec, SPECARR, nchemspec, &
                                           nparams, nparams_old
    USE messy_clamsmix_global,       ONLY: switch_mixing
    USE messy_clams_tools_utils,     ONLY: uppercase
    USE messy_main_mpi_bi,           ONLY: p_pe
!!#D clamschem +
    USE messy_clamschem_defs_mod,    ONLY: chch_defs ! op_pj_20170110
!!#D clamschem -

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsbmix_init_coupling'
    integer :: status, i

    IF (p_pe==0) write(*,*) uppercase(substr)

    ! op_sb_20191021+
    ALLOCATE (MIXSPECARR(nspec))
    ALLOCATE (MIXSPECARR_SHUFFLED(nspec))
    ! op_sb_20191021-


    ! Get arrays from CLAMS submodel:
    CALL get_channel_object(status, 'clams', 'LAT', p1=LAT)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON', p1=LON)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LEV', p1=LEV)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'JULTIME', p0=JULTIME)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LAT_OLD', p1=LAT_OLD)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON_OLD', p1=LON_OLD)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LEV_OLD', p1=LEV_OLD)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'JULTIME_OLD', p0=JULTIME_OLD)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LAT_OLD_MIX', p1=LAT_OLD_MIX)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON_OLD_MIX', p1=LON_OLD_MIX)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'STATE', p1=STATE)
    CALL channel_halt(substr, status)

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

!!$    CALL get_channel_object(status, 'clams', 'TEMP_OLD', p1=TEMP_OLD)
!!$    CALL channel_halt(substr, status)
!!$    IF (TRIM(init_vertcoorname) == 'press') THEN 
!!$       CALL get_channel_object(status, 'clams', 'LEV_OLD', p1=PRESS_OLD)
!!$       CALL channel_halt(substr, status)
!!$    ELSE
!!$       CALL get_channel_object(status, 'clams', 'PRESS_OLD', p1=PRESS_OLD)
!!$       CALL channel_halt(substr, status)
!!$    END IF

    ! THETA_OLD_MIX, BVF_OLD_MIX (if vertical mixing is switched on)
    if (USE_CLAMSMIX .and. switch_mixing==2) then
       CALL get_channel_object(status, 'clams', 'THETA_OLD_MIX', p1=THETA_OLD_MIX)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'clams', 'BVF_OLD_MIX', p1=BVF_OLD_MIX)
       CALL channel_halt(substr, status)
    endif
   
    ! Get chemical species (from CLAMS):
    DO i = 1, nspec
       CALL get_channel_object(status, 'clams', &
            trim(SPECARR(i)%name), p1=MIXSPECARR(i)%values)
       CALL channel_halt(substr, status)
       MIXSPECARR(i)%name     = SPECARR(i)%name
       MIXSPECARR(i)%longname = SPECARR(i)%longname
       MIXSPECARR(i)%units    = SPECARR(i)%units
       MIXSPECARR(i)%ctype    = SPECARR(i)%ctype
    ENDDO

    ! Shuffled arrays
    CALL get_channel_object(status, 'clamsmix', 'ZETA_SHUFFLED', p1=LEV_SHUFFLED)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clamsmix', 'LAT_SHUFFLED', p1=LAT_SHUFFLED)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clamsmix', 'LON_SHUFFLED', p1=LON_SHUFFLED)
    CALL channel_halt(substr, status)
!!!!! ???
!!$    CALL get_channel_object(status, 'clamsmix', 'LAT_OLD_MIX_SHUFFLED', p1=LAT_OLD_MIX_SHUFFLED)
!!$    CALL channel_halt(substr, status)
!!$    CALL get_channel_object(status, 'clamsmix', 'LON_OLD_MIX_SHUFFLED', p1=LON_OLD_MIX_SHUFFLED)
!!$    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clamsmix', 'STATE_SHUFFLED', p1=STATE_SHUFFLED)
    CALL channel_halt(substr, status)
    DO i = 1, nspec
       CALL get_channel_object(status, 'clamsmix', trim(SPECARR(i)%name)//'_SHUFFLED',&  
            p1=MIXSPECARR_SHUFFLED(i)%values)
       CALL channel_halt(substr, status)
       MIXSPECARR_SHUFFLED(i)%name     = SPECARR(i)%name
       MIXSPECARR_SHUFFLED(i)%longname = SPECARR(i)%longname
       MIXSPECARR_SHUFFLED(i)%units    = SPECARR(i)%units
       MIXSPECARR_SHUFFLED(i)%ctype    = SPECARR(i)%ctype
    ENDDO

    ! op_pj_20170110+
    CALL get_channel_info(status,'clamscheme5')
    L_CLAMSCHEME5 = (status == 0)
    ! op_pj_20170110-

#ifdef ECHAM5
!!$ IF (USE_CLAMSCHEME5) THEN ! op_pj_20170110
    IF (L_CLAMSCHEME5) THEN   ! op_pj_20170110
       CALL get_channel_object(status, 'clams', 'E5_ZETA', p3=E5_ZETA_D)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'clams', 'E5_LAT3D', p3=E5_LAT3D_D)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'clams', 'E5_LON3D', p3=E5_LON3D_D)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'ECHAM5', 'press', p3=E5_PRESS_D)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'ECHAM5', 'qm1', p3=E5_SH_D)
       CALL channel_halt(substr, status)
    END IF

    ! op_pj_20170110+
!!#D clamschem +
    IF (L_CLAMSCHEME5) THEN
       ALLOCATE (E5CHEMSPECARR_tte(nchemspec))
       DO i = 1, nchemspec
          E5CHEMSPECARR_tte(i)%name = chch_defs(i)%speci
          E5CHEMSPECARR_tte(i)%ctype = chch_defs(i)%ctype
          E5CHEMSPECARR_tte(i)%longname = chch_defs(i)%speci
          !E5CHEMSPECARR_tte(i)%units = 'm^3/m^3'
       ENDDO
    ENDIF
!!#D clamschem -
    ! op_pj_20170110-

    IF (L_CLAMSCHEME5) THEN
       DO i=1, nchemspec
          IF (E5CHEMSPECARR_tte(i)%ctype == 'TR') THEN
             CALL get_channel_object(status, 'clamscheme5', TRIM(E5CHEMSPECARR_tte(i)%name)//'_tte',&
               p3=E5CHEMSPECARR_tte(i)%values)
             CALL channel_halt(substr, status)
          END IF
       END DO
    ENDIF
#endif
 
  END SUBROUTINE clamsbmix_init_coupling
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsbmix_global_end

    ! BMIL
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_mpi_bi,           ONLY: p_sum, p_pe
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_tracer_mem_bi,    ONLY: xtm1, xtte
    USE messy_main_grid_def_mem_bi,  ONLY: nlev, ngpblks, nproma, npromz

    ! SMIL
    USE messy_clamsmix_si,           ONLY: clamsmix_reshuffle

    ! SMCL
    USE messy_main_switch,           ONLY: USE_CLAMSCIRRUS, &
!!$ USE messy_main_switch,         ONLY: USE_CLAMSCIRRUS, USE_CLAMSCHEME5, &
                                         USE_CLAMSMIX
    USE messy_main_timer,            ONLY: YEAR, MONTH, DAY, HOUR, &
                                           MINUTE, SECOND, delta_time, time_step_len
    USE messy_main_channel,          ONLY: set_channel_output

    USE messy_clams_global,          ONLY: nparts, dnparts, dnparts_max, mdi, ldiagout, &
                                         timetype, dates30, irdsec, irdday, lcoupled, &
                                         fut_year, fut_month, fut_day, fut_sec, &
                                         lbmixevent, lmixevent, nchemspec, &
                                         nparams, paramnames, nparams_old, paramnames_old, &
                                         eps, dnparts_co, nspec, H2O_index,&
                                         H2O_100_index, IWC_index, IWC_100_index, CLWC_index
    USE messy_clamsbmix_global,      ONLY: interpol_from_init, replace_low, replace_up, &
                                         lev_down, lev_in_down, lev_up, lev_in_up, &
                                         nclamsbounds,  &
                                         switch_EMAC_H2O, EMAC_H2O_z
    USE messy_clamsmix_global,       ONLY: levelrange,l_min_act,l_max_act,l_delta_act, &
                                         timestep_mix, adapt_par, switch_mixing, &
                                         vert_mix_param
                                         
    USE messy_clams_tools_interpolreg,  ONLY: interpolate_param
    USE messy_clams_tools_dateconv,     ONLY: incdat_interface
    USE messy_clams_tools_utils,        ONLY: str_pos, uppercase

    USE messy_clamscirrus,              ONLY: CIRRUS
    USE messy_clamscirrus_global,       ONLY: use_traj

    USE messy_clamsbmix_global,         ONLY: nlevs
    USE messy_clamsbmix_replace_bounds, ONLY: replace_boundaries

    USE netcdf

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsbmix_global_end'

    TYPE(timetype) :: itime, ipasttime
    INTEGER        :: status, i, j, k, ilon, ilev, ilat, icnt, ncnt
#ifdef ECHAM5
    REAL(DP)       :: h2o_temp(dnparts_max)
    INTEGER        :: ipart, counter
#endif
    INTEGER        :: rcode, ncid, varid
    CHARACTER(150) :: ncdf_file                      
    INTEGER        :: ispec
!    LOGICAL        :: asc_lat_clim, loglev_clim, logpress_clim, use_modellev_clim
    INTEGER        :: l1, l2, l3, l4
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: temp_array
    LOGICAL,  DIMENSION(:)      , ALLOCATABLE :: calc
    LOGICAL,  DIMENSION(:)      , ALLOCATABLE :: calc0
    INTEGER,DIMENSION(NF90_MAX_VAR_DIMS)      :: dim_array

    IF (lbmixevent) THEN

       IF (p_pe==0) WRITE(*,*) 'ACTIVE BMIX EVENT'

       nlevs = adapt_par%nlevs

       !--------------------------------------------------------------------
       ! define lev- and lat-boundaries
       !--------------------------------------------------------------------
       call define_boundaries


       !--------------------------------------------------------------------
       ! Save lowest and top level for interpolation in BMIX
       ! if interpolation from init-file is chosen
       !--------------------------------------------------------------------
!!!!! Das unterste/oberste Level muss als Ganzes gespeichert werden, wenn
!!!!! Interpolation des untersten/obersten Levels aus dem Init-File gewuenscht ist !!!
       if ((interpol_from_init==1 .or. interpol_from_init==3) .and. replace_low) then
          ! save lowest layer for interpolation
          ! (all points between lev_down and lev_in_down on all ranks)
          call save_layer (lev_down, lev_in_down, &
                           LAT_SHUFFLED, LON_SHUFFLED, LEV_SHUFFLED, MIXSPECARR_SHUFFLED,  &
                           LAT_LOW, LON_LOW, LEV_LOW, SPECARR_LOW)
       else
          ! allocate dummy arrays 
          allocate(LAT_LOW(1),LON_LOW(1),LEV_LOW(1),SPECARR_LOW(1))
       endif

       if ((interpol_from_init==2 .or. interpol_from_init==3) .and. replace_up) then
          ! save upper layer for interpolation
          ! (all points between lev_in_up and lev_up on all ranks)
          call save_layer (lev_in_up, lev_up, &
                           LAT_SHUFFLED, LON_SHUFFLED, LEV_SHUFFLED, MIXSPECARR_SHUFFLED,  &
                           LAT_UPPER, LON_UPPER, LEV_UPPER, SPECARR_UPPER)
       else
          ! allocate dummy arrays 
          allocate(LAT_UPPER(1),LON_UPPER(1),LEV_UPPER(1),SPECARR_UPPER(1))
       endif

       !write (*,*) p_pe,'size(lat_low)=',size(lat_low)
       
       IF (.not. USE_CLAMSMIX) THEN
          ! Set STATE to 0
          STATE      = 0
          STATE_SHUFFLED      = 0
       ENDIF

       !--------------------------------------------------------------------
       ! Call BMIX
       !--------------------------------------------------------------------
       call bmix (status, LAT_SHUFFLED, LON_SHUFFLED, LEV_SHUFFLED, &
                  STATE_SHUFFLED, MIXSPECARR_SHUFFLED, &
                  LAT_LOW, LON_LOW, LEV_LOW, SPECARR_LOW, &
                  LAT_UPPER, LON_UPPER, LEV_UPPER, SPECARR_UPPER, lcoupled)
       IF (status /= 0) CALL error_bi('Error in BMIX !!!',substr)


!!!!!
       dnparts_co(1) = dnparts
       nparts = p_sum(dnparts)
       if (p_pe==0) write (*,*) 'nach bmix: rank, dnparts, nparts =',p_pe,dnparts,nparts
!!!!!
       !--------------------------------------------------------------------
       ! Reshuffle
       !--------------------------------------------------------------------
       if (p_pe==0) write (*,*) 'call clamsmix_reshuffle'
       call clamsmix_reshuffle

       nparts = p_sum(dnparts)
       dnparts_co(1) = dnparts
       if (p_pe==0) write (*,*) 'nach reshuffle: rank, dnparts, nparts=',p_pe,dnparts,nparts

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
          CALL interpolate_param(LON, LAT, LEV, itime, ipasttime, &
               i, paramnames(i), PARAM(i)%values, lcoupled)
       END DO


       !--------------------------------------------------------------------
       ! Set boundaries
       !--------------------------------------------------------------------
       if (nclamsbounds > 0) then
          call replace_boundaries (status, LAT, LON, LEV, PARAM, MIXSPECARR)
          IF (status /= 0) CALL error_bi('Error in replace_boundaries !!!',substr)
       endif


       !--------------------------------------------------------------------
       ! call CIRRUS
       !--------------------------------------------------------------------
       IF (USE_CLAMSCIRRUS .AND. lmixevent) THEN
          IF (dnparts > 0) THEN
             use_traj = .FALSE. 
              CALL CIRRUS(LAT, LON, &
                  PARAM(str_pos(nparams,paramnames,'THETA'))%values, &
                  PARAM(str_pos(nparams,paramnames,'TEMP'))%values,&
                  PARAM(str_pos(nparams,paramnames,'PRESS'))%values, &
                  PARAM(str_pos(nparams,paramnames,'PV'))%values, &
                  MIXSPECARR(H2O_index)%values, &
                  MIXSPECARR(H2O_100_index)%values, MIXSPECARR(IWC_index)%values, &
                  MIXSPECARR(IWC_100_index)%values, MIXSPECARR(CLWC_index)%values, &
                  timestep_mix*3600._PREC, dnparts)
         END IF
       END IF

       !--------------------------------------------------------------------
       ! Tropospheric H2O from EMAC 
       !--------------------------------------------------------------------
#ifdef ECHAM5
       IF (switch_EMAC_H2O) THEN
          ! CLaMS airparcels:
          ALLOCATE(calc(dnparts_max))
          WHERE (LEV .LE. EMAC_H2O_z)
             calc = .TRUE.
          END WHERE

          h2o_temp = mdi
         CALL interpolate_param (LON, LAT, LEV, itime, ipasttime, &
               str_pos(nparams,paramnames,'H2O'), 'H2O', h2o_temp, lcoupled)
          
          DO ipart = 1, dnparts
             IF (calc(ipart)) THEN
                MIXSPECARR(H2O_100_index)%values(ipart)  = h2o_temp(ipart)
                MIXSPECARR(H2O_index)%values(ipart)      = h2o_temp(ipart)
             END IF
          END DO
          DEALLOCATE(calc)
       END IF
#endif

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
       LEV_OLD = LEV
       JULTIME_OLD = JULTIME

       ! TEMP_OLD and PRESS_OLD
       DO i = 1, nparams_old
          PARAM_OLD(i)%values = PARAM(str_pos(nparams,paramnames,paramnames_old(i)))%values
          !write (*,*) 'iparam_old, iparam:', i, str_pos(nparams,paramnames,paramnames_old(i))
       END DO



!!!!!       
!!$       !--------------------------------------------------------------------
!!$       ! write output
!!$       !--------------------------------------------------------------------
!!$       IF (loutput_bmix) THEN
!!$          CALL set_channel_output(status, 'clamsbmix', .TRUE.)
!!$          CALL channel_halt(substr, status)
!!$       ENDIF
!!$       IF (USE_CLAMSCIRRUS .and. loutput_cirrus) THEN
!!$          CALL set_channel_output(status, 'clamscirrus', .TRUE.)
!!$          CALL channel_halt(substr, status)
!!$       END IF
!!$       IF (L_CLAMSCHEME5) THEN
!!$          CALL set_channel_output(status, 'tracer_gp', .TRUE.)
!!$          CALL channel_halt(substr, status)
!!$       END IF
!!$       IF (loutput_bmix .or. (USE_CLAMSCIRRUS .and. loutput_cirrus) .or. L_CLAMSCHEME5) THEN
!!$          CALL set_channel_output(status, 'clams', .FALSE.)
!!$          CALL channel_halt(substr, status)
!!$          CALL messy_write_output
!!$          if (lclamsoutevent) then
!!$             CALL set_channel_output(status, 'clams', .TRUE.)
!!$             CALL channel_halt(substr, status)
!!$          endif
!!$       END IF
    
       ! Clean up
       DEALLOCATE(levelrange)
       DEALLOCATE(l_min_act, l_max_act, l_delta_act)

       if ((interpol_from_init==1 .or. interpol_from_init==3) .and. replace_low) then
          DO i = 1, nspec
             DEALLOCATE (SPECARR_LOW(i)%values)
          ENDDO
       endif
       if ((interpol_from_init==2 .or. interpol_from_init==3) .and. replace_up) then
          DO i = 1, nspec
             DEALLOCATE (SPECARR_UPPER(i)%values)
          ENDDO
       endif
       DEALLOCATE (LAT_LOW, LON_LOW, LEV_LOW, SPECARR_LOW)
       DEALLOCATE (LAT_UPPER, LON_UPPER, LEV_UPPER, SPECARR_UPPER)

    END IF
  
#ifdef ECHAM5
!!$ IF (USE_CLAMSCHEME5) THEN
    IF (L_CLAMSCHEME5) THEN
!!$       DO ispec = 1, nchemspec
!!$          ! Upper boundary: 
!!$          IF (int_up_clim(ispec) .EQ. 1 .OR. int_up_clim(ispec) .EQ. 2) THEN
!!$             SELECT CASE (coor_up_clim(ispec))
!!$             CASE ('Z')
!!$                WHERE (E5_ZETA_D .GE. bound_up_clim(ispec))
!!$                   xtte(:,:,ispec,:) = E5CHEMSPECARR_tte(ispec)%values(:,:,:)
!!$                END WHERE
!!$             CASE('P')
!!$                WHERE (E5_PRESS_D/100. .LE. bound_up_clim(ispec) .AND.&
!!$                     ABS((E5_PRESS_D/100.-mdi)/mdi)>eps)
!!$                   xtte(:,:,ispec,:) = E5CHEMSPECARR_tte(ispec)%values(:,:,:)
!!$                END WHERE
!!$             CASE DEFAULT
!!$                WRITE(*,*) 'Error: Wrong coor_up_clim!'
!!$                STOP
!!$             END SELECT
!!$          END IF
!!$          ! Lower boundary: 
!!$          IF (int_low_clim(ispec) .EQ. 1 .OR. int_low_clim(ispec) .EQ. 2) THEN
!!$             SELECT CASE (coor_low_clim(ispec))
!!$             CASE ('Z')
!!$                WHERE (E5_ZETA_D .LE. bound_low_clim(ispec).AND.&
!!$                     ABS((E5_ZETA_D-mdi)/mdi)>eps)
!!$                   xtte(:,:,ispec,:) = E5CHEMSPECARR_tte(ispec)%values(:,:,:)
!!$                END WHERE
!!$             CASE('P')
!!$                WHERE (E5_PRESS_D/100. .GE. bound_low_clim(ispec))
!!$                   xtte(:,:,ispec,:) = E5CHEMSPECARR_tte(ispec)%values(:,:,:)
!!$                END WHERE
!!$             CASE DEFAULT
!!$                WRITE(*,*) 'Error: Wrong coor_low_clim!'
!!$                STOP
!!$             END SELECT
!!$          END IF
!!$       END DO

       ! Tropospheric water vapor from EMAC:
       WHERE (E5_ZETA_D .LE. EMAC_H2O_z)
          xtte(:,:,H2O_index,:) =  &
               (E5_SH_D(:,:,:)*28.9644/18.015 - xtm1(:,:,H2O_index,:))/time_step_len
          xtte(:,:,H2O_100_index,:) =  &
               (E5_SH_D(:,:,:)*28.9644/18.015 - xtm1(:,:,H2O_100_index,:))/time_step_len
       END WHERE

    END IF

#endif

  END SUBROUTINE clamsbmix_global_end
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsbmix_free_memory

!!$ USE messy_main_switch,       ONLY: USE_CLAMSCHEME5
    USE messy_clams_global,      ONLY: nparams_old

    IMPLICIT NONE

    DEALLOCATE (PARAM)
    IF (nparams_old > 0)  DEALLOCATE (PARAM_OLD)
    DEALLOCATE (MIXSPECARR)
    DEALLOCATE (MIXSPECARR_SHUFFLED)

!!$ IF (USE_CLAMSCHEME5) THEN
    IF (L_CLAMSCHEME5) THEN
       DEALLOCATE (E5CHEMSPECARR_tte)
    ENDIF

  END SUBROUTINE clamsbmix_free_memory
!--------------------------------------------------------------------

  !****************************************************************************
  ! Save layer for interpolation of species
  !****************************************************************************
  SUBROUTINE save_layer (lev_min, lev_max, lat, lon, lev, specarr,  &
                         lat_layer, lon_layer, lev_layer, specarr_layer)
    ! SMCL
    use messy_clams_global,    only: prec, rank, nspec, species_type
    use messy_clamsbmix_tools, only: get_nparts_layer

    implicit none

    REAL(PREC)                                :: lev_min, lev_max
    REAL(PREC),         DIMENSION(:), POINTER :: lat,lon,lev
    TYPE(species_type), DIMENSION(:), POINTER :: specarr
    REAL(PREC),         DIMENSION(:), POINTER :: lat_layer, lon_layer, lev_layer
    TYPE(SPECIES_TYPE), DIMENSION(:), POINTER :: specarr_layer
 
    integer,     dimension(:),   pointer  :: indices

    integer :: nparts_layer
    integer :: i
    
    ! get indices of points between lev_min and lev_max (on current rank)
    call get_nparts_layer (lev, indices, lev_min, lev_max, nparts_layer)

    ! gather all points on layer from all ranks
    call gather_layer (lat_layer, lat, indices, nparts_layer)
    call gather_layer (lon_layer, lon, indices, nparts_layer)
    call gather_layer (lev_layer, lev, indices, nparts_layer)
    allocate (specarr_layer(nspec))
    do i = 1, nspec
       call gather_layer (specarr_layer(i)%values, specarr(i)%values, indices, nparts_layer)
    enddo

    ! clean up
    if (nparts_layer>0) deallocate (indices)

  END SUBROUTINE save_layer

  !****************************************************************************
  ! Gather all points on lowest layer from all ranks
  !****************************************************************************
  subroutine gather_layer (all_array, array, indices, my_nparts)

    ! BMIL
    USE messy_main_mpi_bi,         ONLY: p_sum, p_bcast, p_send, p_recv
    USE messy_main_blather_bi,     ONLY: error_bi

    ! SMCL
    USE messy_clams_global,        ONLY: rank, ntasks

    implicit none

    REAL(PREC),  DIMENSION(:), POINTER :: all_array, array
    INTEGER,     DIMENSION(:), POINTER :: indices
    INTEGER                            :: my_nparts

    REAL(PREC),  DIMENSION(:), POINTER :: my_array
    INTEGER :: all_nparts, nparts_layer
    INTEGER :: irank, startpos

    integer, save :: tag = 100.

    ! my_nparts = number of points on layer (on current rank)
 
    ! number of points on layer (on all ranks)
    all_nparts = p_sum (my_nparts)

    !if (rank==0) write (*,*) 'in gather_layer: all_nparts=', rank, all_nparts

    ! No points on layer:
    if (all_nparts==0) &
         CALL error_bi('Not enough APs in layer','gather_layer')

    ! allocate array for all points on layer (on all ranks)
    allocate (all_array(all_nparts))
    
    ! allocate array for points on layer on current rank
    if (my_nparts>0) then
       allocate (my_array(my_nparts))
       my_array = array(indices)
    endif

    ! set tag for p_send/p_recv
    tag = tag + 2*ntasks

    ! writing position in array "all_array"
    startpos = 1

    if (rank==0) then

       ! write own data on layer to all_array
       if (my_nparts>0) then
          all_array (1:my_nparts) = my_array
          startpos = startpos+my_nparts
       endif

       ! gather data on layer from all ranks 
       do irank = 1, ntasks-1

          ! get number of points on layer (on rank irank)
          call p_recv (nparts_layer,irank,tag+irank)

          ! get points on layer from rank irank
          if (nparts_layer>0) then             
             call p_recv (all_array(startpos:startpos+nparts_layer-1),irank,tag+irank+ntasks)
             startpos = startpos + nparts_layer
          endif

       enddo

    else
       
       ! send number of points on layer to rank 0
       call p_send (my_nparts,0,tag+rank)

       ! send points on layer to rank 0
       if (my_nparts>0) then
          call p_send (my_array,0,tag+rank+ntasks)
       endif

    endif

    ! Broadcast all_array to all ranks
    call p_bcast (all_array,0)

    if (my_nparts>0) deallocate (my_array)

  end subroutine gather_layer

#endif
!**********************************************************************
END MODULE messy_clamsbmix_si
!**********************************************************************
