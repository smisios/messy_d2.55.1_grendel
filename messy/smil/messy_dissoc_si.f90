#if defined (ECHAM5) || defined (COSMO) || defined (BLANK) || defined(MESSYDWARF)
#include "messy_main_ppd_bi.inc"
#endif

!**********************************************************************
MODULE messy_dissoc_si
!**********************************************************************
!  Submodel interface for dissoc -- calculation of photolysis rates
! Authors:
! Thomas Breuer,     IEK-7, Jul 2012
! Jens-Uwe Grooss, IEK-7, Forschungszentrum Juelich, Januar 2014
! Patrick Joeckel, DLR, July 2015, first draft of coupling to EMAC ...
!**********************************************************************

!**********************************************************************
! This SMIL module serves as universal plug-in for various cases:
! (1) coupling to grid-point models (ECHAM5, COSMO, ...):
!     statements specific for this case need to be encapsulated
!     into preprocessor directives #if defined(ECHAM5) || defined(...)
! (2) coupling to the MBM CLaMS and / or to EMAC/CLaMS:
!     statements specific for this case need to be encapsulted
!     into (IF (l_clams) THEN ... ELSE ... ENDIF)
!     The logical l_clams can be set by checking for the CLAMS channel ...
! (3) coupling to the MBM DISSOC (stand allone application)
!     statements specific for this case need to be encapsulated into
!     preprocessor directives #ifdef MBM_DISSOC
!
!**********************************************************************
! TODO: a) the calculation of the solar zenith angle (time) needs to be
!          checked ...
!       b) the import of O3 climatology should be replaced by a proper
!          on-line coupling to an external ozone field which can be selected
!          by a CPL-namelist switch (TYPE(t_chaobj_cpl) :: dissoc_O3), see
!          below; this can then be either the on-line ozone tracer or
!          a climatology provided via IMPORT_GRID (i.e., on the native
!          model grid)
!       c) coupling to MBM CLaMS / EMAC/CLaMS needs to be finalised and tested
!       d) restart test with EMAC must be performed
!**********************************************************************

  USE messy_main_blather_bi,           ONLY: start_message_bi, end_message_bi &
                                           , error_bi
  USE messy_main_channel_dimensions,   ONLY: DIMID_UNDEF
  USE messy_main_channel_repr,         ONLY: REPR_UNDEF
  USE messy_main_channel,              ONLY: t_chaobj_cpl
  USE messy_main_tools,                ONLY: PTR_3D_ARRAY
  USE messy_dissoc
  USE messy_main_timer_event,          ONLY: time_event, io_time_event

  IMPLICIT NONE
  PRIVATE
  SAVE
  INTRINSIC :: NULL

  ! DISSOC DIMENSIONS
  INTEGER :: DIMID_WAVE        = DIMID_UNDEF
  INTEGER :: DIMID_DISSOC_LEV  = DIMID_UNDEF
  INTEGER :: DIMID_DISSOC_LEVC = DIMID_UNDEF
  INTEGER :: DIMID_DISSOC_LAT  = DIMID_UNDEF
  INTEGER :: DIMID_DISSOC_SZA  = DIMID_UNDEF

  ! DISSOC REPRESENTATIONS
  INTEGER :: REPR_DISSOC_WAVE  = REPR_UNDEF
  INTEGER :: REPR_DISSOC_LEV   = REPR_UNDEF
  INTEGER :: REPR_DISSOC_LEVC  = REPR_UNDEF
  INTEGER :: REPR_DISSOC_LAT   = REPR_UNDEF
  INTEGER :: REPR_DISSOC_SZA   = REPR_UNDEF
  INTEGER :: REPR_DISSOC_TABS  = REPR_UNDEF
  INTEGER :: REPR_DISSOC_TABS_DAVG = REPR_UNDEF
  INTEGER :: REPR_DISSOC_TWOD  = REPR_UNDEF
  INTEGER :: REPR_DISSOC_TWODC = REPR_UNDEF

  ! POINTERS TO CHANNEL OBJECTS FOR LOOKUP-TABLES ETC.
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: dissoc_tabs => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER :: dissoc_tabs_davg => NULL()

  REAL(DP), DIMENSION(:), POINTER :: dissoc_szagrid => NULL()
  REAL(DP), DIMENSION(:), POINTER :: dissoc_latgrid => NULL()
  REAL(DP), DIMENSION(:), POINTER :: dissoc_lambda => NULL()
  REAL(DP), DIMENSION(:), POINTER :: dissoc_plevs => NULL()
  REAL(DP), DIMENSION(:), POINTER :: dissoc_plevsc => NULL()

  REAL(DP), DIMENSION(:,:), POINTER :: dissoc_alt => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: dissoc_altc => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: dissoc_temp => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: dissoc_tempc => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: dissoc_do2 => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: dissoc_do2c => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: dissoc_do3 => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: dissoc_do3c => NULL()

  ! TIMER: here this timer is used to trigger (hard-wired!) updates
  !        of read-in O3 climatology at the beginning of each month
  TYPE(io_time_event) :: TIMER_MONTHLY = io_time_event(1, 'months','first',0)
  TYPE(time_event)    :: XTIMER_MONTHLY
  LOGICAL             :: LTRIG_MONTHLY = .FALSE.

  ! TIME INFORAMTION
  REAL(DP) :: jdfrac ! julian day plus fraction
  REAL(DP) :: jsec   ! ... in seconds

! op_j_20150727: currently not yet used;
!                requires flexible O3 input for DISSOC SMCL
!!$  ! CPL namelist
!!$  ! - name of ozone tracer/channel object
!!$  TYPE(t_chaobj_cpl) :: dissoc_O3

  ! op_pj_20150710: new channel objects
  ! POINTERS FOR GRID-POINT REPRESENTATION
  TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER  :: jdiss_gp => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_sza => NULL() ! solar zenith angle

  ! POINTERS FOR LAGRANGIAN REPRESENTATION
  LOGICAL :: l_clams = .FALSE.  ! is CLaMS running ?
  REAL(DP), DIMENSION(:), POINTER :: LAT        => NULL()
  REAL(DP), DIMENSION(:), POINTER :: LON        => NULL()
  REAL(DP), DIMENSION(:), POINTER :: TEMP       => NULL()
  REAL(DP), DIMENSION(:), POINTER :: PRESS      => NULL()
  REAL(DP), POINTER :: JULTIME     => NULL()
  TYPE(rate_type), DIMENSION(IP_MAX) :: jresult

  ! ju_nt_20160217: new channel object SZA for CLaMS
  REAL(DP), DIMENSION(:), POINTER :: SZA        => NULL() ! solar zenith angle

  ! ju_nt_2016_0217: dnparts, dnparts_max coupled from channel "clams"
  REAL(DP), DIMENSION(:), POINTER :: dnparts_co
  REAL(DP), DIMENSION(:), POINTER :: dnparts_max_co
  INTEGER :: dnparts, dnparts_max


#ifdef MBM_DISSOC
  REAL(DP), DIMENSION(:,:,:), POINTER :: tm1_3d => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: tte_3d => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: press_3d => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: altitude_gnd => NULL()
  REAL(DP), DIMENSION(:,:),   POINTER :: philat_2d => NULL()
  REAL(DP), DIMENSION(:,:),   POINTER :: philon_2d => NULL()
!!$REAL(DP), DIMENSION(:,:),   POINTER :: albedo_2d => NULL() ! op_pj_20160617
#endif
  REAL(DP), DIMENSION(:,:),   POINTER :: albedo_2d => NULL()  ! op_j_20160617

  ! ju_nt_20180620+: add dissoc event
  ! ub_ak_20190709+
  ! remove local SAVE, global SAVE exists (g95 error)
  !TYPE(time_event), PUBLIC, SAVE :: dissocevent, dissocevent_gp
  !TYPE(io_time_event), SAVE :: io_dissocevent, io_dissocevent_gp
  TYPE(time_event), PUBLIC :: dissocevent, dissocevent_gp
  TYPE(io_time_event)      :: io_dissocevent, io_dissocevent_gp
  ! ub_ak_20190709-
  LOGICAL :: ldissocevent, ldissocevent_gp
  ! ju_nt_20180620-

  PUBLIC :: dissoc_initialize
  PUBLIC :: dissoc_init_memory
  PUBLIC :: dissoc_init_coupling
  PUBLIC :: dissoc_global_start
  PUBLIC :: dissoc_physc
  PUBLIC :: dissoc_global_end
  PUBLIC :: dissoc_free_memory
  !PRIVATE :: dissoc_initialize_dims
  !PRIVATE :: dissoc_initialize_reprs

!--------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE dissoc_initialize

    USE messy_main_mpi_bi,     ONLY: p_io, p_parallel_io, p_bcast, p_pe
    USE messy_main_tools,      ONLY: find_next_free_unit
    USE messy_main_timer_bi,   ONLY: timer_event_init
    USE messy_main_timer,      ONLY: delta_time

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'dissoc_initialize'
    INTEGER :: status, iou

    ! Read namelist and set default values:
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL dissoc_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('Error in dissoc_read_nml_ctrl ',substr)
    END IF
    CALL p_bcast(o3dat, p_io)
    CALL p_bcast(davg, p_io)
    CALL p_bcast(datdir, p_io)

    ! ju_nt_2016_0217: read mean albedo from namelist
    CALL p_bcast(mean_albedo, p_io)
    albedo = mean_albedo

    ! ju_nt_20180620+
    ! read timesteps from namelist and define events
    CALL p_bcast(timestep_dissoc, p_io)
    CALL p_bcast(timestep_dissoc_gp, p_io)
    if (timestep_dissoc < delta_time) then
       if (p_pe==0) then
          write (*,*)
          write (*,*) 'timestep_dissoc < delta_time'
          write (*,*) 'timestep_dissoc:', timestep_dissoc
          write (*,*) 'delta_time:', delta_time
          write (*,*) 'timestep_dissoc is set to:', delta_time
          write (*,*)
       endif
       timestep_dissoc = delta_time
    elseif (mod(timestep_dissoc,int(delta_time)) /= 0) then
       call error_bi ("Wrong dissoc timestep !!!",substr)
    endif
    if (timestep_dissoc > delta_time) then
       io_dissocevent%counter = timestep_dissoc
       io_dissocevent%unit = 'seconds'
       io_dissocevent%adjustment = 'exact'
       io_dissocevent%offset = -delta_time
       CALL timer_event_init (dissocevent, io_dissocevent, 'Dissoc_Event', 'present')
    endif
    if (timestep_dissoc_gp < delta_time) then
       if (p_pe==0) then
          write (*,*)
          write (*,*) 'timestep_dissoc_gp < delta_time'
          write (*,*) 'timestep_dissoc_gp:', timestep_dissoc_gp
          write (*,*) 'delta_time:', delta_time
          write (*,*) 'timestep_dissoc_gp is set to:', delta_time
          write (*,*)
       endif
       timestep_dissoc_gp = delta_time
    elseif (mod(timestep_dissoc_gp,int(delta_time)) /= 0) then
       call error_bi ("Wrong dissoc timestep !!!",substr)
    endif
    if (timestep_dissoc_gp > delta_time) then
       io_dissocevent_gp%counter = timestep_dissoc_gp
       io_dissocevent_gp%unit = 'seconds'
       io_dissocevent_gp%adjustment = 'exact'
       io_dissocevent_gp%offset = -delta_time
       CALL timer_event_init (dissocevent_gp, io_dissocevent_gp, 'Dissoc_Event_gp', 'present')
    endif
    ! ju_nt_20180620-

! op_pj_20150727: currently not used, see below
!!$    IF (p_parallel_io) THEN
!!$       iou = find_next_free_unit(100,200)
!!$       CALL dissoc_read_nml_cpl(status, iou)
!!$       IF (status /= 0) CALL error_bi(' ',substr)
!!$    END IF
!!$    ! broadcast results
!!$    CALL p_bcast(dissoc_O3%cha, p_io)
!!$    CALL p_bcast(dissoc_O3%obj, p_io)

    ! initialise (hard-wired) trigger for monthly updates of external
    ! ozone climatology ... this should later become obsolete, if
    ! proper ozone coupling has been implemented ...
    CALL timer_event_init(XTIMER_MONTHLY, TIMER_MONTHLY &
         , 'dissoc_monthly', 'present')

  END SUBROUTINE dissoc_initialize
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE dissoc_init_memory

    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         ,  new_attribute, STRLEN_OBJECT   &
                                         , get_channel_info
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, REPR_LG_CLAMS
#ifdef MBM_DISSOC
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL
    USE messy_main_grid_def_mem_bi,  ONLY: nlon, nlat, nlev
    USE messy_main_constants_mem,    ONLY: g
#endif
    USE messy_main_blather_bi,       ONLY: info_bi
    USE messy_cmn_photol_mem,        ONLY: IP_MAX, jname

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'dissoc_init_memory'
    INTEGER :: jt
    INTEGER :: status
!!$ CHARACTER(40) :: cr_date
!!$ CHARACTER(8)  :: ydate
!!$ CHARACTER(10) :: ytime, channel_name
    CHARACTER(LEN=STRLEN_OBJECT) :: obj_name = ''
#ifdef MBM_DISSOC
    INTEGER :: jx, jy, jz
    REAL(DP), dimension(10):: helparr
#endif

    !WRITE(*,*) substr

    CALL start_message_bi(modstr, 'MEMORY INITIALISATION', substr)

    CALL dissoc_initialize_dims

    CALL dissoc_initialize_reprs

    ! Note: channel object, in particular the look-up tables
    !       must be saved for restart ...
    CALL new_channel(status, modstr, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)

    ! jug 01/2014: define channel objects for DISSOC
    CALL new_channel_object(status, modstr, 'DISSOC_TABS', &
            p4=dissoc_tabs, reprid=REPR_DISSOC_TABS)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_TABS_DAVG', &
            p3=dissoc_tabs_davg, reprid=REPR_DISSOC_TABS_DAVG)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_SZA', &
            p1=dissoc_szagrid, reprid=REPR_DISSOC_SZA)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_WAVE', &
            p1=dissoc_lambda, reprid=REPR_DISSOC_WAVE)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_LAT', &
            p1=dissoc_latgrid, reprid=REPR_DISSOC_LAT)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_PRESS', &
            p1=dissoc_plevs, reprid=REPR_DISSOC_LEV)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_PRESC', &
            p1=dissoc_plevsc, reprid=REPR_DISSOC_LEVC)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_ALT', &
            p2=dissoc_alt, reprid=REPR_DISSOC_TWOD)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_ALTC', &
            p2=dissoc_altc, reprid=REPR_DISSOC_TWODC)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_TEMP', &
            p2=dissoc_temp, reprid=REPR_DISSOC_TWOD)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_TEMPC', &
            p2=dissoc_tempc, reprid=REPR_DISSOC_TWODC)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_DO2', &
            p2=dissoc_do2, reprid=REPR_DISSOC_TWOD)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_DO2C', &
            p2=dissoc_do2c, reprid=REPR_DISSOC_TWODC)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_DO3', &
            p2=dissoc_do3, reprid=REPR_DISSOC_TWOD)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DISSOC_DO3C', &
            p2=dissoc_do3c, reprid=REPR_DISSOC_TWODC)
    CALL channel_halt(substr, status)

    ! assign pointers in messy_dissoc to channel objects
    wavenm    => dissoc_lambda
    tabs      => dissoc_tabs
    tabs_davg => dissoc_tabs_davg
    lats      => dissoc_latgrid
    angdeg    => dissoc_szagrid
    pres      => dissoc_plevs
    presc     => dissoc_plevsc
    dtemp     => dissoc_temp
    dtempc    => dissoc_tempc
    alt       => dissoc_alt
    altc      => dissoc_altc
    do2       => dissoc_do2
    do2c      => dissoc_do2c
    do3       => dissoc_do3
    do3c      => dissoc_do3c

    CALL iniphoto(0)

#if defined(ECHAM5) || defined(COSMO) || defined(MBM_DISSOC)
    ! create grid point channel objects
    CALL new_channel_object(status, modstr, 'sza' &
         , p3=p_sza, reprid=GP_3D_MID )
    CALL new_attribute(status, modstr, 'sza' &
         , 'long_name', c='solar zenith angle' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sza' &
         , 'units', c='degrees' )
    CALL channel_halt(substr, status)

    ALLOCATE(jdiss_gp(IP_MAX))

    DO jt=1, IP_MAX
       IF (.NOT. jcalc(jt)) CYCLE
       !
       CALL new_channel_object(status, modstr &
            , 'J_'//TRIM(jname(jt)) &
            , p3=jdiss_gp(jt)%PTR, reprid=GP_3D_MID )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr  &
            , 'J_'//TRIM(jname(jt))             &
            , 'long_name', c='J('//TRIM(jname(jt))//')' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr  &
            , 'J_'//TRIM(jname(jt))             &
            , 'units', c='1/s')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'J_'//TRIM(jname(jt))      &
            , 'missing_value', r=mdi )
       CALL channel_halt(substr, status)
       CALL info_bi('channel / object '//modstr//' / J_'// &
            TRIM(jname(jt))//' was created')
    ENDDO
#endif

    ! check, if CLaMS is running ...
    ! NOTE: make sure that clams_init_memory is called BEFORE
    !       dissoc_init_memory in messy_main_control_e5.f90 ...
    CALL get_channel_info(status, 'clams')
    l_clams = (status == 0)
#ifdef MBM_CLAMS
    ! DISSOC in MBM CLAMS without CLaMS does not work ...
    CALL channel_halt(substr, status)
#endif

    ! IF CLaMS is running, create LG channel objects
    IF (l_clams) THEN

       ! ju_nt_20160217: new channel object SZA
       CALL new_channel_object(status, modstr, 'SZA' &
            , p1=SZA, reprid=REPR_LG_CLAMS)
       CALL new_attribute(status, modstr, 'SZA' &
            , 'long_name', c='solar zenith angle' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'SZA' &
            , 'units', c='degrees' )
       CALL channel_halt(substr, status)

       DO jt = 1, IP_MAX
          IF (.NOT. jcalc(jt)) CYCLE
          jresult(jt)%jname_messy = jname(jt)
          obj_name='J'//jname(jt)
          ! define channel objects for each photolysis rate
          CALL new_channel_object(status, modstr, trim(obj_name), &
               p1=jresult(jt)%values, reprid=REPR_LG_CLAMS)
          CALL channel_halt(substr, status)
! op_pj_20150709: obsolete: these are global attributes anyway
!!$       CALL new_attribute(status, modstr, trim(obj_name), &
!!$            'creator_of_parameter', c = trim(username))
!!$       CALL new_attribute(status, modstr, trim(obj_name), &
!!$            'param_creation_time',  c = trim(cr_date))
!!$       CALL new_attribute(status, modstr, trim(obj_name), &
!!$            'param_modification_time', c = trim(cr_date))
          CALL new_attribute(status, modstr, trim(obj_name), &
               'long_name', c= trim(obj_name))
          CALL new_attribute(status, modstr, trim(obj_name), &
               'units', c= 's^-1')
          CALL new_attribute(status, modstr, trim(obj_name), &
               'missing_value', r=mdi )
!!$       CALL new_attribute(status, modstr, trim(obj_name), &
!!$            'flag', c='NONE')
          CALL channel_halt(substr, status)
       ENDDO
    END IF

#ifdef MBM_DISSOC
    ! set dummy level data for 10 layers, pressure levels 1 to 1000 hPa
    CALL new_channel_object(status, modstr,  'tm1', &
         p3=tm1_3d, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1', &
         'long_name', c='dry air temperature (tm1)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1', 'units', c='K')
    CALL channel_halt(substr, status)
    helparr = (/271., 256., 239., 228., 223., 218., 217., 220., 249., 288./)
    DO jz=1, 10
       tm1_3d(:,:,jz) = helparr(jz)
    END DO

    CALL new_channel_object(status, modstr,  'tte', &
         p3=tte_3d, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tte', &
         'long_name', c='dry air temperature tendency (tte)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tte', 'units', c='K s-1')
    CALL channel_halt(substr, status)
    tte_3d(:,:,:) = 0.0_dp

    CALL new_channel_object(status, modstr,  'press', &
         p3=press_3d, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'press', 'long_name', c='pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'press', 'units', c='Pa')
    CALL channel_halt(substr, status)
    DO jz=1, 10
       press_3d(:,:,jz) = 100._dp * 10.**((jz-1)*1.0_dp/3.0_dp)
    END DO
    CALL new_channel_object(status, modstr,  'altitude_gnd', &
         p3=altitude_gnd, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'altitude_gnd', &
         'long_name', c='altitude above ground')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'altitude_gnd', 'units', c='m')
    CALL channel_halt(substr, status)
    helparr = (/48., 42., 36., 31., 26., 21., 16., 11., 6., 0./) *1000.0_dp
    DO jz=1, 10
       altitude_gnd(:,:,jz) = helparr(jz)
    END DO

    CALL new_channel_object(status, modstr,  'philat_2d', &
         p2=philat_2d, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'philat_2d' &
         , 'long_name', c='geographical latitude')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'philat_2d', 'units', c='degree')
    CALL channel_halt(substr, status)
    DO jx=1, nlon
       DO jy=1, nlat
          philat_2d(jx,jy) = -85.0_dp + REAL((jy-1)*10, DP)
       END DO
    END DO

    CALL new_channel_object(status, modstr,  'philon_2d', &
         p2=philon_2d, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'philon_2d' &
         , 'long_name', c='geographical longitude')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'philon_2d', 'units', c='degree')
    CALL channel_halt(substr, status)
    DO jy = 1, nlat
       DO jx=1, nlon
          philon_2d(jx,jy) = REAL((jx-1)*10, DP) + 5.0_dp
       END DO
    END DO

    CALL new_channel_object(status, modstr,  'alb', &
         p2=albedo_2d, reprid=GP_2D_HORIZONTAL)
    CALL new_attribute(status, modstr, 'alb', 'long_name', &
         c='surface background albedo')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'alb', 'units', c='-')
    CALL channel_halt(substr, status)
    albedo_2d(:,:) = 0.8_dp
#endif

    CALL end_message_bi(modstr, 'MEMORY INITIALISATION', substr)

  END SUBROUTINE dissoc_init_memory
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE dissoc_init_coupling

    USE messy_main_channel,          ONLY: get_channel_object, get_channel_info
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_blather_bi,       ONLY: info_bi
    USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM
    USE messy_main_grid_def_mem_bi,  ONLY: kproma, ngpblks
! op_pj_20150728: not yet used; required for proper ozone coupling, see below
!!$#if defined (ECHAM5) || defined (COSMO)
!!$    USE messy_main_tracer_mem_bi,  ONLY: GPTRSTR, gp_channel
!!$    USE messy_main_channel_tracer, ONLY: set_channel_or_tracer
!!$#endif

! op_pj_20150728: CLaMS dependency, needs to be resolved by coupling
!                 via channel objects ...
!!$    USE messy_clams_global,        ONLY: initfile
!!$    USE messy_clams_tools_utils,   ONLY: uppercase
!!$    USE messy_clams_tools_ncutils, ONLY: nc_get_vertcoorname

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'dissoc_init_coupling'
    INTEGER :: status

! op_pj_20150728: CLaMS dependency, needs to be resolved by coupling
!                 via channel objects ...
!!$    CHARACTER(30) :: vertcoorname
!!$    CALL nc_get_vertcoorname (trim(initfile), vertcoorname)

! op_pj_20150727: later, in case the look-up table calculation is more
!                 flexible, O3, etc. should be coupled via namelist-entry ...
!!$    ! CHECK IF OZONE IS AVAILABLE
!!$    CALL start_message_bi(modstr, 'COUPLING TO Ozone', substr)
!!$    CALL info_bi('Looking for OZONE ... ')
!!$    CALL info_bi('       channel: '//dissoc_o3%cha)
!!$    CALL info_bi('       object : '//dissoc_o3%obj)
!!$
!!$    CALL set_channel_or_tracer(status, GPTRSTR, gp_channel &
!!$         , dissoc_o3%cha, dissoc_o3%obj, ptr_O3, ptr_O3te)
!!$    CALL channel_halt(substr, status)
!!$    CALL end_message_bi(modstr, 'COUPLING TO Ozone', substr)

    !WRITE(*,*) substr

!!$! op_pj_20160617+
!!$#ifndef MBM_DISSOC
!!$    CALL get_channel_object(status, 'rad', 'albedo', p2=albedo_2d)
!!$    CALL channel_halt(substr, status)
!!$#endif
!!$! op_pj_20160617-
#ifndef MBM_DISSOC
    CALL get_channel_info(status, 'rad')
    if (status == 0) then
       CALL get_channel_object(status, 'rad', 'albedo', p2=albedo_2d)
       CALL channel_halt(substr, status)
    else
       allocate (albedo_2d(kproma,ngpblks))
       albedo_2d(:,:) = 0.8_dp
    endif
#endif


    IF (l_clams) THEN
       CALL start_message_bi(modstr, 'COUPLING TO CLaMS DRIVER FIELDS', substr)

       ! ju_nt_20160217: couple dnparts and dnparts_max
       CALL get_channel_object(status, 'clams', 'dnparts', p1=dnparts_co)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'clams', 'dnparts_max', p1=dnparts_max_co)
       CALL channel_halt(substr, status)


       CALL get_channel_object(status, 'clams', 'LAT_OLD', p1=LAT)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'clams', 'LON_OLD', p1=LON)
       CALL channel_halt(substr, status)

       CALL get_channel_object(status, 'clams', 'TEMP_OLD', p1=TEMP)
       IF (status /= 0) CALL error_bi &
            ('Channelobject TEMP_OLD does not exist -> check, if CLAMSCHEM is switched on ! ',substr)
       CALL channel_halt(substr, status)

! ju_nt_20160209: only if PRESS is vertical coordinate,
!                 it is stored on channelobject LEV
!!$       CALL get_channel_object(status, 'clams', 'LEV', p1=PRESS)
!!$       IF (status /= 0) THEN
! ju_nt_20160209: if vertical coordinate is THETA/ZETA:
          CALL get_channel_object(status, 'clams', 'PRESS_OLD', p1=PRESS)
          IF (status /= 0) CALL error_bi &
               ('Channelobject PRESS_OLD does not exist -> check, if CLAMSCHEM is switched on ! ',substr)
!!$       END IF
       CALL channel_halt(substr, status)

       ! requires? better use central TIMER information
       CALL get_channel_object(status, 'clams', 'JULTIME_OLD', p0=JULTIME)
       CALL channel_halt(substr, status)

       CALL end_message_bi(modstr, 'COUPLING TO CLaMS DRIVER FIELDS', substr)
    END IF

  END SUBROUTINE dissoc_init_coupling
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE dissoc_global_start

    USE messy_main_timer_bi,  ONLY: event_state
    USE messy_main_timer,     ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND &
                                   , time_step_len, lstart, lresume &
                                   , current_date
    USE messy_main_blather_bi, ONLY: info_bi
    USE messy_main_timer,      ONLY: gregor2julian, julian_day, delta_time
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_bcast, p_io, p_pe
    USE messy_main_constants_mem, ONLY: STRLEN_VLONG

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'dissoc_global_start'
    INTEGER :: jday
    REAL(DP):: js, jd
    CHARACTER(LEN=STRLEN_VLONG)  :: str = ''

    !WRITE(*,*) substr

    ! ju_nt_20180620+
    if (timestep_dissoc == delta_time) then
       ldissocevent = .TRUE.
    else
       ldissocevent = event_state(dissocevent, current_date)
    endif
    if (timestep_dissoc_gp == delta_time) then
       ldissocevent_gp = .TRUE.
    else
       ldissocevent_gp = event_state(dissocevent_gp, current_date)
    endif
    ! ju_nt_20180620-

    ! TRIGGER CLIMATOLOGY UPDATE?
    LTRIG_MONTHLY = event_state(XTIMER_MONTHLY, current_date)

    if (p_pe==0) write (*,*) substr,': lstart,lresume,ltrig_monthly=',lstart,lresume,ltrig_monthly

    update: IF (lstart .OR. lresume .OR. LTRIG_MONTHLY) THEN

       ! UPDATE TIME INFORMATION
       jdfrac = gregor2julian(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND) &
            - julian_day(REAL(1, DP), 1, 2000) & ! alternative t-zero ???
            - 0.5_dp                             ! 00:00:00 as reference ....
       jsec   = jdfrac*86400._dp

       jday = 0 !INT(jdfrac)
       CALL iniphoto(jday)

       ! CLEAN MEMORY FIRST
       !IF (ALLOCATED(pres_inp)) DEALLOCATE(pres_inp)
       !IF (ALLOCATED(alt_inp))  DEALLOCATE(alt_inp)
       !IF (ALLOCATED(temp_inp)) DEALLOCATE(temp_inp)
       !IF (ALLOCATED(o3_inp))   DEALLOCATE(o3_inp)

       IF (p_parallel_io) THEN
          ! read ozone climatology on I/O pe and broadcast results to all
          ! othters
          call reado3(o3dat,MONTH,nlats)
       ENDIF
       CALL p_bcast(nlats, p_io)
       ! ju_nt_20170217: broadcast lats
       CALL p_bcast(lats, p_io)
       CALL p_bcast(dim_pres_inp, p_io)
       IF (.NOT. p_parallel_io) THEN
          ALLOCATE(pres_inp(dim_pres_inp))
          ALLOCATE(alt_inp(dim_pres_inp,nlats))
          ALLOCATE(temp_inp(dim_pres_inp,nlats))
          ALLOCATE(o3_inp(dim_pres_inp,nlats))
       ENDIF
       CALL p_bcast(pres_inp, p_io)
       CALL p_bcast(alt_inp, p_io)
       CALL p_bcast(temp_inp, p_io)
       CALL p_bcast(o3_inp, p_io)

       CALL p_ascend (dim_pres_inp, nlats, pres_inp, alt_inp, temp_inp, o3_inp)
       CALL setp
       DEALLOCATE(pres_inp)
       DEALLOCATE(alt_inp)
       DEALLOCATE(temp_inp)
       DEALLOCATE(o3_inp)
       CALL settab

       if (davg) then
!!$    js = ymds2js_interface(YEAR, MONTH, 15, 43200, .FALSE.)
! op_pj_20150709: this needs to be checked for consistency
          jd = julian_day(REAL(15, DP), MONTH, YEAR) - &
               julian_day(REAL(1, DP), 1, 2000)   ! alternative T-zero
          ! ju_nt_20160217: js changed (julian seconds at noon)
!!$       js = jd*86400._dp
          js = jd*86400._dp + 43200._dp
          WRITE(str,*) js
          CALL info_bi('init davg photolysis, js = '//TRIM(str), substr)
          ! ju_nt_20160217: set jd for calc_davg
          jd = gregor2julian(YEAR, MONTH, 15, 12, 0, 0)
          if (p_pe==0) write(*,*) substr,' init davg photolysis: js,jd', js,jd
          call calc_davg(js, jd)
       endif

    END IF update

    ! ju_nt_20160217: set dnparts and dnparts_max
    if (l_clams) then
       dnparts = dnparts_co(1)
       dnparts_max = dnparts_max_co(1)
       !nparts  = p_sum(dnparts)
    endif

  END SUBROUTINE dissoc_global_start
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE dissoc_physc

    ! NOTE: This routine can solely be used for the grid-point models!
    !       Lagrangian calculations can only be performed in messy_global_end.

#if defined (ECHAM5) || defined (COSMO) || defined (MBM_DISSOC)

    USE messy_main_timer,         ONLY: time_step_len
    USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, kproma, jrow
#ifndef MBM_DISSOC
    USE messy_main_grid_def_bi,   ONLY: philon_2d, philat_2d, altitude_gnd
    USE messy_main_data_bi,       ONLY: tm1_3d, tte_3d &
                                      , press_3d !!$, albedo_2d=>albedo ! op_pj_20160617
#endif
    USE messy_main_constants_mem, ONLY: g

    IMPLICIT NONE
    INTRINSIC :: MAX, ASSOCIATED

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'dissoc_physc'
    INTEGER :: jk, jp, jt
    REAL(dp) :: temp_2d(nproma,nlev)
    REAL(kind=DP), DIMENSION(IP_MAX) :: photoarr

    !WRITE(*,*) substr

    ! ju_nt_20180620: ldissocevent_gp added
    if (ldissocevent_gp) then

       temp_2d(:,:) = tm1_3d(_RI_XYZ__(:,jrow,:))  + tte_3d(_RI_XYZ__(:,jrow,:))  * time_step_len

       levels: DO jk=1,nlev

       ! ju_nt_20160217: julian day and julian seconds correct ???
!!!!! jdfrac, jsec in richtiger Reihenfolge ?!
!!!!! jdfrac statt jd ?
          CALL calc_zenith(philat_2d(1:kproma,jrow) &
                         , philon_2d(1:kproma,jrow) &
                         , kproma &
                         , jdfrac, jsec &
                         , p_sza(_RI_XYZ__(1:kproma,jrow,jk)) &
                         , altitude_gnd(_RI_XYZ__(1:kproma,jrow,jk))/1000.0_dp &  ! [km]
                         )

          vector: DO jp=1, kproma
             photoarr(:) = 0.0_dp
             albedo = albedo_2d(jp,jrow)
             CALL dissoc(temp_2d(jp,jk) &
                       , press_3d(_RI_XYZ__(jp,jrow,jk))/100.0_dp  &    ! Pa -> hPa
                       , philat_2d(jp,jrow)             &
                       , p_sza(_RI_XYZ__(jp,jrow,jk))*dtr          &    ! degree -> radian
                       , photoarr(:) &
                       )

             DO jt=1, IP_MAX
                IF (jcalc(jt)) THEN
                   jdiss_gp(jt)%PTR(_RI_XYZ__(jp,jrow,jk)) = photoarr(jt)
                END IF
             END DO
          END DO vector

       END DO levels

    endif

#endif

#ifdef MBM_DISSOC

#endif

!!$    USE messy_main_constants_mem, ONLY: pi, dtr
!!$    USE messy_main_tools,         ONLY: PTR_2D_ARRAY
!!$    USE messy_dissoc_global
!!$    USE messy_dissoc_si
!!$    USE messy_cmn_photol_mem
!!$    USE mo_netcdf,                ONLY: open_jval_nc_file,  &
!!$                                        write_jval_nc_file, &
!!$                                        close_jval_nc_file
!!$
!!$    IMPLICIT NONE
!!$
!!$    INTEGER, PARAMETER :: nsza = 10
!!$    INTEGER, PARAMETER :: nlev = 19 ! number of levels
!!$    REAL, DIMENSION(nsza) :: sza = &
!!$         (/ 0., 20., 40., 60., 70., 80., 85., 90., 92., 94./)
!!$
!!$    INTEGER :: i, j, ispec
!!$
!!$    !  njs: number of calculated photolysis rates
!!$    integer, parameter :: njs = 39
!!$    character(10), dimension(njs) :: photonames
!!$    integer, dimension(njs) :: ip_specs
!!$    INTEGER :: ncid, jt,status
!!$
!!$    REAL(kind=DP), DIMENSION(nlev)    :: press_calc
!!$    REAL(kind=DP), DIMENSION(nlev)    :: temp_calc
!!$    REAL(kind=DP), DIMENSION(IP_MAX)  :: photoarr
!!$
!!$    TYPE(PTR_2D_ARRAY), DIMENSION(IP_MAX) :: jarray
!!$
!!$    REAL(kind=DP) :: lat_calc = 17.5000
!!$
!!$    DO ispec = 1, IP_MAX
!!$       allocate(jarray(ispec)%ptr(nsza,nlev))
!!$    END DO
!!$
!!$    albedo  = 0.07
!!$
!!$    davg=.false.
!!$
!!$    press_calc = (/ &
!!$         1000.,     1292.,     1668.,     2154.,     2783.,&
!!$         3594.,     4642.,     5995.,     7743.,    10000.,&
!!$         12915.,    16681.,    21544.,    27826.,    35938.,&
!!$         46416.,    59948.,    77426.,    90000./)
!!$
!!$    temp_calc = (/ &
!!$         228.36,    224.99,    221.91,    219.11,    215.84,&
!!$         211.72,    206.74,    200.57,    197.04,    196.14,&
!!$         203.06,    211.95,    222.80,    235.88,    249.06,&
!!$         262.34,    274.18,    285.72,    292.35/)
!!$
!!$
!!$    do i=1,nsza
!!$      do j=1,nlev
!!$     CALL dissoc(temp_calc(j), press_calc(j)/100., lat_calc, sza(i)*dtr, photoarr)
!!$     do ispec=1, IP_MAX
!!$           if (jcalc(ispec)) jarray(ispec)%ptr(i,j) = photoarr(ispec)
!!$     enddo
!!$      enddo
!!$    enddo
!!$
!!$
!!$    ! output to netcdf file
!!$    CALL open_jval_nc_file(ncid, nlev, nsza, sza)
!!$
!!$    photonames = (/ 'jO2', 'jO3', 'jO3a', 'jH2O2', 'jCl2', 'jCl2O2', &
!!$         'jHOCl', 'jClNO2', 'jClONO2', 'jHNO3', 'jNO2', 'jN2O5', &
!!$         'jHO2NO2a', 'jNO3a', 'jNO3b', 'jCH2Ob', 'jCH2Oa', 'jCH3O2H', &
!!$         'jBrNO3a', 'jBrCl', 'jOClO', 'jH2O', 'jHCl', 'jNO', &
!!$         'jN2O', 'jCH3OCl', 'jHOBr', 'jBr2', 'jMEO2NO2', 'jBrO', &
!!$         'jClONO2a','jBrNO3b', 'jHO2NO2b', 'jF11', 'jF12', 'jF22', &
!!$         'jF113', 'jCH3Cl','jCCl4' /)
!!$
!!$    ip_specs = (/ip_O2, ip_O3P, ip_O1D, ip_H2O2, ip_Cl2, ip_Cl2O2, &
!!$         ip_HOCl, ip_ClNO2, ip_ClNO3, ip_HNO3, ip_NO2, ip_N2O5, &
!!$         ip_HO2NO2, ip_NO2O, ip_NOO2, ip_CHOH, ip_COH2, ip_CH3OOH, &
!!$         ip_BrNO3, ip_BrCl, ip_OClO, ip_H2O, ip_HCl, ip_NO, &
!!$         ip_N2O, ip_CH3OCl, ip_HOBr, ip_Br2, ip_MEO2NO2, ip_BrO, &
!!$         ip_ClONO2, ip_BrONO2, ip_OHNO3, ip_CFCl3, ip_CF2CL2, ip_CHF2Cl, &
!!$         ip_F113, ip_CH3Cl, ip_CCl4/)
!!$
!!$    DO jt = 1, njs
!!$
!!$       IF (jcalc(ip_specs(jt))) THEN
!!$          CALL write_jval_nc_file(ncid, TRIM(photonames(jt)), &
!!$            jarray(ip_specs(jt))%ptr(:,:))
!!$       ENDIF
!!$    ENDDO
!!$
!!$
!!$    CALL close_jval_nc_file(ncid)
!!$
!!$
!!$    DO ispec = 1, IP_MAX
!!$      deallocate(jarray(ispec)%ptr)
!!$    END DO

  END SUBROUTINE dissoc_physc
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE dissoc_global_end

    USE messy_main_timer,         ONLY: YEAR_NEXT, MONTH_NEXT, DAY_NEXT, &
                                        HOUR_NEXT, MINUTE_NEXT, SECOND_NEXT
    USE messy_main_mpi_bi,        ONLY: p_pe

    ! SMCL: MESSy
    USE messy_main_constants_mem, ONLY: dtr
    USE messy_cmn_photol_mem,     ONLY: IP_MAX
    USE messy_main_timer,         ONLY: gregor2julian

    IMPLICIT NONE

    integer                     :: i, irate
    CHARACTER(LEN=*), PARAMETER :: substr = 'dissoc_global_end'
    REAL(DP), DIMENSION(IP_MAX) :: photoarr
    REAL(DP)                    :: jd


    !WRITE(*,*) substr

    ! ju_nt_20180620: ldissocevent added
    if (l_clams .and. ldissocevent) then

       if (p_pe==0) write(*,*) substr

       jd = gregor2julian (YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT, &
                           MINUTE_NEXT, SECOND_NEXT)
       if (p_pe==0) write(*,*) 'call calc_zenith: jultime, jd=', jultime, jd
       call calc_zenith(LAT, LON, dnparts, JULTIME, jd, SZA)


       ! Loop over used air parcels
       do i=1,dnparts

          if(ABS(LAT(i)) < 100.  .and. press(i) > 0.) then

             CALL dissoc(TEMP(i), PRESS(i), LAT(i), sza(i)*dtr, photoarr )

             do irate=1,IP_MAX
                if (jcalc(irate)) jresult(irate)%values(i) = photoarr(irate)
             enddo
!!$             dissoc_rate(i,:) = pack(photoarr, jcalc)
          else
             ! set the photolysis rates of not valid APs to 0.
             do irate=1,iP_MAX
                if (jcalc(irate)) jresult(irate)%values(i) = 0.0
             enddo
!!$             dissoc_rate(i,:) = 0.0
          endif

       enddo

       !  Loop over invalid air parcels
       do irate=1,iP_MAX
          if (jcalc(irate)) jresult(irate)%values(dnparts+1:dnparts_max) = mdi
       enddo
!!$        dissoc_rate(dnparts+1:dnparts_max,:) = mdi


    endif

  END SUBROUTINE dissoc_global_end
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE dissoc_free_memory

    IMPLICIT NONE

#if defined (ECHAM5) || defined (COSMO)
    DEALLOCATE(jdiss_gp)
#endif

  END SUBROUTINE dissoc_free_memory
  !--------------------------------------------------------------------

! op_pj_20150728: currently not used; later required for proper coupling
!                 to external ozone fields
!!$  ! ========================================================================
!!$  SUBROUTINE dissoc_read_nml_cpl(status, iou)
!!$
!!$    ! DISSOC MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
!!$    !
!!$    ! read namelist for 'coupling' to ECHAM5
!!$    !
!!$
!!$    ! MESSy
!!$    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
!!$
!!$    IMPLICIT NONE
!!$    INTRINSIC :: TRIM
!!$
!!$    ! I/O
!!$    INTEGER, INTENT(OUT) :: status     ! error status
!!$    INTEGER, INTENT(IN)  :: iou        ! I/O unit
!!$
!!$    NAMELIST /CPL/ dissoc_o3
!!$
!!$    ! LOCAL
!!$    CHARACTER(LEN=*), PARAMETER :: substr='dissoc_read_nml_cpl'
!!$    LOGICAL              :: lex      ! file exists ?
!!$    INTEGER              :: fstat    ! file status
!!$
!!$    status = 1
!!$
!!$    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
!!$    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist
!!$
!!$    READ(iou, NML=CPL, IOSTAT=fstat)
!!$    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
!!$    IF (fstat /= 0) RETURN  ! error while reading namelist
!!$
!!$    CALL read_nml_close(substr, iou, modstr)
!!$    status = 0 ! NO ERROR
!!$
!!$  END SUBROUTINE dissoc_read_nml_cpl
!!$  ! ========================================================================

  !--------------------------------------------------------------------
  SUBROUTINE dissoc_initialize_dims

    USE messy_main_channel_dimensions,   ONLY: new_dimension            &
                                             , write_dimension          &
                                             , add_dimension_variable   &
                                             , add_dimension_variable_att
    USE messy_main_channel_error_bi,     ONLY: channel_halt

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'dissoc_initialize_dims'
    INTEGER                     :: status

    CALL new_dimension(status, DIMID_WAVE, 'WAVELENGTH', jpwave)
    CALL channel_halt(substr, status)

    CALL new_dimension(status, DIMID_DISSOC_LEV, 'DISSOC_LEV', jpslevall)
    CALL channel_halt(substr, status)

    CALL new_dimension(status, DIMID_DISSOC_LEVC, 'DISSOC_LEVC', jpslev)
    CALL channel_halt(substr, status)

    CALL new_dimension(status, DIMID_DISSOC_LAT, 'DISSOC_LAT', jplats)
    CALL channel_halt(substr, status)

    CALL new_dimension(status, DIMID_DISSOC_SZA, 'DISSOC_SZA', jpschi)
    CALL channel_halt(substr, status)

  END SUBROUTINE dissoc_initialize_dims
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE dissoc_initialize_reprs

    USE messy_main_channel_error_bi,     ONLY: channel_halt
    USE messy_main_channel_bi,           ONLY: DC_IX
    USE messy_main_channel_repr,         ONLY: new_representation &
         , set_representation_decomp, write_representation &
         , AUTO, IRANK               &
         , PIOTYPE_SGL, PIOTYPE_IND  &
         , PIOTYPE_COL

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'dissoc_initialize_reprs'
    INTEGER                        :: status
    INTEGER, DIMENSION(:), POINTER :: scdim => NULL()

! jug: DISSOC representations
    CALL new_representation(status, REPR_DISSOC_WAVE, 'REPR_DISSOC_WAVE'  &
         , rank = 1, link = 'x---', dctype = 0                           &
         , dimension_ids = (/ DIMID_WAVE /)                              &
         , ldimlen       = (/ jpwave /)                                    &
         , axis = 'N---'                                                 &
         )
    CALL channel_halt(substr, status)

!?? jug level edges and level centers
    CALL new_representation(status, REPR_DISSOC_LEV, 'REPR_DISSOC_LEV'  &
         , rank = 1, link = 'x---', dctype = 0                           &
         , dimension_ids = (/ DIMID_DISSOC_LEV /)                        &
         , ldimlen       = (/ jpslevall /)                                &
         , axis = 'Z---'                                                 &
         )
    CALL channel_halt(substr, status)

    CALL new_representation(status, REPR_DISSOC_LEVC, 'REPR_DISSOC_LEVC'  &
         , rank = 1, link = 'x---', dctype = 0                            &
         , dimension_ids = (/ DIMID_DISSOC_LEVC /)                        &
         , ldimlen       = (/ jpslev /)                                   &
         , axis = 'Z---'                                                  &
         )
    CALL channel_halt(substr, status)

    CALL new_representation(status, REPR_DISSOC_SZA, 'REPR_DISSOC_SZA'  &
         , rank = 1, link = 'x---', dctype = 0                           &
         , dimension_ids = (/ DIMID_DISSOC_SZA /)                        &
         , ldimlen       = (/ jpschi /)                                  &
         , axis = 'x---'                                                 &
         )
    CALL channel_halt(substr, status)

    CALL new_representation(status, REPR_DISSOC_LAT, 'REPR_DISSOC_LAT'  &
         , rank = 1, link = 'x---', dctype = 0                           &
         , dimension_ids = (/ DIMID_DISSOC_LAT /)                       &
         , ldimlen       = (/ jplats /)                                  &
         , axis = 'Y---'                                                 &
         )
    CALL channel_halt(substr, status)

    CALL new_representation(status, REPR_DISSOC_TABS, 'REPR_DISSOC_TABS' &
         , rank = 4, link = 'xxxx', dctype = 0                           &
         , dimension_ids = (/ DIMID_DISSOC_LEVC, DIMID_DISSOC_SZA,        &
                              DIMID_WAVE, DIMID_DISSOC_LAT /)            &
         , ldimlen = (/ jpslev, jpschi, jpwave, jplats /)                &
         , axis = 'ZxxY'                                                 &
         )
    CALL channel_halt(substr, status)

    CALL new_representation(status, REPR_DISSOC_TABS_DAVG, 'REPR_DISSOC_TABS_DAVG'    &
         , rank = 3, link = 'xxx-', dctype = 0                           &
         , dimension_ids = (/ DIMID_DISSOC_LEVC, DIMID_WAVE, DIMID_DISSOC_LAT /) &
         , ldimlen = (/ jpslev, jpwave, jplats /)                        &
         , axis = 'ZxY-'                                                 &
         )
    CALL channel_halt(substr, status)

    CALL new_representation(status, REPR_DISSOC_TWOD, 'REPR_DISSOC_TWOD' &
         , rank = 2, link = 'xx--', dctype = 0                           &
         , dimension_ids = (/ DIMID_DISSOC_LEV, DIMID_DISSOC_LAT /)      &
         , ldimlen = (/ jpslevall, jplats /)                             &
         , axis = 'ZY--'                                                 &
         )
    CALL channel_halt(substr, status)

    CALL new_representation(status, REPR_DISSOC_TWODC, 'REPR_DISSOC_TWODC' &
         , rank = 2, link = 'xx--', dctype = 0                           &
         , dimension_ids = (/ DIMID_DISSOC_LEVC, DIMID_DISSOC_LAT /)     &
         , ldimlen = (/ jpslev, jplats /)                                &
         , axis = 'ZY--'                                                 &
         )
    CALL channel_halt(substr, status)

  END SUBROUTINE dissoc_initialize_reprs
  !--------------------------------------------------------------------

!**********************************************************************
END MODULE messy_dissoc_si
!**********************************************************************
