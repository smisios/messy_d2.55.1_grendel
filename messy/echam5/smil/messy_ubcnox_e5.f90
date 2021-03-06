!#include "messy_main_ppd_bi.inc"

! DESCRIPTION:
! This module is the interface (called by messy_main_control) between the
! ECHAM5 base model and the ubcnox module.
!
! AUTHOR:
! S. Versick, Steinbuch Centre for Computing, KIT, Germany
!             Institute for Meteorology and Climate Research, KIT, Germany
!   (2013: Creation of this submodel) (SV)
!      Version 1.0: only top level (SV)
!   (2015: Code improvements (SV)
!      Version 1.1: ubcnox-flux added (SV)
!      Version 1.2: production rate from Holger Nieder added (SV)
!      Version 1.3: multi level support for amount added (SV)
!                   some code cleanup (SV)
!                   set NO2 = 0 for amountbased parameterizations (SV)
!      Version 1.4: added support for elevated stratopause events (SV) -> obsolete
!   (2016)
!      Version 1.5: online calculation for Funke et al 2016 added (SV)
!      Version 1.5.1: bugfix for ESEs in (close to) May (SV)
!---------------------------------------------------------------------------

MODULE messy_ubcnox_e5

  ! ECHAM5/MESSy
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, dp
  USE messy_main_mpi_bi,       ONLY: message
  USE messy_main_channel,       ONLY: t_chaobj_cpl
  USE messy_main_timer_event,  ONLY: time_event, io_time_event ! op_pj_20190410
    ! MESSy
  USE messy_main_grid_def_bi,  ONLY: grvol

  IMPLICIT NONE
  SAVE

  !-----------------------------------------------------------------
  ! Everything is PRIVATE, except when explicitely stated otherwise:
  !-----------------------------------------------------------------
  PRIVATE

  INTRINSIC NULL

  ! op_pj_20190410+
  ! ub_ak_20190709 SAVE deleted , global SAVE defined (g95 ERROR)
  TYPE(io_time_event) :: TIMER_DAILY = &
       io_time_event(1, 'days','first',0)
  TYPE(time_event) :: XTIMER_DAILY
  ! op_pj_20190410-

  REAL(dp), DIMENSION(:,:,:), POINTER :: &
    noxubc => NULL()

  REAL(dp), DIMENSION(:,:,:), POINTER :: &
!    debug  => NULL(), &
    tend   => NULL(), &
    ep     => NULL(), &
    bg     => NULL(), &
    ubc    => NULL(), &
    ese    => NULL()

! ka_sv_20170523+
  REAL(dp), POINTER :: daysofese, daysofesecalc  ! real number of days of elevated stratopause and number that goes into calculation (effects persist longer than ese)
! ka_sv_20170523-


  INTEGER :: idt_PTNOX1, idt_NOX, idt_NO2
  INTEGER, DIMENSION(:), POINTER     :: idt_list_n
  INTEGER                            :: nlntrac_n    ! no of tracers
  !-----
  ! global variable for regrid events
  !-----
!!$  TYPE(RGTEVENT), DIMENSION(:), POINTER :: rgt
! REAL(dp), POINTER :: grvol(:,:,:)   => NULL()

  !-----
  ! switch for tracer initialisation
  !-----
  LOGICAL :: tracer_init_required = .false.

  !-----
  ! 2d-field for tropopause index (details via namelist)
  !-----

  REAL(dp), DIMENSION(:,:,:), POINTER :: ptnox_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: nox_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: tm1_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: tm1_save_ptr => NULL() ! op_pj_20190410
  REAL(dp), DIMENSION(:), POINTER :: kp_data => NULL()
  REAL(dp), DIMENSION(:), POINTER :: f107_data => NULL()
  REAL(dp), DIMENSION(:), POINTER :: ap_data => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: lat_ptr  => NULL()

  CHARACTER (LEN=STRLEN_MEDIUM), PUBLIC        ::  &
       nox_channel= '', &
       nox_object='', &
       target_species=''

  integer, parameter :: nkp = 14, nf107 = 3, nmonth=12, nlat=18
  REAL(dp) :: noxlut(nkp, nmonth, nf107, nlat)   ! NOX Look-Up Table by Holger Nieder
  REAL(dp):: ubc_kp(nkp)
  REAL(dp) :: ubc_f107(nf107)
  REAL(dp) :: ubc_lat_bin(nlat+1)

  REAL(dp), DIMENSION(:,:,:), POINTER :: fac => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: face => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: facdec => NULL() ! decomposed field
  REAL(dp), DIMENSION(:,:,:), POINTER :: facedec => NULL() ! decomposed field

    ! CPL namelist
  ! - name of tracer/channel object
  TYPE(t_chaobj_cpl) :: ubcnox_kp, ubcnox_f107, ubcnox_ap

  PUBLIC :: ubcnox_initialize
  PUBLIC :: ubcnox_new_tracer
  PUBLIC :: ubcnox_init_memory
  PUBLIC :: ubcnox_init_coupling
  PUBLIC :: ubcnox_init_tracer
  PUBLIC :: ubcnox_physc
  PUBLIC :: ubcnox_global_start
  PUBLIC :: ubcnox_free_memory  ! op_pj_20161206

  PUBLIC :: daysofese, daysofesecalc  ! ka_sv_20170523
  PUBLIC :: fac, facdec, face, facedec

  CHARACTER(len=*), PARAMETER :: modstr = 'ubcnox'

CONTAINS

!=============================================================================

  SUBROUTINE ubcnox_initialize

    ! ECHAM5/MESSy
    USE messy_main_tools ,    ONLY: find_next_free_unit
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast, finish
    USE messy_main_timer_bi,   ONLY: timer_event_init ! op_pj_20190410
    USE messy_ubcnox,          ONLY: ubcnox_read_nml_ctrl, nox_switch &
                        , ubcnox_read_noxlut, ubcnox_provide_data, &
                          ese_switch, top_levels, tempthresh, hn_file

    ! MESSy
    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: substr = 'ubcnox_initialize'
    INTEGER                     :: status, iou


    !-----
    ! read namelist CTRL
    !-----
      IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL ubcnox_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL finish(substr,'call ubcnox_read_nml_ctrl failed')
     END IF
     CALL p_bcast(nox_switch, p_io)
     CALL p_bcast(ese_switch, p_io)
    !-----
    ! read namelist CPL
    !-----
    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      CALL ubcnox_read_nml_cpl(status, iou)   ! read /CPL/
      IF (status /= 0) CALL finish(substr,'call ubcnox_read_nml_cpl failed')
    END IF
    CALL p_bcast (nox_channel,      p_io)
    CALL p_bcast (nox_object,       p_io)
    CALL p_bcast (target_species,   p_io)

    IF (nox_switch==3) THEN
        CALL p_bcast(hn_file, p_io)
        CALL ubcnox_read_noxlut(noxlut)
        CALL p_bcast(noxlut, p_io)

        CALL ubcnox_provide_data(ubc_lat_bin,ubc_kp,ubc_f107)
        CALL p_bcast(ubc_lat_bin, p_io)
        CALL p_bcast(ubc_kp, p_io)
        CALL p_bcast(ubc_f107, p_io)

        CALL p_bcast(ubcnox_kp%cha, p_io)
        CALL p_bcast(ubcnox_kp%obj, p_io)
        CALL p_bcast(ubcnox_f107%cha, p_io)
        CALL p_bcast(ubcnox_f107%obj, p_io)

    END IF

    IF (nox_switch==5) THEN
        CALL p_bcast(ubcnox_ap%cha, p_io)
        CALL p_bcast(ubcnox_ap%obj, p_io)
        CALL p_bcast(top_levels, p_io)
        CALL p_bcast(tempthresh, p_io)
    END IF

    ! op_pj_20190410+
    CALL timer_event_init(XTIMER_DAILY,   TIMER_DAILY &
                         ,   'ubcnox_daily',   'present')
    ! op_pj_20190410-

  END SUBROUTINE ubcnox_initialize

!=============================================================================

   SUBROUTINE ubcnox_new_tracer
! ---------------------------------------------------------------
     USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
     USE messy_main_tracer_tools_bi, ONLY: tracer_halt
     USE messy_main_blather_bi,      ONLY: start_message_bi
     USE messy_main_tracer,          ONLY: get_tracer

     IMPLICIT NONE

     CHARACTER(LEN=*), PARAMETER :: substr='ubcnox_new_tracer'
     INTEGER :: i_err
     INTEGER :: status
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: subnames => NULL()


   END SUBROUTINE ubcnox_new_tracer

!=============================================================================

   SUBROUTINE ubcnox_init_memory
! ---------------------------------------------------------------

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, SCALAR
    ! MESSy
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute
     IMPLICIT NONE
     INTEGER :: status
     CHARACTER(LEN=*), PARAMETER           :: substr = 'ubcnox_init_memory'


     CALL message('ubcnox_init_memory','defining streams for ubcnox')

     ! -----
     ! ubcnox specific molecule information
     ! -----
     CALL new_channel(status, modstr, reprid=GP_3D_MID)
     CALL channel_halt(substr, status)

!     CALL new_channel_object(status, modstr, 'debug' &
!          , p3 = debug)
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'debug' &
!          , 'long_name', c='Just for debugging')
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'debug' &
!          , 'units', c='debug')
!     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'tend' &
          , p3 = tend)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'tend' &
          , 'long_name', c='Tendency from UBCNOX')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'tend' &
          , 'units', c='mol/mol/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'ep' &
          , p3 = ep)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'ep' &
          , 'long_name', c='EPP-noy field (without background)')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'ep' &
          , 'units', c='mol/mol')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'bg' &
          , p3 = bg)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'bg' &
          , 'long_name', c='EPP-noy field (only background)')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'bg' &
          , 'units', c='mol/mol')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'ubc' &
          , p3 = ubc)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'ubc' &
          , 'long_name', c='EPP-noy field (complete)')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'ubc' &
          , 'units', c='mol/mol')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'ese' &
          , p3 = ese)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'ese' &
          , 'long_name', c='EPP-noy field (only ESE)')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'ese' &
          , 'units', c='mol/mol')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'daysofese' &
            , p0=daysofese &
            , reprid=SCALAR &
            , lrestreq = .TRUE. )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'daysofese' &
          , 'long_name', c='Number of days since start of ESE')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'daysofese' &
          , 'units', c='days')
     CALL channel_halt(substr, status)

! ka_sv_201705023+
      CALL new_channel_object(status, modstr, 'daysofesecalc' &
            , p0=daysofesecalc &
            , reprid=SCALAR &
            , lrestreq = .TRUE. )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'daysofesecalc' &
          , 'long_name', c='Number of days since start of ESE')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'daysofesecalc' &
          , 'units', c='days')
     CALL channel_halt(substr, status)
! ka_sv_201705023-

     ! op_pj_20190410+
     CALL new_channel_object(status, modstr, 'tm1s' &
            , p3=tm1_save_ptr &
            , reprid=GP_3D_MID &
            , lrestreq = .TRUE. )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'tm1s' &
          , 'long_name', c='temperature at first time step of day')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'tm1s' &
          , 'units', c='K')
     CALL channel_halt(substr, status)
     ! op_pj_20190410-

  END SUBROUTINE ubcnox_init_memory


  SUBROUTINE ubcnox_init_coupling

    USE messy_main_mpi_bi,           ONLY: finish, message, p_parallel_io

    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: get_tracer, get_tracer_list
    USE messy_main_channel,       ONLY: get_channel_object, get_channel_info &
                                      , new_channel_object  &
                                      , new_attribute
   USE messy_main_data_bi, ONLY: modstr_base=>modstr
   USE messy_ubcnox, ONLY: nox_switch

    IMPLICIT NONE

    INTEGER :: status, i_err, jt

    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: subnames => NULL()
    CHARACTER(LEN=*), PARAMETER :: substr = 'ubcnox_init_coupling'

    CALL start_message_bi(modstr,'INIT COUPLING',substr)


    SELECT CASE (nox_switch)

      CASE (1,2,4)      ! precalculated values with Funke Parameterization
        CALL get_channel_object(status, TRIM(nox_channel), TRIM(nox_object), p3=ptnox_ptr)
        IF (status /= 0) &
          CALL finish(substr,'channel object not found ptnox')

      CASE (3)          ! Holger Nieder
        CALL get_channel_object(status &
         , TRIM(ubcnox_kp%cha), TRIM(ubcnox_kp%obj), p1=kp_data)
        CALL channel_halt(substr, status)

        CALL get_channel_object(status &
         , TRIM(ubcnox_f107%cha), TRIM(ubcnox_f107%obj), p1=f107_data)
        CALL channel_halt(substr, status)

      CASE (5)          ! internally calculated CMIP6
        CALL get_channel_object(status &
         , TRIM(ubcnox_ap%cha), TRIM(ubcnox_ap%obj), p1=ap_data)
        CALL channel_halt(substr, status)

      CASE DEFAULT
        write(*,*) 'NOx-switch not supported'

    END SELECT

    CALL get_channel_object(status, 'tracer_gp', TRIM(target_species), p3=nox_ptr)
    IF (status /= 0) &
         CALL finish(substr,'channel object not found nox')

    CALL get_channel_object(status, 'ECHAM5', 'tm1', p3=tm1_ptr)
    IF (status /= 0) &
         CALL finish(substr,'channel object not found tm1')

    CALL get_channel_object(status, 'grid_def', 'sinlat', p2=lat_ptr)
    IF (status /= 0) &
         CALL finish(substr,'channel object not found lat')

      CALL get_tracer(i_err, GPTRSTR, TRIM(target_species), idx=idt_NOX)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, TRIM('NO2'), idx=idt_NO2)
      CALL tracer_halt(substr, i_err)


    CALL get_tracer_list(status, GPTRSTR, TRIM(target_species), idt_list_n, subnames)
    CALL tracer_halt(substr, status)
    nlntrac_n = SIZE(idt_list_n)
    IF (p_parallel_io) THEN
       DO jt=1, nlntrac_n
          IF (TRIM(subnames(jt)) == '') THEN
             WRITE(*,*) ' ... N'
          ELSE
             WRITE(*,*) ' ... N_'//TRIM(subnames(jt))
          END IF
       END DO
    END IF


    CALL end_message_bi(modstr,'INIT COUPLING',substr)

  END SUBROUTINE ubcnox_init_coupling

! !=============================================================================

SUBROUTINE ubcnox_init_tracer

!!   USE messy_main_tracer_bi, ONLY: main_tracer_initialize

   IMPLICIT NONE

!!   IF (tracer_init_required) CALL tracer_init(modstr)

END SUBROUTINE ubcnox_init_tracer

! !=============================================================================

 SUBROUTINE ubcnox_global_start
   ! ----------------------------------------------------------------------
   USE messy_main_transform_bi,  ONLY: trp_gpdc_gpgl
   USE messy_main_grid_def_mem_bi, ONLY: nlev, nlon, ngl, ngpblks, nproma
   USE messy_main_grid_def_bi,     ONLY: gboxarea => gboxarea_2d
   USE messy_main_data_bi,       ONLY:   pmid => press_3d
   USE messy_main_timer,         ONLY: lresume, current_date, lstart
   USE messy_main_timer_bi,      ONLY: event_state
   USE messy_ubcnox, ONLY: tempthresh, top_levels, ubcnox_linear_interp, ese_switch, nox_switch
   USE messy_main_constants_mem, ONLY: N_A, pi, radius_earth
   ! ----------------------------------------------------------------------
   IMPLICIT NONE

   REAL(dp), DIMENSION(:,:,:), POINTER :: temp_glob => NULL()
   REAL(dp), DIMENSION(:,:), POINTER :: lat_glob => NULL()
   REAL(dp), DIMENSION(:,:), POINTER :: area_glob => NULL()
   INTEGER :: i,j, nlat
   INTEGER :: indp,inde,inde2          ! Index for 70N and 30N (closest point)
   INTEGER :: indlev
   REAL(dp), DIMENSION(nlev,ngl) :: zm_temp
   REAL(dp), DIMENSION(ngl) :: zm_temp_1hpa, latitude, lat, xld, xlde

   REAL(dp) :: zm_np, zm_eq           ! zonal mean temperatures in northern high latitudes and in the equator region
   REAL(dp) :: tm_diff
   LOGICAL, SAVE :: startofday = .TRUE.
   LOGICAL, SAVE :: zlstart = .TRUE.
!!$REAL(dp), SAVE :: oldday = 0._dp ! op_pj_20190410
   INTEGER :: ind1, ind2, k, indlat1, indlat2
   REAL(dp) :: pref1, pref2, lat1, lat2
   REAL(dp), DIMENSION(12) :: pref
   REAL(dp), DIMENSION(9) :: lrefs, lrefn, hs, hn, he
   REAL(dp), DIMENSION(12,9) :: lds, ldn, lde
   REAL(dp), DIMENSION(ngl) :: area           !!! difference in calculation to original original code


   ! ----------------------------------------------------------------------

   IF (nox_switch==5) THEN

      pref=(/1.00, 0.70, 0.50, 0.30, 0.20, 0.15, 0.10, 0.07, 0.05, 0.03, 0.02, 0.01/)
      lrefs=(/-85.,-75.,-65.,-55.,-45.,-35.,-25.,-15., -5./)
      lrefn=(/  5., 15., 25., 35., 45., 55., 65., 75., 85./)
      lds(:, 1)=(/ 0.290, 0.285, 0.288, 0.299, 0.311, 0.319, 0.327, 0.340, 0.357, 0.374, 0.393, 0.412/)
        lds(:, 2)=(/ 0.259, 0.248, 0.241, 0.240, 0.246, 0.260, 0.281, 0.302, 0.309, 0.306, 0.297, 0.294/)
        lds(:, 3)=(/ 0.210, 0.201, 0.195, 0.193, 0.196, 0.202, 0.207, 0.208, 0.206, 0.202, 0.197, 0.187/)
        lds(:, 4)=(/ 0.153, 0.160, 0.161, 0.155, 0.147, 0.136, 0.121, 0.104, 0.092, 0.086, 0.083, 0.078/)
        lds(:, 5)=(/ 0.070, 0.084, 0.092, 0.091, 0.082, 0.069, 0.053, 0.038, 0.028, 0.024, 0.022, 0.021/)
        lds(:, 6)=(/ 0.016, 0.020, 0.021, 0.020, 0.017, 0.013, 0.010, 0.008, 0.007, 0.006, 0.006, 0.006/)
        lds(:, 7)=(/ 0.002, 0.002, 0.002, 0.001, 0.001, 0.000, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001/)
        lds(:, 8)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        lds(:, 9)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        WHERE (lds<1.e-20) lds = 1.e-20

        ldn(:, 1)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        ldn(:, 2)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        ldn(:, 3)=(/ 0.003, 0.004, 0.004, 0.004, 0.003, 0.003, 0.003, 0.003, 0.002, 0.002, 0.002, 0.002/)
        ldn(:, 4)=(/ 0.011, 0.016, 0.021, 0.023, 0.023, 0.020, 0.016, 0.011, 0.008, 0.006, 0.005, 0.005/)
        ldn(:, 5)=(/ 0.036, 0.049, 0.060, 0.066, 0.067, 0.057, 0.044, 0.030, 0.023, 0.019, 0.018, 0.018/)
        ldn(:, 6)=(/ 0.096, 0.106, 0.120, 0.130, 0.139, 0.131, 0.115, 0.099, 0.088, 0.083, 0.081, 0.080/)
        ldn(:, 7)=(/ 0.185, 0.184, 0.191, 0.206, 0.229, 0.245, 0.250, 0.246, 0.238, 0.230, 0.228, 0.226/)
        ldn(:, 8)=(/ 0.307, 0.291, 0.273, 0.263, 0.261, 0.278, 0.299, 0.318, 0.325, 0.326, 0.325, 0.325/)
        ldn(:, 9)=(/ 0.362, 0.351, 0.332, 0.308, 0.278, 0.266, 0.273, 0.293, 0.316, 0.333, 0.342, 0.344/)
        WHERE (ldn<1.e-20) ldn = 1.e-20

        lde(:, 1)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        lde(:, 2)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        lde(:, 3)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        lde(:, 4)=(/ 0.001, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.001, 0.001, 0.001, 0.001/)
        lde(:, 5)=(/ 0.009, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.007, 0.007, 0.006, 0.006/)
        lde(:, 6)=(/ 0.044, 0.042, 0.039, 0.036, 0.034, 0.035, 0.037, 0.040, 0.042, 0.042, 0.042, 0.041/)
        lde(:, 7)=(/ 0.146, 0.147, 0.148, 0.139, 0.135, 0.132, 0.137, 0.145, 0.150, 0.153, 0.153, 0.150/)
        lde(:, 8)=(/ 0.333, 0.331, 0.330, 0.334, 0.339, 0.341, 0.344, 0.343, 0.341, 0.336, 0.332, 0.330/)
        lde(:, 9)=(/ 0.467, 0.471, 0.472, 0.481, 0.481, 0.482, 0.471, 0.462, 0.459, 0.462, 0.467, 0.473/)
        WHERE (lde<1.e-20) lde = 1.e-20

! op_pj_20190410+: Note: this is also .TRUE. at each restart (lresume),
!                 because oldday(SAVE) is initialised above with zero.
! This said, the resulting temp_glob below depends on the restart
! time. It is either the temperature field at the
! first time step at a new day (regular), or,
! after an irregular (QTIMER) restart, the actual
! temperature field ...
!
!!$   IF (abs(DAY-oldday)>0.5_dp) THEN
!!$     startofday = .TRUE.
!!$   END IF
!
      ! CHECK EVENT: daily
      startofday   = event_state(XTIMER_DAILY,   current_date)

      IF (startofday .OR. lstart) THEN
         !  SAVE tm1 at beginning of day for restart
         tm1_save_ptr(:,:,:) = tm1_ptr(:,:,:)
      END IF

! op_pj_20190410-

! op_pj_20190410+
! This needs to be done also immediately after restart ...
!!$   IF (startofday) THEN
      IF (startofday .OR. lresume .OR. lstart) THEN
! op_pj_20190410-
      ! decompose fields
! op_pj_20190410+
! Use temperature from start of the day to be consistent after restart ...
!!$    CALL trp_gpdc_gpgl(1, tm1_ptr,   temp_glob)
       CALL trp_gpdc_gpgl(1, tm1_save_ptr,   temp_glob)
! luckily sin(latitude) and the grid-box area are time independent ...
! op_pj_20190410-
       CALL trp_gpdc_gpgl(1, lat_ptr, lat_glob)
       CALL trp_gpdc_gpgl(1, gboxarea, area_glob)
       latitude(:)=asin(lat_glob(1,:))/pi*180._dp
       DO i=1,ngl
         area(i)=sum(area_glob(:,i),1)*1.e4          ! conversion from m?? to cm??
       END DO
   ! calc zonal mean temp
       DO i=1,nlev
         DO j=1,ngl
           zm_temp(i,j)=sum(temp_glob(:,i,j),1)/(max(1,size(temp_glob(:,i,j),1)))
         END DO
       END DO
       indlev=minloc(abs(pmid(1,:,1)-100._dp),1)
       indp=minloc(abs(latitude-70),1)
       inde=minloc(abs(latitude-30),1)
       inde2=ngl/2

       zm_temp_1hpa(:)=zm_temp(indlev,:)
       zm_np=sum(zm_temp_1hpa(1:indp),1)/(max(1,size(zm_temp_1hpa(1:indp),1)))
       zm_eq=sum(zm_temp_1hpa(inde:inde2),1)/(max(1,size(zm_temp_1hpa(inde:inde2),1)))
       tm_diff=zm_eq-zm_np

! ka_sv_201705023+
! bugfix for counting of days of ese since beginning of ese

       ! op_pj_20190410+
       ! op_pj_20190416+
       ! fixes potential memory leak:
       IF (ASSOCIATED(temp_glob)) THEN
          DEALLOCATE(temp_glob) ; NULLIFY(temp_glob)
       END IF
       IF (ASSOCIATED(lat_glob)) THEN
          DEALLOCATE(lat_glob) ; NULLIFY(lat_glob)
       END IF
       IF (ASSOCIATED(area_glob)) THEN
          DEALLOCATE(area_glob) ; NULLIFY(area_glob)
       END IF
       ! op_pj_20190416-
    END IF ! startofday .OR. lresume .OR. lstart

    ! The following should happen independently of the restart, only
    ! at the beginning of a new day ...
    IF (startofday) THEN
    ! op_pj_20190410-

       IF (daysofesecalc>0.9_dp) THEN  ! increment, if already started ...
          ! What about very short ESE?
          daysofesecalc=daysofesecalc + 1._dp
       END IF
!!! CRASH!!!! IF nox_switch==3 ...
       ! op_pj_20190410: check condition
       IF (tm_diff>tempthresh) THEN
! op_pj_20190410+: was wrong and is now obsolete ...
!!$     IF (HOUR==0) THEN    ! don't count a new day if start of model is not at the beginning of one day, e.g. after restart
! op_pj_20190410-
           daysofese=daysofese + 1.0_dp
           IF (daysofesecalc < 0.5_dp) THEN ! first time
              daysofesecalc = 1._dp
           END IF
!!$     END IF ! op_pj_20190410
        ELSE
           daysofese = 0.0_dp
        END IF

        IF (daysofesecalc > 200._dp) THEN ! stop after 200 days
           daysofesecalc = 0._dp
        END IF
! ka_sv_201705023-

!!$       startofday=.FALSE. ! op_pj_20190410
     END IF ! startofday

!!$     oldday=DAY ! op_pj_20190410

   ! Things only needed once (at start and after restart)

     IF (zlstart) THEN
       ALLOCATE(fac(nlon,nlev,ngl))
       ALLOCATE(face(nlon,nlev,ngl))
       CALL trp_gpdc_gpgl(1, lat_ptr, lat_glob)
       lat(:)=asin(lat_glob(1,:))/3.141592654_dp*180._dp
       ! op_pj_20190416+
       IF (ASSOCIATED(lat_glob)) THEN
          DEALLOCATE(lat_glob); NULLIFY(lat_glob)
       ENDIF
       ! op_pj_20190416-
       DO i=1,top_levels
! QQQ+
! op_pj_20190410: Well, this bears the high risk that the results
!                 become dependent on the
!                 - restart frequency (because pmid changes with time
!                   and therefore the field differs after an irregular
!                   restart (QTIMER),
!                 - parallel decompositiom, because element
!                   (1,i,1) is used, which is the first one in the regional
!                   domain decomposition.
!                 It works presumably, as long as
!                   hybm(1:to_levels) == 0.0_dp
!                 because in this case pmid(:,i,:) = const.
!
          ind1=minloc(abs(pref-pmid(1,i,1)),DIM=1)
          pref1=pref(ind1)
          ind2=MAX(ind1-1,1)
          IF (pref1 > pmid(1,i,1)) THEN
            ind2=MIN(ind1+1,12)
          END IF
          pref2=pref(ind2)

          DO j=1,9
            CALL ubcnox_linear_interp(hs(j),pmid(1,i,1),pref1,pref2,lds(ind1,j),lds(ind2,j))
            CALL ubcnox_linear_interp(hn(j),pmid(1,i,1),pref1,pref2,ldn(ind1,j),ldn(ind2,j))
            IF (ese_switch==1) THEN          ! only if ESEs are used
              CALL ubcnox_linear_interp(he(j),pmid(1,i,1),pref1,pref2,lde(ind1,j),lde(ind2,j))
            END IF
          END DO

         DO k=1,ngl
          IF (lat(k)<0._dp) THEN
            indlat1=minloc(abs(lrefs-lat(k)),DIM=1)
             IF (lat(k)<lrefs(indlat1)) THEN
               indlat2=MAX(indlat1-1,1)
             END IF
             IF (lat(k)>=lrefs(indlat1)) THEN
               indlat2=MIN(indlat1+1,9)
             END IF
             lat1=lrefs(indlat1)
             lat2=lrefs(indlat2)

            CALL ubcnox_linear_interp(xld(k),lat(k),lat1,lat2,hs(indlat1),hs(indlat2))
            IF (lat(k)<lrefs(1)) THEN
              xld(k)=hs(1)
            END IF
            IF (lat(k)>lrefs(9)) THEN
              xld(k)=hs(9)
            END IF

          END IF
          IF (lat(k)>=0._dp) THEN

            indlat1=minloc(abs(lrefn-lat(k)),DIM=1)
            IF (lat(k)<lrefn(indlat1)) THEN
               indlat2=MAX(indlat1-1,1)
             END IF
             IF (lat(k)>=lrefn(indlat1)) THEN
               indlat2=MIN(indlat1+1,9)
             END IF
             lat1=lrefn(indlat1)
             lat2=lrefn(indlat2)

            CALL ubcnox_linear_interp(xld(k),lat(k),lat1,lat2,hn(indlat1),hn(indlat2))
            IF (lat(k)>lrefn(9)) THEN
              xld(k)=hn(9)
            END IF
            IF (lat(k)<lrefn(1)) THEN
              xld(k)=hn(1)
            END IF

          END IF


          zlstart=.FALSE.
         END DO ! ngl

         DO k=1,ngl
           IF (k<=ngl/2) THEN
             fac(:,i,k)=xld(k)/(sum(xld(1:ngl/2)*area(1:ngl/2)/N_A))*1.e4
           END IF
           IF (k>ngl/2) THEN
             fac(:,i,k)=xld(k)/(sum(xld(ngl/2+1:ngl)*area(ngl/2+1:ngl)/N_A))*1.e4
           END IF
         END DO ! ngl

         ! for ESE only
         IF (ese_switch==1) THEN
           DO k=1,ngl
             IF (lat(k)>=0._dp) THEN
               indlat1=minloc(abs(lrefn-lat(k)),DIM=1)
               IF (lat(k)<lrefn(indlat1)) THEN
                 indlat2=MAX(indlat1-1,1)
               END IF
               IF (lat(k)>=lrefn(indlat1)) THEN
                 indlat2=MIN(indlat1+1,9)
               END IF
               lat1=lrefn(indlat1)
               lat2=lrefn(indlat2)

               CALL ubcnox_linear_interp(xlde(k),lat(k),lat1,lat2,he(indlat1),he(indlat2))
             ELSE
               xlde(k)=0._dp
             END IF
           END DO
           DO k=1,ngl
             face(:,i,k)=xlde(k)/(sum(xlde(:)*area(:)/N_A))*1.e4
           END DO
         END IF   ! ESE

        END DO ! top_levels_in

     END IF  ! only start
! QQQ- ! op_pj_20190410 (see above, block ends here)
   END IF  ! nox_switch==5

 END SUBROUTINE ubcnox_global_start



! !=============================================================================
 SUBROUTINE ubcnox_physc
  !--------------------------------------------------------------------------
  USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, pxtm1 => qxtm1
!!$  USE messy_main_tracer_mem_bi, ONLY: temp_glob => pxtm1
  USE messy_main_grid_def_mem_bi, ONLY: nlev, nlon, jrow, kproma, nproma, ngpblks
  USE messy_main_grid_def_bi, ONLY: philat_2d, gboxarea_2d
  USE messy_main_data_bi,     ONLY:pmid => press_3d, temp => t_scb & ! mid-level pressures [Pa]
                                    , pint => pressi_3d
  USE messy_ubcnox, ONLY: nox_switch, top_levels, ese_switch,  &
           ubcnox_calc_current_ubc

  USE messy_main_constants_mem, ONLY: M_air, R_gas, N_A, g
  USE messy_main_timer, ONLY: time_span_d, MONTH, DAY, YEAR, time_step_len
  USE mo_geoloc, ONLY: philat_2d
  USE messy_main_transform_bi,  ONLY: trp_gpdc_gpgl


  IMPLICIT none

  INTEGER :: jp, jk, n,  i, status
  INTEGER :: idt, jt

  REAL(DP)                    :: zpb, zpt
  REAL(dp), DIMENSION(nproma,nlev) :: grheight ! height of box [m]

  REAL(dp), DIMENSION(nproma,nlev) :: PTNOX1_t0
  REAL(dp), DIMENSION(nproma,nlev) :: PTNOX1_g
  REAL(dp), DIMENSION(nproma,nlev) :: NOX_t0
  REAL(dp), DIMENSION(nproma,nlev) :: NOX_g


  REAL(dp), DIMENSION(nproma,nlev) :: cair
  REAL(dp), DIMENSION(nproma)      :: local_nox
  INTEGER :: ilat, ikp   ! Indices for LUT
  INTEGER :: if10a, if10b
  REAL(dp) :: f10fak

   REAL(dp), DIMENSION(nlev) :: diff  ! Difference in Pa between pressuregrid and 100 Pa
    INTEGER :: ind_lev_1hpa  ! nearest level to 100 Pa

   REAL(dp), DIMENSION(top_levels) :: ep_temp, ubc_temp, bg_temp, ese_temp  ! EPP-noy field (without background)
   REAL(dp) :: dn , ds


! get level-number for 1 hPa
! assume that level number is everywhere the same
    diff=abs(pmid(1,:,1)-100)
    ind_lev_1hpa = minloc(diff,DIM=1)

  IF (nox_switch==1 .or. nox_switch==2 .or. nox_switch==4) then

  level_loop0: DO jk=1,nlev
    vector_loop0: DO jp=1,kproma
        PTNOX1_t0(jp,jk)  = ptnox_ptr(jp,jk,jrow)
        idt = idt_NOX
        NOX_t0(jp,jk) = pxtm1(jp,jk,idt_NOX) &
                   +pxtte(jp,jk,idt_NOX)*time_step_len

    END DO vector_loop0
  END DO level_loop0

! calculate gas densitys
level_loop1d: DO jk=1,nlev
    vector_loop1d: DO jp=1,kproma
        cair(jp,jk)  = (N_A) * pmid(jp,jk,jrow) / (R_gas*temp(jp,jk,jrow)) / 1.e06  ! cm^(-3)
      END DO vector_loop1d
    END DO level_loop1d

  end if

  IF (nox_switch==1) then

! start ESE detection
  IF (ese_switch==1) THEN   ! Bernd Funke: (T(0-30deg)-T(70-90deg) @ 1 hPa > 45K)
    ! level-number for 1 hPa: ind_lev_1hpa
!    FORALL (jp=1:kproma,philat_2d(jp,jrow)>0._dp .and. philat_2d(jp,jrow)<30._dp)
!      sum_eq=sum_eq+temp(jp,ind_lev_1hpa,jrow)
!      num_eq=num_eq+1
!    END FORALL

  END IF ! ese_switch == 1

! end ESE detection
  NOX_g=NOX_t0
  tend=0._dp
  NOX_g(1:kproma,1) = PTNOX1_t0(1:kproma,1)/cair(1:kproma,1)

           tend(1:kproma,:,jrow) =              &
             + (NOX_g(1:kproma,:)-NOX_t0(1:kproma,:))                       &
             /time_step_len

 IF (nlntrac_n > 0) THEN
       DO jt=1, nlntrac_n
          pxtm1(1:kproma,1,idt_list_n(jt)) = NOX_g(1:kproma,1)
          pxtte(1:kproma,1,idt_list_n(jt)) = 0._dp
          pxtm1(1:kproma,1,idt_NO2) = 0._dp
          pxtte(1:kproma,1,idt_NO2) = 0._dp
       END DO
    END IF


END IF ! nox_switch

IF (nox_switch==4) then
! upper boundary amount; top 4 level
  NOX_g=NOX_t0
  tend=0._dp

  DO jk=1,4
     NOX_g(1:kproma,jk) = PTNOX1_t0(1:kproma,jk)/cair(1:kproma,jk)
  END DO
  tend(1:kproma,:,jrow) =              &
       + (NOX_g(1:kproma,:)-NOX_t0(1:kproma,:))                       &
       /time_step_len

!  write(*,*) " DURING SECOND LOOP"

 IF (nlntrac_n > 0) THEN
       DO jt=1, nlntrac_n
          DO jk=1,4
            pxtm1(1:kproma,jk,idt_list_n(jt)) = NOX_g(1:kproma,jk)
            pxtte(1:kproma,jk,idt_list_n(jt)) = 0._dp
            pxtm1(1:kproma,jk,idt_NO2) = 0._dp
            pxtte(1:kproma,jk,idt_NO2) = 0._dp
          END DO
       END DO
    END IF
END IF

IF (nox_switch==2) then
! upper boundary flux in #/cm??/s
   NOX_g=0._dp
   tend=0._dp
   NOX_g(1:kproma,1) = PTNOX1_t0(1:kproma,1)*1.e04   &  ! convert cm?? to m??
                         * gboxarea_2d(1:kproma,jrow)*time_step_len    &
                         /(cair(1:kproma,1)*1.e06*grvol(1:kproma,1,jrow))


   DO jp=1,kproma
             tend(jp,1,jrow) = NOX_g(jp,1)/time_step_len
   END DO


    !debug(1:kproma,:,jrow) = NOX_g(1:kproma,:)
    IF (nlntrac_n > 0) THEN
       ! ADD N TO TENDENCY
       DO jt=1, nlntrac_n
          pxtte(1:kproma,:,idt_list_n(jt)) = pxtte(1:kproma,:,idt_list_n(jt)) + &
               tend(1:kproma,:,jrow)
       END DO
    END IF

END IF

IF (nox_switch==3) then
   tend=0._dp


  ! get correct kp bin
  IF (kp_data(1)<0._dp) then
     ikp=1
  END IF
  IF (kp_data(1)>ubc_kp(14)) then
     ikp=14
  END IF
  DO i=2,nkp
    IF (kp_data(1)<=ubc_kp(i) .AND. kp_data(1)>ubc_kp(i-1)) THEN
      ikp=i-1
    END IF
  END DO
  DO i=2,nf107
    IF (f107_data(1)<=ubc_f107(i) .AND. f107_data(1)>ubc_f107(i-1)) THEN
      if10a=i-1
      if10b=i
      f10fak=(f107_data(1)-ubc_f107(if10a))/(ubc_f107(if10b)-ubc_f107(if10a))
    END IF
  END DO
  IF (f107_data(1)<=ubc_f107(1)) THEN
     if10a=1
     if10b=1
     f10fak=0.5
  END IF
  IF (f107_data(1)>=ubc_f107(3)) THEN
     if10a=3
     if10b=3
     f10fak=0.5
  END IF

  DO jp=1,kproma
  ! get correct latitude bin
     DO i=1,nlat
        IF (philat_2d(jp,jrow)<=ubc_lat_bin(i) .AND. philat_2d(jp,jrow)>ubc_lat_bin(i+1)) THEN
          ilat=i
        END IF
     END DO
     local_nox(jp)=f10fak*noxlut(ikp, MONTH, if10a, ilat)+(1._dp-f10fak)*noxlut(ikp, MONTH, if10b, ilat)*0.2

  END DO
  tend(1:kproma,1,jrow)=local_nox(1:kproma)
  IF (nlntrac_n > 0) THEN
       ! ADD SPE N TO TENDENCY
       DO jt=1, nlntrac_n
          pxtte(1:kproma,:,idt_list_n(jt)) = pxtte(1:kproma,:,idt_list_n(jt)) + &
               tend(1:kproma,:,jrow)/86400._dp
       END DO
    END IF

END IF

 IF (nox_switch==5) THEN
   ALLOCATE(facdec(nproma,nlev,ngpblks))
   ALLOCATE(facedec(nproma,nlev,ngpblks))
   CALL trp_gpdc_gpgl(-1, facdec, fac)
   CALL trp_gpdc_gpgl(-1, facedec, face)

   IF (MONTH<7) THEN
      CALL time_span_d(dn,YEAR,1,1,0,0,0,YEAR,MONTH,DAY,0,0,0)
      dn=dn+185
   ELSE
      CALL time_span_d(dn,YEAR,7,1,0,0,0,YEAR,MONTH,DAY,0,0,0)
      dn=dn+1
   END IF
   CALL time_span_d(ds,YEAR,1,1,0,0,0,YEAR,MONTH,DAY,0,0,0)
   ds=ds+1

   DO jp=1,kproma

     CALL ubcnox_calc_current_ubc(ubc_temp,ep_temp,bg_temp,ese_temp,&
        top_levels,ap_data,pmid(jp,1:top_levels,jrow)/100._dp,philat_2d(jp,jrow),gboxarea_2d(jp,jrow)*10000._dp,nlon,facdec(jp,1:top_levels,jrow), &              ! pressure conversion from Pa to hPa
! ka_sv_201705023+
        ese_switch, nint(daysofesecalc), facedec(jp,1:top_levels,jrow), dn, ds)
! ka_sv_201705023-
       cair(jp,1:top_levels)  = (N_A) * pmid(jp,1:top_levels,jrow) / (R_gas*temp(jp,1:top_levels,jrow)) / 1.e06  ! cm^(-3)
       ubc(jp,1:top_levels,jrow) = ubc_temp(1:top_levels)/cair(jp,1:top_levels)
       ep(jp,1:top_levels,jrow) = ep_temp(1:top_levels)/cair(jp,1:top_levels)
       bg(jp,1:top_levels,jrow) = bg_temp(1:top_levels)/cair(jp,1:top_levels)
       ese(jp,1:top_levels,jrow) = ese_temp(1:top_levels)/cair(jp,1:top_levels)

   END DO
   IF (nlntrac_n > 0) THEN
       DO jt=1, nlntrac_n
          DO jk=1,top_levels
            pxtm1(1:kproma,jk,idt_list_n(jt)) = ubc(1:kproma,jk,jrow)
            pxtte(1:kproma,jk,idt_list_n(jt)) = 0._dp
            pxtm1(1:kproma,jk,idt_NO2) = 0._dp
            pxtte(1:kproma,jk,idt_NO2) = 0._dp
          END DO
       END DO
    END IF
    DEALLOCATE(facdec)
    DEALLOCATE(facedec)
 END IF



END SUBROUTINE ubcnox_physc

!=============================================================================

!=============================================================================
SUBROUTINE ubcnox_free_memory

  IMPLICIT NONE

  IF (ASSOCIATED(fac)) THEN
     DEALLOCATE(fac); NULLIFY(fac)
  END IF

  IF (ASSOCIATED(face)) THEN
     DEALLOCATE(face); NULLIFY(face)
  END IF

END SUBROUTINE ubcnox_free_memory
!=============================================================================

!=============================================================================

SUBROUTINE ubcnox_read_nml_cpl(status, iou)
  !------------------------------------------
  ! read coupling namelist /CPL/ from maoam.nml
  !------------------------------------------

  ! MESSy
  USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

  IMPLICIT NONE

  INTEGER, INTENT(out) :: status
  INTEGER, INTENT(in) :: iou

  CHARACTER(len=*), PARAMETER :: substr='ubcnox_read_nml_cpl'
  LOGICAL :: lex     ! file exists?
  INTEGER :: fstat   ! file status

  NAMELIST /CPL/ nox_channel, nox_object, target_species, ubcnox_kp, ubcnox_f107, &
                 ubcnox_ap

  status = 1             ! initialise status flag with error code

  CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
  IF (.not.lex) RETURN   ! error: psc.nml does not exist
  read(iou, nml=cpl, iostat=fstat)
  CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
  IF (fstat/=0) RETURN   ! error while reading namelist


  CALL read_nml_close(substr, iou, modstr)
  status = 0   ! no error

END SUBROUTINE ubcnox_read_nml_cpl



!=============================================================================

END MODULE messy_ubcnox_e5
