#include "messy_main_ppd_bi.inc"

! **********************************************************************
MODULE messy_orbit_si
! **********************************************************************

  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_channel,       ONLY: t_chaobj_cpl
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  USE messy_main_timer_bi,      ONLY: nmonth
  USE messy_orbit,              ONLY: modstr

  IMPLICIT NONE
  PRIVATE
  SAVE
  INTRINSIC :: NULL

  ! CPL namelist
  ! switch off (l_quiet=.TRUE.) messages to log-file during time loop
  LOGICAL :: l_quiet = .FALSE.

  ! CHANNEL OBJECTS FOR ORBIT PARAMETERS AT CURRENT TIME STEP
  REAL(DP), POINTER :: cdisse => NULL()  ! distance sun - earth [AU]
  REAL(DP), POINTER :: dec    => NULL()  ! declination of sun
  REAL(DP), POINTER :: ra     => NULL()  ! right ascension of sun
  ! cos(solar zenith angle)
  REAL(DP), DIMENSION(:,:), POINTER :: cossza => NULL()
  ! relative day length
  REAL(DP), DIMENSION(:,:), POINTER :: rdayl => NULL()
  ! cos(solar zenith angle); cut-off at zero
  REAL(DP), DIMENSION(:,:), POINTER :: cosszac => NULL()

  ! FOR ORBIT PARAMETERS AT RADIATION TIME STEP
  TYPE(t_chaobj_cpl) :: c_rad_offset ! op_pj_20110715
  LOGICAL            :: l_calc4rad = .FALSE. ! special calculation for radiation
  REAL(DP), POINTER  :: dt_offset_rad => NULL()  ! offset (in seconds)
  !
  ! CHANNEL OBJECTS
  REAL(DP), POINTER :: cdissem => NULL()  ! distance sun - earth [AU]
  REAL(DP), POINTER :: decm    => NULL()  ! declination of sun
  REAL(DP), POINTER :: ram     => NULL()  ! right ascension of sun
  ! cos(solar zenith angle)
  REAL(DP), DIMENSION(:,:), POINTER :: cosszam => NULL()
  ! relative day length
  REAL(DP), DIMENSION(:,:), POINTER :: rdaylm => NULL()
  ! cos(solar zenith angle); cut-off at zero
  REAL(DP), DIMENSION(:,:), POINTER :: cosszacm => NULL() 

  ! TIME SHIFTED ORBIT PARAMETERS
  TYPE(t_chaobj_cpl) :: c_offset
  REAL(DP), POINTER  :: dt_offset => NULL()  ! offset (in seconds)
  LOGICAL            :: l_alloc_dt_offset = .FALSE.
  !
  ! CHANNEL OBJECTS
  REAL(DP), POINTER :: cdisse_off => NULL()  ! distance sun - earth [AU]
  REAL(DP), POINTER :: dec_off    => NULL()  ! declination of sun
  REAL(DP), POINTER :: ra_off     => NULL()  ! right ascension of sun
  ! cos(solar zenith angle)
  REAL(DP), DIMENSION(:,:), POINTER :: cossza_off => NULL()
  ! relative day length
  REAL(DP), DIMENSION(:,:), POINTER :: rdayl_off => NULL()
  ! cos(solar zenith angle); cut-off at zero
  REAL(DP), DIMENSION(:,:), POINTER :: cosszac_off => NULL()


  PUBLIC :: orbit_initialize
  PUBLIC :: orbit_init_memory
  PUBLIC :: orbit_init_coupling
  PUBLIC :: orbit_global_start
  PUBLIC :: orbit_free_memory

  !PRIVATE :: calc_orbit_param
  !PRIVATE :: orbit_read_nml_cpl

CONTAINS

  ! ---------------------------------------------------------------------------
  SUBROUTINE orbit_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit
    ! MESSy
    USE messy_orbit,         ONLY: orbit_read_nml_ctrl &
                                 , cecc, cobld, clonp, l_orbvsop87

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER     :: substr = 'orbit_initialize'
    INTEGER                         :: iou    ! I/O unit
    INTEGER                         :: status ! error status

  !--- Read CTRL namelist
  IF (p_parallel_io) THEN
    iou = find_next_free_unit(100,200)
    CALL orbit_read_nml_ctrl(status, iou)
    IF (status /= 0) CALL error_bi(' ',substr)
  END IF

  !--- Broadcast over PEs:
  CALL p_bcast (cecc,  p_io)
  CALL p_bcast (cobld, p_io)
  CALL p_bcast (clonp, p_io)
  CALL p_bcast (l_orbvsop87, p_io)

  !--- Read CPL namelist
  IF (p_parallel_io) THEN
    iou = find_next_free_unit(100,200)
    CALL orbit_read_nml_cpl(status, iou)
    IF (status /= 0) CALL error_bi(' ',substr)
  END IF

  !--- Broadcast over PEs:
  CALL p_bcast(l_quiet, p_io)
  CALL p_bcast(c_rad_offset%cha, p_io)
  CALL p_bcast(c_rad_offset%obj, p_io)
  CALL p_bcast(c_offset%cha, p_io)
  CALL p_bcast(c_offset%obj, p_io)

  END SUBROUTINE orbit_initialize
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE orbit_init_memory

    ! MESSy/BMIL
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: SCALAR, GP_2D_HORIZONTAL

    ! MESSy
    USE messy_main_channel,    ONLY: new_channel, new_channel_object &
                                   , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'orbit_init_memory'
    INTEGER :: status

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ! DEFINE NEW CHANNEL
    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)

    ! OBJECTS
    IF (p_parallel_io) WRITE(*,*) ' ... cdisse'
    CALL new_channel_object(status, modstr, 'cdisse' &
         , p0 = cdisse, reprid = SCALAR )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdisse'   &
         , 'long_name', c='distance sun - earth' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdisse'   &
         , 'units', c='AU' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) WRITE(*,*) ' ... dec'
    CALL new_channel_object(status, modstr, 'dec' &
         , p0 = dec, reprid = SCALAR )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dec'   &
         , 'long_name', c='declination of sun' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dec'   &
         , 'units', c='rad' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) WRITE(*,*) ' ... ra'
    CALL new_channel_object(status, modstr, 'ra' &
         , p0 = ra, reprid = SCALAR )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ra'   &
         , 'long_name', c='right ascension of sun' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ra'   &
         , 'units', c='rad' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) WRITE(*,*) ' ... cossza'
    CALL new_channel_object(status, modstr, 'cossza' &
         , p2 = cossza, reprid = GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cossza'   &
         , 'long_name', c='cos(solar zenith angle)' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) WRITE(*,*) ' ... rdayl'
    CALL new_channel_object(status, modstr, 'rdayl' &
         , p2 = rdayl, reprid = GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rdayl'   &
         , 'long_name', c='relative day length' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) WRITE(*,*) ' ... cosszac'
    CALL new_channel_object(status, modstr, 'cosszac' &
         , p2 = cosszac, reprid = GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cosszac'   &
         , 'long_name', c='cos(solar zenith angle); cut-off' )
    CALL channel_halt(substr, status)

    ! ORBIT PARAMETERS WITH TIME OFFSET
    IF (p_parallel_io) WRITE(*,*) ' ... cdisse_off'
    CALL new_channel_object(status, modstr, 'cdisse_off' &
         , p0 = cdisse_off, reprid = SCALAR )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdisse_off'   &
         , 'long_name', c='distance sun - earth' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdisse_off'   &
         , 'units', c='AU' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) WRITE(*,*) ' ... dec_off'
    CALL new_channel_object(status, modstr, 'dec_off' &
         , p0 = dec_off, reprid = SCALAR )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dec_off'   &
         , 'long_name', c='declination of sun' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dec_off'   &
         , 'units', c='rad' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) WRITE(*,*) ' ... ra_off'
    CALL new_channel_object(status, modstr, 'ra_off' &
         , p0 = ra_off, reprid = SCALAR )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ra_off'   &
         , 'long_name', c='right ascension of sun' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ra_off'   &
         , 'units', c='rad' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) WRITE(*,*) ' ... cossza_off'
    CALL new_channel_object(status, modstr, 'cossza_off' &
         , p2 = cossza_off, reprid = GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cossza_off'   &
         , 'long_name', c='cos(solar zenith angle)' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) WRITE(*,*) ' ... rdayl_off'
    CALL new_channel_object(status, modstr, 'rdayl_off' &
         , p2 = rdayl_off, reprid = GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rdayl_off'   &
         , 'long_name', c='relative day length' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) WRITE(*,*) ' ... cosszac_off'
    CALL new_channel_object(status, modstr, 'cosszac_off' &
         , p2 = cosszac_off, reprid = GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cosszac_off'   &
         , 'long_name', c='cos(solar zenith angle); cut-off' )
    CALL channel_halt(substr, status)    

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE orbit_init_memory
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE orbit_init_coupling

    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: info_bi, error_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: SCALAR, GP_2D_HORIZONTAL
    USE messy_main_channel,          ONLY: get_channel_object &
                                         , new_channel_object, new_attribute
    USE messy_main_timer,            ONLY: delta_time

    IMPLICIT NONE
    INTRINSIC :: TRIM, ADJUSTL

    INTEGER :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'orbit_init_coupling'

    CALL start_message_bi(modstr,'RADIATION COUPLING',substr)

    ! CHECK RADIATION OFFSET
    IF ( (TRIM(c_rad_offset%cha) /= '') .AND. &
         (TRIM(c_rad_offset%obj) /= '') ) THEN
       CALL info_bi('Looking for radiation offset [s] ... ')
       CALL info_bi('       channel: '//c_rad_offset%cha)
       CALL info_bi('       object : '//c_rad_offset%obj)
       CALL get_channel_object(status &
            , TRIM(c_rad_offset%cha), TRIM(c_rad_offset%obj), p0=dt_offset_rad)
       l_calc4rad = (status == 0)
       IF (.NOT. l_calc4rad) THEN
          CALL error_bi('channel / object not found', substr)
       ELSE
          CALL info_bi(' ... OK.', substr)
       END IF
    ELSE
       l_calc4rad = .FALSE.
       CALL info_bi('empty c_rad_offset in CPL namelist',substr)
    END IF

    CALL end_message_bi(modstr,'RADIATION COUPLING',substr)
    
    radiation: IF (l_calc4rad) THEN

       CALL start_message_bi(modstr, 'ADDITIONAL CHANNEL OBJECTS', substr)

       ! OBJECTS
       ! NOTE: These objects need to have lrestreq=.TRUE. since they are
       !       only calculated in "radiation time steps".
       !       For "radiation time steps", dt_offset_rad >= 0 (see below).
       IF (p_parallel_io) WRITE(*,*) ' ... cdissem'
       CALL new_channel_object(status, modstr, 'cdissem' &
            , p0 = cdissem, reprid = SCALAR, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'cdissem'   &
            , 'long_name', c='distance sun - earth (radiation calc.)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'cdissem'   &
            , 'units', c='AU' )
       CALL channel_halt(substr, status)

       IF (p_parallel_io) WRITE(*,*) ' ... decm'
       CALL new_channel_object(status, modstr, 'decm' &
            , p0 = decm, reprid = SCALAR, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'decm'   &
            , 'long_name', c='declination of sun (radiation calc.)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'decm'   &
         , 'units', c='rad' )
       CALL channel_halt(substr, status)

       IF (p_parallel_io) WRITE(*,*) ' ... ram'
       CALL new_channel_object(status, modstr, 'ram' &
            , p0 = ram, reprid = SCALAR, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'ram'   &
            , 'long_name', c='right ascension of sun (radiation calc.)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'ram'   &
            , 'units', c='rad' )
       CALL channel_halt(substr, status)

       IF (p_parallel_io) WRITE(*,*) ' ... cosszam'
       CALL new_channel_object(status, modstr, 'cosszam' &
            , p2 = cosszam, reprid = GP_2D_HORIZONTAL, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'cosszam'   &
            , 'long_name', c='cos(solar zenith angle) (radiation calc.)' )
       CALL channel_halt(substr, status)

       IF (p_parallel_io) WRITE(*,*) ' ... rdaylm'
       CALL new_channel_object(status, modstr, 'rdaylm' &
            , p2 = rdaylm, reprid = GP_2D_HORIZONTAL, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'rdaylm'   &
            , 'long_name', c='relative day length (radiation calc.)' )
       CALL channel_halt(substr, status)

       IF (p_parallel_io) WRITE(*,*) ' ... cosszacm'
       CALL new_channel_object(status, modstr, 'cosszacm' &
            , p2 = cosszacm, reprid = GP_2D_HORIZONTAL, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'cosszacm'   &
            , 'long_name' &
            , c='cos(solar zenith angle); cut-off (radiation calc.)' )
       CALL channel_halt(substr, status)
       
       CALL end_message_bi(modstr, 'ADDITIONAL CHANNEL OBJCTS', substr)

    END IF radiation

    CALL start_message_bi(modstr,'DETERMINE OFFSET',substr)

    SELECT CASE (TRIM(ADJUSTL(c_offset%cha)))
    CASE ('')
       ! no offset
       ALLOCATE(dt_offset)
       l_alloc_dt_offset = .TRUE.
       dt_offset = 0.0_dp
       CALL info_bi('offset is 0.0 seconds',substr)
    CASE('#')
       ! offset in seconds
       ALLOCATE(dt_offset)
       l_alloc_dt_offset = .TRUE.
       dt_offset = 0.0_dp
       CALL info_bi('offset is '//TRIM(c_offset%obj)//' seconds',substr)
       READ(unit=c_offset%obj,fmt=*) dt_offset
    CASE('#*')
       ! offset in time-steps
       ALLOCATE(dt_offset)
       l_alloc_dt_offset = .TRUE.
       dt_offset = 0.0_dp
       CALL info_bi('offset is '//TRIM(c_offset%obj)//' time steps',substr)
       READ(unit=c_offset%obj,fmt=*) dt_offset
       dt_offset = dt_offset * delta_time
    CASE DEFAULT
       l_alloc_dt_offset = .FALSE.
       CALL info_bi('offset is provided by object'//TRIM(c_offset%obj)//&
            &'from channel '//TRIM(c_offset%cha),substr)
       CALL get_channel_object(status &
            , TRIM(c_offset%cha), TRIM(c_offset%obj), p0=dt_offset)
       CALL channel_halt(substr, status)
    END SELECT

    CALL end_message_bi(modstr,'DETERMINE OFFSET',substr)

  END SUBROUTINE orbit_init_coupling
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE orbit_global_start

    USE messy_main_timer,         ONLY: time_days, add_date, current_date &
                                      , print_date_components
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_constants_mem, ONLY: STRLEN_ULONG

    IMPLICIT NONE
    INTRINSIC :: INT

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'orbit_global_start'
    TYPE (time_days)            :: rad_date
    TYPE (time_days)            :: shf_date
    INTEGER                     :: status
    CHARACTER(len=STRLEN_ULONG) :: str1

    ! orbit parameters at current time step
    CALL calc_orbit_param(status, current_date, cdisse, dec, ra &
         , cossza, cosszac, rdayl)

    ! orbit parameters with offset in time
    shf_date = current_date
    CALL add_date(0,INT(dt_offset),shf_date)
    CALL calc_orbit_param(status, shf_date, cdisse_off, dec_off, ra_off &
         , cossza_off, cosszac_off, rdayl_off)

    IF (status /=0) &
         CALL error_bi('calc_orbit_param reported an error',substr)

    calc4rad: IF (l_calc4rad) THEN
       trigger_now: IF (dt_offset_rad >= 0.0_dp) THEN
          !--  Compute orbital parameters for full radiation time step.
          rad_date = current_date
          CALL add_date(0,INT(dt_offset_rad),rad_date)

          IF (.NOT. l_quiet) THEN ! op_pj_20141017
             CALL print_date_components(rad_date, status, mess=str1)
             IF (p_parallel_io) &
                  WRITE(*,*) &
                  TRIM('Orbit parameters (for radiation) calculated for : ')//&
                  &' '// TRIM(str1)
          END IF                  ! op_pj_20141017

          CALL calc_orbit_param(status, rad_date, cdissem, decm, ram &
               , cosszam, cosszacm, rdaylm)
          IF (status /=0) &
               CALL error_bi( &
               'calc_orbit_param reported an error (radiation time step)' &
               ,substr)

       END IF trigger_now
    END IF calc4rad

  END SUBROUTINE orbit_global_start
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE orbit_free_memory

    IMPLICIT NONE

    IF (l_alloc_dt_offset) THEN
       DEALLOCATE(dt_offset)
    END IF
    NULLIFY(dt_offset)

  END SUBROUTINE orbit_free_memory
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE calc_orbit_param(status, date, zcdisse, zdec, zra &
       , zcossza, zcosszac, zrdayl)

    USE messy_main_timer,         ONLY: time_days, timer_get_date, date_get &
                                      , julian_day, YearLength, MonthLength
    USE messy_main_blather_bi,    ONLY: warning_bi
    USE messy_main_constants_mem, ONLY: api => pi, OneDay
    USE messy_orbit,              ONLY: l_orbvsop87, orbit, solang

    USE messy_main_grid_def_bi,   ONLY: sinlat_2d, coslat_2d    &
                                      , sinlon_2d, coslon_2d
    USE messy_main_grid_def_mem_bi,ONLY: nproma, npromz, ngpblks

    IMPLICIT NONE
    INTRINSIC :: COS, MOD, REAL, SIN

    ! I/O
    INTEGER        , INTENT(OUT) :: status
    TYPE(time_days), INTENT(IN)  :: date
    REAL(DP),        INTENT(OUT) :: zcdisse
    REAL(DP),        INTENT(OUT) :: zdec
    REAL(DP),        INTENT(OUT) :: zra
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: zcossza
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: zcosszac
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: zrdayl

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'calc_orbit_param'
    !
    REAL(dp) :: czen1        !   *orbital parameters.
    REAL(dp) :: czen2
    REAL(dp) :: czen3
    !
    INTEGER  :: imlen, imlenh, itime   ! perpetual
    REAL(DP) :: zclock
    REAL(DP) :: zdoy
    REAL(DP) :: zdis, zgha
    REAL(DP) :: zvemar0, zyearve, zdoyve
    INTEGER  :: yr, mo, dy, hr, mn, se
    REAL(DP) :: julday, julfrac
    REAL(DP) :: julday_perp, julfrac_perp
    REAL(DP) :: julday_pal, julfrac_pal
    INTEGER  :: csecond, cday
    ! FOR SOALR ZENITH ANGLE AND RELATIVE DAY LENGTH CALCULATION
    INTEGER     :: jjrow, zproma
    REAL(dp)    :: zsin(nproma)
    REAL(dp)    :: ztim1(nproma), ztim2(nproma), ztim3(nproma)

    ! INIT
    status = 0

    !--  Compute middle of the month for perpetual runs
    IF (nmonth.ne.0) THEN
       imlen  = MonthLength(1987,nmonth)
       imlenh = imlen/2 + 1
       IF (MOD(imlen,2).EQ.0) THEN
          itime = 0
       ELSE
          itime = 12
       END IF
       julday_perp  = julian_day(REAL(imlenh,dp),nmonth,1987)
       julfrac_perp = REAL(itime,dp) / OneDay
    END IF

    CALL date_get(date, cday, csecond)
    zclock =  2.0_dp*api*REAL(csecond,dp)/OneDay
    CALL timer_get_date(status, date, yr, mo, dy, hr, mn, se)
    julday  = julian_day(REAL(dy,dp),mo,yr)
    julfrac = REAL(hr*3600+mn*60+se, DP) / OneDay
    
    ! PCMDI-Orbit
    IF (.NOT. l_orbvsop87) THEN
       
       CALL timer_get_date(status, 'start', yr, mo, dy, hr, mn, se)
       julday_pal  = julian_day(REAL(1,dp),1,1900)
       julfrac_pal = 0._dp
       IF(nmonth == 0) THEN
          zdoy = (julday + julfrac) - (julday_pal + julfrac_pal)
       ELSE                   ! perpetual month
          zdoy = (julday_perp + julfrac_perp) &
               - (julday_pal  + julfrac_pal)
       END IF
       
       ! Determine day of vernal equinox in March 1900
       
       zvemar0 = 78.41_dp - 0.0078_dp*REAL((1900-1987),DP) &
            + 0.25_dp*REAL(MOD(1900,4),DP)
       zyearve = zdoy + YearLength() - zvemar0
       zdoyve = MOD(zyearve/YearLength(),1.0_dp)*2.*api
       CALL orbit (zdoyve, zcdisse, zdec, zra, status)
       IF (status /=0) &
            CALL warning_bi('orbit_pcmdi reported an error',substr)
       !
       ! vsop87-orbit
    ELSE

       zdoy = julday + julfrac
       CALL orbit (zdoy, zra, zdec, zdis, zgha)
       zcdisse = 1.0_dp/zdis**2

    END IF
    
    czen1 = SIN(zdec)
    czen2 = COS(zdec)*COS(zclock)
    czen3 = COS(zdec)*SIN(zclock)

    ! CALCULATE SOLAR ZENIT ANGLE AND RELATIVE DAY LENGTH
    DO jjrow=1, ngpblks
#ifndef CESM1
       IF (jjrow == ngpblks) THEN
          zproma = npromz
       ELSE
          zproma = nproma
       END IF
#else
       zproma = npromz(jjrow)
#endif
       zsin (1:zproma) = sinlat_2d(1:zproma,jjrow)
       ztim1(1:zproma) =  czen1*zsin(1:zproma)
       ztim2(1:zproma) = -czen2*coslat_2d(1:zproma,jjrow)
       ztim3(1:zproma) =  czen3*coslat_2d(1:zproma,jjrow)

       CALL solang(coslon_2d(1:zproma,jjrow) &
            , sinlon_2d(1:zproma,jjrow) &
            , ztim1(1:zproma), ztim2(1:zproma), ztim3(1:zproma) &
            , zcossza(1:zproma,jjrow), zcosszac(1:zproma,jjrow)   &
            , zrdayl(1:zproma,jjrow))
    END DO

  END SUBROUTINE calc_orbit_param
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE orbit_read_nml_cpl(status, iou)

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    ! ECHAM5/MESSy

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status     ! error status
    INTEGER,          INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ l_quiet, c_rad_offset, c_offset

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'orbit_read_nml'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1
    
    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist
    
    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE orbit_read_nml_cpl
  ! ---------------------------------------------------------------------------

! **********************************************************************
END MODULE messy_orbit_si
! **********************************************************************
