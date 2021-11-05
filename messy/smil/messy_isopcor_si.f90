#include "messy_main_ppd_bi.inc"

! **********************************************************************
MODULE messy_isopcor_si
! **********************************************************************

#define _EXTERNAL_LT

! TODO:
!   - units of q(i)%ptr: atoms/g/s -> *M_air/N_a -> mol/mol/s (=tendency)
!   - option to add production (tendency) to tracer(s), with scale-factor!
!     (or convert to correct flux units and use OFFEMIS ...)
!     --> OFFEMIS: Volume emissions (3D) must be in  molecules m-3 s-1
!   - strat / trop selection

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      error_bi, info_bi
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY
  USE messy_main_channel,       ONLY: REPR_UNDEF, t_chaobj_cpl
#ifdef _EXTERNAL_LT
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY
#endif
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
  USE messy_main_timer_event,   ONLY:  io_time_event, TRIG_FIRST, time_event
#ifdef MESSYTENDENCY
  ! tendency budget
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,  &
                                      mtend_get_start_l, &
                                      mtend_add_l,       &
                                      mtend_register,    &
                                      mtend_id_tracer
#endif
  USE messy_isopcor

  IMPLICIT NONE
  SAVE
  PRIVATE

  ! CPL namelist parameter
  LOGICAL, DIMENSION(NISO) :: l_iso = .TRUE.   ! which to calculate?
  TYPE(t_chaobj_cpl)       :: c_phi            ! heliospheric potential
  TYPE(t_chaobj_cpl)       :: c_tp             ! tropopause

  ! REPRESENTATION ID FOR ENERGY SPECTRUM
  INTEGER :: REPRID_ENERGY = REPR_UNDEF

  !
  ! POINTER TO OWN CHANNEL OBJECTS
  !
  ! Ec is energy corresponding to the local geomagnetic rigidity cutoff Pc
  ! [GeV] Ec of protons
  REAL(DP), DIMENSION(:,:), POINTER :: Ecp => NULL()
  ! [GeV] Ec of alpha particles
  REAL(DP), DIMENSION(:,:), POINTER :: Eca => NULL()
  ! [GeV] geomagnetic rigidity cutoff
  REAL(DP), DIMENSION(:,:), POINTER :: Pc => NULL()
  ! particle spectra [sr sec cm^2 ]^(-1)
  REAL(DP), DIMENSION(:), POINTER :: Jsa  ! alpha particles
  REAL(DP), DIMENSION(:), POINTER :: Jsp  ! protons
  ! production rates
  TYPE(PTR_3D_ARRAY), DIMENSION(NISO) :: Q
  !
  ! heliospheric potential [MV]
  LOGICAL :: l_phi_const
  REAL(DP), DIMENSION(:), POINTER :: zphi => NULL()

  ! POINTER TO COUPLED CHANNEL OBJECTS
  !
  ! TIME SERIES OF HELIOSPHERIC POTENTIAL
  REAL(DP), DIMENSION(:), POINTER :: phi => NULL() ! heliospheric potential [MV]
  REAL(DP), DIMENSION(:), POINTER :: phi_flag => NULL()  ! missing value ?
  !
  ! TROPOPAUSE
  LOGICAL :: l_tp = .FALSE.
  REAL(DP), DIMENSION(:,:), POINTER :: tp_i => NULL() ! level index
  REAL(DP), DIMENSION(:,:), POINTER :: tp_f => NULL() ! fraction of box below

#ifdef _EXTERNAL_LT
  TYPE(PTR_1D_ARRAY),           DIMENSION(:), POINTER :: dimaxe  => NULL()
  INTEGER,                      DIMENSION(:), POINTER :: dimlen  => NULL()
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: dimname => NULL()
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: dimunit => NULL()
#endif

  ! TIME CONTROL
  TYPE(io_time_event), PUBLIC :: trigger = &
                                       io_time_event (1,'steps',TRIG_FIRST,0)
  TYPE(time_event),    PUBLIC :: ev_trigger
  LOGICAL                     :: l_trigger
  ! check event with present date
  CHARACTER(LEN=*), PARAMETER :: EV_TLEV_PRES = 'present'

  ! COUPLING TO TRACER(S)
#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif
  TYPE T_TRACSRC_IO
     CHARACTER(LEN=STRLEN_MEDIUM) :: tr_basename = '' ! name of tracer
     CHARACTER(LEN=STRLEN_MEDIUM) :: tr_subname  = '' ! OPTIONAL subname of tr
     INTEGER                      :: iiso = 0         ! which isotope
     INTEGER                      :: dom = -1         ! where ?
     REAL(DP)                     :: scalf = 1.0_dp   ! scaling factor
  END type T_TRACSRC_IO
  TYPE T_TRACSRC
     TYPE(T_TRACSRC_IO)               :: io
     CHARACTER(LEN=2*STRLEN_MEDIUM+1) :: fullname
     INTEGER                          :: idt = 0 ! tracer id
  END type T_TRACSRC
  INTEGER, PARAMETER                      :: NTRACMAX = 20
  TYPE(T_TRACSRC_IO), DIMENSION(NTRACMAX) :: TRAC
  INTEGER                                 :: NTRAC
  TYPE(T_TRACSRC), DIMENSION(:), POINTER  :: XTRAC => NULL()
  INTEGER, PARAMETER :: D_BOTH  = 0
  INTEGER, PARAMETER :: D_TROP  = 1
  INTEGER, PARAMETER :: D_STRAT = 2

  ! ADJUSTMENTS
  TYPE T_ADJUST_IO
     CHARACTER(LEN=2*STRLEN_MEDIUM+1) :: s  = ''  ! sum tracer
     CHARACTER(LEN=2*STRLEN_MEDIUM+1) :: f1 = ''  ! fraction 1
     CHARACTER(LEN=2*STRLEN_MEDIUM+1) :: f2 = ''  ! fraction 2
  END type T_ADJUST_IO
  TYPE T_ADJUST
     TYPE(T_ADJUST_IO) :: io
     INTEGER           :: s_idt  = 0
     INTEGER           :: f1_idt = 0
     INTEGER           :: f2_idt = 0
  END type T_ADJUST
  INTEGER, PARAMETER                      :: NADJMAX = 20
  TYPE(T_ADJUST_IO), DIMENSION(NADJMAX)   :: TADJ
  INTEGER                                 :: NADJ
  TYPE(T_ADJUST), DIMENSION(:), POINTER   :: XTADJ => NULL()

  ! ENTRY POINTS
  PUBLIC :: isopcor_initialize
  PUBLIC :: isopcor_init_memory
  PUBLIC :: isopcor_init_coupling
  PUBLIC :: isopcor_global_start
  PUBLIC :: isopcor_physc
  PUBLIC :: isopcor_free_memory
  ! PRIVATE :: isopcor_read_nml_cpl
  ! PRIVATE :: isopcor_new_representation

CONTAINS

  ! ====================================================================
  SUBROUTINE isopcor_initialize

    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,         ONLY: find_next_free_unit
#ifdef _EXTERNAL_LT
    USE messy_main_import_lt,     ONLY: get_lookup_table
#endif
    USE messy_main_timer_bi,      ONLY: p_bcast_event, timer_event_init

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'isopcor_initialize'
    INTEGER :: iou
    INTEGER :: status
#ifdef _EXTERNAL_LT
    INTEGER :: ltrank
    INTEGER :: jr
#endif
    INTEGER :: jt

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL isopcor_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    CALL p_bcast(l_iso, p_io)
    CALL p_bcast(c_phi%cha, p_io)
    CALL p_bcast(c_phi%obj, p_io)
    CALL p_bcast_event(trigger,p_io)
    CALL p_bcast(c_tp%cha, p_io)
    CALL p_bcast(c_tp%obj, p_io)

    ! TRACER
    IF (p_parallel_io) THEN
       NTRAC = 0
       DO jt=1, NTRACMAX
          IF (TRIM(TRAC(jt)%tr_basename) == '') CYCLE
          NTRAC = NTRAC + 1
       END DO
    END IF
    CALL p_bcast(NTRAC, p_io)
    ALLOCATE(XTRAC(NTRAC))
    IF (p_parallel_io) THEN
       NTRAC = 0
       DO jt=1, NTRACMAX
          IF (TRIM(TRAC(jt)%tr_basename) == '') CYCLE
          NTRAC = NTRAC + 1
          XTRAC(NTRAC)%io%tr_basename = TRAC(jt)%tr_basename
          XTRAC(NTRAC)%io%tr_subname  = TRAC(jt)%tr_subname
          XTRAC(NTRAC)%io%iiso        = TRAC(jt)%iiso
          XTRAC(NTRAC)%io%dom         = TRAC(jt)%dom
          XTRAC(NTRAC)%io%scalf       = TRAC(jt)%scalf
       END DO
    ENDIF
    DO jt=1, NTRAC
       CALL p_bcast(XTRAC(jt)%io%tr_basename, p_io)
       CALL p_bcast(XTRAC(jt)%io%tr_subname, p_io)
       CALL p_bcast(XTRAC(jt)%io%iiso, p_io)
       CALL p_bcast(XTRAC(jt)%io%dom, p_io)
       CALL p_bcast(XTRAC(jt)%io%scalf, p_io)
       IF (TRIM(XTRAC(jt)%io%tr_subname) == '') THEN
          XTRAC(jt)%fullname = TRIM(XTRAC(jt)%io%tr_basename)
       ELSE
          XTRAC(jt)%fullname = TRIM(XTRAC(jt)%io%tr_basename)//'_'//&
               &TRIM(XTRAC(jt)%io%tr_subname)
       END IF
    END DO

    ! ADJSUTMENTS
    IF (p_parallel_io) THEN
       NADJ = 0
       DO jt=1, NADJMAX
          IF (TRIM(TADJ(jt)%s) == '') CYCLE
          IF (TRIM(TADJ(jt)%f1) == '') CYCLE
          IF (TRIM(TADJ(jt)%f2) == '') CYCLE
          NADJ = NADJ + 1
       END DO
    END IF
    CALL p_bcast(NADJ, p_io)
    ALLOCATE(XTADJ(NADJ))
    IF (p_parallel_io) THEN
       NADJ = 0
       DO jt=1, NADJMAX
          IF (TRIM(TADJ(jt)%s) == '') CYCLE
          IF (TRIM(TADJ(jt)%f1) == '') CYCLE
          IF (TRIM(TADJ(jt)%f2) == '') CYCLE
          NADJ = NADJ + 1
          XTADJ(NADJ)%io%s  = TADJ(jt)%s
          XTADJ(NADJ)%io%f1 = TADJ(jt)%f1
          XTADJ(NADJ)%io%f2 = TADJ(jt)%f2
       END DO
    END IF
    DO jt=1, NADJ
       CALL p_bcast(XTADJ(NADJ)%io%s,  p_io)
       CALL p_bcast(XTADJ(NADJ)%io%f1, p_io)
       CALL p_bcast(XTADJ(NADJ)%io%f2, p_io)
    END DO

    ! initialize event
    CALL timer_event_init(ev_trigger, trigger, &
         'isotope production rate computation', EV_TLEV_PRES)

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

#ifdef _EXTERNAL_LT
    CALL start_message_bi(modstr, 'LOOKUP TABLE INITIALISATION', substr)

    CALL get_lookup_table(status             &
         , 'yield', p4=Y, rank=ltrank        &
         , dimlen=dimlen, dimaxe=dimaxe, dimname=dimname, dimunit=dimunit)

    IF (status /= 0) CALL error_bi('yield lookup table not available',substr)
    IF (ltrank /= 4) CALL error_bi('yield lookup table has rank /= 4',substr)

    DO jr = 1, ltrank
       ! NOTE: units are NOT tested !
       SELECT CASE(TRIM(dimname(jr)))
       CASE('energy')
          IF (jr /= 1) CALL error_bi('yield LT; energy not rank 1',substr)
          NE = dimlen(jr)
          energy => dimaxe(jr)%ptr(:)
       CASE('depth')
          IF (jr /= 2) CALL error_bi('yield LT; depth not rank 2',substr)
          NH = dimlen(jr)
          depth => dimaxe(jr)%ptr(:)
       CASE('NCR')
          IF (jr /= 3) CALL error_bi('yield LT; NCR not rank 3',substr)
          NCR = dimlen(jr)
       CASE('NISO')
          IF (jr /= 4) CALL error_bi('yield LT; NISO not rank 4',substr)
          IF (dimlen(jr) /= NISO) &
               CALL error_bi('yield lokkup table has wrong'//&
               &' number of isotopes',substr)
       CASE DEFAULT
          CALL error_bi('yield lookup table with unknown axis: '//&
               &TRIM(dimname(jr)),substr)
       END SELECT
    END DO

    CALL end_message_bi(modstr, 'LOOKUP TABLE INITIALISATION', substr)
#endif

  END SUBROUTINE isopcor_initialize
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE isopcor_init_memory

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL, GP_3D_MID
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'isopcor_init_memory'
    INTEGER                     :: status
    INTEGER                     :: i

    CALL isopcor_new_representation

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ! DEFINE NEW CHANNEL
    CALL new_channel(status, modstr, lrestreq = .TRUE.) ! (!) trigger
    CALL channel_halt(substr, status)

    ! ----------------------
    CALL new_channel_object(status, modstr, 'Ecp' &
         , p2 = Ecp &
         , reprid = GP_2D_HORIZONTAL, lstatic=.TRUE. )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Ecp'   &
         , 'long_name' &
         , c='energy corresponding to the local geomagnetic rigidity'//&
         &' cutoff (protons)' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Ecp'   &
         , 'units', c='GeV' )
    CALL channel_halt(substr, status)
    ! ----------------------

    ! ----------------------
    CALL new_channel_object(status, modstr, 'Eca' &
         , p2 = Eca &
         , reprid = GP_2D_HORIZONTAL, lstatic=.TRUE. )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Eca'   &
         , 'long_name' &
         , c='energy corresponding to the local geomagnetic rigidity'//&
         &' cutoff (alpa particles)' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Eca'   &
         , 'units', c='GeV' )
    CALL channel_halt(substr, status)
    ! ----------------------

    ! ----------------------
    CALL new_channel_object(status, modstr, 'Pc' &
         , p2 = Pc &
         , reprid = GP_2D_HORIZONTAL, lstatic=.TRUE. )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Pc'   &
         , 'long_name' &
         , c='geomagnetic rigidity cutoff' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Pc'   &
         , 'units', c='GeV' )
    CALL channel_halt(substr, status)
    ! ----------------------

    ! ----------------------
    DO i=1, NISO
       !
       IF (.NOT. l_iso(i)) CYCLE
       !
       CALL new_channel_object(status, modstr, 'Q_'//TRIM(isoname(i)) &
            , p3 = Q(i)%ptr &
            , reprid = GP_3D_MID )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'Q_'//TRIM(isoname(i)) &
            , 'long_name' &
            , c='production rate of '//TRIM(isoname(i)) )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'Q_'//TRIM(isoname(i)) &
            , 'units', c='atoms/g/s' )
       CALL channel_halt(substr, status)
    END DO
    ! ----------------------

    ! ----------------------
    CALL new_channel_object(status, modstr, 'Jsa' &
         , p1 = Jsa &
         , reprid = REPRID_ENERGY )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Jsa' &
            , 'long_name' &
            , c='spectrum of alpha particles' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'Jsa' &
            , 'units', c='(sr sec cm^2 )^(-1)' )
       CALL channel_halt(substr, status)
    ! ----------------------

    ! ----------------------
    CALL new_channel_object(status, modstr, 'Jsp' &
         , p1 = Jsp &
         , reprid = REPRID_ENERGY )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Jsp' &
            , 'long_name' &
            , c='spectrum of protons' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'Jsp' &
            , 'units', c='(sr sec cm^2 )^(-1)' )
       CALL channel_halt(substr, status)
    ! ----------------------

#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, mtend_id_tracer)
#endif

  END SUBROUTINE isopcor_init_memory
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE isopcor_init_coupling

    USE messy_main_channel,          ONLY: get_channel_object &
                                         , new_channel_object &
                                         , new_attribute
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: SCALAR
    USE messy_main_tools,            ONLY: str2num
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_tracer,           ONLY: get_tracer, full2base_sub

    IMPLICIT NONE
    INTRINSIC :: TRIM, ADJUSTL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'isopcor_init_coupling'
    INTEGER :: status
    CHARACTER(LEN=20) :: rstr = ''
    INTEGER :: jt
    CHARACTER(LEN=STRLEN_MEDIUM) :: sub = ''
    CHARACTER(LEN=STRLEN_MEDIUM) :: base = ''

    CALL start_message_bi(modstr,'INITIALIZE COUPLING', substr)

    IF (TRIM(ADJUSTL(c_phi%cha)) == '#' ) THEN
       !
       ! constant phi
       !
       l_phi_const = .TRUE.
       !
       CALL new_channel_object(status, modstr, 'phi' &
         , p1 = zphi, reprid = SCALAR, lstatic=.TRUE.)
       CALL channel_halt(substr, status)
       !
       CALL str2num(c_phi%obj, zphi(1), status)
       IF (status /= 0) THEN
          CALL error_bi('read error of CPL parameter phi',substr)
       ELSE
          write(rstr,'(f12.4)') zphi(1)
          CALL info_bi('constant phi: '//TRIM(rstr)//' MV', substr)
       ENDIF
    ELSE
       !
       ! time series of phi
       !
       l_phi_const = .FALSE.
       !
       CALL get_channel_object(status &
            , TRIM(c_phi%cha), TRIM(c_phi%obj), p1=phi)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status &
            , TRIM(c_phi%cha), TRIM(c_phi%obj)//'_flg', p1=phi_flag)
       CALL info_bi('time dependent phi', substr)
       !
       CALL new_channel_object(status, modstr, 'phi' &
         , p1 = zphi, reprid = SCALAR)
       CALL channel_halt(substr, status)
       !
    END IF

    CALL new_attribute(status, modstr, 'phi' &
         , 'long_name' &
         , c='heliospheric potential' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'phi' &
         , 'units', c='MV' )
    CALL channel_halt(substr, status)

    ! TROPOPAUSE
    IF ( (TRIM(c_tp%cha) /= '') .AND. (TRIM(c_tp%obj) /= '') ) THEN
       CALL get_channel_object(status &
            , TRIM(c_tp%cha), TRIM(c_tp%obj)//'_i', p2 = tp_i)
       CALL channel_halt(substr, status)
       !
       CALL get_channel_object(status &
            , TRIM(c_tp%cha), TRIM(c_tp%obj)//'_f', p2 = tp_f)
       CALL channel_halt(substr, status)
       !
       l_tp = .TRUE.
       !
    END IF

    ! TRACER
    DO jt = 1, NTRAC
       IF (.NOT. l_iso(XTRAC(jt)%io%iiso)) &
            CALL error_bi('corresponding prod.rate for tracer '//&
            &TRIM(XTRAC(jt)%fullname)//'not enabled',substr)

       IF ( (.NOT. l_tp) .AND. (XTRAC(jt)%io%dom /= D_BOTH) ) &
            CALL error_bi('tropopause information for tracer '//&
            &TRIM(XTRAC(jt)%fullname)//'required',substr)

       CALL get_tracer(status, GPTRSTR       &
            , TRIM(XTRAC(jt)%io%tr_basename) &
            , TRIM(XTRAC(jt)%io%tr_subname)  &
            , idx = XTRAC(jt)%idt            &
            )
       CALL tracer_halt(substr, status)
    END DO

    ! ADJUSTMENTS
    DO jt=1, NADJ
       CALL full2base_sub(status, XTADJ(jt)%io%s, base, sub)
       CALL tracer_halt(substr, status)
       CALL get_tracer(status, GPTRSTR &
            , TRIM(base), TRIM(sub), idx = XTADJ(jt)%s_idt)
       CALL tracer_halt(substr, status)

       CALL full2base_sub(status, XTADJ(jt)%io%f1, base, sub)
       CALL tracer_halt(substr, status)
       CALL get_tracer(status, GPTRSTR &
            , TRIM(base), TRIM(sub), idx = XTADJ(jt)%f1_idt)
       CALL tracer_halt(substr, status)

       CALL full2base_sub(status, XTADJ(jt)%io%f2, base, sub)
       CALL tracer_halt(substr, status)
       CALL get_tracer(status, GPTRSTR &
            , TRIM(base), TRIM(sub), idx = XTADJ(jt)%f2_idt)
       CALL tracer_halt(substr, status)
    END DO

    CALL end_message_bi(modstr,'INITIALIZE COUPLING', substr)

  END SUBROUTINE isopcor_init_coupling
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE isopcor_global_start

    USE messy_main_grid_def_mem_bi, ONLY: nproma, npromz, ngpblks
    USE messy_main_grid_def_bi,     ONLY:philat_2d, philon_2d
    USE messy_main_timer,           ONLY: lstart, lresume, current_date
    USE messy_main_timer_bi,      ONLY: event_state

    IMPLICIT NONE
    INTRINSIC :: MERGE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'isopcor_global_start'
    INTEGER  :: jjrow, zproma
    REAL(dp) :: lon(nproma)

    l_trigger = event_state(ev_trigger, current_date) .OR. lstart

    once: IF (lstart .OR. lresume) THEN

       DO jjrow = 1, ngpblks
#ifndef CESM1
          IF (jjrow == ngpblks) THEN
             zproma = npromz
          ELSE
             zproma = nproma
          END IF
#else
          zproma = npromz(jjrow)
#endif

          lon(1:zproma) = MERGE(philon_2d(1:zproma, jjrow) - 360.0_dp  &
               , philon_2d(1:zproma, jjrow)             &
               , philon_2d(1:zproma, jjrow) > 180.0_dp)

          ! Notes:
          !  - Pc is only output for diagnostic purposes
          !  - it needs to be checked, which range of LON is accepted
          !    [-180,180] or [0,360] or both?
          !  - results are NOT dependent on altitude!
          !  - the soubroutine is ELEMENTAL and can directly be used for arrays

          CALL geomag( &
               lon(1:zproma)                & ! INTENT(IN)
               , philat_2d(1:zproma, jjrow) & ! INTENT(IN)
               , Ecp(1:zproma, jjrow)       & ! INTENT(OUT)
               , Eca(1:zproma, jjrow)       & ! INTENT(OUT)
               , Pc(1:zproma, jjrow)        & ! INTENT(OUT)
               )

       END DO

    END IF once

    ! ONLY FOR TIME SERIES ...
    IF (ASSOCIATED(phi_flag)) THEN
       IF (phi_flag(1) < 1.0) &
            CALL error_bi('mssing value in time series of phi',substr)
    END IF

    IF (.NOT. l_trigger) RETURN

    ! Notes:
    !  - unit conversion of PHI: MV->GV
    !  - force-field approximation: the spectrum depends only on the
    !    heliospheric potential PHI
    IF (.NOT. l_phi_const) zphi(1) = phi(1)
    CALL gcr_spectrum(zphi(1) * 1.e-3_dp, Jsp(:), Jsa(:))

  END SUBROUTINE isopcor_global_start
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE isopcor_physc

    USE messy_main_data_bi,         ONLY: press_3d
    USE messy_main_grid_def_mem_bi, ONLY: jrow, nlev, kproma
    USE messy_main_constants_mem,  ONLY: M_air, N_A
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi,  ONLY: qxtte, qxtm1
#else

#endif
    USE messy_main_timer,          ONLY: ztmst=>time_step_len

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: jp, jk, i
    ! atmospheric depth [g cm^(-2)]  is (p/1013.25) * 1033.0; p in [hPa]
    REAL(DP), PARAMETER          :: p2d = 1033.0_dp/101325.0_dp
    REAL(DP), DIMENSION(NE,NISO) :: Fa, Fp
    REAL(DP)                     :: depth ! [g/cm3]
    REAL(DP)                     :: QISO
    INTEGER                      :: jt, idt
    ! conversion from atoms/(g s) into mol/mol/s:  *(g/mol)/(molec/mol)
    REAL(DP), PARAMETER          :: conv = M_air/N_A
    REAL(DP), DIMENSION(kproma,nlev) :: zflag
    REAL(DP), DIMENSION(kproma,nlev) :: zxtte
    !
    REAL(DP), DIMENSION(kproma,nlev) :: sm1, ste, f1m1, f1te, f2m1, f2te

    IF (l_trigger) THEN
       ! update production rate
       DO jk=1, nlev
          DO jp=1, kproma
             depth = press_3d(_RI_XYZ__(jp,jrow,jk)) * p2d
             CALL yield(depth, Jsp(:), Jsa(:), Fp(:,:), Fa(:,:))
             DO i=1, NISO
                IF (.NOT. l_iso(i)) CYCLE
                CALL integrate(Fp(:,i), Fa(:,i)            &
                     , Ecp(kproma,jrow), Eca(kproma, jrow) &
                     , QISO)
                Q(i)%ptr(_RI_XYZ__(jp,jrow,jk)) = QISO
             END DO
          END DO
       END DO

    END IF

    ! TRACER
    DO jt=1, NTRAC

       zflag(:,:) = 0.0_dp
       SELECT CASE(XTRAC(jt)%io%dom)
       CASE(D_BOTH)
          zflag(:,:) = 1.0_dp
       CASE(D_TROP)
          DO jk=1, nlev
             DO jp=1, kproma
                IF (jk >  NINT(tp_i(jp,jrow))) zflag(jp,jk) = 1.0_dp
                IF (jk == NINT(tp_i(jp,jrow))) zflag(jp,jk) = tp_f(jp,jrow)
             END DO
          END DO
       CASE(D_STRAT)
          DO jk=1, nlev
             DO jp=1, kproma
                IF (jk <  NINT(tp_i(jp,jrow))) zflag(jp,jk) = 1.0_dp
                IF (jk == NINT(tp_i(jp,jrow))) zflag(jp,jk) = 1.0_dp &
                     - tp_f(jp,jrow)
             END DO
          END DO
       END SELECT

       zxtte(:,:) = Q(XTRAC(jt)%io%iiso)%ptr(_RI_XYZ__(1:kproma,jrow,:)) &
            * zflag(:,:) * conv * XTRAC(jt)%io%scalf

       idt = XTRAC(jt)%idt

#ifndef MESSYTENDENCY
       qxtte(_RI_X_ZN_(1:kproma,:,idt)) = qxtte(_RI_X_ZN_(1:kproma,:,idt)) &
            + zxtte(:,:)
#else
       CALL mtend_add_l(my_handle, idt, px=zxtte)
#endif

    END DO

    ! ADJUST TRACER TENDENCIES TO FORCE S = F1 + F2
    DO jt=1, NADJ
       ste(:,:)  = 0.0_dp
       f1te(:,:) = 0.0_dp
       f2te(:,:) = 0.0_dp
#ifndef MESSYTENDENCY
       idt = XTADJ(jt)%s_idt
       sm1(:,:) = qxtm1(_RI_X_ZN_(1:kproma,:,idt)) + &
            qxtte(_RI_X_ZN_(1:kproma,:,idt))*ztmst

       idt = XTADJ(jt)%f1_idt
       f1m1(:,:) = qxtm1(_RI_X_ZN_(1:kproma,:,idt)) + &
            qxtte(_RI_X_ZN_(1:kproma,:,idt))*ztmst

       idt = XTADJ(jt)%f2_idt
       f2m1(:,:) = qxtm1(_RI_X_ZN_(1:kproma,:,idt)) + &
            qxtte(_RI_X_ZN_(1:kproma,:,idt))*ztmst
#else
       CALL mtend_get_start_l(XTADJ(jt)%s_idt,  v0=sm1)
       CALL mtend_get_start_l(XTADJ(jt)%f1_idt, v0=f1m1)
       CALL mtend_get_start_l(XTADJ(jt)%f2_idt, v0=f2m1)
#endif
       DO jp=1, kproma
          CALL adj_tend(sm1(jp,:), ste(jp,:), &
               f1m1(jp,:), f1te(jp,:), f2m1(jp,:), f2te(jp,:), ztmst)
       END DO
#ifndef MESSYTENDENCY
       idt = XTADJ(jt)%f1_idt
       qxtte(_RI_X_ZN_(1:kproma,:,idt)) = qxtte(_RI_X_ZN_(1:kproma,:,idt)) &
            + f1te(:,:)
       idt = XTADJ(jt)%f2_idt
       qxtte(_RI_X_ZN_(1:kproma,:,idt)) = qxtte(_RI_X_ZN_(1:kproma,:,idt)) &
            + f2te(:,:)
#else
       CALL mtend_add_l(my_handle, XTADJ(jt)%f1_idt, px=f1te)
       CALL mtend_add_l(my_handle, XTADJ(jt)%f2_idt, px=f2te)
#endif
    END DO

  END SUBROUTINE isopcor_physc
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE isopcor_free_memory

    IMPLICIT NONE

#ifdef _EXTERNAL_LT
    NULLIFY(dimaxe)
    NULLIFY(dimlen)
    NULLIFY(dimname)
    NULLIFY(dimunit)
    NULLIFY(Y)
#endif
    IF (ASSOCIATED(XTRAC)) THEN
       DEALLOCATE(XTRAC) ; NULLIFY(XTRAC)
    END IF
    IF (ASSOCIATED(XTADJ)) THEN
       DEALLOCATE(XTADJ) ; NULLIFY(XTADJ)
    END IF

  END SUBROUTINE isopcor_free_memory
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE isopcor_read_nml_cpl(status, iou)

    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check &
                                      , read_nml_close
    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'isopcor_read_nml_cpl'

    NAMELIST /CPL/ c_phi, c_tp, trigger, l_iso, trac, tadj

    ! LOCAL
    LOGICAL :: lex      ! file exists ?
    INTEGER :: fstat    ! file status

    status = 1

    ! DEFAULTS
    l_iso(:) = .TRUE.

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE isopcor_read_nml_cpl
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE isopcor_new_representation

    USE messy_main_channel_dimensions,   ONLY: new_dimension            &
                                             , add_dimension_variable   &
                                             , add_dimension_variable_att
    USE messy_main_channel_error_bi,     ONLY: channel_halt
    USE messy_main_channel_bi,           ONLY: DC_BC
    USE messy_main_channel_repr,         ONLY: new_representation, AUTO &
                                             , set_representation_decomp &
                                             , IRANK, PIOTYPE_COL

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'isopcor_new_representation'
    INTEGER :: status
    INTEGER :: DIMID_SPEC
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    CALL start_message_bi(modstr, 'NEW REPRESENTATIONS', substr)

    ! NEW DIMENSION
    CALL new_dimension(status, DIMID_SPEC, 'energy', NE)
    CALL channel_halt(substr, status)

    CALL add_dimension_variable(status, 'energy', 'energy', energy)
    CALL channel_halt(substr, status)

    CALL add_dimension_variable_att(status, 'energy', 'energy', &
         'long_name', c='energy bins')
    CALL channel_halt(substr, status)

    CALL add_dimension_variable_att(status, 'energy', 'energy', &
         'units', c='GeV/nucleon')
    CALL channel_halt(substr, status)

    ! NEW REPRESENTATION
    CALL new_representation(status, REPRID_ENERGY, 'ENERGY' &
         , rank = 1, link = 'x---', dctype = DC_BC                   &
         , dimension_ids = (/ DIMID_SPEC /) &
         , ldimlen       = (/ AUTO  /)     &
         , axis = 'N---'                   &
         )
    CALL channel_halt(substr, status)

    ! -----------------------------------------------------------

    nseg = 1
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:) = 1
    meml(:,:) = 1
    memu(:,:) = 1

    start(:,1) = 1
    cnt(:,1)   = NE
    meml(:,1)  = 1
    memu(:,1)  = NE

    CALL set_representation_decomp(status, REPRID_ENERGY &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    CALL end_message_bi(modstr, 'NEW REPRESENTATIONS', substr)

  END SUBROUTINE isopcor_new_representation
  ! ====================================================================

! **********************************************************************
END MODULE messy_isopcor_si
! **********************************************************************
