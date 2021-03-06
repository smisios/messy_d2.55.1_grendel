#include "messy_main_ppd_bi.inc"
! ***********************************************************************
MODULE messy_orogw_si
! ***********************************************************************

#if defined(ECHAM5) || defined(CESM1)

  ! MESSy-SMIL FOR SUBMODEL OROGRAPHIC GRAVITY WAVES
  ! FORMERLY KNOWN AS SSODRAG AND SSORTNS
  !
  !
  ! Description:
  !
  ! SSO drag following LOTT & MILLER (1997)
  !
  ! Method:
  ! Authors:
  !
  ! m.miller + b.ritter   e.c.m.w.f.     15/06/1986.
  ! f.lott + m. miller    e.c.m.w.f.     22/11/1994
  ! e.manzini             mpi            05.09.2000
  !
  ! for more details see file AUTHORS
  !
  ! Adapted to EMAC2.52:
  ! R. EICHINGER, DLR Oberpfaffenhofen, 2016
  !

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      error_bi, warning_bi ! op_pj_20160618
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,         &
                                      mtend_register,           &   
                                      mtend_add_l,              &
                                      mtend_id_t,               &
                                      mtend_id_u,               &
                                      mtend_id_v
#endif
  ! SMCL
  USE messy_orogw

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! --------------------------
  ! fetched with get_channel_object in init-coupling
  REAL(dp), DIMENSION(:,:), POINTER :: oromea   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: orosig   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: orogam   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: orothe   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: oropic   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: oroval   => NULL()

  ! --------------------------
  ! channel objects
  REAL(dp), DIMENSION(:,:), POINTER :: ustrgw   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: vstrgw   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: vdisgw   => NULL()

  REAL(dp), DIMENSION(:,:,:), POINTER :: &
       gworo_du   => NULL(), &
       gworo_dv   => NULL(), &
       gworo_dt   => NULL(), &
       gwlif_du   => NULL(), &
       gwlif_dv   => NULL(), &
       gwlif_dt   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: w_gwd_kpr   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ampl_gwd    => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: mask_orogw  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: l_z         => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: zeff_fld    => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: llangle     => NULL()
  ! --------------------------

  ! CPL-NAMELIST PARAMETERS
  INTEGER :: i_orogw_scheme = 1 ! SSO gravity wave scheme, to be changed in orogw.nml.
  !                             ! 1: Default Lott & Miller 1997
                                ! 2: Other scheme, to be implemented!!

  LOGICAL :: lgwdrag = .FALSE. ! global switch, truncation dependent

#ifdef MESSYTENDENCY
  INTEGER   :: my_handle
#endif

  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: orogw_initialize
  PUBLIC :: orogw_init_memory
  PUBLIC :: orogw_init_coupling
  PUBLIC :: orogw_physc

  ! PRIVATE SUBROUTINES
!  PRIVATE :: orogw_read_nml_cpl
!  PRIVATE :: orogw_LandM_scheme

CONTAINS

  ! --------------------------------------------------------------------------
  SUBROUTINE orogw_initialize

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_bcast, p_io
    USE messy_main_tools,     ONLY: find_next_free_unit
    USE messy_main_grid_def_mem_bi, ONLY: nn

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'orogw_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

    lgwdrag = (nn > 21) ! see setphys.f90
    IF (.NOT. lgwdrag) THEN
       CALL warning_bi(&
            'GRAVITY WAVE DRAG OFF FOR TRUNCATION LOWER THAN T21 ' &
            ,substr)
       CALL end_message_bi(modstr,'INITIALISATION',substr)
       RETURN
    END IF

    ! READ CPL namelist
    IF (p_parallel_io) THEN                    ! read only on I/O-PE
       iou = find_next_free_unit(100,200)      ! find next free I/O unit
       CALL orogw_read_nml_cpl(status, iou)    ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error reading CPL namelist',substr)
    END IF

    ! BROADCAST CPL namelist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(i_orogw_scheme, p_io)

    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE orogw_initialize
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE orogw_init_memory

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL, GP_3D_MID
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'orogw_init_memory'
    INTEGER                     :: status ! error status

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
    CALL mtend_register(my_handle, mtend_id_t)
    CALL mtend_register(my_handle, mtend_id_u)
    CALL mtend_register(my_handle, mtend_id_v)
#endif

    IF (.NOT. lgwdrag) RETURN

    CALL new_channel(status, modstr, lrestreq=.TRUE., reprid = GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ustrgw', p2 = ustrgw)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ustrgw' &
         , 'long_name', c='u-gravity wave stress')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ustrgw' &
         , 'units', c='Pa')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'vstrgw', p2 = vstrgw)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vstrgw' &
         , 'long_name', c='v-gravity wave stress'                  )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vstrgw' &
         , 'units', c='Pa')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'vdisgw', p2 = vdisgw)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vdisgw' &
         , 'long_name', c='gravity wave dissipation'               )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vdisgw' &
         , 'units', c='W/m**2'   )
    CALL channel_halt(substr, status)

    ! ---------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'gworo_du', &
         p3=gworo_du, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gworo_du', 'long_name' &
         , c='u tendency due to ORO GW DRAG')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gworo_du', 'units', c='m/s^2')
    CALL channel_halt(substr, status)
    ! ---------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'gworo_dv', &
         p3=gworo_dv, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gworo_dv', 'long_name' &
         , c='v tendency due to ORO GW DRAG')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gworo_dv', 'units', c='m/s^2')
    CALL channel_halt(substr, status)
    ! ---------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'gworo_dt', &
         p3=gworo_dt, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gworo_dt', 'long_name' &
         , c='t tendency due to ORO GW DRAG')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gworo_dt', 'units', c='K/s')
    CALL channel_halt(substr, status)
    ! ---------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'gwlif_du', &
         p3=gwlif_du, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gwlif_du', 'long_name' &
         , c='u tendency due to MOUNTAIN LIFT')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gwlif_du', 'units', c='m/s^2')
    CALL channel_halt(substr, status)
    ! ---------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'gwlif_dv', &
         p3=gwlif_dv, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gwlif_dv', 'long_name' &
         , c='v tendency due to MOUNTAIN LIFT')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gwlif_dv', 'units', c='m/s^2')
    CALL channel_halt(substr, status)
    ! ---------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'gwlif_dt', &
         p3=gwlif_dt, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gwlif_dt', 'long_name' &
         , c='t tendency due to MOUNTAIN LIFT')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gwlif_dt', 'units', c='K/s')
    CALL channel_halt(substr, status)
    ! --------------------------------------------------------------------

    ! Extra diagnostics for orographic cirrus in messy_cloud_kuebbeler.f90
    ! --------------------------------------------------------------------
    ! FIX-ME: add the l_orocirrus switch from K14
    CALL new_channel_object(status, modstr,  'w_gwd_kpr', &
         p3=w_gwd_kpr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'w_gwd_kpr', 'long_name' &
         , c='vertical velocity induced by gravity waves')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'w_gwd_kpr', 'units', c='m/s')
    CALL channel_halt(substr, status)
    ! --------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'mask_orogw', p2=mask_orogw)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'mask_orogw', 'long_name' &
         , c='gridpoints where the scheme is active (=1)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'mask_orogw', 'units', c='-')
    CALL channel_halt(substr, status)
    ! --------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'l_z', p2=l_z)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'l_z', 'long_name' &
         , c='half-wavelength seen by the incident flow')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'l_z', 'units', c='m')
    CALL channel_halt(substr, status)
    ! --------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'ampl_gwd', &
         p3=ampl_gwd, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ampl_gwd', 'long_name' &
         , c='wave amplitude')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ampl_gwd', 'units', c='m')
    CALL channel_halt(substr, status)
    ! --------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'zeff_fld', p2=zeff_fld)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zeff_fld', 'long_name' &
         , c='effective mountain height above the blocked flow')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zeff_fld', 'units', c='m')
    CALL channel_halt(substr, status)
    ! --------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'llangle', p2=llangle)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'llangle', 'long_name' &
         , c='angle between low level wind and SS0 main axis')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'llangle', 'units', c='rad')
    CALL channel_halt(substr, status)
    ! --------------------------------------------------------------------

  END SUBROUTINE orogw_init_memory
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE orogw_init_coupling

    USE messy_main_channel,          ONLY: get_channel_object
    USE messy_main_channel_error_bi, ONLY: channel_halt

    IMPLICIT NONE

    INTEGER                     :: status ! error status

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'orogw_init_coupling'

    IF (.NOT. lgwdrag) RETURN

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ! -------------
    ! Get additional fields for sub-grid scale orography (SSO)
    ! op_pj_20160617: these objects can currently not be moved (from g3b)
    !                 to here, because they are initialised in ioinitial.f90
#ifdef ECHAM5
    CALL get_channel_object(status,'g3b','oromea',p2=oromea)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'g3b','orogam',p2=orogam)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'g3b','orothe',p2=orothe)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'g3b','oropic',p2=oropic)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'g3b','oroval',p2=oroval)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'g3b','orosig',p2=orosig)
    CALL channel_halt(substr, status)
#endif
#ifdef CESM1
    CALL get_channel_object(status, 'import_grid', 'SURFACEDAT_OROMEA', p2=oromea)
    CALL channel_halt(&
         substr//': object OROMEA in channel import_grid not found!', status)
    CALL get_channel_object(status, 'import_grid', 'SURFACEDAT_OROGAM', p2=orogam)
    CALL channel_halt(&
         substr//': object OROGAM in channel import_grid not found!', status)
    CALL get_channel_object(status, 'import_grid', 'SURFACEDAT_OROTHE', p2=orothe)
    CALL channel_halt(&
         substr//': object OROTHE in channel import_grid not found!', status)
    CALL get_channel_object(status, 'import_grid', 'SURFACEDAT_OROPIC', p2=oropic)
    CALL channel_halt(&
         substr//': object OROPIC in channel import_grid not found!', status)
    CALL get_channel_object(status, 'import_grid', 'SURFACEDAT_OROVAL', p2=oroval)
    CALL channel_halt(&
         substr//': object OROVAL in channel import_grid not found!', status)
    CALL get_channel_object(status, 'import_grid', 'SURFACEDAT_OROSIG', p2=orosig)
    CALL channel_halt(&
         substr//': object OROSIG in channel import_grid not found!', status)
#endif

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE orogw_init_coupling
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE orogw_physc

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'orogw_physc'

    ! The chosen SSO gravity wave scheme is called here.
    ! So far, there is only the default Lott & Miller (1997) scheme.
    ! More schemes are supposed to be implemented, so knock yourself out!

    IF (.NOT. lgwdrag) RETURN

       IF (i_orogw_scheme == 1) THEN
          ! Default ECHAM5 ssodrag/ssortns following Lott & Miller (1997)
          CALL orogw_LandM_scheme
       ELSEIF (i_orogw_scheme == 2) THEN
          ! To be implemented
          CALL error_bi('Please implement a new SSO gravity wave drag scheme',substr)
       ELSE
          CALL error_bi('No valid SSO gravity wave drag scheme chosen!',substr)
       END IF

  END SUBROUTINE orogw_physc
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE orogw_LandM_scheme


    USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev, nlevp1 &
                                        , nvclev, vct, nn
    USE messy_main_data_bi,       ONLY: tm1, tte_3d                   & ! temperature (t-dt) and tendency
                                      , um1, vom_3d                   & ! zonal wind (t-dt) and tendency
                                      , vm1, vol_3d                   & ! meridional wind (t-dt) and tendency
                                      , orostd                        & ! orographic standard deviation (m)
                                      , aphm1, apm1                   & ! half and full level pressure (t-dt)
                                      , geopot_3d                     ! geopotential above surface (t-dt)
    USE messy_main_constants_mem,   ONLY: g, cpd=>cp_air
    USE messy_main_timer,           ONLY: time_step_len

    IMPLICIT NONE


    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'orogw_LandM_scheme'
 
  ! Pointers for array arguments
  REAL(dp), DIMENSION(:), POINTER :: pmea   ! Mean Orography (m)
  REAL(dp), DIMENSION(:), POINTER :: pstd   ! SSO standard deviation (m)
  REAL(dp), DIMENSION(:), POINTER :: psig   ! SSO slope
  REAL(dp), DIMENSION(:), POINTER :: pgam   ! SSO Anisotropy
  REAL(dp), DIMENSION(:), POINTER :: pthe   ! SSO Angle
  REAL(dp), DIMENSION(:), POINTER :: ppic   ! SSO Peacks elevation (m)
  REAL(dp), DIMENSION(:), POINTER :: pval   ! SSO Valleys elevation (m)
  ! array arguments with intent(INOUT):
  ! Input 1D
  REAL(dp), DIMENSION(:), POINTER :: pustrgw    ! u-gravity wave stress
  REAL(dp), DIMENSION(:), POINTER :: pvstrgw    ! v-gravity wave stress
  REAL(dp), DIMENSION(:), POINTER :: pvdisgw    ! dissipation by gravity wave drag

  ! Local scalars:
  INTEGER :: igwd, jk, jl, ji
  LOGICAL :: l_truncation_error = .FALSE.

  ! Local arrays:
  INTEGER :: idx(kproma), itest(kproma)
  ! To save the tendencies due to orographic gravity wave drag.
  REAL(dp), DIMENSION(:,:), POINTER :: zdu_oro ! tendency due to ORO GW DRAG  (m/s)
  REAL(dp), DIMENSION(:,:), POINTER :: zdv_oro ! tendency due to ORO GW DRAG  (m/s)
  REAL(dp), DIMENSION(:,:), POINTER :: zdt_oro ! tendency due to ORO GW DRAG  (K)
  REAL(dp), DIMENSION(:,:), POINTER :: zdu_lif ! tendency due to MOUNTAIN LIFT(m/s)
  REAL(dp), DIMENSION(:,:), POINTER :: zdv_lif ! tendency due to MOUNTAIN LIFT(m/s)
  REAL(dp), DIMENSION(:,:), POINTER :: zdt_lif ! tendency due to MOUNTAIN LIFT(K)
#ifdef MESSYTENDENCY
  ! variable for tendency budget
  REAL(dp),DIMENSION(kproma,nlev)       :: lo_tte, lo_vom, lo_vol
#endif

  ! Diagnostic output for cirrus
  REAL(dp), DIMENSION(:,:), POINTER :: pgwd, pampl
  REAL(dp), DIMENSION(:),   POINTER :: pmask, plength, peff, pangle

  pmea => oromea(1:kproma,jrow)
  pstd => orostd(1:kproma,jrow)
  psig => orosig(1:kproma,jrow)
  pgam => orogam(1:kproma,jrow)
  pthe => orothe(1:kproma,jrow)
  ppic => oropic(1:kproma,jrow)
  pval => oroval(1:kproma,jrow)

  pustrgw => ustrgw(1:kproma,jrow)
  pvstrgw => vstrgw(1:kproma,jrow)
  pvdisgw => vdisgw(1:kproma,jrow)

  zdu_oro => gworo_du(1:kproma,:,jrow)
  zdv_oro => gworo_dv(1:kproma,:,jrow)
  zdt_oro => gworo_dt(1:kproma,:,jrow)
  zdu_lif => gwlif_du(1:kproma,:,jrow)
  zdv_lif => gwlif_dv(1:kproma,:,jrow)
  zdt_lif => gwlif_dt(1:kproma,:,jrow)

! Diagnostic output for cirrus
  pgwd    => w_gwd_kpr(1:kproma,:,jrow)
  pampl   => ampl_gwd(1:kproma,:,jrow)
  pmask   => mask_orogw(1:kproma,jrow)
  plength => l_z(1:kproma,jrow)
  peff    => zeff_fld(1:kproma,jrow)
  pangle  => llangle(1:kproma,jrow)

  !
  !*         1.    initialization
  !                --------------
  !  INITIALIZE CONSTANT FOR THE GWD SCHEME
  !  ON-LINE SHOULD ONLY BE CALLED AT THE MODEL START UP.

  CALL sugwd(kproma,nlev, nvclev, vct, nn, l_truncation_error)
  IF(l_truncation_error) CALL error_bi('Truncation not supported.', substr)

  pustrgw(:) = 0.0_dp
  pvstrgw(:) = 0.0_dp
  pvdisgw(:) = 0.0_dp

  zdu_oro(:,:) = 0.0_dp
  zdv_oro(:,:) = 0.0_dp
  zdt_oro(:,:) = 0.0_dp

  zdu_lif(:,:) = 0.0_dp
  zdv_lif(:,:) = 0.0_dp
  zdt_lif(:,:) = 0.0_dp

  idx(:) = 0

#ifdef MESSYTENDENCY
   lo_tte = 0._dp
   lo_vom = 0._dp
   lo_vol = 0._dp
#endif

   pgwd(:,:)  = 0.0_dp
   pampl(:,:) = 0.0_dp
   pmask(:)   = 0.0_dp
   plength(:) = 0.0_dp
   peff(:)    = 0.0_dp
   pangle(:)  = 0.0_dp

  !
  !*         2.    orographic gravity wave drag
  !                -----------------------------

  !  SELECTION  POINTS WHERE THE SCHEME IS ACTIVE

  igwd=0
  DO jl=1,kproma
     itest(jl)=0
     IF (((ppic(jl)-pmea(jl)) > gpicmea).AND.(pstd(jl) > gstd)) THEN
        itest(jl)=1
        igwd=igwd+1
        idx(igwd)=jl
     ENDIF
  ENDDO

  pmask(:) = REAL(itest(:), dp)
  !
  CALL orodrag( kproma,  kproma,   nlev,                                             &
                igwd,             idx,   itest,                                      &
                aphm1(1:kproma,1:nlevp1), apm1(1:kproma,:), geopot_3d(1:kproma,:,jrow),       &
                tm1(1:kproma,:,jrow), um1(1:kproma,:,jrow), vm1(1:kproma,:,jrow),    &
                pmea,    pstd,  psig, pgam, pthe, ppic, pval,                        &
                zdu_oro, zdv_oro, zdt_oro,                                           &
                pgwd, plength, pangle, pampl, peff,                                  &
                time_step_len                                                        )
!  !
!  !*         3.    mountain lift
!  !                --------------
!  CALL orolift( kproma,    nlev,                                                     &
!                itest,                                                               &
!                aphm1(1:kproma,1:nlevp1),   geom1(1:kproma,:),                       &
!                tm1(1:kproma,:,jrow), um1(1:kproma,:,jrow), vm1(1:kproma,:,jrow),    &
!                pmea,    pstd,    ppic,                                              &
!                zdu_lif, zdv_lif, zdt_lif,                                           &
!                delta_time, twomu_2d(1:kproma,jrow)                                  )
!
  ! STRESS FROM TENDENCIES

  DO jk = 1, nlev
!CDIR NODEP
     DO jl = 1, igwd
        ji=idx(jl)
        pustrgw(ji) = pustrgw(ji)                                       &
             +(zdu_oro(ji,jk)+zdu_lif(ji,jk))                           &
             *(aphm1(ji,jk+1)-aphm1(ji,jk))/g
        pvstrgw(ji) = pvstrgw(ji)                                       &
             +(zdv_oro(ji,jk)+zdv_lif(ji,jk))                           &
             *(aphm1(ji,jk+1)-aphm1(ji,jk))/g
        pvdisgw(ji) = pvdisgw(ji)                                       &
             +(zdt_oro(ji,jk)+zdt_lif(ji,jk))                           &
             *(aphm1(ji,jk+1)-aphm1(ji,jk))/g*cpd
     ENDDO
  ENDDO
  !
  !*         4.    total quantities
  !                ----------------

#ifndef MESSYTENDENCY
  do jk=1,nlev
!CDIR NODEP
    do jl=1,igwd
      ji=idx(jl)
      tte_3d(_RI_XYZ__(ji,jrow,jk)) = tte_3d(_RI_XYZ__(ji,jrow,jk))+zdt_oro(ji,jk)+zdt_lif(ji,jk)
      vol_3d(_RI_XYZ__(ji,jrow,jk)) = &
           vol_3d(_RI_XYZ__(ji,jrow,jk)) +zdv_oro(ji,jk)+zdv_lif(ji,jk)
      vom_3d(_RI_XYZ__(ji,jrow,jk)) = &
           vom_3d(_RI_XYZ__(ji,jrow,jk))+zdu_oro(ji,jk)+zdu_lif(ji,jk)
    enddo
  enddo
#else
  do jk=1,nlev
!CDIR NODEP
     do jl=1,igwd
        ji=idx(jl)
        lo_tte(ji,jk) =  zdt_oro(ji,jk)+zdt_lif(ji,jk)
        lo_vol(ji,jk) =  zdv_oro(ji,jk)+zdv_lif(ji,jk)
        lo_vom(ji,jk) =  zdu_oro(ji,jk)+zdu_lif(ji,jk)
     enddo
  enddo

  call mtend_add_l (my_handle, mtend_id_t, px = lo_tte)
  call mtend_add_l (my_handle, mtend_id_u, px = lo_vom)
  call mtend_add_l (my_handle, mtend_id_v, px = lo_vol)
#endif

  RETURN
  END SUBROUTINE orogw_LandM_scheme

  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------

  SUBROUTINE orogw_read_nml_cpl(status, iou)
   
    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    
    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ i_orogw_scheme

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='orogw_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    CALL start_message_bi(modstr,'couple namelist',substr)  ! log-output

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

    CALL end_message_bi(modstr,'couple namelist',substr)  ! log-output

  END SUBROUTINE orogw_read_nml_cpl

  ! --------------------------------------------------------------------------

#else

IMPLICIT NONE

#endif

! ***********************************************************************
END MODULE messy_orogw_si
! ***********************************************************************

