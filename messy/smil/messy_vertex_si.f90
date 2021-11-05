MODULE messy_vertex_si

!  AUTHOR:   Huug Ouwersloot, MPI Chemie, Mainz
!            Last modified - July 2016
!
!  Currently only suitable for ECHAM/MESSy.
!  Input/output variables are missing for other core models.
!  Where applicable, this is indicated in the comments for future expansion.
!
!  For EMAC either VDIFF or VERTEX should be enabled

  USE messy_vertex
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,  mtend_register,  &
                                      mtend_get_start_l, mtend_id_t,      &
                                      mtend_id_q,        mtend_id_xl,     &
                                      mtend_id_xi,       mtend_id_u,      &
                                      mtend_id_v,        mtend_id_tracer, &
                                      mtend_add_l
#endif
  USE messy_main_data_bi,       ONLY: ledith   ! op_pj_20180724
  USE messy_main_channel,       ONLY: t_chaobj_cpl  ! ju_te_20190212

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! CPL-NAMELIST PARAMETERS
  TYPE(t_chaobj_cpl)               :: imp_lai     ! ju_te_20190212

  !Locally defined
  !3D
  REAL(dp), DIMENSION(:,:,:), POINTER :: tke       => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: tkem      => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: tkem1     => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: gweddy    => NULL() ! ka_sv_20171219

  !2D
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfl      => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfli     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfll     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ahflw     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfs      => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfsi     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfsl     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfsw     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: az0hi     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: az0hl     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: az0hw     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: az0i      => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: az0l      => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: az0w      => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: dew2      => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: evap      => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: evapi     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: evapl_2d  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: evapot_2d => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: evapw     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: temp2     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ustr      => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ustri     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ustrl     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ustrw     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: vdis      => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: vstr      => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: vstri     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: vstrl     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: vstrw     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: wet_tmp   => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: wind10    => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: wind10w   => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: lai       => NULL() ! ju_te_20190212
!
! op_pj_20181105+
! In the old ECHAM5 times, a lot of accumulated fields were present in many
! parts of the model. These were additional variables which had to be
! accumulated (summed in every time step) "manually". With the introduction of
! CHANNEL these variables became obsolete, because it is possible to output also
! the average (over the output time interval) of each variable (channel object)
! by simply requesting it in the channel.nml &CTRL namelist. As a consequence,
! in the course of modularisation of the ECHAM5 physics package, these
! "accumulated" variables have been eliminated. However, in two parts of the
! physics package the accumulated variables have not just been accumulated, but
! in the course of accumulation additionally "masked" by land, water and ice
! masks: radiation and vdiff. In the submodel RAD this was treated correctly
! (see variables with "(masked)" in the long_name). In E5VDIFF (and likewise in
! VERTEX), however, the variables got lost and need to be re-introduced as
! "(masked)" variables. The (time averaged) output of those "masked" variables
! is required to correctly calculate the flux correction required for the
! submodel MLO (mixed layer ocean).
!
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfslac => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfswac => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfsiac => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfllac => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ahflwac => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfliac => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: evaplac => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: evapwac => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: evapiac => NULL()
! op_pj_20181105-

  !Externally defined
  !3D
  REAL(dp), DIMENSION(:,:,:), POINTER :: emter     => NULL()

  !2D
  REAL(dp), DIMENSION(:,:),   POINTER :: albedo    => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: alsol     => NULL()
#ifdef ECHAM5
  REAL(dp), DIMENSION(:,:),   POINTER :: az0       => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ocu       => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ocv       => NULL()
#endif

#ifndef MESSYTENDENCY
  REAL(dp), POINTER, DIMENSION(:,:) :: ptte   => NULL()
  REAL(dp), POINTER, DIMENSION(:,:) :: pqte   => NULL()
  REAL(dp), POINTER, DIMENSION(:,:) :: pxlte  => NULL()
  REAL(dp), POINTER, DIMENSION(:,:) :: pxite  => NULL()
  REAL(dp), POINTER, DIMENSION(:,:) :: pvol   => NULL()
  REAL(dp), POINTER, DIMENSION(:,:) :: pvom   => NULL()
#endif
  REAL(dp), POINTER, DIMENSION(:,:) :: pgeom1 => NULL()

#ifdef MESSYTENDENCY
  INTEGER                             :: my_handle
#endif

  PUBLIC  :: vertex_initialize, vertex_init_memory, vertex_init_coupling
  PUBLIC  :: vertex_global_start, vertex_vdiff, vertex_global_end, vertex_free_memory

  CONTAINS

!==============================================================================
  SUBROUTINE vertex_initialize

    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_tools,         ONLY: find_next_free_unit

    IMPLICIT NONE

    ! ju_te_20180626+
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='vertex_initialize'
    INTEGER                     :: status
    INTEGER                     :: iou    ! I/O unit

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL vertex_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast(irstom,p_io)                ! Flag for choose rstom calculation
    ! ju_te_20180626-
    CALL p_bcast(izwet,p_io) ! ju_te_20180921  Flag for switch on reduction factors applied to evapotranspiration
    CALL p_bcast(ifws,p_io)  ! ju_te_20190920  Flag for new parametrization of fws

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

    zh2  = 1._dp/3._dp-2._dp*zm1/zm2       ! Initialize constant to calculate S_N_H - gamma_1
    zshn = 3._dp*zh1*zh2*SQRT(2._dp)       ! Initialize neutral exchange coefficient for heat, S_N_H, used to calculate S_H
    zsmn = 3._dp*zm1*(zh2-zm4)*SQRT(2._dp) ! Initialize neutral exchange coefficient for momentum, S_N_M, used to calculate S_M

    ! ju_te_20190212+
    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL VERTEX CORE ROUTINE:
       CALL vertex_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('READ ERROR in &CPL namelist ',substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast(imp_lai%CHA,p_io)
    CALL p_bcast(imp_lai%OBJ,p_io)
    ! ju_te_20190212-

  END SUBROUTINE vertex_initialize
!==============================================================================

!==============================================================================
  SUBROUTINE vertex_init_memory
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL, GP_3D_MID
    USE messy_main_channel,          ONLY: new_channel, new_channel_object, &
                                           new_attribute

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'vertex_init_memory'
    INTEGER                     :: status

#ifdef MESSYTENDENCY
    CALL mtend_register (my_handle,mtend_id_t)
    CALL mtend_register (my_handle,mtend_id_q)
    CALL mtend_register (my_handle,mtend_id_xl)
    CALL mtend_register (my_handle,mtend_id_xi)
    CALL mtend_register (my_handle,mtend_id_u)
    CALL mtend_register (my_handle,mtend_id_v)
    CALL mtend_register (my_handle,mtend_id_tracer)
#endif

    ! DEFINE NEW CHANNEL
    CALL new_channel(status, modstr, lrestreq=.TRUE., reprid = GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)

    ! 2D OBJECTS
    CALL new_channel_object(status, modstr, 'ahfl',      p2 = ahfl)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfl'        &
         , 'long_name', c='latent heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfl',      'units', c='W/m**2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ahfli',     p2 = ahfli)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfli'       &
         , 'long_name', c='latent heat flux over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfli',     'units', c='W/m**2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ahfll',     p2 = ahfll)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfll'       &
         , 'long_name', c='latent heat flux over land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfll',     'units', c='W/m**2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ahflw',     p2 = ahflw)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahflw'       &
         , 'long_name', c='latent heat flux over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahflw',     'units', c='W/m**2')
    CALL channel_halt(substr, status)

    ! op_pj_20181105+
    CALL new_channel_object(status, modstr, 'ahfliac', p2 = ahfliac)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfliac' &
         , 'long_name', c='latent heat flux over ice (masked)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfliac' &
         , 'units', c='W/m**2'   )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ahflwac', p2 = ahflwac)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahflwac' &
         , 'long_name', c='latent heat flux over water (masked)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahflwac' &
         , 'units', c='W/m**2'   )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ahfllac', p2 = ahfllac)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfllac' &
         , 'long_name', c='latent heat flux over land (masked)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfllac' &
         , 'units', c='W/m**2'   )
    CALL channel_halt(substr, status)
    ! op_pj_20181105-

    CALL new_channel_object(status, modstr, 'ahfs',      p2 = ahfs)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfs'        &
         , 'long_name', c='sensible heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfs',      'units', c='W/m**2'   )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ahfsi',     p2 = ahfsi)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfsi'       &
         , 'long_name', c='sensible heat flux over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfsi',     'units', c='W m-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ahfsl',     p2 = ahfsl)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfsl'       &
         , 'long_name', c='sensible heat flux over land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfsl',     'units', c='W m-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ahfsw',     p2 = ahfsw)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfsw'       &
         , 'long_name', c='sensible heat flux over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfsw',     'units', c='W m-2')
    CALL channel_halt(substr, status)

    ! op_pj_20181105+
    CALL new_channel_object(status, modstr,  'ahfsiac', p2=ahfsiac)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfsiac', &
         'long_name', c='sensible heat flux over ice (masked)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfsiac', 'units', c='W m-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'ahfswac', p2=ahfswac)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfswac', &
         'long_name', c='sensible heat flux over water (masked)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfswac', 'units', c='W m-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'ahfslac', p2=ahfslac)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfslac', &
         'long_name', c='sensible heat flux over land (masked)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfslac', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! op_pj_20181105-

    CALL new_channel_object(status, modstr, 'az0hi',     p2 = az0hi)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0hi'       &
         , 'long_name', c='roughness length for heat over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0hi',     'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'az0hl',     p2 = az0hl)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0hl'       &
         , 'long_name', c='roughness length for heat over land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0hl',     'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'az0hw',     p2 = az0hw)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0hw'       &
         , 'long_name', c='roughness length for heat over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0hw',     'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'az0i',      p2 = az0i)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0i'        &
         , 'long_name', c='roughness length over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0i',      'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'az0l',      p2 = az0l)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0l'        &
         , 'long_name', c='roughness length over land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0l',      'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'az0w',      p2 = az0w)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0w'        &
         , 'long_name', c='roughness length over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0w',      'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'dew2',      p2 = dew2)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dew2'        &
         , 'long_name', c='2m dew point temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dew2',      'units', c='K'        )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'evap',      p2 = evap)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evap'        &
         , 'long_name', c='evaporation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evap',      'units', c='kg/m**2s' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'evapi',     p2 = evapi)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapi',      &
         'long_name', c='evaporation over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapi',     'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'evapl_2d',  p2 = evapl_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapl_2d',   &
         'long_name', c='total evaporation, including sublimation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapl_2d',  'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'evapot_2d', p2 = evapot_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapot_2d',  &
         'long_name', c='Potential evaporation/sublimation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapot_2d', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'evapw',     p2 = evapw)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapw',      &
         'long_name', c='evaporation over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapw',     'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    ! op_pj_20181105+
    CALL new_channel_object(status, modstr,  'evapiac', p2=evapiac)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapiac', &
         'long_name', c='evaporation over ice (masked)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapiac', 'units', c='')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'evapwac', p2=evapwac)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapwac', &
         'long_name', c='evaporation over water (masked)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapwac', 'units', c='')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'evaplac', p2=evaplac)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evaplac', &
         'long_name', c='evaporation over land including sublimation (masked)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evaplac', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)
    ! op_pj_20181105-

    CALL new_channel_object(status, modstr, 'temp2',     p2 = temp2)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'temp2'       &
         , 'long_name', c='2m temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'temp2',     'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ustr',      p2 = ustr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ustr'        &
         , 'long_name', c='u-stress')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ustr',      'units', c='Pa')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ustri',     p2 = ustri)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ustri'       &
         , 'long_name', c='zonal wind stress over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ustri',     'units', c='Pa')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ustrl',     p2 = ustrl)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ustrl'       &
         , 'long_name', c='zonal wind stress over land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ustrl',     'units', c='Pa')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ustrw',     p2 = ustrw)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ustrw'       &
         , 'long_name', c='zonal wind stress over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ustrw',     'units', c='Pa')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'vdis',      p2 = vdis)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vdis'        &
         , 'long_name', c='boundary layer dissipation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vdis',      'units', c='W/m**2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'vstr',      p2 = vstr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vstr'        &
         , 'long_name', c='v-stress')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vstr',      'units', c='Pa')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'vstri',     p2 = vstri)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vstri'       &
         , 'long_name', c='meridional wind stress over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vstri',     'units', c='Pa')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'vstrl',     p2 = vstrl)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vstrl'       &
         , 'long_name', c='meridional wind stress over land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vstrl',     'units', c='Pa')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'vstrw',     p2 = vstrw)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vstrw'       &
         , 'long_name', c='meridional wind stress over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vstrw',     'units', c='Pa')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'wet_tmp',   p2 = wet_tmp)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wet_tmp'     &
         , 'long_name', c='temp result of zwet evapotrans in e5vdiff')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wet_tmp',   'units', c='s m-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'wind10',    p2 = wind10)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wind10',     &
         'long_name', c='10 m wind speed')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wind10',    'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'wind10w',   p2 = wind10w)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wind10w',    &
         'long_name', c='10 m wind speed over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wind10w',   'units', c='m s-1')
    CALL channel_halt(substr, status)

    ! 3D OBJECTS
    CALL new_channel_object(status, modstr, 'tke'    &
         , reprid = GP_3D_MID, p3 = tke)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tke'         &
         , 'long_name', c='turbulent kinetic energy (t+1)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tke',       'units', c='m**2/s**2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'tkem'   &
         , reprid = GP_3D_MID, p3 = tkem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tkem'        &
         , 'long_name', c='turbulent kinetic energy')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tkem',      'units', c='m**2/s**2')
    CALL channel_halt(substr, status)
    tkem(:,:,:)    = 1.0e-4_dp ! see ioinitial.f90

    CALL new_channel_object(status, modstr, 'tkem1'  &
         , reprid = GP_3D_MID, p3 = tkem1)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tkem1'       &
         , 'long_name', c='turbulent kinetic energy (t-1)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tkem1',     'units', c='m**2/s**2')
    CALL channel_halt(substr, status)
    tkem1(:,:,:)    = 1.0e-4_dp ! see ioinitial.f90

  END SUBROUTINE vertex_init_memory
!==============================================================================

!==============================================================================
  SUBROUTINE vertex_init_coupling                            ! Currently ECHAM/MESSy specific
    USE messy_main_channel,          ONLY: get_channel_object
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'vertex_init_coupling'
    INTEGER                     :: status

    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)
#ifdef ECHAM5

    CALL get_channel_object(status, 'g3b',   'ocu',    p2=ocu)
    CALL channel_halt(substr//' (ocu)',    status)

    CALL get_channel_object(status, 'g3b',   'ocv',    p2=ocv)
    CALL channel_halt(substr//' (ocv)',    status)

    CALL get_channel_object(status, 'g3b',   'az0',    p2=az0)
    CALL channel_halt(substr//' (az0)',    status)
#endif

    CALL get_channel_object(status, 'rad',   'albedo', p2=albedo)
    CALL channel_halt(substr//' (albedo)', status)

    CALL get_channel_object(status, 'rad',   'alsol',  p2=alsol)
    CALL channel_halt(substr//' (alsol)',  status)

    CALL get_channel_object(status, 'rad01', 'flxt',  p3=emter)
    CALL channel_halt(substr//' (emter)',  status)

    IF (ledith) THEN
       CALL get_channel_object(status, 'gwave', 'gweddy', p3=gweddy)
       CALL channel_halt(substr, status)
    END IF

    ! get leaf area index
    CALL get_channel_object(status, imp_lai%CHA, imp_lai%OBJ, p2=lai)
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)

  END SUBROUTINE vertex_init_coupling
!==============================================================================

!==============================================================================
  SUBROUTINE vertex_global_start

    USE messy_main_timer,         ONLY: lstart
#ifndef ECHAM5
    USE messy_main_data_bi,       ONLY: az0             ! Not defined in BLANK, MBM_CLAMS
#endif

    IMPLICIT NONE

    IF (lstart) THEN   ! see ioinitial.f90
       az0w(:,:)      = az0(:,:)
       az0i(:,:)      = az0(:,:)
       az0l(:,:)      = az0(:,:)
    END IF

  END SUBROUTINE vertex_global_start
!==============================================================================

!==============================================================================
  SUBROUTINE vertex_vdiff ( flag )
!
!**** *vertex_vdiff* - does the vertical exchange of u, v, t, q, xl, xi
!               and xt by turbulence.
!
!       This routine computes the physical tendencies of the seven
!   prognostic variables u, v, t, q, xl, xi and xt due to the vertical
!   exchange by turbulent (= non-moist convective) processes.
!   These tendencies are obtained as the difference between
!   the results after an over-implicit time-step, which starts from
!   values at t-1, and these t-1 values.
!   all the diagnostic computations (exchange coefficients, ...) are
!   done from the t-1 values. As a by-product the roughness length
!   over sea is updated accordingly to the *charnock formula. heat and
!   moisture surface fluxes and their derivatives against ts and ws,
!   later to be used for soil processes treatment, are also
!   computed as well as a stability value to be used as a diagnostic
!   of the depth of the well mixed layer in convective computations.
!
!**   Interface.
!     ----------
!
!          *vertex_vdiff* is called from *messy_vdiff* (messy_main_control).
!
!     Used data.
!      ----------
!
!  - 3d from mo_memory_g1a
!
!  pxtm1    : tracer variables (t-dt)
!
!  - 2d from mo_memory_g1a
!
!  pqm1     : humidity (t-dt)
!  ptm1     : temperature (t-dt)
!  pum1     : zonal wind (t-dt)
!  pvm1     : meridional wind (t-dt)
!  pxlm1    : cloud water (t-dt)
!  pxim1    : cloud ice (t-dt)
!  pxvar    : distribution width (b-a) (t-dt)

! - 2d from mo_memory_g3
!
!  paclc    : cloud cover
!
!  ptke     : turbulent kinetic energy at t+dt (unfiltered)
!  ptkem    :            "             at t    (unfiltered)
!  ptkem1   :            "             at t-dt   (filtered)
!
!  - 1d from mo_memory_g3
!
!  ptsl     : surface temperature over land
!  ptsw     :             "       over water
!  ptsi     :             "       over ice
!
!  pocu     : ocean u-velocity
!  pocv     : ocean v-velocity
!
!  pahfs    : surface sensible heat flux
!  pahfsl   :             "              over land
!  pahfsw   :             "              over water
!  pahfsi   :             "              over ice
!
!  pahfl    : surface latent heat flux
!  pahfll   :             "              over land
!  pahflw   :             "              over water
!  pahfli   :             "              over ice
!
!  pevap    : surface evaporation
!  pevapl   :             "        over land
!  pevapw   :             "        over water
!  pevapi   :             "        over ice
!
!  paz0     : roughness length for momentum
!  paz0l    :      "            over land
!  paz0w    :      "            over water
!  paz0i    :      "            over ice
!
!  pustr    : u-stress
!  pustrl   :     "     over land
!  pustrw   :     "     over sea
!  pustri   :     "     over ice
!
!  pvstr    : v-stress
!  pvstrl   :     "     over land
!  pvstrw   :     "     over water
!  pvstri   :     "     over ice
!
!  pdew2    : dew point temperature at 2 meter
!  peforest : forest coverage
!  psn      : snow depth
!  psnc     : snow depth on canopy
!  ptemp2   : temperature at 2 meter
!  ptsm1    : surface temperature (t-dt)
!  pwind10w : 10m wind over water
!  pu10     : u-wind at 10 meter
!  pv10     : v-wind at 10 meter
!  pwind10  : wind speed at 10 meter
!  pvdis    : boundary layer dissipation
!  pws      : surface soil wetness
!  pwsmx    : field capacity of soil
!  pwlmx    : skin reservoir
!  zlai     : leaf area index                            ! ju_te_20190212
!  pvgrat   : vegetation ratio
!
! - 2d within physics only
!
!  paphm1   : half level pressure (t-dt)
!  papm1    : full level pressure (t-dt)
!  ptvm1    : virtual temperature at t-dt
!  pvdiffp  : rate of change of qv due to vdiff routine for cover
!  pvmixtau : vdiff mixing timescale for variance and skewness
!
! - 1d within physics only
!
!  pgeom1   : geopotential above surface (t-dt)
!  psrfl    : net solar radiative flux at the surface
!  pqhfla   : moisture flux at the surface
!  pevapot  : potential evaporation
!  pcvs     : fractional snow cover (defined in *physc*)
!  pcvw     : wet skin fraction
!  loland   : land-sea flag
!
!        Tendencies
!
!  - 3d
!
!  pxtte    : tendencies of tracer variables
!
!  - 2d
!  pvol     : tendency of meridional wind
!  pvom     : tendency of zonal wind
!  pqte     : tendency of humidity
!  ptte     : tendency of temperature
!  pxlte    : tendency of cloud water
!  pxite    : tendency of cloud ice
!
!
!     Method.
!     -------
!
!        First an auxialiary variable cp(q)t+gz is created on which
!   the vertical diffusion process will work like on u,v and q. then
!   along the vertical and at the surface, exchange coefficients (with
!   the dimension of a pressure thickness) are computed for momentum
!   and for heat (sensible plus latent). the letters m and h are used
!   to distinguish them. the diffusioncoefficents depend on the
!   turbulent kinetic energy (tke) calculated by an additional
!   prognostic equation, which considers advection of tke.
!        In the second part of the routine the implicit linear
!   systems for u,v first and t,q second are solved by a *gaussian
!   elimination back-substitution method. for t and q the lower
!   boundary condition depends on the surface state.
!   for tke the lower boundary condition depends on the square of
!   the frictional velocity.
!   over land, two different regimes of evaporation prevail:
!   a stomatal resistance dependent one over the vegetated part
!   and a soil relative humidity dependent one over the
!   bare soil part of the grid mesh.
!   potential evaporation takes place over the sea, the snow
!   covered part and the liquid water covered part of the
!   grid mesh as well as in case of dew deposition.
!        Finally one returns to the variable temperature to compute
!   its tendency and the later is modified by the dissipation's effect
!   (one assumes no storage in the turbulent kinetic energy range) and
!   the effect of moisture diffusion on cp. z0 is updated and the
!   surface fluxes of t and q and their derivatives are prepared and
!   stored like the difference between the implicitely obtained
!   cp(q)t+gz and cp(q)t at the surface.
!
!
!     Reference.
!
!          See vertical diffusion's part of ECHAM's documentation
!     for details about the mathematics of this routine.
!
!     Authors.
!
!     u. schlese     dkrz-hamburg  feb-93
!       modified     e. roeckner  - 1994
!
!     j.-p. schulz   mpi - 1997 : implementation of implicit
!                                 coupling between land surface
!                                 and atmosphere.
!     m. esch, mpi, june 1999, echam5-modifications
!
!
!     based  on  original ecmwf version by j.f. geleyn  - 1982
!                              modified by c.b. blondin - 1986
!                                          h. feichter  - 1991
!                                          s. brinkop   - 1992
!                                          m. claussen  - 1993

#ifdef ECHAM5
  USE messy_main_timer,         ONLY: zdtime=>delta_time                     ! Time step
  USE messy_main_timer,         ONLY: ztmst=>time_step_len                   ! Usually twice time step; leap frog
  USE messy_main_grid_def_mem_bi, ONLY: kproma, kbdim=>nproma, krow=>jrow
  USE messy_main_grid_def_mem_bi, ONLY: klev=>nlev
#ifdef ECHAM5
  USE messy_main_grid_def_mem_bi, ONLY: klevm1=>nlevm1, klevp1=>nlevp1
#endif
  USE messy_main_tracer_mem_bi, ONLY: ktrac=>ntrac_gp
  USE messy_main_data_bi,       ONLY: wsmx, tsw, tslm1,cvs, cvw              ! Not defined in BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: ws                                     ! Not defined in VERTICO, BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: tsi, sn                                ! Only defined in ECHAM5 & CESM1
  USE messy_main_data_bi,       ONLY: paphm1=>aphm1,  papm1=>apm1            ! Only defined in ECHAM5 & CESM1
  USE messy_main_data_bi,       ONLY: geopot_3d                          ! Not defined in VERTICO, BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: loland_2d                              ! Only defined in ECHAM5 & CESM1
#ifndef MESSYTENDENCY
  USE messy_main_data_bi,       ONLY: vol_3d, vom_3d, qte_3d           ! Not defined in VERTICO, BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: xlte_3d, xite_3d                       ! Not defined in VERTICO, BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: tte_3d                                 ! Not defined in BLANK, MBM_CLAMS
#endif
  USE messy_main_grid_def_bi,   ONLY: philon_2d, philat_2d                   ! Not defined in BLANK, MBM_CLAMS
  USE messy_main_blather_bi,    ONLY: error_bi
  USE messy_main_timer,         ONLY: lstart
  USE messy_main_constants_mem, ONLY: vtmpc1, cpd=>cp_air, rd
  USE messy_main_constants_mem, ONLY: g, vtmpc2, tmelt, alv, als
  USE messy_main_tools,         ONLY: tlucua, jptlucu1, jptlucu2
  USE messy_main_data_bi,       ONLY: pxtems                                 ! Not defined in VERTICO, BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: cfml,   cfmw,  cfmi                    ! Not defined in BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: cfncl,  cfncw, cfnci                   ! Not defined in BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: cdnl,   cdnw,  cdni                    ! Not defined in BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: riw,    rii,   ril,  srfl              ! Not defined in BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: qsl,    qsw,   qsi,  phum              ! Only defined in ECHAM5 & COSMO
  USE messy_main_data_bi,       ONLY: chl,    cfhl,  cfhw, cfhi, cfm           ! Only defined in ECHAM5  ! op_mk_20180221 (added cfm)
  USE messy_main_data_bi,       ONLY: rho_surf                               ! Only defined in ECHAM5

  IMPLICIT NONE
  SAVE
#endif

  INTEGER, INTENT(IN) :: flag
#ifdef ECHAM5

#ifndef ECHAM5
  INTEGER  :: klevm1=klev-1, klevp1=klev+1
#endif

  ! Variables to be stored between the two separate vertex_vdiff calls
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pqm1, ptm1, pum1,   pvm1,    pxim1, pxlm1,  &
                                           zcfv, zcptgz, zqshear, zx

  INTEGER,  DIMENSION(:),   ALLOCATABLE :: ihpbl

  REAL(dp), DIMENSION(:),   ALLOCATABLE :: pfri,  pfrl,   pfrw,   zbmi,  zbmw,  zbml,  &
                                           zbhi,  zbhl,   zbhw,   zbni,  zbnl,  zbnw,  &
                                           zbhnl, zchi,   zchw,   zcpti, zcptw, zdqsl, &
                                           zdu2,  zdu2oc, zhsoil, zwet

  ! Pointers that are used in both separate vertex_vdiff calls
  REAL(dp), DIMENSION(:,:), POINTER     :: ptkem1  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER     :: zcfh    => NULL()
  REAL(dp), DIMENSION(:,:), POINTER     :: zxtems  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER     :: zcfm    => NULL()

  REAL(dp), DIMENSION(:),   POINTER     :: paz0hl  => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: paz0l   => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: pcvs    => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: pcvw    => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: pocu    => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: pocv    => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: psn     => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: ptsi    => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: ptslm1  => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: ptsw    => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: pws     => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: pwsmx   => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zcfmi   => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zcfml   => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zcfmw   => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zcfhi   => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zcfhl   => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zcfhw   => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zcfnci  => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zcfncl  => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zcfncw  => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zchl    => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zdens   => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zhum    => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zrii    => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zril    => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zriw    => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zsrfll  => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zqsi    => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zqsl    => NULL()
  REAL(dp), DIMENSION(:),   POINTER     :: zqsw    => NULL()

!******************************************************************************

  IF (flag .EQ. 1) CALL vertex_vdiff_1
  IF (flag .EQ. 2) CALL vertex_vdiff_2

  CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE vertex_vdiff_1

  USE messy_main_data_bi,       ONLY: qm1, tm1, xlm1, xim1                   ! Not defined in BLANK, MBM_CLAMS
#ifndef MESSYTENDENCY
  USE messy_main_data_bi,       ONLY: um1, vm1                               ! Not defined in BLANK, MBM_CLAMS
#endif
  USE messy_main_data_bi,       ONLY: srfl_2d                                ! Only defined in ECHAM5 & CESM1
  USE messy_main_data_bi,       ONLY: coriol_2d                              ! Not defined in VERTICO, BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: aclc                                   ! Not defined in BLANK, MBM_CLAMS
#ifndef ECHAM5
  USE messy_main_data_bi,       ONLY: az0                                    ! Not defined in BLANK, MBM_CLAMS
#endif
  USE messy_main_data_bi,       ONLY: landcov                                ! Only defined in ECHAM5
  USE messy_main_data_bi,       ONLY: seacov                                 ! Only correctly defined in ECHAM5, defined as "sea cover fraction (covered by ice)" in CESM1
  USE messy_main_data_bi,       ONLY: icecov                                 ! Only defined in ECHAM5 & CESM1
  USE messy_main_constants_mem, ONLY: api=>pi
  USE messy_main_data_bi,       ONLY: tvir, tvl,  tvw, tvi, fws              ! Not defined in BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: rco_leaf                               ! Not defined in BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: tpot_3d                                ! Not defined in VERTICO, BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: zust_2d                                ! Not defined in BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: zcdh_2d,    rinum_3d                   ! Only defined in ECHAM5
  USE messy_main_data_bi,       ONLY: cfh                                    ! Only defined in ECHAM5
  USE messy_main_data_bi,       ONLY: eps
  USE messy_main_data_bi,       ONLY: rh_2m

  IMPLICIT NONE

!******************************************************************************

  ! Temporary variables declared
  INTEGER :: it, it1, jk, jl

  LOGICAL :: lo, lookupoverflow

  REAL(dp):: zblend, zlam, zmix, zqdp, zdisc
  REAL(dp):: zchsnow,zchland,zchmean
  REAL(dp):: zfac,   zcons,  zktest, ztest,  zeps
  REAL(dp):: zdthv,  zepdu2, zepsr,  zepsec, zes
  !REAL(dp):: zln1,   zln2
  REAL(dp):: zplmax
  REAL(dp):: zqlwi1, zqlwi2, zqmitte
  REAL(dp):: zqsmit, zqst1,  zztvm
  REAL(dp):: zsdep2, zsdep1
  !REAL(dp):: zsoil,  zrsi, zsrfl,
  REAL(dp)::  zsrfld
  REAL(dp):: ztemitte, ztmit, zshear
  REAL(dp):: zusus1, zcor, zds, zdz, zzb, zdisl
  !REAL(dp):: zabcs
  REAL(dp):: zalo, zaloh
  REAL(dp):: zvirmitte
  !REAL(dp):: zwcrit, zwpwp
  REAL(dp):: zwslev, zwstop
  REAL(dp):: zepz0o, ztkemin, zepevap
  REAL(dp):: zcdn2m, zcdnr, zcfm2m, zchneu
  REAL(dp):: zucf, zust, zustarm, zustf, zusti, zustl, zustw
  REAL(dp):: zwstf, zz2geo, zfux, zfox, zalf
  REAL(dp):: zghabl, zsh, zsm
  REAL(dp):: zmult1, zmult2, zmult3, zmult4, zmult5
  REAL(dp):: zdus1,  zdus2,  zbuoy,  zdivv,  zdivv1
  REAL(dp):: zteldif, zqddif, zhexp, zconvs, zmonob, zstabf
  REAL(dp):: ztkesq, ztkevi, ztkevl, ztkevw
  REAL(dp):: zqtmit, zrdrv,  zrvrd,  zdqtot, zri

#ifndef MESSYTENDENCY
  REAL(dp):: pqmm1  (kbdim,klev),   ptmm1  (kbdim,klev),   pumm1   (kbdim,klev),   &
             pvmm1  (kbdim,klev),   pximm1 (kbdim,klev),   pxlmm1  (kbdim,klev)
#endif

  REAL(dp):: paclc  (kbdim,klev),   zcdum  (kbdim,klev),   zebsm   (kbdim,klev),   &
             zedif  (kbdim,klev),   zfaxe  (kbdim,klev),   zlteta1 (kbdim,klev),   &
             zqss   (kbdim,klev),   ztkevn (kbdim,klev),   ztvir1  (kbdim,klev),   &
             ptvm1  (kbdim,klev)

  REAL(dp):: zccover(kbdim,klevm1), zfaxen (kbdim,klevm1), zhh     (kbdim,klevm1), &
             zlwcmit(kbdim,klevm1), zqssm  (kbdim,klevm1), zqmit   (kbdim,klevm1), &
             ztemit (kbdim,klevm1), ztmitte(kbdim,klevm1), ztvirmit(kbdim,klevm1)

  INTEGER :: ihpblc (kbdim),        ihpbld (kbdim)

  REAL(dp):: palbedo(kbdim),        palsol (kbdim),        psrfl   (kbdim),        &
                                    zcfnchl(kbdim),        zcfnchw (kbdim),        &
             zchnl  (kbdim),        zchnw  (kbdim),        zcmi    (kbdim),        &
             zcml   (kbdim),        zcmw   (kbdim),        zcr     (kbdim),        &
             zhdyn  (kbdim),        zscfi  (kbdim),        zscfl   (kbdim),        &
             zscfw  (kbdim),        ztcoe  (kbdim),        ztesi   (kbdim),        &
             ztesl  (kbdim),        ztesw  (kbdim),        zucfhl  (kbdim),        &
             zucfi  (kbdim),        zucfl  (kbdim),        zucfw   (kbdim),        &
             zustari(kbdim),        zustarl(kbdim),        zustarw (kbdim),        &
             zwsti  (kbdim),        zwstl  (kbdim),        zwstw   (kbdim)

  REAL(dp):: pgwdiffco(kbdim,klev)

  ! Pointers to channel objects for use in first call to vertex_vdiff
  REAL(dp), POINTER, DIMENSION(:,:) :: ptke     => NULL()
  REAL(dp), POINTER, DIMENSION(:,:) :: ptkem    => NULL()
  REAL(dp), POINTER, DIMENSION(:,:) :: zteta1   => NULL()

  REAL(dp), POINTER, DIMENSION(:)   :: paz0     => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: paz0hi   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: paz0hw   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: paz0i    => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: paz0w    => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pustr    => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pustri   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pustrl   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pustrw   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pvstr    => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pvstri   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pvstrl   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pvstrw   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: zcdni    => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: zcdnl    => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: zcdnw    => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: ztvi     => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: ztvl     => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: ztvw     => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: zwet_tmp => NULL()

!******************************************************************************

  ! Allocation of variables to be stored between the two separate vertex_vdiff calls
  ALLOCATE(  pqm1   (kbdim,klev), ptm1 (kbdim,klev), pum1  (kbdim,klev),  &
             pvm1   (kbdim,klev), pxim1(kbdim,klev), pxlm1 (kbdim,klev)  )

  ALLOCATE(  zcfv (kbdim,klev), zcptgz(kbdim,klev),  &
             zqshear(kbdim,klev), zx   (kbdim,klev)                      )

  ALLOCATE(  pfri   (kbdim),      pfrl (kbdim),      pfrw  (kbdim)       )

  ALLOCATE(  ihpbl  (kbdim),      zbmi (kbdim),      zbmw  (kbdim),       &
             zbml   (kbdim),      zbhi (kbdim),      zbhl  (kbdim),       &
             zbhw   (kbdim),      zbni (kbdim),      zbnl  (kbdim),       &
             zbnw   (kbdim),      zbhnl(kbdim),      zchi  (kbdim),       &
             zchw   (kbdim),      zcpti(kbdim),      zcptw (kbdim),       &
             zdqsl  (kbdim),      zdu2 (kbdim),      zdu2oc(kbdim),       &
             zhsoil (kbdim),      zwet (kbdim)                           )

  ! Setting of pointers that are used in both separate vertex_vdiff calls
  ptkem1  => tkem1   (:,       :,krow)
  zcfh    => cfh     (1:kproma,:,krow)
  zxtems  => pxtems  (:,1,     :,krow)

  paz0hl  => az0hl   (:,         krow)
  paz0l   => az0l    (:,         krow)
  pcvs    => cvs     (:,         krow)
  pcvw    => cvw     (:,         krow)
  pocu    => ocu     (:,         krow)
  pocv    => ocv     (:,         krow)
  psn     => sn      (:,         krow)
  ptsi    => tsi     (:,         krow)
  ptslm1  => tslm1   (:,         krow)
  ptsw    => tsw     (:,         krow)
  pws     => ws      (:,         krow)
  pwsmx   => wsmx    (:,         krow)
  zcfmi   => cfmi    (:,         krow)
  zcfml   => cfml    (:,         krow)
  zcfmw   => cfmw    (:,         krow)
  zcfhi   => cfhi    (1:kproma,  krow)
  zcfhl   => cfhl    (1:kproma,  krow)
  zcfhw   => cfhw    (1:kproma,  krow)
  zcfnci  => cfnci   (:,         krow)
  zcfncl  => cfncl   (:,         krow)
  zcfncw  => cfncw   (:,         krow)
  zchl    => chl     (1:kproma,  krow)
  zdens   => rho_surf(:,         krow)
  zhum    => phum    (1:kproma,  krow)
  zrii    => rii     (:,         krow)
  zril    => ril     (:,         krow)
  zriw    => riw     (:,         krow)
  zsrfll  => srfl    (:,         krow)
  zqsi    => qsi     (1:kproma,  krow)
  zqsl    => qsl     (1:kproma,  krow)
  zqsw    => qsw     (1:kproma,  krow)

#ifndef MESSYTENDENCY
  ptte    => tte_3d   (1:kproma, :, krow)
  pqte    => qte_3d   (1:kproma, :, krow)
  pxlte   => xlte_3d  (1:kproma, :, krow)
  pxite   => xite_3d  (1:kproma, :, krow)
  pvol    => vol_3d   (1:kproma, :, krow)
  pvom    => vom_3d   (1:kproma, :, krow)
#endif
  pgeom1  => geopot_3d(1:kproma, :, krow)
!******************************************************************************

  ! Set (initial) values of arrays
  paclc   =  aclc    (:,:,     krow)
#ifndef MESSYTENDENCY
  pqmm1   =  qm1     (:,:,     krow)
  ptmm1   =  tm1     (:,:,     krow)
  pumm1   =  um1     (:,:,     krow)
  pvmm1   =  vm1     (:,:,     krow)
  pximm1  =  xim1    (:,:,     krow)
  pxlmm1  =  xlm1    (:,:,     krow)
#endif

  palbedo =  albedo  (:,       krow)
  palsol  =  alsol   (:,       krow)
  pfri    =  icecov  (:,       krow)
  pfrl    =  landcov (:,       krow)
  pfrw    =  seacov  (:,       krow)
  psrfl   =  srfl_2d (:,       krow)
  ptvm1   =  0.0_dp

  ! Setting of pointers that are used in first call to vertex_vdiff
  ptke    => tke     (:,:,     krow)
  ptkem   => tkem    (:,:,     krow)
  zteta1  => tpot_3d (:,:,     krow)

  paz0    => az0     (:,       krow)
  paz0hi  => az0hi   (:,       krow)
  paz0hw  => az0hw   (:,       krow)
  paz0i   => az0i    (:,       krow)
  paz0w   => az0w    (:,       krow)
  pustr   => ustr    (:,       krow)
  pustri  => ustri   (:,       krow)
  pustrl  => ustrl   (:,       krow)
  pustrw  => ustrw   (:,       krow)
  pvstr   => vstr    (:,       krow)
  pvstri  => vstri   (:,       krow)
  pvstrl  => vstrl   (:,       krow)
  pvstrw  => vstrw   (:,       krow)
  zcdni   => cdni    (:,       krow)
  zcdnl   => cdnl    (:,       krow)
  zcdnw   => cdnw    (:,       krow)
  ztvi    => tvi     (:,       krow)
  ztvl    => tvl     (:,       krow)
  ztvw    => tvw     (:,       krow)
  zwet_tmp => wet_tmp(1:kproma,krow)
  zcfm    => cfm     (1:kproma,:,krow)

!******************************************************************************

  lookupoverflow    = .FALSE.

  zxtems            = 0._dp

  ptvm1(1:kproma,:) = tm1(1:kproma,:,krow)*(1._dp+vtmpc1*qm1(1:kproma,:,krow) &
                      -(xlm1(1:kproma,:,krow)+xim1(1:kproma,:,krow)))

!*    PHYSICAL CONSTANTS.
!
  zrvrd   = vtmpc1+1._dp                    ! Rv/Rd = 1.608, since vtmpc1 = Rv/Rd-1 = 0.608
  zrdrv   = 1._dp/zrvrd                     ! Rd/Rv = 0.622

!*    SECURITY PARAMETERS.
!
  zepdu2  = 1.0_dp                          ! Minimum difference in |U|^2 between surface and lowest grid [m^2 s^-2]
  zepsr   = 1.e-10_dp                       ! Minimum scaled heat flux or radiation for which convective calculations are performed: 10^-10 m/s or W/m2 (surface TKE & stomatal resistance)
  zepevap = 1.e-10_dp                       ! Minimum conductivity of leaf stomata: 10^-10 m/s
  zepsec  = 1.e-2_dp                        ! Minimum value of C_M,H / (kappa * sqrt(C_M)) [-]
  zepz0o  = 2._dp                           ! Maximum roughness length for the calculation of u*: 2 m
  ztkemin = 1.e-10_dp                       ! Minimum TKE: 10^-10 m^2 s^-2

!*    COMPUTATIONAL CONSTANTS.
!
  zplmax  = 0.75_dp                         ! Fraction of field capacity above which the soil moisture stress correction function = 1
  zblend  = 100._dp                         ! Blending height to calculate z0 over partially snow-covered land
  zchneu  = .3_dp                           ! Factor used to estimate minimum zi using u*/fc

! Factor preceeding u*^2 as part of calculating TKE at the surface, E0
! changed in ECHAM from the documented 3.75 (Mailhot and Benoit, 1982) to
! 1/S_N_M^2, so that for neutral conditions E0 corresponds with
!  a K_M -> kappa z u* close to the surface
  zustf   = 1._dp/zsmn**2
  zwstf   = 0.2_dp                          ! Factor preceeding w*^2 as part of calculating TKE at the surface, E0

! The actual value of temperature, humidity, wind, cloud water and
! cloud ice are calculated and used all through the routine
! ATTENTION:  the code reads still as if the minus 1 values (t-dt)
! are used as in the original ECHAM5 code! However, with this addition
! it is consistent with the leapfrog operator splitting integration scheme

#ifndef MESSYTENDENCY
  do jk=1,klev
    ptm1(1:kproma,jk)  = ptmm1(1:kproma,jk)  + ptte(1:kproma,jk)  * ztmst
    pqm1(1:kproma,jk)  = pqmm1(1:kproma,jk)  + pqte(1:kproma,jk)  * ztmst
    pum1(1:kproma,jk)  = pumm1(1:kproma,jk)  + pvom(1:kproma,jk)  * ztmst
    pvm1(1:kproma,jk)  = pvmm1(1:kproma,jk)  + pvol(1:kproma,jk)  * ztmst
    pxlm1(1:kproma,jk) = pxlmm1(1:kproma,jk) + pxlte(1:kproma,jk) * ztmst
    pxim1(1:kproma,jk) = pximm1(1:kproma,jk) + pxite(1:kproma,jk) * ztmst
  enddo
#else
  ! get compute start values
  call mtend_get_start_l (mtend_id_t,  v0 = ptm1)
  call mtend_get_start_l (mtend_id_q,  v0 = pqm1)
  call mtend_get_start_l (mtend_id_xl, v0 = pxlm1)
  call mtend_get_start_l (mtend_id_xi, v0 = pxim1)
  call mtend_get_start_l (mtend_id_u,  v0 = pum1)
  call mtend_get_start_l (mtend_id_v,  v0 = pvm1)
#endif

!
!     ------------------------------------------------------------------
!
!*         2.     NEW THERMODYNAMIC VARIABLE AND BOUNDARY CONDITIONS.
!
!*         2.1     REPLACE T BY CP(Q)*T+GZ IN THE ATMOSPHERE.
!
  DO jk=1,klev
     DO jl=1,kproma
        zx(jl,jk)=pxlm1(jl,jk)+pxim1(jl,jk)                                  ! Total cloud water
        zcptgz(jl,jk)=pgeom1(jl,jk)+ptm1(jl,jk)                            & ! Dry static energy
                                   *cpd*(1._dp+vtmpc2*pqm1(jl,jk))
        zteta1(jl,jk)=ptm1(jl,jk)*(100000._dp/papm1(jl,jk))**(rd/cpd)        ! Potential temperature
        ztvir1(jl,jk)=zteta1(jl,jk)*(1._dp+vtmpc1*pqm1(jl,jk)-zx(jl,jk))     ! Virtual potential temperature
        lo=ptm1(jl,jk).GE.tmelt                                              ! Flag - true if ice is melting (T >= 0 degrees Celsius)
        zfaxe(jl,jk)=MERGE(alv,als,lo)                                       ! Factor for latent heat release: Lv when melting (vaporisation), Ls when freezing (sublimation)
        zusus1=(zfaxe(jl,jk)/cpd)*zteta1(jl,jk)/ptm1(jl,jk)*zx(jl,jk)        ! Intermediate factor - L/cp_d * theta/T * ql
        zlteta1(jl,jk)=zteta1(jl,jk)-zusus1                                  ! Approximation for liquid water potential temperature
        it = NINT(ptm1(jl,jk)*1000._dp)                                      ! Integer temperature in mK for lookup tables
        IF ((it<jptlucu1 .OR. it>jptlucu2) .AND. &
             (papm1(jl,jk) >= 1._dp)) lookupoverflow = .TRUE.            ! T should lie between 50 and 400 K
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/papm1(jl,jk)                                          ! Rd/Rv * es/p
        zes=MIN(zes,0.5_dp)                                                  ! Upper limit of 0.5 => upper limit to es of 80 % of total pressure
        zqss(jl,jk)=zes/(1._dp-vtmpc1*zes)                                   ! Saturation specific humidity
     END DO
  END DO

  tvir(1:kproma,krow) = ztvir1(1:kproma,klev)

  IF (lookupoverflow) THEN
     do jk=1,1,klev
        do jl=1,kproma
           if ( INT(ptm1(jl,jk)*1000.) <jptlucu1 .OR.                      &
                INT(ptm1(jl,jk)*1000.) >jptlucu2)                          &
                print*, ' vertex(1) ptm1jk: ',ptm1(jl,jk)                  &
                ,philon_2d(jl,krow),philat_2d(jl,krow)
        enddo
     enddo
     CALL error_bi('LOOKUP TABLE OVERFLOW - 1','vertex')
     lookupoverflow = .FALSE.
  ENDIF

  ! Interpolation to values at interface levels, based on relative pressure drop
  ! Note - interface level k is located between grid levels k & k+1; shifted from MESSy standard
  DO jk=1,klevm1
     DO jl=1,kproma
        zhh(jl,jk)=(pgeom1(jl,jk)-pgeom1(jl,jk+1))/g                         ! Vertical distance between adjacent grid levels (k & k+1)
        zsdep1=(paphm1(jl,jk)-paphm1(jl,jk+1))                             & ! Part of pressure drop over grid levels k & k+1 that occurs at level k
              /(paphm1(jl,jk)-paphm1(jl,jk+2))
        zsdep2=(paphm1(jl,jk+1)-paphm1(jl,jk+2))                           & ! Part of pressure drop over grid levels k & k+1 that occurs at level k+1
              /(paphm1(jl,jk)  -paphm1(jl,jk+2))
        zqssm(jl,jk)=zsdep1*zqss(jl,jk)+zsdep2*zqss(jl,jk+1)                 ! Saturation specific humidity, averaged between adjacent full levels
        ztmitte(jl,jk)=zsdep1*ptm1(jl,jk)+zsdep2*ptm1(jl,jk+1)               ! Absolute temperature, averaged between adjacent full levels
        ztvirmit(jl,jk)=zsdep1*ztvir1(jl,jk)+zsdep2*ztvir1(jl,jk+1)          ! Virtual potential temperature, averaged between adjacent full levels
        lo=ztmitte(jl,jk).GE.tmelt                                           ! Flag - true if ice is melting, averaged between adjacent full levels
        zfaxen(jl,jk)=MERGE(alv,als,lo)                                      ! L, averaged between adjacent full levels
        zlwcmit(jl,jk)=zsdep1*zx(jl,jk)+zsdep2*zx(jl,jk+1)                   ! Total cloud water, averaged between adjacent full levels
        zqmit(jl,jk)=zsdep1*pqm1(jl,jk)+zsdep2*pqm1(jl,jk+1)                 ! Specific humidity, averaged between adjacent full levels
        ztemit(jl,jk)=zsdep1*zteta1(jl,jk)+zsdep2*zteta1(jl,jk+1)            ! Potential temperature, averaged between adjacent full levels
        zccover(jl,jk)=paclc(jl,jk)*zsdep1+paclc(jl,jk+1)*zsdep2             ! Cloud cover, averaged between adjacent full levels
     END DO
  END DO

!
!*      2.2   surface humidity and virtual temperature
!                  for land, water and ice
!
  DO jl=1,kproma
!
!    land ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     it = NINT(ptslm1(jl)*1000._dp)                                          ! Integer land-surface temperature in mK for lookup tables
     IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
     it = MAX(MIN(it,jptlucu2),jptlucu1)
     zes=tlucua(it)/paphm1(jl,klevp1)                                        ! Rd/Rv * es/p at land surface
     zqsl(jl)=zes/(1._dp-vtmpc1*zes)                                         ! Saturation specific humidity at land surface
     it1=it+1
     it1= MAX(MIN(it1,jptlucu2),jptlucu1)                                    ! T + 1 mK; to determine dqs/dT at land surface
     zqst1=tlucua(it1)/paphm1(jl,klevp1)                                     ! Rd/Rv * es/p at land surface 1 mK warmer
     zqst1=zqst1/(1._dp-vtmpc1*zqst1)                                        ! Saturation specific humidity at land surface 1 mK warmer
     zdqsl(jl)=(zqst1-zqsl(jl))*1000._dp                                     ! dqs/dT at land surface
     pws(jl)=MIN(pws(jl),pwsmx(jl))                                          ! Surface soil wetness, limited by field capacity
     zwstop=MIN(0.1_dp,pwsmx(jl))                                            ! Maximum difference between surface soil wetness and field capacity where still part of land surface is moist
     zwslev=pwsmx(jl)-zwstop                                                 ! Minimum surface soil wetness where still part of land surface is moist
     IF(pws(jl).GT.zwslev.AND.pws(jl).GT.zplmin*pwsmx(jl)) THEN              ! Only part of the land surface is moist when the soil moisture content is higher than the level above and than wilting point
        zhum(jl)=0.5_dp*(1._dp-COS((pws(jl)-zwslev)*api/zwstop))             ! Moist land surface scales from 0 - 1 between minimum level and field capacity
     ELSE
        zhum(jl)=0._dp
     END IF
     zhsoil(jl)=pcvs(jl)+(1._dp-pcvs(jl))                                  & ! Saturated air at the surface where there is snow, or otherwise wet skin; otherwise relative humidity linked to high soil moisture content
                                  *(pcvw(jl)+(1._dp-pcvw(jl))*zhum(jl))
     lo=pqm1(jl,klev).GT.zqsl(jl)                                            ! When q in the lower air is higher than the saturation specific humidity at land, all land is covered with water
     zhsoil(jl)=MERGE(1._dp,zhsoil(jl),lo)
     ztesl(jl)=ptslm1(jl)*(1.e5_dp/paphm1(jl,klevp1))**(rd/cpd)              ! Potential temperature of land surface
     ztvl(jl)=ztesl(jl)*(1._dp+vtmpc1*zhsoil(jl)*zqsl(jl))                   ! Virtual potential temperature of land surface
!
!    water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     it = NINT(ptsw(jl)*1000._dp)                                            ! Integer water-surface temperature in mK for lookup tables
     IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
     it = MAX(MIN(it,jptlucu2),jptlucu1)
     zes=tlucua(it)/paphm1(jl,klevp1)                                        ! Rd/Rv * es/p at water surface
     zqsw(jl)=zes/(1._dp-vtmpc1*zes)                                         ! Saturation specific humidity at water surface
     zcptw(jl)=ptsw(jl)*cpd*(1._dp+vtmpc2*zqsw(jl))                          ! Dry static energy at water surface
     ztesw(jl)=ptsw(jl)*(1.e5_dp/paphm1(jl,klevp1))**(rd/cpd)                ! Potential temperature at water surface
     ztvw(jl)=ztesw(jl)*(1._dp+vtmpc1*zqsw(jl))                              ! Virtual potential temperature at water surface
!
!    ice   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     it = NINT(ptsi(jl)*1000._dp)                                            ! Integer ice-surface temperature in mK for lookup tables
     IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
     it = MAX(MIN(it,jptlucu2),jptlucu1)
     zes=tlucua(it)/paphm1(jl,klevp1)                                        ! Rd/Rv * es/p at ice surface
     zqsi(jl)=zes/(1._dp-vtmpc1*zes)                                         ! Saturation specific humidity at ice surface
     zcpti(jl)=ptsi(jl)*cpd*(1._dp+vtmpc2*zqsi(jl))                          ! Dry static energy at ice surface
     ztesi(jl)=ptsi(jl)*(1.e5_dp/paphm1(jl,klevp1))**(rd/cpd)                ! Potential temperature at ice surface
     ztvi(jl)=ztesi(jl)*(1._dp+vtmpc1*zqsi(jl))                              ! Virtual potential temperature at ice surface
!
  END DO

  IF (lookupoverflow) THEN
     do jl=1,kproma
        if ( INT(ptslm1(jl)*1000.) <jptlucu1 .OR.                          &
             INT(ptslm1(jl)*1000.) >jptlucu2)                              &
             print*, 'vertex(2) ptslm1: ',ptslm1(jl)                       &
                ,philon_2d(jl,krow),philat_2d(jl,krow)
        if ( INT(ptsw(jl)*1000.) <jptlucu1 .OR.                            &
             INT(ptsw(jl)*1000.) >jptlucu2)                                &
             print*, 'vertex(2) ptsw: ',ptsw(jl)                           &
                ,philon_2d(jl,krow),philat_2d(jl,krow)
        if ( INT(ptsi(jl)*1000.) <jptlucu1 .OR.                            &
             INT(ptsi(jl)*1000.) >jptlucu2)                                &
             print*, 'vertex(2) ptsi: ',ptsi(jl)                           &
                ,philon_2d(jl,krow),philat_2d(jl,krow)
     enddo
     CALL error_bi('LOOKUP TABLE OVERFLOW - 2','vertex')
     lookupoverflow = .FALSE.
  ENDIF

! surface net solar radiation flux over land
!
  DO jl=1,kproma
     zsrfld=psrfl(jl)/(1._dp-palbedo(jl))                                    ! Total incoming solar radiative flux at the surface
     zsrfll(jl)=zsrfld*(1._dp-palsol(jl))                                    ! Net solar radiative flux over land
  END DO
!
!*         2.3     COMPUTATION OF THE STOMATAL RESISTANCE
!
! ju_te_20180625+
!!$  DO jl=1,kproma                                                             ! Calculating stomatal resistance, r_stom, according to Eq. (8) of Ganzeveld and Lelieveld (1995)
!!$     zwcrit=zplmax*pwsmx(jl)                                                 ! Soil moisture content above which the soil moisture stress correction function = 1
!!$     zwpwp=zplmin*pwsmx(jl)                                                  ! Wilting point
!!$     zsoil=MAX(0._dp,MIN(1._dp,(pws(jl)-zwpwp)/(zwcrit-zwpwp)))              ! Soil moisture stress correction function
!!$     zsrfl=MAX(zepsr,zsrfll(jl)*cvrad)                                       ! Net PAR over land
!!$     zabcs=(cva+cvb*cvc)/(cvc*zsrfl)                                         ! Coefficient d * PAR
!!$     zlai=pvlt(jl)                                                           ! LAI
!!$     zln1=LOG((zabcs*EXP(cvk*zlai)+1._dp)/(zabcs+1._dp))                     ! First logarithm in Eq. (8)
!!$     zln2=LOG((zabcs+EXP(-cvk*zlai))/(zabcs+1._dp))                          ! Second logarithm in Eq. (8)
!!$     zrsi=(cvb*zln1/cvabc-zln2)/(cvk*cvc)                                    ! 1 / (r_stom * F(ws) )
!!$     zwet(jl)=1._dp/(zrsi*zsoil+zepevap)                                     ! r_stom

  ! calculate vegetation resistance
  CALL vertex_calc_rstom(kproma, zwet, fws(:,krow), rco_leaf(:,krow) &
       , pws, pwsmx  &
       , zsrfll, lai(1:kproma,krow),  ptslm1(:), rh_2m(:,krow))

  DO jl=1, kproma
! ju_te_20180625-
! Check for dew (0 resistance) and store r_stom
     lo=pqm1(jl,klev).GT.zqsl(jl)                                            ! Check if air is saturated over land
     zwet(jl)=MERGE(0._dp,zwet(jl),lo)                                       ! If air is saturated, r_stom = 0 (since liquid water)
! ju_te_20180625-
!!$     ! mz_lg_20020129+
!!$     fws(jl,krow) = zsoil                                                    ! F(ws) - to be passed on to DDEP
!!$     ! calculation of leaf stomatal resistance using the echam
!!$     ! equation applying an LAI of 1, included by Laurens Ganzeveld, 18/10/01
!!$     rco_leaf(jl,krow)=cvk*cvc/(cvb*LOG((zabcs*EXP(cvk)+1._dp)             & ! r_stom * F(ws) - to be passed on to DDEP
!!$           /(zabcs+1._dp))/cvabc-LOG((zabcs+EXP(-cvk))/(zabcs+1._dp)))
!!$     ! mz_lg_20020129-
! ju_te_20180625-
     zwet_tmp(jl)=zwet(jl) !store temporary result for H2OISO op_re_20140904
  END DO

!     ------------------------------------------------------------------
!
!*         3.     COMPUTATION OF THE EXCHANGE COEFFICIENTS.
!
!        THE SURFACE LAYER IS NOW COMPUTED BEFORE THE OTHER LEVELS
!
!        3.1       COMPUTATION OF BASIC QUANTITIES: WIND SHEAR,
!                  RICHARDSON NUMBER,SQUARED MIXING LENGTHS, UNSTABLE
!                  AND STABLE CASE COMMON FACTORS AND NEUTRAL CASE
!                  COMMON PART OF THE DRAG COEFFICIENTS.
!
  DO jl=1,kproma
!
!    land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zdu2(jl)=MAX(zepdu2,pum1(jl,klev)**2+pvm1(jl,klev)**2)                  ! Square of Delta U (between ground (0 m/s) and center lowest grid cell)
     zqmitte=(pqm1(jl,klev)+zqsl(jl)*zhsoil(jl))/2._dp                       ! Average of water vapor between surface and center lowest grid cell
     zqtmit=zx(jl,klev)*0.5_dp+zqmitte                                       ! Average of total atmospheric water content between surface and center lowest grid cell
     ztmit=(ptm1(jl,klev)+ptslm1(jl))/2._dp                                  ! Average temperature between surface and center lowest grid cell
     zqsmit=(zqss(jl,klev)+zqsl(jl))/2._dp                                   ! Average saturation specific humidity between surface and center lowest grid cell
     ztemitte=(zteta1(jl,klev)+ztesl(jl))/2._dp                              ! Average potential temperature between surface and center lowest grid cell
     zvirmitte=(ztvir1(jl,klev)+ztvl(jl))/2._dp                              ! Average virtual potential temperature between surface and center lowest grid cell
     zteldif=zlteta1(jl,klev)-ztesl(jl)                                      ! Difference in liquid potential temperature between surface and center lowest grid cell
     zqlwi1=pqm1(jl,klev)+zx(jl,klev)                                        ! Total atmospheric water content in lowest grid cell
     zqlwi2=zqsl(jl)*zhsoil(jl)                                              ! Atmospheric water content (vapor) at surface
     zqddif=zqlwi1-zqlwi2                                                    ! Difference in total atmospheric water content between surface and center lowest grid cell
     zfux=zfaxe(jl,klev)/(cpd*ztmit)                                         ! L / (cp_d T) between surface and center lowest grid cell
     zfox=zfaxe(jl,klev)/(rd*ztmit)                                          ! L / (Rd T) between surface and center lowest grid cell
     zmult1=1._dp+vtmpc1*zqtmit                                              ! Calculating intermediates for th_v' for sat. air: 1 + (Rv/Rd - 1) q_t between surface and center lowest grid cell
     zmult2=zfux*zmult1-zrvrd                                                ! Calculating intermediates for th_v' for sat. air: L (1 + (Rv/Rd - 1) q_t) / (cp_d T) - Rv/Rd between surface and center lowest grid cell
     zmult3=zrdrv*zfox*zqsmit/(1._dp+zrdrv*zfox*zfux*zqsmit)                 ! Calculating intermediates for th_v' for sat. air: L q_s /(Rv T (1 + L / (Rv T) L / (cp_d T) q_s ) )
     zmult5=zmult1-zmult2*zmult3                                             ! c1 for sat. air in the equation th_v' = c1 th_l' + c_2 th q'
     zmult4=zfux*zmult5-1._dp                                                ! c2 for sat. air in equation above: c2 = c1 L / (cp_d T) - 1
     zdus1=paclc(jl,klev)*zmult5+(1._dp-paclc(jl,klev))*zmult1               ! c1 for entire land area
     zdus2=paclc(jl,klev)*zmult4+(1._dp-paclc(jl,klev))*vtmpc1               ! c2 for entire land area
     zbuoy=zdus1*zteldif+zdus2*ztemitte*zqddif                               ! Difference in virtual potential temperature between surface and center lowest grid cell
     zril(jl)=pgeom1(jl,klev)*zbuoy/(zvirmitte*zdu2(jl))                     ! Bulk richardson number between land surface and center lowest grid cell
     paz0hl(jl)=MIN(1._dp,paz0l(jl))                                         ! Roughness length for heat over land; upper limit of 1 m
     IF(pcvs(jl).GT.0._dp) THEN                                              ! Blending of z0 for heat between that over land and that over ice - according to Claussen (1991), Eq. 2.2.5
       zchsnow=(LOG(zblend/paz0i(jl)))**2
       zchland=(LOG(zblend/paz0hl(jl)))**2
       zchsnow=pcvs(jl)/zchsnow
       zchland=(1._dp-pcvs(jl))/zchland
       zchmean=1._dp/SQRT(zchsnow+zchland)
       paz0hl(jl)=zblend*EXP(-zchmean)
     END IF
     zalo     =LOG(1._dp+pgeom1(jl,klev)/(g*paz0l(jl)))                      ! Measure of height relative to z0m, used for C_N_M and C_N_H
     zaloh    =LOG(1._dp+pgeom1(jl,klev)/(g*paz0hl(jl)))                     ! Measure of height relative to z0h, used for C_N_H
     zcdnl(jl)=(ckap/zalo)**2                                                ! C_N_M - neutral transfer coefficient for momentum
     zchnl(jl)=ckap**2/(zalo*zaloh)                                          ! C_N_H - neutral transfer coefficient for heat
     zdens(jl)=paphm1(jl,klevp1)/rd/                                       & ! Air density at the surface
               (ptm1(jl,klev)*(1._dp+vtmpc1*pqm1(jl,klev)-zx(jl,klev)))
     zcons=cvdifts*ztmst*g*zdens(jl)                                         ! Computational factor 1.5 * leap frog time step * g * air density
     zcfncl(jl)=zcons*SQRT(zdu2(jl))*zcdnl(jl)                               ! Computational factor 1.5 * leap frog time step * g * air density * C_N_M * |U|
     zcfnchl(jl)=zcons*SQRT(zdu2(jl))*zchnl(jl)                              ! Computational factor 1.5 * leap frog time step * g * air density * C_N_H * |U|
     zdthv=MAX(0._dp,(ztvl(jl)-ztvir1(jl,klev)))                             ! Difference in virtual potential temperature between surface and center lowest grid cell
     zwstl(jl)=zdthv*SQRT(zdu2(jl))/zvirmitte                                ! - Delta{th_v} |U| / th_v = w'th_v'/ C_H / th_v
!
!    water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!    correction for water and ice points
!
     zdu2oc(jl)=MAX(zepdu2,(pum1(jl,klev)-pocu(jl))**2                     & ! Square of Delta U (between surface and center lowest grid cell) over ocean
                    +(pvm1(jl,klev)-pocv(jl))**2)
!
     zqmitte=(pqm1(jl,klev)+zqsw(jl))/2._dp                                  ! Average of water vapor between surface and center lowest grid cell
     zqtmit=zx(jl,klev)*0.5_dp+zqmitte                                       ! Average of total atmospheric water content between surface and center lowest grid cell
     ztmit=(ptm1(jl,klev)+ptsw(jl))/2._dp                                    ! Average temperature between surface and center lowest grid cell
     zqsmit=(zqss(jl,klev)+zqsw(jl))/2._dp                                   ! Average saturation specific humidity between surface and center lowest grid cell
     ztemitte=(zteta1(jl,klev)+ztesw(jl))/2._dp                              ! Average potential temperature between surface and center lowest grid cell
     zvirmitte=(ztvir1(jl,klev)+ztvw(jl))/2._dp                              ! Average virtual potential temperature between surface and center lowest grid cell
     zqlwi1=pqm1(jl,klev)+zx(jl,klev)                                        ! Total atmospheric water content in lowest grid cell
     zqlwi2=zqsw(jl)                                                         ! Atmospheric water content (vapor) at surface
     zqddif=zqlwi1-zqlwi2                                                    ! Difference in total atmospheric water content between surface and center lowest grid cell
     zfux=zfaxe(jl,klev)/(cpd*ztmit)                                         ! L / (cp_d T) between surface and center lowest grid cell
     zfox=zfaxe(jl,klev)/(rd*ztmit)                                          ! L / (Rd T) between surface and center lowest grid cell
     zmult1=1._dp+vtmpc1*zqtmit                                              ! Calculating intermediates for th_v' for sat. air: 1 + (Rv/Rd - 1) q_t between surface and center lowest grid cell
     zmult2=zfux*zmult1-zrvrd                                                ! Calculating intermediates for th_v' for sat. air: L (1 + (Rv/Rd - 1) q_t) / (cp_d T) - Rv/Rd between surface and center lowest grid cell
     zmult3=zrdrv*zfox*zqsmit/(1._dp+zrdrv*zfox*zfux*zqsmit)                 ! Calculating intermediates for th_v' for sat. air: L q_s /(Rv T (1 + L / (Rv T) L / (cp_d T) q_s ) )
     zmult5=zmult1-zmult2*zmult3                                             ! c1 for sat. air in the equation th_v' = c1 th_l' + c_2 th q'
     zmult4=zfux*zmult5-1._dp                                                ! c2 for sat. air in equation above: c2 = c1 L / (cp_d T) - 1
     zdus1=paclc(jl,klev)*zmult5+(1._dp-paclc(jl,klev))*zmult1               ! c1 for entire ocean area
     zdus2=paclc(jl,klev)*zmult4+(1._dp-paclc(jl,klev))*vtmpc1               ! c2 for entire ocean area
     zteldif=zlteta1(jl,klev)-ztesw(jl)                                      ! Difference in liquid potential temperature between surface and center lowest grid cell
     zbuoy=zdus1*zteldif+zdus2*ztemitte*zqddif                               ! Difference in virtual potential temperature between sea surface and center lowest grid cell
     zriw(jl)=pgeom1(jl,klev)*zbuoy/(zvirmitte*zdu2oc(jl))                   ! Bulk richardson number between sea surface and center lowest grid cell
     paz0hw(jl)=paz0w(jl)*EXP(2._dp-86.276_dp*paz0w(jl)**0.375_dp)           ! Roughness length for heat over sea
     zalo=LOG(1._dp+pgeom1(jl,klev)/(g*paz0w(jl)))                           ! Measure of height relative to z0m, used for C_N_M and C_N_H
     zaloh=LOG(1._dp+pgeom1(jl,klev)/(g*paz0hw(jl)))                         ! Measure of height relative to z0h, used for C_N_H
     zcdnw(jl)=(ckap/zalo)**2                                                ! C_N_M - neutral transfer coefficient for momentum
     zchnw(jl)=ckap**2/(zalo*zaloh)                                          ! C_N_H - neutral transfer coefficient for heat
     zcfncw(jl)=zcons*SQRT(zdu2oc(jl))*zcdnw(jl)                             ! Computational factor 1.5 * leap frog time step * g * air density * C_N_M * |U|
     zcfnchw(jl)=zcons*SQRT(zdu2oc(jl))*zchnw(jl)                            ! Computational factor 1.5 * leap frog time step * g * air density * C_N_H * |U|
     zdthv=MAX(0._dp,(ztvw(jl)-ztvir1(jl,klev)))                             ! Difference in virtual potential temperature between surface and center lowest grid cell
     zwstw(jl)=zdthv*SQRT(zdu2oc(jl))/zvirmitte                              ! - Delta{th_v} |U| / th_v = - w'th_v'/ C_v / th_v   over sea
     zcr(jl)=(cfreec/(zchnw(jl)*SQRT(zdu2oc(jl))))*ABS(zbuoy)**(1._dp/3._dp) ! C_R, used to calculate C_H for unstable conditions
!
!    ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zqmitte=(pqm1(jl,klev)+zqsi(jl))/2._dp                                  ! Average of water vapor between surface and center lowest grid cell
     zqtmit=zx(jl,klev)*0.5_dp+zqmitte                                       ! Average of total atmospheric water content between surface and center lowest grid cell
     ztmit=(ptm1(jl,klev)+ptsi(jl))/2._dp                                    ! Average temperature between surface and center lowest grid cell
     zqsmit=(zqss(jl,klev)+zqsi(jl))/2._dp                                   ! Average saturation specific humidity between surface and center lowest grid cell
     ztemitte=(zteta1(jl,klev)+ztesi(jl))/2._dp                              ! Average potential temperature between surface and center lowest grid cell
     zvirmitte=(ztvir1(jl,klev)+ztvi(jl))/2._dp                              ! Average virtual potential temperature between surface and center lowest grid cell
     zqlwi1=pqm1(jl,klev)+zx(jl,klev)                                        ! Total atmospheric water content in lowest grid cell
     zqlwi2=zqsi(jl)                                                         ! Atmospheric water content (vapor) at surface
     zqddif=zqlwi1-zqlwi2                                                    ! Difference in total atmospheric water content between surface and center lowest grid cell
     zfux=zfaxe(jl,klev)/(cpd*ztmit)                                         ! L / (cp_d T) between surface and center lowest grid cell
     zfox=zfaxe(jl,klev)/(rd*ztmit)                                          ! L / (Rd T) between surface and center lowest grid cell
     zmult1=1._dp+vtmpc1*zqtmit                                              ! Calculating intermediates for th_v' for sat. air: 1 + (Rv/Rd - 1) q_t between surface and center lowest grid cell
     zmult2=zfux*zmult1-zrvrd                                                ! Calculating intermediates for th_v' for sat. air: L (1 + (Rv/Rd - 1) q_t) / (cp_d T) - Rv/Rd between surface and center lowest grid cell
     zmult3=zrdrv*zfox*zqsmit/(1._dp+zrdrv*zfox*zfux*zqsmit)                 ! Calculating intermediates for th_v' for sat. air: L q_s /(Rv T (1 + L / (Rv T) L / (cp_d T) q_s ) )
     zmult5=zmult1-zmult2*zmult3                                             ! c1 for sat. air in the equation th_v' = c1 th_l' + c_2 th q'
     zmult4=zfux*zmult5-1._dp                                                ! c2 for sat. air in equation above: c2 = c1 L / (cp_d T) - 1
     zdus1=paclc(jl,klev)*zmult5+(1._dp-paclc(jl,klev))*zmult1               ! c1 for entire ice/snow area
     zdus2=paclc(jl,klev)*zmult4+(1._dp-paclc(jl,klev))*vtmpc1               ! c2 for entire ice/snow area
     zteldif=zlteta1(jl,klev)-ztesi(jl)                                      ! Difference in liquid potential temperature between surface and center lowest grid cell
     zbuoy=zdus1*zteldif+zdus2*ztemitte*zqddif                               ! Difference in virtual potential temperature between sea surface and center lowest grid cell
     zrii(jl)=pgeom1(jl,klev)*zbuoy/(zvirmitte*zdu2oc(jl))                   ! Bulk richardson number between snow/ice surface and center lowest grid cell
     paz0hi(jl)=paz0i(jl)                                                    ! z0h = z0m = 1 mm over ice
     zalo=LOG(1._dp+pgeom1(jl,klev)/(g*paz0i(jl)))                           ! Measure of height relative to z0m (= z0h over snow/ice), used for C_N
     zcdni(jl)=(ckap/zalo)**2                                                ! C_N - neutral transfer coefficient for momentum and heat (same over snow/ice)
     zcfnci(jl)=zcons*SQRT(zdu2oc(jl))*zcdni(jl)                             ! Computational factor 1.5 * leap frog time step * g * air density * C_N * |U|
     zdthv=MAX(0._dp,(ztvi(jl)-ztvir1(jl,klev)))                             ! Difference in virtual potential temperature between surface and center lowest grid cell
     zwsti(jl)=zdthv*SQRT(zdu2oc(jl))/zvirmitte                              ! - Delta{th_v} |U| / th_v = - w'th_v'/ C_v / th_v   over snow/ice
  END DO
!
!     3.2  DIMENSIONLESS HEAT TRANSFER COEFFICIENTS MULTIPLIED
!          BY PRESSURE THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE
!
  DO jl=1,kproma
!
!    land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     IF(zril(jl).GT.0._dp) THEN                                              ! Stability functions for stable conditions
        zscfl(jl) =SQRT(1._dp+ABS(zril(jl)))                                 ! Part of stability functions f_M,H for stable conditions
        zcfml(jl) = zcfncl(jl)/(1._dp+2._dp*cc*zril(jl)/zscfl(jl))           ! Computational factor 1.5 * leap frog time step * g * air density * C_M * |U|
        zcfhl(jl) =zcfnchl(jl)/(1._dp+2._dp*cc*zril(jl)*zscfl(jl))           ! Computational factor 1.5 * leap frog time step * g * air density * C_H * |U|
        zcml(jl)  =  zcdnl(jl)/(1._dp+2._dp*cc*zril(jl)/zscfl(jl))           ! C_M - transfer coefficient for momentum
        zchl(jl)  =  zchnl(jl)/(1._dp+2._dp*cc*zril(jl)*zscfl(jl))           ! C_H - transfer coefficient for heat
     ELSE                                                                    ! Stability functions for unstable conditions
        zucfl(jl) =1._dp/(1._dp+3._dp*(cc**2)*zcdnl(jl)*SQRT(              & ! Part of stability function f_M for unstable conditions over land
                  ABS(zril(jl))*(1._dp +pgeom1(jl,klev)/(g*paz0l(jl)))))
        zucfhl(jl)=1._dp/(1._dp+3._dp*(cc**2)*zchnl(jl)*SQRT(              & ! Part of stability function f_H for unstable conditions over land
                  ABS(zril(jl))*(1._dp +pgeom1(jl,klev)/(g*paz0hl(jl)))))
        zcfml(jl) = zcfncl(jl)*(1._dp-2._dp*cc*zril(jl)*zucfl(jl))           ! Computational factor 1.5 * leap frog time step * g * air density * C_M * |U|
        zcfhl(jl) =zcfnchl(jl)*(1._dp-3._dp*cc*zril(jl)*zucfhl(jl))          ! Computational factor 1.5 * leap frog time step * g * air density * C_H * |U|
        zcml(jl)  =  zcdnl(jl)*(1._dp-2._dp*cc*zril(jl)*zucfl(jl))           ! C_M - transfer coefficient for momentum
        zchl(jl)  =  zchnl(jl)*(1._dp-3._dp*cc*zril(jl)*zucfhl(jl))          ! C_H - transfer coefficient for heat
     END IF
!
!    water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     IF(zriw(jl).GT.0._dp) THEN                                              ! Stability functions for stable conditions
        zscfw(jl)=SQRT(1._dp+ABS(zriw(jl)))                                  ! Part of stability functions f_M,H for stable conditions over ocean
        zcfmw(jl)= zcfncw(jl)/(1._dp+2._dp*cc*zriw(jl)/zscfw(jl))            ! Computational factor 1.5 * leap frog time step * g * air density * C_M * |U|
        zcfhw(jl)=zcfnchw(jl)/(1._dp+2._dp*cc*zriw(jl)*zscfw(jl))            ! Computational factor 1.5 * leap frog time step * g * air density * C_H * |U|
        zcmw(jl) =  zcdnw(jl)/(1._dp+2._dp*cc*zriw(jl)/zscfw(jl))            ! C_M - transfer coefficient for momentum
        zchw(jl) =  zchnw(jl)/(1._dp+2._dp*cc*zriw(jl)*zscfw(jl))            ! C_H - transfer coefficient for heat
     ELSE                                                                    ! Stability functions for unstable conditions
        zucfw(jl)=1._dp/(1._dp+3._dp*(cc**2)*zcdnw(jl)*SQRT(               & ! Part of stability function f_M for unstable conditions over ocean
                 ABS(zriw(jl))*(1._dp+pgeom1(jl,klev)/(g*paz0w(jl)))))
        zcfmw(jl)= zcfncw(jl)*(1._dp-2._dp*cc*zriw(jl)*zucfw(jl))            ! Computational factor 1.5 * leap frog time step * g * air density * C_M * |U|
        zcfhw(jl)=zcfnchw(jl)*(1._dp+zcr(jl)**cgam)**(1._dp/cgam)            ! Computational factor 1.5 * leap frog time step * g * air density * C_H * |U|
        zcmw(jl) =  zcdnw(jl)*(1._dp-2._dp*cc*zriw(jl)*zucfw(jl))            ! C_M - transfer coefficient for momentum
        zchw(jl) =  zchnw(jl)*(1._dp+zcr(jl)**cgam)**(1._dp/cgam)            ! C_H - transfer coefficient for heat
     END IF
!
!    ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     IF(zrii(jl).GT.0._dp) THEN                                              ! Stability functions for stable conditions
        zscfi(jl)=SQRT(1._dp+ABS(zrii(jl)))                                  ! Part of stability functions f_M,H for stable conditions over ocean
        zcfmi(jl)=zcfnci(jl)/(1._dp+2._dp*cc*zrii(jl)/zscfi(jl))             ! Computational factor 1.5 * leap frog time step * g * air density * C_M * |U|
        zcfhi(jl)=zcfnci(jl)/(1._dp+2._dp*cc*zrii(jl)*zscfi(jl))             ! Computational factor 1.5 * leap frog time step * g * air density * C_H * |U|
        zcmi(jl) = zcdni(jl)/(1._dp+2._dp*cc*zrii(jl)/zscfi(jl))             ! C_M - transfer coefficient for momentum
        zchi(jl) = zcdni(jl)/(1._dp+2._dp*cc*zrii(jl)*zscfi(jl))             ! C_H - transfer coefficient for heat
     ELSE                                                                    ! Stability functions for unstable conditions
        zucfi(jl)=1._dp/(1._dp+3._dp*(cc**2)*zcdni(jl)*SQRT(               & ! Part of stability functions f_M,H for unstable conditions over snow/ice
                 ABS(zrii(jl))*(1._dp+pgeom1(jl,klev)/(g*paz0i(jl)))))
        zcfmi(jl)=zcfnci(jl)*(1._dp-2._dp*cc*zrii(jl)*zucfi(jl))             ! Computational factor 1.5 * leap frog time step * g * air density * C_M * |U|
        zcfhi(jl)=zcfnci(jl)*(1._dp-3._dp*cc*zrii(jl)*zucfi(jl))             ! Computational factor 1.5 * leap frog time step * g * air density * C_H * |U|
        zcmi(jl) = zcdni(jl)*(1._dp-2._dp*cc*zrii(jl)*zucfi(jl))             ! C_M - transfer coefficient for momentum
        zchi(jl) = zcdni(jl)*(1._dp-3._dp*cc*zrii(jl)*zucfi(jl))             ! C_H - transfer coefficient for heat
     END IF
!
!    aggregated exchange coefficient for momentum
!
!    Factors A to be used for the Richtmeyer and Morton scheme (Schulz et al., 2001: Appendix)
     zcfm(jl,klev)=pfrl(jl)*zcfml(jl)+pfrw(jl)*zcfmw(jl)+pfri(jl)*zcfmi(jl)  ! Computational factor 1.5 * leap frog time step * g * air density * C_M * |U|
     zcdum(jl,klev)=zcfm(jl,klev)                                            ! Computational factor 1.5 * leap frog time step * g * air density * C_M * |U|
!
!    interpolation functions for diagnostics
!
!    land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zbnl(jl)=ckap/SQRT(zcdnl(jl))                                           ! ln(1 + h/z0m) = kappa / sqrt(C_N_M)
     zbhnl(jl)=ckap/SQRT(zchnl(jl))                                          ! sqrt( ln(1 + h/z0h) ln(1 + h/z0m) ) = kappa / sqrt(C_N_H)
     zbml(jl)=MAX(zepsec,SQRT(zcml(jl))/ckap)                                ! sqrt(C_M) / kappa
     zbhl(jl)=MAX(zepsec,zchl(jl)/zbml(jl)/ckap**2)                          ! C_H / ( sqrt(C_M) kappa )
     zbml(jl)=1._dp/zbml(jl)                                                 ! kappa / sqrt(C_M)
     zbhl(jl)=1._dp/zbhl(jl)                                                 ! sqrt(C_M) kappa / C_H
!
!    water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zbnw(jl)=ckap/SQRT(zcdnw(jl))                                           ! ln(1 + h/z0m) = kappa / sqrt(C_N_M)
     zbmw(jl)=MAX(zepsec,SQRT(zcmw(jl))/ckap)                                ! sqrt(C_M) / kappa
     zbhw(jl)=MAX(zepsec,zchw(jl)/zbmw(jl)/ckap**2)                          ! C_H / ( sqrt(C_M) kappa )
     zbmw(jl)=1._dp/zbmw(jl)                                                 ! kappa / sqrt(C_M)
     zbhw(jl)=1._dp/zbhw(jl)                                                 ! sqrt(C_M) kappa / C_H
!
!    ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zbni(jl)=ckap/SQRT(zcdni(jl))                                           ! ln(1 + h/z0) = kappa / sqrt(C_N)
     zbmi(jl)=MAX(zepsec,SQRT(zcmi(jl))/ckap)                                ! sqrt(C_M) / kappa
     zbhi(jl)=MAX(zepsec,zchi(jl)/zbmi(jl)/ckap**2)                          ! C_H / ( sqrt(C_M) kappa )
     zbmi(jl)=1._dp/zbmi(jl)                                                 ! kappa / sqrt(C_M)
     zbhi(jl)=1._dp/zbhi(jl)                                                 ! sqrt(C_M) kappa / C_H
!
     zcdh_2d(jl,krow)=pfrl(jl)*zchl(jl)+pfrw(jl)*zchw(jl)+pfri(jl)*zchi(jl)  ! C_H
  END DO
!
!*       3.4       COMPUTATION OF SURFACE SHEAR AND TKE, AND THE PBL EXTENSION.
!
  DO jl=1,kproma
!
!    land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     lo=paz0l(jl).GT.zepz0o
     zcdn2m=MERGE((ckap/LOG(1._dp+pgeom1(jl,klev)/(g*zepz0o)))**2,         & ! C_N_M(2m) - Equal to C_N_M, but adapted if z0m > 2m
                     zcdnl(jl),lo)
     zcdnr=zcdn2m/zcdnl(jl)                                                  ! C_N_M(2m)/C_N_M
     zcfm2m=MERGE(zcdn2m*SQRT(zdu2(jl))*(1._dp-2._dp*cc*zril(jl)           & ! C_M * (C_N_M(2m)/C_N_M) * |U| if Ri>= 0 or z0m <= 2 m
                  /(1._dp+3._dp*(cc**2)*zcdn2m*SQRT(ABS(zril(jl))          & ! C_N_M(2m) * f_M(2m) * |U| if Ri < 0 and z0m > 2 m
                  *(1._dp+pgeom1(jl,klev)/(g*zepz0o))))),                  &
                   zcdnr*zcml(jl)*SQRT(zdu2(jl)),lo.AND.zril(jl).LT.0._dp)
     zustl=zcfm2m*SQRT(zdu2(jl))                                             ! C_M(2m) * |U|^2 - where (2m) stands for z0 limited to 2m
     zustarl(jl)=SQRT(zustl)                                                 ! u* = sqrt(C_M(2m)) * |U|
!
!    water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     lo=paz0w(jl).GT.zepz0o
     zcdn2m=MERGE((ckap/LOG(1._dp+pgeom1(jl,klev)/(g*zepz0o)))**2,         & ! C_N_M(2m) - Equal to C_N_M, but adapted if z0m > 2m
                   zcdnw(jl),lo)
     zcdnr=zcdn2m/zcdnw(jl)                                                  ! C_N_M(2m)/C_N_M
     zcfm2m=MERGE(zcdn2m*SQRT(zdu2oc(jl))*(1._dp-2._dp*cc*zriw(jl))        & ! C_M * (C_N_M(2m)/C_N_M) * |U| if Ri>= 0 or z0m <= 2 m
                 /(1._dp+3._dp*(cc**2)*zcdn2m*SQRT(ABS(zriw(jl))           & ! C_N_M(2m) * f_M(2m) * |U| if Ri < 0 and z0m > 2 m
                 *(1._dp+pgeom1(jl,klev)/(g*zepz0o)))),                    &
                  zcdnr*zcmw(jl)*SQRT(zdu2oc(jl)),lo.AND.zriw(jl).LT.0._dp)
     zustw=zcfm2m*SQRT(zdu2oc(jl))                                           ! C_M(2m) * |U|^2 - where (2m) stands for z0 limited to 2m
     zustarw(jl)=SQRT(zustw)                                                 ! u* = sqrt(C_M(2m)) * |U|
!
!    ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     lo=paz0i(jl).GT.zepz0o
     zcdn2m=MERGE((ckap/LOG(1._dp+pgeom1(jl,klev)/(g*zepz0o)))**2,         & ! C_N_M(2m) - Equal to C_N_M, but adapted if z0m > 2m
                   zcdni(jl),lo)
     zcdnr=zcdn2m/zcdni(jl)                                                  ! C_N_M(2m)/C_N_M
     zcfm2m=MERGE(zcdn2m*SQRT(zdu2oc(jl))*(1._dp-2._dp*cc*zrii(jl))        & ! (C_M * C_N_M(2m)/C_N_M) * |U| if Ri>= 0 or z0m <= 2 m
                 /(1._dp+3._dp*(cc**2)*zcdn2m*SQRT(ABS(zrii(jl))           & ! C_N_M(2m) * f_M(2m) * |U| if Ri < 0 and z0m > 2 m
                 *(1._dp+pgeom1(jl,klev)/(g*zepz0o)))),                    &
                  zcdnr*zcmi(jl)*SQRT(zdu2oc(jl)),lo.AND.zrii(jl).LT.0._dp)
     zusti=zcfm2m*SQRT(zdu2oc(jl))                                           ! C_M(2m) * |U|^2 - where (2m) stands for z0 limited to 2m
     zustari(jl)=SQRT(zusti)                                                 ! u* = sqrt(C_M(2m)) * |U|
!
     zust=pfrl(jl)*zustl+pfrw(jl)*zustw+pfri(jl)*zusti                       ! Weighted u*^2
!
     zustarm=SQRT(zust)                                                      ! sqrt( weighted u*^2) => effective u*
     zust_2d(jl,krow) = zustarm
!
!    updating z0
!
     paz0w(jl)=MAX(1.5e-05_dp,cchar*zcmw(jl)*zdu2oc(jl)/g)                   ! z0 over sea set as max(0.018 u*^2 / g, 1.5e-5 m), using u*^2 = C_M |U_{relative to surface}|^2
     paz0i(jl)=cz0ice                                                        ! z0 over ice is set to the constant 1 mm every time step
     paz0(jl)=pfrl(jl)*paz0l(jl)+pfrw(jl)*paz0w(jl)+pfri(jl)*paz0i(jl)       ! Average z0
!
!    windstress
!
!     land   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     pustrl(jl)=zdens(jl)*zcml(jl)*sqrt(zdu2(jl))*pum1(jl,klev)              ! air density * C_M * |U| * u = rho u*^2 u / |U| = - rho u'w'
     pvstrl(jl)=zdens(jl)*zcml(jl)*sqrt(zdu2(jl))*pvm1(jl,klev)              ! air density * C_M * |U| * v = rho u*^2 v / |U| = - rho v'w'
!
!     water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     pustrw(jl)=zdens(jl)*zcmw(jl)*sqrt(zdu2oc(jl))*(pum1(jl,klev)-pocu(jl)) ! air density * C_M * |U_r| * u_r = rho u*^2 u_r / |U_r| = - rho u'w', _r denotes relative to ocean surface movement
     pvstrw(jl)=zdens(jl)*zcmw(jl)*sqrt(zdu2oc(jl))*(pvm1(jl,klev)-pocv(jl)) ! air density * C_M * |U_r| * v_r = rho u*^2 v_r / |U_r| = - rho v'w', _r denotes relative to ocean surface movement
!
!     ice   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     pustri(jl)=zdens(jl)*zcmi(jl)*sqrt(zdu2oc(jl))*(pum1(jl,klev)-pocu(jl)) ! air density * C_M * |U_r| * u_r = rho u*^2 u_r / |U_r| = - rho u'w', _r denotes relative to ocean surface movement
     pvstri(jl)=zdens(jl)*zcmi(jl)*sqrt(zdu2oc(jl))*(pvm1(jl,klev)-pocv(jl)) ! air density * C_M * |U_r| * v_r = rho u*^2 v_r / |U_r| = - rho v'w', _r denotes relative to ocean surface movement
!
!     average   +++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     pustr(jl)=pfrl(jl)*pustrl(jl)+pfrw(jl)*pustrw(jl)+pfri(jl)*pustri(jl)   ! Surface stress in zonal direction
     pvstr(jl)=pfrl(jl)*pvstrl(jl)+pfrw(jl)*pvstrw(jl)+pfri(jl)*pvstri(jl)   ! Surface stress in meridonal direction
!
     zcor=MAX(ABS(coriol_2d(jl,krow)),5.e-05_dp)                             ! |fc| - Absolute Coriolis coefficient
     zhdyn(jl)=MIN(pgeom1(jl,1)/g,zchneu*zustarm/zcor)                       ! Minimum boundary-layer height, based on dynamics (h_dyn = 0.3 u*/|f_c|)
!
     ihpblc(jl)=klev                                                         ! Initial setting for convective boundary-layer top at (physically) lowest level
     ihpbld(jl)=klev                                                         ! Initial setting for dynamic boundary-layer top at (physically) lowest level
  END DO
!
  DO jk=klevm1,1,-1
     DO jl=1,kproma
        zds=zcptgz(jl,jk)-zcptgz(jl,klev)                                    ! Difference in dry static energy between evaluated level and surface
        zdz=pgeom1(jl,jk)/g-zhdyn(jl)                                        ! Evaluated height minus the dynamic boundary-layer height
        ihpblc(jl)=MERGE(jk,ihpblc(jl),ihpblc(jl).EQ.klev.AND.zds.GT.0._dp)  ! If this is the lowest (physical) level (i.e. highest index) with s higher than near the surface, convective BL top is at this level
        ihpbld(jl)=MERGE(jk,ihpbld(jl),ihpbld(jl).EQ.klev.AND.zdz.GE.0._dp)  ! If this is the lowest (physical) level (i.e. highest index) for which z>=h_dyn, dynamic BL top is at this level
     END DO
  END DO
!
!      CONVECTIVE VELOCITY SCALE, MONIN-OBUKHOV LENGTH AND
!      TKE BOUNDARY CONDITION (MAILHOT/BENOIT, 1982)
!
  DO jl=1,kproma
     ihpbl(jl)=MIN(ihpblc(jl),ihpbld(jl))                                    ! Boundary-layer top at (physically) highest grid level between dynamic and convective BL top
     zghabl=MIN(50000._dp,pgeom1(jl,ihpbl(jl)))                              ! Geopotential height of boundary-layer top, g * zi
!    Calculating E0; needed for calculation of K_M and K_H (Mailhot and Benoit, 1982)
!
!    land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     IF(zwstl(jl).GT.zepsr) THEN                                             ! w'th_v'/ C_H / th_v bigger > 0 => convective
        zconvs=(zwstl(jl)*zchl(jl)*zghabl)**(1._dp/3._dp)                    ! w* = (w'th_v'* g/th_v * zi)^(1/3)
        zmonob=(zustarl(jl)**3)/(ckap*g*zwstl(jl)*zchl(jl))                  ! Negative Obukhov length, -L  = u*^3 / (kappa * g/th_v * w'th_v')
        zstabf=(pgeom1(jl,klev)/(g*zmonob))**(2._dp/3._dp)                   ! Part of calculation of E0: (-z/L)^(2/3)
        zstabf=MIN(zustf*3._dp,zstabf)                                       ! Limit (-z/L)^(2/3) to not exceed a value of roughly 10 (-z/L < 30)
     ELSE
        zconvs=0._dp
        zstabf=0._dp
     END IF
     ztkevl=(zustf+zstabf)*(zustarl(jl)**2)+zwstf*(zconvs**2)                ! E0, TKE at surface
!
!    water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     IF(zwstw(jl).GT.zepsr) THEN                                             ! w'th_v'/ C_H / th_v bigger > 0
        zconvs=(zwstw(jl)*zchw(jl)*zghabl)**(1._dp/3._dp)                    ! w* = (w'th_v'* g/th_v * zi)^(1/3)
        zmonob=(zustarw(jl)**3)/(ckap*g*zwstw(jl)*zchw(jl))                  ! Negative Obukhov length, -L  = u*^3 / (kappa * g/th_v * w'th_v')
        zstabf=(pgeom1(jl,klev)/(g*zmonob))**(2._dp/3._dp)                   ! Part of calculation of E0: (-z/L)^(2/3)
        zstabf=MIN(zustf*3._dp,zstabf)                                       ! Limit (-z/L)^(2/3) to not exceed a value of roughly 10 (-z/L < 30)
     ELSE
        zconvs=0._dp
        zstabf=0._dp
     END IF
     ztkevw=(zustf+zstabf)*(zustarw(jl)**2)+zwstf*(zconvs**2)                ! E0, TKE at surface
!
!    ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     IF(zwsti(jl).GT.zepsr) THEN                                             ! w'th_v'/ C_H / th_v bigger > 0
        zconvs=(zwsti(jl)*zchi(jl)*zghabl)**(1._dp/3._dp)                    ! w* = (w'th_v'* g/th_v * zi)^(1/3)
        zmonob=(zustari(jl)**3)/(ckap*g*zwsti(jl)*zchi(jl))                  ! Negative Obukhov length, -L  = u*^3 / (kappa * g/th_v * w'th_v')
        zstabf=(pgeom1(jl,klev)/(g*zmonob))**(2._dp/3._dp)                   ! Part of calculation of E0: (-z/L)^(2/3)
        zstabf=MIN(zustf*3._dp,zstabf)                                       ! Limit (-z/L)^(2/3) to not exceed a value of roughly 10 (-z/L < 30)
     ELSE
        zconvs=0._dp
        zstabf=0._dp
     END IF
     ztkevi=(zustf+zstabf)*(zustari(jl)**2)+zwstf*(zconvs**2)                ! E0, TKE at surface
     ztkevn(jl,klev)=pfrl(jl)*ztkevl+pfrw(jl)*ztkevw+pfri(jl)*ztkevi         ! Weighted TKE at the surface
     ztkevn(jl,klev)=MAX(ztkemin,ztkevn(jl,klev))                            ! Prevent 0 TKE
  END DO
!
  IF(lstart) THEN
     DO jl=1,kproma
        ptkem1(jl,klev)=ztkevn(jl,klev)                                      ! surface TKE at previous timestep
        ptkem(jl,klev)=ztkevn(jl,klev)                                       ! surface TKE at current timestep
     END DO
  END IF
!
!     ==================================================================
!
!*       3.5   Vertical loop: Computation of basic quantities:
!              wind shear, buoyancy, Ri-number, mixing length
!
  DO jk=1,klevm1
     DO jl=1,kproma
        zqtmit=zlwcmit(jl,jk)+zqmit(jl,jk)                                   ! Total atmospheric water content at interface level
        zfux=zfaxen(jl,jk)/(cpd*ztmitte(jl,jk))                              ! L / (cp_d T) at interface level
        zfox=zfaxen(jl,jk)/(rd*ztmitte(jl,jk))                               ! L / (Rd T) at interface level
        zmult1=1._dp+vtmpc1*zqtmit                                           ! Calculating intermediates for th_v' for sat. air: 1 + (Rv/Rd - 1) q_t at interface level
        zmult2=zfux*zmult1-zrvrd                                             ! Calculating intermediates for th_v' for sat. air: L (1 + (Rv/Rd - 1) q_t) / (cp_d T) - Rv/Rd at interface level
        zmult3=zrdrv*zfox*zqssm(jl,jk)/(1._dp+zrdrv*zfux*zfox*zqssm(jl,jk))  ! Calculating intermediates for th_v' for sat. air: L q_s /(Rv T (1 + L / (Rv T) L / (cp_d T) q_s ) )
        zmult5=zmult1-zmult2*zmult3                                          ! c1 for sat. air in the equation th_v' = c1 th_l' + c_2 th q'
        zmult4=zfux*zmult5-1._dp                                             ! c2 for sat. air in equation above: c2 = c1 L / (cp_d T) - 1
        zdus1=zccover(jl,jk)*zmult5+(1._dp-zccover(jl,jk))*zmult1            ! c1 for entire interface level (saturated where cloud)
        zdus2=zccover(jl,jk)*zmult4+(1._dp-zccover(jl,jk))*vtmpc1            ! c2 for entire interface level (saturated where cloud)
        zteldif=(zlteta1(jl,jk)-zlteta1(jl,jk+1))/zhh(jl,jk)                 ! Gradient in liquid potential temperature between adjacent full levels
        zdqtot=(pqm1(jl,jk)+zx(jl,jk))-(pqm1(jl,jk+1)+zx(jl,jk+1))           ! Difference in total atmospheric water content between adjacent full levels
        zqddif=zdqtot/zhh(jl,jk)                                             ! Gradient in total atmospheric water content between adjacent full levels
        zqshear(jl,jk)=zqddif                                                ! Store for variance production
        zbuoy=(zteldif*zdus1+ztemit(jl,jk)*zdus2*zqddif)*g/ztvirmit(jl,jk)   ! Gradient in virtual potential temperature between adjacent levels times g/th_v
        zdivv=(pum1(jl,jk)-pum1(jl,jk+1))**2                                 ! Square of difference in u between adjacent levels
        zdivv1=(pvm1(jl,jk)-pvm1(jl,jk+1))**2                                ! Square of difference in v between adjacent levels
        zshear=(zdivv+zdivv1)/zhh(jl,jk)**2                                  ! Sum of squares of horiz. velocity (vertical) gradients between adjacent levels
        zri=zbuoy/MAX(zshear,1.e-5_dp)                                       ! Gradient Richardson number
        rinum_3d(jl,jk,krow) = zri
!
!       ASYMPTOTIC MIXING LENGTH FOR MOMENTUM AND
!       HEAT (ZLAM) ABOVE THE PBL AS A FUNCTION OF HEIGHT
!       ACCORDING TO HOLTSLAG AND BOVILLE (1992), J. CLIMATE.
!
        IF(jk.GE.ihpbl(jl)) THEN
           zlam  = clam                                                      ! Asymptotic mixing length, lambda, within the boundary layer - 150 m
        ELSE
           zhexp = EXP(1._dp-pgeom1(jl,jk)/pgeom1(jl,ihpbl(jl)))             ! Exponential factor in expression for asymptotic mixing length, lambda, above the boundary layer
           zlam  = 1._dp+(clam-1._dp)*zhexp                                  ! lambda(z>zi) = 1 + 149 * exp(1 - z/zi) m - determined at upper adjacent full level
        END IF
!
!       MIXING LENGTH (BLACKADAR) + STABILITY DEPENDENT FUNCTION
!
        zz2geo = (0.5_dp*ckap/g)*(pgeom1(jl,jk)+pgeom1(jl,jk+1))             ! kappa z at interface level
        zmix   = 1._dp / ( (1._dp/zz2geo) + (1._dp/zlam) )                   ! mixing length - l = 1 / ( 1/(kappa*z) + 1/lambda )
!
!       STABILITY FUNCTIONS (LOUIS, 1979)
!
        IF(zri.LT.0._dp) THEN                                                ! Stability function, g_H,M, for unstable conditions - S_H,M = g_H,M S_N_H,M
           zucf=1._dp/(1._dp+3._dp*((cc*g*zmix)**2)*SQRT(ABS(zri)          & ! Intermediate factor for g_H,M, based on height of lower adjacent full level
                *(((pgeom1(jl,jk)/pgeom1(jl,jk+1))**(1._dp/3._dp)          &
                -1._dp)/(pgeom1(jl,jk)-pgeom1(jl,jk+1)))**3/pgeom1(jl,jk+1)))
           zsh=zshn*(1._dp-3._dp*cc*zri*zucf)*zmix                           ! Lambda_H = S_H * l = S_N_H * g_H * l
           zsm=zsmn*(1._dp-2._dp*cc*zri*zucf)*zmix                           ! Lambda_M = S_M * l = S_N_M * g_M * l
        ELSE                                                                 ! Stable conditions
           zsh=zshn/(1._dp+2._dp*cc*zri*SQRT(1._dp+zri))*zmix                ! Lambda_H
           zsm=zsmn/(1._dp+2._dp*cc*zri/SQRT(1._dp+zri))*zmix                ! Lambda_M
        END IF
!
!       Dimensionless coefficients multiplied by pressure
!            thicknesses for momentum and heat exchange
!
        zzb=zshear*zsm-zbuoy*zsh                                             ! 1/sqrt(E) * {dE/dt due to shear production and buoyancy production/destruction} (Lambda_M,H = K_M,H/sqrt(E))
        zdisl=(zmix/zsmn**3)/ztmst                                           ! Lambda_1 / leap frog time step: Dissipation rate, epsilon, is E**1.5 / Lambda_1
        zktest=1._dp+(zzb*ztmst+SQRT(ptkem1(jl,jk))*2._dp)/zdisl             ! Part of implicit calculation of sqrt(E(t+1)) (Brinkop and Roeckner, 1995: Eq. (A10))
        IF (zktest.LE.1._dp) THEN                                            ! Calculated sqrt(E) would become negative
           ztkevn(jl,jk)=ztkemin
        ELSE
           ztkevn(jl,jk)=MAX(ztkemin,(zdisl*(SQRT(zktest)-1._dp))**2)        ! Square of Eq. (A10)
        END IF
        IF(lstart) THEN
           ptkem1(jl,jk)=ztkevn(jl,jk)                                       ! TKE at previous timestep
           ptkem(jl,jk)=ztkevn(jl,jk)                                        ! TKE at current timestep
        END IF
        ztkesq=SQRT(MAX(ztkemin,ptkem1(jl,jk)))                              ! sqrt(E(t-1))
        zztvm=(ptvm1(jl,jk)+ptvm1(jl,jk+1))*0.5_dp                           ! Interpolated T_v at interface level
        zalf=paphm1(jl,jk+1)/(zztvm*zhh(jl,jk)*rd)                           ! Air density at interface level divided by distance between adjacent full levels
!       Factors A to be used for the Richtmeyer and Morton scheme (Schulz et al., 2001: Appendix)
        IF (.NOT.ledith) THEN
           zcfm(jl,jk)=cvdifts*ztmst*g*zalf*zsm*ztkesq                       !  Computational factor 1.5 * leap frog time step * g * air density * K_M / Delta{z}, K_M = Lambda_M * sqrt(E(t-1))
           zcfh(jl,jk) =cvdifts*ztmst*g*zalf*zsh*ztkesq                      ! Computational factor 1.5 * leap frog time step * g * air density * K_H / Delta{z}, K_H = Lambda_H * sqrt(E(t-1))
        ELSE

           ! including turbulent diffusion by hines scheme
           pgwdiffco(jl,jk)=gweddy(jl,jk,krow)
           zcfm(jl,jk) = cvdifts*ztmst*g*zalf*(zsm*ztkesq+pgwdiffco(jl,jk))
           zcfh(jl,jk) = cvdifts*ztmst*g*zalf*(zsm*ztkesq+pgwdiffco(jl,jk)*0.72_dp) ! zprinv: inverse of prandtl number; for testing fixed here

        ENDIF
        zcfv(jl,jk) =0.5_dp*zcfh(jl,jk)                                      ! Just half of the above, probably since it is used for transport of a standard deviation (proxy) instead of variance
        zcdum(jl,jk)=cvdifts*ztmst*g*zalf*zsm*SQRT(ztkevn(jl,jk))            ! Computational factor 1.5 * leap frog time step * g * air density * K_M~ / Delta{z}, K_M~ = Lambda_M * sqrt(E(t+1))
     END DO
  END DO
!
!     ==================================================================
!
!*       3.8        DIFFUSION IMPLICIT COMPUTATIONS FOR TKE
!
!
! TKE defined at interface levels, thus exchange at full level
! Equations, based on Richtmeyer and Morton (1967), are specified in Schulz et al., 2001
! This is a special adaptation of the Gauss elimination procedure
! Note that X_(t-1) here corresponds to X_(t+1)** of Brinkop and Roeckner (1995)
! i.e. including all tendencies except for vertical diffusion
  DO jk=1,klev
     DO jl=1,kproma
        zedif(jl,jk)=ztkevn(jl,jk)/cvdifts                                   ! E(t+1)** / computational factor 1.5 = X_(t-1) / alpha in Schulz et al. (2001)
     END DO
  END DO
!
  DO jl=1,kproma
     ztcoe(jl)=(zcdum(jl,1)+zcdum(jl,2))*0.5_dp                              ! A_(k+1/2) with respect to highest evaluated interface level (messy interface layer 2, here index 1)
     zqdp=1._dp/(papm1(jl,2)-papm1(jl,1))                                    ! 1 / Delta{p}
     zdisc=1._dp/(1._dp+ztcoe(jl)*zqdp)                                      ! 1 / ( 1 + A_(k+1/2) / Delta{p} ) = 1 / ( 1 + A^_(k+1/2) )
     zebsm(jl,1)=zdisc*ztcoe(jl)*zqdp                                        ! A^_(k+1/2) / ( 1 + A^_(k+1/2) ) = E_(k+1/2)
     zedif(jl,1)=zdisc*zedif(jl,1)                                           ! F_(k+1/2)
  END DO
!
  DO jk=2,klev-2
     DO jl=1,kproma
        zqdp=1._dp/(papm1(jl,jk+1)-papm1(jl,jk))                             ! 1 / Delta{p}
        zfac=ztcoe(jl)*zqdp                                                  ! A_(k-1/2) / Delta{p} = C^
        ztcoe(jl)=(zcdum(jl,jk+1)+zcdum(jl,jk))*0.5_dp                       ! A_(k+1/2) compared to interface layer k (messy: interface layer k+1)
        zdisc=1._dp/(1._dp+zfac*(1._dp-zebsm(jl,jk-1))                     & ! 1 / ( 1 + C^ (1 - E_(k-1/2)) + A^_(k+1/2) ) = 1 / ( B^ - C^  E_(k-1/2) )
                  +ztcoe(jl)*zqdp)                                           !      with B^ = 1 + C^ + A^_(k+1/2)
        zebsm(jl,jk)=zdisc*ztcoe(jl)*zqdp                                    ! A^_(k+1/2) / ( B^ - C^  E_(k-1/2) ) = E_(k+1/2)
        zedif(jl,jk)=zdisc*(zedif(jl,jk)+zfac*zedif(jl,jk-1))                ! ( X_(t-1) / alpha + C^ F_(k-1/2) ) / ( B^ - C^  E_(k-1/2) ) = F_(k+1/2)
     END DO
  END DO
!
  DO jl=1,kproma                                                             ! At level klevm1: lowest inteface level above the surface
     zqdp=1._dp/(papm1(jl,klev)-papm1(jl,klevm1))                            ! 1 / Delta{p}
     zfac=ztcoe(jl)*zqdp                                                     ! C^
     ztcoe(jl)=(zcdum(jl,klev)+zcdum(jl,klevm1))*0.5_dp                      ! A_(k+1/2)
     zdisc=1._dp/(1._dp+zfac*(1._dp-zebsm(jl,klev-2))+ztcoe(jl)*zqdp)        ! 1 / ( B^ - C^  E_(k-1/2) )
     zedif(jl,klevm1)=zdisc*(ztcoe(jl)*zqdp*zedif(jl,klev)                 & ! Not F_(k+1/2), but X*/alpha
                      +zedif(jl,klevm1)                                    & !   E_(k+1/2) * X_(t-1,k+1) / alpha + F_(k+1/2) = X*_k/alpha,
                      +zfac*zedif(jl,klev-2))                                !   since TKE at surface is unaltered by diffusion at surface X*_(k+1) =  X_(t-1,k+1)
  END DO
!
  DO jk=klev-2,1,-1
     DO jl=1,kproma
        zedif(jl,jk)=zedif(jl,jk)+zebsm(jl,jk)*zedif(jl,jk+1)                ! X*_k/alpha = E_(k+1/2) * X*_(k+1)/alpha + F_(k+1/2)
     END DO
  END DO
!
!*    TIME INTEGRATION OF TURBULENT KINETIC ENERGY AND CHECK
!
  DO jk=1,klev
     ztest=0._dp
     DO jl=1,kproma                                                          ! Using X as shorthand notation for E
        ptke(jl,jk)=zedif(jl,jk)+(1._dp-(1._dp/cvdifts))*ztkevn(jl,jk)       !   X_(k,t+1) = X*_k/alpha + (1 - 1/alpha) * X_(k,t-1)
        ztest=ztest+MERGE(1._dp,0._dp,ptke(jl,jk)<0._dp)                     ! Check for negative resulting TKE
     END DO
     IF(ztest.NE.0._dp) CALL error_bi('TKE IS NEGATIVE','vertex')
  END DO
!
!*    TIME FILTER FOR TURBULENT KINETIC ENERGY
!
  IF(.NOT.lstart) THEN
    zeps=eps
  ELSE
    zeps=0._dp
  END IF
  DO jk=1,klev
    DO jl=1,kproma
       ptkem1(jl,jk)=ptkem(jl,jk)                                          & ! Robert-Asselin filter for TKE at t-1
                +zeps*(ptkem1(jl,jk)-2._dp*ptkem(jl,jk)+ptke(jl,jk))
       ptkem(jl,jk)=ptke(jl,jk)                                              ! Setting TKE at t
     END DO
  END DO

!******************************************************************************

  ! Reset pointers used only in first call to vertex_vdiff
  NULLIFY(ptke    )
  NULLIFY(ptkem   )
  NULLIFY(zteta1  )

  NULLIFY(paz0    )
  NULLIFY(paz0hi  )
  NULLIFY(paz0hw  )
  NULLIFY(paz0i   )
  NULLIFY(paz0w   )
  NULLIFY(pustr   )
  NULLIFY(pustri  )
  NULLIFY(pustrl  )
  NULLIFY(pustrw  )
  NULLIFY(pvstr   )
  NULLIFY(pvstri  )
  NULLIFY(pvstrl  )
  NULLIFY(pvstrw  )
  NULLIFY(zcdni   )
  NULLIFY(zcdnl   )
  NULLIFY(zcdnw   )
  NULLIFY(ztvi    )
  NULLIFY(ztvl    )
  NULLIFY(ztvw    )
  NULLIFY(zwet_tmp)

  END SUBROUTINE vertex_vdiff_1
!------------------------------------------------------------------------------
  SUBROUTINE vertex_vdiff_2

#ifndef MESSYTENDENCY
  USE messy_main_tracer_mem_bi, ONLY: xtm1
#endif
  USE messy_main_data_bi,       ONLY: xvar                                   ! Only defined in ECHAM5 & CESM1
  USE messy_main_data_bi,       ONLY: u10, v10                               ! Not defined in BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: grndcapc, grndhflx,tsl, snc            ! Only defined in ECHAM5 & CESM1
  USE messy_main_data_bi,       ONLY: vgrat                                  ! Not defined in BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: tslnew,   loglac_2d                    ! Only defined in ECHAM5 & CESM1
  USE messy_main_data_bi,       ONLY: qflux                                  ! Not defined in VERTICO, BLANK, MBM_CLAMS
  USE messy_main_constants_mem, ONLY: rhoh2o=>rho_H2O, stbo, cemiss, cwlmax
  USE messy_main_tracer_mem_bi, ONLY: ti_gp, ON, I_Vdiff
  USE messy_main_data_bi,       ONLY: rh_2m                                  ! Not defined in BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: zsenkf_2d,  zlatkf_2d                  ! Not defined in VERTICO, BLANK, MBM_CLAMS
  USE messy_main_data_bi,       ONLY: vmixtau,    vdiffp                     ! Only defined in ECHAM5 & CESM1
  USE messy_main_data_bi,       ONLY: ebsh,  qslnew, cair                    ! Only defined in ECHAM5
  USE messy_main_data_bi,       ONLY: pxtte => xtte                          ! Only correctly defined in ECHAM5

  IMPLICIT NONE

!******************************************************************************

  ! Temporary variables declared
  INTEGER :: it, jk, jl, jt

  LOGICAL :: lo, lo1, lookupoverflow

  REAL(dp):: z2geomf
  REAL(dp):: zaph2m
  REAL(dp):: zcbn, zcbs, zcbu
  REAL(dp):: zepot, zsnfac, zwlfac
  REAL(dp):: zcoefi, zcoefl, zcoefw
  REAL(dp):: zcpt, zcvm3
  REAL(dp):: zcvm4, zdew2i, zdew2l, zdew2w
  REAL(dp):: zdisc, zdisci, zdiscl, zdiscw
  REAL(dp):: zdisql, zdisx, zdisxt
  REAL(dp):: zdqdt, zdqtot, zdtdt, zdudt
  REAL(dp):: zdvdt, zdximdt, zdxlmdt, zdxtdt
  REAL(dp):: zephum
  REAL(dp):: zfac, zfrac
  REAL(dp):: zh2m, zhexp, zhtq
  REAL(dp):: zgtl, zgtw, zgti, zgtsum, zgql, zgqw, zgqi, zgqsum
  REAL(dp):: zhuv, zlam
  REAL(dp):: zmix
  REAL(dp):: zq2m, zqdp, zqklevi
  REAL(dp):: zqklevl, zqklevw, zqnlev
  REAL(dp):: zqs1, zqs2, zqvhfl, zrat
  REAL(dp):: zspeedi
  REAL(dp):: zred, zrh2m
  REAL(dp):: zspeedl, zspeedw, zt2i, zt2l
  REAL(dp):: zt2w, ztkemin, zepevap
  REAL(dp):: ztkesq, ztklevi, ztklevl
  REAL(dp):: ztklevw, ztnlev
  REAL(dp):: zqnew, zu10i, zu10l, zu10w
  REAL(dp):: zv10i, zv10l, zv10w
  REAL(dp):: zzcpts, zzqs
  REAL(dp):: zcptlcorr
  REAL(dp):: zqsurf, zrhodz
  REAL(dp):: ztvlan, ztvsea, ztvice, ztvh
  REAL(dp):: zdisv, zfav

#ifndef MESSYTENDENCY
  REAL(dp):: pxtm1    (kbdim,klev,ktrac)
#endif

  REAL(dp):: zxtvn    (kbdim,klev,ktrac), zxtdif   (kbdim,klev,ktrac)

  REAL(dp):: pemter   (kbdim,klevp1),     zqflux   (kbdim,klevp1),  &
             zrho     (kbdim,klevp1),     zvarpr   (kbdim,klevp1)

  REAL(dp):: zdis     (kbdim,klev),       zebsm    (kbdim,klev),    &
             zebsv    (kbdim,klev),       zqdif    (kbdim,klev),    &
             ztdif    (kbdim,klev),       zudif    (kbdim,klev),    &
             zvardif  (kbdim,klev),       zvdif    (kbdim,klev),    &
             zxidif   (kbdim,klev),       zxldif   (kbdim,klev)

  REAL(dp):: pgrndcapc(kbdim),            pgrndhflx(kbdim),         &
             pvgrat   (kbdim),            pwlmx    (kbdim),         &
             zcptl    (kbdim),            zcsat    (kbdim),         &
             zqhfli   (kbdim),            zqhfll   (kbdim),         &
             zqhflw   (kbdim),            zthfli   (kbdim),         &
             zthfll   (kbdim),            zthflw   (kbdim),         &
             zvidis   (kbdim)


  ! Pointers to channel objects for use in second call to vertex_vdiff
  REAL(dp), POINTER, DIMENSION(:,:) :: pvdiffp   => NULL()
  REAL(dp), POINTER, DIMENSION(:,:) :: pvmixtau  => NULL()
  REAL(dp), POINTER, DIMENSION(:,:) :: pxvar     => NULL()
  REAL(dp), POINTER, DIMENSION(:,:) :: zebsh     => NULL()

  REAL(dp), POINTER, DIMENSION(:)   :: pahfl    => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pahfli   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pahfll   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pahflw   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pahfs    => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pahfsi   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pahfsl   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pahfsw   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pdew2    => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pevap    => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pevapi   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pevapl   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pevapot  => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pevapw   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pqhfla   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: psnc     => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: ptemp2   => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: ptsl     => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: ptslnew  => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pu10     => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pvdis    => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pv10     => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pwind10  => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pwind10w => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: zcair    => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: zqslnew  => NULL()

  REAL(dp), POINTER, DIMENSION(:)   :: pahfslac => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pahfswac => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pahfsiac => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pahfllac => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pahflwac => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pahfliac => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pevaplac => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pevapwac => NULL()
  REAL(dp), POINTER, DIMENSION(:)   :: pevapiac => NULL()


!******************************************************************************

! Special groups of local variables

#ifdef MESSYTENDENCY
  ! Variables needed to store tendencies if MESSYTENDENCY is enabled
  REAL(dp),DIMENSION(kbdim,klev)       :: zdqdt_2d
  REAL(dp),DIMENSION(kbdim,klev)       :: zdtdt_2d
  REAL(dp),DIMENSION(kbdim,klev)       :: zdudt_2d
  REAL(dp),DIMENSION(kbdim,klev)       :: zdvdt_2d
  REAL(dp),DIMENSION(kbdim,klev)       :: zdximdt_2d
  REAL(dp),DIMENSION(kbdim,klev)       :: zdxlmdt_2d
  REAL(dp),DIMENSION(kbdim,klev,ktrac) :: zdxtdt_3d
#endif

  ! Variables needed for the subroutine SURFTEMP
  REAL(dp):: zcdrag  (kbdim) ! Drag coefficient for heat and moisture
  REAL(dp):: zcpq    (kbdim) ! Specific heat of air as used in *vdiff*
  REAL(dp):: zcptlnew(kbdim) ! New surface dry static energy
  REAL(dp):: zeqni   (kbdim) ! Richtmyer-morton-coefficients
  REAL(dp):: zeqnl   (kbdim) !   for moisture.
  REAL(dp):: zeqnw   (kbdim) !
  REAL(dp):: zetni   (kbdim) ! Richtmyer-morton-coefficients
  REAL(dp):: zetnl   (kbdim) !   for dry static energy.
  REAL(dp):: zetnw   (kbdim) !
  REAL(dp):: zfqni   (kbdim) ! Richtmyer-morton-coefficients
  REAL(dp):: zfqnl   (kbdim) !   for moisture.
  REAL(dp):: zfqnw   (kbdim) !
  REAL(dp):: zftni   (kbdim) ! Richtmyer-morton-coefficients
  REAL(dp):: zftnl   (kbdim) !   for dry static energy.
  REAL(dp):: zftnw   (kbdim) !
  REAL(dp):: znetr   (kbdim) ! Surface net radiation (old)

  ! Variables only used for WARNINGs
  REAL(dp):: t2i_1d  (kbdim)
  REAL(dp):: t2l_1d  (kbdim)
  REAL(dp):: t2w_1d  (kbdim)

!******************************************************************************

  ! Set (initial) values of arrays
#ifndef MESSYTENDENCY
  pxtm1     =  xtm1     (:,:,:,krow)
#endif

  pemter    =  emter    (:,:,krow)

  pgrndcapc =  grndcapc (:,krow)
  pgrndhflx =  grndhflx (:,krow)
  pvgrat    =  vgrat    (:,krow)

  ! Setting of pointers that are used in second call to vertex_vdiff
  pvdiffp   => vdiffp   (:,       :,krow)
  pvmixtau  => vmixtau  (:,       :,krow)
  pxvar     => xvar     (:,       :,krow)
  zebsh     => ebsh     (1:kproma,:,krow)

  pahfl     => ahfl     (:,         krow)
  pahfli    => ahfli    (:,         krow)
  pahfll    => ahfll    (:,         krow)
  pahflw    => ahflw    (:,         krow)
  pahfs     => ahfs     (:,         krow)
  pahfsi    => ahfsi    (:,         krow)
  pahfsl    => ahfsl    (:,         krow)
  pahfsw    => ahfsw    (:,         krow)
  pdew2     => dew2     (:,         krow)
  pevap     => evap     (:,         krow)
  pevapi    => evapi    (:,         krow)
  pevapl    => evapl_2d (:,         krow)
  pevapot   => evapot_2d(:,         krow)
  pevapw    => evapw    (:,         krow)
  pqhfla    => qflux    (:,         krow)
  psnc      => snc      (:,         krow)
  ptemp2    => temp2    (:,         krow)
  ptsl      => tsl      (:,         krow)
  ptslnew   => tslnew   (:,         krow)
  pu10      => u10      (:,         krow)
  pvdis     => vdis     (:,         krow)
  pv10      => v10      (:,         krow)
  pwind10   => wind10   (:,         krow)
  pwind10w  => wind10w  (:,         krow)
  zcair     => cair     (1:kproma,  krow)
  zqslnew   => qslnew   (1:kproma,  krow)

  pahfslac  => ahfslac  (:,krow)
  pahfswac  => ahfswac  (:,krow)
  pahfsiac  => ahfsiac  (:,krow)
  pahfllac  => ahfllac  (:,krow)
  pahflwac  => ahflwac  (:,krow)
  pahfliac  => ahfliac  (:,krow)
  pevaplac  => evaplac  (:,krow)
  pevapwac  => evapwac  (:,krow)
  pevapiac  => evapiac  (:,krow)

!******************************************************************************

  lookupoverflow = .FALSE.

#ifdef MESSYTENDENCY
  zdtdt_2d(:,:)    = 0.0_dp
  zdqdt_2d(:,:)    = 0.0_dp
  zdudt_2d(:,:)    = 0.0_dp
  zdvdt_2d(:,:)    = 0.0_dp
  zdxlmdt_2d(:,:)  = 0.0_dp
  zdximdt_2d(:,:)  = 0.0_dp
  zdxtdt_3d(:,:,:) = 0.0_dp
#endif

!*      PARAMETERS FOR BOUNDARY LAYER DIAGNOSTICS
!
  zhuv      = 10._dp*g                        ! Geopotential height @ 10 m
  zhtq      = 2._dp*g                         ! Geopotential height @ 2 m

!*    SECURITY PARAMETERS.
!
  zepevap   = 1.e-10_dp                       ! Minimum relative humidity to be considered over land: 10^-10
  zephum    = 5.e-2_dp                        ! Minimum RH at 2 m height: 5 %
  ztkemin   = 1.e-10_dp                       ! Minimum TKE: 10^-10 m^2 s^-2

!
!     ------------------------------------------------------------------
!
!*       4.     DIFFUSION IMPLICIT COMPUTATIONS FOR MOMENTUM.
!
! Equations, based on Richtmeyer and Morton (1967), are specified in Schulz et al., 2001
!
!*       4.1     SETTING OF RIGHT HAND SIDES.
!
  DO jk=1,klev
     DO jl=1,kproma
        zudif(jl,jk)=pum1(jl,jk)/cvdifts                                     ! u(t+1)** / computational factor 1.5 = X_(t-1) / alpha (for u) in Schulz et al. (2001)
        zvdif(jl,jk)=pvm1(jl,jk)/cvdifts                                     ! v(t+1)** / computational factor 1.5 = X_(t-1) / alpha (for v) in Schulz et al. (2001)
     END DO
  END DO
!
!*       4.2     TOP LAYER ELIMINATION.
!
  DO jl=1,kproma                                                             ! Here k = 1; zcfm(k) = A_(k+1/2)
     zqdp=1._dp/(paphm1(jl,2)-paphm1(jl,1))                                  ! 1 / Delta{p} around level 1
     zdisc=1._dp/(1._dp+zcfm(jl,1)*zqdp)                                     ! 1 / ( 1 + A_(k+1/2) / Delta{p} ) = 1 / ( 1 + A^_(k+1/2) )
     zebsm(jl,1)=zdisc*(zcfm(jl,1)*zqdp)                                     ! A^_(k+1/2) / ( 1 + A^_(k+1/2) ) = E_(k+1/2) (for k = 1)
     zudif(jl,1)=zdisc*zudif(jl,1)                                           ! For u: X_(t-1) / alpha / ( 1 + A^_(k+1/2) ) = F_(k+1/2)
     zvdif(jl,1)=zdisc*zvdif(jl,1)                                           ! For v: X_(t-1) / alpha / ( 1 + A^_(k+1/2) ) = F_(k+1/2)
  END DO
!
!*       4.3     ELIMINATION FOR MIDDLE LAYERS.
!
  DO jk=2,klevm1
     DO jl=1,kproma
        zqdp=1._dp/(paphm1(jl,jk+1)-paphm1(jl,jk))                           ! 1 / Delta{p} around level k
        zfac=zcfm(jl,jk-1)*zqdp                                              ! A_(k-1/2) / Delta{p} = C^_(k+1/2)
        zdisc=1._dp/(1._dp+zfac*(1._dp-zebsm(jl,jk-1))+zcfm(jl,jk)*zqdp)     ! 1 / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) )
        zebsm(jl,jk)=zdisc*(zcfm(jl,jk)*zqdp)                                ! A^_(k+1/2) / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) ) = E_(k+1/2)
        zudif(jl,jk)=zdisc*(zudif(jl,jk)+zfac*zudif(jl,jk-1))                ! For u: (X_(t-1) / alpha + C^_(k+1/2) F_(k-1/2)) / (...) = F_(k+1/2)
        zvdif(jl,jk)=zdisc*(zvdif(jl,jk)+zfac*zvdif(jl,jk-1))                ! For v: (X_(t-1) / alpha + C^_(k+1/2) F_(k-1/2)) / (...) = F_(k+1/2)
     END DO
  END DO
!
!*       4.4     BOTTOM LAYER ELIMINATION.
!
  DO jl=1,kproma
     zqdp=1._dp/(paphm1(jl,klevp1)-paphm1(jl,klev))                          ! 1 / Delta{p} around level N
     zfac=zcfm(jl,klevm1)*zqdp                                               ! A_(k-1/2) / Delta{p} = C^_(k+1/2)
     zdisc=1._dp/(1._dp+zfac*(1._dp-zebsm(jl,klevm1))+zcfm(jl,klev)*zqdp)    ! 1 / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) )
     zudif(jl,klev)=zdisc*(zudif(jl,klev)+zfac*zudif(jl,klevm1)            & ! ( X_(t-1) / alpha + C^_(k+1/2) F_(k-1/2) + A^_(k+1/2) X_surf / alpha ) / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) )
                     + zcfm(jl,klev)*zqdp*pocu(jl)*(1._dp-pfrl(jl))/cvdifts) !   = F_(k+1/2) + E_(k+1/2) X_surf / alpha = X*_k/alpha
     zvdif(jl,klev)=zdisc*(zvdif(jl,klev)+zfac*zvdif(jl,klevm1)            & !   X_surf = surface velocity => 0 over land, water surface velocity elsewhere
                     + zcfm(jl,klev)*zqdp*pocv(jl)*(1._dp-pfrl(jl))/cvdifts)
  END DO
!
!*       4.5     BACK-SUBSTITUTION.
!
  DO jk=klevm1,1,-1
     DO jl=1,kproma
        zudif(jl,jk)=zudif(jl,jk)+zebsm(jl,jk)*zudif(jl,jk+1)                ! For u: X*_k/alpha = F_(k+1/2) + E_(k+1/2) * X*_(k+1)/alpha
        zvdif(jl,jk)=zvdif(jl,jk)+zebsm(jl,jk)*zvdif(jl,jk+1)                ! For v: X*_k/alpha = F_(k+1/2) + E_(k+1/2) * X*_(k+1)/alpha
     END DO
  END DO
!
!*       4.6     INCREMENTATION OF U AND V TENDENCIES AND STORAGE OF
!*               THE DISSIPATION.
!
  DO jl=1,kproma
     zvidis(jl)=0._dp                                                        ! Initialisation of dissipation to 0
  END DO
!
  DO jk=1,klev
     DO jl=1,kproma
        zdudt=(zudif(jl,jk)-pum1(jl,jk)/cvdifts)/ztmst                       ! Tendency in u = ( u(t+1) - u(t-1) ) / 2dt = ( X*/alpha - X(t-1)/alpha ) / 2dt
#ifndef MESSYTENDENCY
        pvom(jl,jk)=pvom(jl,jk)+zdudt
#else
        zdudt_2d(jl,jk)=zdudt
#endif
        zdvdt=(zvdif(jl,jk)-pvm1(jl,jk)/cvdifts)/ztmst                       ! Tendency in v = ( v(t+1) - v(t-1) ) / 2dt = ( X*/alpha - X(t-1)/alpha ) / 2dt
#ifndef MESSYTENDENCY
        pvol(jl,jk)=pvol(jl,jk)+zdvdt
#else
        zdvdt_2d(jl,jk)=zdvdt
#endif
        zdis(jl,jk)=0.5_dp*((pum1(jl,jk)/cvdifts-zudif(jl,jk))             & ! Dissipation calculated as 0.5 * (-Delta{u}*(u(t-1)+u(t+1)) + -Delta{v}*(v(t-1)+v(t+1)))
                      *(2._dp-(1._dp/cvdifts)*pum1(jl,jk)+zudif(jl,jk))    & !   = - ( u du/dt + v dv/dt ) * 2dt (using central differences)
                      +(pvm1(jl,jk)/cvdifts-zvdif(jl,jk))                  & !   = - Delta{E~} due to diffusion, E~ = TKE divided by mass
                      *(2._dp-(1._dp/cvdifts)*pvm1(jl,jk)+zvdif(jl,jk)))
        zvidis(jl)=zvidis(jl)+                                             & ! Dissipation (over 2 dt) vertically integrated with (negative of) pressure / g
                   zdis(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))/g             !   Equal to contribution to -Delta{TKE} / A, A = area
     END DO
  END DO
#ifdef MESSYTENDENCY
  CALL mtend_add_l (my_handle, mtend_id_u, px = zdudt_2d)
  CALL mtend_add_l (my_handle, mtend_id_v, px = zdvdt_2d)
#endif
!
  DO jl=1,kproma
     pvdis(jl)=zvidis(jl)/ztmst                                              ! -Delta{TKE} / (A * Delta{t})
  END DO
!
!     ------------------------------------------------------------------
!
!*       5.     DIFFUSION IMPLICIT COMPUTATIONS FOR HEAT (S.+L.).
!
! Equations, based on Richtmeyer and Morton (1967), are specified in Schulz et al., 2001
!
  DO jk=1,klev
     DO jl=1,kproma                                                          ! Unnecessary initializing at 0
        ztdif(jl,jk)=0._dp                                                   !   for s (dry static energy)
        zqdif(jl,jk)=0._dp                                                   !   for q (specific humidity)
        zxldif(jl,jk)=0._dp                                                  !   for q_l (cloud water)
        zxidif(jl,jk)=0._dp                                                  !   for q_i (cloud ice)
        zvardif(jl,jk)=0._dp                                                 !   for b-a (width of q distribution)
     END DO
  END DO
!

  DO jt = 1, ktrac
    IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN                        !   for tracer jt, c, if vertical diffusion is active
      DO jk=1,klev
        DO jl=1,kproma
          zxtdif(jl,jk,jt)=0._dp
        END DO
      END DO
    END IF
  END DO

!
!*       5.1     SETTING OF RIGHT HAND SIDES.
!
  DO jk=1,klev
     DO jl=1,kproma
        ztdif(jl,jk)=zcptgz(jl,jk)/cvdifts                                   ! s(t+1)** / computational factor 1.5 = X_(t-1) / alpha (for u) in Schulz et al. (2001)
        zqdif(jl,jk)=pqm1(jl,jk)/cvdifts                                     ! q(t+1)** / computational factor 1.5 = X_(t-1) / alpha (for u) in Schulz et al. (2001)
        zxldif(jl,jk)=pxlm1(jl,jk)/cvdifts                                   ! q_l(t+1)** / computational factor 1.5 = X_(t-1) / alpha (for u) in Schulz et al. (2001)
        zxidif(jl,jk)=pxim1(jl,jk)/cvdifts                                   ! q_i(t+1)** / computational factor 1.5 = X_(t-1) / alpha (for u) in Schulz et al. (2001)
        zvardif(jl,jk)=pxvar(jl,jk)/cvdifts                                  ! b-a(t+1)** / computational factor 1.5 = X_(t-1) / alpha (for u) in Schulz et al. (2001)
     END DO
  END DO
!
#ifdef MESSYTENDENCY
  call  mtend_get_start_l(mtend_id_tracer, v0t = zxtvn)                      ! Actual value of tracers , consistent with the leapfrog operator splitting integration scheme
#endif
  DO jt = 1, ktrac
     IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN                       ! Only calculate if vertical diffusion is active for this tracer
        DO jk=1,klev
           DO jl=1,kproma
#ifndef MESSYTENDENCY
              zxtvn(jl,jk,jt)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst          ! Actual value of tracers , consistent with the leapfrog operator splitting integration scheme
#endif
              zxtdif(jl,jk,jt)=zxtvn(jl,jk,jt)/cvdifts                       ! c(t+1)** / computational factor 1.5 = X_(t-1) / alpha (for u) in Schulz et al. (2001)
           END DO
        END DO
     END IF
  END DO
!
!*       5.2     TOP LAYER ELIMINATION.
!
  DO jl=1,kproma                                                             ! Here k = 1; A_(k+1/2) = zcfh(k), except for b-a: there it is zcfv(k)
     zqdp=1._dp/(paphm1(jl,2)-paphm1(jl,1))                                  ! 1 / Delta{p} around level 1
     zdisc=1._dp/(1._dp+zcfh(jl,1)*zqdp)                                     ! 1 / ( 1 + A_(k+1/2) / Delta{p} ) = 1 / ( 1 + A^_(k+1/2) )
     zebsh(jl,1)=zdisc*(zcfh(jl,1)*zqdp)                                     ! A^_(k+1/2) / ( 1 + A^_(k+1/2) ) = E_(k+1/2) (for k = 1)
     zdisv=1._dp/(1._dp+zcfv(jl,1)*zqdp)                                     ! For b-a: 1 / ( 1 + A^_(k+1/2) )
     zebsv(jl,1)=zdisv*(zcfv(jl,1)*zqdp)                                     ! For b-a: E_(k+1/2) (for k = 1)
     ztdif(jl,1)=zdisc*ztdif(jl,1)                                           ! For s:   X_(t-1) / alpha / ( 1 + A^_(k+1/2) ) = F_(k+1/2)
     zqdif(jl,1)=zdisc*zqdif(jl,1)                                           ! For q:   F_(k+1/2)
     zxldif(jl,1)=zdisc*zxldif(jl,1)                                         ! For q_l: F_(k+1/2)
     zxidif(jl,1)=zdisc*zxidif(jl,1)                                         ! For q_i: F_(k+1/2)
     zvardif(jl,1)=zdisv*zvardif(jl,1)                                       ! For b-a: F_(k+1/2)
  END DO
!
  DO jt = 1, ktrac
     IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN
        DO jl=1,kproma
           zqdp=1._dp/(paphm1(jl,2)-paphm1(jl,1))                            ! 1 / Delta{p} around level 1
           zdisc=1._dp/(1._dp+zcfh(jl,1)*zqdp)                               ! For c: 1 / ( 1 + A^_(k+1/2) )
           zxtdif(jl,1,jt)=zdisc*zxtdif(jl,1,jt)                             ! For c: F_(k+1/2)
        END DO
     END IF
  END DO
!
!*       5.3     ELIMINATION FOR MIDDLE LAYERS.
!
  DO jk=2,klevm1
     DO jl=1,kproma
        zqdp=1._dp/(paphm1(jl,jk+1)-paphm1(jl,jk))                           ! 1 / Delta{p} around level k
        zfac=zcfh(jl,jk-1)*zqdp                                              ! A_(k-1/2) / Delta{p} = C^_(k+1/2)
        zdisc=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,jk-1))+zcfh(jl,jk)*zqdp)     ! 1 / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) )
        zebsh(jl,jk)=zdisc*(zcfh(jl,jk)*zqdp)                                ! A^_(k+1/2) / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) ) = E_(k+1/2)
        zfav=zcfv(jl,jk-1)*zqdp                                              ! For b-a: C^_(k+1/2)
        zdisv=1._dp/(1._dp+zfav*(1._dp-zebsv(jl,jk-1))+zcfv(jl,jk)*zqdp)     ! For b-a: 1 / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) )
        zebsv(jl,jk)=zdisv*(zcfv(jl,jk)*zqdp)                                ! For b-a: E_(k+1/2)
        ztdif(jl,jk)=zdisc*(ztdif(jl,jk)+zfac*ztdif(jl,jk-1))                ! For s:   (X_(t-1) / alpha + C^_(k+1/2) F_(k-1/2)) / (...) = F_(k+1/2)
        zqdif(jl,jk)=zdisc*(zqdif(jl,jk)+zfac*zqdif(jl,jk-1))                ! For q:   F_(k+1/2)
        zxldif(jl,jk)=zdisc*(zxldif(jl,jk)+zfac*zxldif(jl,jk-1))             ! For q_l: F_(k+1/2)
        zxidif(jl,jk)=zdisc*(zxidif(jl,jk)+zfac*zxidif(jl,jk-1))             ! For q_i: F_(k+1/2)
        zvardif(jl,jk)=zdisv*(zvardif(jl,jk)+zfav*zvardif(jl,jk-1))          ! For b-a: F_(k+1/2)
     END DO
  END DO
!
  DO jt = 1, ktrac
     IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN
        DO jk=2,klevm1
           DO jl=1,kproma
              zqdp=1._dp/(paphm1(jl,jk+1)-paphm1(jl,jk))                     ! 1 / Delta{p} around level k
              zfac=zcfh(jl,jk-1)*zqdp                                        ! For c: C^_(k+1/2)
              zdisc=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,jk-1))               & ! For c: 1 / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) )
                                                 +zcfh(jl,jk)*zqdp)
              zxtdif(jl,jk,jt)=zdisc*                                      & ! For c: F_(k+1/2)
                (zxtdif(jl,jk,jt)+zfac*zxtdif(jl,jk-1,jt))
           END DO
        END DO
     END IF
  END DO
!
!        5.4A    EQUIVALENT EVAPOTRANSPIRATION EFFICIENCY COEFF. OVER LAND.
!
  DO jl=1,kproma
     IF (.NOT.loglac_2d(jl,krow)) THEN
       pwlmx(jl)=cwlmax*(1._dp+lai(jl,krow))
     ELSE
       pwlmx(jl)= 0._dp
     END IF
     IF(loland_2d(jl,krow).AND..NOT.loglac_2d(jl,krow)) THEN   ! ub_ak_20190306 ! Correct fractional snow cover and wet skin fraction - only over non-glacier covered land
       zepot=zdens(jl)*zchl(jl)*SQRT(zdu2(jl))*(zqsl(jl)-pqm1(jl,klev))      ! Air density * C_H * |U| * -Delta{q} = LE / L : kg water evaporating / m2 / s
       ! WARNING: the 'corrections' below are not to improve the physically sound snow and liquid water covers, but to prevent that
       ! with the instantaneously calculated LE's more water is removed from the reservoirs in 1 time step than they contain at the
       ! start of that time step. When water from those reservoirs is depleted, standard LE calculated over land applies. This
       ! technical adaptation of the covers yields that result.
       IF(pcvs(jl).GT.0._dp) THEN                                            ! If snow present
         zsnfac=pcvs(jl)*zepot*zdtime/(rhoh2o*(psn(jl)+psnc(jl)))            ! Fractional snow cover * mass of water able to evaporate per m2 in one timestep / available mass of snow per m2
         IF(zsnfac.GT.1._dp) pcvs(jl)=pcvs(jl)/zsnfac                        ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than equivalent snow water content from soil and canopy
       END IF
       IF(pcvw(jl).GT.0._dp) THEN                                            ! If wet skin reservoir present
         zwlfac=(1._dp-pcvs(jl))*zepot*zdtime/(rhoh2o*pwlmx(jl))             ! Fractional non-snow cover * mass of water able to evaporate per m2 in one timestep / available mass in skin reservoir per m2
         IF(zwlfac.GT.1._dp) pcvw(jl)=pcvw(jl)/zwlfac                        ! Modify wet skin fraction if water loss during the time step due to potential evaporation is larger than water content in reservoir
       END IF
     END IF
  END DO
  DO jl=1,kproma
     IF (pws(jl).GT.zplmin*pwsmx(jl)) THEN                                   ! More soil water than wilting point
        zwet(jl)=pcvs(jl)+(1._dp-pcvs(jl))*(pcvw(jl)                       & ! Full evapotranspiration over snow and wet skin fractions; over rest of vegetated surface, the total resistance is equal to r_a + r_stom instead of solely r_a
                         +(1._dp-pcvw(jl))/                                & !     Therefore over vegetation the flux as calculated using r_a is multiplied by r_a/(r_a+r_stom) = 1/(1+(1/r_a)*r_stom); r_a = 1/(C_H |U|)
                          (1._dp+zchl(jl)*SQRT(zdu2(jl))*zwet(jl)))
     ELSE
        zwet(jl)=pcvs(jl)+(1._dp-pcvs(jl))*pcvw(jl)                          ! Infinite resistance for vegetation -> only evapotranspiration over snow and wet skin fractions
     END IF
     lo=zhum(jl).LE.pqm1(jl,klev)/zqsl(jl) .OR. zhum(jl).LT.zepevap
     zcsat(jl)=pcvs(jl)+(1._dp-pcvs(jl))                                   & ! Over bare soil, full evapotranspiration over snow and wet skin fractions; over rest of (bare) surface, it scales with moistness of soil
             *(pcvw(jl)+(1._dp-pcvw(jl))*MERGE(0._dp,zhum(jl),lo))           ! beta * h in Schulz et al. (2001)
     ! Area fraction where air comes into contact with moisture (either snow,
     ! wet skin or a moist soil) beta in Schulz et al. (2001)
     zcair(jl)=pcvs(jl)+(1._dp-pcvs(jl))                                   &
             *(pcvw(jl)+(1._dp-pcvw(jl))*MERGE(0._dp,1._dp,lo))
     lo=pqm1(jl,klev).GT.zqsl(jl)                                            ! Check where q in atmosphere is higher than sat. spec. humidity of land surface
     zcsat(jl)=MERGE(1._dp,zcsat(jl),lo)                                     ! If too much moisture in the air, everywhere there is full evapotranspiration
     zcair(jl)=MERGE(1._dp,zcair(jl),lo)                                     ! If too much moisture in the air, everywhere air comes into contact with moisture
     zcsat(jl)=pvgrat(jl)*zwet(jl)+(1._dp-pvgrat(jl))*zcsat(jl)              ! Total efficiency of evapotranspiration - fraction of q_sl that contributes to q @ surface
     ! Total efficiency for air to come into contact with soil moisture
     ! (not for too dry soil or, if vegetated, partly due to canopy resistance);
     ! where not the case, q @ surface is governed by atmosphere instead of soil
     zcair(jl)=pvgrat(jl)*zwet(jl)+(1._dp-pvgrat(jl))*zcair(jl)
     zcpq(jl) =cpd*(1._dp+vtmpc2*(zcsat(jl)*zqsl(jl)                       & ! Effective specific heat capacity, c_p
                             +(1._dp-zcair(jl))*pqm1(jl,klev)))
     zcptl(jl)=ptslm1(jl)*zcpq(jl)                                           ! Surface dry static energy
  END DO
!
!*       5.4B    BOTTOM LAYER ELIMINATION.
!
  DO jl=1,kproma
     zqdp=1._dp/(paphm1(jl,klevp1)-paphm1(jl,klev))                          ! 1 / Delta{p} around level N=klev
     zfac=zcfh(jl,klevm1)*zqdp                                               ! A_(k-1/2) / Delta{p} = C^_(k+1/2)
     ! A^(N+1/2) = 0 for ice, water and b-a, since those are not quantities are not directly exchanged between air and surface
     ! As a result, E_(k+1/2) = 0 and B^_(k+1/2) = 1 + C^_(k+1/2)
     zdisx=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1)))                       ! 1 / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) )
     zfav=zcfv(jl,klevm1)*zqdp                                               ! For b-a: C^_(k+1/2)
     zdisv=1._dp/(1._dp+zfav*(1._dp-zebsv(jl,klevm1)))                       ! For b-a: 1 / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) )
     zxldif(jl,klev)=zdisx*(zxldif(jl,klev)+zfac*zxldif(jl,klevm1))          ! For q_l: F_(k+1/2) = X*_k/alpha
     zxidif(jl,klev)=zdisx*(zxidif(jl,klev)+zfac*zxidif(jl,klevm1))          ! For q_i: F_(k+1/2) = X*_k/alpha
     zvardif(jl,klev)=zdisv*(zvardif(jl,klev)+zfav*zvardif(jl,klevm1))       ! For b-a: F_(k+1/2) = X*_k/alpha
!
!*  CALCULATION OF THE EN AND FN COEFFICIENTS OF THE RICHTMYER-
!*  MORTON-SCHEME CONCERNING THE EQUATION:
!
!*  XN = EN * XS + FN
!
!*  WITH XN = S_ATM  OR  XN = QATM : ATM. VALUE OF S OR Q
!*  AND  XS = SSURF  OR  XS = QSAT : SURFACE VALUE OF S OR SAT. SPEC.
!*                                   HUM. AT THE SURFACE
!
!    land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zdiscl=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1))                     & ! For s: 1 / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) )
                                      +zcfhl(jl)*zqdp)                       ! Note that for the surface (A_(N+1/2)) the eddy viscosity is calculated as C_H |U| instead of K / Delta{z}
     zdisql=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1))                     & ! For q: 1 / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) ), adapted A^(k+1/2) (in B^(k+1/2)) to only have exchange over limited area
                                      +zcair(jl)*zcfhl(jl)*zqdp)             ! Accomplished by using zcair, i.e. area fraction of land where air comes into contact with moisture
     zetnl(jl)=zdiscl*zcfhl(jl)*zqdp                                         ! For s: E_(k+1/2)
     zftnl(jl)=zdiscl*(ztdif(jl,klev)+zfac*ztdif(jl,klevm1))*cvdifts         ! For s: F_(k+1/2) * alpha
     zeqnl(jl)=zdisql*zcsat(jl)*zcfhl(jl)*zqdp                               ! For q: E_(k+1/2) * fraction of q_sat that determines q at surface (0 where air does not come into contact with soil moisture), since q_sat is used as X*_surf
     zfqnl(jl)=zdisql*(zqdif(jl,klev)+zfac*zqdif(jl,klevm1))*cvdifts         ! For q: F_(k+1/2) * alpha
!
!    water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zdiscw=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1))+zcfhw(jl)*zqdp)       ! 1 / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) )
     zetnw(jl)=zdiscw*zcfhw(jl)*zqdp                                         ! For s: E_(k+1/2)
     zftnw(jl)=zdiscw*(ztdif(jl,klev)+zfac*ztdif(jl,klevm1))*cvdifts         ! For s: F_(k+1/2) * alpha
     zeqnw(jl)=zetnw(jl)                                                     ! For q: E_(k+1/2), identical to for s
     zfqnw(jl)=zdiscw*(zqdif(jl,klev)+zfac*zqdif(jl,klevm1))*cvdifts         ! For q: F_(k+1/2) * alpha
!
!    ice   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zdisci=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1))+zcfhi(jl)*zqdp)       ! 1 / ( B^_(k+1/2) - C^_(k+1/2) E_(k-1/2) )
     zetni(jl)=zdisci*zcfhi(jl)*zqdp                                         ! For s: E_(k+1/2)
     zftni(jl)=zdisci*(ztdif(jl,klev)+zfac*ztdif(jl,klevm1))*cvdifts         ! For s: F_(k+1/2) * alpha
     zeqni(jl)=zetni(jl)                                                     ! For q: E_(k+1/2), identical to for s
     zfqni(jl)=zdisci*(zqdif(jl,klev)+zfac*zqdif(jl,klevm1))*cvdifts         ! For q: F_(k+1/2) * alpha
  END DO
!
  ! mz_rs_20040329+
  DO jt = 1, ktrac
    IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN                        ! Only calculate if vertical diffusion is active for this tracer
      DO jl=1,kproma
        zqdp=1._dp/(paphm1(jl,klevp1)-paphm1(jl,klev))                       ! 1 / Delta{p} around level N=klev
        zfac=zcfh(jl,klevm1)*zqdp                                            ! A_(k-1/2) / Delta{p} = C^_(k+1/2)
!
! Using prescribed fluxes, the Richtmyer-Morton scheme slightly changes for level N
! X*_N / alpha = E_(N+1/2) X*_surf / alpha + F_(N+1/2) can be rewritten to
! X*_N / alpha * ( 1 + A^_(N+1/2) + C^_(N+1/2) ( 1 - E_(N-1/2) ) ) =
!                A^_(N+1/2) X*_surf / alpha + X_N / alpha + C^_(N+1/2) F_(N-1/2)
! X*_N / alpha = A^_(N+1/2) / ( 1 + C^_(N+1/2) ( 1 - E_(N-1/2) ) )
!                ( X*_surf - X*_N ) / alpha + ( X_N / alpha + C^_(N+1/2) F_(N-1/2) )
!                / ( 1 + C^_(N+1/2) ( 1 - E_(N-1/2) ) )
! X*_N / alpha = E~_(N+1/2) ( X*_surf - X*_N ) / alpha + F~_(N+1/2)
!
! zxtems = rho C_H |U| ( X*_surf - X*_N ) = emission flux, corrected for dry deposition
!
        zdisxt=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1)))                   ! For c: 1 / ( 1 + C^_(N+1/2) ( 1 - E_(N-1/2) ) )
        zxtdif(jl,klev,jt)=zdisxt*(ztmst*g*zqdp*zxtems(jl,jt)+             & ! For c: X*_N / alpha
                                   zxtdif(jl,klev,jt)+                     &
                                   zfac*zxtdif(jl,klevm1,jt))
      END DO
    END IF
  END DO
  ! mz_rs_20040329-
!
  DO jl = 1,kproma
     znetr(jl) =zsrfll(jl)+pemter(jl,klevp1)                                 ! R_net = SW_net + LW_net
     zcdrag(jl)=zdens(jl) * zchl(jl) * SQRT(zdu2(jl))                        ! rho * C_H * |U|
  END DO
!
!    Land surface temperature (cp*Ts=zcptlnew) and
!    saturation specific humidity (qsat=zqslnew)
!
  IF (.NOT.lstart) THEN
     CALL vertex_surftemp(kproma, cvdifts*ztmst,                           & ! Amount of columns (vectorization), alpha * leap frog time step
        cemiss, stbo, zcpq, cpd*vtmpc2, alv, als,                          & ! Emissivity, Stephan-Boltzmann constant, effective cp, cpv - cpd, Lv, Ls
        zftnl, zetnl, zfqnl, zeqnl,                                        & ! alpha * F_(N+1/2) for s, E_(N+1/2) for s, alpha * F_(N+1/2) for q, E_(N+1/2) * csat for q
        zcptl, zqsl, zdqsl,                                                & ! Old values at the surface: s, q_s, d{q_s}/d{T}
        znetr, pgrndhflx,                                                  & ! R_net, G
        zcdrag, zcair, zcsat, pcvs, pgrndcapc,                             & ! rho * C_H * |U|, cair (= beta), csat (= beta * h), snow cover, soil heat capacity (Cs)
        loland_2d(1:kproma,krow),                                          & ! Logical land mask
        zcptlnew, zqslnew)                                                   ! Output: updated s and q_s as X* (= alpha X_(t+1) + (1-alpha) X_(t-1) ) values
  ELSE
     DO jl = 1,kproma
        zcptlnew(jl)=zcptl(jl)
        zqslnew(jl)=zqsl(jl)
     END DO
  END IF
!
!   Land surface temperature and 'zcptlnew' correction for snowmelt
!
  DO jl = 1,kproma
    zcpt=zcptlnew(jl)/cvdifts+(1._dp-(1._dp/cvdifts))*zcptl(jl)              ! Surface dry static energy at t+1, s_(S,t+1) = s*_S/alpha + (1 - 1/alpha) * s_(S,t-1)
    IF (loland_2d(jl,krow)) THEN
      ptsl(jl)=zcpt/zcpq(jl)                                                 ! Over land: T_S = s / c_p
    ELSE
      ptsl(jl)=tmelt                                                         ! Elsewhere: T_S is snow melt temperature
    END IF
    IF (psn(jl).GT.csncri.OR.loglac_2d(jl,krow)) THEN                                ! If snow depth > 5.85 mm or over glacier
       zcptlcorr=MIN(ptsl(jl),tmelt)*zcpq(jl)                                !   @ t+1: T_S, used to calculate s_S, has snow melt temperature as upper limit
       zcptlnew(jl)=cvdifts*zcptlcorr+(1._dp-cvdifts)*zcptl(jl)              !   s*_S = alpha s_(S,t+1) + (1 - alpha) s_(S,t-1)
    ENDIF
!
!   New land surface temperature for sensible heat flux and *radheat*
!
    ptslnew(jl)=zcptlnew(jl)/zcpq(jl)                                        ! T*_S as (implicitly) used in the surface energy balance
  END DO
!
!*  CALCULATION OF SKLEV AND QKLEV USING THE NEW SURFACE VALUES
!*  ZSNEW AND ZQSNEW WHICH WERE CALCULATED IN SUBROUTINE SURFTEMP
!
  DO jl = 1,kproma
!
!    land   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     ztklevl=zetnl(jl)*zcptlnew(jl)+zftnl(jl)                                ! Over land: s*_N = E_(s,N+1/2) s*_surf + alpha F_(s,N+1/2)
     zqklevl=zeqnl(jl)*zqslnew(jl)+zfqnl(jl)                                 ! Over land: q*_N = E_(q,N+1/2) q*_surf + alpha F_(q,N+1/2)
!
!    water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     ztklevw=zetnw(jl)*zcptw(jl)+zftnw(jl)                                   ! Over water: s*_N = E_(s,N+1/2) s*_surf + alpha F_(s,N+1/2)
     zqklevw=zeqnw(jl)*zqsw(jl)+zfqnw(jl)                                    ! Over water: q*_N = E_(q,N+1/2) q*_surf + alpha F_(q,N+1/2)
!
!    ice   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     ztklevi=zetni(jl)*zcpti(jl)+zftni(jl)                                   ! Over ice: s*_N = E_(s,N+1/2) s*_surf + alpha F_(s,N+1/2)
     zqklevi=zeqni(jl)*zqsi(jl)+zfqni(jl)                                    ! Over ice: q*_N = E_(q,N+1/2) q*_surf + alpha F_(q,N+1/2)
!
!    Grid-mean dry static energy and specific humidity at the
!    'blending height' (here: lowest model level 'klev').
!
     zgtl=pfrl(jl)*zcfhl(jl)                                                 ! fraction land * Computational factor 1.5 * leap frog time step * g * air density * C_H * |U| over land
     zgtw=pfrw(jl)*zcfhw(jl)                                                 ! fraction water * Computational factor 1.5 * leap frog time step * g * air density * C_H * |U| over water
     zgti=pfri(jl)*zcfhi(jl)                                                 ! fraction ice * Computational factor 1.5 * leap frog time step * g * air density * C_H * |U| over ice
     zgtsum=(zgtl+zgtw+zgti)*cvdifts                                         ! alpha * sum of the factors above
     IF (pfrl(jl).LT.1._dp) THEN                                             ! if not the entire grid cell is land, account for limited moisture exchange over part of the area (zcair = beta)
        zgql=zgtl*zcair(jl)                                                  ! fraction land * Computational factor 1.5 * leap frog time step * g * air density * C_H * |U| over land
     ELSE
        zgql=zgtl                                                            ! fraction land * Computational factor 1.5 * leap frog time step * g * air density * C_H * |U| over land
     ENDIF
     zgqw=zgtw                                                               ! fraction water * Computational factor 1.5 * leap frog time step * g * air density * C_H * |U| over water
     zgqi=zgti                                                               ! fraction ice * Computational factor 1.5 * leap frog time step * g * air density * C_H * |U| over ice
     zgqsum=(zgql+zgqw+zgqi)*cvdifts                                         ! alpha * sum of the factors above
     ztdif(jl,klev)=(zgtl*ztklevl+zgtw*ztklevw+zgti*ztklevi)/zgtsum          ! Weighted (by cover fraction and exchange efficiency, C_H * |U|) s*_N / alpha
     zqdif(jl,klev)=(zgql*zqklevl+zgqw*zqklevw+zgqi*zqklevi)/zgqsum          ! Weighted (by cover fraction and exchange efficiency, C_H * |U|) q*_N / alpha
  END DO
!
!*         5.5     BACK-SUBSTITUTION.
!
  DO jk=klevm1,1,-1
     DO jl=1,kproma
        ztdif(jl,jk)=ztdif(jl,jk)+zebsh(jl,jk)*ztdif(jl,jk+1)                ! For s:   X*_k/alpha = F_(k+1/2) + E_(k+1/2) * X*_(k+1)/alpha
        zqdif(jl,jk)=zqdif(jl,jk)+zebsh(jl,jk)*zqdif(jl,jk+1)                ! For q:   X*_k/alpha = F_(k+1/2) + E_(k+1/2) * X*_(k+1)/alpha
        zxldif(jl,jk)=zxldif(jl,jk)+zebsh(jl,jk)*zxldif(jl,jk+1)             ! For q_l: X*_k/alpha = F_(k+1/2) + E_(k+1/2) * X*_(k+1)/alpha
        zxidif(jl,jk)=zxidif(jl,jk)+zebsh(jl,jk)*zxidif(jl,jk+1)             ! For q_i: X*_k/alpha = F_(k+1/2) + E_(k+1/2) * X*_(k+1)/alpha
        zvardif(jl,jk)=zvardif(jl,jk)+zebsv(jl,jk)*zvardif(jl,jk+1)          ! For b-a: X*_k/alpha = F_(k+1/2) + E_(k+1/2) * X*_(k+1)/alpha
     END DO
  END DO
!
  DO jt = 1, ktrac
    IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN
      DO jk=klevm1,1,-1
        DO jl=1,kproma
          zxtdif(jl,jk,jt)=zxtdif(jl,jk,jt)+zebsh(jl,jk)*zxtdif(jl,jk+1,jt)  ! For c: X*_k/alpha = F_(k+1/2) + E_(k+1/2) * X*_(k+1)/alpha
        END DO
      END DO
    END IF
  END DO
!
!*         5.6     INCREMENTATION OF T AND Q TENDENCIES.
!
!
  DO jk=1,klev
     DO jl=1,kproma
        zqdif(jl,jk)=zqdif(jl,jk)+(1._dp-(1._dp/cvdifts))*pqm1(jl,jk)        ! q(t+1) = q*/alpha + (1 - 1/alpha) * q(t-1)
        zdqdt=(zqdif(jl,jk)-pqm1(jl,jk))/ztmst                               ! Tendency in q = ( q(t+1) - q(t-1) ) / 2dt
#ifndef MESSYTENDENCY
        pqte(jl,jk)=pqte(jl,jk)+zdqdt
#else
        zdqdt_2d(jl,jk)=zdqdt
#endif
        ztdif(jl,jk)=ztdif(jl,jk)+(1._dp-(1._dp/cvdifts))*zcptgz(jl,jk)      ! s(t+1) = s*/alpha + (1 - 1/alpha) * s(t-1)
        zdtdt=((ztdif(jl,jk)+zdis(jl,jk)-pgeom1(jl,jk))                    & ! Tendency in T = ( ( ( s(t+1) - g z ) + Diss ) / c_p - T(t-1) ) / 2dt
              /(cpd*(1._dp+vtmpc2*zqdif(jl,jk)))-ptm1(jl,jk))/ztmst          !   Diss = heat added due to dissipation of TKE
#ifndef MESSYTENDENCY
        ptte(jl,jk)=ptte(jl,jk)+zdtdt
#else
        zdtdt_2d(jl,jk)=zdtdt
#endif
        zxldif(jl,jk)=zxldif(jl,jk)+(1._dp-(1._dp/cvdifts))*pxlm1(jl,jk)     ! q_l(t+1) = q_l*/alpha + (1 - 1/alpha) * q_l(t-1)
        zxidif(jl,jk)=zxidif(jl,jk)+(1._dp-(1._dp/cvdifts))*pxim1(jl,jk)     ! q_i(t+1) = q_i*/alpha + (1 - 1/alpha) * q_i(t-1)
        zdxlmdt=(zxldif(jl,jk)-pxlm1(jl,jk))/ztmst                           ! Tendency in q_l = ( q_l(t+1) - q_l(t-1) ) / 2dt
        zdximdt=(zxidif(jl,jk)-pxim1(jl,jk))/ztmst                           ! Tendency in q_i = ( q_i(t+1) - q_i(t-1) ) / 2dt
#ifndef MESSYTENDENCY
        pxlte(jl,jk)=pxlte(jl,jk)+zdxlmdt
        pxite(jl,jk)=pxite(jl,jk)+zdximdt
#else
        zdxlmdt_2d(jl,jk)=zdxlmdt
        zdximdt_2d(jl,jk)=zdximdt
#endif
        pxvar(jl,jk)=zvardif(jl,jk)+(1._dp-(1._dp/cvdifts))*pxvar(jl,jk)     ! {b-a}(t+1) = {b-a}*/alpha + (1 - 1/alpha) * {b-a}(t-1)
        pvdiffp(jl,jk)=zdqdt+zdxlmdt+zdximdt !store for production           ! Delta{q_tot}/Delta{t}
     END DO
  END DO
#ifdef MESSYTENDENCY
  CALL mtend_add_l (my_handle, mtend_id_t, px = zdtdt_2d)
  CALL mtend_add_l (my_handle, mtend_id_q, px = zdqdt_2d)
  CALL mtend_add_l (my_handle, mtend_id_xl, px = zdxlmdt_2d)
  CALL mtend_add_l (my_handle, mtend_id_xi, px = zdximdt_2d)
#endif
!
  DO jt = 1, ktrac
    IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN
      DO jk=1,klev
        DO jl=1,kproma
          zxtdif(jl,jk,jt)=zxtdif(jl,jk,jt) +                              & ! c(t+1) = c*/alpha + (1 - 1/alpha) * c(t-1)
            (1._dp-(1._dp/cvdifts))*zxtvn(jl,jk,jt)
          zdxtdt=(zxtdif(jl,jk,jt)-zxtvn(jl,jk,jt))/ztmst                    ! Tendency in c = ( c(t+1) - c(t-1) ) / 2dt
          ! discard negative values
          IF (zxtvn(jl,jk,jt) > 0.0_DP) THEN                                 ! Ensuring that tendencies cannot lead to (more) negative concentrations
             zdxtdt = MAX(-zxtvn(jl,jk,jt)/ztmst, zdxtdt)
          ELSE
             zdxtdt = MAX(0.0_DP, zdxtdt)
          END IF
#ifndef MESSYTENDENCY
          pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
#else
          zdxtdt_3d(jl,jk,jt)=zdxtdt
#endif
        END DO
      END DO
    END IF
  END DO

#ifdef MESSYTENDENCY
  IF (ktrac > 0) &
       CALL mtend_add_l(my_handle, mtend_id_tracer, pxt=zdxtdt_3d)
#endif
!
! back out moisture flux
!
  DO jl=1,kproma
     zqflux(jl,1)=0._dp                                                      ! No flux at top of the domain
     zvarpr(jl,1)=0._dp                                                      ! No q variance production at top of the domain
     ztvlan=ptslm1(jl)*(1._dp+vtmpc1*zhsoil(jl)*zqsl(jl))                    ! T_v at surface-atmosphere interface over land
     ztvsea=ptsw(jl)*(1._dp+vtmpc1*zqsw(jl))                                 ! T_v at surface-atmosphere interface over water
     ztvice=ptsi(jl)*(1._dp+vtmpc1*zqsi(jl))                                 ! T_v at surface-atmosphere interface over ice
     ztvh=pfrl(jl)*ztvlan+pfrw(jl)*ztvsea+pfri(jl)*ztvice                    ! Average T_v
     zrho(jl,klevp1)=paphm1(jl,klevp1)/(rd*ztvh)                             ! Air density directly at the surface
     zqsurf=pfrl(jl)*zqsl(jl)*zcsat(jl)+pfrw(jl)*zqsw(jl)+pfri(jl)*zqsi(jl)  ! Average q at the surface
     zdqtot=(pqm1(jl,klev)+zx(jl,klev))-zqsurf                               ! Difference in total atmospheric water content between surface and center lowest grid cell
     zqshear(jl,klev)=zdqtot*g/pgeom1(jl,klev)                               ! Gradient in total atmospheric water content between surface and center lowest grid cell
  ENDDO !jl
  DO jk=2,klevp1
     DO jl=1,kproma
        IF (jk<klevp1) THEN
           ztvh=(ptm1(jl,jk)*(1._dp+vtmpc1*pqm1(jl,jk)-zx(jl,jk))          & ! Interpolation of T_v at interface levels
                +ptm1(jl,jk-1)*(1._dp+vtmpc1*pqm1(jl,jk-1)                 &
                 -zx(jl,jk-1)))/2._dp
           zrho(jl,jk)=paphm1(jl,jk)/(rd*ztvh)                               ! Air density at interface levels
        ENDIF
        zrhodz=-(paphm1(jl,jk)-paphm1(jl,jk-1))/g                            ! - rho Delta{z} between interface levels
        zqflux(jl,jk)=zrhodz*pvdiffp(jl,jk-1)+zqflux(jl,jk-1)                ! - {rho w'q'}(k+1/2) = - {rho w'q'}(k-1/2) - {rho Delta{z}}(k) * d{q}/d{t}
        zvarpr(jl,jk)=zqshear(jl,jk-1)*zqflux(jl,jk)/zrho(jl,jk)             ! - {rho w'q'}(k+1/2) / rho(k+1/2)  * Delta{q}/Delta{z} = q-variance production by vertical gradient at interface levels
     ENDDO !jl
  ENDDO !jk
  DO jk=1,klev
     DO jl=1,kproma
        pvdiffp(jl,jk)=(zvarpr(jl,jk)+zvarpr(jl,jk+1))/2._dp                 ! Interpolated q-variance production by vertical gradient in grid center
        IF(jk.GE.ihpbl(jl)) THEN
           zlam  = clam                                                      ! Asymptotic mixing length, lambda, within the boundary layer - 150 m
        ELSE
           zhexp = EXP(1._dp-pgeom1(jl,jk)/pgeom1(jl,ihpbl(jl)))             ! Exponential factor in expression for asymptotic mixing length, lambda, above the boundary layer
           zlam  = 1._dp+(clam-1._dp)*zhexp                                  ! lambda(z>zi) = 1 + 149 * exp(1 - z/zi) m
        END IF
        z2geomf  = (ckap/g)*pgeom1(jl,jk)                                    ! kappa z at grid center
        zmix     = 1._dp / ( (1._dp/z2geomf) + (1._dp/zlam) )                ! mixing length - l = 1 / ( 1/(kappa*z) + 1/lambda )
        IF(jk.EQ.1) THEN
           ztkesq=SQRT(MAX(ztkemin,ptkem1(jl,1)))                            ! sqrt(E(t))  (Note that this take place after applying the Robert-Asselin filter for E)
        ELSE
           ztkesq=SQRT(MAX(ztkemin,0.5_dp*(ptkem1(jl,jk-1)+ptkem1(jl,jk))))  ! Square root of the interpolated E(t) at grid center
        END IF
        pvmixtau(jl,jk)=ztkesq/(zmix/zsmn**3)                                ! Dissipation time scale = E / d{E}_diss./d{t} = E / (E^(3/2) / Lambda_1 ) =  sqrt(E) / (Lambda_1) , Lambda_1 = l / S_N_M^3
     ENDDO !jl
  ENDDO !jk
!
!*       5.8     Surface fluxes of heat and moisture
!
  DO jl=1,kproma
!
     zcoefl=cvdifts*zdens(jl)*zchl(jl)*SQRT(zdu2(jl))                        ! alpha rho C_H |U| over land
     zcoefw=cvdifts*zdens(jl)*zchw(jl)*SQRT(zdu2oc(jl))                      ! alpha rho C_H |U| over water
     zcoefi=cvdifts*zdens(jl)*zchi(jl)*SQRT(zdu2oc(jl))                      ! alpha rho C_H |U| over ice
!
!* Moisture fluxes
!
     zqnlev=zqdif(jl,klev)-(1._dp-(1._dp/cvdifts))*pqm1(jl,klev)             ! q*/alpha
!
!    land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zzqs=zqslnew(jl)/cvdifts                                                ! Over land: q_(sat,surf)*/alpha
     zqhfll(jl)=zcoefl*(zcair(jl)*zqnlev-zcsat(jl)*zzqs)                     ! Over land: LE (downwards) / L = rho C_H |U| beta (q* - h * q_(sat,surf)*)
!
!    water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zqhflw(jl)=zcoefw*(zqnlev-zqsw(jl)/cvdifts)                             ! Over water: LE (downwards) / L = rho C_H |U| (q* - q_(sat,surf)*) ; q_(sat,surf)* = q_(sat,surf)

!    ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     zqhfli(jl)=zcoefi*(zqnlev-zqsi(jl)/cvdifts)                             ! Over ice: LE (downwards) / L = rho C_H |U| (q* - q_(sat,surf)*) ; q_(sat,surf)* = q_(sat,surf)

!    Area mean moisture flux (=evaporation)

     pqhfla(jl)=pfrl(jl)*zqhfll(jl)+pfrw(jl)*zqhflw(jl)+pfri(jl)*zqhfli(jl)  ! Averaged LE (downwards) / L
!
!* Sensible heat fluxes

     ztnlev=ztdif(jl,klev)-(1._dp-(1._dp/cvdifts))*zcptgz(jl,klev)           ! s*/alpha
!
!    land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zqnew=pqm1(jl,klev)+zcair(jl)*(zqnlev*cvdifts-pqm1(jl,klev))            ! q(t-1) + beta ( q* - q(t-1) )
     zzcpts=(ptslnew(jl)/cvdifts)*cpd*(1._dp+vtmpc2*zqnew)                   ! s_surf*/alpha
     zthfll(jl)=zcoefl*(ztnlev-zzcpts)                                       ! Over land: SH (downwards) = rho C_H |U| (s* - s_surf) ; probably no correction applied because cp* T* is already used
!
!    water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zthflw(jl)=zcoefw*(ztnlev-zcptw(jl)/cvdifts)                            ! Over water: H (downwards) = rho C_H |U| (s* - s_surf)
     zthflw(jl)=zthflw(jl)-ptsw(jl)*cpd*vtmpc2*zqhflw(jl)                    !   Correction: SH (downwards) = H (downwards) - T_surf (cpv - cpd) LE (downwards) / L
!
!    ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zthfli(jl)=zcoefi*(ztnlev-zcpti(jl)/cvdifts)                            ! Over ice: H (downwards) = rho C_H |U| (s* - s_surf)
     zthfli(jl)=zthfli(jl)-ptsi(jl)*cpd*vtmpc2*zqhfli(jl)                    !   Correction: SH (downwards) = H (downwards) - T_surf (cpv - cpd) LE (downwards) / L
!
! Sensible heat flux and evaporation
!

     pevapl(jl)=zqhfll(jl)                                                   ! Storing LE / L over land in channel
     pevapw(jl)=zqhflw(jl)                                                   ! Storing LE / L over water in channel
     pevapi(jl)=zqhfli(jl)                                                   ! Storing LE / L over ice in channel

     pahfsl(jl)=zthfll(jl)                                                   ! Storing SH over land in channel
     pahfsw(jl)=zthflw(jl)                                                   ! Storing SH over water in channel
     pahfsi(jl)=zthfli(jl)                                                   ! Storing SH over ice in channel

     pahfs(jl)=pfrl(jl)*zthfll(jl)+pfrw(jl)*zthflw(jl)+pfri(jl)*zthfli(jl)   ! Weighted (by cover fraction) SH [W/m2]
     pevap(jl)=pfrl(jl)*zqhfll(jl)+pfrw(jl)*zqhflw(jl)+pfri(jl)*zqhfli(jl)   ! Weighted (by cover fraction) LE / L [kg/m2/s]
!
!* Potential evaporation over land (for snow and skin reservoir)
!
     pevapot(jl)=zcoefl*(zqnlev-zzqs)                                        ! PET / L = rho C_H |U| (q* - q_(sat,surf)*)
!

     ! CALCULATION OF SURFACE KINEMATIC HEAT FLUX [Km/s]
     zsenkf_2d(jl,krow)= pahfs(jl) /(zdens(jl)*cpd)                          ! Kinematic heat flux (downward)

!  Latent heat fluxes
!
!    land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     zqvhfl=zqhfll(jl)-pcvs(jl)*pevapot(jl)                                  ! LE/L over land without snow sublimation
     pahfll(jl)=alv*zqvhfl+als*pcvs(jl)*pevapot(jl)                          ! LE over land with correct L applied
!
!    water
!
     pahflw(jl)=alv*zqhflw(jl)                                               ! LE over water (evaporation)
!
!    ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     pahfli(jl)=als*zqhfli(jl)                                               ! LE over ice (sublimation)
!
     pahfl(jl)=pfrl(jl)*pahfll(jl)+pfrw(jl)*pahflw(jl)+pfri(jl)*pahfli(jl)   ! Weighted (by cover fraction) LE [W/m2]

! op_pj_20181105+
     pahfslac(jl)=pfrl(jl)*pahfsl(jl)
     pahfswac(jl)=pfrw(jl)*pahfsw(jl)
     pahfsiac(jl)=pfri(jl)*pahfsi(jl)
!
     pahfllac(jl)=pfrl(jl)*pahfll(jl)
     pahflwac(jl)=pfrw(jl)*pahflw(jl)
     pahfliac(jl)=pfri(jl)*pahfli(jl)
!
     pevaplac(jl)=pfrl(jl)*pevapl(jl)
     pevapwac(jl)=pfrw(jl)*pevapw(jl)
     pevapiac(jl)=pfri(jl)*pevapi(jl)

     ! CALCULATION OF SURFACE KINEMATIC MOISTURE FLUX [kg/kg m/s]
     zlatkf_2d(jl,krow)= pahfl(jl) / (zdens(jl)*alv)                         ! Kinematic moisture flux (downward) = ( LE + snow sublimation ) / ( rho L )
  END DO
!
!        5.9.1   COMPUTE T2M, T2M_MAX T2M_MIN (based on conditions before diffusion)
!
  DO jl=1,kproma
!
!    land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     lo1=zril(jl).GT.0._dp                                                   ! Bulk Richardson number > 0 over land?
     zrat=zhtq/pgeom1(jl,klev)                                               ! 2 m / z - ratio of 2 m vs height of grid center
     zcbn=SQRT(LOG(1._dp+2._dp/paz0l(jl))*LOG(1._dp+2._dp/paz0hl(jl)))       ! sqrt( ln(1 + 2m/z0m) ln(1 + 2m/z0h) ) = kappa / sqrt(C_N_H(2m))
     zcbs=-(zbhnl(jl)-zbhl(jl))*zrat                                         ! For stable: sqrt(C_M(2m)) kappa / C_H(2m) - kappa / sqrt(C_N_H(2m)) is assumed to be: ( sqrt(C_M) kappa / C_H - kappa / sqrt(C_N_H) ) 2 m / z
     zcbu=-LOG(1._dp+(EXP(zbhnl(jl)-zbhl(jl))-1._dp)*zrat)                   ! For unstable: using X = sqrt(C_M) kappa / C_H - kappa / sqrt(C_N_H), X(2m) is assumed to be: - ln( (exp(-X) - 1) * 2 m / z + 1)
     zred=(zcbn+MERGE(zcbs,zcbu,lo1))/zbhl(jl)                               ! { sqrt(C_M(2m)) kappa / C_H(2m) } / { sqrt(C_M) kappa / C_H } = ( s(2m) - s_S ) / ( s_a - s_S )
     zh2m=zcptl(jl)+zred*(zcptgz(jl,klev)-zcptl(jl))                         ! s(2m) = s_S + ( s(2m) - s_S ) / ( s_a - s_S ) * ( s_a - s_S )
     zt2l=(zh2m-zhtq)/(cpd*(1._dp+vtmpc2*pqm1(jl,klev)))                     ! T(2m) = ( s(2m) - g * 2 m ) / cp; cp approximated using grid cell averaged q
!
!    water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     lo1=zriw(jl).GT.0._dp                                                   ! Bulk Richardson number > 0 over water?
     zcbn=LOG(1._dp+(EXP(zbnw(jl))-1._dp)*zrat)                              ! ln(1 + 2 m / z0m) = kappa / sqrt(C_N_M(2m))
     zcbs=-(zbnw(jl)-zbhw(jl))*zrat                                          ! For stable: sqrt(C_M(2m)) kappa / C_H(2m) - kappa / sqrt(C_N_M(2m)) is assumed to be: ( sqrt(C_M) kappa / C_H - kappa / sqrt(C_N_M) ) 2 m / z
     zcbu=-LOG(1._dp+(EXP(zbnw(jl)-zbhw(jl))-1._dp)*zrat)                    ! For unstable: using X = sqrt(C_M) kappa / C_H - kappa / sqrt(C_N_M), X(2m) is assumed to be: - ln( (exp(-X) - 1) * 2 m / z + 1)
     zred=(zcbn+MERGE(zcbs,zcbu,lo1))/zbhw(jl)                               ! { sqrt(C_M(2m)) kappa / C_H(2m) } / { sqrt(C_M) kappa / C_H } = ( s(2m) - s_S ) / ( s_a - s_S )
     zh2m=zcptw(jl)+zred*(zcptgz(jl,klev)-zcptw(jl))                         ! s(2m) = s_S + ( s(2m) - s_S ) / ( s_a - s_S ) * ( s_a - s_S )
     zt2w=(zh2m-zhtq)/(cpd*(1._dp+vtmpc2*pqm1(jl,klev)))                     ! T(2m) = ( s(2m) - g * 2 m ) / cp; cp approximated using grid cell averaged q
!
!    ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     lo1=zrii(jl).GT.0._dp                                                   ! Bulk Richardson number > 0 over ice?
     zcbn=LOG(1._dp+(EXP(zbni(jl))-1._dp)*zrat)                              ! 2 m / z - ratio of 2 m vs height of grid center
     zcbs=-(zbni(jl)-zbhi(jl))*zrat                                          ! For stable: sqrt(C_M(2m)) kappa / C_H(2m) - kappa / sqrt(C_N(2m)) is assumed to be: ( sqrt(C_M) kappa / C_H - kappa / sqrt(C_N) ) 2 m / z
     zcbu=-LOG(1._dp+(EXP(zbni(jl)-zbhi(jl))-1._dp)*zrat)                    ! For unstable: using X = sqrt(C_M) kappa / C_H - kappa / sqrt(C_N), X(2m) is assumed to be: - ln( (exp(-X) - 1) * 2 m / z + 1)
     zred=(zcbn+MERGE(zcbs,zcbu,lo1))/zbhi(jl)                               ! { sqrt(C_M(2m)) kappa / C_H(2m) } / { sqrt(C_M) kappa / C_H } = ( s(2m) - s_S ) / ( s_a - s_S )
     zh2m=zcpti(jl)+zred*(zcptgz(jl,klev)-zcpti(jl))                         ! s(2m) = s_S + ( s(2m) - s_S ) / ( s_a - s_S ) * ( s_a - s_S )
     zt2i=(zh2m-zhtq)/(cpd*(1._dp+vtmpc2*pqm1(jl,klev)))                     ! T(2m) = ( s(2m) - g * 2 m ) / cp; cp approximated using grid cell averaged q
     ptemp2(jl)=pfrl(jl)*zt2l+pfrw(jl)*zt2w+pfri(jl)*zt2i                    ! Weighted T(2m)
     t2l_1d(jl)=zt2l                                                         ! Storing in channel objects
     t2w_1d(jl)=zt2w                                                         ! Storing in channel objects
     t2i_1d(jl)=zt2i                                                         ! Storing in channel objects
!
!        5.9.2   2M DEW POINT
!
     it = NINT(ptm1(jl,klev)*1000._dp)                                       ! Integer temperature in mK for lookup tables
     IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.               ! T should lie between 50 and 400 K
     it = MAX(MIN(it,jptlucu2),jptlucu1)
     zqs1=tlucua(it)/papm1(jl,klev)                                          ! Rd/Rv * es/p
     zqs1=zqs1/(1._dp-vtmpc1*zqs1)                                           ! Saturation specific humidity in lowest grid cell
     zrh2m=MAX(zephum,pqm1(jl,klev)/zqs1)                                    ! Rel. humidity in lowest grid cell; assumed to be constant throughout the cell
     rh_2m(jl,krow) = MIN(1._dp,zrh2m)                                       ! Set rel. humidity channel object, maximum of 1
!
!    Dewpoint temperature calculated from vapour pressure according to Murray (1967)
!
!    land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     lo=zt2l.GT.tmelt                                                        ! Phase of water at 2 m: Ice or liquid water?
     zcvm3=MERGE(c3les,c3ies,lo)                                             ! a in Murray's Eq. (6)
     zcvm4=MERGE(c4les,c4ies,lo)                                             ! b in Murray's Eq. (6)
     zaph2m=paphm1(jl,klevp1)*(1._dp-zhtq/(rd*zt2l                         & ! P(2m) ~ P_S - rho(2m) g 2m
                        *(1._dp+vtmpc1*pqm1(jl,klev)-zx(jl,klev))))          !   rho(2m) ~ P_S / ( Rd T (1 + 0.61 q - q_i) )
     it = NINT(zt2l*1000._dp)                                                ! Integer temperature in mK for lookup tables
     IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.               ! T should lie between 50 and 400 K
     it = MAX(MIN(it,jptlucu2),jptlucu1)
     zqs2=tlucua(it)/zaph2m                                                  ! Rd/Rv * es/p at 2 m height
     zqs2=zqs2/(1._dp-vtmpc1*zqs2)                                           ! Saturation specific humidity at 2 m height
     zq2m=zrh2m*zqs2                                                         ! Specific humidity at 2 m height
     zfrac=LOG(zaph2m*zq2m/(c2es*(1._dp+vtmpc1*zq2m)))/zcvm3                 ! ln(e_s) / 610.78 at dewpoint; e_s = e = p r / ( r + Rd/Rv ) = p q / ( Rd/Rv + (1-Rd/Rv) q )
     zdew2l=MIN(zt2l,(tmelt-zfrac*zcvm4)/(1._dp-zfrac))                      ! T at dewpoint according to Murray's Eq. (6)
!
!    water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     lo=zt2w.GT.tmelt                                                        ! Phase of water at 2 m: Ice or liquid water?
     zcvm3=MERGE(c3les,c3ies,lo)                                             ! a in Murray's Eq. (6)
     zcvm4=MERGE(c4les,c4ies,lo)                                             ! b in Murray's Eq. (6)
     zaph2m=paphm1(jl,klevp1)*(1._dp-zhtq/(rd*zt2w                         & ! P(2m) ~ P_S - rho(2m) g 2m
                        *(1._dp+vtmpc1*pqm1(jl,klev)-zx(jl,klev))))          !   rho(2m) ~ P_S / ( Rd T (1 + 0.61 q - q_i) )
     it = NINT(zt2w*1000._dp)                                                ! Integer temperature in mK for lookup tables
     IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.               ! T should lie between 50 and 400 K
     it = MAX(MIN(it,jptlucu2),jptlucu1)
     zqs2=tlucua(it)/zaph2m                                                  ! Rd/Rv * es/p at 2 m height
     zqs2=zqs2/(1._dp-vtmpc1*zqs2)                                           ! Saturation specific humidity at 2 m height
     zq2m=zrh2m*zqs2                                                         ! Specific humidity at 2 m height
     zfrac=LOG(zaph2m*zq2m/(c2es*(1._dp+vtmpc1*zq2m)))/zcvm3                 ! ln(e_s) / 610.78 at dewpoint; e_s = e = p r / ( r + Rd/Rv ) = p q / ( Rd/Rv + (1-Rd/Rv) q )
     zdew2w=MIN(zt2w,(tmelt-zfrac*zcvm4)/(1._dp-zfrac))                      ! T at dewpoint according to Murray's Eq. (6)
!
!    ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     lo=zt2i.GT.tmelt                                                        ! Phase of water at 2 m: Ice or liquid water?
     zcvm3=MERGE(c3les,c3ies,lo)                                             ! a in Murray's Eq. (6)
     zcvm4=MERGE(c4les,c4ies,lo)                                             ! b in Murray's Eq. (6)
     zaph2m=paphm1(jl,klevp1)*(1._dp-zhtq/(rd*zt2i                         & ! P(2m) ~ P_S - rho(2m) g 2m
                        *(1._dp+vtmpc1*pqm1(jl,klev)-zx(jl,klev))))          !   rho(2m) ~ P_S / ( Rd T (1 + 0.61 q - q_i) )
     it = NINT(zt2i*1000._dp)                                                ! Integer temperature in mK for lookup tables
     IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.               ! T should lie between 50 and 400 K
     it = MAX(MIN(it,jptlucu2),jptlucu1)
     zqs2=tlucua(it)/zaph2m                                                  ! Rd/Rv * es/p at 2 m height
     zqs2=zqs2/(1._dp-vtmpc1*zqs2)                                           ! Saturation specific humidity at 2 m height
     zq2m=zrh2m*zqs2                                                         ! Specific humidity at 2 m height
     zfrac=LOG(zaph2m*zq2m/(c2es*(1._dp+vtmpc1*zq2m)))/zcvm3                 ! ln(e_s) / 610.78 at dewpoint; e_s = e = p r / ( r + Rd/Rv ) = p q / ( Rd/Rv + (1-Rd/Rv) q )
     zdew2i=MIN(zt2i,(tmelt-zfrac*zcvm4)/(1._dp-zfrac))                      ! T at dewpoint according to Murray's Eq. (6)
     pdew2(jl)=pfrl(jl)*zdew2l+pfrw(jl)*zdew2w+pfri(jl)*zdew2i               ! Weighted 2 m dewpoint temperature

!
!*       5.9.3   10M WIND COMPONENTS, MAX 10M WINDSPEED
!
!    land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     lo1=zril(jl).GT.0._dp                                                   ! Bulk Richardson number > 0 over land?
     zrat=zhuv/pgeom1(jl,klev)                                               ! 10 m / z - ratio of 10 m vs height of grid center
     zcbn=LOG(1._dp+(EXP(zbnl(jl))-1._dp)*zrat)                              ! ln(1 + 10 m / z0m) = kappa / sqrt(C_N_M(10m))
     zcbs=-(zbnl(jl)-zbml(jl))*zrat                                          ! For stable: kappa / sqrt(C_M(10m)) - kappa / sqrt(C_N_M(10m)) is assumed to be: ( kappa / sqrt(C_M) - kappa / sqrt(C_N_M) ) 10 m / z
     zcbu=-LOG(1._dp+(EXP(zbnl(jl)-zbml(jl))-1._dp)*zrat)                    ! For unstable: using X = kappa / sqrt(C_M) - kappa / sqrt(C_N_M), X(10m) is assumed to be: - ln( (exp(-X) - 1) * 10 m / z + 1)
     zred=(zcbn+MERGE(zcbs,zcbu,lo1))/zbml(jl)                               ! { kappa / sqrt(C_M(10m)) } / { kappa / sqrt(C_M) } = ( V(10m) - V_S ) / ( V_a - V_S ) = V(10m) / V_a
     zu10l=zred*pum1(jl,klev)                                                ! u(10m) = u(10m) / u_a * u_a
     zv10l=zred*pvm1(jl,klev)                                                ! v(10m) = V(10m) / v_a * v_a
     zspeedl=SQRT(zu10l**2+zv10l**2)                                         ! V(10m) = |(u(10m),v(10m))|
!
!    water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     lo1=zriw(jl).GT.0._dp                                                   ! Bulk Richardson number > 0 over water?
     zcbn=LOG(1._dp+(EXP(zbnw(jl))-1._dp)*zrat)                              ! ln(1 + 10 m / z0m) = kappa / sqrt(C_N_M(10m))
     zcbs=-(zbnw(jl)-zbmw(jl))*zrat                                          ! For stable: kappa / sqrt(C_M(10m)) - kappa / sqrt(C_N_M(10m)) is assumed to be: ( kappa / sqrt(C_M) - kappa / sqrt(C_N_M) ) 10 m / z
     zcbu=-LOG(1._dp+(EXP(zbnw(jl)-zbmw(jl))-1._dp)*zrat)                    ! For unstable: using X = kappa / sqrt(C_M) - kappa / sqrt(C_N_M), X(10m) is assumed to be: - ln( (exp(-X) - 1) * 10 m / z + 1)
     zred=(zcbn+MERGE(zcbs,zcbu,lo1))/zbmw(jl)                               ! { kappa / sqrt(C_M(10m)) } / { kappa / sqrt(C_M) } = ( V(10m) - V_S ) / ( V_a - V_S )
     zu10w=pocu(jl)+zred*(pum1(jl,klev)-pocu(jl))                            ! u(10m) = u_S + ( u(10m) - u_S) / ( u_a - u_S) * ( u_a - u_S)
     zv10w=pocv(jl)+zred*(pvm1(jl,klev)-pocv(jl))                            ! v(10m) = v_S + ( v(10m) - v_S) / ( v_a - v_S) * ( v_a - v_S)
     zspeedw=SQRT(zu10w**2+zv10w**2)                                         ! V(10m) = |(u(10m),v(10m))|
     pwind10w(jl)=zred*SQRT((pum1(jl,klev)-pocu(jl))**2                    & ! V(10m) - V_S, i.e. relative wind speed over water at 10 m height
                           +(pvm1(jl,klev)-pocv(jl))**2)
!
!    ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     lo1=zrii(jl).GT.0._dp                                                   ! Bulk Richardson number > 0 over ice?
     zcbn=LOG(1._dp+(EXP(zbni(jl))-1._dp)*zrat)                              ! 10 m / z - ratio of 10 m vs height of grid center
     zcbs=-(zbni(jl)-zbmi(jl))*zrat                                          ! For stable: kappa / sqrt(C_M(10m)) - kappa / sqrt(C_N(10m)) is assumed to be: ( kappa / sqrt(C_M) - kappa / sqrt(C_N) ) 10 m / z
     zcbu=-LOG(1._dp+(EXP(zbni(jl)-zbmi(jl))-1._dp)*zrat)                    ! For unstable: using X = kappa / sqrt(C_M) - kappa / sqrt(C_N), X(10m) is assumed to be: - ln( (exp(-X) - 1) * 10 m / z + 1)
     zred=(zcbn+MERGE(zcbs,zcbu,lo1))/zbmi(jl)                               ! { kappa / sqrt(C_M(10m)) } / { kappa / sqrt(C_M) } = ( V(10m) - V_S ) / ( V_a - V_S )
     zu10i=pocu(jl)+zred*(pum1(jl,klev)-pocu(jl))                            ! u(10m) = u_S + ( u(10m) - u_S) / ( u_a - u_S) * ( u_a - u_S)
     zv10i=pocv(jl)+zred*(pvm1(jl,klev)-pocv(jl))                            ! v(10m) = v_S + ( v(10m) - v_S) / ( v_a - v_S) * ( v_a - v_S)
     zspeedi=SQRT(zu10i**2+zv10i**2)                                         ! V(10m) = |(u(10m),v(10m))|
     pu10(jl)=pfrl(jl)*zu10l+pfrw(jl)*zu10w+pfri(jl)*zu10i                   ! Weighted u(10m)
     pv10(jl)=pfrl(jl)*zv10l+pfrw(jl)*zv10w+pfri(jl)*zv10i                   ! Weighted v(10m)
     pwind10(jl)=pfrl(jl)*zspeedl+pfrw(jl)*zspeedw+pfri(jl)*zspeedi          ! Weighted V(10m) (absolute velocity, not amplitude of averaged wind vector)
  END DO

  IF (lookupoverflow) THEN
     do jl=1,kproma
        if ( INT(ptm1(jl,klev)*1000.) <jptlucu1 .OR.                       &
             INT(ptm1(jl,klev)*1000.) >jptlucu2)                           &
             print*, 'vertex(3) ptm1: ',ptm1(jl,klev)                      &
                ,philon_2d(jl,krow),philat_2d(jl,krow)
        if ( INT(t2l_1d(jl)*1000.) <jptlucu1 .OR.                          &
             INT(t2l_1d(jl)*1000.) >jptlucu2)                              &
             print*, 'vertex(3) t2l: ',t2l_1d(jl)                          &
                ,philon_2d(jl,krow),philat_2d(jl,krow)
        if ( INT(t2w_1d(jl)*1000.) <jptlucu1 .OR.                          &
             INT(t2w_1d(jl)*1000.) >jptlucu2)                              &
             print*, 'vertex(3) t2w: ',t2w_1d(jl)                          &
                ,philon_2d(jl,krow),philat_2d(jl,krow)
        if ( INT(t2i_1d(jl)*1000.) <jptlucu1 .OR.                          &
             INT(t2i_1d(jl)*1000.) >jptlucu2)                              &
             print*, 'vertex(3) t2i: ',t2i_1d(jl)                          &
                ,philon_2d(jl,krow),philat_2d(jl,krow)
     enddo
     CALL error_bi('LOOKUP TABLE OVERFLOW - 3','vertex')
     lookupoverflow = .FALSE.
  ENDIF

!******************************************************************************

  ! Reset pointers used in second call to vertex_vdiff
  NULLIFY(pvdiffp )
  NULLIFY(pvmixtau)
  NULLIFY(pxvar   )
  NULLIFY(zebsh   )

  NULLIFY(pahfl   )
  NULLIFY(pahfli  )
  NULLIFY(pahfll  )
  NULLIFY(pahflw  )
  NULLIFY(pahfs   )
  NULLIFY(pahfsi  )
  NULLIFY(pahfsl  )
  NULLIFY(pahfsw  )
  NULLIFY(pdew2   )
  NULLIFY(pevap   )
  NULLIFY(pevapi  )
  NULLIFY(pevapl  )
  NULLIFY(pevapot )
  NULLIFY(pevapw  )
  NULLIFY(pqhfla  )
  NULLIFY(psnc    )
  NULLIFY(ptemp2  )
  NULLIFY(ptsl    )
  NULLIFY(ptslnew )
  NULLIFY(pu10    )
  NULLIFY(pvdis   )
  NULLIFY(pv10    )
  NULLIFY(pwind10 )
  NULLIFY(pwind10w)
  NULLIFY(zcair   )
  NULLIFY(zqslnew )
! op_pj_20181105+
  NULLIFY(pahfslac)
  NULLIFY(pahfswac)
  NULLIFY(pahfsiac)
  NULLIFY(pahfllac)
  NULLIFY(pahflwac)
  NULLIFY(pahfliac)
  NULLIFY(pevaplac)
  NULLIFY(pevapwac)
  NULLIFY(pevapiac)

!******************************************************************************

  ! Reset pointers used in both separate vertex_vdiff calls
  NULLIFY(ptkem1  )
  NULLIFY(zcfh    )
  NULLIFY(zcfm    )
  NULLIFY(zxtems  )

  NULLIFY(paz0hl  )
  NULLIFY(paz0l   )
  NULLIFY(pcvs    )
  NULLIFY(pcvw    )
  NULLIFY(pocu    )
  NULLIFY(pocv    )
  NULLIFY(psn     )
  NULLIFY(ptsi    )
  NULLIFY(ptslm1  )
  NULLIFY(ptsw    )
  NULLIFY(pws     )
  NULLIFY(pwsmx   )
  NULLIFY(zcfmi   )
  NULLIFY(zcfml   )
  NULLIFY(zcfmw   )
  NULLIFY(zcfhi   )
  NULLIFY(zcfhl   )
  NULLIFY(zcfhw   )
  NULLIFY(zcfnci  )
  NULLIFY(zcfncl  )
  NULLIFY(zcfncw  )
  NULLIFY(zchl    )
  NULLIFY(zdens   )
  NULLIFY(zhum    )
  NULLIFY(zrii    )
  NULLIFY(zril    )
  NULLIFY(zriw    )
  NULLIFY(zsrfll  )
  NULLIFY(zqsi    )
  NULLIFY(zqsl    )
  NULLIFY(zqsw    )

  ! Deallocation of variables that were stored between the two separate vertex_vdiff calls
  DEALLOCATE(pqm1,   ptm1,   pum1,    pvm1,    pxim1,    pxlm1  )

  DEALLOCATE(zcfv,   zcptgz,  zqshear, zx               )

  DEALLOCATE(pfri,   pfrl,   pfrw                               )

  DEALLOCATE(ihpbl,  zbmi,   zbmw,    zbml,    zbhi,     zbhl,   &
             zbhw,   zbni,   zbnl,    zbnw,    zbhnl,    zchi,   &
             zchw,   zcpti,  zcptw,   zdqsl,   zdu2,     zdu2oc, &
             zhsoil, zwet                                       )

  END SUBROUTINE vertex_vdiff_2
!------------------------------------------------------------------------------
#endif
  END SUBROUTINE vertex_vdiff
!==============================================================================

!==============================================================================
  SUBROUTINE vertex_global_end

    IMPLICIT NONE

  END SUBROUTINE vertex_global_end
!==============================================================================
  SUBROUTINE vertex_free_memory

    IMPLICIT NONE

  END SUBROUTINE vertex_free_memory
!==============================================================================
  SUBROUTINE vertex_read_nml_cpl(status, iou)

    ! read namelist for 'coupling'
    !
    ! Author: Tamara Emmerichs, FZJ, Jan 2019

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
    ! switch for skipping calculation of Lagrangian
    ! rate coefficients..it is local,not broadcasted

    NAMELIST /CPL/ imp_lai

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='vertex__read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE vertex_read_nml_cpl

! ===========================================================================

END MODULE messy_vertex_si
