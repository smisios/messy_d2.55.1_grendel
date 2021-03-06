! ***********************************************************************
#include "messy_main_ppd_bi.inc"

MODULE MESSY_RAD_SI
! ***********************************************************************
  ! New structure radiation code RAD
  ! Author: Simone Dietmueller, DLR

#if defined(ECHAM5) || defined(CESM1) || defined(MBM_RAD)

  ! BMIL
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                    , error_bi
#if defined(MBM_RAD)
  USE messy_main_grid_def_bi,   ONLY: nn, nlev, ladd_tte, lv_echam, lzonal_mean
#endif

#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,   ONLY:  mtend_get_handle,             &
                                       mtend_add_l,                  &
                                       mtend_register,               &
                                       mtend_id_t
#endif
  ! SMCL
  USE messy_main_channel,       ONLY:  t_chaobj_cpl
  USE messy_main_constants_mem, ONLY:  STRLEN_ULONG, DP, STRLEN_MEDIUM, solc
  USE messy_main_timer_event,   ONLY:  io_time_event, TRIG_FIRST, time_event
  USE messy_rad,                ONLY:  modstr, rsun_scale, t_rad_work &
                                    ,  NRADCALL, l_switch

  ! SCALING FACTORS FOR THE SHORTWAVE RADIATION PARAMETERISATION.
  USE messy_main_tools,         ONLY:  t_reset_par

  ! SPECIAL FOR FUBRAD SUB-SUBMODEL >>>
  ! SMIL
  USE messy_rad_fubrad_si,      ONLY: lfubrad, rad_fubrad_initialize, &
                                      rad_fubrad_init_memory,         &
                                      rad_fubrad_global_start,        &
                                      rad_fubrad_preprad,             &
                                      rad_fubrad_radheat,             &
                                      rad_fubrad_init_coupling
  ! SPECIAL FOR FUBRAD SUB-SUBMODEL <<<

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC :: TRIM, NULL

  TYPE(io_time_event), PUBLIC :: trigrad = &
                                       io_time_event (3,'steps',TRIG_FIRST,0)
  TYPE(time_event),    PUBLIC :: ev_trigrad
  LOGICAL                     :: l_trigrad
  REAL(DP), POINTER :: dt_offset => NULL()

  CHARACTER(LEN=STRLEN_MEDIUM)  :: sname = ''
  CHARACTER(LEN=2)              :: idx   = ''

  ! IDENTIFIER FOR INPUT OBJECTS
  INTEGER,PARAMETER::id_h2o=1
  INTEGER,PARAMETER::id_co2=2
  INTEGER,PARAMETER::id_ch4=3
  INTEGER,PARAMETER::id_o3=4
  INTEGER,PARAMETER::id_n2o=5
  INTEGER,PARAMETER::id_cfc11=6
  INTEGER,PARAMETER::id_cfc12=7
  INTEGER,PARAMETER::id_aot_lw=8
  INTEGER,PARAMETER::id_aot_sw=9
  INTEGER,PARAMETER::id_gamma_sw=10
  INTEGER,PARAMETER::id_omega_sw=11
  INTEGER,PARAMETER::id_clc=12
  INTEGER,PARAMETER::id_cld_lw=13
  INTEGER,PARAMETER::id_cld_sw=14
  INTEGER,PARAMETER::id_cld_gamma=15
  INTEGER,PARAMETER::id_cld_omega=16
  INTEGER,PARAMETER::id_cld_clcv=17
  INTEGER,PARAMETER::id_cld_cld=18
  INTEGER,PARAMETER::id_o2=19
  INTEGER, PARAMETER:: NMAXINP=19     ! no. of input objects

  ! IDENTIFIER FOR OUTPUT OBJECTS
  ! xradout with the indices ido_emter and ido_emtef are not the
  ! LW emissivities, as the output from rad_rad_smcl are the
  ! net LW fluxes (ido_flxt and ido_flxtf) output of
  ! LW emissivities is obsolete.

  INTEGER,PARAMETER::ido_trsol=1
  INTEGER,PARAMETER::ido_trsof=2
  INTEGER,PARAMETER::ido_flxs=3
  INTEGER,PARAMETER::ido_flxt=4
  INTEGER,PARAMETER::ido_flxsf=5
  INTEGER,PARAMETER::ido_flxtf=6
  INTEGER,PARAMETER::ido_dtsw=7
  INTEGER,PARAMETER::ido_dtlw=8
  INTEGER,PARAMETER::ido_dtswc=9
  INTEGER,PARAMETER::ido_dtlwc=10
  INTEGER,PARAMETER::ido_addst=11
  INTEGER,PARAMETER::ido_tte=12
  INTEGER,PARAMETER::ido_addsth=13
  INTEGER,PARAMETER::ido_trnir=14
  INTEGER,PARAMETER::ido_trnif=15
  INTEGER,PARAMETER::ido_trsw1=16
  INTEGER,PARAMETER::ido_trs1f=17
  INTEGER,PARAMETER::ido_flxnir=18
  INTEGER,PARAMETER::ido_flxsw1=19
  INTEGER,PARAMETER::ido_heats1=20
  INTEGER,PARAMETER::ido_heatni=21
  INTEGER,PARAMETER::ido_srad0u=22
  INTEGER,PARAMETER::ido_sradsu=23
  INTEGER,PARAMETER::ido_tradsu=24
  INTEGER,PARAMETER::ido_soflw=25
  INTEGER,PARAMETER::ido_soflwac=26
  INTEGER,PARAMETER::ido_sofli=27
  INTEGER,PARAMETER::ido_sofliac=28
  INTEGER,PARAMETER::ido_sofll=29
  INTEGER,PARAMETER::ido_sofllac=30
  INTEGER,PARAMETER::ido_trflw=31
  INTEGER,PARAMETER::ido_trflwac=32
  INTEGER,PARAMETER::ido_trfli=33
  INTEGER,PARAMETER::ido_trfliac=34
  INTEGER,PARAMETER::ido_trfll=35
  INTEGER,PARAMETER::ido_trfllac=36

  INTEGER,PARAMETER::ido_flxus=37
  INTEGER,PARAMETER::ido_flxut=38
  INTEGER,PARAMETER::ido_flxusf=39
  INTEGER,PARAMETER::ido_flxutf=40
  INTEGER,PARAMETER::ido_flxuni=41
  INTEGER,PARAMETER::ido_flxunif=42
  INTEGER,PARAMETER::ido_trus=43
  INTEGER,PARAMETER::ido_trusf=44
  INTEGER,PARAMETER::ido_truni=45
  INTEGER,PARAMETER::ido_trunif=46

  INTEGER, PARAMETER:: NMAXOUT=46     ! no. of output objects

  CHARACTER(LEN=12), DIMENSION(NMAXOUT):: setname= &
     (/'trsol       ','trsof       ','flxs        ','flxt        ',&
       'flxsf       ','flxtf       ','dtdt_sw     ','dtdt_lw     ',&
       'dtdt_swc    ','dtdt_lwc    ','addst       ','tte         ',&
       'addsth      ','trnir       ','trnif       ','trsw1       ',&
       'trs1f       ','flxnir      ','flxsw1      ','heats1      ',&
       'heatni      ','srad0u      ','sradsu      ','tradsu      ',&
       'soflw       ','soflwac     ','sofli       ','sofliac     ',&
       'sofll       ','sofllac     ','trflw       ','trflwac     ',&
       'trfli       ','trfliac     ','trfll       ','trfllac     ',&
       'flxus       ','flxut       ','flxusf      ','flxutf      ',&
       'flxuni      ','flxunif     ','trus        ','trusf       ',&
       'truni       ','trunif      ' &
      /)

  CHARACTER(LEN=16), DIMENSION(NMAXOUT):: repr_name= &
       (/'GP_3D_INT       ','GP_3D_INT       ','GP_3D_INT       ','GP_3D_INT       ',&
         'GP_3D_INT       ','GP_3D_INT       ','GP_3D_MID       ','GP_3D_MID       ',&
         'GP_3D_MID       ','GP_3D_MID       ','GP_3D_MID       ','GP_3D_MID       ' &
       , 'GP_3D_INT       ','GP_3D_INT       ','GP_3D_INT       ','GP_3D_INT       ' &
       , 'GP_3D_INT       ','GP_3D_INT       ','GP_3D_INT       ','GP_3D_MID       ' &
       , 'GP_3D_MID       ','GP_2D_HORIZONTAL','GP_2D_HORIZONTAL','GP_2D_HORIZONTAL' &
       , 'GP_2D_HORIZONTAL','GP_2D_HORIZONTAL','GP_2D_HORIZONTAL','GP_2D_HORIZONTAL' &
       , 'GP_2D_HORIZONTAL','GP_2D_HORIZONTAL','GP_2D_HORIZONTAL','GP_2D_HORIZONTAL' &
       , 'GP_2D_HORIZONTAL','GP_2D_HORIZONTAL','GP_2D_HORIZONTAL','GP_2D_HORIZONTAL' &
       , 'GP_3D_INT       ','GP_3D_INT       ','GP_3D_INT       ','GP_3D_INT       ' &
       , 'GP_3D_INT       ','GP_3D_INT       ','GP_3D_INT       ','GP_3D_INT       ' &
       , 'GP_3D_INT       ','GP_3D_INT       '   &
       /)

  INTEGER, DIMENSION(NMAXOUT):: repr_idx

  CHARACTER(LEN=29), DIMENSION(NMAXOUT):: LONGNAME=(/&
    'transmissivity               ', 'transmissivity clearsky      ', &
    'shortwave flux               ', 'longwave flux                ', &
    'shortwave flux clearsky      ', 'longwave flux clearsky       ', &
    'shortwave heating rate       ', 'longwave heating rate        ', &
    'shortwave heating rate cl.sky', 'longwave heating rate cl.sky ', &
    'adjusted temp.               ', &
    'temperature tendency         ', &
    'temperature half level       ', &
    'transmissivity NIR           ', &
    'transmis. NIR clearsky       ', &
    'transmissivity UVVIS         ', &
    'transmis UVVIS clearsky      ', &
    'NIR flux all sky             ', &
    'SW1 flux all sky             ', &
    'UV-Vis heating all sky       ', &
    'NIR heating all sky          ', &
    'top solar rad. upward        ', &
    'surf. solar rad. upward      ', &
    'surf. thermal rad. upward    ', &
    'sw flux over water           ', &
    'sw flux over water (masked)  ', &
    'sw flux over ice             ', &
    'sw flux over ice (masked)    ', &
    'sw flux over land            ', &
    'sw flux over land (masked)   ', &
    'lw flux over water           ', &
    'lw flux over water (masked)  ', &
    'lw flux over ice             ', &
    'lw flux over ice (masked)    ', &
    'lw flux over land            ', &
    'lw flux over land (masked)   ', &
    'shortwave flux upwrd         ', 'longwave flux upwrd          ', &
    'shortwave flux upwrd clearsky', 'longwave flux upwrd clear sky', &
    'NIR flux upwrd all sky       ', 'NIR flux upwrd clear sky     ', &
    'transmissivity SW upward     ', 'transmissivity SW up. clear s', &
    'transmissivity NIR upward    ', 'transmissivity NIR up. cl. s.'  &
    /)

  CHARACTER(LEN=6), DIMENSION(NMAXOUT):: unit=(/&
     '-     ', '-     ', 'W/m**2', &
     'W/m**2', 'W/m**2', 'W/m**2', 'K/s   ', 'K/s   ', &
     'K/s   ', 'K/s   ', 'K     ', 'K/s   ', 'K     ', &
     '-     ', '-     ', &
     '-     ', &
     '-     ', &
     'W/m**2', &
     'W/m**2', &
     'K/s   ', &
     'K/s   ', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', &
     'W/m**2', 'W/m**2', 'W/m**2', 'W/m**2', &
     'W/m**2', 'W/m**2', '-     ', '-     ', &
     '-     ', '-     '  &
     /)

  ! fixed tropopause used for adjusted RF diagnostics
  TYPE(t_chaobj_cpl)                                  :: tp_fixed
  TYPE(t_chaobj_cpl), DIMENSION(NMAXINP,NRADCALL)     :: r_inp
  TYPE(t_chaobj_cpl), DIMENSION(NMAXINP,NRADCALL)     :: s_inp
  TYPE(t_rad_work),   DIMENSION(NMAXINP,NRADCALL)     :: xradin
  TYPE(t_rad_work),   DIMENSION(NMAXOUT,NRADCALL)     :: xradout
  INTEGER, DIMENSION(NRADCALL) :: i_rad = 1
  INTEGER, DIMENSION(NRADCALL) :: i_sw  = 1
  INTEGER, DIMENSION(NRADCALL) :: i_lw  = 1

  TYPE(t_reset_par) :: rset_solc    = t_reset_par(.FALSE.,solc)
  TYPE(t_reset_par) :: rset_calbmns = t_reset_par(.FALSE.,0._dp)
  TYPE(t_reset_par) :: rset_calbmxs = t_reset_par(.FALSE.,0._dp)
  TYPE(t_reset_par) :: rset_calbmni = t_reset_par(.FALSE.,0._dp)

  !pointer for coupling
  REAL(DP), DIMENSION(:,:), POINTER:: alsol_x => NULL(),alsoi_x  => NULL(),&
                                      alsow_x => NULL(),albedo_x => NULL()
  REAL(dp), POINTER :: rdayl_x(:,:)    => NULL()
  REAL(dp), POINTER :: cosszac_x(:,:)  => NULL()
  REAL(dp), POINTER :: cosszacm_x(:,:) => NULL()
  REAL(dp), POINTER :: cdisse_x        => NULL()
  REAL(dp), POINTER :: cdissem_x       => NULL()
  REAL(dp), POINTER :: tp_fixed_x(:,:) => NULL()
  REAL(dp), POINTER :: tp_fixed_i(:,:) => NULL()
  REAL(dp), POINTER :: zi0_2d(:,:)     => NULL()

  ! unperturbed temperture
  REAL(dp), POINTER, DIMENSION(:,:,:) ::  padutm1  => NULL(), &
                                          padutm1o => NULL()

  CHARACTER(LEN=STRLEN_ULONG), DIMENSION(NMAXINP,NRADCALL) :: att_unit = ''

  ! feedback to basemodel (TRUE) or purely diagnostic (FALSE)
  LOGICAL :: l_feedback = .TRUE.

  ! switch off (l_quiet=.TRUE.) messages to log-file during time loop
  LOGICAL :: l_quiet = .FALSE.

#ifdef MESSYTENDENCY
  ! variable for tendency budget
  integer                         :: my_handle
#endif

  ! SUBROUTINES
  PUBLIC  :: rad_initialize
  PUBLIC  :: rad_init_memory
  PUBLIC  :: rad_init_coupling
  PUBLIC  :: rad_global_start
  PUBLIC  :: rad_radiation
  PUBLIC  :: rad_radheat
#if defined(MBM_RAD)
  PUBLIC  :: rad_mbm_initialize
#endif
  !PRIVATE :: rad_read_nml_cpl
  !PRIVATE :: rad_read_nml_cpl_mbm

CONTAINS

  ! ====================================================================
  SUBROUTINE rad_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    ! BMIL
    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_timer_bi,      ONLY: p_bcast_event, timer_event_init
    USE messy_main_grid_def_mem_bi,ONLY: nn
    USE messy_main_data_bi,       ONLY: lcouple

    ! SMCL
    USE messy_main_tools,         ONLY: find_next_free_unit
    USE messy_rad,                ONLY: rad_sw_initialize, rad_lw_initialize
    USE messy_rad,                ONLY: su_albedo_rad ! op_pj_20130407
    USE messy_rad,                ONLY: calbmns, calbmxs, calbmni ! op_mk_20190118
    USE messy_rad,                ONLY: &
         absa_1, absb_1, fracrefa_1, fracrefb_1, forref_1, selfref_1, &
         absa_2, absb_2, fracrefa_2, fracrefb_2, forref_2, selfref_2, &
         refparam_2, &
         absa_3, absb_3, fracrefa_3, fracrefb_3, forref_3, selfref_3, &
         absn2oa_3, &
         absn2ob_3, etaref_3, h2oref_3, n2oref_3, co2ref_3, strrat_3, &
         absa_4, absb_4, fracrefa_4, fracrefb_4, selfref_4, strrat1_4, &
         strrat2_4, &
         absa_5, absb_5, ccl4_5, fracrefa_5, fracrefb_5, selfref_5, &
         strrat1_5, strrat2_5, &
         absa_6, absco2_6, cfc11adj_6, cfc12_6, fracrefa_6, selfref_6, &
         absa_7, absb_7, absco2_7, fracrefa_7, fracrefb_7, selfref_7, &
         strrat_7, &
         absa_8, absb_8, fracrefa_8, fracrefb_8, selfref_8, absco2a_8, &
         absco2b_8, &
         absn2oa_8, absn2ob_8, cfc12_8, cfc22adj_8, h2oref_8, n2oref_8, &
         o3ref_8, &
         absa_9, absb_9, fracrefa_9, fracrefb_9, selfref_9, absn2o_9, &
         ch4ref_9, &
         etaref_9, h2oref_9, n2oref_9, strrat_9, &
         absa_10, absb_10, fracrefa_10, fracrefb_10, &
         absa_11, absb_11, fracrefa_11, selfref_11, fracrefb_11, &
         absa_12, fracrefa_12, selfref_12, strrat_12, &
         absa_13, fracrefa_13, selfref_13, strrat_13, absa_14, absb_14, &
         fracrefa_14, fracrefb_14, selfref_14, &
         absa_15, fracrefa_15, selfref_15, strrat_15, &
         absa_16, fracrefa_16, selfref_16, strrat_16,        &
         NGC, NGS, NGM, NGN, NGB, WT,                        &
         corr1, corr2,                                       &
         PREF, PREFLOG, TREF,                                &
         NG, NSPA, NSPB, WAVENUM1, WAVENUM2, DELWAVE,        &
         TOTPLNK, TOTPLK16

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'rad_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: j1,j2
    ! check event with present date
    CHARACTER(LEN=*), PARAMETER :: EV_TLEV_PRES = 'present'

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

#if defined(MBM_RAD)
    ! READ CPL_MBM namelist
    CALL rad_mbm_initialize
#endif

    ! READ CPL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find next free I/O unit
       CALL rad_read_nml_cpl(status, iou)    ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
    END IF

    ! longwave
    IF (p_parallel_io) THEN
       CALL rad_lw_initialize
    END IF

    !
    ! p_bcast of rad_lon_SURRTFTR
    !
    CALL p_bcast(NGC, p_io)
    CALL p_bcast(NGS, p_io)
    CALL p_bcast(NGM, p_io)
    CALL p_bcast(NGN, p_io)
    CALL p_bcast(NGB, p_io)
    CALL p_bcast(WT, p_io)

    !
    ! p_bcast of rad_lon_surrtbg2
    !
    CALL p_bcast(corr1, p_io)
    CALL p_bcast(corr2, p_io)

    !
    ! p_bcast of rad_lon_SURRTRF
    !
    CALL p_bcast(PREF, p_io)
    CALL p_bcast(PREFLOG, p_io)
    CALL p_bcast(TREF, p_io)

    !
    ! p_bcast of rad_lon_SURRTPK
    !
    CALL p_bcast(NG, p_io)
    CALL p_bcast(NSPA, p_io)
    CALL p_bcast(NSPB, p_io)
    CALL p_bcast(WAVENUM1, p_io)
    CALL p_bcast(WAVENUM2, p_io)
    CALL p_bcast(DELWAVE, p_io)
    CALL p_bcast(TOTPLNK, p_io)
    CALL p_bcast(TOTPLK16, p_io)

    !
    ! p_bcast of read_rrta1
    !
    CALL p_bcast(absa_1, p_io)
    CALL p_bcast(absb_1, p_io)
    CALL p_bcast(fracrefa_1, p_io)
    CALL p_bcast(fracrefb_1, p_io)
    CALL p_bcast(forref_1, p_io)
    CALL p_bcast(selfref_1, p_io)

    !
    ! p_bcast of read_rrta2
    !
    CALL p_bcast(absa_2, p_io)
    CALL p_bcast(absb_2, p_io)
    CALL p_bcast(fracrefa_2, p_io)
    CALL p_bcast(fracrefb_2, p_io)
    CALL p_bcast(forref_2, p_io)
    CALL p_bcast(selfref_2, p_io)
    CALL p_bcast(refparam_2, p_io)

    !
    ! p_bcast of read_rrta3
    !
    CALL p_bcast(absa_3, p_io)
    CALL p_bcast(absb_3, p_io)
    CALL p_bcast(fracrefa_3, p_io)
    CALL p_bcast(fracrefb_3, p_io)
    CALL p_bcast(forref_3, p_io)
    CALL p_bcast(selfref_3, p_io)
    CALL p_bcast(absn2oa_3, p_io)
    CALL p_bcast(absn2ob_3, p_io)
    CALL p_bcast(etaref_3, p_io)
    CALL p_bcast(h2oref_3, p_io)
    CALL p_bcast(n2oref_3, p_io)
    CALL p_bcast(co2ref_3, p_io)
    CALL p_bcast(strrat_3, p_io)

    !
    ! p_bcast of read_rrta4
    !
    CALL p_bcast(absa_4, p_io)
    CALL p_bcast(absb_4, p_io)
    CALL p_bcast(fracrefa_4, p_io)
    CALL p_bcast(fracrefb_4, p_io)
    CALL p_bcast(selfref_4, p_io)
    CALL p_bcast(strrat1_4, p_io)
    CALL p_bcast(strrat2_4, p_io)

    !
    ! p_bcast of read_rrta5
    !
    CALL p_bcast(absa_5, p_io)
    CALL p_bcast(absb_5, p_io)
    CALL p_bcast(ccl4_5, p_io)
    CALL p_bcast(fracrefa_5, p_io)
    CALL p_bcast(fracrefb_5, p_io)
    CALL p_bcast(selfref_5, p_io)
    CALL p_bcast(strrat1_5, p_io)
    CALL p_bcast(strrat2_5, p_io)

    !
    ! p_bcast of read_rrta6
    !
    CALL p_bcast(absa_6, p_io)
    CALL p_bcast(absco2_6, p_io)
    CALL p_bcast(cfc11adj_6, p_io)
    CALL p_bcast(cfc12_6, p_io)
    CALL p_bcast(fracrefa_6, p_io)
    CALL p_bcast(selfref_6, p_io)

    !
    ! p_bcast of read_rrta7
    !
    CALL p_bcast(absa_7, p_io)
    CALL p_bcast(absb_7, p_io)
    CALL p_bcast(absco2_7, p_io)
    CALL p_bcast(fracrefa_7, p_io)
    CALL p_bcast(fracrefb_7, p_io)
    CALL p_bcast(selfref_7, p_io)
    CALL p_bcast(strrat_7, p_io)

    !
    ! p_bcast of read_rrta8
    !
    CALL p_bcast(absa_8, p_io)
    CALL p_bcast(absb_8, p_io)
    CALL p_bcast(fracrefa_8, p_io)
    CALL p_bcast(fracrefb_8, p_io)
    CALL p_bcast(selfref_8, p_io)
    CALL p_bcast(absco2a_8, p_io)
    CALL p_bcast(absco2b_8, p_io)
    CALL p_bcast(absn2oa_8, p_io)
    CALL p_bcast(absn2ob_8, p_io)
    CALL p_bcast(cfc12_8, p_io)
    CALL p_bcast(cfc22adj_8, p_io)
    CALL p_bcast(h2oref_8, p_io)
    CALL p_bcast(n2oref_8, p_io)
    CALL p_bcast(o3ref_8, p_io)

    !
    ! p_bcast of read_rrta9
    !
    CALL p_bcast(absa_9, p_io)
    CALL p_bcast(absb_9, p_io)
    CALL p_bcast(fracrefa_9, p_io)
    CALL p_bcast(fracrefb_9, p_io)
    CALL p_bcast(selfref_9, p_io)
    CALL p_bcast(absn2o_9, p_io)
    CALL p_bcast(ch4ref_9, p_io)
    CALL p_bcast(etaref_9, p_io)
    CALL p_bcast(h2oref_9, p_io)
    CALL p_bcast(n2oref_9, p_io)
    CALL p_bcast(strrat_9, p_io)

    !
    ! p_bcast of read_rrta10
    !
    CALL p_bcast(absa_10, p_io)
    CALL p_bcast(absb_10, p_io)
    CALL p_bcast(fracrefa_10, p_io)
    CALL p_bcast(fracrefb_10, p_io)

    !
    ! p_bcast of read_rrta11
    !
    CALL p_bcast(absa_11, p_io)
    CALL p_bcast(absb_11, p_io)
    CALL p_bcast(fracrefa_11, p_io)
    CALL p_bcast(selfref_11, p_io)
    CALL p_bcast(fracrefb_11, p_io)

    !
    ! p_bcast of read_rrta12
    !
    CALL p_bcast(absa_12, p_io)
    CALL p_bcast(fracrefa_12, p_io)
    CALL p_bcast(selfref_12, p_io)
    CALL p_bcast(strrat_12, p_io)

    !
    ! p_bcast of read_rrta13
    !
    CALL p_bcast(absa_13, p_io)
    CALL p_bcast(fracrefa_13, p_io)
    CALL p_bcast(selfref_13, p_io)
    CALL p_bcast(strrat_13, p_io)

    !
    ! p_bcast of read_rrta14
    !
    CALL p_bcast(absa_14, p_io)
    CALL p_bcast(absb_14, p_io)
    CALL p_bcast(fracrefa_14, p_io)
    CALL p_bcast(fracrefb_14, p_io)
    CALL p_bcast(selfref_14, p_io)

    !
    ! p_bcast of read_rrta15
    !
    CALL p_bcast(absa_15, p_io)
    CALL p_bcast(fracrefa_15, p_io)
    CALL p_bcast(selfref_15, p_io)
    CALL p_bcast(strrat_15, p_io)

    !
    ! p_bcast of read_rrta16
    !
    CALL p_bcast(absa_16, p_io)
    CALL p_bcast(fracrefa_16, p_io)
    CALL p_bcast(selfref_16, p_io)
    CALL p_bcast(strrat_16, p_io)

    ! shortwave
    CALL rad_sw_initialize

    ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(l_feedback,p_io)
    CALL p_bcast(l_quiet,p_io)
    CALL p_bcast(lfubrad,p_io)
    CALL p_bcast_event(trigrad,p_io)
    CALL p_bcast (rsun_scale, p_io)
    CALL p_bcast (rset_solc%l, p_io)
    CALL p_bcast (rset_solc%v, p_io)

    CALL p_bcast (rset_calbmns%l, p_io)
    CALL p_bcast (rset_calbmns%v, p_io)
    CALL p_bcast (rset_calbmxs%l, p_io)
    CALL p_bcast (rset_calbmxs%v, p_io)
    CALL p_bcast (rset_calbmni%l, p_io)
    CALL p_bcast (rset_calbmni%v, p_io)
    CALL p_bcast(tp_fixed%cha, p_io)
    CALL p_bcast(tp_fixed%obj, p_io)

    DO  j2=1,NRADCALL
       CALL p_bcast(l_switch(j2),p_io)

       CALL p_bcast(i_rad(j2),p_io)
       IF ( (i_rad(j2) < 1) .OR. (i_rad(j2) > 2)) THEN
          CALL error_bi('i_rad(.) < 1 OR i_rad(.) > 2 in rad.nml', substr)
       ENDIF

       CALL p_bcast(i_sw(j2),p_io)
       IF ( (i_sw(j2) < 1) .OR. (i_sw(j2) > 2)) THEN
          CALL error_bi('i_sw(.) < 1 OR i_sw(.) > 2 in rad.nml', substr)
       ENDIF

       CALL p_bcast(i_lw(j2),p_io)
       IF ( (i_lw(j2) < 1) .OR. (i_lw(j2) > 1)) THEN
          CALL error_bi('i_lw(.) < 1 OR i_lw(.) > 1 in rad.nml', substr)
       ENDIF
    END DO

    IF (l_feedback .AND. (.NOT. l_switch(1))) THEN
       CALL error_bi('l_feedback = T requires l_switch(1) = T in rad.nml' &
            , substr)
    ELSE
       IF (i_rad(1) /= 1) THEN
          CALL error_bi('l_feedback = T requires i_rad(1) = 1 in rad.nml' &
               , substr)
       END IF
    END IF

    DO j2=1,NRADCALL
       DO j1=1,NMAXINP
          CALL p_bcast(r_inp(j1,j2)%cha, p_io)
          CALL p_bcast(r_inp(j1,j2)%obj, p_io)
       END DO
    END DO

    DO j2=2,NRADCALL
       DO j1=1,NMAXINP
          ! reference is always first radiation call
          IF( TRIM(r_inp(j1,j2)%cha) == '' ) THEN
             r_inp(j1,j2)%cha = r_inp(j1,1)%cha
          ENDIF
          IF( TRIM(r_inp(j1,j2)%obj) == '' ) THEN
             r_inp(j1,j2)%obj = r_inp(j1,1)%obj
          END IF
       END DO
    END DO

    ! initialize radiation event
    CALL timer_event_init(ev_trigrad, trigrad, &
         'radiation computation', EV_TLEV_PRES)

    ! initialize resolution dependent albedos
    CALL su_albedo_rad(status, nn, lcouple)
    IF (status /= 0) CALL error_bi('su_albedo_rad reported an error',substr)

    IF (rset_calbmns%l) calbmns = rset_calbmns%v
    IF (rset_calbmxs%l) calbmxs = rset_calbmxs%v
    IF (rset_calbmni%l) calbmni = rset_calbmni%v
    IF (p_parallel_io) THEN
       WRITE(*,*) "rset_calbmns = ",rset_calbmns," calbmns = ",calbmns
       WRITE(*,*) "rset_calbmxs = ",rset_calbmxs," calbmxs = ",calbmxs
       WRITE(*,*) "rset_calbmni = ",rset_calbmni," calbmni = ",calbmni
    END IF

    IF (lfubrad) CALL rad_fubrad_initialize

    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE rad_initialize
  ! ====================================================================
#if defined(MBM_RAD)
  SUBROUTINE rad_mbm_initialize
    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CPL_MBM-namelist
    ! ------------------------------------------------------------------

    ! BMIL
    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast

    ! SMCL
    USE messy_main_tools,         ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'rad_mbm_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    ! READ CPL_MBM namelist
    IF (p_parallel_io) THEN                     ! read only on I/O-PE
       iou = find_next_free_unit(100,200)       ! find next free I/O unit
       CALL rad_read_nml_cpl_mbm(status, iou)   ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CPL_MBM namelist',substr)
    END IF

    CALL p_bcast(ladd_tte,    p_io)
    CALL p_bcast(lv_echam,    p_io)
    CALL p_bcast(lzonal_mean, p_io)
    CALL p_bcast(nn,          p_io)
    CALL p_bcast(nlev,        p_io)

  END SUBROUTINE rad_mbm_initialize
#endif
  ! ====================================================================
  SUBROUTINE rad_init_memory

    ! BMIL
    USE messy_main_mpi_bi,            ONLY: p_parallel_io
    USE messy_main_blather_bi,        ONLY: start_message_bi, end_message_bi
    USE messy_main_channel_error_bi,  ONLY: channel_halt
    USE messy_main_channel_bi,        ONLY: GP_2D_HORIZONTAL, SCALAR, GP_3D_MID
    USE messy_main_channel,           ONLY: new_channel, new_channel_object &
                                          , new_attribute, REPR_UNDEF
    USE messy_main_channel_repr,      ONLY: get_representation_id &
                                          , get_representation_info
    USE messy_main_tools,             ONLY: int2str

    IMPLICIT NONE
    INTRINSIC :: ANY

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER           :: substr = 'rad_init_memory'
    INTEGER                               :: status
    INTEGER                               :: i,j2
    INTEGER, DIMENSION(NMAXOUT)           :: irank

    !following channel objects map ECHAM radiation diagnostics
    !(shortend, some variables can be reached by postprocessing 3D-fluxes)
    !ac -> laccu=true

    ! moved from initialize
#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
    CALL mtend_register (my_handle, mtend_id_t)
#endif
     CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)

    IF (p_parallel_io) WRITE(*,*) 'add new channel rad ...'

    CALL new_channel(status, modstr, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'dt_offset' &
         , p0 = dt_offset, reprid=SCALAR)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dt_offset' &
           , 'long_name', c='offset for radiation calculation (orbit)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dt_offset', 'units', c='s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'alsol' &
         , p2 = alsol_x,reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'alsol' &
           , 'long_name', c='albedo over land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'alsol', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'alsow' &
         , p2 = alsow_x,reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'alsow' &
           , 'long_name', c='albedo over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'alsow', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'alsoi' &
         , p2 = alsoi_x,reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'alsoi' &
           , 'long_name', c='albedo over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'alsoi', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'albedo' &
         , p2 = albedo_x,reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'albedo' &
           , 'long_name', c='surface albedo')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'albedo', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'zi0', &
         p2=zi0_2d, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zi0', &
         'long_name', c='solar incidence')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zi0', 'units', c='W/m**2')
    CALL channel_halt(substr, status)

    IF ( ANY(i_rad(:) > 1) ) THEN

       CALL new_channel_object(status, modstr, 'adutm1' &
            , p3 =padutm1 ,reprid=GP_3D_MID, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'adutm1' &
            , 'long_name', c='unperturbed temp. (=tm1)')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'adutm1', 'units', c='K')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'adutm1o' &
            , p3 =padutm1o ,reprid=GP_3D_MID, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'adutm1o' &
            , 'long_name', c='unperturbed temperature at t-1')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'adutm1o', 'units', c='K')
       CALL channel_halt(substr, status)

    END IF

    ! get representation indices and ranks
    repr_idx(:) = REPR_UNDEF
    DO i=1, NMAXOUT
       CALL get_representation_id(status, TRIM(repr_name(i)), repr_idx(i))
       CALL get_representation_info(status, inpname='' &
            , id=repr_idx(i), rank=irank(i))
       CALL channel_halt(substr, status)
    END DO

    rad_calls: DO j2=1,NRADCALL
       IF (.NOT. l_switch(j2)) CYCLE

       CALL int2str(idx, j2, '0')
       sname=modstr//idx

       IF (p_parallel_io) WRITE(*,*) 'add new channel '//sname//'....'

       CALL new_channel(status, TRIM(sname),lrestreq=.TRUE.)
       CALL channel_halt(substr, status)

       out_objs: DO i=1,NMAXOUT
          SELECT CASE (irank(i))
          CASE(0)
             CALL new_channel_object(status, TRIM(sname), TRIM(setname(i)) &
                  , p0=xradout(i,j2)%ptr0,reprid=repr_idx(i))
          CASE(2)
             CALL new_channel_object(status, TRIM(sname), TRIM(setname(i)) &
                  , p2=xradout(i,j2)%ptr2,reprid=repr_idx(i))
          CASE(3)
             CALL new_channel_object(status, TRIM(sname), TRIM(setname(i)) &
                  , p3=xradout(i,j2)%ptr3,reprid=repr_idx(i))
          CASE(4)
             CALL new_channel_object(status, TRIM(sname), TRIM(setname(i)) &
                  , p4=xradout(i,j2)%ptr4,reprid=repr_idx(i))
          END SELECT
          CALL channel_halt(substr//'channel object '//setname(i)//&
               &' not found', status)

          CALL new_attribute(status, TRIM(sname), TRIM(setname(i)) &
               ,'long_name', c=TRIM(LONGNAME(i)))
          CALL channel_halt(substr, status)
          CALL new_attribute(status, TRIM(sname), TRIM(setname(i)) &
               ,'units', c=TRIM(unit(i)))
          CALL channel_halt(substr, status)
       END DO out_objs

    END DO rad_calls

    ! create channel objects for each radiation call
    IF (lfubrad) CALL rad_fubrad_init_memory

    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)

  END SUBROUTINE rad_init_memory
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE rad_init_coupling

    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the
    ! basemodel and to other submodes.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, GP_2D_HORIZONTAL
    USE messy_main_channel,          ONLY: get_channel_object &
                                         , new_channel_object &
                                         , get_channel_object_info &
                                         , get_attribute ,new_attribute &
                                         , get_channel_info
    USE messy_main_channel_repr,     ONLY: get_representation_info &
                                         , get_representation_id
    USE messy_main_blather_bi,       ONLY: warning_bi
#ifdef CESM1
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, ngpblks, nlev
#endif

    ! SMCL
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, M_H2O, M_air &
                                      , STRLEN_ULONG
    USE messy_main_tools,         ONLY: strcrack,int2str
#ifdef CESM1
    USE messy_rad,                ONLY: nsw, jpband
#endif

    IMPLICIT NONE
    INTRINSIC :: ANY

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER           :: substr = 'rad_init_coupling'
    INTEGER                               :: status
    INTEGER                               :: j1,j2
    INTEGER                               :: m
    INTEGER                               :: repr_id
    INTEGER                               :: irank
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER :: outstring(:) => NULL()
    LOGICAL                               :: l_assoc
    CHARACTER(LEN=STRLEN_ULONG)           :: msg = ''

    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output

    IF (l_feedback) THEN
       CALL get_channel_info(status, 'rad4all')
       IF (status == 0) THEN
          CALL error_bi('rad4all is also running (always with feedback);'//&
               &'either switch off rad4all or set l_feedback=F in rad.nml' &
               , substr)
       END IF
    END IF

    ! distance Sun - Earth in AU
    ! ... in the middle of radiation time interval
    CALL get_channel_object(status, 'orbit', 'cdissem', p0=cdissem_x)
    CALL channel_halt(substr//': object cdissem in channel orbit not found!',&
        status)
    ! ... instantaneous
    CALL get_channel_object(status, 'orbit', 'cdisse', p0=cdisse_x)
    CALL channel_halt(substr//': object cdisse in channel orbit not found!',&
         status)

    ! cos zenith angle
    ! ... in the middle of radiation time interval
    CALL get_channel_object(status, 'orbit', 'cosszacm', p2=cosszacm_x)
    CALL channel_halt(substr//': object cosszacm in channel orbit not found!',&
         status)

    ! ... instantaneous
    CALL get_channel_object(status, 'orbit', 'cosszac', p2=cosszac_x)
    CALL channel_halt(substr//': object cosszac in channel orbit not found!',&
         status)

    !
    CALL get_channel_object(status, 'orbit', 'rdayl', p2=rdayl_x)
    CALL channel_halt(substr//': object rdayl in channel orbit not found!',&
         status)

    IF (ANY(i_rad(:) == 2)) THEN
       IF ((TRIM(tp_fixed%cha) == '') .OR. (TRIM(tp_fixed%obj) == '')) THEN
          CALL error_bi('tp_fixed required in CPL namelist for adjusted RF' &
               ,substr)
       ENDIF
       CALL get_channel_object(status, TRIM(tp_fixed%cha), TRIM(tp_fixed%obj) &
            , p2=tp_fixed_x)
       CALL channel_halt(&
            substr//': object '//TRIM(tp_fixed%obj)//' in channel '//&
            &TRIM(tp_fixed%cha)//' not found!',&
            status)

       CALL new_channel_object(status,modstr,'tp_fixed_i', p2=tp_fixed_i, &
            reprid=GP_2D_HORIZONTAL)
       CALL channel_halt(substr//' tp_fixed_i  ...', status)
    ENDIF

    rad_calls: DO j2=1,NRADCALL
       IF (.NOT. l_switch(j2)) CYCLE

       inp_objs: DO j1=1,NMAXINP
          CALL int2str(idx, j2, '0')
          sname=modstr//idx

          IF ( (TRIM(r_inp(j1,j2)%cha) .EQ. '#const') .OR. &
               (TRIM(r_inp(j1,j2)%cha) .EQ. '#vgrad') ) THEN

             CALL strcrack(TRIM(r_inp(j1,j2)%obj), '=', outstring, m)
             IF (m /=2) THEN
                CALL error_bi('syntax error in #const or #vgrad declaration' &
                     ,substr)
             ENDIF
             ! save names for later check of availability from restart files
             s_inp(j1,j2)%cha = TRIM(sname)
             s_inp(j1,j2)%obj = TRIM(outstring(1))

             SELECT CASE(j1)
             CASE(:id_cfc12,id_o2)
                CALL new_channel_object(status, TRIM(sname) &
                     , TRIM(outstring(1)) &
                     , p3=xradin(j1,j2)%ptr3, reprid=GP_3D_MID)
             CASE(id_aot_lw)
                CALL get_representation_id(status, "REPR_AEROPT_4D_JPBAND" &
                     , repr_id)
                CALL new_channel_object(status &
                     , s_inp(j1, j2)%cha, s_inp(j1, j2)%obj, &
                     p4 = xradin(j1, j2)%ptr4, reprid = repr_id)
             CASE (id_aot_sw, id_gamma_sw, id_omega_sw)
                CALL get_representation_id(status, "REPR_AEROPT_4D_NSW" &
                     , repr_id)
                CALL new_channel_object(status &
                     , s_inp(j1, j2)%cha, s_inp(j1, j2)%obj, &
                     p4 = xradin(j1, j2)%ptr4, reprid = repr_id)
             CASE DEFAULT
                CALL error_bi(&
        '#const/#vgrad only implemented for CH4, N2O, F11, F12, aerosols, O2' &
            , substr)
             END SELECT
             IF (status == 3102) THEN
                CALL channel_halt(substr//': channel object '//&
                     &TRIM(outstring(1))//&
                     &' exists already in channel '//TRIM(sname), status)
             ENDIF
             SELECT CASE(j1)
             CASE(:id_cfc12,id_o2)
                CALL new_attribute(status, TRIM(sname), TRIM(outstring(1)) &
                     , 'units', c='mol/mol' )
                CALL channel_halt(substr, status)
             END SELECT

             DEALLOCATE(outstring) ; NULLIFY(outstring)

          ELSE

             CALL get_channel_object_info(status,TRIM(r_inp(j1,j2)%cha) &
                  , TRIM(r_inp(j1,j2)%obj), reprid=repr_id)
             CALL channel_halt(substr//': info for object '//&
                  &TRIM(r_inp(j1,j2)%obj)//' in channel '//&
                  TRIM(r_inp(j1,j2)%cha)//' not found!',status)

             CALL get_representation_info(status,inpname='',id=repr_id &
                  , rank=irank)
             CALL channel_halt(substr//': representation info for object '//&
                  &TRIM(r_inp(j1,j2)%obj)//' in channel '//&
                  TRIM(r_inp(j1,j2)%cha)//' not found!',status)

             which_rank: SELECT CASE(irank)
             CASE(0)
                CALL get_channel_object(status &
                     , TRIM(r_inp(j1,j2)%cha), TRIM(r_inp(j1,j2)%obj) &
                     , p0=xradin(j1,j2)%ptr0)

             CASE(2)
                CALL get_channel_object(status &
                     , TRIM(r_inp(j1,j2)%cha), TRIM(r_inp(j1,j2)%obj) &
                     , p2=xradin(j1,j2)%ptr2)

             CASE(3)
                CALL get_channel_object(status &
                     , TRIM(r_inp(j1,j2)%cha), TRIM(r_inp(j1,j2)%obj) &
                     , p3=xradin(j1,j2)%ptr3)

                CALL channel_halt(substr//': object '//&
                     &TRIM(r_inp(j1,j2)%obj)//' in channel '//&
                     TRIM(r_inp(j1,j2)%cha)//' not found!', status)

                ! CHECK UNIT
                CALL get_attribute(status  &
                     , TRIM(r_inp(j1,j2)%cha), TRIM(r_inp(j1,j2)%obj) &
                     , 'units', c=att_unit(j1,j2))
                !
                IF(status /= 0) THEN
                   SELECT CASE(j1)
                   CASE(id_h2o,id_co2,id_ch4,id_o3,id_n2o,id_cfc11,id_cfc12)
                      CALL channel_halt(substr//': object '//&
                           &TRIM(r_inp(j1,j2)%obj)//' in channel '//&
                           TRIM(r_inp(j1,j2)%cha)//&
                           &' has no attribute (check unit!)', status)
                   CASE DEFAULT
                       CALL warning_bi('channel object '//&
                              &TRIM(r_inp(j1,j2)%cha)//' - '//&
                              &TRIM(r_inp(j1,j2)%obj)//&
                              ' has no unit!', &
                              substr)
                       status=0
                   END SELECT
                ELSE
                   ! must have SI
                   SELECT CASE(j1)
                   CASE(id_co2,id_ch4,id_o3,id_n2o,id_cfc11,id_cfc12)
! This does currently not yet work, since NCREGRID->IMPORT_GRID does
! not deliver the correct unit ...
!!$                      IF(att_unit(j1,j2) /= 'mol/mol') then
!!$                         CALL error_bi('unknown unit of '//&
!!$                              &TRIM(r_inp(j1,j2)%cha)//' - '//&
!!$                              &TRIM(r_inp(j1,j2)%obj)//&
!!$                              ' distribution! Unit must be mol/mol!)', &
!!$                              substr)
!!$                      ENDIF
                      CALL warning_bi('unit check of '//&
                              &TRIM(r_inp(j1,j2)%cha)//' - '//&
                              &TRIM(r_inp(j1,j2)%obj)//&
                              ' not possible! Unit must be mol/mol!', &
                              substr)
                   ! special case for H2O
                   CASE(id_h2o)
                      IF ((att_unit(j1,j2) /= 'mol/mol') &
                           .AND. (att_unit(j1,j2) /= 'kg kg-1') &
                           .AND. (att_unit(j1,j2) /= 'kg/kg')) then
                         CALL channel_halt(substr//': unknown unit for '//&
                                 &TRIM(r_inp(j1,j2)%obj)//&
                                 ' H2O distribution!)', status)
                      END IF
                   CASE DEFAULT
                        CALL warning_bi('unit check of '//&
                              &TRIM(r_inp(j1,j2)%cha)//' - '//&
                              &TRIM(r_inp(j1,j2)%obj)//&
                              ' not implemented!', &
                              substr)
                        status=0
                   END SELECT
                END IF

             CASE(4)
                CALL get_channel_object(status &
                     , TRIM(r_inp(j1,j2)%cha), TRIM(r_inp(j1,j2)%obj) &
                     , p4=xradin(j1,j2)%ptr4)

             CASE DEFAULT
                  CALL error_bi('rank of channel object '//&
                       &TRIM(r_inp(j1,j2)%cha)//' - '//&
                       &TRIM(r_inp(j1,j2)%obj)//&
                       ' not supported!', &
                       substr)

             END SELECT which_rank

             CALL channel_halt(substr//': object '//&
                  &TRIM(r_inp(j1,j2)%obj)//' in channel '//&
                  TRIM(r_inp(j1,j2)%cha)//' not found!',status)

          ENDIF

          ! CHECK, IF POINTER IS ASSOCIATED
          l_assoc = ASSOCIATED(xradin(j1,j2)%ptr0) .OR. &
                    ASSOCIATED(xradin(j1,j2)%ptr2) .OR. &
                    ASSOCIATED(xradin(j1,j2)%ptr3) .OR. &
                    ASSOCIATED(xradin(j1,j2)%ptr4)
          IF (.NOT. l_assoc) THEN
             WRITE(msg,*) 'ERROR: pointer of input element ',j1, ' in ' &
                  , 'radiation call ',j2,' is not associated'
             CALL error_bi(TRIM(msg),substr)
          END IF

       END DO inp_objs

    END DO rad_calls

    IF (lfubrad) CALL rad_fubrad_init_coupling

    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output

  END SUBROUTINE rad_init_coupling
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE rad_global_start(flag)

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_data_bi,          ONLY: press_3d
    USE messy_main_grid_def_mem_bi,  ONLY: npromz, ngpblks
    USE messy_main_timer_bi,         ONLY: event_state
    USE messy_main_timer,            ONLY: current_date, time_days &
                                         , add_date, print_date_components&
                                         , lstart, lresume, delta_time

    USE messy_main_tools,            ONLY: str2num, strcrack
    USE messy_main_channel,          ONLY: get_channel_object_info
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_blather_bi,       ONLY: info_bi, warning_bi


    IMPLICIT NONE
    INTRINSIC :: INT, REAL, TANH, LOG

    ! I/O
    INTEGER, INTENT(IN)         :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'rad_global_start'
    INTEGER                     :: iradhlen
    TYPE (time_days)            :: rad_date
    CHARACTER(len=STRLEN_ULONG) :: date_mess1, date_mess2
    INTEGER                     :: j1, j2, m, status
    INTEGER                     :: TC_PRN_NATIVE
    REAL(DP)                              :: const_vmr
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER :: outstring(:) => NULL()
    REAL(DP)                    :: z_m, z_d, z_po, z_r, z_da
    LOGICAL                     :: linit
#ifdef CESM1
    INTEGER                     :: zjrow, i
#endif

  ! trigger rad, as it is done in time_set which is
  ! called in stepon and has been moved here

    IF (flag == 1) THEN
       ! called before ORBIT
       l_trigrad= event_state(ev_trigrad, current_date) .OR. lstart

       IF (l_trigrad) THEN
          rad_date = current_date
          iradhlen = INT(0.5*event_state(ev_trigrad,.TRUE.))
          dt_offset = REAL(iradhlen, DP)
          ! orbit parameters are calculated without time offset, if radiation is
          ! calculated every time step
          ! (i.e., if iradhlen is smaller than/equal to 0.5*delta_time)
          IF (dt_offset <= delta_time/2.0_dp) THEN
             dt_offset = 0.0_dp
          END IF
          CALL add_date(0,INT(dt_offset),rad_date)
          IF (.NOT. l_quiet) THEN
             CALL print_date_components(rad_date,TC_PRN_NATIVE, mess=date_mess1)
             date_mess2 = TRIM('Radiation (RAD) calculated for : ')//&
                  &' '// TRIM(date_mess1)
             date_mess1 = TRIM(date_mess2)
             IF (p_parallel_io) WRITE(*,*) TRIM(date_mess1)
          END IF
       ELSE
          dt_offset = -1.0_dp
       ENDIF

       ! SET CONSTANT MIXING RATIOS WITH OR WITHOUT VERTICAL GRADIENT
       IF (lstart .OR. lresume) THEN
          rad_calls: DO j2=1,NRADCALL
             IF (.NOT. l_switch(j2)) CYCLE
             inp_objs: DO j1=1,NMAXINP
                IF ( (TRIM(r_inp(j1,j2)%cha) .EQ. '#const') .OR. &
                     (TRIM(r_inp(j1,j2)%cha) .EQ. '#vgrad') ) THEN

                   CALL strcrack(TRIM(r_inp(j1,j2)%obj), '=', outstring, m)
                   CALL str2num(outstring(2), const_vmr, status)
                   IF (status /= 0) THEN
                      CALL error_bi(&
                           'syntax error in #const or #vgrad declaration' &
                           ,substr)
                   END IF
                   DEALLOCATE(outstring) ; NULLIFY(outstring)

                   IF (lresume) THEN
                      CALL get_channel_object_info(status &
                           , TRIM(s_inp(j1,j2)%cha), TRIM(s_inp(j1,j2)%obj) &
                           , lrestart_read=linit)
                      CALL channel_halt(substr, status)
                      IF (linit) THEN
                         CALL info_bi(&
                              'object '//TRIM(s_inp(j1,j2)%obj)//' in '&
                              &'channel '//TRIM(s_inp(j1,j2)%cha)//&
                              &' initialised from restart file ...' &
                              , substr)
                         CYCLE
                      ELSE
                         CALL warning_bi(&
                              'object '//TRIM(s_inp(j1,j2)%obj)//' in '&
                              &'channel '//TRIM(s_inp(j1,j2)%cha)//&
                              &' NOT initialised from restart file ...'&
                              &' (will be re-initialised)' &
                              , substr)
                      ENDIF
                   END IF

                   const_or_vgrad: SELECT CASE (TRIM(r_inp(j1,j2)%cha))
                   CASE ('#const')
                      SELECT CASE (j1)
                      CASE (id_aot_lw, id_aot_sw, id_gamma_sw, id_omega_sw)
                         xradin(j1,j2)%ptr4(:,:,:,:) = const_vmr
                      CASE DEFAULT
                         xradin(j1,j2)%ptr3(:,:,:) = const_vmr
                      END SELECT
                   CASE ('#vgrad')
                      SELECT CASE(j1)
                      CASE(id_ch4)
                         z_r  = 1.25E-01_dp
                         z_po = 683._dp
                         z_da = -1.43_dp
                      CASE(id_n2o)
                         z_r  = 1.2E-02_dp
                         z_po = 1395._dp
                         z_da = -1.43_dp
                      CASE(id_cfc11)
                         z_r  = 1.E-04_dp
                         z_po = 4159._dp
                         z_da = -0.73_dp
                      CASE(id_cfc12)
                         z_r  = 1.E-04_dp
                         z_po = 3177.4_dp
                         z_da = -0.73_dp
                      CASE DEFAULT
                         CALL error_bi(&
                              '#vgrad only implemented for CH4, N2O, F11, F12' &
                              ,substr)
                      END SELECT
                      z_m  = (const_vmr + z_r * const_vmr) * 0.5_dp
                      z_d  = (const_vmr - z_r * const_vmr) * 0.5_dp
#if defined(ECHAM5) || defined (MBM_RAD)
                      xradin(j1,j2)%ptr3(_RI_XYZ__(:,1:ngpblks-1,:)) = &
                           (1._dp - (z_d/z_m) &
                           * tanh(log(press_3d(_RI_XYZ__(:,1:ngpblks-1,:))/z_po)    &
                           /z_da))*z_m
                      xradin(j1,j2)%ptr3(_RI_XYZ__(1:npromz,ngpblks,:)) = &
                           (1._dp - (z_d/z_m) &
                           * tanh(log(press_3d(_RI_XYZ__(1:npromz,ngpblks,:))/z_po) &
                           /z_da))*z_m
#endif
#ifdef CESM1
                      DO zjrow=1,ngpblks
                         DO i=1,npromz(zjrow)
                            xradin(j1,j2)%ptr3(i,:,zjrow) = &
                                 (1._dp - (z_d/z_m) &
                                 * tanh(log(press_3d(i,:,zjrow)/z_po)    &
                                 /z_da))*z_m
                         END DO
                      END DO
#endif
                   END SELECT const_or_vgrad


                END IF
             END DO inp_objs
          END DO rad_calls
       END IF

       RETURN ! retrun here, if called with flag=1
    END IF

    ! called after ORBIT (with flag==2)

    ! save cdisse for use in FUBRAD; update total solar irradiance (at 1 AU)
    IF (lfubrad) CALL rad_fubrad_global_start

  END SUBROUTINE rad_global_start
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE rad_radiation

    ! BMIL
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, kproma, nlev,nlevp1,jrow
    USE messy_main_data_bi,          ONLY: &
         loland_2d,loglac_2d,forest,seaice,&
         cvs,sni,cvsc,vlt,                 &
         slm,icecov,seacov,                &
         tslm1,tsi,tsw,alb,                &
         aphm1,apm1,tm1,                   &
         srfl_2d, press_3d

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_timer,            ONLY: delta_time, lstart
    USE messy_rad_fubrad_si,         ONLY: solc_fubrad

    ! SMCL
    USE messy_main_constants_mem,    ONLY: solc,crae &
                                         , M_O3, MC, MO, MH, MN, M_air, M_H2O &
                                         , OneDay,vtmpc1, M_O2
    USE messy_main_tools,            ONLY: tlucua, jptlucu1, jptlucu2 &
                                         , iso2ind
    USE messy_rad,                   ONLY: rad_rad_smcl, nsw, jpband &
                                         , albedos_rad

    IMPLICIT NONE
    INTRINSIC :: EPSILON, SPREAD

    CHARACTER(LEN=*), PARAMETER :: substr = 'rad_radiation'
    INTEGER                            :: jl,jk,j2
    REAL(dp)                           :: zsct, zsctm     !solar irrad.
    REAL(dp)                           :: zcrae
    REAL(dp), DIMENSION(nproma)        :: zcosszac !zenith angle

    REAL(dp), DIMENSION(nproma,nlev)   :: zdp
    REAL(dp), DIMENSION(nproma,nlev)   :: ztf
    REAL(dp), DIMENSION(nproma,nlevp1) :: zth
    !
    REAL(dp),PARAMETER                 :: zh2ovmrtokgkg=M_H2O/M_air
    REAL(dp),PARAMETER                 :: zco2vmrtokgkg=(MC+2._dp*MO)/M_air
    REAL(dp),PARAMETER                 :: zo3vmrtokgkg=M_O3/M_air
    REAL(dp),PARAMETER                 :: zo2vmrtokgkg=M_O2/M_air
    REAL(dp),PARAMETER                 :: zch4vmrtokgkg=(MC+4._dp*MH)/M_air
    REAL(dp),PARAMETER                 :: zn2ovmrtokgkg=(2._dp*MN+MO)/M_air
    !
    REAL(dp), DIMENSION(nproma,nlev)   :: zq
    REAL(dp), DIMENSION(nproma,nlev)   :: zqs
    REAL(dp), DIMENSION(nproma,nlev)   :: zclc
    REAL(dp), DIMENSION(nproma,nlev)   :: zco2
    REAL(dp), DIMENSION(nproma,nlev)   :: zo3
    REAL(dp), DIMENSION(nproma,nlev)   :: zo2
    REAL(dp), DIMENSION(nproma,nlev)   :: zch4
    REAL(dp), DIMENSION(nproma,nlev)   :: zn2o
    REAL(dp), DIMENSION(nproma,nlev,2) :: zcfcs
    !
    REAL(dp)::dtemp
    INTEGER, DIMENSION(nproma) :: ztp_fixed_i
    !
    INTEGER                      :: it  !temperature for lookup table entry
    LOGICAL                      :: lookupoverflow=.FALSE.
    !
    INTEGER                     :: status

    ! Intrinsic functions
    INTRINSIC :: INT, MAX, MIN, SQRT

    ! zsct calculation moved before l_trigrad,
    ! since it is needed for all time steps.
    ! zsctm is used within l_trigrad below for consistency
    !
    ! solar insolation = solar irradiation (at 1 AU) at !
    ! current orbit position
    IF (lfubrad) THEN
       zsct=cdisse_x*solc_fubrad
       zsctm=cdissem_x*solc_fubrad
    ELSE
       zsct=cdisse_x*solc
       IF (rset_solc%l) zsct=cdisse_x * rset_solc%v
       zsctm=cdissem_x*solc
       IF (rset_solc%l) zsctm=cdissem_x * rset_solc%v
    ENDIF

    zi0_2d(1:kproma,jrow) = zsct*cosszac_x(1:kproma,jrow) * &
         rdayl_x(1:kproma,jrow)

    !
    trig_rad: IF (l_trigrad) THEN

       ! 0) Initialize for kproma < nproma
       zq(:,:) = 0.0_dp
       zqs(:,:) = 0.0_dp
       zclc(:,:) = 0.0_dp
       zco2(:,:) = 0.0_dp
       zo3(:,:) = 0.0_dp
       zo2(:,:) = 0.0_dp
       zch4(:,:) = 0.0_dp
       zn2o(:,:) = 0.0_dp
       zcfcs(:,:,:) = 0.0_dp
       zth(:,:) = 0.0_dp
       ztf(:,:) = 0.0_dp
       zcosszac(:) = 0.0_dp

       ! 1) Solar insolation
       ! zsct is calculated outside the IF (l_trigrad) condition
       zcrae = crae*(crae+2._dp)
       zcosszac(1:kproma)  = crae/(SQRT(cosszacm_x(1:kproma,jrow)**2+zcrae) &
            - cosszacm_x(1:kproma,jrow))

       ! 2) pre-calculations

       ! layer pressure thickness
       zdp(:,:) = 0.0_dp ! op_sd_20130405 init for kproma < nproma
       zdp(1:kproma,:)=aphm1(1:kproma,2:nlev+1)-aphm1(1:kproma,1:nlev)

       IF (ANY(i_rad(:) == 2)) THEN
          DO jk=1,nlev
             DO  jl=1,kproma
                !save temp. one time step before
                padutm1o(_RI_XYZ__(jl,jrow,jk))=padutm1(_RI_XYZ__(jl,jrow,jk))
                padutm1(_RI_XYZ__(jl,jrow,jk))=tm1(_RI_XYZ__(jl,jrow,jk))
             END DO
          END DO

          IF(lstart)THEN
           padutm1o(_RI_XYZ__(:,jrow,:))=0.0_dp
          ENDIF

       END IF

       rad_calls: DO j2=1,NRADCALL

          IF (.NOT. l_switch(j2)) CYCLE

                 call_type: SELECT CASE(i_rad(j2))

                 CASE(1)  ! instantaneous forcing

                    xradout(ido_addst,j2)%ptr3(_RI_XYZ__(:,jrow,:)) = tm1(_RI_XYZ__(:,jrow,:))

                 CASE(2) ! stratospheric temperature adjustment

                    DO jk=1,nlev
                       DO jl = 1, kproma

                          ! dtemp: fixed dynamical heating rate dT/dt(dyn)
                          !        dT/dt(dyn)= dT/dt-dT/dt(rad)

                          dtemp=(1./delta_time) * (tm1(_RI_XYZ__(jl,jrow,jk)) - &
                               padutm1o(_RI_XYZ__(jl,jrow,jk))) &
                               - xradout(ido_tte,1)%ptr3(_RI_XYZ__(jl,jrow,jk))

                          ! addst: pert. temperature T*
                          ! dT*/dt=dT/dt(dyn)+dT*/dt(rad) ->T*(=addst)
                          xradout(ido_addst,j2)%ptr3(_RI_XYZ__(jl,jrow,jk)) =      &
                               xradout(ido_addst,j2)%ptr3(_RI_XYZ__(jl,jrow,jk)) + &
                               delta_time * (dtemp + &
                               xradout(ido_tte,j2)%ptr3(_RI_XYZ__(jl,jrow,jk)))

                          ! reset to original temp. below trop.
                          IF (aphm1(jl,jk) .GE. tp_fixed_x(jl,jrow)    .AND. &
                              aphm1(jl,jk+1) .GT. tp_fixed_x(jl,jrow)) THEN
                              xradout(ido_addst,j2)%ptr3(_RI_XYZ__(jl,jrow,jk)) = &
                                   tm1(_RI_XYZ__(jl,jrow,jk))
                          END IF

                       END DO
                    END DO

                    CALL iso2ind(kproma,                                   &
                         press_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)),        &
                         tp_fixed_x(1:kproma,jrow), ztp_fixed_i(1:kproma), &
                         lrev=.true.)
                    tp_fixed_i(1:kproma,jrow) = REAL(ztp_fixed_i(1:kproma),dp)

                 END SELECT call_type

       END DO rad_calls

       ! saturation specific humidity
       ! EPSILON(1.) serves to avoid saturated water vapour
       ! content in a layer of less than 2*EPSILON(1.)
       DO jk = 1, nlev
          DO jl = 1, kproma
             it = INT(tm1(_RI_XYZ__(jl,jrow,jk))*1000._dp)
             IF ( (it<jptlucu1 .OR. it>jptlucu2) .AND. &
                  (apm1(jl,jk) >= 1.0_dp) ) lookupoverflow = .TRUE.
             it = MAX(MIN(it,jptlucu2),jptlucu1)
             zqs(jl,jk) = tlucua(it)/apm1(jl,jk)
          END DO
       END DO
       IF (lookupoverflow) &
            CALL error_bi('lookuperror',substr)
       zqs(1:kproma,:)= MIN(zqs(1:kproma,:),0.5_dp)
       zqs(1:kproma,:)= zqs(1:kproma,:)/(1._dp-vtmpc1*zqs(1:kproma,:))
       zqs(1:kproma,:)= MAX(2._dp*EPSILON(1._dp),zqs(1:kproma,:))

       ! 3) surface albedo
       CALL albedos_rad(nproma, kproma,                 &       ! In
            loland_2d(:,jrow), loglac_2d(:,jrow),       &
            forest(:,jrow),seaice(:,jrow),              &
            cvs(:,jrow),sni(:,jrow),                    &
            cvsc(:,jrow),vlt(:,jrow),                   &
            slm(:,jrow),seacov(:,jrow),icecov(:,jrow),  &
            tslm1(:,jrow),tsi(:,jrow),alb(:,jrow),      &
            alsol_x(:,jrow),alsow_x(:,jrow),alsoi_x(:,jrow), &  ! Out
            albedo_x(:,jrow))

       ! 4) input via name list -- changeable

       rad_calls2: DO j2=1,NRADCALL

          IF (.NOT. l_switch(j2)) CYCLE

          ! THIS CALCULATION IS INDEPENDENT OF j2 ...
          ! BUT DEPENDENT ON ztf(j2) and therefore must remain here!

          ztf(:,:)=xradout(ido_addst,j2)%ptr3(_RI_XYZ__(:,jrow,:))

          DO jk=2,nlev
             DO jl = 1, kproma
                zth(jl,jk)= ( ztf(jl,jk-1) * apm1(jl,jk-1) * &
                     (apm1(jl,jk)-aphm1(jl,jk)) +   &
                     ztf(jl,jk)   * apm1(jl,jk)  * &
                     (aphm1(jl,jk)-apm1(jl,jk-1))   &
                     ) &
                     / ( aphm1(jl,jk) * &
                     (apm1(jl,jk)-apm1(jl,jk-1)) )
             END DO
          END DO
          !
          DO jl = 1, kproma
             ! lowermost level
             zth(jl,nlevp1)= &
                  ( slm(jl,jrow)    * tslm1(jl,jrow)**4   &
                  + icecov(jl,jrow) * tsi(jl,jrow)**4     &
                  + seacov(jl,jrow) * tsw(jl,jrow)**4     &
                  )**0.25_dp
             ! uppermost level
             zth(jl,1)=ztf(jl,1)-apm1(jl,1)*(ztf(jl,1) &
                  -zth(jl,2)) / (apm1(jl,1)-aphm1(jl,2))
          END DO

          ! COPY TO CHANNEL OBJECT FOR OUTPUT
          xradout(ido_addsth,j2)%ptr3(_RI_XYZ__(:,jrow,:)) = zth(:,:)

          ! - avoid water vapour content in a layer less than EPSILON(1.)

          ! - query units, if vmr convert to kgkg-1
          SELECT CASE(att_unit(id_h2o,j2))
             ! H2O vmr as input: convert to specific humidity
             ! (mass_H2O / mass_moistair)
          CASE('mol/mol')
             zq(1:kproma,:) = MAX(    &
                  (xradin(id_h2o,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) &
                  * zh2ovmrtokgkg     &
                  / ( 1 + xradin(id_h2o,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) &
                  * zh2ovmrtokgkg )), &
                  EPSILON(1._dp))
             !specific humidity as input
          CASE('kg kg-1', 'kg/kg')
             zq(1:kproma,:) = MAX( &
                  xradin(id_h2o,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)), &
                  EPSILON(1._dp))
          CASE DEFAULT
             CALL channel_halt(substr//': unknown unit for '//&
                  &TRIM(r_inp(id_h2o,j2)%obj)//&
                  ' H2O distribution!)', status)
          END SELECT

          zclc(1:kproma,:)=xradin(id_clc,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:))
          zco2(1:kproma,:)=xradin(id_co2,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) &
               *zco2vmrtokgkg
          zo3(1:kproma,:)=xradin(id_o3,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:))   &
               *zo3vmrtokgkg
          zo2(1:kproma,:)=xradin(id_o2,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:))   &
               *zo2vmrtokgkg
          zch4(1:kproma,:)=xradin(id_ch4,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) &
               *zch4vmrtokgkg
          zn2o(1:kproma,:)=xradin(id_n2o,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) &
               *zn2ovmrtokgkg

          ! NOTE: CFCs need to be in vmr not in mmr
          ! (see comment in rad_int.f90 of ECHAM5)
          zcfcs(1:kproma,:,1) = &
               xradin(id_cfc11,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:))
          zcfcs(1:kproma,:,2) = &
               xradin(id_cfc12,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:))

          ! this really required foreach NRADCALL, due to zo3 ...
          IF (lfubrad) &
               CALL rad_fubrad_preprad(kproma,nlev,jrow,zo3,zo2,apm1,aphm1 &
               ,zcosszac,j2)

          CALL rad_rad_smcl(kproma,nproma,nlev,    &
               & i_sw(j2),                         &
               & INT(OneDay),                      &
               & zsctm,zcosszac,albedo_x(:,jrow),   &
               & apm1(:,:),aphm1(:,:),             &
               & ztf(:,:),zth(:,:),                &
               & zq(:,:),zqs(:,:),zclc(:,:),       &
               & xradin(id_cld_clcv,j2)%ptr2(:,jrow),   &
               & xradin(id_cld_cld, j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
               & zo3(:,:),zco2(:,:),zch4(:,:),zn2o(:,:),zcfcs(:,:,:), &
               & xradin(id_aot_lw,j2)%ptr4(_RI_XYZN_(:,jrow,:,1:jpband)),  &
               & xradin(id_aot_sw,j2)%ptr4(_RI_XYZN_(:,jrow,:,1:nsw)),     &
               & xradin(id_omega_sw,j2)%ptr4(_RI_XYZN_(:,jrow,:,1:nsw)),   &
               & xradin(id_gamma_sw,j2)%ptr4(_RI_XYZN_(:,jrow,:,1:nsw)),   &
               & xradin(id_cld_lw,j2)%ptr4(_RI_XYZN_(:,jrow,:,1:jpband)),  &
               & xradin(id_cld_sw,j2)%ptr4(_RI_XYZN_(:,jrow,:,1:nsw)),     &
               & xradin(id_cld_gamma,j2)%ptr4(_RI_XYZN_(:,jrow,:,1:nsw)),  &
               & xradin(id_cld_omega,j2)%ptr4(_RI_XYZN_(:,jrow,:,1:nsw)),  &
               ! OUTPUT ---------
               ! SW output has changed to transmissivity, therefore
               ! xradout arrays can be used directly, and division by
               ! zsct*zcosszac is no longer necessary.
               ! Transmissivities directed upward are added.
               ! ido_emtef and ido_emter are obsolete, as the output
               ! of rad_rad_smcl is the net LW flux
               & xradout(ido_flxt,j2)%ptr3(_RI_XYZ__(:,jrow,:)),  & ! net flux LW
               & xradout(ido_trsol,j2)%ptr3(_RI_XYZ__(:,jrow,:)), & ! transm. SW (all)
               & xradout(ido_flxtf,j2)%ptr3(_RI_XYZ__(:,jrow,:)),  & ! net flux LW cl.sky
               & xradout(ido_trsof,j2)%ptr3(_RI_XYZ__(:,jrow,:)),  & ! transm. SW (all) cl.sky
               & xradout(ido_trnir,j2)%ptr3(_RI_XYZ__(:,jrow,:)),  & ! transm. NIR
               & xradout(ido_trnif,j2)%ptr3(_RI_XYZ__(:,jrow,:)),  & ! transm. NIR cl.sky
               & xradout(ido_trsw1,j2)%ptr3(_RI_XYZ__(:,jrow,:)),  & ! transm. UVVIS
               & xradout(ido_trs1f,j2)%ptr3(_RI_XYZ__(:,jrow,:)),  & ! transm. UVVIS cl.sky
               & xradout(ido_flxut, j2)%ptr3(_RI_XYZ__(:,jrow,:)), & ! flux LW up
               & xradout(ido_trus,  j2)%ptr3(_RI_XYZ__(:,jrow,:)), & ! transm. SW up
               & xradout(ido_flxutf,j2)%ptr3(_RI_XYZ__(:,jrow,:)), & ! flux LW up cl.sky
               & xradout(ido_trusf, j2)%ptr3(_RI_XYZ__(:,jrow,:)), & ! transm. SW up cl.sky
               & xradout(ido_truni, j2)%ptr3(_RI_XYZ__(:,jrow,:)), & ! transm. NIR up
               & xradout(ido_trunif,j2)%ptr3(_RI_XYZ__(:,jrow,:))  & ! transm. NIR up cl.sky
              &)


       ENDDO rad_calls2

    END IF trig_rad

    ! NOTE: srfl must be calculated here, since it is required for
    !       vdiff, which is called from physc after radiation and before
    !       radheat!

    ! Note: first radiation call is defining the feedback => 1
    IF (l_feedback) THEN
       ! required in physc for vdiff (origin: DATA)
       srfl_2d(1:kproma,jrow)= &
            zi0_2d(1:kproma,jrow) * &
            xradout(ido_trsol,1)%ptr3(_RI_XYZ__(1:kproma,jrow,nlevp1))
    END IF

    ! additional diagnostics ... (requires zi0_2d from above!)
    rad_calls3: DO j2=1,NRADCALL

       IF (.NOT. l_switch(j2)) CYCLE

       xradout(ido_flxus,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) =   &
          xradout(ido_trus,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) * &
          SPREAD(zi0_2d(1:kproma,jrow),2,nlevp1)

       xradout(ido_flxusf,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) =   &
          xradout(ido_trusf,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) * &
          SPREAD(zi0_2d(1:kproma,jrow),2,nlevp1)

       xradout(ido_flxuni,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) =   &
          xradout(ido_truni,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) * &
          SPREAD(zi0_2d(1:kproma,jrow),2,nlevp1)

       xradout(ido_flxunif,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) =   &
          xradout(ido_trunif,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) * &
          SPREAD(zi0_2d(1:kproma,jrow),2,nlevp1)

       ! Save the fluxes in the 3 NIR bands and the UVvis band of the rad bands.
       xradout(ido_flxnir,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) =   &
          xradout(ido_trnir,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) * &
          SPREAD(zi0_2d(1:kproma,jrow),2,nlevp1)

       xradout(ido_flxsw1,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) =   &
          xradout(ido_trsw1,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:)) * &
          SPREAD(zi0_2d(1:kproma,jrow),2,nlevp1)

    ENDDO rad_calls3

  END SUBROUTINE rad_radiation
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE rad_radheat
    ! This routine is called from messy_radheat in messy_main_control

    USE messy_main_grid_def_mem_bi, ONLY: jrow, nlev, nlevp1 &
                                     , nproma, kproma
    USE messy_main_data_bi,      ONLY: aphm1,apm1,tm1,qm1 &
#ifndef MESSYTENDENCY
                                     , tte_3d             &
#endif
                                     , slm,icecov,seacov  &
                                     , tslm1,tsi,tsw      &
                                     , tslnew             &
                                     , radflxw_2d         &
                                     , ledith

    ! SMCL
    USE messy_main_constants_mem,   ONLY: M_O3, M_O2, M_air
    USE messy_rad,                  ONLY: rad_radheat_smcl

    IMPLICIT NONE

    ! local
    INTEGER :: j2
    REAL(dp), DIMENSION(nproma,nlev)   :: zo3  ! for fubrad
    REAL(dp), DIMENSION(nproma,nlev)   :: zo2  ! for fubrad
    REAL(dp),PARAMETER                 :: zo3vmrtokgkg=M_O3/M_air
    REAL(dp),PARAMETER                 :: zo2vmrtokgkg=M_O2/M_air
    REAL(dp), DIMENSION(nproma, nlev)  :: tte_fub

    DO j2=1,NRADCALL
       IF (.NOT. l_switch(j2)) CYCLE

       ! Note: tte is INTENT(OUT) (=> radiative temperature tendency)
       CALL rad_radheat_smcl (kproma, nproma, nlev,  nlevp1, zi0_2d(:,jrow) &
            , tm1(_RI_XYZ__(:,jrow,:))                        &
            , qm1(_RI_XYZ__(:,jrow,:))                        &
            , xradout(ido_trsof,j2)%ptr3(_RI_XYZ__(:,jrow,:)) &
            , xradout(ido_trsol,j2)%ptr3(_RI_XYZ__(:,jrow,:)) &
            , xradout(ido_trsw1,j2)%ptr3(_RI_XYZ__(:,jrow,:)) &
            , xradout(ido_trnir,j2)%ptr3(_RI_XYZ__(:,jrow,:)) &
            , xradout(ido_trs1f,j2)%ptr3(_RI_XYZ__(:,jrow,:)) &
            , xradout(ido_trnif,j2)%ptr3(_RI_XYZ__(:,jrow,:)) &
            , xradout(ido_trfll,j2)%ptr2(:,jrow)   &
            , xradout(ido_trflw,j2)%ptr2(:,jrow)   &
            , xradout(ido_trfli,j2)%ptr2(:,jrow)   &
            , xradout(ido_sofll,j2)%ptr2(:,jrow)   &
            , xradout(ido_soflw,j2)%ptr2(:,jrow)   &
            , xradout(ido_sofli,j2)%ptr2(:,jrow)   &
            , xradout(ido_trfllac,j2)%ptr2(:,jrow) &
            , xradout(ido_trflwac,j2)%ptr2(:,jrow) &
            , xradout(ido_trfliac,j2)%ptr2(:,jrow) &
            , xradout(ido_sofllac,j2)%ptr2(:,jrow) &
            , xradout(ido_soflwac,j2)%ptr2(:,jrow) &
            , xradout(ido_sofliac,j2)%ptr2(:,jrow) &
            , xradout(ido_srad0u,j2)%ptr2(:,jrow)  &
            , xradout(ido_sradsu,j2)%ptr2(:,jrow)  &
            , xradout(ido_tradsu,j2)%ptr2(:,jrow)  &
            , tslm1(:,jrow), tsi(:,jrow), tsw(:,jrow), albedo_x(:,jrow) &
            , alsol_x(:,jrow), alsow_x(:,jrow), alsoi_x(:,jrow) &
            , aphm1(:,:), apm1(:,:), tslnew(:,jrow)             &
            , xradout(ido_tte,j2)%ptr3(_RI_XYZ__(:,jrow,:))     &
            , slm(:,jrow), seacov(:,jrow), icecov(:,jrow)       &
            , xradout(ido_flxs,j2)%ptr3(_RI_XYZ__(:,jrow,:))    &
            , xradout(ido_flxt,j2)%ptr3(_RI_XYZ__(:,jrow,:))    &
            , xradout(ido_flxsf,j2)%ptr3(_RI_XYZ__(:,jrow,:))   &
            , xradout(ido_flxtf,j2)%ptr3(_RI_XYZ__(:,jrow,:))   &
            , xradout(ido_dtsw,j2)%ptr3(_RI_XYZ__(:,jrow,:))    &
            , xradout(ido_dtlw,j2)%ptr3(_RI_XYZ__(:,jrow,:))    &
            , xradout(ido_flxsw1,j2)%ptr3(_RI_XYZ__(:,jrow,:))  &
            , xradout(ido_flxnir,j2)%ptr3(_RI_XYZ__(:,jrow,:))  &
            , xradout(ido_heats1,j2)%ptr3(_RI_XYZ__(:,jrow,:))  &
            , xradout(ido_heatni,j2)%ptr3(_RI_XYZ__(:,jrow,:))  &
            , xradout(ido_dtswc,j2)%ptr3(_RI_XYZ__(:,jrow,:))   & ! SW heating rates clear sky
            , xradout(ido_dtlwc,j2)%ptr3(_RI_XYZ__(:,jrow,:))   & ! LW heating rates clear sky
            , ledith                                  &
            )

       ! save the longwave upward flux of the lowest layer
       xradout(ido_flxut, j2)%ptr3(_RI_XYZ__(:,jrow,nlevp1)) = xradout(ido_tradsu,j2)%ptr2(:,jrow)
       xradout(ido_flxutf,j2)%ptr3(_RI_XYZ__(:,jrow,nlevp1)) = xradout(ido_tradsu,j2)%ptr2(:,jrow)

    END DO

    fubrad: IF (lfubrad) THEN

       zo3(:,:) = 0.0_dp

       DO j2=1,NRADCALL
          IF (.NOT. l_switch(j2)) CYCLE

          zo3(1:kproma,:)=xradin(id_o3,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:))*zo3vmrtokgkg
          zo2(1:kproma,:)=xradin(id_o2,j2)%ptr3(_RI_XYZ__(1:kproma,jrow,:))*zo2vmrtokgkg
          CALL rad_fubrad_preprad(kproma, jrow, j2)

          ! Note: tte is INTENT(OUT) (=> radiative temperature tendency)
          CALL rad_fubrad_radheat(kproma, nproma, jrow, l_trigrad, j2,  &
               nlevp1, zi0_2d(:,jrow),                   &
               zo3(:,:),                                 &
               zo2(:,:),                                 &
               xradout(ido_dtsw,j2)%ptr3(_RI_XYZ__(:,jrow,:)),  &
               xradout(ido_dtswc,j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
               rdayl_x(:,jrow),                          &
               tm1(_RI_XYZ__(:,jrow,:)),                 &
               tte_fub(:,:),                             &
               xradout(ido_trnir,j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
               xradout(ido_trnif,j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
               xradout(ido_trsol,j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
               xradout(ido_trsof,j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
               xradout(ido_flxs,j2)%ptr3(_RI_XYZ__(:,jrow,:)),  &
               xradout(ido_srad0u,j2)%ptr2(:,jrow),      &
               xradout(ido_flxsf,j2)%ptr3(_RI_XYZ__(:,jrow,:)),   &
               xradout(ido_flxnir,j2)%ptr3(_RI_XYZ__(:,jrow,:)),  &
               xradout(ido_flxuni,j2)%ptr3(_RI_XYZ__(:,jrow,:)),  &
               xradout(ido_flxunif,j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
               xradout(ido_flxus,j2)%ptr3(_RI_XYZ__(:,jrow,:)),   &
               xradout(ido_flxusf,j2)%ptr3(_RI_XYZ__(:,jrow,:)),  &
               xradout(ido_truni,j2)%ptr3(_RI_XYZ__(:,jrow,:)),   &
               xradout(ido_trunif,j2)%ptr3(_RI_XYZ__(:,jrow,:)),  &
               xradout(ido_trus,j2)%ptr3(_RI_XYZ__(:,jrow,:)),    &
               xradout(ido_trusf,j2)%ptr3(_RI_XYZ__(:,jrow,:))    &
               )

          ! add FUBRAD tendency to total radiation tendency
          xradout(ido_tte,j2)%ptr3(_RI_XYZ__(:,jrow,:)) = &
               xradout(ido_tte,j2)%ptr3(_RI_XYZ__(:,jrow,:)) + tte_fub(:,:)

       ENDDO

    END IF fubrad

    IF (l_feedback) THEN
       ! Note: radheat call with 1st radiation setup delivers temperature
       ! tendency for feedback; add tendency here
#ifdef MESSYTENDENCY
       !tendency budget
       call mtend_add_l (my_handle, mtend_id_t &
            , px  = xradout(ido_tte,1)%ptr3(_RI_XYZ__(:,jrow,:)))
#else
       tte_3d(_RI_XYZ__(:,jrow,:)) = tte_3d(_RI_XYZ__(:,jrow,:)) + &
         xradout(ido_tte,1)%ptr3(_RI_XYZ__(:,jrow,:))
#endif
       ! save this sum for ocean coupling (see physc.f90)
       radflxw_2d(:,jrow) = xradout(ido_soflw,1)%ptr2(:,jrow) + &
                            xradout(ido_trflw,1)%ptr2(:,jrow)
    END IF

  END SUBROUTINE rad_radheat
  ! ====================================================================

  ! ---------------------------------------------------------------------
  ! ---------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  !----------------------------------------------------------------------
  ! ----------------------------------------------------------------------

  ! ====================================================================
  SUBROUTINE rad_read_nml_cpl(status, iou)

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'rad_read_nml_cpl'

    NAMELIST /CPL/ l_feedback, lfubrad, trigrad, rset_solc    &
         , rset_calbmns, rset_calbmxs, rset_calbmni           &
         , tp_fixed, rsun_scale, l_switch, i_rad, r_inp, i_sw &
         , i_lw, l_quiet

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    ! FOR BACKWARD COMPATIBILITY, SET DEFAULT
    r_inp(id_o2,:)%cha = '#const'
    r_inp(id_o2,:)%obj = 'O2=0.20955775081737177' ! O2  [mol/mol]

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    WRITE(*,*) 'rsun_scale:   ', rsun_scale
    WRITE(*,*) 'rset_solc:    ', rset_solc

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE rad_read_nml_cpl
  ! ====================================================================

#if defined(MBM_RAD)
  ! ====================================================================
  SUBROUTINE rad_read_nml_cpl_mbm(status, iou)

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'rad_read_nml_cpl_mbm'

    NAMELIST /CPL_MBM/  ladd_tte, lv_echam, lzonal_mean, nn, nlev

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    ! NOTE: already at definition

    CALL read_nml_open(lex, substr, iou, 'CPL_MBM', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL_MBM, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL_MBM', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    WRITE(*,*) substr,' ladd_tte: ',ladd_tte,' lv_echam: ',lv_echam,' lzonal_mean: ',lzonal_mean
    WRITE(*,*) substr,' nn = ',nn,' nlev = ',nlev

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE rad_read_nml_cpl_mbm
  ! ====================================================================
#endif

#endif
! ***********************************************************************
END MODULE MESSY_RAD_SI
! ***********************************************************************
