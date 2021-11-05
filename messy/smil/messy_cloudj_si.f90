#include "messy_main_ppd_bi.inc"

!*****************************************************************************

! define ozone climatology as channel objects in ECHAM
! to calculate and use ozone climatology in COSMO instead of ECHAM ozone
! composite
!#ifdef ECHAM5
!#define OZONECLIM
!#endif

! submodel maintainer: Rolf Sander
! major contributions by Patrick Joeckel and Astrid Kerkweg
! see CHANGELOG for a detailed list of modifications

!*****************************************************************************

MODULE messy_cloudj_si

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      info_bi, error_bi, &
                                      warning_bi
  USE messy_main_channel,       ONLY: t_chaobj_cpl
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY, PTR_3D_ARRAY, str
  USE messy_main_constants_mem, ONLY: sp, dp, g, N_A, M_air, STRLEN_MEDIUM
  USE messy_cmn_photol_mem      ! IP_MAX, ip_*, jname
  USE messy_cloudj
  USE messy_cloudj_fjx_init_mod, ONLY: NAMFIL_FJX_spec, NAMFIL_FJX_scat_cld, &
    NAMFIL_FJX_scat_aer, NAMFIL_FJX_scat_UMa, NAMFIL_atmos_std, NAMFIL_FJX_j2j

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: cloudj_initialize
  PUBLIC :: cloudj_init_memory
  PUBLIC :: cloudj_init_coupling
  PUBLIC :: cloudj_global_start
  PUBLIC :: cloudj_physc
  PUBLIC :: cloudj_global_end
  PUBLIC :: cloudj_free_memory
  ! PRIVATE :: cloudj_read_nml_cpl

  INTRINSIC :: NULL

  ! CPL namelist
  ! - name of ozone tracer/channel object
  TYPE(t_chaobj_cpl) :: cloudj_O3
!!$  ! - name of solar cycle data time series channel/object
!!$  TYPE(t_chaobj_cpl) :: cloudj_solar
  ! - name of input data v3_h, o3_h, press_h channel/object
  TYPE(t_chaobj_cpl) :: cloudj_v3h
  TYPE(t_chaobj_cpl) :: cloudj_o3h
  TYPE(t_chaobj_cpl) :: cloudj_pressh
  ! - name of distance earth-sun channel/object
  !TYPE(t_chaobj_cpl) :: cloudj_cdisse
  ! - name of object with cos(zenith angle)
  TYPE(t_chaobj_cpl) :: cloudj_cossza
  !
  TYPE(t_chaobj_cpl) :: cloudj_alb
  ! - switch for skipping calculation of Lagrangian rate coefficients
  LOGICAL            :: l_skip_lg = .FALSE.
  ! - switch to force calculation for all species (independent of tracers)
  LOGICAL            :: l_force = .FALSE.

#if defined(ECHAM5) || defined(CESM1) || defined(COSMO)
  REAL(dp), DIMENSION(:,:), POINTER  :: albedo   => NULL()
#endif

!!$  TYPE(t_chaobj_cpl)           :: jv_aer_sca
!!$  TYPE(t_chaobj_cpl)           :: jv_aer_abs
!!$  TYPE(t_chaobj_cpl)           :: jv_aer_ga
!!$  REAL(dp), DIMENSION(:,:,:,:), POINTER :: aer_sca         => NULL()
!!$  REAL(dp), DIMENSION(:,:,:,:), POINTER :: aer_abs         => NULL()
!!$  REAL(dp), DIMENSION(:,:,:,:), POINTER :: aer_ga          => NULL()
!!$
!!$  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: jv_asca_5d    => NULL()
!!$  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: jv_aabs_5d    => NULL()
!!$  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: jv_ga_5d      => NULL()
!!$  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: jv_asca       => NULL()
!!$  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: jv_aabs       => NULL()
!!$  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: jv_ga         => NULL()

  ! PTR TO OZONE TRACER/CHANNEL OBJECT
  REAL(dp), DIMENSION(:,:,:), POINTER  :: ptr_O3   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER  :: ptr_O3te => NULL()

  ! cos(zenith angle)
  REAL(dp), POINTER, DIMENSION(:,:), PUBLIC :: cossza_2d => NULL()

  ! pressure field of input data (haloe or other)
  REAL(dp), DIMENSION(:,:,:), POINTER :: press_h => NULL()
  ! transmogrified input data (haloe ot other)
  REAL(dp), DIMENSION(:,:,:), POINTER :: o3_h => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: v3_h => NULL()

  REAL(dp), DIMENSION(:,:,:), POINTER :: o3_3d, v3_3d

  LOGICAL, DIMENSION(IP_MAX,2) :: lps = .FALSE.

  ! pointer to grid point channel objects :
  ! ---------------------------------------
  ! pointer to photolysis rate coeff.
  TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER  :: cloudj_gp => NULL()
  ! heating rates
  REAL(dp), DIMENSION(:,:,:), POINTER :: rh_o2_3d   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: rh_o3_3d   => NULL()

!!#D attila +
  ! pointer to Lagrangian channel objects (photolysis rates)
  TYPE(PTR_1D_ARRAY), PUBLIC, DIMENSION(:), POINTER :: cloudj_lg
!!#D attila -

  ! calcluate Lagrangian rate coefficients?
  LOGICAL :: l_calc_lg = .FALSE.

  INTEGER, PARAMETER :: GP=1, LG=2

  LOGICAL            :: lmtskip = .FALSE.
  REAL(DP), POINTER  :: rmode_mtskip => NULL()
#ifdef MESSYTENDENCY
  INTEGER            :: my_handle
#endif

  ! required for distinction between channel object and allocated space
  ! in free memory
  LOGICAL :: l_v3is_chaobj = .FALSE.

CONTAINS

  !***************************************************************************

  SUBROUTINE cloudj_initialize

    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_bcast, p_io
    USE messy_main_tools,        ONLY: find_next_free_unit

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! local
    CHARACTER(LEN=*), PARAMETER :: substr = 'cloudj_initialize'
    INTEGER :: io_unit ! logical I/O unit
    INTEGER :: status  ! status flag

    ! intitialize global switches/parameters
    IF (p_parallel_io) THEN
       io_unit = find_next_free_unit(10,99)
       CALL cloudj_read_nml_ctrl(status, io_unit)
       IF (status /= 0) CALL error_bi('error in cloudj_read_nml_ctrl',substr)
    ENDIF
    CALL p_bcast(NAMFIL_FJX_spec,     p_io)
    CALL p_bcast(NAMFIL_FJX_scat_cld, p_io)
    CALL p_bcast(NAMFIL_FJX_scat_aer, p_io)
    CALL p_bcast(NAMFIL_FJX_scat_UMa, p_io)
    CALL p_bcast(NAMFIL_atmos_std,    p_io)
    CALL p_bcast(NAMFIL_FJX_j2j,      p_io)
    ! intitialize global coupling switches/parameters
    IF (p_parallel_io) THEN
       io_unit = find_next_free_unit(100,200)
       CALL cloudj_read_nml_cpl(status, io_unit)
       IF (status /= 0) CALL error_bi('error in cloudj_read_nml_cpl',substr)
    ENDIF
    ! broadcast results
    CALL p_bcast(cloudj_O3%cha, p_io)
    CALL p_bcast(cloudj_O3%obj, p_io)
    CALL p_bcast(cloudj_o3h%cha, p_io)
    CALL p_bcast(cloudj_o3h%obj, p_io)
    CALL p_bcast(cloudj_v3h%cha, p_io)
    CALL p_bcast(cloudj_v3h%obj, p_io)
    CALL p_bcast(cloudj_pressh%cha, p_io)
    CALL p_bcast(cloudj_pressh%obj, p_io)
    CALL p_bcast(cloudj_cossza%cha, p_io)
    CALL p_bcast(cloudj_cossza%obj, p_io)
    CALL p_bcast(cloudj_alb%cha, p_io)
    CALL p_bcast(cloudj_alb%obj, p_io)
!    CALL p_bcast(cloudj_solar%cha, p_io)
!    CALL p_bcast(cloudj_solar%obj, p_io)
    CALL p_bcast(l_skip_lg, p_io)
    CALL p_bcast(l_force, p_io)
    CALL p_bcast(l_calc_lg, p_io)
!    CALL p_bcast(jv_aer_sca%cha, p_io)
!    CALL p_bcast(jv_aer_sca%obj, p_io)
!    CALL p_bcast(jv_aer_abs%cha, p_io)
!    CALL p_bcast(jv_aer_abs%obj, p_io)
!    CALL p_bcast(jv_aer_ga%cha, p_io)
!    CALL p_bcast(jv_aer_ga%obj, p_io)

  END SUBROUTINE cloudj_initialize

  !***************************************************************************
  SUBROUTINE cloudj_init_memory

    USE messy_main_grid_def_mem_bi, ONLY: ngpblks, nproma
    USE messy_main_tracer_mem_bi,   ONLY: ntrac_gp, ti_gp, &
                                          ntrac_lg, ti_lg, &
                                          t_trinfo_tp
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, LG_ATTILA       &
                                         , GP_3D_INT, gp_cnt, gp_memu &
                                         , gp_meml, gp_nseg, gp_start &
                                         , DC_GP, DIMID_LON, DIMID_LAT
    USE messy_main_channel_repr,         ONLY: new_representation        &
                                             , set_representation_decomp &
                                             , AUTO, IRANK, PIOTYPE_COL
    USE messy_main_channel_dimensions,   ONLY: new_dimension
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! local
    CHARACTER(LEN=STRLEN_MEDIUM) :: basename
    CHARACTER(LEN=*), PARAMETER :: substr = 'cloudj_init_memory'
    TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti => NULL()
    INTEGER :: j, jt, ntrac, status

    ! PARALLEL DECOMPOSITION
    INTEGER :: dimid_levh
    INTEGER :: GP_3D_LEVH
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ! NOTE: cloudj_gp memory must be allocated here, since pointers
    !       should be associated with channel objects
    ALLOCATE(cloudj_gp(IP_MAX))
    ALLOCATE(cloudj_2d(IP_MAX))
!!#D attila +
#ifdef ECHAM5
    IF (l_calc_lg) ALLOCATE(cloudj_lg(IP_MAX))
#endif
!!#D attila -
    DO jt=1, IP_MAX
      NULLIFY(cloudj_gp(jt)%PTR)
      NULLIFY(cloudj_2d(jt)%PTR)
!!#D attila +
#ifdef ECHAM5
      IF (l_calc_lg) NULLIFY(cloudj_lg(jt)%PTR)
#endif
!!#D attila -
    ENDDO

    ! J-values and UV-heating rates
    CALL new_channel(status, modstr//'_gp', reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

!!#D attila +
#ifdef ECHAM5
    IF (l_calc_lg) THEN
       CALL new_channel(status, modstr//'_lg', reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
    ENDIF
#endif
!!#D attila -

    ! diagnostic output
    CALL new_channel(status, modstr//'_diag', reprid=GP_3D_INT)
    CALL channel_halt(substr, status)

    if_lforce: IF (l_force) THEN
      lps(:,:) = .TRUE.
    ELSE
      gp_lg_loop: DO j = 1, 2
        SELECT CASE(j)
        CASE(GP)
          ntrac = ntrac_gp
          ti => ti_gp
        CASE(LG)
          IF (.NOT.l_calc_lg) CYCLE
          ntrac = ntrac_lg
          ti => ti_lg
        END SELECT
        tracer_loop: DO jt = 1, ntrac
          basename = ti(jt)%tp%ident%basename
          IF (TRIM(basename) == 'O2')       lps(ip_O2,j)       = .TRUE. !   1
          IF (TRIM(basename) == 'O3')   THEN
            lps(ip_O3p,j)  = .TRUE.                                     !   2
            lps(ip_O1d,j)  = .TRUE.                                     !   3
          ENDIF
          IF (TRIM(basename) == 'H2O2')     lps(ip_H2O2,j)     = .TRUE. !   4
          IF (TRIM(basename) == 'NO2')      lps(ip_NO2,j)      = .TRUE. !   5
          IF (TRIM(basename) == 'NO3')  THEN
            lps(ip_NO2O,j) = .TRUE.                                     !   6
            lps(ip_NOO2,j) = .TRUE.                                     !   7
          ENDIF
          IF (TRIM(basename) == 'N2O5')  THEN
             lps(ip_N2O5,j)   = .TRUE.                                  !   8
             lps(ip_NO3NOO,j) = .TRUE.                                  !  66
          ENDIF
          IF (TRIM(basename) == 'HNO3')     lps(ip_HNO3,j)     = .TRUE. !   9
          IF (TRIM(basename) == 'HNO4')     lps(ip_HNO4,j)     = .TRUE. !  10
          IF (TRIM(basename) == 'PAN')      lps(ip_PAN,j)      = .TRUE. !  11
          IF (TRIM(basename) == 'HONO')     lps(ip_HONO,j)     = .TRUE. !  12
          IF (TRIM(basename) == 'CH3OOH')   lps(ip_CH3OOH,j)   = .TRUE. !  13
          IF (TRIM(basename) == 'HCHO') THEN
            lps(ip_COH2,j) = .TRUE.                                     !  14
            lps(ip_CHOH,j) = .TRUE.                                     !  15
          ENDIF
          IF (TRIM(basename) == 'CH3CO3H')  lps(ip_CH3CO3H,j)  = .TRUE. !  16
          IF (TRIM(basename) == 'CH3CHO') THEN
            lps(ip_CH3CHO,j)      = .TRUE.                              !  17
            lps(ip_CH3CHO2VINY,j) = .TRUE.                              ! 122
          ENDIF
          IF (TRIM(basename) == 'CH3COCH3') lps(ip_CH3COCH3,j) = .TRUE. !  18
          IF (TRIM(basename) == 'MGLYOX')   lps(ip_MGLYOX,j)   = .TRUE. !  19
          IF (TRIM(basename) == 'HOCl')     lps(ip_HOCl,j)     = .TRUE. !  20
          IF (TRIM(basename) == 'OClO')     lps(ip_OClO,j)     = .TRUE. !  21
          IF (TRIM(basename) == 'Cl2O2')    lps(ip_Cl2O2,j)    = .TRUE. !  22
          IF (TRIM(basename) == 'ClNO3') THEN
            lps(ip_ClNO3,j)    = .TRUE.                                 !  23
            lps(ip_ClONO2,j)   = .TRUE.                                 !  67
          ENDIF
          IF (TRIM(basename) == 'ClNO2')    lps(ip_ClNO2,j)    = .TRUE. !  24
          IF (TRIM(basename) == 'Cl2')      lps(ip_Cl2,j)      = .TRUE. !  25
          IF (TRIM(basename) == 'BrO')      lps(ip_BrO,j)      = .TRUE. !  26
          IF (TRIM(basename) == 'HOBr')     lps(ip_HOBr,j)     = .TRUE. !  27
          IF (TRIM(basename) == 'BrCl')     lps(ip_BrCl,j)     = .TRUE. !  28
          IF (TRIM(basename) == 'BrNO3')    lps(ip_BrNO3,j)    = .TRUE. !  29
          IF (TRIM(basename) == 'BrNO2')    lps(ip_BrNO2,j)    = .TRUE. !  30
          IF (TRIM(basename) == 'Br2')      lps(ip_Br2,j)      = .TRUE. !  31
          IF (TRIM(basename) == 'CCl4')     lps(ip_CCl4,j)     = .TRUE. !  32
          IF (TRIM(basename) == 'CH3Cl')    lps(ip_CH3Cl,j)    = .TRUE. !  33
          IF (TRIM(basename) == 'CH3CCl3')  lps(ip_CH3CCl3,j)  = .TRUE. !  34
          IF (TRIM(basename) == 'CFCl3')    lps(ip_CFCl3,j)    = .TRUE. !  35
          IF (TRIM(basename) == 'CF2Cl2')   lps(ip_CF2Cl2,j)   = .TRUE. !  36
          IF (TRIM(basename) == 'CH3Br')    lps(ip_CH3Br,j)    = .TRUE. !  37
          IF (TRIM(basename) == 'CF2ClBr')  lps(ip_CF2ClBr,j)  = .TRUE. !  38
          IF (TRIM(basename) == 'CF3Br')    lps(ip_CF3Br,j)    = .TRUE. !  39
          IF (TRIM(basename) == 'CH3I')     lps(ip_CH3I,j)     = .TRUE. !  40
          IF (TRIM(basename) == 'C3H7I')    lps(ip_C3H7I,j)    = .TRUE. !  41
          IF (TRIM(basename) == 'CH2ClI')   lps(ip_CH2ClI,j)   = .TRUE. !  42
          IF (TRIM(basename) == 'CH2I2')    lps(ip_CH2I2,j)    = .TRUE. !  43
          IF (TRIM(basename) == 'IO')       lps(ip_IO,j)       = .TRUE. !  44
          IF (TRIM(basename) == 'HOI')      lps(ip_HOI,j)      = .TRUE. !  45
          IF (TRIM(basename) == 'I2')       lps(ip_I2,j)       = .TRUE. !  46
          IF (TRIM(basename) == 'ICl')      lps(ip_ICl,j)      = .TRUE. !  47
          IF (TRIM(basename) == 'IBr')      lps(ip_IBr,j)      = .TRUE. !  48
          IF (TRIM(basename) == 'INO2')     lps(ip_INO2,j)     = .TRUE. !  49
          IF (TRIM(basename) == 'INO3')     lps(ip_INO3,j)     = .TRUE. !  50
          IF (TRIM(basename) == 'SO2')      lps(ip_SO2,j)      = .TRUE. !  51
          IF (TRIM(basename) == 'SO3')      lps(ip_SO3,j)      = .TRUE. !  52
          IF (TRIM(basename) == 'OCS')      lps(ip_OCS,j)      = .TRUE. !  53
          IF (TRIM(basename) == 'CS2')      lps(ip_CS2,j)      = .TRUE. !  54
          ! for definition of lps(ip_H2O,j) see below
          !IF(TRIM(basename) == 'H2O')      lps(ip_H2O,j)      = .TRUE. !  55
          IF (TRIM(basename) == 'N2O')      lps(ip_N2O,j)      = .TRUE. !  56
          IF (TRIM(basename) == 'NO')       lps(ip_NO,j)       = .TRUE. !  57
          IF (TRIM(basename) == 'CO2')      lps(ip_CO2,j)      = .TRUE. !  58
          IF (TRIM(basename) == 'HCl')      lps(ip_HCl,j)      = .TRUE. !  59
          IF (TRIM(basename) == 'CHCl2Br')  lps(ip_CHCl2Br,j)  = .TRUE. !  60
          IF (TRIM(basename) == 'CHClBr2')  lps(ip_CHClBr2,j)  = .TRUE. !  61
          IF (TRIM(basename) == 'CH2ClBr')  lps(ip_CH2ClBr,j)  = .TRUE. !  62
          IF (TRIM(basename) == 'CH2Br2')   lps(ip_CH2Br2,j)   = .TRUE. !  63
          IF (TRIM(basename) == 'CHBr3')    lps(ip_CHBr3,j)    = .TRUE. !  64
          IF (TRIM(basename) == 'SF6')      lps(ip_SF6,j)      = .TRUE. !  65
          IF (TRIM(basename) == 'MACR')     lps(ip_MACR,j)     = .TRUE. !  68
          IF (TRIM(basename) == 'MVK')      lps(ip_MVK,j)      = .TRUE. !  69
          IF (TRIM(basename) == 'GLYOX')    lps(ip_GLYOX,j)    = .TRUE. !  70
          IF (TRIM(basename) == 'HOCH2CHO') lps(ip_HOCH2CHO,j) = .TRUE. !  71
          IF (TRIM(basename) == 'CH4')      lps(ip_CH4,j)      = .TRUE. !  72
          IF (TRIM(basename) == 'H2SO4')    lps(ip_H2SO4,j)    = .TRUE. ! 102
          IF (TRIM(basename) == 'CH3NO3')   lps(ip_CH3NO3,j)   = .TRUE. ! 104
          IF (TRIM(basename) == 'CH3O2NO2') lps(ip_CH3O2NO2,j) = .TRUE. ! 105
          IF (TRIM(basename) == 'CH3ONO')   lps(ip_CH3ONO,j)   = .TRUE. ! 106
          IF (TRIM(basename) == 'CH3O2')    lps(ip_CH3O2,j)    = .TRUE. ! 107
          IF (TRIM(basename) == 'HCOOH')    lps(ip_HCOOH,j)    = .TRUE. ! 108
          !IF (TRIM(basename) == 'HO2NO2')   lps(ip_HO2NO2,j)   = .TRUE. ! 109
          !IF (TRIM(basename) == 'OHNO3')    lps(ip_OHNO3,j)    = .TRUE. ! 110
          !qqqdummy                                                      ! 111
          !IF (TRIM(basename) == 'CH3OCl')   lps(ip_CH3OCl,j)   = .TRUE. ! 112
          !IF (TRIM(basename) == 'MEO2NO2')  lps(ip_MEO2NO2,j)  = .TRUE. ! 113
          IF (TRIM(basename) == 'CHF2Cl')   lps(ip_CHF2Cl,j)   = .TRUE. ! 114
          !IF (TRIM(basename) == 'F113')     lps(ip_F113,j)     = .TRUE. ! 115
          IF (TRIM(basename) == 'C2H5NO3')  lps(ip_C2H5NO3,j)   = .TRUE. ! 116
          IF (TRIM(basename) == 'NOA')      lps(ip_NOA,j)       = .TRUE. ! 117
          IF (TRIM(basename) == 'LMEKNO3')  lps(ip_MEKNO3,j)    = .TRUE. ! 118
          IF (TRIM(basename) == 'BENZAL')   lps(ip_BENZAL,j)    = .TRUE. ! 119
          IF (TRIM(basename) == 'TOL1OHNO2') lps(ip_HOPh3Me2NO2,j) = .TRUE. ! 120
          IF (TRIM(basename) == 'HOC6H4NO2') lps(ip_HOC6H4NO2,j) = .TRUE. ! 121
          IF (TRIM(basename) == 'CH3COCO2H') lps(ip_CH3COCO2H,j) = .TRUE. ! 123
          IF (TRIM(basename) == 'IPRCHO')    lps(ip_IPRCHO2HCO,j) = .TRUE. ! 124
          IF (TRIM(basename) == 'C2H5CHO') THEN
             lps(ip_C2H5CHO2HCO,j)   = .TRUE. ! 125
             lps(ip_C2H5CHO2ENOL,j)  = .TRUE. ! 126
          ENDIF
          IF (TRIM(basename) == 'C3H7CHO') THEN
            lps(ip_C3H7CHO2HCO,j)   = .TRUE.  ! 127
            lps(ip_C3H7CHO2VINY,j)  = .TRUE.  ! 128
          ENDIF
          IF (TRIM(basename) == 'HVMK')      lps(ip_PeDIONE24,j) = .TRUE. ! 129
          IF (TRIM(basename) == 'PINAL') THEN
             lps(ip_PINAL2HCO,j)    = .TRUE.  ! 130
            lps(ip_PINAL2ENOL,j)    = .TRUE.  ! 131
         ENDIF
         IF (TRIM(basename) == 'CF2ClCFCl2')  lps(ip_CF2ClCFCl2,j) = .TRUE. ! 132
         IF (TRIM(basename) == 'CH3CFCl2')    lps(ip_CH3CFCl2,j)   = .TRUE. ! 133
         IF (TRIM(basename) == 'CF3CF2Cl')    lps(ip_CF3CF2Cl,j)   = .TRUE. ! 134
         IF (TRIM(basename) == 'CF2ClCF2Cl')  lps(ip_CF2ClCF2Cl,j) = .TRUE. ! 135
         IF (TRIM(basename) == 'CHCl3')       lps(ip_CHCl3,j)      = .TRUE. ! 136
         IF (TRIM(basename) == 'CH2Cl2')      lps(ip_CH2Cl2,j)     = .TRUE. ! 137

        ENDDO tracer_loop

        ! H2O photolysis is always required due to its presence in MECCA
        ! (q is used if there is no H2O tracer)
        lps(ip_H2O,j)      = .TRUE. !  55

      ENDDO gp_lg_loop

    ENDIF if_lforce

    ! for core file:
    lp(:) = lps(:,GP).OR.lps(:,LG)

    ! make grid point channel objects
    DO jt=1, IP_MAX
      IF (lps(jt,GP).OR.lps(jt,LG)) THEN
         CALL new_channel_object(status, modstr//'_gp'&
              , 'J_'//TRIM(jname(jt)) &
              , p3=cloudj_gp(jt)%PTR  )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr//'_gp' &
              , 'J_'//TRIM(jname(jt))             &
              , 'long_name', c='J('//TRIM(jname(jt))//')' )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr//'_gp' &
              , 'J_'//TRIM(jname(jt))             &
              , 'units', c='1/s')
         CALL channel_halt(substr, status)
         CALL info_bi('channel/object '//modstr//'_gp/J_'// &
           TRIM(jname(jt))//' was created')
      ENDIF
    ENDDO

!!#D attila +
#ifdef ECHAM5
    IF (l_calc_lg) THEN
      ! make Lagrangian channel objects
      DO jt=1, IP_MAX
        IF (lps(jt,LG)) THEN
         CALL new_channel_object(status, modstr//'_lg'&
              , 'J_'//TRIM(jname(jt)) &
              , p1=cloudj_lg(jt)%PTR  )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr//'_lg' &
              , 'J_'//TRIM(jname(jt))             &
              , 'long_name', c='J('//TRIM(jname(jt))//')' )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr//'_lg' &
              , 'J_'//TRIM(jname(jt))             &
              , 'units', c='1/s')
         CALL channel_halt(substr, status)
         CALL info_bi('channel/object '//modstr//'_lg/J_'// &
           TRIM(jname(jt))//' was created')
        ENDIF
      ENDDO
    ENDIF
#endif
!!#D attila -
    ! ----------------------------------------------------------

    ! make channel objects for diagnostic output
    CALL new_channel_object(status, modstr//'_diag', 'v3', p3=v3_3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_diag', 'v3' &
         , 'long_name', c='ozone column')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_diag', 'v3', 'units', c='mcl/cm2')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object '//modstr//'_diag/v3 was created')

    CALL new_channel_object(status, modstr//'_diag', 'O3', p3=o3_3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_diag', 'O3' &
         , 'long_name', c='ozone')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_diag', 'O3', 'units', c='mol/mol')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object '//modstr//'_diag/O3 was created')

   ! make channel objects for diagnostic output
#ifdef OZONECLIM
    CALL new_dimension(status, DIMID_LEVH, 'leveldim36', 36)
    CALL channel_halt(substr, status)

     ! -----------------------------------------------------

    nseg = gp_nseg

    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    CALL new_representation(status, GP_3D_LEVH, 'GP_3D_LEVH' &
         , rank = 3, link = 'xxx-', dctype = DC_GP                &
         , dimension_ids = (/ DIMID_LON, DIMID_LEVH, DIMID_LAT /) &
         , ldimlen       = (/ nproma   , AUTO      , ngpblks   /) &
         , output_order  = (/ 1,3,2 /)                            &
         , axis = 'XZY-'                                          &
         )
    CALL channel_halt(substr, status)

    start(:,:)= gp_start(:,:)
    meml(:,:) = gp_meml(:,:)

    cnt(:,:)  = gp_cnt(:,:)
    memu(:,:) = gp_memu(:,:)

    start(:,4) = 1
    cnt(:,4) = 1
    meml(:,4) = 1
    memu(:,4) = 1

    start(:,3) = gp_start(:,4)
    cnt(:,3)   = gp_cnt(:,4)
    meml(:,3)  = gp_meml(:,4)
    memu(:,3)  = gp_memu(:,4)


    cnt(:,2)  = 36
    memu(:,2) = 36

    CALL set_representation_decomp(status, GP_3D_LEVH &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    !----- DEALLOCATE------------------------------------------
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)
    ! ---------------------------------------------------------
    ! -----------------------------------------------------

    ! lrestreq = .TRUE. required for MECO(n)
    CALL new_channel_object(status, modstr//'_diag', 'v3h', p3=v3_h &
         , reprid = GP_3D_LEVH, lrestreq = .TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_diag', 'v3h' &
         , 'long_name', c='ozone column')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_diag', 'v3h', 'units', c='mcl/cm2')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object '//modstr//'_diag/v3h was created')

    ! lrestreq = .TRUE. required for MECO(n)
    CALL new_channel_object(status, modstr//'_diag', 'pressih', p3=press_h &
         , reprid = GP_3D_LEVH, lrestreq = .TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_diag', 'pressih' &
         , 'long_name', c='ozone')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_diag', 'pressih', 'units', c='mol/mol')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object '//modstr//'_diag/O3h was created')
    l_v3is_chaobj = .TRUE.
#endif
    CALL end_message_bi(modstr,   'CHANNEL DEFINITION', substr)

    CALL init_cloudj

  END SUBROUTINE cloudj_init_memory

  !***************************************************************************

  SUBROUTINE cloudj_init_coupling

    USE messy_main_mpi_bi,             ONLY: p_parallel_io
    USE messy_main_tracer_mem_bi,      ONLY: GPTRSTR, gp_channel
    USE messy_main_channel_error_bi,   ONLY: channel_halt
    USE messy_main_channel_bi,         ONLY: GP_3D_MID, DC_GP  &
                                           , DIMID_LON, DIMID_LAT, DIMID_LEV &
                                           , DC_BC &
                                           , gp_nseg, gp_start, gp_cnt &
                                           , gp_meml, gp_memu
    USE messy_main_channel_tracer,     ONLY: set_channel_or_tracer
    USE messy_main_channel,            ONLY: get_channel_object &
                                           , get_channel_object_dimvar &
                                           , get_channel_info &
                                           , new_channel_object, new_attribute &
                                           , new_channel &
                                           , get_channel_object_dimvalue
    USE messy_main_channel_dimensions, ONLY: new_dimension
    USE messy_main_channel_repr,       ONLY: new_representation, AUTO &
                                           , set_representation_decomp &
                                           , IRANK, PIOTYPE_COL
    USE messy_main_grid_def_mem_bi,    ONLY: nproma, nlev, ngpblks
    USE messy_main_tools,              ONLY: PTR_1D_ARRAY
    USE messy_main_constants_mem,      ONLY: STRLEN_ULONG, iouerr

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER   :: substr = 'cloudj_init_coupling'
    INTEGER                       :: status    ! error flag
    TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER          :: dvs => NULL()
    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: units => NULL()
    INTEGER :: i
    REAL(DP), POINTER                 :: rcloudj  => NULL()

    REAL(dp), DIMENSION(:,:,:,:),   POINTER ::  mem => NULL()
    INTEGER                               :: DIMID_JV
    INTEGER                               :: REPR_JV_AER_4D
    CHARACTER(LEN=1)                      :: ichar
    CHARACTER(LEN=2)                      :: str_num
    CHARACTER(LEN=32)                     :: cname2
    INTEGER                               :: ji

    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    REAL(dp), DIMENSION(:), POINTER :: plevs => NULL()
    INTEGER                         :: jk
    INTEGER                         :: nlev_h

    CALL start_message_bi(modstr, 'COUPLING INITIALIZING', substr)

    !Logic for MTSKIP is set here
    lmtskip = .FALSE.
    CALL get_channel_object(status,'mtskip', 'rcloudj',p0=rcloudj)
    lmtskip = (status == 0)
    ! Note: not with AND in same line, because rcloudj might not be associated
    IF (lmtskip) THEN
       lmtskip = (NINT(rcloudj)==1)
    END IF
    IF (lmtskip) THEN
       ! actual mode of mtskip (0,1,2)
       CALL get_channel_object(status,'mtskip', 'rmode', p0=rmode_mtskip)
       CALL channel_halt(substr, status)
    ENDIF

    ! CHECK IF OZONE IS AVAILABLE
    CALL info_bi('Looking for OZONE ... ')
    CALL info_bi('       channel: '//cloudj_o3%cha)
    CALL info_bi('       object : '//cloudj_o3%obj)

    CALL set_channel_or_tracer(status, GPTRSTR, gp_channel &
         , cloudj_o3%cha, cloudj_o3%obj, ptr_O3, ptr_O3te)
    CALL channel_halt(substr, status)

    ! CHECK IF OZONE IS AVAILABLE
    CALL info_bi('Looking for INPUT DATA UBC OZONE ...')
    CALL info_bi('... O3_h ...')
    CALL info_bi('       channel: '//cloudj_o3h%cha)
    CALL info_bi('       object : '//cloudj_o3h%obj)

    CALL get_channel_object(status  &
         ,  TRIM(cloudj_o3h%cha), TRIM(cloudj_o3h%obj), p3=O3_h)
    CALL channel_halt(substr, status)

    ! get pressure levels via dimension variable
    CALL get_channel_object_dimvalue(status       &
         , TRIM(cloudj_o3h%cha), TRIM(cloudj_o3h%obj) &
         , data= plevs, axis ='Z')

    IF (status /=  953) THEN  !E 953 : DIMVAR does not exist
       CALL channel_halt(substr, status)
       ! analyse attribute / define pressure levels
       nlev_h = SIZE(O3_h, _IZ_XYZ__ )
       IF (.NOT.ASSOCIATED(press_h)) &
            ALLOCATE(press_h(_RI_XYZ__(nproma,ngpblks,nlev_h)))

       IF (SIZE(plevs) /= nlev_h) THEN
          write(iouerr,*) 'CLOUDJ: ', nlev_h, SIZE(plevs)
          CALL error_bi(&
            'number of vertical pressure levels inconsistent', substr)
       END IF
       DO jk = 1, nlev_h
          press_h(_RI_XYZ__(:,:,jk)) = plevs(jk)
       END DO
       IF (ASSOCIATED(plevs)) THEN
          DEALLOCATE(plevs) ; NULLIFY(plevs)
       END IF
    ELSE
       CALL info_bi('... V3_h ...')
       CALL info_bi('       channel: '//cloudj_v3h%cha)
       CALL info_bi('       object : '//cloudj_v3h%obj)

       CALL get_channel_object(status  &
            ,  TRIM(cloudj_v3h%cha), TRIM(cloudj_v3h%obj), p3=v3_h)
       CALL channel_halt(substr, status)

       CALL info_bi('... press_h ...')
       CALL info_bi('       channel: '//cloudj_pressh%cha)
       CALL info_bi('       object : '//cloudj_pressh%obj)

       CALL get_channel_object(status  &
            ,  TRIM(cloudj_pressh%cha), TRIM(cloudj_pressh%obj), p3=press_h)
       CALL channel_halt(substr, status)
       l_v3is_chaobj = .TRUE.
    END IF

    CALL info_bi('Looking for cos(solar zenith angle) (cossza)')
    CALL info_bi('       channel: '//cloudj_cossza%cha)
    CALL info_bi('       object : '//cloudj_cossza%obj)

    CALL get_channel_object(status  &
         ,  TRIM(cloudj_cossza%cha), TRIM(cloudj_cossza%obj), p2=cossza_2d)
    CALL channel_halt(substr, status)

    CALL info_bi('       channel: '//cloudj_alb%cha)
    CALL info_bi('       object : '//cloudj_alb%obj)
    CALL get_channel_object(status  &
         ,  TRIM(cloudj_alb%cha), TRIM(cloudj_alb%obj), p2=albedo)
    CALL channel_halt(substr, status)

  END SUBROUTINE cloudj_init_coupling

  !***************************************************************************

  ! new routine for HALOE input, V3_H column above lev jk
  SUBROUTINE cloudj_global_start

    USE messy_main_grid_def_mem_bi, ONLY: nproma, npromz, ngpblks
    USE messy_main_data_bi,         ONLY: press_3d
    USE messy_main_timer,           ONLY: lstart
    USE messy_main_tools,           ONLY: combine_o3_fields

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE, TRIM

    ![ part./cm^2 * 1/pa]
    REAL(dp), PARAMETER :: sp = N_A * 1000./(M_air*g) * 1.e-4
    CHARACTER(LEN=*), PARAMETER :: substr = 'cloudj_global_start'
    INTEGER                     :: jp,jk,zkproma,zjrow

    ! vertical interpolation of is performed online
    ! (O3, temperature above highest model level)
    INTEGER :: nlev_h ! number of vertical levels in input data

    INTEGER :: jrow, kproma
    INTEGER :: status

    nlev_h = SIZE(o3_h,_IZ_XYZ__)
    IF (.NOT.ASSOCIATED(v3_h)) THEN
       ALLOCATE(v3_h(_RI_XYZ__(nproma,ngpblks,nlev_h)))
       v3_h = 0._dp
    END IF

    DO  zjrow =1,ngpblks
#ifndef CESM1
       IF ( zjrow   == ngpblks ) THEN
          zkproma = npromz
       ELSE
          zkproma = nproma
       ENDIF
#else
       zkproma = npromz(zjrow)
#endif

       DO jp = 1,zkproma
          V3_H(_RI_XYZ__(jp,zjrow,1)) = &
               sp * press_h(_RI_XYZ__(jp,zjrow,1)) * o3_h(_RI_XYZ__(jp,zjrow,1))
       ENDDO
       DO jk = 2,nlev_h
          DO jp = 1,zkproma
             V3_H(_RI_XYZ__(jp,zjrow,jk)) = V3_H(_RI_XYZ__(jp,zjrow,jk-1)) + &
                  sp * (press_h(_RI_XYZ__(jp,zjrow,jk)) &
                  - press_h(_RI_XYZ__(jp,zjrow,jk-1))) * &
                  0.5 * (O3_h(_RI_XYZ__(jp,zjrow,jk)) + &
                  O3_h(_RI_XYZ__(jp,zjrow,jk-1)))
          ENDDO
       ENDDO
    ENDDO

    IF (lstart) THEN
       !-----------------------------------------------------------------------
       ! combine O3 fields and calculate column density of O3
       ! ozone tracer field
       DO jrow =1,ngpblks
#ifndef CESM1
          IF ( jrow   == ngpblks ) THEN
             kproma = npromz
          ELSE
             kproma = nproma
          ENDIF
#else
          kproma = npromz(jrow)
#endif
          CALL combine_o3_fields( &
               ptr_O3(_RI_XYZ__(1:kproma,jrow,:)),    &
               press_3d(_RI_XYZ__(1:kproma,jrow,:)),  &
               press_h(_RI_XYZ__(1:kproma,jrow,:)),   &
               o3_h(_RI_XYZ__(1:kproma,jrow,:)),      &
               v3_h(_RI_XYZ__(1:kproma,jrow,:)),      &
               o3_3d(_RI_XYZ__(1:kproma,jrow,:)),     &
               v3_3d(_RI_XYZ__(1:kproma,jrow,:)))
       ENDDO
    ENDIF
    !------------------------------------------------------------------------

  END SUBROUTINE cloudj_global_start

  !***************************************************************************

  SUBROUTINE cloudj_physc

    USE messy_main_timer,           ONLY: time_step_len
    USE messy_main_grid_def_mem_bi, ONLY: nlev, ngpblks                    &
                                       , jrow, kproma, nproma
    USE messy_main_data_bi,    ONLY: press_3d, pressi_3d                   &
                                   , tm1_3d, tte_3d, rhum_3d               &
                                   , xlm1_3d, xlte_3d, xim1_3d, xite_3d    &
                                   , aclc, slf
    USE messy_main_tools,      ONLY: combine_o3_fields
#if !(defined(ECHAM5) || defined(CESM1) || defined (COSMO))
    ! this is a temporary workaround
    USE messy_main_data_bi,    ONLY: albedo
#endif

    IMPLICIT NONE
    INTRINSIC :: MAX, ASSOCIATED

    CHARACTER(LEN=*), PARAMETER :: substr = 'cloudj_physc'
    INTEGER :: jk, jp, jt

    REAL(dp) :: temp_2d(nproma,nlev)
    REAL(dp) :: xl_2d(nproma,nlev)
    REAL(dp) :: xi_2d(nproma,nlev)
    REAL(dp) :: clp_2d(nproma,nlev)
    REAL(dp) :: zptr_O3(_RI_XYZ__(nproma,ngpblks,nlev))
    INTEGER  :: isw_mtskip

    IF (lmtskip) THEN
       isw_mtskip = NINT(rmode_mtskip)
    ELSE
       isw_mtskip = 0
    ENDIF

    skip_calc: IF (isw_mtskip <= 1 ) THEN

       temp_2d(:,:) = tm1_3d(_RI_XYZ__(:,jrow,:))  &
            + tte_3d(_RI_XYZ__(:,jrow,:))  * time_step_len
       xl_2d(:,:)   = xlm1_3d(_RI_XYZ__(:,jrow,:)) &
            + xlte_3d(_RI_XYZ__(:,jrow,:)) * time_step_len
       xi_2d(:,:)   = xim1_3d(_RI_XYZ__(:,jrow,:)) &
            + xite_3d(_RI_XYZ__(:,jrow,:)) * time_step_len

       DO jk=1,nlev
          DO jp=1,kproma
             ! for cloudiness < 1% assume clear sky
             IF (aclc(_RI_XYZ__(jp,jrow,jk)) < 0.01_dp) THEN
                clp_2d(jp,jk) = 0._dp
             ELSE
                clp_2d(jp,jk) = &
                     MAX( 0._dp , (xl_2d(jp,jk)+xi_2d(jp,jk)) &
                     / aclc(_RI_XYZ__(jp,jrow,jk)) ) &
                     * (1000._dp/g) * &
                     (pressi_3d(_RI_XYZ__(jp,jrow,jk+1)) &
                     - pressi_3d(_RI_XYZ__(jp,jrow,jk)))
             ENDIF
             ! cloud ice path added, will be treated separately later
          ENDDO
       ENDDO

       ! For some reason the intel compiler (16, 17) optimises away the
       ! MAX(0._dp, ...) statement in the above loop structure ...?
       ! This is a workaround:
       clp_2d(1:kproma,:) = MAX(0.0_dp, clp_2d(1:kproma,:))

       !------------------------------------------------------------------------
       IF (ASSOCIATED(ptr_O3te)) THEN
          zptr_O3(_RI_XYZ__(:,jrow,:)) = ptr_O3(_RI_XYZ__(:,jrow,:)) &
               + time_step_len*ptr_O3te(_RI_XYZ__(:,jrow,:))
       ELSE
          zptr_O3(_RI_XYZ__(:,jrow,:)) = ptr_O3(_RI_XYZ__(:,jrow,:))
       ENDIF
       !------------------------------------------------------------------------

       ! combine O3 fields and calculate column density of O3
       ! ozone tracer field
       CALL combine_o3_fields( &
            zptr_O3(_RI_XYZ__(1:kproma,jrow,:)),   &
            press_3d(_RI_XYZ__(1:kproma,jrow,:)),  &
            press_h(_RI_XYZ__(1:kproma,jrow,:)),   &
            o3_h(_RI_XYZ__(1:kproma,jrow,:)),      &
            v3_h(_RI_XYZ__(1:kproma,jrow,:)),      &
            o3_3d(_RI_XYZ__(1:kproma,jrow,:)),     &
            v3_3d(_RI_XYZ__(1:kproma,jrow,:)))
       !------------------------------------------------------------------------

       ! set 2d-pointer for use in core
       DO jt = 1, IP_MAX
          ! test only gp
          IF (lp(jt)) &
               cloudj_2d(jt)%ptr => cloudj_gp(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:))
       ENDDO
       !------------------------------------------------------------------------

       ! calculate J-values and UV-heating rates
       CALL cloudjvalues(                                  &
            REAL(v3_3d(_RI_XYZ__(1:kproma,jrow,:)),SP),    &
            REAL(cossza_2d(1:kproma,jrow),SP),             &
            REAL(press_3d(_RI_XYZ__(1:kproma,jrow,:)),SP), &
            REAL(o3_3d(_RI_XYZ__(1:kproma,jrow,:)),SP),    &
            REAL(rhum_3d(_RI_XYZ__(1:kproma,jrow,:)),SP),  &
            REAL(temp_2d(1:kproma,:),SP),                  &
            REAL(albedo(1:kproma,jrow),SP),                &
            REAL(aclc(_RI_XYZ__(1:kproma,jrow,:)),SP),     &
            REAL(slf(1:kproma,jrow),SP),                   &
            REAL(clp_2d(1:kproma,:),SP))

    ENDIF skip_calc

  END SUBROUTINE cloudj_physc

  !***************************************************************************

  SUBROUTINE cloudj_global_end

#ifdef ECHAM5
!!#D attila +
    USE messy_attila_tools_e5, ONLY: gp2lg_e5
!!#D attila -

    IMPLICIT NONE
!!#D attila +
    INTEGER :: jt

    ! Lagrangian calculation
    IF (l_calc_lg) THEN
      DO jt=1, IP_MAX
        IF (lps(jt,LG)) THEN
           CALL gp2lg_e5(cloudj_gp(jt)%PTR, cloudj_lg(jt)%PTR, lmcons=.FALSE.)
        ENDIF
      ENDDO
    ENDIF
!!#D attila -
#endif
  END SUBROUTINE cloudj_global_end

  !***************************************************************************

  SUBROUTINE cloudj_free_memory

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

#ifndef OZONECLIM
    IF (.NOT. l_v3is_chaobj) THEN
       IF (ASSOCIATED(v3_h)) THEN
          DEALLOCATE(v3_h) ; NULLIFY(v3_h)
       END IF
       IF (ASSOCIATED(press_h)) THEN
          DEALLOCATE(press_h) ; NULLIFY(press_h)
       END IF
    END IF
#endif
    IF (ASSOCIATED(cloudj_gp)) THEN
       DEALLOCATE(cloudj_gp) ; NULLIFY(cloudj_gp)
    END IF
    IF (ASSOCIATED(cloudj_2d)) THEN
       DEALLOCATE(cloudj_2d); NULLIFY(cloudj_2d)
    END IF
!!#D attila +
    IF (ASSOCIATED(cloudj_lg)) THEN
       DEALLOCATE(cloudj_lg); NULLIFY(cloudj_lg)
    END IF
!!#D attila -

  END SUBROUTINE cloudj_free_memory

  !***************************************************************************

  SUBROUTINE cloudj_read_nml_cpl(status, iou)

    ! photolysis module routine
    ! read photolysis namelist, check it, and initialize global variables
    ! Author: Patrick Joeckel, MPICH, Nov 2002

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_main_tracer_mem_bi, ONLY: NGCELL

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status  ! error status
    INTEGER, INTENT(IN)  :: iou     ! logical I/O unit

    NAMELIST /CPL/ cloudj_O3, l_skip_lg, l_force, cloudj_cossza &
          , cloudj_o3h, cloudj_v3h, cloudj_pressh, cloudj_alb

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr = 'cloudj_read_nml_cpl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! init
    status = 1 ! default: error

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    IF (TRIM(cloudj_o3%cha) == '') THEN
       CALL info_bi('ERROR: empty channel name for OZONE')
       RETURN
    ELSE
       CALL info_bi('OZONE channel :'//cloudj_o3%cha)
    END IF
    IF (TRIM(cloudj_o3%obj) == '') THEN
       CALL info_bi('ERROR: empty channel object name for OZONE')
       RETURN
    ELSE
       CALL info_bi('OZONE object  :'//cloudj_o3%obj)
    END IF

   IF (TRIM(cloudj_alb%cha) == '') THEN
       CALL info_bi('Albedo from radiation channel will be used')
#if defined(ECHAM5) || defined(CESM1)
       cloudj_alb%cha = 'rad'
#endif
#if defined(COSMO)
       cloudj_alb%cha = 'COSMO_ORI'
#endif
    ELSE
       CALL info_bi('Albedo channel :'//cloudj_alb%cha)
    END IF
    IF (TRIM(cloudj_alb%obj) == '') THEN
       CALL info_bi('Albedo from radiation channel will be used')
#if defined(ECHAM5) || defined(CESM1)
       cloudj_alb%obj = 'albedo'
#endif
#if defined(COSMO)
       cloudj_alb%obj = 'ALB_RAD'
#endif
    ELSE
       CALL info_bi('Albedo object  :'//cloudj_alb%obj)
    END IF

   IF (TRIM(cloudj_o3h%cha) == '') THEN
       CALL info_bi('ERROR: empty channel name for O3_h')
       RETURN
    ELSE
       CALL info_bi('O3_h channel :'//cloudj_o3h%cha)
    END IF
    IF (TRIM(cloudj_o3h%obj) == '') THEN
       CALL info_bi('ERROR: empty channel object name for O3_h')
       RETURN
    ELSE
       CALL info_bi('O3_h object  :'//cloudj_o3h%obj)
    END IF

#if defined(COSMO) || defined(VERTICO)
    IF (TRIM(cloudj_v3h%cha) == '') THEN
       CALL info_bi('INFO: empty channel name for v3_h')
    ELSE
       CALL info_bi('v3_h channel :'//cloudj_v3h%cha)
    END IF
    IF (TRIM(cloudj_v3h%obj) == '') THEN
       CALL info_bi('INFO: empty channel object name for v3_h')
    ELSE
       CALL info_bi('v3_h object  :'//cloudj_v3h%obj)
    END IF

    IF (TRIM(cloudj_pressh%cha) == '') THEN
       CALL info_bi('INFO: empty channel name for press_h')
    ELSE
       CALL info_bi('press_h channel :'//cloudj_pressh%cha)
    END IF
    IF (TRIM(cloudj_pressh%obj) == '') THEN
       CALL info_bi('INFO: empty channel object name for press_h')
    ELSE
       CALL info_bi('press_h object  :'//cloudj_pressh%obj)
    END IF
#endif

    IF (TRIM(cloudj_cossza%cha) == '') THEN
       CALL info_bi('ERROR: empty channel name for cossza')
       RETURN
    ELSE
       CALL info_bi('cossza channel :'//cloudj_cossza%cha)
    END IF
    IF (TRIM(cloudj_cossza%obj) == '') THEN
       CALL info_bi('ERROR: empty channel object name for cossza')
       RETURN
    ELSE
       CALL info_bi('cossza object  :'//cloudj_cossza%obj)
    END IF

    IF ((.NOT. l_skip_lg) .AND. (NGCELL > 0)) THEN
       l_calc_lg = .TRUE.
!!#D attila +
       CALL info_bi('Lagrangian: ON')
!!#D attila -
    ELSE
       IF (.NOT. l_skip_lg) THEN
!!#D attila +
         CALL info_bi('l_skip_lg = F in namelist')
         CALL info_bi('However no Lagrangian scheme activated ...')
         CALL info_bi(' ... setting l_calc_lg = F')
!!#D attila -
       ENDIF
       l_calc_lg = .FALSE.
!!#D attila +
       CALL info_bi('Lagrangian: OFF')
!!#D attila -
    ENDIF

    IF (l_force) THEN
       CALL info_bi('Force calculation and output of all j-values: ON')
    ELSE
       CALL info_bi('Force calculation and output of all j-values: OFF')
    ENDIF

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! no error

  END SUBROUTINE cloudj_read_nml_cpl

!*****************************************************************************

END MODULE messy_cloudj_si

!*****************************************************************************
