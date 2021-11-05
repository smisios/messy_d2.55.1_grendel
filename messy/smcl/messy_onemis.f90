!*****************************************************************************
!
!*****************************************************************************

! online emissions (onemis)
! Authors:
!   miscellaneous (see individual subroutines)
!   Rolf Sander, 2004: rewritten as MESSy submodel
!   Astrid Kerkweg, 2004: expanded by all than the sea salt and DMS emissions
!   Joerg  Steinkamp 2009: expanded by Yienger and Levy 2 NOx emission scenario
!   Susannah Burrows, 2009: expanded by bioaerosol emissions (Olson and MODIS)
!   Gregor Glaeser, 2010: expanded by Tegen dust emission scheme
!   Marina Astitha, 2009-2012: expanded by a new dust emissions scheme
!   Stefanie Falk, 2017: include bromine emission scheme from snow and sea ice

MODULE messy_onemis

  USE messy_main_constants_mem, ONLY: dp, pi, N_A, M_air
  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: ABS, ACOS, ADJUSTL, ASSOCIATED, COS, EXP, INDEX, INT &
       , MAX, MIN, MOD, REAL, SIN, SIZE , SQRT, TRIM

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'onemis'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.1.5'

  ! online emis types
  !##############################
  ! ### add new emission here ###
  !##############################
  INTEGER, PUBLIC, PARAMETER :: max_emis = 25

  CHARACTER(LEN=12), DIMENSION(max_emis), PUBLIC :: EMIS_TYPE = ''
  ! ALL INCLUDED EMISSION TYPES
  !  CHARACTER(LEN=12), DIMENSION(max_emis), PARAMETER :: &
  !       emis_type = (/'DMS         ','OC/BC       ','SS_lsce     '&
  !                    ,'SS_monahan  ','SS_aerocom  ','O3ice       '&
  !                    ,'CH4         ','VOC         ','NO          '&
  !                    ,'NOpls       ','DU          ','SO2_ant     '&
  !                    ,'            ','NO_yl95sl10 ','terr13C     '&
  !                    ,'BIOO        ','BIOM        ','DU_tegen    '&
  !                    ,'SS_POC_AQUA ','SS_POC_SWIFS','SS_WIOC_AQUA'&
  !                    ,'SS_WIOC_BLEN','DU_Astitha1 ','DU_Astitha2 '&
  !                    ,'AirSnow'/)

  namelist /CTRL/ EMIS_TYPE

  ! cy_ma_20090901+
  ! parameter needed for new dust emissions (M.ASTITHA)
  INTEGER, PUBLIC,PARAMETER :: nbiomes   = 6   !
  INTEGER, PUBLIC,PARAMETER :: nclasses  = 10  !
  INTEGER, PUBLIC,PARAMETER :: maxsizes  = 250 !
  INTEGER, PUBLIC,PARAMETER :: DU_A2_npt = 8   ! Number of transport size bins
  INTEGER, PUBLIC,PARAMETER :: DU_A2_nps = 4   ! Number of source    size bins
  ! cy_ma_20090901-

  ! parameter needed for VOC/NO emissions
  INTEGER, PUBLIC,PARAMETER :: nveglay_hr=4
  ! parameters needed for NO emissions
  ! 12 emission classes in YL95 inventory
  INTEGER, PUBLIC, PARAMETER :: ncl_yl95=12
  ! number of dry days needed to get pul
  INTEGER, PUBLIC, PARAMETER :: ndrydays=14
  ! tropical rainforest class
  INTEGER, PUBLIC            :: trop_class=0

  ! diurnal cycle parameterisation
  INTEGER, PUBLIC, PARAMETER :: ndiurn = 4

  ! NO emission factor for the twelve ecosystems, wet conditions
  REAL(dp), PUBLIC, DIMENSION(ncl_yl95) :: noemfact_wet = &
       (/0._dp,0._dp,0._dp,0._dp,0.05_dp,0.36_dp,0.17_dp  &
       ,0.03_dp,0.03_dp,0.06_dp,2.6_dp,0._dp/)
  ! NO emission factor for the twelve ecosystems, dry conditions
  REAL(dp), PUBLIC, DIMENSION(ncl_yl95) :: noemfact_dry = &
       (/0._dp,0._dp,0._dp,0._dp,0.37_dp,2.65_dp,1.44_dp &
       ,0.22_dp,0.22_dp,0.40_dp,8.6_dp,0._dp/)

  ! mz_js_20081021+
  INTEGER, PUBLIC, PARAMETER :: ncl_yl95sl10 = 24 ! MODIS+Koeppen
  INTEGER, PUBLIC            :: smoist_method = 2 ! 0: precipitation history
                                                  ! 1: soil water column
                                                  ! 2: volumetric soil moisture
  !                                               !    (default now)
  ! NO emission factor for the 24 ecosystems, wet conditions
  REAL (dp), PUBLIC, DIMENSION(ncl_yl95sl10) :: noemfact_wet_yl95sl10 = &
       (/0., 0.,          &
        0., 0., 0.,       &
        0.,               &
        0., 0.,           &
        0.05, 0.05, 0.05, &
        0.36, 0.36,       &
        0.17,             &
        0.03, 0.03, 0.03, &
        0.03, 0.03,       &
        0.06,             &
        2.6,              &
        0., 0., 0./)
  ! NO emission factor for the 24 ecosystems, dry conditions
  REAL (dp), PUBLIC, DIMENSION(ncl_yl95sl10) :: noemfact_dry_yl95sl10 = &
       (/0., 0.,          &
        0., 0., 0.,       &
        0.,               &
        0., 0.,           &
        0.37, 0.37, 0.37, &
        2.65, 2.65,       &
        1.44,             &
        0.22, 0.22, 0.22, &
        0.22, 0.22,       &
        0.40,             &
        8.6,              &
        0., 0., 0./)
  ! mz_js_20081021-

  namelist /CTRL_NOsl10/ noemfact_wet_yl95sl10, noemfact_dry_yl95sl10 &
         , smoist_method

  ! um_gg_20090915+
  !
  ! ******* CLASSIFICATION DEPENDENT VEGETATION PARAMETERS *******
  !
  ! combimax    value   int.    number of different biome types
  !
!  INTEGER, PARAMETER :: combimax=29

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !   The 28 different biome types output by BIOME4
  !-----------------------------------------------------------------
  !     1       Tropical evergreen broadleaf forest
  !     2       Tropical semi-evergreen broadleaf forest
  !     3       Tropical deciduous broadleaf forest and woodland
  !     4       Temperate deciduous broadleaf forest
  !     5       Temperate evergreen needleleaf forest
  !     6       Warm-temperate evergreen broadleaf and mixed forest
  !     7       Cool mixed forest
  !     8       Cool evergreen needleleaf forest
  !     9       Cool-temperate evergreen needleleaf and mixed forest
  !     10      Cold evergreen needleleaf forest
  !     11      Cold deciduous forest
  !     12      Tropical savanna
  !     13      Tropical xerophytic shrubland
  !     14      Temperate xerophytic shrubland
  !     15      Temperate sclerophyll woodland and shrubland
  !     16      Temperate deciduous broadleaf savanna
  !     17      Temperate evergreen needleleaf open woodland
  !     18      Cold parkland
  !     19      Tropical grassland
  !     20      Temperate grassland
  !     21      Desert
  !     22      Graminoid and forb tundra
  !     23      Low and high shrub tundra
  !     24      Erect dwarf-shrub tundra
  !     25      Prostrate dwarf-shrub tundra
  !     26      Cushion-forb tundra
  !     27      Barren
  !     28      Ice
  !    (29)     Water (this is implied)

  !
  ! active vegetation types:
  ! ------------------------
  !
!!$  INTEGER :: active(combimax) = (/ &
!!$       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
!!$       1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, &
!!$       1, 1, 1, 0, 0 /)
!!$
!!$  ! Biomes including shrubs (active=1 only)
!!$  ! ------------------------
!!$  !
!!$  REAL(dp) :: shrub(combimax) = (/ &
!!$       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
!!$       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
!!$       1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
!!$       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
!!$       1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /)

  !-------------------------------------------------------------------------
  ! solspe --> SOIL CARACTERISTICS:
  ! -------------------------------
  !
  !  ZOBLER texture classes |
  !
  ! SOLSPE: for 4 populations : values = 3*(Dmed sig p); ratio of fluxes;
  !                                      residual moisture
  !
  !   Populations: Coarse sand, medium/fine sand, Silt, Clay
  !
  !     soil type 1 : Coarse
  !     soil type 2 : Medium
  !     soil type 3 : Fine
  !     soil type 4 : Coarse Medium
  !     soil type 5 : Coarse Fine
  !     soil type 6 : Medium Fine
  !     soil type 7 : Coarse_dp, Medium_dp, Fine
  !     soil type 8 : Organic
  !     soil type 9 : Ice
  !     soil type 10 : Potential Lakes (additional)
  !     soil type 11 : Potential Lakes (clay)
  !     soil type 12 : Potential Lakes Australia
  !
  !-------------------------------------------------------------------------

  REAL(dp) :: solspe(12,14) = &
       RESHAPE ( (/ &
       0.0707_dp, 2.0_dp, 0.43_dp, 0.0158_dp, 2.0_dp, 0.40_dp, 0.0015_dp, &
       2.0_dp, 0.17_dp, 0.0002_dp, 2.0_dp, 0.00_dp, 2.1e-06_dp, 0.20_dp,  &
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.37_dp, 0.0015_dp, &
       2.0_dp, 0.33_dp, 0.0002_dp, 2.0_dp, 0.30_dp, 4.0e-06_dp, 0.25_dp,  &
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.00_dp, 0.0015_dp, &
       2.0_dp, 0.33_dp, 0.0002_dp, 2.0_dp, 0.67_dp, 1.0e-07_dp, 0.50_dp,  &
       0.0707_dp, 2.0_dp, 0.10_dp, 0.0158_dp, 2.0_dp, 0.50_dp, 0.0015_dp, &
       2.0_dp, 0.20_dp, 0.0002_dp, 2.0_dp, 0.20_dp, 2.7e-06_dp, 0.23_dp,  &
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.50_dp, 0.0015_dp, &
       2.0_dp, 0.12_dp, 0.0002_dp, 2.0_dp, 0.38_dp, 2.8e-06_dp, 0.25_dp,  &
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.27_dp, 0.0015_dp, &
       2.0_dp, 0.25_dp, 0.0002_dp, 2.0_dp, 0.48_dp, 1.0e-07_dp, 0.36_dp,  &
       0.0707_dp, 2.0_dp, 0.23_dp, 0.0158_dp, 2.0_dp, 0.23_dp, 0.0015_dp, &
       2.0_dp, 0.19_dp, 0.0002_dp, 2.0_dp, 0.35_dp, 2.5e-06_dp, 0.25_dp,  &
       0.0707_dp, 2.0_dp, 0.25_dp, 0.0158_dp, 2.0_dp, 0.25_dp, 0.0015_dp, &
       2.0_dp, 0.25_dp, 0.0002_dp, 2.0_dp, 0.25_dp, 0.0e-00_dp, 0.50_dp,  &
       0.0707_dp, 2.0_dp, 0.25_dp, 0.0158_dp, 2.0_dp, 0.25_dp, 0.0015_dp, &
       2.0_dp, 0.25_dp, 0.0002_dp, 2.0_dp, 0.25_dp, 0.0e-00_dp, 0.50_dp,  &
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.00_dp, 0.0015_dp, &
       2.0_dp, 1.00_dp, 0.0002_dp, 2.0_dp, 0.00_dp, 1.0e-05_dp, 0.25_dp,  &
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.00_dp, 0.0015_dp, &
       2.0_dp, 0.00_dp, 0.0002_dp, 2.0_dp, 1.00_dp, 1.0e-05_dp, 0.25_dp,  &
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.00_dp, 0.0027_dp, &
       2.0_dp, 1.00_dp, 0.0002_dp, 2.0_dp, 0.00_dp, 1.0e-05_dp, 0.25_dp   &
       /), (/ 12, 14 /), ORDER=(/2,1/))


  !**** parameters as defined in mo_physc2 in ECHAM5:
  !  *inverse of equivalent water height when snow is considered to cover
  !   completely the ground in the box.
  REAL(dp), PARAMETER :: cqsncr = 0.95_dp
  !  *maximum moisture content of the skin reservoir
  REAL(dp), PARAMETER :: cwlmax = 2.E-4_dp


  ! *** dust parameter
  REAL(dp), PARAMETER     :: Dmin=0.00002_dp  ! minimum particules diameter (cm)
  REAL(dp), PARAMETER     :: Dmax=0.130_dp    ! maximum particules diameter (cm)
  ! diameter increment (cm)
  REAL(dp), PARAMETER     :: Dstep=0.0460517018598807_dp

  INTEGER, PARAMETER      :: ntrace=8            ! number of tracers
  INTEGER, PARAMETER      :: nbin=24             ! number of bins per tracer
  INTEGER, PARAMETER      :: nclass=ntrace*nbin  ! number of particle classes
  INTEGER, PARAMETER      :: nats =12            ! number of soil types
  INTEGER, PARAMETER      :: nmode=4

  REAL(dp), PARAMETER     :: roa=0.001227_dp     ! air density (g/cm-3)
  REAL(dp), PARAMETER     :: rop=2.65_dp         ! particle density (g/cm-3)

  ! roughness length parameters...(see: Thesis of B. Marticorena)
  REAL(dp), PARAMETER     :: z01=0.001_dp
  REAL(dp), PARAMETER     :: z02=Z01
  ! Scale factor for wind stress threshold to adjust dust emission occurence
  ! (Tegen et al., 2004)
  REAL(dp), PUBLIC        :: cuscale_in = 0.86_dp ! um_gg_20130502
  LOGICAL,  PUBLIC        :: l_nudging  = .FALSE. ! um_gg_20130625

  namelist /CTRL_DU/ cuscale_in, l_nudging ! um_gg_20130625

  ! Logical to calculate dust composition in Astitha 1/2 emissions scheme
  LOGICAL,  PUBLIC        :: l_ducomp = .TRUE.    !  mz_ap_20180214

  !kit_sf_20170206
  ! Bromine explosion default values according to Toyota et al. 2011 parameterization
  INTEGER, PUBLIC         :: r_temp_crit = -15      ! Critical temperature [deg celsius]
  INTEGER, PUBLIC         :: r_sun_theta_crit = 85  ! Critical sun zenith angle [deg]
  ! Efficiency of bromine release due to ozone deposit ('dark','sunlit','land')
  REAL(dp), PUBLIC, DIMENSION(3)     :: r_trigger_1 = (/ 0.001_dp, 0.075_dp, 0.0_dp /)

  namelist /CTRL_AirSnow/ r_temp_crit, r_sun_theta_crit, r_trigger_1 ! kit_sf_20170206

  REAL(dp), DIMENSION(nats,nclass)   :: srel, srelV, su_srelV
  REAL(dp), DIMENSION(nclass)        :: Uth
  REAL(dp)                           :: rdp, c_eff

  PUBLIC :: onemis_read_nml_ctrl
  PUBLIC :: onemis_read_nml_ctrl_NOsl10
  PUBLIC :: onemis_read_nml_ctrl_DU
  PUBLIC :: onemis_read_nml_ctrl_AirSnow ! kit_sf_20170206

  PUBLIC :: parse_f2tstr
  PUBLIC :: parse_rgtstr
  !
  PUBLIC :: dms_emissions
  PUBLIC :: seasalt_emissions_lsce
  PUBLIC :: seasalt_emissions_monahan
  PUBLIC :: seasalt_emissions_aerocom
  PUBLIC :: carbon_emissions
  PUBLIC :: O3ice_emissions
  PUBLIC :: CH4_emissions
  PUBLIC :: VOC_emissions
  PUBLIC :: NOpls_emissions
  PUBLIC :: NO_emissions
  PUBLIC :: dust_emissions
  PUBLIC :: dust_tegen_init               ! um_gg_20100326
  PUBLIC :: dust_emissions_tegen          ! um_gg_20090910
  PUBLIC :: so2_emissions
  PUBLIC :: dust_emissions_DU_Astitha1    ! cy_ma_20090901
  PUBLIC :: dust_emissions_DU_Astitha2    ! cy_ma_20090901
  PUBLIC :: dust_emissions_KKDU_Astitha1  ! mz_kk_20170614
  PUBLIC :: NOemis_yl95sl10_pulsing       ! mz_js_20081006
  PUBLIC :: NOemis_yl95sl10               ! mz_js_20081006
  PUBLIC :: NOemis_yl95sl10_crf           ! mz_js_20081006

  PUBLIC :: bioaer_emissions_olson        ! mz_sb_20091217
  PUBLIC :: bioaer_emissions_modis        ! mz_sb_20091217

  PUBLIC :: bio_emissions_lai             ! mz_mt_20160324
  PUBLIC :: terr13C_frac                  ! mz_sg_20091125
  PUBLIC :: POC_emis_SS                   ! mz_ht_20110512
  PUBLIC :: WIOC_EMIS_SS                  ! mz_ht_20110512
  PUBLIC :: airsnow_emissions             ! kit_sf_20170207

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE onemis_read_nml_ctrl(status, iou)

    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Sep 2003 (for tropop)
    !         Astrid Kerkweg, Uni-Mainz 2010

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'onemis_read_nml_ctrl'
    LOGICAL                     :: lex
    INTEGER                     :: fstat

    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! OK

  END SUBROUTINE onemis_read_nml_ctrl

  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------

  SUBROUTINE onemis_read_nml_ctrl_NOsl10(status, iou)

    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Sep 2003 (for tropop)
    !         Astrid Kerkweg, Uni-Mainz 2010

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'onemis_read_nml_ctrl'
    LOGICAL                     :: lex
    INTEGER                     :: fstat

    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL_NOsl10', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_NOsl10, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_NOsl10', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! OK

  END SUBROUTINE onemis_read_nml_ctrl_NOsl10

  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------

  SUBROUTINE onemis_read_nml_ctrl_DU(status, iou)

    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Sep 2003 (for tropop)
    !         Gregor Glaeser, Uni-Mainz 2013

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'onemis_read_nml_ctrl'
    LOGICAL                     :: lex
    INTEGER                     :: fstat

    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL_DU', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_DU, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_DU', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! OK

  END SUBROUTINE onemis_read_nml_ctrl_DU

! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------

  SUBROUTINE onemis_read_nml_ctrl_AirSnow(status, iou)

    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Sep 2003 (for tropop)
    !         Gregor Glaeser, Uni-Mainz 2013

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'onemis_read_nml_ctrl'
    LOGICAL                     :: lex
    INTEGER                     :: fstat

    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL_AirSnow', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_DU, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_AirSnow', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! OK

  END SUBROUTINE onemis_read_nml_ctrl_AirSnow
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE parse_f2tstr(status, strlen, str &
       , name, method, factor)


    ! AUTHOR:  Patrick Joeckel, MPICH, 2004
    !          Astrid  Kerkweg, MPICH, 2004 (extended for needs of onemis)
    USE messy_main_tools,         ONLY: strcrack
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT)   :: status   ! status information
    INTEGER,          INTENT(IN)    :: strlen   ! max. length of strings
    CHARACTER(LEN=*), INTENT(IN)    :: str      ! string to parse
    CHARACTER(LEN=2*STRLEN_MEDIUM+1), DIMENSION(:), POINTER :: name ! (OUT)
    INTEGER,                      DIMENSION(:), POINTER :: method ! INTENT(OUT)
    REAL(DP),                     DIMENSION(:), POINTER :: factor ! INTENT(OUT)
    ! LOCAL
    CHARACTER(LEN=strlen),               POINTER     :: sl1(:)
    CHARACTER(LEN=strlen),               POINTER     :: sl2(:)
    CHARACTER(LEN=strlen),               POINTER     :: sl3(:)
    CHARACTER(LEN=strlen),               POINTER     :: sl4(:)
    INTEGER :: n, m, l, k
    INTEGER :: i, j
    INTEGER :: ix
    INTEGER :: iostat

    status = 1 ! ERROR

    NULLIFY(sl1)
    NULLIFY(sl2)
    NULLIFY(sl3)
    NULLIFY(sl4)

    CALL strcrack(str, ';', sl1, n)

    ALLOCATE(name(n))
    name(:) = ''
    ALLOCATE(method(n))
    method(:) = 0
    ALLOCATE(factor(n))
    factor = 1.0_DP

    DO i=1, n

       CALL strcrack(sl1(i), ':', sl2, m)
       SELECT CASE(m)
       CASE(0)
          ! CANNOT BE REACHED
       CASE(1)
          ! NO VALID ':' -> ONLY TRACER NAME OR ONLY SPECIFICATION
          ix = INDEX(TRIM(sl2(1)), ',') + INDEX(TRIM(sl2(1)), '=')
          IF ((TRIM(ADJUSTL(sl2(1))) == '' ) .OR. (ix > 0)) THEN
             status = 1 ! EMPTY TRACER NAME
          ELSE
             status = 0
             name(i) = TRIM(ADJUSTL(sl2(1)))
          END IF
          RETURN
       CASE(2)
          ! ONE ':' -> FULL FEATURED
          name(i) = TRIM(ADJUSTL(sl2(1)))
          ! -> GO ON BELOW
       CASE DEFAULT
          ! MORE THAN ONE ':' -> ERROR
          IF (m > 2 ) THEN
             status = 2
             RETURN
          END IF
       END SELECT

       CALL strcrack(sl2(2), ',', sl3, l)
       SELECT CASE(l)
       CASE(0)
          ! CANNOT BE REACHED
       CASE(1,2)
          ! ONE OR TWO ',' -> FULL FEATURED (GO ON BELOW)
       CASE DEFAULT
          ! MORE THAN ONE ',' (PLUS ONE OPTIONAL AT THE END)
          status = 3
          RETURN
       END SELECT

       DO j=1, l
          CALL strcrack(sl3(j), '=', sl4, k)
          SELECT CASE(k)
          CASE(0)
             ! CANNOT BE REACHED
          CASE(1)
             ! NO '=' OR ONLY ONE '=' AT THE END -> ERROR
             status = 4
             RETURN
          CASE(2)
             ! -> FULL FEATURED: GO ON BELOW
          CASE DEFAULT
             ! MORE THAN ONE '=' -> ERROR
             status = 5
             RETURN
          END SELECT

          SELECT CASE(TRIM(ADJUSTL(sl4(1))))
          CASE ('M')
             READ(sl4(2),*,IOSTAT=iostat) method(i)
             IF (iostat /= 0) THEN
                status = 6  ! ERROR IN READING INTEGER
                RETURN
             END IF
          CASE ('SC')
             READ(sl4(2),*,IOSTAT=iostat) factor(i)
             IF (iostat /= 0) THEN
                status = 7  ! ERROR IN READING REAL
                RETURN
             END IF
          CASE DEFAULT
             status = 8
             RETURN
          END SELECT
       END DO

    END DO

    ! CLEAN UP
    IF (ASSOCIATED(sl1)) DEALLOCATE(sl1)
    IF (ASSOCIATED(sl2)) DEALLOCATE(sl2)
    IF (ASSOCIATED(sl3)) DEALLOCATE(sl3)
    IF (ASSOCIATED(sl4)) DEALLOCATE(sl4)

    status = 0 ! NO ERROR

  END SUBROUTINE parse_f2tstr
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE parse_rgtstr(status, strlen, str, nml, var, file, &
       type, lrgt)

    ! AUTHOR:  Patrick Joeckel, MPICH, 2004
    !          Astrid  Kerkweg, MPICH, 2004 (extended for needs of onemis)
    USE messy_main_tools, ONLY: strcrack

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT)   :: status   ! status information
    INTEGER,          INTENT(IN)    :: strlen   ! max. length of strings
    CHARACTER(LEN=*), INTENT(IN)    :: str      ! string to parse
    CHARACTER(LEN=*), INTENT(INOUT) :: nml      ! namelist file
    CHARACTER(LEN=*), INTENT(INOUT) :: var      ! netCDF variable
    CHARACTER(LEN=*), INTENT(INOUT) :: file     ! netCDF file
    CHARACTER(LEN=*), INTENT(INOUT) :: type     ! emission type e.g OC, SS etc.
    LOGICAL,          INTENT(OUT)   :: lrgt     ! realy regrid event ?
    ! LOCAL
    CHARACTER(LEN=strlen),               POINTER     :: sl1(:)
    CHARACTER(LEN=strlen),               POINTER     :: sl2(:)
    INTEGER :: n, m
    INTEGER :: i

    status = 1 ! ERROR
    lrgt = .true.

    NULLIFY(sl1)
    NULLIFY(sl2)

    CALL strcrack(str, ';', sl1, n)
    DO i=1, n

       CALL strcrack(sl1(i), '=', sl2, m)
       IF (SIZE(sl2) == 2) THEN
          IF (TRIM(ADJUSTL(sl2(2))) == '') THEN
             status = 50    ! EMPTY SPECIFICATION
             RETURN
          END IF
       END IF
       IF (m > 2 ) THEN
          status = 5    ! to much '=' specifications
          RETURN
       END IF

       SELECT CASE(TRIM(ADJUSTL(sl2(1))))
          CASE('NML')
             nml = TRIM(ADJUSTL(sl2(2)))
          CASE('VAR')
             var = TRIM(ADJUSTL(sl2(2)))
          CASE('FILE')
             file = TRIM(ADJUSTL(sl2(2)))
          CASE('TYPE')
             type = TRIM(ADJUSTL(sl2(2)))
          CASE('NO_RGT')
             lrgt = .FALSE.
          CASE DEFAULT
             status = 8 ! UNKNOWN SPECIFIER
             RETURN
       END SELECT

    END DO

    ! CLEAN UP
    IF (ASSOCIATED(sl1)) DEALLOCATE(sl1)
    IF (ASSOCIATED(sl2)) DEALLOCATE(sl2)

    status = 0 ! NO ERROR

  END SUBROUTINE parse_rgtstr
  ! -------------------------------------------------------------------------

  ! --------------------------------------------------------------------------

  PURE SUBROUTINE dms_emissions &
    (emis_dms_sea, seawater_dms, wind10, tsw, slm, seaice)

    ! Calculate DMS emissions following Liss & Merlivat (1986)
    ! based on code from Philip Stier, MPI-MET

    USE messy_main_constants_mem, ONLY: N_A

    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(OUT) :: emis_dms_sea(:) ! DMS emission flux [mcl m-2 s-1]
    REAL(dp), INTENT(IN)  :: seawater_dms(:) ! DMS seawater conc [nmol/L]
    REAL(dp), INTENT(IN)  :: wind10(:)       ! wind speed at 10 m [m/s]
    REAL(dp), INTENT(IN)  :: tsw(:)          ! surface T over water [K]
    REAL(dp), INTENT(IN)  :: slm(:)          ! sea-land mask [0=water,1=land]
    REAL(dp), INTENT(IN)  :: seaice(:)       ! sea ice fraction [0=liq,1=ice]

    ! LOCAL
    INTEGER :: jp
    REAL(dp):: sst_C   ! sea-surface temperature [Celsius]
    REAL(dp):: schmidt ! Schmidt number after Andreae [1]
    REAL(dp):: zwind10 ! wind speed at 10 m [m/s]
    REAL(dp):: zvdms   ! piston velocity [cm/h]

    DO jp=1,SIZE(emis_dms_sea)
      sst_C=tsw(jp)-273.15_dp
      sst_C = MIN(sst_C, 35.0_dp)
      schmidt=3652.047271_dp-246.99*sst_C+8.536397*sst_C*sst_C     &
        -0.124397*sst_C*sst_C*sst_C
      ! Calculate ocean atmosphere exchange (Liss & Merlivat, 1986)
      zwind10=wind10(jp)
      IF ( (zwind10 > 3.6_dp) .AND. (zwind10 <= 13._dp) ) THEN
        zvdms = (2.85*zwind10-9.65_dp) * (schmidt/600.)**(-0.5)
      ELSE IF(zwind10 <= 3.6_dp) THEN
        zvdms = (0.17*zwind10)      * (schmidt/600.)**(-2._dp/3.)
      ELSE
        zvdms = (5.9*zwind10-49.3_dp)  * (schmidt/600.)**(-0.5)
      END IF
      ! DMS emission flux [mcl m-2 s-1]
      emis_dms_sea(jp) = seawater_dms(jp) * (zvdms*1.E-2/3600.) * &
        1.E-6*N_A * &
        (1._dp-slm(jp)) * (1._dp-seaice(jp))
    END DO

  END SUBROUTINE dms_emissions

  ! --------------------------------------------------------------------------

  PURE SUBROUTINE seasalt_emissions_lsce( &
    wind10, slf, seaice, alake, mss_as, mss_cs, nss_as, nss_cs)

    ! Author:
    ! -------
    ! Michael Schulz
    ! Laboratoire des Sciences du Climat et de l'Environnement / Saclay
    ! 10.1.2002
    !
    ! Modifications:
    ! --------------
    ! Philip Stier, MPI-MET  (Adaption to the ECHAM/HAM structure)       2002
    ! Michael Schulz, LSCE   (Modified source coefficients)        02/08/2002
    ! Rolf Sander, MPICH, 2004: rewritten as MESSy submodel
    ! Astrid Kerkweg, MPICH, 2004 : change in lake mask added
    !
    ! Purpose:
    ! --------
    ! Describe source flux of sea salt aerosol mass and number flux
    ! as a function of wind speed
    ! for two aerosol modes: coarse soluble and accumulation soluble
    !
    ! Interface:
    ! ----------
    ! input
    ! wind10: wind speed at 10 m [m/s]
    ! slf:    sea land fraction  0 = sea, 1 = land
    ! seaice: sea ice fraction   0 = liquid surface, 1 = ice-covered surface
    ! alake:  lake fraction of grid box
    !
    ! output
    ! mss_as, mss_cs: emission fluxes of sea salt mass   [kg m-2 s-1]
    ! nss_as, nss_cs: emission fluxes of sea salt number [# m-2 s-1]
    !
    ! Method:
    ! -------
    ! Tabulated mass and number fluxes following Monahan 86 and
    ! Smith&Harrison 98 for two aerosol modes; interpolated according to
    ! actual wind speed; sea salt mass and number fluxes are fitted to two
    ! lognormal distributions following work published by Guelle et al. (2001)
    ! JGR 106, pp. 27509-27524


    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(in)  :: wind10(:), slf(:), seaice(:), alake(:)
    REAL(dp), INTENT(out) :: mss_as(:), mss_cs(:), nss_as(:), nss_cs(:)

    ! LOCAL
    INTEGER, PARAMETER :: NWCL=40 ! number of windclasses
    INTEGER :: wcl(SIZE(wind10))
    REAL(dp):: landfrac(SIZE(wind10))
    REAL(dp):: dv(SIZE(wind10))
    REAL(dp):: waterfrac(SIZE(wind10))

    ! precalculated mass fluxes per wind class [kg m-2 s-1]:

    REAL, PARAMETER :: mass1flux(0:NWCL) = (/                           &
      0.000E+00, 2.483E-15, 2.591E-14, 1.022E-13, 2.707E-13, 5.761E-13, &
      1.068E-12, 1.800E-12, 2.829E-12, 4.215E-12, 6.023E-12, 8.317E-12, &
      1.117E-11, 1.464E-11, 1.882E-11, 2.378E-11, 2.959E-11, 3.633E-11, &
      4.409E-11, 5.296E-11, 6.301E-11, 7.433E-11, 8.693E-11, 1.012E-10, &
      1.168E-10, 1.342E-10, 1.532E-10, 1.741E-10, 1.970E-10, 2.219E-10, &
      2.489E-10, 2.781E-10, 3.097E-10, 3.437E-10, 3.803E-10, 4.195E-10, &
      4.616E-10, 5.065E-10, 5.544E-10, 6.054E-10, 6.711E-10             /)

    REAL, PARAMETER :: mass2flux(0:NWCL) = (/                           &
      0.000E+00, 2.319E-13, 2.411E-12, 9.481E-12, 2.505E-11, 5.321E-11, &
      9.850E-11, 1.658E-10, 2.602E-10, 3.874E-10, 5.529E-10, 7.628E-10, &
      1.023E-09, 1.341E-09, 1.722E-09, 2.175E-09, 2.704E-09, 3.319E-09, &
      4.026E-09, 4.832E-09, 5.746E-09, 6.776E-09, 7.925E-09, 9.214E-09, &
      1.064E-08, 1.221E-08, 1.394E-08, 1.584E-08, 1.791E-08, 2.016E-08, &
      2.261E-08, 2.526E-08, 2.812E-08, 3.120E-08, 3.451E-08, 3.806E-08, &
      4.186E-08, 4.592E-08, 5.025E-08, 5.486E-08, 6.014E-08             /)

    !--- Precalculated number fluxes per wind class [m-2 s-1]:

    REAL, PARAMETER :: numb1flux(0:NWCL) = (/                           &
      0.000E+00, 3.004E+01, 3.245E+02, 1.306E+03, 3.505E+03, 7.542E+03, &
      1.410E+04, 2.394E+04, 3.787E+04, 5.674E+04, 8.147E+04, 1.130E+05, &
      1.523E+05, 2.005E+05, 2.586E+05, 3.278E+05, 4.091E+05, 5.037E+05, &
      6.129E+05, 7.379E+05, 8.800E+05, 1.041E+06, 1.220E+06, 1.422E+06, &
      1.646E+06, 1.893E+06, 2.166E+06, 2.466E+06, 2.794E+06, 3.152E+06, &
      3.541E+06, 3.962E+06, 4.419E+06, 4.911E+06, 5.441E+06, 6.011E+06, &
      6.621E+06, 7.274E+06, 7.972E+06, 8.716E+06, 8.801E+06             /)

    REAL, PARAMETER :: numb2flux(0:NWCL) = (/                           &
      0.000E+00, 1.934E+01, 2.068E+02, 8.271E+02, 2.211E+03, 4.741E+03, &
      8.841E+03, 1.497E+04, 2.363E+04, 3.534E+04, 5.066E+04, 7.017E+04, &
      9.447E+04, 1.242E+05, 1.600E+05, 2.025E+05, 2.525E+05, 3.106E+05, &
      3.776E+05, 4.542E+05, 5.413E+05, 6.395E+05, 7.501E+05, 8.726E+05, &
      1.009E+06, 1.160E+06, 1.327E+06, 1.509E+06, 1.709E+06, 1.927E+06, &
      2.163E+06, 2.420E+06, 2.697E+06, 2.996E+06, 3.318E+06, 3.664E+06, &
      4.034E+06, 4.430E+06, 4.852E+06, 5.303E+06, 5.740E+06             /)

    !--- Precalculated mass flux gradient for each wind class for
    !    interpolation [kg m-3]: dm/dv(i) where m = m(i) + dm/dv(i) * dv(i)

    REAL, PARAMETER :: dmass1flux(0:NWCL) = (/                          &
      2.483E-15, 2.343E-14, 7.630E-14, 1.684E-13, 3.054E-13, 4.919E-13, &
      7.319E-13, 1.029E-12, 1.386E-12, 1.807E-12, 2.294E-12, 2.850E-12, &
      3.477E-12, 4.174E-12, 4.960E-12, 5.810E-12, 6.745E-12, 7.762E-12, &
      8.863E-12, 1.005E-11, 1.132E-11, 1.260E-11, 1.423E-11, 1.569E-11, &
      1.733E-11, 1.907E-11, 2.090E-11, 2.284E-11, 2.487E-11, 2.700E-11, &
      2.924E-11, 3.158E-11, 3.403E-11, 3.659E-11, 3.925E-11, 4.202E-11, &
      4.490E-11, 4.790E-11, 5.100E-11, 6.578E-11, 6.578E-11             /)

    REAL, PARAMETER :: dmass2flux(0:NWCL) = (/                          &
      2.319E-13, 2.179E-12, 7.070E-12, 1.557E-11, 2.817E-11, 4.528E-11, &
      6.727E-11, 9.446E-11, 1.271E-10, 1.655E-10, 2.100E-10, 2.606E-10, &
      3.177E-10, 3.812E-10, 4.522E-10, 5.297E-10, 6.145E-10, 7.068E-10, &
      8.066E-10, 9.142E-10, 1.030E-09, 1.149E-09, 1.289E-09, 1.425E-09, &
      1.573E-09, 1.730E-09, 1.896E-09, 2.070E-09, 2.254E-09, 2.447E-09, &
      2.649E-09, 2.860E-09, 3.080E-09, 3.311E-09, 3.550E-09, 3.800E-09, &
      4.060E-09, 4.329E-09, 4.609E-09, 5.279E-09, 5.279E-09             /)

    !--- Precalculated number flux gradient for each wind class for
    !    interpolation [m-3]: dn/dv(i) where n = n(i) + dn/dv(i) * dv(i)

    REAL, PARAMETER :: dnumb1flux(0:NWCL) = (/                          &
      3.004E+01, 2.945E+02, 9.811E+02, 2.200E+03, 4.036E+03, 6.562E+03, &
      9.839E+03, 1.393E+04, 1.887E+04, 2.473E+04, 3.153E+04, 3.935E+04, &
      4.818E+04, 5.808E+04, 6.914E+04, 8.130E+04, 9.465E+04, 1.092E+05, &
      1.250E+05, 1.421E+05, 1.605E+05, 1.798E+05, 2.017E+05, 2.237E+05, &
      2.476E+05, 2.729E+05, 2.997E+05, 3.280E+05, 3.577E+05, 3.890E+05, &
      4.219E+05, 4.564E+05, 4.923E+05, 5.301E+05, 5.695E+05, 6.105E+05, &
      6.531E+05, 6.976E+05, 7.437E+05, 8.550E+04, 8.550E+04             /)

    REAL, PARAMETER :: dnumb2flux(0:NWCL) = (/                          &
      1.934E+01, 1.875E+02, 6.203E+02, 1.384E+03, 2.530E+03, 4.100E+03, &
      6.132E+03, 8.659E+03, 1.171E+04, 1.532E+04, 1.951E+04, 2.430E+04, &
      2.972E+04, 3.582E+04, 4.251E+04, 4.997E+04, 5.812E+04, 6.700E+04, &
      7.663E+04, 8.702E+04, 9.820E+04, 1.106E+05, 1.225E+05, 1.366E+05, &
      1.511E+05, 1.664E+05, 1.826E+05, 1.997E+05, 2.177E+05, 2.366E+05, &
      2.565E+05, 2.773E+05, 2.990E+05, 3.218E+05, 3.455E+05, 3.702E+05, &
      3.959E+05, 4.226E+05, 4.504E+05, 4.376E+05, 4.376E+05             /)

    ! calculate fraction of non salty water (land + lakes)
    WHERE (slf(:)+alake(:) >= 0.95_dp)
       landfrac(:) = 1._dp
    ELSEWHERE
       landfrac(:) = slf(:) + alake(:)
    ENDWHERE

    ! windclass in 1 m/s steps as a function of wind speed:
    wcl(:) = MAX(0,MIN(INT(wind10),NWCL))

    ! dv for interpolation in the windclass interval:
    dv(:) = wind10(:) - REAL(wcl(:),dp)

    ! fraction of the grid cell of non ice-covered water:
    waterfrac(:) = MAX(0._dp,MIN((1._dp-landfrac(:))*(1._dp-seaice(:)),1._dp))

    ! mass flux,   accumulation mode soluble [kg m-2 s-1]:
    mss_as(:) = (REAL(mass1flux(wcl),dp) + REAL(dmass1flux(wcl),dp)*dv) &
                  * waterfrac(:)
    ! number flux, accumulation mode soluble [m-2 s-1]:
    nss_as(:) = (REAL(numb1flux(wcl),dp) + REAL(dnumb1flux(wcl),dp)*dv) &
                  * waterfrac(:)

    ! mass flux,   coarse mode soluble       [kg m-2 s-1]:
    mss_cs(:) = (REAL(mass2flux(wcl),dp) + REAL(dmass2flux(wcl),dp)*dv) &
                  * waterfrac(:)
    ! number flux, coarse mode soluble       [m-2 s-1]:
    nss_cs(:) = (REAL(numb2flux(wcl),dp) + REAL(dnumb2flux(wcl),dp)*dv) &
                  * waterfrac(:)

  END SUBROUTINE seasalt_emissions_lsce

  ! --------------------------------------------------------------------------

  PURE SUBROUTINE seasalt_emissions_monahan( &
    wind10, slf, seaice, alake, mss_as, mss_cs, nss_as, nss_cs)

    ! calculate sea salt flux from the 10 m wind speed following
    ! equation (25) in Monahan et al., 1986 using the sea salt
    ! accumulation and coarse mode.

    ! Authors:
    ! Philip Stier, MPIMet
    ! Rolf Sander,  MPICH, 2004: completely rewritten as MESSy submodel
    ! Astrid Kerkweg, MPICH, 2004: change of slm => slf + alake
    !
    ! Interface:
    ! ----------
    ! wind10    : 10 m wind speed [m s-1]
    ! slf       : sea land fraction  0 = sea, 1 = land
    ! alake     : lake fraction  in grid box
    ! seaice    : seaice fraction [0-1]

    USE messy_main_constants_mem, ONLY:  pi

    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(in)  :: slf(:), wind10(:), seaice(:), alake(:)
    REAL(dp), INTENT(out) :: mss_as(:), mss_cs(:), nss_as(:), nss_cs(:)

    !--- Local Variables:
    REAL(dp) :: f_speed(SIZE(wind10))
    REAL(dp) :: waterfrac(SIZE(wind10))
    REAL(dp) :: landfrac(SIZE(wind10))

    ! qqq the values of mass_as and mass_cs (from
    ! echam5.1.07_m7_0.18_ori/modules/mo_chem_forcing.f90) are
    ! probably specific for r_as and r_cs
    REAL, PARAMETER :: mass_as= 0.37E-15 ! particle mass [kg]
    REAL, PARAMETER :: mass_cs= 0.37E-12 ! particle mass [kg]

    REAL(dp), PARAMETER :: r_as=0.416_dp, r_cs=3.49_dp ! radius [um]
    ! dr is probably the radius range [um]
    REAL(dp), PARAMETER :: dr_as=0.5_dp,  dr_cs=4.5_dp
    ! the values of work_as and work_cs are pre-calculated and then defined
    ! as PARAMETERs:
    ! B_as=0.58-1.54*LOG10(r_as) ! similar to (0.380 - log10(r)) / 0.650
    ! B_cs=0.58-1.54*LOG10(r_cs) ! similar to (0.380 - log10(r)) / 0.650
    ! work_as=10**(1.19*EXP(-B_as**2))
    ! work_cs=10**(1.19*EXP(-B_cs**2))
    REAL(dp), PARAMETER :: work_as=2.01900451790032_dp
    REAL(dp), PARAMETER :: work_cs=13.0178923709599_dp
    ! 1.15e3 is probably particle density [kg/m3]
    ! 1e-18 is probably conversion um->m cubed
    REAL(dp), PARAMETER :: fac=4./3.*pi*1.15e3_dp*1.e-18_dp
    REAL(dp) :: ssfl_as, ssfl_cs

    !qqq TODO: define ssfl_* as PARAMETERs (or before time loop)
    ssfl_as = fac*(1._dp+0.057_dp*r_as**1.05_dp)*work_as*dr_as
    ssfl_cs = fac*(1._dp+0.057_dp*r_cs**1.05_dp)*work_cs*dr_cs

    ! calculate fraction of non salty water (land + lakes)
    WHERE (slf(:)+alake(:) > 0.99_dp)
       landfrac(:) = 1._dp
    ELSEWHERE
       landfrac(:) = slf(:) + alake(:)
    ENDWHERE

    ! fraction of the grid cell of non ice-covered water:
    waterfrac(:) = MAX(0._dp,MIN((1._dp-landfrac(:))*(1._dp-seaice(:)),1._dp))
    ! f_speed = function of wind speed (cutoff at 20 m/s)
    ! qqq why 1.373? Monahan uses 1.37!
    f_speed(:)= 1.373 * MIN(wind10(:),20._dp)**3.41
    ! mass flux,   accumulation mode soluble [kg m-2 s-1]:
    mss_as(:) = ssfl_as * f_speed(:) * waterfrac(:)
    ! number flux, accumulation mode soluble [m-2 s-1]:
    nss_as(:) = mss_as(:) / mass_as
    ! mass flux,   coarse mode soluble       [kg m-2 s-1]:
    mss_cs(:) = ssfl_cs * f_speed(:) * waterfrac(:)
    ! number flux, coarse mode soluble       [m-2 s-1]:
    nss_cs(:) = mss_cs(:) / mass_cs

  END SUBROUTINE seasalt_emissions_monahan

! -------------------------------------------------------------------------

PURE SUBROUTINE seasalt_emissions_aerocom(mss_as, mss_cs, nss_as, nss_cs, &
     numflx_as, numflx_cs, massflx_as, massflx_cs, boxarea)


! incoming mass fluxes in kg/gridbox/day
REAL(dp), INTENT(IN)  :: massflx_as(:), massflx_cs(:)
! incoming number fluxes in 1/gridbox/day
REAL(dp), INTENT(IN)  :: numflx_as(:), numflx_cs(:)
! area of gridbox
REAL(dp), INTENT(IN)  :: boxarea(:)
! outgoing mass fluxes in kg m-2 s-1
REAL(dp), INTENT(OUT) :: mss_as(:), mss_cs(:)
! outgoing number fluxes in m-2 s-1
REAL(dp), INTENT(OUT) :: nss_as(:), nss_cs(:)

! LOCAL
REAL(dp) :: factor(SIZE(boxarea))

! calculate factor for unit change: (gridbox*day) => m2 s
factor(:) = boxarea(:) * 86400._dp

mss_as(:) = massflx_as(:) / factor(:)
mss_cs(:) = massflx_cs(:) / factor(:)
nss_as(:) = numflx_as(:)  / factor(:)
nss_cs(:) = numflx_cs(:)  / factor(:)

END SUBROUTINE seasalt_emissions_aerocom

! -------------------------------------------------------------------------

PURE SUBROUTINE POC_EMIS_SS(SS_flux, slf, POC, MPOC_flux)
! Author:
    ! -------
    ! Susannah Burrows
    ! Max Planck Institute for Chemistry, Department of Atmospheric Chemistry
    ! February 2011
    ! Holger Tost
    ! University Mainz, Institute for Physics of the Atmosphere, May 2011
! Purpose:
    ! --------
    ! Describe source flux of sea spray water-insoluble organic carbon as a
    ! function of wind speed and ocean particulate organic carbon,
    ! based on scaling to accumulation mode sea spray emissions.
! Interface:
    ! ----------
    ! input
    ! ss_flux: seasalt emission flux (accumulation mode) [kg m-2 s-1]
    ! slf:    sea land fraction  0 = sea, 1 = land
    ! poc:  particulate organic carbon [mg m-3] in the ocean (satellite product)
    !
    ! output
    ! mpoc: emission flux of submicron POC mass in sea salt   [kg m-2 s-1]
! Method:
    ! --------
    ! The parameterization is based on the assumption that the water-insoluble
    ! organic carbon (WIOC) component of the sea spray aerosol is contributed by
    ! primary sea spray emissions of oceanic particulates, and that these are
    ! essentially proportional to the ocean POC concentration.  It is further
    ! assumed that POC emissions are never more than 76% of total sea spray
    ! emissions (Vignati et al., 2010), but that the concentration of organics in
    ! the submicron sea spray is highly enriched with respect to the bulk sea
    ! water.

  IMPLICIT NONE
  ! I/O
  REAL(dp), INTENT(in)  :: ss_flux(:), slf(:)
  REAL(dp), INTENT(in)  :: poc(:) ! mg m-3 (Satellite data)
  REAL(dp), INTENT(out) :: mpoc_flux(:)

  REAL(dp) :: poc_g(SIZE(SLF))

  WHERE (slf(:) >= 0.95_dp)
    poc_g(:) = poc(:) * 1.e-3_dp
  ELSEWHERE
    ! POC data are set to zero where no data exist (land)
    ! before regridding, therefore along the coasts, the concentrations
    ! have to be adjusted for the fact that these areas are averaged in
    ! during regridding.
    poc_g(:) = poc(:) * 1.e-3_dp / MIN((1._dp - slf(:)),1._dp)
  ENDWHERE

  ! POC flux [g m-2 s-1] =
  mpoc_flux(:) = &
    !                  poc concentration in water [g m-3]
            poc_g(:) &
    !                  * enrichment factor (range ~10--40)
          * 100._dp &
    !                  * sea spray emissions [kg m-2 s-1]
          * ss_flux(:) / 36._dp &
    ! ( --> Density of seawater ca. 1030 kg m-3, density of POC ca. 1000 kg m-3
    !   --> Seasalt concentration approx. 35 g kg-1 ~= 36 kg m-3)
    !                    convert POC to POM
          * 1.8_dp

    ! organic accumulation mode sea spray mass is not more than 76% of total
    ! accumulation mode mass
  mpoc_flux(:) = MIN(ss_flux(:)*.76_dp, mpoc_flux(:))
END SUBROUTINE POC_EMIS_SS

!-----------------------------------------------------------------------------

PURE SUBROUTINE WIOC_EMIS_SS( &
    ss_flux, slf, chlor_a, emis_wioc)
! Author:
    ! -------
    ! Susannah Burrows
    ! Max Planck Institute for Chemistry, Department of Atmospheric Chemistry
    ! February 2011
    ! Holger Tost
    ! University Mainz, Institute for Physics of the Atmosphere, May 2011
! Purpose:
    ! --------
    ! Describe source flux of sea spray water-insoluble organic carbon as a
    ! function of wind speed and ocean particulate organic carbon,
    ! based on scaling to accumulation mode sea spray emissions.
! Interface:
    ! ----------
    ! input
    ! ss_flux: seasalt emission flux (accumulation mode) [kg m-2 s-1]
    ! slf:     sea land fraction  0 = sea, 1 = land
    ! chlor_a: Chlorophyll-a concentration in the ocean water [mg m-3] (satellite product)
    !
    ! output
    ! emis_wioc: emission flux of submicron WIOC mass in sea salt   [kg m-2 s-1]
! Method:
    ! --------
    ! The parameterization is based on the assumption that the water-insoluble
    ! organic carbon (WIOC) component of the sea spray aerosol is contributed by
    ! primary sea spray emissions of oceanic particulates, and that these are
    ! essentially proportional to the ocean POC concentration.  It is further
    ! assumed that POC emissions are never more than 76% of total sea spray
    ! emissions (Vignati et al., 2010), but that the concentration of organics in
    ! the submicron sea spray is highly enriched with respect to the bulk sea
    ! water.
  IMPLICIT NONE
  ! I/O
  REAL(dp), INTENT(in)  :: ss_flux(:), slf(:)
  REAL(dp), INTENT(in)  :: chlor_a(:) ! mg m-3 (Satellite data)
  REAL(dp), INTENT(out) :: emis_wioc(:)

  REAL(dp)  :: scale_wioc(SIZE(slf))
  ! % organic mass =
  !    Chlorophyll-a concentration [mg m-3] * 43.5 + 13.805
  !  (Chl < 1.43 ug m-3)
  !
  ! Vignati et al. (2010): Global scale emission and distribution of
  ! sea-spray aerosol: Sea-salt and organic enrichment. Atmospheric
  ! Environment.
  ! AND (corrected from)
  !
  ! O'Dowd et al. (2008): A combined organic-inorganic sea spray
  ! source function. GRL.

  scale_wioc(:) = (MAX(chlor_a(:),1.43_dp) * 4.3_dp + 13.805_dp) / 100._dp
  ! WIOC flux [kg m-2 s-1] =
  emis_wioc(:) = scale_wioc(:)  &
!                 * sea salt emissions [kg m-2 s-1]
               * ss_flux(:)

END SUBROUTINE WIOC_EMIS_SS

! -------------------------------------------------------------------------
SUBROUTINE carbon_emissions(kproma, pOC_ag,  pOC_ant, pOC_bge, pOC_wf, &
                            pBC_ag,  pBC_ant, pBC_wf, pBC_sum_insol,   &
                            pOC_sum_insol ,pOC_sum_sol, pNemis_ks,     &
                            pNemis_ki, pCseason,                       &
                            ! mz_ht_20110117+
                            pNemis_ki_bc, pNemis_ki_oc,                &
                            pOC_soa_sol, pOC_bb_sol, pOC_ff_sol,       &
                            pOC_ff_insol, pOC_soa_insol, pOC_bb_insol, &
                            pBC_ff_insol, pBC_bb_insol)
                            ! mz_ht_20110117-

! original code P. Stier, MPI-MET
! aug, 2004 completly rewritten A. Kerkweg, MPICH, Mainz

  INTEGER :: kproma
  REAL(dp),INTENT(IN) :: pOC_ag(1:kproma)
  REAL(dp),INTENT(IN) :: pOC_ant(1:kproma)
  REAL(dp),INTENT(IN) :: pOC_bge(1:kproma)
  REAL(dp),INTENT(IN) :: pOC_wf(1:kproma)
  REAL(dp),INTENT(IN) :: pBC_ag(1:kproma)
  REAL(dp),INTENT(IN) :: pBC_ant(1:kproma)
  REAL(dp),INTENT(IN) :: pBC_wf(1:kproma)
  REAL(dp),INTENT(IN) :: pCseason(1:kproma)
  !carbon emission mass fluxes insoluble part in kg/kg
  REAL(dp),INTENT(OUT) :: pBC_sum_insol(1:kproma)
  REAL(dp),INTENT(OUT) :: pOC_sum_insol(1:kproma)
  !carbon emission mass fluxes soluble part in kg/kg
  REAL(dp),INTENT(OUT) :: pOC_sum_sol(1:kproma)
  !number emission fluxes resulting from carbon emissions in 1 /kg
  REAL(dp),INTENT(OUT) :: pNemis_ks(1:kproma)
  REAL(dp),INTENT(OUT) :: pNemis_ki(1:kproma)

  ! mz_ht_20110301+
  REAL(dp),INTENT(OUT) :: pNemis_ki_oc(1:kproma)
  REAL(dp),INTENT(OUT) :: pNemis_ki_bc(1:kproma)
  REAL(dp),INTENT(OUT) :: pOC_soa_sol(1:kproma)
  REAL(dp),INTENT(OUT) :: pOC_bb_sol(1:kproma)
  REAL(dp),INTENT(OUT) :: pOC_ff_sol(1:kproma)
  REAL(dp),INTENT(OUT) :: pOC_soa_insol(1:kproma)
  REAL(dp),INTENT(OUT) :: pOC_bb_insol(1:kproma)
  REAL(dp),INTENT(OUT) :: pOC_ff_insol(1:kproma)
  REAL(dp),INTENT(OUT) :: pBC_bb_insol(1:kproma)
  REAL(dp),INTENT(OUT) :: pBC_ff_insol(1:kproma)
  ! mz_ht_20110301-

  ! parameter definition
  ! Mass ratio organic species to organic carbon
  ! (Seinfeld and Pandis, 1998, p709;  Ferek et al., JGR, 1998)
  REAL(dp), PARAMETER :: zom2oc         = 1.4_dp
  ! Biom. Burn. Percentage of  Water Soluble OC (WSOC) [1]
  ! (M.O. Andreae; Talk: Smoke and Climate)
  REAL(dp), PARAMETER :: zbb_wsoc_perc  = 0.65_dp
  ! Assume same Percentage of WSOC for biogenic OC
  REAL(dp), PARAMETER :: zbge_wsoc_perc = 0.65_dp


  ! coefficients of mass to number conversion for OC/BC emissions
  REAL(dp) :: zm2n_carb_bb
  REAL(dp) :: zm2n_carb_ff
  REAL(dp) :: zm2n_carb_bg

  ! Fossil fuel emissions: assumed number median radius of the emitted
  ! particles [m]. Has to lie within the Aitken mode for the current setup!
  REAL(dp), PARAMETER ::        cmr_ff         = 0.03E-6_dp

  ! Biomass burning emissions: Assumed number median radius of the emitted
  ! particles with the standard deviation given in[m]. Has to lie within the
  ! Aitken mode for the current setup!
  REAL(dp), PARAMETER ::    cmr_bb         = 0.075E-6_dp

  ! Biogenic secondary particle formation: Assumed number median radius of
  ! the emitted  particles with the standard deviation given in [m].
  ! Has to lie within the Aitken mode for the current setup!

  REAL(dp), PARAMETER ::      cmr_bg         = 0.03E-6_dp

  zm2n_carb_bb = 3./(4.*pi*(cmr_bb)**3.)
  zm2n_carb_bg = 3./(4.*pi*(cmr_bg)**3.)
  zm2n_carb_ff = 3./(4.*pi*(cmr_ff)**3.)

  ! black carbon insoluble fraction
  ! mz_ht_20110301+
  pBC_ff_insol(1:kproma) = pBC_ant(1:kproma) * pCseason(1:kproma)
  pBC_bb_insol(1:kproma) = pBC_wf(1:kproma) + pBC_ag(1:kproma)
  ! mz_ht_20110301-
  pBC_sum_insol(1:kproma) =  &
       pBC_ant(1:kproma)  * pCseason(1:kproma)  &
       + pBC_wf(1:kproma)                       & !wildfire
       + pBC_ag(1:kproma)                         !agriculture

  ! organic carbon insoluble fraction
  ! mz_ht_20110301+
  pOC_ff_insol(1:kproma) = pOC_ant(1:kproma)   * pCseason(1:kproma) * zom2oc
  pOC_bb_insol(1:kproma) = pOC_wf(1:kproma)  * zom2oc * (1._dp-zbb_wsoc_perc) &
                         + pOC_ag(1:kproma)  * zom2oc * (1._dp-zbb_wsoc_perc)
  pOC_soa_insol(1:kproma) = pOC_bge(1:kproma) * zom2oc * (1._dp-zbge_wsoc_perc)
  ! mz_ht_20110301-
  pOC_sum_insol(1:kproma) = &
       pOC_ant(1:kproma)   * pCseason(1:kproma) * zom2oc     &
       + pOC_wf(1:kproma)  * zom2oc * (1._dp-zbb_wsoc_perc)  &
       + pOC_ag(1:kproma)  * zom2oc * (1._dp-zbb_wsoc_perc)  &
       + pOC_bge(1:kproma) * zom2oc * (1._dp-zbge_wsoc_perc)

  ! organic carbon soluble fraction
  ! mz_ht_20110301+
  pOC_ff_sol(1:kproma) = 0._dp
  pOC_bb_sol(1:kproma) = (pOC_wf(1:kproma) + pOC_ag(1:kproma)) &
                       * zom2oc * zbb_wsoc_perc
  pOC_soa_sol(1:kproma) = pOC_bge(1:kproma) * zom2oc * zbge_wsoc_perc
  ! mz_ht_20110301-
  pOC_sum_sol(1:kproma) = &
       (pOC_wf(1:kproma) + pOC_ag(1:kproma)) * zom2oc * zbb_wsoc_perc  &
                         + pOC_bge(1:kproma) * zom2oc * zbge_wsoc_perc

  pNemis_ks(1:kproma)   = &
       ! biomass burning
       (pOC_wf (1:kproma) + pOC_ag (1:kproma)) * zom2oc &
                      * zbb_wsoc_perc  * zm2n_carb_bb

  pNemis_ki_bc(1:kproma) = &
                    pBC_ant(1:kproma) * pCseason(1:kproma)    * zm2n_carb_ff &
                 + (pBC_wf (1:kproma) + pBC_ag (1:kproma))    * zm2n_carb_bb

  pNemis_ki_oc(1:kproma) = &
                   pOC_ant(1:kproma)                                         &
                           * zom2oc * pCseason(1:kproma)      * zm2n_carb_ff &
                 + (pOC_wf (1:kproma) + pOC_ag (1:kproma))                   &
                           * zom2oc * (1._dp - zbb_wsoc_perc) * zm2n_carb_bb &
                 + pOC_bge(1:kproma)                                         &
                           * zom2oc *(1._dp - zbge_wsoc_perc) * zm2n_carb_bg

  pNemis_ki(1:kproma) = pNemis_ki_bc(1:kproma) + pNemis_ki_oc(1:kproma)

END SUBROUTINE carbon_emissions
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE O3ice_emissions(kproma, pslm, pseaice, pcvs, O3flux)

! original code Laurens Ganzeveld

  INTEGER,  INTENT(IN)  :: kproma
  REAL(dp), INTENT(IN)  :: pseaice(1:kproma) !seaice mask
  REAL(dp), INTENT(IN)  :: pslm(1:kproma)    ! land sea mask
  REAL(dp), INTENT(IN)  :: pcvs(1:kproma)
  REAL(dp), INTENT(OUT) :: O3flux(1:kproma)

  O3flux(:) =0._dp
  O3flux(:) = ((1._dp-pslm(:))*pseaice(:)  &  ! seaice frac.
             + pslm(:)*pcvs(:))*           &  ! snow frac.
             (0.02_dp*1.e-6/48.)*N_A          ! molecules m-2 s-1

END SUBROUTINE O3ice_emissions
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE CH4_emissions(kproma, pdensair, pdz, deltat, cch4,  ch4init, CH4flux)

  INTEGER,  INTENT(IN)  :: kproma
  REAL(dp), INTENT(IN)  :: deltat

  REAL(dp), INTENT(IN)  :: pdensair(1:kproma)
  REAL(dp), INTENT(IN)  :: pdz(1:kproma)
  REAL(dp), INTENT(IN)  :: cch4(1:kproma)
  REAL(dp), INTENT(IN)  :: CH4init(1:kproma)
  REAL(dp), INTENT(OUT) :: CH4flux(1:kproma)

  CH4flux(:) = 0._dp

  CH4flux(:) = (CH4init(:) - cch4(:)) * N_A /(M_air*1.e-3) &
              * pdensair(:) * pdz(:) /deltat

END SUBROUTINE CH4_emissions
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE VOC_emissions(kproma, psrfl, pcossza, plai &
     ,plad, pdm, pemisfac, ptslm1, pemisflux  )

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: kproma
  REAL(dp), INTENT(IN) :: psrfl(1:kproma)
  REAL(dp), INTENT(IN) :: pcossza(1:kproma)
  REAL(dp), INTENT(IN) :: plai(1:kproma)
  REAL(dp), INTENT(IN) :: plad(1:kproma,1:nveglay_hr)
  REAL(dp), INTENT(IN) :: pdm(1:kproma)
  REAL(dp), INTENT(IN) :: ptslm1(1:kproma)
  REAL(dp), INTENT(IN) :: pemisfac(1:kproma,1:3)
  REAL(dp), INTENT(OUT):: pemisflux(1:kproma,1:3)

  ! LOCAL
  REAL(dp)  :: rbvd(1:kproma)              ! direct radiation
  REAL(dp)  :: rvd(1:kproma,1:nveglay_hr)  ! diff. rad. in canopy
  REAL(dp)  :: fsl(1:kproma,1:nveglay_hr)  ! fraction sunlit leaves


    ! ----------------------------------------------------------------------
    !     This program calculates the emission of volatile organic compounds
    !     from vegetation as a function of: biome (Leaf Area Index) and
    !     temperature and Photosynthetically Active Radiation (PAR). The model
    !     considers the extinction of PAR within the canopy as a function
    !     of the Leaf Area Index which is derived from the Olson ecosystems
    !     database (1992), which discerns 72 ecosystems and
    !     their characteristics (see CDROM and paper by Guenther et al., 1995).
    !     This Olson database is also applied to distinguish between different
    !     biomes which show distinct different standard emission factors
    !     (the emission rate taken at a standard temperature and for a standard
    !     amount of PAR, e.g. 30 degrees C and 1000 umol m-2 s-1). The
    !     units of the GEIA database which contains monthly average
    !     emission fluxes is mg C m-2 month-1 and in order to compare
    !     this model with these data the same units are applied where
    !     possible. This model version applies the model of Weiss and Norman
    !     (1985) to calculate the extinction of PAR as a function of the the
    !     Leaf Area Index, the distribution of the LAI (Leaf Area Density),
    !     the fraction of leaves and the orientation of these leaves. This in
    !     contrast with the original model used for the GEIA emission inventory
    !     which applies the formulas by Norman, 1982.
    !
    !     Laurens Ganzeveld 1997, modified October 2001 for implementation in
    !     ECHAM4/5 f90 versions !
    ! ----------------------------------------------------------------------

    ! Interface:
    ! ----------
    ! input
    ! nveglay_hr: number of canopy layers
    ! iisop     : isoprene index, and similar for monoterpenes and other VOC's
    ! dm        : foliar density [g m-2]
    ! lad       : leaf area density profiles [fraction]
    ! voc_emfact: VOC emission factors [ug C g-1 hr-1]
    ! tslm1     : surface temperature [K]
    ! rbvd      : direct incoming radiation [ W m-2]
    ! rvd       : diffusive radiation [W m-2]
    ! fls       : fraction of sunlit leaves [-]
    !
    ! output
    ! voc_emflux: VOC emission flux [molecules m-2 s-1]
    !

  ! LOCAL
  INTEGER ::  ii

  REAL(dp) :: fluxshade(1:kproma,nveglay_hr), fluxsun(1:kproma,nveglay_hr),    &
              clshade(1:kproma,nveglay_hr), foldenslay(1:kproma,nveglay_hr),   &
              pardif(1:kproma,nveglay_hr)
  REAL(dp) :: clsun(1:kproma), pardir(1:kproma), ct(1:kproma)

  REAL(dp), PARAMETER :: xmc=12._dp        ! molar weight of C



    !  -- assigning of values of the used constants, (see Guenther et al.,
    !     1993, JGR). TSC is the leaf temperature at standard conditions,
    !     The term RECALC is the recalculation factor for getting the
    !     net short wave radiation/PAR im umol m-2 s-1 instead of W m-2.
    !     This term is taken as the average of the recalculation factor
    !     for clear sky (4.24) and diffuse conditions (4.57).
    !     See Ecological Physics by J. Hage, D96-6, IMAU and the official
    !     reference is: Grace, J., Plant-Atmosphere relationships, Chapman &
    !     Hall
    !

    REAL(dp), PARAMETER :: alpha=0.0027_dp
    REAL(dp), PARAMETER :: cl1=1.066_dp
    REAL(dp), PARAMETER :: ct1=95000._dp
    REAL(dp), PARAMETER :: ct2=230000._dp
    REAL(dp), PARAMETER :: tsc=303._dp
    REAL(dp), PARAMETER :: tm=314._dp
    REAL(dp), PARAMETER :: r=8.314_dp
    REAL(dp), PARAMETER :: beta=0.09_dp
    REAL(dp), PARAMETER :: recalc=4.405_dp

    INTEGER, PARAMETER :: iisop =1 , imono=2, iovoc =3
    !  --  initialisation of emission flux of sunlit and shaded leaves

    ! INITIALISE:
    rbvd(:)  = 0._dp
    rvd(:,:) = 0._dp
    fsl(:,:) = 0._dp

    CALL calc_profile(kproma, psrfl(1:kproma), pcossza(1:kproma) &
         , plai(1:kproma), plad(1:kproma,1:nveglay_hr) &
         , rbvd(1:kproma), rvd(1:kproma,1:nveglay_hr)  &
         , fsl(1:kproma,1:nveglay_hr))


    fluxsun(:,:)        = 0._dp
    fluxshade(:,:)      = 0._dp
    pemisflux(:,:)       = 0._dp

    !  --  make flux voc_emfact dependent on the radiation (PAR)
    !      and temperature (see Guenther et al., 1993, JGR)
    !      The PAR is calculated from the net short wave radiation
    !      at the surface

    !      In contrast to Guenther et al. 1995 who used the monthly
    !      mean air temperature (Leemans and Cramer), we use the surface
    !      temperature. This surface temperature represents the temperature
    !      of all the four surface cover fractions. Especially for the
    !      semi-arid regions with some vegetation, this can introduce some
    !      bias since the surface temperature will be mainly controlled
    !      by the bare soil temperature. One solution is to introduce some
    !      diagnostically derived leaf temperature from the energy balance
    !      parameters and resistances (21-09-1999)

    !  calculation of temperature attunation function

    ct(:)=exp((ct1*(ptslm1(:)-tsc))/ &
           (r*tsc*ptslm1(:)))/ &
           (1._dp+exp((ct2*(ptslm1(:)-tm))/ &
           (r*tsc*ptslm1(:))))

    !  --  Four layers within the canopy are distinguished
    !      and for each layer the extinction of PAR is determined.
    !      The amount of total biomass is distributed over these canopy
    !      layers, expressed by the Leaf Area Density (LAD) and
    !      combined with the LAI to calculate the emission from each
    !      layer and the total emission from the biome. NLEVV is the
    !      top layer!!! The direct PAR is calculated from the direct visible
    !      radiation and the zenith angle and combined with the fraction
    !      of sunlit leaves and the total biomass yielding the emission flux
    !      of the fraction directly effected by the sun. The diffuse PAR is
    !      a function of the location within the canopy and the fraction of
    !      shaded leaves (1.-FSL)

    pardir(:)=rbvd(:)*recalc    ! direct incoming PAR

    DO ii=1,nveglay_hr          ! loop vertical layers (radiation profile)
       pardif(:,ii)=rvd(:,ii)*recalc        ! diffuse PAR
       foldenslay(:,ii)=pdm(:)*plad(:,ii)   ! foliar density profile
       ! light attenuation function sunlit leaves
       clsun(:)=(alpha*cl1*pardir(:))/ (sqrt(1._dp+alpha**2*pardir(:)**2))
       ! light att. funct. shaded leaves
       clshade(:,ii)=(alpha*cl1*pardif(:,ii))/ &
            (sqrt(1._dp+alpha**2*pardif(:,ii)**2))
       ! flux from sunlit leaves
       fluxsun(:,ii)=pemisfac(:,iisop)*clsun(:)*ct* &
            foldenslay(:,ii)*fsl(:,ii)
       ! flux from shaded leaves
       fluxshade(:,ii)=pemisfac(:,iisop)*clshade(:,ii)*ct* &
            foldenslay(:,ii)*(1._dp-fsl(:,ii))

       ! calculation of integrated isoprene emission rate for
       ! bulk approach, the emission is in ug C m-2 hr-1 and is
       ! recalculated to the emission in [kg C m-2 s-1] (1.E9/3600)
       ! and from that to molecules m-2 s-1. The term 1.E3 it
       ! to recalculate from kg to g, 1/XMC to recalculate to mol C
       ! and the term 1/5 is to correct for the 5 C molecules.
       ! In order to recalc from mol isoprene m-2 s-1 to molecules its
       ! multiplied with the avogadro number

       ! total flux for use in "big leaf" model
       pemisflux(:,iisop)=pemisflux(:,iisop)+         &
            (fluxsun(:,ii)+fluxshade(:,ii))*1.e-9/(3600.)*      &
            1.e3*(1./xmc)*(1._dp/5.)*N_A

    ENDDO ! end loop vertical layers

    ! emission of monoterpenes and OVOC's

    ! including the recalculation to the units  molecules m-2 s-1, with for
    !  monoterpenes assuming 10 C and for  OVOC's 15 C per molecule

    pemisflux(:,imono)=pemisfac(:,imono)*exp(beta*(ptslm1(:)-tsc))* &
         pdm(:)*1.e-9/(3600.)*1.e3*(1./xmc)*(1._dp/10.)*N_A
    pemisflux(:,iovoc)=pemisfac(:,iovoc)*exp(beta*(ptslm1(:)-tsc))* &
         pdm(:)*1.e-9/(3600.)*1.e3*(1./xmc)*(1._dp/15.)*N_A

CONTAINS

  SUBROUTINE calc_profile(kproma,ppsrfl, ppcossza_2d, pplai &
       , pplad, prbvd, prvd, pfsl)
    ! ---------------------------------------------------------------
    !     Calculation of distribution of PAR (diffuse and direct)
    !     within the canopy as a function of the solar zenith angle
    !     the LAI and the radiation above the canopy. The code
    !     is taken from the DDIM model (Dry deposition Inferential
    !     Model) and developed by Norman and Weiss, 1985. The
    !     windprofile is also calculated acc. to Cionco.
    !
    !     Laurens Ganzeveld, 1998, modified for implementation
    !     in echam5, October, 2001
    ! --------------------------------------------------------------------

    ! Interface:
    ! ----------
    ! input
    ! srfl      : net surface radiation [W M-2]
    ! cossza_2d : cosine of the zenith angle [0-1]
    ! lai       : Leaf Area Index [m2 m-2]
    ! lad       : leaf area density profile
    !
    ! output
    ! rbvd      : direct beam irradiance [W m-2]
    ! rvd       : diffusive irradiance in canopy [W m-2]
    ! fsl       : fraction of sunlit leaves [0-1]

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)   :: kproma
    REAL(dp), INTENT(in)  :: ppsrfl(1:kproma)
    REAL(dp), INTENT(in)  :: ppcossza_2d(1:kproma)
    REAL(dp), INTENT(in)  :: pplai(1:kproma)
    REAL(dp), INTENT(in)  :: pplad(1:kproma,1:nveglay_hr)
    REAL(dp), INTENT(out) :: prbvd(1:kproma)
    REAL(dp), INTENT(out) :: pfsl(1:kproma,1:nveglay_hr)
    REAL(dp), INTENT(out) :: prvd(1:kproma,1:nveglay_hr)

    INTEGER :: ii, iday, jl

    REAL(dp) :: zrvd(nveglay_hr),  zfsl(nveglay_hr), &
                mlai(nveglay_hr,2),mlaitot(nveglay_hr)

    REAL(dp) :: rg,zen,parbeam,zrbvd,twopi

    !
    !  *********************************************************************
    !  *                                                                   *
    !  *        K E Y         V A R I A B L E S                            *
    !  *                                                                   *
    !  *                                                                   *
    !  *   RG = GLOBAL RADIATION (WATTS/M^2)                               *
    !  *                                                                   *
    !  *   LAI = LEAF AREA INDEX                                           *
    !  *                                                                   *
    !  *   JDAY = JULIAN DAY                    LAIW = WINTER LAI          *
    !  *                                                                   *
    !  *   H! = CANOPY HEIGHT (M)                                          *
    !  *                                                                   *
    !  *   FSL = FRACTION OF SUNLIT LEAVES    CSNL = COSINE OF LEAF NORMAL *
    !  *                                                                   *
    !  *   RVD = DIFFUSE VISIBLE RADIATION    RBVD = VISIBLE BEAM RADIATION*
    !  *                                                                   *
    !  *   PAR = PHOTOSYNTHETICALLY ACTIVE RADIATION (.4 - .7 MICRONS)     *
    !  *                                                                   *
    !  *   FVIS = FRACTION OF THE GLOBAL RADIATION THAT IS PAR             *
    !  *                                                                   *
    !  *   PCNTLF = PERCENTAGE OF MAX LEAF AREA FOR NON-CONIFERS           *
    !  *                                                                   *
    !  *                                                                   *
    !  *   ** UNLESS SPECIFIED, UNITS ARE SI                               *
    !  *                                                                   *
    !  *********************************************************************
    !

    twopi = 2.*pi

    kproma_loop: DO jl=1,kproma
       rg=ppsrfl(jl)
       zen=ACOS(ppcossza_2d(jl))

       !  Be carefull with LAD vertical profile. LAD(1) is the LAD value for
       !  the top vegetation layer

       ! normally the LAD profile for one ecosystem type (mlai(:,1) is
       ! considered but this can be extended to more vegetation types

       DO ii=1,nveglay_hr
          mlai(ii,1) = pplad(jl,ii)*pplai(jl)
          mlaitot(ii) = mlai(ii,1)
       ENDDO

       iday = 1

       IF((rg.LT.10._dp).OR.(zen.GT.(twopi/4.))) then
          iday = 0
          zrbvd = 0._dp
          DO ii=1,nveglay_hr
             zrvd(ii) = 0._dp
             zfsl(ii) = 0._dp
          ENDDO
       ENDIF

       ! call of routine in which the radiation profiles in the
       ! canopy are being calculated

       IF(iday.GT.0)  &
            CALL canopy_radiation(nveglay_hr,mlaitot,zen,rg, &
                                zfsl,zrvd,zrbvd,parbeam)
       DO ii=1,nveglay_hr
          pfsl(jl,ii)= zfsl(ii)
          prvd(jl,ii)= zrvd(ii)
       ENDDO

       prbvd(jl)=zrbvd

    ENDDO kproma_loop

  END SUBROUTINE calc_profile

  !=============================================================

  SUBROUTINE canopy_radiation(n,pai,zen,rg,fsl,rvd,rbvd, parbeam)
    !
    !    ***************************************************
    !    *  S U B R O U T I N E    C A N R A D 2           *
    !    *                                                 *
    !    *  THIS SUBROUTINE COMPUTES THE VERTICAL PROFILE  *
    !    *  OF VISIBLE RADIATION (PAR), BOTH BEAM AND      *
    !    *  DIFFUSE COMPONENTS FOR A MULTILAYER CANOPY     *
    !    *  NLEV IS THE TOPLAYER ABOVE THE CANOPY!!!     *
    !    ***************************************************

    !  __________________________________________________________________
    !
    !    It turned out that there is an error in this code
    !    RBVD is calculated from the RDV by dividing through cos(ZEN),
    !    which is incorrect!!!!. It likely has to do with the definition
    !    of the zenith angle in DDIM
    !  _________________________________________________________________
    !

    IMPLICIT NONE

    INTEGER :: N

    REAL(dp) :: PAI(N),RG,ZEN,FVD,RBV(N),PI, &
         ALPHA,KXM,IB(N),MIB(N),AV(N), &
         TV,PV,X1(11),KXD(11),ID(N),MID(N),Z1,Z2,Z3,FSL(N), &
         RVD(N),RVU(N),ORVU,CKVU,CPAI

    INTEGER :: I,J,K,M,ITER

    REAL(dp) :: OT,RDVIS,RFVIS,WA,RDIR,RFIR,RVT,RIRT,FVIS, &
         FIR,RATIO,FVB,PARBEAM,THETA,RBVD

    !
    !    ******************************************************************
    !    *                                                                *
    !    *   THE NEXT FEW STATEMENTS DETERMINE THE BEAM AND DIFFUSE       *
    !    *   COMPONENTS OF VISIBLE RADIATION (FOR DETAILS SEE WEISS AND   *
    !    *   NORMAN, 1985, AGRICULTURAL METEOROLOGY.......                *
    !    *                                                                *
    !    ******************************************************************
    !

    ZEN=MIN(ZEN,1.56_dp)

    OT=35./(1224.*COS(ZEN)**2.+1._dp)**.5
    RDVIS=600.*2.7182**(-.185*OT)*COS(ZEN)
    RFVIS=0.4*(600._dp-RDVIS)*COS(ZEN)
    WA= 1320._dp*.077*(2.*OT)**0.3
    RDIR=(720.*2.7182**(-0.06*OT)-WA)*COS(ZEN)
    RFIR=0.60*(720._dp-WA-RDIR)*COS(ZEN)
    RVT=RDVIS+RFVIS
    RIRT=RDIR+RFIR
    FVIS=RVT/(RIRT+RVT)
    FIR=RIRT/(RIRT+RVT)
    RATIO=RG/(RVT+RIRT)
    IF(RATIO.GE.0.9_dp) THEN
       RATIO=0.899_dp
       RG=RVT+RIRT
    ENDIF

    FVB=RDVIS/RVT*(1._dp-((0.9_dp-RATIO)/7.)**0.67)
    FVD=1._dp-FVB

    PARBEAM = MAX(1.E-5_dp,RG*FVIS*FVB*COS(ZEN))

    !
    !     **************************
    !     *  INITIALIZE CONSTANTS  *
    !     **************************
    !

    PI=3.1415926_dp

    !
    !     *******************************************************
    !     *  SET FRACTION OF TOTAL THAT IS BEAM RADIATION AND   *
    !     *  THEN SEPARATE INTO INTO IR AND VISIBLE COMPONENTS  *
    !     *******************************************************
    !
    !

    RBV(N)=RG*FVIS*(1._dp-FVD)

    !
    !     *****************************************
    !     * COMPUTATION OF EXTINCTION COEFFICIENT *
    !     *****************************************
    !

    ALPHA=PI/2. - ZEN
    KXM = 0.5/SIN(ALPHA)
    THETA = PI/36.
    DO I=1,9
       THETA = THETA+PI/18.
    ENDDO

    !
    !     ****************************************************
    !     *  COMPUTE PROBABILTIY FUNCTION FOR PENETRATION    *
    !     *  OF THE BEAM COMPONENT , FROM NORMAN, 1979       *
    !     *  MODIFICATION OF THE AERIAL ENVIRONMENT OF CROPS *
    !     ****************************************************
    !

    DO I=1,N-1
       K=N - I

       ! originally it is PAI(K+1) but this has been replaced by  PAI(I+1)

       IB(K+1)=2.7182818**(-KXM*PAI(I+1))
       MIB(K+1)=1.0_dp - IB(K+1)
       RBV(K)=RBV(K+1)*IB(K+1)
    ENDDO

    !
    !     ************************************
    !     *  SET SOIL VISIBLE AND IR ALBEDO  *
    !     ************************************
    !

    AV(1)=0.10_dp
    TV=0.01_dp
    PV=0.08_dp
    ALPHA=PI/20.

    !
    !     **************************************************
    !     *   LOOP FOR COMPUTING OFTENLY USED FACTORS      *
    !     **************************************************
    !

    DO I=2,11
       X1(I)=SIN(ALPHA)*COS(ALPHA)
       KXD(I)=.5/SIN(ALPHA)
       ALPHA=ALPHA+PI/20.
    ENDDO

    !
    !     *******************************************************
    !     * LOOP FOR COMPUTING ID (DIFFUSE RADIATION PENETRATION*
    !     * FUNCTION) AND A (R UP/R DOWN)                       *
    !     *******************************************************
    !

    DO I=2,N
       ID(I)=0.0_dp
       Z1=0.0_dp
       K=1
       DO J=1,5

          ! LG-          originally it is PAI(I) but this has been replaced by
          !              PAI(N+1-I)

          Z2=2.718282**(-PAI(N+1-I)*KXD(K+1))*X1(K+1)

          ! LG-          originally it is PAI(I) but this has been replaced by
          !              PAI(N+1-I)

          Z3=2.718282**(-PAI(N+1-I)*KXD(K+2))*X1(K+2)
          ID(I)=ID(I)+PI*(Z1+4.*Z2+Z3)/(20._dp*3.)
          Z1=Z3
          K=K+2
       ENDDO
       ID(I)=ID(I)*2.0
       IF(ID(I).GT.1._dp) ID(I) = 1._dp
       MID(I)=1.0_dp - ID(I)
       AV(I)=AV(I-1)*(TV*MID(I)+ID(I))*(TV*MID(I)+ID(I))/  &
            (1.0_dp - AV(I-1)*PV*MID(I)) + PV*MID(I)
    ENDDO

    !
    !     *************************************************
    !     *  INITIALIZE DOWNWARD DIFFUSE COMPONENTS OF    *
    !     *  VISIBLE RADIATION THE CANOPY                 *
    !     *************************************************
    !

    RVD(N)=RG*FVD*(1.0_dp - FIR)

    !
    !    ****************************************************
    !    * COMPUTE DOWNWARD DIFFUSE RADIATION AT EACH LEVEL *
    !    ****************************************************
    !
    DO I=1,N-1
       K=N-I
       RVD(K)=RVD(K+1)*(TV*MID(K+1)+ID(K+1))/(1.0_dp-AV(K)*PV &
            *MID(K+1))
    ENDDO

    !
    !     *************************************************
    !     *  COMPUTE UPWARD DIFFUSE FLUXES USING THE SOIL *
    !     *  ALBEDO                                     *
    !     *************************************************
    !

    RVU(1)=AV(1)*(RVD(1)+RBV(1))

    DO I=2,N
       RVU(I)=RVU(I-1)*(TV*MID(I)+ID(I))*AV(I)/(AV(I)-PV*MID(I))
    ENDDO

    !
    !     *************************************
    !     *  START MAIN LOOP FOR ITERATIONS   *
    !     *************************************
    !
    !

    ORVU=0.0_dp
    ITER=0
    DO M=1,100
       ITER=ITER + 1
       DO I=1,N-1
          K = N - I
          RVD(K)=RVD(K+1)*(TV*MID(K+1)+ID(K+1))+RVU(K)*PV*MID(K+1) &
               +RBV(K+1)*MIB(K+1)*TV
       ENDDO

       RVU(1)=AV(1)*(RVD(1)+RBV(1))

       !
       !     ************************************************
       !     *   COMPUTE THE UPWARD DIFFUSE FLUXES FOR      *
       !     *   THE VISIBLE WAVELENGTHS                    *
       !     ************************************************
       !

       DO I=2,N
          RVU(I)=RVU(I-1)*(TV*MID(I)+ID(I))+RVD(I)*PV*MID(I)+ &
               RBV(I)*MIB(I)*PV
       ENDDO

       CKVU=ABS(RVU(1)-ORVU)

       !
       !     *****************************************************
       !     *  CHECK FOR CONVERGENCE OF VALUES                  *
       !     *  CRITERIA FOR CONVERGENCE: DIFFERENCE BETWEEN     *
       !     *  THE OLD AND NEW VALUES AT THE TENTH/FOURTH LEVEL *
       !     *  CHANGES BY NO MORE THAN 2 WATTS/METER SQ.        *
       !     *****************************************************
       !
       !

       IF((CKVU.LE.0.01_dp)) GO TO 100
       ORVU=RVU(1)
    ENDDO

100 CPAI = 0.0_dp

    DO I=1,N
       K = N+1-I

       ! in original code PAI(K) with top layer=N, replaced by PAI(I)

       CPAI=CPAI+PAI(I)
       FSL(K)=2.7183**(-KXM*CPAI)
    ENDDO
    RBVD=RBV(N)*COS(ZEN)

    RETURN

  END SUBROUTINE canopy_radiation

! ============================================================================


END SUBROUTINE VOC_emissions
! -------------------------------------------------------------------------

! mz_sg_20091125+
! -------------------------------------------------------------------------
SUBROUTINE terr13C_frac(kproma, delta13C, qatom, flux, flux12C, flux13C)

  ! --------------------------------------------------------------------
  ! emission flux "fracturizer" according to the delta-13C given
  !
  ! Sergey Gromov, 2009
  ! --------------------------------------------------------------------

  ! Interface:
  ! ----------
  ! input
  ! delta13c   : flux isotopic composition as delta-13C (permil VPDB 13C)
  ! qatom      : quantity of carbon atoms in the species molecule
  ! flux       : total flux
  !
  ! output
  ! frac12/13C : fractured fluxes

    IMPLICIT NONE

  ! I/O
    INTEGER,  INTENT(IN)  :: kproma
    REAL(dp), INTENT(IN)  :: delta13C(1:kproma)
    INTEGER,  INTENT(IN)  :: qatom
    REAL(dp), INTENT(IN)  :: flux(1:kproma)
    REAL(dp), INTENT(OUT) :: flux12C(1:kproma), flux13C(1:kproma)

    REAL(dp)              :: r13C(1:kproma), f12C(1:kproma), f13C(1:kproma)

  ! Vienna PeeDee Belemnite (VPDB) 13C carbon reference ratio
    REAL(dp), PARAMETER   ::  R_VPDB_13C = 1.12372E-2_dp

  ! isotope ratio
    r13C(:) = R_VPDB_13C * ( delta13C(:) / 1000.0_dp + 1.0_dp )

  ! isotopologues fraction
    f13C(:) = ( r13C(:) * qatom ) / ( r13C(:) + 1.0_dp )
    f12C(:) = 1.0_dp - f13C(:)
  ! f12C(:) = ( r13C(:) * ( 1.0_dp - qatom ) + 1.0_dp ) / ( r13C(:) + 1.0 )

  ! fracturing the flux
    flux12C(:) = flux(:) * f12C(:)
    flux13C(:) = flux(:) * f13C(:)

END SUBROUTINE terr13C_frac
! -------------------------------------------------------------------------
! mz_sg_20091125-

! -------------------------------------------------------------------------
SUBROUTINE  NOpls_emissions(kproma, ndaylen, delta_time, init_step, nstep &
     , prc , prl, cpold,  lspold , pulsing, plsday, plsdurat, cp, prectot &
     , lsp, pls)

    ! ------------------------------------------------------------------
    !     This program adds the convective and large scale precipitation
    !     and calculates the pulsing of the NO emission occuring after
    !     a rainfall event after a period of drought. For the pulsing,
    !     a history must be recorded and therefore the precipitation data
    !     of the previous month are also incorporated in the calculation
    !     of the NO emission pulse.
    !
    !     edited by Laurens Ganzeveld, 1998, modified for implementation
    !     in echam4/echam5 f90, Laurens Ganzeveld, October, 2001
    ! ------------------------------------------------------------------

    ! mz_lg_20040426+
    ! Interface:
    ! ----------
    ! input
    ! klon      : number of longitudes
    ! klat      : number of latitudes
    ! jrow      : latitude index
    ! nstep     : timestep
    ! init_step : initial timestep
    ! ndaylen   : length of day in seconds
    ! delta_time: timestep
    ! ndrydays  : number of dry days needed to get a pulse (default 14)
    ! prc       : convective rainfall [m]
    ! prl       : large-scale rainfall [m]
    !
    ! input/output
    ! prectot   : accumulated rainfall of last month [mm]
    ! cpold     : convective rainfall record for ndrydays
    ! lspold    : large-scale rainfall record for ndrydays
    ! pulsing   : pulsing regime    [index, 1-3]
    ! plsday    : timing of pulse   [number of timesteps that pulse is active]
    ! plsdurat  : duration of pulse [days]
    !
    ! output
    ! cp        : daily accumulated convective rainfall
    ! lsp       : daily accumulated large-scale rainfall
    ! cpold     : convective rainfall record previous 14 days
    ! lspold    : large-scale rainfall record previous 14 days
    ! pls       : the actual pulse  [ - ]

    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: kproma
    INTEGER, INTENT(IN)     :: ndaylen
    INTEGER, INTENT(IN)     :: init_step, nstep
    REAL(dp), INTENT(IN)    :: delta_time
    REAL(dp), INTENT(in)    :: prc(1:kproma), prl(1:kproma)
    REAL(dp), INTENT(inout) :: cpold(1:kproma,1:ndrydays)          &
                             , lspold(1:kproma,1:ndrydays)         &
                             , pulsing(1:kproma), plsday(1:kproma) &
                             , plsdurat(1:kproma)
    REAL(dp), INTENT(INOUT) :: prectot(1:kproma)
    REAL(dp), INTENT(inout)   :: cp(1:kproma),  lsp(1:kproma), pls(1:kproma)

    ! LOCAL
    REAL(dp), PARAMETER :: wetdry=10._dp
    INTEGER             :: nstepday
    INTEGER             :: istep
    INTEGER             :: k, jl
    REAL(dp)            :: rainsum(1:kproma)
    REAL(dp)            :: prect(1:kproma)
    !--- 0) Initialisations: ----------------------------------------------

    ! determining the number of timesteps of one day
    nstepday= INT(REAL(ndaylen,dp)/delta_time)

    istep=nstep-init_step+1

    ! start calculation of pulsing as a function of total precipitation,
    ! the summed precipitation during a period of 14 days
    !
    ! CP AND LSP are the daily cumulative precipitation,
    ! so the actual precipitation rates are added

    cp(:) = cp(:) + prc(:)
    lsp(:)= lsp(:)+ prl(:)

    IF (MOD(istep,nstepday).eq.0) THEN
       DO k=1,ndrydays-1
          cpold(:,k)=cpold(:,k+1)
          lspold(:,k)=lspold(:,k+1)
       ENDDO
       cpold(:,ndrydays)=cp(:)
       lspold(:,ndrydays)=lsp(:)
       cp(:)=0._dp
       lsp(:)=0._dp
    ENDIF

    ! determining the total precipitation of the previous ndrydays
    rainsum(:)=0.
    DO k=1,ndrydays
       rainsum(:)=rainsum(:)+(cpold(:,k)+lspold(:,k))*1000.
    ENDDO

    !   -- determining the total precipitation of the previous
    !      30 days used as proxy to determine if we are in the wet or
    !      dry season in the tropics.
    prectot(:)=rainsum(:)*(30./ndrydays)

    ! calculating the total amount of precipitation in  mm day-1
    ! (see Yienger and Levy to determine the pulsing regimes).
    ! The daily cumulative rainfall is considered in mm

    prect(:)=(cp(:)+lsp(:))*1000.

    ! the parameter n is the no. of pulsing in each specific class
    ! 1 < prec < 5, 5 < prec < 15 and prec > 15 (n5, n10 and
    ! n15 respect.) Three pulsing regimes are discerned as a function
    ! of the intensity of the precipitation (pulsing=1,2,or 3). The
    ! pulse last for 3,7 or 14 days for the three regimes (plsdurat)
    ! starting at day=1 (all units in mm !!! )

    DO jl=1,kproma
       IF (rainsum(jl).LT.wetdry .AND. INT(plsday(jl)).EQ.0) THEN
          IF (prect(jl) .GT. 1._dp .AND. prect(jl) .LE. 5._dp) THEN
             pulsing(jl)=1._dp
             ! mz_js_20081210+
             ! values of plsday below 1. cause a much to high pulse
             ! so plsday is set to one here (all pulsing regimes).
             ! The pulse lasts now for a slightly shorter period,
             ! but this is better than a much to high pulse.
             plsday(jl)=1.
             ! mz_js_20081210-
             plsdurat(jl)=3._dp
          ELSE IF (prect(jl) .GT. 5._dp .AND. prect(jl) .LE. 15._dp) THEN
             pulsing(jl)=2._dp
             plsday(jl)=1.
             plsdurat(jl)=7._dp
             ! mz_js_20081210 moved the "heavy rain" in here. Was outside the
             !                check of rainsum and plsday
          ELSE IF (prect(jl) .GT. 15._dp) THEN
             pulsing(jl)=3._dp
             plsday(jl)=1.
             plsdurat(jl)=14._dp
          ENDIF
       ENDIF

       IF (INT(pulsing(jl)) .EQ. 1 .AND. &
            plsday(jl) .LE. plsdurat(jl)) THEN
          pls(jl)=MAX(1._dp,11.19*EXP(-0.805*plsday(jl)))
       ELSE IF (INT(pulsing(jl)) .EQ.2 .AND. &
            plsday(jl) .LE. plsdurat(jl)) THEN
          pls(jl)=MAX(1._dp,14.68*EXP(-0.384*plsday(jl)))
       ELSE IF (INT(pulsing(jl)) .EQ. 3 .AND. &
            plsday(jl) .LE. plsdurat(jl)) THEN
          pls(jl)=MAX(1._dp,18.46*EXP(-0.208*plsday(jl)))
       ELSE
          pls(jl)=1._dp
       ENDIF

       ! increasing of the day, so this needs to be corrected
       ! for the number of timesteps of each day

       IF (plsday(jl) .GT. 0._dp) plsday(jl)=plsday(jl)+1./REAL(nstepday,dp)
       IF (plsday(jl) .GT. plsdurat(jl)) THEN
          pulsing(jl)=0._dp
          plsday(jl)=0._dp
          plsdurat(jl)=0._dp
       ENDIF

    ENDDO

  END SUBROUTINE NOpls_emissions
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE  NO_emissions(kproma,  month  &
       , cultiv,  fertil, tsoil, ws , prectot &!, no_emflux_mm &
       , NOemis_w, NOemis_d, noemclass, lai, latitude, NOflux  &
       , noslflux, noslflux_diag)      ! mz_js_20081021
    ! ---------------------------------------------------------------------
    !     This program calculates the soil-biogenic NO-emission
    !     as a function of: biome, soil wetness, soil temperature,
    !     the distribution of cultivation/agriculture, N-fertilizer loss,
    !     canopy reduction, the pulsing (which is the enhanced emission
    !     due to rainfall), and rice-emission-reduction. The program
    !     is originally developed by Peter van den Broek, 1995, and
    !     edited by Laurens Ganzeveld 1996/1998.
    !
    !     10-2001, Modified for including the code in echam4/echam5 f90. More
    !     information about the model can be found in the paper by Yienger
    !     and Levy, "Empirical model of global soil-biogenic NOx emissions"
    !     JGR 100, 1995 and Ganzeveld et al., "The influence of soil-biogenic
    !     NOx emissions on the global distribution of reactive trace gases:
    !     the role of canopy processes", submitted to JGR, 2001. There is
    !     also the option to use an alternative emission inventory by
    !     Davidson, E., and W. Kingerlee, "A global inventory of nitric oxide
    !     emissions from soils", Nutrient Cycling in Agroecosystems, 48, 37-50,
    !     1997. The two different inventories are referred to in this routine
    !     by YL95 and DK97
    !
    !     Laurens Ganzeveld, October, 2001
    ! ---------------------------------------------------------------------

    ! Interface:
    ! ----------
    ! input
    ! klat      : number of latitudes
    ! trop_class : index number that resembles the tropical forest emis. class
    ! cultiv    : cultivation intensity [old: 0-15, new 2004: 0-1]
    ! fertil    : fertilizer application (synthetic/manure)
    ! tsoil     : soil temperature  [K]
    ! ws        : soil moisture     [m]
    ! prc       : convective rainfall [m]
    ! prl       : large-scale rainfall [m]
    ! prectot   : monthly accumulated precipitation [mm]
    ! no_emflux_mm: monthly mean emission flux [molec m-2 s-1]
    ! noemis_w  : wet soil emission factor [ng N m-2 s-1]
    ! noemis_d  : dry soil emission factor [ng N m-2 s-1]
    ! noemclass : NO emission class [0-12, YL95, 0-17 DK97]
    ! lai       : leaf area index [m2 m-2]
    ! lcrfyl95  : switch to use YL95's Canopy Reduction Factor (CRF)
    ! l_veg_mlay: switch to indicate use of explicit canopy model
    ! latitude  : latitude
    !
    ! output
    ! no_emflux : NO soil-biogenic emission flux from canopy [molecules m-2 s-1]
    ! noslflux  : NO soil-biogenic emission flux from soil   [ng m-2 s-1]
    ! noslflux_diag : diagnostic fluxes of indiv. ecosystems [ng m-2 s-1]

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(IN)    :: kproma
    INTEGER,  INTENT(IN)    :: month
    REAL(dp), INTENT(in)    :: cultiv(1:kproma),  fertil(1:kproma)
    REAL(dp), INTENT(in)    :: tsoil(1:kproma),  ws(1:kproma)
    REAL(dp), INTENT(in)    :: prectot(1:kproma)
    REAL(dp), INTENT(in)    :: latitude(1:kproma)
    REAL(dp), INTENT(in)    :: NOemis_w(1:kproma),NOemis_d(1:kproma)
    REAL(dp), INTENT(in)    :: noemclass(1:kproma,ncl_yl95)
    REAL(dp), INTENT(in)    :: lai(1:kproma)
    REAL(dp), INTENT(out)   :: noflux(1:kproma)

    REAL(dp), INTENT(out)   :: noslflux(1:kproma)
    REAL(dp), INTENT(out)   :: noslflux_diag(1:kproma, ncl_yl95)


    ! LOCAL
    REAL, PARAMETER, DIMENSION(ncl_yl95) :: sai_veg = &
         (/0.,0.,0.,0.010,0.019,0.,0.030,0.025,0.036,0.075,0.120,0.032/)
    !switch for YL95's Canopy Reduction Factor (CRF)
    LOGICAL :: lcrfyl95 = .true.
    ! switch to indicate use of explicit canopy model
    LOGICAL :: l_veg_mlay = .false.

    REAL(dp) :: fert,    fertseas,soiltemp, soilws,                   &
         prect,   fwd,      fwdagri,   cult,                 &
         frc_trp, crf,      ks,        kc,      sai

    ! threshold value to distinguish wet and dry soils
    REAL(dp), PARAMETER :: wtd=0.10_dp      ! and dry soils
    ! modified fraction of NO loss  based on paper by Lex Bouwman, GBC: 2002
    ! "Modeling of N2O and NO emissions...": was 0.025)
    REAL(dp), PARAMETER :: tmp1=10._dp         ! threshold temperature
    REAL(dp), PARAMETER :: tmp2=30._dp         ! threshold temperature
    REAL(dp), PARAMETER :: fctr=0.103_dp

    REAL(dp), PARAMETER :: mtsh = -30._dp
    REAL(dp), PARAMETER :: mtnh = 30._dp
    INTEGER :: ic, jl, vegtype

    REAL(dp),PARAMETER :: amn   = 14.00_dp       ! molecular weight of N



    ! --  Definition of threshold latitude to distinguish midlatitudes
    !     from the tropics, mtsh and mtnh (MidlatitudesTropicsSH and NH)
    !     the border is in this study at 30 N and 30 S

    kproma_loop: DO jl = 1, kproma

       ! assigning the field to local parameters

       cult=cultiv(jl)
       fert=fertil(jl)
       soiltemp=tsoil(jl)-273.15_dp ! in echam4, td3 was being used
       soilws=ws(jl)

       prect=prectot(jl)

       !  --   initialisation of flux

       fwd=0.0_dp
       fwdagri=0.0_dp

       ic=INT(cult)

       IF (ic.EQ.1.OR.ic.EQ.-9999) cult=0._dp

       ! --    Bouman data

       IF (ic.EQ.2.OR.ic.EQ.12) cult=0.20_dp
       IF (ic.EQ.3.OR.ic.EQ.13) cult=0.50_dp
       IF (ic.EQ.4.OR.ic.EQ.14) cult=0.75_dp
       IF (ic.EQ.5.OR.ic.EQ.15) cult=1._dp

       ! fill flux (fwd) with biome-related emission-factor
       ! also dependent on soil wetness (A(w/d)), non-agriculture

       IF (soilws.GT.wtd) THEN              !wet
          fwd=noemis_w(jl)
       ELSE                                 !dry
          fwd=noemis_d(jl)
       ENDIF

       ! Not performing futher corrections of emission factors ofthe DK97
       ! inventory since initial tests show that this yields a much too large
       ! emission flux (15-10-2000)


          ! make flux (fwd) dependent of temperature, soil wetness
          ! This is not done for rain forest (f(w/d))

          ! determining the fraction of coverage with rain forest
          ! and from that the surface cover of the other ecosystems

          frc_trp=noemclass(jl,trop_class)

          ! only correction for the YL95 inventory, however still assiging
          ! the temperature correction function for interpretation

          IF (soilws.GT.wtd) then   ! wet
             IF (soiltemp.LT.0.) THEN
                fwd=0.
                noslflux_diag(jl,1:trop_class-1) = 0.
             ELSE IF (soiltemp.ge.0._dp.and.soiltemp.le.tmp1) THEN
                fwd=0.28*fwd*soiltemp
                noslflux_diag(jl,1:trop_class-1) = &
                     0.28*noemfact_wet(1:trop_class-1)*soiltemp
             ELSE IF (soiltemp.GT.tmp1.and.soiltemp.le.tmp2) THEN
                fwd=fwd*exp(fctr*soiltemp)
                noslflux_diag(jl,1:trop_class-1) = &
                     noemfact_wet(1:trop_class-1)*exp(fctr*soiltemp)
             ELSE IF (soiltemp.GT.tmp2) THEN
                fwd=21.97*fwd
                noslflux_diag(jl,1:trop_class-1) = &
                     21.97*noemfact_wet(1:trop_class-1)
             ENDIF
          ELSE                         ! dry
             IF (soiltemp.LT.0.) THEN
                fwd=0.
                noslflux_diag(jl,1:trop_class-1) = 0.
             ELSE IF (soiltemp.ge.0._dp.and.soiltemp.le.tmp2) THEN
                fwd=fwd*soiltemp/tmp2
                noslflux_diag(jl,1:trop_class-1) = &
                     noemfact_dry(1:trop_class-1)*soiltemp/tmp2
             ELSE IF (soiltemp.GT.tmp2) THEN
                noslflux_diag(jl,1:trop_class-1) = &
                     noemfact_dry(1:trop_class-1)
             ENDIF
          ENDIF

          ! fill flux (fwd) for rain forest, depending on month  and biome
          ! (f(w/d) rain forest), non-agriculture distinguishing between wet
          ! and dry season based on the total amount
          ! of precipitation of 150 mm/month for any grid square in the tropics

          IF (prect.LT.150._dp) THEN ! 5 driest months
             fwd=frc_trp*noemfact_dry(trop_class)+fwd
             noslflux_diag(jl,trop_class) = noemfact_dry(trop_class)
          ELSE
             fwd=frc_trp*noemfact_wet(trop_class)+fwd
             noslflux_diag(jl,trop_class) = noemfact_wet(trop_class)
          ENDIF
          noslflux(jl) = fwd

          ! explicit calculation of canopy reduction factor acc. to Yienger
          ! and Levy
          crf=1._dp
          ks=8.75_dp
          kc=0.24_dp
          sai=0._dp

          ! determining the grid average Stomatal Area Index (SAI) from the
          ! individual SAI values of the twelve ecosystems
          DO vegtype=1,ncl_yl95-1
             sai=sai+REAL(sai_veg(vegtype),dp)*noemclass(jl,vegtype)
          ENDDO

          crf=(exp(-ks*sai)+exp(-kc*lai(jl)))/2.

          ! Yienger and Levy canopy reduction factor, only applying this
          ! term when the bigleaf approach is being used and the switch
          ! LCRFYL95 is set to TRUE
          IF (lcrfyl95) THEN
             IF (.not.l_veg_mlay) fwd=fwd*crf
          ENDIF

          !  ---------------------------------------------------------------
          ! start of calculations for agricultural areas
          !
          !  --   fill flux (fwdagri with biome-factor and fert,
          !       dependent of soilwetness and month (A(w/d)),agriculture
          IF (cult.ge.0.20_dp) THEN
             fwdagri=0.36_dp
          ENDIF

          !  --   makes flux (fwdagri dependent on temperature, soil wetness,
          !       agriculture (f(w/d))
         IF (soiltemp.LT.0.) THEN
             fwdagri=0.
          ELSE IF (soiltemp.ge.0..and.soiltemp.le.tmp1) THEN
             fwdagri=0.28*fwdagri*soiltemp
          ELSE IF (soiltemp.GT.tmp1.and.soiltemp.le.tmp2) THEN
             fwdagri=fwdagri*exp(fctr*soiltemp)
          ELSE IF (soiltemp.GT.tmp2) THEN
             fwdagri=21.97*fwdagri
          ENDIF

          !  --   add effect of fertilizer in the periods of application

          fertseas=0._dp

          IF (latitude(jl).GE.mtnh) THEN
             IF (cult.ge.0.20_dp) THEN
                IF (month.GT.4.and.month.LT.9) THEN
                   fertseas=fert
                ENDIF
             ENDIF
          ENDIF
          IF (latitude(jl).lt.mtnh.and.latitude(jl).ge.mtsh) THEN
             IF (cult.ge.0.20_dp) THEN
                fertseas=fert
             ENDIF
          ENDIF
          IF (latitude(jl).lt.mtsh) THEN
             IF (cult.ge.0.20_dp) THEN
                IF (month.LT.3.OR.month.GT.10) THEN
                   fertseas=fert
                ENDIF
             ENDIF
          ENDIF
          fwdagri=fwdagri+fertseas

          !   overlapping problem with agriculture, individual
          !   ecosystems not corrected with cultivation factor. The sum of
          !   the fluxes must be corrected by cult. Make sure, it is in the
          !   output file.

          noslflux_diag(jl, 12) = fwdagri
          noslflux(jl) = noslflux(jl)*(1._dp - cult) + fwdagri*cult

          !  --   correct flux (fwdagri for canopy reduction, month
          !       agriculture

          sai=0._dp
          sai=sai+REAL(sai_veg(12),dp)*noemclass(jl,12)

          crf=(exp(-ks*sai)+exp(-kc*lai(jl)))/2.

          ! Yienger and Levy canopy reduction factor, only applying this
          ! term when the bigleaf approach is being used

          IF (lcrfyl95) THEN
             IF (.not.l_veg_mlay) fwdagri=fwdagri*crf
          ENDIF

          !  --   combine flux from biomes and flux from agriculture
          !       to overall flux, depending of cultivation-fraction

          fwd=fwd*(1._dp-cult)+fwdagri*cult

          !  --   add effect of pulsing to overall flux at each timestep

          !  --   fwd is in [ng N m-2 s-1] whereas the required input is in
          !       molecules NO m-2 s-1, so multiplying with 1e-9 to get g,
          !       dividing by the molecular mass of N to get moles N and
          !       then times avogadro to get molecules

       noflux(jl)=MAX(0.0_dp,((fwd*1.e-9)/amn)*N_A)

    ENDDO kproma_loop

END SUBROUTINE NO_emissions
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE dust_emissions(kproma, time_step_len, wind10, slm, tslm1, prc &
     , prl, emis_du_cla, emis_du_thr, emis_du_src, zprecipinsoil         &
     , duflux )
    ! -----------------------------------------------------------------------
    !
    !     **** dust source  -  calculation of dust flux
    !
    !     Authors:
    !     --------
    !     Michael Schulz; LSCE,Saclay (original source)   9. June 2002
    !     Philip Stier;   MPI-Met, Hamburg, (stier@dkrz.de)
    !                     (adaption to ECHAM5/HAM)           2002/2003
    !
    !     purpose
    !     -------
    !     compute interactively emission flux of mineral dust
    !
    !     interface
    !     ---------
    !      from aerosol modules
    !       zspeed         wind speed at 10 m                 [m/s]
    !       slm            sea land mask                      []
    !       ztemp          surface air temperature            [Å∞C]
    !       zprecip        precipitation amount               [kg m-2 s-1]
    !       zprecipinsoil  proxy for soil humidity            [mm]
    !       densdust       particle density
    !       rhv            source strength factor             [kg s2 m-5]
    !       wth            threshold velocity                 [m/s]
    !       cly            clay content                       [%]
    !
    !       srcsigma       sigma of aerosol mode
    !       srcsigmaln
    !
    !      variables saved
    !       entered        input fields read flag
    !       mask           dust emission area mask            []
    !
    !      output
    !       pxtems         emission of dust per tracer        [kg m-2 s-1]
    !
    !     method
    !     ------
    !
    !     The dust aerosol emitted is described as a lognormal size
    !     distribution with a mass median radius of 2.5 um and a standard
    !     deviation of s = 2. Such a source size distribution was found to
    !     provide a consistent mineral aerosol transport model, which
    !     provides enough mineral dust for long range transport from the
    !     Sahara to the Atlantic and the Mediterranean. Dust plumes
    !     simulated with such a model provide aerosol optical thickness
    !     values very close to those observed with satellites [Schulz et
    !     al., 1998] [Guelle et al., 2000]. Using a fixed lognormal
    !     distribution at the source implies that number concentrations can
    !     be derived from mass emission flux and imposed source size
    !     distribution [Schulz et al., 1998].
    !
    !     While surface wind speed is produced interactively or read from
    !     ECMWF winds (see routine aerosol_meteo_calc), the source function
    !     incorporates also spatial information on three parameters a) the
    !     source strength factor (rhv); b) the threshold velocity (wth) and
    !     c) the clay fraction of the uppermost soil layer (cly).
    !
    !         eflux = max(rhv*(zspeed-wth)*zspeed**2,0.)*landmask
    !
    !     This work is based on an offline dust source developped by
    !     [Claquin, 1999].  He used a global TM3 dust simulation with a
    !     dustsource in all arid regions to obtain an approximate dust
    !     field. In a second step he compared regional averages of TOMS
    !     observations of aerosol index to this simulation, to obtain a
    !     source strength factor for a given region [see also Balkanski et
    !     al. [submitted , 2002]. The threshold velocity was obtained by
    !     comparing the FAO soil database to the maps of Marticorena [1995]
    !     to derive by this procedure a soil specific threshold velocity for
    !     desert soils. Combined with the high resolution FAO dataset this
    !     results in a global map of threshold velocities for dust rise in
    !     arid and semi-arid regions. The offline source formulation was
    !     converted to an interactive dust source for this work [Schulz und
    !     Timmreck, in prep] including prognostic variables produced by the
    !     GCM, such as wind speed, precipitation and surface temperature.
    !
    !     Dust emission is thought to be inhibited, when the surface soil is
    !     wettened or frozen. Drying of desert soils can be parameterised as
    !     a function of recent accumulated precipitation, surface
    !     temperature. The clay content of the soil is a key parameter for
    !     this process, since clay retains water for a longer period. While
    !     Claquin [1999] used an offline calculation to define the periods
    !     when the soil is available for emission of dust, we prefer an
    !     interactive approach here. Soil wettening and drying rate are
    !     computed from GCM calculated preceipitation and surface
    !     temperatures. The clay content is assembled from the FAO soil
    !     data. During freezing conditions, no drying of the soil and no
    !     emissions are assumed.  The parameterisation uses the conditions
    !     given by Claquin [1999] and computes an evaporation rate acting on
    !     the accumulated precipitation amount [Schulz und Timmreck, in
    !     prep]. This might overestimate the time which is needed to dry off
    !     the top soil, since no run-off or losses to the ground water
    !     reservoirs are assumed.
    !     (parameterised M. Schulz LSCE 13.2.2001)
    ! -----------------------------------------------------------------------
    ! Interface:
    ! ----------
    ! input
    ! time_step_len: lenght of the timestep
    ! wind10     : wind speed at 10 m [m/s]
    ! slm        : land-sea mask [0-1]
    ! tslm1      : surface temperature [K]
    ! prc        : convective rainfall [m]
    ! prl        : large-scale rainfall [m]
    ! emis_du_cla: Dust emis. clay content
    ! emis_du_thr: Dust emis. threshold velocity
    ! emis_du_src: Dust emis. source strenght fact.
    !
    ! input/output
    ! zprecipinsoil: proxy for soil humidity [mm]
    !
    ! output
    ! emis_du   : desert dust emission flux for two modes

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)    :: kproma
    REAL(dp), INTENT(in)   :: time_step_len
    REAL(dp), INTENT(in)   :: wind10(1:kproma)
    REAL(dp), INTENT(in)   :: slm(1:kproma)
    REAL(dp), INTENT(in)   :: tslm1(1:kproma)
    REAL(dp), INTENT(in)   :: prc(1:kproma)
    REAL(dp), INTENT(in)   :: prl(1:kproma)
    REAL(dp), INTENT(in)   :: emis_du_cla(1:kproma)
    REAL(dp), INTENT(in)   :: emis_du_thr(1:kproma)
    REAL(dp), INTENT(in)   :: emis_du_src(1:kproma)
    REAL(dp), INTENT(inout):: zprecipinsoil(1:kproma)
    REAL(dp), INTENT(out)  :: duflux(1:kproma)

    ! Local variables

    REAL(dp)  :: avgdryrate     ! Average drying rate [mm dt-1]

    REAL(dp)  :: drying(1:kproma)
    REAL(dp)  :: clyfac(1:kproma)
    REAL(dp)  :: ztempc(1:kproma)
    REAL(dp)  :: zprecip(1:kproma)
    LOGICAL   :: lomask(1:kproma)

    ! INITIALIZATION

    lomask(:)  = .false.
    clyfac(:)  = 0._dp
    ztempc(:)  = 0._dp
    zprecip(:) = 0._dp

    duflux(:) = 0._dp

    ! total precipitation in mm ! This term must be calculated in mm since
    ! it is assigned to zprecipinpsoil which is is given in mm.

    zprecip(:) = (prc(:) + prl(:))*1000. ! mm

    !--- Temperature in Å∞C:

    ztempc(:) = MAX(-10._dp,MIN(50._dp,tslm1(:)-273.17_dp))

    !--- Average drying rate [mm dt-1]:

    avgdryrate=300._dp/365.*time_step_len/86400.

    ! LOGICAL MASK FOR LAND

    lomask(:)=(emis_du_cla(:) > 0._dp .AND. &
         emis_du_thr(:) > 0._dp .AND. &
         slm(:)         > 0.5_dp         )

    !--- 2) Determine soil moisture:

    WHERE (lomask(:))
       !--- Accumulate precipitation:

       ! The parameter zprecip is defined using the convective and large scale
       ! rainfall rates prc and prl [m] whereas the term is not reset  to zero.
       zprecipinsoil(:)=zprecipinsoil(:)+zprecip(:)

       !--- Wet time of the soil depends on clay content:
       !
       ! (New soil wettening formulation f(clay content and soil temperature)
       ! Dependent on clay content and interpolated to obtain
       ! clyfac, which is the maximim number of dry days after a rain event
       ! and corresponds to a maximum amount of water hold in top soil [mm]
       ! at the lowest temperatures [5Å∞C]. The avg drying rate scales the
       ! drying to keep the area an arid area with an aridity limit assumed
       ! to be 300 mm precipitation per year.
       ! More rapid drying only depends on temperature and fits the
       ! 'dry days table' of Claquin to a function with two fitting parameters.

       !--- Maximal amount of water hold in top soil [mm]:

       clyfac(:) = MIN(16._dp,emis_du_cla(:) * 0.4_dp + 8._dp)

       !--- Amount of water drying in top soil during one timestep [mm]:

       drying(:) = avgdryrate &
            * EXP(0.03905491_dp * EXP(0.17446_dp * ztempc(:)))

       !--- Remaining water in top soil [mm]:
       zprecipinsoil(:) = &
            MIN(MAX(0._dp,zprecipinsoil(:) - drying(:)), clyfac(:))

    ENDWHERE

    ! COMPUTATION OF DUST MASS FLUX in [kg m-2 s-1]:

    WHERE ((ztempc(:) > 0._dp) .AND. lomask(:).AND. (zprecipinsoil(:).LE.0._dp))

       duflux(:) = &
            MAX((emis_du_src(:)* (wind10(:) - emis_du_thr(:)) * wind10(:)**2) &
            ,0._dp)

    ENDWHERE

END SUBROUTINE dust_emissions
! -------------------------------------------------------------------------


SUBROUTINE dust_tegen_init(cuscale)

  USE messy_main_constants_mem,  ONLY: g !,iouerr

  ! calculation of the threshold friction velocity Uth

  IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(in)   :: cuscale

    ! LOCAL
    INTEGER             :: ns, nsi, nd, np, j, kk
    REAL(dp)            :: AAA, BB, CCC, DDD, EE, FF, XK, xl, xm, xn, xnV
    REAL(dp)            :: Stotal, StotalV, su, suV, su_loc, su_locV
    REAL(dp)            :: d1, feff
    REAL(dp), DIMENSION(nats)        :: utest
    REAL(dp), DIMENSION(nclass)      :: su_class, su_classV
    REAL(dp), PARAMETER :: gravi=g*100._dp
    REAL(dp), PARAMETER :: a_rnolds=1331.647_dp     ! Reynolds constant
    REAL(dp), PARAMETER :: b_rnolds=0.38194_dp      ! Reynolds constant
    REAL(dp), PARAMETER :: x_rnolds=1.561228_dp     ! Reynolds constant
    REAL(dp), PARAMETER :: aeff=0.35_dp             ! efficient fraction
    REAL(dp), PARAMETER :: xeff=10._dp              ! efficient fraction
    REAL(dp), PARAMETER :: d_thrsld=0.00000231_dp   ! thresold value
    REAL(dp), PARAMETER :: z0S=0.001_dp         ! roughness length of surface...
                 ! without obstacles (cm) (see: Thesis of B. Marticorena, p.85)

    INTEGER :: i ! loop index

       Uth(:) = 0._dp
       i = 0
       rdp = dmin
       DO WHILE(rdp .LE. dmax + 1.E-5_dp)
          i = i + 1
          BB = a_rnolds * (rdp ** x_rnolds) + b_rnolds
          XK = SQRT(rop * gravi * rdp / roa)
          CCC = SQRT(1._dp + d_thrsld /(rdp ** 2.5_dp))
          IF (BB .LT. 10._dp) THEN
             DDD=SQRT(1.928_dp * (BB ** 0.092_dp) - 1._dp)
             Uth(i) = 0.129_dp * XK * CCC / DDD
          ELSE
             EE = -0.0617_dp * (BB - 10._dp)
             FF = 1._dp -0.0858_dp * EXP(EE)
             Uth(i) = 0.12_dp * XK * CCC * FF
          ENDIF
          rdp = rdp * EXP(dstep)
       END DO

       Uth(:)=Uth(:)*cuscale ! scale wind stress threshold

       srel(:,:) = 0._dp
       srelV(:,:) = 0._dp
       su_srelV(:,:) = 0._dp
       utest(:) = 0._dp

       DO ns = 1,nats ! loop over all soil types

          ! calculation of the soil particle distribution and related surfaces
          rdp = dmin
          kk = 0
          Stotal = 0._dp
          StotalV = 0._dp
          su_class(:) = 0._dp
          su_classV(:) = 0._dp
          DO WHILE (rdp.LE.dmax+1.E-5_dp)        ! surface calculations
             kk = kk + 1
             su = 0._dp
             suV = 0._dp
             DO i = 1, nmode
                nd  = ((i - 1) *3 ) + 1
                nsi = nd + 1
                np  = nd + 2
                IF (solspe(ns,nd).EQ.0._dp .or. solspe(ns,nsi).EQ.0._dp &
                     .or. solspe(ns,nsi).EQ.1._dp) THEN
                   IF (solspe(ns,nsi).EQ.0._dp) THEN
!!$                      write(iouerr,*) 'solspe(ns,nsi)=0. or =1.',solspe(ns,nsi),ns,nsi
                   ENDIF
                   su_loc = 0._dp
                   su_locV= 0._dp
                ELSE
                   xk = solspe(ns,np)/(sqrt(2._dp* pi)*log(solspe(ns,nsi)))
                   xl = ((log(rdp)-log(solspe(ns,nd)))**2)/ &
                        (2._dp*(log(solspe(ns,nsi)))**2)
                   xm = xk * exp(-xl)
                   xn =  rop*(2._dp/3._dp)*(rdp/2._dp) !surface
                   xnV =  1._dp !volume
                   su_loc = (xm*dstep/xn)
                   su_locV = (xm*dstep/xnV)
                ENDIF !
                su = su + su_loc
                suV = suV + su_locV
             END DO !nmode
             su_class(kk) = su
             su_classV(kk) = suV
             Stotal = Stotal + su
             StotalV = StotalV + suV
             rdp = rdp * exp(dstep)
          END DO !rdp

          DO j = 1,nclass
             IF (Stotal .EQ. 0._dp) THEN
                srel(ns,j) = 0._dp
                srelV(ns,j) = 0._dp
             ELSE
                srel(ns,j) = su_class(j)/Stotal
                srelV(ns,j) = su_classV(j)/StotalV
                utest(ns) = utest(ns)+srelV(ns,j)
                su_srelV(ns,j) = utest(ns)
             ENDIF
          END DO !j=1,nclass

       END DO !ns (soil type)


       d1   = 0._dp      ! initializations
       feff = 0._dp
       AAA=0._dp
       BB=0._dp
       CCC=0._dp
       DDD=0._dp
       EE=0._dp
       FF=0._dp

       IF (Z01 .EQ. 0._dp) THEN ! calculation of the effective fraction
          feff = 0._dp
       ELSE
          ! partition of energy between the surface and the
          !elements of rugosity (Marticorena thesis, pp 111-112)
          AAA = log(Z01/Z0S)
          BB = log(aeff*(xeff/Z0S)**0.8_dp)
          CCC = 1._dp- AAA/BB
          IF (d1 .EQ. 0._dp) THEN ! partition between Z01 and Z02
             FF = 1._dp
          ELSE
             DDD = log(Z02/Z01)
             EE = log(aeff * (d1/Z01)**0.8_dp)
             FF = 1._dp- DDD/EE
          ENDIF
          feff = FF*CCC                     ! total effective fraction
          IF (feff .LT. 0._dp) feff=0._dp
          IF (feff .GT. 1._dp) feff=1._dp
       ENDIF

       c_eff = feff

END SUBROUTINE dust_tegen_init

! -------------------------------------------------------------------------

SUBROUTINE dust_emissions_tegen(kproma, slf, glac, pcvs, ws, wsmx,         &
                mat_s2, mat_s3, mat_s4, mat_s6, mat_psrc, k_fpar_eff,      &
                wind10_2d, cuscale, duflux, duflux_ai)

    ! -----------------------------------------------------------------------
    !
    !     **** dust source  -  calculation of dust flux
    !
    !
    ! output
    ! emis_du   : desert dust emission flux for two modes

  USE messy_main_constants_mem,  ONLY: g

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(in)    :: kproma
    REAL(dp), INTENT(in)   :: slf(1:kproma)
    REAL(dp), INTENT(in)   :: glac(1:kproma)
    REAL(dp), INTENT(in)   :: pcvs(1:kproma)
    REAL(dp), INTENT(in)   :: ws(1:kproma)
    REAL(dp), INTENT(in)   :: wsmx(1:kproma)
    REAL(dp), INTENT(in)   :: mat_s2(1:kproma)
    REAL(dp), INTENT(in)   :: mat_s3(1:kproma)
    REAL(dp), INTENT(in)   :: mat_s4(1:kproma)
    REAL(dp), INTENT(in)   :: mat_s6(1:kproma)
    REAL(dp), INTENT(in)   :: mat_psrc(1:kproma)
    REAL(dp), INTENT(in)   :: k_fpar_eff(1:kproma)
    REAL(dp), INTENT(in)   :: wind10_2d(1:kproma)
    REAL(dp), INTENT(in)   :: cuscale ! um_gg_20130502
    REAL(dp), INTENT(out)  :: duflux(1:kproma)
    REAL(dp), INTENT(out)  :: duflux_ai(1:kproma)

    ! Local parameters and variables
    ! minimum threshold friction windspeed (cm/s)
    REAL(dp), PARAMETER  :: umin=21._dp
    ! Von Karman constant: 0.4 (0.35 <--> 0.42 see Stull)
    REAL(dp), PARAMETER  :: vk=0.4_dp
    REAL(dp), PARAMETER  :: w0=0.99_dp     ! threshold of relative soil humidity
    REAL(dp), PARAMETER  :: zz=1000._dp    ! wind measurment height (cm)

    ! Bin boundaries:
    INTEGER, PARAMETER :: min_ai = 1
    INTEGER, PARAMETER :: max_ai = 1
    INTEGER, PARAMETER :: min_ci = 2
    INTEGER, PARAMETER :: max_ci = 4

    REAL(dp)           :: cd               ! flux dimensioning parameter
    REAL(dp)           :: mat_s1(kproma)

    ! ECHAM parameter for snow cover calculation
    REAL(dp), PARAMETER  :: zepsec=1.E-12_dp
    ! ECHAM parameter for snow cover calculation
    REAL(dp), PARAMETER  :: zsigfac=0.15_dp

    REAL(dp)           :: alpha

    !REAL(dp)           :: dbmin(ntrace)       ! bin size limits
    !REAL(dp)           :: dbmax(ntrace)       ! bin size limits
    !REAL(dp)           :: dlast               ! bin size limits
    !REAL(dp)           :: rdp
    INTEGER            :: dn(kproma)     ! number of dislocated particle classes
    INTEGER            :: dk(kproma,nclass)! dislocated particle classes

    REAL(dp)           :: du_snow(kproma) ! ECHAM snow coverage interpolated to
                                          ! 0.5x0.5 grid resolution
    REAL(dp)           :: du_W1r(kproma)  ! ECHAM relative soil humidity
                                          !interpolated to
                                          ! 0.5x0.5 grid resolution
    REAL(dp)           :: du_wind(kproma) ! ECHAM wind field interpolated to
                                          ! 0.5x0.5 grid resolution

    INTEGER            :: dust_mask(kproma) ! grid point mask for dust emissions
                                            ! (0: no emission)

    REAL(dp)           :: fluxtyp(kproma,nclass)
    REAL(dp)           :: fluxbin(kproma,ntrace)
    REAL(dp)           :: flux_ai(kproma)  ! dust flux sum for accumulation mode
    REAL(dp)           :: flux_ci(kproma)  ! dust flux sum for coarse mode
    REAL(dp)           :: fdp1, fdp2

    REAL(dp)           :: fluxdiam1(kproma,nclass)  ! flux for soil type #1
    REAL(dp)           :: fluxdiam2(kproma,nclass)  ! flux for soil type #2
    REAL(dp)           :: fluxdiam3(kproma,nclass)  ! flux for soil type #3
    REAL(dp)           :: fluxdiam4(kproma,nclass)  ! flux for soil type #4
    REAL(dp)           :: fluxdiam6(kproma,nclass)  ! flux for soil type #6
    REAL(dp)           :: fluxdiam_pf(kproma,nclass)! flux for preferential
                                                    !     sources soil type

    INTEGER            :: i                        ! loop index
    INTEGER            :: i_soil                   ! soil type

    INTEGER            :: kk,kkk                   ! loop index
    INTEGER            :: kkmin

    INTEGER            :: n,nn                     ! loop index

    REAL(dp)           :: Ustar,Ustar_d(kproma)    ! threshold friction velocity
                                                   ! for saltation
    REAL(dp)           :: uthp

    REAL(dp)           :: zw1r(kproma)             ! relative soil moisture
    REAL(dp)           :: zcvs(kproma)             ! snow cover fraction

    REAL(dp)           :: flux_6h(kproma,ntrace)

    ! INITIALIZATION
    cd = 1.00_dp*roa/(g*100._dp)
    mat_s1(:) = 1._dp


    du_wind(:) = wind10_2d(1:kproma)

    WHERE (du_wind(:) .LT. 0._dp) du_wind(:) = 0._dp
    ! set missing values to zero

    ! calculate ECHAM snow cover fraction *zcvs*
    ! (algorithm taken from routine *physc.f90*)
    ! and ECHAM relative soil moisture *zw1r* (ratio *wl* to *wlmax*)
    DO i=1,kproma
       ! land surface, but not a glacier
      IF ((slf(i).GT.0.5_dp).AND.(glac(i).LT.0.5_dp)) THEN
        zw1r(i)= MIN(ws(i)/wsmx(i),1._dp)
        zcvs(i) = pcvs(i)
      ELSEIF((slf(i).GT.0.5_dp).AND.(glac(i).GT.0.5_dp)) THEN
                    ! glacier on land surface
        zw1r(i)= 1._dp
        zcvs(i)=1._dp
      ELSE          ! ocean grid point
        zw1r(i)= -9.e9_dp
        zcvs(i)=-9.e9_dp
      END IF

    END DO !i

    du_W1r(:) = zw1r(:)
    du_snow(:) = zcvs(:)

    ! set missing values to zero
    WHERE (du_W1r(:)  .LT. 0._dp) du_W1r(:)  = 0._dp
    WHERE (du_snow(:) .LT. 0._dp) du_snow(:) = 0._dp
    IF (c_eff.GT.0._dp) THEN

!!$     rdp=dmin
     !dlast=dmin
!!$     nn=1
!!$     !dbmin(:)=0._dp
!!$     !dbmax(:)=0._dp
!!$     DO kk=1,nclass ! assign fluxes to bins
!!$      IF (mod(kk,nbin).eq.0) THEN
!!$       !dbmax(nn)=rdp*10000._dp*0.5_dp
!!$       ! calculate bin minimum/maximum radius in um
!!$       !dbmin(nn)=dlast*10000._dp*0.5_dp
!!$       nn=nn+1
!!$       !dlast=rdp
!!$      ENDIF
!!$      rdp = rdp * exp(Dstep)
!!$     ENDDO !kk

     fluxbin(:,:)=0._dp
     flux_6h(:,:)=0._dp

     dust_mask(:)=0

     DO i = 1,kproma
       Ustar = (VK * du_wind(i) *100._dp)/(log(ZZ/Z02))  ! set wind speed
       IF (Ustar.gt.0.and.(Ustar.GE. umin*cuscale/c_eff)) THEN
          ! check critical wind speed (wind stress threshold scaled)
          IF (k_fpar_eff(i).GT.1.E-10_dp) THEN
             ! check if the grid cell is a potential dust source
             dust_mask(i)=1
             ! set grid point as a potential dust source point
             Ustar_d(i)=Ustar ! store wind speed of dust grid
          ENDIF ! k_fpar_eff.gt.0.
       ENDIF   ! Ustar
    END DO !i

    fluxtyp(:,:)=0._dp
    fluxdiam1(:,:)=0._dp
    fluxdiam2(:,:)=0._dp
    fluxdiam3(:,:)=0._dp
    fluxdiam4(:,:)=0._dp
    fluxdiam6(:,:)=0._dp
    fluxdiam_pf(:,:)=0._dp
    dn(:)=0

    DO kk=1,nclass
       Uthp=Uth(kk)
       DO i = 1,kproma
          IF (dust_mask(i).EQ.1) THEN
             ! Marticorena:
             fdp1 = (1._dp-(Uthp/(c_eff * UStar_d(i))))
             fdp2 = (1._dp+(Uthp/(c_eff * UStar_d(i))))**2
             ! Shao:
             ! fdp1 = (1.-(Uthp/(c_eff * UStar_d(i)))**2)
             ! fdp2 = 1.
             IF (fdp1.gt.0._dp) THEN
                i_soil=1
                alpha= solspe(i_soil,nmode*3+1)
                fluxdiam1(i,kk) = srel(i_soil,kk) * fdp1 * fdp2 * cd * &
                     UStar_d(i)**3 *alpha ! flux for soil type #1
                i_soil=2
                alpha= solspe(i_soil,nmode*3+1)
                fluxdiam2(i,kk) = srel(i_soil,kk) * fdp1 * fdp2 * cd * &
                     UStar_d(i)**3 *alpha ! flux for soil type #2
                i_soil=3
                alpha= solspe(i_soil,nmode*3+1)
                fluxdiam3(i,kk) = srel(i_soil,kk) * fdp1 * fdp2 * cd * &
                     UStar_d(i)**3 *alpha ! flux for soil type #3
                i_soil=4
                alpha= solspe(i_soil,nmode*3+1)
                fluxdiam4(i,kk) = srel(i_soil,kk) * fdp1 * fdp2 * cd * &
                     UStar_d(i)**3 *alpha ! flux for soil type #4
                i_soil=6
                alpha= solspe(i_soil,nmode*3+1)
                fluxdiam6(i,kk) = srel(i_soil,kk) * fdp1 * fdp2 * cd * &
                     UStar_d(i)**3 *alpha ! flux for soil type #6
                i_soil=10
                alpha= solspe(i_soil,nmode*3+1)
                ! fluxdiam_pf = potential flux if soil type would be a
                !               preferential source
                fluxdiam_pf(i,kk)= srel(i_soil,kk) * fdp1 * fdp2 * cd *  &
                     UStar_d(i)**3 *alpha
                IF (kk.eq.1) THEN
                   fluxtyp(i,kk) = fluxtyp(i,kk)                             &
                        + fluxdiam1(i,kk)*(1._dp-mat_psrc(i))*               &
                        (mat_s1(i)-mat_s2(i)-mat_s3(i) -mat_s4(i)-mat_s6(i)) &
                        + fluxdiam2(i,kk)*(1._dp-mat_psrc(i))*mat_s2(i)      &
                        + fluxdiam3(i,kk)*(1._dp-mat_psrc(i))*mat_s3(i)      &
                        + fluxdiam4(i,kk)*(1._dp-mat_psrc(i))*mat_s4(i)      &
                        + fluxdiam6(i,kk)*(1._dp-mat_psrc(i))*mat_s6(i)      &
                        + fluxdiam_pf(i,kk)*mat_psrc(i)
                ELSE
                   dn(i)=dn(i)+1
                   ! increase number of dislocated particle classes
                   dk(i,dn(i))=kk ! store dislocated particle class
                ENDIF ! kk.eq.1
             ENDIF ! fdp1.gt.0.
          ENDIF ! dust_mask.eq.1
       ENDDO ! i
    ENDDO ! kk

    kkmin=1
    DO i = 1,kproma
       IF (dust_mask(i).EQ.1) THEN
          ! loop over dislocated dust particle classes
          DO n=1,dn(i)
             ! scaling with relative contribution of dust size  fraction
             DO kkk=1,dk(i,n)
                i_soil=1
                fluxtyp(i,kkk) = fluxtyp(i,kkk) + (1._dp-mat_psrc(i))* &
                     (mat_s1(i)-mat_s2(i)-mat_s3(i)-mat_s4(i)-mat_s6(i))  &
                     * fluxdiam1(i,dk(i,n))*srelV(i_soil,kkk) / &
                     ((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
                i_soil=2
                fluxtyp(i,kkk) = fluxtyp(i,kkk) + (1._dp-mat_psrc(i))*mat_s2(i) &
                     * fluxdiam2(i,dk(i,n))*srelV(i_soil,kkk) / &
                     ((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
                i_soil=3
                fluxtyp(i,kkk) = fluxtyp(i,kkk) + (1._dp-mat_psrc(i))*mat_s3(i) &
                     * fluxdiam3(i,dk(i,n))*srelV(i_soil,kkk) / &
                     ((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
                i_soil=4
                fluxtyp(i,kkk) = fluxtyp(i,kkk) + (1._dp-mat_psrc(i))*mat_s4(i) &
                     * fluxdiam4(i,dk(i,n))*srelV(i_soil,kkk) / &
                     ((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
                i_soil=6
                fluxtyp(i,kkk) = fluxtyp(i,kkk) + (1._dp-mat_psrc(i))*mat_s6(i) &
                     * fluxdiam6(i,dk(i,n))*srelV(i_soil,kkk) / &
                     ((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
                IF(du_wind(i).gt.10._dp) THEN
                   i_soil=11 ! flux from preferential source at high wind speeds
                ELSE
                   i_soil=10 ! flux from preferential source at low wind speeds
                ENDIF
                fluxtyp(i,kkk) = fluxtyp(i,kkk) + mat_psrc(i)   &
                     * fluxdiam_pf(i,dk(i,n))*srelV(i_soil,kkk) / &
                     ((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
             ENDDO ! kkk
          ENDDO ! n
       ENDIF ! dust_mask.eq.1
    ENDDO ! i

    DO nn=1,ntrace
       DO i = 1,kproma
          IF (dust_mask(i).EQ.1) THEN
            fluxbin(i,nn) = fluxbin(i,nn)+SUM(fluxtyp(i,(nn-1)*nbin+1:nn*nbin))
             IF (du_W1r(i).gt.w0) fluxbin(i,nn)=0._dp
             ! mask out dust fluxes, where soil moisture threshold is reached
             flux_6h(i,nn)=fluxbin(i,nn)*10000._dp* &
                  ! fluxbin: g/cm2/sec; flux_6h: g/m2/sec
                  (1._dp-du_snow(i))*k_fpar_eff(i)
          ENDIF ! dust_mask.eq.1
       ENDDO ! i
    ENDDO ! nn

    ENDIF ! c_eff

    !--- Sum over internal tracer classes for HAM modal scheme:

    flux_ai(:) = 0._dp
    flux_ci(:) = 0._dp

    ! Accumulation mode:

    DO nn=min_ai, max_ai
       flux_ai(:) = flux_ai(:) + flux_6h(:kproma,nn)
    END DO

    ! Coarse mode:

    DO nn=min_ci, max_ci
       flux_ci(:)=flux_ci(:)+flux_6h(:kproma,nn)
    ENDDO !nn

    WHERE (flux_ai.le.0) flux_ai = 0._dp  ! clean up fluxes
    WHERE (flux_ci.le.0) flux_ci = 0._dp  ! clean up fluxes

    WHERE (glac(1:kproma) > 0.5_dp .OR. slf(1:kproma) < 0.5_dp)
       ! mask out glacier and ocean values on ECHAM grid
       flux_ai(:) = 0._dp
       flux_ci(:) = 0._dp
    ENDWHERE

    ! DUST MASS FLUX in [kg m-2 s-1]:
    duflux(:)    = flux_ci(:kproma)/1000._dp
    duflux_ai(:) = flux_ai(:kproma)/1000._dp

END SUBROUTINE dust_emissions_tegen
! -------------------------------------------------------------------------

SUBROUTINE dust_emissions_KKDU_Astitha1(kproma, airdens, zust_2d, slf, &
      cvs, kkdu_clay, kkdu_topo, kkdu_nap, kkdu_kp, kkdu_capp,         &
      kkdu_mgpp, kkdu_misc, kkdu_mask, kkdu_lai, duflux_ai, duflux_ci, &
      kkdu_nap_emflux_ai, kkdu_nap_emflux_ci, kkdu_kp_emflux_ai,       &
      kkdu_kp_emflux_ci, kkdu_capp_emflux_ai, kkdu_capp_emflux_ci,     &
      kkdu_mgpp_emflux_ai, kkdu_mgpp_emflux_ci, kkdu_misc_emflux_ai,   &
      kkdu_misc_emflux_ci, horflux, ustarthr)

! ----------------------------------------------------------------------------
  !
  !     *************** Online calculation of dust fluxes  ********************
  !
  !     Authors:
  !     --------
  !     Written by Marina Astitha for MESSy1.8; The Cyprus Institute-EEWRC
  !                                             September 2009++
  !     Implemented in MESSy2.41 by Mohamed Abdel Kader (2011)
  !     Revised by Klaus Klingmueller 2017 (Klingmueller et al. 2018)
  !
  !     purpose
  !     -------
  !     Emission Fluxes of Mineral Dust depending on soil,friction velocity and
  !     soil moisture. More Information on the emission scheme can be found in
  !     Astitha et al. (2012) in ACPD (DU1 scheme)
  !     and Klingmueller et al. (2018)
  !
  !
  !     INPUT
  !     ---------
  !       airdens        air density                                     [kg/m3]
  !       kkdu_clay      clay content                                    [%]
  !       kkdu_topo      topography factor
  !       kkdu_mask      dust emission mask
  !       kkdu_lai       leaf/vegetation area index
  !       kkdu_nap       Na+ fraction in desert soil
  !       kkdu_kp        K+ fraction in desert soil
  !       kkdu_capp      Ca++ fraction in desert soil
  !       kkdu_mgpp      Mg++ fraction in desert soil
  !       kkdu_misc      misc. fraction in desert soil
  !       zust_2d        surface friction velocity                       [m/sec]
  !       slf            fractional land-sea mask [0=sea and 1=land]
  !       qm1            specific humidity at 1st model layer            [kg/kg]
  !       qte            specific humidity tendency                  [kg/kg*sec]
  !       tm1            temperature at 1st model layer                  [K]
  !       cvs            fraction of snow cover
  !
  !      OUTPUT
  !       duflux_ai      Vertical flux (emissions) of dust - accumulation
  !                      [kg m-2 s-1]
  !       duflux_ci      Vertical flux (emissions) of dust - coarse
  !                      [kg m-2 s-1]
  !       kkdu_nap_emflux_ai
  !                      Na+ emissions, accumulation mode [kg m-2 s-1]
  !       kkdu_nap_emflux_ci
  !                      Na+ emissions, coarse mode [kg m-2 s-1]
  !       kkdu_kp_emflux_ai
  !                      K+ emissions, accumulation mode [kg m-2 s-1]
  !       kkdu_kp_emflux_ci
  !                      K+ emissions, coarse mode [kg m-2 s-1]
  !       kkdu_capp_emflux_ai
  !                      Ca++ emissions, accumulation mode [kg m-2 s-1]
  !       kkdu_capp_emflux_ci
  !                      Ca++ emissions, coarse mode [kg m-2 s-1]
  !       kkdu_mgpp_emflux_ai
  !                      Mg++ emissions, accumulation mode [kg m-2 s-1]
  !       kkdu_mgpp_emflux_ci
  !                      Mg++ emissions, coarse mode [kg m-2 s-1]
  !       kkdu_misc_emflux_ai
  !                      misc. dust, accumulation mode [kg m-2 s-1]
  !       kkdu_misc_emflux_ci
  !                      misc. dust, coarse mode [kg m-2 s-1]
  !
  !      DIAGNOSTIC OUTPUT
  !       alpha         !sandblasting efficiency [cm-1]
  !       ustarthr      !threshold wind friction velocity [m/s]
  !       baresoil      !fraction of bare soil exposed in a grid cell [-]
  !       horflux       !horizontal flux of dust from soil sources (saltation)
  !                     [kg m-1 s-1]
  !
  ! -----------------------------------------------------------------------
  !
  !    Transport size bins (Perez et al. 2006)
  !    Source Size distribution (d'Almeida 1987;Zender et al. 2003)
  !      RADIUS(min-max)   Effective radius DIAMETER
  !    1) 0.1  -  0.18 um - reff=0.15 um    Dn(um)   Dv(um)   sigma    Mfraction
  !    2) 0.18 -  0.3  um - reff=0.25 um    0.16     0.832     2.10      0.036
  !    3) 0.3  -  0.6  um - reff=0.45 um    1.40     4.820     1.90      0.957
  !    4) 0.6  -  1.0  um - reff=0.78 um   10.00    19.380     1.60      0.007
  !    5) 1.0  -  1.8  um - reff=1.3  um
  !    6) 1.8  -  3.0  um - reff=2.2  um
  !    7) 3.0  -  6.0  um - reff=3.8  um
  !    8) 6.0  - 10.0  um - reff=7.1  um
  ! ----------------------------------------------------------------------------
  USE messy_main_constants_mem,  ONLY: g
  USE messy_main_tools,          ONLY: ERRFUNC

  IMPLICIT NONE

  ! The following parameters are hard coded in GMXE on smcl level. They should
  ! be written to a GMXE channel object to be accesible here. Until then make
  ! sure the following values are consistent with the GMXE values.
  REAL(dp) :: crdiv(4) = &
       (/ 0.0005e-4_dp, 0.006e-4_dp, 0.06e-4_dp, 1.0e-4_dp /) ! in cm
  INTEGER, PARAMETER :: as = 3, cs = 4

  !g       = 9.80665_dp   ! gravity acceleration [m/s2]

  ! I/O
  INTEGER, PARAMETER     :: NPS=3                !Number of source size bins
  INTEGER, PARAMETER     :: NPT=2                !Number of EMAC modes
  INTEGER, INTENT(IN)    :: kproma
  REAL(dp), INTENT(IN)   :: airdens(1:kproma)
  REAL(dp), INTENT(IN)   :: kkdu_clay(1:kproma)
  REAL(dp), INTENT(IN)   :: kkdu_topo(1:kproma)
  REAL(dp), INTENT(IN)   :: kkdu_nap(1:kproma)
  REAL(dp), INTENT(IN)   :: kkdu_kp(1:kproma)
  REAL(dp), INTENT(IN)   :: kkdu_capp(1:kproma)
  REAL(dp), INTENT(IN)   :: kkdu_mgpp(1:kproma)
  REAL(dp), INTENT(IN)   :: kkdu_misc(1:kproma)
  REAL(dp), INTENT(IN)   :: kkdu_lai(1:kproma)
  REAL(dp), INTENT(IN)   :: kkdu_mask(1:kproma,1:2)
  REAL(dp), INTENT(IN)   :: zust_2d(1:kproma)
  REAL(dp), INTENT(IN)   :: cvs(1:kproma)
  REAL(dp), INTENT(IN)   :: slf(1:kproma)
!
  REAL(dp)               :: kappa1(1:kproma),kappa2
  REAL(dp)               :: zos,zo,fdrag,factor(1:kproma),pdiam
  REAL(dp)               :: mij(NPS,NPT),massfraction(NPT),sumtest
!
  REAL(dp)               :: duflux2(1:kproma,NPT)

  REAL(dp), INTENT(OUT)  :: duflux_ai(1:kproma),duflux_ci(1:kproma)
  REAL(dp), INTENT(OUT)  :: kkdu_nap_emflux_ai(1:kproma)
  REAL(dp), INTENT(OUT)  :: kkdu_nap_emflux_ci(1:kproma)
  REAL(dp), INTENT(OUT)  :: kkdu_kp_emflux_ai(1:kproma)
  REAL(dp), INTENT(OUT)  :: kkdu_kp_emflux_ci(1:kproma)
  REAL(dp), INTENT(OUT)  :: kkdu_capp_emflux_ai(1:kproma)
  REAL(dp), INTENT(OUT)  :: kkdu_capp_emflux_ci(1:kproma)
  REAL(dp), INTENT(OUT)  :: kkdu_mgpp_emflux_ai(1:kproma)
  REAL(dp), INTENT(OUT)  :: kkdu_mgpp_emflux_ci(1:kproma)
  REAL(dp), INTENT(OUT)  :: kkdu_misc_emflux_ai(1:kproma)
  REAL(dp), INTENT(OUT)  :: kkdu_misc_emflux_ci(1:kproma)
  REAL(dp), INTENT(OUT)  :: horflux(1:kproma)
  REAL(dp)               :: alpha(1:kproma)
  REAL(dp), INTENT(INOUT):: ustarthr(1:kproma)
  REAL(dp)               :: betar(1:kproma)
  REAL(dp)               :: kappa(1:kproma)
  REAL(dp)               :: baresoil(1:kproma)
  REAL                   :: X1,X2,ERFTERM1,ERFTERM2
  ! bulk density of soil particle (kg/m3)
  REAL(dp), PARAMETER    :: bulkdens=2.65E+03_dp
  ! kinematic air viscosity (m2/sec)
  REAL(dp), PARAMETER    :: airvisc=0.157E-04_dp
  ! leaf area index threshold for complete supression of
  ! dust emissions[m2/m2](Mahowald et al. 1999)
  REAL(dp), PARAMETER    :: vlait=0.35_dp
  INTEGER                :: I,J,bio

!*******************************************************************************
!  Olson Biomes used for dust sources: 8,50,51,52
!    --> renamed to 1,2,3,4 respectively (0 for all others)
!
!  Dustmask is 1= biomes  8, 50 and 71  and  0.5= biomes 51, 52
!*******************************************************************************
    REAL(dp), DIMENSION(6) :: dustmask = &              !fraction
       (/0._dp,1._dp,1._dp,0.5_dp,0.5_dp,1._dp/)

!*******************************************************************************
!   2 RELEVANT EMAC MODES (ACCUMULATION & COARSE)
!*******************************************************************************
    REAL(dp), DIMENSION(2) :: mindiam                       ! micro-m
    REAL(dp), DIMENSION(2) :: maxdiam                       ! micro-m
    REAL(dp), DIMENSION(2) :: pdens = (/2.65_dp, 2.65_dp/)  ! g / cm^3

!*******************************************************************************
!   SOIL (SOURCE) SIZE DISTRIBUTION PARAMETERS
!   (d'Almeida 1987,Zender et al. 2003))
!*******************************************************************************
!
    REAL(dp), DIMENSION(3) :: dv = &                   !Mass median Diameter (um)
       (/0.832_dp,4.82_dp,19.38_dp/)
    REAL(dp), DIMENSION(3) :: massfr = &               !Mass Fraction (-)
       (/0.036_dp,0.957_dp,0.007_dp/)
    REAL(dp), DIMENSION(3) :: sigma = &                !Standard deviation
       (/2.1_dp,1.9_dp,1.6_dp/)
!
    ! for the calculation of the standard error function
    INTEGER, PARAMETER :: jint=0

  REAL(dp) :: zust_limited(1:kproma)

!   INITIALIZATION
    mindiam(1) = 2._dp * crdiv(as) * 1.e4_dp            ! cm to micro-m, r to d
    mindiam(2) = 2._dp * crdiv(cs) * 1.e4_dp            ! cm to micro-m, r to d
    maxdiam(1) = 2._dp * crdiv(as + 1) * 1.e4_dp        ! cm to micro-m, r to d
    maxdiam(2) = 100._dp

    DO J=1,NPT
       duflux2(:,J) = 0._dp
       massfraction(J)= 0._dp
    ENDDO
    baresoil(:)  = 0._dp
    alpha(:)     = 0._dp
    horflux(:)   = 0._dp
    duflux_ai(:) = 0._dp
    duflux_ci(:) = 0._dp
    kkdu_nap_emflux_ai(:) = 0._dp
    kkdu_nap_emflux_ci(:) = 0._dp
    kkdu_kp_emflux_ai(:) = 0._dp
    kkdu_kp_emflux_ci(:) = 0._dp
    kkdu_capp_emflux_ai(:) = 0._dp
    kkdu_capp_emflux_ci(:) = 0._dp
    kkdu_mgpp_emflux_ai(:) = 0._dp
    kkdu_mgpp_emflux_ci(:) = 0._dp
    kkdu_misc_emflux_ai(:) = 0._dp
    kkdu_misc_emflux_ci(:) = 0._dp
!
!   ---------------------------------------------------------------------------
!   CALCULATION OF THE SANDBLASTING EFFICIENCY
!      (a=Vertical Flux/Horizontal Flux) - [cm-1]
!
    DO i = 1, kproma
        alpha(i) = sb_efficiency(kkdu_clay(i))
    END DO

!   ----------------------------------------------------------------------------
!   CALCULATION OF THE FRICTION REYNOLDS NUMBER
!   [ B=USTARTHR*DIAMETER/AIRVISC [dimensionless] ]
!
    IF (MAXVAL(ustarthr(:)) <= 0._dp) THEN
       pdiam=60._dp   !Optimum soil size for saltation [um]
       betar(:)=1331._dp*((pdiam*1.e-04)**1.56)+0.38_dp  !Beta=0.835
    ELSE
       pdiam=60._dp   ![um]
       betar(:)= (ustarthr(:)*(pdiam*1.e-06))/airvisc
    ENDIF
!
!   ----------------------------------------------------------------------------
!   CALCULATION OF THE KAPPA COEFFICIENTS THAT ENTER THE FORMULA OF THE
!   THRESHOLD FRICTION VELOCITY USTARTHR
!   (Marticorena and Bergametti (1995);Marticorena et al. (1997))
!
    kappa1(:)=SQRT((pdens(1)*g*pdiam*1.e-06)/(airdens(:)*1.e-03)) ![m/s]
    kappa2=SQRT(1.+(0.006/(pdens(1)*g*100.*(pdiam*1.e-04)**2.5)))![dimensionless]
    kappa(:)=kappa1(:)*kappa2                                     ![m/s]
!
!   CALCULATION OF THE THRESHOLD FRICTION VELOCITY
!
    WHERE((betar(:) > 0.03_dp).AND.(betar(:) < 10._dp))          ! (0.03<B<10)

       ustarthr(:) = 0.129_dp*kappa(:)/SQRT((1.928_dp*betar(:)**0.092)-1.)
       ! [m/s]
    ELSEWHERE (betar(:) >= 10._dp)                               ! (B>10)
       ustarthr(:) = 0.129_dp*kappa(:)*&
            (1.-0.0858_dp*EXP(-0.0617_dp*(betar(:)-10.))) ![m/s]
    ENDWHERE
!
!  -----------------------------------------------------------------------------
!   DRAG PARTITION SCHEME TO QUANTIFY THE FRACTION OF THE TOTAL WIND SHEAR
!   STRESS ON THE ERODIBLE SURFACES TO MOBILIZE THE SOIL PARTICLES
!
!     fdrag = Ratio of local to total friction velocity
!
    zos = 0.00333_dp   !Local roughness length of the uncovered surface
                       !(smooth roughness length) [cm]
!                       (Zender etal 03;MaB95 suggest 1.e-03cm)
!
    zo  = 0.01_dp   !Surface aeolian Roughness Length  [cm]   (Zender etal 2003)
!
    fdrag=1.-(LOG(zo/zos)/LOG(0.35_dp*((10/zos)**0.8)))

    ustarthr(:)=ustarthr(:)/fdrag    ![m/s]
!
!
!  -----------------------------------------------------------------------------
!  CALCULATION OF THE HORIZONTAL DUST EMISSIONS FLUX [HORFLUX(:)]
!  (SALTATION IS INITIATED WHEN U* > U*thresh)
!
!   The value 2.61 comes from the wind tunnel tests from White (1979)
!   Substitute 2.61 with 1 after Darmenova et al. (2009)  - M.Astitha 14sep2011
!
   zust_limited(:) = MIN(zust_2d(:), .4_dp)
   WHERE (zust_limited(:) > ustarthr(:))
      factor(:)= (airdens(:)/g)*(zust_limited(:)**3)    ![kg/m*sec]
      horflux(:)= factor(:)*(1+(ustarthr(:)/zust_limited(:))) &
         *(1-(ustarthr(:)**2/zust_limited(:)**2))  !Marticorena and Bergametti (1995)
       horflux(:) = horflux(:) * kkdu_topo(:) ! topography factor
       horflux(:) = horflux(:) * 7.9_dp ! normalisation/tuning
   ELSEWHERE
      horflux(:)=0.
   ENDWHERE
!
!  -----------------------------------------------------------------------------
!  CALCULATION OF THE MASSFRACTION (Source Sizes To Transport Sizes)
!  [MASSFRACTION(J)]
!
!  Note: The array mij(I,J) is the fraction of the size distribution I in the
!        source groups that overlaps the transport bin J.
!        In that way we compute the mapping from the lognormal source
!        distribution to the transport size bins.
!
!
   sumtest=0.0
   DO I=1,NPS  !Loop over the source modes
      DO J=1,NPT !Loop over the transport bins
!
         X1=REAL((LOG((maxdiam(J))/dv(I)))/(SQRT(2.)*LOG(sigma(I))))
         X2=REAL((LOG((mindiam(J))/dv(I)))/(SQRT(2.)*LOG(sigma(I))))

         CALL ERRFUNC(X1,ERFTERM1,JINT) ! Calculate standard error function
         CALL ERRFUNC(X2,ERFTERM2,JINT)
!
         mij(I,J)=0.5 * (ERFTERM1 - ERFTERM2) !Eq. 12 of Zender et al. 2003
         mij(I,J)=mij(I,J) * massfr(I)

         X1 = 0.0
         X2 = 0.0
         ERFTERM1 = 0.0
         ERFTERM2 = 0.0

         sumtest=sumtest+mij(I,J)

      ENDDO
   ENDDO

   IF(sumtest.lt.0.9.or.sumtest.gt.1.) THEN
      WRITE(*,*)'Problem with the total mass overlap= ', sumtest
   ENDIF

   DO J=1,NPT
      DO I=1,NPS
         massfraction(J)= massfraction(J) + mij(I,J)
      ENDDO
   ENDDO
   !
   !  --------------------------------------------------------------------------
   !  CALCULATION OF THE VERTICAL FLUX OF DUST EMISSIONS [DUFLUX2(:)]
   !

   DO bio=1, 2
      baresoil(:) = baresoil(:) + kkdu_mask(:,bio)*dustmask(bio)
   ENDDO

   WHERE((slf(:) > 0.5_dp).and.(cvs(:).eq.0._dp))
      baresoil(:) = &
           baresoil(:) * slf(:) * (1._dp-min(1._dp,min(kkdu_lai(:),vlait))/vlait)
   ELSEWHERE
      baresoil(:)=0._dp
   ENDWHERE
!
    DO J=1,NPT      !Number of transport size bins
       duflux2(:,J) = 1.0E-04 * (alpha(:)*100.) * baresoil(:) &
            * horflux(:) * massfraction(J) ![kg/m2*sec]
    ENDDO           !NPT
!
!
    !Accumulation insoluble mode of dust
    duflux_ai(:) = duflux2(:, 1)
!
    !Coarse insoluble mode of dust
    duflux_ci(:) = duflux2(:, 2)
!
    kkdu_nap_emflux_ai(:) = duflux_ai(:) * kkdu_nap(:)
    kkdu_nap_emflux_ci(:) = duflux_ci(:) * kkdu_nap(:)
    kkdu_kp_emflux_ai(:) = duflux_ai(:) * kkdu_kp(:)
    kkdu_kp_emflux_ci(:) = duflux_ci(:) * kkdu_kp(:)
    kkdu_capp_emflux_ai(:) = duflux_ai(:) * kkdu_capp(:)
    kkdu_capp_emflux_ci(:) = duflux_ci(:) * kkdu_capp(:)
    kkdu_mgpp_emflux_ai(:) = duflux_ai(:) * kkdu_mgpp(:)
    kkdu_mgpp_emflux_ci(:) = duflux_ci(:) * kkdu_mgpp(:)
    kkdu_misc_emflux_ai(:) = duflux_ai(:) * kkdu_misc(:)
    kkdu_misc_emflux_ci(:) = duflux_ci(:) * kkdu_misc(:)

  CONTAINS
    ! Sand blasting efficiency in 1 / cm
    ! Input: clay fraction in %
    REAL(dp) FUNCTION sb_efficiency(clay_fraction)
    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: clay_fraction
    ! bandwidth 10
    REAL(dp), PARAMETER, DIMENSION(0:100) :: lut = (/ 3.34954e-06_DP, 3.96634e-06_DP,&
    4.80776e-06_DP, 5.96505e-06_DP, 7.56554e-06_DP, 9.78478e-06_DP, 1.28349e-05_DP,&
    1.70014e-05_DP, 2.26247e-05_DP, 3.00716e-05_DP, 3.96797e-05_DP, 5.16538e-05_DP,&
    6.59188e-05_DP, 8.19545e-05_DP, 9.86693e-05_DP, 0.000114397_DP, 0.000127081_DP,&
    0.000134664_DP, 0.000135607_DP, 0.000129358_DP, 0.000116586_DP, 9.90698e-05_DP,&
    7.9252e-05_DP, 5.96277e-05_DP, 4.21897e-05_DP, 2.81066e-05_DP, 1.7694e-05_DP,&
    1.06138e-05_DP, 6.17332e-06_DP, 3.59925e-06_DP, 2.21797e-06_DP, 1.53081e-06_DP,&
    1.21316e-06_DP, 1.07583e-06_DP, 1.01891e-06_DP, 9.96876e-07_DP, 9.93174e-07_DP,&
    9.86076e-07_DP, 9.73418e-07_DP, 9.52422e-07_DP, 9.20024e-07_DP, 8.73522e-07_DP,&
    8.11436e-07_DP, 7.34327e-07_DP, 6.45246e-07_DP, 5.49516e-07_DP, 4.5382e-07_DP,&
    3.64835e-07_DP, 2.87866e-07_DP, 2.25937e-07_DP, 1.79586e-07_DP, 1.47317e-07_DP,&
    1.26419e-07_DP, 1.1383e-07_DP, 1.06776e-07_DP, 1.03098e-07_DP, 1.01316e-07_DP,&
    1.00511e-07_DP, 1.00174e-07_DP, 1.00042e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP,&
    1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP,&
    1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP,&
    1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP,&
    1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP,&
    1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP /)
    ! bandwidth 20
    REAL(dp), PARAMETER, DIMENSION(0:100) :: lut2 = (/ 1.77259e-05_DP, &
    2.02044e-05_DP, 2.3021e-05_DP, 2.61952e-05_DP, 2.97375e-05_DP, &
    3.3645e-05_DP, 3.78976e-05_DP, 4.24537e-05_DP, 4.7247e-05_DP, &
    5.21832e-05_DP, 5.71405e-05_DP, 6.19709e-05_DP, 6.65054e-05_DP, &
    7.05621e-05_DP, 7.39575e-05_DP, 7.65194e-05_DP, 7.81009e-05_DP, &
    7.85936e-05_DP, 7.79386e-05_DP, 7.61326e-05_DP, 7.32309e-05_DP, &
    6.93432e-05_DP, 6.46263e-05_DP, 5.92714e-05_DP, 5.34898e-05_DP, &
    4.74977e-05_DP, 4.15013e-05_DP, 3.56846e-05_DP, 3.02003e-05_DP, &
    2.51641e-05_DP, 2.0653e-05_DP, 1.67066e-05_DP, 1.33313e-05_DP, &
    1.05067e-05_DP, 8.19194e-06_DP, 6.33264e-06_DP, 4.8674e-06_DP, &
    3.73309e-06_DP, 2.86902e-06_DP, 2.21983e-06_DP, 1.73722e-06_DP, &
    1.3806e-06_DP, 1.1171e-06_DP, 9.20949e-07_DP, 7.72617e-07_DP, &
    6.57743e-07_DP, 5.66112e-07_DP, 4.90677e-07_DP, 4.26717e-07_DP, &
    3.71134e-07_DP, 3.24793e-07_DP, 2.8804e-07_DP, 2.55087e-07_DP, &
    2.26072e-07_DP, 2.00985e-07_DP, 1.79684e-07_DP, 1.61924e-07_DP, &
    1.47382e-07_DP, 1.35691e-07_DP, 1.2646e-07_DP, 1.19303e-07_DP, &
    1.13854e-07_DP, 1.09779e-07_DP, 1.06788e-07_DP, 1.04632e-07_DP, &
    1.03105e-07_DP, 1.02043e-07_DP, 1.01319e-07_DP, 1.00833e-07_DP, &
    1.00513e-07_DP, 1.00306e-07_DP, 1.00175e-07_DP, 1.00093e-07_DP, &
    1.00043e-07_DP, 1.00013e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, &
    1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, &
    1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, &
    1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, 1e-07_DP, &
    1e-07_DP /)

    sb_efficiency = lut(INT((clay_fraction + .5_dp)))

    END FUNCTION sb_efficiency
END SUBROUTINE dust_emissions_KKDU_Astitha1

! ----------------------------------------------------------------------------
SUBROUTINE dust_emissions_DU_Astitha1(kproma, airdens, zust_2d              &
     , slf, ws, cvs, emis_du_cla2, rootdepth                                &
     , dustsrc,lai_in,duflux_ai,duflux_ci, horflux,ustarthr                 &
     , du_nap, du_kp, du_capp, du_mgpp, du_misc                             &
     , du_nap_emflux_ai, du_nap_emflux_ci, du_kp_emflux_ai, du_kp_emflux_ci &
     , du_capp_emflux_ai, du_capp_emflux_ci, du_mgpp_emflux_ai              &
     , du_mgpp_emflux_ci, du_misc_emflux_ai, du_misc_emflux_ci)

! ----------------------------------------------------------------------------
  !
  !     *************** Online calculation of dust fluxes  ********************
  !
  !     Authors:
  !     --------
  !     Written by Marina Astitha for MESSy1.8; The Cyprus Institute-EEWRC
  !                                             September 2009++
  !     Implemented in MESSy2.41 by Mohamed Abdel Kader (2011)
  !
  !     purpose
  !     -------
  !     Emission Fluxes of Mineral Dust depending on soil,friction velocity and
  !     soil moisture. More Information on the emission scheme can be found in
  !     Astitha et al. (2012) in ACPD (DU1 scheme)
  !
  !
  !     INPUT
  !     ---------
  !       airdens        air density                                     [kg/m3]
  !       emis_du_cla2   clay content                                    [%]
  !       du_nap         Na+ fraction in desert soil
  !       du_kp          K+ fraction in desert soil
  !       du_capp        Ca++ fraction in desert soil
  !       du_mgpp        Mg++ fraction in desert soil
  !       du_misc        misc. fraction in desert soil
  !       rootdepth      rooting depth dependent on the land use category[m]
  !       zust_2d        surface friction velocity                       [m/sec]
  !       slf            fractional land-sea mask [0=sea and 1=land]
  !       ws             soil moisture                                   [m]
  !       qm1            specific humidity at 1st model layer            [kg/kg]
  !       qte            specific humidity tendency                  [kg/kg*sec]
  !       tm1            temperature at 1st model layer                  [K]
  !       cvs            fraction of snow cover
  !
  !      OUTPUT
  !       duflux_ai      Vertical flux (emissions) of dust - accumulation
  !                      [kg m-2 s-1]
  !       duflux_ci      Vertical flux (emissions) of dust - coarse
  !                      [kg m-2 s-1]
  !       du_nap_emflux_ai Vertical flux (emissions) of dust Na+, acc. mode
  !                      [kg m-2 s-1]
  !       du_nap_emflux_ci Vertical flux (emissions) of dust Na+, coarse mode
  !                      [kg m-2 s-1]
  !       du_kp_emflux_ai Vertical flux (emissions) of dust K+, acc. mode
  !                      [kg m-2 s-1]
  !       du_kp_emflux_ci Vertical flux (emissions) of dust K+, coarse mode
  !                      [kg m-2 s-1]
  !       du_capp_emflux_ai Vertical flux (emissions) of dust Ca++, acc. mode
  !                      [kg m-2 s-1]
  !       du_capp_emflux_ci Vertical flux (emissions) of dust Ca++, coarse mode
  !                      [kg m-2 s-1]
  !       du_mgpp_emflux_ai Vertical flux (emissions) of dust Mg++, acc. mode
  !                      [kg m-2 s-1]
  !       du_mgpp_emflux_ci Vertical flux (emissions) of dust Mg++, coarse mode
  !                      [kg m-2 s-1]
  !       du_misc_emflux_ai Vertical flux (emissions) of misc. dust, acc. mode
  !                      [kg m-2 s-1]
  !       du_misc_emflux_ci Vertical flux (emissions) of misc. dust, coarse mode
  !                      [kg m-2 s-1]
  !
  !      DIAGNOSTIC OUTPUT
  !       alpha         !sandblasting efficiency [cm-1]
  !       ustarthr      !threshold wind friction velocity [m/s]
  !       baresoil      !fraction of bare soil exposed in a grid cell [-]
  !       horflux       !horizontal flux of dust from soil sources (saltation)
  !                     [kg m-1 s-1]
  !
  ! -----------------------------------------------------------------------
  !
  !    Transport size bins (Perez et al. 2006)
  !    Source Size distribution (d'Almeida 1987;Zender et al. 2003)
  !      RADIUS(min-max)   Effective radius DIAMETER
  !    1) 0.1  -  0.18 um - reff=0.15 um    Dn(um)   Dv(um)   sigma    Mfraction
  !    2) 0.18 -  0.3  um - reff=0.25 um    0.16     0.832     2.10      0.036
  !    3) 0.3  -  0.6  um - reff=0.45 um    1.40     4.820     1.90      0.957
  !    4) 0.6  -  1.0  um - reff=0.78 um   10.00    19.380     1.60      0.007
  !    5) 1.0  -  1.8  um - reff=1.3  um
  !    6) 1.8  -  3.0  um - reff=2.2  um
  !    7) 3.0  -  6.0  um - reff=3.8  um
  !    8) 6.0  - 10.0  um - reff=7.1  um
  ! ----------------------------------------------------------------------------
  USE messy_main_constants_mem,  ONLY: g, rho_H2O
  USE messy_main_tools,          ONLY: ERRFUNC

  !rho_H2O = 999.97_dp    ! density of H2O [kg/m3]
  !g       = 9.80665_dp   ! gravity acceleration [m/s2]

  IMPLICIT NONE

  ! I/O
  INTEGER, PARAMETER     :: NPS=3                !Number of source size bins
  INTEGER, PARAMETER     :: NPT=8                !Number of transport size bins
  INTEGER, INTENT(IN)    :: kproma
  REAL(dp), INTENT(IN)   :: airdens(1:kproma)
  REAL(dp), INTENT(IN)   :: emis_du_cla2(1:kproma)
  REAL(dp), INTENT(IN)   :: lai_in(1:kproma)
  REAL(dp), INTENT(IN)   :: dustsrc(1:kproma,1:nbiomes)
  ! Rooting depth for water and cells with no data is set to -1
  REAL(dp), INTENT(IN)   :: rootdepth(1:kproma)
  REAL(dp), INTENT(IN)   :: zust_2d(1:kproma)
  REAL(dp), INTENT(IN)   :: ws(1:kproma)
  REAL(dp), INTENT(IN)   :: cvs(1:kproma)
  REAL(dp), INTENT(IN)   :: slf(1:kproma)
!
  REAL(dp), INTENT(OUT)  :: horflux(1:kproma)
  REAL(dp), INTENT(INOUT):: ustarthr(1:kproma)
  REAL(dp), INTENT(OUT)  :: duflux_ai(1:kproma),duflux_ci(1:kproma)
  REAL(dp), OPTIONAL, INTENT(IN)   :: du_nap(1:kproma)
  REAL(dp), OPTIONAL, INTENT(IN)   :: du_kp(1:kproma)
  REAL(dp), OPTIONAL, INTENT(IN)   :: du_capp(1:kproma)
  REAL(dp), OPTIONAL, INTENT(IN)   :: du_mgpp(1:kproma)
  REAL(dp), OPTIONAL, INTENT(IN)   :: du_misc(1:kproma)
  REAL(dp), OPTIONAL, INTENT(OUT)  :: du_nap_emflux_ai(1:kproma)
  REAL(dp), OPTIONAL, INTENT(OUT)  :: du_nap_emflux_ci(1:kproma)
  REAL(dp), OPTIONAL, INTENT(OUT)  :: du_kp_emflux_ai(1:kproma)
  REAL(dp), OPTIONAL, INTENT(OUT)  :: du_kp_emflux_ci(1:kproma)
  REAL(dp), OPTIONAL, INTENT(OUT)  :: du_capp_emflux_ai(1:kproma)
  REAL(dp), OPTIONAL, INTENT(OUT)  :: du_capp_emflux_ci(1:kproma)
  REAL(dp), OPTIONAL, INTENT(OUT)  :: du_mgpp_emflux_ai(1:kproma)
  REAL(dp), OPTIONAL, INTENT(OUT)  :: du_mgpp_emflux_ci(1:kproma)
  REAL(dp), OPTIONAL, INTENT(OUT)  :: du_misc_emflux_ai(1:kproma)
  REAL(dp), OPTIONAL, INTENT(OUT)  :: du_misc_emflux_ci(1:kproma)
  REAL(dp)               :: alpha(1:kproma)
  REAL(dP)               :: kappa1(1:kproma),kappa2
  REAL(dp)               :: zos,zo,fdrag,factor(1:kproma),pdiam
  REAL(dp)               :: mij(NPS,NPT),massfraction(NPT),sumtest
!
  REAL(dp)               :: duflux2(1:kproma,NPT)

  REAL(dp)               :: wsgrav(1:kproma)
  REAL(dp)               :: wres(1:kproma)
  REAL(dp)               :: betar(1:kproma)
  REAL(dp)               :: kappa(1:kproma)
  REAL(dp)               :: baresoil(1:kproma)
  REAL                   :: X1,X2,ERFTERM1,ERFTERM2
  ! bulk density of soil particle (kg/m3)
  REAL(dp), PARAMETER    :: bulkdens=2.65E+03_dp
  ! kinematic air viscosity (m2/sec)
  REAL(dp), PARAMETER    :: airvisc=0.157E-04_dp
  ! leaf area index threshold for complete supression of
  ! dust emissions[m2/m2](Mahowald et al. 1999)
  REAL(dp), PARAMETER    :: vlait=0.35_dp
  INTEGER                :: I,J,bio

!*******************************************************************************
!  Olson Biomes used for dust sources: 8,50,51,52
!    --> renamed to 1,2,3,4 respectively (0 for all others)
!
!  Dustmask is 1= biomes  8, 50 and 71  and  0.5= biomes 51, 52
!*******************************************************************************
    REAL, DIMENSION(6) :: dustmask = &              !fraction
       (/0.,1.,1.,0.5,0.5,1./)

!*******************************************************************************
!   8 TRANSPORT SIZE BINS (Perez et al. 2006)
!*******************************************************************************
    REAL, DIMENSION(8) :: mindiam = &              !um
       (/0.2,0.36,0.6,1.2,2.0,3.6,6.0,12.0/)
    REAL, DIMENSION(8) :: maxdiam = &              !um
       (/0.36,0.6,1.2,2.0,3.6,6.0,12.0,20.0/)
    REAL, DIMENSION(8) :: pdens = &                !g/cm-3
       (/2.65,2.65,2.65,2.65,2.65,2.65,2.65,2.65/)

!*******************************************************************************
!   SOIL (SOURCE) SIZE DISTRIBUTION PARAMETERS
!   (d'Almeida 1987,Zender et al. 2003))
!*******************************************************************************
!
    REAL, DIMENSION(3) :: dv = &                   !Mass median Diameter (um)
       (/0.832,4.82,19.38/)
    REAL, DIMENSION(3) :: massfr = &               !Mass Fraction (-)
       (/0.036,0.957,0.007/)
    REAL, DIMENSION(3) :: sigma = &                !Standard deviation
       (/2.1,1.9,1.6/)
!
    ! for the calculation of the standard error function
    INTEGER, PARAMETER :: jint=0

!   INITIALIZATION
    DO J=1,NPT
       duflux2(:,J) = 0._dp
       massfraction(J)= 0._dp
    ENDDO
    baresoil(:)  = 0._dp
    alpha(:)     = 0._dp
    horflux(:)   = 0._dp
    duflux_ai(:) = 0._dp
    duflux_ci(:) = 0._dp

    ! mz_ap_20180214+
    IF (l_ducomp) THEN
       du_nap_emflux_ai(:) = 0._dp
       du_nap_emflux_ci(:) = 0._dp
       du_kp_emflux_ai(:) = 0._dp
       du_kp_emflux_ci(:) = 0._dp
       du_capp_emflux_ai(:) = 0._dp
       du_capp_emflux_ci(:) = 0._dp
       du_mgpp_emflux_ai(:) = 0._dp
       du_mgpp_emflux_ci(:) = 0._dp
       du_misc_emflux_ai(:) = 0._dp
       du_misc_emflux_ci(:) = 0._dp
    ENDIF
!
!   ---------------------------------------------------------------------------
!   CALCULATION OF THE SANDBLASTING EFFICIENCY
!      (a=Vertical Flux/Horizontal Flux) - [cm-1]
!
    WHERE (emis_du_cla2(:) < 20._dp)
       alpha(:)=10**(0.134*emis_du_cla2(:)-6) ! Marticorena and Bergametti (1995)
    ELSEWHERE ((emis_du_cla2(:)>= 20._dp) .AND. (emis_du_cla2(:) < 45._dp))
       alpha(:)=1.e-06_dp                     ! Tegen et al. (2002)
    ELSEWHERE (emis_du_cla2(:)>= 45._dp)
       alpha(:)=1.e-07_dp                     ! Tegen et al. (2002)
    ENDWHERE

!   ----------------------------------------------------------------------------
!   CALCULATION OF THE FRICTION REYNOLDS NUMBER
!   [ B=USTARTHR*DIAMETER/AIRVISC [dimensionless] ]
!
    IF (MAXVAL(ustarthr(:)) <= 0._dp) THEN
       pdiam=60._dp   !Optimum soil size for saltation [um]
       betar(:)=1331._dp*((pdiam*1.e-04)**1.56)+0.38_dp  !Beta=0.835
    ELSE
       pdiam=60._dp   ![um]
       betar(:)= (ustarthr(:)*(pdiam*1.e-06))/airvisc
    ENDIF
!
!   ----------------------------------------------------------------------------
!   CALCULATION OF THE KAPPA COEFFICIENTS THAT ENTER THE FORMULA OF THE
!   THRESHOLD FRICTION VELOCITY USTARTHR
!   (Marticorena and Bergametti (1995);Marticorena et al. (1997))
!
    kappa1(:)=SQRT((pdens(1)*g*pdiam*1.e-06)/(airdens(:)*1.e-03)) ![m/s]
    kappa2=SQRT(1.+(0.006/(pdens(1)*g*100.*(pdiam*1.e-04)**2.5)))![dimensionless]
    kappa(:)=kappa1(:)*kappa2                                     ![m/s]
!
!   CALCULATION OF THE THRESHOLD FRICTION VELOCITY
!
    WHERE((betar(:) > 0.03_dp).AND.(betar(:) < 10._dp))          ! (0.03<B<10)

       ustarthr(:) = 0.129_dp*kappa(:)/SQRT((1.928_dp*betar(:)**0.092)-1.)
       ! [m/s]
    ELSEWHERE (betar(:) >= 10._dp)                               ! (B>10)
       ustarthr(:) = 0.129_dp*kappa(:)*&
            (1.-0.0858_dp*EXP(-0.0617_dp*(betar(:)-10.))) ![m/s]
    ENDWHERE
!
!  -----------------------------------------------------------------------------
!   DRAG PARTITION SCHEME TO QUANTIFY THE FRACTION OF THE TOTAL WIND SHEAR
!   STRESS ON THE ERODIBLE SURFACES TO MOBILIZE THE SOIL PARTICLES
!
!     fdrag = Ratio of local to total friction velocity
!
    zos = 0.00333_dp   !Local roughness length of the uncovered surface
                       !(smooth roughness length) [cm]
!                       (Zender etal 03;MaB95 suggest 1.e-03cm)
!
    zo  = 0.01_dp   !Surface aeolian Roughness Length  [cm]   (Zender etal 2003)
!
    fdrag=1.-(LOG(zo/zos)/LOG(0.35_dp*((10/zos)**0.8)))

    ustarthr(:)=ustarthr(:)/fdrag    ![m/s]
!
!   -----------------------------------------------------------------------------
!   INCREASE OF THE THRESHOLD VELOCITY DUE TO SOIL MOISTURE FOLLOWING
!    Fecan et al. (1999)
!
!   Soil residual moisture WRES (%(mass of water/mass of dry soil))
!   Soil Moisture WS from ECHAM is in meters --> divide by rooting depth to
!   calculate the gravimetric soil moisture (Hillel (1980))
!
    wres(:)=0.0014_dp*(emis_du_cla2(:)**2) + 0.17_dp*emis_du_cla2(:)
!
    WHERE (rootdepth(:) > 0.)
       !Rooting depth for water and cells with no data is set to -1
       !Hillel (1980), Environment. Soil Physics
       wsgrav(:)=100.0_dp*(ws(:)*rho_H2O)/(rootdepth(:)*bulkdens)
    ELSEWHERE
       wsgrav(:)= 0.
    ENDWHERE

    WHERE (wsgrav(:) >= wres(:))
       !Wet threshold friction velocity
       ustarthr(:)=ustarthr(:)*SQRT(1.+1.21_dp*(wsgrav(:) - wres(:))**0.68)
    ELSEWHERE
       ustarthr(:)=ustarthr(:)*1.0_dp
    ENDWHERE
!
!  -----------------------------------------------------------------------------
!  CALCULATION OF THE HORIZONTAL DUST EMISSIONS FLUX [HORFLUX(:)]
!  (SALTATION IS INITIATED WHEN U* > U*thresh)
!
!   The value 2.61 comes from the wind tunnel tests from White (1979)
!   Substitute 2.61 with 1 after Darmenova et al. (2009)  - M.Astitha 14sep2011
!
   WHERE (zust_2d(:) > ustarthr(:))
      factor(:)= (airdens(:)/g)*(zust_2d(:)**3)    ![kg/m*sec]
      horflux(:)= factor(:)*(1+(ustarthr(:)/zust_2d(:))) &
         *(1-(ustarthr(:)**2/zust_2d(:)**2))  !Marticorena and Bergametti (1995)
   ELSEWHERE
      horflux(:)=0.
   ENDWHERE
!
!  -----------------------------------------------------------------------------
!  CALCULATION OF THE MASSFRACTION (Source Sizes To Transport Sizes)
!  [MASSFRACTION(J)]
!
!  Note: The array mij(I,J) is the fraction of the size distribution I in the
!        source groups that overlaps the transport bin J.
!        In that way we compute the mapping from the lognormal source
!        distribution to the transport size bins.
!
!
   sumtest=0.0
   DO I=1,NPS  !Loop over the source modes
      DO J=1,NPT !Loop over the transport bins
!
         X1=REAL((LOG((maxdiam(J))/dv(I)))/(SQRT(2.)*LOG(sigma(I))))
         X2=REAL((LOG((mindiam(J))/dv(I)))/(SQRT(2.)*LOG(sigma(I))))

         CALL ERRFUNC(X1,ERFTERM1,JINT) ! Calculate standard error function
         CALL ERRFUNC(X2,ERFTERM2,JINT)
!
         mij(I,J)=0.5 * (ERFTERM1 - ERFTERM2) !Eq. 12 of Zender et al. 2003
         mij(I,J)=mij(I,J) * massfr(I)

         X1 = 0.0
         X2 = 0.0
         ERFTERM1 = 0.0
         ERFTERM2 = 0.0

         sumtest=sumtest+mij(I,J)

      ENDDO
   ENDDO

   IF(sumtest.lt.0.9.or.sumtest.gt.1.) THEN
      WRITE(*,*)'Problem with the total mass overlap= ', sumtest
   ENDIF

   DO J=1,NPT
      DO I=1,NPS
         massfraction(J)= massfraction(J) + mij(I,J)
      ENDDO
   ENDDO
   !
   !  --------------------------------------------------------------------------
   !  CALCULATION OF THE VERTICAL FLUX OF DUST EMISSIONS [DUFLUX2(:)]
   !

   DO bio=1,nbiomes
      baresoil(:) = baresoil(:) + dustsrc(:,bio)*dustmask(bio)
   ENDDO

   WHERE((slf(:) > 0.5_dp).and.(cvs(:).eq.0._dp))
      baresoil(:) = &
           baresoil(:) * slf(:) * (1._dp-min(1._dp,min(lai_in(:),vlait))/vlait)
   ELSEWHERE
      baresoil(:)=0._dp
   ENDWHERE
!
    DO J=1,NPT      !Number of transport size bins
       duflux2(:,J) = 1.0E-04 * (alpha(:)*100.) * baresoil(:) &
            * horflux(:) * massfraction(J) ![kg/m2*sec]
    ENDDO           !NPT
!
!
    !Accumulation insoluble mode of dust
    DO J=1,3
       duflux_ai(:)=duflux_ai(:)+duflux2(:,J)
    ENDDO
!
    !Coarse insoluble mode of dust
    DO J=4,8
       duflux_ci(:)=duflux_ci(:)+duflux2(:,J)
    ENDDO
!
    IF (l_ducomp) THEN
       du_nap_emflux_ai(:) = duflux_ai(:) * du_nap(:)
       du_nap_emflux_ci(:) = duflux_ci(:) * du_nap(:)
       du_kp_emflux_ai(:) = duflux_ai(:) * du_kp(:)
       du_kp_emflux_ci(:) = duflux_ci(:) * du_kp(:)
       du_capp_emflux_ai(:) = duflux_ai(:) * du_capp(:)
       du_capp_emflux_ci(:) = duflux_ci(:) * du_capp(:)
       du_mgpp_emflux_ai(:) = duflux_ai(:) * du_mgpp(:)
       du_mgpp_emflux_ci(:) = duflux_ci(:) * du_mgpp(:)
       du_misc_emflux_ai(:) = duflux_ai(:) * du_misc(:)
       du_misc_emflux_ci(:) = duflux_ci(:) * du_misc(:)
    ENDIF

  END SUBROUTINE dust_emissions_DU_Astitha1

!-------------------------------------------------------------------------------
SUBROUTINE dust_emissions_DU_Astitha2(kproma, airdens, zust_2d              &
     , slf, ws, cvs, emis_du_cla2                                           &
     , rootdepth, lai_in,dustsrc, soiltext, duflux_ai, duflux_ci , ustarthr &
     , du_nap, du_kp, du_capp, du_mgpp, du_misc                             &
     , du_nap_emflux_ai, du_nap_emflux_ci, du_kp_emflux_ai, du_kp_emflux_ci &
     , du_capp_emflux_ai, du_capp_emflux_ci, du_mgpp_emflux_ai              &
     , du_mgpp_emflux_ci, du_misc_emflux_ai, du_misc_emflux_ci )

  ! --------------------------------------------------------------------------
  !
  !     *************** Online calculation of dust fluxes  ********************
  !
  !     Authors:
  !     --------
  !     Marina Astitha ; The Cyprus Institute-EEWRC ; September 2011++
  !     Implemented in MESSy2.41 by M. Abdel Kader(EEWRC, CyI) (2011)
  !
  !     Purpose
  !     -------
  !     Emission Fluxes of Mineral Dust depending on soil,friction velocity
  !     and soil moisture: More Information on the emission scheme can be found
  !     in Astitha et al. (2012) in ACPD (DU2 scheme)
  !
  !     INPUT
  !     ---------
  !     time_step_len  lenght of the timestep                            [sec]
  !     airdens        air density                                       [kg/m3]
  !     emis_du_cla2   clay content                                      [%]
  !       du_nap         Na+ fraction in desert soil
  !       du_kp          K+ fraction in desert soil
  !       du_capp        Ca++ fraction in desert soil
  !       du_mgpp        Mg++ fraction in desert soil
  !       du_misc        misc. fraction in desert soil
  !     rootdepth      rooting depth dependent on the land use category  [m]
  !     zust_2d        surface friction velocity                         [m/sec]
  !     lai_in         leaf area index                                   [m2/m2]
  !     slf            fractional land-sea mask [0=sea and 1=land]
  !     ws             soil moisture                                     [m]
  !     qm1            specific humidity at 1st model layer              [kg/kg]
  !     qte            specific humidity tendency                    [kg/kg*sec]
  !     tm1            temperature at 1st model layer                    [K]
  !     cvs            fraction of snow cover
  !     soiltext       Percentage of Soil texture classes from Zobler in
  !                    each grid cell (0-9)  [percentage]
  !     dustsrc        Olson world ecosystem biomes for desert areas
  !                                                   [percentage of each biome]
  !
  !
  !    OUTPUT
  !     duflux_ai  Vertical flux (emissions) of dust - accumulation  [kg m-2 s-1]
  !     duflux_ci  Vertical flux (emissions) of dust - coarse        [kg m-2 s-1]
  !       du_nap_emflux_ai Vertical flux (emissions) of dust Na+, acc. mode
  !                      [kg m-2 s-1]
  !       du_nap_emflux_ci Vertical flux (emissions) of dust Na+, coarse mode
  !                      [kg m-2 s-1]
  !       du_kp_emflux_ai Vertical flux (emissions) of dust K+, acc. mode
  !                      [kg m-2 s-1]
  !       du_kp_emflux_ci Vertical flux (emissions) of dust K+, coarse mode
  !                      [kg m-2 s-1]
  !       du_capp_emflux_ai Vertical flux (emissions) of dust Ca++, acc. mode
  !                      [kg m-2 s-1]
  !       du_capp_emflux_ci Vertical flux (emissions) of dust Ca++, coarse mode
  !                      [kg m-2 s-1]
  !       du_mgpp_emflux_ai Vertical flux (emissions) of dust Mg++, acc. mode
  !                      [kg m-2 s-1]
  !       du_mgpp_emflux_ci Vertical flux (emissions) of dust Mg++, coarse mode
  !                      [kg m-2 s-1]
  !       du_misc_emflux_ai Vertical flux (emissions) of misc. dust, acc. mode
  !                      [kg m-2 s-1]
  !       du_misc_emflux_ci Vertical flux (emissions) of misc. dust, coarse mode
  !                      [kg m-2 s-1]
  !
  !    DIAGNOSTIC OUTPUT
  !     alpha         !sandblasting efficiency [cm-1]
  !     ustarthr      !threshold wind friction velocity [m/s]
  !     baresoil      !fraction of bare soil exposed in a grid cell [-]
  !     horflux       !horizontal flux of dust from soil sources (saltation)
  !                   [kg m-1 s-1]
  !
  ! -----------------------------------------------------------------------
  !
  !  Transport size bins (Perez et al. 2006)
  !  Source Size distribution (d'Almeida 1987;Zender et al. 2003)
  !    RADIUS(min-max)   Effective radius     DIAMETER
  !  1) 0.1  -  0.18 um - reff=0.15 um        Dn(um)   Dv(um)   sigma  Mfraction
  !  2) 0.18 -  0.3  um - reff=0.25 um         0.16     0.832   2.10    0.036
  !  3) 0.3  -  0.6  um - reff=0.45 um         1.40     4.820   1.90    0.957
  !  4) 0.6  -  1.0  um - reff=0.78 um        10.00    19.380   1.60    0.007
  !  5) 1.0  -  1.8  um - reff=1.3  um
  !  6) 1.8  -  3.0  um - reff=2.2  um
  !  7) 3.0  -  6.0  um - reff=3.8  um
  !  8) 6.0  - 10.0  um - reff=7.1  um
  ! ----------------------------------------------------------------------------
  USE messy_main_constants_mem,  ONLY: g, rho_H2O
  USE messy_main_tools,          ONLY: ERRFUNC

    !rho_H2O = 999.97_dp    ! density of H2O [kg/m3]
    !g       = 9.80665_dp   ! gravity acceleration [m/s2]

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)    :: kproma
    REAL(dp), INTENT(IN)   :: airdens(1:kproma) ! um_ak_20120625
    ! um_ak_20120625 REAL(dp), INTENT(IN)   :: time_step_len
    REAL(dp), INTENT(IN)   :: emis_du_cla2(1:kproma)
    REAL(dp), INTENT(IN)   :: lai_in(1:kproma)
    REAL(dp), INTENT(IN)   :: dustsrc(1:kproma,1:nbiomes)
    REAL(dp), INTENT(IN)   :: soiltext(1:kproma,1:nclasses)
    ! Rooting depth for water and cells with no data is set to -1
    REAL(dp), INTENT(IN)   :: rootdepth(1:kproma)
    REAL(dp), INTENT(IN)   :: zust_2d(1:kproma)
    REAL(dp), INTENT(IN)   :: ws(1:kproma)
    REAL(dp), INTENT(IN)   :: cvs(1:kproma)
    REAL(dp), INTENT(IN)   :: slf(1:kproma)
    REAL(dp), INTENT(INOUT):: ustarthr(1:kproma,1:maxsizes) ! um_ak_20120625
    REAL(dp), INTENT(OUT)  :: duflux_ai(1:kproma),duflux_ci(1:kproma)

    REAL(dp), OPTIONAL, INTENT(IN)   :: du_nap(1:kproma)
    REAL(dp), OPTIONAL, INTENT(IN)   :: du_kp(1:kproma)
    REAL(dp), OPTIONAL, INTENT(IN)   :: du_capp(1:kproma)
    REAL(dp), OPTIONAL, INTENT(IN)   :: du_mgpp(1:kproma)
    REAL(dp), OPTIONAL, INTENT(IN)   :: du_misc(1:kproma)

    REAL(dp), OPTIONAL, INTENT(OUT)  :: du_nap_emflux_ai(1:kproma)
    REAL(dp), OPTIONAL, INTENT(OUT)  :: du_nap_emflux_ci(1:kproma)
    REAL(dp), OPTIONAL, INTENT(OUT)  :: du_kp_emflux_ai(1:kproma)
    REAL(dp), OPTIONAL, INTENT(OUT)  :: du_kp_emflux_ci(1:kproma)
    REAL(dp), OPTIONAL, INTENT(OUT)  :: du_capp_emflux_ai(1:kproma)
    REAL(dp), OPTIONAL, INTENT(OUT)  :: du_capp_emflux_ci(1:kproma)
    REAL(dp), OPTIONAL, INTENT(OUT)  :: du_mgpp_emflux_ai(1:kproma)
    REAL(dp), OPTIONAL, INTENT(OUT)  :: du_mgpp_emflux_ci(1:kproma)
    REAL(dp), OPTIONAL, INTENT(OUT)  :: du_misc_emflux_ai(1:kproma)
    REAL(dp), OPTIONAL, INTENT(OUT)  :: du_misc_emflux_ci(1:kproma)
!
    REAL(dp), DIMENSION(:),     ALLOCATABLE :: dS,Dn,Dsur,Stotal,kappa2
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: duflux2,total,kappa1,betar
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: horflux, horflux1,kappa,soilscale
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: mijA
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: mij
!
    REAL(dp)               :: zos,zo,fdrag,factor(1:kproma)
    REAL(dp)               :: massfraction(1:kproma,1:DU_A2_npt)
    REAL(dp)               :: fact,logstep
    REAL(dp)               :: wres(1:kproma),sumtest(1:kproma)
    REAL(dp)               :: horflux_tot(1:kproma)
    REAL(dp)               :: alpha(1:kproma)
    REAL(dp)               :: wsgrav(1:kproma)
    REAL(dp)               :: baresoil(1:kproma)
    REAL                   :: X1,X2,ERFTERM1,ERFTERM2
    ! bulk density of soil particle (kg/m3)
    REAL(dp), PARAMETER    :: bulkdens=2.65E+03_dp
    ! kinematic air viscosity (m2/sec)
    REAL(dp), PARAMETER    :: airvisc=0.157E-04_dp
    ! leaf area index threshold for complete supression of
    ! dust emissions [m2/m2](Mahowald et al. 1999)
    REAL(dp), PARAMETER    :: vlait=0.35_dp
    INTEGER                :: I,J,bio,isrc,nc,kk,kmax
    REAL(dp)               :: D
!
!*******************************************************************************
!  [dustsrc]:Olson Biomes used for dust sources: 0,8,50,51,52,71
!    --> renamed to 1,2,3,4,5,6 respectively
!
!  Dustmask is  1  = biomes  8, 50 and 71
!       and     0.5= biomes 51, 52
!       and     0= for water (1st value)
!
    REAL, DIMENSION(nbiomes) :: dustmask = &              !fraction
       (/0.,1.,1.,0.5,0.5,1./)
!*******************************************************************************
!   Zobler Soil Texture Classes (soiltext(i,j,nclasses)=percentage of
!      soil texture in each grid cell):
!
!   1=Water, 2=Coarse, 3=Medium, 4=Fine, 5=coarse-medium, 6=coarse-fine
!   7=medium-fine,  8=coarse-medium-fine, 9=organic, 10=land ice
!
!*******************************************************************************
!*******************************************************************************
!   8 TRANSPORT SIZE BINS (Perez et al. 2006)
!*******************************************************************************
    REAL, DIMENSION(8) :: mindiam = &              !um
       (/0.2,0.36,0.6,1.2,2.0,3.6,6.0,12.0/)
    REAL, DIMENSION(8) :: maxdiam = &              !um
       (/0.36,0.6,1.2,2.0,3.6,6.0,12.0,20.0/)
    REAL, DIMENSION(8) :: pdens = &                !g/cm-3
       (/2.65,2.65,2.65,2.65,2.65,2.65,2.65,2.65/)
!
!*******************************************************************************
!   SOIL SIZE DISTRIBUTION PARAMETERS (Zobler - Tegen)
!*******************************************************************************
!
    REAL, DIMENSION(DU_A2_nps)  :: Dm = &     !mass median diameter (um)
       (/707.,158.,15.,2.0/)
    REAL,DIMENSION(DU_A2_nps)  :: sigma = &
       (/2.,2.,2.,2./)

    !mass fraction for soil mode 1 (for the 10 soil types)
   REAL, DIMENSION(nclasses)  :: massfr1 = &
    (/0.00,0.43, 0.00,0.00,0.10,0.00,0.00,0.23,0.25,0.25/)
   !mass fraction for soil mode 2 (for the 10 soil types)
   REAL, DIMENSION(nclasses)  :: massfr2 = &
    (/0.00,0.40,0.37,0.00,0.50,0.50,0.27,0.23,0.25,0.25/)
   !mass fraction for soil mode 3 (for the 10 soil types)
   REAL, DIMENSION(nclasses)  :: massfr3 = &
    (/0.00,0.17,0.33,0.33,0.20,0.12,0.25,0.19,0.25,0.25/)
   !mass fraction for soil mode 4 (for the 10 soil types)
   REAL, DIMENSION(nclasses)  :: massfr4 = &
    (/0.00,0.00,0.30,0.67,0.20,0.38,0.48,0.35,0.25,0.25/)

!   sandblasting efficiency [cm-1] for each soil type from Tegen
!   REAL, DIMENSION(nclasses)  :: alpha1 = &
!    (/0.00,2.1E-06,4.0E-06,1.E-07,2.7E-06,2.8E-06,1.E-07,2.5E-06,0.00,0.00/)
!
   ! for the calculation of the standard error function
   INTEGER, PARAMETER :: jint=0

   ALLOCATE(dS(1:maxsizes),Dn(1:DU_A2_nps),Dsur(1:DU_A2_nps),Stotal(1:DU_A2_nps)&
        ,kappa2(1:maxsizes))
   ALLOCATE(duflux2(1:kproma,1:DU_A2_npt),total(1:kproma,1:DU_A2_nps)           &
        ,kappa1(1:kproma,1:maxsizes),betar(1:kproma,1:maxsizes))
   ALLOCATE(horflux(1:kproma,1:DU_A2_nps), horflux1(1:kproma,1:maxsizes)     &
        ,kappa(1:kproma,1:maxsizes),soilscale(1:kproma,1:DU_A2_nps)          &
        ,mijA(1:DU_A2_nps,1:DU_A2_npt))
   ALLOCATE(mij(1:kproma,1:DU_A2_nps,1:DU_A2_npt))

!   INITIALIZATION
    DO J=1,DU_A2_npt
       duflux2(:,J)     = 0._dp
       massfraction(:,J)= 0._dp
    ENDDO
    DO isrc=1,DU_A2_nps
       horflux(:,isrc)  = 0._dp
       soilscale(:,isrc)= 0._dp
       total(:,isrc)    = 0._dp
    ENDDO
    horflux1(:,:)    = 0._dp
    betar(:,:)       = 0._dp
    ustarthr(:,:)    = 0._dp

    horflux_tot(:)   = 0._dp
    baresoil(:)      = 0._dp
    alpha(:)         = 0._dp
    duflux_ai(:)     = 0._dp
    duflux_ci(:)     = 0._dp

    IF (l_ducomp) THEN
       du_nap_emflux_ai(:) = 0._dp
       du_nap_emflux_ci(:) = 0._dp
       du_kp_emflux_ai(:) = 0._dp
       du_kp_emflux_ci(:) = 0._dp
       du_capp_emflux_ai(:) = 0._dp
       du_capp_emflux_ci(:) = 0._dp
       du_mgpp_emflux_ai(:) = 0._dp
       du_mgpp_emflux_ci(:) = 0._dp
       du_misc_emflux_ai(:) = 0._dp
       du_misc_emflux_ci(:) = 0._dp
    ENDIF

    DO I=1,DU_A2_nps
       DO J=1,DU_A2_npt
          mijA(I,J)        =0.
          mij(:,I,J)       =0.
       ENDDO
    ENDDO

!   ----------------------------------------------------------------------------
!   CALCULATION OF THE SANDBLASTING EFFICIENCY
!    (a=Vertical Flux/Horizontal Flux) - [cm-1]
!
    WHERE (emis_du_cla2(:) < 20._dp)
       alpha(:)=10**(0.134*emis_du_cla2(:)-6) ! Marticorena and Bergametti (1995)
    ELSEWHERE ((emis_du_cla2(:)>= 20._dp).AND.(emis_du_cla2(:) < 45._dp))
       alpha(:)=1.e-06_dp                     ! Tegen et al. (2002)
    ELSEWHERE (emis_du_cla2(:)>= 45._dp)
       alpha(:)=1.e-07_dp                     ! Tegen et al. (2002)
    ENDWHERE

!   ----------------------------------------------------------------------------
!   CALCULATION OF THE FRICTION REYNOLDS NUMBER
!   [ B=USTARTHR*DIAMETER/AIRVISC [dimensionless] ]
!
!
    ! logstep = (log(Dmax)-log(D1))/maxsizes !Dln(Dp)
    logstep = (log(1000.)-log(0.1))/201
!
    ! [kg/m*sec] The value 2.61 is replaced by 1. after Darmenova et al. (2009)
    factor(:)= 1.0_dp*(airdens(:)/g)*(zust_2d(:)**3)
    ! [kg/m*sec] The value 2.61 comes from the wind tunnel tests from White(1979)
!   factor(:)= 2.61_dp*(airdens(:)/g)*(zust_2d(:)**3)
!
    DO isrc=1,DU_A2_nps
       ! Number median diameter (167,37,3.5,0.47)
       Dn(isrc) = exp(log(Dm(isrc))-3*log(sigma(isrc))**2)
       ! Surface median diameter (437,97.7,9.3,1.2)
       Dsur(isrc) = exp(log(Dn(isrc))+2*log(sigma(isrc))**2)
!
       Stotal(isrc)=0._dp
!
       D = 0.00001   !Dmin (cm) =0.1um
       kk = 0
       DO WHILE (D.LE.0.1+1.E-05_dp)      !diameter (cm) Dmax=0.1cm=1000um
          kk = kk + 1

          betar(:,kk)=1331._dp*(D**1.56)+0.38_dp

          !   ------------------------------------------------------------------
          !   ------------------------------------------------------------------
          !   CALCULATION OF THE KAPPA COEFFICIENTS THAT ENTER THE FORMULA OF THE
          !   THRESHOLD FRICTION VELOCITY USTARTHR
          !   (Marticorena and Bergametti (1995);Marticorena et al. (1997))
!
          kappa1(:,kk)=SQRT((pdens(1)*g*(D*1.e-02))/(airdens(:)*1.e-03))  ![m/s]
          kappa2(kk)=SQRT(1.+(0.006/(pdens(1)*g*100.*(D**2.5))))
          ![dimensionless] (0.006 g cm0.5 s2)
          kappa(:,kk)=kappa1(:,kk)*kappa2(kk)                             ![m/s]
!
          !   CALCULATION OF THE THRESHOLD FRICTION VELOCITY
          !
          WHERE((betar(:,kk) > 0.03_dp).AND.(betar(:,kk) < 10._dp)) !(0.03<B<10)

             ustarthr(:,kk) = &
             0.129_dp*kappa(:,kk)/SQRT((1.928_dp*(betar(:,kk)**0.092))-1.)![m/s]

          ELSEWHERE (betar(:,kk) >= 10._dp)                             ! (B>10)

             ustarthr(:,kk) = 0.129_dp*kappa(:,kk) &
                  *(1.-0.0858_dp*EXP(-0.0617_dp*(betar(:,kk)-10.))) ![m/s]

          ENDWHERE
!
!  -----------------------------------------------------------------------------
!   DRAG PARTITION SCHEME TO QUANTIFY THE FRACTION OF THE TOTAL WIND SHEAR
!   STRESS ON THE ERODIBLE SURFACES TO MOBILIZE THE SOIL PARTICLES
!
!     fdrag = Ratio of local to total friction velocity
!
     zos = 0.00333_dp   !Local roughness length of the uncovered surface
                        !(smooth roughness length) [cm]
                        ! (Zender etal 03;MaB95 suggest 1.e-03cm)
!
     zo  = 0.01_dp   !Surface aeolian Roughness Length  [cm]   (Zender etal 2003)
!
     fdrag=1.-(LOG(zo/zos)/LOG(0.35_dp*((10/zos)**0.8)))

     ustarthr(:,kk)=ustarthr(:,kk)/fdrag    ![m/s]
!
!   ----------------------------------------------------------------------------
!   INCREASE OF THE THRESHOLD VELOCITY DUE TO SOIL MOISTURE FOLLOWING
!    Fecan et al. (1999)
!
!   Soil residual moisture WRES (%(mass of water/mass of dry soil))
!   Soil Moisture WS from ECHAM is in meters --> divide by rooting depth to
!    calculate the gravimetric soil moisture (Hillel (1980))
!
      wres(:)=0.0014_dp*(emis_du_cla2(:)**2) + 0.17_dp*emis_du_cla2(:)
!
      WHERE (rootdepth(:) > 0.)
         ! Rooting depth for water and cells with no data is set to -1
         wsgrav(:)=100.0_dp*(ws(:)*rho_H2O)/(rootdepth(:)*bulkdens)
         !Hillel (1980), Environment. Soil Physics
      ELSEWHERE
         wsgrav(:)= 0.
      ENDWHERE

      WHERE (wsgrav(:) >= wres(:))
         ! Wet threshold friction velocity
         ustarthr(:,kk)=ustarthr(:,kk) &
              *SQRT(1.+1.21_dp*(wsgrav(:) - wres(:))**0.68)
      ELSEWHERE
         ustarthr(:,kk)=ustarthr(:,kk)*1.0_dp
      ENDWHERE
!
!  -----------------------------------------------------------------------------
!  CALCULATION OF THE HORIZONTAL DUST EMISSIONS FLUX [HORFLUX(:)]
!  (SALTATION IS INITIATED WHEN U* > U*thresh)
!
      fact = 1./(sqrt(2*pi)*log(sigma(isrc)))
      ! Surface area covered from particles of size D
      dS(kk) = fact*exp(-(log(D*1.e+04)-log(Dsur(isrc)))**2 &
           /(2*log(sigma(isrc))**2))*logstep
      ! Total surface covered by particles of each source mode
      Stotal(isrc)=Stotal(isrc)+dS(kk)
      D = D * exp(0.04605)  !cm
   ENDDO !D (diameter)
   kmax=kk
!
!
!
   DO kk=1,kmax

      WHERE (zust_2d(:) > ustarthr(:,kk))
         horflux1(:,kk)= factor(:)*(1+(ustarthr(:,kk)/zust_2d(:)))&
              *(1-(ustarthr(:,kk)**2/zust_2d(:)**2))*(dS(kk)/Stotal(isrc))
      ENDWHERE
      ! total horiz flux for 1 source size
      total(:,isrc)=total(:,isrc)+horflux1(:,kk)
!
     ENDDO  !kk loop
!
      horflux(:,isrc)=total(:,isrc)

!
   ENDDO !isrc(source modes)
!
!  -----------------------------------------------------------------------------
!  CALCULATION OF THE MASSFRACTION (Source Sizes To Transport Sizes)
!  [MASSFRACTION(J)]
!
!  Note: The array mij(I,J) is the fraction of the size distribution I in the
!        source groups that overlaps the transport bin J.
!        In that way we compute the mapping from the lognormal source
!        distribution to the transport size bins.
!
   sumtest(:)=0.0
!
   DO I=1,DU_A2_nps  !Loop over the source modes (1..4)
      DO NC=1,nclasses  !Source types (1..10)
         IF (I.EQ.1) THEN
            soilscale(:,I)=soilscale(:,I)+soiltext(:,NC)*massfr1(NC)
         ELSEIF(I.EQ.2) THEN
            soilscale(:,I)=soilscale(:,I)+soiltext(:,NC)*massfr2(NC)
         ELSEIF(I.EQ.3) THEN
            soilscale(:,I)=soilscale(:,I)+soiltext(:,NC)*massfr3(NC)
         ELSEIF(I.EQ.4) THEN
            soilscale(:,I)=soilscale(:,I)+soiltext(:,NC)*massfr4(NC)
         ENDIF
      ENDDO
   ENDDO

   DO I=1,DU_A2_nps  !Loop over the source modes (1..4)
      DO J=1,DU_A2_npt !Loop over the transport bins (1..8)

         X1=REAL((LOG((maxdiam(J))/Dm(I)))/(SQRT(2.)*LOG(sigma(I))))
         X2=REAL((LOG((mindiam(J))/Dm(I)))/(SQRT(2.)*LOG(sigma(I))))

         CALL ERRFUNC(X1,ERFTERM1,JINT) ! Calculate standard error function
         CALL ERRFUNC(X2,ERFTERM2,JINT)

         mijA(I,J) = 0.5 * (ERFTERM1 - ERFTERM2) !Eq. 12 of Zender et al. 2003
         mij(:,I,J)= mijA(I,J) * soilscale(:,I)

         X1 = 0.0
         X2 = 0.0
         ERFTERM1 = 0.0
         ERFTERM2 = 0.0

         sumtest(:)=sumtest(:)+mij(:,I,J)

         massfraction(:,J)= massfraction(:,J) + mij(:,I,J)

      ENDDO    !DU_A2_npt
   ENDDO     !DU_A2_nps

   DO i=1,kproma
      IF(sumtest(i).gt.1.)THEN
         WRITE(*,*)'The total mass overlap is greater than 1 = ', sumtest(i)
      ENDIF
   ENDDO


   !  --------------------------------------------------------------------------
   !  CALCULATION OF THE VERTICAL FLUX OF DUST EMISSIONS [DUFLUX2(:)]
   !
   DO bio=1,nbiomes
      baresoil(:) = baresoil(:) + dustsrc(:,bio)*dustmask(bio)
   ENDDO
!
!  cvs=snow cover, slf=sea-land fraction
   WHERE((slf(:) > 0.5_dp).and.(cvs(:).eq.0._dp))
      baresoil(:)=baresoil(:) * slf(:) &
           * (1._dp-min(1._dp,min(lai_in(:),vlait))/vlait)
   ELSEWHERE
      baresoil(:)=0._dp
   ENDWHERE
!

   DO isrc=1,DU_A2_nps    !Source sizes (1..4)
      horflux_tot(:)=horflux_tot(:)+horflux(:,isrc)*soilscale(:,isrc)
   ENDDO

!
   DO J=1,DU_A2_npt      !Number of transport size bins
    duflux2(:,J) =  1.E-03 *baresoil(:) &
         * (alpha(:)*100) * horflux_tot(:) * massfraction(:,J) ! [kg/m2*sec]
   ENDDO           !DU_A2_npt
!
    DO J=1,3
     duflux_ai(:)=duflux_ai(:)+duflux2(:,J)  !Accumulation insoluble mode of dust
    ENDDO
!
    DO J=4,8
     duflux_ci(:)=duflux_ci(:)+duflux2(:,J)  !Coarse insoluble mode of dust
    ENDDO
!
    IF (l_ducomp) THEN
        du_nap_emflux_ai(:)  = duflux_ai(:) * du_nap(:)
        du_nap_emflux_ci(:)  = duflux_ci(:) * du_nap(:)
        du_kp_emflux_ai(:)   = duflux_ai(:) * du_kp(:)
        du_kp_emflux_ci(:)   = duflux_ci(:) * du_kp(:)
        du_capp_emflux_ai(:) = duflux_ai(:) * du_capp(:)
        du_capp_emflux_ci(:) = duflux_ci(:) * du_capp(:)
        du_mgpp_emflux_ai(:) = duflux_ai(:) * du_mgpp(:)
        du_mgpp_emflux_ci(:) = duflux_ci(:) * du_mgpp(:)
        du_misc_emflux_ai(:) = duflux_ai(:) * du_misc(:)
        du_misc_emflux_ci(:) = duflux_ci(:) * du_misc(:)
    ENDIF

    DEALLOCATE(dS,Dn,Dsur,Stotal,kappa2)
    DEALLOCATE(duflux2,total,kappa1,betar)
    DEALLOCATE(horflux, horflux1,kappa,soilscale,mijA)
    DEALLOCATE(mij)

END SUBROUTINE dust_emissions_DU_Astitha2

! -------------------------------------------------------------------------
SUBROUTINE NOemis_yl95sl10_pulsing( kproma                                    &
                                , timestep_length, current_timestep, ls_prec  &
                                , conv_prec, ndaylength                       &
                                , rpulseregime, pulseday, pulse, prec_hist)
! -------------------------------------------------------------------------
!  calculate the pulsing factor for NO soil emission, that occurs after a
!  certain period of dryness. The precipitation is accumulated for one day
!  and stored fore the last 14 days. At the first time step of a new day
!  it is checked, if a new pulse will be activated, according to the
!  precipitation history of the last 14 days and the amount of precipitation
!  during the last 24 hours. Also at the first timestep of a new day, the
!  precipitation history is shifted by one day.
!  Each timestep the pulsing factor is decremented by the formula
!  pulse(time) = a * EXP(b * time), where time is in days after the
!  initialisation. Also each timestep the precipitation of that day is
!  accumulated.
!
! Interface:
! ----------
!
! Input for time control:
!  - timestep_length      [s]
!  - current_timestep     []
!  - ndaylength           [s]
!
! Input for precipitation history:
!  - ls_prec              [m]       large scale precipitation
!  - conv_prec            [m]       convective precipitation
!
! Output:
!  - rpulseregime         []        Index (0: no pulse, 1: sprinkle, 2: shower
!                                        , 3: heavy rain)
!  - pulseday             [day]     duration of pulse, starting at 1
!  - pulse                []        pulsing factor
!  - prec_hist            [m]       precipitation history
!
! Local variables:
!  - pulseduration        [day]     duration of pulse, depending on the pulsing
!                                   regime
!  - pulsefact            []        "a" value of the formula above
!  - pulseexp             [day^-1]  "b" value
!  - pulselimit           [m]       lower value of last day's precipitation,
!                                   that initializes
!                                   the corresponding pulsing regime
!  - pulseinitlimit       [m]       maximum amount of precipitation during the
!                                   last 14 days,
!                                   below which the pulse is activated
!
! -------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,  INTENT(IN)    :: ndaylength, kproma, current_timestep
  REAL(dp), INTENT(IN)    :: timestep_length
  REAL(dp), INTENT(IN)    :: ls_prec(1:kproma), conv_prec(1:kproma)
  REAL(dp), INTENT(INOUT) :: rpulseregime(1:kproma)
  REAL(dp), INTENT(INOUT) :: pulseday(1:kproma), &
                             pulse(1:kproma), prec_hist(1:kproma,1:ndrydays+1)

  ! LOCAL
  INTEGER, PARAMETER      :: npulseregime = 3   ! number of pulsing regimes
  ! Constant numbers used in the calculation for the pulsing of each
  ! pulsing regime according to
  !   pls = plsfact * exp(plsexp * plsday)
  ! for the period of plsduration

  INTEGER, PARAMETER, DIMENSION(npulseregime)  :: &
       pulseduration = (/3, 7, 14/)
  REAL(dp), PARAMETER, DIMENSION(npulseregime) :: &
       pulsefact     = (/11.19, 14.68, 18.46/)
  REAL(dp), PARAMETER, DIMENSION(npulseregime) :: &
       pulseexp      = (/-0.805, -0.384, -0.208/)
  REAL(dp), PARAMETER, DIMENSION(npulseregime) :: &
       pulselimit    = (/0.001, 0.005, 0.015/)    ! in meter

  ! limit of precipitaion during the last 14 days,
  ! below which pulsing is activated
  ! can eventually be changed to soil moisture
  REAL(dp), PARAMETER :: pulseinitlimit = 0.01_dp                   ! in meter

  ! integer values of pulsing regime
  INTEGER, DIMENSION(kproma) :: pulseregime

  ! used for cumulation of precipitation and initilization of pulsing
  LOGICAL newday

  ! counter
  INTEGER i

  ! convert real values of rpulseregime to integer
  pulseregime = NINT(rpulseregime)

  newday = .false.
  if (INT(current_timestep * timestep_length/ndaylength) .GT.   &
      INT((current_timestep - 1) * timestep_length/ndaylength)) &
      newday = .true.

  ! count over grid points
  DO i=1, kproma
     ! increment day of pulsing and decrement pulse, if already active
     IF(pulseregime(i) .GE. 1._dp) THEN
        pulseday(i) = pulseday(i) + timestep_length/ndaylength
        pulse(i)    = pulsefact(pulseregime(i)) &
             * EXP(pulseexp(pulseregime(i))*pulseday(i))

        ! reset pulse and pulseregime if the pulse is over
        IF ((pulse(i) .LT. 1._dp)  .OR. (pulseday(pulseregime(i)) &
             .GT. pulseduration(pulseregime(i)))) THEN
           pulse(i)       = 1._dp
           pulseregime(i) = 0._dp
           pulseday(i)    = 0._dp
        END IF

     ! initialize pulsing if not yet active,
     !   a new day started and
     !   14 days cumulative precipitation below limit
     ELSE IF (newday .AND. SUM(prec_hist(i,1:ndrydays)) .LT. pulseinitlimit) THEN

        ! which pulsing regime, start with the strongest
        IF (prec_hist(i, ndrydays+1) .GE. pulselimit(3)) THEN
           pulseregime(i) = 3._dp
           pulseday(i)    = 1._dp
          pulse(i)        = pulsefact(3) * EXP(pulseexp(3))
        ELSE IF (prec_hist(i, ndrydays+1) .GE. pulselimit(2)) THEN
           pulseregime(i) = 2._dp
           pulseday(i)    = 1._dp
           pulse(i)       = pulsefact(2) * EXP(pulseexp(2))
        ELSE IF (prec_hist(i, ndrydays+1) .GE. pulselimit(1)) THEN
           pulseregime(i) = 1._dp
           pulseday(i)    = 1._dp
           pulse(i)       = pulsefact(1) * EXP(pulseexp(1))
        END IF
     END IF
  END DO   ! i=1, kproma

  ! shift precipitation history by one day, if a new day started
  IF (newday) THEN

     DO i=1, ndrydays
        prec_hist(:,i) = prec_hist(:,i+1)
     END DO
     prec_hist(:,ndrydays+1) = ls_prec + conv_prec
  ! accumulate todays precipitation
  ELSE
     prec_hist(:,ndrydays+1) = prec_hist(:,ndrydays+1) + ls_prec + conv_prec
  END IF

  ! convert INTEGER value of pulseregime to REAL value of output rpulseregime
  rpulseregime = REAL(pulseregime)

END SUBROUTINE NOemis_yl95sl10_pulsing
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE NOemis_yl95sl10_crf(kproma, month, latitude, noemclass_yl95sl10 &
     , lai, crf)
! -------------------------------------------------------------------------
! calculate the canopy reduction factor for NO soil emissions
! according to LAI (from ECHAM) and SAI (defined here) of the
! individual "ecosystems". Original algorithm is from Jacob and
! Bakwin (1991).
!
!       EXP(-8.75*sai) + EXP(-0.24*lai)
! crf = ------------------------------
!                    2
!
! Interface:
! ---------
!
! input:
!  - latitude           [Å∞]          to distinguish between tropics and none
!                                    tropics
!  - noemclass_yl95sl10 []           fraction of each "ecosystem" in the grid
!  - lai                [m^2/m^2]    leaf area index
!
! output:
!  - crf                []           canopy reduction factor [0 < crf <= 1]
!
! local variables:
!  - sai_yl95sl10_temp  [m^2/m^2]    stomatal area index for none tropics
!                                    [abs(latitude) >  30]
!  - sai_yl95sl10_trop  [m^2/m^2]    stomatal area index for tropics
!                                    [abs(latitude) <= 30]
!
! -------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,  INTENT(IN)    :: kproma
  INTEGER,  INTENT(IN)    :: month
  REAL(dp), INTENT(IN)    :: latitude(1:kproma)
  REAL(dp), INTENT(IN)    :: noemclass_yl95sl10(1:kproma, ncl_yl95sl10)
  REAL(dp), INTENT(IN)    :: lai(1:kproma)
  REAL(dp), INTENT(INOUT) :: crf(1:kproma)

  ! LOCAL
  ! adopt this to MODIS+Koeppen
  REAL(dp), PARAMETER, DIMENSION(ncl_yl95sl10) :: sai_yl95sl10_temp = &
       (/0., 0.,             &
        0., 0., 0.,          &
        0.,                  &
        0.01, 0.01,          &
        0.018, 0.018, 0.018, &
        0., 0.,              &
        0.02,                &
        0.025, 0.025, 0.025, &
        0.036, 0.036,        &
        0.075,               &
        0.12,                &
        0.032, 0.032, 0.032/)
  REAL(dp), PARAMETER, DIMENSION(ncl_yl95sl10) :: sai_yl95sl10_trop = &
       (/0., 0.,             &
        0., 0., 0.,          &
        0.,                  &
        0.01, 0.01,          &
        0.02, 0.02, 0.02,    &
        0., 0.,              &
        0.04,                &
        0.025, 0.025, 0.025, &
        0.036, 0.036,        &
        0.075,               &
        0.12,                &
        0.032, 0.032, 0.032/)

  REAL(dp), PARAMETER          :: ks = -8.75_dp
  REAL(dp), PARAMETER          :: kc = -0.24_dp
  REAL(dp), DIMENSION(kproma)  :: sai

  INTEGER i

  DO i=1, kproma
     ! temperate or tropical SAI, agriculture is treated special
     ! check if the calculation is ok since this is an exponential
     ! function and you can not sum it up so easily. A test program
     ! seemed to work that way, as it is implemented now.
     !
     ! northern hemispheric growing season
     IF(latitude(i) .GT. 30. .AND. month .GE. 5 .AND. month .LE. 8) THEN
        sai(i) = SUM(sai_yl95sl10_temp(1:ncl_yl95sl10 - 3) &
             * noemclass_yl95sl10(i, 1:ncl_yl95sl10 - 3)) + &
             (REAL(month - 5, dp)/3._dp) &
             * SUM(sai_yl95sl10_temp(ncl_yl95sl10 - 2:ncl_yl95sl10) &
             * noemclass_yl95sl10(i, ncl_yl95sl10 - 2:ncl_yl95sl10))
     ! southern hemispheric growing season (first two months)
     ELSE IF(latitude(i) .LT. -30. .AND. month .GE. 11) THEN
        sai(i) = SUM(sai_yl95sl10_trop(1:ncl_yl95sl10 - 3) &
             * noemclass_yl95sl10(i, 1:ncl_yl95sl10 - 3)) + &
             (REAL(month - 11, dp)/3._dp) &
             * SUM(sai_yl95sl10_temp(ncl_yl95sl10 - 2:ncl_yl95sl10) &
             * noemclass_yl95sl10(i, ncl_yl95sl10 - 2:ncl_yl95sl10))
     ! southern hemispheric growing season (last two months)
     ELSE IF(latitude(i) .LT. -30. .AND. month .LE. 2 ) THEN
        sai(i) = SUM(sai_yl95sl10_trop(1:ncl_yl95sl10-3) &
             * noemclass_yl95sl10(i, 1:ncl_yl95sl10-3)) + &
             (REAL(month + 1, dp)/3._dp) &
             * SUM(sai_yl95sl10_temp(ncl_yl95sl10-2:ncl_yl95sl10) &
             * noemclass_yl95sl10(i, ncl_yl95sl10-2:ncl_yl95sl10))
     ! tropical growing season all year round
     ELSE IF (ABS(latitude(i)) .LE. 30.) THEN
        sai(i) = SUM(sai_yl95sl10_trop(1:ncl_yl95sl10-3) &
             * noemclass_yl95sl10(i, 1:ncl_yl95sl10-3)) + &
             0.5_dp * SUM(sai_yl95sl10_temp(ncl_yl95sl10 - 2:ncl_yl95sl10) &
             * noemclass_yl95sl10(i, ncl_yl95sl10 - 2:ncl_yl95sl10))
     ! rest of the year, for temperate latitudes
     ELSE
        sai(i) = SUM(sai_yl95sl10_trop(1:ncl_yl95sl10 - 3) &
             * noemclass_yl95sl10(i, 1:ncl_yl95sl10 - 3))
     END IF
     crf(i) = (EXP(ks * sai(i)) + EXP(kc * lai(i))) / 2._dp

  ENDDO  ! i=1, kproma

END SUBROUTINE NOemis_yl95sl10_crf
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE NOemis_yl95sl10(kproma, month, latitude, longitude  &
                        , prec_hist, noemisclass, fert_rate    &
                        , smoist, vsm, stempin, noslflux       &
                        , noslflux_diag_yl95sl10)
! -------------------------------------------------------------------------
! calculate the NO emission from soils according to the algorithm
! by Yienger and Levy (1995). The "pulsing factor" and the "canopy
! reduction factor" are calculated in separate subroutines
!
! decide first if the soil is wet or dry
! then calculate the temperature dependent flux (no flux below 0Å∞C)
!
! A is the wet or dry emission factor, calculated here as
! biome-weighted factor of the individual ecosystem emission factors (Aw/d,biome)
! dry conditions (Ad):
!   flux = A * T / 30          0 < T < 30,   T in Å∞C
!   flux = A                   T > 30,       T in Å∞C
! wet conditions (Aw):
!   flux = 0.28 * A * T        0 < T < 10,   T in Å∞C
!   flux = A * EXP(0.103 * T)  10 < T < 30,  T in Å∞C
!   flux = 21.97 * A           T > 30,       T in Å∞C
!
! rainforest has constant emissions (dry for five months, wet the rest),
! -> very crude implementation
!
! agriculture is assumed to be wet all time, wet emission factor of grassland
! is assumed + fertilizer induced emission during the growing season. Growing
! season is all year round in the tropics (between 30Å∞S and 30Å∞N) and lasting
! 4 months in the none tropics.
! Rice growing in southeastern Asia is reduced by a factor of 30.
! Half of the rice growing area in India is reduced by a factor of 30.
!
! Interface:
! ----------
!
! input:
!  - month           number of month
!  - latitude        needed for tropics/none-tropics and rice-groing regions
!  - longitude       needed rice-groing regions
!  - noemisclass     fraction of ecosystem classes in gribbox
!  - fert_rate       fertilizer rate per gridbox
!  - smoist          soil water content (column)
!  - vsm             volumetric soil moisture
!  - stempin     [K] temperature of first soil layer
!
! output:
!  - noslflux        NO soil flux from gridbox
!  - noslflux_diag_yl95sl10   diagnostic, potential NO soil flux
!                             for each ecosystem
!
! local:
!  - itrop           index of rainforest in ecosystem data
!  - iagris          start index of agriculture in ecosystem data
!  - iagrie          end index of agriculture in ecosystem data
!  - tlimits_*       temperature limits for the upper calculation
!  - stemp           soil temperature in deg. Celsius
!  - fact_wet        constants used for calculation in above formula
!  - drylimit_prec   value to destinguish between wet and dry base
!                    on precipitation
!  - drylimit_vsm    value to destinguish between wet and dry base
!                    on soil moisture, based on volumetric soil moisture
!  - noslflux_agri   emission flux + FIE for agricultural areas
!  - trop_lat        latitude, below which the tropical region is
!  - ricearea_*_*    regions of rice cultivation:
!                    sea = south east asia: all reduced by a factor of XX
!                    ind = India (half rice, so half area reduced by a factor of XX)
!  - rice_reduct     rice reduction factor
!
! -------------------------------------------------------------------------
  INTEGER,  INTENT(IN)    :: kproma
  INTEGER,  INTENT(IN)    :: month
  REAL(dp), INTENT(IN)    :: latitude(1:kproma)
  REAL(dp), INTENT(IN)    :: longitude(1:kproma)
  REAL(dp), INTENT(IN)    :: prec_hist(1:kproma, 1:ndrydays+1)
  REAL(dp), INTENT(IN)    :: noemisclass(1:kproma, ncl_yl95sl10)
  REAL(dp), INTENT(IN)    :: fert_rate(1:kproma)
  REAL(dp), INTENT(IN)    :: smoist(1:kproma)
  REAL(dp), INTENT(IN)    :: vsm(1:kproma)
  REAL(dp), INTENT(IN)    :: stempin(1:kproma)
  REAL(dp), INTENT(INOUT) :: noslflux(1:kproma)
  REAL(dp), INTENT(INOUT) :: &
       noslflux_diag_yl95sl10(1:kproma, ncl_yl95sl10)

  ! LOCAL
  INTEGER,  PARAMETER               :: itrop           = 21
  INTEGER,  PARAMETER               :: iagris          = 22
  INTEGER,  PARAMETER               :: iagrie          = 24
  REAL(dp), PARAMETER, DIMENSION(3) :: tlimits_wet     = (/ 0., 10., 30. /)
  REAL(dp), PARAMETER, DIMENSION(3) :: fact_wet        = (/ 0.28, 0.103, 21.97 /)
  REAL(dp), PARAMETER, DIMENSION(2) :: tlimits_dry     = (/ 0., 30. /)
  REAL(dp), PARAMETER               :: drylimit_prec   = 0.01    ! [m]
  REAL(dp), PARAMETER               :: drylimit_smoist = 0.1     ! [m]
  REAL(dp), PARAMETER               :: drylimit_vsm    = 0.15    ! [fraction]
  REAL(dp), PARAMETER               :: trop_lat        = 30.
  REAL(dp), PARAMETER, DIMENSION(2) :: ricearea_sea_lon = (/ 80., 140. /)
  REAL(dp), PARAMETER, DIMENSION(2) :: ricearea_sea_lat = (/ 0., 35. /)
  REAL(dp), PARAMETER, DIMENSION(2) :: ricearea_ind_lon = (/ 60., 80. /)
  REAL(dp), PARAMETER, DIMENSION(2) :: ricearea_ind_lat = (/ 0., 35. /)
  REAL(dp)                          :: rice_reduct     = 30.

  INTEGER i
  LOGICAL                        :: dry = .TRUE.
  REAL(dp)                       :: noslflux_agri
  REAL(dp), DIMENSION(kproma)    :: stemp

  ! K-> deg C (soil temperature)
  stemp = stempin - 273.15_dp

  DO i=1, kproma
     ! check if wet or dry
     dry = .TRUE.
     ! -> precipitation history
     IF (SUM(prec_hist(i,1:ndrydays)) .GT. drylimit_prec .AND. &
          smoist_method .EQ. 0) THEN
        dry = .FALSE.
     ! -> soil water column
     ELSE IF (smoist(i) .GT. drylimit_smoist .AND. smoist_method .EQ. 1) THEN
        dry = .FALSE.
     ! volumetric soil moisture
     ELSE IF (vsm(i) .GT. drylimit_vsm .AND. smoist_method .EQ. 2) THEN
        dry = .FALSE.
     ENDIF

     noslflux(i) = 0._dp
     noslflux_diag_yl95sl10(i,:) = 0._dp
     ! calculation for all ecosystems without rainforest and agriculture,
     ! depending on the soil moisture state and soil temperature
     !
     ! hard coded number of ! 4 ! none tropic and agricultural classes
     ! this may change with other emission classes
     !
     IF (dry) THEN
        IF (stemp(i) .GT. tlimits_dry(1) .AND. stemp(i) .LE. tlimits_dry(2)) THEN

           noslflux(i) = SUM(noemfact_dry_yl95sl10(1:ncl_yl95sl10-4) * &
                noemisclass(i, 1:ncl_yl95sl10-4)) * (stemp(i) / tlimits_dry(2))
           noslflux_diag_yl95sl10(i,1:ncl_yl95sl10-4) =   &
                noemfact_dry_yl95sl10(1:ncl_yl95sl10-4) * &
                (stemp(i) / tlimits_dry(2))

        ELSE IF (stemp(i) .GT. tlimits_dry(2)) THEN

           noslflux(i) = SUM(noemfact_dry_yl95sl10(1:ncl_yl95sl10-4) * &
                noemisclass(i, 1:ncl_yl95sl10-4))
           noslflux_diag_yl95sl10(i, 1:ncl_yl95sl10-4) = &
                noemfact_dry_yl95sl10(1:ncl_yl95sl10-4)

        ENDIF
     ELSE
        IF (stemp(i) .GT. tlimits_wet(1) .AND. stemp(i) .LE. tlimits_wet(2)) THEN

           noslflux(i) = SUM(noemfact_wet_yl95sl10(1:ncl_yl95sl10-4) * &
                noemisclass(i, 1:ncl_yl95sl10-4)) * stemp(i) * fact_wet(1)
           noslflux_diag_yl95sl10(i, 1:ncl_yl95sl10-4) = &
                noemfact_wet_yl95sl10(1:ncl_yl95sl10-4) * stemp(i) * fact_wet(1)

        ELSE IF (stemp(i) .GT. tlimits_wet(2) .AND. &
             stemp(i) .LE. tlimits_wet(3)) THEN

           noslflux(i) = SUM(noemfact_wet_yl95sl10(1:ncl_yl95sl10-4) * &
                noemisclass(i, 1:ncl_yl95sl10-4)) * EXP(fact_wet(2)*stemp(i))
           noslflux_diag_yl95sl10(i, 1:ncl_yl95sl10-4) =  &
                noemfact_wet_yl95sl10(1:ncl_yl95sl10-4) * &
                EXP(fact_wet(2)*stemp(i))

        ELSE IF (stemp(i) .GT. tlimits_wet(3)) THEN
           noslflux(i) = SUM(noemfact_wet_yl95sl10(1:ncl_yl95sl10-4) * &
                noemisclass(i, 1:ncl_yl95sl10-4)) * fact_wet(3)
           noslflux_diag_yl95sl10(i, 1:ncl_yl95sl10-4) = &
                noemfact_wet_yl95sl10(1:ncl_yl95sl10-4) * fact_wet(3)
        ENDIF
     ENDIF ! dry

     ! calculation of emission from the rainforest for the northern hemisphere
     ! the dry season is assumed to be from May to September
     ! and for the southern hemisphere from November to March
     IF (latitude(i) .GE. 0. .AND. month .GE. 5 .AND. month .LE. 9) THEN
        noslflux(i) = noslflux(i) + &
             noemisclass(i, itrop) * noemfact_dry_yl95sl10(itrop)
        noslflux_diag_yl95sl10(i, itrop) = noemfact_dry_yl95sl10(itrop)
     ELSE IF (latitude(i) .GE. 0.) THEN
        noslflux(i) = noslflux(i) + &
             noemisclass(i, itrop) * noemfact_wet_yl95sl10(itrop)
        noslflux_diag_yl95sl10(i, itrop) = noemfact_wet_yl95sl10(itrop)
     ELSE IF (latitude(i) .LT. 0. .AND. (month .LE. 3 .OR. month .GE. 11)) THEN
        noslflux(i) = noslflux(i) + &
             noemisclass(i, itrop) * noemfact_dry_yl95sl10(itrop)
        noslflux_diag_yl95sl10(i, itrop) = noemfact_dry_yl95sl10(itrop)
     ELSE IF (latitude(i) .LT. 0.) THEN
        noslflux(i) = noslflux(i) + &
             noemisclass(i, itrop) * noemfact_wet_yl95sl10(itrop)
        noslflux_diag_yl95sl10(i, itrop) = noemfact_wet_yl95sl10(itrop)
     ENDIF

     noslflux_agri=0.
     ! calculation of agricultural fluxes, is always treated
     ! as wet emission
     ! at the moment all three "anthropogenic" land cover types are treatet equal
     IF (stemp(i) .GT. tlimits_wet(1) .AND. stemp(i) .LE. tlimits_wet(2)) THEN
        noslflux_agri = noemfact_wet_yl95sl10(iagris) * stemp(i) * fact_wet(1)
     ELSE IF (stemp(i) .GT. tlimits_wet(2) .AND. &
          stemp(i) .LE. tlimits_wet(3)) THEN
        noslflux_agri = noemfact_wet_yl95sl10(iagris) * EXP(fact_wet(2)*stemp(i))
     ELSE IF (stemp(i) .GT. tlimits_wet(3)) THEN
        noslflux_agri = noemfact_wet_yl95sl10(iagris) * fact_wet(3)
     ENDIF

     ! fertilizer input
     ! add fertilizer loss to noslflux
     ! below months are the agricultural growing season
     IF (latitude(i) .GT. trop_lat .AND. month .GE. 5 .AND. month .LE. 8) THEN
        noslflux_agri = noslflux_agri + fert_rate(i)
     ELSE IF (latitude(i) .LT. -trop_lat .AND. &
          (month .LE. 2 .OR. month .GE. 11)) THEN
        noslflux_agri = noslflux_agri + fert_rate(i)
     ELSE if (ABS(latitude(i)) .LE. trop_lat) THEN
        noslflux_agri = noslflux_agri + fert_rate(i)
     ENDIF

     ! reduction for the rice growing regions
     ! southeastern Asia
     IF (longitude(i) .GE. ricearea_sea_lon(1) .AND. &
          longitude(i) .LE. ricearea_sea_lon(2) &
          .AND. latitude(i) .GE. ricearea_sea_lat(1) &
          .AND. latitude(i) .LE. ricearea_sea_lat(2)) THEN
        noslflux_agri = noslflux_agri / rice_reduct
     ! India
     ELSE IF (longitude(i) .GE. ricearea_ind_lon(1) .AND.&
          longitude(i) .LE. ricearea_ind_lon(2) &
          .AND. latitude(i) .GE. ricearea_ind_lat(1) .AND. &
          latitude(i) .LE. ricearea_ind_lat(2)) THEN
        noslflux_agri = noslflux_agri * (rice_reduct + 1.) / (2. * rice_reduct)
     ENDIF

     ! add the agricultural flux to the total flux
     ! according to its gridbox fraction
     noslflux(i) = noslflux(i) + &
          noslflux_agri * sum(noemisclass(i, iagris:iagrie))
     noslflux_diag_yl95sl10(i, iagris:iagrie) = noslflux_agri

  ENDDO  ! i=1, kproma

  ! change units: from ng(N)/(m^2*s) to molec/(m^2*s)
  !noslflux      = ((noslflux      * 1.e-9) / 14.007) * N_A
  !noslflux_diag = ((noslflux_diag * 1.e-9) / 14.007) * N_A

END SUBROUTINE NOemis_yl95sl10
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE bioaer_emissions_olson(kproma                &
                     , olson_emis_seas                  &
                     , olson_emis_landice               &
                     , olson_emis_deserts               &
                     , olson_emis_forests               &
                     , olson_emis_grasslands            &
                     , olson_emis_crops                 &
                     , olson_emis_wetlands              &
                     , olson_emis_shrubs                &
                     , olson_emis_coastal               &
                     , olson_emis_urban                 &
                     , olson_emis_tundra                &
                     , olson, nolclass                  &
                     )


    ! -----------------------------------------------------------------------
    !
    !     **** Olson bioaerosol source
    !
    !     Authors:
    !     --------
    !     Susannah Burrows, MPI-Mainz                          2006
    !
    !     purpose
    !     -------
    !     Simulate emissions of bioaerosol particles from Olson
    !     lumped ecosystems.
    !     Emission fluxes are hardcoded below.  When used as passive tracers,
    !     tracer concentrations are independent and results scalable.
    !
    !     ---------
    !     method
    !     ---------
    !
    !     The aerosol is emitted as a constant flux from Olson
    !     lumped ecosystem classes.
    !     Aerosol properties are defined in the PTRAC namelist.
    !
    ! -----------------------------------------------------------------------
    ! Interface:
    ! ----------
    ! Input
    !
    ! olson: Olson ecosystem classification (1-72), reduced to nolclass classes
    !        Olson's World Ecosystems (1992).
    ! nolclass: number of ecosystems in reduced scheme (length of
    !          index dimension of olson_emis_flux, currently 11)
    !
    ! Input/Output
    !
    !  none
    !
    ! Output
    !
    ! Aerosol emissions to be sent directly into passive tracers
    !
    ! -----------------------------------------------------------------------

IMPLICIT NONE

! Input
    INTEGER,  INTENT(IN)  :: kproma
    ! number of reduced olson emission classes
    INTEGER,  INTENT(IN)  :: nolclass
    ! Percent of each ecosystem class in the gridbox
    REAL(dp), INTENT(IN)  :: olson(1:kproma,nolclass)

! Output
! Fluxes in lumped ecosystem classes
    REAL(dp), INTENT(OUT) :: olson_emis_seas(1:kproma)
    REAL(dp), INTENT(OUT) :: olson_emis_landice(1:kproma)
    REAL(dp), INTENT(OUT) :: olson_emis_deserts(1:kproma)
    REAL(dp), INTENT(OUT) :: olson_emis_forests(1:kproma)
    REAL(dp), INTENT(OUT) :: olson_emis_grasslands(1:kproma)
    REAL(dp), INTENT(OUT) :: olson_emis_crops(1:kproma)
    REAL(dp), INTENT(OUT) :: olson_emis_wetlands(1:kproma)
    REAL(dp), INTENT(OUT) :: olson_emis_shrubs(1:kproma)
    REAL(dp), INTENT(OUT) :: olson_emis_coastal(1:kproma)
    REAL(dp), INTENT(OUT) :: olson_emis_urban(1:kproma)
    REAL(dp), INTENT(OUT) :: olson_emis_tundra(1:kproma)

! Local variables
    REAL(dp) :: flux_seas
    REAL(dp) :: flux_landice
    REAL(dp) :: flux_deserts
    REAL(dp) :: flux_forests
    REAL(dp) :: flux_grasslands
    REAL(dp) :: flux_crops
    REAL(dp) :: flux_wetlands
    REAL(dp) :: flux_shrubs
    REAL(dp) :: flux_coastal
    REAL(dp) :: flux_urban
    REAL(dp) :: flux_tundra

!  Constant source strengths
    REAL(dp) :: flux_global
    !REAL(dp), DIMENSION(nolclass) :: flux

! Initialization
    olson_emis_seas(:)     = 0.0_dp
    olson_emis_landice(:)  = 0.0_dp
    olson_emis_deserts(:)  = 0.0_dp
    olson_emis_forests(:)  = 0.0_dp
    olson_emis_grasslands(:) = 0.0_dp
    olson_emis_crops(:)    = 0.0_dp
    olson_emis_wetlands(:) = 0.0_dp
    olson_emis_shrubs(:)   = 0.0_dp
    olson_emis_coastal(:)  = 0.0_dp
    olson_emis_urban(:)    = 0.0_dp
    olson_emis_tundra(:)   = 0.0_dp

! Source strengths (hardcoded here)
    flux_global     = 1.0E-9
    flux_seas       = flux_global
    flux_landice    = flux_global
    flux_deserts    = flux_global
    flux_forests    = flux_global
    flux_grasslands = flux_global
    flux_crops      = flux_global
    flux_wetlands   = flux_global
    flux_shrubs     = flux_global
    flux_coastal    = flux_global
    flux_urban      = flux_global
    flux_tundra     = flux_global

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! COMPUTATION OF BACTERIAL FLUX in [num m-2 s-1]:

! By ecosystem

   ! Oceans, Seas, Inland Waters
   olson_emis_seas(:)     = flux_seas*olson(:,1)
   olson_emis_urban(:)    = flux_urban*olson(:,2)
   olson_emis_shrubs(:)   = flux_shrubs*olson(:,3)
   olson_emis_forests(:)  = flux_forests*olson(:,4)
   olson_emis_deserts(:)  = flux_deserts*olson(:,5)
   olson_emis_landice(:)  = flux_landice*olson(:,6)
   olson_emis_crops(:)    = flux_crops*olson(:,7)
   olson_emis_wetlands(:) = flux_wetlands*olson(:,8)
   olson_emis_coastal(:)  = flux_coastal*olson(:,9)
   olson_emis_grasslands(:) = flux_grasslands*olson(:,10)
   olson_emis_tundra(:)   = flux_tundra*olson(:,11)

END SUBROUTINE bioaer_emissions_olson
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE bioaer_emissions_modis(kproma                &
                     , modis_emis_water                 &
                     , modis_emis_ever_need             &
                     , modis_emis_ever_broad            &
                     , modis_emis_deci_need             &
                     , modis_emis_deci_broad            &
                     , modis_emis_mixed_forest          &
                     , modis_emis_closed_shrubs         &
                     , modis_emis_open_shrubs           &
                     , modis_emis_woody_savannas        &
                     , modis_emis_savannas              &
                     , modis_emis_grasslands            &
                     , modis_emis_perm_wetlands         &
                     , modis_emis_crops                 &
                     , modis_emis_urban                 &
                     , modis_emis_crop_nature           &
                     , modis_emis_snow_ice              &
                     , modis_emis_barren                &
                     , modis_emis_unclass               &
                     , modis, nmodisclass               &
                     )


    ! -----------------------------------------------------------------------
    !
    !     **** MODIS bioaerosol source
    !
    !     Authors:
    !     --------
    !     Susannah Burrows, MPI-Mainz                          2009
    !
    !     purpose
    !     -------
    !     Simulate emissions of bioaerosol particles from MODIS ecosystems.
    !     Emission fluxes are hardcoded below.  When used as passive tracers,
    !     tracer concentrations are independent and results scalable.
    !
    !     ---------
    !     method
    !     ---------
    !
    !     The aerosol is emitted as a constant flux from MODIS
    !     ecosystem classes.
    !     Aerosol properties are defined in the PTRAC namelist.
    !
    ! -----------------------------------------------------------------------
    ! Interface:
    ! ----------
    ! Input
    !
    ! modis      : MODIS ecosystem classes
    ! nmodisclass: number of MODIS classes
    !
    ! Input/Output
    !
    !  none
    !
    ! Output
    !
    ! Aerosol emissions to be sent directly into passive tracers
    !
    ! -----------------------------------------------------------------------



IMPLICIT NONE

! Input
    INTEGER,  INTENT(IN)  :: kproma
    ! number of reduced olson emission classes
    INTEGER,  INTENT(IN)  :: nmodisclass
    ! Percent of each ecosystem class in the gridbox
    REAL(dp), INTENT(IN)  :: modis(1:kproma,nmodisclass)

! Output
! Fluxes in lumped ecosystem classes
    REAL(dp), INTENT(OUT) :: modis_emis_water(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_ever_need(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_ever_broad(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_deci_need(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_deci_broad(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_mixed_forest(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_closed_shrubs(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_open_shrubs(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_woody_savannas(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_savannas(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_grasslands(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_perm_wetlands(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_crops(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_urban(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_crop_nature(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_snow_ice(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_barren(1:kproma)
    REAL(dp), INTENT(OUT) :: modis_emis_unclass(1:kproma)

! Local variables
    REAL(dp) :: flux_water
    REAL(dp) :: flux_ever_need
    REAL(dp) :: flux_ever_broad
    REAL(dp) :: flux_deci_need
    REAL(dp) :: flux_deci_broad
    REAL(dp) :: flux_mixed_forest
    REAL(dp) :: flux_closed_shrubs
    REAL(dp) :: flux_open_shrubs
    REAL(dp) :: flux_woody_savannas
    REAL(dp) :: flux_savannas
    REAL(dp) :: flux_grasslands
    REAL(dp) :: flux_perm_wetlands
    REAL(dp) :: flux_crops
    REAL(dp) :: flux_urban
    REAL(dp) :: flux_crop_nature
    REAL(dp) :: flux_snow_ice
    REAL(dp) :: flux_barren
    REAL(dp) :: flux_unclass

!  Constant source strengths
    REAL(dp) :: flux_global

! Initialization
    modis_emis_water(:) = 0.0_dp
    modis_emis_ever_need(:) = 0.0_dp
    modis_emis_ever_broad(:) = 0.0_dp
    modis_emis_deci_need(:) = 0.0_dp
    modis_emis_deci_broad(:) = 0.0_dp
    modis_emis_mixed_forest(:) = 0.0_dp
    modis_emis_closed_shrubs(:) = 0.0_dp
    modis_emis_open_shrubs(:) = 0.0_dp
    modis_emis_woody_savannas(:) = 0.0_dp
    modis_emis_savannas(:) = 0.0_dp
    modis_emis_grasslands(:) = 0.0_dp
    modis_emis_perm_wetlands(:) = 0.0_dp
    modis_emis_crops(:) = 0.0_dp
    modis_emis_urban(:) = 0.0_dp
    modis_emis_crop_nature(:) = 0.0_dp
    modis_emis_snow_ice(:) = 0.0_dp
    modis_emis_barren(:) = 0.0_dp
    modis_emis_unclass(:) = 0.0_dp

! Source strengths (hardcoded here)
    flux_global = 1
    flux_water = flux_global
    flux_ever_need = flux_global
    flux_ever_broad = flux_global
    flux_deci_need = flux_global
    flux_deci_broad = flux_global
    flux_mixed_forest = flux_global
    flux_closed_shrubs = flux_global
    flux_open_shrubs = flux_global
    flux_woody_savannas = flux_global
    flux_savannas = flux_global
    flux_grasslands = flux_global
    flux_perm_wetlands = flux_global
    flux_crops = flux_global
    flux_urban = flux_global
    flux_crop_nature = flux_global
    flux_snow_ice = flux_global
    flux_barren = flux_global
    flux_unclass = flux_global

! COMPUTATION OF FLUX in [num m-2 s-1]:

! By ecosystem

   ! Oceans, Seas, Inland Waters
    modis_emis_water(:) = flux_water*modis(:,1)
    modis_emis_ever_need(:) = flux_ever_need*modis(:,2)
    modis_emis_ever_broad(:) = flux_ever_broad*modis(:,3)
    modis_emis_deci_need(:) = flux_deci_need*modis(:,4)
    modis_emis_deci_broad(:) = flux_deci_broad*modis(:,5)
    modis_emis_mixed_forest(:) = flux_mixed_forest*modis(:,6)
    modis_emis_closed_shrubs(:) = flux_closed_shrubs*modis(:,7)
    modis_emis_open_shrubs(:) = flux_open_shrubs*modis(:,8)
    modis_emis_woody_savannas(:) = flux_woody_savannas*modis(:,9)
    modis_emis_savannas(:) = flux_savannas*modis(:,10)
    modis_emis_grasslands(:) = flux_grasslands*modis(:,11)
    modis_emis_perm_wetlands(:) = flux_perm_wetlands*modis(:,12)
    modis_emis_crops(:) = flux_crops*modis(:,13)
    modis_emis_urban(:) = flux_urban*modis(:,14)
    modis_emis_crop_nature(:) = flux_crop_nature*modis(:,15)
    modis_emis_snow_ice(:) = flux_snow_ice*modis(:,16)
    modis_emis_barren(:) = flux_barren*modis(:,17)
    modis_emis_unclass(:) = flux_unclass*modis(:,18)

END SUBROUTINE bioaer_emissions_modis
! -------------------------------------------------------------------------

SUBROUTINE bio_emissions_lai(kproma                 &
                     , heald_emis                   &
                     , js_emis                      &
                     , hummel_emis                  &
                     , modis_lai                    &
                     , qm1                          &
                     , tm1                          &
                     , philat_2d                    &
                     , imonth                       &
                     )


    ! -----------------------------------------------------------------------
    !
    !     **** bioaerosol sources proportional to LAI
    !          (tracers are independent and results scaleable)
    !          heald_emis
    !
    !     Authors:
    !     --------
    !     Susannah Burrows, MPIC Mainz                          2009
    !
    !     purpose
    !     -------
    !     test transport of large aerosol particles and investigate of
    !     transport dependence on source region or ecosystem.
    !
    !     ---------
    !     method
    !     ------
    !
    !     The aerosol is emitted from land surfaces as a function of LAI.
    !
    !     Aerosol properties are defined in the PTRAC namelist.
    !
    ! -----------------------------------------------------------------------
    ! Interface:
    ! ----------
    ! Input
    !
    ! modis_lai    : Leaf Area Index for current month (to be read in via
    !                onlem namelist), km^2 km^-2
    ! qm1          : Water vapor concentration (at surface)
    ! philat_2d    : latitude (degree)
    ! imonth       : month index
    !
    ! Input/Output
    !
    !  none
    !
    ! Output
    !
    ! heald_emis: Estimated fungal spore emissions following Heald and Spracklen (2009, GRL)
    !             Given by Hoose as:
    !             F_fungi = 500 m^-2 s^-1 * LAI/5 * q/(1.5 * 10^-2 kg kg^-1)
    ! hummel_emis: Estimated fungal spore emissions following Hummel et al, 2015 (ACP)
    !              F_fungi =
    !
    ! js_emis: Estimated pollen emissions following Jacobson and Streets (2009, JGR)
    !          and Hoose et al. (2009, poster from EUCAARI conference):
    !          F_pollen = 0.5 m^-2 s^-1 * LAI * R_month
    !          with R_month = (0.5,0.5,0.5,2.0,2.0,2.0,1.0,1.0,1.0,0.5,0.5,0.5)
    !               for NH, and offset by six months for SH
    !
    ! -----------------------------------------------------------------------

IMPLICIT NONE

! Input
    INTEGER,  INTENT(IN)  :: kproma                 ! Number of local longitudes
    REAL(dp), INTENT(IN)  :: qm1(1:kproma)          ! Water vapor content (kg kg-1)
    REAL(dp), INTENT(IN)  :: tm1(1:kproma)          ! Water vapor content (kg kg-1)
    REAL(dp), INTENT(IN)  :: philat_2d(1:kproma)    ! latitude (degree)
    REAL(dp), INTENT(IN)  :: modis_lai(1:kproma)    ! Leaf Area Index
    INTEGER,  INTENT(IN)  :: imonth                 ! month index

! Output
    REAL(dp), INTENT(OUT) :: heald_emis(1:kproma)
    REAL(dp), INTENT(OUT) :: js_emis(1:kproma)
    REAL(dp), INTENT(OUT) :: hummel_emis(1:kproma)
! Locally-defined constants
    REAL(dp), DIMENSION(12) :: r_month_nh = (/ 0.5,0.5,0.5,2.0,2.0,2.0,1.0,1.0,1.0,0.5,0.5,0.5 /)
    REAL(dp), DIMENSION(12) :: r_month_sh = (/ 1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.5,0.5,2.0,2.0,2.0 /)

! Initialization


! COMPUTATION OF FLUX in [num m-2 s-1]:
   heald_emis =  50000. / 7.5 * modis_lai(1:kproma) * qm1(1:kproma)
! COMPUTATION OF FLUX (FUNGAL, Hummel, 2015)
   hummel_emis = 20.426 * (tm1(1:kproma) - 275.82) + 39300. * qm1(1:kproma) * modis_lai(1:kproma)
! Test which hemisphere this is and use the appropriate set of weights
WHERE (philat_2d(1:kproma) > 0.0_dp)  ! Northern hemisphere
   js_emis(1:kproma) = 0.5 * modis_lai(1:kproma) * r_month_nh(imonth)
ENDWHERE
WHERE (philat_2d(1:kproma) <= 0.0_dp)  ! Southern hemisphere
   js_emis(1:kproma) = 0.5 * modis_lai(1:kproma) * r_month_sh(imonth)
ENDWHERE

END SUBROUTINE bio_emissions_lai
! -------------------------------------------------------------------------
SUBROUTINE so2_emissions(zdz, so2, so2emflux)

  USE messy_main_constants_mem, ONLY: MS, MO, N_A

  ! incoming fluxes in resp. level  unit kg(SO2) m-2 s-1
  ! ??? is it realy so2 and not S ????
  ! netcdf tells S / P. Stiers code tells SO2
  REAL(dp), INTENT(IN)  :: zdz(:,:)
  REAL(dp), INTENT(IN)  :: so2(:,:)
  ! outgoing fluxes in resp. level  unit molecules(SO2) m-3 s-1
  REAL(dp), INTENT(OUT) :: so2emflux(:,:)

  so2emflux(:,:) = so2(:,:) /((MS+ 2._dp*MO)*1.e-3_dp) * N_A / zdz(:,:)

END SUBROUTINE so2_emissions
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE airsnow_emissions( kproma, snow_air_flux_br2, snow_air_flux_brcl, &
                              tsurf_2d, cvs, seaice, cossza_2d,              &
                              sic_multi_year,ddepflux_HOBr, ddepflux_BrNO3,  &
                              ddepflux_HBr, ddepflux_O3 )

  !------------------------------------------------------------------------
  ! AirSnow emission
  !  compute the emission of Br2/BrCl from sea ice and snow covered land
  !  following the parametrized scheme of Toyota et al. 2011
  !  IN:
  !     kproma         - length of field
  !     tsurf_2d       - surface temperature
  !     cvs            - fraction of snow cover on land
  !     seaice         - fraction of ice cover on ocean
  !     cossza_2d      - cosine of sun's zenith angle
  !     sic_multi_year - fraction of multi-year seaice
  !  OUT:
  !     snow_air_flux_br2  - flux of Br2 from parametrization
  !     snow_air_flux_brcl - flux of BrCl from parametrization
  !
  !------------------------------------------------------------------------

  USE messy_main_constants_mem, ONLY: DTR, T0

  IMPLICIT NONE

  ! Input
  INTEGER,  INTENT(IN)  :: kproma
  REAL(dp), DIMENSION(:), INTENT(IN) :: tsurf_2d
  REAL(dp), DIMENSION(:), INTENT(IN) :: cvs
  REAL(dp), DIMENSION(:), INTENT(IN) :: seaice
  REAL(dp), DIMENSION(:), INTENT(IN) :: sic_multi_year
  REAL(dp), DIMENSION(:), INTENT(IN) :: cossza_2d
  REAL(dp), DIMENSION(:), INTENT(IN) :: ddepflux_HOBr
  REAL(dp), DIMENSION(:), INTENT(IN) :: ddepflux_BrNO3
  REAL(dp), DIMENSION(:), INTENT(IN) :: ddepflux_HBr
  REAL(dp), DIMENSION(:), INTENT(IN) :: ddepflux_O3
  ! Output fluxes
  REAL(dp), DIMENSION(:), INTENT(OUT) :: snow_air_flux_br2
  REAL(dp), DIMENSION(:), INTENT(OUT) :: snow_air_flux_brcl
  ! Internal field
  REAL(dp), DIMENSION(1:kproma) :: sic_first_year

  ! Initialisation
  sic_first_year = seaice-sic_multi_year
  ! Make sure we don't get negative fluxes if there is more mysic than seaice
  WHERE ( sic_first_year < 0 )
     sic_first_year = 0._dp
  ENDWHERE
  ! Emission triggered by dry deposition, on sea ice and snow covered land
  ! below critical temperature
  WHERE ( tsurf_2d <= r_temp_crit + T0 )
     ! Br2 release from first year sea ice
     WHERE ( seaice > 0._dp )
        ! Br2 release from O3 deposition on sea ice dependent on sun light
        ! No release from multi-year sea ice
        WHERE ( cossza_2d <= COS(r_sun_theta_crit*DTR) )
           ! sunlit
           snow_air_flux_br2 = ddepflux_O3*r_trigger_1(2) &
                * sic_first_year
        ELSEWHERE
           ! dark
           snow_air_flux_br2 = ddepflux_O3*r_trigger_1(1) &
                * sic_first_year
        ENDWHERE
        ! Br2 release from HBr, HOBr, and BrNO3
        snow_air_flux_br2 = snow_air_flux_br2 + &
             (ddepflux_HOBr + ddepflux_BrNO3) * sic_first_year
        WHERE ( sic_multi_year > 0._dp )
           WHERE ( ddepflux_HOBr + ddepflux_BrNO3 < ddepflux_HBr )
              snow_air_flux_br2 = snow_air_flux_br2 + &
                   (ddepflux_HOBr + ddepflux_BrNO3) * sic_multi_year
           ELSEWHERE
              snow_air_flux_br2 =  snow_air_flux_br2 + &
                   0.5 * (ddepflux_HOBr + ddepflux_BrNO3 - ddepflux_HBr) &
                   * sic_multi_year
              snow_air_flux_brcl =  snow_air_flux_brcl + &
                   0.5 * (ddepflux_HOBr + ddepflux_BrNO3 - ddepflux_HBr) &
                   * sic_multi_year
           ENDWHERE
        ENDWHERE
     ENDWHERE
     ! Br2 release from snow covered land (no BrCl source)
     WHERE ( cvs > 0._dp )
        WHERE ( ddepflux_HOBr + ddepflux_BrNO3 < ddepflux_HBr )
           snow_air_flux_br2 = snow_air_flux_br2 + &
                (ddepflux_HOBr + ddepflux_BrNO3) * cvs
        ELSEWHERE
           snow_air_flux_br2 = snow_air_flux_br2 + &
                ddepflux_HBr * cvs
        ENDWHERE
        snow_air_flux_br2 = snow_air_flux_br2 + &
             (ddepflux_O3*r_trigger_1(3) * cvs)
     ENDWHERE
  ENDWHERE

END SUBROUTINE airsnow_emissions
! -------------------------------------------------------------------------

! ****************************************************************************
END MODULE messy_onemis
! ****************************************************************************
