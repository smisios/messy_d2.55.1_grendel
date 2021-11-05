MODULE MESSY_SCAV_MEM

!  AUTHOR:   HOLGER TOST, MPI CHEMIE, MAINZ
!            LAST MODIFIED 10.01.2005

! IN THIS MODULE THE SCAV_IDTS AND THE VALUES FOR CALCULATION OF THE
! TRANSFER COEFFICIENT ARE DEFINED


  USE MESSY_SCAV_I_KPP,         ONLY: ISPEC => NSPEC
  USE MESSY_SCAV_L_KPP,         ONLY: DP, LSPEC => NSPEC

  USE MESSY_MAIN_CONSTANTS_MEM, ONLY: I4, STRLEN_MEDIUM
  IMPLICIT NONE
  PRIVATE :: DP!, LSPEC, ISPEC
  SAVE
  INTRINSIC :: NULL

  LOGICAL :: LSCAV        ! GLOBAL SWITCH :             SCAV ON/OFF
  LOGICAL :: LSCAV_GP = .TRUE.    ! GLOBAL SWITCH :  GRID POINT SCAV ON/OFF
  LOGICAL :: LSCAV_LG = .TRUE.    ! GLOBAL SWITCH :  LAGRANGIAN SCAV ON/OFF
  LOGICAL :: LSCAV_LS     ! GLOBAL SWITCH : LARGE-SCALE SCAV ON/OFF
  LOGICAL :: LSCAV_CV     ! GLOBAL SWITCH :  CONVECTIVE SCAV ON/OFF
  LOGICAL :: LSCAV_GAS    ! GLOBAL SWITCH : GAS TRACERS SCAV ON/OFF
  LOGICAL :: LSCAV_AER    ! GLOBAL SWITCH :     AEROSOL SCAV ON/OFF
  LOGICAL :: LSCAV_L      ! GLOBAL SWITCH : SCAVENGING BY LIQUID WATER
  LOGICAL :: LSCAV_I      ! GLOBAL SWITCH : SCAVENGING BY SNOW/ICE
  LOGICAL :: LSCAV_NUC    ! GLOBAL SWITCH : NUCLEATION SCAVENGING
  LOGICAL :: LSCAV_IMP    ! GLOBAL SWITCH : IMPACTION  SCAVENGING
  LOGICAL :: LSCAV_EASY   ! SWITCH : NO_KPP=EASY SCAV ON/OFF FOR LIQUID
  LOGICAL :: ISCAV_EASY   ! SWITCH : NO_KPP=EASY SCAV ON/OFF FOR ICE
  INTEGER :: L_SCAV_EASY  ! SWITCH LIQUID: 1 FOR FIXED COEFFICIENTS
                          !                2 FOR ABSOLUTE HENRY'S LAW
                          !                3 FOR EFFECTIVE HENRY'S LAW
  INTEGER :: I_SCAV_EASY  ! SWITCH ICE:    1 FOR FIXED COEFFICIENTS
                          !                2 FOR PSEUDO HENRY's LAW
                          !                3 FOR ITERATIVE LANGMUIR UPTAKE
  LOGICAL :: LOVERWRITE_HENRY = .FALSE.
  LOGICAL :: LOVERWRITE_ALPHA = .FALSE.
  LOGICAL :: LUSE_EMPIRIC_ALPHA = .FALSE.


  LOGICAL :: L_LG         ! FALSE = GRIDPOINT SCAVENGING (ACTUALLY CALCULATION)
                          ! TRUE  = LAGRANGE SCAVENGING  (ACTUALLY CALCULATION)
    
! SWITCH OF DEGREE OF AEROSOL / GASPHASE COUPLING
  INTEGER(I4) :: CPL_AEROSOL      
   ! SWITCH OF DROPLET TO AEROSOL EVAPORATION
  INTEGER(I4) :: I_EVAP   ! 1 = AEROSOLS GO INTO THE MODE WHERE 
                          !     THEY ORIGINATE FROM
                          ! 2 = AEROSOLS GO TO THE LARGEST SOLUBLE 
                          !     MODE FOR THE AEROSOL MODULE (OR SUBMODEL
                          !     SPECIFIC TREATMENT)

  ! FRACTION OF NUMBER CONC. IN RESIDUAL (I. E. ``SOLUBLE'') MODE WITH LARGEST
  ! INDEX AFTER DROPLET EVAPORATION
  ! NOTE: FOR MADE3 NUMBER CONC. IN ALL RESIDUAL MODES WILL BE REDUCED TO THIS
  !       FRACTION
  REAL(dp)    :: FRAC_RESNUM
  ! ICE NUCLEATION RATES FOR T > 238K AND FOR T < 238 (HETEROGENEOUS AND
  ! HOMOGENEOUS)
  REAL(dp)    :: ISCAV_RATE_HITMP, ISCAV_RATE_LOTMP_HET, ISCAV_RATE_LOTMP_HOM
  INTEGER(I4) :: NMAXKPP = 10000   ! MAXIMUM NUMBER OF KPP CALCULATED TIMESTEPS


! MOLECULES IN LARGE SCALE PRECIPITATION
  REAL(DP) ,PUBLIC, POINTER :: KPP_RAIN_LS(:,:,:,:) => NULL()   
  REAL(DP) ,PUBLIC, POINTER :: KPP_SNOW_LS(:,:,:,:) => NULL()   
! MOLECULES IN CONVECTIVE PRECIPITATION
  REAL(DP) ,PUBLIC, POINTER :: KPP_RAIN_CV(:,:,:,:) => NULL()  
  REAL(DP) ,PUBLIC, POINTER :: KPP_SNOW_CV(:,:,:,:) => NULL()  

  TYPE AERO_FLUX_PTR_ARRAY
! AEROSOL MOLECULES IN LARGE SCALE PRECIPITATION
    REAL(DP), DIMENSION(:,:,:,:),   POINTER :: AER_FLX_LS         => NULL()
    REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: AERO_FIELD_LS      => NULL()
    REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: CLOUD_FIELD_LS_AER => NULL()
! AEROSOL MOLECULES IN CONVECTIVE PRECIPITATION
    REAL(DP), DIMENSION(:,:,:,:),   POINTER :: AER_FLX_CV         => NULL()
    REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: AERO_FIELD_CV      => NULL()
    REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: CLOUD_FIELD_CV_AER => NULL()
  END TYPE AERO_FLUX_PTR_ARRAY
  TYPE(AERO_FLUX_PTR_ARRAY), DIMENSION(:), POINTER :: AERO_FLX => NULL()

  TYPE ATTR_2D_PTR_ARRAY
    INTEGER, DIMENSION(:,:), POINTER :: INT_ATT => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: EVAP_TRAC => NULL()
    LOGICAL, DIMENSION(:,:), POINTER :: LOG_ATT => NULL()
  END TYPE ATTR_2D_PTR_ARRAY
  TYPE(ATTR_2D_PTR_ARRAY), DIMENSION(:), POINTER :: ATTR => NULL()

!   INDEX NUMBERS OF LOGICAL TRACER ATTRIBUTES
  INTEGER, PARAMETER  :: LAEROSOL = 1
  INTEGER, PARAMETER  :: LWETDEP  = 2

  INTEGER, PARAMETER  :: LOG_MAX = 2

!   INDEX NUMBERS OF INTEGER TRACER ATTRIBUTES
  INTEGER, PARAMETER  :: LQUANTITY       = 1
  INTEGER, PARAMETER  :: TMODE           = 2
  INTEGER, PARAMETER  :: AEROSOL_INDEX   = 3
  INTEGER, PARAMETER  :: EVAP_TRAC_IDT   = 4
  INTEGER, PARAMETER  :: EVAP_TRAC_IDT2  = 5
  INTEGER, PARAMETER  :: AER_MOD_NUM     = 6

  INTEGER, PARAMETER  :: INT_MAX = 6


  INTEGER   :: LSPEC_GAS
  INTEGER   :: LSPEC_LIQ
  INTEGER   :: ISPEC_GAS
  INTEGER   :: ISPEC_ICE
  TYPE IDX_L_2D_PTR_ARRAY
    INTEGER, DIMENSION(:,:), POINTER ::  GAS_SPEC => NULL()
    INTEGER, DIMENSION(:,:), POINTER ::  LIQ_SPEC => NULL()
    REAL(DP),DIMENSION(:,:), POINTER ::  GAS_ATTR => NULL()
    REAL(DP),DIMENSION(:,:), POINTER ::  LIQ_ATTR => NULL()
  END TYPE IDX_L_2D_PTR_ARRAY
  TYPE IDX_I_2D_PTR_ARRAY
    INTEGER, DIMENSION(:,:), POINTER ::  GAS_SPEC => NULL()
    INTEGER, DIMENSION(:,:), POINTER ::  ICE_SPEC => NULL()
    REAL(DP),DIMENSION(:,:), POINTER ::  GAS_ATTR => NULL()
    REAL(DP),DIMENSION(:,:), POINTER ::  ICE_ATTR => NULL()
  END TYPE IDX_I_2D_PTR_ARRAY
  TYPE(IDX_L_2D_PTR_ARRAY), DIMENSION(:), POINTER :: KPP_L_IDX => NULL()
  TYPE(IDX_I_2D_PTR_ARRAY), DIMENSION(:), POINTER :: KPP_I_IDX => NULL()
  TYPE(IDX_L_2D_PTR_ARRAY)                        :: KPP0_L_IDX
  TYPE(IDX_I_2D_PTR_ARRAY)                        :: KPP0_I_IDX

! FOR GAS_SPEC
  INTEGER, PARAMETER  :: GAS_IDX  = 1   ! index of gas phase compounds 
                                        !     within the species array
  INTEGER, PARAMETER  :: GAS2TRAC = 2   ! index of tracer for gas phase
                                        !     compound of gas species
  INTEGER, PARAMETER  :: GAS_I2L  = 3   ! index of corresponding ice and 
                                        !     liquid compounds, i.e. to 
                                        !     access liquid henry values
                                        !     for the gas phase compounds
                                        !     of the ice phase
  INTEGER, PARAMETER  :: GAS_G2P  = 4   ! index of corresponding species
                                        !     in gas phase and liquid / ice
                                        !     phase
  INTEGER, PARAMETER  :: GAS_MAX  = 4
! FOR GAS_ATTR
  INTEGER, PARAMETER  :: GAS_MW      = 1   ! molar weight / 100.
!  INTEGER, PARAMETER  :: GAS_PSS     = 2   ! effective henry for drydep
 ! INTEGER, PARAMETER  :: GAS_dryreac = 3   ! dryreac_sf for drydep
  INTEGER, PARAMETER  :: GATT_MAX    = 1
! FOR LIQ_SPEC  
  INTEGER, PARAMETER  :: LIQ_IDX  = 1   ! index of liquid phase compounds
                                        !     within the species array
  INTEGER, PARAMETER  :: LIQ2EVAP = 2   ! index of evaporation for species 
                                        !     into tracer_idt
  INTEGER, PARAMETER  :: LIQ2TRAC = 3   ! index of tracer for liquid phase
                                        !     compound of liquid species
  INTEGER, PARAMETER  :: LIQ2EVA2 = 4   ! index of evaporation for species 
                                        !     into tracer_idt 
                                        !     (if evaporation is splitted 
                                        !     between modes)
  INTEGER, PARAMETER  :: LIQ_MAX  = 4
! FOR LIQ_ATTR
  INTEGER, PARAMETER  :: LIQ_MW     = 1   ! molar weight / 100.
!  INTEGER, PARAMETER  :: LIQ_dens   = 2   ! density
!  INTEGER, PARAMETER  :: LIQ_charge = 3   ! compound charge
  INTEGER, PARAMETER  :: LATT_MAX   = 1
! FOR ICE_SPEC  
  INTEGER, PARAMETER  :: ICE_IDX  = 1   ! index of ice phase compounds
                                        !      within the species array
  INTEGER, PARAMETER  :: ICE2EVAP = 2   ! index of evaporation for species 
                                        !      into tracer_idt
  INTEGER, PARAMETER  :: ICE2LIQ  = 3   ! index of corresponding ice and liquid
                                        !      phase species for melting
  INTEGER, PARAMETER  :: ICE2TRAC = 4   ! index of tracer for liquid phase
                                        !     compound of liquid species
  INTEGER, PARAMETER  :: ICE_MAX  = 4
! FOR ICE_ATTR
  INTEGER, PARAMETER  :: ICE_MW     = 1   ! molar weight / 100.
!  INTEGER, PARAMETER  :: ICE_dens   = 2   ! density
!  INTEGER, PARAMETER  :: ICE_charge = 3   ! compound charge
  INTEGER, PARAMETER  :: IATT_MAX   = 1

  TYPE COUNT_1D_PTR_ARRAY
    INTEGER,          DIMENSION(:), POINTER :: NUMBERS => NULL()
  END TYPE COUNT_1D_PTR_ARRAY
  TYPE(COUNT_1D_PTR_ARRAY), DIMENSION(:), POINTER :: AER_COUNT => NULL()
  INTEGER, PARAMETER :: C_ALL        = 1
  INTEGER, PARAMETER :: C_MASS       = 2
  INTEGER, PARAMETER :: C_NUM        = 3
  INTEGER, PARAMETER :: AERCOUNT_MAX = 3

 
  TYPE AER_2D_PTR_ARRAY
    INTEGER, DIMENSION(:,:), POINTER :: AER_SPEC => NULL()
    CHARACTER(2*STRLEN_MEDIUM + 1), DIMENSION(:), POINTER :: NAME => NULL()
    REAL(dp),                       DIMENSION(:), POINTER :: MW   => NULL()
  END TYPE AER_2D_PTR_ARRAY
  TYPE(AER_2D_PTR_ARRAY), DIMENSION(:), POINTER :: AER_ATTR => NULL()

!FOR AER_SPEC
  INTEGER, PARAMETER  :: AER_IDX     = 1
  INTEGER, PARAMETER  :: AER2L_IDX   = 2
  INTEGER, PARAMETER  :: AER2I_IDX   = 3
  INTEGER, PARAMETER  :: AER2TRACL   = 4
  INTEGER, PARAMETER  :: AER2TRACI   = 5
  INTEGER, PARAMETER  :: AER_MAX     = 5

  TYPE PH_3D_PTR_ARRAY
! LS CLOUD PH[0-14, NA=19.99]
    REAL(DP) ,DIMENSION(:,:,:), POINTER :: CLOUDPH_LS  => NULL() 
! CV CLOUD PH[0-14, NA=19.99]                 
! IN CONVECTIVE CLOUDS DUE TO FRESHLY FORMED PRECIPITATION
    REAL(DP) ,DIMENSION(:,:,:), POINTER :: CLOUDPH_CV  => NULL() 
! LS CLOUD H+ CONCENTRATION[0-14, NA=19.99]
    REAL(DP) ,DIMENSION(:,:,:), POINTER :: HP_CLOUD_LS => NULL() 
! CV CLOUD H+ CONCENTRATION[0-14, NA=19.99] 
! IN CONVECTIVE CLOUDS DUE TO FRESHLY FORMED PRECIPITATION
    REAL(DP) ,DIMENSION(:,:,:), POINTER :: HP_CLOUD_CV => NULL() 
! RAIN PH [0-14, NA=19.99]
    REAL(DP) ,DIMENSION(:,:,:), POINTER :: RAINPH_LS   => NULL() 
! RAIN PH [0-14, NA=19.99] IN CONVECTIVE PRECIPITATION
    REAL(DP) ,DIMENSION(:,:,:), POINTER :: RAINPH_CV   => NULL()
  END TYPE PH_3D_PTR_ARRAY
  TYPE(PH_3D_PTR_ARRAY), DIMENSION(:), POINTER :: PH => NULL()

  TYPE KPP_STEPS_3D
    REAL(dp), DIMENSION(:,:,:), POINTER :: ls_cloud_steps => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: ls_rain_steps  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: cv_cloud_steps => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: cv_rain_steps  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: ls_cloud_rsteps => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: ls_rain_rsteps  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: cv_cloud_rsteps => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: cv_rain_rsteps  => NULL()
  END TYPE KPP_STEPS_3D
  TYPE(KPP_STEPS_3D), DIMENSION(:), POINTER :: KPP_INFO   => NULL()

  TYPE Prod_C_3D_PTR_ARRAY
    ! LS Cloud Production rates
    REAL(DP), DIMENSION(:,:,:), POINTER :: cl_ls => NULL()
    ! LS Rain Production Rates
    REAL(DP), DIMENSION(:,:,:), POINTER :: ra_ls => NULL()
    ! CV Cloud Production rates
    REAL(DP), DIMENSION(:,:,:), POINTER :: cl_cv => NULL()
    ! CV Rain Production rates
    REAL(DP), DIMENSION(:,:,:), POINTER :: ra_cv => NULL()
  END TYPE Prod_C_3D_PTR_ARRAY
  TYPE PROD_B_1D_PTR_ARRAY
    TYPE(Prod_C_3D_PTR_ARRAY), DIMENSION(:), POINTER :: Prod => NULL()
    INTEGER, DIMENSION(:),                   POINTER :: PROD_IDX => NULL()
    CHARACTER(LEN=2), DIMENSION(:),          POINTER :: PROD_CHAR => NULL()
  END TYPE PROD_B_1D_PTR_ARRAY
  TYPE(Prod_B_1D_PTR_ARRAY), DIMENSION(:), POINTER :: Pr => NULL()
  INTEGER :: nprod


  LOGICAL   :: L_ANYNCONVECT
    
  REAL(DP)  :: HS0(0:LSPEC)        ! SOLUTION CONSTANT
  REAL(DP)  :: DHT(0:LSPEC)        ! REACTION ENTHALPIE FOR HENRY
  REAL(DP)  :: ALPHA0(0:LSPEC)       ! ACCOMODATIION coefficient
  REAL(DP)  :: ALPHA_T(0:LSPEC)      ! temperature dependency of alpha0
  REAL(DP)  :: SCAV_VALUE(0:LSPEC)
  REAL(DP)  :: K1_T(0:LSPEC),K2_T(0:LSPEC)
  REAL(DP)  :: K1(0:LSPEC),K2(0:LSPEC)

  REAL(DP)  :: HI(0:ISPEC)         ! FACTOR FOR PSEUDO HENRY OVER ICE
  REAL(DP)  :: K_ICE(0:ISPEC)      ! ADSORPTION CONSTANT (LANGMUIR)
  REAL(DP)  :: H_ADS(0:ISPEC)      ! ADSORPTION ENTHALPIE (LANGMUIR)
  REAL(DP)  :: SCAV_VALUE_I(0:ISPEC)

! 2D SPECIES FIELD
  REAL (DP), PUBLIC, POINTER :: L_SPEC(:,:)         => NULL()
  REAL (DP), PUBLIC, POINTER :: I_SPEC(:,:)         => NULL()


  INTEGER   :: nucl_mode
  INTEGER :: n_resmod = 2

END MODULE MESSY_SCAV_MEM
