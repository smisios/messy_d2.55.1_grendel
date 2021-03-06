MODULE MESSY_CONVECT_DONNER_ADDITIONS

  USE MESSY_MAIN_CONSTANTS_MEM,       ONLY: DP, pie=>PI, RDGAS=>rd, RVGAS=>rv, &
                                            DENS_H2O=>rho_H2O, CP_AIR, GRAV=>g,&
                                            KELVIN=>TMELT
  USE MESSY_CONVECT_DONNER_TYPES_MOD, ONLY: DONNER_INITIALIZED_TYPE,           &
                                            DONNER_SAVE_TYPE, DONNER_RAD_TYPE, &
                                            DONNER_NML_TYPE, DONNER_PARAM_TYPE,&
                                            DONNER_COLUMN_DIAG_TYPE,           &
                                            DONNER_CONV_TYPE, DONNER_CAPE_TYPE

  IMPLICIT NONE
!---------------------------------------------------------------------
!---NAMELIST----

!----------------------------------------------------------------------
!  THE FOLLOWING NML VARIABLES ARE STORED IN A DONNER_NML_TYPE DERIVED
!  TYPE VARIABLE SO THEY MAY BE CONVENIENTLY PASSED TO KERNEL SUBROUTINES
!  AS NEEDED:
!----------------------------------------------------------------------

INTEGER             :: MODEL_LEVELS_IN_SFCBL = 2 
                             ! NUMBER OF LEVELS AT WHICH THE TEMPERATURE 
                             ! AND VAPOR PROFILES ARE NOT ALLOWED TO 
                             ! CHANGE FROM LAG VALUE WHEN CALCULATING THE
                             ! TIME TENDENCY OF CAPE
INTEGER             :: PARCEL_LAUNCH_LEVEL = 2 
                             ! LARGE-SCALE MODEL LEVEL FROM WHICH A PAR-
                             ! CEL IS LAUNCHED TO DETERMINE THE LIFTING 
                             ! CONDENSATION LEVEL (LEVEL 1 NEAREST THE 
                             ! SURFACE)
LOGICAL             :: ALLOW_MESOSCALE_CIRCULATION = .TRUE.
                             ! A MESOSCALE CIRCULATION WILL BE INCLUDED 
                             ! IN THOSE COLUMNS WHICH SATISFY THE 
                             ! REQUIRED CONDITIONS ?
INTEGER             :: DONNER_DEEP_FREQ = 4800 
!INTEGER             :: DONNER_DEEP_FREQ = 1800 
                             ! FREQUENCY OF CALLING DONNER_DEEP [ SEC ]; 
                             ! MUST BE <= 86400 
CHARACTER(LEN=16)   :: CELL_LIQUID_SIZE_TYPE = 'BOWER' 
                             ! CHOOSE EITHER 'INPUT' OR 'BOWER' 
REAL                :: CELL_LIQUID_EFF_DIAM_INPUT = -1.0 
                             ! INPUT CELL DROPLET EFFECTIVE DIAMETER 
                             ! [ MICRONS ];
                             ! NEEDED WHEN CELL_LIQUID_SIZE_TYPE == 
                             ! 'INPUT'
CHARACTER(LEN=16)   :: CELL_ICE_SIZE_TYPE = 'DEFAULT' 
                             ! CHOOSE EITHER 'INPUT' OR 'DEFAULT'
REAL                :: CELL_ICE_GENEFF_DIAM_INPUT = -1.0 
                             ! INPUT CELL ICE GENERALIZED EFFECTIVE DIAM-
                             ! ETER [ MICRONS ]; NEEDED WHEN 
                             ! CELL_ICE_SIZE_TYPE == 'INPUT'
REAL                :: MESO_LIQUID_EFF_DIAM_INPUT = -1.0 
                             ! INPUT MESOSCALE DROPLET EFFECTIVE DIAMETER
                             ! [ MICRONS ]; CURRENTLY NO LIQUID ALLOWED 
                             ! IN MESOSCALE CLOUD
LOGICAL             :: DO_AVERAGE = .FALSE.   
                             ! TIME-AVERAGE DONNER CLOUD PROPERTIES FOR 
                             ! USE BY RADIATION PACKAGE?
CHARACTER(LEN=32)   :: ENTRAINMENT_CONSTANT_SOURCE = 'gate'
                             ! SOURCE OF CLOUD ENTRAINMENT CONSTANTS FOR
                             !  CUMULUS ENSEMBLE; EITHER 'GATE' OR 'KEP'

!----------------------------------------------------------------------
!   THE FOLLOWING NML VARIABLES ARE NOT NEEDED IN ANY KERNEL SUBROUTINES
!   AND SO ARE NOT INCLUDED IN THE DONNER_NML_TYPE VARIABLE:
!----------------------------------------------------------------------

LOGICAL             :: DO_NETCDF_RESTART = .TRUE.
                             ! RESTART FILE WRITTEN IN NETCDF FORMAT 
                             ! (OPTION IS NATIVE MODE) 
                             ! NOTE: CURRENT CODE WILL PRODUCE ONLY
                             ! A NETCDF RESTART; IF NATIVE MODE RESTART 
                             ! IS DESIRED, USER MUST UPDATE THE SOURCE 
                             ! CODE APPROPRIATELY.
LOGICAL             :: WRITE_REDUCED_RESTART_FILE = .FALSE.
                             ! BY SETTING THIS VARIABLE TO .TRUE., THE 
                             ! USER IS ASSERTING THAT THE DONNER DEEP 
                             ! CALCULATION WILL BE MADE ON THE FIRST STEP
                             ! OF ANY JOB READING THE RESTART FILE, SO 
                             ! THAT THOSE VARIABLES NEEDED ON STEPS WHEN
                             ! DONNER_DEEP IS NOT CALCULATED MAY BE 
                             ! OMITTED FROM THE FILE, THUS SAVING ARCHIVE
                             ! SPACE. A CHECK IS PROVIDED IN THE CODE 
                             ! WHICH WILL FORCE THE WRITING OF A FULL 
                             ! FILE (EVEN WITH THIS VARIABLE SET TO 
                             !.TRUE.), IF IT IS DETERMINED THAT THE CUR-
                             ! RENT JOB DOES NOT END JUST PRIOR TO A 
                             ! DONNER STEP. 
                             ! USER SHOULD SET THIS VARIABLE TO .FALSE. 
                             ! IF IT IS KNOWN THAT DONNER_DEEP_FREQ IS TO
                             ! BE CHANGED AT THE NEXT RESTART TO AVOID A
                             ! FATAL ERROR IN THAT JOB.
INTEGER, PARAMETER  :: MAX_PTS = 20
                             ! MAX NUNBER OF DIAGNOSTIC COLUMNS
REAL                :: DIAGNOSTICS_PRESSURE_CUTOFF =  50.E+02   
                             ! COLUMN DATA WILL BE OUTPUT ON MODEL LEVELS
                             ! WITH PRESSURES GREATER THAN THIS VALUE 
                             ! [ PA ]
INTEGER,                        &
    DIMENSION(6)    :: DIAGNOSTICS_START_TIME= (/ 0,0,0,0,0,0 /)
                             ! INTEGER SPECIFICATION OF TIME AT WHICH 
                             ! COLUMN DIAGNOSTICS IS TO BE ACTIVATED
                             ! [YEAR, MONTH, DAY, HOUR, MIN, SEC ]
INTEGER             :: NUM_DIAG_PTS_IJ = 0   
                             ! NUMBER OF DIAGNOSTIC COLUMNS SPECIFIED BY 
                             ! GLOBAL(I,J) COORDINATES
INTEGER             :: NUM_DIAG_PTS_LATLON = 0 
                             ! NUMBER OF DIAGNOSTIC COLUMNS SPECIFIED BY
                             ! LAT-LON COORDINATES 
INTEGER,                       &
 DIMENSION(MAX_PTS) :: I_COORDS_GL = -100
                             ! GLOBAL I COORDINATES FOR IJ DIAGNOSTIC 
                             ! COLUMNS
INTEGER,                        &
 DIMENSION(MAX_PTS) :: J_COORDS_GL = -100
                             ! GLOBAL J COORDINATES FOR IJ DIAGNOSTIC 
                             ! COLUMNS
REAL,                           &
 DIMENSION(MAX_PTS) :: LAT_COORDS_GL = -999.
                             ! LATITUDES FOR LAT-LON  DIAGNOSTIC COLUMNS 
                             ! [ DEGREES, -90. -> 90. ]
REAL,                            &
 DIMENSION(MAX_PTS) :: LON_COORDS_GL = -999.
                             ! LONGITUDES FOR LAT-LON  DIAGNOSTIC COLUMNS
                             ! [ DEGREES, 0. -> 360. ]

NAMELIST / DONNER_DEEP_NML /      &

! CONTAINED IN DONNER_RAD_TYPE VARIABLE:
                            MODEL_LEVELS_IN_SFCBL, &
                            PARCEL_LAUNCH_LEVEL, &
                            ALLOW_MESOSCALE_CIRCULATION, &
                            DONNER_DEEP_FREQ, &
                            CELL_LIQUID_SIZE_TYPE,   &
                            CELL_LIQUID_EFF_DIAM_INPUT, &
                            CELL_ICE_SIZE_TYPE, &
                            CELL_ICE_GENEFF_DIAM_INPUT, &
                            MESO_LIQUID_EFF_DIAM_INPUT, &
                            DO_AVERAGE,  &
                            ENTRAINMENT_CONSTANT_SOURCE, &

! NOT CONTAINED IN DONNER_RAD_TYPE VARIABLE:
                            DO_NETCDF_RESTART, &
                            WRITE_REDUCED_RESTART_FILE, &
                            DIAGNOSTICS_PRESSURE_CUTOFF, &
                            DIAGNOSTICS_START_TIME, &
                            NUM_DIAG_PTS_IJ, NUM_DIAG_PTS_LATLON, &
                            I_COORDS_GL, J_COORDS_GL, &
                            LAT_COORDS_GL, LON_COORDS_GL


!--------------------------------------------------------------------
!--- PUBLIC DATA ----------




!--------------------------------------------------------------------
!----PRIVATE DATA-----------



!---------------------------------------------------------------------
!  PARAMETERS STORED IN THE DONNER_PARAM DERIVED TYPE VARIABLE TO FACILI-
!  TATE PASSAGE TO KERNEL SUBROUTINES:
!
real, parameter :: kappa = 2./7.
real, parameter :: seconds_per_day = 86400.
real, parameter :: hlv = 2.500e+06
real, parameter :: hlf = 3.34e+05
real, parameter :: hls = hlv + hlf


REAL,                          &
  PARAMETER                     &
             ::  CP_VAPOR= 4.0*RVGAS  
                       ! SPECIFIC HEAT OF WATER VAPOR AT CONSTANT PRES-
                       ! SURE [ J/(KG K) ]
REAL,                             &
  PARAMETER                    &
             ::  D622 = RDGAS/RVGAS 
                       ! RATIO OF MOLECULAR WEIGHTS OF WATER VAPOR AND 
                       ! DRY AIR
REAL,                        &
  PARAMETER                 &
             ::  D608 = (RVGAS/RDGAS - 1.0)
                        ! FACTOR IN VIRTUAL TEMPERATURE DEFINITION
INTEGER,                   &
  PARAMETER                  &
             ::  KPAR=7         
                        ! NUMBER OF MEMBERS IN CUMULUS ENSEMBLE
INTEGER,                   &
  PARAMETER                  &
             ::  NLEV_HIRES=100       
                        ! NUMBER OF LEVELS IN CLOUD MODEL
REAL,                       &
  PARAMETER                 &
             ::  PDEEP_CV = 500.E02 
                        ! MINIMUM PRESSURE DIFFERENCE BETWEEN LEVEL OF 
                        ! FREE CONVECTION AND LEVEL OF ZERO BUOYANCY 
                        ! NEEDED FOR DEEP CONVECTION TO OCCUR [ PA ].
REAL,                       &
  PARAMETER                 &
             ::  CDEEP_CV = 100.   
                        ! MAXIMUM VALUE OF CONVECTIVE INHIBITION (J/KG) 
                        ! THAT ALLOWS CONVECTION. VALUE OF 10 SUGGESTED 
                        ! BY TABLE 2 IN THOMPSON ET AL. (1979, JAS).
REAL,                         &
  PARAMETER,                   &
  DIMENSION(KPAR)              &
             ::  ARAT  =  (/  1.0, 0.26, 0.35, 0.32, 0.3, 0.54, 0.66/)
                        ! RATIO AT CLOUD BASE OF THE FRACTIONAL AREA OF 
                        ! ENSEMBLE MEMBER I RELATIVE TO ENSEMBLE MEMBER 
                        ! 1. (TAKEN FROM GATE DATA).
REAL,                       &
  PARAMETER                 &
             ::  MAX_ENTRAINMENT_CONSTANT_GATE = 0.0915
                        ! ENTRAINMENT CONSTANT BASED ON GATE DATA FOR 
                        ! MOST ENTRAINING ENSEMBLE MEMBER
REAL,                       &
  PARAMETER                 &
             ::  MAX_ENTRAINMENT_CONSTANT_KEP  = 0.0915
                        ! ENTRAINMENT CONSTANT BASED ON KEP DATA FOR MOST
                        ! ENTRAINING ENSEMBLE MEMBER
REAL,                       &
  PARAMETER,                 &
  DIMENSION(KPAR)              &
             ::  ENSEMBLE_ENTRAIN_FACTORS_GATE =  (/ 1.0, 1.3, 1.8, 2.5,&
                                                     3.3, 4.5, 10. /)
                        ! RATIO OF ENTRAINMENT CONSTANT BETWEEN ENSEMBLE
                        ! MEMBER 1 AND ENSEMBLE MEMBER I FOR GATE-BASED
                        ! ENSEMBLE
REAL,                       &
  PARAMETER,                 &
  DIMENSION(KPAR)              &
             ::  ENSEMBLE_ENTRAIN_FACTORS_KEP  = (/ 1.0, 1.22, 1.56,  &
                                                    2.05, 2.6, 3.21,   &
                                                    7.84 /)
                        ! RATIO OF ENTRAINMENT CONSTANT BETWEEN ENSEMBLE
                        ! MEMBER 1 AND ENSEMBLE MEMBER I FOR KEP-BASED
                        ! ENSEMBLE
REAL,                       &
  PARAMETER                 &
             ::  CLD_BASE_VERT_VEL = 0.5                          
                        ! VERTICAL VELOCITY ASSUMED PRESENT AT CLOUD BASE
                        ! [ M / SEC ]
REAL,                       &
  PARAMETER                 &
             ::  PSTOP = 40.0E02   
                        ! LOWEST POSSIBLE PRESSURE TO WHICH A CLOUD MAY 
                        ! EXTEND IN THE CLOUD MODEL [ PA ]
REAL,                       &
  PARAMETER                 &
             ::  PARCEL_DP = -1.0E02
                        ! PRESSURE INCREMENT USED FOR PARCEL CALCULATIONS
                        ! [ PA ]
REAL,                       &
  PARAMETER                 &
             ::  UPPER_LIMIT_FOR_LFC = 500.E02 
                        ! LOWEST PRESSURE ALLOWED FOR LEVEL OF FREE CONV-
                        ! ECTION [ PA ]
REAL,                       &
  PARAMETER                 &
             ::  DP_OF_CLOUD_MODEL = -10.E02
                        ! PRESSURE THICKNESS (PA) OF THE LAYERS IN THE
                        ! DONNER PARAMETERIZATION'S CLOUD MODEL.
REAL,                       &
  PARAMETER                 &
             ::  CLOUD_BASE_RADIUS = 1000.
                        ! RADIUS ASSUMED FOR CLOUD ENSEMBLE MEMBER #1 AT
                        ! CLOUD BASE [ M ]
REAL,                       &
  PARAMETER                 &
             ::  WDET = .1   
                        ! VERTICAL VELOCITY AT WHICH DETRAINMENT FROM THE
                        ! CLOUDS BEGINS [ M/S ]
REAL,                       &
  PARAMETER                 &
             ::  RBOUND = 0.01    
                        ! VALUE OF CUMULUS RADIUS AT WHICH CLOUD EFFECT-
                        ! IVELY DISAPPEARS AND CLOUD MODEL CALCULATION 
                        ! STOPS [ M ]
REAL,                       &
  PARAMETER                 &
             ::  WBOUND = 0.01  
                        ! VALUE OF CUMULUS VERTICAL VELOCITY AT WHICH 
                        ! CLOUD MODEL CALCULATION STOPS [ M / SEC ]
REAL,                       &
  PARAMETER                 &
             ::  FREEZE_FRACTION = 0.52
                        ! FRACTION OF LIQUID IN CLOUD UPDRAFT WHICH MAY 
                        ! BE FROZEN. (LEARY AND HOUZE (JAS,1980)) 
                        ! [ DIMENSIONLESS ]
REAL,                       &
  PARAMETER                 &
             ::  VIRT_MASS_CO = 0.5
                        ! VIRTUAL MASS COEFFICIENT [ DIMENSIONLESS ]
REAL,                       &
  PARAMETER                 &
             ::  PDEEP_MC = 200.E02 
                        ! PRESSURE THICKNESS [ PA ] REQUIRED FOR MESO-
                        ! SCALE CIRCULATION. IT REFERS TO THE LEAST
                        ! PENETRATIVE ENSEMBLE MEMBER. FOR THIS CHECK 
                        ! TO FUNCTION PROPERLY, THE ENTRAINMENT COEFFIC-
                        ! IENT IN CLOUD_MODEL FOR KOU=1 MUST BE THE 
                        ! LARGEST ENTRAINMENT COEFFICIENT.
REAL,                       &
  PARAMETER                 &
             ::  TR_INSERT_TIME = 0.0
                        ! FRACTIONAL POINT (BASED ON MASS INCREASE) 
                        ! DURING A TIMESTEP AT WHICH AN ENTRAINING PARCEL
                        ! TAKES ON INTERNALLY-GENERATED TRACER 
                        ! [ DIMENSIONLESS, VALUE BETWEEN 0.0 AND 1.0 ]
REAL,                       &
  PARAMETER                 &
             ::  AUTOCONV_RATE = 1.0E-03
                        ! RATE OF AUTOCONVERSION OF CLOUD TO RAINWATER 
                        ! [ SEC**(-1) ]
REAL,                       &
  PARAMETER                 &
             ::  AUTOCONV_THRESHOLD =  0.5    
                        ! THRESHOLD OF CLOUD WATER AT WHICH AUTOCONVER-
                        ! SION OF CLOUD TO RAINWATER BEGINS  [ G / M**3 ]
REAL,                       &
  PARAMETER                 &
             ::  TFRE = 258.  
                        ! TEMPERATURE AT WHICH CLOUD LIQUID BEGINS TO 
                        ! FREEZE [ DEG K ]
REAL,                       &
  PARAMETER                 &
             ::  DFRE = 10.   
                        ! RANGE OF TEMPERATURE BETWEEN THE ONSET AND 
                        ! COMPLETION OF FREEZING  [ DEG K ]
REAL,                       &
  PARAMETER                 &
             ::  EVAP_IN_DOWNDRAFTS = 0.25
                        ! FRACTION OF CONDENSATE AVAILABLE TO THE MESO-
                        ! SCALE WHICH IS EVAPORATED IN CONVECTIVE DOWN-
                        ! DRAFTS (FROM LEARY AND LOUZE, 1980)
                        ! [ DIMENSIONLESS ]
REAL,                       &
  PARAMETER                 &
             ::  EVAP_IN_ENVIRON = 0.13
                        ! FRACTION OF CONDENSATE AVAILABLE TO THE MESO-
                        ! SCALE WHICH IS EVAPORATED IN THE CELL ENVIRON-
                        ! MENT (FROM LEARY AND LOUZE, 1980)
                        ! [ DIMENSIONLESS ]
REAL,                       &
  PARAMETER                 &
             ::  ENTRAINED_INTO_MESO = 0.62
                        ! FRACTION OF CONDENSATE AVAILABLE TO THE MESO-
                        ! SCALE WHICH IS ENTRAINED INTO THE MESOSCALE 
                        ! CIRCULATION (FROM LEARY AND LOUZE, 1980)
                        ! [ DIMENSIONLESS ]
REAL,                       &
  PARAMETER                 &
             ::  UPPER_LIMIT_FOR_LCL = 500.0E02
                        ! LOWEST PRESSURE ALLOWABLE FOR LIFTING CONDENS-
                        ! ATION LEVEL; DEEP CONVECTION WILL NOT BE PRES-
                        ! ENT IF LCL NOT REACHED BEFORE THIS PRESSURE 
                        ! [ PA ]
INTEGER,                   &
  PARAMETER         &
             ::  ISTART = 1    
                        ! INDEX OF LEVEL IN CAPE GRID FROM WHICH THE 
                        ! PARCEL ORIGINATES FOR THE CAPE CALCULATIONS
REAL,                       &
  PARAMETER                 &
             ::  TMIN = 154.       
                        ! CAPE CALCULATIONS ARE TERMINATED WHEN PARCEL 
                        ! TEMPERATURE GOES BELOW TMIN [ DEG K ]
REAL,                       &
  PARAMETER                 &
             ::  ANVIL_PRECIP_EFFICIENCY = 0.5
                        ! FRACTION OF TOTAL CONDENSATE IN ANVIL (TRANSFER
                        ! FROM CELL PLUS IN SITU CONDENSATION) THAT 
                        ! PRECIPITATES OUT (FROM LEARY AND LOUZE, 1980)
                        ! [ DIMENSIONLESS ]
REAL,                       &
  PARAMETER                 &
             ::  MESO_LIFETIME = 64800.
                        ! ASSUMED LIFETIME OF MESOSCALE CIRCULATION 
                        ! (FROM LEARY AND LOUZE, 1980) [ SEC ]
REAL,                       &
  PARAMETER                 &
             ::  MESO_REF_OMEGA = -0.463
                        ! ASSUMED REFERENCE OMEGA FOR MESOSCALE UPDRAFT 
                        ! (FROM LEARY AND LOUZE, 1980) [ PA / SEC ]
REAL,                       &
  PARAMETER                 &
             ::  TPRIME_MESO_UPDRFT = 1.0    
                        ! ASSUMED TEMPERATURE EXCESS OF MESOSCALE UPDRAFT
                        ! OVER ITS ENVIRONMENT [ DEG K ]
REAL,                       &
  PARAMETER                 &
             ::  MESO_SEP = 200.0E+02
                        ! PRESSURE SEPARATION BETWEEN BASE OF MESOSCALE
                        ! UPDRAFT AND TOP OF MESOSCALE DOWNDRAFT [ PA ]
REAL,                       &
  PARAMETER                 &
             ::  REF_PRESS = 1.0E05
                        ! REFERENCE PRESSURE USED IN CALCULATION OF EXNER
                        ! FUMCTION [ PA ]
REAL,                       &
  PARAMETER                 &
             ::  MESO_DOWN_EVAP_FRACTION = 0.4
                        ! FRACTION OF TOTAL ANVIL CONDENSATE ASSUMED 
                        ! EVAPORATED IN THE MESOSCALE DOWNDRAFT
                        ! [ FRACTION ]
REAL,                       &
  PARAMETER                 &
             ::  MESO_UP_EVAP_FRACTION = 0.1
                        ! FRACTION OF TOTAL ANVIL CONDENSATE ASSUMED 
                        ! EVAPORATED IN OUTFLOW FROM THE MESOSCALE 
                        ! UPDRAFT [ FRACTION ]    
REAL,                       &
  PARAMETER                 &
             ::  R_CONV_LAND  = 10.0    
                        ! ASSUMED CONVECTIVE CLOUD DROPLET RADIUS OVER 
                        ! LAND [ MICRONS ]   
REAL,                       &
  PARAMETER                 &
             ::  R_CONV_OCEAN = 16.0  
                        ! ASSUMED CONVECTIVE CLOUD DROPLET RADIUS OVER 
                        ! OCEAN [ MICRONS ]   
REAL,                       &
  PARAMETER                 &
             ::  N_LAND = 600*1.0E6 
                        ! ASSUMED DROPLET NUMBER CONC OVER LAND (M**-3)
REAL,                       &
  PARAMETER                 &
             ::  N_OCEAN = 150*1.0E6 
                        ! ASSUMED DROPLET NUMBER CONC OVER OCEAN (M**-3)
REAL,                       &
  PARAMETER                 &
             ::  DELZ_LAND = 500.0   
                        ! ASSUMED CLOUD DEPTH OVER LAND (M) 
REAL,                       &
  PARAMETER                 &
             ::  DELZ_OCEAN = 1500.0   
                        ! ASSUMED CLOUD DEPTH OVER OCEAN (M)
REAL,                       &
  PARAMETER                 &
             ::  CELL_LIQUID_EFF_DIAM_DEF = 15.0    
                        ! DEFAULT CELL LIQUID EFF DIAMETER [ MICRONS ]
REAL,                       &
  PARAMETER                 &
             ::  CELL_ICE_GENEFF_DIAM_DEF = 18.6   
                        ! DEFAULT CELL ICE GENERALIZED EFFECTIVE DIAMETER
                        ! [ MICRONS ]
INTEGER,                   &
  PARAMETER         &
             ::  ANVIL_LEVELS = 6  
                        ! NUMBER OF LEVELS ASSUMED TO BE IN ANVIL CLOUDS
REAL,                       &
  PARAMETER,                &
  DIMENSION(ANVIL_LEVELS)   &
             ::  DGEICE  = (/ 38.5, 30.72, 28.28, 25.62, 24.8, 13.3 /)
                        ! GENERALIZED EFFECTIVE SIZE OF HEXAGONAL ICE 
                        ! CRYSTALS, DEFINED AS IN FU (1996, J. CLIM.) 
                        ! VALUES FROM TABLE 2 OF MCFARQUHAR ET AL. 
                        ! (1999, JGR) ARE AVERAGED OVER ALL GRID BOXES 
                        ! FOR WHICH D_GE IS DEFINED FOR ALL ALTITUDES 
                        ! BETWEEN 9.9 AND 13.2 KM. INDEX 1 AT BOTTOM OF 
                        ! ANVIL
REAL,                       &
  PARAMETER,                &
  DIMENSION(ANVIL_LEVELS)   &
             ::  RELHT  =  (/0.0, 0.3, 0.45, 0.64, 0.76, 1.0/)
                        ! DISTANCE FROM ANVIL BASE, NORMALIZED BY TOTAL 
                        ! ANVIL THICKNESS. FROM TABLE 2 OF MCFARQUHAR ET
                        ! AL. (1999, JGR) FOR GRID BOXES WITH DATA 
                        ! BETWEEN 9.9 AND 13.2 KM. INDEX 1 AT ANVIL 
                        ! BOTTOM


!--------------------------------------------------------------------
!   LIST OF NATIVE MODE RESTART VERSIONS USABLE BY THIS MODULE:
!
!   NOTE: NONE OF THE EARLIER VERSIONS OF RESTART FILES CAN BE USED TO
!         INITIATE AN EXPERIMENT WITH THIS CODE VERSION DUE TO A CHANGE 
!         IN THE CALCULATION ALGORITHM. EXPERIMENTS BEGUN WITH THIS CODE
!         MUST BE COLDSTARTED, OR USE A NATIVE MODE RESTART FILE GENER-
!         ATED BY AN EXPERIMENT USING THIS CODE VERSION (RESTART VERSION
!         #8), OR A NETCDF RESTART FILE.
!          
!   VERSION 8 HAS THE LAG TEMP, VAPOR AND PRESSURE FIELDS NEEDED TO CAL-
!             CULATE THE LAG TIME VALUE OF CAPE. TEMPBL AND RATPBL
!             REMOVED. 
!
!   VERSION 9 IS RESERVED FOR THE NATIVE MODE RESTART FILE VERSION COR-
!             RESPONDING TO THE CURRENT NETCDF RESTART FILE. IT IS UP TO 
!             THE USER TO GENERATE THE CODE NEEDED TO READ AND WRITE THIS
!             VERSION, IF NEEDED, USING THE SUBROUTINES READ_RESTART AND 
!             WRITE_RESTART THAT ARE PROVIDED AS STARTING POINTS, SINCE 
!             ONLY NETCDF RESTARTS ARE CURRENTLY SUPPORTED.
!

INTEGER, DIMENSION(2)  :: RESTART_VERSIONS = (/ 8, 9 /)


!--------------------------------------------------------------------
!   VARIABLES ASSOCIATED WITH NETCDF DIAGNOSTIC OUTPUT FROM THIS MODULE:
!
!   ID_XXXX         INDICES ASSOCIATED WITH EACH POTENTIAL NETCDF 
!                   DIAGNOSTIC FIELD:
!   MISSING VALUE   VALUE USED BY NETCDF ROUTINES IF DATA NOT PRESENT
!   MOD_NAME        MODULE NAME ASSOCIATED WITH THESE DIAGNOSTICS; USED
!                   TO CONNECT THESE DIAGNOSTICS TO THE DIAG_TABLE
!

INTEGER    :: ID_CEMETF_DEEP, ID_CEEFC_DEEP, ID_CECON_DEEP, &
              ID_CEMFC_DEEP, ID_CEMEMF_DEEP, ID_CEMEMF_MOD_DEEP, &
              ID_CUAL_DEEP, ID_FRE_DEEP, ID_ELT_DEEP, &
              ID_CMUS_DEEP, ID_ECDS_DEEP, ID_ECES_DEEP, &
              ID_EMDS_DEEP, ID_EMES_DEEP, ID_QMES_DEEP,&
              ID_WMPS_DEEP, ID_WMMS_DEEP, ID_TMES_DEEP,&
              ID_DMEML_DEEP, ID_UCEML_DEEP, ID_UMEML_DEEP, &
              ID_XICE_DEEP,  ID_DGEICE_DEEP, ID_DGELIQ_DEEP,  &
              ID_XLIQ_DEEP,    &
              ID_CUQI_DEEP, ID_CUQL_DEEP, &
              ID_PLCL_DEEP, ID_PLFC_DEEP, ID_PLZB_DEEP, &
              ID_XCAPE_DEEP, ID_COIN_DEEP,  &
              ID_DCAPE_DEEP, ID_QINT_DEEP, ID_A1_DEEP, &
              ID_AMAX_DEEP, ID_AMOS_DEEP, &
              ID_TPREA1_DEEP, ID_AMPTA1_DEEP, &
              ID_OMINT_DEEP, ID_RCOA1_DEEP, ID_DETMFL_DEEP

INTEGER, DIMENSION(:), ALLOCATABLE :: ID_QTREN1, ID_QTMES1, &
                                      ID_WTP1, ID_QTCEME
INTEGER, DIMENSION(:), ALLOCATABLE :: ID_QTREN1_COL, ID_QTMES1_COL, &
                                      ID_WTP1_COL, ID_QTCEME_COL

REAL              :: MISSING_VALUE = -999.
CHARACTER(LEN=16) :: MOD_NAME = 'DONNER_DEEP'


!!$!--------------------------------------------------------------------
!!$!   VARIABLES FOR COLUMN DIAGNOSTICS OPTION
!!$!
!!$!   ARRAYS CONTAINING INFORMATION FOR ALL REQUESTED DIAGNOSTIC COLUMNS
!!$!   (1:NUM_DIAG_PTS):
!!$!    COL_DIAG_UNIT         UNIT NUMBERS FOR EACH COLUMN'S OUTPUT FILE 
!!$!    COL_DIAG_LON          EACH COLUMN'S LONGITUDE 
!!$!                          [ DEGREES, 0 < LON < 360 ]
!!$!    COL_DIAG_LAT          EACH COLUMN'S LATITUDE 
!!$!                          [DEGREES, -90 < LAT < 90 ]
!!$!    COL_DIAG_J            EACH COLUMN'S J INDEX (PROCESSOR COORDINATES)
!!$!    COL_DIAG_I            EACH COLUMN'S I INDEX (PROCESSOR COORDINATES) 
!!$!
!!$!    TIME_COL_DIAGNOSTICS  TIME IN MODEL SIMULATION AT WHICH TO ACTIVATE
!!$!                          COLUMN DIAGNOSTICS 
!!$!
!!$
!!$INTEGER, DIMENSION(:), ALLOCATABLE :: COL_DIAG_UNIT
!!$REAL   , DIMENSION(:), ALLOCATABLE :: COL_DIAG_LON, COL_DIAG_LAT   
!!$INTEGER, DIMENSION(:), ALLOCATABLE :: COL_DIAG_J, COL_DIAG_I        
!!$TYPE(TIME_TYPE)                    :: TIME_COL_DIAGNOSTICS  

!---------------------------------------------------------------------
!    DERIVED TYPE VARIABLES PRESENT FOR DURATION OF JOB:
!    (SEE DONNER_TYPES.H FOR DOCUMENTATION OF THEIR CONTENTS)
!

TYPE(DONNER_PARAM_TYPE), SAVE       :: PARAM
TYPE(DONNER_COLUMN_DIAG_TYPE), SAVE :: COL_DIAG
TYPE(DONNER_NML_TYPE), SAVE         :: NML
TYPE(DONNER_SAVE_TYPE), SAVE        :: DON_SAVE
TYPE(DONNER_INITIALIZED_TYPE), SAVE :: INITIALIZED


!-----------------------------------------------------------------------
!   MISCELLANEOUS VARIABLES
!
!     MODULE_IS_INITIALIZED       MODULE HAS BEEN INITIALIZED ?
!

LOGICAL :: MODULE_IS_INITIALIZED = .FALSE. 

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


  PRIVATE

  PUBLIC :: LOOKUP_ES_K, DONNER_INIT
  PUBLIC :: PARAM, COL_DIAG, NML, DON_SAVE, INITIALIZED
  PUBLIC :: NLEV_HIRES, MODULE_IS_INITIALIZED
!---------------------------------------
CONTAINS
!---------------------------------------

  SUBROUTINE DONNER_INIT(idf, jdf, nlev)
!  SUBROUTINE DONNER_INIT(lonb, latb, pref, axes, Time,  &
!                         tracers_in_donner)

!---------------------------------------------------------------------
!    donner_deep_init is the constructor for donner_deep_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!!$real,            dimension(:), intent(in)   :: lonb, latb, pref
!!$integer,         dimension(4), intent(in)   :: axes
!!$real,                          intent(in)   :: Time
!!$!type(time_type),               intent(in)   :: Time
!!$logical,         dimension(:), intent(in)   :: tracers_in_donner
    INTEGER, INTENT(IN)             :: idf, jdf, nlev

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      lonb         array of model longitudes on cell boundaries  
!                   [ radians ]
!      latb         array of model latitudes on cell boundaries
!                   [ radians ]
!      pref         array of reference pressures at full levels (plus 
!                   surface value at nlev+1), based on 1013.25 hPa pstar
!                   [ Pa ]
!      axes         data axes for diagnostics
!      Time         current time [ time_type ]
!      tracers_in_donner 
!                   logical array indicating which of the activated 
!                   tracers are to be transported by donner_deep_mod
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables:

      integer                             :: unit, ierr, io
      integer                             :: ntracers
      integer                             :: secs, days
      logical, dimension(1)               :: do_column_diagnostics
      integer                             :: k, n, nn
  
!-------------------------------------------------------------------
!  local variables:
!
!     unit                   unit number for nml file
!     ierr                   error return flag
!     io                     error return code
!     idf                    number of columns in the x dimension on the
!                            processors domain
!     jdf                    number of columns in the y dimension on the
!                            processors domain
!     nlev                   number of model layers 
!     ntracers               number of tracers to be transported by
!                            the donner deep convection parameterization
!     secs                   seconds component of time_type variable Time
!     days                   days component of time_type variable Time
!     do_column_diagnostics  logical array indicating which latitude rows
!                            in the processor domain contain diagnostic
!                            columns
!     k, n                   do-loop indices
!     nn                     counter of tracers transported by 
!                            donner_deep_mod
!                         
!-------------------------------------------------------------------
!!$!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$!
!!$!    1. READ NAMELIST AND WRITE IT TO LOG FILE.
!!$!
!!$!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$
!!$!---------------------------------------------------------------------
!!$!    read namelist.
!!$!---------------------------------------------------------------------
!!$      if (file_exist('input.nml')) then
!!$        unit =  open_namelist_file ()
!!$        ierr=1; do while (ierr /= 0)
!!$        read (unit, nml=donner_deep_nml, iostat=io, end=10)
!!$        ierr = check_nml_error (io, 'donner_deep_nml')
!!$        enddo
!!$10      call close_file (unit)
!!$      endif
!!$
!!$!---------------------------------------------------------------------
!!$!    write version number and namelist to logfile.
!!$!---------------------------------------------------------------------
!!$      call write_version_number (version, tagname)
!!$      if (mpp_pe() == mpp_root_pe() )    &
!!$                                 write (stdlog(), nml=donner_deep_nml)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    2. DO CONSISTENCY / VALIDITY TESTS ON NML AND PARAMETER VARIABLES.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$
!!$!---------------------------------------------------------------------
!!$!    check for a valid value of donner_deep_freq. 
!!$!---------------------------------------------------------------------
!!$      if (donner_deep_freq > 86400) then
!!$        call error_mesg ( 'donner_deep_mod', 'donner_deep_init: &
!!$         & donner convection must be called at least once per day', &
!!$                                                           FATAL)
!!$      else if (donner_deep_freq <= 0) then
!!$        call error_mesg ( 'donner_deep_mod', 'donner_deep_init: &
!!$          & a positive value must be assigned to donner_deep_freq', &
!!$                                                            FATAL)
!!$      endif
!!$
!!$!---------------------------------------------------------------------
!!$!    check for valid value of entrainment_constant_source.
!!$!---------------------------------------------------------------------
!!$      if (trim(entrainment_constant_source) == 'gate' .or. &
!!$          trim(entrainment_constant_source) == 'kep' ) then
!!$      else
!!$        call error_mesg ('donner_deep_mod', 'donner_deep_init: &
!!$        & invalid string for nml variable entrainment_constant_source', &
!!$                                                             FATAL)
!!$      endif
!!$
!!$!---------------------------------------------------------------------
!!$!    test that PSTOP is smaller than UPPER_LIMIT_FOR_LFC.
!!$!---------------------------------------------------------------------
!!$      if (PSTOP > UPPER_LIMIT_FOR_LFC) then
!!$        call error_mesg ('donner_deep_mod', 'donner_deep_init: &
!!$           & pstop must be above the upper limit of &
!!$                                &the level of free convection', FATAL)
!!$      endif
!!$
!!$!---------------------------------------------------------------------
!!$!    make sure the defined ice size is acceptable to the radiative 
!!$!    properties parameterizations.
!!$!---------------------------------------------------------------------
!!$      if (CELL_ICE_GENEFF_DIAM_DEF < 18.6) then
!!$        call error_mesg ('donner_deep_mod', 'donner_deep_init: &
!!$         & cell_ice_geneff_diam_def must be >= 18.6 microns', FATAL)
!!$      endif
!!$
!!$!---------------------------------------------------------------------
!!$!    test that cell_liquid_size_type has been validly specified, and if
!!$!    it is specified as 'input', an appropriate input value has been
!!$!    supplied.
!!$!---------------------------------------------------------------------
!!$      if (trim(cell_liquid_size_type) == 'input') then
!!$        Initialized%do_input_cell_liquid_size = .true.
!!$        Initialized%do_bower_cell_liquid_size = .false.
!!$        if (cell_liquid_eff_diam_input < 0.0) then
!!$          call error_mesg ('donner_deep_mod', 'donner_deep_init: &
!!$            & cell liquid size must be input, but no value supplied', &
!!$                                                               FATAL)
!!$        endif
!!$      else if (trim(cell_liquid_size_type) == 'bower') then
!!$        Initialized%do_input_cell_liquid_size = .false.
!!$        Initialized%do_bower_cell_liquid_size = .true.
!!$      else
!!$        call error_mesg ( 'donner_deep_mod', 'donner_deep_init: &
!!$           & cell_liquid_size_type must be either input or bower', &
!!$                                                                FATAL)
!!$      endif
!!$
!!$!---------------------------------------------------------------------
!!$!    test that cell_ice_size_type has been validly specified, and if
!!$!    specified as 'input', that cell_ice_geneff_diam_input has also 
!!$!    been appropriately defined.
!!$!---------------------------------------------------------------------
!!$      if (trim(cell_ice_size_type) == 'input') then
!!$        Initialized%do_input_cell_ice_size = .true.
!!$        Initialized%do_default_cell_ice_size = .false.
!!$        if (cell_ice_geneff_diam_input <= 0.0) then
!!$          call error_mesg ('donner_deep_mod', 'donner_deep_init: &
!!$               & must define a nonnegative generalized effective '//&
!!$                'diameter for ice when cell_ice_size_type is input', &
!!$                                                                 FATAL)
!!$        endif
!!$        if (cell_ice_geneff_diam_input < 18.6) then
!!$          cell_ice_geneff_diam_input = 18.6
!!$          call error_mesg ('donner_deep_mod', 'donner_deep_init: &
!!$               & resetting cell_ice_geneff_diam_input to 18.6 microns,&
!!$               & the lowest acceptable value for the radiation&
!!$               & parameterization', NOTE)                        
!!$        endif
!!$      else if (trim(cell_ice_size_type) == 'default') then
!!$        Initialized%do_input_cell_ice_size = .false.
!!$        Initialized%do_default_cell_ice_size = .true.
!!$      else
!!$        call error_mesg ( 'donner_deep_init', 'donner_deep_init: &
!!$             & cell_ice_size_type must be input or default',  FATAL)
!!$      endif


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    3. PROCESS TRACERS THAT ARE TO BE TRANSPORTED BY THE DONNER DEEP
!       CONVECTION PARAMETERIZATION.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    determine how many tracers are to be transported by donner_deep 
!    convection. allocate arrays to contain their names and units for use
!    with diagnostics and restarts. define a logical variable indicating
!    if any tracers are to be so transported. obtain the tracer names and
!    units.
!---------------------------------------------------------------------
!!$      ntracers = count(tracers_in_donner)
!!$      allocate ( Don_save%tracername   (ntracers) )
!!$      allocate ( Don_save%tracer_units (ntracers) )
!!$      if (ntracers > 0) then
!!$        Initialized%do_donner_tracer = .true.
!!$        nn = 1
!!$        do n=1,size(tracers_in_donner(:))
!!$          if (tracers_in_donner(n)) then
!!$            call get_tracer_names (MODEL_ATMOS, n,  &
!!$                                   name = Don_save%tracername(nn), &
!!$                                   units = Don_save%tracer_units(nn))
!!$            nn = nn + 1
!!$          endif
!!$        end do
!!$      else
!!$        Initialized%do_donner_tracer = .false.
!!$      endif
      ntracers = 0
      Initialized%do_donner_tracer = .false.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    4. DEFINE PROCESSOR DIMENSIONS AND ALLOCATE SPACE FOR MODULE 
!       VARIABLES.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-------------------------------------------------------------------
!    define the grid dimensions. idf and jdf are the (i,j) dimensions of
!    the domain on this processor, nlev is the number of model layers.
!-------------------------------------------------------------------
!!$      nlev = size(pref(:)) - 1
!!$      idf  = size(lonb(:)) - 1
!!$      jdf  = size(latb(:)) - 1

!---------------------------------------------------------------------
!    initialize the points processed counter. define the total number 
!    of columns present on the processor. 
!---------------------------------------------------------------------
      Initialized%pts_processed_conv   = 0
      Initialized%total_pts = idf*jdf

!--------------------------------------------------------------------
!    allocate module variables that will be saved across timesteps.
!    these are stored in the derived-type variable Don_save. see 
!    donner_types.h for description of these variables.
!--------------------------------------------------------------------
      allocate ( Don_save%cemetf             (idf, jdf, nlev ) )
      allocate ( Don_save%lag_temp           (idf, jdf, nlev ) )
      allocate ( Don_save%lag_vapor          (idf, jdf, nlev ) )
      allocate ( Don_save%lag_press          (idf, jdf, nlev ) )
      allocate ( Don_save%cememf             (idf, jdf, nlev ) )
      allocate ( Don_save%mass_flux          (idf, jdf, nlev ) )
      allocate ( Don_save%cell_up_mass_flux  (idf, jdf, nlev+1 ) )
      allocate ( Don_save%det_mass_flux      (idf, jdf, nlev ) )
      allocate ( Don_save%dql_strat          (idf, jdf, nlev ) )
      allocate ( Don_save%dqi_strat          (idf, jdf, nlev ) )
      allocate ( Don_save%dqa_strat          (idf, jdf, nlev ) )
      allocate ( Don_save%humidity_area      (idf, jdf, nlev ) )
      allocate ( Don_save%humidity_ratio     (idf, jdf, nlev ) )
      allocate ( Don_save%tracer_tends       (idf, jdf, nlev, ntracers) )
      allocate ( Don_save%parcel_disp        (idf, jdf ) )
      allocate ( Don_save%tprea1             (idf, jdf ) )


!!$!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$!
!!$!    4. INITIALIZE THE NETCDF OUTPUT VARIABLES.
!!$!
!!$!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$
!!$!--------------------------------------------------------------------
!!$!    activate the netcdf diagnostic fields.
!!$!-------------------------------------------------------------------
!!$      call register_fields (Time, axes)
!!$
!!$
!!$!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$!
!!$!    5. PROCESS THE RESTART FILE.
!!$!
!!$!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$
!!$!--------------------------------------------------------------------
!!$!    if a netcdf restart file is present, call read_restart_nc to read 
!!$!    it.
!!$!--------------------------------------------------------------------
!!$      if (file_exist ('INPUT/donner_deep.res.nc') ) then
!!$        Initialized%coldstart= .false.
!!$        call read_restart_nc (ntracers)
!!$
!!$!--------------------------------------------------------------------
!!$!    if a native mode restart file is present, call read_restart 
!!$!    to read it.
!!$!--------------------------------------------------------------------
!!$      else if (file_exist ('INPUT/donner_deep.res') ) then
!!$        Initialized%coldstart= .false.
!!$        call read_restart (ntracers, Time)
!!$
!!$!--------------------------------------------------------------------
!!$!    if no restart file is present, call subroutine process_coldstart
!!$!    to define the needed variables.
!!$!--------------------------------------------------------------------
!!$      else
!!$        call process_coldstart (Time)
!!$      endif

!!$!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$!
!!$!    6. INITIALIZE VARIABLES NEEDED FOR COLUMN_DIAGNOSTICS_MOD OUTPUT.
!!$!
!!$!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$
!!$!---------------------------------------------------------------------
!!$!    define the total number of columns for which diagnostics
!!$!    are desired.
!!$!---------------------------------------------------------------------
!!$      Col_diag%num_diag_pts = num_diag_pts_ij + num_diag_pts_latlon
!!$
!!$!---------------------------------------------------------------------
!!$!    initialize the value of the k index associated with diagnostics
!!$!    cutoff.
!!$!---------------------------------------------------------------------
!!$      Col_diag%kstart = -99
!!$
!!$!---------------------------------------------------------------------
!!$!    if any diagnostics are requested, perform various consistency
!!$!    checks.
!!$!---------------------------------------------------------------------
!!$      if (Col_diag%num_diag_pts > 0) then
!!$
!!$!---------------------------------------------------------------------
!!$!    check that array dimensions are sufficiently large for the number 
!!$!    of columns requested.
!!$!---------------------------------------------------------------------
!!$        if (Col_diag%num_diag_pts > MAX_PTS) then
!!$          call error_mesg ('donner_deep_mod', 'donner_deep_init: &
!!$         &must reset MAX_PTS or reduce number of diagnostic points', &
!!$                                                           FATAL)  
!!$        endif
!!$
!!$!---------------------------------------------------------------------
!!$!    check that the specified time at which diagnostics are to be 
!!$!    activated has been specified.
!!$!---------------------------------------------------------------------
!!$        do n=1,3
!!$          if (diagnostics_start_time(n) == 0) then
!!$            call error_mesg ('donner_deep_mod', 'donner_deep_init:&
!!$             &year, month and/or day invalidly specified for column '//&
!!$                  'diagnostics starting time', FATAL)
!!$          endif
!!$        end do
!!$
!!$!---------------------------------------------------------------------
!!$!    define a time_type variable indicating the requested time to begin
!!$!    outputting diagnostics.
!!$!---------------------------------------------------------------------
!!$        Time_col_diagnostics = set_date (diagnostics_start_time(1), &
!!$                                         diagnostics_start_time(2), &   
!!$                                         diagnostics_start_time(3), &   
!!$                                         diagnostics_start_time(4), &   
!!$                                         diagnostics_start_time(5), &   
!!$                                         diagnostics_start_time(6) )    
!!$
!!$!---------------------------------------------------------------------
!!$!    allocate space for the arrays used to specify the diagnostics 
!!$!    columns and the output units. initialize the arrays with bogus
!!$!    values.
!!$!---------------------------------------------------------------------
!!$        allocate (col_diag_unit    (Col_diag%num_diag_pts) )
!!$        allocate (col_diag_lon     (Col_diag%num_diag_pts) )
!!$        allocate (col_diag_lat     (Col_diag%num_diag_pts) )
!!$        allocate (col_diag_i       (Col_diag%num_diag_pts) )
!!$        allocate (col_diag_j       (Col_diag%num_diag_pts) )
!!$        col_diag_unit  = -1
!!$        col_diag_lon   = -1.0
!!$        col_diag_lat   = -1.0
!!$        col_diag_i     = -1
!!$        col_diag_j     = -1
!!$
!!$!---------------------------------------------------------------------
!!$!    call initialize_diagnostic_columns to determine the locations 
!!$!    (i,j,lat and lon) of any diagnostic columns in this processor's
!!$!    space and to open output files for the diagnostics.
!!$!---------------------------------------------------------------------
!!$        call initialize_diagnostic_columns   &
!!$                     (mod_name, num_diag_pts_latlon, num_diag_pts_ij, &
!!$                      i_coords_gl, j_coords_gl, lat_coords_gl, &
!!$                      lon_coords_gl, lonb, latb, do_column_diagnostics, &
!!$                      col_diag_lon, col_diag_lat, col_diag_i,  &
!!$                      col_diag_j, col_diag_unit)
!!$
!!$!---------------------------------------------------------------------
!!$!    verify that requested pressure cutoff for column diagnostics output
!!$!    is valid. define the model k index which corresponds (kstart).
!!$!---------------------------------------------------------------------
!!$        do k=1,size(pref(:))
!!$          if (pref(k) >= diagnostics_pressure_cutoff) then
!!$            Col_diag%kstart = k
!!$            exit
!!$          endif
!!$        end do
!!$
!!$!----------------------------------------------------------------------
!!$!    if the specified pressure is larger than any pressure level in the
!!$!    model grid, write an error message.
!!$!----------------------------------------------------------------------
!!$        if (Col_diag%kstart == -99) then
!!$          call error_mesg ( 'donner_deep_mod', 'donner_deep_init: &
!!$           &diagnostics_pressure_cutoff is higher than pressure at '//&
!!$                                     'any model level', FATAL)
!!$        endif
!!$
!!$!----------------------------------------------------------------------
!!$!   if column diagnostics is not requested, define the components of
!!$!   Col_diag that will be needed.
!!$!----------------------------------------------------------------------
!!$      else
!!$        Col_diag%in_diagnostics_window = .false.
!!$        Col_diag%ncols_in_window = 0
!!$      endif
!!$
!!$!----------------------------------------------------------------------
!!$!    allocate space for the array elements of the donner_column_diag_type
!!$!    variable Col_diag. These arrays remain for the life of the job and
!!$!    will be defined for each physics window as it is entered.
!!$!----------------------------------------------------------------------
!!$      allocate (Col_diag%i_dc(Col_diag%num_diag_pts))
!!$      allocate (Col_diag%j_dc(Col_diag%num_diag_pts))
!!$      allocate (Col_diag%unit_dc(Col_diag%num_diag_pts))
!!$      allocate (Col_diag%igl_dc(Col_diag%num_diag_pts))
!!$      allocate (Col_diag%jgl_dc(Col_diag%num_diag_pts))

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    7. FILL THE DONNER_PARAM_TYPE VARIABLE WITH VALUES THAT HAVE BEEN 
!       DEFINED HERE.                  
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!----------------------------------------------------------------------
!    define the components of Param that come from constants_mod. see 
!    donner_types.h for their definitions.
!----------------------------------------------------------------------
      Param%dens_h2o        = DENS_H2O
      Param%rdgas           = RDGAS
      Param%grav            = GRAV
      Param%cp_air          = CP_AIR  
      Param%pie             = PIE
      Param%kappa           = KAPPA
      Param%rvgas           = RVGAS
      Param%seconds_per_day = SECONDS_PER_DAY
      Param%hlv             = HLV
      Param%hlf             = HLF
      Param%hls             = HLS
      Param%kelvin          = KELVIN

!----------------------------------------------------------------------
!    store the parameters defined in this module into the 
!    donner_parameter_type variables Param. these variables are defined
!    above.
!----------------------------------------------------------------------
      Param%cp_vapor                = CP_VAPOR
      Param%parcel_dp               = PARCEL_DP
      Param%upper_limit_for_lfc     = UPPER_LIMIT_FOR_LFC
      Param%pstop                   = PSTOP
      Param%cld_base_vert_vel       = CLD_BASE_VERT_VEL
      Param%dp_of_cloud_model       = DP_OF_CLOUD_MODEL
      Param%cloud_base_radius       = CLOUD_BASE_RADIUS
      Param%wdet                    = WDET
      Param%rbound                  = RBOUND
      Param%wbound                  = WBOUND
      Param%freeze_fraction         = FREEZE_FRACTION
      Param%virt_mass_co            = VIRT_MASS_CO
      Param%pdeep_mc                = PDEEP_MC
      Param%tr_insert_time          = TR_INSERT_TIME
      Param%autoconv_rate           = AUTOCONV_RATE
      Param%autoconv_threshold      = AUTOCONV_THRESHOLD
      Param%tfre                    = TFRE
      Param%dfre                    = DFRE
      Param%evap_in_downdrafts      = EVAP_IN_DOWNDRAFTS
      Param%evap_in_environ         = EVAP_IN_ENVIRON      
      Param%entrained_into_meso     = ENTRAINED_INTO_MESO
      Param%d622                    = D622
      Param%d608                    = D608
      Param%upper_limit_for_lcl     = UPPER_LIMIT_FOR_LCL
      Param%tmin                    = TMIN
      Param%anvil_precip_efficiency = ANVIL_PRECIP_EFFICIENCY
      Param%meso_lifetime           = MESO_LIFETIME
      Param%meso_ref_omega          = MESO_REF_OMEGA
      Param%tprime_meso_updrft      = TPRIME_MESO_UPDRFT
      Param%meso_sep                = MESO_SEP
      Param%ref_press               = REF_PRESS
      Param%meso_down_evap_fraction = MESO_DOWN_EVAP_FRACTION
      Param%meso_up_evap_fraction   = MESO_UP_EVAP_FRACTION
      Param%istart                  = ISTART

      Param%max_entrainment_constant_gate = MAX_ENTRAINMENT_CONSTANT_GATE
      Param%max_entrainment_constant_kep  = MAX_ENTRAINMENT_CONSTANT_KEP
      Param%pdeep_cv                      = PDEEP_CV
      Param%cdeep_cv                      = CDEEP_CV
      Param%kpar                          = KPAR
      Param%r_conv_land                   = R_CONV_LAND
      Param%r_conv_ocean                  = R_CONV_OCEAN 
      Param%n_land                        = N_LAND
      Param%n_ocean                       = N_OCEAN
      Param%delz_land                     = DELZ_LAND
      Param%delz_ocean                    = DELZ_OCEAN
      Param%cell_liquid_eff_diam_def      = CELL_LIQUID_EFF_DIAM_DEF 
      Param%cell_ice_geneff_diam_def      = CELL_ICE_GENEFF_DIAM_DEF
      Param%anvil_levels                  = ANVIL_LEVELS 

      allocate (Param%arat(kpar))
      allocate (Param%ensemble_entrain_factors_gate(kpar))
      allocate (Param%ensemble_entrain_factors_kep(kpar))
      Param%arat                          = ARAT
      Param%ensemble_entrain_factors_gate = ENSEMBLE_ENTRAIN_FACTORS_GATE
      Param%ensemble_entrain_factors_kep  = ENSEMBLE_ENTRAIN_FACTORS_KEP

      allocate (Param%dgeice (ANVIL_LEVELS))
      allocate (Param%relht  (ANVIL_LEVELS))
      Param%dgeice  = DGEICE               
      Param%relht   = RELHT                


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    8. STORE THE NAMELIST VARIABLES THAT NEED TO BE MADE AVAILABLE 
!       OUTSIDE OF THIS MODULE INTO THE DONNER_NML_TYPE VARIABLE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      Nml%parcel_launch_level         = parcel_launch_level
      Nml%allow_mesoscale_circulation = allow_mesoscale_circulation
      Nml%entrainment_constant_source = entrainment_constant_source
      Nml%donner_deep_freq            = donner_deep_freq             
      Nml%model_levels_in_sfcbl       = model_levels_in_sfcbl        
      Nml%cell_liquid_size_type       = cell_liquid_size_type 
      Nml%cell_ice_size_type          = cell_ice_size_type
      Nml%cell_liquid_eff_diam_input  = cell_liquid_eff_diam_input
      Nml%cell_ice_geneff_diam_input  = cell_ice_geneff_diam_input
      Nml%meso_liquid_eff_diam_input  = meso_liquid_eff_diam_input
      Nml%do_average                  = do_average


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    9. END OF SUBROUTINE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!--------------------------------------------------------------------
!    set flag to indicate that donner_deep_mod has been initialized.
!--------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------




  END SUBROUTINE DONNER_INIT


!---------------------------------------
  SUBROUTINE LOOKUP_ES_K(TEMP, ES, STATUS)

    USE MESSY_MAIN_TOOLS,             ONLY: TLUCUA, JPTLUCU2, JPTLUCU1

    REAL     :: TEMP, ES
    INTEGER  :: STATUS
    
    INTEGER  :: IT

    STATUS = 0.
    ES     = 0.
    IT     = INT(TEMP*1000.)
    IF ( (IT < JPTLUCU1) .OR. (IT > JPTLUCU2)) THEN
      STATUS = 1
    ELSE
!      ES     = TLUCUA(IT)*RVGAS/RDGAS
      IF (TEMP > 273.15_dp) THEN
        ES      = 610.78_dp *exp( (temp - 273.15_dp) / &
                ( (temp - 273.15_dp) + 238.3_dp ) * 17.2694_dp )
      ELSE
        ES      = exp( -6140.4_dp / ( 273._dp + (temp - 273.15_dp) ) &
                       + 28.916_dp )
      ENDIF
    ENDIF

    RETURN

  END SUBROUTINE LOOKUP_ES_K

!-----------------------------------------

END MODULE MESSY_CONVECT_DONNER_ADDITIONS
