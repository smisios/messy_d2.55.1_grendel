MODULE MESSY_ATTILA

  ! VERSION: see "modver" below
  ! AUTHOR(S):
  !  Original code for ECHAM4:
  !      CH. REITHMEIER, DLR OBERPFAFFENHOFEN, 18.05.98
  !  Separation from ECHAM4 and transition to F90:
  !      Volker Grewe, Gabriele Erhardt, DLR
  !      Michael Traub, Patrick Joeckel, MPICH, May 2003
  !  Further development:
  !      Michael Traub, Patrick Joeckel, MPICH, Sep/Oct 2003
  !      Andrea Pozzer, MPICH, Sep/Oct 2005
  !      Patrick Joeckel, MPICH, DLR, Jul 2007 - ...
  !      Sabine Brinkop, MIM,DLR, Dec 2007 - ...

  USE messy_main_constants_mem, ONLY: DP, cpd=>cp_air, rd, g, a=>radius_earth
  ! op_pj_20140212+
  ! new, to resolve circular dependency:
  !  - messy_main_tendency_bi uses messy_main_tracer_mem_bi uses ...
  !         ... messy_attila_e5 uses messy_main_tendency_bi
  USE messy_attila_mem,         ONLY: NCELL, NGCELL, AMCELL
  ! op_pj_20140212-

  IMPLICIT NONE
  PRIVATE ! default for all
  SAVE

  ! PUBLIC to be USEd by messy_attila AND
  !                      messy_main_tracer_mem_bi, messy_main_tracer_bi
  PUBLIC :: NCELL, NGCELL, AMCELL, DP ! op_pj_20140212

  ! PUBLIC SUBROUTINES CALLED BY GCM/CTM
  ! ### SMIL: ATTILA_INITIALIZE
  PUBLIC  :: ATTILA_READ_NML_CTRL
  PUBLIC  :: ATTILA_MESSAGE
  PUBLIC  :: ATTILA_GLOBAL_INIT
  !           -> ATTILA_INIT_GETA
  PUBLIC  :: ATTILA_INICOM_1
  PUBLIC  :: ATTILA_INICOM_2
  !           -> ATTILA_INICOMPH
  !
  ! ### SMIL: ATTILA_INIT_MEMORY
  PUBLIC  :: ATTILA_ALLOC
  !
  ! ### SMIL: ATTILA_INIT_COUPLING
  !
  ! ### SMIL: ATTILA_INIT_TRACER
  ! -> SMIL: attila_initialize_positions (.NOT. LTRAJEC)
  PUBLIC  :: ATTILA_INIPOS
  !           -> ATTILA_PRESH
  !           -> ATTILA_THEOCELL
  !           -> ATTILA_CELLPOS
  PUBLIC  :: ATTILA_REALLOC
  ! -> SMIL: attila_update_celldist
  PUBLIC  :: ATTILA_COUNT_CELLS
  !
  ! ### SMIL: ATTILA_GLOBAL_START
  ! -> SMIL: attila_initialize_positions_t (LTRAJEC)
  PUBLIC  :: ATTILA_RESETPOS_TRAJEC
  ! -> SMIL: attila_update_celldist
  !
  PUBLIC  :: ATTILA_FIRSTL           ! either ...
  PUBLIC  :: ATTILA_FIRSTL_EXT       ! ... or
  !
  PUBLIC  :: ATTILA_FILL
  ! -> SMIL: attila_update_celldist
  !
  ! -> SMIL: attila_integrate
  PUBLIC  :: ATTILA_DRIVE
  !           -> ATTILA_CAT_MOVE
  !           -> ATTILA_RKO4
  ! -> SMIL: attila_update_celldist
  !
  ! ### SMIL: ATTILA_GLOBAL_END
  ! -> ATTILA_CONVECT
  PUBLIC  :: ATTILA_CONVGPC
  PUBLIC  :: ATTILA_CONVTRAJ
  PUBLIC  :: ATTILA_SAVE_POS ! SAVE POSITION INFORMATION FOR OUTPUT
  !
  ! ### SMIL: ATTILA_FREE_MEMORY
  PUBLIC  :: ATTILA_GLOBAL_EXIT

  ! PRIVATE SUBROUTINES CALLED INTERNALLY BY ATTILA
  !PUBLIC  :: ATTILA_INICOMPH
  !PUBLIC  :: ATTILA_CELLPOS
  !           -> ATTILA_CALIND
  !PUBLIC  :: ATTILA_CALIND
  !PUBLIC  :: ATTILA_PRESH
  !PUBLIC  :: ATTILA_RKO4
  !           -> ATTILA_LOGGRID
  !           -> ATTILA_VELO
  !           -> ATTILA_REPOS
  !           -> ATTILA_CALINDV
  !PUBLIC  :: ATTILA_VELO
  !           -> ATTILA_GINDEXV_Z
  !           -> ATTILA_3DINTER_V
  !           -> ATTILA_3DINTER_2
  !PUBLIC  :: ATTILA_REPOS
  !PUBLIC  :: ATTILA_CALINDV
  !           -> ATTILA_GINDEXV
  !           -> ATTILA_GINDEXV_Z
  !PUBLIC  :: ATTILA_GINDEXV
  !PUBLIC  :: ATTILA_3DINTER_V
  !PUBLIC  :: ATTILA_3DINTER_2
  !PUBLIC  :: ATTILA_THEOCELL
  !           -> ATTILA_PRESH
  !PUBLIC  :: ATTILA_CONVTRAJ
  !          -> ATTILA_GINDEXV
  ! SPECIAL FOR CAT
  !PUBLIC  :: ATTILA_CAT_MOVE
  !          -> ATTILA_GINDEXV

  ! PRIVATE FUNCTIONS CALLED INTERNALLY BY ATTILA
  !PRIVATE :: GINDEX         ! INTEGER
  !PRIVATE :: GINDEX_V       ! INTEGER ! op_sb_20140515
  !PRIVATE :: NRANWEI        ! INTEGER

  ! INTRINSIC PROCEDURES
  INTRINSIC :: INT, REAL, SIZE, ASSOCIATED        &
       , MAX, ASIN, SIN, SQRT                     &
       , MIN, MERGE, MOD, SUM, NINT, COS, PRESENT &
       , ALLOCATED, LOG                           &
       , ABS, NULL

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'attila'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '4.0a'

  REAL(DP),  PARAMETER, PUBLIC :: kappa     = rd/cpd  ! op_sb_20140210
  REAL(dp),  PARAMETER         :: theta_top = 3300._dp! Kelvin ! op_sb_20140826

  INTEGER, PARAMETER, PUBLIC :: NR1 = 4 ! no. of random vectors for init.
  INTEGER, PARAMETER, PUBLIC :: NR2 = 2 ! no. of spare-random vect. for init.

  ! #####################################################################
  ! CTRL-NAMELIST PARAMETERS
  ! GLOBAL SETTINGS
  INTEGER,  PUBLIC    :: NCHUNK    ! NUMBER OF CHUNKS IN PARALLEL LOOPS
  REAL(dp), PUBLIC    :: CPGBAVE   ! CELLS PER (GRID-BOX) ON AVERAGE
  LOGICAL,  PUBLIC    :: LLTINFO   ! TRUE FOR BEING VERBOSE
  INTEGER,  PUBLIC    :: I_PBLH_METHOD   ! HOW TO CALC. PBLH
  ! PROCESS SPECIFIC SETTINGS
  REAL(dp), PUBLIC   :: ADICO(3)  ! (1) DIFFUSION COEFF
                                  !     (HORIZ.,IN FREE ATM.,IN M^2/S)
                                  ! (2) DIFFUSION COEFF
                                  !     (HORIZ.,IN BOUN LAY.,IN M^2/S)
                                  ! (3) DIFFUSION COEFF (VERTICAL 1/S)

  LOGICAL, PUBLIC    :: LLTBLTURB ! SWITCH FOR BL TURBULENCE
                                  !     (RANDOM REASSIGNMENT)
  LOGICAL, PUBLIC    :: LLCONV    ! TRUE FOR CONVECTION
  LOGICAL, PUBLIC    :: LLCAT     ! TRUE FOR CLEAR AIR TURBULENCE
  LOGICAL, PUBLIC    :: L_LG_FMM    ! TRUE FOR FINIT-MASS METHOD
  LOGICAL, PUBLIC    :: LVDIAG    ! VELOCITY DIAGNOSTICS
  ! SPECIAL MODI
  ! - RESOLUTION INDEPENDENT NUMBER OF CELLS:
  INTEGER, PUBLIC    :: I_NCELL   ! >0 NUMBER OF CELLS INDEPENDENT OF
                                  ! GRIDPOINT RESOLUTION
  ! - TRAJECTORY MODE:
  LOGICAL, PUBLIC    :: LTRAJEC   ! SWITCH FOR TRAJECTORY MODE
  INTEGER, PUBLIC    :: LTRAJEC_DATE(5) ! YEAR,MONTH,DAY,HOUR, MINUTE
                                        ! FOR DEFAULT-INIT
  ! OVERWRITE INDIVIDUAL START DATES ?
  LOGICAL, PUBLIC    :: LTRAJEC_SAME_DATE
  ! op_sb_20140121+
  ! 1:eta,2:theta,3:sigma (vertical coordinate)
  INTEGER, PUBLIC    :: I_VERT = 1
  ! reference pressure for theta vertical coordinate
  ! !!! not for use as pressure reference level;
  !     please use sigma_ref_g * apzero instead  !!!
  REAL(dp), PUBLIC   :: press_ref = -1.0_dp
  ! op_sb_20140121-
  ! #####################################################################

  ! #####################################################################
  ! BASE-MODEL SPECIFIC PARAMETERS, NEED TO BE SET VIA SUBROUTINE
  INTEGER, PUBLIC :: NGL    ! NUMBER OF GAUSSIAN LATIUDES
  INTEGER, PUBLIC :: NLON   ! NUMBER OF LONGITUDES
  INTEGER, PUBLIC :: NLEV   ! NUMBER OF LEVELS
  INTEGER         :: NLEVP1 ! NUMBER OF LEVELS + 1
  INTEGER         :: NLEVM1 ! NUMBER OF LEVELS - 1
  REAL(dp)        :: DTIME  ! TIME STEP LENGTH
  INTEGER         :: NN     ! MAX MERIDIONAL WAVE NUMBER FOR M=0
  INTEGER         :: NPLVP1 ! NUMBER OF PRESSURE LEVELS + 1
  INTEGER         :: NLMSGL ! NLEV - (NUMBER OF SIGMA LEVELS)
  INTEGER         :: NPLVP2 ! NUMBER OF PRESSURE LEVELS + 2
  INTEGER         :: NLMSLP ! NLMSGL + 1
  ! mu = sin(Gaussian latitudes)
  REAL(dp), DIMENSION(:), ALLOCATABLE :: gl_gmu
  ! HYBRID COEFFICIENTS (a1, ..., an, b1, ..., bn)
  REAL(dp), DIMENSION(:), ALLOCATABLE :: VCT
  INTEGER  :: NVCLEV ! NUMBER OF HYBRID COEFF. PAIRS
  REAL(dp) :: APSURF ! FIXED GLOBAL MEAN OF SURFACE PRESSURE
  REAL(dp) :: APZERO ! GLOBAL REFERENCE SURFACE PRESSURE
  ! PARALLEL MODE: BASE MODEL TAKES CARE OF DECOMPOSITION BETWEEN CPUs
  LOGICAL                  :: L_PARALLEL_MODE = .false.
  ! #####################################################################

  ! #####################################################################
  ! GLOBAL ATTILA VARIABLES --------------------------------------------
  ! (COMMON LTPARAM; ARRAY SIZES AND OTHER STUFF)
  INTEGER,  PARAMETER         :: JPNAPOS =   3  ! FIRST DIM OF ARRAY APOS
  INTEGER,  PARAMETER         :: JPNNPOS =   5  ! FIRST DIM OF ARRAY NPOS

  ! (COMMON LTCOM ; LAGRANGIAN SCHEME)
  REAL(dp), PUBLIC     :: AMTOT         ! GLOBAL AIR MASS
  !
  INTEGER  :: NAPOS          ! FIRST DIM OF ARRAY APOS
  INTEGER  :: NNPOS          ! FIRST DIM OF ARRAY NPOS
  INTEGER  :: NHLON          ! NLON/2
  INTEGER  :: NGLP1          ! NGL+1
  INTEGER  :: NGLP2          ! NGL+2
  INTEGER  :: NCHUNKM1       ! NCHUNK-1
  INTEGER  :: NCEVL          ! NCELL/NCHUNK (VECTOR LENGTH OF CELL LOOPS)
  REAL(dp) :: ARAD           ! RAD PER DEG = API/180
  REAL(dp) :: ADEGPM         ! DEG PER METRE = 180/(A*API)
  REAL(dp) :: GDHLON         ! LONGITUDE GRID (GRID POINTS = BOX BOUNDARIES)
  REAL(dp) :: GDLATFT2       ! LATITUDE  GRID (GRID POINTS = BOX BOUNDARIES)
  REAL(dp) :: GDLATLT2       ! GDLON/2
  REAL(dp) :: DTIMED2        ! DTIME/2, HALF  TIME STEP, FOR *RUNGE-*KUTTA
  REAL(dp) :: DTIMED6        ! DTIME/6, SIXTH TIME STEP, FOR *RUNGE-*KUTTA
  LOGICAL, PUBLIC :: LLTMOCA !  TRUE, IF *ADICO* > 0 (I.E. MONTE CARLO DIFF.)

  ! COMMON LTCOM3DINT (INTERPOLATION ROUTINE)
  ! COMMON LTCOMPH (PHYSICS, MONTE-CARLO-DIFFUSION)
  !REAL(dp) :: ADICO(3) ! -> see CTRL-NAMELIST ABOVE
  REAL(dp) :: ADISP(3)  ! (1) RANDOM DISPLACEMENT (HORIZ., IN FREE ATM.)
                        ! (2) RANDOM DISPLACEMENT (HORIZ., IN BOUND. LAY.)
                        ! (3) RANDOM DISPLACEMENT (VERTICAL)
  REAL(dp) :: ABLPLUS   ! ADDITIONAL BL HEIGHT, TO AVOID ACCUMULATION
                        !  OF CELLS IN BL (IN CASE OF BL TURBULENCE.)

  ! MISC
  LOGICAL :: LALLOC     ! LOGICAL, TRUE IF APOS IS NOT POINTED TO LSPAC
                        ! BUT HAS TO BE ALLOCATED (NAPOS+1,NCELL) AND
                        ! THEREFOR DEALLOCATED
  LOGICAL, PUBLIC :: LINIPOS ! FALSE, IF NO INITIALIZATION UP TO NOW
                             ! TRUE , IF INITIAL. HAS ALREADY HAPPENED

  ! POSITIONS OF AIR PARCEL
  !(NAPOS+1) x NCELL
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: APOS => NULL()  ! 4 x NCELL
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: NPOS => NULL()  ! 5 x NCELL
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: NPOSM1 => NULL()! 5 x NCELL
  ! REQUIRED FOR LGFMM
  REAL(dp),  DIMENSION(:,:,:), POINTER, PUBLIC :: APOS_FMM => NULL()

  REAL(dp)                            :: API     ! CONSTANT 3.1415...
  REAL(dp), DIMENSION(:), ALLOCATABLE :: AARBOX  ! AREA OF GRIDBOX
  ! LONGITUDE GRID (0..360-360/NLON)
  REAL(dp), DIMENSION(:), ALLOCATABLE :: GLON
  ! GAUSSIAN LATITUDES (0..180)
  REAL(dp), DIMENSION(:), ALLOCATABLE :: GLAT  ! = N-POLE..S-POLE
  REAL(dp)                            :: GDLON ! LONG. GRID SPACING, =360/NLON
  REAL(dp), DIMENSION(:), ALLOCATABLE :: GDLAT ! LAT. GRID SPACING,
                                               ! = GLAT(I+1)-GLAT(I)

  ! ETA VALUES ON HALF LEVELS
  REAL(dp), DIMENSION(:), POINTER :: GETAH  => NULL()
  ! ETA VALUES ON FULL LEVELS
  REAL(dp), DIMENSION(:), POINTER :: GETAF  => NULL()
  ! GETAH GRID SPACING, = GETAH(I+1)-GETAH(I)
  REAL(dp), DIMENSION(:), POINTER :: GDETAH => NULL()
  ! GETAF GRID SPACING, = GETAF(I+1)-GETAF(I)
  REAL(dp), DIMENSION(:), POINTER :: GDETAF => NULL()
  ! op_sb_20161202+
  ! GETA GRID GLOBAL
  REAL(dp), DIMENSION(:,:,:), POINTER :: GETA_3DF_P   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: GETA_3DH_P   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: GETA_3DF_L_P => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: GETA_3DH_L_P => NULL()

   ! op_sb_20161202-
  ! ... all 4 must be globally available
  REAL(dp), DIMENSION(:,:,:), POINTER :: GETAH_3D  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: GETAF_3D  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: GDETAH_3D => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: GDETAF_3D => NULL()
!!$  ! op_sb_20131209+
!!$  REAL(dp), DIMENSION(:,:,:), POINTER :: TM1_3D  => NULL()
!!$  ! op_sb_20131209+
  REAL(dp), DIMENSION(:), ALLOCATABLE :: GLON2
  REAL(dp), DIMENSION(:), ALLOCATABLE :: GLAT2

!!$  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE         :: PWU  ! NLON x NLEV
!!$  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE         :: PWV  ! NLON x NLEV
!!$  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE         :: PWW  ! NLON x NLEVP1
!!$  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE         :: PPH   ! HL pressure
!!$  REAL(dp), DIMENSION(:,:),   ALLOCATABLE         :: PPS   ! SF pressure
  ! GLOBAL (ON EACH PE!) FIELDS IN GRID-POINT SPACE
  REAL(dp), DIMENSION(:,:,:), POINTER, PUBLIC :: PWU => NULL()  ! U-WIND
  REAL(dp), DIMENSION(:,:,:), POINTER, PUBLIC :: PWV => NULL()  ! V-WIND
  REAL(dp), DIMENSION(:,:,:), POINTER, PUBLIC :: PWW => NULL()  ! W-WIND
  REAL(dp), DIMENSION(:,:,:), POINTER, PUBLIC :: PPH => NULL()  ! HL pressure
  REAL(dp), DIMENSION(:,:),   POINTER, PUBLIC :: surftemp => NULL()

  ! NOTE: surrogate for initialisation phase required (ATTILA_INIT_GETA)
  REAL(dp), DIMENSION(:,:),   POINTER, PUBLIC :: PPS => NULL()    ! SF pressure
  REAL(dp), DIMENSION(:,:,:), POINTER, PUBLIC :: PTPOT => NULL()  ! pot. temp.
  REAL(dp), DIMENSION(:,:,:), POINTER, PUBLIC :: g_TPOT_h => NULL()  ! pot. temp. on half levels
  REAL(dp), DIMENSION(:,:,:), POINTER, PUBLIC :: PTEMP => NULL()  ! temperature
  REAL(dp), DIMENSION(:,:),   POINTER, PUBLIC :: sigma_ref_g => NULL()  ! reference sigma for i_vert=2  op_sb_20160422
  REAL(dp), DIMENSION(:,:),   POINTER, PUBLIC :: i_ref_g => NULL() ! level indices for press_ref op_sb_20160708
   !
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC  :: PMBOX   ! AIR MASS
  REAL(dp), DIMENSION(:,:),   ALLOCATABLE, PUBLIC  :: PMBOXBL ! "  in BL

  ! first layer in free atmosphere
  INTEGER,  DIMENSION(:,:),   ALLOCATABLE, PUBLIC  :: KHPBL

  ! FOLLOWING VARIABLES CAN BE FOUND IN ATTILA_DRIVE
  ! NUMBER OF CELLS PER BOX (NLON x NLEV x NGL) ON THIS CPU
  INTEGER,  PUBLIC, DIMENSION(:,:,:), ALLOCATABLE, TARGET :: NCB
  ! NUMBER OF CELLS PER BOX AT T-1 (NLON x NLEV x NGL) ON THIS CPU
  INTEGER,  PUBLIC, DIMENSION(:,:,:), ALLOCATABLE, TARGET :: NCBM1
  ! NUMBER OF CELLS PER BOX in BL (NLON x NGL) ON THIS CPU
  INTEGER,  PUBLIC, DIMENSION(:,:),   ALLOCATABLE :: NCBL
  ! NUMBER OF CELLS PER BOX ON ALL CPUs (NLON x NLEV x NGL)
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: GNCB => NULL()
  ! NUMBER OF CELLS PER BOX IN BL ON ALL CPUs (NLON x NGL)
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: GNCBL => NULL()
  ! NUMBER OF CELLS TO BE MOVED BY SUBSIDENCE ON ALL PEs
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: GNCB_MOVE => NULL()
  ! #####################################################################

  ! #####################################################################
  ! FOLLOWING VARIABLES ARE NEEDED IN CONVECTION ROUTINES
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::  LZMFU   ! UPWARD MASSFLUX
  REAL(dp), DIMENSION(:,:),   ALLOCATABLE ::  MAXCONV ! TOP LEVEL OF CONVECTION
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::  PENTR   ! PROB. FOR A AP ...
                                                      ! ... BEING ENTRAINED
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::  LPENTR  !           " (BUT LOCAL)
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::  PDETR   ! PROB. FOR A AP ...
                                                      ! ... BEING DETRAINED
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::  PPENDEU ! UPDRAFT: EN-/DETRAINM.
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::  PPENDED ! DOWNDRAFT: EN-/DETRAINM.
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::  PPSUBS  ! CONVECTION: SUBSIDENCE
  REAL(dp), DIMENSION(:,:),   ALLOCATABLE ::  KTYPELT
  INTEGER,  DIMENSION(:,:),   ALLOCATABLE ::  KCTOPLT
  ! #####################################################################

  ! #####################################################################
  ! ADDED FOR ADDITIONAL VELOCITY DIAGNOSTICS
  REAL(dp), DIMENSION(:),   POINTER, PUBLIC :: uvel => NULL() ! [m/s]
  REAL(dp), DIMENSION(:),   POINTER, PUBLIC :: vvel => NULL() ! [m/s]
  REAL(dp), DIMENSION(:),   POINTER, PUBLIC :: wvel => NULL() ! [1/s]
  ! #####################################################################

  ! #####################################################################
  ! ADDED FOR VERTICAL FLUX CALCULATIONS
  REAL(dp), DIMENSION(:),   POINTER, PUBLIC :: APOSM1_PRESS => NULL()
  ! #####################################################################

  REAL(dp), DIMENSION(:), ALLOCATABLE :: A3DRZ3 ! =1./Z3
  REAL(dp), DIMENSION(:), ALLOCATABLE :: A3DDH1 ! CONSTANTS NEEDED ...
  REAL(dp), DIMENSION(:), ALLOCATABLE :: A3DDH2 ! ... FOR CUBIC INTERPOLATION

CONTAINS

! ###########################################################################
! ### PUBLIC SUBROUTINES
! ###########################################################################

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_READ_NML_CTRL(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    ! PURPOSE:
    !   READ ATTILA NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! AUTHOR(S)
    !   Michael Traub, MPICH, July 2003

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! NAMELIST CTRL
    NAMELIST /CTRL/ NCHUNK &
         ,CPGBAVE          &
         ,LLTINFO          &
         ,I_PBLH_METHOD    &
         ,ADICO     &
         ,LLCONV    &
         ,LVDIAG    &
         ,LLTBLTURB &
         ,LLCAT     &
         ,I_NCELL   &
         ,LTRAJEC   &
         ,LTRAJEC_DATE &
         ,LTRAJEC_SAME_DATE &
         ,I_VERT    &       ! op_sb_20140121
         ,press_ref         ! op_sb_20140121

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr = 'attila_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES WITH DEFAULT VALUES
    NCHUNK   = 48
    CPGBAVE  = 2.2_dp         ! CELLS PER GRID-BOX ON AVERAGE
    LLTINFO  = .FALSE.     ! NO VERBOSE MODE
    I_PBLH_METHOD = 0      ! DEFAULT: ATTILA INTERNAL
    ADICO    =  (/ 0.0_dp,  0.0_dp, 0.0_dp /)
                           ! FORMERLY (/  5300./4.,  5300.,  7.E-11 /)
                           ! SEE DISSERTATION CH. REITHMEIER PAGE 22!!
    LLTBLTURB= .FALSE.     ! NO BOUNDARY LAYER TURBULENCE
    LLCONV   = .FALSE.     ! NO CONVECTION
    LVDIAG   = .FALSE.     ! NO VELOCITY DIAGNOSTIC
    LLCAT    = .FALSE.     ! NO CLEAR AIR TURBULENCE
    I_NCELL  = -1          ! >0 LIMIT NUMBER OF CELLS INDEPENDENT OF RES.
    LTRAJEC  = .FALSE.     ! NO TRAJECTORY MODE
    LTRAJEC_DATE  = (/ 1978,  1, 1, 1, 0 /) ! COMMON START DATE (TRAJECT MODE)
    LTRAJEC_SAME_DATE = .FALSE.             ! DEFAULT: INDIVIDUAL DATES (")
! op_sb_20140211+
    I_VERT   = 1           ! eta vertical coordinate with etadot (default)
    press_ref= -1._dp      ! reference pressure for xi hybrid (theta/sigma)
    !                      ! vertical coordinate
! op_sb_20140211-

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    WRITE(*,*) 'GLOBAL SETTINGS'
    WRITE(*,*) '---------------'
    !
    WRITE(*,*) 'CHUNK-SIZE    : ', NCHUNK
    !
    IF (LLTINFO) THEN
       WRITE(*,*) 'VERBOSE       : ON'
    ELSE
       WRITE(*,*) 'VERBOSE       : OFF'
    END IF
    !
    SELECT CASE(I_PBLH_METHOD)
    CASE(0)
       WRITE(*,*) 'PBLH CALCULAT.: ATTILA INTERNAL'
    CASE(1)
       WRITE(*,*) 'PBLH CALCULAT.: EXTERNAL'
    CASE DEFAULT
       WRITE(*,*) '*** ERROR *** UNKNOW METHOD FOR PBLH CALCULATION'
       ! ERROR
       RETURN
    END SELECT
    !
    WRITE(*,*) '--------------'
    WRITE(*,*) 'VERTICAL GRID:'
    WRITE(*,*) '--------------'
    SELECT CASE(I_VERT)
       CASE(1)
          WRITE(*,*) 'ETA (p/p0)'
       CASE(2)
          WRITE(*,*) 'XI (hybrid sigma - theta)'
       CASE(3)
          WRITE(*,*) 'SIGMA (p/ps)'
       CASE DEFAULT
          WRITE(*,*) '*** ERROR *** UNKNOWN VERTICAL GRID'
          RETURN
    END SELECT

    WRITE(*,*) '----------'
    WRITE(*,*) 'PROCESSES:'
    WRITE(*,*) '----------'
    !
    WRITE(*,*) 'MC-DIFFUSION (HF, HPBL, V): ',ADICO
    !
    IF (LLTBLTURB) THEN
       WRITE(*,*) 'BOUNDARY LAYER TURBULENCE : ON'
    ELSE
       WRITE(*,*) 'BOUNDARY LAYER TURBULENCE : OFF'
    END IF
    !
    IF (LLCONV) THEN
       WRITE(*,*) 'CONVECTION                : ON'
    ELSE
       WRITE(*,*) 'CONVECTION                : OFF'
    END IF
    !
    IF (LLCAT) THEN
       WRITE(*,*) 'CLEAR AIR TURBULENCE      : ON'
    ELSE
       WRITE(*,*) 'CLEAR AIR TURBULENCE      : OFF'
    END IF
    !
    WRITE(*,*) '------------'
    WRITE(*,*) 'DIAGNOSTICS:'
    WRITE(*,*) '------------'
    !
    IF (LVDIAG) THEN
       WRITE(*,*) 'VELOCITY DIAGNOSTIC       : ON'
    ELSE
       WRITE(*,*) 'VELOCITY DIAGNOSTIC       : OFF'
    END IF
    !
    WRITE(*,*) '-------------'
    WRITE(*,*) 'SPECIAL MODI:'
    WRITE(*,*) '-------------'
    !
    IF (LTRAJEC) THEN
       WRITE(*,*) 'MODE                      : TRAJECTORY MODE'
       IF (LTRAJEC_SAME_DATE) THEN
          WRITE(*,*) 'START DATE (YYYY,MM,DD,HH,MI): ',LTRAJEC_DATE
       ELSE
          WRITE(*,*) 'START DATE (YYYY,MM,DD,HH,MI): <INDIVIDUAL>'
       END IF
    ELSE
       WRITE(*,*) 'MODE                      : STANDARD MODE'
       IF (I_NCELL <= 0) THEN
          WRITE(*,*) 'NUMBER OF CELLS           : RESOLUTION DEPENDENT'
       ELSE
          WRITE(*,*) 'NUMBER OF CELLS           : ',I_NCELL
       END IF
       !
    END IF

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE ATTILA_READ_NML_CTRL
! ----------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_MESSAGE

    ! PURPOSE:
    !   - DIAGNOSTIC MESSAGE
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTMESSAGE)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_MESSAGE)
    !
    ! INTERFACE:
    !   - LTMESSAGE WAS CALLED FROM LTINICOM
    !   - ATTILA_MESSAGE IS CALLED VIA MESSY-INITIALIZE

    IMPLICIT NONE

    WRITE(*,1) ' '
    WRITE(*,1) '*********************************************'
    WRITE(*,1) '* ORIGINAL:                                 *'
    WRITE(*,1) '* ATTILA - VERSION  1.1N [OCT 2001]         *'
    WRITE(*,1) '* CHRISTIAN REITHMEIER                      *'
    WRITE(*,1) '* DLR OBERPFAFFENHOFEN                      *'
    WRITE(*,1) '*-------------------------------------------*'
    WRITE(*,2) '* CURRENT VERSION: ',modver
    WRITE(*,1) '*********************************************'

1   FORMAT(2X,A)
2   FORMAT(2X,A,A)

  END SUBROUTINE ATTILA_MESSAGE
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_GLOBAL_INIT(ngl_ex, nlon_ex, nlev_ex, dtime_ex &
       , l_pm_ex                                                   &
       , vct_ex, nvclev_ex, apsurf_ex                              &
       , apzero_ex, gl_gmu_ex, nn_ex                               &
       , NPLVP1_ex, NLMSGL_ex, NPLVP2_ex, NLMSLP_ex, status)

    ! PURPOSE:
    !   SET BASE MODEL SPECIFIC PARAMETERS AND ALLOCATE MEMORY
    !   called from attila_initilize
    !
    ! AUTHOR(S):
    !   Patrick Joeckel, MPICH, May 2003

    IMPLICIT NONE

    ! I/O
    INTEGER,            INTENT(IN)  :: ngl_ex    ! number of gaussian latitudes
    INTEGER,            INTENT(IN)  :: nlon_ex   ! number of longitudes
    INTEGER,            INTENT(IN)  :: nlev_ex   ! number of levels
    REAL(dp),           INTENT(IN)  :: dtime_ex  ! length of time step
    LOGICAL,            INTENT(IN)  :: l_pm_ex   ! switch for 'parallel mode'
    REAL(dp), DIMENSION(:), INTENT(IN)  :: vct_ex  ! hybrid coeff. (a, b)
    INTEGER,            INTENT(IN)  :: nvclev_ex ! number of hyb. coeff. pairs
    ! fixed global mean of surface pressure
    REAL(dp),               INTENT(IN)  :: apsurf_ex
    REAL(dp),               INTENT(IN)  :: apzero_ex     ! op_sb_20131208

    ! mu = sin(Gaussian latitudes)
    REAL(dp), DIMENSION(:), INTENT(IN)  :: gl_gmu_ex
    ! nn: max meridional wave number for m=0.
    INTEGER,                INTENT(IN)  :: nn_ex
    ! NPLVP1: number of pressure levels + 1
    INTEGER,                INTENT(IN)  :: NPLVP1_ex
    ! NLMSGL: nlev - (number of sigma levels)
    INTEGER,                INTENT(IN)  :: NLMSGL_ex
    ! NPLVP2: number of pressure levels + 2
    INTEGER,                INTENT(IN)  :: NPLVP2_ex
    ! NLMSLP: NLMSGL + 1
    INTEGER,                INTENT(IN)  :: NLMSLP_ex
    !
    INTEGER,                INTENT(OUT) :: status

    ! LOCAL
    INTEGER :: JPGL, JPGLP1, JPGLP2, JPNLVM1, JPNLVP1, JPNLEV, JPNLON

    ! ERROR STATUS
    status=1    !ERROR

    ! INIT
    NGL    = ngl_ex
    NLON   = nlon_ex
    NLEV   = nlev_ex
    NLEVP1 = NLEV + 1
    NLEVM1 = NLEV - 1
    DTIME  = dtime_ex
    L_PARALLEL_MODE = l_pm_ex

    JPGL     = NGL
    JPNLEV   = NLEV
    JPNLON   = NLON

    JPGLP1   = JPGL+1
    JPGLP2   = JPGL+2

    JPNLVM1  = JPNLEV-1
    JPNLVP1  = JPNLEV+1

    ALLOCATE(vct(SIZE(vct_ex)))
    vct(:) = vct_ex(:)
    nvclev = nvclev_ex
    apsurf = apsurf_ex
    apzero = apzero_ex   ! op_sb_20131208
    ALLOCATE(gl_gmu(SIZE(gl_gmu_ex)))
    gl_gmu(:) = gl_gmu_ex(:)
    nn = nn_ex
    NPLVP1 = NPLVP1_ex
    NLMSGL = NLMSGL_ex
    NPLVP2 = NPLVP2_ex
    NLMSLP = NLMSLP_ex

    ALLOCATE(AARBOX(JPGL))   ! AREA OF GRIDBOX   N --> S
    ALLOCATE(GLON(JPNLON))   ! LONGITUDE GRID (0..360-360/NLON)
    ALLOCATE(GLAT(JPGLP2))   ! GAUSSIAN LATITUDES (0..180) = N-POLE..S-POLE
    ALLOCATE(GDLAT(JPGLP1))  ! LAT. GRID SPACING, =GLAT(I+1)-GLAT(I)
    ALLOCATE(GETAH_3D(NLON,JPNLVP1,NGL))  ! ETA VALUES ON HALF LEVELS
    ALLOCATE(GETAF_3D(NLON,JPNLEV,NGL))   ! ETA VALUES ON FULL LEVELS

!!$    ALLOCATE(GETA_3DH_P(NLON,JPNLVP1,NGL))
!!$    GETA_3DH_P(:,:,:) = 0.0_dp
    ALLOCATE(GETA_3DH_L_P(NLON,JPNLVP1,NGL))
    GETA_3DH_L_P(:,:,:) = 0.0_dp
    ALLOCATE(GETA_3DF_L_P(NLON,JPNLEV,NGL))
    GETA_3DF_L_P(:,:,:) = 0.0_dp
    ALLOCATE(GETA_3DF_P(NLON,JPNLEV,NGL))
    GETA_3DF_P(:,:,:) = 0.0_dp

    ALLOCATE(GDETAH_3D(NLON,JPNLEV,NGL))
    ALLOCATE(GDETAF_3D(NLON,JPNLVM1,NGL))
!!$    ! op_sb_20131208+
!!$    IF (I_VERT == 2) ALLOCATE(TM1_3D(NLON,JPNLEV,NGL))
!!$    ! op_sb_20131208-
    ALLOCATE(GLON2(JPNLON))  ! GETAH GRID SPACING, =GETAH(I+1)-GETAH(I)
    ALLOCATE(GLAT2(JPGLP1))  ! GETAF GRID SPACING, =GETAF(I+1)-GETAF(I)

    ALLOCATE(A3DRZ3(JPNLVM1)) ! =1./Z3
    ALLOCATE(A3DDH1(JPNLVM1)) ! CONSTANTS NEEDED FOR CUBIC INTERPOLATION
    ALLOCATE(A3DDH2(JPNLVM1)) ! CONSTANTS NEEDED FOR CUBIC INTERPOLATION

!!$    ALLOCATE(PWU(NLON,NLEV,NGL))
!!$    ALLOCATE(PWV(NLON,NLEV,NGL))
!!$    ALLOCATE(PWW(NLON,NLEVP1,NGL))
    ALLOCATE(PMBOX(NLON,NLEV,NGL))  ! MASS OF GRID BOX
    PMBOX(:,:,:) = 0.0
    ALLOCATE(PMBOXBL(NLON,NGL))     ! MASS OF GRID BOX BL
    PMBOXBL(:,:) = 0.0
!!$    ALLOCATE(PPH(NLON,NLEVP1,NGL))
    ALLOCATE(KHPBL(NLON,NGL))       ! FIRST LAYER IN FREE ATMOSPH.
    KHPBL(:,:) = NLEV - 2           ! DEFINED INITIAL VALUE

    ALLOCATE(PPS(NLON,NGL))
    ! apzero instead of apsurf is used to guarantee parcels in the surface layer
    ! (highest surface pressure), necessary for initialization!
    PPS(:,:) = apzero

    ALLOCATE(NCB(NLON,NLEV,NGL))
    NCB(:,:,:) = 0
    ALLOCATE(NCBM1(NLON,NLEV,NGL))
    NCBM1(:,:,:) = 0
    ALLOCATE(NCBL(NLON,NGL))
    NCBL(:,:) = 0

    ALLOCATE(GNCB(NLON,NLEV,NGL))
    GNCB(:,:,:) = 0.0
    ALLOCATE(GNCBL(NLON,NGL))
    GNCBL(:,:) = 0.0

    ! op_sb_20160708+
    ALLOCATE(g_tpot_h(NLON,NLEVP1,NGL))
    g_tpot_h(:,:,:) = 0.0_dp
    ! op_sb_20160708-

     IF (LLCONV) THEN
       ALLOCATE(GNCB_MOVE(NLON,NLEV,NGL))
       GNCB_MOVE(:,:,:) = 0.0
       !
       ALLOCATE(PENTR(NLON,NLEV,NGL))
       ALLOCATE(LPENTR(NLON,NLEV,NGL))
       ALLOCATE(LZMFU(NLON,NLEV,NGL))
       ALLOCATE(MAXCONV(NLON,NGL))
       ALLOCATE(PDETR(NLON,NLEV,NGL))
       !
       ALLOCATE(PPENDEU(NLON,NLEV,NGL))
       ALLOCATE(PPENDED(NLON,NLEV,NGL))
       ALLOCATE(PPSUBS(NLON,NLEV,NGL))
       ALLOCATE(KTYPELT(NLON,NGL))
       ALLOCATE(KCTOPLT(NLON,NGL))
    END IF

    CALL ATTILA_INIT_GETA(status, .TRUE.)
    IF (status /= 0) RETURN

    LINIPOS = .FALSE.

    ! ERROR STATUS
    status=0 ! NO ERROR

  END SUBROUTINE ATTILA_GLOBAL_INIT
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_INICOM_1

    ! PURPOSE:
    !   - INITIALISE SOME VARIABLES
    !
    ! INTERFACE:
    !   - LTINICOM WAS CALLED FROM *CONTROL*
    !   - ATTILA_INICOM_1 IS CALLED VIA MESSY-INITIALIZE
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTINICOM)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_INICOM_1)
    !     - SPLIT INTO ATTILA_INICOM_1 AND ATTILA_INICOM_2
    !     - STATUS ADDED TO PARAMETERS
    !     - CALL OF LTEMIPOINT REMOVED  (EMISSIONS)
    !     - CALL OF LTINICOMSM REMOVED
    !     - CALLS OF RSETI, RSETL REMOVED
    !     - AUTOMATIC CALCULATION OF NCELL
    !     - TRAJECT MODE

    IMPLICIT NONE

    ! NCHUNK SHOULD BE MINIMUM 1
    NCHUNK= MAX(1,NCHUNK)

    IF (.NOT. LTRAJEC) THEN
       NGCELL = NINT ( REAL( NGL*NLON*NLEV, DP) * CPGBAVE )
    ENDIF

    NAPOS= JPNAPOS
    NNPOS= JPNNPOS

  END SUBROUTINE ATTILA_INICOM_1
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_INICOM_2(status)

    ! PURPOSE:
    !   - INITIALISE SOME VARIABLES
    !
    ! INTERFACE:
    !   - LTCOM WAS CALLED FROM *CONTROL*
    !   - ATTILA_INICOM_2 IS CALLED VIA MESSY-INITIALIZE
    !
    ! EXTERNALS:
    !   - ATTILA_INICOMPH (LTINICOMPH):  INITIALISE COMMON BLOCK LTCOMPH
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTINICOM)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_INICOM_2)
    !     - SPLIT INTO ATTILA_INICOM_1 AND ATTILA_INICOM_2
    !     - STATUS ADDED TO PARAMETERS
    !     - CALL OF LTEMIPOINT REMOVED  (EMISSIONS)
    !     - CALL OF LTINICOMSM REMOVED
    !     - CALLS OF RSETI, RSETL REMOVED
    !   S. Brinkop, DLR, Oct 2012
    !     - ETA GRID CALCULATINS REMOVED
    !     - CALL TO ATTILA_INICOM3D REMOVED

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: STATUS

    ! LOCAL
    INTEGER :: J
    REAL    :: ZLAT

    ! ERROR STATUS
    status = 1  ! ERROR

    API=2._dp*ASIN(1._dp)

    NHLON= NLON/2
    NGLP1= NGL+1
    NGLP2= NGL+2
    NCHUNKM1= NCHUNK - 1
    NCEVL = INT(NCELL / NCHUNK)
    ARAD = API/180._dp
    ADEGPM= 180._dp/(A*API)
    GDLON= 360._dp/REAL(NLON)
    GDHLON= GDLON/2.0_dp
    DTIMED2= DTIME/2.0_dp
    DTIMED6= DTIME/6.0_dp

    !
    !*    TEST ARRAY SIZE ETC.
    !    ---- ----- ---- ----
    !
    IF((NAPOS /= 3).OR.(NNPOS /= 5)) THEN
       PRINT*,' === ATTILA_INICOM: NAPOS NOT 3 OR NNPOS NOT 5! : ' &
            ,NAPOS,NNPOS
       PRINT*,'        ABORT.'
       RETURN ! STATUS = 1
    ENDIF
    !
    IF(NHLON*2 /= NLON) THEN
       PRINT*,' === ATTILA_INICOM: WARNING: NUMBER OF LONGITUDES (',NLON,')'
       PRINT*,'         SHOULD BE BETTER EVEN! RESULTS MAY BE INEXACT!'
    ENDIF
    !
    IF(NCHUNK*NCEVL /= NCELL) THEN
       PRINT*,' === ATTILA_INICOM: ERROR: NCELL NOT DIVISIBLE BY NCHUNK.'
       PRINT*,'     NCELL,NCHUNK,NCEVL= ',NCELL,NCHUNK,NCEVL,'. ABORT.'
       RETURN ! STATUS = 1
    ENDIF

    !
    !*    INITIALISE LONGITUDE GRIDS.
    !     ---------- --------- ------
    !
    DO J=1,NLON
       GLON (J)=          REAL(J-1, dp)*GDLON
       GLON2(J)= GDHLON + REAL(J-1, dp)*GDLON
    ENDDO
    !
    !*    INITIALISE LATITUDE  GRIDS.
    !     ---------- --------  ------
    !
    DO J=1,(NGL/2)
       ZLAT= ASIN(GL_GMU(J))*180._dp/API
       GLAT(  1  +J)= 90._dp - ZLAT
       GLAT(NGLP2-J)= 90._dp + ZLAT
    ENDDO
    GLAT(1)= 0._dp
    GLAT(NGLP2)= 180._dp

    DO J=1,NGL
       GLAT2(J)= (GLAT(J)+GLAT(J+1)) / 2._dp
    ENDDO
    GLAT2(1)= 0._dp
    GLAT2(NGLP1)= 180._dp

    !
    !*    INITIALISE LATITUDE SPACING GRID.
    !     ---------- -------- ------- -----
    !
    DO J=1,NGLP1
       GDLAT(J)= GLAT(J+1) - GLAT(J)
    ENDDO
    GDLATFT2= GDLAT(1    )*2._dp
    GDLATLT2= GDLAT(NGLP1)*2._dp

    !
    !*       CALCULATE AREA OF GRIDBOXES.
    !        --------- ---- -- ----------
    !
    DO J=1,NGL
       AARBOX(J)= GDLON*ARAD * (A*A) *  &
            (   SIN( ARAD*(90._dp - GLAT2(J  )) ) -  &
                SIN( ARAD*(90._dp - GLAT2(J+1)) )   )
    ENDDO

    !
    !*    PRINT OUT PARAMETERS
    !     ----- --- ----------
    !
    IF(LLTINFO)   THEN
       PRINT*,'NCELL, NAPOS, NNPOS, NHLON, NGLP1, NGLP2:'
       PRINT*, NCELL, NAPOS, NNPOS, NHLON, NGLP1, NGLP2
       PRINT*,'NCHUNK, NCHUNKM1, NCEVL:'
       PRINT*, NCHUNK, NCHUNKM1, NCEVL
       PRINT*,'ARAD, ADEGPM, GDLON, GDHLON:'
       PRINT*, ARAD, ADEGPM, GDLON, GDHLON
       PRINT*,'     - - - - - - - - - - - - - - - - - - - -     '
900    FORMAT(1X,9G10.3)
       PRINT*,'GLON: '
       WRITE(*,900) (GLON(J),J=1,NLON)
       PRINT*,'GLAT: '
       WRITE(*,900) (GLAT(J),J=1,NGLP2)
       PRINT*,'GDLAT: '
       WRITE(*,900) (GDLAT(J),J=1,NGLP1)
       PRINT*,'GDLATFT2, GDLATLT2: ', GDLATFT2, GDLATLT2
       PRINT*,'GLON2: '
       WRITE(*,900) (GLON2(J),J=1,NLON)
       PRINT*,'GLAT2: '
       WRITE(*,900) (GLAT2(J),J=1,NGLP1)
       PRINT*,'AARBOX: '
       WRITE(*,900) (AARBOX(J),J=1,NGL)
       PRINT*,'     - - - - - - - - - - - - - - - - - - - -     '
       PRINT*,'DTIMED2, DTIMED6: ',DTIMED2, DTIMED6
       PRINT*,'     - - - - - - - - - - - - - - - - - - - -     '
    ENDIF  ! END OF LLTINFO BLOCK

    PRINT*,'     - - - - - - - - - - - - - - - - - - - -     '

    PRINT*,'     - - - - - - - - - - - - - - - - - - - -     '
    CALL ATTILA_INICOMPH
    PRINT*,'     - - - - - - - - - - - - - - - - - - - -     '

    PRINT*,'============ END OF ATTILA_INICOM_2 ========================='

    ! ERROR STATUS
    STATUS = 0   ! NO ERROR

  END SUBROUTINE ATTILA_INICOM_2
! -------------------------------------------------------------------

! -------------------------------------------------------------------
   SUBROUTINE ATTILA_ALLOC(status, ZNCELL,LSPAC,LSPACINT,LSPACINTM1 &
        , LSPACPRESSM1)

     ! PURPOSE:
     !   - IN CASE OF RUNNING THIS MODULE WITHOUT ECHAM5 APOS, NPOS,
     !     NPOSM1 IS ALLOCTED IN THIS SR AND A SWITCH IS SET FOR
     !     DEALLOCATIAN.
     !   - IN CASE OF RUNNING WITH ECHAM5 THE POSITION ARE WRITTEN
     !     INTO CHANNELS; THIS IS NECESSARY FOR RESTART
     !
     ! AUTHOR(S):
     !   P. Joeckel, MPI Mainz
     !
     ! INTERFACE:
     !   - ATTILA_ALLOC IS CALLED VIA MESSY-...
     !
     ! PARAMETERS:
     !   O STATUS
     !   I ZNZELL
     !   I LSPAC
     !   I LSPACINT
     !   I LSPACINTM1
     !   I LSPACPRESSM1

     IMPLICIT NONE

     ! I/O
     INTEGER                :: status
     INTEGER,      OPTIONAL :: ZNCELL
     ! OUTPUT CHANNEL FOR APOS
     REAL(dp),POINTER, OPTIONAL ::  LSPAC(:,:,:,:,:)
     ! OUTPUT CHANNEL FOR NPOS
     REAL(dp),POINTER, OPTIONAL ::  LSPACINT(:,:,:,:,:)
     ! OUTPUT CHANNEL FOR NPOSM1
     REAL(dp),POINTER, OPTIONAL ::  LSPACINTM1(:,:,:,:,:)
     ! OUTPUT CHANNEL FOR APOS (PRESSURE) TM1
     REAL(dp),POINTER, OPTIONAL ::  LSPACPRESSM1(:,:,:,:,:)

     status = 1 ! ERROR

     IF ( PRESENT(ZNCELL)     .EQV. &
          (PRESENT(LSPAC)     .AND. &
          PRESENT(LSPACINT)   .AND. &
          PRESENT(LSPACINTM1) .AND. &
          PRESENT(LSPACPRESSM1)) ) THEN
        ! ERROR
        RETURN
     ENDIF

     IF (PRESENT(ZNCELL)) THEN

        ALLOCATE(APOS(NAPOS+1,NCELL))
        ALLOCATE(NPOS(NNPOS,NCELL))
        ALLOCATE(NPOSM1(NNPOS,NCELL))
        ALLOCATE(APOSM1_PRESS(NCELL))

        LALLOC=.TRUE.

     ELSE

        APOS =>  LSPAC(:,:,1,1,1)

        NPOS =>  LSPACINT(:,:,1,1,1)

        NPOSM1 => LSPACINTM1(:,:,1,1,1)

        APOSM1_PRESS => LSPACPRESSM1(1,:,1,1,1)

        LALLOC=.FALSE.

     END IF

     status = 0 ! OK

   END SUBROUTINE ATTILA_ALLOC
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_INIPOS(HARVESTINIPOS1, HARVESTINIPOS2, STATUS )

    ! called in attila_initialize_positions (called from attila_init_tracer)

    ! PURPOSE:
    !   - INITIALISE CELL POSITIONS (HORIZONTAL POSITION AND ETA COORDINATE)
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTINIPOS)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_INIPOS)
    !
    ! INTERFACE:
    !   - LTINIPOS WS CALLED FROM LTDRIVE
    !   - ATTILA_INIPOS IS CALLED FROM attila_initialize_positions
    !   - attila_initialize_positions IS CALLED FROM attila_init_tracer
    !
    ! EXTERNALS:
    !   - NRANWEI         GET A WEIGHTED RANDOM NUMBER
    !   - ATTILA_PRESH    CALCULATE PRESSURE ON HALF LEVELS
    !   - ATTILA_CELLPOS  POSITION A CELL
    !   - ATTILA_THEOCELL PRINT OUT THEORETICAL NUMBER OF CELLS
    !
    ! REAMRKS:
    !             PPOS(1-3) : HORIZ. &  ETA POSITION
    !             KPOS(1-5) : INTEGER POSITION
    ! METHODS:
    !         ALL CELLS ARE SUPPOSED TO HAVE EQUAL MASS, THEREFORE
    !         NUMBER OF CELLS PER GRID BOX IS PROPORTIONAL TO
    !         PRESSURE DIFFERENCE AND INTEGRAL OF COS(LAT).

    IMPLICIT NONE

    ! I/O
    ! RANDOM NUMBER ARRAY  (FOR POSITIONING)
    REAL(dp), DIMENSION(:,:), INTENT(IN)    :: HARVESTINIPOS1
    ! RANDOM NUMBER ARRAY  (FOR ADJUSTMENT)
    REAL(dp), DIMENSION(:,:), INTENT(IN)    :: HARVESTINIPOS2
    ! STATUS
    INTEGER,                  INTENT(OUT)   :: STATUS

    ! LOCAL
    REAL(dp) :: ZP(NLEVP1,NLON,NGL)  ! 3D PRESSURE FIELD
    REAL(dp) :: ZFRAC, ZPP(NLEVP1)   ! FRACTION
    REAL(dp) :: ZRAD
    REAL(dp) :: ZW(NLON)
    REAL(dp) :: ZWSUM
    INTEGER  :: IFRAC(NGL,NLEV), ILEV(NLEV), ILAT(NGL)
    INTEGER  :: ICELL, ICELL1, IFR
    INTEGER  :: JCELL, JFRAC, J, JL, JG, JK, JC, JC0, JC1
    INTEGER  :: JCL, JCELLSTART, JCELLEND
    INTEGER  :: JCR1_3, JCR1_1, JCR2, DCELL

    ! STATUS ERROR
    STATUS = 1   !ERROR

    ! SET SOME PARAMETERS
    ZRAD= API/180._dp
    JCR1_3 = 0 ! # of random number triples used from HARVESTINIPOS1(:,1:3)
    JCR1_1 = 0 ! # of random numbers        used from HARVESTINIPOS1(:,4)
    JCR2   = 0 ! # of random number pairs   used from HARVESTINIPOS2(:,1:2)

    PRINT*,'================ ATTILA_INIPOS =================='

    !
    !*  1. CALCULATE 3-D PRESSURE FIELD.
    !   -- --------- --- -------- ------
    !
    DO JG=1,NGL
       DO JL=1,NLON
          CALL ATTILA_PRESH(PPS(JL,JG),ZPP)
          ZP(:,JL,JG)=ZPP(:)
       END DO
    END DO

    !
    !*  2. CALCULATE NUMBER OF CELLS IN EACH BAND.
    !   -- --------- ------ -- ----- -- ---- -----
    !
    ICELL= 0
    DO J=1,NLEV
       ILEV(J)= 0
    END DO
    DO J=1,NGL
       ILAT(J)=0
    END DO

    !
    !    2.1. CALCULATE FRACTION
    !
    DO JG=1,NGL
       DO JK=1,NLEV
          ZFRAC=0.0_dp
          DO JL=1,NLON
             ZFRAC  = ZFRAC +                                      &
                  ( SIN( (GLAT2(JG+1)-90._dp)*ZRAD ) -             &
                  SIN( (GLAT2(JG  )-90._dp)*ZRAD )   ) * 0.5_dp *  &
                  ( ZP(JK+1,JL,JG  ) - ZP(JK,JL,JG) )      /       &
                  ( ZP(NLEVP1,JL,JG) - ZP(1 ,JL,JG) )      /       &
                  REAL(NLON, DP)
          END DO
          IFRAC(JG,JK)= NINT( REAL( NINT(REAL(NGCELL,dp)*ZFRAC), dp )/ &
               (REAL(NGCELL,dp)/REAL(NCELL,dp)) )
          ILEV(JK)= ILEV(JK) + IFRAC(JG,JK)
          ILAT(JG)= ILAT(JG) + IFRAC(JG,JK)
          ICELL   = ICELL    + IFRAC(JG,JK)
       END DO
    END DO

    !
    !    2.2. ADJUST NUMBER OF CELLS
    !
    ICELL1=ICELL
    IF(ICELL < NCELL) THEN
       PRINT*,' ICELL < NCELL : ', ICELL, NCELL
       DCELL = NCELL-ICELL
       DO JCELL= ICELL+1,NCELL
          JG=NINT(HARVESTINIPOS2(JCELL,1)*NGL)
          JK=NINT(HARVESTINIPOS2(JCELL,2)*NLEV)
          JCR2 = JCR2 + 1            ! 1 random number pair used
          IF(JG == 0) JG=NGL
          IF(JK == 0) JK=NLEV
          IFRAC(JG,JK)= IFRAC(JG,JK) + 1
          ILEV(JK)= ILEV(JK) + 1
          ILAT(JG)= ILAT(JG) + 1
          ICELL1= ICELL1 + 1
       END DO
       !
    ELSEIF(ICELL > NCELL) THEN
       PRINT*,' ICELL > NCELL : ', ICELL, NCELL
       DCELL = ICELL-NCELL
       JC0= (NCELL/DCELL)
       JC= JC0
       DO JCELL= NCELL+1,ICELL
          !
          IFR= 0
          !
          DO WHILE( JC > 0 .AND. IFR <= 0 )
             ! WE NEED TO ADJUST DCELL=ICELL-NCELL (ICELL > NCELL !!!)
             ! POSITIONS WITH NCELL AVAILABLE RANDOM NUMBER PAIRS (INDEPENDENT
             ! OF THE PARALLEL DECOMPOSITION!). WE HOPE THAT
             ! DCELL << NCELL, THUS WE SEARCH FOR EVERY JCELL
             ! BETWEEN NCELL+1 AND ICELL IN THE RANDOM NUMBER PAIR BLOCK
             ! STARTING AT [0 ... (DCELL-1)]*JC0+1 AND OF LENGTH
             ! JC0=(NCELL/DCELL) FOR A RANDOM NUMBER PAIR THAT
             ! RESULTS IN A PAIR JG,JK WITH IFRAC(JG,JK) > 0 ...
             ! IF THIS DOES NOT WORK, WE HAVE TO STOP!
             JC1 = (JCELL-NCELL-1)*JC0 + 1 + (JC0-JC)
             IF (JC1 > SIZE(HARVESTINIPOS2,1)) THEN
                PRINT*, ' FATAL ERROR: OUT OF BOUNDS'
                PRINT*, '              PROGRAM ABORTED.'
                RETURN
             ENDIF
             JG=NINT(HARVESTINIPOS2(JC1,1)*NGL)
             JK=NINT(HARVESTINIPOS2(JC1,2)*NLEV)
             JCR2 = JCR2 + 1         ! 1 random number pair used
             IF(JG == 0) JG=NGL
             IF(JK == 0) JK=NLEV
             IFR= IFRAC(JG,JK)
             JC= JC - 1
          END DO
          !
          IF (IFR <= 0) THEN
             PRINT*,' FATAL ERROR: ONLY EMPTY GRID BOXES FOUND.'
             PRINT*,'              PROGRAM ABORTED.'
             RETURN
          END IF
          IFRAC(JG,JK)= IFRAC(JG,JK) - 1
          ILEV(JK)= ILEV(JK) - 1
          ILAT(JG)= ILAT(JG) - 1
          ICELL1= ICELL1 - 1
       END DO
    ELSE
       PRINT*,' ICELL = NCELL : ', ICELL, NCELL
       DCELL = 0
    END IF
    !
    IF (ICELL1.NE.NCELL) THEN
       PRINT*,' FATAL ERROR: ICELL1 != NCELL : ',ICELL1, NCELL
       PRINT*,'       PROGRAM ABORTED.'
       RETURN
    END IF
    PRINT*,' NUMBER OF CELLS IN EACH LEVEL :'
    WRITE(*,'(1X,10I7)') (ILEV(JK),JK=1,NLEV)
    PRINT*,' NUMBER OF CELLS IN EACH LATITUDE :'
    WRITE(*,'(1X,10I7)') (ILAT(JG),JG=1,NGL)
    PRINT*,'TOTAL: ',ICELL1,' CELLS.'
    !
    !       PRINT OUT THEORETICAL NUMBER OF CELLS.
    !
    CALL ATTILA_THEOCELL
    !
    !*   3. POSITION CELLS
    !    -- -------- -----
    !
    JCR1_3=1   ! position of first random number triple
    JCR1_1=1   ! position of first random number
    ICELL=1    ! CURRENT CELL NO. TO ALLOCATE
    DO JG=1,NGL
       DO JK=1,NLEV
          IF(IFRAC(JG,JK) == 0) THEN
             PRINT*,' THIS BAND IS EMPTY: LAT,LEV= ',JG,JK
          ELSE
             !
             !     3.1. CALCULATE WEIGHTS AND INIT. COUNTER.
             !
             ZWSUM = 0.0_dp
             DO JL=1,NLON
                ZW(JL)= ZP(JK+1,JL,JG) - ZP(JK,JL,JG)
                ZWSUM = ZWSUM + ZW(JL)
             END DO
             IF(ZWSUM <= 0._dp) THEN
                PRINT*,' M1: S.TH. WENT DEFINITELY WRONG HERE! ABORTED.'
                RETURN
             END IF
             ZWSUM= 1._dp/ZWSUM
             DO JL=1,NLON
                ZW(JL)= ZW(JL) * ZWSUM
             END DO
             !
             !     3.2. DISTRIBUTE CELLS
             !
             JCELLSTART= ICELL  ! FIRST CELL IN BAND
             JCL= 0             ! NO. OF CELLS IN BAND ALREADY POSITIONED
             !
             !     3.2.1. EVENLY DISTRIBUTED CELLS
             !
             DO JL=1,NLON
                DO JFRAC=1,INT(REAL(IFRAC(JG,JK),dp)*ZW(JL))

                   IF (ICELL > NCELL) THEN
                      PRINT*,' M2: S.TH. WENT DEFINITELY WRONG HERE! ABORTED.'
                      RETURN
                   END IF

                   CALL ATTILA_CELLPOS(JCR1_3, HARVESTINIPOS1,  &
                        JL,JG,JK,ICELL)

                   ICELL= ICELL +1
                   JCL= JCL+1
                END DO
             END DO
             IF (JCL > IFRAC(JG,JK)) THEN
                PRINT*,' ERROR: TOO MANY CELLS IN BAND ',JG,',',JK, &
                     ' : ',JCL,IFRAC(JG,JK),' , ABORTED.'
                RETURN
             END IF

             !
             !     3.2.2. RECALCULATE WEIGHTS.
             !
             ZWSUM = 0.0
             DO JL=1,NLON
                ZW(JL)= ZW(JL)*REAL(IFRAC(JG,JK),dp) - INT(ZW(JL)*IFRAC(JG,JK))
                ZWSUM=  ZWSUM + ZW(JL)
             END DO
             IF (ZWSUM > 0.0_dp) THEN
                ZWSUM= 1._dp/ZWSUM
                DO JL=1,NLON
                   ZW(JL)= ZW(JL)*ZWSUM
                END DO
             ELSE
                DO JL=1,NLON
                   ZW(JL)= 1.0_dp/REAL(NLON,dp)
                END DO
             END IF

             !
             !     3.2.3. RANDOMLY DISTRIBUTED CELLS
             !
             DO JFRAC= JCL+1, IFRAC(JG,JK)
                IF (ICELL > NCELL) THEN
                   PRINT*,' M3: S.TH. WENT DEFINITELY WRONG HERE! ABORTED.'
                   RETURN
                END IF

                JL= NRANWEI(JCR1_1,HARVESTINIPOS1,ZW)
                JCR1_1 = JCR1_1 + 1 ! 1 random number used in NRANWEI
                CALL ATTILA_CELLPOS(JCR1_3, HARVESTINIPOS1, &
                     JL,JG,JK,ICELL)
                ICELL= ICELL +1
             END DO

             JCELLEND= ICELL-1 ! LAST CELL IN BAND

             IF(LLTINFO) THEN
                WRITE(*,'(1X,A,2I3,A,I4,A,I7,A,I7,A)')   &
                     'FINISHED LAT/LEV BAND ',JG,JK,     &
                     ' (',IFRAC(JG,JK),' CELLS, NO. ',   &
                     JCELLSTART,'-',JCELLEND, '). '
             END IF
          END IF
       END DO                                 !LEVEL LOOP
    END DO                                    !LOOP OVER NGL
    PRINT*,' *** ',ICELL-1,' OUT OF ',NCELL,' CELLS POSITIONED.'
    IF(ICELL-1 == NCELL) THEN
       PRINT*,' *** THIS IS OK, FINISHED SUCCESSFULLY.'
    ELSE
       PRINT*,' *** THIS IS NOT OK!'
    ENDIF

    PRINT*,' NUMBER OF USED RANDOM NUMBER TRIPLES (1): ', &
         JCR1_3-1, NCELL, SIZE(HARVESTINIPOS1,1)

    PRINT*,' NUMBER OF USED RANDOM NUMBERS        (1): ', &
         JCR1_1-1, NCELL, SIZE(HARVESTINIPOS1,1)

    PRINT*,' NUMBER OF USED RANDOM NUMBER PAIRS   (2): ', &
         JCR2, DCELL, SIZE(HARVESTINIPOS2,1)

    ! SET FLAG TO TRUE: POSITIONS ARE NOW INITIALIZED
    LINIPOS = .TRUE.

    ! INITIALIZE GLOBAL FIELDS
    NPOSM1(:,:) = NPOS(:,:)

    !STATUS ERROR
    STATUS=0

  END SUBROUTINE ATTILA_INIPOS
! -------------------------------------------------------------------

! -------------------------------------------------------------------
   SUBROUTINE ATTILA_REALLOC(flag, TMP_APOS, TMP_NPOS, TMP_NPOSM1 &
        , LSPAC, LSPACINT, LSPACINTM1)

     IMPLICIT NONE

     ! I/O
     INTEGER, INTENT(IN) :: flag
     REAL(dp), DIMENSION(:,:), POINTER :: TMP_APOS
     REAL(dp), DIMENSION(:,:), POINTER :: TMP_NPOS
     REAL(dp), DIMENSION(:,:), POINTER :: TMP_NPOSM1
     !
     ! OUTPUT CHANNEL FOR APOS
     REAL(dp),POINTER, OPTIONAL ::  LSPAC(:,:,:,:,:)
     ! OUTPUT CHANNEL FOR NPOS
     REAL(dp),POINTER, OPTIONAL ::  LSPACINT(:,:,:,:,:)
     ! OUTPUT CHANNEL FOR NPOSM1
     REAL(dp),POINTER, OPTIONAL ::  LSPACINTM1(:,:,:,:,:)

     ! LOCAL

     SELECT CASE(flag)
        !
     CASE(1)
        !
        ! ALLOCATE TEMPORARY (GLOBAL) APOS, NPOS, etc. ...
        ALLOCATE(TMP_APOS(NAPOS+1,NCELL))
        ALLOCATE(TMP_NPOS(NNPOS,NCELL))
        ALLOCATE(TMP_NPOSM1(NNPOS,NCELL))
        !
        IF (LALLOC) THEN
           ! DEALLOCATE PRE-ALLOCATED LOCAL FIELDS
           DEALLOCATE(APOS);   NULLIFY(APOS)
           DEALLOCATE(NPOS);   NULLIFY(NPOS)
           DEALLOCATE(NPOSM1); NULLIFY(NPOSM1)
        ENDIF
        !
        ! ... AND USE FOR INITIALISATION
        APOS   => TMP_APOS
        NPOS   => TMP_NPOS
        NPOSM1 => TMP_NPOSM1
        !
     CASE(2)
        !
        NULLIFY(APOS)
        NULLIFY(NPOS)
        NULLIFY(NPOSM1)
        !
        IF (LALLOC) THEN
           ALLOCATE(APOS(NAPOS+1,NCELL))
           ALLOCATE(NPOS(NNPOS,NCELL))
           ALLOCATE(NPOSM1(NNPOS,NCELL))
        ELSE
           IF (PRESENT(LSPAC)) &
                APOS =>  LSPAC(:,:,1,1,1)
           IF (PRESENT(LSPACINT)) &
                NPOS =>  LSPACINT(:,:,1,1,1)
           IF (PRESENT(LSPACINTM1)) &
                NPOSM1 => LSPACINTM1(:,:,1,1,1)
        ENDIF
        !
     CASE(3)
        !
        DEALLOCATE(TMP_APOS);   NULLIFY(TMP_APOS)
        DEALLOCATE(TMP_NPOS);   NULLIFY(TMP_NPOS)
        DEALLOCATE(TMP_NPOSM1); NULLIFY(TMP_NPOSM1)
        !
     END SELECT

   END SUBROUTINE ATTILA_REALLOC
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_COUNT_CELLS(RNCB, RNCBL, L_INIT)

    ! - called by attila_update_celldist (called twice by attila_global_end
    !   (before (if lstart) and after attila_integrate))
    ! - called by attila_convect

    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(OUT) :: RNCB(NLON,NLEV,NGL)
    REAL(dp), INTENT(OUT) :: RNCBL(NLON,NGL)
    LOGICAL,  INTENT(IN), OPTIONAL :: L_INIT

    ! LOCAL
    ! op_sb20160509+
    REAL(DP)              :: ptemp_lin(nlon,nlev,ngl)
    REAL(DP)              :: ptpot_lin(nlon,nlev,ngl)
    REAL(DP)              :: ptemp_h(nlon,nlevp1,ngl)
    REAL(DP)              :: ppf(nlon,nlev,ngl)
    REAL(DP)              :: deltap(nlon,nlev,ngl)
    REAL(DP)              :: delta1(nlon,nlev,ngl)
    REAL(DP)              :: delta2(nlon,nlev,ngl)
    REAL(dp)              :: del,del1,del2,t_crit
    REAL(DP)              :: xi_crit(nlon,ngl)

    ! op_sb20160509-
    INTEGER               :: JCHUNK, JCELL
    INTEGER               :: JCELL1, JCELL2, JCELL3
    INTEGER               :: JL, JG, JK, JJ
    LOGICAL               :: ZLSTART ! op_sb_20140516
    LOGICAL               :: l_ptemp
    REAL(dp)              :: pit2

    pit2 = api / 2._dp

!CDIR BEGIN NOVECTOR

    ! INITIALIZE
    l_ptemp = ASSOCIATED(ptemp)
    ! -> AT EACH TIMESTEP NCB HAS TO BE RESET
    NCB(:,:,:)   = 0
    RNCB(:,:,:)  = 0.0_dp
    NCBL(:,:)    = 0
    RNCBL(:,:)   = 0.0_dp
   ! op_sb_20160503+
    ptemp_h(:,:,:)   = 0.0_dp
    ptemp_lin(:,:,:) = 0.0_dp
    ptpot_lin(:,:,:) = 0._dp        ! op_sb_20160708

    IF (l_ptemp) THEN

    ! calculate temperatures at half levels and pressure at full levels
    do jg=1,ngl
      do jl=1,nlon
          ptemp_h(jl,1,jg)= ptemp(jl,1,jg)
          ptemp_h(jl,nlevp1,jg)= surftemp(jl,jg)
          ppf(jl,1,jg) = (pph(jl,1,jg)+pph(jl,2,jg))*.5_dp
      enddo
    enddo

    do jg=1,ngl
      do jk=2,nlev
        do jl=1,nlon
          ppf(jl,jk,jg)     = (pph(jl,jk,jg) + pph(jl,jk+1,jg))*.5_dp
          deltap(jl,jk,jg)  = (ppf(jl,jk,jg) - ppf(jl,jk-1,jg)) ! between full levels
          ptemp_h(jl,jk,jg) = (ptemp(jl,jk,jg) * (pph(jl,jk,jg)-ppf(jl,jk-1,jg))+ &
                              ptemp(jl,jk-1,jg) * (ppf(jl,jk,jg)-pph(jl,jk,jg)))/ &
                              deltap(jl,jk,jg)
        enddo
      enddo
    enddo

    do jg=1,ngl
        do jl=1,nlon
          jk=int(i_ref_g(jl,jg))
          del = abs(PPH(jl,JK+1,jg)-PPH(jl,jk,jg))  ! between half levels
          del1 = abs((PPH(jl,jk+1,jg)-sigma_ref_g(jl,jg)*apzero))
          del2 = abs(sigma_ref_g(jl,jg)*apzero - PPH(jl,jk,jg))
          t_crit=(ptemp_h(jl,jk+1,jg) * del2 + ptemp_h(jl,jk,jg) * del1 ) / del
          xi_crit(jl,jg) = t_crit*(100000._dp/(sigma_ref_g(jl,jg)*apzero))**kappa
         enddo
    enddo

    endif
    ! op_sb_20160503-
    ! op_sb_20140516+
    IF (PRESENT(L_INIT)) THEN
       ZLSTART = L_INIT
    ELSE
       ZLSTART = .FALSE.
    ENDIF

    DO JCHUNK=0,NCHUNKM1

       JCELL1= JCHUNK*NCEVL
       JCELL2= JCELL1 + 1
       JCELL3= JCELL1 + NCEVL

       !
       !   THE FOLLOWING LOOP CONTAINS NON-VECTORIZABLE STATEMENTS
       !
       DO JCELL= JCELL2,JCELL3

          JL= INT(NPOS(4,JCELL))
          JG= INT(NPOS(5,JCELL))
          JK= INT(NPOS(3,JCELL))

          ! GET THE PRESSURE VALUE FOR EACH CELL
          ! PRESSURE IN Pa

          GETAH => GETAH_3D(JL,:,JG)

          SELECT CASE (I_VERT)
          CASE(1)
             APOS(4,JCELL) = APOS(3,JCELL) * APZERO
          CASE(2)
             IF (.NOT. LTRAJEC) THEN ! op_sb_20141113
                IF (ZLSTART) THEN
                   APOS(4,JCELL) = APOS(3, JCELL) * APZERO
                   IF (l_ptemp) THEN
                      ! pressure of parcel can be higher than surface pressure
                      ! due to initialisation, correction applied
                      if (apos(4,jcell) .gt. pps(jl,jg)) &
                           apos(4,jcell) = pps(jl,jg) - 1.E-4_dp

                      ! correction of jk necessary due to switch from
                      ! eta to xi vertical grid
                      if (apos(4,jcell) .gt. PPH(jl,JK+1,jg) .or. &
                           apos(4,jcell) .lt.  PPH(jl,jk,jg)) then
                         if (apos(4,jcell) .gt. PPH(jl,JK+1,jg)) then
                            jj = jk+1
                            do while (apos(4,jcell) .gt. pph(jl,jj+1,jg))
                               jj=jj+1
                            enddo
                         elseif( apos(4,jcell) .lt.  PPH(jl,jk,jg)) then
                            jj = jk
                            do while (apos(4,jcell) .lt. pph(jl,jj,jg))
                               jj=jj-1
                            enddo
                         endif
                         jk=jj
                      endif

                      ! between half levels
                      deltap (jl,jk,jg) = (PPH(jl,JK+1,jg)-PPH(jl,jk,jg) )
                      delta1(jl,jk,jg)  = (PPH(jl,jk+1,jg)-apos(4,jcell))
                      delta2(jl,jk,jg)  = (apos(4,jcell)-PPH(jl,jk,jg))

                      ! temperature linearized on parcel pressure
                      ptemp_lin(jl,jk,jg) = (ptemp_h(jl,jk+1,jg) &
                           * delta2(jl,jk,jg)+              &
                           ptemp_h(jl,jk,jg) * delta1(jl,jk,jg)) &
                           /deltap(jl,jk,jg)
                      ptpot_lin(jl,jk,jg) = ptemp_lin(jl,jk,jg) &
                           * (100000._dp/apos(4,jcell))**kappa

                      !op_sb_20160708+
                      if (apos(4,jcell) .gt. (sigma_ref_g(jl,jg)*apzero)    &
                           .and. jk .ge. int(i_ref_g(jl,jg)) ) then
                         APOS(3,JCELL) = ptpot_lin(jl,jk,jg) *              &
                              sin(pit2*(1._dp-apos(4,jcell)/pps(jl,jg))     &
                              /(1._dp-sigma_ref_g(jl,jg)))
                      elseif(apos(4,jcell) .le. (sigma_ref_g(jl,jg)*apzero)&
                           .and. jk .le. int(i_ref_g(jl,jg)) ) then
                         APOS(3,JCELL) = ptpot_lin(jl,jk,jg)
                      endif
                      !op_sb_20160708-

                      APOS(3,jcell)= MIN(apos(3,jcell),theta_top-1.)
                      NPOS(3,JCELL)= REAL(GINDEX_V(APOS(3,JCELL), GETAH), dp)

                      ! control check
                      if (npos(3,jcell) < 1._dp .or. &
                           npos(3,jcell) .gt. REAL(nlev) ) then
                         print*,'npos out of range!', &
                              npos(3,jcell),jcell,apos(3,jcell),apos(4,jcell), &
                              pps(jl,jg),ptemp_lin(jl,jk,jg),jl,jk,jg,         &
                              deltap(jl,jk,jg),delta1(jl,jk,jg),               &
                              ptemp_h(jl,jk,jg),ptemp_h(jl,jk+1,jg)
                      END IF
                   ENDIF
                ELSE ! not zlstart
                   ! calculate new temperature on parcel,
                   ! apos(3,jcell) has changed
                   ! op_sb_20160708+

                   del  = abs(getah_3d(jl,jk+1,jg) - getah_3d(jl,jk,jg))
                   del1 = abs(getah_3d(jl,jk+1,jg) - apos(3,jcell))
                   del2 = abs(apos(3,jcell) - getah_3d(jl,jk,jg))

                   if (apos(3,jcell) .lt. xi_crit(jl,jg) &
                        .and. jk .ge. int(i_ref_g(jl,jg))) then
                      ! linearization of T*f in xi-vertical grid
                      ptemp_lin(jl,jk,jg) = (ptemp_h(jl,jk+1,jg)*        &
                           sin(pit2*(1._dp-pph(jl,jk+1,jg)/pps(jl,jg))/  &
                           (1._dp-sigma_ref_g(jl,jg))) * del2 +          &
                           ptemp_h(jl,jk,jg)* sin(pit2*                  &
                           (1._dp-pph(jl,jk,jg)/pps(jl,jg))/             &
                           (1._dp-sigma_ref_g(jl,jg))) * del1) / del

                      apos(4,jcell)=(ptemp_lin(jl,jk,jg)/                &
                           apos(3,jcell))**(1._dp/kappa)*100000._dp
                   elseif(apos(3,jcell) .ge. xi_crit(jl,jg) &
                        .and. jk .le. int(i_ref_g(jl,jg))) then
                      ptemp_lin(jl,jk,jg) = (ptemp_h(jl,jk+1,jg) * del2 +  &
                           ptemp_h(jl,jk,jg) * del1) / del
                      apos(4,jcell)=(ptemp_lin(jl,jk,jg) / &
                           apos(3,jcell))**(1._dp/kappa)            &
                           *100000._dp
                   endif
                   ! op_sb_20160708-

                   apos(4,jcell)=MIN(apos(4,jcell),pps(jl,jg))
                   apos(4,jcell)=MAX(apos(4,jcell),0._dp)
                END IF
             ENDIF
          CASE(3)
             if (ZLSTART) THEN
                APOS(4,JCELL) = APOS(3,JCELL) * apzero
             else
                APOS(4,JCELL) = APOS(3,JCELL) * PPS(JL,JG)
             endif
          END SELECT

          !
          !    B)  NUMBER OF CELLS IN BOX (*NCB*)
          !
          NCB(JL,JK,JG)= NCB(JL,JK,JG) + 1

          !
          !    C)  NUMBER OF CELLS IN BL  (*NCBL*)
          !
          NCBL(JL,JG)= MERGE(NCBL(JL,JG) + 1 , NCBL(JL,JG) &
               ,JK > KHPBL(JL,JG)  )
          !
       ENDDO !CELL-LOOP
       !
    END DO !LOOP OVER CHUNKS

    !  NECESSARRY FOR OUTPUT
    RNCB(:,:,:) = REAL(NCB(:,:,:),dp)
    RNCBL(:,:)  = REAL(NCBL(:,:), dp)

    IF (.NOT. L_PARALLEL_MODE) THEN
       GNCB(:,:,:) = RNCB(:,:,:)
       GNCBL(:,:)  = RNCBL(:,:)
    END IF

    !
    !* LT3.+4. LOOP OVER LEVELS + PBL.
    !          ---- ---- -------------
    !
    !*PDIR PARDO BY=1
    AMTOT  = SUM(PMBOX)
    AMCELL = AMTOT/REAL(NGCELL,dp)

!CDIR END

  END SUBROUTINE ATTILA_COUNT_CELLS
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_FIRSTL(PTSM1M,PGEOM1,ZTETA1,ZRI,TPOT)

    ! PURPOSE:
    !   - CALCAULATE THE FIRST FREE LAYER IN THE ATMOSPHERE (KHPBL)
    !     (REQUIRED BY ATTILA_FILL)
    !     (THIS IS AN ALTERNATIVE TO ATTILA_FIRSTL_EXT)
    !
    ! INTERFACE:
    !   - ATTILA_FIRSTL IS CALLED VIA MESSY-GLOBAL_START
    !
    ! AUTHOR(S):
    !   -  M. Traub, MPI Mainz, 2003


    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(IN)  :: PTSM1M(:,:)
    REAL(dp), INTENT(IN)  :: PGEOM1(:,:,:)
    REAL(dp), INTENT(IN)  :: ZTETA1(:,:,:)
    REAL(dp), INTENT(IN)  :: ZRI(:,:,:)
    REAL(dp), INTENT(IN)  :: TPOT(:,:,:)

    ! LOCAL
    REAL(dp) :: ZTETASURF(NLON,NGL)
    REAL(dp) :: ZTSTART(NLON,NGL)
    REAL(dp) :: ZGPBLMAX
    INTEGER  :: IHPBLRI(NLON,NGL)
    INTEGER  :: IHPBLAD(NLON,NGL)
    LOGICAL  :: LOPBLMAX
    INTEGER  :: JL, JR, JK

    ZGPBLMAX= 2500 * G  ! MAX BL HEIGHT OF 2500M (LAG. SCHEME)

    IHPBLRI(:,:)=NLEV
    IHPBLAD(:,:)=NLEV

    DO JR=1,NGL
       DO JL=1,NLON
          ZTETASURF(JL,JR)= PTSM1M(JL,JR)*                 &
               (100000._dp/PPH(JL,NLEVP1,JR))**(RD/CPD)

          ZTSTART(JL,JR)= TPOT(JL,NLEV,JR) +  &
               MERGE(1.2_dp,0.5_dp,ZTETASURF(JL,JR) >= ZTETA1(JL,NLEV,JR))
       END DO
    END DO

    DO JK=1,NLEV-1
       DO JL=1,NLON
          DO JR=1,NGL
             LOPBLMAX=(PGEOM1(JL,JK,JR)+PGEOM1(JL,JK+1,JR))*.5_dp  >  ZGPBLMAX
             IHPBLRI(JL,JR)= MERGE(JK,IHPBLRI(JL,JR), &
                  ZRI(JL,JK,JR) > 1.3_dp .OR. LOPBLMAX)
             IHPBLAD(JL,JR)= MERGE(JK,IHPBLAD(JL,JR),                        &
                  TPOT(JL,JK,JR) > ZTSTART(JL,JR)-PGEOM1(JL,JK,JR)/(RD/CPD) &
                  .OR. LOPBLMAX)
          END DO
       END DO
    END DO

    !  FINAL CALCULATIONS FOR LAGRANGIAN SCHEME
    !
    DO JL=1,NLON
       DO JR=1,NGL
          KHPBL(JL,JR) = MIN(IHPBLAD(JL,JR),IHPBLRI(JL,JR),NLEV-2)
       END DO
    END DO

  END SUBROUTINE ATTILA_FIRSTL
!--------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_FIRSTL_EXT(status, PBLH_IDX)

    ! PURPOSE:
    !   - SET FIRST LAYER IN FREE ATMOSPHERE (KHPBL) TO EXTERNALLY
    !     CALCULATED PLANETARY BOUNDARY LAYER HEIGHT INDEX
    !     (REQUIRED BY ATTILA_FILL)
    !     (THIS IS AN ALTERNATIVE TO ATTILA_FIRSTL)
    !
    ! AUTHOR(S):
    !   P. Joeckel, MPI Mainz, 2003
    !
    ! INTERFACE:
    !   - ATTILA_FIRSTL_EXT IS CALLED VIA MESSY-GLOBAL_START
    !

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(OUT) :: status
    ! INDEX OF PLANETARY BOUNDARY LAYER HEIGHT
    REAL(dp), INTENT(IN)  :: PBLH_IDX(:,:)

    ! LOCAL
    INTEGER :: JL, JR

    status = 1

    IF (SIZE(PBLH_IDX,1) /= NLON) RETURN
    IF (SIZE(PBLH_IDX,2) /= NGL)  RETURN

    DO JL=1,NLON
       DO JR=1,NGL
          KHPBL(JL,JR) = MIN(INT(PBLH_IDX(JL,JR)),NLEV-2)
       END DO
    END DO

    status = 0

  END SUBROUTINE ATTILA_FIRSTL_EXT
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_FILL(status)

    ! PURPOSE:
    !   FILL ATTILA MODULE VARIABLES WITH CURRENT WIND, SURFACE PRESSURE
    !   AND PRESSURE FIELD; CALCULATE MASS OF GRID BOXES
    !   called from attila_global_start
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTFILL)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_FILL)
    !      - VIRTUAL TEMPERATURE (ZTVM1) REMOVED FROM PARAMETERS
    !      - GEOPOPTENTIAL (PGEOM1) REMOVED FROM PARAMETERS

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    ! LOCAL
    INTEGER :: JL, JK, JG, JKK

    DO JG=1,NGL
       DO JK=1,NLEV
          DO JL=1,NLON
             PMBOX(JL,JK,JG)= (PPH(JL,JK+1,JG)-PPH(JL,JK,JG) ) * AARBOX(JG) / G
          END DO
          DO JL=1,NLON
             JKK= KHPBL(JL,JG)
             PMBOXBL(JL,JG)= ( PPH(JL,NLEVP1,JG)-PPH(JL,JKk+1,JG) ) &
                  * AARBOX(JG)/G
          END DO
       END DO

    END DO

    CALL ATTILA_INIT_GETA(status, .FALSE.)
    IF (status /= 0) RETURN

    status = 0

  END SUBROUTINE ATTILA_FILL
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_DRIVE(&
        HARVESTTURB, HARVESTMC &    ! PBL-TURB, MC
       ,HARVESTCAT , TI2, lstart)           ! CAT

    ! PURPOSE:
    !   - CONTROLS THE ADVECTION FOR THE LAGRANGIAN SCHEME, INCLUDING
    !     PLANETARY BOUNDARY LAYER TURBULENCE, MONTE CARLO DIFFUSION
    !     AND CLEAR AIR TURBULENCE (THE LATTER IS EXPERIMENTAL!)
    !
    ! INTERFACE:
    !   - LTDRIVE WAS INLINED INTO SCAN1SL
    !   - ATTILA_DRIVE IS CALLED VIA attila_integrate
    !   - attila_integrate IS CALLED VIA attila_global_start
    !
    ! EXTERNALS:
    !   - ATTILA_CAT_MOVE
    !   - ATTILA_RKO4     ADVECTION (RUNGE-KUTTA, 4TH ORDER)
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTDRIVE)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_DRIVE)
    !      - PARAMETERS CHANGED
    !      - CALL TO LTINICON REMOVED (INITI. CONSTITUENT CONCENTRATIONS
    !                                  IN CELLS)
    !      - CALLS TO LTINIRAN, LTSMOOTH, RESETI, RESETL, RESETR, COPYI
    !        RANDNUM REMOVED
    !
    ! REMARKS:
    ! - ORIGINAL PARAMETERS OF LTDRIVE:
    !            I HARVESTINI               ! RANDOM NUMBER ARRAY - OPTIONAL
    !            I HARVESTIRANPOINT         ! RANDOM NUMBER ARRAY - OPTIONAL
    !            I HARVESTINIPOS            ! RANDOM NUMBER ARRAY - OPTIONAL
    !            I HARVESTINIPOS2           ! RANDOM NUMBER ARRAY - OPTIONAL
    !            O RNCB                     ! NUMBER OF CELLS IN GRID BOXES
    !            O RNCBL                    ! NUMBER OF CELLS IN BL
    !
    ! - ORIGINAL LOCAL FIELDS:
    !            ZMBOX* (     NLON,NGL     ) MASS OF GRID BOXES
    !            IRANPOINT*(4)         RANOM NUMBER FROM RANDOM NUMBER ARRAY
    !            ZRAN*     (NCEVL,4)       CORRESPONDING RANDOM NUMBER
    !            ZPOSV* (NCEVL,NAPOS)      A CHUNK OF APOS
    !            IPOSV* (NCEVL,NNPOS)      A CHUNK OF NPOS
    !
    ! - ORIGINAL METHODS:
    !      0. GENERAL REMARKS
    !         A) FOR MONTE-CARLO DIFFUSION (RANDOM DISPLACEMENT) AND
    !            FOR BOUNDARY LAYER TURBULENCE A LOT OF RANDOM NUMBERS ARE
    !            NEEDED. TO SAVE CPU TIME ARRAYS WITH RANDOM NUMBERS ARE
    !            COMPUTED ONCE EVERY RUN.
    !            EVERY TIME STEP A RANDOM START POSITION IN THE ARRAY IS
    !            DETERMINED (RANDOM POINTER) AND THE CORRESPONDING NUMBER
    !            FOR A CELL IS FOUND IN THE ARRAY AT POSITION
    !            "START POS. + CELL NUMBER".
    !            THIS IS DONE FOR 4 RANDOM NUMBERS:
    !             - 3 NUMBERS RUN FROM -1...1 (ARRAY *ARANDISP*), AND ARE
    !               NEEDED FOR MONTE-CARLO DIFF. (U-/V-/W-DIRECTION)
    !             - 1 NUMBER RUNS FROM 0...1  (ARRAY *ARANBLTURB*), AND IS
    !               NEEDED FOR BL TURBULENCE.
    !            NB: SINCE THE ARRAYS WITH RANDOM NUMBERS ARE COMPUTED
    !                ANEW AT EVERY RESTART, IT DEPENDS ON THE FREQUENCY OF
    !                A RESTART WHAT NUMBERS ARE ACTUALLY USED. THIS IS
    !                ADMITTEDLY NOT IDEAL, BUT THE RESULTS SHOULD NOT BE
    !                SENSITIVE TO WHAT NUMBERS ARE ACTUALLY USED, ANYWAY.
    !                WHEN NOT ALTERING THE FREQUENCY OF RESTART, RESULTS
    !                ARE REPRODUCABLE.
    !
    !      1. INITIALISATION
    !         IF NSTEP=NSTART      ARRAYS FROM (A.1.)+(A.2.) ARE INITIALISED
    !         IN ANY CASE     SOME ARRAYS FROM (B.2.)        ARE INITIALISED
    !
    !      2. 1ST LOOP OVER CELLS
    !         THE FOLLOWING CALCULATIONS ARE DONE:
    !         # NUMBER OF CELLS IN BOX (*NCB*)
    !         # NUMBER OF CELLS IN BL  (*NCBL*)
    !
    !   3.+4. LOOP OVER LEVELS PLUS PBL
    !
    !      5. 2ND LOOP OVER CELLS
    !         THE FOLLOWING CALCULATIONS ARE DONE:
    !         A) PHYSICS II
    !            I.   MIXING BY DIFFUSION
    !         B) DYNAMICS (CALCULATION OF NEW *APOS*,*NPOS*)
    !                     (DONE IN *LTRKO4*)
    !            I.   ADVECTION BY LARGE SCALE WIND
    !            II.  MONTE-CARLO DIFFUSION
    !            III. BL TURBULENCE
    !
    !     REMARK:
    !          THE TWO LOOPS OVER THE CELLS (2. & 5.) ARE SPLITTED INTO
    !          *NCHUNK* CHUNKS WHICH ARE HANDLED IN PARALLEL. IN EACH CHUNK
    !          CALCULATIONS FOR NCELL/NCHUNK = NCEVL CELLS ARE PERFORMED.
    !          THE LOOPS OVER THESE *NCEVL* CELLS ARE VECTORIZED
    !          WHENEVER POSSIBLE.

    IMPLICIT NONE

    !I/O
    ! RANDOM NUMBER ARRAY FOR PBL-TURB
    REAL(dp), DIMENSION(:),     INTENT(IN)       :: HARVESTTURB
    ! RANDOM NUMBER ARRAY FOR MONTE CARLO DIFFUSION
    REAL(dp), DIMENSION(:,:),   INTENT(IN)       :: HARVESTMC
    ! RANDOM NUMBER ARRAY FOR CLEAR AIR TURBULENCE
    REAL(dp), DIMENSION(:),     INTENT(IN)       :: HARVESTCAT
    ! TURBULENCE INDEX
    REAL(dp), DIMENSION(:,:,:), POINTER          :: TI2
    LOGICAL,  INTENT(IN)                         :: lstart

    ! LOCAL
    REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZRAN       ! NCEVL x 4
    REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZPOSV      ! NAPOS x NCELL
    INTEGER,  DIMENSION(:,:), ALLOCATABLE  :: IPOSV      ! NNPOS x NCELL
    INTEGER                                :: JC,JCHUNK, JCELL
    INTEGER                                :: JCELL1, JCELL2, JCELL3
    INTEGER                                :: JCEVL, jl,jg,jk
    REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZVELV      ! NCELL x 3
    REAL(dp), DIMENSION(NLON,NLEVM1,NGL)   :: cfl

    ALLOCATE(ZRAN(NCEVL,4))
    ALLOCATE(ZPOSV(NCEVL,NAPOS))
    ALLOCATE(IPOSV(NCEVL,NNPOS))
    IF (LVDIAG) ALLOCATE(ZVELV(NCEVL,3))

    ! SAVE NCB OF t-1
    NCBM1(:,:,:) = NCB(:,:,:)
    cfl  (:,:,:) = 0._dp

! op_sb_20170110+
    if (I_VERT == 2 .and. .not. LSTART) then ! control CFL criterium
     do jl = 1, nlon
       do jk= 1, nlevm1
         do jg = 1, ngl
         cfl(jl,jk,jg) = dtime/(getaf_3d(jl,jk,jg)-getaf_3d(jl,jk+1,jg)) &
                            *pww(jl,jk,jg)
         if (cfl(jl,jk,jg) > 1._dp) then      ! correct velocity
           pww(jl,jk,jg) = (getaf_3d(jl,jk,jg)-getaf_3d(jl,jk+1,jg)) / dtime
         elseif (cfl(jl,jk,jg) < -1._dp) then
           pww(jl,jk,jg) = -(getaf_3d(jl,jk,jg)-getaf_3d(jl,jk+1,jg)) / dtime
         endif
         enddo
       enddo
     enddo
    endif
! op_sb_20170110-
    !
    !* LT5. 2ND LOOP OVER CELLS.
    !       --- ---- ---- ------
    !
    !*PDIR PARDO BY=1
    DO JCHUNK=0,NCHUNKM1
       JCELL1= JCHUNK*NCEVL
       JCELL2= JCELL1 + 1
       JCELL3= JCELL1 + NCEVL

       !    B) DYNAMICS (CALCULATION OF NEW *APOS*,*NPOS*)
       !       I.   ADVECTION BY LARGE SCALE WIND
       !       II.  MONTE-CARLO DIFFUSION
       !       III. BL TURBULENCE
       !

       DO JCEVL=1,NCEVL
          ZRAN(JCEVL,:) = 0.0_dp
       END DO

       IF (LLTBLTURB) THEN
          DO JCEVL=1,NCEVL
             JCELL= JCELL1+JCEVL
             ZRAN(JCEVL,4)= HARVESTTURB(JCELL)  ! e [0,1]
          END DO
       END IF

       IF (LLTMOCA) THEN
          DO JCEVL=1,NCEVL
             JCELL= JCELL1+JCEVL
             ZRAN(JCEVL,1)= 2.0_dp * HARVESTMC(JCELL,1) - 1.0_dp   ! e [-1,1]
             ZRAN(JCEVL,2)= 2.0_dp * HARVESTMC(JCELL,2) - 1.0_dp
             ZRAN(JCEVL,3)= 2.0_dp * HARVESTMC(JCELL,3) - 1.0_dp
          END DO
       END IF

       DO JC= 1,NAPOS
          DO JCEVL=1,NCEVL
             JCELL= JCELL1+JCEVL
             ZPOSV(JCEVL,JC)= APOS(JC,JCELL)
          ENDDO
       ENDDO

       DO JC= 1,NNPOS
          DO JCEVL=1,NCEVL
             JCELL= JCELL1+JCEVL
             NPOSM1(JC,JCELL)= NPOS(JC,JCELL)
             IPOSV(JCEVL,JC)= INT(NPOS(JC,JCELL))
          ENDDO
       ENDDO

       ! SAVE THE PRESSURE HEIGHT OF EACH CELL AT T-1
       APOSM1_PRESS(:)=APOS(4,:)

       ! TURBULENT MOVEMENT OF AIR PARCELS DUE TO CLEAR AIR TURBULENCE
       IF (LLCAT) &
            CALL ATTILA_CAT_MOVE(ZPOSV(1,3),IPOSV,HARVESTCAT,TI2)

       IF (LVDIAG) THEN
          CALL ATTILA_RKO4(ZPOSV,IPOSV,ZRAN,ZVELV)
       ELSE
          CALL ATTILA_RKO4(ZPOSV,IPOSV,ZRAN)
       END IF

       DO JC= 1,NAPOS
          DO JCEVL=1,NCEVL
             JCELL= JCELL1+JCEVL
             APOS(JC,JCELL)= ZPOSV(JCEVL,JC)
          ENDDO
       ENDDO

       DO JC= 1,NNPOS
          DO JCEVL=1,NCEVL
             JCELL= JCELL1+JCEVL
             NPOS(JC,JCELL)= REAL(IPOSV(JCEVL,JC),dp)
          ENDDO
       ENDDO
       !
       IF (LVDIAG) THEN
          DO JCEVL=1,NCEVL
             JCELL= JCELL1+JCEVL
             uvel(JCELL) = ZVELV(JCEVL,1)
             vvel(JCELL) = ZVELV(JCEVL,2)
             wvel(JCELL) = ZVELV(JCEVL,3)
          ENDDO
       END IF
    END DO       !LOOP OVER CHUNKS

    DEALLOCATE(ZRAN)
    DEALLOCATE(ZPOSV)
    DEALLOCATE(IPOSV)
    IF (ALLOCATED(ZVELV)) DEALLOCATE(ZVELV)

  END SUBROUTINE ATTILA_DRIVE
! -------------------------------------------------------------------

! -------------------------------------------------------------------
   SUBROUTINE ATTILA_CONVGPC(UMASSFL, DMASSFL, TYPECONV, STATUS)

     IMPLICIT NONE

     ! IN/OUT
     REAL(dp), INTENT(IN)  :: UMASSFL(:,:,:)         ! UPWARD MASS FLUX
     REAL(dp), INTENT(IN)  :: DMASSFL(:,:,:)         ! DPWARD MASS FLUX
     REAL(dp), INTENT(IN)  :: TYPECONV(:,:)          ! TYPE OF CONVECTION
     INTEGER,  INTENT(OUT) :: STATUS

     ! LOCAL
     INTEGER :: JK, JL, JM
     REAL(dp)              :: ZFDELTAU, ZFDELTAD

     STATUS=1

     ! SET INITIAL VALUES TO ZERO
     DO JK= 1,NLEV
        DO JL= 1,NLON
           DO JM= 1, NGL
              PPENDEU(JL,JK,JM)= 0._dp
              PPENDED(JL,JK,JM)= 0._dp
              PPSUBS(JL,JK,JM)= 0._dp
           ENDDO
        ENDDO
     ENDDO

     DO JL= 1,NLON
        DO JM= 1,NGL
           KTYPELT(JL,JM)= TYPECONV(JL,JM)
           KCTOPLT(JL,JM)= NLEV
        ENDDO
     ENDDO

     ! CONSISTENCY CHECK
     DO JK= 1,NLEV
        DO JL= 1,NLON
           DO JM= 1, NGL
              IF(UMASSFL(JL,JK,JM)  <  0._dp) RETURN
              IF(DMASSFL(JL,JK,JM)  >  0._dp) RETURN
           ENDDO
        ENDDO
     ENDDO

     !* EN-/DETRAINMENT CALCULATION
     DO JK= 1,NLEV-1
        DO JL= 1,NLON
           DO JM= 1,NGL
              ZFDELTAU= UMASSFL(JL,JK,JM)-UMASSFL(JL,JK+1,JM)
              ZFDELTAD= DMASSFL(JL,JK,JM)-DMASSFL(JL,JK+1,JM)
              ZFDELTAD= MERGE(0._dp,ZFDELTAD, ABS(ZFDELTAD)  <  1.E-10_dp)

              ! TEST FOR TOP LEVEL OF CONVECTION
              KCTOPLT(JL,JM)= MERGE(JK, KCTOPLT(JL,JM), KCTOPLT(JL,JM) &
                    ==  NLEV .AND. ZFDELTAU .NE. 0._dp)
              !
              ! PROBABILITY
              !
              IF(ZFDELTAU  >=  0._dp) THEN
                 PPENDEU(JL,JK,JM)= (ZFDELTAU*DTIME*G)/ &
                      (PPH(JL,JK+1,JM)-PPH(JL,JK,JM))
              ELSE
                 PPENDEU(JL,JK,JM)= ZFDELTAU/UMASSFL(JL,JK+1,JM)
              ENDIF
              !
              IF(ZFDELTAD  >=  0._dp) THEN
                 PPENDED(JL,JK,JM)= (ZFDELTAD*DTIME*G)/ &
                      (PPH(JL,JK+1,JM)-PPH(JL,JK,JM))
              ELSE
                 PPENDED(JL,JK,JM)= ZFDELTAD/(-1._dp*DMASSFL(JL,JK,JM))
              ENDIF
           ENDDO
        ENDDO
     ENDDO

     ! GROUND LEVEL
     DO JL= 1,NLON
        DO JM= 1,NGL
           PPENDEU(JL,NLEV,JM)= (UMASSFL(JL,NLEV,JM)*DTIME*G)/ &
                (PPH(JL,NLEVP1,JM)-PPH(JL,NLEV,JM))
           PPENDED(JL,NLEV,JM)= -1._dp
        ENDDO
     ENDDO

     !
     ! SUBSIDENCE
     !
     PPSUBS(:,:,:)=0._dp

     ! ERROR STATUS
     STATUS=0

   END SUBROUTINE ATTILA_CONVGPC
! -------------------------------------------------------------------

! -------------------------------------------------------------------
   SUBROUTINE ATTILA_CONVTRAJ(IUPDO, JTOFFSET, HARVESTCONV &
        , ZPOS_BEG_UP, ZPOS_END_UP, ZPOS_END_DO, ZPOS_END_SUB, L_EXIT)

     ! NOTE: THE VECTORS ZPOS_* ARE USED TO SAVE INFORMATION ABOUT THE
     !       CONVECTIVE MOVEMENT OF THE AIR PARCELS WITHIN ONE TIME STEP
     !       (FOR LATER USE IN OTHER SUBMODELS, SUCH AS FOR INSTANCE THE
     !        CALCULATION OF CONVECTIVE H2O-PHASE TRANSITIONS ...)
     !       THE ORDER IS AS FOLLOWS:
     !            ZPOS_BEG_UP: VERTICAL POSITION AT THE BEGINNING OF THE
     !                         UPDRAFT
     !            -> UPDRAFT CALCULATION
     !            ZPOS_END_UP: VERTICAL POSITION AT THE END OF THE UPDRAFT
     !            -> DOWNDRAFT CALCULATION
     !            ZPOS_END_DO: VERTICAL POSITION AT THE END OF THE DOWNDRAFT
     !            -> LARGE SCALE SUBSIDENCE
     !            ZPOS_END_SUB: VERTICAL POSITION AT THE END OF THE SUBSIDENCE
     !
     !      ZPOS_END_SUB - ZPOS_BEG_UP SHOULD THEN BE EQUAL TO THE TOTAL
     !      VERTICAL MOVEMENT ...
     !

     IMPLICIT NONE

     ! I/O
     INTEGER,  INTENT(IN)                        :: IUPDO    ! PROCESS SWITCH
     INTEGER,  INTENT(IN)                        :: JTOFFSET ! CELL-# OFFSET
     LOGICAL,  INTENT(OUT)                       :: L_EXIT
     REAL(dp), INTENT(IN)                        :: HARVESTCONV(:)
     REAL(dp), INTENT(INOUT)                     :: ZPOS_BEG_UP(:)
     REAL(dp), INTENT(INOUT)                     :: ZPOS_END_UP(:)
     REAL(dp), INTENT(INOUT)                     :: ZPOS_END_DO(:)
     REAL(dp), INTENT(INOUT)                     :: ZPOS_END_SUB(:)

     ! LOCAL
     INTEGER,  PARAMETER                           :: I_CONV_ITER_MAX = 2000
     INTEGER,                                 SAVE :: NITER
     INTEGER,  DIMENSION(:,:),   ALLOCATABLE, SAVE :: KPOS  ! CASE 1 -> 2
     INTEGER,  DIMENSION(:),     ALLOCATABLE, SAVE :: JRCNT ! CASE 1 -> 2
     INTEGER,  DIMENSION(:),     ALLOCATABLE, SAVE :: JVPOS ! CASE 1 -> 2
     REAL(dp), DIMENSION(:),     ALLOCATABLE, SAVE :: CPOS  ! CASE 1 -> 2
     REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: PPENDES  ! CASE 2
     REAL(dp), DIMENSION(:,:),   ALLOCATABLE, SAVE :: ZTEST2   ! CASE 2
     REAL(dp), DIMENSION(:,:),   ALLOCATABLE, SAVE :: ZTEST3   ! CASE 2
     REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: NCB_SUBS ! CASE 2
     LOGICAL,  DIMENSION(:,:,:), ALLOCATABLE, SAVE :: LONEG    ! CASE 2
     LOGICAL,  DIMENSION(:,:),   ALLOCATABLE, SAVE :: L_SUB    ! CASE 2
     INTEGER :: JT, JL, JG, JK, JC, JTG
     LOGICAL :: LOCO(NCELL), LOGOONU(NCELL), LOGOOND(NCELL)
     REAL(dp):: ZTEST
     REAL(dp):: sc ! op_sb_20160808

     ! op_sb_20160808+
     SELECT CASE(I_VERT)
     CASE(1,3)
        sc = 1.0_dp
     CASE(2)
        sc = -1.0_dp
     END SELECT
     ! op_sb_20160808-

     SELECT CASE (IUPDO)

     CASE(1)
        !
        ! --- PREAPARE LOCAL FIELDS -------------------------------------------
        !
        ! MEMORY
        ! (1) -> (2)
        ALLOCATE(KPOS(NCELL,NNPOS))
        ALLOCATE(JRCNT(NCELL))
        ALLOCATE(JVPOS(NCELL))
        ALLOCATE(CPOS(NCELL))            ! VERTICAL POSITION
        !
        !
        ! INIT
        L_EXIT = .FALSE.
        NITER = 0
        !
        DO JC= 1,NNPOS
           DO JT=1,NCELL
              KPOS(JT,JC)= INT(NPOS(JC,JT))
           END DO
        END DO
        !
        CPOS(:)= APOS(3,:)
        !
        ! INIT DIAGNOSTIC
        ZPOS_BEG_UP(:)  = NPOS(3,:)
        !
        ! --- UPDRAFT/DOWNDRAFT ----------------------------------------------
        !
        ! UPDRAFT
        !
        DO JT= 1,NCELL
           JTG=JT+JTOFFSET
           JL= KPOS(JT,4)
           JG= KPOS(JT,5)
           LOCO(JT)= (KTYPELT(JL,JG) /= 0._dp)
           JRCNT(JT)= JTG
           JK= KPOS(JT,3)
           ZPOS_END_UP(JT) = REAL(JK,dp)
           IF(PPENDEU(JL,JK,JG) > 0 .AND. &
                HARVESTCONV(JTG) <= PPENDEU(JL,JK,JG) .AND. LOCO(JT)) THEN
              JVPOS(JT)= JK-1
              LOGOONU(JT)= .TRUE.
              ZPOS_END_UP(JT) = REAL(JVPOS(JT),dp)
           ELSE
              JVPOS(JT)= JK
              LOGOONU(JT)= .FALSE.
           END IF
        END DO

        DO JK= NLEV,1,-1
           DO JT= 1,NCELL
              !JTG=JT+JTOFFSET
              JL= KPOS(JT,4)
              JG= KPOS(JT,5)
              JRCNT(JT)= MERGE(JRCNT(JT)+1, 1, JRCNT(JT) < NGCELL)
              IF (JK == JVPOS(JT)) THEN
                 IF(LOGOONU(JT) .AND. ((PPENDEU(JL,JK,JG) >= 0._dp) .OR. &
                      ((PPENDEU(JL,JK,JG) < 0._dp) .AND.                 &
                      -1.0_dp * HARVESTCONV(JRCNT(JT)) <=                &
                      PPENDEU(JL,JK,JG)))) THEN
                    JVPOS(JT)= JVPOS(JT)-1
                    ZPOS_END_UP(JT) = REAL(JVPOS(JT),dp)
                 ELSE
                    LOGOONU(JT)= .FALSE.
                 END IF
              END IF
           END DO
        ENDDO

        !
        ! DOWNDRAFT
        !
        DO JT= 1,NCELL
           JTG=JT+JTOFFSET
           JL= KPOS(JT,4)
           JG= KPOS(JT,5)
           LOCO(JT)= (KTYPELT(JL,JG) /= 0._dp)
           JRCNT(JT)= JTG
           JK=JVPOS(JT)
           ZPOS_END_DO(JT) = REAL(JK,dp)
           IF((PPENDED(JL,JK,JG) > 0) .AND. &
                (HARVESTCONV(JTG) <= PPENDED(JL,JK,JG)) &
                .AND. LOCO(JT) .AND. (JK < NLEV)) THEN
              JVPOS(JT)= JK+1
              LOGOOND(JT)= .TRUE.
              ZPOS_END_DO(JT) = REAL(JVPOS(JT),dp)
           ELSE
              JVPOS(JT)= JK
              LOGOOND(JT)= .FALSE.
           END IF
        END DO
        !
        DO JK= 1,NLEV-1
           DO JT= 1,NCELL
              !JTG=JT+JTOFFSET
              JL= KPOS(JT,4)
              JG= KPOS(JT,5)
              JRCNT(JT)= MERGE(JRCNT(JT)+1, 1, JRCNT(JT) < NGCELL)
              IF (JK == JVPOS(JT)) THEN
                 IF(LOGOOND(JT) .AND. ((PPENDED(JL,JK,JG) >= 0._dp) .OR.   &
                      ((PPENDED(JL,JK,JG) < 0._dp) .AND.                   &
                   (-1.0_dp * HARVESTCONV(JRCNT(JT)) <= PPENDED(JL,JK,JG))))) &
                      THEN
                    JVPOS(JT)= JVPOS(JT)+1
                    ZPOS_END_DO(JT) = REAL(JVPOS(JT),dp)
                 ELSE
                    LOGOOND(JT)= .FALSE.
                 END IF
              END IF
           END DO
        END DO

        !
        ! APPLYING NEW POSITION
        !
        DO JT= 1,NCELL
           !JTG=JT+JTOFFSET
           GETAH  => GETAH_3D(KPOS(JT,4),:,KPOS(JT,5))
           GDETAH => GDETAH_3D(KPOS(JT,4),:,KPOS(JT,5))
           IF (JVPOS(JT) /= KPOS(JT,3)) THEN
              CPOS(JT)= GETAH(JVPOS(JT))+ &
                   sc*HARVESTCONV(JRCNT(JT))*GDETAH(JVPOS(JT))
           END IF
           JRCNT(JT)= MERGE(JRCNT(JT)+1, 1, JRCNT(JT) < NGCELL)
           KPOS(JT,3)=GINDEX_V(CPOS(JT),GETAH) ! op_sb_2014021
           NPOS(3,JT)=REAL(KPOS(JT,3),dp)
           APOS(3,JT)=CPOS(JT)
        END DO

     CASE(2)
        !
        ! --- PREAPARE LOCAL FIELDS -------------------------------------------
        !
        ! MEMORY
        ! allocate before first iteration
        IF (NITER == 0) THEN
           ALLOCATE(PPENDES(NLON,NLEV,NGL))
           ALLOCATE(ZTEST2(NLON,NGL))
           ALLOCATE(ZTEST3(NLON,NGL))
           ALLOCATE(NCB_SUBS(NLON,NLEV,NGL))
           ALLOCATE(LONEG(NLON,NLEV,NGL))
           ALLOCATE(L_SUB(NLON,NGL))
        END IF
        !
        ! CHECK IF INITIAL PARCEL DISTRIBUTION IS REACHED AGAIN -----------
        L_EXIT = .TRUE.
        L_SUB(:,:) = .FALSE.
        DO JG = 1,NGL
           DO JL = 1,NLON
              DO JK = 1,NLEV
                 ZTEST = GNCB_MOVE(JL,JK,JG) - GNCB(JL,JK,JG)
                 L_EXIT = L_EXIT .AND. (ZTEST == 0._dp)
                 L_SUB(JL,JG) = L_SUB(JL,JG) .OR. (ZTEST /= 0._dp)
              END DO
           END DO
        END DO
        !
        L_EXIT = L_EXIT .OR. (NITER > I_CONV_ITER_MAX)
        !
        IF (L_EXIT) THEN
           !
           ! --- CLEAN UP ----------------------------------------------------
           !
           ! CASE(1) -> CASE(2)
           DEALLOCATE(CPOS)
           DEALLOCATE(KPOS)
           DEALLOCATE(JVPOS)
           DEALLOCATE(JRCNT)

           ! ONLY CASE(2)
           DEALLOCATE(PPENDES)
           DEALLOCATE(ZTEST2)
           DEALLOCATE(ZTEST3)
           DEALLOCATE(NCB_SUBS)
           DEALLOCATE(LONEG)
           DEALLOCATE(L_SUB)

           RETURN
        END IF

        ! INIT
        NITER = NITER + 1
        !
        DO JC= 1,NNPOS
           DO JT=1,NCELL
              KPOS(JT,JC)= INT(NPOS(JC,JT))
           END DO
        END DO
        !
        NCB_SUBS(:,:,:) = 0._dp
        PPENDES(:,:,:)  = 0._dp
        LONEG(:,:,:)    = .FALSE.
        ZTEST2(:,:) = 0._dp
        ZTEST3(:,:) = 0._dp

        !
        ! --- SUBSIDENCE --------------------------------------------------
        !
        DO JK = 1,NLEV
           DO JG = 1,NGL
              DO JL = 1,NLON
                 IF (KTYPELT(JL,JG) /= 0._dp) THEN
                    ! diff parcels old-new
                    PPSUBS(JL,JK,JG)   = &
                         GNCB(JL,JK,JG) - GNCB_MOVE(JL,JK,JG)
                 ELSE
                    PPSUBS(JL,JK,JG)   = 0._dp
                 END IF
              END DO
           END DO
        END DO

        ! 1. number of parcels sinking per level
        NCB_SUBS(:,1,:) = PPSUBS(:,1,:)
        !
        DO JK = 2,NLEV
           DO JG = 1,NGL
              DO JL = 1,NLON
                 IF (KTYPELT(JL,JG) /= 0._dp) THEN
                    NCB_SUBS(JL,JK,JG) = &
                         NCB_SUBS(JL,JK-1,JG) + PPSUBS(JL,JK,JG)
                 ELSE
                    NCB_SUBS(JL,JK,JG) = 0._dp
                 END IF
                 ! test if at a column any parcel was in updraft
                 IF (KTYPELT(JL,JG) /= 0._dp) THEN
                    ZTEST2(JL,JG) = ZTEST2(JL,JG) + NCB_SUBS(JL,JK,JG)
                 END IF
              END DO
           END DO
        END DO

        ! 2. number of parcels raise per level (downdraft without updraft)
        DO JG = 1,NGL
           DO JL = 1,NLON
              IF (ZTEST2(JL,JG) > 0._dp) THEN
                 NCB_SUBS(JL,NLEV,JG) = PPSUBS(JL,NLEV,JG)
                 ZTEST3(JL,JG)   = NCB_SUBS(JL,NLEV,JG)
              END IF
           END DO
        END DO

        DO JK = NLEV-1,1,-1
           DO JG = 1,NGL
              DO JL = 1,NLON
                 IF (ZTEST2(JL,JG) > 0._dp)THEN
                    NCB_SUBS(JL,JK,JG) = &
                         NCB_SUBS(JL,JK+1,JG) + PPSUBS(JL,JK,JG)
                    ZTEST3(JL,JG) = ZTEST3(JL,JG) + NCB_SUBS(JL,JK,JG)
                 END IF
              END DO
           END DO
        END DO

        ! 3. number of parcels in subsidence (up and down at a column)
        DO JG = 1,NGL
           DO JL= 1,NLON
              IF (ZTEST2(JL,JG) > 0._dp .AND. ZTEST3(JL,JG) > 0._dp) THEN
                 IF (PPSUBS(JL,1,JG) < 0._dp) THEN
                    NCB_SUBS(JL,1,JG) = PPSUBS(JL,1,JG)
                 ELSE
                    NCB_SUBS(JL,1,JG) = 0._dp
                 END IF
              END IF
           END DO
        END DO

        DO JK =2,NLEV
           DO JG = 1,NGL
              DO JL = 1,NLON
                 IF (ZTEST2(JL,JG) > 0._dp .AND. ZTEST3(JL,JG) > 0._dp) THEN
                    IF (PPSUBS(JL,JK,JG) < 0._dp .OR. LONEG(JL,JK-1,JG)) THEN
                       LONEG(JL,JK,JG) = .TRUE.
                       NCB_SUBS(JL,JK,JG) = &
                            NCB_SUBS(JL,JK-1,JG) + PPSUBS(JL,JK,JG)
                    ELSE
                       NCB_SUBS(JL,JK,JG) = 0._dp
                    END IF
                 END IF
              END DO
           END DO
        END DO

        ! 1. and 3. probabilities for subsidence
        DO JK = 2,NLEV
           DO JG = 1,NGL
              DO JL = 1,NLON
                 IF ((ZTEST2(JL,JG) <= 0._dp) .OR.                    &
                      (ZTEST2(JL,JG) > 0._dp  .AND. ZTEST3(JL,JG) > 0._dp)) THEN
                    ! JK => JK-1
                    ZTEST=GNCB_MOVE(JL,JK,JG)+ABS(NCB_SUBS(JL,JK-1,JG))
                    IF ( (KTYPELT(JL,JG) /= 0._dp) .AND. (ZTEST /= 0._dp) ) THEN
                       PPENDES(JL,JK,JG) = NCB_SUBS(JL,JK,JG)/ZTEST
                    ELSE
                       PPENDES(JL,JK,JG) =  0._dp
                    END IF
                 END IF
              END DO
           END DO
        END DO

        ! 2. probabilities for subsidence (downdraft without updraft)
        DO JG = 1,NGL
           DO JL = 1,NLON
              IF (ZTEST2(JL,JG) > 0._dp .AND. ZTEST3(JL,JG) < 0._dp) THEN
                 IF (KTYPELT(JL,JG) /= 0._dp &
                      .AND. GNCB_MOVE(JL,NLEV,JG) /= 0._dp) THEN
                    PPENDES(JL,NLEV,JG) = &
                         NCB_SUBS(JL,NLEV,JG)/GNCB_MOVE(JL,NLEV,JG)
                 ELSE
                    PPENDES(JL,NLEV,JG) =  0._dp
                 END IF
              ENDIF
           END DO
        END DO

        DO JK = NLEV-1,1,-1 ! 1,NLEV => NLEV-1,1
           DO JG = 1,NGL
              DO JL = 1,NLON
                 IF (ZTEST2(JL,JG) > 0._dp .AND. ZTEST3(JL,JG) < 0._dp)THEN
                    ! JK => JK+1
                    ZTEST=GNCB_MOVE(JL,JK,JG)+ABS(NCB_SUBS(JL,JK+1,JG))
                    IF ( (KTYPELT(JL,JG) /= 0._dp) .AND. (ZTEST /= 0._dp) ) THEN
                       PPENDES(JL,JK,JG) = NCB_SUBS(JL,JK,JG)/ZTEST
                    ELSE
                       PPENDES(JL,JK,JG) =  0._dp
                    END IF
                 ENDIF
              END DO
           END DO
        END DO

        ZPOS_END_SUB(:) = REAL(JVPOS(:),dp)

        ! parcels move 1. and 3. -> ZTEST2 < 0 (OR ZTEST2 > 0 AND ZTEST3 > 0)
        DO JK= 1,NLEV-1
           DO JT= 1,NCELL
              !JTG=JT+JTOFFSET
              JL= KPOS(JT,4)
              JG= KPOS(JT,5)
              LOCO(JT)= (KTYPELT(JL,JG) /= 0._dp) .AND. L_SUB(JL,JG)
              JRCNT(JT)= MERGE(JRCNT(JT)+1, 1, JRCNT(JT) < NGCELL)
              IF((JK == JVPOS(JT)) .AND. LOCO(JT)) THEN
                 IF( ZTEST2(JL,JG) <= 0._dp .OR.                    &
                      (ZTEST2(JL,JG) >  0._dp .AND. ZTEST3(JL,JG) > 0._dp)) THEN
                    IF( (PPENDES(JL,JK,JG) < 0._dp) .AND.   &
                       (-1.0_dp*HARVESTCONV(JRCNT(JT)) >= PPENDES(JL,JK,JG)) ) &
                         THEN
                       JVPOS(JT) = JVPOS(JT) + 1
                       ZPOS_END_SUB(JT) = REAL(JVPOS(JT),dp)
                    END IF
                 END IF
              END IF
           END DO
        END DO

        ! parcels move 2. -> ZTEST2 > 0 .AND. ZTEST3 < 0
        DO JK = NLEV,2,-1
           DO JT= 1,NCELL
              !JTG=JT+JTOFFSET
              JL= KPOS(JT,4)
              JG= KPOS(JT,5)
              LOCO(JT)= (KTYPELT(JL,JG) /= 0._dp) .AND. L_SUB(JL,JG)
              JRCNT(JT)= MERGE(JRCNT(JT)+1, 1, JRCNT(JT) < NGCELL)
              IF ( JK == JVPOS(JT) .AND. LOCO(JT)) THEN
                 IF( ZTEST2(JL,JG) > 0._dp .AND. ZTEST3(JL,JG) < 0._dp ) THEN
                    IF( (PPENDES(JL,JK,JG) < 0._dp) .AND.   &
                      (-1.0_dp*HARVESTCONV(JRCNT(JT)) >= PPENDES(JL,JK,JG)) ) &
                         THEN
                       JVPOS(JT) = JVPOS(JT) - 1
                       ZPOS_END_SUB(JT) = REAL(JVPOS(JT),dp)
                    END IF
                 END IF
              END IF
           END DO
        END DO
        !
        CPOS(:)= APOS(3,:)

        ! APPLYING NEW POSITION
        !
        DO JT= 1,NCELL
           !JTG=JT+JTOFFSET
           GETAH  => GETAH_3D(KPOS(JT,4),:,KPOS(JT,5))
           GDETAH => GDETAH_3D(KPOS(JT,4),:,KPOS(JT,5))
           IF (JVPOS(JT) /= KPOS(JT,3)) THEN
              CPOS(JT)= GETAH(JVPOS(JT))+ &
                   sc*HARVESTCONV(JRCNT(JT))*GDETAH(JVPOS(JT))
           END IF
           JRCNT(JT)= MERGE(JRCNT(JT)+1, 1, JRCNT(JT) < NGCELL)
           KPOS(JT,3)=GINDEX_V(CPOS(JT),GETAH) ! op_sb_2014021
           NPOS(3,JT)=REAL(KPOS(JT,3),dp)
           APOS(3,JT)=CPOS(JT)
        END DO

     END SELECT

   END SUBROUTINE ATTILA_CONVTRAJ
! -------------------------------------------------------------------

! -------------------------------------------------------------------
   SUBROUTINE ATTILA_SAVE_POS(RPOS, N)

     REAL(dp), DIMENSION(:), INTENT(OUT) :: RPOS
     INTEGER,                INTENT(IN)  :: N      ! INDEX IN NPOS

     RPOS(:) = NPOS(N,:)

   END SUBROUTINE ATTILA_SAVE_POS
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_GLOBAL_EXIT

    ! PURPOSE:
    !   DEALLOCATE MEMORY
    !
    ! AUTHOR(S):
    !   Patrick Joeckel, MPICH, May 2003

    IMPLICIT NONE

    DEALLOCATE(AARBOX)    ! AREA OF GRIDBOX

    DEALLOCATE(GLON)      ! LONGITUDE GRID (0..360-360/NLON)
    DEALLOCATE(GLAT)      ! GAUSSIAN LATITUDES (0..180) = N-POLE..S-POLE
    DEALLOCATE(GDLAT)     ! LAT. GRID SPACING, =GLAT(I+1)-GLAT(I)
    DEALLOCATE(GETAH_3D)  ! ETA VALUES ON HALF LEVELS
    DEALLOCATE(GETAF_3D)  ! ETA VALUES ON FULL LEVELS

    DEALLOCATE(GETA_3DF_L_P); NULLIFY(GETA_3DF_L_P)
    DEALLOCATE(GETA_3DH_L_P); NULLIFY(GETA_3DH_L_P)
    DEALLOCATE(GETA_3DF_P)  ; NULLIFY(GETA_3DF_P)
!!$    DEALLOCATE(GETA_3DH_P)  ; NULLIFY(GETA_3DH_P)
    NULLIFY(GETA_3DH_P)
    DEALLOCATE(GDETAH_3D) ! DELTA-ETA BETWEEN HALF LEVELS
    DEALLOCATE(GDETAF_3D) ! DELTA-ETA BETWEEN FULL LEVELS
!!$    IF (ASSOCIATED(TM1_3D)) DEALLOCATE(TM1_3D)    ! op_sb_20140211
    DEALLOCATE(GLON2)     ! LONGITUDE GRID (GRID POINTS = BOX BOUNDARIES)
    DEALLOCATE(GLAT2)     ! LATITUDE  GRID (GRID POINTS = BOX BOUNDARIES)

    DEALLOCATE(A3DRZ3)  ! =1./Z3
    DEALLOCATE(A3DDH1)  ! CONSTANTS NEEDED FOR CUBIC INTERPOLATION
    DEALLOCATE(A3DDH2)  ! CONSTANTS NEEDED FOR CUBIC INTERPOLATION

!!$    DEALLOCATE(PWU)
!!$    DEALLOCATE(PWV)
!!$    DEALLOCATE(PWW)
    DEALLOCATE(PMBOX)
    DEALLOCATE(PMBOXBL)
!!$    DEALLOCATE(PPH)
    DEALLOCATE(KHPBL)
    DEALLOCATE(PPS)

    !ATTILA_DRIVE
    DEALLOCATE(NCB)
    DEALLOCATE(NCBL)
    IF (ASSOCIATED(GNCB))        DEALLOCATE(GNCB)
    IF (ASSOCIATED(GNCBL))       DEALLOCATE(GNCBL)
    IF (ASSOCIATED(GNCB_MOVE))   DEALLOCATE(GNCB_MOVE)

    DEALLOCATE(vct)
    DEALLOCATE(gl_gmu)

    ! LOGICAL, SEE *ATTILA_ALLOC*
    IF (LALLOC) THEN
       DEALLOCATE(APOS)
       DEALLOCATE(NPOS)
       DEALLOCATE(NPOSM1)
       DEALLOCATE(APOSM1_PRESS)
    END IF

    ! ONLY FOR CONVECTION
    IF (ALLOCATED(PENTR))    DEALLOCATE(PENTR)
    IF (ALLOCATED(LPENTR))   DEALLOCATE(LPENTR)
    IF (ALLOCATED(LZMFU))    DEALLOCATE(LZMFU)
    IF (ALLOCATED(MAXCONV))  DEALLOCATE(MAXCONV)
    IF (ALLOCATED(PDETR))    DEALLOCATE(PDETR)

    IF (ALLOCATED(PPENDEU))  DEALLOCATE(PPENDEU)
    IF (ALLOCATED(PPENDED))  DEALLOCATE(PPENDED)
    IF (ALLOCATED(PPSUBS))   DEALLOCATE(PPSUBS)
    IF (ALLOCATED(KTYPELT))  DEALLOCATE(KTYPELT)
    IF (ALLOCATED(KCTOPLT))  DEALLOCATE(KCTOPLT)

  END SUBROUTINE ATTILA_GLOBAL_EXIT
! -------------------------------------------------------------------

! ==========================================================================
! === TRAJECTORY MODE
! ==========================================================================

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_RESETPOS_TRAJEC(ICELL,LAT,LON,PRESS)

    IMPLICIT NONE

    ! called from  attila_initialize_positions_t(flag) (SMIL)
    ! in attila_global_end
    ! flag=1: (tpot,temp,pps) not known!
    ! flag=2: (tpot,temp,pps) calculated

    ! I/O
    INTEGER,   INTENT(IN)  :: ICELL
    REAL(dp),  INTENT(IN)  :: LAT
    REAL(dp),  INTENT(IN)  :: LON
    REAL(dp),  INTENT(IN)  :: PRESS

    ! LOCAL
    INTEGER :: JL, JG ,JK
    REAL(dp):: pit2, deltap_loc
    ! op_sb_20160708+
    REAL(DP)              :: deltap(nlon,nlev,ngl)
    REAL(DP)              :: delta1(nlon,nlev,ngl)
    REAL(DP)              :: delta2(nlon,nlev,ngl)
    REAL(DP)              :: ptpot_lin(nlon,nlev,ngl)
    REAL(DP)              :: ptemp_h(nlon,nlevp1,ngl)
    REAL(DP)              :: ppf(nlon,nlev,ngl)

    pit2 = api / 2._dp
    ! op_sb_20160708-

    ! POSITION CELL
    APOS(1,ICELL)=LON
    APOS(2,ICELL)=LAT

    do jg=1,ngl
      do jk=1,nlev
       do jl=1,nlon
        deltap (jl,jk,jg)= (PPH(jl,jk+1,jg)-PPH(jl,jk,jg) ) ! between half levels
        delta1(jl,jk,jg) = (PPH(jl,jk+1,jg)-press)
        delta2(jl,jk,jg) = (press-PPH(jl,jk,jg))
       enddo
      enddo
    enddo

    JL = MOD( INT((APOS(1,ICELL)+GDHLON)/GDLON), NLON ) + 1
    JG = GINDEX(APOS(2,ICELL),GLAT2)

    ! SET SURFACE PRESSURE
    ! NOTE: use constant surface pressure for initialization
    !       if real surface pressure is not available
    !       (e.g., at very first time step)
    !       surface pressure is available if value is not zero
!$$    APOS(3,ICELL) = PRESS / PPS(JL,JG)
    APOS(4,ICELL) = PRESS

    SELECT CASE (I_VERT)
    CASE(1)
       APOS(3,ICELL) = PRESS / APZERO
       IF (APOS(3,ICELL) > 1._dp) THEN
          GETAF => GETAF_3D(JL,:,JG)
          APOS(3,ICELL)=GETAF(NLEV)
          WRITE(*,*) ' NEW INITIALIZATION OF TRAJECTORY ', ICELL,  &
               ' INITIALIZATION PRESSURE HEIGHT > SFC PRESSURE ', &
               ' SET TO LEVEL NLEV !!'
       ENDIF

    CASE(2)
       ! calculate temperatures at half levels and pressure at full levels
       do jg=1,ngl
          do jl=1,nlon
             ptemp_h(jl,1,jg)= ptemp(jl,1,jg)
             ptemp_h(jl,nlevp1,jg)= surftemp(jl,jg)
             ppf(jl,1,jg) = (pph(jl,1,jg)+pph(jl,2,jg)) * 0.5_dp
          enddo
       enddo

       do jg=1,ngl
          do jk=2,nlev
             do jl=1,nlon
                ppf(jl,jk,jg)     = (pph(jl,jk,jg) + pph(jl,jk+1,jg))*.5_dp
                ! between full levels
                deltap_loc        = (ppf(jl,jk,jg) - ppf(jl,jk-1,jg))
                ptemp_h(jl,jk,jg) = (ptemp(jl,jk,jg) * &
                     (pph(jl,jk,jg)-ppf(jl,jk-1,jg))+ &
                     ptemp(jl,jk-1,jg) * (ppf(jl,jk,jg)-pph(jl,jk,jg)))/ &
                     deltap_loc
             enddo
          enddo
       enddo

       JK=GINDEX(PRESS,PPH(JL,:,JG))
       JK = MIN(JK,NLEV)
       ! op_sb_20160708+
       ptpot_lin(jl,jk,jg) = (ptemp_h(jl,jk+1,jg) * delta2(jl,jk,jg)+          &
                        ptemp_h(jl,jk,jg) * delta1(jl,jk,jg))/deltap(jl,jk,jg) &
                      *(100000._dp/press)**kappa

       if (press .gt. (sigma_ref_g(jl,jg)*apzero)) then
         APOS(3,ICELL)=ptpot_lin(jl,jk,jg) * sin(pit2*(1._dp-press/pps(jl,jg))/&
                      (1._dp-sigma_ref_g(jl,jg)))
       else
         APOS(3,ICELL)=ptpot_lin(jl,jk,jg)
       endif
       ! op_sb_20160708-

       apos(3,icell) = MIN(apos(3,icell),theta_top-1.)
       IF (APOS(3,ICELL) <= 0._dp) THEN
          GETAF => GETAF_3D(JL,:,JG)
          WRITE(*,*) ' NEW INITIALIZATION OF TRAJECTORY ', ICELL,  &
               ' INITIALIZATION PRESSURE HEIGHT > SFC PRESSURE ', &
               ' SET TO LEVEL NLEV !!'
          APOS(3,ICELL)=GETAF(NLEV)
       END IF

    CASE(3)
       APOS(3,ICELL) = PRESS / PPS(JL,JG)
       IF (APOS(3,ICELL) > 1._dp) THEN
          GETAF => GETAF_3D(JL,:,JG)
          APOS(3,ICELL)=GETAF(NLEV)
          WRITE(*,*) ' NEW INITIALIZATION OF TRAJECTORY ', ICELL,  &
               ' INITIALIZATION PRESSURE HEIGHT > SFC PRESSURE ', &
               ' SET TO LEVEL NLEV !!'
       END IF
    END SELECT

    CALL ATTILA_CALIND(ICELL,1)

    ! INITIALIZE GLOBAL FIELDS
    NPOSM1(:,:) = NPOS(:,:)

  END SUBROUTINE ATTILA_RESETPOS_TRAJEC
! -------------------------------------------------------------------

! ###########################################################################
! ### PRIVATE SUBROUTINES
! ###########################################################################

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_INIT_GETA(status, l_init)

    ! PURPOSE:
    !   CALCULATE VERTICAL GRID
    !   (ETA FULL- AND HALF-LEVELS or sigma full- and half levels)
    !
    ! AUTHOR(S):
    !   Sabine Brinkop, DLR, Oct 2012
    !
    ! CALLED BY ATTILA_GLOBAL_INIT AND ATTILA_FILL
    !

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    LOGICAL, INTENT(IN)  :: l_init

    ! LOCAL
    INTEGER :: JL,JK,JG
    ! only for I_VERT == 2
    REAL(DP), DIMENSION(nlevp1) :: eta
    REAL(DP), DIMENSION(nlon,nlevp1,ngl) :: FUNC
!!    REAL(DP), DIMENSION(nlon,nlevp1,ngl) :: g_tpot_h ! op_sb_20160708
    REAL(DP), DIMENSION(nlon,nlevp1,ngl) :: pp
    ! op_sb_20160510+
    REAL(DP), DIMENSION(nlon,nlev,ngl)   :: ppf
    REAL(DP)              :: deltap(nlon,nlev,ngl)
    REAL(DP)              :: delta1(nlon,nlev,ngl)
    REAL(DP)              :: delta2(nlon,nlev,ngl)
    REAL(dp) :: pit2
    ! op_sb_20160510-

    pit2 = api/2._dp

    IF (l_init) THEN ! called by attila_global_init (attila_initialize)

       DO JG=1,NGL
          DO JK=1,NLEVP1
             DO JL=1,NLON
                GETAH_3D(JL,JK,JG)= VCT(JK)/PPS(JL,JG) + VCT(NVCLEV+JK)
             END DO
          END DO
       END DO

    ELSE ! time loop (attila_global_end) and lresume=true

       SELECT CASE(I_VERT)

       CASE(1)
          !
          DO JG=1,NGL
             DO JK=1,NLEVP1
                DO JL=1,NLON
                   GETAH_3D(JL,JK,JG)= VCT(JK)/APZERO + VCT(NVCLEV+JK)
                END DO
             END DO
          END DO

       CASE(2)
          !
          IF (.NOT. ASSOCIATED(PTEMP)) THEN
             status = 1
             RETURN
          ENDIF
          IF (.NOT. ASSOCIATED(PTPOT)) THEN
             status = 2
             RETURN
          ENDIF

          do jk = 1, nlevp1
             eta(jk) = VCT(JK)/apzero + VCT(NVCLEV+JK)
             pp(:,jk,:)  = VCT(JK) + VCT(NVCLEV+JK) * pps(:,:)
          enddo

          ! upper and lower boundary
          g_tpot_h(:,1,:) = theta_top
          g_tpot_h(:,nlevp1,:) = ptemp(:,nlev,:)*(100000._dp/pps(:,:))**kappa
          ppf(:,1,:)           = (pp(:,1,:) +  pp(:,2,:)) * 0.5
          ppf(:,nlev,:)        = (pp(:,nlev,:) +  pp(:,nlevp1,:)) * 0.5

          ! pressure at full levels
          do jk= 2,nlev-1
           ppf(:,jk,:)        = (pp(:,jk,:) +  pp(:,jk+1,:)) * 0.5
          enddo

          ! potential temperature at half levels
          do jk = 2,nlev
!!             g_tpot_h(:,jk,:) = (PTPOT(:,jk-1,:) + PTPOT(:,jk,:))*0.5_dp
               deltap(:,jk,:) = (ppf(:,jk,:) - ppf(:,jk-1,:))
               delta1(:,jk,:) = (PP(:,jk,:)-ppf(:,jk-1,:))
               delta2(:,jk,:) = (ppf(:,jk,:)-PP(:,jk,:))

               g_tpot_h(:,jk,:) = (PTPOT(:,jk-1,:) * delta2(:,jk,:)    &
                   + PTPOT(:,jk,:) * delta1(:,jk,:))/deltap(:,jk,:)
          enddo

          DO JG=1,NGL
          DO JK=1,NLEVP1
           DO JL=1,NLON
!             if (eta(jk) >= eta_ref) then
             if ( pp(jl,jk,jg) > sigma_ref_g(jl,jg)*apzero) then
                ! troposphere
                FUNC(jl,jk,jg)=sin(pit2*(1._dp-pp(jl,jk,jg)/pps(jl,jg))/(1._dp-sigma_ref_g(jl,jg)))
!                FUNC(:,JK,:)=((1._dp-pp(:,jk,:)/pps(:,:))/     &
!                         (1._dp-sigma_ref_g(:,:)))**(kappa) !op_sb_20160502
             else
                ! stratosphere
                FUNC(jl,jk,jg) = 1._dp
             endif
          END DO
          enddo
          enddo

          ! new grid
          DO jg= 1,NGL
             DO JK=1,NLEVP1
                DO JL=1,NLON
                     GETAH_3D(JL,JK,JG)= G_TPOT_H(JL,JK,JG) * FUNC(JL,JK,JG)
                END DO
             END DO
          END DO

       CASE(3)
          !
          IF (.NOT. ASSOCIATED(PPS)) THEN
             status = 3
             RETURN
          ENDIF

          DO JG=1,NGL
             DO JK=1,NLEVP1
                DO JL=1,NLON
                   GETAH_3D(JL,JK,JG)= VCT(JK)/PPS(JL,JG) + VCT(NVCLEV+JK)
                END DO
             END DO
          END DO

       END SELECT
    ENDIF

    DO JG=1,NGL
       DO JK=1,NLEV
          DO JL=1,NLON
             GETAF_3D(JL,JK,JG)  = &
                  (GETAH_3D(JL,JK,JG) + GETAH_3D(JL,JK+1,JG)) * 0.5_dp
             GDETAH_3D(JL,JK,JG) = ABS(GETAH_3D(JL,JK+1,JG) &
                  - GETAH_3D(JL,JK,JG))
          END DO
       END DO

       DO JK=1,NLEV-1
          DO JL=1,NLON
             GDETAF_3D(JL,JK,JG) = GETAF_3D(JL,JK+1,JG) - GETAF_3D(JL,JK,JG)
          END DO
       END DO
    END DO

!       PRINT*,'GETAH: '
!       WRITE(*,'(1X,9G10.3)') (GETAH_3D(1,JK,1),JK=1,NLEVP1)
!       PRINT*,'GETAF: '
!       WRITE(*,'(1X,9G10.3)') (GETAF_3D(1,JK,1),JK=1,NLEV)

    status = 0

  END SUBROUTINE ATTILA_INIT_GETA
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_INICOMPH

    ! PURPOSE:
    !   - INITIALISE CONSTANTS FOR PHYSICS
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTINICOMPH)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_INICOMPH)
    !
    ! INTERFACE:
    !   - LTINICOMPH WAS CALLED FROM LTINICOM
    !   - ATTILA_INICOMPH IS CALLED FROM ATTILA_INICOM_2

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: j

    LLTMOCA=.FALSE.

    DO J=1,3
       IF (ADICO(J) > 0._dp) THEN
          ADISP(J)= SQRT(2*DTIME*ADICO(J))
          LLTMOCA=.TRUE.
       ELSE
          ADICO(J)= 0._dp
          ADISP(J)= 0._dp
       END IF
    END DO
    !
    !    CONVERT HORIZONTAL DISPLACEMENT FROM M TO DEG.
    !
    DO J=1,2
       ADISP(J)= ADISP(J) * ADEGPM
    END DO

    ABLPLUS = ADISP(3) / 4._dp

    PRINT*,'=================== ATTILA_INICOMPH ========================'
    PRINT*,' ADICO: ', ADICO
    PRINT*,' ADISP: ', ADISP
    PRINT*,' ABLPLUS: ',ABLPLUS,', LLTMOCA: ',LLTMOCA

  END SUBROUTINE ATTILA_INICOMPH
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_PRESH(PPS,PP)

    ! PURPOSE:
    !  - CALCULATE PRESSURE ON HALF LEVELS
    !    DEPENDING ON SURFACE PRESSURE (PPS) CALCULATE PRESSURE
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTPRESH)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_PRESH)
    !
    ! INTERFACE:
    !   - LTPRESH WAS CALLED FROM LTINIPOS AND LTTHEOCELL
    !   - ATTILA_PRESH IS CALLED FROM ATTILA_INIPOS AND ATTILA_THEOCELL
    !
    ! PARAMETERS:
    !   I PPS  SURFACE PRESSURE
    !   O PP   PRESSURE ON HALF LEVELS
    !
    ! METHODS:
    !   CALCULATIONS ARE PERFORMED SEPARATELY FOR PRESSURE, HYBRID
    !   AND SIGMA LEVELS.

    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(IN)  :: PPS
    REAL(dp), INTENT(OUT) :: PP(NLEVP1)

    ! LOCAL
    INTEGER J

    !
    !*    1. CALCULATE PRESSURE LEVEL VALUES.
    !        --------- -------- ----- -------
    !
    DO J=1,NPLVP1
       PP(J)= VCT(J)
    ENDDO

    !
    !*    2. CALCULATE HYBRID LEVEL VALUES.
    !        --------- ------ ----- -------
    !
    DO J= NPLVP2,NLMSGL
       PP(J)= VCT(J) + VCT(J+NVCLEV)*PPS
    ENDDO
    !
    !*    3. CALCULATE SIGMA LEVEL VALUES.
    !        --------- ----- ----- -------
    !
    DO J= NLMSLP, NLEVP1
       PP(J)= VCT(J+NVCLEV)*PPS
    ENDDO

  END SUBROUTINE ATTILA_PRESH
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_THEOCELL

    ! PURPOSE:
    !   - PRINT THEORETICAL NUMBER OF CELLS IN EACH LATITUDE BAND OR LEVEL
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTTHEOCELL)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_THEOCELL)
    !
    ! INTERFACE:
    !   - LTTHEOCELL WAS CALLED FROM LTINIPOS
    !   - ATTILA_THEOCELL IS CALLED FROM ATTILA_INIPOS
    !
    ! METHODS:
    !       CELLS ARE SUPPOSED TO HAVE EQUAL MASS, THERFORE:
    !       NUMBER OF CELLS IN LATITUDE BAND IS PROPORTIONAL TO
    !       INTEGRAL OF COS(LAT)
    !       NUMBER OF CELLS IN LEVEL IS PROP. TO PRESSURE DIFFERENCE

     IMPLICIT NONE

     ! LOCAL
     REAL(dp)                           :: ZF,ZRAD
     INTEGER, DIMENSION(:), ALLOCATABLE :: ILEV, ILAT
     REAL(dp)                           :: ZP(NLEVP1)
     INTEGER                            :: J

     ! INIT
     ALLOCATE(ILEV(NLEV), ILAT(NGL))

     ! SET PARAMETER
     ZRAD= API/180._dp

     CALL ATTILA_PRESH(APSURF,ZP)

     DO J=1,NLEV
        ZF= (ZP(J+1)-ZP(J))/(ZP(NLEVP1)-ZP(1))
        ILEV(J)= INT(ZF*REAL(NCELL,dp))
     ENDDO
     DO J=1,NGL
        ZF= (SIN((GLAT2(J+1)-90._dp)*ZRAD) - SIN((GLAT2(J)-90._dp)*ZRAD)) *.5_dp
        ILAT(J)= INT(ZF*REAL(NCELL,dp))
     ENDDO

     PRINT*,' --- THEORETICAL NUMBER OF CELLS ---'
     PRINT*,'GLOBAL MEAN OF SURFACE PRESSURE: ',ZP(NLEVP1) !APSURF
     PRINT*,' NUMBER IN LEV: '
     WRITE(*,'(1X,10I7)') (ILEV(J),J=1,NLEV)
     PRINT*,' NUMBER IN LAT: '
     WRITE(*,'(1X,10I7)') (ILAT(J),J=1,NGL)

     ! CLEAN UP
     DEALLOCATE(ILEV, ILAT)

   END SUBROUTINE ATTILA_THEOCELL
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_CELLPOS(NRP,HARVESTINIPOS,KL,KG,KK,KCELL)

    ! PURPOSE:
    !   - POSITION A CELL (HORIZONTAL POSITION AND ETA COORDINATE)
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTCELLPOS)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_CELLPOS)
    !      - CALL TO RANDNUM REMOVED
    !
    ! INTERFACE:
    !   - LTCELLPOS WAS CALLED FROM LTINIPOS
    !   - ATTILA_CELLPOS IS CALLED FROM ATTILA_INIPOS
    !
    ! PARAMETERS:
    !   I NRP               RANDOM NUMBER TRIPLE POSITION
    !   I HARVESTINIPOS     RANOM NUMBER ARRAY (NCELL, NR1=4)
    !   I KL                LONGITUDE INDEX (1..NLON)
    !   I KG                LATITUDE INDEX  (1..NGL)
    !   I KK                ALTITUDE INDEX  (1..NLEV)
    !   I KCELL             NUMBER OF CELL TO BE POSITIONED
    !
    ! REMARKS:
    !   PPOS(1-3) : HORIZ. &  ETA POSITION
    !   KPOS(1-5) : INTEGER POSITION
    !
    ! EXTERNALS:
    !   - ATTILA_CALIND CALCULATE INDICES (INTEGER POSITION)
    !
    ! METHODS:
    !   - CELLS ARE DISTRIBUTED RANDOMLY WITHIN BOX

    IMPLICIT NONE

    ! I/O
    INTEGER,INTENT(INOUT)                :: NRP
    INTEGER,                  INTENT(IN) :: KL,KG,KK,KCELL
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: HARVESTINIPOS

    ! LOCAL
    REAL(dp) :: ZL,ZG,ZK
    INTEGER  :: IL
    LOGICAL  :: l_ptemp

    l_ptemp = ASSOCIATED(ptemp)

    APOS(1,KCELL)=0._dp

    ZL= HARVESTINIPOS(NRP,1)

    ZG= HARVESTINIPOS(NRP,2)

    ZK= HARVESTINIPOS(NRP,3)

    NRP = NRP + 1 ! 1 random number triple used

    IL= MOD( KL+NLON-2 , NLON ) + 1

    GETAH => GETAH_3D(KL,:,KG)

    !CALCULATING THE POSITIONS
    APOS(1,KCELL)= MOD( GLON2(IL) + ZL*GDLON , 360._dp )
    APOS(2,KCELL)= GLAT2(KG) + ZG* (GLAT2(KG+1)-GLAT2(KG))
    APOS(3,KCELL)= GETAH(KK) + ZK* (GETAH(KK+1)-GETAH(KK))

    !PROVING THE POSITIONS
    IF(APOS(1,KCELL) >= 360._dp) THEN
       PRINT*,' POS(1) OF CELL ',KCELL,' >=360: ', &
            APOS(1,KCELL),', ADJUSTED.'
       APOS(1,KCELL)= APOS(1,KCELL) - 360._dp
    ENDIF
    IF(APOS(2,KCELL) > 179.9_dp) THEN
       PRINT*,' POS(2) OF CELL ',KCELL,' >179.9 : ',APOS(2,KCELL) &
            ,', SET TO 179.9.'
       APOS(2,KCELL)= 179.9_dp
    ELSEIF(APOS(2,KCELL) < 0.1_dp) THEN
       PRINT*,' POS(2) OF CELL ',KCELL,' <0.1 : ',APOS(2,KCELL) &
            ,', SET TO 0.1.'
       APOS(2,KCELL)= 0.1_dp
    ENDIF

    SELECT CASE (I_VERT)

    CASE(1, 3)
       IF(APOS(3,KCELL) >= GETAH(NLEVP1)) THEN
          PRINT*,' POS(3) OF CELL ',KCELL,' >= ETA_SURF : ', &
               APOS(3,KCELL), GETAH(NLEVP1), ', SET NEAR SURFACE.'
          APOS(3,KCELL)= 0.99_dp*GETAH(NLEVP1) + (1._dp-0.99_dp)*GETAH(NLEV)
       ELSEIF(APOS(3,KCELL) <= GETAH(1)) THEN
          PRINT*,' POS(3) OF CELL ',KCELL,' <= TOP : ', &
               APOS(3,KCELL), GETAH(1), ', SET NEAR TOP.'
          APOS(3,KCELL)= 0.99_dp*GETAH(1) + (1._dp-0.99_dp)*GETAH(2)
       ENDIF

    CASE(2)

       IF (.NOT. l_ptemp) THEN

          IF(APOS(3,KCELL) >= GETAH(NLEVP1)) THEN
             PRINT*,' POS(3) OF CELL ',KCELL,' >= ETA_SURF (theta_ini) : ', &
                  APOS(3,KCELL), GETAH(NLEVP1), ', SET NEAR SURFACE.'
             APOS(3,KCELL)= 0.99_dp*GETAH(NLEVP1) + (1._dp-0.99_dp)*GETAH(NLEV)
          ELSEIF(APOS(3,KCELL) <= GETAH(1)) THEN
             PRINT*,' POS(3) OF CELL ',KCELL,' <= TOP : ', &
                  APOS(3,KCELL), GETAH(1), ', SET NEAR TOP.'
             APOS(3,KCELL)= 0.99_dp*GETAH(1) + (1._dp-0.99_dp)*GETAH(2)
          ENDIF
       ELSE
          IF(APOS(3,KCELL) <= GETAH(NLEVP1)) THEN
             PRINT*,' POS(3) OF CELL ',KCELL,' <= THETA_SURF : ', &
                  APOS(3,KCELL), GETAH(NLEVP1), ', SET NEAR SURFACE.'
             APOS(3,KCELL)= 0.99_dp*GETAH(NLEVP1) + (1._dp-0.99_dp)*GETAH(NLEV)
          ELSEIF(APOS(3,KCELL) >= GETAH(1)) THEN
             PRINT*,' POS(3) OF CELL ',KCELL,' >= TOP : ', &
                  APOS(3,KCELL), GETAH(1), ', SET NEAR TOP.'
             APOS(3,KCELL)= 0.99_dp*GETAH(1) + (1._dp-0.99_dp)*GETAH(2)
          ENDIF
       ENDIF

    END SELECT

    !
    !*   CALCULATE INDICES.
    !
    CALL ATTILA_CALIND(KCELL,1)

    ! SET APOS(4) - PRESSURE HEIGHT OF EACH CELL - EQUAL ZERO AT
    ! INITIALIZATION
    APOS(4,:)=0.0_dp

  END SUBROUTINE ATTILA_CELLPOS
! -------------------------------------------------------------------

! -------------------------------------------------------------------
   SUBROUTINE ATTILA_CAT_MOVE(PPOS,KPOS,HARVESTCAT,TI2)

     IMPLICIT NONE

     ! PURPOSE:
     !   - TURBULENT MOVEMENT IN VERTICAL DIRECTION IN CLEAR AIR TURBULENCE
     !     REGIONS
     !
     ! AUTHOR(S):
     !   M. Traub, MPI Mainz
     !
     ! INTERFACE:
     !   - ATTILA_CAT_MOVE IS CALLED FROM ATTILA_DRIVE
     !
     ! PARAMETERS:
     !   IO PPOS
     !   IO KPOS
     !   IO HARVESTCAT   RANDOM NUMBER ARRAY
     !   IO TI2          TURBULENCE INDEX

     ! IN/OUT
     REAL(dp)     :: PPOS(NCEVL)        ! VERTICAL POSITION
     INTEGER      :: KPOS(NCEVL,NNPOS)
     REAL(dp)     :: HARVESTCAT(:)
     REAL(dp)     :: TI2(:,:,:)

     ! LOCAL
     INTEGER  :: JK, JL, JG, JT
     REAL(dp) :: deta1, deta2

     DO JT= 1,NCEVL    ! LOOP OVER ALL CELLS
        JL= KPOS(JT,4)
        JG= KPOS(JT,5)
        JK= KPOS(JT,3) ! ACTUAL LEVEL INDEX, FOR TURBULENT MOVEMENT

        GETAH => GETAH_3D(JL,:,JG)

        ! LEVEL INDEX MUST BE MINIMUM 1! OTHERWISE BUG IN PROGRAM
        ! MAXIMUM LEVEL INDEX JK = (NLEV-1), BECAUSE FOR CALCULATION
        ! OF DETA1, GETAH(JK+2) --> GETAH(NLEVP1) IS NEEDED
        ! CHECK THE TI2 INDEX --> WHEN STRONG TURBULENCE??

        IF ((JK > 1) .AND. (TI2(JL,JG,JK) >= 3._dp).AND. (JK <= (NLEV-1))) THEN

           ! DELTA TO LOWER BOUNDARY
           deta1=ABS(APOS(3,JT)-GETAH(JK+2))

           ! DELTA TO UPPER BOUNDARY
           deta2=ABS(APOS(3,JT)-GETAH(JK))

           IF ((2*HARVESTCAT(JT)-1._dp) < 0) THEN

              ! DOWNWARD MOVEMENT
              PPOS(JT)= PPOS(JT)-ABS((2*HARVESTCAT(JT)-1._dp))*deta1

           ELSE

              ! UPWARD MOVEMENT
              PPOS(JT)= PPOS(JT)+ABS((2*HARVESTCAT(JT)-1._dp))*deta2

           ENDIF

        ENDIF

     END DO

     ! CALCULATE NEW VERTICAL INTEGER POSITION
     CALL ATTILA_GINDEXV(PPOS(:),GETAH,KPOS(:,3))

   END SUBROUTINE ATTILA_CAT_MOVE
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_RKO4(PPOS,KPOS,PRAN,PVEL)

    ! PURPOSE:
    !   - ADVECT CELL IN WIND FIELD
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTRKO4)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_RKO4)
    !      - WULT -> PWU , WVLT -> PWV , WELT -> PWW
    !
    ! INTERFACE:
    !   - LTRKO4 WAS CALLED FROM LTDRIVE
    !   - ATTILA_RKO4 IS CALLED FROM ATTILA_DRIVE
    !
    ! PARAMETERS:
    !   IO PPOS    (NCEVL,NAPOS):         POSITION
    !   IO KPOS    (NCEVL,NNPOS): INTEGER POSITION
    !                             NB: KPOS(4-5) IS NOT NEEDED
    !   I  PRAN    (NCEVL,4)    : RANDOM NUMBERS:
    !                         PRAN(1-3) (-1 .. +1) FOR MONTE CARLO DIFFUSION
    !                         PRAN(4)   ( 0 ..  1) FOR BL TURBULENCE
    !
    ! REMARKS:
    !       INPUT:
    !            PPOS(1-3) = POSITION
    !            KPOS(1-3) = INTEGER POSITION
    !       OUTPUT:
    !            PPOS(1-3) = NEW POSITION
    !            KPOS(1-5) = NEW INTEGER POSITION
    !
    ! METHODS:
    !      APPLY *RUNGE - *KUTTA METHOD (4TH ORDER ACCURACY) TO WIND FIELD
    !      (AT TIME) IN ORDER TO CALCULATE NEW POSITION (AT TIME+DTIME).
    !      ASSUMPTION: WINDS ARE CONSTANT IN TIME INTERVALL
    !                  [TIME;TIME+DTIME]
    !      SMALL SCALE TURBULENCE BY RANDOM DISPLACEMENT.
    !      MIXING IN BL BY RANDOM REASSIGNMENT OF VERTICAL POSITION.
    !
    ! EXTERNALS:
    !      ATTILA_LOGGRID  CALCULATE LOG(GRID) FOR INTERPOLATION
    !      ATTILA_VELO     INTERPOLATE WIND FIELD
    !      ATTILA_REPOS    REPOSITION CELL
    !      ATTILA_CALINDV  CALCULATE INDICES KPOS(1-5)  (NO=1)
    !      ATTILA_CALINDV  CALCULATE INDICES KPOS(1-3)  (NO=2)
    !      ATTILA_CALINDV  CALCULATE INDICES KPOS(3-5)  (NO=3)
    !
    ! REFERENCE:
    !      RUNGE-KUTTA METHOD:
    !        BOOK: *W. *H. *PRESS ET AL., *NUMERICAL *RECIPES, P 550.

    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(INOUT)  :: PPOS(NCEVL,NAPOS)   ! NCEVL x NAPOS
    INTEGER , INTENT(INOUT)  :: KPOS(NCEVL,NNPOS)   ! NCEVL x NNPOS
    REAL(dp), INTENT(IN) :: PRAN(NCEVL,4)       ! NCEVL x 4
    REAL(dp), OPTIONAL, INTENT(OUT) :: PVEL(NCEVL,3) ! NCEVL x 3

    ! LOCAL
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: ZPOS
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZV
    INTEGER,  DIMENSION(:,:),   ALLOCATABLE :: IPOS
    INTEGER                                 :: JL, JG, JHPBL, JCEVL
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: ZDISP
    LOGICAL                                 :: LOINPBL, LOMOCA, LOBLTURB
    LOGICAL,  DIMENSION(:),     ALLOCATABLE :: LOREP

    ! INIT
    ALLOCATE(ZPOS(NCEVL,NAPOS), ZV(NCEVL,3,5))
    ALLOCATE(IPOS(NCEVL,NNPOS))
    ALLOCATE(ZDISP(NCEVL,3))
    ALLOCATE(LOREP(NCEVL))

    !
    !*    1. FIND INTEGER POSITION ACCORDING TO FULL LEVEL ETA GRID.
    !        ---- ------- -------- --------- -- ---- ----- --- -----
    !
    !CC       --- OBSOLETE, DONE IN ATTILA_VELO NOW ---
    !
    !*    2. FIND WINDS AT TIME, USE TO FIND NEW POSITION AT TIME+DTIME/2.
    !        ---- ----- -- ----- --- -- ---- --- -------- -- -------------
    !
    CALL ATTILA_LOGGRID
    ! op_sb_20161017+
    CALL ATTILA_GINDEXV_Z(PPOS(:,3),KPOS(:,4:5),GETAH_3D,KPOS(:,3))
    ! op_sb_20161017-
    CALL ATTILA_VELO(ZV(:,:,1),PPOS,KPOS)  ! GET WIND

    IF (L_LG_FMM) THEN
       DO JCEVL=1,NCEVL
          APOS_FMM(JCEVL, :, 1) = PPOS(JCEVL,:)
       ENDDO
    ENDIF

    DO JCEVL=1,NCEVL  ! NEW POSITION
       LOREP(JCEVL)=  .TRUE.
       ZPOS(JCEVL,1)= PPOS(JCEVL,1) + ZV(JCEVL,1,1)*DTIMED2
       ZPOS(JCEVL,2)= PPOS(JCEVL,2) + ZV(JCEVL,2,1)*DTIMED2
       ZPOS(JCEVL,3)= PPOS(JCEVL,3) + ZV(JCEVL,3,1)*DTIMED2
    ENDDO

    CALL ATTILA_REPOS(ZPOS,LOREP)      ! REPOSITION
    CALL ATTILA_CALINDV(ZPOS,IPOS,1)   ! CALC IPOS(1-5)

    IF (L_LG_FMM) THEN
       DO JCEVL=1,NCEVL
          APOS_FMM(JCEVL, :, 2) = ZPOS(JCEVL,:)
       ENDDO
    ENDIF

    !
    !*    3. ESTIMATE WIND AT TIME+DTIME/2 USING POSITION FROM 2., USE TO
    !        -------- ---- -- ------------ ----- -------- -------- --- --
    !        FIND NEW POSITION AT TIME+DTIME/2.
    !        ---- --- -------- -- -------------
    !
    CALL ATTILA_VELO(ZV(:,:,2),ZPOS,IPOS)  ! GET WIND

    DO JCEVL=1,NCEVL  ! NEW POSITION
       ZPOS(JCEVL,1)= PPOS(JCEVL,1) + ZV(JCEVL,1,2)*DTIMED2
       ZPOS(JCEVL,2)= PPOS(JCEVL,2) + ZV(JCEVL,2,2)*DTIMED2
       ZPOS(JCEVL,3)= PPOS(JCEVL,3) + ZV(JCEVL,3,2)*DTIMED2
    ENDDO

    CALL ATTILA_REPOS(ZPOS,LOREP)      ! REPOSITION
    CALL ATTILA_CALINDV(ZPOS,IPOS,2)   ! CALC IPOS(1-5)

    IF (L_LG_FMM) THEN
       DO JCEVL=1,NCEVL
          APOS_FMM(JCEVL, :, 3) = ZPOS(JCEVL,:)
       ENDDO
    ENDIF

    !
    !*    4. ESTIMATE WIND AT TIME+DTIME/2 USING POSITION FROM 3., USE TO
    !        -------- ---- -- ------------ ----- -------- -------- --- --
    !        FIND NEW POSITION AT TIME+DTIME.
    !        ---- --- -------- -- -----------
    !
    CALL ATTILA_VELO(ZV(:,:,3),ZPOS,IPOS)  ! GET WIND

    DO JCEVL=1,NCEVL  ! NEW POSITION
       ZPOS(JCEVL,1)= PPOS(JCEVL,1) + ZV(JCEVL,1,3)*DTIME
       ZPOS(JCEVL,2)= PPOS(JCEVL,2) + ZV(JCEVL,2,3)*DTIME
       ZPOS(JCEVL,3)= PPOS(JCEVL,3) + ZV(JCEVL,3,3)*DTIME
    ENDDO
    CALL ATTILA_REPOS(ZPOS,LOREP)      ! REPOSITION
    CALL ATTILA_CALINDV(ZPOS,IPOS,1)   ! CALC IPOS(1-5)

    IF (L_LG_FMM) THEN
       DO JCEVL=1,NCEVL
          APOS_FMM(JCEVL, :, 4) = ZPOS(JCEVL,:)
       ENDDO
    ENDIF

    !
    !*    5. ESTIMATE WIND AT TIME+DTIME USING POSITION FROM 4.
    !        -------- ---- -- ---------- ----- -------- ---- --
    !
    CALL ATTILA_VELO(ZV(:,:,4),ZPOS,IPOS)  ! GET WIND

    !*    6. ESTIMATE FINAL POSITION.
    !        -------- ----- ---------
    !
    DO JCEVL=1,NCEVL  ! NEW POSITION
       ZV(JCEVL,1,5)=    ZV(JCEVL,1,1) + ZV(JCEVL,1,4) +       &
            2*(ZV(JCEVL,1,2) + ZV(JCEVL,1,3))
       ZV(JCEVL,2,5)=    ZV(JCEVL,2,1) + ZV(JCEVL,2,4) +       &
            2*(ZV(JCEVL,2,2) + ZV(JCEVL,2,3))
       ZV(JCEVL,3,5)=    ZV(JCEVL,3,1) + ZV(JCEVL,3,4) +       &
            2*(ZV(JCEVL,3,2) + ZV(JCEVL,3,3))
       PPOS(JCEVL,1)= PPOS(JCEVL,1) + ZV(JCEVL,1,5) * DTIMED6
       PPOS(JCEVL,2)= PPOS(JCEVL,2) + ZV(JCEVL,2,5) * DTIMED6
       PPOS(JCEVL,3)= PPOS(JCEVL,3) + ZV(JCEVL,3,5) * DTIMED6
       IF (LVDIAG) THEN
          ! for conversion see ATTILA_VELO
          PVEL(JCEVL,1) = ( ZV(JCEVL,1,5) / 6.0_dp ) / ADEGPM * &
               COS( (PPOS(JCEVL,2)-90._dp)*ARAD )                     ! [m/s]
          PVEL(JCEVL,2) = ( - ZV(JCEVL,2,5) / 6.0_dp ) / ADEGPM       ! [m/s]
          PVEL(JCEVL,3) = ( ZV(JCEVL,3,5) / 6.0_dp )                  ! [1/s] or. [K/s]
       END IF
    ENDDO
    CALL ATTILA_REPOS(PPOS,LOREP)    ! REPOSITION
    !
    !     7. DIFFUSION OUTSIDE OF POLAR REGION AND RANDOM REASSIGNMENT
    !        --------- ------- -- ----- ------ --- ------ ------------
    !        IN PBL LAYER.
    !        -- --- ------
    !
    IF(LLTMOCA.OR.LLTBLTURB) THEN

       CALL ATTILA_CALINDV(PPOS,KPOS,3)    ! CALC KPOS(3-5)

       DO JCEVL= 1,NCEVL

          JL=KPOS(JCEVL,4)
          JG=KPOS(JCEVL,5)
          GETAH => GETAH_3D(JL,:,JG)
          JHPBL=KHPBL(JL,JG) + 1
          LOINPBL= (KPOS(JCEVL,3) >= JHPBL)   ! IN PBL
          LOMOCA = (LLTMOCA.AND.(JG > 1).AND.(JG < NGL)) !DO MONTE-CARLO DIFF
          LOBLTURB= LOINPBL .AND. LLTBLTURB
          LOREP(JCEVL)= LOMOCA .OR. LOBLTURB
          ZDISP(JCEVL,1)= MERGE(ADISP(2),ADISP(1),LOINPBL)
          ZDISP(JCEVL,2)= ZDISP(JCEVL,1)
          ZDISP(JCEVL,3)= MERGE(ADISP(3)*PRAN(JCEVL,3),0._dp,LOMOCA)
          PPOS(JCEVL,3) = MERGE(                 &
               1._dp - PRAN(JCEVL,4) *              &
               ( 1._dp - GETAH(JHPBL) + ABLPLUS )   &
               , PPOS(JCEVL,3) + ZDISP(JCEVL,3)  &
               , LOBLTURB )

          IF(LOMOCA) THEN
             PPOS(JCEVL,1)= PPOS(JCEVL,1) + PRAN(JCEVL,1) * ZDISP(JCEVL,1)   &
                  / COS( (PPOS(JCEVL,2)-90._dp)*ARAD )
             PPOS(JCEVL,2)= PPOS(JCEVL,2) + PRAN(JCEVL,2) * ZDISP(JCEVL,2)
          ENDIF

       ENDDO

       CALL ATTILA_REPOS(PPOS,LOREP)

    ENDIF

    CALL ATTILA_CALINDV(PPOS,KPOS,1)     ! CALC KPOS(1-5)

    ! CLEAN UP
    DEALLOCATE(ZPOS, ZV)
    DEALLOCATE(IPOS)
    DEALLOCATE(ZDISP)
    DEALLOCATE(LOREP)

  END SUBROUTINE ATTILA_RKO4
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_CALIND(KCELL, NO)

    ! PURPOSE:
    !   - CALCULATION OF INDICES/INTEGER POSITIONS OF CELLS
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTCALIND[,2,3])
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_CALIND)
    !
    ! INTERFACE:
    !   - LTCALIND[,2,3] WERE CALLED FROM LTCELLPOS
    !   - ATTILA_CALIND IS CALLED FROM ATTILA_CELLPOS
    !
    ! EXTERNALS:
    !   - GINDEX    GET INTEGER POSITION IN AN ASCENDING GRID.
    !
    ! PARAMETERS:
    !   - KCELL  NUMBER OF CELL
    !   - NO     SWITCH (1,2,3)
    !
    ! REAMRKS:
    !   KPOS(1-5) INTEGER POSITION (IF NO=1)
    !   KPOS(1-3) INTEGER POSITION (IF NO=2)
    !   KPOS(3-5) INTEGER POSITION (IF NO=3)
    !

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)  :: NO       ! 1,2,3 (REPLACE ENTRY)
    INTEGER, INTENT(IN)  :: KCELL

    IF (NO == 1) THEN
       NPOS(4,KCELL)= REAL(MOD( INT((APOS(1,KCELL)+GDHLON)/GDLON), NLON ) &
            + 1, dp)
       NPOS(5,KCELL)= REAL(GINDEX(APOS(2,KCELL),GLAT2), dp)
    END IF

    !------------------

    IF ((NO == 1).OR.(NO == 2)) THEN
       GETAH => GETAH_3D(INT(NPOS(4,KCELL)),:,INT(NPOS(5,KCELL)))
       NPOS(1,KCELL)= REAL(INT( APOS(1,KCELL)/GDLON ) + 1, dp)
       NPOS(2,KCELL)= REAL(GINDEX(APOS(2,KCELL), GLAT), dp)
       NPOS(3,KCELL)= REAL(GINDEX_V(APOS(3,KCELL), GETAH), dp) ! op_sb_20140515
       RETURN
    END IF

    !------------------
    GETAH => GETAH_3D(INT(NPOS(4,KCELL)),:,INT(NPOS(5,KCELL)))
    NPOS(3,KCELL)= REAL(GINDEX_V(APOS(3,KCELL), GETAH), dp)   ! op_sb_20140515
    NPOS(4,KCELL)= REAL( MOD( INT((APOS(1,KCELL)+GDHLON)/GDLON) , NLON ) &
         + 1, dp)
    NPOS(5,KCELL)= REAL(GINDEX(APOS(2,KCELL),GLAT2), dp)

    !------------------

  END SUBROUTINE ATTILA_CALIND
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_LOGGRID

    IMPLICIT NONE

    ! local
    INTEGER :: jk

    !   Fields are necessary for interpolation of wind field on grid, where
    !   the parcel exist

    GETA_3DH_P   => pph(:,:,:)

    do jk = 1, nlev
       GETA_3DF_P(:,jk,:)= (pph(:,jk+1,:) + pph(:,jk,:)) * 0.5
    enddo

    GETA_3DF_L_P(:,:,:) = log(GETA_3DF_P(:,:,:))
    GETA_3DH_L_P(:,:,:) = log(GETA_3DH_P(:,:,:))

  END SUBROUTINE ATTILA_LOGGRID
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_VELO(PVELO,PPOS,KPOS)

    ! PURPOSE:
    !   - RETURN WIND VELOCITIES FOR ADVECTION
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTVELO)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_VELO)
    !
    ! INTERFACE:
    !   - LTVELO WAS CALLED FROM LTRKO4
    !   - ATTILA_VELO IS CALLED FROM ATTILA_RKO4
    !
    ! PARAMETERS:
    !   O PVELO  (NCEVL,3): RETURNED VELOCITY, IN DEGREE/SEC RESP. 1/SEC
    !                       NB: POSITIVE *PVELO(2)* MEANS NORTH-WIND
    !   I PPOS   (NCEVL,NAPOS)      : POSITION
    !   I KPOS   KPOS*(NCEVL,NNPOS) : INTEGER POSITION
    !                             NB: KPOS(3) IS ACCORDING HALF LEV ETA GRID
    !                             NB: KPOS(4-5) IS NOT NEEDED
    !
    ! METHODS:
    !    - BOUNDARY CONDITION: WINDS AT TOP/BOTTOM HAVE TO BE ZERO.
    !       THIS HAS TO BE DONE EXPLICITLY FOR *WULT*,*WVLT*,
    !       WHEREAS *WELT* IS ZERO ANYWAY ON TOP/BOTTOM.
    !    !! MODIFICATION, 2.11.98, CH. REITHMEIER:
    !       TO AVOID IRREGULAR CELL DISTRIBUTIONS IN LONG RANGE INTEGRATIONS
    !       IT IS BETTER TO KEEP THE WINDS CONSTANT IN TOP/BOTTOM LAYER.
    !    - CROSS-POLE INTERPOLATION:
    !      + ZONAL WIND: DECLINES LINEARLY WHEN APPROACHING POLE
    !      + MERID.WIND: LINEAR INTERPOL. WITH POINT ON OTHER SIDE OF POLE
    !                   (NB: SIGN OF WIND ON OTHER SIDE HAS TO BE CHANGED)
    !      + VERT. WIND: AS MERID. WIND (WITHOUT CHANGE OF SIGN)
    !    - INTERPOLATION: DONE BY SUBROUTINES
    !        *ATTILA_3DINTERF* AND *ATTILA_3DINTERH*: LINEAR INTERP.
    !        IN HORIZ. DIRECTION, CUBIC INTERP. IN VERTICAL DIRECTION.
    !
    ! EXTERNALS:
    !   - GINDEXV           GET INDEX IN ASCENDING GRID
    !   - ATTILA_3DINTERFV  DO 3-D INTERPOL FOR FIELDS DEF.
    !                       ON FULL LEV ETA GRID
    !   - ATTILA_3DINTERHV  DO 3-D INTERPOL FOR FIELDS DEF.
    !                       ON HALF LEV ETA GRID

    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(OUT) :: PVELO(:,:)  ! NCEVL x 3
    REAL(dp), INTENT(IN)  :: PPOS(:,:)   ! NCEVL x NAPOS
    INTEGER,  INTENT(IN)  :: KPOS(:,:)   ! NCEVL x NNPOS

    ! LOCAL
    INTEGER,  DIMENSION(:), ALLOCATABLE :: I1,I2,I3F,I3H
    INTEGER,  DIMENSION(:), ALLOCATABLE :: I1OTHER, I2INTER
    REAL(dp), DIMENSION(:), ALLOCATABLE :: ZD1,ZD2,ZD3F,ZD3H,ZD2INTER
    REAL(dp), DIMENSION(:), ALLOCATABLE :: ZU,ZV,ZW, ZVOTHER, ZWOTHER
    INTEGER,  DIMENSION(:), ALLOCATABLE :: ICPOL
    REAL(dp), DIMENSION(:), ALLOCATABLE :: ZUFAC, ZVFAC
    REAL(dp), DIMENSION(:), ALLOCATABLE :: ZVOFAC
    INTEGER                             :: INPOL
    REAL(dp)                            :: ZDL
    INTEGER                             :: JCEVL, JNPOL

    ! INIT
    ALLOCATE(I1(NCEVL),I2(NCEVL),I3F(NCEVL),I3H(NCEVL))
    ALLOCATE(I1OTHER(NCEVL), I2INTER(NCEVL))
    ALLOCATE(ZD1(NCEVL),ZD2(NCEVL),ZD3F(NCEVL),ZD3H(NCEVL),ZD2INTER(NCEVL))
    ALLOCATE(ZU(NCEVL),ZV(NCEVL),ZW(NCEVL), ZVOTHER(NCEVL), ZWOTHER(NCEVL))
    ALLOCATE(ICPOL(NCEVL))
    ALLOCATE(ZUFAC(NCEVL), ZVFAC(NCEVL))
    ALLOCATE(ZVOFAC(NCEVL))

    !
    !*    1. INPUT FOR INTERPOLATION
    !        ----- --- -------------
    !
    INPOL= 0
    ZDL= GDLAT(NGL)
    CALL ATTILA_GINDEXV_Z(PPOS(:,3),KPOS(:,4:5),GETAF_3D,I3F)

    DO JCEVL= 1,NCEVL

       I1(JCEVL)= KPOS(JCEVL,1)
       I2(JCEVL)= KPOS(JCEVL,2)
       ZD1(JCEVL)= PPOS(JCEVL,1) - GLON(I1(JCEVL))
       ZD2(JCEVL)= PPOS(JCEVL,2) - GLAT(I2(JCEVL))

       I3H(JCEVL)=KPOS(JCEVL,3)
       IF(I3H(JCEVL) == 1) THEN         ! IF IN TOP LAYER,
          I3H(JCEVL)=2                  ! POSITION AT BOTTOM OF TOP LAYER.
          ZD3H(JCEVL)=0.0_dp
       ELSEIF(I3H(JCEVL) == NLEV) THEN  ! IF IN BOTTOM LAYER,
          I3H(JCEVL)=NLEV-1             ! POSITION AT TOP OF BOTTOM LAYER.
          GDETAH => GDETAH_3D(KPOS(JCEVL,4),:,KPOS(JCEVL,5))
          ZD3H(JCEVL)=GDETAH(NLEV-1)
       ELSE
          GETAH => GETAH_3D(KPOS(JCEVL,4),:,KPOS(JCEVL,5))
          ZD3H(JCEVL)=PPOS(JCEVL,3)-GETAH(I3H(JCEVL))
       ENDIF

       IF(I3F(JCEVL) == 0) THEN         ! IF ABOVE TOP LEVEL,
          I3F(JCEVL)=1                  ! POSITION AT TOP LEVEL.
          ZD3F(JCEVL)=0.0_dp
       ELSEIF(I3F(JCEVL) == NLEV) THEN ! IF BENEATH BOTTOM LEVEL,
          I3F(JCEVL)=NLEV-1             ! POSITION AT BOTTOM LEVEL.
          GDETAF => GDETAF_3D(KPOS(JCEVL,4),:,KPOS(JCEVL,5))
          ZD3F(JCEVL)=GDETAF(NLEV-1)
       ELSE
          GETAF => GETAF_3D(KPOS(JCEVL,4),:,KPOS(JCEVL,5))
          ZD3F(JCEVL)= PPOS(JCEVL,3) - GETAF(I3F(JCEVL))
       ENDIF

       IF(I2(JCEVL) == 1) THEN  ! NORTH-POLE
          ZD2(JCEVL)= ZD2(JCEVL)/GDLATFT2   ! IN [0;0.5)
          I2INTER(JCEVL)= 2
          ZD2INTER(JCEVL)= 0.0_dp
          I1OTHER(JCEVL)= MOD(I1(JCEVL)-1+NHLON,NLON) + 1
          !POINT AT OTHER SIDE OF POLE
          ZUFAC(JCEVL) =  2._dp*ZD2(JCEVL)
          ZVFAC(JCEVL) =0.5_dp+ZD2(JCEVL)
          ZVOFAC(JCEVL)=0.5_dp-ZD2(JCEVL)
       ELSEIF(I2(JCEVL) == NGLP1) THEN     ! SOUTH-POLE
          ZD2(JCEVL)= ZD2(JCEVL)/GDLATLT2       ! IN [0;0.5)
          I2INTER(JCEVL)= NGL
          ZD2INTER(JCEVL)= ZDL
          I1OTHER(JCEVL)= MOD(I1(JCEVL)-1+NHLON,NLON) + 1
          !POINT AT OTHER SIDE OF POLE
          ZUFAC(JCEVL) =1.0_dp-2.0_dp*ZD2(JCEVL)
          ZVFAC(JCEVL) =1.0_dp-ZD2(JCEVL)
          ZVOFAC(JCEVL)=      ZD2(JCEVL)
       ELSE
          I2INTER(JCEVL)= I2(JCEVL)
          ZD2INTER(JCEVL)= ZD2(JCEVL)
          I1OTHER(JCEVL)= -1
       ENDIF

    END DO

    DO JCEVL= 1,NCEVL
       IF((I2(JCEVL) == 1).OR.(I2(JCEVL) == NGLP1)) THEN
          INPOL= INPOL + 1
          ICPOL(INPOL)= JCEVL
       ENDIF
    ENDDO
    !
    !    2. INTERPOLATE WIND DATA, SPECIAL IN POLAR REGION.
    !       ----------- ---- ----- ------- -- ----- -------
    !
    CALL ATTILA_3DINTER_V(ZD1,ZD2INTER,ZD3F,I1,I2INTER,I3F,PWU,ZU, NLEV, KPOS)

    CALL ATTILA_3DINTER_V(ZD1,ZD2INTER,ZD3F,I1,I2INTER,I3F,PWV,ZV, NLEV, KPOS)

    CALL ATTILA_3DINTER_V(ZD1,ZD2INTER,ZD3H,I1,I2INTER,I3H,PWW,ZW, NLEVP1,KPOS)

    CALL ATTILA_3DINTER_2(ZD1,ZD2INTER,ZD3F,I1OTHER,I2INTER,I3F,PWV, &
         ZVOTHER, INPOL,ICPOL, NLEV, KPOS)

    CALL ATTILA_3DINTER_2(ZD1,ZD2INTER,ZD3H,I1OTHER,I2INTER,I3H,PWW, &
         ZWOTHER, INPOL,ICPOL, NLEVP1, KPOS)

    DO JNPOL=1,INPOL
       JCEVL= ICPOL(JNPOL)
       ZU(JCEVL)=   ZU     (JCEVL) * ZUFAC(JCEVL)
       ZV(JCEVL)=   ZV     (JCEVL) * ZVFAC(JCEVL) &
            - ZVOTHER(JCEVL) * ZVOFAC(JCEVL)
       ZW(JCEVL)=   ZW     (JCEVL) * ZVFAC(JCEVL) &
            + ZWOTHER(JCEVL) * ZVOFAC(JCEVL)
    END DO
    !
    !*    3. CONVERT FROM M/S TO DEGREE/S.
    !        ------- ---- --- -- ---------
    !
    DO JCEVL= 1,NCEVL
       PVELO(JCEVL,1)=  ZU(JCEVL) * ADEGPM &
            /COS( (PPOS(JCEVL,2)-90._dp)*ARAD )
       PVELO(JCEVL,2)= -ZV(JCEVL) * ADEGPM     ! POSITIVE= N-WIND
       PVELO(JCEVL,3)=  ZW(JCEVL)
    END DO

    ! CLEAN UP
    DEALLOCATE(I1,I2,I3F,I3H)
    DEALLOCATE(I1OTHER, I2INTER)
    DEALLOCATE(ZD1,ZD2,ZD3F,ZD3H,ZD2INTER)
    DEALLOCATE(ZU,ZV,ZW, ZVOTHER, ZWOTHER)
    DEALLOCATE(ICPOL)
    DEALLOCATE(ZUFAC, ZVFAC)
    DEALLOCATE(ZVOFAC)

  END SUBROUTINE ATTILA_VELO
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_REPOS(PPOS, LPREPOS)

    ! PURPOSE:
    !   - REPOSITION PPOS(1-3) INTO AN ALLOWED RANGE
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTREPOS)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_REPOS)
    !
    ! INTERFACE:
    !   - LTREPOS WAS CALLED FROM LTRKO4
    !   - ATTILA_REPOS IS CALLED FROM ATTILA_RKO4
    !
    ! PARAMETERS:
    !   IO PPOS      (NCEVL,NAPOS): POSITION
    !   I  LPREPOS   (NCEVL)      : TRUE,IF CELL HAS TO BE REPOSITIONED
    !
    ! REMARKS:
    !   - OUTPUT: PPOS(1-3) = HORIZONTAL + VERTIKAL POSITION
    !
    ! METHODS:
    !        RANGE OF PPOS AFTER CALLING ATTILA_REPOS:
    !          PPOS(1) IS IN [ 0 ; 360 )
    !          PPOS(2) IS IN ( 0 ; 180 )
    !            IF PPOS(2) HAPPENS TO BE EQUAL TO 0.0 OR 180.0 IT IS
    !            RESET TO 0.1 RESP. 179.9.
    !          PPOS(3) IS IN (GETAH(1) ; GETAH(NLEVP1)) = ( 0 ; 1 )

    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(INOUT) :: PPOS(:,:)   ! NCEVL x NAPOS
    LOGICAL,  INTENT(IN)    :: LPREPOS(:)  ! NCEVL

    ! LOCAL
    REAL(dp) :: P1,P2,P3
    LOGICAL  :: LOEQ, LOMES1,LOMES2,LOMES3
    INTEGER  :: JCEVL
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: IPOS

    LOMES1= .FALSE.
    LOMES2= .FALSE.
    LOMES3= .FALSE.
    LOEQ  = .FALSE.

    DO JCEVL=1,NCEVL

       P1= PPOS(JCEVL,1)
       P2= PPOS(JCEVL,2)
       !
       !*    1. RESET LONGITUDE.
       !        ----- ----------
       !
       IF((P1 >= 720.0_dp).OR.(P1 < -360.0_dp)) THEN
          LOMES1= .TRUE.
          P1= MOD(P1,360.0_dp)
       ENDIF
       IF(P1 < 0.0_dp) THEN
          P1= P1 + 360._dp
       ENDIF
       IF(P1 >= 360.0_dp) THEN
          P1= P1 - 360._dp
       ENDIF
       !           P1 IS NOW IN [0;360)
       !
       !*    2. RESET LATITUDE.
       !        ----- ---------
       !
       IF((P2 <= -180.0_dp).OR.(P2 >= 360.0_dp)) THEN
          LOMES2= .TRUE.
          P2= 90.0_dp
       ENDIF
       IF(P2 < 0.0_dp) THEN
          P2= -P2
          P1= MOD(P1+180._dp,360._dp)
       ENDIF
       IF(P2 > 180.0_dp) THEN
          P2= 360._dp - P2
          P1= MOD(P1+180._dp,360._dp)
       ENDIF
       IF(P2 <= 0.0_dp) THEN   ! EQUAL IS POSSIBLE
          P2= 0.1_dp
          LOEQ=.TRUE.
       ENDIF
       IF(P2 >= 180.0_dp) THEN ! EQUAL IS POSSIBLE
          P2= 179.9_dp
          LOEQ=.TRUE.
       ENDIF
       !            P2 IS NOW IN (0;180)
       !
       IF(LPREPOS(JCEVL)) THEN
          PPOS(JCEVL,1)= P1
          PPOS(JCEVL,2)= P2
       ENDIF

    ENDDO

    ALLOCATE(IPOS(NCEVL,NNPOS))
    CALL ATTILA_CALINDV(PPOS,IPOS,1)

    SELECT CASE(I_VERT)
    CASE(1,3)
       DO JCEVL=1,NCEVL

          GETAH  => GETAH_3D(IPOS(JCEVL,4),:,IPOS(JCEVL,5))
          GDETAH => GDETAH_3D(IPOS(JCEVL,4),:,IPOS(JCEVL,5))

          P3= PPOS(JCEVL,3)
          !
          P3=MERGE( 2*GETAH(1) - P3 + 0.01_dp*GDETAH(1), P3, (P3 <= GETAH(1)) )
          P3=MERGE( 2*GETAH(NLEV+1) - P3 - 0.01_dp*GDETAH(NLEV), P3,         &
               (P3 >= GETAH(NLEV+1)) )
          ! it is unclear, why this additional check is required here ...
          IF(P3 <= GETAH(1)) THEN
             P3= GETAH((NLEV+1)/2)
             LOMES3= .TRUE.
          ENDIF
          !            P3 IS NOW IN (GETAH(1);GETAH(NLEVP1)) = (0;1)
          !
          IF(LPREPOS(JCEVL)) THEN
             PPOS(JCEVL,3)= P3
          ENDIF

       ENDDO

    CASE(2)

       DO JCEVL=1,NCEVL

          GETAH  => GETAH_3D(IPOS(JCEVL,4),:,IPOS(JCEVL,5))
          GDETAH => GDETAH_3D(IPOS(JCEVL,4),:,IPOS(JCEVL,5))

          P3= PPOS(JCEVL,3)
          !
          P3=MERGE( 2*GETAH(1) - P3 - 0.01_dp*GDETAH(1), P3, (P3 >= GETAH(1)) )
          P3=MERGE( 2*GETAH(NLEV+1) - P3 + 0.01_dp*GDETAH(NLEV), P3,         &
               (P3 <= GETAH(NLEV+1)) )
          ! it is unclear, why this additional check is required here ...
          IF(P3 >= GETAH(1)) THEN
             P3= GETAH((NLEV+1)/2)
             LOMES3= .TRUE.
          ENDIF
          !            P3 IS NOW IN (GETAH(1);GETAH(NLEVP1)) = (3300;0)
          !
          IF(LPREPOS(JCEVL)) THEN
             PPOS(JCEVL,3)= P3
          ENDIF

       ENDDO

    END SELECT

    DEALLOCATE(IPOS)

    IF(LOMES1) THEN
       PRINT*,' === ATTILA_REPOS: POS1 >= 720 OR < -360 !'
    ENDIF
    IF(LOMES2) THEN
       PRINT*,' === ATTILA_REPOS: ERROR: PPOS2 <=-180 OR >=360 !'
       PRINT*,'     --> SET TO 90.0.'
    ENDIF
    IF(LOMES3) THEN
       PRINT*,' === ATTILA_REPOS: ERROR: VERTICAL OVERSHOOT!'
       PRINT*,'              REPOSITIONED TO MIDDLE OF ATMOSPHERE.'
    ENDIF
    IF(LOEQ) PRINT*,' === ATTILA_REPOS: PPOS(2) WAS EQUAL 0.0 OR 180.0, ', &
         'SLIGHTLY MOVED.'

  END SUBROUTINE ATTILA_REPOS
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_CALINDV(PPOS, KPOS, NO)

    ! PURPOSE:
    !   - CALCULATION OF INDICES/INTEGER POSITIONS OF CELLS
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LTCALINDV[,2,3])
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_CALINDV)
    !
    ! INTERFACE:
    !   - LTCALINDV[,2,3] WERE CALLED FROM LTRKO4
    !   - ATTILA_CALINDV IS CALLED FROM ATTILA_RKO4
    !
    ! PARAMETERS:
    !   I  PPOS    (NCEVL,NAPOS)    POSITION OF CELL
    !   O  KPOS    (NCEVL,NNPOS)    INTEGER POSITION OF CELL
    !   I  NO                       SWITCH (1,2,3)
    !
    ! REMARKS:
    !   - OUTPUT:
    !        KPOS(1-5) INTEGER POSITION (IF NO=1)
    !        KPOS(1-3) INTEGER POSITION (IF NO=2)
    !        KPOS(3-5) INTEGER POSITION (IF NO=3)
    !   - INPUT:
    !        PPOS(1-3) HORIZONTAL AND ETA POSITION
    !
    ! EXTERNALS:
    !   - GINDEXV    GET INTEGER POSITION IN AN ASCENDING GRID

    IMPLICIT NONE

    ! I/O
    REAL(dp),    INTENT(IN) :: PPOS(:,:) ! NCEVL x NAPOS
    INTEGER, INTENT(OUT)    :: KPOS(:,:) ! NCEVL x NAPOS
    INTEGER, INTENT(IN)     :: NO        ! 1,2,3 (REPLACE ENTRY)

    ! LOCAL
    INTEGER :: JCEVL

    IF (NO ==1) THEN
       DO JCEVL= 1,NCEVL
          KPOS(JCEVL,4)= MOD( INT((PPOS(JCEVL,1)+GDHLON)/GDLON), NLON) + 1
       ENDDO
       CALL ATTILA_GINDEXV(PPOS(:,2), GLAT2, KPOS(:,5))
    END IF

    ! ------------

    IF ((NO == 1).OR.(NO == 2)) THEN
       DO JCEVL= 1,NCEVL
          KPOS(JCEVL,1)= INT( PPOS(JCEVL,1)/GDLON ) + 1
       ENDDO
       CALL ATTILA_GINDEXV(PPOS(:,2),GLAT, KPOS(:,2))
       CALL ATTILA_GINDEXV_Z(PPOS(:,3),KPOS(:,4:5),GETAH_3D,KPOS(:,3))
       RETURN
    END IF

    !------------------

    DO JCEVL= 1,NCEVL
       KPOS(JCEVL,4)= MOD( INT((PPOS(JCEVL,1)+GDHLON)/GDLON), NLON) + 1
    ENDDO
    CALL ATTILA_GINDEXV(PPOS(:,2),GLAT2,KPOS(:,5))
    CALL ATTILA_GINDEXV_Z(PPOS(:,3),KPOS(:,4:5),GETAH_3D,KPOS(:,3))

    !------------------

  END SUBROUTINE ATTILA_CALINDV
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_GINDEXV(PPOS,PGRID,KPOS) ! KDIM = SIZE(PGRID)
                                             ! K = SIZE(PPOS)

    ! PURPOSE:
    !   - DETERMINE THE MAXIMAL INDEX OF AN ASCENDING GRID FOR WHICH
    !     PPOS IS  >=  GRID VALUE.
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (GINDEXV[,2])
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_GINDEXV)
    !     - KDIM REMOVED FROM PARAMETERS
    !     - K REMOVED FROM PARAMETERS
    !
    ! INTERFACE:
    !   - GINDEXV[,2] WERE CALLED FROM LTCALINDV AND LTVELO
    !   - ATTILA_GINDEXV IS CALLED FROM ATTILA_CALINDV AND ATTILA_VELO
    !
    ! PARAMETERS:
    !   I  PPOS  (K)     POSITION WITHIN GRID
    !   I  PGRID (KDIM)  GRID
    !   O  KPOS  (K)     CALCULATED INDEX
    !
    ! METHODS:
    !   - SIMPLE SEARCH LOOP
    !
    ! REMARKS:
    !      RESULT: PPOS LIES IN THE INTERVAL
    !                             [ PGRID(KPOS);PGRID(KPOS+1) )
    !              KPOS  =0    IF PPOS< PGRID(1)
    !                    =KDIM IF PPOS>=PGRID(KDIM)
    !
    ! REFERENCES:
    !      BOOK: *W. *H. *PRESS ET AL., *NUMERICAL *RECIPES, P 90.

    IMPLICIT NONE

    ! I/O
    REAL(dp),    INTENT(IN)  :: PGRID(:)
    REAL(dp),    INTENT(IN)  :: PPOS(:)
    INTEGER,     INTENT(OUT) :: KPOS(:)

    ! LOCAL
    INTEGER                            :: KDIM, K
    INTEGER                            :: J, JK
    LOGICAL, DIMENSION(:), ALLOCATABLE :: LOCKED
    LOGICAL                            :: LO1

    ! INIT
    KDIM = SIZE(PGRID)
    K    = SIZE(PPOS)
    !
    ALLOCATE(LOCKED(K))
    LOCKED(:)=.FALSE.

    DO JK= 1,K
       LO1= (PPOS(JK) >= PGRID(KDIM))
       KPOS(JK)=   MERGE( KDIM ,0      ,LO1)
!!$       LOCKED(JK)= MERGE(.TRUE.,.FALSE.,LO1)
       LOCKED(JK) = LO1
    ENDDO

    DO J=1,KDIM
       DO JK= 1,K
          LO1= (PPOS(JK) < PGRID(J))
          KPOS(JK)=  MERGE(J-1   ,KPOS(JK)  ,(LO1.AND..NOT.LOCKED(JK)))
          LOCKED(JK)=MERGE(.TRUE.,LOCKED(JK),LO1)
       ENDDO
    ENDDO

    ! CLEAN UP
    DEALLOCATE(LOCKED)

  END SUBROUTINE ATTILA_GINDEXV
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_GINDEXV_Z(PPOS,ZKPOS,PGRID,KPOS)

    IMPLICIT NONE

    ! I/O
    REAL(dp),    INTENT(IN)  :: PPOS(:)
    INTEGER,     INTENT(IN)  :: ZKPOS(:,:)
    REAL(dp),    INTENT(IN)  :: PGRID(:,:,:)
    INTEGER,     INTENT(OUT) :: KPOS(:)

    ! LOCAL
    INTEGER                            :: KDIM, C
    INTEGER                            :: J, JC
    LOGICAL, DIMENSION(:), ALLOCATABLE :: LOCKED
    LOGICAL                            :: LO1

    ! INIT
    KDIM = SIZE(PGRID,2)
    C    = SIZE(PPOS)
    !
    ALLOCATE(LOCKED(C))
    LOCKED(:)=.FALSE.

    SELECT CASE(I_VERT)
    CASE(1,3)

       DO JC= 1,C
          LO1= (PPOS(JC) >= PGRID(ZKPOS(JC,1),KDIM,ZKPOS(JC,2)))
          KPOS(JC)=   MERGE( KDIM ,0      ,LO1)
!!$          LOCKED(JC)= MERGE(.TRUE.,.FALSE.,LO1)
          LOCKED(JC)= LO1
       ENDDO

       DO J=1,KDIM
          DO JC= 1,C
             LO1= (PPOS(JC) < PGRID(ZKPOS(JC,1),J,ZKPOS(JC,2)))
!!$             KPOS(JC)=  MERGE(J-1   ,KPOS(JC)  ,LO1.AND..NOT.LOCKED(JC))
             IF (LO1.AND..NOT.LOCKED(JC)) KPOS(JC)=J-1
             LOCKED(JC)=MERGE(.TRUE.,LOCKED(JC),LO1)
          ENDDO
       ENDDO

    CASE(2)

       DO JC= 1,C
          LO1= (PPOS(JC) > PGRID(ZKPOS(JC,1),1,ZKPOS(JC,2)))
          KPOS(JC)=   MERGE( 1 ,0      ,LO1)
          LOCKED(JC)= MERGE(.TRUE.,.FALSE.,LO1)
       END DO

       DO J=KDIM,1,-1
          DO JC= 1,C
             LO1= (PPOS(JC) <= PGRID(ZKPOS(JC,1),J,ZKPOS(JC,2)))
             KPOS(JC)=  MERGE(J   ,KPOS(JC)  ,LO1.AND..NOT.LOCKED(JC))
             LOCKED(JC)=MERGE(.TRUE.,LOCKED(JC),LO1)
          END DO
       END DO

    END SELECT

    ! CLEAN UP
    DEALLOCATE(LOCKED)

  END SUBROUTINE ATTILA_GINDEXV_Z
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_3DINTER_V(PD1, PD2, PD3  &
       , K1, K2, K3, PF, PFINTER             &
       , KLEV, KPOS)                     ! op_sb_20161202

    ! PURPOSE:
    !   -  DO 3-D INTERPOLATION ON FIELD *PF*, WHICH IS SUPPOSED TO BE
    !      A 3-D MET. FIELD DEFINED ON FULL ETA LEVELS.
    !      NB: REQUIRED POSITION MUST NOT BE IN TOP/BOTTOM LAYER
    !          OR IN POLAR REGION!
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LT3DINTER[F,H]V)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_3DINTER_V)
    !
    ! INTERFACE:
    !   - LT3DINTERFV WAS CALLED FROM LTVELO
    !   - LT3DINTERHV WAS CALLED FROM LTVELO
    !   - ATTILA_3DINTER_V IS CALLED FROM ATTILA_VELO
    !
    ! PARAMETERS:
    !   I  PD1,PD2,PD3: DISTANCE FROM K1,K2,K3 (IN DEG./ETA)
    !   I  K1,K2,K3   : "NORTH-WEST-UPPER POSITION" OF VALUES TO
    !                    INTERPOLATE (= INTEGER POSITION OF REQUIRED POINT)
    !               NB: K1 IS INTEGER POSITION ACCORDING TO GRID GLON.
    !               NB: K2 IS INTEGER POSITION ACCORDING TO GRID GLAT.
    !               NB: K3 IS INTEGER POSITION ACCORDING TO GRID GETAF.
    !                                       (FULL LEVEL GRID)
    !   I  PF (NLON,NLEV,NGL): FIELD TO BE INTERPOLATED
    !   O  PFINTER           : INTERPOLATED VALUES
    !
    ! REMARKS:
    !     ALLOWED RANGES:
    !       K1 = 1..NLON
    !       K2 = 2..NGL
    !       K3 = 1..NLEV-1
    !       PD1 IN [0;GDLON]
    !       PD2 IN [0;GDLAT(K2)]
    !       PD3 IN [0;GDETAF(K3)]
    !
    ! METHODS:
    !       LINEAR INTERPOLATION IN HORIZONTAL DIRECTION,
    !       CUBIC  INTERPOLATION IN VERTICAL   DIRECTION.
    !        THE CUBIC INTERPOLATION POLYNOMIAL F(X) HAS THE FOLLOWING
    !        PROPERTIES:
    !         - F (X(2))= P(2)
    !         - F (X(3))= P(3)
    !         - F'(X(2))= D(2)
    !         - F'(X(3))= D(3)
    !         WHERE P(I) IS THE GIVEN VALUE AT X(I),
    !               D(I)= (P(I+1)-P(I-1)) / (X(I+1)-X(I-1))
    !                    IS THE DERIVATIVE AT X(I),
    !               I=2,3 .

    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(IN)  :: PD1(:),PD2(:),PD3(:) ! NCEVL
    INTEGER,  INTENT(IN)  :: K1(:),K2(:),K3(:)    ! NCEVL
    REAL(dp), INTENT(IN)  :: PF(:,:,:)            ! NLON x NLEV x NGL
    REAL(dp), INTENT(OUT) :: PFINTER(:)           ! NCEVL
    INTEGER,  INTENT(IN)  :: KLEV                 ! =NLEV   (ATTILA_3DINTERFV)
                                                  ! =NLEVP1 (ATTILA_3DINTERHV)
    INTEGER,  INTENT(IN)  :: KPOS(:,:)            ! NCEVL x NNPOS

    ! LOCAL
    INTEGER  :: JCEVL
    INTEGER  :: J, I1P1,I2M1,I3, I3M2
    REAL(dp) :: ZP(4,4), ZR(4), ZD1, Z1,Z2
    !  VARIABLES FOR CUBIC INTERPOLATION
    REAL(dp) :: ZA(4),Z2D3,ZD3T3,ZDH1,ZDH2 !,Z3
    INTEGER  :: ZK3
    REAL(dp) :: ZX(4), Z3
    REAL(dp), DIMENSION(:,:,:), POINTER :: geta_3d_l
    REAL(dp), DIMENSION(:,:,:), POINTER :: geta_3d
    REAL(dp), DIMENSION(:),     POINTER :: geta_l
    REAL(dp), DIMENSION(:),     POINTER :: geta

    ! p-grid for interpolation of velocity for all i_vert
    if (klev == nlev) then       ! grid on full levels
       geta_3d_l =>  GETA_3DF_L_P(:,:,:)
       geta_3d   =>  GETA_3DF_P(:,:,:)
    else                         ! grid on half levels
       geta_3d_l =>  GETA_3DH_L_P(:,:,:)
       geta_3d   =>  GETA_3DH_P(:,:,:)
    endif

    DO JCEVL=1,NCEVL           ! cell loop

          ! p-grid where the parcel lies
          GETA    => GETA_3D(KPOS(JCEVL,4),:,KPOS(JCEVL,5))
          GETA_L  => GETA_3D_L(KPOS(JCEVL,4),:,KPOS(JCEVL,5))

       DO ZK3=1,NLEV-1  ! THIS IS THE ALLOWED RANGE OF THE DUMMY
                        ! ARGUMENT K3 IN ROUTINE ATTILA_3DINTERF

          I3M2= ZK3 - 2
          DO J=1,4   ! LOOP OVER LEVELS
            I3= MIN(NLEV,MAX(1,I3M2+J))
!             ZX(J)= GETA(I3)
             ! for interpolation on parcel (I_VERT dependent grid)
             if (klev == nlev) then       ! grid on full levels
              ZX(J) = GETAF_3D(KPOS(JCEVL,4),I3,KPOS(JCEVL,5))
             ELSE
              ZX(J) = GETAH_3D(KPOS(JCEVL,4),I3,KPOS(JCEVL,5))
             ENDIF
          ENDDO
          !
          Z3       = ZX(3) - ZX(2)
          !A3DZ3(ZK3)= Z3
          A3DRZ3(ZK3)= 1._dp/Z3
          A3DDH1(ZK3)= Z3 / (ZX(3) - ZX(1))
          A3DDH2(ZK3)= Z3 / (ZX(4) - ZX(2))
          !
       END DO
       !
       !*   1. GET SOME INDICES.
       !       --- ---- --------
       !
       I1P1= MOD(K1(JCEVL),NLON) + 1
       I2M1= K2(JCEVL) - 1
       I3M2= K3(JCEVL) - 2

       !
         DO J=1,4   ! LOOP OVER LEVELS
          I3= MIN(KLEV,MAX(1,I3M2+J))
          if (I3 > 1 .and. I3 < KLEV) then
            if (geta_3d(K1(JCEVL),I3,I2M1 ) < geta(I3) ) then
                ZP(J,1) = PF(K1(JCEVL),I3,I2M1) + (PF(K1(JCEVL),I3+1,I2M1)-   &
                        PF(K1(JCEVL),I3,I2M1))/(geta_3d_l(K1(JCEVL),I3+1,I2M1)&
                        -geta_3d_l(K1(JCEVL),I3,I2M1)) *                      &
                         abs(geta_l(I3) - geta_3d_l(K1(JCEVL),I3,I2M1))
             if (i3 == 2) then  ! extra treatment, because log(geta) not defined at I3=1
                ZP(J,1) = PF(K1(JCEVL),I3,I2M1) + (PF(K1(JCEVL),I3+1,I2M1)-   &
                        PF(K1(JCEVL),I3,I2M1))/(geta_3d(K1(JCEVL),I3+1,I2M1)  &
                        -geta_3d(K1(JCEVL),I3,I2M1)) *  &
                         abs(geta(I3) - geta_3d(K1(JCEVL),I3,I2M1))
             endif
            elseif  (geta_3d(K1  (JCEVL),I3,I2M1 ) > geta(I3) ) then
                 ZP(J,1) = PF(K1(JCEVL),I3,I2M1) + (PF(K1(JCEVL),I3-1,I2M1)-   &
                        PF(K1(JCEVL),I3,I2M1))/(geta_3d_l(K1(JCEVL),I3-1,I2M1) &
                        -geta_3d_l(K1(JCEVL),I3,I2M1)) *                       &
                          abs(geta_l(I3) - geta_3d_l(K1(JCEVL),I3,I2M1))
              if (i3 == 2) then  ! extra treatment, because log(geta) not defined at I3=1
                 ZP(J,1) = PF(K1(JCEVL),I3,I2M1) + (PF(K1(JCEVL),I3-1,I2M1)-   &
                        PF(K1(JCEVL),I3,I2M1))/(geta_3d(K1(JCEVL),I3-1,I2M1)   &
                        -geta_3d(K1(JCEVL),I3,I2M1)) *                         &
                          abs(geta(I3) - geta_3d(K1(JCEVL),I3,I2M1))
              endif
            else
                ZP(J,1)= PF(K1(JCEVL),I3,I2M1)
            endif

            if (geta_3d(K1(JCEVL),I3, K2(JCEVL)) < geta(I3) ) then
                 ZP(J,2) = PF(K1(JCEVL),I3,K2(JCEVL))+(PF(K1(JCEVL),I3+1,K2(JCEVL))-    &
                        PF(K1(JCEVL),I3,K2(JCEVL)))/(geta_3d_l(K1(JCEVL),I3+1,K2(JCEVL))&
                        -geta_3d_l(K1(JCEVL),I3,K2(JCEVL)) )  *                         &
                          abs(geta_l(I3) - geta_3d_l(K1(JCEVL),I3,K2(JCEVL)))

             if (i3 == 2) then  ! extra treatment, because log(geta) not defined at I3=1
               ZP(J,2) = PF(K1(JCEVL),I3,K2(JCEVL)) + (PF(K1(JCEVL),I3+1,K2(JCEVL))-    &
                        PF(K1(JCEVL),I3,K2(JCEVL)))/(geta_3d(K1(JCEVL),I3+1,K2(JCEVL))  &
                        -geta_3d(K1(JCEVL),I3,K2(JCEVL)) )  *                           &
                          abs(geta(I3) - geta_3d(K1(JCEVL),I3,K2(JCEVL)))
             endif
             elseif  (geta_3d(K1  (JCEVL),I3,K2(JCEVL) ) > geta(I3) ) then
                  ZP(J,2) = PF(K1(JCEVL),I3,K2(JCEVL)) + (PF(K1(JCEVL),I3-1,K2(JCEVL))-   &
                        PF(K1(JCEVL),I3,K2(JCEVL)))/(geta_3d_l(K1(JCEVL),I3-1,K2(JCEVL))  &
                        -geta_3d_l(K1(JCEVL),I3,K2(JCEVL))) *                             &
                           abs(geta_l(I3) - geta_3d_l(K1(JCEVL),I3,K2(JCEVL)))

              if (i3 == 2) then ! extra treatment, because log(geta) not defined at I3=1
                   ZP(J,2) = PF(K1(JCEVL),I3,K2(JCEVL)) + (PF(K1(JCEVL),I3-1,K2(JCEVL))-   &
                        PF(K1(JCEVL),I3,K2(JCEVL)))/(geta_3d(K1(JCEVL),I3-1,K2(JCEVL))     &
                        -geta_3d(K1(JCEVL),I3,K2(JCEVL))) *                                &
                           abs(geta(I3) - geta_3d(K1(JCEVL),I3,K2(JCEVL)))
              endif
             else
                ZP(J,2)= PF(K1(JCEVL),I3, K2(JCEVL))
             endif

             if (geta_3d(I1P1,I3,I2M1) < geta(I3) ) then
                ZP(J,3) = PF(I1P1,I3,I2M1) + (PF(I1P1,I3+1,I2M1)-   &
                        PF(I1P1,I3,I2M1))/(geta_3d_l(I1P1,I3+1,I2M1)&
                        -geta_3d_l(I1P1,I3,I2M1)) *                 &
                         abs(geta_l(I3) - geta_3d_l(I1P1,I3,I2M1))

             if (i3 == 2) then  ! extra treatment, because log(geta) not defined at I3=1
                ZP(J,3) = PF(I1P1,I3,I2M1) + (PF(I1P1,I3+1,I2M1)-   &
                        PF(I1P1,I3,I2M1))/(geta_3d(I1P1,I3+1,I2M1)  &
                        -geta_3d(I1P1,I3,I2M1)) *                   &
                         abs(geta(I3) - geta_3d(I1P1,I3,I2M1))
             endif
              elseif  (geta_3d(I1P1,I3,I2M1) > geta(I3) ) then
                 ZP(J,3) = PF(I1P1,I3,I2M1) + (PF(I1P1,I3-1,I2M1)-   &
                        PF(I1P1,I3,I2M1))/(geta_3d_l(I1P1,I3-1,I2M1) &
                        -geta_3d_l(I1P1,I3,I2M1)) *                  &
                         abs(geta_l(I3) - geta_3d_l(I1P1,I3,I2M1))
             if (i3 == 2) then  ! extra treatment, because log(geta) not defined at I3=1
                 ZP(J,3) = PF(I1P1,I3,I2M1) + (PF(I1P1,I3-1,I2M1)-   &
                        PF(I1P1,I3,I2M1))/(geta_3d(I1P1,I3-1,I2M1)   &
                        -geta_3d(I1P1,I3,I2M1)) * &
                         abs(geta(I3) - geta_3d(I1P1,I3,I2M1))
             endif
              else
                ZP(J,3)= PF(I1P1,I3,I2M1)
              endif

              if (geta_3d(I1P1,I3,K2(JCEVL)) < geta(I3) ) then
                ZP(J,4) = PF(I1P1,I3,K2(JCEVL)) + (PF(I1P1,I3+1,K2(JCEVL))-   &
                        PF(I1P1,I3,K2(JCEVL)))/(geta_3d_l(I1P1,I3+1,K2(JCEVL))&
                        -geta_3d_l(I1P1,I3,K2(JCEVL))) * &
                          abs(geta_l(I3) - geta_3d_l(I1P1,I3,K2(JCEVL)))

              if (i3 == 2) then   ! extra treatment, because log(geta) not defined at I3=1
                ZP(J,4) = PF(I1P1,I3,K2(JCEVL)) + (PF(I1P1,I3+1,K2(JCEVL))-   &
                        PF(I1P1,I3,K2(JCEVL)))/(geta_3d(I1P1,I3+1,K2(JCEVL))  &
                        -geta_3d(I1P1,I3,K2(JCEVL))) *                   &
                          abs(geta(I3) - geta_3d(I1P1,I3,K2(JCEVL)))
             endif

              elseif  (geta_3d(I1P1,I3,K2(JCEVL) ) > geta(I3) ) then
                 ZP(J,4) = PF(I1P1,I3,K2(JCEVL)) + (PF(I1P1,I3-1,K2(JCEVL))-   &
                        PF(I1P1,I3,K2(JCEVL)))/(geta_3d_l(I1P1,I3-1,K2(JCEVL)) &
                        -geta_3d_l(I1P1,I3,K2(JCEVL))) *                       &
                         abs(geta_l(I3) - geta_3d_l(I1P1,I3,K2(JCEVL)))

              if (i3 == 2) then  ! extra treatment, because log(geta) not defined at I3=1
                 ZP(J,4) = PF(I1P1,I3,K2(JCEVL)) + (PF(I1P1,I3-1,K2(JCEVL))-   &
                        PF(I1P1,I3,K2(JCEVL)))/(geta_3d(I1P1,I3-1,K2(JCEVL))   &
                        -geta_3d(I1P1,I3,K2(JCEVL))) *                         &
                         abs(geta(I3) - geta_3d(I1P1,I3,K2(JCEVL)))
              endif
              else
                ZP(J,4)= PF(I1P1,I3,K2(JCEVL))
              endif

        else
          I3= MIN(KLEV,MAX(1,I3M2+J))
          ZP(J,1)= PF(K1  (JCEVL),I3,I2M1       )
          ZP(J,2)= PF(K1  (JCEVL),I3,K2  (JCEVL))
          ZP(J,3)= PF(I1P1       ,I3,I2M1       )
          ZP(J,4)= PF(I1P1       ,I3,K2  (JCEVL))
        ENDIF
        ENDDO ! LOOP OVER LEVELS(4)

       !
       !*   3. DO CUBIC INTERPOLATION IN VERTICAL DIRECTION.
       !       -- ----- ------------- -- -------- ----------
       !
!       Z3  = A3DZ3(K3(JCEVL))                   ! = ZX(3) - ZX(2)
       Z2D3= PD3(JCEVL) * A3DRZ3(K3(JCEVL))     ! = PD3/Z3
       ZDH1= A3DDH1(K3(JCEVL))                  ! = Z3/(ZX(3)-ZX(1))
       ZDH2= A3DDH2(K3(JCEVL))                  ! = Z3/(ZX(4)-ZX(2))
       DO J=1,4   ! LOOP OVER HORIZONTAL POSITIONS
          ZA(4)= ZP(2,J)
          ZA(3)= (ZP(3,J)-ZP(1,J))*ZDH1     !DERIV AT POINT 2 X Z3
          ZD3T3= (ZP(4,J)-ZP(2,J))*ZDH2     !DERIV AT POINT 3 X Z3
          ZA(2)= 3*ZP(3,J) - ZD3T3 - 2*ZA(3) - 3*ZA(4)
          ZA(1)= ZP(3,J) - ZA(2) - ZA(3) - ZA(4)
          ZR(J)= ZA(4) + Z2D3*(ZA(3)+Z2D3*(ZA(2)+Z2D3*ZA(1)))
       END DO
       !
       !*   4. DO LINEAR INTERPOLATION IN HORIZONTAL DIRECTION.
       !       -- ------ ------------- -- ---------- ----------
       !
       ZD1= PD1(JCEVL)/GDLON
       Z1= ZR(1) + (ZR(3)-ZR(1))*ZD1
       Z2= ZR(2) + (ZR(4)-ZR(2))*ZD1
       PFINTER(JCEVL)= Z1 + (Z2-Z1)*PD2(JCEVL)/GDLAT(K2(JCEVL))

    END DO

  END SUBROUTINE ATTILA_3DINTER_V
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  SUBROUTINE ATTILA_3DINTER_2(PD1, PD2, PD3 &
       , K1, K2, K3, PF, PFINTER            &
       , KNPOL, KCPOL                       &
       , KLEV, KPOS)

    ! PURPOSE:
    !   - DO 3-D INTERPOLATION ON FIELD *PF*, WHICH IS SUPPOSED TO BE
    !     A 3-D MET. FIELD DEFINED ON ETA LEVELS.
    !     NB: REQUIRED POSITION MUST NOT BE IN TOP/BOTTOM LAYER
    !         OR IN POLAR REGION!
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98, (LT3DINTER[F,H]2)
    !   M. Traub, P. Joeckel, MPI Mainz                 (ATTILA_3DINTER_2)
    !
    ! INTERFACE:
    !   - LT3DINTERF2 WAS CALLED FROM ATTILA_VELO
    !   - LT3DINTERH2 WAS CALLED FROM ATTILA_VELO
    !   - ATTILA_3DINTER_2 IS CALLED FROM ATTILA_VELO
    !
    ! PARAMETERS:
    !    I  PD1,PD2,PD3: DISTANCE FROM K1,K2,K3 (IN DEG./ETA)
    !    I  K1,K2,K3 : "NORTH-WEST-UPPER POSITION" OF VALUES TO
    !                   INTERPOLATE (= INTEGER POSITION OF REQUIRED POINT)
    !         NB: K1 IS INTEGER POSITION ACCORDING TO GRID GLON.
    !         NB: K2 IS INTEGER POSITION ACCORDING TO GRID GLAT.
    !         NB: K3 IS INTEGER POSITION ACCORDING TO GRID GETAF.
    !                                       (LEVEL GRID)
    !    I  PF(NLON,NLEV,NGL): FIELD TO BE INTERPOLATED
    !    O PFINTER           : INTERPOLATED VALUES
    !    I  KCPOL(NCEVL)     : NUMBERS OF CELLS WHICH ARE INTERPOLATED
    !
    ! REMARKS:
    !     ALLOWED RANGES:
    !       K1 = 1..NLON
    !       K2 = 2..NGL
    !       K3 = 1..NLEV-1
    !       PD1 IN [0;GDLON]
    !       PD2 IN [0;GDLAT(K2)]
    !       PD3 IN [0;GDETAF(K3)]
    !
    ! METHODS:
    !       LINEAR INTERPOLATION IN HORIZONTAL DIRECTION,
    !       CUBIC  INTERPOLATION IN VERTICAL   DIRECTION.
    !        THE CUBIC INTERPOLATION POLYNOMIAL F(X) HAS THE FOLLOWING
    !        PROPERTIES:
    !         - F (X(2))= P(2)
    !         - F (X(3))= P(3)
    !         - F'(X(2))= D(2)
    !         - F'(X(3))= D(3)
    !         WHERE P(I) IS THE GIVEN VALUE AT X(I),
    !               D(I)= (P(I+1)-P(I-1)) / (X(I+1)-X(I-1))
    !                    IS THE DERIVATIVE AT X(I),
    !               I=2,3 .

    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(IN)  :: PD1(:),PD2(:),PD3(:) ! NCEVL
    INTEGER,  INTENT(IN)  :: K1(:),K2(:),K3(:)    ! NCEVL
    REAL(dp), INTENT(IN)  :: PF(:,:,:)            ! NLON x NLEV x NGL
    REAL(dp), INTENT(OUT) :: PFINTER(:)           ! NCEVL
    INTEGER,  INTENT(IN)  :: KNPOL
    INTEGER,  INTENT(IN)  :: KCPOL(:)             ! NCEVL
    INTEGER,  INTENT(IN)  :: KLEV                 ! =NLEV   (ATTILA_3DINTERF2)
                                                  ! =NLEVP1 (ATTILA_3DINTERH2)
    INTEGER,  INTENT(IN)  :: KPOS(:,:)            ! NCEVL x NNPOS

    ! LOCAL
    INTEGER  :: JNPOL, JCEVL
    INTEGER  :: J, I1P1,I2M1,I3, I3M2
    REAL(dp) :: ZP(4,4), ZR(4), ZD1, Z1,Z2
    !  VARIABLES FOR CUBIC INTERPOLATION
    REAL(dp) :: ZA(4),Z2D3,ZD3T3,ZDH1,ZDH2!,Z3
    INTEGER  :: ZK3
    REAL(dp) :: ZX(4), Z3
    REAL(dp), DIMENSION(:,:,:), POINTER :: geta_3d_l
    REAL(dp), DIMENSION(:,:,:), POINTER :: geta_3d
    REAL(dp), DIMENSION(:),     POINTER :: geta_l
    REAL(dp), DIMENSION(:),     POINTER :: geta

    ! p-grid for interpolation of velocity for all i_vert
    if (klev == nlev) then       ! grid on full levels
       geta_3d_l =>  GETA_3DF_L_P(:,:,:)
       geta_3d   =>  GETA_3DF_P(:,:,:)
    else
       geta_3d_l =>  GETA_3DH_L_P(:,:,:)
       geta_3d   =>  GETA_3DH_P(:,:,:)
    endif

    DO JNPOL=1,KNPOL
       JCEVL= KCPOL(JNPOL)

          ! p-grid where the parcel lies
          GETA    => GETA_3D(KPOS(JCEVL,4),:,KPOS(JCEVL,5))
          GETA_L  => GETA_3D_L(KPOS(JCEVL,4),:,KPOS(JCEVL,5))

       DO ZK3=1,NLEV-1  ! THIS IS THE ALLOWED RANGE OF THE DUMMY
                       ! ARGUMENT K3 IN ROUTINE ATTILA_3DINTERF
          I3M2= ZK3 - 2
          DO J=1,4   ! LOOP OVER LEVELS
             I3= MIN(NLEV,MAX(1,I3M2+J))
!             ZX(J)= GETA(I3)
          ! eta or theta grid dependent on I_VERT, set in
          ! attila_init_geta
             if (klev == nlev) then       ! grid on full levels
              ZX(J) = GETAF_3D(KPOS(JCEVL,4),I3,KPOS(JCEVL,5))
             else                         ! grid on half levels
              ZX(J) = GETAH_3D(KPOS(JCEVL,4),I3,KPOS(JCEVL,5))
             endif
          ENDDO
          !
          Z3       = ZX(3) - ZX(2)
          !A3DZ3(ZK3)= Z3
          A3DRZ3(ZK3)= 1._dp/Z3
          A3DDH1(ZK3)= Z3 / (ZX(3) - ZX(1))
          A3DDH2(ZK3)= Z3 / (ZX(4) - ZX(2))
          !
       END DO
       !
       !*   1. GET SOME INDICES.
       !       --- ---- --------
       !
       I1P1= MOD(K1(JCEVL),NLON) + 1
       I2M1= K2(JCEVL) - 1
       I3M2= K3(JCEVL) - 2
       !
       !*   2. SET UP THE FOUR SETS OF PROFILES FOR INTERPOLATION.
       !       --- -- --- ---- ---- -- -------- --- --------------
       !        REPEAT BOTTOM OR TOP IF NEEDED.

        DO J=1,4   ! LOOP OVER LEVELS
          I3= MIN(KLEV,MAX(1,I3M2+J))
          if (I3 > 1 .and. I3 < KLEV) then
            if (geta_3d(K1(JCEVL),I3,I2M1 ) < geta(I3) ) then
                ZP(J,1) = PF(K1(JCEVL),I3,I2M1) + (PF(K1(JCEVL),I3+1,I2M1)-   &
                        PF(K1(JCEVL),I3,I2M1))/(geta_3d_l(K1(JCEVL),I3+1,I2M1)&
                        -geta_3d_l(K1(JCEVL),I3,I2M1))*                       &
                        abs(geta_l(I3) - geta_3d_l(K1(JCEVL),I3,I2M1))
             if (i3 == 2) then
                ZP(J,1) = PF(K1(JCEVL),I3,I2M1) + (PF(K1(JCEVL),I3+1,I2M1)-   &
                        PF(K1(JCEVL),I3,I2M1))/(geta_3d(K1(JCEVL),I3+1,I2M1)  &
                        -geta_3d(K1(JCEVL),I3,I2M1))*                         &
                        abs(geta(I3) - geta_3d(K1(JCEVL),I3,I2M1))
             endif
            elseif  (geta_3d(K1  (JCEVL),I3,I2M1 ) > geta(I3) ) then
                 ZP(J,1) = PF(K1(JCEVL),I3,I2M1) + (PF(K1(JCEVL),I3-1,I2M1)-   &
                        PF(K1(JCEVL),I3,I2M1))/(geta_3d_l(K1(JCEVL),I3-1,I2M1) &
                        -geta_3d_l(K1(JCEVL),I3,I2M1)) *                       &
                         abs(geta_l(I3) - geta_3d_l(K1(JCEVL),I3,I2M1))
              if (i3 == 2) then
                 ZP(J,1) = PF(K1(JCEVL),I3,I2M1) + (PF(K1(JCEVL),I3-1,I2M1)-   &
                        PF(K1(JCEVL),I3,I2M1))/(geta_3d(K1(JCEVL),I3-1,I2M1)   &
                        -geta_3d(K1(JCEVL),I3,I2M1)) *                         &
                         abs(geta(I3) - geta_3d(K1(JCEVL),I3,I2M1))
              endif
            else
                ZP(J,1)= PF(K1(JCEVL),I3,I2M1)
            endif

            if (geta_3d(K1(JCEVL),I3, K2(JCEVL)) < geta(I3) ) then
                 ZP(J,2) = PF(K1(JCEVL),I3,K2(JCEVL)) + (PF(K1(JCEVL),I3+1,K2(JCEVL))-   &
                        PF(K1(JCEVL),I3,K2(JCEVL)))/(geta_3d_l(K1(JCEVL),I3+1,K2(JCEVL)) &
                        -geta_3d_l(K1(JCEVL),I3,K2(JCEVL)) ) *                           &
                        abs(geta_l(I3) - geta_3d_l(K1(JCEVL),I3,K2(JCEVL)))

             if (i3 == 2) then
                 ZP(J,2) = PF(K1(JCEVL),I3,K2(JCEVL)) + (PF(K1(JCEVL),I3+1,K2(JCEVL))-   &
                        PF(K1(JCEVL),I3,K2(JCEVL)))/(geta_3d(K1(JCEVL),I3+1,K2(JCEVL))   &
                        -geta_3d(K1(JCEVL),I3,K2(JCEVL)) ) *                             &
                        abs(geta(I3) - geta_3d(K1(JCEVL),I3,K2(JCEVL)))
              endif
             elseif  (geta_3d(K1  (JCEVL),I3,K2(JCEVL) ) > geta(I3) ) then
                 ZP(J,2) = PF(K1(JCEVL),I3,K2(JCEVL)) + (PF(K1(JCEVL),I3-1,K2(JCEVL))-    &
                        PF(K1(JCEVL),I3,K2(JCEVL)))/(geta_3d_l(K1(JCEVL),I3-1,K2(JCEVL))  &
                        -geta_3d_l(K1(JCEVL),I3,K2(JCEVL))) *                             &
                          abs(geta_l(I3) - geta_3d_l(K1(JCEVL),I3,K2(JCEVL)))
              if (i3 == 2) then
                 ZP(J,2) = PF(K1(JCEVL),I3,K2(JCEVL)) + (PF(K1(JCEVL),I3-1,K2(JCEVL))-   &
                       PF(K1(JCEVL),I3,K2(JCEVL)))/(geta_3d(K1(JCEVL),I3-1,K2(JCEVL))    &
                        -geta_3d(K1(JCEVL),I3,K2(JCEVL))) *                              &
                          abs(geta(I3) - geta_3d(K1(JCEVL),I3,K2(JCEVL)))
              endif
             else
                ZP(J,2)= PF(K1(JCEVL),I3, K2(JCEVL))
             endif

             if (geta_3d(I1P1,I3,I2M1) < geta(I3) ) then
                ZP(J,3) = PF(I1P1,I3,I2M1) + (PF(I1P1,I3+1,I2M1)-   &
                        PF(I1P1,I3,I2M1))/(geta_3d_l(I1P1,I3+1,I2M1)&
                        -geta_3d_l(I1P1,I3,I2M1)) *                 &
                            abs(geta_l(I3) - geta_3d_l(I1P1,I3,I2M1))

             if (i3 == 2) then
                ZP(J,3) = PF(I1P1,I3,I2M1) + (PF(I1P1,I3+1,I2M1)- &
                        PF(I1P1,I3,I2M1))/(geta_3d(I1P1,I3+1,I2M1)&
                        -geta_3d(I1P1,I3,I2M1)) *                 &
                            abs(geta(I3) - geta_3d(I1P1,I3,I2M1))
             endif
             elseif  (geta_3d(I1P1,I3,I2M1) > geta(I3) ) then
                 ZP(J,3) = PF(I1P1,I3,I2M1) + (PF(I1P1,I3-1,I2M1)-   &
                        PF(I1P1,I3,I2M1))/(geta_3d_l(I1P1,I3-1,I2M1) &
                        -geta_3d_l(I1P1,I3,I2M1)) *                  &
                             abs(geta_l(I3) - geta_3d_l(I1P1,I3,I2M1))

             if (i3 == 2) then
                 ZP(J,3) = PF(I1P1,I3,I2M1) + (PF(I1P1,I3-1,I2M1)-   &
                        PF(I1P1,I3,I2M1))/(geta_3d(I1P1,I3-1,I2M1)   &
                        -geta_3d(I1P1,I3,I2M1)) *                    &
                             abs(geta(I3) - geta_3d(I1P1,I3,I2M1))
             endif
            else
                ZP(J,3)= PF(I1P1,I3,I2M1)
            endif

            if (geta_3d(I1P1,I3,K2(JCEVL)) < geta(I3) ) then
                ZP(J,4) = PF(I1P1,I3,K2(JCEVL)) + (PF(I1P1,I3+1,K2(JCEVL))-   &
                        PF(I1P1,I3,K2(JCEVL)))/(geta_3d_l(I1P1,I3+1,K2(JCEVL))&
                        -geta_3d_l(I1P1,I3,K2(JCEVL)) )   *                   &
                          abs(geta_l(I3) - geta_3d_l(I1P1,I3,K2(JCEVL)))

              if (i3 == 2) then
                ZP(J,4) = PF(I1P1,I3,K2(JCEVL)) + (PF(I1P1,I3+1,K2(JCEVL))-   &
                        PF(I1P1,I3,K2(JCEVL)))/(geta_3d(I1P1,I3+1,K2(JCEVL))  &
                        -geta_3d(I1P1,I3,K2(JCEVL)) )   *                     &
                          abs(geta(I3) - geta_3d(I1P1,I3,K2(JCEVL)))
              endif

             elseif  (geta_3d(I1P1,I3,K2(JCEVL) ) > geta(I3) ) then
                 ZP(J,4) = PF(I1P1,I3,K2(JCEVL)) + (PF(I1P1,I3-1,K2(JCEVL))-   &
                        PF(I1P1,I3,K2(JCEVL)))/(geta_3d_l(I1P1,I3-1,K2(JCEVL)) &
                        -geta_3d_l(I1P1,I3,K2(JCEVL))) *                       &
                         abs(geta_l(I3) - geta_3d_l(I1P1,I3,K2(JCEVL)))

              if (i3 == 2) then
                 ZP(J,4) = PF(I1P1,I3,K2(JCEVL)) + (PF(I1P1,I3-1,K2(JCEVL))-   &
                        PF(I1P1,I3,K2(JCEVL)))/(geta_3d(I1P1,I3-1,K2(JCEVL))   &
                        -geta_3d(I1P1,I3,K2(JCEVL))) *                         &
                         abs(geta(I3) - geta_3d(I1P1,I3,K2(JCEVL)))
              endif
             else
                ZP(J,4)= PF(I1P1,I3,K2(JCEVL))
             endif

        else
          I3= MIN(KLEV,MAX(1,I3M2+J))
          ZP(J,1)= PF(K1  (JCEVL),I3,I2M1       )
          ZP(J,2)= PF(K1  (JCEVL),I3,K2  (JCEVL))
          ZP(J,3)= PF(I1P1       ,I3,I2M1       )
          ZP(J,4)= PF(I1P1       ,I3,K2  (JCEVL))
        ENDIF
        ENDDO ! LOOP OVER LEVELS(4)

       !
       !*   3. DO CUBIC INTERPOLATION IN VERTICAL DIRECTION.
       !       -- ----- ------------- -- -------- ----------
       !
!       Z3  = A3DZ3(K3(JCEVL))                   ! = ZX(3) - ZX(2)
       Z2D3= PD3(JCEVL) * A3DRZ3(K3(JCEVL))     ! = PD3/Z3
       ZDH1= A3DDH1(K3(JCEVL))                  ! = Z3/(ZX(3)-ZX(1))
       ZDH2= A3DDH2(K3(JCEVL))                  ! = Z3/(ZX(4)-ZX(2))
       DO J=1,4   ! LOOP OVER HORIZONTAL POSITIONS
          ZA(4)= ZP(2,J)
          ZA(3)= (ZP(3,J)-ZP(1,J))*ZDH1     !DERIV AT POINT 2 X Z3
          ZD3T3= (ZP(4,J)-ZP(2,J))*ZDH2     !DERIV AT POINT 3 X Z3
          ZA(2)= 3*ZP(3,J) - ZD3T3 - 2*ZA(3) - 3*ZA(4)
          ZA(1)= ZP(3,J) - ZA(2) - ZA(3) - ZA(4)
          ZR(J)= ZA(4) + Z2D3*(ZA(3)+Z2D3*(ZA(2)+Z2D3*ZA(1)))
       ENDDO
       !
       !*   4. DO LINEAR INTERPOLATION IN HORIZONTAL DIRECTION.
       !       -- ------ ------------- -- ---------- ----------
       !
       ZD1= PD1(JCEVL)/GDLON
       Z1= ZR(1) + (ZR(3)-ZR(1))*ZD1
       Z2= ZR(2) + (ZR(4)-ZR(2))*ZD1
       PFINTER(JCEVL)= Z1 + (Z2-Z1)*PD2(JCEVL)/GDLAT(K2(JCEVL))

    ENDDO

  END SUBROUTINE ATTILA_3DINTER_2
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  INTEGER FUNCTION NRANWEI(NRP,HARVESTINIPOS,PW)

    ! PURPOSE:
    !   - GET A WEIGHTED RANDOM NUMBER FROM 1..K
    !     THE PROBABILITY OF *NRANWEI* BEING *I* IS *PW*(*I*), THEREFORE
    !     THE ARRAY *PW* SHOULD ADD UP TO (APPROX.) 1.
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98,
    !   M. Traub, P. Joeckel, MPI Mainz
    !
    ! INTERFACE:
    !   - NRANWEI WAS CALLED FROM LTINIPOS
    !   - NRANWEI IS CALLED FROM ATTILA_INIPOS
    !
    ! PARAMETERS:
    !   I NRP             RANDOM NUMBER POSITION
    !   I HARVESTINIPOS   RANDOM NUMBER ARRAY
    !   I PW              WEIGHTS
    !
    ! METHODS:
    !   - SIMPLE LOOP

    IMPLICIT NONE

    ! I/O
    INTEGER,INTENT(IN)                  :: NRP
    REAL(dp),                INTENT(IN) :: PW(:)
    REAL(dp), DIMENSION(:,:),INTENT(IN) :: HARVESTINIPOS

    ! LOCAL
    INTEGER  :: K
    INTEGER  :: J
    REAL(dp) :: ZWA,Z

    K = SIZE(PW)
    Z=HARVESTINIPOS(NRP,4)

    J=1
    ZWA=PW(1)
    DO WHILE((Z > ZWA).AND.(J < K))
       J=J+1
       ZWA=ZWA+PW(J)
    ENDDO
    NRANWEI = J

  END FUNCTION NRANWEI
! -------------------------------------------------------------------

! -------------------------------------------------------------------
  INTEGER FUNCTION GINDEX(PPOS,PGRID)

    ! PURPOSE:
    !   - DETERMINE THE MAXIMAL INDEX OF AN ASCENDING GRID FOR WHICH
    !     PPOS IS  >=  GRID VALUE.
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98
    !   M. Traub, P. Joeckel, MPI Mainz
    !     - KDIM REMOVED FROM PARAMETERS
    !     - COMMENTED BISECTION LOOP REMOVED (NO VECTORIZATION)
    !
    ! INTERFACE:
    !   GINDEX IS CALLED FROM SEVERAL ROUTINES OF THE LAGRANGIAN SCHEME
    !
    ! PARAMETERS:
    !   I PPOS  POSITION WITHIN GRID
    !   I PGRID GRID
    !
    ! METHODS:
    !       A) BISECTION. THIS LOOP NEEDS APPROX. CEIL(LOG_2(KDIM))
    !                     ITERATIONS, BUT CANNOT BE VECTORIZED.
    !       B) SIMPLE SEARCH LOOP. THIS LOOP NEEDS IN THE AVERAGE
    !          KDIM/2 ITERATIONS, BUT CAN BE VECTORIZED.
    !      ASSUMING A SPEEDUP OF 10 IN VECTOR CODE, KDIM HAS TO BE
    !      AT LEAST 160 TO JUSTIFY BISECTION. THIS IS CURRENTLY NOT
    !      THE CASE.
    !        -> BISECTION HAS BEEN REMOVED
    !
    !      RESULT: PPOS LIES IN THE INTERVAL
    !                             [ PGRID(GINDEX);PGRID(GINDEX+1) )
    !              GINDEX=0    IF PPOS< PGRID(1)
    !                    =KDIM IF PPOS>=PGRID(KDIM)
    !
    ! REFERENCES:
    !      BOOK: *W. *H. *PRESS ET AL., *NUMERICAL *RECIPES, P 90.

    IMPLICIT NONE

    ! I/O
    REAL(dp),    INTENT(IN) :: PPOS
    REAL(dp),    INTENT(IN) :: PGRID(:)  ! KDIM

    ! LOCAL
    INTEGER J, KDIM

    ! INIT
    KDIM = SIZE(PGRID)

       IF ( PPOS >= PGRID(KDIM) ) THEN
          GINDEX = KDIM
       ELSE
          DO J=1,KDIM
             IF(PPOS < PGRID(J)) EXIT
          END DO
          GINDEX=J-1
       ENDIF

  END FUNCTION GINDEX
! op_sb_20140219-
! -------------------------------------------------------------------

! -------------------------------------------------------------------
! op_sb_20140219+
  INTEGER FUNCTION GINDEX_V(PPOS,PGRID)

    ! PURPOSE:
    !   - DETERMINE THE MAXIMAL INDEX OF AN ASCENDING GRID FOR WHICH
    !     PPOS IS  >=  GRID VALUE.
    !
    ! AUTHOR(S):
    !   Ch. Reithmeier, DLR Oberpfaffenhofen,  18.5.98
    !   M. Traub, P. Joeckel, MPI Mainz
    !     - KDIM REMOVED FROM PARAMETERS
    !     - COMMENTED BISECTION LOOP REMOVED (NO VECTORIZATION)
    !
    ! INTERFACE:
    !   GINDEX IS CALLED FROM SEVERAL ROUTINES OF THE LAGRANGIAN SCHEME
    !
    ! PARAMETERS:
    !   I PPOS  POSITION WITHIN GRID
    !   I PGRID GRID
    !
    ! METHODS:
    !       A) BISECTION. THIS LOOP NEEDS APPROX. CEIL(LOG_2(KDIM))
    !                     ITERATIONS, BUT CANNOT BE VECTORIZED.
    !       B) SIMPLE SEARCH LOOP. THIS LOOP NEEDS IN THE AVERAGE
    !          KDIM/2 ITERATIONS, BUT CAN BE VECTORIZED.
    !      ASSUMING A SPEEDUP OF 10 IN VECTOR CODE, KDIM HAS TO BE
    !      AT LEAST 160 TO JUSTIFY BISECTION. THIS IS CURRENTLY NOT
    !      THE CASE.
    !        -> BISECTION HAS BEEN REMOVED
    !
    !      RESULT: PPOS LIES IN THE INTERVAL
    !                             [ PGRID(GINDEX);PGRID(GINDEX+1) )
    !              GINDEX=0    IF PPOS< PGRID(1)
    !                    =KDIM IF PPOS>=PGRID(KDIM)
    !
    ! REFERENCES:
    !      BOOK: *W. *H. *PRESS ET AL., *NUMERICAL *RECIPES, P 90.

    IMPLICIT NONE

    ! I/O
    REAL(dp),    INTENT(IN) :: PPOS
    REAL(dp),    INTENT(IN) :: PGRID(:)  ! KDIM

    ! LOCAL
    INTEGER J, KDIM, I_VERT_CORR

    ! INIT
    KDIM = SIZE(PGRID)

    IF ( .NOT. ASSOCIATED(ptemp) .AND. (I_VERT == 2) ) THEN
       I_VERT_CORR = 1
    ELSE
       I_VERT_CORR = I_VERT
    END IF

    SELECT CASE(I_VERT_CORR)

    CASE(1,3)

       IF ( PPOS >= PGRID(KDIM) ) THEN
          GINDEX_V = KDIM
       ELSE
          DO J=1,KDIM
             IF(PPOS < PGRID(J)) EXIT
          END DO
          GINDEX_V=J-1
       ENDIF

    CASE(2)

       IF ( PPOS > PGRID(1) ) THEN
          GINDEX_V = 0
       ELSE
          DO J=KDIM,1,-1
             IF(PPOS <= PGRID(J)) EXIT
          END DO
          GINDEX_V=J
       ENDIF

    END SELECT

  END FUNCTION GINDEX_V
! op_sb_20140219-
! -------------------------------------------------------------------

END MODULE MESSY_ATTILA
! ===================================================================
