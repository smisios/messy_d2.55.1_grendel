#include "messy_main_ppd_bi.inc"

! **************************************************************************
MODULE messy_mmd2way_child_si
! **************************************************************************
#ifndef I2CINC
#define _VI_SPLINE_A_
#define _VI_SPLINE_B_
#define _I2C03bDEBUG_
#else
#ifndef I2COLD
#define _VI_SPLINE_A_ 'ctspline',
#define _VI_SPLINE_B_ 'ctspline',
#define _I2C03bDEBUG_  ,izdebug
#else
#define _VI_SPLINE_A_
#define _VI_SPLINE_B_
#define _I2C03bDEBUG_
#endif
#endif

#if defined(COSMO) && defined(MESSYMMD) && defined(I2CINC)
  ! MODULE FOR exchanging data with another model in both directions
  !
  ! Author: Klaus Ketelsen,  MPICH,  Oct 2008
  ! Author: Patrick Joeckel, MPICH,  Dec 2008
  ! Author: Astrid Kerkweg,  UNI-MZ, Dec 2008-2016

  ! MESSy/BMIL
  USE messy_main_mpi_bi,        ONLY: p_parallel_io
  USE messy_main_blather_bi,    ONLY: START_MESSAGE_BI, END_MESSAGE_BI &
                                    , DEBUG_BI, ERROR_BI, INFO_BI
  USE messy_main_data_bi,       ONLY: L_IS_CHILD
  USE messy_main_qtimer_bi,     ONLY: main_qtimer_measure, MEAS_ON, MEAS_OFF

  ! MESSy/SMCL
  USE messy_main_constants_mem, ONLY: dp, STRLEN_MEDIUM, sp, STRLEN_ULONG &
                                    , iouerr

  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT &
                                    , t_chaobj_cpl
  USE messy_main_tracer,        ONLY: STRLEN_FNAME
  USE messy_main_timer_event,   ONLY: time_event, io_time_event
  USE messy_main_tools,         ONLY: PTR_4D_ARRAY, PTR_2D_ARRAY

  USE messy_main_grid,             ONLY: t_geohybgrid
  USE messy_main_grid_netcdf,      ONLY: t_ncvar
  USE messy_main_grid_trafo_scrp,  ONLY: t_scrip_data

  ! MMD2WAY
  USE messy_mmd2way_child
  USE messy_mmd2way
  ! MMD
  USE mmd_child,                  ONLY: MMD_C_GetParentType  &
                                      , MMD_ParentIsEcham

  ! INT2COSMO
  USE data_grid_in,               ONLY: pingrid, pinhgrid
  USE data_grid_lm,               ONLY: i2cgrid, i2chgrid, I2C_SD, I2C_SD_ID   &
                                      , i2cUgrid,i2chUgrid,I2C_U_SD,I2C_U_SD_ID&
                                      , i2cVgrid,i2chVgrid,I2C_V_SD,I2C_V_SD_ID

  IMPLICIT NONE

  INTRINSIC :: NULL

  PRIVATE

  ! Local Setup
  ! Lon and lat coarse grid
  REAL(KIND=DP),DIMENSION(:,:), POINTER, SAVE   :: my_lon => NULL()
  REAL(KIND=DP),DIMENSION(:,:), POINTER, SAVE   :: my_lat => NULL()

  ! VERTICAL COORDINATE HYBRID COEFFICIENTS
  REAL(KIND=DP),DIMENSION(:),ALLOCATABLE, SAVE  :: vct

  ! POINTER OF TEST ARRAY
  REAL(DP), DIMENSION(:,:,:), POINTER           :: Test_Ar => NULL()
  LOGICAL, SAVE                                 :: l_test    = .FALSE.

  ! POINTER for DATA EXCHANGE TIMING
  REAL(DP), POINTER                             :: Waittime => NULL()

  ! DEFINE STRUCTURE WHICH combines INPUT, INTERMEDIATE AND OUTPUT POINTERS
  ! TOGETHER WITH THE NAMES OF VARTABs in INT2COSMO and COSMO and the parent
  !  (i.e., ECHAM or COSMO at the moment)

  ! A) CHANNEL IDENTIFICATION OF  COUPING FIELDS
  TYPE T_C_EXCH_IO
     TYPE(t_chaobj_cpl) :: PARENT
     TYPE(t_chaobj_cpl) :: CHILD
     CHARACTER(LEN=4)   :: C_INTERPOL = ''      ! INTERPOLATION METHOD
     ! Specify target field
     LOGICAL            :: L_INITIAL  = .FALSE. ! INITIAL FIELD
     LOGICAL            :: L_BOUND    = .FALSE. ! BOUNDARY FIELD
     LOGICAL            :: L_INPUT    = .FALSE. ! INPUT FIELD
     CHARACTER(LEN=STRLEN_MEDIUM) :: C_REPR ='' ! REPRESENTATION STRING
  END TYPE T_C_EXCH_IO

  ! MAXIMAL NUMBER OF EXCHANGE FIELDS
  INTEGER, PARAMETER                            :: NMAX_EXCH = 1000
  TYPE(T_C_EXCH_IO), DIMENSION(NMAX_EXCH), SAVE :: FIELD

  ! B) COMPLETE DATA STRUCT
  TYPE T_C_COUPLE_DATA
     ! CHANNEL AND CHANNEL OBJECT NAMES IN PARENT AND CHILD
     TYPE(t_chaobj_cpl)                    :: PARENT
     TYPE(t_chaobj_cpl)                    :: CHILD
     ! ORDER OF AXES IN REPRESENTATION ('X','Y','Z','N')
     CHARACTER(LEN=4)                      :: AXIS= ''
     ! DIMENSION LENGTH
     INTEGER, DIMENSION(4)                 :: ldimlen=0
     ! INTERPOLATION METHOD (only valid for arrays not included in vartab)
     ! first three entries like in INT2COSMO
     ! 'Q'  quadratic;  'L': linear; 'M' match interpolation;
     ! 'C'  conservative (via SCRIP additional fields only)
     ! 2.CHAR if 'T' positive definiteness  is required
     ! 3.CHAR if 'T' monotonicity           is required
     ! 4.CHAR if 'V' vertical interpolation is required
     !        if 'W' vertical interpolation via NCREGRID is required (only
     !               possible with 'C' horizontal interpolation and only for
     !               additional fields
     CHARACTER(LEN=4)                        :: C_INTERPOL
     ! INPUT FIELD DELIVERED BY MMD
     REAL(DP), POINTER, DIMENSION(:,:,:,:)   :: ptr_in  => NULL()
     ! INTERMEDIATE FIELD OF INT2COSMO
     REAL(DP), POINTER, DIMENSION(:,:,:,:)   :: ptr_i2c => NULL()
     ! POINTER(S) TO COSMO/MESSy FIELD(S): DIMENSION == number of time levels
     TYPE(PTR_4D_ARRAY), DIMENSION(:), POINTER :: cosmo => NULL()
     ! POINTER TO COSMO/MESSy BOUNDARY FIELD:
     !  (DIMENSION IS ALWAYS TWO FOR THE TWO BOUNDARY LAYER TIME LEVELS)
     TYPE(PTR_4D_ARRAY), POINTER, DIMENSION(:) :: cosmo_bd => NULL()
     ! RANK OF CHILD FIELD (WITHOUT TIME LEVEL DIMENSION)
     INTEGER                                 :: rank         = 0
     ! INDICATOR, IF FIELD IS IN VARTAB
     LOGICAL                                 :: lvartab      = .FALSE.
     ! INDICATOR, IF FIELD IS IN VARTAB and TILES
     LOGICAL                                 :: ltiles       = .FALSE.
     ! NAME OF VARIABLE IN VARTAB
     CHARACTER(LEN=STRLEN_FNAME)             :: vartab_name  = ''
     ! INDEX OF FIELD in var_lm
     INTEGER                                 :: vartab_idx   = 0
     ! INITIAL FIELDS REQUIRED ?
     LOGICAL                                 :: L_INITIAL    = .FALSE.
     ! BOUNDARY FIELDS REQUIRED ?
     LOGICAL                                 :: L_BOUND      = .FALSE.
     ! INPUT FIELD ?
     LOGICAL                                 :: L_INPUT      = .FALSE.
     ! NAME OF REPRESENTATION
     CHARACTER(LEN=STRLEN_MEDIUM)            :: C_REPR       = ''
     ! Number of parent field (requested for shortcut test)
     INTEGER                                 :: scn = -99
     ! Unit of field from Parent required
     LOGICAL                                 :: L_SendUnit = .FALSE.
   END TYPE T_C_COUPLE_DATA

  TYPE (T_C_COUPLE_DATA), DIMENSION(:), ALLOCATABLE :: CplData
  ! NUMBER OF COUPLE FIELDS
  INTEGER, SAVE                                   :: NEXCH
  ! NUMBER OF COPY FIELDS
  INTEGER, SAVE                                   :: NCOPY

  ! Additional interpolation via SCRIP conservative required?
  LOGICAL, SAVE  :: l_i2cscrip = .FALSE.

  REAL(dp), DIMENSION(:,:,:,:), POINTER :: psgl4d => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: pslm4d => NULL()

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! PARENT COUPLING+

  TYPE T_P_COUPLE_DATA
     ! CHANNEL AND CHANNEL OBJECT NAMES IN Child
     TYPE(t_chaobj_cpl)                    :: name
     ! POINTER TO CHILD FIELD
     REAL(DP), POINTER, DIMENSION(:,:,:,:) :: ptr_ori  => NULL()
     ! POINTER TO INTERPOLATED FIELD
     REAL(DP), POINTER, DIMENSION(:,:,:,:) :: ptr_int  => NULL()
     ! POINTER TO horizontally INTERPOLATED FIELD
     REAL(DP), POINTER, DIMENSION(:,:,:,:) :: ptr_hint => NULL()
     ! Interpolation Method
     INTEGER                               :: interpM
     LOGICAL                               :: l_SendUnit
     INTEGER                               :: idt   = -1
     ! Representation String
     CHARACTER(LEN=STRLEN_MEDIUM)          :: repr
     ! ORDER OF AXES IN REPRESENTATION ('X','Y','Z','N')
     CHARACTER(LEN=4)                      :: AXIS= ''
     ! rank of data
     INTEGER                               :: rank = 0
     ! DIMENSION LENGTH
     INTEGER, DIMENSION(4)                 :: ldimlen=0
     ! time dependent prognostic variable?
     ! 0 = not time dependent, 2=time dependent+ horizontal field (2 space dims)
     ! 3 = time dependent, 3 space dimensions
     INTEGER                               :: itimedep = 0
     ! SCRIP DATA
     TYPE(t_scrip_data), POINTER           :: PSD  => NULL()
     ! SCRIP DATA ID
     INTEGER                               :: SD_ID
  END type T_P_COUPLE_DATA
  TYPE(T_P_COUPLE_DATA), DIMENSION(:), ALLOCATABLE :: ParData

  INTEGER, SAVE         :: PD_NUM     = 0 ! number of requested parent data set
  INTEGER, SAVE         :: PD_ADD     = 0
  INTEGER, PARAMETER    :: PD_ADD_MAX = 7

  LOGICAL, SAVE, PUBLIC :: lcpl_global_start = .FALSE.
  ! VARIABLES REQUIRED FOR GENERALISED HUMIDITY:
  ! use generalised humidity calculation
  LOGICAL, SAVE         :: lgrhum  = .FALSE.
  ! PARENT CplData indices
  INTEGER, SAVE :: idx_p_qv = -99
  INTEGER, SAVE :: idx_p_qc = -99
  INTEGER, SAVE :: idx_p_qi = -99

  INTEGER, SAVE :: idx_p_t  = -99
  INTEGER, SAVE :: idx_p_u  = -99
  INTEGER, SAVE :: idx_p_v  = -99
  INTEGER, SAVE :: idx_p_ps = -99

  ! ADDITIONAL VARIABLES ...
  !  ... required for ECHAM5 - COSMO 2-way coupling
  INTEGER, SAVE :: idx_p_pres = -99
  INTEGER, SAVE :: idx_p_fis  = -99
  INTEGER, SAVE :: idx_p_fic  = -99

  ! ... required for COSMO - COSMO 2-way coupling
  INTEGER, SAVE :: idx_p_hsurf = -99
  INTEGER, SAVE :: idx_p_w     = -99
  INTEGER, SAVE :: idx_p_pp    = -99
  INTEGER, SAVE :: idx_p_ts    = -99

  ! SAVE INDEX of parent element horiz. interpolated interpolM = 1
  INTEGER, SAVE :: idx_int1  = -99

  ! save one index for 3d parent element
  INTEGER, SAVE :: idx_p_3d  = -99

  ! CHILD INPUT DATA INDICES
  ! SAVE index of field containing the surface pressure and temperature
  ! (required for vertical interpolation of parent data)
  INTEGER, SAVE :: idx_cl_t   = -99
  INTEGER, SAVE :: idx_cl_ps  = -99
  INTEGER, SAVE :: idx_cl_fis = -99

  ! ARRAYs for required for parent coupling
  REAL(dp), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: p_fis

  ! FOR CHANNEL OBJECT DEFINITION
  ! DIMENSIONS
  INTEGER, SAVE :: DIMID_POUT_IE
  INTEGER, SAVE :: DIMID_POUT_JE
  INTEGER, SAVE :: DIMID_POUT_KMIN
  INTEGER, SAVE :: DIMID_POUT_HINT
  INTEGER, SAVE :: DIMID_POUT_KMINp1
  INTEGER, SAVE :: DIMID_POUT_HINTp1

  ! REPRESENTATIONS
  INTEGER, SAVE :: REPR_POUT_2D
  INTEGER, SAVE :: REPR_POUT_3D_MID
  INTEGER, SAVE :: REPR_POUT_3D_INT
  INTEGER, SAVE :: REPR_PHOUT_3D_MID
  INTEGER, SAVE :: REPR_PHOUT_3D_INT

  INTEGER, SAVE :: REPR_POUTACC_2D
  INTEGER, SAVE :: REPR_POUTACC_3D_MID
  INTEGER, SAVE :: REPR_PHOUTACC_3D_MID
  INTEGER, SAVE :: REPR_POUTACC_3D_INT
  INTEGER, SAVE :: REPR_PHOUTACC_3D_INT

  INTEGER, SAVE :: par_kmin
  INTEGER, SAVE :: chi_kmin
  INTEGER, SAVE :: ke_int   ! number of actual vertical layers sent to Parent
  INTEGER, SAVE :: ke_hint

  ! rotated  Lon and lat parent out-grid
  REAL(KIND=DP),DIMENSION(:),   POINTER, SAVE   :: rlons => NULL()
  REAL(KIND=DP),DIMENSION(:),   POINTER, SAVE   :: rlats => NULL()
  LOGICAL           , SAVE    :: L_gridrotParenteqChild = .FALSE.

  TYPE(t_geohybgrid), SAVE    :: cgrid   ! coupling child grid
                                         ! the horizontal domain of this grid
                                         ! is smaller than the full child grid
  TYPE(t_geohybgrid), SAVE    :: chgrid  ! pure horizontal cgrid
  TYPE(t_geohybgrid), SAVE    :: pgrid   ! parent grid
  TYPE(t_geohybgrid), SAVE    :: phgrid  ! parent grid pure horizontal

  INTEGER, SAVE  :: pgrid_id
  ! INDICES FOR REDUCED CHILD GRID
  ! length of defined region
  INTEGER, SAVE  :: iis     = 0
  INTEGER, SAVE  :: iie     = 0
  INTEGER, SAVE  :: jjs     = 0
  INTEGER, SAVE  :: jje     = 0

  LOGICAL, SAVE  :: loverlap = .TRUE.

  INTEGER, SAVE  :: conserv_idx = 0

  INTEGER, SAVE  :: i_rmy_px = 0
  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: rmy_loc

  ! control level for vertical integration (parent namelist parameter!)
  REAL(dp), SAVE :: pcontrol_fi  = 30000._dp

  ! weighting function for relaxation on parent side
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: mask      => NULL()
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: fw        => NULL()
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: fw_int    => NULL()
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: cofw_int  => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: fracs     => NULL()
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: fractions => NULL()
  REAL(DP), DIMENSION(:,:),     POINTER :: hsurf_full=> NULL()

  ! namelist parameter (parent) for weighing function:
  INTEGER  :: itype_fw = 2
  INTEGER  :: icosexp  = 14
  REAL(dp) :: damprel  = 0.2_dp

  INTEGER  :: itype_VI = 1

  ! SCALING FACTORS FOR GEOHYBGRID CALCULATIONS
  REAL(dp)  :: RCF    = 10000._dp
  REAL(dp)  :: RCF_IN = 10000._dp
  ! Avoid iteration of vertical profiles for back interpolation
  LOGICAL   :: ldiagonly = .FALSE.

  ! maximal number of iterations for adapt pressure 2 (calc_alps_2)
  INTEGER, SAVE  :: nitermax = 2

  ! correct setup
  LOGICAL :: l_shortcut = .FALSE.
  REAL(dp), DIMENSION(:), POINTER :: npsitermax  => NULL()
  REAL(dp), DIMENSION(:), POINTER :: npsiter1    => NULL()
  REAL(dp), DIMENSION(:), POINTER :: npsiter01   => NULL()
  REAL(dp), DIMENSION(:), POINTER :: npsiter001  => NULL()
  REAL(dp), DIMENSION(:), POINTER :: npsiter0001 => NULL()
  ! PARENT COUPLING-

  REAL(dp),           DIMENSION(:,:), POINTER :: boxarea_out => NULL()
  TYPE(PTR_2D_ARRAY), DIMENSION(:),   POINTER :: iterps      => NULL()
  ! PARENT COUPLING -
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! TIMER
  TYPE(io_time_event), SAVE :: &
       CPL_IOEVENT = io_time_event(1, 'steps','first',0)
  TYPE(io_time_event), SAVE :: &
       BREAK_IOEVENT !io_time_event(1,'seconds','first',0)
  TYPE(io_time_event), SAVE :: &
       READEXT_IOEVENT = io_time_event(1, 'months','none',0)
  TYPE(time_event),    SAVE :: CPL_EVENT
  TYPE(time_event),    SAVE :: BREAK_EVENT
  TYPE(time_event),    SAVE :: READEXT_EVENT

  LOGICAL                     :: l_forcevars = .FALSE.
  CHARACTER(LEN=STRLEN_ULONG) :: forcevars   = ' '
  CHARACTER(LEN=STRLEN_ULONG) :: forcefile   = ' '

  LOGICAL                     :: l_I2Cori_output = .FALSE.
  TYPE(io_time_event), SAVE   :: &
       WRITEI2C_IOEVENT = io_time_event(1, 'months','none',0)
  TYPE(time_event),    SAVE   :: WRITEI2C_EVENT

  ! make lcpltimer globally accessible
  LOGICAL, SAVE             :: lcpltimer
  LOGICAL                   :: ldev = .FALSE.
  !LOGICAL                   :: ldev = .TRUE.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! TIME MEASUREMENT: Handles:
  INTEGER, SAVE :: i_initm, i_bufC, i_bufP, i_timeloop, i_interpol
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! INTERFACE section
  INTERFACE mmd2way_child_setup
    MODULE PROCEDURE mmd2way_child_setup
  END INTERFACE mmd2way_child_setup

  INTERFACE mmd2way_child_init_memory
    MODULE PROCEDURE mmd2way_child_init_memory
  END INTERFACE mmd2way_child_init_memory

  INTERFACE mmd2way_child_init_loop
    MODULE PROCEDURE mmd2way_child_init_loop
  END INTERFACE mmd2way_child_init_loop

  INTERFACE mmd2way_child_global_end
    MODULE PROCEDURE mmd2way_child_global_end
  END INTERFACE mmd2way_child_global_end

  INTERFACE mmd2way_child_write_output
    MODULE PROCEDURE mmd2way_child_write_output
  END INTERFACE mmd2way_child_write_output

  INTERFACE mmd2way_child_write_restart
    MODULE PROCEDURE mmd2way_child_write_restart
  END INTERFACE mmd2way_child_write_restart

  INTERFACE mmd2way_child_read_restart
    MODULE PROCEDURE mmd2way_child_read_restart
  END INTERFACE mmd2way_child_read_restart

  INTERFACE mmd2way_child_free_memory
    MODULE PROCEDURE mmd2way_child_free_memory
  END INTERFACE mmd2way_child_free_memory

  PUBLIC :: mmd2way_child_setup
  PUBLIC :: mmd2way_child_init_memory
  PUBLIC :: mmd2way_child_init_loop
  PUBLIC :: mmd2way_child_global_end
  PUBLIC :: mmd2way_child_write_output
  PUBLIC :: mmd2way_child_write_restart
  PUBLIC :: mmd2way_child_read_restart
  PUBLIC :: mmd2way_child_free_memory

  ! Private INTERFACEs
  INTERFACE Setup_data_exchange_with_parent
    MODULE PROCEDURE Setup_data_exchange_with_parent
  END INTERFACE Setup_data_exchange_with_parent

  INTERFACE Define_data_arrays
    MODULE PROCEDURE Define_data_arrays
  END INTERFACE Define_data_arrays

  INTERFACE flag_fields
    MODULE PROCEDURE flag_fields
  END INTERFACE flag_fields

  INTERFACE move_initial_arrays_to_Cosmo
    MODULE PROCEDURE move_initial_arrays_to_Cosmo
  END INTERFACE move_initial_arrays_to_Cosmo

  INTERFACE move_boundary_arrays_to_Cosmo
    MODULE PROCEDURE move_boundary_arrays_to_Cosmo
  END INTERFACE move_boundary_arrays_to_Cosmo

  INTERFACE Interpol_AddiArrays
    MODULE PROCEDURE Interpol_AddiArrays
  END INTERFACE Interpol_AddiArrays

  INTERFACE interpol_coarse_OneLayer
    MODULE PROCEDURE interpol_coarse_OneLayer
  END INTERFACE interpol_coarse_OneLayer

  INTERFACE interpol_vert_AddiArray
    MODULE PROCEDURE interpol_vert_AddiArray
  END INTERFACE interpol_vert_AddiArray

  INTERFACE exchange_grids
    MODULE PROCEDURE exchange_grids
  END INTERFACE exchange_grids

  INTERFACE mmd2way_child_setup_int2cosmo
    MODULE PROCEDURE mmd2way_child_setup_int2cosmo
  END INTERFACE mmd2way_child_setup_int2cosmo

  INTERFACE mmd2way_child_set_CplData
    MODULE PROCEDURE mmd2way_child_set_CplData
  END INTERFACE mmd2way_child_set_CplData

  INTERFACE mmd2way_prepare_external_data
    MODULE PROCEDURE mmd2way_prepare_external_data
  END INTERFACE mmd2way_prepare_external_data

  INTERFACE mmd2way_child_interpolation
    MODULE PROCEDURE mmd2way_child_interpolation
  END INTERFACE mmd2way_child_interpolation

  INTERFACE get_ParentTiming
    MODULE PROCEDURE get_ParentTiming
  END INTERFACE get_ParentTiming

  INTERFACE mmd2way_child_read_nml_cpl_serv
    MODULE PROCEDURE mmd2way_child_read_nml_cpl_serv
  END INTERFACE mmd2way_child_read_nml_cpl_serv

  INTERFACE mmd2way_child_read_nml_cpl
    MODULE PROCEDURE mmd2way_child_read_nml_cpl
  END INTERFACE mmd2way_child_read_nml_cpl

  INTERFACE interpret_namelist
     MODULE PROCEDURE interpret_namelist
  END INTERFACE interpret_namelist

CONTAINS

  ! ---------------------------------------------------------------------
  SUBROUTINE mmd2way_child_setup

    ! MESSy/BMIL
    USE messy_main_mpi_bi,     ONLY: p_io, p_parallel_io, p_bcast
    USE messy_main_data_bi,    ONLY: L_IS_CHILD
    USE messy_main_timer_bi,   ONLY: p_bcast_event
    USE messy_main_qtimer_bi,  ONLY: main_qtimer_measure_init
    ! MESSy/SMCL
    USE messy_main_tools,      ONLY: find_next_free_unit
    USE messy_main_timer,      ONLY: lforcedtime
    ! MMD
    USE mmd_child,             ONLY: MMD_C_Init           &
                                   , MMD_C_GetParentType  &
                                   , MMD_ParentIsEcham
    IMPLICIT NONE

    ! Local Variables
    CHARACTER(LEN=*), PARAMETER :: substr = 'mmd2way_child_setup'
    INTEGER                     :: iou        ! LOGICAL UNIT
    INTEGER                     :: status     ! error status
    CHARACTER(LEN=1)            :: PARENT_TAG

    CALL start_message_bi(submodstr, 'READ NAMELIST', substr)

    ! SET CHILD IDENTIFIER in DATA_BI TRUE
    IF (MMD_C_GetParentType() /= -1) L_IS_CHILD = .TRUE.
    lforcedtime = .TRUE.

    IF (.NOT. L_IS_CHILD) RETURN

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL MMD2WAY CORE ROUTINE:
       CALL mmd2way_child_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast(lvertsnrgd,      p_io)

    ! READ AND INTERPRET CPL-NAMELIST
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL mmd2way_child_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('error in read CPL_CHILD namelist ',substr)
    ENDIF
    CALL p_bcast(l_forcevars,     p_io)
    CALL p_bcast(forcevars,       p_io)
    CALL p_bcast(forcefile,       p_io)
    CALL p_bcast(l_I2Cori_output, p_io)
    CALL p_bcast_event(WRITEI2C_IOEVENT, p_io)
    CALL p_bcast_event(CPL_IOEVENT,      p_io)
    CALL p_bcast_event(READEXT_IOEVENT,  p_io)

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
       IF(MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN
          PARENT_TAG='E'
       ELSE
          PARENT_TAG='C'
       ENDIF

       CALL mmd2way_child_read_nml_cpl_serv(status, iou, PARENT_TAG)
       IF (status /= 0) CALL error_bi( &
            'error in read CPL_CHILD_'//PARENT_TAG//' namelist ',substr)
    ENDIF
    CALL end_message_bi(submodstr, 'READ NAMELIST', substr)

    ! INITIALIZE CHILD IN MMD
    CALL start_message_bi(submodstr, 'SETUP MMD', substr)
    IF (ldev) THEN
       CALL MMD_C_Init(l2way=.TRUE.)
    ELSE
       CALL MMD_C_Init(l2way=.FALSE.)
    ENDIF
    CALL end_message_bi(submodstr, 'SETUP MMD', substr)

    CALL start_message_bi(submodstr, 'SETUP DATES', substr)
    CALL get_ParentTiming
    CALL end_message_bi(submodstr, 'SETUP DATES', substr)

    IF(MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN
       nitermax = nitermax
    ELSE
       ! iteration not required for COSMO - COSMO coupling
       nitermax = 0
    ENDIF

    CALL main_qtimer_measure_init(i_initm,    'MMDC_initm',    submodstr)
    CALL main_qtimer_measure_init(i_bufC,     'MMDC_bufC',     submodstr)
    CALL main_qtimer_measure_init(i_bufP,     'MMDC_bufP',     submodstr)
    CALL main_qtimer_measure_init(i_timeloop, 'MMDC_timeloop', submodstr)
    CALL main_qtimer_measure_init(i_interpol, 'MMDC_interpol', submodstr)

    RETURN

  END SUBROUTINE mmd2way_child_setup
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE mmd2way_child_init_memory

    ! MESSy/BMIL
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, SCALAR, GP_2D_HORIZONTAL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_mpi_bi,           ONLY: switch_par_utilities
    USE messy_main_timer_bi,         ONLY: timer_event_init
    USE messy_main_grid_def_bi,      ONLY: hyam, hybm, hyai, hybi
    USE messy_main_grid_def_mem_bi,  ONLY: rdheight, ke
    USE messy_main_qtimer_bi,        ONLY: main_qtimer_measure

    ! MESSy/SMCL
    USE messy_main_timer_event,   ONLY: event_delta
    USE messy_main_timer,         ONLY: lstart
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute
    USE messy_main_constants_mem, ONLY: g, rd
    ! MMD
    USE mmd_child,                ONLY: MMD_STATUS_OK               &
                                      , MMD_C_Get_ParDataArray_Name &
                                      , MMD_C_Set_DataArray_Name    &
                                      , MMD_C_Set_DataArray_EndList &
                                      , MMD_C_setInd_and_AllocMem   &
                                      , MMD_Inter_Bcast             &
                                      , MMD_Send_to_Parent
    ! INT2COSMO
    USE data_int2lm_control,      ONLY: lcomp_bound, linitial,lbd_frame &
                                      , lbd_frame_cur, i2c_dt => dt!     &
                                      !, lgrh
    USE data_grid_lm,             ONLY: ak_lm, bk_lm

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
#ifdef COSMO
    CHARACTER(LEN=*), PARAMETER :: substr = 'cosmo_child_init_memory'
#else
    CHARACTER(LEN=*), PARAMETER :: substr = 'echam_child_init_memory'
#endif
    INTEGER                     :: status, status_tmp
    INTEGER                     :: ii
    INTEGER, DIMENSION(1)       :: ikmin
    INTEGER, DIMENSION(11)      :: intbuffer
    REAL(KIND=DP), DIMENSION(1) :: p_rdheight
    REAL(KIND=DP), PARAMETER    :: t0sl=288.15_dp, dt0lp=42._dp, p0sl=1.e5_dp
    REAL(dp)                    :: p0hl
    INTEGER                     :: k

    CALL main_qtimer_measure(i_initm, MEAS_ON)

    IF (.NOT. L_IS_CHILD) RETURN
    ! INITIALIZE EVENT FOR READING OF EXTERNAL DATA
    CALL start_message_bi(submodstr, 'SET READEXT EVENT', substr)
    CALL timer_event_init(READEXT_EVENT, READEXT_IOEVENT,&
         'READEXT_EVENT', 'present')
    CALL end_message_bi(submodstr, 'SET READEXT EVENT', substr)

    ! INITIALIZE COUPLING EVENT
    CALL start_message_bi(submodstr, 'SET COUPLING EVENT', substr)
    CALL timer_event_init(CPL_EVENT, CPL_IOEVENT, 'CPL_EVENT', 'next')
    CALL end_message_bi(submodstr, 'SET COUPLING EVENT', substr)

    i2c_dt = event_delta(CPL_EVENT)
    ! INITIALIZE I2C ORIGINAL OUTPUT
    CALL start_message_bi(modstr, 'SET WRITE I2C EVENT', substr)
    CALL timer_event_init(WRITEI2C_EVENT, WRITEI2C_IOEVENT &
         , 'WRITEI2C_EVENT', 'present')
    CALL end_message_bi(modstr, 'SET WRITE I2C EVENT', substr)

    ! INITIALIZE EVENT FOR INFORMATION EXCHANGE ABOUT MODEL INTERRUPTION
    ! (must be equal to the parent time step)
    CALL start_message_bi(submodstr, 'SET TIMER BREAK EVENT', substr)
    CALL timer_event_init(BREAK_EVENT, BREAK_IOEVENT, 'BREAK_EVENT', 'next')
    CALL end_message_bi(submodstr, 'SET TIMER BREAK EVENT', substr)

    CALL start_message_bi(submodstr, 'INIT MEMORY', substr)

    CALL start_message_bi(submodstr, 'INTERPRET NAMELIST', substr)
    CALL interpret_namelist
    CALL end_message_bi(submodstr, 'INTERPRET NAMELIST', substr)

    ! DEFINE SUBMODEL-CHANNEL
    CALL new_channel(status, submodstr, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)

    IF (l_test) THEN
       ! ALLOCATE TEST ARRAY
       CALL new_channel_object(status, submodstr, 'Test_Ar', p3=Test_Ar &
            , reprid = GP_3D_MID, lrestreq=.FALSE.)
       CALL channel_halt(substr, status)
    ENDIF

    ! ALLOCATE TIME ARRAY
    CALL new_channel_object(status, submodstr, 'WaitTime', p0=WaitTime &
         , reprid = SCALAR, lrestreq=.FALSE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'WaitTime' &
         , 'long_name', c='Child waits for Parent')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'WaitTime' &
         , 'units', c='s')
    CALL channel_halt(substr, status)

    IF (ldev) THEN
       CALL new_channel_object(status, submodstr, 'fw', p4=fw &
            , reprid = GP_2D_HORIZONTAL, lrestreq=.FALSE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, submodstr, 'fw' &
            , 'long_name', c='weight function cos**2')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, submodstr, 'fw' &
            , 'units', c='-')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, submodstr, 'mask', p4=mask &
            , reprid = GP_2D_HORIZONTAL, lrestreq=.FALSE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, submodstr, 'mask' &
            , 'long_name', c='mask for nudging region')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, submodstr, 'mask' &
            , 'units', c='-')
       CALL channel_halt(substr, status)
    ENDIF ! ldev

    ALLOCATE(ak_lm(ke+1))
    ALLOCATE(bk_lm(ke+1))
    ak_lm = hyai
    bk_lm = hybi
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF (ldev) THEN
       ! PARENT COUPLING +
       CALL MMD_C_Get_ParDataArray_Name(PD_NUM)

       CALL MMD_Inter_Bcast (intbuffer, .FALSE.)
       IF (intbuffer(1) == 1 ) THEN
          lgrhum = .TRUE.
!!$          lgrh   = .TRUE.
       ELSE
          lgrhum = .FALSE.
!!$          lgrh   = .FALSE.
       END IF
       i_rmy_px = intbuffer(2)

       pcontrol_fi = REAL(intbuffer(3),dp)

       itype_fw = intbuffer(4)

       icosexp  = intbuffer(5)

       damprel  = REAL(intbuffer(6),dp) / 100._dp

       itype_VI = intbuffer(7)

       IF (intbuffer(8) == 1) lcpl_global_start = .TRUE.

       RCF      = REAL(intbuffer(9),dp)
       RCF_IN   = REAL(intbuffer(10),dp)

       IF (intbuffer(11) == 1) THEN
          ldiagonly = .TRUE.
       ELSE
          ldiagonly = .FALSE.
       END IF

       ! determine height for exchange only below damping layer
       ! calculate rdheight in pressure coordinate : p_rdheight
       p_rdheight(1) = p0sl * EXP(-t0sl/dt0lp * (1._dp - &
            SQRT(1._dp - (2._dp*dt0lp*g/(t0sl**2._dp*rd)) * rdheight)))

       chi_kmin=1
       DO k=SIZE(hyam),1,-1
          p0hl = hyam(k) + hybm(k) * p0sl
          IF(p0hl <= p_rdheight(1)) THEN
             chi_kmin = k
             EXIT
          END IF
       END DO

       CALL MMD_Send_to_Parent (p_rdheight, 1, 0, 110, status)
       IF (status /= 0) CALL error_bi('MMD_Send_to_Parent ERROR 110',substr)

       CALL MMD_Inter_Bcast (ikmin, .FALSE.)
       par_kmin = ikmin(1)

    ! PARENT COUPLING -
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ENDIF ! ldev

    DO ii=1, NEXCH
       ! CHECK IF  vertical interpolation is requested for GP_3D_MID
       IF (TRIM(CplData(ii)%C_REPR) == 'GP_3D_MID') THEN
          ! TEST if vertical interpolation is required
          ! a) usual int2lm 'V'
          IF ( (CplData(ii)%C_INTERPOL(4:4) /= 'V') .AND. &
          ! b) SCRIP / NCREGID interpolation  'W'
               (CplData(ii)%C_INTERPOL(4:4) /= 'W')) THEN
             write(iouerr,*) substr, 'NO VERTICAL INTERPOLATION REQUESTED FOR' &
               ,  TRIM(CplData(ii)%CHILD%CHA),' ',TRIM(CplData(ii)%CHILD%OBJ)&
               , 'INTERPOL: ', TRIM(CplData(ii)%C_INTERPOL)
             CALL error_bi(&
                  'no vertical interpolation requested for GP_3D_MID object' &
                  , substr)
          END IF
       END IF
    ENDDO

    ! Set Names of Data Arrays for coupling IN MMD
    status = MMD_STATUS_OK
    DO ii=1, NEXCH
       IF (TRIM(CplData(ii)%CHILD%CHA)=='mmd2way_child' .AND. &
            CplData(ii)%L_INPUT) THEN
          CplData(ii)%L_SendUnit=.TRUE.
       ELSE
          CplData(ii)%L_SendUnit=.FALSE.
       END IF

       CALL MMD_C_Set_DataArray_Name( &
            TRIM(CplData(ii)%PARENT%CHA) , TRIM(CplData(ii)%PARENT%OBJ) &
            , TRIM(CplData(ii)%CHILD%CHA), TRIM(CplData(ii)%CHILD%OBJ)  &
            , TRIM(CplData(ii)%C_REPR), CplData(ii)%L_SendUnit, status_tmp)
       IF(status_tmp /= MMD_STATUS_OK) status = status_tmp
    ENDDO

    IF(status /= MMD_STATUS_OK)   THEN
       CALL error_bi ('Error in MMD_Set_DataArray_Name',substr)
    ENDIF

    ! CLOSE LIST
    CALL MMD_C_Set_DataArray_EndList

    ! EXCHANGE OF REQUIRED LONGITUDES / LATITUDES BETWEEN PARENT AND CHILD
    CALL exchange_grids

    ! +++++++++++++++++++++++++++++++++++
    ! FOR PARENT COUPLING:
    ! this needs to take place before the dat
    IF (ldev) CALL parent_assign_ParData
    ! +++++++++++++++++++++++++++++++++++

    ! SETUP OF INT2COSMO
    CALL mmd2way_child_setup_int2cosmo
    CALL CALC_lonlat_lm
    IF (l_i2cscrip) CALL CALC_forward_weights

    IF (ldev) THEN
       ! +++++++++++++++++++++++++++++++++++
       ! FOR PARENT COUPLING:
       ! determine target parent grid
       CALL match_parent_grid

       ! determine index list
       CALL Setup_ParData_exchange
       ! initialise memory
       ! 1.) make representations for parent out-grid
       CALL parent_make_representations

       ! 2.) define data objects according to input list
       ! and fill library information list
       CALL parent_set_ParData
       ! +++++++++++++++++++++++++++++++++++
    ENDIF ! ldev

    ! Setup data exchange with parent
    CALL Setup_data_exchange_with_parent

    ! SET POINTER OF CplData STRUCT AND TEST IF CHANNEL OBJECTS ARE AVAILABLE
    CALL mmd2way_child_set_CplData

    ! PUT CplData POINTER INTO MMD
    CALL Define_data_arrays

    IF (ldev) CALL init_parent_coupling

    ! DEFINE C-core of library for Child and Parent coupling
    ! and allocate memory
    CALL MMD_C_setInd_and_AllocMem

    IF (lstart) THEN
       ! exchange, interpol and copy data
      CALL exchange_interpol_data
    ELSE
       ! SWITCH BACK TO COSMO PARALLELISATION
       CALL switch_par_utilities(2)
    ENDIF

    ! SET INT2COSMO specific switches
    ! set here as they remain constant during the integration and need not
    ! to be set again each coupling timestep.
    IF (lbd_frame) THEN
       lbd_frame_cur=.TRUE.
    ELSE
       lbd_frame_cur=.FALSE.
    ENDIF

    linitial    = .FALSE.
    lcomp_bound = .TRUE.

    CALL end_message_bi(submodstr, 'INIT MEMORY', substr)
    CALL main_qtimer_measure(i_initm, MEAS_OFF)

  END SUBROUTINE mmd2way_child_init_memory
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  SUBROUTINE init_parent_coupling

    ! MESSy/BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: DIMID_TLV
    ! MESSY/SMCL
    USE messy_main_channel,       ONLY: get_channel_object      &
                                      , get_channel_object_info &
                                      , new_channel_object      &
                                      , get_attribute
    USE messy_main_channel_repr,  ONLY: get_representation_info
    ! MMD
    USE mmd_mpi_wrapper,          ONLY: MMD_Inter_Bcast

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER     :: substr = 'init_parent_coupling'

    CHARACTER(LEN=STRLEN_CHANNEL)   :: myChannel
    CHARACTER(LEN=STRLEN_OBJECT)    :: myName

    INTEGER                         :: ii
    INTEGER                         :: status
    INTEGER                         :: dimids(4)
    INTEGER                         :: rank
    INTEGER                         :: reprid
    INTEGER                         :: i
    CHARACTER(LEN=STRLEN_MEDIUM)    :: unit

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! PARENT COUPLING +
    DO ii = 1 , PD_NUM

        MyChannel = ParData(ii)%name%CHA
        MyName    = ParData(ii)%name%OBJ
        IF (TRIM(MyName)== 'Test_Ar') CYCLE
        IF (TRIM(MyChannel) /= 'mmd2way_child' .AND. .NOT. l_shortcut) THEN
           CALL get_channel_object (status, TRIM(MyChannel), TRIM(myName)&
                , p4=ParData(ii)%ptr_ori)
           IF (status /=0) write(iouerr,*) 'MMD2WAY_CHILD ERROR ' &
                , TRIM(MyChannel),' ',TRIM(myName)
           CALL channel_halt(substr, status)

           ! determine if timedependent
           ! a) get representation id
           CALL get_channel_object_info(status, TRIM(MyChannel), TRIM(myName) &
                ,reprid = reprid)
           CALL channel_halt(substr, status)

           IF (ParData(ii)%l_SendUnit) THEN
              CALL get_attribute(status,TRIM(MyChannel), TRIM(myName) &
                   , 'units', c=unit)
              IF (status /= 0) unit = ' '
              CALL MMD_Inter_Bcast(unit, .TRUE.)
           END IF
        ELSE ! make own memory
           IF (.NOT. l_shortcut) THEN
              CALL get_representation_info(status, ParData(ii)%repr, id=reprid)
              CALL channel_halt(substr, status)

              CALL new_channel_object(status &
                   ,  TRIM(MyChannel), 'Par_'//TRIM(myName)&
                   , p4=ParData(ii)%ptr_ori, reprid=reprid)
              CALL channel_halt(substr, status)
           ELSE
              CALL get_representation_info(status, ParData(ii)%repr, id=reprid)
              CALL channel_halt(substr, status)

              CALL new_channel_object(status, submodstr, 'Par_'//TRIM(myName)&
                   , p4=ParData(ii)%ptr_ori, reprid=reprid)
              CALL channel_halt(substr, status)
           END IF
        END IF

        ! b) get dimension id info
        CALL get_representation_info(status,' ',ID=reprid, dim_ids=dimids &
             , rank=rank)
        CALL channel_halt(substr, status)

        ParData(ii)%itimedep = 0
        DO i=1,4
           IF (dimids(i) ==  DIMID_TLV) THEN
              ParData(ii)%itimedep = rank - 1 ! minus time dimension
              EXIT
           ENDIF
        END DO
     END DO

    ! PARENT COUPLING -
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END SUBROUTINE init_parent_coupling
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  SUBROUTINE mmd2way_child_init_loop

    ! MESSY
    USE messy_main_timer,         ONLY: lstart, lresume

    IMPLICIT NONE

    ! LOCAL
    ! CHARACTER(LEN=*), PARAMETER :: substr='mmd2way_child_init_loop'
    IF (.NOT. L_IS_CHILD) RETURN

    CALL main_qtimer_measure(i_timeloop, MEAS_ON)
    ! CHILD GETS DATA
    IF (.NOT. lstart) THEN
       Waittime = 0._dp
       lcouple2: IF (lcpltimer .OR. lresume) THEN
          ! EXCHANGE, INTERPOLATE and COPY DATA
          CALL exchange_interpol_data
       END IF lcouple2
    END IF
    CALL main_qtimer_measure(i_timeloop, MEAS_OFF)

    RETURN

  END SUBROUTINE mmd2way_child_init_loop
  ! ----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE mmd2way_child_global_end

    ! MESSY
    USE messy_main_qtimer_bi,     ONLY: main_qtimer_measure &
                                      , MEAS_ON, MEAS_OFF, MEAS_ADD
    USE messy_main_timer_bi,      ONLY: event_state
    USE messy_main_timer,         ONLY: next_date
    IMPLICIT NONE

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr='mmd2way_child_global_end'
    LOGICAL                      :: lbreaktimer

    IF (.NOT. L_IS_CHILD) RETURN

    ! CHILD SENDS DATA
    CALL main_qtimer_measure(i_timeloop, MEAS_ON)
    ! INQUIRE COUPLING INTERVAL
    lbreaktimer = event_state(BREAK_EVENT, next_date)
    lcpltimer   = event_state(CPL_EVENT, next_date)

    IF (ldev) THEN
       lcouple1: IF (lcpltimer) THEN
          ! interpolate regional fields to global grid
          CALL main_qtimer_measure(i_interpol, MEAS_ON)
          CALL interpol_parent_data
          CALL main_qtimer_measure(i_interpol, MEAS_OFF)

          CALL main_qtimer_measure(i_bufP, MEAS_ON)
         ! PARENT COUPLING
          CALL exchange_parent_data
          CALL main_qtimer_measure(i_bufP, MEAS_OFF)
       END IF lcouple1
    ENDIF

   ! EXCHANGE BREAK INFO in each PARENT TIME STEP
    IF (lbreaktimer) THEN
       CALL exchange_breakinfo
    ENDIF

    CALL main_qtimer_measure(i_timeloop, MEAS_ADD)
    RETURN

  END SUBROUTINE mmd2way_child_global_end
  ! ----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE mmd2way_child_write_output(flag)

    USE data_grid_in, ONLY: ie_in, je_in

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='mmd2way_child_write_output'

    INTEGER  :: ii, kx, mx
    REAL(dp) :: invfrac(1:ie_in,1:je_in)

    IF (.NOT. L_IS_CHILD) RETURN
    ! Here only output for parent coupling is converted
    IF (PD_NUM== 0) RETURN

    SELECT CASE(flag)
    CASE(-1)
       CALL debug_bi(modstr, 'CONVERT BACKWARD ', '', substr)
       DO ii = 1,PD_NUM + PD_ADD
          DO mx = 1, SIZE(ParData(ii)%ptr_int,4)
             DO kx = 1, SIZE(ParData(ii)%ptr_int,3)
                WHERE (fractions(:,:,1,1) > 0._dp)
                   ParData(ii)%ptr_int(1:ie_in,1:je_in,kx,mx) = &
                     ParData(ii)%ptr_int(1:ie_in,1:je_in,kx,mx) &
                       / fractions(1:ie_in,1:je_in,1,1)
                END WHERE
             END DO
             IF (ParData(ii)%rank > 2 .OR. ii==idx_p_ps) THEN
                DO kx = 1, SIZE(ParData(ii)%ptr_hint,3)
                   WHERE (fractions(:,:,1,1) > 0._dp)
                      ParData(ii)%ptr_hint(1:ie_in,1:je_in,kx,mx) =    &
                           ParData(ii)%ptr_hint(1:ie_in,1:je_in,kx,mx) &
                           / fractions(1:ie_in,1:je_in,1,1)
                   END WHERE
                END DO
             END IF
          END DO
       END DO
       invfrac = 0.

       IF(MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN

          WHERE (fractions(:,:,1,1) > 0._dp)
             invfrac = invfrac / fractions(:,:,1,1)
          END WHERE
          do kx = 1, nitermax
             iterps(kx)%ptr = iterps(kx)%ptr * invfrac(1:ie_in,1:je_in)
          end do
       END IF
       WHERE (fractions(:,:,1,1) > 0._dp)
          cofw_int(:,:,1,1) = cofw_int(:,:,1,1) / fractions(:,:,1,1)
       ENDWHERE

    CASE(1)
       CALL debug_bi(modstr, 'CONVERT FORWARD ', '', substr)
       DO ii = 1,PD_NUM + PD_ADD
          DO mx = 1, SIZE(ParData(ii)%ptr_int,4)
             DO kx = 1, SIZE(ParData(ii)%ptr_int,3)
                ParData(ii)%ptr_int(1:ie_in,1:je_in,kx,mx) = &
                     fractions(1:ie_in,1:je_in,1,1) &
                     * ParData(ii)%ptr_int(1:ie_in,1:je_in,kx,mx)
             END DO
             IF (ParData(ii)%rank > 2 .OR. ii==idx_p_ps) THEN
                DO kx = 1, SIZE(ParData(ii)%ptr_hint,3)
                   ParData(ii)%ptr_hint(1:ie_in,1:je_in,kx,mx) = &
                        fractions(1:ie_in,1:je_in,1,1) &
                        * ParData(ii)%ptr_hint(1:ie_in,1:je_in,kx,mx)
                END DO
             END IF
          END DO
       END DO
       IF(MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN
          do kx = 1, nitermax
             iterps(kx)%ptr = fractions(1:ie_in,1:je_in,1,1) * iterps(kx)%ptr
          END do
       ENDIF
       cofw_int = fractions * cofw_int

       END SELECT
    RETURN

  END SUBROUTINE mmd2way_child_write_output
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE mmd2way_child_write_restart(flag)

    INTEGER, INTENT(IN) :: flag

    IF (.NOT. L_IS_CHILD) RETURN
    IF (ldev) CALL  mmd2way_child_write_output(flag)

  END SUBROUTINE mmd2way_child_write_restart
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE mmd2way_child_read_restart(flag)

    INTEGER, INTENT(IN) :: flag

    IF (.NOT. L_IS_CHILD) RETURN
    IF (ldev) CALL  mmd2way_child_write_output(flag)

  END SUBROUTINE mmd2way_child_read_restart
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE mmd2way_child_free_memory

    ! INT2COSMO
    USE src_cleanup,       ONLY: org_cleanup
    ! MMD
    USE mmd_child,        ONLY: MMD_C_FreeMem
    USE mmd_test,          ONLY: MMD_testC_FreeMem

    IMPLICIT NONE

    INTRINSIC :: ALLOCATED, ASSOCIATED

    IF (.NOT. L_IS_CHILD) RETURN

    IF (ALLOCATED(p_fis))       DEALLOCATE(p_fis)
    IF (ASSOCIATED(hsurf_full)) DEALLOCATE(hsurf_full)

    CALL org_cleanup

    IF (l_test)    CALL MMD_testC_FreeMem
    CALL MMD_C_FreeMem

    DEALLOCATE (CplData)

     IF (ALLOCATED(ParData)) DEALLOCATE (ParData)

  END SUBROUTINE mmd2way_child_free_memory
  ! ----------------------------------------------------------------------

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! Private Subroutines
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE get_ParentTiming

    ! MMD
    USE mmd_child,                ONLY: MMD_Inter_Bcast, MMD_Send_to_Parent
    ! MESSy/BMIL
    USE messy_main_data_bi,       ONLY: dt, ydate_ini, ntstep            &
                                      , hstart, hstop, nstop
    USE messy_main_mpi_bi,        ONLY: p_io
    USE messy_main_timer_bi,      ONLY: messy_timer_init_manager         &
                                      , timer_message                    &
                                      , messy_timer_COSMO_reinit_time    &
                                      , p_bcast_event
    ! MESSy/SMCL
    USE messy_main_timer,         ONLY: timer_set_date, timer_get_date     &
                                      , time_span_d, delta_time, add_date  &
                                      , lstart, time_days                  &
                                      , YEAR_START, MONTH_START, DAY_START &
                                      , HOUR_START, MINUTE_START, SECOND_START
    USE messy_main_timer_event,   ONLY: convert_unit2seconds

    IMPLICIT NONE

    INTRINSIC :: NINT, REAL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='mmd2way_child_get_ParentTiming'
    !
    INTEGER, DIMENSION(9)       :: timing_array
    INTEGER                     :: status
    INTEGER                     :: syr, smo, sdy, shr, smi, sse
     ! time difference in second (resume-start date)
    REAL(dp)                    :: timediff
    INTEGER                     :: dt_sc
    INTEGER, DIMENSION(1)       :: cplseconds
    REAL(dp)                    :: rcplsecs
    TYPE(time_days)             :: test_date

    !  BROADCAST NUMBER OF COUPLE SECONDS:
    CALL convert_unit2seconds(status, rcplsecs, CPL_IOEVENT)

    CALL timer_message(status, substr)
    cplseconds(1) = NINT(rcplsecs)

    IF (p_parallel_io) &
         CALL MMD_Send_to_Parent (cplseconds, 1, 0, 61, status)

    ! Get Parent Timing data
    CALL MMD_Inter_Bcast (timing_array, .FALSE.)

   IF (lstart) THEN
       CALL timer_set_date(status, 'current', timing_array(1), timing_array(2))
       CALL timer_message(status,substr)
       ! in case of very first start resume date is equal to start date
       ! if the Parent restarts and the CHILD starts for the first time
       ! the resume date of the Parent is the start date of the Child.
       CALL timer_set_date(status, 'start',  timing_array(5), timing_array(6))
       CALL timer_message(status,substr)
       CALL timer_set_date(status, 'resume',  timing_array(5), timing_array(6))
       CALL timer_message(status,substr)
    ELSE
       CALL timer_set_date(status, 'start',   timing_array(3), timing_array(4))
       CALL timer_message(status,substr)
       ! CALCULATE CORRECT CHILD RESUME DATA
       ! 1. determine Parent resume data
       CALL timer_set_date(status, test_date, timing_array(5), timing_array(6))
       CALL timer_message(status,substr)
       ! 2. determine different in time steps between Child and parent (seconds)
       dt_sc = INT(delta_time - REAL(timing_array(9),dp)) ! must be negative !
       ! 3. subtract time step difference from Parent resume date
       CALL add_date(0, dt_sc, test_date, status)
       CALL timer_message(status,substr)
       ! 4. get date
       CALL timer_get_date(status, test_date, syr, smo, sdy, shr, smi, sse)
       CALL timer_message(status,substr)
       CALL timer_set_date(status, 'resume', syr, smo, sdy, shr, smi, sse)
       CALL timer_message(status,substr)
     ENDIF

    CALL timer_set_date(status, 'stop',    timing_array(7), timing_array(8))
    CALL timer_message(status,substr)


    WRITE ( ydate_ini(1:4) , '(I4.4)' ) YEAR_START
    WRITE ( ydate_ini(5:6) , '(I2.2)' ) MONTH_START
    WRITE ( ydate_ini(7:8) , '(I2.2)' ) DAY_START
    WRITE ( ydate_ini(9:10), '(I2.2)' ) HOUR_START

    CALL timer_get_date(status, 'stop', syr, smo, sdy, shr, smi, sse)

    IF (p_parallel_io) &
         write (*,*) 'MMD2WAY_CHILD SET STOP DATE' , syr, smo, sdy, shr, smi, sse

    CALL time_span_d( timediff , YEAR_START, MONTH_START,  DAY_START    &
         , HOUR_START, MINUTE_START, SECOND_START &
         , syr, smo,   sdy, shr,     smi, sse     )
    hstop = timediff * 24.0_dp ! stop hour
    nstop = NINT(timediff*86400.0_dp /delta_time)

    CALL messy_timer_COSMO_reinit_time

    ! SET TIMER FOR MODEL BREAK AND STOP COUPLING
    BREAK_IOEVENT = io_time_event(INT(timing_array(9)),'seconds','first', 0)
    CALL p_bcast_event(BREAK_IOEVENT, p_io)

    IF (hstart <= 0._dp) THEN
       ! INIT TIME MANAGER AT MODEL START
! op_bk_20161024+
!!$    CALL messy_timer_init_manager(INT(dt), INT(ntstep))
       CALL messy_timer_init_manager(dt, INT(ntstep))
! op_bk_20161024-
    ENDIF

    RETURN

  END SUBROUTINE get_ParentTiming
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE interpret_namelist

    ! MESSy/BMIL
    USE messy_main_mpi_bi,           ONLY: p_io, p_bcast, p_parallel_io
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_tracer_mem_bi,    ONLY: ti_gp, GPTRSTR
    USE messy_main_data_bi,          ONLY: yvarini, nyvar_i, yvarbd, nyvar_b
    ! MESSy/SMCL
    USE messy_main_tracer,        ONLY: get_tracer, full2base_sub, TR_NEXIST
    USE messy_main_tools,         ONLY: match_wild, strcrack
    USE messy_main_channel,       ONLY: get_channel_info
    ! MMD
    USE mmd_child,                ONLY: MMD_C_GetParentType  &
                                      , MMD_ParentIsEcham

    IMPLICIT NONE

    INTRINSIC :: INDEX, SIZE, TRIM

    CHARACTER(LEN=*), PARAMETER             :: substr = 'interpret_namelist'
    INTEGER,          PARAMETER             :: LMAX_EXCH = 10000
    TYPE(T_C_EXCH_IO), DIMENSION(LMAX_EXCH) :: LFIELD
    LOGICAL                                 :: l_wild_exist = .FALSE.
    INTEGER                                 :: ij, ii, ik
    INTEGER                                 :: status
    CHARACTER(LEN=STRLEN_OBJECT), DIMENSION(:), POINTER  :: objname => NULL()
    CHARACTER(LEN=STRLEN_OBJECT), DIMENSION(:), POINTER  :: string => NULL()
    INTEGER                                 :: num
    LOGICAL                                 :: lw_so = .FALSE.
    INTEGER                                 :: idt
    LOGICAL                                 :: ladd
    CHARACTER(LEN=STRLEN_MEDIUM)            :: trbasename
    CHARACTER(LEN=STRLEN_MEDIUM)            :: trsubname
    CHARACTER(LEN=STRLEN_OBJECT)            :: channelname

    IF (p_parallel_io) THEN

       ! COUNT  COUPLE FIELDS
       NEXCH = 0

       DO ii=1, NMAX_EXCH
          IF (TRIM(FIELD(ii)%PARENT%CHA) == '') CYCLE
          IF (TRIM(FIELD(ii)%CHILD%CHA) == '') CYCLE
          IF (TRIM(FIELD(ii)%PARENT%OBJ) == '') CYCLE
          IF (TRIM(FIELD(ii)%CHILD%OBJ) == '') CYCLE
          !
          IF (TRIM(FIELD(ii)%CHILD%CHA) == 'tracer_gp') THEN
             FIELD(ii)%C_REPR          = 'GP_3D_MID'
             IF (TRIM(FIELD(ii)%C_INTERPOL(4:4)) == '') &
                  FIELD(ii)%C_INTERPOL(4:4) = 'V'
          ENDIF
          l_wild_exist =  (INDEX(FIELD(ii)%CHILD%OBJ, '?') > 0)      &
               .OR. (INDEX(FIELD(ii)%CHILD%OBJ, '*') > 0)
          IF (l_wild_exist) THEN
             ! WILDCARD NOT POSSIBLE FOR THE MMD2WAY_CHILD CHANNEL ITSELF
             IF (TRIM(FIELD(ii)%CHILD%CHA) == 'mmd2way_child')       &
                  CALL error_bi(substr, &
                  'wildcard for child modules not allowed')

             CALL get_channel_info(status, TRIM(FIELD(ii)%CHILD%CHA) &
                  , ONAMES=objname)
             CALL channel_halt(substr,status)
             do_obj: DO ij = 1, SIZE(objname)
                IF ( match_wild(TRIM(FIELD(ii)%CHILD%OBJ)            &
                     ,TRIM(objname(ij)))) THEN

                   IF (TRIM(FIELD(ii)%CHILD%CHA) == 'tracer_gp') THEN
                   ! Exclude COSMO tracers (individual vartab setting)
                      IF (TRIM(ti_gp(ij)%tp%ident%submodel) == 'COSMO') CYCLE

                      ! exclude SCAV liquid tracers in case of coupling to ECHAM
                      IF (MMD_C_GetParentType() == MMD_ParentIsEcham)  THEN

                         CALL strcrack(TRIM(objname(ij)),'_',string,num)
                         IF (num == 2 .AND. &
                              TRIM(ti_gp(ij)%tp%ident%submodel) == 'scav') THEN
                            IF (TRIM(string(2)) == 'l') CYCLE
                            IF (TRIM(string(2)) == 'i') CYCLE
                         END IF
                      END IF
                   END IF

                   DO ik = 1, NEXCH
                      IF (TRIM(LFIELD(ik)%CHILD%OBJ) == TRIM(objname(ij)) &
                           .AND. TRIM(LFIELD(ik)%CHILD%CHA) == &
                           TRIM(FIELD(ii)%CHILD%CHA)) THEN
                         write (*,*) 'OBJECT ',TRIM(FIELD(ii)%CHILD%CHA),' ' &
                              ,TRIM(objname(ij)),' already exists! Cycle'
                         CYCLE do_obj
                      ENDIF
                   ENDDO

                   NEXCH   = NEXCH + 1
                   ! FIRST COPY SETTINGS
                   LFIELD(NEXCH) = FIELD(ii)
                   ! RESET WITH OBJECT SPECIFIC STUFF
                   LFIELD(NEXCH)%CHILD%OBJ = objname(ij)
                   LFIELD(NEXCH)%PARENT%OBJ = objname(ij)
                   IF (LFIELD(NEXCH)%CHILD%OBJ(1:4) == 'W_SO')  lw_so = .TRUE.
                ENDIF
             ENDDO do_obj
          ELSE
             NEXCH = NEXCH + 1
             LFIELD(NEXCH) = FIELD(ii)
             IF (LFIELD(NEXCH)%CHILD%OBJ(1:4) == 'W_SO')  lw_so = .TRUE.
          ENDIF
       ENDDO
       IF (.NOT. lw_so) CALL error_bi('w_so must be provided by parent',substr)

       ! CHECK L_INITIAL AND L_BOUND with COSMO setup (yvarini and yvarbd)
       ! ADD FIELDS CALCULATED BY INT2COSMO ADDITIONALLY TO THE EXCHANGE
       ! VARIABLES
       NCOPY = NEXCH
       ! CHECK INITIAL
       varini: DO ii = 1, nyvar_i
          IF ( TRIM(yvarini(ii)) == 'P') CYCLE
          fld1: DO ij = 1, NCOPY
             IF (TRIM(LFIELD(ij)%CHILD%OBJ) == TRIM(yvarini(ii))) THEN
                ! CHECK FOR L_INITIAL
                IF (.NOT. LFIELD(ij)%L_INITIAL) THEN
                   IF (TRIM(LFIELD(ij)%CHILD%CHA) == 'tracer_gp' ) THEN
                      ! CHANGE INITIAL CONDITION FOR COSMO TRACERS ONLY
                      CALL full2base_sub(status, TRIM(LFIELD(ij)%CHILD%OBJ)&
                           , trbasename, trsubname)
                      call get_tracer(status, GPTRSTR &
                           , trbasename, trsubname, idx=idt)
                      IF (TRIM(ti_gp(idt)%tp%ident%submodel) == 'COSMO') THEN
                         write (*,*) ' SET L_INITIAL = T FOR ' &
                              , TRIM(yvarini(ii))
                         LFIELD(ij)%L_INITIAL = .TRUE.
                      END IF
                   ELSE
                      write (*,*) ' SET L_INITIAL = T FOR ', TRIM(yvarini(ii))
                      LFIELD(ij)%L_INITIAL = .TRUE.
                   END IF
                ENDIF
                EXIT
             ENDIF
          ENDDO fld1
          IF (ij == NCOPY+1) THEN

             ladd = .FALSE.
             CALL full2base_sub(status, TRIM(yvarini(ii)) &
                  , trbasename, trsubname)
             call get_tracer(status, GPTRSTR &
                  , trbasename, trsubname, idx=idt)
             IF (status == TR_NEXIST) THEN
                ladd = .TRUE.
                channelname = 'COSMO_ORI'
             ELSE
                IF (TRIM(ti_gp(idt)%tp%ident%submodel) == 'COSMO') THEN
                   ladd =.TRUE.
                   channelname = 'tracer_gp'
                ENDIF
             END IF

             add: IF (ladd) THEN

                ! ADD FIELD TO LIST
                write (*,*) 'NEW MEMBER OF COPY LIST (INITIAL): ' &
                     , TRIM(yvarini(ii))
                NCOPY = NCOPY + 1
                LFIELD(NCOPY)%CHILD%OBJ = TRIM(yvarini(ii))
                LFIELD(NCOPY)%CHILD%CHA = TRIM(channelname)
                LFIELD(NCOPY)%L_INITIAL  = .TRUE.
                ! SPECIAL CASES:
                ! A) ROUGHNESS LENGTH z0 in INT2COSMO g*z0 in COSMO
                IF ( TRIM(yvarini(ii)) == 'Z0') THEN
                   LFIELD(NCOPY)%CHILD%OBJ = 'gZ0'
                END IF
             END IF add
            END IF
         END DO varini
       ! BOUNDARY TRANSFER
       varbd: DO ii = 1, nyvar_b
          IF ( TRIM(yvarbd(ii)) == 'P') CYCLE
          fld2: DO ij = 1, NCOPY
             IF (TRIM(LFIELD(ij)%CHILD%OBJ) == TRIM(yvarbd(ii))) THEN
                ! CHECK FOR L_BOUND
                IF (.NOT. LFIELD(ij)%L_BOUND) THEN
                   IF (TRIM(LFIELD(ij)%CHILD%CHA) == 'tracer_gp' ) THEN
                      ! CHANGE BOUNDARY CONDITION FOR COSMO TRACERS ONLY
                      CALL full2base_sub(status, TRIM(LFIELD(ij)%CHILD%OBJ)&
                           , trbasename, trsubname)
                      call get_tracer(status, GPTRSTR &
                           , trbasename, trsubname, idx=idt)
                      IF (TRIM(ti_gp(idt)%tp%ident%submodel) == 'COSMO') THEN
                         write (*,*) ' SET L_BOUND = T FOR ', TRIM(yvarbd(ii))
                         LFIELD(ij)%L_BOUND = .TRUE.
                      END IF
                   ELSE
                      write (*,*) ' SET L_BOUND = T FOR ', TRIM(yvarbd(ii))
                      LFIELD(ij)%L_BOUND = .TRUE.
                   END IF
                ENDIF
                EXIT
             ENDIF
          ENDDO fld2
          IF (ij == NCOPY+1) THEN

             ladd = .FALSE.
             CALL full2base_sub(status, TRIM(yvarbd(ii)) &
                  , trbasename, trsubname)
             call get_tracer(status, GPTRSTR &
                  , trbasename, trsubname, idx=idt)
             IF (status == TR_NEXIST) THEN
                ladd = .TRUE.
                channelname = 'COSMO_ORI'
             ELSE
                IF (TRIM(ti_gp(idt)%tp%ident%submodel) == 'COSMO') THEN
                   ladd =.TRUE.
                   channelname = 'tracer_gp'
                ENDIF
             END IF

             add2: IF (ladd) THEN
                ! ADD FIELD TO LIST
                write (*,*) 'NEW MEMBER OF COPY LIST (BOUND): ' &
                     ,  TRIM(yvarbd(ii))
                NCOPY = NCOPY + 1
                LFIELD(NCOPY)%CHILD%OBJ = TRIM(yvarbd(ii))
                LFIELD(NCOPY)%CHILD%CHA = TRIM(channelname)
                LFIELD(NCOPY)%L_BOUND  = .TRUE.
             END IF add2
          ENDIF
       END DO varbd

       ! ALLOCATE CplData STRUCT TO ACTUAL NUMBER OF COUPLING FIELDS
       ALLOCATE (CplData(NCOPY))

       WRITE(*,*) '.........................................................'
       WRITE(*,*) ' NUMBER OF COUPLING-FIELDS: ',NEXCH
       WRITE(*,*) '.........................................................'
       WRITE(*,*) '.........................................................'
       WRITE(*,*) ' OVERALL NUMBER OF FIELDS TO COUPLE AND COPY: ',NCOPY
       WRITE(*,*) '.........................................................'

       DO ii=1, NCOPY
          CplData(ii)%PARENT%CHA = LFIELD(ii)%PARENT%CHA
          CplData(ii)%PARENT%OBJ = LFIELD(ii)%PARENT%OBJ
          CplData(ii)%CHILD%CHA = LFIELD(ii)%CHILD%CHA
          CplData(ii)%CHILD%OBJ = LFIELD(ii)%CHILD%OBJ
          CplData(ii)%C_INTERPOL = LFIELD(ii)%C_INTERPOL
          CplData(ii)%L_BOUND    = LFIELD(ii)%L_BOUND
          CplData(ii)%L_INITIAL  = LFIELD(ii)%L_INITIAL
          CplData(ii)%L_INPUT    = LFIELD(ii)%L_INPUT
          CplData(ii)%C_REPR     = LFIELD(ii)%C_REPR
          WRITE(*,*) '.........................................................'
          WRITE(*,*) ' NUMBER OF COPY-FIELD: ',ii
          WRITE(*,*) ' CHANNEL - OBJECT PARENT = ' &
               , TRIM(CplData(ii)%PARENT%CHA),'  -  ' &
               , TRIM(CplData(ii)%PARENT%OBJ)
          WRITE(*,*) ' CHANNEL - OBJECT CHILD = ', &
               TRIM(CplData(ii)%CHILD%CHA), '  -  ' &
               ,TRIM(CplData(ii)%CHILD%OBJ)
          WRITE(*,*) ' INTERPOLATION METHOD   = ',  CplData(ii)%C_INTERPOL
          WRITE(*,*) ' CALC INITIAL / BOUND / INPUT   = ' &
               ,  CplData(ii)%L_INITIAL,'/', CplData(ii)%L_BOUND,'/' &
               ,  CplData(ii)%L_INPUT
          WRITE(*,*) ' REPRESENTATION         = ', CplData(ii)%C_REPR
          WRITE(*,*) '.........................................................'

       ENDDO
    ENDIF

    ! DISTRIBUTE VALUES
    CALL p_bcast(NEXCH, p_io)
    CALL p_bcast(NCOPY, p_io)
    IF (.NOT. p_parallel_io) ALLOCATE (CplData(NCOPY))

    DO  ii=1,NCOPY
       CALL p_bcast(CplData(ii)%PARENT%CHA, p_io)
       CALL p_bcast(CplData(ii)%PARENT%OBJ, p_io)
       IF (TRIM(CplData(ii)%PARENT%OBJ)=='Test_Ar') l_test = .TRUE.
       CALL p_bcast(CplData(ii)%CHILD%CHA, p_io)
       CALL p_bcast(CplData(ii)%CHILD%OBJ, p_io)
       CALL p_bcast(CplData(ii)%C_INTERPOL, p_io)
       CALL p_bcast(CplData(ii)%L_BOUND,    p_io)
       CALL p_bcast(CplData(ii)%L_INITIAL,  p_io)
       CALL p_bcast(CplData(ii)%L_INPUT,    p_io)
       CALL p_bcast(CplData(ii)%C_REPR,     p_io)

       IF (CplData(ii)%L_INPUT .AND. &
            (CplData(ii)%L_INITIAL .OR. CplData(ii)%L_BOUND)) THEN
          CALL error_bi('input and initial/boundary field requested', substr)
       ENDIF
    ENDDO

    DO  ii=1,NCOPY
       IF (TRIM(CplData(ii)%CHILD%OBJ) == 'PS' .AND. &
            TRIM(CplData(ii)%CHILD%CHA) == 'COSMO_ORI') THEN
          idx_cl_ps = ii
       END IF
       IF (TRIM(CplData(ii)%CHILD%OBJ) == 'T' .AND. &
            TRIM(CplData(ii)%CHILD%CHA) == 'COSMO_ORI') THEN
          idx_cl_t = ii
       END IF
       IF (TRIM(CplData(ii)%CHILD%OBJ) == 'FIS' .AND. &
            TRIM(CplData(ii)%CHILD%CHA) == '#XXX') THEN
          idx_cl_fis = ii
       ENDIF
       IF (TRIM(CplData(ii)%CHILD%OBJ) == 'HSURF' .AND. idx_cl_fis == -99) THEN
          idx_cl_fis = ii
       END IF
    END DO
    IF (idx_cl_ps == -99 .OR. idx_cl_t == -99 .OR. idx_cl_fis == -99) THEN
       CALL info_bi(&
            'surface pressure, geopotential and temperature','')
       CALL info_bi(&
            'are required for parent coupling!', '')
       CALL error_bi('Define surf. press, geopotential and temperature' &
            , substr)
    ENDIF

    ! CHECK if conservative interpolation via SCRIP is required
    DO ii= 1, NCOPY
       IF (CplData(ii)%C_INTERPOL(1:1) == 'C' .OR. &
            CplData(ii)%C_INTERPOL(4:4) == 'W') THEN
          l_i2cscrip = .TRUE.
          EXIT
       ENDIF
    END DO
    IF (p_parallel_io) THEN

       IF (l_i2cscrip) THEN
          WRITE(*,*) ' INTEPOLATION VIA SCRIP REQUIRED FOR I2C'
       ELSE
          WRITE(*,*) ' NO CONSERVATIVE VARIABLE VIA SCRIP INTERPOLATION REQUIRED FOR I2C'
       ENDIF
    END IF

  END SUBROUTINE interpret_namelist
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  SUBROUTINE exchange_grids

    ! INT2COSMO
    USE data_grid_in,             ONLY: pollat_in,pollon_in, polgam_in       &
                                      , dlat_in, dlon_in                     &
                                      , startlat_in_tot,startlon_in_tot      &
                                      , ie_in_tot,je_in_tot, ke_in_tot       &
                                      , ke1in, ke_in                         &
                                      , endlat_in_tot, endlon_in_tot         &
                                      , latitudes_in,longitudes_in           &
                                      , ak_in,   bk_in, akh_in, bkh_in       &
                                      , dak_in, dbk_in                       &
                                      , ke_soil_in, czmls_in
    USE vgrid_refatm_utils,       ONLY: vcoord_in, refatm_in, svc1_in, svc2_in

    USE data_int2lm_control,      ONLY: itype_t_cl, itype_w_so_rel &
                                      , yinput_model, l_bicub_spl
    ! MMD
    USE mmd_child,               ONLY:  MMD_Send_to_Parent   &
                                      , MMD_Recv_from_Parent &
                                      , MMD_C_GetParentType  &
                                      , MMD_ParentIsEcham
    ! MESSy/BMIL
    USE messy_main_grid_def_bi,     ONLY: rlat, rlon
    USE messy_main_grid_def_mem_bi, ONLY: ie_tot, je_tot, ie, je   &
                                        , pollon, pollat, polgam
    USE messy_main_data_bi,         ONLY: rmy
    USE messy_main_mpi_bi,          ONLY: p_bcast, num_compute, gather_field   &
                                        , p_parallel_io, p_io, distribute_field
    ! SMCL
    USE messy_main_tools,         ONLY: calc_hybrid_coeff
    USE messy_main_constants_mem, ONLY: pi !,iouerr

    IMPLICIT NONE

    INTRINSIC :: COS, SIZE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER               :: substr = 'child_exchange_grids'
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: rlon_glob,rlat_glob
    REAL(KIND=DP), DIMENSION(17)               :: Define_coarse_grid_real
    INTEGER, DIMENSION(8)                      :: Define_coarse_grid_int
    INTEGER, DIMENSION(3)                      :: val
    INTEGER                                    :: status
    INTEGER(I4)                                :: stat_i4
    INTEGER                                    :: i, k
    INTEGER                                    :: j, is_p3, js_p3, ie_p3, je_p3
    LOGICAL                                    :: l_ijf
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: rmy_glob
    INTEGER                                    :: is_w, ie_w, js_w, je_w
    INTEGER                                    :: iwlen, jwlen
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: fw_glob ! weighting function
    REAL(dp)                                   :: cosi, cosj, x ,y
    ! length of relaxation zone (number of COSMO gridboxes)
    REAL(dp)                                   :: damplen(2)
    INTEGER                                    :: idamplen(2)

    !******************************************************************
    ! PARENT COUPLING +
    !******************************************************************
    IF (ldev) THEN

       ALLOCATE(rmy_loc(ie,je))

       IF (num_compute > 1) THEN
          ALLOCATE(rmy_glob(ie_tot,je_tot))
          ALLOCATE(fw_glob(ie_tot,je_tot))
          CALL gather_field (rmy(:,:,1),ie,je, rmy_glob,ie_tot,je_tot &
               , p_io, stat_i4)
          IF (stat_i4 /= 0) CALL error_bi('ERROR: gather rmy', substr)
          pio: IF (p_parallel_io) THEN
             l_ijf = .FALSE.
             DO i=1,ie_tot
                DO j=1,je_tot
                   IF(rmy_glob(i,j) <= 0._dp) THEN
                      is_p3=i
                      js_p3=j
                      l_ijf = .TRUE.
                      exit
                   END IF
                END DO
                IF(l_ijf) exit
             END DO
             l_ijf = .FALSE.
             DO i=ie_tot,1,-1
                DO j=je_tot,1,-1
                   IF(rmy_glob(i,j) <= 0._dp) THEN
                      l_ijf = .TRUE.
                      ie_p3=i
                      je_p3=j
                      exit
                   END IF
                END DO
                IF(l_ijf) exit
             END DO
             IF (i_rmy_px >= 0) THEN
                DO i=is_p3,ie_p3
                   DO j=js_p3,js_p3+(i_rmy_px-1)
                      rmy_glob(i,j) = 0.0001_dp
                   END DO
                   DO j=je_p3-(i_rmy_px-1),je_p3
                      rmy_glob(i,j) = 0.0001_dp
                   END DO
                END DO
                DO j=js_p3+i_rmy_px,je_p3-i_rmy_px
                   DO i=is_p3,is_p3+(i_rmy_px-1)
                      rmy_glob(i,j) = 0.0001_dp
                   END DO
                   DO i=ie_p3-(i_rmy_px-1),ie_p3
                      rmy_glob(i,j) = 0.0001_dp
                   END DO
                END DO
                ! set indices for weight field
                is_w = is_p3 + i_rmy_px
                ie_w = ie_p3 - i_rmy_px
                js_w = js_p3 + i_rmy_px
                je_w = je_p3 - i_rmy_px
             ELSE
                is_w = is_p3
                ie_w = ie_p3
                js_w = js_p3
                je_w = je_p3
             ENDIF
             ! calculate length of coupled region
             ! iwlen = (ie_tot - is_w + 1) - (ie_tot - ie_w)
             iwlen = ie_w - is_w + 1
             jwlen = je_w - js_w + 1
             idamplen(1) = INT(iwlen * damprel)
             idamplen(2) = INT(jwlen * damprel)
             damplen = REAL(idamplen,dp)

             fw_glob = 0._dp
             IF (itype_fw == 0) THEN
                fw_glob = 1._dp
             ELSE IF (itype_fw == 1) THEN
                DO i=1, ie_tot
                   x = pi * i / ( ie_tot+1)
                   DO j = 1, je_tot
                      y = pi * j / ( je_tot+1)
                      fw_glob(i,j) = 1._dp - (COS(x)**icosexp + COS(y)**icosexp)
                   END DO
                END DO
                WHERE (fw_glob(:,:) < 0._dp)
                   fw_glob(:,:) = 0._dp
                END WHERE
             ELSE IF (itype_fw == 2) THEN
                fw_glob(is_w:ie_w,js_w:je_w) = 1._dp
                ! left border
                DO i = 1, idamplen(1)

                   cosi = 0.5 - COS(REAL(i,dp)/damplen(1)*pi)/2.
                   ! left / right lower / upper corner
                   DO j = 1, idamplen(2)
                      cosj = 0.5 -COS(REAL(j,dp)/damplen(2)*pi)/2.

                     ! left lower corner
                     fw_glob(is_w + i-1,js_w+j-1) = cosi * cosj
                     ! left upper corner
                     fw_glob(is_w + i-1,je_w-j+1) = fw_glob(is_w + i-1,js_w+j-1)
                     ! right upper corner
                     fw_glob(ie_w - i+1,je_w-j+1) = fw_glob(is_w + i-1,js_w+j-1)
                     ! right lower corner
                     fw_glob(ie_w - i+1,js_w+j-1) = fw_glob(is_w + i-1,js_w+j-1)
                   END DO
                   !
                   cosj = 1._dp
                   ! left border
                   fw_glob(is_w + i-1,js_w+idamplen(2):je_w-idamplen(2)+1) &
                        = cosi*cosj
                   ! right border
                   fw_glob(ie_w - i +1,js_w+idamplen(2):je_w-idamplen(2)+1) &
                        = cosi * cosj
                END DO
                ! lower border
                cosi = 1._dp
                DO j = 1, idamplen(2)
                   cosj = 0.5-COS(REAL(j,dp)/damplen(2)*pi)/2.
                   fw_glob(is_w + idamplen(1): ie_w -idamplen(1)+1, js_w+j-1) &
                        = cosi * cosj
                   fw_glob(is_w + idamplen(1): ie_w -idamplen(1)+1, je_w-j+1) &
                        = cosi * cosj
                END DO
             END IF
          END IF pio

          CALL distribute_field( rmy_glob, ie_tot, je_tot, rmy_loc &
               , ie, je, p_io, stat_I4)
          IF (stat_i4 /= 0) CALL error_bi('ERROR: distribute rmy', substr)

          CALL distribute_field( fw_glob, ie_tot, je_tot, fw(:,:,1,1) &
               , ie, je, p_io, stat_I4)
          IF (stat_i4 /= 0) CALL error_bi('ERROR: distribute weight function' &
               , substr)

          DEALLOCATE(rmy_glob)
          DEALLOCATE(fw_glob)

          WHERE (fw > 0.9999)
             mask = 1.0
          ENDWHERE
       ENDIF
    !******************************************************************
    ! PARENT COUPLING -
    !******************************************************************
    ENDIF ! ldev

    ALLOCATE(rlat_glob(ie_tot,je_tot))
    ALLOCATE(rlon_glob(ie_tot,je_tot))

    IF (num_compute > 1) THEN
      CALL gather_field (rlat,ie,je, rlat_glob,ie_tot,je_tot, -1_I4, stat_i4)
      IF (stat_i4 /= 0) CALL error_bi('gather_field ERROR 01', substr)
      CALL gather_field (rlon,ie,je, rlon_glob,ie_tot,je_tot, -1_I4, stat_i4)
      IF (stat_i4 /= 0) CALL error_bi('gather_field ERROR 02', substr)
    ELSE
      rlat_glob(:,:) = rlat(:,:)
      rlon_glob(:,:) = rlon(:,:)
    ENDIF

    IF (p_parallel_io) THEN
       ! Send Child grid to parent
       val(1) = ie_tot
       val(2) = je_tot
       ! Send also number of coupling fields
       val(3) = NEXCH
       CALL MMD_Send_to_Parent (val, size(val), 0, 123, status)
       IF (status /= 0) CALL error_bi('MMD_Send_to_Parent ERROR 01',substr)
       CALL MMD_Send_to_Parent (rlat_glob, INT(ie_tot*je_tot), 0, 11, status)
       IF (status /= 0) CALL error_bi('MMD_Send_to_Parent ERROR 02',substr)
       CALL MMD_Send_to_Parent (rlon_glob, INT(ie_tot*je_tot), 0, 12, status)
       IF (status /= 0) CALL error_bi('MMD_Send_to_Parent ERROR 03',substr)

       ! Receiver Coarse grid information
       CALL MMD_Recv_from_Parent (Define_coarse_grid_real, 17, 0, 21,  status)
       IF (status /= 0) CALL error_bi('MMD_Recv_from_Parent ERROR 01', substr)
       CALL MMD_Recv_from_Parent (Define_coarse_grid_int,   8, 0, 22,  status)
       IF (status /= 0) CALL error_bi('MMD_Recv_from_Parent ERROR 02', substr)

!!$       write(iouerr,*) 'Coarse grid from Parent '
!!$       write(iouerr,*) 'startlat_in_tot  = ',Define_coarse_grid_real(1)
!!$       write(iouerr,*) 'startlon_in_tot  = ',Define_coarse_grid_real(2)
!!$       write(iouerr,*) 'endlat_in_tot    = ',Define_coarse_grid_real(8)
!!$       write(iouerr,*) 'endlon_in_tot    = ',Define_coarse_grid_real(9)
!!$       write(iouerr,*) 'pollat_in        = ',Define_coarse_grid_real(3)
!!$       write(iouerr,*) 'pollon_in        = ',Define_coarse_grid_real(4)
!!$       write(iouerr,*) 'polgam_in        = ',Define_coarse_grid_real(5)
!!$       write(iouerr,*) 'dlat_in          = ',Define_coarse_grid_real(6)
!!$       write(iouerr,*) 'dlon_in          = ',Define_coarse_grid_real(7)
!!$       IF (llm2lm) THEN
!!$          write(iouerr,*) 'vcflat_in      = ',Define_coarse_grid_real(10)
!!$          write(iouerr,*) 'p0sl_in        = ',Define_coarse_grid_real(11)
!!$          write(iouerr,*) 't0sl_in        = ',Define_coarse_grid_real(12)
!!$          write(iouerr,*) 'dt0lp_in       = ',Define_coarse_grid_real(13)
!!$          write(iouerr,*) 'delta_t_in     = ',Define_coarse_grid_real(14)
!!$          write(iouerr,*) 'h_scal_in      = ',Define_coarse_grid_real(15)
!!$          write(iouerr,*) 'svc1_in        = ',Define_coarse_grid_real(16)
!!$          write(iouerr,*) 'svc2_in        = ',Define_coarse_grid_real(17)
!!$       ENDIF
!!$       write(iouerr,*) 'ie_coarse     = ',Define_coarse_grid_int(1)
!!$       write(iouerr,*) 'je_coarse     = ',Define_coarse_grid_int(2)
!!$       write(iouerr,*) 'ke_coarse     = ',Define_coarse_grid_int(3)
!!$       IF (llm2lm) THEN
!!$          write(iouerr,*) 'ivctype     = ',Define_coarse_grid_int(4)
!!$          write(iouerr,*) 'irefatm     = ',Define_coarse_grid_int(5)
!!$       ENDIF
!!$       write(iouerr,*) 'ke_soil_coarse = ',Define_coarse_grid_int(6)
!!$       write(iouerr,*) 'itype_w_so_rel = ',Define_coarse_grid_int(7)
!!$       write(iouerr,*) 'itype_t_cl     = ',Define_coarse_grid_int(8)
    ENDIF

    ! clean
    DEALLOCATE(rlat_glob)
    DEALLOCATE(rlon_glob)

    ! Decompose coarse grid
    CALL p_bcast (Define_coarse_grid_real, 0_I4)
    CALL p_bcast (Define_coarse_grid_int , 0_I4)

    startlat_in_tot   = Define_coarse_grid_real(1)
    startlon_in_tot   = Define_coarse_grid_real(2)
    pollat_in         = Define_coarse_grid_real(3)
    pollon_in         = Define_coarse_grid_real(4)
    polgam_in         = Define_coarse_grid_real(5)
    dlat_in           = Define_coarse_grid_real(6)
    dlon_in           = Define_coarse_grid_real(7)
    endlat_in_tot     = Define_coarse_grid_real(8)
    endlon_in_tot     = Define_coarse_grid_real(9)
    vcoord_in%vcflat  = Define_coarse_grid_real(10)
    refatm_in%p0sl    = Define_coarse_grid_real(11)
    refatm_in%t0sl    = Define_coarse_grid_real(12)
    refatm_in%dt0lp   = Define_coarse_grid_real(13)
    refatm_in%delta_t = Define_coarse_grid_real(14)
    refatm_in%h_scal  = Define_coarse_grid_real(15)
    svc1_in           = Define_coarse_grid_real(16)
    svc2_in           = Define_coarse_grid_real(17)
    ie_in_tot         = Define_coarse_grid_int(1)
    je_in_tot         = Define_coarse_grid_int(2)
    ke_in_tot         = Define_coarse_grid_int(3)
    vcoord_in%ivctype = Define_coarse_grid_int(4)
    refatm_in%irefatm = Define_coarse_grid_int(5)
    ke_soil_in        = Define_coarse_grid_int(6)
    itype_w_so_rel    = Define_coarse_grid_int(7)
    itype_t_cl        = Define_coarse_grid_int(8)
    ke_in             = ke_in_tot        ! Has to be set here to get correct
    ke1in             = ke_in_tot+1      ! memory allocation

    IF(MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN
       yinput_model = 'CM'
       l_bicub_spl  = .TRUE.
    ELSE
       yinput_model = 'COSMO'
       l_bicub_spl  = .FALSE.
    ENDIF

    ! get vertical coordinate and soil depth information
    IF (vcoord_in%ivctype == 1) THEN
       ALLOCATE (vcoord_in%sigm_coord(ke_in_tot+1))
       vcoord_in%sigm_coord(:)=0.0_dp
    ELSE
       ALLOCATE (vcoord_in%vert_coord(ke_in_tot+1))
       vcoord_in%vert_coord(:)=0.0_dp
    ENDIF
    ALLOCATE(czmls_in(ke_soil_in+1)); czmls_in= 0.0_dp

    IF (p_parallel_io) THEN
       IF(MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN
          ALLOCATE (vct(2*(ke_in_tot+1))); vct=0.0
          CALL MMD_Recv_from_Parent (vct, INT(2*(ke_in_tot+1)), 0, 23, status)
       ELSE
          IF (vcoord_in%ivctype == 1) THEN
             CALL MMD_Recv_from_Parent ( &
                  vcoord_in%sigm_coord, INT(ke_in_tot+1), 0, 23, status)
          ELSE
             CALL MMD_Recv_from_Parent ( &
                  vcoord_in%vert_coord, INT(ke_in_tot+1), 0, 23, status)
          ENDIF
       ENDIF
       CALL MMD_Recv_from_Parent (czmls_in, INT(ke_soil_in+1), 0, 26, status)
    ENDIF
    IF(MMD_C_GetParentType() /= MMD_ParentIsEcham)   THEN
       IF (vcoord_in%ivctype == 1) THEN
          CALL p_bcast (vcoord_in%sigm_coord , 0)
       ELSE
          CALL p_bcast (vcoord_in%vert_coord , 0)
       ENDIF
    ENDIF
    CALL p_bcast (czmls_in  , 0)

    ! Get Parent latitudes and longitudes on coarse grid
    ALLOCATE (latitudes_in(je_in_tot))
    ALLOCATE (longitudes_in(ie_in_tot))

    IF (p_parallel_io) THEN
      CALL MMD_Recv_from_Parent (latitudes_in, INT(je_in_tot), 0, 24, status)
      CALL MMD_Recv_from_Parent (longitudes_in,INT(ie_in_tot), 0, 25, status)
    ENDIF

    CALL p_bcast (latitudes_in,  0)
    CALL p_bcast (longitudes_in, 0)

    ALLOCATE ( ak_in(ke_in+1), bk_in(ke_in+1)  &
             , akh_in(ke_in),  bkh_in(ke_in)   &
             , dak_in(ke_in),  dbk_in(ke_in))
    ak_in       = 0.0_dp
    bk_in       = 0.0_dp
    akh_in      = 0.0_dp
    bkh_in      = 0.0_dp
    dak_in      = 0.0_dp
    dbk_in      = 0.0_dp

    ! Distribute ak_in and bk_in, compute akh_in ans bkh_in
    IF(MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN
       IF (p_parallel_io)  THEN
          ak_in = vct(1:ke_in_tot+1)
          bk_in = vct(ke_in_tot+2:2*(ke_in_tot+1))
          DEALLOCATE(vct)
       ENDIF
       CALL p_bcast (ak_in , 0)
       CALL p_bcast (bk_in , 0)

       CALL compute_akh_bkh
    ELSE
       IF (vcoord_in%ivctype == 1) THEN
          CALL calc_hybrid_coeff(ke_in, vcoord_in%sigm_coord, vcoord_in%vcflat &
            , refatm_in%irefatm, vcoord_in%ivctype, akh_in, ak_in, bkh_in, bk_in)
       ELSE
          CALL calc_hybrid_coeff(ke_in, vcoord_in%vert_coord, vcoord_in%vcflat &
            , refatm_in%irefatm, vcoord_in%ivctype, akh_in, ak_in, bkh_in, bk_in)
       END IF
       DO k=1, ke_in
          dak_in(k) =  ak_in(k+1) - ak_in(k)
          dbk_in(k) =  bk_in(k+1) - bk_in(k)
       END DO
   ENDIF

   IF (MMD_C_GetParentType() == MMD_ParentIsEcham) THEN
      L_gridrotParenteqChild = .FALSE.
   ELSE

      IF (pollon_in == pollon .AND. pollat_in == pollat .AND. &
           polgam_in == polgam)  THEN
         L_gridrotParenteqChild = .TRUE.
      ELSE  IF  (((pollon_in == 180._dp .AND. polgam_in == 0._dp &
           .AND. pollon == 0._dp .AND. polgam == 180._dp) .OR.  &
           (pollon_in == 0._dp .AND. polgam_in == 180._dp &
           .AND. pollon == 180._dp .AND. polgam == 0._dp)) .AND. &
           pollat_in == pollat)  THEN
         L_gridrotParenteqChild = .TRUE.
      ELSE
         L_gridrotParenteqChild = .FALSE.
      ENDIF
      ! L_gridrotParenteqChild = .FALSE.
   END IF

   IF (p_parallel_io) &
        write (*,*) 'L_gridrotParenteqChild = ',L_gridrotParenteqChild &
        , 'POL LON: ', pollon_in, pollon, 'POL LAT: ',pollat_in, pollat &
        ,  'POL GAMMA', polgam_in , polgam

  END SUBROUTINE exchange_grids
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE mmd2way_child_setup_int2cosmo

    ! MESSy/BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy/SMCL
    USE messy_main_channel,    ONLY: get_channel_object

    ! INT2COSMO
    USE data_int2lm_control,   ONLY: llm2lm, lcm2lm, lgme2lm, lec2lm, lhm2lm &
                                   , lum2lm, lcomp_bound, linitial,lbd_frame &
                                   , lhir2lm, lhm2lm, lgfs2lm, lgsm2lm       &
                                   , lbd_frame_cur, pcontrol_fi
    USE data_grid_in,          ONLY: slatitudes_in, slongitudes_in           &
                                   , ie_in_tot, je_in_tot                    &
                                   , dlat_in, dlon_in, ie_in, je_in          &
                                   , latitudes_in, longitudes_in             &
                                   , lushift_in, lvshift_in, east_add_in     &
                                   , west_add_in, south_add_in, north_add_in
    USE data_grid_in,          ONLY:                                         &
                                   ! Input data has hybrid height coordinates
                                     lcm_hgt_coor                            &
                                   ! Input data has pressure coordinates
                                   , lcm_pres_coor

    ! MMD
    USE mmd_child,            ONLY: MMD_C_GetParentType &
                                   , MMD_ParentIsEcham
    USE mmd_test,              ONLY: MMD_testC_Setup

    ! MESSY/SMCL
    USE messy_main_timer,      ONLY: lstart
    USE messy_main_constants_mem, ONLY: iouerr

    IMPLICIT NONE

    EXTERNAL  :: setup_int2lm
    INTRINSIC :: INT, TRIM

    ! LOCAL
    INTEGER                      :: status
    CHARACTER(LEN=80)            :: yerror
    CHARACTER(LEN=*), PARAMETER  :: substr='mmd2way_child_setup_int2cosmo'
    INTEGER                      :: i

    ! Intialize int2cosmo
    lgme2lm = .FALSE.
    lec2lm  = .FALSE.
    lhm2lm  = .FALSE.
    lum2lm  = .FALSE.
    lgfs2lm = .FALSE.
    lgsm2lm = .FALSE.
    lhir2lm = .FALSE.
    lhm2lm  = .FALSE.

    IF(MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN
      IF (p_parallel_io) write(*,*) 'Parent is ECHAM '
      llm2lm = .FALSE.
      lcm2lm = .TRUE.
      lushift_in(:) =.FALSE.
      lvshift_in(:) =.FALSE.
      lcm_hgt_coor  = .FALSE.
      lcm_pres_coor = .FALSE.
    ELSE
      IF (p_parallel_io) write(*,*) 'Parent is COSMO '
      llm2lm = .TRUE.
      lcm2lm = .FALSE.
      lushift_in(1) =.TRUE.
      lushift_in(2) =.FALSE.
      lvshift_in(1) =.FALSE.
      lvshift_in(2) =.TRUE.
      lcm_hgt_coor  = .FALSE.
      lcm_pres_coor = .FALSE.
    ENDIF

       linitial     =.TRUE.
       lcomp_bound  =.FALSE.
       lbd_frame_cur=.FALSE.

    status = 0

    ALLOCATE (slatitudes_in(je_in_tot))
    ALLOCATE (slongitudes_in(ie_in_tot))
    IF (llm2lm .OR. lum2lm ) THEN
       DO i = 1+west_add_in, ie_in_tot-east_add_in
          slongitudes_in(i) = longitudes_in(i) + 0.5_dp * dlon_in
       ENDDO
       DO i = 1+south_add_in, je_in_tot-north_add_in
          slatitudes_in(i) = latitudes_in(i) + 0.5_dp * dlat_in
       ENDDO
    ELSE
       slatitudes_in  = latitudes_in
       slongitudes_in = longitudes_in
    ENDIF

    CALL   setup_int2lm (status, yerror)
    IF(status /= 0)  THEN
      write(iouerr,*)  'Error in setup_int2lm ',status,trim(yerror)
      CALL error_bi ('Error in setup_int2lm ',substr)
    ENDIF

    IF (.NOT. lstart) THEN
       IF (lbd_frame) THEN
          lbd_frame_cur=.TRUE.
       ELSE
          lbd_frame_cur=.FALSE.
       ENDIF
    ENDIF

    IF (l_test) CALL MMD_testC_Setup (INT(ie_in),INT(je_in))

    ! Distribute ak_in and bk_in, compute akh_in ans bkh_in
    IF(MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN
       ! avoid error in pp_lm calculation:
       IF (pcontrol_fi < 0.) &
            CALL error_bi('setting for pcontrol_fi is required for E2C' &
            , substr)
    END IF

    ! acquire ps_gl  / ps_lm as 4D pointer (required for vertical interpolation
    ! of additional fields via NCREGRID
    CALL get_channel_object(status, 'MMDC4', 'PS_GL', p4=psgl4d)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'MMDC4', 'PS',    p4=pslm4d)
    CALL channel_halt(substr, status)

   RETURN

 END SUBROUTINE MMD2WAY_CHILD_SETUP_INT2COSMO

   ! ------------------------------------------------------------------

 SUBROUTINE compute_akh_bkh

   USE data_grid_in,  ONLY: ke_in_tot,  ak_in,  bk_in, akh_in, bkh_in &
                          , dak_in, dbk_in

   IMPLICIT NONE

   INTEGER :: k

   DO k = 1, ke_in_tot
      akh_in(k) = (ak_in(k) + ak_in(k+1)) * 0.5_dp
      bkh_in(k) = (bk_in(k) + bk_in(k+1)) * 0.5_dp
      dak_in(k) =  ak_in(k+1) - ak_in(k)
      dbk_in(k) =  bk_in(k+1) - bk_in(k)
   ENDDO

   RETURN

 END SUBROUTINE compute_akh_bkh

 ! ----------------------------------------------------------------------
 ! ----------------------------------------------------------------------
  SUBROUTINE Setup_data_exchange_with_Parent

    ! INT2COSMO
    USE data_grid_in,             ONLY: latitudes_in,longitudes_in        &
                                      , ie_in_max,je_in_max, ie_in, je_in
#ifdef MESSYTWOWAY
    USE data_fields_in,           ONLY: lon_in, lat_in
#endif
    USE data_int2lm_parallel,     ONLY: isubpos_coarse
    ! MMD
    USE mmd_child,               ONLY: MMD_C_Get_Indexlist   &
                                      , MMD_Send_to_Parent
    ! MESSy/BMIL
    USE messy_main_mpi_bi,        ONLY: num_compute, icomm_compute &
                                      , imp_reals, gather_values, my_cart_id

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! Local Data
    CHARACTER(LEN=*), PARAMETER :: substr ='Setup_data_exchange_with_Parent'
    INTEGER                     :: i,j
    INTEGER                     :: status
    INTEGER(I4)                 :: stat_i4
    INTEGER,DIMENSION(3)        :: size_of_array
    CHARACTER(LEN=80)           :: yerrmsg

    REAL(kind=dp),DIMENSION(ie_in_max, je_in_max,num_compute) :: &
         all_lon, all_lat

    ALLOCATE(my_lon(ie_in_max, je_in_max))
    ALLOCATE(my_lat(ie_in_max, je_in_max))
    my_lon  = -9999.0_dp
    my_lat  = -9999.0_dp
    all_lon = -9999.0_dp
    all_lat = -9999.0_dp

    ! Compute coordinates of mid point of all local cells of the coarse grid
    DO j=1,je_in
      DO i=1,ie_in
        my_lon(i,j) = longitudes_in(isubpos_coarse(my_cart_id,1)+i-1)
        my_lat(i,j) = latitudes_in(isubpos_coarse(my_cart_id,2)+j-1)
        IF (my_lat(i,j) < -400._dp .or. my_lon(i,j) < -400._dp ) &
             CALL error_bi('longitudes_in undefined','setup data exchange')
      ENDDO
    ENDDO
#ifdef MESSYTWOWAY
    lon_in(1:ie_in,1:je_in) = my_lon(1:ie_in,1:je_in)
    lat_in(1:ie_in,1:je_in) = my_lat(1:ie_in,1:je_in)
#endif
    CALL gather_values (my_lon, all_lon, ie_in_max, je_in_max, num_compute, &
                        imp_reals, 0_I4, icomm_compute, yerrmsg, stat_i4)
    IF (stat_i4 /= 0) CALL error_bi(yerrmsg,substr)
    CALL gather_values (my_lat, all_lat, ie_in_max, je_in_max, num_compute, &
                        imp_reals, 0_I4, icomm_compute, yerrmsg, stat_i4)
    IF (stat_i4 /= 0) CALL error_bi(yerrmsg,substr)

    ! send coordinates to parent
    IF (p_parallel_io) THEN
      size_of_array(1) = size(all_lon,1)
      size_of_array(2) = size(all_lon,2)
      size_of_array(3) = size(all_lon,3)
      CALL MMD_Send_to_Parent (size_of_array, 3, 0, 40, status)
      CALL MMD_Send_to_Parent (all_lon, size(all_lon), 0, 41, status)
      IF (status /= 0) CALL error_bi('MMD_Send_to_Parent ERROR 01', substr)
      CALL MMD_Send_to_Parent (all_lat, size(all_lat), 0, 42, status)
      IF (status /= 0) CALL error_bi('MMD_Send_to_Parent ERROR 02', substr)
    ENDIF

    ! Get Data from Parent model
    CALL MMD_C_Get_Indexlist

    RETURN

  END SUBROUTINE Setup_data_exchange_with_parent
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE mmd2way_child_set_CplData

    ! MESSy/BMIL
    USE messy_main_grid_def_mem_bi,  ONLY: ie, je
    USE messy_main_data_bi,          ONLY: l2tls
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: DC_I2C_IN, DC_I2C                &
                                         , gp_meml, gp_nseg, gp_start       &
                                         , DC_GP, DIMID_LON, DIMID_LAT      &
                                         , REPR_I2C_2D, REPR_I2C_3D_MID     &
                                         , REPR_I2C_2D_IN, REPR_3D_MID_IN   &
                                         , DIMID_I2C_IE_IN, DIMID_I2C_JE_IN &
                                         , DIMID_I2C_KE_IN, DIMID_I2C_KEDIM &
                                         , DIMID_I2C_IE, DIMID_I2C_JE
    USE messy_main_tracer_mem_bi,    ONLY: xt, xtm1, xtf, ti_gp, GPTRSTR
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    ! MESSY/SMCL
    USE messy_main_channel,       ONLY: get_channel_object               &
                                      , new_channel_object               &
                                      , get_channel_object_info          &
                                      , new_attribute, get_attribute

    USE messy_main_channel_dimensions,   ONLY: new_dimension             &
                                             , get_dimension_info        &
                                             , ADD_DIMENSION_VARIABLE    &
                                             , ADD_DIMENSION_VARIABLE_ATT
    USE messy_main_channel_repr,         ONLY: new_representation        &
                                             , set_representation_decomp &
                                             , AUTO, IRANK, PIOTYPE_COL  &
                                             , get_representation_info   &
                                             , repr_def_axes
    USE messy_main_tracer,        ONLY: get_tracer, full2base_sub, I_MMD_INIT &
                                      , ON, set_tracer, I_LBC, T_LBC_CST
    USE messy_main_timer,         ONLY: lstart
    USE messy_main_tools,         ONLY: int2str, str, strcrack,  str2num
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, STRLEN_ULONG
    ! INT2COSMO
    USE data_int2lm_control,      ONLY: lcm2lm
    USE data_int2lm_io,           ONLY: nvar_in, var_in, nvar_lm, var_lm
    USE data_grid_lm,             ONLY: je2lm, ie2lm, kedim
    USE data_grid_in,             ONLY: je_in, ie_in, ke_in

    ! MMD
    USE mmd_test,                 ONLY: MMD_testC_GetTestPtr
    USE mmd_child,                ONLY: MMD_C_Get_Repr      &
                                      , MMD_C_GetParentType &
                                      , MMD_ParentIsEcham
    USE mmd_utilities,            ONLY: STRLEN_ATT
    USE mmd_mpi_wrapper,          ONLY: MMD_Inter_Bcast

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, SIZE, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER        :: substr='mmd2way_child_set_CplData'

    INTEGER                            :: status
    INTEGER                            :: itab, ii, il, jj, l, ic
    INTEGER                            :: idt    ! tracer idt
    CHARACTER(LEN=STRLEN_MEDIUM)       :: trbasename, chbasename
    CHARACTER(LEN=STRLEN_MEDIUM)       :: trsubname, chsubname
    CHARACTER(LEN=STRLEN_MEDIUM)       :: objname
    CHARACTER(LEN=STRLEN_MEDIUM)       :: COSMOOBJ
    INTEGER                            :: reprid     = 0
    INTEGER                            :: repr_input = 0
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: bdptr => NULL()  ! boundary pointer
    CHARACTER(LEN=4)                   :: par_axis
    INTEGER, DIMENSION(4)              :: par_gdimlen
    CHARACTER(LEN=STRLEN_MEDIUM)       :: par_repr
    CHARACTER(LEN=STRLEN_ATT)          :: par_att
    CHARACTER(LEN=2)                   :: dimstr
    INTEGER                            :: zdim_length
    INTEGER                            :: ndim_length
    INTEGER, DIMENSION(4)              :: dim_ids
    CHARACTER(LEN=4)                   :: dim_axis
    INTEGER, DIMENSION(4)              :: dim_len
    INTEGER                            :: length, tmp
    CHARACTER(LEN=1)                   :: stmp
    INTEGER                            :: rank
    INTEGER                            :: yind, zind, nind
    INTEGER                            :: REPR_I2C = 0
    INTEGER                            :: REPR_IN  = 0
    INTEGER                            :: timelev
    INTEGER                            :: ntiles
    CHARACTER(LEN=STRLEN_ULONG)        :: unit
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: chlevs => NULL()
    REAL(dp),                     DIMENSION(:), POINTER ::   levs => NULL()

    CplData_loop: DO ii = 1 , NCOPY
       ! SET DEFAULTS
       CplData(ii)%lvartab     =.FALSE.
       CplData(ii)%vartab_name = ''
       CplData(ii)%vartab_idx  = 0
       CplData(ii)%rank        = 0
       NULLIFY(CplData(ii)%ptr_in)
       NULLIFY(CplData(ii)%ptr_i2c)

       ! GET COSMO CHANNEL OBJECT AND COSMO BOUNDARY CHANNEL OBJECT
       !write(0,*) substr, ' ',TRIM(CplData(ii)%CHILD%CHA) &
       !                  , ' ',TRIM(CplData(ii)%CHILD%OBJ)
       par_att = ''
       if_ccha: IF (TRIM(CplData(ii)%CHILD%CHA)/='#XXX') THEN
          if_ccha2: IF ((TRIM(CplData(ii)%CHILD%CHA) == 'mmd2way_child')  &
               .AND. CplData(ii)%L_INPUT) THEN
             ! FIND REPRESENTATION IF UNKNOWN
             CALL MMD_Inter_Bcast(unit, .FALSE.)

             IF (TRIM(ADJUSTL(CplData(ii)%C_REPR)) == '#UNKNOWN') THEN
                CALL MMD_C_Get_Repr(par_axis, par_gdimlen &
                                  , par_repr, par_att)
                !write(iouerr,*) 'LOOKING FOR REPRESENTATION ' &
                !     , TRIM(CplData(ii)%CHILD%OBJ),'' , par_axis &
                !     , par_gdimlen, TRIM(par_repr)
                ! AS MMD only transport field with XY DIMENSION WE
                ! ASSUME THEY EXIST HERE AS WELL
                IF ((TRIM(par_repr) /= 'GP_3D_MID') .AND. &
                     (TRIM(par_repr) /= 'GP_2D_HORIZONTAL') .AND. &
                     (TRIM(par_repr) /= 'GP_3D_1LEV')) THEN
                   dim_ids(:) = 0
                   dim_axis = '----'
                   dim_len(:) = 0
                   zdim_length = 0
                   ndim_length = 0
                   rank = 0
                   yind = -1
                   zind = -1
                   nind = -1
                   DO il=1,4
                      IF (par_axis(il:il) == 'X') THEN
                         rank = rank+1
                         dim_ids(il)=DIMID_LON
                         dim_axis(il:il)='X'
                         dim_len(il)=ie
                      ELSE IF (par_axis(il:il) == 'Y') THEN
                         rank = rank+1
                         yind = il
                         dim_ids(il)=DIMID_LAT
                         dim_axis(il:il)='Y'
                         dim_len(il)=je
                      ELSE IF (par_axis(il:il) == 'Z') THEN
                         rank = rank+1
                         zind = il
                         zdim_length = par_gdimlen(il)
                         CALL int2str(dimstr,zdim_length)
                         CALL get_dimension_info(status      &
                              , 'DIM_'//dimstr//'LEV'   &
                              , id=dim_ids(il), len=length)
                         IF (status /=0) THEN
                            CALL new_dimension(status, dim_ids(il) &
                                , 'DIM_'//dimstr//'LEV', zdim_length)
                            CALL channel_halt(substr,status)
                            CALL STRCRACK(TRIM(par_att), ',', chlevs, l)
                            ALLOCATE(levs(SIZE(chlevs)))
                            DO ic = 2, SIZE(chlevs)
                               CALL str2num(chlevs(ic), levs(ic))
                            END DO
                            CALL add_dimension_variable(status, dim_ids(il) &
                                 , 'DIM_'//dimstr//'LEV', levs)
                            CALL CHANNEL_HALT(substr, status)
                            CALL add_dimension_variable_att(status,dim_ids(il)&
                                 , 'DIM_'//dimstr//'LEV', 'type' &
                                 , c='heights')
                            CALL CHANNEL_HALT(substr, status)
                            DEALLOCATE(chlevs); NULLIFY(chlevs)
                            DEALLOCATE(levs);   NULLIFY(levs)
                         ELSE
                            IF (length /= zdim_length) THEN
                               CALL error_bi(' OVERLAP OF VERITCAL DIMENSIONS' &
                                    , substr)
                            ENDIF
                         ENDIF
                         dim_axis(il:il)='Z'
                         dim_len(il)=AUTO
                      ELSE IF (par_axis(il:il) == 'N') THEN
                         rank = rank+1
                         nind = il
                         ndim_length = par_gdimlen(il)
                         CALL int2str(dimstr,ndim_length)
                         CALL get_dimension_info(status &
                              , 'DIM_'//dimstr//'N'     &
                              , id=dim_ids(il), len=length)
                         IF (status /=0) THEN
                            CALL new_dimension(status, dim_ids(il) &
                                 , 'DIM_'//dimstr//'N', ndim_length)
                            CALL channel_halt(substr,status)
                         ELSE
                            IF (length /= ndim_length) THEN
                               CALL error_bi( &
                                    ' OVERLAP OF FREE DIMENSIONS', substr)
                            ENDIF
                         ENDIF
                         dim_axis(il:il)='N'
                         dim_len(il)=AUTO
                      ELSE IF (par_axis(il:il) == '-') THEN
                         dim_ids(il) = 1
                      ELSE
                         CALL error_bi('UNDEFINED DIMENSION', substr)
                      END IF
                   ENDDO

                   IF (rank == 2 .AND. zdim_length== 0 .AND. &
                        ndim_length == 0) THEN
                         CplData(ii)%C_REPR='GP_2D_HORIZONTAL'
                   ELSE IF (rank > 2 .AND. rank <=4 .AND. &
                        dim_axis(1:1) /= '-' .AND. dim_axis(2:2) /= '-') THEN
                      ! switch Y and Z/N axis, if PARENT is ECHAM
                      IF(MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN
                         IF (_RI_TRRANK == 3) THEN
                            IF (zind /= -1) THEN  ! Z dimension exists
                               tmp            = dim_ids(zind)
                               dim_ids(zind)  = dim_ids(yind)
                               dim_ids(yind)  = tmp
                               tmp            = dim_len(zind)
                               dim_len(zind)  = dim_len(yind)
                               dim_len(yind)  = tmp
                               stmp                = dim_axis(zind:zind)
                               dim_axis(zind:zind) = dim_axis(yind:yind)
                               dim_axis(yind:yind) = stmp
                            ELSE IF(nind /= -1) THEN !???
                               tmp            = dim_ids(nind)
                               dim_ids(nind)  = dim_ids(yind)
                               dim_ids(yind)  = tmp
                               tmp            = dim_len(nind)
                               dim_len(nind)  = dim_len(yind)
                               dim_len(yind)  = tmp
                               stmp                = dim_axis(nind:nind)
                               dim_axis(nind:nind) = dim_axis(yind:yind)
                               dim_axis(yind:yind) = stmp
                            ENDIF
                         ELSE IF (_RI_TRRANK ==  4) THEN
                            IF (zind /= -1) THEN  ! Z dimension exists
                               tmp            = dim_ids(zind)
                               dim_ids(zind)  = dim_ids(yind)
                               dim_ids(yind)  = tmp
                               tmp            = dim_len(zind)
                               dim_len(zind)  = dim_len(yind)
                               dim_len(yind)  = tmp
                               stmp                = dim_axis(zind:zind)
                               dim_axis(zind:zind) = dim_axis(yind:yind)
                               dim_axis(yind:yind) = stmp
                            END IF
                            IF(nind /= -1) THEN !???
                               tmp            = dim_ids(nind)
                               dim_ids(nind)  = dim_ids(yind)
                               dim_ids(yind)  = tmp
                               tmp            = dim_len(nind)
                               dim_len(nind)  = dim_len(yind)
                               dim_len(yind)  = tmp
                               stmp                = dim_axis(nind:nind)
                               dim_axis(nind:nind) = dim_axis(yind:yind)
                               dim_axis(yind:yind) = stmp
                            ENDIF

                         ELSE
                            CALL error_bi(substr, ' UNKNOWN TRACER RANK')
                         ENDIF
                      ENDIF
                      CALL make_cosmo_representation(rank,zdim_length &
                           , ndim_length, dim_ids, dim_axis, dim_len  &
                           , reprid)
                      ! get name of representation
                      CALL get_representation_info(status, ' '   &
                           , name=CplData(ii)%C_REPR, id=reprid  &
                           , axis=CplData(ii)%AXIS               &
                           , ldimlen=CplData(ii)%ldimlen         &
                           , rank= CplData(ii)%rank              )
                      CALL channel_halt(substr, status)
                   ELSE
                      CALL error_bi('NO MEANINGFUL REPRESENTATION FOUND' &
                           , substr)
                   ENDIF
                ELSE
                   CplData(ii)%C_REPR=TRIM(par_repr)
                   IF (TRIM(par_repr) == 'GP_3D_1LEV') &
                        CplData(ii)%C_REPR='GP_2D_HORIZONTAL'
                ENDIF
             ENDIF ! REPR UNKNOWN

             CALL get_representation_info(status, CplData(ii)%C_REPR   &
                  , id=repr_input, rank=CplData(ii)%rank               )
             CALL channel_halt(substr, status)

             ! NOTE: L_INPUT FIELDS CANNOT BE PROGNOSTIC
             ALLOCATE(CplData(ii)%cosmo(1))
             CALL new_channel_object(status, submodstr                      &
                  , CplData(ii)%CHILD%OBJ, p4=CplData(ii)%cosmo(1)%ptr      &
                  , reprid=repr_input)
             CALL channel_halt(substr, status)

             IF (CplData(ii)%l_SendUnit) THEN
                CALL new_attribute(status, submodstr, CplData(ii)%CHILD%OBJ &
                      , 'units', c=TRIM(unit) )
                CALL channel_halt(substr, status)
              ENDIF
          ELSE ! if_ccha2
             !
             ! ALLOCATE PTR ARRAY FOR COSMO FIELD
             !
             ! INQUIRE IF FIELD IS PROGNOSTIC, i.e. search for attribute
             ! "number_of_timelevels". If it exists, the object is time
             ! dependent. IF status /= 0 the object is diagnostic
             !
             COSMOOBJ = TRIM(CplData(ii)%CHILD%OBJ)
             CALL get_channel_object(status      &
                  , TRIM(CplData(ii)%CHILD%CHA)  &
                  , TRIM(COSMOOBJ))
             IF (status == 3103) THEN
                COSMOOBJ = TRIM(CplData(ii)%CHILD%OBJ)//'_tile0'
                CALL get_channel_object(status      &
                     , TRIM(CplData(ii)%CHILD%CHA)  &
                     , TRIM(COSMOOBJ))
             END IF
             CALL channel_halt(substr, status) 

             CALL get_attribute(status                &
                  , TRIM(CplData(ii)%CHILD%CHA)       &
!                  , TRIM(CplData(ii)%CHILD%OBJ)       &
                  , TRIM(COSMOOBJ)                    &
                  , 'number_of_timelevels', i=timelev )

             nottimedep: IF (status /= 0) THEN !
                ALLOCATE(CplData(ii)%cosmo(1))
                ! COSMO CHANNEL OBJECT
                CALL get_channel_object(status      &
                     , TRIM(CplData(ii)%CHILD%CHA)  &
                     , TRIM(COSMOOBJ)  &
                     , p4=CplData(ii)%cosmo(1)%ptr  )
                IF (status /=0) THEN
                   write(iouerr,*) substr,' COSMO_set_CplData 1', status &
                        , TRIM(CplData(ii)%CHILD%CHA),' '           &
                        , TRIM(COSMOOBJ)
                   IF (status == 3103) THEN

                   END IF
                   CALL channel_halt(substr, status)
                ENDIF
                ! get rank of object
                CALL get_rank(status, TRIM(CplData(ii)%CHILD%CHA) &
!                     , TRIM(CplData(ii)%CHILD%OBJ), CplData(ii)%rank)
                     , TRIM(COSMOOBJ), CplData(ii)%rank)
             ELSE ! nottimdep
                ! FIELD IS PROGNOSTIC / HAS DIFFERENT TIME LEVELS
                if_trac: IF (TRIM(CplData(ii)%CHILD%CHA) == 'tracer_gp') THEN
                   ! SEPCIAL TREATMENTS OF TRACERS REQUIRED
                   ALLOCATE(CplData(ii)%cosmo(2))
                   ! get tracer idt
                   CALL full2base_sub(status, TRIM(CplData(ii)%CHILD%CHA)     &
                        , chbasename, chsubname)
                   CALL full2base_sub(status, TRIM(CplData(ii)%CHILD%OBJ)     &
                        , trbasename, trsubname)
                   CALL get_tracer( status, TRIM(chsubname), TRIM(trbasename) &
                        , TRIM(trsubname), idx=idt)
                   CALL tracer_halt(substr,status)

                   ! INDEX 1 always nnew
                   CplData(ii)%cosmo(1)%ptr => xt(_RI_XYZN_(:,:,:,idt:idt))
                   ! INDEX 2 always nnow
                   IF (l2tls) THEN
                      CplData(ii)%cosmo(2)%ptr => xtm1(_RI_XYZN_(:,:,:,idt:idt))
                   ELSE
                      CplData(ii)%cosmo(2)%ptr => xtf(_RI_XYZN_(:,:,:,idt:idt))
                   ENDIF
                   CplData(ii)%rank = 3

                   ! SET INITIALISATION FLAG
                   IF (lstart .AND. CplData(ii)%L_INITIAL)  &
                        ti_gp(idt)%tp%meta%cask_i(I_MMD_INIT) = ON

                ELSE ! if tracer
                   ALLOCATE(CplData(ii)%cosmo(timelev))
                   DO jj=1,timelev
                      objname=TRIM(CplData(ii)%CHILD%OBJ)//'_tl'//str(jj)

                      ! COSMO CHANNEL OBJECT
                      CALL get_channel_object(status                      &
                           , TRIM(CplData(ii)%CHILD%CHA), objname         &
                           , p4=CplData(ii)%cosmo(jj)%ptr )
                      IF (status /=0 .AND. status /= 2020) THEN
                         write(iouerr,*) substr,' COSMO_set_CplData 2', status &
                              , TRIM(CplData(ii)%CHILD%CHA),' ', TRIM(OBJNAME)
                         CALL channel_halt(substr, status)
                      ENDIF
                      IF (jj==1) &
                           CALL get_rank(status,TRIM(CplData(ii)%CHILD%CHA) &
                           , TRIM(objname), CplData(ii)%rank)
                   ENDDO
                END IF if_trac
             ENDIF nottimedep

             if_lbound: IF (CplData(ii)%L_BOUND) THEN
                ALLOCATE(CplData(ii)%cosmo_bd(2))
                IF (CplData(ii)%rank == 2 .AND.  &
                     TRIM(CplData(ii)%C_repr) /='GP_3D_1LEV') THEN
                   CALL get_channel_object(status             &
                        , TRIM(CplData(ii)%CHILD%CHA)         &
                        , TRIM(CplData(ii)%CHILD%OBJ)//'_BD'  &
                        , p4=bdptr                            )
                   IF (status /=0) write(iouerr,*) substr,' COSMO_BD '  &
                        , TRIM(CplData(ii)%CHILD%CHA),' '           &
                        , TRIM(CplData(ii)%CHILD%OBJ)//'_BD'
                   CALL channel_halt(substr, status)
                   CplData(ii)%cosmo_bd(1)%ptr =>bdptr (:,:,1:1,1:1)
                   CplData(ii)%cosmo_bd(2)%ptr =>bdptr (:,:,2:2,1:1)
                   NULLIFY(bdptr)
                ELSE IF ( (CplData(ii)%rank == 3 .OR.                 &
                     TRIM(CplData(ii)%C_repr) /='GP_3D_1LEV' ) .AND.  &
                     (TRIM(CplData(ii)%CHILD%CHA) /= 'tracer_gp')) THEN
                   CALL get_channel_object(status             &
                        , TRIM(CplData(ii)%CHILD%CHA)         &
                        , TRIM(CplData(ii)%CHILD%OBJ)//'_BD'  &
                        , p4=bdptr                            )
                   IF (status /=0) write(iouerr,*) substr,' COSMO_BD '  &
                        , TRIM(CplData(ii)%CHILD%CHA),' '           &
                        , TRIM(CplData(ii)%CHILD%OBJ)//'_BD'
                   CALL channel_halt(substr, status)
#ifdef COSMOv509
                   ntiles = -99
                   CALL get_attribute(status                  &
                        , TRIM(CplData(ii)%CHILD%CHA)         &
                        , TRIM(CplData(ii)%CHILD%OBJ)         &
                        , 'number_of_tiles', i=ntiles)
                   IF (status == 0 .AND. ntiles >= 0) THEN
                      CplData(ii)%ltiles = .TRUE.
                      CplData(ii)%cosmo_bd(1)%ptr =>bdptr (:,:,1:1,1:1)
                      CplData(ii)%cosmo_bd(2)%ptr =>bdptr (:,:,2:2,1:1)
                   ELSE
#endif
                   CplData(ii)%cosmo_bd(1)%ptr =>bdptr (:,:,:,1:1)
                   CplData(ii)%cosmo_bd(2)%ptr =>bdptr (:,:,:,2:2)
#ifdef COSMOv509
                   ENDIF
#endif
                   NULLIFY(bdptr)
                ELSE IF (TRIM(CplData(ii)%CHILD%CHA) == 'tracer_gp') THEN
                   CALL get_channel_object(status              &
                        , TRIM(CplData(ii)%CHILD%CHA)//'_x001' &
                        , TRIM(CplData(ii)%CHILD%OBJ)          &
                        , p4=CplData(ii)%cosmo_bd(1)%ptr        )
                   IF (status /=0) write(iouerr,*) substr,' COSMO_BD '  &
                        ,TRIM(CplData(ii)%CHILD%CHA)//'_x001',' '   &
                        ,TRIM(CplData(ii)%CHILD%OBJ)
                   CALL channel_halt(substr, status)
                   CALL get_channel_object(status              &
                        , TRIM(CplData(ii)%CHILD%CHA)//'_x002' &
                        , TRIM(CplData(ii)%CHILD%OBJ)          &
                        , p4=CplData(ii)%cosmo_bd(2)%ptr        )
                   IF (status /=0) write(iouerr,*) substr,' COSMO_BD '  &
                        ,TRIM(CplData(ii)%CHILD%CHA)//'_x002',' '   &
                        ,TRIM(CplData(ii)%CHILD%OBJ)
                   CALL channel_halt(substr, status)

                ELSE
                   CALL error_bi(&
              'channel name of 4D prognostic element needs to be implemented'&
                        , substr)
                ENDIF
             ELSE
                ! reset tracer property for boundary condition
                IF (TRIM(CplData(ii)%CHILD%CHA) == 'tracer_gp') THEN
                   CALL set_tracer(status, GPTRSTR, idt, I_LBC, T_LBC_CST)
                   CALL tracer_halt(substr, status)
                END IF

             ENDIF if_lbound

          ENDIF if_ccha2
       ENDIF if_ccha

       ! SET VARTAB NAME and INIDCATOR OF IN2COSMO variable "lvartab"
       ! SPECIAL CASE ROUGHNESS LENGTH in INT2COSMO is g*z0 in COSMO
        vartab_name: DO itab = 1, nvar_lm
          IF (TRIM(CplData(ii)%CHILD%OBJ) == TRIM(var_lm(itab)%name)) THEN
             CplData(ii)%vartab_name =  TRIM(var_lm(itab)%name)
             CplData(ii)%lvartab     = .TRUE.
             CplData(ii)%vartab_idx  = itab
#ifdef COSMOv509
             IF (TRIM(CplData(ii)%CHILD%OBJ) == 'W_SO' .OR. &
                  TRIM(CplData(ii)%CHILD%OBJ) == 'T_SO') THEN
                CplData(ii)%ltiles    = .TRUE.
             END IF
#endif
             EXIT
          ELSE IF (TRIM(CplData(ii)%CHILD%OBJ) == 'gZ0' .AND. &
               TRIM(var_lm(itab)%name) == 'Z0') THEN
             CplData(ii)%vartab_name =  'Z0'
             CplData(ii)%lvartab     = .TRUE.
             CplData(ii)%vartab_idx  = itab
             EXIT
          ENDIF
       ENDDO vartab_name

       if_lvartab: IF (CplData(ii)%lvartab) THEN
          ! "IN" POINTER IS NOT REQUIRED FOR COSMO FIELDS ADDITIONALLY
          ! CALCULATED WITHIN INT2COSMO
          IF (ii <= NEXCH ) THEN
             ! FIND "IN" POINTER
             CALL get_channel_object(status, 'MMDC4_IN'  &
                  , TRIM(CplData(ii)%vartab_name)//'_IN' &
                  , p4 = CplData(ii)%ptr_in)
             IF (status /= 0) THEN
                IF (TRIM(CplData(ii)%vartab_name) == 'W_SO') THEN
                   CALL get_channel_object(status, 'MMDC4_IN'       &
                        , TRIM(CplData(ii)%vartab_name)//'_REL_IN'  &
                        , p4 = CplData(ii)%ptr_in)
                   IF (status /= 0) write(iouerr,*) substr,' MMDC4_IN ' &
                        , TRIM(CplData(ii)%vartab_name)//'_REL_IN'
                   CALL channel_halt(substr, status)
                   CALL get_channel_object_info(status, 'MMDC4_IN'  &
                        , TRIM(CplData(ii)%vartab_name)//'_REL_IN'  &
                        , reprid=reprid)
                   CALL channel_halt(substr, status)
                ELSE
                   write(iouerr,*) substr,' MMDC4_IN ' &
                        , TRIM(CplData(ii)%vartab_name)//'_IN'
                   CALL channel_halt(substr, status)
                ENDIF
             ELSE
                CALL get_channel_object_info(status, 'MMDC4_IN' &
                     , TRIM(CplData(ii)%vartab_name)//'_IN', reprid=reprid)
                CALL channel_halt(substr, status)
             ENDIF

             ! DEFINE AXIS STRING
             ! 1. get representation id
             ! 2. get axis string
             CALL get_representation_info(status, ' ', reprid &
                  , axis=CplData(ii)%AXIS, ldimlen=CplData(ii)%ldimlen)
             CALL channel_halt(substr, status)

             IF (SIZE(CplData(ii)%ptr_in,3) == 1 ) THEN
                CplData(ii)%rank = 2
             ELSE
                CplData(ii)%rank = 3
             ENDIF
             ! SPECIAL CASE W_SO set lreadin correctly
             vtab: DO itab = 1, nvar_in
                IF (TRIM(var_in(itab)%name) == &
                     TRIM(CplData(ii)%vartab_name)) THEN
                   IF (lcm2lm .AND. TRIM(CplData(ii)%vartab_name) == 'W_SO') THEN
                      var_in(35)%lreadin = .TRUE. ! qqq 35 hard coded
                      EXIT
                   ELSE
                      var_in(itab)%lreadin = .TRUE.
                      EXIT
                   ENDIF
                END IF
             ENDDO vtab
             IF (itab == nvar_in + 1 ) THEN
                IF (TRIM(CplData(ii)%vartab_name) == 'W_SO') THEN
                     var_in(35)%lreadin = .TRUE. ! qqq 35 hard coded
                ENDIF
             ENDIF

             ! OVERWRITE INTERPOLATION FLAGS SET IN GRIBTRABS
             IF (TRIM(CplData(ii)%C_INTERPOL(1:3)) /= '') THEN
                DO itab = 1, nvar_in
                   IF (TRIM(var_in(itab)%name) == &
                        TRIM(CplData(ii)%vartab_name)) THEN
                      IF (p_parallel_io) write (*,*) 'WARNING: '     &
                           , 'CHANGE INTERPOLATION OF INT2LM FIELD ' &
                           , TRIM(var_in(itab)%name), ' FROM ',var_in(itab)%ipc
                      IF (TRIM(CplData(ii)%C_INTERPOL(1:1)) /= '') THEN
                         var_in(itab)%ipc(1:1) = CplData(ii)%C_INTERPOL(1:1)
                      END IF
                      IF (TRIM(CplData(ii)%C_INTERPOL(2:2)) /= '')  THEN
                         var_in(itab)%ipc(2:2) = CplData(ii)%C_INTERPOL(2:2)
                      END IF
                      IF (TRIM(CplData(ii)%C_INTERPOL(3:3)) /= '')  THEN
                         var_in(itab)%ipc(3:3) = CplData(ii)%C_INTERPOL(3:3)
                      END IF
                      IF (p_parallel_io) write (*,*) &
                           ' ... TO ', var_in(itab)%ipc
                      EXIT
                   ENDIF
                ENDDO
             END IF
          ENDIF

          ! FIND INTERMEDIATE POINTER OF INT2COSMO
          CALL get_channel_object(status, 'MMDC4'                 &
               , CplData(ii)%vartab_name, p4 = CplData(ii)%ptr_i2c)
          IF (status /= 0) THEN
             write(iouerr,*) substr,' INT2COSMO POINTER ' &
                  ,TRIM(CplData(ii)%CHILD%CHA),' '    &
                  ,TRIM(CplData(ii)%CHILD%OBJ)
          ENDIF
          CALL channel_halt(substr, status)
       ELSE ! if_lvartab
          ! A) MAKE/FIND REPRESENTATION
          ifrank: IF  (CplData(ii)%rank == 2) THEN
             REPR_IN  = REPR_I2C_2D_IN
             REPR_I2C = REPR_I2C_2D
          ELSE IF (CplData(ii)%rank == 3) THEN
             ! SPECIAL CASE TEST ARRAY, SET PTR_IN BUT NOT PTR_I2C
             if_test: IF (TRIM(CplData(ii)%CHILD%CHA) == submodstr .AND. &
                  TRIM(CplData(ii)%CHILD%OBJ) == 'Test_Ar' .AND. l_test) THEN
                ! CHANNEL OBJECT IS MADE BY THIS MMD2WAY_CHILD
                CALL MMD_testC_GetTestPtr(CplData(ii)%ptr_in &
                     , CplData(ii)%AXIS, CplData(ii)%ldimlen)
                NULLIFY(CplData(ii)%ptr_i2c)
                CYCLE
             ENDIF if_test
             IF (TRIM(CplData(ii)%C_REPR) == 'GP_3D_MID') THEN
                ! IN THIS CASE THE REPRESENTATIONS ARE ALREADY KNOWN
                REPR_IN  = REPR_3D_MID_IN
                REPR_I2C = REPR_I2C_3D_MID
             ELSE
                ! FIND/ MAKE REPRESENTATION
                IF (CplData(ii)%axis(3:3) == 'Z') THEN
                   CALL make_i2c_representations(REPR_I2C, REPR_IN &
                        , vsize=SIZE(CplData(ii)%cosmo(1)%ptr,_IZ_XYZ__) )
                ELSE IF (CplData(ii)%axis(3:3) == 'N') THEN
                   CALL make_i2c_representations(REPR_I2C, REPR_IN &
                        , nsize=SIZE(CplData(ii)%cosmo(1)%ptr,_IN_XY_N_) )
                ELSE
                   CALL error_bi('CANNOT IDENTIFY REPRESENTATION'//&
                        &' ('//TRIM(CplData(ii)%CHILD%CHA)//' - '//&
                        &TRIM(CplData(ii)%CHILD%OBJ)//')'          &
                        , substr)
                ENDIF
             ENDIF
          ELSE IF (CplData(ii)%rank == 4) THEN
             CALL make_i2c_representations( REPR_I2C, REPR_IN &
                  , nsize=SIZE(CplData(ii)%cosmo(1)%ptr,_IN_XYZN_)    &
                  , vsize=SIZE(CplData(ii)%cosmo(1)%ptr,_IZ_XYZN_)    )
          ENDIF ifrank
          !----------------------------------------------------------
          ! B) MAKE CHANNEL OBJECTS
          ! MAKE MEMORY FOR INPUT FIELD FILLED BY MMD
          CALL new_channel_object(status, 'MMDC4_IN'                       &
               , TRIM(CplData(ii)%CHILD%OBJ)//'_IN', p4=CplData(ii)%ptr_in &
               , reprid=REPR_IN)
          CALL channel_halt(substr, status)
          ! GET INFO ABOUT AXES
          CALL get_representation_info(status, ' ', REPR_IN        &
               , axis=CplData(ii)%AXIS, ldimlen=CplData(ii)%ldimlen)
          CALL channel_halt(substr, status)

          ! MAKE INTERMEDIATE MEMORY NEEDED BY INT2COSMO
          ! FOR ADDITIONAL ARRAYS
          IF (CplData(ii)%rank == 4) THEN
             ! begin object name with marker to switch off output of 4D Objects
             CALL new_channel_object(status, 'MMDC4'       &
                  , '&4D_'//CplData(ii)%CHILD%OBJ          &
                  , p4=CplData(ii)%ptr_i2c, reprid=REPR_I2C)
             CALL channel_halt(substr, status)
          ELSE
             CALL new_channel_object(status, 'MMDC4'      &
                  , CplData(ii)%CHILD%OBJ                 &
                  , p4=CplData(ii)%ptr_i2c, reprid=REPR_I2C)
             CALL channel_halt(substr, status)
          END IF
       END IF if_lvartab
    ENDDO CplData_loop

    RETURN

    ! -----------------------------------------------------------------

    CONTAINS

      SUBROUTINE get_rank(status, cname, oname, rank)

        IMPLICIT NONE

        ! I/O
        INTEGER,                      INTENT(OUT)          :: status
        ! CHANNEL NAME
        CHARACTER(LEN=*),             INTENT(IN)           :: cname
        ! OBJECT NAME
        CHARACTER(LEN=*),             INTENT(IN)           :: oname
        INTEGER,                      INTENT(OUT)          :: rank

        ! LOCAL
        INTEGER :: rep_id ! representation id
        CHARACTER(LEN=*), PARAMETER :: substr = 'get_rank'

        CALL get_channel_object_info(status, cname, oname, reprid=rep_id)
        CALL channel_halt(substr,status)

        CALL get_representation_info(status, ' ', rep_id, rank=rank)
        CALL channel_halt(substr,status)

        status = 0

      END SUBROUTINE get_rank

      !--------------------------------------------------------------

      SUBROUTINE make_i2c_representations( REPR_MMD_I2C,  REPR_MMD_IN &
                                         , nsize, vsize)

        IMPLICIT NONE

        INTRINSIC :: INT, PRESENT

        INTEGER, INTENT(OUT)             :: REPR_MMD_IN
        INTEGER, INTENT(OUT)             :: REPR_MMD_I2C
        INTEGER, INTENT(IN),    OPTIONAL :: nsize
        INTEGER, INTENT(IN),    OPTIONAL :: vsize

        ! PARALLEL DECOMPOSITION
        INTEGER                          :: nseg = 0
        INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
        INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
        INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
        INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

        CHARACTER(LEN=*), PARAMETER :: substr='make_i2c_representations'
        CHARACTER(LEN=2)            :: zchar, nchar, char
        INTEGER                     :: zsize,size, DIMID
        INTEGER                     :: DIMID_Z  = 0
        INTEGER                     :: DIMID_N  = 0
        CHARACTER(LEN=3)            :: tag=''
        INTEGER                     :: rank

        nseg = gp_nseg

        ALLOCATE(start(nseg,IRANK))
        ALLOCATE(cnt(nseg,IRANK))
        ALLOCATE(meml(nseg,IRANK))
        ALLOCATE(memu(nseg,IRANK))

        start(:,:)=1
        cnt(:,:)  =1
        meml(:,:) =1
        memu(:,:) =1

        IF ((PRESENT(vsize)) .AND. (PRESENT(nsize))) THEN
           rank = 4
        ELSE IF (.NOT. PRESENT(vsize) .AND. (.NOT. PRESENT(nsize))) THEN
           CALL error_bi('visze and nsize missing for representation', substr)
        ELSE
           rank = 3
        ENDIF

        IF (PRESENT(nsize)) THEN
           CALL int2str(nchar,nsize)
           CALL get_dimension_info(status, 'DIM_'//nchar//'N' &
                , id=DIMID_N, len=length)
           IF (status /=0) THEN
              CALL new_dimension(status, DIMID_N, 'DIM_'//nchar//'N', nsize)
              CALL channel_halt(substr,status)
           ELSE
              IF (length /= nsize) THEN
                 CALL error_bi(' OVERLAP OF NUMBER DIMENSIONS', substr)
              ENDIF
           ENDIF
        ENDIF

        IF (PRESENT(vsize)) THEN
           IF (vsize == 0) THEN
              DIMID_Z=DIMID_I2C_KE_IN
              zsize=ke_in
              zchar=''
           ELSE
              zsize=vsize
              CALL int2str(zchar,zsize)
              CALL get_dimension_info(status, 'DIM_'//zchar//'LEV'  &
                   , id=DIMID_Z, len=length)
              IF (status /=0) THEN
                 CALL new_dimension(status, DIMID_Z, 'DIM_'//zchar//'LEV', zsize)
                 CALL channel_halt(substr,status)
              ELSE
                 IF (length /= zsize) THEN
                    CALL error_bi(' OVERLAP OF VERTICAL DIMENSIONS', substr)
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
        IF (rank ==4) THEN

           CALL get_representation_info(status                   &
                , 'REPR_MMD_IN_'//TRIM(zchar)//'LEV'//nchar//'N' &
                , id= REPR_MMD_IN                       )

           IF (status /=0 ) THEN
              ! NEW REPRESENTATIONS
              CALL new_representation(status, REPR_MMD_IN                    &
                   , 'REPR_MMD_IN_'//TRIM(zchar)//'LEV'//nchar//'N'          &
                   , rank = 4, link = 'xxxx', dctype = DC_I2C_IN             &
                   , dimension_ids = (/&
                    _RI_XYZN_(DIMID_I2C_IE_IN,DIMID_I2C_JE_IN,DIMID_Z,DIMID_N)&
                     /) &
                   , ldimlen       = (/ &
                      _RI_XYZN_(INT(ie_in), INT(je_in), AUTO, AUTO)/)         &
                   , output_order  = (/&
                       _IX_XYZN_, _IY_XYZN_, _IZ_XYZN_, _IN_XYZN_ /)           &
                   , axis          = repr_def_axes(_RI_XYZN_('X','Y','N','Z')) &
                   )
              CALL channel_halt(substr, status)

              start(:,:) = gp_start(:,:)
              cnt(:,_IX_XYZN_)   = ie_in
              cnt(:,_IY_XYZN_)   = je_in
              cnt(:,_IN_XYZN_)   = nsize
              cnt(:,_IZ_XYZN_)   = zsize
              meml(:,:)  = gp_meml(:,:)
              memu(:,_IX_XYZN_)  = ie_in
              memu(:,_IY_XYZN_)  = je_in
              memu(:,_IN_XYZN_)  = nsize
              memu(:,_IZ_XYZN_)  = zsize

              CALL set_representation_decomp(status, REPR_MMD_IN &
                   , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
              CALL channel_halt(substr, status)
           ENDIF
           !------------------------------------------------------------
           IF (PRESENT(vsize)) THEN
              IF (vsize==0) THEN
                 DIMID_Z=DIMID_I2C_KEDIM
                 zsize=kedim
              ENDIF
           ENDIF

           CALL get_representation_info(status              &
                , 'REPR_MMD_I2C_'//zchar//'LEV'//nchar//'N' &
                , id= REPR_MMD_I2C                          )

           IF (status /=0 ) THEN
              CALL new_representation(status, REPR_MMD_I2C                   &
                   , 'REPR_MMD_I2C_'//zchar//'LEV'//nchar//'N'               &
                   , rank = 4, link = 'xxxx', dctype = DC_I2C                &
                    , dimension_ids = (/&
                    _RI_XYZN_(DIMID_I2C_IE,DIMID_I2C_JE,DIMID_Z,DIMID_N)&
                     /) &
                   , ldimlen       = (/ &
                      _RI_XYZN_(INT(ie2lm), INT(je2lm), AUTO, AUTO)/)         &
                   , output_order  = (/&
                       _IX_XYZN_,_IY_XYZN_, _IZ_XYZN_, _IN_XYZN_ /)           &
                   , axis          = repr_def_axes(_RI_XYZN_('X','Y','Z','N')) &
                   )
              CALL channel_halt(substr, status)

              start(:,:) = gp_start(:,:)
              cnt(:,_IX_XYZN_)   = ie2lm
              cnt(:,_IY_XYZN_)   = je2lm
              cnt(:,_IN_XYZN_)   = nsize
              cnt(:,_IZ_XYZN_)   = zsize
              meml(:,:)  = gp_meml(:,:)
              memu(:,_IX_XYZN_)  = ie2lm
              memu(:,_IY_XYZN_)  = je2lm
              memu(:,_IN_XYZN_)  = nsize
              memu(:,_IZ_XYZN_)  = zsize

              CALL set_representation_decomp(status, REPR_MMD_I2C &
                   , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
              CALL channel_halt(substr, status)

           ENDIF
        ELSE
           ! rank=3
           IF (PRESENT(vsize)) THEN
              IF (vsize==0) THEN
                 size = ke_in
                 char = ''
                 DIMID=DIMID_I2C_KE_IN
              ELSE
                 size=zsize
                 char=zchar
                 DIMID=DIMID_Z
              ENDIF
              tag='LEV'
           ELSE
              !nsize is present
              size=nsize
              char=nchar
              DIMID=DIMID_N
              tag='N'
           ENDIF


           CALL get_representation_info(status                     &
                , 'REPR_MMD_IN_'//TRIM(char)//TRIM(tag), id= REPR_MMD_IN)

           IF (status /= 0) THEN

              ! NEW REPRESENTATIONS
              CALL new_representation(status, REPR_MMD_IN                    &
                   , 'REPR_MMD_IN_'//TRIM(char)//TRIM(tag)                   &
                   , rank = 3, link = 'xxx-', dctype = DC_I2C_IN             &
                   , dimension_ids = (/&
                      _RI_XYZ__(DIMID_I2C_IE_IN, DIMID_I2C_JE_IN, DIMID) /)  &
                   , ldimlen       = (/ &
                     _RI_XYZ__(INT(ie_in), INT(je_in), AUTO)/)         &
                   , output_order  = (/ _IX_XYZ__,_IY_XYZ__,_IZ_XYZ__ /)     &
                   , axis = repr_def_axes(_RI_XYZ__('X','Y',tag),'-')        &
                   )
              CALL channel_halt(substr, status)

              start(:,:) = gp_start(:,:)
              cnt(:,_IX_XYZ__)   = ie_in
              cnt(:,_IY_XYZ__)   = je_in
              cnt(:,_IZ_XYZ__)   = size
              cnt(:,4)   = 1
              meml(:,:)  = gp_meml(:,:)
              memu(:,_IX_XYZ__)  = ie_in
              memu(:,_IY_XYZ__)  = je_in
              memu(:,_IZ_XYZ__)  = size
              memu(:,4)  = 1

              CALL set_representation_decomp(status, REPR_MMD_IN &
                   , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
              CALL channel_halt(substr, status)

           ENDIF
           !------------------------------------------------------------

           IF (PRESENT(vsize)) THEN
              IF (vsize==0) THEN
                 size = kedim
                 char = ''
                 DIMID=DIMID_I2C_KEDIM
              ENDIF
           ENDIF
           CALL get_representation_info(status                     &
                , 'REPR_MMD_I2C_'//char//TRIM(tag), id= REPR_MMD_I2c)

           IF (status /= 0) THEN
              CALL new_representation(status, REPR_MMD_I2C                   &
                   , 'REPR_MMD_I2C_'//char//TRIM(tag)                        &
                   , rank = 3, link = 'xxx-', dctype = DC_I2C                &
                   , dimension_ids = (/&
                      _RI_XYZ__(DIMID_I2C_IE, DIMID_I2C_JE, DIMID) /)  &
                   , ldimlen       = (/ &
                     _RI_XYZ__(INT(ie2lm), INT(je2lm), AUTO)/)         &
                   , output_order  = (/ _IX_XYZ__,_IY_XYZ__,_IZ_XYZ__ /)     &
                   , axis = repr_def_axes(_RI_XYZ__('X','Y',tag),'-')        &
                   )
              CALL channel_halt(substr, status)

              start(:,:) = gp_start(:,:)
              cnt(:,_IX_XYZ__) = ie2lm
              cnt(:,_IY_XYZ__) = je2lm
              cnt(:,_IZ_XYZ__) = size
              cnt(:,4) = 1
              meml(:,:) = gp_meml(:,:)
              memu(:,_IX_XYZ__) = ie2lm
              memu(:,_IY_XYZ__) = je2lm
              memu(:,_IZ_XYZ__) = size
              memu(:,4) = 1

              CALL set_representation_decomp(status, REPR_MMD_I2C &
                   , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
              CALL channel_halt(substr, status)
           ENDIF
        ENDIF

        !----- DEALLOCATE------------------------------------------
        DEALLOCATE(start) ; NULLIFY(start)
        DEALLOCATE(cnt)   ; NULLIFY(cnt)
        DEALLOCATE(meml)  ; NULLIFY(meml)
        DEALLOCATE(memu)  ; NULLIFY(memu)
        ! ---------------------------------------------------------
      END SUBROUTINE make_i2c_representations
      !------------------------------------------------------------

      !------------------------------------------------------------
      SUBROUTINE make_cosmo_representation (rank, vsize, nsize &
           , dim_ids, dim_axis, dim_len, reprid)

        IMPLICIT NONE

        INTEGER, INTENT(IN)  :: rank
        INTEGER, INTENT(IN)  :: vsize
        INTEGER, INTENT(IN)  :: nsize
        INTEGER, DIMENSION(4), INTENT(IN) :: dim_ids
        CHARACTER(LEN=4),      INTENT(IN) :: dim_axis
        INTEGER, DIMENSION(4), INTENT(IN) :: dim_len

        INTEGER, INTENT(OUT) :: reprid

        ! PARALLEL DECOMPOSITION
        INTEGER                          :: nseg = 0
        INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
        INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
        INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
        INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

        CHARACTER(LEN=2) :: vchar, nchar

        nseg = gp_nseg

        ALLOCATE(start(nseg,IRANK))
        ALLOCATE(cnt(nseg,IRANK))
        ALLOCATE(meml(nseg,IRANK))
        ALLOCATE(memu(nseg,IRANK))

        start(:,:)=1
        cnt(:,:)  =1
        meml(:,:) =1
        memu(:,:) =1

        CALL int2str(vchar,vsize)
        CALL int2str(nchar,nsize)

        IF (rank == 4) THEN

           CALL get_representation_info(status                    &
                ,  'REPR_4D_'//vchar//'LEV_'//nchar//'N', id=reprid)

           IF (status /= 0) THEN
              ! NEW REPRESENTATIONS
              CALL new_representation(status, reprid             &
                   , 'REPR_4D_'//vchar//'LEV_'//nchar//'N'       &
                   , rank = rank, link = 'xxxx', dctype = DC_GP  &
                   , dimension_ids = dim_ids                     &
                   , ldimlen     = dim_len                       &
                   , output_order= (/ _IX_XYZN_,_IY_XYZN_,_IZ_XYZN_,_IN_XYZN_ /)&
                   , axis        = dim_axis                      &
                   )
              CALL channel_halt(substr, status)

              start(:,:) = gp_start(:,:)
              cnt(:,_IX_XYZN_) = ie
              cnt(:,_IY_XYZN_) = je
              cnt(:,_IN_XYZN_) = nsize
              cnt(:,_IZ_XYZN_) = vsize
              meml(:,:) = gp_meml(:,:)
              memu(:,_IX_XYZN_) = ie
              memu(:,_IY_XYZN_) = je
              memu(:,_IN_XYZN_) = nsize
              memu(:,_IZ_XYZN_) = vsize

              CALL set_representation_decomp(status, reprid      &
                   , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
              CALL channel_halt(substr, status)

           ENDIF
        !------------------------------------------------------------

     ELSE ! rank =3

        IF (vsize > 0) THEN
           CALL get_representation_info(status                 &
                , 'REPR_3D_'//vchar//'LEV', id=reprid          )
           IF (status /= 0) THEN
              CALL new_representation(status, reprid &
                   , 'REPR_3D_'//vchar//'LEV' &
                   , rank = 3, link = 'xxx-', dctype = DC_GP              &
                   , dimension_ids = dim_ids, ldimlen= dim_len            &
                   , output_order  = (/ _IX_XYZ__, _IY_XYZ__, _IZ_XYZ__/) &
                   , axis      = dim_axis )
              CALL channel_halt(substr, status)

              start(:,:) = gp_start(:,:)
              cnt(:,_IX_XYZ__)   = ie
              cnt(:,_IY_XYZ__)   = je
              cnt(:,_IZ_XYZ__)   = vsize
              cnt(:,4)   = 1
              meml(:,:)  = gp_meml(:,:)
              memu(:,_IX_XYZ__)  = ie
              memu(:,_IY_XYZ__)  = je
              memu(:,_IZ_XYZ__)  = vsize
              memu(:,4)  = 1

              CALL set_representation_decomp(status, reprid      &
                   , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
              CALL channel_halt(substr, status)
           ENDIF

        ELSE IF (nsize > 0) THEN
           CALL get_representation_info(status        &
                , 'REPR_3D_'//nchar//'N', id = reprid )

           IF (status /= 0 ) THEN
              CALL new_representation(status, reprid, 'REPR_3D_'//nchar//'N' &
                   , rank = 3, link = 'xxx-',      dctype   = DC_GP          &
                   , dimension_ids  = dim_ids,     ldimlen  = dim_len        &
                   , output_order  = (/ _IX_XY_N_, _IY_XY_N_, _IN_XY_N_/)    &
                   , axis     = dim_axis )
              CALL channel_halt(substr, status)

              start(:,:) = gp_start(:,:)
              cnt(:,_IX_XY_N_)   = ie
              cnt(:,_IY_XY_N_)   = je
              cnt(:,_IN_XY_N_)   = nsize
              cnt(:,4)   = 1
              meml(:,:)  = gp_meml(:,:)
              memu(:,_IX_XY_N_)  = ie
              memu(:,_IY_XY_N_)  = je
              memu(:,_IN_XY_N_)  = nsize
              memu(:,4)  = 1

              CALL set_representation_decomp(status, reprid      &
                   , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
              CALL channel_halt(substr, status)
           ENDIF

        ENDIF
     ENDIF
     !----- DEALLOCATE------------------------------------------
     DEALLOCATE(start) ; NULLIFY(start)
     DEALLOCATE(cnt)   ; NULLIFY(cnt)
     DEALLOCATE(meml)  ; NULLIFY(meml)
     DEALLOCATE(memu)  ; NULLIFY(memu)
     ! ---------------------------------------------------------
   END SUBROUTINE make_cosmo_representation
   ! --------------------------------------------------------------------

 END SUBROUTINE mmd2way_child_set_CplData
 ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE Define_data_arrays

    ! MMD
    USE mmd_child,  ONLY: MMD_C_GetNextArray  &
                         , MMD_C_Set_DataArray

    IMPLICIT NONE

    ! local variables
    INTEGER                         :: ii
    CHARACTER(LEN=*), PARAMETER     :: &
                             substr = 'mmd2way_child:Define_data_arrays'
    CHARACTER(LEN=STRLEN_CHANNEL)   :: myChannel
    CHARACTER(LEN=STRLEN_OBJECT)    :: myName
    INTEGER                         :: status

    ii = 0
    DO WHILE (MMD_C_GetNextArray (myChannel, myName))

       ii=ii+1
       CALL MMD_C_Set_DataArray(status, CplData(ii)%ldimlen &
            , CplData(ii)%axis, p4=CplData(ii)%ptr_in)
       IF (status /=0)  &
            CALL error_bi('ERROR in MMD_C_Set_DataArray',substr)
    ENDDO

    RETURN

  END SUBROUTINE Define_data_arrays
  ! ---------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE exchange_interpol_data

    ! MESSY
    USE messy_main_mpi_bi,        ONLY: switch_par_utilities
    USE messy_main_timer,         ONLY: lstart, lfirst_cycle
    ! MMD
    USE mmd_child,                ONLY: MMD_C_GetBuffer, MMD_STATUS_OK
    USE mmd_test,                 ONLY: MMD_testC_Compare

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='mmd2way_child_exchg_interpol_data'
    INTEGER                     :: status

    ! START CplData EXCHANGE
    CALL debug_bi(submodstr, 'START DATA EXCHANGE', '', substr)

    ! Get Parent / Coarse Grid DATA
    CALL main_qtimer_measure(i_bufC, MEAS_ON)
    CALL MMD_C_GetBuffer(WaitTime = WaitTime)
    CALL main_qtimer_measure(i_bufC, MEAS_OFF)

    IF (lfirst_cycle .AND. l_test) THEN
       CALL MMD_testC_Compare (my_lon,my_lat,status)
       IF(status /= MMD_STATUS_OK)   THEN
          CALL error_bi ('ERROR : MMD_test_Compare','mmd2way_child')
       ENDIF
    ENDIF

    ! switch ON INT2COSMO parallelisation
    CALL switch_par_utilities (1)

    ! READ IN THE EXTERNAL PARAMETER AND SET SOME FLAGS
    CALL mmd2way_prepare_external_data

    ! INTERPOLATION OF INT2COSMO ARRAYs
    CALL mmd2way_child_interpolation

    ! INTERPOLATION OF ADDITIONAL ARRAYs
    CALL Interpol_AddiArrays

    IF (lstart .AND. l_forcevars) CALL force_vars_from_file

    ! switch OFF INT2COSMO parallelisation
    CALL switch_par_utilities (2)

    ! COPY RESULTS TO COSMO VARIABLES
    ! boundary arrays only
    ! a) INITIAL DATA
    IF (lstart) CALL move_initial_arrays_to_COSMO

    ! b) BOUNDARY DATA
    CALL move_boundary_arrays_to_COSMO

    CALL debug_bi(submodstr, 'STOP DATA EXCHANGE', '', substr)

  END SUBROUTINE exchange_interpol_data
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE mmd2way_prepare_external_data

    ! MESSy/BMIL
    USE messy_main_timer_bi,   ONLY: event_state
    USE messy_main_mpi_bi,     ONLY: p_pe, switch_par_utilities &
                                   , gather_field
    USE messy_main_blather_bi, ONLY: error_bi
    ! MESSy/SMCL
    USE messy_main_timer,      ONLY: current_date, lfirst_cycle
    USE messy_main_tools,      ONLY: find_next_free_unit

    ! INT2COSMO
    USE data_int2lm_control,   ONLY: noutput, lgme2lm
    USE src_read_coarse_grid,  ONLY: org_read_coarse_grid
    USE data_fields_in,        ONLY: fr_land_in, fland_in_tot
    USE data_grid_in,          ONLY: ie_in, je_in, ie_in_tot, je_in_tot
    USE src_lm_output,         ONLY: init_lm_output


    IMPLICIT NONE

    EXTERNAL  :: external_data

    CHARACTER(LEN=80)          :: yerror
    CHARACTER(LEN=*),PARAMETER :: substr="mmd2way_prepare_external_data"
    INTEGER                    :: status
    INTEGER(I4)                :: stat_i4
    LOGICAL, SAVE              :: lread  = .TRUE.
    LOGICAL, SAVE              :: lfirst = .TRUE.

    lread = event_state(READEXT_EVENT, current_date) .OR. lread

    IF (lfirst_cycle) THEN
       ! Open file OUTPUT for compute PEs
       IF (p_pe == 0) THEN
          ! get again a unit number for OUTPUT
          noutput = find_next_free_unit(100,200)

          OPEN(noutput, FILE='OUTPUT', FORM=  'FORMATTED', STATUS='OLD',      &
               POSITION='APPEND', IOSTAT=status)
          IF(status /= 0) THEN
             yerror  = ' ERROR    *** Error while opening file OUTPUT *** '
             CALL error_bi(substr, yerror )
          ENDIF
       ENDIF
    ENDIF

    IF (.NOT. lgme2lm) THEN
       CALL switch_par_utilities(5)
       CALL gather_field( fr_land_in,   ie_in,     je_in                &
                        , fland_in_tot, ie_in_tot, je_in_tot, -1_I4, stat_I4)
       IF (stat_i4 /= 0) CALL error_bi('ERROR: gather fland_in_tot', 'mmd2way')
       CALL switch_par_utilities(1)
    ENDIF

    ! Read external data after receiving initial coarse data from parent
    CALL external_data (lread, status, yerror)
    IF (status /= 0) THEN
       write(iouerr,*) status, 'ERROR in mmd2way child external_data: ', yerror
       CALL error_bi('ERROR : external_data','mmd2way_child')
    ENDIF

    IF (lread .and. l_I2Cori_output) CALL init_lm_output
    lread = .FALSE.

    CALL org_read_coarse_grid(lfirst)

  END SUBROUTINE mmd2way_prepare_external_data
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE mmd2way_child_interpolation

    ! INT2COSMO
    USE data_int2lm_control,      ONLY: llm2lm, lcomp_bound
    USE data_fields_lm,           ONLY: grh_lm
    USE data_grid_lm,             ONLY: kedim, ie2lm,je2lm
    USE src_coarse_interpol,      ONLY: org_coarse_interpol
    USE src_vert_interpol,        ONLY: org_vert_interpol
    USE src_vert_inter_lm,        ONLY: org_vert_inter_lm
    USE src_2d_fields,            ONLY: org_2d_fields
    USE src_lm_fields,            ONLY: org_lm_fields
    USE src_lm_output,            ONLY: org_lm_output
    USE data_int2lm_io,           ONLY: numlist_ini,  numlist_bd  &
                                      , youtlist_ini, youtlist_bd
    USE data_runcontrol,          ONLY: ntstep, yakdat1
    ! MESSy/BMIL
    USE messy_main_timer_bi,      ONLY: event_state
    ! MESSy/SMCL
    USE messy_main_timer,         ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND &
                                      , DAYOFYEAR, current_date, lstart

    IMPLICIT NONE

    INTRINSIC :: INT, REAL

    ! LOCAL
    CHARACTER(LEN=14)          :: cdate
    LOGICAL                    :: lout

    ! ALLOCATE local arrays  needed for vertical interpolation
    ALLOCATE(grh_lm(ie2lm,je2lm,kedim))
    grh_lm = 0._dp

    ! Horizontal interpolation
    CALL org_coarse_interpol

    write(cdate,'(i4.4,5i2.2)') YEAR,MONTH,DAY,HOUR, MINUTE, SECOND

    ! Vertical interpolation, Do the interpolation in two steps
    IF (llm2lm) THEN
       CALL org_vert_inter_lm
       CALL org_2d_fields (cdate,INT(DAYOFYEAR,I4),REAL(HOUR,dp))
    ELSE
       CALL org_vert_interpol
       CALL org_2d_fields (cdate,INT(DAYOFYEAR,I4),REAL(HOUR,dp))
       CALL org_lm_fields
    ENDIF

    IF (l_I2Cori_output) THEN

       lout = event_state(WRITEI2C_EVENT, current_date) .OR. lstart
       IF (lout) THEN

          IF (.NOT. lcomp_bound) THEN
             CALL org_lm_output(numlist_ini,youtlist_ini,ntstep &
                  , yakdat1 &
#ifdef I2COLD
                  , .TRUE. &
#endif
                  )
          ELSE
             CALL org_lm_output(numlist_bd, youtlist_bd, ntstep &
                  , yakdat1 &
#ifdef I2COLD
                  , .FALSE. &
#endif
                  )
          ENDIF
       END IF
    END IF

    DEALLOCATE(grh_lm)

  END SUBROUTINE mmd2way_child_interpolation
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE Interpol_AddiArrays

    ! MESSy/SMCL
    USE messy_main_grid,            ONLY: GRID_ERROR, INIT_GEOHYBGRID
    USE messy_main_grid_netcdf,     ONLY: t_ncvar, INIT_NCVAR
    USE messy_main_grid_tools,      ONLY: RGTOOL_CONVERT_DAT2VAR, RGTOOL_CONVERT
    USE messy_main_grid_trafo,      ONLY: RG_INT
    USE messy_main_grid_trafo_scrp, ONLY: SCRIP_CONTROL
    ! INT2COSMO
    USE data_grid_in,               ONLY: ke_in_tot

    IMPLICIT NONE

    INTRINSIC :: SIZE, TRIM

    ! Local Data
    CHARACTER(LEN=*), PARAMETER :: substr = 'Interpol_AddiArrays'
    CHARACTER(LEN=80)           :: yerrmsg = ''
    INTEGER(I4)                 :: stat_I4
    INTEGER                     :: my_lev
    INTEGER                     :: ii, k, iN

    TYPE(t_geohybgrid)                     :: intgrid
    TYPE(t_ncvar), DIMENSION(:),   POINTER :: vari => NULL()
    TYPE(t_ncvar), DIMENSION(:),   POINTER :: varo => NULL()
    REAL(dp), DIMENSION(:,:,:,:),  POINTER :: dat  => NULL()
    INTEGER                                :: ix
    INTEGER                                :: status
    INTEGER                                :: RGT(1)

    DO ii = 1 , NEXCH
       ! IF  ARRAY IS NO ADDITIONAL ARRAY CYCLE
       IF (CplData(ii)%lvartab) CYCLE
       IF (TRIM(CplData(ii)%CHILD%OBJ) == 'Test_Ar') CYCLE

       IF (SIZE(CplData(ii)%ptr_in,4) == 1) THEN ! 2D/3D-ARRAY
          ! Horizontal interpolation
          ! CplData(ii)%C_INTERPOL => Flag for interpolation type
          my_lev = SIZE(CplData(ii)%ptr_in,_IZ_XYZ__)

          IF (CplData(ii)%C_INTERPOL(1:1) /= 'C') THEN
             DO k=1,my_lev
                CALL interpol_coarse_OneLayer(&
                     CplData(ii)%ptr_in(_RI_XYZ__(:,:,k),1)     &
                   , CplData(ii)%ptr_i2c(_RI_XYZ__(:,:,k),1)    &
                   , CplData(ii)%C_INTERPOL(1:3) &
                   , yerrmsg, stat_I4, CplData(ii)%CHILD%OBJ)
                IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)
             ENDDO
          ELSE
             ! USE Conservative remapping via SCRIP
             ALLOCATE(vari(1))

             IF (my_lev == 1) THEN
                ! 2D-ARRAY
                CALL RGTOOL_CONVERT_DAT2VAR(vari(1), CplData(ii)%ptr_in &
                     , CplData(ii)%CHILD%OBJ, pinhgrid, 'xyzn')
             ELSE
                 ! 3D-ARRAY
                CALL RGTOOL_CONVERT_DAT2VAR(vari(1), CplData(ii)%ptr_in &
                     , CplData(ii)%CHILD%OBJ, pingrid, 'xyzn')
             ENDIF

             RGT(1) = RG_INT
             CALL SCRIP_CONTROL(status, I2C_SD_ID, pingrid, i2cgrid &
                  , RGT, .TRUE., vari, varo, intgrid, llrgz=.FALSE.)
             IF (status /= 0) CALL error_bi(grid_error(status), substr)

             CALL RGTOOL_CONVERT(varo(1), dat, i2chgrid, 'xyzn')
             CplData(ii)%ptr_i2c(:,:,:,1) = dat(:,:,1,:)
             DEALLOCATE(dat)
             NULLIFY(dat)

             DO ix= 1, SIZE(vari)
                CALL INIT_NCVAR(vari(ix))
                CALL INIT_NCVAR(varo(ix))
             END DO
             DEALLOCATE(vari)
             NULLIFY(vari)
             DEALLOCATE(varo)
             NULLIFY(varo)
             CALL INIT_GEOHYBGRID(intgrid)
          END IF

          ! Vertical interpolation an lm_field
          ! Vertical interpolation Only if NLEV fields
          IF( my_lev == ke_in_tot) THEN
             IF (CplData(ii)%C_INTERPOL(4:4) =='V') THEN
                CALL interpol_vert_AddiArray (CplData(ii)%ptr_i2c(:,:,:,1) &
                     , CplData(ii)%CHILD%OBJ,CplData(ii)%C_INTERPOL(3:3))
             ELSE IF (CplData(ii)%C_INTERPOL(4:4) =='W') THEN
                CALL vert_ncinterpol(CplData(ii)%CHILD%OBJ   &
                     , CplData(ii)%ptr_i2c                   &
                     , psgl4d, pslm4d, pingrid, i2cgrid)
             END IF
          ENDIF
       ELSE ! 4D-ARRAY
          my_lev = SIZE(CplData(ii)%ptr_in,_IZ_XYZN_)

          IF (CplData(ii)%C_INTERPOL(1:1) /= 'C') THEN

             DO iN = 1, SIZE(CplData(ii)%ptr_in,_IN_XYZN_)
                DO k=1,my_lev
                   CALL interpol_coarse_OneLayer(&
                        CplData(ii)%ptr_in(_RI_XYZN_(:,:,k,iN))          &
                        , CplData(ii)%ptr_i2c(_RI_XYZN_(:,:,k,iN))       &
                        , CplData(ii)%C_INTERPOL(1:3), yerrmsg, stat_I4  &
                        , CplData(ii)%CHILD%OBJ)
                   IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)
                ENDDO
             END DO
          ELSE

             ! USE Conservative remapping via SCRIP
             ALLOCATE(vari(1))

             ! 3D/4D-ARRAY
             CALL RGTOOL_CONVERT_DAT2VAR(vari(1), CplData(ii)%ptr_in &
                  , CplData(ii)%CHILD%OBJ, pingrid, 'xyzn')

             RGT(1) = RG_INT
             CALL SCRIP_CONTROL(status, I2C_SD_ID, pingrid, i2cgrid &
                  , RGT, .TRUE., vari, varo, llrgz=.FALSE.)
             IF (status /= 0) CALL error_bi(grid_error(status), substr)

             CALL RGTOOL_CONVERT(varo(1), dat, i2chgrid, 'xyzn')
             CplData(ii)%ptr_i2c(:,:,:,:) = dat(:,:,:,:)
             DEALLOCATE(dat)
             NULLIFY(dat)

             DO ix= 1, SIZE(vari)
                CALL INIT_NCVAR(vari(ix))
                CALL INIT_NCVAR(varo(ix))
             END DO
             DEALLOCATE(vari)
             NULLIFY(vari)
             DEALLOCATE(varo)
             NULLIFY(varo)
             CALL INIT_GEOHYBGRID(intgrid)
          END IF
          IF (my_lev == ke_in_tot) THEN
             IF (CplData(ii)%C_INTERPOL(4:4) =='V')  THEN
                DO iN = 1, SIZE(CplData(ii)%ptr_in,_IN_XYZN_)
                   ! Vertical interpolation Only if NLEV fields
                   CALL interpol_vert_AddiArray( &
                        CplData(ii)%ptr_i2c(_RI_XYZN_(:,:,:,iN)) &
                        , CplData(ii)%CHILD%OBJ,CplData(ii)%C_INTERPOL(3:3) )
                END DO
             ELSE IF (CplData(ii)%C_INTERPOL(4:4) =='W')  THEN
                CALL vert_ncinterpol(CplData(ii)%CHILD%OBJ   &
                     , CplData(ii)%ptr_i2c &
                     , psgl4d, pslm4d, pingrid, i2cgrid)
             ENDIF
          END IF
       ENDIF
    ENDDO

  END SUBROUTINE Interpol_AddiArrays
  !-------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE interpol_coarse_OneLayer (field_in, field_lm, yzitype &
                                      , yerrmsg, stat_I4, name)

    ! INT2COSMO
    USE data_int2lm_control, ONLY: llm2lm, lcm2lm, lec2lm, lum2lm,lbd_frame_cur &
                                 , l_bicub_spl, l_cressman
    USE data_int2lm_io,      ONLY: undef
    USE data_fields_lm,      ONLY: i_index, j_index, x_wght, y_wght             &
                                 , lmask_lm, lolp_lm
    USE data_grid_lm,        ONLY: ie2lm, je2lm
    USE data_grid_in,        ONLY: grdpt_rel_in, ie_in_tot, je_in_tot           &
                                 , startlat_in, startlon_in, ie_in, je_in       &
                                 , latitudes_in,longitudes_in
    USE data_fields_in,      ONLY: lolp_in, lat_coarse_m, lon_coarse_m
    USE interp_utilities,    ONLY: interp_l, interp_q, interp_q_lm, interp_q_bs

    IMPLICIT NONE

    INTRINSIC :: TRIM

    REAL(KIND=dp),DIMENSION(:,:),INTENT(IN)   :: field_in
    REAL(KIND=dp),DIMENSION(:,:),INTENT(OUT)  :: field_lm
    CHARACTER(LEN=3), INTENT(IN)              :: yzitype
    CHARACTER(LEN=*), INTENT(IN)              :: name
    CHARACTER(LEN=*), INTENT(OUT)             :: yerrmsg
    INTEGER(I4),      INTENT(OUT)             :: stat_I4

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='interpol_coarse_OneLayer'
    LOGICAL                     :: lzmono
    LOGICAL                     :: lzposdef

    grdpt_rel_in = 0
    lzmono   = .FALSE.
    lzposdef = .FALSE.
    IF (yzitype(2:2) == 'T') lzmono   = .TRUE.
    IF (yzitype(3:3) == 'T') lzposdef = .TRUE.

    IF (yzitype(1:1) == 'Q') THEN
       ! quadratic interpolation to LM-field
       IF (lec2lm .OR. lcm2lm) THEN
          ! Introduce option of bicubic spline-interpolation instead of quadratic
          IF (l_bicub_spl) THEN
             CALL interp_q_bs                                                  &
                  (field_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),     &
                  lzmono, lzposdef, lbd_frame_cur, lolp_in, lolp_lm,           &
                  undef, lmask_lm, x_wght(:,:,1), y_wght(:,:,1),               &
                  field_lm, 1_I4, ie2lm, 1_I4, je2lm, startlat_in, startlon_in,&
                  latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,     &
                  grdpt_rel_in, ie_in_tot, je_in_tot, lcm2lm, l_cressman,      &
                  yerrmsg, stat_I4)
             IF (stat_I4 /= 0) CALL error_bi(yerrmsg,substr)
          ELSE
             CALL interp_q(field_in, ie_in,je_in, i_index(:,:,1),j_index(:,:,1),&
                  lzmono, lzposdef, lbd_frame_cur, undef, lmask_lm,             &
                  x_wght(:,:,1), y_wght(:,:,1),                                 &
                  field_lm, 1_I4, ie2lm, 1_I4, je2lm, yerrmsg, stat_I4)
             IF (stat_I4 /= 0) CALL error_bi(yerrmsg,substr)
          ENDIF
       ELSEIF (llm2lm .OR. lum2lm) THEN
          CALL interp_q_lm(field_in, ie_in, je_in,          &
               field_lm, i_index(:,:,1), j_index(:,:,1),    &
               x_wght(:,:,1), y_wght(:,:,1), ie2lm, je2lm,  &
               lzmono, lzposdef, yerrmsg, stat_I4)
          IF (stat_I4 /= 0) CALL error_bi(yerrmsg,substr)
       ENDIF
    ELSE IF (yzitype(1:1)=='M') THEN
       ! normal linear (or match) interpolation to LM-field
       IF (l_cressman) grdpt_rel_in=1
       CALL interp_l(field_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
            lzmono, lzposdef, 'M', lbd_frame_cur, lolp_in,                   &
            lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),          &
            field_lm, 1_I4, ie2lm, 1_I4, je2lm, startlat_in, startlon_in,    &
            latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,         &
            grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                  &
            yerrmsg, stat_I4)
       IF (stat_I4 /= 0) CALL error_bi(yerrmsg,substr)
    ELSE IF (yzitype(1:1)=='N') THEN
       ! normal linear (or match) interpolation to LM-field
       CALL interp_l(field_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
            lzmono, lzposdef, 'N', lbd_frame_cur, lolp_in,                   &
            lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),          &
            field_lm, 1_I4, ie2lm, 1_I4, je2lm, startlat_in, startlon_in,    &
            latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,         &
            grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                  &
            yerrmsg, stat_I4)
       IF (stat_I4 /= 0) CALL error_bi(yerrmsg,substr)
    ELSE IF (yzitype(1:1)=='L') THEN
       CALL interp_l(field_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
              .FALSE., .FALSE., 'L', lbd_frame_cur, lolp_in,                 &
              lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),        &
              field_lm, 1_I4, ie2lm, 1_I4, je2lm, startlat_in, startlon_in,  &
              latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,       &
              grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                &
              yerrmsg, stat_I4)
       IF (stat_I4 /= 0) CALL error_bi(yerrmsg,substr)
    ELSE
       CALL error_bi(substr, &
            'INTERPOLATION METHOD '//yzitype(1:1)//' FOR '//&
            &TRIM(name)//' IS NOT IMPLEMENTED')
    ENDIF

    RETURN

  END SUBROUTINE interpol_coarse_OneLayer
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE interpol_vert_AddiArray (field_lm, varname, cposdef)

    ! INT2COSMO
    USE data_int2lm_control, ONLY: lcm2lm
    USE data_fields_lm,      ONLY: zps1_lm, zkzgr, zfi_fl, zhi_fl
    USE data_grid_lm,        ONLY: kedim
    USE vgrid_refatm_utils,  ONLY: vcoord
    USE data_grid_in,        ONLY: ke_in, akh_in, bkh_in
    USE src_vert_interpol,   ONLY: vert_interpol
    USE src_lm_fields,       ONLY: vert_int_lm, vert_z_lm
    USE src_vert_inter_lm,   ONLY: vert_interp

    IMPLICIT NONE

    INTRINSIC :: INT, MAX

    CHARACTER(LEN=STRLEN_OBJECT)                  :: varname
    REAL(KIND=DP),DIMENSION(:,:,:),INTENT(INOUT)  :: field_lm
    CHARACTER(LEN=1), INTENT(IN)                  :: cposdef

    ! local variables
    INTEGER(I4)                                   :: izdebug = 0
#if !defined(I2COLD) && defined(I2CINC)
    REAL(KIND=DP)    :: zh_vertgrad(2) = (/500.0_dp, 1000.0_dp/)
#endif
    IF (lcm2lm) THEN   ! ECHAM is Parent
       CALL vert_interpol (field_lm, varname, _VI_SPLINE_A_ ke_in, akh_in, bkh_in &
            , zps1_lm, INT(zkzgr,I4))

      ! Force positive values (done hard-coded in vert_interpol for qv, qi etc.)
       IF (cposdef == 'T') field_lm = MAX(0._dp, field_lm)

       IF (vcoord%vcflat > 0.0_dp) THEN
          CALL vert_int_lm (field_lm, varname, _VI_SPLINE_B_  kedim, zfi_fl, ke_in, izdebug)
       ELSE
          CALL vert_z_lm   (field_lm, varname,  _VI_SPLINE_B_ kedim, zfi_fl, ke_in _I2C03bDEBUG_ )
       ENDIF
    ELSE                                    !COSMO is Parent
       CALL vert_interp (field_lm, varname,   _VI_SPLINE_A_ kedim, zhi_fl, ke_in &
#if !defined(I2COLD) && defined(I2CINC)
            , izdebug, INT(zkzgr,I4), zh_vertgrad(1), zh_vertgrad(2), &
            ypbl_method_for_shrink='linear-transition', &
            ypbl_method_for_valleyextrapol='gradient-lintrans' &
            , lpbl_slope_corr=.FALSE. &
#else
            , INT(zkzgr,I4), izdebug &
#endif
            )
    ENDIF

    IF (cposdef == 'T') field_lm = MAX(0._dp, field_lm)

    RETURN

  END SUBROUTINE Interpol_vert_AddiArray
  ! ----------------------------------------------------------------------

  !========================================================================

  SUBROUTINE force_vars_from_file

    ! This subroutine overwrites the INT2LM results by values read from
    ! a netcdf-file,
    ! NOTE: the dimensions in the netcdf file and that of the COSMO domain
    !       have to be the same !!!

    ! BMIL
    USE messy_main_mpi_bi,        ONLY: scatter_gp
    ! SMCL
    USE messy_main_tools,         ONLY: strcrack, find_next_free_unit
    USE messy_main_grid_netcdf
    ! INT2COSMO
    USE data_grid_lm,             ONLY: ie2lm_tot, je2lm_tot

    IMPLICIT NONE
    !
    CHARACTER(LEN=*), PARAMETER :: substr = 'force_vars_from_file'
    ! names of variable to be forced by offline field
    CHARACTER(LEN=STRLEN_OBJECT), DIMENSION(:), POINTER :: varname  => NULL()
    ! number of offline forced variables
    INTEGER         :: nvar

    ! Unit for input file
    INTEGER         :: iou
    LOGICAL         :: lex
    ! IMPORTED VARIABLE (1D !)
    TYPE(t_ncvar) :: fvar
    ! GLOBAL VARIABLE FIELD
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: globvar  => NULL()
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: gvarsort => NULL()
    REAL(dp), POINTER, DIMENSION(:,:,:)   :: gptr_3d  => NULL()
    REAL(dp), POINTER, DIMENSION(:,:,:)   :: lptr_3d  => NULL()
    INTEGER, DIMENSION(4)      :: dimlen
    INTEGER                    :: i,j            ! loop variable
    INTEGER                    :: idim           ! loop variable
    INTEGER                    :: j1, j2, j3, j4 ! loop variables
    INTEGER                    :: idx            ! index help variable

    IF (.NOT. l_forcevars) RETURN

    ! 1.) break variable string
    CALL strcrack(forcevars, ';' ,varname, nvar)

    ! 2) test of the given file is readable
    iou = find_next_free_unit(100,200)
    INQUIRE(file=TRIM(forcefile), exist=lex, number=iou)
    IF (.NOT. lex) THEN
       CALL error_bi('forcing file '''//TRIM(forcefile)//''' does not exist' &
            ,substr)
    ENDIF
    !
    forcevar_loop: DO i=1,nvar
       CplData_loop: DO j=1, NCOPY
          IF (TRIM(CplData(j)%CHILD%OBJ) == TRIM(varname(i))) THEN

         ! write(iouerr,*) ' IN ',substr,' processing variable: ' &
!             ,TRIM(varname(i)), CplData(j)%rank, CplData(j)%ltiles
             CALL debug_bi(modstr,'processing variable',varname(i),substr)

             ! ALLOCATE SPACE FOR "GLOBAL" VARIABLE
             IF (CplData(j)%rank == 2) THEN
                ALLOCATE(globvar(ie2lm_tot, je2lm_tot, 1, 1))
                ALLOCATE(gvarsort(ie2lm_tot-2, je2lm_tot-2, 1, 1))
             ELSE IF (CplData(j)%rank == 3  .OR. &
                  (CplData(j)%rank == 4 .AND. CplData(j)%ltiles)) THEN
                ! INT2LM does not know of tiles, therefore the input data
                ! of rank 4 variables in COSMO is still 3D
                ALLOCATE(globvar(ie2lm_tot, je2lm_tot &
                                   ,SIZE(CplData(j)%ptr_i2c,3),1))
                ALLOCATE(gvarsort(ie2lm_tot-2, je2lm_tot-2 &
                                   ,SIZE(CplData(j)%ptr_i2c,3),1))
             ELSE
                CALL error_bi( &
                     'offline forcing of 4D-fields not yet implemented' &
                     , substr)
             END IF
             globvar = -999._dp

             IF (p_parallel_io) THEN
                ! IMPORT VARIABLE
                CALL IMPORT_NCVAR(fvar, ustep=1, varname=varname(i) &
                     , file=TRIM(forcefile))
                dimlen(:) = 1
                DO idim = 1, fvar%ndims
                   dimlen(idim) = fvar%dim(idim)%len
                END DO
                IF (fvar%dat%type == VTYPE_DOUBLE) THEN
                   DO j4 = 1, dimlen(4)
                      DO j3 = 1, dimlen(3)
                         DO j2 = 1, dimlen(2)
                            DO j1 = 1, dimlen(1)
                               idx = j1 + (j2-1) * dimlen(1)           &
                                    + (j3-1) * dimlen(1)*dimlen(2) &
                                    + (j4-1) * dimlen(1)*dimlen(2)*dimlen(3)
                               gvarsort(j1,j2,j3,j4) = fvar%dat%vd(idx)
                            END DO
                         END DO
                      END DO
                   END DO
                ELSE IF (fvar%dat%type == VTYPE_REAL) THEN
                   DO j4 = 1, dimlen(4)
                      DO j3 = 1, dimlen(3)
                         DO j2 = 1, dimlen(2)
                            DO j1 = 1, dimlen(1)
                               idx = j1 + (j2-1) * dimlen(1)           &
                                    + (j3-1) * dimlen(1)*dimlen(2) &
                                    + (j4-1) * dimlen(1)*dimlen(2)*dimlen(3)
                               gvarsort(j1,j2,j3,j4) = REAL(fvar%dat%vr(idx),dp)
                            END DO
                         END DO
                      END DO
                   END DO
                ELSE
                   write(iouerr,*) 'data type ',fvar%dat%type,' not yet implemented'
                   CALL error_bi('data type  not yet implemented', substr)
                ENDIF ! DATA TYPE

                globvar(2:ie2lm_tot-1,2:je2lm_tot-1,:,:) = gvarsort(:,:,:,:)
                ! initialise boundary
                globvar(1,2:je2lm_tot-1,:,:) = gvarsort(1,:,:,:)
                globvar(ie2lm_tot,2:je2lm_tot-1,:,:)= &
                                               gvarsort(ie2lm_tot-2,:,:,:)
                globvar(:,1,:,:)             = globvar(:,2,:,:)
                globvar(:,je2lm_tot,:,:)     = globvar(:,je2lm_tot-1,:,:)

             END IF ! p_parallel_io

             IF (TRIM(varname(i)) /= 'T_SO') THEN
                CALL scatter_gp(globvar, CplData(j)%ptr_i2c)
             ELSE
                DO j4 = 1, SIZE(CplData(j)%ptr_i2c,3)-1
                   gptr_3d => globvar(:,:,j4,:)
                   lptr_3d => CplData(j)%ptr_i2c(:,:,j4+1,:)
                   CALL scatter_gp(gptr_3d, lptr_3d)
                END DO
             END IF

          END IF
          IF (ASSOCIATED(globvar))  DEALLOCATE(globvar)
          IF (ASSOCIATED(gvarsort)) DEALLOCATE(gvarsort)
       END DO CplData_loop
    END DO forcevar_loop

    IF (ASSOCIATED(varname)) DEALLOCATE(varname)

  END SUBROUTINE force_vars_from_file
  !========================================================================

  ! ---------------------------------------------------------------------
  SUBROUTINE flag_fields(itab, p4)

    ! INT2COSMO
    USE data_grid_lm,           ONLY: istartcos, iendcos        &
                                    , jstartcos, jendcos
    USE data_int2lm_io,         ONLY: var_lm, undefncdf
    USE data_fields_lm,         ONLY: lolp_lm, w_snow_lm

    ! MESSy/BMIL
    USE messy_main_grid_def_mem_bi, ONLY: istartpar_c4 => istartpar  &
                                    , iendpar_c4   => iendpar    &
                                    , jstartpar_c4 => jstartpar  &
                                    , jendpar_c4   => jendpar
    ! MESSY/SMCL
    USE messy_main_constants_mem, ONLY: g


    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM

    INTEGER                     , INTENT(IN)     :: itab
    REAL(DP), DIMENSION(:,:,:,:), INTENT(INOUT)  :: p4

    IF (ASSOCIATED( var_lm(itab)%p2) ) THEN
       IF (TRIM(var_lm(itab)%lsm) == 'l') THEN
          WHERE (.NOT. lolp_lm(istartcos:iendcos,jstartcos:jendcos))     &
               p4(istartpar_c4:iendpar_c4,jstartpar_c4:jendpar_c4,1,1) = &
                                                                     undefncdf
       ENDIF
       IF  (TRIM(var_lm(itab)%name) == 'T_SNOW' ) THEN
          WHERE (w_snow_lm(istartcos:iendcos,jstartcos:jendcos) <= 0.0_dp) &
               p4(istartpar_c4:iendpar_c4,jstartpar_c4:jendpar_c4,1,1) =   &
                                                                     undefncdf
       ENDIF

       ! Z0 is g * Z0 in COSMO
       ! Z0 is not needed for any calculations within INT2COSMO
       ! THUS THE ASSIGNMENT BELOW IS OK
       IF  (TRIM(var_lm(itab)%name) == 'Z0' ) THEN
          p4(istartpar_c4:iendpar_c4,jstartpar_c4:jendpar_c4,1,1) = &
               p4(istartpar_c4:iendpar_c4,jstartpar_c4:jendpar_c4,1,1) * g
       ENDIF
    ENDIF

  END SUBROUTINE flag_fields
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE move_initial_arrays_to_Cosmo

    ! NOTE: IN THIS SUBROUTINE TWO TYPE OF ARRAYS ARE MOVED TO COSMO i.e.:
    !       I) THE CONSTANT FIELDS READ IN WITH THE EXTERNAL FIELD
    !       II) THE INTERPOLATED FIELDS PROVIDED BY THE PARENT
    !    In CASE of a RESTART this subroutine should NOT be executed
    !    because all necessary data is included in the RESTART-FILES

    ! MESSy/BMIL
    USE messy_main_grid_def_mem_bi, ONLY: istartpar_c4 => istartpar  &
                                      , iendpar_c4   => iendpar    &
                                      , jstartpar_c4 => jstartpar  &
                                      , jendpar_c4   => jendpar    &
                                      , ie_c4        => ie         &
                                      , je_c4        => je
    USE messy_main_grid_def_bi,     ONLY:lperi_x, lperi_y, l2dim
    USE messy_main_data_bi,         ONLY: l2tls, nnow, nnew

    USE messy_main_mpi_bi,        ONLY: icomm_cart, sendbuf, isendbuflen      &
                                      , imp_reals, nboundlines, my_cart_neigh &
                                      , ncomm_type, ldatatypes, num_compute   &
                                      , exchg_boundaries
    ! INT2COSMO
    USE data_grid_lm,             ONLY: istartcos, iendcos, jstartcos, jendcos
    USE messy_main_constants_mem, ONLY: g

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, INT, SIZE, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='move_initial_arrays_to_Cosmo'
    INTEGER(I4)                 :: ii, ix
    INTEGER                     :: size3, size4
    INTEGER(I4)                 :: kdims(24)
    INTEGER(I4)                 :: stat_I4
    CHARACTER(LEN=80)           :: yerrmsg
    INTEGER(I4)                 :: icount

    icount = 1

    IF (.NOT. ALLOCATED(p_fis)) ALLOCATE (p_fis(ie_c4,je_c4,1,1))
    IF (MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN
       p_fis(istartpar_c4:iendpar_c4,jstartpar_c4:jendpar_c4,1,1) = &
            CplData(idx_cl_fis)%ptr_i2c(istartcos:iendcos,jstartcos:jendcos,1,1)
    ELSE
       ! hsurf multiply by g
       p_fis(istartpar_c4:iendpar_c4,jstartpar_c4:jendpar_c4,1,1) = g * &
            CplData(idx_cl_fis)%ptr_i2c(istartcos:iendcos,jstartcos:jendcos,1,1)
    END IF
    kdims(1:24) = &
         (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                             &
         (100+icount, sendbuf, isendbuflen, imp_reals, icomm_cart     &
         , num_compute, ie_c4, je_c4, kdims, jstartpar_c4, jendpar_c4 &
         , nboundlines, nboundlines, my_cart_neigh                    &
         , lperi_x, lperi_y, l2dim, 50000+icount                      &
         , ldatatypes, ncomm_type, stat_I4, yerrmsg                   &
         , p_fis(:,:,:,1))
    IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)
    icount=icount+1

    ! B) MOVE THE INTERPOLATED FIELDS
    DO ii = 1, NCOPY
       IF (.NOT. CplData(ii)%L_INITIAL) CYCLE
       IF (TRIM(CplData(ii)%CHILD%OBJ) == 'Test_Ar') CYCLE
       lasso: IF (ASSOCIATED(CplData(ii)%cosmo))  THEN
         ! CHECK IF COSMO 4D-POINTER IS TIME DEPENDENT
          ltime: IF (SIZE(CplData(ii)%cosmo) == 1) THEN
             ! NOT TIME DEPENDENT
             size3 = SIZE(CplData(ii)%cosmo(1)%ptr,3)
             size4 = SIZE(CplData(ii)%cosmo(1)%ptr,4)
             DO iX=1,size4
                CplData(ii)%cosmo(1)%ptr(istartpar_c4:iendpar_c4 &
                     ,jstartpar_c4:jendpar_c4,1:size3,iX) =      &
                     CplData(ii)%ptr_i2c(istartcos:iendcos       &
                     ,jstartcos:jendcos,1:size3,iX)
                IF (size4 == 1 .AND. &
                     CplData(ii)%vartab_idx /= 0)  CALL flag_fields(&
                     CplData(ii)%vartab_idx, CplData(ii)%cosmo(1)%ptr)
                kdims(1:24) = &
                     (/size3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
                CALL exchg_boundaries                                       &
                     (INT(100+icount,I4),sendbuf,isendbuflen, imp_reals     &
                     , icomm_cart, num_compute                              &
                     , ie_c4, je_c4, kdims, jstartpar_c4, jendpar_c4        &
                     , nboundlines, nboundlines, my_cart_neigh              &
                     , lperi_x, lperi_y, l2dim, 50000+icount                &
                     , ldatatypes, ncomm_type, stat_I4, yerrmsg             &
                     , CplData(ii)%cosmo(1)%ptr(:,:,:,IX))
                IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)
                icount=icount+1
             END DO
          ELSE
             ! POINTER IS TIME DEPENDENT
             if_trac_c:IF ( TRIM(CplData(ii)%CHILD%CHA) /= 'tracer_gp') THEN
                size3 = SIZE(CplData(ii)%cosmo(1)%ptr,3)
                size4 = SIZE(CplData(ii)%cosmo(1)%ptr,4)
                DO iX=1,size4
                   CplData(ii)%cosmo(nnew)%ptr(istartpar_c4:iendpar_c4 &
                        ,jstartpar_c4:jendpar_c4,1:size3,iX) =         &
                        CplData(ii)%ptr_i2c(istartcos:iendcos          &
                        ,jstartcos:jendcos,1:size3,iX)
                   IF (size4 == 1 .AND. &
                        CplData(ii)%vartab_idx /= 0)  CALL flag_fields (&
                        CplData(ii)%vartab_idx, CplData(ii)%cosmo(nnew)%ptr)
                   kdims(1:24) = &
                        (/size3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
                   CALL exchg_boundaries                                       &
                        (100+icount,sendbuf,isendbuflen,imp_reals,icomm_cart   &
                        , num_compute                                          &
                        , ie_c4, je_c4, kdims, jstartpar_c4, jendpar_c4        &
                        , nboundlines, nboundlines, my_cart_neigh              &
                        , lperi_x, lperi_y, l2dim, 50000+icount                &
                        , ldatatypes, ncomm_type, stat_I4, yerrmsg             &
                        , CplData(ii)%cosmo(nnew)%ptr(:,:,:,iX))
                   IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)
                   icount=icount+1
                END DO

                IF (.NOT. l2tls) THEN
                   DO iX=1,size4
                      CplData(ii)%cosmo(nnow)%ptr(istartpar_c4:iendpar_c4 &
                           ,jstartpar_c4:jendpar_c4,1:size3,iX) =         &
                           CplData(ii)%ptr_i2c(istartcos:iendcos          &
                           ,jstartcos:jendcos,1:size3,iX)
                      IF (size4 == 1 .AND. &
                           CplData(ii)%vartab_idx /= 0) CALL flag_fields (&
                           CplData(ii)%vartab_idx, CplData(ii)%cosmo(nnow)%ptr)
                      kdims(1:24) = &
                         (/size3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
                      CALL exchg_boundaries                                &
                           (100+icount, sendbuf, isendbuflen, imp_reals    &
                           , icomm_cart, num_compute, ie_c4, je_c4, kdims  &
                           , jstartpar_c4, jendpar_c4, nboundlines         &
                           , nboundlines, my_cart_neigh                    &
                           , lperi_x, lperi_y,  l2dim, 50000+icount        &
                           , ldatatypes, ncomm_type, stat_I4, yerrmsg     &
                           , CplData(ii)%cosmo(nnow)%ptr(:,:,:,iX))
                      IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)
                      icount=icount+1
                   END DO

                ENDIF
             ELSE ! if_trac_c
                size4 = SIZE(CplData(ii)%cosmo(1)%ptr,_IZ_XYZN_)
                DO iX = 1,size4
                   CplData(ii)%cosmo(1)%ptr(&
                        _RI_XYZN_(istartpar_c4:iendpar_c4,jstartpar_c4:jendpar_c4,iX,1)) =            &
                        CplData(ii)%ptr_i2c(_RI_XYZN_(istartcos:iendcos,jstartcos:jendcos,iX,1))
                ENDDO
                kdims(1:24) = &
                     (/size4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
                CALL exchg_boundaries                                         &
                     (100+icount, sendbuf, isendbuflen, imp_reals, icomm_cart &
                     , num_compute                                            &
                     , ie_c4, je_c4, kdims, jstartpar_c4, jendpar_c4          &
                     , nboundlines, nboundlines, my_cart_neigh                &
                     , lperi_x, lperi_y, l2dim, 50000+icount                  &
                     , ldatatypes, ncomm_type, stat_I4, yerrmsg               &
                     , CplData(ii)%cosmo(1)%ptr(_RI_XYZN_(:,:,:,1)))
                IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)
                icount=icount+1

                IF (.NOT. l2tls) THEN
                   DO iX =1,size4
                    CplData(ii)%cosmo(2)%ptr(&
                         _RI_XYZN_(istartpar_c4:iendpar_c4,jstartpar_c4:jendpar_c4,iX,1)) =            &
                           CplData(ii)%ptr_i2c(_RI_XYZN_(istartcos:iendcos,jstartcos:jendcos,iX,1))
                   ENDDO
                   CALL exchg_boundaries                                       &
                        (100+icount, sendbuf, isendbuflen, imp_reals           &
                        , icomm_cart, num_compute, ie_c4, je_c4, kdims         &
                        , jstartpar_c4, jendpar_c4, nboundlines, nboundlines   &
                        , my_cart_neigh, lperi_x, lperi_y, l2dim, 50000+icount &
                        , ldatatypes, ncomm_type, stat_I4, yerrmsg             &
                        , CplData(ii)%cosmo(2)%ptr(_RI_XYZN_(:,:,:,1)))
                   IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)
                   icount=icount+1
                ENDIF
             END IF if_trac_c
          END IF ltime

       ELSE ! not associated ptr_c

          CALL error_bi(substr, 'POINTER NOT ASSOCIATED')

       ENDIF lasso
   ENDDO

  END SUBROUTINE move_initial_arrays_to_Cosmo
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE move_boundary_arrays_to_Cosmo

    ! INT2COSMO
    USE data_grid_lm,        ONLY: istartcos, iendcos, jstartcos, jendcos
    ! MESSy/BMIL
    USE messy_main_mpi_bi,   ONLY: icomm_cart, sendbuf, isendbuflen      &
                                 , imp_reals, nboundlines, my_cart_neigh &
                                 , ncomm_type, ldatatypes, num_compute   &
                                 , exchg_boundaries
    USE messy_main_grid_Def_mem_bi, ONLY: istartpar_c4   => istartpar &
                                  , iendpar_c4     => iendpar   &
                                  , jstartpar_c4   => jstartpar &
                                  , jendpar_c4     => jendpar   &
                                  , ie_c4=> ie,  je_c4 =>  je
    USE messy_main_grid_def_bi,     ONLY:lperi_x, lperi_y, l2dim

    USE messy_main_constants_mem,  ONLY: g

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='move_boundary_arrays_to_Cosmo'
    INTEGER(I4)                 :: ii, iX
    INTEGER                     :: size3, size4
    INTEGER(I4)                 :: kdims(24)
    CHARACTER(LEN=80)           :: yerrmsg
    INTEGER(I4)                 :: stat_I4
    INTEGER(I4)                 :: icount

    icount = 1

    IF (.NOT. ALLOCATED(p_fis)) ALLOCATE (p_fis(ie_c4,je_c4,1,1))
    IF (MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN
       p_fis(istartpar_c4:iendpar_c4,jstartpar_c4:jendpar_c4,1,1) = &
            CplData(idx_cl_fis)%ptr_i2c(istartcos:iendcos,jstartcos:jendcos,1,1)
    ELSE
       ! hsurf multiply by g
       p_fis(istartpar_c4:iendpar_c4,jstartpar_c4:jendpar_c4,1,1) = g * &
            CplData(idx_cl_fis)%ptr_i2c(istartcos:iendcos,jstartcos:jendcos,1,1)
    END IF
    kdims(1:24) = &
         (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                             &
         (100+icount, sendbuf, isendbuflen, imp_reals, icomm_cart     &
         , num_compute, ie_c4, je_c4, kdims, jstartpar_c4, jendpar_c4 &
         , nboundlines, nboundlines, my_cart_neigh                    &
         , lperi_x, lperi_y, l2dim, 50000+icount                      &
         , ldatatypes, ncomm_type, stat_I4, yerrmsg                   &
         , p_fis(:,:,:,1))
    IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)
    icount=icount+1

    ! MOVE  FIELDS
    DO ii = 1, NCOPY
       IF (TRIM(CplData(ii)%CHILD%OBJ) == 'Test_Ar') CYCLE
       l_bound_if: IF (CplData(ii)%L_BOUND) THEN

          lasso: IF (ASSOCIATED(CplData(ii)%cosmo_bd)) THEN
             size3 = SIZE(CplData(ii)%cosmo_bd(1)%ptr,3)
             size4 = SIZE(CplData(ii)%cosmo_bd(1)%ptr,4)
             CplData(ii)%cosmo_bd(1)%ptr(istartpar_c4:iendpar_c4 &
                  ,jstartpar_c4:jendpar_c4,1:size3,1:size4) =    &
                  CplData(ii)%ptr_i2c(istartcos:iendcos          &
                  ,jstartcos:jendcos,1:size3,1:size4)

             IF (size4 == 1 .AND. &
                  CplData(ii)%vartab_idx /= 0) CALL flag_fields (&
                  CplData(ii)%vartab_idx, CplData(ii)%cosmo_bd(1)%ptr)

             kdims(1:24) = &
                  (/size3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
             DO iX=1,size4
                CALL exchg_boundaries                                        &
                     (100+icount, sendbuf, isendbuflen, imp_reals, icomm_cart&
                     , num_compute                                           &
                     , ie_c4, je_c4, kdims, jstartpar_c4, jendpar_c4         &
                     , nboundlines, nboundlines, my_cart_neigh               &
                     , lperi_x, lperi_y, l2dim, 60000+icount                 &
                     , ldatatypes, ncomm_type, stat_I4, yerrmsg              &
                     , CplData(ii)%cosmo_bd(1)%ptr(:,:,:,iX)                 )
                IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)
                icount =icount+1
             ENDDO
             CplData(ii)%cosmo_bd(2)%ptr = CplData(ii)%cosmo_bd(1)%ptr
          ELSE
             CALL error_bi(substr, 'POINTER NOT ASSOCIATED')
          ENDIF lasso

          IF (CplData(ii)%scn /= -99) THEN
             ParDATA(CplData(ii)%scn)%ptr_ori = &
                  CplData(ii)%cosmo_bd(1)%ptr(:,:,:,1:1)
          END IF
          !**************************************************************
       ELSE IF (CplData(ii)%L_INPUT) THEN !l_bound
          !**************************************************************

          lasso_cp: IF (ASSOCIATED(CplData(ii)%cosmo) ) THEN
             size3 = SIZE(CplData(ii)%cosmo(1)%ptr,3)
             size4 = SIZE(CplData(ii)%cosmo(1)%ptr,4)
             CplData(ii)%cosmo(1)%ptr(istartpar_c4:iendpar_c4 &
                  ,jstartpar_c4:jendpar_c4,1:size3,1:size4) = &
                  CplData(ii)%ptr_i2c(istartcos:iendcos       &
                  ,jstartcos:jendcos,1:size3,1:size4)

             IF (size4 == 1 .AND. &
                  CplData(ii)%vartab_idx /= 0) CALL flag_fields (&
                  CplData(ii)%vartab_idx, CplData(ii)%cosmo(1)%ptr)

             kdims(1:24) = &
                  (/size3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
             DO ix=1,size4
                CALL exchg_boundaries                                        &
                     (100+icount, sendbuf, isendbuflen, imp_reals, icomm_cart&
                     , num_compute                                           &
                     , ie_c4, je_c4, kdims, jstartpar_c4, jendpar_c4         &
                     , nboundlines, nboundlines, my_cart_neigh               &
                     , lperi_x, lperi_y, l2dim, 60000+icount                 &
                     , ldatatypes, ncomm_type, stat_I4, yerrmsg              &
                     , CplData(ii)%cosmo(1)%ptr(:,:,:,ix)                    )
                IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)
                icount =icount+1
             ENDDO
          ELSE
             CALL error_bi(substr, 'POINTER NOT ASSOCIATED')
          ENDIF lasso_cp

          ! PARENT coupling l_shortcut =.TRUE.
          IF (CplData(ii)%scn /= -99) THEN
             ParDATA(CplData(ii)%scn)%ptr_ori = CplData(ii)%cosmo(1)%ptr
          END IF

       ENDIF l_bound_if

    ENDDO
    IF (l_shortcut .AND. idx_p_ps /= -99) THEN
       ParDATA(idx_p_ps)%ptr_ori(istartpar_c4:iendpar_c4 &
                  ,jstartpar_c4:jendpar_c4,:,:) =  CplData(idx_cl_ps)%ptr_i2c(&
            istartcos:iendcos,jstartcos:jendcos,:,:)
       kdims(1:24) = &
            (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       CALL exchg_boundaries                                              &
            (100+icount, sendbuf, isendbuflen, imp_reals, icomm_cart      &
            , num_compute, ie_c4, je_c4, kdims, jstartpar_c4, jendpar_c4  &
            , nboundlines, nboundlines, my_cart_neigh                     &
            , lperi_x, lperi_y, l2dim, 60000+icount                       &
            , ldatatypes, ncomm_type, stat_I4, yerrmsg                    &
            , ParDATA(idx_p_ps)%ptr_ori(:,:,:,1)                          )

       IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)
       icount =icount+1
    ENDIF

  END SUBROUTINE move_boundary_arrays_to_Cosmo
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE exchange_breakinfo

    ! MESSy
    USE messy_main_timer,         ONLY: lbreak, l_rerun &
                                      , L_TRIGGER_RESTART
    ! MMD
    USE mmd_child,               ONLY: MMD_Inter_Bcast

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='mmd2way_child_exchange_breakinfo'
    INTEGER,DIMENSION(3)        :: Timeflags

    CALL debug_bi(submodstr, 'EXCHANGE BREAK INFO', '', substr)

    ! EXCHANGE FLAGS FOR TIMER (MODEL RESTART AND/OR STOP)
    CALL MMD_Inter_Bcast (Timeflags, .FALSE.)
    lbreak       = (Timeflags(1) == 1)
    l_rerun      = (Timeflags(2) == 1) .OR. l_rerun

    IF (lbreak)  THEN
       l_rerun           = .TRUE.
       L_TRIGGER_RESTART = .TRUE.
    ENDIF

    CALL debug_bi(submodstr, 'EXCHANGE BREAK INFO', '', substr)

  END SUBROUTINE exchange_breakinfo
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE mmd2way_child_read_nml_cpl(status, iou)

    ! MESSy/SMCL
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check &
                                      , read_nml_close
    USE messy_mmd2way,            ONLY: modstr

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT)         :: status     ! error status
    INTEGER, INTENT(IN)          :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'mmd2way_child_read_nml_cpl'


    NAMELIST /CPL_CHILD/ CPL_IOEVENT, READEXT_IOEVENT&
         , l_forcevars, forcevars, forcefile &
         , l_I2Cori_output, WRITEI2C_IOEVENT

    ! LOCAL
    LOGICAL :: lex      ! file exists
    INTEGER :: fstat

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL_CHILD', modstr)
    IF (.not.lex) RETURN    ! <submodstr>.nml does not exist

    READ(iou, NML=CPL_CHILD, IOSTAT=fstat)

    CALL read_nml_check(fstat, substr, iou, 'CPL_CHILD', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, submodstr)

    status = 0  ! no ERROR

  END SUBROUTINE mmd2way_child_read_nml_cpl
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE mmd2way_child_read_nml_cpl_serv(status, iou, PARENT_TAG)

    ! MESSy/SMCL
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check &
                                      , read_nml_close
    USE messy_mmd2way,            ONLY: modstr

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(OUT)         :: status     ! error status
    INTEGER, INTENT(IN)          :: iou        ! I/O unit
    CHARACTER(LEN=1), INTENT(IN) :: PARENT_TAG

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'mmd2way_child_read_nml_cpl_serv'

    NAMELIST /CPL_CHILD_ECHAM/ FIELD
    NAMELIST /CPL_CHILD_COSMO/ FIELD

    ! LOCAL
    LOGICAL           :: lex      ! file exists ?
    INTEGER           :: fstat    ! file status
    INTEGER           :: ii
    CHARACTER(LEN=15) :: NML

    status = 1
    ! SET DEFAULTS
    IF (PARENT_TAG =='E') THEN
       NML = 'CPL_CHILD_ECHAM'
    ELSE
       NML = 'CPL_CHILD_COSMO'
    ENDIF

    CALL read_nml_open(lex, substr, iou, NML, modstr)
    IF (.not.lex) RETURN    ! <submodstr>.nml does not exist

    IF (PARENT_TAG =='E') THEN
       READ(iou, NML=CPL_CHILD_ECHAM, IOSTAT=fstat)
    ELSE
       READ(iou, NML=CPL_CHILD_COSMO, IOSTAT=fstat)
    ENDIF
    CALL read_nml_check(fstat, substr, iou, NML, modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES

    cplfield_loop: DO ii = 1, NMAX_EXCH

       IF (TRIM(FIELD(ii)%PARENT%CHA) == '') CYCLE
       IF (TRIM(FIELD(ii)%PARENT%OBJ) == '') CYCLE
       IF (TRIM(FIELD(ii)%CHILD%CHA) == '') CYCLE
       IF (TRIM(FIELD(ii)%CHILD%OBJ) == '') CYCLE

       WRITE(*,*) '.........................................................'
       WRITE(*,*) ' NUMBER OF CPL-FIELDS: ',ii
       WRITE(*,*) ' CHANNEL PARENT         = ', TRIM(FIELD(ii)%PARENT%CHA)
       WRITE(*,*) ' CHANNEL OBJECT PARENT  = ', TRIM(FIELD(ii)%PARENT%OBJ)

       WRITE(*,*) ' CHANNEL CHILD         = ', TRIM(FIELD(ii)%CHILD%CHA)
       WRITE(*,*) ' CHANNEL OBJECT CHILD  = ', TRIM(FIELD(ii)%CHILD%OBJ)

       WRITE(*,*) ' INTERPOLATION METHOD   = ',  FIELD(ii)%C_INTERPOL
       WRITE(*,*) ' CALCULATE BOUNDARIES   = ',  FIELD(ii)%L_BOUND
       WRITE(*,*) ' CALCULATE INITIAL FIELD= ',  FIELD(ii)%L_INITIAL
       WRITE(*,*) ' INPUT FIELD            = ',  FIELD(ii)%L_INPUT
       WRITE(*,*) ' REPRESENTATION         = ',  FIELD(ii)%C_REPR
       WRITE(*,*) '.........................................................'
    END DO cplfield_loop

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE mmd2way_child_read_nml_cpl_serv
  ! ----------------------------------------------------------------------

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! PARENT COUPLING +
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ----------------------------------------------------------------------
  SUBROUTINE parent_assign_ParData

    ! MMD
    USE mmd_child,               ONLY: MMD_C_GetNextParArray

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    INTEGER :: ii
    CHARACTER(LEN=STRLEN_CHANNEL)   :: myChannel
    CHARACTER(LEN=STRLEN_OBJECT)    :: myName
    CHARACTER(LEN=STRLEN_MEDIUM)    :: repr
    INTEGER                         :: interpM
    LOGICAL                         :: l_SendUnit
    INTEGER                         :: ix

    ! NO PARENT FIELDS
    IF (PD_NUM == 0) RETURN

    ALLOCATE(ParData(PD_NUM+PD_ADD_MAX))

    ii = 0
    DO WHILE (MMD_C_GetNextParArray (myChannel, myName, repr, interpM &
         , l_SendUnit))

       ii = ii + 1
       !write(iouerr,*) 'parent_assign_ParData ' &
       !     ,TRIM(myChannel), ' ',TRIM(myName), interpM, TRIM(repr)
       ParData(ii)%name%CHA = TRIM(myChannel)
       ParData(ii)%name%OBJ = TRIM(myName)
       ParData(ii)%repr     = TRIM(repr)
       ParData(ii)%interpM  = interpM
       ParData(ii)%l_SendUnit = l_SendUnit

       IF (TRIM(myChannel) == 'mmd2way_child' .OR. l_shortcut) THEN
          DO ix = 1, NEXCH
             IF (TRIM(CplData(ix)%CHILD%OBJ) == TRIM(MyName)) THEN
                CplData(ix)%scn = ii
                EXIT
             END IF
          END DO
       END IF
    END DO

  END SUBROUTINE parent_assign_ParData
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  SUBROUTINE parent_make_representations

    ! MESSy/BMIL
    USE messy_main_grid_def_mem_bi,      ONLY: ke
    USE messy_main_channel_bi,           ONLY: DC_MMD_OUT &
                                             , DC_MMD_OUTACC               &
                                             , gp_meml, gp_nseg, gp_start
    USE messy_main_channel_error_bi,     ONLY: channel_halt

    ! MESSy/SMCL
    USE messy_main_channel_dimensions,   ONLY: new_dimension             &
                                             , add_dimension_variable    &
                                             , add_dimension_variable_att
    USE messy_main_channel_repr,         ONLY: new_representation        &
                                             , set_representation_decomp &
                                             , AUTO, IRANK, PIOTYPE_COL
    ! INT2COSMO
    USE data_grid_in,                    ONLY: ke_in, ie_in, je_in
    USE data_int2lm_parallel,            ONLY: ie_in_tot_red, je_in_tot_red



    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'parent_make_representations'
    INTEGER                     :: status
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    INTEGER :: i
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: array


    ! NO PARENT FIELDS
    IF (PD_NUM == 0) RETURN

    ! 1.) make dimensions
    CALL new_dimension(status, DIMID_POUT_IE, 'MMDOUT_ie', ie_in_tot_red)
    CALL channel_halt(substr, status)

    CALL new_dimension(status, DIMID_POUT_JE, 'MMDOUT_je', je_in_tot_red)
    CALL channel_halt(substr, status)

    CALL new_dimension(status, DIMID_POUT_KMIN, 'MMDOUT_kmin', ke_int)
    CALL channel_halt(substr, status)


    ALLOCATE(array(ke_int))
    DO i=par_kmin, ke_in
       array(i-par_kmin+1) = REAL(i, DP)
    END DO
    CALL add_dimension_variable(status, 'MMDOUT_kmin', 'MMDOUT_kmin', array)
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'MMDOUT_kmin', 'MMDOUT_kmin', &
         'positive', c='down')
    CALL channel_halt(substr, status)
    DEALLOCATE(array)

    CALL new_dimension(status, DIMID_POUT_KMINp1, 'MMDOUT_kminp1', ke_int+1)
    CALL channel_halt(substr, status)
    ALLOCATE(array(ke_int+1))
    DO i=par_kmin, ke_in+1
       array(i-par_kmin+1) = REAL(i, DP)
    END DO
    CALL add_dimension_variable(status, 'MMDOUT_kminp1', 'MMDOUT_kminp1', array)
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'MMDOUT_kminp1', 'MMDOUT_kminp1', &
         'positive', c='down')
    CALL channel_halt(substr, status)
    DEALLOCATE(array)

    ke_hint = ke - chi_kmin + 1
    CALL new_dimension(status, DIMID_POUT_HINT, 'MMDOUT_hint', ke_hint)
    CALL channel_halt(substr, status)

    ALLOCATE(array(ke_hint))
    DO i=chi_kmin, ke
       array(i-chi_kmin+1) = REAL(i, DP)
    END DO
    CALL add_dimension_variable(status,'MMDOUT_hint' , 'MMDOUT_hint', array)
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'MMDOUT_hint', 'MMDOUT_hint', &
         'positive', c='down')
    CALL channel_halt(substr, status)
    DEALLOCATE(array)

    CALL new_dimension(status, DIMID_POUT_HINTp1, 'MMDOUT_hintp1', ke_hint+1)
    CALL channel_halt(substr, status)

    ALLOCATE(array(ke_hint+1))
    DO i=chi_kmin, ke+1
       array(i-chi_kmin+1) = REAL(i, DP)
    END DO
    CALL add_dimension_variable(status,'MMDOUT_hintp1' , 'MMDOUT_hintp1', array)
    CALL channel_halt(substr, status)
    CALL add_dimension_variable_att(status, 'MMDOUT_hintp1', 'MMDOUT_hintp1', &
         'positive', c='down')
    CALL channel_halt(substr, status)
    DEALLOCATE(array)

    ! 2.) make representations
    CALL new_representation(status, REPR_POUT_2D, 'REPR_POUT_2D' &
       , rank = 2, link = 'xx--', dctype = DC_MMD_OUT            &
       , dimension_ids = (/ DIMID_POUT_IE, DIMID_POUT_JE  /)     &
       , ldimlen       = (/ ie_in        , je_in          /)     &
       , output_order  = (/ 1,2 /)                               &
       , axis = 'XY--'                                           &
       )
    CALL channel_halt(substr, status)

    nseg = gp_nseg
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:)   = 1
    meml(:,:)  = 1
    memu(:,:)  = 1

    start(:,:) = gp_start(:,:)
    cnt(:,1)   = ie_in
    cnt(:,2)   = je_in
    cnt(:,3)   = 1
    meml(:,:)  = gp_meml(:,:)
    memu(:,1)  = ie_in
    memu(:,2)  = je_in
    memu(:,3)  = 1

    CALL set_representation_decomp(status, REPR_POUT_2D &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------

    CALL new_representation(status, REPR_POUT_3D_MID, 'REPR_POUT_3D_MID'      &
       , rank = 3, link = 'xxx-', dctype = DC_MMD_OUT                         &
       , dimension_ids = (/ DIMID_POUT_IE, DIMID_POUT_JE, DIMID_POUT_KMIN /)  &
       , ldimlen       = (/ ie_in        , je_in        ,  AUTO           /)  &
       , output_order  = (/ 1,2,3 /)                                          &
       , axis = 'XYZ-'                                                        &
       )
    CALL channel_halt(substr, status)
    nseg = gp_nseg

    start(:,:) = gp_start(:,:)
    cnt(:,_IX_XYZ__)   = ie_in
    cnt(:,_IY_XYZ__)   = je_in
    cnt(:,_IZ_XYZ__)   = ke_int

    meml(:,:)  = gp_meml(:,:)
    memu(:,_IX_XYZ__)  = ie_in
    memu(:,_IY_XYZ__)  = je_in
    memu(:,_IZ_XYZ__)  = ke_int

    CALL set_representation_decomp(status, REPR_POUT_3D_MID &
         , start, cnt, memu, meml, .TRUE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    !--------------------------------------------------------------------

    !--------------------------------------------------------------------

    CALL new_representation(status, REPR_POUT_3D_INT, 'REPR_POUT_3D_INT'      &
       , rank = 3, link = 'xxx-', dctype = DC_MMD_OUT                         &
       , dimension_ids = (/ DIMID_POUT_IE, DIMID_POUT_JE, DIMID_POUT_KMINp1/) &
       , ldimlen       = (/ ie_in        , je_in        ,  AUTO            /) &
       , output_order  = (/ 1,2,3 /)                                          &
       , axis = 'XYZ-'                                                        &
       )
    CALL channel_halt(substr, status)
    nseg = gp_nseg

    start(:,:) = gp_start(:,:)
    cnt(:,_IX_XYZ__)   = ie_in
    cnt(:,_IY_XYZ__)   = je_in
    cnt(:,_IZ_XYZ__)   = ke_int + 1
    meml(:,:)  = gp_meml(:,:)
    memu(:,_IX_XYZ__)  = ie_in
    memu(:,_IY_XYZ__)  = je_in
    memu(:,_IZ_XYZ__)  = ke_int + 1

    CALL set_representation_decomp(status, REPR_POUT_3D_INT &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    !--------------------------------------------------------------------
    CALL new_representation(status, REPR_PHOUT_3D_MID, 'REPR_PHOUT_3D_MID'    &
       , rank = 3, link = 'xxx-', dctype = DC_MMD_OUT                         &
       , dimension_ids = (/ DIMID_POUT_IE, DIMID_POUT_JE, DIMID_POUT_HINT /)  &
       , ldimlen       = (/ ie_in        , je_in        ,  AUTO           /)  &
       , output_order  = (/ 1,2,3 /)                                          &
       , axis = 'XYZ-'                                                        &
       )
    CALL channel_halt(substr, status)
    nseg = gp_nseg

    start(:,:) = gp_start(:,:)
    cnt(:,_IX_XYZ__)   = ie_in
    cnt(:,_IY_XYZ__)   = je_in
    cnt(:,_IZ_XYZ__)   = ke_hint

    meml(:,:)  = gp_meml(:,:)
    memu(:,_IX_XYZ__)  = ie_in
    memu(:,_IY_XYZ__)  = je_in
    memu(:,_IZ_XYZ__)  = ke_hint

    CALL set_representation_decomp(status, REPR_PHOUT_3D_MID &
         , start, cnt, memu, meml, .TRUE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    !--------------------------------------------------------------------

    !--------------------------------------------------------------------

    CALL new_representation(status, REPR_PHOUT_3D_INT, 'REPR_PHOUT_3D_INT'    &
       , rank = 3, link = 'xxx-', dctype = DC_MMD_OUT                         &
       , dimension_ids = (/ DIMID_POUT_IE, DIMID_POUT_JE, DIMID_POUT_HINTp1/) &
       , ldimlen       = (/ ie_in        , je_in        ,  AUTO            /) &
       , output_order  = (/ 1,2,3 /)                                          &
       , axis = 'XYZ-'                                                        &
       )
    CALL channel_halt(substr, status)
    nseg = gp_nseg

    start(:,:) = gp_start(:,:)
    cnt(:,_IX_XYZ__)   = ie_in
    cnt(:,_IY_XYZ__)   = je_in
    cnt(:,_IZ_XYZ__)   = ke_hint + 1
    meml(:,:)  = gp_meml(:,:)
    memu(:,_IX_XYZ__)  = ie_in
    memu(:,_IY_XYZ__)  = je_in
    memu(:,_IZ_XYZ__)  = ke_hint + 1

    CALL set_representation_decomp(status, REPR_PHOUT_3D_INT &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    DEALLOCATE(start)
    DEALLOCATE(cnt)
    DEALLOCATE(meml)
    DEALLOCATE(memu)
    !--------------------------------------------------------------------

    !--------------------------------------------------------------------
    CALL new_representation(status, REPR_POUTACC_2D, 'REPR_POUTACC_2D' &
       , rank = 2, link = 'xx--', dctype = DC_MMD_OUTACC               &
       , dimension_ids = (/ DIMID_POUT_IE, DIMID_POUT_JE  /)           &
       , ldimlen       = (/ ie_in        , je_in          /)           &
       , output_order  = (/ 1,2 /)                                     &
       , axis = 'XY--'                                                 &
       )
    CALL channel_halt(substr, status)

    nseg = gp_nseg
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:)   = 1
    meml(:,:)  = 1
    memu(:,:)  = 1

    start(:,:) = gp_start(:,:)
    cnt(:,_IX_XYZ__)   = ie_in
    cnt(:,_IY_XYZ__)   = je_in

    meml(:,:)  = gp_meml(:,:)
    memu(:,_IX_XYZ__)  = ie_in
    memu(:,_IY_XYZ__)  = je_in

    CALL set_representation_decomp(status, REPR_POUTACC_2D &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------

    CALL new_representation(status, REPR_POUTACC_3D_MID, 'REPR_POUTACC_3D_MID'&
       , rank = 3, link = 'xxx-', dctype = DC_MMD_OUTACC                      &
       , dimension_ids = (/ DIMID_POUT_IE, DIMID_POUT_JE, DIMID_POUT_KMIN /)  &
       , ldimlen       = (/ ie_in        , je_in        ,  AUTO           /)  &
       , output_order  = (/ 1,2,3 /)                                          &
       , axis = 'XYZ-'                                                        &
       )
    CALL channel_halt(substr, status)
    nseg = gp_nseg

    start(:,:) = gp_start(:,:)
    cnt(:,_IX_XYZ__)   = ie_in
    cnt(:,_IY_XYZ__)   = je_in
    cnt(:,_IZ_XYZ__)   = ke_int

    meml(:,:)  = gp_meml(:,:)
    memu(:,_IX_XYZ__)  = ie_in
    memu(:,_IY_XYZ__)  = je_in
    memu(:,_IZ_XYZ__)  = ke_int

    CALL set_representation_decomp(status, REPR_POUTACC_3D_MID &
         , start, cnt, memu, meml, .TRUE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    !--------------------------------------------------------------------
    CALL new_representation(status, REPR_PHOUTACC_3D_MID, 'REPR_PHOUTACC_3D_MID'&
       , rank = 3, link = 'xxx-', dctype = DC_MMD_OUTACC                        &
       , dimension_ids = (/ DIMID_POUT_IE, DIMID_POUT_JE, DIMID_POUT_HINT /)    &
       , ldimlen       = (/ ie_in        , je_in        ,  AUTO           /)    &
       , output_order  = (/ 1,2,3 /)                                            &
       , axis = 'XYZ-'                                                          &
       )
    CALL channel_halt(substr, status)
    nseg = gp_nseg

    start(:,:) = gp_start(:,:)
    cnt(:,_IX_XYZ__)   = ie_in
    cnt(:,_IY_XYZ__)   = je_in
    cnt(:,_IZ_XYZ__)   = ke_hint
    meml(:,:)  = gp_meml(:,:)
    memu(:,_IX_XYZ__)  = ie_in
    memu(:,_IY_XYZ__)  = je_in
    memu(:,_IZ_XYZ__)  = ke_hint

    CALL set_representation_decomp(status, REPR_PHOUTACC_3D_MID &
         , start, cnt, memu, meml, .TRUE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------


    !--------------------------------------------------------------------
    CALL new_representation(status, REPR_POUTACC_3D_INT, 'REPR_POUTACC_3D_INT'&
       , rank = 3, link = 'xxx-', dctype = DC_MMD_OUTACC                      &
       , dimension_ids = (/ DIMID_POUT_IE, DIMID_POUT_JE, DIMID_POUT_KMINp1/) &
       , ldimlen       = (/ ie_in        , je_in        ,  AUTO            /) &
       , output_order  = (/ 1,2,3 /)                                          &
       , axis = 'XYZ-'                                                        &
       )
    CALL channel_halt(substr, status)
    nseg = gp_nseg

    start(:,:) = gp_start(:,:)
    cnt(:,_IX_XYZ__)   = ie_in
    cnt(:,_IY_XYZ__)   = je_in
    cnt(:,_IZ_XYZ__)   = ke_int+1

    meml(:,:)  = gp_meml(:,:)
    memu(:,_IX_XYZ__)  = ie_in
    memu(:,_IY_XYZ__)  = je_in
    memu(:,_IZ_XYZ__)  = ke_int+1

    CALL set_representation_decomp(status, REPR_POUTACC_3D_INT &
         , start, cnt, memu, meml, .TRUE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    !--------------------------------------------------------------------
    CALL new_representation(status, REPR_PHOUTACC_3D_INT, 'REPR_PHOUTACC_3D_INT'&
       , rank = 3, link = 'xxx-', dctype = DC_MMD_OUTACC                        &
       , dimension_ids = (/ DIMID_POUT_IE, DIMID_POUT_JE, DIMID_POUT_HINTp1/)   &
       , ldimlen       = (/ ie_in        , je_in        ,  AUTO            /)   &
       , output_order  = (/ 1,2,3 /)                                            &
       , axis = 'XYZ-'                                                          &
       )
    CALL channel_halt(substr, status)
    nseg = gp_nseg

    start(:,:) = gp_start(:,:)
    cnt(:,_IX_XYZ__)   = ie_in
    cnt(:,_IY_XYZ__)   = je_in
    cnt(:,_IZ_XYZ__)   = ke_hint+1
    meml(:,:)  = gp_meml(:,:)
    memu(:,_IX_XYZ__)  = ie_in
    memu(:,_IY_XYZ__)  = je_in
    memu(:,_IZ_XYZ__)  = ke_hint+1

    CALL set_representation_decomp(status, REPR_PHOUTACC_3D_INT &
         , start, cnt, memu, meml, .TRUE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    !--------------------------------------------------------------------
    DEALLOCATE(start)
    DEALLOCATE(cnt)
    DEALLOCATE(meml)
    DEALLOCATE(memu)

  END SUBROUTINE parent_make_representations
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  SUBROUTINE  parent_set_ParData

    ! MESSy/BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, GP_2D_HORIZONTAL
    USE messy_main_grid_def_mem_bi,  ONLY: ke
    ! MESSy/SMCL
    USE messy_main_channel,       ONLY: new_channel             &
                                      , new_channel_object      &
                                      , get_channel_object
    USE messy_main_channel_repr,  ONLY: get_representation_info &
                                      , get_representation_id
    USE messy_main_tools,         ONLY: str, int2str

    ! MMD
    USE mmd_child,                ONLY: MMD_C_GetNextParArray  &
                                      , MMD_C_Set_ParDataArray


    ! INT2COSMO
    USE data_int2lm_constants,    ONLY: r_earth, pi ! for consistency
    USE data_grid_in,             ONLY: je_in, dlon_in, dlat_in

    IMPLICIT NONE

    INTRINSIC :: TRIM

    CHARACTER(LEN=*), PARAMETER   :: substr = 'parent_set_ParData'
    CHARACTER(LEN=*), PARAMETER   :: chname = 'MMDC4_OUT'
    CHARACTER(LEN=STRLEN_CHANNEL) :: myChannel
    CHARACTER(LEN=STRLEN_OBJECT)  :: myName
    CHARACTER(LEN=STRLEN_MEDIUM)  :: repr
    INTEGER                       :: interpM
    LOGICAL                       :: lsu    ! dummy use
    INTEGER                       :: ii     ! loop index
    INTEGER                       :: status ! error status
    INTEGER                       :: reprId = 0
    INTEGER                       :: iter
    CHARACTER(LEN=2)              :: istr
    INTEGER                       :: jx

    ! NO PARENT FIELDS
    IF (PD_NUM == 0) RETURN

    ! DEFINE SUBMODEL-CHANNEL
    CALL new_channel(status, chname)
    CALL channel_halt(substr, status)

    ! define channel object for weight function
    CALL new_channel_object(status, chname, 'fw_int' &
                           , p4=cofw_int, reprid=REPR_POUTACC_2D)
    CALL channel_halt(substr, status)
    cofw_int = fw_int

    DEALLOCATE(fw_int)
    ! define channel object for weight fractions
    CALL new_channel_object(status, chname, 'fractions' &
         , p4=fractions, reprid=REPR_POUTACC_2D)
    CALL channel_halt(substr, status)
    fractions(:,:,1,1) = fracs
    DEALLOCATE(fracs)

    ! debug+
    CALL get_representation_id(status, 'REPR_TIMER' , reprid)
    CALL channel_halt(substr, status)

    IF (MMD_C_GetParentType() == MMD_ParentIsEcham) THEN

       CALL new_channel_object(status, chname, 'NPSITERMAX'  &
            , p1=npsitermax, reprid=reprid)
       CALL channel_halt(substr, status)
       CALL new_channel_object(status, chname, 'NPSITER1'    &
            , p1=npsiter1, reprid=reprid)
       CALL channel_halt(substr, status)
       CALL new_channel_object(status, chname, 'NPSITER01'   &
            , p1=npsiter01, reprid=reprid)
       CALL channel_halt(substr, status)
       CALL new_channel_object(status, chname, 'NPSITER001'  &
            , p1=npsiter001, reprid=reprid)
       CALL channel_halt(substr, status)
       CALL new_channel_object(status, chname, 'NPSITER0001' &
            , p1=npsiter0001, reprid=reprid)
       CALL channel_halt(substr, status)

       ! Calculate boxarea for diagnosis in interpol_parent_data
       CALL new_channel_object(status,chname,'boxarea_out' &
            , p2=boxarea_out, reprid=REPR_POUTACC_2D)
       CALL channel_halt(substr, status)

       DO jx = 1, je_in
          boxarea_out(:,jx) = COS (rlats(jx) * pi /180._dp) &
            * r_earth**2 * (pi/180.0_dp)**2 * dlon_in * dlat_in
       END DO

       ALLOCATE(iterps(nitermax))
       DO iter = 1,nitermax
          CALL int2str(istr, iter)
          CALL new_channel_object(status, chname, 'iteraps'//istr &
               , p2=iterps(iter)%ptr &
               , reprid = REPR_POUTACC_2D, lrestreq=.FALSE.)
          CALL channel_halt(substr, status)
       END DO
    ENDIF
    ! debug-

    ii = 0
    DO WHILE (MMD_C_GetNextParArray (myChannel, myName, repr, interpM &
         , lsu))

       ii = ii + 1

          SELECT CASE (TRIM(myName))
          CASE('QV')
             idx_p_qv   = ii
          CASE('QC')
             idx_p_qc   = ii
          CASE('QI')
             idx_p_qi   = ii
          CASE('U')
             idx_p_u    = ii
          CASE ('V')
             idx_p_v    = ii
          CASE ('T')
             idx_p_t    = ii
          CASE ('PS')
             idx_p_ps   = ii
          CASE ('HSURF')
             idx_p_hsurf  = ii
          CASE ('W')
             idx_p_w    = ii
          CASE ('PP')
             idx_p_pp   = ii
          CASE ('T_S')
             idx_p_ts   = ii
          CASE DEFAULT
             ! NOTHING TO DO
          END SELECT


          SELECT CASE (TRIM(ParData(ii)%repr))
          CASE("GP_2D_HORIZONTAL")
             ! TARGET VARIABLE is STANDARD 2D
             CALL new_channel_object(status,chname,TRIM(ParData(ii)%name%obj) &
                  , p4=ParData(ii)%ptr_int, reprid=REPR_POUTACC_2D)
             IF (status == 3102) THEN ! channel object exists already
                CALL new_channel_object(status, chname               &
                     , TRIM(ParData(ii)%name%obj)//'_'//str(ii)      &
                     , p4=ParData(ii)%ptr_int, reprid=REPR_POUTACC_2D)
             ENDIF
             CALL channel_halt(substr, status)
             reprId = REPR_POUTACC_2D

             IF (ii == idx_p_ps) THEN
               CALL new_channel_object(status, chname                &
                    , 'H_'//TRIM(ParData(ii)%name%obj)               &
                    , p4=ParData(ii)%ptr_hint, reprid=REPR_POUTACC_2D)
               IF (status == 3102) THEN ! channel object exists already
                  CALL new_channel_object(status, chname                &
                       , 'H_'//TRIM(ParData(ii)%name%obj)//'_'//str(ii) &
                       , p4=ParData(ii)%ptr_hint, reprid=REPR_POUTACC_2D)
               ENDIF
             ELSE
                ParData(ii)%ptr_hint => ParData(ii)%ptr_int
             END IF

          CASE ("GP_3D_MID")
             ! TARGET VARIABLE is STANDARD 3D
             CALL new_channel_object(status,chname,TRIM(ParData(ii)%name%obj) &
                  , p4=ParData(ii)%ptr_int, reprid=REPR_POUTACC_3D_MID)
             IF (status == 3102) THEN ! channel object exists already
                CALL new_channel_object(status, chname                   &
                     , TRIM(ParData(ii)%name%obj)//'_'//str(ii)          &
                     , p4=ParData(ii)%ptr_int, reprid=REPR_POUTACC_3D_MID)
             ENDIF
             CALL channel_halt(substr, status)
             reprId = REPR_POUTACC_3D_MID

             CALL new_channel_object(status, chname                     &
                  ,'H_'//TRIM(ParData(ii)%name%obj)                     &
                  , p4=ParData(ii)%ptr_hint, reprid=REPR_PHOUTACC_3D_MID)
             IF (status == 3102) THEN ! channel object exists already
                CALL new_channel_object(status ,chname                     &
                     , TRIM(ParData(ii)%name%obj)//'_'//str(ii)            &
                     , p4=ParData(ii)%ptr_hint, reprid=REPR_PHOUTACC_3D_MID)
             ENDIF
             CALL channel_halt(substr, status)

             ! set 3d index
             idx_p_3d = ii

          CASE ("GP_3D_INT")
             ! TARGET VARIABLE is STANDARD 3D
             CALL new_channel_object(status,chname,TRIM(ParData(ii)%name%obj) &
                  , p4=ParData(ii)%ptr_int, reprid=REPR_POUTACC_3D_INT)
             IF (status == 3102) THEN ! channel object exists already
                CALL new_channel_object(status, chname                   &
                     , TRIM(ParData(ii)%name%obj)//'_'//str(ii)          &
                     , p4=ParData(ii)%ptr_int, reprid=REPR_POUTACC_3D_INT)
             ENDIF
             CALL channel_halt(substr, status)
             reprId = REPR_POUTACC_3D_INT

             CALL new_channel_object(status, chname                     &
                  ,'H_'//TRIM(ParData(ii)%name%obj)                     &
                  , p4=ParData(ii)%ptr_hint, reprid=REPR_PHOUTACC_3D_INT)
             IF (status == 3102) THEN ! channel object exists already
                CALL new_channel_object(status, chname                     &
                     , TRIM(ParData(ii)%name%obj)//'_'//str(ii)            &
                     , p4=ParData(ii)%ptr_hint, reprid=REPR_PHOUTACC_3D_INT)
             ENDIF
             CALL channel_halt(substr, status)
             ! TARGET VARIABLE is STANDARD 3D

          CASE DEFAULT
             CALL error_bi(substr,                           &
                  'REPRESENTATION '//TRIM(ParData(ii)%repr)//&
                  &'NOT YET IMPLEMENTED FOR PARENT COUPLING')
          END SELECT
          ! DEFINE AXIS STRING
          CALL get_representation_info(status, ' ', reprId &
               , rank=ParDATA(ii)%rank                     &
               , axis=ParDATA(ii)%AXIS                     &
               , ldimlen=ParDATA(ii)%ldimlen)
          CALL channel_halt(substr, status)

          IF (idx_int1 == -99) THEN
             IF (ParDATA(ii)%interpM == 1) idx_int1 = ii
          END IF

       CALL MMD_C_Set_ParDataArray(status, ParData(ii)%ldimlen &
            , ParData(ii)%axis, p4=ParData(ii)%ptr_int)
       IF (status /=0)  &
            CALL error_bi('ERROR in MMD_C_Set_ParDataArray',substr)
    END DO

    IF (idx_int1 == -99) CALL error_bi(&
         'we need at least ONE element with interpM = 1', substr)

    IF (ldiagonly) THEN
       IF (idx_p_ps == -99) &
            CALL add_parent_list_element('PS', 1, 2, idx_p_ps, 'COSMO_ORI' &
                                              , lhor=.TRUE., itimedep = 2)
        RETURN
    END IF

    IF (MMD_C_GetParentType() == MMD_ParentIsEcham) THEN

       IF (idx_p_ps == -99) THEN
          IF (l_shortcut) THEN
             CALL add_parent_list_element('PS', 1, 2, idx_p_ps, lhor=.TRUE.)
          ELSE
             CALL add_parent_list_element('PS', 1, 2, idx_p_ps, 'COSMO_ORI' &
                                              , lhor=.TRUE., itimedep = 2)
          ENDIF
       END IF
       CALL add_parent_list_element('FIS', 1, 2, idx_p_fis)
       ! FIC is required for hydrostatic pressure adaption
       ! FOR HYDROSTATIC PRESSURE ADAPTION fis NEEDS TO BE CREATED AND
       ! INTERPOLATED
       CALL add_parent_list_element('FIC', 1, 2, idx_p_fic)

       ! ADD 3D PRESSURE
       CALL add_parent_list_element('PRES', 1, 3, idx_p_pres)

       ! IF  missing add temperature (required for lgrhum and ps adaption)
       IF (idx_p_t == -99) THEN
          CALL add_parent_list_element('T', 1, 3, idx_p_t, 'COSMO_ORI' &
                                          , itimedep = 3)
          IF (idx_p_3d == -99) idx_p_3d = idx_p_t
       ENDIF

       ! IF  missing add qv and qc (required for lgrhum and ps adaption)
       IF (idx_p_qv == -99) &
            CALL add_parent_list_element('QV', 1, 3, idx_p_qv, 'tracer_gp')

       IF (idx_p_qc == -99) &
         CALL add_parent_list_element('QC', 1, 3, idx_p_qc, 'tracer_gp')



    ELSE ! COSMO

       ! IF  missing add temperature (required for lgrhum)
       IF (lgrhum .AND. idx_p_t == -99) &
            CALL add_parent_list_element('T', 1, 3, idx_p_t, 'COSMO_ORI'   &
                                            , itimedep = 3)

       IF (lgrhum .AND. idx_p_ps == -99) &
            CALL add_parent_list_element('PS', 1, 2, idx_p_ps, 'COSMO_ORI' &
                                            , itimedep = 2)

       ! surface temperature required for vertical interpolation (vert_interp)
       IF (idx_p_ts == -99) &
            CALL add_parent_list_element('T_S', 1, 2, idx_p_ts, 'COSMO_ORI'&
                                               , itimedep = 2)

       ! vertical height information required for vertical interp.(vert_interp)
       IF (idx_p_hsurf == -99) &
           CALL add_parent_list_element('HSURF', 1, 2, idx_p_hsurf, 'COSMO_ORI')

       ! pressure deviation required for vertical interpolation
       IF (idx_p_pp == -99) THEN
            CALL add_parent_list_element('PP', 1, 3, idx_p_pp, 'COSMO_ORI' &
                                             , itimedep = 3)
            IF (idx_p_3d == -99) idx_p_3d = idx_p_pp
       ENDIF

    ENDIF

    RETURN

  CONTAINS

    SUBROUTINE add_parent_list_element(name, ipol, rank, idx, chaname, lhor &
                                       , itimedep)

!      IMPLICIT NONE

      CHARACTER(LEN=*),           INTENT(IN)  :: name
      INTEGER,                    INTENT(IN)  :: ipol ! interpolation method
      INTEGER,                    INTENT(IN)  :: rank ! rank of object
      INTEGER,                    INTENT(OUT) :: idx
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: chaname
      LOGICAL,          OPTIONAL, INTENT(IN)  :: lhor
      INTEGER,          OPTIONAL, INTENT(IN)  :: itimedep

      ! LOCAL
      INTEGER :: reprid, hreprid
      LOGICAL :: lh

      ! pressure is required for hydrostatic pressure adaption and calculation
      ! of generalised humidity
      PD_ADD = PD_ADD+1
      idx    = PD_NUM+PD_ADD

      ParData(idx)%name%CHA = '#XXX'
      ParData(idx)%name%OBJ = TRIM(name)
      ParData(idx)%interpM  = ipol
      IF (PRESENT(itimedep)) THEN
         ParData(ii)%itimedep  = itimedep
      ELSE
         ParData(ii)%itimedep  = 0
      ENDIF

      IF (PRESENT(chaname)) THEN
         CALL get_channel_object(status,TRIM(chaname), TRIM(name) &
              ,p4=ParData(idx)%ptr_ori)
         CALL channel_halt(substr, status)
      ELSE
        IF (l_shortcut) THEN
            IF (rank == 3) THEN
               CALL new_channel_object(status,chname &
                    , TRIM(ParData(idx)%name%obj)&
                    , p4=ParData(idx)%ptr_ori, reprid=GP_3D_MID)
            ELSE
               CALL new_channel_object(status,chname &
                    , TRIM(ParData(idx)%name%obj)&
                    , p4=ParData(idx)%ptr_ori, reprid=GP_2D_HORIZONTAL)
            ENDIF
            CALL channel_halt(substr, status)

         ELSE
            ! set dummy value
            IF (rank > 2) THEN
               ALLOCATE(ParData(idx)%ptr_ori(1,1,ke,1))
            ELSE
               ALLOCATE(ParData(idx)%ptr_ori(1,1,1,1))
            ENDIF
         ENDIF
      ENDIF

      IF (rank == 3 .AND. idx_p_w /= ii) THEN
         reprid  = REPR_POUTACC_3D_MID
         hreprid = REPR_PHOUTACC_3D_MID
         lh      = .TRUE.
      ELSE IF (rank == 3 ) THEN
         reprid  = REPR_POUTACC_3D_INT
         hreprid = REPR_PHOUTACC_3D_INT
         lh      = .TRUE.

      ELSE
         reprid  = REPR_POUTACC_2D
         hreprid = REPR_POUTACC_2D
         IF (PRESENT(lhor)) THEN
            lh = lhor
         ELSE
            lh = .FALSE.
         END IF
      END IF

      IF (ldiagonly) THEN
         IF (TRIM(ADJUSTL(name)) == 'PS') THEN
            CALL get_channel_object(status, 'MMDC4_IN', 'PS_IN' &
                 , p4=ParData(idx)%ptr_int)
            CALL channel_halt(substr, status)
         ENDIF
      ELSE
         CALL new_channel_object(status, chname, TRIM(ParData(idx)%name%obj)&
              , p4=ParData(idx)%ptr_int, reprid=reprid)
         IF (status == 3102) THEN ! channel object exists already
            CALL new_channel_object(status, chname            &
                 , TRIM(ParData(idx)%name%obj)//'_'//str(idx) &
                 , p4=ParData(idx)%ptr_int, reprid=reprid)
         ENDIF
         CALL channel_halt(substr, status)
      END IF
      IF (lh) THEN
         ! horiz interpolated
         CALL new_channel_object(status,chname,'H_'//TRIM(ParData(idx)%name%obj)&
              , p4=ParData(idx)%ptr_hint, reprid=hreprid)
         IF (status == 3102) THEN ! channel object exists already
            CALL new_channel_object(status,chname             &
                 , TRIM(ParData(idx)%name%obj)//'_'//str(idx) &
                 , p4=ParData(idx)%ptr_hint, reprid=hreprid)
         ENDIF
         CALL channel_halt(substr, status)
      ELSE
         ParData(idx)%ptr_hint => ParData(idx)%ptr_int
      END IF

      ! DEFINE AXIS STRING
      CALL get_representation_info(status, ' ', reprId, rank=ParDATA(idx)%rank &
           , axis=ParDATA(idx)%AXIS, ldimlen=ParDATA(idx)%ldimlen)
      CALL channel_halt(substr, status)

      ParData(idx)%SD_ID = ParData(idx_int1)%SD_ID
      ParData(idx)%PSD   => ParData(idx_int1)%PSD

    END SUBROUTINE add_parent_list_element

  END SUBROUTINE parent_set_ParData
  !----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE match_parent_grid

    ! MESSy/BMIL
    USE messy_main_mpi_bi,       ONLY:  my_cart_id
    ! INT2COSMO
    USE data_int2lm_parallel,    ONLY: isubpos_in, nboundlines_in
    USE data_grid_in,            ONLY: ie_in_tot, je_in_tot                   &
                                     , startlon_in_tot, startlat_in_tot       &
                                     , ie_in, je_in &
                                     , dlon_in, dlat_in

    IMPLICIT NONE

    INTRINSIC :: REAL

    ! LOCAL
    INTEGER                     :: i,j         ! loop indices
    INTEGER                     :: itot, jtot
    CHARACTER(LEN=*), PARAMETER :: substr= 'match_parent_grid'

    ! NO PARENT FIELDS
    IF (PD_NUM == 0) RETURN

    ! 1.) locate parent grid boxes to that child Pe where it is
    ALLOCATE( rlons(ie_in_tot))
    ALLOCATE( rlats(je_in_tot))

    DO i = 1, ie_in
       itot = isubpos_in(my_cart_id,1)-nboundlines_in-1
       rlons(i) = startlon_in_tot + dlon_in * REAL(itot+i-1,dp)
    END DO
    DO j = 1, je_in
       jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
       rlats(j) = startlat_in_tot + dlat_in * REAL(jtot+j-1,dp)
    END DO

    CALL start_message_bi(submodstr, ' calculating weights', substr)
    CALL CALC_backward_WEIGHTS
    CALL end_message_bi(submodstr, ' calculating weights', substr)

  END SUBROUTINE MATCH_PARENT_GRID
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE CALC_backward_WEIGHTS

    ! MESSy/SMCL
    USE messy_main_grid,                 ONLY: GRID_ERROR, COPY_GEOHYBGRID     &
                                             , NEW_GEOHYBGRID
    USE messy_main_grid_netcdf,          ONLY: MAIN_GRID_SET_MESSAGEMODE       &
                                             , MSGMODE_E &
                                             !, MSGMODE_VL, MSGMODE_S &
                                             !, MSGMODE_W, MSGMODE_VM, MSGMODE_I&
                                             , INIT_NCVAR
    USE messy_main_grid_tools,           ONLY: RGTOOL_CONVERT     &
                                             , RGTOOL_CONVERT_DAT2VAR
    USE messy_main_grid_trafo,           ONLY: RG_INT, SWITCH_GEOHYBGRID
    USE messy_main_grid_trafo_scrp,      ONLY: CALC_SCRIPDATA     &
                                             , CALC_SCRIP_WEIGHTS &
                                             , SCRIP_CONTROL
    USE messy_main_grid_trafo_scrp_base, ONLY: norm_opt_dstarea

    ! INT2COSMO
    USE data_grid_in,                    ONLY: ke_in, ie_in, je_in

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'calc_backward_weights'
    INTEGER                     :: status
    LOGICAL                     :: lok
    INTEGER, DIMENSION(1)       :: RGTYPE
    INTEGER                     :: ii

    ! used for interpolation of horizontal weight function
    TYPE(t_scrip_data), POINTER           :: fw_PSD   => NULL()
    INTEGER                               :: fw_SD_ID = -99
    TYPE(t_ncvar), DIMENSION(:),  POINTER :: vari     => NULL()
    TYPE(t_ncvar), DIMENSION(:),  POINTER :: sovar    => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: dat      => NULL()
    INTEGER                               :: ix

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!    CALCULATE WEIGHTS   !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    CALL MAIN_GRID_SET_MESSAGEMODE(MSGMODE_S + MSGMODE_E + MSGMODE_VL &!)
!                                 + MSGMODE_W + MSGMODE_VM + MSGMODE_I)
    CALL MAIN_GRID_SET_MESSAGEMODE(MSGMODE_E)
    ! define reduced level number
    ke_int = ke_in-par_kmin+1
    ! 1a) DEFINE (TOTAL) CHILD GRID
    IF (L_gridrotParenteqChild) THEN
       CALL define_rotred_BM_grid(cgrid, lok)
    ELSE
       CALL define_reduced_BM_grid(cgrid, lok)
    ENDIF

    ! 1b) determine horizontal child grid
    CALL COPY_GEOHYBGRID(chgrid,cgrid)
    CALL SWITCH_GEOHYBGRID(chgrid, .TRUE., .TRUE., .FALSE.)

    ! 2a) DEFINE (TOTAL) PARENT GRID
    IF (L_gridrotParenteqChild) THEN
       CALL define_rotlocal_parent_grid(pgrid, lok)
    ELSE
       CALL define_local_parent_grid(pgrid, lok)
    ENDIF

    CALL NEW_GEOHYBGRID(status,pgrid_id,pgrid)

    ! 2b) DEFINE horizontal PARENT GRID
    CALL COPY_GEOHYBGRID(phgrid, pgrid)
    CALL SWITCH_GEOHYBGRID(phgrid, .TRUE., .TRUE., .FALSE.)
    ! 3.) convert to SCRIP GRID
    RGTYPE(1) = RG_INT ! dummy

    IF (loverlap) THEN
       DO ii= 1, PD_NUM
          CALL CALC_SCRIPDATA(status, cgrid, pgrid, RGTYPE, ParData(ii)%SD_ID &
               , PSD=ParData(ii)%PSD &
               , norm_opt_in=norm_opt_dstarea, map_type_in=ParData(ii)%interpM)
          IF (status /= 0) THEN
             IF (status /= 01 )  CALL error_bi(grid_error(status), substr)
          ELSE
             ! CALCULATE WEIGHTS
             CALL CALC_SCRIP_WEIGHTS(status, ParData(ii)%PSD)
             IF (fw_SD_ID < 0 .AND. ParData(ii)%interpM == 1) THEN
                fw_SD_ID = ParData(ii)%SD_ID
             END IF
          ENDIF
       END DO
    END IF

    !*************************************************************
    !*************************************************************
    !*************************************************************
    ! interpolate horizontal weight function
    !*************************************************************
    IF (.NOT. loverlap) THEN
       ALLOCATE(fw_int(ie_in, je_in,1,1))
       fw_int = 0._dp
       RETURN
    END IF

    IF (fw_SD_ID < 0) THEN ! map_type == 1 not yet calculated
       CALL CALC_SCRIPDATA(status, cgrid, pgrid, RGTYPE, fw_SD_ID    &
            , PSD=fw_PSD, norm_opt_in=norm_opt_dstarea, map_type_in=1)
       IF (status /= 0) THEN
          IF (status /= 01 )  RETURN
       ELSE
          ! CALCULATE WEIGHTS
          CALL CALC_SCRIP_WEIGHTS(status, fw_PSD)
       ENDIF
    END IF
    CALL MAIN_GRID_SET_MESSAGEMODE(MSGMODE_E)

    ALLOCATE(vari(1))

    dat => fw(iis:iie, jjs:jje,:,:)
    CALL RGTOOL_CONVERT_DAT2VAR(vari(1), dat, 'fw', chgrid, 'xyzn')
    NULLIFY(dat)

    CALL SCRIP_CONTROL(status, fw_SD_ID, cgrid, pgrid &
         , RGTYPE, .TRUE., vari, sovar, llrgz=.FALSE.)
    IF (status /= 0) CALL ERROR_BI(grid_error(status), substr)

    CALL RGTOOL_CONVERT(sovar(1), fw_int, phgrid, 'xyzn')

    ! CLEAN
    DO ix= 1, SIZE(vari)
       CALL INIT_NCVAR(vari(ix))
       CALL INIT_NCVAR(sovar(ix))
    END DO
    DEALLOCATE(vari)
    NULLIFY(vari)
    DEALLOCATE(sovar)
    NULLIFY(sovar)

    NULLIFY(fw_PSD)

    !*************************************************************
    !*************************************************************

 CONTAINS

   SUBROUTINE define_local_parent_grid(g,ok)
    !
    ! remove halo to avoid double counting of halo area

     ! INT2COSMO
     USE data_int2lm_parallel,    ONLY: isubpos_in, nboundlines_in, my_cart_id
     USE data_grid_in,            ONLY: startlon_in_tot, startlat_in_tot       &
                                      , dlon_in, dlat_in                       &
                                      , polgam_in, pollon_in, pollat_in        &
                                      , ie_in, je_in, ke_in, ak_in, bk_in
     ! MESSy/SMCL
     USE messy_main_timer,        ONLY: YEAR_START, MONTH_START, DAY_START     &
                                      , HOUR_START, MINUTE_START, SECOND_START &
                                      , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

     USE messy_main_constants_mem, ONLY: sp
     USE messy_main_timer,         ONLY: time_span_d

     ! GRID MODULES
     USE messy_main_grid_trafo,    ONLY: RGEMPTY, COMPLETE_GEOHYBGRID &
                                       , PHIROT2PHI, RLAROT2RLA
     USE messy_main_grid_netcdf,   ONLY: ERRMSG, COPY_NCDIM       &
                                       , null_dimid, null_varid   &
                                       , VTYPE_REAL, VTYPE_DOUBLE &
                                       , nf90_float, nf90_double  &
                                       , POSITION, INIT_NARRAY, ADD_NCATT
     USE messy_main_grid,          ONLY: INIT_GEOHYBGRID

     IMPLICIT NONE

     ! I/O
     TYPE(t_geohybgrid),  INTENT(OUT) :: g
     LOGICAL            , INTENT(OUT) :: ok

     !LOCAL
     REAL(DP)            :: dts
     INTEGER             :: i,j, n, itot, jtot
     INTEGER             :: status
     CHARACTER(LEN=100)  :: tunit
     REAL(dp)            :: clon, clat
     INTEGER             :: itmp
     INTEGER             :: istartlon, istartlat, idlon, idlat

     istartlon = NINT(startlon_in_tot * RCF_IN)
     istartlat = NINT(startlat_in_tot * RCF_IN)
     idlon     = NINT(dlon_in  * RCF_IN)
     idlat     = NINT(dlat_in  * RCF_IN)

     ! INIT
     CALL INIT_GEOHYBGRID(g)

     g%name = 'PARENT GRID'

     g%file = ' '       ! Filename
     g%t    = 0                               ! time step

     ! Curvilinear LONGITUDE (MID) ...
     g%clonm%name  = 'lon'
     g%clonm%id    = NULL_VARID
     g%clonm%xtype = NF90_DOUBLE
     g%clonc = .FALSE.
     ! ... dimensions
     g%clonm%ndims = 2
     ALLOCATE(g%clonm%dim(g%clonm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%clonm%dim(1)%name  = 'lon'
     g%clonm%dim(1)%id    = NULL_DIMID
     g%clonm%dim(1)%len   = ie_in
     g%clonm%dim(1)%fuid  = .false.
     g%clonm%dim(1)%varid = NULL_VARID
     g%clonm%dim(2)%name  = 'lat'
     g%clonm%dim(2)%id    = NULL_DIMID
     g%clonm%dim(2)%len   = je_in
     g%clonm%dim(2)%fuid  = .false.
     g%clonm%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%clonm%dat, g%clonm%ndims &
          , (/g%clonm%dim(1)%len,g%clonm%dim(2)%len/) &
          ,VTYPE_DOUBLE)

     !  CALCULATE clonm BELOW together with clatm

     ! ... attributes
     CALL ADD_NCATT(g%clonm, 'long_name', vs='longitude')
     CALL ADD_NCATT(g%clonm, 'units', vs='degrees_east')
     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%clatm%name  = 'lat'
     g%clatm%id    = NULL_VARID
     g%clatm%xtype = NF90_DOUBLE
     ! ... dimensions
     g%clatm%ndims = 2
     ALLOCATE(g%clatm%dim(g%clatm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%clatm%dim(1)%name  = 'lon'
     g%clatm%dim(1)%id    = NULL_DIMID
     g%clatm%dim(1)%len   = ie_in
     g%clatm%dim(1)%fuid  = .false.
     g%clatm%dim(1)%varid = NULL_VARID
     g%clatm%dim(2)%name  = 'lat'
     g%clatm%dim(2)%id    = NULL_DIMID
     g%clatm%dim(2)%len   = je_in
     g%clatm%dim(2)%fuid  = .false.
     g%clatm%dim(2)%varid = NULL_VARID

     ! ... attributes
     CALL ADD_NCATT(g%clatm, 'long_name', vs='mid latitude')
     CALL ADD_NCATT(g%clatm, 'units', vs='degrees_north')

     g%ranges(2,1) = RGEMPTY
     g%ranges(2,2) = RGEMPTY

     ! ... data
     CALL INIT_NARRAY(g%clatm%dat, g%clatm%ndims &
          , (/g%clatm%dim(1)%len, g%clatm%dim(2)%len/) &
!          ,VTYPE_REAL)
          ,VTYPE_DOUBLE)
     DO i=1, ie_in
        itot = isubpos_in(my_cart_id,1)-nboundlines_in-1
        itmp = istartlon + (itot+i-1) *  idlon
        clon = REAL(itmp,dp) / RCF_IN
        DO j = 1 , je_in
           jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
           itmp = istartlat + (jtot+j-1) * idlat
           clat = REAL(itmp,dp) / RCF_IN
           n = (j-1) * g%clonm%dim(1)%len + i
           g%clonm%dat%vd(n) = &
                rlarot2rla(clat, clon, pollat_in, pollon_in, polgam_in)
           g%clatm%dat%vd(n) = &
                 phirot2phi(clat, clon, pollat_in, polgam_in)
        END DO
     END DO

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     g%cloni%name  = 'lon_I'
     g%cloni%id    = NULL_VARID
     g%cloni%xtype = NF90_DOUBLE
     ! ... dimensions
     g%cloni%ndims = 2
     ALLOCATE(g%cloni%dim(g%cloni%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%cloni%dim(1)%name  = 'lon_I'
     g%cloni%dim(1)%id    = NULL_DIMID
     g%cloni%dim(1)%len   = ie_in+1
     g%cloni%dim(1)%fuid  = .false.
     g%cloni%dim(1)%varid = NULL_VARID
     g%cloni%dim(2)%name  = 'lat_I'
     g%cloni%dim(2)%id    = NULL_DIMID
     g%cloni%dim(2)%len   = je_in+1
     g%cloni%dim(2)%fuid  = .false.
     g%cloni%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%cloni%dat, g%cloni%ndims &
          , (/g%cloni%dim(1)%len,g%cloni%dim(2)%len/) &
          ,VTYPE_DOUBLE)

     ! ... attributes
     CALL ADD_NCATT(g%cloni, 'long_name', vs='interface longitude')
     CALL ADD_NCATT(g%cloni, 'units', vs='degrees_east')

     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! CALCULATE CLONI together with clati below

     ! LATITUDE (MID) ...
     g%clati%name  = 'lat_I'
     g%clati%id    = NULL_VARID
     g%clati%xtype = NF90_DOUBLE
     ! ... dimensions
     g%clati%ndims = 2
     ALLOCATE(g%clati%dim(g%clati%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%clati%dim(1)%name  = 'lon_I'
     g%clati%dim(1)%id    = NULL_DIMID
     g%clati%dim(1)%len   = ie_in+1
     g%clati%dim(1)%fuid  = .false.
     g%clati%dim(1)%varid = NULL_VARID
     g%clati%dim(2)%name  = 'lat_I'
     g%clati%dim(2)%id    = NULL_DIMID
     g%clati%dim(2)%len   = je_in+1
     g%clati%dim(2)%fuid  = .false.
     g%clati%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%clati%dat, g%clati%ndims &
          , (/g%clati%dim(1)%len, g%clati%dim(2)%len/), VTYPE_DOUBLE)
     DO i=1, ie_in+1
        itot = isubpos_in(my_cart_id,1)-nboundlines_in-1
        itmp = istartlon + (itot+i) * idlon - NINT(1.5_dp * dlon_in * RCF_IN)
        clon = REAL(itmp,dp) / RCF_IN
        DO j = 1 , je_in+1
           jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
           itmp = istartlat + (jtot+j) * idlat - NINT(1.5_dp * dlat_in * RCF_IN)
           clat = REAL(itmp,dp) / RCF_IN
           n = (j-1) * g%clati%dim(1)%len + i
           g%clati%dat%vd(n) = &
                phirot2phi(clat, clon, pollat_in, polgam_in)
           g%cloni%dat%vd(n) = &
                rlarot2rla(clat, clon, pollat_in, pollon_in, polgam_in)
        END DO
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%clati, 'long_name', vs='interface latitude')
     CALL ADD_NCATT(g%clati, 'units', vs='degrees_north')

     ! .. set ranges
     g%ranges(2,1) = RGEMPTY !-90.0_sp
     g%ranges(2,2) = RGEMPTY

     ! DETERMIN MIN / MAX LATITUDE
     g%minmaxlonlat(1,1) = MINVAL(g%cloni%dat%vd)
     g%minmaxlonlat(1,2) = MAXVAL(g%cloni%dat%vd)

     g%minmaxlonlat(2,1) = MINVAL(g%clati%dat%vd)
     g%minmaxlonlat(2,2) = MAXVAL(g%clati%dat%vd)

     !**************************************************************************
     ! Rotated  Curvilinear LONGITUDE (MID) ...
     g%rlonm%name  = 'lon'
     g%rlonm%id    = NULL_VARID
     g%rlonm%xtype = NF90_DOUBLE
     g%rlonc = .FALSE.
     ! ... dimensions
     g%rlonm%ndims = 1
     ALLOCATE(g%rlonm%dim(g%rlonm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%rlonm%dim(1)%name  = 'lon'
     g%rlonm%dim(1)%id    = NULL_DIMID
     g%rlonm%dim(1)%len   = ie_in
     g%rlonm%dim(1)%fuid  = .false.

     ! ... data
     CALL INIT_NARRAY(g%rlonm%dat, g%rlonm%ndims &
          , (/g%rlonm%dim(1)%len/),VTYPE_DOUBLE)
     DO i=1, ie_in
        itot = isubpos_in(my_cart_id,1)-nboundlines_in-1
        itmp = istartlon + (itot+i-1) * idlon
        g%rlonm%dat%vd(i) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rlonm, 'long_name', vs='longitude')
     CALL ADD_NCATT(g%rlonm, 'units', vs='degrees_east')
     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%rlatm%name  = 'lat'
     g%rlatm%id    = NULL_VARID
     g%rlatm%xtype = NF90_DOUBLE
     ! ... dimensions
     g%rlatm%ndims = 1
     ALLOCATE(g%rlatm%dim(g%rlatm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%rlatm%dim(1)%name  = 'lat'
     g%rlatm%dim(1)%id    = NULL_DIMID
     g%rlatm%dim(1)%len   = je_in
     g%rlatm%dim(1)%fuid  = .false.
     g%rlatm%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%rlatm%dat, g%rlatm%ndims &
       , (/g%rlatm%dim(1)%len/),VTYPE_DOUBLE)

     DO j=1,je_in
        jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
        itmp = istartlat + (jtot+j-1) * idlat
        g%rlatm%dat%vd(j) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rlatm, 'long_name', vs='mid latitude')
     CALL ADD_NCATT(g%rlatm, 'units', vs='degrees_north')

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     g%rloni%name  = 'lon_I'
     g%rloni%id    = NULL_VARID
     g%rloni%xtype = NF90_DOUBLE
     ! ... dimensions
     g%rloni%ndims = 1
     ALLOCATE(g%rloni%dim(g%rloni%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%rloni%dim(1)%name  = 'lon_I'
     g%rloni%dim(1)%id    = NULL_DIMID
     g%rloni%dim(1)%len   = ie_in+1
     g%rloni%dim(1)%fuid  = .false.
     g%rloni%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%rloni%dat, g%rloni%ndims &
          , (/g%rloni%dim(1)%len/)  ,VTYPE_DOUBLE)

     DO i=1, ie_in+1
       itot = isubpos_in(my_cart_id,1)-nboundlines_in-1
       itmp = istartlon +  (itot+i) *idlon - NINT(1.5_dp * dlon_in * RCF_IN)
       g%rloni%dat%vd(i) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rloni, 'long_name', vs='interface longitude')
     CALL ADD_NCATT(g%rloni, 'units', vs='degrees_east')

     ! LATITUDE (MID) ...
     g%rlati%name  = 'lat_I'
     g%rlati%id    = NULL_VARID
     g%rlati%xtype = NF90_DOUBLE
     ! ... dimensions
     g%rlati%ndims = 1
     ALLOCATE(g%rlati%dim(g%rlati%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%rlati%dim(1)%name  = 'lat_I'
     g%rlati%dim(1)%id    = NULL_DIMID
     g%rlati%dim(1)%len   = je_in+1
     g%rlati%dim(1)%fuid  = .false.
     g%rlati%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%rlati%dat, g%rlati%ndims &
          , (/g%rlati%dim(1)%len/) ,VTYPE_DOUBLE)

     DO j = 1 , je_in+1
        jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
        itmp = istartlat + (jtot+j) * idlat - NINT(1.5_dp * dlat_in * RCF_IN)
        g%rlati%dat%vd(j) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rlati, 'long_name', vs='interface latitude')
     CALL ADD_NCATT(g%rlati, 'units', vs='degrees_north')

     ! *************************************************************************
     ! HYBRID-A-COEFFICIENTS (INTERFACES) ...
     g%hyai%name  = 'hyai'
     g%hyai%id    = NULL_VARID
     g%hyai%xtype = NF90_FLOAT
     ! ... dimensions
     g%hyai%ndims = 1
     ALLOCATE(g%hyai%dim(g%hyai%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
     g%hyai%dim(1)%name  = 'ilev'
     g%hyai%dim(1)%id    = NULL_DIMID
     g%hyai%dim(1)%len   = ke_int +1
     g%hyai%dim(1)%fuid  = .false.
     g%hyai%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%hyai%dat, g%hyai%ndims &
          , (/g%hyai%dim(1)%len/),VTYPE_REAL)
     g%hyai%dat%vr = REAL(ak_in(par_kmin:ke_in+1),sp)

     ! ... attributes
     CALL ADD_NCATT(g%hyai, 'long_name'                  &
          ,vs='hybrid-A-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hyai, 'units', vs='1')
     g%ranges(3,1) = RGEMPTY
     g%ranges(3,2) = RGEMPTY

     ! ************************************************************************
     ! HYBRID-B-COEFFICIENTS (INTERFACES) ...
     g%hybi%name  = 'hybi'
     g%hybi%id    = NULL_VARID
     g%hybi%xtype = NF90_FLOAT
     ! ... dimensions
     g%hybi%ndims = 1
     ALLOCATE(g%hybi%dim(g%hybi%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
     CALL COPY_NCDIM(g%hybi%dim(1),g%hyai%dim(1))
     ! ... data
     CALL INIT_NARRAY(g%hybi%dat, g%hybi%ndims &
          , (/g%hybi%dim(1)%len/),VTYPE_REAL)

     g%hybi%dat%vr = REAL(bk_in(par_kmin:ke_in+1),sp)
     ! ... attributes
     CALL ADD_NCATT(g%hybi, 'long_name'                  &
          ,vs='hybrid-B-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hybi, 'units', vs='1')

     ! ranges
     g%ranges(4,1) = 0.0_sp
     g%ranges(4,2) = 1.0_sp

     ! ************************************************************************
     ! SURFACE PRESSURE
     g%ps%name  = 'ps'
     g%ps%id    = NULL_VARID
     g%ps%xtype = NF90_FLOAT
     ! ... dimensions
     g%ps%ndims = 2
     ALLOCATE(g%ps%dim(g%ps%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

     CALL COPY_NCDIM(g%ps%dim(1), g%rlonm%dim(1))
     CALL COPY_NCDIM(g%ps%dim(2), g%rlatm%dim(1))
     ! ... data
     CALL INIT_NARRAY(g%ps%dat, g%ps%ndims      &
          ,(/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
          ,VTYPE_REAL)
     DO i=1, g%ps%dim(1)%len
        DO j=1, g%ps%dim(2)%len
           g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
                ,(/i,j/)))  = 101325.0_sp
        END DO
     END DO
     ! ... attributes
     CALL ADD_NCATT(g%ps, 'long_name', vs='surface pressure')
     CALL ADD_NCATT(g%ps, 'units', vs='Pa')

     ! **************************************************************************
     ! REFERENCE PRESSURE ...
     g%p0%name  = 'p0'
     g%p0%id    = NULL_VARID
     g%p0%xtype = NF90_FLOAT
     ! ... dimensions
     g%p0%ndims = 0
     ! ... data
     CALL INIT_NARRAY(g%p0%dat, 1, (/ 1 /), VTYPE_REAL)
     g%p0%dat%vr(1) = 1.0_sp
     ! ... attributes
     CALL ADD_NCATT(g%p0, 'long_name', vs='reference pressure')
     CALL ADD_NCATT(g%p0, 'units', vs='Pa')

     ! TIME (MID) ...
     g%timem%name  = 'time'
     g%timem%id    = NULL_VARID
     g%timem%xtype = NF90_DOUBLE
     ! ... dimensions
     g%timem%ndims = 1
     ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
     g%timem%dim(1)%name  = 'time'
     g%timem%dim(1)%id    = NULL_DIMID
     g%timem%dim(1)%len   = 1
     g%timem%dim(1)%fuid  = .true.
     g%timem%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
          ,VTYPE_DOUBLE)
     ! TIME: SECONDS SINCE MODEL START
     CALL time_span_d(dts   &
          , YEAR_START, MONTH_START, DAY_START &
          , HOUR_START, MINUTE_START, SECOND_START  &
          , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
     g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds

     ! ... attributes
     CALL ADD_NCATT(g%timem, 'long_name'                              &
          ,vs='time in seconds since model start')
     WRITE(tunit, &
          '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
          YEAR_START, MONTH_START, DAY_START &
          , HOUR_START, MINUTE_START, SECOND_START

     CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

     ! CALCULATE INTs from MIDs
     CALL COMPLETE_GEOHYBGRID(g)

     ! INTERFACE OK
     ok = .true.

   END SUBROUTINE define_local_parent_grid

   ! ------------------------------------------------------------------
   SUBROUTINE DEFINE_REDUCED_BM_GRID(g, ok)

    !
    ! remove halo to avoid double counting of halo area

    ! MESSY/BMIL
     USE messy_main_mpi_bi,       ONLY: isubpos, my_cart_id, nboundlines
     USE messy_main_grid_def_mem_bi, ONLY: nlev => ke_tot &
                                      , startlon_tot, startlat_tot             &
                                      , dlon, dlat, polgam, pollon, pollat     &
                                      , istart, iend, jstart, jend
     USE messy_main_grid_def_bi,  ONLY:hyai, hybi
     USE messy_main_data_bi,      ONLY: aps

     USE messy_main_timer,        ONLY: YEAR_START, MONTH_START, DAY_START     &
                                      , HOUR_START, MINUTE_START, SECOND_START &
                                      , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

     ! MESSy
     USE messy_main_constants_mem, ONLY: sp
     USE messy_main_timer,         ONLY: time_span_d

     ! NCREGRID MODULES
     USE messy_main_grid_trafo,  ONLY: RGEMPTY, COMPLETE_GEOHYBGRID &
                                     , PHIROT2PHI, RLAROT2RLA
     USE messy_main_grid_netcdf, ONLY: ERRMSG, COPY_NCDIM       &
                                     , null_dimid, null_varid   &
                                     , vtype_real, vtype_double &
                                     , nf90_float, nf90_double  &
                                     , POSITION, INIT_NARRAY, ADD_NCATT
     USE messy_main_grid,        ONLY: INIT_GEOHYBGRID

     IMPLICIT NONE

     ! I/O
     TYPE(t_geohybgrid), INTENT(OUT) :: g
     LOGICAL           , INTENT(OUT) :: ok

     !LOCAL
     REAL(DP)             :: dts
     INTEGER              :: i,j, n, id , jd, itot, jtot
     INTEGER              :: status
     CHARACTER(LEN=100)   :: tunit
     REAL(dp)             :: clon, clat
     INTEGER              :: itmp
     LOGICAL              :: l_ijf
     INTEGER              :: clonlen = 0
     INTEGER              :: clatlen = 0
     INTEGER             :: istartlon, istartlat, idlon, idlat

     istartlon = NINT(startlon_tot * RCF)
     istartlat = NINT(startlat_tot * RCF)
     idlon     = NINT(dlon * RCF)
     idlat     = NINT(dlat * RCF)

     ! DEFINITION OF START, END AND LENGTH OF REDUCED CHILD GRID
     ! start and end of inner region
     ! define locally to be able to exclude damping regions etc.

     IF(i_rmy_px < 0) THEN
        iis     = istart
        iie     = iend
        jjs     = jstart
        jje     = jend
     ELSE
        l_ijf = .FALSE.
        DO i=istart,iend
           DO j=jstart,jend
              IF(rmy_loc(i,j) <= 0._dp) THEN
                 iis=i
                 jjs=j
                 l_ijf = .TRUE.
                 exit
              END IF
           END DO
           IF(l_ijf) exit
        END DO
        IF(.NOT. l_ijf) THEN
           iis=0
           iie=-1
           jjs=0
           jje=-1
        ELSE
           l_ijf = .FALSE.
           DO i=iend,istart,-1
              DO j=jend,jstart,-1
                 IF(rmy_loc(i,j) <= 0._dp) THEN
                    iie=i
                    jje=j
                    l_ijf = .TRUE.
                    exit
                 END IF
              END DO
              IF(l_ijf) exit
           END DO
        END IF
        DEALLOCATE(rmy_loc)
     ENDIF
     clonlen = iie - iis + 1
     clatlen = jje - jjs + 1
     IF (clonlen == 0 .OR. clatlen == 0) THEN
        write(iouerr,*) 'CLON / CLAT no overlap ', clonlen , clatlen
        loverlap = .FALSE.
     ENDIF

     ! INIT
     CALL INIT_GEOHYBGRID(g)

     g%name = 'COSMOohalo'

     g%file = ' '   ! Filename
     g%t    = 0                               ! time step

     ! Curvilinear LONGITUDE (MID) ...
     g%clonm%name  = 'lon'
     g%clonm%id    = NULL_VARID
     g%clonm%xtype = NF90_DOUBLE
     g%clonc = .FALSE.
     ! ... dimensions
     g%clonm%ndims = 2
     ALLOCATE(g%clonm%dim(g%clonm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%clonm%dim(1)%name  = 'lon'
     g%clonm%dim(1)%id    = NULL_DIMID
     g%clonm%dim(1)%len   = clonlen
     g%clonm%dim(1)%fuid  = .false.
     g%clonm%dim(1)%varid = NULL_VARID
     g%clonm%dim(2)%name  = 'lat'
     g%clonm%dim(2)%id    = NULL_DIMID
     g%clonm%dim(2)%len   = clatlen
     g%clonm%dim(2)%fuid  = .false.
     g%clonm%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%clonm%dat, g%clonm%ndims      &
          , (/g%clonm%dim(1)%len,g%clonm%dim(2)%len/) &
          , VTYPE_DOUBLE)

     ! ... attributes
     CALL ADD_NCATT(g%clonm, 'long_name', vs='longitude')
     CALL ADD_NCATT(g%clonm, 'units', vs='degrees_east')
     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%clatm%name  = 'lat'
     g%clatm%id    = NULL_VARID
     g%clatm%xtype = NF90_DOUBLE
     ! ... dimensions
     g%clatm%ndims = 2
     ALLOCATE(g%clatm%dim(g%clatm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%clatm%dim(1)%name  = 'lon'
     g%clatm%dim(1)%id    = NULL_DIMID
     g%clatm%dim(1)%len   = clonlen
     g%clatm%dim(1)%fuid  = .false.
     g%clatm%dim(1)%varid = NULL_VARID
     g%clatm%dim(2)%name  = 'lat'
     g%clatm%dim(2)%id    = NULL_DIMID
     g%clatm%dim(2)%len   = clatlen
     g%clatm%dim(2)%fuid  = .false.
     g%clatm%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%clatm%dat, g%clatm%ndims       &
          , (/g%clatm%dim(1)%len, g%clatm%dim(2)%len/) &
          , VTYPE_DOUBLE)
     DO i = iis, iie
        itot = isubpos(my_cart_id,1)-nboundlines-1
        itmp = istartlon + (itot+i-1) * idlon
        clon = REAL(itmp,dp) / RCF
        DO j= jjs, jje
           id = i - iis + 1
           jd = j - jjs + 1
           n = (jd-1) * g%clatm%dim(1)%len + id
           jtot = isubpos(my_cart_id,2)-nboundlines-1
           itmp = istartlat  + (jtot+j-1) * idlat
           clat = REAL(itmp,dp) / RCF
           g%clonm%dat%vd(n) = &
                rlarot2rla(clat, clon, pollat, pollon, polgam)
           g%clatm%dat%vd(n) = &
                phirot2phi(clat, clon, pollat, polgam)
        END DO
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%clatm, 'long_name', vs='mid latitude')
     CALL ADD_NCATT(g%clatm, 'units', vs='degrees_north')

     g%ranges(2,1) = RGEMPTY !-90.0_dp
     g%ranges(2,2) = RGEMPTY

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     g%cloni%name  = 'lon_I'
     g%cloni%id    = NULL_VARID
     g%cloni%xtype = NF90_DOUBLE
     ! ... dimensions
     g%cloni%ndims = 2
     ALLOCATE(g%cloni%dim(g%cloni%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%cloni%dim(1)%name  = 'lon_I'
     g%cloni%dim(1)%id    = NULL_DIMID
     g%cloni%dim(1)%len   = clonlen + 1
     g%cloni%dim(1)%fuid  = .false.
     g%cloni%dim(1)%varid = NULL_VARID
     g%cloni%dim(2)%name  = 'lat_I'
     g%cloni%dim(2)%id    = NULL_DIMID
     g%cloni%dim(2)%len   = clatlen + 1
     g%cloni%dim(2)%fuid  = .false.
     g%cloni%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%cloni%dat, g%cloni%ndims      &
          , (/g%cloni%dim(1)%len,g%cloni%dim(2)%len/) &
          , VTYPE_DOUBLE)

     ! ... attributes
     CALL ADD_NCATT(g%cloni, 'long_name', vs='interface longitude')
     CALL ADD_NCATT(g%cloni, 'units', vs='degrees_east')

     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%clati%name  = 'lat_I'
     g%clati%id    = NULL_VARID
     g%clati%xtype = NF90_DOUBLE
     ! ... dimensions
     g%clati%ndims = 2
     ALLOCATE(g%clati%dim(g%clati%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%clati%dim(1)%name  = 'lon_I'
     g%clati%dim(1)%id    = NULL_DIMID
     g%clati%dim(1)%len   = clonlen + 1
     g%clati%dim(1)%fuid  = .false.
     g%clati%dim(1)%varid = NULL_VARID
     g%clati%dim(2)%name  = 'lat_I'
     g%clati%dim(2)%id    = NULL_DIMID
     g%clati%dim(2)%len   = clatlen + 1
     g%clati%dim(2)%fuid  = .false.
     g%clati%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%clati%dat, g%clati%ndims       &
          , (/g%clati%dim(1)%len, g%clati%dim(2)%len/) &
          , VTYPE_DOUBLE)

     DO i=iis, iie+1
        itot = isubpos(my_cart_id,1)-nboundlines-1
        itmp = istartlon + (itot+i) * idlon - NINT(1.5_dp * dlon * RCF)
        clon = REAL(itmp,dp) / RCF
        DO j = jjs, jje+1
           jtot = isubpos(my_cart_id,2)-nboundlines-1
           itmp = istartlat + (jtot+j) * idlat- NINT(1.5_dp * dlat * RCF)
           clat = REAL(itmp,dp) / RCF
           id = i - iis + 1
           jd = j - jjs + 1
           n = (jd-1) * g%cloni%dim(1)%len + id
           g%cloni%dat%vd(n) = &
                rlarot2rla(clat, clon, pollat, pollon, polgam)
           g%clati%dat%vd(n) = &
                phirot2phi(clat, clon, pollat, polgam)
        END DO
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%clati, 'long_name', vs='interface latitude')
     CALL ADD_NCATT(g%clati, 'units', vs='degrees_north')

     g%ranges(2,1) = RGEMPTY
     g%ranges(2,2) = RGEMPTY

     ! DETERMIN MIN / MAX LATITUDE
     g%minmaxlonlat(1,1) = MINVAL(g%cloni%dat%vd)
     g%minmaxlonlat(1,2) = MAXVAL(g%cloni%dat%vd)

     g%minmaxlonlat(2,1) = MINVAL(g%clati%dat%vd)
     g%minmaxlonlat(2,2) = MAXVAL(g%clati%dat%vd)

     !**************************************************************************
     ! Rotated  Curvilinear LONGITUDE (MID) ...
     g%rlonm%name  = 'lon'
     g%rlonm%id    = NULL_VARID
     g%rlonm%xtype = NF90_DOUBLE
     g%rlonc = .FALSE.
     ! ... dimensions
     g%rlonm%ndims = 1
     ALLOCATE(g%rlonm%dim(g%rlonm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%rlonm%dim(1)%name  = 'lon'
     g%rlonm%dim(1)%id    = NULL_DIMID
     g%rlonm%dim(1)%len   = clonlen
     g%rlonm%dim(1)%fuid  = .false.

     ! ... data
     CALL INIT_NARRAY(g%rlonm%dat, g%rlonm%ndims &
          , (/g%rlonm%dim(1)%len/),VTYPE_DOUBLE)

     DO i=iis, iie
        itot = isubpos(my_cart_id,1)-nboundlines-1
        itmp = istartlon + (itot+i-1) * idlon
        g%rlonm%dat%vd(i-iis+1) = REAL(itmp,dp) / RCF
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rlonm, 'long_name', vs='longitude')
     CALL ADD_NCATT(g%rlonm, 'units', vs='degrees_east')
     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%rlatm%name  = 'lat'
     g%rlatm%id    = NULL_VARID
     g%rlatm%xtype = NF90_DOUBLE
     ! ... dimensions
     g%rlatm%ndims = 1
     ALLOCATE(g%rlatm%dim(g%rlatm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%rlatm%dim(1)%name  = 'lat'
     g%rlatm%dim(1)%id    = NULL_DIMID
     g%rlatm%dim(1)%len   = clatlen
     g%rlatm%dim(1)%fuid  = .false.
     g%rlatm%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%rlatm%dat, g%rlatm%ndims &
          , (/g%rlatm%dim(1)%len/), VTYPE_DOUBLE)

     DO j=jjs,jje
        jtot = isubpos(my_cart_id,2)-nboundlines-1
        itmp = istartlat + (jtot+j-1) * idlat
        g%rlatm%dat%vd(j-jjs+1) = REAL(itmp,dp) / RCF
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rlatm, 'long_name', vs='mid latitude')
     CALL ADD_NCATT(g%rlatm, 'units', vs='degrees_north')

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     g%rloni%name  = 'lon_I'
     g%rloni%id    = NULL_VARID
     g%rloni%xtype = NF90_DOUBLE
     ! ... dimensions
     g%rloni%ndims = 1
     ALLOCATE(g%rloni%dim(g%rloni%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%rloni%dim(1)%name  = 'lon_I'
     g%rloni%dim(1)%id    = NULL_DIMID
     g%rloni%dim(1)%len   = clonlen+1
     g%rloni%dim(1)%fuid  = .false.
     g%rloni%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%rloni%dat, g%rloni%ndims &
          , (/g%rloni%dim(1)%len/)  ,VTYPE_DOUBLE)

     DO i=iis, iie+1
        itot = isubpos(my_cart_id,1)-nboundlines-1
        itmp = istartlon + (itot+i) * idlon - NINT(1.5_dp * dlon * RCF)
        g%rloni%dat%vd(i-iis+1) = REAL(itmp,dp) / RCF
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rloni, 'long_name', vs='interface longitude')
     CALL ADD_NCATT(g%rloni, 'units', vs='degrees_east')

     ! LATITUDE (MID) ...
     g%rlati%name  = 'lat_I'
     g%rlati%id    = NULL_VARID
     g%rlati%xtype = NF90_DOUBLE
     ! ... dimensions
     g%rlati%ndims = 1
     ALLOCATE(g%rlati%dim(g%rlati%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%rlati%dim(1)%name  = 'lat_I'
     g%rlati%dim(1)%id    = NULL_DIMID
     g%rlati%dim(1)%len   = clatlen+1
     g%rlati%dim(1)%fuid  = .false.
     g%rlati%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%rlati%dat, g%rlati%ndims &
          , (/g%rlati%dim(1)%len/) ,VTYPE_DOUBLE)

     DO j = jjs , jje+1
        jtot = isubpos(my_cart_id,2)-nboundlines-1
        itmp = istartlat +  (jtot+j) * idlat - NINT(1.5_dp * dlat * RCF)
        g%rlati%dat%vd(j-jjs+1) = REAL(itmp,dp) / RCF
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rlati, 'long_name', vs='interface latitude')
     CALL ADD_NCATT(g%rlati, 'units', vs='degrees_north')

     ! *************************************************************************
     ! HYBRID-A-COEFFICIENTS (INTERFACES) ...
     g%hyai%name  = 'hyai'
     g%hyai%id    = NULL_VARID
     g%hyai%xtype = NF90_FLOAT
     ! ... dimensions
     g%hyai%ndims = 1
     ALLOCATE(g%hyai%dim(g%hyai%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
     g%hyai%dim(1)%name  = 'ilev'
     g%hyai%dim(1)%id    = NULL_DIMID
     g%hyai%dim(1)%len   = nlev+1-chi_kmin+1
     g%hyai%dim(1)%fuid  = .false.
     g%hyai%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%hyai%dat, g%hyai%ndims, (/g%hyai%dim(1)%len/) &
          ,VTYPE_REAL)
     g%hyai%dat%vr(:) = REAL(hyai(chi_kmin:nlev+1),sp)

     ! ... attributes
     CALL ADD_NCATT(g%hyai, 'long_name'                              &
          ,vs='hybrid-A-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hyai, 'units', vs='1')
     g%ranges(3,1) = RGEMPTY
     g%ranges(3,2) = RGEMPTY

     ! *************************************************************************
     ! HYBRID-B-COEFFICIENTS (INTERFACES) ...
     g%hybi%name  = 'hybi'
     g%hybi%id    = NULL_VARID
     g%hybi%xtype = NF90_FLOAT
     ! ... dimensions
     g%hybi%ndims = 1
     ALLOCATE(g%hybi%dim(g%hybi%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
     CALL COPY_NCDIM(g%hybi%dim(1),g%hyai%dim(1))
     ! ... data
     CALL INIT_NARRAY(g%hybi%dat, g%hybi%ndims, (/g%hybi%dim(1)%len/) &
          ,VTYPE_REAL)
     g%hybi%dat%vr(:) = REAL(hybi(chi_kmin:nlev+1),sp)

     ! ... attributes
     CALL ADD_NCATT(g%hybi, 'long_name'                               &
          ,vs='hybrid-B-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hybi, 'units', vs='1')
     g%ranges(4,1) = 0.0_sp
     g%ranges(4,2) = 1.0_sp

     ! *************************************************************************
     ! SURFACE PRESSURE
     g%ps%name  = 'ps'
     g%ps%id    = NULL_VARID
     g%ps%xtype = NF90_FLOAT
     ! ... dimensions
     g%ps%ndims = 2
     ALLOCATE(g%ps%dim(g%ps%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

     g%ps%dim(1) = g%rlonm%dim(1)
     g%ps%dim(2) = g%rlatm%dim(1)

     ! ... data
     CALL INIT_NARRAY(g%ps%dat, g%ps%ndims      &
          ,(/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
          ,VTYPE_REAL)
     IF (ASSOCIATED(aps)) THEN
        DO i=1, g%ps%dim(1)%len
           DO j=1, g%ps%dim(2)%len
              g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                   ,(/i,j/)))                           &
                   !           = aps(i,j)
                   = 101325.0_sp
           END DO
        END DO
     ELSE
        DO i=1, g%ps%dim(1)%len
           DO j=1, g%ps%dim(2)%len
              g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                   ,(/i,j/)))                           &
                   = 101325.0_sp
           END DO
        END DO
     END IF
     ! ... attributes
     CALL ADD_NCATT(g%ps, 'long_name', vs='surface pressure')
     CALL ADD_NCATT(g%ps, 'units', vs='Pa')

     ! *************************************************************************
     ! REFERENCE PRESSURE ...
     g%p0%name  = 'p0'
     g%p0%id    = NULL_VARID
     g%p0%xtype = NF90_FLOAT
     ! ... dimensions
     g%p0%ndims = 0
     ! ... data
     CALL INIT_NARRAY(g%p0%dat, 1, (/ 1 /), VTYPE_REAL)
     g%p0%dat%vr(1) = 1.0_sp
     ! ... attributes
     CALL ADD_NCATT(g%p0, 'long_name', vs='reference pressure')
     CALL ADD_NCATT(g%p0, 'units', vs='Pa')

     ! TIME (MID) ...
     g%timem%name  = 'time'
     g%timem%id    = NULL_VARID
     g%timem%xtype = NF90_DOUBLE
     ! ... dimensions
     g%timem%ndims = 1
     ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
     g%timem%dim(1)%name  = 'time'
     g%timem%dim(1)%id    = NULL_DIMID
     g%timem%dim(1)%len   = 1
     g%timem%dim(1)%fuid  = .true.
     g%timem%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
          ,VTYPE_DOUBLE)
     ! TIME: SECONDS SINCE MODEL START
     CALL time_span_d(dts   &
          , YEAR_START, MONTH_START, DAY_START &
          , HOUR_START, MINUTE_START, SECOND_START  &
          , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
     g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds

     ! ... attributes
     CALL ADD_NCATT(g%timem, 'long_name'                              &
          ,vs='time in seconds since model start')
     WRITE(tunit, &
          '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
          YEAR_START, MONTH_START, DAY_START &
          , HOUR_START, MINUTE_START, SECOND_START

     CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

     ! CALCULATE INTs from MIDs
     CALL COMPLETE_GEOHYBGRID(g)

     ! INTERFACE OK
     ok = .true.

   END SUBROUTINE DEFINE_REDUCED_BM_GRID

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   SUBROUTINE define_rotlocal_parent_grid(g,ok)
    !
    ! remove halo to avoid double counting of halo area

     ! INT2COSMO
     USE data_int2lm_parallel,    ONLY: isubpos_in, nboundlines_in, my_cart_id
     USE data_grid_in,            ONLY: startlon_in_tot, startlat_in_tot       &
                                      , dlon_in, dlat_in                       &
                                      , ie_in, je_in, ke_in, ak_in, bk_in
     ! MESSy/SMCL
     USE messy_main_timer,        ONLY: YEAR_START, MONTH_START, DAY_START     &
                                      , HOUR_START, MINUTE_START, SECOND_START &
                                      , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

     USE messy_main_constants_mem, ONLY: sp
     USE messy_main_timer,         ONLY: time_span_d

     ! GRID MODULES
     USE messy_main_grid_trafo,    ONLY: RGEMPTY, COMPLETE_GEOHYBGRID
     USE messy_main_grid_netcdf,   ONLY: ERRMSG, COPY_NCDIM       &
                                       , null_dimid, null_varid   &
                                       , VTYPE_REAL, VTYPE_DOUBLE &
                                       , nf90_float, nf90_double  &
                                       , POSITION, INIT_NARRAY, ADD_NCATT
     USE messy_main_grid,          ONLY: INIT_GEOHYBGRID

     IMPLICIT NONE

     ! I/O
     TYPE(t_geohybgrid),  INTENT(OUT) :: g
     LOGICAL            , INTENT(OUT) :: ok

     !LOCAL
     REAL(DP)            :: dts
     INTEGER             :: i,j, itot, jtot
     INTEGER             :: status
     CHARACTER(LEN=100)  :: tunit
     INTEGER             :: itmp
     INTEGER             :: idlon, idlat, istartlon, istartlat

     idlon = NINT(dlon_in * RCF_IN)
     idlat = NINT(dlat_in * RCF_IN)
     istartlon = NINT(startlon_in_tot *RCF_IN)
     istartlat = NINT(startlat_in_tot *RCF_IN)

     ! INIT
     CALL INIT_GEOHYBGRID(g)

     g%name = 'PARENT GRID'

     g%file = ' '       ! Filename
     g%t    = 0         ! time step

     !**************************************************************************
     ! Rotated  Curvilinear LONGITUDE (MID) ...
     g%lonm%name  = 'lon'
     g%lonm%id    = NULL_VARID
     g%lonm%xtype = NF90_DOUBLE
     g%lonc = .FALSE.
     ! ... dimensions
     g%lonm%ndims = 1
     ALLOCATE(g%lonm%dim(g%lonm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%lonm%dim(1)%name  = 'lon'
     g%lonm%dim(1)%id    = NULL_DIMID
     g%lonm%dim(1)%len   = ie_in
     g%lonm%dim(1)%fuid  = .false.

     ! ... data
     CALL INIT_NARRAY(g%lonm%dat, g%lonm%ndims &
          , (/g%lonm%dim(1)%len/),VTYPE_DOUBLE)
     DO i=1, ie_in
        itot = isubpos_in(my_cart_id,1)-nboundlines_in-1
        itmp = istartlon + (itot+i-1) * idlon
        g%lonm%dat%vd(i) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%lonm, 'long_name', vs='longitude')
     CALL ADD_NCATT(g%lonm, 'units', vs='degrees_east')
     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%latm%name  = 'lat'
     g%latm%id    = NULL_VARID
     g%latm%xtype = NF90_DOUBLE
     ! ... dimensions
     g%latm%ndims = 1
     ALLOCATE(g%latm%dim(g%latm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%latm%dim(1)%name  = 'lat'
     g%latm%dim(1)%id    = NULL_DIMID
     g%latm%dim(1)%len   = je_in
     g%latm%dim(1)%fuid  = .false.
     g%latm%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%latm%dat, g%latm%ndims &
       , (/g%latm%dim(1)%len/),VTYPE_DOUBLE)
     DO j=1,je_in
        jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
        itmp = istartlat + (jtot+j-1) * idlat
        g%latm%dat%vd(j) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%latm, 'long_name', vs='mid latitude')
     CALL ADD_NCATT(g%latm, 'units', vs='degrees_north')

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     g%loni%name  = 'lon_I'
     g%loni%id    = NULL_VARID
     g%loni%xtype = NF90_DOUBLE
     ! ... dimensions
     g%loni%ndims = 1
     ALLOCATE(g%loni%dim(g%loni%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%loni%dim(1)%name  = 'lon_I'
     g%loni%dim(1)%id    = NULL_DIMID
     g%loni%dim(1)%len   = ie_in+1
     g%loni%dim(1)%fuid  = .false.
     g%loni%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%loni%dat, g%loni%ndims &
          , (/g%loni%dim(1)%len/)  ,VTYPE_DOUBLE)

     DO i=1, ie_in+1
       itot = isubpos_in(my_cart_id,1)-nboundlines_in-1
       itmp = istartlon +  (itot+i) * idlon - NINT(1.5_dp * dlon_in * RCF)
       g%loni%dat%vd(i) = REAL(itmp,dp)/ RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%loni, 'long_name', vs='interface longitude')
     CALL ADD_NCATT(g%loni, 'units', vs='degrees_east')

     ! LATITUDE (MID) ...
     g%lati%name  = 'lat_I'
     g%lati%id    = NULL_VARID
     g%lati%xtype = NF90_DOUBLE
     ! ... dimensions
     g%lati%ndims = 1
     ALLOCATE(g%lati%dim(g%lati%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%lati%dim(1)%name  = 'lat_I'
     g%lati%dim(1)%id    = NULL_DIMID
     g%lati%dim(1)%len   = je_in+1
     g%lati%dim(1)%fuid  = .false.
     g%lati%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%lati%dat, g%lati%ndims &
          , (/g%lati%dim(1)%len/) ,VTYPE_DOUBLE)
     DO j = 1 , je_in+1
        jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
        itmp = istartlat + (jtot+j) * idlat - NINT(1.5_dp * dlat_in * RCF_IN)
        g%lati%dat%vd(j) = REAL(itmp,dp) /RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%lati, 'long_name', vs='interface latitude')
     CALL ADD_NCATT(g%lati, 'units', vs='degrees_north')

     ! *************************************************************************
     ! HYBRID-A-COEFFICIENTS (INTERFACES) ...
     g%hyai%name  = 'hyai'
     g%hyai%id    = NULL_VARID
     g%hyai%xtype = NF90_FLOAT
     ! ... dimensions
     g%hyai%ndims = 1
     ALLOCATE(g%hyai%dim(g%hyai%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
     g%hyai%dim(1)%name  = 'ilev'
     g%hyai%dim(1)%id    = NULL_DIMID
     g%hyai%dim(1)%len   = ke_int +1
     g%hyai%dim(1)%fuid  = .false.
     g%hyai%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%hyai%dat, g%hyai%ndims &
          , (/g%hyai%dim(1)%len/),VTYPE_REAL)
     g%hyai%dat%vr = REAL(ak_in(par_kmin:ke_in+1),sp)

     ! ... attributes
     CALL ADD_NCATT(g%hyai, 'long_name'                  &
          ,vs='hybrid-A-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hyai, 'units', vs='1')
     g%ranges(3,1) = RGEMPTY
     g%ranges(3,2) = RGEMPTY

     ! ************************************************************************
     ! HYBRID-B-COEFFICIENTS (INTERFACES) ...
     g%hybi%name  = 'hybi'
     g%hybi%id    = NULL_VARID
     g%hybi%xtype = NF90_FLOAT
     ! ... dimensions
     g%hybi%ndims = 1
     ALLOCATE(g%hybi%dim(g%hybi%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
     CALL COPY_NCDIM(g%hybi%dim(1),g%hyai%dim(1))
     ! ... data
     CALL INIT_NARRAY(g%hybi%dat, g%hybi%ndims &
          , (/g%hybi%dim(1)%len/),VTYPE_REAL)

     g%hybi%dat%vr = REAL(bk_in(par_kmin:ke_in+1),sp)
     ! ... attributes
     CALL ADD_NCATT(g%hybi, 'long_name'                  &
          ,vs='hybrid-B-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hybi, 'units', vs='1')
     g%ranges(4,1) = 0.0_sp
     g%ranges(4,2) = 1.0_sp

     ! ************************************************************************
     ! SURFACE PRESSURE
     g%ps%name  = 'ps'
     g%ps%id    = NULL_VARID
     g%ps%xtype = NF90_FLOAT
     ! ... dimensions
     g%ps%ndims = 2
     ALLOCATE(g%ps%dim(g%ps%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

     CALL COPY_NCDIM(g%ps%dim(1), g%lonm%dim(1))
     CALL COPY_NCDIM(g%ps%dim(2), g%latm%dim(1))
     ! ... data
     CALL INIT_NARRAY(g%ps%dat, g%ps%ndims      &
          ,(/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
          ,VTYPE_REAL)
        DO i=1, g%ps%dim(1)%len
           DO j=1, g%ps%dim(2)%len
              g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
                   ,(/i,j/)))  = 101325.0_sp
           END DO
        END DO
     ! ... attributes
     CALL ADD_NCATT(g%ps, 'long_name', vs='surface pressure')
     CALL ADD_NCATT(g%ps, 'units', vs='Pa')

     ! *************************************************************************
     ! REFERENCE PRESSURE ...
     g%p0%name  = 'p0'
     g%p0%id    = NULL_VARID
     g%p0%xtype = NF90_FLOAT
     ! ... dimensions
     g%p0%ndims = 0
     ! ... data
     CALL INIT_NARRAY(g%p0%dat, 1, (/ 1 /), VTYPE_REAL)
     g%p0%dat%vr(1) = 1.0_sp
     ! ... attributes
     CALL ADD_NCATT(g%p0, 'long_name', vs='reference pressure')
     CALL ADD_NCATT(g%p0, 'units', vs='Pa')

     ! TIME (MID) ...
     g%timem%name  = 'time'
     g%timem%id    = NULL_VARID
     g%timem%xtype = NF90_DOUBLE
     ! ... dimensions
     g%timem%ndims = 1
     ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
     g%timem%dim(1)%name  = 'time'
     g%timem%dim(1)%id    = NULL_DIMID
     g%timem%dim(1)%len   = 1
     g%timem%dim(1)%fuid  = .true.
     g%timem%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
          ,VTYPE_DOUBLE)
     ! TIME: SECONDS SINCE MODEL START
     CALL time_span_d(dts   &
          , YEAR_START, MONTH_START, DAY_START &
          , HOUR_START, MINUTE_START, SECOND_START  &
          , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
     g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds

     ! ... attributes
     CALL ADD_NCATT(g%timem, 'long_name'                              &
          ,vs='time in seconds since model start')
     WRITE(tunit, &
          '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
          YEAR_START, MONTH_START, DAY_START &
          , HOUR_START, MINUTE_START, SECOND_START

     CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

     ! CALCULATE INTs from MIDs
     CALL COMPLETE_GEOHYBGRID(g)

     ! INTERFACE OK
     ok = .true.

   END SUBROUTINE define_rotlocal_parent_grid

   ! ------------------------------------------------------------------
   SUBROUTINE DEFINE_ROTRED_BM_GRID(g, ok)

    !
    ! remove halo to avoid double counting of halo area

    ! MESSY/BMIL
     USE messy_main_mpi_bi,       ONLY: isubpos, my_cart_id, nboundlines
     USE messy_main_data_bi,      ONLY: aps
     USE messy_main_grid_def_mem_bi, ONLY: nlev => ke_tot             &
                                      , startlon_tot, startlat_tot &
                                      , dlon, dlat                 &
                                      , istart, iend, jstart, jend
     USE messy_main_grid_def_bi,   ONLY: hyai, hybi

     USE messy_main_timer,        ONLY: YEAR_START, MONTH_START, DAY_START     &
                                      , HOUR_START, MINUTE_START, SECOND_START &
                                      , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

     ! MESSy
     USE messy_main_constants_mem, ONLY: sp
     USE messy_main_timer,         ONLY: time_span_d

     ! NCREGRID MODULES
     USE messy_main_grid_trafo,  ONLY: RGEMPTY, COMPLETE_GEOHYBGRID
     USE messy_main_grid_netcdf, ONLY: ERRMSG, COPY_NCDIM       &
                                     , null_dimid, null_varid   &
                                     , vtype_real, vtype_double &
                                     , nf90_float, nf90_double  &
                                     , POSITION, INIT_NARRAY, ADD_NCATT
     USE messy_main_grid,        ONLY: INIT_GEOHYBGRID

     IMPLICIT NONE

     ! I/O
     TYPE(t_geohybgrid), INTENT(OUT) :: g
     LOGICAL           , INTENT(OUT) :: ok

     !LOCAL
     REAL(DP)             :: dts
     INTEGER              :: i,j, itot, jtot, itmp
     INTEGER              :: status
     CHARACTER(LEN=100)   :: tunit
     LOGICAL              :: l_ijf
     INTEGER              :: clonlen = 0
     INTEGER              :: clatlen = 0
     INTEGER              :: istartlon, istartlat, idlon, idlat


     istartlon = NINT(startlon_tot * RCF)
     istartlat = NINT(startlat_tot * RCF)
     idlon     = NINT(dlon * RCF)
     idlat     = NINT(dlat * RCF)

     ! DEFINITION OF START, END AND LENGTH OF REDUCED CHILD GRID
     ! start and end of inner region
     ! define locally to be able to exclude damping regions etc.

     IF(i_rmy_px < 0) THEN
        iis     = istart
        iie     = iend
        jjs     = jstart
        jje     = jend
     ELSE
        l_ijf = .FALSE.
        DO i=istart,iend
           DO j=jstart,jend
              IF(rmy_loc(i,j) <= 0._dp) THEN
                 iis=i
                 jjs=j
                 l_ijf = .TRUE.
                 exit
              END IF
           END DO
           IF(l_ijf) exit
        END DO
        IF(.NOT. l_ijf) THEN
           iis=0
           iie=-1
           jjs=0
           jje=-1
        ELSE
           l_ijf = .FALSE.
           DO i=iend,istart,-1
              DO j=jend,jstart,-1
                 IF(rmy_loc(i,j) <= 0._dp) THEN
                    iie=i
                    jje=j
                    l_ijf = .TRUE.
                    exit
                 END IF
              END DO
              IF(l_ijf) exit
           END DO
        END IF
        DEALLOCATE(rmy_loc)
     ENDIF
     clonlen = iie - iis + 1
     clatlen = jje - jjs + 1
     IF (clonlen == 0 .OR. clatlen == 0) THEN
        write(iouerr,*) 'CLON / CLAT no overlap ', clonlen , clatlen
        loverlap = .FALSE.
     ENDIF

     ! INIT
     CALL INIT_GEOHYBGRID(g)

     g%name = 'COSMOohalo'

     g%file = ' '   ! Filename
     g%t    = 0                               ! time step

     !**************************************************************************
     ! Rotated  Curvilinear LONGITUDE (MID) ...
     g%lonm%name  = 'lon'
     g%lonm%id    = NULL_VARID
     g%lonm%xtype = NF90_DOUBLE
     g%lonc = .FALSE.
     ! ... dimensions
     g%lonm%ndims = 1
     ALLOCATE(g%lonm%dim(g%lonm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%lonm%dim(1)%name  = 'lon'
     g%lonm%dim(1)%id    = NULL_DIMID
     g%lonm%dim(1)%len   = clonlen
     g%lonm%dim(1)%fuid  = .false.

     ! ... data
     CALL INIT_NARRAY(g%lonm%dat, g%lonm%ndims &
          , (/g%lonm%dim(1)%len/),VTYPE_DOUBLE)
     DO i=iis, iie
        itot = isubpos(my_cart_id,1)-nboundlines-1
        itmp = istartlon + (itot+i-1) * idlon
        g%lonm%dat%vd(i-iis+1) = REAL(itmp,dp) / RCF
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%lonm, 'long_name', vs='longitude')
     CALL ADD_NCATT(g%lonm, 'units', vs='degrees_east')
     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%latm%name  = 'lat'
     g%latm%id    = NULL_VARID
     g%latm%xtype = NF90_DOUBLE
     ! ... dimensions
     g%latm%ndims = 1
     ALLOCATE(g%latm%dim(g%latm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%latm%dim(1)%name  = 'lat'
     g%latm%dim(1)%id    = NULL_DIMID
     g%latm%dim(1)%len   = clatlen
     g%latm%dim(1)%fuid  = .false.
     g%latm%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%latm%dat, g%latm%ndims &
          , (/g%latm%dim(1)%len/), VTYPE_DOUBLE)
     DO j=jjs,jje
        jtot = isubpos(my_cart_id,2)-nboundlines-1
        itmp = istartlat + (jtot+j-1) * idlat
        g%latm%dat%vd(j-jjs+1) = REAL(itmp,dp) / RCF
    END DO

     ! ... attributes
     CALL ADD_NCATT(g%latm, 'long_name', vs='mid latitude')
     CALL ADD_NCATT(g%latm, 'units', vs='degrees_north')

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     g%loni%name  = 'lon_I'
     g%loni%id    = NULL_VARID
     g%loni%xtype = NF90_DOUBLE
     ! ... dimensions
     g%loni%ndims = 1
     ALLOCATE(g%loni%dim(g%loni%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%loni%dim(1)%name  = 'lon_I'
     g%loni%dim(1)%id    = NULL_DIMID
     g%loni%dim(1)%len   = clonlen+1
     g%loni%dim(1)%fuid  = .false.
     g%loni%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%loni%dat, g%loni%ndims &
          , (/g%loni%dim(1)%len/)  ,VTYPE_DOUBLE)

     DO i=iis, iie+1
        itot = isubpos(my_cart_id,1)-nboundlines-1
        itmp = istartlon + (itot+i) * idlon - NINT(1.5_dp * dlon * RCF)
        g%loni%dat%vd(i-iis+1) =  REAL(itmp,dp) / RCF
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%loni, 'long_name', vs='interface longitude')
     CALL ADD_NCATT(g%loni, 'units', vs='degrees_east')

     ! LATITUDE (MID) ...
     g%lati%name  = 'lat_I'
     g%lati%id    = NULL_VARID
     g%lati%xtype = NF90_DOUBLE
     ! ... dimensions
     g%lati%ndims = 1
     ALLOCATE(g%lati%dim(g%lati%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%lati%dim(1)%name  = 'lat_I'
     g%lati%dim(1)%id    = NULL_DIMID
     g%lati%dim(1)%len   = clatlen+1
     g%lati%dim(1)%fuid  = .false.
     g%lati%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%lati%dat, g%lati%ndims &
          , (/g%lati%dim(1)%len/) ,VTYPE_DOUBLE)
     DO j = jjs , jje+1
        jtot = isubpos(my_cart_id,2)-nboundlines-1
        itmp = istartlat + (jtot+j) * idlat - NINT(1.5_dp*dlat * RCF)
        g%lati%dat%vd(j-jjs+1) = REAL(itmp,dp)/RCF
      END DO

     ! ... attributes
     CALL ADD_NCATT(g%lati, 'long_name', vs='interface latitude')
     CALL ADD_NCATT(g%lati, 'units', vs='degrees_north')

     ! *************************************************************************
     ! HYBRID-A-COEFFICIENTS (INTERFACES) ...
     g%hyai%name  = 'hyai'
     g%hyai%id    = NULL_VARID
     g%hyai%xtype = NF90_FLOAT
     ! ... dimensions
     g%hyai%ndims = 1
     ALLOCATE(g%hyai%dim(g%hyai%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
     g%hyai%dim(1)%name  = 'ilev'
     g%hyai%dim(1)%id    = NULL_DIMID
     g%hyai%dim(1)%len   = nlev+1-chi_kmin+1
     g%hyai%dim(1)%fuid  = .false.
     g%hyai%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%hyai%dat, g%hyai%ndims, (/g%hyai%dim(1)%len/) &
          ,VTYPE_REAL)
     g%hyai%dat%vr(:) = REAL(hyai(chi_kmin:nlev+1),sp)

     ! ... attributes
     CALL ADD_NCATT(g%hyai, 'long_name'                              &
          ,vs='hybrid-A-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hyai, 'units', vs='1')
     g%ranges(3,1) = RGEMPTY
     g%ranges(3,2) = RGEMPTY

     ! *************************************************************************
     ! HYBRID-B-COEFFICIENTS (INTERFACES) ...
     g%hybi%name  = 'hybi'
     g%hybi%id    = NULL_VARID
     g%hybi%xtype = NF90_FLOAT
     ! ... dimensions
     g%hybi%ndims = 1
     ALLOCATE(g%hybi%dim(g%hybi%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
     CALL COPY_NCDIM(g%hybi%dim(1),g%hyai%dim(1))
     ! ... data
     CALL INIT_NARRAY(g%hybi%dat, g%hybi%ndims, (/g%hybi%dim(1)%len/) &
          ,VTYPE_REAL)
     g%hybi%dat%vr(:) = REAL(hybi(chi_kmin:nlev+1),sp)

     ! ... attributes
     CALL ADD_NCATT(g%hybi, 'long_name'                              &
          ,vs='hybrid-B-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hybi, 'units', vs='1')
     g%ranges(4,1) = 0.0_sp
     g%ranges(4,2) = 1.0_sp

     ! *************************************************************************
     ! SURFACE PRESSURE
     g%ps%name  = 'ps'
     g%ps%id    = NULL_VARID
     g%ps%xtype = NF90_FLOAT
     ! ... dimensions
     g%ps%ndims = 2
     ALLOCATE(g%ps%dim(g%ps%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

     g%ps%dim(1) = g%lonm%dim(1)
     g%ps%dim(2) = g%latm%dim(1)

     ! ... data
     CALL INIT_NARRAY(g%ps%dat, g%ps%ndims      &
          ,(/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
          ,VTYPE_REAL)
     IF (ASSOCIATED(aps)) THEN
        DO i=1, g%ps%dim(1)%len
           DO j=1, g%ps%dim(2)%len
              g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                   ,(/i,j/)))                                            &
                   !           = aps(i,j)
                   = 101325.0_sp
           END DO
        END DO
     ELSE
        DO i=1, g%ps%dim(1)%len
           DO j=1, g%ps%dim(2)%len
              g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                   ,(/i,j/)))                                            &
                   = 101325.0_sp
           END DO
        END DO
     END IF
     ! ... attributes
     CALL ADD_NCATT(g%ps, 'long_name', vs='surface pressure')
     CALL ADD_NCATT(g%ps, 'units', vs='Pa')

     ! *************************************************************************
     ! REFERENCE PRESSURE ...
     g%p0%name  = 'p0'
     g%p0%id    = NULL_VARID
     g%p0%xtype = NF90_FLOAT
     ! ... dimensions
     g%p0%ndims = 0
     ! ... data
     CALL INIT_NARRAY(g%p0%dat, 1, (/ 1 /), VTYPE_REAL)
     g%p0%dat%vr(1) = 1.0_sp
     ! ... attributes
     CALL ADD_NCATT(g%p0, 'long_name', vs='reference pressure')
     CALL ADD_NCATT(g%p0, 'units', vs='Pa')

     ! TIME (MID) ...
     g%timem%name  = 'time'
     g%timem%id    = NULL_VARID
     g%timem%xtype = NF90_DOUBLE
     ! ... dimensions
     g%timem%ndims = 1
     ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
     g%timem%dim(1)%name  = 'time'
     g%timem%dim(1)%id    = NULL_DIMID
     g%timem%dim(1)%len   = 1
     g%timem%dim(1)%fuid  = .true.
     g%timem%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
          ,VTYPE_DOUBLE)
     ! TIME: SECONDS SINCE MODEL START
     CALL time_span_d(dts                           &
          , YEAR_START, MONTH_START, DAY_START      &
          , HOUR_START, MINUTE_START, SECOND_START  &
          , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
     g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds

     ! ... attributes
     CALL ADD_NCATT(g%timem, 'long_name'         &
          ,vs='time in seconds since model start')
     WRITE(tunit, &
          '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
          YEAR_START, MONTH_START, DAY_START &
          , HOUR_START, MINUTE_START, SECOND_START

     CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

     ! CALCULATE INTs from MIDs
     CALL COMPLETE_GEOHYBGRID(g)

     ! INTERFACE OK
     ok = .true.

   END SUBROUTINE DEFINE_ROTRED_BM_GRID

 END SUBROUTINE CALC_BACKWARD_WEIGHTS
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE Setup_Pardata_exchange

    ! MESSy/BMIL
    USE messy_main_mpi_bi,        ONLY: num_compute, icomm_compute &
                                      , imp_reals, gather_values
    ! MMD
    USE mmd_child,                ONLY: MMD_Send_to_Parent         &
                                      , MMD_Recv_from_Parent       &
                                      , MMD_C_Set_ParIndexlist

    ! INT2COSMO
    USE data_grid_in,             ONLY: ie_in, je_in                     &
                                      , ie_in_max, je_in_max

    IMPLICIT NONE

    INTRINSIC :: INT, SIZE

    ! Local Data
    CHARACTER(LEN=*), PARAMETER :: substr='Setup_Pardata_exchange'
    INTEGER                     :: i,j
    INTEGER(I4)                 :: stat_I4
    INTEGER                     :: status
    INTEGER,DIMENSION(3)        :: size_of_array
    CHARACTER(len=80)           :: yerrmsg

    INTEGER                               :: icount(1)
    INTEGER,  DIMENSION(:,:), ALLOCATABLE :: IndexList
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: my_plon, my_plat
    REAL(dp), DIMENSION(:,:), POINTER     :: tmpfracs => NULL()
    INTEGER                               :: itmp

    REAL(kind=dp),DIMENSION(ie_in_max, je_in_max,0:num_compute-1) :: &
         all_lon, all_lat

    ! NO PARENT FIELDS
    IF (PD_NUM == 0) RETURN

    ALLOCATE(my_plon(ie_in_max, je_in_max))
    ALLOCATE(my_plat(ie_in_max, je_in_max))

    my_plon  = -9999.0_dp
    my_plat  = -9999.0_dp
    all_lon = -9999.0_dp
    all_lat = -9999.0_dp

    ! Compute coordinates of mid point of all local cells of the coarse grid
    DO i = 1, ie_in
       DO j = 1, je_in
          my_plon(i,j) = rlons(i) *1000._dp
          itmp         = NINT(my_plon(i,j))
          my_plon(i,j) = REAL(itmp,dp)/1000._dp
          my_plat(i,j) = rlats(j)*1000._dp
          itmp         = NINT(my_plat(i,j))
          my_plat(i,j) = REAL(itmp,dp)/1000._dp
          ! calculate geographical longitude and latitude
          IF (my_plat(i,j) < -400._dp  .or. my_plon(i,j) < -400._dp ) &
             CALL error_bi(substr, 'longitudes or latitudes undefined')
      ENDDO
    ENDDO

    CALL gather_values (my_plon, all_lon,INT(ie_in_max,I4),INT(je_in_max,I4) &
         ,num_compute, imp_reals, 0_I4, icomm_compute, yerrmsg, stat_I4)
    IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)

    CALL gather_values (my_plat, all_lat,INT(ie_in_max,I4),INT(je_in_max,I4) &
         , num_compute, imp_reals, 0_I4, icomm_compute, yerrmsg, stat_I4)
    IF (stat_I4 /= 0) CALL error_bi(yerrmsg, substr)

    DEALLOCATE (my_plon, my_plat)

    ! send coordinates to parent
    IF (p_parallel_io) THEN
       size_of_array(1) = size(all_lon,1)
       size_of_array(2) = size(all_lon,2)
       size_of_array(3) = size(all_lon,3)
       CALL MMD_Send_to_Parent (size_of_array, 3, 0, 140, status)
       IF (status /= 0) CALL error_bi('MMD_Send_to_Parent ERROR 01', substr)
       CALL MMD_Send_to_Parent (all_lon, size(all_lon), 0, 141, status)
       IF (status /= 0) CALL error_bi('MMD_Send_to_Parent ERROR 02', substr)
       CALL MMD_Send_to_Parent (all_lat, size(all_lat), 0, 142, status)
       IF (status /= 0) CALL error_bi('MMD_Send_to_Parent ERROR 03', substr)

       ! receive Index_list
       CALL MMD_Recv_from_Parent (icount, 1, 0, 144, status)
       IF (status /= 0)  CALL error_bi('MMD_Recv_from_Parent', substr)

       ALLOCATE (IndexList(6,icount(1)))
       CALL MMD_Recv_from_Parent(IndexList, SIZE(IndexList),0, 145, status)

    ELSE
       ALLOCATE (IndexList(6,1))
    ENDIF

    CALL determine_fractions(tmpfracs)

    ALLOCATE(fracs(ie_in, je_in))
    fracs = 0._dp

    DO i = 1, ie_in
       DO j =1, je_in
          IF (fw_int(i,j,1,1) > 0._dp) THEN
             fracs(i,j) = tmpfracs(i,j)
          ENDIF
       END DO
    END DO
    CALL MMD_C_Set_ParIndexlist(IndexList, fracs, fw_int(:,:,1,1))

    DEALLOCATE(IndexList)

    RETURN

 CONTAINS

   ! -------------------------------------------------------------
   SUBROUTINE determine_fractions(fractions)

     ! MESSy/SMCL
     USE messy_main_grid,            ONLY: grid_error
     USE messy_main_grid_tools,      ONLY: dealine_array

     ! INT2COSMO
     USE data_grid_in,               ONLY: ie_in, je_in

     IMPLICIT NONE

     ! I/O
     REAL(dp), DIMENSION(:,:), POINTER :: fractions
     ! LOCAL
     INTEGER                           :: ii
     INTEGER                           :: status

     IF (ASSOCIATED(fractions)) DEALLOCATE(fractions)
     ALLOCATE(fractions(ie_in,je_in))
     DO ii=1,PD_NUM
        IF (ParData(ii)%interpM == 1) THEN
           conserv_idx = ii
           EXIT
        ENDIF
     END DO
     IF (conserv_idx== 0 .OR. .NOT. loverlap) THEN
        ! no conservative interpolation required
        fractions(:,:) = 1._dp
     ELSE
        CALL DEALINE_ARRAY(status, ie_in, je_in &
             , ParData(conserv_idx)%PSD%wghts%dstfrac, fractions)
        IF (status /= 0) CALL error_bi(grid_error(status), substr)

        IF (MAXVAL(fractions) > 1.000001_dp) THEN
           write(iouerr,*) 'WARNING FRACS ', &
             MAXVAL(fractions), MAXVAL(ParData(conserv_idx)%PSD%wghts%dstfrac)

           CALL error_bi(substr, 'FRACs > 1')
        ENDIF
     ENDIF

   END SUBROUTINE determine_fractions
   ! -------------------------------------------------------------

  END SUBROUTINE Setup_Pardata_exchange
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE exchange_parent_data

    USE mmd_child,       ONLY: MMD_C_FillBuffer

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'exchange_parent_data'

    CALL debug_bi(submodstr, 'START exchange_parent_data', '', substr)

    ! Send  Data to Parent
    CALL MMD_C_FillBuffer

    CALL debug_bi(submodstr, 'END exchange_parent_data', '', substr)

   END SUBROUTINE exchange_parent_data
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE interpol_parent_data

    ! INT2COSMO
    USE data_int2lm_parallel,   ONLY: isubpos_in, nboundlines_in
    USE data_grid_in,           ONLY: ie_in, je_in, ke_in, akh_in, bkh_in  &
                                    , dlat_in, startlat_in_tot
    USE meteo_utilities,        ONLY: moist_split

    USE mmd_child,              ONLY: MMD_C_GetParentType  &
                                    , MMD_ParentIsEcham

    ! MESSy/BMIL
    USE messy_main_mpi_bi,           ONLY: my_cart_id
    USE messy_main_data_bi,          ONLY: p0, pp, nnew
    USE messy_main_grid_def_bi,      ONLY: hsurf
    USE messy_main_grid_def_mem_bi,  ONLY: ie, je, ke

    ! MESSy/SMCL
    USE messy_main_grid_tools,       ONLY: RGTOOL_CONVERT &
                                         , RGTOOL_CONVERT_DAT2VAR
    USE messy_main_grid_trafo,       ONLY: RG_INT
    USE messy_main_grid_trafo_nrgd,  ONLY: REGRID_CONTROL
    USE messy_main_grid_trafo_scrp,  ONLY: SCRIP_CONTROL
    USE messy_main_grid_netcdf,      ONLY: t_ncvar, INIT_NCVAR, COPY_NCVAR    &
                                         , QDEF_NCVAR                         &
                                         , MAIN_GRID_SET_MESSAGEMODE          &
                                         , MSGMODE_E
    USE messy_main_grid,             ONLY: INIT_GEOHYBGRID, COPY_GEOHYBGRID &
                                         , GRID_ERROR
    USE messy_main_constants_mem,    ONLY: rd, rv, pi, g

    IMPLICIT NONE

    INTRINSIC :: TRIM

    TYPE(t_ncvar), DIMENSION(:),  POINTER :: vari  => NULL()
    TYPE(t_ncvar)                         :: varo
    TYPE(t_ncvar), DIMENSION(:),  POINTER :: sovar => NULL()
    TYPE(t_ncvar), DIMENSION(:),  POINTER :: ovar  => NULL()

    REAL(dp), DIMENSION(:,:,:,:), POINTER :: dat   => NULL()

    CHARACTER(LEN=*), PARAMETER           :: substr = 'interpol_parent_data'
    INTEGER                               :: ii, ix     ! loop indices
    INTEGER                               :: status ! error status
    TYPE(t_geohybgrid)                    :: intgrid
    TYPE(t_geohybgrid)                    :: intzogrid
    TYPE(t_geohybgrid)                    :: intpsgrid
    INTEGER, DIMENSION(1)                 :: RGT
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: ptr  => NULL()

    REAL(dp), PARAMETER                   :: o_m_rdv = 1._dp-Rd/Rv
    REAL(dp), PARAMETER                   :: Rdv  = Rd/Rv
    REAL(dp), PARAMETER                   :: Rvd_m_o = Rv/Rd - 1.0_dp
    REAL(dp)                              :: zaqi, zbqi
    REAL(dp)                              :: qv_ijk

    REAL(dp), TARGET                      :: grhum(ie,je,ke,1)
    REAL(dp), TARGET                      :: grhum_in(ie_in,je_in,ke_int,1)
    ! time level for usage of COSMO_ORI time-level containing prognostic fields
    ! NOTE: tracer are treated differently
    INTEGER                               :: tlev
    INTEGER                               :: kkmin, kkmax

    REAL(dp), TARGET                      :: geou(ie,je,ke,1)
    REAL(dp), TARGET                      :: geov(ie,je,ke,1)
    REAL(dp), TARGET                      :: geops  (ie,je,1,1)

    REAL(dp), TARGET                      :: fic_cosmo(ie,je,1,1)
    REAL(dp), TARGET                      :: fis_cosmo(ie,je,1,1)
    REAL(dp), TARGET                      :: pres_cosmo(ie,je,ke,1)

    REAL(dp)                              :: ps_iter(ie_in,je_in,1,1)
    REAL(dp)                              :: maxps_iterdiff

    ! horizontal grid with ke vertical layers
    TYPE(t_geohybgrid)                    :: phgrid_3d
    LOGICAL                               :: lcalc (ie_in,je_in)
    INTEGER                               :: jtot
    REAL(dp)                              :: clat

    REAL(dp)                              :: hhl_hint(ie_in, je_in,ke+1)
    LOGICAL                               :: luv_calc = .FALSE.
    LOGICAL                               :: lvertint = .FALSE.

    CALL debug_bi(submodstr, 'interpol_parent_data', '', substr)

    ! NO PARENT FIELDS
    IF (.NOT. loverlap) RETURN
    IF (PD_NUM == 0)    RETURN

    IF (.NOT. l_shortcut) THEN
       tlev = nnew
    ELSE
       tlev = 1
    ENDIF

    CALL MAIN_GRID_SET_MESSAGEMODE(MSGMODE_E)
!    CALL MAIN_GRID_SET_MESSAGEMODE()

    IF (MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN
       IF (idx_p_u /= -99 .AND. idx_p_v /= -99) THEN
          ! calculate wind at grid point mids and rotate to
          ! geographical (ECHAM) coordinate system.
          ! parent coordinate system
          CALL CALC_WIND
          luv_calc= .TRUE.
       END IF
    ELSE
       ! moved from above, only required for COSMO parent
       IF (.NOT. ASSOCIATED(hsurf_full)) CALL calc_p_hsurf_full
    ENDIF

    IF (idx_p_ps /= -99) THEN
       geops(:,:,1,1) = ParData(idx_p_ps)%ptr_ori(:,:,tlev,1)
    END IF

    DO ii= 1, PD_NUM+PD_ADD
       IF (ASSOCIATED(ptr)) THEN
          write(iouerr,*) 'WARNING PTR ASSOCIATED !!'
          DEALLOCATE(ptr)
       ENDIF
       NULLIFY(ptr)
       ! TEST array needs not to be interpolated
       IF (TRIM(ParData(ii)%name%obj) == 'Test_Ar') CYCLE

       IF (ParData(ii)%rank /= 2) THEN
          kkmin = chi_kmin
          kkmax = SIZE(ParData(ii)%ptr_ori,3)
       ELSE
          kkmin = 1
          kkmax = 1
       END IF

       ! u and v are only re-adjusted to grid-midpoints for EMAC as parent
       IF (ii == idx_p_u .AND. luv_calc) THEN
          ! calculcate wind components on mass grid point
          ! U
          PTR => geou(iis:iie, jjs:jje,kkmin:kkmax,:)
       ELSE IF (ii == idx_p_v .AND. luv_calc) THEN
          ! calculcate wind components on mass grid point
          ! V
          PTR => geov(iis:iie, jjs:jje,kkmin:kkmax,:)
       ELSE IF (ii == idx_p_ps) THEN
          PTR => geops(iis:iie, jjs:jje,:,:)
       ELSE IF (ii==idx_p_fis) THEN
          ! FIS
          IF (l_shortcut) THEN
             PTR => p_fis(iis:iie, jjs:jje,:,:)
          ELSE
             fis_cosmo(:,:,1,1) = hsurf * g
             PTR => fis_cosmo(iis:iie, jjs:jje,:,:)
          END IF
       ELSE IF (ii==idx_p_fic) THEN
          ! FIC
          CALL CALC_FIC
          PTR => fic_cosmo(iis:iie, jjs:jje,:,:)
       ELSE IF (ii==idx_p_pres) THEN
          pres_cosmo(:,:,:,1) = p0(:,:,:) +  pp(:,:,:,tlev)
          PTR => pres_cosmo(iis:iie, jjs:jje,kkmin:kkmax,:)
       ELSE
          SELECT CASE(ParData(ii)%itimedep)
          CASE(0)
             ptr => ParData(ii)%ptr_ori(iis:iie, jjs:jje,kkmin:kkmax,:)
          CASE(2)
             ! time dependent horizontal field
             ptr => ParData(ii)%ptr_ori(iis:iie, jjs:jje,tlev:tlev,:)
          CASE (3)
             ! 3d field
             ptr => ParData(ii)%ptr_ori(iis:iie, jjs:jje,kkmin:kkmax,tlev:tlev)
          CASE DEFAULT
             CALL error_bi('UNKOWN itimedep / rank of timedependent field' &
                  , substr)
          END SELECT

       END IF

       IF (ASSOCIATED(ptr)) THEN

          ALLOCATE(vari(1))

          IF (ParData(ii)%rank == 2) THEN
             CALL RGTOOL_CONVERT_DAT2VAR(vari(1) &
                  , ptr, ParData(ii)%name%obj, chgrid, 'xyzn')
          ELSE
             CALL RGTOOL_CONVERT_DAT2VAR(vari(1) &
                  , ptr, ParData(ii)%name%obj, cgrid, 'xyzn')
          ENDIF
          NULLIFY(ptr)

          ! INTERPOLATE
          ! A) HORIZONTAL INTERPOLATION
          RGT(1) = ParData(ii)%interpM

          CALL SCRIP_CONTROL(status, ParData(ii)%SD_ID, cgrid, pgrid &
               , RGT, .TRUE., vari, sovar, intgrid, llrgz=.FALSE.)
          IF (status /= 0) CALL error_bi(grid_error(status), substr)

          ! CONVERT BACK

          IF (ParData(ii)%rank == 2) THEN
             CALL RGTOOL_CONVERT(sovar(1), dat, phgrid, 'xyzn')
          ELSE
             CALL RGTOOL_CONVERT(sovar(1), dat, intgrid, 'xyzn')
             CALL COPY_GEOHYBGRID(phgrid_3d, intgrid)
          ENDIF
          CALL INIT_GEOHYBGRID(intgrid)

          ParData(ii)%ptr_hint(:,:,1:SIZE(dat,4),1) = dat(:,:,1,:)

          DEALLOCATE(dat)
          NULLIFY(dat)

          DO ix= 1, SIZE(vari)
             CALL INIT_NCVAR(vari(ix))
             CALL INIT_NCVAR(sovar(ix))
          END DO
          DEALLOCATE(vari)
          NULLIFY(vari)
          DEALLOCATE(sovar)
          NULLIFY(sovar)

       END IF ! asso ptr

    END DO

    lvertint = .FALSE.
    DO ii= 1, PD_NUM
       IF (ParData(ii)%rank /= 2) THEN
          lvertint =.TRUE.
          EXIT
       END IF
    END DO

    IF (.NOT. lvertint) RETURN

    IF (MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN

       CALL vert_interpol_lm2echam

    ELSE

       CALL vert_interpol_lm2lm

    END IF

  CONTAINS

  SUBROUTINE vert_interpol_lm2echam

    USE messy_main_grid,       ONLY: INIT_GEOHYBGRID             &
                                   , COPY_GEOHYBGRID             &
                                   , GEOHYBGRID_UPDATE_SURFPRESS
    USE messy_main_grid_tools, ONLY: RGTOOL_CONVERT_DAT2PREDEFVAR
    IMPLICIT NONE

    INTEGER            :: ix, jx, iz, ii
    TYPE(t_ncvar)      :: psurf
    TYPE(t_geohybgrid) :: basegrid
    TYPE(t_geohybgrid) :: igrid

    IF (lvertsnrgd) THEN
       CALL COPY_GEOHYBGRID(igrid,cgrid)
       CALL DEFINE_PRESSURE_COMPONENTS(igrid, iis, iie, jjs,jje, kkmin, kkmax)
    END IF

    DO jx = 1, je_in
       DO ix = 1, ie_in
          IF (cofw_int(ix,jx,1,1) > 0._dp) THEN
             lcalc(ix,jx) = .TRUE.
          ELSE
             lcalc(ix,jx) = .FALSE.
          END IF
       END DO
    END DO
    IF (.NOT. ldiagonly) THEN

    CALL CALC_ALPS_1

    CALL COPY_NCVAR(psurf,pgrid%ps)
    CALL RGTOOL_CONVERT_DAT2PREDEFVAR(status, psurf &
            , ParData(idx_p_ps)%ptr_int, 'xy--', 'xyzn')
    CALL GEOHYBGRID_UPDATE_SURFPRESS(status,pgrid_id,psurf)

    DO iz=1, nitermax

       ps_iter = ParData(idx_p_ps)%ptr_int
       IF (lvertsnrgd) THEN
          CALL vert_c2p_ncinterpol(                                  &
               ParData(idx_p_t)%ptr_hint,  ParData(idx_p_t)%ptr_int  &
               , ParData(idx_p_t)%SD_ID,igrid)

          CALL vert_c2p_ncinterpol(                                  &
               ParData(idx_p_qv)%ptr_hint, ParData(idx_p_qv)%ptr_int &
               , ParData(idx_p_qv)%SD_ID, igrid)

          CALL vert_c2p_ncinterpol(                                  &
               ParData(idx_p_qc)%ptr_hint, ParData(idx_p_qc)%ptr_int &
               , ParData(idx_p_qc)%SD_ID, igrid)

          IF (idx_p_qi /= -99)&
               CALL vert_c2p_ncinterpol(                             &
               ParData(idx_p_qi)%ptr_hint, ParData(idx_p_qi)%ptr_int &
               , ParData(idx_p_qi)%SD_ID, igrid)

       ELSE
          CALL vert_c2p_ncinterpol_old(idx_p_t,                        &
               ParData(idx_p_t)%ptr_hint,  ParData(idx_p_t)%ptr_int    &
               , ParData(idx_p_ps)%ptr_hint, ParData(idx_p_ps)%ptr_int &
               , phgrid_3d)

          CALL vert_c2p_ncinterpol_old(idx_p_qv ,                      &
               ParData(idx_p_qv)%ptr_hint, ParData(idx_p_qv)%ptr_int   &
               , ParData(idx_p_ps)%ptr_hint, ParData(idx_p_ps)%ptr_int &
               , phgrid_3d)

          CALL vert_c2p_ncinterpol_old(idx_p_qc,                       &
               ParData(idx_p_qc)%ptr_hint, ParData(idx_p_qc)%ptr_int   &
               , ParData(idx_p_ps)%ptr_hint, ParData(idx_p_ps)%ptr_int &
               , phgrid_3d)

          IF (idx_p_qi /= -99)&
               CALL vert_c2p_ncinterpol_old(idx_p_qi,                  &
               ParData(idx_p_qi)%ptr_hint, ParData(idx_p_qi)%ptr_int   &
               , ParData(idx_p_ps)%ptr_hint, ParData(idx_p_ps)%ptr_int &
               , phgrid_3d)
       ENDIF

       CALL CALC_ALPS_2
       CALL RGTOOL_CONVERT_DAT2PREDEFVAR(status, psurf &
            , ParData(idx_p_ps)%ptr_int, 'xy--', 'xyzn')
       CALL GEOHYBGRID_UPDATE_SURFPRESS(status,pgrid_id,psurf)

       maxps_iterdiff = MAXVAL(ParData(idx_p_ps)%ptr_int-ps_iter)
       IF (maxps_iterdiff > 1.0_dp) THEN
          npsiter1    = iz + 1
          npsiter01   = iz + 1
          npsiter001  = iz + 1
          npsiter0001 = iz + 1
       ELSE IF (maxps_iterdiff > 0.1_dp) THEN
          npsiter01   = iz + 1
          npsiter001  = iz + 1
          npsiter0001 = iz + 1
       ELSE IF (maxps_iterdiff > 0.01_dp) THEN
          npsiter001  = iz + 1
          npsiter0001 = iz + 1
       ELSE IF (maxps_iterdiff > 0.001_dp) THEN
          npsiter0001 = iz + 1
       ENDIF

       iterps(iz)%ptr = ParData(idx_p_ps)%ptr_int(:,:,1,1)
     END DO
 END IF !  ldiagonly
    IF (lgrhum) CALL CALC_GRHUM(1)

    ! b) VERTICAL INTERPOLATION

    IF (itype_VI == 1) THEN

       DO ii = 1, PD_NUM+PD_ADD

          IF (TRIM(ParData(ii)%name%obj) == 'Test_Ar') CYCLE

          IF (ParData(ii)%rank == 3) THEN
             IF (ii == idx_p_t .OR. ii== idx_p_qc .OR. ii == idx_p_qv &
                  .OR. ii == idx_p_qi) CYCLE

             IF (lvertsnrgd) THEN
                CALL vert_c2p_ncinterpol(                          &
                   ParData(ii)%ptr_hint,       ParData(ii)%ptr_int &
                   , ParData(ii)%SD_ID, igrid)
             ELSE
                CALL vert_c2p_ncinterpol_old(ii, &
                     ParData(ii)%ptr_hint,       ParData(ii)%ptr_int &
                     , ParData(idx_p_ps)%ptr_hint, ParData(idx_p_ps)%ptr_int &
                     , phgrid_3d)

             ENDIF
          ELSE
             ! This overwrites ps and it should not be necessary for others
             ! because int points to hint
             IF (ii /= idx_p_ps) THEN
                ParData(ii)%ptr_int = ParData(ii)%ptr_hint
             END IF
          ENDIF

       END DO

       IF (lgrhum) CALL CALC_GRHUM(-1)

       CALL INIT_GEOHYBGRID(phgrid_3d)
       CALL INIT_GEOHYBGRID(igrid)
       CALL INIT_GEOHYBGRID(basegrid)
       CALL INIT_NCVAR(psurf)

    ELSE IF (itype_VI == 2) THEN
       CALL error_bi('itype_VI does currently not work', substr)
    ELSE
       CALL error_bi('Unknown vertical integration type', substr)
    END IF

    IF(MMD_C_GetParentType() == MMD_ParentIsEcham)   THEN
       jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
       IF (idx_p_u /= -99) THEN
          ! CONVERT to u * COS(lat)
          DO jx = 1, je_in
             clat = startlat_in_tot + REAL(jtot+jx-1,dp) * dlat_in
             ParData(idx_p_u)%ptr_int(:,jx,:,:) = &
                  ParData(idx_p_u)%ptr_int(:,jx,:,:) * ABS(COS(clat*pi/180))
          END DO
       END IF
       IF (idx_p_v /= -99) THEN
          ! CONVERT to u * COS(lat)
          DO jx = 1, je_in
             clat = startlat_in_tot + REAL(jtot+jx-1,dp) * dlat_in
             ParData(idx_p_v)%ptr_int(:,jx,:,:) = &
                  ParData(idx_p_v)%ptr_int(:,jx,:,:)* ABS(COS(clat*pi/180))
          END DO
       END IF
    END IF

    CALL debug_bi(submodstr, 'interpol parent data', '', substr)

  END SUBROUTINE vert_interpol_lm2echam

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------

  REAL(dp) FUNCTION sf_psat_w(x,y,z,v,w)

    REAL(dp), INTENT(IN) :: x, y ,z, v, w

    sf_psat_w = y * EXP(z*(x-v)/(x-w))

  END FUNCTION sf_psat_w

  !-----------------------------------------------------------------------

  REAL(dp) FUNCTION sf_psat_i(x,y,zi,v,wi)

    REAL(dp), INTENT(IN) :: x, y ,zi, v, wi

    sf_psat_i = y * EXP(zi*(x-v)/(x-wi))

  END FUNCTION sf_psat_i

  !-----------------------------------------------------------------------

  REAL(dp) FUNCTION sf_qsat (x,y,z,v)

    REAL(dp), INTENT(IN) :: x, y ,z, v

    sf_qsat = z * x / MAX( (y-v*x), 1.0_dp)

  END FUNCTION sf_qsat

  !-----------------------------------------------------------------------

  SUBROUTINE CALC_WIND

    USE messy_main_grid_def_mem_bi, ONLY: pollon, pollat
    USE messy_main_grid_def_bi,     ONLY:philon_2d, philat_2d
    USE messy_main_grid_trafo,      ONLY: uvrot2uv_vec, uv2uvrot_vec

    ! MMD
    USE mmd_child,                ONLY: MMD_C_GetParentType  &
                                      , MMD_ParentIsEcham
    ! INT2COSMO
    USE data_grid_in,             ONLY:  pollat_in, pollon_in

    IMPLICIT NONE

    INTEGER  :: ix,jx,kx

    ! calculcate wind components on mass grid point
    ! U
    geou = 0._dp

    DO kx = 1, ke
       DO jx = jjs, jje
          DO ix = MAX(2,iis), iie
             geou(ix,jx,kx,1) = (ParData(idx_p_u)%ptr_ori(ix-1,jx,kx,tlev) &
                  +  ParData(idx_p_u)%ptr_ori(ix,jx,kx,tlev)) *0.5_dp
          END DO
          geou(1,jx,kx,1) = geou(2,jx,kx,1)
       END DO
    END DO

    geov = 0._dp
    DO kx = 1, ke
       DO ix = iis, iie
          DO jx = MAX(2,jjs), jje
             geov(ix,jx,kx,1) = (ParData(idx_p_v)%ptr_ori(ix,jx-1,kx,tlev) &
                  +  ParData(idx_p_v)%ptr_ori(ix,jx,kx,tlev)) *0.5_dp
          END DO
          geov(ix,1,kx,1) = geov(ix,2,kx,1)
       END DO
    END DO

    IF ( MMD_C_GetParentType() == MMD_ParentIsEcham) THEN
       DO kx = 1, ke
          CALL uvrot2uv_vec(geou(iis:iie,jjs:jje,kx,1)                  &
               , geov(iis:iie,jjs:jje,kx,1)                             &
               , philat_2d(iis:iie,jjs:jje), philon_2d(iis:iie,jjs:jje) &
               , pollat, pollon, iie-iis+1, jje-jjs+1)
       END DO

    ELSE IF (L_gridrotParenteqChild) THEN
       ! both grids are rotated equally, nothing to do
    ELSE
       ! rotate to geographical system first
       DO kx = 1, ke
          CALL uvrot2uv_vec(geou(iis:iie,jjs:jje,kx,1)                  &
               , geov(iis:iie,jjs:jje,kx,1)                             &
               , philat_2d(iis:iie,jjs:jje), philon_2d(iis:iie,jjs:jje) &
               , pollat, pollon, iie-iis+1, jje-jjs+1)
       END DO
       ! rotate velocities from geographical system to parent grid
       DO kx = 1, ke
          CALL uv2uvrot_vec(geou(iis:iie,jjs:jje,kx,1)                  &
               , geov(iis:iie,jjs:jje,kx,1)                             &
               , philat_2d(iis:iie,jjs:jje), philon_2d(iis:iie,jjs:jje) &
               , pollat_in, pollon_in, iie-iis+1, jje-jjs+1)
       END DO

    ENDIF

  END SUBROUTINE calc_wind
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE calc_alps_1

    ! MESSy/BMIL
    USE messy_main_grid_def_bi,   ONLY: hyai, hybi, hyam, hybm
    ! MESSy/SMCL
    USE messy_main_constants_mem, ONLY: g, Rd

    IMPLICIT NONE

    ! LOCAL
    INTEGER  :: ix,jx,kx
    REAL(dp) :: dfis(ie_in,je_in), zfio(ie_in,je_in), zfiu(ie_in,je_in)
    REAL(dp) :: ztvpq(ie_in,je_in), ztvq(ie_in,je_in)
    REAL(dp) :: ztv(ie_in,je_in,ke_hint)
    REAL(dp) :: zbtv(ie_in,je_in),zctv(ie_in,je_in),zdtv(ie_in,je_in)
    REAL(dp) :: zpq(ie_in,je_in), zpq2(ie_in,je_in)
    REAL(dp) :: zpn(ie_in,je_in,ke_hint+1)
    REAL(dp) :: zpdfis, zpht, zpresm
    INTEGER  :: kgr, kdiff
    REAL(dp), PARAMETER :: zpgr = 850.0E2_dp

    ! A] make hydrastatic adjustment of surface presse (TODO)

    ! difference of surface geopotentials (input-output = child - parent)
    IF (MMD_C_GetParentType() == MMD_ParentIsEcham) THEN
       dfis(:,:)  =  ParData(idx_p_fis)%ptr_hint(:,:,1,1) &
            - CplData(idx_cl_fis)%ptr_in(:,:,1,1)
    ELSE
       dfis(:,:)  =  ParData(idx_p_fis)%ptr_hint(:,:,1,1) &
            - CplData(idx_cl_fis)%ptr_in(:,:,1,1) * g
    ENDIF

    zpn(:,:,ke_hint+1)   = ParData(idx_p_ps)%ptr_hint(:,:,1,1)

    ! use the pressure on COSMO half levels (pressi)
    ! and calculate the virtual temperature (tv) on COSMO main levels
    ztv = 0._dp

    kdiff = ke - ke_hint
    DO kx = 1, ke_hint
       DO jx = 1, je_in
          DO ix = 1, ie_in
             IF (.NOT. lcalc(ix,jx)) CYCLE
             zpn(ix,jx,kx) = hyai(kx+kdiff) +                        &
                  hybi(kx+kdiff)*ParData(idx_p_ps)%ptr_hint(ix,jx,1,1)
             IF (idx_p_qi /= -99) THEN
                ztv(ix,jx,kx) = ParData(idx_p_t)%ptr_hint(ix,jx,kx,1)*(1. + &
                     Rvd_m_o*ParData(idx_p_qv)%ptr_hint(ix,jx,kx,1)         &
                     - ParData(idx_p_qc)%ptr_hint(ix,jx,kx,1)               &
                     - ParData(idx_p_qi)%ptr_hint(ix,jx,kx,1))
             ELSE
                ztv(ix,jx,kx) = ParData(idx_p_t)%ptr_hint(ix,jx,kx,1)*(1. + &
                     Rvd_m_o*ParData(idx_p_qv)%ptr_hint(ix,jx,kx,1)         &
                     - ParData(idx_p_qc)%ptr_hint(ix,jx,kx,1))
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    ! boundary layer and gradient of virtual temperature (zdtv)
    ! Find the boundary layer height (kgr)
    DO kx = ke, 1, -1
       zpht = hyam(kx) + hybm(kx)*1.0E5_dp
       IF (zpht > zpgr) kgr = kx
    ENDDO
    kgr = kgr - 1 - kdiff

    ! Regression of the 4 layers above k = kgr for linear approximation of
    ! the virtual temperature (ztv(p) = zbtv*p + zctv)
    zpq  (:,:)  = 0.0_dp
    zpq2 (:,:)  = 0.0_dp
    ztvq (:,:)  = 0.0_dp
    ztvpq(:,:)  = 0.0_dp

    DO kx = kgr-3, kgr
       DO jx = 1,je_in
          DO ix = 1,ie_in
             IF (.NOT. lcalc(ix,jx)) CYCLE
             zpresm = 0.5_dp * (zpn(ix,jx,kx+1) + zpn(ix,jx,kx))
             zpq(ix,jx)   = zpq  (ix,jx) + zpresm
             zpq2 (ix,jx) = zpq2 (ix,jx) + zpresm * zpresm
             ztvq (ix,jx) = ztvq (ix,jx) + ztv(ix,jx,kx)
             ztvpq(ix,jx) = ztvpq(ix,jx) + ztv(ix,jx,kx) * zpresm
          ENDDO
       ENDDO
    ENDDO

    ! The gradient of the virtual temperature (zdtv) is computed from zbtv
    ! This is used for calculating the difference of the temperature
    ! regarding the change of the geopotential difference fis_gl - fis_lm.
    DO jx = 1,je_in
       DO ix = 1,ie_in
          IF (.NOT. lcalc(ix,jx)) CYCLE
          zbtv(ix,jx) = (ztvpq(ix,jx) - ztvq(ix,jx) * zpq(ix,jx) / 4.0_dp) /&
               (zpq2(ix,jx) -  zpq(ix,jx) * zpq(ix,jx) / 4.0_dp)
          zctv(ix,jx) = (ztvq (ix,jx) - zbtv(ix,jx) * zpq(ix,jx)) / 4.0_dp

          zdtv(ix,jx) =  zbtv (ix,jx) * zpgr * dfis  (ix,jx) /              &
               (rd * (zbtv(ix,jx)*zpgr + zctv(ix,jx)))
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    ! Compute pressure adaptation (see int2lm src_vert_interpol)
    !-----------------------------------------------------------------------

    zfio (:,:)  = 0.0_dp
    zfiu (:,:)  = 0.0_dp
    zpq  (:,:)  = 0.0_dp
    ztvq (:,:)  = 0.0_dp
    ztvpq(:,:)  = 0.0_dp

    do_height: DO kx = ke_hint, 2, -1
       DO jx = 1,je_in
          DO ix = 1,ie_in
             IF (.NOT. lcalc(ix,jx)) CYCLE

             zfio(ix,jx) = zfiu(ix,jx) + rd*ztv(ix,jx,kx) &
                  *LOG(zpn(ix,jx,kx+1)/zpn(ix,jx,kx))

             IF     (dfis(ix,jx) < -10.0_dp) THEN
                ! parent orography is higher than child grid orography
                ! Look for the proper layer in the child model,
                ! where the parent model orography is in
                IF ( (ABS(dfis(ix,jx)) >  zfiu(ix,jx)) .AND.        &
                     (ABS(dfis(ix,jx)) <= zfio(ix,jx)) ) THEN
                   ParData(idx_p_ps)%ptr_int(ix,jx,1,1) = zpn(ix,jx,kx+1)*&
                        EXP( (ABS(dfis(ix,jx)) - zfiu(ix,jx))         &
                        / (zfio (ix,jx)  - zfiu(ix,jx))            &
                        * LOG(zpn(ix,jx,kx)/zpn(ix,jx,kx+1))  )
                ENDIF
             ELSEIF (ABS(dfis(ix,jx)) <= 10.0_dp) THEN

                ! orographies are about the same, so the pressure is also
                ! this holds for all k, but has only to be set once
                ParData(idx_p_ps)%ptr_int(ix,jx,1,1) = &
                     ParData(idx_p_ps)%ptr_hint(ix,jx,1,1)

             ELSEIF (dfis(ix,jx) > 10.0_dp) THEN

                ! Compute geops for dfis   > 0
                ! (child grid orography higher than parent orography)
                IF (dfis(ix,jx) > zfiu(ix,jx) .AND. &
                     dfis(ix,jx) <= zfio(ix,jx)) THEN
                   ! Compute the average temperature (tvq) in the layer
                   ! from fis_gl to fis_gl + zdfis   for zdfis   > 0
                   zpdfis = (dfis(ix,jx) - zfiu(ix,jx))/(zfio(ix,jx) &
                        - zfiu(ix,jx))* (zpn(ix,jx,kx)               &
                        - zpn(ix,jx,kx+1)) + zpn(ix,jx,kx+1)
                   ztvpq(ix,jx) = ztv(ix,jx,kx)*   &
                        (zpn(ix,jx,kx+1) - zpdfis)  + ztvpq(ix,jx)
                   zpq (ix,jx)  = zpn(ix,jx,kx+1) - zpdfis + zpq(ix,jx)
                   ztvq(ix,jx)  = ztvpq(ix,jx)/zpq(ix,jx)
                   ParData(idx_p_ps)%ptr_int(ix,jx,1,1) =      &
                        ParData(idx_p_ps)%ptr_hint(ix,jx,1,1)* &
                        EXP(dfis(ix,jx)/(rd*(ztvq(ix,jx) + zdtv(ix,jx))))
                ELSE
                   ztvpq(ix,jx) = ztv(ix,jx,kx)* &
                        (zpn(ix,jx,kx+1) - zpn(ix,jx,kx)) + ztvpq(ix,jx)
                   zpq(ix,jx) =  zpn(ix,jx,kx+1) - zpn(ix,jx,kx) + zpq(ix,jx)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       zfiu(:,:) = zfio(:,:)
    ENDDO do_height

  END SUBROUTINE calc_alps_1
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE calc_alps_2

    ! INT2COSMO
    USE data_grid_in,             ONLY: ak_in, bk_in, ie_in, je_in
    USE messy_main_constants_mem, ONLY: Rd, g

    IMPLICIT NONE

    INTEGER  :: ix,jx,kx
    INTEGER  :: iziter, kediff
    REAL(dp) :: zpsold(ie_in,je_in),zpsnew(ie_in,je_in)
    REAL(dp) :: ztv(ie_in,je_in,ke_hint)
    REAL(dp) :: zfio(ie_in,je_in), zfiu(ie_in,je_in), zpc_fi(ie_in,je_in)
    REAL(dp) :: zfic(ie_in,je_in), zdfidps(ie_in,je_in)
    REAL(dp) :: zpno, zpnu

    ! INITIALISE
    zpsold  = 0._dp
    zpsnew  = 0._dp
    ztv     = 0._dp
    zfio    = 0._dp
    zfiu    = 0._dp
    zpc_fi  = 0._dp
    zfic    = 0._dp
    zdfidps = 0._dp

    kediff = ke_in - ke_int
    ! Compute virtual temperature on main levels
    IF (idx_p_qi /= -99) THEN
       DO kx = 2, ke_int
          DO jx = 1, je_in
             DO ix = 1, ie_in
                IF (.NOT. lcalc(ix,jx)) CYCLE
                   ztv(ix,jx,kx) = ParData(idx_p_t)%ptr_int(ix,jx,kx,1)     &
                        *(1. + Rvd_m_o*ParData(idx_p_qv)%ptr_int(ix,jx,kx,1)&
                        - ParData(idx_p_qc)%ptr_int(ix,jx,kx,1)             &
                        - ParData(idx_p_qi)%ptr_int(ix,jx,kx,1))
             ENDDO
          ENDDO
       ENDDO
    ELSE
       DO kx = 2, ke_int
          DO jx = 1, je_in
             DO ix = 1, ie_in
                IF (.NOT. lcalc(ix,jx)) CYCLE
                   ztv(ix,jx,kx) = ParData(idx_p_t)%ptr_int(ix,jx,kx,1)     &
                        *(1. + Rvd_m_o*ParData(idx_p_qv)%ptr_int(ix,jx,kx,1)&
                        - ParData(idx_p_qc)%ptr_int(ix,jx,kx,1))
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    !-----------------------------------------------------------------------
    ! Section 2: Compute surface pressure with Newton iterations
    !-----------------------------------------------------------------------

    iziter      = 0

    zpsold(:,:) = ParData(idx_p_ps)%ptr_int(:,:,1,1)

    iteration: DO WHILE (iziter < 5)
       iziter = iziter + 1
       IF (MMD_C_GetParentType() == MMD_ParentIsEcham) THEN
           zfiu (:,:)   = CplData(idx_cl_fis)%ptr_in(:,:,1,1)
       ELSE
           zfiu (:,:)   = CplData(idx_cl_fis)%ptr_in(:,:,1,1) * g
       ENDIF

       zfic (:,:)   = 0.0_dp
       zdfidps(:,:) = 0.0_dp

       zpc_fi(:,:)  = pcontrol_fi

       DO kx = ke_int, 2, - 1
          DO jx = 1, je_in
             DO ix = 1, ie_in
                IF (.NOT. lcalc(ix,jx)) CYCLE
                zpno = ak_in(kx+kediff) + bk_in(kx+kediff) * zpsold(ix,jx)
                zpnu = ak_in(kx+1+kediff) &
                     + bk_in(kx+1+kediff) * zpsold(ix,jx)
                zfio(ix,jx) = zfiu(ix,jx)  + Rd * ztv(ix,jx,kx) * LOG(zpnu/zpno)
                zdfidps(ix,jx) = zdfidps(ix,jx) + Rd * ztv(ix,jx,kx) * &
                     (bk_in(kx+1+kediff)*zpno                          &
                     - bk_in(kx+kediff)*zpnu)/(zpnu*zpno)
                IF (zpnu > zpc_fi(ix,jx) .AND. zpno <= zpc_fi(ix,jx)) THEN
                   zfic(ix,jx) = zfiu(ix,jx) + Rd*ztv(ix,jx,kx)&
                        *LOG(zpnu/zpc_fi(ix,jx))
                   zdfidps(ix,jx) = zdfidps(ix,jx) &
                        + Rd*ztv(ix,jx,kx)*bk_in(kx+1+kediff)/zpnu
                ENDIF
             ENDDO
          ENDDO
          zfiu(:,:) = zfio(:,:)
       END DO
        DO ix = 1, ie_in
          Do jx = 1, je_in
             IF (zdfidps(ix,jx) /= 0._dp) THEN
                zpsnew(ix,jx) = zpsold(ix,jx) &
                     - (zfic(ix,jx)- ParData(idx_p_fic)%ptr_hint(ix,jx,1,1)) &
                     /zdfidps(ix,jx)
                zpsold(ix,jx) = zpsnew(ix,jx)
             ENDIF
          END DO
       END DO

    ENDDO iteration

    ParData(idx_p_ps)%ptr_int(:,:,1,1) = zpsnew(:,:)

  END SUBROUTINE calc_alps_2
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE calc_fic

    ! MESSy/BMIL
    USE messy_main_grid_def_bi,       ONLY: hyai, hybi
    ! MESSy/SMCL
    USE messy_main_constants_mem, ONLY: g, Rd

    IMPLICIT NONE

    ! LOCAL
    REAL(dp) :: ztv(ie,je)
    REAL(dp) :: zpo (ie,je), zpu(ie,je)
    REAL(dp) :: zfiu(ie,je), zfio(ie,je)
    INTEGER  :: kx

    zpu    (:,:)   = geops  (:,:,1,1)
    zfiu   (:,:)   = hsurf * g
    fic_cosmo(:,:,1,1) = 0.0_dp

    DO kx = ke, 1, - 1

       IF (kx == 1) THEN
          zpo (:,:) = 0.0
          IF (idx_p_qi /= -99) THEN
             ztv (:,:) = ParData(idx_p_t)%ptr_ori(:,:,kx,tlev)*(1.0_dp + &
                  Rvd_m_o*ParData(idx_p_qv)%ptr_ori(:,:,kx,1) &
                  - ParData(idx_p_qc)%ptr_ori(:,:,kx,1)       &
                  - ParData(idx_p_qi)%ptr_ori(:,:,kx,1))
          ELSE
             ztv (:,:) = ParData(idx_p_t)%ptr_ori(:,:,kx,tlev)*(1.0_dp + &
                  Rvd_m_o*ParData(idx_p_qv)%ptr_ori(:,:,kx,1) &
                  - ParData(idx_p_qc)%ptr_ori(:,:,kx,1))
          ENDIF

          zfio(:,:) = zfiu(:,:) + Rd*ztv(:,:)*LOG(2.0_dp)
       ELSE
          zpo (:,:) = hyai(kx)    + hybi(kx) * geops(:,:,1,1)
          IF (idx_p_qi /= -99) THEN
             ztv (:,:) = ParData(idx_p_t)%ptr_ori(:,:,kx,tlev)*(1.0_dp + &
                  Rvd_m_o * ParData(idx_p_qv)%ptr_ori(:,:,kx,1)          &
                  - ParData(idx_p_qc)%ptr_ori(:,:,kx,1)                  &
                  - ParData(idx_p_qi)%ptr_ori(:,:,kx,1))
          ELSE
             ztv (:,:) = ParData(idx_p_t)%ptr_ori(:,:,kx,tlev)*(1.0_dp + &
                  Rvd_m_o*ParData(idx_p_qv)%ptr_ori(:,:,kx,1)            &
                  - ParData(idx_p_qc)%ptr_ori(:,:,kx,1))
          ENDIF

          zfio(:,:) = zfiu(:,:) + Rd*ztv(:,:)*LOG(zpu(:,:)/zpo(:,:))
       ENDIF

       WHERE (zpu(:,:) > pcontrol_fi .AND. zpo(:,:) <= pcontrol_fi)
          fic_cosmo(:,:,1,1) = &
               zfiu(:,:) + Rd*ztv(:,:)*LOG(zpu(:,:)/pcontrol_fi)
       END WHERE

       zpu (:,:) = zpo (:,:)
       zfiu(:,:) = zfio(:,:)

    ENDDO

  END SUBROUTINE calc_fic
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE calc_grhum(flag)

    ! INT2COSMO
    USE data_int2lm_constants,    ONLY: B1, B2_w, B2_i, B3, B4_w, B4_i

    IMPLICIT NONE

    ! IN
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    REAL(dp) :: temp(ie_in,je_in,ke_hint)
    REAL(dp) :: press(ie_in,je_in,ke_hint)
    INTEGER  :: kediff, ix, jx, kx, i, j
    REAL(dp) :: zph(ie_in,je_in)
    REAL(dp) :: zaq, zbq
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: qi_dummy  => NULL()

    SELECT CASE(flag)

    CASE(1)
       ! CALCULATE GENERALISED HUMIDITY, if required
       CALL debug_bi(submodstr, 'START humidity coupling', '', substr)

       temp(:,:,:)  = ParData(idx_p_t)%ptr_hint(:,:,:,1)    ! new
       press(:,:,:) = ParData(idx_p_pres)%ptr_hint(:,:,:,1) ! new

       grhum = 0._dp

       kkmin = chi_kmin
       !
       kkmax = SIZE(ParData(idx_p_qv)%ptr_hint,3)
       !
       DO ix = 1, SIZE(ParData(idx_p_qv)%ptr_hint,1)
          DO jx = 1, SIZE(ParData(idx_p_qv)%ptr_hint,2)
             DO kx = kkmin, kkmax
                IF (.NOT. lcalc(ix,jx)) CYCLE
                ! calculate saturation pressure
                zaq  = sf_psat_w (temp(ix,jx,kx), B1, B2_w, B3, B4_w)
                zbq  = sf_qsat   (zaq, press(ix,jx,kx), rdv, O_m_rdv)
                zaqi = sf_psat_i (temp(ix,jx,kx), B1, B2_i, B3, B4_i)
                zbqi = sf_qsat   (zaqi,press(ix,jx,kx), rdv, O_m_rdv)

                qv_ijk = ParData(idx_p_qv)%ptr_hint(ix,jx,kx,1)
                IF (qv_ijk > zbq) THEN
                   grhum(ix,jx,kx,1) = &
                        ParData(idx_p_qv)%ptr_hint(ix,jx,kx,1) / zbq

                ELSE IF (qv_ijk < zbq .AND. qv_ijk > zbqi) THEN
                   grhum(ix,jx,kx,1)=(ParData(idx_p_qv)%ptr_hint(ix,jx,kx,1)&
                        +  ParData(idx_p_qc)%ptr_hint(ix,jx,kx,1) ) / zbq
                ELSE IF (idx_p_qi > 0 .AND. itype_VI /= 2) THEN
                   grhum(ix,jx,kx,1)=(ParData(idx_p_qv)%ptr_hint(ix,jx,kx,1)&
                        +  ParData(idx_p_qc)%ptr_hint(ix,jx,kx,1)           &
                        +  ParData(idx_p_qi)%ptr_hint(ix,jx,kx,1)) / zbq
                ELSE
                   grhum(ix,jx,kx,1)=(ParData(idx_p_qv)%ptr_hint(ix,jx,kx,1)&
                        +  ParData(idx_p_qc)%ptr_hint(ix,jx,kx,1) ) / zbq
                END IF
             END DO
          END DO
       END DO

       ParData(idx_p_qv)%ptr_hint = grhum

       CALL debug_bi(submodstr, 'END humidity coupling', '', substr)

    CASE(-1)
       ! REDISTRIBUTE generalised humidity
       CALL debug_bi(submodstr, 'START moist split', '', substr)
       grhum_in = ParData(idx_p_qv)%ptr_int
       ParData(idx_p_qv)%ptr_int = 0._dp
       ParData(idx_p_qc)%ptr_int = 0._dp

       kediff = ke_in-ke_int

       IF (idx_p_qi /= -99 .AND. itype_VI /= 2) THEN

          ParData(idx_p_qi)%ptr_int = 0._dp

          DO kx = 1, ke_int
             zph(:,:) = akh_in(kediff+kx) + bkh_in(kediff+kx)&
                  * ParData(idx_p_ps)%ptr_int(:,:,1,1)
             DO i=1,ie_in
                DO j= 1, je_in
                   IF (.NOT. lcalc(i,j)) CYCLE
                   CALL moist_split( &
                        ParData(idx_p_t)%ptr_int(i:i,j:j,kx,1)&
                        , zph(i:i,j:j), grhum_in(i:i,j:j,kx,1)         &
                        , 0._dp,0._dp,0._dp &
                        , pi,b1,b2_w,b2_i,b3,b4_w,b4_i,Rdv,O_m_rdv     &
                        , ParData(idx_p_qv)%ptr_int(i:i,j:j,kx,1)      &
                        , ParData(idx_p_qc)%ptr_int(i:i,j:j,kx,1)      &
                        , ParData(idx_p_qi)%ptr_int(i:i,j:j,kx,1), 1,1)
                END DO
             END DO
          END DO
       ELSE IF (itype_VI /= 2) THEN
          ALLOCATE(qi_dummy(ie_in,je_in,ke_int,1))
          qi_dummy = 0._dp

          DO kx = 1, ke_int
             zph(:,:) = akh_in(kediff+kx) + bkh_in(kediff+kx)&
                  * ParData(idx_p_ps)%ptr_int(:,:,1,1)
             DO i=1,ie_in
                DO j= 1, je_in
                   IF (.NOT. lcalc(i,j)) CYCLE
                   CALL moist_split( &
                        ParData(idx_p_t)%ptr_int(i:i,j:j,kx,1)&
                        , zph(i:i,j:j), grhum_in(i:i,j:j,kx,1)         &
                        , 0._dp,0._dp,0._dp &
                        , pi,b1,b2_w,b2_i,b3,b4_w,b4_i,Rdv,O_m_rdv     &
                        , ParData(idx_p_qv)%ptr_int(i:i,j:j,kx,1)      &
                        , ParData(idx_p_qc)%ptr_int(i:i,j:j,kx,1)      &
                        , qi_dummy(i:i,j:j,kx,1), 1,1)
                   IF (qi_dummy(i,j,kx,1)> 0._dp) THEN
                      write(iouerr,*) 'WARNING moist_split qi > 0. adding to qc' &
                           , i,j,kx
                      ParData(idx_p_qc)%ptr_int(i,j,kx,1) = &
                        ParData(idx_p_qc)%ptr_int(i,j,kx,1) + qi_dummy(i,j,kx,1)
                   END IF
                 END DO
             END DO
          END DO

          DEALLOCATE(qi_dummy)
          NULLIFY(qi_dummy)
       ELSE
          CALL error_bi('MesoTEl grhum version, should not yet be called','')
       ENDIF

       CALL debug_bi(submodstr, 'END calc grhum', '', substr)

    END SELECT

  END SUBROUTINE calc_grhum
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE vert_c2p_ncinterpol_old(idx, var_in, var_out, ps_s, ps_d &
                                    , grid_s,grid_d)

    USE messy_main_grid,          ONLY: GRID_ERROR
    USE messy_main_grid_tools,    ONLY: SET_SURFACE_PRESSURE
    USE messy_main_grid_trafo,    ONLY: COMPLETE_GEOHYBGRID

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: idx
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: var_in
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: var_out
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: ps_s
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: ps_d
    TYPE(t_geohybgrid), INTENT(IN)            :: grid_s
    TYPE(t_geohybgrid), INTENT(OUT), OPTIONAL :: grid_d
    INTEGER                                   :: ix

    ! SURFACE PRESSURE needs to be available for both grid
    ! on the rotated grid

    ! --------------------
    ! make input grid
    ! --------------------

    CALL COPY_GEOHYBGRID(intpsgrid, grid_s)
    ! -- correct also lon/lat definition in ps field

    IF ((QDEF_NCVAR(pgrid%rlonm) .AND. QDEF_NCVAR(pgrid%rlatm))  .OR. &
         (QDEF_NCVAR(pgrid%rloni) .AND. QDEF_NCVAR(pgrid%rlati))) THEN
       ! -- correct also lon/lat definition
       intpsgrid%lonc=.FALSE.
       CALL COPY_NCVAR(intpsgrid%lonm, pgrid%rlonm)
       CALL COPY_NCVAR(intpsgrid%loni, pgrid%rloni)
       CALL COPY_NCVAR(intpsgrid%latm, pgrid%rlatm)
       CALL COPY_NCVAR(intpsgrid%lati, pgrid%rlati)

       CALL INIT_NCVAR(intpsgrid%rlonm)
       CALL INIT_NCVAR(intpsgrid%rloni)
       CALL INIT_NCVAR(intpsgrid%rlatm)
       CALL INIT_NCVAR(intpsgrid%rlati)
    ELSE
       ! ??? required ??
       ! -- correct also lon/lat definition
       CALL COPY_NCVAR(intpsgrid%lonm, pgrid%lonm)
       CALL COPY_NCVAR(intpsgrid%loni, pgrid%loni)
       CALL COPY_NCVAR(intpsgrid%latm, pgrid%latm)
       CALL COPY_NCVAR(intpsgrid%lati, pgrid%lati)
    ENDIF
    ! remove curvilinear information
    ! vertical interpolation proceeds on target grid
    CALL INIT_NCVAR(intpsgrid%clonm)
    CALL INIT_NCVAR(intpsgrid%cloni)
    CALL INIT_NCVAR(intpsgrid%clatm)
    CALL INIT_NCVAR(intpsgrid%clati)

    ! add vertical grid (deleted for horizontal interpolation)
    CALL COPY_NCVAR(intpsgrid%hyai, cgrid%hyai)
    CALL COPY_NCVAR(intpsgrid%hybi, cgrid%hybi)
    intpsgrid%ranges(4,1) = cgrid%ranges(4,1)
    intpsgrid%ranges(4,2) = cgrid%ranges(4,2)

    ! construct intermediate ps grid
    CALL SET_SURFACE_PRESSURE(status, intpsgrid, ps_s(:,:,1,1), lcalc)
    IF (status /= 0) CALL error_bi(grid_error(status), substr)

    CALL COMPLETE_GEOHYBGRID(intpsgrid)

    ! --------------------
    ! make output grid
    ! --------------------
    CALL COPY_GEOHYBGRID(intzogrid, pgrid)
    ! assume congruency of grids ogrid and intgrid
    ! take rotated coordinates for curvi-linear grids
    ! instead of intermediate grid. NCREGRID can only
    ! handle lat-lon-grids

    IF ((QDEF_NCVAR(pgrid%rlonm) .AND. QDEF_NCVAR(pgrid%rlatm))  .OR. &
         (QDEF_NCVAR(pgrid%rloni) .AND. QDEF_NCVAR(pgrid%rlati))) THEN
       ! set non-rotated lat-lon grid
       intzogrid%lonc=.FALSE.
       CALL COPY_NCVAR(intzogrid%lonm, pgrid%rlonm)
       CALL COPY_NCVAR(intzogrid%loni, pgrid%rloni)
       CALL COPY_NCVAR(intzogrid%latm, pgrid%rlatm)
       CALL COPY_NCVAR(intzogrid%lati, pgrid%rlati)
       CALL INIT_NCVAR(intzogrid%rlonm)
       CALL INIT_NCVAR(intzogrid%rloni)
       CALL INIT_NCVAR(intzogrid%rlatm)
       CALL INIT_NCVAR(intzogrid%rlati)
    ELSE
       ! ??? required ??
       ! -- correct also lon/lat definition
       CALL COPY_NCVAR(intzogrid%lonm, pgrid%lonm)
       CALL COPY_NCVAR(intzogrid%loni, pgrid%loni)
       CALL COPY_NCVAR(intzogrid%latm, pgrid%latm)
       CALL COPY_NCVAR(intzogrid%lati, pgrid%lati)
    ENDIF

    ! remove curvilinear information
    ! vertical interpolation proceeds on target grid
    CALL INIT_NCVAR(intzogrid%clonm)
    CALL INIT_NCVAR(intzogrid%cloni)
    CALL INIT_NCVAR(intzogrid%clatm)
    CALL INIT_NCVAR(intzogrid%clati)

    CALL SET_SURFACE_PRESSURE(status, intzogrid, ps_d(:,:,1,1), lcalc)
    IF (status /= 0) CALL error_bi(grid_error(status), substr)

    ! ----------------------------
    ! make vertical interpolation
    ! ----------------------------
    ALLOCATE(dat(SIZE(var_in,1), SIZE(var_in,2), ke_hint, 1))
    dat(:,:,:,1) = var_in(:,:,:,1)

    ALLOCATE(sovar(1))
    ! construct grid with vertical input grid, but horizontal
    CALL RGTOOL_CONVERT_DAT2VAR(sovar(1), dat &
         , ParData(idx)%name%obj, intpsgrid, 'xyzn')
    DEALLOCATE(dat)
    NULLIFY(dat)

    RGT(1) = RG_INT
    CALL REGRID_CONTROL(intpsgrid, intzogrid, sovar, ovar, RGT, .TRUE. &
         , lrgx=.FALSE., lrgy=.FALSE., lrgz=.TRUE., lpresaxis=.FALSE.  &
         , lwork=.TRUE., lstatout=.FALSE.)

    CALL COPY_NCVAR(varo, ovar(1))   ! ... RETURN VARIABLES

    ! CONVERT BACK
    CALL RGTOOL_CONVERT(varo, dat, intzogrid, 'xyzn')
    IF (PRESENT(grid_d)) CALL COPY_GEOHYBGRID(grid_d, intzogrid)
    var_out = dat

    CALL INIT_GEOHYBGRID(intpsgrid)
    CALL INIT_GEOHYBGRID(intzogrid)
    CALL INIT_NCVAR(varo)
    DO ix= 1, SIZE(sovar)
       CALL INIT_NCVAR(sovar(ix))
       CALL INIT_NCVAR(ovar(ix))
    END DO
    DEALLOCATE(ovar)
    NULLIFY(ovar)
    DEALLOCATE(sovar)
    NULLIFY(sovar)
    DEALLOCATE(dat)
    NULLIFY(dat)

  END SUBROUTINE vert_c2p_ncinterpol_old
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  SUBROUTINE vert_c2p_ncinterpol(var_in, var_out, SCRIP_ID, grid_s,grid_d)

    USE messy_main_grid_trafo,    ONLY: RG_INT

    USE MESSY_MAIN_GRID_TRAFO_NRGD_BASE, ONLY: SNREGRID
    USE MESSY_MAIN_GRID_TRAFO_SCRP,      ONLY: CONSTRUCT_VERTICAL_AXIS
    USE MESSY_MAIN_GRID_NETCDF,          ONLY: ERRMSG

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:,:),    POINTER  :: var_in
    REAL(dp), DIMENSION(:,:,:,:),    POINTER  :: var_out
    INTEGER           , INTENT(IN)            :: SCRIP_ID
    TYPE(t_geohybgrid), INTENT(IN)            :: grid_s
    TYPE(t_geohybgrid), INTENT(OUT), OPTIONAL :: grid_d

    INTEGER                               :: nx, ny, nn, nz
    REAL(dp), DIMENSION(:,:,:), POINTER   :: vax_in  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER   :: vax_out => NULL()
    INTEGER                               :: xsize, ysize
    LOGICAL                               :: lrot = .FALSE.
    INTEGER                               :: zdim
    REAL(dp), DIMENSION(:), ALLOCATABLE   :: help
    ! pressure or sigma axis
    LOGICAL                               :: lpresax=.FALSE.
    ! input time
    LOGICAL                               :: lint   =.FALSE.
    INTEGER, DIMENSION(:), POINTER        :: RGT => NULL() ! regridding type

    ALLOCATE(RGT(1))
    RGT = RG_INT

    xsize = SIZE(var_in,1)
    ysize = SIZE(var_in,2)
    var_out = 0

    ! a DEFINITION OF VERTICAL AXIS: IN-GRID
    CALL CONSTRUCT_VERTICAL_AXIS(status, xsize, ysize, vax_in &
         , grid_s, lpresax, SCRIP_ID, pgrid, RGT, lint)
    IF (status /= 0) &
         CALL ERRMSG('CONSTRUCT_VERTICAL_AXIS: ' ,status,1)

    ! a DEFINITION OF VERTICAL AXIS: IN-GRID
    CALL CONSTRUCT_VERTICAL_AXIS(status, xsize, ysize, vax_out &
         , pgrid, lpresax)
    IF (status /= 0) &
         CALL ERRMSG('CONSTRUCT_VERTICAL_AXIS: ' ,status,1)
    ! CHECK if both vertical axis are orientated in the same way
    ! assume axis orientation is equal for all columns, check only
    ! point (1,1)

    IF ( ( (vax_in(1,1,1)-vax_in(1,1,SIZE(vax_in,3))) * &
         (vax_out(1,1,1)-vax_out(1,1,SIZE(vax_out,3)))) < 0) THEN
       ! V-AXIS orientation differs
       ! => AXIS ROTATION in VAX_IN and DAT required
       lrot = .TRUE.
    ELSE
       lrot = .FALSE.
    END IF

    ! ROTATE VAX-IN and DATA IF REQUIRED
    IF (lrot) THEN
       zdim = SIZE(VAX_IN,3)
       ALLOCATE(help(zdim))
       DO nx = 1, SIZE(var_in,1)
          DO ny = 1, SIZE(var_in,2)
             help(:) = VAX_IN(nx,ny,:)
             DO nz = 1, SIZE(help)
                VAX_IN(nx,ny,nz) = help(zdim-nz+1)
             END DO
             DO nn = 1, SIZE(var_in,3)
                help(1:zdim-1) = var_in(nx,ny,:,1)
                DO nz = 1, zdim-1
                   var_in(nx,ny,nz,1) = help(zdim-1-nz+1)
                END DO
             END DO
          END DO
       END DO
       DEALLOCATE(help)
    END IF

    DO nx = 1, SIZE(var_in,1)
       DO ny = 1, SIZE(var_in,2)
             !
          IF (MAXVAL(vax_in(nx,ny,:)) > 0.) THEN

             CALL SNREGRID(status                             &
                  , vax_in(nx,ny,:), vax_out(nx,ny,:)         &
                  , var_in(nx,ny,:,1), var_out(nx,ny,:,1), .FALSE.)
             IF (status /= 0) &
                  CALL ERRMSG('SNREGRID ERROR: ' ,status,24)

          END IF
       END DO
    END DO

    ! free memory
    DEALLOCATE(vax_in)
    NULLIFY(vax_in)
    DEALLOCATE(vax_out)
    NULLIFY(vax_out)

    IF (PRESENT(grid_d)) CALL COPY_GEOHYBGRID(grid_d, pgrid)
    DEALLOCATE(RGT)
    NULLIFY(RGT)

  END SUBROUTINE vert_c2p_ncinterpol
  !-----------------------------------------------------------------------

  SUBROUTINE DEFINE_PRESSURE_COMPONENTS(grid, iis, iie, jjs,jje, kkmin, kkmax)

    USE MESSY_MAIN_GRID_DEF_MEM_BI,  ONLY: nlev
    USE MESSY_MAIN_CHANNEL_ERROR_BI, ONLY: CHANNEL_HALT
    USE MESSY_MAIN_CHANNEL,          ONLY: GET_CHANNEL_OBJECT
    USE MESSY_MAIN_GRID_NETCDF,      ONLY: ERRMSG, INIT_NCVAR         &
                                         , NULL_DIMID, NULL_VARID     &
                                         , VTYPE_DOUBLE   &
                                         , NF90_DOUBLE, INIT_NARRAY
    USE MESSY_MAIN_GRID_TOOLS,       ONLY: RGTOOL_CONVERT_DAT2PREDEFVAR
    USE MESSY_MAIN_GRID_TRAFO,       ONLY: RGMAX
    USE MESSY_MAIN_GRID,             ONLY: GRID_ERROR


    IMPLICIT NONE

    TYPE(t_geohybgrid), INTENT(INOUT) :: grid
    INTEGER,            INTENT(IN)    :: iis, iie, jjs,jje, kkmin, kkmax

    ! LCOAL
    INTEGER :: status
    REAL(dp), DIMENSION(:,:,:,:),  POINTER :: press4d  => NULL()
    REAL(dp), DIMENSION(:,:,:,:),  POINTER :: pressi4d => NULL()

    REAL(dp), DIMENSION(:,:,:,:),  POINTER :: p4d  => NULL()


    ! prepare space for the up-to-the-minte 3D pressure field
    ! required for regridding of pressure coordinates to height coordinates
    ! 3d PRESSURE FIELD INTERFACES
    CALL INIT_NCVAR(grid%hyai)
    CALL INIT_NCVAR(grid%hybi)
    CALL INIT_NCVAR(grid%hyam)
    CALL INIT_NCVAR(grid%hybm)
    CALL INIT_NCVAR(grid%ps)

    IF (.NOT. QDEF_NCVAR(grid%pressi)) THEN
       CALL INIT_NCVAR(grid%pressi)

       grid%pressi%name  = 'pressi'
       grid%pressi%id    = NULL_VARID
       grid%pressi%xtype = NF90_DOUBLE
       ! ... dimensions
       grid%pressi%ndims = 3
       ALLOCATE(grid%pressi%dim(grid%pressi%ndims), STAT=status)
       CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

       grid%pressi%dim(1) = grid%clonm%dim(1)
       grid%pressi%dim(2) = grid%clonm%dim(2)
       grid%pressi%dim(3)%name  = 'ilev'
       grid%pressi%dim(3)%id    = NULL_DIMID
       grid%pressi%dim(3)%len   = nlev+1-chi_kmin+1
       grid%pressi%dim(3)%fuid  = .false.
       grid%pressi%dim(3)%varid = NULL_VARID

       ! ... data
       CALL INIT_NARRAY(grid%pressi%dat, grid%pressi%ndims      &
            ,(/grid%pressi%dim(1)%len, grid%pressi%dim(2)%len   &
            ,  grid%pressi%dim(3)%len/)                      &
            , VTYPE_DOUBLE)
    END IF

    CALL GET_CHANNEL_OBJECT(status, TRIM('COSMO'),TRIM('pressi'),p4=pressi4d)
    CALL CHANNEL_HALT(substr, status)

    p4d => pressi4d(iis:iie, jjs:jje,kkmin:kkmax+1,:)
    CALL RGTOOL_CONVERT_DAT2PREDEFVAR(status, grid%pressi &
         , p4d, 'xyz-', 'xyzn')
    IF (status /= 0) CALL ERROR_BI(grid_error(status), substr)

    ! 3d PRESSURE FIELD MIDS
    IF (.NOT.  QDEF_NCVAR(grid%pressm)) THEN
       CALL INIT_NCVAR(grid%pressm)

       grid%pressm%name  = 'pressm'
       grid%pressm%id    = NULL_VARID
       grid%pressm%xtype = NF90_DOUBLE
       ! ... dimensions
       grid%pressm%ndims = 3
       ALLOCATE(grid%pressm%dim(grid%pressm%ndims), STAT=status)
       CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

       grid%pressm%dim(1) = grid%clonm%dim(1)
       grid%pressm%dim(2) = grid%clonm%dim(2)
       grid%pressm%dim(3)%name  = 'ilev'
       grid%pressm%dim(3)%id    = NULL_DIMID
       grid%pressm%dim(3)%len   =  nlev  -chi_kmin+1
       grid%pressm%dim(3)%fuid  = .false.
       grid%pressm%dim(3)%varid = NULL_VARID

       ! ... data
       CALL INIT_NARRAY(grid%pressm%dat, grid%pressm%ndims          &
            ,(/grid%pressm%dim(1)%len,   grid%pressm%dim(2)%len     &
            ,  grid%pressm%dim(3)%len/)                             &
            , VTYPE_DOUBLE)
    END IF

    CALL GET_CHANNEL_OBJECT(status, TRIM('COSMO'),TRIM('press'),p4=press4d)
    CALL CHANNEL_HALT(substr, status)
    p4d => press4d(iis:iie, jjs:jje,kkmin:kkmax,:)
    CALL RGTOOL_CONVERT_DAT2PREDEFVAR(status, grid%pressm &
         , p4d, 'xyz-', 'xyzn')
    IF (status /= 0) CALL ERROR_BI(grid_error(status), substr)

    NULLIFY(pressi4d)
    NULLIFY(press4d)
    NULLIFY(p4d)

    grid%ranges(3,1) = 0._dp
    grid%ranges(3,2) = RGMAX

  END SUBROUTINE DEFINE_PRESSURE_COMPONENTS
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  ! ----------------------------------------------------------------------


  SUBROUTINE vert_interpol_lm2lm

    !---------------------------------------------------------------------------
    !
    ! Description:
    !   org_vert_inter_lm organizes the interpolation of LM-fields that are
    !   necessary for initial or boundary files for the nonhydrostatic LM.
    !   These computations include:
    !     - for the variables qc_lm and qv_lm, the generalized relative humidity
    !       is interpolated
    !     - special vertical interpolation for pressure pp_lm
    !     - vertical interpolation to LM levels
    !     - pressure deviation and splitting grh_lm in qc_lm and qv_lm
    !
    ! Method:
    !
    !---------------------------------------------------------------------------

    !
    ! MESSy/SMCL
    USE messy_main_constants_mem, ONLY: Rd

    ! INT2COSMO
    USE data_fields_in,        ONLY: hhl_in, p0_in

    USE vgrid_refatm_utils,    ONLY: vcoord_in, refatm_in   &
                                   , vcoord,    refatm, svc1, svc2   &
                                   , reference_atmosphere   &
                                   , reference_atmosphere_2

    IMPLICIT NONE

    ! Subroutine arguments

    ! Local arrays
    REAL(dp)         ::  &
         zhi_hl(ie_in,je_in,ke_hint+1), & ! height on INPUT half levels
         zhi_fl(ie_in,je_in,ke_hint  )    ! height on INPUT full levels

    INTEGER                    ::  &
         zpgr, kloc,                  & ! top of boundary layer for input levels
         k

#ifdef COSMOv5s5
    INTEGER  :: izerror ! status and error status variable
    CHARACTER (LEN=200)        ::  &
         yzerrmsg    ! error message for error handling
#endif

    ! Definition of statement functions
    REAL    (KIND=dp)    ::    &
         zpm, zgdrt, ztdbe, zbetf, zt00

    REAL    (KIND=dp)    ::    &
         zp0_lm_refmod(ie,je,ke_int+1 ), &
         zplm    (ie,je,ke_in+1)

    CHARACTER(LEN=*), PARAMETER :: substr = "vert_interpol_lm2lm"
    INTEGER                     :: kedim, kediff, keindiff
    INTEGER                     :: ii
    LOGICAL                     :: lanalyt_calc_t0p0 = .FALSE.
    REAL(dp)                    :: zhsurfs(ie_in, je_in,2) ! dummy
    REAL(dp)                    :: p0_hint(ie_in, je_in,ke)
    REAL(dp)                    :: zp0hl_hint(ie_in, je_in,ke+1)
    REAL(dp)                    :: rho0_hint(ie_in, je_in,ke)
    REAL(dp)                    :: zt0_hint(ie_in, je_in,ke)
    REAL(dp)                    :: zt0hl_hint(ie_in, je_in,ke)
    REAL(dp)                    :: dp0_hint(ie_in, je_in,ke)

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 0: Initializations
!------------------------------------------------------------------------------

  kloc      = 0
#ifdef COSMOv5s5
  izerror   = 0
  yzerrmsg  = '  '
#endif
!------------------------------------------------------------------------------
! Section 1: Compute generalized relative humidity on horiz. interp. grid
!------------------------------------------------------------------------------

  ! CALCULATE hhl_hint
  IF (refatm%irefatm == 1) THEN
     CALL reference_atmosphere                                   &
          ( hhl_hint, p0_hint, zp0hl_hint, rho0_hint, zt0_hint,  &
          zt0hl_hint, dp0_hint, hsurf_full,                      &
          zhsurfs, ie_in, je_in, ke, refatm,                     &
          vcoord, svc1, svc2, rd, g, lanalyt_calc_t0p0, .TRUE.   &
#ifdef COSMOv5s5
          , yzerrmsg, izerror                                    &
#endif
          )
  ELSE IF (refatm%irefatm == 2) THEN
     CALL reference_atmosphere_2                                 &
          ( hhl_hint, p0_hint, zp0hl_hint, rho0_hint, zt0_hint,  &
            zt0hl_hint, dp0_hint, hsurf_full(:,:),               &
            zhsurfs, ie_in, je_in, ke, refatm,                   &
            vcoord, svc1, svc2, rd, g, .TRUE.                    &
#ifdef COSMOv5s5
          , yzerrmsg, izerror                                    &
#endif
          )
  ENDIF

    IF (lgrhum) CALL CALC_GRHUM(1)

!------------------------------------------------------------------------------
! Section 2: Compute boundary layer height and interpolate pressure
!------------------------------------------------------------------------------

  ! Find the boundary layer height (kzgr)
  zpgr = 850.0E2
  kediff = ke - ke_hint
  keindiff = ke_in - ke_int
  kedim  = MAX(ke_hint,ke_int)
  DO k = ke_hint, 1, -1
     zpm = (vcoord_in%sigm_coord(k+1+kediff)  &
          + vcoord_in%sigm_coord(k+kediff)) * 0.5_dp * refatm_in%p0sl
     IF (zpm > zpgr) kloc = k
  ENDDO
  kloc = kloc - 1

  ! Compute height on full LM INPUT levels
  DO k = 1, ke_hint
     zhi_fl(:,:,k) = 0.5_dp * hhl_hint(:,:,k+kediff+1)  &
                       + 0.5_dp * hhl_hint(:,:,k+kediff)
  ENDDO

  IF (idx_p_pp /= -99) THEN
     ! interpol pressure deviation
     CALL vert_interp (ParData(idx_p_pp)%ptr_hint, ParData(idx_p_pp)%ptr_int &
          , 'pp', kedim, ke_int, zhi_fl, ke_hint, kloc, keindiff)
  END IF
!------------------------------------------------------------------------------
! Section 3: Compute or interpolate the vertical velocity
!------------------------------------------------------------------------------

  DO k = 1, ke_hint+1
     zhi_hl(:,:,k) = hhl_hint(:,:,k+kediff)
  END DO

  IF (idx_p_w /= -99) THEN
     CALL vert_interp (ParData(idx_p_w)%ptr_hint(:,:,:,1)&
          , ParData(idx_p_w)%ptr_int(:,:,:,1)&
          , 'W', kedim+1, ke_int+1, zhi_hl, ke_hint+1, kloc, keindiff)
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Vertical interpolations for the other variables
!------------------------------------------------------------------------------
  DO ii = 1, PD_NUM
     IF (TRIM(ParData(ii)%name%obj) == 'Test_Ar') CYCLE
     IF (lgrhum .AND. (ii == idx_p_qc .OR. ii == idx_p_qi)) CYCLE
     IF (ii == idx_p_pp) CYCLE  ! already interpolated
     IF (ii == idx_p_w) CYCLE   ! already interpolated


     IF (ParData(ii)%rank == 3) THEN
        CALL vert_interp (ParData(ii)%ptr_hint, ParData(ii)%ptr_int &
             , TRIM(ParData(ii)%name%obj), kedim, ke_int, zhi_fl &
             , ke_hint, kloc, keindiff)

     ELSE
        ParData(ii)%ptr_int = ParData(ii)%ptr_hint
     ENDIF

  END DO

!------------------------------------------------------------------------------
! Section 5: Switch from incoming reference atmosphere to outgoing
!------------------------------------------------------------------------------

  IF (refatm_in%irefatm /= refatm%irefatm) THEN

    IF (my_cart_id == 0) THEN
      PRINT *, ' *** Switching from incoming reference atmosphere:  ' &
           , refatm_in%irefatm
      PRINT *, ' *** to outgoing reference atmosphere:              ' &
           , refatm%irefatm
    ENDIF

    ! Compute (incoming) reference atmosphere on outgoing vertical levels
    ! (the formulas for ivctype=2 from SR reference_atmosphere
    ! (for refatm_in%irefatm=1)
    !  or reference_atmosphere_2 (for refatm_in%irefatm=2) are used)

    IF     (refatm%irefatm == 1) THEN

      zgdrt = g/rd/refatm%t0sl
      IF (refatm%dt0lp /= 0.0_dp) THEN
        ztdbe = refatm%t0sl/refatm%dt0lp
      ELSE
        ztdbe = 0.0_dp
      ENDIF
      zbetf = 2.0_dp*refatm%dt0lp*zgdrt/refatm%t0sl

      DO k = 1, ke_int+1
        IF (refatm%dt0lp == 0.0_dp) THEN
          zp0_lm_refmod(:,:,k) = refatm%p0sl  &
               * EXP ( - zgdrt*hhl_in(:,:,k+keindiff) )
        ELSE
          zp0_lm_refmod(:,:,k) = refatm%p0sl * EXP ( - ztdbe*(1.0_dp        &
                      - SQRT(1.0_dp - zbetf*hhl_in(:,:,k+keindiff))) )
        ENDIF

        IF (k > 1) THEN
          ! averaging to main levels
          zp0_lm_refmod(:,:,k-1) = &
               0.5_dp * (zp0_lm_refmod(:,:,k-1) + zp0_lm_refmod(:,:,k))
        ENDIF
      ENDDO

    ELSEIF (refatm%irefatm == 2) THEN

      zt00 = refatm%t0sl - refatm%delta_t

      DO k = 1, ke_int
         zp0_lm_refmod(:,:,k) = &
              refatm%p0sl*EXP ( - g/rd*refatm%h_scal/zt00 * LOG(  &
              (EXP(0.5_dp*(hhl_in(:,:,k+keindiff) + hhl_in(:,:,k+keindiff+1) ) &
               /refatm%h_scal)*zt00 + refatm%delta_t)/(zt00 + refatm%delta_t)) )
      ENDDO

    ENDIF

    ! Compute full pressure on outgoing vertical levels
    DO k = 1, ke_int
      zplm(:,:,k) = zp0_lm_refmod(:,:,k) + ParData(idx_p_pp)%ptr_int(:,:,k,1)
    ENDDO

    ! Compute pressure deviation from (outgoing) reference atmosphere pressure
    ! on outgoing vertical levels
    DO k = 1, ke_in
      ParData(idx_p_pp)%ptr_int(:,:,k,1) = zplm(:,:,k) - p0_in(:,:,k+keindiff)
    ENDDO

  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE vert_interpol_lm2lm


!===============================================================================
!+ Vertical Interpolation to LM levels.
!------------------------------------------------------------------------------

SUBROUTINE vert_interp (xlm_in, xlm_out, yname, idim_max, idim_out, height &
                      , idim_in, kgr, keindiff)

!------------------------------------------------------------------------------
!
! Description:
!   Does the vertical interpolation on the sigma reference levels of the
!   variables: u_lm, v_lm, w_lm, t_lm, pp_lm, grh_lm and qi_lm
!   (grh = generalized relative humidity)
!
! Method:
!   The routine is entered with the variable xlm on the height levels of the
!   input profile (coarse grid). An intermediate profile is constructed by
!   shifting the input profile in order to take into account the topographical
!   difference between the coarse and the fine grid models (variable zxexp_vec
!   on the levels zpexp_vec).
!   For each gridpoint the tension spline routine tautsp is called for
!   the interpolation.
!   For break(i) <= x <= break(i+1) the interpolation function has the form
!      F(X) = COEF(1,I)+DX(COEF(2,I)+DX/2(COEF(3,I)+DX/3(COEF(4,I)))
!      using DX=X-BREAK(I) and i=1,...,izl
!   This interpolation leads to the output profile with the variable xlm on the
!   fine grid levels zhhl.
!
!------------------------------------------------------------------------------

  USE data_fields_in,    ONLY: hhl_in

  USE messy_main_tools,  ONLY: tautsp2D


  IMPLICIT NONE



! Subroutine arguments
  INTEGER,   INTENT(IN)   ::  &
       idim_max,     &
       idim_out,     & ! vertical dimension of outgoing LM fields
       idim_in,      & ! vertical dimension of incoming LM fields
       kgr,          & ! height level of boundary layer top
       keindiff
  ! height on incoming LM levels used as abscissas
  REAL (KIND=dp),    INTENT(IN)  :: height (ie_in,je_in,idim_in)

  ! field to be interpolated
  REAL (KIND=dp),    INTENT(IN)  :: xlm_in (ie_in,je_in,idim_in)
  ! interpolated field
  REAL (KIND=dp),    INTENT(OUT) :: xlm_out (ie_in,je_in,idim_out)

  CHARACTER (LEN=*), INTENT(IN)  :: yname           ! name of the variable

!------------------------------------------------------------------------------

! Local scalars:
INTEGER       :: i, j, k, kk, k1, k2, k3, kstart, kend, izerror, idone,    &
                 izln_vec(ie_in),        & ! number of abscissas for spline
                 nztau_vec(ie_in),       & ! number of abscissas
                 kzgrn_vec(ie_in),       & !
                 izn      (ie_in),       & !
                 izind_vec(ie_in,ke_in+1)

LOGICAL       :: ldone(ie_in)

REAL(KIND=dp) :: zdx, zp1(ie_in), zp2(ie_in), zp3(ie_in), &
                 zcheck(ie_in),              & !
                 zdelh(ie_in), zdelhchk,     & ! height difference
                 zgamma                        ! tension-parameter = 5.5

! Local arrays:
REAL(KIND=dp) :: zhhl       (ie_in,ke_in+1),   & !
                 zpexp_vec  (ie_in,idim_in+6), & !
                 zxexp_vec  (ie_in,idim_in+6), & !
                 zbreak_vec (ie_in,3*idim_max)    ! abscissas of spline

REAL(KIND=dp) :: zcoef_vec (ie_in,4,3*idim_max),& ! coefficients  of spline
                 zs_vec    (ie_in,idim_in+6,6), & ! work array for tautsp
                 zbx       (ie_in,je_in),       & !
                 zdelx     (ie_in,je_in),       & !
                 zpq       (ie_in,je_in),       & ! work arrays for regression
                 zpq2      (ie_in,je_in),       & !
                 zxq       (ie_in,je_in),       & !
                 zxpq      (ie_in,je_in)

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 0: Initializations
!------------------------------------------------------------------------------

  izerror   = 0
  zgamma    = 5.5
  zdelhchk  = 1.0E-4_dp

  ! Initialize intervals
  izind_vec(:,:) = 0.0

  ! Initializations
  nztau_vec(:)=0.0
  zbx(:,:)    =0.0
  zdelx(:,:)  =0.0
  zdelh       =0.0
  zpexp_vec(:,:)=-23456789
  zxexp_vec(:,:)=-23456789

  kstart = 1
  kend = idim_out

  CALL debug_bi(submodstr, &
       'vertical interpolation of:  '//TRIM(yname), '', substr)

!------------------------------------------------------------------------------
! Section 1: Loop over all grid points
!------------------------------------------------------------------------------

  ! Construct intermediate profiles that can be used for interpolation:
  !   variable zxexp_vec on the levels zpexp_vec.
  ! The levels in zpexp_vec are numbered from the ground to the top of the
  ! atmosphere (just the inverse counting compared to the atmospheric
  ! COSMO-Model variables). The first levels to be set are:
  !   zpexp_vec(:,1): under the ground
  !   zpexp_vec(:,2): at the soil surface
  !   zpexp_vec(:,3): 1st atmospheric level
  ! All other levels are taken from the coarse grid input model
  !   zpexp_vec(:,4:idim_in+2):  height(:,j;idim_in+1-2:1)

  DO j = 1, je_in

    ! 1.1 Lower boundary values
    ! -------------------------
    ! The first value that is set is the surface value; for w, an artificial
    ! level of 20 m below the fine grid surface is introduced
    ! For u and v-gridpoints, the height is averaged to the C-grid.
    ! On the eastern and northern boundary the values from the last mass grid
    ! point are used. For subdomains with a right or upper neighbor, this does
    ! not matter, for the rightmost subdomain, it is the only value available

    IF     (yname == 'W') THEN
       zpexp_vec(:,2) =  hhl_in(:,j,ke_in+1)-20._dp
    ELSEIF (yname == 'U') THEN
      DO i = 1, ie_in-1
         zpexp_vec(i,2) = &
              0.5_dp * (hhl_in(i,j,ke_in+1) + hhl_in(i+1,j,ke_in+1))
      ENDDO
      zpexp_vec(ie_in,2) = hhl_in(ie_in,j,ke_in+1)
    ELSEIF (yname == 'V') THEN
      IF (j < je_in) THEN
        zpexp_vec(:,2) = 0.5_dp * (hhl_in(:,j,ke_in+1) &
                                   + hhl_in(:,j+1,ke_in+1))
      ELSE
        zpexp_vec(:,2) = hhl_in(:,je_in,ke_in+1)
      ENDIF
    ELSE   ! t, grh, qi, pp and chemistry variables
      zpexp_vec(:,2) = hhl_in(:,j,ke_in+1)
    ENDIF

    kzgrn_vec(:) = 0

    ! For all variables, an artificial value of 50 m below the surface is
    ! introduced for the intermediate profile
    IF     (yname == 'W') THEN
      zpexp_vec(:,1) = zpexp_vec(:,2)- 30.0_dp   ! overall: 30 m below surface
    ELSE
      zpexp_vec(:,1) = zpexp_vec(:,2)- 50.0_dp   ! also: 50 m below surface
    ENDIF

    IF     (yname == 'U') THEN
      DO i = 1, ie_in-1
        zpexp_vec(i,3) = 0.5_dp * (height(i,    j,idim_in) &
                                   + height(i+1,  j,idim_in))
      ENDDO
      zpexp_vec(ie_in,3) = height(ie_in,j,idim_in)
    ELSEIF (yname == 'V') THEN
      IF (j < je_in) THEN
        zpexp_vec(:,3) = 0.5_dp * (height(:,    j,idim_in) &
                                  + height(:,  j+1,idim_in))
      ELSE
        zpexp_vec(:,3) = height(:,je_in,idim_in)
      ENDIF
    ELSE   ! w, t, grh, qi, pp and chemistry variables
      zpexp_vec(:,3) =  height(:,    j,idim_in)
    ENDIF

    IF (yname == 'T') THEN
      zxexp_vec(:,1) = ParData(idx_p_ts)%ptr_hint(:,j,1,1)
      zxexp_vec(:,2) = ParData(idx_p_ts)%ptr_hint(:,j,1,1)
      zxexp_vec(:,3) = xlm_in(:,j,idim_in)
    ELSE IF (yname == 'U') THEN
      zxexp_vec(:,1) = 0.0
      zxexp_vec(:,2) = 0.0
      zxexp_vec(:,3) = xlm_in(:,j,idim_in)
    ELSE IF (yname == 'V') THEN
      zxexp_vec(:,1) = 0.0
      zxexp_vec(:,2) = 0.0
      zxexp_vec(:,3) = xlm_in(:,j,idim_in)
    ELSE
      zxexp_vec(:,1) = xlm_in(:,j,idim_in)
      zxexp_vec(:,2) = xlm_in(:,j,idim_in)
      zxexp_vec(:,3) = xlm_in(:,j,idim_in)
    END IF

    ! 1.2 Levels in the atmosphere
    ! ----------------------------
    DO k = 2, idim_in
      ! For u and v-gridpoints, the height is averaged to the C-grid.
      ! On the eastern and northern boundary the values from the last mass grid
      ! point are used. For subdomains with a right or upper neighbor, this does
      ! not matter, for the rightmost subdomain, it is the only value available
      IF     (yname == 'U') THEN
        DO i = 1, ie_in-1
          zpexp_vec(i,k+2) = 0.5_dp *(height(i,    j,idim_in+1-k) &
                                     + height(i+1,  j,idim_in+1-k))
        ENDDO
        zpexp_vec(ie_in,k+2) = height(ie_in,j,idim_in+1-k)
      ELSEIF (yname == 'V') THEN
        IF (j < je_in) THEN
          zpexp_vec(:,k+2) = 0.5_dp *(height(:,    j,idim_in+1-k) &
                                     + height(:,  j+1,idim_in+1-k))
        ELSE
          zpexp_vec(:,k+2) = height(:,je_in,idim_in+1-k)
        ENDIF
      ELSE
        zpexp_vec  (:,k+2) = height(:,    j,idim_in+1-k)
      ENDIF
      zxexp_vec    (:,k+2) = xlm_in(:,    j,idim_in+1-k)
    ENDDO

    ! 1.3 Set the height vector for the output model at the correct grid points
    !     (output profiles)
    ! -----------------------------------------------------------------------
    IF     (yname == 'W') THEN
      DO k = 1, kend
        zhhl(:,k) = hhl_in(:,j,k+keindiff)
      ENDDO
    ELSEIF (yname == 'U') THEN
      DO k = 1, kend
        DO i = 1, ie_in-1
          zhhl(i,k)=0.25_dp*(hhl_in(i,j,k+keindiff) + hhl_in(i+1,j,k+keindiff) &
                         + hhl_in(i,j,k+keindiff+1) + hhl_in(i+1,j,k+keindiff+1))
        ENDDO
        zhhl(ie_in,k)= &
             0.5_dp * (hhl_in(ie_in,j,k+keindiff) + hhl_in(ie_in,j,k+keindiff+1))
      ENDDO
    ELSEIF (yname == 'V') THEN
      DO k = 1, kend
        IF (j < je_in) THEN
          zhhl(:,k) = 0.25_dp*(hhl_in(:,j,k+keindiff) + hhl_in(:,j+1,k+keindiff)&
                         + hhl_in(:,j,k+keindiff+1) + hhl_in(:,j+1,k+keindiff+1))
        ELSE
          zhhl(:,k) = &
               0.5_dp*(hhl_in(:,je_in,k+keindiff) + hhl_in(:,je_in,k+keindiff+1))
        ENDIF
      ENDDO
    ELSE
      DO k = 1, kend
        zhhl(:,k) = 0.5_dp * (hhl_in(:,j,k+keindiff) + hhl_in(:,j,k+keindiff+1))
      ENDDO
    ENDIF

    nztau_vec(:) = idim_in + 2

    !--------------------------------------------------------------------------
    ! For the temperature a linear approximation of the vertical profile above
    ! the boundary layer is computed.
    ! This linear profile is used for the extrapolation of the fields.
    !
    ! The index 'kgr' of the boundary layer top.
    !
    ! The linear approximation is of the form:
    !  x(h) = bx*h + cx  ,  where h is the height of input levels.
    !
    ! We need only the coefficient bx and delx . We do not need cx.
    !
    ! Version 1.9:
    !  The computation of this regression had a severe bug before, where
    !  the index kgr was used for the intermediate profile, but kgr refers
    !  to the indexing in the initial profile. This has been corrected by
    !  Anne Roches, MCH.
    !
    !  Also, the regression is now computed only for temperature, not for the
    !  horizontal wind speeds any more.
    !
    !--------------------------------------------------------------------------
     IF (yname == 'T') THEN
      zpq  (:,j) = 0.0
      zpq2 (:,j) = 0.0
      zxq  (:,j) = 0.0
      zxpq (:,j) = 0.0
      DO k = 3+idim_in-kgr,3+idim_in-kgr+3
        zpq (:,j) = zpq (:,j) + zpexp_vec(:,k)
        zpq2(:,j) = zpq2(:,j) + zpexp_vec(:,k) * zpexp_vec(:,k)
        zxq (:,j) = zxq (:,j) + zxexp_vec(:,k)
        zxpq(:,j) = zxpq(:,j) + zxexp_vec(:,k) * zpexp_vec(:,k)
      ENDDO
      zbx  (:,j) =( zxpq(:,j) - zxq(:,j)*zpq(:,j)/4.0_dp) /    &
                   (zpq2(:,j) - zpq(:,j)*zpq(:,j)/4.0_dp)
      zdelx(:,j) =  zbx (:,j) * (zhhl(:,kend) - zpexp_vec(:,3))
    ENDIF

    ! 1.4 Shift boundary layer depending on height differences
    ! --------------------------------------------------------

    DO i = 1, ie_in
      ! Height difference LM minus COARSE LM (horizontaly interpolated)
      zdelh(i)  = zhhl(i,kend) - zpexp_vec(i,3)
      zcheck(i) = zpexp_vec(i,3)
    ENDDO

    ! Shifting within the boundary layer:
    ! Extrapolation of the profiles for hhl_gl > zhhl :
    !  for yname= 't',         with a linear regression
    !  for all other variables with a constant value.
    IF (yname == 'T') THEN
      DO k = idim_in, kgr+1, -1
        kk = 3+idim_in-k
        DO i = 1, ie_in
          IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
            zpexp_vec(i,kk) = zpexp_vec(i,kk) + zdelh(i)
            zxexp_vec(i,kk) = xlm_in (i,j,k) + zdelx(i,j) ! shift all PBL
          ENDIF
        ENDDO
      ENDDO
    ELSE  ! 'rh', 'qi', 'pp', 'p' and 'w'
      DO k = idim_in, kgr+1, -1
        kk = 3+idim_in-k
        DO i = 1, ie_in
          IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
            zpexp_vec(i,kk) = zpexp_vec(i,kk) + zdelh(i)
            zxexp_vec(i,kk) = xlm_in (i,j,k)        ! shift PBL height ONLY
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    ! Shift ground values
    IF (yname == 'T') THEN
        DO i = 1, ie_in
          zxexp_vec(i,2)= &
               zxexp_vec(i,2)+zbx(i,j)*(zpexp_vec(i,2)-hhl_hint(i,j,ke_in+1))
          zxexp_vec(i,1)= zxexp_vec(i,2)
        ENDDO
    ENDIF

    ! At the top of the boundary layer
    k1 = 3+idim_in-kgr+2
    k2 = 3+idim_in-kgr+1
    k3 = 3+idim_in-kgr
    DO i = 1, ie_in
       IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
          zp1(i) = zpexp_vec(i,k1)
          zp2(i) = zpexp_vec(i,k2)
          zp3(i) = zpexp_vec(i,k3)
          zpexp_vec(i,k1) = zp1(i) + zdelh(i)/4.0_dp * 1.0_dp
          zpexp_vec(i,k2) = zp2(i) + zdelh(i)/4.0_dp * 2.0_dp
          zpexp_vec(i,k3) = zp3(i) + zdelh(i)/4.0_dp * 3.0_dp
       ENDIF
    ENDDO

    IF (yname == 'T') THEN
       DO i = 1, ie_in
          IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
             ! add 2 regression values
             zxexp_vec(i,k1) = xlm_in(i,j,kgr-2) &
                               + zbx(i,j)*zdelh(i)/4.0_dp * 1.0_dp
             zxexp_vec(i,k2) = xlm_in(i,j,kgr-1) &
                               + zbx(i,j)*zdelh(i)/4.0_dp * 2.0_dp
             zxexp_vec(i,k3) = xlm_in(i,j,kgr  ) &
                               + zbx(i,j)*zdelh(i)/4.0_dp * 3.0_dp
          ENDIF
       ENDDO
    ELSE   ! 'rh', 'qi', 'pp' and 'w'
       DO i = 1, ie_in
          IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
             ! add transition zone values
             zxexp_vec(i,k1) = xlm_in(i,j,kgr-2)
             zxexp_vec(i,k2) = xlm_in(i,j,kgr-1)
             zxexp_vec(i,k3) = xlm_in(i,j,kgr  )
          ENDIF
       ENDDO
    ENDIF

    ! upper part of the profile

    DO k = kgr-3, 1, -1              ! TOP
       kk = 3+idim_in-k
       DO i = 1, ie_in
          IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
             zp1(i) = zpexp_vec(i,kk)   ! no change
             zpexp_vec(i,kk) = zp1(i)
             zxexp_vec(i,kk) = xlm_in(i,j,k)
          ENDIF
       ENDDO
    ENDDO

    DO i = 1, ie_in
      IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
        zpexp_vec(i,idim_in + 4) = zp1(i) + 50.0_dp
        zpexp_vec(i,idim_in + 3) = zp1(i) + 20.0_dp
        zxexp_vec(i,idim_in + 4) = xlm_in(i,j,1)
        zxexp_vec(i,idim_in + 3) = xlm_in(i,j,1)
        nztau_vec(i) = idim_in + 4
      ENDIF
    ENDDO

    ! Elevate (shift) the boundary layer for hhl_gl < zhhl :
    DO i = 1, ie_in
      IF (zcheck(i) < zhhl(i,kend) - zdelhchk) THEN
        k = 3+idim_in-kgr
        DO WHILE (zpexp_vec(i,k) < zpexp_vec(i,3+idim_in-kgr) + zdelh(i))
          k = k + 1
        ENDDO
        kzgrn_vec(i) = k-3-idim_in+kgr
      ENDIF
    ENDDO

    IF (yname == 'T') THEN
      DO k = idim_in, kgr, -1
        kk = 3+idim_in-k
        DO i = 1, ie_in
          IF (zcheck(i) < zhhl(i,kend) - zdelhchk) THEN
            zpexp_vec(i,kk) = zpexp_vec(i,kk) + zdelh(i)
            zxexp_vec(i,kk) = xlm_in  (i,j,k) + zdelx(i,j) ! shift all PBL
          ENDIF
        ENDDO
      ENDDO
    ELSE   ! 'rh', 'qi', 'pp' and 'w'
       DO k = idim_in, kgr, -1
          kk = 3+idim_in-k
          DO i = 1, ie_in
             IF (zcheck(i) < zhhl(i,kend) - zdelhchk) THEN
                zpexp_vec(i,kk) = zpexp_vec(i,kk) + zdelh(i)
                zxexp_vec(i,kk) = xlm_in (i,j,k)        ! shift PBL height ONLY
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    DO k = kgr-1, 1, -1     ! LEAVE OUT zdelh (kzgrn_vec(i)) LEVELS
      DO i = 1, ie_in
        IF (k >= kzgrn_vec(i)+1) THEN
          kk = 3+idim_in-k
          IF (zcheck(i) < zhhl(i,kend) - zdelhchk) THEN
            zxexp_vec(i,kk) = xlm_in   (i,j,k - kzgrn_vec(i))
            zpexp_vec(i,kk) = zpexp_vec(i,kk + kzgrn_vec(i))
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    DO i = 1, ie_in
      IF (zcheck(i) < zhhl(i,kend) - zdelhchk) THEN
        nztau_vec(i)= idim_in + 2 - kzgrn_vec(i)
      ENDIF
    ENDDO

    ! 1.5. Upper boundary values
    ! --------------------------
    DO i = 1, ie_in
      k = nztau_vec(i) + 2
      nztau_vec(i) = k
      zpexp_vec(i,k  ) = zpexp_vec(i,k-2) + 2000.0_dp
      zpexp_vec(i,k-1) = zpexp_vec(i,k-2) + 1000.0_dp
      zxexp_vec(i,k  ) = xlm_in(i,j,1)
      zxexp_vec(i,k-1) = xlm_in(i,j,1)
    ENDDO


    izln_vec(:) = 3*idim_max

    CALL tautsp2D(zpexp_vec, zxexp_vec, nztau_vec, ie_in, 1, ie_in, idim_in+6, &
                  zgamma, zs_vec, zbreak_vec, zcoef_vec, izln_vec, izerror)

    IF (izerror == 0) THEN

      ! 2.1 Locate intervals and put in izind
      izn(:) = izln_vec(:) - 1
      DO k = kstart, kend
        ldone(:) = .FALSE.
        idone    = 0
        DO WHILE (idone < ie_in)
          DO i = 1, ie_in
            IF (.NOT. ldone(i)) THEN
             IF ( (zbreak_vec(i,izn(i)) <= zhhl(i,k))             .AND. &
                               (zhhl(i,k) <= zbreak_vec(i,izn(i)+1)) ) THEN
                izind_vec(i,k) = izn(i)
                ldone(i)       = .TRUE.
                idone          = idone + 1
              ELSE
                izn(i) = izn(i) - 1
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ! 2.2 Compute interpolated values
      !
      DO k = kstart, kend
        DO i = 1, ie_in
          zdx        = zhhl(i,k) - zbreak_vec(i,izind_vec(i,k))

          xlm_out(i,j,k) = zcoef_vec(i,1,Izind_vec(i,k)) +          &
                           zdx*(zcoef_vec(i,2,izind_vec(i,k)) +     &
                           zdx*0.5_dp*(zcoef_vec(i,3,izind_vec(i,k)) + &
                           zdx/3.0_dp*zcoef_vec(i,4,izind_vec(i,k))))
        ENDDO
      ENDDO

      ! Limit values of qi
      IF (yname=='QI' .OR. yname=='QR' .OR. yname=='QS' .OR. yname=='QG' &
            .OR. yname=='QV' .OR. yname=='QC') THEN
         DO k = 1, kend
            DO i = 1, ie_in
               IF (xlm_out(i,j,k) < 0.0_dp) THEN
                  xlm_out(i,j,k) = 0.0_dp
               ENDIF
            ENDDO
         ENDDO
     ENDIF

      ! Set w_lm to 0.0 at the top and the bottom
      IF (yname == 'W') THEN
        DO i = 1, ie_in
          xlm_out(i,j,    1)    = 0.0_dp
          xlm_out(i,j,idim_out) = 0.0_dp
        ENDDO
      ENDIF

   ELSE
      PRINT *, '*** ERROR in tautsp2D: while processing j-index:  ', j
      CALL error_bi ('Error in tautsp', substr)
   ENDIF

ENDDO ! j = 1, je_in

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE vert_interp

  END SUBROUTINE interpol_parent_data
!===============================================================================
  ! PARENT COUPLING -
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE CALC_lonlat_lm

    ! MESSy/BMIL
    USE messy_main_mpi_bi,     ONLY: my_cart_id
    ! MESSy/SMCL
    USE messy_main_grid_trafo, ONLY: rlarot2rla, phirot2phi
    ! INT2COSMO
    USE data_int2lm_parallel,  ONLY: isubpos, nboundlines
    USE data_grid_lm,          ONLY: ie2lm, je2lm, startlon_tot, startlat_tot &
                                   , pollon, pollat, polgam, dlon, dlat
    USE data_fields_lm,        ONLY: geolon_lm,  geolat_lm,  lon_lm,  lat_lm &
                                   , geoloni_lm, geolati_lm, loni_lm, lati_lm

    IMPLICIT NONE

    INTEGER                     :: ix, jx, itot, jtot
    INTEGER                     :: itmp
    INTEGER                     :: istartlon, istartlat, idlon, idlat


    istartlon = NINT(startlon_tot * RCF)
    istartlat = NINT(startlat_tot * RCF)
    idlon     = NINT(dlon         * RCF)
    idlat     = NINT(dlat         * RCF)

    itot = isubpos(my_cart_id,1) - nboundlines - 1
    jtot = isubpos(my_cart_id,2) - nboundlines - 1

    DO ix = 1, ie2lm
       DO jx = 1, je2lm
          ! Note: startlon_tot is the longitude of the first COSMO grid box,
          !       i.e., the second INT2LM grid box (=> additionally -1)
          !       for lstartlat_tot respectively
          !lon_lm(ix,jx) = startlon_tot + (itot+ix-2) * dlon
          !lat_lm(ix,jx) = startlat_tot + (jtot+jx-2) * dlat
          itmp = istartlon + (itot+ix-2) * idlon
          lon_lm(ix,jx) = REAL(itmp,dp) / RCF

          itmp = istartlat + (jtot+jx-2) * idlat
          lat_lm(ix,jx) = REAL(itmp,dp) / RCF

          geolon_lm(ix,jx) = rlarot2rla ( lat_lm(ix,jx) , lon_lm(ix,jx) &
               , pollat, pollon, polgam)
          geolat_lm(ix,jx) = phirot2phi ( lat_lm(ix,jx) , lon_lm(ix,jx) &
               , pollat, polgam)

       END DO
    END DO

    DO ix = 1, ie2lm+1
       DO jx = 1, je2lm+1
          itmp = istartlon + (itot+ix-2) * idlon - NINT(0.5_dp * dlon * RCF)
          loni_lm(ix,jx) = REAL(itmp,dp) / RCF

          itmp = istartlat + (jtot+jx-2) * idlat - NINT(0.5_dp * dlat * RCF)
          lati_lm(ix,jx) = REAL(itmp,dp) / RCF

          geoloni_lm(ix,jx) = rlarot2rla ( lati_lm(ix,jx) , loni_lm(ix,jx) &
               , pollat, pollon, polgam)
          geolati_lm(ix,jx) = phirot2phi ( lati_lm(ix,jx) , loni_lm(ix,jx) &
               , pollat, polgam)

       END DO
    END DO

  END SUBROUTINE CALC_lonlat_lm
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE CALC_forward_weights

    ! MESSy/SMCL
    USE messy_main_grid,                 ONLY: GRID_ERROR, COPY_GEOHYBGRID
    USE messy_main_grid_netcdf,          ONLY: MAIN_GRID_SET_MESSAGEMODE! &
                                             !, MSGMODE_E
    USE messy_main_grid_trafo,           ONLY: RG_INT, SWITCH_GEOHYBGRID
    USE messy_main_grid_trafo_scrp,      ONLY: CALC_SCRIPDATA    &
                                             , CALC_SCRIP_WEIGHTS
    USE messy_main_grid_trafo_scrp_base, ONLY: norm_opt_dstarea  &
                                             , map_type_conserv

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'calc_forward_weights'
    INTEGER                     :: status
    LOGICAL                     :: lok
    INTEGER, DIMENSION(1)       :: RGTYPE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!    CALCULATE FORWARD WEIGHTS   !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    CALL MAIN_GRID_SET_MESSAGEMODE(MSGMODE_E)
    CALL MAIN_GRID_SET_MESSAGEMODE()

    ! 1a) DEFINE (TOTAL) CHILD GRID
    IF (L_gridrotParenteqChild) THEN
       CALL define_rotINT2COSMO_grid(i2cgrid, lok)
    ELSE
       CALL define_INT2COSMO_grid(i2cgrid, lok)
    ENDIF

    ! 1b) determine horizontal child grid
    CALL COPY_GEOHYBGRID(i2chgrid,i2cgrid)
    CALL SWITCH_GEOHYBGRID(i2chgrid, .TRUE., .TRUE., .FALSE.)

    ! 1c) define staggered u grid
    CALL COPY_GEOHYBGRID(i2cUgrid, i2cgrid)
    CALL redefine_staggered_Ugrid(i2cUgrid, lok)
    IF (.NOT. lok) CALL error_bi(' ERROR in i2cUgrid definition', substr)
    CALL COPY_GEOHYBGRID(i2chUgrid,i2cUgrid)
    CALL SWITCH_GEOHYBGRID(i2chUgrid, .TRUE., .TRUE., .FALSE.)

    ! 1d) define staggered v grid
    CALL COPY_GEOHYBGRID(i2cVgrid, i2cgrid)
    CALL redefine_staggered_Vgrid(i2cVgrid, lok)
    IF (.NOT. lok) CALL error_bi(' ERROR in i2cVgrid definition', substr)
    CALL COPY_GEOHYBGRID(i2chVgrid,i2cVgrid)
    CALL SWITCH_GEOHYBGRID(i2chVgrid, .TRUE., .TRUE., .FALSE.)

    ! 2a) DEFINE (TOTAL) PARENT GRID
    IF (L_gridrotParenteqChild) THEN
       CALL define_rot_parentin_grid(pingrid, lok)
    ELSE
       CALL define_parentin_grid(pingrid, lok)
    ENDIF

    ! 2b) DEFINE horizontal PARENT GRID
    CALL COPY_GEOHYBGRID(pinhgrid, pingrid)
    CALL SWITCH_GEOHYBGRID(pinhgrid, .TRUE., .TRUE., .FALSE.)

    ! 3.) convert to SCRIP GRID
    RGTYPE(1) = RG_INT ! dummy

    ! 3a) mid point grid (unstaggered grid)
    CALL CALC_SCRIPDATA(status, pingrid, i2cgrid, RGTYPE            &
         , I2C_SD_ID, PSD=I2C_SD                                    &
         , norm_opt_in=norm_opt_dstarea, map_type_in=map_type_conserv)
    IF (status /= 0) THEN
       IF (status /= 01 )  CALL error_bi(grid_error(status), substr)
    ELSE
       ! CALCULATE WEIGHTS
       CALL CALC_SCRIP_WEIGHTS(status, I2C_SD)
    END IF

    ! 3b) u grid (staggered grid)
    CALL CALC_SCRIPDATA(status, pingrid, i2cUgrid, RGTYPE            &
         , I2C_U_SD_ID, PSD=I2C_U_SD                                 &
         , norm_opt_in=norm_opt_dstarea, map_type_in=map_type_conserv)
    IF (status /= 0) THEN
       IF (status /= 01 )  CALL error_bi(grid_error(status), substr)
    ELSE
       ! CALCULATE WEIGHTS
       CALL CALC_SCRIP_WEIGHTS(status, I2C_U_SD)
    END IF

    ! 3b) v grid (staggered grid)
    CALL CALC_SCRIPDATA(status, pingrid, i2cVgrid, RGTYPE           &
         , I2C_V_SD_ID, PSD=I2C_V_SD                                &
         , norm_opt_in=norm_opt_dstarea, map_type_in=map_type_conserv)
    IF (status /= 0) THEN
       IF (status /= 01 )  CALL error_bi(grid_error(status), substr)
    ELSE
       ! CALCULATE WEIGHTS
       CALL CALC_SCRIP_WEIGHTS(status, I2C_V_SD)
    END IF

    !*************************************************************
    !*************************************************************

 CONTAINS

   SUBROUTINE define_parentin_grid(g,ok)
    !
    ! remove halo to avoid double counting of halo area

     ! INT2COSMO
     USE data_int2lm_parallel,    ONLY: isubpos_in, nboundlines_in, my_cart_id
     USE data_grid_in,            ONLY: startlon_in_tot, startlat_in_tot       &
                                      , dlon_in, dlat_in                       &
                                      , polgam_in, pollon_in, pollat_in        &
                                      , ie_in, je_in, ke_in, ak_in, bk_in
     ! MESSy/SMCL
     USE messy_main_timer,        ONLY: YEAR_START, MONTH_START, DAY_START     &
                                      , HOUR_START, MINUTE_START, SECOND_START &
                                      , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

     USE messy_main_constants_mem, ONLY: sp
     USE messy_main_timer,         ONLY: time_span_d

     ! GRID MODULES
     USE messy_main_grid_trafo,    ONLY: RGEMPTY, COMPLETE_GEOHYBGRID &
                                       , PHIROT2PHI, RLAROT2RLA
     USE messy_main_grid_netcdf,   ONLY: ERRMSG, COPY_NCDIM       &
                                       , null_dimid, null_varid   &
                                       , VTYPE_REAL, VTYPE_DOUBLE &
                                       , nf90_float, nf90_double  &
                                       , POSITION, INIT_NARRAY, ADD_NCATT
     USE messy_main_grid,          ONLY: INIT_GEOHYBGRID

     IMPLICIT NONE

     ! I/O
     TYPE(t_geohybgrid),  INTENT(OUT) :: g
     LOGICAL            , INTENT(OUT) :: ok

     !LOCAL
     REAL(DP)            :: dts
     INTEGER             :: i,j, n, itot, jtot
     INTEGER             :: status
     CHARACTER(LEN=100)  :: tunit
     REAL(dp)            :: clon, clat
     REAL(dp)            :: tmp
     INTEGER             :: itmp
     INTEGER             :: istartlon, istartlat, idlon, idlat

     istartlon = NINT(startlon_in_tot * RCF_IN)
     istartlat = NINT(startlat_in_tot * RCF_IN)
     idlon     = NINT(dlon_in  *RCF_IN)
     idlat     = NINT(dlat_in  *RCF_IN)

     ! INIT
     CALL INIT_GEOHYBGRID(g)

     g%name = 'Parent INGRID'

     g%file = ' '       ! Filename
     g%t    = 0                               ! time step

     ! Curvilinear LONGITUDE (MID) ...
     g%clonm%name  = 'lon'
     g%clonm%id    = NULL_VARID
     g%clonm%xtype = NF90_DOUBLE
     g%clonc = .FALSE.
     ! ... dimensions
     g%clonm%ndims = 2
     ALLOCATE(g%clonm%dim(g%clonm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%clonm%dim(1)%name  = 'lon'
     g%clonm%dim(1)%id    = NULL_DIMID
     g%clonm%dim(1)%len   = ie_in
     g%clonm%dim(1)%fuid  = .false.
     g%clonm%dim(1)%varid = NULL_VARID
     g%clonm%dim(2)%name  = 'lat'
     g%clonm%dim(2)%id    = NULL_DIMID
     g%clonm%dim(2)%len   = je_in
     g%clonm%dim(2)%fuid  = .false.
     g%clonm%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%clonm%dat, g%clonm%ndims      &
          , (/g%clonm%dim(1)%len,g%clonm%dim(2)%len/) &
          ,VTYPE_DOUBLE)

     !  CALCULATE clonm BELOW together with clatm

     ! ... attributes
     CALL ADD_NCATT(g%clonm, 'long_name', vs='longitude')
     CALL ADD_NCATT(g%clonm, 'units', vs='degrees_east')
     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%clatm%name  = 'lat'
     g%clatm%id    = NULL_VARID
     g%clatm%xtype = NF90_DOUBLE
     ! ... dimensions
     g%clatm%ndims = 2
     ALLOCATE(g%clatm%dim(g%clatm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%clatm%dim(1)%name  = 'lon'
     g%clatm%dim(1)%id    = NULL_DIMID
     g%clatm%dim(1)%len   = ie_in
     g%clatm%dim(1)%fuid  = .false.
     g%clatm%dim(1)%varid = NULL_VARID
     g%clatm%dim(2)%name  = 'lat'
     g%clatm%dim(2)%id    = NULL_DIMID
     g%clatm%dim(2)%len   = je_in
     g%clatm%dim(2)%fuid  = .false.
     g%clatm%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%clatm%dat, g%clatm%ndims &
          , (/g%clatm%dim(1)%len, g%clatm%dim(2)%len/) &
          ,VTYPE_DOUBLE)

     DO i=1, ie_in
        itot = isubpos_in(my_cart_id,1)-nboundlines_in-1
        itmp = istartlon + (itot+i-1) * idlon
        clon = REAL(itmp,dp) / RCF_IN

        DO j = 1 , je_in
           jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
           itmp = istartlat + (jtot+j-1) * idlat
           clat = REAL(itmp,dp) / RCF_IN
           n = (j-1) * g%clonm%dim(1)%len + i
           tmp = rlarot2rla(clat, clon, pollat_in, pollon_in, polgam_in)
           g%clonm%dat%vd(n) = REAL(tmp,dp)
           tmp = phirot2phi(clat, clon, pollat_in, polgam_in)
           g%clatm%dat%vd(n) = REAL(tmp,dp)
        END DO
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%clatm, 'long_name', vs='mid latitude')
     CALL ADD_NCATT(g%clatm, 'units', vs='degrees_north')

     g%ranges(2,1) = RGEMPTY !-90.0_dp
     g%ranges(2,2) = RGEMPTY

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     g%cloni%name  = 'lon_I'
     g%cloni%id    = NULL_VARID
     g%cloni%xtype = NF90_DOUBLE
     ! ... dimensions
     g%cloni%ndims = 2
     ALLOCATE(g%cloni%dim(g%cloni%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%cloni%dim(1)%name  = 'lon_I'
     g%cloni%dim(1)%id    = NULL_DIMID
     g%cloni%dim(1)%len   = ie_in+1
     g%cloni%dim(1)%fuid  = .false.
     g%cloni%dim(1)%varid = NULL_VARID
     g%cloni%dim(2)%name  = 'lat_I'
     g%cloni%dim(2)%id    = NULL_DIMID
     g%cloni%dim(2)%len   = je_in+1
     g%cloni%dim(2)%fuid  = .false.
     g%cloni%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%cloni%dat, g%cloni%ndims &
          , (/g%cloni%dim(1)%len,g%cloni%dim(2)%len/) &
          ,VTYPE_DOUBLE)

     ! CALCULATE CLONI together with clati below

     ! ... attributes
     CALL ADD_NCATT(g%cloni, 'long_name', vs='interface longitude')
     CALL ADD_NCATT(g%cloni, 'units', vs='degrees_east')

     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%clati%name  = 'lat_I'
     g%clati%id    = NULL_VARID
     g%clati%xtype = NF90_DOUBLE
     ! ... dimensions
     g%clati%ndims = 2
     ALLOCATE(g%clati%dim(g%clati%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%clati%dim(1)%name  = 'lon_I'
     g%clati%dim(1)%id    = NULL_DIMID
     g%clati%dim(1)%len   = ie_in+1
     g%clati%dim(1)%fuid  = .false.
     g%clati%dim(1)%varid = NULL_VARID
     g%clati%dim(2)%name  = 'lat_I'
     g%clati%dim(2)%id    = NULL_DIMID
     g%clati%dim(2)%len   = je_in+1
     g%clati%dim(2)%fuid  = .false.
     g%clati%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%clati%dat, g%clati%ndims &
          , (/g%clati%dim(1)%len, g%clati%dim(2)%len/) &
          ,VTYPE_DOUBLE)

     DO i=1, ie_in+1
        itot = isubpos_in(my_cart_id,1)-nboundlines_in-1
        itmp = istartlon + (itot+i) * idlon -NINT(1.5_dp * dlon_in * RCF_IN)
        clon = REAL(itmp,dp) / RCF_IN
        DO j = 1 , je_in+1
           jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
           itmp= istartlat + (jtot+j) * idlat - NINT(1.5_dp * dlat_in * RCF_IN)
           clat = REAL(itmp,dp) / RCF_IN
           n = (j-1) * g%clati%dim(1)%len + i
           tmp = phirot2phi(clat, clon, pollat_in, polgam_in)
           g%clati%dat%vd(n) = REAL(tmp,dp)
           tmp = rlarot2rla(clat, clon, pollat_in, pollon_in, polgam_in)
           g%cloni%dat%vd(n) = REAL(tmp,dp)
        END DO
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%clati, 'long_name', vs='interface latitude')
     CALL ADD_NCATT(g%clati, 'units', vs='degrees_north')

     g%ranges(2,1) = RGEMPTY !-90.0_sp
     g%ranges(2,2) = RGEMPTY

     !**************************************************************************
     ! Rotated  Curvilinear LONGITUDE (MID) ...
     g%rlonm%name  = 'lon'
     g%rlonm%id    = NULL_VARID
     g%rlonm%xtype = NF90_DOUBLE
     g%rlonc = .FALSE.
     ! ... dimensions
     g%rlonm%ndims = 1
     ALLOCATE(g%rlonm%dim(g%rlonm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%rlonm%dim(1)%name  = 'lon'
     g%rlonm%dim(1)%id    = NULL_DIMID
     g%rlonm%dim(1)%len   = ie_in
     g%rlonm%dim(1)%fuid  = .false.

     ! ... data
     CALL INIT_NARRAY(g%rlonm%dat, g%rlonm%ndims &
          , (/g%rlonm%dim(1)%len/),VTYPE_DOUBLE)

     DO i=1, ie_in
        itot = isubpos_in(my_cart_id,1)-nboundlines_in-1
        itmp = istartlon + (itot+i-1) * idlon
        g%rlonm%dat%vd(i) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rlonm, 'long_name', vs='longitude')
     CALL ADD_NCATT(g%rlonm, 'units', vs='degrees_east')
     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%rlatm%name  = 'lat'
     g%rlatm%id    = NULL_VARID
     g%rlatm%xtype = NF90_DOUBLE
     ! ... dimensions
     g%rlatm%ndims = 1
     ALLOCATE(g%rlatm%dim(g%rlatm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%rlatm%dim(1)%name  = 'lat'
     g%rlatm%dim(1)%id    = NULL_DIMID
     g%rlatm%dim(1)%len   = je_in
     g%rlatm%dim(1)%fuid  = .false.
     g%rlatm%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%rlatm%dat, g%rlatm%ndims &
       , (/g%rlatm%dim(1)%len/),VTYPE_DOUBLE)

     DO j=1,je_in
        jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
        itmp = istartlat + (jtot+j-1) * idlat
        g%rlatm%dat%vd(j) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rlatm, 'long_name', vs='mid latitude')
     CALL ADD_NCATT(g%rlatm, 'units', vs='degrees_north')

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     g%rloni%name  = 'lon_I'
     g%rloni%id    = NULL_VARID
     g%rloni%xtype = NF90_DOUBLE
     ! ... dimensions
     g%rloni%ndims = 1
     ALLOCATE(g%rloni%dim(g%rloni%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%rloni%dim(1)%name  = 'lon_I'
     g%rloni%dim(1)%id    = NULL_DIMID
     g%rloni%dim(1)%len   = ie_in+1
     g%rloni%dim(1)%fuid  = .false.
     g%rloni%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%rloni%dat, g%rloni%ndims &
          , (/g%rloni%dim(1)%len/)  ,VTYPE_DOUBLE)

     DO i=1, ie_in+1
       itot = isubpos_in(my_cart_id,1)-nboundlines_in-1
       itmp = istartlon + (itot+i) * idlon - NINT(1.5_dp * dlon_in * RCF_IN)
       g%rloni%dat%vd(i) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rloni, 'long_name', vs='interface longitude')
     CALL ADD_NCATT(g%rloni, 'units', vs='degrees_east')

     ! LATITUDE (MID) ...
     g%rlati%name  = 'lat_I'
     g%rlati%id    = NULL_VARID
     g%rlati%xtype = NF90_DOUBLE
     ! ... dimensions
     g%rlati%ndims = 1
     ALLOCATE(g%rlati%dim(g%rlati%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%rlati%dim(1)%name  = 'lat_I'
     g%rlati%dim(1)%id    = NULL_DIMID
     g%rlati%dim(1)%len   = je_in+1
     g%rlati%dim(1)%fuid  = .false.
     g%rlati%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%rlati%dat, g%rlati%ndims &
          , (/g%rlati%dim(1)%len/) ,VTYPE_DOUBLE)

     DO j = 1 , je_in+1
        jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
        itmp = istartlat + (jtot+j) * idlat - NINT(1.5_dp * dlat_in * RCF_IN)
        g%rlati%dat%vd(j) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rlati, 'long_name', vs='interface latitude')
     CALL ADD_NCATT(g%rlati, 'units', vs='degrees_north')

     ! *************************************************************************
     ! HYBRID-A-COEFFICIENTS (INTERFACES) ...
     g%hyai%name  = 'hyai'
     g%hyai%id    = NULL_VARID
     g%hyai%xtype = NF90_FLOAT
     ! ... dimensions
     g%hyai%ndims = 1
     ALLOCATE(g%hyai%dim(g%hyai%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
     g%hyai%dim(1)%name  = 'ilev'
     g%hyai%dim(1)%id    = NULL_DIMID
     g%hyai%dim(1)%len   = ke_in+1
     g%hyai%dim(1)%fuid  = .false.
     g%hyai%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%hyai%dat, g%hyai%ndims &
          , (/g%hyai%dim(1)%len/),VTYPE_REAL)
     g%hyai%dat%vr(:) = REAL(ak_in,sp)

     ! ... attributes
     CALL ADD_NCATT(g%hyai, 'long_name'                  &
          ,vs='hybrid-A-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hyai, 'units', vs='1')
     g%ranges(3,1) = RGEMPTY
     g%ranges(3,2) = RGEMPTY

     ! ************************************************************************
     ! HYBRID-B-COEFFICIENTS (INTERFACES) ...
     g%hybi%name  = 'hybi'
     g%hybi%id    = NULL_VARID
     g%hybi%xtype = NF90_FLOAT
     ! ... dimensions
     g%hybi%ndims = 1
     ALLOCATE(g%hybi%dim(g%hybi%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
     CALL COPY_NCDIM(g%hybi%dim(1),g%hyai%dim(1))
     ! ... data
     CALL INIT_NARRAY(g%hybi%dat, g%hybi%ndims &
          , (/g%hybi%dim(1)%len/),VTYPE_REAL)

     g%hybi%dat%vr = REAL(bk_in,sp)
     ! ... attributes
     CALL ADD_NCATT(g%hybi, 'long_name'                  &
          ,vs='hybrid-B-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hybi, 'units', vs='1')
     g%ranges(4,1) = 0.0_sp
     g%ranges(4,2) = 1.0_sp

     ! ************************************************************************
     ! SURFACE PRESSURE
     g%ps%name  = 'ps'
     g%ps%id    = NULL_VARID
     g%ps%xtype = NF90_FLOAT
     ! ... dimensions
     g%ps%ndims = 2
     ALLOCATE(g%ps%dim(g%ps%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

     CALL COPY_NCDIM(g%ps%dim(1), g%rlonm%dim(1))
     CALL COPY_NCDIM(g%ps%dim(2), g%rlatm%dim(1))
     ! ... data
     CALL INIT_NARRAY(g%ps%dat, g%ps%ndims      &
          ,(/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
          ,VTYPE_REAL)
        DO i=1, g%ps%dim(1)%len
           DO j=1, g%ps%dim(2)%len
              g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
                   ,(/i,j/)))  = 101325.0_sp
           END DO
        END DO
     ! ... attributes
     CALL ADD_NCATT(g%ps, 'long_name', vs='surface pressure')
     CALL ADD_NCATT(g%ps, 'units', vs='Pa')

     ! *************************************************************************
     ! REFERENCE PRESSURE ...
     g%p0%name  = 'p0'
     g%p0%id    = NULL_VARID
     g%p0%xtype = NF90_FLOAT
     ! ... dimensions
     g%p0%ndims = 0
     ! ... data
     CALL INIT_NARRAY(g%p0%dat, 1, (/ 1 /), VTYPE_REAL)
     g%p0%dat%vr(1) = 1.0_sp
     ! ... attributes
     CALL ADD_NCATT(g%p0, 'long_name', vs='reference pressure')
     CALL ADD_NCATT(g%p0, 'units', vs='Pa')

     ! TIME (MID) ...
     g%timem%name  = 'time'
     g%timem%id    = NULL_VARID
     g%timem%xtype = NF90_DOUBLE
     ! ... dimensions
     g%timem%ndims = 1
     ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
     g%timem%dim(1)%name  = 'time'
     g%timem%dim(1)%id    = NULL_DIMID
     g%timem%dim(1)%len   = 1
     g%timem%dim(1)%fuid  = .true.
     g%timem%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
          ,VTYPE_DOUBLE)
     ! TIME: SECONDS SINCE MODEL START
     CALL time_span_d(dts   &
          , YEAR_START, MONTH_START, DAY_START &
          , HOUR_START, MINUTE_START, SECOND_START  &
          , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
     g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds

     ! ... attributes
     CALL ADD_NCATT(g%timem, 'long_name'                              &
          ,vs='time in seconds since model start')
     WRITE(tunit, &
          '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
          YEAR_START, MONTH_START, DAY_START &
          , HOUR_START, MINUTE_START, SECOND_START

     CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

     ! CALCULATE INTs from MIDs
     CALL COMPLETE_GEOHYBGRID(g)

     ! INTERFACE OK
     ok = .true.

   END SUBROUTINE define_parentin_grid

   ! ------------------------------------------------------------------
   SUBROUTINE DEFINE_INT2COSMO_GRID(g, ok)

     ! MESSY/BMIL
     USE messy_main_grid_def_mem_bi, ONLY: nlev => ke_tot
     USE messy_main_grid_def_bi,  ONLY: hyai, hybi
     USE messy_main_timer,        ONLY: YEAR_START, MONTH_START, DAY_START     &
                                      , HOUR_START, MINUTE_START, SECOND_START &
                                      , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
     ! MESSy
     USE messy_main_timer,         ONLY: time_span_d
     ! NCREGRID MODULES
     USE messy_main_grid_trafo,  ONLY: RGEMPTY, COMPLETE_GEOHYBGRID !&
     USE messy_main_grid_netcdf, ONLY: ERRMSG, COPY_NCDIM       &
                                     , null_dimid, null_varid   &
                                     , vtype_real, vtype_double &
                                     , nf90_float, nf90_double  &
                                     , POSITION, INIT_NARRAY, ADD_NCATT
     USE messy_main_grid,        ONLY: INIT_GEOHYBGRID
     ! INT2COSMO
     USE data_grid_lm,           ONLY: ie2lm, je2lm
     USE data_fields_lm,         ONLY: geolon_lm,  geolat_lm,  lon_lm,  lat_lm &
                                     , geoloni_lm, geolati_lm, loni_lm, lati_lm

     IMPLICIT NONE

     ! I/O
     TYPE(t_geohybgrid), INTENT(OUT) :: g
     LOGICAL           , INTENT(OUT) :: ok

     !LOCAL
     REAL(DP)             :: dts
     INTEGER              :: i,j, n
     INTEGER              :: status
     CHARACTER(LEN=100)   :: tunit

     ! INIT
     CALL INIT_GEOHYBGRID(g)

     g%name = 'INT2COSMOohalo'

     g%file = ' '   ! Filename
     g%t    = 0     ! time step

     ! Curvilinear LONGITUDE (MID) ...
     g%clonm%name  = 'lon'
     g%clonm%id    = NULL_VARID
     g%clonm%xtype = NF90_DOUBLE
     g%clonc = .FALSE.
     ! ... dimensions
     g%clonm%ndims = 2
     ALLOCATE(g%clonm%dim(g%clonm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%clonm%dim(1)%name  = 'lon'
     g%clonm%dim(1)%id    = NULL_DIMID
     g%clonm%dim(1)%len   = ie2lm! - 2*nboundlines
     g%clonm%dim(1)%fuid  = .false.
     g%clonm%dim(1)%varid = NULL_VARID
     g%clonm%dim(2)%name  = 'lat'
     g%clonm%dim(2)%id    = NULL_DIMID
     g%clonm%dim(2)%len   = je2lm! - 2*nboundlines
     g%clonm%dim(2)%fuid  = .false.
     g%clonm%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%clonm%dat, g%clonm%ndims      &
          , (/g%clonm%dim(1)%len,g%clonm%dim(2)%len/) &
          , VTYPE_DOUBLE)

     ! ... attributes
     CALL ADD_NCATT(g%clonm, 'long_name', vs='longitude')
     CALL ADD_NCATT(g%clonm, 'units', vs='degrees_east')
     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%clatm%name  = 'lat'
     g%clatm%id    = NULL_VARID
     g%clatm%xtype = NF90_DOUBLE
     ! ... dimensions
     g%clatm%ndims = 2
     ALLOCATE(g%clatm%dim(g%clatm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%clatm%dim(1)%name  = 'lon'
     g%clatm%dim(1)%id    = NULL_DIMID
     g%clatm%dim(1)%len   = ie2lm
     g%clatm%dim(1)%fuid  = .false.
     g%clatm%dim(1)%varid = NULL_VARID
     g%clatm%dim(2)%name  = 'lat'
     g%clatm%dim(2)%id    = NULL_DIMID
     g%clatm%dim(2)%len   = je2lm
     g%clatm%dim(2)%fuid  = .false.
     g%clatm%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%clatm%dat, g%clatm%ndims       &
          , (/g%clatm%dim(1)%len, g%clatm%dim(2)%len/) &
          , VTYPE_DOUBLE)

     ! ... attributes
     CALL ADD_NCATT(g%clatm, 'long_name', vs='mid latitude')
     CALL ADD_NCATT(g%clatm, 'units', vs='degrees_north')

     g%ranges(2,1) = RGEMPTY !-90.0_dp
     g%ranges(2,2) = RGEMPTY

     ! FILL CLONM / CLATM
     DO i = 1, ie2lm
        DO j = 1, je2lm
           n = (j-1) * g%clatm%dim(1)%len + i
           g%clonm%dat%vd(n) = geolon_lm(i,j)
           g%clatm%dat%vd(n) = geolat_lm(i,j)
        END DO
     END DO

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     g%cloni%name  = 'lon_I'
     g%cloni%id    = NULL_VARID
     g%cloni%xtype = NF90_DOUBLE
     ! ... dimensions
     g%cloni%ndims = 2
     ALLOCATE(g%cloni%dim(g%cloni%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%cloni%dim(1)%name  = 'lon_I'
     g%cloni%dim(1)%id    = NULL_DIMID
     g%cloni%dim(1)%len   = ie2lm +1
     g%cloni%dim(1)%fuid  = .false.
     g%cloni%dim(1)%varid = NULL_VARID
     g%cloni%dim(2)%name  = 'lat_I'
     g%cloni%dim(2)%id    = NULL_DIMID
     g%cloni%dim(2)%len   = je2lm +1
     g%cloni%dim(2)%fuid  = .false.
     g%cloni%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%cloni%dat, g%cloni%ndims      &
          , (/g%cloni%dim(1)%len,g%cloni%dim(2)%len/) &
          , VTYPE_DOUBLE)

     ! ... attributes
     CALL ADD_NCATT(g%cloni, 'long_name', vs='interface longitude')
     CALL ADD_NCATT(g%cloni, 'units', vs='degrees_east')

     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%clati%name  = 'lat_I'
     g%clati%id    = NULL_VARID
     g%clati%xtype = NF90_DOUBLE
     ! ... dimensions
     g%clati%ndims = 2
     ALLOCATE(g%clati%dim(g%clati%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%clati%dim(1)%name  = 'lon_I'
     g%clati%dim(1)%id    = NULL_DIMID
     g%clati%dim(1)%len   = ie2lm + 1
     g%clati%dim(1)%fuid  = .false.
     g%clati%dim(1)%varid = NULL_VARID
     g%clati%dim(2)%name  = 'lat_I'
     g%clati%dim(2)%id    = NULL_DIMID
     g%clati%dim(2)%len   = je2lm +1
     g%clati%dim(2)%fuid  = .false.
     g%clati%dim(2)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%clati%dat, g%clati%ndims       &
          , (/g%clati%dim(1)%len, g%clati%dim(2)%len/) &
          , VTYPE_DOUBLE)
     DO i= 1, ie2lm + 1
        DO j = 1, je2lm + 1
           n = (j-1) * g%cloni%dim(1)%len + i
           g%cloni%dat%vd(n) = geoloni_lm(i,j)
           g%clati%dat%vd(n) = geolati_lm(i,j)
        END DO
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%clati, 'long_name', vs='interface latitude')
     CALL ADD_NCATT(g%clati, 'units', vs='degrees_north')

     g%ranges(2,1) = RGEMPTY !-90.0_sp
     g%ranges(2,2) = RGEMPTY

     !**************************************************************************
     ! Rotated  Curvilinear LONGITUDE (MID) ...
     g%rlonm%name  = 'lon'
     g%rlonm%id    = NULL_VARID
     g%rlonm%xtype = NF90_DOUBLE
     g%rlonc = .FALSE.
     ! ... dimensions
     g%rlonm%ndims = 1
     ALLOCATE(g%rlonm%dim(g%rlonm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%rlonm%dim(1)%name  = 'lon'
     g%rlonm%dim(1)%id    = NULL_DIMID
     g%rlonm%dim(1)%len   = ie2lm !- 2*nboundlines
     g%rlonm%dim(1)%fuid  = .false.

     ! ... data
     CALL INIT_NARRAY(g%rlonm%dat, g%rlonm%ndims &
         , (/g%rlonm%dim(1)%len/),VTYPE_DOUBLE)
     DO i=1, ie2lm
        g%rlonm%dat%vd(i) = lon_lm(i,1)
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rlonm, 'long_name', vs='longitude')
     CALL ADD_NCATT(g%rlonm, 'units', vs='degrees_east')
     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%rlatm%name  = 'lat'
     g%rlatm%id    = NULL_VARID
     g%rlatm%xtype = NF90_DOUBLE
     ! ... dimensions
     g%rlatm%ndims = 1
     ALLOCATE(g%rlatm%dim(g%rlatm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%rlatm%dim(1)%name  = 'lat'
     g%rlatm%dim(1)%id    = NULL_DIMID
     g%rlatm%dim(1)%len   = je2lm !- 2*nboundlines
     g%rlatm%dim(1)%fuid  = .false.
     g%rlatm%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%rlatm%dat, g%rlatm%ndims &
          , (/g%rlatm%dim(1)%len/), VTYPE_DOUBLE)
     DO j=1, je2lm
        g%rlatm%dat%vd(j) = lat_lm(1,j)
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rlatm, 'long_name', vs='mid latitude')
     CALL ADD_NCATT(g%rlatm, 'units', vs='degrees_north')

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     g%rloni%name  = 'lon_I'
     g%rloni%id    = NULL_VARID
     g%rloni%xtype = NF90_DOUBLE
     ! ... dimensions
     g%rloni%ndims = 1
     ALLOCATE(g%rloni%dim(g%rloni%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%rloni%dim(1)%name  = 'lon_I'
     g%rloni%dim(1)%id    = NULL_DIMID
     g%rloni%dim(1)%len   = ie2lm + 1
     g%rloni%dim(1)%fuid  = .false.
     g%rloni%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%rloni%dat, g%rloni%ndims &
          , (/g%rloni%dim(1)%len/)  ,VTYPE_DOUBLE)

     DO i= 1, ie2lm + 1
        g%rloni%dat%vd(i) = loni_lm(i,1)
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rloni, 'long_name', vs='interface longitude')
     CALL ADD_NCATT(g%rloni, 'units', vs='degrees_east')

     ! LATITUDE (MID) ...
     g%rlati%name  = 'lat_I'
     g%rlati%id    = NULL_VARID
     g%rlati%xtype = NF90_DOUBLE
     ! ... dimensions
     g%rlati%ndims = 1
     ALLOCATE(g%rlati%dim(g%rlati%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%rlati%dim(1)%name  = 'lat_I'
     g%rlati%dim(1)%id    = NULL_DIMID
     g%rlati%dim(1)%len   = je2lm + 1
     g%rlati%dim(1)%fuid  = .false.
     g%rlati%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%rlati%dat, g%rlati%ndims &
          , (/g%rlati%dim(1)%len/) ,VTYPE_DOUBLE)
     DO j = 1 , je2lm + 1
        g%rlati%dat%vd(j) = lati_lm(1,j)
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%rlati, 'long_name', vs='interface latitude')
     CALL ADD_NCATT(g%rlati, 'units', vs='degrees_north')

     ! *************************************************************************
     ! HYBRID-A-COEFFICIENTS (INTERFACES) ...
     g%hyai%name  = 'hyai'
     g%hyai%id    = NULL_VARID
     g%hyai%xtype = NF90_FLOAT
     ! ... dimensions
     g%hyai%ndims = 1
     ALLOCATE(g%hyai%dim(g%hyai%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
     g%hyai%dim(1)%name  = 'ilev'
     g%hyai%dim(1)%id    = NULL_DIMID
     g%hyai%dim(1)%len   = nlev+1
     g%hyai%dim(1)%fuid  = .false.
     g%hyai%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%hyai%dat, g%hyai%ndims, (/g%hyai%dim(1)%len/) &
          ,VTYPE_REAL)

     g%hyai%dat%vr(:) = REAL(hyai(1:nlev+1),sp)

     ! ... attributes
     CALL ADD_NCATT(g%hyai, 'long_name'                              &
          ,vs='hybrid-A-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hyai, 'units', vs='1')
     g%ranges(3,1) = RGEMPTY
     g%ranges(3,2) = RGEMPTY

     ! *************************************************************************
     ! HYBRID-B-COEFFICIENTS (INTERFACES) ...
     g%hybi%name  = 'hybi'
     g%hybi%id    = NULL_VARID
     g%hybi%xtype = NF90_FLOAT
     ! ... dimensions
     g%hybi%ndims = 1
     ALLOCATE(g%hybi%dim(g%hybi%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
     CALL COPY_NCDIM(g%hybi%dim(1),g%hyai%dim(1))
     ! ... data
     CALL INIT_NARRAY(g%hybi%dat, g%hybi%ndims, (/g%hybi%dim(1)%len/) &
          ,VTYPE_REAL)
     g%hybi%dat%vr(:) = REAL(hybi(1:nlev+1),sp)

     ! ... attributes
     CALL ADD_NCATT(g%hybi, 'long_name'                              &
          ,vs='hybrid-B-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hybi, 'units', vs='1')
     g%ranges(4,1) = 0.0_sp
     g%ranges(4,2) = 1.0_sp

     ! *************************************************************************
     ! SURFACE PRESSURE
     g%ps%name  = 'ps'
     g%ps%id    = NULL_VARID
     g%ps%xtype = NF90_FLOAT
     ! ... dimensions
     g%ps%ndims = 2
     ALLOCATE(g%ps%dim(g%ps%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

     g%ps%dim(1) = g%rlonm%dim(1)
     g%ps%dim(2) = g%rlatm%dim(1)

     ! ... data
     CALL INIT_NARRAY(g%ps%dat, g%ps%ndims      &
          ,(/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
          ,VTYPE_REAL)

     DO i=1, g%ps%dim(1)%len
        DO j=1, g%ps%dim(2)%len
           g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                ,(/i,j/)))                           &
                = 101325.0_sp
        END DO
     END DO
     ! ... attributes
     CALL ADD_NCATT(g%ps, 'long_name', vs='surface pressure')
     CALL ADD_NCATT(g%ps, 'units', vs='Pa')

     ! *************************************************************************
     ! REFERENCE PRESSURE ...
     g%p0%name  = 'p0'
     g%p0%id    = NULL_VARID
     g%p0%xtype = NF90_FLOAT
     ! ... dimensions
     g%p0%ndims = 0
     ! ... data
     CALL INIT_NARRAY(g%p0%dat, 1, (/ 1 /), VTYPE_REAL)
     g%p0%dat%vr(1) = 1.0_sp
     ! ... attributes
     CALL ADD_NCATT(g%p0, 'long_name', vs='reference pressure')
     CALL ADD_NCATT(g%p0, 'units', vs='Pa')

     ! TIME (MID) ...
     g%timem%name  = 'time'
     g%timem%id    = NULL_VARID
     g%timem%xtype = NF90_DOUBLE
     ! ... dimensions
     g%timem%ndims = 1
     ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
     g%timem%dim(1)%name  = 'time'
     g%timem%dim(1)%id    = NULL_DIMID
     g%timem%dim(1)%len   = 1
     g%timem%dim(1)%fuid  = .true.
     g%timem%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
          ,VTYPE_DOUBLE)
     ! TIME: SECONDS SINCE MODEL START
     CALL time_span_d(dts                           &
          , YEAR_START, MONTH_START, DAY_START      &
          , HOUR_START, MINUTE_START, SECOND_START  &
          , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
     g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds

     ! ... attributes
     CALL ADD_NCATT(g%timem, 'long_name'         &
          ,vs='time in seconds since model start')
     WRITE(tunit, &
          '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
          YEAR_START, MONTH_START, DAY_START &
          , HOUR_START, MINUTE_START, SECOND_START

     CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

     ! CALCULATE INTs from MIDs
     CALL COMPLETE_GEOHYBGRID(g)

     ! INTERFACE OK
     ok = .true.

   END SUBROUTINE DEFINE_INT2COSMO_GRID

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   SUBROUTINE define_rot_parentin_grid(g,ok)
    !
    ! remove halo to avoid double counting of halo area

     ! INT2COSMO
     USE data_int2lm_parallel,    ONLY: isubpos_in, nboundlines_in, my_cart_id
     USE data_grid_in,            ONLY: startlon_in_tot, startlat_in_tot       &
                                      , dlon_in, dlat_in                       &
                                      , ie_in, je_in, ke_in, ak_in, bk_in
     ! MESSy/SMCL
     USE messy_main_timer,        ONLY: YEAR_START, MONTH_START, DAY_START     &
                                      , HOUR_START, MINUTE_START, SECOND_START &
                                      , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

     USE messy_main_constants_mem, ONLY: sp
     USE messy_main_timer,         ONLY: time_span_d

     ! GRID MODULES
     USE messy_main_grid_trafo,    ONLY: RGEMPTY, COMPLETE_GEOHYBGRID
     USE messy_main_grid_netcdf,   ONLY: ERRMSG, COPY_NCDIM       &
                                       , null_dimid, null_varid   &
                                       , VTYPE_REAL, VTYPE_DOUBLE &
                                       , nf90_float, nf90_double  &
                                       , POSITION, INIT_NARRAY, ADD_NCATT
     USE messy_main_grid,          ONLY: INIT_GEOHYBGRID

     IMPLICIT NONE

     ! I/O
     TYPE(t_geohybgrid),  INTENT(OUT) :: g
     LOGICAL            , INTENT(OUT) :: ok

     !LOCAL
     REAL(DP)            :: dts
     INTEGER             :: i,j, itot, jtot
     INTEGER             :: status
     CHARACTER(LEN=100)  :: tunit
     INTEGER             :: itmp
     INTEGER             :: istartlon, istartlat, idlon, idlat

     istartlon = NINT(startlon_in_tot * RCF_IN)
     istartlat = NINT(startlat_in_tot * RCF_IN)
     idlon     = NINT(dlon_in  *RCF_IN)
     idlat     = NINT(dlat_in  *RCF_IN)

     ! INIT
     CALL INIT_GEOHYBGRID(g)

     g%name = 'Parent INGRID'

     g%file = ' '       ! Filename
     g%t    = 0         ! time step

     !**************************************************************************
     ! Rotated  Curvilinear LONGITUDE (MID) ...
     g%lonm%name  = 'lon'
     g%lonm%id    = NULL_VARID
     g%lonm%xtype = NF90_DOUBLE
     g%lonc = .FALSE.
     ! ... dimensions
     g%lonm%ndims = 1
     ALLOCATE(g%lonm%dim(g%lonm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%lonm%dim(1)%name  = 'lon'
     g%lonm%dim(1)%id    = NULL_DIMID
     g%lonm%dim(1)%len   = ie_in
     g%lonm%dim(1)%fuid  = .false.

     ! ... data
     CALL INIT_NARRAY(g%lonm%dat, g%lonm%ndims &
          , (/g%lonm%dim(1)%len/),VTYPE_DOUBLE)
     itot = isubpos_in(my_cart_id,1)-nboundlines_in-1
     DO i=1, ie_in
        itmp = istartlon + (itot+i-1) * idlon
        g%lonm%dat%vd(i) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%lonm, 'long_name', vs='longitude')
     CALL ADD_NCATT(g%lonm, 'units', vs='degrees_east')
     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%latm%name  = 'lat'
     g%latm%id    = NULL_VARID
     g%latm%xtype = NF90_DOUBLE
     ! ... dimensions
     g%latm%ndims = 1
     ALLOCATE(g%latm%dim(g%latm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%latm%dim(1)%name  = 'lat'
     g%latm%dim(1)%id    = NULL_DIMID
     g%latm%dim(1)%len   = je_in
     g%latm%dim(1)%fuid  = .false.
     g%latm%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%latm%dat, g%latm%ndims &
       , (/g%latm%dim(1)%len/),VTYPE_DOUBLE)
     jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
     DO j=1,je_in
        itmp = istartlat + (jtot+j-1) * idlat
        g%latm%dat%vd(j) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%latm, 'long_name', vs='mid latitude')
     CALL ADD_NCATT(g%latm, 'units', vs='degrees_north')

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     g%loni%name  = 'lon_I'
     g%loni%id    = NULL_VARID
     g%loni%xtype = NF90_DOUBLE
     ! ... dimensions
     g%loni%ndims = 1
     ALLOCATE(g%loni%dim(g%loni%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%loni%dim(1)%name  = 'lon_I'
     g%loni%dim(1)%id    = NULL_DIMID
     g%loni%dim(1)%len   = ie_in+1
     g%loni%dim(1)%fuid  = .false.
     g%loni%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%loni%dat, g%loni%ndims &
          , (/g%loni%dim(1)%len/)  ,VTYPE_DOUBLE)

     itot = isubpos_in(my_cart_id,1)-nboundlines_in-1
     DO i=1, ie_in+1
        itmp = istartlon + (itot+i) * idlon - NINT(1.5_dp * dlon_in * RCF_IN)
        g%loni%dat%vd(i) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%loni, 'long_name', vs='interface longitude')
     CALL ADD_NCATT(g%loni, 'units', vs='degrees_east')

     ! LATITUDE (MID) ...
     g%lati%name  = 'lat_I'
     g%lati%id    = NULL_VARID
     g%lati%xtype = NF90_DOUBLE
     ! ... dimensions
     g%lati%ndims = 1
     ALLOCATE(g%lati%dim(g%lati%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%lati%dim(1)%name  = 'lat_I'
     g%lati%dim(1)%id    = NULL_DIMID
     g%lati%dim(1)%len   = je_in+1
     g%lati%dim(1)%fuid  = .false.
     g%lati%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%lati%dat, g%lati%ndims &
          , (/g%lati%dim(1)%len/) ,VTYPE_DOUBLE)
     jtot = isubpos_in(my_cart_id,2)-nboundlines_in-1
     DO j = 1 , je_in+1
        itmp = istartlat + (jtot+j) * idlat - NINT(1.5_dp * dlat_in * RCF_IN)
        g%lati%dat%vd(j) = REAL(itmp,dp) / RCF_IN
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%lati, 'long_name', vs='interface latitude')
     CALL ADD_NCATT(g%lati, 'units', vs='degrees_north')

     ! *************************************************************************
     ! HYBRID-A-COEFFICIENTS (INTERFACES) ...
     g%hyai%name  = 'hyai'
     g%hyai%id    = NULL_VARID
     g%hyai%xtype = NF90_FLOAT
     ! ... dimensions
     g%hyai%ndims = 1
     ALLOCATE(g%hyai%dim(g%hyai%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
     g%hyai%dim(1)%name  = 'ilev'
     g%hyai%dim(1)%id    = NULL_DIMID
     g%hyai%dim(1)%len   = ke_in+1
     g%hyai%dim(1)%fuid  = .false.
     g%hyai%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%hyai%dat, g%hyai%ndims &
          , (/g%hyai%dim(1)%len/),VTYPE_REAL)
     g%hyai%dat%vr(:) = REAL(ak_in,sp)

     ! ... attributes
     CALL ADD_NCATT(g%hyai, 'long_name'                  &
          ,vs='hybrid-A-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hyai, 'units', vs='1')
     g%ranges(3,1) = RGEMPTY
     g%ranges(3,2) = RGEMPTY

     ! ************************************************************************
     ! HYBRID-B-COEFFICIENTS (INTERFACES) ...
     g%hybi%name  = 'hybi'
     g%hybi%id    = NULL_VARID
     g%hybi%xtype = NF90_FLOAT
     ! ... dimensions
     g%hybi%ndims = 1
     ALLOCATE(g%hybi%dim(g%hybi%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
     CALL COPY_NCDIM(g%hybi%dim(1),g%hyai%dim(1))
     ! ... data
     CALL INIT_NARRAY(g%hybi%dat, g%hybi%ndims &
          , (/g%hybi%dim(1)%len/),VTYPE_REAL)

     g%hybi%dat%vr = REAL(bk_in,sp)
     ! ... attributes
     CALL ADD_NCATT(g%hybi, 'long_name'                  &
          ,vs='hybrid-B-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hybi, 'units', vs='1')
     g%ranges(4,1) = 0.0_sp
     g%ranges(4,2) = 1.0_sp

     ! ************************************************************************
     ! SURFACE PRESSURE
     g%ps%name  = 'ps'
     g%ps%id    = NULL_VARID
     g%ps%xtype = NF90_FLOAT
     ! ... dimensions
     g%ps%ndims = 2
     ALLOCATE(g%ps%dim(g%ps%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

     CALL COPY_NCDIM(g%ps%dim(1), g%lonm%dim(1))
     CALL COPY_NCDIM(g%ps%dim(2), g%latm%dim(1))
     ! ... data
     CALL INIT_NARRAY(g%ps%dat, g%ps%ndims      &
          ,(/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
          ,VTYPE_REAL)
        DO i=1, g%ps%dim(1)%len
           DO j=1, g%ps%dim(2)%len
              g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
                   ,(/i,j/)))  = 101325.0_sp
           END DO
        END DO
     ! ... attributes
     CALL ADD_NCATT(g%ps, 'long_name', vs='surface pressure')
     CALL ADD_NCATT(g%ps, 'units', vs='Pa')

     ! *************************************************************************
     ! REFERENCE PRESSURE ...
     g%p0%name  = 'p0'
     g%p0%id    = NULL_VARID
     g%p0%xtype = NF90_FLOAT
     ! ... dimensions
     g%p0%ndims = 0
     ! ... data
     CALL INIT_NARRAY(g%p0%dat, 1, (/ 1 /), VTYPE_REAL)
     g%p0%dat%vr(1) = 1.0_sp
     ! ... attributes
     CALL ADD_NCATT(g%p0, 'long_name', vs='reference pressure')
     CALL ADD_NCATT(g%p0, 'units', vs='Pa')

     ! TIME (MID) ...
     g%timem%name  = 'time'
     g%timem%id    = NULL_VARID
     g%timem%xtype = NF90_DOUBLE
     ! ... dimensions
     g%timem%ndims = 1
     ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
     g%timem%dim(1)%name  = 'time'
     g%timem%dim(1)%id    = NULL_DIMID
     g%timem%dim(1)%len   = 1
     g%timem%dim(1)%fuid  = .true.
     g%timem%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
          ,VTYPE_DOUBLE)
     ! TIME: SECONDS SINCE MODEL START
     CALL time_span_d(dts                           &
          , YEAR_START, MONTH_START, DAY_START      &
          , HOUR_START, MINUTE_START, SECOND_START  &
          , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
     g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds

     ! ... attributes
     CALL ADD_NCATT(g%timem, 'long_name'         &
          ,vs='time in seconds since model start')
     WRITE(tunit, &
          '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
          YEAR_START, MONTH_START, DAY_START &
          , HOUR_START, MINUTE_START, SECOND_START

     CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

     ! CALCULATE INTs from MIDs
     CALL COMPLETE_GEOHYBGRID(g)

     ! INTERFACE OK
     ok = .true.

   END SUBROUTINE define_rot_parentin_grid

   ! ------------------------------------------------------------------
   SUBROUTINE DEFINE_ROTINT2COSMO_GRID(g, ok)

    !
    ! remove halo to avoid double counting of halo area

    ! MESSY/BMIL
     USE messy_main_grid_def_mem_bi, ONLY: nlev => ke_tot
     USE messy_main_grid_def_bi,     ONLY: hyai, hybi

     USE messy_main_timer,        ONLY: YEAR_START, MONTH_START, DAY_START     &
                                      , HOUR_START, MINUTE_START, SECOND_START &
                                      , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

     ! MESSy
     USE messy_main_timer,         ONLY: time_span_d

     ! NCREGRID MODULES
     USE messy_main_grid_trafo,  ONLY: RGEMPTY, COMPLETE_GEOHYBGRID
     USE messy_main_grid_netcdf, ONLY: ERRMSG, COPY_NCDIM       &
                                     , null_dimid, null_varid   &
                                     , vtype_real, vtype_double &
                                     , nf90_float, nf90_double  &
                                     , POSITION, INIT_NARRAY, ADD_NCATT
     USE messy_main_grid,        ONLY: INIT_GEOHYBGRID
     ! INT2COSMO
     USE data_grid_lm,           ONLY: ie2lm, je2lm
     USE data_fields_lm,         ONLY: lon_lm, lat_lm, loni_lm, lati_lm

     IMPLICIT NONE

     ! I/O
     TYPE(t_geohybgrid), INTENT(OUT) :: g
     LOGICAL           , INTENT(OUT) :: ok

     !LOCAL
     REAL(DP)             :: dts
     INTEGER              :: i,j
     INTEGER              :: status
     CHARACTER(LEN=100)   :: tunit

     ! INIT
     CALL INIT_GEOHYBGRID(g)

     g%name = 'INT2COSMOohalo'

     g%file = ' '   ! Filename
     g%t    = 0     ! time step

     !**************************************************************************
     ! Rotated  Curvilinear LONGITUDE (MID) ...
     g%lonm%name  = 'lon'
     g%lonm%id    = NULL_VARID
     g%lonm%xtype = NF90_DOUBLE
     g%lonc = .FALSE.
     ! ... dimensions
     g%lonm%ndims = 1
     ALLOCATE(g%lonm%dim(g%lonm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%lonm%dim(1)%name  = 'lon'
     g%lonm%dim(1)%id    = NULL_DIMID
     g%lonm%dim(1)%len   = ie2lm !- 2*nboundlines
     g%lonm%dim(1)%fuid  = .false.

     ! ... data
     CALL INIT_NARRAY(g%lonm%dat, g%lonm%ndims &
          , (/g%lonm%dim(1)%len/),VTYPE_DOUBLE)
     DO i= 1, ie2lm
        g%lonm%dat%vd(i) = lon_lm(i,1)
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%lonm, 'long_name', vs='longitude')
     CALL ADD_NCATT(g%lonm, 'units', vs='degrees_east')
     g%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
     g%ranges(1,2) = RGEMPTY

     ! LATITUDE (MID) ...
     g%latm%name  = 'lat'
     g%latm%id    = NULL_VARID
     g%latm%xtype = NF90_DOUBLE
     ! ... dimensions
     g%latm%ndims = 1
     ALLOCATE(g%latm%dim(g%latm%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%latm%dim(1)%name  = 'lat'
     g%latm%dim(1)%id    = NULL_DIMID
     g%latm%dim(1)%len   = je2lm
     g%latm%dim(1)%fuid  = .false.
     g%latm%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%latm%dat, g%latm%ndims &
          , (/g%latm%dim(1)%len/), VTYPE_DOUBLE)

     DO j=1, je2lm
        g%latm%dat%vd(j) = lat_lm(1,j)
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%latm, 'long_name', vs='mid latitude')
     CALL ADD_NCATT(g%latm, 'units', vs='degrees_north')

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     g%loni%name  = 'lon_I'
     g%loni%id    = NULL_VARID
     g%loni%xtype = NF90_DOUBLE
     ! ... dimensions
     g%loni%ndims = 1
     ALLOCATE(g%loni%dim(g%loni%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
     g%loni%dim(1)%name  = 'lon_I'
     g%loni%dim(1)%id    = NULL_DIMID
     g%loni%dim(1)%len   = ie2lm + 1
     g%loni%dim(1)%fuid  = .false.
     g%loni%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%loni%dat, g%loni%ndims &
          , (/g%loni%dim(1)%len/)  ,VTYPE_DOUBLE)
     DO i=1, ie2lm +1
        g%loni%dat%vd(i) =  loni_lm(i,1)
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%loni, 'long_name', vs='interface longitude')
     CALL ADD_NCATT(g%loni, 'units', vs='degrees_east')

     ! LATITUDE (MID) ...
     g%lati%name  = 'lat_I'
     g%lati%id    = NULL_VARID
     g%lati%xtype = NF90_DOUBLE
     ! ... dimensions
     g%lati%ndims = 1
     ALLOCATE(g%lati%dim(g%lati%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
     g%lati%dim(1)%name  = 'lat_I'
     g%lati%dim(1)%id    = NULL_DIMID
     g%lati%dim(1)%len   = je2lm +  1
     g%lati%dim(1)%fuid  = .false.
     g%lati%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%lati%dat, g%lati%ndims &
          , (/g%lati%dim(1)%len/) ,VTYPE_DOUBLE)
     DO j = 1 , je2lm + 1
        g%lati%dat%vd(j) = lati_lm(1,j)
     END DO

     ! ... attributes
     CALL ADD_NCATT(g%lati, 'long_name', vs='interface latitude')
     CALL ADD_NCATT(g%lati, 'units', vs='degrees_north')

     ! *************************************************************************
     ! HYBRID-A-COEFFICIENTS (INTERFACES) ...
     g%hyai%name  = 'hyai'
     g%hyai%id    = NULL_VARID
     g%hyai%xtype = NF90_FLOAT
     ! ... dimensions
     g%hyai%ndims = 1
     ALLOCATE(g%hyai%dim(g%hyai%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
     g%hyai%dim(1)%name  = 'ilev'
     g%hyai%dim(1)%id    = NULL_DIMID
     g%hyai%dim(1)%len   = nlev+1
     g%hyai%dim(1)%fuid  = .false.
     g%hyai%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%hyai%dat, g%hyai%ndims, (/g%hyai%dim(1)%len/) &
          ,VTYPE_REAL)
     g%hyai%dat%vr(:) = REAL(hyai(1:nlev+1),sp)

     ! ... attributes
     CALL ADD_NCATT(g%hyai, 'long_name'                  &
          ,vs='hybrid-A-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hyai, 'units', vs='1')
     g%ranges(3,1) = RGEMPTY
     g%ranges(3,2) = RGEMPTY

     ! *************************************************************************
     ! HYBRID-B-COEFFICIENTS (INTERFACES) ...
     g%hybi%name  = 'hybi'
     g%hybi%id    = NULL_VARID
     g%hybi%xtype = NF90_FLOAT
     ! ... dimensions
     g%hybi%ndims = 1
     ALLOCATE(g%hybi%dim(g%hybi%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
     CALL COPY_NCDIM(g%hybi%dim(1),g%hyai%dim(1))
     ! ... data
     CALL INIT_NARRAY(g%hybi%dat, g%hybi%ndims, (/g%hybi%dim(1)%len/) &
          ,VTYPE_REAL)
     g%hybi%dat%vr(:) = REAL(hybi(1:nlev+1),sp)
     ! ... attributes
     CALL ADD_NCATT(g%hybi, 'long_name'                  &
          ,vs='hybrid-B-coefficients at layer interfaces')
     CALL ADD_NCATT(g%hybi, 'units', vs='1')
     g%ranges(4,1) = 0.0_sp
     g%ranges(4,2) = 1.0_sp

     ! *************************************************************************
     ! SURFACE PRESSURE
     g%ps%name  = 'ps'
     g%ps%id    = NULL_VARID
     g%ps%xtype = NF90_FLOAT
     ! ... dimensions
     g%ps%ndims = 2
     ALLOCATE(g%ps%dim(g%ps%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)

     g%ps%dim(1) = g%lonm%dim(1)
     g%ps%dim(2) = g%latm%dim(1)

     ! ... data
     CALL INIT_NARRAY(g%ps%dat, g%ps%ndims      &
          ,(/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
          ,VTYPE_REAL)
     DO i=1, g%ps%dim(1)%len
        DO j=1, g%ps%dim(2)%len
           g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                ,(/i,j/)))                                            &
                !           = aps(i,j)
                = 101325.0_sp
        END DO
     END DO
     ! ... attributes
     CALL ADD_NCATT(g%ps, 'long_name', vs='surface pressure')
     CALL ADD_NCATT(g%ps, 'units',     vs='Pa')

     ! *************************************************************************
     ! REFERENCE PRESSURE ...
     g%p0%name  = 'p0'
     g%p0%id    = NULL_VARID
     g%p0%xtype = NF90_FLOAT
     ! ... dimensions
     g%p0%ndims = 0
     ! ... data
     CALL INIT_NARRAY(g%p0%dat, 1, (/ 1 /), VTYPE_REAL)
     g%p0%dat%vr(1) = 1.0_sp
     ! ... attributes
     CALL ADD_NCATT(g%p0, 'long_name', vs='reference pressure')
     CALL ADD_NCATT(g%p0, 'units', vs='Pa')

     ! TIME (MID) ...
     g%timem%name  = 'time'
     g%timem%id    = NULL_VARID
     g%timem%xtype = NF90_DOUBLE
     ! ... dimensions
     g%timem%ndims = 1
     ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
     CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
     g%timem%dim(1)%name  = 'time'
     g%timem%dim(1)%id    = NULL_DIMID
     g%timem%dim(1)%len   = 1
     g%timem%dim(1)%fuid  = .true.
     g%timem%dim(1)%varid = NULL_VARID
     ! ... data
     CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
          ,VTYPE_DOUBLE)
     ! TIME: SECONDS SINCE MODEL START
     CALL time_span_d(dts   &
          , YEAR_START, MONTH_START, DAY_START &
          , HOUR_START, MINUTE_START, SECOND_START  &
          , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
     g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds

     ! ... attributes
     CALL ADD_NCATT(g%timem, 'long_name'                              &
          ,vs='time in seconds since model start')
     WRITE(tunit, &
          '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
          YEAR_START, MONTH_START, DAY_START &
          , HOUR_START, MINUTE_START, SECOND_START

     CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

     ! CALCULATE INTs from MIDs
     CALL COMPLETE_GEOHYBGRID(g)

     ! INTERFACE OK
     ok = .true.

   END SUBROUTINE DEFINE_ROTINT2COSMO_GRID

  END SUBROUTINE CALC_forward_weights
  ! ---------------------------------------------------------------------------

   ! ------------------------------------------------------------------
   SUBROUTINE REDEFINE_STAGGERED_UGRID(g, ok)

     ! MESSY/BMIL
     USE messy_main_mpi_bi,      ONLY: my_cart_id
     ! MESSy/SMCL
     ! GRID MODULES
     USE messy_main_grid_trafo,  ONLY: rlarot2rla, phirot2phi
     USE messy_main_grid_netcdf, ONLY: VTYPE_DOUBLE     &
                                     , INIT_NARRAY, QDEF_NCVAR
     ! INT2COSMO
     USE data_grid_lm,           ONLY: ie2lm, je2lm, startlon_tot,startlat_tot &
                                     , pollon, pollat, polgam, dlon, dlat
     USE data_fields_lm,         ONLY: lon_lm, loni_lm
     USE data_int2lm_parallel,   ONLY: isubpos, nboundlines

     IMPLICIT NONE

     ! I/O
     TYPE(t_geohybgrid), INTENT(INOUT) :: g
     LOGICAL           , INTENT(OUT) :: ok

     !LOCAL
     INTEGER              :: i,j, n, itot, jtot
     REAL(dp)             :: tmplon, tmplat
     INTEGER              :: itmp
     INTEGER              :: istartlon, istartlat, idlon, idlat

     istartlon = NINT(startlon_tot * RCF)
     istartlat = NINT(startlat_tot * RCF)
     idlon     = NINT(dlon         * RCF)
     idlat     = NINT(dlat         * RCF)

     ! INIT
     IF (QDEF_NCVAR(g%clonm)) THEN

        IF (.NOT. QDEF_NCVAR(g%clatm)) THEN
           ok = .FALSE.
           RETURN
        END IF

        ! ...reset data to staggered grid , ie mids = interfaces
        CALL INIT_NARRAY(g%clonm%dat, g%clonm%ndims      &
             , (/g%clonm%dim(1)%len,g%clonm%dim(2)%len/) &
             , VTYPE_DOUBLE)
        CALL INIT_NARRAY(g%clatm%dat, g%clatm%ndims      &
             , (/g%clatm%dim(1)%len,g%clatm%dim(2)%len/) &
             , VTYPE_DOUBLE)

        ! Note: the geographical midpoints for the staggered longitude
        !       are calculated by latitude ad midpoints and longitude at
        !       interfaces
        itot = isubpos(my_cart_id,1) - nboundlines - 1
        jtot = isubpos(my_cart_id,2) - nboundlines - 1
        DO i = 1, ie2lm
           itmp  = istartlon + (itot+i-2) * idlon - NINT(0.5_dp * dlon * RCF)
           tmplon = REAL(itmp,dp) / RCF
           DO j = 1, je2lm
              itmp = istartlat + (jtot+j-2) * idlat
              tmplat = REAL(itmp,dp) / RCF
              n = (j-1) * g%clonm%dim(1)%len + i
              g%clonm%dat%vd(n) =  &
                   rlarot2rla(tmplat , tmplon, pollat, pollon, polgam)
              g%clatm%dat%vd(n) =  &
                   phirot2phi(tmplat , tmplon, pollat, polgam)
           END DO
        END DO
     ENDIF

     IF (QDEF_NCVAR(g%cloni)) THEN

        ! ...reset data to staggered  u grid
        ! ----------------------------------------------------------
        ! define interfaces:
        ! Curvilinear LONGITUDE (INTERFACES) ...

        IF (.NOT. QDEF_NCVAR(g%clati)) THEN
           ok = .FALSE.
           RETURN
        END IF

        CALL INIT_NARRAY(g%cloni%dat, g%cloni%ndims      &
             , (/g%cloni%dim(1)%len,g%cloni%dim(2)%len/) &
             , VTYPE_DOUBLE)
        CALL INIT_NARRAY(g%clati%dat, g%clati%ndims      &
             , (/g%clati%dim(1)%len,g%clati%dim(2)%len/) &
             , VTYPE_DOUBLE)

        ! Note: the geographical interface points for the staggered longitude
        !       are calculated by interfaces of latitude and midpoints (shifted
        !       by -1 for longitudes)
        itot = isubpos(my_cart_id,1) - nboundlines - 1
        jtot = isubpos(my_cart_id,2) - nboundlines - 1
        DO i = 1, ie2lm + 1
           ! "interface" of interfacer => " -dlon"
           itmp = istartlon + (itot+i-2) * idlon - idlon
           tmplon = REAL(itmp,dp) / RCF
           DO j = 1, je2lm + 1
              itmp = istartlat + (jtot+j-2) * idlat - NINT(0.5_dp * dlat * RCF)
              tmplat = REAL(itmp,dp) / RCF
              n = (j-1) * g%cloni%dim(1)%len + i
              g%cloni%dat%vd(n) = &
                   rlarot2rla(tmplat , tmplon, pollat, pollon, polgam)
              g%clati%dat%vd(n) = &
                   phirot2phi(tmplat , tmplon, pollat, polgam)
           END DO
        END DO

     END IF


     !**************************************************************************
     ! Rotated  Curvilinear LONGITUDE (MID) ...
     IF (QDEF_NCVAR(g%rlonm)) THEN
        ! ... reset data
        CALL INIT_NARRAY(g%rlonm%dat, g%rlonm%ndims &
             , (/g%rlonm%dim(1)%len/),VTYPE_DOUBLE)
        DO i=1, ie2lm
           g%rlonm%dat%vd(i) = loni_lm(i,1)
        END DO
     END IF

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     IF (QDEF_NCVAR(g%rloni)) THEN
        ! ... reset data
        CALL INIT_NARRAY(g%rloni%dat, g%rloni%ndims &
             , (/g%rloni%dim(1)%len/)  ,VTYPE_DOUBLE)

        DO i= 2, ie2lm + 1
           g%rloni%dat%vd(i) = lon_lm(i-1,1)
        END DO
        ! i= 1
        g%rloni%dat%vd(1) =  g%rloni%dat%vd(2) - dlon
     END IF
     ! INTERFACE OK
     ok = .true.

   END SUBROUTINE REDEFINE_STAGGERED_UGRID
   ! ------------------------------------------------------------------

   ! ------------------------------------------------------------------
   SUBROUTINE REDEFINE_STAGGERED_VGRID(g, ok)

     ! MESSY/BMIL
     USE messy_main_mpi_bi,      ONLY: my_cart_id
     ! MESSy/SMCL
     ! GRID MODULES
     USE messy_main_grid_trafo,  ONLY: rlarot2rla, phirot2phi
     USE messy_main_grid_netcdf, ONLY: VTYPE_DOUBLE     &
                                     , INIT_NARRAY, QDEF_NCVAR
     ! INT2COSMO
     USE data_grid_lm,           ONLY: ie2lm, je2lm, startlon_tot,startlat_tot &
                                     , pollon, pollat, polgam, dlon, dlat
     USE data_fields_lm,         ONLY: lat_lm, lati_lm
     USE data_int2lm_parallel,   ONLY: isubpos, nboundlines

     IMPLICIT NONE

     ! I/O
     TYPE(t_geohybgrid), INTENT(INOUT) :: g
     LOGICAL           , INTENT(OUT)   :: ok

     !LOCAL
     INTEGER              :: i,j, n, itot, jtot
     REAL(dp)             :: tmplon, tmplat
     INTEGER              :: itmp
     INTEGER              :: istartlon, istartlat, idlon, idlat

     istartlon = NINT(startlon_tot * RCF)
     istartlat = NINT(startlat_tot * RCF)
     idlon     = NINT(dlon         * RCF)
     idlat     = NINT(dlat         * RCF)

     ! INIT
     IF (QDEF_NCVAR(g%clatm)) THEN

        IF (.NOT. QDEF_NCVAR(g%clonm)) THEN
           ok = .FALSE.
           RETURN
        END IF

        ! ...reset data to staggered grid , ie mids = interfaces
        CALL INIT_NARRAY(g%clonm%dat, g%clonm%ndims      &
             , (/g%clonm%dim(1)%len,g%clonm%dim(2)%len/) &
             , VTYPE_DOUBLE)
        CALL INIT_NARRAY(g%clatm%dat, g%clatm%ndims      &
             , (/g%clatm%dim(1)%len,g%clatm%dim(2)%len/) &
             , VTYPE_DOUBLE)

        ! Note: the geographical midpoints for the staggered latitude
        !       are calculated by longitude ad midpoints and latitude at
        !       interfaces
        itot = isubpos(my_cart_id,1) - nboundlines - 1
        jtot = isubpos(my_cart_id,2) - nboundlines - 1
        DO i = 1, ie2lm
           itmp = istartlon + (itot+i-2) * idlon
           tmplon = REAL(itmp,dp) / RCF
           DO j = 1, je2lm
              itmp = istartlat + (jtot+j-2) * idlat - NINT(0.5_dp * dlat *RCF)
              tmplat = REAL(itmp,dp) / RCF
              n = (j-1) * g%clonm%dim(1)%len + i
              g%clonm%dat%vd(n) =  &
                   rlarot2rla(tmplat , tmplon, pollat, pollon, polgam)
              g%clatm%dat%vd(n) =  &
                   phirot2phi(tmplat , tmplon, pollat, polgam)
           END DO
        END DO
     ENDIF

     IF (QDEF_NCVAR(g%clati)) THEN

        ! ...reset data to staggered  u grid
        ! ----------------------------------------------------------
        ! define interfaces:
        ! Curvilinear LONGITUDE (INTERFACES) ...

        IF (.NOT. QDEF_NCVAR(g%cloni)) THEN
           ok = .FALSE.
           RETURN
        END IF

        CALL INIT_NARRAY(g%cloni%dat, g%cloni%ndims      &
             , (/g%cloni%dim(1)%len,g%cloni%dim(2)%len/) &
             , VTYPE_DOUBLE)
        CALL INIT_NARRAY(g%clati%dat, g%clati%ndims      &
             , (/g%clati%dim(1)%len,g%clati%dim(2)%len/) &
             , VTYPE_DOUBLE)

        ! Note: the geographical interface points for the staggered latitude
        !       are calculated by interfaces of longitude and midpoints (shifted
        !       by -1 for latitudes)
        itot = isubpos(my_cart_id,1) - nboundlines - 1
        jtot = isubpos(my_cart_id,2) - nboundlines - 1
        DO i = 1, ie2lm + 1
           itmp = istartlon + (itot+i-2) * idlon - NINT(0.5_dp * dlon * RCF)
           tmplon = REAL(itmp,dp) / RCF
           DO j = 1, je2lm + 1
              n = (j-1) * g%cloni%dim(1)%len + i
              ! "interface" of interfacer => " -dlon"
              itmp = istartlat + (jtot+j-2) * idlat - idlat
              tmplat = REAL(itmp,dp) / RCF
              g%cloni%dat%vd(n) =  &
                   rlarot2rla(tmplat , tmplon, pollat, pollon, polgam)
              g%clati%dat%vd(n) =  &
                   phirot2phi(tmplat , tmplon, pollat, polgam)
           END DO
        END DO

     END IF


     !**************************************************************************
     ! Rotated  Curvilinear LONGITUDE (MID) ...
     IF (QDEF_NCVAR(g%rlatm)) THEN
        ! ... reset data
        CALL INIT_NARRAY(g%rlatm%dat, g%rlatm%ndims &
             , (/g%rlatm%dim(1)%len/),VTYPE_DOUBLE)
        DO j=1, je2lm
           g%rlatm%dat%vd(j) = lati_lm(1,j)
        END DO
     END IF

     ! ----------------------------------------------------------
     ! define interfaces:
     ! Curvilinear LONGITUDE (INTERFACES) ...
     IF (QDEF_NCVAR(g%rlati)) THEN
        ! ... reset data
        CALL INIT_NARRAY(g%rlati%dat, g%rlati%ndims &
             , (/g%rlati%dim(1)%len/)  ,VTYPE_DOUBLE)

        DO j= 2, je2lm + 1
           g%rlati%dat%vd(j) = lat_lm(1,j-1)
        END DO
        ! i= 1
        g%rlati%dat%vd(1) =  g%rlati%dat%vd(2) - dlat
     END IF
     ! INTERFACE OK
     ok = .true.

   END SUBROUTINE REDEFINE_STAGGERED_VGRID
   ! ------------------------------------------------------------------
   ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! ---------------------------------------------------------------------------
  SUBROUTINE vert_ncinterpol(name, var, ps_s, ps_d, grid_s, grid_d, grid_o)

    USE data_grid_in,               ONLY: ke_in
    USE messy_main_grid,            ONLY: GRID_ERROR &
                                        , COPY_GEOHYBGRID          &
                                        , INIT_GEOHYBGRID
    USE messy_main_grid_netcdf,     ONLY: INIT_NCVAR, COPY_NCVAR, QDEF_NCVAR
    USE messy_main_grid_tools,      ONLY: RGTOOL_CONVERT         &
                                        , RGTOOL_CONVERT_DAT2VAR &
                                        , SET_SURFACE_PRESSURE
    USE messy_main_grid_trafo,      ONLY: COMPLETE_GEOHYBGRID, RG_INT
    USE messy_main_grid_trafo_nrgd, ONLY: REGRID_CONTROL

    IMPLICIT NONE

    CHARACTER(LEN=*),             INTENT(IN)  :: name
    REAL(dp), DIMENSION(:,:,:,:), POINTER     :: var
    REAL(dp), DIMENSION(:,:,:,:), POINTER     :: ps_s
    REAL(dp), DIMENSION(:,:,:,:), POINTER     :: ps_d
    ! full 3d source grid of vertical interpolation
    TYPE(t_geohybgrid), INTENT(IN)            :: grid_s
    ! source grid as defined after horizontal interpolation (grid of var_in)
    ! destination/target grid of vertical interpolation
    TYPE(t_geohybgrid), INTENT(IN)            :: grid_d
    ! output grid after vertical interpolation (grid of var_out)
    TYPE(t_geohybgrid), INTENT(OUT), OPTIONAL :: grid_o

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr ='vert_ncinterpol'

    TYPE(t_geohybgrid)                    :: intzogrid
    TYPE(t_geohybgrid)                    :: intpsgrid

    TYPE(t_ncvar), DIMENSION(:),  POINTER :: sovar => NULL()
    TYPE(t_ncvar), DIMENSION(:),  POINTER :: ovar  => NULL()
    TYPE(t_ncvar)                         :: varo

    REAL(dp), DIMENSION(:,:,:,:), POINTER :: dat   => NULL()
    INTEGER,       DIMENSION(1)           :: RGT
    INTEGER                               :: ix
    INTEGER                               :: status

    ! Note:
    ! a) both variables / grids need to have the same horizontal rectangular
    !    grid, which might be differently defined by the grids
    !    Thus the curvilinear and rotated definitions are deleted and the
    !    rectangular lon/lat struct elements are filled.
    !
    ! b) SURFACE PRESSURE needs to be available for both grid
    !    on the same horizontal (target) grid

    ! --------------------
    ! make input grid
    ! --------------------

    CALL COPY_GEOHYBGRID(intpsgrid, grid_d)
    ! -- correct also lon/lat definition in ps field

    IF ((QDEF_NCVAR(grid_d%lonm) .AND. QDEF_NCVAR(grid_d%latm))  .OR. &
         (QDEF_NCVAR(grid_d%loni) .AND. QDEF_NCVAR(grid_d%lati))) THEN
       ! NOTHING TO DO
    ELSE IF ((QDEF_NCVAR(grid_d%rlonm) .AND. QDEF_NCVAR(grid_d%rlatm))  .OR. &
         (QDEF_NCVAR(grid_d%rloni) .AND. QDEF_NCVAR(grid_d%rlati))) THEN
       ! -- correct also lon/lat definition
       intpsgrid%lonc=.FALSE.
       CALL COPY_NCVAR(intpsgrid%lonm, grid_d%rlonm)
       CALL COPY_NCVAR(intpsgrid%loni, grid_d%rloni)
       CALL COPY_NCVAR(intpsgrid%latm, grid_d%rlatm)
       CALL COPY_NCVAR(intpsgrid%lati, grid_d%rlati)
    ELSE
       CALL error_bi('ncregrid vertical integration requires rectangular grid' &
            , substr)
    ENDIF

    ! remove curvilinear / rotated information
    ! vertical interpolation proceeds on target grid
    CALL INIT_NCVAR(intpsgrid%clonm)
    CALL INIT_NCVAR(intpsgrid%cloni)
    CALL INIT_NCVAR(intpsgrid%clatm)
    CALL INIT_NCVAR(intpsgrid%clati)

    CALL INIT_NCVAR(intpsgrid%rlonm)
    CALL INIT_NCVAR(intpsgrid%rloni)
    CALL INIT_NCVAR(intpsgrid%rlatm)
    CALL INIT_NCVAR(intpsgrid%rlati)

    ! add vertical grid (deleted for horizontal interpolation)
    CALL COPY_NCVAR(intpsgrid%hyam, grid_s%hyam)
    CALL COPY_NCVAR(intpsgrid%hybm, grid_s%hybm)
    CALL COPY_NCVAR(intpsgrid%hyai, grid_s%hyai)
    CALL COPY_NCVAR(intpsgrid%hybi, grid_s%hybi)
    intpsgrid%ranges(4,1) = grid_s%ranges(4,1)
    intpsgrid%ranges(4,2) = grid_s%ranges(4,2)

    ! construct intermediate ps grid
    CALL SET_SURFACE_PRESSURE(status, intpsgrid, ps_s(:,:,1,1))
    IF (status /= 0) CALL error_bi(grid_error(status), substr)
    CALL COMPLETE_GEOHYBGRID(intpsgrid)

    ! --------------------
    ! make output grid
    ! --------------------
    CALL COPY_GEOHYBGRID(intzogrid, grid_d)
    ! assume congruency of grids (grid_sh and grid_s)
    ! take rotated coordinates of curvi-linear grids
    ! instead of intermediate grid. NCREGRID can only
    ! handle lat-lon-grids

    IF ((QDEF_NCVAR(grid_d%rlonm) .AND. QDEF_NCVAR(grid_d%rlatm))  .OR. &
         (QDEF_NCVAR(grid_d%rloni) .AND. QDEF_NCVAR(grid_d%rlati))) THEN
       ! set non-rotated lat-lon grid
       intzogrid%lonc=.FALSE.
       CALL COPY_NCVAR(intzogrid%lonm, grid_d%rlonm)
       CALL COPY_NCVAR(intzogrid%loni, grid_d%rloni)
       CALL COPY_NCVAR(intzogrid%latm, grid_d%rlatm)
       CALL COPY_NCVAR(intzogrid%lati, grid_d%rlati)
    ELSE
       ! ??? required ??
       ! -- correct also lon/lat definition
       CALL COPY_NCVAR(intzogrid%lonm, grid_d%lonm)
       CALL COPY_NCVAR(intzogrid%loni, grid_d%loni)
       CALL COPY_NCVAR(intzogrid%latm, grid_d%latm)
       CALL COPY_NCVAR(intzogrid%lati, grid_d%lati)
    ENDIF

    ! remove curvilinear / rotated information
    ! vertical interpolation proceeds on target grid
    CALL INIT_NCVAR(intzogrid%clonm)
    CALL INIT_NCVAR(intzogrid%cloni)
    CALL INIT_NCVAR(intzogrid%clatm)
    CALL INIT_NCVAR(intzogrid%clati)

    CALL INIT_NCVAR(intzogrid%rlonm)
    CALL INIT_NCVAR(intzogrid%rloni)
    CALL INIT_NCVAR(intzogrid%rlatm)
    CALL INIT_NCVAR(intzogrid%rlati)

    ! construct destination surface pressure (ps) for destination grid
    CALL SET_SURFACE_PRESSURE(status, intzogrid, ps_d(:,:,1,1))
    IF (status /= 0) CALL error_bi(grid_error(status), substr)

    ! ----------------------------
    ! make vertical interpolation
    ! ----------------------------
    ALLOCATE(dat(SIZE(var,1), SIZE(var,2), ke_in, 1))
    dat(:,:,:,1) = var(:,:,1:ke_in,1)

    ALLOCATE(sovar(1))
    ! construct grid with vertical input grid, but horizontal
    CALL RGTOOL_CONVERT_DAT2VAR(sovar(1), dat, name, intpsgrid, 'xyzn')
    DEALLOCATE(dat)

    RGT(1) = RG_INT
    CALL REGRID_CONTROL(intpsgrid, intzogrid, sovar, ovar, RGT, .TRUE. &
         , lrgx=.FALSE., lrgy=.FALSE., lrgz=.TRUE., lpresaxis=.FALSE. &
         , lwork=.TRUE., lstatout=.FALSE.)

    CALL COPY_NCVAR(varo, ovar(1))   ! ... RETURN VARIABLES
    CALL INIT_NCVAR(ovar(1))
    DEALLOCATE(ovar)
    NULLIFY(ovar)

    ! CONVERT BACK
    CALL RGTOOL_CONVERT(varo, dat, intzogrid, 'xyzn')
    IF (PRESENT(grid_o)) CALL COPY_GEOHYBGRID(grid_o, intzogrid)

    var = 0._dp
    var(1:SIZE(dat,1),1:SIZE(dat,2),1:SIZE(dat,3),1:SIZE(dat,4)) = dat

    CALL INIT_GEOHYBGRID(intpsgrid)
    CALL INIT_GEOHYBGRID(intzogrid)
    CALL INIT_NCVAR(varo)
    DO ix= 1, SIZE(sovar)
       CALL INIT_NCVAR(sovar(ix))
    END DO
    DEALLOCATE(sovar)
    NULLIFY(sovar)
    DEALLOCATE(dat)
    NULLIFY(dat)

  END SUBROUTINE vert_ncinterpol
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE calc_p_hsurf_full

    ! MESSY/BMIL
    USE messy_main_grid_def_mem_bi,      ONLY: BASEGRID_ID
    ! MESSy/SMCL
    USE messy_main_grid,                 ONLY: LOCATE_GEOHYBGRID               &
                                             , GRID_ERROR, COPY_GEOHYBGRID     &
                                             , INIT_GEOHYBGRID
    USE messy_main_grid_netcdf,          ONLY: MAIN_GRID_SET_MESSAGEMODE       &
                                             , MSGMODE_E, INIT_NCVAR
    USE messy_main_grid_tools,           ONLY: RGTOOL_CONVERT                  &
                                             , RGTOOL_CONVERT_DAT2VAR
    USE messy_main_grid_trafo,           ONLY: RG_INT, SWITCH_GEOHYBGRID
    USE messy_main_grid_trafo_scrp,      ONLY: CALC_SCRIPDATA, SCRIP_CONTROL   &
                                             , CALC_SCRIP_WEIGHTS

    USE messy_main_grid_trafo_scrp_base, ONLY: norm_opt_dstarea
    USE messy_main_channel,              ONLY: get_channel_object
    USE messy_main_channel_mem,          ONLY: dom_current
    USE messy_main_channel_error_bi,     ONLY: channel_halt

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'calc_p_hsurf_full'
    INTEGER                     :: status
    INTEGER, DIMENSION(1)       :: RGTYPE = RG_INT

    ! used for interpolation of horizontal weight function
    TYPE(t_ncvar), DIMENSION(:),  POINTER :: vari  => NULL()
    TYPE(t_ncvar), DIMENSION(:),  POINTER :: sovar => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: dat   => NULL()
    INTEGER                               :: ix
    TYPE(t_geohybgrid)                    :: basegrid
    TYPE(t_geohybgrid)                    :: basehgrid
    TYPE(t_scrip_data), POINTER           :: hsurf_PSD   => NULL()
    INTEGER                               :: hsurf_SD_ID = -99
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: hsurf_loc   => NULL()

    ! get basemodel grid
    CALL INIT_GEOHYBGRID(basegrid)
    CALL LOCATE_GEOHYBGRID (status, BASEGRID_ID(dom_current), grid= basegrid)
    IF (status /= 0) CALL ERROR_BI(grid_error(status), substr)

    CALL COPY_GEOHYBGRID(basehgrid,basegrid)
    CALL SWITCH_GEOHYBGRID(basehgrid, .TRUE., .TRUE., .FALSE.)

    CALL CALC_SCRIPDATA(status, basegrid, pgrid, RGTYPE, hsurf_SD_ID &
         , PSD=hsurf_PSD, norm_opt_in=norm_opt_dstarea, map_type_in=1)

    IF (status /= 0) THEN
       IF (status /= 01 )  RETURN
    ELSE
       ! CALCULATE WEIGHTS
       CALL CALC_SCRIP_WEIGHTS(status, hsurf_PSD)
    ENDIF

    CALL MAIN_GRID_SET_MESSAGEMODE(MSGMODE_E)
    !CALL MAIN_GRID_SET_MESSAGEMODE()

    ALLOCATE(vari(1))

    CALL get_channel_object(status,'COSMO_ORI','HSURF', p4 = hsurf_loc)
    CALL channel_halt(substr, status)

    dat => hsurf_loc(:,:,:,:)

    CALL RGTOOL_CONVERT_DAT2VAR(vari(1), dat, 'hsurf', basehgrid, 'xyzn')
    NULLIFY(dat)

    CALL SCRIP_CONTROL(status, hsurf_SD_ID, basehgrid, phgrid &
         , RGTYPE, .TRUE., vari, sovar, llrgz=.FALSE.)
    IF (status /= 0) CALL ERROR_BI(grid_error(status), substr)

    CALL RGTOOL_CONVERT(sovar(1), dat, phgrid, 'xyzn')

    IF (.NOT. ASSOCIATED(hsurf_full)) &
         ALLOCATE(hsurf_full(SIZE(dat,1),SIZE(dat,2)))
    hsurf_full(:,:) = dat(:,:,1,1)
    DEALLOCATE(dat)
    NULLIFY(dat)

    ! CLEAN
    DO ix= 1, SIZE(vari)
       CALL INIT_NCVAR(vari(ix))
       CALL INIT_NCVAR(sovar(ix))
    END DO
    DEALLOCATE(vari)
    NULLIFY(vari)
    DEALLOCATE(sovar)
    NULLIFY(sovar)

    NULLIFY(hsurf_PSD)

    CALL INIT_GEOHYBGRID(basegrid)
    CALL INIT_GEOHYBGRID(basehgrid)

  END SUBROUTINE calc_p_hsurf_full
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------


#else

  IMPLICIT NONE
  PRIVATE

  INTERFACE mmd2way_child_setup
    MODULE PROCEDURE mmd2way_child_setup
  END INTERFACE mmd2way_child_setup

  INTERFACE mmd2way_child_init_memory
    MODULE PROCEDURE mmd2way_child_init_memory
  END INTERFACE mmd2way_child_init_memory

  INTERFACE mmd2way_child_init_loop
    MODULE PROCEDURE mmd2way_child_init_loop
  END INTERFACE mmd2way_child_init_loop

  INTERFACE mmd2way_child_global_end
    MODULE PROCEDURE mmd2way_child_global_end
  END INTERFACE mmd2way_child_global_end

  INTERFACE mmd2way_child_write_output
    MODULE PROCEDURE mmd2way_child_write_output
  END INTERFACE mmd2way_child_write_output

  INTERFACE mmd2way_child_write_restart
    MODULE PROCEDURE mmd2way_child_write_restart
  END INTERFACE mmd2way_child_write_restart

  INTERFACE mmd2way_child_read_restart
    MODULE PROCEDURE mmd2way_child_read_restart
  END INTERFACE mmd2way_child_read_restart

  INTERFACE mmd2way_child_free_memory
    MODULE PROCEDURE mmd2way_child_free_memory
  END INTERFACE mmd2way_child_free_memory

  PUBLIC :: mmd2way_child_setup
  PUBLIC :: mmd2way_child_init_memory
  PUBLIC :: mmd2way_child_init_loop
  PUBLIC :: mmd2way_child_global_end
  PUBLIC :: mmd2way_child_write_output
  PUBLIC :: mmd2way_child_write_restart
  PUBLIC :: mmd2way_child_read_restart
  PUBLIC :: mmd2way_child_free_memory


CONTAINS

  SUBROUTINE mmd2way_child_setup
  END SUBROUTINE mmd2way_child_setup

  SUBROUTINE mmd2way_child_init_memory
  END SUBROUTINE mmd2way_child_init_memory

  SUBROUTINE mmd2way_child_init_loop
  END SUBROUTINE mmd2way_child_init_loop

  SUBROUTINE  mmd2way_child_global_end
  END  SUBROUTINE mmd2way_child_global_end

  SUBROUTINE  mmd2way_child_write_output(flag)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

    IF (flag == 0) RETURN

  END  SUBROUTINE mmd2way_child_write_output

  SUBROUTINE  mmd2way_child_write_restart(flag)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

    IF (flag == 0) RETURN

  END  SUBROUTINE mmd2way_child_write_restart

  SUBROUTINE  mmd2way_child_read_restart(flag)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

    IF (flag == 0) RETURN

  END  SUBROUTINE mmd2way_child_read_restart

  SUBROUTINE mmd2way_child_free_memory
  END SUBROUTINE mmd2way_child_free_memory

#endif

END MODULE messy_mmd2way_child_si
