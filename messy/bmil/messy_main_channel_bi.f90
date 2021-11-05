#include "messy_main_ppd_bi.inc"

! **********************************************************************
MODULE messy_main_channel_bi
! **********************************************************************

  ! -------------------------------------------------------------------
  !  CHANNEL INTERFACE SHARED BETWEEN DIFFERENT BASEMODELS
  ! -------------------------------------------------------------------

  USE messy_main_channel_error_bi, ONLY: channel_halt
#if defined(ECHAM5)
  USE mo_grib,                  ONLY: L_BM_ORIG_OUTPUT
#endif
  USE messy_main_timer_event,   ONLY: &
                                      time_event, io_time_event  &
                                    , TIME_INC_MONTHS, TRIG_EXACT
  USE messy_main_constants_mem, ONLY: DP
  USE messy_main_channel,       ONLY: NMAXCHANNELS, STRLEN_CHANNEL, NCHANNEL &
                                    , modstr, IOMODE_OUT, IOMODE_RST         &
                                    , AMODE_WRITE, AMODE_READ, AMODE_INIT    &
                                    , REPR_UNDEF, DIMID_UNDEF, FTYPE_MAXIMUM
  USE messy_main_channel_mem,   ONLY: n_dom
  USE messy_main_channel_attributes, ONLY: t_attribute_list
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  USE messy_main_tools,         ONLY: PTR_4D_ARRAY, PTR_1D_ARRAY

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE

  PUBLIC :: DIMID_UNDEF, REPR_UNDEF, IOMODE_OUT, IOMODE_RST
  PUBLIC :: n_dom

  ! DIMENSION IDs ----------------------------------------------------------
  ! - TIME
  INTEGER, SAVE, PUBLIC :: DIMID_TIME     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_TBNDS    = DIMID_UNDEF
  ! - GRIDPOINT
  INTEGER, SAVE, PUBLIC :: DIMID_LEV      = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_ILEV     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_ILEVp1   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_LON      = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_LAT      = DIMID_UNDEF
  ! - SPECTRAL
  INTEGER, SAVE, PUBLIC :: DIMID_COMPLEX  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_NSP      = DIMID_UNDEF
  ! - FOURIER (ANTI-SYMMETRIC, SYMMETRIC)
  INTEGER, SAVE, PUBLIC :: DIMID_NHGL     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_NMP1     = DIMID_UNDEF
  ! - SPECIAL
  INTEGER, SAVE, PUBLIC :: DIMID_BELOWSF  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_2LEV     = DIMID_UNDEF
  !
  ! SPECIFIC BASEMODELS / SUBMODELS:
  ! - ATTILA
  INTEGER, SAVE, PUBLIC :: DIMID_NCELL_ATTILA = DIMID_UNDEF
  !
  ! - COSMO
  INTEGER, SAVE, PUBLIC :: DIMID_SRLON  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SRLAT  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_BNDS   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_DAV    = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_PRESS  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_ALTIT  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_FIXHT  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_H2m    = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_H10m   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_HTOA   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_WBT13  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SOIL   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SOIL1  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SOIL2  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SNOW   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SNOW1  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_TLV    = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_TKE    = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_CANOPY = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_CANOPY1= DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SAT8   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SAT32  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_2N     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_N5     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_N9     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_NHORI  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_WINDSEC = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_ECHOTOP = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_NEL_COMP= DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_TILES   = DIMID_UNDEF
  !
#if defined(BLANK) || defined(VERTICO)
  ! - BLANK, VERTICO
  INTEGER, SAVE, PUBLIC :: DIMID_AL     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_MID_BND= REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: ARRAY        = REPR_UNDEF
#endif

#ifdef I2CINC
  ! - I2CINC
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_IE        = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_JE        = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2Ci_IE       = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2Ci_JE       = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_KEDIM     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_KE1DIM    = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_KE1LM     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_KELM      = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_NHORI     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_MONTH     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_MONTHx9   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_9         = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_SOIL1     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_SOIL2     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_SPEC      = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_SPECMON   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_IE_IN     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_JE_IN     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_IE_LIN    = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_JE_LIN    = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_KE_IN     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_KE1_IN    = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_KEDIM_IN  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_KEDIM1_IN = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_SOIL1_IN  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_I2C_SOIL2_IN  = DIMID_UNDEF
#endif

!!#D clams +
  ! - CLaMS
  INTEGER, SAVE, PUBLIC :: DIMID_TRAJ         = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_NTASKS       = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SHUFFLE      = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_THETA        = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_MIX_GRID     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_POS          = DIMID_UNDEF
  ! - CLaMS - SEDI
  INTEGER, SAVE, PUBLIC :: DIMID_SEDI_PARTICLE = DIMID_UNDEF
!!#D clams -

#ifdef ICON
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_LEVS(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_LEVSP1(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_NCELLS(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_NVERTS(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_NEDGES(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_SOILLEVS(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_SOILLEVSP1(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_SNOWLEVS(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_SNOWLEVSP1(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_ICE(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_ONE(:)
#endif

  ! REPRESENTATION IDs -------------------------------------------------------
  ! - TIME
  INTEGER, SAVE, PUBLIC :: REPR_TIMEBNDS    = REPR_UNDEF
  ! - GRIDPOINT
  INTEGER, SAVE, PUBLIC :: GP_3D_MID        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_INT        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_2D_HORIZONTAL = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_1LEV       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_1D_LEV        = REPR_UNDEF
  ! - SPECTRAL
  INTEGER, SAVE, PUBLIC :: SP_3D_MID        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: SP_3D_INT        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: SP_2D_HORIZONTAL = REPR_UNDEF
  ! - FOURIER
  INTEGER, SAVE, PUBLIC :: FAS_MID          = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: FAS_INT          = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: FAS_MID_ZM       = REPR_UNDEF
  ! - SPECIAL
  INTEGER, SAVE, PUBLIC :: SCALAR           = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_BELOWSF    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_2LEV       = REPR_UNDEF

  ! SPECIFIC BASEMODELS / SUBMODELS:
  ! - COSMO
  INTEGER, SAVE, PUBLIC :: GP_4D_MID_TLV    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_INT_TLV    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_SOIL1_TLV  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_SOIL2_TLV  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_SNOW_TLV   = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_SNOW1_TLV  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_MID_BND    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_INT_BND    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_MID_SLON   = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_MID_SLAT   = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_MID_P      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_MID_Z      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_SOIL       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_SOIL1      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_SOIL2      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_SNOW       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_SNOW1      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_NHORI      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_HORIZ_T    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_HORIZ_BND  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_HORIZ_DAV  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_2D_SLONLAT    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_2D_LONSLAT    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_2D_HTOTLON    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_2D_TOT        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_2D_LAT_BND    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_1D_ILEV       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_1D_LAT        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_2D_LAT        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_SAT8       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_SAT32      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_1D_ILEVp1     = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_ECHOTOP    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_NEL_COMP   = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_HOR_N5       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_HOR_N9       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_WINDSEC      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_1D_WINDSEC      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_1D_WINDSEC_BNDS = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_1D_SOIL1_BNDS   = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_TILES        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_MID_TILES    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_INT_TILES    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_SOIL1_TILES  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_SOIL2_TILES  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_SNOW_TILES   = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_SNOW1_TILES  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_TILES_TLV    = REPR_UNDEF
#ifdef I2CINC
  ! - I2CINC
  INTEGER, SAVE, PUBLIC :: REPR_I2C_2D       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_2Dx9     = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2Ci_2D      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_3D_MID   = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_NHORI    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_3D_INT   = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_MONTH    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_MONTHx9  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_3D_MID_LM    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_3D_INT_LM    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_SOIL1    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_SOIL2    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_SPEC     = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_SPECMON  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_2D_IN    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_2D_IN_TOT= REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_3D_MID_IN    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_3D_MIDK_IN   = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_3D_INT_IN    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_3D_INTK_IN   = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_3D_MIX       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_3D_MIX1      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_SOIL1_IN = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_SOILIN1  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_I2C_SOIL2_IN = REPR_UNDEF
#endif
  !
  ! - ATTILA
  INTEGER, SAVE, PUBLIC :: LG_ATTILA        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: LG_FMM_POS       = REPR_UNDEF
  !
  ! - MPIOM
  INTEGER, SAVE, PUBLIC :: GP_3D_MPIOM           = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_MPIOM_INT       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_2D_MPIOM           = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_LEV_MPIOM       = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_LEV_MPIOM_INT   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_DEPTH_MPIOM     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_DEPTH_MPIOM_INT = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_LON_MPIOM       = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_LAT_MPIOM       = DIMID_UNDEF
  !
  ! - CESM1, for HOMME-SE dynamical core
  INTEGER, SAVE, PUBLIC :: DIMID_LAT_SE          = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_LON_SE          = DIMID_UNDEF

  ! - CLaMS
  INTEGER, SAVE, PUBLIC :: REPR_LG_CLAMS         = REPR_UNDEF
!!#D clams +
  INTEGER, SAVE, PUBLIC :: REPR_CLAMS_SHUFFLE    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_3DINP_CLAMS      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_3DINP_CLAMSTHETA = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_NTASKS           = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_MIX_GRID         = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPR_POS              = REPR_UNDEF
  ! - CLaMS - SEDI
  INTEGER, SAVE, PUBLIC :: REPR_SEDI_PARTICLE = REPR_UNDEF
!!#D clams -

#ifdef ICON
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UCELL_2D_SURFACE(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UCELL_3D_1LEV(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UCELL_3D_MID(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UCELL_3D_INT(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UCELL_3D_INT_HHL(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UCELL_3D_SNOW(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UCELL_3D_SNOW_HALF(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UCELL_3D_DEPTH_BELOW_LAND(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UCELL_3D_DEPTH_BELOW_LAND_P1(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UCELL_3D_GENERIC_ICE(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UEDGE_2D_SURFACE(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UEDGE_3D_MID(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UEDGE_3D_INT(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UVERT_2D_SURFACE(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UVERT_3D_MID(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UVERT_3D_INT(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: ICON_SCALAR(:)

  TYPE t_piommap
     INTEGER, DIMENSION(:,:), POINTER, PUBLIC :: gp_start => NULL()
     INTEGER, DIMENSION(:,:), POINTER, PUBLIC :: gp_cnt   => NULL()
     INTEGER, DIMENSION(:,:), POINTER, PUBLIC :: gp_meml  => NULL()
     INTEGER, DIMENSION(:,:), POINTER, PUBLIC :: gp_memu  => NULL()
  END TYPE t_piommap
  TYPE(t_piommap), DIMENSION(:), POINTER, PUBLIC :: piommap => NULL()
#endif

  ! DECOMPOSITION TYPES
  INTEGER, PARAMETER, PUBLIC :: DC_BC = 0  ! BROADCAST (IO-PE ONLY)
  INTEGER, PARAMETER, PUBLIC :: DC_GP = 1  ! GRIDPOINT
  INTEGER, PARAMETER, PUBLIC :: DC_SP = 2  ! SPECTRAL
  INTEGER, PARAMETER, PUBLIC :: DC_SA = 3  ! FOURIER (SYM. AND ANTISYM.)
  INTEGER, PARAMETER, PUBLIC :: DC_IX = 4  ! INDEX
  INTEGER, PARAMETER, PUBLIC :: DC_GP_MPIOM   =  5 ! GRIDPOINT MPIOM
  INTEGER, PARAMETER, PUBLIC :: DC_I2C        =  6 ! GRIDPOINT INT2COSMO
  INTEGER, PARAMETER, PUBLIC :: DC_I2C_IN     =  7 ! GRIDPOINT INT2COSMO INPUT
  INTEGER, PARAMETER, PUBLIC :: DC_MMD_OUT    =  8 ! GRIDPOINT CHILD OUTPUT
                                                   ! FIELD TRANSFERED TO PARENT
  INTEGER, PARAMETER, PUBLIC :: DC_AG         = 10 ! ALL_GATHER
  INTEGER, PARAMETER, PUBLIC :: DC_MMD_OUTACC = 11 ! GRIDPOINT ACCUMULATED
                                                   ! CHILD OUTPUT
  ! ICON GRIDPOINT AGGREGATION ON REGULAR LON-LAT GRIDS
  INTEGER, PARAMETER, PUBLIC :: DC_AGG_SUM     = 12 ! SUM e.g for average
  INTEGER, PARAMETER, PUBLIC :: DC_AGG_PDF     = 13 ! handling PDFs
  INTEGER, PARAMETER, PUBLIC :: DC_LATLON_NONE = 14 ! alreday (manually) synced
  INTEGER, PARAMETER, PUBLIC :: DC_TRIX        = 15 ! LaMETTA trajectory index

  ! GP-DECOMPOSITION INFORMATION
#if defined(ECHAM5)
  INTEGER, PARAMETER, PUBLIC             :: gp_nseg = 2
#endif
#ifdef COSMO
  INTEGER, PARAMETER, PUBLIC             :: gp_nseg = 1
  INTEGER                                :: notl        ! number of time levels
#endif
#ifdef CESM1
  INTEGER, PARAMETER, PUBLIC             :: gp_nseg = 1 ! is this correct?
#endif
#ifdef ICON
  INTEGER, PARAMETER, PUBLIC             :: gp_nseg = 1 ! is this correct?
#endif
#if defined(VERTICO) || defined(MBM_CLAMS) || defined(MBM_MPIOM)  || defined(BLANK)
  INTEGER, PARAMETER, PUBLIC             :: gp_nseg = 1
#endif
#if defined(MESSYDWARF)
  INTEGER, PARAMETER, PUBLIC             :: gp_nseg = 1
#endif
  !
  INTEGER, DIMENSION(:,:), POINTER, PUBLIC :: gp_start => NULL()
  INTEGER, DIMENSION(:,:), POINTER, PUBLIC :: gp_cnt   => NULL()
  INTEGER, DIMENSION(:,:), POINTER, PUBLIC :: gp_meml  => NULL()
  INTEGER, DIMENSION(:,:), POINTER, PUBLIC :: gp_memu  => NULL()

  ! OUTPUT TIME
  ! qqq need for PUBLIC needs to be checked
  TYPE(time_event), DIMENSION(:), ALLOCATABLE, PUBLIC :: OUTPUT_EVENT
  LOGICAL,          DIMENSION(:), ALLOCATABLE, PUBLIC :: LOUTPUT_NOW
  LOGICAL, SAVE                               :: LFORCE_NEW_OUTPUT = .FALSE.
  INTEGER, SAVE,    DIMENSION(:), ALLOCATABLE         :: I_PATCH

  TYPE(PTR_1D_ARRAY), DIMENSION(:), ALLOCATABLE :: TIME_BNDS
  INTEGER, PARAMETER :: tbnds = 2

#ifdef ECHAM5
  ! SPECIAL INDICES FOR ECHAM5
  ! NOTE: CHANNEL OUTPUT OF tdiag, tdiag_gp, nudg, nudg_gp
  !       WILL BE SYNCHRONIZED TO THE SHORTEST OUTPUT INTERVAL (IF PRESENT !)
  INTEGER, SAVE :: js_tdiag     = 0
  INTEGER, SAVE :: js_tdiag_gp  = 0
  INTEGER, SAVE :: js_nudg      = 0
  INTEGER, SAVE :: js_nudg_gp   = 0
#endif
  ! TRIGGER NEW FILE
  TYPE(time_event), DIMENSION(:), ALLOCATABLE         :: TNF_EVENT
  LOGICAL,          DIMENSION(:), ALLOCATABLE, PUBLIC :: LTNF_NOW

  ! SPECIAL INDICES FOR COSMO
  ! NOTE: CHANNEL OUTPUT OF COSMOm, COSMOp and COSMOz
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC :: js_COSMOm
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC :: js_COSMOz
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC :: js_COSMOp
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC :: js_COSMOs
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC :: js_COSMOc

  ! SPECIAL RESTART ATTRIBUTES
  TYPE(t_attribute_list), POINTER, SAVE :: restart_att => NULL()

  ! SPECIAL FOR PARALLEL I/O OPTIMISATION
  INTEGER, DIMENSION(:), ALLOCATABLE :: nodeid

#if defined(ECHAM5)
  ! <EXP_NAME> + _YYYYMMDD_HHMM_  = 14 + 15 (+1)
  INTEGER, PARAMETER :: BASENAME_LEN = 30
  ! YYYYMMDD hhmm[ss]
  INTEGER, PARAMETER :: DATE_TIME_STR_LEN = 15
#endif
#if (defined COSMO) || defined(BLANK) || defined(MBM_CLAMS) || defined(CESM1) || defined(VERTICO) || defined(MBM_MPIOM)
  ! <EXP_NAME> + _YYYYMMDD_HHMMSS_  = 14 + 17 (+1)
  INTEGER, PARAMETER :: BASENAME_LEN = 32
  ! YYYYMMDD hhmmss
  INTEGER, PARAMETER :: DATE_TIME_STR_LEN = 17
#endif
#if defined(ICON)
  ! <EXP_NAME> + _YYYYMMDD_HHMMSSsss_ + <CHANNEL>  = 14 + 20 (+1)
  ! xy = domain/patch number
  INTEGER, PARAMETER :: BASENAME_LEN = 35
  ! YYYYMMDD hhmmss.sss
  INTEGER, PARAMETER :: DATE_TIME_STR_LEN = 20
#endif

  ! ====================================================================
  ! FOR CPL-NAMELIST
  !
  TYPE t_channel_timer
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname    = ''       ! CHANNEL NAME
     TYPE(io_time_event)          :: io_event = &
          io_time_event(1, TIME_INC_MONTHS,TRIG_EXACT,0) ! DEFAULT
  END TYPE t_channel_timer
  !
  TYPE(t_channel_timer),                          SAVE :: TIMER_DEFAULT
  TYPE(t_channel_timer), DIMENSION(NMAXCHANNELS), SAVE :: TIMER_CHANNEL
  !
  TYPE(t_channel_timer),                          SAVE :: TIMER_TNF_DEFAULT
  TYPE(t_channel_timer), DIMENSION(NMAXCHANNELS), SAVE :: TIMER_TNF_CHANNEL
  !
  ! METHOD FOR PSEUDO PARALLEL OUTPUT
  INTEGER, DIMENSION(FTYPE_MAXIMUM),              SAVE :: CHANNEL_PPIO_TYPE
  ! ====================================================================

#if defined(ECHAM5)
  ! ====================================================================
  ! SPECIAL HANDLING FOR ECHAM5-STREAM-ELEMENTS WITH laccu = .TRUE.
  TYPE(PTR_4D_ARRAY), DIMENSION(:), POINTER, SAVE :: ptr_accu
  INTEGER,                                   SAVE :: naccu = 0
  ! ====================================================================
#endif

#ifdef COSMO
  ! POINTER FIELDS FOR COSMO OUTPUT VARIABLES per CHANNEL
  TYPE cosmo_out_fields
     CHARACTER(LEN=9)                                :: label
     TYPE(PTR_4D_ARRAY), DIMENSION(:), POINTER       :: vars
  END TYPE cosmo_out_fields
  !
  TYPE cosmo_output_list
     TYPE(cosmo_out_fields)           :: this
     TYPE(cosmo_output_list), POINTER :: next => NULL()
  END TYPE cosmo_output_list
  !
  PUBLIC :: cosmo_output_list,  cosmo_out_fields
  !
  TYPE(cosmo_output_list), POINTER, SAVE, PUBLIC :: COSMOOUT => NULL()
  !
  ! FORCE OUTPUT CALCULATIONS in COSMO
  LOGICAL, PUBLIC :: L_FORCE_calcout = .TRUE.
#endif

#if defined(COSMO) || defined(BLANK) || defined(MBM_CLAMS) || defined(CESM1) || defined(VERTICO) || defined(MBM_MPIOM) || defined(ICON)
  ! SWITCH OFF BASEMODEL OUTPUT AND RESTARTS
  LOGICAL, PUBLIC :: L_BM_ORIG_OUTPUT
#endif

  ! CALCULATE COSMO DIAGNOSIC VARIABLE EACH TIME STEP
  ! (independent of output interval)
  LOGICAL, PUBLIC :: L_CALCOUT_EACH_STEP = .FALSE.

  ! CALLED FROM BMIL (messy_main_control_<basemodel>.f90)
  PUBLIC   :: main_channel_setup
                            !-> channel_init_restart_bi
  PUBLIC   :: main_channel_initialize
  !                         |-> main_channel_read_nml_ctrl (CORE !)
  !                         |-> main_channel_initialize_gatts
  !                         |-> main_channel_initialize_dims
  !                         |-> main_channel_initialize_reprs
  !                         |-> main_channel_initialize_parallel_io
  PUBLIC   :: main_channel_init_memory
  !                         |-> associate_streams_to_channels
  !                         |-> set_COSMO_ORI_attributes
  PUBLIC   :: main_channel_init_coupling
  !                         | case 0 (COSMO only!)
  !                         |-> messy_channel_cosmo_output
  !                                    |-> make_cosmo_channel
  !                         | case 1
  !                         |-> modify_attributes
  !                         | case 2
  !                         |-> messy_channel_cosmo_vcatt
  !                         |-> fixate_channels (CORE !)
  !                         !-> distribute_channels_io
  !                         |-> main_channel_init_timer
  !                                    |-> main_channel_read_nml_cpl
  !                         |-> write_attribute
  !                         |-> write_dimension
  !                         |-> write_representation
  !                         |-> write_channel
  PUBLIC   :: main_channel_global_start
  !                         |-> set_channel_output (COSMOc)
  !                         |-> main_channel_update_timer
  !                             |-> trigger_channel_output(_part)
  PUBLIC   :: main_channel_write_output
  !                         |-> trigger_channel_output(_full)
  !                         !-> channel_write_output_bi
  PUBLIC   :: main_channel_write_restart
  !                         |-> channel_write_output_bi
  PUBLIC   :: main_channel_free_memory
  !                         |-> channel_finish_io
  !                         |-> clean_representations (CORE !)
  !                         |-> clean_dimensions (CORE !)
  !                         |-> clean_channels (CORE !)
  !                         |-> write_channel
  PUBLIC   :: main_channel_read_restart
  !                         !-> initialize_restart_attributes
  !
#ifdef COSMO
  PUBLIC   :: messy_COSMO_create_channel
  !
  PUBLIC   :: messy_COSMO_dealloc_meteofields
  PUBLIC   :: messy_COSMO_channel_readrst
#ifdef I2CINC
  PUBLIC   :: messy_I2CINC_channel_alloc_lm
  !
  PUBLIC   :: messy_I2CINC_channel_alloc_cg
  !
  PUBLIC   :: messy_I2CINC_channel_make_MMDC4
#endif
#endif
#ifdef COSMOv5s5
  PUBLIC   ::  messy_channel_alloc_block
  PUBLIC   ::  messy_channel_dealloc_block
  PUBLIC   ::  messy_channel_register_block_fields
  PUBLIC   ::  messy_channel_register_copy
#endif
#ifdef ICON
  PUBLIC   :: main_channel_set_domain
#endif
  !
  ! CALLED FROM SMIL
  !
  PUBLIC   :: p_bcast_attribute
  !
  ! SHARED BETWEEN DIFFERENT BASEMODELS
  !
  !PRIVATE :: initialize_restart_attributes  ! INIT RESTART ATTRIBUTES
  !PRIVATE :: main_channel_init_timer        ! INITIALIZE OUTPUT TIMERS
  !PRIVATE :: main_channel_update_timer
  !PRIVATE :: main_channel_read_nml_cpl
  !PRIVATE :: bi_decompose
  !PRIVATE :: bi_vector
  !PRIVATE :: distribute_channels_io
  !
  !PRIVATE :: channel_write_output_bi
  !                         !-> reset_accu_stream_elements(1)
  !                             (reset accu elements to instantaneous values)
  !                         !-> reset_accu_cosmo_elements(1)
  !                             (reset accu elements to instantaneous values)
  !                         |-> update_channels(1) (CORE !)
  !                             ( ACCUMULATE 2ndary DATA )
  !                         |-> update_channels(2) (CORE !)
  !                             ( PREPARE FOR OUTPUT )
  !                         !-> channel_finish_io
  !                         !-> reset_accu_stream_elements(2)
  !                             (set accu elements to zero)
  !                         !-> reset_accu_cosmo_elements(2)
  !                             (set accu elements to zero)
  !                         |-> update_channels(3) (CORE !)
  !                             ( RESET AFTER OUTPUT )
  !                         !-> initialize_restart_attributes
  !
  !PRIVATE :: channel_init_restart_bi
  !                         !-> initialize_parallel_io
  !                         !-> initialize_restart_attributes
  !
  ! BASEMODEL SPECIFIC: SEE messy_main_channel_xx.inc
  !PRIVATE :: main_channel_initialize_gatts  ! SET GLOBAL ATTRIBUTES
  !PRIVATE :: main_channel_initialize_dims   ! DEFINE DIMENSIONS
  !PRIVATE :: main_channel_initialize_reprs  ! DEFINE REPRESENTATIONS
  !
  ! ... ECHAM5
  !PRIVATE :: associate_streams_to_channels
  !PRIVATE :: reset_accu_stream_elements
  !PRIVATE :: timer_sync
  !
  ! ... COSMO
  !PRIVATE :: set_cosmo_channel_attributes
  !PRIVATE :: messy_channel_cosmo_output
  !PRIVATE :: reset_accu_cosmo_elements
  !PRIVATE :: set_COSMO_ORI_attributes
  !PRIVATE :: set_cosmo_attributes

CONTAINS

#if defined(ECHAM5)
#include "messy_main_channel_echam5.inc"
#endif

#ifdef COSMO
  !INCLUDE 'messy_main_channel_cosmo.inc'
#include "messy_main_channel_cosmo.inc"
#ifdef I2CINC
#include "messy_main_channel_int2cosmo.inc"
#endif
#endif

#ifdef CESM1
#include "messy_main_channel_cesm1.inc"
#endif

#ifdef MESSYDWARF
#include "messy_main_channel_dwarf.inc"
#endif

#ifdef BLANK
#if defined(MBM_BLANK) || defined (MBM_GWAVE) || defined (MBM_QBO) || defined (MBM_SPE)
#include "messy_main_channel_blank.inc"
#endif
#if defined(MBM_DISSOC)
#include "messy_main_channel_dissoc.inc"
#endif
#if defined(MBM_RAD)
#include "messy_main_channel_rad.inc"
#endif
#endif

#ifdef MBM_MPIOM
#include "messy_main_channel_mpiom.inc"
#endif

#ifdef MBM_CLAMS
#include "messy_main_channel_clams.inc"
#endif
#ifdef VERTICO
#include "messy_main_channel_vertico.inc"
#endif

#ifdef ICON
#include "messy_main_channel_icon.inc"
#endif

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES (MAIN ENTRY POINTS)
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
#ifndef ICON
  SUBROUTINE main_channel_setup

    IMPLICIT NONE

    CALL channel_init_restart_bi

  END SUBROUTINE main_channel_setup

#else

  SUBROUTINE main_channel_setup(flag, ndom)

    USE messy_main_channel_mem, ONLY: l_dom, n_dom_max, n_dom
    USE messy_main_bmluse_bi,   ONLY: max_dom_icon => max_dom

    IMPLICIT NONE

    INTEGER, INTENT(IN)           :: flag
    INTEGER, INTENT(IN), OPTIONAL :: ndom

    SELECT CASE(flag)

    CASE(1)
       l_dom = .TRUE.
       n_dom_max = max_dom_icon
       IF (PRESENT(ndom)) THEN
          IF (ndom > 0) THEN
             n_dom = ndom
          END IF
       END IF

    CASE(2)

       CALL channel_init_restart_bi

    END SELECT

END SUBROUTINE main_channel_setup
#endif
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize

    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_mpi_bi,   ONLY: p_parallel_io, p_io, p_bcast &
                                 , p_pe, p_all_comm, p_nprocs
    USE messy_main_channel,  ONLY: main_channel_read_nml_ctrl    &
                                 , NMAXCHANNELS, NMAXOBJECTS     &
                                 , NMAXADDCHANNEL, NMAXADDREF    &
                                 , NMAXADDATT                    &
                                 , ADD_CHANNEL, ADD_REF, OUT_DEFAULT &
                                 , ADD_ATT &
                                 , OUT_PREC &
                                 , OUT_CHANNEL, OUT_OBJECT, EXP_NAME &
                                 , EXEC_CHECKSUM &
                                 , L_FLUSH_IOBUFFER &
                                 , I_VERBOSE_LEVEL
    USE messy_main_tools,    ONLY: find_next_free_unit

#ifndef NOMPI
    USE messy_main_channel_io, ONLY: initialize_parallel_io
#ifdef HAVE_PNETCDF
    USE messy_main_channel_pnetcdf, ONLY: ch_pnetcdf_read_nml_ctrl &
                                        , NMAX_MPI_IO_HINTS        &
                                        , MPI_IO_HINT
#endif
#endif

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize'
    INTEGER :: status
    INTEGER :: iou
    INTEGER :: i

    CALL start_message_bi(modstr,'INITIALIZE CHANNELS',substr)

    ! READ CTRL-NAMELIST FOR FIXATION ( ... AT END OF INIT_COUPLING)
    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_channel_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    ! BROADCAST RESULTS
    CALL p_bcast(EXP_NAME, p_io)
    CALL p_bcast(EXEC_CHECKSUM, p_io)
    CALL p_bcast(L_FLUSH_IOBUFFER, p_io)
    CALL p_bcast(I_VERBOSE_LEVEL, p_io)
    DO i=1, NMAXADDCHANNEL
       CALL p_bcast(ADD_CHANNEL(i)%cname, p_io)
    END DO
    !
    DO i=1, NMAXADDREF
       CALL p_bcast(ADD_REF(i)%cname1, p_io)
       CALL p_bcast(ADD_REF(i)%oname1, p_io)
       CALL p_bcast(ADD_REF(i)%cname2, p_io)
       CALL p_bcast(ADD_REF(i)%oname2, p_io)
    END DO
    !
    DO i=1, NMAXADDATT
       CALL p_bcast(ADD_ATT(i)%cname,   p_io)
       CALL p_bcast(ADD_ATT(i)%oname,   p_io)
       CALL p_bcast(ADD_ATT(i)%attname, p_io)
       CALL p_bcast(ADD_ATT(i)%atttype, p_io)
       CALL p_bcast(ADD_ATT(i)%attval,  p_io)
       CALL p_bcast(ADD_ATT(i)%lforce,  p_io)
       CALL p_bcast(ADD_ATT(i)%lout,    p_io)
    END DO
    !
    CALL p_bcast(OUT_DEFAULT%cname, p_io)
    CALL p_bcast(OUT_DEFAULT%cio%ftype(:), p_io)
    CALL p_bcast(OUT_DEFAULT%cio%ntpf, p_io)
    CALL p_bcast(OUT_DEFAULT%oio%lrestart, p_io)
    CALL p_bcast(OUT_DEFAULT%oio%lignore, p_io)
    CALL p_bcast(OUT_DEFAULT%oio%lout(:), p_io)
    CALL p_bcast(OUT_DEFAULT%oio%range(:), p_io)
    !
    CALL p_bcast(OUT_PREC(:), p_io)
    !
    DO i=1, NMAXCHANNELS
       CALL p_bcast(OUT_CHANNEL(i)%cname, p_io)
       CALL p_bcast(OUT_CHANNEL(i)%cio%ftype(:), p_io)
       CALL p_bcast(OUT_CHANNEL(i)%cio%ntpf, p_io)
       CALL p_bcast(OUT_CHANNEL(i)%oio%lrestart, p_io)
       CALL p_bcast(OUT_CHANNEL(i)%oio%lignore, p_io)
       CALL p_bcast(OUT_CHANNEL(i)%oio%lout(:), p_io)
       CALL p_bcast(OUT_CHANNEL(i)%oio%range(:), p_io)
    END DO
    !
    DO i=1, NMAXOBJECTS
       CALL p_bcast(OUT_OBJECT(i)%cname, p_io)
       CALL p_bcast(OUT_OBJECT(i)%oname, p_io)
       CALL p_bcast(OUT_OBJECT(i)%io%lrestart, p_io)
       CALL p_bcast(OUT_OBJECT(i)%io%lignore, p_io)
       CALL p_bcast(OUT_OBJECT(i)%io%lout(:), p_io)
       CALL p_bcast(OUT_OBJECT(i)%io%range(:), p_io)
    END DO

    ! INTIALIZE GLOBAL ATTRIBUTES
    CALL main_channel_initialize_gatts

    ! INTIALIZE DIMENSIONS
    CALL main_channel_initialize_dims

    ! INTIALIZE REPRESENTATIONS
    CALL main_channel_initialize_reprs

    CALL end_message_bi(modstr,'INITIALIZE CHANNELS',substr)

#ifndef NOMPI
    CALL start_message_bi(modstr,'INITIALIZE PARALLEL I/O',substr)
#ifdef HAVE_PNETCDF
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL ch_pnetcdf_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! BROADCAST RESULTS
    DO i=1, NMAX_MPI_IO_HINTS
       CALL p_bcast(MPI_IO_HINT(i)%hint, p_io)
       CALL p_bcast(MPI_IO_HINT(i)%value, p_io)
    END DO
#endif

    CALL initialize_parallel_io(status, p_pe, p_io, p_all_comm, p_nprocs)
    CALL channel_halt(substr, status)
    CALL end_message_bi(modstr,'INITIALIZE PARALLEL I/O',substr)
#endif

    ! moved here from init_coupling; otherwise L_BM_ORIG_OUTPUT is not defined
    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_channel_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast(L_BM_ORIG_OUTPUT, p_io)
    CALL p_bcast(L_CALCOUT_EACH_STEP, p_io)
    CALL p_bcast(CHANNEL_PPIO_TYPE, p_io)

  END SUBROUTINE main_channel_initialize
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_init_memory

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_init_memory'

    CALL start_message_bi(modstr,'INITIALIZE CHANNEL MEMORY',substr)

#if defined(ECHAM5)
    ! ASSOCIATE ECHAM5 STREAMS TO MESSy CHANNELS
    CALL associate_streams_to_channels
#endif

#ifdef COSMO
    CALL set_COSMO_ORI_attributes
#endif

#ifdef CESM1
    CALL messy_CESM_create_channel
#endif

#ifdef ICON
    CALL associate_var_lists_to_channels
#endif

    CALL end_message_bi(modstr,'INITIALIZE CHANNEL MEMORY',substr)

  END SUBROUTINE main_channel_init_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_init_coupling(flag)

    USE messy_main_mpi_bi,   ONLY: p_parallel_io, p_pe, p_io
    USE messy_main_channel,  ONLY: fixate_channels, write_channel &
                                 , write_attribute &
                                 , I_VERBOSE_LEVEL &
                                 , modify_attributes
    USE messy_main_channel_dimensions,  ONLY: write_dimension
    USE messy_main_channel_repr,        ONLY: write_representation &
                                            , write_representation_dc

#ifdef ICON
    USE messy_main_bmluse_bi,              ONLY: patch_info            &
#ifndef ICON_2_1_01
                                               , icell, iedge, ivert   &
#endif
                                               , output_mode           &
                                               , number_of_grid_used   &
                                               , p_patch, p_phys_patch &
                                               , t_phys_patch          &
                                               , collect_all_grid_info &
                                               , GRID_INFO_NONE        &
                                               , GRID_INFO_FILE        &
                                               , GRID_INFO_BCAST       &
                                               , l_output_phys_patch   &
                                               , t_reorder_info        &
                                               , nproma, my_process_is_io
#endif

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)         :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_init_coupling'
    INTEGER                     :: status
#ifdef ICON
    INTEGER                     :: jg
    INTEGER                     :: jl
    INTEGER                     :: n_dom_out
#endif
    INTEGER          :: js
    CHARACTER(LEN=STRLEN_CHANNEL) :: cname
    INTEGER          :: dom_id


    SELECT CASE(flag)

#ifdef COSMO
    CASE(0)
       CALL start_message_bi(modstr,'COSMO OUTPUT TO CHANNELS',substr)
       ! ASSOCIATE COSMO OUTPUT WITH CHANNELS
       ! MUST BE CALLED BEFORE FIXATE CHANNELS AND before init_output
       CALL messy_channel_cosmo_output
       CALL end_message_bi(modstr,'COSMO OUTPUT TO CHANNELS',substr)
#endif

    CASE(1)
       CALL modify_attributes(status)
       CALL channel_halt(substr, status)

    CASE(2)
#ifdef COSMO
       CALL start_message_bi(modstr,'COSMO SET VCOORD ATTRIBUTS',substr)
       CALL messy_channel_cosmo_vcatt
       ! moved to init_cpl (0) CALL messy_channel_cosmo_output
       CALL end_message_bi(modstr,'COSMO SET VCOORD ATTRIBUTES',substr)
#endif
#ifdef ICON
       CALL main_channel_icon_attributes
#endif
       CALL start_message_bi(modstr,'FIXATE CHANNELS PART 1',substr)
       CALL fixate_channels(status,1)
       CALL channel_halt(substr, status)
       CALL end_message_bi(modstr,'FIXATE CHANNELS PART 1',substr)

       CALL main_channel_make_cf_conform

       CALL start_message_bi(modstr,'FIXATE CHANNELS PART 2',substr)
       CALL fixate_channels(status,2)
       CALL channel_halt(substr, status)
       CALL end_message_bi(modstr,'FIXATE CHANNELS PART 2',substr)

       CALL start_message_bi(modstr,'DISTRIBUTE CHANNEL OUTPUT',substr)
       CALL distribute_channels_io
       CALL end_message_bi(modstr,'DISTRIBUTE CHANNEL OUTPUT',substr)

#ifdef ICON
! QQQ: why is this needed here?
! option 1: if the same thing is done (under certain conditions) already in
!           icon, then add as addtional (.OR.) condition that #ifdef MESSY
! ?1: if in icon more than this is done, does that harm icon (#ifdef MESSY)
!     or can it be "cut" (e.g. with #ifdef MESSY and no icon output then return)
! ?2: can this logic be connected directly with the L_BM_ORIG_OUTPUT switch?

       ! set some important patch informations,
       ! in case original ICON "name list output" is disabled.
       ! in that case patch_info is not allocated...
       IF (.NOT. output_mode%l_nml) THEN
          IF (.NOT. ALLOCATED(patch_info)) THEN
             ALLOCATE(patch_info(n_dom))

             n_dom_out = n_dom

             ! Set number of global cells/edges/verts and logical patch ID
             patch_loop: DO jg = 1, n_dom_out
                IF(l_output_phys_patch) THEN
                   patch_info(jg)%log_patch_id = p_phys_patch(jg)%logical_id
                   IF (.NOT. my_process_is_io()) THEN
                      patch_info(jg)%p_pat_c    => &
                           p_phys_patch(jg)%comm_pat_gather_c
                      patch_info(jg)%nblks_glb_c = &
                           (p_phys_patch(jg)%n_patch_cells-1)/nproma + 1
                      patch_info(jg)%p_pat_e    => &
                           p_phys_patch(jg)%comm_pat_gather_e
                      patch_info(jg)%nblks_glb_e = &
                           (p_phys_patch(jg)%n_patch_edges-1)/nproma + 1
                      patch_info(jg)%p_pat_v    => &
                           p_phys_patch(jg)%comm_pat_gather_v
                      patch_info(jg)%nblks_glb_v = &
                           (p_phys_patch(jg)%n_patch_verts-1)/nproma + 1
                      patch_info(jg)%max_cell_connectivity = &
                           p_patch(patch_info(jg)%log_patch_id)%cells%max_connectivity
                   END IF
                ELSE
                   patch_info(jg)%log_patch_id = jg
                   IF (.NOT. my_process_is_io()) THEN
                      patch_info(jg)%p_pat_c    => &
                           p_patch(jg)%comm_pat_gather_c
                      patch_info(jg)%nblks_glb_c = &
                           (p_patch(jg)%n_patch_cells_g-1)/nproma + 1
                      patch_info(jg)%p_pat_e    => &
                           p_patch(jg)%comm_pat_gather_e
                      patch_info(jg)%nblks_glb_e = &
                           (p_patch(jg)%n_patch_edges_g-1)/nproma + 1
                      patch_info(jg)%p_pat_v    => &
                           p_patch(jg)%comm_pat_gather_v
                      patch_info(jg)%nblks_glb_v = &
                           (p_patch(jg)%n_patch_verts_g-1)/nproma + 1
                      patch_info(jg)%max_cell_connectivity = &
                           p_patch(patch_info(jg)%log_patch_id)%cells%max_connectivity
                   END IF
                ENDIF
             ENDDO patch_loop

             ! Pure I/O PEs may skip this...
             ! Go over all output domains
             DO jg = 1, n_dom_out
                IF (patch_info(jg)%grid_info_mode == GRID_INFO_BCAST) THEN
                   CALL collect_all_grid_info(&
                        p_patch(patch_info(jg)%log_patch_id), patch_info(jg))
                END IF
             END DO

             DO jg = 1, n_dom_out

                jl = patch_info(jg)%log_patch_id

                IF(.NOT.my_process_is_io()) THEN
                   ! Set reorder_info on work and test PE
#ifdef ICON_2_1_01
                   CALL set_reorder_info(jg &
                        , p_patch(jl)%n_patch_cells_g &
                        , p_patch(jl)%n_patch_cells   &
                        , p_patch(jl)%cells%decomp_info%owner_mask &
                        , p_patch(jl)%cells%phys_id   &
                        , p_patch(jl)%cells%decomp_info%glb_index &
                        , patch_info(jg)%cells )

                   CALL set_reorder_info(jg &
                        , p_patch(jl)%n_patch_edges_g &
                        , p_patch(jl)%n_patch_edges   &
                        , p_patch(jl)%edges%decomp_info%owner_mask &
                        , p_patch(jl)%edges%phys_id   &
                        , p_patch(jl)%edges%decomp_info%glb_index &
                        , patch_info(jg)%edges )

                   CALL set_reorder_info(jg &
                        , p_patch(jl)%n_patch_verts_g &
                        , p_patch(jl)%n_patch_verts   &
                        , p_patch(jl)%verts%decomp_info%owner_mask &
                        , p_patch(jl)%verts%phys_id   &
                        , p_patch(jl)%verts%decomp_info%glb_index  &
                        , patch_info(jg)%verts )
#else
                   CALL set_reorder_info(jg &
                        , p_patch(jl)%n_patch_cells_g &
                        , p_patch(jl)%n_patch_cells   &
                        , p_patch(jl)%cells%decomp_info%owner_mask &
                        , p_patch(jl)%cells%phys_id   &
                        , p_patch(jl)%cells%decomp_info%glb_index &
                        , patch_info(jg)%ri(icell)    &
                        , patch_info(jg)%grid_info(icell))

                   CALL set_reorder_info(jg &
                        , p_patch(jl)%n_patch_edges_g &
                        , p_patch(jl)%n_patch_edges   &
                        , p_patch(jl)%edges%decomp_info%owner_mask &
                        , p_patch(jl)%edges%phys_id   &
                        , p_patch(jl)%edges%decomp_info%glb_index &
                        , patch_info(jg)%ri(iedge)    &
                        , patch_info(jg)%grid_info(iedge))

                   CALL set_reorder_info(jg &
                        , p_patch(jl)%n_patch_verts_g &
                        , p_patch(jl)%n_patch_verts   &
                        , p_patch(jl)%verts%decomp_info%owner_mask &
                        , p_patch(jl)%verts%phys_id   &
                        , p_patch(jl)%verts%decomp_info%glb_index  &
                        , patch_info(jg)%ri(ivert)    &
                        , patch_info(jg)%grid_info(ivert))
#endif
                   ! Set grid_filename on work and test PE
                   patch_info(jg)%grid_filename = &
                        TRIM(p_patch(jl)%grid_filename)
                   ! Set UUID on work and test PE
                   patch_info(jg)%grid_uuid = p_patch(jl)%grid_uuid
                   ! Set information about numberOfGridUsed on work and test PE
                   patch_info(jg)%number_of_grid_used = number_of_grid_used(jl)

                   patch_info(jg)%max_cell_connectivity = &
                        p_patch(jl)%cells%max_connectivity
                ENDIF
             END DO
          END IF
       END IF
#endif

       ! INITIALIZE CHANNEL OUTPUT TIMERS (via CPL-NAMELIST)
       CALL main_channel_init_timer

       IF (I_VERBOSE_LEVEL > 0) THEN
          IF (p_parallel_io) THEN
             CALL write_attribute(status)
             CALL channel_halt(substr, status)
          END IF

          IF (p_parallel_io) THEN
             CALL write_dimension(status)
             CALL channel_halt(substr, status)
          END IF

          IF (p_parallel_io) THEN
             CALL write_representation(status)
             CALL channel_halt(substr, status)
          END IF

          CALL write_representation_dc(status, p_pe)
          CALL channel_halt(substr, status)

          IF (p_parallel_io) THEN
             CALL write_channel(status)
             CALL channel_halt(substr, status)
          END IF
       END IF

       CALL end_message_bi(modstr,'FIXATE CHANNELS',substr)

    END SELECT

  END SUBROUTINE main_channel_init_coupling
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_global_start

#ifdef COSMO
    !  COSMO
    USE data_io,              ONLY: ngribout
    ! MESSy/SMCL
    USE messy_main_tools,     ONLY: int2str
    USE messy_main_timer,     ONLY: lfirst_cycle
    USE messy_main_channel,   ONLY: set_channel_output
    IMPLICIT NONE
#endif

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_channel_global_start'
#ifdef COSMO
    CHARACTER(LEN=3) :: str
    INTEGER          :: icc, status

       IF (lfirst_cycle) THEN
          ! force output of channels with constants at the beginning
          DO iCc = 1, ngribout
             IF (js_COSMOc(iCc) > 0) THEN
                CALL int2str(str,iCc)
                CALL set_channel_output(status,'COSMOc'//str, .TRUE.)
                CALL channel_halt(substr,status)
             ENDIF
          ENDDO
       ENDIF
#endif

    CALL main_channel_update_timer

  END SUBROUTINE main_channel_global_start
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_write_output

    USE messy_main_channel,    ONLY: trigger_channel_output

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_write_output'
    INTEGER :: status

    CALL trigger_channel_output(status, LOUTPUT_NOW, LTNF_NOW &
         , LFORCE_NEW_OUTPUT) ! full
    CALL channel_halt(substr, status)
    LFORCE_NEW_OUTPUT = .FALSE.

    CALL channel_write_output_bi(IOMODE_OUT)

  END SUBROUTINE main_channel_write_output
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_write_restart

    IMPLICIT NONE

    CALL channel_write_output_bi(IOMODE_RST)

  END SUBROUTINE main_channel_write_restart
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_free_memory

    USE messy_main_mpi_bi,              ONLY: p_parallel_io, p_pe
    USE messy_main_channel_io,          ONLY: channel_finish_io
    USE messy_main_channel_dimensions,  ONLY: clean_dimensions
    USE messy_main_channel_repr,        ONLY: clean_representations
    USE messy_main_channel,             ONLY: clean_channels, write_channel

    IMPLICIT NONE

#ifdef COSMO
    INTRINSIC :: ALLOCATED
#endif

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_free_memory'
    INTEGER :: status
#ifdef ICON
    INTEGER :: jg
#endif

    CALL start_message_bi(modstr,'FINISH CHANNELS',substr)

    ! ------------------------------------------
    ! CLOSE ALL OUTPUT FILES
    ! ------------------------------------------
    CALL channel_finish_io(status,  IOMODE_OUT, .TRUE., p_pe)
    CALL channel_halt(substr, status)

    CALL clean_representations(status)
    CALL channel_halt(substr, status)

    CALL clean_dimensions(status)
    CALL channel_halt(substr, status)

    CALL clean_channels(status)
    CALL channel_halt(substr, status)

    IF (p_parallel_io) THEN
       CALL write_channel(status)
       CALL channel_halt(substr, status)
    END IF

#ifdef COSMO
    ! CLEAN COSMO output channel indices
    IF (ALLOCATED(js_COSMOm)) DEALLOCATE(js_COSMOm)
    IF (ALLOCATED(js_COSMOp)) DEALLOCATE(js_COSMOp)
    IF (ALLOCATED(js_COSMOz)) DEALLOCATE(js_COSMOz)
    IF (ALLOCATED(js_COSMOs)) DEALLOCATE(js_COSMOs)
    IF (ALLOCATED(js_COSMOc)) DEALLOCATE(js_COSMOc)
#endif

    DEALLOCATE(LOUTPUT_NOW)
    DEALLOCATE(TIME_BNDS)
    DEALLOCATE(LTNF_NOW)
    DEALLOCATE(I_PATCH)

#ifndef ICON
    IF (ASSOCIATED(gp_start)) THEN
       DEALLOCATE(gp_start) ; NULLIFY(gp_start)
    END IF
    IF (ASSOCIATED(gp_cnt)) THEN
       DEALLOCATE(gp_cnt) ; NULLIFY(gp_cnt)
    END IF
    IF (ASSOCIATED(gp_meml)) THEN
       DEALLOCATE(gp_meml) ; NULLIFY(gp_meml)
    END IF
    IF (ASSOCIATED(gp_memu)) THEN
       DEALLOCATE(gp_memu) ; NULLIFY(gp_memu)
    END IF
#else
    NULLIFY(gp_start)
    NULLIFY(gp_cnt)
    NULLIFY(gp_meml)
    NULLIFY(gp_memu)
    DO jg = 1, n_dom
       DEALLOCATE(piommap(jg)%gp_start); NULLIFY(piommap(jg)%gp_start)
       DEALLOCATE(piommap(jg)%gp_cnt);   NULLIFY(piommap(jg)%gp_cnt)
       DEALLOCATE(piommap(jg)%gp_meml);  NULLIFY(piommap(jg)%gp_meml)
       DEALLOCATE(piommap(jg)%gp_memu);  NULLIFY(piommap(jg)%gp_memu)
    END DO
    DEALLOCATE(piommap); NULLIFY(piommap)
#endif

    IF (ALLOCATED(nodeid)) DEALLOCATE(nodeid)

    CALL end_message_bi(modstr,'FINISH CHANNELS',substr)

  END SUBROUTINE main_channel_free_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE (PART I)
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_write_output_bi(IOMODE)

    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_pe
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_timer,      ONLY: &
                                     YEAR, MONTH, DAY, HOUR, MINUTE, SECOND &
                                   , YEAR_NEXT, MONTH_NEXT, DAY_NEXT     &
                                   , HOUR_NEXT, MINUTE_NEXT, SECOND_NEXT &
                                   , YEAR_START, MONTH_START, DAY_START     &
                                   , HOUR_START, MINUTE_START, SECOND_START &
                                   , MILLISECOND_NEXT &
                                   , MILLISECOND_START &
                                   , delta_time, current_time_step &
                                   , time_step_len
    USE messy_main_channel_io, ONLY: channel_init_io           &
                                   , channel_write_header      &
                                   , channel_write_time        &
                                   , channel_write_data        &
                                   , channel_finish_io
    USE messy_main_channel,    ONLY: EXP_NAME, update_channels &
                                   , I_VERBOSE_LEVEL &
                                   , get_channel_info
    USE messy_main_channel_dimensions, ONLY: update_dimension_variable
    USE messy_main_timer,      ONLY: time_span_d   &
                                   , gregor2julian
#if defined(ICON) || defined(COSMO) || defined(MBM_RAD)
    USE messy_main_timer,      ONLY: MILLISECOND
#endif

    IMPLICIT NONE

    INTRINSIC :: ABS, ASSOCIATED, LEN, LEN_TRIM, NULL, REAL

    ! I/O
    INTEGER, INTENT(IN)         :: IOMODE  ! OUTPUT MODE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'channel_write_output_bi'
    INTEGER                     :: status
    CHARACTER(LEN=BASENAME_LEN) :: fnamebase = ''
    INTEGER                     :: i
    REAL(DP)                    :: yyyymmdd
    REAL(DP)                    :: now
    LOGICAL                     :: lexit
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: ptr  => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: gptr => NULL()
    INTEGER                               :: reprid
    INTEGER, SAVE    :: nrstcount     = 0
    CHARACTER(LEN=4) :: nrstcount_str = ''
    LOGICAL          :: lp
    INTEGER          :: YEAR_DATE, MONTH_DATE, DAY_DATE
    INTEGER          :: HOUR_DATE, MINUTE_DATE, SECOND_DATE
    INTEGER          :: MILLISECOND_DATE
    LOGICAL          :: l_clams = .FALSE.  ! is CLaMS running ?
    REAL(DP)         :: julsec  ! time in julian seconds
    INTEGER           :: js
    INTEGER          :: p_io_c ! IO PE per CHANNEL, required in bi_decompose

    IF (I_VERBOSE_LEVEL >= 2) &
         CALL start_message_bi(modstr,'WRITE OUTPUT',substr)

    SELECT CASE(IOMODE)
       !
    CASE(IOMODE_OUT)   ! ### ---------------- OUTPUT -------------------- ###
       !
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! SPECIAL HANDLING FOR ECHAM5-STREAM-ELEMENTS WITH laccu = .TRUE.
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(ECHAM5)
       CALL reset_accu_stream_elements(1)
#endif
#ifdef COSMO
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! SPECIAL HANDLING FOR FILTERED COSMO FIELDS
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL filter_cosmo_elements
#endif
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       CALL update_channels(status, 1, time_step_len) ! ACCUMULATE 2ndary
       CALL channel_halt(substr, status)
       !
       CALL update_channels(status, 2, time_step_len) ! PREPARE FOR OUTPUT
       CALL channel_halt(substr, status)
       !
       ! UPDATE FILENAME BASE
       fnamebase = EXP_NAME
       DO i = LEN_TRIM(EXP_NAME)+1, LEN(EXP_NAME)
          WRITE(fnamebase(i:i),'(a1)') '_'
       END DO

       YEAR_DATE   = YEAR_NEXT
       MONTH_DATE  = MONTH_NEXT
       DAY_DATE    = DAY_NEXT
       HOUR_DATE   = HOUR_NEXT
       MINUTE_DATE = MINUTE_NEXT
       SECOND_DATE = SECOND_NEXT
       MILLISECOND_DATE = MILLISECOND_NEXT

#if defined(ECHAM5) || defined(CESM1)
       WRITE(fnamebase(LEN(EXP_NAME)+1:),'(a1,i4,i2.2,i2.2,a1,i2.2,i2.2,a1)') &
            '_', YEAR_DATE, MONTH_DATE, DAY_DATE, '_' &
            , HOUR_DATE, MINUTE_DATE, '_'
#endif

#if defined(BLANK) || defined(MBM_CLAMS) || defined(MBM_MPIOM)
       WRITE(fnamebase(LEN(EXP_NAME)+1:)  &
            , '(a1,i4,i2.2,i2.2,a1,i2.2,i2.2,i2.2,a1)')&
            '_', YEAR_DATE, MONTH_DATE,  DAY_DATE, '_' &
            , HOUR_DATE, MINUTE_DATE, SECOND_DATE,'_'
#endif

#if defined(VERTICO)
       WRITE(fnamebase(LEN(EXP_NAME)+1:)  &
            , '(a1,i4,i2.2,i2.2,a1,i2.2,i2.2,i2.2,a1)')&
            '_', YEAR_DATE, MONTH_DATE, DAY_DATE, '_'
#endif

#if defined(MBM_RAD)
       YEAR_DATE   = YEAR
       MONTH_DATE  = MONTH
       DAY_DATE    = DAY
       HOUR_DATE   = HOUR
       MINUTE_DATE = MINUTE
       SECOND_DATE = SECOND
       MILLISECOND_DATE = MILLISECOND
       WRITE(fnamebase(LEN(EXP_NAME)+1:) &
            , '(a1,i4,i2.2,i2.2,a1,i2.2,i2.2,a1)') &
            '_', YEAR_DATE, MONTH_DATE, DAY_DATE, '_' &
            , HOUR_DATE, MINUTE_DATE, '_'
#endif

#if defined(COSMO)
       YEAR_DATE   = YEAR
       MONTH_DATE  = MONTH
       DAY_DATE    = DAY
       HOUR_DATE   = HOUR
       MINUTE_DATE = MINUTE
       SECOND_DATE = SECOND
       MILLISECOND_DATE = MILLISECOND
       WRITE(fnamebase(LEN(EXP_NAME)+1:)  &
            , '(a1,i4,i2.2,i2.2,a1,i2.2,i2.2,i2.2,a1)')&
            '_', YEAR_DATE, MONTH_DATE,  DAY_DATE, '_' &
            , HOUR_DATE, MINUTE_DATE, SECOND_DATE,'_'
#endif

#if defined (ICON)
       YEAR_DATE   = YEAR
       MONTH_DATE  = MONTH
       DAY_DATE    = DAY
       HOUR_DATE   = HOUR
       MINUTE_DATE = MINUTE
       SECOND_DATE = SECOND
       MILLISECOND_DATE = MILLISECOND
       WRITE(fnamebase(LEN(EXP_NAME)+1:)  &
            , '(a1,i4,i2.2,i2.2,a1,i2.2,i2.2,i2.2,i3.3,a1)')&
            '_', YEAR_DATE, MONTH_DATE,  DAY_DATE, '_' &
            , HOUR_DATE, MINUTE_DATE, SECOND_DATE, MILLISECOND_DATE,'_'
#endif

       ! UPDATE TIME STEP - DATA
       ! - TIME
       CALL time_span_d(now   &
            , YEAR_START, MONTH_START, DAY_START &
            , HOUR_START, MINUTE_START, SECOND_START, MILLISECOND_START  &
            , YEAR_DATE, MONTH_DATE, DAY_DATE         &
            , HOUR_DATE, MINUTE_DATE, SECOND_DATE, MILLISECOND_DATE)

       CALL update_dimension_variable(status, 'time', 'time', (/ now /))
       CALL channel_halt(substr, status)
       ! - YYYYMMDD
       yyyymmdd = REAL(ABS(YEAR_DATE)*10000 + MONTH_DATE*100 +DAY_DATE, DP) &
        + REAL((HOUR_DATE*3600 + MINUTE_DATE*60 + SECOND_DATE), DP)/86400.0_DP &
        + REAL(MILLISECOND_DATE, DP)/86400000.0_DP
       IF (YEAR<0) yyyymmdd = -yyyymmdd
       CALL update_dimension_variable(status, 'time', 'YYYYMMDD' &
            , (/ yyyymmdd /))
       CALL channel_halt(substr, status)
       ! - DT
       CALL update_dimension_variable(status, 'time', 'dt', (/ delta_time /))
       CALL channel_halt(substr, status)
       ! - TIME STEP
       CALL update_dimension_variable(status, 'time', 'nstep' &
            , (/ REAL(current_time_step, DP) /))
       CALL channel_halt(substr, status)
       !
       DO js = 1, NCHANNEL
          TIME_BNDS(js)%ptr(2) = now
       END DO

       CALL get_channel_info(status, 'clams')
       l_clams = (status == 0)
       if (l_clams) then
          julsec = REAL((HOUR*3600 + MINUTE*60 + SECOND),DP) + &         ! s
               ( (gregor2julian(YEAR, MONTH, DAY, 0, 0, 0) + 0.5_dp) - & ! ymd2jd()
               2451545.0_dp ) * &   ! 1.1.2000
               86400.0_dp
          julsec = julsec + delta_time
          CALL update_dimension_variable(status, 'time', 'JULIAN_SECONDS' &
               , (/ julsec /))
          CALL channel_halt(substr, status)
       endif
       !
    CASE(IOMODE_RST)   ! ### ---------------- RESTART -------------------- ###
       !
       ! FORCE NEW OUTPUT FILE IN NEXT STEP
       LFORCE_NEW_OUTPUT = .TRUE.
       ! UPDATE FILENAME BASE
       nrstcount = nrstcount + 1
       write(nrstcount_str,'(i4.4)') nrstcount
       ! UPDATE FILENAME BASE
       fnamebase = 'restart_'//nrstcount_str//'_'
       !
       ! SET RESTART ATTRIBUTES
       CALL initialize_restart_attributes(AMODE_WRITE)
       !
    END SELECT

    ! PREPARE OUTPUT / RESTART FILE
    ! NEW FILE, SAVE I/O-UNITS, FILE-IDs etc.
    ! p_pe is passed to determine I/O PE based on CHANNEL-settings, similar for
    ! all IO-routines
    CALL channel_init_io(status, IOMODE, fnamebase, AMODE_WRITE, p_pe)
    CALL channel_halt(substr, status)
    !
    ! WRITE HEADER
    SELECT CASE(IOMODE)
    CASE(IOMODE_OUT)
       CALL channel_write_header(status, IOMODE, DIMID_TIME, p_pe)
    CASE(IOMODE_RST)
       CALL channel_write_header(status, IOMODE, DIMID_TIME, p_pe &
            , restart_att)
    END SELECT
    CALL channel_halt(substr, status)
    !
    ! WRITE TIME STEP DATA TO OUTPUT FILE
    CALL channel_write_time(status, IOMODE, DIMID_TIME, p_pe)
    CALL channel_halt(substr, status)
    ! WRITE DATA
    ! p_io_c is passed from channel_write_data to bi_decompose to determine PE
    ! for gathering fields (IO PE)
    DO
       ! NEXT POINTER
       CALL channel_write_data(status, lp  &
            , IOMODE, lexit, ptr, reprid, p_io_c, p_pe)
       CALL channel_halt(substr, status)
       IF (lexit) EXIT

       IF (lp) THEN ! parallel I/O
          CALL bi_vector(status, -1, reprid, gptr, ptr)
          IF (status /= 0) &
               CALL error_bi('bi_vector reported an error', substr)
       ELSE
          CALL bi_decompose(status, -1, reprid, gptr, ptr, p_io_c)
          IF (status /= 0) &
               CALL error_bi('bi_decompose reported an error',substr)
       END IF
       ! OUTPUT
       ! NOTE: gptr only associated on I/O-PE
       CALL channel_write_data(status, lp &
            , IOMODE, lexit, gptr, reprid, p_io_c, p_pe)
       CALL channel_halt(substr, status)
    END DO
    ! FLUSH (OUTPUT) / CLOSE (RESTART) BUFFER
    CALL channel_finish_io(status, IOMODE, (IOMODE==IOMODE_RST),  p_pe)
    CALL channel_halt(substr, status)

    ! CLEAN MEMORY
    IF (ASSOCIATED(gptr)) THEN
       DEALLOCATE(gptr)
       NULLIFY(gptr)
    END IF

    IF (IOMODE == IOMODE_OUT) THEN
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! SPECIAL HANDLING FOR ECHAM5-STREAM-ELEMENTS WITH laccu = .TRUE.
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(ECHAM5)
       CALL reset_accu_stream_elements(2)
#endif
       DO js = 1, NCHANNEL
          IF (LOUTPUT_NOW(js)) TIME_BNDS(js)%ptr(:) = now
       END DO
       !
       CALL update_channels(status, 3, time_step_len) ! RESET AFTER OUTPUT
       CALL channel_halt(substr, status)
       !
    END IF

    IF (I_VERBOSE_LEVEL >= 2) &
         CALL end_message_bi(modstr,'WRITE OUTPUT',substr)

  END SUBROUTINE channel_write_output_bi
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_init_restart_bi

    USE messy_main_mpi_bi,  ONLY: p_io, p_pe, p_all_comm, p_nprocs
    USE messy_main_timer,   ONLY: lresume
    USE messy_main_channel_io, ONLY: initialize_parallel_io

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'channel_init_restart_bi'
    INTEGER :: status

    IF (.NOT. lresume) RETURN

    CALL initialize_parallel_io(status, p_pe, p_io, p_all_comm, p_nprocs)
    CALL channel_halt(substr, status)

    CALL initialize_restart_attributes(AMODE_INIT)

  END SUBROUTINE channel_init_restart_bi
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_read_restart(chname, dom_id)

    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast, p_pe
    USE messy_main_blather_bi, ONLY: error_bi

    USE messy_main_channel_io, ONLY: channel_init_io   &
                                   , channel_read_data &
                                   , channel_finish_io
    USE messy_main_channel,    ONLY: t_channel_list, t_channel, GCHANNELLIST &
                                   , l_cheat
#ifdef ICON
    USE messy_main_channel,    ONLY: STRLEN_OBJECT
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, NULL

    ! I/O
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: chname
    INTEGER,          INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'main_channel_read_restart'
    INTEGER                      :: status
    CHARACTER(LEN=BASENAME_LEN)  :: fnamebase = ''
    TYPE(t_channel_list),  POINTER        :: ls
    TYPE(t_channel),       POINTER        :: channel
    LOGICAL                               :: lexit
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: ptr  => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: gptr => NULL()
    INTEGER                               :: reprid
    LOGICAL                               :: lok
    LOGICAL                               :: lp
    INTEGER                               :: ife
#ifdef ICON
    CHARACTER(LEN=34+STRLEN_CHANNEL+8)    :: filename = ''
    CHARACTER(LEN=STRLEN_OBJECT+4)        :: objname = ''
    LOGICAL                               :: lignore  = .FALSE.
    LOGICAL                               :: lrestreq = .FALSE.
#endif

    CALL start_message_bi(modstr,'READ RESTART',substr)

    IF (PRESENT(chname)) THEN
       IF (p_parallel_io) WRITE(*,*) 'chname: ', TRIM(chname)
       IF (TRIM(chname) == 'tracer_gp') l_cheat = .TRUE.
    ENDIF

    ! UPDATE FILENAME BASE
    fnamebase = 'restart_'

    ! INITIALIZE RESTART ATTRIBUTES
    CALL initialize_restart_attributes(AMODE_READ)

    ! OPEN RESTART FILES AND CHECK HEADER INFORMATION
    ! OPEN FILE FOR READ
    CALL channel_init_io(status,  IOMODE_RST &
         ,  fnamebase, AMODE_READ, p_pe, restart_att, chname, dom_id &
#ifdef ICON
         , lclose = .TRUE. &
#endif
         , llp_io = p_parallel_io)
    CALL channel_halt(substr, status)

    ! BROADCAST %tslo IN ALL CHANNELS
    ! this should only be required for non-parallel I/O, however ...
    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT
       channel => ls%this
       ! ------------------------------------------------
       CALL p_bcast(channel%int%tslo, p_io)
       ! ------------------------------------------------
       ! broadcast deleted filename for missing restart files
       IF (p_parallel_io) THEN
          IF (TRIM(channel%int%fname(IOMODE_RST)) == '') THEN
             ife = 0
          ELSE
             ife = 1
          END IF
       END IF
       CALL p_bcast(ife, p_io)
       IF (ife == 0) channel%int%fname(IOMODE_RST) = ' '
       ls => ls%next
    END DO channel_loop

    ! READ DATA
    DO
       ! NEXT POINTER
       CALL channel_read_data(status, p_parallel_io &
            , IOMODE_RST, lexit, gptr, reprid, lp, chname, dom_id &
#ifdef ICON
            , lskipinput=.TRUE., filename=filename, objname=objname &
            , lrestreq = lrestreq,lignore = lignore                 &
#endif
            )
       CALL channel_halt(substr, status)
       IF (lexit) EXIT
       !
#ifndef ICON
       IF (lp) THEN ! parallel I/O
          IF (ASSOCIATED(gptr)) THEN
             CALL bi_vector(status, 1, reprid, gptr, ptr)
             IF (status /= 0) CALL error_bi( &
                  'bi_vector reported an error', substr)
          END IF
       ELSE
          ! NOTE: gptr only associated on I/O-PE
          IF (p_parallel_io) lok = ASSOCIATED(gptr)
          CALL p_bcast(lok, p_io)
          !
          IF (lok) THEN
             CALL bi_decompose(status, 1, reprid, gptr, ptr, p_io)
             IF (status /= 0) CALL error_bi( &
                  'bi_decompose reported an error',substr)
             ! NOTE: ptr now associated on all PEs
             !       or =>NULL if .NOT. lok
          END IF
       END IF
       !
#else
          CALL channel_readdistrib_ICON_data(status, ptr, reprid &
               , filename, objname, lrestreq, lignore)
#endif
       ! DISTRIBUTE
       CALL channel_read_data(status, p_parallel_io, &
            IOMODE_RST, lexit, ptr, reprid, lp       &
            )
       CALL channel_halt(substr, status)
       !
       ! RESET MEMORY
       IF (ASSOCIATED(ptr)) THEN
          DEALLOCATE(ptr)
          NULLIFY(ptr)
       END IF
       IF (ASSOCIATED(gptr)) THEN
          DEALLOCATE(gptr)
          NULLIFY(gptr)
       END IF
    END DO

#ifndef ICON
    ! CLOSE ALL RESTART FILES
    CALL channel_finish_io(status, IOMODE_RST, .TRUE., p_pe, chname &
         , dom_id, llp_io=p_parallel_io)
    CALL channel_halt(substr, status)
#endif

    l_cheat = .FALSE.

    CALL end_message_bi(modstr,'READ RESTART',substr)

#ifdef ICON

  CONTAINS

    ! ---------------------------------------------------------------------
    SUBROUTINE channel_readdistrib_ICON_data(status         &
         , ptr, reprid, filename, objname, lrestreq, lignore)

      USE messy_main_channel_repr,  ONLY: t_representation, get_representation &
                                        , repr_reorder                         &
                                        , repr_rank_location, IRANK
      USE messy_main_channel,       ONLY: get_channel_object, STRLEN_OBJECT
      USE messy_main_mpi_bi,        ONLY: p_pe, p_io, p_bcast
      USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, INT_UNDEF
      USE messy_main_bmluse_bi,     ONLY: wp, p_patch                     &
                                        , openInputFile, closeFile        &
                                        , t_stream_id                     &
                                        , on_cells, on_edges, on_vertices &
                                        , read_netcdf_broadcast_method    &
                                        , read_netcdf_distribute_method   &
                                        , read_3D_extdim, read_2d_extdim  &
                                        , read_3D, read_2d                &
                                        , read_bcast_REAL_2D              &
                                        , read_1D, read_0D_real           &
                                        , read_inq_varexists

      IMPLICIT NONE
      INTRINSIC :: TRIM

      INTEGER,                            INTENT(OUT):: status
      REAL(DP), DIMENSION(:,:,:,:), POINTER          :: ptr
      INTEGER,                            INTENT(IN) :: reprid
      CHARACTER(LEN=34+STRLEN_CHANNEL+8), INTENT(IN) :: filename
      CHARACTER(LEN=STRLEN_OBJECT+4),     INTENT(IN) :: objname
      LOGICAL,                            INTENT(IN) :: lrestreq
      LOGICAL,                            INTENT(IN) :: lignore

      ! LOCAL
      CHARACTER(LEN=*), PARAMETER :: substr = "channel_readdistrib_ICON_data"
      INTEGER                                :: iz
      TYPE(t_stream_id), SAVE                :: fileID
      INTEGER                                :: location
      INTEGER ,  SAVE                        :: ifileID
      INTEGER                                :: pID
      INTEGER                                :: read_method
      REAL(wp),                     POINTER  :: lptr0d => NULL()
      REAL(wp), DIMENSION(:),       POINTER  :: lptr1d => NULL()
      REAL(wp), DIMENSION(:,:),     POINTER  :: lptr2d => NULL()
      REAL(wp), DIMENSION(:,:,:),   POINTER  :: lptr3d => NULL()
      REAL(wp), DIMENSION(:,:,:,:), POINTER  ::   lptr => NULL()
      REAL(wp), DIMENSION(:,:,:,:), POINTER  ::   zptr => NULL()
      REAL(wp), DIMENSION(:,:),     POINTER  :: zptr2d => NULL()

      TYPE(t_representation),       POINTER  :: repr
      INTEGER                                :: istat
      INTEGER                                :: rankidx(IRANK)
      LOGICAL                                :: ret

      ! CLEAN POINTER
      IF (ASSOCIATED(ptr)) THEN
         DEALLOCATE(ptr)
         NULLIFY(ptr)
      END IF

      IF (TRIM(filename) == '') THEN
         status = 0
         RETURN
      END IF

      ! DETERMINE READ METHOD from REPRID
      CALL get_representation(status, reprid, repr)
      CALL channel_halt(substr,status)

      pID = repr%dom_id

      SELECT CASE(repr%dctype)
      CASE(DC_BC)
         read_method = read_netcdf_broadcast_method
         IF (repr%rank > 2) CALL error_bi(&
              'restart broadcast for rank > 2 not implemented', substr)
      CASE(DC_GP)
         read_method = read_netcdf_distribute_method
         IF (repr%rank >= 2) THEN
            ! DETERMINE location: cells, vertices or edges
            IF (repr%dim(1)%ptr%id == DIMID_NCELLS(pID)) THEN
               location = on_cells
            ELSE IF (repr%dim(1)%ptr%id == DIMID_NEDGES(pID)) THEN
               location = on_edges
            ELSE IF (repr%dim(1)%ptr%id == DIMID_NVERTS(pID)) THEN
               location = on_vertices
            ELSE
               CALL error_bi("LOCATION NOT KNOWN for Repr: "//TRIM(repr%name) &
                    , substr)
            END IF
         ELSE
            CALL error_bi("DC_GP for rank < 2 does not make sense: " &
                 , substr)
         END IF
      CASE (DC_AG)
         read_method = read_netcdf_broadcast_method
      CASE (DC_TRIX)
         read_method = read_netcdf_distribute_method
      CASE DEFAULT
         CALL error_bi("DECOMPOSITION METHOD NOT YET IMPLEMENTED", substr)
         RETURN
      END SELECT

      ! OPEN FILE
      SELECT CASE(repr%rank)
      CASE (0)
         ifileID = openInputFile(TRIM(filename))
         ret = read_inq_varexists(ifileID, TRIM(objname))
         IF (read_inq_varexists(ifileID, TRIM(objname))) THEN
            ALLOCATE(ptr(1,1,1,1))
            ptr(1,1,1,1) = &
                 read_0D_real(ifileID, variable_name=TRIM(objname))
         ELSE
            CALL handle_status(lrestreq, lignore)
         END IF
         CALL closeFile(ifileID)
      CASE(1)
         SELECT CASE (repr%dctype)

         CASE (DC_AG)
            read_method = read_netcdf_broadcast_method
            ifileID = openInputFile(TRIM(filename))
            IF (read_inq_varexists(ifileID, TRIM(objname))) THEN
               CALL read_1D(ifileID, variable_name=TRIM(objname) &
                    , return_pointer=lptr1d)
               ALLOCATE(ptr(repr%ldimlen(repr%order_mem2out(1))  &
                    , repr%ldimlen(repr%order_mem2out(2))        &
                    , repr%ldimlen(repr%order_mem2out(3))        &
                    , repr%ldimlen(repr%order_mem2out(4))))
               SELECT CASE (repr%link)
               CASE ('---x')
                  !qqq this fits for QTIMER case, could be different ..
                  ptr(1,1,1,1)=lptr1d(p_pe+1)
               CASE DEFAULT
                  CALL error_bi('rank1 DC_AG, link not implemented'//repr%link &
                       , substr)
               END SELECT
               NULLIFY(lptr1d)
            ELSE
               CALL handle_status(lrestreq, lignore)
            ENDIF
            CALL closeFile(ifileID)

         CASE DEFAULT
            ifileID = openInputFile(TRIM(filename))
            ret= read_inq_varexists(ifileID,TRIM(objname))
            IF (read_inq_varexists(ifileID, TRIM(objname))) THEN
               CALL read_1D(ifileID, variable_name=TRIM(objname) &
                    , return_pointer=lptr1d)
               ALLOCATE(ptr(repr%ldimlen(repr%order_mem2out(1))  &
                    , repr%ldimlen(repr%order_mem2out(2))        &
                    , repr%ldimlen(repr%order_mem2out(3))        &
                    , repr%ldimlen(repr%order_mem2out(4))))
               SELECT CASE (repr%link)
               CASE ('x---')
                  ptr(:,1,1,1)=lptr1d(:)
               CASE ('-x--')
                  ptr(1,:,1,1)=lptr1d(:)
               CASE ('--x-')
                  ptr(1,1,:,1)=lptr1d(:)
               CASE ('---x')
                  ptr(1,1,1,:)=lptr1d(:)
               CASE DEFAULT
                  CALL error_bi('rank1 DC_BC, wrong link: '//repr%link &
                       , substr)
               END SELECT
               NULLIFY(lptr1d)
            ELSE
               CALL handle_status(lrestreq, lignore)
            ENDIF
            CALL closeFile(ifileID)
         END SELECT

      CASE(2)
         SELECT CASE (repr%dctype)

         CASE (DC_AG)
            read_method = read_netcdf_broadcast_method
            ifileID = openInputFile(TRIM(filename))
            IF (read_inq_varexists(ifileID, TRIM(objname))) THEN
               CALL read_bcast_REAL_2D(ifileID, TRIM(objname) &
                    , return_pointer=lptr2d)

               ALLOCATE(ptr(repr%ldimlen(repr%order_mem2out(1))  &
                    , repr%ldimlen(repr%order_mem2out(2))        &
                    , repr%ldimlen(repr%order_mem2out(3))        &
                    , repr%ldimlen(repr%order_mem2out(4))))
               SELECT CASE (repr%link)
               CASE ('x--x')
                  !qqq this fits for RND case, could be different ..
                  ptr(:,1,1,1)=lptr2d(:,p_pe+1)
               CASE DEFAULT
                  CALL error_bi('rank2 DC_AG, link not implemented'//repr%link &
                       , substr)
               END SELECT
               NULLIFY(lptr2d)
            ELSE
               CALL handle_status(lrestreq, lignore)
            ENDIF
            CALL closeFile(ifileID)

         CASE (DC_GP)
            read_method = read_netcdf_distribute_method
            fileID = openInputFile(TRIM(filename), p_patch(pID), read_method)

            IF (read_inq_varexists(fileID, TRIM(objname))) THEN
               ALLOCATE(zptr(repr%ldimlen(repr%order_mem2out(1)) &
                    , repr%ldimlen(repr%order_mem2out(2))        &
                    , repr%ldimlen(repr%order_mem2out(3))        &
                    , repr%ldimlen(repr%order_mem2out(4))))
               zptr(:,:,:,:) = 0._dp

               CALL read_2D(fileID, location, TRIM(objname) &
                    , return_pointer=lptr2d)
               zptr(:,:,1,1) = lptr2d(:,:)

               IF (p_parallel_io) write (*,*) 'OBJECT: ', TRIM(objname)&
                    ,' read from file: ', TRIM(filename)

               CALL repr_reorder(status, -1, .TRUE., repr, ptr, zptr)
               IF (status /=0) CALL error_bi('REORDER 2D ERROR', substr)

               NULLIFY (lptr2d)
               DEALLOCATE(zptr); NULLIFY(zptr)
            ELSE
               CALL handle_status(lrestreq, lignore)
            END IF
            CALL closeFile(fileID)

         CASE DEFAULT
           CALL error_bi("DECOMPOSITION METHOD NOT YET IMPLEMENTED FOR rank 2" &
                 , substr)
         END SELECT
      CASE(3)
         SELECT CASE(repr%dctype)
         CASE(DC_GP)
            read_method = read_netcdf_distribute_method

            fileID = openInputFile(TRIM(filename), p_patch(pID), read_method)

            IF (read_inq_varexists(fileID, TRIM(objname))) THEN
               ! ALLOCATE return pointer dependent on representation
               ! NOTE: We have to DISTINGUISH BETWEEN 'XZY-' and 'XYN-' fields
               CALL  repr_rank_location(repr%axis, rankidx)

               ! ERROR IF XY are not defined
               IF (rankidx(1) == INT_UNDEF .OR.  rankidx(2) == INT_UNDEF) THEN
                  CALL error_bi (&
                       'restart field cannot be read by ICON routines'&
                       , substr)
               ELSE IF (rankidx(3) /= INT_UNDEF) THEN
                  ! FIELD OF AXIS TYPE 'XZY-'
                  CALL read_3D(fileID, location, objname &
                       , return_pointer=lptr3d &
                       , levelsDimName=TRIM(repr%dim(2)%ptr%name))
               ELSE IF  (rankidx(4) /= INT_UNDEF) THEN
                  ! FIELD OF AXIS TYPE 'XYN'
                  CALL read_2D_extdim(fileID, location, TRIM(objname)         &
                       , return_pointer=lptr3d                                &
                       , start_extdim=1,end_extdim=repr%dim(rankidx(4))%ptr%len&
                       , extdim_name=TRIM(repr%dim(rankidx(4))%ptr%name))
               ENDIF

               IF (p_parallel_io) write (*,*) 'OBJECT :', TRIM(objname)&
                    ,' read from file: ', TRIM(filename)

               ALLOCATE(ptr(repr%ldimlen(1), repr%ldimlen(2)     &
                    , repr%ldimlen(3), repr%ldimlen(4)))
               ptr = 0._dp
               ptr(:,:,:,1) = lptr3d(:,:,:)

               DEALLOCATE(lptr3d)
               NULLIFY(lptr3d)
            ELSE
               CALL handle_status(lrestreq, lignore)
            END IF
            CALL closeFile(fileID)
            CASE DEFAULT
               CALL error_bi(&
                    "DECOMPOSITION METHOD NOT YET IMPLEMENTED FOR RANK 3"&
                    , substr)
            END SELECT
      CASE(4)

         read_method = read_netcdf_distribute_method
         fileID = openInputFile(TRIM(filename), p_patch(pID), read_method)

         IF (read_inq_varexists(fileID, TRIM(objname))) THEN
            ! NOTE: ICON read routine assumes extdim in third place
            !       level on second place
            !  currently we assume the same for MESSy /could be changed if
            !  necessary
            CALL read_3D_extdim(fileID, location, TRIM(objname)          &
                 , return_pointer=lptr                                   &
                 , start_extdim=1, end_extdim= repr%dim(3)%ptr%len       &
                 , levelsDimName=TRIM(repr%dim(2)%ptr%name)              &
                 , extdim_name=TRIM(repr%dim(3)%ptr%name))

            ALLOCATE(ptr(repr%ldimlen(1), repr%ldimlen(2)     &
                 , repr%ldimlen(3), repr%ldimlen(4)))
            ptr = lptr
            NULLIFY(lptr)
         ELSE
            CALL handle_status(lrestreq, lignore)
         END IF
         CALL closeFile(fileID)

      CASE DEFAULT
         write (*,*) substr, 'ERROR: rank: ', TRIM(objname),' ',repr%rank
         CALL error_bi("RANK out of range / not yet implemented", substr)
      END SELECT

    END SUBROUTINE channel_readdistrib_ICON_data

    ! ----------------------------------------------------------------------

  SUBROUTINE handle_status(lreq, lign)

    USE messy_main_channel, ONLY: I_VERBOSE_LEVEL

    IMPLICIT NONE

    LOGICAL, INTENT(IN)  :: lreq
    LOGICAL, INTENT(IN)  :: lign

    IF (lreq) THEN
       IF (lign) THEN
          IF (I_VERBOSE_LEVEL >= 1 .AND. p_parallel_io) &
               WRITE(*,*) substr &
               ,'    WARNING: REQUIRED RESTART VARIABLE ', &
               'NOT PRESENT! ' &
               ,'HOWEVER: IGNORE = T FOR THIS OBJECT'
       ELSE
          IF (p_parallel_io) &
               write (*,*) 'OBJECT :', TRIM(objname)&
               ,' could not be read from file: ', TRIM(FILENAME)
          !          stat = 3211 ! RESTART VARIABLE REQUIRED BUT NOT PRESENT
          CALL channel_halt(substr,3211)
       END IF
    ELSE
       IF (I_VERBOSE_LEVEL >= 1 .AND. p_parallel_io) &
            WRITE(*,*) '  WARNING: OBJECT: '//TRIM(objname)//&
            &' NOT PRESENT IN RESTART FILE '//TRIM(filename)
    END IF

  END SUBROUTINE handle_status

#endif

  END SUBROUTINE main_channel_read_restart
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PUBLIC HELPER ROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE p_bcast_attribute(att, p)

    USE messy_main_channel_attributes, ONLY: t_attribute
    USE messy_main_mpi_bi,             ONLY: p_bcast

    IMPLICIT NONE

    ! I/O
    TYPE(t_attribute) :: att
    INTEGER           :: p

    CALL p_bcast(att%name, p)
    CALL p_bcast(att%type, p)
    CALL p_bcast(att%iflag, p)
    CALL p_bcast(att%i, p)
    CALL p_bcast(att%c, p)
    CALL p_bcast(att%r, p)

  END SUBROUTINE p_bcast_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE (PART II)
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE initialize_restart_attributes(AMODE)

    USE messy_main_data_bi,    ONLY: basemodstr => modstr
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_timer_bi,   ONLY: get_time_step => timer_get_time_step &
                                   , messy_timer_init_manager             &
                                   , timer_message
    USE messy_main_timer,      ONLY: INIT_STEP                              &
                                   , timer_get_date,timer_set_date          &
                                   , delta_time                             &
                                   , YEAR_START,MONTH_START,DAY_START       &
                                   , HOUR_START, MINUTE_START, SECOND_START
    USE messy_main_channel_io, ONLY: channel_init_restart
    USE messy_main_channel,    ONLY: new_attribute, get_attribute &
                                   , AF_RST_CMP, AF_RST_INP, AF_RST_NONE
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM

#if defined(ECHAM5)
    USE mo_time_control,     ONLY: ec_manager_init &
                                 , E5_resume_date => resume_date &
                                 , E5_start_date  => start_date &
                                 , lresume
    USE mo_time_conversion,  ONLY: TC_set, time_native, TC_convert, time_intern
    USE messy_main_timer,    ONLY: timer_get_lresume
#endif

#ifdef COSMO
    ! COSMO
    USE data_modelconfig,           ONLY: dt
    USE data_runcontrol,            ONLY: nnow, nnew, nold    &
                                        , nbd1, nbd2          &
                                        , hnextrad, nextrad, ntke
    ! MESSy/BMIL
    USE messy_main_data_bi,         ONLY: L_IS_CHILD
    USE messy_main_grid_def_mem_bi, ONLY: vcoord, refatm, nfltvc
#endif

#ifdef ICON
    USE messy_main_timer,      ONLY: MILLISECOND_START
#endif

    IMPLICIT NONE
    INTRINSIC :: INT, MOD, SIGN, TRIM

    ! I/O
    INTEGER, INTENT(IN) :: AMODE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='initialize_restart_attributes'
    CHARACTER(LEN=DATE_TIME_STR_LEN) :: str_date_time = ''
    !
    INTEGER           :: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
    INTEGER           :: MILLISECOND
    INTEGER           :: iflag
    INTEGER           :: iflag_dt
    INTEGER           :: yyyymmdd
    REAL(DP)          :: hhmmss
    ! qqq why SAVE?
    INTEGER, SAVE     :: nstep
    REAL(DP)          :: timestep
    INTEGER           :: status
    LOGICAL           :: lp = .FALSE.
    INTEGER           :: kyr, kmo, kdy, khr, kmn, kse, iymd
    INTEGER           :: kms
    CHARACTER(LEN=STRLEN_MEDIUM) :: zdstr = ''
#if defined(ECHAM5)
    TYPE(time_intern) :: io_date
#endif

#if defined(ECHAM5)
    ! set lresume of original ECHAM5 timer
    CALL timer_get_lresume(lresume)
#endif

    SELECT CASE(AMODE)
       CASE(AMODE_WRITE)
          iflag    = AF_RST_NONE
          iflag_dt = AF_RST_NONE
#ifndef ICON
          zdstr    = 'next'
#else
          ! 'current' is required in ICON in contrast to all other models, as
          ! timer_time(2) is called at the beginning and not at the end of the
          ! time loop
          zdstr    = 'current'
#endif
          nstep    = get_time_step()
       CASE(AMODE_READ)
          iflag    = AF_RST_CMP
          iflag_dt = AF_RST_NONE
          zdstr    = 'resume'
       CASE(AMODE_INIT)
          iflag    = AF_RST_INP
          iflag_dt = AF_RST_INP
          ! use start-date here as restart date ist not yet defined (for TIMER)
          ! and will be overwritten anyway
          zdstr = 'start'
          nstep    = INIT_STEP
    END SELECT

    ! NOTE: The time information in *_START for the MESSY-TIMER must be
    !       set here, therfore the local variables (YEAR, MONTH, ...)
    !       cannot be used.
#ifndef ICON
    WRITE(str_date_time,'(i4.4,i2.2,i2.2,a1,3(i2.2))') &
          YEAR_START, MONTH_START,  DAY_START,' '    &
          ,HOUR_START, MINUTE_START, SECOND_START
#else
    WRITE(str_date_time,'(i4,i2.2,i2.2,a1,3(i2.2),a1,i3.3)') &
          YEAR_START, MONTH_START,  DAY_START,' '    &
          ,HOUR_START, MINUTE_START, SECOND_START, '.', MILLISECOND_START
#endif
    !
    CALL new_attribute(status, restart_att &
         , 'start_date_time', c=TRIM(str_date_time) &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)

    ! - RESTART DATE AND TIME
    IF (AMODE /= AMODE_INIT) THEN
       CALL timer_get_date(status, TRIM(zdstr)  &
            ,YEAR,MONTH,DAY,HOUR,MINUTE,SECOND,MILLISECOND)
       CALL timer_message(status, substr)
    ELSE
       ! in case of restart start_date is not initialised, use
       ! start date components directly
       YEAR        =   YEAR_START
       MONTH       =  MONTH_START
       DAY         =    DAY_START
       HOUR        =   HOUR_START
       MINUTE      = MINUTE_START
       SECOND      = SECOND_START
#ifdef ICON
       MILLISECOND = MILLISECOND_START
#endif
    END IF
#ifndef ICON

    WRITE(str_date_time,'(i4.4,i2.2,i2.2,a1,3(i2.2))') &
         YEAR, MONTH, DAY,' ',HOUR, MINUTE, SECOND
#else
    WRITE(str_date_time,'(i4.4,i2.2,i2.2,a1,3(i2.2),a1,i3.3)') &
         YEAR, MONTH, DAY,' ',HOUR, MINUTE, SECOND, '.', MILLISECOND
#endif
    !
    CALL new_attribute(status, restart_att &
         , 'restart_date_time', c=TRIM(str_date_time) &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)

    ! - CURRENT TIME STEP
    CALL new_attribute(status, restart_att &
         , 'nstep', i=nstep                &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)

    ! - TIME STEP LENGTH
    CALL new_attribute(status, restart_att &
         , 'timestep', r=delta_time        &
         , loverwrite=.TRUE., iflag = iflag_dt)
    CALL channel_halt(substr, status)

#ifdef COSMO
    ! - counter for radiation calculation
    CALL new_attribute(status, restart_att &
         , 'hnextrad', r=hnextrad   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    !
    ! - counter for radiation calculation
    CALL new_attribute(status, restart_att &
         , 'nextrad', i=nextrad   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    !
    ! - COORDINATE TYPE
    CALL new_attribute(status, restart_att &
         , 'irefatm', i=refatm%irefatm   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    !
    ! - COORDINATE TYPE
    CALL new_attribute(status, restart_att &
         , 'ivctype', i=vcoord%ivctype   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    !
    ! - current index of time_levels
    CALL new_attribute(status, restart_att &
         , 'nnow', i=nnow   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, restart_att &
         , 'nnew', i=nnew   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, restart_att &
         , 'nold', i=nold   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, restart_att &
         , 'nfltvc', i=nfltvc   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, restart_att &
         , 'ntke', i=ntke   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)

    IF (L_IS_CHILD) THEN
       CALL new_attribute(status, restart_att &
            , 'nbd1', i=nbd1   &
            , loverwrite=.TRUE., iflag = iflag)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, restart_att &
            , 'nbd2', i=nbd2   &
            , loverwrite=.TRUE., iflag = iflag)
       CALL channel_halt(substr, status)
    END IF
#endif

    ! CONTINUE ONLY DURING INITIALIZATION AFTER RESTART
    IF (AMODE /= AMODE_INIT) RETURN

#if defined(ECHAM5)
    CALL channel_init_restart(status, lp, p_parallel_io &
         , 'restart_g1a', restart_att)
#elif defined(COSMO)
    CALL channel_init_restart(status, lp, p_parallel_io &
         , 'restart_COSMO_ORI', restart_att)
#elif defined(CESM1) || defined(MBM_MPIOM)
    CALL channel_init_restart(status, lp, p_parallel_io &
         , 'restart_'//TRIM(basemodstr), restart_att)
#elif defined(ICON)
    ! at least a restart file for patch D01 should be available ...
    CALL channel_init_restart(status, lp, p_parallel_io &
         , 'restart_tracer_gp_D01', restart_att)
#else
    ! BLANK, MBM_CLAMS, VERTICO
    CALL channel_init_restart(status, lp, .TRUE. &
         , 'restart_'//TRIM(basemodstr), restart_att)
#endif
    CALL channel_halt(substr, status)

    ! READ DATA FROM ATTRIBUTES:
    ! - START DATE AND TIME
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'start_date_time' &
            , c=str_date_time)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(str_date_time, p_io)
    READ(str_date_time,*) yyyymmdd, hhmmss
#if defined(ECHAM5)
    CALL TC_set (yyyymmdd, INT(hhmmss), io_date)
    CALL TC_convert(io_date, E5_start_date)
#endif

    iymd = SIGN(1,INT(yyyymmdd))*INT(yyyymmdd)
    kyr = INT(iymd/10000)
    kmo = MOD(iymd,10000)/100
    kdy = MOD(iymd,100)
    khr = INT(hhmmss)/10000
    kmn = MOD(INT(hhmmss),10000)/100
    kse = MOD(INT(hhmmss),100)
    kms = INT((hhmmss - REAL(INT(hhmmss),dp)) * 1000._dp)
    CALL timer_set_date(status,'start',kyr, kmo, kdy, khr, kmn, kse, kms)

    ! - RESTART DATE
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'restart_date_time' &
            , c=str_date_time)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(str_date_time, p_io)
    READ(str_date_time,*) yyyymmdd, hhmmss
#if defined(ECHAM5)
    CALL TC_set (yyyymmdd, INT(hhmmss), io_date)
    CALL TC_convert(io_date, E5_resume_date)
#endif

    iymd = SIGN(1,INT(yyyymmdd))*INT(yyyymmdd)
    kyr = INT(iymd/10000)
    kmo = MOD(iymd,10000)/100
    kdy = MOD(iymd,100)
    khr = INT(hhmmss)/10000
    kmn = MOD(INT(hhmmss),10000)/100
    kse = MOD(INT(hhmmss),100)
    kms = INT((hhmmss - REAL(INT(hhmmss),dp)) * 1000._dp)
    CALL timer_set_date(status,'resume',kyr, kmo, kdy, khr, kmn, kse, kms)

    ! - CURRENT TIME STEP
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'nstep' &
            , i=nstep)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(nstep, p_io)

    ! - TIME STEP LENGTH
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'timestep' &
            , r=timestep)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(timestep, p_io)

#ifdef COSMO
    IF (INT(timestep) /= INT(dt)) THEN
       IF (p_parallel_io) THEN
          write (*,*) 'COSMO TIME STEP (',dt &
               ,') DOES NOT MATCH RESTART TIME STEP (',timestep,')'
          write (*,*) 'CHANGE OF TIMESTEP NOT POSSIBLE IN COSMO'
          CALL channel_halt(substr,3210)
       ENDIF
    ENDIF
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'hnextrad' &
            , r=hnextrad)
       CALL channel_halt(substr, status)
    END IF

    ! radiation trigger
    IF (.NOT. lp) CALL p_bcast(hnextrad, p_io)

    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'nextrad' &
            , i=nextrad)
       CALL channel_halt(substr, status)
    END IF

    ! radiation trigger
    IF (.NOT. lp) CALL p_bcast(nextrad, p_io)

    ! - coordinate type
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'ivctype' &
            , i=vcoord%ivctype)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(vcoord%ivctype, p_io)
    ! - coordinate type
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'irefatm' &
            , i=refatm%irefatm)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(refatm%irefatm, p_io)

    ! - coordinate type
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'nfltvc' &
            , i=nfltvc)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(nfltvc, p_io)

    ! the indices for the timelevels nnow, nnew, nold are rotating
    ! thus the information about the values of these three varaibles
    ! must be saved in the restart files in check in "init_restart"

    ! GET INDEX OF TIME LEVEL NNOW
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'nnow' &
            , i=nnow)
       CALL channel_halt(substr, status)
    END IF
    IF (.NOT. lp) CALL p_bcast(nnow, p_io)

    ! GET INDEX OF TIME LEVEL NNEW
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'nnew' &
            , i=nnew)
       CALL channel_halt(substr, status)
    END IF
    IF (.NOT. lp) CALL p_bcast(nnew, p_io)

    ! GET INDEX OF TIME LEVEL NOLD
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'nold' &
            , i=nold)
       CALL channel_halt(substr, status)
    END IF
    IF (.NOT. lp) CALL p_bcast(nold, p_io)

    ! GET INDEX OF TKE TIME LEVEL NTKE
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'ntke' &
            , i=ntke)
       CALL channel_halt(substr, status)
    END IF
    IF (.NOT. lp) CALL p_bcast(ntke, p_io)

    IF (L_IS_CHILD) THEN
       IF (p_parallel_io .OR. lp) THEN
          CALL get_attribute(status, restart_att, 'nbd1' &
               , i=nbd1)
          CALL channel_halt(substr, status)
       END IF
       IF (.NOT. lp) CALL p_bcast(nbd1, p_io)
       IF (p_parallel_io .OR. lp) THEN
          CALL get_attribute(status, restart_att, 'nbd2' &
               , i=nbd2)
          CALL channel_halt(substr, status)
       END IF
       IF (.NOT. lp) CALL p_bcast(nbd2, p_io)
    END IF

#endif

#if defined(ECHAM5)
    CALL timer_sync
    ! (RE-)INITIALIZE ECHAM5 TIME MANAGER
    CALL ec_manager_init(INT(timestep), nstep)
#endif

#ifndef ICON
    CALL messy_timer_init_manager(timestep, nstep)
#endif

  END SUBROUTINE initialize_restart_attributes
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE bi_decompose(status, flag, reprid, gptr, ptr, p_io_c)

    USE messy_main_channel_repr, ONLY: t_representation, get_representation
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast, p_pe
    USE messy_main_mpi_bi,       ONLY: p_nprocs, p_allgather, p_scatter
    USE messy_main_blather_bi,   ONLY: error_bi

#if defined(ECHAM5)
    USE messy_main_mpi_bi,       ONLY: gather_sa,  gather_sp   &
                                     , scatter_sa, scatter_sp  &
                                     , dcg
    USE messy_main_mpi_bi,       ONLY: gather_field
#endif

#if defined(ECHAM5) || defined(MBM_CLAMS)
    USE messy_main_transform_bi, ONLY: gather_glix, scatter_glix
#endif

!!#D mpiom +
#if defined(ECHAM5) || defined(MBM_MPIOM)
    USE messy_main_mpi_bi,       ONLY: gather_arr, scatter_arr    ! MPIOM
#endif
!!#D mpiom -

#if defined (COSMO) || defined(BLANK) || defined(VERTICO) || defined(MBM_RAD)
    USE messy_main_mpi_bi,       ONLY: gather_gp, scatter_gp
#endif
#if defined(ECHAM5)
    USE messy_main_mpi_bi,       ONLY: scatter_gp
#endif

#if defined(VERTICO) || defined(MBM_RAD)
    USE messy_main_mpi_bi,       ONLY: dcg
#endif

#if defined(I2CINC)
    USE messy_main_mpi_bi,       ONLY: switch_par_utilities
    USE messy_main_blather_bi,   ONLY: info_bi
#endif

#if defined(CESM1)
    USE messy_main_mpi_bi,       ONLY: gather_gp, scatter_gp, dcg
#endif

#ifdef ICON
    USE messy_main_mpi_bi,       ONLY: p_comm_work, p_sum, p_all_comm, p_io0
    USE messy_main_bmluse_bi,    ONLY: wp, sp &
         , my_process_is_mpi_workroot, my_process_is_mpi_seq &
         , p_patch, p_phys_patch, t_phys_patch &
         , n_phys_dom, nproma &
         , exchange_data, t_comm_gather_pattern, idx_no, blk_no &
         , l_output_phys_patch, t_reorder_info &
         , patch_info, t_ScatterPattern &
         , icell, iedge, ivert
    USE messy_main_transform_bi,  ONLY: gather_trix, scatter_trix
    USE messy_main_constants_mem, ONLY: INT_UNDEF
    USE messy_main_channel_repr,  ONLY: IRANK
    USE messy_main_channel_repr,  ONLY: repr_rank_location
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE

    ! I/O
    INTEGER,                      INTENT(OUT) :: status
    INTEGER,                      INTENT(IN)  :: flag
    INTEGER,                      INTENT(IN)  :: reprid
    REAL(DP), DIMENSION(:,:,:,:), POINTER     :: gptr
    REAL(DP), DIMENSION(:,:,:,:), POINTER     :: ptr
    ! I/O PE of current channel, currently only used in ECHAM's decompositions
    ! other decompositions still gather on p_parallel_io
    INTEGER,                      INTENT(IN)  :: p_io_c

    ! LOCAL
    CHARACTER(LEN=*),       PARAMETER :: substr = 'bi_decompose'
    TYPE(t_representation), POINTER   :: repr
    INTEGER :: stat
#if defined(ECHAM5) || defined(MBM_MPIOM)
    INTEGER :: k ! needed for MPIOM
#endif
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: tmpbuf => NULL()
#ifdef ICON
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: tmp_ptr
    INTEGER :: npoints, nblks, n_own, nlev, n_points, nlevs
    INTEGER :: i, ib, il, jk, jp, n, i_dom
    LOGICAL, ALLOCATABLE :: phys_owner_mask(:) ! owner mask for physical patch
    INTEGER, ALLOCATABLE :: glbidx_own(:), glbidx_glb(:)
    INTEGER, ALLOCATABLE :: own_idx(:), own_blk(:) &
                          , own_dst_idx(:), own_dst_blk(:)
    REAL(wp), ALLOCATABLE :: r_tmp(:,:,:)
    REAL(wp), ALLOCATABLE :: r_ptr(:,:,:)
    TYPE(t_phys_patch), POINTER   :: p_phys
    TYPE(t_comm_gather_pattern), POINTER :: p_pat
    CLASS(t_ScatterPattern), POINTER :: scatter_pattern
    TYPE(t_reorder_info),  POINTER :: p_ri
    INTEGER :: jb, jc, lon, lat, rankidx(IRANK), n1, n2, j1, j2
#endif

    ! INIT
    CALL get_representation(status, reprid, repr)
    CALL channel_halt(substr, status)

    SELECT CASE(flag)
    CASE(-1)
       ! ################ RE-COMPOSE ##################################
       !
       ! INIT
       IF (ASSOCIATED(gptr)) THEN
          DEALLOCATE(gptr)
          NULLIFY(gptr)
       END IF
       !
       IF (p_pe == p_io_c) THEN
          ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          IF (stat /= 0) THEN
             status = 1000 ! memory allocation failed
             CALL channel_halt(substr, status)
          END IF
          gptr = 0._dp
       END IF
       !
       SELECT CASE(repr%dctype)

       CASE(DC_BC)
          ! NOTE: MUST BE SYNCHRONIZED ON ALL PEs
          IF (p_pe == p_io_c) &
               gptr(:,:,:,:) = ptr(:,:,:,:)

       CASE(DC_AG)
          ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          IF (stat /= 0) THEN
             status = 1000 ! memory allocation failed
             CALL channel_halt(substr, status)
          END IF
          ALLOCATE (tmpbuf(repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), 0:p_nprocs-1))
          CALL p_allgather(tmpbuf(:,:,:,:), ptr(:,:,:,1))
          gptr(:,:,:,1:p_nprocs) = tmpbuf(:,:,:,0:p_nprocs-1)
          DEALLOCATE(tmpbuf)

#ifdef COSMO
       CASE(DC_GP)
          ! The allocation of the global pointer is required here also
          ! on the non-IO PEs, since the underlying COSMO routine
          ! gather_field expects this.
          IF (.NOT.p_parallel_io) THEN
             ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
                  , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          CALL gather_gp(gptr, ptr)
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
#ifdef I2CINC
       CASE(DC_I2C)
          CALL switch_par_utilities (1)
          ! The allocation of the global pointer is required here also
          ! on the non-IO PEs, since the underlying COSMO routine
          ! gather_field expects this.
          IF (.NOT.p_parallel_io) THEN
             ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
                  , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          CALL gather_gp(gptr, ptr)
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
          CALL switch_par_utilities (2)
       CASE(DC_I2C_IN)
          CALL switch_par_utilities (3)
          ! The allocation of the global pointer is required here also
          ! on the non-IO PEs, since the underlying COSMO routine
          ! gather_field expects this.
          IF (.NOT.p_parallel_io) THEN
             ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
                  , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          CALL gather_gp(gptr, ptr)
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
          CALL switch_par_utilities (2)
       CASE(DC_MMD_OUT,DC_MMD_OUTACC)
          CALL switch_par_utilities (4)
          ! The allocation of the global pointer is required here also
          ! on the non-IO PEs, since the underlying COSMO routine
          ! gather_field expects this.
          IF (.NOT.p_parallel_io) THEN
             ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
                  , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          CALL gather_gp(gptr, ptr, lsum=(repr%dctype==DC_MMD_OUTACC))
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
          CALL switch_par_utilities (2)
#endif
#endif

#if defined(VERTICO) || defined(MBM_RAD)
       CASE(DC_GP)
          CALL gather_gp(gptr, ptr, dcg)
#endif

#if defined(ECHAM5)
       CASE(DC_GP)
          CALL gather_field(gptr, repr%gdimlen, ptr, p_io_c=p_io_c)
       CASE(DC_SP)
          CALL gather_sp(gptr, ptr, dcg, p_io_c=p_io_c)
       CASE(DC_SA)
          CALL gather_sa(gptr, ptr, dcg, p_io_c=p_io_c)
#endif

#if defined(ECHAM5) || defined(MBM_CLAMS)
       CASE(DC_IX)
          CALL gather_glix(gptr, ptr, repr%dcindex &
               , XNG=repr%gdimlen(repr%dcindex))
#endif

!!#D mpiom +
#if defined(ECHAM5) || defined(MBM_MPIOM)
       CASE(DC_GP_MPIOM)
          ! NOTE: The following allocation / deallocation on non-IO PEs
          !       is required, because a section of the array is
          !       explicitely selected: (:,:,k,1). This is not allowed
          !       for a nullified pointer ...
          IF (.NOT.p_parallel_io) THEN
              ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          DO k=1,SIZE(ptr,3)
            CALL gather_arr(ptr(:,:,k,1), gptr(:,:,k,1), p_io)
          ENDDO
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
#endif
!!#D mpiom -

#if defined(CESM1)
       CASE(DC_GP)
          ! The allocation of the global pointer is required here also
          ! on the non-IO PEs, since the underlying CESM1 routine
          ! expects this.
          IF (.NOT.p_parallel_io) THEN
             ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
                  , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          CALL gather_gp(gptr, ptr, dcg,lg32=(reprid == GP_2D_HORIZONTAL))
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
#endif

#if defined(ICON)
       CASE(DC_GP)
          i_dom = repr%dom_id
          IF (i_dom > 0) THEN
             DO jp=1, n_phys_dom
                IF (p_phys_patch(jp)%logical_id == i_dom) THEN
                   p_phys => p_phys_patch(jp)
                END IF
             END DO
             IF (repr%dim(1)%ptr%id == DIMID_NCELLS(i_dom)) THEN
#ifdef ICON_2_1_01
                p_ri => patch_info(i_dom)%cells
#else
                p_ri => patch_info(i_dom)%ri(icell)
#endif
                p_pat => patch_info(i_dom)%p_pat_c
             ELSE IF (repr%dim(1)%ptr%id == DIMID_NEDGES(i_dom)) THEN
#ifdef ICON_2_1_01
                p_ri => patch_info(i_dom)%edges
#else
                p_ri => patch_info(i_dom)%ri(iedge)
#endif
                p_pat => patch_info(i_dom)%p_pat_e
             ELSE IF (repr%dim(1)%ptr%id == DIMID_NVERTS(i_dom)) THEN
#ifdef ICON_2_1_01
                p_ri => patch_info(i_dom)%verts
#else
                p_ri => patch_info(i_dom)%ri(ivert)
#endif
                p_pat => patch_info(i_dom)%p_pat_v
             ELSE
                status = 1
                WRITE(*,*) "MESSY DEBUG: Unkown grid type"
                WRITE(*,'(A18,I3)') "MESSY DEBUG: id = ", repr%dim(1)%ptr%id
                WRITE(*,'(A)') "MESSY DEBUG: name = "//TRIM(repr%name)
                !CALL warning_bi('unknown grid type', substr)
                RETURN
             END IF

             n_points = p_ri%n_glb

             CALL repr_rank_location(repr%axis, rankidx)

             ! Gather data on root
             IF (p_parallel_io) THEN
                gptr(:,:,:,:) = 0._wp
             END IF

             IF (rankidx(1) == 1 .AND. rankidx(2) == 2) THEN
                ! among other  normal 2D case
                n1 = SIZE(ptr,3)
                n2 = SIZE(ptr,4)
                ALLOCATE(r_tmp(MERGE(n_points, 0, &
                     &      my_process_is_mpi_workroot()), n1, n2))
                r_tmp(:,:,:) = 0._wp
                DO j1 = 1, n1
                   DO j2 = 1, n2
                      CALL exchange_data(ptr(:,:,j1,j2), r_tmp(:,j1, j2) &
                           , gather_pattern=p_pat)
                   END DO
                END DO
             ELSE  IF (rankidx(1) == 1 .AND. rankidx(2) == 3) THEN
                n1 = SIZE(ptr,2)
                n2 = SIZE(ptr,4)
                ALLOCATE(r_tmp(MERGE(n_points, 0, &
                     &      my_process_is_mpi_workroot()), n1, n2))
                r_tmp(:,:,:) = 0._wp
                DO j1 = 1, n1
                   DO j2 = 1, n2
                      CALL exchange_data(ptr(:,j1,:,j2), r_tmp(:,j1, j2) &
                           , gather_pattern=p_pat)
                   END DO
                END DO

             ELSE
                CALL error_bi (&
                     'only cases for X on first place an y on 2nd'//&
                     &' or 3rd place implemented' &
                     , substr)
             END IF

             IF(my_process_is_mpi_workroot()) THEN
                IF (rankidx(4) /= INT_UNDEF) THEN ! N rank exists
                   gptr(:,1,:,:) = r_tmp(:,:,:)
                ELSE
                   gptr(:,:,:,1) = r_tmp(:,:,:)
                END IF
             END IF

             DEALLOCATE(r_tmp)
          ELSE
             status = 1
                WRITE(*,*) "MESSY DEBUG: idom < 0"
                !CALL warning_bi('idom < 0 !',substr)
          END IF

!QQQ: GENERAL: DO NOT USE IMPLICIT ASSUMPTIONS ABOUT REPRESENTATIONS,
!     WHICH ARE AVAILABLE FROM THE REPRESENTATION STRUCTURE. THIS IS ERROR
!     PRONE, BECAUSE AT LEAST TWO CODE-POSITIONS NEED TO BE CHANGED, IF THE
!     REPRESENTAION IS MODIFIED.
!     NEXT: DO NOT USE BM-variables (like nproma) WHICH ARE AVAILABLE FROM
!           THE REPRESENTATION STRUCTURE ... (FOR THE SAME REASON!)
       CASE(DC_AGG_SUM)
!qqq      ! NOTE: The following allocation / deallocation on non-IO PEs
          !       is a workaround to bypass Lahey/Fujitsu runtime checks,
          !       which obviously do not allow to pass a nullified pointer
          !       as parameter to a subroutine.
          IF (.NOT.p_parallel_io) THEN
             ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
                  , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          ALLOCATE(tmp_ptr(repr%ldimlen(1),repr%ldimlen(2) &
               ,repr%ldimlen(3),repr%ldimlen(4)))
          tmp_ptr = p_sum(ptr, p_comm_work)

          IF (repr%ldimlen(3) /= 1) THEN
             DO jk = 1, repr%gdimlen(2)
                DO i=1, repr%gdimlen(1) * repr%gdimlen(3)
                   jb = (i-1)/nproma + 1
                   jc = i - (jb-1)*nproma

                   lat = (i - 1)/repr%gdimlen(1) + 1
                   lon = MOD((i - 1),repr%gdimlen(1)) + 1

                   gptr(lon,jk,lat,1) = tmp_ptr(jc,jk,jb,1)
                END DO
             END DO
          ELSE
             DO i=1, repr%gdimlen(1) * repr%gdimlen(2)
                jb = (i-1)/nproma + 1
                jc = i - (jb-1)*nproma

                lat = (i - 1)/repr%gdimlen(1) + 1
                lon = MOD((i - 1),repr%gdimlen(1)) + 1

                gptr(lon,lat,1,1) = tmp_ptr(jc,jb,1,1)
             END DO
          END IF

          DEALLOCATE(tmp_ptr)

          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF

       CASE(DC_AGG_PDF)
!qqq      ! NOTE: The following allocation / deallocation on non-IO PEs
          !       is a workaround to bypass Lahey/Fujitsu runtime checks,
          !       which obviously do not allow to pass a nullified pointer
          !       as parameter to a subroutine.
          IF (.NOT.p_parallel_io) THEN
             ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
                  , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          ALLOCATE(tmp_ptr(repr%ldimlen(1),repr%ldimlen(2) &
               , repr%ldimlen(3),repr%ldimlen(4)))
          tmp_ptr = p_sum(ptr, p_comm_work)

          DO jk = 1, repr%gdimlen(2)
             DO i=1, repr%gdimlen(1) * repr%gdimlen(3)
                jb = (i-1)/nproma + 1
                jc = i - (jb-1)*nproma

                lat = (i - 1)/repr%gdimlen(1) + 1
                lon = MOD((i - 1),repr%gdimlen(1)) + 1

                gptr(lon,jk,lat,:) = tmp_ptr(jc,jk,jb,:)
             END DO
          END DO

          DEALLOCATE(tmp_ptr)

          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF

       CASE(DC_LATLON_NONE)
!qqq      ! NOTE: The following allocation / deallocation on non-IO PEs
          !       is a workaround to bypass Lahey/Fujitsu runtime checks,
          !       which obviously do not allow to pass a nullified pointer
          !       as parameter to a subroutine.
          IF (.NOT.p_parallel_io) THEN
             ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
                  , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          DO jk = 1, repr%gdimlen(2)
             DO i=1, repr%gdimlen(1) * repr%gdimlen(3)
                jb = (i-1)/nproma + 1
                jc = i - (jb-1)*nproma

                lat = (i - 1)/repr%gdimlen(1) + 1
                lon = MOD((i - 1),repr%gdimlen(1)) + 1

                gptr(lon,jk,lat,:) = ptr(jc,jk,jb,:)
             END DO
          END DO
          IF (.NOT.p_parallel_io) THEN
             DEALLOCATE(gptr)
             NULLIFY(gptr)
          ENDIF

       CASE(DC_IX)
          DO i=1, IRANK
             IF (repr%link(i:i) =='x') EXIT
          END DO
          CALL gather_trix(gptr, ptr, i)

       CASE(DC_TRIX)
          CALL gather_trix(gptr, ptr, repr%dcindex) !LaMETTA
#endif

       CASE DEFAULT
          CALL error_bi( &
               'UNKNOWN DECOMPOSITION TYPE FOR RE-COMPOSITION',substr)

       END SELECT
       !
       ! ##############################################################
    CASE(1)
       ! ################ DE-COMPOSE ##################################
       !
       ! INIT
       IF (ASSOCIATED(ptr)) THEN
          DEALLOCATE(ptr)
          NULLIFY(ptr)
       END IF
       !
       ! subtract bounds
       ALLOCATE( ptr( repr%ldimlen(1)-2*repr%bounds%nbounds(1) &
                    , repr%ldimlen(2)-2*repr%bounds%nbounds(2) &
                    , repr%ldimlen(3)-2*repr%bounds%nbounds(3) &
                    , repr%ldimlen(4)-2*repr%bounds%nbounds(4) ), STAT=stat )
       IF (stat /= 0) THEN
          status = 1000 ! memory allocation failed
          CALL channel_halt(substr, status)
       END IF
       !
       SELECT CASE(repr%dctype)

       CASE(DC_BC)
          ! BROADCAST FROM IO-PE TO ALL OTHERS
          IF (p_parallel_io) ptr(:,:,:,:) = gptr(:,:,:,:)
          CALL p_bcast(ptr, p_io)

       CASE(DC_AG)
          ! this is just a workaround for the Intel compiler ...
          IF (.NOT. p_parallel_io) THEN
             IF (ASSOCIATED(gptr)) DEALLOCATE(gptr)
             NULLIFY(gptr)
          END IF
          CALL p_scatter(gptr, ptr(:,:,:,1))

#ifdef COSMO
       CASE(DC_GP)
          ! The allocation of the global pointer is required here also
          ! on the non-IO PEs, since the underlying COSMO routine
          ! expects this.
          IF (.NOT.p_parallel_io) THEN
              ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          CALL scatter_gp(gptr, ptr)
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
#ifdef I2CINC
       CASE (DC_I2C)
          CALL switch_par_utilities (1)
          ! The allocation of the global pointer is required here also
          ! on the non-IO PEs, since the underlying COSMO routine
          ! expects this.
          IF (.NOT.p_parallel_io) THEN
              ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          CALL scatter_gp(gptr, ptr)
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
          CALL switch_par_utilities (2)
       CASE (DC_I2C_IN)
          CALL info_bi('SORRY, NO INPUT OF INT2COSMO INPUT FIELDS POSSIBLE' &
               , substr)
          IF (ASSOCIATED(gptr)) DEALLOCATE(gptr)
          NULLIFY(gptr)
#endif
#endif

#if defined (VERTICO) || defined(MBM_RAD)
       CASE(DC_GP)
          CALL scatter_gp(gptr, ptr, dcg)
#endif

#if defined(ECHAM5)
       CASE(DC_GP)
          CALL scatter_gp(gptr, ptr, dcg)
       CASE(DC_SP)
          CALL scatter_sp(gptr, ptr, dcg)
       CASE(DC_SA)
          CALL scatter_sa(gptr, ptr, dcg)
#endif

#if defined(ECHAM5) || defined(MBM_CLAMS)
       CASE(DC_IX)
          CALL scatter_glix(gptr, ptr, repr%dcindex, xishpg=repr%gdimlen(:))
#endif

!!#D mpiom +
#if defined(ECHAM5) || defined(MBM_MPIOM)
       CASE(DC_GP_MPIOM)
          ! NOTE: The following allocation / deallocation on non-IO PEs
          !       is required, because a section of the array is
          !       explicitely selected: (:,:,k,1). This is not allowed
          !       for a nullified pointer ...
          IF (.NOT.p_parallel_io) THEN
              ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          DO k=1,SIZE(gptr, 3)
            CALL scatter_arr(gptr(:,:,k,1), ptr(:,:,k,1), p_io)
          ENDDO
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
#endif
!!#D mpiom -

#if defined(CESM1)
       CASE(DC_GP)
#ifdef LF
!qqq      ! NOTE: The following allocation / deallocation on non-IO PEs
          !       is a workaround to bypass Lahey/Fujitsu runtime checks,
          !       which obviously do not allow to pass a nullified pointer
          !       as parameter to a subroutine.
          IF (.NOT.p_parallel_io) THEN
             ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
                  , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
#endif
          CALL scatter_gp(gptr, ptr, dcg,lg32=(reprid == GP_2D_HORIZONTAL))
#ifdef LF
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
#endif
#endif

#if defined(ICON)
       CASE(DC_GP)
          i_dom = repr%dom_id
          IF (i_dom >= 0) THEN
             DO jp=1, n_phys_dom
                IF (p_phys_patch(jp)%logical_id == i_dom) THEN
                   p_phys => p_phys_patch(jp)
                END IF
             END DO
             IF (repr%dim(1)%ptr%id == DIMID_NCELLS(i_dom)) THEN
                scatter_pattern => p_patch(i_dom)%comm_pat_scatter_c
             ELSE IF (repr%dim(1)%ptr%id == DIMID_NEDGES(i_dom)) THEN
                scatter_pattern => p_patch(i_dom)%comm_pat_scatter_e
             ELSE IF (repr%dim(1)%ptr%id == DIMID_NVERTS(i_dom)) THEN
                scatter_pattern => p_patch(i_dom)%comm_pat_scatter_v
             ELSE
                status = 1
                WRITE(*,*) "MESSY DEBUG: Unkown grid type"
                WRITE(*,'(A18,I3)') "MESSY DEBUG: id = ", repr%dim(1)%ptr%id
                WRITE(*,'(A)') "MESSY DEBUG: name = "//TRIM(repr%name)
                !CALL warning_bi('unknown grid type', substr)
                RETURN
             END IF

             !WRITE(*,*) "MESSY DEBUG 001: read number of global points"
             n_points = repr%gdimlen(1)
             !WRITE(*,*) "MESSY DEBUG 002: test representation ranks"
             IF (repr%rank == 2) THEN
                nlevs = 1
             ELSE
                nlevs = repr%ldimlen(2)
             END IF

             !WRITE(*,*) "MESSY DEBUG 003: allocate r_ptr"
             IF (nlevs == 1) THEN
                ALLOCATE(r_ptr(repr%ldimlen(1),1,repr%ldimlen(2)))
                r_ptr(:,1,:) = ptr(:,:,1,1)
             ELSE
                ALLOCATE(r_ptr(repr%ldimlen(1),repr%ldimlen(2),repr%ldimlen(3)))
                r_ptr = ptr(:,:,:,1)
             END IF
             !WRITE(*,*) "MESSY DEBUG 004: allocate gptr"
             ALLOCATE(gptr(MERGE(n_points, 0, &
                  my_process_is_mpi_workroot()), nlevs, 1, 1))
             !WRITE(*,*) "MESSY DEBUG 005: scatter data"
             DO jk = 1, nlevs
                CALL scatter_pattern%distribute(gptr(:,jk,1,1) &
                     , r_ptr(:,jk,:), .FALSE.)
             END DO
          ELSE
             status = 1
                WRITE(*,*) "MESSY DEBUG: idom < 0"
                !CALL warning_bi('idom < 0 !',substr)
          END IF

       CASE(DC_TRIX)
          CALL scatter_trix(gptr, ptr, repr%dcindex, xpe=p_io0) !LAMETTA
#endif

       CASE DEFAULT
          CALL error_bi( &
               'UNKNOWN DECOMPOSITION TYPE FOR DE-COMPOSITION',substr)

       END SELECT
       !
       ! ##############################################################
    END SELECT

    status = 0

  END SUBROUTINE bi_decompose
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE bi_vector(status, flag, reprid, ptr, vptr)

    USE messy_main_mpi_bi,       ONLY: dcl, reorder
    USE messy_main_channel_repr, ONLY: t_representation, get_representation

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE

    ! I/O
    INTEGER,                      INTENT(OUT) :: status
    INTEGER,                      INTENT(IN)  :: flag
    INTEGER,                      INTENT(IN)  :: reprid
    REAL(DP), DIMENSION(:,:,:,:), POINTER     :: ptr
    REAL(DP), DIMENSION(:,:,:,:), POINTER     :: vptr

    ! LOCAL
    CHARACTER(LEN=*),       PARAMETER :: substr = 'bi_vector'
    TYPE(t_representation), POINTER   :: repr
    INTEGER :: stat
    INTEGER :: i1, i2, n1, n2

    ! INIT
    CALL get_representation(status, reprid, repr)
    CALL channel_halt(substr, status)

    SELECT CASE(flag)
    CASE(-1)
       ! ################ DE-VECTOR ##################################
       !
       ! INIT
       IF (ASSOCIATED(ptr)) THEN
          DEALLOCATE(ptr)
          NULLIFY(ptr)
       END IF
       !
       IF (.NOT. repr%pdecomp%lpdecomp) THEN
          ! REPRESENTATION DECOMPOSITION TABLE DOES NOT EXIST
          status = 2025
          CALL channel_halt(substr, status)
       END IF
       !
       ALLOCATE( ptr( repr%pdecomp%shape_mem(1), repr%pdecomp%shape_mem(2) &
            , repr%pdecomp%shape_mem(3), repr%pdecomp%shape_mem(4) ) &
            , STAT=stat )
       IF (stat /= 0) THEN
          status = 1000 ! memory allocation failed
          CALL channel_halt(substr, status)
       END IF
       !
       SELECT CASE(repr%dctype)
       CASE(DC_BC)
          !ptr(:,:,:,:) = vptr(:,:,:,:)
          ptr = vptr
       CASE(DC_GP)
          IF (dcl%lreg) THEN
             !ptr(:,:,:,:) = vptr(:,:,:,:)
             ptr = vptr
          ELSE
             SELECT CASE (repr%rank)
             CASE(1)
                ! UNKNOWN REPRESENTATION DECOMPOSITION TYPE / RANK
                status = 2022
                CALL channel_halt(substr, status)
             CASE(2)
                IF (repr%id == GP_3D_1LEV) THEN
                   CALL reorder(ptr(:,1,:,1), vptr(:,1,:,1))
                ELSE
                   CALL reorder(ptr(:,:,1,1), vptr(:,:,1,1))
                END IF
             CASE(3)
                n1 = SIZE(vptr, 2)
                DO i1=1, n1
                   CALL reorder(ptr(:,i1,:,1), vptr(:,i1,:,1))
                END DO
             CASE(4)
                n1 = SIZE(vptr, 2)
                n2 = SIZE(vptr, 3)
                DO i1=1, n1
                   DO i2=1, n2
                      CALL reorder(ptr(:,i1,i2,:), vptr(:,i1,i2,:))
                   END DO
                END DO
             END SELECT
          END IF
       CASE(DC_SP)
          !ptr(:,:,:,:) = vptr(:,:,:,:)
          ptr = vptr
       CASE(DC_SA)
          !ptr(:,:,:,:) = vptr(:,:,:,:)
          ptr = vptr
       CASE(DC_IX)
          !ptr(:,:,:,:) = vptr(:,:,:,:)
          ptr = vptr
       CASE DEFAULT
          status = 1
       END SELECT
       !
       ! ##############################################################
    CASE(1)
       ! ################ VECTOR ######################################
       !
       ! INIT
       IF (ASSOCIATED(vptr)) THEN
          DEALLOCATE(vptr)
          NULLIFY(vptr)
       END IF
       !
       ALLOCATE( vptr( repr%ldimlen(1), repr%ldimlen(2) &
            , repr%ldimlen(3), repr%ldimlen(4) ), STAT=stat )
       IF (stat /= 0) THEN
          status = 1000 ! memory allocation failed
          CALL channel_halt(substr, status)
       END IF
       !
       SELECT CASE(repr%dctype)
       CASE(DC_BC)
          !vptr(:,:,:,:) = ptr(:,:,:,:)
          vptr = ptr
       CASE(DC_GP)
          IF (dcl%lreg) THEN
             !vptr(:,:,:,:) = ptr(:,:,:,:)
             vptr = ptr
          ELSE
             SELECT CASE (repr%rank)
             CASE(1)
                ! UNKNOWN REPRESENTATION DECOMPOSITION TYPE / RANK
                status = 2022
                CALL channel_halt(substr, status)
                ! NOT EXISTENT
             CASE(2)
                IF (repr%id == GP_3D_1LEV) THEN
                   CALL reorder(vptr(:,1,:,1), ptr(:,1,:,1))
                ELSE
                   CALL reorder(vptr(:,:,1,1), ptr(:,:,1,1))
                END IF
             CASE(3)
                n1 = SIZE(ptr, 2)
                DO i1=1, n1
                   CALL reorder(vptr(:,i1,:,1), ptr(:,i1,:,1))
                END DO
             CASE(4)
                n1 = SIZE(ptr, 2)
                n2 = SIZE(ptr, 3)
                DO i1=1, n1
                   DO i2=1, n2
                      CALL reorder(vptr(:,i1,i2,:), ptr(:,i1,i2,:))
                   END DO
                END DO
             END SELECT
          END IF
       CASE(DC_SP)
          !vptr(:,:,:,:) = ptr(:,:,:,:)
          vptr = ptr
       CASE(DC_SA)
          !vptr(:,:,:,:) = ptr(:,:,:,:)
          vptr = ptr
       CASE(DC_IX)
          !vptr(:,:,:,:) = ptr(:,:,:,:)
          vptr = ptr
       CASE DEFAULT
          status = 1
       END SELECT
       !
       ! ##############################################################
    END SELECT

    status = 0

  END SUBROUTINE bi_vector
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE distribute_channels_io

    USE messy_main_mpi_bi,     ONLY: p_all_comm, p_nprocs, get_node_ids, p_io
    USE messy_main_channel_io, ONLY: channel_dist_io
    USE messy_main_blather_bi, ONLY: error_bi

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'distribute_channels_io'
    INTEGER                     :: status
    INTEGER                     :: i

    CALL get_node_ids(status, p_all_comm, nodeid)
    IF (status /= 0) CALL error_bi(substr, 'get_node_ids reported an error')

    CALL channel_dist_io(status, nodeid, p_io, CHANNEL_PPIO_TYPE(:))
    CALL channel_halt(substr, status)

  END SUBROUTINE distribute_channels_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_init_timer

    USE messy_main_mpi_bi,       ONLY: p_io, p_bcast
    USE messy_main_timer_bi,     ONLY: timer_event_init, p_bcast_event
    USE messy_main_channel,      ONLY: get_channel_name
    USE messy_main_tools,        ONLY: match_wild &
                                     , int2str
    USE messy_main_channel_mem,  ONLY: dom_unbound
#ifdef COSMO
    USE data_io,                 ONLY: ngribout
    USE messy_main_tools,        ONLY: str2num
#endif

    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_init_timer'
    INTEGER          :: status
    INTEGER          :: i, js
    CHARACTER(LEN=STRLEN_CHANNEL) :: cname
    LOGICAL          :: lexplicit
    CHARACTER(LEN=8) :: evaldate
    INTEGER          :: dom_id
    CHARACTER(LEN=2) :: domstr
#ifdef COSMO
    INTEGER          :: icn

    ! dimension indices for COSMO output channels
    ALLOCATE(js_COSMOm(ngribout)); js_COSMOm(:) = -1
    ALLOCATE(js_COSMOp(ngribout)); js_COSMOp(:) = -1
    ALLOCATE(js_COSMOz(ngribout)); js_COSMOz(:) = -1
    ALLOCATE(js_COSMOs(ngribout)); js_COSMOs(:) = -1
    ALLOCATE(js_COSMOc(ngribout)); js_COSMOc(:) = -1
#endif

    ! needs to be done here because of TIMER MANAGER initialization
    ! in case of RESTART
    CALL p_bcast(TIMER_DEFAULT%cname, p_io)
    CALL p_bcast_event(TIMER_DEFAULT%io_event, p_io)
    DO i=1, NMAXCHANNELS
       CALL p_bcast(TIMER_CHANNEL(i)%cname, p_io)
       CALL p_bcast_event(TIMER_CHANNEL(i)%io_event, p_io)
    END DO

    CALL p_bcast(TIMER_TNF_DEFAULT%cname, p_io)
    CALL p_bcast_event(TIMER_TNF_DEFAULT%io_event, p_io)
    DO i=1, NMAXCHANNELS
       CALL p_bcast(TIMER_TNF_CHANNEL(i)%cname, p_io)
       CALL p_bcast_event(TIMER_TNF_CHANNEL(i)%io_event, p_io)
    END DO

    ! SPACE FOR ACTIVE OUTPUT EVENTS (ONE PER CHANNEL)
    ALLOCATE(OUTPUT_EVENT(NCHANNEL))
    ALLOCATE(LOUTPUT_NOW(NCHANNEL))
    LOUTPUT_NOW(:) = .FALSE.

    ALLOCATE(I_PATCH(NCHANNEL))
    I_PATCH(:) = dom_unbound

    ALLOCATE(TNF_EVENT(NCHANNEL))
    ALLOCATE(LTNF_NOW(NCHANNEL))
    LTNF_NOW(:) = .FALSE.

    ! PRESET NON-EXPLICIT TO DEFAULT
    channel_loop: DO js=1, NCHANNEL

       CALL get_channel_name(status, js, cname, dom_id)
       CALL channel_halt(substr, status)
       I_PATCH(js) = dom_id

       ! SET EXPLICIT OR DEFAULT
       lexplicit = .FALSE.
       CALL int2str(domstr, dom_id)

       event_loop: DO i=1, NMAXCHANNELS
          IF (TRIM(ADJUSTL(TIMER_CHANNEL(i)%cname)) == '') CYCLE
          IF ( (match_wild(TRIM(ADJUSTL(TIMER_CHANNEL(i)%cname)),&
               TRIM(cname)) ) .OR. &
               (TRIM(ADJUSTL(TIMER_CHANNEL(i)%cname)) == &
                              TRIM(cname)//"_D"//TRIM(domstr)) &
               ) THEN
             lexplicit  = .TRUE.
             EXIT
          END IF
       END DO event_loop

       ! INITIALIZE EVENT
       ! Note: to trigger the output to exactly the time choosen in the namelist
       !       the evaluation date depends on the time integration scheme
#if defined(ICON) || defined(COSMO)
          evaldate =  'present'
#else
          evaldate =  'next'
#endif

       IF (lexplicit) THEN
          CALL timer_event_init(OUTPUT_EVENT(js), TIMER_CHANNEL(i)%io_event &
               , TRIM(cname), TRIM(evaldate) &
#ifdef ICON
               , clock_id=dom_id           &
#endif
               )
       ELSE
          CALL timer_event_init(OUTPUT_EVENT(js), TIMER_DEFAULT%io_event &
               , TRIM(cname),  TRIM(evaldate) &
#ifdef ICON
               , clock_id = dom_id          &
#endif
               )
       END IF

#if defined(ECHAM5)
       ! SPECIAL FOR ECHAM5
       ! NOTE: CHANNEL OUTPUT OF tdiag, tdiag_gp, nudg, nudg_gp
       !       WILL BE SYNCHRONIZED TO THE SHORTEST OUTPUT INTERVAL
       !       (IF PRESENT !)
       SELECT CASE(TRIM(cname))
       CASE('tdiag')
          js_tdiag = js
       CASE('tdiag_gp')
          js_tdiag_gp = js
       CASE('nudg')
          js_nudg = js
       CASE('nudg_gp')
          js_nudg_gp = js
       CASE DEFAULT
          ! DO NOTHING
       END SELECT
#endif

#ifdef COSMO
       ! SPECIAL FOR COSMO
       ! calculation of OUTPUT fields in COSMO will be triggered,
       ! if ANY of the output channels is .TRUE.
       IF (cname(1:5) == 'COSMO') THEN
          CALL str2num(cname(7:9),iCn)
          SELECT CASE(cname(6:6))
          CASE('m')
             js_COSMOm(iCn) = js
          CASE('p')
             js_COSMOp(iCn) = js
          CASE('z')
             js_COSMOz(iCn) = js
          CASE('s')
             js_COSMOs(iCn) = js
          CASE('c')
             js_COSMOs(iCn) = js
          CASE DEFAULT
          ! DO NOTHING
          END SELECT
       ENDIF
#endif

       ! SET EXPLICIT OR DEFAULT
       lexplicit = .FALSE.
       event_loop_tnf: DO i=1, NMAXCHANNELS
          IF (TRIM(ADJUSTL(TIMER_TNF_CHANNEL(i)%cname)) == '') CYCLE
          IF (match_wild(TRIM(ADJUSTL(TIMER_TNF_CHANNEL(i)%cname)),&
               TRIM(cname))) THEN
             lexplicit  = .TRUE.
             EXIT
          END IF
       END DO event_loop_tnf

       ! INITIALIZE EVENT
       IF (lexplicit) THEN
          CALL timer_event_init(TNF_EVENT(js), TIMER_TNF_CHANNEL(i)%io_event &
               , 'TNF'//TRIM(cname), TRIM(evaldate) &
#ifdef ICON
               , clock_id=dom_id           &
#endif
               )
       ELSE
          CALL timer_event_init(TNF_EVENT(js), TIMER_TNF_DEFAULT%io_event   &
               , 'TNF'//TRIM(cname), TRIM(evaldate) &
#ifdef ICON
               , clock_id=dom_id           &
#endif
               )
       END IF

    END DO channel_loop

  END SUBROUTINE main_channel_init_timer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_update_timer

    USE messy_main_timer_bi,    ONLY: event_state
    USE messy_main_channel,     ONLY: trigger_channel_output &
                                    , set_channel_output
    USE messy_main_channel_mem, ONLY: dom_current
#if defined (COSMO) || defined (ICON)
    USE messy_main_timer,       ONLY: current_date, time_days
#else
    USE messy_main_timer,       ONLY: next_date, time_days
#endif

#if defined(ECHAM5)
    USE mo_diag_tendency,       ONLY: dio_index
    USE mo_nudging_buffer,      ONLY: nio_index
    USE mo_nudging_constants,   ONLY: lnudg_out
    USE mo_time_control,        ONLY: l_putdata
#endif

#ifdef COSMO
    USE data_io,                ONLY: ngribout
    USE messy_main_tools,       ONLY: int2str
#endif

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_update_timer'
    INTEGER :: status
    INTEGER :: js
#ifdef ECHAM5
    LOGICAL :: lsync
#endif
    TYPE(time_days)  :: date
#ifdef COSMO
    CHARACTER(LEN=3) :: str
    INTEGER          :: icc
#endif

    ! Note: to trigger the output to exactly the time choosen in the namelist
    !       the evaluation date depends on the time integration scheme
#if defined(ICON) || defined(COSMO)
    date = current_date
#else
    date = next_date
#endif

    DO js = 1, NCHANNEL
       IF (I_PATCH(js) /= dom_current) CYCLE
       LOUTPUT_NOW(js) = event_state(OUTPUT_EVENT(js), date)
    END DO

    DO js = 1, NCHANNEL
       IF (I_PATCH(js) /= dom_current) CYCLE
       LTNF_NOW(js) = event_state(TNF_EVENT(js), date)
    END DO

#ifdef COSMO
    ! SPECIAL FOR COSMO
    ! set L_FORCE_calcout to TRUE, if any COSMO output channel is true
    L_FORCE_calcout = L_CALCOUT_EACH_STEP
#endif

#if defined(ECHAM5)
    ! SPECIAL FOR ECHAM5
    ! NOTE: CHANNEL OUTPUT OF tdiag, tdiag_gp, nudg, nudg_gp
    !       WILL BE SYNCHRONIZED TO THE SHORTEST OUTPUT INTERVAL (IF PRESENT !)
    lsync = .FALSE.
    IF (js_tdiag > 0)    lsync = lsync .OR. LOUTPUT_NOW(js_tdiag)
    IF (js_tdiag_gp > 0) lsync = lsync .OR. LOUTPUT_NOW(js_tdiag_gp)
    IF (js_nudg > 0)     lsync = lsync .OR. LOUTPUT_NOW(js_nudg)
    IF (js_nudg_gp > 0)  lsync = lsync .OR. LOUTPUT_NOW(js_nudg_gp)
    !
    IF (js_tdiag + js_tdiag_gp > 0) l_putdata(dio_index) = lsync
    IF (js_nudg  + js_nudg_gp  > 0) l_putdata(nio_index) = lsync
    !
    IF (js_tdiag > 0)    LOUTPUT_NOW(js_tdiag)    = lsync
    IF (js_tdiag_gp > 0) LOUTPUT_NOW(js_tdiag_gp) = lsync
    IF (js_nudg > 0)     LOUTPUT_NOW(js_nudg)     = lsync
    IF (js_nudg_gp > 0)  LOUTPUT_NOW(js_nudg_gp)  = lsync

    ! SPECIAL FOR ECHAM5
    ! NUDGING OUTPUT ONLY, IF NUDGING IS CURRENTLY ACTIVE
    ! (after nudging stop date has been reached, the nudging memory
    !  is deallocated)
    IF ((js_nudg > 0) .AND. .NOT. lnudg_out) THEN
       CALL set_channel_output(status, 'nudg', .FALSE., .FALSE.)
       CALL channel_halt(substr, status)
    ENDIF
    IF ((js_nudg_gp > 0) .AND. .NOT. lnudg_out) THEN
       CALL set_channel_output(status, 'nudg_gp', .FALSE., .FALSE.)
       CALL channel_halt(substr, status)
    END IF
    ! PREVENT ALSO FROM WRITING RESTART FILES ...
    IF ((js_nudg > 0) .AND. .NOT. lnudg_out) THEN
       CALL set_channel_output(status, 'nudg', .FALSE., .TRUE.)
       CALL channel_halt(substr, status)
    ENDIF
    IF ((js_nudg_gp > 0) .AND. .NOT. lnudg_out) THEN
       CALL set_channel_output(status, 'nudg_gp', .FALSE., .TRUE.)
       CALL channel_halt(substr, status)
    END IF
#endif

#ifdef COSMO
    ! suppress output of channels with constants at the beginning
    DO iCc = 1, ngribout
       IF (js_COSMOc(iCc) > 0) THEN
          CALL int2str(str,iCc)
          CALL set_channel_output(status,'COSMOc'//str, .FALSE.)
          CALL channel_halt(substr,status)
       END IF
    END DO
#endif

    CALL trigger_channel_output(status, LOUTPUT_NOW) ! part
    CALL channel_halt(substr, status)

  END SUBROUTINE main_channel_update_timer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_read_nml_cpl(status, iou)

    ! MODULE ROUTINE (SMIL)
    !
    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2004

    USE messy_main_tools,   ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_main_channel, ONLY: modstr

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CPL/ L_BM_ORIG_OUTPUT     &
         , L_CALCOUT_EACH_STEP          &
         , TIMER_DEFAULT, TIMER_CHANNEL &
         , TIMER_TNF_DEFAULT, TIMER_TNF_CHANNEL &
         , CHANNEL_PPIO_TYPE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_read_nml_cpl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! DO NOT ENABLE ADDITIONAL ECHAM5 STANDARD STREAM OUTPUT PER DEFAULT
    ! (only needed for GRIB-template generation)
    L_BM_ORIG_OUTPUT = .FALSE.
    CHANNEL_PPIO_TYPE(:) = 0

    ! -> OTHER DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_channel_read_nml_cpl
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------

  SUBROUTINE main_channel_make_cf_conform

#ifndef COSMO
    USE messy_main_data_bi,        ONLY: basemodstr => modstr
#endif
    USE messy_main_channel_repr,   ONLY: IRANK
    USE messy_main_channel,        ONLY: t_channel_list, t_channel &
                                       , t_channel_object_list     &
                                       , t_channel_object          &
                                       , new_attribute, get_attribute &
                                       , SND_MAXLEN, GCHANNELLIST &
                                       , new_channel_object                &
                                       , new_channel_object_reference      &
                                       , set_channel_object_inst           &
                                       , NCHANNEL
    USE messy_main_constants_mem,  ONLY: STRLEN_ULONG
    USE messy_main_timer,          ONLY: YEAR_START, MONTH_START, DAY_START &
                                       , HOUR_START, MINUTE_START, SECOND_START
#ifdef ICON
    USE messy_main_bmluse_bi,      ONLY: p_patch, uuid_unparse
#endif
#ifdef COSMO
    USE messy_main_grid_def_mem_bi, ONLY: pollon, pollat, polgam
#endif

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr ='main_channel_make_cf_conform'
    TYPE(t_channel_list),        POINTER, SAVE :: ls
    TYPE(t_channel),             POINTER, SAVE :: channel
    TYPE(t_channel_object_list), POINTER, SAVE :: le
    TYPE(t_channel_object),      POINTER, SAVE :: object
    ! output statistics ? => add timebnds variable
    LOGICAL                                    :: loutstat = .FALSE.
    ! output data on regular lon / lat grid DIMID_LON / DIMID_LAT
    ! add lon / lat variables to channel
    CHARACTER(LEN=STRLEN_ULONG)                :: coordstring = ' '
    CHARACTER(LEN=STRLEN_ULONG)                :: levstring   = ' '
    INTEGER                                    :: ix
    LOGICAL                                    :: llon, llat
    INTEGER                                    :: status
    CHARACTER(LEN=STRLEN_ULONG)                :: ydate = ' '
    INTEGER                                    :: js
    LOGICAL                                    :: lcell
#if defined(COSMO) || defined (ICON)
    LOGICAL                                    :: llonlat   = .FALSE.
#endif
#ifdef COSMO
    LOGICAL                                    :: lslonlat  = .FALSE.
    LOGICAL                                    :: llonslat  = .FALSE.
    LOGICAL                                    :: lslat, lslon
    LOGICAL                                    :: lsoil
#endif
#ifdef ICON
    CHARACTER(LEN=STRLEN_ULONG)                :: griduuid= ' '
    LOGICAL                                    :: ledge, lgrid
    LOGICAL                                    :: lvert
#endif

    ! **************************
    ! ENDLESS LOOP over all CHANNELs
    !
    ! allocate pointers for time_bnds object
    ALLOCATE(TIME_BNDS(NCHANNEL))
    DO js = 1, NCHANNEL
       NULLIFY(TIME_BNDS(js)%ptr)
    END DO

    ls => GCHANNELLIST
    !
    DO
       NULLIFY(le)
       NULLIFY(channel)
       NULLIFY(object)
       loutstat  = .FALSE.
#if defined(COSMO) || defined (ICON)
       llonlat   = .FALSE.
#endif
#ifdef COSMO
       lslonlat  = .FALSE. ! staggered U-wind grid
       llonslat  = .FALSE. ! staggered V-wind_grid
#endif
#ifdef ICON
       lgrid     = .FALSE.
#endif
       !
       IF (.NOT. ASSOCIATED(ls)) THEN
          RETURN
       ELSE
          channel => ls%this
          le => channel%list
       END IF
       ! ********************************
       ! LOOK FOR NEXT OBJECT WITH OUTPUT
       DO
          IF (.NOT. ASSOCIATED(le)) THEN  ! NO MORE OBJECT
             EXIT
          ELSE
             object => le%this
             ! check for lon / lat variable
             llon  = .FALSE.
             llat  = .FALSE.
             lcell = .FALSE.
#ifdef COSMO
             lsoil = .FALSE.
             lslon = .FALSE. ! staggered lon
             lslat = .FALSE. ! staggered lat
#endif
#ifdef ICON
             ledge = .FALSE.
             lvert = .FALSE.
#endif
             coordstring = ' '
             levstring = ' '
             DO ix = 1, IRANK
                IF (ASSOCIATED(object%repr%dim(ix)%ptr)) THEN
#ifndef ICON
                   IF (object%repr%dim(ix)%ptr%id == DIMID_LON) THEN
                      llon = .TRUE.
                   ELSE IF (object%repr%dim(ix)%ptr%id == DIMID_LAT) THEN
                      llat = .TRUE.
#ifdef COSMO
                   ELSE IF (object%repr%dim(ix)%ptr%id == DIMID_SRLON) THEN
                      lslon = .TRUE.
                   ELSE IF (object%repr%dim(ix)%ptr%id == DIMID_SRLAT) THEN
                      lslat = .TRUE.
                   ELSE IF (object%repr%dim(ix)%ptr%id == DIMID_SOIL1) THEN
                      lsoil = .TRUE.
#endif
                   END IF
#else
! ICON
                   IF (channel%dom_id > 0) THEN
                      IF (object%repr%dim(ix)%ptr%id == &
                           DIMID_NCELLS(channel%dom_id)) THEN
                         lcell = .TRUE.
                      ELSE IF (object%repr%dim(ix)%ptr%id == &
                           DIMID_NVERTS(channel%dom_id)) THEN
                         lvert = .TRUE.
                      ELSE IF (object%repr%dim(ix)%ptr%id == &
                           DIMID_NEDGES(channel%dom_id)) THEN
                         ledge = .TRUE.
                      END IF
                   END IF
#endif
                END IF
             END DO
             ! COMPOSE COORDINATE STRING
             ! add lon-lat variables to channel
             IF (llon .AND. llat .OR. lcell ) THEN
                coordstring = TRIM(ADJUSTL(coordstring))//' lon lat'
#if defined(COSMO) || defined (ICON)
                llonlat = .TRUE.
#endif
#ifdef COSMO
             ELSE IF (lslon .AND. llat) THEN
                coordstring =  &
                  TRIM(ADJUSTL(coordstring))//' slonu slatu'
                lslonlat = .TRUE.
             ELSE IF (llon .AND. lslat) THEN
                coordstring = &
                  TRIM(ADJUSTL(coordstring))//' slonv slatv'
                llonslat = .TRUE.
#endif
             END IF
             CALL get_attribute(status, channel%name, object%name &
                  , 'levelindex', c=levstring, dom_id= channel%dom_id)

             IF (status /= 805 ) THEN ! 805 = ATTRIBUTE DOES NOT EXIST
                CALL channel_halt(substr, status)
                coordstring=&
                     TRIM(ADJUSTL(coordstring))//' '//TRIM(ADJUSTL(levstring))
                CALL  new_channel_object_reference(status       &
#ifndef COSMO
                     , basemodstr, TRIM(ADJUSTL(levstring))     &
#else
                     , 'COSMO_ORI', TRIM(ADJUSTL(levstring))    &
#endif
                     , channel%name, TRIM(ADJUSTL(levstring))   &
                     , dom_id1 =  channel%dom_id                &
                     , dom_id2 =  channel%dom_id                &
                     )
                CALL set_channel_object_inst(status           &
                     , channel%name, TRIM(ADJUSTL(levstring)) &
                     , dom_id=channel%dom_id)
                CALL channel_halt(substr//'2',status)
             END IF

             ! ADD attribute
             IF (TRIM(ADJUSTL(coordstring)) /= '') THEN
                CALL new_attribute(status, channel%name, object%name &
                     , 'coordinates', c=TRIM(ADJUSTL(coordstring))   &
                     , loverwrite=.TRUE., dom_id= channel%dom_id)
                IF (status /= 804) & ! attribute exists already
                     CALL channel_halt(substr, status)
#ifdef COSMO
                CALL new_attribute(status, channel%name, object%name &
                     , 'grid_mapping', c='rotated_pole', loverwrite=.TRUE.)
                CALL channel_halt(substr, status)
#endif
             END IF
#ifdef ICON
             IF (lcell) THEN
                CALL new_attribute(status, channel%name, object%name &
                    , 'coordinates (grid_file)'                      &
                    , c='lon_cell_centre lat_cell_centre'            &
                    , dom_id= channel%dom_id)
                CALL channel_halt(substr, status)
                lgrid = .TRUE.
             ELSE IF (ledge) THEN
                CALL new_attribute(status, channel%name, object%name &
                    , 'coordinates (grid_file)'                      &
                    , c='lon_edge_centre lat_edge_centre'            &
                    , dom_id= channel%dom_id)
                CALL channel_halt(substr, status)
                lgrid = .TRUE.
             ELSE IF (lvert) THEN
                CALL new_attribute(status, channel%name, object%name &
                    , 'coordinates (grid_file)'                      &
                    , c='longitude_vertices latitude_vertices'       &
                    , dom_id= channel%dom_id)
                CALL channel_halt(substr, status)
                lgrid = .TRUE.
             END IF
#endif
             ! ******************************
             ! LOOK FOR OTHER OUTPUT than instantaneous
             IF (ANY(object%io%lout(2:SND_MAXLEN))) loutstat = .TRUE.
          END IF
          le => le%next
       END DO
       ! END OBJECT LOOP

       ! ADD required attributes and variables
       ! LAT / LON
       ! FOR EMAC and CESM1 lon / lat are already dimension variables
       ! add dimension variables for ICON
#if defined(COSMO) || defined(ICON)
       IF (llonlat) THEN
          CALL new_channel_object_reference(status &
               , 'grid_def', 'philon_2d'           &
               , channel%name, 'lon'               &
               , dom_id1= channel%dom_id           &
               , dom_id2= channel%dom_id)
          IF (status /= 3102 ) THEN  ! CHANNEL OBJECT EXISTS ALREADY
             CALL channel_halt(substr,status)

             CALL new_attribute(status, channel%name, 'lon' &
                  , 'standard_name', c='longitude', loverwrite=.TRUE.&
                  , dom_id= channel%dom_id)
             CALL channel_halt(substr,status)
             CALL new_attribute(status, channel%name, 'lon' &
                  , 'long_name', c='longitude', loverwrite=.TRUE.&
                  , dom_id= channel%dom_id)
             CALL channel_halt(substr,status)
             CALL new_attribute(status, channel%name, 'lon' &
                  , 'units', c='degrees_east', loverwrite=.TRUE. &
                  , dom_id= channel%dom_id)
             CALL channel_halt(substr,status)
          END IF

          CALL new_channel_object_reference(status &
               , 'grid_def', 'philat_2d'           &
               , channel%name, 'lat'               &
               , dom_id1= channel%dom_id           &
               , dom_id2= channel%dom_id)
          IF (status /= 3102 ) THEN ! CHANNEL OBJECT EXISTS ALREADY
             CALL channel_halt(substr,status)

             CALL new_attribute(status, channel%name, 'lat' &
                  , 'standard_name', c='latitude', loverwrite=.TRUE.&
                  , dom_id= channel%dom_id)
             CALL channel_halt(substr,status)
             CALL new_attribute(status, channel%name, 'lat' &
                  , 'long_name', c='latitude', loverwrite=.TRUE. &
                  , dom_id= channel%dom_id)
             CALL channel_halt(substr,status)
             CALL new_attribute(status, channel%name, 'lat' &
                  , 'units', c='degrees_north', loverwrite=.TRUE. &
                  , dom_id= channel%dom_id)
             CALL channel_halt(substr,status)
          END IF
       END IF
#endif

#ifdef COSMO
       IF (lslonlat) THEN
          CALL new_channel_object_reference(status &
               , 'grid_def', 'philon_U_2d'         &
               , channel%name, 'slonu')
          IF (status /= 3102 ) THEN  ! CHANNEL OBJECT EXISTS ALREADY
             CALL channel_halt(substr,status)

             CALL new_attribute(status, channel%name, 'slonu' &
                  , 'standard_name', c='longitude', loverwrite=.TRUE.)
             CALL channel_halt(substr,status)
             CALL new_attribute(status, channel%name, 'slonu'   &
                  , 'long_name', c='staggered U-wind longitude' &
                  , loverwrite=.TRUE.)
             CALL channel_halt(substr,status)
             CALL new_attribute(status, channel%name, 'slonu' &
                  , 'units', c='degrees_east', loverwrite=.TRUE.)
             CALL channel_halt(substr,status)
          END IF

          CALL new_channel_object_reference(status &
               , 'grid_def', 'philat_U_2d'           &
               , channel%name, 'slatu')
          IF (status /= 3102 ) THEN ! CHANNEL OBJECT EXISTS ALREADY
             CALL channel_halt(substr,status)

             CALL new_attribute(status, channel%name, 'slatu' &
                  , 'standard_name', c='latitude', loverwrite=.TRUE.)
             CALL channel_halt(substr,status)
             CALL new_attribute(status, channel%name, 'slatu'  &
                  , 'long_name', c='staggered U-wind latitude' &
                  , loverwrite=.TRUE.)
             CALL channel_halt(substr,status)
             CALL new_attribute(status, channel%name, 'slatu' &
                  , 'units', c='degrees_north', loverwrite=.TRUE.)
             CALL channel_halt(substr,status)
          END IF
       END IF
       ! ----------------------------------------------------------
       IF (llonslat) THEN
          CALL new_channel_object_reference(status &
               , 'grid_def', 'philon_V_2d'         &
               , channel%name, 'slonv')
          IF (status /= 3102 ) THEN  ! CHANNEL OBJECT EXISTS ALREADY
             CALL channel_halt(substr,status)

             CALL new_attribute(status, channel%name, 'slonv' &
                  , 'standard_name', c='longitude', loverwrite=.TRUE.)
             CALL channel_halt(substr,status)
             CALL new_attribute(status, channel%name, 'slonv'   &
                  , 'long_name', c='staggered V-wind longitude' &
                  , loverwrite=.TRUE.)
             CALL channel_halt(substr,status)
             CALL new_attribute(status, channel%name, 'slonv' &
                  , 'units', c='degrees_east', loverwrite=.TRUE.)
             CALL channel_halt(substr,status)
          END IF

          CALL new_channel_object_reference(status &
               , 'grid_def', 'philat_V_2d'           &
               , channel%name, 'slatv')
          IF (status /= 3102 ) THEN ! CHANNEL OBJECT EXISTS ALREADY
             CALL channel_halt(substr,status)

             CALL new_attribute(status, channel%name, 'slatv' &
                  , 'standard_name', c='latitude', loverwrite=.TRUE.)
             CALL channel_halt(substr,status)
             CALL new_attribute(status, channel%name, 'slatv'  &
                  , 'long_name', c='staggered V-wind latitude' &
                  , loverwrite=.TRUE.)
             CALL channel_halt(substr,status)
             CALL new_attribute(status, channel%name, 'slatv' &
                  , 'units', c='degrees_north', loverwrite=.TRUE.)
             CALL channel_halt(substr,status)
          END IF
       END IF

       IF (lsoil) THEN
          ! SOIL_BNDS
          CALL new_channel_object_reference(status &
               , 'COSMO_ORI', 'soil1_bnds', channel%name, 'soil1_bnds')
          CALL channel_halt(substr,status)
       END IF

       IF (llonlat .OR.lslonlat .OR.llonslat) THEN
          CALL new_channel_object(status,channel%name &
               ,'rotated_pole',reprid=SCALAR, lstatic=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, channel%name,  'rotated_pole'   &
               , 'grid_mapping_name', c='rotated_latitude_longitude' &
               , loverwrite=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, channel%name,  'rotated_pole'   &
               , 'grid_north_pole_latitude', r=pollat                &
               , loverwrite=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, channel%name,  'rotated_pole'   &
               , 'grid_north_pole_longitude', r=pollon               &
               , loverwrite=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, channel%name,  'rotated_pole'   &
               , 'north_pole_grid_longitude', r=polgam               &
               , loverwrite=.TRUE.)
          CALL channel_halt(substr, status)
       END IF

#endif

#ifdef ICON
       ! HERE we take advantage of the knowledge, that each channel is
       ! dedicated to one domain only (I.e. icell == iedge == ivert ...)
       ! if > 0 ...
       IF (lgrid) THEN
          CALL new_attribute(status, channel%name, 'grid_file'  &
               , c=TRIM(p_patch(channel%dom_id)%grid_filename)  &
               , dom_id=channel%dom_id)
          CALL channel_halt(substr, status)
          CALL uuid_unparse(p_patch(channel%dom_id)%grid_uuid,griduuid)
          CALL new_attribute(status, channel%name, 'uuidOfHGrid'  &
               , c=TRIM(griduuid) , dom_id=channel%dom_id)
          CALL channel_halt(substr, status)
       END IF
#endif

       ! TIME_BNDS
       CALL new_channel_object(status, channel%name, "time_bnds"  &
            , reprid=REPR_TIMEBNDS, p1=TIME_BNDS(channel%id)%ptr &
            , lrestreq=loutstat,ldp=.TRUE., dom_id=channel%dom_id)
       CALL channel_halt(substr,status)
       CALL new_attribute(status, channel%name, "time_bnds" &
            , "long_name", c= "time bounds", dom_id=channel%dom_id)
       CALL channel_halt(substr, status)
       WRITE(ydate, '(a11,i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1)')&
            'days since ',YEAR_START,'-',MONTH_START,'-', DAY_START, 'T' &
         , HOUR_START,':',MINUTE_START,':',SECOND_START,'Z'
       CALL new_attribute(status, channel%name, "time_bnds" &
            , "units",c= TRIM(ydate), dom_id=channel%dom_id)
       CALL channel_halt(substr, status)
       CALL set_channel_object_inst(status         &
            , channel%name, 'time_bnds' , dom_id=channel%dom_id)
       CALL channel_halt(substr//'11',status)
       !

#ifdef COSMO
       IF (loutstat) THEN
          IF (channel%name(1:5) == 'COSMO') THEN
             SELECT CASE(channel%name(6:6))
             CASE ('m', 'p', 'z', 's', 'c')
                write (*,*) 'CALCULATION OF SRC_OUTPUT required every timestep!'
                L_CALCOUT_EACH_STEP=.TRUE.
             CASE DEFAULT
                ! COSMO or COSMO_ORI => nothing to do
             END SELECT
          END IF
       END IF
#endif
       ! --------------------
       ! CHECK next channel
       ls => ls%next
    END DO

  END SUBROUTINE main_channel_make_cf_conform
  ! -------------------------------------------------------------------

#ifdef ICON
#ifdef ICON_2_1_01
  !---------------------------------------------------------------------------
  !> Sets the reorder_info for cells/edges/verts
  !  ATTENTION: This routine must only be called on work and test PE
  !  (i.e. not on IO PEs) ...
  !  ... the arguments don't make sense on the IO PEs anyways
  !
  SUBROUTINE set_reorder_info(phys_patch_id, n_points_g, n_points &
       , owner_mask, phys_id  &
       , glb_index, p_ri)

    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_bmluse_bi, ONLY: l_output_phys_patch, t_reorder_info &
         , patch_info &
         , GRID_INFO_NONE, GRID_INFO_FILE, GRID_INFO_BCAST &
         , finish &
         , idx_no, blk_no

  USE mo_mpi,           ONLY: p_bcast, get_my_mpi_work_id, p_max,             &
    &                         get_my_mpi_work_communicator,                   &
    &                         p_comm_work, p_comm_work_2_io,                  &
    &                         p_comm_io, p_comm_work_io,                      &
    &                         p_int, p_real_dp, p_real_sp,                    &
    &                         my_process_is_stdio, my_process_is_mpi_test,    &
    &                         my_process_is_mpi_workroot,                     &
    &                         my_process_is_mpi_seq, my_process_is_io,        &
    &                         my_process_is_mpi_ioroot,                       &
    &                         process_mpi_stdio_id, process_work_io0,         &
    &                         process_mpi_io_size, num_work_procs, p_n_work,  &
    &                         p_pe_work, p_io_pe0, p_pe

    INTEGER, INTENT(IN) :: phys_patch_id   ! Physical patch ID
    INTEGER, INTENT(IN) :: n_points_g      ! Global number of cells/edges/verts in logical patch
    INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
    LOGICAL, INTENT(IN) :: owner_mask(:,:) ! owner_mask for logical patch
    INTEGER, INTENT(IN) :: phys_id(:,:)    ! phys_id for logical patch
    INTEGER, INTENT(IN) :: glb_index(:)    ! glb_index for logical patch
    TYPE(t_reorder_info), INTENT(INOUT) :: p_ri ! Result: reorder info
    ! local variables
    INTEGER :: i, n, il, ib, mpierr
    LOGICAL, ALLOCATABLE :: phys_owner_mask(:) ! owner mask for physical patch
    INTEGER, ALLOCATABLE :: glbidx_own(:), glbidx_glb(:), reorder_index_log_dom(:)

    CHARACTER(LEN=*), PARAMETER :: routine = modstr//"::set_reorder_info"

    ! Just for safety
    IF(my_process_is_io()) CALL finish(routine, 'Must not be called on IO PEs')

    p_ri%n_log = n_points_g ! Total points in logical domain

    ! Set the physical patch owner mask
    ALLOCATE(phys_owner_mask(n_points))
    DO i = 1, n_points
      il = idx_no(i)
      ib = blk_no(i)
      phys_owner_mask(i) = owner_mask(il,ib)
      IF(l_output_phys_patch) &
        phys_owner_mask(i) = phys_owner_mask(i) .AND. &
        (phys_id(il,ib) == phys_patch_id)
    ENDDO

    ! Get number of owned cells/edges/verts (without halos, physical patch only)
    p_ri%n_own = COUNT(phys_owner_mask(:))

    ! Set index arrays to own cells/edges/verts
    ALLOCATE(p_ri%own_idx(p_ri%n_own))
    ALLOCATE(p_ri%own_blk(p_ri%n_own))
    ALLOCATE(glbidx_own(p_ri%n_own)) ! Global index of my own points

    n = 0
    DO i = 1, n_points
      IF(phys_owner_mask(i)) THEN
        n = n+1
        p_ri%own_idx(n) = idx_no(i)
        p_ri%own_blk(n) = blk_no(i)
        glbidx_own(n)   = glb_index(i)
      ENDIF
    ENDDO

    ! Gather the number of own points for every PE into p_ri%pe_own
    ALLOCATE(p_ri%pe_own(0:p_n_work-1))
    ALLOCATE(p_ri%pe_off(0:p_n_work-1))
#ifndef NOMPI
    CALL MPI_Allgather(p_ri%n_own,  1, p_int, &
                       p_ri%pe_own, 1, p_int, &
                       p_comm_work, mpierr)
#else
    p_ri%pe_own(0) = p_ri%n_own
#endif

    ! Get offset within result array
    p_ri%pe_off(0) = 0
    DO i = 1, p_n_work-1
      p_ri%pe_off(i) = p_ri%pe_off(i-1) + p_ri%pe_own(i-1)
    ENDDO

    ! Get global number of points for current (physical!) patch
    p_ri%n_glb = SUM(p_ri%pe_own(:))

    ! Get the global index numbers of the data when it is gathered on PE 0
    ! exactly in the same order as it is retrieved later during I/O
    ALLOCATE(glbidx_glb(p_ri%n_glb))
#ifndef NOMPI
    CALL MPI_Allgatherv(glbidx_own, p_ri%n_own, p_int, &
                        glbidx_glb, p_ri%pe_own, p_ri%pe_off, p_int, &
                        p_comm_work, mpierr)
#else
    glbidx_glb(:) = glbidx_own(:)
#endif

    ! Get reorder_index
    ALLOCATE(p_ri%reorder_index(p_ri%n_glb))
    IF (patch_info(phys_patch_id)%grid_info_mode == GRID_INFO_FILE) THEN
      ALLOCATE(p_ri%grid_info%log_dom_index(p_ri%n_glb))
    END IF
    ! spans the complete logical domain
    ALLOCATE(reorder_index_log_dom(n_points_g))
    reorder_index_log_dom(:) = 0

    DO i = 1, p_ri%n_glb
       ! reorder_index_log_dom stores where a global point in logical
       ! domain comes from.
       ! It is nonzero only at the physical patch locations
       reorder_index_log_dom(glbidx_glb(i)) = i
    ENDDO

    ! Gather the reorder index for the physical domain
    n = 0
    DO i = 1, n_points_g
      IF(reorder_index_log_dom(i)>0) THEN
        n = n+1
        p_ri%reorder_index(reorder_index_log_dom(i)) = n
      ENDIF
    ENDDO

    IF (patch_info(phys_patch_id)%grid_info_mode == GRID_INFO_FILE) THEN
      n = 0
      DO i = 1, n_points_g
        IF(reorder_index_log_dom(i)>0) THEN
          n = n+1
          p_ri%grid_info%log_dom_index(n) = i
        ENDIF
      ENDDO
    END IF

    ! Safety check
    IF(n/=p_ri%n_glb) CALL finish(routine,'Reordering failed')

    ! set trivial destination indices:
    IF(my_process_is_mpi_seq()) THEN
      ALLOCATE(p_ri%own_dst_idx(p_ri%n_own), &
        &      p_ri%own_dst_blk(p_ri%n_own))
      DO i=1,p_ri%n_own
        p_ri%own_dst_idx(i) = idx_no(i)
        p_ri%own_dst_blk(i) = blk_no(i)
      END DO ! i
    END IF

    DEALLOCATE(phys_owner_mask)
    DEALLOCATE(glbidx_own)
    DEALLOCATE(glbidx_glb)
    DEALLOCATE(reorder_index_log_dom)

  END SUBROUTINE set_reorder_info
  !^°-- ICON_2_1_01
#else
  !v-- ICON_2_4_00
  !---------------------------------------------------------------------------
  !> Sets the reorder_info for cells/edges/verts
  !  ATTENTION: This routine must only be called on work and test PE
  ! (i.e. not on IO PEs)
  !             The arguments don't make sense on the IO PEs anyways
  !
  SUBROUTINE set_reorder_info(phys_patch_id, n_points_g, n_points &
       , owner_mask, phys_id, glb_index, p_ri, grid_info)

    USE messy_main_blather_bi,     ONLY: error_bi
    USE messy_main_constants_mem,  ONLY: i8
    USE messy_main_mpi_bi,         ONLY: p_comm_work
    USE messy_main_bmluse_bi,      ONLY: patch_info, l_output_phys_patch &
                                       , GRID_INFO_NONE, GRID_INFO_FILE  &
                                       , GRID_INFO_BCAST &
                                       , idx_no, blk_no
    USE mo_reorder_info,           ONLY: mask2reorder_info
    USE mo_name_list_output_types, ONLY: t_reorder_info, t_grid_info
    USE mo_exception,              ONLY: message_text

    INTEGER, INTENT(IN) :: phys_patch_id   ! Physical patch ID
    INTEGER, INTENT(IN) :: n_points_g      ! Global number of cells/edges/verts in logical patch
    INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
    LOGICAL, INTENT(IN) :: owner_mask(:,:) ! owner_mask for logical patch
    INTEGER, INTENT(IN) :: phys_id(:,:)    ! phys_id for logical patch
    INTEGER, INTENT(IN) :: glb_index(:)    ! glb_index for logical patch
    TYPE(t_reorder_info), INTENT(INOUT) :: p_ri ! Result: reorder info
    TYPE(t_grid_info), INTENT(INOUT) :: grid_info
    ! local variables
    INTEGER :: i, n, il, ib
    LOGICAL, ALLOCATABLE :: phys_owner_mask(:) ! owner mask for physical patch
    INTEGER(i8), ALLOCATABLE :: occupation_mask(:)
    CHARACTER(LEN=*), PARAMETER :: routine = modstr//"::set_reorder_info"

    ! Set the physical patch owner mask
    ALLOCATE(phys_owner_mask(n_points))
    n = 0
    DO i = 1, n_points
      il = idx_no(i)
      ib = blk_no(i)
      phys_owner_mask(i) =       owner_mask(il,ib) &
        &                  .AND. (     .NOT. l_output_phys_patch &
        &                         .OR. phys_id(il,ib) == phys_patch_id)
    ENDDO

    CALL mask2reorder_info(p_ri, phys_owner_mask, n_points_g, glb_index, &
         p_comm_work, occupation_mask)
    DEALLOCATE(phys_owner_mask)

    grid_info%n_log = n_points_g ! Total points in logical domain

    IF (patch_info(phys_patch_id)%grid_info_mode == GRID_INFO_FILE) THEN
      CALL bitmask2start_count_blks(occupation_mask, INT(n_points_g, i8), &
           grid_info%log_dom_starts, grid_info%log_dom_counts)
      ! Safety check
      n = SUM(grid_info%log_dom_counts)
      IF (n /= p_ri%n_glb) THEN
        WRITE(message_text, '(2(a,i0))')  'Reordering failed, n=', n, &
             ', p_ri%n_glb = ', p_ri%n_glb
        CALL error_bi(message_text, routine)
      END IF
    END IF
  END SUBROUTINE set_reorder_info

  SUBROUTINE bitmask2start_count_blks(mask, nb, starts, counts)

    USE messy_main_constants_mem,  ONLY: i8

    INTEGER(i8), INTENT(in) :: mask(0:), nb
    INTEGER, ALLOCATABLE, INTENT(out) :: starts(:), counts(:)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: mask
#endif
    INTEGER(i8) :: i
    INTEGER(i8), PARAMETER :: nbits_i8 = BIT_SIZE(i)
    INTEGER :: num_cblk, n
    LOGICAL :: prev_is_set, current_is_set
    INTEGER(i8) :: pos, apos, bpos, bmask

    num_cblk = 0
    prev_is_set = .FALSE.
    bmask = 1_i8
    DO i = 0_i8, nb-1_i8
      ! extract logical from single bit
      apos = i/nbits_i8
      current_is_set = IAND(mask(apos), bmask) /= 0_i8
      num_cblk = num_cblk + MERGE(1, 0, (.NOT.prev_is_set) .AND. current_is_set)
      bmask = ISHFTC(bmask, 1_i8)
      prev_is_set = current_is_set
    ENDDO
    ALLOCATE(starts(num_cblk), counts(num_cblk))
    bmask = 1_i8
    DO i = 0_i8, nb-1_i8
      apos = i/nbits_i8
      current_is_set = IAND(mask(apos), bmask) /= 0_i8
      bmask = ISHFTC(bmask, 1_i8)
      IF (current_is_set) EXIT
    END DO
    n = 0
    DO WHILE (i < nb)
      ! at this point bit i is set because of the previous loop
      n = n + 1
      ! i is zero-based but starts are Fortran indices
      starts(n) = INT(i+1_i8)
      DO WHILE (current_is_set)
        i = i + 1_i8
        IF (i >= nb) EXIT
        apos = i/nbits_i8
        current_is_set = IAND(mask(apos), bmask) /= 0_i8
        bmask = ISHFTC(bmask, 1_i8)
      END DO
      counts(n) = INT(i+1_i8) - starts(n)
      DO WHILE (.NOT. current_is_set)
        i = i + 1_i8
        IF (i >= nb) EXIT
        apos = i/nbits_i8
        current_is_set = IAND(mask(apos), bmask) /= 0_i8
        bmask = ISHFTC(bmask, 1_i8)
      END DO
    END DO
  END SUBROUTINE bitmask2start_count_blks
#endif
#endif

! **********************************************************************
END MODULE messy_main_channel_bi
! **********************************************************************
