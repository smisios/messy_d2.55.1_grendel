! *************************************************************************
MODULE messy_main_switch
! *************************************************************************

  ! DEFINITION AND NAMELIST INPUT OF GLOBAL SUBMODEL SWITCHES
  ! NOTE: TO ADD A NEW SWITCH LOOK FOR
  !       '### ADD HERE'

  IMPLICIT NONE
  PUBLIC

  CHARACTER(LEN=*), PARAMETER :: modstr = 'switch'
  CHARACTER(LEN=*), PARAMETER :: modver = '1.0'

  ! GLOBAL SWITCHES
  !
  LOGICAL :: L_TIME_INFO  = .true.
  !
  LOGICAL :: USE_A2O      = .false.
  LOGICAL :: USE_ACCF     = .false.
  LOGICAL :: USE_AEROPT   = .false.
  LOGICAL :: USE_AIRSEA   = .false.
  LOGICAL :: USE_AIRTRAC  = .false.
  LOGICAL :: USE_AIRTRAF  = .false.
  LOGICAL :: USE_ATTILA   = .false.
  LOGICAL :: USE_AVEOUT   = .false.
  LOGICAL :: USE_BIOBURN  = .false.
  LOGICAL :: USE_BUFLY    = .false.
  LOGICAL :: USE_CAT      = .false.
  LOGICAL :: USE_CH4      = .false.
  LOGICAL :: USE_CHEMGLUE = .false.
  LOGICAL :: USE_CLOUD    = .false.
  LOGICAL :: USE_CLOUDJ   = .false.
  LOGICAL :: USE_CLOUDOPT = .false.
  LOGICAL :: USE_CONTRAIL = .false.
  LOGICAL :: USE_CONVECT  = .false.
  LOGICAL :: USE_CRM      = .false.
  LOGICAL :: USE_COSMOTOY = .false.
  LOGICAL :: USE_CVTRANS  = .false.
  LOGICAL :: USE_D14CO    = .false.
  LOGICAL :: USE_DDEP     = .false.
  LOGICAL :: USE_DISSOC   = .false.
  LOGICAL :: USE_DIUMOD   = .false.
  LOGICAL :: USE_DOMINT   = .false.
  LOGICAL :: USE_DRADON   = .false.
  LOGICAL :: USE_E4CHEM   = .false.
  LOGICAL :: USE_E5VDIFF  = .false.
  LOGICAL :: USE_EDITH    = .false.
  LOGICAL :: USE_EC2COSMO = .false.
  LOGICAL :: USE_EVER     = .false.
  LOGICAL :: USE_GMXE     = .false.
  LOGICAL :: USE_GWAVE    = .false.
  LOGICAL :: USE_H2O      = .false.
  LOGICAL :: USE_H2OEMIS  = .false.
  LOGICAL :: USE_H2OISO   = .false.
  LOGICAL :: USE_HD       = .false.
  LOGICAL :: USE_HAMOCC   = .false.
  LOGICAL :: USE_IONS     = .false.
  LOGICAL :: USE_ISOPCOR  = .false.
  LOGICAL :: USE_JVAL     = .false.
  LOGICAL :: USE_JVST     = .false.
  LOGICAL :: USE_LGGP     = .false.
  LOGICAL :: USE_LGTMIX   = .false.
  LOGICAL :: USE_LGVFLUX  = .false.
  LOGICAL :: USE_LNOX     = .false.
  LOGICAL :: USE_M7       = .false.
  LOGICAL :: USE_MADE     = .false.
  LOGICAL :: USE_MADE3    = .false.
  LOGICAL :: USE_MECCA    = .false.
  LOGICAL :: USE_MEGAN    = .false.
  LOGICAL :: USE_MESOENERGY = .false.
  LOGICAL :: USE_MLOCEAN  = .false.
  LOGICAL :: USE_MMFORCE  = .false.
  LOGICAL :: USE_MPIOM    = .false.
  LOGICAL :: USE_MSBM     = .false.
  LOGICAL :: USE_MTSKIP   = .false.
  LOGICAL :: USE_NAN      = .false.
  LOGICAL :: USE_O3ORIG   = .false.
  LOGICAL :: USE_OASIS3MCT= .false. ! ub_ak_20171221
  LOGICAL :: USE_OFFEMIS  = .false.
  LOGICAL :: USE_ONEMIS   = .false.
  LOGICAL :: USE_ORACLE   = .false.
  LOGICAL :: USE_ORBIT    = .false.
  LOGICAL :: USE_OROGW    = .false.
  LOGICAL :: USE_OTPHYSC  = .false.
  LOGICAL :: USE_PLUMEGAS = .false.
  LOGICAL :: USE_PTRAC    = .false.
  LOGICAL :: USE_PTRACINI = .false.
  LOGICAL :: USE_QBO      = .false.
  LOGICAL :: USE_RAD      = .false.
  LOGICAL :: USE_RELAX    = .false.
  LOGICAL :: USE_RNDTEST  = .false.
  LOGICAL :: USE_SATSIMS  = .false.
  LOGICAL :: USE_SCALC    = .false.
  LOGICAL :: USE_SCAV     = .false.
  LOGICAL :: USE_SCOUT    = .false.
  LOGICAL :: USE_SEDI     = .false.
  LOGICAL :: USE_S4D      = .false.
  LOGICAL :: USE_SF6      = .false.
  LOGICAL :: USE_SORBIT   = .false.
  LOGICAL :: USE_SPACENOX = .false.
  LOGICAL :: USE_SVOC     = .false.
  LOGICAL :: USE_GEC      = .false.
  LOGICAL :: USE_SPE      = .false.
  LOGICAL :: USE_SURFACE  = .false.
  LOGICAL :: USE_TAGGING  = .false.
  LOGICAL :: USE_TBUDGET  = .false.
  LOGICAL :: USE_TIMEPOS  = .false.
  LOGICAL :: USE_TNUDGE   = .false.
  LOGICAL :: USE_TPULSE   = .false.
  LOGICAL :: USE_TREXP    = .false.
  LOGICAL :: USE_TROPOP   = .false.
  LOGICAL :: USE_TRSYNC   = .false.
  LOGICAL :: USE_UBCNOX   = .false.
  LOGICAL :: USE_VAHR     = .false.
  LOGICAL :: USE_VAXTRA   = .false.
  LOGICAL :: USE_VERTDIFF = .false.
  LOGICAL :: USE_VERTEX   = .false.
  LOGICAL :: USE_VISO     = .false.
  LOGICAL :: USE_VISOP    = .false.
  LOGICAL :: USE_MMD2WAY  = .false.
  LOGICAL :: USE_SUBMOD1  = .false.
  LOGICAL :: USE_SUBMOD2  = .false.
  LOGICAL :: USE_SUBMOD3  = .false.
  LOGICAL :: USE_TESTEVENT= .false.
! ju_ch_20110429+
  LOGICAL :: USE_CLAMS       = .false.
  LOGICAL :: USE_CLAMSTRAJ   = .false.
  LOGICAL :: USE_CLAMSCHEM   = .false.
  LOGICAL :: USE_CLAMSCHEME5 = .false.
  LOGICAL :: USE_CLAMSMIX    = .false.
  LOGICAL :: USE_CLAMSBMIX   = .false.
  LOGICAL :: USE_CLAMSCIRRUS = .false.
  LOGICAL :: USE_CLAMSRDFRC  = .false.
  LOGICAL :: USE_CLAMSSEDI   = .false.
  LOGICAL :: USE_CLAMSTRACER = .false.
  LOGICAL :: USE_CLAMSDEEPCONV=.false.
! ju_ch_20110429-
  LOGICAL :: USE_MXL         = .false. ! mz_rj_20150206

  ! ### ADD HERE

  NAMELIST /CTRL/ L_TIME_INFO, &
     L_TIME_INFO,  USE_A2O,      USE_ACCF,     USE_AEROPT,     &
     USE_AIRSEA,   USE_AIRTRAC,  USE_AIRTRAF,  USE_ATTILA,     &
     USE_AVEOUT,   USE_BIOBURN,  USE_BUFLY,    USE_CAT,        &
     USE_CH4,      USE_CHEMGLUE, USE_CLAMS,    USE_CLAMSBMIX,  &
     USE_CLAMSCHEM,USE_CLAMSCHEME5,USE_CLAMSCIRRUS,USE_CLAMSDEEPCONV,  &
     USE_CLAMSMIX, USE_CLAMSRDFRC,USE_CLAMSSEDI,USE_CLAMSTRACER,  &
     USE_CLAMSTRAJ,USE_CLOUD,    USE_CLOUDJ,   USE_CLOUDOPT,   &
     USE_CONTRAIL, USE_CONVECT,  USE_COSMOTOY, USE_CRM,        &
     USE_CVTRANS,  USE_D14CO,    USE_DDEP,     USE_DISSOC,     &
     USE_DIUMOD,   USE_DOMINT,   USE_DRADON,   USE_E4CHEM,     &
     USE_E5VDIFF,  USE_EC2COSMO, USE_EDITH,    USE_EVER,       &
     USE_GEC,      USE_GMXE,     USE_GWAVE,    USE_H2O,        &
     USE_H2OEMIS,  USE_H2OISO,   USE_HAMOCC,   USE_HD,         &
     USE_IONS,     USE_ISOPCOR,  USE_JVAL,     USE_JVST,       &
     USE_LGGP,     USE_LGTMIX,   USE_LGVFLUX,  USE_LNOX,       &
     USE_M7,       USE_MADE,     USE_MADE3,    USE_MECCA,      &
     USE_MEGAN,    USE_MESOENERGY,USE_MLOCEAN,  USE_MMD2WAY,    &
     USE_MMFORCE,  USE_MPIOM,    USE_MSBM,     USE_MTSKIP,     &
     USE_MXL,      USE_NAN,      USE_O3ORIG,   USE_OASIS3MCT,  &
     USE_OFFEMIS,  USE_ONEMIS,   USE_ORACLE,   USE_ORBIT,      &
     USE_OROGW,    USE_OTPHYSC,  USE_PLUMEGAS, USE_PTRAC,      &
     USE_PTRACINI, USE_QBO,      USE_RAD,      USE_RELAX,      &
     USE_RNDTEST,  USE_S4D,      USE_SATSIMS,  USE_SCALC,      &
     USE_SCAV,     USE_SCOUT,    USE_SEDI,     USE_SF6,        &
     USE_SORBIT,   USE_SPACENOX, USE_SPE,      USE_SUBMOD1,    &
     USE_SUBMOD2,  USE_SUBMOD3,  USE_SURFACE,  USE_SVOC,       &
     USE_TAGGING,  USE_TBUDGET,  USE_TESTEVENT,USE_TIMEPOS,    &
     USE_TNUDGE,   USE_TPULSE,   USE_TREXP,    USE_TROPOP,     &
     USE_TRSYNC,   USE_UBCNOX,   USE_VAHR,     USE_VAXTRA,     &
     USE_VERTDIFF, USE_VERTEX,   USE_VISO,     USE_VISOP
  ! ### ADD HERE

  !PUBLIC :: messy_main_read_nml_ctrl
  !PUBLIC :: switch_init
  ! op_pj_20110621+
  ! temprarily moved to messy_main_switch_bi.f90 to reduce dependencies
  ! (via USE_MLOCEAN) of BML to SMCL
!!$  PRIVATE :: put_submodel_att
  ! op_pj_20110621-

CONTAINS

  ! -------------------------------------------------------------
  SUBROUTINE messy_main_read_nml_ctrl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'messy_main_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1  ! DEFAULT: ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! NO ERROR

  END SUBROUTINE messy_main_read_nml_ctrl
  ! -------------------------------------------------------------

! *************************************************************************
END MODULE messy_main_switch
! *************************************************************************
