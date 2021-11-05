MODULE MESSY_GMXE_MEM


  USE messy_main_constants_mem, ONLY: strlen_medium, dp

  IMPLICIT NONE
  PUBLIC
  SAVE

  INTRINSIC :: NULL
  
  INTEGER   :: nmod
  INTEGER   :: nsoluble, ndiff

  ! TD species structure
  TYPE gas_info
     INTEGER,  DIMENSION(:), POINTER :: aerml_idx    => NULL()
     REAL(dp)                        :: molmass      = 0._dp
     REAL(dp)                        :: density      = 0._dp
     REAL(dp), DIMENSION(:), POINTER :: caccgas      => NULL()
     INTEGER                         :: iso_idx      = 0
     INTEGER                         :: eqs_idx      = 0
     CHARACTER(LEN=STRLEN_MEDIUM)    :: name         = ''
     LOGICAL,  DIMENSION(:), POINTER :: ltreat       => NULL()
     CHARACTER(LEN=STRLEN_MEDIUM)    :: ionname      = ''
  END type gas_info
  TYPE anion_info
     INTEGER,  DIMENSION(:), POINTER :: aerml_idx  => NULL()
     REAL(dp)                        :: molmass    = 0._dp    
     REAL(dp)                        :: density    = 0._dp
     INTEGER                         :: iso_idx    = 0
     INTEGER                         :: eqs_idx    = 0
     INTEGER                         :: gas_idx    = 0
     REAL(dp)                        :: charge     = 0._dp
     CHARACTER(LEN=STRLEN_MEDIUM)    :: name         = ''
     LOGICAL,  DIMENSION(:), POINTER :: ltreat       => NULL()
  END type anion_info
  TYPE cation_info
     INTEGER,  DIMENSION(:), POINTER :: aerml_idx  => NULL()
     REAL(dp)                        :: molmass    = 0._dp
     REAL(dp)                        :: density    = 0._dp
     INTEGER                         :: iso_idx    = 0
     INTEGER                         :: eqs_idx    = 0
     INTEGER                         :: gas_idx    = 0
     REAL(dp)                        :: charge     = 0._dp
     CHARACTER(LEN=STRLEN_MEDIUM)    :: name       = ''
     LOGICAL,  DIMENSION(:), POINTER :: ltreat     => NULL()
  END type cation_info
  TYPE solute_info
     INTEGER,  DIMENSION(:), POINTER :: aerml_idx  => NULL()
     REAL(dp)                        :: molmass    = 0._dp
     REAL(dp)                        :: density    = 0._dp
     INTEGER                         :: iso_idx    = 0
     INTEGER                         :: eqs_idx    = 0
     INTEGER                         :: gas_idx    = 0
     CHARACTER(LEN=STRLEN_MEDIUM)    :: name         = ''
     LOGICAL,  DIMENSION(:), POINTER :: ltreat       => NULL()
  END type solute_info
  TYPE bulk_info
     INTEGER,  DIMENSION(:), POINTER :: aerml_idx  => NULL()
     REAL(dp)                        :: molmass    = 0._dp
     REAL(dp)                        :: density    = 0._dp
     REAL(dp)                        :: kappa      = 0._dp
     CHARACTER(LEN=STRLEN_MEDIUM)    :: name       = ''
     LOGICAL,  DIMENSION(:), POINTER :: ltreat     => NULL()
     LOGICAL                         :: l_oc       =.FALSE.
  END type bulk_info

  TYPE td_info
     TYPE(gas_info),    DIMENSION(:), POINTER :: gas           => NULL()
     TYPE(anion_info),  DIMENSION(:), POINTER :: anion         => NULL()
     TYPE(cation_info), DIMENSION(:), POINTER :: cation        => NULL()
     TYPE(solute_info), DIMENSION(:), POINTER :: solute        => NULL()
     INTEGER                                  :: gas_sulph_idx = 0
     INTEGER                                  :: hp_idx        = 0
  END type td_info
  TYPE(td_info)                    :: td   
  TYPE(bulk_info), DIMENSION(:), POINTER :: bulk    => NULL()

  TYPE species_info
    INTEGER,  DIMENSION(:), POINTER  :: tracidx     => NULL()
    INTEGER,  DIMENSION(:), POINTER  :: aermlidx    => NULL()
    INTEGER,  DIMENSION(:), POINTER  :: aernlidx    => NULL()
    INTEGER,  DIMENSION(:), POINTER  :: kppidx      => NULL()
    INTEGER,  DIMENSION(:), POINTER  :: SOAIDX      => NULL()
    REAL(dp), DIMENSION(:), POINTER  :: zcoup       => NULL()
    CHARACTER(LEN=STRLEN_MEDIUM)     :: name        = ""
    REAL(dp)                         :: molmass     = 0.0_dp 
    REAL(dp)                         :: density     = 0.0_dp 
    REAL(dp)                         :: kappa       = 0.0_dp 
    LOGICAL,  DIMENSION(:), POINTER  :: npassive    => NULL() ! mz_dk_20120119
    INTEGER                          :: charge      = 0
  END TYPE species_info
  TYPE(species_info), DIMENSION(:), POINTER :: species => NULL()

  INTEGER                                   :: spec_number = 0
  INTEGER                                   :: nbulk       = 0
  INTEGER                                   :: nanions     = 0
  INTEGER                                   :: ncations    = 0
  INTEGER                                   :: ngas        = 0
  INTEGER                                   :: nsolutes    = 0

  INTEGER                                   :: spec_idx_oh    = 0
  INTEGER                                   :: spec_idx_ohm   = 0
  INTEGER                                   :: spec_idx_hp    = 0
  INTEGER                                   :: spec_idx_nh3   = 0
  INTEGER                                   :: spec_idx_nh4p  = 0
  INTEGER                                   :: spec_idx_h2so4 = 0
  INTEGER                                   :: spec_idx_hso4m = 0
  INTEGER                                   :: spec_idx_so4mm = 0
  INTEGER                                   :: spec_idx_h2o   = 0
  INTEGER                                   :: spec_idx_oc    = 0
!  INTEGER                                   :: spec_idx_bc    = 0
!  INTEGER                                   :: spec_idx_ss    = 0
!  INTEGER                                   :: spec_idx_du    = 0
  INTEGER                                   :: number_idx     = 0

  REAL(dp), PUBLIC, POINTER :: sigma(:)      => NULL()   ! width of modes

  REAL(dp), PUBLIC, POINTER :: crdiv(:)      => NULL()     ! radius boundary per mode
  REAL(dp), PUBLIC, POINTER :: crdiv_r3(:)   => NULL()  ! radius boundary per mode to the third power
  REAL(dp), PUBLIC, POINTER :: crdiv_mid(:)  => NULL() ! mean radius per mode

END MODULE MESSY_GMXE_MEM
