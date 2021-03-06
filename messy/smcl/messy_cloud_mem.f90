MODULE MESSY_CLOUD_MEM

  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE

  PRIVATE :: dp

  INTEGER, PUBLIC :: ncdnc = 0
  INTEGER, PUBLIC :: nicnc = 0

  INTEGER, PUBLIC :: idt_cdnc

  INTEGER, PUBLIC :: idt_icnc

  INTEGER, PUBLIC :: idt_cover1
  INTEGER, PUBLIC :: idt_cover2
  INTEGER, PUBLIC :: idt_cc_cov_1dt
  INTEGER, PUBLIC :: idt_cc_cov_2dt
  INTEGER, PUBLIC :: idt_cc_cov_3dt
  INTEGER, PUBLIC :: idt_cc_cov_4dt
  INTEGER, PUBLIC :: idt_cc_cov_5dt
  INTEGER, PUBLIC :: idt_cc_cov
  INTEGER, PUBLIC :: idt_cc_vol_1dt
  INTEGER, PUBLIC :: idt_cc_vol_2dt
  INTEGER, PUBLIC :: idt_cc_vol_3dt
  INTEGER, PUBLIC :: idt_cc_vol_4dt
  INTEGER, PUBLIC :: idt_cc_vol_5dt
  INTEGER, PUBLIC :: idt_cc_vol
  INTEGER, PUBLIC :: idt_cc_len_1dt
  INTEGER, PUBLIC :: idt_cc_len_2dt
  INTEGER, PUBLIC :: idt_cc_len_3dt
  INTEGER, PUBLIC :: idt_cc_len_4dt
  INTEGER, PUBLIC :: idt_cc_len_5dt
  INTEGER, PUBLIC :: idt_cc_len
  INTEGER, PUBLIC :: idt_cc_iwc
  INTEGER, PUBLIC :: idt_cc_n_i
!
! Define parameters as actual variables, that they can be set in 
! cloud_init_coupling depending on the coupled aerosol submodel
! Distinguish between mixed and soluble indexes
  INTEGER, PUBLIC :: inucs, iaits, iaccs, icoas, &
                            iaitm, iaccm, icoam, &
                            iaiti, iacci, icoai

! Variables for the coupling with the aerosol submodels, moved from 
! messy_cloud_lohmann10 to be accessible from all cloud modules
  INTEGER, PUBLIC :: nmod    = 0    ! number of aerosol modes
  INTEGER, PUBLIC :: nsol    = 0    ! number of soluble aerosol modes
  INTEGER, PUBLIC :: nfrzmod = 0    ! number of aerosol modes for freezing

! Number fraction of dust in coarse mixed mode above which dust
! dominance is assumed (for IN calculation in K14)
  REAL(dp), PUBLIC :: du_nfrac_thr = 0.7_dp

  TYPE aer_species
    INTEGER                           :: idt
    REAL(dp)                          :: density
    REAL(dp)                          :: molarmass
  END TYPE aer_species
  PUBLIC :: aer_species

  TYPE FREEZ_SPEC_ARRAY
    REAL(dp), DIMENSION(:,:), POINTER       :: ptr => NULL()
    INTEGER                                 :: idt
    REAL(dp)                                :: density
    REAL(dp)                                :: molarmass
    CHARACTER(LEN=32)                       :: name
  END TYPE FREEZ_SPEC_ARRAY
  PUBLIC :: freez_spec_array

  TYPE aer_mode_cloud
    INTEGER                           :: nr
    LOGICAL                           :: sol
    LOGICAL                           :: nucl
    LOGICAL                           :: lskip
    REAL(dp), POINTER, DIMENSION(:,:) :: aernum => NULL()
    REAL(dp), POINTER, DIMENSION(:,:) :: wetrad => NULL()
    REAL(dp), POINTER, DIMENSION(:,:) :: dryrad => NULL()
    REAL(dp)                          :: sigma
    REAL(dp)                          :: crdiv
    INTEGER                           :: no_aerspec
    INTEGER                           :: no_freezspec
    TYPE(aer_species),      DIMENSION(:), POINTER   :: aerspec => NULL()
    TYPE(FREEZ_SPEC_ARRAY), DIMENSION(:), POINTER   :: freezspec => NULL()
  END TYPE aer_mode_cloud
  PUBLIC :: aer_mode_cloud

  TYPE(aer_mode_cloud), DIMENSION(:), POINTER, PUBLIC  :: mode => NULL()

 
  INTEGER, PARAMETER, PUBLIC                                :: freez_spec = 4
  CHARACTER(Len=*),DIMENSION(freez_spec), PARAMETER, PUBLIC :: COMP_str = &
       (/'DU   ','BC   ','BCtag', 'OC   '/)


  REAL(dp),PARAMETER :: alv = 2.5008e6_dp ! latent heat for vaporisation in J/kg
  REAL(dp),PARAMETER :: als = 2.8345e6_dp ! latent heat for sublimation in J/kg
                                          ! tropopause calc.

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: swat              => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: qnuc              => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: qmel              => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: qfre              => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: qeva              => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: qacc              => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: qaut              => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: reffl             => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: reffl_acc         => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cloud_tm1         => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cdnc              => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cdnc_acc          => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:)   :: cdnc_burden       => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:)   :: cdnc_burden_acc   => NULL()
  
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: sice              => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: reffi             => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: reffi_acc         => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: icnc              => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: icnc_acc          => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:)   :: icnc_burden       => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:)   :: icnc_burden_acc   => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: lwc_acc           => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: iwc_acc           => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:)   :: cdnc_ct           => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:)   :: reffl_ct          => NULL()


  ! HAM conv_stream

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cdncact_cv        => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: nbcsol_strat      => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: ndusol_strat      => NULL()

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: nbcsol_cirrus     => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: nbctagsol_cirrus  => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: ndusol_cirrus     => NULL()
  ! variables needed for BN09 scheme
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: nocsolks          => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: nocsolas          => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: nocsolcs          => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: nocinsol          => NULL()

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: nbcinsol          => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: nbctaginsol       => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: nduinsolai        => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: nduinsolci        => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: naerinsol         => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: naersol           => NULL()

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: twc_conv          => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: conv_time         => NULL()

  ! Different components of the vertical velocity
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: vervel_ls         => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: vervel_tke        => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: vervel_gw         => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: vervel_p18        => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:)   :: random_2d         => NULL()

  ! IC in the different freezing modes
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: Nice_preex     => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: Nice_DUdep     => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: Nice_DUimm     => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: Nice_BC        => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: Nice_BCtag     => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: Nice_homog     => NULL()

  ! For comparison with CDNC observational data
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cdnc_insitu    => NULL()

  ! Cirrus-specific diagnostics for aircraft data by M. Kraemer (FZJ)
  INTEGER, PARAMETER :: N_CIRRUS_TEMP     = 62

  INTEGER, PARAMETER :: N_CIRRUS_BIN_IWC  = 30
  INTEGER, PARAMETER :: N_CIRRUS_BIN_Nice = 35
  INTEGER, PARAMETER :: N_CIRRUS_BIN_Rice = 25
  INTEGER, PARAMETER :: N_CIRRUS_BIN_RHi  = 24
  INTEGER, PARAMETER :: N_CIRRUS_BIN_VEL  = 30

  REAL(dp), DIMENSION(N_CIRRUS_TEMP)         :: CIRRUS_TEMP

  REAL(dp), DIMENSION(N_CIRRUS_BIN_IWC)      :: CIRRUS_BIN_IWC
  REAL(dp), DIMENSION(N_CIRRUS_BIN_Nice)     :: CIRRUS_BIN_Nice
  REAL(dp), DIMENSION(N_CIRRUS_BIN_Rice)     :: CIRRUS_BIN_Rice
  REAL(dp), DIMENSION(N_CIRRUS_BIN_RHi)      :: CIRRUS_BIN_RHi
  REAL(dp), DIMENSION(N_CIRRUS_BIN_VEL)      :: CIRRUS_BIN_VEL

  REAL(dp), DIMENSION(N_CIRRUS_BIN_IWC + 1)  :: CIRRUS_IBIN_IWC
  REAL(dp), DIMENSION(N_CIRRUS_BIN_Nice + 1) :: CIRRUS_IBIN_Nice
  REAL(dp), DIMENSION(N_CIRRUS_BIN_Rice + 1) :: CIRRUS_IBIN_Rice
  REAL(dp), DIMENSION(N_CIRRUS_BIN_RHi + 1)  :: CIRRUS_IBIN_RHi
  REAL(dp), DIMENSION(N_CIRRUS_BIN_VEL + 1)  :: CIRRUS_IBIN_VEL

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:,:) :: CIRRUS_IWC       => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:,:) :: CIRRUS_Nice      => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:,:) :: CIRRUS_Nice_ML   => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:,:) :: CIRRUS_Rice      => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:,:) :: CIRRUS_RHi_cloud => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:,:) :: CIRRUS_RHi_clear => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:,:) :: CIRRUS_vervel    => NULL()

! for CCMod
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: ccfh2o         => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: ccfkme         => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: ccfh2o_invent  => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: ccfkme_invent  => NULL()

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: B_co           => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: B_cc           => NULL()

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:)   :: gboxarea_2d    => NULL()

  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: ccvol          => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cccov          => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:)   :: cccov_tot      => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:)   :: cccov_tot2     => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: ccicnc         => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cciwc          => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: ccreffi        => NULL()
  REAL(dp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cctau          => NULL()

END MODULE MESSY_CLOUD_MEM
