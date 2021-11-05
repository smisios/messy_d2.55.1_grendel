MODULE MESSY_AEROPT_MEM

  USE MESSY_MAIN_CONSTANTS_MEM,   ONLY: dp, STRLEN_MEDIUM, STRLEN_ULONG
  USE messy_main_tools,           ONLY: PTR_3D_ARRAY
  USE messy_main_channel,         ONLY: STRLEN_CHANNEL, STRLEN_OBJECT

  IMPLICIT NONE

  SAVE

  CHARACTER(LEN=*), PUBLIC, PARAMETER :: MODSTR='aeropt'
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: MODVER='2.0.2'
  TYPE nml_lut_struct
    INTEGER                      :: lut_number             = 0
    CHARACTER(LEN=STRLEN_ULONG)  :: sw_filename = ""
    CHARACTER(LEN=STRLEN_ULONG)  :: lw_filename = ""
  END TYPE nml_lut_struct
  INTEGER, PARAMETER                                :: MAX_AEROPT_CALLS = 25
  TYPE(nml_lut_struct), DIMENSION(MAX_AEROPT_CALLS) :: read_lut_sets
  NAMELIST /CTRL/ read_lut_sets

!--------------------------------------------------------------------------
!---------------- B E G I N   R A D I A T I O N   S T U F F ---------------
!--------------------------------------------------------------------------
! aerosol optical properties mapped to ECHAM5 radiation scheme 
! (adopted from MADE)
!

  INTEGER, PUBLIC :: NSW = 4  ! number of shortwave bands of SW scheme
                                        ! (need to match the defintion in 
                                        !  messy_rad_short)
                                        !  4 bands: as in ECMWF model
                                    ! visible : BAND 1 : 0.25 - 0.69 micrometer
                                    ! near IR : BAND 2 : 0.69 - 1.19 micrometer
                                    !           BAND 3 : 1.19 - 2.38 micrometer
                                    !           BAND 4 : 2.38 - 4.00 micrometer
  INTEGER, PUBLIC :: JPBAND = 16  ! number of longwave bands of 
                                          ! LW scheme (need to match the 
                                          ! defintion in messy_rad_long)
                 !  BAND 1:    10- 250 cm-1 (low - H2O; high - H2O)
                 !  BAND 2:   250- 500 cm-1 (low - H2O; high - H2O)
                 !  BAND 3:   500- 630 cm-1 (low - H2O,CO2; high - H2O,CO2)
                 !  BAND 4:   630- 700 cm-1 (low - H2O,CO2; high - O3,CO2)
                 !  BAND 5:   700- 820 cm-1 (low - H2O,CO2; high - O3,CO2)
                 !  BAND 6:   820- 980 cm-1 (low - H2O; high - nothing)
                 !  BAND 7:   980-1080 cm-1 (low - H2O,O3; high - O3)
                 !  BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)
                 !  BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)
                 !  BAND 10: 1390-1480 cm-1 (low - H2O; high - H2O)
                 !  BAND 11: 1480-1800 cm-1 (low - H2O; high - H2O)
                 !  BAND 12: 1800-2080 cm-1 (low - H2O,CO2; high - nothing)
                 !  BAND 13: 2080-2250 cm-1 (low - H2O,N2O; high - nothing)
                 !  BAND 14: 2250-2380 cm-1 (low - CO2; high - CO2)
                 !  BAND 15: 2380-2600 cm-1 (low - N2O,CO2; high - nothing)
                 !  BAND 16: 2600-3000 cm-1 (low - H2O,CH4; high - nothing)
  

  INTEGER,  PARAMETER,PUBLIC                    :: max_diag_wavelens = 16

  INTEGER :: ncID_sw   ! NetCDF ID of file containing lookuptables (shortwave)
  INTEGER :: ncID_lw   ! NetCDF ID of file containing lookuptables (longwave)

  ! Dimensions of look-up tables for optical properties of aerosol modes

!  INTEGER, DIMENSION(2) :: idim_nr  = 0 ! 1st dim (real refractive index)
!  INTEGER, DIMENSION(2) :: idim_ni  = 0 ! 2nd dim (imaginary refractive index)
!  INTEGER, DIMENSION(2) :: idim_msp = 0 ! 3rd dim (mie size parameter)

  INTEGER, PARAMETER :: max_sw_intervals = 32 ! max. number of sub-intervals
                                              ! for shortwave bands
  INTEGER, PARAMETER :: max_lw_intervals = 32 ! max. number of sub-intervals
                                              ! for longwave bands

  INTEGER, PARAMETER :: max_aerosols = 16     ! max. number of aerosol
                                              ! components

  INTEGER, POINTER :: num_opt_wavelens => NULL()  ! number of wavelengths for optional
                                                  ! diagnostic output


  ! indices for arrays containing refractive indices of different aerosol types

  INTEGER, PARAMETER :: num_rad_types = 8 ! number of aerosol types for
                                          ! aerosol optical properties
  INTEGER, PARAMETER :: ri_bc      = 1 ! BC
  INTEGER, PARAMETER :: ri_oc      = 2 ! OC
  INTEGER, PARAMETER :: ri_du      = 3 ! mineral dust
  INTEGER, PARAMETER :: ri_ss      = 4 ! sea salt
  INTEGER, PARAMETER :: ri_so4     = 5 ! SO4
  INTEGER, PARAMETER :: ri_ammsulf = 6 ! (NH4)2SO4
!!$  INTEGER, PARAMETER :: ri_waso    = 7 ! OPAC "water soluble"
  INTEGER, PARAMETER :: ri_waso    = 6 ! OPAC "water soluble"
!cb  changed from 7 to 5 because refractive index should not be like dust
  INTEGER, PARAMETER :: ri_h2o     = 8 ! H2O

  ! total number of entries in aerosol optical thickness array

  INTEGER, PARAMETER,PUBLIC :: num_diag_elements = 7

  ! indices for array containing aerosol optical thickness

  INTEGER, PARAMETER, PUBLIC :: diag_bc    = 1 ! BC
  INTEGER, PARAMETER, PUBLIC :: diag_oc    = 2 ! OC
  INTEGER, PARAMETER, PUBLIC :: diag_du    = 3 ! DU
  INTEGER, PARAMETER, PUBLIC :: diag_ss    = 4 ! sea salt
  ! OPAC "water soluble" = SO4/NH4/NO3
  INTEGER, PARAMETER, PUBLIC :: diag_waso  = 5 
  INTEGER, PARAMETER, PUBLIC :: diag_h2o   = 6 ! aerosol liquid water
  INTEGER, PARAMETER, PUBLIC :: diag_tot   = 7 ! total
  
  !--------------------------------------------------------------------------
  !------------------ E N D   R A D I A T I O N   S T U F F -----------------
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !--------------------- B E G I N    J V A L   S T U F F -------------------
  !--------------------------------------------------------------------------
  ! JVAL has 8 bands, but band 0 does not use aerosol information and 
  ! can therefore be neglected -> is removed from the wavelength boundary parameters
  INTEGER, PARAMETER, PUBLIC :: n_jv      = 7       ! number of bands in JVAL
  INTEGER, PARAMETER, PUBLIC :: n_jv_calc = 7       ! number of bands for
                                                    !   calculations for JVAL
  ! JVAL lower wavelength boundary for the seven bands (µm)
  REAL(dp), PARAMETER, PUBLIC :: jval_wavelen_l(n_jv) = &
!    (/0.1786_dp,0.202_dp,0.241_dp,0.2899_dp,0.3055_dp,0.3135_dp,0.3375_dp,0.4225_dp/)
    (/0.202_dp,0.241_dp,0.2899_dp,0.3055_dp,0.3135_dp,0.3375_dp,0.4225_dp/)
  ! JVAL upper wavelength boundary for the seven bands (µm)
  REAL(dp), PARAMETER, PUBLIC :: jval_wavelen_u(n_jv) = &
!    (/0.202_dp,0.241_dp,0.2899_dp,0.3055_dp,0.3135_dp,0.3375_dp,0.4225_dp,0.6825_dp/)
    (/0.241_dp,0.2899_dp,0.3055_dp,0.3135_dp,0.3375_dp,0.4225_dp,0.6825_dp/)
  ! one can calculate from the band boundaries a mean wavelength per band
  ! alternative approach: determine specific efective wavelnegth per band (according to C. Bruehl 
  !                       based on solar spectrum and ozone / o2 absorption)
  REAL(dp), PARAMETER, PUBLIC :: jval_wavelen_lc(n_jv) = &
   (/0.2051_dp,0.2879_dp,0.302_dp,0.309_dp,0.320_dp,0.370_dp,0.580_dp/)
  LOGICAL :: L_mean_wavelength = .FALSE.

  REAL(dp), PUBLIC :: jval_wavelen(n_jv_calc)

  !--------------------------------------------------------------------------
  !--------------------- E N D    J V A L   S T U F F -----------------------
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  !----------- B E G I N    A E R O S O L    C O U P L I N G  S T U F F -----
  !--------------------------------------------------------------------------
  
  TYPE t_aero_mode
    ! number of compounds per mode
    INTEGER                         :: count      =  0
    ! number of tracer compounds per mode
    INTEGER                         :: count0     =  0
    ! index number of number tracer associated with this compound
    INTEGER                         :: idx_num    =  0
    ! radius change for mass
    REAL(dp)                        :: ram2cmr    =  0._dp
    ! index number of sodium (Na+) tracer
    INTEGER                         :: Na_idx     =  0
    ! index number of Chloride (Cl-) tracer
    INTEGER                         :: Cl_idx     =  0
    ! index number of Seasalt (SS-) tracer
    INTEGER                         :: SS_idx     =  0

    ! logical switch for water calculation in this mode from volumes
    LOGICAL                         :: l_calc_water = .false.
    ! index number of Water (H2O) compound
    INTEGER                         :: H2O_idx    =  0
    
    ! LOGICAL for the treatment of seasalt
    LOGICAL                         :: l_calc_seasalt = .TRUE.
    ! index number of Seasalt (SS) compound
    INTEGER                         :: seasalt_idx2 =  0
    ! index number of sodium (Na+) compound
    INTEGER                         :: Na_idx2      =  0
    ! index number of Chloride (Cl-) compound
    INTEGER                         :: Cl_idx2      =  0
    
    ! index number of radiation species associated with this compound
    INTEGER,  DIMENSION(:),       POINTER :: rad_spec   => NULL()
    ! index number of tracer associated with this compound
    INTEGER,  DIMENSION(:),       POINTER :: idx_trac   => NULL()
    ! molar mass of this compound
    REAL(dp), DIMENSION(:),       POINTER :: molarmass  => NULL()
    ! density of this compound
    REAL(dp), DIMENSION(:),       POINTER :: density    => NULL()
    ! names of this compound
    CHARACTER(LEN=32), DIMENSION(:), POINTER :: name    => NULL()
    ! mixing ratios of compounds
    REAL(dp), DIMENSION(:,:,:),   POINTER :: paerml     => NULL()

  END TYPE t_aero_mode
  TYPE(t_aero_mode), DIMENSION(:), POINTER :: aerspec   => NULL()

  TYPE t_inp_data_sw
     REAL(dp), DIMENSION(:,:,:), POINTER :: extinct => NULL()
     REAL(dp), DIMENSION(:,:,:), POINTER :: omega   => NULL()
     REAL(dp), DIMENSION(:,:,:), POINTER :: gamma   => NULL()
  END type t_inp_data_sw

  TYPE t_inp_data_lw
     REAL(dp), DIMENSION(:,:,:), POINTER :: extinct => NULL()
  END type t_inp_data_lw

  INTEGER, PARAMETER             :: max_modes = 9

  TYPE t_aero_set
    ! array information of species
    TYPE(t_aero_mode), DIMENSION(:), POINTER :: aerspec_e5   => NULL()
    ! aerosol mode width
    REAL(dp), POINTER, DIMENSION(:)     :: sigma
    ! total aerosol compounds for this aero_set
    INTEGER                             :: naertot
    ! total aerosol modes and number of soluble modes for this aero_set
    INTEGER                             :: nmod = 0
    INTEGER                             :: nsol = 0
    ! name for this aero_set
    CHARACTER(LEN=STRLEN_MEDIUM)        :: aero_set_name
    ! aerosol optical thicknes (AOT), shortwave (sw) / longwave (lw)
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: aot_sw   => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: aot_lw   => NULL()
    ! aerosol single scattering albedo, shortwave (sw)
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: omega_sw => NULL()
    ! aerosol asymmetry factor, shortwave (sw)
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: gamma_sw => NULL()

    REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: aot_sw_5d   => NULL()
    REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: aot_lw_5d   => NULL()
    ! aerosol single scattering albedo, shortwave (sw)
    REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: omega_sw_5d => NULL()
    ! aerosol asymmetry factor, shortwave (sw)
    REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: gamma_sw_5d => NULL()
    ! aerosol optical thickness (optional output, 2D fields)
    TYPE (PTR_3D_ARRAY), DIMENSION (:,:,:), POINTER     :: aot_opt=> NULL()
    ! aerosol extinction coefficient (optional output, 2D fields)
    TYPE (PTR_3D_ARRAY), DIMENSION (:,:,:), POINTER     :: extcoeff_opt=> NULL()

    ! CPL elements for this aero_set
    ! pointers for coupling to aerosol stream
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: wetradius=> NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: dryradius=> NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: aernumber=> NULL()
    LOGICAL         :: lcalc_seasalt   =.FALSE.
    LOGICAL         :: l_extmixt       =.FALSE.
    LOGICAL         :: l_tracer        =.FALSE.

    ! links to jval 
    ! calculating jval values for this set
    LOGICAL         :: lcalc_jval          =.FALSE.
    ! number of jval calculation bands for this set
    INTEGER         :: n_jv_bands = 0            
    ! aerosol scattering extinction
    REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: jv_asca_5d   => NULL()
    REAL(dp), DIMENSION(:,:,:,:),   POINTER :: jv_asca      => NULL()
    ! aerosol absorbption extinction
    REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: jv_aabs_5d   => NULL()
    REAL(dp), DIMENSION(:,:,:,:),   POINTER :: jv_aabs      => NULL()
    ! aerosol asymmetry factor
    REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: jv_ga_5d     => NULL()
    REAL(dp), DIMENSION(:,:,:,:),   POINTER :: jv_ga        => NULL()

    LOGICAL, DIMENSION(:), POINTER       :: spec_exclude => NULL()
    CHARACTER(LEN=STRLEN_ULONG)          :: exclude_string_full=""
    CHARACTER(LEN=STRLEN_ULONG)          :: exclude_string_spec=""
    CHARACTER(LEN=STRLEN_MEDIUM)         :: aermodelname     = ""
    CHARACTER(LEN=STRLEN_MEDIUM)         :: tracerset        = ""

    TYPE(t_inp_data_sw), DIMENSION(:), POINTER :: inp_sw => NULL()
    TYPE(t_inp_data_lw), DIMENSION(:), POINTER :: inp_lw => NULL()
    INTEGER                       :: inp_dim = 4
    INTEGER                       :: inp_nsw = 0
    INTEGER                       :: inp_nlw = 0
    LOGICAL                       :: l_flip  = .FALSE.

    INTEGER                       :: merge_idx1 = 0
    INTEGER                       :: merge_idx2 = 0
    REAL(dp)                      :: weight1 = 1.0_dp
    REAL(dp)                      :: weight2 = 1.0_dp
    REAL(dp)                      :: press_bound1 = 0.0_dp
    REAL(dp)                      :: press_bound2 = 0.0_dp
    LOGICAL                       :: lconv_sw = .FALSE.
    LOGICAL                       :: lconv_lw = .FALSE.

    INTEGER          :: num_opt_wavelens
    INTEGER          :: lut_number
    INTEGER          :: znwave           = 0

    REAL(dp), DIMENSION(max_diag_wavelens) :: diag_wavelen = 0.0_dp
    REAL(dp), DIMENSION(:),   POINTER :: lambda            => NULL()
    REAL(dp), DIMENSION(:),   POINTER :: lambda_squared    => NULL()
    REAL(dp), DIMENSION(:),   POINTER :: weight            => NULL()
    REAL(dp), DIMENSION(:,:), POINTER :: ref_re            => NULL()
    REAL(dp), DIMENSION(:,:), POINTER :: ref_im            => NULL()
    INTEGER,  DIMENSION(:),   POINTER :: int2band          => NULL()
    INTEGER          :: NS 
    INTEGER          :: KS 
    INTEGER          :: AS 
    INTEGER          :: CS 
    INTEGER, DIMENSION(max_modes) :: radmod = 0
    LOGICAL          :: TANRE =.FALSE.

  END TYPE t_aero_set

  INTEGER, POINTER :: NS => NULL()
  INTEGER, POINTER :: KS => NULL()
  INTEGER, POINTER :: AS => NULL()
  INTEGER, POINTER :: CS => NULL()

  INTEGER, DIMENSION(:), POINTER :: radmod => NULL()

  CHARACTER(LEN=STRLEN_ULONG), POINTER :: rad_sw_filename   => NULL()
  CHARACTER(LEN=STRLEN_ULONG), POINTER :: rad_lw_filename   => NULL()

  TYPE lookuptable_set
    INTEGER            :: nmodrad = 3        ! number of aerosol radiation modes 
                            ! (e.g. 3 for aitken, accumulation, coarse)
    INTEGER            :: table_number      = 0
    CHARACTER(LEN=STRLEN_ULONG) :: rad_sw_filename   = ''
    CHARACTER(LEN=STRLEN_ULONG) :: rad_lw_filename   = ''
    ! Dimensions of look-up tables for optical properties of aerosol modes

    INTEGER, DIMENSION(:), POINTER :: idim_nr  =>NULL() ! 1st dim (real refractive index)
    INTEGER, DIMENSION(:), POINTER :: idim_ni  =>NULL() ! 2nd dim (imaginary refractive index)
    INTEGER, DIMENSION(:), POINTER :: idim_msp =>NULL() ! 3rd dim (mie size parameter)


    INTEGER          :: num_sw_intervals = 0         ! number of sub-intervals for shortwave bands
    INTEGER          :: num_lw_intervals = 0         ! number of sub-intervals for longwave bands

    INTEGER          :: num_aerosols     = 0         ! number of aerosol components
    ! --------------------------------------------------------------------------
    ! temporary fields
    INTEGER, POINTER, DIMENSION(:)     :: sw_intervals ! number of sub-intervals per
                                                       ! ECHAM5 sw band
    INTEGER, POINTER, DIMENSION(:)     :: lw_intervals ! number of sub-intervals per
                                                       ! ECHAM5 lw band
    ! --------------------------------------------------------------------------
    ! originally DIMENSION(2,nmodrad)
    ! min/max real refractive index
    REAL(dp), POINTER, DIMENSION(:,:) :: nr_min
    REAL(dp), POINTER, DIMENSION(:,:) :: nr_max 
    ! min/max imaginary refractive index
    REAL(dp), POINTER, DIMENSION(:,:) :: ni_min
    REAL(dp), POINTER, DIMENSION(:,:) :: ni_max
    ! min/max of mie-size-param.                                   
    REAL(dp), POINTER, DIMENSION(:,:) :: msp_min
    REAL(dp), POINTER, DIMENSION(:,:) :: msp_max
    ! increment of nr in look-up table
    REAL(dp), POINTER, DIMENSION(:,:) :: nr_step
    ! increment of ni in look-up table
    REAL(dp), POINTER, DIMENSION(:,:) :: ni_step
    ! increment of msp in look-up table
    REAL(dp), POINTER, DIMENSION(:,:) :: msp_step
    REAL(dp), POINTER, DIMENSION(:,:) :: log_ni_min
    REAL(dp), POINTER, DIMENSION(:,:) :: log_msp_min

    ! Look-up tables for optical properties of aerosols

    ! wavelength of sub-intervals (center) for E5 shortwave+longwave bands [m]
    ! DIMENSION(max_sw_intervals+max_lw_intervals+max_diag_wavelens)
    REAL(dp), DIMENSION(:), POINTER  :: lambda
    REAL(dp), DIMENSION(:), POINTER  :: lambda_squared
    REAL(dp), DIMENSION(:), POINTER  :: weight
    REAL(dp), DIMENSION(:), POINTER  :: lambda_sw ! temporary
    REAL(dp), DIMENSION(:), POINTER  :: lambda_lw ! temporary  
    REAL(dp), DIMENSION(:), POINTER  :: weight_sw ! temporary
    REAL(dp), DIMENSION(:), POINTER  :: weight_lw ! temporary

    ! sum of weights of sub-intervals for mapping to ECHAM5-bands
    REAL(dp), DIMENSION(:), POINTER                :: sumweight_sw
    REAL(dp), DIMENSION(:), POINTER                :: sumweight_lw
    ! refractive index, real part
    REAL(dp), DIMENSION(:,:), POINTER :: ref_re
    REAL(dp), DIMENSION(:,:), POINTER :: ref_re_sw => NULL() ! temporary
    REAL(dp), DIMENSION(:,:), POINTER :: ref_re_lw => NULL() ! temporary
    
    ! refractive index, imaginary part
    REAL(dp), DIMENSION(:,:), POINTER :: ref_im
    REAL(dp), DIMENSION(:,:), POINTER :: ref_im_sw => NULL() ! temporary
    REAL(dp), DIMENSION(:,:), POINTER :: ref_im_lw => NULL() ! temporary

    ! extinction cross section divided by wavelength squared [part-1]
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: lut_sw_sigma => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: lut_lw_sigma => NULL()
    ! single scattering albedo [1]
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: lut_sw_omega => NULL()
    ! assymetry factor [1]
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: lut_sw_gamma => NULL()

   ! indices for mapping sub-intervals to ECHAM5-bands
    INTEGER,  DIMENSION(:),   POINTER :: int2band => NULL()
    REAL(dp), DIMENSION(:,:), POINTER :: sigma_g  => NULL()
  END TYPE lookuptable_set

  TYPE(lookuptable_set), DIMENSION(:), POINTER :: tab_set => NULL()

  ! locally used Pointers that are directed to the entries 
  ! of the corresponding tab_sets

  REAL(dp), DIMENSION(:),   POINTER :: sigma           => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: sigma_g         => NULL()
  LOGICAL,                  POINTER :: l_extmixt       => NULL()
  LOGICAL,                  POINTER :: lcalc_seasalt   => NULL()
  LOGICAL,                  POINTER :: lcalc_jval      => NULL()

  INTEGER, POINTER               :: nmodrad          => NULL()              
  INTEGER, POINTER               :: nmod             => NULL()              
  INTEGER, POINTER               :: nsol             => NULL()              
  INTEGER, POINTER               :: naertot          => NULL()              

  INTEGER, POINTER               :: num_sw_intervals => NULL()             
  INTEGER, POINTER               :: num_lw_intervals => NULL()             
  INTEGER, POINTER               :: num_aerosols     => NULL()              
  INTEGER, POINTER               :: znwave           => NULL()              


  INTEGER, DIMENSION(:), POINTER :: idim_nr         => NULL()
  INTEGER, DIMENSION(:), POINTER :: idim_ni         => NULL()
  INTEGER, DIMENSION(:), POINTER :: idim_msp        => NULL()

  INTEGER, POINTER, DIMENSION(:)    :: sw_intervals => NULL() 
  INTEGER, POINTER, DIMENSION(:)    :: lw_intervals => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:) :: nr_min       => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:) :: nr_max       => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:) :: ni_min       => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:) :: ni_max       => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:) :: msp_min      => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:) :: msp_max      => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:) :: nr_step      => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:) :: ni_step      => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:) :: msp_step     => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:) :: log_ni_min   => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:) :: log_msp_min  => NULL() 

  REAL(dp), DIMENSION(:), POINTER   :: lambda         => NULL() 
  REAL(dp), DIMENSION(:), POINTER   :: lambda_squared => NULL() 
  REAL(dp), DIMENSION(:), POINTER   :: weight         => NULL() 
  REAL(dp), DIMENSION(:), POINTER   :: lambda_sw      => NULL() 
  REAL(dp), DIMENSION(:), POINTER   :: lambda_lw      => NULL() 
  REAL(dp), DIMENSION(:), POINTER   :: weight_sw      => NULL() 
  REAL(dp), DIMENSION(:), POINTER   :: weight_lw      => NULL() 

  REAL(dp), DIMENSION(:), POINTER   :: sumweight_sw => NULL()
  REAL(dp), DIMENSION(:), POINTER   :: sumweight_lw => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: ref_re       => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: ref_re_sw    => NULL() 
  REAL(dp), DIMENSION(:,:), POINTER :: ref_re_lw    => NULL() 
  REAL(dp), DIMENSION(:,:), POINTER :: ref_im       => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: ref_im_sw    => NULL() 
  REAL(dp), DIMENSION(:,:), POINTER :: ref_im_lw    => NULL() 

  REAL(dp), DIMENSION(:,:,:,:), POINTER :: lut_sw_sigma => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: lut_lw_sigma => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: lut_sw_omega => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: lut_sw_gamma => NULL()
  
  INTEGER, DIMENSION(:), POINTER :: int2band => NULL()


  CHARACTER(LEN=STRLEN_MEDIUM), POINTER :: aermodelname     => NULL()
  CHARACTER(LEN=STRLEN_MEDIUM), POINTER :: tracerset        => NULL()
  CHARACTER(LEN=STRLEN_MEDIUM), POINTER :: aero_set_name    => NULL()
  CHARACTER(LEN=STRLEN_ULONG),  POINTER :: exclude_string_full => NULL()
  CHARACTER(LEN=STRLEN_ULONG),  POINTER :: exclude_string_spec => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: aot_sw   => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: aot_lw   => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: omega_sw => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: gamma_sw => NULL()
  REAL(dp), DIMENSION(:),       POINTER :: rad_diag_wavelen   => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: wetradius=> NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: dryradius=> NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: aernumber=> NULL()

  TYPE (PTR_3D_ARRAY), DIMENSION (:,:,:), POINTER :: aot_opt => NULL()
  TYPE (PTR_3D_ARRAY), DIMENSION (:,:,:), POINTER :: extcoeff_opt => NULL()

  LOGICAL,                      POINTER :: l_tracer   => NULL()

  INTEGER,                      POINTER :: lut_number => NULL()

  INTEGER,                      POINTER :: n_jv_bands => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: jv_asca    => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: jv_aabs    => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: jv_ga      => NULL()
END MODULE MESSY_AEROPT_MEM
