!_______________________________________________________________________________
MODULE messy_gmxe
!_______________________________________________________________________________
! DESCRIPTION
! -----------
! GMXe CORE MODULE FOR ECHAM5/MESSY AND MESSY BOX-MODEL
!
! AUTHORS
! ------
! Holger Tost, JGU, Mainz, Germany
! Swen Metzger, Max Planck Institute for Chemistry, Mainz, Germany
! Kirsty Pringle, Max Planck Institute for Chemistry, Mainz, Germany
! questions/suggestions: holger.tost@mpic.de, swen.metzger@mpic.de
! Copyright 2006-2010+. All rights reserved.
!_______________________________________________________________________________
! GMXe = Generalized version of M7 microphysics
!        added thermodynamics by ISORROPIA or EQSAM4CLIM
!              explicit SOA treatment
!              detailed aqueous phase chemistry
!              paramterised organic aerosol ageing
!
!
! M7 was originally implemented in ECHAM by Philip Stier, MPI-M, Hamburg, 2001.
! Original M7 source code (box model) by J. Wilson & E. Vignati, JRC/EI,  2000.
!
! For details see:
! http://www.geosci-model-dev.net/3/391/2010/gmd-3-391-2010.html
!
! For license see:
! http://www.messy-interface.org/
!_______________________________________________________________________________
!_______________________________________________________________________________
! USE messy_main_constants_mem, ONLY: dp
  USE messy_main_blather,         ONLY: start_message, end_message

  USE MESSY_GMXE_MEM
  USE MESSY_GMXE_SOA,             ONLY: umode_soa, lmode_soa, l_gasprec, l_oxi, &
                                        excl_str_soa
  USE messy_main_constants_mem,   ONLY: pi, avo => N_a, R_gas, T0, MW => M_H2O, &
                                        Dw => rho_h2o, p0 => atm2Pa, M_air, g,  &
                                        mwh2o => M_H2O
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  ! SUBROUTINES
  PUBLIC :: gmxe_main
  PUBLIC :: gmxe_read_nml_ctrl
  PUBLIC :: gmxe_initialize_core, gmxe_initialize_species
  PUBLIC :: init_paerosol
  !
  INTRINSIC EXP, LOG, MAX, MIN, SQRT, ALOG10, MINVAL

  PUBLIC :: pi, g, avo, mwh2o, m_air
  !
  ! GLOBAL PARAMETERS ==========================================================
  !
  ! SUBMODEL
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr  = 'gmxe'     ! name of module
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver  = '2.2.0'    ! module version
  !
  ! General definitions ========================================================
  !
  INTEGER, PUBLIC,PARAMETER  :: IO = 6  ! debugging unit
!_______________________________________________________________________________
  !

 !_______________________________________________________________________________
  !
  !--- 1) Pre-set switches for the processes of GMXe:
  !
  ! CTRL-NAMELIST PARAMETERS (to be defined in gmxe.nml)
  LOGICAL, PUBLIC :: &
             ! <--- GMXe general options
             lgmxe      = .FALSE.,& ! Aerosol composition scheme GMXe
             loutput    = .FALSE.,& ! Write GMXe stream to diagnose output
             l_io       = .FALSE.,& ! write(*,*) only if p_parallel_io
             lmass_diag = .FALSE.,& ! Mass balance check in gmxe_interface
             lstrat     = .FALSE.,& ! Exclude stratosphere (requires tropop stream)
             lpsc       = .FALSE.,& ! Exclude psc region (requires psc stream)
             lnucl      = .FALSE.,& ! Calculate nucleation of aerosol particles
             lcoat      = .FALSE.,& ! Calculate coating of primary particles
             lsize      = .FALSE.,& ! Reshape the size distribution due to hygroscopic growth
             lcoag      = .FALSE.,& ! Reshape the size distribution due to particle coagulation
             lcond      = .FALSE.,& ! Calculate condensation rate of aerosol precursor gases (otherwise equilibrium is assumed)
             lkappa     = .TRUE.,& ! Calculate aerosol kappa, water uptake due to kappa and CCN at fixed Supersaturation
             ladyn      = .FALSE.,& ! Calculate aerosol dynamics (equilibration time limited by transport time)
             lah2o      = .FALSE.,& ! Calculate aerosol water concentration limited by specific humidity (otherwise eq. is assumed)
             lgh2o      = .FALSE.,& ! Update water vapor concentration (H2O tracer only if defined in H2O tracer coupling)
             lshum      = .FALSE.,& ! Update specific humidity (ECHAM5)
             lclc       = .FALSE.,& ! Update cloud cover
             lacc       = .FALSE.,& ! Aerosol-cloud-coupling (diag only, if not LCLC,LCLWC,LCIWC,LCDNC,LICNC)
             lclwc      = .FALSE.,& ! Update cloud liquid water concentration (sets lcdnc = t)
             lciwc      = .FALSE.,& ! Update cloud ice water concentration    (sets licnc = t)
             lcdnc      = .FALSE.,& ! Calculate Cloud Droplet Number Concentration
             licnc      = .FALSE.,& ! Calculate Ice Crystal Number Concentration
             lgas       = .FALSE.,& ! Update gas phase chemisty (H2SO4, HNO3, HCl, NH3)
             laerosol   = .FALSE.,& ! Update aerosol phase chemisty (NH4+, NO3-, Cl-)
                                    ! SO4, SS, DU, OC, BC are not updated
                                    ! (concentrations of non-volatile, particulate matter remain unchanged)
             lnumber    = .FALSE.,& ! Update/calculate aerosol numbers
             lwetrad    = .FALSE.,& ! Update/calculate wet aerosol radius
             ldryrad    = .FALSE.,& ! Update/calculate dry aerosol radius
             ldrydens   = .FALSE.   ! Update/calculate dry aerosol density
             ! <--- Interface specific options
             !
             ! 3. Calculate aerosol composition and hygroscopic growth for selected modes:
  INTEGER, PUBLIC :: &  !  (to be defined in gmxe.nml)
             nlowermode = 1,      & ! Lowest mode (min. mode number = 1)
             nuppermode = 7         ! Uppest mode (max. mode number = 7) for an M7 setup
             !
             ! 4. Choice of aerosol composition module
             !
  INTEGER, PUBLIC :: &  !  (to be defined in gmxe.nml)
             neqm       = 0,      & ! Thermodynamic module: 0 = None (bulk hygroscopic growth only)
                                    !                       1 = EQSAM4CLIM
                                    !                       2 = ISORROPIA
                                    !                     [-1 = no (gmxe_thermo), but gmxe core]
                                    ! <--- Choice of the nucleation scheme
             nnucl      = 0         ! 0 = nucleation test, 1 = Vehkamaeki (2002), 2 = Kulmala (1998)

  LOGICAL, PUBLIC :: l_aerchem = .FALSE. ! Calculate explicit aerosol chemistry
                                     ! using a kpp mechanism and reaction system
  INTEGER, PUBLIC :: lmode = 0 ! lower mode boundary for aerchem
  INTEGER, PUBLIC :: umode = 0 ! upper mode boundary for aerchem
  LOGICAL, PUBLIC, POINTER :: AERCHEM(:) => NULL()

!mz_ap_20160929+
  LOGICAL, PUBLIC :: l_oracle =.FALSE.  ! calculate organics with ORACLE (interactions of oracle defined aerosols)
!mz_ap_20160929-

  LOGICAL, PUBLIC :: l_oc_aging =.FALSE.  ! calculate OC aging (oxidation by OH)

  LOGICAL, PUBLIC :: l_soa = .FALSE.   ! use the SOA sub-submodel for
                                       ! explicit SOA calculations
  PUBLIC :: lmode_soa  ! lower mode boundary for soa
  PUBLIC :: umode_soa  ! upper mode boundary for soa
  LOGICAL, PUBLIC, POINTER :: LSOA(:) => NULL()

  ! mz_dk_20120119+
  ! PASSIVE AEROSOLS
  LOGICAL, PUBLIC :: l_passive_aer = .FALSE. ! use passive aerosols
  INTEGER, PUBLIC :: num_pa = 0 ! # of passive aerosols
  INTEGER, PUBLIC :: pamode1 = 0 ! lower mode boundary for passive aerosol
  INTEGER, PUBLIC :: pamode2 = 0 ! upper mode boundary for passive aerosol
  ! mz_dk_20120119-
!_________________________________________________________________________________________________________________________________
  LOGICAL,PUBLIC  :: &

      ! namelist switches for thermodynamics (bith ISORROPIA2 and EQSAM4CLIM

      ldry    = .FALSE.,& ! Force aerosol particles to be dry (no aerosol water)
      lhyster = .FALSE.   ! Inclusion of hysteresis effect (includes solids/metastable)

  LOGICAL,PARAMETER, PUBLIC :: &
      ! Estimate aerosol water for bulk species
      lbulk_aerosol_water = .TRUE.,&
      ! If condensation is not limited by gas phase diffusion,
      ! distribute gases equally between both aerosol and gas
      ! phase (lcond_gas_phase=F) or allow complete condensation
      ! (if thermodynamically possible)
      lcond_gas_phase     = .FALSE., &
      ! calculate subgridscale relative humidity (experimental)
      lsubrh              = .FALSE., &
      ! calculate dconc as in M7
      ldconc_m7           = .FALSE.


!_________________________________________________________________________________________________________________________________
  !
  !--- 3.3) Computational constants:
  !
  REAL(dp), PUBLIC,PARAMETER :: FACp0=1._dp
  REAL(dp), PUBLIC,PARAMETER :: FACp1=10._dp
  REAL(dp), PUBLIC,PARAMETER :: FACp2=100._dp
  REAL(dp), PUBLIC,PARAMETER :: FACp3=1000._dp
  REAL(dp), PUBLIC,PARAMETER :: FACp4=10000._dp
  REAL(dp), PUBLIC,PARAMETER :: FACp5=1.e+5_dp
  REAL(dp), PUBLIC,PARAMETER :: FACp6=1.e+6_dp
  REAL(dp), PUBLIC,PARAMETER :: FACp9=1.e+9_dp
  REAL(dp), PUBLIC,PARAMETER :: FACp10=1.e+10_dp
  REAL(dp), PUBLIC,PARAMETER :: FACp12=1.e+12_dp
  REAL(dp), PUBLIC,PARAMETER :: FACm1=0.1_dp
  REAL(dp), PUBLIC,PARAMETER :: FACm2=0.01_dp
  REAL(dp), PUBLIC,PARAMETER :: FACm3=0.001_dp
  REAL(dp), PUBLIC,PARAMETER :: FACm4=0.0001_dp
  REAL(dp), PUBLIC,PARAMETER :: FACm5=1.e-5_dp
  REAL(dp), PUBLIC,PARAMETER :: FACm6=1.e-6_dp
  REAL(dp), PUBLIC,PARAMETER :: FACm7=1.e-7_dp
  REAL(dp), PUBLIC,PARAMETER :: FACm8=1.e-8_dp
  REAL(dp), PUBLIC,PARAMETER :: FACm9=1.e-9_dp
  REAL(dp), PUBLIC,PARAMETER :: FACm10=1.e-10_dp
  REAL(dp), PUBLIC,PARAMETER :: FACm12=1.e-12_dp
  REAL(dp), PUBLIC,PARAMETER :: FACm15=1.e-12_dp
  REAL(dp), PUBLIC,PARAMETER :: sqrt2=1.414213562373_dp
  !
  !--- 3.4) Physical constants:
  !
REAL(dp),PUBLIC,PARAMETER ::  ZERO = 0._dp,      & ! REAL zero
                              REALZERO=1.E-40_dp,& ! REAL values below REALZERO are assumed zero
                              TINYAW=1.E-10_dp,  & ! limit for aw parameterization (0 < aw > 1)
                              TINYCO=tiny(0._dp),& ! threshold value for concentrations [mol]
!                             TINYCO=1.E-15_dp,  & ! threshold value for concentrations [mol]
                              TINYRH=1.E-5_dp,   & ! RH limit
                              NaN=-1.E+34_dp       ! Undefined (missing value of ferret)
REAL(dp), PARAMETER,PUBLIC :: zeps=1.e-20_dp       ! EPSILON(1._dp) (determined during initialization)
REAL(dp), PARAMETER        :: bk     = 1.38e-16_dp  ! Bolzman constant [erg / K]
REAL(dp), PARAMETER        :: rerg   = 8.314E+7_dp      ! Ideal gas constant [erg. K-1 mole-1]
REAL(dp), PARAMETER        :: r_kcal = 1.986E-3_dp      ! Ideal gas constant [kcal K-1 mole-1]
REAL(dp), PARAMETER        :: Tm=273.15_dp         ! melting   temperature [K]
REAL(dp), PARAMETER        :: Ti=90.00_dp          ! max. freezing depression [K]
REAL(dp), PARAMETER        :: Tc=1.86_dp           ! cryoscopic constant for water [K l mol-1]
REAL(dp), PARAMETER        :: RT0=R_gas*T0*FACm3   ! Gas constant and standard state temperature [J mol-1 K-1] * [K] -> [KJ mol-1]
REAL(dp), PARAMETER        :: MVAIR=0.02477412_dp  ! molar volume air [m3/mol] at STD (1013.15 hPa, 298.15K)
REAL(dp), PARAMETER        :: Kw=1.010e-14_dp      ! autodiss.  constant for water [mol2 kg-2]
REAL(dp), PARAMETER        :: rwatermass=1000._dp  ! reference mass of water [g]
REAL(dp), PARAMETER        :: Pw=3167.0_dp         ! saturation vapor pressure of water [Pa] at T=298K
REAL(dp), PARAMETER        :: N2=2._dp             ! cation-anion pair [-]
REAL(dp), PARAMETER        :: WMOL=rwatermass/Mw   ! moles of H2O [mol kg-1(H2O)]
REAL(dp), PUBLIC,PARAMETER :: RHMIN=0.001_dp       ! minimum fractional relative humidity [-]
REAL(dp), PUBLIC           :: RHMAX=1._dp-TINYRH   ! maximum fractional relative humidity [-]
!box REAL(dp), PUBLIC,PARAMETER :: RHMAX=1._dp-TINYRH   ! maximum fractional relative humidity
REAL(dp), PUBLIC,PARAMETER :: TMIN=150._dp,TMAX=350._dp ! maximum temperature (lookup tables) [K]
REAL(dp), PUBLIC,PARAMETER :: FACp9divMw = FACp9/Mw     ! EQSAM3 core [mol ng-1]
!_________________________________________________________________________________________________________________________________
#ifndef _BOX
INTEGER,  PARAMETER        :: ndat=10002             ! 10002 = default 3D (SX6/IBM)
#else
INTEGER,  PARAMETER        :: ndat=1002              ! 1002  = Box model (test)
!INTEGER,  PARAMETER        :: ndat=100002           ! extended for huge box model data sets
                                                     ! (742 => T = 298.1K)
#endif
!_________________________________________________________________________________________________________________________________

  !
  !--- 3.5) Assumed parameters:
  !
  !        Assumed mass of an nucleated sulfate particle [molecules]
  REAL(dp) :: critn=100._dp
  !        Factor that limits the condensation of sulfate to fmax
  !        times the available sulfate in the gas phase [1]. (gmxe_dgas)
  REAL(dp), PARAMETER :: fmax=0.95_dp
  !         Assumed required layer thickness of
  !         sulfate to transfer an insoluble
  !         particle to a soluble mode. It is
  !         given in units of layers of
  !         monomolecular sulfate. Determines the
  !         transfer rate from insoluble to soluble modes.
!  REAL(dp), PARAMETER :: cLayerThickness = 1.0_dp
!  REAL(dp), PARAMETER :: cLayerThickness = 3.0_dp
!  REAL(dp), PARAMETER :: cLayerThickness = 5.0_dp
  REAL(dp), PARAMETER :: cLayerThickness = 10.0_dp

  REAL(dp), PARAMETER :: csurf_molec = 2.39E-15_dp   ! Average cross-section
                                                     ! of a single H2SO4 molecule [cm+2]

! GMXe help parameters
  INTEGER,PUBLIC, POINTER    ::  nwh2o (:)    => NULL()       ! H2O   indices
  INTEGER,PUBLIC, POINTER    ::  nh2so4(:)    => NULL()       ! H2SO4 indices
  INTEGER,PUBLIC, POINTER    ::  nso42m(:)    => NULL()       ! SO42m indices
  INTEGER,PUBLIC, POINTER    ::  nhso4m(:)    => NULL()       ! SO42m indices
  INTEGER,PUBLIC, POINTER    ::  nhp(:)       => NULL()       ! Hp    indices
  INTEGER,PUBLIC, POINTER    ::  nohm(:)      => NULL()       ! OHm   indices
  INTEGER,PUBLIC, POINTER    ::  nnum(:)      => NULL()       ! number index


! setting up species structure for thermodynamics
! first string contains the name,
! second the modes in which the species are allowed to occur

  INTEGER, PUBLIC, PARAMETER :: NGAS_MAX = 30
  INTEGER, PUBLIC, PARAMETER :: Ncations_MAX = 30
  INTEGER, PUBLIC, PARAMETER :: Nanions_MAX = 30
  INTEGER, PUBLIC, PARAMETER :: Nsolutes_MAX = 30
  INTEGER, PUBLIC, PARAMETER :: Nbulk_MAX = 30
  CHARACTER(LEN=1024), PUBLIC :: CASK_gases(ngas_max,2)      = ''
  CHARACTER(LEN=1024), PUBLIC :: CASK_anions(nanions_max,2) = ''
  CHARACTER(LEN=1024), PUBLIC :: CASK_cations(ncations_max,2) = ''
  CHARACTER(LEN=1024), PUBLIC :: CASK_solutes(nsolutes_max,2)= ''
  CHARACTER(LEN=1024), PUBLIC :: CASK_bulk(nbulk_max,2)      = ''
  CHARACTER(LEN=1024), PUBLIC :: sigma_nml, crdiv_nml, cmodes_nml

  INTEGER, POINTER, PUBLIC           :: nmodes(:)   => NULL()
  CHARACTER(LEN=2), POINTER, PUBLIC  :: cmodes(:)   => NULL()
  LOGICAL, POINTER                   :: locoagmask(:,:) => NULL()
!_________________________________________________________________________________________________________________________________
!--- 4) Tracer mapping and indices defintions:
!_________________________________________________________________________________________________________________________________
!--- 4.1) Aerosol compounds (array paerml):
  INTEGER, PUBLIC   ::  &
                              nbulk_eq      = 7                , &              ! number of bulk (unkown/insoluble) species (7)
                              ncomp= 0,     &! number of species per aerosol mode (56)
!                              ncomp=ngas+ncations+nanions+nbulk_eq+nsol,     &! number of species per aerosol mode (56)
                              naertot_all= 0
!                              naertot_all= ncomp * (nmod+1)                  ! 448 maximum number of all species
!INTEGER,PUBLIC             :: mcomp(0:nmod), naertot, naero                  ! variable array boundary (for run time memory)
  INTEGER, PUBLIC, POINTER   :: mcomp(:) => NULL()
  INTEGER, PUBLIC            :: naertot
!--- 4.2) Diagnostic aerosol properties (array paerosol):
  INTEGER,PUBLIC,PARAMETER   :: iTT        =  1         , &  ! T  [K]
                              iRH        =  2           , &  ! RH [0-1]
                              iPX        =  3           , &  ! aerosol water history            [0=solid, else=wet]
                              iZIONIC    =  4           , &  ! ionic strength (aq)              [mol/kg]
                              irho       =  5           , &
                              iVOL       =  6           , &  ! bulk volume                      [ucm3 m-3 (air)]
                              iPH        =  7           , &  ! aerosol pH (from TD)             [log H+]
                              isPM       =  8           , &  ! total solid  matter              [umol m-3 (air)]
                              iaPM       =  9           , &  ! total liquid matter              [umol m-3 (air)]
                              iPMt       = 10           , &  ! total PM (liquids & solids)      [ug   m-3 (air)]
                              iPMs       = 11           , &  ! total PM (solids)                [ug   m-3 (air)]
                              irhdcr_sum = 12           , &  ! crystallization RH of mixed solution [-]
                              irhdms_sum = 13           , &  ! deliquescence   RH of mixed solution [-]
                              iRHcr      = 14           , &  ! lowest crystallization RH in mixed solution [-]
                              iRHD       = 15           , &  ! lowest deliquescence   RH in mixed solution [-]
                              iWH2O      = 16           , &  ! aerosol Water  (aq)              [ug m-3 (air)]
                              iGF        = 17           , &  ! hygroscopic growth factor        [-]
                              iTice      = 18           , &  !- freezing depression by aerosols [K]
!   additional diagnostic output - example ...
                              inewph     = 19           , &  ! new calculation of pH (after all aerosol processes)
                              iKappa     = 20           , &  ! Kappa vaule for water uptake
                              iKa_water  = 21           , &  ! Volume aerosol water (calculated from Kappa)
                              iSc        = 22           , &  ! Critical supersaturation (from Kappa)
                              iCCN2      = 23           , &  ! CCN at 0.2% supersaturation (from Kappa)
                              iCCN4      = 24           , &  ! CCN at 0.4% supersaturation (from Kappa)
                              iKappa_insol = 25         , &  ! Kappa value including insolubles
                              iKa_vol      = 26         , &  ! Volume aerosol water (calculated from Kappa)
                              iKa_vol_insol = 27        , &  ! Volume aerosol water (calculated from Kappa_insol)
                              ioc_water = 28            , &  ! water uptake by organic aerosols
                              ioc_kappa = 29            , &  ! kappa value for organics (in case of ageing)
                              ifsulf     = 30            , &  ! fraction of sulphate from dissociated H2SO4
                              naerodiag  = ifsulf               ! max. number of diagnostic aerosol properties

  CHARACTER(LEN=60),DIMENSION(:,:), ALLOCATABLE,PUBLIC :: cdiagaer             ! GMXe diagnostic output

  INTEGER, PARAMETER, PUBLIC   :: naero = naerodiag

!_______________________________________________________________________________________________
!_________________________________________________________________________________________________________________________________
!_________________________________________________________________________________________________________________________________
  !
  !
  !--- 6) Conversion factors for lognormal particle size distributions: -------------
  !       Calculated in gmxe_initialize.
  !
  ! Conversion factor: count median radius to radius of average surface
!  REAL(dp),  PUBLIC          :: cmr2ras(nmod)
  REAL(dp),  PUBLIC, POINTER :: cmr2ras(:)        => NULL()
  !
  ! Conversion factor: count median radius to mass mean radius
!  REAL(dp),  PUBLIC          :: cmr2mmr(nmod)
  REAL(dp),  PUBLIC, POINTER :: cmr2mmr(:)        => NULL()
  !
  ! Conversion factor: count median radius to mass median radius
!  REAL(dp),  PUBLIC          :: cmedr2mmedr(nmod)
  REAL(dp),  PUBLIC, POINTER :: cmedr2mmedr(:)        => NULL()
  !
  ! Conversion factor: count median radius to radius of average mass
!  REAL(dp),  PUBLIC          :: cmr2ram(nmod)
  REAL(dp),  PUBLIC, POINTER :: cmr2ram(:)        => NULL()
  !
  ! Conversion factor: radius of average mass to count median radius
!  REAL(dp),  PUBLIC          :: ram2cmr(nmod)
  REAL(dp),  PUBLIC, POINTER :: ram2cmr(:)        => NULL()
  ! Conversion factor:  LOG(sigma(nmod))
!  REAL(dp),  PUBLIC          :: sigma_ln(nmod)
  REAL(dp),  PUBLIC, POINTER :: sigma_ln(:)        => NULL()
  ! Conversion factor: EXP(1.5*LOG(sigma(nmod))**2)
!  REAL(dp),  PUBLIC          :: sigma_exp_ln(nmod)
  REAL(dp),  PUBLIC, POINTER :: sigma_exp_ln(:)        => NULL()
!_______________________________________________________________________________
  !--- 7) Assumed thresholds for occurence of specific quantities: -------------
  !
  REAL(dp) :: cmin_aerml     = 1.e-12_dp, &! Aerosol bulk mass [umol m-3 (air)]
              cmin_aernl     = 1.e-7_dp,  &! Aerosol number    [N   cm-3 (air)]
              cmin_avol      = 1.e-14_dp, &! Aerosol volume    [cm3 cm-3 (air)]
              cmin_nucl      = 1.e-15_dp, &! H2SO4/H2O concen. [umol m-3 (air)]
              cmin_epsilon   = 1.e-18_dp, &! Minium value
              cmin_nuclmolec = 1.e+4_dp    ! H2SO4 concen. [molecules cm-3(air)]
!_____________________________________________________________________________

  ! matrix for moving of particles in coagulation
!  REAL(dp) :: xmov(nmod,nmod)
!  INTEGER  :: idest(nmod,nmod), ndiff

  REAL(dp), POINTER :: xmov(:,:)  => NULL()
  INTEGER,  POINTER :: idest(:,:) => NULL()

!
! logical switch if an equilibrium model is used for that mode
  LOGICAL, POINTER  :: lleqm(:)   => NULL()

  ! level index of 300 hPa, determined in init_cpl
  INTEGER, PUBLIC :: k300
!=============================================================================


!=============================================================================
CONTAINS

!=============================================================================

 SUBROUTINE gmxe_read_nml_ctrl(status, iou)

    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    SAVE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    !--- Local variables:

    NAMELIST /CTRL/lgmxe,     loutput,    lmass_diag, lstrat,     lpsc,       lnucl,     &
                   lcoat,     lsize,      lcoag,      lcond,      ladyn,      lah2o,     &
                   lgh2o,     lshum,      lclc,       lacc,       lclwc,      lciwc,     &
                   lcdnc,     licnc,      lgas,       laerosol,                          &
                   lnumber,   lwetrad,    ldryrad,    ldrydens,                          &
                   nlowermode, nuppermode, neqm,      nnucl,      l_aerchem,  l_oracle, l_oc_aging, &
                   l_passive_aer, l_soa !mz_dk_20120119

    NAMELIST /CTRL_GMXE_DISC/ &
                   nmod, sigma_nml, crdiv_nml, cmodes_nml

    NAMELIST /CTRL_GMXE_SPECIES/ &
                   CASK_gases, CASK_anions, CASK_cations, CASK_solutes, CASK_bulk

    NAMELIST /CTRL_GMXE_TD/ &
                   ldry,  lhyster

    NAMELIST /CTRL_GMXE_AERCHEM/ &
                   lmode , umode

    ! mz_dk_20120119+
    ! PASSIVE AEROSOL NAMELIST
    NAMELIST /CTRL_GMXE_PASSIVE/ &
                   num_pa, pamode1, pamode2
    ! mz_dk_20120119-
    NAMELIST /CTRL_GMXE_SOA/ &
                   lmode_soa, umode_soa, L_GASPREC, L_OXI, EXCL_STR_SOA

    CHARACTER(LEN=*), PARAMETER :: substr = 'gmxe_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    CALL start_message(TRIM(modstr),'INITIALISATION', substr)

    ! INITIALIZE

    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    !--- 1) Read CTRL namelist:
    ! for processes
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

    !  for size discretisation
    CALL read_nml_open(lex, substr, iou, 'CTRL_GMXE_DISC', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_GMXE_DISC, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_GMXE_DISC', modstr)
    IF (fstat /= 0) status = 1
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL_GMXE_SPECIES', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_GMXE_SPECIES, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_GMXE_SPECIES', modstr)
    IF (fstat /= 0) status = 1
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

    !--- 2) Read CTRL_GMXE_TD namelist:

    CALL read_nml_open(lex, substr, iou, 'CTRL_GMXE_TD', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_GMXE_TD, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_GMXE_TD', modstr)
    IF (fstat /= 0) status = 1
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

    !--- 2) Read CTRL_GMXE_AERCHEM namelist:

    CALL read_nml_open(lex, substr, iou, 'CTRL_GMXE_AERCHEM', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_GMXE_AERCHEM, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_GMXE_AERCHEM', modstr)
    IF (fstat /= 0) status = 1
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

    !--- 4) Read CTRL_GMXE_SOA namelist:

    CALL read_nml_open(lex, substr, iou, 'CTRL_GMXE_SOA', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_GMXE_SOA, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_GMXE_SOA', modstr)
    IF (fstat /= 0) status = 1
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

    !--- 5) Read CTRL_GMXE_PASSIVE namelist:

    CALL read_nml_open(lex, substr, iou, 'CTRL_GMXE_PASSIVE', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_GMXE_PASSIVE, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_GMXE_PASSIVE', modstr)
    IF (fstat /= 0) status = 1
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

    OPEN(unit=iou,file='GMXE_TDMODEL_VERSION.txt',status='unknown',form='formatted')

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    ! CHECK TRACER INITIALIZATION

    WRITE(*,*) '    Initialization of aerosol dynamics and composition module GMXe'
    WRITE(*,*) ' '
    WRITE(*,*) 'DIAGNOSE NAMELISTS'
    WRITE(*,*) ' '
    WRITE(*,*) ' 1) CTRL - GMXe general options:'
    WRITE(*,*) ' '
    WRITE(*,*) '              lgmxe      = ', lgmxe
    WRITE(*,*) '              loutput    = ', loutput
!    WRITE(*,*) '              ldebug     = ', ldebug
    WRITE(*,*) '              lmass_diag = ', lmass_diag
    WRITE(*,*) '              lstrat     = ', lstrat
    WRITE(*,*) '              lpsc       = ', lpsc
    WRITE(*,*) '              lnucl      = ', lnucl
    WRITE(*,*) '              lcoat      = ', lcoat
    WRITE(*,*) '              lsize      = ', lsize
    WRITE(*,*) '              lcoag      = ', lcoag
    WRITE(*,*) '              lcond      = ', lcond
    WRITE(*,*) '              ladyn      = ', ladyn
    WRITE(*,*) '              lah2o      = ', lah2o
    WRITE(*,*) '              lgh2o      = ', lgh2o
    WRITE(*,*) '              lshum      = ', lshum
    WRITE(*,*) '              lclc       = ', lclc
    WRITE(*,*) '              lacc       = ', lacc
    WRITE(*,*) '              lclwc      = ', lclwc
    WRITE(*,*) '              lciwc      = ', lciwc
    WRITE(*,*) '              lcdnc      = ', lcdnc
    WRITE(*,*) '              licnc      = ', licnc
    WRITE(*,*) '              lgas       = ', lgas
    WRITE(*,*) '              laerosol   = ', laerosol
    WRITE(*,*) '<--- GMXe stream:'
    WRITE(*,*) '              lnumber    = ', lnumber
    WRITE(*,*) '              lwetrad    = ', lwetrad
    WRITE(*,*) '              ldryrad    = ', ldryrad
    WRITE(*,*) '              ldrydens   = ', ldrydens
    WRITE(*,*) ' '
    WRITE(*,*) '<--- Calculate aerosol composition and hygroscopic growth for selected modes:'
    WRITE(*,*) ' '
    WRITE(*,*) '  NLOWERMODE = ', NLOWERMODE  ,'= Lowest mode (min. mode number = 1; nucleation mode,   soluble)'
    WRITE(*,*) '  NUPPERMODE = ', NUPPERMODE  ,'= Uppest mode (max. mode number = 7; coarse     mode, insoluble)'
    WRITE(*,*) ' '
    WRITE(*,*) '<--- Choice of aerosol composition module'
    WRITE(*,*) '              neqm       = ', neqm
    WRITE(*,*) ' '
    IF (lgmxe .AND. neqm==-1) THEN
    WRITE(*,*) ' >>> Caution: (gmxe_thermo) switched off <<<'
    ELSE IF (lgmxe .AND. neqm==0) THEN
    WRITE(*,*) ' >>> Caution: No thermodynamic module will be used! <<<'
    WRITE(*,*)               'bulk hygroscopic growth only'
    ELSE IF (lgmxe .AND. neqm==1)  THEN
    WRITE(*,*) '              EQSAM4CLIM     = ', modver,' (vectorized MESSy Version)'
    WRITE(*,*) '                         => Metzger et al., 2016'
    WRITE(*,*)               'EQSAM4CLIM'
    WRITE(iou,*)             'EQSAM4CLIM'
    ELSE IF (lgmxe .AND. neqm==2)  THEN
    WRITE(*,*) '              ISORROPIA   = ', 2.1,' (MESSy Version)'
    WRITE(*,*)               'ISORROPIA2'
    WRITE(iou,*)             'ISORROPIA2'
    ELSE IF (lgmxe .AND. (neqm/=-1.OR. neqm/=0  .OR. neqm/=1 .OR. neqm/=2)) THEN
             lgmxe=.false.
    WRITE(*,*) ' '
    WRITE(*,*) ' >>> Caution: No proper thermodynamic module selected <<< '
    WRITE(*,*) ' >>> Aaerosol composition interface switched off !    <<< '
    WRITE(*,*) ' '
    ELSE
    WRITE(*,*) ' '
    WRITE(*,*) ' >>> Caution: GMXe interface (gmxe_main) switched off <<<'
    WRITE(*,*) ' '
    END IF
    WRITE(*,*) '<--- Choice of nucleation module'
    IF (LNUCL .AND. nnucl==0) THEN
       WRITE(*,*) ' '
       WRITE(*,*) '              nnucl      = ', nnucl
       WRITE(*,*) '                         => nucleation test'
    ELSE IF (LNUCL .AND. nnucl==1) THEN
       WRITE(*,*) ' '
       WRITE(*,*) '              nnucl      = ', nnucl
       WRITE(*,*) '                         => Vehkamaeki et al., 2002'
    ELSE IF (LNUCL .AND. nnucl==2)  THEN
       WRITE(*,*) ' '
       WRITE(*,*) '              nnucl      = ', nnucl
       WRITE(*,*) '                         => Kulmala et al., 1998'
    ELSE IF (LNUCL .AND. (nnucl/=0 .OR. nnucl/=1 .OR. nnucl/=2)) THEN
       LNUCL=.FALSE.
       WRITE(*,*) ' '
       WRITE(*,*) 'Nucleation requested but none selected. Nucleation switched off !'
       WRITE(*,*) ' '
    END IF
    IF(lcoat .AND. neqm < 1)  THEN
       lcoat=.FALSE.
       WRITE(*,*) 'Coating switched off (neqm < 1)!'
       WRITE(*,*) ' '
    END IF
    IF(.NOT.lacc) THEN
       WRITE(*,*) 'Aerosol-cloud-coupling switched off !'
       WRITE(*,*) 'Maximum relative humidity for aerosol water uptake limited to RH=95%'
       WRITE(*,*) ' '
    END IF
    CLOSE(iou)
    WRITE(*,*) ' '


    CALL end_message(TRIM(modstr),'INITIALISATION', substr)

  END SUBROUTINE gmxe_read_nml_ctrl

!==============================================================================

  SUBROUTINE gmxe_initialize_species(pe,p_io)

    USE messy_gmxe_soa,      ONLY: names, nsoa, mw_soa => MW
    USE messy_gmxe_oc_aging, ONLY: num_wsoc, kappa
    USE messy_main_tools,    ONLY: str2num, strcrack

    INTEGER  :: pe, p_io, status
    INTEGER  :: jm, ji, jt, num
    INTEGER  :: counter, counter2
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER     :: outstring(:) => NULL()
    LOGICAL  :: found
    CHARACTER(LEN=2)            :: str_wsoc
! gases
    counter = 0
    DO ji=1,ngas_max
      IF ( TRIM(ADJUSTL(CASK_gases(ji,1))) /= '' ) THEN
        counter = counter + 1
      ENDIF
    END DO
    ngas = counter
    ALLOCATE(td%gas(ngas))

    counter = 0
    DO ji=1,ngas_max
      IF ( TRIM(CASK_gases(ji,1)) /='' ) THEN
        counter = counter + 1
        td%gas(counter)%name=TRIM(CASK_gases(ji,1))
        ALLOCATE(td%gas(counter)%ltreat(0:nmod))
        ALLOCATE(td%gas(counter)%aerml_idx(0:nmod))
        ALLOCATE(td%gas(counter)%caccgas(nmod))
        td%gas(counter)%ltreat(0:nmod)    = .FALSE.
        td%gas(counter)%aerml_idx(0:nmod) = 0

        call strcrack(CASK_gases(ji,2), ';', outstring,  &
          counter2)
        DO jm=0,nmod
          DO jt=1,counter2
            call str2num(TRIM(outstring(jt)),num)
            IF (jm == num) THEN
              td%gas(counter)%ltreat(jm) = .TRUE.
              EXIT
            ENDIF
          END DO
        END DO
      ENDIF
    END DO

!anions
    counter = 0
    DO ji=1,nanions_max
      IF ( TRIM(CASK_anions(ji,1)) /='' ) THEN
        counter = counter + 1
      ENDIF
    END DO
    nanions = counter
    ALLOCATE(td%anion(nanions))

    counter = 0
    DO ji=1,nanions_max
      IF ( TRIM(CASK_anions(ji,1)) /='' ) THEN
        counter = counter + 1
        td%anion(counter)%name=TRIM(CASK_anions(ji,1))
        ALLOCATE(td%anion(counter)%ltreat(nmod))
        ALLOCATE(td%anion(counter)%aerml_idx(nmod))
        td%anion(counter)%ltreat(1:nmod)    = .FALSE.
        td%anion(counter)%aerml_idx(1:nmod) = 0

        call strcrack(CASK_anions(ji,2), ';', outstring,  &
          counter2)
        DO jm=1,nmod
          DO jt=1,counter2
            call str2num(TRIM(outstring(jt)),num)
            IF (jm == num) THEN
              td%anion(counter)%ltreat(jm) = .TRUE.
              EXIT
            ENDIF
          END DO
        END DO
      ENDIF
    END DO

!cations
    counter = 0
    DO ji=1,ncations_max
      IF ( TRIM(CASK_cations(ji,1)) /='' ) THEN
        counter = counter + 1
      ENDIF
    END DO
    ncations = counter
    ALLOCATE(td%cation(ncations))

    counter = 0
    DO ji=1,ncations_max
      IF ( TRIM(CASK_cations(ji,1)) /='' ) THEN
        counter = counter + 1
        td%cation(counter)%name=TRIM(CASK_cations(ji,1))
        ALLOCATE(td%cation(counter)%ltreat(nmod))
        ALLOCATE(td%cation(counter)%aerml_idx(nmod))
        td%cation(counter)%ltreat(1:nmod)    = .FALSE.
        td%cation(counter)%aerml_idx(1:nmod) = 0

        call strcrack(CASK_cations(ji,2), ';', outstring,  &
          counter2)
        DO jm=1,nmod
          DO jt=1,counter2
            call str2num(TRIM(outstring(jt)),num)
            IF (jm == num) THEN
              td%cation(counter)%ltreat(jm) = .TRUE.
              EXIT
            ENDIF
          END DO
        END DO
      ENDIF
    END DO

!solutes
    counter = 0
    DO ji=1,nsolutes_max
      IF ( TRIM(CASK_solutes(ji,1)) /='' ) THEN
        counter = counter + 1
      ENDIF
    END DO
    nsolutes = counter
    ALLOCATE(td%solute(nsolutes))

    counter = 0
    DO ji=1,nsolutes_max
      IF ( TRIM(CASK_solutes(ji,1)) /='' ) THEN
        counter = counter + 1
        td%solute(counter)%name=TRIM(CASK_solutes(ji,1))
        ALLOCATE(td%solute(counter)%ltreat(0:nmod))
        ALLOCATE(td%solute(counter)%aerml_idx(0:nmod))
        td%solute(counter)%ltreat(1:nmod)    = .FALSE.
        td%solute(counter)%aerml_idx(0:nmod) = 0
        call strcrack(CASK_solutes(ji,2), ';', outstring,  &
          counter2)
        DO jm=0,nmod
          DO jt=1,counter2
            call str2num(TRIM(outstring(jt)),num)
            IF (jm == num) THEN
              td%solute(counter)%ltreat(jm) = .TRUE.
              EXIT
            ENDIF
          END DO
        END DO
      ENDIF
    END DO

!bulk
    counter = 0
    DO ji=1,nbulk_max
      IF ( TRIM(CASK_bulk(ji,1)) /='' ) THEN
        counter = counter + 1
      ENDIF
    END DO
    nbulk = counter

    IF (L_OC_AGING) THEN
      found = .FALSE.
      DO ji=1,nbulk_max
        IF ( TRIM(CASK_bulk(ji,1)) =='OC' ) THEN
          found = .TRUE.
          EXIT
        END IF
      END DO
      IF (FOUND) THEN
        counter = counter + num_wsoc
      ELSE
        L_OC_AGING = .FALSE.
        IF (pe == p_io) &
          write (*,*) "Warning! No OC bulk species found, therefore OC aging switched off!"
      END IF
      nbulk = counter
    END IF

    IF (L_SOA) THEN
      counter = counter + nsoa
      nbulk = counter
    END IF

    ALLOCATE(bulk(nbulk))

    DO ji=1,nbulk
      ALLOCATE(bulk(ji)%ltreat(0:nmod))
      ALLOCATE(bulk(ji)%aerml_idx(nmod))
      bulk(ji)%ltreat(0:nmod)    = .FALSE.
      bulk(ji)%aerml_idx(1:nmod) = 0
    ENDDO


    counter = 0
    DO ji=1,nbulk_max
      IF ( TRIM(CASK_bulk(ji,1)) /='' ) THEN
        counter = counter + 1
        bulk(counter)%name=TRIM(CASK_bulk(ji,1))

        call strcrack(CASK_bulk(ji,2), ';', outstring,  &
          counter2)
        DO jm=1,nmod
          DO jt=1,counter2
            call str2num(TRIM(outstring(jt)),num)
            IF (jm == num) THEN
              bulk(counter)%ltreat(jm) = .TRUE.
              EXIT
            ENDIF
          END DO
        END DO
      ENDIF
    END DO

    IF (L_OC_AGING) THEN
      DO ji=1,num_wsoc
        counter = counter + 1
        IF (ji < 10) then
          write (str_wsoc,'(A,I1)') "0",ji
        ELSE
          write (str_wsoc,'(I2)') ji
        END IF
        bulk(counter)%name='WSOC'//str_wsoc
        bulk(counter)%molmass = 12.01_dp
        bulk(counter)%density = 2.0_dp
        bulk(counter)%kappa   = kappa(ji)
        bulk(counter)%l_oc    = .TRUE.
        DO jm=1,nsoluble
          found = .FALSE.
          DO jt=1,nbulk
            IF ( TRIM(bulk(jt)%name)=='OC' ) THEN
              IF (bulk(jt)%ltreat(jm)) found = .TRUE.
              EXIT
            ENDIF
          END DO
          IF (FOUND) THEN
            bulk(counter)%ltreat(jm) = .TRUE.
          ENDIF
        END DO
      END DO
    END IF

    IF (L_SOA) THEN
      DO ji=1,NSOA
        counter = counter + 1
        bulk(counter)%name=TRIM(names(ji))
        bulk(counter)%molmass = MW_SOA(ji)
        bulk(counter)%density = 2.0_dp
        bulk(counter)%kappa   = 0.1_dp
        bulk(counter)%l_oc    = .TRUE.
        DO jm=1,nmod
          IF (.NOT. LSOA(jm)) CYCLE
          bulk(counter)%ltreat(jm) = .TRUE.
        END DO
      END DO
    END IF

    CALL initialise_bulk_gas

  END SUBROUTINE gmxe_initialize_species
!==============================================================================
  SUBROUTINE initialise_bulk_gas

    USE messy_gmxe_mem,           ONLY: nbulk, bulk, ngas, td
    USE messy_main_constants_mem, ONLY: MNa, MCl

    TYPE bulk_inf
      CHARACTER(LEN=STRLEN_MEDIUM) :: name    = ''
      REAL(dp)                     :: kappa   = 0._dp
      REAL(dp)                     :: density = 0._dp
      REAL(dp)                     :: molmass = 0._dp
      LOGICAL                      :: l_oc    = .FALSE.
    END type bulk_inf

    TYPE(BULK_inf)    :: bulk__info(nbulk_max)
    INTEGER  :: jt, jc
    TYPE gas_inf
      CHARACTER(LEN=STRLEN_MEDIUM) :: name    = ''
      REAL(dp), DIMENSION(2)       :: caccgas = 0._dp
      CHARACTER(LEN=STRLEN_MEDIUM) :: ionname = ''
    END type gas_inf
    TYPE(gas_inf)    :: gases_info(ngas_max)


    BULK__INFO(1)%name    = "DU"
    BULK__INFO(1)%kappa   = 0._dp
    BULK__INFO(1)%density = 2.650_dp
    BULK__INFO(1)%molmass = 40.08_dp

    BULK__INFO(2)%name    = "OC"
    BULK__INFO(2)%kappa   = 0.1_dp
    BULK__INFO(2)%density = 2.0_dp
    BULK__INFO(2)%molmass = 12.01_dp
    BULK__INFO(2)%l_oc    = .TRUE.

    BULK__INFO(3)%name    = "BC"
    BULK__INFO(3)%kappa   = 0._dp
    BULK__INFO(3)%density = 2.0_dp
    BULK__INFO(3)%molmass = 12.01_dp

    BULK__INFO(4)%name    = "SS"
    BULK__INFO(4)%kappa   = 1.12_dp
    BULK__INFO(4)%density = 2.170_dp
    BULK__INFO(4)%molmass = MNa + MCl

    BULK__INFO(5)%name    = "SOA"
    BULK__INFO(5)%kappa   = 0.1_dp
    BULK__INFO(5)%density = 2.0_dp
    BULK__INFO(5)%molmass = 12.01_dp
    BULK__INFO(5)%l_oc    = .TRUE.

    DO jt = 1, nbulk_max
      DO jc = 1, nbulk
        IF (TRIM(BULK__INFO(jt)%name) == TRIM(BULK(jc)%name) ) THEN
          bulk(jc)%kappa   = BULK__INFO(jt)%kappa
          bulk(jc)%density = BULK__INFO(jt)%density
          bulk(jc)%molmass = BULK__INFO(jt)%molmass
          bulk(jc)%l_oc    = BULK__INFO(jt)%l_oc
        END IF
      END DO
    END DO


    DO jc = 1, ngas_max
      gases_info(jc)%caccgas(1)  = 1.0_dp
      gases_info(jc)%caccgas(2)  = 0.3_dp
    END DO

    gases_info(1)%name       = "HNO3"
    gases_info(1)%caccgas(1) = 0.1_dp
    gases_info(1)%caccgas(2) = 0.1_dp
    gases_info(1)%ionname    = "NO3m"

    gases_info(2)%name       = "H2SO4"
    gases_info(2)%caccgas(1) = 1.0_dp
    gases_info(2)%caccgas(2) = 0.3_dp
    gases_info(2)%ionname    = "SO4mm"

    gases_info(3)%name       = "HCl"
    gases_info(3)%caccgas(1) = 0.064_dp
    gases_info(3)%caccgas(2) = 0.064_dp
    gases_info(3)%ionname    = "Clm"

    gases_info(4)%name       = "HBr"
    gases_info(4)%caccgas(1) = 1.0_dp
    gases_info(4)%caccgas(2) = 0.3_dp
    gases_info(4)%ionname    = "Brm"

    gases_info(5)%name       = "HI"
    gases_info(5)%caccgas(1) = 1.0_dp
    gases_info(5)%caccgas(2) = 0.3_dp
    gases_info(5)%ionname    = "Im"

    gases_info(6)%name       = "NH3"
    gases_info(6)%caccgas(1) = 0.097_dp
    gases_info(6)%caccgas(2) = 0.097_dp
    gases_info(6)%ionname    = "NH4p"


    DO jc = 1, ngas
      td%gas(jc)%caccgas(1:nsoluble)      = 1.0_dp
      td%gas(jc)%caccgas(nsoluble+1:nmod) = 0.3_dp
      DO jt=1,ngas_max
        IF (TRIM(gases_info(jt)%name) == TRIM(td%gas(jc)%name) ) THEN
          td%gas(jc)%caccgas(1:nsoluble)      = gases_info(jt)%caccgas(1)
          td%gas(jc)%caccgas(nsoluble+1:nmod) = gases_info(jt)%caccgas(2)
          td%gas(jc)%ionname                  = gases_info(jt)%ionname
        END IF
      END DO
    END DO

    DO jc=1,nanions
      DO jt=1,ngas
        IF (TRIM(td%anion(jc)%name) == TRIM(td%gas(jt)%ionname) ) &
          td%anion(jc)%gas_idx=jt
      END DO
    END DO
    DO jc=1,ncations
      DO jt=1,ngas
        IF (TRIM(td%cation(jc)%name) == TRIM(td%gas(jt)%ionname) ) &
          td%cation(jc)%gas_idx=jt
      END DO
    END DO


  END SUBROUTINE initialise_bulk_gas

!==============================================================================
  SUBROUTINE gmxe_initialize_core(status)

    USE messy_main_tools,    ONLY: str2num, strcrack
    USE messy_main_constants_mem, ONLY: iouerr

    INTEGER :: jm, ji, counter, status
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER :: outstring(:) => NULL()
    CHARACTER(LEN=*), PARAMETER :: substr = 'gmxe_initialize_core'

    IF(l_io) &
    CALL start_message(TRIM(modstr),'CORE INITIALISATION', substr)

    Status = -1

    ALLOCATE(aerchem(nmod))
    AERCHEM(:) = .FALSE.
    IF (L_SOA) THEN
      ALLOCATE(LSOA(nmod))
      LSOA(:) = .FALSE.
    END IF

    ALLOCATE(sigma(nmod))
    ALLOCATE(sigma_ln(nmod))
    ALLOCATE(sigma_exp_ln(nmod))

    ALLOCATE(crdiv(nmod))
    ALLOCATE(crdiv_r3(nmod))
    ALLOCATE(crdiv_mid(nmod))


    ALLOCATE(nmodes(0:nmod))
    ALLOCATE(cmodes(nmod))
    ALLOCATE(lleqm(nmod))

    ALLOCATE(cmr2ras(nmod))
    ALLOCATE(cmr2mmr(nmod))
    ALLOCATE(cmedr2mmedr(nmod))
    ALLOCATE(cmr2ram(nmod))
    ALLOCATE(ram2cmr(nmod))

    ALLOCATE(xmov(nmod,nmod))
    ALLOCATE(idest(nmod,nmod))
    ALLOCATE(locoagmask(nmod,nmod))

    ALLOCATE(nwh2o(0:nmod))
    ALLOCATE(nh2so4(0:nmod))
    ALLOCATE(nso42m(0:nmod))
    ALLOCATE(nhso4m(0:nmod))
    ALLOCATE(nHp(0:nmod))
    ALLOCATE(nohm(0:nmod))
    ALLOCATE(nnum(nmod))

    nwh2o  = 0
    nh2so4 = 0
    nso42m = 0
    nhso4m = 0
    nhp    = 0
    nohm   = 0
    nnum   = 0

    LLEQM(:) = .FALSE.
    DO jm=nlowermode, nuppermode
      LLEQM(jm) = .TRUE.
    END DO

    CALL strcrack(sigma_nml, ';', outstring, counter)
    IF (counter /= nmod) THEN
      write(iouerr,*) "Number of modes and entries for sigma do not correspond! EXITING!"
      RETURN
    END IF
    DO ji = 1, counter
      CALL str2num(outstring(ji),sigma(ji))
    END DO
    CALL strcrack(crdiv_nml, ';', outstring, counter)
    IF (counter /= nmod) THEN
      write(iouerr,*) "Number of modes and entries for crdiv do not correspond! EXITING!"
      RETURN
    END IF
    DO ji = 1, counter
      CALL str2num(outstring(ji),crdiv(ji))
    END DO
    CALL strcrack(cmodes_nml, ';', outstring, counter)
    IF (counter /= nmod) THEN
      write(iouerr,*) "Number of modes and entries for cmodes do not correspond! EXITING!"
      RETURN
    END IF
    DO ji = 1, counter
      cmodes(ji) = TRIM(outstring(ji))
    END DO

    NSOLUBLE = NINT(NMOD/2._dp + 0.1)
    ndiff = INT(nsoluble - INT(nmod / 2._dp))

! order of modes in their treatment
! gas phase first
    NMODES(0) = 0
! GMXE previously used a size ordering
! preferring insolbule over soluble modes
    IF (NDIFF /= 0) THEN
      NMODES(1) = NDIFF
      DO jm=2,nmod,2
        NMODES(jm)   = (jm / 2) + nsoluble
        NMODES(jm+1) = ndiff + (jm / 2)
      END DO
    ELSE
      DO jm=1,nmod,2
        NMODES(jm)   = (jm/2) + 1 + nsoluble
        NMODES(jm+1) = (jm/2) + 1
      END DO
    END IF

    DO jm=1,nsoluble - 1
      crdiv_mid(jm) = sqrt(crdiv(jm+1) * crdiv(jm))
    END DO
    DO jm=nsoluble+1, nmod - 1
      crdiv_mid(jm) = sqrt(crdiv(jm+1) * crdiv(jm))
    END DO
    crdiv_mid(nsoluble) = sqrt(10._dp * crdiv(nsoluble) * crdiv(nsoluble))
    crdiv_mid(nmod) = sqrt(10._dp * crdiv(nmod) * crdiv(nmod))

    DO jm=1, nmod

       !--- 1) Calculate conversion factors for lognormal distributions:----
       !       Radius of average mass (ram) to count median radius (cmr) and
       !       vice versa. Count median radius to radius of average
       !       mass (ram).
       !       These factors depend on the standard deviation (sigma)
       !       of the lognormal distribution.
       !       (Based on the Hatch-Choate Conversins Equations;
       !       see Hinds, Chapter 4.5, 4.6 for more details.
       !       In particular equation 4.53.)
       !       Adopted from ECHAM5/M7, Philip Stier, MPI, May 2001.

       !--- Count Median Radius to Mass Median Radius:

       cmedr2mmedr(jm) = EXP(3.0_dp*(LOG(sigma(jm)))**2)

       !--- Count Median Radius to Mass Mean Radius:

!       cmr2mmr(jm) = EXP(3.5_dp*(LOG(sigma(jm)))**2)
       cmr2mmr(jm) = EXP(3._dp*(LOG(sigma(jm)))**2)

       !--- Count Median Radius to Radius of Average Mass:

       cmr2ram(jm) = EXP(1.5_dp*(LOG(sigma(jm)))**2)

       crdiv_r3(jm) = (1.e-2_dp * crdiv_mid(jm)*cmr2ram(jm))**3._dp

       !--- Radius of Average Mass to Count Median Radius:

       ram2cmr(jm) = 1._dp / cmr2ram(jm)

       !--- Count Median Radius to Radius of Average Surface:

       cmr2ras(jm) = EXP(1._dp*(LOG(sigma(jm)))**2)

       !--- 2.1) Calculate the natural logarithm of the standard deviation:

       sigma_ln(jm) = LOG(sigma(jm))

       !--- 2.2) Calculate factor for radius of average mass to the third power
       !         for emission mass to number conversion
       !         cmr2ram to the third power
        sigma_exp_ln(jm) = EXP(4.5*LOG(sigma(jm))*LOG(sigma(jm)))

    END DO

    !    Calculate coagulation moving matrix xmov(nmod,nmod)
    !    and destination mode matrix in coagulation idest(nmod,nmod)
    XMOV(:,:) = 0._dp
    idest(:,:) = 0

    ! Define coagulation matrix
    LOCOAGMASK(:,:) = .FALSE.
    DO jm=1,nmod
      DO ji=jm,nmod
        LOCOAGMASK(ji,jm) = .TRUE.
        IF ( (crdiv(jm) > 0.3e-4_dp) .AND. (ji > nsoluble) )  &
          LOCOAGMASK(ji,jm) = .FALSE.
        IF ( (jm > nsoluble) .AND. (crdiv(ji) > 0.055e-4_dp) ) &
          LOCOAGMASK(ji,jm) = .FALSE.
      END DO
    END DO

    DO jm=1,nmod
      DO ji=1,nmod
        IF ( (JM <= nsoluble .AND. ji <= nsoluble) .OR. &
             (JM >  nsoluble .AND. ji >  nsoluble) ) THEN
          IF (jm < ji ) xmov(jm,ji) = 1._dp
          idest(jm,ji) = MAX(jm,ji)
        END IF
        IF (jm >  nsoluble .AND. ji <= nsoluble) THEN
          IF (jm - nsoluble + ndiff <= ji) THEN
            xmov(jm,ji) = 1._dp
            idest(jm,ji) = ji
          ELSE
            idest(jm,ji) = jm - nsoluble + ndiff
          END IF
        END IF
        IF (jm <= nsoluble .AND. ji >  nsoluble) THEN
          IF ( (ji - nsoluble + ndiff) > jm) THEN
            xmov(jm,ji) = 1._dp
            idest(jm,ji) = ji - nsoluble + ndiff
          ELSE
            idest(jm,ji) = jm
          END IF
        ENDIF
      ENDDO
    ENDDO

!!$      print*, "mask 1:", locoagmask(1,:)
!!$      print*, "mask 2:", locoagmask(2,:)
!!$      print*, "mask 3:", locoagmask(3,:)
!!$      print*, "mask 4:", locoagmask(4,:)
!!$      print*, "mask 5:", locoagmask(5,:)
!!$      print*, "mask 6:", locoagmask(6,:)
!!$      print*, "mask 7:", locoagmask(7,:)
!!$      print*, "mov 1:", xmov(1,:)
!!$      print*, "mov 2:", xmov(2,:)
!!$      print*, "mov 3:", xmov(3,:)
!!$      print*, "mov 4:", xmov(4,:)
!!$      print*, "mov 5:", xmov(5,:)
!!$      print*, "mov 6:", xmov(6,:)
!!$      print*, "mov 7:", xmov(7,:)
!!$      print*, "dest 1:", idest(1,:)
!!$      print*, "dest 2:", idest(2,:)
!!$      print*, "dest 3:", idest(3,:)
!!$      print*, "dest 4:", idest(4,:)
!!$      print*, "dest 5:", idest(5,:)
!!$      print*, "dest 6:", idest(6,:)
!!$      print*, "dest 7:", idest(7,:)

    IF(.NOT.lacc) THEN
       RHMAX=0.95_dp
!       RHMAX=0.98_dp
    END IF

    !--- Set the lower mode boundary for the nucleation mode:
    !       (Depends on the choice of the nucleation scheme.)

    SELECT CASE (nnucl)
    CASE(0,1,2)
      ! molar mass 98.080
      ! density 1.830
       crdiv(1)=( critn * 98.080_dp / avo / &
         1.830_dp * 0.75_dp / pi )**(1._dp/3._dp)
    END SELECT

    status = 0

    IF(l_io) &
    CALL end_message(TRIM(modstr),'CORE INITIALISATION', substr)

  END SUBROUTINE gmxe_initialize_core

!==============================================================================

  SUBROUTINE gmxe_main(kproma,    klev,     dt,                     &  ! ECHAM indices
                       prhum,     ptemp,    ppress,   psat,         &  !   "   thermodynamics
                       paclc,     paopt,    pcih2o,   pclh2o,       &  !   "   aerosol-cloud properties
                       pcdnc,     picnc,    pph,      pcrain,       &  !   "   aerosol-cloud properties
                       pccn,      paerml,   paernl,   paerosol,     &  !   "   aerosol-cloud properties
                       pconv,     pvap,     pvap2,    pprod,        &  !   "   aerosol-cloud properties
                       pminc,     prho,     ppressi,  pvervel,      &
                       prdry,     prwet,    pddry )

  !
  !**** *gmxe_main* aerosol composition and dynamics interface
  !
  !  Purpose:
  !  --------
  !  This routine calculates mixed particle properties.
  !
  !  Method:
  !  -------
  !  aerosol dynamics including equilibrium composition of mixed particles
  !
  !**Interface:
  !  ----------
  !  *gmxe_main* is called from *gmxe_driver*
  !
  !   Externals
  !   ---------
  !
  !   *gmxe_averageproperties*
  !       calculates the average mass for all modes and the particle
  !       dry radius and density for the insoluble modes.
  !
  !   *gmxe_thermo*
  !       calculates the ambient radius of mixed particles
  !
  !   *gmxe_dgas*
  !       calculates the sulfate condensation on existing particles
  !
  !   *gmxe_dnum*
  !       calculates new gas phase sulfate and aerosol numbers and masses
  !       after condensation, nucleation and coagulation over one timestep
  !
  !   *gmxe_dconc*
  !       repartitions aerosol number and mass between the
  !       the modes to account for condensational growth and the formation
  !       of an accumulation mode from the upper tail of the aitken mode and
  !       of a coarse mode from the upper tail of the accumulation mode
  !
  !   *gmxe_evap*
  !        evaporation of sulfate particles in the stratosphere
  !
  !   *gmxe_rad*
  !       calculates radiative properties of aerosol particles
  !
  !--- Parameter list:
  !
  !  ppress     = atmospheric pressure at time t+1 [Pa]
  !  prhum      = atmospheric relative humidity [0-1]
  !  ptemp      = atmospheric temperature at time t+1 [K]
  !  paerml     = total aerosol mass for each compound [umol m-3 (air)]
  !  paernl     = aerosol number for each mode [N cm-3]
  !  pgasi      = conc. of gas phase species i [molec. cm-3]
  !  prwet      = aerosol wet radius  [cm]
  !  prdry      = aerosol dry radius  [cm]
  !  pddry      = aerosol dry density [g cm-3]
  !

  USE MESSY_GMXE_AERCHEM,    ONLY: aerchem_driver
  USE MESSY_GMXE_OC_AGING,   ONLY: gmxe_oc_aging
  USE MESSY_GMXE_SOA,        ONLY: nspec, soa_driver, &
                                   excl_str_soa, exclude_max_soa

  IMPLICIT NONE

  INTEGER :: kproma,klev,jk,jl,jc,jm,kcomp,jcomp,ic,it, jt

  REAL(dp):: dt, ztmst, zcondg

  REAL(dp), DIMENSION(kproma,klev,0:nmod,0:naero) :: paerosol
  REAL(dp), DIMENSION(kproma,klev,0:naertot)      :: paerml
  REAL(dp), DIMENSION(kproma,klev,0:naertot)      :: paerml1
  REAL(dp), DIMENSION(kproma,klev,0:naertot)      :: paerml2
  REAL(dp), DIMENSION(kproma,klev,ngas,nmod)      :: pgas
  REAL(dp), DIMENSION(kproma,klev,ngas)           :: pgasi
  REAL(dp), DIMENSION(kproma,klev,nmod)           :: paernl, &
                                                     prwet, prdry , &
                                                     pddry,         &
                                                     pccn,  pph
  REAL(dp), DIMENSION(kproma,klev)   ::  ppress, prho,   ptemp,  prhum, psat,  &
                                         paclc,  pcih2o, pclh2o, pcdnc, paopt, &
                                         picnc,  pcrain, pprod,  pvap2, pvap,  &
                                         zncrit,         pminc,  pso4g, pvervel
  REAL(dp)                           ::  pconv
  REAL(dp), DIMENSION(kproma,klev+1) ::  ppressi
  REAL(dp), DIMENSION(kproma,klev)   ::  zmass_post2, zmass_pre2
  REAL(dp):: pre(spec_number), post(spec_number)
  INTEGER :: jmod, idx1
  ! mz_ht_20100827+
  !For evaporation
  REAL(dp), DIMENSION(kproma,klev,0:nsoluble) :: input_h2so4,  &
                                                 output_h2so4, &
                                                 delta_h2so4
  REAL(dp), DIMENSION(kproma,klev,0:nsoluble) :: frac_h2so4,   &
                                                 frac_hso4m,   &
                                                 frac_so4mm
  ! mz_ht_20100827-
  ! additional variables for SOA
  REAL(dp), DIMENSION(kproma,klev,nmod)            :: zbulkm
  LOGICAL                                          :: npass
  REAL(dp), DIMENSION(kproma,klev,0:NSPEC, 0:nmod) :: xt_aer


#ifdef _FTRACE
CALL FTRACE_REGION_BEGIN("gmxe_main_region_7")
#endif

  !--- 0) Initialisations: -------------------------------------------------
  !
  pgasi (:,:,:)     = zero
  pgas  (:,:,:,:)   = zero
  !
  zcondg            = 1._dp
  IF(lcond_gas_phase) &
    zcondg            = ZERO
  !
  ztmst = dt
  !
#ifdef _FTRACE
CALL FTRACE_REGION_BEGIN("gmxe_main_region_1")
#endif
  !
  !--------------------------------------------------------------------
  !
  !--- 1) Calculation of particle properties under ambient conditions: -----
  !
!CDIR NOIEXPAND
  CALL gmxe_averageproperties(kproma, klev,  paernl, paerml, prwet,  prdry,    &
                              pddry,  ptemp, prhum,  psat,   ppress, prho,     &
                              paerosol)

  !
  !--------------------------------------------------------------------
  !

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_main_region_1")
#endif
  !
  ! if lcond is true the kinetic limitation for all gas phase compounds
  !  over ALL modes is calculated once at the beginning (HERE!!!!)
!!$
  paerml(1:kproma,1:klev,:) =  paerml(1:kproma,1:klev,:) / pconv

  !
  ! mz_kp_20080112+    ! Changed order so pgas is set (for any microphysics)
  jm = 0
  DO jc=1,ngas
     jcomp = td%gas(jc)%aerml_idx(jm)
     IF (jcomp == nwh2o(jm) .OR. jcomp == 0) CYCLE ! exclude H2O vapor
     pgasi(1:kproma,1:klev,jc) = paerml(1:kproma,1:klev,jcomp)
  ENDDO

  IF (LCOND) THEN
     !
!CDIR NOIEXPAND
!    print*, "Before dgas: ", paernl(1,31,1:nmod)

    CALL gmxe_dgas(kproma, klev,   pgasi,  paernl,         &
                   ptemp,  ppress, prwet,  pgas,   ztmst,  &
                   1, nmod)
!    print*, "after dgas: ", paernl(1,31,1:nmod)

     !
     ! condensable gases (by pure diffusion limited processes)
     ! pgas (1:kproma,1:klev,jc,jm)
     !

  ELSE


    DO jc=1,ngas
      DO jm=1,nmod
        jcomp = td%gas(jc)%aerml_idx(jm)
        IF (jcomp == nwh2o(0) .OR. jcomp == 0) CYCLE ! exclude H2O vapor
        pgas(1:kproma,1:klev,jc,jm) = pgasi(1:kproma,1:klev,jc) / &
          ( REAL(nmod,dp) + zcondg )
      ENDDO
!!$      pgasi(1:kproma,1:klev,jc)  = pgasi(1:kproma,1:klev,jc) / &
!!$        ( REAL(nmod,dp) + zcondg)
      DO jm=1,nmod
        pgasi(1:kproma,1:klev,jc)  = pgasi(1:kproma,1:klev,jc) - &
                                     pgas(1:kproma,1:klev,jc,jm)
      END DO
    END DO
  ENDIF ! lcond

  DO jc=1,ngas
     jm = 0
     jcomp = td%gas(jc)%aerml_idx(jm)
     IF (jcomp == nwh2o(jm) .OR. jcomp == 0) CYCLE ! exclude H2O vapor
     paerml(1:kproma,1:klev,jcomp) = pgasi(1:kproma,1:klev,jc)
  END DO

  paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) * pconv
  pgas(1:kproma,1:klev,:,1:nmod) = pgas(1:kproma,1:klev,:,1:nmod) * pconv

  !--------------------------------------------------------------------
  IF (neqm > -1) THEN
  !--------------------------------------------------------------------

#ifdef _FTRACE
CALL FTRACE_REGION_BEGIN("gmxe_main_region_2")
#endif


  !--- 2) Calculate aerosol composition, update dry and ambient radii and density
  !       gas and aerosol species [umol m-3 (air)]:
  !

!CDIR NOIEXPAND
!  print*, "Before thermo: ", paernl(1,31,1:nmod)

  CALL gmxe_thermo(kproma, klev,   prhum, paerml,  paernl,                  &
                   psat,   paclc,  paopt, pcih2o,  pclh2o, pcdnc,    picnc, &
                   pph,    pcrain, pccn,                                    &
                   prwet,  prdry ,        pddry,   ptemp,  paerosol,        &
                   pvap,   pvap2,  pprod, pminc,   ppress, prho,            &
                   ztmst,  ppressi,        pvervel, pgas )

!  print*, "after thermo: ", paernl(1,31,1:nmod)

  !
#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_main_region_2")
#endif


  CALL gmxe_averageproperties(kproma, klev,  paernl, paerml, prwet,  prdry,    &
                              pddry,  ptemp, prhum,  psat,   ppress, prho,     &
                              paerosol)
! print*, "ave prop 2:", MINVAL(prdry), MAXVAL(prdry), MINVAL(prwet), MAXVAL(prwet)
  !--------------------------------------------------------------------
  END IF ! (neqm > -1)
  !--------------------------------------------------------------------
  !

!  call sum_mass3(kproma,klev,paerml,zmass_pre2)
    !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--- All species [umol m-3 (air)] => [molecules cm-3 (air)]:
  !
  paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) / pconv
  pgas(1:kproma,1:klev,:,1:nmod) = pgas(1:kproma,1:klev,:,1:nmod) / pconv

  ! pgas is set (for any microphysics)
  DO jc=1,ngas
    jm = 0
    jcomp = td%gas(jc)%aerml_idx(jm)
    IF (jcomp == nwh2o(jm) .OR. jcomp == 0) CYCLE ! exclude H2O vapor
    pgasi(1:kproma,1:klev,jc) = paerml(1:kproma,1:klev,jcomp)
  END DO
  !
  !--------------------------------------------------------------------

  !
  !--------------------------------------------------------------------
  !!! IF(lcoag .OR. (lnucl .AND. nnucl > 0))  THEN         !! KP Origional
  IF(lcond .OR. lcoag .OR. (lnucl .AND. nnucl > 0))  THEN  !! KP Need to redistribute aersol number / mass if lcond=TRUE
  !--------------------------------------------------------------------

#ifdef _FTRACE
    CALL FTRACE_REGION_BEGIN("gmxe_main_region_4")
#endif
  !
  !--- 4) Calculate change in particle number concentrations ---------------
  !       due to nucleation and coagulation:
  !       Change particle mass/number relationships
  !

  !--- Assign H2SO4 for nucleation  !
    jcomp = td%gas_sulph_idx

    IF (lmass_diag) THEN
      paerml1 = paerml
      paerml1(1:kproma,1:klev,nh2so4(0)) = &
        paerml1(1:kproma,1:klev,nh2so4(0)) + pgasi(1:kproma,1:klev,jcomp)
      call sum_mass3(kproma,klev,paerml1,zmass_pre2)
    END IF
    pso4g(1:kproma,1:klev) = pgasi(1:kproma,1:klev,jcomp)


    !CDIR NOIEXPAND

!    print*, "before dnum aernl:" , paernl(1,31,1:nmod)

    CALL gmxe_dnum(kproma, klev,  pso4g,  paerml, paernl, ptemp,        &
        ppress, prhum, prwet,  pddry,  pgas,   zncrit, ztmst )

!    print*, "after dnum aernl:" , paernl(1,31,1:nmod)

  !! set pgasi (for h2so4) to updated h2so4 concentration
    pgasi(1:kproma,1:klev,jcomp) =  pso4g(1:kproma,1:klev)

    If (lmass_diag) THEN
      paerml2 = paerml
      paerml2(1:kproma,1:klev,nh2so4(0)) = &
        paerml2(1:kproma,1:klev,nh2so4(0)) + pgasi(1:kproma,1:klev,jcomp)

      call sum_mass3(kproma,klev,paerml2,zmass_post2)
      DO jk=1,klev
        DO jl=1,kproma
          pre(:)  = 0._dp
          post(:) = 0._dp
          DO jmod=0,nmod
            DO jc=1,spec_number
              idx1 = species(jc)%aermlidx(jmod)
              if (idx1 == 0 .OR. idx1 == species(jc)%aernlidx(jmod)) cycle
              pre(jc)  = pre(jc) + paerml1(jl,jk,idx1)
              post(jc) = post(jc) + paerml2(jl,jk,idx1)
            END DO
          END DO
          IF (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) >                &
            (10._dp * SPACING(zmass_post2(jl,jk))) ) THEN
            print*, "mass changes itemwise; dnum box: ",jl,jk," Error: (%)",  jm, &
              (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) /               &
              zmass_post2(jl,jk)*100._dp),                               &
              zmass_pre2(jl,jk), zmass_post2(jl,jk)
            DO jc=1,spec_number
              IF (ABS(pre(jc) - post(jc)) > 1.e-30_dp ) THEN
                DO jmod=0,nsoluble
                  idx1 = species(jc)%aermlidx(jmod)
                  print*, trim(species(jc)%name), "  per mode: ", jmod, jc, jl, jk, idx1, &
                    paerml1(jl,jk,idx1)*species(jc)%molmass, &
                    paerml2(jl, jk, idx1)*species(jc)%molmass, &
                    (paerml1(jl,jk,idx1)- paerml2(jl, jk, idx1))*species(jc)%molmass
                END DO
              ENDIF
            END DO
          END IF
        END DO
      END DO

    ENDIF

  END IF

  DO jc=1,ngas
     DO jm=0,nmod
        jcomp = td%gas(jc)%aerml_idx(jm)
        IF (jcomp == nwh2o(0) .OR. jcomp == 0) CYCLE ! exclude H2O vapor
        IF (jm == 0) THEN
           paerml(1:kproma,1:klev,jcomp) =  pgasi(1:kproma,1:klev,jc)
        ELSE
           paerml(1:kproma,1:klev,jcomp) =  paerml(1:kproma,1:klev,jcomp) + &
                                            pgas(1:kproma, 1:klev, jc, jm)
        END IF
     END DO
  END DO

  IF(lcond .OR. lcoag .OR. (lnucl .AND. nnucl > 0))  THEN
    paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) * pconv

    CALL gmxe_averageproperties(kproma, klev,  paernl, paerml, prwet,  prdry,  &
                                pddry,  ptemp, prhum,  psat,   ppress, prho,   &
                                paerosol)

    paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) / pconv
  ENDIF

! mz_ht_20100826+
! 4.5 ) Potential Evaporation of H2SO4
!       (mainly relevant in the upper stratosphere)

! Select H2SO4 and SO4mm concentration in all 8 modes (7 + gas)
! in molecules cm-3
  ! for Evaporation in the stratosphere and mesosphere

  input_h2so4(:,:,:)  = 0.0_dp
  output_h2so4(:,:,:) = 0.0_dp
  frac_h2so4(:,:,:)   = 0.0_dp
  frac_hso4m(:,:,:)   = 0.0_dp
  frac_so4mm(:,:,:)   = 0.0_dp
  delta_h2so4(:,:,:)  = 0.0_dp

  DO jm=0,nsoluble
    kcomp = nh2so4(jm)
    IF (kcomp /= 0) &  !Only select species that are "switched on"
    input_h2so4(1:kproma,1:k300,jm) = input_h2so4(1:kproma,1:k300,jm) &
                                    + paerml(1:kproma,1:k300,kcomp)
    kcomp = nhso4m(jm)
    IF (kcomp /= 0) &  !Only select species that are "switched on"
    input_h2so4(1:kproma,1:k300,jm) = input_h2so4(1:kproma,1:k300,jm) &
                                    + paerml(1:kproma,1:k300,kcomp)
    kcomp = nso42m(jm)
    IF (kcomp /= 0) &  !Only select species that are "switched on"
    input_h2so4(1:kproma,1:k300,jm) = input_h2so4(1:kproma,1:k300,jm) &
                                    + paerml(1:kproma,1:k300,kcomp)

    output_h2so4(1:kproma,1:k300,jm) = input_h2so4(1:kproma,1:k300,jm)

    kcomp = nh2so4(jm)
    IF (kcomp /= 0) THEN  !Only select species that are "switched on"
      DO jk=1,k300
        DO jl=1,kproma
          IF (input_h2so4(jl,jk,jm) > zeps)             &
            frac_h2so4(jl,jk,jm) = paerml(jl,jk,kcomp)  &
                                 / input_h2so4(jl,jk,jm)
        END DO
      END DO
    END IF

    kcomp = nhso4m(jm)
    IF (kcomp /= 0) THEN
      DO jk=1,k300
        DO jl=1,kproma
          IF (input_h2so4(jl,jk,jm) > zeps)             &
            frac_hso4m(jl,jk,jm) = paerml(jl,jk,kcomp)  &
                                 / input_h2so4(jl,jk,jm)
        END DO
      END DO
    END IF

    kcomp = nso42m(jm)
    IF (kcomp /= 0) THEN
      DO jk=1,k300
        DO jl=1,kproma
          IF (input_h2so4(jl,jk,jm) > zeps)             &
            frac_so4mm(jl,jk,jm) = paerml(jl,jk,kcomp)  &
                                 / input_h2so4(jl,jk,jm)
        END DO
      END DO
    END IF

  END DO         ! nsoluble

  CALL gmXe_evap2(kproma,k300,input_h2so4(1:kproma,1:k300,0:nsoluble),  &
                              output_h2so4(1:kproma,1:k300,0:nsoluble), &
                              prwet(1:kproma,1:k300,1:nsoluble),        &
                              ptemp(1:kproma,1:k300),                   &
                              ppress(1:kproma,1:k300),                  &
                              paernl(1:kproma,1:k300,1:nsoluble), ztmst )


  DO jm=0,nsoluble
    delta_h2so4(1:kproma,1:k300,jm) = output_h2so4(1:kproma,1:k300,jm)  &
                                    -  input_h2so4(1:kproma,1:k300,jm)
    kcomp = nh2so4(jm)
    IF (kcomp /= 0) &
    paerml(1:kproma,1:k300,kcomp) = paerml(1:kproma,1:k300,kcomp)   &
                                  + delta_h2so4(1:kproma,1:k300,jm) &
                                  * frac_h2so4(1:kproma,1:k300,jm)
    kcomp = nhso4m(jm)
    IF (kcomp /= 0) &
    paerml(1:kproma,1:k300,kcomp) = paerml(1:kproma,1:k300,kcomp)   &
                                  + delta_h2so4(1:kproma,1:k300,jm) &
                                  * frac_hso4m(1:kproma,1:k300,jm)
    kcomp = nso42m(jm)
    IF (kcomp /= 0) &
    paerml(1:kproma,1:k300,kcomp) = paerml(1:kproma,1:k300,kcomp)   &
                                  + delta_h2so4(1:kproma,1:k300,jm) &
                                  * frac_so4mm(1:kproma,1:k300,jm)
    kcomp = nhp(jm)
    IF (kcomp /= 0) &
    paerml(1:kproma,1:k300,kcomp) = paerml(1:kproma,1:k300,kcomp)   &
                                  + delta_h2so4(1:kproma,1:k300,jm) &
                                  * (frac_so4mm(1:kproma,1:k300,jm) &
                                  * 2._dp                           &
                                  + frac_hso4m(1:kproma,1:k300,jm) )
  END DO   ! nsoluble
! mz_ht_20100826-


  !--- 5) Call to non-equilibrium aerosol chemistry
  !       in sub-submodel aerchem
  IF (L_aerchem) THEN
    DO jm=1,nsoluble
      IF (.NOT. AERCHEM(JM)) CYCLE
    ! paerml is already in molecules/cm^3
      kcomp = species(spec_idx_h2o)%aermlidx(jm)
      CALL AERCHEM_DRIVER(kproma, klev, jm, kcomp, naertot,          &
                          paerml, prwet(1:kproma,1:klev,jm),         &
                          ptemp, ppress, ztmst)                    !    &
    ENDDO
    paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) * pconv

    CALL gmxe_averageproperties(kproma, klev,  paernl, paerml, prwet,  prdry,  &
                                pddry,  ptemp, prhum,  psat,   ppress, prho,   &
                                paerosol)

    paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) / pconv
  ENDIF


!!$  paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) * pconv
!!$  call sum_mass3(kproma,klev,paerml,zmass_pre2)
!!$  paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) * pconv

  IF (L_oc_aging) THEN
    CALL GMXE_OC_AGING(kproma, klev, naertot, nmod, nsoluble, ztmst, &
                       paerml, paerosol(1:kproma,1:klev,1:nmod,ioc_kappa) )

    paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) * pconv

    CALL gmxe_averageproperties(kproma, klev,  paernl, paerml, prwet,  prdry,  &
                                pddry,  ptemp, prhum,  psat,   ppress, prho,   &
                                paerosol)

    paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) / pconv
  ENDIF


  IF (L_SOA) THEN
    ! Determine mass of seed aerosol, on which the SOA condenses, i.e. the
    !   wet aerosol mass
    zbulkm(:,:,:) = 0._dp
    DO jm=1,nmod
       IF (.NOT. LSOA(jm)) CYCLE
       DO jt=1,spec_number
          DO jk=1,EXCLUDE_MAX_SOA
             IF (TRIM(species(jt)%name) == TRIM(EXCL_STR_SOA(jk))) &
                  kcomp = 0
          END DO
          npass = species(jt)%npassive(jm) ! mz_dk_20120120
          IF(kcomp == 0 .or. kcomp == species(jt)%aernlidx(jm) .or. npass) CYCLE
          DO jk=1,klev
             DO jl=1,kproma
                            ! zbulkm is in [?g  m-3(air)]
                zbulkm (jl,jk,jm) = zbulkm (jl,jk,jm) +        &
                                    paerml(jl,jk,kcomp) * pconv * &
                                    species(jt)%molmass
             END DO
          ENDDO
       END DO
    END DO
!!$    DO jm=1,nmod
!!$       print*, "zbulkm before SOA", jm, MAXVAL(zbulkm(:,:,jm))
!!$    END DO
    ! pack soa species into a specific array and call SOA subroutine
    xt_aer(:,:,:,:) = 0._dp
    DO jm=0,nmod
       DO jt=1,spec_number
          kcomp = species(jt)%aermlidx(jm)
          jcomp = species(jt)%SOAidx(jm)
          IF (jcomp == 0) CYCLE
!!$          print*, "species in time loop in", jt, jcomp, kcomp, species(jt)%name, jm, &
!!$               MAXVAL(paerml(:,:,kcomp))
          IF (kcomp == 0) CYCLE
          DO jk=1,klev
             DO jl=1,kproma
                xt_aer(jl,jk,jcomp,jm) = paerml(jl,jk,kcomp)
             END DO
          ENDDO
       END DO
    ENDDO
!!$    print*, "before SOA", MAXVAL(xt_aer(:,:,:,0)), &
!!$         MAXVAL(xt_aer(:,:,:,1)),MAXVAL(xt_aer(:,:,:,2)), MAXVAL(xt_aer(:,:,:,3)), &
!!$         MAXVAL(xt_aer(:,:,:,4)),MAXVAL(xt_aer(:,:,:,5)), MAXVAL(xt_aer(:,:,:,6))
    DO jk=1,klev
       DO jl=1,kproma
          CALL SOA_DRIVER(ptemp(jl,jk), ztmst, nmod, &
               xt_aer(jl,jk,0:nspec,0:nmod),         &
               zbulkm(jl,jk,1:nmod))
       END DO
    END DO
!!$    print*, "after SOA", MAXVAL(xt_aer(:,:,:,0)), &
!!$         MAXVAL(xt_aer(:,:,:,1)),MAXVAL(xt_aer(:,:,:,2)), MAXVAL(xt_aer(:,:,:,3)), &
!!$         MAXVAL(xt_aer(:,:,:,4)),MAXVAL(xt_aer(:,:,:,5)), MAXVAL(xt_aer(:,:,:,6))
    DO jm=0,nmod
       DO jt=1,spec_number
          kcomp = species(jt)%aermlidx(jm)
          jcomp = species(jt)%SOAidx(jm)
          IF (jcomp == 0) CYCLE
          IF (kcomp == 0) CYCLE
          DO jk=1,klev
             DO jl=1,kproma
                paerml(jl,jk,kcomp) = xt_aer(jl,jk,jcomp,jm)
             END DO
          ENDDO
!!$          print*, "species in time loop out", jt, jcomp, kcomp, species(jt)%name, jm, &
!!$               MAXVAL(paerml(:,:,kcomp))
       END DO
    ENDDO
    ! water uptake on SOA is calculated in the thermo subroutine
    paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) * pconv

    CALL gmxe_averageproperties(kproma, klev,  paernl, paerml, prwet,  prdry,  &
                                pddry,  ptemp, prhum,  psat,   ppress, prho,   &
                                paerosol)

    paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) / pconv
  ENDIF

!!$  paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) * pconv
!!$  call sum_mass3(kproma,klev,paerml,zmass_pre2)
!!$  paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) * pconv
!!$
!!$  DO jk=1,klev
!!$    DO jl=1,kproma
!!$      IF (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) >                &
!!$        (3._dp * SPACING(zmass_post2(jl,jk))) ) THEN
!!$        print*, "mass changes itemwise; box: ",jl,jk," Error: (%)", &
!!$          (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) /               &
!!$          zmass_post2(jl,jk)*100._dp),                               &
!!$          zmass_pre2(jl,jk), zmass_post2(jl,jk)
!!$      END IF
!!$    END DO
!!$  END DO
!!$

  !
  !--- 6) Repartitition particles among the modes:
  !
  !CDIR NOIEXPAND
  IF(.NOT. lsize)THEN
!! mz_kp_20080112 lsize reshapes size distribution => double couting processes.
 ! redistributing after ANY microphysical process (in GMXE_MAIN or EQSAM)
!    print*, "Before dconc: ", paernl(1,31,1:nmod)

    IF (lmass_diag) THEN
      paerml1 = paerml
      call sum_mass3(kproma,klev,paerml1,zmass_pre2)
    ENDIF

    IF(lcoag .OR. lcond .OR. lcoat .OR. lnucl .OR. l_aerchem .OR. neqm > 0) &
    CALL gmxe_dconc(kproma, klev, paerml, paernl, prdry)

!    print*, "after dconc: ", paernl(1,31,1:nmod)
    IF (lmass_diag) THEN
      paerml2 = paerml
      call sum_mass3(kproma,klev,paerml2,zmass_post2)
      DO jk=1,klev
        DO jl=1,kproma
          pre(:)  = 0._dp
          post(:) = 0._dp
          DO jmod=0,nmod
            DO jc=1,spec_number
              idx1 = species(jc)%aermlidx(jmod)
              if (idx1 == 0 .OR. idx1 == species(jc)%aernlidx(jmod)) cycle
              pre(jc)  = pre(jc) + paerml1(jl,jk,idx1)
              post(jc) = post(jc) + paerml2(jl,jk,idx1)
            END DO
          END DO
          IF (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) >                &
            (10._dp * SPACING(zmass_post2(jl,jk))) ) THEN
            print*, "mass changes itemwise; dconc box: ",jl,jk," Error: (%)",  jm, &
              (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) /               &
              zmass_post2(jl,jk)*100._dp),                               &
              zmass_pre2(jl,jk), zmass_post2(jl,jk)
            DO jc=1,spec_number
              IF (ABS(pre(jc) - post(jc)) > 1.e-30_dp ) THEN
                DO jmod=0,nmod
                  idx1 = species(jc)%aermlidx(jmod)
                  print*, trim(species(jc)%name), "  per mode: ", jmod, jc, jl, jk, idx1, &
                    paerml1(jl,jk,idx1)*species(jc)%molmass, &
                    paerml2(jl, jk, idx1)*species(jc)%molmass, &
                    (paerml1(jl,jk,idx1)- paerml2(jl, jk, idx1))*species(jc)%molmass
                END DO
              ENDIF
            END DO
          END IF
        END DO
      END DO
    END IF

  ENDIF ! .NOT. lsize

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_main_region_4")
#endif


  !--- 7) Assign output values:

#ifdef _FTRACE
CALL FTRACE_REGION_BEGIN("gmxe_main_region_5")
#endif
  !
  !--- Trace gases + aerosols [molecules cm-3 (air)] => [umol m-3 (air)]:
  !
    paerml(1:kproma,1:klev,:) = paerml(1:kproma,1:klev,:) * pconv

!!$  call sum_mass3(kproma,klev,paerml,zmass_pre2)
!!$    DO jk=1,klev
!!$       DO jl=1,kproma
!!$         IF (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) >                &
!!$           (3._dp * SPACING(zmass_post2(jl,jk))) ) THEN
!!$           print*, "mass changes itemwise; box: ",jl,jk," Error: (%)", &
!!$             (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) /               &
!!$             zmass_post2(jl,jk)*100._dp),                               &
!!$             zmass_pre2(jl,jk), zmass_post2(jl,jk)
!!$         END IF
!!$       END DO
!!$     END DO#

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_main_region_5")
#endif


  IF(lcond .OR. lcoag .OR. lnucl .OR. lcoat)  THEN
  !--- 8) Re - Calculation of particle properties under ambient conditions: ---
  !            since they have been modified by the redistribution
!CDIR NOIEXPAND
    CALL gmxe_averageproperties(kproma, klev,  paernl, paerml, prwet,  prdry,  &
                                pddry,  ptemp, prhum,  psat,   ppress, prho,   &
                                paerosol)
!    print*, "ave prop 4:", MINVAL(prdry), MAXVAL(prdry), MINVAL(prwet), MAXVAL(prwet)
  !
  !--------------------------------------------------------------------
  ENDIF

  DO jm=1,nmod
   kcomp = nwh2o(jm) ! H2O aerosol
   DO jk = 1,klev
     DO jl = 1,kproma
        ! air temperature   [K]
        paerosol(jl,jk,jm,iTT)        = ptemp     (jl,jk)
        ! relative humidity [%]
        paerosol(jl,jk,jm,iRH)        = prhum     (jl,jk)
        paerosol(jl,jk,jm,iRH)        = paerosol  (jl,jk,jm,iRH) * FACp2
        paerosol(jl,jk,jm,iWH2O)      = paerml    (jl,jk,kcomp) * mwh2o

 ! hygroscopic growth factor [-]
          paerosol(jl,jk,jm,iGF)        = 1._dp
        IF (prdry(jl,jk,jm) > REALZERO) &
          paerosol(jl,jk,jm,iGF)        = prwet(jl,jk,jm) / prdry(jl,jk,jm)
     END DO ! jl
   END DO ! jk
  END DO ! jm


#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_main_region_7")
#endif


! mz_ht_20090310+
! subroutine checking for negative H+/OH- balance
  CALL Hp_OHm_check(kproma, klev, paerml,paerosol,1)
! mz_ht_20090310-

!  print*, "finalise gmxe_main"

END SUBROUTINE gmxe_main
! -------------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! PRIVATE ROUTINES
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


!------------------------------------------------------------------------------!

SUBROUTINE gmxe_averageproperties(kproma, klev,  paernl, paerml, prwet,  prdry,&
                                  pddry,  ptemp, prhum,  psat,   ppress, prho, &
                                  paerosol)
  !
  !  Author:
  !  --------
  !  E. Vignati, JRC/EI (original source)                10/2000
  !  P. Stier, MPI      (f90-version, changes, comments)    2001
  !
  !  Modifications:
  !  --------------
  !  H. Tost, UMZ generalised species structure
  !  S. Metzger, MPI-CHEM (modified for use with EQSAM3),   July 2008
  !
  !  Purpose:
  !  ---------
  !  Calculation of the (dry) radius and the density
  !  of the particles of the insoluble modes.
  !
  !  Interface:
  !  ----------
  !  gmxe_averageproperties is called from gmxe
  !
  !  Externals:
  !  ----------
  !  none
  !
  !--- Parameter list:
  !
  ! paerml(kproma,klev,0:naertot) = total aerosol mass for each
  !                                 compound [umol m-3]
  ! paernl(kproma,klev,nmod)      = aerosol number for each mode [N cm-3]
  ! prwet (kproma,klev,nmod)      = aerosol wet radius  [cm]
  ! prdry (kproma,klev,nmod)      = aerosol dry radius  [cm]
  ! pddry (kproma,klev,nmod)      = aerosol dry density [g cm-3]
  !
  !   !
  !--- Local variables:
  !
  ! zvbulk                        = average volume for single particle in
  !                                 each mode [cm3]
  ! zmbulk                        = average   mass for single particle in
  !                                 each mode [g]
  !
  ! kproma, klev, it, jc, kcomp   = loop indices
  !

  USE MESSY_MAIN_TOOLS,          ONLY: tlucuaw, jptlucu1, jptlucu2

  !--- Parameters:

  IMPLICIT NONE

  INTEGER :: kproma, klev, it

  REAL(dp), DIMENSION(kproma,klev,0:nmod,0:naero) :: paerosol
  REAL(dp), DIMENSION(kproma,klev,0:naertot)      :: paerml
  REAL(dp), DIMENSION(kproma,klev,nmod)           :: paernl, prwet, prdry , pddry,  &
                                                     zdrybulk, zdrybulkm, zdrybulkv,&
                                                     zwatbulkv, zdrybulk2

  REAL(dp), DIMENSION(kproma,klev)   :: ppress, prho, ptemp,  prhum, psat

  REAL(dp), DIMENSION(kproma,klev)   ::  &
             znumber,      zhelp,        zm2n,         zcmm,      &
             znu,          zmol,         zmbulk,       zvbulk,    zwatbulk
  REAL(dp), DIMENSION(kproma,klev)   ::  ztotwat

  !--- Local variables:

  INTEGER :: jm, jk, jl, jc, kcomp, jcomp, lcomp, jt, min_kcomp
  LOGICAL :: npass ! mz_dk_20120120

#ifdef _FTRACE
CALL FTRACE_REGION_BEGIN("gmxe_averageproperties_region_1")
#endif

  jcomp = nwh2o(0) ! H2O vapor
  !--------------------------------------------------------------------

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_averageproperties_region_1")
CALL FTRACE_REGION_BEGIN("gmxe_averageproperties_region_2")
#endif

  !--- Calculate bulk aerosol mass and volume:

  zdrybulk  (:,:,:) = ZERO ![mol m-3 (air)]
  zdrybulk2 (:,:,:) = ZERO ![mol m-3 (air)]
  zdrybulkm (:,:,:) = ZERO ![ g  m-3 (air)]
  zdrybulkv (:,:,:) = ZERO ![cm3 m-3 (air)]
  zwatbulkv (:,:,:) = ZERO ![cm3 m-3 (air)]

  DO jm=1,nmod
    DO jt=1,spec_number
      kcomp = species(jt)%aermlidx(jm)
      npass = species(jt)%npassive(jm) ! mz_dk_20120120
      IF(kcomp == 0 .or. kcomp == species(jt)%aernlidx(jm) .or. npass) CYCLE ! mz_dk_20120124, added npass
      IF(kcomp /= nwh2o(jm)) THEN ! exclude H2O
        DO jk=1,klev
          DO jl=1,kproma
            ! [mol m-3(air)]
            zdrybulk (jl,jk,jm) = zdrybulk (jl,jk,jm) +        &
                                   paerml(jl,jk,kcomp) * FACm6
            ! [g   m-3(air)]
            zdrybulkm (jl,jk,jm) = zdrybulkm (jl,jk,jm) +        &
                                   paerml(jl,jk,kcomp) * FACm6 * &
                                   species(jt)%molmass
            ! [cm3 m-3(air)]
            zdrybulkv (jl,jk,jm) = zdrybulkv (jl,jk,jm) +        &
                                   paerml(jl,jk,kcomp) * FACm6 * &
                                   species(jt)%molmass /         &
                                   species(jt)%density
          END DO
        END DO
      ELSE ! H2O
        DO jk=1,klev
          DO jl=1,kproma
            ! [cm3 m-3(air)]
            zwatbulkv (jl,jk,jm) = zwatbulkv (jl,jk,jm) +        &
                                   paerml(jl,jk,kcomp) * FACm6 * &
                                   species(jt)%molmass /         &
                                   species(jt)%density
          ENDDO
        END DO
      END IF
    END DO
  END DO

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_averageproperties_region_2")
#endif

  !
  !--- Calculate dry and ambient count median radii and density:
  !

  DO jm=1,nmod

#ifdef _FTRACE
CALL FTRACE_REGION_BEGIN("gmxe_averageproperties_region_3")
#endif

     DO jk=1,klev
        DO jl=1,kproma
           IF (paernl(jl,jk,jm)  > cmin_aernl .AND. &
               zdrybulkv(jl,jk,jm) > cmin_avol) THEN

               znumber(jl,jk)    = paernl(jl,jk,jm) * FACp6

!If mass but no number, calculate a new number assuming a const. radius
           ELSE IF (paernl(jl,jk,jm) <= cmin_aernl .AND. &
                    zdrybulkv(jl,jk,jm) > cmin_avol) THEN

               zm2n   (jl,jk)    = FACp0 / (4._dp/3._dp * pi *               &
                                   zdrybulkm(jl,jk,jm)/zdrybulkv(jl,jk,jm) * &
                                   FACp3*crdiv_r3(jm))
               zhelp  (jl,jk)    = zdrybulkm(jl,jk,jm) * FACm3 * zm2n(jl,jk)
               znumber(jl,jk)    = zhelp    (jl,jk)
               paernl (jl,jk,jm) = znumber(jl,jk) * FACp6

!If number but no mass, prescribe radius
           ELSE IF (paernl(jl,jk,jm) > cmin_aernl .AND. &
                    zdrybulkv(jl,jk,jm) <= cmin_avol) THEN

               znumber(jl,jk)    = ZERO

           ELSE

               znumber(jl,jk)    = ZERO

           END IF
        END DO ! jl
     END DO ! jk

     DO jk=1,klev
        DO jl=1,kproma

           IF (znumber(jl,jk) * FACm6 > cmin_aernl .AND. &
               zdrybulkv(jl,jk,jm) > cmin_avol) THEN ! [N m-3]

             !zmbulk = mass per particle = total mass / total number
             zmbulk(jl,jk)=zdrybulkm(jl,jk,jm)/znumber(jl,jk) ! [g]

             zvbulk(jl,jk)=zdrybulkv(jl,jk,jm)/znumber(jl,jk) ! [cm3]
             !--- dry radius:  ! [cm]
             prdry (jl,jk,jm)=(0.75_dp/pi*zvbulk(jl,jk))**(1._dp/3._dp) * &
                              ram2cmr(jm)
             !--- dry density:! [g cm-3]
             pddry(jl,jk,jm)=zmbulk(jl,jk)/zvbulk(jl,jk)
             ! now a wet volume : ! [cm3]
             zvbulk(jl,jk) = ( zdrybulkv(jl,jk,jm) + zwatbulkv(jl,jk,jm) ) / &
                             znumber(jl,jk)
             !--- wet radius: ! [cm]
             prwet(jl,jk,jm)=(0.75_dp/pi*zvbulk(jl,jk))**(1._dp/3._dp) * &
                             ram2cmr(jm)

           ELSE

             prdry  (jl,jk,jm) = 0._dp
             prwet  (jl,jk,jm) = 0._dp
             pddry  (jl,jk,jm) = 0._dp

           END IF

          !! mz_kp_20080112 Error checks for radius going out of bounds
!!$          IF(jm == KS)THEN
!!$            IF(prdry(jl,jk,jm) > 150.0E-7_dp )THEN
!!$               print*,' Error Aitken Mode radius out of reasonable bounds!!! '
!!$               print*,jl,jk,jm,'Rad=',prdry(jl,jk,jm)*1.0E7,' zvbulk=',zvbulk(jl,jk),ram2cmr(jm)
!!$               print*,' number=',znumber(jl,jk)*FACm6,' tot vol=',zdrybulkm(jl,jk,jm)
!!$            ENDIF
!!$          ENDIF
!!$          IF(jm == AS)THEN
!!$            IF(prdry(jl,jk,jm) > 800.0E-7_dp )THEN
!!$               print*,' Error Accumulation Mode radius out of reasonable bounds!!! '
!!$               print*,jl,jk,jm,'Rad=',prdry(jl,jk,jm)*1.0E7,' zvbulk=',zvbulk(jl,jk),ram2cmr(jm)
!!$               print*,' number=',znumber(jl,jk)*FACm6,' tot vol=',zdrybulkm(jl,jk,jm)
!!$            ENDIF
!!$          ENDIF

        END DO ! jl
     END DO ! jk

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_averageproperties_region_3")
CALL FTRACE_REGION_BEGIN("gmxe_averageproperties_region_4")
#endif

  !--- Calculate aerosol number:

  IF(lnumber) THEN
  DO jk=1,klev
  DO jl=1,kproma
     zmbulk(jl,jk)=zero
     IF(zdrybulkv(jl,jk,jm) > cmin_avol) &
     zmbulk(jl,jk)=zdrybulkm(jl,jk,jm)/znumber(jl,jk) ! [g]
     !
     ! [kg]      =                 [g cm-3] -> [kg m-3]   *   [m3]
     zcmm(jl,jk) =  4._dp/3._dp*pi*pddry(jl,jk,jm)*FACp3*crdiv_r3(jm)

     ! [kg-1]    = 1.   /  [kg]
     zm2n(jl,jk) = 1._dp/zcmm(jl,jk)
     !
     !--- Convert mass to number for given number median radii and standard deviation
     !         (implicitly by the conversion factor cmr2ram) of the modes.
     !
     !    M2N  = 1/M0 = 1/(4/3 * pi * dens * R0(averageMass)**3)
     !    N(t) = M(t) * M2N * N(t-1)
     !
     ! [N m-3]      = [g]           =>   [kg]  * [kg-1]   * [N m-3(air)]
     zhelp  (jl,jk) = zmbulk(jl,jk) * FACm3 * zm2n(jl,jk) * znumber(jl,jk)
     !
     znumber(jl,jk) = zhelp (jl,jk)
     !
     ! [N cm-3]
  END DO
  END DO
  END IF ! lnumber

  DO jk=1,klev
  DO jl=1,kproma
     paernl (jl,jk,jm) = znumber(jl,jk) * FACm6
  END DO
  END DO

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_averageproperties_region_4")
#endif

  END DO ! jm

  !--------------------------------------------------------------------

 DO jm=1,nmod
   kcomp = nwh2o(jm) ! H2O aerosol
   DO jk = 1,klev
     DO jl = 1,kproma

       ! density (bulk mass) [g/m3]
       IF(zdrybulkv (jl,jk,jm) > ZERO) &
         paerosol(jl,jk,jm,iRHO)     = zdrybulkm (jl,jk,jm) / &
                                       zdrybulkv (jl,jk,jm)
       ! bulk volume         [ucm3 m-3 (air)]
       paerosol(jl,jk,jm,iVOL)       = zdrybulkv (jl,jk,jm) * FACp6
       ! total liquid matter [umol m-3 (air)]
       paerosol(jl,jk,jm,iaPM)       = zdrybulk  (jl,jk,jm) * FACp6
       ! total PM (aq & s)   [ug m-3 (air)]
       paerosol(jl,jk,jm,iPMt)       = zdrybulkm (jl,jk,jm) * FACp6
       ! aerosol Water (aq)  [ug m-3 (air)]
       paerosol(jl,jk,jm,iWH2O)      = paerml    (jl,jk,kcomp) * mwh2o

       ! hygroscopic growth factor [-]
         paerosol(jl,jk,jm,iGF)  = 1._dp
        IF(prdry(jl,jk,jm) > REALZERO) &
         paerosol(jl,jk,jm,iGF)  = prwet(jl,jk,jm) / prdry(jl,jk,jm)
     END DO ! jl
   END DO ! jk
 END DO ! jm

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_averageproperties_region_5")
#endif

  !--------------------------------------------------------------------

END SUBROUTINE gmxe_averageproperties
!------------------------------------------------------------------------------!
SUBROUTINE gmxe_thermo(kproma, klev,   prhum,  paerml, paernl,                  &
                       psat,   paclc,  paopt,  pcih2o, pclh2o,   pcdnc, picnc,  &
                       pph,    pcrain, pccn,                                    &
                       prwet,  prdry , pddry,  ptemp,  paerosol,                &
                       pvap,   pvap2,  pprod,  pminc,  ppress,   prho,          &
                       ztmst,  ppressi,pvervel, pgas  )
  !
  !**** *gmxe_thermo* aerosol composition interface
  !
  !  Purpose:
  !  --------
  !  This routine calculates mixed particle properties.
  !
  !  Method:
  !  -------
  !  equilibrium composition of mixed particles
  !
  !**Interface:
  !  ----------
  !  *gmxe_thermo* is called from *gmxe_main*
  !

  USE MESSY_MAIN_TOOLS,      ONLY: tlucuaw, jptlucu1, jptlucu2
  USE MESSY_GMXE_OC_AGING,   ONLY: num_wsoc, spec_idx_wsoc
  USE messy_gmxe_kappa,      ONLY: calc_kappa, calc_wateruptake_bulk, salt, nss
  USE messy_gmxe_isorropia2, ONLY: isorropia_interface
  USE messy_gmxe_eqsam4clim, ONLY: eqsam4clim_interface

  IMPLICIT NONE
  INTRINSIC :: SQRT ! op_pj_20131105

  CHARACTER(LEN=*), PARAMETER :: substr = 'gmxe_thermo'

  INTEGER :: Is, it, ic,jc, jl, jk, jm, km, imod, jmod, jcomp, kcomp
  INTEGER :: kproma, klev
  INTEGER, DIMENSION(kproma,klev,nmod) :: mask
  INTEGER, DIMENSION(kproma,klev)      :: mask2, ncompi
  REAL(dp),DIMENSION(kproma,klev,nmod) :: zmask

  LOGICAL,  SAVE              :: entered = .FALSE.
  LOGICAL                     :: leqm, lcoati

  REAL(dp), DIMENSION(kproma,klev,0:nmod,0:naero) :: paerosol
  REAL(dp), DIMENSION(kproma,klev,0:naertot)      :: paerml

  REAL(dp), DIMENSION(kproma,klev,0:naertot)      :: paerml1
  REAL(dp), DIMENSION(kproma,klev,0:naertot)      :: paerml2
  REAL(dp), DIMENSION(kproma,klev,0:naertot)      :: paerml3
  REAL(dp):: zmass_post2(kproma,klev), zmass_pre2(kproma,klev)


  REAL(dp), DIMENSION(kproma,klev)  ::  &
              ptemp,  prhum,  paclc, paopt, pcih2o, pclh2o,         &
              znumber,prho,   ppress,pcdnc, picnc,                  &
              pvap,   pvap2,  pprod, pminc, psat,   pcrain !,pcsnow

  REAL(dp), DIMENSION(kproma,klev)  ::  &
              zcoat, zsurf, zsurfm,       zmax,         zcoati,      &
              zdrymass,     zwetmass,     zdryden,      zwatmass,    &
              zdryvol,      zwetvol,      zdryrad,      zwetrad,     &
              zwetden,      zdryvol_mean, zwetvol_mean, zionbulk

 REAL(dp), DIMENSION(kproma,klev) ::    &
               zfac,         zhelp,        zmbulk,      zmbulki,     &
               zice,         zm2n,         zcmm,        zmcoat,      &
               zdrybulk,     zwatbulk,     zdrybulkm,   zwatbulkm,   &
               zplus,        zmin,         zkw,         ztot, zcoef, &
               cmin_fog,     cmin_water,   cmin_cloud,  cmin_rain,   &
               znu,          zmol,         ztotwat

 REAL(dp), DIMENSION(kproma,klev+1)           :: ppressi
 REAL(dp), DIMENSION(kproma,klev)           :: pvervel, psphum, TT, RH
 REAL(dp), DIMENSION(kproma,klev,ngas)      :: pgasi
 REAL(dp), DIMENSION(kproma,klev,ngas,nmod) :: pgas

 REAL(dp), DIMENSION(kproma,klev,nmod,1,ngas)     :: xgases
 REAL(dp), DIMENSION(kproma,klev,nmod,2,nanions)  :: xanions
 REAL(dp), DIMENSION(kproma,klev,nmod,2,ncations) :: xcations
 REAL(dp), DIMENSION(kproma,klev,nmod,2,nsolutes) :: xsolutes

 REAL(dp), DIMENSION(ngas)                  :: ygases
 REAL(dp) ::   ztmst
 REAL(dp), DIMENSION(kproma,klev,nmod) ::  &
               prwet, prdry , pddry, paernl, pccn, phplus, pph
 REAL(dp)  :: xhelp1, znu_oc

 INTEGER   :: idx1, idx2, nhy
 REAL(dp)  :: zcharge

 REAL(dp), DIMENSION(kproma,klev)           :: frac, RH_rest
 REAL(dp), DIMENSION(kproma,klev,nmod)      :: WH2O, ZHP, ZDIFF
 REAL(dp)  :: V_s, V_w, M_wat
 REAL(dp), DIMENSION(kproma,klev,8)         :: diagout
 REAL(dp), DIMENSION(9)                     :: xdiag


 REAL(dp):: pre(spec_number), post(spec_number)
  !--------------------------------------------------------------------
  IF(l_io .AND. .NOT. entered) &
  CALL start_message(modstr,'calculate aerosol composition',substr)

  !--------------------------------------------------------------------
  !--- 1. Assign input values per mode:
  !--------------------------------------------------------------------
  ! NOTE: It is implicit in this treatment that gas/aerosol partitioning
  !       starts from small to large particles (using the residual gases),
  !       since smaller particles equilibrate faster, because of the larger
  !       surface to volume ratio.
  !       It is further assumed that the coating of primary/insoluble
  !       particles is preferred compared to secondary/soluble particles,
  !       since the vapor pressures over insoluble/dry surfaces is zero
  !       and hence lower than the corresponding vapor pressures over aqueous
  !       solutions (that can be present only in the soluble modes).
  !
  !       We further assume that primary organic, black carbon (OC/BC) and
  !       dust particles are (partly) transferred into the soluble modes
  !       when coated, with a chemical speciation of the transferred mass
  !
  !       The coating/ageing itself depends on the condensation of selected
  !       aerosol precursor gases (those in messy_gmxe_parameters.inc with
  !       '.T.' in the first column 'GP'), including water vapor.
  !       The condensed mass of gases is calculated within gmxe_dgas, and
  !       the actual mass transfer (from insoluble to soluble modes) depends
  !       on the total surface area of the primary particles and the number
  !       of moles that are required for a 1:1 stoichiometric surface reaction.
  !--------------------------------------------------------------------
  !
  zmask   (:,:,:)     = ZERO
  mask    (:,:,:)     = 0
  mask2   (:,:)       = 0
  xcations(:,:,:,:,:) = ZERO
  xanions (:,:,:,:,:) = ZERO
  xsolutes(:,:,:,:,:) = ZERO
  xgases  (:,:,:,:,:) = ZERO
  ygases (:) = 0._dp

  LLEQM(:) = .FALSE.

  DO jc=1,nss
     salt(jc)%mass(:,:,:) = 0._dp
  END DO

  !--------------------------------------------------------------------

#ifdef _FTRACE
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_1")
#endif

     !--------------------------------------------------------------------
  DO imod = 1, nmod

 !--------------------------------------------------------------------

#ifdef _FTRACE
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_2")
#endif

     !--------------------------------------------------------------------

     ! initialize switch and help array for coating
     lcoati        = .FALSE.
     zwatbulk(:,:) = ZERO

         ! Start with insoluble modes (vapor pressure is zero)
         jm  = NMODES(imod)
     IF (jm >= nlowermode .AND. jm <= nuppermode) THEN
         leqm = .TRUE.
     ELSE
         leqm = .FALSE.
     END IF
     IF (NEQM <= 0) LEQM = .FALSE.
     paerml1 = paerml
     DO jc = 1,ngas       ! ngas = nacids + 1 (ammonia)
        jcomp = td%gas(jc)%aerml_idx(0)
        IF (jcomp == nwh2o(0) .OR. jcomp == 0) CYCLE ! exclude H2O vapor
        paerml1(1:kproma,1:klev,jcomp) = paerml1(1:kproma,1:klev,jcomp) + &
             pgas(1:kproma,1:klev,jc,jm)

     END DO ! jc

     IF (lmass_diag) call sum_mass3(kproma,klev,paerml1,zmass_pre2)
     !--------------------------------------------------------------------

     !--- Assign T [K] and RH [0-1]:
     DO jk=1,klev
     DO jl=1,kproma
        TT(jl,jk) = MAX(MIN(ptemp(jl,jk),TMAX),TMIN)
        RH(jl,jk) = MAX(MIN(prhum(jl,jk),RHMAX),RHMIN)
     END DO
     END DO

     !--- Calculate saturation water vapor mass [g m-3 (air)] and RH [0-1]:
     jcomp = nwh2o(0) ! H2O vapor
     IF(lgh2o .AND. jcomp > 0) THEN
     DO jk=1,klev
        ! Update RH - SM check why RH becomes zero after one mode !
        DO jl=1,kproma
           zmcoat(jl,jk)   = ZERO
           ! Saturation water mass [g m-3 (air)]
           it = NINT(ptemp(jl,jk)*FACp3)
           it = MAX(MIN(it,jptlucu2),jptlucu1)
#ifdef _DEBUG
           WRITE(*,*) 'input: jl,jk,jm,it,tlucuaw(it),ptemp(jl,jk),ppress(jl,jk),prho(jl,jk),psat(jl,jk),prhum(jl,jk) =',&
                              jl,jk,jm,it,tlucuaw(it),ptemp(jl,jk),ppress(jl,jk),prho(jl,jk),psat(jl,jk),prhum(jl,jk)
#endif
           psat(jl,jk) = tlucuaw(it)/ppress(jl,jk) * prho(jl,jk)
           IF(psat(jl,jk) < REALZERO)  psat(jl,jk) = REALZERO
           IF(lah2o) &
           ! [0-1]               [umol/m3 (air)]                 / [g/m3 (air)]
           prhum   (jl,jk) = paerml(jl,jk,jcomp) * FACm6 * mwh2o / psat(jl,jk)
           IF(prhum(jl,jk) < ZERO ) prhum(jl,jk) = ZERO
        END DO
        DO jl=1,kproma
            !  meteorology
            ! limit T  [K]
            TT(jl,jk) = MAX(MIN(ptemp(jl,jk),TMAX),TMIN)
            ! limit RH [0-1]
            RH(jl,jk) = MAX(MIN(prhum(jl,jk),RHMAX),RHMIN)
#ifdef _DEBUG
            WRITE(*,*) 'output: jl,jk,jm,TT(jl,jk),RH(jl,jk),psat(jl,jk),prhum(jl,jk),paerml(jl,jk,jcomp),jcomp =',&
                                jl,jk,jm,TT(jl,jk),RH(jl,jk),psat(jl,jk),prhum(jl,jk),paerml(jl,jk,jcomp),jcomp
#endif
        END DO
     END DO

!!$     ! mz_ht_20090601+
     DO jk=1,klev
       DO jl=1,kproma
         psphum(jl,jk) = paerml(jl,jk,jcomp) * FACm6 * mwh2o / prho(jl,jk)
       ENDDO
     ENDDO
     IF(lsubrh) &
     CALL SUBGRID_SCALE_RHUM(kproma, klev,                                &
                             ppress, ppressi, ptemp, psphum, pvervel, &
                             frac(:,:), prhum, RH_rest )
!!$     DO jk=nlev1,nlev2
!!$       DO jl=1,kproma
!!$         IF (RH(jl,jk,jm) > 0.90_dp) print*, "RH value",jl,jk,jm,prhum(jl,jk), frac(jl,jk,1), ptemp(jl,jk), RH_rest(jl,jk)
!!$       ENDDO
!!$     ENDDO
     ! mz_ht_20090601-
     END IF ! lgh2o/jcomp

!#ifdef _DEBUG
#ifdef _INFO
  ! debugging only
  ! check all species  [umol m-3 (air)]

     DO jc=1,spec_number
        kcomp = species(jc)%aermlidx(jm)
        IF (kcomp == 0) CYCLE
        DO jk = nlev1,nlev2
           DO jl = 1,kproma
              IF(l_io .AND. jk == klev) WRITE(*,*)  ' 1.1 gmxe_thermo - ', & ! remove
                   ' input gmxe_thermo - paerml(jl,jk,kcomp) = ',jl,jk,kcomp,paerml(jl,jk,kcomp),jc,jm,species(jc)%name
           END DO ! jl
        END DO ! jk
     END DO ! jc
#endif

     !
     !--------------------------------------------------------------------
     !

     !
     km = jm
!CDIR NOIEXPAND
     !

     IF (LEQM) THEN
       DO jc=1,ngas
        ! condensed gases (by pure diffusion limited processes)
         xgases(1:kproma,1:klev,km,1,jc) = pgas (1:kproma,1:klev,jc,km)
         pgas (1:kproma,1:klev,jc,km)     = 0._dp

       END DO
     ENDIF
     !
      !--------------------------------------------------------------------
     !

     !
     mask (:,:,jm) = 0
     zmask(:,:,jm) = ZERO
     !
     !--------------------------------------------------------------------
     ! set mask for nucleation mode
     !

     IF (leqm .AND. jm == 1 .AND. lnucl .AND. nnucl == 0) THEN

     ! Dynamical and thermodynamical apporach:
     ! Allow nucleation of gas phase species only, if total aerosol number
     ! is below cmin_aernl; otherwise condensation is preferred (assumed).
     ! Currently only the nucleation thresholds of sulfuric acid is given.
     ! All other gas phase species are thus only allowed to condense during
     ! the nucleation process of H2SO4 by which the number of moles that
     ! condense depends on the neutralization reactions.
     ! Implicitly we assume that model dynamics determine first order effects.
     !
     ! Note that this part is experimental. Results should be treated with care!
     !
     ![N cm-3 (air)]
     znumber(:,:) = zero
     Do jmod=1,nmod
       DO jk = 1,klev
         DO jl = 1,kproma
           znumber(jl,jk) = znumber(jl,jk) + paernl(jl,jk,jmod)
         END DO
       END DO
     END DO
     DO jk = 1,klev
       DO jl = 1,kproma
         ! If nucleation mode particles already exist,
         ! this mode should be calculateed as well, since the
         ! condensation on these particles is aleo preferred.
         !
         IF (paernl(jl,jk,jm) > FACp0) znumber(jl,jk) = 0._dp
         !
         IF (znumber(jl,jk) < FACp0                     &     ! [N   cm-3 (air)]
           .AND. paerml (jl,jk,nh2so4(0)) > ZERO       &     ! [umol m-3 (air)]
           .AND. paerml (jl,jk,nwh2o (0)) > ZERO )  THEN  ! Water vapor

           mask (jl,jk,jm) = 1
           zmask(jl,jk,jm) = FACp0
         END IF
       END DO ! jl
     END DO ! jk
     !
     END IF ! (leqm .AND. jm == 1 .AND. lnucl .AND. nnucl == 0)

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_thermo_region_2")
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_3")
#endif

     !
     !--------------------------------------------------------------------
     ! set mask for soluble modes (above nucleation)
     !
     IF (jm >= 1 .AND. jm <= nsoluble) THEN
!     IF (jm > NS .AND. jm <= nsoluble) THEN
        mask (:,:,jm) = 1
        zmask(:,:,jm) = FACp0
     END IF
     !
     !--------------------------------------------------------------------

     !--------------------------------------------------------------------
     ! set mask for insoluble modes (mass transfer to soluble modes)
     !--------------------------------------------------------------------
     IF (leqm .AND. lcoat .AND. jm > nsoluble) THEN

     ! Allow coating/ageing of primary species only, if total aerosol mass
     ! is above cmin_aerml; otherwise condensation is neglected.

     ! skip calculations in case of no primary particles
     lcoati = .FALSE.
     ! sum mass of primary particles [ug m-3 (air)]
     zmcoat(:,:) = ZERO

     DO jc=1,spec_number
       kcomp = species(jc)%aermlidx(jm)
       IF (kcomp == 0 .OR. kcomp == species(jc)%aernlidx(jm)) CYCLE
       DO jk = 1,klev
         DO jl = 1,kproma
           zmcoat(jl,jk) = zmcoat(jl,jk) &
                         + paerml(jl,jk,kcomp) * species(jc)%molmass
         ENDDO
       ENDDO
     ENDDO



     DO jk = 1,klev
     DO jl = 1,kproma
        IF (zmcoat(jl,jk) < cmin_aerml) CYCLE
        ! at least one possible coating event within this vector loop
        mask (jl,jk,jm) = 1
        zmask(jl,jk,jm) = FACp0
        lcoati         = .TRUE.
     END DO ! jl
     END DO ! jk

     END IF ! (leqm .AND. lcoat .AND. jm > nsoluble)

     !------------------------------------------------------------------------

     ! add aerosols         [umol m-3 (air)]
     ! anions+acids
     DO jc=1,nanions
        kcomp = td%anion(jc)%aerml_idx(jm)
        IF (kcomp == 0 ) CYCLE
        DO jk = 1,klev
           DO jl = 1,kproma
              xanions(jl,jk,jm,1,jc) = paerml(jl,jk,kcomp)
           END DO ! jl
        END DO ! jk
     END DO ! jc


     ! cations
    DO jc = 1,ncations
       kcomp = td%cation(jc)%aerml_idx(jm)
       IF (kcomp == 0) CYCLE
       DO jk = 1,klev
         DO jl = 1,kproma
           xcations(jl,jk,jm,1,jc) = paerml(jl,jk,kcomp)
         END DO ! jl
       END DO ! jk
     END DO ! jc

     ! neutral solutes
     DO jc = 1,nsolutes
       kcomp = td%solute(jc)%aerml_idx(jm)
       IF (kcomp == 0) CYCLE
       DO jk = 1,klev
         DO jl = 1,kproma
           xsolutes(jl,jk,jm,1,jc) = paerml(jl,jk,kcomp)
         END DO ! jl
       END DO ! jk
     END DO ! jc

     ! aerosol water
     kcomp=nwh2o(jm)    ! aerosol water
     IF (kcomp /= 0) THEN
       DO jk = 1,klev
         DO jl = 1,kproma
           zwatbulk(jl,jk)  = paerml(jl,jk,kcomp)
           ! SM assign aerosol water of the previous time step
           WH2O(jl,jk,jm) = paerosol(jl,jk,jm,iPX) ! aerosol H2O [ug m-3 (air)]
         END DO ! jl
       END DO ! jk
     END IF

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_thermo_region_3")
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_4")
#endif



     !--------------------------------------------------------------------
     !--- 2. Call aerosol composition routine ----------------------------
     !--------------------------------------------------------------------

     IF(neqm < 1) lcoati = .FALSE.

     !--------------------------------------------------------------------
     IF (leqm) THEN
     !--------------------------------------------------------------------
     ! gas phase species    [umol m-3 (air)]
     ! anions/acids
! mz_ht_20090305+
       ZHp (:,:,jm)  = 0._dp
       ZPLUS(:,:)    = 0._dp
       ZMIN(:,:)     = 0._dp
       ZDIFF(:,:,jm) = 0._dp


! adding up the gas phase compounds to the aerosol
! mz_ht_20090305-
       DO jc=1,nanions
          jcomp = td%anion(jc)%gas_idx
          IF (jcomp == 0) CYCLE
          IF (.NOT. td%gas(jcomp)%ltreat(jm)) CYCLE
          zcharge = td%anion(jc)%charge
          nhy = td%Hp_idx
          IF (.NOT. td%cation(nhy)%ltreat(jm)) nhy = 0
          DO jk = 1,klev
             DO jl = 1,kproma
                IF(mask (jl,jk,jm) == 0) CYCLE
                xanions(jl,jk,jm,1,jc) = xanions(jl,jk,jm,1,jc) + &
                                         xgases(jl,jk,km,1,jcomp)
                ! charge is negative for anions, to add H+ a minus is required
                IF (NHY /= 0) &
                xcations(jl,jk,jm,1,nhy) = xcations(jl,jk,jm,1,nhy) - &
                                           zcharge * xgases(jl,jk,km,1,jcomp)
                xgases(jl,jk,km,1,jcomp) = 0._dp
             END DO
          END DO
       END DO

       DO jc = 1,ncations
          jcomp = td%cation(jc)%gas_idx
          IF (jcomp == 0) CYCLE
          IF (.NOT. td%gas(jcomp)%ltreat(jm)) CYCLE
          zcharge = td%cation(jc)%charge
          nhy = td%Hp_idx
          IF (.NOT. td%cation(nhy)%ltreat(jm)) nhy = 0
          DO jk = 1,klev
             DO jl = 1,kproma
                IF(mask (jl,jk,jm) == 0) CYCLE
                xcations(jl,jk,jm,1,jc) = xcations(jl,jk,jm,1,jc) + &
                     xgases(jl,jk,km,1,jcomp)
                IF (NHY /= 0) &
                     xcations(jl,jk,jm,1,nhy) = xcations(jl,jk,jm,1,nhy) - &
                     zcharge * xgases(jl,jk,km,1,jcomp)
                xgases(jl,jk,km,1,jcomp) = 0._dp
             END DO ! jl
          END DO ! jk
       END DO ! jc

     ! convert  [umol m-3(air)] -> [mol m-3(air)]
     ! FACm6  = 1.e-6_dp

     ! input total gas+aerosol    [mol m-3(air)]
     ! anions/acids
     DO jc = 1,nanions
       jcomp = td%anion(jc)%aerml_idx(jm)
       IF( jcomp == 0 ) CYCLE
        DO jk = 1,klev
        DO jl = 1,kproma
           IF(mask (jl,jk,jm) == 0) CYCLE
           xanions(jl,jk,jm,1,jc) = xanions (jl,jk,jm,1,jc) * FACm6
        END DO ! jl
        END DO ! jk
     END DO ! jc
     ! cations
     DO jc = 1,ncations
       jcomp = td%cation(jc)%aerml_idx(jm)
       IF( jcomp == 0 ) CYCLE
       DO jk = 1,klev
         DO jl = 1,kproma
           IF(mask (jl,jk,jm) == 0) CYCLE
           xcations(jl,jk,jm,1,jc) = xcations(jl,jk,jm,1,jc) * FACm6
         END DO ! jl
       END DO ! jk
     END DO ! jc
     ! neutral solutes
     DO jc = 1,nsolutes
       jcomp = td%solute(jc)%aerml_idx(jm)
       IF( jcomp == 0 ) CYCLE
       DO jk = 1,klev
         DO jl = 1,kproma
           IF(mask (jl,jk,jm) == 0) CYCLE
           xsolutes(jl,jk,jm,1,jc) = xsolutes(jl,jk,jm,1,jc) * FACm6
         END DO ! jl
       END DO ! jk
     END DO ! jc

     SELECT CASE (NEQM)
     CASE(1)
        ! EQSAM4CLIM
        do jk=1,klev
           do jl=1,kproma
              IF(mask (jl,jk,jm) == 0) CYCLE
              CALL EQSAM4CLIM_INTERFACE(xanions(jl,jk,jm,:,:),     &
                                        xcations(jl,jk,jm,:,:),    &
                                        ygases(:),                 &
                                        RH(jl,jk), TT(jl,jk),      &
                                        WH2O(jl,jk,jm),            &
                                        paerosol(jl,jk,jm,ifsulf), &
                                        prdry(jl,jk,jm), xdiag, 9, &
                                        jl, jk, jm, ldry, lhyster)
              paerosol(jl,jk,jm,ispm) = xdiag(1)
              paerosol(jl,jk,jm,iapm) = xdiag(2)
              paerosol(jl,jk,jm,ipms) = xdiag(3)
              paerosol(jl,jk,jm,ipmt) = xdiag(4)
              paerosol(jl,jk,jm,irho) = xdiag(5)
              paerosol(jl,jk,jm,ivol) = xdiag(6)
              paerosol(jl,jk,jm,ipH)  = xdiag(7)
              paerosol(jl,jk,jm,iGF)  = xdiag(8)
!              paerosol(jl,jk,jm,iHp) = xdiag(9)  ! does not exist yet
           END do
        END do



     CASE(2)
        ! ISORROPIA 2

!        print*, "before Iso_inter: ",jm, MINVAL(XANIONS(:,:,jm,:,:)), MAXVAL(XANIONS(:,:,jm,:,:)), &
!             MINVAL(XCATIONS(:,:,jm,:,:)), MAXVAL(XCATIONS(:,:,jm,:,:))
       do jk=1,klev
         do jl=1,kproma
           IF(mask (jl,jk,jm) == 0) CYCLE
           CALL ISORROPIA_INTERFACE(xanions(jl,jk,jm,:,:),     &
                                    xcations(jl,jk,jm,:,:),    &
                                    ygases(:),                 &
                                    RH(jl,jk), TT(jl,jk),      &
                                    WH2O(jl,jk,jm),            &
                                    paerosol(jl,jk,jm,ifsulf), &
                                    jl, jk, jm, ldry, lhyster)
           xgases(jl,jk,km,1,:) = xgases(jl,jk,km,1,:) + ygases(:)
           xsolutes(jl,jk,jm,1,:) = xsolutes(jl,jk,jm,1,:) / FACm6
         end do
       end do
!       print*, "after Iso_inter: ", jm,MINVAL(XANIONS(:,:,jm,:,:)), MAXVAL(XANIONS(:,:,jm,:,:)), &
!            MINVAL(XCATIONS(:,:,jm,:,:)), MAXVAL(XCATIONS(:,:,jm,:,:)), MINVAL(WH2O(:,:,jm)), MAXVAL(WH2O(:,:,jm))
     END SELECT

     ! output gas, liquids, solids [umol m-3(air)]
     ! gases

     ! xgases should already contain all gaseous material
     ! for eqsam the remaining gases should be added to xgases in the eqsam interface

     ! anions/acids     (lumped liquids + solids)
     DO jc = 1,nanions
       jcomp = td%anion(jc)%aerml_idx(jm)
       IF( jcomp == 0 ) CYCLE
       DO jk = 1,klev
         DO jl = 1,kproma
           IF(mask (jl,jk,jm) == 0) CYCLE
           xanions (jl,jk,jm,1,jc) = xanions (jl,jk,jm,1,jc) + xanions (jl,jk,jm,2,jc)
        END DO ! jl
        END DO ! jk
     END DO ! jc
     ! cations          (lumped liquids + solids)
     DO jc = 1,ncations
       jcomp = td%cation(jc)%aerml_idx(jm)
       IF( jcomp == 0 ) CYCLE
!       IF (trim(td%cation(jc)%name) == 'NH4p') print*, "out", MAXVAL( xcations(:,19,jm,1,jc) ), &
!         MAXVAL( xcations(:,19,jm,2,jc) ), MAXVAL(xgases(:,19,jm,1,td%cation(jc)%gas_idx))
       DO jk = 1,klev
         DO jl = 1,kproma
           IF(mask (jl,jk,jm) == 0) CYCLE
           xcations(jl,jk,jm,1,jc) = xcations(jl,jk,jm,1,jc) + xcations(jl,jk,jm,2,jc)
        END DO ! jl
        END DO ! jk
     END DO ! jc
     ! neutral solutes  (lumped liquids + solids)
     DO jc = 1,nsolutes
       jcomp = td%solute(jc)%aerml_idx(jm)
       IF( jcomp == 0 ) CYCLE
       DO jk = 1,klev
         DO jl = 1,kproma
           IF(mask (jl,jk,jm) == 0) CYCLE
           xsolutes(jl,jk,jm,1,jc) = xsolutes(jl,jk,jm,1,jc) + xsolutes(jl,jk,jm,2,jc)
         END DO ! jl
       END DO ! jk
     END DO ! jc

     IF(.NOT. ldry) THEN
     ! aerosol water of the current time step [umol m-3 (air)]
        DO jk = 1,klev
           DO jl = 1,kproma
              IF(mask (jl,jk,jm) == 0) CYCLE
              zwatbulk(jl,jk) = WH2O(jl,jk,jm) / mwh2o
           END DO ! jl
        END DO ! jk
     END IF ! ldry

! Add contribution of aerosol water caused by the bulk species not
! considered in thermodynamics, namely bulk seasalt and potentially OC

     IF(lbulk_aerosol_water) THEN
       IF (jm <= nsoluble) &
       CALL calc_wateruptake_bulk(kproma, klev, jm, naertot,     &
                                  paerml, zwatbulk, WH2O, prhum, &
                                  paerosol(1:kproma,1:klev,jm,ioc_water))

     END IF ! lbulk_aerosol_water

     !--------------------------------------------------------------------
     ELSE   ! .NOT. leqm
     !--------------------------------------------------------------------
     ! aerosol water of the current time step
     DO jk = 1,klev
     DO jl = 1,kproma
        WH2O(jl,jk,jm) = zwatbulk(jl,jk) * mwh2o ! aerosol H2O [ug m-3 (air)]
     END DO ! jl
     END DO ! jk

     !--------------------------------------------------------------------
     END IF ! (leqm)
     !--------------------------------------------------------------------

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_thermo_region_4")
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_5")
#endif

     !--------------------------------------------------------------------

     ! Limit calculated aerosol water concentration by available water vapor
     ! concentration that is derived from specific humidity. The water vapor
     ! concentration is corrected by the condensed water mass and optionally
     ! the specific humidity is updated (otherwise equilibrium is assumed).
     ztotwat(:,:) = 0._dp
     jcomp = nwh2o(0) ! H2O vapor
     kcomp = nwh2o(jm) ! H2O aerosol
     IF(lah2o .AND. jcomp > 0) THEN
       IF (kcomp > 0) THEN
         DO jk = 1,klev
           DO jl=1,kproma
           ! [umol m-3 (air)]
           ! max. water mass at saturation [umol m-3 (air)] = (absolute hum. [umol m-3 (air)] / rel. hum.)
             ztotwat(jl,jk)  = paerml(jl,jk,jcomp) + paerml(jl,jk,kcomp)
             zwatbulk(jl,jk) = MAX(zero,MIN(zwatbulk(jl,jk),ztotwat(jl,jk)))
           END DO
         END DO
       ELSE
         DO jk = 1,klev
           DO jl=1,kproma
             ztotwat(jl,jk)  = paerml(jl,jk,jcomp)
             zwatbulk(jl,jk) = MAX(zero,MIN(zwatbulk(jl,jk),ztotwat(jl,jk)))
           END DO
         ENDDO
       END IF
     END IF


#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_thermo_region_5")
#endif

     !--------------------------------------------------------------------
     !--- 3. Assign output values per mode:
     !--------------------------------------------------------------------
     !
#ifdef _FTRACE
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_6")
#endif
     ! gas phase species [umol m-3 (air)]
     DO jc = 1,ngas       ! ngas = nacids + 1 (ammonia)
       jcomp = td%gas(jc)%aerml_idx(0)
       IF (jcomp == nwh2o(0) .OR. jcomp == 0) CYCLE ! exclude H2O vapor
       paerml(1:kproma,1:klev,jcomp) = paerml(1:kproma,1:klev,jcomp) + &
                                       xgases(1:kproma,1:klev,km,1,jc)
     END DO ! jc
     !
     !--------------------------------------------------------------------
     !
     IF (jm <= nsoluble) THEN
     !
     ! aerosol species   [umol m-3 (air)]
     ! anions/acids
     DO jc = 1,nanions
        kcomp = td%anion(jc)%aerml_idx(jm)
        IF (kcomp == 0) CYCLE
        DO jk = 1,klev
           DO jl = 1,kproma
              paerml (jl,jk,kcomp) = xanions(jl,jk,jm,1,jc)
           END DO ! jl
        END DO ! jk
     END DO ! jc
     ! cations
     DO jc = 1,ncations
        kcomp = td%cation(jc)%aerml_idx(jm)
        IF (kcomp == 0) CYCLE
        DO jk = 1,klev
           DO jl = 1,kproma
              paerml (jl,jk,kcomp) = xcations(jl,jk,jm,1,jc)
           END DO ! jl
        END DO ! jk
     END DO ! jc
     ! neutral solutes
     DO jc = 1,nsolutes
        kcomp = td%solute(jc)%aerml_idx(jm)
        IF (kcomp == 0) CYCLE
        DO jk = 1,klev
           DO jl = 1,kproma
              paerml (jl,jk,kcomp) = xsolutes(jl,jk,jm,1,jc)
           END DO ! jl
        END DO ! jk
     END DO ! jc
     ! aerosol water
     kcomp = nwh2o(jm) ! H2O aerosol
     IF (kcomp > 0) THEN
        DO jk = 1,klev
           DO jl = 1,kproma
              paerml(jl,jk,kcomp) = zwatbulk(jl,jk)
           END DO ! jl
        END DO ! jk
     END IF
     !

     END IF ! (jm <= nsoluble)
     !
     !--------------------------------------------------------------------
     !
     IF (lcoat .AND. jm > nsoluble) THEN
     !
     ! Update tracers from coating events (numbers are updated below)

     ! aerosol species [umol m-3 (air)]
     !
        km = 1 + jm - nsoluble  ! add to resepctive soluble modes
     !
        zdrybulk(:,:) = 0._dp
     ! anions/acids
        DO jc = 1,nanions
           kcomp=td%anion(jc)%aerml_idx(km) !   soluble
           IF (kcomp == 0) CYCLE
           DO jk = 1,klev
              DO jl = 1,kproma
                 paerml (jl,jk,kcomp) = &
                      paerml (jl,jk,kcomp) + xanions(jl,jk,jm,1,jc)
                 zdrybulk(jl,jk) = zdrybulk(jl,jk) &
                      + xanions(jl,jk,jm,1,jc) * td%anion(jc)%molmass
              END DO ! jl
           END DO ! jk
        END DO ! jc
        ! cations
        DO jc = 1,ncations
           kcomp=td%cation(jc)%aerml_idx(km) !   soluble
           IF (kcomp == 0) CYCLE
           DO jk = 1,klev
              DO jl = 1,kproma
                 paerml (jl,jk,kcomp) = &
                      paerml (jl,jk,kcomp) + xcations(jl,jk,jm,1,jc)
                 zdrybulk(jl,jk) = zdrybulk(jl,jk) + &
                      xcations(jl,jk,jm,1,jc) * td%cation(jc)%molmass
              END DO ! jl
           END DO ! jk
        END DO ! jc
        ! neutral solutes
        DO jc = 1,nsolutes
           kcomp=td%solute(jc)%aerml_idx(km) !   soluble
           IF (kcomp == 0 .AND. kcomp /= nwh2o(km)) CYCLE
           DO jk = 1,klev
              DO jl = 1,kproma
                 paerml (jl,jk,kcomp) = &
                      paerml (jl,jk,kcomp) + xsolutes(jl,jk,jm,1,jc)
                 zdrybulk(jl,jk) = zdrybulk(jl,jk) + &
                      xsolutes(jl,jk,jm,1,jc) * td%solute(jc)%molmass
              END DO ! jl
           END DO ! jk
        END DO ! jc
        ! aerosol water
        kcomp = nwh2o(km)
        DO jk = 1,klev
           DO jl = 1,kproma
              paerml(jl,jk,kcomp)     = paerml (jl,jk,kcomp) + zwatbulk(jl,jk)
           END DO ! jl
        END DO ! jk
        !
     END IF ! (lcoat .AND. jm > nsoluble)

     !--------------------------------------------------------------------

! Updating water vapour concentration, just after updating the aerosol water
     jcomp = nwh2o(0)
     kcomp = nwh2o(jm)
     DO jk = 1,klev
        DO jl = 1,kproma
           paerml(jl,jk,jcomp) = MAX(zero, ztotwat(jl,jk) - zwatbulk(jl,jk))
        ENDDO
     ENDDO


     CALL Hp_OHm_check(kproma, klev, paerml,paerosol,1)



     DO jk=1,klev
       DO jl=1,kproma

!!$        DO jc=1,spec_number
!!$          IF (ABS( pre(jc) - post(jc) ) > (10._dp * SPACING(post(jc))) ) THEN
!!$            print*, "change specieswise: ", jm, jc, pre(jc), post(jc),       &
!!$              jl, jk,ABS( pre(jc) - post(jc))
!!$            DO jmod=1,nmod
!!$              idx1 = species(jc)%aermlidx(jmod)
!!$              if (idx1 == 0 .OR. idx1 == species(jc)%aernlidx(jmod)) cycle
!!$              print*, "per mode: ", jmod, jc, jl, jk, idx1, &
!!$                paerml1(jl,jk,idx1), paerml2(jl, jk, idx1)
!!$            END DO
!!$          END IF
!!$        END DO
       END DO
     END DO



#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_thermo_region_6")
#endif

#ifdef _FTRACE
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_7")
#endif

     !--------------------------------------------------------------------
     !--- 3.1 Calculate Coating of primary particles
     !--------------------------------------------------------------------
     IF (lcoati .AND. jm > nsoluble) THEN
     !--------------------------------------------------------------------

#ifdef _FTRACE
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_8")
#endif

       km = 1 + jm - nsoluble  ! add to resepctive soluble modes

       zcoat(:,:) = 0._dp
       DO jk = 1,klev
         DO jl = 1,kproma
           ! determmine surface area of existing particles [cm^2/cm^3(air)]
           zsurf(jl,jk) = paernl(jl,jk,jm) * 4._dp * pi          &
                        * prdry(jl,jk,jm)*prdry(jl,jk,jm)        &
                        * cmr2ras(jm) * cmr2ras(jm)
           if (zsurf(jl,jk) > 0._dp ) THEN
             ! determine surface area of coating particles
             xhelp1 = 0._dp
             ! zdrybulk is in ?g m-3(air) = 1e-12 g cm-3(air)
             ! zdrybulk is multiplied with molar mass per compound (in the
             !   calculation of zdrybulk above) and then divided by the
             !   molar mass of sulphate to yield equivalent moles of
             !   sulphate
             ! csurf_molec = average cross-section of a single H2SO4 molecule
             ! use zdrybulk as equivalent moles of H2SO4
             !   to determine the surface area which can be covered by the coating
             !   material
             xhelp1 = 1.e-12_dp * zdrybulk(jl,jk) / &
               species(spec_idx_so4mm)%molmass * avo * csurf_molec

             zcoat(jl,jk) = MIN( xhelp1/(zsurf(jl,jk) &
                               * clayerthickness), 1._dp)

           END if
         END DO
       END DO

       ! moving stuff accordingly from the insoluble mode to the soluble mode
       jmod = jm - nsoluble + ndiff
       ! particle numbers
!       print*, "thermo coating paernl before :", paernl(1,31,:)
       DO jk=1,klev
         DO jl=1,kproma
           paernl(jl,jk,jmod)  = paernl(jl,jk,jmod)  + &
                                 paernl(jl,jk,jm) * zcoat(jl,jk)
           paernl(jl,jk,jm)    = paernl(jl,jk,jm) * (1._dp - zcoat(jl,jk))
         END DO
       END DO
!       print*, "thermo coating paernl after :", paernl(1,31,:)

       DO jc=1,spec_number
          kcomp = species(jc)%aermlidx(jm)
          IF (kcomp == 0 .OR. kcomp == species(jc)%aernlidx(jm)) CYCLE
          jcomp = species(jc)%aermlidx(jmod)
          IF (jcomp == 0 .OR. jcomp == species(jc)%aernlidx(jmod)) CYCLE
          DO jk=1,klev
             DO jl=1,kproma
! think whether zcoat is ok here or if it should be modified from
! particle numbers to the mass distributions
! this will reduce insolube species lifetime!!!!
                paerml(jl,jk,jcomp) = paerml(jl,jk,jcomp) + &
                                      paerml(jl,jk,kcomp) * zcoat(jl,jk)

                paerml(jl,jk,kcomp) = paerml(jl,jk,kcomp) * (1._dp - zcoat(jl,jk))
             ENDDO
          END DO
       END DO

     !--------------------------------------------------------------------
     END IF ! (lcoati .AND. jm > nsoluble)
     !--------------------------------------------------------------------

     IF (lmass_diag) THEN
       paerml2 = paerml
       call sum_mass3(kproma,klev,paerml2,zmass_post2)
       DO jk=1,klev
         DO jl=1,kproma
           pre(:)  = 0._dp
           post(:) = 0._dp
           DO jmod=0,nmod
             DO jc=1,spec_number
               idx1 = species(jc)%aermlidx(jmod)
               if (idx1 == 0 .OR. idx1 == species(jc)%aernlidx(jmod)) cycle
               pre(jc)  = pre(jc) + paerml1(jl,jk,idx1)
               post(jc) = post(jc) + paerml2(jl,jk,idx1)
             END DO
           END DO
           IF (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) >                &
             (10._dp * SPACING(zmass_post2(jl,jk))) ) THEN
             print*, "mass changes itemwise; thermo box: ",jl,jk," Error: (%)",  jm, &
               (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) /               &
               zmass_post2(jl,jk)*100._dp),                               &
               zmass_pre2(jl,jk), zmass_post2(jl,jk)
             DO jc=1,spec_number
               IF (ABS(pre(jc) - post(jc)) > 1.e-30_dp ) THEN
                 DO jmod=0,nsoluble
                   idx1 = species(jc)%aermlidx(jmod)
                   print*, trim(species(jc)%name), "  per mode: ", jmod, jc, jl, jk, idx1, &
                     paerml1(jl,jk,idx1)*species(jc)%molmass, &
                     paerml2(jl, jk, idx1)*species(jc)%molmass, &
                     (paerml1(jl,jk,idx1)- paerml2(jl, jk, idx1))*species(jc)%molmass
                 END DO
               ENDIF
             END DO
           END IF
         END DO
       END DO
     END IF

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_thermo_region_7")
#endif


     DO jk = 1,klev
     DO jl = 1,kproma
!        paerosol(jl,jk,jm,iTT)        = TT        (jl,jk,jm) ! air temperature [K]
!        paerosol(jl,jk,jm,iRH)        = RH        (jl,jk,jm) ! relative humidity [%]
!        paerosol(jl,jk,jm,iRH)        = paerosol  (jl,jk,jm,iRH) * FACp2
!        paerosol(jl,jk,jm,iPX)        = PX        (jl,jk,jm) ! aerosol water history  [0=solid, else=wet]
!        paerosol(jl,jk,jm,iawano)     = thmss     (jl,jk,jm,fe,awano) ! activity coefficient of NH4NO3 [-]
!        paerosol(jl,jk,jm,iZIONIC)    = ZIONIC    (jl,jk,jm) ! ionic strength (aq) [mol/kg]
!        paerosol(jl,jk,jm,iRHO)       = RHO       (jl,jk,jm) ! density (bulk mass) [g/m3]
!        paerosol(jl,jk,jm,iVOL)       = VOL       (jl,jk,jm) ! bulk volume [ucm3 m-3 (air)]

!        paerosol(jl,jk,jm,iPH)        = PH        (jl,jk,jm) ! aerosol pH [log H+]
!        paerosol(jl,jk,jm,isPM)       = sPM       (jl,jk,jm) ! total solid  matter [umol m-3 (air)]
!        paerosol(jl,jk,jm,iaPM)       = aPM       (jl,jk,jm) ! total liquid matter [umol m-3 (air)]
!        paerosol(jl,jk,jm,iPMt)       = PMt       (jl,jk,jm) ! total PM (liquids & solids) [ug   m-3 (air)]
!        paerosol(jl,jk,jm,iPMs)       = PMs       (jl,jk,jm) ! total PM (solids, treated by EQSAM) [ug   m-3 (air)]
!        paerosol(jl,jk,jm,irhdcr_sum) = rhdcr_sum (jl,jk,jm) ! crystallization RH of mixed solution [-]
!        paerosol(jl,jk,jm,irhdms_sum) = rhdms_sum (jl,jk,jm) ! deliquescence   RH of mixed solution [-]
!        paerosol(jl,jk,jm,iRHcr)      = RHcr      (jl,jk,jm) ! lowest crystallization RH in mixed solution [-]
!        paerosol(jl,jk,jm,iRHD)       = RHD       (jl,jk,jm) ! lowest deliquescence   RH in mixed solution [-]
        paerosol(jl,jk,jm,iWH2O)      = WH2O      (jl,jk,jm) ! aerosol Water (aq) [ug m-3 (air)]

        ! hygroscopic growth factor [-]
        paerosol(jl,jk,jm,iGF)        = 1._dp
        IF(prdry(jl,jk,jm) > REALZERO) &
             paerosol(jl,jk,jm,iGF)        = prwet(jl,jk,jm) / prdry(jl,jk,jm)     ! hygroscopic growth factor [-]

!        paerosol(jl,jk,jm,iNCawano)   = thmss     (jl,jk,jm,nc,awano) ! solid ammonium nitrate [umol m-3 (air)]
!        paerosol(jl,jk,jm,iNAawano)   = thmss     (jl,jk,jm,na,awano) ! aqueous ammonium nitrate [umol m-3 (air)]
!        paerosol(jl,jk,jm,iNAawasu)   = thmss     (jl,jk,jm,na,awasu) ! aqueous ammonium sulfate [umol m-3 (air)]

     END DO ! jl
     END DO ! jk

! Calculate the Kappa (water uptake) vaule, following Petters and Kreidenweis
! Calcualate dry volume fraction per compound.
   END DO

   DO jm=1,nmod
     IF(lkappa .AND. neqm >= 1) THEN

        call calc_kappa(kproma,klev,jm, naertot, diagout, paerml, &
                        paernl(:,:,jm), RH, TT, prdry(:,:,jm))
        paerosol(:,:,jm,ikappa)        = diagout(:,:,1)
        paerosol(:,:,jm,ikappa_insol)  = diagout(:,:,2)
        paerosol(:,:,jm,ika_vol)       = diagout(:,:,3)
        paerosol(:,:,jm,ika_vol_insol) = diagout(:,:,4)
        paerosol(:,:,jm,ika_water)     = diagout(:,:,5)
        paerosol(:,:,jm,iSc)           = diagout(:,:,6)
        paerosol(:,:,jm,iCCN2)         = diagout(:,:,7)
        paerosol(:,:,jm,iCCN4)         = diagout(:,:,8)
    END IF ! lkappa


     !--------------------------------------------------------------------
  END DO ! jm
     !--------------------------------------------------------------------


#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_thermo_region_1")
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_10")
#endif

     !--------------------------------------------------------------------
  DO jm = 1,nmod
     !--------------------------------------------------------------------

#ifdef _FTRACE
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_11")
#endif

     ! re-assign/calculate aerosol properties (radius, density, number)

     ! Initialization:

     DO jk = 1,klev
     DO jl = 1,kproma
       zfac(jl,jk) = zero
       zm2n(jl,jk) = zero
       zcmm(jl,jk) = zero
       zhelp (jl,jk) = zero
       zdrymass(jl,jk) = zero
       zwetmass(jl,jk) = zero
       zwatmass(jl,jk) = zero
       zdryvol(jl,jk)  = zero
       zwetvol(jl,jk)  = zero
       zdryrad(jl,jk)  = zero
       zwetrad(jl,jk)  = zero
       zdryden(jl,jk)  = FACp0 ! pure water [g m-3]
       zwetden(jl,jk)  = zdryden(jl,jk)
       zdryvol_mean(jl,jk) = zero
       zwetvol_mean(jl,jk) = zero
       zionbulk    (jl,jk) = zero
       zdrybulkm   (jl,jk) = zero
       zdrybulk    (jl,jk) = zero
       zwatbulkm   (jl,jk) = zero
       zwatbulk    (jl,jk) = zero
       zmbulk      (jl,jk) = zero
       zmbulki     (jl,jk) = zero
     END DO
     END DO
!!! QQQQ test comment
!!$     IF(leqm .AND. neqm >= 1 .AND. jm <= nsoluble) THEN
!!$        DO jk = 1,klev
!!$        DO jl=1,kproma
!!$           zdryden  (jl,jk) = RHO (jl,jk,jm)              ! aerosol dry density  [g   cm-3]
!!$           zdrybulkm(jl,jk) = PMt (jl,jk,jm)              ! total PM (excl. H2O) [ug   m-3 (air)]
!!$           zdrybulk (jl,jk) = sPM (jl,jk,jm)+aPM(jl,jk,jm)! total PM (excl. H2O) [umol m-3 (air)]
!!$           zionbulk (jl,jk) = aPM (jl,jk,jm)              ! aqueous PM           [umol m-3 (air)]
!!$           zdryvol  (jl,jk) = VOL (jl,jk,jm)              ! bulk volume          [ucm3 m-3 (air)]
!!$           zwatbulk (jl,jk) = WH2O(jl,jk,jm)/mwh2o        ! aerosol Water (aq)   [umol m-3 (air)]
!!$        END DO ! jl
!!$        END DO ! jk
!!$     ELSE

        kcomp = species(spec_idx_h2o)%aermlidx(jm) ! H2O aerosol
        IF (kcomp > 0) THEN
        DO jk = 1,klev
        DO jl=1,kproma
!           RHO(jl,jk,jm)    = FACp0
           zwatbulk (jl,jk) = paerml(jl,jk,kcomp)
        END DO
        END DO
        END IF
        DO jc=1,spec_number
          kcomp=species(jc)%aermlidx(jm)
           IF(kcomp == nwh2o(jm) .OR. kcomp == 0. .OR.&
             kcomp == species(jc)%aernlidx(jm)) CYCLE ! exclude H2O aerosol
           DO jk = 1,klev
           DO jl=1,kproma
              ![umol m-3 (air)]= [umol m-3 (air)]
              zdrybulk (jl,jk) = zdrybulk (jl,jk) &
                               + MAX(0._dp,paerml(jl,jk,kcomp))
              ![ug   m-3 (air)]= [umol m-3 (air)]   * [g mol-1]<=>[ug umol-1]
              zdrybulkm(jl,jk) = zdrybulkm(jl,jk) &
                               + MAX(0._dp,paerml(jl,jk,kcomp)) &
                               * species(jc)%molmass
      !        PMt(jl,jk,jm)    = zdrybulkm(jl,jk)
      !        PMs(jl,jk,jm)    = zdrybulkm(jl,jk)
              IF(jm <= nsoluble) THEN
                IF (species(jc)%charge /= 0) &
                  zionbulk (jl,jk) = zionbulk (jl,jk) &
                                   + MAX(0._dp,paerml(jl,jk,kcomp))
              END IF
              ! total   PM  [umol m-3 (air)]
              ![ucm3  m-3 (air)]= [umol m-3 (air)] * [g mol-1]   / [g cm-3]
              zdryvol(jl,jk)   = zdryvol(jl,jk)   &
                               + MAX(0._dp,paerml(jl,jk,kcomp))  &
                               * species(jc)%molmass / species(jc)%density
              IF(zdryvol(jl,jk)> cmin_avol) &
              zdryden(jl,jk)   = zdrybulkm(jl,jk) / zdryvol(jl,jk)

              ! security check
              IF(zdryden(jl,jk) < 0.5_dp) zdryden(jl,jk) = FACp0
              IF(zdryden(jl,jk) > FACp1 ) zdryden(jl,jk) = FACp1
!              RHO(jl,jk,jm)     = zdryden(jl,jk)

           END DO
           END DO
        END DO
!!$     END IF
!QQQ test comment

     DO jk = 1,klev
     DO jl=1,kproma
        ! [ug   m-3 (air)]
        zwatbulkm(jl,jk) = zwatbulk (jl,jk) * mwh2o
        ! total condensed matter, including aerosol water
        ! [umol m-3 (air)]
        zmbulk (jl,jk)   = zdrybulk (jl,jk) + zwatbulk (jl,jk)
        ! [ug   m-3 (air)]
        zmbulki(jl,jk)   = zdrybulkm(jl,jk) + zwatbulkm(jl,jk)
     END DO ! jl
     END DO ! jk

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_thermo_region_11")
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_12")
#endif

     !--------------------------------------------------------------------
     !--- 5. Calculate thermodynamic properties per mode
     !--------------------------------------------------------------------

     !--------------------------------------------------------------------
     mask(:,:,jm) = 0
     DO jk = 1,klev
     DO jl = 1,kproma
       !--------------------------------------------------------------------
       IF (paernl (jl,jk,jm) > cmin_aernl .AND. zmbulk(jl,jk) > cmin_aerml) THEN
           mask(jl,jk,jm)    = 1
       ELSE IF (paernl(jl,jk,jm) <= cmin_aernl .AND. &
                zmbulk(jl,jk) > cmin_aerml) THEN
           mask(jl,jk,jm)    = 2
       END IF

       znumber(jl,jk) = paernl(jl,jk,jm) * FACp6
       zdryrad(jl,jk) = prdry (jl,jk,jm)
       zwetrad(jl,jk) = prwet (jl,jk,jm)

     END DO ! jl
     END DO ! jk

     DO jk = 1,klev
     DO jl = 1,kproma
       !--------------------------------------------------------------------
       IF (mask(jl,jk,jm) > 0) THEN
       !--------------------------------------------------------------------

         !--- 5.1 Calculate aerosol number [N m-3]:

         ! [kg]      =         [g cm-3] -> [kg m-3]     *         [m3]
         zcmm(jl,jk)    =  4._dp/3._dp*pi*zdryden(jl,jk)*1.e3_dp*crdiv_r3(jm)
         !!! zcmm(jl,jk)    =  4._dp/3._dp*pi*RHO(jl,jk,jm)*1.e3_dp*crdiv_r3(jm)  !!KP Need mode info? No as zdryden set within loop

         ! [kg-1]    = 1.   /  [kg]
         zm2n(jl,jk)    = FACp0/zcmm(jl,jk)

!!$         IF(zcmm(jl,jk) == 0.0)THEN
!!$            print*,' zm2n(jl,jk) =',zm2n(jl,jk) ,FACp0,zcmm(jl,jk),zdryden(jl,jk),crdiv_r3(jm),RHO(jl,jk,jm)
!!$         ENDIF

       !--------------------------------------------------------------------
       END IF !
       !--------------------------------------------------------------------
       !--------------------------------------------------------------------
       IF (mask(jl,jk,jm) == 2 .OR. lnumber) THEN
       !--------------------------------------------------------------------
         !--- Convert mass to number for given number median radii and standard deviation
         !         (implicitly by the conversion factor cmr2ram) of the modes.
         !
         !    M2N  = 1/M0 = 1/(4/3 * pi * dens * R0(averageMass)**3)
         !    N(t) = M(t) * M2N * N(t-1)
         !
         ! [N m-3]      = [ug m-3]         => [kg]    * [kg-1]
         znumber(jl,jk) = zdrybulkm(jl,jk)  * FACm9    * zm2n(jl,jk)
         paernl (jl,jk,jm) = znumber(jl,jk) * FACm6

       !--------------------------------------------------------------------
       END IF ! (mask(jl,jk,jm) == 2 .OR. lnumber)
       !--------------------------------------------------------------------
       IF (mask(jl,jk,jm) > 0 .AND. paernl (jl,jk,jm) > cmin_aernl) THEN
       !--------------------------------------------------------------------

         !--- 5.1.3 Total dry aerosol mass (liquids & solids) [ug]:

         ! [ug]          =  [ug m-3 (air)   / [N m-3 (air)]
         zdrymass(jl,jk) = zdrybulkm(jl,jk) / znumber(jl,jk)

         !--- 5.1.4 Total aerosol water mass [ug]:

         IF (.NOT. ldry) &
         ! [ug]          = [ug m-3 (air)    / [N m-3 (air)]
         zwatmass(jl,jk) = zwatbulkm(jl,jk) / znumber(jl,jk)

         !--- 5.2.2  Mean aerosol dry volume [cm3]:
         !   [cm3]             [ucm3 m-3(air)]        / [N m-3(air)]
         zdryvol_mean(jl,jk) = zdryvol(jl,jk) * FACm6 / znumber(jl,jk)

         !--- 5.2.3 Dry radius [cm]:

         ![cm] =                    [cm3] -> [cm]
         zdryrad(jl,jk)=(0.75_dp/pi*zdryvol_mean(jl,jk))**(1._dp/3._dp) * &
                         ram2cmr(jm)

         !--------------------------------------------------------------------
         !--- 5.3.1 Wet (dry mass + water) volume [cm3]:
         !          (assuming volume additivity)

         !    [cm3]                     [g]            / [g cm-3]  +     [cm3]
         zwetvol_mean(jl,jk) = zwatmass(jl,jk) * FACm6 /    Dw     + &
                               zdryvol_mean(jl,jk)

         !--- 5.3.2 Wet radius [cm]:

         zwetrad(jl,jk) = zdryrad(jl,jk)

         IF (zwatmass(jl,jk) > ZERO ) &
         zwetrad(jl,jk)=(0.75_dp/pi*zwetvol_mean(jl,jk))**(1._dp/3._dp) * &
                         ram2cmr(jm)

        !--------------------------------------------------------------------
        END IF ! (mask(jl,jk,jm) > 0)
        !--------------------------------------------------------------------

     END DO ! jl
     END DO ! jk

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_thermo_region_12")
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_13")
#endif

     !--------------------------------------------------------------------
     !--- 5.5 Recalculate aerosol mean radius for the particles with the
     !        standard deviation given by the tracer property.
     !        If the mean radius does not lie in the mode, a mass difference
     !        is transfered to next higher mode (see also comment on 5.4).
     zhelp  (:,:) = ZERO
     zmbulki(:,:) = ZERO
     !--------------------------------------------------------------------
     IF (lsize .AND. jm < nsoluble) THEN
     !--------------------------------------------------------------------
        DO jk = 1,klev
        DO jl = 1,kproma
           !--- Wet radius [cm]:
           !--- 5.5.1 Respective bulk mass of the mode at given aerosol number:
           ![ug m-3 (air)] =        [kg] -> [ug]   *     [N m-3]
           ! [kg]      =         [g cm-3] -> [kg m-3]     *         [m3]
           zhelp(jl,jk)    = zcmm(jl,jk) * FACp9   * znumber(jl,jk)
       END DO
       END DO

       DO jc = 1,spec_number
         kcomp = species(jc)%aermlidx(jm)
         IF (kcomp == 0 .OR. kcomp == species(jc)%aernlidx(jm)) CYCLE
         DO jk = 1,klev
           DO jl = 1,kproma
             !--- 5.5.2 Total bulk aerosol mass (total PM & water):
             ![ug m-3 (air)]= [umol m-3 (air)] * [g mol-1]<=>[ug umol-1]
             zmbulki(jl,jk) = zmbulki(jl,jk) &
                            +  paerml(jl,jk,kcomp) * species(jc)%molmass
           END DO
         ENDDO
       ENDDO



     !--- 5.5.3 Aerosol mass change:

       mask(:,:,jm) = 0
       DO jk = 1,klev
       DO jl = 1,kproma
          IF(zmbulki(jl,jk) > zhelp(jl,jk) ) THEN
             mask(jl,jk,jm) = 1
             zfac (jl,jk)   = (zmbulki(jl,jk)-zhelp(jl,jk))         ! mass difference [ug m-3 (air)]
             zhelp(jl,jk)   =  zfac   (jl,jk)
          END IF
       END DO
       END DO

       !--- 5.5.3 Aerosol mass change:
       DO jc=1,spec_number
         kcomp = species(jc)%aermlidx(jm)
         jcomp = species(jc)%aermlidx(jm+1)
         IF (kcomp == 0 .OR. jcomp == 0) CYCLE
         IF (kcomp == species(jc)%aernlidx(jm) .OR. &
             jcomp == species(jc)%aernlidx(jm+1) ) CYCLE
         DO jk = 1,klev
           DO jl = 1,kproma
             IF(mask(jl,jk,jm) == 0) CYCLE
             zfac(jl,jk) = zhelp  (jl,jk)                                   & ! [ug m-3 (air)]
               * paerml (jl,jk,kcomp) * species(jc)%molmass       & ! weight with indiviual PMs ->[ug m-3]
               / zmbulki(jl,jk)                                     ! scale with total PM & water [ug m-3]
             paerml(jl,jk,jcomp) = paerml(jl,jk,jcomp) + zfac(jl,jk)/species(jc)%molmass  !  [umol m-3 (air)]
             paerml(jl,jk,kcomp) = paerml(jl,jk,kcomp) - zfac(jl,jk)/species(jc)%molmass  !  [umol m-3 (air)]
           END DO
         END DO
       END DO
     !--------------------------------------------------------------------
     END IF ! (lsize .AND. jm < CS .AND. paernl (jl,jk,jm) > cmin_aernl)
     !--------------------------------------------------------------------

#ifdef _FTRACE
     CALL FTRACE_REGION_END  ("gmxe_thermo_region_13")
     CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_14")
#endif


       !--------------------------------------------------------------------
       !--- 7. Update aerosol quantities:

       !--- Aerosol Number [N cm-3]:

     !      IF (lnumber) THEN
     DO jk = 1,klev
       DO jl = 1,kproma
         paernl(jl,jk,jm) = znumber(jl,jk) * FACm6
       END DO ! jl
     END DO ! jk
     !      END IF

     !--- 7.1 Dry Count Mean Radius [cm]:

     IF (ldryrad) THEN
       DO jk = 1,klev
         DO jl = 1,kproma
           prdry   (jl,jk,jm) = zdryrad(jl,jk)
         END DO ! jl
       END DO ! jk
     END IF

       !--- 7.2 Wet Count Mean Radius [cm]:

     IF (lwetrad) THEN
       DO jk = 1,klev
         DO jl = 1,kproma
           prwet   (jl,jk,jm) = zwetrad(jl,jk)
         END DO ! jl
       END DO ! jk
     END IF

       !--- 7.3 Dry mode density [g/cm3]:

     IF (ldrydens) THEN
       DO jk = 1,klev
         DO jl = 1,kproma
           pddry   (jl,jk,jm) = zdryden(jl,jk)
         END DO ! jl
       END DO ! jk
     END IF



     DO jk = 1,klev
       DO jl = 1,kproma
!         paerosol(jl,jk,jm,iPX)        = PX        (jl,jk,jm) ! aerosol water history  [0=solid, else=wet]
!         paerosol(jl,jk,jm,iZIONIC)    = ZIONIC    (jl,jk,jm) ! ionic strength (aq) [mol/kg]
!         paerosol(jl,jk,jm,iRHO)       = RHO       (jl,jk,jm) ! density (bulk mass) [g/m3]
!         paerosol(jl,jk,jm,iVOL)       = VOL       (jl,jk,jm) ! bulk volume [ucm3 m-3 (air)]
!
!         paerosol(jl,jk,jm,iPH)        = PH        (jl,jk,jm) ! aerosol pH [log H+]
         paerosol(jl,jk,jm,iWH2O)      = WH2O      (jl,jk,jm) ! aerosol Water (aq) [ug m-3 (air)]

       END DO ! jl
     END DO ! jk

#ifdef _FTRACE
     CALL FTRACE_REGION_END  ("gmxe_thermo_region_14")
#endif

     !--------------------------------------------------------------------
     END DO ! jm
     !--------------------------------------------------------------------

#ifdef _FTRACE
     CALL FTRACE_REGION_END  ("gmxe_thermo_region_10")
     CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_15")
#endif

     DO jk = 1,klev
       DO jl=1,kproma
         pph   (jl,jk,:) = NaN
         paopt (jl,jk) = zero
         paclc (jl,jk) = zero
         picnc (jl,jk) = zero
         pcdnc (jl,jk) = zero
         pclh2o(jl,jk) = zero
         pcih2o(jl,jk) = zero
         pcrain(jl,jk) = zero
         zmbulki(jl,jk)= zero
       END DO
     END DO

     jcomp = nwh2o(0) ! H2O vapor
     !--------------------------------------------------------------------
     IF(lacc .AND. jcomp > 0) THEN
     !--------------------------------------------------------------------

#ifdef _FTRACE
       CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_16")
#endif

     !--- 9. Determine aerosol cloud & optical properties:
     !
     ! Correct cloud cover fraction for high aerosol load (including water)
     ! to account for optical thin clouds at very low visibility (dusty air,
     ! fog and/or thin cirrus clouds).
     !
       mask2(:,:) = 0
       DO jk = 1,klev
         DO jl=1,kproma
           ! Fog formation,   if aerosol water mass > 1e-3  *  psat(jl,jk) [g m-3 (air)], which is about [1 mg  H2O/m3]
           cmin_fog(jl,jk)   = psat(jl,jk) * 1.e-5_dp  !                      [g m-3 (air)]
           ! Cloud formation, if aerosol water mass > 1e-1  *  psat(jl,jk) [g m-3 (air)], which is about [0.1g  H2O/m3]
           cmin_cloud(jl,jk) = psat(jl,jk) * 1.e-4_dp  !                      [g m-3 (air)]
           ! Cloud formation, if water vapor mass   > 0.95  *  psat(jl,jk) [g m-3 (air)], which is about [0.95g H2O/m3]
           cmin_water(jl,jk) = psat(jl,jk) * 0.95_dp   !                      [g m-3 (air)]
           ! Rain formation, if cloud water mass    > 1e+0  *  psat(jl,jk) [g m-3 (air)], which is about [10g   H2O/m3]
           cmin_rain(jl,jk)  = psat(jl,jk) * 0.99_dp   !                      [g m-3 (air)]
           ! Cloud water formation from condensing water vapor
           IF((paerml (jl,jk,jcomp) * FACm6 * mwh2o) >= cmin_water(jl,jk)) mask2(jl,jk) = 1
           pvap (jl,jk) = paerml (jl,jk,jcomp) * FACm6 * mwh2o
           pvap2(jl,jk) = paerml (jl,jk,jcomp) * FACm6 * mwh2o
           pprod(jl,jk) = 0._dp
           pminc(jl,jk) = cmin_water(jl,jk)
         END DO
       END DO
       DO jk = 1,klev
         DO jl=1,kproma
           ! Cloud water formation from condensing water vapor
           IF(mask2(jl,jk) == 0) CYCLE
           ! Cloud water mass
           ! [g m-3 (air)]=  [g m-3 (air)]  - [g m-3 (air)]
           zhelp  (jl,jk)    = (paerml(jl,jk,jcomp) * FACm6 * mwh2o) - cmin_water(jl,jk)
           pprod(jl,jk) = zhelp(jl,jk)
           ! [g m-3 (air)]=  [g m-3 (air)]  + [g m-3 (air)]
           pclh2o (jl,jk) = pclh2o  (jl,jk) + zhelp(jl,jk)
           ! Correct water vapor mass for condensed cloud water
           IF(lgh2o) &
             paerml (jl,jk,jcomp) = cmin_water(jl,jk) * FACp6 / mwh2o
           pvap2(jl,jk)         = paerml (jl,jk,jcomp) * FACm6 * mwh2o
         END DO
       END DO

       DO jm=1,nsoluble ! nmod
         !--- Total bulk aerosol water mass: [g m-3 (air)]
         kcomp = nwh2o(jm) ! H2O aerosol
         IF(kcomp > 0) THEN
           DO jk = 1,klev
             DO jl=1,kproma
               zmbulki(jl,jk) = zmbulki(jl,jk)  + paerml(jl,jk,kcomp) * FACm6 * mwh2o
             END DO
           END DO
         END IF
       END DO ! jm

#ifdef _FTRACE
       CALL FTRACE_REGION_END  ("gmxe_thermo_region_16")
       CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_17")
#endif

       mask2(:,:) = 0
       DO jk = 1,klev
         DO jl=1,kproma
           ! Cloud water formation from excess aerosol water
           IF(zmbulki(jl,jk) > cmin_cloud(jl,jk)) THEN
             mask2(jl,jk)   = 1
             ! Cloud water formation from excess aerosol water
             ! [g m-3 (air)]=  [g m-3 (air)]  - [g m-3 (air)]
             zhelp  (jl,jk) = zmbulki (jl,jk) - cmin_cloud(jl,jk)
             ! [g m-3 (air)]=  [g m-3 (air)]  + [g m-3 (air)]
             pclh2o (jl,jk) = pclh2o  (jl,jk) + zhelp(jl,jk)
             ! Correct bulk aerosol water for condensed cloud water
             zmbulki(jl,jk) = cmin_cloud(jl,jk)
           END IF
         END DO
       END DO
       pccn(:,:,:) = 0._dp
       DO jm=1,nsoluble
         ! [umol m-3 (air)]         = [umol m-3 (air)]                     , [g m-3 (air)]
         ! New total aerosol mass in this mode (wet & dry)
         zwetmass(:,:)=zero
         DO jc=1,spec_number
           kcomp = species(jc)%aermlidx(jm)
           IF(kcomp == 0 .OR. kcomp==species(jc)%aernlidx(jm) ) CYCLE
           DO jk = 1,klev
             DO jl=1,kproma
               IF(mask2(jl,jk) == 0) CYCLE
               ! [g m-3 (air)] = [g m-3 (air)]    + [umol m-3 (air)]    *   [g mol-1]
               zwetmass(jl,jk) = zwetmass(jl,jk)  + paerml(jl,jk,kcomp) * &
                 species(jc)%molmass * FACm6 ! tPM & H2O
             END DO
           END DO
         END DO ! jc
         ! Excess total aerosol mass (= new cloud mass) as fraction of total aerosol mass (tPM & H2O)
         kcomp = nwh2o(jm) ! H2O aerosol
         DO jk = 1,klev
           DO jl=1,kproma
             zhelp(jl,jk)=zero
             IF(mask2(jl,jk) == 0) CYCLE
             IF(zwetmass(jl,jk)   > cmin_cloud(jl,jk)) THEN
               zhelp(jl,jk)=1._dp - (cmin_cloud(jl,jk)/(zwetmass(jl,jk)))
             END IF
             ! Fraction of aerosol number which acts as CCN [N cm-3 (air)]
             pccn (jl,jk,jm)=zhelp(jl,jk) * paernl(jl,jk,jm)
             ! CDNC from this mode
             pcdnc(jl,jk)=pcdnc(jl,jk) + pccn(jl,jk,jm)
             ! New aerosol water mass in this mode
             ! [umol m-3 (air)]        = [umol m-3 (air)]                      , [g m-3 (air)]
!!$              IF(.NOT. LADYN .AND. kcomp > 0) &
!!$              paerml  (jl,jk,kcomp)     = MAX(zero,MIN(paerml(jl,jk,kcomp)      , cmin_cloud(jl,jk)*FACp6/mwh2o))
!!$              paerosol(jl,jk,jm,iWH2O)  = MAX(zero,MIN(paerosol(jl,jk,jm,iWH2O) , cmin_cloud(jl,jk)*FACp6/mwh2o))
           END DO
         END DO
       END DO ! jm

#ifdef _FTRACE
       CALL FTRACE_REGION_END  ("gmxe_thermo_region_17")
       CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_18")
#endif

       DO jk = 1,klev
         DO jl=1,kproma
           ! [g m-3 (air)]
           zhelp(jl,jk)=zmbulki(jl,jk)+pclh2o(jl,jk)
           mask2   (jl,jk) = 0
           IF(zhelp(jl,jk) > zero) mask2(jl,jk) = 1
           zmbulk  (jl,jk) = zero
         END DO
       END DO
       DO jm=1,nsoluble
         DO jk = 1,klev
           DO jl=1,kproma
             pph  (jl,jk,jm) = NaN
             zplus(jl,jk) = zero
             zmin (jl,jk) = zero
           END DO
         END DO
         DO jc=1,spec_number
           kcomp=species(jc)%aermlidx(jm)
           IF(kcomp == nwh2o(jm) .OR. kcomp == 0 .OR. &
             kcomp == species(jc)%aernlidx(jm) ) CYCLE ! exclude H2O aerosol
           DO jk = 1,klev
             DO jl=1,kproma
               ![mol m-3 (air)]
               zmbulk(jl,jk)    = zmbulk (jl,jk) + zionbulk(jl,jk)     * FACm6                ! aqueous PM
               ! zmbulk (jl,jk)   = zmbulk (jl,jk) + paerml(jl,jk,kcomp) * FACm6                ! total   PM
               ! total cation charge
               IF (species(jc)%charge > 0) &
                 zplus  (jl,jk)   = zplus  (jl,jk) + &
                 paerml(jl,jk,kcomp) * species(jc)%charge
               ! total anion charge
               IF (species(jc)%charge < 0) &
                 zmin  (jl,jk)   = zmin  (jl,jk) + &
                 paerml(jl,jk,kcomp) * species(jc)%charge
             END DO
           END DO
         END DO ! jc
         DO jk = 1,klev
           DO jl=1,kproma
             IF(mask2(jl,jk) == 0) CYCLE
             ! AUTODISSOCIATION CONSTANT (KW) FOR WATER H2O <==> H(aq) + OH(aq) [mol^2/kg^2]
             ztot (jl,jk)=T0/ptemp(jl,jk)
             zcoef(jl,jk)=1.0_dp+LOG(ztot(jl,jk)) - ztot(jl,jk)
             zkw  (jl,jk)=Kw*EXP(-22.52_dp*(ztot(jl,jk)-1.0_dp)+26.920_dp*zcoef(jl,jk))
             ! Aerosol water [kg m-3 (air)]
             zhelp(jl,jk)=paerosol(jl,jk,jm,iPMs)*mwh2o*1.e-9_dp
             ! Use total aerosol and cloud water in case of cloud [kg m-3 (air)]
             IF(pclh2o(jl,jk) > zero) &
               zhelp(jl,jk)=zhelp(jl,jk)+pclh2o(jl,jk)*1.e-3_dp
             ![mol^2/m-3 (air)/m-3 (air)]
             zkw  (jl,jk)=zkw(jl,jk)*prhum(jl,jk)*zhelp(jl,jk)*zhelp(jl,jk)
             ! OH- = H+ [mol/m-3 (air)]
             zkw  (jl,jk) = zkw(jl,jk)**0.5_dp
             ![mol/m-3 (air)]
             pph   (jl,jk,jm) = 0.5_dp * (-(-(zmin(jl,jk)-zplus(jl,jk))) + &
               SQRT((-(zmin(jl,jk)-zplus(jl,jk)))**2._dp - 4.0_dp*(-zkw(jl,jk))))
             phplus(jl,jk,jm) = zero
             IF(zhelp(jl,jk)  > zero) &
               ! Hydrogen ion concentration [mol/kg H2O]
               phplus(jl,jk,jm) = pph(jl,jk,jm) / zhelp(jl,jk)
             ! aerosol pH
             pph(jl,jk,jm) = NaN
             IF(pclh2o(jl,jk) > zero) &
               pph(jl,jk,jm) = 5.6_dp  ! pure wtaer & CO2
             IF(phplus(jl,jk,jm) > zero) &
               pph(jl,jk,jm) = -DLOG10(phplus(jl,jk,jm))
           END DO
         END DO
       END DO ! jm

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_thermo_region_18")
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_19")
#endif

     DO jk = 1,klev
     DO jl=1,kproma

        ! Assumed freezing point of cloud water with mixed salt solutions [K]
        ! Set to freezing point of water in case of no ions
        zice(jl,jk) = Tm
        ! Set to max. freezing depression in case of ions but no water
        IF(zmbulk(jl,jk) > zero) &
        zice(jl,jk) = Tm - Ti
        ! Calculate freezing depression in case of ions and water
        ! Tc=Cryoscopic constant for water = 1.86 [K l mol-1]
        ! 1  mol solute in 1 l water reduces the freezing point by -1.86 K.
        ! 1 l = 1000 g H2O = 55.5 mol/kg_H2O ~ 55.51 mol/l_H2O
        ! T = -1.86 K per 1  mol per 55.51 mol H2O
        ! dT [K] = PM[mol m-3 (air)]/H2O[mol m-3 (air)] * 55.51[mol/l] * 1.86 [K l mol-1]
        ! limit to -70 deg C. <=> Ti
        zhelp(jl,jk)     = (zmbulki(jl,jk)+pclh2o(jl,jk))/mwh2o
        ! Ions [mol m-3 (air)]/ (kg H2O = H2O [mol m-3 (air)] * Mw[g/mol]*1e-3[kg/mol])
        IF(zhelp (jl,jk) > zero) &
        zice (jl,jk)     = Tm - MAX(zero,MIN(zmbulk(jl,jk)/zhelp(jl,jk)*FACp3/mwh2o*Tc,Ti))
!       zice (jl,jk)     = 243.15_dp
        paerosol(jl,jk,nlowermode,iTice) = zice (jl,jk) - Tm

!!$        ! Rain formation [g m-3 (air)]
!!$        IF(pclh2o (jl,jk) > cmin_rain(jl,jk)) THEN
!!$           ! Rain water mass
!!$           zhelp  (jl,jk) = pclh2o (jl,jk) - cmin_rain(jl,jk)
!!$           pcrain (jl,jk) = pcrain (jl,jk) + zhelp(jl,jk)
!!$           ! Correct cloud water for rain water
!!$           pclh2o (jl,jk) = cmin_rain(jl,jk)
!!$           ! Correct CDNC for removal by rain water
!!$           IF(zhelp(jl,jk) > zero) &
!!$           zhelp(jl,jk)=(1._dp-cmin_rain(jl,jk)/zhelp(jl,jk))
!!$           ! Fraction of CDNC removed by rain [N cm-3 (air)]
!!$           pcdnc(jl,jk)=pcdnc(jl,jk)-zhelp(jl,jk)*pcdnc(jl,jk)
!!$        END IF
!!$        ! Evaopration of rain water [g m-3 (air)]
!!$        IF(pcrain (jl,jk) < cmin_rain(jl,jk) .AND.  pcrain (jl,jk) > zero) THEN
!!$           pclh2o (jl,jk) = pclh2o (jl,jk) + pcrain (jl,jk)
!!$           pcrain (jl,jk) = zero
!!$        END IF

!!$        ! Freezing of cloud droplets [g m-3 (air)]
!!$        IF(ptemp  (jl,jk) < zice  (jl,jk) .AND. pclh2o(jl,jk) > zero) THEN
!!$           pcih2o (jl,jk) = pcih2o(jl,jk) +  pclh2o(jl,jk)
!!$           pclh2o (jl,jk) = zero
!!$           ! ICNC from this mode
!!$           picnc  (jl,jk) = picnc (jl,jk) +  pcdnc(jl,jk)
!!$           pcdnc  (jl,jk) = zero
!!$        ! Melting of cloud droplets [g m-3 (air)]
!!$        ELSE IF(ptemp(jl,jk) > zice(jl,jk).AND. pcih2o(jl,jk) > zero) THEN
!!$           pclh2o (jl,jk) = pclh2o(jl,jk) +  pcih2o(jl,jk)
!!$           pcih2o (jl,jk) = zero
!!$           pcdnc  (jl,jk) = pcdnc (jl,jk) +  picnc(jl,jk)
!!$           picnc  (jl,jk) = zero
!!$           ! Evaporation of cloud droplets [g m-3 (air)]
!!$           IF(pclh2o(jl,jk)    < cmin_cloud(jl,jk) .AND. pclh2o(jl,jk) > zero) THEN
!!$              paernl(jl,jk,CS) = paernl(jl,jk,CS) + pcdnc(jl,jk)
!!$              pcdnc (jl,jk)    = zero
!!$              IF(lgh2o) &
!!$              paerml (jl,jk,jcomp) = paerml (jl,jk,jcomp) + pclh2o(jl,jk) * FACp6 / mwh2o
!!$              pclh2o(jl,jk)        = zero
!!$           END IF
!!$        END IF
     END DO ! jl
     END DO ! jk

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_thermo_region_19")
CALL FTRACE_REGION_BEGIN("gmxe_thermo_region_20")
#endif

     DO jk = 1,klev
     DO jl=1,kproma
        ! Cloud cover - fog and high aerosol load:
        ! Fog type: 0=no fog, 0-50=brown fog (dense smog), 50-100=fog [%]
        paopt (jl,jk) = zero
        paclc (jl,jk) = zero
        ! Fog cloud, i.e. high aerosol water load
        IF(zmbulki(jl,jk) > cmin_fog(jl,jk)) then
        paopt  (jl,jk) = 100._dp - cmin_fog(jl,jk)/zmbulki(jl,jk) * 50._dp      ! 50-100 [%]
        END IF
        ! Brown fog (dense smog), i.e. low aerosol water but very high dry aerosol load
        IF(zmbulk (jl,jk) > cmin_fog(jl,jk)) &
        paopt  (jl,jk) = 49.0_dp - cmin_fog(jl,jk)/zmbulk(jl,jk)  * 49._dp      !  0-50 [%]
        ! Cloud cover (ice and warm cloud, rain/snow, fog and dense smog):
        zhelp(jl,jk)      = zmbulk(jl,jk)+zmbulki(jl,jk)+pclh2o(jl,jk)+pcih2o(jl,jk)+pcrain(jl,jk)
        IF(zhelp (jl,jk) >= cmin_cloud(jl,jk)) &
        paclc  (jl,jk) = 100._dp - cmin_cloud(jl,jk) / zhelp(jl,jk) * 100._dp  ! [%]
     END DO ! jl
     END DO ! jk

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_thermo_region_20")
#endif

     !--------------------------------------------------------------------
     END IF ! lacc
     !--------------------------------------------------------------------

#ifdef _FTRACE
CALL FTRACE_REGION_END  ("gmxe_thermo_region_15")
#endif

  !--------------------------------------------------------------------
  IF(l_io .AND. .NOT. entered ) &
  CALL end_message(modstr,'calculate aerosol composition',substr)
  entered = .TRUE.
  !--------------------------------------------------------------------

END SUBROUTINE gmxe_thermo
!
!------------------------------------------------------------------------------
!
SUBROUTINE gmxe_dgas(kproma, klev,   pgasi,  paernl,         &
                     ptemp,  ppress, prwet,  pgas,   ztmst,  &
                     mod1,   mod2 )
  !
  !**** *gmxe_dgas*  calculates the transfer of mass due to
  !                 condensation of gas phase species on soluble and
  !                 insoluble modes.
  !
  !    Authors:
  !    -----------
  !    J. Wilson, E. Vignati, JRC/EI (original source)                05/2000
  !    P. Stier, MPI                 (f90-version, changes, comments)    2001
  !
  !    Modifications:
  !    --------------
  !    H. Tost, UMZ generalised species structure
  !    S. Metzger, MPI-CHEM (extended and generalized for use with EQSAM3), June 2008
  !
  !    Purpose:
  !    -----------
  !    This routine calculates the changes in aerosol mass and
  !    gas phase species due to condensation.
  !
  !**  Interface:
  !    -----------
  !    *gmxe_dgas* is called from *gmxe*
  !
  !    Method:
  !    -----------------
  !    The transfer of gas phase species to the particles is based on
  !    Fuchs (1959). Soluble and insoluble particles are distinguished
  !    by a different accomodation coefficient "caccgas" defined
  !    in messy_gmxe_mem.
  !    (Currently 1.0 for soluble and 0.3 for insoluble modes).
  !
  !    Externals
  !    -----------
  !    none
  !
  !    References:
  !    -----------
  !    Fuchs, N.A. (1959). Evaporation and droplet growth in gaseous media;
  !       Pergamon, New York, pp72.
  !
  !--- Parameter list:
  !
  ! pgasi     = conc. of gas phase species i [molec. cm-3]
  ! pgas      = mass of condensable gas phase species on any mode [molec. cm-3]
  ! mod1,mod2 = mode boundaries for loops
  !
  !--- Local Variables:
  !
  ! zde       = molecular diffusion []
  ! zvelb     = velocity []
  ! zcondo    = condensation coefficient []
  ! zc2(nmod) = flux of gas phase species condensing on the respective mode
  !             per gas phase concentration []
  ! zcondo    = total flux of condensing gas phase species
  !             per gas phase concentration []
  ! zfcond    = total mass of condensing gases for one timestep []

  ! mz_ht_20081210+
  ! treating only the modes that have not already been treated
  ! within gmxe_thermo



  IMPLICIT NONE

  INTEGER :: kproma, klev, kcomp, jk, jl, jgas, jmod, jcomp
  INTEGER, INTENT(IN) :: mod1, mod2
  REAL(dp):: ptemp(kproma,klev),       ppress(kproma,klev),              &
             pgasi(kproma,klev,ngas),  pgas  (kproma,klev,ngas,nmod),  &
             zxcond(kproma,klev,ngas)
  !
  REAL(dp):: paernl(kproma,klev,nmod), prwet(kproma,klev,nmod)
  !
  ! Local variables:

  REAL(dp):: zfcond,     zftot,      zpbyone,    zde2,                   &
             zvelb,      zxibc,      zrwet,      zf1,                    &
             ztmst,      zqtmst

  REAL(dp):: zcondo(kproma,klev,ngas),zc2(kproma,klev,ngas,nmod)

  zc2   (:,:,:,:) = zero
  zcondo(:,:,:)   = zero
  zqtmst          = 1._dp/ztmst
  zxcond(:,:,:)   = 0._dp

  !--- 1) Calculate condensation rate for cm diameter aerosols: ---
  !
  DO jmod = mod1,mod2
    DO jgas=1,ngas      ! loop over all gas phase species
      DO jk=1,klev
        DO jl=1,kproma
          IF (prwet(jl,jk,jmod) > REALZERO) THEN

            !--- Diffusion coefficient (Reference???):

            zpbyone=1000.0_dp / (ppress(jl,jk)/100.0_dp) ! hPa?

            zde2=0.073_dp * zpbyone * (ptemp(jl,jk) / 298.15_dp)**1.5_dp

            !--- Mean molecule velocity (Moore, 1962 (S+P equ. 8.2)):

            zvelb=SQRT(8.0_dp * rerg * ptemp(jl,jk) / &
              pi / td%gas(jgas)%molmass)

            !--- ???Fuchs???

            zxibc=8.0_dp * zde2 / pi / zvelb

            !
            ! Use count median radius:

            zrwet=prwet(jl,jk,jmod)

            !--- Distance from particle up to which the kinetic regime applies:

            zf1=( (zrwet + zxibc)**3.0_dp - &
                (zrwet**2.0_dp + zxibc**2.0_dp)**1.5_dp ) / &
                (3.0_dp * zrwet * zxibc) - zrwet

            !--- Diffusive flux to single particle surface:
            !    (Elisabetta's thesis: fraction in equ. 2.26)

            zc2(jl,jk,jgas,jmod)=(4.0_dp * pi * zde2 * zrwet ) /        &
                                 ((4.0_dp * zde2) /                     &
                                 (zvelb * zrwet * td%gas(jgas)%caccgas(jmod)) + &
                                 (zrwet/(zrwet+zf1)) )

            !--- Total diffusive flux to all particles in the respective mode:
            !    (per concentration of gas phase species)

            zc2(jl,jk,jgas,jmod) = zc2(jl,jk,jgas,jmod) * paernl(jl,jk,jmod)

            !--- Total diffusive flux to all particles in all modes:
            !    (per concentration of gas phase species)

            zcondo(jl,jk,jgas)=zcondo(jl,jk,jgas) + zc2(jl,jk,jgas,jmod)

          END IF
        END DO
      END DO
    END DO ! jgas
  END DO ! jmod
  !
  !--- 2) Calculation of the new aerosol masses and of the ---------
  !       mass of gas phase species condensing on the respective modes:
  !
  DO jmod = mod1,mod2
    DO jgas=1,ngas ! loop over all gas phase species

      DO jk=1,klev
        DO jl=1,kproma
          IF(zcondo(jl,jk,jgas) > REALZERO .AND. &
             pgasi(jl,jk,jgas) > REALZERO) THEN

            !--- Total diffusive flux to all particles in all modes:

            zfcond=zcondo(jl,jk,jgas)*pgasi(jl,jk,jgas)

            !--- Total mass of species i condensing during 1 timestep:

            zftot=zfcond*ztmst

            !--- Limit condensing gases to
            !    fmax[%] x (available gas-phase conc.) :

            zfcond=MIN(zftot,(pgasi(jl,jk,jgas)*fmax))
            zxcond(jl,jk,jgas) = zfcond

            !--- Add mass to respective aerosol species:
            !    zc2(:,:,:,jmod)/zcondo = fraction of total gas phase flux
            !                           that condenses on mode jmod
            !    => (   "    )*zfcond = mass of gases condensing on
            !                           the respective mode

            !--- Mass of gas phase sepcies condensing on the insoluble modes:
            !    (Transfer from insoluble to soluble modes
            !     calculated in gmxe_concoag.)

            pgas(jl,jk,jgas,jmod)=zc2(jl,jk,jgas,jmod)/zcondo(jl,jk,jgas)*zfcond


          END IF
        END DO
      END DO
    END DO ! jgas
  END DO ! jmod

  do jgas=1,ngas
    do jk=1,klev
      do jl=1,kproma
    !--- Remaining gas phase conc.:
        pgasi(jl,jk,jgas) = pgasi(jl,jk,jgas)-zxcond(jl,jk,jgas)
      enddo
    enddo
  enddo


END SUBROUTINE gmxe_dgas
!------------------------------------------------------------------------------!
SUBROUTINE gmxe_dnum(kproma, klev,                  &
                     pso4g,  paerml, paernl, ptemp, &
                     ppress, prhum,  prwet,  pddry, &
                     pgas,   pncrit, ztmst          )
  !
  !  *gmxe_dnum*  changes gas phase conc., aerosol numbers and masses
  !               due to nucleation and coagulation
  !
  !  Authors:
  !  ---------
  !  J. Wilson  and E. Vignati, JRC/EI (original source)                09/2000
  !  P. Stier, MPI                     (f90-version, changes, comments)    2001
  !
  !  Modifications:
  !  --------------
  !  H. Tost, UMZ generalised species structure
  !  S. Metzger, MPI-CHEM (extended and generalized for use with EQSAM3), June 2008
  !
  !  Version:
  !  ---------
  !  This version is equivalent to the version dnum2 of the gmxe boxmodel.
  !
  !  Purpose
  !  ---------
  !  This routine calculates new gas phase conc. and aerosol
  !  numbers and masses after timestep ztmst.
  !
  !  Interface
  !  -----------
  !  *gmxe_dnum* is called from *gmxe*
  !
  !  Externals
  !  -----------
  !  *gmxe_coaset* calculates the coagulation kernels
  !  *gmxe_nuck*   calculates the integral mass nucleated sulfate and
  !                the number of nucleated particles over one timestep
  !  *gmxe_delcoa* integrates equations for the changes in aerosol numbers
  !                dn/dt=c -an2 -bn over one timestep and calculates the
  !                corresponding changes in aerosol masses
  !  *gmxe_concoag*calculates particle numbers and mass moved from the
  !                insoluble to the mixed modes due to coagulation with
  !                smaller mixed particles and condensation of H2SO4.
  !
  !  Warning:
  !  --------
  !  For optimization purposes currently only "physically reasonable" elements of the
  !  coagulation kernel zcom are calculated in gmxe_concoag.
  !  These elements are specified
  !  in the matrix locoagmask in messy_gmxe.
  !  Check carefully and adapt locoagmask
  !  accordingly  before changes in the code below.
  !
  !--- Parameters:
  !
  !  pso4g      = mass of gas phase sulfate [molec. cm-3] for nucleation
  !  paerml     = total aerosol mass for each compound
  !  paernl     = aerosol number for each mode [cm-3]
  !  ptemp      = atmospheric temperature at time t+1 [K]
  !  ppress     = atmospheric pressure at time t+1 [Pa]
  !  prhum      = atmospheric relative humidity [0-1]
  !  prwet      = aerosol wet radius  [cm]
  !  prdry      = aerosol dry radius  [cm]
  !  pddry      = aerosol dry density [g cm-3]
  !  pncrit     = number of molecules in the critical cluster of nucleation [1]
  !
  !--- Local variables:
  !
  ! zcom            = general coagulation coefficient []
  !                   (calculated in gmxe_coaset)
  ! za              = effectively used coagulation coefficient for
  !                   unimodal coagulation []
  ! zb              = effectively used coagulation coefficient for
  !                   inter-modal coagulation []
  ! zbfractx(:,:,y) = fraction of the total number of particles removed by
  !                   coagulation from mode x (finally calculated in gmxe_delcoa)
  !                   that is moved to mode y / y+1 (modes 5,6,7 and mode 1,2 resp.) [1]
  ! za4delt(:,:,:)  = change in H2SO4 mass of the respective mode over one timstep
  !                   due to:
  !                      - nucleation of H2SO4 (calculated in gmxe_nuck)
  !                      - coagulation (calculated in gmxe_concoag)
  ! zanew           = number of nucleated particles
  !                   zanew=za4delt/critn i.e. mass of formed sulfate
  !                   divided by an assumed mass of a nucleus. []
  !                   (calculated in gmxe_nuck)
  ! zxxavy          = average mass of species xx in mode y []
  !                   where xx is ss, du, bc, oc, or a4 for sulfate
  !                   [molecules for sulfate and ug for others]

  !--- Parameters:

  IMPLICIT NONE

  INTEGER :: kproma, klev

  REAL(dp):: pso4g(kproma,klev),          ptemp(kproma,klev),               &
             ppress(kproma,klev),         prhum(kproma,klev),               &
             pncrit(kproma,klev),         pgas(kproma,klev,ngas,nmod)

  REAL(dp):: paerml(kproma,klev,0:naertot), paernl(kproma,klev,nmod),       &
             prwet(kproma,klev,nmod),       pddry(kproma,klev,nmod)

  ! Local Variables:

  INTEGER :: jl, jk, ji, jm

  REAL(dp):: ztmst

  REAL(dp):: zanew(kproma,klev)

  REAL(dp):: za(kproma,klev,nmod),        zb(kproma,klev,nmod),              &
             za4delt(kproma,klev,0:naertot)

  REAL(dp):: zbfract1(kproma,klev,nmod-1),zbfract2(kproma,klev,nmod-1),      &
             zbfract5(kproma,klev,3),     zbfract6(kproma,klev,2),           &
             zbfract7(kproma,klev,2)


  REAL(dp):: zcom(kproma,klev,nmod,nmod)

  REAL(dp) :: xbfract(kproma,klev,nmod,nmod)

  xbfract(:,:,:,:) = 0._dp

  !--- 0) Initialisations: ----------------------------------------------

  za4delt(:,:,:) = zero  ! Mode 1 changed by gmxe_nuck if lnucl=TRUE .AND. nnucl==1.
                         ! Has to be initialized for the other modes.
  zanew(:,:)     = zero  ! Changed by gmxe_nuck if lnucl=TRUE.AND. nnucl==1.

  !--- 1) Calculate  coagulation coefficients: --------------------------
  !
!CDIR NOIEXPAND
  IF (lcoag) CALL gmxe_coaset(kproma, klev,  paernl, ptemp,  &
                              ppress, prwet, pddry,  zcom    )

  !
  !--- 2) Calculate nucleation rate, number of nucleated particles ------
  !       and changes in gas phase sulfate mass.
  !

!CDIR NOIEXPAND

  IF (lnucl.AND.nnucl>=1) CALL gmxe_nuck(kproma, klev,  pso4g,          &
                                         ptemp,  prhum, zanew, za4delt, &
                                         pncrit, ztmst                  )
  !
  !--- 3) Assign coagulation coefficients (unimodal coagulation)---------
  !       and the normalised fraction of particle numbers moved
  !       between the modes (inter-modal coagulation):
  !
  !       The general equation for dn/dt for each mode is:
  !
  !       dn/dt=-za*n^2 - zb*n + zc
  !
  !       where za=unimodal coagulation coefficient (zcom(mod))
  !             zb=inter-modal coagulation with higher modes
  !                (zcom(mod) * n(jmod+1))
  !             zc=particle formation rate
  !                (=zanew/ztmst if jmod=1, zero for higher modes)
  !
  !             zb is zero when n(jmod+1)=zero, or jmod=naertot
  !
  za(:,:,:)       = zero
  zb(:,:,:)       = zero
  zbfract1(:,:,:) = zero
  zbfract2(:,:,:) = zero
  zbfract5(:,:,:) = zero
  zbfract6(:,:,:) = zero
  zbfract7(:,:,:) = zero

  IF (lcoag) THEN

  ! todo : generalize coagulation (jmod/kmod)

    DO jm=1,nmod
      IF (locoagmask(jm,jm)) THEN
        DO jk=1,klev
          DO jl=1,kproma

          !---  3.1) Unimodal coagulation coefficients:
          !          only if allowed in the coagulation matrix locoagmask

            za(jl,jk,jm) = zcom(jl,jk,jm,jm) / 2._dp
          END DO
        END DO
      ENDIF
    END DO
    DO jm=1,nsoluble
      DO ji=1,nmod-1
        DO jk=1,klev
          DO jl=1,kproma
            ! number of particles that are moved from mode jm
            ! to mode ji
            ! Since not all coagulation moves particles away
            ! from mode jm (i.e. coagulation of hydrophilic
            ! with hydrophobic particles) a mode transfer matrix
            ! xmov(nmod,nmod) is multiplied with the coagulation
            ! matrix and the numbers
            xbfract(jl,jk,jm,ji) = zcom(jl,jk,ji+1,jm) * xmov(jm,ji+1) * &
                                   paernl(jl,jk,ji+1)
            zb(jl,jk,jm) = zb(jl,jk,jm) + xbfract(jl,jk,jm,ji)
          ENDDO
        ENDDO
      ENDDO
    END DO
    DO jm=nsoluble+1,nmod
      DO ji=1,nmod-1
        DO jk=1,klev
          DO jl=1,kproma
            ! number of particles that are moved from mode jm
            ! to mode ji
            ! Since not all coagulation moves particles away
            ! from mode jm (i.e. coagulation of hydrophilic
            ! with hydrophobic particles) a mode transfer matrix
            ! xmov(nmod,nmod) is multiplied with the coagulation
            ! matrix and the numbers
            xbfract(jl,jk,jm,ji) = zcom(jl,jk,ji,jm) * xmov(jm,ji) * &
                                   paernl(jl,jk,ji)
            zb(jl,jk,jm) = zb(jl,jk,jm) + xbfract(jl,jk,jm,ji)
          ENDDO
        ENDDO
      ENDDO
    END DO
    DO jm=1,nmod
      DO ji=1,nmod
        DO jk=1,klev
          DO jl=1,kproma
            IF (zb(jl,jk,jm) > zero) THEN
              xbfract(jl,jk,jm,ji) = xbfract(jl,jk,jm,ji) / zb(jl,jk,jm)
            END IF
          END DO
        END DO
      END DO
    END DO

     !
   END IF !(lcoag)
  !
  !
!CDIR NOIEXPAND

  CALL gmxe_delcoa(kproma,   klev,    paerml,            &
                   paernl,   prwet,  za4delt,  zanew,    &
                   za,       zb,     zbfract1, zbfract2, &
                   zbfract5, pgas,   ztmst               )

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

CONTAINS

SUBROUTINE gmxe_coaset(kproma, klev,  paernl, ptemp, &
                       ppress, prwet, pddry,  pcom   )
  !
  ! *gmxe_coaset*  calculates the coagulation kernels between the modes
  !
  ! Authors:
  ! ---------
  ! J. Wilson  and E. Vignati, JRC/EI (original source)                09/2000
  ! P. Stier, MPI                     (f90-version, changes, comments)    2001
  !
  ! Modifications:
  ! --------------
  ! Philip Stier, MPI                             2001
  !
  ! Purpose
  ! ---------
  ! This routine calculates the coaglation kernels between particles
  ! with the count median radii of the three modes.
  ! Coagulation allowed between the following modes:
  ! soluble modes:   1+1=1, 2+2=2, 1+2=2, 1+3=3, 1+4=4, 2+3=3, 2+4=4
  ! insoluble modes: 2i+2i=2i
  ! mixed modes:     1+2i=2, 1+4i=4, 2+2i=2, 2+4i=4, 3+2i=3, 4+2i=4
  !
  ! Interface:
  ! -----------
  ! *gmxe_coaset* is called from *gmxe_dnum*
  !
  ! Externals:
  ! -----------
  ! none
  !
  ! Reference:
  ! -----------
  ! The calculations are based on:
  ! Fuchs, N.A. (1964). The Mechanics of Aerosols. Pergamon Press. Oxford.
  ! (Chapter VII, Section 49)
  !
  !  Warning:
  !  --------
  !  For optimization purposes currently only "physically reasonable"
  !  elements of the coagulation kernel pcom are calculated in gmxe_concoag.
  !  These elements are specified in the matrix locoagmask in messy_gmxe.
  !  Check carefully and adapt locoagmask accordingly  before changes
  !  in the code below.
  !
  !--- Parameter list:
  !
  !  paernl            = aerosol number for each mode [cm-3]
  !  ptemp             = atmospheric temperature at time t+1 [K]
  !  ppress            = atmospheric pressure at time t+1 [Pa]
  !  prwet             = aerosol wet radius  [cm]
  !  prdry             = aerosol dry radius  [cm]
  !  pddry             = aerosol dry density [g cm-3]
  !  pcom(:,:,jm1,jm2) = Coagulation coefficient for modes jm1 and jm2 []
  !
  !--- List of local variables:
  !
  ! zwlc              = mean free pathlength []
  ! zairvisc          = air viscosity []
  ! zrpav             = average radius of the two interacting modes
  ! zpvx              = volume of the xth interacting mode
  ! zpmx              = mass of the xth interacting mode
  !
  ! zrknudx           = knudsen number of the xth interacting mode
  ! zpd2x             = particle diffusion of the xth interacting mode
  !
  !--- Parameters:
  !
  IMPLICIT NONE

  INTEGER :: kproma, klev

  REAL(dp):: ptemp(kproma,klev),         ppress(kproma,klev)

  REAL(dp):: prwet(kproma,klev,nmod),    pddry(kproma,klev,nmod),                 &
             paernl(kproma,klev,nmod)

  REAL(dp):: pcom(kproma,klev,nmod,nmod)

  !--- Local variables:

  INTEGER :: jm2, jm1, jl, jk
  !
  REAL(dp):: zpbyone,     zwlc,        zairvisc,    zbtketc,   &
             zrpvm1,      zrpvm2,      zrpav,       zpv1,      &
             zpm1,        zcv21,       zpv2,        zpm2,      &
             zcv22,       zcv2av,      zrknud1,     zrknud2,   &
             ze1,         ze2,         zpd21,       zpd22,     &
             zpdav,       zxd1,        zh21,        zxd2,      &
             zh22,        zh2av,       zcoc,        zhu1,      &
             zhu2,        zhu

  !--- 1) Calculation of the coagulation coefficient: ---------------------------
  !
  pcom(:,:,:,:)=zero
  DO jm2=1,nmod
     DO jm1=jm2,nmod

        IF (locoagmask(jm1,jm2)) THEN

           DO jk=1,klev
              DO jl=1,kproma

                 IF (paernl(jl,jk,jm1) > cmin_aernl .AND.  &
                     paernl(jl,jk,jm2) > cmin_aernl .AND.  &
                     prwet (jl,jk,jm1) > FACm10     .AND.  &
                     prwet (jl,jk,jm2) > FACm10     .AND.  &
                     pddry (jl,jk,jm1) > FACm10     .AND.  &
                     pddry (jl,jk,jm2) > FACm10        ) THEN

                    !--- 1.1) Calculate ambient properties:

                    !--- Mean free pathlength ? (from Knudsen Number below):
                    !    Parametrisation?
                    zpbyone=1000.0_dp / (ppress(jl,jk)/100.0_dp) ! hPa?
                    zwlc=6.6e-6_dp * ptemp(jl,jk) / 293.15_dp * zpbyone

                    !--- Viscosity:
                    zairvisc=1.827e-4_dp * (ptemp(jl,jk) / 291.15_dp)**0.74_dp

                    !---
                    zbtketc=bk * ptemp(jl,jk) / 6.0_dp / pi / zairvisc

                    !--- Count median radii of the respective modes:
                    zrpvm1=prwet(jl,jk,jm1)
                    zrpvm2=prwet(jl,jk,jm2)

                    !--- Average radius of the modes:
                    zrpav=(zrpvm1 + zrpvm2) / 2.0_dp

                    !--- Volume and mass of mode jm1:
                    zpv1=4.0_dp / 3.0_dp * pi * zrpvm1**3.0_dp
                    zpm1=zpv1 * pddry(jl,jk,jm1)

                    !--- Squared mean particle velocity of mode jm1:
                    zcv21=8.0_dp * bk * ptemp(jl,jk) / (pi * zpm1)

                    !--- Volume and mass of particles in mode jm2:
                    zpv2=4.0_dp / 3.0_dp * pi * zrpvm2**3.0_dp
                    zpm2=zpv2 * pddry(jl,jk,jm2)

                    !--- Squared mean particle velocity of mode jm2:
                    zcv22=8.0_dp * bk * ptemp(jl,jk) / (pi * zpm2)

                    !--- Fuchs: G_r (below Eq. 49.27):
                    zcv2av=SQRT(zcv21 + zcv22)

                    !--- Knudsen numbers of the modes:
                    zrknud1=zwlc/zrpvm1/2.0_dp
                    zrknud2=zwlc/zrpvm2/2.0_dp

                    !--- Diffusivities of the respective modes:
                    ze1=EXP(-0.43_dp/zrknud1)
                    ze2=EXP(-0.43_dp/zrknud2)
                    ! the bracket is (?) the Cunningham Slip correction
                    zpd21=zbtketc * (1.0_dp + zrknud1*2.492_dp + &
                      zrknud1*0.84_dp*ze1) / zrpvm1
                    zpd22=zbtketc * (1.0_dp + zrknud2*2.492_dp + &
                      zrknud2*0.84_dp*ze2) / zrpvm2
!                      zrknud1*0.84_dp*ze2) / zrpvm2  ! corrected due to unlikely combination; however not sure !!!!

                    !--- Average diffusivity of the modes:
                    zpdav=(zpd21 + zpd22) / 2.0_dp

                    !--- Average mean free path of particles in jm1:
                    zxd1=8.0_dp * zpd21 / (pi*SQRT(zcv21))

                    !--- Mean distance from surface after mean
                    !    free path (Eq. 49.13):
                    zh21=(((2.0_dp*zrpvm1 + zxd1)**3.0_dp -                   &
                          (4.0_dp*zrpvm1*zrpvm1 + zxd1*zxd1)**1.5_dp) /       &
                          (6.0_dp*zrpvm1*zxd1) - 2._dp*zrpvm1   ) * sqrt2

                    !--- Average mean free path of particles in jm2:
                    zxd2=8.0_dp * zpd22 / (pi*SQRT(zcv22))

                    !--- Mean distance from surface after mean
                    !    free path (Eq. 49.13):

                    zh22=(((2.0_dp*zrpvm2 + zxd2)**3 -                        &
                          (4.0_dp*zrpvm2*zrpvm2 + zxd2*zxd2)**1.5_dp) /       &
                          (6.0_dp*zrpvm2*zxd2) - 2._dp*zrpvm2        ) * sqrt2

                    !--- Fuchs: delta_r
                    zh2av=SQRT(zh21*zh21 + zh22*zh22) / sqrt2

                    !- 1.2) Calculation of the coagulation coefficient pcom
                    !       (Eq. 49.26):
                    !       Factor 16 instead factor 8 as in Fuchs as his
                    !       formulation applies for the inter-modal coagulation.
                    !       This is taken into account in the assignment of the
                    !       inter-modal coagulation coefficient.

                    zcoc=16.0_dp * pi * zpdav * zrpav

                    !--- Calculation of beta=1/zhu (Eq. 49.27):
                    zhu1=4.0_dp * zpdav / (zcv2av * zrpav)
                    zhu2=zrpav / (zrpav + zh2av /2.0_dp)
                    zhu=zhu1 +  zhu2

                    !--- Coagulation coefficient following (Eq.49.26):
                    pcom(jl,jk,jm1,jm2)=zcoc / zhu
                 END IF

                 !--- 2) Mirror the coagulation matrix (symmetric): -------
                 pcom(jl,jk,jm2,jm1)=pcom(jl,jk,jm1,jm2)
              END DO
           END DO

         END IF ! locoagmask

       END DO
     END DO


END SUBROUTINE gmxe_coaset

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
SUBROUTINE gmxe_nuck(kproma, klev,  pso4g,           &
                     ptemp,  prhum, panew, pa4delt,  &
                     pncrit, ztmst                   )
  !
  !  Authors:
  !  --------
  !  J. Wilson, E. Vignati, JRC/EI, original source                 09/2000
  !  P. Stier, MPI                  f90-version,
  !                                 changes,
  !                                 comments,
  !                                 modularisation and
  !                                 implementation of
  !                                 Vehkamaeki (2002)             2001-2003
  !
  !  Purpose:
  !  --------
  !  This routine calls the routines for the computation of the
  !  nucleation rate znucrate [molec. cm-3 s-1] and, for
  !  Vehkamaeki (2002), of the number of molecules in the critical
  !  cluster from a given gas phase H2SO4 concentration
  !  pso4g [molec. cm-3]. It also calculates the integrated change of
  !  H2SO4 gas phase mass over one timestep due to nucleation
  !  pa4delt(:,:,nso42m(NS)) [molec. cm-3] as well as the number of nucleated
  !  particles panew [1] during the timestep. Whilst this is done
  !  analytically for the old Kulmala (1998) parameterization, it has
  !  to be done numerically for the new Vehkamaeki (2002) scheme.
  !
  ! The old parameterization of Kulmala (1998), that apparently contains
  ! some errors is kept for consistency. It is recommended to use the
  ! new scheme of Vehkamaeki et al. (2002). However, this parameterization
  ! is no longer analytically integrable, the effect of this is to be
  ! scrutinized.
  ! The modularized version of the Kulmala parameterization has been tested
  ! to give bit-identical results with the old version.
  !
  ! References:
  ! -----------
  ! Vehkamaeki et al. (2002), An improved parameterization for sulfuric
  !    acid/water nucleation rates for tropospheric and stratospheric
  !    conditions, J. Geophys. Res, 107, D22, 4622
  ! Kulmala et al. (1998), Parameterizations for sulfuric acid/water
  !    nucleation rates. J. Geophys. Res., 103, No D7, 8301-8307
  ! Vignatti, E. (1999), Modelling Interactions between Aerosols and
  !    Gaseous Compounds in the Polluted Marine Atmosphere. PhD-Thesis,
  !    RISO National Laborartory Copenhagen, Riso-R-1163(EN)
  !
  !  Interface:
  !  ----------
  !  *gmxe_nuck* is called from *gmxe_dnum*
  !
  !  Method:
  !  -------
  !
  !  1) Kulmala et al. (1998):
  !
  !  For the Kulmala et al. (1998) scheme the formula for binary
  !  nucleation is rewritten to take the form
  !  znucrate = exp[zalpha+ln(pso4g)*beta].
  !  Equation numbers are taken from Kulmala et al. (1998).
  !  After the calculation of the nucleation rate znucrate, it is
  !  integrated in 2) analytically over one timestep, i.e.:
  !
  !  Integration of:
  !
  !  znucrate=d(critn*znav)/dt=exp[zalpha + zbeta*ln(znav)]
  !
  !  gives znav(t0+dt, znav(t0)) where znav(t0) is pso4g.
  !  znav is temporarily stored in zso4g_new and confined between
  !  0 and pso4g.
  !  The number of nucleated particles is then calculated by
  !  assuming a fixed critical mass critn of newly nucleated
  !  particles and dividing the total mass of nucleated sulfate
  !  by this critical mass of one nucleated particle.
  !
  !
  !  2) Vehkamaeki et al. (2002):
  !
  !  An analytical integration of the nucleation equation is not
  !  possible, therefor the nucleation rate is simply multiplied
  !  by dt, implying a fixed gas-phase H2SO4 concentration over
  !  the timestep.
  !  The number of nucleated particles is then calculated by
  !  taking the calculated critical mass critn of newly nucleated
  !  particles and dividing the total mass of nucleated sulfate
  !  by this critical mass of one nucleated particle.
  !
  !  Externals:
  !  ----------
  !  gnucl_kulmala
  !  gnucl_vehkamaeki
  !
  !  References:
  !  -----------
  !  Vehkamaeki et al. (2002), An improved parameterization for sulfuric
  !     acid/water nucleation rates for tropospheric and stratospheric
  !     conditions, J. Geophys. Res, 107, D22, 4622
  !  Kulmala et al. (1998), Parameterizations for sulfuric acid/water
  !     nucleation rates. J. Geophys. Res., 103, No D7, 8301-8307
  !  Vignatti, E. (1999), Modelling Interactions between Aerosols and
  !     Gaseous Compounds in the Polluted Marine Atmosphere. PhD-Thesis,
  !     RISO National Laborartory Copenhagen, Riso-R-1163(EN)
  !
  !
  !--- Parameters:
  !
  ! pso4g          = mass of gas phase sulfate [molec. cm-3]
  ! ptemp          = atmospheric temperature at time t+1 [K]
  ! prhum          = atmospheric relative humidity [0-1]
  ! pa4delt(:,:,nso42m(NS)) = mass of H2SO4 added to the nucleation mode due
  !                  to nucleation of H2SO4 over ztmst.
  !                  Equilvalent to the integral of H2SO4 gas loss
  !                  due to nucleation over timestep ztmst. [molec. cm-3]
  ! panew          = number of nucleated particles over timestep ztmst
  !                  panew=pa4delt/critn i.e. mass of formed sulfate
  !                  divided by an assumed mass of a nucleus. [cm-3]
  ! pncrit         = number of molecules in the critical cluster [1]
  !
  !--- Local variables:
  !
  ! zso4g_new        = temporay storage of gas phase sulfate [molec. cm-3]
  !
  ! See comments!

  IMPLICIT NONE

  INTEGER :: kproma, klev

  REAL(dp):: pso4g(kproma,klev),        ptemp(kproma,klev),     &
             prhum(kproma,klev),        panew(kproma,klev)

  REAL(dp):: pa4delt(kproma,klev,0:naertot), pncrit(kproma,klev)

  ! Local variables:

  INTEGER :: jk,          jl,           NS

  REAL(dp):: ztmst,       zqtmst,       zf1!,          zeps

  REAL(dp):: znucrate(kproma,klev),  & ! nucleation rate [m-3 s-1]
             zso4g_new(kproma,klev), & ! new gas phase sulfate concentration [molec. cm-3]
             zalpha(kproma,klev),    & ! auxiliary coefficient for the analytical integration of Kulmala
             zbeta(kproma,klev)        ! auxiliary coefficient for the analytical integration of Kulmala

  !--- 0) Initialisations: ------------------------------------------------------

  zqtmst=1._dp/ztmst
  NS = 1
  !zeps=EPSILON(1.0_dp)

  !--- 1) Calculate nucleation rate:

  IF(nnucl==1) THEN

!CDIR NOIEXPAND
     CALL gnucl_vehkamaeki(kproma,   klev,         & ! ECHAM5 dimensions
                           ptemp,    prhum, pso4g, & ! ECHAM5 temperature, relative humidity
                           znucrate, pncrit          )
     ! pncrit is the number of nucleated particles (intent(out)) for vehkamaeki
     !--- Calculate updated gas phase concentration:
     !
     !    N(t)   = N(0)   - znucrate   * pncrit * dt
     !    [cm-3] = [cm-3] - [cm-3 s-1] * [1]    * [s]

     DO jk=1, klev
        DO jl=1, kproma

           zso4g_new(jl,jk)=pso4g(jl,jk)-(znucrate(jl,jk)*pncrit(jl,jk)*ztmst)

           IF(pso4g(jl,jk) > cmin_epsilon .AND. &
             znucrate(jl,jk) > cmin_epsilon) THEN
   !! KP I don't understand why you need pncrit here?  Its the same as M7 though
             zso4g_new(jl,jk) = pso4g(jl,jk) - (znucrate(jl,jk) * &
                                                pncrit(jl,jk) * ztmst)
           ELSE
!             zso4g_new(jl,jk)=MAX(pso4g(jl,jk),0.0_dp)
             zso4g_new(jl,jk)=pso4g(jl,jk)
           ENDIF
        END DO
     END DO

 ELSE IF (nnucl==2) THEN

     pncrit(:,:)=critn

!CDIR NOIEXPAND
     CALL gnucl_kulmala(kproma,   klev,            &
                        pso4g,    ptemp,  prhum,   &
                        znucrate, zalpha, zbeta    )

     !--- 2) Analytical integration of the nucleation rate (Eq. 19) ----------------
     !       over timestep ztmst assuming no new H2SO4(g) production:
     !
     !       d(N_av/critn)/dt=exp(alpha + beta*ln(N_av)) => ... =>
     !
     !       N_av(t0+dt)=
     !       [N_av(t0)**(1-beta) + critn*(1-beta)exp(alpha)*dt]**(1/(1-beta)
     !

     DO jk=1, klev
        DO jl=1, kproma

           IF (znucrate(jl,jk) .GT. 1e-10_dp) THEN
              zf1 = pso4g(jl,jk)**(1.0_dp-zbeta(jl,jk))-&
                          critn*EXP(zalpha(jl,jk))*(1.0_dp-zbeta(jl,jk))*ztmst
              zso4g_new(jl,jk) = EXP(LOG(zf1)/(1.0_dp - zbeta(jl,jk)))
           ELSE
              zso4g_new(jl,jk) = pso4g(jl,jk)
           END IF

        END DO
     END DO

  END IF


  !--- 3) Calculate changes in gas-phase and aerosol mass and aerosol numbers: --

  DO jk=1, klev
     DO jl=1, kproma

       IF( znucrate(jl,jk) > zeps  ) THEN

  !--- 3.1) Security check:

           zso4g_new(jl,jk) = MAX(zso4g_new(jl,jk), zero)
           zso4g_new(jl,jk) = MIN(zso4g_new(jl,jk), pso4g(jl,jk))

  !--- 3.2) Calculate mass of nucleated H2SO4 (equals the net
  !         gas phase H2SO4 loss):
  !
           !pa4delt(jl,jk,nso42m(NS)) = (pso4g(jl,jk)-zso4g_new(jl,jk))
           pa4delt(jl,jk,nso42m(NS)) = MAX((pso4g(jl,jk)-zso4g_new(jl,jk)),0.0_dp)  !mz_kp_20080112
           ! adding also H+ mass to compensate for H2SO4 = 2 H+ + SO4--
           pa4delt(jl,jk,nhp(NS))    = MAX((2._dp * &
             (pso4g(jl,jk)-zso4g_new(jl,jk))),0.0_dp)  ! mz_ht_20090310
  !
  !--- 3.3) Calculate the number of nucleated particles (nucleated mass
  !         divided by the assumed mass of a critical cluster critn):

           panew(jl,jk)=pa4delt(jl,jk,nso42m(NS))/pncrit(jl,jk)

  !--- 3.4) Calculate changes in gas phase H2SO4 due to nucleation:
  !
!           print*, pa4delt(jl,jk,nso42m(NS)), pso4g(jl,jk), jl,jk
           pso4g(jl,jk)=pso4g(jl,jk)-pa4delt(jl,jk,nso42m(NS))
        END IF

     END DO
  END DO

END SUBROUTINE gmxe_nuck
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~?
SUBROUTINE gmxe_delcoa(kproma,  klev,     paerml,   &
                      paernl,   prwet,  pa4delt,  panew,    &
                      pa,       pb,     pbfract1, pbfract2, &
                      pbfract5, pgas,   ztmst               )
  !
  !    Authors:
  !    ---------
  !    E. Vignati and J. Wilson, JRC/EI (original source)                09/2000
  !    P. Stier, MPI                    (f90-version, changes, comments)    2001
  !
  !    Modifications:
  !    --------------
  !    H. Tost, MPI-CHEM (extended and generalized for use with GMXE), Aug 2010
  !
  !    Version/History:
  !    ----------------
  !    equivalent to the version delco_n2 of the gmxe boxmodel
  !    + use of the analytical solution
  !
  !    Purpose
  !    ---------
  !    This routine calculates changes in number concentration of
  !    each aerosol mode over the time step, due to coagulation with
  !    the current mode and all higher ones.
  !
  !    Method:
  !    -----------
  !    *delcoa*  integrates for each mode dn/dt=c -a*n^2 -b*n  over ztmst
  !
  !    The resulting particles are assumed to reside in the
  !    mode of highest mode of the pair of particles colliding.
  !    1+1=>1 1+2=>2 1+3=>3, 2+2=>2 2+3=>3, 3+3=>3.
  !    zc is now non zero for mode 1 only (nucleation).
  !    All formation of higher mode particles is handled in dconc
  !
  !    For climatological studies, 5 day accumulation mode concs are
  !    within a factor of 2 of the full model.
  !
  !    Interface
  !    -----------
  !    *gmxe_delcoa* is called from *gmxe_dnum*
  !
  !    Externals
  !    -----------
  !    none
  !
  !--- Parameter list:
  !
  ! paerml          = total aerosol mass for each compound
  ! paernl          = aerosol number for each mode [cm-3]
  ! pa4delt(:,:,:)  = change in H2SO4 mass of the respective mode over one timstep
  !                   due to:
  !                      - nucleation of H2SO4 (calculated in gmxe_nuck)
  !                      - coagulation (calculated in gmxe_concoag)
  ! panew           = number of nucleated particles (during 1 timestep) [1]
  ! pa              = unimodal coagulation coefficient (zcom(mod)) []
  ! pb              = inter-modal coagulation with higher modes
  !                   (zcom(mod) * n(jmod+1)) []
  ! pbfractx(:,:,y) = fraction of the total number of particles removed by
  !                   coagulation from mode x that is moved to mode y+1 [1]
  !
  !--- Local variables:
  !
  ! zansum          = aerosol number in the respective mode [cm-3]
  ! zxxsum          = aerosol mass for compound xx in the respective
  !                   mode, e.g. xx = bc, oc, a4 (sulfate)
  !                   [g cm-3 for bc,oc and molec. cm-3 for sulfate]
  ! zxxav           = average mass of a sulfate particle in the respective
  !                   mode [molecules]
  ! zxxavy          = average mass of species xx in mode y []
  !                   where xx is ss, du, bc, oc, or a4 for sulfate
  !                   [molecules for sulfate and ug for others]
  ! zanli(:,:,:)    = Number of particles moved by the inter-modal
  !                   coagulation []
  ! zansq(:,:,:)    = Number of particles moved by the intra-modal
  !                   coagulation []
  ! zaernt(:,:,:)   = New particle number n(t+dt) after the integration
  !                   of the aerosol dynamics equation [cm-3]

  !--- Parameters:

  IMPLICIT NONE

  ! Local variables:
  !
  INTEGER :: kmod,         jmod,           jk,          jl,                    &
             jcomp,        jc,             kcomp,       kproma,  klev

  REAL(dp):: paerml(kproma,klev,0:naertot),   paernl(kproma,klev,nmod),        &
             pa4delt(kproma,klev,0:naertot),  panew(kproma,klev),              &
             prwet(kproma,klev,nmod)

  REAL(dp):: pbfract1(kproma,klev,nmod-1),  pbfract2(kproma,klev,nmod-1),      &
             pbfract5(kproma,klev,3)

  REAL(dp):: pa(kproma,klev,nmod),          pb(kproma,klev,nmod)


  REAL(dp):: za4av1(kproma,klev),           za4av2(kproma,klev)
  REAL(dp):: za4av(kproma,klev,nmod)

  REAL(dp):: zanli(kproma,klev,nmod),       zansq(kproma,klev,nmod),           &
             zaernt(kproma,klev,nmod),      zav(kproma,klev,0:naertot),        &
             pgas(kproma,klev,ngas,nmod)

  REAL(dp):: zrwet(nmod), zcrtcst(nmod)
  INTEGER :: imask(kproma,klev,nmod)
  REAL(dp):: zansum_vec(kproma,klev,nsoluble)

  ! Auxiliary variables:
  INTEGER :: idx1, idx2, jt

  REAL(dp):: zansum,                                                        &
             ztop,   zbot,   zatot,  zanse,  zanle,                         &
             ze1,    zf1,    zf2,    zf3,    zf4,    zr1,                   &
             ztmst,  zc

  REAL(dp):: zamloss,       zanloss,       ztloss,        zbfofac,          &
             zbfnfac,       zbftot,        zaerni,        zanytnt,          &
             zanytni,       zanytns,       zanytnm,       ztotn,            &
             zaerns
  LOGICAL :: npass ! mz_dk_20120120
!!$  REAL(dp):: zanli1(kproma,klev,nmod),       zansq1(kproma,klev,nmod),      &
!!$             zaernt1(kproma,klev,nmod),      paernl1(kproma,klev,nmod)
  !
  !--- 0) Initialisations: ------------------------------------------------
  !
  za4av1(:,:)  = zero
  za4av2(:,:)  = zero
  za4av(:,:,:) = zero
  zav (:,:,:)  = zero
  !
  !--- 1) Insoluble modes
  !
     !--- Insoluble modes:
     !    Calculate average mass (used in concoag)
  DO jmod=nsoluble+1,nmod
    DO jt=1,spec_number
      kcomp=species(jt)%aermlidx(jmod)
      IF(kcomp == 0 .OR. kcomp == species(jt)%aernlidx(jmod)) CYCLE
      DO jk=1,klev
        DO jl=1,kproma
          IF(paernl(jl,jk,jmod) > cmin_aernl) &
            zav(jl,jk,kcomp)=paerml(jl,jk,kcomp)/paernl(jl,jk,jmod)
        END DO
      END DO
    END DO
  END DO
  !
  ! mz_ht_20100805+
  ! insoluble modes
  DO jmod=nsoluble+1,nmod
    DO jk=1,klev
      DO jl=1,kproma
        zansum=paernl(jl,jk,jmod)
        zaernt(jl,jk,jmod)=zansum
        zanli(jl,jk,jmod)=zero
        zansq(jl,jk,jmod)=zero
        imask(jl,jk,jmod)=0
        !--- Calculations only in presence of sufficient particles:
        IF (zansum > cmin_aernl) THEN
          imask(jl,jk,jmod)=1
          !--- 1.2) Case with coagulation:
          ! either intra - modal or intermodal coagulation takes place
          IF (pa(jl,jk,jmod) >= cmin_epsilon .OR. &
              pb(jl,jk,jmod) >= cmin_epsilon) THEN
            !--- 1.2.1) Case of no inter-modal coagulation:
            !           dn/dt = -a*n**2 =>
            !           n(t)  = n0/(1 + n0*a*(t-t0))
            IF (pb(jl,jk,jmod) < cmin_epsilon) THEN
              zaernt(jl,jk,jmod)=zansum/(1.0_dp+zansum*pa(jl,jk,jmod)*ztmst)
              zanli (jl,jk,jmod)=zero
              zansq (jl,jk,jmod)=zansum-zaernt(jl,jk,jmod)
            ELSE
            !--- 1.2.2) Case with inter- and intra-modal coagulation:
            !           dn/dt = -a*n**2 - b*n =>
            !           n(t)  = (b*n0*exp(-b(t-t0)))/((n0*a)(1-exp(-b(t-t0)))+b)
            !--- Calculate n(t+dt):
              ze1=EXP(-pb(jl,jk,jmod)*ztmst)
              ztop=pb(jl,jk,jmod)*zansum*ze1
              zbot=zansum*pa(jl,jk,jmod)*(1.0_dp-ze1)+pb(jl,jk,jmod)
              zaernt(jl,jk,jmod)=ztop/zbot
              !--- Limit n(t+dt) to available particle in the mode:
              zaernt(jl,jk,jmod)=MIN(zaernt(jl,jk,jmod), zansum)
              !--- Total change in particle numbers of the mode
              !    due to coagulation:
              zatot=zansum-zaernt(jl,jk,jmod)
              !--- Contribution of the intra-modal coagulation:
              zanse=zansum*zansum*pa(jl,jk,jmod)
              !--- Contribution of the inter-modal coagulation:
              zanle=zansum*pb(jl,jk,jmod)
              !--- Number of particles moved by the inter-modal coagulation:
              zanli(jl,jk,jmod)=zatot*zanle/(zanse+zanle)
              !--- Number of particles moved by the intra-modal coagulation:
              zansq(jl,jk,jmod)=zatot*zanse/(zanse+zanle)
            END IF
          END IF
        END IF
      END DO
    END DO
    DO jt=1,spec_number
      kcomp=species(jt)%aermlidx(jmod)  ! component of aitken insoluble
      IF(kcomp == 0 .OR. kcomp==species(jt)%aernlidx(jmod)) CYCLE
      DO jm = 1,nmod                ! coagulating modes
        ! mode after coagulation larger than the soluble mode
        ! (originating from the insoluble mode)
        IF (idest(jmod,jm) > jm) THEN
          DO jk=1,klev
            DO jl=1,kproma
              IF (imask(jl,jk,jmod) == 0) CYCLE
          !-- 1.2.3) Change masses of the modes due to
          !          intra-modal coagulation and the coagulation with the
          !          smaller soluble modes (transfers according to a coating
          !          calculated later (concoag):
              paerml(jl,jk,kcomp)=(zaernt(jl,jk,jmod) +                        &
                                   zansq(jl,jk,jmod)  +                        &
                                   xbfract(jl,jk,jmod,jm)*zanli(jl,jk,jmod)) * &
                                   zav(jl,jk,kcomp)
            ENDDO
          ENDDO
        ELSE        ! mode after coagulation is the soluble mode

          jcomp = species(jt)%aermlidx(jm)
          IF (jcomp == 0 .OR. jcomp == species(jt)%aernlidx(jm)) CYCLE
          DO jk=1,klev
            DO jl=1,kproma
              IF (imask(jl,jk,jmod) == 0) CYCLE
            !--- 1.2.5) Store changes in masses of compounds in
            !           the insoluble aitken mode due to inter-modal
            !           coagulation:
            ! (zanli(:,:,x)   = total number of particles moved from mode x
            !  pbfract5(:,:,x)= fraction of the total number of particles
            !                   moved from mode 5 that is moved to mode x )
              pa4delt(jl,jk,jcomp) = xbfract(jl,jk,jmod,jm)* &
                                     zanli(jl,jk,jmod)*zav(jl,jk,kcomp)
            END DO
          END DO
        END IF
      END DO
    END DO
    !--- 1.2.4) Change the numbers of the insoluble
    !           modes due to inter-modal coagulation is done in concoag:

!    print*, "before inter-modal coag: ", jmod, paernl(1,31,:)

    DO jm = 1,nmod      ! coagulating modes
      IF (idest(jmod,jm) > jm) THEN
        DO jk=1,klev
          DO jl=1,kproma
            IF (imask(jl,jk,jmod) == 0) CYCLE
            paernl(jl,jk,jmod)= zaernt  (jl,jk,jmod)   +               &
                                xbfract(jl,jk,jmod,jm)*zanli(jl,jk,jmod)
          ENDDO
        END DO
      END IF
    ENDDO
!    print*, "after inter-modal coag: ", jmod, paernl(1,31,:)

  END DO      ! insoluble modes

  ! mz_ht_20100805-





  !
  !--- 2) Soluble modes: --------------------------------------------------
  !

! adding up the mass from nucleation to an updated mass in the soluble mode
! which then is used for coagulation
! new
  jmod = 1
  DO jt=1,spec_number
    kcomp=species(jt)%aermlidx(jmod)
    IF(kcomp == 0 .OR. kcomp == species(jt)%aernlidx(jmod)) CYCLE
    DO jk=1,klev
      DO jl=1,kproma
        paerml(jl,jk,kcomp)  = paerml(jl,jk,kcomp)+pa4delt(jl,jk,kcomp)
        pa4delt(jl,jk,kcomp) = 0._dp
      END DO
    END DO
  END DO


  DO jmod=1,nsoluble

    IF (jmod == 1) THEN
      DO jk=1,klev
        DO jl=1,kproma
          zansum_vec(jl,jk,jmod) = paernl(jl,jk,jmod) + panew(jl,jk)
        ENDDO
      ENDDO
    ELSE ! jmod /= NS
      DO jk=1,klev
        DO jl=1,kproma
          zansum_vec(jl,jk,jmod) = paernl(jl,jk,jmod)
        ENDDO
      ENDDO
    ENDIF
    DO jt=1,spec_number
      kcomp=species(jt)%aermlidx(jmod)
      IF(kcomp == 0 .OR. kcomp == species(jt)%aernlidx(jmod)) CYCLE
      DO jk=1,klev
        DO jl=1,kproma
          IF (zansum_vec(jl,jk,jmod) > cmin_aernl) &
              zav(jl,jk,kcomp) = paerml(jl,jk,kcomp) / zansum_vec(jl,jk,jmod)

        END DO
      END DO
    END DO ! jt - species

    DO jk=1,klev
      DO jl=1,kproma
        zaernt(jl,jk,jmod)=zansum_vec(jl,jk,jmod)
        zanli(jl,jk,jmod)=zero
        zansq(jl,jk,jmod)=zero
        imask(jl,jk,jmod)=0
      END DO
    END DO
  END DO ! jmod
  !

  za4av1(:,:)=zero
  za4av2(:,:)=zero


  DO jmod=1,nsoluble
    DO jt=1,spec_number
      kcomp=species(jt)%aermlidx(jmod)
      npass=species(jt)%npassive(jmod) ! mz_dk_20120120
      IF(kcomp == 0 .or. kcomp == species(jt)%aernlidx(jmod) .or. npass) CYCLE ! mz_dk_20120124, added npass
      !print*, "za4av",species(jt)%name, npass ! mz_dk_20120124
      DO jk=1,klev
        DO jl=1,kproma
          !--- Calculations only in presence of sufficient particles:
          IF(zansum_vec(jl,jk,jmod) > cmin_aernl) THEN

            IF (jmod == 1) THEN
              za4av1(jl,jk)=za4av1(jl,jk)+zav(jl,jk,kcomp)
            ELSE IF (jmod == 2) THEN
              za4av2(jl,jk)=za4av2(jl,jk)+zav(jl,jk,kcomp)
            END IF
            za4av(jl,jk,jmod) = za4av(jl,jk,jmod) + zav(jl,jk,kcomp)

            !--- 2.1) Case of no coagulation:
            !
            IF (pa(jl,jk,jmod) < cmin_epsilon .AND. &
                pb(jl,jk,jmod) < cmin_epsilon) THEN
              !--- Nucleation in the smallest (= first) mode only.
              !    Nothing to be done for other modes.
!!$              IF(jmod == NS) THEN
              IF(jmod == 1) THEN
                paernl(jl,jk,jmod) =zansum_vec(jl,jk,jmod)
              END IF
              !--- 2.2) Case with coagulation:
            ELSE
              imask(jl,jk,jmod) = 1
            END IF
          END IF ! (zansum_vec(jl,jk,jmod) > cmin_aernl)
        END DO ! jl
      END DO ! jk
    END DO ! jt - species
  END DO ! jmod
  !
  DO jmod=1,nsoluble
    DO jk=1,klev
      DO jl=1,kproma
        IF (imask(jl,jk,jmod) == 0) CYCLE
        !--- 2.2.1) Case of no nucleation:
        !--- Not Mode 1 or Nucleation rate below 1/s:
!!$        IF ( (jmod /= NS) .OR. (panew(jl,jk)/ztmst < 1.0_dp) ) THEN
        IF ( (jmod /= 1) .OR. (panew(jl,jk)/ztmst < 1.0_dp) ) THEN
          paernl(jl,jk,jmod)=zansum_vec(jl,jk,jmod)
          !--- 2.2.1a) Case of no inter-modal coagulation:
          !            dn/dt = -a*n**2 =>
          !            n(t)  = n0/(1 + n0*a*(t-t0))
          IF (pb(jl,jk,jmod) < cmin_epsilon) THEN  !! KP New
            zaernt(jl,jk,jmod) = zansum_vec(jl,jk,jmod) / &
              (1.0_dp + zansum_vec(jl,jk,jmod) * pa(jl,jk,jmod) * ztmst)
            zanli(jl,jk,jmod) = zero
            zansq(jl,jk,jmod) = zansum_vec(jl,jk,jmod) - zaernt(jl,jk,jmod)
            !--- 2.2.1b) Case with inter- and intra-modal coagulation:
            !            dn/dt = -a*n**2 - b*n =>
            !            n(t)  =(b*n0*exp(-b(t-t0)))/((n0*a)(1-exp(-b(t-t0)))+b)

          ELSE
            !--- Calculate n(t+dt):
            ze1  = EXP(-pb(jl,jk,jmod) * ztmst)
            ztop = pb(jl,jk,jmod) * zansum_vec(jl,jk,jmod) * ze1
            zbot = zansum_vec(jl,jk,jmod) * pa(jl,jk,jmod) &
                 * (1.0_dp - ze1) + pb(jl,jk,jmod)
            zaernt(jl,jk,jmod) = ztop / zbot
            !--- Limit n(t+dt) to available particle in the mode
            zaernt(jl,jk,jmod) = MIN(zaernt(jl,jk,jmod), zansum_vec(jl,jk,jmod))
            !--- Total change in particle numbers of the mode due to coagulation
            zatot = zansum_vec(jl,jk,jmod) - zaernt(jl,jk,jmod)
            !--- Contribution of the intra-modal coagulation
            zanse = zansum_vec(jl,jk,jmod) * zansum_vec(jl,jk,jmod) * &
                    pa(jl,jk,jmod)
            !--- Contribution of the inter-modal coagulation:
            zanle = zansum_vec(jl,jk,jmod) * pb(jl,jk,jmod)
            !--- Number of particles moved by the inter-modal coagulation:
            zanli(jl,jk,jmod) = zatot * zanle / (zanse + zanle)
            !--- Number of particles moved by the intra-modal coagulation:
            zansq(jl,jk,jmod) = zatot * zanse / (zanse + zanle)
          END IF
          !--- 2.2.2) Case with nucleation:
!        ELSE IF ( (jmod == NS) .AND. (panew(jl,jk)/ztmst >= 1.0_dp) ) THEN
        ELSE IF ( (jmod == 1) .AND. (panew(jl,jk)/ztmst >= 1.0_dp) ) THEN
          !--- 2.2.2a) Nucleation, inter- and intra-modal coagulation:
          !            dn/dt = -a*n**2 - b*n + c =>
          !            n(t)  = -(b/(2a)) +
          !                    R/2a * [((1 - (-2ax0-b+R)/(+2ax0+b+R))exp(-Rt)) /
          !                             ((1 + (-2ax0-b+R)/(+2ax0+b+R))exp(-Rt))]
          !            where:  R=SQRT(b**2+4ac)
          !
          !           If b/=0 then always a/=0. The only case where a would be 0
          !           and b unequal 0 is the case of no pre-existing particles
          !           in the nucleation mode but pre-existing particles in other
          !           modes. For this case a is calculated for an assumed radius
          !           of a critical cluster in gmxe_coaset.

          IF (pb(jl,jk,jmod) >= cmin_epsilon) THEN
            !--- Calculate n(t):
            !--- c:
            zc = panew(jl,jk) / ztmst
            !--- R:
            zf1 = pb(jl,jk,jmod) * pb(jl,jk,jmod) + 4.0_dp * pa(jl,jk,jmod) * zc
            zr1 = SQRT(zf1)
            !--- exp(-Rt):
            ze1 = EXP(-zr1*ztmst)
            !--- 2ax0+b:
            zf2=2.0_dp*pa(jl,jk,jmod)*paernl(jl,jk,jmod)+pb(jl,jk,jmod)
            !--- Term in squared bracket:
            zf3=ze1*(zr1-zf2)/(zr1+zf2)
            zf4=(1.0_dp-zf3)/(1.0_dp+zf3)
            !--- n(t):
            zaernt(jl,jk,jmod)=(zr1*zf4-pb(jl,jk,jmod))/2.0_dp/pa(jl,jk,jmod)
            !--- Limit n(t+dt) to available particle in the mode:
            zaernt(jl,jk,jmod)=MIN(zaernt(jl,jk,jmod), zansum_vec(jl,jk,jmod))
            !--- Total change in particle numbers of the mode due to coagulation:
            zatot=zansum_vec(jl,jk,jmod)-zaernt(jl,jk,jmod)
            !--- Contribution of the intra-modal coagulation:
            zanse=zansum_vec(jl,jk,jmod)*zansum_vec(jl,jk,jmod)*pa(jl,jk,jmod)
            !--- Contribution of the inter-modal coagulation:
            zanle=zansum_vec(jl,jk,jmod)*pb(jl,jk,jmod)
            !--- Number of particles moved by the inter-modal coagulation:
            zanli(jl,jk,jmod)=zatot*zanle/(zanse+zanle)
            !--- Number of particles moved by the intra-modal coagulation:
            zansq(jl,jk,jmod)=zatot*zanse/(zanse+zanle)

            !--- 2.2.2b) Nucleation and intra-modal coagulation:
            !            dn/dt = -a*n**2 - b*n + c with b=0 =>
            !            dn/dt = -a*n**2 + c =>
            !            n(t)  = R/2a * [ ((1 - (-2ax0+R)/(+2ax0+R))exp(-Rt)) /
            !                             ((1 + (-2ax0+R)/(+2ax0+R))exp(-Rt))  ]
            !            where:  R=SQRT(4ac)
            !
            !            Can be shown to be equivalent to:
            !
            !            n(t)  = R1*((x0+R1)/(x0-R1)+exp(-SQRT(-4ac)t)) /
            !                       ((x0+R1)/(x0-R1)-exp(-SQRT(-4ac)t))
            !            where R1=SQRT(c/a)

          ELSE IF (pb(jl,jk,jmod) < cmin_epsilon) THEN
            !--- c:
            zc=panew(jl,jk)/ztmst
            !--- R1:
            zr1=SQRT(zc/pa(jl,jk,jmod))
            !--- exp(-Rt):
            ze1=EXP(-zr1*2.0_dp*pa(jl,jk,jmod)*ztmst)
            !--- n(t):
            zf1=(paernl(jl,jk,jmod)+zr1)/(paernl(jl,jk,jmod)-zr1)
            ztop=zr1*(zf1+ze1)
            zbot=zf1-ze1
            IF (zbot < cmin_epsilon) THEN
              zaernt(jl,jk,jmod)=zansum_vec(jl,jk,jmod)
            ELSE
              zaernt(jl,jk,jmod)=ztop/zbot
            END IF
            !--- Limit n(t+dt) to available particle in the mode:
            zaernt(jl,jk,jmod)=MIN(zaernt(jl,jk,jmod), zansum_vec(jl,jk,jmod))
            !--- Number of particles moved by the inter-modal coagulation:
            zanli(jl,jk,jmod)=zero
            !--- Number of particles moved by the intra-modal coagulation:
            zansq(jl,jk,jmod)=zansum_vec(jl,jk,jmod)-zaernt(jl,jk,jmod)

          END IF
        END IF
      END DO ! jl - kproma
    END DO   ! jk - klev
  END DO     ! jmod
  !
  DO jmod=1,nsoluble
    DO jk=1,klev
      DO jl=1,kproma
        IF (imask(jl,jk,jmod) == 0) CYCLE
        !---2.2.3 New bit for insoluble/souble coagulation
        !--- sum total insoluble+soluble paticles in mode jmod JJNW
!!$        IF (jmod == NS .AND. zanli(jl,jk,jmod) > zero) THEN
!!$          zaerni=paernl(jl,jk,iaiti)+paernl(jl,jk,iacci)+paernl(jl,jk,icoai)
!!$          zaerns=zansum_vec(jl,jk,jmod)+paernl(jl,jk,iaits)
        IF (jmod == 1 .AND. zanli(jl,jk,jmod) > zero) THEN
          zaerni=sum(paernl(jl,jk,nsoluble+1:nmod))
!!$          zaerns=zansum_vec(jl,jk,jmod)+sum(paernl(jl,jk,2:nsoluble))
          zaerns=zansum_vec(jl,jk,jmod)+paernl(jl,jk,jmod+1)
          ztotn=zaerns+zaerni
          IF (zaerns > zaerni .and. zaerni > zero) THEN
            !calculate analytical solution no of mixed particles for coagulation
            !between paernl(jl,jk,jmod) soluble particles and zaerni insouble of
            !the same dimensions
            IF (zaerni > 1.0_dp) then
              zanytni=4.0_dp*zaerni/((2.0_dp+pa(jl,jk,jmod)*ztmst*ztotn)*  &
                (2.0_dp+pa(jl,jk,jmod)*ztmst*(ztotn-zaerni)))
            ELSE
              zanytni = zero
            END IF
            zanytnt=2.0_dp*ztotn/(2.0_dp+pa(jl,jk,jmod)*ztmst*ztotn)
            zanytns=4.0_dp*zaerns/((2.0_dp+pa(jl,jk,jmod)*ztmst*ztotn)*       &
              (2.0_dp+pa(jl,jk,jmod)*ztmst*(ztotn-zaerns)))
            zanytnm=zanytnt-(zanytni+zanytns)

            !scale analytical solution to real aernt
            !                   zanytnm=zatot/(ztotn-zanytnt)*zanytnm
            zanytnm=min(zanytnm,zaerni)
!CDIR UNROLL=nmod
            DO kmod =1,nmod
              zrwet(kmod) = prwet(jl,jk,kmod)
            END DO
            CALL gmxe_coat_0d(zrwet, zcrtcst)

            ztloss=zero
!CDIR UNROLL=nmod
            DO kmod = nsoluble+1,nmod
              zamloss=paernl(jl,jk,kmod)/zaerni*zanytnm*zcrtcst(kmod)
!!$              IF (za4av1(jl,jk)>zero) THEN
!!$                zanloss=zamloss/za4av1(jl,jk)
              IF (za4av(jl,jk,jmod) > zero) THEN
                zanloss=zamloss/za4av(jl,jk,jmod)
              ELSE
                zanloss=zero
              END IF
              ztloss=ztloss+zanloss
            END DO
            ztloss=min(ztloss,zansq(jl,jk,jmod)*0.95_dp)
            zbfofac=zanli(jl,jk,jmod)/(zanli(jl,jk,jmod)+ztloss)
            zbfnfac=ztloss/(zanli(jl,jk,jmod)+ztloss)
            zanli(jl,jk,jmod)=zanli(jl,jk,jmod)+ztloss
            zansq(jl,jk,jmod)=zansq(jl,jk,jmod)-ztloss
            zbftot=zero
!CDIR UNROLL=nmod
            DO kmod=1,nmod
              IF(kmod > jmod) THEN
                pbfract1(jl,jk,kmod-jmod)=pbfract1(jl,jk,kmod-jmod)*zbfofac
                IF (kmod > nsoluble) THEN
                  pbfract1(jl,jk,kmod-jmod)=pbfract1(jl,jk,kmod-jmod)+      &
                    zbfnfac*paernl(jl,jk,kmod)/zaerni
                ENDIF
                xbfract(jl,jk,jmod,kmod-jmod) = &
                  xbfract(jl,jk,jmod,kmod-jmod) * zbfofac
                IF (kmod > nsoluble) THEN
                  xbfract(jl,jk,jmod,kmod-jmod) =   &
                    xbfract(jl,jk,jmod,kmod-jmod) + &
                    zbfnfac*paernl(jl,jk,kmod)/zaerni
                END IF
!!$                zbftot=zbftot+pbfract1(jl,jk,kmod-jmod)
                zbftot=zbftot+xbfract(jl,jk,jmod,kmod-jmod)
              END IF
            END DO
!CDIR UNROLL=nmod
            DO kmod=1,nmod
              IF (kmod > jmod) THEN
                pbfract1(jl,jk,kmod-jmod)=pbfract1(jl,jk,kmod-jmod)/zbftot
                xbfract(jl,jk,jmod,kmod-jmod) = &
                  xbfract(jl,jk,jmod,kmod-jmod) / zbftot
              END IF
            END DO
          END IF
        END IF
        !---- End of new inslouble/soluble caogulation routine JJNW
      END DO
    END DO
  END DO
  !
  DO jmod=1,nsoluble
    DO jt=1,spec_number
      kcomp=species(jt)%aermlidx(jmod)
      IF(kcomp == 0 .OR. kcomp == species(jt)%aernlidx(jmod)) CYCLE
      DO jk=1,klev
        DO jl=1,kproma
          IF (imask(jl,jk,jmod) == 0) CYCLE
          !- 2.3) Change masses and numbers of the respective modes to account--
          !         for intra-modal coagulation (zansq) and coagulation with
          !         higher modes (zaernt):
          !
          !--- 2.3.1. - 2.3.2) Change mass of aerosol compounds:
          paerml(jl,jk,kcomp) = (zaernt(jl,jk,jmod) + zansq(jl,jk,jmod)) * &
                                 zav(jl,jk,kcomp)
        END DO
      END DO
    END DO
  END DO

  !
  DO jmod=1,nsoluble
    DO jk=1,klev
      DO jl=1,kproma
        IF (imask(jl,jk,jmod) == 0) CYCLE
        !--- 2.3.3) Particle numbers:
        paernl(jl,jk,jmod)=zaernt(jl,jk,jmod)
      END DO
    END DO
  END DO
  !

  DO jmod=1,nsoluble

    do jt=1,spec_number
      kcomp = species(jt)%aermlidx(jmod)
      if (kcomp == 0 .OR. kcomp == species(jt)%aernlidx(jmod)) CYCLE
      DO kmod = jmod + 1, nsoluble
        jcomp = species(jt)%aermlidx(kmod)
        if (jcomp == 0 .OR. jcomp == species(jt)%aernlidx(kmod)) CYCLE
        DO jk=1,klev
          DO jl=1,kproma
            IF (imask(jl,jk,jmod) == 0) CYCLE
          !--- 2.4) Calculate changes in particle masses due to inter-modal ---
          !         coagulation:
          !--- 2.4.1) Transfer of mass from mode jmod to higher modes:
          ! Mass from (jmod) to higher solbule modes (kmod):

            pa4delt(jl,jk,jcomp) = pa4delt(jl,jk,jcomp) +   &
                                   xbfract(jl,jk,jmod,kmod-1) * &
                                   zanli(jl,jk,jmod) * zav(jl,jk,kcomp)

          ! Mass from (jmod) to higher soluble modes (kmod) due to coag.
          ! with insoluble modes (kmod+nsoluble-iaits):
            pa4delt(jl,jk,jcomp) = pa4delt(jl,jk,jcomp) +  &
                                   xbfract(jl,jk,jmod,kmod+nsoluble-1-ndiff)*&
                                   zanli(jl,jk,jmod) * zav(jl,jk,kcomp)
!            ENDIF
          ENDDO !kproma
        END DO  ! klev
      END DO    ! kmod
    END DO      ! species

  END DO  ! mode loop



!  print*, "before concoag: ", paernl(1,31,:)
  !
  !--- 3) Calculate transfer from the insoluble to the soluble modes: ---------
!CDIR NOIEXPAND
  CALL gmxe_concoag (kproma,   klev,                                 &
                     paerml,   paernl,   prwet,   pa4delt, zanli,    &
                     za4av1,   za4av2,   zav,     pgas,             &
                     pbfract1, pbfract2                             ,&
                     za4av )
!  print*, "after concoag: ", paernl(1,31,:)

  !
  !--- 4) Final change of the aerosol masses due to nucleation, ---------------
  !       inter-modal coagulation and condensation on the insoluble modes:
  !       (Nucleation mode already done above.)
  !

  DO jmod=1,nmod
    DO jt=1,spec_number
      kcomp=species(jt)%aermlidx(jmod)
      IF(kcomp == 0 .OR. kcomp == species(jt)%aernlidx(jmod)) CYCLE
      paerml(1:kproma,1:klev,kcomp) = paerml(1:kproma,1:klev,kcomp) + &
                                      pa4delt(1:kproma,1:klev,kcomp)
    END DO
  END DO



  !
END SUBROUTINE gmxe_delcoa
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
SUBROUTINE gmxe_concoag (kproma,   klev,                                &
                         paerml,   paernl,   prwet,   pa4delt, panli,   &
                         pa4av1,   pa4av2,   pav,     pgas,             &
                         pbfract1, pbfract2, pa4av             )
  !
  !   *gmxe_concoag*
  !
  !   Author:
  !   ----------
  !   E. Vignati, JRC/EI     (original source)                09/2000
  !   P. Stier, MPI          (f90-version, changes, comments)    2001
  !
  !   Modifications:
  !   --------------
  !   H. Tost, MPI-CHEM (restructured and generalised for GMXE)
  !   S. Metzger, MPI-CHEM (extended and generalized for use with EQSAM3), June 2008
  !
  !   Version:
  !   ----------
  !   This version is equivalent to the version concoa_n of the boxmodel.
  !
  !   Purpose
  !   ----------
  !   gmxe_concoag transfers aerosol mass and numbers from the insoluble
  !   to the soluble modes.
  !
  !   Interface:
  !   ----------
  !   *gmxe_concoag* is called from *gmxe_delcoa*
  !
  !   Externals
  !   ----------
  !   none
  !
  !--- Parameters:
  !
  ! paerml          = total aerosol mass for each compound
  ! paernl          = aerosol number for each mode [cm-3]
  ! pa4delt(:,:,:)  = change in H2SO4 mass of the respective mode over one timstep
  !                   due to:
  !                      - nucleation of H2SO4 (calculated in gmxe_nuck)
  !                      - coagulation (calculated here in gmxe_concoag)
  ! pxxavy          = average mass of species xx in mode y []
  ! panli(:,:,x)    = total number of particles moved by inter-modal
  !                   coagulation from mode x [cm-3]
  ! pbfractx(:,:,y) = fraction of the total number of particles removed by
  !                   coagulation from mode x that is moved to mode y+1 [1]
  !
  !--- Local variables / Constants:
  !
  ! zcrtcst  = Critical constant, i.e. number of sulfate molecules to cover
  !            an average particle of the mode with a layer of the thickness
  !            determined by cLayerThickness in messy_gmxe. Calculated by
  !            gmxe_coat.
  !
  ! zcrit_x  = total available number of particles in mode x that are moved from
  !            insoluble mode x to the corresponding soluble mode.

  IMPLICIT NONE

  INTEGER :: kproma, klev

  REAL(dp):: pa4av1(kproma,klev),          pa4av2(kproma,klev),               &
             zgas(kproma,klev,nmod),       pav(kproma,klev,0:naertot),        &
             pa4av(kproma,klev,nmod)

  REAL(dp):: paerml(kproma,klev,0:naertot),paernl(kproma,klev,nmod),         &
             pbfract1(kproma,klev,nmod-1), pbfract2(kproma,klev,nmod-1),     &
             panli(kproma,klev,nmod),      pa4delt(kproma,klev,0:naertot),   &
             prwet(kproma,klev,nmod),      pgas(kproma,klev,ngas,nmod)


  ! Local variables:
  !
  INTEGER :: kmod, jmod, jk, jl, jcomp, jc, kcomp, idx1, idx2, jt

  REAL(dp):: zcrit(nmod), zcrit_vec(kproma,klev,nmod), zsum_vec(kproma,klev,nmod), zscale

  ! Auxiliary variables:

  REAL(dp):: zcrtcst(kproma,klev,nmod)

  !--- 0) Initializations:

  !zeps=EPSILON(1._dp)
  zcrit(:)=zero
  zcrit_vec(:,:,:)=zero

  !--- 1) Redistribution of mass and numbers after nucleation, coagulation ----
  !       and coagulation calculated in the preceeding subroutines:

  !--- 1.1) Determine number of particles that can be sufficiently coated
  !         by the available sulfate to be transfered to the soluble modes:

!CDIR NOIEXPAND
  CALL gmxe_coat(kproma, klev, prwet, zcrtcst)

  !--- 1.2) Sum masses added to insoluble modes due to
  !         coagulation with modes 1 and 2 (1st term) and the mass
  !         of gases condensed on the insoluble mode x (zgas):

  zgas(1:kproma,1:klev,1:nmod)=zero

  DO jmod=1,nmod
!!$     kmod = jmod
!!$     if (jmod > nsoluble ) kmod = jmod + ndiff - nsoluble
     jc = td%gas_sulph_idx
     IF (jc == 0) CYCLE
     zgas(1:kproma,1:klev,jmod)=zgas(1:kproma,1:klev,jmod) +&
                                pgas(1:kproma,1:klev,jc,jmod)
  ENDDO

  DO kmod = nsoluble + 1, nmod
    zsum_vec(1:kproma,1:klev,kmod) = zgas(1:kproma,1:klev,kmod)
    DO jmod = 1,nsoluble
      DO jk=1,klev
        DO jl=1,kproma
          zsum_vec(jl,jk,kmod) = zsum_vec(jl,jk,kmod) + panli(jl,jk,jmod) * &
                                 xbfract(jl,jk,jmod,kmod-1) * pa4av(jl,jk,jmod)
        END DO
      END DO
    END DO
  END do



  !--- 1.3) Number of particles moved from the insoluble to soluble modes due to
  !         interaction with 1 and due to condensation:
  DO kmod = nsoluble + 1,nmod
    jmod = kmod - nsoluble + ndiff
    DO jk=1,klev
      DO jl=1,kproma
        IF(paernl(jl,jk,kmod) >= cmin_aernl .AND. &
          zcrtcst(jl,jk,kmod) > zeps) THEN
          zcrit_vec(jl,jk,kmod) = MIN(paernl(jl,jk,kmod), &
                                  zsum_vec(jl,jk,kmod)/zcrtcst(jl,jk,kmod))
          paernl(jl,jk,jmod)=paernl(jl,jk,jmod)+zcrit_vec(jl,jk,kmod)
          paernl(jl,jk,kmod)=paernl(jl,jk,kmod)-zcrit_vec(jl,jk,kmod)
        END IF
      END DO
    END DO
  END DO


  DO jmod=1,nmod
     kmod = jmod
     if (jmod > nsoluble ) kmod = jmod + ndiff - nsoluble
     jc = td%gas_sulph_idx
     jcomp = nso42m(kmod)
     IF ((jc == 0) .OR. (jcomp == 0)) CYCLE
     pa4delt(1:kproma,1:klev,jcomp)     = &
          pa4delt(1:kproma,1:klev,jcomp) + &
          pgas(1:kproma,1:klev,jc,jmod)

     jcomp = nhp(kmod)
     pa4delt(1:kproma,1:klev,jcomp)     = &
          pa4delt(1:kproma,1:klev,jcomp) + &
          2._dp * pgas(1:kproma,1:klev,jc,jmod)

     pgas(1:kproma,1:klev,jc,jmod) = 0._dp
  END DO
  !--- 1.4) pa4delt: Mass moved   from mode insoluble to soluble:
  !--- 1.5) paerml:  Mass remaining in mode insoluble:


!!$  DO jmod = iaits, nsoluble      ! loop over soluble modes
!!$    kmod  = jmod + nsoluble - 1  ! insoluble mode index
  DO kmod = nsoluble+1, nmod        ! insoluble mode index
    jmod  = kmod - nsoluble + ndiff
    DO jt = 1,spec_number        ! loop over all species
      kcomp = species(jt)%aermlidx(kmod)    ! insoluble species index
      IF (kcomp == 0 .OR. kcomp == species(jt)%aernlidx(kmod)) CYCLE  ! species does not exist in insoluble mode
      jcomp = species(jt)%aermlidx(jmod)    ! soluble species index
      IF (jcomp == 0 .OR. jcomp == species(jt)%aernlidx(jmod)) CYCLE   ! species does not exist in soluble mode
      DO jk=1,klev
        DO jl=1,kproma
          zscale = zcrit_vec(jl,jk,kmod) * pav(jl,jk,kcomp)
          pa4delt(jl,jk,jcomp) = pa4delt(jl,jk,jcomp) + zscale
          paerml (jl,jk,kcomp) = paerml (jl,jk,kcomp) - zscale
        ENDDO
      END DO
    END DO
  END DO
! mz_ht_20080112-

END SUBROUTINE gmxe_concoag
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  SUBROUTINE gmxe_coat_0d(prwet_loc, pcrtcst_loc)

    ! Purpose:
    ! ---------
    ! *gmxe_coat* calculates the number of sulfate
    !           molecules required to coat a particle
    !           with cLayerThickness of sulfate
    !
    ! Author:
    ! ---------
    ! Philip Stier, MPI                          2001
    !
    ! Modifications:
    ! --------------
    ! H. Tost, UMZ generalised species structure
    ! S. Metzger, MPI-CHEM (extended and generalized for use with EQSAM3), June 2008
    !
    ! Interface:
    ! ---------
    ! *gmxe_coat* is called from *gmxe_concoag*
    !

    IMPLICIT NONE

    INTEGER         :: jmod_loc

    REAL(dp)        :: prwet_loc(nmod)          ! Ambient radii [cm]
    REAL(dp)        :: pcrtcst_loc(nmod)        ! Critical constant, i.e. number of
                                                ! sulfate to cover an average particle
                                                ! of the mode with a layer of the
                                                ! thickness determined by cLayerThickness.
    REAL(dp)        :: zras(nmod)               ! Radius of average surface
                                                ! for a single particle [cm]
    REAL(dp)        :: zas(nmod)                ! Average surface
                                                ! for single particle [cm+2]


  DO jmod_loc=nsoluble + 1,nmod
     !--- 1) Calculate the radii of average surface for modes 5-7:

     zras(jmod_loc) = prwet_loc(jmod_loc) * cmr2ras(jmod_loc)

     !--- 2) Calculate the average surface of an particle for modes 5-7:

     zas(jmod_loc)    = 4._dp * zras(jmod_loc)**2 * pi

     !--- 3) Determine the number of sulfate molecules needed to form
     !       a cLayerThickness thick layer of sulfate on the particles
     !       in modes 5-7:

     pcrtcst_loc(jmod_loc) = (zas(jmod_loc) / csurf_molec) * cLayerThickness

  END DO

  END SUBROUTINE gmxe_coat_0d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  SUBROUTINE gmxe_coat(kproma, klev, prwet, pcrtcst)

    ! Purpose:
    ! ---------
    ! *gmxe_coat* calculates the number of sulfate
    !           molecules required to coat a particle
    !           with cLayerThickness of sulfate
    !
    ! Author:
    ! ---------
    ! Philip Stier, MPI                          2001
    !
    ! Modifications:
    ! --------------
    ! H. Tost, UMZ, generalised species structure
    ! S. Metzger, MPI-CHEM (extended and generalized for use with EQSAM3), June 2008
    !
    ! Interface:
    ! ---------
    ! *gmxe_coat* is called from *gmxe_concoag*
    !

    IMPLICIT NONE

    INTEGER         :: kproma, klev, jmod, jk, jl

    REAL(dp)        :: prwet(kproma,klev,nmod)   ! Ambient radii [cm]
    REAL(dp)        :: pcrtcst(kproma,klev,nmod) ! Critical constant, i.e. number of
                                                ! sulfate to cover an average particle
                                                ! of the mode with a layer of the
                                                ! thickness determined by cLayerThickness.
    REAL(dp)        :: zras(nmod)               ! Radius of average surface
                                                ! for a single particle [cm]
    REAL(dp)        :: zas(nmod)                ! Average surface
                                                ! for single particle [cm+2]


  DO jmod=nsoluble + 1,nmod
     DO jk=1,klev
        DO jl=1,kproma

           !--- 1) Calculate the radii of average surface for modes 5-7:

           zras(jmod) = prwet(jl,jk,jmod) * cmr2ras(jmod)

           !--- 2) Calculate the average surface of an particle for modes 5-7:

           zas(jmod)    = 4._dp * zras(jmod)**2 * pi

           !--- 3) Determine the number of sulfate molecules needed to form
           !       a cLayerThickness thick layer of sulfate on the particles
           !       in modes 5-7:

           pcrtcst(jl,jk,jmod) = (zas(jmod) / csurf_molec) * cLayerThickness

        END DO
     END DO
  END DO

  END SUBROUTINE gmxe_coat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
END SUBROUTINE gmxe_dnum
!------------------------------------------------------------------------------!
SUBROUTINE gmxe_dconc (kproma, klev, paerml, paernl, prad)
  !
  !    *gmxe_dconc*  changes aerosol numbers and masses to account for
  !                  condensational growth of the mode mean radii
  !
  !    Authors:
  !    --------
  !    J. Wilson and E. Vignati, JRC (original source)            May 2000
  !    P. Stier, MPI-MET (f90 version, changes, comments)             2001
  !
  !    Modifications:
  !    --------------
  !    S. Metzger, MPI-CHEM (extended and generalized for use with EQSAM3), June 2008
  !    H. Tost, U-MZ (gerneralised and updated), Dec 2012
  !    Purpose:
  !    --------
  !    This routine repartitions aerosol number and mass between the
  !    the modes to account for condensational growth and the formation
  !    of an accumulation mode from the upper tail of the aitken mode.
  !
  !    Interface:
  !    ----------
  !    *gmxe_dconc* is called from *gmxe*
  !
  !    Method:
  !    -------
  !    The routine calculates the cumulativ number and mass distribution of the
  !    modes up to the respective mode boundary:
  !
  !                        / x                              _
  !                 N      |       1           1   ln(R)-ln(R)  2
  !    N(0,x) = ---------  |   --------  exp(- - ( ----------- )   ) d ln(R)
  !             ln(sigma)  |   sqrt(2PI)       2    ln(sigma)
  !                        / 0
  !
  !                         /tx                   2
  !                        |        1            t
  !           =     N      |     --------  exp(- - ) d t
  !                        |     sqrt(2PI)       2
  !                        /-inf
  !
  !    where:
  !
  !                        _
  !               ln(R)-ln(R)
  !    t      =   -----------
  !                ln(sigma)
  !
  !    and:
  !                        _
  !               ln(x)-ln(R)
  !    tx     =   -----------
  !                ln(sigma)
  !    _
  !    R is the Count Mean Radius or the Mass Mean Radius.
  !
  !    Practically, the routine gmxe_cumnor calculates the fraction of the number and
  !    mass distribution for each mode lying below the respective upper mode boundary (1).
  !    In a next step the net fraction of each mode lying between the upper and lower
  !    mode boundaries are summed up (2) and the numbers and masses exceeding the mode
  !    boundaries are transfered to the neighboring larger mode (3).
  !    Finally, these quantities are stored in the respective arrays
  !    paernl and paerml (4).
  !    The repartititioning is currently only done for the soluble modes as it is
  !    assumed that insoluble modes are rather transfered to the soluble modes
  !    and grow as soluble particles.
  !
  !    Externals:
  !    ----------
  !    None
  !
  !--- Parameter list:
  !
  !    paerml(kproma,klev,0:naertot)= total aerosol mass for each compound
  !                                   in [molecules cm^-3]
  !    paernl(kproma,klev,nmod)     = aerosol number for each mode [cm-3]  !
  !    sigma(jmod)                  = standard deviation of mode jmod [1]
  !    crdiv              = threshold radii between the different modes [cm]
  !                         crdiv(jmod) is the lower bound and crdiv(jmod+1) is
  !                         the upper bound of the respective mode
  !
  !--- Local Variables:
  !
  !    zfconn(:,:,jnum,jmod)  = absolute fraction of the number of particles
  !                             in mode jmod,
  !                             with CMD=2*prad (jmod) and a geometric standard
  !                             deviation zrcsig, that are smaller than
  !                             crdiv(jnum+1).
  !                             I.e. with 0 < CMD < crdiv(jnum+1) [1]
  !    zfconm(:,:,jnum,jmod)  = absolute fraction of the mass in mode jmod,
  !                             with CMD=2*prad (jmod) and a geometric standard
  !                             deviation zrcsig, that are smaller than
  !                             crdiv(jnum+1).
  !                             I.e. with 0 < CMD < crdiv(jnum+1) [1]

  IMPLICIT NONE

  INTEGER :: kproma, klev, kmod
  INTEGER :: idx01, idx02, idx1, idx2, idx_h2so4
  REAL(dp):: paerml(kproma,klev,0:naertot),   paernl(kproma,klev,nmod),       &
             prad (kproma,klev,nmod)

  REAL(dp):: paerml_pre(kproma,klev,0:naertot)
  REAL(dp):: paerml_post(kproma,klev,0:naertot)
  REAL(dp):: pre(spec_number), post(spec_number)

  ! Local variables:
  !
  INTEGER :: jnum,         jmod,         jk,          jl,       jc,     kcomp

  REAL(dp):: zrcsig,       zarg1,        zarg2,       zdpmm,         zdpcm,   &
             zarg3,        zdpam,        zcongn,      zcongm,        zdummy,  &
             zr1,          zr2,          zttnj,       zavnj,         zmrj,    &
             zmt,          znt,          zavmt,       zmcr,          zfconmj, &
             zntnew,       zmtnew,       zdm,         znewvol,       ztotvol, &
             zdpbm,        zarg4

  REAL(dp):: sumn_vec(kproma,nmod), summ_vec(kproma,nmod,naertot), summ2_vec(kproma,naertot)

  REAL(dp):: ztotmass(kproma,klev)

  REAL(dp):: zfconn(kproma,klev,nmod,nmod), &
             zfconm(kproma,klev,nmod,nmod)
  REAL(dp):: zmass_post2(kproma,klev), zmass_pre2(kproma,klev)

  LOGICAL :: npass ! mz_dk_20120120
  !--- 0) Initialisations: ----------------------------------------------------------

  zfconm(:,:,:,:) = zero
  zfconn(:,:,:,:) = zero
  !
  zmass_post2(:,:) = 0._dp
  zmass_pre2(:,:) = 0._dp
  !--- 1) Identify how much the mode jmod has grown into the next higher mode -------
  !
  idx_h2so4 = species(spec_idx_h2so4)%aermlidx(0)
  DO jmod=1,nsoluble-1

     !--- Total mass of the mode in equivalent molecules of sulphate:

     ztotmass(1:kproma,1:klev) = zero
     DO jc=1,spec_number
       kcomp=species(jc)%aermlidx(jmod)
       npass=species(jc)%npassive(jmod) ! mz_dk_20120124
       IF(kcomp == 0 .OR. kcomp == species(spec_idx_h2o)%aermlidx(jmod) &
            .or. npass) CYCLE ! mz_dk_20120124, added npass
       IF (kcomp == species(jc)%aernlidx(jmod)) CYCLE
       DO jk=1,klev
         DO jl=1,kproma
           ztotmass(jl,jk)  = ztotmass(jl,jk) + paerml(jl,jk,kcomp)
         END DO
       END DO
     END DO

     DO jnum=jmod,nsoluble-1
        DO jk=1,klev
           DO jl=1,kproma
              IF (paernl(jl,jk,jmod) > cmin_aernl .AND. &
                  prad (jl,jk,jmod) > ZERO) THEN

                 !--- 1.1) Calculate necessary parameters:

                 !--- Geometric Standard Deviation:
                 zrcsig=sigma_ln(jmod) ! LOG(sigma(jmod))

                 !--- Mass Median Radius:
                 zarg1=prad (jl,jk,jmod) * cmr2mmr(jmod)

                 !--- Count Median Radius:
                 zarg2=prad (jl,jk,jmod)

                 !--- Threshold radius between the modes:
                 zarg3=crdiv(jnum+1)
                 zarg4=crdiv(jnum+1) * cmr2mmr(jmod)
                 !--- Transfer to logarithmic scale:
                 zdpmm=LOG(zarg1)
                 zdpcm=LOG(zarg2)
                 zdpam=LOG(zarg3)
                 zdpbm=LOG(zarg4)

                 !--- Distance of the CMD of the mode from the threshold mode
                 !    diameter in terms of geometric standard deviations:

                 zcongn=(zdpam-zdpcm)/zrcsig

                 !--- Distance of the MMD of the mode from the threshold mode
                 !    diameter in terms of geometric standard deviations (t):

                 zcongm=(zdpbm-zdpmm)/zrcsig

                 !--- Calculate the cumulative of the log-normal number distribution:

                 CALL gmxe_cumnor(zcongn,zfconn(jl,jk,jnum,jmod),zdummy)

                 !--- Limit transfer only to adjacent modes:

                 IF (jnum > jmod) THEN
                    zfconn(jl,jk,jnum,jmod)= 1.0_dp
                    zfconm(jl,jk,jnum,jmod)= 1.0_dp
                 END IF

                 !--- Set minimum radius and maximum radius:

                 zr1 = crdiv(jmod)
                 zr2 = crdiv(jmod+1)
                 IF (ldconc_m7) THEN

                 !--- Radius of average mass for a lognormal distribution

                   zdm = EXP((LOG(zr1)+LOG(zr2))/2.0_dp)*cmr2ram(jmod)

                 !--- Average mass contained in the mode

                   zttnj = ztotmass(jl,jk)/paernl(jl,jk,jmod)

                 !--- Average number of sulfate molecules or
                 !    equivalent for mixed modes,
                 !    for a particle with radius zdm

                   zavnj= zdm*zdm*zdm                              &
                        * pi*avo*species(spec_idx_h2so4)%density   &
                        / species(spec_idx_h2so4)%molmass/0.75_dp

                 !--- If the average mass contained in the mode is larger
                 !    than the average mass
                 !    for a particle with radius zdm, the transfer of number
                 !    and mass is done, else there is no transfer

                   IF (zttnj > zavnj .AND. jnum == jmod) THEN

                   !--- Mass remaining in the mode
                     zmrj=zfconn(jl,jk,jnum,jmod)*paernl(jl,jk,jmod)*zavnj
!                   zmrj=zfconn(jl,jk,jnum,jmod)*paernl(jl,jk,jmod)*zttnj

                   !--- Mass transferred

                     zmt=ztotmass(jl,jk)-zmrj

                   !--- Numbers transferred

                     znt=(1.0_dp-zfconn(jl,jk,jnum,jmod))*paernl(jl,jk,jmod)

                   !--- Average mass of particles transferred

                     IF(znt > zeps) THEN
                       zavmt=zmt/znt
                     ELSE
                       zavmt=zero
                     END IF

                   !--- Average mass of particles of radius zr2

                     zmcr = (zr2*cmr2ram(jmod)) &
                          * (zr2*cmr2ram(jmod)) &
                          * (zr2*cmr2ram(jmod)) &
                          * pi * avo            &
                          * species(spec_idx_h2so4)%density    &
                          / species(spec_idx_h2so4)%molmass /0.75_dp

                   !--- If the average mass of particle transferred is
                   !    smaller than the average mass
                   !    of particles with radius zr2 then reduce
                   !    the particles transferred
                   !    so that zavmt=zmcr, else calculate the mass
                   !    fraction transferred zfconmj

                     IF (zavmt >= zmcr) THEN
                       zfconmj=zmrj/ztotmass(jl,jk)
                     ELSE
                       zntnew = znt/(1.0_dp + (zmcr-zavmt)/(zavmt-zavnj))
                       zmtnew = zntnew*zmcr
                       zfconmj = 1.0_dp - zmtnew/ztotmass(jl,jk)
                       zfconn(jl,jk,jnum,jmod) = 1.0_dp - zntnew &
                                               / paernl(jl,jk,jmod)
                     END IF

                     zfconm(jl,jk,jnum,jmod)=zfconmj
                   ELSE
                     zfconn(jl,jk,jnum,jmod)=1.0_dp
                     zfconm(jl,jk,jnum,jmod)=1.0_dp
                   END IF

                 ! mz_ht_20090122+
                 ELSE ! ldconc_m7
                   ! calculate the transferred mass directly
                   ! this has the advantage that more mass than numbers
                   ! are transferred, i.e. the bigger particles of the mode
                   ! are shifted into the next mode

                 !--- Mass Median Radius):
                 zarg1 = prad (jl,jk,jmod) * cmr2mmr(jmod)

                 !--- Threshold Radius):
                 zarg4 = crdiv(jnum+1)
                 !--- Transfer to logarithmic scale:
                 zdpmm=LOG(zarg1)
                 zdpbm=LOG(zarg4)

                 zcongm=(zdpbm-zdpmm)/zrcsig

                 CALL gmxe_cumnor(zcongm,zfconm(jl,jk,jnum,jmod),zdummy)
!!$! current volume in that mode
!!$                   ztotvol = paernl(jl,jk,jmod) * prad(jl,jk,jmod) &
!!$                                                * prad(jl,jk,jmod) &
!!$                                                * prad(jl,jk,jmod)
!!$! remaining volume
!!$                   znewvol = ztotvol * zfconn(jl,jk,jnum,jmod)
!!$! maximum volume of remaining particles
!!$                   znewvol = MIN(znewvol, paernl(jl,jk,jmod)  &
!!$                           * zfconn(jl,jk,jnum,jmod)          &
!!$                           * zr2 * zr2 * zr2 )
!!$! minimum transferred volume
!!$                   znewvol = paernl(jl,jk,jmod)                  &
!!$                           * ( 1._dp - zfconn(jl,jk,jnum,jmod) ) &
!!$                           * zr2 * zr2 * zr2
!!$! remaining volume
!!$                   znewvol = MAX(0._dp, ztotvol - znewvol)
!!$! remaining mass factor
!!$                   If (ztotvol > 0._dp)        &
!!$                     zfconm(jl,jk,jnum,jmod) = &
!!$                     znewvol / ztotvol

                 ENDIF ! ldconc_m7
                 ! mz_ht_20090122-
               ELSE
                 zfconn(jl,jk,jnum,jmod)=1.0_dp
                 zfconm(jl,jk,jnum,jmod)=1.0_dp
              END IF
           END DO
        END DO
     END DO
  END DO

!BS-21032007+
  IF (lmass_diag) call sum_mass3(kproma,klev,paerml,zmass_pre2)
  zfconn(:,:,nsoluble,1:nsoluble) = 1.0_dp
  zfconm(:,:,nsoluble,1:nsoluble) = 1.0_dp

  IF (lmass_diag) paerml_pre = paerml


  DO jk=1,klev


    ! mz_ht_20101020+
    DO jnum=nsoluble,2,-1
      DO jmod=1,nsoluble
        DO jl=1,kproma
          zfconn(jl,jk,jnum,jmod) = zfconn(jl,jk,jnum,jmod) - &
                                    zfconn(jl,jk,jnum-1,jmod)
          zfconm(jl,jk,jnum,jmod) = zfconm(jl,jk,jnum,jmod) - &
                                    zfconm(jl,jk,jnum-1,jmod)
        END DO
      END DO
    END DO
  END DO


!!$  print*, "zfconn 1:", zfconn(1,31,:,1)
!!$  print*, "zfconn 2:", zfconn(1,31,:,2)
!!$  print*, "zfconn 3:", zfconn(1,31,:,3)
!!$  print*, "zfconn 4:", zfconn(1,31,:,4)

  DO jk=1,klev
    sumn_vec(:,:)   = 0._dp
    summ2_vec(:,:)  = 0._dp
    summ_vec(:,:,:) = 0._dp

    DO jmod=1,nsoluble

! transfer of numbers of the soluble modes

      !--- Particle numbers:
      DO jnum=1,nsoluble
        DO jl=1,kproma
          sumn_vec(jl,jmod) = sumn_vec(jl,jmod) + &
            paernl(jl,jk,jnum) * zfconn(jl,jk,jmod,jnum)
!!$          if ((jl == 1) .and. (jk == 31) ) THEN
!!$            print*, "in loop: ", jmod, jnum, sumn_vec(jl,jmod), &
!!$              paernl(jl,jk,jnum), zfconn(jl,jk,jmod,jnum)
!!$          END if
        ENDDO
      ENDDO

    END DO
    DO jmod=1,nsoluble
      DO jl=1,kproma
        paernl(jl,jk,jmod)=sumn_vec(jl,jmod)
      ENDDO



        DO jc=1,spec_number
          idx1 = species(jc)%aermlidx(jmod)
          if ((idx1 == 0) .OR. (idx1 == species(jc)%aernlidx(jmod))) cycle
!!$!          DO jnum=jmod,nsoluble
!!$          DO jnum=jmod,jmod+1
!!$            idx2 = species(jc)%aermlidx(jnum)
!!$            if (idx2 == 0 .OR. (idx2 == species(jc)%aernlidx(jnum))) cycle
!!$            DO jl=1,kproma
!!$              zarg1 = paerml(jl,jk,idx1) * zfconm(jl,jk,jnum,jmod)
!!$              summ_vec(jl,jnum,idx2) = summ_vec(jl,jnum,idx2) + zarg1
!!$
!!$
!!$            ENDDO
!!$          END DO


!!$          jnum = jmod
!!$          idx2 = species(jc)%aermlidx(jnum)
!!$          if (idx2 == 0 .OR. (idx2 == species(jc)%aernlidx(jnum))) cycle
!!$          DO jl=1,kproma
!!$            zarg1 = paerml(jl,jk,idx1) * zfconm(jl,jk,jnum,jmod)
!!$            summ_vec(jl,jnum,idx2) = summ_vec(jl,jnum,idx2) + zarg1
!!$          ENDDO
!!$          jnum = jmod + 1
!!$          idx2 = species(jc)%aermlidx(jnum)
!!$          if (idx2 == 0 .OR. (idx2 == species(jc)%aernlidx(jnum))) cycle
!!$          DO jl=1,kproma
!!$            zarg1 = paerml(jl,jk,idx1) * SUM(zfconm(jl,jk,jnum:nsoluble,jmod))
!!$            summ_vec(jl,jnum,idx2) = summ_vec(jl,jnum,idx2) + zarg1
!!$          END DO
!!$          if (jc==5 .AND. jk==31) &
!!$            print*, "zfconm2: ", jmod, zfconm(1,31,jmod,jmod),SUM(zfconm(1,31,jnum:nsoluble,jmod)),&
!!$            zfconm(1,31,jnum,jmod), zfconm(1,31,jnum+1,jmod), "aerml: ",&
!!$            paerml(1,31,idx1), summ_vec(jl,jnum,idx2), paerml(1,31,idx1) * zfconm(1,31,jmod,jmod), &
!!$            paerml(1,31,idx1) * SUM(zfconm(1,31,jnum:nsoluble,jmod))

        END DO
      END DO
!!$      DO jmod=1,nsoluble
!!$        DO jc=1,spec_number
!!$          idx1 = species(jc)%aermlidx(jmod)
!!$          if ((idx1 == 0) .OR. (idx1 == species(jc)%aernlidx(jmod))) cycle
!!$ ! --- Particle mass
!!$          DO jl=1,kproma
!!$            paerml(jl,jk,idx1) = summ_vec(jl,jmod,idx1)
!!$          ENDDO
!!$        END DO

      jmod = nsoluble
      DO jc=1,spec_number
        idx1 = species(jc)%aermlidx(jmod)
        if (idx1 == 0 .OR. (idx1 == species(jc)%aernlidx(jmod))) cycle
        DO jl=1,kproma
          summ2_vec(jl,idx1) = paerml(jl,jk,idx1)
        END DO
      END DO

      DO jmod=1,nsoluble-1
        jnum = jmod + 1
        DO jc=1,spec_number
          idx1 = species(jc)%aermlidx(jmod)
          if (idx1 == 0 .OR. (idx1 == species(jc)%aernlidx(jmod))) cycle
          idx2 = species(jc)%aermlidx(jnum)
          if (idx2 == 0 .OR. (idx2 == species(jc)%aernlidx(jnum))) cycle
          DO jl=1,kproma

            summ2_vec(jl,idx1) =  summ2_vec(jl,idx1) + &
                                  paerml(jl,jk,idx1) * zfconm(jl,jk,jmod,jmod)
            summ2_vec(jl,idx2) =  summ2_vec(jl,idx2) + &
                                  paerml(jl,jk,idx1) * SUM(zfconm(jl,jk,jnum:nsoluble,jmod))

          END DO
        END DO
      END DO

      DO jmod=1,nsoluble
        DO jc=1,spec_number
          idx1 = species(jc)%aermlidx(jmod)
          if (idx1 == 0 .OR. (idx1 == species(jc)%aernlidx(jmod))) cycle
          paerml(1:kproma,jk,idx1) = summ2_vec(1:kproma,idx1)
        END DO

      END DO  ! jmod = 1,nsoluble

    END DO ! klev

    ! mz_ht_20100827+
    ! transfer of shrinked particles into the smaller modes
    ! only one mode per time step possible
    ! If mean radius < mode boundary then shift as many particles into the next
    ! smaller mode as are there according to the modal distribution
    ! The mass transferred should be as big that the remaining particles are
    ! above the mode boundary, i.e. that the new mean radius is 1% above the
    ! mode boundary.
    ! To calculate the particle number fraction the cumnor function is used.
    ! To determine the mass fraction the required mass remaining in the mode is
    ! calculated.
    ! mass_rest = volume_rest / aerosol density
    ! volume_rest = (mode limit radius(+1%) )^3 * 4/3 pi * number_rest
    ! mass_fract = mass_rest / total_mass
    ! total_mass = total_volume / aerosol density
    ! total_volume = radius^3 * 4/3 pi * number
    ! => mass_fract = volme_rest / total_volume
    !               = ( (mode limit radius (+1%))^3 * number_rest ) / &
    !                 ( radius^3 * number )
    DO jmod=nsoluble,2,-1
      DO jk=1,klev
        DO jl=1,kproma
          IF (paernl(jl,jk,jmod) > cmin_aernl &
               .AND. prad (jl,jk,jmod) > ZERO) THEN
            IF (prad(jl,jk,jmod) < crdiv(jmod)) THEN
              zarg2=prad(jl,jk,jmod)
              zdpcm=log(zarg2)
              zdpam=log(crdiv(jmod))
              zdpmm=log(zarg2 * cmr2mmr(jmod))
              zcongn=(zdpcm-zdpam)/sigma_ln(jmod)
              zcongm=(zdpmm-zdpam)/sigma_ln(jmod)

              call gmXe_cumnor(zcongn,zdummy,zfconn(jl,jk,jmod,jmod))
              call gmXe_cumnor(zcongm,zdummy,zfconm(jl,jk,jmod,jmod))
!!$!              print*, jl,jk,jmod, prad(jl,jk,jmod), crdiv(jmod), &
!!$!                   zfconm(jl,jk,jmod,jmod),&
!!$!                   zfconn(jl,jk,jmod,jmod), &
!!$!                   prad(jl,jk,jmod) * cmr2mmr(jmod), zdummy

              paernl(jl,jk,jmod-1) = paernl(jl,jk,jmod-1) &
                                   + paernl(jl,jk,jmod)   &
                                   * zfconn(jl,jk,jmod,jmod)

              ztotvol = paernl(jl,jk,jmod) * ( prad(jl,jk,jmod)   &
                                             * prad(jl,jk,jmod)   &
                                             * prad(jl,jk,jmod)  )

              paernl(jl,jk,jmod)   = paernl(jl,jk,jmod)   &
                                   * (1._dp - zfconn(jl,jk,jmod,jmod))

              znewvol = paernl(jl,jk,jmod) * (sigma(jmod)*crdiv(jmod)) &
                                           * (sigma(jmod)*crdiv(jmod)) &
                                           * (sigma(jmod)*crdiv(jmod))

              zfconm(jl,jk,jmod,jmod) = 1._dp - MIN(1._dp, znewvol / ztotvol)


              DO jc=1,spec_number
                idx1 = species(jc)%aermlidx(jmod)
                if (idx1 == 0 .OR. idx1 == species(jc)%aernlidx(jmod)) cycle
                idx2 = species(jc)%aermlidx(jmod-1)
                if (idx2 == 0 .OR. idx2 == species(jc)%aernlidx(jmod-1)) cycle
                paerml(jl,jk,idx2) = paerml(jl,jk,idx2) &
                                   + paerml(jl,jk,idx1) &
                                   * zfconm(jl,jk,jmod,jmod)
                paerml(jl,jk,idx1) = paerml(jl,jk,idx1) &
                                   * (1._dp - zfconm(jl,jk,jmod,jmod))
              ENDDO
            END IF
          END IF
        END DO
      END DO
    END DO


    ! allow this shrinking due to unequal loss processes of mass and numbers
    ! also for the hydrophobic modes
    ! As more mass is lost by wet/dry deposition or sedimentation than numbers, the remaining particles
    ! are smaller. Hence en effective shrinking of the particles occurs; however, the particles should
    ! not become smaller than the mode boundary. If this occurs a transfer into the next smaller mode
    ! should be performed.
    ! Note: this is only a discretisation issue and has no physical meaning. However, it makes a difference
    ! if the prescribed aerosol modes have a different width (sigma).

    DO jmod=nmod,nsoluble+2,-1
      kmod = jmod - nsoluble + 1
      DO jk=1,klev
        DO jl=1,kproma
          IF (paernl(jl,jk,jmod) > cmin_aernl &
               .AND. prad (jl,jk,jmod) > ZERO) THEN
            IF (prad(jl,jk,jmod) < crdiv(kmod) ) THEN
              zarg2=prad(jl,jk,jmod)
              zdpcm=log(zarg2)
              zdpam=log(crdiv(kmod))
              zdpmm=log(zarg2 * cmr2mmr(jmod))
              zcongn=(zdpcm-zdpam)/sigma_ln(jmod)
              zcongm=(zdpmm-zdpam)/sigma_ln(jmod)
              call gmXe_cumnor(zcongn,zdummy,zfconn(jl,jk,jmod,jmod))
              call gmXe_cumnor(zcongm,zdummy,zfconm(jl,jk,jmod,jmod))
!!$!              print*, jl,jk,jmod, prad(jl,jk,jmod), crdiv(jmod), &
!!$!                   zfconm(jl,jk,jmod,jmod),&
!!$!                   zfconn(jl,jk,jmod,jmod), &
!!$!                   prad(jl,jk,jmod) * cmr2mmr(jmod), zdummy

              paernl(jl,jk,jmod-1) = paernl(jl,jk,jmod-1) &
                                   + paernl(jl,jk,jmod)   &
                                   * zfconn(jl,jk,jmod,jmod)

              ztotvol = paernl(jl,jk,jmod) * ( prad(jl,jk,jmod)   &
                                             * prad(jl,jk,jmod)   &
                                             * prad(jl,jk,jmod)  )

              paernl(jl,jk,jmod)   = paernl(jl,jk,jmod)   &
                                   * (1._dp - zfconn(jl,jk,jmod,jmod))

              znewvol = paernl(jl,jk,jmod) * (sigma(jmod)*crdiv(kmod)) &
                                           * (sigma(jmod)*crdiv(kmod)) &
                                           * (sigma(jmod)*crdiv(kmod))

              zfconm(jl,jk,jmod,jmod) = 1._dp - MIN(1._dp, znewvol / ztotvol)
              DO jc=1,spec_number
                idx1 = species(jc)%aermlidx(jmod)
                if (idx1 == 0 .OR. idx1 == species(jc)%aernlidx(jmod)) cycle
                idx2 = species(jc)%aermlidx(jmod-1)
                if (idx2 == 0 .OR. idx2 == species(jc)%aernlidx(jmod-1)) cycle
                paerml(jl,jk,idx2) = paerml(jl,jk,idx2) &
                                   + paerml(jl,jk,idx1) &
                                   * zfconm(jl,jk,jmod,jmod)
                paerml(jl,jk,idx1) = paerml(jl,jk,idx1) &
                                   * (1._dp - zfconm(jl,jk,jmod,jmod))
              ENDDO
            END IF
          END IF
        END DO
      END DO
    END DO


    ! mz_ht_20100930-
    IF (lmass_diag) THEN
      paerml_post = paerml

      DO jk=1,klev
        DO jl=1,kproma
          pre(:)  = 0._dp
          post(:) = 0._dp
          DO jmod=1,nsoluble
            DO jc=1,spec_number
              idx1 = species(jc)%aermlidx(jmod)
              if (idx1 == 0 .OR. idx1 == species(jc)%aernlidx(jmod)) cycle
              pre(jc)  = pre(jc) + paerml_pre(jl,jk,idx1)
              post(jc) = post(jc) + paerml_post(jl,jk,idx1)
            END DO
          END DO
          DO jc=1,spec_number
            IF (ABS( pre(jc) - post(jc) ) > (10._dp * SPACING(post(jc))) ) THEN
              print*, "change specieswise: ", jc, pre(jc), post(jc),       &
                jl, jk,ABS( pre(jc) - post(jc))
              DO jmod=1,nsoluble
                idx1 = species(jc)%aermlidx(jmod)
                if (idx1 == 0 .OR. idx1 == species(jc)%aernlidx(jmod)) cycle
                print*, "per mode: ", jmod, jc, jl, jk, idx1, &
                  paerml_pre(jl,jk,idx1), paerml_post(jl, jk, idx1), &
                  zfconm(jl,jk,jmod,jmod), SUM(zfconm(jl,jk,jmod+1:nsoluble,jmod))
              END DO
            END IF
          END DO
        END DO
      END DO
    END IF


    IF (lmass_diag) THEN
      call sum_mass3(kproma,klev,paerml,zmass_post2)
      DO jk=1,klev
        DO jl=1,kproma
          IF (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) >                &
            (10._dp * SPACING(zmass_post2(jl,jk))) ) THEN
            print*, "mass changes itemwise; box: ",jl,jk," Error: (%)",  &
              (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) /               &
              zmass_post2(jl,jk)*100._dp),                               &
              zmass_pre2(jl,jk), zmass_post2(jl,jk)
          END IF
        END DO
     END DO
   END IF

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
CONTAINS
SUBROUTINE gmxe_cumnor ( arg, RESULT, ccum )
  !
  !*****************************************************************************
  !
  !! CUMNOR computes the cumulative normal distribution.
  !
  !
  !     the integral from -infinity to x of
  !          (1/sqrt(2*pi)) exp(-u*u/2) du
  !
  !  Author:
  !  -------
  !  Original source:
  !
  !    W. J. Cody    Mathematics and Computer Science Division
  !                  Argonne National Laboratory
  !                  Argonne, IL 60439
  !
  !    DCDFLIB is attributed to Barry Brown, James Lovato, and Kathy Russell
  !            bwb@odin.mda.uth.tmc.edu.
  !
  !    Adopted to GMXe:
  !
  !    Philip Stier  (MPI-MET)                    2001
  !
  !
  !  Reference:
  !  ----------
  !
  !    W D Cody,
  !    "ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of Special
  !    Function Routines and Test Drivers"
  !    ACM Transactions on Mathematical Software,
  !    Volume 19, 1993, pages 22-32.
  !
  !  Parameters:
  !
  !     ARG --> Upper limit of integration.
  !                                        X is REAL(dp)
  !
  !     RESULT <-- Cumulative normal distribution.
  !                                        RESULT is REAL(dp)
  !
  !     CCUM <-- Complement of Cumulative normal distribution.
  !                                        CCUM is REAL(dp)
  !
  !
  ! Original Comments:
  !
  !
  ! This function evaluates the normal distribution function:
  !
  !                              / x
  !                     1       |       -t*t/2
  !          P(x) = ----------- |      e       dt
  !                 sqrt(2 pi)  |
  !                             /-oo
  !
  !   The main computation evaluates near-minimax approximations
  !   derived from those in "Rational Chebyshev approximations for
  !   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
  !   This transportable program uses rational functions that
  !   theoretically approximate the normal distribution function to
  !   at least 18 significant decimal digits.  The accuracy achieved
  !   depends on the arithmetic system, the compiler, the intrinsic
  !   functions, and proper selection of the machine-dependent
  !   constants.
  !
  !  Explanation of machine-dependent constants.
  !
  !   MIN   = smallest machine representable number.
  !
  !   EPS   = argument below which anorm(x) may be represented by
  !           0.5  and above which  x*x  will not underflow.
  !           A conservative value is the largest machine number X
  !           such that   1.0 + X = 1.0   to machine precision.
  !
  !  Error returns
  !
  !  The program returns  ANORM = 0     for  ARG .LE. XLOW.
  !
  !  Author:
  !
  !    W. J. Cody
  !    Mathematics and Computer Science Division
  !    Argonne National Laboratory
  !    Argonne, IL 60439
  !
  !  Latest modification: March 15, 1992
  !
  REAL(dp), PARAMETER, DIMENSION ( 5 ) :: a = (/ &
       2.2352520354606839287d00, &
       1.6102823106855587881d02, &
       1.0676894854603709582d03, &
       1.8154981253343561249d04, &
       6.5682337918207449113d-2 /)
  REAL(dp) arg
  REAL(dp), PARAMETER, DIMENSION ( 4 ) :: b = (/ &
       4.7202581904688241870d01, &
       9.7609855173777669322d02, &
       1.0260932208618978205d04, &
       4.5507789335026729956d04 /)
  REAL(dp), PARAMETER, DIMENSION ( 9 ) :: c = (/ &
       3.9894151208813466764d-1, &
       8.8831497943883759412d00, &
       9.3506656132177855979d01, &
       5.9727027639480026226d02, &
       2.4945375852903726711d03, &
       6.8481904505362823326d03, &
       1.1602651437647350124d04, &
       9.8427148383839780218d03, &
       1.0765576773720192317d-8 /)
  REAL(dp) ccum
  REAL(dp), PARAMETER, DIMENSION ( 8 ) :: d = (/ &
       2.2266688044328115691d01, &
       2.3538790178262499861d02, &
       1.5193775994075548050d03, &
       6.4855582982667607550d03, &
       1.8615571640885098091d04, &
       3.4900952721145977266d04, &
       3.8912003286093271411d04, &
       1.9685429676859990727d04 /)
  REAL(dp) del
  REAL(dp) eps
  INTEGER i
  REAL(dp) min
  REAL(dp), PARAMETER, DIMENSION ( 6 ) :: p = (/ &
       2.1589853405795699d-1, &
       1.274011611602473639d-1, &
       2.2235277870649807d-2, &
       1.421619193227893466d-3, &
       2.9112874951168792d-5, &
       2.307344176494017303d-2 /)
  REAL(dp), PARAMETER, DIMENSION ( 5 ) :: q = (/ &
       1.28426009614491121d00, &
       4.68238212480865118d-1, &
       6.59881378689285515d-2, &
       3.78239633202758244d-3, &
       7.29751555083966205d-5 /)
  REAL(dp) RESULT
  REAL(dp), PARAMETER :: root32 = 5.656854248d0
  REAL(dp), PARAMETER :: sixten = 16.0
  REAL(dp) temp
  REAL(dp), PARAMETER :: sqrpi = 3.9894228040143267794d-1
  REAL(dp), PARAMETER :: thrsh = 0.66291d0
  REAL(dp) x
  REAL(dp) xden
  REAL(dp) xnum
  REAL(dp) y
  REAL(dp) xsq
  !
  !  Machine dependent constants
  !
  eps = EPSILON ( 1.0d0 ) * 0.5d0
  !
  !@@@ Simplified calculation of the smallest machine representable number
  !    (Higher accuracy than needed!)
  !
  !@@@ min = dpmpar(2)

  min = epsilon ( 1.0D0)

  x = arg
  y = ABS ( x )

  IF ( y <= thrsh ) THEN
     !
     !  Evaluate  anorm  for  |X| <= 0.66291
     !
     IF ( y > eps ) THEN
        xsq = x * x
     ELSE
        xsq = zero
     END IF

     xnum = a(5) * xsq
     xden = xsq
     DO i = 1, 3
        xnum = ( xnum + a(i) ) * xsq
        xden = ( xden + b(i) ) * xsq
     END DO
     RESULT = x * ( xnum + a(4) ) / ( xden + b(4) )
     temp = RESULT
     RESULT = 0.5 + temp
     ccum = 0.5 - temp
     !
     !  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
     !
  ELSE IF ( y <= root32 ) THEN

     xnum = c(9) * y
     xden = y
!CDIR UNROLL=7
     DO i = 1, 7
        xnum = ( xnum + c(i) ) * y
        xden = ( xden + d(i) ) * y
     END DO
     RESULT = ( xnum + c(8) ) / ( xden + d(8) )
     xsq = AINT ( y * sixten ) / sixten
     del = ( y - xsq ) * ( y + xsq )
     RESULT = EXP(-xsq*xsq*0.5) * EXP(-del*0.5) * RESULT
     ccum = 1.0_dp - RESULT

     IF ( x > zero ) THEN
        temp = RESULT
        RESULT = ccum
        ccum = temp
     END IF
     !
     !  Evaluate  anorm  for |X| > sqrt(32).
     !
  ELSE

     RESULT = zero
     xsq = 1.0_dp / ( x * x )
     xnum = p(6) * xsq
     xden = xsq
     DO i = 1, 4
        xnum = ( xnum + p(i) ) * xsq
        xden = ( xden + q(i) ) * xsq
     END DO

     RESULT = xsq * ( xnum + p(5) ) / ( xden + q(5) )
     RESULT = ( sqrpi - RESULT ) / y
     xsq = AINT ( x * sixten ) / sixten
     del = ( x - xsq ) * ( x + xsq )
     RESULT = EXP ( - xsq * xsq * 0.5 ) * EXP ( - del * 0.5 ) * RESULT
     ccum = 1.0_dp - RESULT

     IF ( x > zero ) THEN
        temp = RESULT
        RESULT = ccum
        ccum = temp
     END IF

  END IF

  IF ( RESULT < min ) THEN
     RESULT = zero
  END IF

  IF ( ccum < min ) THEN
     ccum = zero
  END IF

END SUBROUTINE gmxe_cumnor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
END SUBROUTINE gmxe_dconc
!-----------------------------------------------------------------------------!
  SUBROUTINE gnucl_kulmala(kproma,  klev,            &
                           pso4g,   ptemp,  prhum,   &
                           pbnrate, palpha, pbeta    )


    !  Authors:
    !  --------
    !  P. Stier, MPI-Met, Hamburg,    from the original f77 code
    !                                 in the routine gmxe_nuck        2001-2003
    !  J. Wilson, E. Vignati, JRC/EI, original source                 09/2000
    !
    !  Purpose:
    !  --------
    !  This routine calculates the instananeous nucleation rate
    !  znucrate [molec. cm-3 s-1] from a given gas phase H2SO4 concentration
    !  pso4g [molec. cm-3]. It also calculates the integrated change of
    !  H2SO4 gas phase mass over one timestep due to nucleation
    !  pa4delt(:,:,1) [molec. cm-3] as well as the number of nucleated
    !  particles panew [1] during the timestep.
    !
    !  Interface:
    !  ----------
    !  *gnucl_kulmala* is called from *gmxe_nuck*
    !
    !  Method:
    !  -------
    !  Kulmala et al. (1998) 's formula for binary nucleation is
    !  rewritten to take the form znucrate = exp[zalpha+ln(pso4g)*beta].
    !  Equation numbers are taken from Kulmala et al. (1998).
    !  After the calculation of the nucleation rate znucrate, it is
    !  integrated in 2) analytically over one timestep, i.e.:
    !
    !  Integration of:
    !
    !  znucrate=d(critn*znav)/dt=exp[zalpha + zbeta*ln(znav)]
    !
    !  gives znav(t0+dt, znav(t0)) where znav(t0) is pso4g.
    !  znav is temporarily stored in zso4g_new and confined between
    !  0 and pso4g.
    !  The number of nucleated particles is then calculated by
    !  assuming a fixed critical mass critn of newly nucleated
    !  particles and dividing the total mass of nucleated sulfate
    !  by this critical mass of one nucleated particle.
    !
    !  Externals:
    !  ----------
    !  None

    IMPLICIT NONE

    INTEGER :: kproma, klev

    REAL(dp):: pso4g(kproma,klev),        ptemp(kproma,klev),    &
               prhum(kproma,klev),        pbnrate(kproma,klev),  &
               palpha(kproma,klev),       pbeta(kproma,klev)

    INTEGER :: jl, jk

    REAL(dp):: znwv,        zln_nac,      ztk,         zsupsat,  &
               zpeh2o,      zpeh2so4,     zra,         zxal,     &
               ztkn,        zssn,         zdelta

    !---1) Calculation of the nucleation rate: ----------------------------

    DO jk=1,klev
       DO jl=1,kproma

         IF (pso4g(jl,jk) .GT. 1e-5 .AND. prhum(jl,jk) > ZERO) THEN  !!KP - Threshold from M7

!         IF (pso4g(jl,jk) .GT. cmin_nuclmolec) THEN
!          IF (pso4g(jl,jk) .GT. cmin_nuclmolec .AND. prhum(jl,jk) > ZERO) THEN
             ztk=ptemp(jl,jk)
             zsupsat=prhum(jl,jk)
             !
             !--- 1.1) Restrict t, and rh to limits where the parameterization ok:
             !
             ztkn = MAX(ztk, 220.0_dp)
             zssn = MIN(zsupsat, 0.90_dp)

             !
             !--- 1.2) Equlibrium vapour pressures (Jaeker-Mirabel (1995), JGR):
             !
             !--- H2O equlibrium vapour pressure (Tabata):
             !
             zpeh2o=0.750064*(10.**(8.42926609-1827.17843/          &
                  ztkn-71208.271/ztkn/ztkn))*1333./bk/ztkn
             !
             !--- H2SO4 equlibrium vapour pressure at 360
             !
             zpeh2so4=EXP(-10156./ztkn+16.259)*7.6e2*1333./bk/ztkn
             !
             !--- H2SO4 equlibrium vapour pressure - correction of ayers
             !    by kulmala - currently not used
             !
             !     payers=exp(-10156/360+16.259)*7.6e2
             !     zpeh2so4=exp(log(payers)+10156*(-1./ztkn+1./360.+0.38/(905-360) * &
             !              (1+log(360./ztkn)-360./ztkn)))*1333/bk/ztkn
             !
             !--- 1.3) Relative acidity (0.0 -1.0):
             !
             zra=pso4g(jl,jk)/zpeh2so4
             !
             !--- 1.4) Water vapour molecule concentration [cm-3]:
             !
             znwv=zsupsat*zpeh2o
             !
             !--- 1.5) Factor delta in Eq. 22:
             !
             zdelta=1.0_dp+(ztkn-273.15)/273.15
             !
             !--- 1.6) Molefraction of H2SO4 in the critical cluster
             !         minus the H2SO4(g) term in Eq. 17:
             !
             zxal = 1.2233-0.0154*zra/(zra+zssn)-0.0415*LOG(znwv)+ 0.0016*ztkn
             !
             !--- 1.7) Exponent of the critical cluster (Eq. 18):
             !
             zln_nac = -14.5125+0.1335*ztkn-10.5462*zssn+1958.4*zssn/ztkn
             !
             !--- 1.8) Sum of all terms in Eq. 20 containing H2SO4(g):
             !
             pbeta(jl,jk) = 25.1289 - 4890.8/ztkn + 7643.4*0.0102/ztkn - &
                            2.2479*zdelta*zssn - 1.9712*0.0102*zdelta/zssn
             !
             !--- 1.9) Sum all terms in Eq. 20 not containing H2SO4(g):
             !
             palpha(jl,jk) = zln_nac*(-25.1289 + 4890.8/ztkn + 2.2479*zdelta*zssn) - &
                             1743.3/ztkn + zxal*(7643.4/ztkn - 1.9712*zdelta/zssn)
             !
             !--- 1.10) Nucleation rate [cm-3 s-1] (Kulmala et al., 1998):
             !
             pbnrate(jl,jk) = EXP(palpha(jl,jk)+LOG(pso4g(jl,jk))*pbeta(jl,jk))

          ELSE

             palpha(jl,jk) =zero
             pbeta(jl,jk)  =zero
             pbnrate(jl,jk)=zero

          END IF ! pso4g(jl,jk) > cmin_nuclmolec
       END DO ! kproma
    END DO !klev

  END SUBROUTINE gnucl_kulmala

!===========================================================================!

  SUBROUTINE gnucl_vehkamaeki(kproma,   klev,               &  ! ECHAM5 dimensions
                              ptemp,    prhd,  pmolecH2SO4, &  ! ECHAM5 temperature, relative humidity
                              pxtrnucr, pntot               )  ! nucleation rate, number of molecules in the
                                                               ! critical cluster
    !
    !   Authors:
    !   ---------
    !   C. TIMMRECK, MPI HAMBURG                                             2002
    !
    !   Purpose
    !   ---------
    !   Calculation of classical nucleation rate
    !
    !   calculation of the nucleation rate after Vehkamaeki et al. (2002)
    !   The calculation of the nucrate ZKNH2SO4 is in cm^-3 s^-1
    !   and a coarse approxmation for the first class
    !
    !   Modifications:
    !   --------------
    !   R. Hommel; rewrite in f90, adopted to ECHAM5; MPI HAMBURG;      Dec. 2002
    !   P. Stier; modularisation and optimization;    MPI HAMBURG;      Jan  2003
    !
    !   H2SO4 still fixed to xxx molc/cm3, no sulfur cycle coupling yet
    !
    !   References:
    !   -----------
    !   Vehkamaeki et al. (2002), An improved parameterization for sulfuric
    !      acid/water nucleation rates for tropospheric and stratospheric
    !      conditions, J. Geophys. Res, 107, D22, 4622
    !
    !   Parameters
    !   ----------
    !   prho = prhop_neu in *sam*
    !
    !   prhd = relative humidity in sam_aeroprop & sam_nucl
    !
    !   pxtrnucr = nucleation rate in [1/m3s]
    !
    !----------------------------------------------------

  IMPLICIT NONE

    INTEGER     :: kproma, klev
    INTEGER     :: jk, jl

    !----------------------------------------------------
    !

    REAL(dp)::   ptemp(kproma,klev),       &
                 prhd(kproma,klev),        &
                 pxtrnucr(kproma,klev),    &
                 pmolecH2SO4(kproma,klev), &
                 pntot(kproma,klev)

    !----------------------------------------------------
    ! Local Arrays

!   REAL(dp)::   zrxc(kproma)

    REAL(dp)::   zrhoa, zrh, zt, x, zjnuc, zrc, zxmole, zntot

    REAL(dp):: zsubtotal  !mz_ak_20031218 introduced to avoid compiler
                          ! error by too many continuation lines
    REAL(dp):: zrh_log, zrhoa_log
    REAL(dp):: zrh_log2, zrh_log3, zt2, zt3, zrhoa_log2, zrhoa_log3

    ! mz_ht_20081001+
    pntot(:,:)   =zero
    pxtrnucr(:,:)=zero
    ! mz_ht_20081001-



  DO jk=1, klev
     DO jl=1,kproma

  !----1.) Parameterization of  nucleation rate after Vehkamaeki et al. (2002)

        ! t: temperature in K (190.15-300.15K)
        ! zrh: saturatio ratio of water (0.0001-1)
        ! zrhoa: sulfuric acid concentration in 1/cm3 (10^4-10^11 1/cm3)
        ! jnuc: nucleation rate in 1/cm3s (10^-7-10^10 1/cm3s)
        ! ntot: total number of molecules in the critical cluster (ntot>4)
        ! x: molefraction of H2SO4 in the critical cluster
        ! rc: radius of the critical cluster in nm

        ! Calculate nucleation only for valid thermodynamic conditions:

!       IF( (pmolecH2SO4(jl,jk)>=1.E+4)                      .AND. &
        IF( (pmolecH2SO4(jl,jk) > cmin_nuclmolec)            .AND. &  !!KP Check threshold
            (prhd(jl,jk) >=1.E-4)                            .AND. &
            (ptemp(jl,jk)>=190.15 .OR. ptemp(jl,jk)<=300.15_dp)    ) THEN

           zrhoa=MIN(pmolecH2SO4(jl,jk),1.e11_dp)
           zrh=MIN(prhd(jl,jk),1.0_dp)
           zt=ptemp(jl,jk)
           zt2=zt*zt
           zt3=zt2*zt
           zrh_log=LOG(zrh)
           zrh_log2=zrh_log*zrh_log
           zrh_log3=zrh_log2*zrh_log
           zrhoa_log=LOG(zrhoa)
           zrhoa_log2=zrhoa_log*zrhoa_log
           zrhoa_log3=zrhoa_log2*zrhoa_log


           ! Equation (11) - molefraction of H2SO4 in the critical cluster

        x=0.7409967177282139_dp - 0.002663785665140117*zt   &
          + 0.002010478847383187*LOG(zrh)    &
          - 0.0001832894131464668*zt*LOG(zrh)    &
          + 0.001574072538464286*LOG(zrh)**2        &
          - 0.00001790589121766952*zt*LOG(zrh)**2    &
          + 0.0001844027436573778*LOG(zrh)**3     &
          -  1.503452308794887e-6*zt*LOG(zrh)**3    &
          - 0.003499978417957668*LOG(zrhoa)   &
          + 0.0000504021689382576*zt*LOG(zrhoa)

        zxmole=x !qqq

        ! Equation (12) - nucleation rate in 1/cm3s

        zjnuc =0.1430901615568665_dp + 2.219563673425199*zt -   &
              0.02739106114964264*zt**2 +     &
              0.00007228107239317088*zt**3 + 5.91822263375044/x +     &
              0.1174886643003278*LOG(zrh) + 0.4625315047693772*zt*LOG(zrh) -   &
              0.01180591129059253*zt**2*LOG(zrh) +     &
              0.0000404196487152575*zt**3*LOG(zrh) +    &
              (15.79628615047088*LOG(zrh))/x -     &
              0.215553951893509*LOG(zrh)**2 -    &
              0.0810269192332194*zt*LOG(zrh)**2 +     &
              0.001435808434184642*zt**2*LOG(zrh)**2 -    &
              4.775796947178588e-6*zt**3*LOG(zrh)**2 -     &
              (2.912974063702185*LOG(zrh)**2)/x -   &
              3.588557942822751*LOG(zrh)**3 +     &
              0.04950795302831703*zt*LOG(zrh)**3 -     &
              0.0002138195118737068*zt**2*LOG(zrh)**3 +    &
              3.108005107949533e-7*zt**3*LOG(zrh)**3 -     &
              (0.02933332747098296*LOG(zrh)**3)/x +     &
              1.145983818561277*LOG(zrhoa) -    &
              0.6007956227856778*zt*LOG(zrhoa) +    &
              0.00864244733283759*zt**2*LOG(zrhoa) -    &
              0.00002289467254710888*zt**3*LOG(zrhoa)! -    &


        zjnuc =zjnuc - &
             (8.44984513869014*LOG(zrhoa))/x +    &
              2.158548369286559*LOG(zrh)*LOG(zrhoa) +   &
              0.0808121412840917*zt*LOG(zrh)*LOG(zrhoa) -    &
              0.0004073815255395214*zt**2*LOG(zrh)*LOG(zrhoa) - &
              4.019572560156515e-7*zt**3*LOG(zrh)*LOG(zrhoa) +    &
              (0.7213255852557236*LOG(zrh)*LOG(zrhoa))/x +    &
              1.62409850488771*LOG(zrh)**2*LOG(zrhoa) -    &
              0.01601062035325362*zt*LOG(zrh)**2*LOG(zrhoa) +   &
              0.00003771238979714162*zt**2*LOG(zrh)**2*LOG(zrhoa) +    &
              3.217942606371182e-8*zt**3*LOG(zrh)**2*LOG(zrhoa) -    &
              (0.01132550810022116*LOG(zrh)**2*LOG(zrhoa))/x +    &
              9.71681713056504*LOG(zrhoa)**2 -    &
              0.1150478558347306*zt*LOG(zrhoa)**2 +    &
              0.0001570982486038294*zt**2*LOG(zrhoa)**2 +    &
              4.009144680125015e-7*zt**3*LOG(zrhoa)**2 +    &
              (0.7118597859976135*LOG(zrhoa)**2)/x -    &
              1.056105824379897*LOG(zrh)*LOG(zrhoa)**2 +    &
              0.00903377584628419*zt*LOG(zrh)*LOG(zrhoa)**2 -    &
              0.00001984167387090606*zt**2*LOG(zrh)*LOG(zrhoa)**2 +    &
              2.460478196482179e-8*zt**3*LOG(zrh)*LOG(zrhoa)**2 -    &
              (0.05790872906645181*LOG(zrh)*LOG(zrhoa)**2)/x -    &
              0.1487119673397459*LOG(zrhoa)**3 +    &
              0.002835082097822667*zt*LOG(zrhoa)**3 -    &
              9.24618825471694e-6*zt**2*LOG(zrhoa)**3 +    &
              5.004267665960894e-9*zt**3*LOG(zrhoa)**3 -    &
              (0.01270805101481648*LOG(zrhoa)**3)/x

        zjnuc=EXP(zjnuc)      !   add. Eq. (12) [1/(cm^3s)]


        ! Equation (13) - total number of molecules in the critical cluster

        zntot =-0.002954125078716302_dp - 0.0976834264241286*zt +   &
               0.001024847927067835*zt**2 - 2.186459697726116e-6*zt**3 -    &
               0.1017165718716887/x - 0.002050640345231486*LOG(zrh) -   &
               0.007585041382707174*zt*LOG(zrh) +    &
               0.0001926539658089536*zt**2*LOG(zrh) -   &
               6.70429719683894e-7*zt**3*LOG(zrh) -    &
               (0.2557744774673163*LOG(zrh))/x +   &
               0.003223076552477191*LOG(zrh)**2 +   &
               0.000852636632240633*zt*LOG(zrh)**2 -    &
               0.00001547571354871789*zt**2*LOG(zrh)**2 +   &
               5.666608424980593e-8*zt**3*LOG(zrh)**2 +    &
               (0.03384437400744206*LOG(zrh)**2)/x +   &
               0.04743226764572505*LOG(zrh)**3 -    &
               0.0006251042204583412*zt*LOG(zrh)**3 +   &
               2.650663328519478e-6*zt**2*LOG(zrh)**3 -    &
               3.674710848763778e-9*zt**3*LOG(zrh)**3 -   &
               (0.0002672510825259393*LOG(zrh)**3)/x -    &
               0.01252108546759328*LOG(zrhoa) !+   &

        zntot =zntot + &
               0.005806550506277202*zt*LOG(zrhoa) -    &
               0.0001016735312443444*zt**2*LOG(zrhoa) +   &
               2.881946187214505e-7*zt**3*LOG(zrhoa) +    &
               (0.0942243379396279*LOG(zrhoa))/x -   &
               0.0385459592773097*LOG(zrh)*LOG(zrhoa) -   &
               0.0006723156277391984*zt*LOG(zrh)*LOG(zrhoa) + &
               2.602884877659698e-6*zt**2*LOG(zrh)*LOG(zrhoa) +    &
               1.194163699688297e-8*zt**3*LOG(zrh)*LOG(zrhoa) -   &
               (0.00851515345806281*LOG(zrh)*LOG(zrhoa))/x -    &
               0.01837488495738111*LOG(zrh)**2*LOG(zrhoa) +   &
               0.0001720723574407498*zt*LOG(zrh)**2*LOG(zrhoa) -   &
               3.717657974086814e-7*zt**2*LOG(zrh)**2*LOG(zrhoa) -    &
               5.148746022615196e-10*zt**3*LOG(zrh)**2*LOG(zrhoa) +    &
               (0.0002686602132926594*LOG(zrh)**2*LOG(zrhoa))/x -   &
               0.06199739728812199*LOG(zrhoa)**2 +    &
               0.000906958053583576*zt*LOG(zrhoa)**2 -   &
               9.11727926129757e-7*zt**2*LOG(zrhoa)**2 -    &
               5.367963396508457e-9*zt**3*LOG(zrhoa)**2 -   &
               (0.007742343393937707*LOG(zrhoa)**2)/x +    &
               0.0121827103101659*LOG(zrh)*LOG(zrhoa)**2 -   &
               0.0001066499571188091*zt*LOG(zrh)*LOG(zrhoa)**2 +    &
               2.534598655067518e-7*zt**2*LOG(zrh)*LOG(zrhoa)**2 -    &
               3.635186504599571e-10*zt**3*LOG(zrh)*LOG(zrhoa)**2 +    &
               (0.0006100650851863252*LOG(zrh)*LOG(zrhoa)**2)/x +   &
               0.0003201836700403512*LOG(zrhoa)**3 -    &
               0.0000174761713262546*zt*LOG(zrhoa)**3 +   &
               6.065037668052182e-8*zt**2*LOG(zrhoa)**3 -    &
               1.421771723004557e-11*zt**3*LOG(zrhoa)**3 +   &
               (0.0001357509859501723*LOG(zrhoa)**3)/x


           pntot(jl,jk)=EXP(zntot)  !  add. Eq. (13)


           ! Equation (14) - radius of the critical cluster in nm

           zrc=EXP(-1.6524245+0.42316402*x+0.33466487*LOG(pntot(jl,jk)))    ! [nm]

           ! Conversion [nm -> m]

!ueberfl.  zrxc(jl)=zrc*1e-9

           !----1.2) Limiter

           IF (pntot(jl,jk) < 4.0 ) THEN
              IF ( zt < 195.15) THEN
                 !!! zjnuc=1.e5
                 zjnuc=zero !!!KP testing - same as m7
              END IF
           END IF

           IF(zjnuc < 1.e-7) THEN
              zjnuc=zero
           END IF

           ! limitation to 1E+10 [1/cm3s]
           zjnuc=MIN(zjnuc,1.e10_dp)          !mz_kp_20080112 cmin_aernl is far too low a threshold
           !!!!zjnuc=MIN(zjnuc,cmin_aernl)

           pxtrnucr(jl,jk) = zjnuc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!KP Insered from M7  - for testing only
!       ! mz_ak_20041102+
!       ! convert total number of molecules in the critical cluster
!       ! to number of sulfate molecules:
!
!       pntot(jl,jk)=zntot*zxmole
!
!        !print*,'pntot(jl,jk)=',pntot(jl,jk),zntot,zxmole
!        ! mz_ak_20041102-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ELSE ! pmolecH2SO4, ptemp , prhd out of range

           pntot(jl,jk)   =zero
           pxtrnucr(jl,jk)=zero

        END IF

     END DO ! kproma
  END DO ! klev

END SUBROUTINE gnucl_vehkamaeki

!_________________________________________________________________________________________________________________________________


SUBROUTINE init_paerosol

    !-------------------------------------------------------------
    ! GMXe DIAGNOSTIC OUTPUT
    !-------------------------------------------------------------

    ALLOCATE(cdiagaer (naerodiag,2)) ; cdiagaer = ''

    cdiagaer(iTT,1)        = 'T = air temperature'
    cdiagaer(iTT,2)        = '[K]'
    cdiagaer(iRH,1)        = 'RH = relative humidity'
    cdiagaer(iRH,2)        = '[0-1]'
    cdiagaer(iPX,1)        = 'PX = aerosol water history '
    cdiagaer(iPX,2)        = '[0=solid, else=wet]'
    cdiagaer(iZIONIC,1)    = 'ZIONIC = ionic strength (aq)'
    cdiagaer(iZIONIC,2)    = '[mol/kg]'
    cdiagaer(iRHO,1)       = 'RHO = density (bulk mass)'
    cdiagaer(iRHO,2)       = '[g/m3]'
    cdiagaer(iVOL,1)       = 'VOL = bulk volume'
    cdiagaer(iVOL,2)       = '[ucm3 m-3 (air)]'
    cdiagaer(iPH,1)        = 'PH = aerosol pH'
    cdiagaer(iPH,2)        = '[log H+]'
    cdiagaer(isPM,1)       = 'sPM = total solid  matter'
    cdiagaer(isPM,2)       = '[umol m-3 (air)]'
    cdiagaer(iaPM,1)       = 'aPM = total liquid matter'
    cdiagaer(iaPM,2)       = '[umol m-3 (air)]'
    cdiagaer(iPMt,1)       = 'PMt = total PM (liquids & solids)'
    cdiagaer(iPMt,2)       = '[ug   m-3 (air)]'
    cdiagaer(iPMs,1)       = 'PMs = total PM (solids treated by EQSAM)'
    cdiagaer(iPMs,2)       = '[ug   m-3 (air)]'
    cdiagaer(irhdcr_sum,1) = 'rhdcr_sum = crystallization RH of mixed solution'
    cdiagaer(irhdcr_sum,2) = '[-]'
    cdiagaer(irhdms_sum,1) = 'rhdms_sum = deliquescence   RH of mixed solution'
    cdiagaer(irhdms_sum,2) = '[-]'
    cdiagaer(iRHcr,1)      = 'RHcr = lowest crystallization RH in mixed solution'
    cdiagaer(iRHcr,2)      = '[-]'
    cdiagaer(iRHD,1)       = 'RHD  = lowest deliquescence   RH in mixed solution'
    cdiagaer(iRHD,2)       = '[-]'
    cdiagaer(iWH2O,1)      = 'WH2O = aerosol Water (aq)'
    cdiagaer(iWH2O,2)      = '[ug m-3 (air)]'
    cdiagaer(iGF,1)        = 'GF = hygroscopic growth factor'
    cdiagaer(iGF,2)        = '[-]'
    cdiagaer(iTice,1)      = 'freezing depression by aerosols'
    cdiagaer(iTice,2)      = '[K]'
!   cdiagaer(iConv,1)      = 'conv fac [umol m-3 (air)] -> [molecules cm-3 (air)]'
!   cdiagaer(iConv,2)      = '[-]'
    cdiagaer(inewph,1)     = 'new aerosol pH'
    cdiagaer(inewph,2)     = '[-]'
    cdiagaer(iKappa,1)     = 'Kappa value'
    cdiagaer(iKappa,2)     = '[-]'
    cdiagaer(iKa_water,1)  = 'Dry Diameter'
    cdiagaer(iKa_water,2)  = '[-]'
    cdiagaer(iSc,1)        = 'Critical Supersaturation'
    cdiagaer(iSc,2)        = '[-]'
    cdiagaer(iCCN2,1)      = 'CCN at 0.2% Supersaturation'
    cdiagaer(iCCN2,2)      = '[-]'
    cdiagaer(iCCN4,1)      = 'CCN at 0.4% Supersaturation'
    cdiagaer(iCCN4,2)      = '[-]'
    cdiagaer(iKappa_insol,1) = 'Kappa Insol value'
    cdiagaer(iKappa_insol,2) = '[-]'
    cdiagaer(iKa_vol,1)   = 'Soluble Volume'
    cdiagaer(iKa_vol,2)   = '[-]'
    cdiagaer(iKa_vol_insol,1) = 'Insoluble + Soluble Volume'
    cdiagaer(iKa_vol_insol,2) = '[-]'
    cdiagaer(ioc_water,1)     = 'aerosol water attached to organic carbon'
    cdiagaer(ioc_water,2)     = '[kg / m^3]'
    cdiagaer(ioc_kappa,1)     = 'organic aerosol Kappa value'
    cdiagaer(ioc_kappa,2)     = '[-]'
    cdiagaer(ifsulf,1)        = 'fraction of sulphate from H2SO4'
    cdiagaer(ifsulf,2)        = '[-]'
!!$    cdiagaer(itest,1)        = 'test quantity'
!!$    cdiagaer(itest,2)        = '[-]'
!!$    cdiagaer(itest2,1)        = 'test quantity 2'
!!$    cdiagaer(itest2,2)        = '[-]'


!   additional diagnostic output - example ...
!!$   cdiagaer(iNCawano,1)   = 'solid ammonium nitrate'
!!$   cdiagaer(iNCawano,2)   = '[umol m-3 (air)]'
!!$   cdiagaer(iNAawano,1)   = 'aqueous ammonium nitrate'
!!$   cdiagaer(iNAawano,2)   = '[umol m-3 (air)]'
!!$   cdiagaer(iNAawasu,1)   = 'aqueous ammonium sulfate'
!!$   cdiagaer(iNAawasu,2)   = '[umol m-3 (air)]'
    !-------------------------------------------------------------
END SUBROUTINE init_paerosol

!_________________________________________________________________________________________________________________________________
!_________________________________________________________________________________________________________________________________
ELEMENTAL FUNCTION power2(abasis,xexpon)
!
IMPLICIT NONE
!
! second order approximation of a**x
!
REAL(dp), INTENT(in) :: abasis, xexpon
REAL(dp)             :: power2
REAL(dp)             :: b, b2, lna, x
REAL(dp), PARAMETER :: rloc0=2._dp,rloc1=1._dp/3._dp, rloc2=1._dp/5._dp
!
b  = (abasis-1._dp)/(abasis+1._dp)
b2 = b*b
lna= rloc0*b*(1._dp+b2*(rloc1+b2*rloc2))
x  = xexpon*lna
power2=1._dp+x+x*x*0.5_dp
!
END FUNCTION power2
!_________________________________________________________________________________________________________________________________
ELEMENTAL FUNCTION power1(abasis,xexpon)
!
IMPLICIT NONE
!
! first order approximation of a**x
!
REAL(dp), INTENT(in) :: abasis, xexpon
REAL(dp)             :: power1
REAL(dp)             :: b, lna
REAL(dp), PARAMETER :: rloc1=2._dp/3._dp
!
b=(abasis-1._dp)/(abasis+1._dp)
lna=rloc1*b*(3._dp+b*b)
power1=1._dp+xexpon*lna
!
END FUNCTION power1
!_________________________________________________________________________________________________________________________________

! mz_ht_20090310+
!===============================================================================
SUBROUTINE Hp_OHm_check(kproma,nlev,paerml,paerosol,ep)


! This subroutine exchanges negative H+ ions by forming corresponding OH- ions
! To compensate for the mass created, the same amount of water is destroyed.
! In case of sufficient aerosol water this is done from the aerosol water
! reservoir, else from the gas phase reservoir.
  USE messy_main_constants_mem,              ONLY: M_H2O, MH, MO


  INTEGER,  INTENT(IN)                            :: kproma,nlev
  INTEGER,  INTENT(IN)                            :: ep   ! entry-point
                                                  ! 1 after thermodynamics
                                                  ! 2 after microphysics
  REAL(dp), DIMENSION(kproma,nlev,0:naertot)      :: paerml
  REAL(dp), DIMENSION(kproma,nlev,0:nmod,0:naero) :: paerosol

  ! Local variables
  REAL(dp) :: zdif
  REAL(dp) :: z_pH, z_hp
  REAL(dp) :: befo, afte, befo_h2og, befo_h2oL,befo_hp,befo_oh,afte_h2og,&
              afte_h2oL,afte_hp,afte_oh
  INTEGER  :: acomp, bcomp, dcomp, ecomp
  INTEGER  :: jm, jk, jl

  acomp = nwh2o(0)
  DO jm=1,nsoluble
    IF ( (ep == 2) .AND. AERCHEM(JM) ) CYCLE
    bcomp = nwh2o(jm)
    dcomp = nhp(jm)
    ecomp = nohm(jm)
    DO jk = 1,nlev
      DO jl = 1,kproma

        befo_h2og = paerml(jl,jk,acomp) * M_H2O
        befo_h2ol = paerml(jl,jk,bcomp) * M_H2O
        befo_hp   = paerml(jl,jk,dcomp) * MH
        befo_oh   = paerml(jl,jk,ecomp) * (MO + MH)
        befo = befo_h2ol + befo_hp + befo_oh + befo_h2og


        ! a) checking for negative H+ to form OH-
        ! only if negative H+ is in that grid box
        IF (paerml(jl,jk,dcomp) < 0._dp) THEN

          ! forming OH (lcomp)
          paerml(jl,jk,ecomp) = paerml(jl,jk,ecomp) - paerml(jl,jk,dcomp)
          ! taking it from the water phases
          IF (paerml(jl,jk,bcomp) > (-1._dp) * paerml(jl,jk,dcomp) ) THEN
            ! a) from the aerosol phase (if sufficient water is available)
            paerml(jl,jk,bcomp) = paerml(jl,jk,bcomp) + paerml(jl,jk,dcomp)
          ELSE
            ! b) from the gas phase (first depleating the little amount of
            !    aerosol water, and then taking the rest from the gas phase
            !    water vapour)
            zdif = paerml(jl,jk,dcomp) + paerml(jl,jk,bcomp)
            paerml(jl,jk,acomp) = paerml(jl,jk,acomp) + zdif
            paerml(jl,jk,bcomp) = 0._dp
          ENDIF
          paerml(jl,jk,dcomp) = 0._dp

        ENDIF      ! negative H+

        ! b) neutralisation of H+ + OH- -> H2O
        IF (paerml(jl,jk,bcomp) > 0._dp) THEN   ! only if there is wet aerosol
          IF ( (paerml(jl,jk,dcomp) > 0._dp) .and. &
               (paerml(jl,jk,ecomp) > 0._dp) ) THEN
            ! amount of water produced by neutralisation
            zdif = MIN( paerml(jl,jk,dcomp),paerml(jl,jk,ecomp) )
            ! reduction of H+
            paerml(jl,jk,dcomp) = paerml(jl,jk,dcomp) - zdif
            ! reduction of OH-
            paerml(jl,jk,ecomp) = paerml(jl,jk,ecomp) - zdif
            ! formation of new aerosol water
            paerml(jl,jk,bcomp) = paerml(jl,jk,bcomp) + zdif
          ENDIF
        ENDIF            ! neutralisation
        ! c) determination of aerosol pH
        ! only if there is wet aerosol
        z_pH = NaN
        IF (paerml(jl,jk,bcomp) > 1.e-20_dp) THEN
          z_pH = 7._dp
          !
          ! a) acidic regime (if H+ exists, OHm does not exist)
          IF (paerml(jl,jk,dcomp) > 1.e-20_dp) THEN
            ! z_hp is H+ in mol H+ / liter (aerosol water)
            z_hp = paerml(jl,jk,dcomp) / (paerml(jl,jk,bcomp) * mwh2o) * &
                   Dw * 1.e3_dp
            ! autodissociation dominates if
            ! H+ by other processes/compounds is lower than
            ! that produced by the autodissocation
            ! calculation of pH (if autodissociation is NOT dominant)
            IF (z_hp > 1.0e-7_dp) z_pH = -LOG10(z_hp)
          ELSE

            ! alkaline regime
            IF (paerml(jl,jk,ecomp) > 1.e-20_dp) THEN
            ! z_hp is OH- in mol OH- / liter (aerosol water)
              z_hp = paerml(jl,jk,ecomp) / (paerml(jl,jk,bcomp) * mwh2o) * &
                     Dw * 1.e3_dp
              ! autodissociation dominates if
              ! OH- by other processes/compounds is lower than
              ! that produced by the autodissocation
              ! calculation of pH (if autodissociation is NOT dominant)
              ! first determining pOH (since alkaline regime)
              ! pOH = -LOG10(z_hp)
              ! pH + pOH = 14
              ! pH = 14 - pOH = 14 - (-LOG10(z_hp)) = 14 + LOG10(z_hp)
              IF (z_hp > 1.0e-7_dp) z_pH = 14._dp + LOG10(z_hp)
            ENDIF
          ENDIF
        ENDIF
        paerosol(jl,jk,jm,inewPH) = z_pH ! aerosol pH [log H+]

        afte_h2og = paerml(jl,jk,acomp) * M_H2O
        afte_h2ol = paerml(jl,jk,bcomp) * M_H2O
        afte_hp   = paerml(jl,jk,dcomp) * MH
        afte_oh   = paerml(jl,jk,ecomp) * (MH+ MO)
        afte = afte_h2ol + afte_hp + afte_oh + afte_h2og

        IF (ABS(afte - befo) > 100._dp * spacing(afte)) &
          print*, ABS(afte - befo), &
          ABS(afte - befo)/befo, jl, jk, jm, afte, befo, &
          "before: ",  befo_h2og, befo_h2ol, befo_hp, befo_oh, &
          "after: ",   afte_h2og, afte_h2ol, afte_hp, afte_oh

      END DO
    END DO
  END DO

END SUBROUTINE Hp_OHm_check
!===============================================================================
SUBROUTINE sum_mass3(kproma,nlev,paerml,pmasssum)

      IMPLICIT NONE

      INTEGER     :: jl, jk, jm, jc, it, kproma,nlev
      REAL(dp)    :: pmasssum(kproma,nlev), zscale
      REAL(dp), DIMENSION(kproma,nlev,0:naertot)  :: paerml

! this one tests the conservation of molecules

!!$      pmasssum = ZERO
!!$      DO jm = 0, nmod
!!$        DO jc = 1, mcomp(jm)
!!$          it = icomp(jc,jm)
!!$          zscale = FACp6
!!$          IF(it == nwh2o(jm)) zscale = FACm4 ! scale down water vapor
!!$          DO jk = 1,klev
!!$            DO jl = 1,kproma
!!$              pmasssum = pmasssum + zaerml(jl,jk,it) * zscale
!!$            END DO
!!$          END DO
!!$        END DO
!!$      END DO

! this one tests the conservation of masses

      pmasssum = ZERO
      DO jl=1,kproma
        DO jk=1,nlev
          DO jm = 0, nmod
            DO jc = 1, spec_number
              it = species(jc)%aermlidx(jm)
              IF (it == species(jc)%aernlidx(jm)) CYCLE
              zscale = 1._dp
              pmasssum(jl,jk) = pmasssum(jl,jk) + paerml(jl,jk,it) * &
                                species(jc)%molmass * zscale
            END DO
          END DO
        END DO
      END DO

    END SUBROUTINE sum_mass3
!===============================================================================
    SUBROUTINE SUBGRID_SCALE_RHUM(kproma, klev,                                &
                                  ppress, ppress_half, ptemp, psphum, pvervel, &
                                  pfrac,  prh,         prh_rest )


      USE MESSY_MAIN_CONSTANTS_MEM, ONLY: vtmpc1
      USE MESSY_MAIN_TOOLS,         ONLY: tlucuaw, tlucua, jptlucu1, jptlucu2
      IMPLICIT NONE

      REAL(dp), PARAMETER :: crs     = 0.9_dp
      REAL(dp), PARAMETER :: crt     = 0.7_dp
      REAL(dp), PARAMETER :: nex     = 4._dp
      REAL(dp), PARAMETER :: crhsc   = 0.6_dp
      REAL(dp), PARAMETER :: csatsc  = 0.8_dp

      INTEGER,  INTENT(IN) :: kproma, klev

      REAL(dp), INTENT(IN)    :: ppress(kproma,klev)
      REAL(dp), INTENT(IN)    :: ppress_half(kproma,klev+1)
      REAL(dp), INTENT(IN)    :: ptemp(kproma,klev),prh(kproma,klev)
      REAL(dp), INTENT(IN)    :: psphum(kproma,klev)
      REAL(dp), INTENT(IN)    :: pvervel(kproma,klev)
      REAL(dp), INTENT(INOUT) :: pfrac(kproma,klev), prh_rest(kproma,klev)

      INTEGER  :: jk, jl, it
      REAL(dp) :: zrhc, zsat, zqsm1, zqr

      DO jk=1,klev
        DO jl=1,kproma
!!$          !  need pressure
!!$          zrhc = crt + (crs-crt) * EXP(1._dp - &
!!$                (ppress_half(jl,jk+1) / ppress(jl,jk))**nex)
!!$          zsat=1._dp
!!$          ! need temperature
!!$          it = NINT(ptemp(jl,jk)*FACp3)
!!$          it = MAX(MIN(it,jptlucu2),jptlucu1)
!!$          !  need pressure
!!$!          zqsm1       = MERGE(tlucua(it),tlucuaw(it),lo2)/ppress(jl,jk)
!!$          zqsm1       = tlucuaw(it)/ppress(jl,jk)
!!$          zqsm1       = MIN(zqsm1,0.5_dp)
!!$          zqsm1       = zqsm1/(1.0_dp-vtmpc1*zqsm1)
!!$
!!$!!!$          jb=knvb(jl)
!!$!!!$          lo=(jb.GE.jbmin .AND. jb.LE.jbmax .AND. pvervel(jl,jk).GT.0._dp)
!!$!!!$          lo1=(jk.EQ.jb .OR. jk.EQ.jb+1)
!!$!!!$          IF (lo .AND. lo1) THEN
!!$
!!$          ! need vervel
!!$          IF (pvervel(jl,jk) > 0._dp) THEN
!!$            zsat=csatsc*zrhc
!!$            zrhc=crhsc*zrhc
!!$          END IF
!!$
!!$          zqr=zqsm1*zsat*zrhc
!!$          ! need moisture (kg/kg)
!!$          pfrac(jl,jk)=(psphum(jl,jk)-zqr)/(zqsm1*zsat-zqr)
!!$          pfrac(jl,jk)=MAX(pfrac(jl,jk),0.0_dp)
!!$          pfrac(jl,jk)=MIN(pfrac(jl,jk),1.0_dp)
!!$          pfrac(jl,jk)=1._dp-SQRT(1._dp-pfrac(jl,jk))
!!$!          IF (pxim1(jl,jk)<=ccwmin .AND. pxlm1(jl,jk)<=ccwmin) THEN
!!$!            paclc(jl,jk)=0.0_dp
!!$!          ENDIF
!!$          pfrac(jl,jk) = MAX(MIN(pfrac(jl,jk),1.0_dp),0.0_dp)
!!$        END DO !jl
!!$      END DO  !jk
!      print*, "fractions: " ,pfrac
          IF(pvervel(jl,jk) < 0._dp .and. pRH(jl,jk) > 0.9_dp) THEN
!            pfrac(jl,jk)    = MIN(1._dp,SQRT(-1._dp * pvervel(jl,jk)) * 1._dp)
!            pfrac(jl,jk)    = MIN(1._dp,SQRT(-1._dp * pvervel(jl,jk)) * 3._dp)
            pfrac(jl,jk)    = MIN(0.9999_dp,SQRT(-1._dp * pvervel(jl,jk)) * 6._dp)

            pRH_rest(jl,jk) = (pRH(jl,jk) - pfrac(jl,jk)) / &
                              (1._dp - pfrac(jl,jk))

!            pfrac(jl,jk)     = 1._dp
!            pRH_rest(jl,jk)  = 0._dp


            IF (pRH(jl,jk) > 1._dp) THEN
              pfrac(jl,jk) = 0._dp
              pRH_rest(jl,jk) = pRH(jl,jk)
            ENDIF
          ELSE
            pfrac(jl,jk)    = 0._dp
            pRH_rest(jl,jk) = pRH(jl,jk)
          ENDIF
        ENDDO
      ENDDO

    END SUBROUTINE SUBGRID_SCALE_RHUM
!------------------------------------------------------------------------------



!-----------------------------------------------------------------------------
    SUBROUTINE  gmxe_evap2(kproma,k300, input_h2so4, output_h2so4, &
                        radius, ptemp, ppress, paernl, ztmst)
   IMPLICIT NONE

   INTEGER :: kproma, k300
   REAL(dp), DIMENSION(kproma,k300,0:nsoluble) :: input_h2so4,  &
                                                  output_h2so4
   REAL(dp), DIMENSION(kproma,k300,nsoluble)   :: radius
   REAL(dp), DIMENSION(kproma,k300,nsoluble)   :: paernl
   REAL(dp), DIMENSION(kproma,k300)        :: ptemp,ppress
   REAL(dp)                                :: ztmst


   REAL(dp), DIMENSION(kproma,k300,nmod)   :: zc2
   REAL(dp), DIMENSION(nmod) :: aern
   REAL(dp) :: ptemku, psatku, zpbyone, zde2, zvelb, zxibc, zrwet, zf1
   REAL(dp) :: pd, eva, evah, psatk0, psatve, ptemve, qtemve, psatv0
   INTEGER  :: jl, jk, jm

   aern(1)=1000._dp    !10000
   aern(2)=200._dp   !500
   aern(3)=20._dp   !50
   aern(4)=0.02_dp   !0.01, 0.3 not stable
   aern(5)=500._dp
   aern(6)=50._dp
   aern(7)=0.01_dp

   zc2(:,:,:)=0.0_dp

 ! generalize, this is only for L39MA

   ! vehkamaeki-saturation vapour pressure (JGR, 2002)
   ! p_sat in Pa; TEMP in K
   ! p_sat = 101325 * exp { L + 10156 *
   !         [ 1 / 360.15 - 1/TEMP +
   !         0.38/545 * (1 + ln (360.15/TEMP) - 360.15/TEMP) ] }
   ! fitting parameter L = - 11.695
   ! conversion from Pa to 1 /cm^3
   ! using k_b (in 10^-7 J / K) * TEMP  -> 10^-7 / m^3
   ! results in a factor of 10 from 1 / cm^3
   ! ref at 200K
   ! psatv0 = exp ( -11.695_dp + 10156.0_dp * (0.00277662_dp - 1._dp/200.0_dp &
   !       + 0.000697248_dp * (1._dp + log(1.80075_dp) - 1.80075_dp)))       &
   !       * 101325._dp * 10._dp / (bk * 200.0_dp)
   DO jk=1,k300  !25
     DO jl=1,kproma
       output_h2so4(jl,jk,0:nsoluble) = input_h2so4(jl,jk,0:nsoluble)

!k       ptemku=max(220.0_dp,ptemp(jl,jk))
!k       psatku=exp(16.259_dp-10156.0_dp/ptemku)*1013000._dp &
!k       /(bk*ptemp(jl,jk))         !see old version for Kulmala

       ptemve = MAX(190.0_dp,ptemp(jl,jk))
       qtemve = 360.15_dp/ptemve

       psatve = exp ( -11.695_dp + 10156.0_dp * (0.00277662_dp &
              - 1._dp/ptemve + 0.000697248_dp * (1._dp         &
              + log(qtemve) - qtemve))) * 1013250._dp / (bk * ptemp(jl,jk) )

       ! loop over soluble modes
       DO jm=nsoluble,1,-1
         IF ((output_h2so4(jl,jk,0) < psatve) .AND. &
             (ptemp(jl,jk) > 200.0_dp)) THEN

! about 1 for acc and T=210  (test also with 0.02 for T=200 and 0.002)
!           pd=(psatve-output_h2so4(jl,jk,0))/psatv0*0.0009
           pd = (psatve - output_h2so4(jl,jk,0))
           IF (radius(jl,jk,jm) > 1.e-10) THEN

           !--- Diffusion coefficient (Reference???):
             zpbyone=1000.0_dp / (ppress(jl,jk)/100.0_dp) ! hPa?

             zde2=0.073_dp * zpbyone * (ptemp(jl,jk) / 298.15_dp)**1.5_dp

           !--- Mean molecule velocity (Moore, 1962 (S+P equ. 8.2)):

             zvelb=SQRT(8.0_dp * rerg * ptemp(jl,jk) / pi / 98.08_dp )

           !--- ???Fuchs???

             zxibc=8.0_dp * zde2 / pi / zvelb

           !
           ! Use count median radius:

             zrwet=radius(jl,jk,jm)

           !--- Distance from particle up to which the kinetic regime applies:

             zf1=( (zrwet + zxibc)**3 - &
                 (  zrwet**2 + zxibc**2)**1.5_dp ) / &
                 (  3.0_dp * zrwet * zxibc) - zrwet

             ! mz_ht_20100827+
!!$! from chb
!!$             zc2(jl,jk,jm) = min(1.0_dp, (4.0_dp * pi * zde2 * zrwet ) /     &
!!$                                 ((4.0_dp * zde2) / (zvelb * zrwet ) +       &
!!$                                 (zrwet/(zrwet+zf1)) )*ztmst * pd * aern(jm) )
             zc2(jl,jk,jm) = min(1.0_dp, (4.0_dp * pi * zde2 * zrwet ) /     &
                                 ((4.0_dp * zde2) / (zvelb * zrwet ) +       &
                                 (zrwet/(zrwet+zf1)) ) * paernl(jl,jk,jm) *  &
                                 ztmst * pd )
             ! mz_ht_20100827-

             ! changed to avoid supersaturation
             eva  = MIN ( psatve - output_h2so4(jl,jk,0), &
                          input_h2so4(jl,jk,jm) * zc2(jl,jk,jm) )

             output_h2so4(jl,jk,jm) = input_h2so4(jl,jk,jm) - eva
             output_h2so4(jl,jk,0)  = output_h2so4(jl,jk,0) + eva

           ELSE
             output_h2so4(jl,jk,0)  = output_h2so4(jl,jk,0) + &
                                      output_h2so4(jl,jk,jm)
             output_h2so4(jl,jk,jm) = 0._dp
             paernl(jl,jk,jm)       = 0._dp


           endif   !radius

         endif     !sub - saturation and TEMP
       enddo       ! mode loop
     enddo         ! kproma
   enddo           ! klev


 END SUBROUTINE  gmxe_evap2

!=============================================================================

END MODULE messy_gmxe
