! =====================================================================
!
! DESCRIPTION
!   Fortran 95 module for simulation of polar stratospheric clouds (PSCs)
!   within a global chemistry-climate model. 
!   Contains functions and subroutines for calculation of
!   - PSC occurrence
!   - Solid PSC particle sedimentation
!   - reaction rates for heterogeneous reactions on PSC surfaces
!
! AUTHORS
!   S. Meilinger, MPI for Chemistry, Mainz, Germany
!     (Fortran 90 version for use in ECHAM 5, based on the heterogeneous
!     chemistry in the "Mainz Stratospheric Box Model - Version 6" by
!     K. Carslaw)
!   J. Buchholz (major revision + psc sedimentation)
!   B. Steil (functions mz_psc_surface_liquid and mz_psc_dim_liquid)
!   O. Kirner (kinetic nat parameterisation)
!
! BIBLIOGRAPHY
! K. S. Carslaw, B. Luo, T. Peter: "An analytic expression for the
!   composition of aqueous HNO3-H2SO4 stratospheric aerosols including
!   gas phase removal of HNO3"; GEOPHYSICAL RESEARCH LETTERS, Vol. 14,
!   No. 14, pp. 1877-1880, 1995
! K. S. Carslaw, T. Peter: "Uncertainties in reactive uptake coefficients
!   for solid stratospheric particles. 1. Surface chemistry"; Geophysical
!   Research Letters, Vol. 24, No. 14, pp. 1743-1746, 1997
! K. S. Carslaw, T. Peter, R. Muller: "Uncertainties in reactive uptake
!   coefficients for solid stratospheric particles. 2. Effect on ozone
!   depletion"; Geophysical Research Letters, Vol. 24, No. 14, pp. 1747-1750,
!   1997
! W. B. DeMore, S. P. Sander, D. M. Golden, R. F. Hampson, M. J. Kurylo, C. J.
!   Howard, A. R. Ravishankara, C. E. Kolb, M. J. Molina: "Chemical Kinetics
!   and Photochemical Data for Use in Stratospheric Modeling, Evaluation Number
!   12"; JPL Publication 97-4, Jet Propulsion Laboratory, California Institute
!   of Technology, Pasadena, CA, 1997; http://jpldataeval.jpl.nasa.gov
! K. Drdla, R. P. Turco, S. Elliott: "Heterogeneous chemistry on antarctic
!   polar stratospheric clouds - a microphysical estimate of the extent of
!   chemical processing"; Journal of Geophysical Research - Atmospheres,
!   Vol. 98, No. D5, pp. 8965-8981, 1993
! R. G. Grainger, A. Lambert, C. D. Rodgers, F. W. Taylor, T. Deshler:
!   "Stratospheric aerosol effective radius, surface area and volume estimated
!   from infrared measurements", Journal of Geophysical Research, Vol. 100,
!   No. D8, pp. 16507--16518, August 1995
! D. Hanson, K. Mauersberger: "Laboratory studies of the nitric-acid
!   trihydrate - implications for the south polar stratosphere"; Geophysical
!   Research Letters, Vol. 15, No. 8, pp. 855-858, 1988
! D. R. Hanson, A. R. Ravishankara, S. Solomon: "Heterogeneous reactions in 
!   sulfuric acid aerosols: A framework for model calculations", Journal of 
!   Geophysical Research, Vol. 99, No. D2, pp. 3615--3629, February 1994
! T. Huthwelker, T. Peter, B. P. Luo, S. L. Clegg, K. S. Carslaw, P.
!   Brimblecombe: "Solubility of HOCl in Water and aqueous H2SO4 to
!   stratospheric temperatures"; Journal of Atmospheric Chemistry, Vol. 21,
!   No. 1, pp. 81-95, 1995
! P. Joeckel, R. Sander, A. Kerkweg, H. Tost, J. Lelieveld: "Technical Note: 
!   The Modular Earth Submodel System (MESSy) - a new approach towards Earth
!   System Modeling; Atmospheric Chemistry and Physics, Vol. 5, pp. 433-444, 
!   2005, www.atmos-chem-phys.org/acp/5/433/
! B. Luo, K. S. Carslaw, T. Peter, S. L. Clegg: "Vapour pressures of
!   H2SO4/HNO3/HCl/HBr/H2O solutions to low stratospheric temperatures";
!   Geophysical Research Letters, Vol. 22, No.3, pp. 247-250, 1995
! B. P. Luo, U. K. Krieger, T. Peter: "Densities and refractive indices of
!   H2SO4/HNO3/H2O solutions to stratospheric temperatures"; Geophysical
!   Research Letters, Vol. 23, No. 25, pp. 3707-3710, 1996
! J. Marti, K. Mauersberger: "A survey and new measurements of ice vapor
!   pressure at temperatures between 170 and 250 K"; Geophysical Research
!   Letters, Vol. 20, No. 5, pp. 363-366, 1993
! S. P. Sander, R. R. Friedl, W. B. DeMore, D. M. Golden, M. J. Kurylo, R. F.
!   Hampson, R. E. Huie, G. K. Moortgat, A. R. Ravishankara, C. E. Kolb, M. J.
!   Molina: \glqq Chemical Kinetics and Photochemial Data for Use in
!   Stratospheric Modeling, Evaluation Number 13\grqq ; JPL Publication 00-3,
!   Jet Propulsion Laboratory, California Institute of Technology, Pasadena,
!   CA, 2000; http://jpldataeval.jpl.nasa.gov
! A. Tabazadeh, R. P. Turco, K. Drdla, M. Z. Jacobson, O. B. Toon: "A study
!   of type-I polar stratospheric cloud formation"; Geophysical Research
!   Letters, Vol. 21, No. 15, pp. 1619-1622, 1994
! A. Waibel: "Anomalien ozonchemisch relevanter Spurengase" (Dissertation
!   at the University Heidelberg); Shaker, Aachen, 1997
! C. J. Walcek: "Minor flux adjustment near mixing ratio extremes for
!   simplified yet highly accurate monotonic calculatio of tracer advection"
!   Journal of Geophysical Research, Vol. 105, No. D7, pp. 9335-9348, 2000
! J. Buchholz, Simulations of physics and chemistry of polar stratospheric
!   clouds with a general circulation model, Dissertation, Johannes Gutenberg 
!   University Mainz, 2005.
!   (http://nbn-resolving.de/urn/resolver.pl?urn=urn:nbn:de:hebis:77-8187)
! O. Kirner, Prozessstudien der stratosphaerischen Chemie und Dynamik mit
!   Hilfe des Chemie-Klima-Modells ECHAM5/MESSy1 Dissertation, Universitaet
!   Karsruhe, 2008.
!   (http://digbib.ubka.uni-karlsruhe.de/volltexte/1000010199)
! M. M. P. van den Broek, J. E. Williams, and A. Bregman,
!   Implementing growth and sedimentation of NAT particles in a
!   global Eulerian model, Atmos. Chem. Phys., 4, 1869-1883, 2004.
!   (http://www.atmos-chem-phys.org/acp/4/1869/)

!=======================================================================

MODULE messy_msbm

  USE messy_main_constants_mem, ONLY: &
    dp,                               & ! kind parameter for real
    i4,                               & ! kind parameter for integer
    atm2Pa,                           & ! pressure unit conversion factor
    pi,                               & !=pi
    N_Avogadro => N_A,                & !=Avogadro constant / (1/mol) 
    MolMassH2O => M_H2O,              & !=molar mass of H2O / (g/mol)
    g_acc => g,                       & !=gravity acceleration / (m/s**2)
    Rgas => R_gas                       !=universal gas constant / (J/(mol*K))

  IMPLICIT NONE

!-----------------------------------------------------------------
! Everything is PRIVATE, except when explicitely stated otherwise:
!-----------------------------------------------------------------
  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modver = '2.4'
  CHARACTER(LEN=*), PARAMETER :: modstr = 'msbm'
  
  !gas constant / (atm*m**3/(K*mol)): 
  REAL(dp), PARAMETER :: R1 = Rgas/atm2Pa
  ! mass density of NAT / (kg/m**3): 
  ! valid for 190K (Drdla et al., 1993):
  REAL(dp), PARAMETER :: dens_nat= 1620.0_dp
  ! mass density of ice / (kg/m**3): 
  REAL(dp), PARAMETER :: dens_ice= 990.0_dp
  ! molar masses / (g/mol):
  REAL(dp), PARAMETER ::            &
    MolMassH2SO4= 98.076_dp,        &
    MolMassHNO3 = 63.012_dp,        &
    MolMassHCl  = 36.461_dp,        &
    MolMassHOCl = 52.46_dp,         &
    MolMassHBr  = 80.91_dp,         &
    MolMassHOBr = 96.91_dp

  ! mz_rs_20060123+
  ! NOTE: THIS LIST (IHS_*) MUST BE THE SAME AS IN MECCA_KHET (gas.eqn)
  INTEGER, PUBLIC, PARAMETER :: IHS_MAX = 13
  ! ihs_ = index of stratospheric heterogeneous reactions
  INTEGER, PUBLIC, PARAMETER :: &
    ihs_N2O5_H2O  =  1, ihs_HOCl_HCl  =  2, ihs_ClNO3_HCl =  3, &
    ihs_ClNO3_H2O =  4, ihs_N2O5_HCl  =  5, ihs_ClNO3_HBr =  6, &
    ihs_BrNO3_HCl =  7, ihs_HOCl_HBr  =  8, ihs_HOBr_HCl  =  9, &
    ihs_HOBr_HBr  = 10, ihs_BrNO3_H2O = 11, ihs_Hg        = 12, &
    ihs_RGM       = 13
  CHARACTER(LEN=9), PUBLIC, DIMENSION(IHS_MAX) :: khet_St_name = (/ &
    'N2O5_H2O ', 'HOCl_HCl ', 'ClNO3_HCl', &
    'ClNO3_H2O', 'N2O5_HCl ', 'ClNO3_HBr', &
    'BrNO3_HCl', 'HOCl_HBr ', 'HOBr_HCl ', &
    'HOBr_HBr ', 'BrNO3_H2O', 'Hg       ', 'RGM      ' /)
  ! mz_rs_20060123-

! ka_ok_20100115+
  INTEGER,  PARAMETER   :: NSB = 9  ! number of size bins
  REAL(dp), PARAMETER   :: MolMassNAT  = 3.0_dp*MolMassH2O+MolMassHNO3
  REAL(dp), PARAMETER   :: rad_ini     = 1.e-7_dp
  REAL(dp), PARAMETER   :: ctoa        = 7.336E21_dp
  ! crystal mass density of NAT in g/m**3
  REAL(dp), PARAMETER   :: rho_NAT     = 1.626e6_dp
  LOGICAL, SAVE         :: KinPar
! ka_ok_20100115-

  ! advection influence on ice/NAT nucleation: 
  LOGICAL, SAVE :: LAdvectIceNat
  ! switch for homogeneous nucleation of NAT: 
  LOGICAL, SAVE :: LHomNucNAT
  ! supercooling for NAT formation: 
  REAL(dp), SAVE :: NatFormThreshold
  ! minimum reaction rate / (cm**3/s): 
  REAL(dp), SAVE :: minKhet
  ! maximum reaction rate / (cm**3/s): 
  REAL(dp), SAVE :: maxKhet
  ! supersaturation required for ice formation: 
  REAL(dp), SAVE :: SupSatIce
  ! minimum radius of solid aerosol particles / m:
  REAL(dp), SAVE :: r_min
  ! maximum solid particle number concentration / (1/m**3):
  REAL(dp), SAVE :: N_solid_max
  ! switch for sedimentation scheme: 
  INTEGER(i4), SAVE :: SedScheme


!---------------
!     PUBLIC
!---------------
  PUBLIC :: dp, modver, modstr
  PUBLIC :: LAdvectIceNat, LHomNucNAT, NatFormThreshold, &
            minKhet, maxKhet, &
            SupSatIce, r_min, N_solid_max, &
            SedScheme
  PUBLIC :: MolMassH2O, MolMassH2SO4, &
            MolMassHCl, MolMassHOCl, &
            MolMassHBr, MolMassHOBr, &
            Rgas
  ! ka_ok_20100115+
  PUBLIC :: KinPar, MolMassNAT, rad_ini, ctoa, rho_NAT, NSB
  PUBLIC :: K_PARA_NAT         ! calculates growth of NAT particles
  ! ka_ok_20100115-
  PUBLIC :: msbm_read_nml_ctrl       ! read CTRL namelist
  PUBLIC :: mz_psc_phase       ! checks the phase of strat. particles (PRESS)
  PUBLIC :: mz_psc_liq_check   ! checks wether cond. meet param. for liquids
  PUBLIC :: mz_psc_ice_H2O     ! calculates ice phase H2O
  PUBLIC :: mz_psc_nat_HNO3    ! calculates nat phase HNO3
  PUBLIC :: mz_psc_liq_bH2SO4b ! calculates auxiliary variable
  PUBLIC :: mz_psc_liq_bHNO3b  ! calculates auxiliary variable
  PUBLIC :: mz_psc_liq_bH2SO4  ! calculates liquid particle composition
  PUBLIC :: mz_psc_liq_bHNO3   ! calculates liquid particle composition
  PUBLIC :: mz_psc_liq_partHNO3! calculates liquid particle composition
  PUBLIC :: mz_psc_liq_wnen    ! calculates auxiliary variable
  PUBLIC :: mz_psc_liq_wHNO3   ! calculates liquid particle composition
  PUBLIC :: mz_psc_liq_wH2SO4  ! calculates liquid particle composition
  PUBLIC :: mz_psc_liq_hHCl    ! calculates liquid particle composition
  PUBLIC :: mz_psc_liq_hHBr    ! calculates liquid particle composition
  PUBLIC :: mz_psc_liq_hHOCl   ! calculates liquid particle composition
  PUBLIC :: mz_psc_liq_hHOBr   ! calculates liquid particle composition
  PUBLIC :: mz_psc_liq_partHBr ! calculates liquid particle composition
  PUBLIC :: mz_psc_liq_wHCl    ! calculates liquid particle composition
  PUBLIC :: mz_psc_liq_wHOCl   ! calculates liquid particle composition
  PUBLIC :: mz_psc_liq_wHOBr   ! calculates liquid particle composition
  PUBLIC :: mz_psc_liq_H2O     ! calculates liquid particle composition 
  PUBLIC :: mz_psc_liq_HCl     ! calculates liquid particle composition 
  PUBLIC :: mz_psc_liq_HOCl    ! calculates liquid particle composition 
  PUBLIC :: mz_psc_liq_HOBr    ! calculates liquid particle composition 
  PUBLIC :: mz_psc_N_solid     ! calculates solid particle number density
  PUBLIC :: mz_psc_r_solid     ! calculates solid particle radius
! ka_ok_20100115+
  PUBLIC :: mz_psc_N_ice       ! calculates ice particle number density
  PUBLIC :: mz_psc_r_ice       ! calculates ice particle radius
! ka_ok_20100115-
  PUBLIC :: mz_psc_dim_liquid  ! calc. log-normal distribution of liquid ptcls
  PUBLIC :: mz_psc_surface_liquid ! liq. particle surface area density
  PUBLIC :: mz_psc_het_nat1    ! calc. het. reaction rate on nat particles
  PUBLIC :: mz_psc_het_nat2    ! calc. het. reaction rate on nat particles
  PUBLIC :: mz_psc_het_nat3    ! calc. het. reaction rate on nat particles
  PUBLIC :: mz_psc_het_nat4    ! calc. het. reaction rate on nat particles
  PUBLIC :: mz_psc_het_nat5    ! calc. het. reaction rate on nat particles
  PUBLIC :: mz_psc_het_nat6    ! calc. het. reaction rate on nat particles
  PUBLIC :: mz_psc_het_nat7    ! calc. het. reaction rate on nat particles
  PUBLIC :: mz_psc_het_nat8    ! calc. het. reaction rate on nat particles
  PUBLIC :: mz_psc_het_nat9    ! calc. het. reaction rate on nat particles
  PUBLIC :: mz_psc_het_nat10   ! calc. het. reaction rate on nat particles
  PUBLIC :: mz_psc_het_nat11   ! calc. het. reaction rate on nat particles
  PUBLIC :: mz_psc_het_ice12   ! calc. het. reaction rate on ice particles
  PUBLIC :: mz_psc_het_ice13   ! calc. het. reaction rate on ice particles
  PUBLIC :: mz_psc_het_ice14   ! calc. het. reaction rate on ice particles
  PUBLIC :: mz_psc_het_ice15   ! calc. het. reaction rate on ice particles
  PUBLIC :: mz_psc_het_ice16   ! calc. het. reaction rate on ice particles
  PUBLIC :: mz_psc_het_ice17   ! calc. het. reaction rate on ice particles
  PUBLIC :: mz_psc_het_ice18   ! calc. het. reaction rate on ice particles
  PUBLIC :: mz_psc_het_ice19   ! calc. het. reaction rate on ice particles
  PUBLIC :: mz_psc_het_ice20   ! calc. het. reaction rate on ice particles
  PUBLIC :: mz_psc_het_ice21   ! calc. het. reaction rate on ice particles
  PUBLIC :: mz_psc_het_ice22   ! calc. het. reaction rate on ice particles
  PUBLIC :: mz_psc_het_liq23   ! calc. het. reaction rate on liquid particles
  PUBLIC :: mz_psc_het_liq24   ! calc. het. reaction rate on liquid particles
  PUBLIC :: mz_psc_het_liq25   ! calc. het. reaction rate on liquid particles
  PUBLIC :: mz_psc_het_liq26   ! calc. het. reaction rate on liquid particles
  PUBLIC :: mz_psc_het_liq27   ! calc. het. reaction rate on liquid particles
  PUBLIC :: mz_psc_het_liq28   ! calc. het. reaction rate on liquid particles
  PUBLIC :: mz_psc_het_liq29   ! calc. het. reaction rate on liquid particles
  PUBLIC :: mz_psc_het_liq30   ! calc. het. reaction rate on liquid particles
  PUBLIC :: mz_psc_het_liq_gen ! calc. het. reaction rate on liquid particles
  PUBLIC :: mz_psc_het_sol_gen ! calc. het. reaction rate on ice/nat particles
  PUBLIC :: mix2conc           ! converts mixing ratio to concentration
  PUBLIC :: mz_psc_density
  PUBLIC :: mz_psc_diff
  PUBLIC :: mz_psc_vel         ! sedimentation velocity
  ! ka_ok_20100115+
  !                            ! ... (of ice if KinPar=.true.)
  ! ka_ok_20100115-
  PUBLIC :: mz_psc_SedStep     ! calculates sedimentation steps
  PUBLIC :: mz_psc_sed         ! particle redistribution due to sedimentation
  PUBLIC :: T_nat, T_ice
  ! ka_ok_20100115+
  ! elemental function for kinetic PSC parametrisation
  PUBLIC :: diff_H2O           ! diffusion coefficient of H2O in air
  PUBLIC :: diff_HNO3          ! diffusion coefficient of HNO3 in air
  PUBLIC :: speed_HNO3         ! mean molecular speed of HNO3 
  PUBLIC :: press_HNO3         ! HNO3 partial pressure
  PUBLIC :: press_H2O          ! H2O partial pressure
  PUBLIC :: press_HNO3_over_NAT ! HNO3 partial pressure over NAT
  PUBLIC :: press_H2OS         ! saturation vapor pressure of H2O
  PUBLIC :: press_HNO3S        ! saturation vapor pressure of HNO3
  PUBLIC :: svapn_fkt    
  PUBLIC :: growth_factor      ! growth factor for NAT particles
  PUBLIC :: mz_NAT_vel         ! sedimentation velocity of NAT particles 
  !                              for every sizebin
  PUBLIC :: mean_free_path     ! mean free path of a molecule
  PUBLIC :: viscosity_air      ! viscosity of air
  ! ka_ok_20100115-

!------------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------

!=========================================================================
SUBROUTINE msbm_read_nml_ctrl(status, iou)
  !--------------------------------------------------------------------
  ! This routine reads and checks the msbm CTRL namelist. It is designed
  ! according to the MESSy (Mainz Earth Submodel System) standard. 
  !--------------------------------------------------------------------

  USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
  
  IMPLICIT NONE
  
  !-----
  ! input/output variables
  !-----
  INTEGER, INTENT(out) :: status   ! error status
  INTEGER, INTENT(in) :: iou       ! logical I/O unit
  
  
  !-----
  ! local variables
  !-----
  CHARACTER(len=*), PARAMETER :: substr = 'msbm_read_nml_ctrl'
  LOGICAL                     :: lex     ! file existence flag
  INTEGER                     :: fstat   ! file status

  NAMELIST /CTRL/ LAdvectIceNat, LHomNucNAT, NatFormThreshold, &
                  minKhet, maxKhet, &
                  SupSatIce, r_min, N_solid_max, SedScheme,    &
                  KinPar ! ka_ok_20100115

  !-----
  ! initialisation
  !-----
  status = 1
  LAdvectIceNat = .false. 
  LHomNucNAT = .false.
  NatFormThreshold = -3.0_dp
  minKhet = 0.0_dp
  maxKhet = 1.0e-13_dp
  SupSatIce = 1.5_dp
  r_min = 1.0e-7_dp
  N_solid_max = 0.01e6_dp
  SedScheme = 3
  KinPar=.false. ! ka_ok_20100115
  
  CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
  IF (.not.lex) RETURN   ! msbm.nml does not exist
  
  read(iou, nml=CTRL, iostat=fstat)
  CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
  IF (fstat /= 0) RETURN   ! error while reading namelist

  ! ka_ok_20100115+
  write (*,*)'switch for kinetic PSC parameterisation=',KinPar
  ! ka_ok_20100115-

  write (*,*) &
    'switch for advection influence on ice/NAT formation: LAdvectIceNat = ', &
    LAdvectIceNat
  write (*,*) 'switch for homogeneous NAT nucleation: LHomNucNAT = ', LHomNucNAT
  write (*,*) 'supercooling for NAT formation: NatFormThreshold = ', &
    NatFormThreshold
  write (*,*) 'minimum reaction rate / (cm**3/s) = ', minKhet
  write (*,*) 'maximum reaction rate / (cm**3/s) = ', maxKhet
  write (*,*) 'supersaturation required for ice formation = ', SupSatIce
  write (*,*) 'minimum radius of solid aerosol particles / m = ', r_min
  write (*,*) 'maximum solid particle number concentration / (1/m**3) = ', &
    N_solid_max
  write (*,*) 'switch for sedimentation scheme: SedScheme = ', SedScheme
  
  CALL read_nml_close(substr, iou, modstr)
  status = 0   ! no error

END SUBROUTINE msbm_read_nml_ctrl
!=========================================================================



!=========================================================================
ELEMENTAL FUNCTION mz_psc_phase(TEMP, PRESS, H2O_tot, HNO3_tot, phase, &
                                H2O_ice, HNO3_nat)
  !-------------------------------------------------------------------------
  ! This function checks, which phases exist:
  ! 1 = binary or ternary liquid aerosol only
  ! 2 = mixture of liquid and solid NAT particles
  ! 3 = mixture of liquid and solid NAT and ice particles
  !-------------------------------------------------------------------------
  ! The function does NOT consider SAT and is only valid if T_nat>T_ice
  ! It assumes, that NAT nucleates on ice and ice forms at supersaturation
  ! S=pH2O_tot/p_ice=SupSatIce. 
  ! The additional possibility of homogeneous NAT nucleation is controlled
  ! via the switch LHomNucNAT. 
  ! The question whether ice or nat are already present can be answered based
  ! on the input variable phase only (switch LAdvectIceNat==.false.) or based
  ! on the presence of H2O_ice and HNO3_nat, too (switch LAdvectIceNat==.true.).
  !-------------------------------------------------------------------------
  ! INPUT:
  ! ------
  ! TEMP    = temperature / K
  ! PRESS   = pressure / hPa
  ! HX_tot  = amount-of-substance ratio of total HX to dry air / (mol/mol)
  !           Note: HX is H2O or HNO3.
  ! phase   = phase indicator (1=liquid, 2=liquid+nat, 3=liquid+nat+ice)
  ! H2O_ice = amount-of-substance ratio of ice phase H2O to dry air / (mol/mol)
  ! HNO3_nat= amount-of-substance ratio of NAT HNO3 to dry air / (mol/mol)
  ! ----------
  ! PARAMETER:
  ! ----------
  ! local_phase = local phase indicator
  ! CRITERIA  = Supersaturation (p_vap/p_part) necessary for phase transitions
  ! HX        = amount-of-substance ratio of gas phase HX to dry air / (mol/mol)
  ! pH2O_tot  = H2O pressure if all H2O molecules were gaseous / atm
  ! pHNO3_tot = HNO3 pressure if all HNO3 molecules were gaseous / atm
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_phase = phase indicator (1=liquid, 2=liquid+nat, 3=liquid+nat+ice)
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: H2O_tot, HNO3_tot
  INTEGER, INTENT(in) :: phase
  REAL(dp), INTENT(in) :: H2O_ice, HNO3_nat

  INTEGER :: mz_psc_phase

  REAL(dp) :: H2O, HNO3
  REAL(dp) :: pH2O_tot, pHNO3_tot
  INTEGER :: local_phase
  REAL(dp), DIMENSION(3,3) :: CRITERIA
  
  INTRINSIC tiny

  ! -----
  ! Supersaturations (CRITERIA=p_vap/p_part) necessary for phase transition
  ! -----
  CRITERIA(1,2)=9999.0_dp  !liquid -> nat    (phase transition not allowed)
  CRITERIA(1,3)=SupSatIce  !liquid -> ice    (supersaturation required)
  CRITERIA(2,3)=1.0_dp     !nat    -> ice    (no supersaturation required)
  CRITERIA(2,1)=1.0_dp     !nat    -> liquid (no subsaturation required)
  CRITERIA(3,1)=9999.0_dp  !ice    -> liquid (only allowed via nat)
  CRITERIA(3,2)=1.0_dp     !ice    -> nat    (no subsaturation required)
  ! -----

  local_phase = phase

  H2O = 1013.25_dp*p_ice(TEMP)/PRESS
  pH2O_tot = H2O_tot*PRESS/1013.25_dp
  pHNO3_tot = HNO3_tot*PRESS/1013.25_dp
  HNO3 = 1013.25_dp*p_nat(TEMP, pH2O_tot)/PRESS

  IF (LAdvectIceNat) THEN
    IF (H2O_ice>tiny(1.0_dp)) THEN
      local_phase=3
    ELSEIF (HNO3_nat>tiny(1.0_dp)) THEN
      local_phase=2
    END IF
  END IF

  IF (local_phase==3) THEN
    IF (H2O_tot<(CRITERIA(3,2)*H2O) .AND. &
        HNO3_tot<(CRITERIA(2,1)*HNO3)) THEN
      ! ice and nat melt, only liquid left
      mz_psc_phase = 1
    ELSEIF (H2O_tot<(CRITERIA(3,2)*H2O) .AND. &
            HNO3_tot>=(CRITERIA(2,1)*HNO3)) THEN
      ! ice melts, nat and liquid left
      mz_psc_phase = 2
    ELSE
      ! ice-nat-liquid-mixture, unchanged
      mz_psc_phase = 3
    END IF
  ELSEIF (local_phase==2) THEN
    IF (HNO3_tot<(CRITERIA(2,1)*HNO3)) THEN
      !liquid only
      mz_psc_phase = 1
    ELSEIF (H2O_tot>(CRITERIA(2,3)*H2O)) THEN
      !ice-nat-liquid-mixture
      mz_psc_phase = 3
    ELSE
      ! nat-liquid-mixture, unchanged
      mz_psc_phase = 2
    END IF
  ELSE
    IF (H2O_tot>=(CRITERIA(1,3)*H2O)) THEN
      ! liquid --> ice-nat-liquid-mixture
      mz_psc_phase = 3
    ELSEIF (LHomNucNAT .AND. &
            (TEMP<=T_nat(pH2O_tot,pHNO3_tot)+NatFormThreshold)) THEN
      ! nat-liquid mixture
      mz_psc_phase = 2
    ELSE
      ! liquid, unchanged
      mz_psc_phase = 1
    END IF
  END IF

!-------------------------------------------------------------------------
END FUNCTION mz_psc_phase
!=========================================================================



!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_check(TEMP, PRESS, H2O, HNO3, H2SO4)
  !-------------------------------------------------------------------------
  ! This function checks, wether the actual conditions meet the conditions
  ! for which the parameterisation for H2O-HNO3-HCl-HBr-HOCl-HOBr liquid
  ! stratospheric aerosol is valid (i_val: 1=valid, 0=not valid)
  !-------------------------------------------------------------------------
  ! INPUT:
  ! ------
  ! TEMP = temperature / K
  ! PRESS = total pressure / hPa
  ! H2O   = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  !         Note: Gas phase depletion due to H2O uptake in droplets is
  !         neglected.
  ! HNO3  = amount-of-substance ratio of gaseous HNO3 to dry air / (mol/mol)
  ! H2SO4 = amount-of-substance ratio of total H2SO4 to dry air / (mol/mol)
  ! ----------
  ! PARAMETER:
  ! ----------
  ! pH2O   = water partial pressure / hPa
  ! Tice   = ice frost point temperature / K
  ! i_val  = 1, if strat.aerosol parameterisation is valid
  !        = 0, if strat.aerosol parameterisation is NOT valid
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_liq_check = i_val
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS, H2O, HNO3, H2SO4
  INTEGER :: i_val, mz_psc_liq_check
  REAL(dp) :: pH2O, Tice


  pH2O = PRESS * H2O
  Tice = T_ice(pH2O)
  i_val = 1

  IF ((TEMP  <= 185.0_dp)  .OR. (TEMP  <= Tice-3.0_dp) .OR.   &
      (TEMP  >= 240.0_dp)  .OR.                               &
      (pH2O  >= 2.0e-3_dp) .OR. (pH2O  <= 2.0e-5_dp)   .OR.   &
      (H2O   >= 1.0e-5_dp) .OR. (H2O   <= 1.0e-6_dp)   .OR.   &
      (HNO3  >= 2.0e-8_dp) .OR. (HNO3  <= 0.0_dp)      .OR.   &
      (H2SO4 >= 1.0e-7_dp) .OR. (H2SO4 <= 1.0e-10_dp)  .OR.   &
      (PRESS >= 200.0_dp)  .OR. (PRESS <= 20.0_dp)) THEN
     i_val = 0
  ENDIF

  mz_psc_liq_check = i_val
!-------------------------------------------------------------------------
END FUNCTION mz_psc_liq_check
!=========================================================================



!==========================================================================
ELEMENTAL FUNCTION mz_psc_ice_H2O(phase, TEMP, PRESS, H2O_tot)
  !--------------------------------------------------------------------------
  ! This function calculates ice phase H2O as a function of temperature,
  ! pressure, and total water.
  !--------------------------------------------------------------------------
  ! INPUT:
  ! ------
  ! phase   = phase indicator (1=liquid, 2=liquid+nat, 3=liquid+nat+ice)
  ! TEMP    = temperature / K
  ! PRESS   = pressure / hPa
  ! H2O_tot = amount-of-substance ratio of total H2O to dry air / (mol/mol)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_ice_H2O
  !         = amount-of-substance ratio of ice phase H2O to dry air / (mol/mol)
  !--------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in)     :: phase
  REAL(dp), INTENT(in)    :: TEMP, PRESS, H2O_tot

  REAL(dp) :: mz_psc_ice_H2O
  
  INTRINSIC max

  IF (phase==3) THEN
    !-----
    ! ice, nat, liquid present
    ! ice formation leads to dehydration
    !-----
    mz_psc_ice_H2O = max(H2O_tot - p_ice(TEMP)*1013.25_dp/PRESS, 0.0_dp)
  ELSE
    !-----
    ! no ice
    !-----
    mz_psc_ice_H2O  = 0.0_dp
  END IF
END FUNCTION mz_psc_ice_H2O
!==========================================================================




!==========================================================================
ELEMENTAL FUNCTION mz_psc_nat_HNO3(phase, TEMP, PRESS, &
                                   H2O_tot, HNO3_tot)
  !--------------------------------------------------------------------------
  ! This function calculates gas phase HNO3 as a function of temperature,
  ! pressure, total water, and total HNO3.
  !--------------------------------------------------------------------------
  ! INPUT:
  ! ------
  ! phase  = phase indicator (1=liquid, 2=liquid+nat, 3=liquid+nat+ice)
  ! TEMP   = temperature / K
  ! PRESS  = pressure / hPa
  ! HX_tot = amount-of-substance ratio of total HX to dry air / (mol/mol)
  !          Note: HX is H2O or HNO3.
  ! ----------
  ! PARAMETER:
  ! ----------
  ! pice     = equilibrium H2O vapour pressure over ice / atm
  !            Note: The H2O equilibrium vapour pressure equals the H2O
  !            partial pressure.
  ! pnat     = equilibrium HNO3 vapour pressure over NAT / atm
  !            Note: The HNO3 equilibrium vapour pressure equals the HNO3
  !            partial pressure.
  ! pH2O_tot = H2O pressure if all H2O molecules were gaseous / atm
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_nat_HNO3
  !        = amount-of-substance ratio of NAT phase HNO3 to dry air / (mol/mol)
  !--------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in)     :: phase
  REAL(dp), INTENT(in)    :: TEMP, PRESS, H2O_tot, HNO3_tot

  REAL(dp) :: mz_psc_nat_HNO3
  REAL(dp) :: pice, pnat, pH2O_tot
  
  INTRINSIC max
  
  IF (phase==3) THEN
    !-----
    ! ice, nat, liquid present
    ! some HNO3 in nat particles
    !-----
    pice = p_ice(TEMP)
    pnat = p_nat(TEMP,pice)
    mz_psc_nat_HNO3 = max(HNO3_tot - pnat*1013.25_dp/PRESS, 0.0_dp)
  ELSEIF (phase==2) THEN
    !-----
    ! nat, liquid present
    ! some HNO3 in nat particles
    !-----
    pH2O_tot = PRESS/1013.25_dp * H2O_tot
    pnat = p_nat(TEMP, pH2O_tot)
    mz_psc_nat_HNO3 = max(HNO3_tot - pnat*1013.25_dp/PRESS, 0.0_dp)
  ELSE
    !-----
    ! only liquid present
    ! all HNO3 gaseous
    !-----
    mz_psc_nat_HNO3 = 0.0_dp
  END IF

END FUNCTION mz_psc_nat_HNO3
!==========================================================================



!============================================================================
ELEMENTAL FUNCTION mz_psc_liq_bH2SO4b(TEMP, PRESS, H2O_gl)
  !-------------------------------------------------------------------------
  ! This elemental function is one of the functions that calculate the
  ! composition of liquid stratospheric aerosol droplets according to
  ! Carslaw et al., GRL, 1995
  !-------------------------------------------------------------------------
  ! The HNO3/H2SO4 composition is based on a thermodynamic model of the
  ! HCl-HNO3-H2SO4-H2O-System (Carslaw et al, JPC, 1995).
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP     = temperature / K
  ! PRESS    = total pressure / hPa
  ! H2O_gl   = amount-of-substance ratio of gaseous+liquid H2O to dry air
  !            / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! rTEMP    = temperature (restricted) / K
  ! pH2O_tot = H2O partial pressure if all H2O was gaseous / atm
  ! xsb = ??? fraction of binary H2SO4/H2O solution
  ! ks = auxiliary variable
  ! -------
  ! OUTPUT:
  ! -------
  ! bH2SO4b  = molality of H2SO4 in water for binary solution / (mol/kg)
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: H2O_gl

  REAL(dp) bH2SO4b, mz_psc_liq_bH2SO4b

  REAL(dp) :: pH2O_tot, rTEMP
  REAL(dp) :: xsb
  REAL(dp), PARAMETER :: &
       ks1 = -21.661_dp, &
       ks2 = 2724.2_dp, &
       ks3 = 51.81_dp, &
       ks4 = -15732.0_dp, &
       ks5 = 47.004_dp, &
       ks6 = -6969.0_dp, &
       ks7 = -4.6183_dp

  INTRINSIC log, max, min, sqrt
  
  rTEMP = min(max(185.0_dp,TEMP),240.0_dp)

  pH2O_tot  = min(max(PRESS*H2O_gl,2e-5_dp),2e-3_dp)/1013.25_dp

  xsb = 1.0_dp / (2.0_dp * (ks3 + ks4/rTEMP))                 &
       * ( - ks1 - ks2/rTEMP                                  &
       - SQRT( (ks1 + ks2/rTEMP)**2 - 4.0_dp*(ks3           &
       +ks4/rTEMP)                                              &
       * (ks5 + ks6/rTEMP + ks7*LOG(rTEMP)-LOG(pH2O_tot)) ) )
  bH2SO4b = 55.51_dp * xsb / (1.0_dp - xsb)

  mz_psc_liq_bH2SO4b = bH2SO4b
END FUNCTION mz_psc_liq_bH2SO4b
!============================================================================



!============================================================================
ELEMENTAL FUNCTION mz_psc_liq_bHNO3b(TEMP, PRESS, H2O_gl)
  !-------------------------------------------------------------------------
  ! This elemental function is one of the functions that calculate the
  ! composition of liquid stratospheric aerosol droplets according to
  ! Carslaw et al., GRL, 1995
  !-------------------------------------------------------------------------
  ! The HNO3/H2SO4 composition is based on a thermodynamic model of the
  ! HCl-HNO3-H2SO4-H2O-System (Carslaw et al, JPC, 1995).
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP = temperature / K
  ! PRESS    = total pressure / hPa
  ! H2O_gl   = amount-of-substance ratio of gaseous+liquid H2O to dry air
  !            / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! rTEMP    = temperature (restricted) / K
  ! pH2O_tot = H2O partial pressure if all H2O was gaseous / atm
  ! xnb = ??? fraction of binary HNO3/H2O solution
  ! kn
  ! -------
  ! OUTPUT:
  ! -------
  ! bHNO3b   = molality of HNO3 in water for binary solution / (mol/kg)
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: H2O_gl

  REAL(dp) :: bHNO3b, mz_psc_liq_bHNO3b

  REAL(dp) :: pH2O_tot, rTEMP
  REAL(dp) :: xnb
  REAL(dp), PARAMETER :: &
       kn1 = -39.136_dp, &
       kn2 = 6358.4_dp, &
       kn3 = 83.29_dp, &
       kn4 = -17650.0_dp, &
       kn5 = 198.53_dp, &
       kn6 = -11948.0_dp, &
       kn7 = -28.469_dp

  INTRINSIC log, max, min, sqrt, tiny

  IF (TEMP>215.0_dp) THEN
    mz_psc_liq_bHNO3b = 0.0_dp
    RETURN
  END IF
  
  rTEMP = max(185.0_dp,TEMP)

  pH2O_tot  = min(max(PRESS*H2O_gl,2e-5_dp),2e-3_dp)/1013.25_dp

  IF ((kn3+kn4/rTEMP)<TINY(xnb)) THEN
    xnb = -(kn5 + kn6/rTEMP + kn7*LOG(rTEMP)-LOG(pH2O_tot)) &
          /(kn1 + kn2/rTEMP)
  ELSE
    xnb = 1.0_dp / (2.0_dp * (kn3 + kn4/rTEMP))               &
         * ( - kn1 - kn2/rTEMP                                &
         - SQRT( (kn1 + kn2/rTEMP)**2                         &
         -4.0_dp*(kn3+kn4/rTEMP)                              &
         * (kn5 + kn6/rTEMP + kn7*LOG(rTEMP)-LOG(pH2O_tot))))
  END IF
  bHNO3b = 55.51_dp * xnb / (1.0_dp - xnb)

  mz_psc_liq_bHNO3b = bHNO3b
END FUNCTION mz_psc_liq_bHNO3b
!============================================================================




!============================================================================
 ELEMENTAL FUNCTION mz_psc_liq_bH2SO4(TEMP, PRESS,               &
                                       H2SO4_gl, H2O_gl, HNO3_gl, &
                                       bH2SO4b, bHNO3b)
  !-------------------------------------------------------------------------
  ! This elemental function is one of the functions that calculate the
  ! composition of liquid stratospheric aerosol droplets according to
  ! Carslaw et al., GRL, 1995
  !-------------------------------------------------------------------------
  ! The HNO3/H2SO4 composition is based on a thermodynamic model of the
  ! HCl-HNO3-H2SO4-H2O-System (Carslaw et al, JPC, 1995).
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP  = temperature / K
  ! PRESS = total pressure / hPa
  ! HX_gl = amount-of-substance ratio of gaseous+liquid HX to dry air 
  !         / (mol/mol)
  !         Note: HX_gl is H2SO4_gl, H2O_gl, HNO3_gl
  ! -------
  ! OUTPUT:
  ! -------
  ! bH2SO4   = molality of H2SO4 in water / (mol/kg)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! rTEMP  = temperature (restricted) / K
  ! bHNO3    = molality of HNO3 in water  / (mol/kg)
  ! ns = auxiliary variable
  ! pH2O_tot = H2O partial pressure if all H2O was gaseous / atm
  ! pHNO3_tot = HNO3 partial pressure if all HNO3 was gaseous / atm
  ! pHNO3 = HNO3 equilibrium vapour pressure over liquid solution / atm
  !         Note: pHNO3 is calculated from the partitioning of HNO3 between
  !         the liquid aerosol and the gas phase
  ! a, b, c = auxiliary variables
  ! pr, pr2, pr3, phi = auxiliary variables
  ! sqrtfunc = auxiliary variable
  ! tt, tr, tr2, tr3 = auxiliary variables
  ! hnb = solubility of binary HNO3/H2O solution
  ! hsb = solubility of binary H2SO4/H2O solution
  ! qn, qs = auxiliary variables
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: H2SO4_gl, H2O_gl, HNO3_gl
  REAL(dp), INTENT(in) :: bH2SO4b, bHNO3b

  REAL(dp) :: bH2SO4, mz_psc_liq_bH2SO4

  REAL(dp) :: bHNO3, rTEMP
  REAL(dp) :: ns
  REAL(dp) :: pH2O_tot, pHNO3_tot
  REAL(dp) :: pHNO3
  REAL(dp) :: a,b,c
  REAL(dp) :: pr, pr2, pr3, phi
  REAL(dp) :: sqrtfunc
  REAL(dp) :: tt, tr, tr2, tr3
  REAL(dp) :: hnb, hsb
  REAL(dp), PARAMETER :: &
       qn1 = 14.5734_dp, &
       qn2 = 0.0615994_dp, &
       qn3 = -1.14895_dp, &
       qn4 = 0.691693_dp, &
       qn5 = -0.098863_dp, &
       qn6 = 0.0051579_dp, &
       qn7 = 0.123472_dp, &
       qn8 = -0.115574_dp, &
       qn9 = 0.0110113_dp, &
       qn10 = 0.0097914_dp, &
       qs1 = 14.4700_dp, &
       qs2 = 0.0638795_dp, &
       qs3 = -3.29597_dp, &
       qs4 = 1.778224_dp, &
       qs5 = -0.223244_dp, &
       qs6 = 0.0086486_dp, &
       qs7 = 0.536695_dp, &
       qs8 = -0.335164_dp, &
       qs9 = 0.0265153_dp, &
       qs10 = 0.0157550_dp
  
  INTRINSIC atan, cos, epsilon, exp, log, max, min, sqrt

  IF (TEMP>215.0_dp .OR. HNO3_gl<epsilon(1.0_dp)) THEN
    !-----
    ! assume solution is pure H2SO4/H2O
    !-----
    mz_psc_liq_bH2SO4 = bH2SO4b
    RETURN
  END IF

  rTEMP = max(TEMP,185.0_dp)

  !-----
  ! convert input
  !-----
  ns= min(max(H2SO4_gl,1.0e-10_dp),1.0e7_dp)*PRESS/1013.25_dp/R1/rTEMP

  pH2O_tot  = min(max(PRESS*H2O_gl,2e-5_dp),2e-3_dp)/1013.25_dp
  pHNO3_tot = PRESS * min(HNO3_gl,2.0e-8_dp) / 1013.25_dp


  tr   = 1.0e4_dp / rTEMP - 43.4782608_dp
  pr   = LOG(pH2O_tot) + 18.4_dp
  tr2  = tr**2
  tr3  = tr**3
  pr2  = pr**2
  pr3  = pr**3

  !-----
  ! the HNO3/H2SO4/H2O solution composition
  !-----
  hsb = qs1                                            &
       + qs2*tr2                                       &
       + (qs3 + qs4*tr + qs5*tr2 + qs6*tr3)*pr   &
       + (qs7 + qs8*tr + qs9*tr2)*pr2              &
       + qs10*tr*pr3
  hsb = EXP(hsb)
  hnb = qn1                                          &
       + qn2*tr2                                     &
       + (qn3 + qn4*tr + qn5*tr2 + qn6*tr3)*pr &
       + (qn7 + qn8*tr + qn9*tr2)*pr2            &
       + qn10*tr*pr3
  hnb = EXP(hnb)
  tt   = r1 * rTEMP * ns
  !-----
  ! The following code is taken from Carslaw et al., GRL, 1995
  ! Variables have been renamed, however, as follows:
  ! T    |--> rTEMP
  ! NS   |--> ns
  ! MNB  |--> bHNO3b
  ! MSB  |--> bH2SO4b
  ! PN0  |--> pHNO3_tot
  ! PHI  |--> phi
  ! PI   |--> pi
  ! MS   |--> bH2SO4
  ! MN   |--> bHNO3
  ! PN   |--> pHNO3
  ! WS   |--> wH2SO4
  ! Other changes include the use of sqrtfunc as auxiliary variable.
  !-----
  a = (tt*hnb*bHNO3b**2 - tt*hsb*bHNO3b*bH2SO4b &
       - 2.0_dp*bHNO3b**2*bH2SO4b + bHNO3b*bH2SO4b**2 &
       + hnb*bHNO3b*bH2SO4b*pHNO3_tot - hsb*bH2SO4b**2*pHNO3_tot) &
       / (bHNO3b**2 - bHNO3b * bH2SO4b)
  b = bH2SO4b * (-2.0_dp*tt*hnb*bHNO3b         &
                 +tt*hsb*bH2SO4b           &
                 +bHNO3b*bH2SO4b           &
                 -hnb*bH2SO4b*pHNO3_tot)   &
              / (bHNO3b - bH2SO4b)
  c = (tt*hnb*bHNO3b*bH2SO4b**2)     &
       / (bHNO3b - bH2SO4b)
  sqrtfunc = -2.0_dp*a**3 +9.0_dp*a*b-27.0_dp*c
  phi = ATAN(SQRT(4.0_dp * (a**2 - 3.0_dp* b)**3         &
       - sqrtfunc**2)                                    &
       / sqrtfunc )
  IF (phi < 0.0_dp) phi = phi + pi
  bH2SO4 = - (a + 2.0_dp * SQRT(a**2 - 3.0_dp * b)       &
          * COS((pi + phi) / 3.0_dp)) / 3.0_dp
  IF (bH2SO4 < 0.0_dp) THEN
    bH2SO4 = bH2SO4b / 1.0e5_dp
  ENDIF
  bHNO3 = bHNO3b * (1.0_dp - bH2SO4 / bH2SO4b)
  pHNO3 = bHNO3 / ((hnb*bHNO3 + hsb*bH2SO4)              &
         / (bHNO3+bH2SO4))
  !-----
  ! End of quote from Carslaw et al., GRL, 1995
  !-----
  IF (bH2SO4 > bH2SO4b .OR. pHNO3 > pHNO3_tot) THEN
    bH2SO4 = bH2SO4b
  ENDIF

  mz_psc_liq_bH2SO4 = bH2SO4
END FUNCTION mz_psc_liq_bH2SO4
!=========================================================================




!============================================================================
ELEMENTAL FUNCTION mz_psc_liq_bHNO3(TEMP, PRESS,           &
                                     H2O_gl, HNO3_gl,       &
                                     bHNO3b, bH2SO4b, bH2SO4)
  !-------------------------------------------------------------------------
  ! This elemental function is one of the functions that calculate the
  ! composition of liquid stratospheric aerosol droplets according to
  ! Carslaw et al., GRL, 1995
  !-------------------------------------------------------------------------
  ! The HNO3/H2SO4 composition is based on a thermodynamic model of the
  ! HCl-HNO3-H2SO4-H2O-System (Carslaw et al, JPC, 1995).
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP = temperature / K
  ! PRESS    = total pressure / hPa
  ! H2O_gl   = amount-of-substance ratio of gaseous+liquid H2O to dry air
  !            / (mol/mol)
  ! HNO3_gl   = amount-of-substance ratio of gaseous+liquid HNO3 to dry air
  !            / (mol/mol)
  ! bHNO3b   = molality of HNO3 in water for binary solution / (mol/kg)
  ! bH2SO4b  = molality of H2SO4 in water for binary solution / (mol/kg)
  ! bH2SO4   = molality of H2SO4 in water / (mol/kg)
  ! -------
  ! OUTPUT:
  ! -------
  ! bHNO3b   = molality of HNO3 in water for binary solution / (mol/kg)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! pH2O_tot = H2O partial pressure if all H2O was gaseous / atm
  ! pHNO3_tot = HNO3 partial pressure if all HNO3 was gaseous / atm
  ! pHNO3 = HNO3 equilibrium vapour pressure over liquid solution / atm
  !         Note: pHNO3 is calculated from the partitioning of HNO3 between
  !         the liquid aerosol and the gas phase
  ! pr, pr2, pr3 = auxiliary variables
  ! tr, tr2, tr3 = auxiliary variables
  ! hnb = solubility of binary HNO3/H2O solution
  ! hsb = solubility of binary H2SO4/H2O solution
  ! qn, qs = auxiliary variables
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: H2O_gl, HNO3_gl
  REAL(dp), INTENT(in) :: bHNO3b, bH2SO4b, bH2SO4

  REAL(dp) :: bHNO3, mz_psc_liq_bHNO3

  REAL(dp) :: pH2O_tot, pHNO3_tot
  REAL(dp) :: pHNO3
  REAL(dp) :: pr, pr2, pr3
  REAL(dp) :: tr, tr2, tr3
  REAL(dp) :: hsb, hnb
  REAL(dp), PARAMETER :: &
       qn1 = 14.5734_dp, &
       qn2 = 0.0615994_dp, &
       qn3 = -1.14895_dp, &
       qn4 = 0.691693_dp, &
       qn5 = -0.098863_dp, &
       qn6 = 0.0051579_dp, &
       qn7 = 0.123472_dp, &
       qn8 = -0.115574_dp, &
       qn9 = 0.0110113_dp, &
       qn10 = 0.0097914_dp, &
       qs1 = 14.4700_dp, &
       qs2 = 0.0638795_dp, &
       qs3 = -3.29597_dp, &
       qs4 = 1.778224_dp, &
       qs5 = -0.223244_dp, &
       qs6 = 0.0086486_dp, &
       qs7 = 0.536695_dp, &
       qs8 = -0.335164_dp, &
       qs9 = 0.0265153_dp, &
       qs10 = 0.0157550_dp

  INTRINSIC epsilon, exp, log, max, min

  IF (TEMP>215.0_dp .OR. HNO3_gl<epsilon(1.0_dp)) THEN
    !-----
    ! assume solution is pure H2SO4/H2O
    !-----
    mz_psc_liq_bHNO3 = 0.0_dp
    RETURN
  END IF

  !-----
  ! convert input
  !-----
  pH2O_tot  = min(max(PRESS*H2O_gl,2e-5_dp),2e-3_dp)/1013.25_dp
  pHNO3_tot = PRESS * min(HNO3_gl,2.0e-8_dp) / 1013.25_dp

  tr   = 1.0e4_dp / max(TEMP,185.0_dp) - 43.4782608_dp
  pr   = LOG(pH2O_tot) + 18.4_dp
  tr2  = tr**2
  tr3  = tr**3
  pr2  = pr**2
  pr3  = pr**3

  !-----
  ! the HNO3/H2SO4/H2O solution composition
  !-----
  hsb = qs1                                            &
       + qs2*tr2                                       &
       + (qs3 + qs4*tr + qs5*tr2 + qs6*tr3)*pr   &
       + (qs7 + qs8*tr + qs9*tr2)*pr2              &
       + qs10*tr*pr3
  hsb = EXP(hsb)
  hnb = qn1                                          &
       + qn2*tr2                                     &
       + (qn3 + qn4*tr + qn5*tr2 + qn6*tr3)*pr &
       + (qn7 + qn8*tr + qn9*tr2)*pr2            &
       + qn10*tr*pr3
  hnb = EXP(hnb)

  bHNO3 = bHNO3b * (1.0_dp - bH2SO4 / bH2SO4b)
  pHNO3 = bHNO3 / ((hnb*bHNO3 + hsb*bH2SO4)              &
         / (bHNO3+bH2SO4))
  IF (bH2SO4 > bH2SO4b .OR. pHNO3 > pHNO3_tot) THEN
    bHNO3 = 0.0_dp
  ENDIF

  mz_psc_liq_bHNO3 = bHNO3
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_liq_bHNO3
!=========================================================================




!============================================================================
 ELEMENTAL FUNCTION mz_psc_liq_partHNO3(TEMP, PRESS,          &
                                         H2O_gl, HNO3_gl,      &
                                         bHNO3, bH2SO4, bH2SO4b)
  !-------------------------------------------------------------------------
  ! This elemental function is one of the functions that calculate the
  ! composition of liquid stratospheric aerosol droplets according to
  ! Carslaw et al., GRL, 1995
  !-------------------------------------------------------------------------
  ! The HNO3/H2SO4 composition is based on a thermodynamic model of the
  ! HCl-HNO3-H2SO4-H2O-System (Carslaw et al, JPC, 1995).
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP = temperature / K
  ! PRESS    = total pressure / hPa
  ! H2O_gl   = amount-of-substance ratio of gaseous+liquid H2O to dry air
  !            / (mol/mol)
  ! HNO3_gl  = amount-of-substance ratio of gaseous+liquid HNO3 to dry air
  !            / (mol/mol)
  ! bHNO3    = molality of HNO3 in water  / (mol/kg)
  ! bH2SO4   = molality of H2SO4 in water / (mol/kg)
  ! bH2SO4b  = molality of H2SO4 in water for binary solution / (mol/kg)
  ! -------
  ! OUTPUT:
  ! -------
  ! partHNO3 = fraction of HNO3 remaining in gas phase after partitioning
  !            into liquid (1 = no removal from gas phase to the aerosol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! pr, pr2, pr3 = auxiliary variables
  ! tr, tr2, tr3 = auxiliary variables
  ! hnb = solubility of binary HNO3/H2O solution
  ! hsb = solubility of binary H2SO4/H2O solution
  ! qn, qs = auxiliary variables
  !-------------------------------------------------------------------------

   IMPLICIT NONE

   REAL(dp), INTENT(in) :: TEMP, PRESS
   REAL(dp), INTENT(in) :: H2O_gl, HNO3_gl
   REAL(dp), INTENT(in) :: bHNO3, bH2SO4, bH2SO4b
   REAL(dp) :: partHNO3, mz_psc_liq_partHNO3

   REAL(dp) :: pH2O_tot, pHNO3_tot
   REAL(dp) :: pHNO3

   ! help variables
   REAL(dp) :: pr, pr2, pr3
   REAL(dp) :: tr, tr2, tr3

   REAL(dp) :: hsb, hnb

   ! array declaration
   REAL(dp), PARAMETER :: &
        qn1 = 14.5734_dp, &
        qn2 = 0.0615994_dp, &
        qn3 = -1.14895_dp, &
        qn4 = 0.691693_dp, &
        qn5 = -0.098863_dp, &
        qn6 = 0.0051579_dp, &
        qn7 = 0.123472_dp, &
        qn8 = -0.115574_dp, &
        qn9 = 0.0110113_dp, &
        qn10 = 0.0097914_dp, &
        qs1 = 14.4700_dp, &
        qs2 = 0.0638795_dp, &
        qs3 = -3.29597_dp, &
        qs4 = 1.778224_dp, &
        qs5 = -0.223244_dp, &
        qs6 = 0.0086486_dp, &
        qs7 = 0.536695_dp, &
        qs8 = -0.335164_dp, &
        qs9 = 0.0265153_dp, &
        qs10 = 0.0157550_dp

   INTRINSIC exp, log, max, min

  IF (TEMP>215.0_dp) THEN
    !-----
    ! assume solution is pure H2SO4/H2O
    !-----
    mz_psc_liq_partHNO3 = 1.0_dp
    RETURN
  END IF

  !-----
  ! convert input
  !-----
  pH2O_tot  = min(max(PRESS*H2O_gl,2e-5_dp),2e-3_dp)/1013.25_dp
  pHNO3_tot = PRESS * min(HNO3_gl,2.0e-8_dp) / 1013.25_dp

  !-----
  ! initialise partitioning
  !-----
   partHNO3 = 1.0_dp

   tr   = 1.0e4_dp / max(TEMP,185.0_dp) - 43.4782608_dp
   pr   = LOG(pH2O_tot) + 18.4_dp
   tr2  = tr**2
   tr3  = tr**3
   pr2  = pr**2
   pr3  = pr**3

  !-----
  ! the HNO3/H2SO4/H2O solution composition
  !-----
    hsb = qs1                                            &
         + qs2*tr2                                       &
         + (qs3 + qs4*tr + qs5*tr2 + qs6*tr3)*pr   &
         + (qs7 + qs8*tr + qs9*tr2)*pr2              &
         + qs10*tr*pr3
    hsb = EXP(hsb)
    hnb = qn1                                          &
         + qn2*tr2                                     &
         + (qn3 + qn4*tr + qn5*tr2 + qn6*tr3)*pr &
         + (qn7 + qn8*tr + qn9*tr2)*pr2            &
         + qn10*tr*pr3
    hnb = EXP(hnb)
    pHNO3 = bHNO3 / ((hnb*bHNO3 + hsb*bH2SO4)              &
         / (bHNO3+bH2SO4))
    IF (bH2SO4 > bH2SO4b .OR. pHNO3 >= pHNO3_tot) THEN
      partHNO3 = 1.0_dp
    ELSE
      partHNO3 = pHNO3 / pHNO3_tot
    END IF

   mz_psc_liq_partHNO3 = min(max(partHNO3, 0.0_dp), 1.0_dp)
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_liq_partHNO3
!=========================================================================



!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_wnen(TEMP,                   &
                                    bHNO3, bH2SO4, bH2SO4b, &
                                    partHNO3)
  !-------------------------------------------------------------------------
  ! This elemental function is one of the functions that calculate the
  ! composition of liquid stratospheric aerosol droplets according to
  ! Carslaw et al., GRL, 1995
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP = temperature / K
  ! bHNO3    = molality of HNO3 in water  / (mol/kg)
  ! bH2SO4   = molality of H2SO4 in water / (mol/kg)
  ! bH2SO4b  = molality of H2SO4 in water for binary solution / (mol/kg)
  ! partHNO3 = fraction of HNO3 remaining in gas phase after partitioning
  !            into liquid (1 = no removal from gas phase to the aerosol)
  ! -------
  ! OUTPUT:
  ! -------
  ! wnen     = auxiliary variable
  !          = m(H2O)/m(H2O) + m(H2SO4)/m(H2O) + m(HNO3)/m(H2O)
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP
  REAL(dp), INTENT(in) :: bHNO3, bH2SO4, bH2SO4b
  REAL(dp), INTENT(in) :: partHNO3

  REAL(dp) :: wnen, mz_psc_liq_wnen

  IF (TEMP <= 215.0_dp) THEN
    !-----
    ! the HNO3/H2SO4/H2O solution composition
    !-----
    wnen = 1.0_dp + bH2SO4*MolMassH2SO4/1000.0_dp &
          + bHNO3*MolMassHNO3/1000.0_dp
    IF (bH2SO4 > bH2SO4b .OR. partHNO3 > 1.0_dp) THEN
      wnen = 1.0_dp + bH2SO4*MolMassH2SO4/1000.0_dp
    ENDIF
  ELSE
    !-----
    ! assume solution is pure H2SO4/H2O
    !-----
    wnen = 1.0_dp + bH2SO4*MolMassH2SO4/1000.0_dp
  ENDIF

  mz_psc_liq_wnen = wnen
END FUNCTION mz_psc_liq_wnen
!=========================================================================





!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_wHNO3(TEMP,                   &
                                     bHNO3, bH2SO4, bH2SO4b, &
                                     partHNO3, wnen)
  !-------------------------------------------------------------------------
  ! This elemental function is one of the functions that calculate the
  ! composition of liquid stratospheric aerosol droplets according to
  ! Carslaw et al., GRL, 1995
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP = temperature / K
  ! bHNO3    = molality of HNO3 in water  / (mol/kg)
  ! bH2SO4   = molality of H2SO4 in water / (mol/kg)
  ! bH2SO4b  = molality of H2SO4 in water for binary solution / (mol/kg)
  ! partHNO3 = fraction of HNO3 remaining in gas phase after partitioning
  !            into liquid (1 = no removal from gas phase to the aerosol)
  ! wnen     = auxiliary variable
  ! -------
  ! OUTPUT:
  ! -------
  ! wHNO3    = mass fraction of HNO3 in the liquid aerosol / (kg/kg)
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP
  REAL(dp), INTENT(in) :: bHNO3, bH2SO4, bH2SO4b
  REAL(dp), INTENT(in) :: partHNO3, wnen

  REAL(dp) :: wHNO3, mz_psc_liq_wHNO3

  IF (TEMP <= 215.0_dp) THEN
    !-----
    ! the HNO3/H2SO4/H2O solution composition
    !-----
    wHNO3 = bHNO3 * MolMassHNO3 / (1000.0_dp*wnen)
    IF (bH2SO4 > bH2SO4b .OR. partHNO3 > 1.0_dp) THEN
      wHNO3 = 0.0_dp
    ENDIF
  ELSE
    !-----
    ! assume solution is pure H2SO4/H2O
    !-----
    wHNO3 = 0.0_dp
  ENDIF

  mz_psc_liq_wHNO3 = wHNO3
END FUNCTION mz_psc_liq_wHNO3
!=========================================================================



!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_wH2SO4(TEMP,            &
                                      bH2SO4, bH2SO4b, &
                                      partHNO3, wnen)
  !-------------------------------------------------------------------------
  ! This elemental function is one of the functions that calculate the
  ! composition of liquid stratospheric aerosol droplets according to
  ! Carslaw et al., GRL, 1995
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP = temperature / K
  ! bH2SO4   = molality of H2SO4 in water / (mol/kg)
  ! bH2SO4b  = molality of H2SO4 in water for binary solution / (mol/kg)
  ! partHNO3 = fraction of HNO3 remaining in gas phase after partitioning
  !            into liquid (1 = no removal from gas phase to the aerosol)
  ! wnen     = auxiliary variable
  ! -------
  ! OUTPUT:
  ! -------
  ! wH2SO4   = mass fraction of H2SO4 in the liquid aerosol / (kg/kg)
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP
  REAL(dp), INTENT(in) :: bH2SO4, bH2SO4b
  REAL(dp), INTENT(in) :: partHNO3, wnen

  REAL(dp) :: wH2SO4, mz_psc_liq_wH2SO4

  IF (TEMP <= 215.0_dp) THEN
    !-----
    ! the HNO3/H2SO4/H2O solution composition
    !-----
    wH2SO4 = bH2SO4 * MolMassH2SO4/1000.0_dp / wnen
    IF (bH2SO4 > bH2SO4b .OR. partHNO3 > 1.0_dp) THEN
      wH2SO4 = bH2SO4b * MolMassH2SO4/1000.0_dp / wnen
    ENDIF
  ELSE
    !-----
    ! assume solution is pure H2SO4/H2O
    !-----
    wH2SO4 = bH2SO4b * MolMassH2SO4/1000.0_dp / wnen
  ENDIF

  mz_psc_liq_wH2SO4 = wH2SO4
END FUNCTION mz_psc_liq_wH2SO4
!=========================================================================




!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_hHCl(TEMP, wHNO3, wH2SO4, wnen)
  !-------------------------------------------------------------------------
  ! The HCl solubility parameterisation is taken from Luo et al., GRL, 1995.
  ! It is assumed that HCl is a trace component of the aerosol.
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP = temperature / K
  ! wHNO3    = mass fraction of HNO3 in the liquid aerosol / (kg/kg)
  ! wH2SO4   = mass fraction of H2SO4 in the liquid aerosol / (kg/kg)
  ! wnen     = auxiliary variable
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_liq_hHCl = effective HCl henrys law constant in 
  !                   HNO3-H2SO4-H2O solution / (mol/(kg*atm))
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! sqrtwH2SO4 = sqrt(wH2SO4)
  ! sqrtwHNO3 = sqrt(wHNO3)
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP
  REAL(dp), INTENT(in) :: wHNO3, wH2SO4, wnen

  REAL(dp) :: ln_hHCl, mz_psc_liq_hHCl

  REAL(dp) :: sqrtwH2SO4, sqrtwHNO3

  INTRINSIC log, exp, max, min, sqrt

  sqrtwH2SO4 = SQRT(wH2SO4)
  sqrtwHNO3 = SQRT(wHNO3)

  ln_hHCl =-(21.0_dp + 46.61_dp*wHNO3 + 4.069_dp*wH2SO4                     &
           - 4.837_dp*sqrtwHNO3 + 2.186_dp*sqrtwH2SO4                       &
           - 63.0_dp*wHNO3**2 - 40.17_dp*wHNO3*wH2SO4 - 1.571_dp*wH2SO4**2) &
           - 1.0_dp/min(max(185.0_dp,TEMP),235.0_dp)                        &
           * (-7437.0_dp - 8327.8_dp*wHNO3 + 1300.9_dp*wH2SO4               &
           + 1087.2_dp*sqrtwHNO3 - 242.71_dp*sqrtwH2SO4                     &
           + 18749.0_dp*wHNO3**2 + 18500.0_dp*wHNO3*wH2SO4                  &
           + 5632.0_dp*wH2SO4**2)                                           &
           - LOG(wHNO3 + 0.61_dp*wH2SO4)                                    &
           - LOG(MolMassHCl / (1000.0_dp*wnen))

  mz_psc_liq_hHCl = EXP(ln_hHCl) * 1.013e3_dp
END FUNCTION mz_psc_liq_hHCl
!=========================================================================




!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_hHBr(TEMP, wHNO3, wH2SO4, wnen)
  !-------------------------------------------------------------------------
  ! The HBr solubility parameterisation is taken from Luo et al., GRL, 1995.
  ! It is assumed that HBr is a trace component of the aerosol.
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP = temperature / K
  ! wHNO3    = mass fraction of HNO3 in the liquid aerosol / (kg/kg)
  ! wH2SO4   = mass fraction of H2SO4 in the liquid aerosol / (kg/kg)
  ! wnen     = auxiliary variable
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_liq_hHBr = effective HBr henrys law constant in 
  !                   HNO3-H2SO4-H2O solution / (mol/(kg*atm))
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! sqrtwH2SO4 = sqrt(wH2SO4)
  ! sqrtwHNO3 = sqrt(wHNO3)
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP
  REAL(dp), INTENT(in) :: wHNO3, wH2SO4, wnen

  REAL(dp) :: ln_hHBr, mz_psc_liq_hHBr

  REAL(dp) :: sqrtwH2SO4, sqrtwHNO3

  INTRINSIC log, exp, max, min, sqrt

  sqrtwH2SO4 = SQRT(wH2SO4)
  sqrtwHNO3 = SQRT(wHNO3)

  ln_hHBr = - (17.83_dp + 1.02_dp*wHNO3 - 1.08_dp*wH2SO4                  &
               + 3.9_dp*sqrtwHNO3 + 4.38_dp*sqrtwH2SO4                    &
               - 8.87_dp*wHNO3**2 - 17.0_dp*wHNO3*wH2SO4                  &
               + 3.73_dp*wH2SO4**2)                                       &
         - 1.0_dp/min(max(185.0_dp,TEMP),235.0_dp)                        &
                      *(- 8220.5_dp - 362.76_dp*wHNO3 + 658.93_dp*wH2SO4  &
                        - 914.0_dp*sqrtwHNO3 - 955.3_dp*sqrtwH2SO4        &
                        + 9976.6_dp*wHNO3**2 + 19778.5_dp*wHNO3*wH2SO4    &
                        + 7680.0_dp*wH2SO4**2)                            &
         - LOG(wHNO3 + 0.41_dp*wH2SO4)                                    &
         - LOG(MolMassHBr/(1000.0_dp*wnen))

  mz_psc_liq_hHBr = EXP(ln_hHBr) * 1.013e3_dp
END FUNCTION mz_psc_liq_hHBr
!=========================================================================




!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_hHOCl(TEMP, bH2SO4, bHNO3)
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  ! The HOCl solubility parameterisation is taken from Huthwelker et al.,
  ! JAC, 1995. 
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP = temperature / K
  ! bH2SO4   = molality of H2SO4 in water / (mol/kg)
  ! bHNO3    = molality of HNO3 in water  / (mol/kg)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_liq_hHOCl = effective HOCl henrys law constant in 
  !                    HNO3-H2SO4-H2O solution / (mol/(kg*atm))
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP
  REAL(dp), INTENT(in) :: bH2SO4, bHNO3

  REAL(dp) :: ln_hHOCl, mz_psc_liq_hHOCl

  INTRINSIC exp, max

  ln_hHOCl = 6.4946_dp                                            &
            - (-0.04107_dp + 54.56_dp/TEMP) * max(bH2SO4, bHNO3)  &
            - 5862.0_dp * (1.0_dp/298.15_dp - 1.0_dp/TEMP)

  mz_psc_liq_hHOCl = exp(ln_hHOCl)
END FUNCTION mz_psc_liq_hHOCl
!=========================================================================



!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_hHOBr(hHOCl)
  !-------------------------------------------------------------------------
  ! The solubility of HOBr is not known for al stratospheric conditions. 
  ! Limited Data (Hanson and Ravishankara, 1995, at 210 K, mass fraction of 
  ! H2SO4 = 60 %, indicate that hHOBr = ca. 18*hHOCl. The solubility of 
  ! HOBr is low enough to ignore gas phase removal.
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! hHOCl    = effective HOCl henrys law constant in HNO3-H2SO4-H2O solution
  !            / (mol/(kg*atm))
  ! -------
  ! OUTPUT:
  ! -------
  ! hHOBr    = effective HOBr henrys law constant in HNO3-H2SO4-H2O solution
  !            / (mol/(kg*atm))
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: hHOCl

  REAL(dp) :: mz_psc_liq_hHOBr

  mz_psc_liq_hHOBr = 18.0_dp * hHOCl

END FUNCTION mz_psc_liq_hHOBr
!=========================================================================



!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_partHBr(PRESS,            &
                                       HBr_gl, H2SO4_gl, &
                                       hHBr, bH2SO4)
  !-------------------------------------------------------------------------
  ! The equation for partHBr can be derived with the approximation that
  ! all H2SO4 is in the liquid phase. 
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! PRESS    = total pressure / hPa
  ! HBr_gl   = amount-of-substance ratio of gaseous+liquid HBr to dry air 
  !            / (mol/mol)
  ! H2SO4_gl = amount-of-substance ratio of gaseous+liquid H2SO4 to dry air 
  !            / (mol/mol)
  ! hHBr     = effective HBr henrys law constant in HNO3-H2SO4-H2O solution
  !            / (mol / (kg*atm))
  ! bH2SO4   = molality of H2SO4 in water / (mol/kg)
  ! -------
  ! OUTPUT:
  ! -------
  ! partHBr  = fraction of HBr remaining in gas phase after partitioning
  !            into liquid (1 = no removal from gas phase to aerosols)
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: PRESS
  REAL(dp), INTENT(in) :: HBr_gl, H2SO4_gl
  REAL(dp), INTENT(in) :: hHBr, bH2SO4

  REAL(dp) :: partHBr, mz_psc_liq_partHBr
  
  INTRINSIC max, min

  IF (hHBr>0.0_dp .AND. HBr_gl>0.0_dp) THEN
    partHBr = 1013.25_dp*bH2SO4 &
             / (H2SO4_gl*hHBr*PRESS + 1013.25_dp*bH2SO4)
  ELSE
      partHBr = 1.0_dp
  END IF

  mz_psc_liq_partHBr = min(max(partHBr, 0.0_dp), 1.0_dp)
END FUNCTION mz_psc_liq_partHBr
!=========================================================================



!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_wHCl(PRESS,            &
                                    HCl_gl, H2SO4_gl, &
                                    hHCl, bH2SO4, wnen)
  !-------------------------------------------------------------------------
  ! The equation for wHCl can be derived with the approximation that
  ! all H2SO4 is in the liquid phase. 
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! PRESS    = total pressure / hPa
  ! HCl_gl   = amount-of-substance ratio of gaseous+liquid HCl to dry air 
  !            / (mol/mol)
  ! H2SO4_gl = amount-of-substance ratio of gaseous+liquid H2SO4 to dry air 
  !            / (mol/mol)
  ! hHCl     = effective HCl henrys law constant in HNO3-H2SO4-H2O solution
  !            / (mol/(kg*atm))
  ! bH2SO4   = molality of H2SO4 in water / (mol/kg)
  ! wnen     = auxiliary variable
  ! -------
  ! OUTPUT:
  ! -------
  ! wHCl     = mass fraction of HCl in the liquid aerosol / (kg/kg)
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: PRESS
  REAL(dp), INTENT(in) :: HCl_gl, H2SO4_gl
  REAL(dp), INTENT(in) :: hHCl, bH2SO4, wnen

  REAL(dp) :: wHCl, mz_psc_liq_wHCl

  IF (hHCl>0.0_dp .AND. HCl_gl>0.0_dp) THEN
    wHCl = PRESS*HCL_gl*hHCl*bH2SO4*MolMassHCl   &
          /(1000.0_dp*wnen*(H2SO4_gl*hHCl*PRESS + 1013.25_dp*bH2SO4))
  ELSE
    wHCl = 0.0_dp
  END IF

  mz_psc_liq_wHCl = wHCl
END FUNCTION mz_psc_liq_wHCl
!=========================================================================



!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_wHOCl(PRESS,             &
                                     HOCl_gl, H2SO4_gl, &
                                     hHOCl, bH2SO4, wnen)
  !-------------------------------------------------------------------------
  ! The equation for wHOCl can be derived with the approximation that
  ! all H2SO4 is in the liquid phase. 
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! PRESS    = total pressure / hPa
  ! HOCl_gl  = amount-of-substance ratio of gaseous+liquid HOCl to dry air
  !            / (mol/mol)
  ! H2SO4_gl = amount-of-substance ratio of gaseous+liquid H2SO4 to dry air 
  !            / (mol/mol)
  ! hHOCl    = effective HOCl henrys law constant in HNO3-H2SO4-H2O solution
  !            / (mol/(kg*atm))
  ! bH2SO4   = molality of H2SO4 in water / (mol/kg)
  ! wnen     = auxiliary variable
  ! -------
  ! OUTPUT:
  ! -------
  ! wHOCl    = mass fraction of HOCl in the liquid aerosol / (kg/kg)
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: PRESS
  REAL(dp), INTENT(in) :: HOCl_gl, H2SO4_gl
  REAL(dp), INTENT(in) :: hHOCl, bH2SO4, wnen

  REAL(dp) :: wHOCl, mz_psc_liq_wHOCl

  IF (hHOCl>0.0_dp .AND. HOCl_gl>0.0_dp) THEN
    wHOCl = PRESS*HOCl_gl*hHOCl*bH2SO4*MolMassHOCl   &
           /((1000.0_dp*wnen)*(H2SO4_gl*PRESS*hHOCl + 1013.25_dp*bH2SO4))
  ELSE
    wHOCl = 0.0_dp
  END IF

  mz_psc_liq_wHOCl = wHOCl
END FUNCTION mz_psc_liq_wHOCl
!=========================================================================



!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_wHOBr(PRESS,             &
                                     HOBr_gl, H2SO4_gl, &
                                     hHOBr, bH2SO4, wnen)
  !-------------------------------------------------------------------------
  ! The equation for wHOBr can be derived with the approximation that
  ! all H2SO4 is in the liquid phase. 
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! PRESS    = total pressure / hPa
  ! HOBr_gl  = amount-of-substance ratio of gaseous+liquid HOBr to dry air
  !            / (mol/mol)
  ! H2SO4_gl = amount-of-substance ratio of gaseous+liquid H2SO4 to dry air 
  !            / (mol/mol)
  ! hHOBr    = effective HOBr henrys law constant in HNO3-H2SO4-H2O solution
  !            / (mol/(kg*atm))
  ! bH2SO4   = molality of H2SO4 in water / (mol/kg)
  ! wnen     = auxiliary variable
  ! -------
  ! OUTPUT:
  ! -------
  ! wHOBr    = mass fraction of HOBr in the liquid aerosol / (kg/kg)
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: PRESS
  REAL(dp), INTENT(in) :: HOBr_gl, H2SO4_gl
  REAL(dp), INTENT(in) :: hHOBr, bH2SO4, wnen

  REAL(dp) :: wHOBr, mz_psc_liq_wHOBr

  IF (hHOBr>0.0_dp .AND. HOBr_gl>0.0_dp) THEN
    wHOBr = PRESS*HOBr_gl*hHOBr*bH2SO4*MolMassHOBr   &
           /((1000.0_dp*wnen)*(H2SO4_gl*PRESS*hHOBr + 1013.25_dp*bH2SO4))
  ELSE
    wHOBr = 0.0_dp
  END IF

  mz_psc_liq_wHOBr = wHOBr
END FUNCTION mz_psc_liq_wHOBr
!=========================================================================


!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_H2O(wH2SO4, wHNO3, H2SO4_l, H2O_t)
  !-------------------------------------------------------------------------
  ! The function calculates the amount-of-substance ratio of H2O in 
  ! liquid aerosol particles to dry air. It is explicitely avoided that 
  ! more H2O than available ends up in the liquid phase. 
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! wH2SO4  = mass fraction of H2SO4 in the liquid aerosol / (kg/kg)
  ! wHNO3   = mass fraction of HNO3 in the liquid aerosol / (kg/kg)
  ! H2SO4_l = amount-of-substance ratio of liquid H2SO4 to dry air
  !          / (mol/mol)
  ! H2O_t  = amount-of-substance ratio of total H2O to dry air
  !          / (mol/mol)
  ! -------
  ! OUTPUT:
  ! -------
  ! H2O_l   = amount-of-substance ratio of liquid phase H2O to dry air 
  !          / (mol/mol)
  !-------------------------------------------------------------------------
  
  IMPLICIT NONE
  
  REAL(dp), INTENT(in) :: wH2SO4, wHNO3, H2SO4_l, H2O_t
  
  REAL(dp) :: H2O_l, mz_psc_liq_H2O
  
  INTRINSIC min 
  
  H2O_l   = (1.0_dp-wH2SO4-wHNO3)          &
           *(H2SO4_l/wH2SO4                &
           *MolMassH2SO4)/MolMassH2O
  
  mz_psc_liq_H2O = min(H2O_l, H2O_t)
END FUNCTION mz_psc_liq_H2O
!=========================================================================


!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_HCl(wHCl, H2SO4_l, wH2SO4, HCl_t)
  !-------------------------------------------------------------------------
  ! The function calculates the amount-of-substance fraction of HCl in 
  ! liquid aerosol particles to dry air. It is explicitely avoided that 
  ! more HCl than available ends up in the liquid phase. 
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! wHCl    = mass fraction of HCl in the liquid aerosol / (kg/kg)
  ! H2SO4_l = amount-of-substance ratio of liquid H2SO4 to dry air
  !          / (mol/mol)
  ! wH2SO4  = mass fraction of H2SO4 in the liquid aerosol / (kg/kg)
  ! HCl_t   = amount-of-substance ratio of total H2O to dry air
  !          / (mol/mol)
  ! -------
  ! OUTPUT:
  ! -------
  ! HCl_t   = amount-of-substance ratio of liquid phase HCl to dry air 
  !          / (mol/mol)
  !-------------------------------------------------------------------------

  IMPLICIT NONE
  
  REAL(dp), INTENT(in) :: wHCl, H2SO4_l, wH2SO4, HCl_t
  
  REAL(dp) :: HCl_l, mz_psc_liq_HCl
  
  INTRINSIC min 
  
  HCl_l   = wHCl*MolMassH2SO4*H2SO4_l/(wH2SO4*MolMassHCl)
  
  mz_psc_liq_HCl = min(HCl_l, HCl_t)
END FUNCTION mz_psc_liq_HCl
!=========================================================================


!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_HOCl(wHOCl, H2SO4_l, wH2SO4, HOCl_t)
  !-------------------------------------------------------------------------
  ! The function calculates the amount-of-substance ratio of HOCl in 
  ! liquid aerosol particles to dry air. It is explicitely avoided that 
  ! more HOCl than available ends up in the liquid phase. 
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! wHOCl   = mass fraction of HOCl in the liquid aerosol / (kg/kg)
  ! H2SO4_l = amount-of-substance ratio of liquid H2SO4 to dry air
  !          / (mol/mol)
  ! wH2SO4  = mass fraction of H2SO4 in the liquid aerosol / (kg/kg)
  ! HOCl_t  = amount-of-substance ratio of total HOCl to dry air
  !          / (mol/mol)
  ! -------
  ! OUTPUT:
  ! -------
  ! HOCl_t  = amount-of-substance ratio of liquid phase HOCl to dry air 
  !          / (mol/mol)
  !-------------------------------------------------------------------------

  IMPLICIT NONE
  
  REAL(dp), INTENT(in) :: wHOCl, H2SO4_l, wH2SO4, HOCl_t
  
  REAL(dp) :: HOCl_l, mz_psc_liq_HOCl
  
  INTRINSIC min 
  
  HOCl_l   = wHOCl*MolMassH2SO4*H2SO4_l/(wH2SO4*MolMassHOCl)
  
  mz_psc_liq_HOCl = min(HOCl_l, HOCl_t)
END FUNCTION mz_psc_liq_HOCl
!=========================================================================


!=========================================================================
ELEMENTAL FUNCTION mz_psc_liq_HOBr(wHOBr, H2SO4_l, wH2SO4, HOBr_t)
  !-------------------------------------------------------------------------
  ! The function calculates the amount-of-substance ratio of HOBr in 
  ! liquid aerosol particles to dry air. It is explicitely avoided that 
  ! more HOBr than available ends up in the liquid phase. 
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! wHOBr    = mass fraction of HOBr in the liquid aerosol / (kg/kg)
  ! H2SO4_l = amount-of-substance ratio of liquid H2SO4 to dry air
  !          / (mol/mol)
  ! wH2SO4  = mass fraction of H2SO4 in the liquid aerosol / (kg/kg)
  ! HOBr_t   = amount-of-substance ratio of total H2O to dry air
  !          / (mol/mol)
  ! -------
  ! OUTPUT:
  ! -------
  ! HOBr_t   = amount-of-substance ratio of liquid phase HOBr to dry air 
  !          / (mol/mol)
  !-------------------------------------------------------------------------

  IMPLICIT NONE
  
  REAL(dp), INTENT(in) :: wHOBr, H2SO4_l, wH2SO4, HOBr_t
  
  REAL(dp) :: HOBr_l, mz_psc_liq_HOBr
  
  INTRINSIC min 
  
  HOBr_l   = wHOBr*MolMassH2SO4*H2SO4_l/(wH2SO4*MolMassHOBr)
  
  mz_psc_liq_HOBr = min(HOBr_l, HOBr_t)
END FUNCTION mz_psc_liq_HOBr
!=========================================================================


!==========================================================================
 ELEMENTAL FUNCTION mz_psc_N_solid(phase,             &
                                   TEMP, PRESS,       &
                                   H2O_ice, HNO3_nat)
  !-------------------------------------------------------------------------
  ! This function calculates the particle number density of solid, spherical
  ! particles containing ice and/or NAT.
  ! If there is much ice/nat in a grid box, the function returns a constant 
  ! minimum value N_solid_max. However, if the input particle number density
  ! together with a rather small amount of ice/nat in the grid box leads to
  ! very small particle radi, the particle number density is decreased, so
  ! that the particle radius equals a default minimum radius.
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! phase    = phase indicator (1=liquid, 2=liquid+nat, 3=liquid+nat+ice)
  ! TEMP = temperature / K
  ! PRESS    = pressure / hPa
  ! H2O      = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! HNO3_nat = amount-of-substance ratio of nat phase HNO3 to dry air
  !            / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! m_ice    = mass of ice per volume of air / (kg/m**3)
  ! m_nat    = mass of NAT per volume of air / (kg/m**3)
  ! v_ice    = volume fraction of ice in air / (m**3/m**3)
  ! v_nat    = volume fraction of NAT in air / (m**3/m**3)
  ! m_solid  = total mass of solid particles per volume of air / (kg/m**3)
  ! v_solid  = volume fraction of particles in air / (m**3/m**3)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_N_solid = number of solid particles per m**3 air
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in) :: phase
  REAL(dp), INTENT(in) :: TEMP, PRESS, H2O_ice, HNO3_nat

  REAL(dp) :: mz_psc_N_solid

  REAL(dp) :: m_solid, v_solid
  REAL(dp) :: m_ice, v_ice
  REAL(dp) :: m_nat, v_nat
  REAL(dp) :: factor

  INTRINSIC MIN

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/m**3)
  ! -----
  factor   = mix2conc(PRESS, TEMP)*1.0E6

  IF (phase==3) THEN
    m_ice   =(H2O_ice*factor) * MolMassH2O/(N_Avogadro*1000.0_dp)
    v_ice   = m_ice / dens_ice
    m_nat   =(HNO3_nat*factor) * (MolMassHNO3+3.0_dp*MolMassH2O) &
             /(N_Avogadro*1000.0_dp)
    v_nat   = m_nat / dens_nat
    m_solid = m_ice+m_nat
    v_solid = v_ice+v_nat
    mz_psc_N_solid = MIN(3.0_dp*v_solid/(4.0_dp*pi*r_min**3), N_solid_max)
  ELSEIF (phase==2) THEN
    m_solid =(HNO3_nat*factor) * (MolMassHNO3+3.0_dp*MolMassH2O) &
              /(N_Avogadro*1000.0_dp)
    v_solid = m_solid / dens_nat
    mz_psc_N_solid = MIN(3.0_dp*v_solid/(4.0_dp*pi*r_min**3), N_solid_max)
  ELSE   ! should not happen
    mz_psc_N_solid = N_solid_max
  END IF
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_N_solid
!=========================================================================



!==========================================================================
ELEMENTAL FUNCTiON mz_psc_r_solid(phase,             &
                                  TEMP, PRESS,       &
                                  N_solid,           &
                                  H2O_ice, HNO3_nat)
  !-------------------------------------------------------------------------
  ! This function calculates the radius of solid particles containing
  ! ice and/or NAT.
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! phase    = phase indicator (1=liquid, 2=liquid+nat, 3=liquid+nat+ice)
  ! TEMP = temperature / K
  ! PRESS    = pressure / hPa
  ! N_solid  = number of solid particles per m**3 air
  ! H2O      = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! HNO3_nat = amount-of-substance ratio of nat phase HNO3 to dry air
  !            / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! m_ice = mass of ice per volume of air / (kg/m**3)
  ! m_nat = mass of NAT per volume of air / (kg/m**3)
  ! v_ice = volume fraction of ice in air / (m**3/m**3)
  ! v_nat = volume fraction of NAT in air / (m**3/m**3)
  ! m_solid = total mass of solid particles per volume of air / (kg/m**3)
  ! v_solid = volume fraction of particles in air / (m**3/m**3)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_r_solid = r_solid = radius of one single solid particle / m
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in) :: phase
  REAL(dp), INTENT(in) :: TEMP, PRESS, H2O_ice, HNO3_nat
  REAL(dp), INTENT(in) :: N_solid

  REAL(dp) :: mz_psc_r_solid, r_solid

  REAL(dp) :: m_solid, v_solid
  REAL(dp) :: m_ice, v_ice
  REAL(dp) :: m_nat, v_nat
  REAL(dp) :: factor
  
  INTRINSIC max, tiny

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/m**3)
  ! -----
  factor   = mix2conc(PRESS, TEMP)*1.0E6

  IF (phase==3) THEN
    m_ice   =(H2O_ice*factor) * MolMassH2O/(N_Avogadro*1000.0_dp)
    v_ice   = m_ice / dens_ice
    m_nat   =(HNO3_nat*factor) * (MolMassHNO3+3.0_dp*MolMassH2O) &
             /(N_Avogadro*1000.0_dp)
    v_nat   = m_nat / dens_nat
    m_solid = m_ice+m_nat
    v_solid = v_ice+v_nat
    r_solid = (3.0_dp*v_solid &
             /max(4.0_dp*N_solid*pi,tiny(1.0_dp)))**(1.0_dp/3.0_dp)
  ELSEIF (phase==2) THEN
    m_solid =(HNO3_nat*factor) * (MolMassHNO3+3.0_dp*MolMassH2O) &
              /(N_Avogadro*1000.0_dp)
    v_solid = m_solid / dens_nat
    r_solid = (3.0_dp*v_solid &
             /max(4.0_dp*N_solid*pi,tiny(1.0_dp)))**(1.0_dp/3.0_dp)
  ELSE ! should not happen
    r_solid = 0.0_dp
  END IF
  mz_psc_r_solid = r_solid
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_r_solid
!=========================================================================


!=========================================================================
ELEMENTAL FUNCTION mz_psc_surface_liquid(TEMP, PRESS, &
                                         H2O_l, HNO3_l, H2SO4_l)
  !-------------------------------------------------------------------------
  ! This function calculates the surface density in cm**2/cm**3
  ! at a certain temperature (TEMP / K) and pressure (PRESS / hPa) 
  ! and at a given amount-of-substance ratio 
  ! of liquid HX to dry air (HX_l / (mol/mol))
  !
  ! From Grainger et al. (1995) a relation between the numerical value of 
  ! the aerosol surface area density in square micrometer per cubic centimeter
  ! and the numerical value of the aerosol volume density in cubic micrometer
  ! per cubic centimeter is used:
  !   S = 8.406 * V**0.751
  ! Note that great care has to be taken in the usage of units (as always if
  ! numerical-value equations are used instead of quantity equations).
  !-------------------------------------------------------------------------
  ! INPUT:
  ! ------
  ! TEMP       = temperature / K
  ! PRESS      = pressure / hPa
  ! HX         = amount-of-substance ratio of liquid phase HX to dry air 
  !              / (mol/mol)
  !              Note: HX is H2O, HNO3, or H2SO4.
  ! ----------
  ! PARAMETER:
  ! ----------
  ! factor        = conversion factor 
  !                 (from amount-of-substance ratio / (mol/mol)
  !                  to particle number concentration / (1/cm**3))
  ! zhx_l         = particle number concentration of HX_l / (1/cm**3)
  !                 Note: zhx_l=factor*HX_l
  ! wH2SO4, wHNO3 = mass fraction of H2SO4, HNO3 / (kg/kg)
  ! dens          = mass density of the liquid solution / (g/cm**3)
  ! mtot          = mass of liquid per volume of air / (g/cm**3)
  ! v_liq         = volume fraction of liquid aerosol in air / (cm**3/cm**3)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_surface_liquid = mean surface of liquid droplets / (cm**2/cm**3)
  !-------------------------------------------------------------------------

   IMPLICIT NONE

   REAL(dp), INTENT(in) :: TEMP, PRESS
   REAL(dp), INTENT(in) :: H2O_l
   REAL(dp), INTENT(in) :: HNO3_l
   REAL(dp), INTENT(in) :: H2SO4_l

   REAL(dp) :: mz_psc_surface_liquid

   REAL(dp) :: mtot, v_liq
   REAL(dp) :: factor
   REAL(dp) :: dens
   REAL(dp) :: wH2SO4, wHNO3
   REAL(dp) :: zH2O_l, zHNO3_l, zH2SO4_l
   REAL(dp) :: factcm2

   factcm2=8.406_dp*1.0E-08_dp*(1.0E+12_dp)**0.751_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
   factor  = mix2conc(PRESS, TEMP)

   zH2O_l   = H2O_l   * factor
   zHNO3_l  = HNO3_l  * factor
   zH2SO4_l = H2SO4_l * factor

   ! -----
   ! total mass of all droplet per volume of air
   ! -----
   mtot=(zH2O_l*MolMassH2O+zHNO3_l*MolMassHNO3+zH2SO4_l*MolMassH2SO4)/ &
        N_Avogadro

   ! mz_rs_20060110+
!qqq+ PSC: no IF ... ELSE ... ENDIF, always
!     wH2SO4=(zH2SO4_l*MolMassH2SO4/N_Avogadro)/mtot
!     wHNO3=(zHNO3_l *MolMassHNO3 /N_Avogadro)/mtot
!qqq-
   IF (mtot>0._dp) THEN
     wH2SO4=(zH2SO4_l*MolMassH2SO4/N_Avogadro)/mtot
     wHNO3=(zHNO3_l *MolMassHNO3 /N_Avogadro)/mtot
   ELSE
     wH2SO4 = 0._dp
     wHNO3  = 0._dp
   ENDIF
   ! mz_rs_20060110-

   ! -----
   ! mass density of liquid solution calculated according to Luo et al. (1996)
   ! -----
   dens  = mz_psc_density(wH2SO4,wHNO3,TEMP)

   ! -----
   ! volume fraction of liquid aerosol in air / (cm**3/cm**3)
   ! -----

   v_liq = mtot/dens

   ! -----
   ! surface of liquid stratospheric aerosol / (cm**2/cm**3)
   ! -----

   mz_psc_surface_liquid = factcm2*v_liq**0.751

!-------------------------------------------------------------------------
 END FUNCTION mz_psc_surface_liquid
!=========================================================================



!=========================================================================
 ELEMENTAL FUNCTION mz_psc_dim_liquid(surface, sigmaaero)
  !-----------------------------------------------------------------------
  ! From Grainger et al. (1995) a relation between the numerical value of 
  ! the aerosol effective radius in micrometer and the  numerical value of 
  ! the aerosol volume density in cubic micrometer per cubic centimeter 
  ! is used: 
  !   reff=0.357*V^0.249
  ! Note that great care has to be taken in the usage of units (as always if
  ! numerical-value equations are used instead of quantity equations).
  ! 
  ! The surface median is then calculated via the equation
  !   r_surfmed=reff*exp{-0.5*(ln sigmaaero)**2}
  !
  ! For the input parameter sigmaaero a value of 1.8 is recommended. 
  !-------------------------------------------------------------------------
  ! INPUT:
  ! ------
  ! surface    = total liquid aerosol surface / (cm**2/cm**3)
  ! sigmaaero  = width of lognormal distribution
  ! ----------
  ! PARAMETER:
  ! ----------
  ! factexp       = conversion exponent
  ! factcm2       = conversion factor 
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_dim_liquid = surface median radius of liquid droplets / cm
  !-------------------------------------------------------------------------

   IMPLICIT NONE

   REAL(dp), INTENT(in) :: surface, sigmaaero

   REAL(dp) :: mz_psc_dim_liquid

   REAL(dp), PARAMETER :: factexp=0.249_dp/0.751_dp
   REAL(dp) :: factcm2
   
   INTRINSIC exp, log

   factcm2=0.357_dp*1.E-04_dp*(1.E+08_dp/8.406_dp)**factexp

   ! -----
   ! surface median of liquid stratospheric aerosol [cm]
   ! -----

   mz_psc_dim_liquid = &
     factcm2*surface**factexp*exp(-0.5_dp*log(sigmaaero)**2)
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_dim_liquid
!=========================================================================



!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_nat1(TEMP, PRESS,          &
                                   RnatMeter, NnatMeter, &
                                   HCl)
  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on nat particles.
  !-------------------------------------------------------------------------
  !--------------- 1) HCl+ClONO2->HNO3+Cl2 --------------
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RnatMeter = Radius of NAT particles / m
  ! NnatMeter = number density of particles in air / (1/m**3)
  ! HCl       = amount-of-substance ratio of gas phase HCl to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rnat     = Radius of NAT particles / cm
  ! Nnat     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHCl     = number concentration of gaseous HCl molecules / (1/cm**3)
  ! gamma1   = Uptake coefficient (unit 1)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_nat1 = khet1  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RnatMeter, NnatMeter
  REAL(dp), INTENT(in) :: HCl
  REAL(dp) :: khet1, mz_psc_het_nat1, gamma1

  REAL(dp) :: Rnat, Nnat
  REAL(dp) :: factor
  REAL(dp) :: cHCl

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rnat=RnatMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nnat=NnatMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHCl  = factor * HCl

  gamma1 = 0.2_dp     !JPL2000
  khet1  = gamma1 *4.56e4_dp*SQRT(TEMP/97.46_dp)*Rnat**2*Nnat/ &
           (1.0_dp + 3.3e4_dp*gamma1*Rnat*PRESS/TEMP)/(cHCl+1.0_dp)
  khet1  = max(khet1, minKhet)
  khet1  = min(khet1, maxKhet)

  mz_psc_het_nat1 = khet1
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_het_nat1
!=========================================================================




!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_nat2(TEMP, PRESS,          &
                                   RnatMeter, NnatMeter, &
                                   H2O)
  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on nat particles.
  !-------------------------------------------------------------------------
  !--------------- 2) H2O+ClONO2->HNO3+HOCl -------------
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RnatMeter = Radius of NAT particles / m
  ! NnatMeter = number density of particles in air / (1/m**3)
  ! H2O       = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rnat     = Radius of NAT particles / cm
  ! Nnat     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cH2O     = number concentration of gaseous H2O molecules / (1/cm**3)
  ! gamma2   = Uptake coefficient (unit 1)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_nat2 = khet2  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RnatMeter, NnatMeter
  REAL(dp), INTENT(in) :: H2O

  REAL(dp) :: khet2, mz_psc_het_nat2, gamma2
  REAL(dp) :: Rnat, Nnat
  REAL(dp) :: factor
  REAL(dp) :: cH2O

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rnat=RnatMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nnat=NnatMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cH2O  = factor * H2O

  gamma2 = 0.004_dp   !JPL2000
  khet2  = gamma2 *4.56e4_dp*SQRT(TEMP/97.46_dp)*Rnat**2*Nnat/ &
             (1.0_dp + 3.3e4_dp*gamma2*Rnat*PRESS/TEMP)/(cH2O+1.0_dp)
  khet2  = max(khet2, minKhet)
  khet2  = min(khet2, maxKhet)

  mz_psc_het_nat2 = khet2
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_het_nat2
!=========================================================================




!==========================================================================
 ELEMENTAL FUNCTION mz_psc_het_nat3(TEMP, PRESS,          &
                                    RnatMeter, NnatMeter, &
                                    HCl)
  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on nat particles.
  !-------------------------------------------------------------------------
  !--------------- 3) HOCl+ HCl->H2O+Cl2 ----------------
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RnatMeter = Radius of NAT particles / m
  ! NnatMeter = number density of particles in air / (1/m**3)
  ! HCl       = amount-of-substance ratio of gas phase HCl to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rnat     = Radius of NAT particles / cm
  ! Nnat     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHCl     = number concentration of gaseous HCl molecules / (1/cm**3)
  ! gamma3   = Uptake coefficient (unit 1)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_nat3 = khet3  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RnatMeter, NnatMeter
  REAL(dp), INTENT(in) :: HCl

  REAL(dp) :: khet3, mz_psc_het_nat3, gamma3
  REAL(dp) :: Rnat, Nnat
  REAL(dp) :: factor
  REAL(dp) :: cHCl

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rnat=RnatMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nnat=NnatMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHCl  = factor * HCl

  gamma3 = 0.1_dp                              !JPL 1997/2000
  khet3  = gamma3*4.56e4_dp*SQRT(TEMP/52.46_dp)*Rnat**2*Nnat/ &
             (1.0_dp + 3.3e4_dp*gamma3*Rnat*PRESS/TEMP)/(cHCl+1.0_dp)
  khet3  = max(khet3, minKhet)
  khet3  = min(khet3, maxKhet)

  mz_psc_het_nat3 = khet3
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_het_nat3
!=========================================================================



!==========================================================================
 ELEMENTAL FUNCTION mz_psc_het_nat4(TEMP, PRESS,            &
                                    RnatMeter, NnatMeter, &
                                    HCl)
  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on nat particles.
  !-------------------------------------------------------------------------
  !--------------- 4) N2O5+HCl->ClNO2+HNO3 --------------
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RnatMeter = Radius of NAT particles / m
  ! NnatMeter = number density of particles in air / (1/m**3)
  ! HCl       = amount-of-substance ratio of gas phase HCl to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rnat     = Radius of NAT particles / cm
  ! Nnat     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHCl     = number concentration of gaseous HCl molecules / (1/cm**3)
  ! gamma4   = Uptake coefficient (unit 1)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_nat4 = khet4  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RnatMeter, NnatMeter
  REAL(dp), INTENT(in) :: HCl

  REAL(dp) :: khet4, mz_psc_het_nat4, gamma4
  REAL(dp) :: Rnat, Nnat
  REAL(dp) :: factor
  REAL(dp) :: cHCl

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rnat=RnatMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nnat=NnatMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHCl  = factor * HCl

  gamma4 = 0.003_dp  !JPL1997
  khet4  = gamma4*4.56e4_dp*SQRT(TEMP/108.0_dp)*Rnat**2*Nnat/ &
             (1.0_dp + 3.3e4_dp*gamma4*Rnat*PRESS/TEMP)/(cHCl+1.0_dp)
  khet4  = max(khet4, minKhet)
  khet4  = min(khet4, maxKhet)

  mz_psc_het_nat4 = khet4
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_het_nat4
!=========================================================================



!==========================================================================
 ELEMENTAL FUNCTION mz_psc_het_nat5(TEMP, PRESS,          &
                                    RnatMeter, NnatMeter, &
                                    H2O)
  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on nat particles.
  !-------------------------------------------------------------------------
  !--------------- 5) N2O5+H2O->2HNO3 -------------------
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RnatMeter = Radius of NAT particles / m
  ! NnatMeter = number density of particles in air / (1/m**3)
  ! H2O       = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rnat     = Radius of NAT particles / cm
  ! Nnat     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cH2O     = number concentration of gaseous H2O molecules / (1/cm**3)
  ! gamma5   = Uptake coefficient (unit 1)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_nat5 = khet5  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RnatMeter, NnatMeter
  REAL(dp), INTENT(in) :: H2O

  REAL(dp) :: khet5, mz_psc_het_nat5, gamma5
  REAL(dp) :: Rnat, Nnat
  REAL(dp) :: factor
  REAL(dp) :: cH2O

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rnat=RnatMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nnat=NnatMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cH2O  = factor * H2O

  gamma5 = 0.0004_dp !JPL2000
  khet5  = gamma5*4.56e4_dp*SQRT(TEMP/108.0_dp)*Rnat**2*Nnat/ &
             (1.0_dp + 3.3e4_dp*gamma5*Rnat*PRESS/TEMP)/(cH2O+1.0_dp)
  khet5  = max(khet5, minKhet)
  khet5  = min(khet5, maxKhet)

  mz_psc_het_nat5 = khet5
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_het_nat5
!=========================================================================



!==========================================================================
 ELEMENTAL FUNCTION mz_psc_het_nat6(TEMP, PRESS,            &
                                    RnatMeter, NnatMeter, &
                                    HBr)
  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on nat particles.
  !-------------------------------------------------------------------------
  !--------------- 6) ClONO2+HBr->HNO3+BrCl -------------
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RnatMeter = Radius of NAT particles / m
  ! NnatMeter = number density of particles in air / (1/m**3)
  ! HBr       = amount-of-substance ratio of gas phase HBr to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rnat     = Radius of NAT particles / cm
  ! Nnat     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHBr     = number concentration of gaseous HBr molecules / (1/cm**3)
  ! gamma6   = Uptake coefficient (unit 1)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_nat6 = khet6  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RnatMeter, NnatMeter
  REAL(dp), INTENT(in) :: HBr

  REAL(dp) :: khet6, mz_psc_het_nat6, gamma6
  REAL(dp) :: Rnat, Nnat
  REAL(dp) :: factor
  REAL(dp) :: cHBr

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rnat=RnatMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nnat=NnatMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHBr  = factor * HBr

  gamma6 = 0.3_dp !JPL1997 gamma>0.3
  khet6  = gamma6*4.56e4_dp*SQRT(TEMP/97.46_dp)*Rnat**2*Nnat/ &
             (1.0_dp + 3.3e4_dp*gamma6*Rnat*PRESS/TEMP)/(cHBr+1.0_dp)
  khet6  = max(khet6, minKhet)
  khet6  = min(khet6, maxKhet)

  mz_psc_het_nat6 = khet6
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_het_nat6
!=========================================================================



!==========================================================================
 ELEMENTAL FUNCTION mz_psc_het_nat7(TEMP, PRESS,          &
                                    RnatMeter, NnatMeter, &
                                    HCl)
  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on nat particles.
  !-------------------------------------------------------------------------
  !--------------- 7) BrONO2+HCl->HNO3+BrCl -------------
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RnatMeter = Radius of NAT particles / m
  ! NnatMeter = number density of particles in air / (1/m**3)
  ! HCl       = amount-of-substance ratio of gas phase HCl to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rnat     = Radius of NAT particles / cm
  ! Nnat     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHCl     = number concentration of gaseous HCl molecules / (1/cm**3)
  ! gamma7   = Uptake coefficient (unit 1)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_nat7 = khet7  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RnatMeter, NnatMeter
  REAL(dp), INTENT(in) :: HCl

  REAL(dp) :: khet7, mz_psc_het_nat7, gamma7
  REAL(dp) :: Rnat, Nnat
  REAL(dp) :: factor
  REAL(dp) :: cHCl

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rnat=RnatMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nnat=NnatMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHCl  = factor * HCl

  gamma7 = 0.3_dp   !Carslaw original model (guess)
  khet7  = gamma7*4.56e4_dp*SQRT(TEMP/142.0_dp)*Rnat**2*Nnat/ &
             (1.0_dp + 3.3e4_dp*gamma7*Rnat*PRESS/TEMP)/(cHCl+1.0_dp)
  khet7  = max(khet7, minKhet)
  khet7  = min(khet7, maxKhet)

  mz_psc_het_nat7 = khet7
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_het_nat7
!=========================================================================



!==========================================================================
 ELEMENTAL FUNCTION mz_psc_het_nat8(TEMP, PRESS,          &
                                    RnatMeter, NnatMeter, &
                                    HBr)
  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on nat particles.
  !-------------------------------------------------------------------------
  !--------------- 8) HBr+HOCl->H2O+BrCl ----------------
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RnatMeter = Radius of NAT particles / m
  ! NnatMeter = number density of particles in air / (1/m**3)
  ! HBr       = amount-of-substance ratio of gas phase HBr to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rnat     = Radius of NAT particles / cm
  ! Nnat     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHBr     = number concentration of gaseous HBr molecules / (1/cm**3)
  ! gamma8   = Uptake coefficient (unit 1)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_nat8 = khet8  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RnatMeter, NnatMeter
  REAL(dp), INTENT(in) :: HBr

  REAL(dp) :: khet8, mz_psc_het_nat8, gamma8
  REAL(dp) :: Rnat, Nnat
  REAL(dp) :: factor
  REAL(dp) :: cHBr

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rnat=RnatMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nnat=NnatMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHBr  = factor * HBr

  gamma8 = 0.3_dp !Carslaw original code (guess)
  khet8  = gamma8*4.56e4_dp*SQRT(TEMP/52.46_dp)*Rnat**2*Nnat/ &
             (1.0_dp + 3.3e4_dp*gamma8*Rnat*PRESS/TEMP)/(cHBr+1.0_dp)
  khet8  = max(khet8, minKhet)
  khet8  = min(khet8, maxKhet)

  mz_psc_het_nat8 = khet8
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_het_nat8
!=========================================================================



!==========================================================================
 ELEMENTAL FUNCTION mz_psc_het_nat9(TEMP, PRESS,          &
                                    RnatMeter, NnatMeter, &
                                    HCl)
  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on nat particles.
  !-------------------------------------------------------------------------
  !--------------- 9) HOBr+HCl->H2O+BrCl ----------------
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RnatMeter = Radius of NAT particles / m
  ! NnatMeter = number density of particles in air / (1/m**3)
  ! HCl       = amount-of-substance ratio of gas phase HCl to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rnat     = Radius of NAT particles / cm
  ! Nnat     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHCl     = number concentration of gaseous HCl molecules / (1/cm**3)
  ! gamma9   = Uptake coefficient (unit 1)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_nat9 = khet9  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RnatMeter, NnatMeter
  REAL(dp), INTENT(in) :: HCl
  
  REAL(dp) :: khet9, mz_psc_het_nat9, gamma9
  REAL(dp) :: Rnat, Nnat
  REAL(dp) :: factor
  REAL(dp) :: cHCl

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rnat=RnatMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nnat=NnatMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHCl  = factor * HCl

  gamma9 = 0.1_dp !Carslaw original code (guess: analogy to HOBr+HBr)
  khet9  = gamma9*4.56e4_dp*SQRT(TEMP/96.91_dp)*Rnat**2*Nnat/ &
             (1.0_dp + 3.3e4_dp*gamma9*Rnat*PRESS/TEMP)/(cHCl+1.0_dp)
  khet9  = max(khet9, minKhet)
  khet9  = min(khet9, maxKhet)

  mz_psc_het_nat9 = khet9
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_het_nat9
!=========================================================================


!==========================================================================
 ELEMENTAL FUNCTION mz_psc_het_nat10(TEMP, PRESS,          &
                                     RnatMeter, NnatMeter, &
                                     HBr)
  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on nat particles.
  !-------------------------------------------------------------------------
  !--------------- 10) HOBr+HBr->H2O+Br2 ----------------
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RnatMeter = Radius of NAT particles / m
  ! NnatMeter = number density of particles in air / (1/m**3)
  ! HBr       = amount-of-substance ratio of gas phase HBr to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rnat     = Radius of NAT particles / cm
  ! Nnat     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHBr     = number concentration of gaseous HBr molecules / (1/cm**3)
  ! gamma10  = Uptake coefficient (unit 1)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_nat10 = khet10  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RnatMeter, NnatMeter
  REAL(dp), INTENT(in) :: HBr

  REAL(dp) :: khet10, mz_psc_het_nat10, gamma10
  REAL(dp) :: Rnat, Nnat
  REAL(dp) :: factor
  REAL(dp) :: cHBr

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rnat=RnatMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nnat=NnatMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHBr  = factor * HBr

  gamma10= 0.1_dp !Carslaw original code (guess: for ice gamma=0.1)
  khet10 = gamma10*4.56e4_dp*SQRT(TEMP/96.91_dp)*Rnat**2*Nnat/ &
           (1.0_dp + 3.3e4_dp*gamma10*Rnat*PRESS/TEMP)/(cHBr+1.0_dp)
  khet10 = max(khet10, minKhet)
  khet10 = min(khet10, maxKhet)

  mz_psc_het_nat10 = khet10
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_het_nat10
!=========================================================================



!==========================================================================
 ELEMENTAL FUNCTION mz_psc_het_nat11(TEMP, PRESS,          &
                                     RnatMeter, NnatMeter, &
                                     H2O)
  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on nat particles.
  !-------------------------------------------------------------------------
  !--------------- 11) BrONO2+H2O->HNO3+HOBr ------------
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RnatMeter = Radius of NAT particles / m
  ! NnatMeter = number density of particles in air / (1/m**3)
  ! H2O       = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rnat     = Radius of NAT particles / cm
  ! Nnat     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cH2O     = number concentration of gaseous H2O molecules / (1/cm**3)
  ! gamma11  = Uptake coefficient (unit 1)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_nat11 = khet11  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RnatMeter, NnatMeter
  REAL(dp), INTENT(in) :: H2O

  REAL(dp) :: khet11, mz_psc_het_nat11, gamma11
  REAL(dp) :: Rnat, Nnat
  REAL(dp) :: factor
  REAL(dp) :: cH2O

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rnat=RnatMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nnat=NnatMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cH2O  = factor * H2O

  gamma11= 0.001_dp !Carslaw original code
  khet11 = gamma11*4.56e4_dp*SQRT(TEMP/142.0_dp)*Rnat**2*Nnat/ &
           (1.0_dp + 3.3e4_dp*gamma11*Rnat*PRESS/TEMP)/(cH2O+1.0_dp)
  khet11 = max(khet11, minKhet)
  khet11 = min(khet11, maxKhet)

  mz_psc_het_nat11 = khet11
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_het_nat11
!=========================================================================


!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_ice12(TEMP, PRESS,             &
                                    RiceMeter, NiceMeter,    &
                                    HCl)

  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on ice particles.
  !-------------------------------------------------------
  !--------------- 12) HCl+ClONO2->HNO3+Cl2 --------------
  !-------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RiceMeter = Radius of ice particles / m
  ! NiceMeter = number density of particles in air / (1/m**3)
  ! HCl       = amount-of-substance ratio of gas phase HCl to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rice     = Radius of ice particles / cm
  ! Nice     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHCl     = number concentration of gaseous HCl molecules / (1/cm**3)
  ! gamma12  = uptake coefficient / 1
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_ice12 = khet12  = heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RiceMeter, NiceMeter
  REAL(dp), INTENT(in) :: HCl
  REAL(dp) :: mz_psc_het_ice12, khet12, gamma12

  REAL(dp) :: factor, Rice, Nice
  REAL(dp) :: cHCl

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rice=RiceMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nice=NiceMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHCl  = factor * HCl

  gamma12 = 0.3_dp !JPL1997/2000
  khet12  = gamma12*4.56e4_dp*SQRT(TEMP/97.46_dp)*Rice**2*Nice/ &
            (1.0_dp + 3.3e4_dp*gamma12*Rice*PRESS/TEMP)/         &
            (cHCl+1.0_dp)
  khet12 = max(khet12, minKhet)
  khet12 = min(khet12, maxKhet)

  mz_psc_het_ice12 = khet12
!------------------------------------------------------------------------
END FUNCTION mz_psc_het_ice12
!=========================================================================



!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_ice13(TEMP, PRESS,             &
                                    RiceMeter, NiceMeter,    &
                                    H2O)

  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on ice particles.
  !-------------------------------------------------------
  !--------------- 13) H2O+ClONO2->HNO3+HOCl -------------
  !-------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RiceMeter = Radius of ice particles / m
  ! NiceMeter = number density of particles in air / (1/m**3)
  ! H2O       = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rice     = Radius of ice particles / cm
  ! Nice     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cH2O     = number concentration of gaseous H2O molecules / (1/cm**3)
  ! gamma13  = uptake coefficient / 1
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_ice13 = khet13  = heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RiceMeter, NiceMeter
  REAL(dp), INTENT(in) :: H2O
  REAL(dp) :: mz_psc_het_ice13, khet13, gamma13

  REAL(dp) :: factor, Rice, Nice
  REAL(dp) :: cH2O

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rice=RiceMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nice=NiceMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cH2O  = factor * H2O

  gamma13 = 0.3_dp !JPL1997/2000
  khet13  = gamma13*4.56e4_dp*SQRT(TEMP/97.46_dp)*Rice**2*Nice/ &
              (1.0_dp + 3.3e4_dp*gamma13*Rice*PRESS/TEMP)/(cH2O+1.0_dp)
  khet13 = max(khet13, minKhet)
  khet13 = min(khet13, maxKhet)

  mz_psc_het_ice13 = khet13
!------------------------------------------------------------------------
END FUNCTION mz_psc_het_ice13
!=========================================================================



!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_ice14(TEMP, PRESS,            &
                                    RiceMeter, NiceMeter,   &
                                    HCl)

  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on ice particles.
  !-------------------------------------------------------
  !--------------- 14) HOCl+HCl->H2O+Cl2 -----------------
  !-------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RiceMeter = Radius of ice particles / m
  ! NiceMeter = number density of particles in air / (1/m**3)
  ! HCl       = amount-of-substance ratio of gas phase HCl to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rice     = Radius of ice particles / cm
  ! Nice     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHCl     = number concentration of gaseous HCl molecules / (1/cm**3)
  ! gamma14  = uptake coefficient / 1
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_ice14 = khet14  = heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RiceMeter, NiceMeter
  REAL(dp), INTENT(in) :: HCl
  REAL(dp) :: mz_psc_het_ice14, khet14, gamma14

  REAL(dp) :: factor, Rice, Nice
  REAL(dp) :: cHCl

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rice=RiceMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nice=NiceMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHCl  = factor * HCl

  gamma14 = 0.2_dp !JPL2000
  khet14  = gamma14*4.56e4_dp*SQRT(TEMP/52.46_dp)*Rice**2*Nice/ &
              (1.0_dp + 3.3e4_dp*gamma14*Rice*PRESS/TEMP)/         &
              (cHCl+1.0_dp)
  khet14 = max(khet14, minKhet)
  khet14 = min(khet14, maxKhet)

  mz_psc_het_ice14 = khet14
!-------------------------------------------------------------------------
END FUNCTION mz_psc_het_ice14
!=========================================================================



!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_ice15(TEMP, PRESS,            &
                                    RiceMeter, NiceMeter,   &
                                    HCl)

  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on ice particles.
  !-------------------------------------------------------
  !--------------- 15) N2O5+HCl->ClNO2+HNO3 --------------
  !-------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RiceMeter = Radius of ice particles / m
  ! NiceMeter = number density of particles in air / (1/m**3)
  ! HCl       = amount-of-substance ratio of gas phase HCl to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rice     = Radius of ice particles / cm
  ! Nice     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHCl     = number concentration of gaseous HCl molecules / (1/cm**3)
  ! gamma15  = uptake coefficient / 1
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_ice15 = khet15  = heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RiceMeter, NiceMeter
  REAL(dp), INTENT(in) :: HCl
  REAL(dp) :: mz_psc_het_ice15, khet15, gamma15

  REAL(dp) :: factor, Rice, Nice
  REAL(dp) :: cHCl

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rice=RiceMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nice=NiceMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHCl  = factor * HCl

  gamma15 = 0.03_dp  !JPL1997
  khet15  = gamma15*4.56e4_dp*SQRT(TEMP/108.0_dp)*Rice**2*Nice/ &
              (1.0_dp + 3.3e4_dp*gamma15*Rice*PRESS/TEMP)/(cHCl+1.0_dp)
  khet15 = max(khet15, minKhet)
  khet15 = min(khet15, maxKhet)

  mz_psc_het_ice15 = khet15
!-------------------------------------------------------------------------
END FUNCTION mz_psc_het_ice15
!=========================================================================



!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_ice16(TEMP, PRESS,            &
                                    RiceMeter, NiceMeter,   &
                                    H2O)

  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on ice particles.
  !-------------------------------------------------------
  !--------------- 16) N2O5+H2O->2HNO3 -------------------
  !-------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RiceMeter = Radius of ice particles / m
  ! NiceMeter = number density of particles in air / (1/m**3)
  ! H2O       = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rice     = Radius of ice particles / cm
  ! Nice     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cH2O     = number concentration of gaseous H2O molecules / (1/cm**3)
  ! gamma16  = uptake coefficient / 1
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_ice16 = khet16  = heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RiceMeter, NiceMeter
  REAL(dp), INTENT(in) :: H2O
  REAL(dp) :: mz_psc_het_ice16, khet16, gamma16

  REAL(dp) :: factor, Rice, Nice
  REAL(dp) :: cH2O

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rice=RiceMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nice=NiceMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cH2O  = factor * H2O

  gamma16 = 0.02_dp !JPL2000
  khet16  = gamma16*4.56e4_dp*SQRT(TEMP/108.0_dp)*Rice**2*Nice/ &
              (1.0_dp + 3.3e4_dp*gamma16*Rice*PRESS/TEMP)/(cH2O+1.0_dp)
  khet16 = max(khet16, minKhet)
  khet16 = min(khet16, maxKhet)

  mz_psc_het_ice16 = khet16
!-------------------------------------------------------------------------
END FUNCTION mz_psc_het_ice16
!=========================================================================



!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_ice17(TEMP, PRESS,            &
                                    RiceMeter, NiceMeter,   &
                                    HBr)

  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on ice particles.
  !------------------------------------------------------
  !--------------- 17) ClONO2+HBr->HNO3+BrCl -------------
  !------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RiceMeter = Radius of ice particles / m
  ! NiceMeter = number density of particles in air / (1/m**3)
  ! HBr       = amount-of-substance ratio of gas phase HBr to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rice     = Radius of ice particles / cm
  ! Nice     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHBr     = number concentration of gaseous HBr molecules / (1/cm**3)
  ! gamma17  = uptake coefficient / 1
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_ice17 = khet17  = heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RiceMeter, NiceMeter
  REAL(dp), INTENT(in) :: HBr
  REAL(dp) :: mz_psc_het_ice17, khet17, gamma17

  REAL(dp) :: factor, Rice, Nice
  REAL(dp) :: cHBr

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rice=RiceMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nice=NiceMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHBr  = factor * HBr

  gamma17 = 0.3_dp  !JPL1997 gamma>0.3_dp
  khet17  = gamma17*4.56e4_dp*SQRT(TEMP/97.46_dp)*Rice**2*Nice/ &
              (1.0_dp + 3.3e4_dp*gamma17*Rice*PRESS/TEMP)/          &
              (cHBr+1.0_dp)
  khet17 = max(khet17, minKhet)
  khet17 = min(khet17, maxKhet)

  mz_psc_het_ice17 = khet17
!-------------------------------------------------------------------------
END FUNCTION mz_psc_het_ice17
!=========================================================================



!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_ice18(TEMP, PRESS,             &
                                    RiceMeter, NiceMeter,    &
                                    HCl)

  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on ice particles.
  !-------------------------------------------------------
  !--------------- 18) BrONO2+HCl->HNO3+BrCl -------------
  !-------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RiceMeter = Radius of ice particles / m
  ! NiceMeter = number density of particles in air / (1/m**3)
  ! HCl       = amount-of-substance ratio of gas phase HCl to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rice     = Radius of ice particles / cm
  ! Nice     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHCl     = number concentration of gaseous HCl molecules / (1/cm**3)
  ! gamma18  = uptake coefficient / 1
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_ice18 = khet18  = heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RiceMeter, NiceMeter
  REAL(dp), INTENT(in) :: HCl
  REAL(dp) :: mz_psc_het_ice18, khet18, gamma18

  REAL(dp) :: factor, Rice, Nice
  REAL(dp) :: cHCl

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rice=RiceMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nice=NiceMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHCl  = factor * HCl

  gamma18 = 0.3_dp !Carslaw original code (guess)
  khet18  = gamma18*4.56e4_dp*SQRT(TEMP/142.0_dp)*Rice**2*Nice/ &
              (1.0_dp + 3.3e4_dp*gamma18*Rice*PRESS/TEMP)/        &
              (cHCl+1.0_dp)
  khet18 = max(khet18, minKhet)
  khet18 = min(khet18, maxKhet)

  mz_psc_het_ice18 = khet18
!-------------------------------------------------------------------------
END FUNCTION mz_psc_het_ice18
!=========================================================================



!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_ice19(TEMP, PRESS,            &
                                    RiceMeter, NiceMeter,   &
                                    HBr)

  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on ice particles.
  !------------------------------------------------------
  !--------------- 19) HBr+HOCl->H2O+BrCl ----------------
  !------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RiceMeter = Radius of ice particles / m
  ! NiceMeter = number density of particles in air / (1/m**3)
  ! HBr       = amount-of-substance ratio of gas phase HBr to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rice     = Radius of ice particles / cm
  ! Nice     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHBr     = number concentration of gaseous HBr molecules / (1/cm**3)
  ! gamma19  = uptake coefficient / 1
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_ice19 = khet19  = heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RiceMeter, NiceMeter
  REAL(dp), INTENT(in) :: HBr
  REAL(dp) :: mz_psc_het_ice19, khet19, gamma19

  REAL(dp) :: factor, Rice, Nice
  REAL(dp) :: cHBr

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rice=RiceMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nice=NiceMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHBr  = factor * HBr

  gamma19 = 0.3_dp !in analogy to HOBr+HCl
  khet19  = gamma19*4.56e4_dp*SQRT(TEMP/52.46_dp)*Rice**2*Nice/ &
              (1.0_dp + 3.3e4_dp*gamma19*Rice*PRESS/TEMP)/         &
              (cHBr+1.0_dp)
  khet19 = max(khet19, minKhet)
  khet19 = min(khet19, maxKhet)

  mz_psc_het_ice19 = khet19
!-------------------------------------------------------------------------
END FUNCTION mz_psc_het_ice19
!=========================================================================



!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_ice20(TEMP, PRESS,            &
                                    RiceMeter, NiceMeter,   &
                                    HCl)

  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on ice particles.
  !------------------------------------------------------
  !--------------- 20) HOBr+HCl->H2O+BrCl ----------------
  !------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RiceMeter = Radius of ice particles / m
  ! NiceMeter = number density of particles in air / (1/m**3)
  ! HCl       = amount-of-substance ratio of gas phase HCl to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rice     = Radius of ice particles / cm
  ! Nice     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHCl     = number concentration of gaseous HCl molecules / (1/cm**3)
  ! gamma20  = uptake coefficient / 1
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_ice20 = khet20  = heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RiceMeter, NiceMeter
  REAL(dp), INTENT(in) :: HCl
  REAL(dp) :: mz_psc_het_ice20, khet20, gamma20

  REAL(dp) :: factor, Rice, Nice
  REAL(dp) :: cHCl

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rice=RiceMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nice=NiceMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHCl  = factor * HCl

  gamma20 = 0.3_dp !JPL1997/2000
  khet20  = gamma20*4.56e4_dp*SQRT(TEMP/96.91_dp)*Rice**2*Nice/ &
              (1.0_dp + 3.3e4_dp*gamma20*Rice*PRESS/TEMP)/         &
              (cHCl+1.0_dp)
  khet20 = max(khet20, minKhet)
  khet20 = min(khet20, maxKhet)

  mz_psc_het_ice20 = khet20
!-------------------------------------------------------------------------
END FUNCTION mz_psc_het_ice20
!=========================================================================



!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_ice21(TEMP, PRESS,            &
                                    RiceMeter, NiceMeter,   &
                                    HBr)

  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on ice particles.
  !------------------------------------------------------
  !---------------21) HOBr+HBr->H2O+Br2 ----------------
  !------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RiceMeter = Radius of ice particles / m
  ! NiceMeter = number density of particles in air / (1/m**3)
  ! HBr       = amount-of-substance ratio of gas phase HBr to dry air 
  !             / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rice     = Radius of ice particles / cm
  ! Nice     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cHBr     = number concentration of gaseous HBr molecules / (1/cm**3)
  ! gamma21  = uptake coefficient / 1
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_ice21 = khet21  = heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RiceMeter, NiceMeter
  REAL(dp), INTENT(in) :: HBr
  REAL(dp) :: mz_psc_het_ice21, khet21, gamma21

  REAL(dp) :: factor, Rice, Nice
  REAL(dp) :: cHBr

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rice=RiceMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nice=NiceMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cHBr  = factor * HBr

  gamma21= 0.1_dp !JPL1997
  khet21 = gamma21*4.56e4_dp*SQRT(TEMP/96.91_dp)*Rice**2*Nice/ &
             (1.0_dp + 3.3e4_dp*gamma21*Rice*PRESS/TEMP)/         &
             (cHBr+1.0_dp)
  khet21 = max(khet21, minKhet)
  khet21 = min(khet21, maxKhet)

  mz_psc_het_ice21 = khet21
!-------------------------------------------------------------------------
END FUNCTION mz_psc_het_ice21
!=========================================================================



!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_ice22(TEMP, PRESS,           &
                                    RiceMeter, NiceMeter,  &
                                    H2O)

  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on ice particles.
  !------------------------------------------------------
  !--------------- 22) BrONO2+H2O->HNO3+HOBr ------------
  !------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RiceMeter = Radius of ice particles / m
  ! NiceMeter = number density of particles in air / (1/m**3)
  ! H2O       = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rice     = Radius of ice particles / cm
  ! Nice     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cH2O     = number concentration of gaseous H2O molecules / (1/cm**3)
  ! gamma22  = uptake coefficient / 1
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_ice22 = khet22  = heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RiceMeter, NiceMeter
  REAL(dp), INTENT(in) :: H2O
  REAL(dp) :: mz_psc_het_ice22, khet22, gamma22

  REAL(dp) :: factor, Rice, Nice
  REAL(dp) :: cH2O

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rice=RiceMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nice=NiceMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cH2O  = factor * H2O

  gamma22= 0.3_dp !JPL1997/2000 (gamma>0.3_dp)
  khet22 = gamma22*4.56e4_dp*SQRT(TEMP/142.0_dp)*Rice**2*Nice/ &
             (1.0_dp + 3.3e4_dp*gamma22*Rice*PRESS/TEMP)/(cH2O+1.0_dp)
  khet22 = max(khet22, minKhet)
  khet22 = min(khet22, maxKhet)

  mz_psc_het_ice22 = khet22
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_het_ice22
!=========================================================================



!=========================================================================
ELEMENTAL FUNCTION mz_psc_het_liq23(TEMP, PRESS,     &
                                    HCl, HOCl,       &
                                    bH2SO4, bHNO3,   &
                                    wHCl,            &
                                    dens,            &
                                    hHOCl,           &
                                    a_liq, rmean, dl)
  !-------------------------------------------------------------------------
  ! This function calculates a seond order reaction rate for a heterogeneous
  ! reaction on liquid stratospheric aerosols
  ! (H2O-HNO3-H2SO4-HCl-HBr-HOCl-HOBr).
  ! Example: The reduction of A due to the heterogeneous reaction
  !   A + B -> C + D
  ! is
  !   d(A)/dt = -khet*[A]*[B]
  ! where khet is a second order reactin rate measured in 1/(s*cm**3) and [A],
  ! [B] are particle number concentrations, i.e. molecules per cm**3.
  ! -----------------------------------------------------------------
  ! ---------------------- 23) HOCL + HCL ---------------------------
  ! --------------------- ON LIQUID AEROSOL -------------------------
  ! -----------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = temperature / K
  ! PRESS     = pressure / hPa
  ! HCl       = amount-of-substance ratio of gas phase HCl to dry air 
  !             / (mol/mol)
  ! HOCl      = amount-of-substance ratio of gas phase HOCl to dry air 
  !             / (mol/mol)
  ! bH2SO4    = molality of H2SO4 in water / (mol/kg)
  ! bHNO3     = molality of HNO3  in water / (mol/kg)
  ! wHCl      = mass fraction of HCl in liquid / (kg/kg)
  ! dens      = density of ternary solution / (g/cm**3)
  ! hHOCl     = henrys law constant of HOCl solved in water / (mol/(kg*atm))
  ! a_liq     = total area of liquid aerosols per volume of air / (cm**2/cm**3)
  ! rmean     = mean radius of liquid droplets / cm
  !           Note: This varies with the total liquid volume per unit volume of
  !           air (called vliq)
  ! dl        = dl liquid phase diffusion constant / (cm**2/s)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! gamma23   = uptake coefficient / 1
  ! cHCl_l    = amount-of-substance concentration of HCl in liquid / (mol/dm**3)
  ! factor    = conversion factor 
  !             (from amount-of-substance ratio / (mol/mol)
  !              to particle number concentration / (1/cm**3))
  ! fhHOCl    = Henry's law constant / (mol/(atm*cm**3))
  ! k1HOClHCl = first order reaction rate / (dm**3/(mol*s))
  ! zHCl      = particle number concentration of HCl / (1/cm**3)
  ! zHOCl     = particle number concentration of HOCl / (1/cm**3)
  ! cbar, fq, q = auxiliary variables
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_liq23 = khet23  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: HCl, HOCl
  REAL(dp), INTENT(in) :: bH2SO4, bHNO3
  REAL(dp), INTENT(in) :: wHCl
  REAL(dp), INTENT(in) :: dens
  REAL(dp), INTENT(in) :: hHOCl
  REAL(dp), INTENT(in) :: a_liq, rmean, dl

  REAL(dp) :: gamma23, khet23, mz_psc_het_liq23

  REAL(dp) :: cHCl_l
  REAL(dp) :: factor
  REAL(dp) :: fhHOCl
  REAL(dp) :: k1HOClHCl
  REAL(dp) :: zHCl, zHOCl
  REAL(dp) :: cbar, fq, q

  ! first order rate coeff. for HOCl+HCl: ! mz_rs_20060115: 2nd order?
  REAL(dp), PARAMETER :: k2HOClHCl = 1.0e5_dp

  INTRINSIC sqrt, max, min

  !-------------------------
  ! END VARIABLE DECLARATION
  !-------------------------

  ! -----
  ! amount-of-substance concentration in one droplet / (mol/dm**3)
  ! -----
  cHCl_l  = max(wHCl,0.0_dp)  * dens / (MolMassHCl/1000.0_dp)

  ! -----
  ! first order rate constants for reaction in liquid / (dm**3/(mol*s))
  ! -----
  k1HOClHCl = k2HOClHCl * cHCl_l

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor = mix2conc(PRESS, TEMP) ! conversion factor
  zHCl   = HCl  * factor
  zHOCl  = HOCl * factor

  ! -----
  ! convert henrys law constants hHX / (mol/(kg*atm))
  ! to other units fHX / (mol/(atm*cm**3))
  ! -----
  fhHOCl = hHOCl * dens / (1.0_dp                            &
                           + bH2SO4 * MolMassH2SO4/1000.0_dp &
                           + bHNO3                           &
                           + MolMassHNO3/1000.0_dp)

  IF ((fhHOCl > 0.0_dp)                  &
       .AND. (max(zHCl,zHOCl) > 0.0_dp)  &
      ! mz_rs_20060115: why max and not min? Is it okay if only one is >0?
       .AND. k1HOClHCl > 0.0_dp) THEN
    cbar = 2008.0_dp * sqrt(TEMP)
    q = qfunc(rmean,k1HOClHCl,dl)
    fq = fqfunc(q)

    gamma23 = gamfunc(fq,cbar,fhHOCl,TEMP,k1HOClHCl,dl)
    khet23= gamma23 * cbar/4.0_dp * a_liq  &
            / (zHCl+1.0_dp)
  ELSE
    khet23 = 0.0_dp
  END IF
  khet23 = max(khet23, minKhet)
  khet23 = min(khet23, maxKhet)

  mz_psc_het_liq23 = khet23
!---------------------------------------------------------------------------
END FUNCTION mz_psc_het_liq23
!===========================================================================



!=========================================================================
ELEMENTAL FUNCTION mz_psc_het_liq24(TEMP, PRESS, H2O,  &
                                    HCl,               &
                                    ClNO3, wHCl,       &
                                    dens,              &
                                    a_liq, rmean)
  !-------------------------------------------------------------------------
  ! This function calculates a seond order reaction rate for a heterogeneous
  ! reaction on liquid stratospheric aerosols
  ! (H2O-HNO3-H2SO4-HCl-HBr-HOCl-HOBr).
  ! Example: The reduction of A due to the heterogeneous reaction
  !   A + B -> C + D
  ! is
  !   d(A)/dt = -khet*[A]*[B]
  ! where khet is a second order reactin rate measured in 1/(s*cm**3) and [A],
  ! [B] are particle number concentrations, i.e. molecules per cm**3.
  ! -----------------------------------------------------------------
  ! ------------------24) ClONO2+HCl -> Cl2+HNO3 --------------------
  ! ----------------------- ON LIQUID AEROSOL -----------------------
  ! -----------------------------------------------------------------
  ! taken directly from Hanson and Ravishankara, J.Phys.Chem, 98,
  ! 5728, 1994, except for function f (see following comments) and
  ! for the HCl solubility, which is calculated according to Luo et al.
  ! -----------------------------------------------------------------
  ! Function f: The form of g used by Hanson and Ravishankara
  ! can "explode" under certain conditions. It has been replace here
  ! a stable function that is accurate within about 4%. This is also
  ! the case for other reactions that follow
  ! -----------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = temperature / K
  ! PRESS     = pressure / hPa
  ! H2O       = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! HCl       = amount-of-substance ratio of gas phase HCl to dry air 
  !             / (mol/mol)
  ! ClNO3     = amount-of-substance ratio of gas phase ClNO3 to dry air 
  !             / (mol/mol)
  ! wHCl      = mass fraction of HCl in liquid / (kg/kg)
  ! dens      = density of ternary solution / (g/cm**3)
  ! a_liq     = total area of liquid aerosols per volume of air / (cm**2/cm**3)
  ! rmean     = mean radius of liquid droplets / cm
  !           Note: This varies with the total liquid volume per unit volume of
  !           air (called vliq)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! gamma24   = uptake coefficient / 1
  ! cHCl_l    = amount-of-substance concentration of HCl in liquid / (mol/dm**3)
  ! factor    = conversion factor 
  !             (from amount-of-substance ratio / (mol/mol)
  !              to particle number concentration / (1/cm**3))
  ! zHCl      = particle number concentration of HCl / (1/cm**3)
  ! zClNO3    = particle number concentration of ClNO3 / (1/cm**3)
  ! aH2O, cbar, g0, gs, gcalc, ge, fq, q = auxiliary variables
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_liq24 = khet24  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS, H2O
  REAL(dp), INTENT(in) :: HCl
  REAL(dp), INTENT(in) :: ClNO3
  REAL(dp), INTENT(in) :: wHCl
  REAL(dp), INTENT(in) :: dens
  REAL(dp), INTENT(in) :: a_liq, rmean

  REAL(dp) :: gamma24, khet24, mz_psc_het_liq24

  REAL(dp) :: cHCl_l
  REAL(dp) :: factor
  REAL(dp) :: zHCl
  REAL(dp) :: zClNO3
  REAL(dp) :: aH2O, cbar, g0, gs, gcalc, ge, fq, p, q

  REAL(dp), PARAMETER :: ksur  = 576.0_dp  ! ? surface tension [??]
  REAL(dp), PARAMETER :: alpha = 0.30_dp   ! ? mass accomodation coeff. (unit 1)
  REAL(dp), PARAMETER :: rho   = 2.0e3_dp   ! approximative particle denisty

  INTRINSIC sqrt, max, min

  !-------------------------
  ! END VARIABLE DECLARATION
  !-------------------------

  !-----
  ! check for unrealistically low H2O
  !-----
  IF (H2O<1.0e-10_dp) THEN
    mz_psc_het_liq24 = minKhet
    return
  END IF

  ! -----
  ! amount-of-substance concentration in one droplet / (mol/dm**3)
  ! -----
  cHCl_l  = max(wHCl,0.0_dp)  * dens / (MolMassHCl/1000.0_dp)

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor = mix2conc(PRESS, TEMP) ! conversion factor
  zHCl   = HCl  * factor
  zClNO3 = ClNO3* factor

  IF (max(zHCl,zClNO3)>0.0_dp) THEN
    cbar = 1474.0_dp * sqrt(TEMP)
    aH2O=PRESS * H2O/(10**(9.217_dp-2190.0_dp/(TEMP-12.70_dp)))
    g0=1.18E-4_dp+9.1E-3_dp*aH2O+0.50_dp*aH2O**2
    gs=aH2O*ksur*cHCl_l
    p=RHO*cHCl_l/aH2O
    gcalc=g0*SQRT(1.0_dp+p)
    q = rmean * SQRT(aH2O)/ 1.4e-6_dp
    fq = fqfunc(q)
    ge=1.0_dp/(1.0_dp/(gs+fq*gcalc)+1.0_dp/ALPHA)
    gamma24 = ge * (gs+fq*gcalc*p/(1.0_dp+p)) / (gs+fq*gcalc)
    khet24  = gamma24 * cbar/4.0_dp * a_liq  &
             / (zHCl+1.0_dp)
  ELSE
    khet24 = 0.0_dp
  ENDIF
  khet24 = max(khet24, minKhet)
  khet24 = min(khet24, maxKhet)

   mz_psc_het_liq24 = khet24
!---------------------------------------------------------------------------
 END FUNCTION mz_psc_het_liq24
!===========================================================================




!=========================================================================
ELEMENTAL FUNCTION mz_psc_het_liq25(TEMP, PRESS, H2O,  &
                                    wHCl,              &
                                    dens,              &
                                    a_liq, rmean)
  !-------------------------------------------------------------------------
  ! This function calculates a seond order reaction rate for a heterogeneous
  ! reaction on liquid stratospheric aerosols
  ! (H2O-HNO3-H2SO4-HCl-HBr-HOCl-HOBr).
  ! Example: The reduction of A due to the heterogeneous reaction
  !   A + B -> C + D
  ! is
  !   d(A)/dt = -khet*[A]*[B]
  ! where khet is a second order reactin rate measured in 1/(s*cm**3) and [A],
  ! [B] are particle number concentrations, i.e. molecules per cm**3.
  ! -----------------------------------------------------------------
  ! ------------------25) ClONO2+H2O -> HOCl+HNO3 -------------------
  ! ----------------------- ON LIQUID AEROSOL -----------------------
  ! -----------------------------------------------------------------
  ! taken directly from Hanson and Ravishankara, J.Phys.Chem, 98,
  ! 5728, 1994, except for function f (see following comments) and
  ! for the HCl solubility, which is calculated according to Luo et al.
  ! -----------------------------------------------------------------
  ! Function f: The form of g used by Hanson and Ravishankara
  ! can "explode" under certain conditions. It has been replace here
  ! a stable function that is accurate within about 4%. This is also
  ! the case for other reactions that follow
  ! -----------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = temperature / K
  ! PRESS     = pressure / hPa
  ! H2O       = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! wHCl      = mass fraction of HCl in liquid / (kg/kg)
  ! dens      = density of ternary solution / (g/cm**3)
  ! a_liq     = total area of liquid aerosols per volume of air / (cm**2/cm**3)
  ! rmean     = mean radius of liquid droplets / cm
  !           Note: This varies with the total liquid volume per unit volume of
  !           air (called vliq)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! gamma25   = uptake coefficient / 1
  ! cHCl_l    = amount-of-substance concentration of HCl in liquid / (mol/dm**3)
  ! zH2O      = particle number concentration of H2O / (1/cm**3)
  ! aH2O, cbar, g0, gs, gcalc, ge, fq, p, q = auxiliary variables
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_liq25 = khet25  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS, H2O
  REAL(dp), INTENT(in) :: wHCl
  REAL(dp), INTENT(in) :: dens
  REAL(dp), INTENT(in) :: a_liq, rmean

  REAL(dp) :: gamma25, khet25, mz_psc_het_liq25

  REAL(dp) :: cHCl_l
  REAL(dp) :: zH2O
  REAL(dp) :: aH2O, cbar, g0, gs, gcalc, ge, fq, p, q

  REAL(dp), PARAMETER :: ksur  = 576.0_dp  ! ? surface tension [??]
  REAL(dp), PARAMETER :: alpha = 0.30_dp   ! ? mass accomodation coeff. (unit 1)
  REAL(dp), PARAMETER :: rho   = 2.0e3_dp   ! approximative particle denisty

  INTRINSIC sqrt, max, min

  !-------------------------
  ! END VARIABLE DECLARATION
  !-------------------------

  !-----
  ! check for unrealistically low H2O
  !-----
  IF (H2O<1.0e-10_dp) THEN
    mz_psc_het_liq25 = minKhet
    return
  END IF

  ! -----
  ! amount-of-substance concentration in one droplet / (mol/dm**3)
  ! -----
  cHCl_l  = max(wHCl,0.0_dp)  * dens / (MolMassHCl/1000.0_dp)

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  zH2O   = H2O  * mix2conc(PRESS, TEMP)

  cbar = 1474.0_dp * sqrt(TEMP)
  aH2O=PRESS*H2O/(10**(9.217_dp-2190.0_dp/(TEMP-12.70_dp)))
  g0=1.18E-4_dp+9.1E-3_dp*aH2O+0.50_dp*aH2O**2
  gs=aH2O*ksur*cHCl_l
  p=RHO*cHCl_l/aH2O
  gcalc=g0*SQRT(1.0_dp+p)
  q = rmean / (1.4e-6_dp * SQRT(1.0_dp/aH2O))
  fq = fqfunc(q)
  ge=1.0_dp/(1.0_dp/(gs+fq*gcalc)+1.0_dp/ALPHA)
  gamma25 = ge - ge * (gs+fq*gcalc*p/(1.0_dp+p)) / (gs+fq*gcalc)
  khet25  = gamma25 * cbar/4.0_dp * a_liq / (zH2O+1.0_dp)

  khet25 = max(khet25, minKhet)
  khet25 = min(khet25, maxKhet)

  mz_psc_het_liq25 = khet25
!---------------------------------------------------------------------------
 END FUNCTION mz_psc_het_liq25
!===========================================================================



!=========================================================================
ELEMENTAL FUNCTION mz_psc_het_liq26(TEMP, PRESS, H2O,  &
                                    a_liq)
  !-------------------------------------------------------------------------
  ! This function calculates a seond order reaction rate for a heterogeneous
  ! reaction on liquid stratospheric aerosols
  ! (H2O-HNO3-H2SO4-HCl-HBr-HOCl-HOBr).
  ! Example: The reduction of A due to the heterogeneous reaction
  !   A + B -> C + D
  ! is
  !   d(A)/dt = -khet*[A]*[B]
  ! where khet is a second order reactin rate measured in 1/(s*cm**3) and [A],
  ! [B] are particle number concentrations, i.e. molecules per cm**3.
  ! -----------------------------------------------------------------
  ! ------------------------26) N2O5+H2O -> 2HNO3 -------------------
  ! --------------------------- ON LIQUID AEROSOL -------------------
  ! -----------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = temperature / K
  ! PRESS     = pressure / hPa
  ! H2O       = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! a_liq     = total area of liquid aerosols per volume of air / (cm**2/cm**3)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! gamma26   = uptake coefficient / 1
  ! zH2O      = particle number concentration of H2O / (1/cm**3)
  ! cbar      = auxiliary variables
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_liq26 = khet26  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS, H2O
  REAL(dp), INTENT(in) :: a_liq

  REAL(dp) :: gamma26, khet26, mz_psc_het_liq26

  REAL(dp) :: zH2O
  REAL(dp) :: cbar

  INTRINSIC sqrt, max, min

  !-------------------------
  ! END VARIABLE DECLARATION
  !-------------------------

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  zH2O   = H2O  * mix2conc(PRESS, TEMP)

  cbar = 1400.1_dp * sqrt(TEMP)

  gamma26 = 0.1_dp
  khet26  = gamma26 * cbar/4.0_dp * a_liq / (zH2O+1.0_dp)

  khet26 = max(khet26, minKhet)
  khet26 = min(khet26, maxKhet)

  mz_psc_het_liq26 = khet26
!---------------------------------------------------------------------------
 END FUNCTION mz_psc_het_liq26
!===========================================================================




!=========================================================================
ELEMENTAL FUNCTION mz_psc_het_liq27(TEMP, PRESS,     &
                                    HCl, HOBr,       &
                                    bH2SO4, bHNO3,   &
                                    wHCl,            &
                                    dens,            &
                                    hHOBr,           &
                                    a_liq, rmean, dl)
  !-------------------------------------------------------------------------
  ! This function calculates a seond order reaction rate for a heterogeneous
  ! reaction on liquid stratospheric aerosols
  ! (H2O-HNO3-H2SO4-HCl-HBr-HOCl-HOBr).
  ! Example: The reduction of A due to the heterogeneous reaction
  !   A + B -> C + D
  ! is
  !   d(A)/dt = -khet*[A]*[B]
  ! where khet is a second order reactin rate measured in 1/(s*cm**3) and [A],
  ! [B] are particle number concentrations, i.e. molecules per cm**3.
  ! -----------------------------------------------------------------
  ! -------------------------27) HOBR + HCL -------------------------
  ! ------------------------ ON LIQUID AEROSOL ----------------------
  ! -----------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = temperature / K
  ! PRESS     = pressure / hPa
  ! HCl       = amount-of-substance ratio of gas phase HCl to dry air 
  !             / (mol/mol)
  ! HOBr      = amount-of-substance ratio of gas phase HOBr to dry air 
  !             / (mol/mol)
  ! bH2SO4    = molality of H2SO4 in water / (mol/kg)
  ! bHNO3     = molality of HNO3  in water / (mol/kg)
  ! wHCl      = mass fraction of HCl in liquid / (kg/kg)
  ! dens      = density of ternary solution / (g/cm**3)
  ! hHOBr     = henrys law constant of HOBr solved in water / (mol/(kg*atm))
  ! a_liq     = total area of liquid aerosols per volume of air / (cm**2/cm**3)
  ! rmean     = mean radius of liquid droplets / cm
  !           Note: This varies with the total liquid volume per unit volume of
  !           air (called vliq)
  ! dl        = dl liquid phase diffusion constant / (cm**2/s)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! gamma27   = uptake coefficient / 1
  ! cHCl_l    = amount-of-substance concentration of HCl in liquid / (mol/dm**3)
  ! factor    = conversion factor 
  !             (from amount-of-substance ratio / (mol/mol)
  !              to particle number concentration / (1/cm**3))
  ! fhHOBr    = Henry's law constant / (mol/(atm*cm**3))
  ! k1HOBrHCl = first order reaction rate / (dm**3/(mol*s))
  ! zHCl      = particle number concentration of HCl / (1/cm**3)
  ! zHOBr     = particle number concentration of HOBr / (1/cm**3)
  ! cbar, fq, q = auxiliary variables
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_liq27 = khet27  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: HCl, HOBr
  REAL(dp), INTENT(in) :: bH2SO4, bHNO3
  REAL(dp), INTENT(in) :: wHCl
  REAL(dp), INTENT(in) :: dens
  REAL(dp), INTENT(in) :: hHOBr
  REAL(dp), INTENT(in) :: a_liq, rmean, dl

  REAL(dp) :: gamma27, khet27, mz_psc_het_liq27

  REAL(dp) :: cHCl_l
  REAL(dp) :: factor
  REAL(dp) :: fhHOBr
  REAL(dp) :: k1HOBrHCl
  REAL(dp) :: zHCl, zHOBr
  REAL(dp) :: cbar, q, fq

  ! first order rate coeff. for HOBr+HCl:
  REAL(dp), PARAMETER :: k2HOBrHCl = 1.0e5_dp

  INTRINSIC sqrt, max, min

  !-------------------------
  ! END VARIABLE DECLARATION
  !-------------------------

  ! -----
  ! amount-of-substance concentration in one droplet / (mol/dm**3)
  ! -----
  cHCl_l  = max(wHCl,0.0_dp)  * dens / (MolMassHCl/1000.0_dp)

  ! -----
  ! first order rate constants for reaction in liquid / (dm**3/(mol*s))
  ! -----
  k1HOBrHCl = k2HOBrHCl * cHCl_l

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor = mix2conc(PRESS, TEMP) ! conversion factor
  zHCl   = HCl  * factor
  zHOBr  = HOBr * factor

  ! -----
  ! convert henrys law constants hHX / (mol/(kg*atm))
  ! to other units fHX / (mol/(atm*cm**3))
  ! -----
  fhHOBr = hHOBr * dens / (1.0_dp                    &
                           + bH2SO4 * MolMassH2SO4/1000.0_dp &
                           + bHNO3                           &
                           + MolMassHNO3/1000.0_dp)

  IF ((fhHOBr > 0.0_dp)                  &
      .AND. (max(zHCl,zHOBr) > 0.0_dp)   &
      .AND. k1HOBrHCl > 0.0_dp) THEN
    cbar = 1477.0_dp * sqrt(TEMP)
    q = qfunc(rmean,k1HOBrHCl,dl)
    fq = fqfunc(q)

    gamma27 = gamfunc(fq,cbar,fhHOBr,TEMP,k1HOBrHCl,dl)
    khet27  = gamma27 * cbar/4.0_dp * a_liq  &
              / (zHCl+1.0_dp)
  ELSE
    khet27 = 0.0_dp
  END IF

  khet27 = max(khet27, minKhet)
  khet27 = min(khet27, maxKhet)

  mz_psc_het_liq27 = khet27
!---------------------------------------------------------------------------
 END FUNCTION mz_psc_het_liq27
!===========================================================================




!=========================================================================
ELEMENTAL FUNCTION mz_psc_het_liq28(TEMP, PRESS,     &
                                    HBr, HOBr,       &
                                    partHBr,         &
                                    bH2SO4, bHNO3,   &
                                    wHOBr, dens,     &
                                    hHBr, hHOBr,     &
                                    a_liq, rmean, dl)
  !-------------------------------------------------------------------------
  ! This function calculates a seond order reaction rate for a heterogeneous
  ! reaction on liquid stratospheric aerosols
  ! (H2O-HNO3-H2SO4-HCl-HBr-HOCl-HOBr).
  ! Example: The reduction of A due to the heterogeneous reaction
  !   A + B -> C + D
  ! is
  !   d(A)/dt = -khet*[A]*[B]
  ! where khet is a second order reactin rate measured in 1/(s*cm**3) and [A],
  ! [B] are particle number concentrations, i.e. molecules per cm**3.
  ! -----------------------------------------------------------------
  ! ----------------------28) HBR + HOBr ----------------------------
  ! --------------------- ON LIQUID AEROSOL -------------------------
  ! -----------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = temperature / K
  ! PRESS     = pressure / hPa
  ! HBr       = amount-of-substance ratio of gas phase HBr to dry air 
  !             / (mol/mol)
  ! HOBr      = amount-of-substance ratio of gas phase HOBr to dry air 
  !             / (mol/mol)
  ! partHBr   = partitioning (fraction in the gas phase) of HBr
  ! bH2SO4    = molality of H2SO4 in water / (mol/kg)
  ! bHNO3     = molality of HNO3  in water / (mol/kg)
  ! wHOBr     = mass fraction of HOBr in liquid / (kg/kg)
  ! dens      = density of ternary solution / (g/cm**3)
  ! hHBr      = henrys law constant of HBr solved in water / (mol/(kg*atm))
  ! hHOBr     = henrys law constant of HOBr solved in water / (mol/(kg*atm))
  ! a_liq     = total area of liquid aerosols per volume of air / (cm**2/cm**3)
  ! rmean     = mean radius of liquid droplets / cm
  !           Note: This varies with the total liquid volume per unit volume of
  !           air (called vliq)
  ! dl        = dl liquid phase diffusion constant / (cm**2/s)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! gamma28   = uptake coefficient / 1
  ! cHOBr_l   = amount-of-substance concentration of HOBr in liquid 
  !             / (mol/dm**3)
  ! factor    = conversion factor 
  !             (from amount-of-substance ratio / (mol/mol)
  !              to particle number concentration / (1/cm**3))
  ! fhHBr     = Henry's law constant / (mol/(atm*cm**3))
  ! fhHOBr    = Henry's law constant / (mol/(atm*cm**3))
  ! k1HOBrHBr = first order reaction rate / (dm**3/(mol*s))
  ! zHBr      = particle number concentration of HBr / (1/cm**3)
  ! zHOBr     = particle number concentration of HOBr / (1/cm**3)
  ! cbar, fq, q = auxiliary variables
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_liq28 = khet28  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: HBr, HOBr
  REAL(dp), INTENT(in) :: partHBr
  REAL(dp), INTENT(in) :: bH2SO4, bHNO3
  REAL(dp), INTENT(in) :: wHOBr, dens
  REAL(dp), INTENT(in) :: hHBr, hHOBr
  REAL(dp), INTENT(in) :: a_liq, rmean, dl

  REAL(dp) :: gamma28, khet28, mz_psc_het_liq28

  REAL(dp) :: cHOBr_l
  REAL(dp) :: factor
  REAL(dp) :: fhHBr, fhHOBr
  REAL(dp) :: k1HOBrHBr
  REAL(dp) :: zHBr, zHOBr
  REAL(dp) :: cbar, fq, q

  ! first order rate coeff. for HOBr+HBr:
  REAL(dp), PARAMETER :: k2HOBrHBr = 1.0e7_dp

  INTRINSIC sqrt, max, min

  !-------------------------
  ! END VARIABLE DECLARATION
  !-------------------------

  ! -----
  ! amount-of-substance concentration in one droplet / (mol/dm**3)
  ! -----
  cHOBr_l = max(wHOBr,0.0_dp) * dens / (MolMassHOBr/1000.0_dp)

  ! -----
  ! first order rate constants for reaction in liquid / (dm**3/(mol*s))
  ! -----
  k1HOBrHBr = k2HOBrHBr * cHOBr_l

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor = mix2conc(PRESS, TEMP) ! conversion factor
  zHBr   = HBr  * factor
  zHOBr  = HOBr * factor

  ! -----
  ! convert henrys law constants hHX / (mol/(kg*atm))
  ! to other units fHX / (mol/(atm*cm**3))
  ! -----
  fhHBr = hHBr / (1.0_dp                            &
                  + bH2SO4 * MolMassH2SO4/1000.0_dp &
                  + bHNO3                           &
                  + MolMassHNO3/1000.0_dp)
  fhHOBr = hHOBr / (1.0_dp                            &
                    + bH2SO4 * MolMassH2SO4/1000.0_dp &
                    + bHNO3                           &
                    + MolMassHNO3/1000.0_dp)

  IF ((fhHOBr > 0.0_dp)                  &
      .AND. (max(zHBr,zHOBr) > 0.0_dp)   &
      .AND. k1HOBrHBr > 0.0_dp) THEN
    cbar = 1616.0_dp * sqrt(TEMP)
    q = qfunc(rmean,k1HOBrHBr,dl)
    fq = fqfunc(q)

    gamma28 = gamfunc(fq,cbar,fhHBr,TEMP,k1HOBrHBr,dl)
    khet28  = gamma28 * cbar/4.0_dp * a_liq* partHBr / (zHBr+1.0_dp)
  ELSE
    khet28 = 0.0_dp
  END IF

  khet28 = max(khet28, minKhet)
  khet28 = min(khet28, maxKhet)

  mz_psc_het_liq28 = khet28
!---------------------------------------------------------------------------
 END FUNCTION mz_psc_het_liq28
!===========================================================================




!=========================================================================
ELEMENTAL FUNCTION mz_psc_het_liq29(TEMP, PRESS,     &
                                    HBr, HOCl,       &
                                    partHBr,         &
                                    bH2SO4, bHNO3,   &
                                    wHOCl,           &
                                    dens,            &
                                    hHOCl, hHBr,     &
                                    a_liq, rmean, dl)
  !-------------------------------------------------------------------------
  ! This function calculates a seond order reaction rate for a heterogeneous
  ! reaction on liquid stratospheric aerosols
  ! (H2O-HNO3-H2SO4-HCl-HBr-HOCl-HOBr).
  ! Example: The reduction of A due to the heterogeneous reaction
  !   A + B -> C + D
  ! is
  !   d(A)/dt = -khet*[A]*[B]
  ! where khet is a second order reactin rate measured in 1/(s*cm**3) and [A],
  ! [B] are particle number concentrations, i.e. molecules per cm**3.
  ! -----------------------------------------------------------------
  ! ----------------------29) HBR + HOCL ----------------------------
  ! --------------------- ON LIQUID AEROSOL -------------------------
  ! -----------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = temperature / K
  ! PRESS     = pressure / hPa
  ! HBr       = amount-of-substance ratio of gas phase HBr to dry air 
  !             / (mol/mol)
  ! HOCl      = amount-of-substance ratio of gas phase HOCl to dry air 
  !             / (mol/mol)
  ! partHBr   = partitioning (fraction in the gas phase) of HBr
  ! bH2SO4    = molality of H2SO4 in water / (mol/kg)
  ! bHNO3     = molality of HNO3  in water / (mol/kg)
  ! wHOCl     = mass fraction of HOCl in liquid / (kg/kg)
  ! dens      = density of ternary solution / (g/cm**3)
  ! hHOCl     = henrys law constant of HOCl solved in water / (mol/(kg*atm))
  ! hHBr      = henrys law constant of HBr solved in water / (mol/(kg*atm))
  ! a_liq     = total area of liquid aerosols per volume of air / (cm**2/cm**3)
  ! rmean     = mean radius of liquid droplets / cm
  !           Note: This varies with the total liquid volume per unit volume of
  !           air (called vliq)
  ! dl        = dl liquid phase diffusion constant / (cm**2/s)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! gamma29   = uptake coefficient / 1
  ! cHOCl_l   = amount-of-substance concentration of HOCl in liquid 
  !             / (mol/dm**3)
  ! factor    = conversion factor 
  !             (from amount-of-substance ratio / (mol/mol)
  !              to particle number concentration / (1/cm**3))
  ! fhHOCl    = Henry's law constant / (mol/(atm*cm**3))
  ! fhHBr     = Henry's law constant / (mol/(atm*cm**3))
  ! k1HOClHBr = first order reaction rate / (dm**3/(mol*s))
  ! zHBr      = particle number concentration of HBr / (1/cm**3)
  ! zHOCl     = particle number concentration of HOCl / (1/cm**3)
  ! cbar, fq, q = auxiliary variables
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_liq29 = khet29  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: HBr, HOCl
  REAL(dp), INTENT(in) :: partHBr
  REAL(dp), INTENT(in) :: bH2SO4, bHNO3
  REAL(dp), INTENT(in) :: wHOCl
  REAL(dp), INTENT(in) :: dens
  REAL(dp), INTENT(in) :: hHOCl, hHBr
  REAL(dp), INTENT(in) :: a_liq, rmean, dl

  REAL(dp) :: gamma29, khet29, mz_psc_het_liq29

  REAL(dp) :: cHOCl_l
  REAL(dp) :: factor
  REAL(dp) :: fhHOCl, fhHBr
  REAL(dp) :: k1HOClHBr
  REAL(dp) :: zHBr, zHOCl
  REAL(dp) :: cbar, fq, q

  ! first order rate coeff. for HOCl+HBr:
  REAL(dp), PARAMETER :: k2HOClHBr = 1.0e6_dp

  INTRINSIC sqrt, max, min

  !-------------------------
  ! END VARIABLE DECLARATION
  !-------------------------

  ! -----
  ! amount-of-substance concentration in one droplet / (mol/dm**3)
  ! -----
  cHOCl_l = max(wHOCl,0.0_dp) * dens / (MolMassHOCl/1000.0_dp)

  ! -----
  ! first order rate constants for reaction in liquid / (dm**3/(mol*s))
  ! -----
  k1HOClHBr = k2HOClHBr * cHOCl_l

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor = mix2conc(PRESS, TEMP) ! conversion factor
  zHOCl  = HOCl * factor
  zHBr   = HBr  * factor

  ! -----
  ! convert henrys law constants hHX / (mol/(kg*atm))
  ! to other units fHX / (mol/(atm*cm**3))
  ! -----
  factor = dens / (1.0_dp                            &
                   + bH2SO4 * MolMassH2SO4/1000.0_dp &
                   + bHNO3                           &
                   + MolMassHNO3/1000.0_dp)
  fhHOCl = factor * hHOCl
  fhHBr  = factor * hHBr

  IF ((fhHOCl > 0.0_dp)                  &
      .AND. (max(zHBr,zHOCl) > 0.0_dp)   &
      .AND. k1HOClHBr > 0.0_dp) THEN
    cbar = 1616.0_dp * sqrt(TEMP)
    q = qfunc(rmean,k1HOClHBr,dl)
    fq = fqfunc(q)

    gamma29 = gamfunc(fq,cbar,fhHBr,TEMP,k1HOClHBr,dl)
    khet29  = gamma29 * cbar/4.0_dp * a_liq * partHBr / (zHBr+1.0_dp)
  ELSE
    khet29 = 0.0_dp
  END IF

  khet29 = max(khet29, minKhet)
  khet29 = min(khet29, maxKhet)

  mz_psc_het_liq29 = khet29
!---------------------------------------------------------------------------
 END FUNCTION mz_psc_het_liq29
!===========================================================================




!=========================================================================
ELEMENTAL FUNCTION mz_psc_het_liq30(TEMP, PRESS, H2O,  &
                                    a_liq)
  !-------------------------------------------------------------------------
  ! This function calculates a seond order reaction rate for a heterogeneous
  ! reaction on liquid stratospheric aerosols
  ! (H2O-HNO3-H2SO4-HCl-HBr-HOCl-HOBr).
  ! Example: The reduction of A due to the heterogeneous reaction
  !   A + B -> C + D
  ! is
  !   d(A)/dt = -khet*[A]*[B]
  ! where khet is a second order reactin rate measured in 1/(s*cm**3) and [A],
  ! [B] are particle number concentrations, i.e. molecules per cm**3.
  ! -----------------------------------------------------------------
  ! ----------------------30) BRONO2 + H2O --------------------------
  ! --------------------- ON LIQUID AEROSOL -------------------------
  ! -----------------------------------------------------------------
  ! FROM HANSON ET AL., JGR, 101, 9063-9069, 1996.
  ! -----------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = temperature / K
  ! PRESS     = pressure / hPa
  ! H2O       = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! a_liq     = total area of liquid aerosols per volume of air / (cm**2/cm**3)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! gamma30   = uptake coefficient / 1
  ! zH2O      = particle number concentration of H2O / (1/cm**3)
  ! aH2O, cbar, grxn = auxiliary variables
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_liq30 = khet30  = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS, H2O
  REAL(dp), INTENT(in) :: a_liq

  REAL(dp) :: gamma30, khet30, mz_psc_het_liq30

  REAL(dp) :: zH2O
  REAL(dp) :: aH2O
  REAL(dp) :: cbar
  REAL(dp) :: grxn

  INTRINSIC sqrt, max, min

  !-------------------------
  ! END VARIABLE DECLARATION
  !-------------------------

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  zH2O   = H2O  * mix2conc(PRESS, TEMP)

  cbar = 1221.4_dp * sqrt(TEMP)
  aH2O = PRESS*H2O/(10**(9.217_dp-2190.0_dp/(TEMP-12.70_dp)))
  grxn = 211.0_dp * aH2O**1.37_dp
  gamma30 = (0.84_dp * grxn) / (grxn + 0.84_dp)
  khet30  = gamma30 * cbar/4. * a_liq / (zH2O+1.0_dp)

  khet30 = max(khet30, minKhet)
  khet30 = min(khet30, maxKhet)

  mz_psc_het_liq30 = khet30
!---------------------------------------------------------------------------
 END FUNCTION mz_psc_het_liq30
!===========================================================================

!=========================================================================
ELEMENTAL FUNCTION mz_psc_het_liq_gen(TEMP, PRESS, H2O, M, gamma,  &
                                    a_liq)
  !-------------------------------------------------------------------------
  ! This function calculates a seond order reaction rate for a heterogeneous
  ! reaction on liquid stratospheric aerosols
  ! (H2O-HNO3-H2SO4-HCl-HBr-HOCl-HOBr).
  ! Example: The reduction of A due to the heterogeneous reaction
  !   A + B -> C + D
  ! is
  !   d(A)/dt = -khet*[A]*[B]
  ! where khet is a second order reactin rate measured in 1/(s*cm**3) and [A],
  ! [B] are particle number concentrations, i.e. molecules per cm**3.
  ! -----------------------------------------------------------------
  ! ------------------------gen) *_gas + H2O -> *_aq ----------------
  ! --------------------------- ON LIQUID AEROSOL -------------------
  ! -----------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = temperature / K
  ! PRESS     = pressure / hPa
  ! H2O       = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! M         = molar mass of species (kg/mol)
  ! gamma     = uptake coefficient (1)
  ! a_liq     = total area of liquid aerosols per volume of air / (cm**2/cm**3)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! zH2O      = particle number concentration of H2O / (1/cm**3)
  ! cbar      = Maxwell-Boltzmann mean velocity (cm/s)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_liq_gen = khet = Heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS, H2O, M, gamma
  REAL(dp), INTENT(in) :: a_liq

  REAL(dp) :: khet, mz_psc_het_liq_gen

  REAL(dp) :: zH2O
  REAL(dp) :: cbar

  INTRINSIC sqrt, max, min

  !-------------------------
  ! END VARIABLE DECLARATION
  !-------------------------

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  zH2O   = H2O  * mix2conc(PRESS, TEMP)

  ! vmean=sqrt(8*R_gas*T/(M*pi))      (M in kg/mol)
  ! sqrt(8*R_gas/pi)=4.60138
  ! * 100 (m -> cm)
  cbar = SQRT(TEMP/M)*460.138

  khet    = gamma * cbar/4.0_dp * a_liq / (zH2O+1.0_dp)

  khet = max(khet, minKhet)
  khet = min(khet, maxKhet)

  mz_psc_het_liq_gen = khet
!---------------------------------------------------------------------------
END FUNCTION mz_psc_het_liq_gen
!===========================================================================

!==========================================================================
ELEMENTAL FUNCTION mz_psc_het_sol_gen(TEMP, PRESS,          &
                                    RiceMeter, NiceMeter,   &
                                    H2O, M, gamma)

  !-------------------------------------------------------------------------
  ! This function calculates a second order reaction rate for a heterogeneous
  ! reaction on ice particles.
  !-------------------------------------------------------
  !--------------- gen) *_gas+ICE/NAT->*_ice -------------
  !-------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! TEMP      = Temperature / K
  ! PRESS     = Pressure / hPa
  ! RiceMeter = Radius of ice particles / m
  ! NiceMeter = number density of particles in air / (1/m**3)
  ! H2O       = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! M         = molar mass of species (kg/mol)
  ! gamma     = uptake coefficient (1)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! Rice     = Radius of ice particles / cm
  ! Nice     = number density of particles in air / (1/cm**3)
  ! factor   = conversion factor 
  !            (from amount-of-substance ratio / (mol/mol)
  !             to particle number concentration / (1/cm**3))
  ! cH2O     = number concentration of gaseous H2O molecules / (1/cm**3)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_het_ice = khet  = heterogeneous reaction rate / (cm**3/s)
  !   where number concentrations in 1/cm**3 of reactants are given
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TEMP, PRESS
  REAL(dp), INTENT(in) :: RiceMeter, NiceMeter
  REAL(dp), INTENT(in) :: H2O, gamma, M
  REAL(dp) :: mz_psc_het_sol_gen, khet 

  REAL(dp) :: factor, Rice, Nice
  REAL(dp) :: cH2O

  INTRINSIC sqrt, max, min

  !-----
  ! convert radius unit from m to cm
  !-----
  Rice=RiceMeter*100.0_dp

  !-----
  ! convert number density unit from 1/m**3 to 1/cm**3
  !-----
  Nice=NiceMeter*1.0E-6_dp

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/cm**3)
  ! -----
  factor  = mix2conc(PRESS, TEMP)
  cH2O  = factor * H2O

  ! Turco et al,(1989) as cited by Buchholz (PhD thesis, 2005)
  khet  = gamma*4.56e4_dp*SQRT(TEMP/(M*1.e3_dp))*Rice**2*Nice/ &
              (1.0_dp + 3.3e4_dp*gamma*Rice*PRESS/TEMP)/(cH2O+1.0_dp)
  khet = max(khet, minKhet)
  khet = min(khet, maxKhet)

  mz_psc_het_sol_gen = khet
!-------------------------------------------------------------------------
END FUNCTION mz_psc_het_sol_gen
!=========================================================================

!===========================================================================
ELEMENTAL FUNCTION qfunc(rm, kx, dl)
  !-----
  ! Auxiliary function used to calculate the uptake coefficients and
  ! heterogeneous reaction rates on liquid stratospheric aerosols, i.e.
  ! by the functions mz_psc_het_liq23, ...mz_psc_het_liq30.
  !-----
  IMPLICIT NONE

  REAL(dp), INTENT(in) :: rm, dl, kx
  REAL(dp) :: qfunc

  INTRINSIC sqrt

  qfunc = rm * sqrt(kx/dl)
END FUNCTION qfunc
!===========================================================================



!===========================================================================
ELEMENTAL FUNCTION fqfunc(qx)
  !-----
  ! Auxiliary function used to calculate the uptake coefficients and
  ! heterogeneous reaction rates on liquid stratospheric aerosols, i.e.
  ! by the functions mz_psc_het_liq23, ...mz_psc_het_liq30.
  !-----
  IMPLICIT NONE

  REAL(dp), INTENT(in) :: qx
  REAL(dp) :: fqfunc
  fqfunc = (qx + 0.312_dp*qx**2) / (3.0_dp + qx + 0.312_dp*qx**2)
END FUNCTION fqfunc
!===========================================================================



!===========================================================================
ELEMENTAL FUNCTION gamfunc(a1, a2, a3, a4, a5, a6)
  !-----
  ! Auxiliary function used to calculate the uptake coefficients and
  ! heterogeneous reaction rates on liquid stratospheric aerosols, i.e.
  ! by the functions mz_psc_het_liq23, ...mz_psc_het_liq30.
  !-----
  IMPLICIT NONE

  REAL(dp), INTENT(in) :: a1, a2, a3, a4, a5, a6
  REAL(dp) :: gamfunc

  INTRINSIC SQRT

  gamfunc = &
       a1 / (a1 + a2 /(4.0_dp * a3 * 0.082_dp * a4 * sqrt(a5*a6)))
END FUNCTION gamfunc
!===========================================================================



!=========================================================================

  ELEMENTAL FUNCTION T_nat(pH2O_tot, pHNO3_tot)
    ! -----
    ! calculates equilibrium temperature of NAT in K according to
    ! Hanson and Mauersberger, GRL, 15, 1988
    ! valid for 180-200 K
    ! INPUT :
    ! pH2O_tot = H2O pressure if all H2O molecules were gaseous / atm
    ! pHNO3_tot = HNO3 pressure of all HNO3 molecules were gaseous / atm
    ! Note: 760 Torr = 1 atm)
    ! OUTPUT:
    ! T_nat / K
    ! -----
    REAL(dp), INTENT(in) :: pH2O_tot, pHNO3_tot
    REAL(dp) :: rpH2O_tot, rpHNO3_tot ! mz_rs_20060125
    REAL(dp) :: aa, bb, cc
    REAL(dp), DIMENSION(5) :: kNAT
    REAL(dp) :: T_nat

    INTRINSIC LOG10, TINY

    kNAT(1)=2.7836_dp
    kNAT(2)=0.00088_dp
    kNAT(3)=38.9855_dp
    kNAT(4)=11397.0_dp
    kNAT(5)=0.009179_dp

    ! mz_rs_20060110+
!qqq+ PSC: 
!    aa=kNAT(2)*log10(pH2O_tot*760.0_dp)-kNAT(5)
!    bb=log10(pHNO3_tot*760.0_dp)+kNAT(1)*log10(pH2O_tot*760.0_dp)-kNAT(3)
!qqq-
    ! set pH2O_tot and pHNO3_tot to at least TINY to avoid problems
    ! with log10:
    rpH2O_tot = MAX(pH2O_tot, TINY(1.0_dp))
    rpHNO3_tot = MAX(pHNO3_tot, TINY(1.0_dp))
    aa=kNAT(2)*log10(rpH2O_tot*760.0_dp)-kNAT(5)
    bb=log10(rpHNO3_tot*760.0_dp)+kNAT(1)*log10(rpH2O_tot*760.0_dp)-kNAT(3)
    ! mz_rs_20060110-
    cc=kNAT(4)

    T_nat=(-bb-(bb**2-4.0_dp*aa*cc)**0.5_dp )/(2.0_dp*aa)

  END FUNCTION T_nat

!=========================================================================

  ELEMENTAL FUNCTION p_nat(TEMP, pH2O)
    ! -----
    ! Calculates equilibrium HNO3 vapour pressure over NAT in atm acc. to
    ! Hanson et al., GRL, 1988
    ! valid for 180-200 K
    ! -------
    ! INPUT :
    ! -------
    ! TEMP = temperature / K
    ! pH2O = H2O partial pressure / atm
    ! -------
    ! OUTPUT:
    ! -------
    ! p_nat = HNO3 vapour pressure over NAT particles / atm
    ! Notes: 760 Torr = 1 atm
    !        1 mbar = 1 hPa = atm/1013.25_dp
    ! -----
    REAL(dp), INTENT(in) :: TEMP, pH2O
    REAL(dp) :: aa, bb, cc, rpH2O
    REAL(dp), DIMENSION(5) :: kNAT   ! empirical constants according to
                                 ! Hanson et al. (1988)
    REAL(dp) :: p_nat

    INTRINSIC epsilon, log10, max

    rpH2O = MAX(pH2O, epsilon(1.0_dp))

    kNAT(1)=2.7836_dp
    kNAT(2)=0.00088_dp
    kNAT(3)=38.9855_dp
    kNAT(4)=11397.0_dp
    kNAT(5)=0.009179_dp

    aa=kNAT(1)+kNAT(2)*TEMP
    bb=kNAT(3)-kNAT(4)/TEMP+kNAT(5)*TEMP

    cc=-aa*log10(rpH2O*760.0_dp)+bb
    p_nat=(10.0_dp**cc)/760.0_dp

  END FUNCTION p_nat

!=========================================================================

  ELEMENTAL FUNCTION T_ice(pH2O_tot)
    !--------------------------------------------------------------------
    ! calculates the equilibrium temperature of ice ccording to Marti and 
    ! Mauersberger, GRL, 1993
    ! ------
    ! INPUT:
    ! ------
    ! pH2O_tot = H2O pressure if all H2O molecules were gaseous / hPa
    ! -------
    ! OUTPUT:
    ! -------
    ! T_ice = equilibrium temperature of ice / K
    !--------------------------------------------------------------------
    IMPLICIT NONE

    REAL(dp), INTENT(in) :: pH2O_tot
    REAL(dp) :: T_ice

    INTRINSIC epsilon, log10, max

    T_ice = 2663.5_dp/(12.537_dp &
                      -log10(max(pH2O_tot*100.0_dp, epsilon(1.0_dp))))

  END FUNCTION T_ice

!=========================================================================

  ELEMENTAL FUNCTION p_ice(TEMP)
    !---------------------------------------------------------------------
    ! calculation of equilibrium water vapour pressure over ice according 
    ! to Marti and Mauersberger, GRL, 1993
    ! ------
    ! INPUT: 
    ! ------
    ! TEMP = temperature / K
    !
    ! -------
    ! OUTPUT:
    ! -------
    ! p_ice = water vapour pressure over ice / atm
    !---------------------------------------------------------------------

    IMPLICIT NONE

    REAL(dp), INTENT(in) :: TEMP
    REAL(dp) :: p_ice

    p_ice= 10.0_dp**(-2663.5_dp/TEMP + 12.537_dp)/101325.0_dp

  END FUNCTION p_ice

!=========================================================================

 ELEMENTAL FUNCTION mix2conc(PRESS, TEMP)
   ! -----
   ! this function converts for given temperature in K and pressure in hPa
   ! amount-of-substance ratios of particles to dry air / (mol/mol) 
   ! to 
   ! particle number concentrations / (1/cm**3)
   ! -----
   IMPLICIT NONE

   REAL(dp), INTENT(in) :: PRESS, TEMP
   REAL(dp) :: mix2conc

   mix2conc = 7.2427E18_dp * PRESS/TEMP;
 END FUNCTION mix2conc

!=========================================================================
 ELEMENTAL FUNCTION mz_psc_diff(TEMP, dens, bHNO3, bH2SO4)
   ! -----
   ! Function to calculate the liquid phase diffusion constant (in cm**2/s)
   ! "diff" = diffusion constant of HOCl in H2SO4/HNO3 solution
   !          this is used also for HBr and HOBr.
   ! -----
   ! INPUT:
   ! mass fraction (wH2SO4,wHNO3) of H2SO4 and HNO3
   ! molality (bH2SO4,bHNO3) of H2SO4 and HNO3
   ! TEMP = temperature / K
   ! dens = density of the ternary solution / (g/cm**3)
   ! OUTPUT:
   ! diffusion coeff. (diff) for HOCl based on Houghton cubic cell model
   ! with cell dimension 3.65 Angstroems
   ! Note: this is in good agreement (+-10% with a composition dependent
   ! cell dimension as given in Huthwelker et al. (J.At.Sci., 1995)
   ! -----
   IMPLICIT NONE

   REAL(dp), INTENT(in) :: TEMP, dens, bH2SO4, bHNO3
   REAL(dp) :: wH2SO4, wnen
   REAL(dp) :: c    !wt% of H2SO4
   REAL(dp) :: visc !viscosity in SI units (mole ration of  H2SO4 and HNO3)
   REAL(dp) :: viss ,visn !viscosity of H2SO4 (s), HNO3 (n) in water
   REAL(dp) :: a, xb, t0, xn
   REAL(dp) :: mz_psc_diff   ! diffusion coefficient of HOCl / (cm**2/s)

   INTRINSIC EXP

   wnen = 1.0_dp + bH2SO4*MolMassH2SO4/1000.0_dp &
         + bHNO3*MolMassHNO3/1000.0_dp
   wH2SO4   = bH2SO4*MolMassH2SO4/1000.0_dp/ wnen

   c = wH2SO4*100.0_dp

   a = exp(-7.722133_dp + C*3.773159e-2_dp)
   xb= 623.8082_dp + 5.221606_dp*C - 8.085769e-2_dp*C**2 &
      + 2.1769575e-4_dp*C**3
   t0= 154.3466_dp - 0.9521694_dp*C - 2.6749929e-3_dp*C**2 &
      + 1.984055e-4_dp*C**3
   xn=0.5186040_dp

   viss=1.0e-3*a*TEMP**xn*EXP(xb/(TEMP-t0))

   a =0.02656_dp+0.001971_dp*bHNO3+0.0002376_dp*bHNO3**2
     !bHNO3 = molality of HNO3
   xb=735.7_dp
   t0=92.89_dp+0.6848_dp*bHNO3
   xn=-0.01275_dp*bHNO3

   visn=1.0e-3_dp*a*TEMP**xn*EXP(xb/(TEMP-t0))

   visc=bH2SO4/(bH2SO4+bHNO3)*viss + bHNO3/(bH2SO4+bHNO3)*visn

   mz_psc_diff=3.37e-14_dp * TEMP *1.0e3_dp * dens/visc
 END FUNCTION mz_psc_diff

!=========================================================================

  ELEMENTAL FUNCTION mz_psc_density(wH2SO4, wHNO3, TEMP)
     ! -----
     ! INPUT:
     ! wH2SO4 and wHNO3 are mass fraction of H2SO4 and HNO3
     ! TEMP = temperature / K
     ! OUTPUT:
     ! density of ternary solution in g/cm3
     ! taken from Luo et al., GRL, 23,3707-, 1996.
     ! (fitted to 0.05<wH2SO4+wHNO3<0.70, but extrapolates well 185 < TEMP
     ! -----
     IMPLICIT NONE

     REAL(dp), INTENT(in) :: wH2SO4, wHNO3, TEMP
     REAL(dp) :: w, wh, v1, vs, vn, vmcal
     REAL(dp), DIMENSION(22) :: X(22)
     REAL(dp) :: mz_psc_density

       !-----
       ! Note: Initialising the array X as a whole, like X = (/.../), would
       ! be more simple but caused difficulties on a NEC SX6 vector computer.
       !-----
       X(1) = 2.393284e-02_dp
       X(2) = -4.359335e-05_dp
       X(3) = 7.961181e-08_dp
       X(4) = 0.0_dp
       X(5) = -0.198716351_dp
       X(6) = 1.39564574e-03_dp
       X(7) = -2.020633e-06_dp
       X(8) = 0.51684706_dp
       X(9) = -3.0539e-03_dp
       X(10)= 4.505475e-06_dp
       X(11)= -0.30119511_dp
       X(12)= 1.840408e-03_dp
       X(13)= -2.7221253742e-06_dp
       X(14)= -0.11331674116_dp
       X(15)= 8.47763e-04_dp
       X(16)= -1.22336185e-06_dp
       X(17)= 0.3455282_dp
       X(18)= -2.2111e-03_dp
       X(19)= 3.503768245e-06_dp
       X(20)= -0.2315332_dp
       X(21)= 1.60074e-03_dp
       X(22)= -2.5827835e-06_dp

      w=wH2SO4+wHNO3
      wh=1.0_dp-w
      v1=X(1)+X(2)*TEMP+X(3)*TEMP**2+X(4)*TEMP**3
      vs=X(5)+X(6)*TEMP+X(7)*TEMP**2+(X(8)+X(9)*TEMP+X(10)*TEMP**2)*w       &
           +(X(11)+X(12)*TEMP+X(13)*TEMP**2)*w*w
      vn=X(14)+X(15)*TEMP+X(16)*TEMP**2+(X(17)+X(18)*TEMP+X(19)*TEMP**2)*w  &
           +(X(20)+X(21)*TEMP+X(22)*TEMP**2)*w*w
      vmcal=wh/18.0160_dp*v1 + vs*wH2SO4/98.080_dp + vn*wHNO3/63.0160_dp

      mz_psc_density=0.001_dp/vmcal

   END FUNCTION mz_psc_density

!=========================================================================

ELEMENTAL FUNCTION mz_psc_vel(PRESS, TEMP, radius)
  !-------------
  ! DESCRIPTION:
  !-------------
  ! The function calculates the velocity of fall of solid particles according 
  ! to Waibel (1997), pp. 104/105. 
  ! The parameters in this function are (empirical) constants necessary 
  ! for the calculation of the velocity; v_s has no physical meaning, it is 
  ! just a part of the velocity formula.
  !
  !-----------------
  ! INPUT VARIABLES:
  !-----------------
  ! PRESS    = air pressure within grid box / Pa
  ! TEMP     = air temperature / K
  ! radius   = aerosol particle radius / m
  !
  !--------
  ! OUTPUT:
  !--------
  ! mz_psc_vel      = velocity of fall of aerosol particles / (m/s)
  !------------------------------------------------------------------------
  IMPLICIT NONE

  REAL(dp), INTENT(in) :: PRESS, TEMP, radius
  REAL(dp) :: mz_psc_vel

  REAL(dp), PARAMETER :: &
    etaFactor=6.45E-8_dp,             & ! [etaFactor]=1 kg/(m*s*K)
    f1=1.12_dp,                       & ! [f1]=1
    f2=0.58248_dp,                    & ! [f1]=1
    alpha1=f2*1.246_dp*0.23E-4_dp/f1, & ! [alpha1]=1 m*Pa/K
    alpha2=f2*0.42_dp*0.23E-4_dp/f1,  & ! [alpha2]=1 m*Pa/K
    alpha3=0.23E-4_dp/0.87_dp           ! [alpha3]=1 m*Pa/K
  REAL(dp) :: v_s                       ! [v_s]=1 m/s

  INTRINSIC EXP

  IF (radius>=r_min) THEN
    v_s=2.0_dp*g_acc*dens_ice*radius**2/(9.0_dp*etaFactor*Temp)
    mz_psc_vel=(v_s/f1)*( 1.0_dp                                       &
                         +alpha1*Temp/(press*radius)                   &
                         +alpha2*Temp*exp(-press*radius/(alpha3*Temp)) &
                          /(press*radius))
  ELSE
    mz_psc_vel=0.0_dp
  END IF
  
END FUNCTION mz_psc_vel

!=========================================================================

ELEMENTAL FUNCTION mz_psc_SedStep(TimeStep, PRESS, TEMP, vel)
  !-------------
  ! DESCRIPTION:
  !-------------
  ! The function calculates the vertical distance which a falling particle
  ! travels during one time step. As the vertical coordinate is the pressure,
  ! this distance is a pressure difference.
  !
  !-----------------
  ! INPUT VARIABLES:
  !-----------------
  ! TimeStep = duration of time step / s
  ! PRESS    = air pressure within grid box / Pa
  ! TEMP     = air temperature / K
  ! vel      = velocity of fall of aerosol particles / (m/s)
  !
  !--------
  ! OUTPUT:
  !--------
  ! mz_psc_SedStep  = sedimentation distance within one TimeStep / Pa
  !------------------------------------------------------------------------
  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TimeStep, PRESS, TEMP, vel
  REAL(dp) :: mz_psc_SedStep

  REAL(dp), PARAMETER :: MolMassAir=0.02897_dp !molar mass of dry air / (kg/mol)

  mz_psc_SedStep=(g_acc*MolMassAir*PRESS*vel*TimeStep)/(Rgas*TEMP)
END FUNCTION mz_psc_SedStep

!=========================================================================

PURE FUNCTION mz_psc_sed(kproma, klev,             &
                         pTop,                     &
                         pBottom,                  &
                         SedStep,                  &
                         InFrac,                   &
                         val_sed)
  !-------------
  ! DESCRIPTION:
  !-------------
  ! This function calculates the changes in amount-of-substance ratios of 
  ! H2O or HNO3 to air due to ice or NAT sedimentation. 
  ! Three different sedimentation schemes are available. 
  !
  !-----------------
  ! INPUT VARIABLES:
  !-----------------
  ! klev          = number of grid boxes within one column
  ! kproma        = number of columns
  ! pTop(:)       = air pressure at the top of grid box / Pa
  !   Note: The air pressure at the top of grid box i equals the air pressure
  !   at the bottom of the grid box above, which is grid box (i-1).
  ! pBottom(:)    = air pressure at the bottom of grid box / Pa
  !   Note: The air pressure at the bottom of grid box i equals the air
  !   pressure at the top of the grid box below, which is grid box (i+1).
  ! SedStep(:)    = sedimentation distance within one TimeStep / Pa
  ! InFrac(:)     = amount-of-substance ratio of X to air / (mol/mol)
  !   Note: X is ice phase H2O or NAT phase HNO3
  ! val_sed(:)    = flag indication grid boxes where sedimentation takes place
  !
  !------------------
  ! OUTPUT VALUES:
  !------------------
  ! mz_psc_sed(:) = change of amount-of-substance ratio of X to air 
  !                / (mol/mol)
  !   Note: X is ice phase H2O or NAT phase HNO3
  !------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(in) :: kproma, klev
  REAL(dp), DIMENSION(kproma,klev), INTENT(in) :: &
    pTop, pBottom, SedStep, InFrac
  LOGICAL, DIMENSION(kproma,klev), INTENT(in) :: val_sed

  REAL(dp), DIMENSION(kproma,klev) :: mz_psc_sed

  IF (SedScheme==1) THEN
    mz_psc_sed=SimpleUpwind(kproma, klev,    &
                            pTop,            &
                            pBottom,         &
                            SedStep,         &
                            InFrac,          &
                            val_sed)
  ELSEIF (SedScheme==2) THEN
    mz_psc_sed=Walcek2000(kproma, klev,      &
                          pTop,              &
                          pBottom,           &
                          SedStep,           &
                          InFrac,            &
                          val_sed)
  ELSEIF (SedScheme==3) THEN
    mz_psc_sed=TrapezoidScheme(kproma, klev, &
                               pTop,         &
                               pBottom,      &
                               SedStep,      &
                               InFrac,       &
                               val_sed)
  ELSE
    mz_psc_sed=0.0_dp
  END IF


CONTAINS

  !---------------------------------------------------------------------------
  PURE FUNCTION SimpleUpwind(kproma, klev,             &
                             pTop,                     &
                             pBottom,                  &
                             SedStep,                  &
                             InFrac,                   &
                             val_sed)
    !-------------
    ! DESCRIPTION:
    !------------- 
    ! Calculates sedimentation the easiest way.  
    ! Note that no sedimentation takes place from the highest grid box. This 
    ! neglection is a quick-and-dirty trick to avoid problems that occur when
    ! the upper end of the highest grid box has the pressure value pTop(1)=0.
    !
    !-----------------
    ! INPUT VARIABLES:
    !-----------------
    ! klev          = number of grid boxes within one column
    ! kproma        = number of columns
    ! pTop(:)       = air pressure at the top of grid box / Pa
    !   Note: The air pressure at the top of grid box i equals the air pressure
    !   at the bottom of the grid box above, which is grid box (i-1).
    ! pBottom(:)    = air pressure at the bottom of grid box / Pa
    !   Note: The air pressure at the bottom of grid box i equals the air
    !   pressure at the top of the grid box below, which is grid box (i+1).
    ! SedStep(:)    = sedimentation distance within one TimeStep / Pa
    ! InFrac(:)     = amount-of-substance ratio of X to air / (mol/mol)
    !   Note: X is ice phase H2O or NAT phase HNO3
    ! val_sed(:)    = flag indication grid boxes where sedimentation takes place
    !
    !------------------ 
    ! OUTPUT VALUES:
    !------------------ 
    ! SimpleUpwind(:) = change of amount-of-substance ratio of X to air 
    !                  / (mol/mol)
    !   Note: X is ice phase H2O or NAT phase HNO3
    !------------------------------------------------------------------------

    INTEGER, INTENT(in) :: kproma, klev
    REAL(dp), DIMENSION(kproma,klev), INTENT(in) :: &
      pTop, pBottom, SedStep, InFrac
    LOGICAL, DIMENSION(kproma,klev), INTENT(in) :: val_sed

    INTEGER :: i,j
    REAL(dp), DIMENSION(kproma,klev) :: frac, change, SimpleUpwind
    LOGICAL, DIMENSION(kproma,klev) :: val_loc

    INTRINSIC max

    frac=max(InFrac, 0.0_dp)
    change=0.0_dp

    val_loc = val_sed
    val_loc(1:kproma,1) = .false.
    val_loc(1:kproma,klev) = .false.

    FORALL (i=2:klev,j=1:kproma, val_loc(j,i) .AND. (.NOT. val_loc(j,i-1)))
      change(j,i)=-frac(j,i)*SedStep(j,i)             &
                        /(pBottom(j,i)-pTop(j,i))
    END FORALL

    FORALL (i=2:klev,j=1:kproma, val_loc(j,i) .AND. val_loc(j,i-1))
      change(j,i)=+frac(j,i-1)*SedStep(j,i-1)         &
                          /(pBottom(j,i)-pTop(j,i)) &
                -frac(j,i)*SedStep(j,i)             &
                        /(pBottom(j,i)-pTop(j,i))
    END FORALL

    FORALL (i=2:klev,j=1:kproma, (.NOT. val_loc(j,i)) .AND. val_loc(j,i-1))
      change(j,i)=+frac(j,i-1)*SedStep(j,i-1)   &
                          /(pBottom(j,i)-pTop(j,i))
    END FORALL

    SimpleUpwind=MAX(change,-InFrac)
  END FUNCTION SimpleUpwind
  !---------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  PURE FUNCTION Walcek2000(kproma, klev,             &
                           pTop,                     &
                           pBottom,                  &
                           SedStep,                  &
                           InFrac,                   &
                           val_sed)

    !-------------
    ! DESCRIPTION:
    !------------- 
    ! See literature: 
    ! Chris J. Walcek, "Minor flux adjustment near mixing ratio extremes for
    ! simplified yet highly accurate monotonic calculation of tracer advection";
    ! Journal of Geophysical Research, Vol. 105 (D7), pp. 9335-9348, 2000
    !
    ! The original program code published in the above paper has been modified 
    ! for the current purpose. The most important modifications are: 
    ! - Code transformed from Fortran 77 to Fortran 95
    ! - Subroutine replaced by function
    ! - Particle movement only in one direction (u >= 0)
    ! - Variable q0 replaced by variable InFrac
    ! - Variable qn replaced by variable frac
    ! - Courant number calculation
    !     x1= dt*u(I)/dxx(I)
    !   replaced by
    !     x1 = log(pBottom(i)/(pBottom(i)-SedStep(i)))/log(pBottom(i)/pTop(i))
    ! - Use only one density instead of den0, den1, and dd0
    ! - Replace
    !     dxx(I)*den0(i)
    !   by
    !     (pBottom(i)-pTop(i))*N_Avogadro/(g_acc*MolMassAir)
    ! - Use array dimension 1:klev
    ! - Make sedimentation calculation dependent on val_sed
    !
    !-----------------
    ! INPUT VARIABLES:
    !-----------------
    ! klev          = number of grid boxes within one column
    ! kproma        = number of columns
    ! pTop(:)       = air pressure at the top of grid box / Pa
    !   Note: The air pressure at the top of grid box i equals the air pressure
    !   at the bottom of the grid box above, which is grid box (i-1).
    ! pBottom(:)    = air pressure at the bottom of grid box / Pa
    !   Note: The air pressure at the bottom of grid box i equals the air
    !   pressure at the top of the grid box below, which is grid box (i+1).
    ! SedStep(:)    = sedimentation distance within one TimeStep / Pa
    ! InFrac(:)     = amount-of-substance ratio of X to air / (mol/mol)
    !   Note: X is ice phase H2O or NAT phase HNO3
    ! val_sed(:)    = flag indication grid boxes where sedimentation takes place
    !
    !------------------ 
    ! OUTPUT VALUES:
    !------------------ 
    ! Walcek2000(:) = change of amount-of-substance ratio of X to air 
    !                / (mol/mol)
    !   Note: X is ice phase H2O or NAT phase HNO3
    !------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kproma, klev
    REAL(dp), DIMENSION(kproma,klev), INTENT(in) :: &
      pTop, pBottom, SedStep, InFrac
    LOGICAL, DIMENSION(kproma,klev), INTENT(in) :: val_sed

    ! molar mass of dry air / (kg/mol): 
    REAL(dp), PARAMETER :: MolMassAir=0.02897_dp 

    INTEGER :: i,j
    LOGICAL, DIMENSION(kproma,klev) :: imxmn
    REAL(dp) :: cf, cf1, x1, x1n
    REAL(dp), DIMENSION(kproma,klev) :: frac, flux, vcmax, vcmin, Walcek2000
    REAL(dp) :: rhelp

    INTRINSIC log, max, min

    ! Identify local max and min, specify mixing ratio limits at new time
    ! (equation 7)
    imxmn(1:kproma,1)       = .false.
    imxmn(1:kproma,2:klev-1)=InFrac(1:kproma,2:klev-1)>= &
                     MAX(InFrac(1:kproma,1:klev-2),InFrac(1:kproma,3:klev)) &
                     .OR. InFrac(1:kproma,2:klev-1)<= &
                     MIN(InFrac(1:kproma,1:klev-2),InFrac(1:kproma,3:klev))
    imxmn(1:kproma,klev)    = .false.

    vcmax(1:kproma,1)       = InFrac(1:kproma,1)
    vcmin(1:kproma,1)       = InFrac(1:kproma,1)
    vcmax(1:kproma,2:klev)  = MAX(InFrac(1:kproma,2:klev), &
         InFrac(1:kproma,1:klev-1))
    vcmin(1:kproma,2:klev)  = MIN(InFrac(1:kproma,2:klev), &
         InFrac(1:kproma,1:klev-1))

    ! Initialize mixing ratios and fluxes with zero
    frac= 0.0_dp
    flux= 0.0_dp
    rhelp=g_acc*MolMassAir/N_Avogadro
    DO i=2,klev-1
      DO j=1,kproma
        IF (val_sed(j,i)) THEN
          ! Courant number: 
          x1 = log(pBottom(j,i)/(pBottom(j,i)-SedStep(j,i)))/ &
               log(pBottom(j,i)/pTop(j,i))
          x1n= (1.0_dp-x1)*(InFrac(j,i+1)-InFrac(j,i-1))/4.0_dp

          ! estimate mixing ratio in outgoing fluid
          IF (imxmn(j,i-1)) THEN       ! see equation 10b
            cf = InFrac(j,i)+MAX(1.5_dp, 1.2_dp+0.6_dp*x1)*x1n
          ELSEIF (imxmn(j,i+1)) THEN   ! see equation 10a
            cf = InFrac(j,i)+(1.75_dp -0.45_dp*x1)*x1n
          ELSE                       ! see equation 4a
            cf = InFrac(j,i) + x1n
          END IF

          ! Limit cf to be between mixing ratio on either side of edge where flux 
          ! is being calculated
          cf1= MIN(MAX(cf,MIN(InFrac(j,i),InFrac(j,i+1))), &
               MAX(InFrac(j,i),InFrac(j,i+1)))

          ! Calculate mixing ratio at new time
          frac(j,i)= InFrac(j,i) &
                  -x1*cf1    &
                  +flux(j,i-1)*rhelp/(pBottom(j,i)-pTop(j,i))
!         +flux(j,i-1)*g_acc*MolMassAir/((pBottom(j,i)-pTop(j,i))*N_Avogadro)

          ! Now use vcmax and vcmin as physical limits for frac. 
          ! This limitation is necessary to avoid oscillatory behaviour of the
          ! solution, including negative values; however, it limits the
          ! applicability of the Walcek (2000) scheme for sedimentation. 
          frac(j,i)=MAX(vcmin(j,i), MIN(vcmax(j,i), frac(j,i)))

          ! Re-calculate OUTFLOWING flux before moving on to next cell
!     flux(j,i)= (pBottom(j,i)-pTop(j,i))*N_Avogadro*(InFrac(j,i)-frac(j,i)) &
!                /(g_acc*MolMassAir) &
          flux(j,i)= (pBottom(j,i)-pTop(j,i))*(InFrac(j,i)-frac(j,i)) &
                   /rhelp &
                  +flux(j,i-1)
        ELSEIF (val_sed(j,i-1)) THEN
          frac(j,i)= InFrac(j,i) &
                  +flux(j,i-1)*rhelp/(pBottom(j,i)-pTop(j,i))
!          +flux(j,i-1)*g_acc*MolMassAir/((pBottom(j,i)-pTop(j,i))*N_Avogadro)
        END IF
      END DO
    END DO

    DO j=1,kproma
      IF (val_sed(j,klev-1)) THEN
        frac(j,klev)= InFrac(j,klev) &
                   +flux(j,klev-1)*rhelp &
                    /(pBottom(j,klev)-pTop(j,klev))
!                   +flux(j,klev-1)*g_acc*MolMassAir &
!                    /((pBottom(j,klev)-pTop(j,klev))*N_Avogadro)
      END IF
    END DO

    Walcek2000 = MAX(frac-InFrac, -InFrac)
  END FUNCTION Walcek2000
  !---------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ELEMENTAL FUNCTION Trapezoid(frac_above, pTop_above,    &
                               frac_i, pTop_i, pBottom_i, &
                               frac_below, pBottom_below, &
                               SedStep)
    !-------------
    ! DESCRIPTION:
    !-------------
    ! The elemental function Trapezoid is used by the pure function
    ! TrapezoidScheme. 
    ! 
    ! The simplest algorithm to calculate the changes in the amount-of-
    ! substance ratios of H2O or HNO3 to dry air due to sedimentation
    ! ("simple upwind scheme") includes in its equations the product of the
    ! amount-of-substance ratio and a pressure difference as measure for the
    ! height of the vertical air layer from which particles fall into the next
    ! box.
    ! This simple algorithm has the disadvantage of "numerical diffusion".
    ! The problem of numerical diffusion can be reduced if the product of
    ! the amount-of-substance ratio and the sedimentation step
    ! is replaced by a more sophisticated calculation. The function Trapezoid
    ! uses integrals over straight line approximations of the amount-of-
    ! substance ratio distribution for that purpose. Definite integrals
    ! over straight lines have trapezoidal shape, hence the function name.
    ! Where SedStep is greater than half the box height, the calculation
    ! becomes more difficult; the function does not return just the area
    ! of a trapezoid but rather the sum of the areas of a trapezoid and a
    ! rectangle.
    !
    !-------
    ! INPUT:
    !-------
    ! frac_above    = amount-of-substance ratio in grid box i-1 / (mol/mol)
    ! pTop_above    = pressure at the top of grid box i-1 / Pa
    ! frac_i        = amount-of-substance ratio in grid box i / (mol/mol)
    ! pTop_i        = pressure at the top of grid box i / Pa
    !               = pressure at the bottom of grid box i-1 / Pa
    ! pBottom_i     = pressure at the bottom of grid box i / Pa
    !               = pressure at the top of grid box i+1 / Pa
    ! frac_below    = amount-of-substance ratio in grid box i+1 / (mol/mol)
    ! pBottom_below = pressure at the bottom of grid box i+1 / Pa
    ! SedStep       = sedimentation distance within one time step / Pa
    ! 
    !--------
    ! OUTPUT:
    !--------
    ! Trapezoid = product of an amount-of-substance ratio and a pressure
    !             difference
    !--------- 
    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: &
      frac_above, pTop_above,    &
      frac_i, pTop_i, pBottom_i, &
      frac_below, pBottom_below, &
      SedStep
    REAL(dp) :: &
      Trapezoid, Trapezoid_tmp, &
      SlopeAbove, InterceptAbove, SlopeBelow, InterceptBelow, &
      pHeight_i

    INTRINSIC min

    pHeight_i=pBottom_i-pTop_i
    Trapezoid=0.0_dp

    IF (frac_above<=frac_i .AND. frac_i<=frac_below &
        .AND. SedStep<=0.5_dp*pHeight_i) THEN
      !----------------------------------------------------------------------
      ! If the amount-of-substance ratio "frac" increases from box i-1 to 
      ! box i and from box i to box i+1, the straight line approximation method
      ! increases the amount of transported particles compared to
      ! simple volume shifting.
      !----------------------------------------------------------------------
      SlopeBelow=(frac_i-frac_below)/(0.5_dp*(pTop_i-pBottom_below))
! op_pj_20100714+
!!$      InterceptBelow=frac_i-0.5_dp*(pBottom_i-pTop_i)*SlopeBelow
      InterceptBelow=frac_i-0.5_dp*(pBottom_i+pTop_i)*SlopeBelow
! op_pj_20100714-
      Trapezoid=SedStep                        &
                *(SlopeBelow*pBottom_i         &
                  -0.5_dp*SlopeBelow*SedStep   &
                  +InterceptBelow)
      IF (Trapezoid>pHeight_i*frac_i) Trapezoid=pHeight_i*frac_i
    ELSEIF (frac_above<=frac_i .AND. frac_i<=frac_below &
            .AND. SedStep>0.5_dp*pHeight_i) THEN
      !----------------------------------------------------------------------
      ! If the amount-of-substance ratio "frac" increases from box i-1 to 
      ! box i and from box i to box i+1, the straight line approximation method
      ! increases the amount of transported particles compared to
      ! simple volume shifting.
      ! However, where the straight line values are lower than the average
      ! frac-values in box i, the average frac-values are integrated instead
      ! of the line.
      !----------------------------------------------------------------------
      SlopeBelow=(frac_i-frac_below)/(0.5_dp*(pTop_i-pBottom_below))
! op_pj_20100714+
!!$      InterceptBelow=frac_i-0.5_dp*(pBottom_i-pTop_i)*SlopeBelow
      InterceptBelow=frac_i-0.5_dp*(pBottom_i+pTop_i)*SlopeBelow
! op_pj_20100714-
      Trapezoid= 0.5_dp*pHeight_i                       &
                 *(SlopeBelow*pBottom_i                 &
                   -0.25_dp*SlopeBelow*pHeight_i        &
                   +InterceptBelow)                     &
                +frac_i*(SedStep-0.5_dp*pHeight_i)
      IF (Trapezoid>pHeight_i*frac_i) Trapezoid=pHeight_i*frac_i
    ELSEIF (frac_above>frac_i .AND. frac_i<=frac_below) THEN
      !---------------------------------------------------------------------
      ! If box i is a local minimum, sedimentation is calculated like in
      ! the simple upwind scheme
      !---------------------------------------------------------------------
      Trapezoid=frac_i*MIN(SedStep,pHeight_i)
    ELSEIF (frac_above<=frac_i .AND. frac_i>frac_below) THEN
      !---------------------------------------------------------------------
      ! If box i is a local maximum, sedimentation is calculated like in
      ! the simple upwind scheme
      !---------------------------------------------------------------------
      Trapezoid=frac_i*MIN(SedStep,pHeight_i)
    ELSEIF (frac_above>frac_i .AND. frac_i>frac_below &
            .AND. SedStep<=0.5_dp*pHeight_i) THEN
      !---------------------------------------------------------------------
      ! If the amount-of-substance ratio "frac" decreases between the 
      ! boxes i-1 and i as well as between the boxes i and i+1, numerical
      ! diffusion correction by integrating straight line approximation can be
      ! applied in two ways. The one that yields the smaller Trapezoid is
      ! choosen.
      !---------------------------------------------------------------------
      SlopeAbove=(frac_above-frac_i)/(0.5_dp*(pTop_above-pBottom_i))
! op_pj_20100714+
!!$      InterceptAbove=frac_i-0.5_dp*(pBottom_i-pTop_i)*SlopeAbove
      InterceptAbove=frac_i-0.5_dp*(pBottom_i+pTop_i)*SlopeAbove
! op_pj_20100714-
      Trapezoid_tmp=SedStep                        &
                    *(SlopeAbove*pBottom_i         &
                      -0.5_dp*SlopeAbove*SedStep   &
                      +InterceptAbove)
      ! op_pj_20100716+
      ! Note that Trapezoid_tmp can become negative, if the straight line
      ! intercepts the frac=0 line within the grid box i. 
      ! In this case the calculations below will result in a zero 
      ! sedimentation. The fix here is to approximate the mixing ratio
      ! distribution within the box by a triangle and to sediment the
      ! lower part (of height SedStep).
      IF (Trapezoid_tmp < 0.0_dp) &
           Trapezoid_tmp = 0.5_dp * SedStep * SedStep*(frac_i/(0.5*pHeight_i))
      ! op_pj_20100716-
      SlopeBelow=(frac_i-frac_below)/(0.5_dp*(pTop_i-pBottom_below))
! op_pj_20100714+
!!$      InterceptBelow=frac_i-0.5_dp*(pBottom_i-pTop_i)*SlopeBelow
      InterceptBelow=frac_i-0.5_dp*(pBottom_i+pTop_i)*SlopeBelow
! op_pj_20100714-
      Trapezoid=SedStep                        &
                *(SlopeBelow*pBottom_i         &
                  -0.5_dp*SlopeBelow*SedStep   &
                  +InterceptBelow)
      IF (Trapezoid_tmp<Trapezoid) Trapezoid=Trapezoid_tmp
    ELSEIF (frac_above>frac_i .AND. frac_i>frac_below &
            .AND. SedStep>0.5_dp*pHeight_i) THEN
      !---------------------------------------------------------------------
      ! If the amount-of-substance ratio "frac" decreases between the boxes 
      ! i-1 and i as well as between the boxes i and i+1, numerical diffusion
      ! correction by integrating straight line approximation can be applied
      ! in two ways. The one that yields the smaller Trapezoid is
      ! choosen.
      ! Again, where the straigt line values are higher than the
      ! average concentration values in box i, the average concentration
      ! values are integrated instead of the line.
      !---------------------------------------------------------------------
      SlopeAbove=(frac_above-frac_i)/(0.5_dp*(pTop_above-pBottom_i))
! op_pj_20100714+
!!$      InterceptAbove=frac_i-0.5_dp*(pBottom_i-pTop_i)*SlopeAbove
      InterceptAbove=frac_i-0.5_dp*(pBottom_i+pTop_i)*SlopeAbove
! op_pj_20100714-
      Trapezoid_tmp= 0.5_dp*pHeight_i                  &
                     *(SlopeAbove*pBottom_i            &
                       -0.25_dp*SlopeAbove*pHeight_i   &
                       +InterceptAbove)                &
                    +frac_i*(SedStep-0.5_dp*pHeight_i)
      ! op_pj_20100716+
      ! Note that Trapezoid_tmp can become negative, if the straight line
      ! intercepts the frac=0 line within the grid box i. 
      ! In this case the calculations below will result in a zero 
      ! sedimentation. The fix here is to approximate the mixing ratio
      ! distribution within the box by a triangle and to sediment the
      ! lower part (of height pHeight_i/2) plus the rectangle above
      ! (of height SedStep-pHeight_i).
      IF (Trapezoid_tmp < 0.0_dp) &
           Trapezoid_tmp = (SedStep - 0.25_dp*pHeight_i) * frac_i
      ! op_pj_20100716-
      SlopeBelow=(frac_i-frac_below)/(0.5_dp*(pTop_i-pBottom_below))
! op_pj_20100714+
!!$      InterceptBelow=frac_i-0.5_dp*(pBottom_i-pTop_i)*(frac_i-frac_below) &
!!$                            /(0.5_dp*(pTop_i-pBottom_below))
      InterceptBelow=frac_i-0.5_dp*(pBottom_i+pTop_i)*SlopeBelow
! op_pj_20100714-
      Trapezoid= 0.5_dp*pHeight_i                  &
                 *(SlopeBelow*pBottom_i            &
                   -0.25_dp*SlopeBelow*pHeight_i   &
                   +InterceptBelow)                &
                +frac_i*(SedStep-0.5_dp*pHeight_i)
      IF (Trapezoid_tmp<Trapezoid) Trapezoid=Trapezoid_tmp
    END IF
    IF (Trapezoid<0.0_dp) Trapezoid=0.0_dp
  END FUNCTION Trapezoid
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  PURE FUNCTION TrapezoidScheme(kproma, klev,             &
                                pTop,                     &
                                pBottom,                  &
                                SedStep,                  &
                                InFrac,                   &
                                val_sed)
    !-------------
    ! DESCRIPTION:
    !-------------
    ! This function calculates the changes in amount-of-substance ratios of H2O
    ! or HNO3 molecules to dry air due to sedimentation. 
    ! The idea of the sedimentation algorithm is to calculate the amount of
    ! sedimenting particles not by integrating the vertical amount-of-substance
    ! ratio profile profile but local straight line approximations to that step
    ! function.
    ! No sedimentation takes place from the highest grid box and from the lowest
    ! grid box. 
    !
    !-----------------
    ! INPUT VARIABLES:
    !-----------------
    ! klev          = number of grid boxes within one column
    ! kproma        = number of columns
    ! pTop(:)       = air pressure at the top of grid box / Pa
    !   Note: The air pressure at the top of grid box i equals the air pressure
    !   at the bottom of the grid box above, which is grid box (i-1).
    ! pBottom(:)    = air pressure at the bottom of grid box / Pa
    !   Note: The air pressure at the bottom of grid box i equals the air
    !   pressure at the top of the grid box below, which is grid box (i+1).
    ! SedStep(:)    = sedimentation distance within one TimeStep / Pa
    ! InFrac(:)     = amount-of-substance ratio of X to air / (mol/mol)
    !   Note: X is ice phase H2O or NAT phase HNO3
    ! val_sed(:)    = flag indication grid boxes where sedimentation takes place
    !
    !------------------
    ! OUTPUT VALUES:
    !------------------
    ! TrapezoidScheme(:) = change of amount-of-substance ratio of X to air 
    !                     / (mol/mol)
    !   Note: X is ice phase H2O or NAT phase HNO3
    !------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kproma, klev
    REAL(dp), DIMENSION(kproma, klev), INTENT(in) :: &
      pTop, pBottom, SedStep, InFrac
    LOGICAL, DIMENSION(kproma, klev), INTENT(in) :: val_sed

    INTEGER :: i,j
    REAL(dp), DIMENSION(kproma,klev) :: frac, change, TrapezoidScheme
    LOGICAL, DIMENSION(kproma,klev) :: val_loc

    INTRINSIC MAX

    frac=max(InFrac, 0.0_dp)
    change=0.0_dp

    val_loc = val_sed
    val_loc(1:kproma,1) = .false.
    val_loc(1:kproma,klev) = .false.

    FORALL (i=2:klev-1,j=1:kproma, val_loc(j,i) .AND. (.NOT. val_loc(j,i-1)))
      change(j,i)=-Trapezoid(frac(j,i-1), pTop(j,i-1),               &
                           frac(j,i), pTop(j,i), pBottom(j,i),       &
                           frac(j,i+1), pBottom(j,i+1),              &
                           SedStep(j,i))                             &
                 /(pBottom(j,i)-pTop(j,i))
    END FORALL

    FORALL (i=3:klev-1,j=1:kproma, val_loc(j,i) .AND. val_loc(j,i-1))
      change(j,i)= Trapezoid(frac(j,i-2), pTop(j,i-2),               &
                           frac(j,i-1), pTop(j,i-1), pBottom(j,i-1), &
                           frac(j,i), pBottom(j,i),                  &
                           SedStep(j,i-1))                           &
                 /(pBottom(j,i)-pTop(j,i))                           &
                -Trapezoid(frac(j,i-1), pTop(j,i-1),                 &
                           frac(j,i), pTop(j,i), pBottom(j,i),       &
                           frac(j,i+1), pBottom(j,i+1),              &
                           SedStep(j,i))                             &
                 /(pBottom(j,i)-pTop(j,i))
    END FORALL

    FORALL (i=3:klev,j=1:kproma, (.NOT. val_loc(j,i)) .AND. val_loc(j,i-1))
      change(j,i)= Trapezoid(frac(j,i-2), pTop(j,i-2),               &
                           frac(j,i-1), pTop(j,i-1), pBottom(j,i-1), &
                           frac(j,i), pBottom(j,i),                  &
                           SedStep(j,i-1))                           &
                 /(pBottom(j,i)-pTop(j,i))
    END FORALL

    TrapezoidScheme=MAX(change,-InFrac)
  END FUNCTION TrapezoidScheme
  !--------------------------------------------------------------------------


END FUNCTION mz_psc_sed

!=========================================================================

! ka_ok_20100115+
!=========================================================================
SUBROUTINE K_PARA_NAT(HNO3_nsarr_00, HNO3_nsarr_s, DHNO3, &
     vhno3, pHNO3, pHNO3_over_NAT, &
     pHNO3S, svapn, TEMP, PRESS, dtime, radius_NAT, &
     dens_nat, r_NAT_arr, G)

  IMPLICIT NONE

  REAL(dp), INTENT (in) ::  DHNO3, vhno3, pHNO3, pHNO3_over_NAT, &
       pHNO3S, svapn, TEMP, PRESS, dtime 
  !!!REAL(dp), INTENT (inout) :: HNO3
  
  REAL(dp), DIMENSION(NSB), INTENT (in) :: HNO3_nsarr_00
  REAL(dp), DIMENSION(NSB), INTENT (out) :: HNO3_nsarr_s
  REAL(dp), DIMENSION(NSB), INTENT (out) :: dens_nat
  REAL(dp), DIMENSION(NSB), INTENT (out) :: r_NAT_arr
  REAL(dp), DIMENSION(NSB), INTENT (out) :: G
  REAL(dp), INTENT (out) :: radius_NAT
  ! Minimum radii and particle density of the NAT size bins

  REAL(dp), DIMENSION(NSB), PARAMETER :: &
       rbin_min=(/0.0,0.2E-6,1.E-6,2.E-6,6.E-6,9.E-6,12.E-6,&
       16.E-6,20.E-6/) ! min radius in m (9bins)
  REAL(dp), DIMENSION(NSB), PARAMETER :: &
       rbin_max=(/0.2E-6,1.E-6,2.E-6,6.E-6,9.E-6,12.E-6,16.E-6,&
       20.E-6,25.E-6/) ! max radius in m (9bins)
  REAL(dp), DIMENSION(NSB), PARAMETER :: &
       rbin_av =(/0.1E-6,0.6E-6,1.5E-6,4.E-6,7.5E-6,10.5E-6,14.0E-6,&
       18.0E-6,22.5E-6/) !  average radius in m (9bins) 
  REAL(dp), DIMENSION(NSB), PARAMETER :: &
       no_density_limit = (/3.2857e-5,3.2857e-5,3.2857e-5,3.2857e-5,3.2857e-5, &
       3.2857e-5,1.64785e-5,1.64785e-5,1.64785e-5/) ! 9 bins
  REAL(dp), DIMENSION(NSB) :: m_nat, m_nat_old, dens_diff, &
       mass_diff ,mass_ini , nat_diff
  REAL(dp), DIMENSION(NSB) :: new_mass, rad_0, nat_rad, x, HNO3_nsarr_0 
  REAL(dp) :: arg, radius, factor, total_dens_nat
  INTEGER :: n
  
  INTRINSIC :: SUM
  
  REAL(dp), parameter :: r_mol_HNO3=3.0e-10_dp 
  REAL(dp) :: mfp_HNO3, visc_air ! only test
  REAL(dp) :: Cc, SedFactor ! only test
  REAL(dp), Dimension(NSB) :: NAT_vel ! only test
  
  factor = mix2conc(PRESS, TEMP)
  HNO3_nsarr_0(:)=factor*HNO3_nsarr_00(:)
  
  DO n=1,NSB ! loop over NSB bins
     m_nat(n) = 0.0_dp ! this is re-calculated after the modification in radius due to growth/shrinkage. 
     G(n) = 0.0_dp ! sets growth value to zero for size bin (n) for each time step.
     dens_diff(n) = 0.0_dp ! initialization of xs particle no.density in size bin (n)
     mass_diff(n) = 0.0_dp ! initialization of mass to be transferred to next size bin
     mass_ini(n)=0.0_dp ! initialization of initial density and initial mass in size bin after advection.
     dens_nat(n) = 0.0_dp ! initialization of particle number density before re-calculation
     new_mass = 0.0_dp

     IF (HNO3_nsarr_0(n)<0.0_dp) THEN 
        !     write(*,*) 'ERROR -> negative initial particle number of ', HNO3_nsarr_0(n)
        !     write(*,*) 'SET it to 0.0_dp'
        HNO3_nsarr_0(n) = 0.0_dp
        !     write(*,*) 'HNO3_nsarr_0(n)', HNO3_nsarr_0(n)
     ENDIF

     rad_0(n)=rbin_av(n) ! particle radius at t=0

     nat_rad(n)=rbin_av(n) ! NAT-radius

     ! only to write out the amounts of sedimentation calculation

     m_nat_old(n) = 4.0_dp/3.0_dp*pi*rho_NAT*nat_rad(n)**3.0_dp  ! mass per particle [g] ; van den Broek et al. 1871-2

     ! For size bins > 1, update mass with transfered mass from previous size bin
     mass_ini(n) = HNO3_nsarr_0(n)*MolMassNAT/N_Avogadro ! initial mass in size bin at start of de-nit step for bin > 1 [g]
     dens_nat(n) = (HNO3_nsarr_0(n)*MolMassNAT)/(m_nat_old(n)*N_Avogadro) ! calculates no.density of particles for bin > 1       

     ! Initialization of NAT particlea at start of simulation
     IF ((n==1) .and. (dens_nat(1) .lt. 1.5e-5_dp) .and. (svapn .ge. 1.)) THEN
        dens_nat(1) = 1.5e-5_dp ! 1/cm**3
        nat_rad(1) = rbin_av(1)
     END IF

     G(n)=growth_factor(DHNO3, vhno3, nat_rad(n), pHNO3, pHNO3_over_NAT, TEMP)

     !--------------------------------
     !In case that the size decrease is larger than the initial radius then the particle has completely evaporated
     !---------------------------------
     IF (dens_nat(n) .ge. 1e-10_dp) THEN !only grow particles if there are any present
        arg = nat_rad(n)**2.0_dp + 2.0_dp*G(n)*dtime
        IF (arg .ge. 0.) radius = sqrt(arg)
        IF (arg .lt. 0.) radius = 0.0_dp 
        nat_rad(n) = radius  
        ! determine new particle mass after growth
        m_nat(n) = 4.0_dp/3.0_dp*pi*rho_NAT*nat_rad(n)**3.0_dp               
        ! calculates mass per particle [kg] for new radius
        ! here the particle no. density remains the same (i.e. ) all particles experience identical radius change
        ! updates the number of particles
     ELSE                       
        m_nat(n) = m_nat_old(n) 
     END IF ! no. density filter

     x(n) = 0.0_dp 
     HNO3_nsarr_s(n) = N_Avogadro*m_nat(n)*dens_nat(n)/MolMassNAT

     x(n) = HNO3_nsarr_s(n) - HNO3_nsarr_0(n) ! difference new to old concentration 
     !----------------------------------------------------------------
     ! Each sizebin is 'refilled' in the next section
     !----------------------------------------------------------------
     dens_nat(n) = (HNO3_nsarr_s(n)*MolMassNAT)/(m_nat_old(n)*N_Avogadro)
     new_mass(n) = HNO3_nsarr_s(n)*MolMassNAT/N_Avogadro

     !---------------------------------------------------------------
     ! Transfer mass from the previous size bin if necessary
     !----------------------------------------------------------------
     IF (n>1) THEN
        IF (mass_diff(n-1)>0.0_dp) THEN
           new_mass(n) = new_mass(n)+mass_diff(n-1)
           dens_nat(n) = new_mass(n)/m_nat_old(n) ! calculate new particles no. density from added mass
           HNO3_nsarr_s(n) = N_Avogadro*m_nat_old(n)*dens_nat(n)/MolMassNAT                                    
        END IF
     END IF
     !------------------------------------------------------------------------------------------               
     ! Checks to see if the particle number density exceeds the threshold set for each size bin      
     ! If so, the xs particles are essentially stored and added, as mass, to the next size bin up.
     ! No growth of stored particles occurs till they are added to the next bin.
     IF (dens_nat(n)>no_density_limit(n)) THEN
        dens_diff(n) = dens_nat(n)-no_density_limit(n) ! how much particles over the limit
        nat_diff(n) = (N_Avogadro * dens_diff(n) * m_nat_old(n))/MolMassNAT ! convert to  
        mass_diff(n) = nat_diff(n)*MolMassNAT/N_Avogadro ! [kg]
        IF (nat_diff(n)>1.e-30_dp) THEN! filter to ensure -ve NAT concentrations don't occur
           HNO3_nsarr_s(n) = (N_Avogadro * no_density_limit(n) * m_nat_old(n))/MolMassNAT  !modify no.particles in bin (n)
           ! Calculate a new number density after XS particles stored
           dens_nat(n) = (HNO3_nsarr_s(n)*MolMassNAT)/(m_nat_old(n) * N_Avogadro) 
        ELSE
           mass_diff(n) = 0.0_dp
        END IF
     END IF

     !--------------------------------------------------------------------------------------------------------------------
     ! Warning given if an overflow of the last size bin occurs (during prolonged low temps)         
     IF ((n==NSB).and.(dens_nat(NSB)>no_density_limit(NSB))) write (*,*) ' BIN 9 in overflowing !!!!!! with density',dens_nat(NSB)
     !--------------------------------------------------------------------------------------------------------------------

     !----------------------------------------------------------------------------------------
     ! Warning if the maximum particle radius becomes too big (i.e) greater than 25e-6 m
     !----------------------------------------------------------------------------------------

     IF (nat_rad(n) > rbin_max(NSB)) then
        write (*,*) 'error in radius'
        write (*,*) 'nat_rad,radius',n,nat_rad(n),rbin_max(n)
     END IF
  END DO ! loop over 9 bins 

  ! calculation of the total number density and mean radius of NAT particles
  total_dens_nat=sum(dens_nat(:))
  IF (total_dens_nat > 5.0e-7_dp) THEN 
     radius_NAT=sqrt((dens_nat(1)*rbin_av(1)**2+dens_nat(2)*rbin_av(2)**2   &
          +dens_nat(3)*rbin_av(3)**2+dens_nat(4)*rbin_av(4)**2        &
          +dens_nat(5)*rbin_av(5)**2+dens_nat(6)*rbin_av(6)**2        &
          +dens_nat(7)*rbin_av(7)**2+dens_nat(8)*rbin_av(8)**2        &
          +dens_nat(9)*rbin_av(9)**2)/total_dens_nat)
  ELSE
     radius_NAT=0.0_dp
  ENDIF

  r_NAT_arr(:)=rbin_av(:)
  HNO3_nsarr_s(:)=HNO3_nsarr_s(:)/factor
END SUBROUTINE K_PARA_NAT
!=========================================================================

!==========================================================================
 ELEMENTAL FUNCTION mz_psc_N_ice(phase, TEMP, PRESS, H2O_ice)
  ! O.Kirner 06/11/23 from orginal mz_psc_N_solid 
  !-------------------------------------------------------------------------
  ! This function calculates the particle number density of ice particles.
  ! If there is much ice in a grid box, the function returns a constant 
  ! minimum value N_solid_max. However, if the input particle number density
  ! together with a rather small amount of ice in the grid box leads to
  ! very small particle radi, the particle number density is decreased, so
  ! that the particle radius equals a default minimum radius.
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! phase    = phase indicator (1=liquid, 2=liquid+nat, 3=liquid+nat+ice)
  ! TEMP     = temperature / K
  ! PRESS    = pressure / hPa
  ! H2O      = amount-of-substance ratio of gaseous H2O to dry air / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! m_ice    = mass of ice per volume of air / (kg/m**3)
  ! v_ice    = volume fraction of ice in air / (m**3/m**3)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_N_ice = number of solid particles per m**3 air
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in) :: phase
  REAL(dp), INTENT(in) :: TEMP, PRESS, H2O_ice

  REAL(dp) :: mz_psc_N_ice

  REAL(dp) :: m_ice, v_ice
  REAL(dp) :: factor

  INTRINSIC MIN

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/m**3)
  ! -----
  factor   = mix2conc(PRESS, TEMP)*1.0E6

  m_ice   =(H2O_ice*factor) * MolMassH2O/(N_Avogadro*1000.0_dp)
  v_ice   = m_ice / dens_ice
  mz_psc_N_ice = MIN(3.0_dp*v_ice/(4.0_dp*pi*r_min**3), N_solid_max)

!-------------------------------------------------------------------------
 END FUNCTION mz_psc_N_ice
!=========================================================================

!==========================================================================
ELEMENTAL FUNCTiON mz_psc_r_ice(phase, TEMP, PRESS, N_ice, H2O_ice)
  !-------------------------------------------------------------------------
  ! This function calculates the radius of solid particles containing
  ! ice and/or NAT.
  !-------------------------------------------------------------------------
  ! ------
  ! INPUT:
  ! ------
  ! phase    = phase indicator (1=liquid, 2=liquid+nat, 3=liquid+nat+ice)
  ! TEMP = temperature / K
  ! PRESS    = pressure / hPa
  ! N_ice    = number of ice particles per m**3 air
  ! H2O_ice  = amount-of-substance ratio of ice phase to dry air / (mol/mol)
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! m_ice = mass of ice per volume of air / (kg/m**3)
  ! v_ice = volume fraction of ice in air / (m**3/m**3)
  ! -------
  ! OUTPUT:
  ! -------
  ! mz_psc_r_ice = r_ice = radius of one single ice particle / m
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in) :: phase
  REAL(dp), INTENT(in) :: TEMP, PRESS, H2O_ice
  REAL(dp), INTENT(in) :: N_ice

  REAL(dp) :: mz_psc_r_ice, r_ice

  REAL(dp) :: m_ice, v_ice
  REAL(dp) :: factor
  
  INTRINSIC max, tiny

  ! -----
  ! convert from amount-of-substance ratio of HX to dry air / (mol/mol)
  ! to particle number concentration / (1/m**3)
  ! -----
  factor   = mix2conc(PRESS, TEMP)*1.0E6
  m_ice   = (H2O_ice*factor) * MolMassH2O/(N_Avogadro*1000.0_dp)
  v_ice   = m_ice / dens_ice
  r_ice = (3.0_dp*v_ice/max(4.0_dp*N_ice*pi,tiny(1.0_dp)))**(1.0_dp/3.0_dp)
  mz_psc_r_ice = r_ice
!-------------------------------------------------------------------------
 END FUNCTION mz_psc_r_ice
!=========================================================================

!=========================================================================
  ELEMENTAL FUNCTION diff_H2O(TEMP,PRESS)
     ! O.Kirner 06/11/08   
     ! Function to calculate the diffusion coefficient of H2O in air
     ! INPUT :
     ! TEMP = temperature; K
     ! PRESS = pressure; hPa
     ! OUTPUT :
     ! [diff_H2O]=m**2/s
     ! Reif, 1965
     
     IMPLICIT NONE
     
     REAL(dp), INTENT(in) :: TEMP,PRESS
     REAL(dp) :: diff_H2O
     
     diff_H2O=0.22_dp*((TEMP/273.15_dp)**1.94_dp)*(1013.25_dp/(PRESS/100_dp))*1.0E-4_dp 
  
  END FUNCTION diff_H2O
  
  ELEMENTAL FUNCTION diff_HNO3(DH2O)
     ! O.Kirner 06/11/08
     ! FunctIon to calculate the diffusion coefficient of HNO3 in air
     ! INPUT :
     ! DH2O = diffusion coefficient of H2O in air; m**2/s
     ! OUTPUT:
     ! [diff_HNO3]=m**2/s
     ! Reif, 1965
     
     IMPLICIT NONE
     
     REAL(dp), INTENT(in) :: DH2O
     REAL(dp) :: diff_HNO3
     
     diff_HNO3=DH2O*sqrt(MolMassH2O/MolMassHNO3)
  
  END FUNCTION diff_HNO3
  
  ELEMENTAL FUNCTION speed_HNO3(TEMP)
      ! O.Kirner 06/11/08
      ! Function to calculate the mean molevular speed of HNO3
      ! INPUT:
      ! TEMP = temperature; K
      ! OUTPUT:
      ! [speed_HNO3]=m/s
      
      IMPLICIT NONE
      
      REAL(dp), INTENT(in) :: TEMP
      REAL(dp) :: speed_HNO3
      
      speed_HNO3=sqrt((8.0_dp*Rgas*TEMP)/(pi*MolMassHNO3*1.0E-3_dp))
      
  END FUNCTION speed_HNO3
  
  ELEMENTAL FUNCTION press_HNO3(HNO3_gl,TEMP,PRESS)
      ! O.Kirner 06/11/08
      ! Function to calculate the ambient HNO3 partial Pressure
      ! INPUT:
      ! HNO3_gl=gas phase of HNO3; mol/mol 
      ! TEMP = temperature; K
      ! PRESS = pressure; hPa
      ! OUTPUT:
      ! [press_HNO3]=Pa
      
      IMPLICIT NONE
      
      REAL(dp), INTENT(in) :: HNO3_gl, TEMP, PRESS
      REAL(dp) :: factor, press_HNO3
      
      factor   = mix2conc(PRESS, TEMP)*1e6 ! convert in 1/m**2
      press_HNO3=HNO3_gl*factor*TEMP*Rgas/N_Avogadro
  END FUNCTION press_HNO3
  
  ELEMENTAL FUNCTION press_H2O(H2O_g,TEMP,PRESS)
      ! O.Kirner 06/11/08
      ! Function to calculate the ambient H2O partial Pressure
      ! INPUT:
      ! HNO3_gl=gas phase of H2O; mol/mol 
      ! TEMP = temperature; K
      ! PRESS = pressure; hPa
      ! OUTPUT:
      ! [press_H2O]=Pa
      
      IMPLICIT NONE
      
      REAL(dp), INTENT(in) :: H2O_g, TEMP, PRESS
      REAL(dp) :: factor, press_H2O
      
      factor   = mix2conc(PRESS, TEMP)*1e6 ! convert in 1/m**2
      press_H2O=H2O_g*factor*TEMP*Rgas/N_Avogadro
  END FUNCTION press_H2O
  
  ELEMENTAL FUNCTION press_HNO3_over_NAT(TEMP,pH2O)
      ! O.Kirner 06/11/08
      ! Function to calculate the vapor pressure of HNO3 over NAT
      ! INPUT:
      ! TEMP = temperature; K
      ! pH2O = ambient H2O partial pressure; Pa
      ! OUTPUT:
      ! [press_HNO3_over_NAT]=Pa
      ! Hanson and Mauersberger, 1988
      
      IMPLICIT NONE
      
      REAL(dp), INTENT(in) :: TEMP, pH2O
      REAL(dp) :: press_HNO3_over_NAT, DUM1, DUM2
      
      INTRINSIC LOG10
      
      DUM1=-2.7836_dp-0.00088_dp*TEMP
      DUM2=38.9855_dp-11397.0_dp/TEMP+0.009179*TEMP
      
      press_HNO3_over_NAT=10.0_dp**(DUM1*LOG10((pH2O/1.01325E5_dp)*760.0_dp)+DUM2) &
                          /760.0_dp*1.01325E5_dp
  
  END FUNCTION press_HNO3_over_NAT
  
  ELEMENTAL FUNCTION press_H2OS(H2O_g,TEMP,PRESS)
      ! O.Kirner 06/11/08
      ! Function to calculate saturation vapor pressure of H2O
      ! INPUT:
      ! H2O_g = gas phase of H2O; mol/mol
      ! TEMP = temperature; K
      ! PRESS = pressure; hPa
      ! OUTPUT:
      ! [press_H2OS]=Pa
      ! Hanson and Mauersberger, 1988
      
      IMPLICIT NONE
      
      REAL(dp), INTENT(in) :: PRESS, TEMP, H2O_g
      REAL(dp) :: factor, press_H2OS
      
      factor   = mix2conc(PRESS, TEMP) ! convert in 1/cm**2
  
      press_H2OS = H2O_g*factor*TEMP/ctoa
  
  END FUNCTION press_H2OS
  
  ELEMENTAL FUNCTION press_HNO3S(pH2OS,TEMP)
      ! O.Kirner 06/11/08
      ! Function to calculate the saturation vapor pressure
      ! INPUT:
      ! TEMP = temperature; K
      ! pH2O = saturation vapor pressure of H2O; Pa
      ! OUTPUT:
      ! [press_HNO3S]=Pa
      ! Hanson and Mauersberger, 1988
      
      IMPLICIT NONE
      
      REAL(dp), INTENT(in) :: TEMP, pH2OS
      REAL(dp) :: press_HNO3S, DUM1, DUM2
      
      INTRINSIC LOG10
      
      DUM1=-2.7836_dp-0.00088_dp*TEMP
      DUM2=38.9855_dp-11397.0_dp/TEMP+0.009179*TEMP
      
      press_HNO3S=10.0_dp**(DUM1*LOG10(pH2OS*760.0_dp)+DUM2)/760.0_dp
  
  END FUNCTION press_HNO3S  
  
  ELEMENTAL FUNCTION svapn_fkt(HNO3_gl,pH2OS,PRESS,TEMP)
      ! O.Kirner 06/11/08
      ! Function to calculate svapn
      ! INPUT:
      ! TEMP = gas phase of HNO3; mol/mol
      ! pHNO3S = saturation vapor pressure of HNO3; Pa
      ! PRESS = pressure; hPa
      ! TEMP = temperature; K
      ! TEMP_svapn = temperature include supercooling
      
      IMPLICIT NONE
      
      REAL(dp), INTENT(in) :: HNO3_gl, pH2OS, PRESS, TEMP
      REAL(dp) :: factor, svapn_fkt, TEMP_svapn, pHNO3S_svapn
            
      TEMP_svapn=TEMP-NatFormThreshold ! supercooling 
          
      factor   = mix2conc(PRESS, TEMP_svapn) ! convert in 1/cm**2
      
      pHNO3S_svapn=press_HNO3S(pH2OS,TEMP_svapn)
        
      svapn_fkt=HNO3_gl*factor*TEMP_svapn/(ctoa*pHNO3S_svapn)
  
  END FUNCTION svapn_fkt
  
  ELEMENTAL FUNCTION  growth_factor(DHNO3, vhno3, nat_rad, pHNO3, pHNO3_over_NAT, TEMP)
     ! O.Kirner 06/11/09
      ! Function to calculate the growth_factor for NAT-particles
      ! INPUT:
      ! DHNO3 = diffusion coefficient of HNO3 in air in m**2/s
      ! vhno3 = molecular speed of HNO3 in m/s
      ! nat_rad = Radius of particles in m
      ! pHNO3 = HNO3 partial Pressure in Pa
      ! pNAT_HNO3 = vapor pressure of HNO3 over NAT in Pa
      ! TEMP = Temperature in K
      ! OUTPUT:
      ! [growth_factor]=m**2/s
      ! Carslaw et al. ,2002
      
      IMPLICIT NONE
      
      REAL(dp), INTENT(in) :: DHNO3, vhno3, nat_rad, pHNO3, pHNO3_over_NAT, TEMP
      REAL(dp) :: growth_factor 
      REAL(dp) :: D_star_HNO3 
      
      D_star_HNO3 = DHNO3/(1.0_dp+4.0_dp*DHNO3/(vhno3*nat_rad))
         !diff. coeff. of HNO3 in air, [m**2/s] accounting for mass transfer continuum
         !effect (when size ~ mean free path)
        
      growth_factor = D_star_HNO3*MolMassHNO3/(rho_NAT*TEMP*Rgas)*(pHNO3-pHNO3_over_NAT)     
  
  END FUNCTION growth_factor
!=========================================================================  

!=========================================================================
ELEMENTAL FUNCTION mz_NAT_vel(TimeStep,PRESS, TEMP, r_NAT, G, N_NAT)
  !------------
  ! DESCRIPTION:
  !------------
  ! O.Kirner 06/12/06
  ! Carslaw et. al, 2002
  ! -----------
  ! INPUT:
  ! -----------
  ! TimeStep       = duration of time step in s
  ! PRESS          = air pressure within grid box in hPa
  ! TEMP           = air temperature in K
  ! r_NAT          = mean radius of NAT particle in current sizebin in m
  ! N_NAT        = number density of particle in current size bin in 1/m**3
  !  -----------
  ! Output:
  ! -----------
  ! mz_NAT_vel     = sedimentation velocity of NAT in m/s
  ! -----------
  ! PARAMETERS:
  ! -----------
  ! r_NAT          = radius of the NAT molecules
  ! mfp_HNO3       = mean free path of one HNO3 molecule in m
  ! visc_air       = viscosity of air in g/(ms)
  ! Cc             = Cunninham slip flow correction factor, dimensionless
  ! SedFactor      = Sedimentation factor in 1/(ms)
  ! G              = growth factor in m**2/s
  
  IMPLICIT NONE
  
  REAL(dp), INTENT(in) :: PRESS, TEMP
  REAL(dp), INTENT(in) :: TimeStep, r_NAT 
  REAL(dp), INTENT(in) :: G
  REAL(dp), INTENT(in) :: N_NAT

  REAL(dp), parameter :: r_mol_HNO3=3.0e-10_dp 
  REAL(dp) :: mfp_HNO3, visc_air
  REAL(dp) :: Cc, SedFactor 
  REAL(dp) :: mz_NAT_vel 
  
  INTRINSIC MAX
  
  mfp_HNO3=mean_free_path(PRESS,TEMP,r_mol_HNO3)
  visc_air=viscosity_air(TEMP)
  
  Cc=1.0_dp+mfp_HNO3/r_NAT*(1.257_dp+0.4_dp*exp(-1.1*r_NAT/mfp_HNO3))
  SedFactor=2.0_dp*g_acc*rho_NAT*Cc/(9.0_dp*visc_air)
  IF (N_NAT> 0.1_dp) THEN   
        mz_NAT_vel=max(SedFactor*(r_NAT**2.0_dp+G*TimeStep),0.0_dp)
  ELSE 
        mz_NAT_vel=0.0_dp
  END IF

END FUNCTION mz_NAT_vel

ELEMENTAL FUNCTION mean_free_path(PRESS,TEMP,r_mol)
  !-----------
  ! DESCRIPTION
  !-----------
  ! O.Kirner 06/12/07
  ! This function calculates the mean free path dependent on r_mol
  ! Reif, 1965
  ! ---------
  ! INPUT:
  ! ---------
  ! PRESS          = air pressure within grid box in hPa
  ! TEMP           = air temperature in K
  ! r_mol          = radius of one molecule in m
  ! ---------
  ! Output:
  ! ---------
  ! mean_free_path = mean free path of a molecule in m
  ! --------- 
  ! Other variables:
  ! --------
  ! sigm_mol     = cross section of 1 air molecules in m**2 
  ! air_density  = in molecules/cm3
  IMPLICIT NONE
  
  REAL(dp), INTENT(in) :: r_mol   
  REAL(dp), INTENT(in) :: PRESS, TEMP 
  REAL(dp) ::  air_density, mean_free_path, factor, sigm_mol
  
  factor=mix2conc(PRESS, TEMP)
  air_density=1.0_dp*factor
  
  sigm_mol=pi*r_mol**2.0_dp
  
  mean_free_path=1.0_dp/(sqrt(2.0_dp)*air_density*1.0e6_dp*sigm_mol)

END FUNCTION mean_free_path
  
ELEMENTAL FUNCTION viscosity_air(TEMP)
  !-----------
  ! DESCRIPTION
  !-----------
  ! O.Kirner 06/12/07
  ! This function calculates the viscosity of air
  ! ---------
  ! INPUT:
  ! ---------
  ! TEMP           = air temperature in K
  ! ---------
  ! Output:
  ! ---------
  ! viscosity_air in g/ms
  IMPLICIT NONE
      
  REAL(dp), INTENT(in) :: TEMP
  REAL(dp) :: viscosity_air
  
  IF ((TEMP-273.15) .ge. 0.0_dp) THEN
      viscosity_air=(1.718_dp+0.0049_dp*(TEMP-273.15_dp))*1.0e-2_dp
  ELSE
      viscosity_air=(1.718_dp+0.0049_dp*(TEMP-273.15_dp)-1.2e-5_dp*(TEMP-273.15_dp)**2.0_dp)*1.0e-2_dp 
  END IF

END FUNCTION viscosity_air
!=========================================================================
! ka_ok_20100115-

END MODULE messy_msbm
