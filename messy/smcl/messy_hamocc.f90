MODULE messy_hamocc
  !----------------------------------------------------------------------------
  ! MESSY ocean chemistry. This code is based on the HAMOCC code
  !
  ! Please contact Maier-Reimer (MPIM) for any questions and 
  ! aknowledgement.
  !
  ! A. Pozzer - 22-01-2008
  ! Updates to hamocc version used in mpiesm-1.0.00 - B.Kern 26 Jun 2012
  ! Iron in pore water not included! (yet)
  ! --> maybe problems when reading fields from netcdf file ?? (messy_hamocc_e5)
  ! mz_bk_20120707 detritus now advected tracer!
  !
  !   [M]  = [mol/L] = [kmol/m^3]
  !   [nM] = [ 10^-9 mol/L] = [10^-9 mol/dm^3] = [10^-6 mol/m^3]
  !        = [10^-9 kmol/m^3] 
  !   [M]  = [mol/dm^3] = [mol/L] =  [kmol/m^3]
  !----------------------------------------------------------------------------

  ! MESSy
  USE messy_main_constants_mem,  ONLY: DP, WP, SP, STRLEN_MEDIUM, STRLEN_LONG

  IMPLICIT NONE 

  INTRINSIC MAX, SIGN
  
  PUBLIC

  SAVE ! mz_bk_20120613

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'hamocc'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.2'

  ! CTRL-NAMELIST PARAMETER
  LOGICAL, PUBLIC :: L_CHECK 
  LOGICAL, PUBLIC :: L_AGG
  ! mz_bk_20120606+
  LOGICAL, PUBLIC :: L_BGC_DIAG
  ! mz_bk_20120606-
  INTEGER  :: isac          ! acceleration factor for sediment, from namelist
  REAL(DP) :: rmasks = 0.0  ! value at wet cells in sediment.
  CHARACTER(LEN=200), PUBLIC :: hamocc_ini_files_path = '' ! path of initialised
  !                                                                        field

  TYPE PTR_3D_ARRAY
     REAL(DP), DIMENSION(:,:,:), POINTER :: PTR
  END TYPE PTR_3D_ARRAY
  TYPE PTR_2D_ARRAY
     REAL(DP), DIMENSION(:,:), POINTER :: PTR
  END TYPE PTR_2D_ARRAY

  INTEGER  :: n90depth
  !____________________________________________________________________________!
  !
  !   PARAMETERS : MO_PARAM1_BGC
  !____________________________________________________________________________!

  INTEGER, PARAMETER :: ks = 12, ksp = ks + 1                ! sediment layers
  INTEGER, PARAMETER :: kwrbioz=8                            ! euphotic layers

  ! advected tracers
  INTEGER, PARAMETER :: i_base_adv = 15,                                      &
       &                isco212    =  1,                                      &
       &                ialkali    =  2,                                      &
       &                iphosph    =  3,                                      &
       &                ioxygen    =  4,                                      &
       &                igasnit    =  5,                                      &
       &                iano3      =  6,                                      &
       &                isilica    =  7,                                      &
       &                idoc       =  8,                                      &
       &                iphy       =  9,                                      &
       &                izoo       = 10,                                      &
       &                ian2o      = 11,                                      &
       &                idms       = 12,                                      &
       &                iiron      = 13,                                      &
       &                ifdust     = 14,                                      &
  ! mz_bk_20120706+
       &                idet       = 15
  ! mz_bk_20120706-

  INTEGER, PARAMETER :: i_iso_adv = 0

  INTEGER, PARAMETER :: i_cfc_adv= 0

  INTEGER, PARAMETER ::                                                       &
       ! AGG always active here. In case the pointer are not associated
       ! (i.e. fields eq 0)
!!$!!#ifdef AGG
       &                i_agg_adv = 2,                                        &
       &                inos      = i_base_adv + i_iso_adv + i_cfc_adv + 1,   &
       &                iadust    = i_base_adv + i_iso_adv + i_cfc_adv + 2
!!$!!#else 
  !                      i_agg_adv= 0
!!$!!#endif

  ! total number of advected tracers
  INTEGER, PARAMETER :: ntraad = i_base_adv + i_iso_adv + i_cfc_adv + i_agg_adv

  ! non-advected (fast sinking) tracers
  ! mz_bk_20120706+
  ! INTEGER, PARAMETER :: idet     = ntraad + 1,                              &
  !      &                icalc    = ntraad + 2,                              &
  !      &                iopal    = ntraad + 3,                              &
  !      &                i_base   = 3
  INTEGER, PARAMETER :: icalc    = ntraad + 1,                                &
       &                iopal    = ntraad + 2,                                &
       &                i_base   = 2
  ! mz_bk_20120706-
                            
  INTEGER, PARAMETER :: i_iso    = 0

  INTEGER, PARAMETER :: nocetra  = ntraad + i_base + i_iso

  ! sediment
  INTEGER, PARAMETER :: issso12  = 1,                                         &
       &                isssc12  = 2,                                         &
       &                issssil  = 3,                                         &
       &                issster  = 4,                                         &
       &                nsss_base= 4
  INTEGER, PARAMETER :: nsss_iso = 0 !for isotope
  INTEGER, PARAMETER :: nsedtra  = nsss_base + nsss_iso

  ! pore water tracers
  ! index must be the same as for ocetra otherwise problems in dipowa.f90!
  INTEGER, PARAMETER :: npowa_base = 7,                                       &
       &                ipowaic    = 1,                                       &
       &                ipowaal    = 2,                                       &
       &                ipowaph    = 3,                                       &
       &                ipowaox    = 4,                                       &
       &                ipown2     = 5,                                       &
       &                ipowno3    = 6,                                       &
       &                ipowasi    = 7
  INTEGER, PARAMETER :: npowa_iso  = 0 ! for isotope
  INTEGER, PARAMETER :: npowtra    = npowa_base + npowa_iso

  ! mz_bk_20120606+
  ! Parameter for BGC diagnostics
  INTEGER, PARAMETER :: ibgcprod = 12,                                        &
       &                kphosy   =  1,                                        &
       &                kgraz    =  2,                                        &
       &                kexport  =  3,                                        &
       &                kdelcar  =  4,                                        &
       &                kdelsil  =  5,                                        &
       &                kdmsprod =  6,                                        &
       &                kdms_bac =  7,                                        &
       &                kdms_uv  =  8,                                        &
       ! mz_bk_20120813+
       ! &                krain    =  9,                                      &
       ! mz_bk_20120813-
       &                kflim    = 10,                                        &
       &                kplim    = 11,                                        &
       &                knlim    = 12
  ! mz_bk_20120606-

  ! BGC spinup
  ! not implemented in MESSy system ! mz_bk_20120626
  LOGICAL :: lspinbgc = .false.

  ! Logical unit number for I/O.
  INTEGER :: io_stdo_bgc = 6       !  standard output stream

  ! Control variables
  REAL(DP) :: dtbgc                !  time step length [sec].
  REAL(DP) :: dtb                  !  time step length [days].


  !____________________________________________________________________________!
  !
  !   VARIABLES FOR MARINE BIOLOGY
  !____________________________________________________________________________!

  INTEGER,  DIMENSION (:,:), POINTER :: kbo
  REAL(DP), DIMENSION (:,:), POINTER :: bolay

  REAL(DP), DIMENSION (:,:), POINTER :: strahl
  
  REAL(DP), DIMENSION (:,:), POINTER :: alar1max
  REAL(DP), DIMENSION (:,:), POINTER :: TSFmax
  REAL(DP), DIMENSION (:,:), POINTER :: TMFmax

  REAL(DP), DIMENSION (:,:,:), POINTER :: abs_oce


  ! mz_bk_20120724+
  ! Pre-define biogeochemical variables (beleg.f90), can be overwritten by
  ! namelist values
  ! parameters depending on time are in d^-1 and multiplied by dtb in
  ! hamocc_beleg to get time-step^-1
  !
  ! 1) Redfield ratios and constants
  REAL(DP)::                                                                  &
       ! extended redfield ratio declaration
       ! Note: stoichiometric ratios are based on Takahashi etal. (1985)
       ! P:N:C:-O2 + 1:16:122:172
       & ro2ut    = 172._dp,      & ! -O2:P ratio
       & rcar     = 122._dp,      & ! C:P ratio
       & rnit     = 16._dp,       & ! N:P ratio
       & nitdem   = 121.6_dp,     & ! nitrate demand to remin. 1 mol P
       !                            ! in suboxic water
       & n2prod   = 68.8_dp,      & ! N2 production for 1 mol P remineralized
       !                            ! in suboxic water
       & rcalc    = 35._dp,       & ! calcium carbonate to organic phosphorous
       !                            ! production ratio ! iris 40!
       & ropal    = 25._dp,       & ! opal to organic phosphorous production
       !                            ! ratio !iris 25 !
       & n2_fixation = 0.005_dp,  & ! nitrogen fixation by blue algae rate
       !                            ! in 1/d (hamocc_cyano)
       ! & o2ut = 172._dp,        & ! -O2:P ratio in pore water ! never used!!
       & rno3     =  16._dp,      & ! N:P ratio in pore water in sediment
       & perc_diron = 0.035_dp * 0.01_dp / 55.85_dp,& ! weight percent iron in
       !                            ! dust deposition (0.035) times Fe
       !                            ! solubility (0.01) / 55.85 g--> mol
       ! the three values below are from Johnson et al., 1997
       ! Mar. Chemistry 57, 137-161
       !  riron = 5._dp * rcar * 1.e-6_dp ! 'Q'=5.*122.*1.e-6 = 6.1e-4
       !  !              ! (iron to phosphate stoichiometric ratio * umol->mol)
       !  riron = 2._dp * rcar * 1.e-6_dp ! 14.2.06: 5.->2. (2.5umol P/nmol Fe)
       !  !                               ! emr: 3.*rcar*1.e-6
       !& riron = 3._dp * rcar * 1.e-6_dp,& ! 06.3.06: 2.->3. coex90 = 12.2GTC/a
       !!                                   !          before change
       & riron = 3._dp*122._dp*1.e-6_dp,& ! 06.3.06: 2.->3. coex90 = 12.2GTC/a
       !                                   !          before change
       & fesoly  = 0.6_dp * 1.e-9_dp,&     ! global mean/max (?) dissolved iron
       !                ! concentration in deep water (for relaxation) [mol/m3]
       !              ! 'fesoly' stands for 'apparent solubility value' of iron
       & relaxfe = 0.05_dp / 365._dp ! relaxation time for iron to fesoly
       !                             ! in 1/d; corresponds to 20 yrs
       !  relaxfe = 0.005_dp / 365._dp * dtb  ! changed as suggested by emr
       !  !                                   ! (should probably be named
       !  !                                   ! scavenge_fe or remove_fe)
       ! (relaxfe is only used to remove iron, corresponding to Johnson's 200yr
       ! residence time)
       ! but 200 yrs might be too long, as iron concentrations in the Atlantic
       ! are too high.
  !
  ! 2) Radiation and Pothosythesis
  REAL(DP)::                                                                  &
       & pi_alpha = 0.02_dp,      & ! initial slope of production vs irradiance
       !                            ! curve (alpha) (0.002 for 10 steps per day)
       & fPAR     = 0.4_dp          ! fraction of Photosynthetic Active
  !                                 ! Radiation
#if defined MPIOM_13B
  REAL(DP)::                                                                  &
       ! parameters for sw-radiation attenuation
       ! Analogue to Moore et al., Deep-Sea Research II 49 (2002), 403-462
       ! 1 kmolP = (122*12/60)*10^6 mg[Chlorophyll] 
       ! js: in Moore et al. atten_w and atten_c are mean values 
       ! for the mixed layer. (the values just appear in their appendix,
       ! without further reference. can they be used here?)
       ! here we might need smaller values, e.g. scale by layer depth / mld
       & ctochl   = 60._dp,       & ! C to Chlorophyll ratio ! HadOCC: 40. 
       & atten_w  = 0.04_dp,      & ! Gelbstoff attenuation in 1/m
       & atten_f  = 0.4_dp,       & ! fraction of sw-radiation directly
       !                            ! absorbed in surface layer (only used
       !                            ! if FB_BGC_OCE) [feedback bgc-ocean]
       & atten_c                    ! phytoplankton attenuation in 1/m
  !                                 ! "self-shading" 7.32e5 mg Chl/m3;
  !                                 ! not set here, set in hamocc_beleg to:
  !                                 ! 0 .03 * rcar * (12. / ctochl) * 1.e6
#endif
  !
  ! 3) Phytoplankton
  REAL(DP)::                                                                  &
       & phytomi  = 1.e-11_dp,    & ! minimum concentration of phytoplankton
       !                            ! i.e. 1.e-11 kmol P/m3 = 1e-5 mmol P/m3 
       & bkphy    = 1.e-8_dp,     & ! half saturation constant for
       !                            ! phytoplankton (kmol/m3)
       & bkopal   = 1.e-6_dp,     & ! half saturation constant for Si(OH)4
       !                            ! i.e. 1.0  mmol Si/m3
       & remido   = 0.008_dp,     & ! 1/d -remineralization rate (of DOM)
       !                            ! KS, JS, EMR 12/07/2007
       & dyphy    = 0.008_dp,     & ! 1/d -mortality rate of phytoplankton 
       & gammap   = 0.03_dp         ! 1/d -exudation rate
  !
  ! 4) Zooplankton
  REAL(DP)::                                                                  &
       & bkzoo    = 4.e-8_dp,     & ! half saturation constant for grazing
       !                            ! i.e. 0.04 mmol P/m3
       & grami    = 1.e-11_dp,    & ! minimum concentration of zooplankton
       !                            ! i.e. 1.e-11 kmol P/m3 = 1e-5 mmol P/m3 
       & zinges   = 0.6_dp,       & ! dimensionless fraction -
       !                            ! assimilation efficiency of zooplankton
       & epsher   = 0.8_dp,       & ! dimensionless fraction - 
       !                            ! (1-epsher)=fraction of grazing egested
       & grazra   = 1.0_dp,       & ! 1/d -grazing rate [emr: 0.6-0.9]
       & spemor   = 3.e6_dp,      & ! 1/d -mortality rate of zooplankton
       & gammaz   = 0.06_dp,      & ! 1/d -excretion rate
       & ecan     = 0.95_dp         ! fraction of mortality as PO_4
  !
  ! 5) Deep Ocean Remineralisation, Dissolution, Mortality, Sinking
  REAL(DP)::                                                                  &
       ! note: this sinking speeds have no effect in case of L_AGG=.TRUE.
       & sinkspeed_poc  = 5.0_dp, & ! daily sinking speed of poc
       & sinkspeed_opal = 30.0_dp,& ! daily sinking speed of opal
       & sinkspeed_cal  = 30.0_dp,& ! daily sinking speed of cal
       ! water column remineralisation constants
       & drempoc  = 0.025_dp,     & ! detritus remineralisation rate in 1/d
       !                            ! 0.75/month.
       !                            ! H3.1: 0.2/month k=1, 0.1/month k>1
       & dremdoc  = 0.004_dp,     & ! DOC remineralisation rate in 1/d
       & dremn2o  = 0.01_dp,      & ! remineralisation using N2O in 1/d
       & denitrification = 0.05_dp,&! fraction of remineralisation from
       !                            ! denitrification in oxygen low waters
       & sulfate_reduction  = 0.005_dp,& ! remineralisation from sulfate
       !                            ! reduction in oxygen low waters in 1/d
       & dremopal = 0.01_dp,      & ! Opal dissolution rate in 1/d
       !                            ! 0.01 -> 0.001 js10072006 : 
       !                            ! slightly overdone -->0.0075
       & dremcalc = 0.075_dp,     & !1/d ! 0.2 -> 0.02 js10072006 : 
       !                            ! slightly overdone --> 0.075
       & dphymor  = 0.1_dp,       & ! Phytoplankton mortality rate in 1/d
       & dzoomor  = 0.02_dp         ! Zooplankton mortality rate in 1/d
  !
  ! 7) DMS
       ! Set constants for calculation of DMS ( mo_carbch )
       ! Parameters are a result from Kettle optimisation
       ! (fit to Kettle data, version emr, not tinka)
  REAL(DP), DIMENSION(6)::                                                    &
       & dmspar = (/ 10._dp,  & !temperature dependent release by phytoplankton
       & 0.0075_dp,               & ! photolysis (uv-destruction)
       & 0.0096_dp,               & ! microbial consumption
       ! dmspar(4)=1.25_dp*0.107638502E+00_dp ! production with delcar
       ! dmspar(5)=1.25_dp*0.109784522E-01_dp ! production with delsil (diatoms)
       & 0.5_dp*0.107638502E+00_dp,&! production with delcar
       & 0.5_dp*0.109784522E-01_dp,&! production with delsil (diatoms)
       & 0.100000000E-07_dp /)      ! half saturation microbial
  !
  ! 8) AGG (Aggregation sinking)
  REAL(DP)::                                                                  &
       & calmax   = 0.20_dp         ! max fraction of CaCO3 production (AGG)
  !                                 ! (max rain ratio)
  ! further AGG parameter defined below...


  ! values are set in hamocc_initialize after call to hamocc_read_nml_ctrl
  ! in messy_hamocc_e5.f90 
  REAL(DP) :: &
         remido_dtb,   &
         dyphy_dtb,    &
         grazra_dtb,   &
         spemor_dtb,   &
         gammap_dtb,   &
         gammaz_dtb,   &
         drempoc_dtb,  &
         dremdoc_dtb,  &
         dremopal_dtb, &
         dremcalc_dtb, &
         dremn2o_dtb,  &
         dphymor_dtb,  &
         dzoomor_dtb,  &
         relaxfe_dtb  
    ! mz_bk_20120724-


  REAL(DP) :: rrrcl
  REAL(DP) ::wdust  
  ! mz_bk_20120725+
  REAL(DP) :: rnoi
  ! mz_bk_20120725-


#if defined MPIOM_2000
  ! parameters for sw-radiation fraction
  ! analogous to Zielinski et al., Deep-Sea Research II 49 (2002), 3529-3542

  REAL(DP), PARAMETER :: redfrac = 0.4_dp ! red fraction of the spectral domain
  !                                       ! (> 580nm)

  REAL(DP), PARAMETER :: c_to_chl = 12.0_dp / 60.0_dp ! Carbon to Chlorophyll
  !                                                   ! ratio
  REAL(DP), PARAMETER :: r_car = 122.0_dp             ! Redfield ratio
  REAL(DP), PARAMETER :: pho_to_chl = r_car * c_to_chl * 1.e6_dp
  !                                 ! 1 kmolP = (122*12/60)*10^6 mg[Chlorophyll]

  REAL(DP), PARAMETER :: atten_r = 0.35_dp ! attenuation of red light [m-1]
  REAL(DP), PARAMETER :: atten_w = 0.03_dp ! attenuation of blue/green light
  !                               ! in clear water between 400nm and 580nm [m-1]
  REAL(DP), PARAMETER :: atten_c = 0.04_dp ! attenuation of blue/green light
  !                                        ! by chlorophyll [m-1]
  ! mz_bk_20120724+
  ! REAL(DP):: ctochl = c_to_chl
  REAL(DP):: ctochl = 12._dp / c_to_chl
  ! mz_bk_20120724-
      
#endif

!   REAL(DP):: phytomi, grami, grazra, rrrcl
!   REAL(DP):: remido, dyphy, zinges, epsher, spemor, gammap, gammaz, ecan
!   ! mz_bk_20120627+
!   REAL(DP):: ro2ut, rcar, rnit, rnoi, rnit23, rnit13, rcalc, ropal
!   ! REAL(DP):: ro2ut, rcar, rnit, rnoi, rnit23, rnit13, rcalc, ropal, bluefix
!   ! mz_bk_20120627-
!   REAL(DP):: bkphy, bkzoo, bkopal, bifr13, bifr14, plafr13, plafr14
!   ! mz_bk_20120626+
!   ! REAL(DP):: wpoc, wcal, wopal, drempoc, dremdoc, dremcalc, dremn2o
!   REAL(DP):: drempoc, dremdoc, dremcalc, dremn2o

!   REAL(DP) :: n2_fixation
!   REAL(DP) :: sulfate_reduction
!   REAL(DP) :: denitrification

!   REAL(DP) :: sinkspeed_poc  ! daily sinking speed of poc (namelist parameter)
!   REAL(DP) :: sinkspeed_opal ! daily sinking speed of opal (namelist parameter)
!   REAL(DP) :: sinkspeed_cal  ! daily sinking speed of cal (namelist parameter)
!   ! mz_bk_20120626-
!   REAL(DP):: dphymor, dzoomor, dremopal, calmax, gutc
!   REAL(DP):: psedi, csedi, ssedi
!   REAL(DP):: perc_diron, riron, fesoly, relaxfe, wdust  
!   ! mz_bk_20120627+
!   REAL(wp) :: nitdem,n2prod
!   ! mz_bk_20120627-

!   ! RADIATION COUPLING
! #if defined MPIOM_13B
!   REAL(DP) :: atten_w, atten_c, atten_f
!   REAL(DP):: ctochl 
! #elif defined MPIOM_2000
!   ! parameters for sw-radiation fraction
!   ! analogous to Zielinski et al., Deep-Sea Research II 49 (2002), 3529-3542

!   REAL(DP), PARAMETER :: redfrac = 0.4_dp ! red fraction of the spectral domain
!   !                                       ! (> 580nm)

!   REAL(DP), PARAMETER :: c_to_chl = 12.0_dp / 60.0_dp ! Carbon to Chlorophyll
!   !                                                   ! ratio
!   REAL(DP), PARAMETER :: r_car = 122.0_dp             ! Redfield ratio
!   REAL(DP), PARAMETER :: pho_to_chl = r_car * c_to_chl * 1.e6_dp
!   !                                 ! 1 kmolP = (122*12/60)*10^6 mg[Chlorophyll]

!   REAL(DP), PARAMETER :: atten_r = 0.35_dp ! attenuation of red light [m-1]
!   REAL(DP), PARAMETER :: atten_w = 0.03_dp ! attenuation of blue/green light
!   !                               ! in clear water between 400nm and 580nm [m-1]
!   REAL(DP), PARAMETER :: atten_c = 0.04_dp ! attenuation of blue/green light
!   !                                        ! by chlorophyll [m-1]
!   ! mz_bk_20120724+
!   ! REAL(DP):: ctochl = c_to_chl
!   REAL(DP):: ctochl = 12._dp / c_to_chl
!   ! mz_bk_20120724-
        
! #endif

!   REAL(DP):: pi_alpha
!   ! mz_bk_20120515+
!   ! included fPAR fraction of radiation available for biological production
!   ! from HAMOCC mpiesm-1.0.00
!   REAL(DP):: fPAR
!   ! mz_bk_20120515-
  ! mz_bk_20120724-

!!$!! #ifdef AGG      
  REAL(DP):: SinkExp, FractDim, Stick, cellmass, cellsink, fsh, fse
  REAL(DP):: alow1, alow2, alow3, alar1, alar2, alar3, TSFac, TMFac
  REAL(DP):: vsmall, safe, pupper, plower, zdis
  REAL(DP):: dustd1, dustd2, dustd3, dustsink
  REAL(DP):: snow
!!$!! #endif

  !____________________________________________________________________________!
  !
  !   VARIABLES FOR SEDIMENT MODULES
  !____________________________________________________________________________!

!!$!!      REAL(DP), DIMENSION (:,:,:,:), POINTER :: sedlay
!!$!!      REAL(DP), DIMENSION (:,:,:,:), ALLOCATABLE :: powtra
!!$!!      REAL(DP), DIMENSION (:,:,:), ALLOCATABLE :: burial
!!$!!      REAL(DP), DIMENSION (:,:,:), ALLOCATABLE :: sedhpl

  TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER  ::  sedlay => NULL()
  TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER  ::  powtra => NULL()
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER  ::  burial => NULL()
  REAL(DP), DIMENSION (:,:,:), POINTER :: sedhpl => NULL()
  ! mz_bk_20120613+
  ! TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER, SAVE  ::  sedlay => NULL()
  ! TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER, SAVE  ::  powtra => NULL()
  ! TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  ::  burial => NULL()
  ! REAL(DP), DIMENSION (:,:,:), POINTER, SAVE :: sedhpl => NULL()
  ! mz_bk_20120613-

  REAL(DP), DIMENSION (:), POINTER :: seddw
  REAL(DP), DIMENSION (:), POINTER :: porsol
  REAL(DP), DIMENSION (:), POINTER :: porwah
  REAL(DP), DIMENSION (:), POINTER :: porwat

  REAL(DP), DIMENSION (:), POINTER :: dzs
  REAL(DP), DIMENSION (:), POINTER :: seddzi

  REAL(DP), DIMENSION (:,:), POINTER :: silpro
  REAL(DP), DIMENSION (:,:), POINTER :: prorca
  REAL(DP), DIMENSION (:,:), POINTER :: prcaca
  REAL(DP), DIMENSION (:,:), POINTER :: pror13
  REAL(DP), DIMENSION (:,:), POINTER :: prca13
  REAL(DP), DIMENSION (:,:), POINTER :: pror14
  REAL(DP), DIMENSION (:,:), POINTER :: prca14
  REAL(DP), DIMENSION (:,:), POINTER :: produs

  REAL(DP):: sedict, calcon, ansed, sedac, sedifti
  REAL(DP):: calcwei, opalwei, orgwei
  REAL(DP):: calcdens, opaldens, orgdens, claydens
  REAL(DP):: calfa, oplfa, orgfa, clafa, solfu

  !____________________________________________________________________________!
  !
  !   VARIABLES FOR INORGANIC CARBON CYCLE
  !____________________________________________________________________________!

!!$!!      REAL(DP), DIMENSION (:,:,:,:), POINTER :: ocetra
  TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER  ::  ocetra => NULL()
  TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER  ::  bgcprod => NULL()
  ! mz_bk_20120613+
  ! TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER, SAVE  ::  ocetra => NULL()
  ! TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER, SAVE  ::  bgcprod => NULL()
  ! mz_bk_20120613-
  REAL(DP), DIMENSION (:,:,:),   POINTER :: co3 => NULL()
  REAL(DP), DIMENSION (:,:,:),   POINTER :: hi  => NULL()
  REAL(DP), DIMENSION (:,:),   POINTER :: co2 => NULL()

  REAL(DP), DIMENSION (:,:,:), POINTER :: chemcm
  REAL(DP), DIMENSION (:,:,:), POINTER :: akw3
  REAL(DP), DIMENSION (:,:,:), POINTER :: akb3
  REAL(DP), DIMENSION (:,:,:), POINTER :: ak13
  REAL(DP), DIMENSION (:,:,:), POINTER :: ak23
  REAL(DP), DIMENSION (:,:,:), POINTER :: aksp
  REAL(DP), DIMENSION (:,:,:), POINTER :: satoxy
  REAL(DP), DIMENSION (:,:),   POINTER :: satn2o
!!$!!      REAL(DP), DIMENSION (:,:,:), ALLOCATABLE :: sedfluxo

  REAL(DP), DIMENSION (:,:), POINTER :: dusty ! in kg/m2/year
  !                                           ! (0 if ddpo(i,j,1) .gt. 0.5)
  REAL(DP), DIMENSION (:,:,:), POINTER :: c14pool

  !____________________________________________________________________________!

  PUBLIC :: dp, sp, STRLEN_MEDIUM

  PUBLIC :: hamocc_read_nml_ctrl
  PUBLIC :: hamocc_bodensed
  PUBLIC :: hamocc_beleg
  PUBLIC :: hamocc_init_var
  PUBLIC :: hamocc_chemcon
  PUBLIC :: hamocc_ocprod
  PUBLIC :: hamocc_cyano
  PUBLIC :: hamocc_carchm
  PUBLIC :: hamocc_powach
  PUBLIC :: hamocc_dipowa
  PUBLIC :: hamocc_powadi
  PUBLIC :: hamocc_sedshi

  PUBLIC :: hamocc_photoprod

CONTAINS

  !----------------------------------------------------------------------------
  SUBROUTINE hamocc_bodensed(ie,je,ke,ddpo)
    !
    ! correspondes to bodensed.f90
    ! set up of sediment layers
    !--------------------------------------------------------------------------

    INTRINSIC REAL

    INTEGER, INTENT(IN) :: ie, je, ke
    REAL(DP),INTENT(IN) :: ddpo(ie,je,ke)
    REAL(DP):: sumsed
    INTEGER :: i, j, k

!!$!!       dzs(1) = 0.001
!!$!!       dzs(2) = 0.003
!!$!!       dzs(3) = 0.005
!!$!!       dzs(4) = 0.007
!!$!!       dzs(5) = 0.009
!!$!!       dzs(6) = 0.011
!!$!!       dzs(7) = 0.013
!!$!!       dzs(8) = 0.015
!!$!!       dzs(9) = 0.017
!!$!!       dzs(10) = 0.019
!!$!!       dzs(11) = 0.021
!!$!!       dzs(12) = 0.023
!!$!!       dzs(13) = 0.025
!!$!!  
!!$!!       WRITE(*,*)  ' '
!!$!!       WRITE(*,*)  'Sediment layer thickness [m] : '
!!$!!       WRITE(*,'(5F9.3)') dzs
!!$!!       WRITE(*,*)  ' '

    porwat(1) = 0.85
    porwat(2) = 0.83
    porwat(3) = 0.8
    porwat(4) = 0.79
    porwat(5) = 0.77
    porwat(6) = 0.75
    porwat(7) = 0.73
    porwat(8) = 0.7
    porwat(9) = 0.68
    porwat(10) = 0.66
    porwat(11) = 0.64
    porwat(12) = 0.62

!      WRITE(io_stdo_bgc,*)  ' '
!      WRITE(io_stdo_bgc,*)  'Pore water in sediment: '
!      WRITE(io_stdo_bgc,'(5F9.3)') porwat
!      WRITE(io_stdo_bgc,*)  ' '

      
    seddzi(1)=500._dp
    DO k=1,ks
       porsol(k)=1._dp - porwat(k)
       IF(k .EQ. 1) porwah(k) = 0.5_dp * (1._dp + porwat(1))
       IF(k .GE. 2) porwah(k) = 0.5_dp * (porwat(k) + porwat(k-1))
       seddzi(k+1) = 1._dp / dzs(k+1)     ! inverse thickness of sediment layer
       seddw(k)    = 0.5_dp * (dzs(k) + dzs(k+1))
    ENDDO

    ! ******************************************************************
    ! the following section is to include the SEDIMENT ACCELERATION
    ! mechanism to accelerate the sediment:

    sedac = 1._dp / REAL(isac,dp)

    ! determine total solid sediment thickness sumsed
    ! and reduced sediment volume
    sumsed = 0._dp
    do k = 1, ks
       porwat(k) = porwat(k) * sedac
       porsol(k) = porsol(k) * sedac
       sumsed = sumsed + seddw(k)
    enddo

    ! determine reduced diffusion rate sedict,
    ! and scaling for bottom layer ventilation, sedifti

    sedict = 1.e-9_dp * dtbgc                           ! units? m/s2 ?
    sedifti = sedict / (sumsed**2)                      ! not used
    sedict=sedict*sedac

    !  WRITE(io_stdo_bgc,*)  'sediment acceleration factor: ',isac
    !  WRITE(io_stdo_bgc,*)  'sediment area reduction: ',sedac
    !  WRITE(io_stdo_bgc,*)  'new diffusion is: ',sedict
    !  WRITE(io_stdo_bgc,*)  'total sediment thickness [m]: ',sumsed
    !  WRITE(io_stdo_bgc,*)  'sedict / sumsed**2: ',sedifti
    !  WRITE(io_stdo_bgc,*)  'new pore water fraction: ',porwat


    ! ******************************************************************

    ! ******************************************************************
    ! densities etc. for SEDIMENT SHIFTING

    ! define (wei)ght of calcium carbonate, opal, and poc [kg/kmol]
    calcwei = 100._dp          ! 40+12+3*16 kg/kmol C
    opalwei =  60._dp          ! 28 + 2*16  kg/kmol Si
    orgwei  =  30._dp          ! from 12 kg/kmol * 2.5 POC[kg]/DW[kg]
    !                          ! after Alldredge, 1998: POC(g)/DW(g) 
    !                          ! = 0.4 of diatom marine snow, size 1mm3

    ! define densities of caco3, opal, caco3 [kg/m3]
    calcdens = 2600._dp
    opaldens = 2200._dp
    orgdens  = 1000._dp
    claydens = 2600._dp !quartz

    ! define volumes occupied by solid constituents [m3/kmol]
    calfa = calcwei / calcdens         ! 3.85e-2
    oplfa = opalwei / opaldens         ! 2.73e-2
    orgfa = orgwei / orgdens           ! 3.0e-2
    clafa = 1._dp / claydens           ! 3.85e-4 (clay is calculated in kg/m3)

    ! determine total solid column integrated sediment volume  (1-d)
    solfu = 0._dp
    DO k=1,ks
       solfu = solfu + seddw(k) * porsol(k)
    ENDDO

    ! ******************************************************************

    k = ke
    DO j=1,je
       DO i=1,ie
          kbo(i,j) = 1
          bolay(i,j) = 0._dp
          IF(ddpo(i,j,k) .GT. 0.5_dp) THEN
             bolay(i,j) = ddpo(i,j,k)               ! thickness of bottom layer
             kbo(i,j) = k
          ENDIF
       ENDDO
    ENDDO

    DO k=ke-1,1,-1
       DO j=1,je
          DO i=1,ie
             IF(ddpo(i,j,k).GT.0.5 .AND. ddpo(i,j,k+1).LE.0.5) THEN
                bolay(i,j)=ddpo(i,j,k)
                kbo(i,j)=k
             ENDIF
          ENDDO
       ENDDO
    ENDDO


  END SUBROUTINE hamocc_bodensed
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE hamocc_beleg(ie,je,ke,ddpo,zmini)
    !
    ! corresponds to mo_beleg_bgc.f90
    ! initialization of parameters for ocean BGC
    !--------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: ie, je, ke
    REAL(DP),INTENT(IN) :: ddpo(ie,je,ke)
    REAL(DP),INTENT(IN) :: zmini

    !LOCAL
    INTEGER :: i, j, k, l, m

!!$!! #ifdef AGG
    REAL(DP):: shear, talar1, snow, checksink
!!$!! #endif 
!!$!!#ifndef AGG
    REAL(DP):: dustd1, dustd2, dustsink
!!$!!#endif
    ! mz_bk_20120626+
    ! REAL(DP):: xpi, rad, radi, rmissing


    ! xpi       = 4. * ATAN(1.)
    ! rad       = xpi / 180.
    ! radi      = 1. / rad
    ! mz_bk_20120626-

    ! mz_bk_20120724+
!     !
!     ! Biology
!     !
!     !ik note that the unit is kmol/m3!
!     phytomi = 1.e-11_dp          ! i.e. 1e-5 mmol P/m3 minimum concentration of
!     !                            ! phytoplankton
!     grami   = 1.e-11_dp          ! i.e. 1e-5 mmol P/m3 minimum concentration of
!     !                            ! zooplankton | test e-11 ->e-8 very slow decay

!     !ik addded parameter definition; taken from OCPROD.F (dtb= 1./ steps/day)

!     ! mz_bk_20120515+
!     ! updated from HAMOCC mpiesm-1.0.00
!     remido = 0.008_dp * dtb      !KS, JS, EMR 12/07/2007
!     ! ! "old" value
!     ! remido=0.025*dtb          !1/d -remineralization rate (of DOM)
!     ! mz_bk_20120515-

!     ! mz_bk_20120718+
!     ! testing dyphy = 0.007_dp * dtb
!     dyphy = 0.008_dp * dtb       !1/d -mortality rate of phytoplankton 
!     ! dyphy = 0.007_dp * dtb       !1/d -mortality rate of phytoplankton 
!     ! mz_bk_20120718-

!     ! mz_bk_20120515+
!     ! updated from HAMOCC mpiesm-1.0.00
!     zinges = 0.6_dp              ! dimensionless fraction -
!     ! !                            ! assimilation efficiency of zooplankton
!     ! ! "old" value
!     ! zinges=0.5_dp              ! dimensionless fraction - 
!     !                          ! assimilation efficiency of zooplankton
!     ! mz_bk_20120515-

!     epsher = 0.8_dp              ! dimensionless fraction - 
!     !                            ! (1-epsher)=fraction of grazing egested

!     ! mz_bk_20120515+
!     ! updated from HAMOCC mpiesm-1.0.00
!     grazra = 1.0_dp * dtb        !1/d -grazing rate [emr: 0.6-0.9]
!     ! ! "old" value
!     ! grazra=0.8_dp * dtb        !1/d -grazing rate [emr: 0.6-0.9]
!     ! mz_bk_20120515-

!     ! mz_bk_20120515+
!     ! updated from HAMOCC mpiesm-1.0.00
!     ! mz_bk_20120709+
!     ! testing spemor = 5.e6 * dtb
!     spemor = 3.e6_dp * dtb       !1/d -mortality rate of zooplankton
!     ! spemor = 5.e6_dp * dtb       !1/d -mortality rate of zooplankton
!     ! mz_bk_20120709-
!     ! "old" value
!     !! DON'T USE this old value anymore,
!     !! because calculation of zoomor was updated:
!     !! old: zoomor=spemor*zoothresh
!     !! new: zoomor=spemor*zoothresh*zoothresh
!     !! with zoothresh=MAX(0., (ocetra(i,j,k,izoo) - 2. * grami))
!     !  spemor=0.008*dtb  !1/d -mortality rate of zooplankton
!     ! mz_bk_20120515-

!     gammap = 0.03_dp * dtb       !1/d -exudation rate
!     ! mz_bk_20120718+
!     ! testing gammaz = 0.03_dp * dtb
!     gammaz = 0.06_dp * dtb       !1/d -excretion rate
!     ! gammaz = 0.03_dp * dtb       !1/d -excretion rate
!     ! mz_bk_20120718-
!     ecan = 0.95_dp               ! fraction of mortality as PO_4
!     ! m_bk_20120710+
!     ! test for DMS scaling pi_alpha = 0.005_dp
!     ! pi_alpha = 0.005_dp         ! initial slope of production vs irradiance
!     ! mz_bk_20120710-
!     pi_alpha = 0.02_dp         ! initial slope of production vs irradiance
!                                ! curve (alpha) (0.002 for 10 steps per day)
!                                ! note: in mpiesm-1.0.00
!                                ! dtb is not included in pi_alpha !

!     ! mz_bk_20120515+
!     ! added from HAMOCC mpiesm-1.0.00
!     fPAR     = 0.4_dp            ! fraction of Photosynthetic Active Radiation
!     ! mz_bk_20120515-

!     ! half sat. constants, note that the units are kmol/m3
!     ! (conc0 in hamocc3.1)
!     ! mz_bk_20120718+
!     ! testing bkphy  = 4.e-8_dp
!     bkphy  = 1.e-8_dp  !i.e. 0.04 mmol P/m3 |js: 0.01 vs. 0.04? check 0.4 16.9.
!     ! bkphy  = 4.e-8_dp!i.e. 0.04 mmol P/m3 |js: 0.01 vs. 0.04? check 0.4 16.9.
!     ! mz_bk_20120718-
!     bkzoo  = 4.e-8_dp  !i.e. 0.04 mmol P/m3
!     bkopal = 1.e-6_dp  !i.e. 1.0  mmol Si/m3
!     !js: no bkiron? (fitzwater et al 1996) et al. 200 120 nMol/m3,
!     ! moore et al 80nMol/m3)

!     !sinking speed
!     ! mz_bk_20120626+
!     sinkspeed_poc  = 5.0_dp      ! daily sinking speed of poc
!     sinkspeed_opal = 30.0_dp     ! daily sinking speed of opal
!     ! "old" value:
!     ! sinkspeed_opal = 60.0_dp     ! daily sinking speed of opal
!     sinkspeed_cal  = 30.0_dp     ! daily sinking speed of cal
!     ! "old" values:
!     ! wpoc  = 5.*dtb  !m/d  iris : 5. 
!     !                 !  (for dtb=.1 (10 time steps / day) -> wpoc=1) pw:10
!     ! wcal  = 30.*dtb !m/d 
!     ! ! mz_bk_20120515+
!     ! ! updated from HAMOCC mpiesm-1.0.00
!     ! !  wopal = 30.0_dp*dtb ! daily sinking speed of opal 
!     ! ! "old" value
!     ! wopal = 60.*dtb !m/d  iris : 60 pw: 30
!     ! ! mz_bk_20120515-
!     ! mz_bk_20120626-

      
!     ! water column remineralisation constants

!     drempoc  = 0.025_dp * dtb    !1/d ! 0.75/month.
!     !                                 ! H3.1: 0.2/month k=1, 0.1/month k>1
!     ! mz_bk_20120515+
!     ! updated from HAMOCC mpiesm-1.0.00
!     dremdoc  = 0.004_dp * dtb    !1/d
!     ! "old" value
!     ! dremdoc  = 0.03_dp * dtb     !1/d
!     ! mz_bk_20120515-

!     dphymor  = 0.1_dp * dtb      !1/d
!     dzoomor  = 0.02_dp * dtb     !1/d

!     ! mz_bk_20120515+
!     ! updated from HAMOCC mpiesm-1.0.00
!     dremopal  = 0.01_dp * dtb    !1/d
!     ! "old" value
!     ! dremopal = 0.0075_dp * dtb   !1/d ! 0.01 -> 0.001 js10072006 : 
!     !                                 ! slightly overdone -->0.0075
!     ! mz_bk_20120515-

!     ! mz_bk_20120515+
!     ! updated from HAMOCC mpiesm-1.0.00
!     dremn2o  = 0.01_dp * dtb     !1/d
!     ! "old" value
!     ! dremn2o  = 0.1_dp * dtb      !1/d
!     ! mz_bk_20120515-

!     ! mz_bk_20120629+
!     ! mz_bk_20120709+
!     sulfate_reduction  = 0.005_dp    
!     ! sulfate_reduction  = 0.0_dp    
!     ! mz_bk_20120709-
!     ! mz_bk_20120629-
!     dremcalc = 0.075_dp * dtb    !1/d ! 0.2 -> 0.02 js10072006 : 
!     !                                 ! slightly overdone --> 0.075

! !!$!!#ifdef AGG
!     IF (L_AGG) THEN 
!        drempoc  = 0.1_dp * dtb      !1/d       re-introduced 09062006
!        dremopal = 3.3e-3_dp * dtb   ! js 4.7.2006 0.0033 from .01/3. (60/20 m/d)
!        dphymor  = 0.2_dp * dtb      !1/d
!     ENDIF
! !!$!!#endif
      
!     ! nitrogen fixation by blue green algae (cyano.f90)

!     ! mz_bk_20120515+
!     ! updated from HAMOCC mpiesm-1.0.00
!     ! bluefix = 0.005_dp * dtb     !1/d
!     ! mz_bk_20120709+
!     n2_fixation = 0.005_dp
!     ! n2_fixation = 0.02_dp
!     ! mz_bk_20120709-
!     ! "old" value
!     ! bluefix=0.02*dtb     !1/d 
!     ! n2_fixation = 0.02_dp
!     ! mz_bk_20120515-

!     ! mz_bk_20120626+
!     ! total denitrification rate is a fraction of aerob remineralisation rate
!     ! drempoc

!     !  dremn3o = 0.05_dp*drempoc           ! 1/d
!     ! mz_bk_20120709+
!     denitrification = 0.05_dp           ! 1/d
!     ! denitrification = 0.5_dp           ! 1/d
!     ! mz_bk_20120709-
!     ! mz_bk_20120626-

!     ! extended redfield ratio declaration
!     ! Note: stoichiometric ratios are based on Takahashi etal. (1985)
!     ! P:N:C:-O2 + 1:16:122:172

!     ro2ut = 172._dp
!     rcar  = 122._dp
!     rnit  = 16._dp
!     rnoi  = 1._dp / rnit

!     ! mz_bk_20120626+
!     ! not longer needed, in new version:
!     ! rnit23 -> nitdem
!     ! rnit13 -> n2prod
!     ! rnit23 = ro2ut*2./3. !114 rno3*2/3 for nitrate reduction
!     ! !                    ! during denitrification
!     ! rnit13 = ro2ut*1./3. !57  rno3*1/3 for n2 production
!     ! !                    !during denitrification
!     ! N consumption of denitrification corrected after Paulmier etal, 2009)
!     nitdem = 121.6_dp      ! nitrate demand to remin. 1 mol P in suboxic water
!     n2prod = 68.8_dp       ! N2 production for 1 mol P remineralized
!     !                      ! in suboxic water
!     ! mz_bk_20120626-

!     rcalc  = 35._dp   ! iris 40 !calcium carbonate to organic phosphorous
!     !                 ! production ratio
!     ropal  = 25._dp   ! iris 25 !opal to organic phosphorous production ratio
!     calmax = 0.20_dp  ! max fraction of CaCO3 production (AGG) (max rain ratio)
!     gutc   = 0._dp    ! fraction of caco3 that dissolves during passage
!     !                 ! through zooplankton guts (not used yet)

! #if defined MPIOM_13B
!     ! parameters for sw-radiation attenuation
!     ! Analogue to Moore et al., Deep-Sea Research II 49 (2002), 403-462
!     ! 1 kmolP = (122*12/60)*10^6 mg[Chlorophyll] 
!     ! js: in Moore et al. atten_w and atten_c are mean values 
!     ! for the mixed layer. (the values just appear in their appendix,
!     ! without further reference. can they be used here?)
!     ! here we might need smaller values, e.g. scale by layer depth / mld

!     ctochl  = 60._dp        ! C to Chlorophyll ratio     ! HadOCC: 40.
!     atten_w = 0.04_dp       ! Gelbstoff attenuation in 1/m
!     atten_c = 0.03_dp * rcar * (12._dp / ctochl) * 1.e6_dp  ! phytoplankton
!     !                       ! attenuation in 1/m "self-shading" 7.32e5 mg Chl/m3
!     atten_f = 0.4_dp        ! fraction of sw-radiation directly absorbed in
!     !                       ! surface layer 
!     !                       ! (only used if FB_BGC_OCE) [feedback bgc-ocean]
! #endif

!     !ik for interaction with sediment module
!     ! mz_bk_20120724+
!     ! never used!!
!     ! o2ut = 172._dp
!     ! mz_bk_20120724-
!     rno3 =  16._dp

!     !ik weight percent iron in dust deposition (0.035) times Fe solubility
!     ! (0.01) /55.85 g--> mol
!     perc_diron = 0.035_dp * 0.01_dp / 55.85_dp
!     ! the three values below are from Johnson et al., 1997
!     ! Mar. Chemistry 57, 137-161
!     !  riron = 5._dp * rcar * 1.e-6_dp  ! 'Q'=5.*122.*1.e-6 = 6.1e-4
!     !  !                  ! (iron to phosphate stoichiometric ratio * umol->mol)
!     !  riron = 2._dp * rcar * 1.e-6_dp  ! 14.2.06: 5.->2. (2.5umol P/nmol Fe)
!     !  !                                ! emr: 3.*rcar*1.e-6
!     riron = 3._dp * rcar * 1.e-6_dp  ! 06.3.06: 2.->3. coex90 = 12.2GTC/a
!     !                                !          before change
!     fesoly  = 0.6_dp * 1.e-9_dp      ! global mean/max (?) dissolved iron
!     !                  ! concentration in deep water (for relaxation) [mol/m3]
!     !                  ! 'fesoly' stands for 'apparent solubility value' of iron

!     ! mz_bk_20120626+
!     relaxfe = 0.05_dp / 365._dp * dtb   ! relaxation time for iron to fesoly
!                                         ! corresponds to 20 yrs
!     !  relaxfe = 0.005_dp / 365._dp * dtb  ! changed as suggested by emr
!     !  !                                   ! (should probably be named
!     !  !                                   ! scavenge_fe or remove_fe)
!     ! (relaxfe is only used to remove iron, corresponding to Johnson's 200yr
!     ! residence time)
!     ! but 200 yrs might be too long, as iron concentrations in the Atlantic
!     ! are too high.
!     !  relaxfe = 0.5/365.*dtb              ! relaxation time for iron to fesoly
!     !  !                 ! ->1.37e-3*dtb(iris' paper states 0.005!?)back15206js 
!     ! mz_bk_20120626-

!     !
!     ! Set constants for calculation of DMS ( mo_carbch )
!     ! Parameters are a result from Kettle optimisation
!     ! (fit to Kettle data, version emr, not tinka)

!     dmspar(6)=0.100000000E-07_dp        !0 half saturation microbial

!     ! mz_bk_20120515+
!     ! updated from HAMOCC mpiesm-1.0.00
!     ! seems both paramter values were swaped
!     ! mz_bk_20120719+
!     ! dmspar(5)=1.25*0.109784522E-01_dp  ! production with delsil (diatoms)
!     ! dmspar(4)=1.25*0.107638502E+00_dp  ! production with delcar
!     dmspar(5)=0.5*0.109784522E-01_dp  ! production with delsil (diatoms)
!     dmspar(4)=0.5*0.107638502E+00_dp  ! production with delcar
!     ! mz_bk_20120719-
!     !                                  ! (coccolithoporides)
!     ! "old" values
!     ! dmspar(5)=1.25*0.107638502E+00  !2*1.3e-5 production with delsil (diatoms)
!     ! dmspar(4)=1.25*0.109784522E-01  !2*0.02   production with delcar
!     ! !                               !         (coccolithoporides)
!     ! mz_bk_20120515-

!     ! mz_bk_20120612+
!     ! These are in units/day, should we multiply with dtb to get
!     ! units/timestep??
!     dmspar(3)=0.0096  !4.8e-5       !2*1.6e-3 microbial consumption
!     dmspar(2)=0.0075  !0.0003       !2*0.005  photolysis (uv-destruction)
!     ! mz_bk_20120711 test values Kloster et al. 2006
!     ! dmspar(3)=0.0011_dp  !4.8e-5       !2*1.6e-3 microbial consumption
!     ! dmspar(2)=0.1728_dp  !0.0003       !2*0.005  photolysis (uv-destruction)
!     ! dmspar(3)=0.0096*dtb  !4.8e-5       !2*1.6e-3 microbial consumption
!     ! dmspar(2)=0.0075*dtb  !0.0003       !2*0.005  photolysis (uv-destruction)
!     ! mz_bk_20120612-

!     dmspar(1)=10._dp     !2*5.   temperature dependent release by phytoplankton
    ! mz_bk_20120724-

    ! mz_bk_20120725+
    rnoi     = 1._dp / rnit   ! inverse of rnit
    atten_c = 0.03_dp * rcar * (12._dp / ctochl) * 1.e6_dp  ! attenuation by
    !                                                       ! phytoplankton
    ! mz_bk_20120725-

    ! mz_bk_20120724+
    WRITE(*,*)  '***************************************************'
    WRITE(*,*)  '* Values of BELEG_BGC variables : '
    WRITE(*,*)  '***************************************************'
    WRITE(*,*)  '*    ro2ut             = ',ro2ut   
    WRITE(*,*)  '*    rcar              = ',rcar 
    WRITE(*,*)  '*    rnit              = ',rnit
    WRITE(*,*)  '*    rnoi              = ',rnoi
    WRITE(*,*)  '*    nitdem            = ',nitdem
    WRITE(*,*)  '*    n2prod            = ',n2prod
    WRITE(*,*)  '*    rcalc             = ',rcalc
    WRITE(*,*)  '*    ropal             = ',ropal
    WRITE(*,*)  '*    n2_fixation       = ',n2_fixation
    WRITE(*,*)  '*    rno3              = ',rno3
    WRITE(*,*)  '*    perc_diron        = ',perc_diron
    WRITE(*,*)  '*    riron             = ',riron
    WRITE(*,*)  '*    fesoly            = ',fesoly
    WRITE(*,*)  '*    relaxfe           = ',relaxfe
    WRITE(*,*)  '*    pi_alpha          = ',pi_alpha
    WRITE(*,*)  '*    fPAR              = ',fPAR
#if defined MPIOM_13B
    WRITE(*,*)  '*    ctochl            = ',ctochl
    WRITE(*,*)  '*    atten_w           = ',atten_w
    WRITE(*,*)  '*    atten_f           = ',atten_f
    WRITE(*,*)  '*    atten_c           = ',atten_c
#endif
    WRITE(*,*)  '*    phytomi           = ',phytomi
    WRITE(*,*)  '*    bkphy             = ',bkphy
    WRITE(*,*)  '*    bkopal            = ',bkopal    
    WRITE(*,*)  '*    remido            = ',remido
    WRITE(*,*)  '*    dyphy             = ',dyphy
    WRITE(*,*)  '*    gammap            = ',gammap
    WRITE(*,*)  '*    bkzoo             = ',bkzoo    
    WRITE(*,*)  '*    grami             = ',grami
    WRITE(*,*)  '*    zinges            = ',zinges
    WRITE(*,*)  '*    epsher            = ',epsher
    WRITE(*,*)  '*    grazra            = ',grazra
    WRITE(*,*)  '*    spemor            = ',spemor
    WRITE(*,*)  '*    gammaz            = ',gammaz
    WRITE(*,*)  '*    ecan              = ',ecan    
    WRITE(*,*)  '*    sinkspeed_poc     = ',sinkspeed_poc
    WRITE(*,*)  '*    sinkspeed_opal    = ',sinkspeed_opal   
    WRITE(*,*)  '*    sinkspeed_cal     = ',sinkspeed_cal    
    WRITE(*,*)  '*    drempoc           = ',drempoc    
    WRITE(*,*)  '*    dremdoc           = ',dremdoc   
    WRITE(*,*)  '*    dremn2o           = ',dremn2o   
    WRITE(*,*)  '*    denitrification   = ',denitrification   
    WRITE(*,*)  '*    sulfate_reduction = ',sulfate_reduction   
    WRITE(*,*)  '*    dremopal          = ',dremopal   
    WRITE(*,*)  '*    dphymor           = ',dphymor   
    WRITE(*,*)  '*    dzoomor           = ',dzoomor   
    WRITE(*,*)  '*    dmspar(1)         = ',dmspar(1)
    WRITE(*,*)  '*    dmspar(2)         = ',dmspar(2)
    WRITE(*,*)  '*    dmspar(3)         = ',dmspar(3)
    WRITE(*,*)  '*    dmspar(4)         = ',dmspar(4)
    WRITE(*,*)  '*    dmspar(5)         = ',dmspar(5)    
    WRITE(*,*)  '*    dmspar(6)         = ',dmspar(5)    
    WRITE(*,*)  '***************************************************'
    
    ! mz_bk_20120724-

    !WRITE(io_stdo_bgc,*)  '***************************************************'
    !WRITE(io_stdo_bgc,*)  '* Values of BELEG_BGC variables : '
    !WRITE(io_stdo_bgc,*)  '***************************************************'
    !WRITE(io_stdo_bgc,*)  '*    phytomi      = ',phytomi
    !WRITE(io_stdo_bgc,*)  '*    grami        = ',grami
    !WRITE(io_stdo_bgc,*)  '*    remido       = ',remido
    !WRITE(io_stdo_bgc,*)  '*    dyphy        = ',dyphy
    !WRITE(io_stdo_bgc,*)  '*    zinges       = ',zinges
    !WRITE(io_stdo_bgc,*)  '*    epsher       = ',epsher
    !WRITE(io_stdo_bgc,*)  '*    grazra       = ',grazra
    !WRITE(io_stdo_bgc,*)  '*    spemor       = ',spemor
    !WRITE(io_stdo_bgc,*)  '*    gammap       = ',gammap
    !WRITE(io_stdo_bgc,*)  '*    gammaz       = ',gammaz
    !WRITE(io_stdo_bgc,*)  '*    ecan         = ',ecan    
    !WRITE(io_stdo_bgc,*)  '*    bkphy        = ',bkphy
    !WRITE(io_stdo_bgc,*)  '*    bkzoo        = ',bkzoo    
    !WRITE(io_stdo_bgc,*)  '*    bkopal       = ',bkopal    
    !WRITE(io_stdo_bgc,*)  '*    wpoc         = ',wpoc
    !WRITE(io_stdo_bgc,*)  '*    wcal         = ',wcal    
    !WRITE(io_stdo_bgc,*)  '*    wopal        = ',wopal   
    !WRITE(io_stdo_bgc,*)  '*    drempoc      = ',drempoc    
    !WRITE(io_stdo_bgc,*)  '*    dremdoc      = ',dremdoc   
    !WRITE(io_stdo_bgc,*)  '*    dremopal     = ',dremopal   
    !WRITE(io_stdo_bgc,*)  '*    dphymor      = ',dphymor   
    !WRITE(io_stdo_bgc,*)  '*    dzoomor      = ',dzoomor   
    !WRITE(io_stdo_bgc,*)  '*    bluefix      = ',bluefix   
    !WRITE(io_stdo_bgc,*)  '*    ro2ut        = ',ro2ut   
    !WRITE(io_stdo_bgc,*)  '*    rcar         = ',rcar 
    !WRITE(io_stdo_bgc,*)  '*    rnit         = ',rnit
    !WRITE(io_stdo_bgc,*)  '*    rnoi         = ',rnoi
    !WRITE(io_stdo_bgc,*)  '*    rnit23       = ',rnit23
    !WRITE(io_stdo_bgc,*)  '*    rnit13       = ',rnit13
    !WRITE(io_stdo_bgc,*)  '*    rcalc        = ',rcalc
    !WRITE(io_stdo_bgc,*)  '*    ropal        = ',ropal
    !WRITE(io_stdo_bgc,*)  '*    gutc         = ',gutc
    !WRITE(io_stdo_bgc,*)  '*    ctochl       = ',ctochl
    !WRITE(io_stdo_bgc,*)  '*    atten_w      = ',atten_w
    !WRITE(io_stdo_bgc,*)  '*    atten_c      = ',atten_c
    !WRITE(io_stdo_bgc,*)  '*    atten_f      = ',atten_f
    !WRITE(io_stdo_bgc,*)  '*    o2ut         = ',o2ut
    !WRITE(io_stdo_bgc,*)  '*    rno3         = ',rno3
    !WRITE(io_stdo_bgc,*)  '*    perc_diron   = ',perc_diron
    !WRITE(io_stdo_bgc,*)  '*    riron        = ',riron
    !WRITE(io_stdo_bgc,*)  '*    fesoly       = ',fesoly
    !WRITE(io_stdo_bgc,*)  '*    relaxfe      = ',relaxfe
    !WRITE(io_stdo_bgc,*)  '*    dmspar(1)    = ',dmspar(1)
    !WRITE(io_stdo_bgc,*)  '*    dmspar(2)    = ',dmspar(2)
    !WRITE(io_stdo_bgc,*)  '*    dmspar(3)    = ',dmspar(3)
    !WRITE(io_stdo_bgc,*)  '*    dmspar(4)    = ',dmspar(4)
    !WRITE(io_stdo_bgc,*)  '*    dmspar(5)    = ',dmspar(5)    
    !WRITE(io_stdo_bgc,*)  '***************************************************'

!!$!!#ifndef AGG
    IF (.not.L_AGG) THEN
       dustd1 = 0.0001_dp            !cm = 1 um, boundary between clay and silt
       dustd2 = dustd1 * dustd1
       dustsink = (9.81_dp * 86400._dp / 18._dp                               &
            ! g * sec per day / 18. | js: Stoke's law for small particles
            &     * (claydens - 1025._dp) / 1.567_dp * 1000._dp               &
            ! excess density / dyn. visc. | -> cm/s to m/day
            &     * dustd2 * 1.e-4_dp) * dtb
            ! *diameter**2 |*1000 *1.e-4?
       wdust = dustsink

       WRITE(io_stdo_bgc,*)                                                   &
            &'*                              dustd1       = ',dustd1
       WRITE(io_stdo_bgc,*)                                                   &
            &'*                              dustd2       = ',dustd2 
       WRITE(io_stdo_bgc,*)                                                   &
            &'*                              dustsink     = ',dustsink
       WRITE(io_stdo_bgc,*)                                                   &
            &'*                              wdust        = ',wdust
    ENDIF
!!$!!#endif

    WRITE(io_stdo_bgc,*)                                                      &
         &'****************************************************************'

    !--------------------------------------------------------------------------
    ! Initialization of aggregation
    ! from subroutine ini_aggregation in mo_beleg_bgc.f90
    !--------------------------------------------------------------------------
!!$!! #ifdef AGG
    IF (L_AGG) THEN
       ! parameters needed for the aggregation module
       ! (see Kriest 2002, DSR I vol.49, p. 2133-2162)

       SinkExp = 0.62_dp           ! exponent of the sinking speed vs. diameter
       !                           ! relationship
       FractDim = 1.62_dp          ! exponent of the diameter vs. phosphorous
       !                           ! content relationship
       Stick = 0.40_dp             ! maximum stickiness
       cellmass = 0.012_dp / rnit  ! [nmol P]; minimum mass of a particle in
       !                           ! phosphorous units (rnit=16)
       !ik  cellmass = 0.0039_dp / rnit ! [nmol P] for 10 um diameter
       cellsink = 1.40_dp *dtb     ! [m/d]; see Kriest 2002, Table 2 Ref 8 (from
       !                           ! Stokes' formula, delta rho 0.052 g/cm3)
       !ik  cellsink = 0.911_dp * dtb   ! [m/d]  for a 10 um diameter
       shear = 86400._dp           ! wind induced shear in upper 100m , 1 d^-1
       fsh = 0.163_dp * shear *dtb ! turbulent shear (used for aggregation)
       fse = 0.125_dp * 3.1415927_dp * cellsink * 100._dp ! differential
       !                    ! settling (used for aggregation) (100=10**2 [d**2])
       alow1 = 0.002_dp            ! diameter of smallest particle [cm]
       !ik  alow1 = 0.001_dp            ! diameter of smallest particle [cm]
       alow2 = alow1 * alow1
       alow3 = alow2 * alow1
       alar1 = 1.0_dp              ! diameter of largest particle for size
       !                           ! dependend aggregation and sinking [cm]
       vsmall = 1.e-9_dp
       safe = 1.e-6_dp
       pupper = safe / ((FractDim + safe) * cellmass) ! upper boundary for
       !                                              ! cells per aggregate (?)
       plower = 1._dp / (1.1_dp * cellmass)           ! lower boundary for
       !                                              ! cells per aggregate (?)
       zdis = 0.01_dp / ((FractDim + 0.01_dp) * cellmass)

       !ik check max possible sinking speed in relation to min.
       !ik layer thinkness and time step for all standard layers, except
       !ik the bottom layer.
       !ik if max possible sinking speed (per time step) is greater
       !ik than min layer thickness, decrease max. length for sinking and
       !ik aggregation
      
       !DONE in messy_hamocc_e5
       ! zmini = 8000._dp
       !    DO j=1,je
       !       DO i=1,ie
       !          DO k=1,kbo(i,j)-1
       !             if(ddpo(i,j,k) .gt. 0.5_dp) then
       !                zmini = min(ddpo(i,j,k),zmini)
       !             endif 
       !          ENDDO
       !       ENDDO
       !    ENDDO
       !
       ! CALL global_min(zmini)
             
       checksink =(zmini / cellsink)**(1._dp / SinkExp) * alow1 
 
       if (checksink .lt. alar1) then
 
          write(io_stdo_bgc,*) 'Allowed max. length for sinking'              &
               & , 'with min. depth of '                                      &
               & , zmini, ' m for layers 1-(kbo-1) and time step of ', dtb    &
               & , ' days is' , checksink                                     &
               & ,'cm, which is smaller than prescribed value of', alar1, ' cm'

          talar1 = alar1
          alar1 = checksink
          write(io_stdo_bgc,*) 'Set max. length for sinking and aggregation'  &
               &, 'from ',talar1,' to ', alar1

       endif

       alar2 = alar1 * alar1
       alar3 = alar2 * alar1
       TSFac = (alar1 / alow1)**SinkExp
       TMFac = (alar1 / alow1)**FractDim
 
       !ik check the maximum possible sinking speed for the bottom layer (which
       !ik may be smaller than zmini, and write to array alar1max, tsfmax, 
       !ik tmfmax

       DO j=1,je
          DO i=1,ie
             alar1max(i,j) = alar1
             TSFmax(i,j) = TSFac
             TMFmax(i,j) = TMFac
             if (ddpo(i,j,kbo(i,j)) .gt. 0.5_dp) then

                !ik evaluate safe length scale for size dependent sinking and
                !ik aggregation, and the resulting sinking rate and
                !ik aggregation rate.

                checksink = (ddpo(i,j,kbo(i,j)) / cellsink)**(1._dp / SinkExp)&
                     &      * alow1
                if (checksink .lt. alar1) then
                   alar1max(i,j) = checksink
                   TSFmax(i,j) = (checksink / alow1)**SinkExp
                   TMFmax(i,j) = (checksink / alow1)**FractDim
                   write(io_stdo_bgc,*) 'resetting alar1 to',checksink        &
                        & , 'at i =', i,' j = ', j, ' k = ', kbo(i,j)         &
                        & , ' with dz = ', ddpo(i,j,kbo(i,j))  
                endif
             ENDIF
          ENDDO
       ENDDO
 
       ! for shear aggregation of dust:
       dustd1 = 0.0001_dp          ![cm] = 1 um, boundary between clay and silt
       dustd2 = dustd1 * dustd1
       dustd3 = dustd2 * dustd1
       dustsink = (9.81_dp * 86400._dp / 18._dp                               &
            ! g * sec per day / 18.                 
            &     * (claydens - 1025._dp) / 1.567_dp * 1000._dp               &
            !excess density / dyn. visc.
            &     * dustd2 * 1.e-4_dp) * dtb                  ! --> 4.73e-2 m/d
       write(io_stdo_bgc,*) 'dust diameter (cm)', dustd1
       write(io_stdo_bgc,*) 'dust sinking speed (m/d)', dustsink / dtb
       if (dustsink .gt. cellsink) then 
          write(io_stdo_bgc,*) 'dust sinking speed greater than cellsink'
          dustsink = cellsink
          write(io_stdo_bgc,*) 'set dust sinking speed to cellsink'
       endif
    ENDIF
!!$!! #endif /*AGG*/  

  END SUBROUTINE hamocc_beleg
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE hamocc_init_var(ie,je,ke,ddpo)
    !
    ! corresponds to subroutines ini_aquatic_tracers and
    ! ini_pore_water_tracersin mo_beleg_bgc.f90
    ! Initialization of aquatic (advected) ocean tracers (new start)
    ! Initialization of sediment pore water tracers (new start)
    ! Initialization of fluxes to sediment
    ! Initialization of fields for diagnostic output
    !-------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: ie,je,ke
    REAL(DP),INTENT(IN) :: ddpo(ie,je,ke)

    !LOCAL
    INTEGER :: i,j,k


    DO j=1,je
       DO i=1,ie

          DO k=1,8
             chemcm(i,j,k)=rmasks
          ENDDO

          DO  k=1,ke
             aksp(i,j,k)=rmasks
          ENDDO

          !
          ! Initial values for aquatic (advected) ocean tracers (for new start)
          ! 
          DO k=1,ke
             IF (ddpo(i,j,k) .GT. 0.5_dp) THEN                 ! wet points
                ocetra(isco212)%ptr(i,j,k) = 2.27e-3_dp        ! [kmol/m3]
                ocetra(ialkali)%ptr(i,j,k) = 2.37e-3_dp
                ocetra(iphosph)%ptr(i,j,k) = 2.17e-6_dp
                ocetra(ioxygen)%ptr(i,j,k) = 3.e-4_dp
                ocetra(igasnit)%ptr(i,j,k) = 0._dp
                ocetra(iano3)%ptr(i,j,k)   = 2.17e-6_dp * rnit ! old 32.e-6
                ocetra(isilica)%ptr(i,j,k) = 1.2e-4_dp
                ocetra(idoc)%ptr(i,j,k)    = 1.e-10_dp
                ocetra(iphy)%ptr(i,j,k)    = 1.e-8_dp
                ocetra(izoo)%ptr(i,j,k)    = 1.e-8_dp
                ocetra(idet)%ptr(i,j,k)    = 1.e-8_dp
                ocetra(icalc)%ptr(i,j,k)   = 0._dp
                ocetra(iopal)%ptr(i,j,k)   = 1.e-8_dp
                ocetra(ian2o)%ptr(i,j,k)   = 0._dp
                ocetra(idms)%ptr(i,j,k)    = 0._dp
                ocetra(ifdust)%ptr(i,j,k)  = 0._dp
                ocetra(iiron)%ptr(i,j,k)   = 0.6_dp * 1.e-9_dp
                hi(i,j,k)                  = 3.e-9_dp
                co3(i,j,k)                 = 0._dp   ! this good for
                !                                    ! initialisation -> 2.e-4?
!!$!! #ifdef AGG
                IF (L_AGG) THEN
                   ! calculate initial numbers from mass, 
                   ! to start with appropriate size distribution
                   snow = (ocetra(iphy)%ptr(i,j,k)                            &
                        & + ocetra(idet)%ptr(i,j,k)) * 1.e+6_dp
                   ocetra(inos)%ptr(i,j,k) = snow / cellmass / (FractDim+1._dp)
                   ocetra(iadust)%ptr(i,j,k) = 0._dp
                ENDIF
!!$!! #endif /*AGG*/
             ELSE                              ! dry points
                ocetra(iphosph)%ptr(i,j,k)   = rmasks
                ocetra(isilica)%ptr(i,j,k)   = rmasks
                ocetra(ioxygen)%ptr(i,j,k)   = rmasks
                ocetra(ialkali)%ptr(i,j,k)   = rmasks
                ocetra(isco212)%ptr(i,j,k)   = rmasks
                ocetra(iano3)%ptr(i,j,k)     = rmasks
                ocetra(igasnit)%ptr(i,j,k)   = rmasks
                ocetra(idoc)%ptr(i,j,k)      = rmasks
                ocetra(iphy)%ptr(i,j,k)      = rmasks
                ocetra(izoo)%ptr(i,j,k)      = rmasks
                ocetra(idet)%ptr(i,j,k)      = rmasks
                ocetra(icalc)%ptr(i,j,k)     = rmasks
                ocetra(iopal)%ptr(i,j,k)     = rmasks
                ocetra(ian2o)%ptr(i,j,k)     = rmasks
                ocetra(idms)%ptr(i,j,k)      = rmasks
                ocetra(ifdust)%ptr(i,j,k)    = rmasks
                ocetra(iiron)%ptr(i,j,k)     = rmasks
                hi(i,j,k)                    = rmasks
                co3(i,j,k)                   = rmasks
                co2(i,j)                     = rmasks 
!!$!! #ifdef AGG
                IF(L_AGG) THEN
                   ocetra(inos)%ptr(i,j,k)   = rmasks
                   ocetra(iadust)%ptr(i,j,k) = rmasks 
                ENDIF
!!$!! #endif /*AGG*/
             ENDIF
          ENDDO

          !
          ! Initial values for sediment pore water tracers. (solid components?)
          !
          DO  k=1,ks
             IF (bolay(i,j) .GT. 0._dp) THEN
                powtra(ipowaic)%ptr(i,j,k) = ocetra(isco212)%ptr(i,j,kbo(i,j))
                powtra(ipowaal)%ptr(i,j,k) = ocetra(ialkali)%ptr(i,j,kbo(i,j))
                powtra(ipowaph)%ptr(i,j,k) = ocetra(iphosph)%ptr(i,j,kbo(i,j))
                powtra(ipowaox)%ptr(i,j,k) = ocetra(ioxygen)%ptr(i,j,kbo(i,j))
                powtra(ipown2)%ptr(i,j,k)  = 0._dp
                powtra(ipowno3)%ptr(i,j,k) = ocetra(iano3)%ptr(i,j,kbo(i,j))
                powtra(ipowasi)%ptr(i,j,k) = ocetra(isilica)%ptr(i,j,kbo(i,j))
                sedlay(issso12)%ptr(i,j,k) = 1.e-8_dp
                sedlay(isssc12)%ptr(i,j,k) = 1.e-8_dp
                sedlay(issster)%ptr(i,j,k) = 30._dp
                sedlay(issssil)%ptr(i,j,k) = 3._dp
                sedhpl(i,j,k)              = hi(i,j,kbo(i,j))
             ELSE
                powtra(ipowno3)%ptr(i,j,k) = rmasks   ! pore water
                powtra(ipown2)%ptr(i,j,k)  = rmasks
                powtra(ipowaic)%ptr(i,j,k) = rmasks
                powtra(ipowaal)%ptr(i,j,k) = rmasks
                powtra(ipowaph)%ptr(i,j,k) = rmasks
                powtra(ipowaox)%ptr(i,j,k) = rmasks
                powtra(ipowasi)%ptr(i,j,k) = rmasks
                sedlay(issso12)%ptr(i,j,k) = rmasks   ! solid sediment 
                sedlay(isssc12)%ptr(i,j,k) = rmasks
                sedlay(issssil)%ptr(i,j,k) = rmasks
                sedlay(issster)%ptr(i,j,k) = rmasks
                sedhpl(i,j,k)              = rmasks
             ENDIF
          ENDDO
          !
          ! values for sediment fluxes
          !
!!$!!        DO k=1,npowtra
!!$!!          sedfluxo(i,j,k)=0.     
!!$!!        ENDDO

          !
          ! fluxes of organic carbon, caco3, opal, and dust to the sediment
          ! (js: still old [bad] nomenclature from HAMOCC3....)
          !
          ! change default from 0. to rmasko
          prorca(i,j)=0._dp
          prcaca(i,j)=0._dp
          silpro(i,j)=0._dp
          produs(i,j)=0._dp

          ! mz_bk_20120607+
          IF (L_BGC_DIAG) THEN
             DO k=1,ke
                bgcprod(kphosy)%ptr(i,j,k)    = 0._dp
                bgcprod(kgraz)%ptr(i,j,k)     = 0._dp
                bgcprod(kexport)%ptr(i,j,k)   = 0._dp
                bgcprod(kdelcar)%ptr(i,j,k)   = 0._dp
                bgcprod(kdelsil)%ptr(i,j,k)   = 0._dp
                bgcprod(kdmsprod)%ptr(i,j,k)  = 0._dp
                bgcprod(kdms_bac)%ptr(i,j,k)  = 0._dp
                bgcprod(kdms_uv)%ptr(i,j,k)   = 0._dp
                ! mz_bk_20120813+
                ! bgcprod(krain)%ptr(i,j,k)     = 0._dp
                ! mz_bk_20120813-
                bgcprod(kflim)%ptr(i,j,k)     = 0._dp
                bgcprod(kplim)%ptr(i,j,k)     = 0._dp
                bgcprod(knlim)%ptr(i,j,k)     = 0._dp
             ENDDO
          ENDIF
          ! mz_bk_20120607-

       ENDDO
    ENDDO

    DO k=1,ke
       DO j=1,je
          DO i=1,ie
             abs_oce(i,j,k) = 1._dp
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE hamocc_init_var
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE hamocc_chemcon(ie,je,ke,sao,tho,ddpo,tiestu)
    ! Written by A.Pozzer , MPICH, 2007/11/08
    ! Based on CHEMCON written by Ernst Maier-Reimer,*MPI-Met, HH* 10.04.01
    !**   Interface to ocean model (parameter list):
    !     -----------------------------------------
    !
    !     *INTEGER* *IE*    - 1st dimension of model grid.
    !     *INTEGER* *JE*    - 2nd dimension of model grid.
    !     *INTEGER* *KE*    - 3rd (vertical) dimension of model grid.
    !     *REAL*    *SAO*    - salinity [psu.].
    !     *REAL*    *THO*    - potential temperature [deg C].
    !     *REAL*    *ddpo*   - size of scalar grid cell (3rd dimension) [m].
    !     *REAL*    *tiestu* - level depths [m].
    !
    !     Externals
    !     ---------
    !     .
    !**********************************************************************
    ! Calculation of chemical constsnts in surface layer (chemcm)
    ! and in the water column (ak13, ak23, akb3, aksp)
    !--------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: ie,je,ke

    REAL(DP), INTENT(IN) :: ddpo(IE,JE,KE)
    REAL(DP), INTENT(IN) :: tho (IE,JE,KE)
    REAL(DP), INTENT(IN) :: sao (IE,JE,KE)
    REAL(DP), INTENT(IN) :: tiestu(KE+1)

    ! LOCAL
    INTEGER :: i,j,k
    REAL(DP):: ZERO,TENM7,SMICR,THOUSI,PERC,FOURTH,THIRD
    REAL(DP):: HALF,ONE,TWO,TEN
    REAL(DP):: devkb,devk1t,devk2t,devkbt,devkst,devks
    REAL(DP):: rgas,bor1,bor2,oxyco
    REAL(DP):: c00,c01,c02,c03,c04,c05,c10,c11,c12,c13,c20,c21,c22,c23
    REAL(DP):: cb0,cb1,cb2,cb3,cw0,cw1,cw2
    REAL(DP):: ox0,ox1,ox2,ox3,ox4,ox5,ox6
    REAL(DP):: an0,an1,an2,an3,an4,an5,an6
    REAL(DP):: akcc1, akcc2, akcc3, akcc4
    REAL(DP):: salchl, temzer, arafra, calfra,sucall,aracal
    REAL(DP):: devk1, devk2, t,q,s,cl
    REAL(DP):: cek0, ckb, ck1, ck2, ckw, oxy, ani
    REAL(DP):: ak1, ak2, akb, akw, ak0, aksp0
    REAL(DP):: rrr, rho,bor,p,cp,tc
    REAL(DP):: a1,a2,a3,b1,b2,b3,atn2o,rs

!mz_ap_20080114+
!      WRITE(*,*) "CHEMCON IS CALLED"
!mz_ap_20080114-

    !     -----------------------------------------------------------------
    !*         1. SET HALF PRECISION CONSTANTS
    !             --- ---- --------- ---------
    !
    ZERO=0._dp
    TENM7=10._dp**(-7.0_dp)
    SMICR=1.E-6_dp
    THOUSI=1._dp/1000._dp
    PERC=0.01_dp
    FOURTH=0.25_dp
    THIRD=1._dp/3._dp
    HALF=0.5_dp
    ONE=1._dp
    TWO=2._dp
    TEN=10._dp
    !
    !     -----------------------------------------------------------------
    !*         3. SET CONVERSION FACTOR SALINITY -> CHLORINITY
    !             ------ ------- -- ---- ------- -- ----------
    !             (AFTER WOOSTER ET AL., 1969)
    !
    SALCHL=1._dp/1.80655_dp
    !
    !     -----------------------------------------------------------------
    !*         4. SET ZERO DEG CENTIGRADE AT KELVIN SCALE
    !             --- ---- --- ---------- -- ------ -----
    !
    TEMZER=273.16_dp
    !
    !     -----------------------------------------------------------------
    !*         5. SET MEAN TOTAL [CA++] IN SEAWATER (MOLES/KG)
    !             (SEE BROECKER A. PENG, 1982, P. 26)
    !             ([CA++](MOLES/KG)=1.026E-2*(S/35.) AFTER
    !             CULKIN(1965), CF. BROECKER ET AL. 1982)
    !             ------------- --- -------- -- --- -----
    !
    CALCON=1.03E-2_dp
    !
    !     -----------------------------------------------------------------
    !*         6. SET COEFFICIENTS FOR APPARENT SOLUBILITY EQUILIBRIUM
    !             OF CALCITE (INGLE, 1800, EQ. 6)
    !             -- ------- ------- ----- --- ----------- -----------
    !
    AKCC1=-34.452_dp
    AKCC2=-39.866_dp
    AKCC3=110.21_dp
    AKCC4=-7.5752E-6_dp
    !
    !     -----------------------------------------------------------------
    !*        6A. SET FRACTION OF ARAGONITE IN BIOGENIC CACO3 PARTICLES
    !             --- -------- --------- -- -------- ----- ---------
    !
    ARAFRA=0._dp   ! js: obsolete
    !
    !     -----------------------------------------------------------------
    !*        6B. FRACTION OF CALCITE IN BIOGENIC CACO3 PARTICLES
    !             -------- ------- -- -------- ----- ---------
    !
    CALFRA=1._dp-ARAFRA
    SUCALL=ARAFRA+CALFRA

    !
    !     -----------------------------------------------------------------
    !*         7. FACTOR TO GET APPARENT SOLUBILITY PRODUCT FOR
    !             ARAGONIT BY MULTIPLICATION WITH APPARENT SOLUBILITY
    !             PRODUCT FOR CALCIT (CF. BERNER, 1976,
    !             OR BROECKER ET AL., 1982)
    !             -- -------- -- ---- ----- -- ------ --- ------- -----
    !
    ARACAL=ARAFRA*1.45_dp+CALFRA
    !     WRITE(io_stdo_bgc,*) 'ARACAL=',ARACAL
    !
    !     -----------------------------------------------------------------
    !*         8. SET COEFFICIENTS FOR SEAWATER PRESSURE CORRECTION OF

    !             (AFTER CULBERSON AND PYTKOWICZ, 1968, CF. BROECKER
    !             ET AL., 1982)
    !             ------- -------- --- ---------- ----- --- -------- -
    !
    DEVK1=24.2_dp
    DEVK2=16.4_dp
    DEVKB=27.5_dp
    DEVK1T=0.085_dp
    DEVK2T=0.04_dp
    DEVKBT=0.095_dp
    !
    !     -----------------------------------------------------------------
    !*         9. SET COEFFICIENTS FOR PRESSURE CORRECTION OF SOLUBILITY
    !             PRODUCT OF CACO3 CORRESPONDING TO ARAGONITE/CALCITE
    !             RATIO AFTER EDMOND AND GIESKES (1970),
    !             P. 1285
    !             -- ---- -- --------- ----- ------ --- ------- -------
    !
    DEVKST=0.23_dp
    DEVKS=32.8_dp*ARAFRA+35.4_dp*CALFRA
    !     WRITE(io_stdo_bgc,*) '***DEVKS=',DEVKS,' DEVKST=',DEVKST,TEN
    !
    !     -----------------------------------------------------------------
    !*        11. SET UNIVERSAL GAS CONSTANT
    !             --- --------- --- --------
    !
    RGAS=83.143_dp
    !
    !     -----------------------------------------------------------------
    !*        12. SET BORON CONCENTRATION IN SEA WATER
    !*            IN G/KG PER O/OO CL ACCORDING
    !             TO RILEY AND SKIRROW, 1965 (P.250)
    !             -- ----- --- -------- ---- ------- -
    !
    BOR1=0.00023_dp
    !
    !     -----------------------------------------------------------------
    !*        13. SET INVERSE OF ATOMIC WEIGHT OF BORON [G**-1]
    !             (USED TO CONVERT SPECIFIC TOTAL BORAT INTO CONCENTRATIONS)
    !             ----- -- ------- -------- ----- ----- ---- ---------------
    !
    BOR2=1._dp/10.82_dp
    !
    !     -----------------------------------------------------------------
    !*        14. SET INVERS OF NORMAL MOLAL VOLUME OF AN IDEAL GAS
    !             [CM**3]
    !             ---------- -- ------ ----- ------ -- -- ----- ---
    !
    OXYCO=1._dp/22414.4_dp
    !
    !     -----------------------------------------------------------------
    !*        15. SET VOLUMETRIC SOLUBILITY CONSTANTS FOR CO2 IN ML/L
    !             WEISS, R. F. (1974)
    !             CARBON DIOXIDE IN WATER AND SEAWATER: THE SOLUBILITY OF A
    !             NON IDEAL GAS. MARINE CHEMISTRY, VOL. 2, 203-215.
    !     -----------------------------------------------------------------

    C00=-58.0931_dp          ! C null null
    C01=90.5069_dp
    C02=22.2940_dp
    C03=0.027766_dp
    C04=-0.025888_dp
    C05=0.0050578_dp
    !
    !     -----------------------------------------------------------------
    !*        16. SET COEFF. FOR 1. DISSOC. OF CARBONIC ACID
    !             (EDMOND AND GIESKES, 1970)
    !             ------- --- -------- ------- -------- ----
    !
    C10=812.27_dp
    C11=3.356_dp
    C12=-0.00171_dp
    C13= 0.000091_dp

    !     -----------------------------------------------------------------
    !*        17. SET COEFF. FOR 2. DISSOC. OF CARBONIC ACID
    !             (EDMOND AND GIESKES, 1970)
    !             ------- --- -------- ------- -------- ----
    !
    C20=1450.87_dp
    C21=4.604_dp
    C22=-0.00385_dp
    C23= 0.000182_dp
    !
    !     -----------------------------------------------------------------
    !*        18. SET COEFF. FOR 1. DISSOC. OF BORIC ACID
    !             (EDMOND AND GIESKES, 1970)
    !             ------- --- -------- ------- ----- ----
    !
    CB0=2291.90_dp
    CB1=0.01756_dp
    CB2=-3.385_dp
    CB3=-0.32051_dp
    !
    !     -----------------------------------------------------------------
    !*        19. SET COEFF. FOR DISSOC. OF WATER
    !             (DICKSON AND RILEY, 1979, EQ. 7, COEFFICIENT
    !             CW2 CORRECTED FROM 0.9415 TO 0.09415 AFTER
    !             PERS. COMMUN. TO B. BACASTOW, 1988)
    !             ----- ------- -- -- --------- ------ -------
    !
    CW0=3441._dp
    CW1=2.241_dp
    CW2=-0.09415_dp
    !
    !     -----------------------------------------------------------------
    !*        20. SET VOLUMETRIC SOLUBILITY CONSTANTS FOR O2 IN ML/L
    !             (WEISS, 1970)
    !             ------- ------ --------- --------- --- -- -- ----
    !
    OX0=-173.4292_dp
    OX1=249.6339_dp
    OX2=143.3483_dp
    OX3=-21.8492_dp
    OX4=-0.033096_dp
    OX5=0.014259_dp
    OX6=-0.0017_dp

    !     -----------------------------------------------------------------
    !*            SET VOLUMETRIC SOLUBILITY CONSTANTS FOR N2 IN ML/L
    !             WEISS, R. F. (1970) THE SOLUBILITY OF NITROGEN
    !             OXYGEN AND ARGON IN WATER AND SEAWATER.
    !             DEEP-SEA RESEARCH, VOL. 17, 721-735.
    !     -----------------------------------------------------------------

    AN0=-172.4965_dp
    AN1=248.4262_dp
    AN2=143.0738_dp
    AN3=-21.7120_dp
    AN4=-0.049781_dp
    AN5=0.025018_dp
    AN6=-0.0034861_dp

    !      Constants for laughing gas solubility 
    !      (WEISS, 1974, MARINE CHEMISTRY)
    !      --------------------------------------  
    a1=-62.7062_dp
    a2=97.3066_dp
    a3=24.1406_dp
    b1=-0.058420_dp
    b2=0.033193_dp
    b3=-0.0051313_dp
    atn2o=3.e-7_dp
       
    ! mz_bk_20120708+
    rrrcl = salchl * 1.025_dp * bor1 * bor2
    ! mz_bk_20120708-
    !
    !     -----------------------------------------------------------------
    !*        21. CHEMICAL CONSTANTS - SURFACE LAYER
    !             -------- --------- - ------- -----
    !

    DO j=1,JE
       DO i=1,IE
          IF (ddpo(i,j,1) .GT. 0.5_dp) THEN                ! wet cell

             !
             !*        21.1 SET ABSOLUTE TEMPERATURE
             !              ------------------------
             T=THO(i,j,1)+TEMZER                     ! degC to K
             Q=T*PERC                                ! perc=0.01
             S=MAX(25._dp,SAO(i,j,1))                ! minimum salinity 25

             !      Laughing gas solubility (WEISS, 1974)
             !      --------------------------------------  
             rs=a1+a2*(100._dp/t)+a3*log(t/100._dp)                           &
                    &    +s*( b1 +b2*(t/100._dp) + b3*(t/100._dp)**2)

             satn2o(i,j)=atn2o*exp(rs)

             !
             !*        21.2 CHLORINITY (WOOSTER ET AL., 1969)
             !              ---------------------------------
             !
             CL=S*SALCHL
             !
             !*        21.3 LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1974)
             !              -------------------------------------------------

             CEK0=C00+C01/Q+C02*LOG(Q)+S*(C03+C04*Q+C05*Q**2)

             !
             !*        21.4 PK1, PK2 OF CARB. ACID, PKB OF BORIC ACID 
             !              -----------------------------------------
             !*             AFTER EDMOND AND GIESKES (1970)
             !              -------------------------------

             CKB=CB0/T+CB1*T+CB2+CB3*CL**THIRD
             CK1=C10/T+C11+C12*S*LOG(T)+C13*S**2
             CK2=C20/T+C21+C22*S*LOG(T)+C23*S**2

             !
             !*        21.5 CKW (H2O) (DICKSON AND RILEY, 1979)
             !              ------------------------------------

             CKW=CW0/T+CW1+CW2*SQRT(S)

             !
             !*****CKW COULD ADDITIONALLY BE EXPRESSED SALIN. DEPENDENT *****
             !

             !
             !*        21.6 LN(K0) OF SOLUBILITY OF O2 (EQ. 4, WEISS, 1970)
             !              -----------------------------------------------

             OXY=OX0+OX1/Q+OX2*LOG(Q)+OX3*Q+S*(OX4+OX5*Q+OX6*Q**2)

             !*      SOLUBILITY OF N2 
             !       WEISS, R. F. (1970), DEEP-SEA RESEARCH, VOL. 17, 721-735.
             !              -----------------------------------------------

             ANI=AN0+AN1/Q+AN2*LOG(Q)+AN3*Q+S*(AN4+AN5*Q+AN6*Q**2)

             !
             !*        21.7 K1, K2 OF CARB. ACID, KB OF BORIC ACID
             !              (EDMOND AND GIESKES,1970)
             !              -------------------------------------------------
             AK1=TEN**(-CK1)
             AK2=TEN**(-CK2)
             AKB=TEN**(-CKB)
             !
             !*        21.8 IONIC PRODUCT OF WATER KW (H2O)
             !              (DICKSON AND RILEY, 1979)
             !              -------------------------------------------------
             AKW=TEN**(-CKW)
             AKW3(I,J,1)=AKW
             !
             !*       21.9 CO2 SOLUBILITY IN SEAWATER
             !             (WEISS, 1974, CF. EQ. 12)
             !              -------------------------------------------------
             AK0=EXP(CEK0)*SMICR
             !
             !*       21.10 DENSITY OF SEAWATER AND TOTAL BORATE IN MOLES/L

             !      RRR=RHO(S,THO(i,j,1),ZERO) *THOUSI
             !      BOR=BOR1*RRR*CL*BOR2

             !      reformulation after R. Bacastow
             BOR=1.22e-5_dp*S 

             !
             !*       21.11 SET CHEMICAL CONSTANTS
             CHEMCM(i,j,5)=AK0
             CHEMCM(i,j,4)=ak1
             CHEMCM(i,j,3)=ak2
             CHEMCM(i,j,1)=akb
             CHEMCM(i,j,2)=AKW
             CHEMCM(i,j,6)=BOR
             !
             !*       21.12 O2/N2 SOLUBILITY IN SEAWATER (WEISS, 1970)
             !              -------------------------------------------------
             CHEMCM(i,j,7)=EXP(OXY)*OXYCO
             CHEMCM(i,j,8)=EXP(ANI)*OXYCO

             !      CHEMC(I,J,7,LAUMON)=EXP(OXY)*OXYCO/196800.
             !      CHEMC(I,J,8,LAUMON)=EXP(ani)*OXYCO/802000.

          ENDIF
       ENDDO
    ENDDO


    !
    !     -----------------------------------------------------------------
    !*        22. CHEMICAL CONSTANTS - DEEP OCEAN
    !        ----------------------------------------------------------------

    DO k=1,KE
       !
       !*        22.1 APPROX. SEAWATER PRESSURE AT U-POINT DEPTH (BAR)
       !              -------------------------------------------------------

       P=1.025E-1_dp*tiestu(k)
       !  WRITE(io_stdo_bgc,*) 'CHEMCON: P=1.025E-1*tiestu(k)', P,tiestu(k),k

       DO i=1,IE
          DO j=1,JE

             !
             !*        22.1.1 Zonal mean temperature/salinity
             !                -------------------------------
             !
             !      tzsum = 0.0
             !      szsum = 0.0
             !      vzsum = 0.0
             !      DO i=1,IE
             !         IF(ddpo(i,j,k).GT.0.5) THEN
             !            tzsum = tzsum+THO(i,j,k)*dlxp(i,j)*dlyp(i,j)*ddpo(i,j,k)
             !            szsum = szsum+SAO(i,j,k)*dlxp(i,j)*dlyp(i,j)*ddpo(i,j,k)
             !            vzsum = vzsum+            dlxp(i,j)*dlyp(i,j)*ddpo(i,j,k)
             !         ENDIF
             !      ENDDO
             !      IF(vzsum.GT.0.5) THEN
             !         tzmean = tzsum/vzsum
             !         szmean = szsum/vzsum
             !      ELSE
             !         tzmean = 0.0
             !         szmean = 0.0
             !      ENDIF
             !
             !
             !*        22.2 SET LIMITS FOR SEAWATER TEMP. AND SALINITY
             !              -------------------------------------------------
             !              (THIS IS DONE TO AVOID COMPUTATIONAL CRASH AT DRY
             !               POINTS DURING CALCULATION OF CHEMICAL CONSTANTS)
             !
             !*        22.3 SET [H+] (FIRST GUESS)
             !              -------------------------------------------------

             !
             !*        22.4 SET ABSOLUTE TEMPERATURE
             !              -------------------------------------------------
             t=THO(i,j,k)+273.16_dp
             q=t*perc
             ! mz_bk_20120708+
             ! s=MAX(34._dp,SAO(i,j,k))
             s=MAX(25._dp,SAO(i,j,k))
             ! mz_bk_20120708-

             !
             !*        22.5 CHLORINITY (WOOSTER ET AL., 1969)
             !              -------------------------------------------------
             CL=S*SALCHL
             !
             !*        22.6 LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1974)
             !              -------------------------------------------------
             CEK0=C00+C01/Q+C02*LOG(Q)+S*(C03+C04*Q+C05*Q**2)
             !
             !*        22.7 PK1, PK2 OF CARBONIC ACID, PKB OF BORIC ACID 
             !              -------------------------------------------------
             !              AFTER EDMOND AND GIESKES (1970)

             CKB=CB0/T+CB1*T+CB2+CB3*CL**THIRD
             !
             CK1=C10/T+C11+C12*S*LOG(T)+C13*S**2
             CK2=C20/T+C21+C22*S*LOG(T)+C23*S**2

             !*        22.8 LN(K0) OF SOLUBILITY OF O2 (EQ. 4, WEISS, 1970)
             !              -------------------------------------------------
             OXY=OX0+OX1/Q+OX2*LOG(Q)+OX3*Q+S*(OX4+OX5*Q+OX6*Q**2)

             satoxy(i,j,k)=exp(oxy)*oxyco
             !
             !*        22.9 K1, K2 OF CARBONIC ACID, KB OF BORIC ACID, KW
             !              (H2O) (LIT.?)
             AK1=TEN**(-CK1)
             AK2=TEN**(-CK2)
             AKB=TEN**(-CKB)
             !
             !*       22.10 APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE IN
             !              SEAWATER
             !              -------------------------------------------------
             !              (S=27-43, T=2-25 DEG C) AT P=0 (ATMOSPH. PRESSURE)
             !              (INGLE, 1800, EQ. 6)

             AKSP0=1.E-7_dp*(AKCC1+AKCC2*S**(THIRD)+AKCC3*LOG10(S)+AKCC4*T**2)
             CKW=CW0/T+CW1+CW2*SQRT(S)
             AKW3(I,J,K)=TEN**(-CKW)
             !
             !*       22.11 FORMULA FOR CP AFTER EDMOND AND GIESKES (1970)
             !              -------------------------------------------------

             !           (REFERENCE TO CULBERSON AND PYTKOQICZ (1968) AS MADE
             !           IN BROECKER ET AL. (1982) IS INCORRECT; HERE RGAS IS
             !           TAKEN TENFOLD TO CORRECT FOR THE NOTATION OF P IN
             !           DBAR INSTEAD OF BAR AND THE EXPRESSION FOR CP IS
             !           MULTIPLIED BY LN(10.) TO ALLOW USE OF EXP-FUNCTION
             !           WITH BASIS E IN THE FORMULA FOR AKSP (CF. EDMOND
             !           AND GIESKES (1970), P. 1285 AND P. 1286 (THE SMALL
             !           FORMULA ON P. 1286 IS RIGHT AND CONSISTENT WITH THE
             !           SIGN IN PARTIAL MOLAR VOLUME CHANGE AS SHOWN ON
             !           P. 1285))
             !    WRITE(io_stdo_bgc,*)  'CHEMCON: CP=P/(RGAS*T)', CP,P,RGAS,T 
             CP=P/(RGAS*T)
             !
             !*       22.12 KB OF BORIC ACID, K1,K2 OF CARBONIC ACID PRESSURE
             !              CORRECTION AFTER CULBERSON AND PYTKOWICZ (1968)
             !              (CF. BROECKER ET AL., 1982)

             TC=THO(i,j,k)

             !      WRITE(io_stdo_bgc,*)  ' CHEMCON: modal&(CP*(DEVKB-DEVKBT*TC)) ',modal,CP,DEVKB,DEVKBT,TC
      
             AKB3(I,J,K)=AKB*EXP(CP*(DEVKB-DEVKBT*TC))
             AK13(I,J,K)=AK1*EXP(CP*(DEVK1-DEVK1T*TC))
             AK23(I,J,K)=AK2*EXP(CP*(DEVK2-DEVK2T*TC))
             !
             !        22.13 APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE
             !              (OR ARAGONITE)
             !              -------------------------------------------------
             !              AS FUNCTION OF PRESSURE FOLLOWING EDMOND AND
             !              GIESKES (1970) (P. 1285) AND BERNER (1976)
             AKSP(I,J,K)=ARACAL*AKSP0*EXP(CP*(DEVKS-DEVKST*TC))



             !
             !*       22.14 DENSITY OF SEAWATER AND TOTAL BORATE CONCENTR.
             !              [MOLES/L]
             !              -------------------------------------------------

             !      rrr=rho(s,tc,p)*thousi                    
             !      bor=bor1*rrr*cl*bor2

             !     reformulation after R. Bacastow
             ! mz_bk_20120708+
             ! BOR = 1.22e-5_dp*S 
             ! mz_bk_20120708-

             ! mz_bk_20120708+
             ! rrrcl=salchl*1.025_dp*bor1*bor2
             ! mz_bk_20120708-



             !
             !     ----------------------------------------------------------
             !*        23. INITIATE [H+] AND [CO3--]
             !             -------- ---- --- -------
             !
             !*       23.1  FIRST GUESSES FOR [CO3--] AND [H+]
             !              -------------------------------------------------

             !      CARALK=ocetra(ialkali)%ptr(i,j,k)-BOR/(ONE+TENM7/AKB3(I,J,K))

             !
             !*       20.15 DENSITY OF SEAWATER AND TOTAL BORATE IN MOLES/L
             !              -------------------------------------------------
             !      RRR=RHO(S,TC,P)*THOUSI
             !      BOR=BOR1*RRR*CL*BOR2

             !     reformulation after R. Bacastow
             ! mz_bk_20120708+
             ! BOR=1.22e-5_dp*S 
             ! mz_bk_20120708-
             !
             !*       20.16 [CO3--], FIRST ESTIMATE
             !              -------------------------------------------------

          ENDDO
       ENDDO
    ENDDO

    !     -----------------------------------------------------------------
    !*        21. ITERATION TO INITIATE [H+] AND [CO3--]
    !             --------- -- -------- ---- --- -------
    !      BT=BOR
    !      BREMS=1.5
    !      DO 43 KI=1,30
    !         zer(ki)=0.
    !         DO 43 k=1,KE
    !         DO 43 j=1,JE
    !         DO 43 i=1,IE
    !
    !            IF(ddpo(i,j,k).GT.0.5) THEN
    !
    !            ak1=ak13(i,j,k)
    !            ak2=ak23(i,j,k)
    !            H=HI(i,j,k)
    !            R=CO3(i,j,k)
    !            ALKA=ocetra(ialkali)%ptr(i,j,k)
    !            C=ocetra(isco212)%ptr(i,j,k)
    !            T1=h/ak13(i,j,k)
    !            T2=h/ak23(i,j,k)
    !            AKW=AKW3(J,K)
    !            BT=rrrcl*SAO(i,j,k)
    !            AKB=AKB3(J,K)
    !            alk=ocetra(ialkali)%ptr(i,j,k)
    !            A=!*(2.+t2)/(1.+t2+t2*t1)  +AKW/H-H+BT/(1.+H/AKB)-ALK
    !            A=!*(2.+t2)/(1.+t2+t2*t1)  +AKW/H-H+BT/(1.+H/AKB)-ALK
    !            zer(ki)=zer(ki)+a**2
    !            DADH=!*(1./(AK2*(1.+T2+T2*T1))-(2.+T2)*(1./AK2+2.*T1/AK2)/
    !     1          (1.+T2+T2*T1)**2)
    !     1          -AKW/H**2-1.-(BT/AKB)/(1.+H/AKB)**2
    !            dddhhh=a/dadh
    !            reduk=MAX(1._dp,1.2*abs(dddhhh/h))
    !            H=H-dddhhh/reduk
    !  
    !            HI(i,j,k)=HI(i,j,k)-dddhhh
    !            co3(i,j,k)
    !     1      =c/(1.+hi(i,j,k)*(1.+hi(i,j,k)/ak13(i,j,k))/ak23(i,j,k))
    !
    !           ENDIF
    !
    !43      CONTINUE

  END SUBROUTINE hamocc_chemcon
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE hamocc_ocprod(ie,je,ke,tho,ddpo,dpio)
    !
    ! mz_bk_20120626 this subroutine was updated on 26 Jun 2012 to match
    ! hamocc version of mpiesm-1.0.00
    !
    ! Correspondes to subroutine ocprod in ocprod.f90
    ! Calculates biological production, settling of debris,
    ! and related biogeochemistry
    !
    ! Original by:
    ! Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    !     Modified
    !     --------
    !     S.Legutke,                                *MPI-MaD, HH*    10.04.01
    !     S.Lorenz/JO.Beismann, OpenMP parallel     *MPI-Met, HH*    24.08.07
    !
    ! surface layer = 1.
    !------------------------------------------------------------------------
    implicit none

    INTEGER, INTENT(IN) :: ie, je, ke
    REAL(DP),INTENT(IN) :: tho (ie,je,ke)
    REAL(DP),INTENT(IN) :: ddpo(ie,je,ke)
    REAL(DP),INTENT(IN) :: dpio(ie,je,ke)
    REAL(DP):: abs_bgc(ie,je,ke)           !js name implies "absorbtion bgc",
    !                                      !   but it is the light fraction
    ! note: in mpiesm-1.0.00 abs_bgc is swr_frac, which is calculated in
    ! mpiom mo_swr_absorption.f90

    INTEGER :: i, j, k, l
    REAL(DP):: dmsp1, dmsp2, dmsp3, dmsp4, dmsp5, dmsp6
    REAL(DP):: atten, avphy, avanut, avanfe, pho, xa, xn, ya, yn, phosy,      &
         &     volcell, avgra, grazing, avsil, graton,                        &
         &     gratpoc, grawa, bacfra, phymor, zoomor, excdoc, exud,          &
         &     export, delsil, delcar, sterph, sterzo, remin,                 &
         &     docrem, opalrem, remin2o, aou, refra
    ! mz_bk_20120210+
    REAL(DP):: phofa, temfa
    ! mz_bk_20120210-
    ! mz_bk_20120626+
    REAL(DP):: zoothresh
    ! mz_bk_20120626-
    ! mz_bk_20120626+
    REAL(DP) :: wopal, wpoc, wcal
    ! mz_bk_20120626-
    ! mz_bk_20120626+
    ! REAL(DP):: epfac, tloc, eppz15    ! not used
    ! mz_bk_20120626-
  
    REAL(DP):: dustinp, dustinp_volcano
    ! mz_bk_20120626+
    ! REAL(DP):: gutscale
    ! mz_bk_20120626-
    REAL(DP):: fopa, fdet, fcal
    REAL(DP):: absorption
    REAL(DP):: dmsprod, dms_bac, dms_uv
    REAL(DP):: detref, detrl
    INTEGER :: volcano
!!$!! #ifdef AGG
    REAL(DP):: wmass(ie,je,ke)     ! sinking speed for 'mass of aggregates'
    REAL(DP):: wnumb(ie,je,ke)     ! sinking speed for 'numbers of aggregates'
    REAL(DP):: aggregate(ie,je,ke) ! aggregation (should be renamed)
    REAL(DP):: dustagg(ie,je,ke)   ! aggregation of dust

    REAL(DP):: avmass, avnos, anosloss
    ! nosin = nos increase, nosde = nos decrease
    REAL(DP):: zmornos, avm, avn, eps, e1, e2, e3, e4, es1, es3
    REAL(DP):: TopM, TopF, fshear, sagg1, sagg2, sagg4
    REAL(DP):: sett_agg, shear_agg, effsti, dfirst, dshagg, dsett
    REAL(DP):: checksize, nacheck, flar1, flar2, flar3
    REAL(DP):: fTSFac, fTMFac, fTopF, fTopM, wphy, wphyup, wnosup, wnos
!!$!! #endif 

    volcano=0           ! volcano = 0: no dust from volcano,
    !                   ! volcano = 1: ashes from volcano in dustinp_volcano
    !
    ! Constant parameters
    !
    ! parameter definition in beleg_bgc.f90

    dmsp6 = dmspar(6)
    dmsp5 = dmspar(5)
    dmsp4 = dmspar(4)
    dmsp3 = dmspar(3)
    dmsp2 = dmspar(2)
    dmsp1 = dmspar(1)

#if defined MPIOM_13B
    ! Calculate swr absorption by water and phytoplankton

    ! Almost half of the SWR is absorbed in the surface layer (0.04*10m).
    ! --> *0.6 is transmitted (js: for atten_c=0!)
    ! 1) Upper 2 layers
    DO j=1,je
       DO i=1,ie
          abs_bgc(i,j,1)=1._dp        ! surface layer has full light available
          abs_oce(i,j,1)=1._dp
          IF (ddpo(i,j,1) .GT. 0.5_dp) THEN
             atten = atten_w + atten_c * ocetra(iphy)%ptr(i,j,1)
             !                           ! atten_w = 0.04 m^-1, 
             !                           ! atten_c=7.32e5 (defined in beleg_bgc)
             !mz_ap_20080118+
             ! This atten_c*ocetra(iphy)%ptr(i,j,1) = 0.03/m [Chl]
             ! ([Chl]=1/80[phy]) see Wetzel et al,2006) 
             !mz_ap_20080118-
             ! for ocetra(iphy)=5e-7 attenuation is 0.36/m : seems too much
             absorption = exp(-atten*ddpo(i,j,1))   ! absorption in surface
             !                                      ! layer: 0.67 for atten_c=0,
             !                                      ! sfc layer 10m
             abs_bgc(i,j,2) = absorption            ! second layer 
             ! (Implicit would be faster:)
             !  abs_bgc(i,j,2)=1./(1+atten*ddpo(i,j,1))

             abs_oce(i,j,2)=atten_f*absorption    ! atten_f =0.4 (beleg_bgc.f90)
             !  abs_oce(i,j,2)=0.5/(1+atten*ddpo(i,j,1))

          ENDIF
       ENDDO
    ENDDO
      
    ! 2) water column
    DO k=3,ke          ! (js) why run over whole water colum? kwrbioz sufficient
       DO j=1,je
          DO i=1,ie
             IF (ddpo(i,j,k-1) .GT. 0.5_dp) THEN      ! wet point
                !  atten=0.04+0.03*2.44e7*ocetra(iphy)%ptr(i,j,k-1)
                !  abs_bgc(i,j,k)=abs_bgc(i,j,k-1)/(1+atten*ddpo(i,j,k-1))
                atten = atten_w + atten_c * ocetra(iphy)%ptr(i,j,k-1)
                absorption = exp(-atten * ddpo(i,j,k-1))
                abs_bgc(i,j,k) = abs_bgc(i,j,k-1) * absorption
                !  abs_oce(i,j,k)=abs_oce(i,j,k-1)/(1+atten*ddpo(i,j,k-1))
                abs_oce(i,j,k) = abs_oce(i,j,k-1) * absorption
             ENDIF
          ENDDO
       ENDDO
    ENDDO
#elif defined MPIOM_2000
    !NB in MPIOM_2000 abs_oce => swr_frac!!!
    abs_oce(:,:,1) = 1.0_dp

    DO k=2,ke
       DO j=1,je
          DO i=1,ie
             IF (ddpo(i,j,k-1) .GT. 0.5_dp) THEN      ! wet point
                abs_oce(i,j,k) = abs_oce(i,j,k-1) * (                         &
                     !!mz_ap_20100903+
                     !  & redfrac * EXP(-dzw(k-1) * atten_r)                  &
                     !  & + (1.0_wp-redfrac) * EXP(-dzw(k-1) * (atten_w +     &
                     & redfrac * EXP(-ddpo(i,j,k-1) *  atten_r)               &
                     & + (1.0_dp-redfrac) * EXP(-ddpo(i,j,k-1) * (atten_w +   &
                     !mz_ap_20100903-
                     & atten_c * pho_to_chl                                   &
                     & * MAX(0.0_dp, ocetra(iphy)%ptr(i,j,k-1)))))
             ENDIF
             abs_bgc(i,j,k) = abs_oce(i,j,k)
          END DO
       END DO
    END DO
#endif

    ! mz_bk_20120626+
    wopal = sinkspeed_opal * dtb
    wcal = sinkspeed_cal * dtb
    wpoc = sinkspeed_poc * dtb
    ! mz_bk_20120626-

    ! dust flux from the atmosphere to the surface layer; 
    ! dust fields are monthly mean values in units of kg/m2/year
    ! dissolved iron is a fixed fraction (3.5%), and immediately released
    ! 1% of the iron input is dissolved [='bio-available']
    ! (see perc_diron in beleg_bgc)

    dustinp_volcano=0._dp

    do j=1,je
       do i=1,ie
          if (ddpo(i,j,1) .gt. 0.5_dp) then   ! wet point
             dustinp = dusty(i,j) / 365._dp * dtb * dpio(i,j,1)
             ocetra(ifdust)%ptr(i,j,1) = ocetra(ifdust)%ptr(i,j,1)            &
                  &                    + dustinp + dustinp_volcano 
             ocetra(iiron)%ptr(i,j,1)  = ocetra(iiron)%ptr(i,j,1)             &
                  &                    + dustinp * perc_diron                 &
                  ! perc_diron =.035*.01 /55.85 ! why *0.8 in IPCC_HAM?
                  &                    + dustinp_volcano * perc_diron         &
                  &                    * 0.8_dp / 3.5_dp  
             ! volcanic ash has lower iron content
             !mz_ap_20090727+
             ! we do the same for silica following Cotrim da Cunha, GBC, 2007 
             ! i.e. the PISCES model
             ! n.b. dustinp_volcano = 0!
             ocetra(isilica)%ptr(i,j,1) = ocetra(isilica)%ptr(i,j,1)          &
                  &                     + dustinp * 0.308_dp * 0.075_dp / 60._dp
             !                                                      (60: g->mol)
             !mz_ap_20090727-
          endif
       enddo
    enddo

!!$!! #ifdef AGG
    IF (L_AGG) THEN
       !***********************************************************************
       !
       ! special resetting for particle numbers, that sets their concentra-
       ! tion (particles per volume, ocetra(inos)) depending on the mass of
       ! marine snow:
       !
       ! Compartments have already been set to 0 in 
       ! ADVECTION_BGC.h and OCTDIFF_BGC.h. js: ???
       !
       ! Ensure that if there is no mass, there are no particles, and 
       ! that the number of particles is in the right range (this is crude,
       ! but is supposed to be needed only due to numerical errors such as
       ! truncation or overshoots during advection)
       !
       ! (1) avnos<<avmass, such that eps = FractDim + 1: increase numbers
       !     such that eps = FractDim + 1 + safe (currently set to 1.e-6 in
       !     BELEG_BGC) 
       !
       ! (2) avnos>>avmass, such that  Nbar (=Mass/Nos/cellmass) <=1:
       !     decrease numbers such that Nbar=1.1 (i.e. 1.1 cells per
       !     aggregate, set in BELEG_BGC) 

       DO k=1,ke
          DO j=1,je
             DO i=1,ie
                IF (ddpo(i,j,k) .GT. 0.5_dp) THEN               ! wet cell
                   avmass = ocetra(iphy)%ptr(i,j,k) + ocetra(idet)%ptr(i,j,k)
                   snow = avmass * 1.e+6_dp                     ! why *1.e6??
                   ! check whether the numbers have to be decreased
                   ! or increased     
                   ocetra(inos)%ptr(i,j,k) =                                  &
                        &             MAX(snow*pupper, ocetra(inos)%ptr(i,j,k)) 
                   ocetra(inos)%ptr(i,j,k) =                                  &
                        &             MIN(snow*plower, ocetra(inos)%ptr(i,j,k))
                   !js (MAX/MIN correct?)
                ENDIF     ! endif wet cell
             ENDDO
          ENDDO
       ENDDO

    ENDIF !if L_AGG
!!$!! #endif  /*AGG*/


    !
    ! Biological productivity in the euphotic zone (upper 90m)
    !
    DO k=1,kwrbioz

       DO j=1,je
          DO i=1,ie

             IF (ddpo(i,j,k) .GT. 0.5_dp) THEN
!!$!! #ifdef AGG
                IF (L_AGG) THEN
                   avmass = ocetra(iphy)%ptr(i,j,k) + ocetra(idet)%ptr(i,j,k)
                ENDIF
!!$!! #endif /*AGG*/

                avphy = MAX(phytomi,ocetra(iphy)%ptr(i,j,k))    ! 'available'
                !                                               ! phytoplankton
                avgra = MAX(grami,ocetra(izoo)%ptr(i,j,k))      ! 'available'
                !                                               ! zooplankton
                avsil = MAX(0._dp,ocetra(isilica)%ptr(i,j,k))   ! available
                !                                               ! silicate
                avanut = MAX(0._dp,MIN(ocetra(iphosph)%ptr(i,j,k),            &
                     !            !available nutrients (phosphate   [kmol P /m3]
                     & rnoi * ocetra(iano3)%ptr(i,j,k)))        !     + nitrate)
                avanfe = MAX(0._dp,MIN(avanut,ocetra(iiron)%ptr(i,j,k)/riron))
                !                                               ! available iron
                ! mz_bk_20120606+
                ! mark limiting element in bgc diagnostics
                IF (L_BGC_DIAG) THEN
                   IF (avanfe .EQ. ocetra(iiron)%ptr(i,j,k) / riron) THEN
                      bgcprod(kflim)%ptr(i,j,k) = 1._dp
                      bgcprod(knlim)%ptr(i,j,k) = 0._dp
                      bgcprod(kplim)%ptr(i,j,k) = 0._dp
                   ELSE IF (avanfe .EQ. avanut) THEN
                      IF (ocetra(iphosph)%ptr(i,j,k) .LE.                     &
                           &    rnoi*ocetra(iano3)%ptr(i,j,k)) THEN
                         bgcprod(kplim)%ptr(i,j,k) = 1._dp
                         bgcprod(knlim)%ptr(i,j,k) = 0._dp
                         bgcprod(kflim)%ptr(i,j,k) = 0._dp
                      ELSE
                         bgcprod(kplim)%ptr(i,j,k) = 0._dp
                         bgcprod(knlim)%ptr(i,j,k) = 1._dp
                         bgcprod(kflim)%ptr(i,j,k) = 0._dp
                      END IF
                   END IF
                ENDIF   ! logical end of bgc_diag
                ! mz_bk_20120606-

                ! mz_bk_20120210+
                ! mz_bk_20120606+
                ! here fPAR is missing (or should be included in pi_alpha,
                ! but is not!); so DON'T USE the following 2 lines, use the
                ! 2 lines below that instead!
                ! dtb is included in the following equations; as in
                ! mpiesm-1.0.00 it is not included in pi_alpha!
                ! !pho=pi_alpha*dtb*strahl(i,j)*abs_bgc(i,j,k)                &
                ! !                                     ! biological production
                ! !   & *(1.+0.06*tho(i,j,k)*(1.+0.03*tho(i,j,k)))
                ! temperature dependency from ernst, 0.06^T -> ln exp...
                ! pho=pi_alpha*dtb*fPAR*strahl(i,j)*abs_bgc(i,j,k)            &
                ! !                                     ! biological production
                !    & *(1.+0.06*tho(i,j,k)*(1.+0.03*tho(i,j,k)))
                ! !       temperature dependency from ernst, 0.06^T -> ln exp...
                ! mz_bk_20120606-
                ! mz_bk_20120210-

                ! mz_bk_20120210+
                ! production calculation updated from 'new HAMOCC'
                ! (mpiom_2000!)
                ! mz_bk_20120515 : included fPAR according HAMOCC
                ! in mpiesm-1.0.0
                ! in abs_oce factor of 0.4 already included? see above
                phofa = pi_alpha * fPAR * strahl(i,j) * abs_bgc(i,j,k)
                !                      ! * abs_bgc: light absorbtion coefficient
                temfa= 0.6_dp * 1.066_dp**tho(i,j,k)
                pho= dtb * phofa * temfa / SQRT(phofa**2 + temfa**2)
                ! mz_bk_20120210-

                xa = avanfe
                xn = xa / (1._dp + pho * avphy / (xa + bkphy))
                !                             ! bkphy = half saturation constant
                phosy = MAX(0._dp, xa-xn)     ! photo synthesis
                xn = MAX(xn, 1.e-10_dp)
                ya = avphy + phosy            ! new phytoplankton concentration
                !                             ! before grazing
                yn = (ya + grazra_dtb*avgra*phytomi/(avphy + bkzoo)) & ! grazing
                     &  / (1._dp + grazra_dtb*avgra / (avphy + bkzoo))
                grazing = MAX(0._dp, ya-yn)          ! what about grazing below
                !                                    ! euphotic zone?
                graton = epsher * (1._dp - zinges) * grazing    ! "grazing to
                !                                    ! (re-dissolved) nutrients"
                gratpoc = (1._dp - epsher) * grazing            ! epsher=0.8
                !                                             ! "grazing to POC"
                grawa = epsher * zinges * grazing         ! grazer 'wachstum(?)'

                ! mz_bk_20120626+
                !                         remineralization of poc12 using oxygen
                IF (ocetra(ioxygen)%ptr(i,j,k) .gt. 5.e-8_dp) THEN
                   remin = MIN(drempoc_dtb * ocetra(idet)%ptr(i,j,k),         &
                        &  0.5_dp * ocetra(ioxygen)%ptr(i,j,k)/ro2ut)
                   !                    ! 'detritus remineralized fraction' (?)
                   detref = remin / (ocetra(idet)%ptr(i,j,k) + 1.e-20_dp)
                ELSE
                   remin = 0._dp
                   ! mz_bk_20120626 detref never used, so it may be ok,
                   ! not to set it to 0._dp here
                ENDIF
                ! mz_bk_20120626-

                bacfra = remido_dtb * ocetra(idoc)%ptr(i,j,k)
                !      ! remido_dtb = remineralization rate of DOM
                phymor = dyphy_dtb * MAX(0._dp,                               &
                     & (ocetra(iphy)%ptr(i,j,k) - 2._dp * phytomi))         
                !              ! phytoplankton mortality dyphy_dtb = .008 * dtb

                ! mz_bk_20120626+
                zoothresh = MAX(0._dp,                                        &
                     &      (ocetra(izoo)%ptr(i,j,k) - 2._dp * grami))
                zoomor = spemor_dtb*zoothresh*zoothresh  ! zooplankton mortality
                excdoc = gammaz_dtb * zoothresh ! excretion to DOC (zooplankton)
                exud = gammap_dtb * MAX(0._dp,                                &
                     (ocetra(iphy)%ptr(i,j,k) - 2._dp * phytomi))
                !                             ! exudation to DOC (phytoplankton)

                !zoomor=spemor_dtb*MAX(0._dp,(ocetra(izoo)%ptr(i,j,k)-2.*grami))
                ! zooplankton mortality  
                !  excdoc=gammaz_dtb*MAX(0._dp,(ocetra(izoo)%ptr(i,j,k)-2.*grami))
                ! excretion to DOC (zooplankton)
                !  exud=gammap_dtb*MAX(0._dp,(ocetra(iphy)%ptr(i,j,k)-2.*phytomi))
                ! exudation to DOC (phytoplankton)
                ! mz_bk_20120626-

                ! mz_bk_20120626+
                ocetra(iphosph)%ptr(i,j,k) =                                  &
                     & ocetra(iphosph)%ptr(i,j,k) + bacfra - phosy + graton   &
                     & + ecan * zoomor + remin
                ocetra(iano3)%ptr(i,j,k) =                                    &
                     & ocetra(iano3)%ptr(i,j,k) + (bacfra - phosy + graton    &
                     & + ecan * zoomor) * rnit + remin * rnit
                ! ocetra(iphosph)%ptr(i,j,k) =                                &
                !     & ocetra(iphosph)%ptr(i,j,k) + bacfra - phosy           &
                !     & + graton + ecan * zoomor
                ! ocetra(iano3)%ptr(i,j,k) =                                  &
                !     & ocetra(iano3)%ptr(i,j,k) + (bacfra - phosy            &
                !     & + graton + ecan * zoomor) * rnit
                ! mz_bk_20120626-

                export = zoomor * (1._dp - ecan) + phymor + gratpoc
                !                      ! ecan=.95, gratpoc= .2*grazing [P-units]

                ! mz_bk_20120626+
                ocetra(idet)%ptr(i,j,k) = ocetra(idet)%ptr(i,j,k) + export    &
                     &                  - remin                          ! k=1,8
                ! ocetra(idet)%ptr(i,j,k) = ocetra(idet)%ptr(i,j,k) + export
                !                                                          k=1,8
                ! mz_bk_20120626-

!!$!!#ifdef AGG       
                IF (L_AGG) THEN
                   delsil = MIN(ropal*phosy*avsil/(avsil+bkopal), 0.5_dp*avsil)
                   delcar = rcalc * MIN(calmax*phosy, (phosy - delsil/ropal))
                ELSE
!!$!!#else
                   delsil = MIN(ropal*export*avsil/(avsil+bkopal),            &
                        &       0.5_dp*avsil)      
                   delcar = rcalc * export * bkopal/(avsil+bkopal)          
                   !             ! 'detritus linked calcium carbonate ' ?P units
                ENDIF
!!$!!#endif

                ! mz_bk_20120707+
                ! Added factor dtb as test, as otherwise we get 1/d ??
                ! ! DMS (js: slightly out of place)
                dmsprod = (dmsp5 * delsil + dmsp4 * delcar)                 &
                     &  * (1._dp + 1._dp / (tho(i,j,k) + dmsp1)**2)
                dms_bac = dmsp3 * abs(tho(i,j,k) + 3._dp)                   &
                     &  * ocetra(idms)%ptr(i,j,k)    & ! bacterial consumption
                     &  * (ocetra(idms)%ptr(i,j,k)   &
                     &  / (dmsp6 + ocetra(idms)%ptr(i,j,k)))
                ! DMS (js: slightly out of place)
                ! dmsprod = dtb * (dmsp5 * delsil + dmsp4 * delcar)           &
                !      &  * (1._dp + 1._dp / (tho(i,j,k) + dmsp1)**2)
                ! dms_bac = dtb * dmsp3 * abs(tho(i,j,k) + 3._dp)             &
                !      &  * ocetra(idms)%ptr(i,j,k)    & ! bacterial consumption
                !      &  * (ocetra(idms)%ptr(i,j,k)   &
                !      &  / (dmsp6 + ocetra(idms)%ptr(i,j,k)))
                ! mz_bk_20120707-

                ! mz_bk_20120210+
                ! dms_uv  = dmsp2*4.*pho*ocetra(idms)%ptr(i,j,k)
                !                                      decay due to UV-radiation
                ! mz_bk_20120210-

                ! mz_bk_20120210+
                dms_uv  = dmsp2 * 4._dp * dtb * phofa                         &
                     &  * ocetra(idms)%ptr(i,j,k)    ! decay due to UV-radiation
                ! mz_bk_20120210-

                ocetra(idms)%ptr(i,j,k) = ocetra(idms)%ptr(i,j,k)             &
                     &                  + dmsprod - dms_bac - dms_uv
                ! end DMS

                ! mz_bk_20120626+
                ocetra(isco212)%ptr(i,j,k) = ocetra(isco212)%ptr(i,j,k)       &
                     &                     - delcar      &  ! - CACO3 production
                     &                     + rcar * ( bacfra - phosy          &
                     &                     + graton + ecan * zoomor + remin)
                !                                     + remineralization C-units
                ocetra(ialkali)%ptr(i,j,k) = ocetra(ialkali)%ptr(i,j,k)       &
                     &                     - 2._dp * delcar                   &
                     &                     - rnit * ( bacfra - phosy          &
                     &                     + graton + ecan * zoomor + remin)
                ! ocetra(isco212)%ptr(i,j,k)=ocetra(isco212)%ptr(i,j,k)       &
                !      & - delcar                        &  ! - CACO3 production
                !      & + rcar * ( bacfra - phosy + graton + ecan * zoomor)
                ! !                                   + remineralization C-units
                ! ocetra(ialkali)%ptr(i,j,k)=ocetra(ialkali)%ptr(i,j,k)       &
                !      & - 2._dp * delcar                                     &
                !      & - rnit * ( bacfra - phosy + graton + ecan * zoomor)
                ! mz_bk_20120626-

                ocetra(iphy)%ptr(i,j,k) = ocetra(iphy)%ptr(i,j,k) + phosy     &
                     &                  - grazing - phymor - exud

                ! mz_bk_20120626+
                ocetra(ioxygen)%ptr(i,j,k) = ocetra(ioxygen)%ptr(i,j,k)       &
                     &                     + ro2ut * (phosy - bacfra)         &
                     &                     - (graton + ecan * zoomor          &
                     &                     + remin) * ro2ut
                ! ocetra(ioxygen)%ptr(i,j,k) = ocetra(ioxygen)%ptr(i,j,k)     &
                !    &                     + ro2ut * (phosy - bacfra)         &
                !    &                     - (graton + ecan * zoomor) * ro2ut
                ! mz_bk_20120626-

                ocetra(izoo)%ptr(i,j,k) = ocetra(izoo)%ptr(i,j,k) + grawa     &
                     &                  - excdoc - zoomor
                ocetra(idoc)%ptr(i,j,k) = ocetra(idoc)%ptr(i,j,k) - bacfra    &
                     &                  + excdoc + exud
                ocetra(icalc)%ptr(i,j,k)=ocetra(icalc)%ptr(i,j,k) + delcar

                ! mz_bk_20120626+
                opalrem = dremopal_dtb * 0.1_dp                                   &
                     &  * (tho(i, j, k) + 3.0_dp) * ocetra(iopal)%ptr(i,j,k)
                ocetra(isilica)%ptr(i,j,k) = ocetra(isilica)%ptr(i,j,k)       &
                     &                     - delsil + opalrem
                ocetra(iopal)%ptr(i,j,k) = ocetra(iopal)%ptr(i,j,k)           &
                     &                   + delsil - opalrem
                ocetra(iiron)%ptr(i,j,k) = ocetra(iiron)%ptr(i,j,k)           &
                     &                   + (bacfra - phosy + graton           &
                     &                   + ecan * zoomor + remin) * riron     &
                     &                   - relaxfe_dtb                            &
                     &          * MAX(ocetra(iiron)%ptr(i,j,k) - fesoly, 0._dp)
                ! ocetra(isilica)%ptr(i,j,k) = ocetra(isilica)%ptr(i,j,k)     &
                !      &             - delsil                                 &
                !      &             + dremopal_dtb * ocetra(iopal)%ptr(i,j,k)
                ! ocetra(iopal)%ptr(i,j,k) = ocetra(iopal)%ptr(i,j,k)         &
                !      &             + delsil                                 &
                !      &             - dremopal_dtb * ocetra(iopal)%ptr(i,j,k)
                ! ocetra(iiron)%ptr(i,j,k) = ocetra(iiron)%ptr(i,j,k)         &
                !      &  +(bacfra-phosy+graton+ecan*zoomor)*riron            &
                !      &- relaxfe_dtb*MAX(ocetra(iiron)%ptr(i,j,k)-fesoly,0._dp)
                ! mz_bk_20120626-

!!$!!  #ifdef AGG
                IF (L_AGG) THEN
                   !*********************************************************
                   ! effects of biological processes on number of particles:
                   ! photosynthesis creates POM
                   ! exudation removes POM
                   ! grazing removes POM; but only the fraction that is not
                   ! egested as fecal pellets again (grawa remains in zoo,
                   ! graton goes to po4)
                   ! none of the processes at the current time is assumed
                   ! to change the size distribution (subject to change)
                   ! NOTE that phosy, exud etc. are in kmol/m3! 
                   ! Thus divide by avmass (kmol/m3)
                   !**********************************************************
                   if (avmass .gt. 0._dp) then
                      avnos = ocetra(inos)%ptr(i,j,k) 
                      anosloss = (phosy-exud-graton-grawa)*avnos / avmass
                      ocetra(inos)%ptr(i,j,k)=ocetra(inos)%ptr(i,j,k)+anosloss
                   endif

                   !*********************************************************
                   ! dead zooplankton corpses come with their own, flat
                   ! distribution
                   ! this flow even takes place if there is neither nos nor
                   ! mass
                   ! NOTE: zoomor is in kmol/m3!! Thus multiply flow by 1.e+6
                   !*********************************************************
                   zmornos = zoomor * (1._dp - ecan) * zdis * 1.e+6_dp
                   ocetra(inos)%ptr(i,j,k) = ocetra(inos)%ptr(i,j,k) + zmornos

                ENDIF
!!$!!  #endif /*AGG*/

                !
                ! write output for bgcmean (sum over kwrbioz for 2d fields)
                !                          (should that be at the bottom of
                !                           euphotic zone? (at least for
                !                           export))
                !

                ! mz_bk_20120606+
                ! bgc diagnostics output
                IF (L_BGC_DIAG) THEN
                   bgcprod(kdmsprod)%ptr(i,j,k)  = dmsprod/dtbgc
                   bgcprod(kdms_bac)%ptr(i,j,k)  = dms_bac/dtbgc
                   bgcprod(kdms_uv)%ptr(i,j,k)   = dms_uv/dtbgc
                   bgcprod(kexport)%ptr(i,j,k)   = export*rcar/dtbgc
                   bgcprod(kdelcar)%ptr(i,j,k)   = delcar/dtbgc
                   bgcprod(kdelsil)%ptr(i,j,k)   = delsil/dtbgc
                   bgcprod(kphosy)%ptr(i,j,k)    = phosy/dtbgc
                   bgcprod(kgraz)%ptr(i,j,k)     = grazing/dtbgc
                ENDIF
                ! mz_bk_20120606-

             ENDIF      ! ddpo(i,j,k).GT.0.5
          ENDDO   ! ie
       ENDDO   ! je
    ENDDO   ! kwrbioz

!!$!! #ifdef AGG
    IF (L_AGG) THEN
       DO  k=1,ke
          DO j=1,je
             DO i=1,ie
                IF(ddpo(i,j,k) .GT. 0.5_dp) THEN
                   avmass = ocetra(iphy)%ptr(i,j,k) + ocetra(idet)%ptr(i,j,k)
                   snow = avmass * 1.e+6_dp
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!!$!! #endif /*AGG*/


    !-----below euphotic zone
    DO k=kwrbioz+1,ke
       DO j=1,je
          DO i=1,ie
             IF(ddpo(i,j,k) .GT. 0.5_dp) THEN
!!$!! #ifdef AGG
                IF (L_AGG) THEN
                   avmass = ocetra(iphy)%ptr(i,j,k) + ocetra(idet)%ptr(i,j,k)
                ENDIF
!!$!! #endif /*AGG*/          
                ! 'sterberate' phytoplankton
                sterph = dphymor_dtb                                          &
                     & * MAX(0._dp, ocetra(iphy)%ptr(i,j,k) - phytomi)
                ! 'sterberate' zooplankton
                sterzo = dzoomor_dtb*MAX(0._dp, ocetra(izoo)%ptr(i,j,k) - grami)
                ocetra(iphy)%ptr(i,j,k) = ocetra(iphy)%ptr(i,j,k) - sterph   
                ocetra(izoo)%ptr(i,j,k) = ocetra(izoo)%ptr(i,j,k) - sterzo   
                ! remineralization of poc12 using oxygen
                IF (ocetra(ioxygen)%ptr(i,j,k) .gt. 5.e-8_dp) THEN
                   remin = MIN(drempoc_dtb * ocetra(idet)%ptr(i,j,k),         &
                        &      0.5_dp * ocetra(ioxygen)%ptr(i,j,k) / ro2ut)
                   ! 'detritus remineralized fraction' (?)
                   detref = remin / (ocetra(idet)%ptr(i,j,k) + 1.e-20_dp)
                   ! remineralization of doc
                   docrem = MIN(dremdoc_dtb * ocetra(idoc)%ptr(i,j,k),        &
                          & 0.5_dp*(ocetra(ioxygen)%ptr(i,j,k)-5.e-8_dp)/ro2ut)
                ELSE
                   !changed from =max(remin,0.) etc. to =0. js3.5.2006
                   remin = 0._dp
                   docrem= 0._dp
                ENDIF
                ocetra(idet)%ptr(i,j,k) = ocetra(idet)%ptr(i,j,k) - remin     &
                     &                  + sterph + sterzo  
                ocetra(ialkali)%ptr(i,j,k) = ocetra(ialkali)%ptr(i,j,k)       &
                     &                     - rnit * (remin + docrem)
                ocetra(isco212)%ptr(i,j,k) = ocetra(isco212)%ptr(i,j,k)       &
                     &                     + rcar * (remin + docrem)
                ocetra(idoc)%ptr(i,j,k) = ocetra(idoc)%ptr(i,j,k) - docrem
                ocetra(ioxygen)%ptr(i,j,k) = ocetra(ioxygen)%ptr(i,j,k)       &
                     &                     - ro2ut * (remin + docrem)
                ocetra(iphosph)%ptr(i,j,k) = ocetra(iphosph)%ptr(i,j,k)       &
                     &                     + remin + docrem
                ocetra(iano3)%ptr(i,j,k) = ocetra(iano3)%ptr(i,j,k)           &
                     &                   + (remin + docrem) * rnit
                ocetra(iiron)%ptr(i,j,k) = ocetra(iiron)%ptr(i,j,k)           &
                     &                   + (remin + docrem) * riron           &
                     & - relaxfe_dtb*MAX(ocetra(iiron)%ptr(i,j,k)-fesoly,0._dp)
                !**************************************************************
                ! as ragueneau (2000) notes, Si(OH)4sat is about 1000 umol,
                ! but Si(OH)4 varies only between 0-100 umol so the
                ! expression dremopal*(Si(OH)4sat-Si(OH)4) would change the 
                ! rate only from 0 to 100%     
                !**************************************************************

                ! mz_bk_20120626+
                opalrem = dremopal_dtb * 0.1_dp * (tho(i,j,k) + 3.0_dp)       &
                     &  * ocetra(iopal)%ptr(i,j,k)
                ! opalrem = dremopal_dtb * ocetra(iopal)%ptr(i,j,k)
                ! mz_bk_20120626-
                ocetra(iopal)%ptr(i,j,k) = ocetra(iopal)%ptr(i,j,k)   - opalrem
                ocetra(isilica)%ptr(i,j,k)=ocetra(isilica)%ptr(i,j,k) + opalrem

                !**************************************************************
                ! There is about 1.e4 O2 on 1 N2O molecule (Broecker&Peng)
                ! refra : Tim Rixen, pers. communication
                !**************************************************************
                aou = satoxy(i,j,k) - ocetra(ioxygen)%ptr(i,j,k)
                refra = 1._dp + 3._dp * (0.5_dp + SIGN(0.5_dp,                &
                     &                                 aou - 1.97e-4_dp))
                ocetra(ian2o)%ptr(i,j,k) = ocetra(ian2o)%ptr(i,j,k)           &
                     &             + (remin + docrem) * 1.e-4_dp * ro2ut * refra
                ocetra(igasnit)%ptr(i,j,k) = ocetra(igasnit)%ptr(i,j,k)       &
                     &             - (remin + docrem) * 1.e-4_dp * ro2ut * refra
                ocetra(ioxygen)%ptr(i,j,k) = ocetra(ioxygen)%ptr(i,j,k)       &
                     &    - (remin + docrem) * 1.e-4_dp * ro2ut * refra * 0.5_dp

                ! careful, pho is very small at large depths
                ! (js: why careful? photolysis of DMS will be small,
                ! which it should)
                !  ocetra(idms)%ptr(i,j,k) = ocetra(idms)%ptr(i,j,k)          &
                !       &      - dmsp2 * 8. * pho * ocetra(idms)%ptr(i,j,k)   &
                !       &      - dmsp3 * abs(tho(i,j,k) + 3.)                 &
                !       &      * ocetra(idms)%ptr(i,j,k)

                ! mz_bk_20120708+
                ! testing factor dtb to get units 1/timestep instead of 1/d
                dms_bac = dmsp3 * abs(tho(i,j,k) + 3._dp)                   &
                     &  * ocetra(idms)%ptr(i,j,k)
                ! dms_bac = dtb * dmsp3 * abs(tho(i,j,k) + 3._dp)             &
                !      &  * ocetra(idms)%ptr(i,j,k)
                ! mz_bk_20120708-

                ocetra(idms)%ptr(i,j,k) = ocetra(idms)%ptr(i,j,k) - dms_bac

!!$!! #ifdef AGG
                IF (L_AGG) THEN
                   !***********************************************************
                   ! loss of snow aggregates (by numbers) due to
                   ! remineralization of poc
                   ! gain of snow aggregates (by numbers) due to zooplankton
                   ! mortality
                   ! NOTE that remin is in kmol/m3. Thus divide by avmass
                   ! (kmol/m3)
                   !***********************************************************
                   if (avmass .gt. 0._dp) then  
                      avnos = ocetra(inos)%ptr(i,j,k)
                      ocetra(inos)%ptr(i,j,k) = ocetra(inos)%ptr(i,j,k)       &
                           &                  - remin * avnos / avmass
                   endif
                   !***********************************************************
                   ! dead zooplankton corpses come with their own, flat
                   ! distribution
                   ! this flow even takes place if there is neither nos nor
                   ! mass
                   ! NOTE: zoomor is in kmol/m3!! Thus multiply flow by 1.e+6
                   !***********************************************************
                   zmornos = sterzo * zdis * 1.e+6_dp
                   ocetra(inos)%ptr(i,j,k) = ocetra(inos)%ptr(i,j,k) + zmornos
                ENDIF
!!$!! #endif /*AGG*/

             ENDIF !ddpo
          ENDDO
       ENDDO
    ENDDO

    !-----below euphotic zone

    !--------------------------------------------------------------------------
    DO k=kwrbioz+1,ke
       DO j=1,je
          DO i=1,ie
             !mz_ap_20080111+
             IF (ddpo(i,j,k) .GT. 0.5_dp) THEN
             !mz_ap_20080111-
                IF (ocetra(ioxygen)%ptr(i,j,k) .LT. 5.e-7_dp) THEN
                   ! denitrification
!!$!! #ifdef AGG
                   IF (L_AGG) THEN
                      avmass = ocetra(iphy)%ptr(i,j,k)+ocetra(idet)%ptr(i,j,k)
                   ENDIF
!!$!! #endif /*AGG*/

                   ! mz_bk_20120626+
                   remin = denitrification * drempoc_dtb                      &
                        & * MIN(ocetra(idet)%ptr(i,j,k),                      &
                        ! remineralization using NO3;
                        ! proportional to remineralization of poc 
                        & 0.5_dp * ocetra(iano3)%ptr(i,j,k) / nitdem)

                   !  remin = 0.5*drempoc_dtb*MIN(ocetra(idet)%ptr(i,j,k),    &
                   !       !                       ! remineralization using NO3
                   !       &  0.5*ocetra(iano3)%ptr(i,j,k)/rnit23)
                   ! mz_bk_20120626-

                   detref = remin / (ocetra(idet)%ptr(i,j,k) + 1.e-60_dp)
                   !                                                   ! P-units

                   !mz_ap_20071109+
                   !  detref=remin/(ocetra(idet)%ptr(i,j,k)+1.e-30)
                   !                                                   ! P-units
                   !mz_ap_20071109-

                   remin2o = dremn2o_dtb * MIN(ocetra(idet)%ptr(i,j,k),       &
                        ! remineralization using N2O
                        & 0.003_dp * ocetra(ian2o)%ptr(i,j,k) / (2._dp * ro2ut))

                   detrl = remin2o / (ocetra(idet)%ptr(i,j,k) + 1.e-60_dp)
                   !                                                    ! detrl?

                   !mz_ap_20071109+
                   !  detrl=remin2o/(ocetra(idet)%ptr(i,j,k)+1.e-30)
                   !mz_ap_20071109-

                   ! mz_bk_20120626+
                   !  ocetra(ialkali)%ptr(i,j,k)=ocetra(ialkali)%ptr(i,j,k)   &
                   !       &                    -rnit*(remin + remin2o)
                   ! mz_bk_20120626-

                   ocetra(isco212)%ptr(i,j,k) = ocetra(isco212)%ptr(i,j,k)    &
                        &                     + rcar * (remin + remin2o)

                   ocetra(idet)%ptr(i,j,k)    = ocetra(idet)%ptr(i,j,k)       &
                        &                     -        (remin + remin2o)

                   ocetra(iphosph)%ptr(i,j,k) = ocetra(iphosph)%ptr(i,j,k)    &
                        &                     +        (remin + remin2o)

                   ! mz_bk_20120626+
                   ocetra(iano3)%ptr(i,j,k) = ocetra(iano3)%ptr(i,j,k)        &
                        &                   - nitdem * remin + rnit * remin2o
                   ! ocetra(iano3)%ptr(i,j,k) = ocetra(iano3)%ptr(i,j,k)      &
                   !      &                   - rnit23 * remin                &
                   !      &                   + rnit * (remin + remin2o)
                   ! mz_bk_20120626-

                   ! mz_bk_20120626+
                   ocetra(igasnit)%ptr(i,j,k) = ocetra(igasnit)%ptr(i,j,k)    &
                        &           + n2prod * remin + 2._dp * ro2ut * remin2o
                   ! ocetra(igasnit)%ptr(i,j,k)=ocetra(igasnit)%ptr(i,j,k)    &
                   !      &         + rnit13 * remin + 2 * ro2ut * remin2o
                   ! mz_bk_20120626-

                   ocetra(iiron)%ptr(i,j,k) = ocetra(iiron)%ptr(i,j,k)        &
                        &                   + riron * (remin + remin2o)
                   ocetra(ian2o)%ptr(i,j,k) = ocetra(ian2o)%ptr(i,j,k)        &
                        &                   - 2._dp * ro2ut * remin2o

                   ! mz_bk_20120626+
                   !alkalinity is increased during denitrification due to
                   !consumption of H+ (see Wolf-Gladrow etal,2007)
                   ocetra(ialkali)%ptr(i,j,k) = ocetra(ialkali)%ptr(i,j,k)    &
                        &                     + nitdem * remin - rnit * remin2o
                   ! mz_bk_20120626-

                   ! mz_bk_20120626+
                   ! n2budget(i, j, k) = n2budget(i, j, k)                   &
                   !      &            + 2._wp * n2prod * remin
                   ! ! denitrification produces water (H2O), the corresponding
                   ! ! O2 uptake is budgeted in h2obudget
                   ! h2obudget(i, j, k) = h2obudget(i, j, k)                 &
                   !      &             + 0.5_wp * n2prod * remin

                   ! bgc_o_pro(i,j,k,kdenit) = nitdem * remin / dtbgc
                   ! mz_bk_20120626-


!!$!! #ifdef AGG
                   IF (L_AGG) THEN
                      !********************************************************
                      ! loss of snow aggregates (numbers) due to
                      ! remineralization of poc
                      ! NOTE that remin is in kmol/m3. Thus divide by avmass
                      ! (kmol/m3)
                      !********************************************************
                      if (avmass .gt. 0._dp) then  
                         avnos = ocetra(inos)%ptr(i,j,k)
                         ocetra(inos)%ptr(i,j,k) = ocetra(inos)%ptr(i,j,k)    &
                              &            - (remin + remin2o) * avnos / avmass
                      endif
                   ENDIF
!!$!! #endif /*AGG*/

                ENDIF
                ! mz_bk_20120626+
                ! sulphate reduction   ! introduced 11.5.2007 to improve
                ! poc-remineralisation in the
                ! oxygen minimum zone in the subsurface equatorial Pacific
                ! assumption of endless pool of SO4 (typical concentration
                ! are on the order of mmol/l)

                !      DO 301 k=kwrbioz+1,kpke
                !         DO 301 j=1,kpje
                !         DO 301 i=1,kpie
                IF (ocetra(ioxygen)%ptr(i,j,k) .LT. 3.e-6_dp) THEN
!!$!! #ifdef AGG
                   IF (L_AGG) THEN
                      avmass = ocetra(iphy)%ptr(i,j,k) + ocetra(idet)%ptr(i,j,k)
                   ENDIF
!!$!! #endif /*AGG*/

                   remin = sulfate_reduction * dtb * ocetra(idet)%ptr(i,j,k)
                   detref = sulfate_reduction * dtb

                   ocetra(idet)%ptr(i,j,k)    = ocetra(idet)%ptr(i,j,k)       &
                        &                     - remin
                   ocetra(ialkali)%ptr(i,j,k) = ocetra(ialkali)%ptr(i,j,k)    &
                        &                     - rnit * remin
                   ocetra(isco212)%ptr(i,j,k) = ocetra(isco212)%ptr(i,j,k)    &
                        &                     + rcar * remin
                   ocetra(iphosph)%ptr(i,j,k) = ocetra(iphosph)%ptr(i,j,k)    &
                        &                     + remin
                   ocetra(iano3)%ptr(i,j,k)   = ocetra(iano3)%ptr(i,j,k)      &
                        &                     + rnit  * remin
                   ocetra(iiron)%ptr(i,j,k)   = ocetra(iiron)%ptr(i,j,k)      &
                        &                     + riron * remin

                   ! ! sulphate reduction indirectly effects O2 bugdet, which
                   ! ! is budgeted in h2obudget
                   ! h2obudget(i,j,k) = h2obudget(i,j,k) - ro2ut * remin
                ENDIF
                !301   CONTINUE
                ! end sulphate reduction
                ! mz_bk_20120626-

                !mz_ap_20080111+
             ENDIF
             !mz_ap_20080111-
          ENDDO
       ENDDO
    ENDDO

!!$!! #ifdef AGG
    IF (L_AGG) THEN
       DO  k=1,ke
          DO j=1,je
             DO i=1,ie
                IF (ddpo(i,j,k) .GT. 0.5_dp) THEN
                   avmass = ocetra(iphy)%ptr(i,j,k) + ocetra(idet)%ptr(i,j,k)
                   snow = avmass * 1.e+6_dp
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!!$!! #endif /*AGG*/

    !js why twice ifdef AGG?
!!$!! #ifdef AGG
    IF (L_AGG) THEN
 
       ! ************************AGGREGATION*****(by Iris Kriest)**************
       ! General:
       ! Sinking speed, size distribution and aggregation are calculated 
       ! as in Kriest and Evans, 2000.
       ! I assume that opal and calcium carbonate sink at the same speed as P
       ! (mass).
       !
       ! Sinking speed and aggregation: I assume that if there is no
       ! phosphorous mass, the sinking speed is the maximal sinking speed of
       ! aggregates. I further assume that then there are no particles, and
       ! that the rate of aggregation is 0. This scheme removes no P in the
       ! absence of P, but still opal and/or calcium carbonate.
       ! This could or should be changed, because silica as well as carbonate
       ! shells will add to the aggregate mass, and should be considered.
       ! Puh. Does anyone know functional relationships between
       ! size and Si or CaCO3? Perhaps in a later version, I have to
       ! take the relationship between mass and size (i.e., density)?
       !
       ! 1. Size distribution and loss of marine snow aggregates due to
       ! aggregation (aggregate(i,j,k)) and sinking speed of mass and numbers
       ! (wmass(i,j,k) and wnumb(i,j,k) are calculated in a loop over 2-ke.
       !
       ! 2. The depth of the first layer may change due to ice drift, etc.
       ! This puts a restriction onto the maximum sinking speed.
       ! I currently set the max. size for size dependent sinking onto
       ! one appropriate for this depth.
       !
       ! 3. The fluxes out of the bottom layer are calculated from sinking
       ! speed and mass concentration, and form the boundary condition for
       ! the sediment.
       !
       ! 4. The fluxes in layer ke->2 are calculated from sinking speed and
       ! mass concentration, sinking speed, aggregation and number
       ! concentration.
       !
       ! 5. The fluxes in layer 1 are calculated from sinking speed and mass
       ! concentration, sinking speed, aggregation and number concentration.
       ! (??)
       !***********************************************************************
 
       do k=2,ke
 
          do i=1,ie
             do j=1,je
                if (ddpo(i,j,k) .gt. 0.5_dp) then
                   avm = ocetra(iphy)%ptr(i,j,k) + ocetra(idet)%ptr(i,j,k)
                   if (avm .gt. 0._dp) then
                      snow = avm * 1.e+6_dp
                      avn  = ocetra(inos)%ptr(i,j,k)
                      eps    = ((1._dp + FractDim) * snow - avn * cellmass)   &
                           & / (snow - avn * cellmass)
 
                      ! prevent epsilon from becoming exactly one of the
                      ! values which are needed for the division 
                      if (abs(eps - 3._dp) .lt. 1.e-15_dp)                    &
                           &   eps = 3._dp + vsmall
                      if (abs(eps - 4._dp) .lt. 1.e-15_dp)                    &
                           &   eps = 4._dp + vsmall
                      if (abs(eps - 3._dp - SinkExp) .lt. 1.e-15_dp)          &
                           &   eps = 3._dp + SinkExp + vsmall
                      if (abs(eps-1._dp-SinkExp-FractDim) .lt. 1.e-15_dp)     &
                           &   eps = 1._dp + SinkExp + FractDim + vsmall

                      e1 = 1._dp - eps
                      e2 = 2._dp - eps
                      e3 = 3._dp - eps
                      e4 = 4._dp - eps
                      es1 = e1 + SinkExp
                      es3 = e3 + SinkExp
                      TopF = (alar1 / alow1)**e1     ! alar1 'largest diameter',
                      !                              ! alow1 'smallest diameter'
                      TopM = TopF * TMFac

                      ! SINKING SPEED FOR THIS LAYER
                      wmass(i,j,k) = cellsink*((FractDim+e1)/(FractDim+es1)   &
                           &       + TopM * TSFac * SinkExp / (FractDim+es1))
                      wnumb(i,j,k) = cellsink*(e1/es1+TopF*TSFac*SinkExp/es1)

                      ! AGGREGATION

                      ! As a first step, assume that shear in the upper 4
                      ! layers is high and zero below. Subject to change.
                      ! js: this is for 20 layer version.
                      ! include 40 layer version
                      ! should be replaced by check for depth.
                      ! done 29072005js  (DONE!!)
                      if (k .le. n90depth) then
                         fshear = fsh
                      else
                         fshear = 0._dp
                      endif

                      ! shear kernel:
                      sagg1  = (TopF-1._dp) * (TopF*alar3-alow3) * e1 / e4    &
                           & + 3._dp*(TopF*alar1-alow1)*(TopF*alar2-alow2)    &
                           & * e1 * e1 / (e2 * e3)
                      sagg2  = TopF                                           &
                           & * ((alar3 + 3._dp                                &
                           & * (alar2*alow1*e1/e2 + alar1*alow2*e1/e3)        &
                           & + alow3 * e1 / e4)                               &
                           & - TopF*alar3*(1._dp + 3._dp*(e1/e2 + e1/e3)      &
                           & + e1 / e4))
                      sagg4 = TopF * TopF * 4._dp * alar3
                      shear_agg = (sagg1 + sagg2 + sagg4) * fshear

                      ! settlement kernel:
                      sagg1  = (TopF * TopF * alar2 * TSFac - alow2)          &
                           & * SinkExp / (es3 * e3 * (es3 + e1))              &
                           & + alow2 * ((1._dp - TopF * TSFac) / (e3 * es1)   &
                           & - (1._dp - TopF) / (es3*e1))
                      sagg2  = TopF*e1*(TSFac * ( alow2 - TopF * alar2)/e3    &
                           & - (alow2 - TopF * alar2 * TSFac) / es3)
                      sett_agg = (e1 * e1 * sagg1 + sagg2) * fse

                      effsti = Stick * (ocetra(iopal)%ptr(i,j,k) * 1.e+6_dp   &
                           & / ropal) / ((ocetra(iopal)%ptr(i,j,k)            &
                           & * 1.e+6_dp / ropal) + snow)

                      aggregate(i,j,k) = (shear_agg + sett_agg) * effsti      &
                           &           * avn * avn

                      ! dust aggregation:
                      ! shear kernel:
                      ! dustd3: dust diameter**3, d2: **2
                      dfirst = dustd3 + 3._dp * dustd2 * alar1                &
                           & + 3._dp * dustd1 * alar2 + alar3
                      dshagg = e1 * fsh * (dfirst * TopF / e1 - (             &
                           &   (TopF - 1._dp) / e1 * dustd3                   &
                           & + 3._dp * (TopF*alar1 - alow1) / e2 * dustd2     &
                           & + 3._dp * (TopF*alar2 - alow2) / e3 * dustd1     &
                           & + (TopF * alar3 - alow3) / e4))

                      ! settlement kernel:
                      dsett  = fse * dustd2 * ((e1 + SinkExp*TopF*TSFac)      &
                           & / es1 - dustsink / cellsink)

                      dustagg(i,j,k) = effsti*avn*ocetra(ifdust)%ptr(i,j,k)   &
                           &         * (dshagg + dsett)

                   else    ! available mass le 0
                      wmass(i,j,k)            = TSFac * cellsink
                      wnumb(i,j,k)            = 0._dp
                      aggregate(i,j,k)        = 0._dp
                      dustagg(i,j,k)          = 0._dp
                      ocetra(inos)%ptr(i,j,k) = 0._dp
                   endif

                endif   ! wet cell

             enddo   ! je
          enddo   ! ie
       enddo   ! ke


       ! EVALUATE SINKING RATE AND AGGREGATION FOR SURFACE LAYER, WHICH MAY BE
       ! LESS DEEP THAN INITIALLY SET BECAUSE OF EVAPORATION, ICE ETC.

       DO j=1,je
          DO i=1,ie
             if (ddpo(i,j,1) .gt. 0.5_dp ) then

                !ik evaluate safe length scale for size dependent sinking and
                !ik aggregation, and the resulting sinking rate and
                !ik aggregation rate.
                !ik zo may reduce the first layer depth to values that are
                !ik small and may cause the sinking length to exceed the
                !ik layers depth.
                !ik to be safe, for this upper layer set the upper size such
                !ik that loss due to sinking is at max the whole inventory of
                !ik this box.
                !ik aggregation will be calculated accordingly.

                checksize = (ddpo(i,j,1)/cellsink)**(1._dp/SinkExp) * alow1
                if (alar1 .gt. checksize) then
                   nacheck = nacheck + 1            ! js: seems not to be used
                endif
                flar1 = MIN(alar1,checksize)        ! reduce diameter of largest
                !                                   ! particle
                flar2 = flar1 * flar1
                flar3 = flar2 * flar1
                fTSFac = (flar1 / alow1)**SinkExp
                fTMFac = (flar1 / alow1)**FractDim

                ! SIZE DITRIBUTION
                avm = ocetra(iphy)%ptr(i,j,1) + ocetra(idet)%ptr(i,j,1)   
                ! available mass
                ! (js: add dust here to account for ballast effect?)
                if (avm .gt. 0._dp) then
                   snow = avm * 1.e+6_dp
                   avn = ocetra(inos)%ptr(i,j,1)             ! available numbers
                   eps = ((1._dp + FractDim) * snow - avn * cellmass)         &
                        ! exponential coefficient of size distribution
                        &  / (snow - avn * cellmass)

                   if (abs(eps - 3._dp) .lt. 1.e-15_dp) eps = 3._dp + vsmall
                   if (abs(eps - 4._dp) .lt. 1.e-15_dp) eps = 4._dp + vsmall
                   if (abs(eps - 3._dp - SinkExp) .lt. 1.e-15_dp)             &
                           eps = 3._dp + SinkExp + vsmall
                   if (abs(eps-1._dp-SinkExp-FractDim) .lt. 1.e-15_dp)        &
                        &  eps = 1._dp + SinkExp + FractDim + vsmall

                   e1 = 1._dp - eps
                   e2 = 2._dp - eps
                   e3 = 3._dp - eps
                   e4 = 4._dp - eps
                   es1 = e1 + SinkExp
                   es3 = e3 + SinkExp

                   fTopF = (flar1 / alow1)**e1
                   fTopM = fTopF * fTMFac

                   ! SINKING SPEEDS
                   wmass(i,j,1) = cellsink * ( (FractDim+e1)/(FractDim+es1)   &
                        &       + fTopM * fTSFac * SinkExp / (FractDim + es1))
                   wnumb(i,j,1) = cellsink * (e1/es1 + fTopF*fTSFac*SinkExp/es1)
                   ! AGGREGATION
                   sagg1  = (fTopF - 1._dp)*(fTopF*flar3 - alow3)*e1/e4       &
                        & + 3._dp * (fTopF*flar1-alow1)                       &
                        & * (fTopF*flar2 - alow2)*e1*e1 / (e2*e3)
                   sagg2  = fTopF * (                                         &
                        &   (flar3 + 3._dp                                    &
                        & * (flar2*alow1*e1/e2+flar1*alow2*e1/e3)             &
                        & + alow3*e1/e4)                                      &
                        & - fTopF*flar3*(1._dp+3._dp*(e1/e2 + e1/e3)+ e1/e4))
                   sagg4 = fTopF * fTopF * 4._dp * flar3
                   shear_agg = (sagg1 + sagg2 + sagg4) * fsh

                   sagg1  = (fTopF * fTopF * flar2 * fTSFac - alow2)          &
                        & * SinkExp / (es3 * e3 * (es3 + e1))                 &
                        & + alow2 * ((1._dp - fTopF * fTSFac) / (e3 * es1)    &
                        & - (1._dp - fTopF) / (es3*e1))
                   sagg2  = fTopF * e1 * (fTSFac*( alow2 - fTopF*flar2)/e3    &
                        & - (alow2 - fTopF * flar2 * fTSFac) / es3)
                   sett_agg = (e1 * e1 * sagg1 + sagg2) * fse

                   effsti = Stick*(ocetra(iopal)%ptr(i,j,1)*1.e+6 / ropal)    &
                        & / ((ocetra(iopal)%ptr(i,j,1)*1.e+6 / ropal) + snow)

                   aggregate(i,j,1) = (shear_agg + sett_agg)*effsti*avn*avn

                   ! dust aggregation:
                   ! shear kernel:
                   dfirst = dustd3 + 3._dp * dustd2 * flar1                   &
                        & + 3._dp * dustd1 * flar2 + flar3
                   dshagg = e1 * fsh * (dfirst * fTopF / e1 - (               &
                        &  (fTopF - 1._dp) / e1 * dustd3                      &
                        & + 3._dp * (fTopF * flar1 - alow1) / e2 * dustd2     &
                        & + 3._dp * (fTopF * flar2 - alow2) / e3 * dustd1     &
                        & + (fTopF * flar3 - alow3) / e4))

                   ! settlement kernel:
                   dsett  = fse * dustd2 * ((e1 + SinkExp * fTopF * fTSFac)   &
                        & / es1 - dustsink / cellsink)

                   dustagg(i,j,1) = effsti*avn*ocetra(ifdust)%ptr(i,j,1)      &
                        &         * (dshagg + dsett)

                else                            ! available mass le 0.
                   wmass(i,j,1)            = fTSFac * cellsink
                   wnumb(i,j,1)            = 0._dp
                   aggregate(i,j,1)        = 0._dp
                   dustagg(i,j,1)          = 0._dp
                   ocetra(inos)%ptr(i,j,1) = 0._dp
                endif

             endif    ! wet cell

          enddo
       enddo

       ! EVALUATE SINKING RATE AND AGGREGATION FOR BOTTOM LAYER, WHICH MAY BE
       ! LESS THICK THAN THE MINIMUM LAYER THICKNESS

       DO j=1,je
          DO i=1,ie
             if (ddpo(i,j,1) .gt. 0.5_dp) then
                if (alar1max(i,j) .lt. alar1) then

                   !ik take safe length scale for size dependent sinking and
                   !ik aggregation, and the resulting sinking rate and
                   !ik aggregation rate.

                   flar1 = alar1max(i,j)
                   flar2 = flar1 * flar1
                   flar3 = flar2 * flar1
                   fTSFac = TSFmax(i,j)
                   fTMFac = TMFmax(i,j)

                   ! SIZE DITRIBUTION
                   avm    = ocetra(iphy)%ptr(i,j,kbo(i,j))                    &
                        & + ocetra(idet)%ptr(i,j,kbo(i,j))
                   if (avm .gt. 0._dp ) then
                      snow = avm * 1.e+6_dp                         ! why *1.e6?
                      avn = ocetra(inos)%ptr(i,j,kbo(i,j))
                      eps    = ((1._dp + FractDim) * snow - avn*cellmass)     &
                           & / (snow - avn*cellmass)

                      if (abs(eps - 3._dp) .lt. 1.e-15_dp)                    &
                           &  eps = 3._dp + vsmall
                      if (abs(eps - 4._dp) .lt. 1.e-15_dp)                    &
                           &  eps = 4._dp + vsmall
                      if (abs(eps - 3._dp - SinkExp) .lt. 1.e-15_dp)          &
                           &  eps = 3._dp + SinkExp + vsmall
                      if (abs(eps-1._dp-SinkExp-FractDim) .lt. 1.e-15_dp)     &
                           &  eps = 1._dp + SinkExp + FractDim + vsmall

                      e1 = 1._dp - eps
                      e2 = 2._dp - eps
                      e3 = 3._dp - eps
                      e4 = 4._dp - eps
                      es1 = e1 + SinkExp
                      es3 = e3 + SinkExp

                      fTopF = (flar1 / alow1)**e1
                      fTopM = fTopF * fTMFac

                      ! SINKING SPEEDS
                      wmass(i,j,kbo(i,j)) = cellsink                          &
                           &            * ( (FractDim+e1) / (FractDim+es1)    &
                           &            + fTopM*fTSFac*SinkExp/(FractDim+es1) )
                      wnumb(i,j,kbo(i,j)) = cellsink                          &
                           &              * (e1/es1+fTopF*fTSFac*SinkExp/es1)

                      ! AGGREGATION
                      sagg1  = (fTopF - 1._dp)*(fTopF*flar3 - alow3)*e1/e4    &
                           & + 3._dp * (fTopF*flar1 - alow1)                  &
                           & * (fTopF * flar2 - alow2) * e1 * e1 / (e2 * e3)
                      sagg2  = fTopF * (                                      &
                           &   (flar3 + 3._dp * (flar2*alow1*e1/e2            &
                           &    + flar1 * alow2 * e1/e3) + alow3 * e1/e4)     &
                           & - fTopF * flar3 * (1._dp + 3._dp                 &
                           & * (e1 / e2 + e1 / e3) + e1 / e4))
                      sagg4  = fTopF * fTopF * 4._dp * flar3
                      shear_agg = (sagg1 + sagg2 + sagg4) * fsh

                      sagg1  = (fTopF * fTopF * flar2 * fTSFac - alow2)       &
                           & * SinkExp / (es3 * e3 * (es3 + e1))              &
                           & + alow2 * ((1._dp - fTopF * fTSFac)/(e3 * es1)   &
                           & - (1._dp - fTopF) / (es3 * e1))
                      sagg2  = fTopF*e1*(fTSFac*(alow2 - fTopF * flar2)/e3    &
                           & - (alow2 - fTopF * flar2 * fTSFac) / es3)
                      sett_agg =  (e1 * e1 * sagg1 + sagg2) * fse

                      effsti = Stick * (ocetra(iopal)%ptr(i,j,kbo(i,j))       &
                           & * 1.e+6_dp / ropal)                              &
                           & / ((ocetra(iopal)%ptr(i,j,kbo(i,j)) * 1.e+6_dp   &
                           &    / ropal) + snow)

                      aggregate(i,j,kbo(i,j)) = (shear_agg + sett_agg)        &
                           &                  * effsti * avn * avn

                      ! dust aggregation:
                      ! shear kernel:
                      dfirst = dustd3 + 3._dp * dustd2 * flar1                &
                           & + 3._dp * dustd1 * flar2 + flar3
                      dshagg = e1 * fsh * (dfirst * fTopF / e1                &
                           & - ((fTopF - 1._dp) / e1 * dustd3                 &
                           & + 3._dp * (fTopF * flar1 - alow1) / e2*dustd2    &
                           & + 3._dp * (fTopF * flar2 - alow2) / e3*dustd1    &
                           & + (fTopF * flar3 - alow3) / e4))

                      ! settlement kernel:
                      dsett  = fse*dustd2*((e1 + SinkExp*fTopF*fTSFac)/es1    &
                           & - dustsink / cellsink)

                      dustagg(i,j,kbo(i,j)) = effsti * avn                    &
                           &             * ocetra(ifdust)%ptr(i,j,kbo(i,j))   &
                           &             * (dshagg + dsett)

                   else
                      wmass(i,j,kbo(i,j))            = fTSFac * cellsink
                      wnumb(i,j,kbo(i,j))            = 0._dp
                      aggregate(i,j,kbo(i,j))        = 0._dp
                      dustagg(i,j,kbo(i,j))          = 0._dp
                      ocetra(inos)%ptr(i,j,kbo(i,j)) = 0._dp
                   endif ! avm

                endif ! alar1max

             endif ! ddpo

          enddo
       enddo

       !IK COMPUTE FLUXES FOR BOUNDARY CONDITION/BOTTOM LAYER
       !js fluxes to sediment   (still AGG)

       ! mz_bk_20120627+
       DO j=1,je
          DO i=1,ie
       ! DO 36 j=1,je
       !    DO 36 i=1,ie
       ! mz_bk_20120627-
             if (ddpo(i,j,kbo(i,j)) .gt. 0.5_dp) then
                wphy = wmass(i,j,kbo(i,j))
                prorca(i,j) = ocetra(iphy)%ptr(i,j,kbo(i,j)) * wphy           &
                     &      + ocetra(idet)%ptr(i,j,kbo(i,j)) * wphy
                prcaca(i,j) = ocetra(icalc)%ptr(i,j,kbo(i,j)) * wphy
                silpro(i,j) = ocetra(iopal)%ptr(i,j,kbo(i,j)) * wphy
                produs(i,j) = ocetra(ifdust)%ptr(i,j,kbo(i,j)) * dustsink     &
                     &      + ocetra(iadust)%ptr(i,j,kbo(i,j)) * wphy
             endif
       ! mz_bk_20120627+
       ! 36    CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-

       ! mz_bk_20120627+
       DO k=ke,2,-1
          DO j=1,je
             DO i=1,ie
       ! DO 2 K=ke,2,-1
       !    DO 34 j=1,je
       !       DO 34 i=1,ie
       ! mz_bk_20120627-
                if (ddpo(i,j,k) .gt. 0.5_dp) then

                   ! SINKING SPEED FOR UPPER LAYER
                   wphyup = wmass(i,j,k-1)      ! settling velocity of mass
                   wnosup = wnumb(i,j,k-1)      ! settling velocity of number
                   !                            ! of marine snow aggregates

                   ! SINKING SPEED FOR ACTUAL LAYER
                   wphy = wmass(i,j,k)
                   wnos = wnumb(i,j,k)

                   ! SUM-UP FLUXES (compute new concentrations)
                   ocetra(iphy)%ptr(i,j,k) = ocetra(iphy)%ptr(i,j,k)          &
                        &               + (ocetra(iphy)%ptr(i,j,k-1) * wphyup &
                        &               - ocetra(iphy)%ptr(i,j,k) * wphy)     &
                        &               * dpio(i,j,k)
                   ocetra(idet)%ptr(i,j,k) = ocetra(idet)%ptr(i,j,k)          &
                        &               + (ocetra(idet)%ptr(i,j,k-1) * wphyup &
                        &               - ocetra(idet)%ptr(i,j,k) * wphy)     &
                        &               * dpio(i,j,k)
                   ocetra(icalc)%ptr(i,j,k) = ocetra(icalc)%ptr(i,j,k)        &
                        &               + (ocetra(icalc)%ptr(i,j,k-1)*wphyup  &
                        &               - ocetra(icalc)%ptr(i,j,k) * wphy)    &
                        &               * dpio(i,j,k)
                   ocetra(iopal)%ptr(i,j,k) = ocetra(iopal)%ptr(i,j,k)        &
                        &               + (ocetra(iopal)%ptr(i,j,k-1)*wphyup  &
                        &               - ocetra(iopal)%ptr(i,j,k) * wphy)    &
                        &               * dpio(i,j,k)
                   ocetra(inos)%ptr(i,j,k) = ocetra(inos)%ptr(i,j,k)          &
                        &               - aggregate(i,j,k)                    &
                        &               + (ocetra(inos)%ptr(i,j,k-1) * wnosup &
                        &               - ocetra(inos)%ptr(i,j,k) * wnos)     &
                        &               * dpio(i,j,k)
                   ! sinking of free dust and loss due to attachment to
                   ! aggregated dust
                   ocetra(ifdust)%ptr(i,j,k) = ocetra(ifdust)%ptr(i,j,k)      &
                        &               - dustagg(i,j,k)                      &
                        &               + (ocetra(ifdust)%ptr(i,j,k-1)        &
                        &               - ocetra(ifdust)%ptr(i,j,k))          &
                        &               * dustsink * dpio(i,j,k)
                   ! sinking of aggregated dust and gain due to attachment of
                   ! free dust to aggregates
                   ocetra(iadust)%ptr(i,j,k) = ocetra(iadust)%ptr(i,j,k)      &
                        &               + dustagg(i,j,k)                      &
                        &               + (ocetra(iadust)%ptr(i,j,k-1)*wphyup &
                        &               - ocetra(iadust)%ptr(i,j,k) * wphy)   &
                        &               * dpio(i,j,k)
                endif

       ! mz_bk_20120627+
       ! 34 CONTINUE   ! end i,j-loop
             ENDDO
          ENDDO
       !  2 CONTINUE   ! end k-loop
       ENDDO
       ! mz_bk_20120627-

       !IK  COMPUTE FLUXES FOR SURFACE LAYER
       ! mz_bk_20120627+
       DO j=1,je
          DO i=1,ie
       ! DO 35 j=1,je
       !    DO 35 i=1,ie
       ! mz_bk_20120627-
             if (ddpo(i,j,1) .gt. 0._dp) then
                wphy = wmass(i,j,1)
                wnos = wnumb(i,j,1)

                ! SUM-UP FLUXES
                ocetra(iphy)%ptr(i,j,1) = ocetra(iphy)%ptr(i,j,1)             &
                     &           - ocetra(iphy)%ptr(i,j,1) * wphy * dpio(i,j,1)
                ocetra(idet)%ptr(i,j,1) = ocetra(idet)%ptr(i,j,1)             &
                     &           - ocetra(idet)%ptr(i,j,1) * wphy * dpio(i,j,1)
                ocetra(icalc)%ptr(i,j,1) = ocetra(icalc)%ptr(i,j,1)           &
                     &           - ocetra(icalc)%ptr(i,j,1) * wphy * dpio(i,j,1)
                ocetra(iopal)%ptr(i,j,1) = ocetra(iopal)%ptr(i,j,1)           &
                     &           - ocetra(iopal)%ptr(i,j,1) * wphy * dpio(i,j,1)
                ocetra(inos)%ptr(i,j,1) = ocetra(inos)%ptr(i,j,1)             &
                     &           - aggregate(i,j,1)                           &
                     &           - ocetra(inos)%ptr(i,j,1) * wnos * dpio(i,j,1)

                ! sinking of free dust and loss of free dust to
                ! aggregated dust
                ocetra(ifdust)%ptr(i,j,1) = ocetra(ifdust)%ptr(i,j,1)         &
                     &      - dustagg(i,j,1)                                  &
                     &      - ocetra(ifdust)%ptr(i,j,1) * dustsink * dpio(i,j,1)

                ! sinking of aggregated dust and gain due to attachment of
                ! free dust to aggregates
                ocetra(iadust)%ptr(i,j,1) = ocetra(iadust)%ptr(i,j,1)         &
                     &          + dustagg(i,j,1)                              &
                     &          - ocetra(iadust)%ptr(i,j,1) * wphy * dpio(i,j,1)

             endif
       ! mz_bk_20120627+
       ! 35    CONTINUE 
          ENDDO
       ENDDO
       ! mz_bk_20120627-

       !-------------------------------------------end aggregation part
    ELSE !(L_AGG)

       ! implicit method:
       ! C(k,T+dt)=C(k,T) + (w*dt/ddpo(k))*(C(k-1,T+1)-C(k,T+1))
       ! -->    
       ! C(k,T+dt)=(ddpo(k)*C(k,T)+w*dt*C(k-1,T+dt))/(ddpo(k)+w*dt)
       ! sedimentation=w*dt*C(ks,T+dt)
       k=1                         ! -----------surface layer
       DO j=1,je
          DO i=1,ie
             IF (ddpo(i,j,k) .GT. 0.5_dp) THEN
                ocetra(idet)%ptr(i,j,k)  = (ocetra(idet)%ptr(i,j,k)           &
                     &                   * ddpo(i,j,k)) / (ddpo(i,j,k) + wpoc)
                ocetra(icalc)%ptr(i,j,k) = (ocetra(icalc)%ptr(i,j,k)          &
                     &                   * ddpo(i,j,k)) / (ddpo(i,j,k) + wcal)
                ocetra(iopal)%ptr(i,j,k) = (ocetra(iopal)%ptr(i,j,k)          &
                     &                   * ddpo(i,j,k)) / (ddpo(i,j,k)+wopal)
                ocetra(ifdust)%ptr(i,j,k)= (ocetra(ifdust)%ptr(i,j,k)         &
                     &                   * ddpo(i,j,k)) / (ddpo(i,j,k)+wdust)
             ENDIF
          enddo
       enddo

       ! mz_bk_20120627+
       DO k=2,ke                ! ------------ water column
          DO j=1,je
             DO i=1,ie
       ! DO 10 k=2,ke                ! ------------ water column
       !    DO 12 j=1,je
       !       DO 12 i=1,ie
       ! mz_bk_20120627-
                IF (ddpo(i,j,k) .GT. 0.5_dp) THEN
                   ocetra(idet)%ptr(i,j,k) = (ocetra(idet)%ptr(i,j,k)         &
                        &                 * ddpo(i,j,k)                       &
                        &                 + ocetra(idet)%ptr(i,j,k-1) * wpoc) &
                        &                 / (ddpo(i,j,k) + wpoc)
                   ocetra(icalc)%ptr(i,j,k) = (ocetra(icalc)%ptr(i,j,k)       &
                        &                * ddpo(i,j,k)                        &
                        &                + ocetra(icalc)%ptr(i,j,k-1) * wcal) &
                        &                / (ddpo(i,j,k) + wcal)
                   ocetra(iopal)%ptr(i,j,k) = (ocetra(iopal)%ptr(i,j,k)       &
                        &               * ddpo(i,j,k)                         &
                        &               + ocetra(iopal)%ptr(i,j,k-1) * wopal) &
                        &               / (ddpo(i,j,k) + wopal)        
                   ocetra(ifdust)%ptr(i,j,k) = (ocetra(ifdust)%ptr(i,j,k)     &
                        &              * ddpo(i,j,k)                          &
                        &              + ocetra(ifdust)%ptr(i,j,k-1) * wdust) &
                        &              / (ddpo(i,j,k) + wdust)        
                ENDIF
       ! mz_bk_20120627+
       ! 12      CONTINUE
             ENDDO
          ENDDO
       ! 10    CONTINUE
       ENDDO
       ! mz_bk_20120627-

       ! ------------------------------flux to sediment
       ! mz_bk_20120627+
       DO j=1,je
          DO i=1,ie
       ! DO 33 j=1,je
       !    DO 33 i=1,ie
       ! mz_bk_20120627-
             IF (ddpo(i,j,1) .GT. 0.5_dp) THEN
                prorca(i,j) = ocetra(idet)%ptr(i,j,kbo(i,j)) * wpoc
                prcaca(i,j) = ocetra(icalc)%ptr(i,j,kbo(i,j)) * wcal
                silpro(i,j) = ocetra(iopal)%ptr(i,j,kbo(i,j)) * wopal
                produs(i,j) = ocetra(ifdust)%ptr(i,j,kbo(i,j)) * wdust
             ENDIF
       ! mz_bk_20120627+
       ! 33    CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-

    ENDIF !(L_AGG)

  END SUBROUTINE hamocc_ocprod
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE hamocc_cyano(ie,je,ke,ddpo)
    !
    ! correspondes to cyano.f90
    ! Nitrate reduction by cyano bacteria (2NO3 + O2 => N2O + O2).
    !
    ! Original:
    !     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    !     Modified
    !     --------
    !     S.Legutke,             *MPI-MaD, HH*    10.04.01
    !     - included : surface reduction of gaseous nitrogen
    !------------------------------------------------------------------------
    implicit none

    INTEGER,  INTENT(IN) :: ie,je,ke 
    REAL(DP), INTENT(IN) :: ddpo(ie,je,ke)

    !LOCAL
    INTEGER ::i,j,k
    REAL(DP):: oldocetra

    !
    ! N-fixation by cyano bacteria from atmospheric nitrogen if nitrate is
    ! below redfield ratio wrt phsophate(this is not a surface flux!)
    ! N2 = nitrogen, NO3 = nitrate
    ! from tech report: 
    ! uptake of atmospheric nitrogen and its immediate release as nitrate
    ! by diazotrophs
    !
    DO j=1,je
       DO i=1,ie
          IF (ddpo(i,j,1) .GT. 0.5_dp) THEN
             IF (ocetra(iano3)%ptr(i,j,1) .LT.                                &
                  &     (rnit * ocetra(iphosph)%ptr(i,j,1))) THEN
                oldocetra = ocetra(iano3)%ptr(i,j,1)
                ! mz_bk_20120627+
                ocetra(iano3)%ptr(i,j,1) = ocetra(iano3)%ptr(i,j,1)           &
                     &                   * (1._dp - n2_fixation * dtb)        &
                     &                   + n2_fixation * dtb * rnit           &
                     &                   * ocetra(iphosph)%ptr(i,j,1)
                ! ocetra(iano3)%ptr(i,j,1) = ocetra(iano3)%ptr(i,j,1)         &
                !      &         * (1-bluefix)                                &
                !      &         + bluefix * rnit * ocetra(iphosph)%ptr(i,j,1)
                ! mz_bk_20120627-
                ocetra(igasnit)%ptr(i,j,1) = ocetra(igasnit)%ptr(i,j,1)       &
                     &        - (ocetra(iano3)%ptr(i,j,1) - oldocetra) * .5_dp
                !            *(1./2.) half of the nitrogen from gaseous pool (?)
                ocetra(ioxygen)%ptr(i,j,1) = ocetra(ioxygen)%ptr(i,j,1)       &
                     &        - (ocetra(iano3)%ptr(i,j,1) - oldocetra) * 1.5_dp
                !       1.5 (?) should be 172./122. redfield ratios O2/C , PO4/C
                ! alkalinity is decreased during n-fixation due to release of H+
                ! change in ALK is identical to change in NO3
                ocetra(ialkali)%ptr(i,j,1) = ocetra(ialkali)%ptr(i,j,1)       & 
                     &                 - (ocetra(iano3)%ptr(i,j,1) - oldocetra)
                ! mz_bk_20120627+
                ! n2budget(i, j, 1) = n2budget(i, j, 1)                       &
                !      &            - (ocetra(i,j,1,iano3) - oldocetra)
                ! mz_bk_20120627-
             ENDIF
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE hamocc_cyano
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE hamocc_carchm(ie,je,ke,sao,tho,ddpo)
    !
    ! Correspondes to carchm.f90
    ! Calculation of inorganic carbon cycle
    !
    ! Original:
    !     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    !     Modified
    !     --------
    !     S.Legutke,        *MPI-MaD, HH*    10.04.01
    !------------------------------------------------------------------------
    implicit none

    INTEGER, INTENT(IN) :: ie,je,ke

    REAL(DP), INTENT(IN) :: sao(ie,je,ke)
    REAL(DP), INTENT(IN) :: ddpo(ie,je,ke)
    REAL(DP), INTENT(IN) :: tho(ie,je,ke)

    INTEGER :: i,j,k,iter
    INTEGER :: kk,ii,jj

    REAL(DP):: supsat, undsa, dissol
    REAL(DP):: dddhhh,dadh,a,h,c,alk,t1,t2
    REAL(DP):: akbi,ak2i,ak1i
    REAL(DP):: AK0,AK1,AK2,AKB,AKW,BT,oxysa,anisa

    ! mz_bk_20120627+
    ! REAL(DP):: thickness
    ! mz_bk_20120627-

    k=1      ! surface layer

    ! mz_bk_20120627+
    DO j=1,je
       DO i=1,ie
    ! DO 1 j=1,je
    ! DO 1 i=1,ie
    ! mz_bk_20120627-

          IF (ddpo(i,j,1) .GT. 0.5_dp) THEN

             !*       21.11 SET CHEMICAL CONSTANTS

             AK0   = CHEMCM(i,j,5)
             AK1   = CHEMCM(i,j,4)
             AK2   = CHEMCM(i,j,3)
             AKB   = CHEMCM(i,j,1)
             AKW   = CHEMCM(i,j,2)
             BT    = CHEMCM(i,j,6)
             oxysa = CHEMCM(i,j,7)
             anisa = CHEMCM(i,j,8)


             ak1i = 1._dp / ak1
             ak2i = 1._dp / ak2
             akbi = 1._dp / akb

             !***************************************************************
             !
             ! Compute the Schmidt number of CO2 in seawater and the transfer
             ! (piston) velocity using the formulation presented
             ! by Wanninkhof (1992, J. Geophys. Res., 97, 7373-7382).
             ! Input is temperature in deg C.
             ! CO2 Schmidt number after
             ! Wanninkhof (1992, J. Geophys. Res., 97, 7373-7382)
             ! DMS Schmidt number after
             ! Saltzmann et al. (1993, J. Geophys. Res. 98, 16,481-16,486)
             ! O2 Schmidt number after
             ! Keeling et al. (1998, Global Biogeochem. Cycles, 12, 141-163)
             ! CFC Schmidt number ref:  Zheng et al (1998), JGR, vol 103,No C1
             !  ---> www.ipsl.jussieu.fr/OCMIP/
             !***************************************************************
!!$!!!
!!$!!       scco2 = 2073.1 - 125.62*tho(i,j,1) + 3.6276*tho(i,j,1)**2  &
!!$!!     &       - 0.043219*tho(i,j,1)**3
!!$!!
!!$!!       scdms = 2674.0-147.12*tho(i,j,1)+3.726*tho(i,j,1)**2       &
!!$!!     &       - 0.038*tho(i,j,1)**3
!!$!!
!!$!!       sco2  = 1638.0 - 81.83*tho(i,j,1) + 1.483*tho(i,j,1)**2    & 
!!$!!     &       - 0.008004*tho(i,j,1)**3  

             !TODO: WE NEED FLUXES OF: CO2, O2, N2,N2O

             ! Calculate new hi concentration  (four iterations)

             h = hi(i,j,k)
             c = ocetra(isco212)%ptr(i,j,k)
             akw = akw3(i,j,k)                         ! IONIC PRODUCT OF WATER
             bt  = rrrcl * sao(i,j,k)                  ! salinity
             akb = akb3(i,j,k)
             alk = ocetra(ialkali)%ptr(i,j,k)
             t1 = h * ak1i
             t2 = h * ak2i
             a      = c * (2._dp + t2) / (1._dp + t2 + t2*t1) + akw/h - h     &
                  & + bt / (1._dp + h * akbi) - alk
             dadh   = c * (1._dp / (ak2 * (1._dp + t2 + t2 * t1))             &
                  & - (2._dp + t2) * ((1._dp + 2._dp * t1) * ak2i)            &
                  & / (1._dp + t2 + t2 * t1)**2)                              &
                  & - akw / h**2 - 1._dp - (bt * akbi) / (1._dp + h * akbi)**2

             dddhhh = a / dadh
             ! mz_bk_20110208+
             ! hi(i,j,k)=hi(i,j,k)-dddhhh
             ! h=hi(i,j,k)
             h = MAX(h - dddhhh, 1.e-10_dp)   ! Prevent overshooting to negative
             !                                ! values at start of iteration
             ! mz_bk_20110208-

             t1 = h * ak1i
             t2 = h * ak2i
             a      = c * (2._dp + t2) / (1._dp + t2 + t2*t1) + akw/h - h     &
                  & + bt / (1._dp + h * akbi) - alk
             dadh   = c * (1._dp / (ak2 * (1._dp + t2 + t2 * t1))             &
                  & - (2._dp + t2) * ((1._dp + 2._dp * t1) * ak2i)            &
                  & / (1._dp + t2 + t2 * t1)**2)                              &
                  & - akw / h**2 - 1._dp - (bt * akbi) / (1._dp + h * akbi)**2
             dddhhh = a / dadh
             ! mz_bk_20110208+
             ! hi(i,j,k)=hi(i,j,k)-dddhhh
             ! h=hi(i,j,k)
             h = MAX(h - dddhhh, 1.e-10_dp)   ! Prevent overshooting to negative
             !                                ! values at start of iteration
             ! mz_bk_20110208-

             t1 = h * ak1i
             t2 = h * ak2i
             a      = c * (2._dp + t2) / (1._dp + t2 + t2*t1) + akw/h - h     &
                  & + bt / (1._dp + h * akbi) - alk
             dadh   = c * (1._dp / (ak2 * (1._dp + t2 + t2 * t1))             &
                  & - (2._dp + t2) * ((1._dp + 2._dp * t1) * ak2i)            &
                  & / (1._dp + t2 + t2 * t1)**2)                              &
                  & - akw / h**2 - 1._dp - (bt * akbi) / (1._dp + h * akbi)**2
             dddhhh = a / dadh
             ! mz_bk_20110208+
             ! hi(i,j,k)=hi(i,j,k)-dddhhh
             ! h=hi(i,j,k)
             h = MAX(h - dddhhh, 1.e-10_dp)   ! Prevent overshooting to negative
             !                                ! values at start of iteration
             hi(i,j,k) = h
             ! mz_bk_20110208-

             !mz_ap_20080206+
             co2(i,j) = (c / ((1._dp + ak1 * (1._dp + ak2/h) / h) * ak0)) * ak0
             !mz_ap_20080206-

          ENDIF ! wet cell

    ! mz_bk_20120627+
    ! 1 CONTINUE    ! i,j loop
       ENDDO
    ENDDO
    ! mz_bk_20120627-

    !     -----------------------------------------------------------------
    !*        22. CHEMICAL CONSTANTS - water column
    ! mz_bk_20120729+
    ! ! mz_bk_20120627+
    ! ! three iterations
    ! DO iter=1,3
    ! ! mz_bk_20120627-
    !    DO k=1,ke
    !       DO j=1,je
    !          DO i=1,ie
    !             IF (ddpo(i,j,k) .GT. 0.5_dp) THEN
    !                h   = hi(i,j,k)
    !                c   = ocetra(isco212)%ptr(i,j,k)
    !                t1  = h / ak13(i,j,k)
    !                t2  = h / ak23(i,j,k)
    !                ak2 = ak23(i,j,k)
    !                akw = akw3(i,j,k)
    !                bt  = rrrcl * sao(i,j,k)
    !                akb = akb3(i,j,k)
    !                alk = ocetra(ialkali)%ptr(i,j,k)
    !                ! Determine hydrogen ion HI so that ALK(DIC,BT,HI) matches
    !                ! given alk by Newton iteration
    !                ! Actual mismatch
    !                a      = c*(2._dp + t2)/(1._dp + t2 + t2*t1) + akw/h - h   &
    !                     & + bt / (1._dp + h/akb) - alk
    !                ! Derivative
    !                dadh   = c * (1._dp / (ak2 * (1._dp + t2 + t2 * t1))       &
    !                     & - (2._dp + t2) * (1._dp / ak2 + 2._dp * t1 / ak2)   &
    !                     & / (1._dp + t2 + t2 * t1)**2)                        &
    !                     & - akw / h**2 - 1._dp - (bt/akb) / (1._dp + h/akb)**2
    !                dddhhh = a / dadh
    !                ! mz_bk_20120627+
    !                h = MAX(h - dddhhh, 1.e-10_dp)   ! Prevent overshooting to
    !                !                     ! negative values at start of iteration
    !                hi(i,j,k) = h
    !                ! h=h-dddhhh
    !                ! hi(i,j,k)=hi(i,j,k)-dddhhh
    !                ! mz_bk_20120627-
    !                co3(i,j,k) = c / (1._dp + h * (1._dp + h / ak13(i,j,k))    &
    !                     &     / ak23(i,j,k))
    !             ENDIF ! wet cell
    !          ENDDO ! i,j,k loop
    !       ENDDO ! i,j,k loop
    !    ENDDO ! i,j,k loop
    ! ! mz_bk_20120627+
    ! ENDDO  ! iteration
    ! ! mz_bk_20120627-
    DO k=1,ke
       DO j=1,je
          DO i=1,ie
             IF (ddpo(i,j,k) .GT. 0.5_dp) THEN
                c   = ocetra(isco212)%ptr(i,j,k)
                alk = ocetra(ialkali)%ptr(i,j,k)
                ak2 = ak23(i,j,k)
                akw = akw3(i,j,k)
                akb = akb3(i,j,k)
                bt  = rrrcl * sao(i,j,k)
                h   = hi(i,j,k)
                ! mz_bk_20120627+
                ! three iterations
                DO iter=1,3
                   ! mz_bk_20120627-
                   t1  = h / ak13(i,j,k)
                   t2  = h / ak23(i,j,k)
                   ! Determine hydrogen ion HI so that ALK(DIC,BT,HI) matches
                   ! given alk by Newton iteration
                   ! Actual mismatch
                   a      = c*(2._dp + t2)/(1._dp + t2 + t2*t1) + akw/h - h   &
                        & + bt / (1._dp + h/akb) - alk
                   ! Derivative
                   dadh   = c * (1._dp / (ak2 * (1._dp + t2 + t2 * t1))       &
                        & - (2._dp + t2) * (1._dp / ak2 + 2._dp * t1 / ak2)   &
                        & / (1._dp + t2 + t2 * t1)**2)                        &
                        & - akw / h**2 - 1._dp - (bt/akb) / (1._dp + h/akb)**2
                   dddhhh = a / dadh
                   ! mz_bk_20120627+
                   h = MAX(h - dddhhh, 1.e-10_dp)   ! Prevent overshooting to
                   !                     ! negative values at start of iteration
                   hi(i,j,k) = h
                   ! h=h-dddhhh
                   ! hi(i,j,k)=hi(i,j,k)-dddhhh
                   ! mz_bk_20120627-
                   ! mz_bk_20120627+
                ENDDO  ! iteration
                ! mz_bk_20120627-
                co3(i,j,k) = c / (1._dp + h * (1._dp + h / ak13(i,j,k))    &
                     &     / ak23(i,j,k))
             ENDIF ! wet cell
          ENDDO ! i,j,k loop
       ENDDO ! i,j,k loop
    ENDDO ! i,j,k loop
    ! mz_bk_20120729-


    !
    ! Dissolution of calcium
    ! Note : mixed layer (k=1) is assumed to be always supersaturated
    !        (saturation in reality depends on temperature/DIC/alkalinity)
    !
    DO k=2,ke
       DO j=1,je
          DO i=1,ie
             IF (ddpo(i,j,k) .GT. 0.5_dp) THEN
                supsat = co3(i,j,k) - 97._dp * aksp(i,j,k)  
                ! 97. = 1./1.03e-2 (MEAN TOTAL [CA++] IN SEAWATER [kmol/m3])
                undsa  = MAX(0._dp, -supsat)
                dissol = MIN(undsa, dremcalc_dtb * ocetra(icalc)%ptr(i,j,k))
                ocetra(icalc)%ptr(i,j,k) = ocetra(icalc)%ptr(i,j,k) - dissol
                ocetra(ialkali)%ptr(i,j,k) = ocetra(ialkali)%ptr(i,j,k)       &
                     &                     + 2._dp * dissol

                ocetra(isco212)%ptr(i,j,k) = ocetra(isco212)%ptr(i,j,k) + dissol
             ENDIF   ! wet cell
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE hamocc_carchm
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE hamocc_powach(ie,je,ke,sao)
    !
    ! Correspondes to powach.f90
    !
    ! Original:
    !     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    !     Modified
    !     --------
    !     S.Legutke,        *MPI-MaD, HH*    10.04.01
    !
    !     S.Lorenz/JO.Beismann, OpenMP parallel    *MPI-Met, HH*  24.08.07
    !------------------------------------------------------------------------
    implicit none

    INTEGER :: i,j,k, iter
    INTEGER, INTENT(IN)  :: ie,je,ke

    REAL(DP), INTENT(IN) :: sao(ie,je,ke)

    REAL(DP):: sedb1(ie,0:ks),sediso(ie,0:ks)
    REAL(DP):: solrat(ie,ks),powcar(ie,ks)
    REAL(DP):: aerob(ie,ks),anaerob(ie,ks),ansulf(ie,ks)

    REAL(DP):: disso, dissot, undsa, silsat, posol, dissot1, dissot2
    REAL(DP):: umfa,denit,hconve,bt,alk,c
    REAL(DP):: ak1,ak2,akb,akw
    REAL(DP):: ratc13,ratc14,rato13,rato14,poso13,poso14
    REAL(DP):: h,t1,t2,a,dadh,dddhhh,reduk,satlev

    ! *****************************************************************
    ! accelerated sediment
    ! needed for boundary layer ventilation in fast sediment routine      
    !           REAL(DP):: wo(ie,je,ke+1)
    REAL(DP):: bolven(ie)

    ! A LOOP OVER J
    ! RJ: This loop must go from 1 to je in the parallel version,
    !     otherways we had to do a boundary exchange
    !e    write(0,737)nsedtra
    !737   format('powach',i6)

    ! Silicate saturation concentration is 1 mol/m3
    silsat = 0.001_dp

    ! Dissolution rate constant of opal (disso) [1/(kmol Si(OH)4/m3)*1/sec]
    disso  = 3.e-8_dp  ! (2011-01-04) EMR
    !  disso = 1.e-6_dp ! test vom 03.03.04 half live sil ca. 20.000 yr
    dissot = disso * dtbgc

    ! Degradation rate constant of POP (disso) [1/(kmol O2/m3)*1/sec]
    disso   = 0.01_dp / 86400._dp  !  disso=3.e-5 was quite high
    dissot1 = disso * dtbgc

    ! Dissolution rate constant of CaCO3 (disso) [1/(kmol CO3--/m3)*1/sec]
    disso   = 1.e-7_dp
    dissot2 = disso * dtbgc

    !ik      denit = 1.e-6*dtbgc
    denit = 0.01_dp / 86400._dp *dtbgc

    ! mz_bk_20120627+
    DO j=1,je
    ! DO 8888 j=1,je
    ! mz_bk_20120627-
     
       ! mz_bk_20120627+
       DO k=1,ks
          DO i=1,ie
       ! DO 1189 k=1,ks
       ! DO 1189 i=1,ie
       ! mz_bk_20120627-

             solrat(i,k) = 0._dp
             powcar(i,k) = 0._dp
             anaerob(i,k)= 0._dp
             aerob(i,k)  = 0._dp
             ansulf(i,k) = 0._dp

       ! mz_bk_20120627+
       ! 1189  CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-

       ! calculate bottom ventilation rate for scaling of sediment-water
       ! exchange
       ! mz_bk_20120627+
       do i=1,ie      
       ! do 1170 i=1,ie      
       ! mz_bk_20120627-
          bolven(i) = 1._dp
       ! mz_bk_20120627+
       ! 1170  continue
       enddo
       ! mz_bk_20120627-

       ! mz_bk_20120627+
       DO k=0,ks
          DO i=1,ie
       ! DO 1171 k=0,ks
       ! DO 1171 i=1,ie
       ! mz_bk_20120627-
             sedb1(i,k) = 0._dp
             sediso(i,k)= 0._dp
       ! mz_bk_20120627+
       ! 1171  CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-


       ! CALCULATE SILICATE-OPAL CYCLE AND SIMULTANEOUS SILICATE DIFFUSION
       !******************************************************************

       ! mz_bk_20120627+
       ! ! Dissolution rate constant of opal (disso) [1/(kmol Si(OH)4/m3)*1/sec]
       ! !      disso=1.e-8 
       !       disso=1.e-6 ! test vom 03.03.04 half live sil ca. 20.000 yr 
       !       dissot=disso*dtbgc

       ! ! Silicate saturation concentration is 1 mol/m3

       !       silsat=0.001
       ! mz_bk_20120627-

       ! Evaluate boundary conditions for sediment-water column exchange.
       ! Current undersaturation of bottom water: sedb(i,0) and 
       ! Approximation for new solid sediment, as from sedimentation flux:
       ! solrat(i,1)
     
       ! mz_bk_20120627+
       DO i=1,ie
       ! DO 3 i=1,ie
       ! mz_bk_20120627-
          IF (bolay(i,j) .GT. 0._dp) THEN
             undsa = silsat - powtra(ipowasi)%ptr(i,j,1)
             sedb1(i,0)  = bolay(i,j)                                         &
                  &      * (silsat - ocetra(isilica)%ptr(i,j,kbo(i,j)))       &
                  &      * bolven(i)       
             solrat(i,1) = (sedlay(issssil)%ptr(i,j,1)                        &
                  &      + silpro(i,j) / (porsol(1)*seddw(1)))                &
                  &      * dissot / (1._dp + dissot*undsa)*porsol(1)/porwat(1)
          ENDIF
       ! mz_bk_20120627+
       ! 3     CONTINUE
       ENDDO
       ! mz_bk_20120627-

       ! Evaluate sediment undersaturation and degradation.
       ! Current undersaturation in pore water: sedb(i,k) and
       ! Approximation for new solid sediment, as from degradation:
       ! solrat(i,k) 

       ! mz_bk_20120627+
       DO k=1,ks
          DO i=1,ie
       ! DO 2 k=1,ks
       ! DO 2 i=1,ie
       ! mz_bk_20120627-
             IF (bolay(i,j) .GT. 0._dp) THEN
                undsa = silsat - powtra(ipowasi)%ptr(i,j,k)
                sedb1(i,k) = seddw(k) * porwat(k)                             &
                     &     * (silsat - powtra(ipowasi)%ptr(i,j,k))
                IF(k .GT. 1) solrat(i,k) = sedlay(issssil)%ptr(i,j,k)         &
                     &                   * dissot / (1._dp + dissot * undsa)  &
                     &                   * porsol(k) / porwat(k)
             ENDIF
       ! mz_bk_20120627+
       ! 2     CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-

       ! Solve for new undersaturation sediso, from current undersaturation
       ! sedb1, and first guess of new solid sediment solrat.     

       CALL hamocc_powadi(j,ie,solrat,sedb1,sediso,bolven)

       ! Update water column silicate, and store the flux for budget.
       ! Add biogenic opal flux to top sediment layer.

       ! mz_bk_20120627+
       DO i=1,ie
       ! DO 4 i=1,ie
       ! mz_bk_20120627-
          IF (bolay(i,j) .GT. 0._dp) THEN
!!$!!               sedfluxo(i,j,ipowasi)=sedfluxo(i,j,ipowasi) +                &
!!$!!           &  (silsat-sediso(i,0)-ocetra(isilica)%ptr(i,j,kbo(i,j)))*bolay(i,j)

             ocetra(isilica)%ptr(i,j,kbo(i,j)) = silsat - sediso(i,0)
             sedlay(issssil)%ptr(i,j,1) = sedlay(issssil)%ptr(i,j,1)          &
                  &                     + silpro(i,j) / (porsol(1) * seddw(1))
             ! mz_bk_20120627+
             silpro(i,j) = 0._dp
             ! mz_bk_20120627-
          ENDIF
       ! mz_bk_20120627+
       ! 4     CONTINUE
       ENDDO
       ! mz_bk_20120627-

       ! Calculate updated degradation rate from updated undersaturation.
       ! Calculate new solid sediment.
       ! Update pore water concentration from new undersaturation.

       ! mz_bk_20120627+
       DO k=1,ks
          DO i=1,ie
       ! DO 5 k=1,ks
       ! DO 5 i=1,ie
       ! mz_bk_20120627-
             IF (bolay(i,j) .GT. 0._dp) THEN
                ! mz_bk_20120627+
                ! umfa=porsol(k)/porwat(k)
                ! mz_bk_20120627-
                solrat(i,k) = sedlay(issssil)%ptr(i,j,k)                      &
                     &      * dissot / (1._dp + dissot * sediso(i,k))
                posol = sediso(i,k) * solrat(i,k)
                sedlay(issssil)%ptr(i,j,k)= sedlay(issssil)%ptr(i,j,k) - posol
                powtra(ipowasi)%ptr(i,j,k)= silsat - sediso(i,k)
             ENDIF
       ! mz_bk_20120627+
       ! 5     CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-

       ! CALCULATE OXYGEN-POC CYCLE AND SIMULTANEOUS OXYGEN DIFFUSION
       !*************************************************************

       ! mz_bk_20120627+
       ! ! Degradation rate constant of POP (disso) [1/(kmol O2/m3)*1/sec]
       !   disso=0.01/86400.  !  disso=3.e-5 was quite high
       !   dissot=disso*dtbgc
       ! mz_bk_20120627-

       ! This scheme is not based on undersaturation, but on O2 itself

       ! Evaluate boundary conditions for sediment-water column exchange.
       ! Current concentration of bottom water: sedb(i,0) and 
       ! Approximation for new solid sediment, as from sedimentation flux:
       ! solrat(i,1)

       ! mz_bk_20120627+
       DO i=1,ie
       ! DO 13 i=1,ie
       ! mz_bk_20120627-
          IF (bolay(i,j) .GT. 0._dp) THEN
             undsa = powtra(ipowaox)%ptr(i,j,1)
             sedb1(i,0) = bolay(i,j) * ocetra(ioxygen)%ptr(i,j,kbo(i,j))      &
                  &     * bolven(i)       
             solrat(i,1)= (sedlay(issso12)%ptr(i,j,1)                         &
                  &     + prorca(i,j) / (porsol(1)*seddw(1)))                 &
                  &     * ro2ut * dissot1 / (1._dp + dissot1 * undsa)         &
                  &     * porsol(1) / porwat(1)
          ENDIF
       ! mz_bk_20120627+
       ! 13    CONTINUE
       ENDDO
       ! mz_bk_20120627-

       ! Evaluate sediment concentration and degradation.
       ! Current concentration in pore water: sedb(i,k) and
       ! Approximation for new solid sediment, as from degradation:
       ! solrat(i,k) 

       ! mz_bk_20120627+
       DO k=1,ks
          DO i=1,ie
       ! DO 12 k=1,ks
       ! DO 12 i=1,ie
       ! mz_bk_20120627-
             IF (bolay(i,j) .GT. 0._dp) THEN
                undsa = powtra(ipowaox)%ptr(i,j,k)
                sedb1(i,k) = seddw(k) * porwat(k) * powtra(ipowaox)%ptr(i,j,k)
                IF (k .GT. 1) solrat(i,k) = sedlay(issso12)%ptr(i,j,k)        &
                     &            * ro2ut * dissot1 / (1._dp + dissot1*undsa) &
                     &            * porsol(k) / porwat(k)
             ENDIF
       ! mz_bk_20120627+
       ! 12     CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-

       ! Solve for new O2 concentration sediso, from current concentration
       ! sedb1, and first guess of new solid sediment solrat.

       CALL hamocc_powadi(j,ie,solrat,sedb1,sediso,bolven)

       ! Update water column oxygen, and store the flux for budget (opwflux). 
       ! js: opwflux not in present model code
       ! Add organic carbon flux 'prorca' to top sediment layer.

       ! mz_bk_20120627+
       DO i=1,ie
       ! DO 14 i=1,ie
       ! mz_bk_20120627-
          IF (bolay(i,j) .GT. 0._dp) THEN
             ocetra(ioxygen)%ptr(i,j,kbo(i,j)) = sediso(i,0)
             sedlay(issso12)%ptr(i,j,1) = sedlay(issso12)%ptr(i,j,1)          &
                  &                     + prorca(i,j) / (porsol(1) * seddw(1))
             prorca(i,j)=0._dp
          ENDIF
       ! mz_bk_20120627+
       ! 14    CONTINUE
       ENDDO
       ! mz_bk_20120627-

       ! Calculate updated degradation rate from updated concentration.
       ! Calculate new solid sediment.
       ! Update pore water concentration.
       ! Store flux in array aerob, for later computation of DIC and
       ! alkalinity.

       ! mz_bk_20120627+
       DO k=1,ks
       ! DO 15 k=1,ks
       ! mz_bk_20120627-
          ! mz_bk_20120627+
          umfa = porsol(k) / porwat(k)
          ! mz_bk_20120627-
          ! mz_bk_20120627+
          DO i=1,ie
          ! DO 15 i=1,ie
          ! mz_bk_20120627-
             IF (bolay(i,j) .GT. 0._dp) THEN
                ! mz_bk_20120627+
                ! umfa=porsol(k)/porwat(k)
                ! mz_bk_20120627-
                solrat(i,k) = sedlay(issso12)%ptr(i,j,k)                      &
                     &      * dissot1 / (1._dp + dissot1 * sediso(i,k))
                posol = sediso(i,k) * solrat(i,k)
                aerob(i,k) = posol * umfa             ! this has P units:
                !                                     ! kmol P/m3 of pore water
                sedlay(issso12)%ptr(i,j,k) = sedlay(issso12)%ptr(i,j,k)       &
                     &                     - posol
                powtra(ipowaph)%ptr(i,j,k) = powtra(ipowaph)%ptr(i,j,k)       &
                     &                     + posol * umfa
                powtra(ipowno3)%ptr(i,j,k) = powtra(ipowno3)%ptr(i,j,k)       &
                     &                     + posol * rnit * umfa
                powtra(ipowaox)%ptr(i,j,k) = sediso(i,k)
             ENDIF
       ! mz_bk_20120627+
       ! 15    CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-

       ! CALCULATE NITRATE REDUCTION UNDER ANAEROBIC CONDITIONS EXPLICITELY
       !*******************************************************************

       ! Denitrification rate constant of POP (disso) [1/sec]
       ! Store flux in array anaerob, for later computation of DIC and
       ! alkalinity.

       ! mz_bk_20120627+
       ! !ik      denit = 1.e-6*dtbgc
       !       denit = 0.01/86400. *dtbgc
       ! mz_bk_20120627-
       ! mz_bk_20120627+
       DO k=1,ks
       ! DO 124 k=1,ks
       ! mz_bk_20120627-
          ! mz_bk_20120627+
          umfa=porsol(k)/porwat(k)
          ! mz_bk_20120627-
          ! mz_bk_20120627+
          DO i=1,ie
          ! DO 124 i=1,ie
          ! mz_bk_20120627-
             if (bolay(i,j) .gt. 0._dp) then
                IF (powtra(ipowaox)%ptr(i,j,k) .LT. 1.e-6_dp) THEN
                   ! mz_bk_20120627+
                   posol = denit * MIN(0.5_dp * powtra(ipowno3)%ptr(i,j,k)    &
                        &              /nitdem, sedlay(issso12)%ptr(i,j,k))
                   ! posol=denit*MIN(0.5*powtra(ipowno3)%ptr(i,j,k)/114.,     &
                   !      &                    sedlay(issso12)%ptr(i,j,k))
                   ! mz_bk_20120627-
                   ! mz_bk_20120627+
                   ! umfa=porsol(k)/porwat(k)
                   ! mz_bk_20120627-
                   anaerob(i,k) = posol * umfa        ! this has P units:
                   !                                  ! kmol P/m3 of pore water
                   sedlay(issso12)%ptr(i,j,k) = sedlay(issso12)%ptr(i,j,k)    &
                        &                     - posol
                   powtra(ipowaph)%ptr(i,j,k) = powtra(ipowaph)%ptr(i,j,k)    &
                        &                     + posol * umfa
                   ! mz_bk_20120627+
                   ! powtra(ipowno3)%ptr(i,j,k) = powtra(ipowno3)%ptr(i,j,k)  &
                   !      &                     - 98._dp * posol * umfa
                   ! powtra(ipown2)%ptr(i,j,k)  = powtra(ipown2)%ptr(i,j,k)   &
                   !      &                     + 57._dp * posol * umfa
                   !tk27012010 changed no3 use and n2 production in
                   ! denitrification
                   powtra(ipowno3)%ptr(i,j,k) = powtra(ipowno3)%ptr(i,j,k)    &
                        &                     - nitdem * posol * umfa
                   powtra(ipown2)%ptr(i,j,k)  = powtra(ipown2)%ptr(i,j,k)     &
                        &                     + n2prod * posol * umfa
                   ! mz_bk_20120627-
                   ! mz_bk_20120627+
                   ! powh2obud(i,j,k) = powh2obud(i,j,k)                      &
                   !      &           + 0.5_dp * n2prod * posol * umfa
                   ! mz_bk_20120627-     
                ENDIF   ! oxygen <1.e-6
             endif   ! bolay
       ! mz_bk_20120627+
       ! 124   CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-

       ! mz_bk_20120627+
       ! new in mpiesm-1.0.00 version of hamocc
       ! sulphate reduction in sediments
       DO k=1,ks
       !  DO 125 k=1,ks
          umfa = porsol(k) / porwat(k)
          DO i=1,ie
          !  DO 125 i=1,kpie
             if (bolay(i,j) .gt. 0._dp) then
                IF (powtra(ipowaox)%ptr(i,j,k) .LT. 1.e-6_dp) THEN
                   posol = denit * 0.01_dp * sedlay(issso12)%ptr(i,j,k)
                   ansulf(i,k) = posol * umfa         ! this has P units:
                   !                                  ! kmol P/m3 of pore water
                   sedlay(issso12)%ptr(i,j,k) = sedlay(issso12)%ptr(i,j,k)    &
                        &                     - posol
                   powtra(ipowaph)%ptr(i,j,k) = powtra(ipowaph)%ptr(i,j,k)    &
                        &                     + posol * umfa
                   powtra(ipowno3)%ptr(i,j,k) = powtra(ipowno3)%ptr(i,j,k)    &
                        &                     + posol * umfa * rno3
                   !  powh2obud(i,j,k) = powh2obud(i,j,k) - ro2ut * posol * umfa
                ENDIF
             endif
       !  125   CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-

       ! CALCULATE CaCO3-CO3 CYCLE AND SIMULTANEOUS CO3-UNDERSATURATION
       ! DIFFUSION
       !*********************************************************************

       ! COMPUTE NEW POWCAR=CARBONATE ION CONCENTRATION IN THE SEDIMENT
       ! FROM CHANGED ALKALINITY (NITRATE PRODUCTION DURING REMINERALISATION)
       ! AND DIC GAIN. ITERATE 5 TIMES. THIS CHANGES PH (SEDHPL) OF SEDIMENT.

       ! mz_bk_20120627+
       DO iter=1,5
       ! DO 10 ITER=1,5
       ! mz_bk_20120627-

          ! mz_bk_20120627+
          ! hconve=0.
          ! mz_bk_20120627-
          ! mz_bk_20120627+
          DO k=1,ks
             DO i=1,ie
          ! DO 1 K=1,KS
          ! DO 1 i=1,ie
          ! mz_bk_20120627-
                IF (bolay(i,j) .GT. 0._dp) THEN
                   bt = rrrcl * sao(i,j,kbo(i,j))
                   ! mz_bk_20120627+
                   ! alkalinity is increased during denitrification due to
                   ! consumption of H+ via NO3 (see Wolf-Gladrow etal,2007)
                   alk    = powtra(ipowaal)%ptr(i,j,k)                        &
                        & - (ansulf(i,k) + aerob(i,k)) * rnit                 &
                        & + nitdem * anaerob(i,k)
                   ! alk    = powtra(ipowaal)%ptr(i,j,k)                      &
                   !      & -(anaerob(i,k) + aerob(i,k)) * 16._dp
                   c      = powtra(ipowaic)%ptr(i,j,k)                        &
                        & + (anaerob(i,k) + aerob(i,k) + ansulf(i,k)) * rcar
                   ! c      = powtra(ipowaic)%ptr(i,j,k)                      &
                   !      & + (anaerob(i,k) + aerob(i,k)) * 122._dp
                   ! mz_bk_20120627-
                   ak1 = ak13(i,j,kbo(i,j))
                   ak2 = ak23(i,j,kbo(i,j))
                   akb = akb3(i,j,kbo(i,j))
                   akw = akw3(i,j,kbo(i,j))
                   h   = sedhpl(i,j,k)
                   t1  = h / ak1
                   t2  = h / ak2
                   a      = c * (2._dp + t2)/(1._dp + t2 + t2*t1) + akw/h     &
                        & - h + bt / (1._dp + h/akb) - alk
                   dadh   = c * (1._dp / (ak2 * (1._dp + t2 + t2 * t1))       &
                        & - (2._dp + t2) * (1._dp / ak2 + 2._dp * t1 / ak2)   &
                        & / (1._dp + t2 + t2 * t1)**2)                        &
                        & - akw/h**2 - 1._dp - (bt / akb)/(1._dp + h/akb)**2
                   dddhhh = a / dadh
                   ! mz_bk_20120627+
                   sedhpl(i,j,k) = MAX(h - dddhhh, 1.e-11_dp)
                   ! reduk=MAX(1._dp,2.*abs(dddhhh/h))
                   ! sedhpl(i,j,k)=h-dddhhh/reduk
                   ! hconve=hconve+dddhhh**2
                   ! mz_bk_20120627-
                   powcar(i,k) = c / (1._dp + t2 * (1._dp + t1))
                ENDIF

       ! mz_bk_20120627+
       ! 1     CONTINUE
             ENDDO       ! loop over i
          ENDDO          ! loop over k
       ! 10    CONTINUE
       ENDDO             ! loop over iter
       ! mz_bk_20120627-

       ! mz_bk_20120627+
       ! !Dissolution rate constant of CaCO3 (disso) [1/(kmol CO3--/m3)*1/sec]
       !       disso=1.e-7
       !       dissot=disso*dtbgc
       ! mz_bk_20120627-

       ! Evaluate boundary conditions for sediment-water column exchange.
       ! Current undersaturation of bottom water: sedb(i,0) and 
       ! Approximation for new solid sediment, as from sedimentation flux:
       ! solrat(i,1)

       ! CO3 saturation concentration is aksp/calcon as in CARCHM 
       ! (calcon defined in BELEG_BGC with 1.03e-2; 1/calcon =~ 97.) 

       ! mz_bk_20120627+
       DO i=1,ie
       ! DO 23 i=1,ie
       ! mz_bk_20120627-
          IF (bolay(i,j) .GT. 0._dp) THEN
             satlev = aksp(i,j,kbo(i,j)) / calcon + 2.e-5_dp
             undsa  = MAX(satlev - powcar(i,1), 0._dp)
             sedb1(i,0) = bolay(i,j) * (satlev - co3(i,j,kbo(i,j)))           &
                  &     * bolven(i)
             solrat(i,1)= (sedlay(isssc12)%ptr(i,j,1)                         &
                  &     + prcaca(i,j) / (porsol(1) * seddw(1)))               &
                  &     * dissot2/(1._dp + dissot2*undsa)*porsol(1)/porwat(1)
          ENDIF
       ! mz_bk_20120627+
       ! 23     CONTINUE
       ENDDO
       ! mz_bk_20120627-

       ! Evaluate sediment undersaturation and degradation.
       ! Current undersaturation in pore water: sedb(i,k) and
       ! Approximation for new solid sediment, as from degradation:
       ! solrat(i,k) 

       ! mz_bk_20120627+
       DO k=1,ks
          DO i=1,ie
       ! DO 22 k=1,ks
       ! DO 22 i=1,ie
       ! mz_bk_20120627-
             IF (bolay(i,j) .GT. 0._dp) THEN
                undsa = MAX(aksp(i,j,kbo(i,j)) / calcon - powcar(i,k), 0._dp)
                sedb1(i,k) = seddw(k) * porwat(k) * undsa
                IF (k .GT. 1) solrat(i,k) = sedlay(isssc12)%ptr(i,j,k)        &
                     &                    * dissot2/(1._dp + dissot2*undsa)   &
                     &                    * porsol(k) / porwat(k)
                IF (undsa .LE. 0._dp) solrat(i,k) = 0._dp
             ENDIF
       ! mz_bk_20120627+
       ! 22     CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-

       ! Solve for new undersaturation sediso, from current undersaturation
       ! sedb1, and first guess of new solid sediment solrat.     

       CALL hamocc_powadi(j,ie,solrat,sedb1,sediso,bolven)

       ! There is no exchange between water and sediment with respect to co3
       ! so far.
       ! Add calcite flux 'prcaca' to uppermost sediment layer.
       ! mz_bk_20120627+
       DO i=1,ie
       ! DO 24 i=1,ie
       ! mz_bk_20120627-
          IF (bolay(i,j) .GT. 0._dp) THEN
             sedlay(isssc12)%ptr(i,j,1)= sedlay(isssc12)%ptr(i,j,1)           &
                  &                    + prcaca(i,j) / (porsol(1) * seddw(1))
             prcaca(i,j) = 0._dp
          ENDIF
       ! mz_bk_20120627+
       ! 24    CONTINUE
       ENDDO
       ! mz_bk_20120627-

       ! Calculate updated degradation rate from updated undersaturation.
       ! Calculate new solid sediment.
       ! No update of powcar pore water concentration from new undersaturation
       ! so far.
       ! Instead, only update DIC, and, of course, alkalinity.
       ! This also includes gains from aerobic and anaerobic decomposition.

       ! mz_bk_20120627+
       DO k=1,ks
       ! DO 25 k=1,ks
       ! mz_bk_20120627-
          ! mz_bk_20120627+
          umfa=porsol(k)/porwat(k)
          ! mz_bk_20120627-
          ! mz_bk_20120627+
          DO i=1,ie
          ! DO 25 i=1,ie
          ! mz_bk_20120627-
             IF (bolay(i,j) .GT. 0._dp) THEN
                ! mz_bk_20120627+
                ! umfa=porsol(k)/porwat(k)
                ! mz_bk_20120627-
                solrat(i,k) = sedlay(isssc12)%ptr(i,j,k)                      &
                     &      * dissot2 / (1._dp + dissot2 * sediso(i,k))
                posol = sediso(i,k) * solrat(i,k)
                sedlay(isssc12)%ptr(i,j,k) = sedlay(isssc12)%ptr(i,j,k)       &
                     &                     - posol
                ! mz_bk_20120627+
                ! powtra(ipowaic)%ptr(i,j,k) = powtra(ipowaic)%ptr(i,j,k)     &
                !      &          +posol*umfa+(aerob(i,k)+anaerob(i,k))*122._dp
                ! powtra(ipowaal)%ptr(i,j,k) = powtra(ipowaal)%ptr(i,j,k)     &
                !      &     +2._dp*posol*umfa-16._dp*(aerob(i,k)+anaerob(i,k))
                powtra(ipowaic)%ptr(i,j,k) = powtra(ipowaic)%ptr(i,j,k)       &
                     &                     + posol * umfa                     &
                     &                     + (aerob(i,k) + anaerob(i,k)       &
                     &                     + ansulf(i,k)) * rcar
                powtra(ipowaal)%ptr(i,j,k) = powtra(ipowaal)%ptr(i,j,k)       &
                     &                     + 2._dp * posol * umfa             &
                     &                     - rnit                             &
                     &                     * (aerob(i,k) + ansulf(i,k))       &
                     &                     + nitdem * anaerob(i,k)
                !  pown2bud(i,j,k) = pown2bud(i,j,k) + 2._dp*n2prod*anaerob(i,k)
                ! include iron 25102010
                ! powtra(ipowafe)%ptr(i,j,k) = powtra(ipowafe)%ptr(i,j,k)     &
                !      &                     + (aerob(i,k) + anaerob(i,k)     &
                !      &                     + ansulf(i,k)) * riron
                ! mz_bk_20120627-
             ENDIF
       ! mz_bk_20120627+
       ! 25    CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-


    ! mz_bk_20120627+
    ! 8888  CONTINUE
    ENDDO
    ! mz_bk_20120627-

    CALL hamocc_dipowa(ie,je,ke)

    !ik add clay sedimentation
    !ik this is currently assumed to depend on total and corg sedimentation:
    !ik f(POC) [kg C] / f(total) [kg] = 0.05
    !ik thus it is 
    do j=1,je
       do i=1,ie
          IF (bolay(i,j).GT.0._dp) THEN
             sedlay(issster)%ptr(i,j,1) = sedlay(issster)%ptr(i,j,1)          &
                  &                     + produs(i,j) / (porsol(1) * seddw(1))
          END IF
       enddo
    enddo

    ! mz_bk_20120627+
    DO j=1,je
       DO i=1,ie
    ! DO 91 j=1,je
    ! DO 91 i=1,ie
    ! mz_bk_20120627-
          silpro(i,j) = 0._dp
          prorca(i,j) = 0._dp
          prcaca(i,j) = 0._dp
          produs(i,j) = 0._dp
    ! mz_bk_20120627+
    ! 91     CONTINUE
       ENDDO
    ENDDO
    ! mz_bk_20120627-

  END SUBROUTINE hamocc_powach
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  SUBROUTINE hamocc_dipowa(ie,je,ke)
    !
    ! Correspondes to dipowa.f90
    !
    !     Purpose
    !     -------
    !     calculate vertical diffusion of sediment pore water properties
    !     and diffusive flux through the ocean/sediment interface.
    !     integration.
    !
    !     Method
    !     -------
    !     implicit formulation;
    !     constant diffusion coefficient : 1.e-9 set in BODENSED.
    !     diffusion coefficient : zcoefsu/zcoeflo for upper/lower
    !     sediment layer boundary.
    !
    ! Original:
    !     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    !     Modified
    !     --------
    !     S.Legutke,        *MPI-MaD, HH*    10.04.01
    !     - all npowtra-1 properties are diffused in 1 go.
    !     js: not mass conserving check c13/powtra/ocetra
    !-------------------------------------------------------------------------
    implicit none

    INTEGER :: i,j,k,l,iv
    INTEGER :: ie,je,ke
    integer :: iv_oc                           ! index of ocetra in powtra loop

    REAL(DP):: sedb1(ie,0:ks,npowtra)          ! ????
    REAL(DP):: zcoefsu(0:ks),zcoeflo(0:ks)     ! diffusion coefficients
    !                                          ! (upper/lower)
    ! mz_bk_20120627+
    ! REAL(DP):: TREDSY(ie,0:ke,3)               ! redsy for 'reduced system'?
    REAL(DP):: TREDSY(ie,0:ks,3)               ! redsy for 'reduced system'?
    ! mz_bk_20120627-

    REAL(DP):: aprior           ! start value of oceanic tracer in bottom layer

    !ik accelerated sediment
    !ik needed for boundary layer ventilation in fast sediment routine      

    REAL(DP):: bolven(ie)                      ! bottom layer ventilation rate

    zcoefsu(0) = 0.0_dp
    DO  k=1,ks
       ! sediment diffusion coefficient * 1/dz * fraction of pore water
       ! at half depths
       zcoefsu(k)   = -sedict * seddzi(k) * porwah(k)
       ! the lowerdiffusive flux of layer k is identical to the upper diff.
       ! flux of layer k+1
       zcoeflo(k-1) = -sedict * seddzi(k) * porwah(k)    
    ENDDO
    ! diffusion coefficient for bottom sediment layer
    zcoeflo(ks)=0.0_dp

    ! mz_bk_20120627+
    DO j=1,je
    ! DO 11000 j=1,je
    ! mz_bk_20120627-

       ! calculate bottom ventilation rate for scaling of
       ! sediment-water exchange
       ! mz_bk_20120627+
       do i=1,ie      
       ! do 1170 i=1,ie      
       ! mz_bk_20120627-
          bolven(i) = 1._dp
       ! mz_bk_20120627+
       ! 1170  continue
       enddo
       ! mz_bk_20120627-


       k=0
       ! mz_bk_20120627+
       DO i=1,ie
       ! DO 1421 i=1,ie
       ! mz_bk_20120627-
          tredsy(i,k,1) = zcoefsu(k)
          tredsy(i,k,3) = zcoeflo(k)
          tredsy(i,k,2) = bolven(i)*bolay(i,j) - tredsy(i,k,1) - tredsy(i,k,3)
          !                            dz(kbo) - diff upper    - diff lower
       ! mz_bk_20120627+
       ! 1421  CONTINUE
       ENDDO
       ! mz_bk_20120627-

       k=0
       ! mz_bk_20120627+
       DO iv=1,npowtra      ! loop over pore water tracers
       ! DO 1422 iv=1,npowtra      ! loop over pore water tracers
       ! mz_bk_20120627-
          iv_oc=iv
          ! mz_bk_20120627+
          ! if(iv .eq. ipowafe) iv_oc = iiron
          ! mz_bk_20120627-
          ! mz_bk_20120627+
          DO i=1,ie
          ! DO 1422 i=1,ie
          ! mz_bk_20120627-
             sedb1(i,k,iv) = 0._dp
             IF(bolay(i,j) .GT. 0._dp)                                       &
                  &      sedb1(i,k,iv) = ocetra(iv_oc)%ptr(i,j,kbo(i,j))     &
                  &                    * bolay(i,j) * bolven(i)
             ! tracer_concentration(kbo) * dz(kbo)
       ! mz_bk_20120627+
       ! 1422 CONTINUE
          ENDDO     ! loop over i
       ENDDO        ! loop over iv
       ! mz_bk_20120627-

       ! mz_bk_20120627+
       DO k=1,ks
          DO i=1,ie     
       ! DO 1321 k=1,ks
       !    DO 1321 i=1,ie     
       ! mz_bk_20120627-
             tredsy(i,k,1) = zcoefsu(k)
             tredsy(i,k,3) = zcoeflo(k)
             tredsy(i,k,2) = seddw(k)*porwat(k) - tredsy(i,k,1) - tredsy(i,k,3)
       ! mz_bk_20120627+
       ! 1321    CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-

       ! mz_bk_20120627+
       DO iv=1,npowtra
          DO k=1,ks
             DO i=1,ie     
       ! DO 1322 iv=1,npowtra
       ! DO 1322 k=1,ks
       ! DO 1322 i=1,ie     
       ! mz_bk_20120627-
                sedb1(i,k,iv) = powtra(iv)%ptr(i,j,k) * porwat(k) * seddw(k) 
                ! tracer_concentration(k[1:ks])*porewater fraction(k)*dz(k)
       ! mz_bk_20120627+
       ! 1322 CONTINUE
             ENDDO     ! loop over i
          ENDDO        ! loop over k
       ENDDO           ! loop over iv
       ! mz_bk_20120627-

       ! mz_bk_20120627+
       DO k=1,ks
          DO i=1,ie
       ! DO 132 k=1,ks
       !    DO 133 i=1,ie
       ! mz_bk_20120627-
             IF (bolay(i,j) .GT. 0._dp) THEN
                ! this overwrites tredsy(k=0) for k=1
                tredsy(i,k-1,1) = tredsy(i,k,1) / tredsy(i,k-1,2)
                !                 diff upper    / conc (k-1)
                tredsy(i,k,2) = tredsy(i,k,2)                               &
                     &        - tredsy(i,k-1,3)*tredsy(i,k,1)/tredsy(i,k-1,2)
                !concentration - diff lower    * diff upper  /conc(k-1)
             ENDIF
       ! mz_bk_20120627+
       ! 133        CONTINUE
          ENDDO
       ! 132     CONTINUE
       ENDDO
       ! mz_bk_20120627-

       ! diffusion from above
       ! mz_bk_20120627+
       DO iv=1,npowtra
          DO k=1,ks
             DO i=1,ie
       ! DO 135 iv=1,npowtra
       !    DO 135 k=1,ks
       !    DO 135 i=1,ie
       ! mz_bk_20120627-
                sedb1(i,k,iv) = sedb1(i,k,iv)                                 &
                     &        - tredsy(i,k-1,1) * sedb1(i,k-1,iv)       
       ! mz_bk_20120627+
       ! 135     CONTINUE
             ENDDO      ! loop over i
          ENDDO         ! loop over k
       ENDDO            ! loop over iv
       ! mz_bk_20120627-

       ! sediment bottom layer
       k=ks
       ! mz_bk_20120627+
       DO iv=1,npowtra
          DO i=1,ie
       ! DO 136 iv=1,npowtra
       !    DO 136 i=1,ie
       ! mz_bk_20120627-
             IF (bolay(i,j) .GT. 0._dp) THEN
                powtra(iv)%ptr(i,j,k) = sedb1(i,k,iv) / tredsy(i,k,2)
             ENDIF
       ! mz_bk_20120627+
       ! 136     CONTINUE
          ENDDO
       ENDDO
       ! mz_bk_20120627-

       ! sediment column
       ! mz_bk_20120627+
       DO iv=1,npowtra
          DO k=1,ks-1
       ! DO 137 iv=1,npowtra
       !    DO 137 k=1,ks-1
       ! mz_bk_20120627-
             l = ks - k
             ! mz_bk_20120627+
             DO i=1,ie
             ! DO 137 i=1,ie
             ! mz_bk_20120627-
                IF (bolay(i,j) .GT. 0._dp) THEN
                   powtra(iv)%ptr(i,j,l) = ( sedb1(i,l,iv)                    &
                        &       - tredsy(i,l,3) * powtra(iv)%ptr(i,j,l+1) )   &
                        &       / tredsy(i,l,2)
                ENDIF
       ! mz_bk_20120627+
       ! 137    CONTINUE
             ENDDO       ! loop over i
          ENDDO          ! loop over k
       ENDDO             ! loop over iv
       ! mz_bk_20120627-

       ! sediment ocean interface
       ! mz_bk_20120627+
       DO iv=1,npowtra         ! caution - the following assumes same indecees
          !                    ! for ocetra and powtra test npowa_base 071106
       ! DO 139 iv=1,npowtra   ! caution - the following assumes same indecees
       !                       ! for ocetra and powtra test npowa_base 071106
       ! mz_bk_20120627-
          ! check mo_param1_bgc.f90 for consistency
          iv_oc = iv

          ! mz_bk_20120627+
          ! if (iv .eq. ipowafe) iv_oc = iiron
          ! mz_bk_20120627-

          ! mz_bk_20120627+
          DO i=1,ie
          ! DO 139 i=1,ie
          ! mz_bk_20120627-
             l=0
             IF (bolay(i,j) .GT. 0._dp) THEN

                aprior = ocetra(iv_oc)%ptr(i,j,kbo(i,j))
                ocetra(iv_oc)%ptr(i,j,kbo(i,j)) =                             &
                     &        ( sedb1(i,l,iv)                                 &
                     &          - tredsy(i,l,3) * powtra(iv)%ptr(i,j,l+1) )   &
                     &        / tredsy(i,l,2)

!!$!!                  sedfluxo(i,j,iv)=sedfluxo(i,j,iv)                          &    !used in inventory_bgc/maschk (diagnostics) 
!!$!!            &                     +ocetra(iv)%ptr(i,j,kbo(i,j))-aprior
      
             ENDIF
       ! mz_bk_20120627+
       ! 139  CONTINUE
          ENDDO      ! loop over i
       ENDDO         ! loop over iv
       ! mz_bk_20120627-


    ! mz_bk_20120627+
    ! 11000 CONTINUE    ! j loop
    ENDDO            ! loop over j
    ! mz_bk_20120627-

  END SUBROUTINE hamocc_dipowa
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE hamocc_powadi(j,ie,solrat,sedb1,sediso,bolven)
    !
    ! Correspondes to powadi.f90
    ! vertical diffusion with simultaneous dissolution.
    !
    !     Input  solrat : dissolution rate
    !     =====       j : zonal grid index
    !             sedb1 : tracer at entry
    !
    ! Original:
    !     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    !     Modified
    !     --------
    !     S.Legutke,        *MPI-MaD, HH*    10.04.01
    !------------------------------------------------------------------------
    implicit none

    INTEGER :: ie,i,j,k,l
    REAL(DP):: sedb1(ie,0:ks),sediso(ie,0:ks)
    REAL(DP):: solrat(ie,ks)
    REAL(DP):: bolven(ie)
    REAL(DP):: TREDSY(ie,0:ks,3)
    REAL(DP):: asu,alo

    DO k=1,ks
       asu = sedict * seddzi(k) * porwah(k)
       alo = 0._dp
       IF (k .LT. ks) alo = sedict * seddzi(k+1) * porwah(k+1)
       DO i=1,ie
          tredsy(i,k,1) = -asu
          tredsy(i,k,3) = -alo 
          tredsy(i,k,2) =  seddw(k) * porwat(k) - tredsy(i,k,1)               &
               &        - tredsy(i,k,3) + solrat(i,k) * porwat(k) * seddw(k)
       ENDDO
    ENDDO

    k = 0
    asu = 0._dp
    alo = sedict * seddzi(1) * porwah(1)
    DO i=1,ie
       IF (bolay(i,j) .GT. 0._dp) THEN
          tredsy(i,k,1) = -asu
          tredsy(i,k,3) = -alo 
          tredsy(i,k,2) = bolven(i) * bolay(i,j)                              &
               &        - tredsy(i,k,1) - tredsy(i,k,3)
       ELSE
          tredsy(i,k,1) = 0._dp
          tredsy(i,k,3) = 0._dp
          tredsy(i,k,2) = 0._dp
       ENDIF
    ENDDO

    DO k=1,ks
       DO i=1,ie
          IF (bolay(i,j) .GT. 0._dp) THEN
             tredsy(i,k-1,1) = tredsy(i,k,1) / tredsy(i,k-1,2)
             tredsy(i,k,2)   = tredsy(i,k,2)                                  &
                  &          - tredsy(i,k-1,3) * tredsy(i,k,1)/tredsy(i,k-1,2)
          ENDIF
       ENDDO
    ENDDO

    DO k=1,ks
       DO i=1,ie
          sedb1(i,k) = sedb1(i,k) - tredsy(i,k-1,1) * sedb1(i,k-1)
       ENDDO
    ENDDO

    k=ks
    DO i=1,ie
       IF (bolay(i,j) .GT. 0._dp) sediso(i,k) = sedb1(i,k) / tredsy(i,k,2)
    ENDDO

    DO k=1,ks
       l = ks - k
       DO i=1,ie
          IF (bolay(i,j) .GT. 0._dp) sediso(i,l) =                            &
               &            ( sedb1(i,l) - tredsy(i,l,3) * sediso(i,l+1) )    &
               &            / tredsy(i,l,2)
       ENDDO
    ENDDO

  END SUBROUTINE hamocc_powadi
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE hamocc_sedshi(ie,je)
    !
    ! Correspondes to sedshi.f90
    !
    ! Original:
    !     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    !     Modified
    !     --------
    !     S.Legutke,        *MPI-MaD, HH*    10.04.01
    !     - rename ssssil(i,j,k)=sedlay(i,j,k,issssil) etc.
    !     I. Kriest         *MPI-Met, HH*,   27.05.03
    !     - change specific weights for opal, CaCO3, POC --> bodensed.f90
    !     - include upward transport
    !------------------------------------------------------------------------
    implicit none

    INTEGER :: ie,je,i,j,k,l,iv
    REAL(DP):: wsed(ie,je)      ! shifting velocity for upward/downward shifts
    REAL(DP):: fulsed(ie,je)
    REAL(DP):: sedlo,uebers
    REAL(DP):: seddef           ! sediment deficiency
    REAL(DP):: spresent, buried
    REAL(DP):: refill,frac

    ! DOWNWARD SHIFTING
    ! shift solid sediment downwards, if layer is full, i.e., if 
    ! the volume filled by the four constituents poc, opal, caco3, clay
    ! is more than porsol*seddw
    ! the outflow of layer i is given by sedlay(i)*porsol(i)*seddw(i), it is
    ! distributed in the layer below over a volume of porsol(i+1)*seddw(i+1)

    !e    write(0,121)orgfa,calfa,oplfa,clafa
    !121   format('sediments',4e12.4)

    do k=1,ks-1

       do j=1,je
          do i=1,ie
             if (bolay(i,j) .gt. 0._dp) then
                sedlo  = orgfa * rcar * sedlay(issso12)%ptr(i,j,k)            &
                     & + calfa * sedlay(isssc12)%ptr(i,j,k)                   &
                     & + oplfa * sedlay(issssil)%ptr(i,j,k)                   &
                     & + clafa * sedlay(issster)%ptr(i,j,k)
                ! "full" sediment has sedlo=1. for sedlo>1., wsed is >0.
                ! downward shifting velocity (?)
                wsed(i,j) = MAX(0._dp, (sedlo-1._dp)/(sedlo + 1.e-10_dp))
             endif
          enddo !end i-loop
       enddo !end j-loop

       ! filling downward  (accumulation)
       do iv=1,nsedtra
          do j=1,je
             do i=1,ie
                if (bolay(i,j) .gt. 0._dp) then
                   ! 'uebersaettigung?'
                   uebers = wsed(i,j) * sedlay(iv)%ptr(i,j,k)
                   sedlay(iv)%ptr(i,j,k)   = sedlay(iv)%ptr(i,j,k) - uebers
                   sedlay(iv)%ptr(i,j,k+1) = sedlay(iv)%ptr(i,j,k+1)          &
                        &      + uebers                                       &
                        &      * (seddw(k)*porsol(k))/(seddw(k+1)*porsol(k+1))
                endif
             enddo !end i-loop
          enddo !end j-loop
       enddo !end iv-loop

    enddo !end k-loop


    ! store amount lost from bottom sediment layer - this is a kind of 
    ! permanent burial in deep consolidated layer, and this stuff is 
    ! effectively lost from the whole ocean+sediment(+atmosphere) system.
    ! Would have to be supplied by river runoff or simple addition e.g. 
    ! to surface layers in the long range. Can be supplied again if a 
    ! sediment column has a deficiency in volume.

    do j=1,je
       do i=1,ie
          if (bolay(i,j) .gt. 0._dp) then
             sedlo  = orgfa * rcar * sedlay(issso12)%ptr(i,j,ks)              &
                  & + calfa * sedlay(isssc12)%ptr(i,j,ks)                     &
                  & + oplfa * sedlay(issssil)%ptr(i,j,ks)                     &
                  & + clafa * sedlay(issster)%ptr(i,j,ks)
             wsed(i,j) = MAX(0._dp, (sedlo - 1._dp) / (sedlo + 1.e-10_dp))
          endif
       enddo !end i-loop
    enddo !end j-loop

    do iv=1,nsedtra
       do j=1,je
          do i=1,ie
             if (bolay(i,j) .gt. 0._dp) then
                uebers = wsed(i,j) * sedlay(iv)%ptr(i,j,k)
                sedlay(iv)%ptr(i,j,ks) = sedlay(iv)%ptr(i,j,ks) - uebers
                burial(iv)%ptr(i,j)    = burial(iv)%ptr(i,j)                  &
                     &                 + uebers * seddw(k) * porsol(k)
             endif
          enddo !end i-loop
       enddo !end j-loop
    enddo !end iv-loop

    return !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< switch off upward shift

    ! now the loading nowhere exceeds 1.

    ! digging from below in case of erosion (js and initialization with 0 or
    ! none, as currently)
    ! UPWARD SHIFTING
    ! shift solid sediment upwards, if total sediment volume is less
    ! than required, i.e., if the volume filled by the four constituents 
    ! poc, opal, caco3, clay (integrated over the total sediment column)
    ! is less than porsol*seddw (integrated over the total sediment column)
    ! first, the deepest box is filled from below with total required volume; 
    ! then, successively, the following layers are filled upwards.
    ! if there is not enough solid matter to fill the column, add clay.
    ! (js: so implicite initial state is all clay)

    do j=1,je
       do i=1,ie
          fulsed(i,j) = 0._dp
       enddo !end i-loop
    enddo !end j-loop

    ! determine how the total sediment column is filled 
    do k=1,ks
       do j=1,je
          do i=1,ie
             if (bolay(i,j) .gt. 0._dp) then
                sedlo  = orgfa * rcar * sedlay(issso12)%ptr(i,j,k)            &
                     & + calfa * sedlay(isssc12)%ptr(i,j,k)                   &
                     & + oplfa * sedlay(issssil)%ptr(i,j,k)                   &
                     & + clafa * sedlay(issster)%ptr(i,j,k)
                fulsed(i,j) = fulsed(i,j) + porsol(k) * seddw(k) * sedlo
             endif
          enddo !end i-loop
       enddo !end j-loop
    enddo !end k-loop

    ! shift the sediment deficiency from the deepest (burial) 
    ! layer into layer ks
    do j=1,je
       do i=1,ie
          if (bolay(i,j) .gt. 0._dp) then

             ! deficiency with respect to fully loaded sediment | packed in
             ! sedlay(i,j,ks) ??
             ! this is the volume of sediment shifted upwards from the
             ! burial layer

             ! 'sediment deficiency', solfu = total column inegrated solid
             ! fraction volume (bodensed)
             seddef = solfu - fulsed(i,j)    

             ! total volume of solid constituents in buried layer
             spresent = orgfa * rcar * burial(issso12)%ptr(i,j)               &
                  &   + calfa * burial(isssc12)%ptr(i,j)                      &
                  &   + oplfa * burial(issssil)%ptr(i,j)                      &
                  &   + clafa * burial(issster)%ptr(i,j)

             ! determine whether an additional amount of clay is needed from
             ! the burial layer to fill the whole sediment; I assume that
             ! there is an infinite supply of clay from below
             burial(issster)%ptr(i,j) = burial(issster)%ptr(i,j)              &
                  &                   + MAX(0._dp, seddef - spresent) / clafa

             ! determine new volume of buried layer
             buried = orgfa * rcar * burial(issso12)%ptr(i,j)                 &
                  & + calfa * burial(isssc12)%ptr(i,j)                        &
                  & + oplfa * burial(issssil)%ptr(i,j)                        &
                  & + clafa * burial(issster)%ptr(i,j)

             ! fill the deepest active sediment layer
             refill = seddef / buried 
             frac = porsol(ks) * seddw(ks) !changed k to ks, ik

             sedlay(issso12)%ptr(i,j,ks) = sedlay(issso12)%ptr(i,j,ks)        &
                  &                  + refill * burial(issso12)%ptr(i,j) / frac
             sedlay(isssc12)%ptr(i,j,ks) = sedlay(isssc12)%ptr(i,j,ks)        &
                  &                  + refill * burial(isssc12)%ptr(i,j) / frac
             sedlay(issssil)%ptr(i,j,ks) = sedlay(issssil)%ptr(i,j,ks)        &
                  &                  + refill * burial(issssil)%ptr(i,j) / frac
             sedlay(issster)%ptr(i,j,ks) = sedlay(issster)%ptr(i,j,ks)        &
                  &                  + refill * burial(issster)%ptr(i,j) / frac

             ! account for losses in buried sediment
             burial(issso12)%ptr(i,j) = burial(issso12)%ptr(i,j)              &
                  &                   - refill * burial(issso12)%ptr(i,j)
             burial(isssc12)%ptr(i,j) = burial(isssc12)%ptr(i,j)              &
                  &                   - refill * burial(isssc12)%ptr(i,j)
             burial(issssil)%ptr(i,j) = burial(issssil)%ptr(i,j)              &
                  &                   - refill * burial(issssil)%ptr(i,j)
             burial(issster)%ptr(i,j) = burial(issster)%ptr(i,j)              &
                  &                   - refill * burial(issster)%ptr(i,j)

          endif ! bolay >0
       enddo !end i-loop
    enddo !end j-loop

    ! redistribute overload of deepest layer ks to layers 2 to ks
    do k=ks,2,-1
       do j=1,je
          do i=1,ie
             if (bolay(i,j) .gt. 0._dp) then
                sedlo  = orgfa * rcar * sedlay(issso12)%ptr(i,j,k)            &
                     & + calfa * sedlay(isssc12)%ptr(i,j,k)                   &
                     & + oplfa * sedlay(issssil)%ptr(i,j,k)                   &
                     & + clafa * sedlay(issster)%ptr(i,j,k)
                wsed(i,j) = MAX(0._dp, (sedlo - 1._dp) / (sedlo + 1.e-10_dp))
             endif
          enddo !end i-loop
       enddo !end j-loop

       do iv=1,4
          do j=1,je
             do i=1,ie
                if (bolay(i,j) .gt. 0._dp) then
                   uebers = sedlay(iv)%ptr(i,j,k) * wsed(i,j)
                   frac = porsol(k) * seddw(k) / (porsol(k-1) * seddw(k-1))
                   sedlay(iv)%ptr(i,j,k) = sedlay(iv)%ptr(i,j,k) - uebers
                   sedlay(iv)%ptr(i,j,k-1) = sedlay(iv)%ptr(i,j,k-1)          &
                        &                  + uebers * frac
                   !                                note k-1 here = upward shift
                endif
             enddo !end i-loop
          enddo !end j-loop
       enddo !end iv-loop

    enddo  !end k-loop

  END SUBROUTINE hamocc_sedshi
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE hamocc_photoprod(ie,je,ke,co,isop,ch3i,sao, tho)
    !
    ! Calculates ocean concentrations of CO, Isoprene and CH3I
    ! changed by production via photolysis
    !------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(DP),INTENT(INOUT) :: co(:,:,:)
    REAL(DP),INTENT(INOUT) :: isop(:,:,:)
    REAL(DP),INTENT(INOUT) :: ch3i(:,:,:)
    REAL(DP),INTENT(IN) :: sao(:,:,:)
    REAL(DP),INTENT(IN) :: tho(:,:,:)
    INTEGER,INTENT(IN)  :: ie,je,ke
    REAL(DP)            :: chl(ie,je,ke)

    INTEGER :: i,j,k
    REAL(DP) :: dark, tau_c, consumption(ie,je,ke) 
    REAL(DP) :: phoprod(ie,je,ke), cdom_abs(ie,je,ke)
    REAL(DP) :: cl_react(ie,je,ke)
    REAL(DP) :: cl_conc(ie,je,ke)

    !----------------CO-------------------!
    !     See Kettle et al. (2005)        !
    !-------------------------------------!

    !DARK PRODUCTION RATE:
    ! optimized for surface CO: 0.116 [nM/day]
    ! optimized for deep water: 0.025 [nM/day]
    ! optimized for deep and surface water: 0.027 [nM/day]
    dark  = 0.027*1E-9_dp / (60._dp*60._dp*24._dp)
    !          ! 0.116 [nM/day] => 0.116 * 1E-9 /(60*60*24) [M/s] =[kmol/m^3/s]
    ! E-FOLDING TIME :
    ! optimized for surface CO: 43.4 [h]
    ! optimized for deep water CO: ... [h]
    ! optimized for deep and surface water CO: 184.1  [h]
    tau_c = 184.4*(60*60)   ! 166 [h] => (60*60)[s]
    !MICROBIAL CONSUMPTION [kmol/m^3/s]
    consumption(:,:,:) = (-1._dp) * co(:,:,:) / tau_c
    ! CDOM absorption(cross section) at 320 nm 
    ! (Retamal-2007)--> taken number with low salinity
    ! and normalized to DOM concentration. 
    ! mgC/L = 1E-3 kgC/m^3. Hence 2.8 mgC/L=(2.8*1E-3)/(MW(C)) kmol/m^3
    ! 1 mol C =12 gC :  MW(C) = 12
    ! cdom_abs(i,j,k) = MAX(0.0_dp, 11.6_dp - 0.2_dp * sao(i,j,k))            &
    !      &          / (2.8_dp * 1.E-3_dp / 12._dp)
    !                                                             ! [m^2/kmolC]

    ! mz_bk_20120721+
    ! DO j=1,je
    !    DO i=1,ie
    !       DO k=1,ke
    !          ! cdom_abs(i,j,k) = MAX(0.0_dp, 11.6_dp - 0.2_dp * sao(i,j,k))   &
    !          !      &          / (2.8E-3_dp / 12._dp)
    !          !                                        ! [m^2/kmolC] <-- Retamal
    !          cdom_abs(i,j,k) = MAX(0.0_dp, 0.87_dp) / (1.15E-3_dp / 12._dp)
    !          !                                        ! [m^2/kmolC] <-- Retamal
    !       ENDDO
    !    ENDDO
    ! ENDDO

    ! DO j=1,je
    !    DO i=1,ie
    !       DO k=1,ke
    !          !                      ! radiation [photons/m^2/s] => [kmol/m^2/s]
    !          phoprod(i,j,k) = (strahl(i,j) * abs_oce(i,j,k) / 0.0642E-17_dp   &
    !               &            / 6.023E23_dp / 1000._dp)                      &
    !               &         * (ocetra(idoc)%ptr(i,j,k) * 12._dp / 30._dp      &
    !               !         ! DOM cross section [kmolC/m^3]*[m^2/kmolC] = [1/m]
    !               &         * cdom_abs(i,j,k))                                &
    !               &         *   4.7E-5               ! apparent quantum yield :
    !          !         ! linear interpolation between 325 and 313 nm wavelenght
    !       ENDDO
    !    ENDDO
    ! ENDDO

    DO k=1,ke
       DO j=1,je
          DO i=1,ie
             ! cdom_abs(i,j,k) = MAX(0.0_dp, 11.6_dp - 0.2_dp * sao(i,j,k))   &
             !      &          / (2.8E-3_dp / 12._dp)
             !                                        ! [m^2/kmolC] <-- Retamal
             cdom_abs(i,j,k) = 0.87_dp / (1.15E-3_dp / 12._dp)
             !                                        ! [m^2/kmolC] <-- Retamal
             !                      ! radiation [photons/m^2/s] => [kmol/m^2/s]
             phoprod(i,j,k) = (strahl(i,j) * abs_oce(i,j,k) / 0.0642E-17_dp   &
                  &            / 6.023E23_dp / 1000._dp)                      &
                  &         * (ocetra(idoc)%ptr(i,j,k) * 12._dp / 30._dp      &
                  !         ! DOM cross section [kmolC/m^3]*[m^2/kmolC] = [1/m]
                  &         * cdom_abs(i,j,k))                                &
                  &         *   4.7E-5               ! apparent quantum yield :
             !         ! linear interpolation between 325 and 313 nm wavelenght
          ENDDO
       ENDDO
    ENDDO
    ! mz_bk_20120721-

    co(:,:,:) = co(:,:,:) + (phoprod(:,:,:) + consumption(:,:,:) + dark) * dtbgc

    !--------------CHLOROPHYLL--------------------------------------------!
    !     See Moore et al., Deep-Sea Research II 49 (2002), 403-462       !
    !---------------------------------------------------------------------!
    ! Analogous to Moore et al., Deep-Sea Research II 49 (2002), 403-462
    ! 1 kmolP = (122*12/60)*10^6 mg[Chlorophyll] 
    !  hence [Chl] mg/m^3 = rcar*(12./ctochl)*1.e6 [Phy] kmolP/m^3 
    !  hence [Chl] mg/L   = rcar*(12./ctochl)*1.e3 [Phy] kmolP/m^3

    chl(:,:,:) = (rcar*(12._dp/ctochl) * 1.E3_dp * (ocetra(iphy)%ptr(:,:,:)))

    !--------------ISOPRENE---------------!
    !     See Broadgate et al (1997)      !
    !-------------------------------------!
    !   [Chl] ug/L = [Chl]*1E3 mg/L  
    !   Broadgate-1997 => [ISOP] [pmol/L] = 6.43*[Chl][ug/L]+1.2 =>
    !   [ISOP] [mol/L] = [ISOP] [kmol/m^3] = [ISOP][pmol/L]*1E-12 =>
    !   [ISOP] [kmol/m^3] = (6.43*([phy]*1/80*893.49*1E6)+1.2)*1E-12

    isop(:,:,:) = (6.43_dp * (chl(:,:,:) * 1.E3_dp) + 1.2_dp) * 1.E-12_dp

    !--------Methyliodide (CH3I)----------!
    !     See Bell et al (2002) JGR       !
    !-------------------------------------!

    ! mz_bk_20120721+
    ! DO j=1,je
    !    DO i=1,ie
    !       DO k=1,ke
    !          !                      ! radiation [photons/m^2/s] => [kmol/m^2/s]
    !          phoprod(i,j,k) = (strahl(i,j) * abs_oce(i,j,k) / 0.0642E-17_dp   &
    !               &            / 6.023E23_dp / 1000._dp)                      &
    !               &         * (ocetra(idoc)%ptr(i,j,k) * 12._dp / 30._dp      &
    !               !         ! DOM cross section [kmolC/m^3]*[m^2/kmolC] = [1/m]
    !               &         * cdom_abs(i,j,k))                                &
    !               &         * 0.1                      ! apparent quantum yield
    !       ENDDO
    !    ENDDO
    ! ENDDO
    DO k=1,ke
       DO j=1,je
          DO i=1,ie
             !                      ! radiation [photons/m^2/s] => [kmol/m^2/s]
             phoprod(i,j,k) = (strahl(i,j) * abs_oce(i,j,k) / 0.0642E-17_dp   &
                  &            / 6.023E23_dp / 1000._dp)                      &
                  &         * (ocetra(idoc)%ptr(i,j,k) * 12._dp / 30._dp      &
                  !         ! DOM cross section [kmolC/m^3]*[m^2/kmolC] = [1/m]
                  &         * cdom_abs(i,j,k))                                &
                  &         * 0.1                      ! apparent quantum yield
          ENDDO
       ENDDO
    ENDDO
    ! mz_bk_20120721-

    ! The United Nations Scientific, Education and Cultural Organization
    ! (UNESCO) definitions (1969) :
    ! S (ppt) = 1.80655 Cl (ppt) (1969)
    ! Cl-(mg/L) concentration ~= salinity(psu)* 553.5412
    ! Cl-  +  CH3I -->  
    ! reaction rates [1/M s]
    cl_react(:,:,:) = 7.78_dp * 1.E13_dp * exp(-13518._dp / tho(:,:,:))
    ! [g/L]/[MW] -->   [mol/L] = [M]
    cl_conc(:,:,:) = sao(:,:,:) * 0.55354_dp / 35.45_dp
    cl_react(:,:,:) = cl_react(:,:,:) * cl_conc(:,:,:) * ch3i(:,:,:)    ! [M/s]
    ch3i(:,:,:) = (phoprod(:,:,:) - cl_react(:,:,:)) * dtbgc

  END SUBROUTINE hamocc_photoprod
  !----------------------------------------------------------------------------

  !---------------------------------------------------------------------------- 
  SUBROUTINE hamocc_read_nml_ctrl(status, iou)
    !
    ! MPIOM MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Pozzer Andrea, MPICH, Oct 2004


    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status      ! error status
    INTEGER, INTENT(IN)  :: iou         ! I/O unit

    ! mz_bk_20120606 added L_BGC_DIAG for diagnostic bio-geochemistry output
    NAMELIST /CTRL/ L_CHECK, L_AGG, L_BGC_DIAG, isac, rmasks                  &
         &        , hamocc_ini_files_path                                     &
    ! mz_bk_20120724+
         ! make all biogeochemistry model parameters available for user control
         ! via namelist
         &        , ro2ut, rcar, rnit, nitdem, n2prod, rcalc, ropal     &
         &        , n2_fixation, rno3, perc_diron, riron, fesoly, relaxfe     &
         &        , pi_alpha, fPAR                                            &
#if defined MPIOM_13B
         &        , ctochl, atten_w, atten_f                                  &
#endif
         &        , phytomi, bkphy, bkopal, remido, dyphy, gammap             &
         &        , bkzoo, grami, zinges, epsher, grazra, spemor, gammaz      &
         &        , ecan                                                      &
         &        , sinkspeed_poc, sinkspeed_opal, sinkspeed_cal              &
         &        , drempoc, dremdoc, dremn2o, denitrification                &
         &        , sulfate_reduction, dremopal, dremcalc, dphymor, dzoomor   &
         &        , dmspar, calmax
    ! mz_bk_20120724-
    !
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='hamocc_read_nml_ctrl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    !----------------------------------------------------------------------
    ! DEFAULT PARAMETER SETTINGS - CAN BE OVERWRITEN BY THE NAMELIST
    L_CHECK = .FALSE.
    L_AGG   = .FALSE.
    ! mz_bk_20120606+
    L_BGC_DIAG = .FALSE.
    ! mz_bk_20120606-
    isac    = 1
    rmasks  = 0._dp
    !----------------------------------------------------------------------

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)

    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! NO ERROR

  END SUBROUTINE hamocc_read_nml_ctrl
  !----------------------------------------------------------------------------

END MODULE messy_hamocc
