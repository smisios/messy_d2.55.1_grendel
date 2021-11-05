! ***********************************************************************
MODULE messy_e4chem
  ! ***********************************************************************
  ! Fast chemistry module based on ECHAM4/CHEM                       *
  ! Original code by B. Steil                                        *
  ! Converted to MESSy submodel by A. Baumgaertner, 2009             *
  !                                                                  *
  ! * References:                                                    *
  ! * -----------                                                    *
  ! * Steil,B.: PhD thesis, Fachbereich Geowissenschaften der        *
  ! *  Universitaet Hamburg, Hamburg, Germany, jan, 1999.       *
  ! * Steil et al: Development of a chemistry module for GCMs:       *
  ! *  first results of a multi-annual integration,                  *
  ! *  Ann. Geophysicae, 1998,                             *
  ! *  also available as Report No. 74, DLR - Institut fuer Physik   *
  ! *  der Atmosphaere, Report No. 74, Oberpfaffenhofen, Germany,    *
  ! *  1997.                                                         *
  ! ***********************************************************************

  USE messy_main_constants_mem, ONLY: DP
  USE messy_main_blather,       ONLY: warning ! mz_bk_20101216

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'e4chem'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.0'

  ! CTRL NAMELIST PARAMETER
  LOGICAL, PUBLIC :: l_fastscav = .FALSE.
  LOGICAL, PUBLIC :: l_Brparam = .FALSE. ! mz_ab_20100505

  INTEGER, PUBLIC :: NSTCHPH = 1 

  INTEGER, PARAMETER,PUBLIC :: NSPEC = 39

  INTEGER, PARAMETER,PUBLIC :: ind_H          = 1
  INTEGER, PARAMETER,PUBLIC :: ind_OH         = 2
  INTEGER, PARAMETER,PUBLIC :: ind_HO2        = 3
  INTEGER, PARAMETER,PUBLIC :: ind_N          = 4
  INTEGER, PARAMETER,PUBLIC :: ind_NO         = 5
  INTEGER, PARAMETER,PUBLIC :: ind_NO2        = 6
  INTEGER, PARAMETER,PUBLIC :: ind_NO3        = 7
  INTEGER, PARAMETER,PUBLIC :: ind_N2O5       = 8
  INTEGER, PARAMETER,PUBLIC :: ind_HNO4       = 9
  INTEGER, PARAMETER,PUBLIC :: ind_CL         = 10
  INTEGER, PARAMETER,PUBLIC :: ind_CLO        = 11
  INTEGER, PARAMETER,PUBLIC :: ind_HOCl       = 12 ! CLOH  --> HOCl
  INTEGER, PARAMETER,PUBLIC :: ind_CL2O2      = 13
  INTEGER, PARAMETER,PUBLIC :: ind_CL2        = 14
  INTEGER, PARAMETER,PUBLIC :: ind_HCHO       = 15 ! CH2O --> HCHO
  INTEGER, PARAMETER,PUBLIC :: ind_CH3O2      = 16
  INTEGER, PARAMETER,PUBLIC :: ind_CH4        = 17
  INTEGER, PARAMETER,PUBLIC :: ind_N2O        = 18
  INTEGER, PARAMETER,PUBLIC :: ind_H2O2       = 19
  INTEGER, PARAMETER,PUBLIC :: ind_HCl        = 20
  INTEGER, PARAMETER,PUBLIC :: ind_CO         = 21
  INTEGER, PARAMETER,PUBLIC :: ind_CH3OOH     = 22 ! CH3O2H --> CH3OOH
  INTEGER, PARAMETER,PUBLIC :: ind_ClNO3      = 23
  INTEGER, PARAMETER,PUBLIC :: ind_CFCl3      = 24
  INTEGER, PARAMETER,PUBLIC :: ind_CF2Cl2     = 25
  INTEGER, PARAMETER,PUBLIC :: ind_CH3CL      = 26
  INTEGER, PARAMETER,PUBLIC :: ind_CCL4       = 27
  INTEGER, PARAMETER,PUBLIC :: ind_CH3CCL3    = 28
  INTEGER, PARAMETER,PUBLIC :: ind_H2         = 29
  INTEGER, PARAMETER,PUBLIC :: ind_HNO3       = 30
  INTEGER, PARAMETER,PUBLIC :: ind_NAT        = 31
  INTEGER, PARAMETER,PUBLIC :: ind_OHAB       = 32
  INTEGER, PARAMETER,PUBLIC :: ind_HO2AB      = 33
  INTEGER, PARAMETER,PUBLIC :: ind_ICE        = 34
  INTEGER, PARAMETER,PUBLIC :: ind_H2O        = 35
  INTEGER, PARAMETER,PUBLIC :: ind_O3P        = 36
  INTEGER, PARAMETER,PUBLIC :: ind_O3         = 37
  INTEGER, PARAMETER,PUBLIC :: ind_O1D        = 38
  INTEGER, PARAMETER,PUBLIC :: ind_CO2        = 39


  INTEGER, PARAMETER,PUBLIC :: JVALS = 4, NUMRAT = 67, IRCTMIN = 150, IRCTMAX = 320
  INTEGER, PARAMETER,PUBLIC :: JDIFTE = IRCTMAX - IRCTMIN
  INTEGER, PARAMETER,PUBLIC :: NUMTEM = JDIFTE * JVALS + 1

  PUBLIC :: CHEMICS
  PUBLIC :: INRCGAS
  PUBLIC :: inisulnew
  PUBLIC :: CLSCAV
  PUBLIC :: e4chem_read_nml_ctrl

  ! Precalculated gas reaction rates
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, PUBLIC, SAVE :: RCGAS
  ! lookup table for sulphur 
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC, SAVE :: sulook 
  INTEGER, PARAMETER,  PUBLIC :: ksul=177 ! sulook length

  PUBLIC :: SPC_NAMES

  CHARACTER(LEN=32), PARAMETER, DIMENSION(39) :: SPC_NAMES =  (/ &
       'H                               ','OH                              ','HO2                             ',&
       'N                               ','NO                              ','NO2                             ',&
       'NO3                             ','N2O5                            ','HNO4                            ',&
       'CL                              ','CLO                             ','HOCl                            ',&
       'CL2O2                           ','CL2                             ','HCHO                            ',&
       'CH3O2                           ','CH4                             ','N2O                             ',&
       'H2O2                            ','HCl                             ',                                   &
       'CO                              ','CH3OOH                          ','ClNO3                           ',&
       'CFCl3                           ','CF2Cl2                          ','CH3CL                           ',&
       'CCL4                            ','CH3CCL3                         ','H2                              ',&
       'HNO3                            ','NAT                             ','OHAB                            ',&
       'HO2AB                           ','ICE                             ','H2O                             ',&
       'O3P                             ','O3                              ','O1D                             ',&
       'CO2                             '/) 

CONTAINS

  !
  !     ZTEMP   : TEMPERATURE IN [K]
  !     ZCON    : CONCENTRATION IN [MOLEC/CM3] OF THE AIR
  !     ZPRES   : PRESSURE IN [PA]
  !     ZTMSTDT : TIMESTEP [SEC]
  !     ZTMSTHDT: HALVE TIMESTEP [SEC]
  !
  !     ZINVCON : inverse of concentration
  !
  !     THE PREFIX ZMOLEC  MEANS [MOLECULES/CM3]
  !     THE PREFIX ZC      MEANS REACTIONCONSTANT
  !     THE PREFIX ZD      FOLLOWED BY THE NAME OF THE SPECIES,
  !                        MEANS PHOTODISSOCIATION FREQUENCY [SEC-1]
  !     THE PREFIX ZSTO    MEANS TEMPORARY STORAGE
  !     THE EXTENSION T0   MEANS VALUES FOR THE MOLECULE
  !                        CONCENTRATION AT T0. THESE ARE THE P1-VALUES
  !                        OF ECHAM WIHOUT CHEMISTRY !
  !     THE EXTENSION M    MEANS UNIT MOLECULES/CM3
  !     THE EXTENSION V    MEANS Vector (?) with index JL=1:KLON, corresponding
  !                        variable without V (scalar)
  !
  !
  !
  ! INPUT
  ! ===================
  ! KLON                : Number of longitudes  --> kproma
  ! KLEV                : Number of levels
  ! JROW                : Current Row
  ! PTEMP               : temperature
  ! PQ                  : Specific humidity
  ! Conc                : Tracer concentrations
  ! PAPP1               : PRESSURE AT FULL LEVELS (T+DT)
  ! PAPHP1              : PRESSURE AT HALF LEVELS (T+DT)
  ! DANI/M              : Day/night flag
  ! RJ_...              : Photolysis Rates
  ! zlat                : latitude [deg]
  ! imonth              : current month
  ! dtime               : time_step_length
  ! tp_i                : tropopause level index

      SUBROUTINE CHEMICS(                                             &
       KLON,      KLEV,     JROW,                                     &
       PTEMP,     PQ,       Conc,                                     &
       PAPP1,     PAPHP1,   DANI,     DANIM,                          &
       RJ_O3P,    RJ_O1D,   RJ_NO2,   RJ_HNO3,                        &
       RJ_COH2,   RJ_CHOH,  RJ_N2O5,  RJ_HNO4,                        &
       RJ_NO2O,   RJ_NOO2,  RJ_H2O2,  RJ_CH3OOH,                      &
       RJ_O2,     RJ_CFC11, RJ_CFC12, RJ_N2O,                         &
       RJ_CLONO2, RJ_CL2O2, RJ_HOCL,                                  &
       RJ_CCL4,   RJ_CH3CL, RJ_CH3CCL3, RJ_HCL,                       &
       RJ_H2O,    RJ_NO,    RJ_CO2,                                   &
       zlat,      IMONTH,   DTIME, tp_i,                              &
       ZDELTAO3_BRV,                                                  &
       ZPRODO2,   ZPRODCO,                                            &
       ZPRODCH4,  ZDESTH12,                                           &
       ZDESTH14,  ZDESTN13,                                           &
       ZDESTC1,   ZDESTCL2O2,                                         &
       ZDESTCLOH, ZDESTH8,                                            &
       ZDESTH4,   ZDESTH11,                                           &
       ZDESTO1)

    USE messy_main_constants_mem, ONLY : dp, &
         zavogadro => N_A,                   & ! Avogadro constant / (1/mol)
         zm_air => M_air,                    & ! molar mass of dry air [g/mol]
         MO, MC,                             &
         vtmpc1,                             &
         rd



    IMPLICIT NONE

    INTEGER,                     INTENT(IN) :: KLON, KLEV, JROW
    REAL(DP), DIMENSION(:,:),    INTENT(IN) :: PTEMP,  PQ
    REAL(DP), DIMENSION(:,:,:),  INTENT(INOUT) :: Conc

    REAL(DP), DIMENSION(:,:),  INTENT(IN)  :: RJ_O3P,                   &
         RJ_O1D,   RJ_NO2,   RJ_HNO3,                                   &
         RJ_COH2,   RJ_CHOH,  RJ_N2O5,  RJ_HNO4,                        &
         RJ_NO2O,   RJ_NOO2,  RJ_H2O2,  RJ_CH3OOH,                      &
         RJ_O2,     RJ_CFC11, RJ_CFC12, RJ_N2O,                         &
         RJ_CLONO2, RJ_CL2O2, RJ_HOCL,                                  &
         RJ_CCL4,   RJ_CH3CL, RJ_CH3CCL3, RJ_HCL,                       &
         RJ_H2O,    RJ_NO,    RJ_CO2


    REAL(dp), DIMENSION(:,:), INTENT(IN) :: PAPP1, PAPHP1 
    REAL(dp), DIMENSION(:),   INTENT(IN) :: DANI, DANIM 
    REAL(dp), DIMENSION(:),   INTENT(IN) :: zlat

    INTEGER, INTENT(IN)  :: IMONTH
    REAL(dp), INTENT(IN) :: DTIME

    REAL(dp), DIMENSION(:), INTENT(IN) :: tp_i ! tropopause level index

    ! mz_bk_20101222+ ozone diagnostics/budget
    REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: zdeltao3_brv
    REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: ZPRODO2,   ZPRODCO      &
         , ZPRODCH4,   ZDESTH12,   ZDESTH14,  ZDESTN13,    ZDESTC1     &
         , ZDESTCL2O2, ZDESTCLOH,  ZDESTH8,   ZDESTH4,     ZDESTH11    &
         , ZDESTO1
    ! mz_bk_20101222-

    REAL(dp), PARAMETER :: RMULTRHO=zavogadro*0.001/zm_air, ZQMIN=0._dp

    !     PARAMETERS FOR ICE-PARTICLES
    !
    REAL(dp), PARAMETER :: RHOICE =0.92_dp
    !      PARAMETER(ZNUMICE=0.03536,RHOICE =0.92,
    !     * ZHELPICE=18./(6.0E+23*RHOICE*ZNUMICE))
    !
    !     TEMPERATURE AND PRESSURE INDEPENDENT RATE CONSTANT
    !
    REAL(dp), PARAMETER :: ZCO2   = 1.2E-10_dp, &
         ZCN3  = 5.3E-11_dp, &
         ZCN4  = 2.0E-10_dp, &
        ! ZCN21 = 6.7E-11_dp, & ! mz_ab_20091221
        ! ZCN22 = 4.9E-11_dp, & ! mz_ab_20091221
         ZCC12 = 1.65E-10_dp, &
         ZCC13 = 3.0E-10_dp, &
         ZCC16 = 1.4E-10_dp, &
         ZCC17 = 2.3E-10_dp, &
         ZCC18 = 3.3E-10_dp, &
         ZCH3  = 1.1E-10_dp, &
         ZCH7  = 1.5E-10_dp, &
        ! ZCH8  = 2.2E-10_dp, & ! mz_ab_20091221
         ZCH10 = 7.2E-11_dp, &
        ! ZCH16 = 1.0E-11_dp,    & ! mz_ab_20091221
         ZCH17 = 6.9E-12_dp,  &
         ZCH18 = 1.6E-12_dp

    REAL(dp), DIMENSION(:),   ALLOCATABLE :: ZTEMPV,ZPRESV,ZCONV,ZINVCONV ! (KLON)

    REAL(dp), DIMENSION(:),   ALLOCATABLE :: ZMOLECO3V, ZMOLECOV, ZMOLECO1DV & !(KLON)
         , ZMOLECOXV, ZMOLECHV, ZMOLECOHV, ZMOLECHO2V &
         , ZMOLECNV, ZMOLECNOV, ZMOLECNO2V            &
         , ZMOLECNO3V, ZMOLECHNO4V, ZMOLECN2O5V       &
         , ZMOLECCLV, ZMOLECCLOV, ZMOLECCLOHV         &
         , ZMOLECCL2O2V, ZMOLECCLNO3V, ZMOLECCL2V     &
         , ZMOLECCH2OV, ZMOLECCH3O2V,  ZMOLECCH3OV    &
         , ZMOLECCFCL3V, ZMOLECCF2CL2V,  ZMOLECCCL4V  &
         , ZMOLECCH3CLV, ZMOLECCH3CCL3V, ZMOLECH2V    &
         , ZMOLECCH4V, ZMOLECN2OV, ZMOLECH2O2V        &
         , ZMOLECHCLV, ZMOLECCOV, ZMOLECCH3O2HV       &
         , ZMOLECNATV, ZMOLECICEV, ZMOLECHNO3V        &
         , ZMOLECCLOXV, ZMOLECNOXV, ZMOLECODDV        &
         , ZFMOLECNOXV, ZFMOLECHNO3V, ZFMOLECHCLV     &
         , ZMOLECH2OV, ZMOLECNATT0V, ZMOLECCLOXT0V    &
         , ZMOLECH2OT0V, ZFMOLECCLOXV

    REAL(dp), DIMENSION(:), ALLOCATABLE :: ZMOLECO2V,ZMOLECN2V,ZMOLECCO2V & ! (KLON)
         , ZSCALEO3

    REAL(dp), DIMENSION(:), ALLOCATABLE :: ZCP1V,ZCP2V,ZCP3V,ZCP4V &
         , ZCS1V, ZCS2V, ZHET1V, ZHET2V, ZHET3V ! (KLON)
    REAL(dp), DIMENSION(:), ALLOCATABLE :: ZUDTAICE,ZUDTBICE &
         , ZUDTANAT,ZUDTBNAT, ZHEILEV ! (KLON)

    INTEGER, DIMENSION(:), ALLOCATABLE :: INIHET0V, JINDV ! (KLON)

    REAL(dp), DIMENSION(:), ALLOCATABLE :: ZNUMICE,ZHELPICE !(KLEV)
     
    REAL(dp) :: ZQOHOUFA, ZRATOHHO2, ZRATHO2HOX , ZQHO2OUFA, ZQHOX, ZEQHOX, ZLIMHOX            &
         , ZMOLECHOX, ZAUX1CLOX, ZEMEQCLNO3 , ZAUX2CLOX, ZAUX3CLOX, ZAUX4CLOX, ZAUX10CLOX      &
         , ZAUX5CLOX, ZAUX6CLOX, ZAUX7CLOX, ZWCLOX, ZUCLOX, ZVCLOX, Z1QHNO4, Z1QNO3            &
         , Z1QN2O5, ZPNO3 , ZEQNO3, ZAUX10NOX, ZQHNO3, ZEQHNO3, ZLIMHNO3, ZEQNOX, ZLIMNOX      &
         , ZP2NOX, ZQ1NOX, ZQNOX, ZSCALATOMSN, ZSTOSCALTOT, ZSTOSCALNOX, ZQCH4, ZSTODENOMCH3O2 &
         , ZSTONUMERCH3O2, ZSTOPRODCH22HO2, ZSTOPRODCH23O2, ZQCH3O2H, ZPCH3O2H, ZQCH2O, ZPCH2O &
         , ZQH2, ZPH2, ZSTO0, ZSTO1, ZSTO2, ZSTO3, ZSTO4, ZSTO5, ZDESTODD, ZQH2O     &
         , ZPROH2O, ZMOLECHOXT0, ZQOHINFA                                 &
         , ZSTOEHOX, ZSTOFHOX, ZSTOKKHOX, ZSTOCCHOX, ZSTOFFHOX, ZAUXHOX, ZSTOPHOX, ZSTOQHOX    &
         , ZQH2O2, ZHETLOSS1, ZHETLOSS2, ZHETLOSSNOXTEST, ZHETLOSSNOX, SCALE, ZHETLOSSHCL      &
         , ZP1HCL, ZP2HCL, ZQ1HCL, ZQCLOX, ZEQCLOX, ZAUX, ZSCY, ZSTOSCALCLOX, ZQN2O, ZMOLECNOYT0 &
         , ZQ1NOY, ZQNOY, ZMOLECNOY, ZP1HNO3, ZP1NOX                                           &
         , ZC1,  ZC2,  ZC3,  ZMAX,  ZMAX1,  ZAUXNNEW1,  ZAUXNNEW2,  ZPCLNO3,  ZPCLOH           &
         , ZPN2O5,  ZSTO1DENOMH,  ZSTO2DENOMH,  ZSTODENOMH,  ZSTOKH,  ZSTOKKH,  ZSTOBH         &
         , ZSTOBBH,  ZSTOK1OH,  ZSTOKOH,  ZSTOB1OH,  ZSTOBOH,  ZSTOCOH,  ZSTOEOH,  ZSTODENOMOH &
         , ZSTOKKOH,  ZSTOCCOH,  ZSTOKHOX,  ZSTOBHOX,  ZSTOCHOX,  ZSTODHOX                     &
         , ZU, ZV, ZW, ZQNO3A, ZQNO3, ZQN2O5, ZEQN2O5, ZAUX1NOX, ZEMEQN2O5            &
         , ZQHNO4, ZEQHNO4, ZAUX2NOX, ZEMEQHNO4, ZAUXK1, ZAUXK2, ZAUXK3, ZAUXK4                &
         , ZKNOX, ZG1, ZG2, ZAUXF, ZF1, ZF2, ZF3, ZANOX, ZBNOX, ZCNOX, ZDNOX, ZC0              &
         , ZDNO, ZDCO2, ZDHCL, ZDCLOH, ZDCH3O2H, ZDHNO4, ZDNO3B, ZDCL2O2A                      &
         , ZPRCFCL3, ZPRCF2CL2, ZPRCCL4, ZPRCH3CL                      &
         , ZPRCH3CCL3, ZPRODCLOY, ZMOLECCLOY, ZRATNONO2, ZRATNNO, ZDECL                        &
         , ZQCLOH, ZEQCLOH, ZAUXU1, ZQCLNO3, ZEQCLNO3, ZAUXU4, ZAUXU2                          &
         , ZAUXU3, ZAUXU10, ZAUXU11                                                            &
         , ZCN18, ZCH13, ZCCLOCLO, ZCN5, ZCN7, ZCEQ, ZCCL2O2M, ZK0                             &
         , ZK3, ZCN19, ZCH1, ZCH2, ZCH5, ZCH8, ZCN2, ZCN23, ZCH20, ZCH21                       &
         , ZDO2, ZDH2O2, ZDH2O, ZDHNO3, ZDN2O, ZDN2O5, ZDCH2OA                                 &
         , ZDCH2OB, ZDCLNO3, ZDCH3CL, ZDCCL4, ZDCFCL3, ZDCF2CL2, ZDCH3CCL3 &
         , ZCN20, ZCC1, ZCC2, ZCC3, ZCC3B, ZCC4, ZCC5, ZCC6, ZCC7, ZCC8 &
         , ZCC9, ZCC10, ZCC11, ZCC14, ZCC15, ZCH4, ZCH6, ZCH9, ZCH11 &
         , ZCH12, ZCH14, ZCH15, ZCH19, ZCH22, ZCH23, ZK1, ZK2, ZCN14 &
         , ZCN16, ZCN17 
    REAL(dp) ::           ZMOLECN2O5T0, ZMOLECHNO4T0, ZMOLECNO2T0 &
         , ZMOLECCH3O2T0 & 
         , ZMOLECHO2T0 &
         , ZCO1, ZCN1, ZCN6, ZCN8, ZCN9 &
         , ZCN10, ZCN11, ZCN12, ZCN13, ZCN15, ZCN21, ZCN22, ZCH16 &
         , ZMOLECOXT0, ZDO3B, ZDO3A, ZDNO3A, ZDNO2, ZCO4, ZCO5 &
         , ZCO3, ZSTO, ZRATO1DO3, ZRATOO3, ZRATO3OX, ZRATOOX &
         , ZRATO1DOX, ZHET1, ZHET2, ZHET3, ZMOLECO2 &
         , zmolecco2, ZMOLECHCLT0, ZMOLECCLOXT0 &
         , ZMOLECCLNO3T0, ZMOLECCLOHT0, ZMOLECCLOT0 &
         , ZMOLECCL2O2T0, ZMOLECNOXT0 &
         , VHNO3, sulup, SUL, SULSURF, VOL0, WEIGHT &
         , RHO, ZHEHCL1, ZHEHCL2, ZHEHCL, ZALPH, ZLOW &
         , ZHIGH, ZNEU, ZEQUA, ZFMOLECHCL, ZFMOLECCLOX, ZFMOLECNOX &
         , ZFMOLECHNO3, ZTEIND, zjind &
         , ZCPN4, ZNATSURF, ZTEMPNUCLI, ZMOLECH2OBFC, ZMOLECTOTH2O &
         , ZMOLECEQH2O, ZDIFFH2O, ZVH2O, ZRADICE, ZICEDENOM1, ZICEDENOM2 &
         , ZTEMPNUCLN, ZMOLECTOTHNO3, ZMOLECEQHNO3, ZDIFFHNO3, ZVHNO3 &
         , ZRADNAT, ZNATDENOM1, ZNATDENOM2, ZMOLECODD
    REAL(dp) ::  ZSULSURF, ZVOL0, ZSCALSURF, ZPH2O, ZWEIGHT, ZCSS, ZRHO &
         , ZGAMMAS1, ZMOLECCF2CL2, ZMOLECCH3CL, ZMOLECCCL4, ZMOLECCH3CCL3 &
         , ZMOLECCLOX, ZMOLECNOX, ZMOLECOX, ZMOLECH2O, ZCP1, ZCP2, ZCP3 &
         , ZCP4, ZCS1, ZCS2, ZCPI1, ZCPI2, ZCPI3, ZCPI4, ZICESURF &
         , ZCPN1, ZCPN2, ZCPN3 &
         , ZMOLECN, ZMOLECNO, ZMOLECNO2, ZMOLECNO3, ZMOLECN2O5, ZMOLECHNO4 &
         , ZMOLECCL, ZMOLECCLO, ZMOLECCLOH, ZMOLECCL2O2, ZMOLECCL2, ZMOLECCH2O &
         , ZMOLECCH3O2, ZMOLECCH3O, ZMOLECCH4 &
         , ZMOLECH2, ZMOLECN2O, ZMOLECH2O2, ZMOLECHCL, ZMOLECCO, ZMOLECCH3O2H &
         , ZMOLECICE, ZMOLECCLNO3, ZMOLECHNO3, ZMOLECNAT &
         , ZWEIGHT0, ZRHO0, ZUDT, ZTEMP &
         , ZQ, ZTV, ZPRES, RHOAIR, ZCON, ZINVCON, ZMOLECO3 &
         , ZMOLECO, ZMOLECO1D, ZMOLECH, ZMOLECOH, ZMOLECHO2 &
         , ZMOLECCFCL3, TWODT, ZTMSTDT, ZTMSTHDT &
         , ZCHLOR_SCALE

    REAL(dp), DIMENSION(:), ALLOCATABLE :: ZRNATUP, ZRICEUP

    REAL(dp), PARAMETER :: RNATUP=2000._dp, RICEUP=RNATUP, &
         rsulup=500._dp,rsuldo=18000._dp


    INTEGER, DIMENSION(:), ALLOCATABLE :: IIHET, INHET, INSED
    INTEGER :: JL, JK,  ISHET, JIND, JSUL, INIHET0


    ! zu_as_20100504+ Variables for Br parameterisation
    REAL(dp) :: zbrtemp     ! temperature [K]
    REAL(dp) :: zbrpres     ! pressure [hPa]
    REAL(dp) :: zproxy      ! scaling factor
!!$    REAL(dp) :: zfac_br     ! another factor ! op_pj_20101220
    REAL(dp) :: zdeltao3_br = 0._dp ! o3 destroyed by bromine [molec*cm-3]
!!$    REAL(dp), ALLOCATABLE :: zdeltao3_brv(:) ! mz_bk_20101222
    REAL(dp) :: zclox       ! volume mixing ratio (ClO + 2*Cl2O2) in ppbv
    REAL(dp), PARAMETER :: zbr1=-0.064996_dp, zbr2=0.8939601_dp &
         , zbr3=-3.760605_dp, zbr4=0.9998547_dp, zbr5=-3.674095_dp &
         , zbr6=1.0190500_dp
    ! zu_as_20100504- 

    INTRINSIC ABS, TINY, MAX, MIN, INT, EXP, LOG, REAL, SQRT

    ! mz_ab_20100504+
!!$    ALLOCATE(zdeltao3_brv(KLON)) ! mz_bk_20101222
    zdeltao3_brv = 0.0_dp ! mz_bk_20101213
    ! mz_ab_20100504-
    ! mz_bk_20101222+
    ZPRODO2    = 0.0_dp
    ZPRODCO    = 0.0_dp
    ZPRODCH4   = 0.0_dp
    ZDESTH12   = 0.0_dp
    ZDESTH14   = 0.0_dp
    ZDESTN13   = 0.0_dp
    ZDESTC1    = 0.0_dp
    ZDESTCL2O2 = 0.0_dp
    ZDESTCLOH  = 0.0_dp
    ZDESTH8    = 0.0_dp
    ZDESTH4    = 0.0_dp
    ZDESTH11   = 0.0_dp
    ZDESTO1    = 0.0_dp
    ! mz_bk_20101222-

    ALLOCATE(ZTEMPV (KLON),ZPRESV (KLON),ZCONV (KLON),ZINVCONV (KLON))

    ALLOCATE(ZMOLECO3V(KLON),ZMOLECOV(KLON), ZMOLECO1DV(KLON)                 & 
         , ZMOLECOXV(KLON), ZMOLECHV(KLON), ZMOLECOHV(KLON), ZMOLECHO2V(KLON) &
         , ZMOLECNV(KLON), ZMOLECNOV(KLON), ZMOLECNO2V(KLON)                  &
         , ZMOLECNO3V(KLON), ZMOLECHNO4V(KLON), ZMOLECN2O5V(KLON)             &
         , ZMOLECCLV(KLON), ZMOLECCLOV(KLON), ZMOLECCLOHV(KLON)               &
         , ZMOLECCL2O2V(KLON), ZMOLECCLNO3V(KLON), ZMOLECCL2V(KLON)           &
         , ZMOLECCH2OV(KLON), ZMOLECCH3O2V(KLON), ZMOLECCH3OV(KLON)           &
         , ZMOLECCFCL3V(KLON), ZMOLECCF2CL2V(KLON), ZMOLECCCL4V(KLON)         &
         , ZMOLECCH3CLV(KLON), ZMOLECCH3CCL3V(KLON), ZMOLECH2V(KLON)          &
         , ZMOLECCH4V(KLON), ZMOLECN2OV(KLON), ZMOLECH2O2V(KLON)              &
         , ZMOLECHCLV(KLON), ZMOLECCOV(KLON), ZMOLECCH3O2HV(KLON)             &
         , ZMOLECNATV(KLON), ZMOLECICEV(KLON), ZMOLECHNO3V(KLON)              &
         , ZMOLECCLOXV(KLON), ZMOLECNOXV(KLON), ZMOLECODDV(KLON)              &
         , ZFMOLECNOXV(KLON), ZFMOLECHNO3V(KLON), ZFMOLECHCLV(KLON)           &
         , ZFMOLECCLOXV(KLON), ZMOLECH2OV(KLON)                               &
         , ZMOLECH2OT0V(KLON), ZSCALEO3(KLON))

    ALLOCATE(ZMOLECO2V(KLON),ZMOLECN2V(KLON),ZMOLECCO2V(KLON),JINDV(KLON))

    ALLOCATE(ZCP1V(KLON),ZCP2V(KLON),ZCP3V(KLON),ZCP4V(KLON),           &
         ZCS1V(KLON),ZCS2V(KLON), ZHET1V(KLON),ZHET2V(KLON),ZHET3V(KLON))
    ALLOCATE(ZUDTAICE(KLON),ZUDTBICE(KLON), ZUDTANAT(KLON),ZUDTBNAT(KLON), &
         ZHEILEV(KLON))

    ALLOCATE(INIHET0V(KLON))

    ALLOCATE(ZNUMICE(KLEV),ZHELPICE(KLEV))

    ALLOCATE(IIHET(KLON),INHET(KLON),INSED(KLON))
    ALLOCATE(ZRNATUP(KLON),ZRICEUP(KLON))

    TWODT       = 2.*DTIME
    ZTMSTDT     = 0.5*TWODT*REAL(NSTCHPH)
    ZTMSTHDT    = 0.25*TWODT*REAL(NSTCHPH)
    !
    ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !
    DO JL=1,KLON
       ZUDTAICE(JL) = 0._dp
       ZUDTANAT(JL) = 0._dp
    ENDDO

    ! ******************************************************************
    !
    ! IIHET = 1 <==> het reactions on ice are permitted
    ! INHET = 1 <==> het reactions on NAT are permitted
    ! ISHET = 1 <==> het reactions on sulfate aerosols are permitted
    !
    ! Het. reactions on ICE and NAT are only permitted for
    ! latitudes poleward of 40deg.
    ! Latitudes are formulated resolution independent.
    ! R. Hein, 28-aug-1997
    !
    ! Vertical restrictions for het reactions are given inside
    ! the loop 1100 (over longitudes).
    ! They depend on tropopause heights and limits set in comdeck COMSUL2
    ! New format for reading sulfate aerosol surfaces has been included
    ! by B. Steil, 03-apr-1998
    ! R. Hein, 08-apr-1998
    !
    ! ******************************************************************
    DO JL=1,KLON
       IF (abs(zlat(JL)) > 40._dp) THEN
          IIHET(JL) = 1
          INHET(JL) = 1
       ELSE
          IIHET(JL) = 0
          INHET(JL) = 0
       ENDIF
    ENDDO

    ISHET     = 1
    ZWEIGHT0  = 0.75_dp
    ZRHO0     = 1.7275378_dp

    !!$    ! FOR DEBUGGING+
!!$        IIHET = 0
!!$        INHET = 0
!!$        ISHET = 0
    !!$    ! FOR DEBUGGING-

    ! sedimentation of nat, only in winter and only in latitudes
    ! poleward of 55 deg, because of poorly defined hygropause!
    !
    ! seasonal variation of rnatup and riceup!
    !
    ZRNATUP(:) = RNATUP
    ZRICEUP(:) = RICEUP
    !
    INSED(:) = 0
    DO JL=1,KLON
       IF ( (abs(zlat(JL)) > 55._dp).AND.  &
            (IMONTH.LE.10).AND.(IMONTH.GE.5) ) INSED(JL)=1
       IF ( (abs(zlat(JL)) >  55._dp).AND. &
            ((IMONTH.LE.4).OR.(IMONTH.EQ.12)) ) INSED(JL)=1
       !
       !
       IF ((zlat(JL)<0.).and.(IMONTH.GT.7)) THEN
          ZRNATUP(JL) = 2800._dp
          ZRICEUP(JL) = ZRNATUP(JL)
       ENDIF
    ENDDO



    ! ******************************************************************
    ! *                                                                *
    ! * start loop over levels                                         *
    ! *                                                                *
    ! ******************************************************************
    DO  JK=1,KLEV
       ! ******************************************************************
       ! *                                                                *
       ! * start loop over longitudes                                     *
       ! *                                                                *
       ! ******************************************************************
       DO JL=1,KLON
          !
          IF (PAPP1(JL,JK)<2744._dp) ZNUMICE(JK)=0.1_dp
          IF ((PAPP1(JL,JK)>=2744._dp).and.(PAPP1(JL,JK)<4709._dp)) ZNUMICE(JK)=0.054_dp
          IF (PAPP1(JL,JK)>=4709._dp) ZNUMICE(JK)=0.035_dp
          ZHELPICE(JK)=18._dp/(6.0E+23_dp*RHOICE*ZNUMICE(JK))

          ZUDTBICE(JL) = ZUDTAICE(JL)
          ZUDTBNAT(JL) = ZUDTANAT(JL)
          ZUDT         = 0._dp
          !
          ZTEMP      = PTEMP(JL,JK) !PTM1(JL,JK)+PTTE(JL,JK)*twodt ! mz_ab_20090617
          ZTEMPV(JL) = ZTEMP
          ZQ         = PQ(JL,JK) !PQM1(JL,JK)+PQTE(JL,JK)*twodt ! mz_ab_20090617
          !          ZQBFCOND   = PQM1(JL,JK)+PQTEBFCOND(JL,JK)*twodt ! ZQBFCOND mz_ab_20090617 ZQBFCOND not needed
          !MAR ---------------------------------------------------
          ZQ         = MAX(ZQ,ZQMIN)
          !          ZQBFCOND   = MAX(ZQBFCOND,ZQMIN) ! mz_ab_20090617 ZQBFCOND not needed
          !MAR ---------------------------------------------------
          ZTV     = ZTEMP*(1._dp+VTMPC1*ZQ)
          ZPRES   = PAPP1(JL,JK)
          ZPRESV(JL) = ZPRES
          RHOAIR  = ZPRES/(RD*ZTV)
          ZCON    = RHOAIR*RMULTRHO
          ZCONV(JL) = ZCON
          ZINVCON = 1._dp/ZCON
          ZINVCONV(JL) = ZINVCON

          ZHEILEV(JL)  = (PAPHP1(JL,JK+1)-PAPHP1(JL,JK))/ &
               (0.0980665_dp*RHOAIR)
          !
          !     SET LOCAL VARIABLES, UNITS [molec/cm3]
          !
          ZMOLECO3         = Conc(JL,JK,ind_O3)
          ZMOLECO3         = MAX(ZMOLECO3,ZINVCON)
          ZMOLECO          = Conc(JL,JK,ind_O3P)
          ZMOLECO          = MAX(ZMOLECO,ZINVCON)
          ZMOLECO1D        = Conc(JL,JK,ind_O1d)
          ZMOLECO1D        = MAX(ZMOLECO1D,ZINVCON)
          ZMOLECH          = Conc(JL,JK,ind_H) 
          ZMOLECH          = MAX(ZMOLECH,ZINVCON)
          !
          ! Morning: use saved OH and HO2
          IF ((DANI(JL).GT.1.).AND.(DANIM(JL).LT.1.)) THEN
             ZMOLECOH      = Conc(JL,JK,ind_OHAB)
             ZMOLECHO2     = Conc(JL,JK,ind_HO2AB)
          ELSE
             ZMOLECOH      = Conc(JL,JK,ind_OH)
             ZMOLECHO2     = Conc(JL,JK,ind_HO2)
          ENDIF
          ZMOLECOH         = MAX(ZMOLECOH,ZINVCON)
          ZMOLECHO2        = MAX(ZMOLECHO2,ZINVCON)
          !
          ! Evening: save OH and HO2
          IF ((DANI(JL).LT.1.).AND.(DANIM(JL).GT.1.)) THEN
             Conc(JL,JK,ind_OHAB)   = Conc(JL,JK,ind_OH) 
             Conc(JL,JK,ind_HO2AB)  = Conc(JL,JK,ind_HO2)
          ENDIF
          ! For Debugging+
!!$          ZMOLECOH      = Conc(JL,JK,ind_OH)
!!$          ZMOLECHO2     = Conc(JL,JK,ind_HO2)
!!$          ZMOLECOH         = MAX(ZMOLECOH,ZINVCON)
!!$          ZMOLECHO2        = MAX(ZMOLECHO2,ZINVCON)
          ! For Debugging-
          !
          ZMOLECN       = Conc(JL,JK,ind_N)
          ZMOLECN       = MAX(ZMOLECN,TINY(0._dp))
          ZMOLECNO      = Conc(JL,JK,ind_NO)
          ZMOLECNO      = MAX(ZMOLECNO,ZINVCON)
          ZMOLECNO2     = Conc(JL,JK,ind_NO2)
          ZMOLECNO2     = MAX(ZMOLECNO2,ZINVCON)
          ZMOLECNO3     = Conc(JL,JK,ind_NO3)
          ZMOLECNO3     = MAX(ZMOLECNO3,ZINVCON)
          ZMOLECN2O5    = Conc(JL,JK,ind_N2O5)
          ZMOLECN2O5    = MAX(ZMOLECN2O5,ZINVCON)
          ZMOLECHNO4    = Conc(JL,JK,ind_HNO4)
          ZMOLECHNO4    = MAX(ZMOLECHNO4,ZINVCON)
          ZMOLECCL      = Conc(JL,JK,ind_CL)
          ZMOLECCL      = MAX(ZMOLECCL,ZINVCON)
          ZMOLECCLO     = Conc(JL,JK,ind_CLO)
          ZMOLECCLO     = MAX(ZMOLECCLO,ZINVCON)
          ZMOLECCLOH    = Conc(JL,JK,ind_HOCl)
          ZMOLECCLOH    = MAX(ZMOLECCLOH,ZINVCON)
          ZMOLECCL2O2   = Conc(JL,JK,ind_CL2O2)
          ZMOLECCL2O2   = MAX(ZMOLECCL2O2,ZINVCON)
          ZMOLECCL2     = Conc(JL,JK,ind_CL2)
          ZMOLECCL2     = MAX(ZMOLECCL2,ZINVCON)
          ZMOLECCH2O    = Conc(JL,JK,ind_HCHO)
          ZMOLECCH2O    = MAX(ZMOLECCH2O,ZINVCON)
          ZMOLECCH3O2   = Conc(JL,JK,ind_CH3O2)
          ZMOLECCH3O2   = MAX(ZMOLECCH3O2,ZINVCON)


          ZMOLECCH3O    = 0._dp
          ZMOLECNAT     = Conc(JL,JK,ind_NAT)
          ZMOLECCH4     = Conc(JL,JK,ind_CH4)
          ZMOLECCH4     = MAX(ZMOLECCH4,1._dp)
          ZMOLECH2      = Conc(JL,JK,ind_H2)
          ZMOLECH2      = MAX(ZMOLECH2,1._dp)
          ZMOLECN2O     = Conc(JL,JK,ind_N2O)
          ZMOLECN2O     = MAX(ZMOLECN2O,1._dp)
          ZMOLECH2O2    = Conc(JL,JK,ind_H2O2)
          ZMOLECH2O2    = MAX(ZMOLECH2O2,1._dp)
          ZMOLECHCL     = Conc(JL,JK,ind_HCL)
          ZMOLECHCL     = MAX(ZMOLECHCL,1._dp)
          ZMOLECCO      = Conc(JL,JK,ind_CO)
          ZMOLECCO      = MAX(ZMOLECCO,1._dp)
          ZMOLECCH3O2H  = Conc(JL,JK,ind_CH3OOH)
          ZMOLECCH3O2H  = MAX(ZMOLECCH3O2H,1._dp)
          ZMOLECICE     = Conc(JL,JK,ind_ICE)
          ZMOLECICE     = MAX(ZMOLECICE,1._dp)
          ZMOLECCLNO3   = Conc(JL,JK,ind_CLNO3)
          ZMOLECCLNO3   = MAX(ZMOLECCLNO3,TINY(0._dp))
          ZMOLECHNO3    = Conc(JL,JK,ind_HNO3)
          ZMOLECHNO3    = MAX(ZMOLECHNO3,1._dp)
          ZMOLECCFCL3   = Conc(JL,JK,ind_CFCL3)
          ZMOLECCFCL3   = MAX(ZMOLECCFCL3,0._dp)
          ZMOLECCF2CL2  = Conc(JL,JK,ind_CF2CL2)
          ZMOLECCF2CL2  = MAX(ZMOLECCF2CL2,0._dp)
          ZMOLECCH3CL   = Conc(JL,JK,ind_CH3CL)
          ZMOLECCH3CL   = MAX(ZMOLECCH3CL,0._dp)
          ZMOLECCCL4    = Conc(JL,JK,ind_CCL4)
          ZMOLECCCL4    = MAX(ZMOLECCCL4 ,0._dp)
          ZMOLECCH3CCL3 = Conc(JL,JK,ind_CH3CCL3)
          ZMOLECCH3CCL3 = MAX(ZMOLECCH3CCL3,0._dp)


          ZMOLECOX      = ZMOLECO3 + ZMOLECO1D + ZMOLECO ! mz_ab_20090912 instead of transport scaling
          ! mz_ab_20090915+ Rescue O3/OX partitioning for partitioning at the end
          ZSCALEO3(JL)      = ZMOLECO3/ZMOLECOX
          ! mz_ab_20090915-


          ZMOLECNOX         = ZMOLECN + ZMOLECNO + ZMOLECNO2 + &
                              ZMOLECNO3 + 2. * ZMOLECN2O5 + ZMOLECHNO4
          ZMOLECNOX         = ZMOLECNOX + ZMOLECCLNO3

          ZMOLECCLOX        = ZMOLECCL + ZMOLECCLO + ZMOLECCLOH + &
                               2. * (ZMOLECCL2O2 + ZMOLECCL2)

          ZMOLECCLOX        = ZMOLECCLOX + ZMOLECCLNO3

          ZMOLECH2O         = Conc(JL,JK,ind_H2O) !

          ! mz_bk_20101216+
!!$          ZMOLECH2O         = MAX(ZMOLECH2O,0._dp)
          IF (ZMOLECH2O < -1.e-12_dp) &
               CALL warning("ZMOLECH2O negative (< -1.E12)! Set to 1.E-12!" &
               , modstr)
          ZMOLECH2O         = MAX(ZMOLECH2O,1.e-12_dp)
          ! mz_bk_20101216-

          ZMOLECH2OBFC = ZMOLECH2O
     
          !     T0-WERTE SETZEN
          !
          ZMOLECH2OT0V(JL)     = ZMOLECH2O

          !---------------------------------------------------------------------
          !
          !     HETEROGENE CHEMIE TEIL1
          !
          !     HETEROGENE REAKTIONEN
          !
          !     FOLGENDES VORGEHEN :
          !               ZUNAECHST SO TUN ALS OB ALLES IN GASPHASE!
          !               DANN, BEVOR DIE CHEMIE GERECHNET WIRD, UEBERGLEICHGE-
          !               WICHTSDAMPFDRUCK BERECHNEN WIEVIEL HNO3 IM NAT.
          !               DIE NAEHERUNG AUCH DAS AUF EIS PRODUZIERTE HNO3 IN
          !               NAT ZU STECKEN IST WOHL ERLAUBT. DA DIE REAKTIONSRATEN
          !               AUF EIS SO GROSS SIND WG. 1. GROESSEREN GAMMAS
          !                                         2. UM 1 GROESSENORDNUNG
          !                                            GROESSERE OBERFLAECHE.
          !               IM 3D-MODEL SOLL NAT SO SCHNELL FALLEN WIE EIS, WENN
          !               EIS DA IST!
          !               BERUECKSICHTIGE BEI REAKTIONSRATEN ABER NICHT DIE
          !               VERAENDERTE OBERFLAECHE (NAT WAECHST LANGSAM).
          !
          ZCP1 = 0._dp
          ZCP2 = 0._dp
          ZCP3 = 0._dp
          ZCP4 = 0._dp
          ZCS1 = 0._dp
          ZCS2 = 0._dp
          INIHET0 = 0
          !
          ZCPI1 = 0._dp
          ZCPI2 = 0._dp
          ZCPI3 = 0._dp
          ZCPI4 = 0._dp
          ZICESURF = 0._dp
          !
          ZCPN1 = 0._dp
          ZCPN2 = 0._dp
          ZCPN3 = 0._dp
          ZCPN4 = 0._dp
          ZNATSURF = 0._dp
          !
          ! ******************************************************************
          ! *                                                                *
          ! * if het. reactions on ice permitted:                            *
          ! *                                                                *
          ! ******************************************************************
          IF ((IIHET(JL).EQ.1).and.(REAL(JK,dp).lt.tp_i(jl)).and. &
               (zpres.ge.zriceup(JL))) THEN
             !
             !     FUER TESTZWECKE :
             !
             !     ICE, SURF=3.0E-07 CM^2/CM^3
             !       ZICESURF  = 3.0E-07
             !
             !     EISBILDUNG
             !
             !     ENTSCHEIDUNG OB UEBERHAUPT EIS DA.
             !     FALLS KEINE EISTEILCHEN DA : TEMPNUCLI = TEMP + 1.8
             !     FALLS       EISTEILCHEN DA : TEMPNUCLI = TEMP
             !     1.8 K ENTSPRECHEN UNGEFAEHR 1.3 FACHER UEBERSAETTIGUNG
             !
             !     KRITERIUM FUER EISTEILCHEN (WG. TRANSPORT) :
             !        EISTEILCHEN DUERFEN NICHT KLEINER ALS IHR KERN WERDEN !
             !        KERN IST NAT UND HAT GROESSE VON 0.3 MUE !
             !               (ENTSPRICHT ETWA 3.5E+08 MOLEC/CM^3 EIS BEI UNTIGER
             !                TEILCHENZAHLDICHTE)
             !     ANGENOMMENE EISTEILCHENZAHLDICHTE [TEILCHEN/CM^3]
             !      ZNUMICE  = 0.03536
             !     EISDICHTE : RHO = 0.92 G/CM^3
             !      RHOICE   = 0.92
             !        3.2608E-22 = 18./(6.0E+23*0.92*0.1)
             !      ZHELPICE = 18./(6.0E+23*RHOICE*ZNUMICE)
             !
             !       IF (ZMOLECICE.LT.3.5E+08) THEN
             !        ZTEMPNUCLI = ZTEMP + 1.8
             ZTEMPNUCLI = ZTEMP
             !       ELSE
             !        ZTEMPNUCLI = ZTEMP
             !       ENDIF
             !
             ZMOLECTOTH2O  = ZMOLECICE + ZMOLECH2OBFC
             ZMOLECEQH2O   = 7.2427E+18_dp/ZTEMPNUCLI* &
                  EXP(24.306_dp-6144.9_dp/ZTEMPNUCLI)
             ZDIFFH2O      = ZMOLECTOTH2O - ZMOLECEQH2O
             IF(ZDIFFH2O.GT.0._dp) THEN
                ZMOLECICE    = ZDIFFH2O
                ZMOLECH2O    = MAX(ZMOLECEQH2O,1._dp)
                !        ZVH2O        = 3.2608E-22*ZMOLECICE
                ZVH2O        = ZHELPICE(JK)*ZMOLECICE
                ZRADICE      = (0.2387324_dp*ZVH2O)**0.333333_dp
                ZICESURF     = 12.566371_dp*ZRADICE*ZRADICE
                !
                ZICEDENOM1   = 3.3E02_dp*ZRADICE*ZPRES/ZTEMP
                ZICEDENOM2   = 1._dp+ZICEDENOM1*0.3_dp
                !        ZCPI1 = 0.3    * 3628.7327 *SQRT(ZTEMP/97.5) * ZICESURF * 0.1
                ZCPI1 = 0.3_dp*3628.7327_dp*SQRT(ZTEMP/97.5_dp)*ZICESURF*ZNUMICE(JK) &
                     /ZICEDENOM2
                !        ZCPI2 = 0.3    * 3628.7327 *SQRT(ZTEMP/97.5) * ZICESURF * 0.1
                ZCPI2 = 0.3_dp*3628.7327_dp*SQRT(ZTEMP/97.5_dp)*ZICESURF*ZNUMICE(JK) &
                     /ZICEDENOM2
                !        ZCPI3 = 0.3    * 3628.7327 *SQRT(ZTEMP/52.5) * ZICESURF * 0.1
                ZCPI3 = 0.3_dp*3628.7327_dp*SQRT(ZTEMP/52.5_dp)*ZICESURF*ZNUMICE(JK) &
                     /ZICEDENOM2
                ZICEDENOM2   = 1._dp+ZICEDENOM1*0.03_dp
                !        ZCPI4 = 0.03   * 3628.7327 *SQRT(ZTEMP/108.) * ZICESURF * 0.1
                ZCPI4 = 0.03_dp*3628.7327_dp*SQRT(ZTEMP/108._dp)*ZICESURF*ZNUMICE(JK) &
                     /ZICEDENOM2
                !
                ZCP1  = ZCPI1
                ZCP2  = ZCPI2
                ZCP3  = ZCPI3
                ZCP4  = ZCPI4
                INIHET0 = 1
                !
                !     SEDIMENTATION !
                !
                !
                ZUDT         = 22085558._dp*ZICESURF/ZTEMP*ZTMSTDT
                ZUDT = MIN(ZUDT,0.5_dp*ZHEILEV(JL))
                ZUDTAICE(JL) = ZUDT*ZMOLECICE
             ELSE
                ZMOLECICE    = 0._dp
                ZMOLECH2O    = ZMOLECTOTH2O
                ZUDTAICE(JL) = 0._dp
                ZUDT         = 0._dp
             ENDIF
          ELSE
             ! ******************************************************************
             ! *                                                                *
             ! * if het. reactions on ice not permitted:                        *
             ! *                                                                *
             ! ******************************************************************
             ! mz_ab_20100305+ Benedikt and Andreas decided that this does not make sense
!             ZMOLECH2O = ZMOLECICE + ZMOLECH2O
!             ZMOLECICE = 0._dp             
             ! mz_ab_20100305-
             ZUDTAICE(JL) = 0._dp
             ZUDT         = 0._dp
          ENDIF
          !
!!$          PXTRSDICESURF(JL,JK) = PXTRSDICESURF(JL,JK)+ & ! mz_ab_20091007 not used
!!$                                 ZTMSTDT*ZICESURF
          !
          !
          ! ******************************************************************
          ! *                                                                *
          ! * if het. reactions on NAT permitted:                            *
          ! *                                                                *
          ! ******************************************************************
          IF ((INHET(JL).EQ.1).and.(REAL(JK,dp).lt.tp_i(jl)).and.&
               (zpres.ge.zrnatup(JL))) THEN
             !
             !     NATBILDUNG
             !
             !     ENTSCHEIDUNG OB UEBERHAUPT NAT DA.
             !     FALLS KEINE NATTEILCHEN DA : TEMPNUCL = TEMP + 3
             !     FALLS       NATTEILCHEN DA : TEMPNUCL = TEMP
             !     3 K ENTSPRECHEN UNGEFAEHR 10 FACHER UEBERSAETTIGUNG
             !
             !     KRITERIUM FUER NATTEILCHEN (WG. TRANSPORT) :
             !        NATTEILCHEN DUERFEN NICHT KLEINER ALS IHR KERN WERDEN !
             !        KERN IST SCHWEFELAEROSOL UND HAT GROESSE VON 0.07 MUE !
             !               (ENTSPRICHT ETWA 11.8E+06 MOLEC/CM^3 NAT BEI UNTIGER
             !                TEILCHENZAHLDICHTE)
             !     ANGENOMMENE NATTEILCHENZAHLDICHTE : 1 TEILCHEN/CM^3
             !     ANGENOMMENE NATDICHTE : RHO = 1.6 G/CM^3 M(NAT)=117
             !     1.2188E-22 = 117/(1.6*6.0E+23*1.)
             !
             !     FUER  TESTZWECKE: ZNATSURF = 1.0E-07 CM^2/CM^3
             !
             !      ZCORENAT = 7.0E-06
             !      ZVHNO3 = ZMOLECNAT * 1.2188E-22
             !      ZRADNAT = (0.2387324*ZVHNO3)**0.3333333
             !       ZNATSURF     = 1.0E-07
             !
             !       IF (ZMOLECNAT.LT.11.8E+06) THEN
             !        ZTEMPNUCLN = ZTEMP + 3.
             ZTEMPNUCLN = ZTEMP
             !       ELSE
             !        ZTEMPNUCLN = ZTEMP
             !       ENDIF
             ZMOLECTOTHNO3 = ZMOLECNAT + ZMOLECHNO3
             ZMOLECEQHNO3  = 7.2427E+18/ZTEMPNUCLN* &
                             EXP((-2.7836_dp-0.00088_dp*ZTEMPNUCLN)* &
                             LOG(ZMOLECH2O*ZTEMPNUCLN/7.2427E+18_dp)+ &
                             90.8556_dp-26242.6_dp/ZTEMPNUCLN+ &
                             0.0213885_dp*ZTEMPNUCLN)
             ZDIFFHNO3     = ZMOLECTOTHNO3 - ZMOLECEQHNO3
             IF(ZDIFFHNO3.GT.0._dp) THEN
                ZMOLECNAT    = ZDIFFHNO3
                ZMOLECHNO3   = MAX(ZMOLECEQHNO3,1._dp)
                ZVHNO3       = 1.2188E-22*ZMOLECNAT
                ZRADNAT      = (0.2387324_dp*ZVHNO3)**0.333333_dp
                ZNATSURF     = 12.566371_dp*ZRADNAT*ZRADNAT
                !
                ZNATDENOM1   = 3.3E02_dp*ZRADNAT*ZPRES/ZTEMP
                ZNATDENOM2   = 1._dp+ZNATDENOM1*0.3_dp
                ZCPN1 = 0.3_dp    * 3628.7327_dp *SQRT(ZTEMP/97.5_dp) * ZNATSURF &
                     /ZNATDENOM2
                ZNATDENOM2   = 1._dp+ZNATDENOM1*0.006
                ZCPN2 = 0.006_dp  * 3628.7327_dp *SQRT(ZTEMP/97.5_dp) * ZNATSURF &
                     /ZNATDENOM2
                ZNATDENOM2   = 1._dp+ZNATDENOM1*0.1_dp
                ZCPN3 = 0.1_dp    * 3628.7327_dp *SQRT(ZTEMP/52.5_dp) * ZNATSURF &
                     /ZNATDENOM2
                ZNATDENOM2   = 1._dp+ZNATDENOM1*0.0006_dp
                ZCPN4 = 0.0006_dp * 3628.7327_dp *SQRT(ZTEMP/108._dp) * ZNATSURF &
                     /ZNATDENOM2
                !
                ZCP1  = ZCP1 + ZCPN1
                ZCP2  = ZCP2 + ZCPN2
                ZCP3  = ZCP3 + ZCPN3
                ZCP4  = ZCP4 + ZCPN4
                INIHET0 = 1
                !
                !     SEDIMENTATION
                !
                IF(INSED(JL).EQ.1) THEN
                   !         ZUDTANAT(JL) = ZUDT*ZMOLECNAT*0.5
                   ZUDTANAT(JL) = ZUDT*ZMOLECNAT
                ELSE
                   ZUDTANAT(JL) = 0._dp
                ENDIF
             ELSE
                ZMOLECNAT    = 0._dp
                ZMOLECHNO3   = ZMOLECTOTHNO3
                ZUDTANAT(JL) = 0._dp
             ENDIF
          ELSE
             ! ******************************************************************
             ! *                                                                *
             ! * if het. reactions on NAT not permitted:                        *
             ! *                                                                *
             ! ******************************************************************
             ZMOLECHNO3   = ZMOLECNAT + ZMOLECHNO3
             ZMOLECNAT    = 0._dp
             ZUDTANAT(JL) = 0._dp
          ENDIF
          !
!!$          PXTRSDNATSURF(JL,JK) = PXTRSDNATSURF(JL,JK)+ & ! mz_ab_20091007 not used
!!$                                 ZTMSTDT*ZNATSURF
          !
          ZMOLECODD=ZMOLECOX+ &
               ZMOLECCLOH+ZMOLECCLO+ &
               2._dp*(ZMOLECCLNO3+ZMOLECCL2O2+ZMOLECNO3)+ &
               ZMOLECNO2+3._dp*ZMOLECN2O5+ZMOLECHNO4+ &
               ZMOLECHNO3
          !
          !     SULFAT, SURF=1.0E-08 CM^2/CM^3
          !
          !sulnew start
          IF ((ISHET.EQ.1).and.(REAL(JK,dp).lt.tp_i(jl)).and. &
               (zpres.ge.rsulup).and.(zpres.le.rsuldo)) THEN
             JSUL      = MAX(MIN(INT(ZPRES/100._dp),ksul),1)
             ZSULSURF  = SULOOK(JSUL,JL,jrow)
             ZVOL0     = SQRT(ZSULSURF*ZSULSURF*ZSULSURF)*0.0437046_dp
             ZSCALSURF = ZWEIGHT0*ZVOL0*ZRHO0
             !sulnew end
             ZPH2O     = LOG(ZMOLECH2O/ZCON * ZPRES * .01_dp)
             ZWEIGHT   = ZTEMP/3674.3_dp*(23.1_dp-ZPH2O)-1.591_dp
             IF(ZWEIGHT.LE.0._dp) THEN
                ZWEIGHT   = 0.247_dp
             ELSE
                ZWEIGHT   = MIN(0.552_dp,ZWEIGHT)
                ZWEIGHT   = 0.247_dp+SQRT(ZWEIGHT)
             ENDIF
             ZCSS      = 10.204081_dp*ZWEIGHT/(1._dp-ZWEIGHT)
             ZRHO      = 1._dp+ZCSS*(0.12364_dp- &
                  ZTEMP*ZTEMP*(5.6E-07_dp+1.324E-08*ZCSS) &
                  +SQRT(ZCSS)*(1.814E-07*ZTEMP*ZTEMP-0.02954_dp) &
                  +ZCSS*(0.002343_dp-1.487E-06*ZTEMP))
             ZRHO      = ABS(ZRHO)
             ZSULSURF  = (ZSCALSURF/ZWEIGHT/ZRHO)**0.666666666*8.0596316
             ZGAMMAS1  = 10._dp**(1.86_dp-7.47_dp*ZWEIGHT)
             ZCS1 = ZGAMMAS1 * 3628.7327_dp *SQRT(ZTEMP/97.5_dp) * ZSULSURF
             ZCS2 = 0.1_dp      * 3628.7327_dp *SQRT(ZTEMP/108._dp) * ZSULSURF
          ENDIF
          !
          !---------------------------------------------------------------------
          !
          !     HETEROGENE CHEMIE TEIL2
          !
          !     FOLGENDE RECHNUNG IMMER UNTER DER ANNAHME, DASS NACH REAKTIONEN
          !     AUF SCHWEFELAEROSOLEN HNO3 IN DIE GASPHASE GEHT!
          !     SCHWEFELAEROSOL WIRD BEI NAT-TEMPERATUR NICHT ABGESCHALTET!
          !
          !     WAS NOX HETEROGEN VERLIERT IST NACH AUFSPALTUNG BEKANNT!
          !     SOMIT AUCH WAS [NAT + HNO3] GEWINNEN. DIESER GEWINN WIRD IM
          !     VERHAELTNIS DER HETEROGENEN REAKTIONSRATEN ZUGETEILT.
          !
          !
          IF (ZCP1.GT.0.1E-15_dp) THEN
             ZHEHCL1     = (1._dp-EXP(-1._dp*ZCP1*ZTMSTDT))*ZMOLECCLNO3
             ZHEHCL2     = (1._dp-EXP(-1._dp*ZCP3*ZTMSTDT))*ZMOLECCLOH
             ZHEHCL      = ZHEHCL1+ZHEHCL2
             IF( ZHEHCL.GT.ZMOLECHCL) THEN
                ZHEHCL     = ZMOLECHCL
                ZALPH      = ZCP3/ZCP1
                ZLOW       = 1._dp
                ZHIGH      = 0._dp
                ZNEU       = (ZLOW+ZHIGH)*0.5_dp
                !
                !        1.
                !
                ZEQUA      = (1._dp-ZNEU)*ZMOLECCLNO3+(1._dp-ZNEU**ZALPH)*ZMOLECCLOH
                IF (ZEQUA.GT.ZMOLECHCL) THEN
                   ZHIGH = ZNEU
                ELSE
                   ZLOW = ZNEU
                ENDIF
                ZNEU       = (ZLOW+ZHIGH)*0.5_dp
                !
                !        2.
                !
                ZEQUA      = (1._dp-ZNEU)*ZMOLECCLNO3+(1._dp-ZNEU**ZALPH)*ZMOLECCLOH
                IF (ZEQUA.GT.ZMOLECHCL) THEN
                   ZHIGH = ZNEU
                ELSE
                   ZLOW = ZNEU
                ENDIF
                ZNEU       = (ZLOW+ZHIGH)*0.5_dp
                !
                !        3.
                !
                ZEQUA      = (1._dp-ZNEU)*ZMOLECCLNO3+(1._dp-ZNEU**ZALPH)*ZMOLECCLOH
                IF (ZEQUA.GT.ZMOLECHCL) THEN
                   ZHIGH = ZNEU
                ELSE
                   ZLOW = ZNEU
                ENDIF
                ZNEU       = (ZLOW+ZHIGH)*0.5_dp
                !
                !        4.
                !
                ZEQUA      = (1._dp-ZNEU)*ZMOLECCLNO3+(1._dp-ZNEU**ZALPH)*ZMOLECCLOH
                IF (ZEQUA.GT.ZMOLECHCL) THEN
                   ZHIGH = ZNEU
                ELSE
                   ZLOW = ZNEU
                ENDIF
                ZNEU       = (ZLOW+ZHIGH)*0.5
                !
                !        5.
                !
                ZEQUA      = (1._dp-ZNEU)*ZMOLECCLNO3+(1._dp-ZNEU**ZALPH)*ZMOLECCLOH
                IF (ZEQUA.GT.ZMOLECHCL) THEN
                   ZHIGH = ZNEU
                ELSE
                   ZLOW = ZNEU
                ENDIF
                ZNEU       = (ZLOW+ZHIGH)*0.5_dp
                !
                !        6.
                !
                ZEQUA      = (1._dp-ZNEU)*ZMOLECCLNO3+(1._dp-ZNEU**ZALPH)*ZMOLECCLOH
                IF (ZEQUA.GT.ZMOLECHCL) THEN
                   ZHIGH = ZNEU
                ELSE
                   ZLOW = ZNEU
                ENDIF
                ZNEU       = (ZLOW+ZHIGH)*0.5_dp
                !
                !        7.
                !
                ZEQUA      = (1._dp-ZNEU)*ZMOLECCLNO3+(1._dp-ZNEU**ZALPH)*ZMOLECCLOH
                IF (ZEQUA.GT.ZMOLECHCL) THEN
                   ZHIGH = ZNEU
                ELSE
                   ZLOW = ZNEU
                ENDIF
                ZNEU       = (ZLOW+ZHIGH)*0.5_dp
                !
                !        8.
                !
                ZEQUA      = (1._dp-ZNEU)*ZMOLECCLNO3+(1._dp-ZNEU**ZALPH)*ZMOLECCLOH
                IF (ZEQUA.GT.ZMOLECHCL) THEN
                   ZHIGH = ZNEU
                ELSE
                   ZLOW = ZNEU
                ENDIF
                ZNEU       = (ZLOW+ZHIGH)*0.5_dp
                !
                !        9.
                !
                ZEQUA      = (1._dp-ZNEU)*ZMOLECCLNO3+(1._dp-ZNEU**ZALPH)*ZMOLECCLOH
                IF (ZEQUA.GT.ZMOLECHCL) THEN
                   ZHIGH = ZNEU
                ELSE
                   ZLOW = ZNEU
                ENDIF
                ZNEU       = (ZLOW+ZHIGH)*0.5_dp
                !
                !        10.
                !
                ZEQUA      = (1._dp-ZNEU)*ZMOLECCLNO3+(1._dp-ZNEU**ZALPH)*ZMOLECCLOH
                IF (ZEQUA.GT.ZMOLECHCL) THEN
                   ZHIGH = ZNEU
                ELSE
                   ZLOW = ZNEU
                ENDIF
                !
                ZCP1 = -1._dp*LOG(ZLOW)/ZTMSTDT
                ZCP3 = ZALPH*ZCP1
             ENDIF
          ELSE
             ZHEHCL = 0._dp
          ENDIF
          !
          ZFMOLECHCL  = ZMOLECHCL - ZHEHCL
          ZFMOLECCLOX = ZMOLECCLOX + ZHEHCL
          ZFMOLECNOX = ZMOLECNOX- &
               (1._dp-EXP(-1._dp*(ZCP1+ZCP2+ZCS1)*ZTMSTDT))*ZMOLECCLNO3- &
               (1._dp-EXP(-1._dp*(ZCP4+ZCS2)*ZTMSTDT))*2._dp*ZMOLECN2O5
          ZFMOLECHNO3= ZMOLECHNO3+(1._dp-EXP(-1._dp*ZTMSTDT*ZCS1))*ZMOLECCLNO3+ &
               2._dp*(1._dp-EXP(-1._dp*ZTMSTDT*ZCS2))*ZMOLECN2O5
          ZMOLECO3V(JL)   =  ZMOLECO3
          ZMOLECOXV(JL)   =  ZMOLECOX
          ZMOLECHV(JL)    =  ZMOLECH
          ZMOLECOHV(JL)   =  ZMOLECOH
          ZMOLECHO2V(JL)  =  ZMOLECHO2
          ZMOLECNV(JL)    =  ZMOLECN
          ZMOLECNOV(JL)   =  ZMOLECNO
          ZMOLECNO2V(JL)  =  ZMOLECNO2
          ZMOLECNO3V(JL)  =  ZMOLECNO3
          ZMOLECHNO4V(JL) =  ZMOLECHNO4
          ZMOLECN2O5V(JL) =  ZMOLECN2O5
          ZMOLECCLV(JL)   =  ZMOLECCL
          ZMOLECCLOV(JL)  =  ZMOLECCLO
          ZMOLECCLOHV(JL) =  ZMOLECCLOH
          ZMOLECCL2O2V(JL)=  ZMOLECCL2O2
          ZMOLECCLNO3V(JL)=  ZMOLECCLNO3
          ZMOLECCL2V(JL)  =  ZMOLECCL2
          ZMOLECCH2OV(JL) =  ZMOLECCH2O
          ZMOLECCH3O2V(JL)=  ZMOLECCH3O2
          ZMOLECCH3OV(JL) =  ZMOLECCH3O
          ZMOLECCFCL3V(JL)=  ZMOLECCFCL3
          ZMOLECCF2CL2V(JL)= ZMOLECCF2CL2
          ZMOLECCCL4V(JL) =  ZMOLECCCL4
          ZMOLECCH3CLV(JL)=  ZMOLECCH3CL
          ZMOLECCH3CCL3V(JL)= ZMOLECCH3CCL3
          ZMOLECH2V(JL)   =  ZMOLECH2
          ZMOLECCH4V(JL)  =  ZMOLECCH4
          !
          ZMOLECN2OV(JL)  =  ZMOLECN2O
          !
          ZMOLECH2O2V(JL)  =  ZMOLECH2O2
          !
          ZMOLECHCLV(JL)  =  ZMOLECHCL
          !      ZMOLECHCLOLDV(JL) = ZMOLECHCL
          !
          ZMOLECCOV(JL)   =  ZMOLECCO
          ZMOLECCH3O2HV(JL)= ZMOLECCH3O2H
          ZMOLECNATV(JL)  =  ZMOLECNAT
          ZMOLECICEV(JL)  =  ZMOLECICE
          !
          ZMOLECHNO3V(JL) =  ZMOLECHNO3
          !      ZMOLECHNO3OLDV(JL) = ZMOLECHNO3
          !
          ZMOLECCLOXV(JL) =  ZMOLECCLOX
          !      ZMOLECCLOXOLDV2(JL) = ZMOLECCLOX
          ZMOLECNOXV(JL)  =  ZMOLECNOX
          !      ZMOLECNOXOLDV2(JL) = ZMOLECNOX
          !
          ZMOLECODDV(JL)  =  ZMOLECODD
          ZFMOLECNOXV(JL) =  ZFMOLECNOX
          ZFMOLECHNO3V(JL)=  ZFMOLECHNO3
          ZFMOLECHCLV(JL) =  ZFMOLECHCL
          ZFMOLECCLOXV(JL)=  ZFMOLECCLOX
          ZMOLECH2OV(JL) = ZMOLECH2O
          !
          ZCP1V(JL)       =  ZCP1
          ZCP2V(JL)       =  ZCP2
          ZCP3V(JL)       =  ZCP3
          ZCP4V(JL)       =  ZCP4
          ZCS1V(JL)       =  ZCS1
          ZCS2V(JL)       =  ZCS2
          !
          INIHET0V(JL)    =  INIHET0
          !
          ! ******************************************************************
          ! *                                                                *
          ! * end loop over longitudes                                       *
          ! *                                                                *
          ! ******************************************************************
       ENDDO

       ! ******************************************************************
       ! *                                                                *
       ! * start loops over longitudes                                    *
       ! * (homogenous chemistry)                                         *
       ! *                                                                *
       ! ******************************************************************
       !
       DO JL=1,KLON
          !
          ZCON          = ZCONV(JL)
          ZTEMP         =  ZTEMPV(JL)
          !
          ZMOLECO2V(JL) = 0.2095_dp*ZCON
          ZMOLECN2V(JL) = ZCON-ZMOLECO2V(JL)
          zmolecco2V(JL)= Conc(JL,JK,ind_CO2)
          !
          !     CALCULATE TEMPERATURE INDEX JIND OF LOOKUP-TABLE
          !     FOR GAS-PHASE REACTION RATES
          !
          !     T is forced to be within (IRCTMIN,IRCTMAX)
          !
          ZTEIND = MIN(ZTEMP,REAL(IRCTMAX,dp))
          ZTEIND = MAX(ZTEIND,REAL(IRCTMIN,dp))
          !
          !     so ist's kuerzer und hoffentlich richtig:
          !
          zjind=1._dp+(zteind-REAL(irctmin))*REAL(jvals)
          jind=int(zjind+0.5_dp)
          JINDV(JL) = JIND
          !
          IF(DANI(JL).GT.1._dp) THEN
             !
             ZMOLECO3      =  ZMOLECO3V(JL)
             ZMOLECOXT0    =  ZMOLECOXV(JL)
             ZMOLECH2O     =  ZMOLECH2OV(JL)
             ZMOLECNO2     =  ZMOLECNO2V(JL)
             ZMOLECNO3     =  ZMOLECNO3V(JL)
             !
             !     SET PHOTOLYSES-RATES!
             !
             ZDO3B      = RJ_O3P(JL,JK)
             ZDO3A      = RJ_O1D(JL,JK)
             ZDNO3A     = RJ_NO2O(JL,JK)
             ZDNO2      = RJ_NO2(JL,JK)
             !
             ZCO4        = RCGAS(JIND,2)
             ZCO5        = RCGAS(JIND,3)
             ZCO3        = RCGAS(JIND,37)*ZCON
             ZCH8        = RCGAS(JIND,65)
             !
             !     FIRST ESTIMATE OF [OX]=[O1D]+[O]+[O3];
             !     CALCULATED FROM STEADY-STATE ASSUMPTIONS FOR [O1D] AND [O]
             !
             ZSTO      = ZCO4*ZMOLECO2V(JL)+ZCO5*ZMOLECN2V(JL)
             ZRATO1DO3 = ZDO3A/(ZSTO+ZCH8*ZMOLECH2O)
             ZRATOO3   = (ZDO3B+(ZDNO2*ZMOLECNO2+ZDNO3A*ZMOLECNO3)/ZMOLECO3+ &
                  ZSTO*ZRATO1DO3)/(ZCO3*ZMOLECO2V(JL))
             ZRATO3OX  = 1._dp/(1._dp+ZRATO1DO3+ZRATOO3)
             ZRATOOX   = ZRATO3OX*ZRATOO3
             ZRATO1DOX = ZRATO3OX*ZRATO1DO3
             !
             ZMOLECO3V(JL)  = ZRATO3OX*ZMOLECOXT0
             ZMOLECO1DV(JL) = ZRATO1DOX*ZMOLECOXT0
             ZMOLECOV(JL)   = ZRATOOX*ZMOLECOXT0
          ENDIF
       ENDDO
       !
       DO JL=1,KLON
          !
          ZTEMP         =  ZTEMPV(JL)
          ZCON          =  ZCONV(JL)
          ZINVCON       =  ZINVCONV(JL)
          ZPRES         =  ZPRESV(JL)
          !
          ZMOLECO3      =  ZMOLECO3V(JL)
          ZMOLECOX      =  ZMOLECOXV(JL)
          ZMOLECH       =  ZMOLECHV(JL)
          ZMOLECOH      =  ZMOLECOHV(JL)
          ZMOLECHO2     =  ZMOLECHO2V(JL)
          ZMOLECN       =  ZMOLECNV(JL)
          ZMOLECNO      =  ZMOLECNOV(JL)
          ZMOLECNO2     =  ZMOLECNO2V(JL)
          ZMOLECNO3     =  ZMOLECNO3V(JL)
          ZMOLECHNO4    =  ZMOLECHNO4V(JL)
          ZMOLECN2O5    =  ZMOLECN2O5V(JL)
          ZMOLECCL      =  ZMOLECCLV(JL)
          ZMOLECCLO     =  ZMOLECCLOV(JL)
          ZMOLECCLOH    =  ZMOLECCLOHV(JL) 
          ZMOLECCL2O2   =  ZMOLECCL2O2V(JL)
          ZMOLECCLNO3   =  ZMOLECCLNO3V(JL)
          ZMOLECCL2     =  ZMOLECCL2V(JL)
          ZMOLECCH2O    =  ZMOLECCH2OV(JL)
          ZMOLECCH3O2   =  ZMOLECCH3O2V(JL)
          ZMOLECCH3O    =  ZMOLECCH3OV(JL)
          ZMOLECCFCL3   =  ZMOLECCFCL3V(JL)
          ZMOLECCF2CL2  =  ZMOLECCF2CL2V(JL)
          ZMOLECCCL4    =  ZMOLECCCL4V(JL)
          ZMOLECCH3CL   =  ZMOLECCH3CLV(JL)
          ZMOLECCH3CCL3 =  ZMOLECCH3CCL3V(JL)
          ZMOLECH2      =  ZMOLECH2V(JL)
          ZMOLECCH4     =  ZMOLECCH4V(JL)
          ZMOLECN2O     =  ZMOLECN2OV(JL)
          ZMOLECH2O2    =  ZMOLECH2O2V(JL)
          ZMOLECHCL     =  ZMOLECHCLV(JL)
          ZMOLECCO      =  ZMOLECCOV(JL)
          ZMOLECCH3O2H  =  ZMOLECCH3O2HV(JL)
          ZMOLECNAT     =  ZMOLECNATV(JL)
          ZMOLECICE     =  ZMOLECICEV(JL)
          ZMOLECHNO3    =  ZMOLECHNO3V(JL)
          ZMOLECCLOX    =  ZMOLECCLOXV(JL)
          ZMOLECNOX     =  ZMOLECNOXV(JL)
          ZMOLECODD     =  ZMOLECODDV(JL)
          ZFMOLECNOX    =  ZFMOLECNOXV(JL)
          ZFMOLECHNO3   =  ZFMOLECHNO3V(JL)
          ZFMOLECHCL    =  ZFMOLECHCLV(JL)
          ZFMOLECCLOX   =  ZFMOLECCLOXV(JL)
          ZMOLECH2O     =  ZMOLECH2OV(JL)
          !
          ZCP1       =  ZCP1V(JL)
          ZCP2       =  ZCP2V(JL)
          ZCP3       =  ZCP3V(JL)
          ZCP4       =  ZCP4V(JL)
          ZCS1       =  ZCS1V(JL)
          ZCS2       =  ZCS2V(JL)
          ZHET1      = ZCP4+ZCS2
          ZHET2      = ZCP1+ZCP2+ZCS1
          ZHET3      = ZCP2+ZCS1
          !
          INIHET0    =  INIHET0V(JL)
          !
          !
          !---------------------------------------------------------------------
          !
          ZMOLECO2   = ZMOLECO2V(JL)
          zmolecco2  = ZMOLECCO2V(JL)
          !
          !
          !----------------------------------------------------------------------
          !
          !
          !----------------------------------------------------------------------
          !
          !     T0-WERTE SETZEN
          !
          ZMOLECOXT0         = ZMOLECOX
          ZMOLECHCLT0        = ZMOLECHCL
          ZMOLECCLOXT0       = ZMOLECCLOX
          ZMOLECCLNO3T0      = ZMOLECCLNO3
          ZMOLECCLOHT0       = ZMOLECCLOH
          ZMOLECCLOT0        = ZMOLECCLO
          ZMOLECCL2O2T0      = ZMOLECCL2O2
          !
          ZMOLECNOXT0        = ZMOLECNOX
          ZMOLECN2O5T0       = ZMOLECN2O5
          ZMOLECHNO4T0       = ZMOLECHNO4
          ZMOLECNO2T0        = ZMOLECNO2
          !
          ZMOLECCH3O2T0      = ZMOLECCH3O2
          ZMOLECHO2T0        = ZMOLECHO2
          !
          !
          !     TEMPERATURE INDEX JIND OF LOOKUP-TABLE
          !     FOR GAS-PHASE REACTION RATES
          !
          jind=JINDV(JL)
          !
          !     CALCULATION OF RATE-CONSTANTS, NAMES, SEE CHR.BRUEHL
          !
          ! Table of chemical reactions
          !
          ! *********************************************************************
          ! *   TEMPERATURE AND PRESSURE INDEPENDENT RATE CONSTANTS             *
          ! *                                                                   *
          ! *              ZCO2   :  O(1D) + O3     ----> 2O2                   *
          ! *              ZCN2   :  N + NO         ----> N2 + O(3P)            *
          ! *              ZCN3   :  N + OH         ----> NO + H                *
          ! *              ZCN4   :  N + HO2        ----> NO + OH               *
          ! *              ZCN21  :  N2O + O(1D)    ----> 2NO                   *
          ! *              ZCN22  :  N2O + O(1D)    ----> N2 + O2               *
          ! *              ZCN23  :  N + NO2        ----> N2 + O2               *
          ! *              ZCC12  :  CH3Cl + O(1D)  ----> CH2Cl + OH            *
          ! *              ZCC13  :  MCF + O(1D)    ----> CH2CCl3 + OH          *
          ! *              ZCC16  :  CF2Cl2 + O(1D) ----> CF2Cl + ClO           *
          ! *              ZCC17  :  CFCl3 + O(1D)  ----> CFCl2 + ClO           *
          ! *              ZCC18  :  CCl4 + O(1D)   ----> CCl3 + ClO            *
          ! *              ZCH3   :  H2 + O(1D)     ----> H + OH                *
          ! *              ZCH7   :  CH4 + O(1D) (+O2)--> CH3O2 + OH            *
          ! *              ZCH8   :  H2O + O(1D)    ----> 2OH                   *
          ! *              ZCH10  :  H + HO2        ----> 2OH                   *
          ! *              ZCH16  :  CH2O + OH      ----> CO + H2O + HO2        *
          ! *              ZCH17  :  H + HO2        ----> H2 + O2               *
          ! *              ZCH18  :  H + HO2        ----> H2O +O(3P)            *
          ! *              ZCH20  :  CH3O2H +OH     ----> CH3O2 + H2O           *
          ! *              ZCH21  :  CH3O2H +OH     ----> CH2O + H2O +OH        *
          ! *                                                                   *
          ! *********************************************************************
          ! * other rate constants:                                     *
          ! *                                                           *
          ! *   ZCO1      :  O(3P) + O3     ----> 2 O2                  *
          ! *   ZCO4      :  O(1D) + O2     ----> O(3P) + O2            *
          ! *   ZCO5      :  O(1D) + N2     ----> O(3P) + N2            *
          ! *   ZCN1      :  N + O2         ----> NO + O(3P)            *
          ! *   ZCN6      :  HNO4 + OH      ----> NO2 + H2O + O2        *
          ! *   ZCN8      :  NO + O3        ----> NO2 + O2              *
          ! *   ZCN9      :  NO + ClO       ----> NO2 + Cl              *
          ! *   ZCN10     :  NO + HO2       ----> NO2 + OH              *
          ! *   ZCN11     :  NO + CH3O2     ----> NO2 + CH3O            *
          ! *   ZCN12     :  NO + NO3       ----> 2 NO2                 *
          ! *   ZCN13     :  NO2 + O(3P)    ----> NO + O2               *
          ! *   ZCN15     :  NO2 + O3       ----> NO3 + O2              *
          ! *   ZCN20     :  ClONO2 +O(3P)  ----> ClO + NO3             *
          ! *   ZCC1      :  ClO + O(3P)    ----> Cl + O2               *
          ! *   ZCC2      :  HCl + OH       ----> Cl + H2O              *
          ! *   ZCC3      :  ClO + OH       ----> Cl + HO2              *
          ! *   ZCC3B     :  ClO + OH       ----> HCl + O2              *
          ! *   ZCC4      :  Cl + O3        ----> ClO + O2              *
          ! *   ZCC5      :  Cl + H2        ----> HCl + H               *
          ! *   ZCC6      :  Cl + CH4 (+O2) ----> HCl + CH3O2           *
          ! *   ZCC7      :  Cl + HO2       ----> HCl + O2              *
          ! *   ZCC8      :  Cl + HO2       ----> ClO + OH              *
          ! *   ZCC9      :  Cl + HCHO (+O2)----> HCl + CO + HO2        *
          ! *   ZCC10     :  ClO + HO2      ----> HOCl + O2             *
          ! *   ZCC11     :  HOCl + OH      ----> ClO + H2O             *
          ! *   ZCC14     :  CH3Cl + OH     ----> CH2Cl + H2O           *
          ! *   ZCC15     :  CH3CCl3 + OH   ----> CH2CCl3 + H2O         *
          ! *   ZCH4      :  OH + O(3P)     ----> O2 + H                *
          ! *   ZCH6      :  H2 + OH        ----> H2O + H               *
          ! *   ZCH9      :  O3 + H         ----> OH + O2               *
          ! *   ZCH11     :  O(3P) + HO2    ----> OH + O2               *
          ! *   ZCH12     :  O3 + HO2       ----> OH + 2 O2             *
          ! *   ZCH14     :  O3 + OH        ----> HO2 + O2              *
          ! *   ZCH15     :  H2O2 + OH      ----> H2O + HO2             *
          ! *   ZCH19     :  CH4 + OH (+O2) ----> CH3O2 + H2O           *
          ! *   ZCH22     :  CH3O2 + HO2    ----> CH3O2H + O2           *
          ! *   ZCH23     :  O2 + CH3O      ----> HCHO + HO2            *
          ! *   ZCO3      :  O(3P) + O2 (+M)----> O3                    *
          ! *   ZCN14     :  NO2 + ClO (+M) ----> ClONO2                *
          ! *   ZCN16     :  NO2 + HO2 (+M) ----> HNO4                  *
          ! *   ZCN17     :  NO2 + OH (+M)  ----> HNO3                  *
          ! *   ZCN18     :  NO2 + NO3 (+M) ----> N2O5                  *
          ! *   ZCH13     :  H + O2 (+M)    ----> HO2                   *
          ! *   ZCCLOCLO  :  2 ClO (+M)     ----> Cl2O2                 *
          ! *   ZCN19     :  HNO3 + OH      ----> H2O + NO3             *
          ! *   ZCH1      :  OH + HO2       ----> H2O + O2              *
          ! *   ZCH2      :  2 HO2          ----> H2O2 + O2             *
          ! *   ZCH5      :  CO + OH        ----> CO2 + H               *
          ! *                                                           *
          ! * equilibrium constants:
          ! *   ZCN5      : HNO4   / ZCN16
          ! *   ZCN7      : N2O5   / ZCN18
          ! *   ZCCL2O2   : Cl2O2  / ZCCLOCLO
          ! *
          ! *********************************************************************
          !     ARRHENIUS-TYPE
          !
          !      ZCO1        = 8.0E-12*EXPHF(-2060./ZTEMP)
          ZCO1        = RCGAS(JIND,1)

          !      ZCO4        = 3.2E-11*EXPHF(70./ZTEMP)
          ZCO4        = RCGAS(JIND,2)

          !      ZCO5        = 1.8E-11*EXPHF(110./ZTEMP)
          ZCO5        = RCGAS(JIND,3)

          !
          !     ZCN1        = 1.5E-11*EXPHF(-3600./ZTEMP)
          ZCN1        = RCGAS(JIND,4)

          !      ZCN6        = 1.3E-12*EXPHF(380./ZTEMP)
          ZCN6        = RCGAS(JIND,5)

          !      ZCN8        = 2.0E-12*EXPHF(-1400./ZTEMP)
          ZCN8        = RCGAS(JIND,6)

          !      ZCN9        = 6.4E-12*EXPHF(290./ZTEMP)
          ZCN9        = RCGAS(JIND,7)

          !      ZCN10       = 3.7E-12*EXPHF(250./ZTEMP)
          ZCN10       = RCGAS(JIND,8)

          !      ZCN11       = 4.2E-12*EXPHF(180./ZTEMP)
          ZCN11       = RCGAS(JIND,9)

          !      ZCN12       = 1.5E-11*EXPHF(170./ZTEMP)
          ZCN12       = RCGAS(JIND,10)

          !      ZCN13       = 6.5E-12*EXPHF(120./ZTEMP)
          ZCN13       = RCGAS(JIND,11)

          !      ZCN15       = 1.2E-13*EXPHF(-2450./ZTEMP)
          ZCN15       = RCGAS(JIND,12)

          !      ZCN20       = 2.9E-12*EXPHF(-800./ZTEMP)
          ZCN20       = RCGAS(JIND,13)

          !
          !      ZCC1        = 3.0E-11*EXPHF(70./ZTEMP)
          ZCC1        = RCGAS(JIND,14)

          !      ZCC2        = 2.6E-12*EXPHF(-350./ZTEMP)
          ZCC2        = RCGAS(JIND,15)

          !      ZCC3        = 1.1E-11*EXPHF(120./ZTEMP)
          ZCC3        = 0.95*RCGAS(JIND,16)
          ZCC3B       = 0.05*RCGAS(JIND,16)

          !      ZCC4        = 2.9E-11*EXPHF(-260./ZTEMP)
          ZCC4        = RCGAS(JIND,17)

          !      ZCC5        = 3.7E-11*EXPHF(-2300./ZTEMP)
          ZCC5        = RCGAS(JIND,18)

          !      ZCC6        = 1.1E-11*EXPHF(-1400./ZTEMP)
          ZCC6        = RCGAS(JIND,19)

          !      ZCC7        = 1.8E-11*EXPHF(170./ZTEMP)
          ZCC7        = RCGAS(JIND,20)

          !      ZCC8        = 4.1E-11*EXPHF(-450./ZTEMP)
          ZCC8        = RCGAS(JIND,21)

          !      ZCC9        = 8.1E-11*EXPHF(-30./ZTEMP)
          ZCC9        = RCGAS(JIND,22)

          !      ZCC10       = 4.8E-13*EXPHF(700./ZTEMP)
          ZCC10       = RCGAS(JIND,23)

          !      ZCC11       = 3.0E-12*EXPHF(-500./ZTEMP)
          ZCC11       = RCGAS(JIND,24)

          !      ZCC14       = 2.1E-12*EXPHF(-1150./ZTEMP)
          ZCC14       = RCGAS(JIND,25)

          !      ZCC15       = 1.8E-12*EXPHF(-1550./ZTEMP)
          ZCC15       = RCGAS(JIND,26)
          !
          !     HELLEIS,CROWLEY,MOORTGAT MPI-MAINZ
          !
          !      ZCCH3O2CLO  = 3.25E-12*EXPHF(-114./ZTEMP)
          !
          !      ZCH4        = 2.2E-11*EXPHF(120./ZTEMP)
          ZCH4        = RCGAS(JIND,27)

          !      ZCH6        = 5.5E-12*EXPHF(-2000./ZTEMP)
          ZCH6        = RCGAS(JIND,28)

          !      ZCH9        = 1.4E-10*EXPHF(-470./ZTEMP)
          ZCH9        = RCGAS(JIND,29)

          !      ZCH11       = 3.0E-11*EXPHF(200./ZTEMP)
          ZCH11       = RCGAS(JIND,30)

          !      ZCH12       = 1.1E-14*EXPHF(-500./ZTEMP)
          ZCH12       = RCGAS(JIND,31)

          !      ZCH14       = 1.6E-12*EXPHF(-940./ZTEMP)
          ZCH14       = RCGAS(JIND,32)

          !      ZCH15       = 2.9E-12*EXPHF(-160./ZTEMP)
          ZCH15       = RCGAS(JIND,33)

          !      ZCH19       = 2.95E-12*EXPHF(-1815./ZTEMP)
          ZCH19       = RCGAS(JIND,34)

          !      ZCH22       = 3.8E-13*EXPHF(800./ZTEMP)
          ZCH22       = RCGAS(JIND,35)

          !      ZCH23       = 3.9E-14*EXPHF(-900./ZTEMP)
          ZCH23       = RCGAS(JIND,36)
          !
          !     RATE CONSTANTS FOR THREE-BODY REACTIONS
          !
          !     DEFINITION OF ZZT,ZK1,ZK2 IN DEMORE ET AL. (1992)
          !
          !      ZZT         = ZTEMP/300.
          !      ZZTLN       = LOG(ZZT)
          !
          !      ZCO3        = 6.0E-34/EXPHF(ZZTLN*2.3)*ZCON
          ZCO3        = RCGAS(JIND,37)*ZCON

          !      ZK1         = 1.8E-31/EXPHF(ZZTLN*3.4)*ZCON
          ZK1         = RCGAS(JIND,38)*ZCON

          !      ZK2         = 1.5E-11/EXPHF(ZZTLN*1.9)
          ZK2         = RCGAS(JIND,39)
          ZCN14       = ZFUNC3BODY(ZK1,ZK2)

          !      ZK1         = 1.8E-31/EXPHF(ZZTLN*3.2)*ZCON
          ZK1         = RCGAS(JIND,40)*ZCON

          !      ZK2         = 4.7E-12/EXPHF(ZZTLN*1.4)
          ZK2         = RCGAS(JIND,41)
          ZCN16       = ZFUNC3BODY(ZK1,ZK2)

          !      ZK1         = 2.6E-30/EXPHF(ZZTLN*3.2)*ZCON
          ZK1         = RCGAS(JIND,42)*ZCON

          !      ZK2         = 2.4E-11/EXPHF(ZZTLN*1.3)
          ZK2         = RCGAS(JIND,43)
          ZCN17       = ZFUNC3BODY(ZK1,ZK2)

          !      ZK1         = 2.2E-30/EXPHF(ZZTLN*3.9)*ZCON
          ZK1         = RCGAS(JIND,44)*ZCON

          !      ZK2         = 1.5E-12/EXPHF(ZZTLN*0.7)
          ZK2         = RCGAS(JIND,45)
          ZCN18       = ZFUNC3BODY(ZK1,ZK2)

          !      ZK1         = 5.7E-32/EXPHF(ZZTLN*1.6)*ZCON
          ZK1         = RCGAS(JIND,46)*ZCON
          ZK2         = RCGAS(JIND,66)*ZCON 
          ZCH13       = ZFUNC3BODY(ZK1,ZK2)

          !      ZK1         = 1.9E-32/EXPHF(ZZTLN*3.9)*ZCON
          ZK1         = RCGAS(JIND,47)*ZCON
          ZK2         = RCGAS(JIND,61)
          ZCCLOCLO    = ZFUNC3BODY(ZK1,ZK2)
          !S1      ZK1         = 1.5E-30/ZZT**4*ZCON
          !S1      ZK2         = 6.5E-12/ZZT**2
          !S1      ZCNO2CH3O2  = ZFUNC3BODY(ZK1,ZK2)
          !
          !     THERMAL DECOMPOSITION RATE CONSTANTS ; EQUILIBRIUM CONSTANTS
          !     FROM DEMORE ET AL. (1997)
          !     (IN DEMORE GILT: HINREAKT./EQUIKON.=REAKTIONSRATE, IN DIE
          !      M SCHON EINMULTIPLIZIERT IST !)

          !
          !
          !      ZCN5        = ZCN16/(2.1E-27*EXPHF(10900./ZTEMP))
          ZCN5        = ZCN16/RCGAS(JIND,48)

          !      ZCN7        = ZCN18/(4.0E-27*EXPHF(10930./ZTEMP))
          ZCN7        = ZCN18/RCGAS(JIND,49)

          !      ZCEQ        = 3.0E-27*EXPHF(8450./ZTEMP)
          ZCEQ        = RCGAS(JIND,50)
          ZCCL2O2M    = ZCCLOCLO/ZCEQ
          !S1      ZCCH3O2NO2M = ZCNO2CH3O2/(1.3E-28*EXPHF(11200./ZTEMP))
          !
          !     BINARY REACTIONS WITH SPECIAL CALCULATIONS; DEMORE ET AL. (1992)
          !
          !      ZK0         = 7.2E-15*EXPHF(785./ZTEMP)
          ZK0         = RCGAS(JIND,51)

          !      ZK3         = 1.9E-33*EXPHF(725./ZTEMP)*ZCON
          ZK3         = RCGAS(JIND,52)*ZCON

          !      ZK2         = 4.1E-16*EXPHF(1440./ZTEMP)
          ZK2         = RCGAS(JIND,53)

          ZCN19       = ZK0+ZK3/(1._dp+ZK3/ZK2)
          ! Try whether this is faster:
          !      ZCN19       = ZK0+ZK2*ZK3/(ZK2+ZK3)
          ! rh, 29-sep-1997

          !      ZCH1        = 4.8E-11*EXPHF(250./ZTEMP)
          ZCH1        = RCGAS(JIND,54)

          !      ZCH2        = 2.3E-13*EXPHF(600./ZTEMP)+ &
          !     *              1.7E-33*ZCON*EXPHF(1000./ZTEMP)
          ZCH2        = (RCGAS(JIND,55)+ &
               RCGAS(JIND,56)*ZCON)* &
               (1._dp+RCGAS(JIND,62)*ZMOLECH2O)
          !
          ! Use 101325. instead of 101400. rh, 29-sep-1997
          !
          ZCH5        = 1.5E-13*(1._dp+0.6_dp*ZPRES/101400._dp)
          !
          ZCN2      = RCGAS(JIND,57)
          ZCN23     = RCGAS(JIND,58)
          ZCH20     = RCGAS(JIND,59)
          ZCH21     = RCGAS(JIND,60)
          !
          ! mz_ab_20091221+
          ZCN21 = RCGAS(JIND,63)
          ZCN22 = RCGAS(JIND,64)
          ZCH8  = RCGAS(JIND,65)
          ZCH16 = RCGAS(JIND,66)
          ! mz_ab_20091221-
          !
          !     DAY NIGHT BRANCHING
          !
          !     DANI = 1.1  DAYTIME
          !     DANI = 0.1  FIRST TIMESTEP IN NIGHTTIME
          !     DANI = -1.1 NIGHTTIME
          !
!!$          ifdani: IF(DANI(JL).GT.1._dp) THEN ! op_bk_20121219 moved to below
             ! ******************************************************************
             ! *                                                                *
             ! * daytime chemistry                                              *
             ! *                                                                *
             ! ******************************************************************
             !
             !     SET PHOTOLYSES-RATES!
             !
             ZDO2       = RJ_O2(JL,JK)
             ZDO3B      = RJ_O3P(JL,JK)
             ZDH2O2     = RJ_H2O2(JL,JK)
             ZDH2O      = RJ_H2O(JL,JK)
             ZDHNO3     = RJ_HNO3(JL,JK)
             ZDN2O      = RJ_N2O(JL,JK)
             ZDNO2      = RJ_NO2(JL,JK)
             ZDN2O5     = RJ_N2O5(JL,JK)
             ZDCH2OA    = RJ_COH2(JL,JK)
             ZDCH2OB    = RJ_CHOH(JL,JK)
             ZDCLNO3    = RJ_CLONO2(JL,JK)
             ZDCH3CL    = RJ_CH3CL(JL,JK)
             ZDCCL4     = RJ_CCL4(JL,JK)
             ZDCFCL3    = RJ_CFC11(JL,JK)
             ZDCF2CL2   = RJ_CFC12(JL,JK)
             ZDCH3CCL3  = RJ_CH3CCL3(JL,JK)
             ZDNO       = RJ_NO(JL,JK)
             ZDO3A      = RJ_O1D(JL,JK)
             ZDCO2      = RJ_CO2(JL,JK)
             ZDHCL      = RJ_HCL(JL,JK)
             ZDCLOH     = RJ_HOCL(JL,JK)
             ZDCH3O2H   = RJ_CH3OOH(JL,JK)
             ZDHNO4     = RJ_HNO4(JL,JK)
             ZDNO3A     = RJ_NO2O(JL,JK)
             ZDNO3B     = RJ_NOO2(JL,JK)
             ZDCL2O2A   = RJ_CL2O2(JL,JK)
             !
             !-------------------------------------------------------------
         ifdani: IF(DANI(JL).GT.1._dp) THEN ! op_bk_20121219 moved from above
             !
             ! Perform chemical integration now:
             !
             ZMOLECCL2  = 0._dp
             !
             !     FIRST GUESS FOR [O],[O1D],[O3]
             !
             !
             !     FIRST ESTIMATE OF [OX]=[O1D]+[O]+[O3];
             !     CALCULATED FROM STEADY-STATE ASSUMPTIONS FOR [O1D] AND [O]
             !
             ZMOLECO1D   = ZMOLECO1DV(JL)
             ZMOLECO     = ZMOLECOV(JL)
             !
             !
             !     PRODUCTION OF CLOX BY HALOCARBONCHEMISTRY
             !     MASS IS STRICTLY CONSERVED!
             !     LOOK FOR EXAMPLE JUST AT F11:
             !     THE FOLLOWING EXPRESSION FOR THE CHANGE OF CLOY
             !     CLOY=CLOYt0+3*DT*(ZDCFCL3+ZCC17*ZMOLECO1D)*F11
             !     IS EQUAL TO CLOY=CLOYt0+3(F11t0-F11)
             !
             !
             ZPRCFCL3     = ZDCFCL3+ZCC17*ZMOLECO1D
             ZPRCFCL3     = ZPRCFCL3*ZTMSTDT
             IF(ZPRCFCL3.LT.1.E-10_dp) THEN
                ZPRCFCL3     = 1._dp
             ELSE
                ZPRCFCL3     = EXP(-1._dp*ZPRCFCL3)
             ENDIF
             !
             ZPRCF2CL2    = ZDCF2CL2+ZCC16*ZMOLECO1D
             ZPRCF2CL2    = ZPRCF2CL2*ZTMSTDT
             IF(ZPRCF2CL2.LT.1.E-10_dp) THEN
                ZPRCF2CL2     = 1._dp
             ELSE
                ZPRCF2CL2     = EXP(-1._dp*ZPRCF2CL2)
             ENDIF
             !
             ZPRCCL4      = ZDCCL4+ZCC18*ZMOLECO1D
             ZPRCCL4      = ZPRCCL4*ZTMSTDT
             IF(ZPRCCL4.LT.1.E-10_dp) THEN
                ZPRCCL4     = 1._dp
             ELSE
                ZPRCCL4     = EXP(-1._dp*ZPRCCL4)
             ENDIF
             !
             ZPRCH3CL     = ZDCH3CL+ZCC12*ZMOLECO1D+ZCC14*ZMOLECOH
             ZPRCH3CL     = ZPRCH3CL*ZTMSTDT
             IF(ZPRCH3CL.LT.1.E-10_dp) THEN
                ZPRCH3CL     = 1._dp
             ELSE
                ZPRCH3CL     = EXP(-1._dp*ZPRCH3CL)
             ENDIF
             !
             ZPRCH3CCL3   = ZDCH3CCL3+ZCC13*ZMOLECO1D+ZCC15*ZMOLECOH
             ZPRCH3CCL3   = ZPRCH3CCL3*ZTMSTDT
             IF(ZPRCH3CCL3.LT.1.E-10_dp) THEN
                ZPRCH3CCL3     = 1._dp
             ELSE
                ZPRCH3CCL3     = EXP(-1._dp*ZPRCH3CCL3)
             ENDIF
             !
             ZPRODCLOY    = 3._dp*ZMOLECCFCL3*(1._dp-ZPRCFCL3)+ &
                  2._dp*ZMOLECCF2CL2*(1._dp-ZPRCF2CL2)+ &
                  4._dp*ZMOLECCCL4*(1._dp-ZPRCCL4)+ &
                  3._dp*ZMOLECCH3CCL3*(1._dp-ZPRCH3CCL3)+ &
                  ZMOLECCH3CL*(1._dp-ZPRCH3CL)
             !
             ZMOLECCFCL3  = ZMOLECCFCL3*ZPRCFCL3
             ZMOLECCF2CL2 = ZMOLECCF2CL2*ZPRCF2CL2
             ZMOLECCCL4   = ZMOLECCCL4*ZPRCCL4
             ZMOLECCH3CL  = ZMOLECCH3CL*ZPRCH3CL
             ZMOLECCH3CCL3= ZMOLECCH3CCL3*ZPRCH3CCL3
             !
             !      ZPRODCLOY    = 0.
             !
             ZMOLECCLOY   = ZMOLECCLOX+ZMOLECHCL
             ZMOLECCLOY   = ZPRODCLOY+ZMOLECCLOY
             ZMOLECCLOX   = ZPRODCLOY+ZMOLECCLOX
             !
             ZRATNONO2    = ZFUNCRATNONO2(ZDNO2,ZCN13,ZMOLECO,ZCN8,ZMOLECO3, &
                  ZCN9,ZMOLECCLO,ZCN10,ZMOLECHO2,ZCN11,ZMOLECCH3O2)
             ZRATNNO      = ZFUNCRATNNO (ZDNO,ZCN1,ZMOLECO2,ZCN23,ZMOLECNO2, &
                  ZCN3,ZMOLECOH,ZCN4,ZMOLECHO2,ZCN2,ZMOLECNO)
             !
             !     SOLVE THE COMBINED SYSTEM OF NOX AND CLOX
             !
             !     A) CLOX
             !
             ZDECL        = ZCC4*ZMOLECO3+ZCC8*ZMOLECHO2
             ZDECL        = 1._dp/ZDECL
             ZQCLOH       = ZDCLOH+ZCC11*ZMOLECOH+ZCP3
             IF(ZQCLOH.LT.1.E-10_dp) THEN
                ZEQCLOH      = 1._dp
                ZAUXU1       = ZTMSTDT
             ELSE
                ZEQCLOH      = EXP(-1._dp*ZQCLOH*ZTMSTDT)
                ZAUXU1       = (1._dp-ZEQCLOH)/ZQCLOH
             ENDIF
             ZQCLNO3      = ZDCLNO3+ZCN20*ZMOLECO+ZHET2
             IF(ZQCLNO3.LT.1.E-10_dp) THEN
                ZEQCLNO3     = 1._dp
                ZAUXU4       = ZTMSTDT
             ELSE
                ZEQCLNO3     = EXP(-1._dp*ZQCLNO3*ZTMSTDT)
                ZAUXU4       = (1._dp-ZEQCLNO3)/ZQCLNO3
             ENDIF
             ZAUXU2       = 1._dp+(ZDCLOH+2._dp*ZCP3)*ZDECL
             ZAUXU3       = 1._dp+(ZDCLNO3+2._dp*ZCP1)*ZDECL
             ZAUXU10      = 1._dp+ZTMSTDT*(ZDCL2O2A+ZCCL2O2M)
             ZAUXU11      = 1._dp+ZDCL2O2A*ZDECL
             ZU           = ZFMOLECCLOX- &
                  ZAUXU2*ZMOLECCLOHT0*ZEQCLOH- &
                  (ZAUXU3+ZAUXU2*ZAUXU1*ZHET3)*ZMOLECCLNO3T0*ZEQCLNO3 &
                  -2._dp*ZMOLECCL2O2T0/ZAUXU10*ZAUXU11
             ZV           = 1._dp+(ZCC1*ZMOLECO+ZCC3*ZMOLECOH)*ZDECL+ &
                  ZAUXU2*ZAUXU1*ZCC10*0._dp+ &  
 !                ZAUXU2*ZAUXU1*ZCC10*ZMOLECHO2+ & ! mz_ab_20090929 comment this line to get ClOx right
                  2._dp*ZTMSTDT*ZCCLOCLO*ZMOLECCLOT0/ZAUXU10*ZAUXU11
             ZW           = (ZAUXU3+ZAUXU2*ZAUXU1*ZHET3)*ZCN14*ZAUXU4+ &
                  ZCN9*ZRATNONO2*ZDECL
             !      ZWV(JL) = ZW
             !      ZVV(JL) = ZV
             !
             !     PARACHUTE
             !
             IF(ZU.LE.0._dp) THEN
                ZU          = ZMOLECCLOT0*0.5*(ZV+ZW*ZMOLECNO2T0*0.5)
             ENDIF
             !      ZUV(JL) = ZU
             !
             !     B) NOX
             !
             ZQNO3A       = ZDNO3A+ZDNO3B+ZCN12*ZRATNONO2*ZMOLECNO2
             ZQNO3        = ZQNO3A+ZCN18*ZMOLECNO2
             ZQN2O5       = ZHET1+(ZCN7+ZDN2O5)*ZQNO3A/ZQNO3
             IF(ZQN2O5.LT.1.E-10_dp) THEN
                ZEQN2O5      = 1._dp
                ZAUX1NOX     = ZTMSTDT
             ELSE
                ZEQN2O5      = EXP(-1._dp*ZQN2O5*ZTMSTDT)
                ZEMEQN2O5    = 1._dp-ZEQN2O5
                ZAUX1NOX     = ZEMEQN2O5/ZQN2O5
             ENDIF
             ZQHNO4       = ZDHNO4+ZCN6*ZMOLECOH+ZCN5
             IF(ZQHNO4.LT.1.E-10_dp) THEN
                ZEQHNO4      = 1._dp
                ZAUX2NOX     = ZTMSTDT
             ELSE
                ZEQHNO4      = EXP(-1._dp*ZQHNO4*ZTMSTDT)
                ZEMEQHNO4    = 1._dp-ZEQHNO4
                ZAUX2NOX     = ZEMEQHNO4/ZQHNO4
             ENDIF
             ZAUXK1       = 2._dp*ZHET1*ZAUX1NOX
             ZAUXK2       = 2._dp+(ZCN7+ZDN2O5)/ZQNO3
             ZAUXK3       = ZHET2*ZAUXU4
             ZAUXK4       = 1._dp+(ZDCLNO3+ZCN20*ZMOLECO)/ZQNO3
             ZKNOX        = ZMOLECNOXT0-ZMOLECHNO4T0*ZEQHNO4- &
                  ZMOLECN2O5T0*(ZAUXK1+ZAUXK2*ZEQN2O5)- &
                  ZMOLECCLNO3T0*(ZAUXK3+ZAUXK4*ZEQCLNO3)
             !
             !      ZKNOXV(JL) = ZKNOX
             !
             !     PARACHUTE
             !
             IF(ZKNOX.GT.0._dp) THEN
                ZG1         = (ZAUXK2-2._dp*ZHET1/ZQN2O5)*ZAUX1NOX+ &
                     2._dp*ZTMSTDT*ZHET1/ZQN2O5
                ZG2         = (ZAUXK4-ZHET2/ZQCLNO3)*ZAUXU4+ZTMSTDT*ZHET2/ZQCLNO3
                ZAUXF       = (ZDCLNO3+ZCN20*ZMOLECO)*ZCN18/ZQNO3
                ZF1         = ZAUXF*ZMOLECCLNO3T0*ZEQCLNO3
                ZF2         = ZCN18*ZCN15*ZMOLECO3/ZQNO3
                ZF3         = ZAUXF*ZAUXU4*ZCN14
                ZANOX       = 1._dp+ZRATNONO2*(1._dp+ZRATNNO)+ &
                     ZAUX2NOX*ZCN16*ZMOLECHO2+ & 
                     ZCN15*ZMOLECO3/ZQNO3+ZG1*ZF1
                ZBNOX       = ZG1*ZF2
                ZCNOX       = ZG2*ZCN14
                ZDNOX       = ZG1*ZF3
                ZC0         = ZKNOX*ZV
                ZC1         = ZANOX*ZV+ZCNOX*ZU-ZKNOX*ZW
                ZC2         = ZANOX*ZW+ZBNOX*ZV+ZDNOX*ZU
                ZC3         = ZBNOX*ZW
                !
                !       ZG1V(JL) = ZG1
                !       ZG2V(JL) = ZG2
                !       ZF1V(JL) = ZF1
                !       ZF2V(JL) = ZF2
                !       ZF3V(JL) = ZF3
                !       ZANOXV(JL) = ZANOX
                !       ZBNOXV(JL) = ZBNOX
                !       ZCNOXV(JL) = ZCNOX
                !       ZDNOXV(JL) = ZDNOX
                !       ZC0V(JL) = ZC0
                !       ZC1V(JL) = ZC1
                !       ZC2V(JL) = ZC2
                !       ZC3V(JL) = ZC3
                !
                !     PARACHUTE
                !
                ZMAX        = ZC2*ZC2
                ZMAX1       = ZMAX-3._dp*ZC1*ZC3
                IF(ZMAX1.GT.ZMAX) THEN
                   ZMOLECNO2  = MAX((-1._dp*ZC2+SQRT(ZMAX1))/ZC3,ZMOLECNO2)
                ENDIF
                !
                !     SOLVE WITH NEWTON-ITERATION!
                !
                !       DO 110 J=1,5
                !       ZAUXNNEW1   = ZC3*ZMOLECNO2**2
                !       ZAUXNNEW2   = ZC2*ZMOLECNO2
                !       ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/
                !     @               (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !      1.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                ZMOLECNO2   = MAX(1._dp,ZMOLECNO2)
                !
                !      2.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !      3.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !      4.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !      5.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                ZMOLECNO2   = MAX(1._dp,ZMOLECNO2)
                !
                !      6.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !      7.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !      8.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !      9.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !      10.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                ZMOLECNO2   = MAX(1._dp,ZMOLECNO2)
                !      ELSE
                !     UEBERNEHME FIRST GUESS WERTE
                !     WOBEI IN FIRST GUESS WERTE SCHON SICHERUNGEN EINGEBAUT
                !     SEIN MUESSEN!
             ENDIF
             !
             !      CLO NOCH BEZUEGLICH ZRATNONO2 NACHITERIEREN?
             !      ZRATNONO2 NEU BERECHNET MIT NEUEM CLO!
             !      ZRATNNO   NEU BERECHNET MIT NEUEM NO2!
             !      QNO3,QN2O5 NEU BERECHNET MIT NEUEM NO NO2
             !
             ZMOLECCLO   = ZU/(ZV+ZW*ZMOLECNO2)
             !-----------------------------------------------------------------
             ZMOLECNO    = ZRATNONO2*ZMOLECNO2
             !
             ZPCLNO3     = ZCN14*ZMOLECNO2*ZMOLECCLO
             ZMOLECCLNO3 = ZMOLECCLNO3T0*ZEQCLNO3+ &
                  ZPCLNO3*ZAUXU4
             ZPCLOH      = ZCC10*ZMOLECHO2*ZMOLECCLO+ZHET3*ZMOLECCLNO3
             ZMOLECCLOH  = ZMOLECCLOHT0*ZEQCLOH+ &
                  ZAUXU1*ZPCLOH
             !
             !     VERSION FUER EUBA IN CL2O2
             !
             ZMOLECCL2O2 = (ZMOLECCL2O2T0+ &
                  ZTMSTDT*ZCCLOCLO*ZMOLECCLOT0*ZMOLECCLO)/ZAUXU10
             ZMOLECCL    = (ZMOLECCLO*(ZCC1*ZMOLECO+ZCC3*ZMOLECOH+ &
                  ZCN9*ZMOLECNO)+ &
                  ZMOLECCLNO3*(ZDCLNO3+2._dp*ZCP1)+ &
                  ZMOLECCLOH*(ZDCLOH+2._dp*ZCP3)+ &
                  2._dp*ZDCL2O2A*ZMOLECCL2O2)*ZDECL
             ZMOLECHNO4  = ZMOLECHNO4T0*ZEQHNO4 &
                  +ZAUX2NOX*ZMOLECNO2*ZMOLECHO2*ZCN16
             ZMOLECN     = ZRATNNO*ZMOLECNO
             ZPN2O5       = ZCN18/ZQNO3*ZMOLECNO2* &
                  (ZCN15*ZMOLECO3*ZMOLECNO2+ &
                  (ZDCLNO3+ZCN20*ZMOLECO)*ZMOLECCLNO3)
             ZMOLECN2O5   = ZMOLECN2O5T0*ZEQN2O5+ &
                  ZAUX1NOX*ZPN2O5
             ZMOLECNO3    = (ZCN15*ZMOLECO3*ZMOLECNO2+(ZDN2O5+ZCN7)*ZMOLECN2O5+ &
                  (ZDCLNO3+ZCN20*ZMOLECO)*ZMOLECCLNO3)/ZQNO3
             !
             !    INTEGRATION AND PARTITIONING OF [HOX]
             !
             ZMOLECHO2T0     = ZMOLECHO2
             ZSTO1DENOMH     = ZCH9*ZMOLECO3+ZCH13*ZMOLECO2
             ZSTO2DENOMH     = ZCH10+ZCH17+ZCH18
             ! T0 FUER HO2
             ZSTODENOMH      = ZSTO1DENOMH+ZSTO2DENOMH*ZMOLECHO2T0
             !
             ZSTOKH          = ZDH2O*ZMOLECH2O+ZDCH2OB*ZMOLECCH2O+ &
                  ZMOLECH2*(ZCC5*ZMOLECCL+ZCH3*ZMOLECO1D)+ &
                  ZDHCL*ZFMOLECHCL
             ZSTOKKH         = ZSTOKH/ZSTODENOMH
             ZSTOBH          = ZCH4*ZMOLECO+ZCH5*ZMOLECCO+ZCH6*ZMOLECH2
             ZSTOBBH         = ZSTOBH/ZSTODENOMH
             !
             ZSTOK1OH        = ZMOLECO1D*(ZCC12*ZMOLECCH3CL+ZCH3*ZMOLECH2+ &
                  ZCH7*ZMOLECCH4+2._dp*ZCH8*ZMOLECH2O)+ &
                  ZDH2O*ZMOLECH2O+2._dp*ZDH2O2*ZMOLECH2O2+ZDHNO3* &
                  ZFMOLECHNO3+ZDCLOH*ZMOLECCLOH+ &
                  ZDCH3O2H*ZMOLECCH3O2H
             ZSTOKOH         = ZSTOK1OH+ &  ! T0 FUER HO2
                  ZSTOKKH*(ZCH9*ZMOLECO3+ZCH10*ZMOLECHO2T0)
             ZSTOB1OH        = ZCH16*ZMOLECCH2O+ZCH4*ZMOLECO+ZCH5*ZMOLECCO+ &
                  ZCH6*ZMOLECH2+ZCH14*ZMOLECO3+ZCH15*ZMOLECH2O2+ &
                  ZCN19*ZFMOLECHNO3+ZCN17*ZMOLECNO2+ &
                  ZCH19*ZMOLECCH4+ZCN6*ZMOLECHNO4+ZCC2*ZFMOLECHCL+ &
                  ZCC11*ZMOLECCLOH+ZCC3*ZMOLECCLO+ &
                  ZCC14*ZMOLECCH3CL+ZCH20*ZMOLECCH3O2H
             ZSTOBOH         = ZSTOB1OH-    & ! T0 FUER HO2
                  ZSTOBBH*(ZCH9*ZMOLECO3+ZCH10*ZMOLECHO2T0)
             ZSTOCOH         = ZCH11*ZMOLECO+ZCH12*ZMOLECO3+ZCN10*ZMOLECNO+ &          ! T0 FUER H
                  ZCC8*ZMOLECCL+ZCH10*ZMOLECH
             ZSTOEOH         = ZCH1
             ! T0 FUER HO2
             ZSTODENOMOH     = ZSTOBOH+ZSTOEOH*ZMOLECHO2T0
             ZSTOKKOH        = ZSTOKOH/ZSTODENOMOH
             ZSTOCCOH        = ZSTOCOH/ZSTODENOMOH
             !
             ZSTOKHOX        = 2._dp*(ZDH2O*ZMOLECH2O+ZDCH2OB*ZMOLECCH2O+ &
                  ZDH2O2*ZMOLECH2O2+ZDCH3O2H*ZMOLECCH3O2H)+ &
                  ZDHNO3*ZFMOLECHNO3+ZDCLOH*ZMOLECCLOH+ &
                  ZDHNO4*ZMOLECHNO4+ZDHCL*ZFMOLECHCL+ &
                  ZMOLECO1D*(ZCC12*ZMOLECCH3CL+2._dp*ZCH3*ZMOLECH2+ &
                  ZCH7*ZMOLECCH4+2._dp*ZCH8*ZMOLECH2O)+ &
                  ZMOLECCL*(ZCC5*ZMOLECH2+ZCC9*ZMOLECCH2O)+ &
                  ZCN5*ZMOLECHNO4+ZCN11*ZMOLECNO*ZMOLECCH3O2
             ZSTOBHOX        = ZCN19*ZFMOLECHNO3+ZCN17*ZMOLECNO2+ &
                  ZCH19*ZMOLECCH4+ZCN6*ZMOLECHNO4+ &
                  ZCC2*ZFMOLECHCL+ZCC11*ZMOLECCLOH+ &
                  ZCC14*ZMOLECCH3CL+ZCH20*ZMOLECCH3O2H
             ZSTOCHOX        = ZCH22*ZMOLECCH3O2+ZCN16*ZMOLECNO2+ &
                  ZCC7*ZMOLECCL+ZCC10*ZMOLECCLO
             ZSTODHOX        = 2._dp*(ZCH17+ZCH18)
             ZSTOEHOX        = 2._dp*ZCH1
             ZSTOFHOX        = 2._dp*ZCH2
             ZSTOKKHOX       = ZSTOKHOX-ZSTOKKOH*ZSTOBHOX
             ZSTOCCHOX       = ZSTODHOX*(ZSTOKKH+ZSTOBBH*ZSTOKKOH)+ &
                  ZSTOCCOH*ZSTOBHOX+ &
                  ZSTOCHOX+ZSTOKKOH*ZSTOEHOX
             ZSTOFFHOX       = ZSTODHOX*ZSTOBBH*ZSTOCCOH+ &
                  ZSTOCCOH*ZSTOEHOX+ZSTOFHOX
             !
             !
             ZAUXHOX         = 1._dp+ZSTOCCOH*(1._dp+ZSTOBBH)
             !      ZAUXHOX         = ZAUXHOX/ZTMSTHDT
             ZAUXHOX         = ZAUXHOX/ZTMSTDT
             ZSTOPHOX        = (ZSTOCCHOX+ZAUXHOX)/(2._dp*ZSTOFFHOX)
             ! T0 FUER HO2
             ZSTOQHOX        = (ZMOLECHO2T0*ZAUXHOX+ZSTOKKHOX)/ZSTOFFHOX
             ZMOLECHO2       = -1._dp*ZSTOPHOX+SQRT(ZSTOPHOX*ZSTOPHOX+ZSTOQHOX)
             !
             !
             ZSTODENOMH      = ZSTO1DENOMH+ZSTO2DENOMH*ZMOLECHO2
             ZSTOKKH         = ZSTOKH/ZSTODENOMH
             ZSTOBBH         = ZSTOBH/ZSTODENOMH
             ZSTOKOH         = ZSTOK1OH+ &
                  ZSTOKKH*(ZCH9*ZMOLECO3+ZCH10*ZMOLECHO2)
             ZSTOBOH         = ZSTOB1OH- &
                  ZSTOBBH*(ZCH9*ZMOLECO3+ZCH10*ZMOLECHO2)
             !
             ZMOLECOH        = (ZSTOKOH+ZSTOCOH*ZMOLECHO2)/ &
                  (ZSTOBOH+ZSTOEOH*ZMOLECHO2)
             !
             !
             ZMOLECH         = (ZSTOKH+ZSTOBH*ZMOLECOH)/ZSTODENOMH
             !
             !----------------------------------------------------------------------
             !
             !     CALCULATION OF H2O2
             !
             ZQH2O2      = ZDH2O2+ZCH15*ZMOLECOH
             ZMOLECH2O2  = (ZMOLECH2O2+ZTMSTDT*ZCH2*ZMOLECHO2**2)/ &
                  (1._dp+ZTMSTDT*ZQH2O2)
             !
             !     AFTERBURNER CLOH
             !
             ZPCLOH      = ZCC10*ZMOLECHO2*ZMOLECCLO+ZHET3*ZMOLECCLNO3
             ZQCLOH       = ZDCLOH+ZCC11*ZMOLECOH+ZCP3
             IF(ZQCLOH.LT.1.E-10_dp) THEN
                ZEQCLOH      = 1._dp
                ZAUXU1       = ZTMSTDT
             ELSE
                ZEQCLOH      = EXP(-1._dp*ZQCLOH*ZTMSTDT)
                ZAUXU1       = (1._dp-ZEQCLOH)/ZQCLOH
             ENDIF
             ZMOLECCLOH  = ZMOLECCLOHT0*ZEQCLOH+ &
                  ZAUXU1*ZPCLOH
             !
             !     AFTERBURNER HNO4
             !
             ZQHNO4       = ZDHNO4+ZCN6*ZMOLECOH+ZCN5
             IF(ZQHNO4.LT.1.E-10_dp) THEN
                ZEQHNO4      = 1._dp
                ZAUX2NOX     = ZTMSTDT
             ELSE
                ZEQHNO4      = EXP(-1._dp*ZQHNO4*ZTMSTDT)
                ZEMEQHNO4    = 1._dp-ZEQHNO4
                ZAUX2NOX     = ZEMEQHNO4/ZQHNO4
             ENDIF
             ZMOLECHNO4  = ZMOLECHNO4T0*ZEQHNO4 &
                  +ZAUX2NOX*ZMOLECNO2*ZMOLECHO2*ZCN16
             !
             !     HETEROGENER VERLUST VON NOX
             !
             IF(ZHET1.GT.1.E-20_dp) THEN
                ZHETLOSS1 = (ZMOLECN2O5T0-ZMOLECN2O5+ZPN2O5*ZTMSTDT)/ZQN2O5
                ZHETLOSS2    = (ZMOLECCLNO3T0-ZMOLECCLNO3+ZPCLNO3*ZTMSTDT)/ &
                     ZQCLNO3
                ZHETLOSSNOXTEST  = 2._dp*ZHET1*ZHETLOSS1+ZHET2*ZHETLOSS2
                ZMOLECNOX    = MAX(2._dp,ZMOLECNOX)
                ZHETLOSSNOX  = MIN(ZMOLECNOX-1._dp,ZHETLOSSNOXTEST)
                SCALE        = ZHETLOSSNOX/ZHETLOSSNOXTEST
                ZHETLOSS1    = SCALE*ZHETLOSS1
                ZHETLOSS2    = SCALE*ZHETLOSS2
                ZMOLECNOX    = ZMOLECNOX-ZHETLOSSNOX
                ZMOLECHNO3   = 2._dp*ZCS2*ZHETLOSS1+ZCS1*ZHETLOSS2+ &
                     ZMOLECHNO3
                ZMOLECNAT    = 2._dp*ZCP4*ZHETLOSS1+(ZCP1+ZCP2)*ZHETLOSS2+ &
                     ZMOLECNAT
             ENDIF
             !
             !     HETEROGENER VERLUST VON HCL
             !
             IF((ZCP1+ZCP3).GT.1.E-20_dp) THEN
                ZHETLOSSHCL  = ZCP1/ZQCLNO3*(ZMOLECCLNO3T0-ZMOLECCLNO3+ &
                     ZPCLNO3*ZTMSTDT)+ &
                     ZCP3/ZQCLOH*(ZMOLECCLOHT0-ZMOLECCLOH+ &
                     ZPCLOH*ZTMSTDT)
                ZHETLOSSHCL  = MIN(ZHETLOSSHCL,ZMOLECHCLT0)
                ZMOLECHCL    = ZMOLECHCLT0-ZHETLOSSHCL
                ZMOLECCLOX   = ZMOLECCLOXT0+ZHETLOSSHCL
                !       IF(ZHETLOSSHCL.LT.0.) THEN
                !       IDIAGV(JL) = 1
                !       ZHETLOSSHCLV(JL,1) = ZCP1
                !       ZHETLOSSHCLV(JL,2) = ZQCLNO3
                !       ZHETLOSSHCLV(JL,3) = ZMOLECCLNO3T0
                !       ZHETLOSSHCLV(JL,4) = ZMOLECCLNO3
                !       ZHETLOSSHCLV(JL,5) = ZPCLNO3
                !       ZHETLOSSHCLV(JL,6) = ZCP3
                !       ZHETLOSSHCLV(JL,7) = ZQCLOH
                !       ZHETLOSSHCLV(JL,8) = ZMOLECCLOHT0
                !       ZHETLOSSHCLV(JL,9) = ZMOLECCLOH
                !       ZHETLOSSHCLV(JL,10) = ZPCLOH
                !       ENDIF
             ENDIF
             IF((ZMOLECHCL.LE.10._dp).AND.(INIHET0.EQ.1)) THEN
                ZMOLECCLOX    = ZMOLECCLOY - ZMOLECHCL
             ELSE
                ZP1HCL        = ZCC6*ZMOLECCH4+ZCC7*ZMOLECHO2+ &
                     ZCC9*ZMOLECCH2O+ZCC5*ZMOLECH2
                ZP2HCL        = ZCC3B*ZMOLECOH
                ZQ1HCL        = ZCC2*ZMOLECOH+ZDHCL
                ZQCLOX        = (ZMOLECCL*ZP1HCL+ZMOLECCLO*ZP2HCL)/ &
                     ZMOLECCLOX+ZQ1HCL
                IF(ZQCLOX.LT.1.E-10_dp) THEN
                   ZEQCLOX       = 1._dp
                   ZAUX          = ZTMSTDT
                ELSE
                   ZEQCLOX       = EXP(-1._dp*ZTMSTDT*ZQCLOX)
                   ZAUX          = (1._dp-ZEQCLOX)/ZQCLOX
                ENDIF
                ZMOLECCLOX    = ZMOLECCLOX*ZEQCLOX+ZAUX*ZMOLECCLOY*ZQ1HCL
                ZMOLECHCL     = ZMOLECCLOY-ZMOLECCLOX
             ENDIF
             !
             ZMOLECCLOX = MAX(1._dp,ZMOLECCLOX)
             ZMOLECHCL  = MAX(1._dp,ZMOLECHCL)
             ZSCY       = ZMOLECCLOY/(ZMOLECCLOX+ZMOLECHCL)
             ZMOLECCLOX = ZSCY*ZMOLECCLOX
             ZMOLECHCL  = ZSCY*ZMOLECHCL
             !
             !     SCALING OF CL-ATOMS
             !     (FIRST CL-ATOMS, BECAUSE OF CLNO3)
             !
             ZSTOSCALCLOX    = ZMOLECCLOX/(ZMOLECCL+ZMOLECCLO+ZMOLECCLOH+ &
                  ZMOLECCLNO3+2._dp*ZMOLECCL2O2)
             ! 
             ZMOLECCL        = ZMOLECCL*ZSTOSCALCLOX
             ZMOLECCLO       = ZMOLECCLO*ZSTOSCALCLOX
             ZMOLECCLOH      = ZMOLECCLOH*ZSTOSCALCLOX
             ZMOLECCLNO3     = ZMOLECCLNO3*ZSTOSCALCLOX
             ZMOLECCL2O2     = ZMOLECCL2O2*ZSTOSCALCLOX
             !
             !     INTEGRATION OF N2O
             !
             ZQN2O       = ZMOLECO1D*(ZCN21+ZCN22)+ZDN2O
             ZMOLECN2O   = ZMOLECN2O/(1._dp+ZTMSTDT*ZQN2O)
             !
             !     INTEGRATION OF NOY
             !
             ZMOLECNOYT0 = ZMOLECNOX+ZMOLECHNO3
             ZQ1NOY      = MAX(ZMOLECN,1._dp)*2._dp* &
                  (ZCN2*ZMOLECNO+ZCN23*ZMOLECNO2)
             ZQNOY       = 1._dp+ZQ1NOY*ZTMSTDT/ZMOLECNOYT0
             ZMOLECNOY   = (ZMOLECNOYT0+2._dp*ZTMSTDT*ZCN21*ZMOLECN2O*ZMOLECO1D) &
                  /ZQNOY
             !
             !     HOMOGENEOUS INTEGRATION OF NOX AND HNO3
             !
             ZP1HNO3       = ZMOLECNO2/ZMOLECNOX*ZMOLECOH*ZCN17
             ZP1NOX        = ZDHNO3+ZCN19*ZMOLECOH
             ZP2NOX        = 2._dp*ZCN21*ZMOLECO1D*ZMOLECN2O
             ZQ1NOX        = ZQ1NOY/ZMOLECNOX
             ZQNOX         = ZQ1NOX+ZP1HNO3
             ZMOLECHNO3    = (ZMOLECHNO3+ZTMSTDT*ZMOLECNOYT0*ZP1HNO3)/ &
                  (1._dp+ZTMSTDT*(ZP1HNO3+ZP1NOX))
             ZMOLECNOX     = (ZMOLECNOX+ZTMSTDT*(ZMOLECNOYT0*ZP1NOX+ZP2NOX))/ &
                  (1._dp+ZTMSTDT*(ZQNOX+ZP1NOX))
             !
             !     MASS-CONSERVATION OF N-ATOMS
             !
             ZSCALATOMSN  = ZMOLECNOY/ &
                  (ZMOLECNOX+ZMOLECHNO3)
             ZMOLECHNO3   = ZMOLECHNO3*ZSCALATOMSN
             ZMOLECNOX    = ZMOLECNOX*ZSCALATOMSN
             !
             !     MASS-CONSERVATION FOR NOX AND CLOX FAMILY-MEMBERS, &
             !     WHEN CLNO3 TOO BIG !
             !
             !      ZSWMBTOTNOXV2(JL) = 0.
             IF ((ZMOLECNOX-ZMOLECCLNO3).LE.0.) THEN
                ZSTOSCALTOT =(ZMOLECNOX+ZMOLECCLOX)/ &
                     (ZMOLECN+ZMOLECNO+ZMOLECNO2+ZMOLECNO3+ &
                     ZMOLECHNO4+2._dp*ZMOLECN2O5+ZMOLECCL+ZMOLECCLO+ &
                     ZMOLECCLOH+2._dp*ZMOLECCLNO3+2._dp*ZMOLECCL2O2)
                ZMOLECCL    = ZMOLECCL*ZSTOSCALTOT
                ZMOLECCLO   = ZMOLECCLO*ZSTOSCALTOT
                ZMOLECCLOH  = ZMOLECCLOH*ZSTOSCALTOT
                ZMOLECCLNO3 = ZMOLECCLNO3*ZSTOSCALTOT
                ZMOLECCL2O2 = ZMOLECCL2O2*ZSTOSCALTOT
                ZMOLECN     = ZMOLECN*ZSTOSCALTOT
                ZMOLECNO    = ZMOLECNO*ZSTOSCALTOT
                ZMOLECNO2   = ZMOLECNO2*ZSTOSCALTOT
                ZMOLECNO3   = ZMOLECNO3*ZSTOSCALTOT
                ZMOLECHNO4  = ZMOLECHNO4*ZSTOSCALTOT
                ZMOLECN2O5  = ZMOLECN2O5*ZSTOSCALTOT
                !       ZSWMBTOTNOXV2(JL) = 1.
                !
                ! NOCH STRENGERE MASSENBILANZ CLOY !
                !
                ZSTOSCALCLOX    = ZMOLECCLOX/(ZMOLECCL+ZMOLECCLO+ZMOLECCLOH+ &
                     ZMOLECCLNO3+2._dp*ZMOLECCL2O2)
                ZMOLECCL        = ZMOLECCL*ZSTOSCALCLOX
                ZMOLECCLO       = ZMOLECCLO*ZSTOSCALCLOX
                ZMOLECCLOH      = ZMOLECCLOH*ZSTOSCALCLOX
                ZMOLECCLNO3     = ZMOLECCLNO3*ZSTOSCALCLOX
                ZMOLECCL2O2     = ZMOLECCL2O2*ZSTOSCALCLOX
             ELSE
                !
                !     MASS-CONSERVATION OF NOX-FAMILY-MEMBERS
                !
                ZSTOSCALNOX    = (ZMOLECNOX-ZMOLECCLNO3)/ &
                     (ZMOLECN+ZMOLECNO+ZMOLECNO2+ZMOLECNO3+ &
                     ZMOLECHNO4+2._dp*ZMOLECN2O5)
                ZMOLECN         = ZMOLECN*ZSTOSCALNOX
                ZMOLECNO        = ZMOLECNO*ZSTOSCALNOX
                ZMOLECNO2       = ZMOLECNO2*ZSTOSCALNOX
                ZMOLECNO3       = ZMOLECNO3*ZSTOSCALNOX
                ZMOLECHNO4      = ZMOLECHNO4*ZSTOSCALNOX
                ZMOLECN2O5      = ZMOLECN2O5*ZSTOSCALNOX
             ENDIF
             !-----------------------------------------------------------------
             !
             !     CALCULATION OF [CH4]
             !
             ZQCH4        = ZCH19*ZMOLECOH+ZCH7*ZMOLECO1D+ZCC6*ZMOLECCL
             ZMOLECCH4    = ZMOLECCH4/(1._dp+ZTMSTDT*ZQCH4)
             !
             !
             !     EUBA-VERSION VON CH3O2
             !
             ZSTODENOMCH3O2  = ZCN11*ZMOLECNO+ZCH22*ZMOLECHO2
             ZSTONUMERCH3O2  = ZMOLECCH4*ZQCH4+ZDCH3CL*ZMOLECCH3CL
             ZSTOPRODCH22HO2 = ZCH22*ZMOLECHO2
             ZSTOPRODCH23O2  = ZCH23*ZMOLECO2
             !
             !     INTEGRATION OF [CH3O2H]
             !
             ZQCH3O2H     = ZDCH3O2H+ZMOLECOH*(ZCH21+ &
                  ZCH20*(1._dp+ZTMSTDT*ZCN11*ZMOLECNO)/ &
                  (1._dp+ZTMSTDT*ZSTODENOMCH3O2))
             ZPCH3O2H     = ZSTOPRODCH22HO2*(ZMOLECCH3O2T0+ &
                  ZTMSTDT*ZSTONUMERCH3O2)/ &
                  (1._dp+ZTMSTDT*ZSTODENOMCH3O2)
             ZMOLECCH3O2H = (ZMOLECCH3O2H+ZTMSTDT*ZPCH3O2H)/ &
                  (1._dp+ZTMSTDT*ZQCH3O2H)
             !
             !     EUBA CALCULATION OF [CH3O2]
             !
             ZMOLECCH3O2  = (ZMOLECCH3O2T0+ZTMSTDT* &
                  (ZSTONUMERCH3O2+ZCH20*ZMOLECOH*ZMOLECCH3O2H))/ &
                  (1._dp+ZTMSTDT*ZSTODENOMCH3O2)
             !
             !     STEADY-STATE CALCULATION OF [CH3O]
             !
             ZMOLECCH3O   = (ZCN11*ZMOLECCH3O2*ZMOLECNO+ZDCH3O2H*ZMOLECCH3O2H)/ &
                  ZSTOPRODCH23O2
             !
             !     INTEGRATION OF [CH2O]
             !
             ZQCH2O       = ZDCH2OA+ZDCH2OB+ZCH16*ZMOLECOH+ZCC9*ZMOLECCL 
             ZPCH2O       = ZSTOPRODCH23O2*ZMOLECCH3O+ZCH21*ZMOLECOH* &
                  ZMOLECCH3O2H
             ZMOLECCH2O   = (ZMOLECCH2O+ZTMSTDT*ZPCH2O)/(1._dp+ZTMSTDT*ZQCH2O)
             !     INTEGRATION OF [CO]
             !
             ZMOLECCO     = ((ZQCH2O*ZMOLECCH2O+ZDCO2*ZMOLECCO2)*ZTMSTDT &
                  +ZMOLECCO)/ &
                  (1._dp+ZCH5*ZMOLECOH*ZTMSTDT)
             !
             !     CALCULATION OF [H2]
             !
             ZQH2     = ZCH3*ZMOLECO1D+ZCH6*ZMOLECOH
             ZPH2     = ZCH17*ZMOLECH*ZMOLECHO2+ZDCH2OA*ZMOLECCH2O
             ZMOLECH2 = (ZMOLECH2+ZTMSTDT*ZPH2)/(1._dp+ZTMSTDT*ZQH2)
             !
             ! zu_as_20100504+
             ! ***********************
             ! * Br parameterisation *
             ! ***********************
             ! Ozone depletion by bromine is proportional to photolysis of cl2o2
             ! 1. Calculation of scaling factor
             ! zbrpres: pressure in hPa
             ! zclox: volume mixing ratio (ClO + 2*Cl2O2) in ppbv
             !
             zclox = (zmolecclo + 2._dp*zmoleccl2o2)*zinvcon*1.E9
             zbrpres = zpres*0.01_dp
             !
             if (l_Brparam .AND. (zclox .ge. 0.25_dp) .and. (zbrpres.le.80._dp) .and. &
                  (zbrpres .ge. 20._dp)) then
                zbrtemp=ztemp+max(zbrpres-60._dp,0._dp)/2._dp
                zproxy=ABS(EXP(zbr1 - zbr2*LOG(zclox))+ &
                     EXP(zbr3 - zbr4*LOG(zclox))*(zbrtemp-195._dp)-&
                     EXP(zbr5 - zbr6*LOG(zclox))*(zbrpres-50._dp))
             else
                zproxy=0._dp
             endif
             !
             ! number of ozone molecules destroyed by BrO+ClO cycle
             ! proportinal to number of photolysed cl2o2 molecules
             ! zdeltao3_br only for budgets, integration below
             !
             ! mz_bk_20101222+
!              zdeltao3_br = zproxy*ZTMSTHDT*2._dp*ZDCL2O2A*ZMOLECCL2O2
              zdeltao3_br = zproxy*2._dp*ZDCL2O2A*ZMOLECCL2O2
             ! mz_bk_20101222-
             ! zu_as_20100504-
             !
             !     INTEGRATION OF [ODD]
             !     NEUER ANSATZ FUER DIE ODD-INTEGRATION
             !     HOMOGENE PRODUKTION UND HOMOGENE ZERSTOERUNG VON NICHT OX
             !     EXPLIZIT BEHANDELN.(BIS AUF H+HO2-->H2O+O)
             !     HETEROGENE ZERSTOERUNG UND REAKTIONEN MIT OX BETEILIGUNG
             !     IMPLIZIT!
             !
             ZSTO0        = ZMOLECCLNO3*(2._dp*ZCP1+ZCP2)+ &
                  ZMOLECCLOH*ZCP3+ &
                  ZMOLECN2O5*(3._dp*ZCP4+ZCS2) &
                  +ZDCLOH*ZMOLECCLOH+2._dp*ZMOLECCL2O2*ZDCL2O2A+ &
                  ZCC3*ZMOLECCLO*ZMOLECOH
             ZSTO1        = ZMOLECO1D*(ZCH3*ZMOLECH2+(ZCN21+ZCN22)*ZMOLECN2O+ &
                  ZCH8*ZMOLECH2O+ZCH7*ZMOLECCH4)
             ZSTO2        = ZMOLECO*(ZCH4*ZMOLECOH+ZCH11*ZMOLECHO2+ &
                  2._dp*(ZCN13*ZMOLECNO2+ZCC1*ZMOLECCLO))- &
                  ZMOLECHO2*ZCH18*ZMOLECH
             ! mz_ab_20100504+
             ! Br parameterisation:
             ! Consider ozone depletion by bromine
!!$          ZSTO3        = ZMOLECO3*(ZCH9*ZMOLECH+ZCH12*ZMOLECHO2+ &
!!$               ZCH14*ZMOLECOH)
             ! mz_bk_20110331+
             ! According to Volker, and Steil et al., 1998, Appendix B, Eq. 1, 
             ! here schould be a + instead of a -. Because we want to solve 
             ! dq/dt = P - Lq, where q is number density (cm^-3) P is
             ! production (cm^-3) and L is loss freqency (s^-1), with an Euler
             ! backward => q_new = (q_old + P_old * dt) / (1 + L_old * dt), 
             ! loss L should be positive (zproxy is positive (ABS))!
!!$          ZSTO3        = ZMOLECO3*(ZCH9*ZMOLECH+ZCH12*ZMOLECHO2+ &
!!$                ZCH14*ZMOLECOH)-zproxy*2._dp*ZDCL2O2A*ZMOLECCL2O2
              ZSTO3        = ZMOLECO3*(ZCH9*ZMOLECH+ZCH12*ZMOLECHO2+ &
                   ZCH14*ZMOLECOH)+zproxy*2._dp*ZDCL2O2A*ZMOLECCL2O2
             ! mz_bk_20110331-
             ! mz_ab_20100504-

             ZSTO4        = 2._dp*ZMOLECO3*(ZCO2*ZMOLECO1D+ZCO1*ZMOLECO)
             ZSTO5        = ZDN2O*ZMOLECN2O+2._dp*ZDO2*ZMOLECO2+ &
                  ZMOLECHO2*ZCN10*ZMOLECNO+ &
                  ZCN19*ZMOLECHNO3*ZMOLECOH+ &
                  ZCN11*ZMOLECCH3O2*ZMOLECNO+ &
                  ZDCO2*ZMOLECCO2
             !     @               -ZDCLOH*ZMOLECCLOH-2._dp*ZMOLECCL2O2*ZDCL2O2A-
             !     @               ZCC3*ZMOLECCLO*ZMOLECOH
             !
             !     EULER BACKWARD LINEAR
             !
             ZMOLECODD    = (ZMOLECODD+ZTMSTDT*ZSTO5)/ &
                  (1._dp+ZTMSTDT/ZMOLECODD* &
                  (ZSTO0+ZSTO1+ZSTO2+ZSTO3+ZSTO4))
             ZDESTODD     = ZTMSTDT*(ZSTO0+ZSTO1+ZSTO2+ZSTO3+ZSTO4)
             !
             !
             !     PARTITIONING OF [OX]
             !
             ZMOLECOX     = ZMOLECODD &
                  -ZMOLECCLO-ZMOLECCLOH &
                  -2._dp*ZMOLECCLNO3-2._dp*ZMOLECCL2O2 &
                  -ZMOLECNO2-ZMOLECHNO4 &
                  -2._dp*ZMOLECNO3-3._dp*ZMOLECN2O5 &
                  -ZMOLECHNO3
             !
             !     CALCULATION OF [H2O]
             !
             ZQH2O        = ZDH2O+ZMOLECO1D*ZCH8
             ZPROH2O      = (ZFUNCPRODH2OLONG (ZCN19,ZMOLECHNO3,ZCC2,ZMOLECHCL, &
                  ZCC14,ZMOLECCH3CL,ZCC15, &
                  ZMOLECCH3CCL3, &
                  ZCH6,ZMOLECH2,ZCH15,ZMOLECH2O2, &
                  ZCH19,ZMOLECCH4)+ &
                  ZFUNCPRODH2OSHORT (ZCN6,ZMOLECHNO4,ZCC11, &
                  ZMOLECCLOH, &
                  ZCH16,ZMOLECCH2O,ZCH20,ZCH21, &
                  ZMOLECCH3O2H))*ZMOLECOH+ &
                  ZFUNCPRODH2OHOX (ZCH18,ZMOLECH,ZMOLECHO2, &
                  ZCH1,ZMOLECOH)
             ZMOLECH2O    = (ZMOLECH2O+ZPROH2O*ZTMSTDT)/(1._dp+ZTMSTDT*ZQH2O)
             Conc(JL,JK,ind_H2O)=ZMOLECH2O 
             !
          ELSE
             ! ******************************************************************
             ! *                                                                *
             ! * nighttime chemistry                                            *
             ! *                                                                *
             ! ******************************************************************
             ! mz_bk_20101213+
             ! ***********************
             ! * Br parameterisation *
             ! ***********************
             zdeltao3_br = 0.0_dp
             ! mz_bk_20101213-

             !
             !     INITIALIZE NIGHTTIME-CHEMISTRY
             !
             ZMOLECO3   = ZMOLECOX
             ZMOLECNO2  = ZMOLECNO+ZMOLECN+ZMOLECNO2
             ZMOLECNO2T0= ZMOLECNO2
             ZMOLECO    = 0._dp
             ZMOLECO1D  = 0._dp
             ZMOLECNO   = 0._dp
             ZMOLECN    = 0._dp
             ZMOLECCLO  = ZMOLECCLO+ZMOLECCL
             ZMOLECCLOT0= ZMOLECCLO
             ZMOLECCL   = 0._dp
             ZMOLECCH3O = 0._dp
             Conc(JL,JK,ind_H2O)=ZMOLECH2O
             !
             !     INTEGRATION OF HOX
             !

             ZMOLECHOXT0  = ZMOLECOH+ZMOLECHO2+ZMOLECH
             ZMOLECH    = 0._dp
             !
             IF(ZMOLECHO2.GT.1._dp) THEN
                !
                ZQOHINFA   = ZCH16*ZMOLECCH2O+ZCH5*ZMOLECCO+ZCH6*ZMOLECH2+ &
                     ZCH14*ZMOLECO3+ZCH15*ZMOLECH2O2
                ZQOHOUFA   = ZCN19*ZFMOLECHNO3+ &
                     ZCN17*ZMOLECNO2+ZCH19*ZMOLECCH4+ZCN6*ZMOLECHNO4+ &
                     ZCC11*ZMOLECCLOH+ZCC14*ZMOLECCH3CL+ &
                     ZCH20*ZMOLECCH3O2H+ZCH1*ZMOLECHO2+ZCC2*ZFMOLECHCL+ &
                     ZCC15*ZMOLECCH3CCL3
                ZRATOHHO2  = (ZCH12*ZMOLECO3)/(ZQOHINFA+ZQOHOUFA)
                ZRATHO2HOX = 1._dp/(1._dp+ZRATOHHO2)
                !
                ZQHO2OUFA  = ZCH22*ZMOLECCH3O2+ZCN16*ZMOLECNO2+ZCC10*ZMOLECCLO
                !
                ZQHOX      = ZRATHO2HOX* &
                     (ZQOHOUFA*ZRATOHHO2+ZQHO2OUFA+2._dp*ZCH2*ZMOLECHO2)
                IF(ZQHOX.LT.1.E-10_dp) THEN
                   ZEQHOX     = 1._dp
                   ZLIMHOX    = ZTMSTDT
                ELSE
                   ZEQHOX     = EXP(-1._dp*ZTMSTDT*ZQHOX)
                   ZLIMHOX    = (1._dp-ZEQHOX)/ZQHOX
                ENDIF
                ZMOLECHOX  = ZMOLECHOXT0*ZEQHOX+ZLIMHOX*ZCN5*ZMOLECHNO4
                !
                ZMOLECHO2  = ZRATHO2HOX*ZMOLECHOX
                ZMOLECOH   = MAX(ZRATOHHO2*ZMOLECHOX,0.00000001_dp)
             ELSE
                ZMOLECHO2 = 0.00000001_dp
                ZMOLECOH  = 0.00000001_dp
             ENDIF
             !
             !     INTEGRATION VON [H2O2]
             !
             ZMOLECH2O2 = (ZMOLECH2O2+ZCH2*ZMOLECHO2*ZMOLECHO2T0*ZTMSTDT)/ &
                  (1._dp+ZCH15*ZMOLECOH*ZTMSTDT)
             !
             !     CALCULATION OF [CH4]
             !
             ZQCH4        = ZCH19*ZMOLECOH
             ZMOLECCH4    = ZMOLECCH4/(1._dp+ZTMSTDT*ZQCH4)
             !
             !     CALCULATION OF [CH3O2]
             !
             ZMOLECCH3O2 = (ZTMSTDT*ZMOLECOH*(ZCH19*ZMOLECCH4+ &
                  ZCH20*ZMOLECCH3O2H)+ZMOLECCH3O2)/ &
                  (1._dp+ZTMSTDT*ZCH22*ZMOLECHO2)
             !
             !     CALCULATION OF [CH3O2H]
             !
             ZMOLECCH3O2H =(ZTMSTDT*ZCH22*ZMOLECCH3O2*ZMOLECHO2+ZMOLECCH3O2H)/ &
                  (1._dp+(ZCH20+ZCH21)*ZMOLECOH*ZTMSTDT)
             !
             !     CALCULATION OF [CH2O]
             !
             ZMOLECCH2O   = (ZMOLECCH2O+ZCH21*ZMOLECCH3O2H*ZMOLECOH*ZTMSTDT)/ &
                  (1._dp+ZTMSTDT*ZCH16*ZMOLECOH)
             !
             !     INTEGRATION OF [CO]
             !

             ZMOLECCO     = (ZCH16*ZMOLECOH*ZMOLECCH2O*ZTMSTDT+ZMOLECCO)/ &
                  (1._dp+ZCH5*ZMOLECOH*ZTMSTDT)
             !
             !     PRODUCTION OF CLOX BY HALOCARBONCHEMISTRY
             !     MASS IS STRICTLY CONSERVED!
             !
             !
             ZPRCH3CL     = ZCC14*ZMOLECOH
             ZPRCH3CL     = ZPRCH3CL*ZTMSTDT
             IF(ZPRCH3CL.LT.1.E-10_dp) THEN
                ZPRCH3CL     = 1._dp
             ELSE
                ZPRCH3CL     = EXP(-1._dp*ZPRCH3CL)
             ENDIF
             !
             ZPRCH3CCL3   = ZCC15*ZMOLECOH
             ZPRCH3CCL3   = ZPRCH3CCL3*ZTMSTDT
             IF(ZPRCH3CCL3.LT.1.E-10_dp) THEN
                ZPRCH3CCL3     = 1._dp
             ELSE
                ZPRCH3CCL3     = EXP(-1._dp*ZPRCH3CCL3)
             ENDIF
             !
             ZPRODCLOY    = 3._dp*ZMOLECCH3CCL3*(1._dp-ZPRCH3CCL3)+ &
                  ZMOLECCH3CL*(1._dp-ZPRCH3CL)
             !
             ZMOLECCH3CL  = ZMOLECCH3CL*ZPRCH3CL
             ZMOLECCH3CCL3= ZMOLECCH3CCL3*ZPRCH3CCL3
             !
             !      ZPRODCLOY    = 0._dp
             ZMOLECCLOY   = ZMOLECCLOX+ZMOLECHCL
             ZMOLECCLOY   = ZPRODCLOY+ZMOLECCLOY
             ZMOLECCLOX   = ZPRODCLOY+ZMOLECCLOX
             !
             !     SOLVE THE COMBINED SYSTEM OF NOX AND CLOX
             !
             !     A) CLOX
             !
             ZQCLNO3      = ZHET2
             IF(ZHET2.LT.1.E-10_dp) THEN
                ZAUX1CLOX   = ZTMSTDT
                ZEQCLNO3    = 1._dp
             ELSE
                ZEQCLNO3    = EXP(-1._dp*ZTMSTDT*ZQCLNO3)
                ZEMEQCLNO3  = 1._dp-ZEQCLNO3
                ZAUX1CLOX   = ZEMEQCLNO3/ZQCLNO3
             ENDIF
             ZQCLOH       = ZMOLECOH*ZCC11+ZCP3
             IF(ZQCLOH.LT.1.E-10_dp) THEN
                ZEQCLOH      = 1._dp
                ZAUX2CLOX    = ZTMSTDT
             ELSE
                ZEQCLOH      = EXP(-1._dp*ZTMSTDT*ZQCLOH)
                ZAUX2CLOX    = (1._dp-ZEQCLOH)/ZQCLOH
             ENDIF
             ZAUX3CLOX    = ZCP3/ZQCLOH
             ZAUX4CLOX    = 1._dp-ZAUX3CLOX
             ZAUX10CLOX   = 1._dp+ZTMSTDT*ZCCL2O2M
             !
             !     WEICHE HETEROGENE CHEMIE ODER NICHT!
             !
             IF(ZCP1.LT.1E-20_dp) THEN
                ZAUX5CLOX   = 0._dp
             ELSE
                ZAUX5CLOX   = ZCP1/ZQCLNO3
             ENDIF
             !
             ZAUX6CLOX   = 1._dp-ZAUX5CLOX+ZAUX4CLOX*ZHET3*ZAUX2CLOX+ &
                  ZAUX3CLOX*ZTMSTDT*ZHET3
             ZAUX7CLOX   = ZAUX5CLOX+ZEQCLNO3*ZAUX6CLOX
             ZWCLOX      = ZCN14*(ZAUX5CLOX*ZTMSTDT+ZAUX1CLOX*ZAUX6CLOX) 
             ZUCLOX       = ZMOLECCLOXT0-2._dp*ZMOLECCL2O2/ZAUX10CLOX- &
                  2._dp*ZMOLECCL2- &
                  ZMOLECCLOH*(ZAUX3CLOX+ZAUX4CLOX*ZEQCLOH)- &
                  ZMOLECCLNO3*ZAUX7CLOX
             IF(ZUCLOX.LE.0.) ZUCLOX = ZMOLECCLO
             ZVCLOX       = 1._dp+2._dp*ZCCLOCLO*ZMOLECCLO*ZTMSTDT/ZAUX10CLOX+ &
                  (ZAUX4CLOX*ZAUX2CLOX+ZAUX3CLOX*ZTMSTDT)* &
                  ZCC10*ZMOLECHO2
             !      ZWV(JL) = ZWCLOX
             !      ZVV(JL) = ZVCLOX
             !      ZUV(JL) = ZUCLOX
             !
             !     B) NOX
             !
             ZQHNO4       = ZCN5+ZCN6*ZMOLECOH
             Z1QHNO4      = 1._dp/(1._dp+ZTMSTDT*ZQHNO4)
             ZQNO3        = ZCN18*ZMOLECNO2
             Z1QNO3       = 1._dp/(1._dp+ZTMSTDT*ZQNO3)
             ZQN2O5       = ZHET1+ZCN7*Z1QNO3
             Z1QN2O5      = 1._dp/(1._dp+ZTMSTDT*ZQN2O5)
             ZAUX1NOX     = (4._dp-2._dp*EXP(-1._dp*ZTMSTDT*ZHET1)+ &
                  ZTMSTDT*ZCN7*Z1QNO3)*Z1QN2O5
             ZKNOX        = ZMOLECNOX-ZMOLECCLNO3-ZMOLECHNO4*Z1QHNO4- &
                  ZMOLECNO3*Z1QNO3- &
                  ZMOLECN2O5*ZAUX1NOX
             !     @               ZMOLECN2O5*(ZAUX1NOX*Z1QN2O5+ &
             !     @               2._dp*(1.-EXP(-1._dp*ZTMSTDT*ZHET1)))
             !
             !      ZKNOXV(JL) = ZKNOX
             !
             !     PARACHUTE
             !
             IF(ZKNOX.GT.0._dp) THEN
                ZANOX        = 1._dp+ZCN15*ZMOLECO3*ZTMSTDT*Z1QNO3+ &
                     ZAUX1NOX*ZCN18*ZMOLECNO3*ZTMSTDT*Z1QNO3+ &
                     ZCN16*ZMOLECHO2*Z1QHNO4*ZTMSTDT
                ZBNOX        = ZCN18*ZCN15*ZMOLECO3*ZTMSTDT*ZTMSTDT*ZAUX1NOX* &
                     Z1QNO3
                ZCNOX        = ZCN14*ZTMSTDT
                ZC0          = ZKNOX*ZVCLOX
                ZC1          = ZANOX*ZVCLOX+ZCNOX*ZUCLOX-ZKNOX*ZWCLOX
                ZC2          = ZANOX*ZWCLOX+ZBNOX*ZVCLOX
                ZC3          = ZBNOX*ZWCLOX
                !       ZANOXV(JL) = ZANOX
                !       ZBNOXV(JL) = ZBNOX
                !       ZCNOXV(JL) = ZCNOX
                !       ZDNOXV(JL) = 3333333333.
                !       ZC0V(JL) = ZC0
                !       ZC1V(JL) = ZC1
                !       ZC2V(JL) = ZC2
                !       ZC3V(JL) = ZC3
                !
                !     PARACHUTE
                !
                ZMAX        = ZC2*ZC2
                ZMAX1       = ZMAX-3._dp*ZC1*ZC3
                IF(ZMAX1.GT.ZMAX) THEN
                   ZMOLECNO2  = MAX((-1._dp*ZC2+SQRT(ZMAX1))/ZC3,ZMOLECNO2)
                ENDIF
                !
                !     SOLVE WITH NEWTON-ITERATION!
                !
                !     1.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                ZMOLECNO2   = MAX(1._dp,ZMOLECNO2)
                !
                !     2.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !     3.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !     4.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !     5.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                ZMOLECNO2   = MAX(1._dp,ZMOLECNO2)
                !
                !     6.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !     7.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !     8.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !     9.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                !
                !     10.ITERATION
                !
                ZAUXNNEW1   = ZC3*ZMOLECNO2**2._dp
                ZAUXNNEW2   = ZC2*ZMOLECNO2
                ZMOLECNO2   = (ZMOLECNO2*(2._dp*ZAUXNNEW1+ZAUXNNEW2)+ZC0)/ &
                     (3._dp*ZAUXNNEW1+2._dp*ZAUXNNEW2+ZC1)
                ZMOLECNO2   = MAX(1._dp,ZMOLECNO2)
             ENDIF
             !
             ZMOLECCLO   = ZUCLOX/(ZVCLOX+ZWCLOX*ZMOLECNO2)
             ZPCLNO3     = ZCN14*ZMOLECNO2*ZMOLECCLO
             ZMOLECCLNO3 = ZMOLECCLNO3T0*ZEQCLNO3+ZPCLNO3*ZAUX1CLOX
             ZPCLOH      = ZCC10*ZMOLECHO2*ZMOLECCLO+ZMOLECCLNO3*ZHET3
             ZMOLECCLOH  = ZMOLECCLOHT0*ZEQCLOH+ZAUX2CLOX*ZPCLOH
             ZMOLECCL2O2 = (ZCCLOCLO*ZMOLECCLOT0*ZMOLECCLO*ZTMSTDT+ &
                  ZMOLECCL2O2)/ZAUX10CLOX
             !
             !     HETEROGENER VERLUST VON HCL
             !
             IF((ZCP1+ZCP3).GT.1.E-20_dp) THEN
                ZHETLOSSHCL  = ZAUX5CLOX*(ZMOLECCLNO3T0-ZMOLECCLNO3+ &
                     ZPCLNO3*ZTMSTDT)+ &
                     ZAUX3CLOX*(ZMOLECCLOHT0-ZMOLECCLOH+ &
                     ZPCLOH*ZTMSTDT)
                ZHETLOSSHCL  = MIN(ZHETLOSSHCL,ZMOLECHCLT0)
                ZMOLECHCL    = ZMOLECHCLT0-ZHETLOSSHCL
                ZMOLECCLOX   = ZMOLECCLOXT0+ZHETLOSSHCL
                ZMOLECCL2    = ZMOLECCL2+ZHETLOSSHCL
                !       IF(ZHETLOSSHCL.LT.0.) THEN
                !       IDIAGV(JL) = 1
                !       ZHETLOSSHCLV(JL,1) = ZAUX5CLOX
                !       ZHETLOSSHCLV(JL,2) = -123456.
                !       ZHETLOSSHCLV(JL,3) = ZMOLECCLNO3T0
                !       ZHETLOSSHCLV(JL,4) = ZMOLECCLNO3
                !       ZHETLOSSHCLV(JL,5) = ZPCLNO3
                !       ZHETLOSSHCLV(JL,6) = ZAUX3CLOX
                !       ZHETLOSSHCLV(JL,7) = -123456.
                !       ZHETLOSSHCLV(JL,8) = ZMOLECCLOHT0
                !       ZHETLOSSHCLV(JL,9) = ZMOLECCLOH
                !       ZHETLOSSHCLV(JL,10) = ZPCLOH
                !       ENDIF
             ENDIF
             !
             ZMOLECCLOX = MAX(1._dp,ZMOLECCLOX)
             ZMOLECHCL  = MAX(1._dp,ZMOLECHCL)
             ZSCY       = ZMOLECCLOY/(ZMOLECCLOX+ZMOLECHCL)
             ZMOLECCLOX = ZSCY*ZMOLECCLOX
             ZMOLECHCL  = ZSCY*ZMOLECHCL
             !
             ZQNO3        = ZCN18*ZMOLECNO2
             Z1QNO3       = 1._dp/(1._dp+ZTMSTDT*ZQNO3)
             ZQN2O5       = ZHET1+ZCN7*Z1QNO3
             Z1QN2O5      = 1._dp/(1._dp+ZTMSTDT*ZQN2O5)
             ZMOLECHNO4  = ZMOLECHNO4T0*Z1QHNO4+ &
                  ZMOLECHO2*ZMOLECNO2*Z1QHNO4*ZTMSTDT &
                  *ZCN16
             ZPN2O5       = ZMOLECNO2*ZCN18*Z1QNO3* &
                  (ZMOLECNO3+ &
                  ZMOLECNO2*ZCN15*ZMOLECO3*ZTMSTDT)
             ZMOLECN2O5   = ZMOLECN2O5*Z1QN2O5+ZTMSTDT*ZPN2O5*Z1QN2O5
             ZPNO3        = ZCN7*ZMOLECN2O5+ZCN15*ZMOLECO3*ZMOLECNO2
             ZMOLECNO3    = ZMOLECNO3*Z1QNO3+ZTMSTDT*ZPNO3*Z1QNO3
             !
             !     HETEROGENER VERLUST VON NOX
             !
             IF(ZHET1.GT.1.E-20_dp) THEN
                IF(ZQNO3.LT.1.E-10) THEN
                   ZEQNO3       = 1._dp
                   ZAUX10NOX    = ZTMSTDT
                ELSE
                   ZEQNO3       = EXP(-1._dp*ZTMSTDT*ZQNO3)
                   ZAUX10NOX    = (1._dp-ZEQNO3)/ZQNO3
                ENDIF
                ZQN2O5       = ZHET1+ZCN7*ZEQNO3
                ZPN2O5       = ZMOLECNO2*ZCN18*(ZMOLECNO3*ZEQNO3+ &
                     ZMOLECNO2*ZCN15*ZMOLECO3*ZAUX10NOX)
                IF(ZQN2O5.LT.1.E-10_dp) THEN
                   ZEQN2O5   = 1._dp
                   ZEMEQN2O5 = ZTMSTDT
                ELSE
                   ZEQN2O5   = EXP(-1._dp*ZTMSTDT*ZQN2O5)
                   ZEMEQN2O5 = (1._dp-ZEQN2O5)/ZQN2O5
                ENDIF
                ZHETLOSS1    = (ZMOLECN2O5T0*(1._dp-ZEQN2O5)+ &
                     ZPN2O5*(ZTMSTDT-ZEMEQN2O5))/ZQN2O5
                ZHETLOSS2    = ZMOLECCLNO3T0-ZMOLECCLNO3+ZPCLNO3*ZTMSTDT
                ZHETLOSSNOXTEST  = 2._dp*ZHET1*ZHETLOSS1+ZHETLOSS2
                ZMOLECNOX    = MAX(2._dp,ZMOLECNOX)
                ZHETLOSSNOX  = MIN(ZMOLECNOX-1.,ZHETLOSSNOXTEST)
                SCALE        = ZHETLOSSNOX/ZHETLOSSNOXTEST
                ZHETLOSS1    = SCALE*ZHETLOSS1
                ZHETLOSS2    = SCALE*ZHETLOSS2
                !
                !       ZMOLECNOX    = 2._dp*ZMOLECN2O5+ZMOLECNO3+ZMOLECHNO4+ &
                !     @                ZMOLECNO2+ZMOLECNO+ZMOLECN+ZMOLECCLNO3
                ZMOLECNOX    = ZMOLECNOX-ZHETLOSSNOX
                ZMOLECHNO3   = 2._dp*ZCS2*ZHETLOSS1+ZCS1/ZQCLNO3*ZHETLOSS2+ &
                     ZMOLECHNO3
                ZMOLECNAT    = 2._dp*ZCP4*ZHETLOSS1+(ZCP1+ZCP2)/ZQCLNO3*ZHETLOSS2+ &
                     ZMOLECNAT
             ENDIF
             !
             !     MASS-CONSERVATION OF CL-ATOMS
             !     (FIRST CL-ATOMS, BECAUSE OF CLNO3)
             !
             ZSTOSCALCLOX    = ZMOLECCLOX/(ZMOLECCLO+ZMOLECCLOH+ &
                               ZMOLECCLNO3+2._dp*ZMOLECCL2O2+2._dp*ZMOLECCL2)
             ZMOLECCL2       = ZMOLECCL2*ZSTOSCALCLOX
             ZMOLECCLO       = ZMOLECCLO*ZSTOSCALCLOX
             ZMOLECCLOH      = ZMOLECCLOH*ZSTOSCALCLOX
             ZMOLECCLNO3     = ZMOLECCLNO3*ZSTOSCALCLOX
             ZMOLECCL2O2     = ZMOLECCL2O2*ZSTOSCALCLOX
             !
             !     MASSCONSERVATION OF N-ATOMS
             !
             ZMOLECNOY     = ZMOLECNOX+ZMOLECHNO3
             !      ZDELATOMSNV(JL) = 0._dp
             !
             !     HOMOGENEOUS INTEGRATION OF NOX AND HNO3
             !
             ZP1HNO3       = ZMOLECNO2/ZMOLECNOX*ZMOLECOH*ZCN17
             ZP1NOX        = ZCN19*ZMOLECOH
             ZQHNO3        = ZP1NOX+ZP1HNO3
             IF(ZQHNO3.LT.1.E-10_dp) THEN
                ZEQHNO3       = 1._dp
                ZLIMHNO3      = ZTMSTDT
             ELSE
                ZEQHNO3       = EXP(-1._dp*ZTMSTDT*ZQHNO3)
                ZLIMHNO3      = (1._dp-ZEQHNO3)/ZQHNO3
             ENDIF
             ZQNOX         = ZQHNO3
             IF(ZQNOX.LT.1.E-10_dp) THEN
                ZEQNOX        = 1._dp
                ZLIMNOX       = ZTMSTDT
             ELSE
                ZEQNOX        = EXP(-1._dp*ZTMSTDT*ZQNOX)
                ZLIMNOX       = (1._dp-ZEQNOX)/ZQNOX
             ENDIF
             ZMOLECHNO3    = ZMOLECHNO3*ZEQHNO3+ &
                  ZLIMHNO3*ZMOLECNOY*ZP1HNO3
             ZMOLECNOX     = ZMOLECNOX*ZEQNOX+ &
                  ZLIMNOX*(ZMOLECNOY*ZP1NOX)
             !
             !     MASS-CONSERVATION OF N-ATOMS
             !
             ZSCALATOMSN  = ZMOLECNOY/(ZMOLECNOX+ZMOLECHNO3)
             ZMOLECHNO3   = ZMOLECHNO3*ZSCALATOMSN
             ZMOLECNOX    = ZMOLECNOX*ZSCALATOMSN
             !
             !     MASS-CONSERVATION FOR NOX AND CLOX FAMILY-MEMBERS
             !     WHEN CLNO3 TOO BIG !
             !
             !      ZSWMBTOTNOXV2(JL) = 0.
             IF ((ZMOLECNOX-ZMOLECCLNO3).LE.0._dp) THEN
                ZSTOSCALTOT =(ZMOLECNOX+ZMOLECCLOX)/ &
                     (ZMOLECNO2+ZMOLECNO3+ &
                     ZMOLECHNO4+2._dp*ZMOLECN2O5+2._dp*ZMOLECCL2+ZMOLECCLO+ &
                     ZMOLECCLOH+2._dp*ZMOLECCLNO3+2._dp*ZMOLECCL2O2)
                ZMOLECCL2   = ZMOLECCL2*ZSTOSCALTOT
                ZMOLECCLO   = ZMOLECCLO*ZSTOSCALTOT
                ZMOLECCLOH  = ZMOLECCLOH*ZSTOSCALTOT
                ZMOLECCLNO3 = ZMOLECCLNO3*ZSTOSCALTOT
                ZMOLECCL2O2 = ZMOLECCL2O2*ZSTOSCALTOT
                ZMOLECNO2   = ZMOLECNO2*ZSTOSCALTOT
                ZMOLECNO3   = ZMOLECNO3*ZSTOSCALTOT
                ZMOLECHNO4  = ZMOLECHNO4*ZSTOSCALTOT
                ZMOLECN2O5  = ZMOLECN2O5*ZSTOSCALTOT
                !       ZSWMBTOTNOXV2(JL) = 1.
                !
                !      MASSENBILANZ CLOY NOCH STRENGER !
                !
                ZSTOSCALCLOX    = ZMOLECCLOX/(ZMOLECCLO+ZMOLECCLOH+ &
                                  ZMOLECCLNO3+2._dp*ZMOLECCL2O2+2._dp*ZMOLECCL2)
                ZMOLECCL2       = ZMOLECCL2*ZSTOSCALCLOX
                ZMOLECCLO       = ZMOLECCLO*ZSTOSCALCLOX
                ZMOLECCLOH      = ZMOLECCLOH*ZSTOSCALCLOX
                ZMOLECCLNO3     = ZMOLECCLNO3*ZSTOSCALCLOX
                ZMOLECCL2O2     = ZMOLECCL2O2*ZSTOSCALCLOX
             ELSE
                !
                !     MASS-CONSERVATION OF NOX-FAMILY-MEMBERS
                !
                ZSTOSCALNOX    = (ZMOLECNOX-ZMOLECCLNO3)/ &
                     (ZMOLECN+ZMOLECNO+ZMOLECNO2+ZMOLECNO3+ &
                     ZMOLECHNO4+2._dp*ZMOLECN2O5)
                ZMOLECN         = ZMOLECN*ZSTOSCALNOX
                ZMOLECNO        = ZMOLECNO*ZSTOSCALNOX
                ZMOLECNO2       = ZMOLECNO2*ZSTOSCALNOX
                ZMOLECNO3       = ZMOLECNO3*ZSTOSCALNOX
                ZMOLECHNO4      = ZMOLECHNO4*ZSTOSCALNOX
                ZMOLECN2O5      = ZMOLECN2O5*ZSTOSCALNOX
             ENDIF
             !
             !     INTEGRATION OF [ODD], CALCULATION OF [O3]
             !
             ZSTO0        = ZMOLECCLNO3*(2._dp*ZCP1+ZCP2)+ &
                  ZMOLECCLOH*ZCP3+ &
                  ZMOLECN2O5*(3._dp*ZCP4+ZCS2)
             ZSTO1        = 1._dp+ZMOLECO3*ZTMSTDT/ZMOLECODD* &
                  (ZCH12*ZMOLECHO2+ZCH14*ZMOLECOH) &
                  +ZTMSTDT/ZMOLECODD*ZSTO0
             ZSTO2        = ZTMSTDT*ZMOLECOH* &
                  ZCN19*ZMOLECHNO3+ZMOLECODD
             ZDESTODD     = ZTMSTDT*(ZSTO0+ZMOLECO3* &
                  (ZCH12*ZMOLECHO2+ZCH14*ZMOLECOH))
             ZMOLECODD    = ZSTO2/ZSTO1
             ZMOLECOX     = ZMOLECODD-ZMOLECCLO-ZMOLECCLOH- &     
                  2._dp*(ZMOLECCLNO3+ZMOLECCL2O2)-ZMOLECNO2- &
                  ZMOLECHNO3-ZMOLECHNO4-3._dp*ZMOLECN2O5-2._dp*ZMOLECNO3
          ENDIF ifdani
          !
          !     PARACHUTE !
          !     BECAUSE OF NOX-EMISSIONS AND DRY DEPOSITION OF O3
          !     TAKE CARE!
          !
          !      ZSCDOXV(JL) = 1._dp
          IF(ZMOLECOX.LE.0._dp) THEN
             ZMOLECOX = 0.2_dp*ZMOLECOXT0 
             !      ZSCDOXV(JL) = -1._dp
          ENDIF
          !
          !     SEDIMENTATION !
          !
          !      ZTESTICE  = ZMOLECICE+(ZUDTBICE(JL)-ZUDTAICE(JL))/ZHEILEV(JL)
          !      IF(ZTESTICE.LT.1.) THEN
          !        ZUDTAICE(JL) = ZMOLECICE*ZHEILEV(JL)+ZUDTBICE(JL)
          !        ZTESTICE     = 1.
          !      ENDIF
          !      ZMOLECICE = ZTESTICE
          ZMOLECICE = ZMOLECICE+(ZUDTBICE(JL)-ZUDTAICE(JL))/ZHEILEV(JL)
          ZMOLECICE = MAX(ZMOLECICE,1._dp)
          !
          !      ZTESTNAT  = ZMOLECNAT+(ZUDTBNAT(JL)-ZUDTANAT(JL))/ZHEILEV(JL)
          !      IF(ZTESTNAT.LT.1.) THEN
          !        ZUDTANAT(JL) = ZMOLECNAT*ZHEILEV(JL)+ZUDTBNAT(JL)
          !        ZTESTNAT     = 1.
          !      ENDIF
          !      ZMOLECNAT = ZTESTNAT
          ZMOLECNAT = ZMOLECNAT+(ZUDTBNAT(JL)-ZUDTANAT(JL))/ZHEILEV(JL)
          ZMOLECNAT = MAX(ZMOLECNAT,1._dp)
          !
          ! mz_ab_20100504+
          ! Br parameterisation:
          ! Save ozone depletion by bromine !!
          zdeltao3_brv(JL,JK) = zdeltao3_br
          ! mz_ab_20100504-

          ! mz_bk_20101222+
          ! original code by TS/VG(DLR)
          ! changed to units molec/sec
          ! Production: ZPRODO2  Photolysis of molecular oxigen
          !             ZPRODCO  Ozoneproduction by Carbonmonoxide oxidation
          !             ZPRODCH4 Ozoneproduction by Methane oxidation
          ZPRODO2(jl,jk)     = 2.*ZDO2*ZMOLECO2
          ZPRODCO(jl,jk)     = ZCN10*ZMOLECHO2*ZMOLECNO
          ZPRODCH4(jl,jk)    = ZCN11*ZMOLECCH3O2*ZMOLECNO
          ! Ozone destruction:
          ZDESTH12(jl,jk)    = ZCH12*ZMOLECO3*ZMOLECHO2
          ZDESTH14(jl,jk)    = ZCH14*ZMOLECO3*ZMOLECOH
          ZDESTN13(jl,jk)    = 2.*ZCN13*ZMOLECO*ZMOLECNO2
          ZDESTC1(jl,jk)     = 2.*ZCC1*ZMOLECCLO*ZMOLECO
          ZDESTCL2O2(jl,jk) = 2.*ZDCL2O2A*ZMOLECCL2O2
          ZDESTCLOH(jl,jk)  = ZDCLOH*ZMOLECCLOH
          ZDESTH8(jl,jk)     = ZCH8*ZMOLECH2O*ZMOLECO1D
          ! little impact; important in upper stratosphere and mesosphere
          ZDESTH4(jl,jk)     = ZCH4*ZMOLECOH*ZMOLECO
          ZDESTH11(jl,jk)    = ZCH11*ZMOLECHO2*ZMOLECO
          ZDESTO1(jl,jk)     = ZCO1*2.*ZMOLECO*ZMOLECO3
          ! mz_bk_20101222-

          ! mz_ab_20090915+ partitioning
          ZMOLECO3V(JL)     =  ZSCALEO3(JL)*ZMOLECOX
          ! mz_ab_20090915+ 
          ZMOLECOV(JL)      = ZMOLECO
          ZMOLECO1DV(JL)    = ZMOLECO1D
          ZMOLECHV(JL)      = ZMOLECH
          ZMOLECOHV(JL)     = ZMOLECOH
          ZMOLECHO2V(JL)    = ZMOLECHO2
          ZMOLECNV(JL)      = ZMOLECN
          ZMOLECNOV(JL)     = ZMOLECNO
          ZMOLECNO2V(JL)    = ZMOLECNO2
          ZMOLECNO3V(JL)    = ZMOLECNO3
          ZMOLECHNO4V(JL)   = ZMOLECHNO4
          ZMOLECN2O5V(JL)   = ZMOLECN2O5
          ZMOLECCLV(JL)     = ZMOLECCL
          ZMOLECCLOV(JL)    = ZMOLECCLO
          ZMOLECCLOHV(JL)   = ZMOLECCLOH
          ZMOLECCL2O2V(JL)  = ZMOLECCL2O2
          ZMOLECCLNO3V(JL)  = ZMOLECCLNO3
          ZMOLECCL2V(JL)    = ZMOLECCL2
          ZMOLECCH2OV(JL)   = ZMOLECCH2O
          ZMOLECCH3O2V(JL)  = ZMOLECCH3O2
          ZMOLECCH3OV(JL)   = ZMOLECCH3O
          ZMOLECCFCL3V(JL)  = ZMOLECCFCL3
          ZMOLECCF2CL2V(JL) = ZMOLECCF2CL2
          ZMOLECCCL4V(JL)   = ZMOLECCCL4
          ZMOLECCH3CLV(JL)  = ZMOLECCH3CL
          ZMOLECCH3CCL3V(JL)= ZMOLECCH3CCL3
          ZMOLECH2V(JL)     = ZMOLECH2
          ZMOLECCH4V(JL)    = ZMOLECCH4
          ZMOLECN2OV(JL)    = ZMOLECN2O
          ZMOLECH2O2V(JL)   = ZMOLECH2O2
          ZMOLECHCLV(JL)    = ZMOLECHCL
          ZMOLECCOV(JL)     = ZMOLECCO
          ZMOLECCH3O2HV(JL) = ZMOLECCH3O2H
          ZMOLECNATV(JL)    = ZMOLECNAT
          ZMOLECICEV(JL)    = ZMOLECICE
          ZMOLECHNO3V(JL)   = ZMOLECHNO3

          ZMOLECOXV(JL)   =  ZMOLECOX
          ZMOLECCLOXV(JL) =  ZMOLECCLOX
          ZMOLECNOXV(JL)  =  ZMOLECNOX
          ZMOLECODDV(JL)  =  ZMOLECODD
          ZMOLECH2OV(JL)  = ZMOLECH2O

          ! ******************************************************************
          ! *                                                                *
          ! * end loop over longitudes                                       *
          ! *                                                                *
          ! ******************************************************************
       ENDDO

       DO JL = 1,KLON

          ! ******************************************************************
          ! *                                                                *
          ! *   Copy values to output argument                               *
          ! *                                                                *
          ! ******************************************************************
          Conc(JL,JK,ind_O3)      = ZMOLECO3V(JL) 
          Conc(JL,JK,ind_O3P)     = ZMOLECOV(JL)                  
          Conc(JL,JK,ind_O1D)     = ZMOLECO1DV(JL)                
          Conc(JL,JK,ind_H)       = ZMOLECHV(JL)                  
          Conc(JL,JK,ind_OH)      = ZMOLECOHV(JL)                 
          Conc(JL,JK,ind_HO2)     = ZMOLECHO2V(JL)                
          Conc(JL,JK,ind_N)       = ZMOLECNV(JL)                  
          Conc(JL,JK,ind_NO)      = ZMOLECNOV(JL)                 
          Conc(JL,JK,ind_NO2)     = ZMOLECNO2V(JL)                
          Conc(JL,JK,ind_NO3)     = ZMOLECNO3V(JL)                
          Conc(JL,JK,ind_HNO4)    = ZMOLECHNO4V(JL)               
          Conc(JL,JK,ind_N2O5)    = ZMOLECN2O5V(JL)               
          Conc(JL,JK,ind_CL)      = ZMOLECCLV(JL)                 
          Conc(JL,JK,ind_CLO)     = ZMOLECCLOV(JL)                
          Conc(JL,JK,ind_HOCl)    = ZMOLECCLOHV(JL)               
          Conc(JL,JK,ind_CL2O2)   = ZMOLECCL2O2V(JL)              
          Conc(JL,JK,ind_ClNO3)   = ZMOLECCLNO3V(JL)              
          Conc(JL,JK,ind_CL2)     = ZMOLECCL2V(JL)                
          Conc(JL,JK,ind_HCHO)    = ZMOLECCH2OV(JL)               
          Conc(JL,JK,ind_CH3O2)   = ZMOLECCH3O2V(JL)              
          Conc(JL,JK,ind_CFCL3)   = ZMOLECCFCL3V(JL)              
          Conc(JL,JK,ind_CF2CL2)  = ZMOLECCF2CL2V(JL)             
          Conc(JL,JK,ind_CCL4)    = ZMOLECCCL4V(JL)               
          Conc(JL,JK,ind_CH3CL)   = ZMOLECCH3CLV(JL)              
          Conc(JL,JK,ind_CH3CCL3) = ZMOLECCH3CCL3V(JL)            
          Conc(JL,JK,ind_H2)      = ZMOLECH2V(JL)                 
          Conc(JL,JK,ind_CH4)     = ZMOLECCH4V(JL)                
          Conc(JL,JK,ind_N2O)     = ZMOLECN2OV(JL)                
          Conc(JL,JK,ind_H2O2)    = ZMOLECH2O2V(JL)               
          Conc(JL,JK,ind_HCL)     = ZMOLECHCLV(JL)                
          Conc(JL,JK,ind_CO)      = ZMOLECCOV(JL)                 
          Conc(JL,JK,ind_CH3OOH)  = ZMOLECCH3O2HV(JL)             
          Conc(JL,JK,ind_NAT)     = ZMOLECNATV(JL)                
          Conc(JL,JK,ind_ICE)     = ZMOLECICEV(JL)                
          Conc(JL,JK,ind_HNO3)    = ZMOLECHNO3V(JL)
       ENDDO

       ! ******************************************************************
       ! *                                                                *
       ! * end loop over levels                                           *
       ! *                                                                *
       ! ******************************************************************
    ENDDO

    ! mz_ab_20100504+
 !!$   DEALLOCATE(zdeltao3_brv) ! mz_bk_20101222
    ! mz_ab_20100504-
    DEALLOCATE(ZTEMPV ,ZPRESV ,ZCONV ,ZINVCONV )
    DEALLOCATE(ZMOLECO3V,ZMOLECOV, ZMOLECO1DV       & 
         , ZMOLECOXV, ZMOLECHV, ZMOLECOHV, ZMOLECHO2V &
         , ZMOLECNV, ZMOLECNOV, ZMOLECNO2V &
         , ZMOLECNO3V, ZMOLECHNO4V, ZMOLECN2O5V &
         , ZMOLECCLV, ZMOLECCLOV, ZMOLECCLOHV &
         , ZMOLECCL2O2V, ZMOLECCLNO3V, ZMOLECCL2V &
         , ZMOLECCH2OV, ZMOLECCH3O2V, ZMOLECCH3OV &
         , ZMOLECCFCL3V, ZMOLECCF2CL2V, ZMOLECCCL4V &
         , ZMOLECCH3CLV, ZMOLECCH3CCL3V, ZMOLECH2V &
         , ZMOLECCH4V, ZMOLECN2OV, ZMOLECH2O2V &
         , ZMOLECHCLV, ZMOLECCOV, ZMOLECCH3O2HV &
         , ZMOLECNATV, ZMOLECICEV, ZMOLECHNO3V &
         , ZMOLECCLOXV, ZMOLECNOXV, ZMOLECODDV &
         , ZFMOLECNOXV, ZFMOLECHNO3V, ZFMOLECHCLV &
         , ZFMOLECCLOXV, ZMOLECH2OV  &
         , ZMOLECH2OT0V)
    DEALLOCATE(ZSCALEO3)
    DEALLOCATE(ZMOLECO2V,ZMOLECN2V,ZMOLECCO2V,JINDV)
    DEALLOCATE(ZCP1V,ZCP2V,ZCP3V,ZCP4V, ZCS1V,ZCS2V, ZHET1V,ZHET2V,ZHET3V)
    DEALLOCATE(ZUDTAICE,ZUDTBICE, ZUDTANAT,ZUDTBNAT, ZHEILEV)
    DEALLOCATE(INIHET0V)
    DEALLOCATE(ZNUMICE,ZHELPICE)
    DEALLOCATE(ZRNATUP,ZRICEUP)
    DEALLOCATE(IIHET,INHET,INSED)

    RETURN

      CONTAINS

    FUNCTION ZFUNCRATNONO2 (ZDNO2,ZCN13,ZMOLECO,ZCN8,ZMOLECO3,ZCN9, &
         ZMOLECCLO,ZCN10,ZMOLECHO2,ZCN11,ZMOLECCH3O2)
      REAL(dp) :: ZFUNCRATNONO2
      REAL(dp), INTENT(IN) :: ZDNO2,ZCN13,ZMOLECO,ZCN8,ZMOLECO3,ZCN9,   &
           ZMOLECCLO,ZCN10,ZMOLECHO2,ZCN11,ZMOLECCH3O2
      ZFUNCRATNONO2 = (ZDNO2+ZCN13*ZMOLECO)/            &
           (ZCN8*ZMOLECO3+ZCN9*ZMOLECCLO+      &
           ZCN10*ZMOLECHO2+ZCN11*ZMOLECCH3O2)
    END FUNCTION ZFUNCRATNONO2


    !
    FUNCTION ZFUNCRATNNO (ZDNO,ZCN1,ZMOLECO2,ZSCNNO2,ZMOLECNO2,    &
         ZSCN3,ZMOLECOH,ZSCN4,ZMOLECHO2,ZSCN2,ZMOLECNO)
      REAL(dp) :: ZFUNCRATNNO
      REAL(dp), INTENT(IN) :: ZDNO,ZCN1,ZMOLECO2,ZSCNNO2,ZMOLECNO2,   &
           ZSCN3,ZMOLECOH,ZSCN4,ZMOLECHO2,ZSCN2,ZMOLECNO
      ZFUNCRATNNO  = ZDNO/(ZCN1*ZMOLECO2+ZSCNNO2*ZMOLECNO2+        &
           ZSCN3*ZMOLECOH+ZSCN4*ZMOLECHO2+ZSCN2*ZMOLECNO)
    END FUNCTION ZFUNCRATNNO
    !
    FUNCTION ZFUNCPRODH2OLONG (ZCN19,ZMOLECHNO3,ZCC2,ZMOLECHCL,      &
         ZCC14,ZMOLECCH3CL,ZCC15,ZMOLECCH3CCL3,     &
         ZCH6,ZMOLECH2,ZCH15,ZMOLECH2O2,          &
         ZCH19,ZMOLECCH4)
      REAL(dp) :: ZFUNCPRODH2OLONG
      REAL(dp), INTENT(IN) :: ZCN19,ZMOLECHNO3,ZCC2,ZMOLECHCL,        &
           ZCC14,ZMOLECCH3CL,ZCC15,ZMOLECCH3CCL3,   &
           ZCH6,ZMOLECH2,ZCH15,ZMOLECH2O2,           &
           ZCH19,ZMOLECCH4
      ZFUNCPRODH2OLONG  = ZCN19*ZMOLECHNO3+ZCC2*ZMOLECHCL+   &
           ZCC14*ZMOLECCH3CL+ZCC15*ZMOLECCH3CCL3+       &
           ZCH6*ZMOLECH2+ZCH15*ZMOLECH2O2+            &
           ZCH19*ZMOLECCH4
    END FUNCTION ZFUNCPRODH2OLONG

    !
    FUNCTION ZFUNCPRODH2OSHORT (ZCN6,ZMOLECHNO4,ZCC11,ZMOLECCLOH,     &
         ZSCH16,ZMOLECCH2O,ZSCH20,ZSCH21,        &
         ZMOLECCH3O2H)
      REAL(dp) :: ZFUNCPRODH2OSHORT
      REAL(dp), INTENT(IN) :: ZCN6,ZMOLECHNO4,ZCC11,ZMOLECCLOH,       &
           ZSCH16,ZMOLECCH2O,ZSCH20,ZSCH21,        &
           ZMOLECCH3O2H
      ZFUNCPRODH2OSHORT   = ZCN6*ZMOLECHNO4+ZCC11*ZMOLECCLOH+      &
           ZSCH16*ZMOLECCH2O+(ZSCH20+ZSCH21)*ZMOLECCH3O2H
    END FUNCTION ZFUNCPRODH2OSHORT
    !
    FUNCTION    ZFUNCPRODH2OHOX (ZSCH18,ZMOLECH,ZMOLECHO2,ZCH1,ZMOLECOH)
      REAL(dp) :: ZFUNCPRODH2OHOX
      REAL(dp), INTENT(IN) ::  ZSCH18,ZMOLECH,ZMOLECHO2,ZCH1,ZMOLECOH
      ZFUNCPRODH2OHOX = (ZSCH18*ZMOLECH+ZCH1*ZMOLECOH)*ZMOLECHO2
    END FUNCTION ZFUNCPRODH2OHOX
    !
    FUNCTION ZFUNC3BODY (ZK1,ZK2)
      REAL(dp) :: ZFUNC3BODY
      REAL(dp), INTENT(IN) :: ZK1,ZK2
      ZFUNC3BODY=ZK1/(1._dp+ZK1/ZK2)/EXP(.5108256_dp/(1._dp+LOG10(ZK1/ZK2)**2._dp)) ! = *0.6^(...), so fc=0.6
    END FUNCTION ZFUNC3BODY

      END SUBROUTINE CHEMICS


  SUBROUTINE INRCGAS 
    !                                                             
    ! purpose:                                                    
    ! --------                                                    
    ! PRECALCULATES GASPHASE REACTION RATES FOR CHEMISTRY         
    ! IN A LOOKUPTABLE.                                           
    !                                                             
    ! method:                                                     
    ! -------                                                     
    ! THE REACTION RATES ARE CALCULATED IN THE                    
    ! TEMPERATURE INTERVAL IRCTMIN - IRCTMAX K IN 1./JDIFTE K STEPS.
    ! TEMPERATURES HIGHER (LOWER) THAN IRCTMAX (IRCTMIN) ARE      
    ! TREATED LIKE IRCTMAX (IRCTMIN).                             
    ! CURRENTLY USED WITH IRCTMIN=170.                            
    !                     IRCTMAX=320.                            
    !                     JVALS=4                                 
    !                                                             
    ! THE ACCURACY OF THE RATES IS THEN +-1/(2*JVALS) (=0.125) K! 
    !                                                             
    ! important parameters:                                       
    ! ---------------------                                       
    !     IRCTMIN: MINIMUM TEMPERATURE FOR WHICH RATES ARE        
    !             CALCULATED                                      
    !     IRCTMAX: MAXIMUM TEMPERATURE FOR WHICH RATES ARE        
    !             CALCULATED                                      
    !     JVALS : NUMBER OF INTERVALS IN WHICH THE DISTANCE       
    !             BETWEEN TWO INTEGER TEMPERATURES IS DEVIDED     
    !     JDIFTE: IRCTMAX-IRCTMIN                                 
    !     NUMTEM: TOTAL NUMBER OF TEMPERATURES ON WHICH THE       
    !             RATES ARE CALCULATED                            
    !     NUMRAT: NUMBER OF REACTION RATES AND                    
    !             PARTS OF REACTION RATES THAT ARE PRECALCULATED  
    !                                                             
    ! interface:                                                  
    ! ----------                                                  
    ! called from global_start                                         
    !                                                             
    ! externals:                                                  
    ! ----------                                                  
    ! none                                                        
    !                                                             
    ! written by B. Steil, 08-sep-1997                            
    ! last modified by R. Hein, 25-sep-1997                       
    !                                                             
    !                                                             
    !*    *PARAMETERS* CONTROLLING ARRAY SIZES.                   
    !                                                             
    !       ----------------------------------------------------------------
    !                      

    IMPLICIT NONE

    INTEGER, PARAMETER :: JPM = 106, JPGL = 160, JPNLON =  &
         320, JPNLEV = 39,     &
         JPTASKS = 16, &
         JPNLVP1 = &
         JPNLEV + 1
    !                                                             
    !                                                             
    !                                                             

    REAL(dp) :: ZTEMPST, ZTEMP, ZZT, ZZTLN
    INTEGER :: J

    INTRINSIC EXP, LOG, REAL

    ZTEMPST = 1._dp / REAL (JVALS) 
    DO J = 1, NUMTEM 
       ZTEMP = REAL (IRCTMIN) + REAL (J - 1) * ZTEMPST 
       !                                                             
       ! *************************************************************** 
       ! * TEMPERATURE DEPENDENT RATE CONSTANTS FOR GASPHASE-CHEMISTRY * 
       ! *                                                             * 
       ! *   ZCO1      :  O(3P) + O3     ----> 2 O2                    * 
       ! *   ZCO4      :  O(1D) + O2     ----> O(3P) + O2              * 
       ! *   ZCO5      :  O(1D) + N2     ----> O(3P) + N2              * 
       ! *   ZCN1      :  N + O2         ----> NO + O(3P)              * 
       ! *   ZCN6      :  HNO4 + OH      ----> NO2 + H2O + O2          * 
       ! *   ZCN8      :  NO + O3        ----> NO2 + O2                * 
       ! *   ZCN9      :  NO + ClO       ----> NO2 + Cl                * 
       ! *   ZCN10     :  NO + HO2       ----> NO2 + OH                * 
       ! *   ZCN11     :  NO + CH3O2     ----> NO2 + CH3O              * 
       ! *   ZCN12     :  NO + NO3       ----> 2 NO2                   * 
       ! *   ZCN13     :  NO2 + O(3P)    ----> NO + O2                 * 
       ! *   ZCN15     :  NO2 + O3       ----> NO3 + O2                * 
       ! *   ZCN20     :  ClONO2 +O(3P)  ----> ClO + NO3               * 
       ! *   ZCC1      :  ClO + O(3P)    ----> Cl + O2                 * 
       ! *   ZCC2      :  HCl + OH       ----> Cl + H2O                * 
       ! *   ZCC3      :  ClO + OH       ----> Cl + HO2                * 
       ! *   ZCC4      :  Cl + O3        ----> ClO + O2                * 
       ! *   ZCC5      :  Cl + H2        ----> HCl + H                 * 
       ! *   ZCC6      :  Cl + CH4 (+O2) ----> HCl + CH3O2             * 
       ! *   ZCC7      :  Cl + HO2       ----> HCl + O2                * 
       ! *   ZCC8      :  Cl + HO2       ----> ClO + OH                * 
       ! *   ZCC9      :  Cl + HCHO (+O2)----> HCl + CO + HO2          * 
       ! *   ZCC10     :  ClO + HO2      ----> HOCl + O2               * 
       ! *   ZCC11     :  HOCl + OH      ----> ClO + H2O               * 
       ! *   ZCC14     :  CH3Cl + OH     ----> CH2Cl + H2O             * 
       ! *   ZCC15     :  CH3CCl3 + OH   ----> CH2CCl3 + H2O           * 
       ! *   ZCH4      :  OH + O(3P)     ----> O2 + H                  * 
       ! *   ZCH6      :  H2 + OH        ----> H2O + H                 * 
       ! *   ZCH9      :  O3 + H         ----> OH + O2                 * 
       ! *   ZCH11     :  O(3P) + HO2    ----> OH + O2                 * 
       ! *   ZCH12     :  O3 + HO2       ----> OH + 2 O2               * 
       ! *   ZCH14     :  O3 + OH        ----> HO2 + O2                * 
       ! *   ZCH15     :  H2O2 + OH      ----> H2O + HO2               * 
       ! *   ZCH19     :  CH4 + OH (+O2) ----> CH3O2 + H2O             * 
       ! *   ZCH22     :  CH3O2 + HO2    ----> CH3O2H + O2             * 
       ! *   ZCH23     :  O2 + CH3O      ----> HCHO + HO2              * 
       ! *   ZCO3      :  O(3P) + O2 (+M)----> O3                      * 
       ! *   ZCN14     :  NO2 + ClO (+M) ----> ClONO2                  * 
       ! *   ZCN16     :  NO2 + HO2 (+M) ----> HNO4                    * 
       ! *   ZCN17     :  NO2 + OH (+M)  ----> HNO3                    * 
       ! *   ZCN18     :  NO2 + NO3 (+M) ----> N2O5                    * 
       ! *   ZCH13     :  H + O2 (+M)    ----> HO2                     * 
       ! *   ZCCLOCLO  :  2 ClO (+M)     ----> Cl2O2                   * 
       ! *   ZCN19     :  HNO3 + OH      ----> H2O + NO3               * 
       ! *   ZCH1      :  OH + HO2       ----> H2O + O2                * 
       ! *   ZCH2      :  2 HO2          ----> H2O2 + O2               * 
       ! *   ZCH5      :  CO + OH        ----> CO2 + H                 * 
       ! *                                                             * 
       ! * equilibrium constants:                                      * 
       ! *   ZCN5      : HNO4   / ZCN16                                * 
       ! *   ZCN7      : N2O5   / ZCN18                                * 
       ! *   ZCCL2O2   : Cl2O2  / ZCCLOCLO                             * 
       ! *                                                             * 
       ! *************************************************************** 
       !                                                             
       !     ARRHENIUS-TYPE                                          
       !                                                             
       ! ZCO1                                                        
       RCGAS (J, 1) = 8.0E-12 * EXP ( - 2060. / ZTEMP) 
       ! ZCO4                                                        
       RCGAS (J, 2) = 3.3E-11 * EXP (55. / ZTEMP) 
       ! ZCO5                                                        
       RCGAS (J, 3) = 2.15E-11 * EXP (110. / ZTEMP) 
       !                                                             
       ! ZCN1                                                        
       RCGAS (J, 4) = 1.5E-11 * EXP ( - 3600. / ZTEMP) 
       ! ZCN6                                                        
       RCGAS (J, 5) = 1.3E-12 * EXP (380. / ZTEMP) 
       ! ZCN8                                                        
       RCGAS (J, 6) = 3.0E-12 * EXP ( - 1500. / ZTEMP) 
       ! ZCN9                                                        
       RCGAS (J, 7) = 6.2E-12 * EXP (295. / ZTEMP) 
       ! ZCN10                                                       
       RCGAS (J, 8) = 3.5E-12 * EXP (250. / ZTEMP) 
       ! ZCN11                                                       
       RCGAS (J, 9) = 2.8E-12 * EXP (300. / ZTEMP) 
       ! ZCN12                                                       
       RCGAS (J, 10) = 1.5E-11 * EXP (170. / ZTEMP) 
       ! ZCN13                                                       
       RCGAS (J, 11) = 5.1E-12 * EXP (210. / ZTEMP) 
       !      RCGAS(J,11) = 6.5E-12*EXP(120./ZTEMP)                  
       ! ZCN15                                                       
       RCGAS (J, 12) = 1.2E-13 * EXP ( - 2450. / ZTEMP) 
       ! ZCN20                                                       
       RCGAS (J, 13) = 4.5E-12 * EXP ( - 900. / ZTEMP) 
       !                                                             
       ! ZCC1                                                        
       RCGAS (J, 14) = 2.5E-11 * EXP (110. / ZTEMP) 
       ! ZCC2                                                        
       RCGAS (J, 15) = 1.7E-12 * EXP ( - 230. / ZTEMP) 
       ! ZCC3                                                        
       RCGAS (J, 16) = 7.3E-12 * EXP (300. / ZTEMP) ! mz_ab_20091221  --> ... + HCl + O2 in MECCA
       ! ZCC4                                                        
       RCGAS (J, 17) = 2.8E-11 * EXP ( - 250. / ZTEMP) 
       ! ZCC5                                                        
       RCGAS (J, 18) = 3.9E-11 * EXP ( - 2310. / ZTEMP) 
       ! ZCC6                                                        
       RCGAS (J, 19) = 6.6E-12 * EXP ( - 1240. / ZTEMP) 
       ! ZCC7                                                        
       RCGAS (J, 20) = 4.4E-11-7.5E-11*EXP(-620./ZTEMP) 
       ! ZCC8                                                        
       RCGAS (J, 21) = 7.5E-11*EXP(-620./ZTEMP) 
       ! ZCC9                                                        
       RCGAS (J, 22) = 8.1E-11 * EXP ( - 34. / ZTEMP) 
       ! ZCC10                                                       
       RCGAS (J, 23) = 2.2E-12 * EXP (340. / ZTEMP) 
       ! ZCC11                                                       
       RCGAS (J, 24) = 3.0E-12 * EXP ( - 500. / ZTEMP) 
       ! ZCC14                                                       
       RCGAS (J, 25) = 2.4E-12 * EXP ( - 1250. / ZTEMP) 
       ! ZCC15                                                       
       RCGAS (J, 26) = 1.64E-12 * EXP ( - 1520. / ZTEMP) 
       !                                                             
       !     HELLEIS,CROWLEY,MOORTGAT MPI-MAINZ                      
       !                                                             
       !      ZCCH3O2CLO  = 3.25E-12*EXP(-114./ZTEMP)                
       !                                                             
       ! ZCH4                                                        
       RCGAS (J, 27) = 2.2E-11 * EXP (120. / ZTEMP) 
       ! ZCH6                                                        
       RCGAS (J, 28) = 2.8E-12 * EXP ( - 1800. / ZTEMP) 
       ! ZCH9                                                        
       RCGAS (J, 29) = 1.4E-10 * EXP ( - 470. / ZTEMP) 
       ! ZCH11                                                       
       RCGAS (J, 30) = 3.0E-11 * EXP (200. / ZTEMP) 
       ! ZCH12                                                       
       RCGAS (J, 31) = 1.0E-14 * EXP ( - 490. / ZTEMP) 
       ! ZCH14                                                       
       RCGAS (J, 32) = 1.7E-12 * EXP ( - 940. / ZTEMP) 
       ! ZCH15                                                       
       RCGAS (J, 33) = 1.8E-12_dp
       ! ZCH19                                                       
       RCGAS (J, 34) = 1.85E-20 * EXP (2.82*LOG(ZTEMP) -987./ZTEMP) 
       ! ZCH22                                                       
       RCGAS (J, 35) = 4.1E-13 * EXP (750. / ZTEMP) 
       ! ZCH23                                                       
       RCGAS (J, 36) = 3.9E-14 * EXP ( - 900. / ZTEMP) 
       !                                                             
       ! RATE CONSTANTS FOR THREE-BODY REACTIONS                     
       ! REF: DEMORE ET AL. (1992)                                   
       !                                                             
       ZZT = ZTEMP / 300. 
       ZZTLN = LOG (ZZT) 
       ! ZCO3                                                        
       RCGAS (J, 37) = 6.0E-34 / EXP (ZZTLN * 2.4) 
       ! ZCN14                                                       
       RCGAS (J, 38) = 1.8E-31 / EXP (ZZTLN * 3.4) 
       RCGAS (J, 39) = 1.5E-11 / EXP (ZZTLN * 1.9) 
       ! ZCN16                                                       
       RCGAS (J, 40) = 2.0E-31 / EXP (ZZTLN * 3.4) 
       RCGAS (J, 41) = 2.9E-12 / EXP (ZZTLN * 1.1) 
       ! ZCN17                                                       
       !     Brown et al., GRL, oct. 26,98                           
       RCGAS (J, 42) = 1.8E-30 / EXP (ZZTLN * 3.0) 
       RCGAS (J, 43) = 2.8E-11 / EXP (ZZTLN * 0.0) 
       !      RCGAS(J,42) = 2.6E-30/EXP(ZZTLN*3.2)                   
       !      RCGAS(J,43) = 2.4E-11/EXP(ZZTLN*1.3)                   
       ! ZCN18                                                       
       RCGAS (J, 44) = 2.0E-30 / EXP (ZZTLN * 4.4) 
       RCGAS (J, 45) = 1.4E-12 / EXP (ZZTLN * 0.7) 
       ! ZCH13                                                       
       RCGAS (J, 46) = 4.4E-32 / EXP (ZZTLN * 1.3) 
       RCGAS (J, 66) = 4.7E-11 / EXP (ZZTLN * 0.2) 
       ! ZCCLOCLO                                                    
       RCGAS (J, 47) = 1.6E-32 / EXP (ZZTLN * 4.5) 
       RCGAS (J, 61) = 2.0E-12 / EXP (ZZTLN * 2.4) 
       !                                                             
       !                                                             
       !S1      ZK1         = 1.5E-30/ZZT**4*ZCON                    
       !S1      ZK2         = 6.5E-12/ZZT**2                         
       !S1      ZCNO2CH3O2  = ZFUNC3BODY(ZK1,ZK2)                    
       !                                                             
       ! THERMAL DECOMPOSITION RATE CONSTANTS ; EQUILIBRIUM CONSTANTS
       !  FROM DEMORE ET AL. (1992)                                  
       !  (IN DEMORE GILT: HINREAKT./EQUIKON.=REAKTIONSRATE, IN DIE  
       !   M SCHON EINMULTIPLIZIERT IST !)                           
       !                                                             
       ! ZCN5 ( = ZCN16/RCGAS(J,48) )                                
       RCGAS (J, 48) = 2.1E-27 * EXP (10900. / ZTEMP) 
       !                                                             
       ! ZCN7 ( = ZCN18/RCGAS(J,49) )                                
       RCGAS (J, 49) = 2.7E-27 * EXP (11000. / ZTEMP) 
       !                                                             
       !      ZCCL2O2M    = ZCCLOCLO/ZCEQ                            
       RCGAS (J, 50) = 9.3E-28 * EXP (8835. / ZTEMP) 
       !                                                             
       !S1      ZCCH3O2NO2M = ZCNO2CH3O2/(1.3E-28*EXP(11200./ZTEMP)) 
       !                                                             
       !     BINARY REACTIONS WITH SPECIAL CALCULATIONS; DEMORE ET AL. (1992)
       !                                                             
       ! ZCN19  ( = RCGAS(J,51)+RCGAS(J,52)/(1+RCGAS(J,52)/RCGAS(J,53)) )
       !      RCGAS(J,51) = 7.2E-15*EXP(785./ZTEMP)                  
       !      RCGAS(J,52) = 1.9E-33*EXP(725./ZTEMP)                  
       !      RCGAS(J,53) = 4.1E-16*EXP(1440./ZTEMP)                 
       !     Brown et al., GRL, oct. 26,98                           
       RCGAS (J, 51) = 2.4E-14 * EXP (460. / ZTEMP) 
       RCGAS (J, 52) = 6.5E-34 * EXP (1335. / ZTEMP) 
       RCGAS (J, 53) = 2.7E-17 * EXP (2199. / ZTEMP) 
       !                                                             
       ! ZCH1                                                        
       RCGAS (J, 54) = 4.8E-11 * EXP (250. / ZTEMP) 
       !                                                             
       ! ZCH2 ( = RCGAS(J,55)+ZCON*RCGAS(J,56) )                     
       RCGAS (J, 55) = 1.5E-12 * EXP (19. / ZTEMP) 
       RCGAS (J, 56) = 1.7E-33 * EXP (1000. / ZTEMP) 
       RCGAS (J, 62) = 1.4E-21 * EXP (2200. / ZTEMP) 
       !      ZCH5        = 1.5E-13*(1+0.6*ZPRES/101400.)            
       !                                                             
       !     ZCN2                                                    
       RCGAS (J, 57) = 2.1E-11 * EXP (100. / ZTEMP) 
       !     ZCN23                                                   
       RCGAS (J, 58) = 5.8E-12 * EXP (220. / ZTEMP) 
       !     ZCH20                                                   
       RCGAS (J, 59) = 2.66E-12 * EXP (200. / ZTEMP) 
       !     ZCH21                                                   
       RCGAS (J, 60) = 1.14E-12 * EXP (200. / ZTEMP) 
 
       ! mz_ab_20091221+
       ! ZCN21
       RCGAS (J, 63) = 6.7E-11 * EXP (20. / ZTEMP) 
       ! ZCN22
       RCGAS (J, 64) = 4.7E-11 * EXP (20. / ZTEMP) 
       ! ZCH8
       RCGAS (J, 65) = 1.63E-10 * EXP (60. / ZTEMP) 
       ! ZCH16
       RCGAS (J, 66) = 9.52E-18 * EXP (2.03*LOG(ZTEMP) + 636./ ZTEMP) 
      
       ! mz_ab_20091221-
      !                                                             
    ENDDO
    !                                                             
    RETURN 
  END SUBROUTINE INRCGAS

  !********************************************************************

! op_pj_20101216+
!!$  SUBROUTINE inisulnew(latitude, ilat, nproma, npromz, ngpblks)
  SUBROUTINE inisulnew(latitude, nproma, ngpblks, npromz, vnpromz)
! op_pj_20101216-

    IMPLICIT NONE

! op_pj_20101216+
!!$    REAL(dp), DIMENSION(:),   INTENT(IN) :: latitude ! [deg]
!!$    INTEGER,  DIMENSION(:,:), INTENT(IN) :: ilat
    REAL(dp), DIMENSION(:,:),   INTENT(IN) :: latitude ! [deg]
! op_pj_20101216-
    INTEGER, INTENT(IN)           :: nproma, ngpblks 
    INTEGER, INTENT(IN), OPTIONAL :: npromz
    INTEGER, INTENT(IN), OPTIONAL :: vnpromz(:)

    ! LWMO    = # of latitudes in WMO aerosol climatology
    ! KWMO    = # of levels    in WMO aerosol climatology
    ! KSUL    = lowest level of   WMO aerosol climatology [hPa]
    ! RBASE   = conversion factor for units (10^-8cm^2/cm^3 ---> cm^-1)
    !           and downscaling from WMO median scenario to baseline
    !           scenario (0.25)
    ! RSULUP/DO = upper/lower boundary for het. chem. on sulfor aerosols
    ! RNATUP/DO = upper/lower boundary for het. chem. on NAT
    ! RICEUP/DO = upper/lower boundary for het. chem. on ICE
    ! SHEIWMO(K)            HEIGHTS
    ! SULWMO(K,JL=1:LWMO)   AEROSOL SURFACE AREA

!!$AEROSOL SURFACE AREA (10^-8 CM^2CM^-3)
!!$Z^*[KM] = 16*LOG10(1000/P[HPA])
!!$JAHRESMITTEL f\"ur 5 BREITENGUERTEL
!!$-------------------------------------------------------
!!$Z^* | 60-90 N | 30-60 N | 30N-30S | 30-60 S | 60-90 S |
!!$    |   75    |   45    |    0    |  -45    |  -75    |
!!$-------------------------------------------------------
!!$32  | 0.10    | 0.10    | 0.40    | 0.25    | 0.10    |
!!$30  | 0.20    | 0.30    | 0.70    | 0.50    | 0.20    |
!!$28  | 0.50    | 0.70    | 1.30    | 0.90    | 0.50    |
!!$26  | 1.00    | 1.00    | 1.70    | 1.35    | 1.00    |
!!$24  | 1.40    | 1.50    | 2.00    | 1.85    | 1.40    |
!!$22  | 1.80    | 2.00    | 2.50    | 2.50    | 2.00    |
!!$20  | 2.50    | 2.50    | 3.00    | 3.10    | 2.75    |
!!$18  | 3.00    | 3.00    | 3.00    | 4.00    | 4.00    |
!!$16  | 3.25    | 3.25    | 2.00    | 4.00    | 4.75    |
!!$14  | 4.00    | 4.00    | 2.00    | 4.50    | 5.50    |
!!$12  | 4.50    | 4.50    | 2.00    | 4.50    | 6.00    |


    INTEGER,  PARAMETER :: kwmo=11, lwmo=5
    REAL(dp), PARAMETER :: rbase=2.5e-9_dp
    REAL(dp), PARAMETER :: rnatup=2000._dp

    INTEGER, PARAMETER :: k0=10
    INTEGER, PARAMETER :: k0h=k0/2

    REAL(dp) :: PR, BB, ZSTR,  AA, DD, ZLAT
    INTEGER :: JL, K, JLL, INMAX, INMIN, jrow

    REAL(dp), DIMENSION(KSUL,LWMO) ::  SHELP1
    REAL(dp), DIMENSION(KWMO,LWMO) ::  SULWMO
    REAL(dp), DIMENSION(KWMO)      ::  SHEIWMO
    REAL(dp), DIMENSION(LWMO)      ::  SLATWMO

    INTEGER :: zkproma

    INTRINSIC INT, LOG10, REAL

    ! Last height is read first, so order changed to above file
    DATA SULWMO /  4.50, 4.00, 3.25, 3.00, 2.50, 1.80, 1.40, 1.00, 0.50, 0.20, 0.10, &
         4.50, 4.00, 3.25, 3.00, 2.50, 2.00, 1.50, 1.00, 0.70, 0.30, 0.10, &
         2.00, 2.00, 2.00, 3.00, 3.00, 2.50, 2.00, 1.70, 1.30, 0.70, 0.40, &
         4.50, 4.50, 4.00, 4.00, 3.10, 2.50, 1.85, 1.35, 0.90, 0.50, 0.25, &
         6.00, 5.50, 4.75, 4.00, 2.75, 2.00, 1.40, 1.00, 0.50, 0.20, 0.10/

    DATA SLATWMO / 75.,   45.,   0.,  -45.,  -75. /
    ! Last height is read first, so order changed to above file
    DATA SHEIWMO / 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 32. /

    !     BASELINE IS REQUIRED!
    SULWMO(:,:) = RBASE*SULWMO(:,:)

    !     FOR LOOKUP-TABLE CALCULATION PERFORM FIRST
    !     THE INTERPOLATION ON 1 HPA STEPS!

    SHELP1(:,:) = 0.

    DD    = SHEIWMO(KWMO)-SHEIWMO(KWMO-1)

    DO JL=1,LWMO
       DO K=k0,KSUL
          PR           = REAL(K,dp)
          ZSTR         = 16.*LOG10(1000./PR)-0.01_dp
          INMAX        = INT(ZSTR/dd)-k0h+1
          INMIN        = INMAX-1
          BB           = SULWMO(INMAX,JL)
          AA           = (SULWMO(INMIN,JL)-BB)/DD
          SHELP1(K,JL) = AA*(REAL(INMAX+k0h)*dd-ZSTR+0.01_dp)+BB
       ENDDO
    ENDDO

    !     PERFORM LATITUDE INTERPOLATION (PINGPONG) !
    SULOOK(:,:,:) = 0._dp

    DO jrow=1,ngpblks

       zkproma = -1
       IF (PRESENT(npromz)) THEN
          IF ( jrow == ngpblks ) THEN
             zkproma = npromz
          ELSE
             zkproma = nproma
          END IF
       ENDIF
       IF (PRESENT(vnpromz)) THEN
          zkproma = vnpromz(jrow)
       END IF

       DO JL=1,zkproma
! op_pj_20101216+
!!$          ZLAT=latitude(ilat(JL,jrow)) ! mz_ab_20090618
          ZLAT=latitude(JL,jrow) ! mz_ab_20090618
! op_pj_20101216-

          IF (ZLAT.GE.SLATWMO(1)) THEN
             DO K=k0,KSUL
                SULOOK(K,JL,jrow) = SHELP1(K,1)
             ENDDO
          ELSE

             IF (ZLAT.LE.SLATWMO(LWMO)) THEN
                DO K=k0,KSUL
                   SULOOK(K,JL,jrow) = SHELP1(K,LWMO)
                ENDDO
             ELSE

                INMAX = LWMO
                DO JLL=LWMO-1,2,-1
                   IF (ZLAT.ge.SLATWMO(JLL)) INMAX = INMAX - 1
                ENDDO
                INMIN = INMAX - 1
                dd=slatwmo(inmin)-slatwmo(inmax)

                DO K=k0,KSUL
                   BB           = SHELP1(K,INMAX)
                   AA           = (SHELP1(K,INMIN)-SHELP1(K,INMAX))/DD
                   SULOOK(K,JL,jrow) = AA*(zlat-slatwmo(inmax))+BB
                ENDDO
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    return

  END SUBROUTINE inisulnew

  !===================================================================

  SUBROUTINE CLSCAV (Conc,PTMST,PRHOA,PDP                             &
       ,zmratep, zfprec, zfevap, zclcover, zmlwc                  &
       ,cvdprec,kconbot,PT,NLON,NLEV                  &
       )                            
!!$  SUBROUTINE CLSCAV (Conc,PTMST,PRHOA,PDP                             &
!!$       ,zmratep, zfprec, zfevap, zclcover, zmlwc                  &
!!$       ,xn3depcv,cvdprec,kconbot,PT,NLON,NLEV,jrow                  &
!!$       )                            

  ! INPUT:
  ! Conc(NLON,NLEV,NSPEC) = Concentrations
  ! PTMST = time_step_len
  ! PRHOA = air density in g cm-3          
  ! PDP(JL, JK) = APHP1 (JL, JK + 1) - APHP1 (JL, JK) 
  ! zmratep(nlon,nlev) = PRECIPITATION FORMATION RATE [KG/(M3*SEC)] 
  ! zfprec(nlon,nlev) = precipitation flux
  ! zfevap(nlon,nlev) = evaporation flux
  ! zclcover(nlon,nlev) = CLOUD COVER
  ! zmlwc(nlon,nlev) =  LIQUID WATER CONTENT [KGH2O/KGAIR] 
  ! xn3depcv(nlon,nlev,3) = deposition HNO3, HCl ???
  ! cvdprec(nlon,nlev) =
  ! kconbot(nlon) =
  ! PT(nlon,nlev) = temperature

  ! OUTPUT:
  ! Conc

    USE messy_main_constants_mem, ONLY : dp, &
         zavogadro => N_A,                   & ! Avogadro constant / (1/mol)
         zm_air => M_air,                    & ! molar mass of dry air [g/mol]
         MO, MC,                             &
         ZGASC => R_gas, g,  &
         api => pi

    IMPLICIT NONE

    INTEGER,  INTENT(IN)    :: NLON,NLEV
    REAL(dp), INTENT(IN)    :: PTMST
    REAL(dp), INTENT(INOUT) :: Conc(NLON,NLEV,NSPEC)
    REAL(dp), INTENT(IN)    :: PRHOA(NLON,NLEV), PDP(NLON,NLEV)                                          
    REAL(dp), INTENT(IN)    :: PT(NLON,NLEV) 
    REAL(dp), INTENT(IN)    :: zmratep(nlon,nlev),zfprec(nlon,nlev)  &
         , zfevap(nlon,nlev),zclcover(nlon,nlev), zmlwc(nlon,nlev)   &
         , cvdprec(nlon,nlev) 
!!$    REAL(dp), INTENT(IN)    :: xn3depcv(nlon,nlev,3)
    REAL(dp), INTENT(IN)    :: kconbot(nlon) 

    REAL(dp), PARAMETER     :: zmoltom2s=zavogadro*1000./zm_air

    REAL(dp) :: WDHNO3(NLON),WDHCL(NLON),OLDCOV(NLON),               &
         KCLTOP(NLON)                                                 
    REAL(dp) :: WDHNO3CV(NLON),WDHCLCV(NLON) 
    REAL(dp) :: RAINCV(NLON) 
    REAL(dp) :: ZMTOF(NLON,NLEV) 

    REAL(dp) :: DGHNO3=0.136_dp &
         , DGHCL=0.136_dp          &
         , DGH2O2=0.184_dp         &
         , DGAIR=0.133_dp 

    REAL(dp) :: zeop,DHNO3L,DHClL,DH2O2L,DHNO3,DHCl,washfrac,dwd,HNO3L,HClL,zmxcov &
         , zclear,fac1,H2O2L,cvcover,henper,PERFRACL, znre, znsc, rn, ru &
         , beta, betaex, znsh, zkg, ZTOTDEP4, ZTOTDEP5 &
         , rlwc, rdrad, rdvol, cl2rain, zinflux, rflx


    INTEGER :: jk, jl

    INTRINSIC SQRT

    WDHNO3(1:NLON)=0._dp
    WDHCL(1:NLON)=0._dp
    OLDCOV(1:NLON)=0._dp
    WDHNO3CV(1:NLON)=0._dp
    WDHCLCV(1:NLON)=0._dp
    RAINCV(1:NLON)=0._dp 
    KCLTOP(1:NLON)=REAL(nlev,dp)

    ! ******************************************************************    
    ! *                                                                *    
    ! * wet dep. due to stratiform precipitation                       *    
    ! *                                                                *    
    ! ******************************************************************    
    !                                                                       
    !   Coefficients for sub-cloud scavenging:                              
    !                                                                       
    !                                                                       
    DO JK=1,NLEV 
       DO JL=1,NLON 
          ! - ZMTOF: mass-to-flux (cm3 m-2 s-1)                                   
          ZMTOF(JL,JK)=PDP(JL,JK)/(PTMST*G*PRHOA(JL,JK)*1.E-3) 
       ENDDO
    ENDDO
    ZTOTDEP5=0._dp
    ZTOTDEP4=0._dp 
    DO JK=1,NLEV 
       !                                                                       
       !   Below cloud scavenging:                                             
       !    adds to total wet deposition flux                                  
       !   Note that ZFPREC includes the additional flux in the cloud!         
       !                                                                       
       DO JL=1,NLON 
          CL2RAIN=ZMRATEP(JL,JK)*PTMST*ZMTOF(JL,JK)*1.E-6 
          ZINFLUX=MAX(ZFPREC(JL,JK)-CL2RAIN,0._dp) 
          ! CAUTION: explicit level number!!! ??????????????????????????????????  
          IF (ZINFLUX.LT.1E-15) KCLTOP(JL)=REAL(nlev,dp)
          IF ((ZINFLUX-ZFEVAP(JL,JK)).GT.1E-15.AND.REAL(JK,dp).GT.KCLTOP(JL)) THEN 
             ! -- calculate rain parameters (Kumar, 1985)                            
             RFLX=ZINFLUX/OLDCOV(JL)*3600._dp
             RLWC=72._dp*RFLX**0.88_dp
             RDRAD=0.3659_dp*RFLX**0.21_dp
             RDVOL=4._dp*api/3._dp*RDRAD**3._dp
             RN=RLWC/RDVOL 
             RU=9.58_dp*(1._dp-EXP(-(RDRAD/0.885)**1.147_dp)) 
             ! -- calculate mass transfer coefficient HNO3 (Durham et. al., 1981)    
             ZNRE=SQRT(20._dp*RDRAD*RU/DGAIR) 
             ZNSC=DGAIR/DGHNO3 
             ZNSH=1._dp+0.3_dp*ZNRE*(ZNSC**(1._dp/3._dp)) 
             ZKG=10._dp*DGHNO3/RDRAD*ZNSH 
             ! -- calculate washout (Kumar, 1985); account for cloud overlap         
             BETA=4._dp*api*RDRAD**2._dp*RN*ZKG*1.E-8 
             BETAEX = EXP(-1._dp*BETA*PTMST)-1._dp 
             DHNO3=Conc(JL,JK,ind_HNO3)*BETAEX 
             DHCL =Conc(JL,JK,ind_HCl)*BETAEX 
             ZMXCOV=MIN(OLDCOV(JL),ZCLCOVER(JL,JK)) 
             WASHFRAC=OLDCOV(JL)-ZMXCOV 
             DHNO3=DHNO3*WASHFRAC 
             DHCL =DHCL *WASHFRAC 
             ! NO WET LARGESCALE BELOW CLOUD SCAVENGING                              
             !         DHNO3 = 0._dp                                                  
             !         DHCL = 0._dp                                                   
             Conc(JL,JK,ind_HNO3)=Conc(JL,JK,ind_HNO3)+DHNO3 
             Conc(JL,JK,ind_HCl)=Conc(JL,JK,ind_HCl)+DHCL 
             DWD=DHNO3*ZMTOF(JL,JK) 
             WDHNO3(JL)=WDHNO3(JL)-DWD 
             DWD=DHCL*ZMTOF(JL,JK) 
             WDHCL(JL)=WDHCL(JL)-DWD 
             !        ZTOTDEP5=ZTOTDEP5-DHNO3*GRVOL(JL,JK)                           
             !        ZTOTDEP4=ZTOTDEP4-DHCL*GRVOL(JL,JK)                            
          ENDIF
       ENDDO
       ! - - -                                                                 
       !   Re-evaporation of incoming flux                                     
       !                                                                       
       DO JL=1,NLON 
          CL2RAIN=ZMRATEP(JL,JK)*PTMST*ZMTOF(JL,JK)*1.E-6 
          ZINFLUX=MAX(ZFPREC(JL,JK)-CL2RAIN,0._dp) 
          ZCLEAR=1._dp-ZCLCOVER(JL,JK) 
          IF (ZCLEAR.GT.1E-20.AND.ZINFLUX.GT.1E-20) THEN 
             ZEOP=ZFEVAP(JL,JK)/ZINFLUX 
             ZEOP=MAX(0._dp,ZEOP) 
             ZEOP=MIN(1._dp,ZEOP) 
             DHNO3=ZEOP*WDHNO3(JL)/ZMTOF(JL,JK) 
             DHCL =ZEOP*WDHCL(JL)/ZMTOF(JL,JK) 
             WDHNO3(JL)=(1._dp-ZEOP)*WDHNO3(JL) 
             WDHCL(JL) =(1._dp-ZEOP)*WDHCL(JL) 
             ! NO WET LARGESCALE REEVAPORATION OF INCOMING FLUX                      
             !         DHNO3 = 0.                                                    
             !         DHCL = 0.                                                     
             Conc(JL,JK,ind_HNO3)= Conc(JL,JK,ind_HNO3)+DHNO3 
             Conc(JL,JK,ind_HCl)=Conc(JL,JK,ind_HCl)+DHCL 
             !        ZTOTDEP5=ZTOTDEP5-DHNO3*GRVOL(JL,JK)                           
             !        ZTOTDEP4=ZTOTDEP4-DHCL*GRVOL(JL,JK)                            
          ENDIF
       END DO
       ! - - -                                                                 
       !   Transfer from cloud to rain water:                                  
       !    adds to wet deposition flux                                        
       !   (for H2O2 this is the only process considered)                      
       !                                                                       
       DO JL=1,NLON 
          IF (ZMRATEP(JL,JK).GT.1E-20) THEN 
             IF(REAL(JK,dp).LT.KCLTOP(JL)) KCLTOP(JL)=REAL(JK,dp)
             HNO3L=ZCLCOVER(JL,JK)*Conc(JL,JK,ind_HNO3) 
             HCLL =ZCLCOVER(JL,JK)*Conc(JL,JK,ind_HCl) 
             FAC1=HNO3L*ZMRATEP(JL,JK)/(PRHOA(JL,JK)*ZMLWC(JL,JK)) 
             DHNO3L=FAC1*PTMST*1.E-3 
             FAC1=HCLL*ZMRATEP(JL,JK)/(PRHOA(JL,JK)*ZMLWC(JL,JK)) 
             DHCLL =FAC1*PTMST*1.E-3 
             DHNO3L=MAX(0._dp,DHNO3L) 
             DHCLL =MAX(0._dp,DHCLL) 
             DHNO3L=MIN(Conc(JL,JK,ind_HNO3),DHNO3L) 
             DHCLL =MIN(Conc(JL,JK,ind_HCl),DHCLL) 
             WDHNO3(JL)=WDHNO3(JL)+DHNO3L*ZMTOF(JL,JK) 
             WDHCL(JL) =WDHCL(JL) +DHCLL*ZMTOF(JL,JK) 
             ! NO WET LARGESCALE TRANSFER FROM CLOUD TO RAINWATER                    
             !         DHNO3L = 0.                                                   
             !         DHCLL = 0.                                                    
             Conc(JL,JK,ind_HNO3)=Conc(JL,JK,ind_HNO3)-DHNO3L 
             Conc(JL,JK,ind_HCl)=Conc(JL,JK,ind_HCl)-DHCLL 
             !        ZTOTDEP5=ZTOTDEP5+DHNO3L*GRVOL(JL,JK)                          
             !        ZTOTDEP4=ZTOTDEP4+DHCLL*GRVOL(JL,JK)                           
             !     H2O2                                                              
             RLWC=ZMLWC(JL,JK)*PRHOA(JL,JK)*1.E9 
             HENPER=5.49E-7*EXP(-6620*(1./PT(JL,JK)-1./298.)) 
             PERFRACL=1.E-3/(HENPER/(RLWC*1E-6)+1.E-3) 
             H2O2L=ZCLCOVER(JL,JK)*PERFRACL*Conc(JL,JK,ind_H2O2) 
             DH2O2L=H2O2L*FAC1*PTMST*1.E-3 
             DH2O2L=MAX(0._dp,DH2O2L) 
             DH2O2L=MIN(Conc(JL,JK,ind_H2O2),DH2O2L) 
             Conc(JL,JK,ind_H2O2)=Conc(JL,JK,ind_H2O2)-DH2O2L 
             !                                                                       
             ! --- Set cover for gridbox below                                       
             OLDCOV(JL)=ZCLCOVER(JL,JK) 
          ENDIF
       END DO
       !ONSCAV  \|/                                                            
       !3) to be inserted in SUBROUTINE CLSCAV:                                
       !  113 CONTINUE              <-- This is already in CLSCAV !!           
       ! - - -                     <-- Insert from here .......                
       ! ******************************************************************    
       ! *                                                                *    
       ! * wet dep. due to convective precipitation                       *    
       ! *                                                                *    
       ! ******************************************************************    
       !  Below cloud scavenging                                               
       !                                                                       
       !ONSCAV  ECHAM-assumed convective cloud-cover                           
       CVCOVER=0.05_dp
       !                                                                       
       DO JL=1,NLON 
          IF (KCONBOT(JL).GT.0) THEN 
             IF (JK.GT.KCONBOT(JL).AND.RAINCV(JL).GT.0._dp) THEN 
                RFLX=RAINCV(JL)/CVCOVER*3600._dp
                RLWC=72._dp*RFLX**0.88_dp
                RDRAD=0.3659_dp*RFLX**0.21_dp
                RDVOL=4._dp*api/3*RDRAD**3._dp
                RN=RLWC/RDVOL 
                RU=9.58*(1.-EXP(-(RDRAD/0.885_dp)**1.147_dp)) 
                ZNRE=SQRT(20*RDRAD*RU/DGAIR) 
                ZNSC=DGAIR/DGHNO3 
                ZNSH=1._dp+0.3_dp*ZNRE*(ZNSC**(1._dp/3._dp)) 
                ZKG=10._dp*DGHNO3/RDRAD*ZNSH 
                BETA=4._dp*api*RDRAD**2._dp*RN*ZKG*1.E-8 
                BETAEX=EXP(-1._dp*BETA*PTMST)-1._dp 
                DHNO3=Conc(JL,JK,ind_HNO3)*BETAEX 
                DHNO3=DHNO3*CVCOVER 
                ! NO WET CONVECTIVE BELOW CLOUD SCAVENGING                              
                !         DHNO3 = 0.                                                    
                Conc(JL,JK,ind_HNO3)=Conc(JL,JK,ind_HNO3)+DHNO3 
                DHCL=Conc(JL,JK,ind_HCl)*BETAEX 
                DHCL=DHCL*CVCOVER 
                ! NO WET CONVECTIVE BELOW CLOUD SCAVENGING                              
                !         DHCL = 0.                                                     
                Conc(JL,JK,ind_HCl)=Conc(JL,JK,ind_HCl)+DHCL 
                DWD=DHNO3*ZMTOF(JL,JK) 
                WDHNO3CV(JL)=WDHNO3CV(JL)-DWD 
                DWD=DHCL*ZMTOF(JL,JK) 
                WDHCLCV(JL)=WDHCLCV(JL)-DWD 
             ENDIF
          ENDIF
       END DO
       ! - - -                                                                 
       !  Below cloud evaporation                                              
       !                                                                       
       DO  JL=1,NLON 
          IF (KCONBOT(JL).GT.0.AND.RAINCV(JL).GT.0.) THEN 
             IF (JK.GT.KCONBOT(JL).AND.CVDPREC(JL,JK).LT.0._dp) THEN 
                ZEOP=-1._dp*CVDPREC(JL,JK)/RAINCV(JL) 
                ZEOP=MAX(0._dp,ZEOP) 
                ZEOP=MIN(1._dp,ZEOP) 
                DHNO3=ZEOP*WDHNO3CV(JL)/ZMTOF(JL,JK) 
                ! NO WET CONVECTIVE BELOW CLOUD EVAPORATION                             
                !         DHNO3 = 0.                                                    
                WDHNO3CV(JL)=(1.-ZEOP)*WDHNO3CV(JL) 
                Conc(JL,JK,ind_HNO3)=Conc(JL,JK,ind_HNO3)+DHNO3 
                DHCL=ZEOP*WDHCLCV(JL)/ZMTOF(JL,JK) 
                ! NO WET CONVECTIVE BELOW CLOUD EVAPORATION                             
                !         DHCL = 0.                                                     
                WDHCLCV(JL)=(1.-ZEOP)*WDHCLCV(JL) 
                Conc(JL,JK,ind_HCl)=Conc(JL,JK,ind_HCl)+DHCL 
             ENDIF
          ENDIF
       END DO
       ! - - -                                                                 
!!$       !  Add up the in-cloud scavenging calculated in Tiedtke                 
!!$       !                                                                       
!!$       DO  JL=1,NLON 
!!$          IF (KCONBOT(JL).GT.0) THEN 
!!$             RAINCV(JL)=RAINCV(JL)+CVDPREC(JL,JK) 
!!$             IF (JK.LE.KCONBOT(JL)) THEN 
!!$                !         WDHNO3CV(JL)=WDHNO3CV(JL)+XN3DEPCV(JL,JK,1)                   
!!$                ! WDHNO3CV [MOLEC/(M2*SEC)]                                             
!!$                WDHNO3CV(JL)=WDHNO3CV(JL)+XN3DEPCV(JL,JK,1)                    &
!!$                     *ZMOLTOM2S                                        
!!$                WDHCLCV(JL)=WDHCLCV(JL)+XN3DEPCV(JL,JK,2)                      &
!!$                     *ZMOLTOM2S                                        
!!$             ENDIF
!!$          ENDIF
!!$       END DO
       !  110 CONTINUE   <--  already in CLSCAV                                
       !ONSCAV  /|                                                            
    ENDDO
    RETURN 
  END SUBROUTINE CLSCAV

  SUBROUTINE e4chem_read_nml_ctrl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/  NSTCHPH, l_fastscav, l_Brparam


    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='e4chem_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    WRITE(*,*) 'NUMBER OF OUTER STEPS PER CHEMISTRY STEP : ', NSTCHPH
    WRITE(*,*) 'FAST SCAVENGING                          : ', l_fastscav


    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE e4chem_read_nml_ctrl
  ! =========================================================================

  ! ***********************************************************************
END MODULE messy_e4chem
! ***********************************************************************
