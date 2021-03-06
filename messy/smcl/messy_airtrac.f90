! **********************************************************************
! AIRTRAC module: air traffic emissions can be tracked
! simplified O3-NOx-CH4-chemistry is calculated on trajectories
! Author : Christine Froemming, Volker Grewe, Sabine Brinkop, DLR, 2011/2012
!          Patrick Joeckel, DLR, 2014
! **********************************************************************
!
! Reference:
!
! * Frömming, C., Grewe, V., Brinkop, S., and Jöckel, P.: Documentation of the
!   EMAC submodels AIRTRAC 1.0 and CONTRAIL 1.0, supplementary material of
!   Grewe et al., 2014a, Geoscientific Model Development, 7, 175–201,
!   https://doi.org/10.5194/gmd-7-175-2014, 2014.
! **********************************************************************

MODULE messy_airtrac

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE
  SAVE
  PUBLIC :: DP

  ! ----------- <

  ! AIRTRAC SMCL ROUTINES (MAIN ENTRY POINTS)
  PUBLIC :: airtrac_read_nml_ctrl      ! read CTRL namelist and initialize
  PUBLIC :: airtrac_integrate

  ! GLOBAL PARAMETER
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'airtrac' ! name of module
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '13.02.2014' ! name of module

  ! GLOBAL NAMELIST SWITCH
   INTEGER,PUBLIC  :: n_emis_points = 0

! ----------------------------------------------------------------------

CONTAINS

! ************************************************************************
! AIRTRAC CORE ROUTINES
! ************************************************************************

! ----------------------------------------------------------------------
  SUBROUTINE airtrac_read_nml_ctrl(status, iou)

    ! airtrac MODULE ROUTINE (SMCL)
    ! READ airtrac NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    ! Author: Christine Froemming, DLR, 2010/2011

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    NAMELIST /CTRL/  n_emis_points
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'airtrac_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! ### ADD HERE DIAGNOSTIC OUPUT FOR LOG-FILE
    WRITE(*,*) 'Number of Emission points : ',n_emis_points


    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE airtrac_read_nml_ctrl
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  ELEMENTAL  SUBROUTINE airtrac_integrate(dt & ! time step length
       , o3_lg                         &  ! O3 background
       , no_lg,no2_lg                  &  ! NO background,NO2 background
       , o3prodn_lg_te                 &  ! ProdO3 tendency
       , o3lossn_lg_te                 &  ! LossO3N tendency
       , o3lossy_lg_te                 &  ! LossO3Y tendency
       , noxloss_lg_te                 &  ! LossNOx tendency
       , hno3loss_lg_te                &  ! LossHNO3 tendency
       , hno3_lg                       &  ! HNO3 background
       , hno3scav_lg_te                &  ! HNO3 Scavenging tendency
       , oh_lg                         &  ! OH
       , ho2_lg                        &  ! HO2
       , ohprod1_lg_te, ohprod2_lg_te,ohprod3_lg_te    & ! ProdOH1,2,3
       , ohloss1_lg_te, ohloss2_lg_te,ohloss3_lg_te    & ! LossOH1,2,3
       , ohloss4_lg_te, ohloss5_lg_te                  & ! LossOH4,5
       , ho2prod1_lg_te, ho2loss1_lg_te,ho2loss2_lg_te & ! ProdHO21,2,3
       , nox_air_lg,o3_air_lg          &  ! airNOx,airO3
       , o3prodn_air_lg,o3lossn_air_lg &  ! airProdO3N,airLossO3N
       , o3lossy_air_lg,noxloss_air_lg &  ! airLossO3Y,airLossNOx
       , hno3_air_lg,oh_air_lg         &  ! airHNO3,airOH
       , ho2_air_lg                    &  ! airHO2
       , nox_air_lg_te,o3_air_lg_te          &  ! airNOx_te,airO3_te
       , o3prodn_air_lg_te,o3lossn_air_lg_te &  ! airProdO3_te,airLossO3N_te
       , o3lossy_air_lg_te,noxloss_air_lg_te &  ! airLossO3Y_te,airLossNOx_te
       , hno3_air_lg_te,oh_air_lg_te         &  ! airHNO3_te,airOH_te
       , ho2_air_lg_te,ch4_air_lg_te         &  ! airHO2_te,airCH4_te
       , hno3scav_air_lg_te, hno3loss_air_lg_te &
       , h2o_lg                               &  ! H2O background
       , h2o_air_lg                           &  ! airH2O
       , h2o_loss_lg                          &  ! H2O Loss
       , h2o_air_lg_te )                         ! H2O air tendency

       ! airtrac MODULE ROUTINE
       ! NOxneu = NOxalt - L * dt
       ! d/dt O3  = P - D*O3
       ! integration method: euler backward
       ! O3neu = (O3alt + P * dt) / (1 + L/O3alt*dt)
       !

    IMPLICIT NONE
    INTRINSIC :: ABS, MIN

    ! I/O
    REAL(DP), INTENT(IN) :: dt  ! time step in seconds
    ! background tracer concentrations from messy (converted to lagrangian tracers)
    ! units: mol/mol  (tendencies in mol/mol/s)
    REAL(DP), INTENT(IN) :: no_lg             ! background no
    REAL(DP), INTENT(IN) :: no2_lg            ! background no2
    REAL(DP), INTENT(IN) :: o3_lg             ! o3 background
    REAL(DP), INTENT(IN) :: o3prodn_lg_te     ! o3prod tendency
    REAL(DP), INTENT(IN) :: o3lossn_lg_te     ! o3loss tendency via n
    REAL(DP), INTENT(IN) :: o3lossy_lg_te     ! o3loss tendency via other reactions
    REAL(DP), INTENT(IN) :: noxloss_lg_te     ! noxloss tendency via hno3 formation
    REAL(DP), INTENT(IN) :: hno3loss_lg_te    ! hno3loss tendency return to nox via NO3 or photolysis
    REAL(DP), INTENT(IN) :: hno3_lg           ! hno3 background
    REAL(DP), INTENT(IN) :: hno3scav_lg_te    ! hno3 tendency from scavenging
    ! airtrac new tracers: background values
    REAL(DP), INTENT(INOUT) :: nox_air_lg        ! new nox related to nox emission
    REAL(DP)                :: no2_air_lg        ! additional no2
    REAL(DP), INTENT(INOUT) :: o3_air_lg         ! new o3 related to nox emission
    REAL(DP), INTENT(IN)    :: o3prodn_air_lg    ! new o3prod related to nox emission
    REAL(DP), INTENT(IN)    :: o3lossn_air_lg    ! new o3lossn related to nox emission
    REAL(DP), INTENT(IN)    :: o3lossy_air_lg    ! new o3lossy related to nox emission
    REAL(DP), INTENT(IN)    :: noxloss_air_lg    ! noxloss via hno3 formation
    REAL(DP), INTENT(INOUT) :: hno3_air_lg       ! hno3_air on trajectory
    ! airtrac new tracer tendencies
    REAL(DP), INTENT(OUT) :: nox_air_lg_te     ! nox due to nox additional emission
    REAL(DP), INTENT(OUT) :: o3_air_lg_te      ! o3 tendency due to nox additional emission
    REAL(DP), INTENT(OUT) :: o3prodn_air_lg_te ! o3prod tendency due to nox additional emission
    REAL(DP), INTENT(OUT) :: o3lossn_air_lg_te ! o3lossn tendency due to nox additional emission
    REAL(DP), INTENT(OUT) :: o3lossy_air_lg_te ! o3lossy tendency due to nox additional emission
    REAL(DP), INTENT(OUT) :: noxloss_air_lg_te ! noxloss tendency via hno3 formation
    REAL(DP), INTENT(OUT) :: hno3_air_lg_te    ! hno3_air tendency on trajectory
    REAL(DP), INTENT(IN)  :: oh_lg, ho2_lg
    REAL(DP), INTENT(IN)  :: ohprod1_lg_te,ohprod2_lg_te,ohprod3_lg_te
    REAL(DP), INTENT(IN)  :: ohloss1_lg_te,ohloss2_lg_te,ohloss3_lg_te,ohloss4_lg_te,ohloss5_lg_te
    REAL(DP), INTENT(IN)  :: ho2prod1_lg_te,ho2loss1_lg_te,ho2loss2_lg_te
    REAL(DP), INTENT(IN)  :: oh_air_lg,ho2_air_lg
    REAL(DP), INTENT(OUT) :: oh_air_lg_te, ho2_air_lg_te,ch4_air_lg_te
    ! H2O
    REAL(DP), INTENT(IN)    :: h2o_lg             ! background h2o
    REAL(DP), INTENT(INOUT) :: h2o_air_lg         ! new h2o
    REAL(DP), INTENT(IN)    :: h2o_loss_lg        ! h2o loss
    REAL(DP), INTENT(OUT)   :: h2o_air_lg_te      ! h2o air new tendency
    !
    REAL(DP), INTENT(OUT)   :: hno3scav_air_lg_te, hno3loss_air_lg_te
    ! local
    REAL(DP)                :: h2o_loss_air_lg, h2o_air_lg_m1
    REAL(DP)                :: o3_air_lg_m1,o3loss_air_help          &
                             , nox_lg, nox_air_lg_m1, hno3_air_lg_m1
    ! methane variables
    REAL(DP)                :: a0,a1,a2,b0,b1,b2,ohfac,ho2fac,noxfac,o3fac
    ! diagnostic purpose
    ! REAL(DP)                :: PR3,PR4,PR5,DR6,DR7,DR8,DR9,DR10

    REAL(DP), Parameter :: &
           eps_o3     = 1.e-18_dp &   ! minimum value for background o3
         , eps_nox    = 1.e-21_dp &   ! minimum value for background nox
         , eps_o3_air = 1.e-25_dp &
         , eps_h2o    = 1.e-21_dp

    ! NOX, O3, HNO3
    o3_air_lg_m1      = o3_air_lg
    nox_air_lg_m1     = nox_air_lg
    hno3_air_lg_m1    = hno3_air_lg
    h2o_air_lg_m1     = h2o_air_lg
    nox_lg            = no_lg + no2_lg

! NOX
      if (abs(nox_lg)>eps_nox) then
! noxloss = hno3 formation on trajectory
         noxloss_air_lg_te = noxloss_lg_te * nox_air_lg_m1 / nox_lg
      else
         noxloss_air_lg_te = 0._dp
      endif
      hno3_air_lg =  hno3_air_lg_m1 + noxloss_air_lg_te * dt
! proportiate scavenging of hno3 (only negative tendencies)
      if (abs(hno3_lg)>eps_nox) then
         hno3scav_air_lg_te = MIN(0._dp,hno3scav_lg_te) * hno3_air_lg / hno3_lg
      else
         hno3scav_air_lg_te = 0._dp
      endif
      hno3_air_lg = hno3_air_lg + hno3scav_air_lg_te * dt
! reconversion from hno3 to nox via no3 or photolysis of hno3 = hno3loss
      if (abs(hno3_lg)>eps_nox) then
         hno3loss_air_lg_te = hno3loss_lg_te * hno3_air_lg / hno3_lg
      else
         hno3loss_air_lg_te = 0._dp
      endif
! new hno3 air
      hno3_air_lg = hno3_air_lg - hno3loss_air_lg_te * dt
! new hno3 air tendency
      hno3_air_lg_te = (hno3_air_lg - hno3_air_lg_m1) / dt
!
! new nox air
      nox_air_lg = nox_air_lg_m1 + (hno3loss_air_lg_te - noxloss_air_lg_te)* dt
! nox air new tendency
      nox_air_lg_te = (nox_air_lg - nox_air_lg_m1) / dt


! O3 Prod + Loss
      if (abs(nox_lg)>eps_nox) then
         o3prodn_air_lg_te = o3prodn_lg_te * nox_air_lg_m1 / nox_lg
         no2_air_lg    = no2_lg * nox_air_lg_m1 / nox_lg
         if (abs(o3_lg)>eps_o3) then
            o3lossn_air_lg_te =  o3lossn_lg_te * 0.5_dp * (no2_air_lg/no2_lg + o3_air_lg/o3_lg)
         else
            o3lossn_air_lg_te = 0._dp
         endif
      else
        o3prodn_air_lg_te = 0._dp
        o3lossn_air_lg_te = 0._dp
      endif
      if (abs(o3_lg)>eps_o3) then
         o3lossy_air_lg_te = (o3lossy_lg_te / o3_lg) * o3_air_lg
      else
         o3lossy_air_lg_te = 0._dp
      endif

      if (abs(o3_air_lg_m1)>eps_o3_air) then
         o3loss_air_help = (o3lossn_air_lg_te + o3lossy_air_lg_te ) * dt / o3_air_lg_m1
      else
         o3loss_air_help = 0._dp
      endif
      o3_air_lg = (o3_air_lg_m1 + o3prodn_air_lg_te * dt)/  &
                  (1._dp + o3loss_air_help)
! o3 air new tendency
      o3_air_lg_te      = (o3_air_lg - o3_air_lg_m1)  / dt


! H2O
     if (abs(h2o_lg)>eps_h2o) then
        h2o_loss_air_lg = h2o_loss_lg * h2o_air_lg_m1 / h2o_lg
     else
        h2o_loss_air_lg = 0._dp
     endif
     h2o_air_lg = h2o_air_lg_m1 - h2o_loss_air_lg
     h2o_air_lg_te = (h2o_air_lg - h2o_air_lg_m1) / dt


! CH4
! Caution nox_lg=0??
      if ((abs(o3_lg)>eps_o3).and.(abs(nox_lg)>eps_nox)) then
         noxfac = nox_air_lg/nox_lg
         o3fac = o3_air_lg/o3_lg
         a0 = 0.5_dp*(2._dp*ohprod1_lg_te+ohprod2_lg_te-ohloss1_lg_te)*o3fac &
            + 0.5_dp*ohprod3_lg_te*noxfac
         a1 = 0.5_dp*(ohprod2_lg_te+ohprod3_lg_te-ohloss5_lg_te)
         a2 = 0.5_dp*(-ohloss1_lg_te-2._dp*ohloss2_lg_te-2._dp*ohloss3_lg_te-2._dp*ohloss4_lg_te-ohloss5_lg_te)
         b0 = 0.5_dp*(ohloss1_lg_te-ohprod2_lg_te)*o3fac           &
            + 0.5_dp*(2._dp*ho2prod1_lg_te-ohprod3_lg_te)*noxfac
         b1 = 0.5_dp*(-ohprod2_lg_te-ohprod3_lg_te-ohloss5_lg_te-2._dp*ho2loss1_lg_te-ho2loss2_lg_te)
         b2 = 0.5_dp*(ohloss1_lg_te+2._dp*ohloss2_lg_te+ohloss5_lg_te)
         ohfac = (a0*b1-a1*b0)/(a1*b2-a2*b1)
         ho2fac = (a2*b0-a0*b2)/(a1*b2-a2*b1)
         oh_air_lg_te = (oh_lg*ohfac-oh_air_lg)/dt
         ho2_air_lg_te = (ho2_lg*ho2fac-ho2_air_lg)/dt
         ch4_air_lg_te = ohloss4_lg_te*ohfac
      else
         oh_air_lg_te = 0._dp
         ho2_air_lg_te = 0._dp
         ch4_air_lg_te = 0._dp
      endif

!! diagnostics
!      PR3 = ohprod1_lg_te*o3fac
!      PR4 = 0.5_dp*ohprod2_lg_te*(ho2fac+o3fac)
!      PR5 = 0.5_dp*ohprod3_lg_te*(ho2fac+noxfac)
!      DR6 = 0.5_dp*ohloss1_lg_te*(ohfac+o3fac)
!      DR7 = ohloss2_lg_te*ohfac
!      DR8 = ohloss3_lg_te*ohfac
!      DR9 = ohloss4_lg_te*ohfac
!      DR10 =0.5_dp*ohloss5_lg_te*(ohfac+ho2fac)

  END SUBROUTINE airtrac_integrate
! ----------------------------------------------------------------------

END MODULE messy_airtrac
! **********************************************************************
