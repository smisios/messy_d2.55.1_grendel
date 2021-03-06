! **********************************************************************
!
! SUBMODEL CORE LAYER (SMCL) routines for MESSy SUBMODEL CONTRAIL
!
! This submodel is used to provide information for 
!      - contrail formation 
!      - potential contrail coverage
! 
! Author : Christine Froemming,Volker Grewe,Sabine Brinkop, DLR-IPA, 2011/2012
!
! References: (see below)
!
! **********************************************************************

! **********************************************************************
MODULE messy_contrail
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'contrail'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.1'

  ! CTRL-NAMELIST PARAMETERS
  ! Emissionindex of water vapour [kg(H2O)/kg(fuel)]
  REAL(DP), PUBLIC :: ei_h2o = 1.25_dp
  ! Combustion heat of fuel [J/kg] (Schumann et al., 2000)
  REAL(DP), PUBLIC :: q_fuel = 43.2e6_dp
  ! overall propulsion efficiency of aircraft [fraction] (Schumann et a., 2000)
  REAL(DP), PUBLIC :: eta_ac = 0.31_dp
  !qqq at the moment, r_sac NOT used
  REAL(DP), PUBLIC :: a_sac = 0.9_dp   ! (Burkhardt et al., 2008)
  REAL(DP), PUBLIC :: r_sac = 1.1_dp   ! (Burkhardt et al., 2008)

  ! PUBLIC SUBROUTINES
  PUBLIC :: contrail_read_nml_ctrl
  PUBLIC :: contrail_pot_cov
  PUBLIC :: contrail_calc
  PUBLIC :: contrail_calc_dev
  PUBLIC :: contrail_uv_grad

  ! PRIVATE SUBROUTINES
  !PRIVATE :: contrail_spread
  !PRIVATE :: contrail_sedi

CONTAINS

  ! #########################################################################
  ! PUBLIC SUBROUTINES
  ! #########################################################################

  ! =========================================================================
  SUBROUTINE contrail_read_nml_ctrl(status, iou)

    ! ------------------------------------------------------------------
    ! This routine is used to read the CTRL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy INTERFACE
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ EI_H2O, Q_fuel, eta_ac, a_SAC, r_SAC

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='contrail_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! ### ADD HERE DIAGNOSTIC OUPUT FOR LOG-FILE
    WRITE(*,*) 'Emission index for H2O: ',EI_H2O,' kg/kg(fuel)'
    WRITE(*,*) 'Combustion heat:        ',Q_fuel,' J/kg(fuel)'
    WRITE(*,*) 'Overall propulsion efficiency ',eta_ac, ' (fraction)'
    WRITE(*,*) 'Ice supersaturation parameter ',a_SAC,' (-)'       
    WRITE(*,*) 'Rel. humidity for Schmitt-Appleman-Criterium ',r_SAC, ' (-)'


    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE contrail_read_nml_ctrl
  ! =========================================================================

  ! =========================================================================
  ELEMENTAL SUBROUTINE contrail_pot_cov(t,p,q,zqte,b_ci,zrhc,dt, & !IN
       b_cc,             & !OUT
       potcov,           & !OUT
       zqsm1, zconpn)      !OUT

    USE messy_main_tools,           ONLY: tlucua,jptlucu1,jptlucu2
    USE messy_main_constants_mem,   ONLY: tmelt, vtmpc1, g, cp_air

    implicit none
    INTRINSIC :: MAX, EXP, LOG, MERGE, MIN, NINT, SQRT

    ! Atmospheric input
    REAL(DP), INTENT(IN) :: T     ! air temperature [K]
    REAL(DP), INTENT(IN) :: p     ! pressure [Pa]
    REAL(DP), INTENT(IN) :: q     ! water vapour [kg/kg]
    REAL(DP), INTENT(IN) :: b_ci  ! cloud cover [fraction]
    REAL(DP), INTENT(IN) :: zrhc  ! critical humidity for nat. clouds
    REAL(DP), INTENT(IN) :: zqte  ! water vapour tendency [kg/kg/s]
    REAL(DP), INTENT(IN) :: dt      ! time step length
    !
    REAL(DP), INTENT(OUT) :: b_cc   ! potential contrail cirrus coverage = ISS (contrail+cirrus)
    REAL(DP), INTENT(OUT) :: potcov ! potential contrail coverage
    REAL(DP), INTENT(OUT) :: zqsm1  !
    REAL(DP), INTENT(OUT) :: zconpn !

    ! UNKNOWN PARAMETERS (TAKEN FROM ECHAM4 CODE)
    REAL(DP), parameter  :: C1ES = 610.78_dp,     & ! taken from E4
         C3LES = 17.269_dp,    &
         C3IES = 21.875_dp,    & 
         C4LES = 35.86_dp,     &
         C4IES = 7.66_dp,      & 
         zcaa = 0.0059_dp,     &
         zcab = 0.003102_dp,   &
         ZCP=1004._dp,         & ! heat capacity (of what?) [J/kg/K] 
         ALV=2.5008E6_dp,      &
         ALS=2.8345E6_dp,      &
         ZEPSI=vtmpc1+1._dp

    ! LOCAL VARIABLES
    REAL(DP) :: ZG  ! slope of the mixing line [Pa/K]

    REAL(DP) :: q_env ! water vapour mixing ratios 
    ! env=environment = no cloud area
    REAL(DP) :: ZQR
    REAL(DP) :: r_nuc  ! nucleation threshold for rel. humidity, and saturation
    REAL(DP) :: l_over_ice
    REAL(DP) :: b_totn, b_tot
    REAL(DP) :: ztc, zalpha,zzfac,zsat,zrwkr,zqswp1
    REAL(DP) :: zsubt,zqhelp,zcor
    REAL(DP) :: zqsatkr,zeswp1,zgfac
    REAL(DP) :: zqconpa,zqcon
    REAL(DP) :: zzqt,zcons1,zrcp,zlsdcp,zlvdcp 
    REAL(DP) :: zcvm4,zcvm5,c5les,c5ies,alf
    INTEGER  :: it
    LOGICAL  :: lo1,lo2
    REAL(DP) :: r_cc   !
    REAL(DP) :: q_cc   !
    REAL(DP) :: T_cont ! critical temperature[K] for contrail formation
    REAL(DP) :: r_cont ! critica rel humidity [?] for contrail formation
    REAL(DP) :: q_cont ! spec. humidity [kg/kg]  for contrail formation

    ! saturation specific humidity
    it=NINT(t*1000._dp)
    it = MAX(MIN(it,jptlucu2),jptlucu1)
    zqsm1=tlucua(it)/p
    zqsm1=MIN(zqsm1,0.5_dp)
    zcor=1._dp/(1._dp-vtmpc1*zqsm1)
    zqsm1=zqsm1*zcor

    ZESWP1=C1ES*EXP(C3LES*(T-TMELT)/(T-C4LES))

    ! critical temperature and humidity for contrail formation 
    ! (Schmidt-Appleman)
    ZGFAC=EI_H2O*ZCP*ZEPSI/(Q_fuel*(1._dp-eta_ac))
    ZG=ZGFAC*P

    ! INIT (ELSE BRANCH MISSING BELOW!)
    ZQSWP1 = 0.0_dp
    ZTC = 0.0_dp
    b_cc = b_ci
    b_totn = b_ci
    Q_CONT = 0.0_dp
    r_cc = 0.0_dp
    t_cont = 0.0_dp

    IF (ZG.gt.0.053_dp) THEN
       ZSUBT=LOG(ZG-0.053_dp)
       T_CONT=226.7_dp+9.43_dp*ZSUBT+0.72_dp*ZSUBT*ZSUBT
       ZQSATKR=C1ES*EXP(C3LES*(T_CONT-TMELT)/(T_CONT-C4LES))
       IF (ZESWP1.GT.0._dp) THEN 
          ZRWKR=(ZG*(T-T_CONT)+ZQSATKR)/ZESWP1
          ZRWKR=MAX(ZRWKR,0._dp)

          ZESWP1=ZESWP1/(ZEPSI*P)
          ZCOR=1._dp/(1._dp-VTMPC1*ZESWP1)
          ZQSWP1=ZESWP1*ZCOR

          ! potential contrail cirrus cover

          ! homogeneous nucelation threshold
          ZSAT=1._dp
          R_NUC=2.349_dp-T/259._dp
          L_OVER_ICE=EXP(-(210368._dp+131.438_dp*T-3.32373E+06_dp&
               /T-41729.1_dp*LOG(T))/(8.31441_dp*T))
          R_CC=ZRHC/(R_NUC*a_sac)
          Q_CC=zqsm1*ZSAT*R_CC
          ZQR=zqsm1*ZSAT*ZRHC

          IF(B_CI.LT.1._dp) THEN
             Q_ENV=(Q-B_CI* zqsm1)/(1._dp-B_CI)
          ELSEIF (B_CI.EQ.1._dp) THEN
             Q_ENV=zqsm1
          ENDIF

          ! potential contrail cirrus cover
          B_TOTN=(Q_ENV-Q_CC)/((zqsm1-(ZQR-Q_CC))*ZSAT-Q_CC)
          B_TOTN=MAX(B_TOTN,0._dp)
          B_TOTN=MIN(B_TOTN,1._dp)
          B_TOTN=1._dp-SQRT(1._dp-B_TOTN)
          B_CC=B_TOTN-B_CI
          B_CC=MAX(B_CC,0._dp)

          ! probability for ice liquidwater at temperature t
          ZTC=T-TMELT
          LO1 = ZTC.LT.0._dp
          ZZFAC=MERGE(1._dp,0._dp,LO1)
          ZALPHA=(1._dp-ZZFAC)+ZZFAC*(ZCAA+(1._dp-ZCAA)*EXP(-ZCAB*ZTC**2._dp))
          B_CC=(1._dp-ZALPHA)*B_CC

          IF(T.LT.T_CONT) THEN

             ! combination of rcrit for cirrus and rcrit appleman
             R_CONT=ZRWKR*ZRHC/(R_NUC*a_sac*L_OVER_ICE)
             Q_CONT=zqsm1*ZSAT*R_CONT

             IF(B_CI.LT.1._dp) THEN
                Q_ENV=(Q-B_CI* zqsm1)/(1._dp-B_CI)
             ELSEIF (B_CI.EQ.1._dp) THEN
                Q_ENV=zqsm1
             ENDIF

             ! potential contrail coverage
             B_TOT=(Q_ENV-Q_CONT)/((zqsm1-(ZQR-Q_CONT))*ZSAT-Q_CONT)
             B_TOT=MAX(B_TOT,0._dp)
             B_TOT=MIN(B_TOT,1._dp)
             B_TOT=1._dp-SQRT(1._dp-B_TOT)
             POTCOV=B_TOT-B_CI
             POTCOV=MAX(POTCOV,0._dp)
             POTCOV=(1._dp-ZALPHA)*POTCOV
          ELSE
             POTCOV=0._dp
             B_TOT=B_CI
          ENDIF
       ELSE
          POTCOV=0._dp
          B_TOT=B_CI
       ENDIF
    ELSE 
       POTCOV=0._dp
       B_TOT=B_CI
    ENDIF
    POTCOV=MIN(POTCOV,B_CC)

    ! no contrails below humidity threshold. only persistent contails
    IF (ZQSWP1.GT.0._dp) THEN
       ZQHELP=ZRHC*zqsm1*zqsm1/ZQSWP1
       IF(Q.LT.ZQHELP) THEN
          POTCOV=0._dp
          B_TOT=B_CI
       ENDIF
    ELSE
       POTCOV=0._dp
       B_TOT=B_CI
    ENDIF

    ! condensation rate
    ALF=ALS-ALV         
    ZCONS1=cp_air/(ALF*G*dt)
    ZRCP=1./(cp_air+ZCONS1*Q)
    ZLVDCP=ALV*ZRCP
    ZLSDCP=ALS*ZRCP
    C5LES=C3LES*(TMELT-C4LES)
    C5IES=C3IES*(TMELT-C4IES)
    LO2 = ZTC.GT.0._dp
    ZCVM4=MERGE(C4LES,C4IES,LO2)
    ZCVM5=MERGE(C5LES*ZLVDCP,C5IES*ZLSDCP,LO2)
    ZZQT=1._dp/(T-ZCVM4)
    ZQCON=1._dp/(1._dp+ZCVM5*zqsm1*ZCOR*ZZQT**2._dp)

    ! ???
    ZQCONPA=MAX(ZQSM1,Q)
    ZCONPN=ZQCONPA+zqte*dt-zqsm1
    ZCONPN=ZCONPN*ZQCON

  END SUBROUTINE contrail_pot_cov
  ! =========================================================================

  ! =========================================================================
  ELEMENTAL SUBROUTINE contrail_calc(    &
       q, b_ci, zconrat, potcov, zconpn, & !IN
       emis, scal, gboxarea,             & ! IN
       pconcov, pconiwc ) ! OUT

    IMPLICIT NONE
    INTRINSIC :: EPSILON, MAX, MIN 

    REAL(DP), INTENT(IN)  :: q
    REAL(DP), INTENT(IN)  :: b_ci   
    REAL(DP), INTENT(IN)  :: zconrat  ! condensation rate from cloud 
    REAL(DP), INTENT(IN)  :: potcov   ! potential contrail coverage
    REAL(DP), INTENT(IN)  :: zconpn   !
    REAL(DP), INTENT(IN)  :: emis     !
    REAL(DP), INTENT(IN)  :: scal     !
    REAL(DP), INTENT(IN)  :: gboxarea !
    REAL(DP), INTENT(OUT) :: pconcov  ! contrail coverage
    REAL(DP), INTENT(OUT) :: pconiwc  ! contrail ice water mixing ratio (kg/kg)

    ! LOCAL VARIABLES
    REAL(DP)              :: zcontra

! Correct scalingfactor for gp contrail cover to be determined later
!    pconcov = potcov * emis * scal

!    if ((emis>EPSILON(1.0_dp)).and.(gboxarea>EPSILON(1.0_dp))) then
    if ((emis>0._dp).and.(gboxarea>EPSILON(1.0_dp))) then
       pconcov=(potcov*SQRT(gboxarea)*200.)/gboxarea
    else
       pconcov=0._dp
    endif

    ! contrail iwc 
    IF(pconcov.GT.0._dp) THEN
       ZCONTRA=PCONCOV*ZCONPN
    ELSE
       ZCONTRA=0._dp
    ENDIF

    ! avoid negative condensation rates and transparent contrails 
    IF(ZCONTRA.LE.0._dp) THEN
       ZCONTRA=0._dp
       PCONCOV=0._dp
    ENDIF
    ! Condensation in contrails must not exceed the residue water vapour after cloud formation
    IF(ZCONTRA.GT.0._dp) THEN
       ZCONTRA=MIN(ZCONTRA,Q-ZCONRAT)
    ENDIF
    ! if no available water for condensation in contrails avoid contrail cover
    IF(ZCONTRA.LE.0._dp) THEN
       ZCONTRA=0._dp
       PCONCOV=0._dp
       pconiwc=0._dp
    ELSE
       pconiwc=ZCONTRA
    ENDIF

    ! concov must not be greater than potcov and coniwc not greater than residue water vapour

    pconcov = MIN(pconcov,potcov)
    pconiwc = MIN(pconiwc,Q-ZCONRAT)

    ! correction of IWC in contrail part
    IF (PCONCOV>EPSILON(1.0_dp)) THEN
       PCONIWC = PCONIWC/PCONCOV
    else
       PCONIWC = 0._dp
    ENDIF
    ! total cover (cirrus+contrails) must not exceed 100%
    PCONCOV=MIN(PCONCOV,(1.-B_CI))
    PCONCOV=MAX(PCONCOV,0._dp)

    ! update of IWC in contrail part
    PCONIWC = PCONIWC*PCONCOV
    PCONIWC = MAX(PCONIWC,0._dp)

  END SUBROUTINE Contrail_calc
  ! =========================================================================

  ! =========================================================================
  ELEMENTAL SUBROUTINE contrail_calc_dev(&
       q, b_ci, zconrat, potcov, gboxarea, zconpn,                & ! IN
       dt, emis, scal,                                            & ! IN
       potcov_m1, pconcov_m1, pconiwc_m1,                         & ! IN
       rho , eta_dot, dudz, dvdz,                                 & ! IN
       pconcov_sum, pconiwc_sum, pconcov_now, pconiwc_now,        & ! OUT
       concov_te_spread,coniwc_te_sedi,coniwc_te_pot )              ! OUT

    IMPLICIT NONE
    INTRINSIC :: EPSILON, MAX, MIN

    REAL(DP), INTENT(IN)  :: q,b_ci   
    REAL(DP), INTENT(IN)  :: zconrat  ! condensation rate from cloud 
    REAL(DP), INTENT(IN)  :: potcov   ! potential contrail coverage
    REAL(DP), INTENT(IN)  :: gboxarea
    REAL(DP), INTENT(IN)  :: zconpn   !
    REAL(DP), INTENT(IN)  :: dt
    REAL(DP), INTENT(IN)  :: emis     !
    REAL(DP), INTENT(IN)  :: scal     !
    ! FOR T-1:
    REAL(DP), INTENT(IN)  :: potcov_m1  ! potential contrail coverage (t-1)
    REAL(DP), INTENT(IN)  :: pconcov_m1 ! contrail coverage (t-1)
    REAL(DP), INTENT(IN)  :: pconiwc_m1 ! iwc (t-1)
    REAL(DP), INTENT(IN)  :: rho
    REAL(DP), INTENT(IN)  :: eta_dot
    REAL(DP), INTENT(IN)  :: dudz
    REAL(DP), INTENT(IN)  :: dvdz
    !
    ! contrail coverage:
    REAL(DP), INTENT(OUT) :: pconcov_sum  ! "sum" of aged plus act. contribution 
    REAL(DP), INTENT(OUT) :: pconcov_now  ! contribution from act. source only
    ! contrail ice water mixing ratio (kg/kg):
    REAL(DP), INTENT(OUT) :: pconiwc_sum  ! "sum" of aged plus act. contribution 
    REAL(DP), INTENT(OUT) :: pconiwc_now  ! contribution from act. source only
    ! diagnostic tendencies:
    REAL(DP), INTENT(OUT) :: concov_te_spread  ! contrail coverage
    REAL(DP), INTENT(OUT) :: coniwc_te_sedi    ! ice water mixing ratio
    REAL(DP), INTENT(OUT) :: coniwc_te_pot     ! ice water mixing ratio

    ! LCCAL VARIABLES
    REAL(DP) :: pconcov_age ! spread,sedi ... coverage ...
    REAL(DP) :: pconiwc_age !  ... from t-1

    ! calculate contrail coverage from emission in actual time step
    CALL contrail_calc( &
         q, b_ci, zconrat, potcov, zconpn, & !IN
         emis, scal, gboxarea, & ! IN
         pconcov_now, pconiwc_now )

    ! initialize pconcov_old and pconiwc_old
    pconcov_age = 0._dp
    pconiwc_age = 0._dp

    ! development of contrail coverage and iwc from last timestep  
    !
    ! change of potential contrail coverage
    if (potcov_m1>EPSILON(1.0_dp)) then
       ! mit coniwc_m1 gewichten
       coniwc_te_pot = (((potcov-potcov_m1)/potcov_m1)*pconiwc_m1)/dt
    else
       coniwc_te_pot = 0._dp
    endif

    ! spreading calculated according to vertical wind gradient

    if (potcov_m1>EPSILON(1.0_dp)) then
       call contrail_spread(concov_te_spread, dudz, dvdz, pconcov_m1, gboxarea)
    else
       concov_te_spread = 0._dp
    endif

    ! sedimentation
    if (potcov_m1>EPSILON(1.0_dp)) then
       ! Trajectory only
       call contrail_sedi (pconiwc_m1, rho, coniwc_te_sedi, eta_dot)
    else
       coniwc_te_sedi = 0._dp
    endif

    ! add tendencies, but pconcov_old not larger than actual potcov
    if (pconcov_m1.gt.0._dp) then
       pconcov_age = pconcov_m1 + (concov_te_spread * dt)
       pconcov_age = MIN(pconcov_age,potcov)
    else
       pconcov_age =0._dp
    endif

    if (pconiwc_m1.gt.0._dp) then
       pconiwc_age = pconiwc_m1 + (coniwc_te_sedi + coniwc_te_pot)* dt
    else
       pconiwc_age = 0._dp
    endif

    pconcov_sum = pconcov_age + pconcov_now
    pconcov_sum = MIN(pconcov_sum,potcov)
    pconiwc_sum = pconiwc_age + pconiwc_now
    pconiwc_sum = MIN(pconiwc_sum,Q-ZCONRAT)

    ! correction of IWC in contrail part
    IF (PCONCOV_sum>EPSILON(1.0_dp)) THEN
       PCONIWC_sum = PCONIWC_sum/PCONCOV_sum
    else
       PCONIWC_sum = 0._dp
    ENDIF
    ! total cover (cirrus+contrails) must not exceed 100%
    PCONCOV_sum=MIN(PCONCOV_sum,(1.-B_CI))
    PCONCOV_sum=MAX(PCONCOV_sum,0._dp)

    ! update of IWC in contrail part
    PCONIWC_sum = PCONIWC_sum*PCONCOV_sum
    PCONIWC_sum = MAX(PCONIWC_sum,0._dp)

  END SUBROUTINE contrail_calc_dev
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE contrail_uv_grad(u_scb,v_scb,sqcst_2d,rho,zpf,dudz,dvdz)

    USE messy_main_constants_mem, ONLY: g

    IMPLICIT NONE
    INTRINSIC :: SIZE

    ! I/O
    REAL(DP), INTENT(IN)   :: u_scb(:,:)
    REAL(DP), INTENT(IN)   :: v_scb(:,:)
    REAL(DP), INTENT(IN)   :: sqcst_2d(:)
    REAL(DP), INTENT(IN)   :: rho(:,:)
    REAL(DP), INTENT(IN)   :: zpf(:,:)
    !
    REAL(DP), INTENT(INOUT):: dudz(:,:)
    REAL(DP), INTENT(INOUT):: dvdz(:,:)

    ! LOCAL
    REAL(DP), DIMENSION(:,:), POINTER :: udcos
    REAL(DP), DIMENSION(:,:), POINTER :: vdcos
    REAL(DP), DIMENSION(:),   POINTER :: zdz
    REAL(DP), DIMENSION(:),   POINTER :: zsqcst

    INTEGER  :: kproma, nlev, jk, jl, jkp1, jkm1

    kproma = SIZE(u_scb,1)
    nlev   = SIZE(u_scb,2)

    NULLIFY(udcos)
    NULLIFY(vdcos)
    NULLIFY(zdz)
    NULLIFY(zsqcst)

    allocate(udcos(kproma,nlev))
    allocate(vdcos(kproma,nlev))
    allocate(zdz  (kproma))
    allocate(zsqcst(kproma))

    zsqcst(:) = 1._dp/sqcst_2d(:)

    do jk = 1, nlev
       do jl = 1, kproma
          ! wind divided by cos(lat)
          udcos(jl,jk) = u_scb(jl,jk) * zsqcst(jl)
          vdcos(jl,jk) = v_scb(jl,jk) * zsqcst(jl)
       enddo
    enddo

    DO JK = 1, NLEV
       do jl=1,kproma
          IF (JK.LT.NLEV) THEN
             JKP1 = JK+1
          ELSE
             JKP1 = JK
          END IF

          IF (JK.GT.1) THEN
             JKM1 = JK-1
          ELSE
             JKM1 = JK
          END IF

          ZDZ(jl) = (ZPF(jl,JKP1) - ZPF(jl,JKM1))/(RHO(jl,JK)*G)

          DUDZ(jl,jk) = (udcos(jl,JKM1)  -  (udcos(jl,JKP1))) / ZDZ(jl)

          DVDZ(jl,jk) = (-vdcos(jl,JKM1)  - (-vdcos(jl,JKP1))) / ZDZ(jl)
       end do
    end do

    ! release space
    DEALLOCATE(udcos)  ; NULLIFY(udcos)
    DEALLOCATE(vdcos)  ; NULLIFY(vdcos)
    DEALLOCATE(zdz)    ; NULLIFY(zdz)
    DEALLOCATE(zsqcst) ; NULLIFY(zsqcst)

  END SUBROUTINE contrail_uv_grad
  ! =========================================================================

  ! #########################################################################
  ! PRIVATE SUBROUTINES
  ! #########################################################################

  ! =========================================================================
  ! op_sb_20120411+
  ELEMENTAL SUBROUTINE contrail_spread(concov_te_spread &
       , dudz,dvdz, zconcov_m1,gboxarea)

    ! Parametrisation of the spread of contrails according to
    ! Burkhardt und K?rcher (2009)


    IMPLICIT NONE
    INTRINSIC :: SQRT

    REAL(DP), INTENT(IN) :: gboxarea
    REAL(DP), INTENT(IN) :: zconcov_m1
    REAL(DP), INTENT(IN) :: dudz
    REAL(DP), INTENT(IN) :: dvdz

    REAL(DP), INTENT(OUT):: concov_te_spread ! tendency contrail cover

    REAL(DP), PARAMETER  :: H = 200._dp  ! thickness of contrail 200m (1h)
    REAL(DP), PARAMETER  :: c = 0.72_dp  ! depends on angle of the vertical wind shear 
    ! vector with the flight direction 
    ! hier muss contrail cover eingehen

    concov_te_spread = c * sqrt(dudz**2+dvdz**2 ) &
         * (H *zconcov_m1/ sqrt(gboxarea))

  END SUBROUTINE contrail_spread
  ! =========================================================================
  ! op_sb_20120411-

  ! op_sb_20120404+
  ! =========================================================================
  ELEMENTAL SUBROUTINE contrail_sedi(lwi_lg, rho_lg, lwi_lg_tte, w_wind)

    ! Parametrisation of sedimentation according to
    ! Heymsfield and Donner (1990)
    !
    !  dri/dt = 1/rho*dF/dz   
    !  v      = alpha *(rho * ri)**beta
    !  F      = v*ri  F_top = 0 (no ice from above)

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: lwi_lg     ! ice water content (kg/kg)
    REAL(DP), INTENT(IN) :: rho_lg     ! air density (kg/m3)
    REAL(DP), INTENT(IN) :: w_wind     ! vertical wind (1/s, etadot)
    REAL(DP), INTENT(OUT):: lwi_lg_tte ! tendency of lwi due to sedimentation

    ! LOCAL VARIABELS
    REAL(dp), PARAMETER :: alpha  = 3.29_dp
    REAL(dp), PARAMETER :: beta   = 0.16_dp
    REAL(dp), PARAMETER :: height = 200._dp    !500._dp  ! height of contrail (m) test: 1000m

    ! sedimentation only if contrail moves upward
    if (w_wind .lt. 0._dp) then              ! upward
       !    lwi_lg_tte = (- alpha) * (rho_lg * lwi_lg)**(beta) * lwi_lg / height
       lwi_lg_tte = -(alpha * (rho_lg * lwi_lg)**(beta) * lwi_lg / height)
    else                                     ! downward or zero
       lwi_lg_tte = 0._dp
    endif

  END SUBROUTINE contrail_sedi
  ! =========================================================================

! **********************************************************************
END MODULE messy_contrail
! **********************************************************************

