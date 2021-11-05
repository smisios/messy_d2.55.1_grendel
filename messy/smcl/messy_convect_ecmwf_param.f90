MODULE MESSY_CONVECT_ECMWF_PARAM

!      Collection of convection specific functions, constants and parameters
!      as well as their declaration 


!    Author of this module: 
!    little adjustements due to implementation into the MESSy system
!    H. Tost,    MPI - Chemie, Mainz Feb. 2005

!    Original code as described in the routines

  USE MESSY_MAIN_CONSTANTS_MEM,     ONLY: dp, RG => g, RPI => PI,           &
                                          R => R_gas, M_air, M_H2O, rho_h2o,&
                                          RKAP => c_vKar
  IMPLICIT NONE
  SAVE
  PRIVATE

!!!   PUBLIC SUBROUTINES
  PUBLIC :: SATUR, VDIV, VEXP, init_convection_constants

!!!   PUBLIC FUNCTIONS
  PUBLIC :: FOEALFA, FOEALFCU, FOEDEM, FOELDCPM, FOELDCPMCU, FOEEWM
  PUBLIC :: FOEEWMCU, FOEDEMCU, FOELHMCU

!!!  PARAMETERS
  REAL(dp), PARAMETER :: RLVTT = 2.5008E+6_dp
  REAL(dp), PARAMETER :: RLSTT = 2.8345E+6_dp
  REAL(dp), PARAMETER :: RLMLT = RLSTT - RLVTT
  REAL(dp), PARAMETER :: RTT   = 273.16_dp
  REAL(dp), PARAMETER :: RD    = 1000._dp * R / M_AIR
  REAL(dp), PARAMETER :: RV    = 1000._dp * R / M_H2O
  REAL(dp), PARAMETER :: RETV  = RV / RD - 1._dp
  
  REAL(dp), PARAMETER :: RCPD  = 3.5_dp * RD
  REAL(dp), PARAMETER :: RCPV  = 4.0_dp * RV

  PUBLIC :: RG, RPI, RCPD, RLVTT, RLSTT, RLMLT, RTT, RETV, RD, RKAP
  PUBLIC :: RVTMP2, R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES,     &
            R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTBER,  &           
            RTBERCU, RTICE, RTICECU, RTWAT_RTICE_R, RTWAT_RTICECU_R 
  PUBLIC :: ENTRPEN, ENTRSCV, ENTRMID, ENTRDD, RMFCTOP, RMFCMAX,      &
            RMFCMIN, RMFDEPS, RDEPTHS, RHCDD, RPRCON, RCPECONS,       &
            RCUCOV, RTAUMEL, RHEBC, RTAU, RMFCFL, RMFSHCFL, RMFSOLUV, &
            RMFSOLCT, NJKT1, NJKT2, NJKT3, LMFSCL_WSTAR, LMFTRAC
  PUBLIC :: LMFDUDV, LMFMID, LMFDD, LEPCLD, LMFPEN, LMFSCV, NSMAX
  PUBLIC :: LPHYLIN, RLPTRC, RLPAL1, RLPAL2
  PUBLIC :: N_VMASS, LHOOK, DR_HOOK
  PUBLIC :: RLMIN

!========================================================================================

!     ------------------------------------------------------------------
!*     *YOETHF* DERIVED CONSTANTS SPECIFIC TO ECMWF THERMODYNAMICS
!     ------------------------------------------------------------------

  REAL(dp), PARAMETER :: RVTMP2          = RCPV / RCPD - 1._dp
  REAL(dp), PARAMETER :: R2ES            = 611.21_dp * RD / RV
  REAL(dp), PARAMETER :: R3LES           = 17.502_dp
  REAL(dp), PARAMETER :: R3IES           = 22.587_dp
  REAL(dp), PARAMETER :: R4LES           = 32.19_dp
  REAL(dp), PARAMETER :: R4IES           = - 0.7_dp
  REAL(dp), PARAMETER :: R5LES           = R3LES *(RTT - R4LES)
  REAL(dp), PARAMETER :: R5IES           = R3IES *(RTT - R4IES)
  REAL(dp), PARAMETER :: R5ALVCP         = R5LES * RLVTT / RCPD
  REAL(dp), PARAMETER :: R5ALSCP         = R5IES * RLSTT / RCPD
  REAL(dp), PARAMETER :: RALVDCP         = RLVTT / RCPD
  REAL(dp), PARAMETER :: RALSDCP         = RLSTT / RCPD
  REAL(dp), PARAMETER :: RALFDCP         = RLMLT / RCPD
  REAL(dp), PARAMETER :: RTWAT           = RTT
  REAL(dp), PARAMETER :: RTBER           = RTT - 5.0_dp
  REAL(dp), PARAMETER :: RTBERCU         = RTT - 5.0_dp
  REAL(dp), PARAMETER :: RTICE           = RTT - 23._dp
  REAL(dp), PARAMETER :: RTICECU         = RTT - 23._dp
  REAL(dp), PARAMETER :: RTWAT_RTICE_R   = 1._dp /(RTWAT - RTICE)
  REAL(dp), PARAMETER :: RTWAT_RTICECU_R = 1._dp /(RTWAT - RTICECU) 
 

!     J.-J. MORCRETTE                   91/07/14  ADAPTED TO I.F.S.
!
!      NAME     TYPE      PURPOSE
!      ----     ----      -------
!
!     *R__ES*   REAL      *CONSTANTS USED FOR COMPUTATION OF SATURATION
!                         MIXING RATIO OVER LIQUID WATER(*R_LES*) OR
!                         ICE(*R_IES*).
!     *RVTMP2*  REAL      *RVTMP2=RCPV/RCPD-1.
!     *RHOH2O*  REAL      *DENSITY OF LIQUID WATER.   (RATM/100.)
!     *R5ALVCP* REAL      *R5LES*RLVTT/RCPD
!     *R5ALSCP* REAL      *R5IES*RLSTT/RCPD
!     *RALVDCP* REAL      *RLVTT/RCPD
!     *RALSDCP* REAL      *RLSTT/RCPD
!     *RALFDCP* REAL      *RLMLT/RCPD
!     *RTWAT*   REAL      *RTWAT=RTT
!     *RTBER*   REAL      *RTBER=RTT-0.05
!     *RTBERCU  REAL      *RTBERCU=RTT-5.0
!     *RTICE*   REAL      *RTICE=RTT-0.1
!     *RTICECU* REAL      *RTICECU=RTT-23.0
!     *RTWAT_RTICE_R*   REAL      *RTWAT_RTICE_R=1./(RTWAT-RTICE)
!     *RTWAT_RTICECU_R* REAL      *RTWAT_RTICECU_R=1./(RTWAT-RTICECU)
!
!     ------------------------------------------------------------------
!*     END OF *YOETHF* DERIVED CONSTANTS SPECIFIC TO ECMWF THERMODYNAMICS
!     ------------------------------------------------------------------
!=======================================================================================
!     ----------------------------------------------------------------
!*    ** *YOECUMF* - PARAMETERS FOR CUMULUS MASSFLUX SCHEME
!     ----------------------------------------------------------------
!     1.           SPECIFY PARAMETERS FOR MASSFLUX-SCHEME
!                  --------------------------------------

!     ENTRPEN: AVERAGE ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
!     -------

!  REAL(dp), PARAMETER :: ENTRPEN    = 1.0E-4_dp
  REAL(dp), PARAMETER :: ENTRPEN    = 1.2E-4_dp
!  REAL(dp), PARAMETER :: ENTRPEN    = 5.2E-5_dp

!     ENTRSCV: AVERAGE ENTRAINMENT RATE FOR SHALLOW CONVECTION
!     -------

  REAL(dp), PARAMETER :: ENTRSCV    = 3.0E-4_dp
!  REAL(dp), PARAMETER :: ENTRSCV    = 1.0E-4_dp

!     ENTRMID: AVERAGE ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
!     -------

  REAL(dp), PARAMETER :: ENTRMID    = 1.0E-4_dp
!  REAL(dp), PARAMETER :: ENTRMID    = 5.0E-5_dp

!     ENTRDD: AVERAGE ENTRAINMENT RATE FOR DOWNDRAFTS
!     ------

  REAL(dp), PARAMETER :: ENTRDD     = 2.0E-4_dp

!     RMFCTOP:   RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVE
!     -------

  REAL(dp), PARAMETER :: RMFCTOP    = 0.33_dp

!     RMFCMAX:   MAXIMUM MASSFLUX VALUE ALLOWED FOR UPDRAFTS ETC
!     -------

  REAL(dp), PARAMETER :: RMFCMAX    = 1.0_dp

!     RMFCMIN:   MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     -------

  REAL(dp), PARAMETER :: RMFCMIN    = 1.E-10_dp

!     RMFDEPS:   FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     -------

  REAL(dp), PARAMETER :: RMFDEPS    = 0.3_dp

!     RDEPTHS:   MAXIMUM ALLOWED SHALLOW CLOUD DEPTH (Pa)
!     -------

  REAL(dp), PARAMETER :: RDEPTHS    = 2.E4_dp

!     RPRCON:    COEFFICIENTS FOR DETERMINING CONVERSION FROM CLOUD WATE
!     ------

! mz_ht_20070911+
  !original
  !REAL(dp), PARAMETER :: RPRCON     = 1.5E-3_dp !increased value to avoid UTropo warm bias
  ! modified for T21L19, gives global OLR close to ERBE
  !  REAL(dp), PARAMETER :: RPRCON     = 5.0E-3_dp 
  ! modified for T63L87, gives global OLR close to ERBE free
  !  REAL(dp), PARAMETER :: RPRCON     = 9.0E-4_dp 
  ! modified for T63L87, gives global OLR close to ERBE nudged
    REAL(dp), PARAMETER :: RPRCON     = 1.5E-3_dp 
! mz_ht_20070911-

!                COEFFICIENTS FOR RAIN EVAPORATION BELOW CLOUD
!                AND MELTING
!                ---------------------------------------------
!     RCPECONS:  KESSLER COEFFICIENT
!     RCUCOV:    ASSUMED CONVECTIVE CLOUD COVER
!     RTAUMEL:   MELTING TIME SCALE
!     RHEBC:     CRITICAL RELATIVE HUMIDITY BELOW CLOUD  FOR EVAPORATION

  REAL(dp), PARAMETER :: RCUCOV      = 0.05_dp
  REAL(dp), PARAMETER :: RCPECONS    = 5.44E-4_dp / RG
  REAL(dp), PARAMETER :: RTAUMEL     = 5._dp * 3.6E3_dp
  REAL(dp), PARAMETER :: RHEBC       = 0.8_dp

!     NEXT VALUE IS RELATIVE SATURATION IN DOWNDRAFRS
!     BUT IS NO LONGER USED ( FORMULATION IMPLIES SATURATION)
!     -------------------------------------------------------

  REAL(dp), PARAMETER :: RHCDD       = 1._dp

!     RMFCFL:     MASSFLUX MULTIPLE OF CFL STABILITY CRITERIUM
!     -------

  REAL(dp), PARAMETER :: RMFCFL      = 1._dp

!     RMFSHCFL:   MASSFLUX MULTIPLE OF CFL FOR SHALLOW IF dt>1800 s,
!     -------     "TEMPORARY SOLUTION", ( DEFAULT=1)

  REAL(dp), PARAMETER :: RMFSHCFL    = 1._dp

!     MASSFLUX SOLVERs FOR MOMEMTUM AND TRACERS
!     0: EXPLICIT 0-1 SEMI-IMPLICIT >=1: IMPLICIT
!     -------------------------------------------

  REAL(dp), PARAMETER :: RMFSOLUV   = 0._dp  ! mass flux solver for momentum
  REAL(dp), PARAMETER :: RMFSOLCT   = 0._dp  ! mass flux solver for chemical tracers 


  REAL(dp) :: RTAU  
  LOGICAL  :: LMFPEN
  LOGICAL  :: LMFSCV
  LOGICAL  :: LMFMID
  LOGICAL  :: LMFDD
  LOGICAL  :: LMFDUDV
  LOGICAL  :: LMFTRAC
!*UPG change to operations 
  LOGICAL  :: LEPCLD
  LOGICAL  :: LMFSCL_WSTAR
!*UPG change to operations
  INTEGER  :: NJKT1, NJKT2, NJKT3

!*    ** *YOECUMF* - PARAMETERS FOR CUMULUS MASSFLUX SCHEME

!     M.TIEDTKE       E. C. M. W. F.      18/1/89

!     NAME      TYPE      PURPOSE
!     ----      ----      -------

!     LMFPEN    LOGICAL  TRUE IF PENETRATIVE CONVECTION IS SWITCHED ON
!     LMFSCV    LOGICAL  TRUE IF SHALLOW     CONVECTION IS SWITCHED ON
!     LMFMID    LOGICAL  TRUE IF MIDLEVEL    CONVECTION IS SWITCHED ON
!     LMFDD     LOGICAL  TRUE IF CUMULUS DOWNDRAFT      IS SWITCHED ON
!     LMFDUDV   LOGICAL  TRUE IF CUMULUS FRICTION       IS SWITCHED ON
!     LMFTRAC   LOGICAL  TRUE IF TRACER TRANSPORT       IS SWITCHED ON
!     LMFSCL_WSTAR LOG   TRUE IF W* SHALLOW CLOSURE     IS SWITCHED ON
!     ENTRPEN   REAL     ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
!     ENTRSCV   REAL     ENTRAINMENT RATE FOR SHALLOW CONVECTION
!     ENTRMID   REAL     ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
!     ENTRDD    REAL     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     RMFCTOP   REAL     RELAT. CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANC
!     RMFCMAX   REAL     MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     RMFCMIN   REAL     MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     RMFDEPS   REAL     FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     RDEPTHS   REAL     MAXIMUM ALLOWED CLOUD THICKNESS FOR SHALLOW
!     RHCDD     REAL     RELATIVE SATURATION IN DOWNDRAFTS
!     RPRCON    REAL     COEFFICIENTS FOR DETERMINING CONVERSION
!                        FROM CLOUD WATER TO RAIN
!     RCPECONS  REAL     COEFFICIENT FOR RAIN EVAPORATION BELOW CLOUD
!     RCUCOV    REAL     CONVECTIVE CLOUD COVER FOR RAIN EVPORATION
!     RTAUMEL   REAL     TIME CONSTANT FOR MELTING
!     RHEBC     REAL     REL. HUMIDITY BELOW CLOUD FOR WHICH EVAPORATION STARTS
!     RTAU      REAL     ADJUSTMENT TIME SCALE IN CAPE CLOSURE
!     RMFCFL    REAL     MASSFLUX MULTIPLE OF CFL STABILITY CRITERIUM 
!     RMFSHCFL  REAL     MASSFLUX MULTIPLE OF CFL CRITERIUM FOR SHALLOW
!     RMFSOLUV  REAL     SOLVER FOR MASSFLUX ADVECTION EQUATION FOR MOMENTUM
!     RMFSOLCT  REAL     SOLVER FOR MASSFLUX ADVECTION EQUATION FOR TRACERS
!                        0 : EXPLICIT  0-1 : SEMI-IMPLICIT >=1 : FULLY IMPLICIT
!     NJKT1, NJKT2, NJKT3 INTEGER  LEVEL LIMITS FOR CUBASEN/CUDDR
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!*    END OF *YOECUMF* - PARAMETERS FOR CUMULUS MASSFLUX SCHEME
!     ----------------------------------------------------------------


!====================================================================================

!     ------------------------------------------------------------------
!*    ** *YOEPHLI* CONTAINS CONSTANTS FOR THE LINEARIZED PHYSICS
!     ------------------------------------------------------------------

  LOGICAL  :: LPHYLIN = .FALSE.
!  LOGICAL LENOPERT
!  LOGICAL LRAISANEN

  REAL(dp) :: RLPTRC
  REAL(dp) :: RLPAL1
  REAL(dp) :: RLPAL2
!  REAL(dp) :: RLPBB
!  REAL(dp) :: RLPCC
!  REAL(dp) :: RLPDD
!  REAL(dp) :: RLPMIXL
!  REAL(dp) :: RLPBETA
!  REAL(dp) :: RLPDRAG
!  REAL(dp) :: RLPEVAP
!  REAL(dp) :: RLPP00

!*     *YOEPHLI* CONTAINS CONSTANTS NEEDED BY 
!     THE LINEARIZED PHYSICS


!     J.F. MAHFOUF        E.C.M.W.F.    23/06/96


!     NAME        TYPE     DESCRIPTION
!     ----        ----     -----------

!     *RLPTRC*    REAL     CRITICAL TEMPERATURE FOR MIXED PHASE PROPERTIES
!                          OF WATER 
!     *RLPAL1*    REAL     SMOOTHING COEFFICIENT
!     *RLPAL2*    REAL     SMOOTHING COEFFICIENT
!     *RLPBB*     REAL     CONSTANT FROM THE LOUIS ET AL. FORMULATION
!     *RLPCC*     REAL     CONSTANT FROM THE LOUIS ET AL. FORMULATION
!     *RLPDD*     REAL     CONSTANT FROM THE LOUIS ET AL. FORMULATION
!     *RLPMIXL*   REAL     PSEUDO DEPTH OF THE PLANETARY BOUNDARY LAYER
!     *RLPBETA*   REAL     REDUCTION FACTOR OF THE ASYMPTOTIC MIXING LENGTH
!     *RLPDRAG*   REAL     COEFFICIENT FOR THE ESTIMATION OF SURFACE DRAG
!     *RLPEVAP*   REAL     FRACTION OF POSSIBLE RAINFALL EVAPORATION
!     *RLPP00*    REAL     PRESSURE ABOVE WHICH RADIATION IS NOT APPLIED
!     *LPHYLIN*   LOGICAL  TRUE WHEN LINEARIZED PHYSICS IS ACTIVATED 
!     *LENOPERT   LOGICAL  TRUE WHEN NO PERTURBATION IS REQUIRED
!                          FOR SURFACE ARRAYS
!     *LRAISANEN   LOGICAL  TRUE WHEN RAISANEN OVERLAP SCHEME IS
!                           ACTIVATED
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!*    *END OF *YOEPHLI* CONTAINS CONSTANTS FOR THE LINEARIZED PHYSICS
!     ------------------------------------------------------------------

  INTEGER :: NSMAX 
  INTEGER :: N_VMASS=0
  LOGICAL :: LHOOK=.FALSE.
  
!====================================================================================
!     -----------------------------------------------------------------
!     ** YOECLDP - CONTROL PARAMETERS FOR PROGNOSTIC CLOUD SCHEME
!     -----------------------------------------------------------------

!     * E.C.M.W.F. PHYSICS PACKAGE *

!     C. JAKOB        E.C.M.W.F.          94/02/07

  REAL(dp), PARAMETER :: RLMIN = 1.E-8_dp        !     *RLMIN*   REAL      LIMIT FOR L

!     -----------------------------------------------------------------
!     END OF YOECLDP - CONTROL PARAMETERS FOR PROGNOSTIC CLOUD SCHEME
!     -----------------------------------------------------------------

!=========================================================================================

CONTAINS

!=========================================================================================
  SUBROUTINE DR_HOOK(CH,K1,P1)

    IMPLICIT NONE

    CHARACTER :: CH
    INTEGER   :: K1
    REAL(dp)  :: P1

  END SUBROUTINE DR_HOOK

!-----------------------------------------------------------------------------------------
  SUBROUTINE VEXP(P1,P2,K1)

    IMPLICIT NONE

    INTEGER,INTENT(IN)   :: K1
    REAL(dp),INTENT(OUT) :: P1(:)
    REAL(dp),INTENT(IN)  :: P2(:)

    INTEGER :: i
    INTRINSIC :: EXP

    do i=1,k1
      P1(i) = exp(P2(i))
    enddo

  END SUBROUTINE VEXP

!-----------------------------------------------------------------------------------------

  SUBROUTINE VDIV(P1,P2,P3,K1)

    IMPLICIT NONE
    INTEGER,INTENT(IN)   :: K1
    REAL(dp),INTENT(IN)  :: P2(:),P3(:)
    REAL(dp),INTENT(OUT) :: P1(:)

    INTEGER :: i

    do i=1,K1
      P1(i) = P2(i) / P3(i)
    enddo

  END SUBROUTINE VDIV

!=========================================================================================
  SUBROUTINE init_convection_constants(nnsmax, nflevg)

    IMPLICIT NONE

    INTEGER :: nnsmax, nflevg
    

!     LOGICAL SWITCHES
!     ----------------
!     SET ADJUSTMENT TIME SCALE FOR CAPE CLOSURE AS A FUNCTION
!     OF MODEL RESOLUTION
!     RTAU IS 20 MINUTES FOR RESOLUTIONS HIGHER THAN TL319
!     RTAU IS 1 HOUR FOR ANY OTHER RESOLUTION
!     --------------------------------------------------------

    NSMAX = nnsmax

    IF(NSMAX < 63) THEN
      RTAU= 14400._dp
    ELSEIF (NSMAX > 63) THEN
      RTAU= 7200.0_dp
    ELSEIF (NSMAX > 159) THEN
      RTAU= 3600.0_dp
    ELSEIF (NSMAX > 319) THEN
      RTAU= 1200.0_dp
      ! op_mm_20140123+
   ELSEIF (NSMAX < -0) THEN  !SPECIAL CASE to distinguish COSMO (nn<0)
      RTAU=120.0_dp
      ! op_mm_20140123-
    ENDIF

!!    LMFPEN  =.TRUE.   ! deep convection
!!    LMFSCV  =.TRUE.   ! shallow convection
!!    LMFMID  =.TRUE.   ! mid-level convection
!!    LMFDD   =.TRUE.   ! use downdrafts
!!    LMFDUDV =.TRUE.   ! use convective momentum transport
!!    LMFTRAC =.TRUE.   ! convective chemical tracer transport

!*UPG add to operations
!!    LEPCLD  =.TRUE.      ! produce detrained cloud water/ice
!!                         ! to reuse in prognostic cloud scheme
!!    LMFSCL_WSTAR=.FALSE. ! use w* shallow convective closure proposed by Grant.
!!                         ! this is useful if your host model cannot provide
!!                         ! T and q tendencies for the boundary-layer
!*UPG add to operations



    NJKT1=2
    NJKT2=2
    NJKT3=NFLEVG-2
    !DO JLEV=NFLEVG,2,-1
    !  IF(STPRE(JLEV) > 350.E2_dp)NJKT1=JLEV
    !  IF(STPRE(JLEV) >  60.E2_dp)NJKT2=JLEV
    !  IF(STPRE(JLEV) > 950.E2_dp)NJKT3=JLEV
    !ENDDO
    !NJKT3=MIN(NFLEVG-2,NJKT3)

  END SUBROUTINE init_convection_constants

!=======================================================================================

  SUBROUTINE SATUR ( KIDIA , KFDIA , KLON  , KTDIA , KLEV,&
                   & PAPRSF, PT    , PQSAT , KFLAG)  

!***

! **   *SATUR* -  COMPUTES SPECIFIC HUMIDITY AT SATURATION

!       J.F. MAHFOUF       E.C.M.W.F.     15/05/96

!       Modified J. HAGUE          13/01/03 MASS Vector Functions       

!       PURPOSE.
!       --------

!       SPECIFIC HUMIDITY AT SATURATION IS USED BY THE
!       DIAGNOSTIC CLOUD SCHEME TO COMPUTE RELATIVE HUMIDITY
!       AND LIQUID WATER CONTENT  

!       INTERFACE
!       ---------

!       THIS ROUTINE IS CALLED FROM *CALLPAR*.

!       PARAMETER     DESCRIPTION                                 UNITS
!       ---------     -----------                                 -----
!       INPUT PARAMETERS (INTEGER):

!      *KIDIA*        START POINT
!      *KFDIA*        END POINT
!      *KLON*         NUMBER OF GRID POINTS PER PACKET
!      *KTDIA*        START OF THE VERTICAL LOOP
!      *KLEV*         NUMBER OF LEVELS

!       INPUT PARAMETERS (REAL):

!      *PAPRSF*        PRESSURE ON FULL LEVELS                      PA
!      *PT*            TEMPERATURE AT T-DT                          K

!       INPUT PARAMETERS (INTEGER):

!      *KFLAG*         FLAG TO DETECT CALL FROM

!                      CONVECTION  KFLAG=1
!                      OTHER       KFLAG=2

!       OUTPUT PARAMETER (REAL):

!      *PQSAT*         SATURATION SPECIFIC HUMIDITY                 KG/KG

!-------------------------------------------------------------------------

IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER,INTENT(IN)    :: KTDIA 
REAL(dp)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PQSAT(KLON,KLEV) 
INTEGER,INTENT(IN)    :: KFLAG 
INTEGER :: JK, JL, JLEN

REAL(dp) :: Z3ES, Z4ES, ZCOR, ZEW, ZFOEEW, ZQMAX, ZQS, ZTARG
REAL(dp) :: Z_EXPARG1(KIDIA:KFDIA+N_VMASS)
REAL(dp) :: Z_EXPARG2(KIDIA:KFDIA+N_VMASS)
REAL(dp) :: Z_EXPOUT1(KIDIA:KFDIA+N_VMASS)
REAL(dp) :: Z_EXPOUT2(KIDIA:KFDIA+N_VMASS)
REAL(dp) :: ZHOOK_HANDLE

INTRINSIC :: EXP, MIN, MOD
!----------------------------------------------------------------------

!*    1.           DEFINE CONSTANTS
!                  ----------------

IF (LHOOK) CALL DR_HOOK('SATUR',0,ZHOOK_HANDLE)
ZQMAX=0.5_dp

!     *
!----------------------------------------------------------------------

!     *    2.           CALCULATE SATURATION SPECIFIC HUMIDITY
!                       --------------------------------------

IF (LPHYLIN) THEN
  DO JK=KTDIA,KLEV
    DO JL=KIDIA, KFDIA
      ZTARG = PT(JL,JK)
      IF (ZTARG > RTT) THEN
        Z3ES=R3LES
        Z4ES=R4LES
      ELSE
        Z3ES=R3IES
        Z4ES=R4IES
      ENDIF
      ZFOEEW = R2ES*EXP(Z3ES*(ZTARG-RTT)/(ZTARG-Z4ES))
      ZQS    = ZFOEEW/PAPRSF(JL,JK)
      IF (ZQS > ZQMAX) THEN
        ZQS=ZQMAX
      ENDIF
      ZCOR = 1.0_dp/(1.0_dp-RETV*ZQS)
      PQSAT(JL,JK)=ZQS*ZCOR
    ENDDO
  ENDDO
ELSE

  IF(N_VMASS <= 0) THEN ! Not using Vector MASS

    DO JK=KTDIA,KLEV
      DO JL=KIDIA, KFDIA
        IF(KFLAG == 1) THEN
          ZEW  = FOEEWMCU(PT(JL,JK))
        ELSE
          ZEW  = FOEEWM(PT(JL,JK))
        ENDIF
        ZQS  = ZEW/PAPRSF(JL,JK)
        ZQS  = MIN(ZQMAX,ZQS)
        ZCOR = 1.0_dp/(1.0_dp-RETV*ZQS)
        PQSAT(JL,JK)=ZQS*ZCOR
      ENDDO
    ENDDO

  ELSE ! Using Vector MASS

    JLEN=KFDIA-KIDIA+N_VMASS-MOD(KFDIA-KIDIA,N_VMASS)

    IF(KFDIA-KIDIA+1 /= JLEN) THEN
      Z_EXPARG1(KFDIA+1:KIDIA+JLEN-1)=1.0_dp     
      Z_EXPARG2(KFDIA+1:KIDIA+JLEN-1)=1.0_dp     
    ENDIF

    DO JK=KTDIA,KLEV
      DO JL=KIDIA, KFDIA
        Z_EXPARG1(JL)=FOELES_V(PT(JL,JK))
        Z_EXPARG2(JL)=FOEIES_V(PT(JL,JK))
      ENDDO
      CALL VEXP(Z_EXPOUT1,Z_EXPARG1,JLEN)
      CALL VEXP(Z_EXPOUT2,Z_EXPARG2,JLEN)
      DO JL=KIDIA, KFDIA
        IF(KFLAG == 1) THEN
          ZEW  = FOEEWMCU_V( PT(JL,JK),Z_EXPOUT1(JL),Z_EXPOUT2(JL) )
        ELSE
          ZEW  = FOEEWM_V  ( PT(JL,JK),Z_EXPOUT1(JL),Z_EXPOUT2(JL) )
        ENDIF
!       ZQS  = ZEW/PAPRSF(JL,JK)
!       ZQS  = MIN(ZQMAX,ZQS)
!!      ZCOR = 1._dp/(1._dp-RETV*ZQS)
!       PQSAT(JL,JK)=ZQS/(1._dp-RETV*ZQS)
        ZQS  = MIN(ZQMAX*PAPRSF(JL,JK),ZEW)
        PQSAT(JL,JK)=ZQS/(PAPRSF(JL,JK)-RETV*ZQS)
      ENDDO
    ENDDO
  
  ENDIF

ENDIF

IF (LHOOK) CALL DR_HOOK('SATUR',1,ZHOOK_HANDLE)
END SUBROUTINE SATUR

!==========================================================================

! content from fcttre.h 
! transformed from statement fuctions into ELEMETAL FUNCTIONS

!*
!     ------------------------------------------------------------------

!     This COMDECK includes the Thermodynamical functions for the cy39
!       ECMWF Physics package.
!       Consistent with YOMCST Basic physics constants, assuming the
!       partial pressure of water vapour is given by a first order
!       Taylor expansion of Qs(T) w.r.t. to Temperature, using constants
!       in YOETHF
!       Two sets of functions are available. In the first set only the
!       cases water or ice are distinguished by temperature.  This set 
!       consists of the functions FOEDELTA,FOEEW,FOEDE and FOELH.
!       The second set considers, besides the two cases water and ice 
!       also a mix of both for the temperature range RTICE < T < RTWAT.
!       This set contains FOEALFA,FOEEWM,FOEDEM,FOELDCPM and FOELHM.

!       Depending on the consideration of mixed phases either the first 
!       set (e.g. surface, post-processing) or the second set 
!       (e.g. clouds, condensation, convection) should be used.

!     ------------------------------------------------------------------
!     *****************************************************************

!                NO CONSIDERATION OF MIXED PHASES

!     *****************************************************************

!----------------------------------------------------------------------------------  
  ELEMENTAL REAL(dp) FUNCTION  FOEDELTA(PTARE)

    REAL(dp), INTENT(IN) :: PTARE
    INTRINSIC :: MAX, SIGN
    
    FOEDELTA  = MAX (0.0_dp,SIGN(1.0_dp,PTARE-RTT))
    

!                  FOEDELTA = 1    water
!                  FOEDELTA = 0    ice
  END FUNCTION FOEDELTA
!     THERMODYNAMICAL FUNCTIONS .
!----------------------------------------------------------------------------------

!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE

  ELEMENTAL REAL(dp) FUNCTION FOEEW(PTARE)
    REAL(dp), INTENT(IN) :: PTARE
    INTRINSIC :: EXP

    FOEEW  = R2ES*EXP ((R3LES * FOEDELTA(PTARE) +                   &
             R3IES * (1.0_dp - FOEDELTA(PTARE))) * (PTARE - RTT) /  &
             (PTARE - (R4LES * FOEDELTA(PTARE) +                    & 
             R4IES * (1.0_dp - FOEDELTA(PTARE)))))
  END FUNCTION FOEEW

  ELEMENTAL REAL(dp) FUNCTION FOEDE(PTARE)
    REAL(dp), INTENT(IN) :: PTARE

    FOEDE  = (FOEDELTA(PTARE) * R5ALVCP +                   &
             (1.0_dp - FOEDELTA(PTARE)) * R5ALSCP) /        &
             (PTARE - (R4LES * FOEDELTA(PTARE) +            &
             R4IES * (1.0_dp - FOEDELTA(PTARE))))**2
  END FUNCTION FOEDE

  ELEMENTAL REAL(dp) FUNCTION FOEDESU(PTARE)
    REAL(dp), INTENT(IN) :: PTARE

    FOEDESU = (FOEDELTA(PTARE) * R5LES +                    &
              (1.0_dp - FOEDELTA(PTARE)) * R5IES) /         &
              (PTARE - (R4LES * FOEDELTA(PTARE) +           &
              R4IES*(1.0_dp-FOEDELTA(PTARE))))**2
  END FUNCTION FOEDESU

  ELEMENTAL REAL(dp) FUNCTION FOELH(PTARE)
    REAL(dp), INTENT(IN) :: PTARE

    FOELH   = FOEDELTA(PTARE) * RLVTT +           &
              (1.0_dp - FOEDELTA(PTARE)) * RLSTT
  END FUNCTION FOELH

  ELEMENTAL REAL(dp) FUNCTION FOELDCP(PTARE)
    REAL(dp), INTENT(IN) :: PTARE

    FOELDCP = FOEDELTA(PTARE) * RALVDCP +         &
              (1.0_dp - FOEDELTA(PTARE)) * RALSDCP
  END FUNCTION FOELDCP

!     *****************************************************************

!           CONSIDERATION OF MIXED PHASES

!     *****************************************************************

!     FOEALFA is calculated to distinguish the three cases:

!                       FOEALFA=1            water phase
!                       FOEALFA=0            ice phase
!                       0 < FOEALFA < 1      mixed phase

!               INPUT : PTARE = TEMPERATURE
  ELEMENTAL REAL(dp) FUNCTION FOEALFA(PTARE)
    REAL(dp), INTENT(IN) :: PTARE
    INTRINSIC :: MIN, MAX

    FOEALFA  = MIN(1.0_dp,((MAX(RTICE, MIN(RTWAT,PTARE)) - RTICE) * &
               RTWAT_RTICE_R)**2) 

  END FUNCTION FOEALFA


!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE
  ELEMENTAL REAL(dp) FUNCTION FOEEWM(PTARE)
    REAL(dp), INTENT(IN) :: PTARE
    INTRINSIC :: EXP

    FOEEWM   = R2ES * (FOEALFA(PTARE) * EXP (R3LES *   &
               (PTARE - RTT) / (PTARE - R4LES)) +      &
               (1.0_dp - FOEALFA(PTARE)) * EXP(R3IES * &
               (PTARE - RTT) / (PTARE - R4IES)))

  END FUNCTION FOEEWM

  ELEMENTAL REAL(dp) FUNCTION FOEDEM(PTARE)
    REAL(dp), INTENT(IN) :: PTARE

    FOEDEM   = FOEALFA(PTARE) * R5ALVCP * (1.0_dp / (PTARE - R4LES)**2) + &
               (1.0_dp - FOEALFA(PTARE)) * R5ALSCP *                      &
               (1.0_dp / (PTARE - R4IES)**2)

  END FUNCTION FOEDEM

  ELEMENTAL REAL(dp) FUNCTION FOELDCPM(PTARE)
    REAL(dp), INTENT(IN) :: PTARE

    FOELDCPM  = FOEALFA(PTARE) * RALVDCP +           &
               (1.0_dp - FOEALFA(PTARE)) * RALSDCP

  END FUNCTION FOELDCPM

  ELEMENTAL REAL(dp) FUNCTION FOELHM(PTARE)
    REAL(dp), INTENT(IN) :: PTARE

    FOELHM    = FOEALFA(PTARE) * RLVTT +             &
                (1.0_dp - FOEALFA(PTARE)) * RLSTT

  END FUNCTION FOELHM

!     Temperature normalization for humidity background change of variable
!        INPUT : PTARE = TEMPERATURE

  ELEMENTAL REAL(dp) FUNCTION FOETB(PTARE)
    REAL(dp), INTENT(IN) :: PTARE

    FOETB     = FOEALFA(PTARE) * R3LES * (RTT - R4LES) *   &
                (1.0_dp / (PTARE - R4LES)**2) +            &
                (1.0_dp - FOEALFA(PTARE)) * R3IES *        &
                (RTT - R4IES) * (1.0_dp / (PTARE - R4IES)**2)

  END FUNCTION FOETB

!     ------------------------------------------------------------------
!     *****************************************************************

!           CONSIDERATION OF DIFFERENT MIXED PHASE FOR CONV

!     *****************************************************************

!     FOEALFCU is calculated to distinguish the three cases:

!                       FOEALFCU=1            water phase
!                       FOEALFCU=0            ice phase
!                       0 < FOEALFCU < 1      mixed phase

!               INPUT : PTARE = TEMPERATURE
    ELEMENTAL REAL(dp) FUNCTION FOEALFCU(PTARE)
      REAL(dp), INTENT(IN) :: PTARE
      INTRINSIC :: MIN, MAX

      FOEALFCU = MIN(1.0_dp,((MAX(RTICECU, &
                 MIN(RTWAT,PTARE)) - RTICECU) * RTWAT_RTICECU_R)**2) 

    END FUNCTION FOEALFCU


!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE

  ELEMENTAL REAL(dp) FUNCTION FOEEWMCU(PTARE)
    REAL(dp), INTENT(IN) :: PTARE
    INTRINSIC :: EXP
    
    FOEEWMCU   = R2ES * (FOEALFCU(PTARE) *                      &
                 EXP(R3LES * (PTARE - RTT) / (PTARE - R4LES)) + &
                 (1.0_dp - FOEALFCU(PTARE)) *                   &
                 EXP(R3IES * (PTARE - RTT) / (PTARE - R4IES)))

  END FUNCTION FOEEWMCU

  ELEMENTAL REAL(dp) FUNCTION FOEDEMCU(PTARE)
    REAL(dp), INTENT(IN) :: PTARE

    FOEDEMCU   = FOEALFCU(PTARE) * R5ALVCP *            &
                 (1.0_dp / (PTARE - R4LES)**2) +        &
                 (1.0_dp - FOEALFCU(PTARE)) * R5ALSCP * &
                 (1.0_dp / (PTARE - R4IES)**2)

  END FUNCTION FOEDEMCU

  ELEMENTAL REAL(dp) FUNCTION FOELDCPMCU(PTARE)
    REAL(dp), INTENT(IN) :: PTARE

    FOELDCPMCU = FOEALFCU(PTARE) * RALVDCP +          &
                 (1.0_dp - FOEALFCU(PTARE)) * RALSDCP

  END FUNCTION FOELDCPMCU

  ELEMENTAL REAL(dp) FUNCTION FOELHMCU(PTARE)
    REAL(dp), INTENT(IN) :: PTARE

    FOELHMCU   = FOEALFCU(PTARE) * RLVTT +            &
                 (1.0_dp - FOEALFCU(PTARE)) * RLSTT

  END FUNCTION FOELHMCU

!     ------------------------------------------------------------------
!     Pressure of water vapour at saturation
!     This one is for the WMO definition of saturation, i.e. always
!     with respect to water.

  ELEMENTAL REAL(dp) FUNCTION FOEEWMO(PTARE)
    REAL(dp), INTENT(IN) :: PTARE
    INTRINSIC :: EXP

    FOEEWMO    = R2ES * EXP(R3LES * (PTARE - RTT) / (PTARE - R4LES))

  END FUNCTION FOEEWMO

  ELEMENTAL REAL(dp) FUNCTION FOELES_V(PTARE)
    REAL(dp), INTENT(IN) :: PTARE

    FOELES_V   = R3LES * (PTARE - RTT) / (PTARE - R4LES)

  END FUNCTION FOELES_V

  ELEMENTAL REAL(dp) FUNCTION FOEIES_V(PTARE)
    REAL(dp), INTENT(IN) :: PTARE

    FOEIES_V   = R3IES * (PTARE - RTT) / (PTARE - R4IES)
    
  END FUNCTION FOEIES_V

  ELEMENTAL REAL(dp) FUNCTION FOEEWM_V(PTARE, EXP1, EXP2)
    REAL(dp), INTENT(IN) :: PTARE, EXP1, EXP2

    FOEEWM_V   = R2ES * (FOEALFA(PTARE) * EXP1 +    &
                 (1.0_dp - FOEALFA(PTARE)) * EXP2)

  END FUNCTION FOEEWM_V

  ELEMENTAL REAL(dp) FUNCTION FOEEWMCU_V(PTARE, EXP1, EXP2)
    REAL(dp), INTENT(IN) :: PTARE, EXP1, EXP2

    FOEEWMCU_V = R2ES * (FOEALFCU(PTARE) * EXP1 +   &
                 (1.0_dp - FOEALFCU(PTARE)) * EXP2)
  END FUNCTION FOEEWMCU_V

!===============================================================================

    END MODULE MESSY_CONVECT_ECMWF_PARAM
