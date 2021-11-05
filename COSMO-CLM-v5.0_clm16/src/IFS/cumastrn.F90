SUBROUTINE CUMASTRN &
 & (  KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & KSTEP,    KSTART,   LDLAND,   PTSPHY,&
 & PTEN,     PQEN,     PUEN,     PVEN,     PLITOT,&
 & PVERVEL,  PQHFL,    PAHFS,&
 & psstru,   psstrv,&
 & PAP,      PAPH,     PGEO,     PGEOH,&
 & PTENT,    PTENQ,    PTENU,    PTENV,&
 & PTENL,    PTENI, &
 & LDCUM,    KTYPE,    KCBOT,    KCTOP,&
 & KBOTSC,   LDSC,&
 & PTU,      PQU,      PLU,&
 & PMFLXR,   PMFLXS,   PRAIN,&
 & PMFU,     PMFD,&
 & PMFUDE_RATE,        PMFDDE_RATE,    PCAPE,&
 & KTRAC,    PCEN,     PTENC )  

!**** *CUMASTR*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME

!     M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
!     D.GREGORY      E.C.M.W.F.     1996
!     P.BECHTOLD     E.C.M.W.F.     2005/2007

!     PURPOSE
!     -------

!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
!     PROGNOSTIC VARIABLES T,Q,U, V AND TRACERS DUE TO CONVECTIVE PROCESSES.
!     PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
!     PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
!     SATURATED CUMULUS DOWNDRAFTS.

!**   INTERFACE.
!     ----------

!          *CUMASTR* IS CALLED FROM *CUCALL*
!     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
!     T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
!     IT RETURNS ITS OUTPUT TO THE SAME SPACE
!      1.MODIFIED TENDENCIES OF MODEL VARIABLES
!      2.RATES OF CONVECTIVE PRECIPITATION
!        (USED IN SUBROUTINE SURF)
!      3.CLOUD BASE, CLOUD TOP AND PRECIP FOR RADIATION
!        (USED IN SUBROUTINE CLOUD)

!     METHOD
!     -------

!     PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
!        (1) DEFINE CONSTANTS AND PARAMETERS
!        (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
!            INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
!        (3) CALCULATE CLOUD BASE IN 'CUBASE'
!            AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
!        (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
!        (5) DO DOWNDRAFT CALCULATIONS:
!              (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
!              (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
!              (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
!                  EFFECT OF CU-DOWNDRAFTS
!        (6) DO FINAL CLOUD ASCENT IN 'CUASC'
!        (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
!            DO EVAPORATION IN SUBCLOUD LAYER
!        (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
!        (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KTRAC*        NUMBER OF CHEMICAL TRACERS
!    *KSTEP*        CURRENT TIME STEP INDEX
!    *KSTART*       FIRST STEP OF MODEL

!     INPUT PARAMETERS (LOGICAL)

!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)

!     INPUT PARAMETERS (REAL)

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PCEN*         PROVISIONAL ENVIRONMENT TRACER CONCENTRATIONS KG/KG 
!    *PLITOT*       GRID MEAN LIQUID WATER+ICE CONTENT            KG/KG
!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS             PA
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PSSTRU*       SURFACE MOMENTUM FLUX U                - not used presently
!    *PSSTRV*       SURFACE MOMENTUM FLUX V                - not used presently 

!    UPDATED PARAMETERS (REAL):

!    *PTENT*        TEMPERATURE TENDENCY                           K/S
!    *PTENQ*        MOISTURE TENDENCY                             KG/(KG S)
!    *PTENL*        LIQUID WATER TENDENCY                         KG/(KG S)
!    *PTENI*        ICE CONDENSATE TENDENCY                       KG/(KG S)
!    *PTENU*        TENDENCY OF U-COMP. OF WIND                    M/S2
!    *PTENV*        TENDENCY OF V-COMP. OF WIND                    M/S2
!    *PTENC*        TENDENCY OF CHEMICAL TRACERS                   1/S

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDSC*         FLAG: .TRUE. FOR SC-POINTS

!    OUTPUT PARAMETERS (INTEGER):

!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP*        CLOUD TOP LEVEL
!    *KBOTSC*       CLOUD BASE LEVEL FOR SC-CLOUDS

!    OUTPUT PARAMETERS (REAL):

!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PLUDE*        DETRAINED LIQUID WATER                        KG/(M2*S)
!    *PENTH*        INCREMENT OF DRY STATIC ENERGY                 J/(KG*S)
!    *PMFLXR*       CONVECTIVE RAIN FLUX                          KG/(M2*S)
!    *PMFLXS*       CONVECTIVE SNOW FLUX                          KG/(M2*S)
!    *PRAIN*        TOTAL PRECIP. PRODUCED IN CONV. UPDRAFTS      KG/(M2*S)
!                   (NO EVAPORATION IN DOWNDRAFTS)
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PMFUDE_RATE*  UPDRAFT DETRAINMENT RATE                      KG/(M3*S)
!    *PMFDDE_RATE*  DOWNDRAFT DETRAINMENT RATE                    KG/(M3*S)
!    *PCAPE*        CONVECTVE AVAILABLE POTENTIAL ENERGY           J/KG
!    *PWMEAN*       VERTICALLY AVERAGED UPDRAUGHT VELOCITY         M/S

!     EXTERNALS.
!     ----------

!       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
!       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
!       CUASC:  CLOUD ASCENT FOR ENTRAINING PLUME
!       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
!       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
!       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
!       CUDQDT: UPDATES TENDENCIES FOR T AND Q
!       CUDUDV: UPDATES TENDENCIES FOR U AND V

!     SWITCHES.
!     --------

!          LMFPEN=.TRUE.   PENETRATIVE CONVECTION IS SWITCHED ON
!          LMFSCV=.TRUE.   SHALLOW CONVECTION IS SWITCHED ON
!          LMFMID=.TRUE.   MIDLEVEL CONVECTION IS SWITCHED ON
!          LMFIT=.TRUE.    UPDRAUGHT ITERATION
!          LMFDD=.TRUE.    CUMULUS DOWNDRAFTS SWITCHED ON
!          LMFDUDV=.TRUE.  CUMULUS FRICTION SWITCHED ON
!          LMFTRAC=.false. TRACER TRANSPORT 

!     MODEL PARAMETERS (DEFINED IN SUBROUTINE CUPARAM)
!     ------------------------------------------------
!     ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
!     ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
!     ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
!     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     RMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     RMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     RMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     RPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN

!     REFERENCE.
!     ----------

!          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
!          DRAFT PAPER ON MASSFLUX SCHEME (NORDENG, 1995)

!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!             96-03-12 : Introduce CAPE closure for deep convection
!                        (based upon previous work by Nordeng, 1995)
!             99-06-21 : Optimisation   D.Salmond
!             03-08-29 : Clean-up, deep/shallow switches  P.Bechtold
!             04-02-11 : Add tracer transport             P.Bechtold
!             05-02-11 : Positive scaling of total Mflux  P.Bechtold
!             M.Hamrud      01-Oct-2003 CY28 Cleaning
!             04-12-03 : Turn off shallow convection over stratocu.
!                                                         M.Ko"hler
!             05-06-27 : Switch off ddraught if idtop<kctop  
!                        correction for detrainment rates P.Bechtold
!             05-11-22 : Mods for coarser/finer physics D.Salmond + M.Hortal
!             06-02-11 : Enable TQ implicit               P.Bechtold
!             07-06-01 : Only single updraught call with  P.Bechtold
!                        scaling, convective turnover time
!                        scale, contain momentum computations
!                        in cumastrn
!----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK_IFS   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RG       ,RD       ,RCPD     ,RETV     ,&
 & RLVTT    ,RLSTT   ,RTT      ,RPI
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
 & RTWAT_RTICECU_R    ,RTWAT_RTICE_R  
USE YOECUMF  , ONLY : ENTRPEN  ,ENTRSCV  ,LMFDD    ,LMFDUDV  ,&
 & RTAU      ,RDEPTHS ,LMFSCV  ,LMFPEN   ,LMFIT    ,NJKT2    ,&
 & RMFCFL    ,RMFLIC  ,RMFLIA  ,RMFSOLUV ,NJKT4    ,NJKT5    ,&
 & RUVPER    ,RMFSOLTQ,RMFSOLCT,RMFCMIN  ,LMFSMOOTH,LMFWSTAR ,&
 & LMFTRAC
!* change to operations
!USE YOEPHY   , ONLY : LMFTRAC

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM)               :: KTDIA ! Argument NOT used
INTEGER(KIND=JPIM)               :: KSTEP ! Argument NOT used
INTEGER(KIND=JPIM)               :: KSTART ! Argument NOT used
LOGICAL           ,INTENT(IN)    :: LDLAND(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLITOT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
!REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQSEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQHFL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAHFS(KLON,KLEV+1) 
REAL(KIND=JPRB)                  :: PSSTRU(KLON) ! Argument NOT used
REAL(KIND=JPRB)                  :: PSSTRV(KLON) ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCEN(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC) 
LOGICAL           ,INTENT(INOUT) :: LDCUM(KLON)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KTYPE(KLON)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCBOT(KLON)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCTOP(KLON)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KBOTSC(KLON) 
LOGICAL           ,INTENT(OUT)   :: LDSC(KLON) 
!LOGICAL           ,INTENT(IN)    :: LDSHCV(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLU(KLON,KLEV) 
!REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLUDE(KLON,KLEV) 
!REAL(KIND=JPRB)   ,INTENT(OUT)   :: PENTH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMFLXR(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMFLXS(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAIN(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMFU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMFD(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMFUDE_RATE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMFDDE_RATE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPE(KLON) 
!REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWMEAN(KLON) 
!*UPG change to operations
REAL(KIND=JPRB) :: PWMEAN(KLON) 
REAL(KIND=JPRB) :: PLUDE(KLON,KLEV) ! only local variable
REAL(KIND=JPRB) :: PENTH(KLON,KLEV) ! only local variable
REAL(KIND=JPRB) :: PQSEN(KLON,KLEV) ! only local variable
!*UPG


REAL(KIND=JPRB) ::     ZTENH(KLON,KLEV),       ZQENH(KLON,KLEV),&
 & ZQSENH(KLON,KLEV),&
 & ZTD(KLON,KLEV),         ZQD(KLON,KLEV),&
 & ZMFUS(KLON,KLEV),       ZMFDS(KLON,KLEV),&
 & ZMFUQ(KLON,KLEV),       ZMFDQ(KLON,KLEV),&
 & ZDMFUP(KLON,KLEV),      ZDMFDP(KLON,KLEV),&
 & ZMFUL(KLON,KLEV),       ZRFL(KLON),&
 & ZUU(KLON,KLEV),         ZVU(KLON,KLEV),&
 & ZUD(KLON,KLEV),         ZVD(KLON,KLEV),&
 & ZKINEU(KLON,KLEV),      ZKINED(KLON,KLEV) 
REAL(KIND=JPRB) ::     ZENTR(KLON),            ZHCBASE(KLON),&
 & ZMFUB(KLON),            ZMFUB1(KLON),&
 & ZDQCV(KLON)  
REAL(KIND=JPRB) ::                  ZDPMEL(KLON,KLEV),ZLGLAC(KLON,KLEV)
REAL(KIND=JPRB) ::     ZDHPBL(KLON),       ZWUBASE(KLON)

REAL(KIND=JPRB) ::     ZDMFEN(KLON,KLEV),      ZDMFDE(KLON,KLEV)

INTEGER(KIND=JPIM) ::  ILAB(KLON,KLEV),        IDTOP(KLON),ICTOP0(KLON),ILWMIN(KLON)
INTEGER(KIND=JPIM) ::  IDPL(KLON) ! departure level for convection
REAL(KIND=JPRB) ::     ZCAPE(KLON), ZHEAT(KLON)
LOGICAL ::  LLDDRAF(KLON), LLDDRAF3(KLON), LLDCUM(KLON)
LOGICAL ::  LLO1, LLO2(KLON)

INTEGER(KIND=JPIM) :: IKB, ITOPM2, JK, IK, JL

REAL(KIND=JPRB) ::   ZCONS2, ZCONS, ZDH,&
 & ZDQMIN, ZDZ, ZEPS, ZFAC, &
 & ZMFMAX, ZPBMPT, ZQUMQE, ZRO, ZMFA, ZERATE, ZDERATE, ZORCPD, ZRDOCPD,&
 & ZALV, ZSFL(KLON)
          
REAL(KIND=JPRB) :: ZTAU(KLON)  ! adjustment time

! scaling factor for momentum and tracer massflux
REAL(KIND=JPRB) :: ZMFS(KLON),  ZMFUUS(KLON,KLEV), ZMFDUS(KLON,KLEV) ,&
                  &             ZMFUDR(KLON,KLEV), ZMFDDR(KLON,KLEV) ,&
                  & ZMFUUB(KLON), ZMFUVB(KLON),&
                  & ZMF_SHAL(KLON), zz

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!*UPG change to operations
!    LOCALS FOR CONSERVATION CHECK
LOGICAL :: LLCONSCHECK=.false.
INTEGER(KIND=JPIM) :: JN
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZTENT, ZTENQ, ZTENU, ZTENV, ZSUMC
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE :: ZTENC

#include "cuascn.intfb.h"
#include "cubasen.intfb.h"
#include "cuddrafn.intfb.h"
#include "cudlfsn.intfb.h"
#include "cudtdqn.intfb.h"
#include "cududv.intfb.h"
#include "cuflxn.intfb.h"
#include "cuinin.intfb.h"
#include "cuctracer.intfb.h"

#include "fcttre.h"
!---------------------------------------------------------------------
!*UPG Change to operations call SATUR routine here

!     0.           Compute Saturation specific humidity
!                  ------------------------------------

    LDCUM(:)=.FALSE.
    PQSEN(:,:)=PQEN(:,:)

    CALL SATUR (KIDIA , KFDIA , KLON  , NJKT2, KLEV,&
               &PAP,    PTEN  , PQSEN , 1  )


!*UPG
!---------------------------------------------------------------------

!     1.           SPECIFY CONSTANTS AND PARAMETERS
!                  --------------------------------

IF (LHOOK) CALL DR_HOOK('CUMASTRN',0,ZHOOK_HANDLE)
ZCONS2=RMFCFL/(RG*PTSPHY)
ZCONS=1.0_JPRB/(RG*PTSPHY)
ZORCPD=1.0_JPRB/RCPD
ZRDOCPD=RD*ZORCPD

!----------------------------------------------------------------------

!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------

CALL CUININ &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,&
 & PVERVEL,  PGEO,     PAPH,     PAP,&
 & ILWMIN,   ILAB,&
 & ZTENH,    ZQENH,    ZQSENH,   PGEOH,&
 & PTU,      PQU,      ZTD,      ZQD,&
 & ZUU,      ZVU,      ZUD,      ZVD,&
 & PLU     )  

!---------------------------------------------------------------------

!*    3.0          CLOUD BASE CALCULATIONS
!                  -----------------------

!*             (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
!                  ---------------------------------------

CALL CUBASEN &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & LDLAND,&
 & ZTENH,    ZQENH,    PGEOH,    PAPH,&
 & PQHFL,    PAHFS,    PSSTRU,   PSSTRV,   PVERVEL,&
 & PTEN,     PQEN,     PGEO,&
 & PUEN,     PVEN,&
 & PTU,      PQU,      PLU,      ZUU,      ZVU,    ZWUBASE,&
 & ILAB,     LDCUM,    LDSC,     KCBOT,    KBOTSC,&
 & ICTOP0,   IDPL,     PCAPE )   


!*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*                 DECIDE ON TYPE OF CUMULUS CONVECTION 
!*                 ONE THE BASIS OF THE DEPTH OF THE CONVECTION
!*                 DEEP IF CLOUD DEPTH > 200MB
!*                 SHALLOW IF CLOUD DEPTH <200MB
!                  -----------------------------------------

! CALCULATE COLUMN AND SUB CLOUD LAYER MOISTURE CONVERGENCE
! AND SUB CLOUD LAYER MOIST STATIC ENERGY CONVERGENCE

DO JL=KIDIA,KFDIA
  ZDQCV(JL) =0.0_JPRB
  ZDHPBL(JL)=0.0_JPRB
  IDTOP(JL)=0
ENDDO
DO JK=NJKT2,KLEV
  DO JL=KIDIA,KFDIA
    ZDQCV(JL)=ZDQCV(JL)+MAX(0.0_JPRB,PTENQ(JL,JK))*(PAPH(JL,JK+1)-PAPH(JL,JK))
    IF(LDCUM(JL).AND.JK >= KCBOT(JL)) THEN
      ZDHPBL(JL)=ZDHPBL(JL)+(RLVTT*PTENQ(JL,JK)+RCPD*PTENT(JL,JK))&
       & *(PAPH(JL,JK+1)-PAPH(JL,JK))  
    ENDIF
  ENDDO
ENDDO

!*                 ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*                 CALCULATIONS IN CUASC AND INITIAL DETERMINATION OF 
!*                 CLOUD TYPE
!*                 (MAX.POSSIBLE CLOUD HEIGHT
!*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
!                  -------------------------------------------------

!DIR$ IVDEP
!OCL NOVREC      
DO JL=KIDIA,KFDIA
  IF (LDCUM(JL)) THEN
    IKB=KCBOT(JL)
    ZHCBASE(JL)=RCPD*PTU(JL,IKB)+PGEOH(JL,IKB)+RLVTT*PQU(JL,IKB)
    IF(.NOT.LDCUM(JL)) ICTOP0(JL)=-1
  ELSE
    ZHCBASE(JL)=0.0_JPRB
  ENDIF
ENDDO

!*                 SPECIFY INITIAL CLOUD TYPE
!*

DO JL=KIDIA,KFDIA
  IF (LDCUM(JL)) THEN
    IKB=KCBOT(JL)
    ITOPM2=ICTOP0(JL)
    ZPBMPT=PAPH(JL,IKB)-PAPH(JL,ITOPM2)
    IF (ZPBMPT >= RDEPTHS) THEN
      KTYPE(JL)=1
    ELSE
      KTYPE(JL)=2
    ENDIF
  ELSE
    KTYPE(JL)=0
  ENDIF
ENDDO

!*             (C) calculate initial updraught mass flux
!*                 and set lateral mixing rates
!*
!*                 for deep convection assume it is 10% of 
!*                 maximum value which is determined by the 
!*                 thickness of the layer and timestep
!*
!*                 for shallow convection calculated assuming
!*                 a balance of moist static energy in the 
!*                 sub-cloud layer (ignores present of downdraughts)
!                  ------------------------------------------

!DIR$ IVDEP
!OCL NOVREC
IF (LMFWSTAR) THEN
  DO JL=KIDIA,KFDIA
    IF (LDCUM(JL)) THEN
      IKB=KCBOT(JL)
      ZDZ=MAX(0.0_JPRB,MIN(1.5E3_JPRB,(PGEOH(JL,IKB)-PGEOH(JL,KLEV+1))/RG))
      ZMF_SHAL(JL)=0.07_JPRB*(RG/PTEN(JL,KLEV)*ZDZ*&
        &MAX(0.0_JPRB,-PAHFS(JL,KLEV+1)*ZORCPD-RETV*PTEN(JL,KLEV)*PQHFL(JL,KLEV+1)))**.3333
      ZMFMAX=(PAPH(JL,IKB)-PAPH(JL,IKB-1))*ZCONS2*RMFLIC+RMFLIA
      ZMF_SHAL(JL)=MIN(ZMF_SHAL(JL),ZMFMAX)
    ENDIF
  ENDDO
ENDIF

DO JL=KIDIA,KFDIA
  IF (LDCUM(JL)) THEN
    IKB=KCBOT(JL)
    ZMFMAX=(PAPH(JL,IKB)-PAPH(JL,IKB-1))*ZCONS2*RMFLIC+RMFLIA

! deep convection

    IF (KTYPE(JL) == 1) THEN
      ZMFUB(JL)=ZMFMAX*0.1_JPRB
      ZENTR(JL)=ENTRPEN

    ELSEIF (KTYPE(JL) == 2) THEN

! shallow convection

      ZQUMQE=PQU(JL,IKB)+PLU(JL,IKB)-ZQENH(JL,IKB)
      ZDQMIN=MAX(0.01_JPRB*ZQENH(JL,IKB),1.E-10_JPRB)
      ZDH=RCPD*(PTU(JL,IKB)-ZTENH(JL,IKB))+RLVTT*ZQUMQE
      ZDH=RG*MAX(ZDH,1.E5_JPRB*ZDQMIN)
      IF (ZDHPBL(JL) > 0.0_JPRB) THEN
        ZMFUB(JL)=ZDHPBL(JL)/ZDH
      !EPS: temporary solution for explicit
        if(ptsphy>1800._JPRB.and.rmfcfl==1.0_jprb) then
          ZMFUB(JL)=MIN(ZMFUB(JL),3._JPRB*ZMFMAX)
        else
          ZMFUB(JL)=MIN(ZMFUB(JL),ZMFMAX)
        endif
      ELSE
        ZMFUB(JL)=ZMFMAX*0.1_JPRB
        LDCUM(JL)=.FALSE.
      ENDIF
      IF(LMFWSTAR) ZMFUB(JL)=ZMF_SHAL(JL)
      ZENTR(JL)=ENTRSCV
    ENDIF

  ELSE

! no buoyancy cloud base from surface
! set cloud base mass flux and mixing rate
! to default value for safety

    ZMFUB(JL)=0.0_JPRB
    ZENTR(JL)=ENTRSCV
  ENDIF
ENDDO

!-----------------------------------------------------------------------

!*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!                  -------------------------------------------

!*             (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*                 CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
!*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
!                  -------------------------------------------------

! CALCULATIONS NOW DONE IS SECTION 3 ABOVE SO THAT
! INITIAL CLOUD DEPTH CAN BE USED TO SPECIFY
! THE TYPE OF CONVECTION

!*             (B) DO ASCENT IN 'CUASC'IN ABSENCE OF DOWNDRAFTS
!                  --------------------------------------------

CALL CUASCN &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & PTSPHY,&
 & ZTENH,    ZQENH,    PUEN,     PVEN,&
 & PTEN,     PQEN,     PQSEN,    PLITOT,&
 & PGEO,     PGEOH,    PAP,      PAPH,&
 & PTENQ,    PVERVEL,  zwubase,  ILWMIN,&
 & LDLAND,   LDCUM,    KTYPE,    ILAB,&
 & PTU,      PQU,      PLU,&
 & PMFU,     ZMFUB,    ZENTR,    ZLGLAC,&
 & ZMFUS,    ZMFUQ,    ZMFUL,    PLUDE,    ZDMFUP,&
 & ZDMFEN,&  
 & KCBOT,    KCTOP,    ICTOP0,   IDPL,     PMFUDE_RATE,   ZKINEU,   PWMEAN )   

!*         (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!              CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!              -----------------------------------------------------
!DIR$ IVDEP
!OCL NOVREC
DO JL=KIDIA,KFDIA
  IF (LDCUM(JL)) THEN
    IKB=KCBOT(JL)
    ITOPM2=KCTOP(JL)
    ZPBMPT=PAPH(JL,IKB)-PAPH(JL,ITOPM2)
    IF(KTYPE(JL) == 1.AND.ZPBMPT < RDEPTHS) KTYPE(JL)=2
    IF(KTYPE(JL) == 2.AND.ZPBMPT >= RDEPTHS) KTYPE(JL)=1
    ICTOP0(JL)=KCTOP(JL)
    IF(KTYPE(JL) == 1) ZENTR(JL)=ENTRPEN
    IF(KTYPE(JL) == 2) ZENTR(JL)=ENTRSCV
  ENDIF
  ZRFL(JL)=ZDMFUP(JL,1)
ENDDO
DO JK=2,KLEV
  DO JL=KIDIA,KFDIA
    ZRFL(JL)=ZRFL(JL)+ZDMFUP(JL,JK)
  ENDDO
ENDDO
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PMFD(JL,JK)=0.0_JPRB
    ZMFDS(JL,JK)=0.0_JPRB
    ZMFDQ(JL,JK)=0.0_JPRB
    ZDMFDP(JL,JK)=0.0_JPRB
    ZDPMEL(JL,JK)=0.0_JPRB
  ENDDO
ENDDO

!-----------------------------------------------------------------------

!*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
!                  ------------------------------

IF(LMFDD) THEN

!*             (A) DETERMINE LFS IN 'CUDLFS'
!                  -------------------------

  CALL CUDLFSN &
   & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
   & KCBOT,    KCTOP,    LDLAND,   LDCUM,&
   & ZTENH,    ZQENH,    PUEN,     PVEN,&
   & PTEN,     PQSEN,    PGEO,&
   & PGEOH,    PAPH,     PTU,      PQU,      PLU,&
   & ZUU,      ZVU,      ZMFUB,    ZRFL,&
   & ZTD,      ZQD,&
   & PMFD,     ZMFDS,    ZMFDQ,    ZDMFDP,&
   & IDTOP,    LLDDRAF )  

!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                  -----------------------------------------------

  CALL CUDDRAFN &
   & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
   & LLDDRAF,&
   & ZTENH,    ZQENH,    PUEN,     PVEN,&
   & PGEO,     PGEOH,    PAPH,     ZRFL,&
   & ZTD,      ZQD,      PMFU,&
   & PMFD,     ZMFDS,    ZMFDQ,    ZDMFDP,&
   & ZDMFDE,   PMFDDE_RATE,        ZKINED )  

ENDIF

!*                 (C)  RECALCULATE CLOUD BASE MASSFLUX FROM A
!*                 CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
!*                 AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
!*                 ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)          
!                  --------------------------------------------

!   DEEP CONVECTION

DO JL=KIDIA,KFDIA
  ZHEAT(JL)=0.0_JPRB
  ZCAPE(JL)=0.0_JPRB
  ZMFUB1(JL)=ZMFUB(JL)
ENDDO
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    LLO1=LDCUM(JL).AND.KTYPE(JL) == 1
    IF(LLO1.AND.JK <= KCBOT(JL).AND.JK > KCTOP(JL)) THEN
      IKB=KCBOT(JL)
      ZRO=PAPH(JL,JK)/(RD*ZTENH(JL,JK)*(1.0_JPRB+RETV*ZQENH(JL,JK)))
      ZDZ=(PGEOH(JL,JK-1)-PGEOH(JL,JK))
      ZHEAT(JL)=ZHEAT(JL) +&
       & (  (PTEN(JL,JK-1)-PTEN(JL,JK) + ZDZ*ZORCPD)/ZTENH(JL,JK)&
       & +  RETV*(PQEN(JL,JK-1)-PQEN(JL,JK))  ) *&
       & (RG*(PMFU(JL,JK)+PMFD(JL,JK)))/ZRO  
      ZCAPE(JL)=ZCAPE(JL) +&
       & ((PTU(JL,JK)-ZTENH(JL,JK))/ZTENH(JL,JK)&
       & +RETV*(PQU(JL,JK)-ZQENH(JL,JK))&
       & -PLU(JL,JK) ) * ZDZ  
    ENDIF
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  IF(LDCUM(JL).AND.KTYPE(JL) == 1) THEN
    IKB=KCBOT(JL)
    IK=KCTOP(JL)
    ZCAPE(JL)=MAX(0.0_JPRB,MIN(ZCAPE(JL),5000.0_JPRB))
    ZHEAT(JL)=MAX(1.E-4_JPRB,ZHEAT(JL))
    ZTAU(JL)=(PGEOH(JL,IK)-PGEOH(JL,IKB))/((2.0_JPRB+MIN(15.0_JPRB,PWMEAN(JL)))*RG)*RTAU
    ZTAU(JL)=MAX(PTSPHY,MIN(10800._JPRB,ZTAU(JL)))
    ZTAU(JL)=MAX(720.0_JPRB,ZTAU(JL))
    ZMFUB1(JL)=(ZCAPE(JL)*ZMFUB(JL))/(ZHEAT(JL)*ZTAU(JL))
    ZMFUB1(JL)=MAX(ZMFUB1(JL),0.001_JPRB)
    ZMFMAX=(PAPH(JL,IKB)-PAPH(JL,IKB-1))*ZCONS2*RMFLIC+RMFLIA
    ZMFUB1(JL)=MIN(ZMFUB1(JL),ZMFMAX)
  ENDIF
ENDDO

!  SHALLOW CONVECTION AND MID_LEVEL

!DIR$ IVDEP
!OCL NOVREC
DO JL=KIDIA,KFDIA
  IF ( LDCUM(JL) .AND. (KTYPE(JL) == 2.OR. KTYPE(JL) == 3) ) THEN
    IKB=KCBOT(JL)
    IF(PMFD(JL,IKB) < 0.0_JPRB) THEN
      ZEPS=-PMFD(JL,IKB)/MAX(ZMFUB(JL),1.E-10_JPRB)
    ELSE
      ZEPS=0.0_JPRB
    ENDIF
    ZQUMQE=PQU(JL,IKB)+PLU(JL,IKB)-&
     & ZEPS*ZQD(JL,IKB)-(1.0_JPRB-ZEPS)*ZQENH(JL,IKB)  
    ZDQMIN=MAX(0.01_JPRB*ZQENH(JL,IKB),1.E-10_JPRB)
! maximum permisable value of ud base mass flux 
    ZMFMAX=(PAPH(JL,IKB)-PAPH(JL,IKB-1))*ZCONS2*RMFLIC+RMFLIA

! shallow convection

    IF(KTYPE(JL) == 2) THEN
      ZDH=RCPD*(PTU(JL,IKB)-ZEPS*ZTD(JL,IKB)-&
       & (1.0_JPRB-ZEPS)*ZTENH(JL,IKB))+RLVTT*ZQUMQE  
      ZDH=RG*MAX(ZDH,1.E5_JPRB*ZDQMIN)
      IF(ZDHPBL(JL) > 0.0_JPRB) THEN
        ZMFUB1(JL)=ZDHPBL(JL)/ZDH
      ELSE
        ZMFUB1(JL)=ZMFUB(JL)
      ENDIF
   !EPS: temporary solution for explicit
      if(ptsphy>1800._JPRB.and.rmfcfl==1.0_jprb) then
        zmfub1(jl)=min(zmfub1(jl),3._JPRB*zmfmax)
      else
        zmfub1(jl)=min(zmfub1(jl),zmfmax)
      endif
      IF(LMFWSTAR) ZMFUB1(JL)=ZMF_SHAL(JL)
    ENDIF

! mid-level convection

    IF(KTYPE(JL) == 3)THEN
      ZMFUB1(JL)=ZMFUB(JL)*(1.0_JPRB+ZEPS)
      ZMFUB1(JL)=MIN(ZMFUB1(JL),ZMFMAX)
    ENDIF

  ENDIF
ENDDO

! rescale DD fluxes if deep and shallow convection

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF ( LLDDRAF(JL) .AND.( KTYPE(JL) == 1.OR. KTYPE(JL) == 2 ) ) THEN
      ZFAC=ZMFUB1(JL)/MAX(ZMFUB(JL),1.E-10_JPRB)
      PMFD(JL,JK)=PMFD(JL,JK)*ZFAC
      ZMFDS(JL,JK)=ZMFDS(JL,JK)*ZFAC
      ZMFDQ(JL,JK)=ZMFDQ(JL,JK)*ZFAC
      ZDMFDP(JL,JK)=ZDMFDP(JL,JK)*ZFAC
!  also rescale detrainment flux for ERA pp
      PMFDDE_RATE(JL,JK)=PMFDDE_RATE(JL,JK)*ZFAC
    ENDIF
  ENDDO
ENDDO

IF(LMFIT) THEN

  DO JK=2,KLEV-1
    DO JL=KIDIA,KFDIA
      ZUU(JL,JK)=PUEN(JL,JK-1)
      ZVU(JL,JK)=PVEN(JL,JK-1)
    ENDDO
  ENDDO

  ! reset updraught mass flux at cloud base

  DO JL=KIDIA,KFDIA
    ZMFUB(JL)=ZMFUB1(JL)
  ENDDO

!-----------------------------------------------------------------------

!*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
!*                 FOR PENETRATIVE CONVECTION (TYPE=1),
!*                 FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
!*                 AND FOR MID-LEVEL CONVECTION (TYPE=3).
!                  -------------------------------------------------

  CALL CUASCN &
   & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
   & PTSPHY,&
   & ZTENH,    ZQENH,    PUEN,     PVEN,&
   & PTEN,     PQEN,     PQSEN,    PLITOT,&
   & PGEO,     PGEOH,    PAP,      PAPH,&
   & PTENQ,    PVERVEL,  zwubase,  ILWMIN,&
   & LDLAND,   LDCUM,    KTYPE,    ILAB,&
   & PTU,      PQU,      PLU,&
   & PMFU,     ZMFUB,    ZENTR,    ZLGLAC,&
   & ZMFUS,    ZMFUQ,    ZMFUL,    PLUDE,    ZDMFUP,&
   & ZDMFEN,&
   & KCBOT,    KCTOP,    ICTOP0,   IDPL,     PMFUDE_RATE,    ZKINEU,   PWMEAN )   
  
  DO JL=KIDIA,KFDIA
    IF (LDCUM(JL)) THEN
      IKB=KCBOT(JL)
      ITOPM2=KCTOP(JL)
      ZPBMPT=PAPH(JL,IKB)-PAPH(JL,ITOPM2)
      IF(KTYPE(JL) == 1.AND.ZPBMPT < RDEPTHS) KTYPE(JL)=2
      IF(KTYPE(JL) == 2.AND.ZPBMPT >= RDEPTHS) KTYPE(JL)=1
    ENDIF
  ENDDO

ELSE

  DO JL=KIDIA,KFDIA
    IF(LDCUM(JL)) THEN
      ZMFS(JL)=ZMFUB1(JL)/MAX(RMFCMIN,ZMFUB(JL))
    ENDIF
  ENDDO
  DO JK=2,KLEV
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL).AND.JK>=KCTOP(JL)-1) THEN
        IKB=KCBOT(JL)
        IF(JK>IKB) THEN
          ZDZ=((PAPH(JL,KLEV+1)-PAPH(JL,JK))/(PAPH(JL,KLEV+1)-PAPH(JL,IKB)))
          PMFU(JL,JK)=PMFU(JL,IKB)*ZDZ
        ENDIF
        ZMFMAX=(PAPH(JL,JK)-PAPH(JL,JK-1))*ZCONS2*RMFLIC+RMFLIA
        IF(PMFU(JL,JK)*ZMFS(JL)>ZMFMAX) &
           & ZMFS(JL)=MIN(ZMFS(JL),ZMFMAX/PMFU(JL,JK))  
      ENDIF
    ENDDO
  ENDDO
  DO JK=2,KLEV
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL).AND.JK<=KCBOT(JL).AND.JK>=KCTOP(JL)-1) THEN
        PMFU(JL,JK)=PMFU(JL,JK)*ZMFS(JL)
        ZMFUS(JL,JK)=ZMFUS(JL,JK)*ZMFS(JL)
        ZMFUQ(JL,JK)=ZMFUQ(JL,JK)*ZMFS(JL)
        ZMFUL(JL,JK)=ZMFUL(JL,JK)*ZMFS(JL)
        ZDMFUP(JL,JK)=ZDMFUP(JL,JK)*ZMFS(JL)
        ZDMFEN(JL,JK)=ZDMFEN(JL,JK)*ZMFS(JL)
        PLUDE(JL,JK)=PLUDE(JL,JK)*ZMFS(JL)
        PMFUDE_RATE(JL,JK)=PMFUDE_RATE(JL,JK)*ZMFS(JL)
      ENDIF
    ENDDO
  ENDDO

ENDIF

!-----------------------------------------------------------------------

!*    6.5          IN CASE THAT EITHER DEEP OR SHALLOW IS SWITCHED OFF
!                  RESET LDCUM TO FALSE-> FLUXES SET TO ZERO IN CUFLXN
!                  ---------------------------------------------------

!                 exclude pathological KTYPE=2 KCBOT=KCTOP=KLEV-1

DO JL=KIDIA,KFDIA
  IF(KTYPE(JL)==2.AND.KCBOT(JL)==KCTOP(JL).AND.KCBOT(JL)>=KLEV-1) THEN
    LDCUM(JL)=.FALSE.
    KTYPE(JL)=0
  ENDIF
ENDDO

IF (.NOT.LMFSCV .OR. .NOT.LMFPEN) THEN
  DO JL=KIDIA,KFDIA
    LLO2(JL)=.FALSE.
    IF((.NOT.LMFSCV .AND. KTYPE(JL)==2).OR.(.NOT.LMFPEN .AND. KTYPE(JL)==1))THEN
      LLO2(JL)=.TRUE.
      LDCUM(JL)=.FALSE.
    ENDIF
  ENDDO
ENDIF

!-----------------------------------------------------------------------

!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------

!- set DD mass fluxes to zero above cloud top
!  (because of inconsistency with second updraught)
DO JL=KIDIA,KFDIA
  IF(LLDDRAF(JL).AND.IDTOP(JL)<=KCTOP(JL)) THEN
    IDTOP(JL)=KCTOP(JL)+1
  ENDIF
ENDDO

DO JK=2,KLEV
  DO JL=KIDIA,KFDIA
    IF (LLDDRAF(JL)) THEN
      IF (JK<IDTOP(JL)) THEN
        PMFD(JL,JK)=0.0_JPRB
        ZMFDS(JL,JK)=0.0_JPRB
        ZMFDQ(JL,JK)=0.0_JPRB
        PMFDDE_RATE(JL,JK)=0.0_JPRB
        ZDMFDP(JL,JK)=0.0_JPRB
      ELSEIF (JK==IDTOP(JL)) THEN
        PMFDDE_RATE(JL,JK)=0.0_JPRB
      ENDIF
    ENDIF
  ENDDO
ENDDO
CALL CUFLXN &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & PTSPHY,&
 & PTEN,     PQEN,     PQSEN,    ZTENH,    ZQENH,&
 & PAPH,     PAP,      PGEOH,    LDLAND,   LDCUM,&
 & KCBOT,    KCTOP,    IDTOP,    ITOPM2,&
 & KTYPE,    LLDDRAF,&
 & PMFU,     PMFD,     ZMFUS,    ZMFDS,&
 & ZMFUQ,    ZMFDQ,    ZMFUL,    PLUDE,&
 & ZDMFUP,   ZDMFDP,   ZDPMEL,   ZLGLAC,&
 & PMFLXR,   PMFLXS,   PRAIN,    PMFDDE_RATE )  
 
!- rescale DD fluxes if total mass flux becomes negative
!- correct DD detrainment rates if entrainment becomes negative
!- correct UD detrainment rates if entrainment becomes negative
!- conservation correction for precip

ZMFS(:)=1.0_JPRB
!DO JK=2,KLEV-1
DO JK=2,KLEV ! change for stability
  DO JL=KIDIA,KFDIA
    IF ( LLDDRAF(JL) .AND. JK>=IDTOP(JL)-1 ) THEN
       ZMFMAX=PMFU(JL,JK)*0.98_JPRB
       IF(PMFD(JL,JK)+ZMFMAX+1.E-15_JPRB<0.0_JPRB) THEN
          ZMFS(JL)=MIN(ZMFS(JL),-ZMFMAX/PMFD(JL,JK))
       ENDIF
    ENDIF
  ENDDO
ENDDO
ZMFUUB(:)=0.0_JPRB
DO JK=2,KLEV
  DO JL=KIDIA,KFDIA
    IF ( ZMFS(JL)<1.0_JPRB .AND. JK>=IDTOP(JL)-1 ) THEN
      PMFD(JL,JK)=PMFD(JL,JK)*ZMFS(JL)
      ZMFDS(JL,JK)=ZMFDS(JL,JK)*ZMFS(JL)
      ZMFDQ(JL,JK)=ZMFDQ(JL,JK)*ZMFS(JL)
      PMFDDE_RATE(JL,JK)=PMFDDE_RATE(JL,JK)*ZMFS(JL)
      ZMFUUB(JL)=ZMFUUB(JL)-(1.0_JPRB-ZMFS(JL))*ZDMFDP(JL,JK)
      PMFLXR(JL,JK+1)=PMFLXR(JL,JK+1)+ZMFUUB(JL)
      ZDMFDP(JL,JK)=ZDMFDP(JL,JK)*ZMFS(JL)
    ENDIF
  ENDDO
ENDDO

DO JK=2,KLEV-1
  DO JL=KIDIA,KFDIA
    IF ( LLDDRAF(JL) .AND. JK>=IDTOP(JL)-1 ) THEN
       ZERATE=-PMFD(JL,JK)+PMFD(JL,JK-1)+PMFDDE_RATE(JL,JK)
       IF(ZERATE<0.0_JPRB) THEN
         PMFDDE_RATE(JL,JK)=PMFDDE_RATE(JL,JK)-ZERATE
       ENDIF
    ENDIF
    IF ( LDCUM(JL) .AND. JK>=KCTOP(JL)-1 ) THEN
       ZERATE=PMFU(JL,JK)-PMFU(JL,JK+1)+PMFUDE_RATE(JL,JK)
       IF(ZERATE<0.0_JPRB) THEN
         PMFUDE_RATE(JL,JK)=PMFUDE_RATE(JL,JK)-ZERATE
       ENDIF
     ! ZDMFUP(JL,JK)=ZDMFUP(JL,JK)+ZDMFDP(JL,JK)
       ZDMFUP(JL,JK)=PMFLXR(JL,JK+1)+PMFLXS(JL,JK+1)&
                   &-PMFLXR(JL,JK)-PMFLXS(JL,JK)
       ZDMFDP(JL,JK)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

! avoid negative humidities at ddraught top
DO JL=KIDIA,KFDIA
 IF(LLDDRAF(JL)) THEN
   JK=IDTOP(JL)
   IK=MIN(JK+1,KLEV)
   IF(ZMFDQ(JL,JK)<0.3_JPRB*ZMFDQ(JL,IK)) THEN
     IF(RMFSOLTQ==0.0_JPRB) THEN
       ZMFDQ(JL,JK)=0.3_JPRB*ZMFDQ(JL,IK)
     ELSE
       PMFD(JL,JK)=0.3_JPRB*PMFD(JL,IK)
     ENDIF
   ENDIF
 ENDIF
ENDDO

! avoid negative humidities near cloud top because gradient of precip flux
! and detrainment / liquid water flux too large
DO JK=2,KLEV
  DO JL=KIDIA,KFDIA
    IF(LDCUM(JL).AND.JK>=KCTOP(JL)-1.AND.JK<KCBOT(JL)) THEN
      ZDZ=PTSPHY*RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
      ZMFA=ZMFUQ(JL,JK+1)+ZMFDQ(JL,JK+1)-ZMFUQ(JL,JK)-ZMFDQ(JL,JK)+&
      &ZMFUL(JL,JK+1)-ZMFUL(JL,JK)+ZDMFUP(JL,JK)
      ZMFA=(ZMFA-PLUDE(JL,JK))*ZDZ
      IF(PQEN(JL,JK)+ZMFA<0.0_JPRB) THEN
        PLUDE(JL,JK)=PLUDE(JL,JK)+2.0_JPRB*(PQEN(JL,JK)+ZMFA)/ZDZ
      ENDIF
      IF(PLUDE(JL,JK)<0.0_JPRB) THEN
        PLUDE(JL,JK)=0.0_JPRB
      ENDIF
    ENDIF
  ENDDO
ENDDO

!*UPG change to operations
IF ( LLCONSCHECK ) THEN
    ALLOCATE(ZTENT(KLON,KLEV))
    ALLOCATE(ZTENQ(KLON,KLEV))
    ALLOCATE(ZTENU(KLON,KLEV))
    ALLOCATE(ZTENV(KLON,KLEV))
    DO JK=2,KLEV
    DO JL=KIDIA,KFDIA
       IF ( LDCUM(JL) ) THEN
            ZTENT(JL,JK)=PTENT(JL,JK)
            ZTENQ(JL,JK)=PTENQ(JL,JK)
            ZTENU(JL,JK)=PTENU(JL,JK)
            ZTENV(JL,JK)=PTENV(JL,JK)
       ENDIF
    ENDDO
    ENDDO
    IF ( LMFTRAC .AND. KTRAC>0 ) THEN
       ALLOCATE(ZTENC(KLON,KLEV,KTRAC))
       ALLOCATE(ZSUMC(KLON,4+KTRAC))
       DO JN=1,KTRAC
          DO JK=2,KLEV
          DO JL=KIDIA,KFDIA
             IF ( LDCUM(JL) ) THEN
                ZTENC(JL,JK,JN)=PTENC(JL,JK,JN)
             ENDIF
          ENDDO
          ENDDO
       ENDDO
    ELSE
       ALLOCATE(ZSUMC(KLON,4))
    ENDIF
ENDIF
!*UPG change to operations

!----------------------------------------------------------------------

!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------

IF( RMFSOLTQ>0.0_JPRB) THEN
! derive draught properties for implicit

  DO JK=KLEV,2,-1
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
        IF(JK>KCBOT(JL)) THEN
          ZMFA=1.0_JPRB/MAX(1.E-15_JPRB,PMFU(JL,JK))
          PQU(JL,JK)=ZQENH(JL,JK)+ZMFUQ(JL,JK)*ZMFA
          PTU(JL,JK)=ZTENH(JL,JK)+ZMFUS(JL,JK)*ZMFA*ZORCPD
          ZMFUS(JL,JK)=PMFU(JL,JK)*(RCPD*PTU(JL,JK)+PGEOH(JL,JK))
          ZMFUQ(JL,JK)=PMFU(JL,JK)*PQU(JL,JK)
          IF(LLDDRAF(JL)) THEN
            ZMFA=1.0_JPRB/MIN(-1.E-15_JPRB,PMFD(JL,JK))
            ZQD(JL,JK)=ZQENH(JL,JK)+ZMFDQ(JL,JK)*ZMFA
            ZTD(JL,JK)=ZTENH(JL,JK)+ZMFDS(JL,JK)*ZMFA*ZORCPD
            ZMFDQ(JL,JK)=PMFD(JL,JK)*ZQD(JL,JK)
            ZMFDS(JL,JK)=PMFD(JL,JK)*(RCPD*ZTD(JL,JK)+PGEOH(JL,JK))
          ENDIF
        ELSEIF(JK<=KCBOT(JL).AND.JK>=KCTOP(JL)) THEN
          ZMFUS(JL,JK)=PMFU(JL,JK)*(RCPD*PTU(JL,JK)+PGEOH(JL,JK))
          ZMFUQ(JL,JK)=PMFU(JL,JK)*PQU(JL,JK)
          ZMFDS(JL,JK)=PMFD(JL,JK)*(RCPD*ZTD(JL,JK)+PGEOH(JL,JK))
          ZMFDQ(JL,JK)=PMFD(JL,JK)*ZQD(JL,JK)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  
ENDIF

CALL CUDTDQN &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & ITOPM2,   KTYPE,    KCTOP,    IDTOP,    LDCUM,    LLDDRAF,   PTSPHY,&
 & PAPH,     PGEOH,    PGEO,&
 & PTEN,     ZTENH,    PQEN,     ZQENH,    PQSEN,&
 & ZLGLAC,   PLUDE,    PMFU,     PMFD,&
 & ZMFUS,    ZMFDS,    ZMFUQ,    ZMFDQ,&
 & ZMFUL,    ZDMFUP,   ZDPMEL,&
 & PTENT,    PTENQ,    PENTH )

!----------------------------------------------------------------------

!*    9.0          COMPUTE MOMENTUM IN UPDRAUGHT AND DOWNDRAUGHT
!                  ---------------------------------------------

IF(LMFDUDV) THEN

  DO JK=KLEV-1,2,-1
    IK=JK+1
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
        IF(JK==KCBOT(JL).AND.KTYPE(JL)<3) THEN
           IKB=IDPL(JL)
           ZUU(JL,JK)=PUEN(JL,IKB-1)
           ZVU(JL,JK)=PVEN(JL,IKB-1)
        ELSEIF(JK==KCBOT(JL).AND.KTYPE(JL)==3) THEN
           ZUU(JL,JK)=PUEN(JL,JK-1)
           ZVU(JL,JK)=PVEN(JL,JK-1)
        ENDIF
        IF( JK<KCBOT(JL).AND.JK>=KCTOP(JL)) THEN
          ZFAC=0.0_JPRB
          IF(KTYPE(JL)==1.OR.KTYPE(JL)==3) ZFAC=2.0_JPRB
          if(ktype(jl)==1.and.jk<=kctop(jl)+2) zfac=3.0_jprb
          ZERATE=PMFU(JL,JK)-PMFU(JL,IK)+(1.0_JPRB+ZFAC)*PMFUDE_RATE(JL,JK)
          ZDERATE=(1.0_JPRB+ZFAC)*PMFUDE_RATE(JL,JK)
          ZMFA=1.0_JPRB/MAX(RMFCMIN,PMFU(JL,JK))
          ZUU(JL,JK)=(ZUU(JL,IK)*PMFU(JL,IK)+ZERATE*PUEN(JL,JK)-ZDERATE*ZUU(JL,IK))*ZMFA
          ZVU(JL,JK)=(ZVU(JL,IK)*PMFU(JL,IK)+ZERATE*PVEN(JL,JK)-ZDERATE*ZVU(JL,IK))*ZMFA
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DO JK=3,KLEV
    IK=JK-1
    DO JL=KIDIA,KFDIA
      IF( LDCUM(JL)) THEN
        IF(JK==IDTOP(JL)) THEN
          ZUD(JL,JK)=0.5_JPRB*(ZUU(JL,JK)+PUEN(JL,IK))
          ZVD(JL,JK)=0.5_JPRB*(ZVU(JL,JK)+PVEN(JL,IK))
        ELSEIF(JK>IDTOP(JL)) THEN
          ZERATE=-PMFD(JL,JK)+PMFD(JL,IK)+PMFDDE_RATE(JL,JK)
          ZMFA=1.0_JPRB/MIN(-RMFCMIN,PMFD(JL,JK))
          ZUD(JL,JK)=(ZUD(JL,IK)*PMFD(JL,IK)-ZERATE*PUEN(JL,IK)+PMFDDE_RATE(JL,JK)*ZUD(JL,IK))*ZMFA
          ZVD(JL,JK)=(ZVD(JL,IK)*PMFD(JL,IK)-ZERATE*PVEN(JL,IK)+PMFDDE_RATE(JL,JK)*ZVD(JL,IK))*ZMFA
        ENDIF
      ENDIF
    ENDDO
  ENDDO

!*    9.1          UPDATE TENDENCIES FOR U AND V IN SUBROUTINE CUDUDV
!                  --------------------------------------------------

! for explicit/semi-implicit rescale massfluxes for stability in Momentum
!------------------------------------------------------------------------

  ZMFS(:)=1.0_JPRB
! IF(RMFSOLUV<=0.5_JPRB) THEN
  IF(RMFSOLUV<=1.0_JPRB) THEN
    DO JK=2,KLEV
      DO JL=KIDIA,KFDIA
        IF(LDCUM(JL).AND.JK>=KCTOP(JL)-1) THEN
          ZMFMAX=(PAPH(JL,JK)-PAPH(JL,JK-1))*ZCONS
          IF(PMFU(JL,JK)>ZMFMAX.AND.JK>=KCTOP(JL)) &
           & ZMFS(JL)=MIN(ZMFS(JL),ZMFMAX/PMFU(JL,JK))  
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZMFUUS(JL,JK)=PMFU(JL,JK)
      ZMFDUS(JL,JK)=PMFD(JL,JK)
      IF(LDCUM(JL).AND.JK>=KCTOP(JL)-1) THEN
        ZMFUUS(JL,JK)=PMFU(JL,JK)*ZMFS(JL)
        ZMFDUS(JL,JK)=PMFD(JL,JK)*ZMFS(JL)
      ENDIF
    ENDDO
  ENDDO

! recompute Draught properties below for Implicit
! based on linear flux profiles

  IF(RMFSOLUV>0.0_JPRB) THEN

    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
        JK=KCBOT(JL)
        IK=JK-1
        ZMFUUB(JL)=ZMFUUS(JL,JK)*(ZUU(JL,JK)-PUEN(JL,IK))
        ZMFUVB(JL)=ZMFUUS(JL,JK)*(ZVU(JL,JK)-PVEN(JL,IK))
      ENDIF
    ENDDO

    DO JK=2,KLEV
      IK=JK-1
      DO JL=KIDIA,KFDIA
        IF ( LDCUM(JL).AND.JK>KCBOT(JL) ) THEN
          IKB=KCBOT(JL)
          ZDZ=((PAPH(JL,KLEV+1)-PAPH(JL,JK))/(PAPH(JL,KLEV+1)-PAPH(JL,IKB)))
          IF(KTYPE(JL) == 3) THEN
            ZDZ=ZDZ*ZDZ
          ENDIF
          ZMFA=1.0_JPRB/MAX(RMFCMIN,ZMFUUS(JL,JK))
          ZUU(JL,JK)=PUEN(JL,IK)+ZMFUUB(JL)*ZDZ*ZMFA
          ZVU(JL,JK)=PVEN(JL,IK)+ZMFUVB(JL)*ZDZ*ZMFA

          ZMFDUS(JL,JK)=ZMFDUS(JL,IKB)*ZDZ
          ZUD(JL,JK)=PUEN(JL,IK)+ZUD(JL,IKB)-PUEN(JL,IKB-1)
          ZVD(JL,JK)=PVEN(JL,IK)+ZVD(JL,IKB)-PVEN(JL,IKB-1)
        ENDIF
    ! add UV perturb to correct wind bias
        IF ( LDCUM(JL).AND.JK>=KCTOP(JL) ) THEN
          ZUU(JL,JK)=ZUU(JL,JK)-RUVPER*SIGN(1.0_JPRB,ZUU(JL,JK))
          ZVU(JL,JK)=ZVU(JL,JK)-RUVPER*SIGN(1.0_JPRB,ZVU(JL,JK))
        ENDIF
      ENDDO
    ENDDO

  ENDIF

!-------------------------------------------------------------------
! End
! Intermediate Solution for stability in EPS: 
! For original code replace line
!  &, PUEN,     PVEN,     ZMFUUS,   ZMFDUS &
!by
!  &, PUEN,     PVEN,     PMFU,     PMFD

  CALL CUDUDV &
   & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
   & ITOPM2,   KTYPE,    KCBOT,    KCTOP,    LDCUM,    PTSPHY,&
   & PAPH,     PUEN,     PVEN,     ZMFUUS,   ZMFDUS,&
   & ZUU,      ZUD,      ZVU,      ZVD,&
   & PTENU,    PTENV     )  

ENDIF

!----------------------------------------------------------------------

!*   10.           IN CASE THAT EITHER DEEP OR SHALLOW IS SWITCHED OFF
!                  NEED TO SET SOME VARIABLES A POSTERIORI TO ZERO
!                  ---------------------------------------------------

IF (.NOT.LMFSCV .OR. .NOT.LMFPEN) THEN
  DO JK=2,KLEV
    DO JL=KIDIA,KFDIA
      IF(LLO2(JL).AND.JK>=KCTOP(JL)-1) THEN
        PTU(JL,JK)  =PTEN(JL,JK)
        PQU(JL,JK)  =PQEN(JL,JK)
        PLU(JL,JK)  =0.0_JPRB
        PENTH(JL,JK) =0.0_JPRB
        PMFUDE_RATE(JL,JK) =0.0_JPRB
        PMFDDE_RATE(JL,JK) =0.0_JPRB
      ENDIF
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    IF(LLO2(JL)) THEN
      KCTOP(JL)=KLEV-1
      KCBOT(JL)=KLEV-1
    ENDIF
  ENDDO
ENDIF

!----------------------------------------------------------------------

!*   11.0          CHEMICAL TRACER TRANSPORT
!                  -------------------------

IF ( LMFTRAC .AND. KTRAC>0 ) THEN

! transport switched off for mid-level convection
  DO JL=KIDIA,KFDIA
    !IF( LDCUM(JL).AND.KTYPE(JL)/=3 ) THEN
     IF( LDCUM(JL).AND.KTYPE(JL)/=3.and.kcbot(jl)-kctop(jl)>=1 ) THEN
         LLDCUM(JL)=.TRUE.
         LLDDRAF3(JL)=LLDDRAF(JL)
     ELSE
         LLDCUM(JL)=.FALSE.
         LLDDRAF3(JL)=.FALSE.
     ENDIF
  ENDDO
  
! check and correct mass fluxes for CFL criterium

  ZMFS(:)=1.0_JPRB
  IF(RMFSOLCT<=3.0_JPRB) THEN
    DO JK=2,KLEV
      DO JL=KIDIA,KFDIA
        IF(LLDCUM(JL).AND.JK>=KCTOP(JL)) THEN
          ZMFMAX=(PAPH(JL,JK)-PAPH(JL,JK-1))*1.0_JPRB*ZCONS
          IF(PMFU(JL,JK)>ZMFMAX) &
           & ZMFS(JL)=MIN(ZMFS(JL),ZMFMAX/PMFU(JL,JK))  
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      IF(LLDCUM(JL).AND.JK>=KCTOP(JL)-1) THEN
        ZMFUUS(JL,JK)=PMFU(JL,JK)*ZMFS(JL)
        ZMFUDR(JL,JK)=PMFUDE_RATE(JL,JK)*ZMFS(JL)
      ELSE
        ZMFUUS(JL,JK)=0._JPRB
        ZMFUDR(JL,JK)=0._JPRB
      ENDIF
      IF ( LLDDRAF3(JL) .AND. JK>=IDTOP(JL)-1) THEN
        ZMFDUS(JL,JK)=PMFD(JL,JK)*ZMFS(JL)
        ZMFDDR(JL,JK)=PMFDDE_RATE(JL,JK)*ZMFS(JL)
      ELSE
        ZMFDUS(JL,JK)=0._JPRB
        ZMFDDR(JL,JK)=0._JPRB
      ENDIF
    ENDDO
  ENDDO

  IF( LMFSMOOTH ) THEN
! smmoothing of mass fluxes (gradients) at top and bottom of draughts
    DO JK=2,KLEV-1
      DO JL=KIDIA,KFDIA
        IF(LLDDRAF3(JL).AND.ZMFDUS(JL,JK)<0.0_JPRB .AND. ZMFDUS(JL,JK+1)==0.0_JPRB) THEN
          ZERATE=MIN(0._JPRB,ZMFDUS(JL,JK)-0.5*ZMFDUS(JL,JK-1))
          ZMFDUS(JL,JK)=ZMFDUS(JL,JK)-ZERATE
          ZMFDDR(JL,JK)=ZMFDDR(JL,JK)-ZERATE
          ZMFDDR(JL,JK+1)=-ZMFDUS(JL,JK)
        ENDIF
        IF(LLDCUM(JL).AND.JK==KCTOP(JL)) THEN
          ZERATE=MAX(0.0_JPRB,ZMFUUS(JL,JK)-0.5_JPRB*ZMFUUS(JL,JK+1))
          ZMFUUS(JL,JK)=ZMFUUS(JL,JK)-ZERATE
          ZMFUDR(JL,JK)=ZMFUDR(JL,JK)+ZERATE
          ZMFUDR(JL,JK-1)=ZMFUUS(JL,JK)
        ENDIF
      ENDDO
    ENDDO
    DO JK=KLEV-1,2,-1
      DO JL=KIDIA,KFDIA
        IF(LLDCUM(JL)) THEN
          IF(ZMFUDR(JL,JK)==0.0_JPRB.AND.ZMFUDR(JL,JK-1)>0.0_JPRB) THEN
            ZMFUDR(JL,JK)=0.5_JPRB*ZMFUDR(JL,JK-1)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  CALL CUCTRACER &
   & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,     KTRAC,&
   & KTYPE,    KCTOP,    KCBOT,    IDPL,     IDTOP,&
   & LLDCUM,   LLDDRAF3,  PTSPHY,   PAPH,&
   & ZMFUUS,   ZMFDUS,   ZMFUDR,   ZMFDDR,&
   & PCEN,     PTENC     )  

ENDIF

!----------------------------------------------------------------------

!*   12.           PUT DETRAINMENT RATES FROM MFLX UNITS IN UNITS MFLX/M 
!                  FOR ERA40
!                  ---------------------------------------------------

DO JK=2,KLEV
  DO JL=KIDIA,KFDIA
    IF ( LDCUM(JL) ) THEN
      ZRO=RG/(PGEOH(JL,JK)-PGEOH(JL,JK+1))  ! 1/dz
      PMFUDE_RATE(JL,JK)=PMFUDE_RATE(JL,JK)*ZRO
      PMFDDE_RATE(JL,JK)=PMFDDE_RATE(JL,JK)*ZRO
      IF(JK<KCTOP(JL)) THEN
         PLU(JL,JK)=0.0_JPRB
         PTU(JL,JK)=PTEN(JL,JK)
         PQU(JL,JK)=PQEN(JL,JK)
      ENDIF
    ENDIF
  ENDDO
ENDDO

!----------------------------------------------------------------------
!*UPG change to operations

IF ( LLCONSCHECK ) THEN

!*   13.0          CONSERVATION CHECK and CORRECTION
!                  ---------------------------------

    DO JL=KIDIA,KFDIA
      ZSUMC(JL,:)=0.
    ENDDO
    DO JK=KLEV,2,-1
    DO JL=KIDIA,KFDIA
     IF ( LDCUM(JL) .AND. JK>=KCTOP(JL)-1) THEN
       ZDZ=(PAPH(JL,JK+1)-PAPH(JL,JK))/RG
       ZSUMC(JL,1)=ZSUMC(JL,1)+(PTENQ(JL,JK)-ZTENQ(JL,JK))*ZDZ+PLUDE(JL,JK)
       ZALV=FOELHMCU(PTEN(JL,JK))
       ZSUMC(JL,2)=ZSUMC(JL,2)+RCPD*(PTENT(JL,JK)-ZTENT(JL,JK))*ZDZ-ZALV*PLUDE(JL,JK)
       ZSUMC(JL,3)=ZSUMC(JL,3)+(PTENU(JL,JK)-ZTENU(JL,JK))*ZDZ
       ZSUMC(JL,4)=ZSUMC(JL,4)+(PTENV(JL,JK)-ZTENV(JL,JK))*ZDZ
     ENDIF
    ENDDO
    ENDDO
    IF ( LMFTRAC .AND. KTRAC>0 ) THEN
      DO JN=1,KTRAC
         DO JK=KLEV,2,-1
         DO JL=KIDIA,KFDIA
            IF ( LDCUM(JL) .AND. JK>=KCTOP(JL)-1) THEN
               ZDZ=(PAPH(JL,JK+1)-PAPH(JL,JK))/RG
               ZSUMC(JL,4+JN)=ZSUMC(JL,4+JN)+(PTENC(JL,JK,JN)-ZTENC(JL,JK,JN))*ZDZ
            ENDIF
          ENDDO
          ENDDO
       ENDDO
     ENDIF

    DO JL=KIDIA,KFDIA
     IF ( LDCUM(JL) ) THEN
       ZALV=FOELHMCU(PTEN(JL,KLEV))
       ZSFL(JL)=PMFLXR(JL,KLEV+1)+PMFLXS(JL,KLEV+1)

       write(61,'(i4,a9,2f15.8,i4,a9,f15.8,a10,2f15.8)')jl,' CONS q: ',&
          &-zsumc(jl,1)*zalv,zsfl(jl)*zalv,ktype(jl),&
          &' CONS h: ',zsumc(jl,2),' CONS uv: ',zsumc(jl,3),zsumc(jl,4)
       if ( lmftrac .and. ktrac>0 ) then
          write(61,*)' Conserv Error Tracers 1-',ktrac,' :'
         do jn=1,ktrac
          write(61,'(i4,e12.4)')jn,zsumc(jl,4+jn)
         enddo
       endif

       IKB=KCTOP(JL)
       ZDZ=(PAPH(JL,KLEV+1)-PAPH(JL,IKB-1))/RG
       ZSUMC(JL,1)=(ZSUMC(JL,1)+ZSFL(JL))/ZDZ
       ZSUMC(JL,2)=(ZSUMC(JL,2)-ZALV*ZSFL(JL))/(ZDZ*RCPD)
     ENDIF
    END DO

  ! DO JK=KLEV,2,-1
  ! DO JL=KIDIA,KFDIA
  !  IKB=KCTOP(JL)
  !  IF ( LDCUM(JL) .AND. JK >= IKB-1) THEN
  !    PTENQ(JL,JK)=PTENQ(JL,JK)-ZSUMC(JL,1)
  !    PTENT(JL,JK)=PTENT(JL,JK)-ZSUMC(JL,2)
  !    PTENU(JL,JK)=PTENU(JL,JK)-ZSUMC(JL,3)
  !    PTENV(JL,JK)=PTENV(JL,JK)-ZSUMC(JL,4)
  !  ENDIF
  ! ENDDO
  ! ENDDO

    DEALLOCATE(ZSUMC)
    IF ( LMFTRAC .AND. KTRAC>0 ) THEN
       DEALLOCATE(ZTENC)
    ENDIF
    DEALLOCATE(ZTENV)
    DEALLOCATE(ZTENU)
    DEALLOCATE(ZTENQ)
    DEALLOCATE(ZTENT)

ENDIF
!----------------------------------------------------------------------

!*    14.0         COMPUTE CONVECTIVE TENDENCIES FOR LIQUID AND SOLID
!                  CLOUD CONDENSATE, CHANGE PRECIP UNITS IN M/S
!                  --------------------------------------------------

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PTENL(JL,JK)=PLUDE(JL,JK)*RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
    PTENI(JL,JK)=(1.0_JPRB-FOEALFA(PTEN(JL,JK)))*PTENL(JL,JK)
    PTENL(JL,JK)=PTENL(JL,JK)-PTENI(JL,JK)
    PMFLXR(JL,JK)=PMFLXR(JL,JK)*1.E-3
    PMFLXS(JL,JK)=PMFLXS(JL,JK)*1.E-3
  ENDDO
ENDDO
  DO JL=KIDIA,KFDIA
    PMFLXR(JL,KLEV+1)=PMFLXR(JL,KLEV+1)*1.E-3
    PMFLXS(JL,KLEV+1)=PMFLXS(JL,KLEV+1)*1.E-3
  ENDDO
!----------------------------------------------------------------------
!*UPG Change to operations


IF (LHOOK) CALL DR_HOOK('CUMASTRN',1,ZHOOK_HANDLE)
END SUBROUTINE CUMASTRN
