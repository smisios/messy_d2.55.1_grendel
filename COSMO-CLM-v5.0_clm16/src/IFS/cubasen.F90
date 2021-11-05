SUBROUTINE CUBASEN &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & LDLAND,&
 & PTENH,    PQENH,    PGEOH,    PAPH,&
 & PQHFL,    PAHFS,    PSSTRU,   PSSTRV,   PWN,&
 & PTEN,     PQEN,     PGEO,&
 & PUEN,     PVEN,&
 & PTU,      PQU,      PLU,      PUU,      PVU ,  PWUBASE,&
 & KLAB,     LDCUM,    LDSC,     KCBOT,    KBOTSC,&
 & KCTOP,    KDPL,     PCAPE )  

!          THIS ROUTINE CALCULATES CLOUD BASE FIELDS
!          CLOUD BASE HEIGHT AND CLOUD TOP HEIGHT

!          A. Pier Siebesma   KNMI ********      
!          modified C Jakob (ECMWF) (01/2001) 
!          modified P Bechtold (ECMWF) (08/2002) 
!          (include cycling over levels to find unstable departure/base level+
!           mixed layer properties +w Trigger)

!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD BASE AND CLOUD TOP VALUES FOR CU-PARAMETRIZATION

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD FIELDS VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=0 FOR STABLE LAYERS
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CLOUD LEVELS LEVEL

!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD TOP
!          (ENTRAINING PLUME, WITH ENTRAINMENT PROPORTIONAL TO (1/Z))

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS

!    INPUT PARAMETERS (REAL):

! not used at the moment because we want to use linear intepolation
! for fields on the half levels.

!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG

!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    *PSSTRU*       KINEMATIC surface U-MOMENTUM FLUX             (M/S)^2
!    *PSSTRV*       KINEMATIC surface V-MOMENTUM FLUX             (M/S)^2
!    *PWN*          NORMALIZED LARGE-SCALE VERTICAL VELOCITY      (M/S)
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2

!    UPDATED PARAMETERS (REAL):

!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PUU*          U-VELOCITY IN UPDRAFTS                         M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                         M/S

!    UPDATED PARAMETERS (INTEGER):

!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CLOUD LEVELS

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDSC*         FLAG: .TRUE. IF BL-CLOUDS EXIST

!    OUTPUT PARAMETERS (INTEGER):

!    *KCBOT*       CLOUD BASE LEVEL !    
!    *KCTOP*       CLOUD TOP LEVEL = HEIGHEST HALF LEVEL 
!                  WITH A NON-ZERO CLOUD UPDRAFT.
!    *KBOTSC*      CLOUD BASE LEVEL OF BL-CLOUDS
!    *KDPL*        DEPARTURE LEVEL
!    *PCAPE*       PSEUDOADIABATIQUE max CAPE (J/KG)

!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT

!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!             02-11-02 : Use fixed last possible departure level and 
!                        last updraft computation level for bit-reproducibility
!                                            D.Salmond &  J. Hague
!             03-07-03 : Tuning for p690     J. Hague
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK_IFS   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RCPD     ,RETV, RD, RG,&
 & RLVTT    ,RLSTT    ,RTT  
USE YOEVDF   , ONLY : RKAP   
USE YOECUMF  , ONLY : LMFDUDV, ENTRPEN, RDEPTHS, NJKT1, NJKT2
USE YOECLDP  , ONLY : RLMIN
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 & RALVDCP  ,RALSDCP  ,RALFDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
 & RTWAT_RTICECU_R    ,RTWAT_RTICE_R  

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM)               :: KTDIA ! Argument NOT used
LOGICAL           ,INTENT(IN)    :: LDLAND(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQENH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQHFL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAHFS(KLON,KLEV+1) 
REAL(KIND=JPRB)                  :: PSSTRU(KLON) ! Argument NOT used
REAL(KIND=JPRB)                  :: PSSTRV(KLON) ! Argument NOT used
REAL(KIND=JPRB)                  :: PWN(KLON,KLEV) ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWUBASE(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KLAB(KLON,KLEV) 
LOGICAL           ,INTENT(INOUT) :: LDCUM(KLON) 
LOGICAL           ,INTENT(OUT)   :: LDSC(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCBOT(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KBOTSC(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCTOP(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KDPL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPE(KLON) 
INTEGER(KIND=JPIM) ::  ICTOP(KLON),            ICBOT(KLON),&
 & IBOTSC(KLON),           ILAB(KLON,KLEV),&
 & IDPL(KLON)  

!             LOCAL STORAGE
!             ----- -------

LOGICAL ::         LL_LDBASE(KLON),&
 & LLGO_ON(KLON),&
 & LLDEEP(KLON),    LLDCUM(KLON), &
 & LLDSC(KLON),     LLFIRST(KLON)  
LOGICAL ::     LLRESET,        LLRESETJL(KLON)

INTEGER(KIND=JPIM) :: ICALL, IK, IKB, IS, JK, JL, JKK, JKT1, JKT2, JKT, JKB

REAL(KIND=JPRB)    :: ZS(KLON,KLEV),&
 & ZSENH(KLON,KLEV+1),&
 & ZQENH(KLON,KLEV+1),&
 & ZSUH (KLON,KLEV),&
 & ZWU2H(KLON,KLEV),&
 & ZBUOH(KLON,KLEV)  
REAL(KIND=JPRB) :: ZQOLD(KLON),ZPH(KLON)
REAL(KIND=JPRB) :: ZMIX(KLON)
REAL(KIND=JPRB) :: ZDZ(KLON),zcbase(klon)

REAL(KIND=JPRB) ::    ZLU(KLON,KLEV),   ZQU(KLON,KLEV),&
 & ZTU(KLON,KLEV), &
 & ZUU(KLON,KLEV),   ZVU(KLON,KLEV)  

REAL(KIND=JPRB) :: ZCAPE(KLON,KLEV) ! local for CAPE at every departure level

REAL(KIND=JPRB) :: ZBUOF, ZZ, ZC2, ZEPSADD
REAL(KIND=JPRB) :: ZRHO      ! DENSITY AT SURFACE (KG/M^3) 
REAL(KIND=JPRB) :: ZKHVFL    ! SURFACE BUOYANCY FLUX (K M/S)
REAL(KIND=JPRB) :: ZWS       ! SIGMA_W AT LOWEST MODEL HALFLEVEL (M/S)
REAL(KIND=JPRB) :: ZQEXC     ! HUMIDITY EXCESS AT LOWEST MODEL HALFLEVEL (KG/KG)
REAL(KIND=JPRB) :: ZTEXC     ! TEMPERATURE EXCESS AT LOWEST MODEL HALFLEVEL (K)
REAL(KIND=JPRB) :: ZEPS      ! FRACTIONAL ENTRAINMENT RATE   [M^-1]
REAL(KIND=JPRB) :: ZTVENH    ! ENVIRONMENT VIRTUAL TEMPERATURE AT HALF LEVELS (K)  
REAL(KIND=JPRB) :: ZTVUH     ! UPDRAFT VIRTUAL TEMPERATURE AT HALF LEVELS     (K)
REAL(KIND=JPRB) :: ZLGLAC    ! UPDRAFT LIQUID WATER FROZEN IN ONE LAYER
REAL(KIND=JPRB) :: zqsu, ZCOR, zdq, zalfaw, zfacw, zfaci, zfac,&
 & zesdp, zdqsdt, zdtdp, zdp,zpdifftop,zpdiffbot,ZSF,ZQF,zaw,zbw  
REAL(KIND=JPRB) :: ZTVEN1, ZTVEN2, ZTVU1, ZTVU2 ! pseudoadiabatique T_v
REAL(KIND=JPRB) :: ZDTVTRIG(KLON) ! virtual temperatures
REAL(KIND=JPRB) :: ZWORK1, ZWORK2 ! work arrays for T and w perturbations
REAL(KIND=JPRB) :: ZRCPD, ZRG, ZTMP
!REAL(KIND=JPRB) :: ZTINY
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "cuadjtq.intfb.h"
#include "fcttre.h"

!----------------------------------------------------------------------
!     0.           INITIALIZE CONSTANTS AND FIELDS
!                  -------------------------------
!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUBASEN',0,ZHOOK_HANDLE)
ZC2    = 0.55_JPRB
ZAW    = 1.0_JPRB
ZBW    = 1.0_JPRB
ZEPSADD= 1.E-4_JPRB

DO JL=KIDIA,KFDIA
  PWUBASE(JL)=0.0_JPRB
  LLGO_ON(JL)=.TRUE.
  LLFIRST(JL)=.TRUE.
  KDPL(JL)=KLEV
ENDDO

! Set last possible departure level and last updraft computation level
! NOT Bit-reproducible
!DO JK=KLEV+1,2,-1
!  IF((PAPH(KIDIA,KLEV+1)-PAPH(KIDIA,JK)) < 350.E2_JPRB ) JKT1=JK
!  IF(PAPH(KIDIA,JK) > 60.E2_JPRB ) JKT2=JK
!END DO

JKT1=NJKT1
JKT2=NJKT2
ZRG=1.0_JPRB/RG
ZRCPD=1.0_JPRB/RCPD
! The function TINY is standard Fortran90, but is not available on all systems
!ZTINY=TINY(ZRG)
!ZTINY=1.0E-35

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZTU(JL,JK) = PTU(JL,JK)
    ZQU(JL,JK) = PQU(JL,JK)
    ZLU(JL,JK) = PLU(JL,JK)
    ZUU(JL,JK) = PUU(JL,JK)
    ZVU(JL,JK) = PVU(JL,JK)
    ILAB(JL,JK)= KLAB(JL,JK)
    ZCAPE(JL,JK)= 0.0_JPRB
  ENDDO
ENDDO

!----------------------------------------------------------------------
!       -----------------------------------------------------------
!       1.1  PREPARE FIELDS ON HALF LEVELS BY LINEAR INTERPOLATION
!             OF SPECIFIC HUMIDITY AND STATIC ENERGY
!       -----------------------------------------------------------

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZWU2H(JL,JK)=0.0_JPRB
    ZS   (JL,JK) = RCPD*PTEN(JL,JK) + PGEO(JL,JK)
    ZQENH(JL,JK) = PQENH(JL,JK)
    ZSENH(JL,JK) = RCPD*PTENH(JL,JK)+PGEOH(JL,JK)
  ENDDO
ENDDO

DO JKK=KLEV,JKT1,-1 ! Big external loop for level testing:
                    ! find first departure level that produces deepest cloud top
                    ! or take surface level for shallow convection and Sc
   !
   !        ---------------------------------------------------------
   !        1.2    INITIALISE FIELDS AT DEPARTURE HALF MODEL LEVEL
   !        ---------------------------------------------------------
   !
  IS=0
  DO JL=KIDIA,KFDIA
    IF (LLGO_ON(JL)) THEN
      IS=IS+1
      IDPL(JL)    =JKK      ! departure level
      ICBOT  (JL) =JKK      ! cloud base level for convection, (-1 if not found)
      IBOTSC (JL) =KLEV-1   ! sc    base level for sc-clouds , (-1 if not found)
      ICTOP(JL)   =KLEV-1   ! cloud top for convection (-1 if not found)
      LLDCUM(JL)  =.FALSE.  ! on exit: true if cloudbase=found
      LLDSC (JL)  =.FALSE.  ! on exit: true if cloudbase=found
      LL_LDBASE(JL)   =.FALSE. ! on exit: true if cloudbase=found
      ZDTVTRIG(JL) =0.0_JPRB
      ZUU(JL,JKK) =PUEN(JL,JKK)*(PAPH(JL,JKK+1)-PAPH(JL,JKK))
      ZVU(JL,JKK) =PVEN(JL,JKK)*(PAPH(JL,JKK+1)-PAPH(JL,JKK))
    ENDIF 
  ENDDO

  IF(IS /= 0) THEN

    IF(JKK == KLEV) THEN

      DO JL=KIDIA,KFDIA
        IF (LLGO_ON(JL)) THEN
          ZRHO  = PAPH(JL,JKK+1)/(RD*(PTEN(JL,JKK)*(1.+RETV*PQEN(JL,JKK))))
          ZKHVFL= (PAHFS(JL,JKK+1)/RCPD+RETV*PTEN(JL,JKK)*PQHFL(JL,JKK+1))/ZRHO
!          ZUST  = MAX(SQRT(PSSTRU(JL)**2 + PSSTRV(JL)**2),ZREPUST)     !u* (repust=10e-4)
!          ZWS=ZUST**3._JPRB- 1.5_JPRB*RKAP*ZKHVFL*PGEOH(JL,KLEV)/PTEN(JL,KLEV)
!!**$$     ZWS=0.001_JPRB - 1.5_JPRB*RKAP*ZKHVFL*PGEOH(JL,KLEV)/PTEN(JL,KLEV)
           ZWS=0.001_JPRB - 1.5_JPRB*RKAP*ZKHVFL &              !!**$$
             & *(PGEOH(JL,KLEV)-PGEOH(JL,KLEV+1))/PTEN(JL,KLEV) !!**$$
          IF( ZKHVFL < 0.0_JPRB ) THEN
            ZWS=1.2_JPRB*ZWS**.3333_JPRB
            ILAB(JL,JKK)= 1
            ZTEXC     = MAX(-1.5_JPRB*PAHFS(JL,JKK+1)/(ZRHO*ZWS*RCPD),0.0_JPRB)
            ZQEXC     = MAX(-1.5_JPRB*PQHFL(JL,JKK+1)/(ZRHO*ZWS),0.0_JPRB)
            ZQU (JL,JKK) = ZQENH(JL,JKK) + ZQEXC
            ZSUH (JL,JKK) = ZSENH(JL,JKK) + RCPD*ZTEXC
            ZTU (JL,JKK) = (ZSENH(JL,JKK)-PGEOH(JL,JKK))/RCPD + ZTEXC
            ZLU (JL,JKK) = 0.0_JPRB
            ZWU2H(JL,JKK) = ZWS**2
            ZLU (JL,JKK) = 0.0_JPRB
        !
        !  determine buoyancy at lowest half level
        !
            ZTVENH            = (1.0_JPRB+RETV*ZQENH(JL,JKK)) &
             & *(ZSENH(JL,JKK)-PGEOH(JL,JKK))/RCPD  
            ZTVUH             = (1.0_JPRB+RETV*ZQU(JL,JKK))*ZTU(JL,JKK)
            ZBUOH(JL,JKK) = (ZTVUH-ZTVENH)*RG/ZTVENH
          ELSE
            LLGO_ON(JL)=.FALSE.      ! non-convective point
          ENDIF
        ENDIF
      ENDDO
   
    ELSE

      DO JL=KIDIA,KFDIA
        IF (LLGO_ON(JL)) THEN
          ZRHO  = PAPH(JL,JKK+1)/(RD*(PTEN(JL,JKK)*(1.+RETV*PQEN(JL,JKK))))
          ILAB(JL,JKK)= 1
          ZTEXC=.2_JPRB
          ZQEXC=1.E-4_JPRB
          ZQU (JL,JKK) = ZQENH(JL,JKK) + ZQEXC
          ZSUH (JL,JKK) = ZSENH(JL,JKK) + RCPD*ZTEXC
          ZTU (JL,JKK) = (ZSENH(JL,JKK)-PGEOH(JL,JKK))*ZRCPD + ZTEXC
          ZLU (JL,JKK) = 0.0_JPRB
         ! construct mixed layer for parcels emanating in lowest 60 hPa
          IF (PAPH(JL,KLEV+1)-PAPH(JL,JKK-1)<60.E2_JPRB) THEN
            ZQU(JL,JKK) =0.0_JPRB
            ZSUH(JL,JKK)=0.0_JPRB
            ZWORK1      =0.0_JPRB
            DO JK=JKK+1,JKK-1,-1
              IF( ZWORK1 < 50.E2_JPRB ) THEN
                ZWORK2=PAPH(JL,JK)-PAPH(JL,JK-1)
                ZWORK1      =ZWORK1+ZWORK2
                ZQU(JL,JKK) =ZQU(JL,JKK) +ZQENH(JL,JK)*ZWORK2
                ZSUH(JL,JKK)=ZSUH(JL,JKK)+ZSENH(JL,JK)*ZWORK2
              ENDIF
            ENDDO
            ZQU(JL,JKK) =ZQU(JL,JKK) /ZWORK1+ZQEXC
            ZSUH(JL,JKK)=ZSUH(JL,JKK)/ZWORK1+RCPD*ZTEXC
            ZTU(JL,JKK) =(ZSUH(JL,JKK)-PGEOH(JL,JKK))/RCPD+ZTEXC
          ENDIF
          ZWU2H(JL,JKK) = 1.0_JPRB
      !
      !  determine buoyancy at lowest half level
      !
          ZTVENH            = (1.0_JPRB+RETV*ZQENH(JL,JKK)) &
           & *(ZSENH(JL,JKK)-PGEOH(JL,JKK))/RCPD  
          ZTVUH             = (1.0_JPRB+RETV*ZQU(JL,JKK))*ZTU(JL,JKK)
          ZBUOH(JL,JKK) = (ZTVUH-ZTVENH)*RG/ZTVENH
        ENDIF
      ENDDO
   
    ENDIF

  ENDIF
   
   !----------------------------------------------------------------------
   !     2.0          DO ASCENT IN SUBCLOUD AND LAYER,
   !                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
   !                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
   !                  CHECK FOR BUOYANCY AND SET FLAGS
   !                  -------------------------------------
   !       ------------------------------------------------------------
   !        1.2  DO THE VERTICAL ASCENT UNTIL VELOCITY BECOMES NEGATIVE
   !       ------------------------------------------------------------
  DO JK=JKK-1,JKT2,-1
    IS=0

    IF(JKK==KLEV) THEN ! 1/z mixing for shallow

      DO JL=KIDIA,KFDIA
        IF (LLGO_ON(JL)) THEN
          IS         = IS+1
          ZDZ(JL)        = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
          ZEPS       = ZC2/(PGEO(JL,JK)*ZRG + ZDZ(jl)) + ZEPSADD
          ZMIX(JL)       = 0.5_JPRB*ZDZ(JL)*ZEPS
          ZQF = (PQENH(JL,JK+1) + PQENH(JL,JK))*0.5_JPRB
          ZSF = (ZSENH(JL,JK+1) + ZSENH(JL,JK))*0.5_JPRB
          ZTMP = 1.0_JPRB/(1.0_JPRB+ZMIX(JL))
          ZQU(JL,JK)= (ZQU(JL,JK+1)*(1.0_JPRB-ZMIX(JL))&
           & +2.0_JPRB*ZMIX(jl)*ZQF) * ZTMP  
          ZSUH (JL,JK)= (ZSUH(JL,JK+1)*(1.0_JPRB-ZMIX(JL))&
           & +2.0_JPRB*ZMIX(jl)*ZSF) * ZTMP  
          ZQOLD(JL)  = ZQU(JL,JK)
          ZTU (JL,JK) = (ZSUH(JL,JK)-PGEOH(JL,JK))*ZRCPD
          ZPH  (JL)    = PAPH(JL,JK)
        ENDIF
      ENDDO
     
    ELSE

      DO JL=KIDIA,KFDIA
        IF (LLGO_ON(JL)) THEN
          IS         = IS+1
          ZDZ(JL)        = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
          ZQF = (PQENH(JL,JK+1) + PQENH(JL,JK))*0.5_JPRB
          ZSF = (ZSENH(JL,JK+1) + ZSENH(JL,JK))*0.5_JPRB
         !ZMIX(JL)=ENTRPEN*ZDZ(JL)
          ZMIX(JL)=2.0_JPRB*ENTRPEN*ZDZ(JL)*(PAPH(JL,JK)/PAPH(JL,KLEV+1))**3
          ZQU(JL,JK)= ZQU(JL,JK+1)*(1.0_JPRB-ZMIX(JL))+ ZQF*ZMIX(JL)
          ZSUH(JL,JK)= ZSUH(JL,JK+1)*(1.0_JPRB-ZMIX(JL))+ ZSF*ZMIX(JL)
          ZQOLD(JL)  = ZQU(JL,JK)
          ZTU (JL,JK) = (ZSUH(JL,JK)-PGEOH(JL,JK))*ZRCPD
          ZPH  (JL)    = PAPH(JL,JK)
        ENDIF
      ENDDO
  
    ENDIF
       
    IF (IS == 0) EXIT
     
    IK=JK
    ICALL=1
     
    CALL CUADJTQ &
     & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
     & IK,&
     & ZPH,      ZTU,      ZQU,      LLGO_ON,   ICALL)  
   
   !DIR$ IVDEP
   !OCL NOVREC
   
    DO JL=KIDIA,KFDIA
      IF(LLGO_ON(JL)) THEN
   
   ! add condensation to water
   
        ZDQ=MAX(ZQOLD(JL)-ZQU(JL,JK),0.0_JPRB)
        ZLU(JL,JK)=ZLU(JL,JK+1)+ZDQ

   ! freezing
   
        ZLGLAC=ZDQ*((1.0_JPRB-FOEALFCU(ZTU(JL,JK)))-&
         & (1.0_JPRB-FOEALFCU(ZTU(JL,JK+1))))  
              
   
   ! pseudo-microphysics
   
        if(jkk==klev) then  ! no precip for shallow
          ZLU(JL,JK)=MIN(ZLU(JL,JK),5.E-3_JPRB)
   !* chose a more pseudo-adiabatic formulation as original overestimates
   !* water loading efect and therefore strongly underestimates cloud thickness
        else 
          ZLU(JL,JK)=0.5_JPRB*ZLU(JL,JK) 
        endif
   
   ! update dry static energy after condensation + freezing
   
        ZSUH(JL,JK)    = RCPD*(ZTU(JL,JK)+RALFDCP*ZLGLAC)+PGEOH(JL,JK)
         
   ! Buoyancy on half and full levels
            
        ZTVUH           = (1.0_JPRB+RETV*ZQU(JL,JK)-ZLU(JL,JK))*ZTU(JL,JK)&
         & +RALFDCP*ZLGLAC  
        ZTVENH          = (1.0_JPRB+RETV*ZQENH(JL,JK)) &
         & *(ZSENH(JL,JK)-PGEOH(JL,JK))*ZRCPD  
        ZBUOH(JL,JK)   = (ZTVUH-ZTVENH)*RG/ZTVENH
        ZBUOF          = (ZBUOH(JL,JK) + ZBUOH(JL,JK+1))*0.5_JPRB
   
   ! solve kinetic energy equation
   
        ZTMP=1.0_JPRB/(1.0_JPRB+2.0_JPRB*ZBW*ZMIX(jl))
        ZWU2H(JL,JK) = (ZWU2H(JL,JK+1)*(1.0_JPRB-2.0_JPRB*ZBW*ZMIX(jl))&
         & +2.0_JPRB*ZAW*ZBUOF*ZDZ(jl)) * ZTMP  
   
   ! compute pseudoadiabatique CAPE for diagnostics
   
        ZTVU2 = ZTU(JL,JK)  *(1.0_JPRB+RETV*ZQU(JL,JK))
        ZTVEN2= PTENH(JL,JK)*(1.0_JPRB+RETV*PQENH(JL,JK))
        IF (JK == JKK-1) THEN
          ZTVU1  = ZTVU2
          ZTVEN1 = ZTVEN2
        ENDIF
        ZBUOF = (ZTVU2+ZTVU1-ZTVEN1-ZTVEN2)/ZTVEN2
        ZBUOF = ZBUOF*ZDZ(JL)*RG
        ZCAPE(JL,JKK)  = ZCAPE(JL,JKK) + MAX(0.0_JPRB,ZBUOF)
        ZTVU1=ZTVU2
        ZTVEN1=ZTVEN2
   
   ! first layer with liquid water - find exact cloud base
   
        IF(ZLU(JL,JK) >0.0_JPRB.AND.ILAB(JL,JK+1)==1) THEN
           
          IK=JK+1
          ZQSU=FOEEWM(ZTU(JL,IK))/PAPH(JL,IK)
          ZQSU=MIN(0.5_JPRB,ZQSU)
          ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSU)
          ZQSU=ZQSU*ZCOR
          ZDQ=MIN(0._JPRB,ZQU(JL,IK)-ZQSU)
          ZALFAW=FOEALFA(ZTU(JL,IK))
          ZFACW=R5LES/((ZTU(JL,IK)-R4LES)**2)
          ZFACI=R5IES/((ZTU(JL,IK)-R4IES)**2)
          ZFAC=ZALFAW*ZFACW+(1.-ZALFAW)*ZFACI
          ZESDP=FOEEWM(ZTU(JL,IK))/PAPH(JL,IK)
          ZCOR=1.0_JPRB/(1.0_JPRB-RETV*ZESDP)
          ZDQSDT=ZFAC*ZCOR*ZQSU
          ZDTDP=RD*ZTU(JL,IK)/(RCPD*PAPH(JL,IK))
          ZDP=ZDQ/(ZDQSDT*ZDTDP)
          ZCBASE(JL)=PAPH(JL,IK)+ZDP
           
   ! chose nearest half level as cloud base
   
          ZPDIFFTOP=ZCBASE(JL)-PAPH(JL,JK)
          ZPDIFFBOT=PAPH(JL,JK+1)-ZCBASE(JL)
           
          IF(ZPDIFFTOP > ZPDIFFBOT.AND.ZWU2H(JL,JK+1)>0.0_JPRB) THEN
            JKB=MIN(KLEV-1,JK+1)
            ILAB(JL,JKB)=2 
            ILAB(JL,JK)=2
            LL_LDBASE(JL) =.TRUE.
            LLDSC(JL)   =.TRUE.
            IBOTSC(JL) =JKB
            ICBOT(JL)  =JKB
            ZLU(JL,JK+1) = RLMIN
          ELSEIF(ZPDIFFTOP <= ZPDIFFBOT.AND.ZWU2H(JL,JK)>0.0_JPRB) THEN
            ILAB(JL,JK)=2
            LL_LDBASE(JL) =.TRUE.
            LLDSC(JL)   =.TRUE.
            IBOTSC(JL) =JK
            ICBOT(JL)  =JK
          ENDIF
          JKB=ICBOT(JL)
   
        ENDIF
   
   ! decide on presence of convection, cloud base and cloud top based on
   ! kinetic energy
   
        IF (ZWU2H(JL,JK) < 0.0_JPRB) THEN
          LLGO_ON(JL) = .FALSE.             
          IF (ZLU(JL,JK+1)>0.0_JPRB) THEN
            ICTOP(JL)   = JK
            LLDCUM(JL)   = .TRUE.
          ELSE
            LLDCUM(JL)   = .FALSE.
          ENDIF
        ELSE
          IF (ZLU(JL,JK)>0.0_JPRB) then
            ILAB(JL,JK) = 2
          ELSE
            ILAB(JL,JK) = 1
          ENDIF
        ENDIF
      ENDIF
    ENDDO
   
    IF(LMFDUDV.AND.JKK==KLEV) THEN
      DO JL=KIDIA,KFDIA
        IF(.NOT.LL_LDBASE(JL).AND.LLGO_ON(JL)) THEN
          ZUU(JL,JKK)=ZUU(JL,JKK)+PUEN(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
          ZVU(JL,JKK)=ZVU(JL,JKK)+PVEN(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
        ENDIF
      ENDDO
    ENDIF
   
!     IF (IS == 0) EXIT
  ENDDO
   
  IF( JKK==KLEV) THEN
      ! set values for departure level for PBL clouds = first model level
    DO JL=KIDIA,KFDIA
      LDSC(JL)  = LLDSC(JL)
      IF(LDSC(JL)) THEN
        KBOTSC(JL)= IBOTSC(JL)
      ELSE
        KBOTSC(JL)=-1
      ENDIF
    
      LLGO_ON(JL) = .FALSE.
      JKT=ICTOP(JL)
      JKB=ICBOT(JL)
      LLDEEP(JL)=PAPH(JL,JKB)-PAPH(JL,JKT)>RDEPTHS
      IF(LLDEEP(JL)) LLDCUM(JL)=.FALSE. ! no deep allowed for KLEV
      lldeep(jl)=.false. ! for deep convection start only at level KLEV-1
                            ! and form mixed layer, so go on
      ! test further for deep convective columns as not yet found
      IF ( LLDEEP(JL) ) LLFIRST(JL)=.FALSE.
      LLGO_ON(JL) = .NOT.LLDEEP(JL)
      IF(LLDCUM(JL)) THEN
        KCBOT(JL)= ICBOT(JL)
        KCTOP(JL)= ICTOP(JL)
        KDPL(JL)  = IDPL(JL)
        LDCUM(JL) = LLDCUM(JL)
        PWUBASE(JL)=SQRT(MAX(ZWU2H(JL,JKB),0.0_JPRB))
      ELSE
        KCTOP(JL)=-1
        KCBOT(JL)=-1
        KDPL(JL) =KLEV-1
        LDCUM(JL)=.FALSE.
        PWUBASE(JL)=0.0_JPRB
      ENDIF
    ENDDO
    DO JK=KLEV,1,-1
      DO JL=KIDIA,KFDIA
        JKT=ICTOP(JL)
        IF ( JK>=JKT ) THEN
          KLAB(JL,JK)=ILAB(JL,JK)
          PTU(JL,JK)=ZTU(JL,JK)
          PQU(JL,JK)=ZQU(JL,JK)
          PLU(JL,JK)=ZLU(JL,JK)
        ENDIF
      ENDDO
    ENDDO
  ! IF(LMFDUDV) THEN
  !   DO JL=KIDIA,KFDIA
  !     IF(LDCUM(JL)) THEN
  !       IKB=KCBOT(JL)
  !       ZZ=1.0_JPRB/(PAPH(JL,JKK+1)-PAPH(JL,IKB))
  !       PUU(JL,JKK)=ZUU(JL,JKK)*ZZ
  !       PVU(JL,JKK)=ZVU(JL,JKK)*ZZ
  !     ENDIF
  !   ENDDO
  ! ENDIF
  ENDIF
   
  IF( JKK < KLEV ) THEN
    LLRESET=.FALSE.
    DO JL=KIDIA,KFDIA
      IF ( .NOT.LLDEEP(JL) ) THEN
        JKT=ICTOP(JL)
        JKB=ICBOT(JL)
           ! test on cloud thickness and buoyancy
        LLDEEP(JL)=PAPH(JL,JKB)-PAPH(JL,JKT)>=RDEPTHS 
       !LLDEEP(JL)=PAPH(JL,JKB)-PAPH(JL,JKT)>=RDEPTHS &
       !   &.AND. ZDTVTRIG(JL)>0._JPRB
      ENDIF
      LLRESETJL(JL)=LLDEEP(JL).AND.LLFIRST(JL)
      LLRESET=LLRESET.OR.LLRESETJL(JL)
    ENDDO


    IF(LLRESET) THEN
      DO JK=KLEV,1,-1
        DO JL=KIDIA,KFDIA
         ! keep first departure level that produces deep cloud
!          IF ( LLDEEP(JL) .AND. LLFIRST(JL) ) THEN 
          IF ( LLRESETJL(JL) ) THEN 
            JKT=ICTOP(JL)
            JKB=IDPL(JL)
            IF ( JK<=JKB .AND. JK>=JKT ) THEN
              KLAB(JL,JK)=ILAB(JL,JK)
              PTU(JL,JK)=ZTU(JL,JK)
              PQU(JL,JK)=ZQU(JL,JK)
              PLU(JL,JK)=ZLU(JL,JK)
            ELSE 
              KLAB(JL,JK)=1
              PTU(JL,JK)=PTENH(JL,JK)
              PQU(JL,JK)=PQENH(JL,JK)
              PLU(JL,JK)=0.0_JPRB
            ENDIF
            IF ( JK<JKT ) KLAB(JL,JK)=0
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    DO JL=KIDIA,KFDIA
      IF ( LLDEEP(JL) .AND. LLFIRST(JL) ) THEN
        KDPL(JL)  = IDPL(JL)
        KCTOP(JL) = ICTOP(JL)
        KCBOT(JL) = ICBOT(JL)
        LDCUM(JL) = LLDCUM(JL)
        LDSC(JL)  = .FALSE.
        KBOTSC(JL)= -1
        JKB=KCBOT(JL)
        PWUBASE(JL)=SQRT(MAX(ZWU2H(JL,JKB),0.0_JPRB))
!  no initialization of wind for deep here, this is done in
!  CUINI and CUASCN
        LLFIRST(JL)=.FALSE.
      ENDIF
      LLGO_ON(JL) = .NOT.LLDEEP(JL)
    ENDDO
  ENDIF

ENDDO ! end of big loop for search of departure level     

      ! chose maximum CAPE value
DO JL=KIDIA,KFDIA
  PCAPE(JL) = MAXVAL(ZCAPE(JL,:))
ENDDO

IF (LHOOK) CALL DR_HOOK('CUBASEN',1,ZHOOK_HANDLE)
END SUBROUTINE CUBASEN
