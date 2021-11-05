SUBROUTINE SWR &
 &( KIDIA , KFDIA , KBDIM , KLEV  , KNU &
 &, PALBD , PCG   , PCLD , POMEGA, PSEC , PTAU &
 &, PCGAZ , PPIZAZ, PRAY1, PRAY2 , PREFZ, PRJ  , PRK , PRMUE &
 &, PTAUAZ, PTRA1 , PTRA2, PTRCLD &
 &)

!**** *SWR* - CONTINUUM SCATTERING COMPUTATIONS

!     PURPOSE.
!     --------
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
!     CONTINUUM SCATTERING

!**   INTERFACE.
!     ----------

!          *SWR* IS CALLED EITHER FROM *SW1S*
!                              OR FROM *SWNI*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES CONTINUUM FLUXES CORRESPONDING TO AEROSOL
!     OR/AND RAYLEIGH SCATTERING (NO MOLECULAR GAS ABSORPTION)

!     EXTERNALS.
!     ----------

!          *SWDE*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!        Ph. DANDIN Meteo-France 05-96 : Effect of cloud layer
!        JJMorcrette 990128 : sunshine duration
    
!     ------------------------------------------------------------------


USE MO_KIND  , ONLY : DP

USE MO_SW    , ONLY : NOVLP    ,NSW


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER :: KFDIA
INTEGER :: KIDIA
INTEGER :: KLEV
INTEGER :: KBDIM
INTEGER :: KNU



!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL(DP):: PALBD(KBDIM,NSW)      , PCG(KBDIM,NSW,KLEV)&
  &,  PCLD(KBDIM,KLEV)&
  &,  POMEGA(KBDIM,NSW,KLEV)&
  &,  PSEC(KBDIM)           , PTAU(KBDIM,NSW,KLEV)

REAL(DP):: PRAY1(KBDIM,KLEV+1)   , PRAY2(KBDIM,KLEV+1)&
  &,  PREFZ(KBDIM,2,KLEV+1) , PRJ(KBDIM,6,KLEV+1)&
  &,  PRK(KBDIM,6,KLEV+1)   , PRMUE(KBDIM,KLEV+1)&
  &,  PCGAZ(KBDIM,KLEV)     , PPIZAZ(KBDIM,KLEV)&
  &,  PTAUAZ(KBDIM,KLEV)&
  &,  PTRA1(KBDIM,KLEV+1)   , PTRA2(KBDIM,KLEV+1)&
  &,  PTRCLD(KBDIM)

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL(DP):: ZC1I(KBDIM,KLEV+1)    , ZCLEQ(KBDIM,KLEV)&
  &,  ZCLEAR(KBDIM)         , ZCLOUD(KBDIM) &
  &,  ZGG(KBDIM)            , ZREF(KBDIM)&
  &,  ZRE1(KBDIM)           , ZRE2(KBDIM)&
  &,  ZRMUZ(KBDIM)          , ZRNEB(KBDIM)&
  &,  ZR21(KBDIM)           , ZR22(KBDIM)&
  &,  ZR23(KBDIM)           , ZSS1(KBDIM)&
  &,  ZTO1(KBDIM)           , ZTR(KBDIM,2,KLEV+1)&
  &,  ZTR1(KBDIM)           , ZTR2(KBDIM)&
  &,  ZW(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IKL, IKLP1, JA, JAJ, JK, JKM1, JL

!     LOCAL REAL SCALARS
REAL(DP):: ZBMU0, ZBMU1, ZCORAE, ZCORCD, ZDEN, ZDEN1,&
          &ZFACOA, ZFACOC, ZGAP, ZMU1, ZMUE, ZRE11, &
          &ZTO, ZWW



!     ------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------


DO JK = 1 , KLEV+1
  DO JA = 1 , 6
    DO JL = KIDIA,KFDIA
      PRJ(JL,JA,JK) = 0._DP
      PRK(JL,JA,JK) = 0._DP
    ENDDO
  ENDDO
ENDDO


!     ------------------------------------------------------------------

!*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
!                ----------------------------------------------


DO JL = KIDIA,KFDIA
  ZR23(JL) = 0._DP
  ZC1I(JL,KLEV+1) = 0._DP
  ZCLEAR(JL) = 1._DP
  ZCLOUD(JL) = 0._DP
ENDDO

JK = 1
IKL = KLEV+1 - JK
IKLP1 = IKL + 1
DO JL = KIDIA,KFDIA
  ZFACOA = 1._DP - PPIZAZ(JL,IKL)*PCGAZ(JL,IKL)*PCGAZ(JL,IKL)
  ZFACOC = 1._DP - POMEGA(JL,KNU,IKL) * PCG(JL,KNU,IKL)* PCG(JL,KNU,IKL)
  ZCORAE = ZFACOA * PTAUAZ(JL,IKL) * PSEC(JL)
  ZCORCD = ZFACOC * PTAU(JL,KNU,IKL) * PSEC(JL)
  ZCORAE = MIN(200._DP,ZCORAE)
  ZCORCD = MIN(200._DP,ZCORCD)
  ZR21(JL) = EXP(-ZCORAE   )
  ZR22(JL) = EXP(-ZCORCD   )
  ZSS1(JL) = PCLD(JL,IKL)*(1._DP-ZR21(JL)*ZR22(JL))&
   &+ (1._DP-PCLD(JL,IKL))*(1._DP-ZR21(JL))
  ZCLEQ(JL,IKL) = ZSS1(JL)

  IF (NOVLP == 1) THEN
!* maximum-random      
    ZCLEAR(JL) = ZCLEAR(JL)&
     &*(1._DP-MAX(ZSS1(JL),ZCLOUD(JL)))&
     &/(1._DP-MIN(ZCLOUD(JL),1._DP-EPSILON(1._DP)))
    ZC1I(JL,IKL) = 1._DP - ZCLEAR(JL)
    ZCLOUD(JL) = ZSS1(JL)
  ELSEIF (NOVLP == 2) THEN
!* maximum
    ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
    ZC1I(JL,IKL) = ZCLOUD(JL)
  ELSEIF (NOVLP == 3) THEN
!* random
    ZCLEAR(JL) = ZCLEAR(JL)*(1._DP - ZSS1(JL))
    ZCLOUD(JL) = 1._DP - ZCLEAR(JL)
    ZC1I(JL,IKL) = ZCLOUD(JL)
  ENDIF
ENDDO

DO JK = 2 , KLEV
  IKL = KLEV+1 - JK
  IKLP1 = IKL + 1
  DO JL = KIDIA,KFDIA
    ZFACOA = 1._DP - PPIZAZ(JL,IKL)*PCGAZ(JL,IKL)*PCGAZ(JL,IKL)
    ZFACOC = 1._DP - POMEGA(JL,KNU,IKL) * PCG(JL,KNU,IKL)* PCG(JL,KNU,IKL)
    ZCORAE = ZFACOA * PTAUAZ(JL,IKL) * PSEC(JL)
    ZCORCD = ZFACOC * PTAU(JL,KNU,IKL) * PSEC(JL)
    ZCORAE = MIN(200._DP,ZCORAE)
    ZCORCD = MIN(200._DP,ZCORCD)
    ZR21(JL) = EXP(-ZCORAE   )
    ZR22(JL) = EXP(-ZCORCD   )
    ZSS1(JL) = PCLD(JL,IKL)*(1._DP-ZR21(JL)*ZR22(JL))&
     &+ (1._DP-PCLD(JL,IKL))*(1._DP-ZR21(JL))
    ZCLEQ(JL,IKL) = ZSS1(JL)

    IF (NOVLP == 1) THEN
!* maximum-random      
      ZCLEAR(JL) = ZCLEAR(JL)&
       &*(1._DP-MAX(ZSS1(JL),ZCLOUD(JL)))&
       &/(1._DP-MIN(ZCLOUD(JL),1._DP-EPSILON(1._DP)))
      ZC1I(JL,IKL) = 1._DP - ZCLEAR(JL)
      ZCLOUD(JL) = ZSS1(JL)
    ELSEIF (NOVLP == 2) THEN
!* maximum
      ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
      ZC1I(JL,IKL) = ZCLOUD(JL)
    ELSEIF (NOVLP == 3) THEN
!* random
      ZCLEAR(JL) = ZCLEAR(JL)*(1._DP - ZSS1(JL))
      ZCLOUD(JL) = 1._DP - ZCLEAR(JL)
      ZC1I(JL,IKL) = ZCLOUD(JL)
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
!                -----------------------------------------------


DO JL = KIDIA,KFDIA
  PRAY1(JL,KLEV+1) = 0._DP
  PRAY2(JL,KLEV+1) = 0._DP
  PREFZ(JL,2,1) = PALBD(JL,KNU)
  PREFZ(JL,1,1) = PALBD(JL,KNU)
  PTRA1(JL,KLEV+1) = 1._DP
  PTRA2(JL,KLEV+1) = 1._DP
ENDDO

DO JK = 2 , KLEV+1
  JKM1 = JK-1
  DO JL = KIDIA,KFDIA
    ZRNEB(JL)= PCLD(JL,JKM1)
    ZRE1(JL)=0._DP
    ZTR1(JL)=0._DP
    ZRE2(JL)=0._DP
    ZTR2(JL)=0._DP


!     ------------------------------------------------------------------

!*         3.1  EQUIVALENT ZENITH ANGLE
!               -----------------------


    ZMUE = (1._DP-ZC1I(JL,JK)) * PSEC(JL)+ ZC1I(JL,JK) * 1.66_DP
    PRMUE(JL,JK) = 1._DP/ZMUE


!     ------------------------------------------------------------------

!*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
!               ----------------------------------------------------


    ZGAP = PCGAZ(JL,JKM1)
    ZBMU0 = 0.5_DP - 0.75_DP * ZGAP / ZMUE
    ZWW = PPIZAZ(JL,JKM1)
    ZTO = PTAUAZ(JL,JKM1)
    ZDEN = 1._DP + (1._DP - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE &
     &+ (1-ZWW) * (1._DP - ZWW +2._DP*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE
    PRAY1(JL,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN
    PTRA1(JL,JKM1) = 1._DP / ZDEN

    ZMU1 = 0.5_DP
    ZBMU1 = 0.5_DP - 0.75_DP * ZGAP * ZMU1
    ZDEN1= 1._DP + (1._DP - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1 &
     &+ (1-ZWW) * (1._DP - ZWW +2._DP*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1
    PRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1
    PTRA2(JL,JKM1) = 1._DP / ZDEN1


!     ------------------------------------------------------------------

!*         3.3  EFFECT OF CLOUD LAYER
!               ---------------------


    ZW(JL) = POMEGA(JL,KNU,JKM1)
    ZTO1(JL) = PTAU(JL,KNU,JKM1)/ZW(JL)+ PTAUAZ(JL,JKM1)/PPIZAZ(JL,JKM1)
    ZR21(JL) = PTAU(JL,KNU,JKM1) + PTAUAZ(JL,JKM1)
    ZR22(JL) = PTAU(JL,KNU,JKM1) / ZR21(JL)
    ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)&
     &+ (1._DP - ZR22(JL)) * PCGAZ(JL,JKM1)
    IF (ZW(JL) == 1._DP .AND. PPIZAZ(JL,JKM1) == 1._DP) THEN
      ZW(JL)=1._DP
    ELSE
      ZW(JL) = ZR21(JL) / ZTO1(JL)
    ENDIF
    ZREF(JL) = PREFZ(JL,1,JKM1)
    ZRMUZ(JL) = PRMUE(JL,JK)
  ENDDO

  CALL SWDE ( KIDIA, KFDIA , KBDIM &
   &, ZGG  , ZREF  , ZRMUZ , ZTO1 , ZW &
   &, ZRE1 , ZRE2  , ZTR1  , ZTR2      )

  DO JL = KIDIA,KFDIA

    PREFZ(JL,1,JK) = (1._DP-ZRNEB(JL)) * (PRAY1(JL,JKM1)&
     &+ PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)&
     &* PTRA2(JL,JKM1)&
     &/ (1._DP-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))&
     &+ ZRNEB(JL) * ZRE2(JL)

    ZTR(JL,1,JKM1) = ZRNEB(JL) * ZTR2(JL) + (PTRA1(JL,JKM1)&
     &/ (1._DP-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))&
     &* (1._DP-ZRNEB(JL))

    PREFZ(JL,2,JK) = (1._DP-ZRNEB(JL)) * (PRAY1(JL,JKM1)&
     &+ PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)&
     &* PTRA2(JL,JKM1) )&
     &+ ZRNEB(JL) * ZRE1(JL)

    ZTR(JL,2,JKM1) = ZRNEB(JL) * ZTR1(JL)+ PTRA1(JL,JKM1) * (1._DP-ZRNEB(JL))

  ENDDO
ENDDO
DO JL = KIDIA,KFDIA
  ZMUE = (1._DP-ZC1I(JL,1))*PSEC(JL)+ZC1I(JL,1)*1.66_DP
  PRMUE(JL,1)=1._DP/ZMUE
  PTRCLD(JL)=1._DP-ZC1I(JL,1)
ENDDO


!     ------------------------------------------------------------------

!*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
!                 -------------------------------------------------


IF (KNU == 1) THEN
  JAJ = 2
  DO JL = KIDIA,KFDIA
    PRJ(JL,JAJ,KLEV+1) = 1._DP
    PRK(JL,JAJ,KLEV+1) = PREFZ(JL, 1,KLEV+1)
  ENDDO

  DO JK = 1 , KLEV
    IKL = KLEV+1 - JK
    IKLP1 = IKL + 1
    DO JL = KIDIA,KFDIA
      ZRE11= PRJ(JL,JAJ,IKLP1) * ZTR(JL,  1,IKL)
      PRJ(JL,JAJ,IKL) = ZRE11
      PRK(JL,JAJ,IKL) = ZRE11 * PREFZ(JL,  1,IKL)
    ENDDO
  ENDDO

ELSE

  DO JAJ = 1 , 2
    DO JL = KIDIA,KFDIA
      PRJ(JL,JAJ,KLEV+1) = 1._DP
      PRK(JL,JAJ,KLEV+1) = PREFZ(JL,JAJ,KLEV+1)
    ENDDO

    DO JK = 1 , KLEV
      IKL = KLEV+1 - JK
      IKLP1 = IKL + 1
      DO JL = KIDIA,KFDIA
        ZRE11= PRJ(JL,JAJ,IKLP1) * ZTR(JL,JAJ,IKL)
        PRJ(JL,JAJ,IKL) = ZRE11
        PRK(JL,JAJ,IKL) = ZRE11 * PREFZ(JL,JAJ,IKL)
      ENDDO
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SWR
