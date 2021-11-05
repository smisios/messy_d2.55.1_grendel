SUBROUTINE SW &
 &( KIDIA, KFDIA , KBDIM  , KLEV &
 &, KAER , KNEWAER,KAERH &
 &, PSCT , PCARDI, PPSOL , PALBD, PALBP , PWV, PQS  &
 &, PRMU0, PCG   , PCLEAR,PCLDFRAC,PDP  , POMEGA, POZ, PPMB &
 &, PTAU , PTAVE , PAER &
 &, PHEAT, PFDOWN, PFUP &
 &, PCEAT, PCDOWN, PCUP &
 &, PFDNN, PFDNV , PFUPN, PFUPV &
 &, PCDNN, PCDNV , PCUPN, PCUPV &
 &, PSUDU &
 &)

!**** *SW* - COMPUTES THE SHORTWAVE RADIATION FLUXES.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
!     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

!**   INTERFACE.
!     ----------

!          *SW* IS CALLED FROM *RADLSW*


!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES ABSORBER AMOUNTS                 (SWU)
!          2. COMPUTES FLUXES IN U.V./VISIBLE  SPECTRAL INTERVAL (SW1S)
!          3. COMPUTES FLUXES IN NEAR-INFRARED SPECTRAL INTERVAL (SWNI)

!     EXTERNALS.
!     ----------

!          *SWU*, *SW1S*, *SWNI*

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
!        95-01-01   J.-J. MORCRETTE  Direct/Diffuse Albedo
!        95-12-07   J.-J. MORCRETTE  Near-Infrared in nsw-1 Intervals
!        990128     JJMorcrette      sunshine duration
!        99-05-25   JJMorcrette      Revised aerosols

!     ------------------------------------------------------------------


USE MO_KIND     , ONLY : DP

USE MO_CONSTANTS,    ONLY : RG=>G,RD
USE MO_TIME_CONTROL, ONLY : NDAYLEN
USE MO_SW       ,    ONLY : NSW


IMPLICIT NONE

REAL(DP) :: RCDAY

!     DUMMY INTEGER SCALARS
INTEGER :: KAER
INTEGER :: KNEWAER
INTEGER :: KFDIA
INTEGER :: KIDIA
INTEGER :: KLEV
INTEGER :: KBDIM

!     DUMMY REAL SCALARS
REAL(DP):: PCARDI
REAL(DP):: PSCT



!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL(DP):: PPSOL(KBDIM), PAER(KBDIM,KLEV,KAER+KNEWAER),PRMU0(KBDIM)&
  &,  PWV(KBDIM,KLEV),PQS(KBDIM,KLEV)

REAL(DP):: PALBD(KBDIM,NSW)      , PALBP(KBDIM,NSW)&
  &,  PCG(KBDIM,NSW,KLEV)   , PCLEAR(KBDIM), PCLDFRAC(KBDIM,KLEV)&
  &,  PDP(KBDIM,KLEV)  &
  &,  POMEGA(KBDIM,NSW,KLEV), POZ(KBDIM,KLEV)&
  &,  PPMB(KBDIM,KLEV+1)&
  &,  PTAU(KBDIM,NSW,KLEV)  , PTAVE(KBDIM,KLEV)

INTEGER :: KAERH(KBDIM,KLEV)

REAL(DP):: PHEAT(KBDIM,KLEV), PFDOWN(KBDIM,KLEV+1), PFUP(KBDIM,KLEV+1),&
     &PFUPV(KBDIM), PFUPN(KBDIM), PFDNV(KBDIM), PFDNN(KBDIM)&
  &,  PCEAT(KBDIM,KLEV), PCDOWN(KBDIM,KLEV+1), PCUP(KBDIM,KLEV+1)&
  &,  PCUPV(KBDIM), PCUPN(KBDIM), PCDNV(KBDIM), PCDNN(KBDIM)&
  &,  PSUDU(KBDIM)

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL(DP):: ZAKI(KBDIM,2,NSW)&
  &,  ZDSIG(KBDIM,KLEV)   , ZFACT(KBDIM)&
  &,  ZFD(KBDIM,KLEV+1)   , ZCD(KBDIM,KLEV+1)&
  &,  ZCDOWN(KBDIM,KLEV+1), ZCDNIR(KBDIM,KLEV+1)&
  &,  ZFDOWN(KBDIM,KLEV+1), ZFDNIR(KBDIM,KLEV+1)&
  &,  ZFU(KBDIM,KLEV+1)   , ZCU(KBDIM,KLEV+1)&
  &,  ZCUP(KBDIM,KLEV+1)  , ZCUNIR(KBDIM,KLEV+1)&
  &,  ZFUP(KBDIM,KLEV+1)  , ZFUNIR(KBDIM,KLEV+1)&
  &,  ZRMU(KBDIM)         , ZSEC(KBDIM)         &
  &,  ZSUDU1(KBDIM)       , ZSUDU2(KBDIM)       , ZSUDU2T(KBDIM)&
  &,  ZUD(KBDIM,5,KLEV+1)

!     LOCAL INTEGER SCALARS
INTEGER :: INU, JK, JKL, JL, JNU

!     LOCAL REAL SCALARS
REAL(DP):: ZDCNET, ZDFNET


!     ------------------------------------------------------------------

!*         1.     ABSORBER AMOUNTS AND OTHER USEFUL QUANTITIES
!                 --------------------------------------------

RCDAY = REAL(NDAYLEN,dp)*RG/3.5_dp/RD

CALL SWU ( KIDIA,KFDIA ,KBDIM ,KLEV &
         &, PSCT ,PCARDI,PPMB ,PPSOL &
         &, PRMU0,PTAVE ,PWV &
         &, ZAKI ,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD )


!     ------------------------------------------------------------------

!*         2.     INTERVAL (0.25-0.68 MICRON): U.V. AND VISIBLE
!                 ---------------------------------------------


INU = 1

CALL SW1S &
  &( KIDIA , KFDIA , KBDIM , KLEV , KAER , KNEWAER, KAERH, INU &
  &,  PAER , PALBD , PALBP, PCG  , PCLDFRAC , PCLEAR       &
  &,  ZDSIG, POMEGA, POZ  , ZRMU , ZSEC , PTAU  , ZUD  &
  &,  ZFD  , ZFU   , ZCD  , ZCU  , ZSUDU1 &
  &)


!     ------------------------------------------------------------------

!*         3.     INTERVAL (0.68-4.00 MICRON): NEAR-INFRARED
!                 ------------------------------------------


DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    ZFDOWN(JL,JK)=0._DP
    ZFUP  (JL,JK)=0._DP
    ZCDOWN(JL,JK)=0._DP
    ZCUP  (JL,JK)=0._DP
    ZSUDU2T(JL)  =0._DP
  ENDDO
ENDDO

DO JNU = 2 , NSW

  CALL SWNI &
   &(  KIDIA ,KFDIA , KBDIM , KLEV , KAER , KNEWAER, KAERH, JNU &
   &,  PAER  ,ZAKI  , PALBD, PALBP, PCG  , PCLDFRAC, PCLEAR &
   &,  ZDSIG ,POMEGA, POZ  , ZRMU , ZSEC , PTAU, ZUD    &
   &,  PWV   ,PQS &
   &,  ZFDNIR,ZFUNIR,ZCDNIR,ZCUNIR,ZSUDU2 &
   &)

  DO JK = 1 , KLEV+1
    DO JL = KIDIA,KFDIA
      ZFDOWN(JL,JK)=ZFDOWN(JL,JK)+ZFDNIR(JL,JK)
      ZFUP  (JL,JK)=ZFUP  (JL,JK)+ZFUNIR(JL,JK)
      ZCDOWN(JL,JK)=ZCDOWN(JL,JK)+ZCDNIR(JL,JK)
      ZCUP  (JL,JK)=ZCUP  (JL,JK)+ZCUNIR(JL,JK)
    ENDDO
  ENDDO
  DO JL = KIDIA,KFDIA
    ZSUDU2T(JL)=ZSUDU2T(JL)+ZSUDU2(JL)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         4.     FILL THE DIAGNOSTIC ARRAYS
!                 --------------------------


DO JL = KIDIA,KFDIA
  PFDNN(JL)=ZFDOWN(JL,1)*ZFACT(JL)
  PFDNV(JL)=ZFD(JL,1)*ZFACT(JL)
  PFUPN(JL)=ZFUP(JL,KLEV+1)*ZFACT(JL)
  PFUPV(JL)=ZFU(JL,KLEV+1)*ZFACT(JL)

  PCDNN(JL)=ZCDOWN(JL,1)*ZFACT(JL)
  PCDNV(JL)=ZCD(JL,1)*ZFACT(JL)
  PCUPN(JL)=ZCUP(JL,KLEV+1)*ZFACT(JL)
  PCUPV(JL)=ZCU(JL,KLEV+1)*ZFACT(JL)

  PSUDU(JL)=(ZSUDU1(JL)+ZSUDU2T(JL))*ZFACT(JL)
ENDDO

DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PFUP(JL,JK)   = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
    PFDOWN(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
    PCUP(JL,JK)   = (ZCUP(JL,JK)   + ZCU(JL,JK)) * ZFACT(JL)
    PCDOWN(JL,JK) = (ZCDOWN(JL,JK) + ZCD(JL,JK)) * ZFACT(JL)
  ENDDO
ENDDO

DO JKL = 1 , KLEV
  JK = KLEV+1 - JKL
  DO JL = KIDIA,KFDIA
    ZDFNET = PFUP(JL,JK+1) - PFDOWN(JL,JK+1)-PFUP(JL,JK  ) + PFDOWN(JL,JK  )
    PHEAT(JL,JK) = RCDAY * ZDFNET / PDP(JL,JKL)
    ZDCNET = PCUP(JL,JK+1) - PCDOWN(JL,JK+1)-PCUP(JL,JK  ) + PCDOWN(JL,JK  )
    PCEAT(JL,JK) = RCDAY * ZDCNET / PDP(JL,JKL)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SW
