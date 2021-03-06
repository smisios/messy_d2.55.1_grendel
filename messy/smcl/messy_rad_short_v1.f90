! *************************************************************************
MODULE messy_rad_short_v1
! *************************************************************************

  USE messy_main_constants_mem, ONLY: DP
  USE messy_rad_fubrad_mem,     ONLY: lfubrad
  USE messy_rad_fubrad,         ONLY: middle_atmosphere_downward_flux  &
                                    , middle_atmosphere_upward_flux
  USE messy_rad_short_cmn ! NSW, NOVLP, rsun_scale, rad_sw_initialize, ...

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! SUBROUTINES
  !
  interface rad_sw_swclr
     module procedure rad_sw_swclr
  end interface
  !
  PUBLIC :: rad_sw_SW

CONTAINS

  ! ------------------------------------------------------------------
  SUBROUTINE rad_sw_SW &
       &( NDAYLEN &
       &, KIDIA, KFDIA , KBDIM  , KLEV                  &
       &, PAOT_GAMMA,PAOT_OMEGA, PAOT_TAU               &
       &, PSCT , PCO2, PPSOL , PALBD, PALBP , PWV, PQS  & 
       &, PRMU0, PCG,  PCLEAR,PCLDFRAC,PDP  , POMEGA, POZ, PPMB &
       &, PTAU , PTAVE &
       &, PFDOWN, PFUP &
       &, PCDOWN, PCUP &
       &, PFDNIR, PFUNIR &           ! solar flux NIR, up, down
       &, PCDNIR, PCUNIR &           ! solar flux NIR, up, down clear sky
       &, PFDSW1, PFUSW1 &           ! solar flux UVVIS, up, down
       &, PCDSW1, PCUSW1 &           ! solar flux UVVIS, up, down clear sky
       )

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
    !
    !        THIS SUBROUTINE COMPUTES THE SHORTWAVE RADIATION 
    !        TRANSMISSIVITIES IN FOUR SPECTRAL INTERVALS FOLLOWING 
    !        FOUQUART AND BONNEL (1980).
    !        THE OUTPUT HAS CHANGED AS TRANSMISSIVITIES ARE NEEDED 
    !        RATHER THAN THE SHORTWAVE FLUXES IN THE SMIL.
    !        ALL OBSOLETE VARIABLES ARE INACTIVE NOW.
    !
    !     ------------------------------------------------------------------

    USE messy_main_constants_mem, ONLY: RG=>g, rd

    IMPLICIT NONE
    INTRINSIC :: REAL

    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KLEV
    INTEGER :: KBDIM

    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    INTEGER, INTENT(IN) :: NDAYLEN

    REAL(DP), DIMENSION(KBDIM,KLEV) :: PCO2
    REAL(DP) :: PSCT

    REAL(DP) :: PPSOL(KBDIM), PRMU0(KBDIM)      &
              , PWV(KBDIM,KLEV),PQS(KBDIM,KLEV)
    REAL(DP) :: PAOT_GAMMA(KBDIM,KLEV,NSW) &
              , PAOT_OMEGA(KBDIM,KLEV,NSW) &
              , PAOT_TAU(KBDIM,KLEV,NSW)

    REAL(DP):: PALBD(KBDIM,NSW), PALBP(KBDIM,NSW)&
         &,  PCG(KBDIM,NSW,KLEV), PCLEAR(KBDIM), PCLDFRAC(KBDIM,KLEV)&
         &,  PDP(KBDIM,KLEV)  &
         &,  POMEGA(KBDIM,NSW,KLEV), POZ(KBDIM,KLEV)&
         &,  PPMB(KBDIM,KLEV+1)&
         &,  PTAU(KBDIM,NSW,KLEV)  , PTAVE(KBDIM,KLEV)

    REAL(DP):: PFDOWN(KBDIM,KLEV+1), PFUP(KBDIM,KLEV+1) &
         &,    PCDOWN(KBDIM,KLEV+1), PCUP(KBDIM,KLEV+1) 

    ! All output variables are transmissivities, as needed in the SMIL.
    REAL(DP), DIMENSION(KBDIM,KLEV+1) :: PFDNIR ! NIR, downward, total sky
    REAL(DP), DIMENSION(KBDIM,KLEV+1) :: PFUNIR ! NIR, upward, total sky
    REAL(DP), DIMENSION(KBDIM,KLEV+1) :: PCDNIR ! NIR, downward, clear sky
    REAL(DP), DIMENSION(KBDIM,KLEV+1) :: PCUNIR ! NIR, upward, clear sky
    REAL(DP), DIMENSION(KBDIM,KLEV+1) :: PFDSW1 ! UVVIS, downward, total sky
    REAL(DP), DIMENSION(KBDIM,KLEV+1) :: PFUSW1 ! UVVIS, upward, total sky
    REAL(DP), DIMENSION(KBDIM,KLEV+1) :: PCDSW1 ! UVVIS, downward, clear sky
    REAL(DP), DIMENSION(KBDIM,KLEV+1) :: PCUSW1 ! UVVIS, upward, clear sky

    !     ------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    REAL(DP):: ZAKI(KBDIM,2,NSW)&
         &,  ZDSIG(KBDIM,KLEV)   , ZFACT(KBDIM)&
         &,  ZCDNIR(KBDIM,KLEV+1)&
         &,  ZFDNIR(KBDIM,KLEV+1)&
         &,  ZCUNIR(KBDIM,KLEV+1)&
         &,  ZFUNIR(KBDIM,KLEV+1)&
         &,  ZRMU(KBDIM)         , ZSEC(KBDIM)         &
         &,  ZSUDU1(KBDIM)       , ZSUDU2(KBDIM)       &
         &,  ZUD(KBDIM,5,KLEV+1)

    ! LOCAL INTEGER SCALARS
    INTEGER :: INU, JK, JNU

    ! LOCAL REAL SCALARS
    REAL(DP) :: RCDAY

    !     ------------------------------------------------------------------

    !*         1.     ABSORBER AMOUNTS AND OTHER USEFUL QUANTITIES
    !                 --------------------------------------------

    RCDAY = REAL(NDAYLEN)*RG/3.5_dp/RD

    CALL rad_sw_SWU ( KIDIA,KFDIA ,KBDIM ,KLEV &
         &, PSCT ,PCO2,PPMB ,PPSOL &
         &, PRMU0,PTAVE ,PWV &
         &, ZAKI ,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD )

    !     ------------------------------------------------------------------

    !*         2.     INTERVAL (0.25-0.68 MICRON): U.V. AND VISIBLE
    !                 ---------------------------------------------

    INU = 1

    CALL rad_sw_SW1S &
         &(  KIDIA , KFDIA , KBDIM , KLEV , INU &
         &,  PTAVE, PAOT_GAMMA, PAOT_OMEGA,  PAOT_TAU &
         &,  PALBD , PALBP, PCG  , PCLDFRAC , PCLEAR       &
         &,  ZDSIG, POMEGA, POZ  , ZRMU , ZSEC , PTAU  , ZUD  &
         &,  PFDSW1, PFUSW1, PCDSW1, PCUSW1, ZSUDU1 &
         &)
    !     ------------------------------------------------------------------

    !*         3.     INTERVAL (0.68-4.00 MICRON): NEAR-INFRARED
    !                 ------------------------------------------

    PFDNIR(KIDIA:KFDIA,:) = 0._DP
    PFUNIR(KIDIA:KFDIA,:) = 0._DP
    PCDNIR(KIDIA:KFDIA,:) = 0._DP
    PCUNIR(KIDIA:KFDIA,:) = 0._DP

    DO JNU = 2 , NSW

       CALL rad_sw_SWNI &
            &(  KIDIA ,KFDIA , KBDIM , KLEV , JNU &
            &,  PAOT_GAMMA,PAOT_OMEGA, PAOT_TAU &
            &,  ZAKI  , PALBD, PALBP, PCG  , PCLDFRAC, PCLEAR &
            &,  ZDSIG ,POMEGA, POZ  , ZRMU , ZSEC , PTAU, ZUD    &
            &,  PWV   ,PQS &
            &,  ZFDNIR,ZFUNIR,ZCDNIR,ZCUNIR,ZSUDU2)

       DO JK=1,KLEV+1
          PFDNIR(KIDIA:KFDIA,JK) = PFDNIR(KIDIA:KFDIA,JK) + ZFDNIR(KIDIA:KFDIA,JK)
          PFUNIR(KIDIA:KFDIA,JK) = PFUNIR(KIDIA:KFDIA,JK) + ZFUNIR(KIDIA:KFDIA,JK)
          PCDNIR(KIDIA:KFDIA,JK) = PCDNIR(KIDIA:KFDIA,JK) + ZCDNIR(KIDIA:KFDIA,JK)
          PCUNIR(KIDIA:KFDIA,JK) = PCUNIR(KIDIA:KFDIA,JK) + ZCUNIR(KIDIA:KFDIA,JK)
       ENDDO
    ENDDO

    !    ------------------------------------------------------------------

    !*         4.     FILL THE DIAGNOSTIC ARRAYS
    !                 --------------------------
    ! conversion from transmissivity to fluxes is not necessary here, as this
    ! conversion is withdrawn in messy_rad_si.f90 (rad_radiation).
    DO JK = 1 , KLEV+1
       ! total transmissivity (NIR + UV) upward
       PFUP  (KIDIA:KFDIA,JK) = PFUNIR(KIDIA:KFDIA,JK) + PFUSW1(KIDIA:KFDIA,JK)
       ! total transmissivity (NIR + UV) downward
       PFDOWN(KIDIA:KFDIA,JK) = PFDNIR(KIDIA:KFDIA,JK) + PFDSW1(KIDIA:KFDIA,JK)
       ! total transmissivity (NIR + UV) upward, clear sky
       PCUP  (KIDIA:KFDIA,JK) = PCUNIR(KIDIA:KFDIA,JK) + PCUSW1(KIDIA:KFDIA,JK)
       ! total transmissivity (NIR + UV) downward, clear sky
       PCDOWN(KIDIA:KFDIA,JK) = PCDNIR(KIDIA:KFDIA,JK) + PCDSW1(KIDIA:KFDIA,JK)
    ENDDO

    !     ------------------------------------------------------------------

    RETURN
  END SUBROUTINE rad_sw_SW
  !     ------------------------------------------------------------------

  !     ------------------------------------------------------------------
!OPTIONS XOPT(HSFUN)
  SUBROUTINE rad_sw_SWU &
       &( KIDIA, KFDIA , KBDIM  , KLEV &
       &, PSCT , PCO2, PPMB , PPSOL, PRMU0, PTAVE, PWV &
       &, PAKI , PDSIG , PFACT, PRMU , PSEC , PUD &
       &)

    !**** *SWU* - SHORTWAVE RADIATION, ABSORBER AMOUNTS

    !     PURPOSE.
    !     --------
    !           COMPUTES THE ABSORBER AMOUNTS USED IN SHORTWAVE RADIATION
    !     CALCULATIONS

    !**   INTERFACE.
    !     ----------
    !          *SWU* IS CALLED BY *SW*


    !        IMPLICIT ARGUMENTS :
    !        --------------------

    !     ==== INPUTS ===
    !     ==== OUTPUTS ===

    !     METHOD.
    !     -------

    !          1. COMPUTES ABSORBER AMOUNTS WITH TEMPERATURE AND PRESSURE
    !     SCALING.

    !     EXTERNALS.
    !     ----------

    !          *SWTT*

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

    !     ------------------------------------------------------------------

    IMPLICIT NONE
    INTRINSIC :: EPSILON, LOG, MAX, SQRT 

    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KLEV
    INTEGER :: KBDIM

    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    REAL(DP):: PPMB(KBDIM,KLEV+1), PPSOL(KBDIM)&
         &,  PRMU0(KBDIM)      , PTAVE(KBDIM,KLEV) , PWV(KBDIM,KLEV)

    REAL(DP):: PAKI(KBDIM,2,NSW)&
         &,  PDSIG(KBDIM,KLEV) , PFACT(KBDIM)      , PRMU(KBDIM)&
         &,  PSEC(KBDIM)       , PUD(KBDIM,5,KLEV+1)

    REAL(DP), DIMENSION(KBDIM,KLEV) :: PCO2
    REAL(DP):: PSCT

    !     ------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    INTEGER :: IIND(2)
    REAL(DP):: ZN175(KBDIM), ZN190(KBDIM), ZO175(KBDIM)&
         &,  ZO190(KBDIM), ZSIGN(KBDIM)&
         &,  ZR(KBDIM,2) , ZSIGO(KBDIM), ZUD(KBDIM,2)

    !     LOCAL INTEGER SCALARS
    INTEGER :: JA, JK, JKL, JKP1, JL, JNU

    !     LOCAL REAL SCALARS
    REAL(DP):: ZDSCO2, ZDSH2O, ZFPPW, ZRTH, ZRTU, ZWH2O

    !     DUMMY REAL SCALARS
    REAL(DP):: PCARDI

    !     ------------------------------------------------------------------

    !*         1.     COMPUTES AMOUNTS OF ABSORBERS
    !                 -----------------------------

    IIND(1)=1
    IIND(2)=2

    !*         1.1    INITIALIZES QUANTITIES
    !                 ----------------------

    DO JL = KIDIA,KFDIA
       PUD(JL,1,KLEV+1)=0._DP
       PUD(JL,2,KLEV+1)=0._DP
       PUD(JL,3,KLEV+1)=0._DP
       PUD(JL,4,KLEV+1)=0._DP
       PUD(JL,5,KLEV+1)=0._DP
       PFACT(JL)= PRMU0(JL) * PSCT
       PRMU(JL)=SQRT(1224._DP* PRMU0(JL) * PRMU0(JL) + 1._DP) / 35._DP
       PSEC(JL)=1._DP/PRMU(JL)
    ENDDO

    !*          1.3    AMOUNTS OF ABSORBERS
    !                  --------------------

    DO JL= KIDIA,KFDIA
       ZUD(JL,1) = 0._DP
       ZUD(JL,2) = 0._DP
       ZO175(JL) = PPSOL(JL)** RPDU1
       ZO190(JL) = PPSOL(JL)** RPDH1
       ZSIGO(JL) = PPSOL(JL)
    ENDDO

    DO JK = 1 , KLEV
       JKP1 = JK + 1
       JKL = KLEV+1 - JK
       DO JL = KIDIA,KFDIA
          ZRTH=(RTH2O/PTAVE(JL,JK))**RTDH2O
          ZRTU=(RTUMG/PTAVE(JL,JK))**RTDUMG
          ZWH2O = MAX (PWV(JL,JKL) , EPSILON(1._DP) )
          ZSIGN(JL) = 100._DP * PPMB(JL,JKP1)
          PDSIG(JL,JK) = (ZSIGO(JL) - ZSIGN(JL))/PPSOL(JL)
          ZN175(JL) = ZSIGN(JL) ** RPDU1
          ZN190(JL) = ZSIGN(JL) ** RPDH1
          ZDSCO2 = ZO175(JL) - ZN175(JL)
          ZDSH2O = ZO190(JL) - ZN190(JL)
          PUD(JL,1,JK) = RPNH * ZDSH2O * ZWH2O  * ZRTH
          PCARDI = MAX (PCO2(JL,JKL) , EPSILON(1._DP) )
          PUD(JL,2,JK) = RPNU * ZDSCO2 * PCARDI * ZRTU
          ZFPPW=1.6078_DP*ZWH2O/(1._DP+0.608_DP*ZWH2O)
          PUD(JL,4,JK)=PUD(JL,1,JK)*ZFPPW
          PUD(JL,5,JK)=PUD(JL,1,JK)*(1._DP-ZFPPW)
          ZUD(JL,1) = ZUD(JL,1) + PUD(JL,1,JK)
          ZUD(JL,2) = ZUD(JL,2) + PUD(JL,2,JK)
          ZSIGO(JL) = ZSIGN(JL)
          ZO175(JL) = ZN175(JL)
          ZO190(JL) = ZN190(JL)
       ENDDO
    ENDDO

    !*         1.4    COMPUTES CLEAR-SKY GREY ABSORPTION COEFFICIENTS
    !                 -----------------------------------------------

    DO JA = 1,2
       DO JL = KIDIA,KFDIA
          ZUD(JL,JA) = ZUD(JL,JA) * PSEC(JL)
       ENDDO
    ENDDO

    DO JNU= 2,NSW

       CALL rad_sw_SWTT1 ( KIDIA,KFDIA,KBDIM, JNU, 2, IIND &
            &, ZUD &
            &, ZR                            )

       DO JA = 1,2
          DO JL = KIDIA,KFDIA
             PAKI(JL,JA,JNU) = -LOG( ZR(JL,JA) ) / ZUD(JL,JA)
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE rad_sw_SWU
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE rad_sw_SWTT1 ( KIDIA,KFDIA,KBDIM,KNU,KABS,KIND, PU, PTR )

    !**** *SWTT1* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

    !     PURPOSE.
    !     --------
    !           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
    !     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
    !     INTERVALS.

    !**   INTERFACE.
    !     ----------
    !          *SWTT1* IS CALLED FROM *SW1S*.


    !        EXPLICIT ARGUMENTS :
    !        --------------------
    ! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
    ! KABS   :                     ; NUMBER OF ABSORBERS
    ! KIND   : (KABS)              ; INDICES OF THE ABSORBERS
    ! PU     : (KBDIM,KABS)         ; ABSORBER AMOUNT
    !     ==== OUTPUTS ===
    ! PTR    : (KBDIM,KABS)         ; TRANSMISSION FUNCTION

    !        IMPLICIT ARGUMENTS :   NONE
    !        --------------------

    !     METHOD.
    !     -------

    !          TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
    !     AND HORNER'S ALGORITHM.

    !     EXTERNALS.
    !     ----------

    !          NONE

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
    !        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 95-01-20

    !-----------------------------------------------------------------------

    IMPLICIT NONE

    !     DUMMY INTEGER SCALARS
    INTEGER :: KABS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KBDIM
    INTEGER :: KNU

    !-----------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    INTEGER :: KIND(KABS)
    REAL(DP):: PU(KBDIM,KABS)
    REAL(DP):: PTR(KBDIM,KABS)

    !-----------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    REAL(DP):: ZR1(KBDIM), ZR2(KBDIM), ZU(KBDIM)

    !     LOCAL INTEGER SCALARS
    INTEGER :: IA, JA, JL

    !-----------------------------------------------------------------------

    !*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION

    DO JA = 1,KABS
       IA=KIND(JA)
       DO JL = KIDIA,KFDIA
          ZU(JL) = PU(JL,JA)
          ZR1(JL) = APAD(KNU,IA,1) + ZU(JL) * (APAD(KNU,IA,2) + ZU(JL)&
               &* ( APAD(KNU,IA,3) + ZU(JL) * (APAD(KNU,IA,4) + ZU(JL)&
               &* ( APAD(KNU,IA,5) + ZU(JL) * (APAD(KNU,IA,6) + ZU(JL)&
               &* ( APAD(KNU,IA,7) ))))))

          ZR2(JL) = BPAD(KNU,IA,1) + ZU(JL) * (BPAD(KNU,IA,2) + ZU(JL)&
               &* ( BPAD(KNU,IA,3) + ZU(JL) * (BPAD(KNU,IA,4) + ZU(JL)&
               &* ( BPAD(KNU,IA,5) + ZU(JL) * (BPAD(KNU,IA,6) + ZU(JL)&
               &* ( BPAD(KNU,IA,7) ))))))

     !*         2.      ADD THE BACKGROUND TRANSMISSION

          PTR(JL,JA) = (ZR1(JL)/ZR2(JL)) * (1._DP-D(KNU,IA)) + D(KNU,IA)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE rad_sw_SWTT1
  ! -----------------------------------------------------------------------

  ! -----------------------------------------------------------------------
  SUBROUTINE rad_sw_SW1S &
       &( KIDIA , KFDIA , KBDIM , KLEV , KNU &
       &, PTAVE, PAOT_GAMMA,PAOT_OMEGA, PAOT_TAU &    ! fb_mk_20150925 added PTAVE
       &, PALBD , PALBP, PCG  , PCLD , PCLEAR &
       &, PDSIG , POMEGA, POZ  , PRMU , PSEC , PTAU  , PUD  &
       &, PFD   , PFU   , PCD  , PCU  , PSUDU1 &
       &)

    !**** *SW1S* - SHORTWAVE RADIATION, FIRST SPECTRAL INTERVAL

    !     PURPOSE.
    !     --------

    !          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
    !     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

    !**   INTERFACE.
    !     ----------

    !          *SW1S* IS CALLED FROM *SW*.


    !        IMPLICIT ARGUMENTS :
    !        --------------------

    !     ==== INPUTS ===
    !     ==== OUTPUTS ===

    !     METHOD.
    !     -------

    !          1. COMPUTES QUANTITIES FOR THE CLEAR-SKY FRACTION OF THE
    !     COLUMN
    !          2. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
    !     CONTINUUM SCATTERING
    !          3. MULTIPLY BY OZONE TRANSMISSION FUNCTION

    !     EXTERNALS.
    !     ----------

    !          *SWCLR*, *SWR*, *SWTT*

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
    !        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
    !        96-01-15   J.-J. MORCRETTE    SW in nsw SPECTRAL INTERVALS 
    !        990128     JJMorcrette        sunshine duration
    !        99-05-25   JJMorcrette        Revised aerosols

    !     ------------------------------------------------------------------

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

    REAL(DP)::  PALBD(KBDIM,NSW)      , PALBP(KBDIM,NSW)&
         &,  PCG(KBDIM,NSW,KLEV)   , PCLD(KBDIM,KLEV) &
         &,  PCLEAR(KBDIM)&
         &,  PDSIG(KBDIM,KLEV)&
         &,  POMEGA(KBDIM,NSW,KLEV), POZ(KBDIM,KLEV)&
         &,  PRMU(KBDIM)           , PSEC(KBDIM)&
         &,  PTAU(KBDIM,NSW,KLEV)  , PUD(KBDIM,5,KLEV+1)

    REAL(DP):: PAOT_GAMMA(KBDIM,KLEV,NSW),PAOT_OMEGA(KBDIM,KLEV,NSW), &
         &  PAOT_TAU(KBDIM,KLEV,NSW)

    REAL(DP):: PTAVE(KBDIM,KLEV)
    REAL(DP):: PFD(KBDIM,KLEV+1)     , PFU(KBDIM,KLEV+1)&
         &,  PCD(KBDIM,KLEV+1)     , PCU(KBDIM,KLEV+1)&
         &,  PSUDU1(KBDIM)

    !     ------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    INTEGER :: IIND(6)

    REAL(DP):: ZCGAZ(KBDIM,KLEV)&
         &,  ZDIFF(KBDIM)        , ZDIRF(KBDIM)        &
         &,  ZDIFT(KBDIM)        , ZDIRT(KBDIM)        &
         &,  ZPIZAZ(KBDIM,KLEV)&
         &,  ZRAYL(KBDIM), ZRAY1(KBDIM,KLEV+1), ZRAY2(KBDIM,KLEV+1)&
         &,  ZREFZ(KBDIM,2,KLEV+1)&
         &,  ZRJ(KBDIM,6,KLEV+1), ZRJ0(KBDIM,6,KLEV+1)&
         &,  ZRK(KBDIM,6,KLEV+1), ZRK0(KBDIM,6,KLEV+1)&
         &,  ZRMUE(KBDIM,KLEV+1), ZRMU0(KBDIM,KLEV+1)&
         &,  ZR(KBDIM,6)&
         &,  ZTAUAZ(KBDIM,KLEV)&
         &,  ZTRA1(KBDIM,KLEV+1), ZTRA2(KBDIM,KLEV+1)&
         &,  ZTRCLD(KBDIM)      , ZTRCLR(KBDIM)&
         &,  ZW(KBDIM,6) &
         &,  ZRISWLEV(KBDIM,6)

    !     LOCAL INTEGER SCALARS
    INTEGER :: IKL, IKM1, JAJ, JK, JL

    REAL(DP) :: pffrac(kbdim) ! fb_kk_20061002
    INTEGER  :: ISWLEV

    !     ------------------------------------------------------------------

    !*         1.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
    !                 ----------------------- ------------------

    !*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
    !                 -----------------------------------------

    DO JL = KIDIA,KFDIA
       ZRAYL(JL) =  RRAY(KNU,1) + PRMU(JL) * (RRAY(KNU,2) + PRMU(JL)&
            &* (RRAY(KNU,3) + PRMU(JL) * (RRAY(KNU,4) + PRMU(JL)&
            &* (RRAY(KNU,5) + PRMU(JL) *  RRAY(KNU,6)       ))))
    ENDDO

    !     ------------------------------------------------------------------

    IF (lfubrad) THEN
       !                CALCULATE MIDDLE ATMOSPHERE DOWNWARD FLUX
       CALL middle_atmosphere_downward_flux (kfdia,kbdim,klev,ptave, &
            &                                      iswlev,pffrac=pffrac)

    END IF

    !*         2.    CONTINUUM SCATTERING CALCULATIONS
    !                ---------------------------------
    !*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
    !

    CALL rad_sw_SWCLR &
         &( KIDIA  , KFDIA , KBDIM  , KLEV , KNU &
         &, PAOT_GAMMA,PAOT_OMEGA, PAOT_TAU&
         &,  PALBP , PDSIG , ZRAYL, PSEC &
         &, ZCGAZ  , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0 &
         &, ZRK0   , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2, ZTRCLR &
         &)

    !*         2.2   CLOUDY FRACTION OF THE COLUMN
    !                -----------------------------

    CALL rad_sw_SWR &
         &( KIDIA ,KFDIA ,KBDIM  ,KLEV  , KNU &
         &, PAOT_GAMMA,PAOT_OMEGA, PAOT_TAU&
         &, PALBD ,PCG   ,PCLD  ,POMEGA, PSEC , PTAU &
         &, ZCGAZ ,ZPIZAZ,ZRAY1 ,ZRAY2 , ZREFZ, ZRJ  ,ZRK , ZRMUE &
         &, ZTAUAZ,ZTRA1 ,ZTRA2 ,ZTRCLD &
         &)

    !     ------------------------------------------------------------------

    !*         3.    OZONE ABSORPTION
    !                ----------------

    IIND(1)=1
    IIND(2)=2
    IIND(3)=3
    IIND(4)=1
    IIND(5)=2
    IIND(6)=3

    !*         3.1   DOWNWARD FLUXES
    !                ---------------

    JAJ = 2

    DO JL = KIDIA,KFDIA
       ZW(JL,1)=0._DP
       ZW(JL,2)=0._DP
       ZW(JL,3)=0._DP
       ZW(JL,4)=0._DP
       ZW(JL,5)=0._DP
       ZW(JL,6)=0._DP
       PFD(JL,KLEV+1)=((1._DP-PCLEAR(JL))*ZRJ(JL,JAJ,KLEV+1)&
            &+ PCLEAR(JL) *ZRJ0(JL,JAJ,KLEV+1)) * RSUN(KNU)
       PCD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1) * RSUN(KNU)
    ENDDO

    fubrad: IF (lfubrad) THEN
       DO JK = 1, KLEV
          IKL = KLEV+1-JK
          DO JL = KIDIA,KFDIA
             ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
             ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKL)/ZRMUE(JL,IKL)
             ZW(JL,3)=ZW(JL,3)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
             ZW(JL,4)=ZW(JL,4)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
             ZW(JL,5)=ZW(JL,5)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
             ZW(JL,6)=ZW(JL,6)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
          ENDDO
          IF (JK .GE. ISWLEV) THEN

             CALL rad_sw_SWTT1 ( KIDIA, KFDIA, KBDIM, KNU, 6 &
                  &, IIND &
                  &, ZW  &
                  &, ZR  )    

             IF (JK .EQ. ISWLEV) ZRISWLEV(:,:)=ZR(:,:) 

             DO JL = KIDIA,KFDIA
                ZDIFF(JL) = (1._DP-(ZRISWLEV(JL,1)-ZR(JL,1)))*&
                     &            (1._DP-(ZRISWLEV(JL,2)-ZR(JL,2)))*&
                     &            (1._DP-(ZRISWLEV(JL,3)-ZR(JL,3)))*&
                     &            PFFRAC(JL)*ZRJ(JL,JAJ,IKL)
                ZDIRF(JL) = (1._DP-(ZRISWLEV(JL,4)-ZR(JL,4)))*&
                     &            (1._DP-(ZRISWLEV(JL,5)-ZR(JL,5)))*&
                     &            (1._DP-(ZRISWLEV(JL,6)-ZR(JL,6)))*&
                     &            PFFRAC(JL)*ZRJ0(JL,JAJ,IKL)
                PFD(JL,IKL) = ((1._DP-PCLEAR(JL)) * ZDIFF(JL)&
                     &              +PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
                PCD(JL,IKL) = ZDIRF(JL) * RSUN(KNU)
             ENDDO
          ENDIF
       ENDDO
       DO JK = 0, ISWLEV-1
          IKL = KLEV+1-JK
          DO JL = KIDIA,KFDIA
             PFD(JL,IKL) = PFD(JL,KLEV+1-ISWLEV)
             PCD(JL,IKL) = PCD(JL,KLEV+1-ISWLEV)
          ENDDO
       ENDDO
    ELSE
       DO JK = 1 , KLEV
          IKL = KLEV+1-JK
          DO JL = KIDIA,KFDIA
             ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
             ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKL)/ZRMUE(JL,IKL)
             ZW(JL,3)=ZW(JL,3)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
             ZW(JL,4)=ZW(JL,4)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
             ZW(JL,5)=ZW(JL,5)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
             ZW(JL,6)=ZW(JL,6)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
          ENDDO

          CALL rad_sw_SWTT1 ( KIDIA, KFDIA, KBDIM, KNU, 6 &
               &, IIND &
               &, ZW  &
               &, ZR                          )

          DO JL = KIDIA,KFDIA
             ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZR(JL,3)*ZRJ(JL,JAJ,IKL)
             ZDIRF(JL) = ZR(JL,4)*ZR(JL,5)*ZR(JL,6)*ZRJ0(JL,JAJ,IKL)
             PFD(JL,IKL) = ((1._DP-PCLEAR(JL)) * ZDIFF(JL)&
                  &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
             PCD(JL,IKL) = ZDIRF(JL) * RSUN(KNU)
          ENDDO
       ENDDO
    ENDIF fubrad

    DO JL=KIDIA,KFDIA
       ZDIFT(JL) = ZR(JL,1)*ZR(JL,2)*ZR(JL,3)*ZTRCLD(JL)
       ZDIRT(JL) = ZR(JL,4)*ZR(JL,5)*ZR(JL,6)*ZTRCLR(JL)
       PSUDU1(JL) = ((1._DP-PCLEAR(JL)) * ZDIFT(JL)&
            &+PCLEAR(JL) * ZDIRT(JL)) * RSUN(KNU)
    ENDDO


    !*         3.2   UPWARD FLUXES
    !                -------------


    DO JL = KIDIA,KFDIA
       PFU(JL,1) = ((1._DP-PCLEAR(JL))*ZDIFF(JL)*PALBD(JL,KNU)&
            &+ PCLEAR(JL) *ZDIRF(JL)*PALBP(JL,KNU))&
            &* RSUN(KNU)
       PCU(JL,1) = ZDIRF(JL) * PALBP(JL,KNU) * RSUN(KNU)
    ENDDO

    DO JK = 2 , KLEV+1
       IKM1=JK-1
       DO JL = KIDIA,KFDIA
          ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKM1)*1.66_DP
          ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKM1)*1.66_DP
          ZW(JL,3)=ZW(JL,3)+POZ(JL,  IKM1)*1.66_DP
          ZW(JL,4)=ZW(JL,4)+PUD(JL,1,IKM1)*1.66_DP
          ZW(JL,5)=ZW(JL,5)+PUD(JL,2,IKM1)*1.66_DP
          ZW(JL,6)=ZW(JL,6)+POZ(JL,  IKM1)*1.66_DP
       ENDDO

       CALL rad_sw_SWTT1 ( KIDIA, KFDIA, KBDIM, KNU, 6 &
            &, IIND &
            &, ZW  &
            &, ZR                          )

       fubrad2: IF (LFUBRAD) THEN 
          DO JL = KIDIA,KFDIA
             ZDIFF(JL) = (1._DP-(ZRISWLEV(JL,1)-ZR(JL,1)))*&
                  &            (1._DP-(ZRISWLEV(JL,2)-ZR(JL,2)))*&
                  &            (1._DP-(ZRISWLEV(JL,3)-ZR(JL,3)))*&
                  &            PFFRAC(JL)*ZRK(JL,JAJ,JK)
             ZDIRF(JL) = (1._DP-(ZRISWLEV(JL,4)-ZR(JL,4)))*&
                  &            (1._DP-(ZRISWLEV(JL,5)-ZR(JL,5)))*&
                  &            (1._DP-(ZRISWLEV(JL,6)-ZR(JL,6)))*&
                  &            PFFRAC(JL)*ZRK0(JL,JAJ,JK)

             PFU(JL,JK) = ((1._DP-PCLEAR(JL)) * ZDIFF(JL)&
                  &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
             PCU(JL,JK) = ZDIRF(JL) * RSUN(KNU)
          ENDDO
       ELSE
          DO JL = KIDIA,KFDIA
             ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZR(JL,3)*ZRK(JL,JAJ,JK)
             ZDIRF(JL) = ZR(JL,4)*ZR(JL,5)*ZR(JL,6)*ZRK0(JL,JAJ,JK)
             PFU(JL,JK) = ((1._DP-PCLEAR(JL)) * ZDIFF(JL)&
                  &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
             PCU(JL,JK) = ZDIRF(JL) * RSUN(KNU)
          ENDDO
       ENDIF fubrad2
    ENDDO

    ! Calculate middle atmosphere upward flux
    IF (lfubrad) &
         CALL middle_atmosphere_upward_flux (kfdia,kbdim,klev,pfu,pfd, &
         pcu,pcd)
    !     ------------------------------------------------------------------

    RETURN
  END SUBROUTINE rad_sw_SW1S
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE rad_sw_SWCLR &
       &( KIDIA , KFDIA , KBDIM  , KLEV  , KNU &
       &, PAOT_GAMMA,PAOT_OMEGA, PAOT_TAU&
       &, PALBP , PDSIG , PRAYL , PSEC &
       &, PCGAZ , PPIZAZ, PRAY1 , PRAY2 , PREFZ , PRJ  &
       &, PRK   , PRMU0 , PTAUAZ, PTRA1 , PTRA2 , PTRCLR &
       &)

    !**** *SWCLR* - CLEAR-SKY COLUMN COMPUTATIONS

    !     PURPOSE.
    !     --------
    !           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
    !     CLEAR-SKY COLUMN

    !**   INTERFACE.
    !     ----------

    !          *SWCLR* IS CALLED EITHER FROM *SW1S*
    !                                OR FROM *SWNI*

    !        IMPLICIT ARGUMENTS :
    !        --------------------

    !     ==== INPUTS ===
    !     ==== OUTPUTS ===

    !     METHOD.
    !     -------


    !     EXTERNALS.
    !     ----------

    !          NONE

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    !        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 94-11-15
    !        Modified : 96-03-19 JJM-PhD (loop 107 in absence of aerosols)
    !        JJMorcrette 990128 : sunshine duration
    !        99-05-25   JJMorcrette    Revised aerosols

    !     ------------------------------------------------------------------

    IMPLICIT NONE
    INTRINSIC :: EPSILON, EXP, MAX, MIN

    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KLEV
    INTEGER :: KBDIM
    INTEGER :: KNU

    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    REAL(DP):: PALBP(KBDIM,NSW)&
         &,  PDSIG(KBDIM,KLEV)&
         &,  PRAYL(KBDIM)&
         &,  PSEC(KBDIM)

    REAL(DP) :: PAOT_GAMMA(KBDIM,KLEV,NSW) &
              , PAOT_OMEGA(KBDIM,KLEV,NSW), PAOT_TAU(KBDIM,KLEV,NSW)

    REAL(DP)::&
         &   PCGAZ(KBDIM,KLEV)     &
         &,  PPIZAZ(KBDIM,KLEV)&
         &,  PRAY1(KBDIM,KLEV+1)  , PRAY2(KBDIM,KLEV+1)&
         &,  PREFZ(KBDIM,2,KLEV+1), PRJ(KBDIM,6,KLEV+1)&
         &,  PRK(KBDIM,6,KLEV+1)  , PRMU0(KBDIM,KLEV+1)&
         &,  PTAUAZ(KBDIM,KLEV)&
         &,  PTRA1(KBDIM,KLEV+1)  , PTRA2(KBDIM,KLEV+1)&
         &,  PTRCLR(KBDIM)
    !     ------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    REAL(DP):: ZC0I(KBDIM,KLEV+1)&
         &,  ZCLE0(KBDIM,KLEV), ZCLEAR(KBDIM) &
         &,  ZR21(KBDIM)&
         &,  ZR23(KBDIM) , ZSS0(KBDIM) , ZSCAT(KBDIM)&
         &,  ZTR(KBDIM,2,KLEV+1)

    !     LOCAL INTEGER SCALARS
    INTEGER :: JA, JAJ, JK, JKL, JKLP1, JKM1, JL

    !     LOCAL REAL SCALARS
    REAL(DP):: ZBMU0, ZBMU1, ZCORAE, ZDEN, ZDEN1, ZFACOA,&
         &ZFF, ZGAP, ZGAR, ZMU1, ZMUE, ZRATIO, ZRE11, &
         &ZTO, ZTRAY, ZWW

    !     ------------------------------------------------------------------

    !*         1.    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH
    !                --------------------------------------------

    DO JK = 1 , KLEV+1
       DO JA = 1 , 6          !qqq better as parameter ???
          DO JL = KIDIA,KFDIA
             PRJ(JL,JA,JK) = 0._DP
             PRK(JL,JA,JK) = 0._DP
          ENDDO
       ENDDO
    ENDDO

    DO JK = 1 , KLEV
       PTAUAZ(KIDIA:KFDIA,JK) = PAOT_TAU(KIDIA:KFDIA,JK,KNU)
       PPIZAZ(KIDIA:KFDIA,JK) = PAOT_OMEGA(KIDIA:KFDIA,JK,KNU)
       PCGAZ(KIDIA:KFDIA,JK)  = PAOT_GAMMA(KIDIA:KFDIA,JK,KNU)
    ENDDO

    DO JK = 1 , KLEV  
       DO JL = KIDIA,KFDIA

          ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
          ZRATIO = ZTRAY / (ZTRAY + PTAUAZ(JL,JK))
          ZGAR = PCGAZ(JL,JK)
          ZFF = ZGAR * ZGAR
          PTAUAZ(JL,JK)=ZTRAY+PTAUAZ(JL,JK)*(1._DP-PPIZAZ(JL,JK)*ZFF)
          PCGAZ(JL,JK) = ZGAR * (1._DP - ZRATIO) / (1._DP + ZGAR)
          PPIZAZ(JL,JK) =ZRATIO+(1._DP-ZRATIO)*PPIZAZ(JL,JK)*(1._DP-ZFF)&
               &/ (1._DP - PPIZAZ(JL,JK) * ZFF)

       ENDDO
    ENDDO

    !     ------------------------------------------------------------------

    !*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
    !                ----------------------------------------------

    DO JL = KIDIA,KFDIA
       ZR23(JL) = 0._DP
       ZC0I(JL,KLEV+1) = 0._DP
       ZCLEAR(JL) = 1._DP
       ZSCAT(JL) = 0._DP
    ENDDO

    JK = 1
    JKL = KLEV+1 - JK

    DO JL = KIDIA,KFDIA
       ZFACOA = 1._DP - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
       ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
       ZR21(JL) = EXP(-ZCORAE   )
       ZSS0(JL) = 1._DP-ZR21(JL)
       ZCLE0(JL,JKL) = ZSS0(JL)

       IF (NOVLP == 1) THEN
          !* maximum-random      
          ZCLEAR(JL) = ZCLEAR(JL)&
               &*(1._DP-MAX(ZSS0(JL),ZSCAT(JL)))&
               &/(1._DP-MIN(ZSCAT(JL),1._DP-EPSILON(1._dp)))
          ZC0I(JL,JKL) = 1._DP - ZCLEAR(JL)
          ZSCAT(JL) = ZSS0(JL)
       ELSEIF (NOVLP == 2) THEN
          !* maximum
          ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
          ZC0I(JL,JKL) = ZSCAT(JL)
       ELSEIF (NOVLP == 3) THEN
          !* random
          ZCLEAR(JL)=ZCLEAR(JL)*(1._DP-ZSS0(JL))
          ZSCAT(JL) = 1._DP - ZCLEAR(JL)
          ZC0I(JL,JKL) = ZSCAT(JL)
       ENDIF
    ENDDO

    DO JK = 2 , KLEV
       JKL = KLEV+1 - JK
       !JKLP1 = JKL + 1
       DO JL = KIDIA,KFDIA
          ZFACOA = 1._DP - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
          ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
          ZR21(JL) = EXP(-ZCORAE   )
          ZSS0(JL) = 1._DP-ZR21(JL)
          ZCLE0(JL,JKL) = ZSS0(JL)

          IF (NOVLP == 1) THEN
             !* maximum-random      
             ZCLEAR(JL) = ZCLEAR(JL)&
                  &*(1._DP-MAX(ZSS0(JL),ZSCAT(JL)))&
                  &/(1._DP-MIN(ZSCAT(JL),1._DP-EPSILON(1._dp)))
             ZC0I(JL,JKL) = 1._DP - ZCLEAR(JL)
             ZSCAT(JL) = ZSS0(JL)
          ELSEIF (NOVLP == 2) THEN
             !* maximum
             ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
             ZC0I(JL,JKL) = ZSCAT(JL)
          ELSEIF (NOVLP == 3) THEN
             !* random
             ZCLEAR(JL)=ZCLEAR(JL)*(1._DP-ZSS0(JL))
             ZSCAT(JL) = 1._DP - ZCLEAR(JL)
             ZC0I(JL,JKL) = ZSCAT(JL)
          ENDIF
       ENDDO
    ENDDO

    !     ------------------------------------------------------------------

    !*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
    !                -----------------------------------------------

    DO JL = KIDIA,KFDIA
       PRAY1(JL,KLEV+1) = 0._DP
       PRAY2(JL,KLEV+1) = 0._DP
       PREFZ(JL,2,1) = PALBP(JL,KNU)
       PREFZ(JL,1,1) = PALBP(JL,KNU)
       PTRA1(JL,KLEV+1) = 1._DP
       PTRA2(JL,KLEV+1) = 1._DP
    ENDDO

    DO JK = 2 , KLEV+1
       JKM1 = JK-1
       DO JL = KIDIA,KFDIA

          !     -------------------------------------------------------------

          !*         3.1  EQUIVALENT ZENITH ANGLE
          !               -----------------------

          ZMUE = (1._DP-ZC0I(JL,JK)) * PSEC(JL)+ ZC0I(JL,JK) * 1.66_DP
          PRMU0(JL,JK) = 1._DP/ZMUE

          !     ------------------------------------------------------------

          !*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
          !               ----------------------------------------------------

          ZGAP = PCGAZ(JL,JKM1)
          ZBMU0 = 0.5_DP - 0.75_DP * ZGAP / ZMUE
          ZWW = PPIZAZ(JL,JKM1)
          ZTO = PTAUAZ(JL,JKM1)
          ZDEN = 1._DP + (1._DP - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE &
               &+ (1._DP-ZWW) * (1._DP - ZWW +2._DP*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE
          PRAY1(JL,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN
          PTRA1(JL,JKM1) = 1._DP / ZDEN

          ZMU1 = 0.5_DP
          ZBMU1 = 0.5_DP - 0.75_DP * ZGAP * ZMU1
          ZDEN1= 1._DP + (1._DP - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1 &
               &+ (1._DP-ZWW) * (1._DP - ZWW +2._DP*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1
          PRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1
          PTRA2(JL,JKM1) = 1._DP / ZDEN1

          PREFZ(JL,1,JK) = (PRAY1(JL,JKM1)&
               &+ PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)&
               &* PTRA2(JL,JKM1)&
               &/ (1._DP-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))

          ZTR(JL,1,JKM1) = (PTRA1(JL,JKM1)&
               &/ (1._DP-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))

          PREFZ(JL,2,JK) = (PRAY1(JL,JKM1)&
               &+ PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)&
               &* PTRA2(JL,JKM1) )

          ZTR(JL,2,JKM1) = PTRA1(JL,JKM1)

       ENDDO
    ENDDO
    DO JL = KIDIA,KFDIA
       ZMUE = (1._DP-ZC0I(JL,1))*PSEC(JL)+ZC0I(JL,1)*1.66_DP
       PRMU0(JL,1)=1._DP/ZMUE
       PTRCLR(JL)=1._DP-ZC0I(JL,1)
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
          JKL = KLEV+1 - JK
          JKLP1 = JKL + 1
          DO JL = KIDIA,KFDIA
             ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,  1,JKL)
             PRJ(JL,JAJ,JKL) = ZRE11
             PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,  1,JKL)
          ENDDO
       ENDDO

    ELSE

       DO JAJ = 1 , 2
          DO JL = KIDIA,KFDIA
             PRJ(JL,JAJ,KLEV+1) = 1._DP
             PRK(JL,JAJ,KLEV+1) = PREFZ(JL,JAJ,KLEV+1)
          ENDDO

          DO JK = 1 , KLEV
             JKL = KLEV+1 - JK
             JKLP1 = JKL + 1
             DO JL = KIDIA,KFDIA
                ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,JAJ,JKL)
                PRJ(JL,JAJ,JKL) = ZRE11
                PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,JAJ,JKL)
             ENDDO
          ENDDO
       ENDDO

    ENDIF

    !     ------------------------------------------------------------------

    RETURN
  END SUBROUTINE rad_sw_SWCLR
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE rad_sw_SWR &
       &( KIDIA , KFDIA , KBDIM , KLEV  , KNU &
       &, PAOT_GAMMA,PAOT_OMEGA, PAOT_TAU&
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

    IMPLICIT NONE
    INTRINSIC :: EPSILON, EXP, MAX, MIN

    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KLEV
    INTEGER :: KBDIM
    INTEGER :: KNU

    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    REAL(DP) :: PALBD(KBDIM,NSW)      , PCG(KBDIM,NSW,KLEV)&
         &,  PCLD(KBDIM,KLEV)&
         &,  POMEGA(KBDIM,NSW,KLEV)&
         &,  PSEC(KBDIM)           , PTAU(KBDIM,NSW,KLEV)
    REAL(DP) :: PAOT_GAMMA(KBDIM,KLEV,NSW) &
             ,  PAOT_OMEGA(KBDIM,KLEV,NSW), PAOT_TAU(KBDIM,KLEV,NSW)

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
       DO JA = 1 , 6               ! qqq better as parameter?
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
               &/(1._DP-MIN(ZCLOUD(JL),1._DP-EPSILON(1.)))
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
                  &/(1._DP-MIN(ZCLOUD(JL),1._DP-EPSILON(1.)))
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

          !     --------------------------------------------------------------

          !*         3.1  EQUIVALENT ZENITH ANGLE
          !               -----------------------

          ZMUE = (1._DP-ZC1I(JL,JK)) * PSEC(JL)+ ZC1I(JL,JK) * 1.66_DP
          PRMUE(JL,JK) = 1._DP/ZMUE

          !     -------------------------------------------------------------

          !*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
          !               ----------------------------------------------------

          ZGAP = PCGAZ(JL,JKM1)
          ZBMU0 = 0.5_DP - 0.75_DP * ZGAP / ZMUE
          ZWW = PPIZAZ(JL,JKM1)
          ZTO = PTAUAZ(JL,JKM1)
          ZDEN = 1._DP + (1._DP - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE &
               &+ (1._DP-ZWW) * (1._DP - ZWW +2._DP*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE
          PRAY1(JL,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN
          PTRA1(JL,JKM1) = 1._DP / ZDEN

          ZMU1 = 0.5_DP
          ZBMU1 = 0.5_DP - 0.75_DP * ZGAP * ZMU1
          ZDEN1= 1._DP + (1._DP - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1 &
               &+ (1._DP-ZWW) * (1._DP - ZWW +2._DP*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1
          PRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1
          PTRA2(JL,JKM1) = 1._DP / ZDEN1

          !     ---------------------------------------------------------------

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

       CALL rad_sw_SWDE ( KIDIA, KFDIA , KBDIM &
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

          ZTR(JL,2,JKM1) = ZRNEB(JL) * ZTR1(JL)+ PTRA1(JL,JKM1) &
               * (1._DP-ZRNEB(JL))

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
  END SUBROUTINE rad_sw_SWR
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE rad_sw_SWNI &
       &( KIDIA , KFDIA , KBDIM  , KLEV , KNU &
       &, PAOT_GAMMA,PAOT_OMEGA, PAOT_TAU &
       &, PAKI  , PALBD , PALBP, PCG , PCLD, PCLEAR &
       &, PDSIG , POMEGA, POZ   , PRMU , PSEC, PTAU &
       &, PUD   , PWV   , PQS &
       &, PFDOWN, PFUP  , PCDOWN, PCUP , PSUDU2 )

    !**** *SWNI* - SHORTWAVE RADIATION, NEAR-INFRARED SPECTRAL INTERVALS

    !     PURPOSE.
    !     --------

    !          COMPUTES THE SHORTWAVE RADIATION FLUXES IN THE NEAR-INFRARED 
    !     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

    !**   INTERFACE.
    !     ----------

    !          *SWNI* IS CALLED FROM *SW*.


    !        IMPLICIT ARGUMENTS :
    !        --------------------

    !     ==== INPUTS ===
    !     ==== OUTPUTS ===

    !     METHOD.
    !     -------

    !          1. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING TO
    !     CONTINUUM SCATTERING
    !          2. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING FOR
    !     A GREY MOLECULAR ABSORPTION
    !          3. LAPLACE TRANSFORM ON THE PREVIOUS TO GET EFFECTIVE AMOUNTS
    !     OF ABSORBERS
    !          4. APPLY H2O AND U.M.G. TRANSMISSION FUNCTIONS
    !          5. MULTIPLY BY OZONE TRANSMISSION FUNCTION

    !     EXTERNALS.
    !     ----------

    !          *SWCLR*, *SWR*, *SWDE*, *SWTT*

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
    !        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
    !        95-12-07   J.-J. MORCRETTE    NEAR-INFRARED SW
    !        990128     JJMorcrette        Sunshine duration
    !        99-05-25   JJMorcrette        Revised aerosols

    !     ------------------------------------------------------------------

    IMPLICIT NONE
    INTRINSIC :: EPSILON, EXP, LOG, MAX, MIN

    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KLEV
    INTEGER :: KBDIM
    INTEGER :: KNU

    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    REAL(DP)::  PAKI(KBDIM,2,NSW)&
         &,  PALBD(KBDIM,NSW)      , PALBP(KBDIM,NSW)&
         &,  PCG(KBDIM,NSW,KLEV)   , PCLD(KBDIM,KLEV)&
         &,  PCLEAR(KBDIM)         , PDSIG(KBDIM,KLEV)&
         &,  POMEGA(KBDIM,NSW,KLEV), POZ(KBDIM,KLEV)&
         &,  PQS(KBDIM,KLEV)&
         &,  PRMU(KBDIM)           , PSEC(KBDIM)&
         &,  PTAU(KBDIM,NSW,KLEV)  , PUD(KBDIM,5,KLEV+1)&
         &,  PWV(KBDIM,KLEV)

    REAL(DP) :: PAOT_GAMMA(KBDIM,KLEV,NSW) &
         , PAOT_OMEGA(KBDIM,KLEV,NSW), PAOT_TAU(KBDIM,KLEV,NSW)

    REAL(DP):: PFDOWN(KBDIM,KLEV+1)  , PFUP(KBDIM,KLEV+1)&
         &,  PCDOWN(KBDIM,KLEV+1)  , PCUP(KBDIM,KLEV+1)&
         &,  PSUDU2(KBDIM)

    !     ------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    INTEGER :: IIND2(2), IIND3(3)
    REAL(DP):: ZCGAZ(KBDIM,KLEV)  , ZDIFF(KBDIM)         , ZDIRF(KBDIM)&
         &,  ZFD(KBDIM,KLEV+1)  , ZFU(KBDIM,KLEV+1) &
         &,  ZG(KBDIM)          , ZGG(KBDIM)&
         &,  ZPIZAZ(KBDIM,KLEV)&
         &,  ZRAYL(KBDIM)       , ZRAY1(KBDIM,KLEV+1)  , ZRAY2(KBDIM,KLEV+1)&
         &,  ZREF(KBDIM)        , ZREFZ(KBDIM,2,KLEV+1)&
         &,  ZRE1(KBDIM)        , ZRE2(KBDIM)&
         &,  ZRJ(KBDIM,6,KLEV+1), ZRJ0(KBDIM,6,KLEV+1)&
         &,  ZRK(KBDIM,6,KLEV+1), ZRK0(KBDIM,6,KLEV+1)&
         &,  ZRL(KBDIM,8)&
         &,  ZRMUE(KBDIM,KLEV+1), ZRMU0(KBDIM,KLEV+1)  , ZRMUZ(KBDIM)&
         &,  ZRNEB(KBDIM)       , ZRUEF(KBDIM,8)       , ZR1(KBDIM) &
         &,  ZR2(KBDIM,2)       , ZR3(KBDIM,3)         , ZR4(KBDIM)&
         &,  ZR21(KBDIM)        , ZR22(KBDIM)&
         &,  ZS(KBDIM)&
         &,  ZTAUAZ(KBDIM,KLEV) , ZTO1(KBDIM)          , ZTR(KBDIM,2,KLEV+1)&
         &,  ZTRA1(KBDIM,KLEV+1), ZTRA2(KBDIM,KLEV+1)&
         &,  ZTRCLD(KBDIM)      , ZTRCLR(KBDIM)&
         &,  ZTR1(KBDIM)        , ZTR2(KBDIM)&
         &,  ZW(KBDIM)          , ZW1(KBDIM)           , ZW2(KBDIM,2)&
         &,  ZW3(KBDIM,3)       , ZW4(KBDIM)           , ZW5(KBDIM)

    !     LOCAL INTEGER SCALARS
    INTEGER :: IABS, IKL, IKM1, JABS, JAJ, JAJP, JK, JKKI,&
         &JKKP4, JKL, JKLP1, JKM1, JL, JN, JN2J, JREF

    !     LOCAL REAL SCALARS
    REAL(DP):: ZAA, ZBB, ZCNEB, ZRE11, ZRKI, ZRMUM1, ZWH2O
    REAL(DP):: ZCHKS, ZCHKG, ZARGJ, ZARGK

    !     ------------------------------------------------------------------

    !*         1.     NEAR-INFRARED SPECTRAL INTERVAL (0.68-4.00 MICRON)
    !                 --------------------------------------------------

    !*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
    !                 -----------------------------------------

    DO JL = KIDIA,KFDIA
       ZRMUM1 = 1._DP - PRMU(JL)
       ZRAYL(JL) =  RRAY(KNU,1) + ZRMUM1   * (RRAY(KNU,2) + ZRMUM1 &
            &* (RRAY(KNU,3) + ZRMUM1   * (RRAY(KNU,4) + ZRMUM1 &
            &* (RRAY(KNU,5) + ZRMUM1   *  RRAY(KNU,6)     ))))
       ZRAYL(JL) =  MAX(ZRAYL(JL),0._DP)
    ENDDO

    !     ------------------------------------------------------------------

    !*         2.    CONTINUUM SCATTERING CALCULATIONS
    !                ---------------------------------

    !*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
    !                --------------------------------

    CALL rad_sw_SWCLR &
         &( KIDIA , KFDIA , KBDIM ,  KLEV , KNU &
         &, PAOT_GAMMA,PAOT_OMEGA,PAOT_TAU &
         &, PALBP , PDSIG , ZRAYL, PSEC &
         &, ZCGAZ , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0 &
         &, ZRK0  , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2, ZTRCLR &
         &)

    !*         2.2   CLOUDY FRACTION OF THE COLUMN
    !                -----------------------------

    CALL rad_sw_SWR &
         &( KIDIA , KFDIA , KBDIM , KLEV  , KNU &
         &, PAOT_GAMMA,PAOT_OMEGA,PAOT_TAU &
         &, PALBD , PCG   , PCLD , POMEGA, PSEC , PTAU &
         &, ZCGAZ , ZPIZAZ, ZRAY1, ZRAY2 , ZREFZ, ZRJ  , ZRK, ZRMUE &
         &, ZTAUAZ, ZTRA1 , ZTRA2, ZTRCLD &
         &)

    !     ------------------------------------------------------------------

    !*         3.    SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
    !                ------------------------------------------------------

    JN = 2

    DO JABS=1,2

       !*         3.1  SURFACE CONDITIONS
       !               ------------------

       DO JL = KIDIA,KFDIA
          ZREFZ(JL,2,1) = PALBD(JL,KNU)
          ZREFZ(JL,1,1) = PALBD(JL,KNU)
       ENDDO

       !*         3.2  INTRODUCING CLOUD EFFECTS
       !               -------------------------

       DO JK = 2 , KLEV+1
          JKM1 = JK - 1
          IKL=KLEV+1-JKM1
          DO JL = KIDIA,KFDIA
             ZRNEB(JL) = PCLD(JL,JKM1)
             IF (JABS == 1 .AND. ZRNEB(JL) > 2._DP*EPSILON(1._dp)) THEN
                ZWH2O=MAX(PWV(JL,IKL),EPSILON(1._DP))
                ZCNEB=MAX(EPSILON(1._DP),MIN(ZRNEB(JL),1._DP-EPSILON(1._DP)))
                ZBB=PUD(JL,JABS,JKM1)*PQS(JL,IKL)/ZWH2O
                ZAA=MAX((PUD(JL,JABS,JKM1)-ZCNEB*ZBB)/(1._DP-ZCNEB),EPSILON(1._DP))
             ELSE
                ZAA=PUD(JL,JABS,JKM1)
                ZBB=ZAA
             ENDIF
             ZRKI = PAKI(JL,JABS,KNU)
             ZCHKS= MIN(200._DP,ZRKI * ZAA * 1.66_DP)
             ZS(JL) = EXP(-ZCHKS)
             ZCHKG= MIN(200._DP,ZRKI * ZAA / ZRMUE(JL,JK))
             ZG(JL) = EXP(-ZCHKG)
             ZTR1(JL) = 0._DP
             ZRE1(JL) = 0._DP
             ZTR2(JL) = 0._DP
             ZRE2(JL) = 0._DP

             ZW(JL)= POMEGA(JL,KNU,JKM1)
             ZTO1(JL) = PTAU(JL,KNU,JKM1) / ZW(JL)&
                  &+ ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1)&
                  &+ ZBB * ZRKI

             ZR21(JL) = PTAU(JL,KNU,JKM1) + ZTAUAZ(JL,JKM1)
             ZR22(JL) = PTAU(JL,KNU,JKM1) / ZR21(JL)
             ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)&
                  &+ (1._DP - ZR22(JL)) * ZCGAZ(JL,JKM1)
             ZW(JL) = ZR21(JL) / ZTO1(JL)
             ZREF(JL) = ZREFZ(JL,1,JKM1)
             ZRMUZ(JL) = ZRMUE(JL,JK)
          ENDDO

          CALL rad_sw_SWDE ( KIDIA, KFDIA, KBDIM &
               &, ZGG  , ZREF , ZRMUZ, ZTO1, ZW &
               &, ZRE1 , ZRE2 , ZTR1 , ZTR2     )

          DO JL = KIDIA,KFDIA

             ZREFZ(JL,2,JK) = (1._DP-ZRNEB(JL)) * (ZRAY1(JL,JKM1)&
                  &+ ZREFZ(JL,2,JKM1) * ZTRA1(JL,JKM1)&
                  &* ZTRA2(JL,JKM1) ) * ZG(JL) * ZS(JL)&
                  &+ ZRNEB(JL) * ZRE1(JL)

             ZTR(JL,2,JKM1)=ZRNEB(JL)*ZTR1(JL)&
                  &+ (ZTRA1(JL,JKM1)) * ZG(JL) * (1._DP-ZRNEB(JL))

             ZREFZ(JL,1,JK)=(1._DP-ZRNEB(JL))*(ZRAY1(JL,JKM1)&
                  &+ZREFZ(JL,1,JKM1)*ZTRA1(JL,JKM1)*ZTRA2(JL,JKM1)&
                  &/(1._DP-ZRAY2(JL,JKM1)*ZREFZ(JL,1,JKM1)))*ZG(JL)*ZS(JL)&
                  &+ ZRNEB(JL) * ZRE2(JL)

             ZTR(JL,1,JKM1)= ZRNEB(JL) * ZTR2(JL)&
                  &+ (ZTRA1(JL,JKM1)/(1._DP-ZRAY2(JL,JKM1)&
                  &* ZREFZ(JL,1,JKM1)))&
                  &* ZG(JL) * (1._DP -ZRNEB(JL))

          ENDDO
       ENDDO

       !*         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
       !               -------------------------------------------------

       DO JREF=1,2

          JN = JN + 1

          DO JL = KIDIA,KFDIA
             ZRJ(JL,JN,KLEV+1) = 1._DP
             ZRK(JL,JN,KLEV+1) = ZREFZ(JL,JREF,KLEV+1)
          ENDDO

          DO JK = 1 , KLEV
             JKL = KLEV+1 - JK
             JKLP1 = JKL + 1
             DO JL = KIDIA,KFDIA
                ZRE11 = ZRJ(JL,JN,JKLP1) * ZTR(JL,JREF,JKL)
                ZRJ(JL,JN,JKL) = ZRE11
                ZRK(JL,JN,JKL) = ZRE11 * ZREFZ(JL,JREF,JKL)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !     ------------------------------------------------------------------

    !*         4.    INVERT GREY AND CONTINUUM FLUXES
    !                --------------------------------

    !*         4.1   UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES
    !                ---------------------------------------------

    DO JK = 1 , KLEV+1
       DO JAJ = 1 , 5 , 2
          JAJP = JAJ + 1
          DO JL = KIDIA,KFDIA
             ZRJ(JL,JAJ,JK)=        ZRJ(JL,JAJ,JK) - ZRJ(JL,JAJP,JK)
             ZRK(JL,JAJ,JK)=        ZRK(JL,JAJ,JK) - ZRK(JL,JAJP,JK)
             ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , EPSILON(1._DP) )
             ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , EPSILON(1._DP) )
          ENDDO
       ENDDO
    ENDDO

    DO JK = 1 , KLEV+1
       DO JAJ = 2 , 6 , 2
          DO JL = KIDIA,KFDIA
             ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , EPSILON(1._DP) )
             ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , EPSILON(1._DP) )
          ENDDO
       ENDDO
    ENDDO

    !*         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
    !                 ---------------------------------------------

    DO JK = 1 , KLEV+1
       JKKI = 1
       DO JAJ = 1 , 2
          IIND2(1)=JAJ
          IIND2(2)=JAJ
          DO JN = 1 , 2
             JN2J = JN + 2 * JAJ
             JKKP4 = JKKI + 4

             !*         4.2.1  EFFECTIVE ABSORBER AMOUNTS
             !                 --------------------------

             DO JL = KIDIA,KFDIA
                ZARGJ     = MAX( ZRJ(JL,JN,JK) / ZRJ(JL,JN2J,JK),1._DP)
                ZARGK     = MAX( ZRK(JL,JN,JK) / ZRK(JL,JN2J,JK),1._DP)
                ZW2(JL,1) = LOG( ZARGJ )/ PAKI(JL,JAJ,KNU)
                ZW2(JL,2) = LOG( ZARGK )/ PAKI(JL,JAJ,KNU)
             ENDDO

             !*         4.2.2  TRANSMISSION FUNCTION
             !                 ---------------------

             CALL rad_sw_SWTT1 ( KIDIA,KFDIA,KBDIM, KNU, 2, IIND2 &
                  &, ZW2 &
                  &, ZR2                              )

             DO JL = KIDIA,KFDIA
                ZRL(JL,JKKI) = ZR2(JL,1)
                ZRUEF(JL,JKKI) = ZW2(JL,1)
                ZRL(JL,JKKP4) = ZR2(JL,2)
                ZRUEF(JL,JKKP4) = ZW2(JL,2)
             ENDDO

             JKKI=JKKI+1
          ENDDO
       ENDDO

       !*         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
       !                 ------------------------------------------------------

       DO JL = KIDIA,KFDIA
          PFDOWN(JL,JK) = ZRJ(JL,1,JK) * ZRL(JL,1) * ZRL(JL,3)&
               &+ ZRJ(JL,2,JK) * ZRL(JL,2) * ZRL(JL,4)
          PFUP(JL,JK)   = ZRK(JL,1,JK) * ZRL(JL,5) * ZRL(JL,7)&
               &+ ZRK(JL,2,JK) * ZRL(JL,6) * ZRL(JL,8)
       ENDDO

    ENDDO

    !     ------------------------------------------------------------------

    !*         5.    MOLECULAR ABSORPTION ON CLEAR-SKY FLUXES
    !                ----------------------------------------

    !*         5.1   DOWNWARD FLUXES
    !                ---------------

    JAJ = 2
    IIND3(1)=1
    IIND3(2)=2
    IIND3(3)=3

    DO JL = KIDIA,KFDIA
       ZW3(JL,1)=0._DP
       ZW3(JL,2)=0._DP
       ZW3(JL,3)=0._DP
       ZW4(JL)  =0._DP
       ZW5(JL)  =0._DP
       ZR4(JL)  =1._DP
       ZFD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1)
    ENDDO
    DO JK = 1 , KLEV
       IKL = KLEV+1-JK
       DO JL = KIDIA,KFDIA
          ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
          ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
          ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
          ZW4(JL)  =ZW4(JL)  +PUD(JL,4,IKL)/ZRMU0(JL,IKL)
          ZW5(JL)  =ZW5(JL)  +PUD(JL,5,IKL)/ZRMU0(JL,IKL)
       ENDDO

       CALL rad_sw_SWTT1 ( KIDIA,KFDIA,KBDIM, KNU, 3, IIND3 &
            &, ZW3 &
            &, ZR3                              )

       DO JL = KIDIA,KFDIA
          ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))
          ZFD(JL,IKL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)* ZRJ0(JL,JAJ,IKL)
       ENDDO
    ENDDO

    DO JL=KIDIA,KFDIA
       ZDIFF(JL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)*ZTRCLD(JL)
       ZDIRF(JL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)*ZTRCLR(JL)
       PSUDU2(JL) = ((1._DP-PCLEAR(JL)) * ZDIFF(JL)&
            &+PCLEAR(JL) * ZDIRF(JL)) * RSUN(KNU)
    ENDDO

    !*         5.2   UPWARD FLUXES
    !                -------------

    DO JL = KIDIA,KFDIA
       ZFU(JL,1) = ZFD(JL,1)*PALBP(JL,KNU)
    ENDDO

    DO JK = 2 , KLEV+1
       IKM1=JK-1
       DO JL = KIDIA,KFDIA
          ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKM1)*1.66_DP
          ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKM1)*1.66_DP
          ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKM1)*1.66_DP
          ZW4(JL)  =ZW4(JL)  +PUD(JL,4,IKM1)*1.66_DP
          ZW5(JL)  =ZW5(JL)  +PUD(JL,5,IKM1)*1.66_DP
       ENDDO

       CALL rad_sw_SWTT1 ( KIDIA,KFDIA,KBDIM, KNU, 3, IIND3 &
            &, ZW3 &
            &, ZR3                              )

       DO JL = KIDIA,KFDIA
          ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))
          ZFU(JL,JK) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)* ZRK0(JL,JAJ,JK)
       ENDDO
    ENDDO

    !     ------------------------------------------------------------------

    !*         6.     INTRODUCTION OF OZONE AND H2O CONTINUUM ABSORPTION
    !                 --------------------------------------------------

    IABS=3

    !*         6.1    DOWNWARD FLUXES
    !                 ---------------

    DO JL = KIDIA,KFDIA
       ZW1(JL)=0._DP
       ZW4(JL)=0._DP
       ZW5(JL)=0._DP
       ZR1(JL)=0._DP
       PFDOWN(JL,KLEV+1) = ((1._DP-PCLEAR(JL))*PFDOWN(JL,KLEV+1)&
            &+ PCLEAR(JL) * ZFD(JL,KLEV+1)) * RSUN(KNU)
       PCDOWN(JL,KLEV+1) = ZFD(JL,KLEV+1) * RSUN(KNU)
    ENDDO

    DO JK = 1 , KLEV
       IKL=KLEV+1-JK
       DO JL = KIDIA,KFDIA
          ZW1(JL) = ZW1(JL)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
          ZW4(JL) = ZW4(JL)+PUD(JL,4,IKL)/ZRMUE(JL,IKL)
          ZW5(JL) = ZW5(JL)+PUD(JL,5,IKL)/ZRMUE(JL,IKL)
          ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))
       ENDDO

       CALL rad_sw_SWTT ( KIDIA,KFDIA,KBDIM, KNU, IABS, ZW1, ZR1 )

       DO JL = KIDIA,KFDIA
          PFDOWN(JL,IKL) = ((1._DP-PCLEAR(JL))*ZR1(JL)*ZR4(JL)*PFDOWN(JL,&
               &IKL)&
               &+PCLEAR(JL)*ZFD(JL,IKL)) * RSUN(KNU)
          PCDOWN(JL,IKL) = ZFD(JL,IKL) * RSUN(KNU)
       ENDDO
    ENDDO

    !*         6.2    UPWARD FLUXES
    !                 -------------

    DO JL = KIDIA,KFDIA
       PFUP(JL,1) = ((1._DP-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,1)&
            &+PCLEAR(JL)*ZFU(JL,1)) * RSUN(KNU)
       PCUP(JL,1) = ZFU(JL,1) * RSUN(KNU)
    ENDDO

    DO JK = 2 , KLEV+1
       IKM1=JK-1
       DO JL = KIDIA,KFDIA
          ZW1(JL) = ZW1(JL)+POZ(JL  ,IKM1)*1.66_DP
          ZW4(JL) = ZW4(JL)+PUD(JL,4,IKM1)*1.66_DP
          ZW5(JL) = ZW5(JL)+PUD(JL,5,IKM1)*1.66_DP
          ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))
       ENDDO

       CALL rad_sw_SWTT ( KIDIA,KFDIA,KBDIM, KNU, IABS, ZW1, ZR1 )

       DO JL = KIDIA,KFDIA
          PFUP(JL,JK) = ((1._DP-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,JK)&
               &+PCLEAR(JL)*ZFU(JL,JK)) * RSUN(KNU)
          PCUP(JL,JK) = ZFU(JL,JK) * RSUN(KNU)
       ENDDO
    ENDDO

    !     ------------------------------------------------------------------

    RETURN
  END SUBROUTINE rad_sw_SWNI
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
!OPTIONS XOPT(HSFUN)
  SUBROUTINE rad_sw_SWDE &
       &( KIDIA, KFDIA, KBDIM &
       &, PGG  , PREF , PRMUZ, PTO1, PW &
       &, PRE1 , PRE2 , PTR1 , PTR2 &
       &)

    !**** *SWDE* - DELTA-EDDINGTON IN A CLOUDY LAYER

    !     PURPOSE.
    !     --------
    !           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLOUDY
    !     LAYER USING THE DELTA-EDDINGTON'S APPROXIMATION.

    !**   INTERFACE.
    !     ----------
    !          *SWDE* IS CALLED BY *SWR*, *SWNI*


    !        EXPLICIT ARGUMENTS :
    !        --------------------
    ! PGG    : (KBDIM)             ; ASSYMETRY FACTOR
    ! PREF   : (KBDIM)             ; REFLECTIVITY OF THE UNDERLYING LAYER
    ! PRMUZ  : (KBDIM)             ; COSINE OF SOLAR ZENITH ANGLE
    ! PTO1   : (KBDIM)             ; OPTICAL THICKNESS
    ! PW     : (KBDIM)             ; SINGLE SCATTERING ALBEDO
    !     ==== OUTPUTS ===
    ! PRE1   : (KBDIM)             ; LAYER REFLECTIVITY ASSUMING NO
    !                             ; REFLECTION FROM UNDERLYING LAYER
    ! PTR1   : (KBDIM)             ; LAYER TRANSMISSIVITY ASSUMING NO
    !                             ; REFLECTION FROM UNDERLYING LAYER
    ! PRE2   : (KBDIM)             ; LAYER REFLECTIVITY ASSUMING
    !                             ; REFLECTION FROM UNDERLYING LAYER
    ! PTR2   : (KBDIM)             ; LAYER TRANSMISSIVITY ASSUMING
    !                             ; REFLECTION FROM UNDERLYING LAYER

    !        IMPLICIT ARGUMENTS :   NONE
    !        --------------------

    !     METHOD.
    !     -------

    !          STANDARD DELTA-EDDINGTON LAYER CALCULATIONS.

    !     EXTERNALS.
    !     ----------

    !          NONE

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
    !        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 88-12-15
    !                   96-05-30 Michel Deque (security in EXP()) 

    !     ------------------------------------------------------------------


    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    IMPLICIT NONE
    INTRINSIC :: EXP, MAX, MIN, SQRT

    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KBDIM

    REAL(DP):: PGG(KBDIM),PREF(KBDIM),PRMUZ(KBDIM),PTO1(KBDIM),PW(KBDIM)
    REAL(DP):: PRE1(KBDIM),PRE2(KBDIM),PTR1(KBDIM),PTR2(KBDIM)

    !     LOCAL INTEGER SCALARS
    INTEGER :: JL

    !     LOCAL REAL SCALARS
    REAL(DP):: ZA11, ZA12, ZA13, ZA21, ZA22, ZA23, ZALPHA,&
         &ZAM2B, ZAP2B, ZARG, ZARG2, ZB21, ZB22, ZB23, &
         &ZBETA, ZC1A, ZC1B, ZC2A, ZC2B, ZDENA, ZDENB, &
         &ZDT, ZEXKM, ZEXKP, ZEXMU0, ZFF, ZGP, ZRI0A, &
         &ZRI0B, ZRI0C, ZRI0D, ZRI1A, ZRI1B, ZRI1C, &
         &ZRI1D, ZRK, ZRM2, ZRP, ZTOP, ZWCP, ZWM, ZX1, &
         &ZX2, ZXM2P, ZXP2P

    !     ------------------------------------------------------------------

    !*         1.      DELTA-EDDINGTON CALCULATIONS

    DO JL   =   KIDIA,KFDIA

       !*         1.1     SET UP THE DELTA-MODIFIED PARAMETERS

       ZFF = PGG(JL)*PGG(JL)
       ZGP = PGG(JL)/(1._DP+PGG(JL))
       ZTOP = (1._DP- PW(JL) * ZFF) * PTO1(JL)
       ZWCP = (1._DP-ZFF)* PW(JL) /(1._DP- PW(JL) * ZFF)
       ZDT = 2._DP/3._DP
       ZX1 = 1._DP-ZWCP*ZGP
       ZWM = 1._DP-ZWCP
       ZRM2 =  PRMUZ(JL) * PRMUZ(JL)
       ZRK = SQRT(3._DP*ZWM*ZX1)
       ZX2 = 4._DP*(1._DP-ZRK*ZRK*ZRM2)
       ZRP=ZRK/ZX1
       ZALPHA = 3._DP*ZWCP*ZRM2*(1._DP+ZGP*ZWM)/ZX2
       ZBETA = 3._DP*ZWCP* PRMUZ(JL) *(1._DP+3._DP*ZGP*ZRM2*ZWM)/ZX2
       ZARG=MAX(-200._DP,MIN(ZTOP/PRMUZ(JL),200._DP))
       ZEXMU0=EXP(-ZARG)
       ZARG2=MIN(ZRK*ZTOP,200._DP)
       ZEXKP=EXP(ZARG2)
       ZEXKM = 1._DP/ZEXKP
       ZXP2P = 1._DP+ZDT*ZRP
       ZXM2P = 1._DP-ZDT*ZRP
       ZAP2B = ZALPHA+ZDT*ZBETA
       ZAM2B = ZALPHA-ZDT*ZBETA

       !*         1.2     WITHOUT REFLECTION FROM THE UNDERLYING LAYER

       ZA11 = ZXP2P
       ZA12 = ZXM2P
       ZA13 = ZAP2B
       ZA22 = ZXP2P*ZEXKP
       ZA21 = ZXM2P*ZEXKM
       ZA23 = ZAM2B*ZEXMU0
       ZDENA = ZA11 * ZA22 - ZA21 * ZA12
       ZC1A = (ZA22*ZA13-ZA12*ZA23)/ZDENA
       ZC2A = (ZA11*ZA23-ZA21*ZA13)/ZDENA
       ZRI0A = ZC1A+ZC2A-ZALPHA
       ZRI1A = ZRP*(ZC1A-ZC2A)-ZBETA
       PRE1(JL) = (ZRI0A-ZDT*ZRI1A)/ PRMUZ(JL)
       ZRI0B = ZC1A*ZEXKM+ZC2A*ZEXKP-ZALPHA*ZEXMU0
       ZRI1B = ZRP*(ZC1A*ZEXKM-ZC2A*ZEXKP)-ZBETA*ZEXMU0
       PTR1(JL) = ZEXMU0+(ZRI0B+ZDT*ZRI1B)/ PRMUZ(JL)

       !*         1.3     WITH REFLECTION FROM THE UNDERLYING LAYER

       ZB21 = ZA21- PREF(JL) *ZXP2P*ZEXKM
       ZB22 = ZA22- PREF(JL) *ZXM2P*ZEXKP
       ZB23 = ZA23- PREF(JL) *ZEXMU0*(ZAP2B - PRMUZ(JL) )
       ZDENB = ZA11 * ZB22 - ZB21 * ZA12
       ZC1B = (ZB22*ZA13-ZA12*ZB23)/ZDENB
       ZC2B = (ZA11*ZB23-ZB21*ZA13)/ZDENB
       ZRI0C = ZC1B+ZC2B-ZALPHA
       ZRI1C = ZRP*(ZC1B-ZC2B)-ZBETA
       PRE2(JL) = (ZRI0C-ZDT*ZRI1C) / PRMUZ(JL)
       ZRI0D = ZC1B*ZEXKM + ZC2B*ZEXKP - ZALPHA*ZEXMU0
       ZRI1D = ZRP * (ZC1B*ZEXKM - ZC2B*ZEXKP) - ZBETA*ZEXMU0
       PTR2(JL) = ZEXMU0 + (ZRI0D + ZDT*ZRI1D) / PRMUZ(JL)

    ENDDO
    RETURN
  END SUBROUTINE rad_sw_SWDE
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE rad_sw_SWTT ( KIDIA, KFDIA, KBDIM, KNU, KA , PU, PTR)

    !**** *SWTT* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

    !     PURPOSE.
    !     --------
    !           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
    !     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
    !     INTERVALS.

    !**   INTERFACE.
    !     ----------
    !          *SWTT* IS CALLED FROM *SW1S*, *SWNI*.


    !        EXPLICIT ARGUMENTS :
    !        --------------------
    ! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
    ! KA     :                     ; INDEX OF THE ABSORBER
    ! PU     : (KBDIM)             ; ABSORBER AMOUNT
    !     ==== OUTPUTS ===
    ! PTR    : (KBDIM)             ; TRANSMISSION FUNCTION

    !        IMPLICIT ARGUMENTS :   NONE
    !        --------------------

    !     METHOD.
    !     -------

    !          TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
    !     AND HORNER'S ALGORITHM.

    !     EXTERNALS.
    !     ----------

    !          NONE

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
    !        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 88-12-15

    !-----------------------------------------------------------------------

    IMPLICIT NONE

    !     DUMMY INTEGER SCALARS
    INTEGER :: KA
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KBDIM
    INTEGER :: KNU

    !-----------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    REAL(DP):: PU(KBDIM), PTR(KBDIM)

    !-----------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    REAL(DP):: ZR1(KBDIM), ZR2(KBDIM)

    !     LOCAL INTEGER SCALARS
    INTEGER :: JL


    !-----------------------------------------------------------------------

    !*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION

    DO JL = KIDIA,KFDIA
       ZR1(JL) = APAD(KNU,KA,1) + PU(JL) * (APAD(KNU,KA,2) + PU(JL)&
            &* ( APAD(KNU,KA,3) + PU(JL) * (APAD(KNU,KA,4) + PU(JL)&
            &* ( APAD(KNU,KA,5) + PU(JL) * (APAD(KNU,KA,6) + PU(JL)&
            &* ( APAD(KNU,KA,7) ))))))

       ZR2(JL) = BPAD(KNU,KA,1) + PU(JL) * (BPAD(KNU,KA,2) + PU(JL)&
            &* ( BPAD(KNU,KA,3) + PU(JL) * (BPAD(KNU,KA,4) + PU(JL)&
            &* ( BPAD(KNU,KA,5) + PU(JL) * (BPAD(KNU,KA,6) + PU(JL)&
            &* ( BPAD(KNU,KA,7) ))))))


       !*         2.      ADD THE BACKGROUND TRANSMISSION

       PTR(JL) = (ZR1(JL) / ZR2(JL)) * (1._DP - D(KNU,KA)) + D(KNU,KA)
    ENDDO

    RETURN
  END SUBROUTINE rad_sw_SWTT
  !-----------------------------------------------------------------------

! ***************************************************************************
END MODULE messy_rad_short_v1
! ***************************************************************************
