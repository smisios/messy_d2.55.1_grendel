! mz_ap_20070823 exacly as growth for coupled system in MESSy
!--------------------------------------------------
SUBROUTINE GROWTH_coupled(HO,AO,HSNO,HN,AN,HSNN,                          &
     RPRECW,RPRECI,AOFLDHW,AOFLNHW,AOFLCHI,AOFLRHI,               &
     QS,QT,DZ,Z,DT,FW,OM,                                         &
     PRECN,PRECH,                                                 &
     AOFLSHW,QNSW,QNLW,QNSE,QNLA,HEATABS,SWSUM)
  !=======================================================================

  !
  !     PROGRAMMED BY:
  !     --------------
  !     B.OWENS, P.LEMKE
  !
  !     MODIFIED BY:
  !     ------------
  !     A. STOESSEL                 MPI,  HAMBURG                     1989
  !     S. LEGUTKE                  DKRZ, HAMBURG               24.12.1995
  !
  !        - ONLY SURFACE RESIDUAL HEAT USED FOR SNOW MELT
  !
  !        - SNOW --> ICE CONVERSION BECOMES SNOW --> WATER CONVERSION
  !          IF ICE THICKNESS > HSNOTOICE 
  !
  !     S. LEGUTKE                  DKRZ, HAMBURG               24.06.1997
  !        - ADJUST COMPACTNESS IN ORDER TO HAVE MINIMUM ICE THICKNESS SICTHMIN
  !                
  !        - COMPILE OPTION FOR FORCING WITH FLUXES INSTEAD OF SURFACE FIELDS 
  !          (FLUXES, USED WITH THE COUPLED MODEL)
  !
  !     UWE MIKOLAJEWICZ     5/99
  !        -INCLUDE EVAPORATION (FROM LATENT HEAT) IN FORING WITH SURFACE FIELDS
  !
  !     S. VENZKE                    MPI, HAMBURG               31.08.1999
  !         -ADJUST COMPILE OPTION FOR FORCING WITH FLUXES INSTEAD OF SURFACE 
  !          FIELDS
  !
  !     UWE                                                     28.06.00
  !         INCLUDE SW-PENETRATION UNDER SEA ICE
  !
  !     S. LEGUTKE                  DKRZ, HAMBURG               04.08.2001
  !        - At exit PRECN [m/s] contains all atmospheric freshwater fluxes
  !          relevant for tracer concentration (i.e. snowmelt, river runoff,
  !          P-E over open water, but no snow sublimation). PRECN also
  !          contains snow->water conversion for cases of very thick ice.
  !        - At exit BRINE [m/s] contains all changes of mean ice thickness
  !          (eq. water depth) relevant for tracer concentration.
  !        - BRINE(,) -> COMMON COMMOAU3 for use in marine bgc tracer dilution.
  !
  !      H.HAAK  MPI, HAMBURG       31.11.05
  !
  !          CORE Forcing added
  !
  !
  !     PURPOSE:
  !     --------
  !     UPDATE OF ICE THICKNESS, COMPACTNESS, SNOW DEPTH, UPPER OCEAN
  !     TEMPERATURE AND SALINITY DUE TO ATMOSPHERIC HEAT AND FRESHWATER FLUX.
  !
  !     METHOD:
  !     -------
  !     GROWTH RATES OF THIN/THICK ICE AND SNOW DEPTH DUE TO ATMOSPHERIC 
  !     FLUXES OF HEAT AND FRESHWATER ARE COMPUTED BASED ON THICKNESS
  !     AND COMPACTNESS BEFORE THE UPDATE BY ADVECTIVE CHANGES IN SBR OCICE,
  !     I.E. H(,,1) AND A(,,1) (I.E. HO AND AO) ARE USED IN BUDGET.
  !     THESE GROWTH RATES ARE THEN USED TO UPDATE THICKNESS AND COMPACTNESS
  !     AS OBTAINED IN SBR OCICE, I.E. H(,,2) AND A(,,2) (HN AND AN)
  !     ARE UPDATED IN GROWTH.
  !     NOTE: SURFACE ELEVATION Z CONTAINS ICE AND SNOW LAYER FOR DYNAMICS,
  !     THUS ALL PRECIPITATION HAS TO BE ADDED TO Z, EVEN SNOWFALL.
  !     FOR MIXING PROCESSES ICE+SNOW DRAFT HAS TO BE SUBTRACTED.
  !
  !SV 21.04.99 INCLUDED OPTION FOR FORCING WITH FLUXES
  !C     INTERFACE:
  !     ----------
  !     INPUT:
  !     *HO*       OLD ICE THICKNESS                      [M]
  !     *AO*       OLD ICE COMPACTNESS                    [FRAC.]
  !     *HSNO*     OLD SNOW DEPTH                         [M]
  !     *RPRECW*   (P-E)WATER+RUNOFF+(P-E)ICE:IF RAIN     [M/S,EQ. WATER COLUMN]
  !     *RPRECI*   (P-E)ICE:IF SNOW FALL OR SUBLIMATION   [M/S, " ]
  !     *PA*       ATMOSPHERIC SURFACE PRESURE            [PA]
  !     *AOFLDHW*  HEAT FLUX CORRECTION                   [W/M**2]
  !     *AOFLNHW*  NET HEAT FLUX OVER WATER               [W/M**2]
  !     *AOFLCHI*  CONDUCTIVE HEAT FLUX > 0 IF UPWARD     [W/M**2]
  !     *AOFLRHI*  RESIDUAL HEAT FLUX (FOR SNOW/ICE MELT) [W/M**2]
  !                > 0 IF DOWNWARD
  !SV 26.04.99 READ SOLAR HEAT FLUX FIELD TO ENABLE UWE'S QNET-DIAGNOSTICS
  !     *AOFLSHW*  SOLAR HEAT FLUX                   [W/M**2]
  !
  !     *DZ*       UNDISTURBED 1.LAYER THICK.             [M]
  !     *DT*       TIME STEP                              [S]
  !     *OM*       LAND/SEA MASK                          [0/1]
  !     *SWSUM*    FRACTION OF SW-RADIATION PENETRATION INTO DEEPER LAYERS
  !
  !     OUTPUT:
  !     *FW*       CHANGE OF ICE THICKNESS (MELT > 0)
  !                INCLUDING SNOW-TO-ICE CONVERSION       [M]
  !     *PRECN*    NET FRESHWATER FLUX MOD. BY SNOW CHANGE[M/S]
  !SV     *BRINE*    NET ICE GROWTH                         [M/S]
  !SV 31.08.99 INCLUDED HEAT FLUX DUMMY FIELDS TO ENABLE UWE'S QNET-DIAGNOSTICS
  !     *QNSW*     ABSORBED SOLAR RADIATION          [W/M**2]
  !     *QNLW*     OUTGOING LONGWAVE HEAT FLUX       [W/M**2]
  !     *QNSE*     SENSIBLE HEAT FLUX                [W/M**2]
  !     *QNLA*     LATENT HEAT FLUX                  [W/M**2]
  !     *HEATABS*  SW-RADIATION PENETRATING INTO DEEPER LAYERS [W/M**2]
  !SVX 31.08.99
  ! 
  !     CHANGED:
  !     *Z*        SEA SURFACE ELEVATION                  [M]
  !     *HN*       NEW ICE THICKNESS                [M]
  !     *AN*       NEW ICE COMPACTNESS              [FRAC.]
  !     *HSNN*     NEW SNOW DEPTH                   [M]
  !     *QS*       OCEANIC UPPER LAYER SALINITY     [PSU]
  !     *QT*       OCEANIC UPPER LAYER TEMP.        [DEG C]
  !
  !SVX 31.08.99
  !
  !     EXTERNALS:
  !     ----------
  !     OBUDGET: CALCULATES OPEN WATER ICE GROWTH RATES
  !      BUDGET: CALCULATES ICE GROWTH RATES OVER ICE
  !
  !=======================================================================
  !

  USE MO_PARAM1
  USE MO_COMMOAU1 
  USE MO_COMMOAU3 
  USE MO_UNITS

  USE MO_COMMO1, ONLY: tice,sicuo,sicve,sictho,sicsno,   &
       uko,vke,tho,weto,sao,dlxp,dlyp,                   &
       tafo,txo,tye,fprec,fswr,ftdew,fu10,fclou,giriv,fslp &
       ,lhfldiag

#ifdef CORE
  USE MO_NCAR_OCEAN_FLUXES
  USE MO_PARALLEL
#else
  USE MO_OMIP
#endif


!#ifdef TESTOUT_HFL
  USE MO_MEAN, only : dqswo,dqlwo,dqseo,dqlao,dqtho   &
                     ,dqswi,dqlwi,dqsei,dqlai,dqthi,dticeo
!#endif
  !SV 31.08.99 INCLUDED OPTION FOR FORCING WITH FLUXES

  DIMENSION HO(IE,JE),AO(IE,JE),HSNO(IE,JE)
  DIMENSION RPRECW(IE,JE),RPRECI(IE,JE),Z(IE,JE),OM(IE,JE)
  DIMENSION AOFLNHW(IE,JE),AOFLDHW(IE,JE),AOFLRHI(IE,JE)
  DIMENSION AOFLCHI(IE,JE)
  !SV 31.08.99 INCLUDE SOLAR HEAT FLUX FIELD TO ENABLE UWE'S QNET-DIAGNOSTICS
  DIMENSION AOFLSHW(IE,JE)
  !
  DIMENSION FW(IE,JE),PRECN(IE,JE),PRECH(IE,JE)
  !
  DIMENSION HN(IE,JE),AN(IE,JE),HSNN(IE,JE)
  DIMENSION QS(IE,JE),QT(IE,JE)
  !SV 31.08.99 INCLUDED HEAT FLUX DUMMY FIELDS TO ENABLE UWE'S QNET-DIAGNOSTICS
  DIMENSION QNSW(IE,JE),QNLW(IE,JE),QNSE(IE,JE),QNLA(IE,JE)
  !
  !  LOCAL VARIABLES
  DIMENSION QTM(IE,JE),SH(IE,JE)
  DIMENSION TICM(IE,JE), TICA(IE,JE)
  DIMENSION FO(IE,JE)
  DIMENSION QH(IE,JE),HS(IE,JE),HT(IE,JE)
  DIMENSION HDRAFT(IE,JE)
  DIMENSION RA(IE,JE),RHS(IE,JE),RHB(IE,JE),RH(IE,JE)
  DIMENSION RHSA(IE,JE),RHBA(IE,JE)
  DIMENSION QHST(IE,JE),SN(IE,JE),QFM(IE,JE)
  DIMENSION H(IE,JE,2),A(IE,JE,2),HSN(IE,JE,2)
  DIMENSION TMP(IE,JE),TMP2(IE,JE),TMP3(IE,JE),TMP4(IE,JE)
  DIMENSION HEATABS(IE,JE)                          

  DIMENSION SWSUM(IE,JE)
   
!UWE     ADDITIONAL LOCAL FIELD  
  REAL :: TMPUWE(IE,JE),TMPUW2(IE,JE),TMPUW3(IE,JE),TMPUW4(IE,JE)

#ifdef CORE
  RPRECW(:,:)=0.0
  RPRECI(:,:)=0.0
  QPRE(:,:)=0.0
  UO(:,:)=0.0
  VO(:,:)=0.0
  TO(:,:)=0.0
  TXOC(:,:)=0.0
  TYOC(:,:)=0.0
  TXIC(:,:)=0.0
  TYIC(:,:)=0.0
  STP(:,:)=0.0
  STPP(:,:)=0.0
  FP(:,:)=0.0
  FPP(:,:)=0.0
  fdiff=0.
#endif

  L = IE
  M = JE
  
  !-----------------------------------------------------------------------
  !  OLD AND NEW TIME LEVEL INDEX
  !-----------------------------------------------------------------------
  
  LRHS=1
  LNEW=2
  
  !-----------------------------------------------------------------------
  !  DETERMINE ICE DRAFT AND TOTAL LAYER THICKNESS:
  !-----------------------------------------------------------------------
  
  DO J=1,JE
     DO I=1,IE

        FO(I,J)=0.0  

        H  (I,J,LNEW) = HN  (I,J)
        A  (I,J,LNEW) = AN  (I,J)
        HSN(I,J,LNEW) = HSNN(I,J)

        FW (I,J)      = HN(I,J)

        H  (I,J,LRHS) = HO  (I,J)
        A  (I,J,LRHS) = AO  (I,J)
        HSN(I,J,LRHS) = HSNO(I,J)

        HDRAFT(I,J)=(RHOSNO*HSN(I,J,LNEW) + RHOICE*H(I,J,LNEW))/RHOWAT

        QH(I,J)= DZ + Z(I,J) - HDRAFT(I,J)

     END DO
  END DO

  !-----------------------------------------------------------------------
  !  NEW-ICE GROWTH RATE:
  !-----------------------------------------------------------------------

  DO J=1,M
     DO I=1,L


        FO(I,J) = -(AOFLNHW(I,J)+AOFLDHW(I,J))/CLB                     &
             +SWSUM(I,J)*AOFLSHW(I,J)/CLB

        !SV 31.08.99 SET HEAT FLUX DUMMY FIELDS FOR UWE'S DIAGNOSTICS ONLY
        !SV   (QNSE IS SET TO QNET-QSOLAR AND QNSW IS SET TO QSOLAR; OTHERS ARE SET TO ZERO)
        QNSE(I,J) = (AOFLNHW(I,J)-AOFLSHW(I,J))*(1.-A(I,J,LRHS))
        QNSW(I,J) = AOFLSHW(I,J)*(1.-A(I,J,LRHS)) 
        QNLW(I,J) = 0.0
        QNLA(I,J) = 0.0

        HEATABS(I,J)=SWSUM(I,J)*AOFLSHW(I,J)*(1.-A(I,J,LRHS))

        
        !-----------------------------------------------------------------------
        !  THICK ICE GROWTH RATES.
        !-----------------------------------------------------------------------
        
        RHS(I,J) = -AOFLRHI(I,J) / CLB
        RHB(I,J) = -(AOFLCHI(I,J)+AOFLDHW(I,J)) / CLB
     END DO
  END DO
  
  
  
  DO J=1,M
     DO I=1,L
        QNSE(I,J) = QNSE(I,J)+(AOFLRHI(I,J)+AOFLCHI(I,J))*A(I,J,LRHS)
     ENDDO
  ENDDO

  
  DO J=1,M
     DO I=1,L

        RA (I,J)  = FO (I,J)*DT
        RHS(I,J)  = RHS(I,J)*DT
        RHB(I,J)  = RHB(I,J)*DT
        RH(I,J)  = RHS(I,J)+RHB(I,J)


     ENDDO
  ENDDO


  DO J=1,M
     DO I=1,L
        

        SH(I,J)       = RH(I,J)*A(I,J,LRHS) + RA(I,J)*(1.-A(I,J,LRHS))
        !
        PRECN(I,J)    = RPRECW(I,J)*DT
        HSN(I,J,LNEW) = HSN(I,J,LNEW)+RPRECI(I,J)                      &
             *DT*RHOWAT/RHOSNO
        
        Z(I,J)        = Z(I,J) + (RPRECW(I,J)+RPRECI(I,J))*DT
        
        !UWE INCLUDE COMPACTNESS   22.12.99
        TMP(I,J)      = RHS(I,J)*A(I,J,LRHS)

        !-----------------------------------------------------------------------
        !  MAKE SURE WE DO NOT END UP WITH NEGATIVE SNOW THICKNESS THROUGH 
        !  EVAPORATION. IN THAT CASE, IN ORDER CONSERE MASS, EVAPORATE 
        !  ICE. THIS DOES NOT CONSERVE HEAT, HOWEVER.
        !-----------------------------------------------------------------------
        !
        TMP2(I,J)=MIN(HSN(I,J,LNEW),0.)
        HSN(I,J,LNEW)=MAX(HSN(I,J,LNEW),0.)
        H(I,J,LNEW) = H(I,J,LNEW) + TMP2(I,J)*RHOSNO/RHOICE
        !

        !-----------------------------------------------------------------------
        !  AT MELTING CONDITIONS AT THE SURFACE (RHS<0), FIRST MELT SNOW 
        !  MAKE SURE WE DO NOT END UP WITH NEGATIVE SNOW THICKNESS.
        !-----------------------------------------------------------------------
        TMP2(I,J)=HSN(I,J,LNEW)+MIN(TMP(I,J),0.)*RHOICE/RHOSNO
        !
        !-----------------------------------------------------------------------
        !  MODIFY ICE GROWTH TO ACCOUNT FOR HEAT LOSS BY SNOW MELT ( WITH CLO).
        !  MODIFY PRECIPITATION TO ACCOUNT FOR SNOW MELT.
        !  COMPUTE ICE MASS + HEAT STORAGE + ATMOS.hEAT FLUX (=-TOTAL HEAT) IN 
        !  UPPER LAYER.
        !  NOTE: IT IS ASSUMED THAT ALL PRECIPITATED WATER (P-E AND SNOW MELT)
        !  HAS 2M-AIR TEMPERATURE (COMPILE OPTION FLUXES ONLY, SINCE IN THE
        !  COUPLED VERSION WE DO NOT WANT TO TRANSFER 2M TEMPERATURE).
        !-----------------------------------------------------------------------
        SN(I,J)=MAX(TMP2(I,J),0.)
        TMP2(I,J)    = HSN(I,J,LNEW)-SN(I,J)
        HSN(I,J,LNEW)= SN(I,J)

        RH(I,J)      = SH(I,J)   +TMP2(I,J)*RHOSNO/RHOICE

        PRECN(I,J)   = PRECN(I,J)+TMP2(I,J)*RHOSNO/RHOWAT
        QHST(I,J)    =  H(I,J,LNEW)                                    &
             -( (QT  (I,J)-TFREZ)*QH(I,J)                     &
             ) * CC/CLB*OM(I,J)                             &
             + RH(I,J)

        !-----------------------------------------------------------------------
        !  WHEN NO ICE IS LEFT (QHST<0) MELT SNOW FIRST.
        !-----------------------------------------------------------------------
        TMP(I,J)=MAX(QHST(I,J),0.)
        TMP3(I,J)=MIN(QHST(I,J),0.)
        SN(I,J)=SN(I,J)+TMP3(I,J)*RHOICE/RHOSNO
        !-----------------------------------------------------------------------
        !  UPDATE FRESH WATER FLUX INTO UPPER OCEAN LAYER (PRECN) WITH 
        !         ADDITIONAL SNOW MELT.
        !  UPDATE HEAT STORAGE IN UPPER OCEAN (QHST) DUE TO ADDITIONAL SNOW MELT.
        !
        !  UPDATE SNOW DEPTH WITH ADDITIONAL SNOW MELT (NOTE: IF SNOW REMAINS
        !         WITHOUT ICE BENEATH, IT WILL BE TRANSFORMED INTO ICE (LOOP 141).
        !UWE REMOVED C  QTM IS HEAT CONTENT IN UPPER OCEAN LAYER.
        ! 
        !  H-QFM IS ICE THICKNESS CHANGE: > 0 ==> MELTING; < 0 ==> FREEZING.
        !
        !  UPDATE OCEANIC LAYER THICKNESS DUE TO ADDITIONAL SNOW MELT.
        !-----------------------------------------------------------------------

        TMP2(I,J)=MIN(SN(I,J),0.)
        SN(I,J)=MAX(SN(I,J),0.)
        PRECN(I,J)= PRECN(I,J)+(HSN(I,J,LNEW)-SN(I,J))*RHOSNO/RHOWAT

     enddo
  enddo

  DO J=1,M
     DO I=1,L
        

        QHST(I,J) = QHST(I,J)+(HSN(I,J,LNEW)-SN(I,J))*RHOSNO/RHOICE


        QFM(I,J)      = H(I,J,LNEW)
        HSN(I,J,LNEW) = SN(I,J)
        RH(I,J)       =-RH(I,J)
        
        !-----------------------------------------------------------------------
        !  UPDATE ICE THICKNESS (EQ.9 IN OWENS & LEMKE 90)
        !-----------------------------------------------------------------------
        !
        !  UPDATE SALT CONTENT AND
        !         HEAT CONTENT OF UPPER OCEANIC LAYER.
        !-----------------------------------------------------------------------
        H(I,J,LNEW)=MAX(QHST(I,J),0.)
        BRINE(I,J) = (QFM(I,J)-H(I,J,LNEW))*OM(I,J)
        !
        HS(I,J)  = QS(I,J)*QH(I,J)                                     &
             +BRINE(I,J)*RHOICWA*MIN(SICE,qs(i,j))
        !
        HT(I,J)  = -TMP2(I,J)*RHOSNO*CLB/(CC*RHOICE)
        !
        !-----------------------------------------------------------------------
        !  UPDATE POTENTIAL TEMPERATURE AND SALINITY OF UPPER OCEANIC LAYER
        !         AND LAYER THICKNESS FOR THERMODYNAMICS.
        !-----------------------------------------------------------------------
        HDRAFT(I,J) =   RHOSNWA*HSN(I,J,LNEW)                          &
             +RHOICWA*H(I,J,LNEW) 
        !
        QH(I,J) = DZ + Z(I,J) - HDRAFT(I,J)
        !
        QT(I,J) = (HT(I,J)/QH(I,J) + TFREZ)*    OM(I,J)                &
             + QT(I,J)                 *(1.-OM(I,J))
        !
        QS(I,J) = HS(I,J)/QH(I,J)
        !
        !-----------------------------------------------------------------------
        !  MAKE SURE WE DON'T TRY TO MELT MORE ICE THAN IS AVAILABLE:
        !  NOTE: POSITIVE RH MEANS MELTING HERE.
        !-----------------------------------------------------------------------
        !
        RH(I,J)=-MIN(RH(I,J),H(I,J,LNEW))
        TMP3(I,J)=MAX(H(I,J,LRHS),HMIN)
        TMP2(I,J)=MIN(RH(I,J),0.)
        TMP(I,J)=MAX(RA(I,J),0.)
        !-----------------------------------------------------------------------
        !  UPDATE ICE COMPACTNESS (EQ.16 IN HIBLER 79)
        !  MAKE SURE WE DO NOT DIVIDE BY 0 
        !  IF MELTING THICK ICE, THEN EVALUATE THE MELTING TERM: TMP2
        !  IF FREEZING THIN ICE, THEN EVALUATE THE FREEZING TERM: TMP.
        !-----------------------------------------------------------------------
        RA(I,J)     = 0.5*TMP2(I,J) * A(I,J,LRHS)/TMP3(I,J)            &
                                ! MODIFY LEADCLOSING IN CASE OF FREEZING
             +    TMP (I,J) *(1.-A(I,J,LRHS))*5.               &
#ifndef LEADCLOSE
        /(H0+4.*MAX(H0,H(I,J,LRHS)/MAX(A(I,J,LRHS),ARMIN)))
#else
        /(H0+3.*MAX(H0,H(I,J,LRHS)/MAX(A(I,J,LRHS),ARMIN)))
#endif 
        !
        !
        A(I,J,LNEW) = A(I,J,LNEW) + RA(I,J)
        !-----------------------------------------------------------------------
        !  ENSURE THAT COMPACTNESS > 0 WHERE THERE IS ICE.
        !  SET COMPACTNESS TO 0 WHERE THERE IS NO ICE.
        !  COMPACTNESS IS NOT ALLOWED TO BECOME LARGER THAN ARMAX.
        !  COMPACTNESS IS NOT ALLOWED TO BECOME LESS THAN 0.
        !-----------------------------------------------------------------------
        IF( H(I,J,LNEW) .GT. 0. .AND. A(I,J,LNEW) .LE. 0.) THEN
           A(I,J,LNEW) = H(I,J,LNEW)/SICTHMIN
        ENDIF
        !
        A(I,J,LNEW) = ( A(I,J,LNEW)*(0.5+SIGN(0.5, A(I,J,LNEW)      )) &
             *(0.5-SIGN(0.5, A(I,J,LNEW)-ARMAX)) &
             +      ARMAX*(0.5+SIGN(0.5, A(I,J,LNEW)-ARMAX)))&
             *(0.5-SIGN(0.5,-H(I,J,LNEW)))
     ENDDO
  ENDDO
  
  !-----------------------------------------------------------------------
  !  SNOW TO ICE CONVERSION:
  !-----------------------------------------------------------------------

  IF ( ISNFLG .EQ. 1 ) THEN
     !
     !-----------------------------------------------------------------------
     !     HSNTOICE IS  MAXIMUM ICE THICKNESS FOR WHICH SNOW --> ICE 
     !     CONVERSION WILL BE DONE.
     !-----------------------------------------------------------------------
     !
     DO J = 1,M
        DO I = 1,L
           !
           TMP3  (I,J)=(0.5-SIGN(0.5,H(I,J,LNEW)-HSNTOICE))
           TMP2(I,J)=MIN(HDRAFT(I,J),H(I,J,LNEW))
           TMP(I,J)=MAX(HDRAFT(I,J),H(I,J,LNEW))
        END DO
     END DO
     !-----------------------------------------------------------------------
     !  SNOW TO ICE CONVERSION:
     !            IN CASE THE ICE SURFACE LIES IN A DEPTH DELTH BELOW THE 
     !            WATER LINE BECAUSE OF THE SNOW LOAD , THE ICE THICKNESS 
     !            IS INCREASED BY THIS DEPTH, AND AN EQUIVALENT AMOUNT (IN
     !            IN HEAT) OF SNOW IS MELTED.
     !            IF THERE IS NET SNOW ACCUMULATION AND NOT ENOUGH ICE MELT,
     !            THIS CAN LEAD TO UNLIMITED ICE GROWTH. THUS IN CASE THE
     !            ICE BECOMES TOO THICK LOCALLY, THE SNOW TO ICE CONVERSION
     !            IS TRANSFORMED INTO A SNOW TO FRESHWATER CONVERSION THERE.
     !  CLOSE THE SALINITY BALANCE.
     !            NOTE:IF ICE IS FORMED, THE HEAT IS TAKEN FROM SNOW MELT,
     !            I.E. NO IMPACT ON OCEAN TEMPERATURE; IF NO ICE IS FORMED
     !            THE HEAT BALANCE IS NOT CLOSED.
     !-----------------------------------------------------------------------
     DTI=1./DT

     DO J=1,M
        DO I=1,L
           !
           DELTH         = HDRAFT(I,J)-TMP2(I,J)
           !
           HSN(I,J,LNEW) = HSN(I,J,LNEW)-DELTH*RHOICE/RHOSNO
           !
           H  (I,J,LNEW) =       TMP3(I,J) *TMP(I,J)                   &
                +(1.0-TMP3(I,J))*H(I,J,LNEW)
           !
           HDRAFT(I,J)   = ( RHOSNO*HSN(I,J,LNEW)                      &
                +RHOICE*H(I,J,LNEW)  )/RHOWAT
           !
           QH(I,J)       = DZ + Z(I,J) - HDRAFT(I,J)
           !
           HS(I,J)       = HS(I,J)                                     &
                -DELTH*SICE*RHOICE/RHOWAT*TMP3(I,J)
           !
           QS(I,J)       = MAX(HS(I,J)/QH(I,J),om(i,j))
           !
           !UWE  CONVERT PRECN TO M/S
           !
           !CHANGES ACC. TO S.L.
           PRECH(I,J)=RPRECW(I,J)+RPRECI(I,J)
           PRECN(I,J) = (PRECN(I,J)+DELTH*RHOICE/RHOWAT)/DT

           BRINE(I,J) = (BRINE(I,J)-DELTH*TMP3(I,J))*RHOICE/RHOWAT/DT

        END DO
     END DO

  ENDIF
  !
  !-----------------------------------------------------------------------
  !  store heat flux (w/m**2) without heat flux correction term in sh.
  !  fw is total thermodynamic change of ice thickness (m of water).
  !  adjust compactness to have minimum ice thickness sicthmin.
  !-----------------------------------------------------------------------
  !
  DO j=1,m
     DO i=1,l
        hn  (i,j) = h  (i,j,lnew)

        an  (i,j) = a(i,j,lnew)                                        &
             *(0.5+SIGN(0.5,h(i,j,lnew)/sicthmin-a(i,j,lnew)))  &
             +h(i,j,lnew)/sicthmin                               &
             *(0.5-SIGN(0.5,h(i,j,lnew)/sicthmin-a(i,j,lnew)))

        hsnn(i,j) = hsn(i,j,lnew)

        fw  (i,j) =( fw(i,j)-h(i,j,lnew) )*rhoice/rhowat

     END DO
  END DO

#ifdef DEBUG

  zzzhh=0.
  zzzth=0.
  zzzsh=0.
  zzzfl=0.
  zzzsw=0.
  zzzse=0.
  zzzla=0.
  zzzlw=0.  

  do i=2,ie-1
    do j=2,je-1
      zzzsw = zzzsw+qnsw(i,j)*dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
      zzzlw = zzzlw+qnlw(i,j)*dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
      zzzse = zzzse+qnse(i,j)*dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
      zzzla = zzzla+qnla(i,j)*dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
     
      zzzhh = zzzhh+prech(i,j)*dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
      zzzth = zzzth+tho(i,j,1)*dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
      zzzsh = zzzsh+sao(i,j,1)*dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
      zzzfl = zzzfl+dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
    END DO
  END DO

  CALL global_sum(zzzhh)
  CALL global_sum(zzzth)
  CALL global_sum(zzzsh)
  CALL global_sum(zzzfl)
  
  CALL global_sum(zzzsw)
  CALL global_sum(zzzlw)
  CALL global_sum(zzzse)
  CALL global_sum(zzzla)

  IF (p_pe==p_io) THEN
    WRITE(IO_STDOUT,*)'global sum p-e+r (Sv): ',(zzzhh)*(1.e-6)
    WRITE(IO_STDOUT,*)'global mean sst (C)   : ',(zzzth/zzzfl)
    WRITE(IO_STDOUT,*)'global mean sss (psu) : ',(zzzsh/zzzfl)
    WRITE(IO_STDOUT,*)'global mean sw (Wm-2) : ',(zzzsw/zzzfl)
    WRITE(IO_STDOUT,*)'global mean lw (Wm-2) : ',(zzzlw/zzzfl)
    WRITE(IO_STDOUT,*)'global mean se (Wm-2) : ',(zzzse/zzzfl)
    WRITE(IO_STDOUT,*)'global mean la (Wm-2) : ',(zzzla/zzzfl)
    WRITE(IO_STDOUT,*)'global mean hf (Wm-2) : ',(zzzsw/zzzfl)+(zzzlw/zzzfl)  &
                                        +(zzzse/zzzfl)+(zzzla/zzzfl)
  ENDIF

#endif
END SUBROUTINE GROWTH_coupled

