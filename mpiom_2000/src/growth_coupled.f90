! mz_ap_20070823 exacly as growth for coupled system in MESSy
!--------------------------------------------------
SUBROUTINE GROWTH_coupled(HO,AO,HSNO,HN,AN,HSNN,                          &
     RPRECW,RPRECI,AOFLNHW,AOFLCHI,AOFLRHI,                       &
     SLN,                                                         &
     QS,QT,DZ,Z,DT,FW,OM,                                         &
     QNSW,QNLW,QNSE,QNLA,PRECN,PRECH,HEATABS,SWSUM)
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
  !     *AOFLNHW*  NET HEAT FLUX OVER WATER               [W/M**2]
  !     *AOFLCHI*  CONDUCTIVE HEAT FLUX > 0 IF UPWARD     [W/M**2]
  !     *AOFLRHI*  RESIDUAL HEAT FLUX (FOR SNOW/ICE MELT) [W/M**2]
  !                > 0 IF DOWNWARD
  !SV 26.04.99 READ SOLAR HEAT FLUX FIELD TO ENABLE UWE'S QNET-DIAGNOSTICS
  !     *SLN*      SOLAR HEAT FLUX /
  !                INCOMING SURFACE SOLAR RAD.            [W/M**2]
  !
  !     *DZ*       UNDISTURBED 1.LAYER THICK.             [M]
  !     *DT*       TIME STEP                              [S]
  !     *OM*       LAND/SEA MASK                          [0/1]
  !     *SWSUM*    FRACTION OF SW-RADIATION PENETRATION INTO DEEPER LAYERS
  !
  !     OUTPUT:
  !     *FW*       CHANGE OF ICE THICKNESS (MELT > 0)
  !                INCLUDING SNOW-TO-ICE CONVERSION       [M]
  !     *PRECN*    NET FRESHWATER FLUX (P-E)
  !                MOD. BY SNOW CHANGE                    [M/S]
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
  !     CHANGED GLOBALLY:
  !     *TICE*     ICE/SNOW SURFACE TEMPERATURE[DEG C]
  !
  !     EXTERNALS:
  !     ----------
  !     OBUDGET: CALCULATES OPEN WATER ICE GROWTH RATES
  !      BUDGET: CALCULATES ICE GROWTH RATES OVER ICE
  !
  !=======================================================================
  !

  USE mo_kind, ONLY: dp, wp
  USE MO_PARAM1, ONLY: ie, je
  USE mo_commo1, ONLY: almzer
  USE MO_COMMOAU1, ONLY: armax, armin, cc, clb, &
       con, consn, h0, hmin, hsntoice, isnflg, sice, sicthmin, subl, &
       tfreeze, vapl
  USE MO_COMMOAU3, ONLY: brine
  USE MO_UNITS, ONLY: io_stdout
  USE MO_COMMO1, ONLY: rleadclose, tice, tho, tafo, fprec, fswr, ftdew, fu10, fclou, &
       lhfldiag !,fslp,giriv
  USE mo_planetary_constants, ONLY : rhoref_water,rhoref_ice,rhoref_snow &
       , rhoicwa, rhosnwa


#ifdef DEBUG
  USE mo_commo1, ONLY: weto, sao, dlxp, dlyp
  USE mo_mpi, ONLY: p_pe, p_io
  USE mo_parallel, ONLY: global_sum
#endif/*def DEBUG */

  USE MO_MEAN, only : dqswo,dqlwo,dqseo,dqlao,dqtho   &
                     ,dqswi,dqlwi,dqsei,dqlai,dqthi,dticeo

  IMPLICIT NONE

  ! in case one wants to use special albedo settings and not those
  ! adjusted to the different forcing sets, consider one of the
  ! following sets:
#ifdef ALBMELTHI
  REAL(wp), PARAMETER :: albsnm = 0.8_wp, albsn = 0.87_wp, albi = 0.78_wp, &
       albm = 0.7_wp
#endif /*ALBMELTHI*/
#ifdef ALBMSN07
  REAL(wp), PARAMETER :: albsn = 0.82_wp, albsnm = 0.7_wp, albm = 0.63_wp
#endif /*ALBMSN07*/
#if 0
  ! old default value:
  REAL, PARAMETER :: albi = 0.75, ALBM   = 0.66
#endif

  !SV 31.08.99 INCLUDED OPTION FOR FORCING WITH FLUXES
  REAL(dp) :: HO(IE,JE),AO(IE,JE),HSNO(IE,JE)
  REAL(dp) :: RPRECW(IE,JE),RPRECI(IE,JE),Z(IE,JE),OM(IE,JE)
  REAL(dp) :: AOFLNHW(IE,JE),AOFLRHI(IE,JE)
  REAL(dp) :: AOFLCHI(IE,JE)
  !SV 31.08.99 INCLUDE SOLAR HEAT FLUX FIELD TO ENABLE UWE'S QNET-DIAGNOSTICS
  REAL(dp) :: SLN(IE,JE)
  !
  REAL(dp) :: FW(IE,JE),PRECN(IE,JE),PRECH(IE,JE)
  !
  REAL(dp) :: HN(IE,JE),AN(IE,JE),HSNN(IE,JE)
  REAL(dp) :: QS(IE,JE),QT(IE,JE)
  !SV 31.08.99 INCLUDED HEAT FLUX DUMMY FIELDS TO ENABLE UWE'S QNET-DIAGNOSTICS
  REAL(dp) :: QNSW(IE,JE),QNLW(IE,JE),QNSE(IE,JE),QNLA(IE,JE)
  !
  !  LOCAL VARIABLES
  REAL(dp) :: QTM(IE,JE),SH(IE,JE)
  REAL(dp) :: TICM(IE,JE), TICA(IE,JE)
  REAL(dp) :: FO(IE,JE)
  REAL(dp) :: QH(IE,JE),HS(IE,JE),HT(IE,JE)
  REAL(dp) :: HDRAFT(IE,JE)
  REAL(dp) :: RA(IE,JE),RHS(IE,JE),RHB(IE,JE),RH(IE,JE)
  REAL(dp) :: RHSA(IE,JE),RHBA(IE,JE)
  REAL(dp) :: QHST(IE,JE),SN(IE,JE),QFM(IE,JE)
  REAL(dp) :: H(IE,JE,2),A(IE,JE,2),HSN(IE,JE,2)
  REAL(dp) :: TMP(IE,JE),TMP2(IE,JE),TMP3(IE,JE),TMP4(IE,JE)
  REAL(dp) :: HEATABS(IE,JE)

  REAL(dp), INTENT(IN) :: dz, dt

  REAL(dp) :: SWSUM(IE,JE)


!UWE     ADDITIONAL LOCAL FIELD
  REAL(dp) :: TMPUWE(IE,JE),TMPUW2(IE,JE),TMPUW3(IE,JE),TMPUW4(IE,JE)
  REAL(dp) :: vapli, anzlevi, subrsni, subrici, rwri, subli, delth
  INTEGER :: l, m, lrhs, lnew, i, j, k, icelev

#ifdef DEBUG

  ! Variables needed for gathering debugging information
   
  ! Global conservation checks
  REAL(wp) :: zzzhh
  REAL(wp) :: zzzth
  REAL(wp) :: zzzsh
  REAL(wp) :: zzzfl
  REAL(wp) :: zzzsw
  REAL(wp) :: zzzse
  REAL(wp) :: zzzla
  REAL(wp) :: zzzlw

#endif/*def DEBUG */

#ifdef CORE
  rprecw(:,:) = 0.0_wp
  rpreci(:,:) = 0.0_wp
  qpre(:,:) = 0.0_wp
  uo(:,:) = 0.0_wp
  vo(:,:) = 0.0_wp
  to(:,:) = 0.0_wp
  txoc(:,:) = 0.0_wp
  tyoc(:,:) = 0.0_wp
  txic(:,:) = 0.0_wp
  tyic(:,:) = 0.0_wp
  stp(:,:) = 0.0_wp
  stpp(:,:) = 0.0_wp
  fp(:,:) = 0.0_wp
  fpp(:,:) = 0.0_wp
  fdiff = 0._wp
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

        fo(i, j) = 0.0_wp


        H  (I,J,LNEW) = HN  (I,J)
        A  (I,J,LNEW) = AN  (I,J)
        HSN(I,J,LNEW) = HSNN(I,J)

        FW (I,J)      = HN(I,J)

        H  (I,J,LRHS) = HO  (I,J)
        A  (I,J,LRHS) = AO  (I,J)
        HSN(I,J,LRHS) = HSNO(I,J)

        HDRAFT(I,J)=(RHOREF_SNOW*HSN(I,J,LNEW) + RHOREF_ICE*H(I,J,LNEW))/RHOREF_WATER

        QH(I,J)= DZ + Z(I,J) - HDRAFT(I,J)

     END DO
  END DO

  !-----------------------------------------------------------------------
  !  NEW-ICE GROWTH RATE:
  !-----------------------------------------------------------------------


  DO J=1,M
     DO I=1,L



        FO(I,J) = -AOFLNHW(I,J)/CLB                     &
             +SWSUM(I,J)*SLN(I,J)/CLB

        !SV 31.08.99 SET HEAT FLUX DUMMY FIELDS FOR UWE'S DIAGNOSTICS ONLY
        !SV   (QNSE IS SET TO QNET-QSOLAR AND QNSW IS SET TO QSOLAR; OTHERS ARE SET TO ZERO)
        QNSE(I,J) = (AOFLNHW(I,J)-SLN(I,J))*(1.-A(I,J,LRHS))
        QNSW(I,J) = SLN(I,J)*(1.-A(I,J,LRHS))
        QNLW(I,J) = 0.0
        QNLA(I,J) = 0.0

        HEATABS(I,J)=SWSUM(I,J)*SLN(I,J)*(1.-A(I,J,LRHS))


        !-----------------------------------------------------------------------
        !  THICK ICE GROWTH RATES.
        !-----------------------------------------------------------------------

        RHS(I,J) = -AOFLRHI(I,J) / CLB
        RHB(I,J) = -AOFLCHI(I,J) / CLB
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


        sh(i, j) = rh(i, j) * a(i, j, lrhs) + ra(i, j) * (1._wp - a(i, j, lrhs))
        !
        PRECN(I,J)    = RPRECW(I,J)*DT
        HSN(I,J,LNEW) = HSN(I,J,LNEW)+RPRECI(I,J)                      &
             *DT*RHOREF_WATER/RHOREF_SNOW

        Z(I,J)        = Z(I,J) + (RPRECW(I,J)+RPRECI(I,J))*DT

        !UWE INCLUDE COMPACTNESS   22.12.99
        TMP(I,J)      = RHS(I,J)*A(I,J,LRHS)

        !-----------------------------------------------------------------------
        !  MAKE SURE WE DO NOT END UP WITH NEGATIVE SNOW THICKNESS THROUGH
        !  EVAPORATION. IN THAT CASE, IN ORDER CONSERE MASS, EVAPORATE
        !  ICE. THIS DOES NOT CONSERVE HEAT, HOWEVER.
        !-----------------------------------------------------------------------
        !
        tmp2(i, j) = MIN(hsn(i, j, lnew), 0._wp)
        hsn(i, j, lnew) = MAX(hsn(i, j, lnew), 0._wp)
        h(i, j, lnew) = h(i, j, lnew) + tmp2(i, j) * (rhoref_snow / rhoref_ice)
        !

        !-----------------------------------------------------------------------
        !  AT MELTING CONDITIONS AT THE SURFACE (RHS<0), FIRST MELT SNOW
        !  MAKE SURE WE DO NOT END UP WITH NEGATIVE SNOW THICKNESS.
        !-----------------------------------------------------------------------
        tmp2(i, j) = hsn(i, j, lnew) + MIN(tmp(i, j), 0._wp) &
             * (rhoref_ice/rhoref_snow)
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
        sn(i, j) = MAX(tmp2(i, j), 0._wp)
        TMP2(I,J)    = HSN(I,J,LNEW)-SN(I,J)
        HSN(I,J,LNEW)= SN(I,J)

        RH(I,J)      = SH(I,J)   +TMP2(I,J)*RHOREF_SNOW/RHOREF_ICE

        PRECN(I,J)   = PRECN(I,J)+TMP2(I,J)*RHOREF_SNOW/RHOREF_WATER

        QHST(I,J)    =  H(I,J,LNEW)                                    &
             - ((qt(i,j) - tfreeze) * qh(i,j)) * cc/clb * om(i, j)     &
             + RH(I,J)

        !-----------------------------------------------------------------------
        !  WHEN NO ICE IS LEFT (QHST<0) MELT SNOW FIRST.
        !-----------------------------------------------------------------------
        tmp(i, j) = MAX(qhst(i, j), 0._wp)
        tmp3(i, j) = MIN(qhst(i, j), 0._wp)
        SN(I,J)=SN(I,J)+TMP3(I,J)*RHOREF_ICE/RHOREF_SNOW
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

        tmp2(i, j) = MIN(sn(i, j), 0._wp)
        sn(i, j) = MAX(sn(i, j), 0._wp)
        PRECN(I,J)= PRECN(I,J)+(HSN(I,J,LNEW)-SN(I,J))*RHOREF_SNOW/RHOREF_WATER

     enddo
  enddo

  DO J=1,M
     DO I=1,L



        QHST(I,J) = QHST(I,J)+(HSN(I,J,LNEW)-SN(I,J))*RHOREF_SNOW/RHOREF_ICE


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
        h(i, j, lnew) = MAX(qhst(i, j), 0._wp)
        BRINE(I,J) = (QFM(I,J)-H(I,J,LNEW))*OM(I,J)
        !
        HS(I,J)  = QS(I,J)*QH(I,J)                                     &
             +BRINE(I,J)*RHOICWA*MIN(SICE,qs(i,j))
        !
        HT(I,J)  = -TMP2(I,J)*RHOREF_SNOW*CLB/(CC*RHOREF_ICE)
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
        qt(i, j) = (ht(i, j)/qh(i, j) + tfreeze) * om(i, j)                &
             + qt(i, j) * (1._wp - om(i, j))
        !
        QS(I,J) = HS(I,J)/QH(I,J)
        !
        !-----------------------------------------------------------------------
        !  MAKE SURE WE DON'T TRY TO MELT MORE ICE THAN IS AVAILABLE:
        !  NOTE: POSITIVE RH MEANS MELTING HERE.
        !-----------------------------------------------------------------------
        !
        rh(i, j) = -MIN(rh(i, j), h(i, j, lnew))
        tmp3(i, j) = MAX(h(i, j, lrhs), hmin)
        tmp2(i, j) = MIN(rh(i, j), 0._wp)
        tmp(i, j) = MAX(ra(i, j), 0._wp)
        !-----------------------------------------------------------------------
        !  UPDATE ICE COMPACTNESS (EQ.16 IN HIBLER 79)
        !  MAKE SURE WE DO NOT DIVIDE BY 0
        !  IF MELTING THICK ICE, THEN EVALUATE THE MELTING TERM: TMP2
        !  IF FREEZING THIN ICE, THEN EVALUATE THE FREEZING TERM: TMP.
        !-----------------------------------------------------------------------

        ! LEADCLOSING parameters (default) can be overwritten in the namelist
        ! rleadclose(1)=0.5
        ! rleadclose(2)=5.
        ! rleadclose(3)=4.
        ! With old LEADCLOSE compile flag, rleadclose(3) would have been 3.

        ra(i,j)= rleadclose(1) * tmp2(i,j) * a(i,j,lrhs) / tmp3(i,j)          &
             ! modify leadclosing in case of freezing
             + tmp(i, j) * (1._wp - a(i,j,lrhs))                              &
             * rleadclose(2)                                                  &
             / (h0+rleadclose(3)*MAX(h0,h(i,j,lrhs)/MAX(a(i,j,lrhs),armin)))

        !
        !
        A(I,J,LNEW) = A(I,J,LNEW) + RA(I,J)
        !-----------------------------------------------------------------------
        !  ENSURE THAT COMPACTNESS > 0 WHERE THERE IS ICE.
        !  SET COMPACTNESS TO 0 WHERE THERE IS NO ICE.
        !  COMPACTNESS IS NOT ALLOWED TO BECOME LARGER THAN ARMAX.
        !  COMPACTNESS IS NOT ALLOWED TO BECOME LESS THAN 0.
        !-----------------------------------------------------------------------
        IF( H(I,J,LNEW) .GT. 0._wp .AND. A(I,J,LNEW) .LE. 0._wp) THEN
           A(I,J,LNEW) = H(I,J,LNEW)/SICTHMIN
        ENDIF
        !
        a(i, j, lnew) = ( a(i, j, lnew) &
             * (0.5_wp + SIGN(0.5_wp, a(i,j,lnew))) &
             * (0.5_wp - SIGN(0.5_wp, a(i,j,lnew) - armax)) &
             +      armax * (0.5_wp + SIGN(0.5_wp, a(i, j, lnew) - armax))) &
             * (0.5_wp - SIGN(0.5_wp, -h(i, j, lnew)))
     END DO
  END DO

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
           tmp3(i, j) = (0.5_wp - SIGN(0.5_wp, h(i, j, lnew) - hsntoice))
           tmp2(i, j) = MIN(hdraft(i, j), h(i, j, lnew))
           tmp(i, j) = MAX(hdraft(i, j), h(i, j, lnew))
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

     DO J=1,M
        DO I=1,L
           !
           DELTH         = HDRAFT(I,J)-TMP2(I,J)
           !
           HSN(I,J,LNEW) = HSN(I,J,LNEW)-DELTH*RHOREF_ICE/RHOREF_SNOW
           !
           H  (I,J,LNEW) =       TMP3(I,J) *TMP(I,J)                   &
                + (1.0_wp - tmp3(i, j)) * h(i, j, lnew)
           !
           HDRAFT(I,J)   = ( RHOREF_SNOW*HSN(I,J,LNEW)                      &
                +RHOREF_ICE*H(I,J,LNEW)  )/RHOREF_WATER
           !
           QH(I,J)       = DZ + Z(I,J) - HDRAFT(I,J)
           !
           HS(I,J)       = HS(I,J)                                     &
                -DELTH*SICE*RHOREF_ICE/RHOREF_WATER*TMP3(I,J)
           !
           QS(I,J)       = MAX(HS(I,J)/QH(I,J),om(i,j))
           !
           !UWE  CONVERT PRECN TO M/S
           !
           !CHANGES ACC. TO S.L.
           PRECH(I,J)=RPRECW(I,J)+RPRECI(I,J)

           PRECN(I,J) = (PRECN(I,J)+DELTH*RHOREF_ICE/RHOREF_WATER)/DT

           BRINE(I,J) = (BRINE(I,J)-DELTH*TMP3(I,J))*RHOREF_ICE/RHOREF_WATER/DT

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
             *(0.5_wp + SIGN(0.5_wp, h(i,j,lnew)/sicthmin-a(i,j,lnew)))  &
             +h(i,j,lnew)/sicthmin                               &
             *(0.5_wp - SIGN(0.5_wp, h(i,j,lnew)/sicthmin-a(i,j,lnew)))

        hsnn(i,j) = hsn(i,j,lnew)

        fw  (i,j) =( fw(i,j)-h(i,j,lnew) )*rhoref_ice/rhoref_water

     END DO
  END DO

#ifdef DEBUG

  zzzhh=0._wp
  zzzth=0._wp
  zzzsh=0._wp
  zzzfl=0._wp
  zzzsw=0._wp
  zzzse=0._wp
  zzzla=0._wp
  zzzlw=0._wp

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
    WRITE(IO_STDOUT,*)'global sum p-e+r (Sv): ', zzzhh * 1.e-6_wp
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

END
