SUBROUTINE GROWTH(HO,AO,HSNO,HN,AN,HSNN,                          &
#ifdef __coupled
     RPRECW,RPRECI,AOFLNHW,AOFLCHI,AOFLRHI,                       &
#else
     RPREC,TAIR,TD,ACL,PA,UG,                                     &
#endif
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
#ifdef __coupled
  !     *RPRECW*   (P-E)WATER+RUNOFF+(P-E)ICE:IF RAIN     [M/S,EQ. WATER COLUMN]
  !     *RPRECI*   (P-E)ICE:IF SNOW FALL OR SUBLIMATION   [M/S, " ]
  !     *AOFLNHW*  NET HEAT FLUX OVER WATER               [W/M**2]
  !     *AOFLCHI*  CONDUCTIVE HEAT FLUX > 0 IF UPWARD     [W/M**2]
  !     *AOFLRHI*  RESIDUAL HEAT FLUX (FOR SNOW/ICE MELT) [W/M**2]
  !                > 0 IF DOWNWARD
#else
  !     *RPREC*    ATMOSPHERIC P-E                        [M/S]
  !     *TAIR*     AIR TEMPERATURES                       [DEG C]
  !     *TD*       DEW POINT TEMPERATURES                 [K]
  !     *ACL*      FRACTIONAL CLOUD COVER                 [FRAC.]
  !     *PA*       ATMOSPHERIC SURFACE PRESURE            [PA]
  !     *UG*       WIND SPEED                             [M/S]
#endif
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

#ifndef __coupled
#ifdef CORE
  ! CORE forcing currently needs some additional calculation
  USE mo_commo1, ONLY: sicuo, sicve, uko, vke, weto, sao, dlxp, dlyp, txo, tye
  USE mo_commoau1, ONLY: tmelt
  USE mo_boundsexch, ONLY : bounds_exch
  USE mo_ncar_ocean_fluxes, ONLY: budget_ocean_core, budget_ice_core, &
       corrunoff, normpem, albi, albm, albsn, albsnm
#else/*def CORE */
  ! OMIP forcing is default
  USE mo_omip, ONLY: budget, obudget, albi, albm, albsn, albsnm
#endif/*else def CORE */
#endif/*ndef __coupled */

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
#ifdef __coupled
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
#else/*def __coupled */


  REAL(wp) :: precn(ie,je),prech(ie,je)

#ifdef CORE
  REAL(wp) :: rprecw(ie,je),rpreci(ie,je),ti(ie,je)
  REAL(wp) :: qpre(ie,je),uo(ie,je),vo(ie,je),to(ie,je)
  REAL(wp) :: txic(ie,je),txoc(ie,je),tyic(ie,je),tyoc(ie,je)
  REAL(wp) :: rai, rao
  REAL(wp) :: stp(ie,je),stpp(ie,je),fp(ie,je),fpp(ie,je),fdiff
  REAL(wp) :: rcorfw, rsumfw, rarea

  INTEGER, PARAMETER :: imax = 10
  INTEGER :: iter
#endif/*def CORE */

  REAL(wp) :: ho(ie,je),ao(ie,je),hsno(ie,je)
  REAL(wp) :: tair(ie,je),rprec(ie,je)
  REAL(wp) :: td(ie,je),acl(ie,je),pa(ie,je)
  REAL(wp) :: sln(ie,je),z(ie,je)
  REAL(wp) :: om(ie,je)
  REAL(wp) :: fw(ie,je),sh(ie,je)
  REAL(wp) :: qnsw(ie,je),qnlw(ie,je),qnse(ie,je),qnla(ie,je)
  REAL(wp) :: hn(ie,je),an(ie,je),hsnn(ie,je)
  REAL(wp) :: qs(ie,je),qt(ie,je)



  !  LOCAL VARIABLES

  REAL(dp) :: TICM(IE,JE), TICA(IE,JE)
  REAL(dp) :: UG(IE,JE),HEATABS(IE,JE)
  REAL(dp) :: FO(IE,JE)
  REAL(dp) :: QH(IE,JE),HS(IE,JE),HT(IE,JE)
  REAL(dp) :: HDRAFT(IE,JE)
  REAL(dp) :: RA(IE,JE),RHS(IE,JE),RHB(IE,JE),RH(IE,JE)
  REAL(dp) :: RHSA(IE,JE),RHBA(IE,JE)
  REAL(dp) :: QHST(IE,JE),SN(IE,JE),QFM(IE,JE)
  REAL(dp) :: H(IE,JE,2),A(IE,JE,2),HSN(IE,JE,2)
  REAL(dp) :: QSW(IE,JE),QLW(IE,JE),QSE(IE,JE),QLA(IE,JE)
  REAL(dp) :: TMP(IE,JE),TMP2(IE,JE),TMP3(IE,JE),TMP4(IE,JE)
#endif/*else def __coupled */

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

#ifndef __coupled

        qsw(i, j) = 0.0_wp
        qlw(i, j) = 0.0_wp
        qse(i, j) = 0.0_wp
        qla(i, j) = 0.0_wp
        qnsw(i, j) = 0.0_wp
        qnlw(i, j) = 0.0_wp
        qnse(i, j) = 0.0_wp
        qnla(i, j) = 0.0_wp
#endif

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

#ifndef __coupled

#ifdef CORE
  DO i=1,ie
     DO j=1,je
        uo(i,j)=uko(i,j,1)
        vo(i,j)=vke(i,j,1)
        to(i, j) = tho(i, j, 1) + tmelt
     END DO
  END DO

    CALL bounds_exch(1,'u',uo,'growth 1')
    CALL bounds_exch(1,'v',vo,'growth 2')
    CALL bounds_exch(1,'p',to,'growth 3')

  CALL budget_ocean_core(uo,vo,to,qla,qse,qsw,qlw,qpre,fo,rprecw,txoc,tyoc)

#else

! in: tair,qt,sln,om,TD,ACL,PA,UG
! out: FO,QSW,QLW,QSE,QLA

   DO i=1,ie
     DO j=1,je
        tair(i,j)=tafo(i,j)
        qt(i,j)=tho(i,j,1)
        TD(i,j)=ftdew(i,j)
        ACL(i,j)=fclou(i,j)
!        PA(i,j)=fslp(i,j)
        UG(i,j)=fu10(i,j)
        SLN(i,j)=fswr(i,j)
        rprec(i,j)=fprec(i,j)
!   +giriv(i,j)
     END DO
   END DO

  CALL OBUDGET(QT,TAIR,TD,ACL,PA,UG,SLN,FO,OM,QSW,QLW,QSE,QLA)

#endif /* CORE */


  vapli = 1._wp / vapl

#endif /*__coupled*/

  DO J=1,M
     DO I=1,L

#ifdef __coupled


        FO(I,J) = -AOFLNHW(I,J)/CLB                     &
             +SWSUM(I,J)*SLN(I,J)/CLB

        !SV 31.08.99 SET HEAT FLUX DUMMY FIELDS FOR UWE'S DIAGNOSTICS ONLY
        !SV   (QNSE IS SET TO QNET-QSOLAR AND QNSW IS SET TO QSOLAR; OTHERS ARE SET TO ZERO)
        QNSE(I,J) = (AOFLNHW(I,J)-SLN(I,J))*(1.-A(I,J,LRHS))
        QNSW(I,J) = SLN(I,J)*(1.-A(I,J,LRHS))
        QNLW(I,J) = 0.0
        QNLA(I,J) = 0.0

        HEATABS(I,J)=SWSUM(I,J)*SLN(I,J)*(1.-A(I,J,LRHS))

#else /*__coupled*/

        FO(I,J)=FO(I,J)+SWSUM(I,J)*QSW(I,J)/CLB
        heatabs(i, j) = swsum(i, j) * qsw(i, j) * (1._wp - a(i, j, lrhs))


        qnsw(i, j) = qsw(i, j) * (1._wp - a(i, j, lrhs))
        qnlw(i, j) = qlw(i, j) * (1._wp - a(i, j, lrhs))
        qnse(i, j) = qse(i, j) * (1._wp - a(i, j, lrhs))
        qnla(i, j) = qla(i, j) * (1._wp - a(i, j, lrhs))

!#ifdef TESTOUT_HFL
 if (LHFLDIAG) then
        DQSWO(I,J)=QSW(I,J)
        DQLWO(I,J)=QLW(I,J)
        DQSEO(I,J)=QSE(I,J)
        DQLAO(I,J)=QLA(I,J)
        DQTHO(I,J)=QSW(I,J)+QLW(I,J)+QSE(I,J)+QLA(I,J)
 endif
!#endif

        !UWE     STORE EVAPORATION IN TMPUW3
        tmpuw3(i, j) = dt * om(i, j) * qla(i, j) * (1._wp - a(i, j, lrhs)) &
             * vapli

#endif /*__coupled*/

        !-----------------------------------------------------------------------
        !  THICK ICE GROWTH RATES.
        !-----------------------------------------------------------------------

#ifdef __coupled
        RHS(I,J) = -AOFLRHI(I,J) / CLB
        RHB(I,J) = -AOFLCHI(I,J) / CLB
#else /*__coupled*/
        !-----------------------------------------------------------------------
        !  CALCULATE EFFECTIVE ICE THICKNESS FOR CONDUCTION TERM IN BUDGET:
        !               FIRST MAKE SURE WE HAVE NON-ZERO COMPACTNESS.
        !               INCLUDE SNOW THICKNESS BECAUSE OF CONDUCTION EFFECT.
        !               MAKE SURE TO HAVE NON-ZERO EFFECTIVE ICE THICKNESS.
        !               INITIALISE MEAN ICE GROWTH RH AND TEMPERATURE (TICM).
        !               STORE EFFECTIVE ICE THICKNESS IN ARRAY SN.
        !-----------------------------------------------------------------------
        !
        TMP (I,J) = ( H(I,J,LRHS) + HSN(I,J,LRHS)*CON/CONSN )          &
             /MAX(A(I,J,LRHS),ARMIN)
        rhs(i, j) = 0.0_wp
        rhb(i, j) = 0.0_wp
        rhsa(i, j) = 0.0_wp
        rhba(i, j) = 0.0_wp
        ticm(i, j) = 0.0_wp
        SN(I,J)=MAX(TMP(I,J),HMIN)
        !
        !-----------------------------------------------------------------------
        !     SET ALBEDO (TMP2) ACCORDING TO PRESENCE OF SNOW  (TMP4)
        !                                 TO MELTING CONDITIONS(TMP3).
        !                                 USEMEAN PREVIOUS ICE TEMP..
        !-----------------------------------------------------------------------
        !
        !UWE    INCLUDE FINITE VALUE FOR SNOW THICKNESS IN ALBEDO CALCULATION!
        !
        tmp4(i, j) = (0.5_wp - SIGN(0.5_wp, 1.E-2_wp - hsn(i, j, lrhs)))
        tmp3(i, j) = (0.5_wp + SIGN(0.5_wp, tice(i, j)))
        tmp2(i, j) = tmp4(i, j) &
             * (tmp3(i, j) * albsnm + (1._wp - tmp3(i, j)) * albsn) &
             + (1._wp - tmp4(i, j)) * (tmp3(i, j) * albm &
             &                         + (1._wp - tmp3(i, j)) * albi)
        !
        !-----------------------------------------------------------------------
        !  CALCULATE GROWTH RATES FOR THICK ICE.
        !  DO ONCE FOR EACH ICE THICKNESS CATEGORIE ( ICELEV).
        !  USE THE SAME MEAN ICE SURFACE TEMPERATURE OF THE PREVIOUS TIME
        !  STEP FOR ALL THICKNESS CATEGORIES (LOOP 6).
        !-----------------------------------------------------------------------
        !
        !UWE   INITIALIZE TEMPORARY FILEDS FOR SNOW EVAPORATION
        !
        tmpuwe(i, j) = 0._wp
        tmpuw2(i, j) = 0._wp
        tmpuw4(i, j) = 0._wp
#endif /*__coupled*/
     END DO
  END DO


#ifdef __coupled

  DO J=1,M
     DO I=1,L
        QNSE(I,J) = QNSE(I,J)+(AOFLRHI(I,J)+AOFLCHI(I,J))*A(I,J,LRHS)
     ENDDO
  ENDDO

#else /*__coupled*/

  ICELEV=1
  anzlevi = 1._wp / REAL(icelev, wp)

  DO K=1,ICELEV
     DO J=1,M
        DO I=1,L
           tmp(i, j) = REAL(2 * K - 1, wp) * sn(i, j) * anzlevi
           TICA(I,J) = TICE(I,J)
        END DO
     END DO

#ifdef CORE

     !-----------------------------------------------------------------------
     !  make first guess for surface temperature and heat balance
     !-----------------------------------------------------------------------

     DO i=1,iE
        DO j=1,jE
          tica(i, j) = tica(i, j) + tmelt
        END DO
     END DO

     CALL budget_ice_core(sicuo,sicve,tica,tmp,           &
          qla,qse,qsw,qlw,rhsa,rhba,rpreci,txic,tyic)

     DO i=1,ie
        DO j=1,je
           stp(i,j)=tica(i,j)
           fp(i,j)=rhsa(i,j)
           tica(i, j) = tica(i, j) + 1._wp
        END DO
     END DO

     !-----------------------------------------------------------------------
     !  calculate the surface temperature (start of iteration procedure)
     !-----------------------------------------------------------------------

     DO iter=1,imax
        DO i=1,ie
           DO j=1,je
              stpp(i,j)=stp(i,j)
              fpp (i,j)=fp(i,j)
              stp (i,j)=tica(i,j)
           END DO
        END DO

        CALL budget_ice_core(sicuo,sicve,tica,tmp,           &
             qla,qse,qsw,qlw,rhsa,rhba,rpreci,txic,tyic)


       DO i=1,ie

           DO j=1,je
              fp(i,j)  = rhsa(i,j)
              fdiff  = fp(i,j)-fpp(i,j)
              tica(i,j)  = MAX(stp(i,j)-(stp(i,j)-stpp(i,j))*fp(i,j)       &
                   /(ABS(fdiff) + almzer) * SIGN(1._wp, fdiff),100._wp)
           END DO
        END DO
     END DO

     DO i=1,ie
        DO j=1,je
           rai=a(i,j,lrhs)
           rao = 1._wp - rai
           txo(i,j)=(rao*txoc(i,j)+rai*txic(i,j))/rhoref_water
           tye(i,j)=(rao*tyoc(i,j)+rai*tyic(i,j))/rhoref_water
           rprec(i,j)=(rao*rprecw(i,j)+rai*rpreci(i,j))
           rprec(i,j)=rprec(i,j)+(corrunoff(i,j)/(dlxp(i,j)*dlyp(i,j)))
       END DO
     END DO



     DO i=1,iE
        DO j=1,jE
          tica(i, j) = tica(i, j) - tmelt
        END DO
     END DO

#else /*CORE*/

!in: LRHS,A,TAIR,TD,ACL,PA,UG,SLN,OM,TMP,TMP2,TMP4

!out: RHSA,RHBA,QSW,QLW,QSE,QLA
!inout: TICA

     CALL BUDGET(RHSA,RHBA,TICE,LRHS,A,                           &
          TAIR,TD,ACL,PA,UG,SLN,OM,TMP,TMP2,TMP4,QSW,QLW,QSE,QLA)

#endif /*CORE*/

     !UWE   INCLUDE SUBLIMATION OF SNOW AND ICE FROM EVAPORATION

     subrsni = 1._wp / (subl * rhoref_snow)
     subrici = 1._wp / (subl * rhoref_ice)
     RWRI=RHOREF_SNOW/RHOREF_ICE
     subli = 1._wp / subl

     DO I=1,IE
        DO J=1,JE
           TICM(I,J) = TICM(I,J)+TICA(I,J)*ANZLEVI
           RHS (I,J) = RHS (I,J)+RHSA(I,J)*ANZLEVI
           RHB (I,J) = RHB (I,J)+RHBA(I,J)*ANZLEVI
           QNSW(I,J) = QNSW(I,J)+QSW(I,J)*A(I,J,LRHS)*ANZLEVI
           QNSE(I,J) = QNSE(I,J)+QSE(I,J)*A(I,J,LRHS)*ANZLEVI
           QNLW(I,J) = QNLW(I,J)+QLW(I,J)*A(I,J,LRHS)*ANZLEVI
           QNLA(I,J) = QNLA(I,J)+QLA(I,J)*A(I,J,LRHS)*ANZLEVI

!#ifdef TESTOUT_HFL
 if (LHFLDIAG) then
          DQSWI(I,J)=QSW(I,J)
           DQLWI(I,J)=QLW(I,J)
           DQSEI(I,J)=QSE(I,J)
           DQLAI(I,J)=QLA(I,J)
           DQTHI(I,J)=QSW(I,J)+QLW(I,J)+QSE(I,J)+QLA(I,J)
           IF (a(i, j, lrhs) .LT. 1.e-3_wp) THEN
              dticeo(i, j) = 99._wp
           ELSE
              DTICEO(i,j)=TICM(I,J)
           ENDIF
endif
!#endif

!     TMPUWE CONTAINS AMOUNT OF SNOW AND ICE EVAPORATED DURING THIS TIME STEP
!     AVERAGED OVER GRID BOX   [KG/M**2]

           TMPUWE(I,J)=TMPUWE(I,J)+A(I,J,LRHS)*QLA(I,J)                    &
                *ANZLEVI*OM(I,J)*SUBLI*DT
        END DO
     END DO

  END DO !ICELEV



!UWE        COMPUTE SUM OF SNOW AND ICE THICKNESS EVAPORATED AWAY.
!        IF NOT SUFFICIENT AVAILABLE, USE THE ENERGY TO EVAPORATE WATER
!        STORED IN TMPUW3

#endif /* __coupled*/

  DO J=1,M
     DO I=1,L

#ifndef __coupled

        TMPUW2(I,J)=-MIN(HSN(I,J,LNEW)*RHOREF_SNOW                         &
             +RHOREF_ICE*H(I,J,LNEW)                                &
             ,-TMPUWE(I,J))

        TMPUW3(I,J)=TMPUW3(I,J)+(TMPUWE(I,J)-TMPUW2(I,J))*SUBL/VAPL
        tmpuw4(i, j) = MIN(0._wp, tmpuw2(i, j) + hsn(i, j, lnew) * rhoref_snow)
        hsn(i, j, lnew) = MAX(0._wp, &
             hsn(i, j, lnew) + tmpuwe(i, j) / rhoref_snow)
        h(i, j, lnew) = MAX(0._wp, h(i, j, lnew) + tmpuw4(i, j) / rhoref_ice)

        !       TMPUW2 CONTAINS SUBLIMATION OF ICE IN KG/M**2
        !       TMPUW3 CONTAINS EVAPORATION OVER ICE FREE WATER IN KG/M**2

        !-----------------------------------------------------------------------
        !  STORE MEAN ICE TEMPERATURE ON TICE.
        !  DETERMINE THERMODYNAMIC TOTAL ICE THICKNESS CHANGE (SH)
        !  FOR CONTINUITY EQUATION.
        !  MULTIPLY THICK ICE GROWTH RATE WITH COMPACTNESS A(,,LRHS) (SNOW DEPTH
        !  IS MEAN OVER GRID CELL AREA).
        !  CONVERT P-E OVER ICE INTO SNOWFALL IF TAIR < 0.
        !  ADD SNOWFALL TO SNOW LAYER OF NEW TIME STEP.
        !  SUBTRACT SNOWFALL FROM PRECIPITATION INTO MIXED LAYER.
        !-----------------------------------------------------------------------

#endif /*__coupled*/

        RA (I,J)  = FO (I,J)*DT
        RHS(I,J)  = RHS(I,J)*DT
        RHB(I,J)  = RHB(I,J)*DT
        RH(I,J)  = RHS(I,J)+RHB(I,J)


     ENDDO
  ENDDO


  DO J=1,M
     DO I=1,L

#ifdef __coupled

        sh(i, j) = rh(i, j) * a(i, j, lrhs) + ra(i, j) * (1._wp - a(i, j, lrhs))
        !
        PRECN(I,J)    = RPRECW(I,J)*DT
        HSN(I,J,LNEW) = HSN(I,J,LNEW)+RPRECI(I,J)                      &
             *DT*RHOREF_WATER/RHOREF_SNOW

        Z(I,J)        = Z(I,J) + (RPRECW(I,J)+RPRECI(I,J))*DT

        !UWE INCLUDE COMPACTNESS   22.12.99
        TMP(I,J)      = RHS(I,J)*A(I,J,LRHS)

#else /*__coupled*/
        ! thermodynamic sea ice thickness change
        sh(i, j) = rh(i, j) * a(i, j, lrhs) + ra(i, j) * (1._wp - a(i, j, lrhs))

        tmp3(i, j)     = (0.5_wp - SIGN(0.5_wp, tair(i, j))) &
             * a(i, j, lrhs) * REAL(isnflg, wp)         ! snowfall flag 1 or 0
        precn(i, j) = ((1._wp - tmp3(i, j)) * rprec(i, j)) * dt * om(i, j)                     ! liquid part
        HSN(I,J,LNEW) = HSN(I,J,LNEW)+TMP3(I,J)*RPREC(I,J)*DT*OM(I,J)*RHOREF_WATER/RHOREF_SNOW ! solid part

        Z(I,J)        = Z(I,J) + RPREC(I,J)*DT*OM(I,J)                              ! precip changes sea level


        !UWE     ADD EVAPORATION TO SEA LEVEL

        Z(I,J)        = Z(I,J)+OM(I,J)*(TMPUW2(I,J)+TMPUW3(I,J))/RHOREF_WATER             ! evap changes to sea level
        PRECN(I,J)    = PRECN(I,J)+TMPUW3(I,J)/RHOREF_WATER                               ! liquid part
        TMP(I,J)      = RHS(I,J)*A(I,J,LRHS)                                        ! heatflux ???

#endif /*__coupled*/

#ifdef __coupled
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
#endif  /*__coupled*/

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

#ifdef __coupled
        QHST(I,J)    =  H(I,J,LNEW)                                    &
             - ((qt(i,j) - tfreeze) * qh(i,j)) * cc/clb * om(i, j)     &
             + RH(I,J)
#else /*__coupled*/
        QHST(I,J)    =  H(I,J,LNEW)                                    &
             - ((qt(i,j) - tfreeze) * qh(i,j)                          &
             &  + (tair(i,j) - tfreeze) * precn(i,j)) * cc/clb*om(i,j) &
             + RH(I,J)
#endif /*__coupled*/

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


#ifndef __coupled

        QHST(I,J) = QHST(I,J)+(HSN(I,J,LNEW)-SN(I,J))*RHOREF_SNOW/RHOREF_ICE      &
             - (HSN(I,J,LNEW)-SN(I,J))*RHOREF_SNOW/RHOREF_WATER             &
                                !UWE USE 0. INSTEAD OF TAIR FOR SNOW MELT
             * (0._wp - tfreeze) * cc/clb

#else /*__coupled*/

        QHST(I,J) = QHST(I,J)+(HSN(I,J,LNEW)-SN(I,J))*RHOREF_SNOW/RHOREF_ICE

#endif /*__coupled*/

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
#ifdef __coupled
           PRECH(I,J)=RPRECW(I,J)+RPRECI(I,J)
#else
           PRECH(I,J)=RPREC(I,J)+(TMPUW2(I,J)+TMPUW3(I,J))/(RHOREF_WATER*DT)
#endif
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
