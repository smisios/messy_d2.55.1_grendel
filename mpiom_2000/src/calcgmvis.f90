      SUBROUTINE CALCGMVIS

      USE mo_kind, ONLY: wp
      USE MO_PARAM1
      USE MO_PARALLEL
      USE mo_boundsexch, ONLY : bounds_exch
      USE MO_COMMO1
      USE mo_planetary_constants, ONLY: g, rhoref_water
      USE MO_UNITS, ONLY: io_stdout
      IMPLICIT NONE
      REAL(wp) :: shelp(ie,je),thelp(ie,je)
      REAL(wp) :: corio, visn, dudo, gor, &
           hho, richi, rinumo, visfac, vislx, visly, vismax
      !:: SET LOWER LIMIT FOR DEPTH INTEGRATED RICHARDSON NUMBER
      REAL(wp), PARAMETER :: depdown = 1000.0_wp
      !:: SET UPPER LIMIT FOR DEPTH INTEGRATED RICHARDSON NUMBER
      REAL(wp), PARAMETER :: depup = 20.0_wp
      REAL(wp), PARAMETER :: visal = 0.005_wp
      REAL(wp), PARAMETER :: vismin = 25.0_wp
      INTEGER :: kdown, ko, kvisd, kvisu, i, j, k
!
!-----------------------------------------------------------------------
!
!       CALCULATION OF RICHARDSON NUMBER DEPENDENT
!       COEFFICIENTS AFTER VISBECK ET AL. JPO 27, 1997
!
!
!         RINUM : ( (G/RHO)*D(RHO)/DZ ) / ( (D(VELOCITY)/DZ)**2 )
!                  RICHARDSON NUMBER (EVEN,ODD ==> RINUME,RINUMO)
!         GOR   : G/RHO
!
!
!---------------------------------------------------------------------
!
      REAL(wp) :: VISRICH_M(IE,JE),VISROSSBY_M(IE,JE)
!=======================================================================
!
      VISMAX = 2.0_wp * CAH00
      KVISD=1
      KVISU=1
      DO k=1,ke
       IF (DEPUP .GE. TIESTU(K)) KVISU= K
       IF (DEPDOWN .GE. TIESTU(K)) KVISD= K
      ENDDO
!
!
      WRITE(IO_STDOUT,*) ' CALCULATING VISBECK COEFFICIENT BETWEEN'
      WRITE(IO_STDOUT,*) ' K= ',KVISU,' AND K= ',KVISD
      GOR=G/RHOREF_WATER
!
!--------------------------------------------------------------------
      DO I=1,IE
       DO J=1,JE
        visrich_m(i, j) = 0.0_wp
        visrossby_m(i, j) = 0.0_wp
        thelp(i, j) = 0.0_wp
        shelp(i, j) = 0.0_wp
        BOLX(I,J)=MAX(BOLX(I,J),VISMIN)
        BOLY(I,J)=MAX(BOLY(I,J),VISMIN)
       ENDDO
      ENDDO
!
! CALCULATE VERTICALLY AVERAGED RICHARDSEN NR. FOR
! CALCULATION OF VISBECK COEFFICIENTS
! AFTER GENT ET AL. J. CLIM. 2002
!
      DO J=2,JE1
       DO I=2,IE-1
        KDOWN=MIN(KBOT(I,J),KVISD)
        IF (KDOWN .GT. KVISU) THEN
        DO K=KVISU,KDOWN
         KO=MAX(K-1,1)
!
         DUDO=ALMZER + WETO(I,J,K) * DI(K)**2                           &
     & * (   ( UKO(I-1,J,K) - UKO(I-1,J,KO) )**2                        &
     &     + ( UKO(I,J,K)   - UKO(I,J,KO)   )**2                        &
     &     + ( VKE(I,J-1,K) - VKE(I,J-1,KO) )**2                        &
     &     + ( VKE(I,J,K)   - VKE(I,J,KO)   )**2   ) * 0.5_wp
!
         HHO=WETO(I,J,K)*(AMSUO(I,J,K)+AMSUO(I-1,J,K)+AMSUE(I,J-1,K)    &
     & +AMSUE(I,J,K)) * 0.25_wp
!
         rinumo = hho * MAX(gor * stabio(i,j,k) / dudo, 0._wp)
! INTEGRATE VERTICALLY
           VISRICH_M(I,J)=VISRICH_M(I,J)+RINUMO*DZ(K)
           visn = SQRT(hho * MAX(gor * stabio(i,j,k), 0._wp))
!:: USE THELP AS TEMP FIELD TO SUM N
           THELP(I,J)=THELP(I,J)+VISN*DZ(K)
!:: USE SHELP AS TEMP FIELD TO SUM EFFECTIVE DEPTH RANGE
           SHELP(I,J)=SHELP(I,J)+DZ(K)
!
        ENDDO
       ENDIF
!
       ENDDO
      ENDDO
!
!OtB_GMVISETAL_TEST ifdef GMVISETAL_TEST
!OtB_GMVISETAL_TEST      ii1=ldt
!OtB_GMVISETAL_TEST      ii4=ie*je
!OtB_GMVISETAL_TEST      ii3=-100
!OtB_GMVISETAL_TEST      ii2=251
!OtB_GMVISETAL_TEST      WRITE(161) ii1,ii2,ii3,ii4
!OtB_GMVISETAL_TEST      WRITE(161) ((VISRICH_M(I,J)/(SHELP(I,J)+ALMZER),I=1,IE),J=1,JE)
!OtB_GMVISETAL_TEST      ii2=252
!OtB_GMVISETAL_TEST      WRITE(161) ii1,ii2,ii3,ii4
!OtB_GMVISETAL_TEST      WRITE(161) ((THELP(I,J)/(SHELP(I,J)+ALMZER),I=1,IE),J=1,JE)
!OtB_GMVISETAL_TEST      ii2=253
!OtB_GMVISETAL_TEST      WRITE(161) ii1,ii2,ii3,ii4
!OtB_GMVISETAL_TEST      WRITE(161) ((SHELP(I,J),I=1,IE),J=1,JE)
!OtB_GMVISETAL_TEST      igm=56
!OtB_GMVISETAL_TEST      jgm=24
!OtB_GMVISETAL_TEST      WRITE(*,*) 'gmvis at 56,24'
!OtB_GMVISETAL_TEST      WRITE(*,*) 'depto',depto(igm,jgm)
!OtB_GMVISETAL_TEST      WRITE(*,*) visrich_m(igm,jgm),thelp(igm,jgm),shelp(igm,jgm)
!OtB_GMVISETAL_TEST endif /*GMVISETAL_TEST*/
!
!:: CALCULATE VISBECK ET AL. DIFFUSION COEFFICIENTS FOR GM
        DO I=2,IE1
         DO J=2,JE1
          IF (lweto(i,j,1)) THEN
!:: CALCULATE ROSSBY RADIUS
          corio = ABS(0.25_wp * (ftwou(i, j) + ftwou(i-1, j) &
               + ftwov(i, j) + ftwov(i, j-1)))
!:: LIMIT F TO 2.5 degree off equator
          corio = MAX(2.5e-6_wp, corio)
          VISROSSBY_M(I,J)=THELP(I,J)/CORIO
          visrossby_m(i,j) = MIN(2.5e6_wp, visrossby_m(i,j))
!:: CALCULATE RICHARDSON_NUMBER
          RICHI=VISRICH_M(I,J)/(SHELP(I,J)+ALMZER)
          VISLX=MAX(VISROSSBY_M(I,J),DLXP(I,J))
          VISLY=MAX(VISROSSBY_M(I,J),DLYP(I,J))
          VISFAC=CORIO*VISAL/SQRT(RICHI+ALMZER)
          IF (shelp(i, j) .LT. almzer) visfac = 0.0_wp
! STORE ROSSBY NR. FOR DIAGNOSTICS
          VISROSSBY_M(I,J)=VISLX
          thelp(i, j) = MIN(vismax, MAX(vismin, vislx**2 * visfac))
          shelp(i, j) = MIN(vismax, MAX(vismin, visly**2 * visfac))
          ELSE
           THELP(I,J)=VISMIN
           SHELP(I,J)=VISMIN
          ENDIF
        ENDDO
       ENDDO
!
       CALL bounds_exch(1,'p',THELP,'calcgmvis 10')
       CALL bounds_exch(1,'p',SHELP,'calcgmvis 10')
       CALL bounds_exch(1,'p',VISRICH_M,'calcgmvis 10')
       CALL bounds_exch(1,'p',VISROSSBY_M,'calcgmvis 10')
!
! TAKE 0.5 times old day bolx/boly
       DO I=1,IE1
        DO J=1,JE1
         bolx(i, j) = 0.5_wp * bolx(i, j) &
              + 0.25_wp * (thelp(i, j) + thelp(i+1, j))
         boly(i, j) = 0.5_wp * boly(i, j) &
              + 0.25_wp * (shelp(i, j) + shelp(i, j+1))
        ENDDO
       ENDDO

       CALL bounds_exch(1,'u',BOLX,'calcgmvis 10')
       CALL bounds_exch(1,'v',BOLY,'calcgmvis 10')

!OtB_GMVISETAL_TEST ifdef GMVISETAL_TEST
!OtB_GMVISETAL_TEST      ii1=ldt
!OtB_GMVISETAL_TEST      ii4=ie*je
!OtB_GMVISETAL_TEST      ii3=-100
!OtB_GMVISETAL_TEST      ii2=254
!OtB_GMVISETAL_TEST      WRITE(161) ii1,ii2,ii3,ii4
!OtB_GMVISETAL_TEST      WRITE(161) ((THELP(I,J),I=1,IE),J=1,JE)
!OtB_GMVISETAL_TEST      ii2=255
!OtB_GMVISETAL_TEST      WRITE(161) ii1,ii2,ii3,ii4
!OtB_GMVISETAL_TEST      WRITE(161) ((SHELP(I,J),I=1,IE),J=1,JE)
!OtB_GMVISETAL_TEST      ii2=256
!OtB_GMVISETAL_TEST      WRITE(161) ii1,ii2,ii3,ii4
!OtB_GMVISETAL_TEST      WRITE(161) ((VISROSSBY_M(I,J),I=1,IE),J=1,JE)
!OtB_GMVISETAL_TEST endif /*GMVISETAL_TEST*/
!
!=====================================================================

      RETURN
      END
