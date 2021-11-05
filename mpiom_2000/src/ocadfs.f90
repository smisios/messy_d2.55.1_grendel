      SUBROUTINE OCADFS
!**********************************************************************
!
!
!      OOO    CCCCC  AAAAAA   DDDDD   FFFFFF   SSS
!     O   O  C       A    A   D    D  F       S
!     O   O  C       AAAAAA   D    D  FFF      SSS
!     O   O  C       A    A   D    D  F           S
!      OOO    CCCCC  A    A   DDDDD   F        SSS
!
!
!=======================================================================
!
!     SBR OCADFS
!
!
!     PURPOSE :
!
!     A) 3-D ADVECTION OF TEMPERATURE AND SALINITY
!        PREDICTOR-CORRECTOR SCHEME
!        IF IOCAD=6 (replaces CPP key QUICK) THEN
!          USING QUICK-SCHEME AS PROPOSED BY FARROW AND STEVENS 1995
!          JPO 25, 1731...
!          3RD ORDER ACCURATE IN SPACE, 2ND IN TIME
!        ELSEIF IOCAD=7 (replaces CPP key QUICK2)
!          USE TEMPERATURE OF THE INFLOWING WATER IN I-1 INSTEAD OF
!          TEMPERATURE IN I-2
!          SEEMS TO BE MORE EFFICIENT IN REDUCING OVERSHOOTINGS.
!          BETTER CONDITIONED AT BOUNDARIES
!        ELSE
!          ORDINARY CENTERED DIFFERENCES
!
!
!      BY UWE MIKOLAJEWICZ    1/99
!      OPTION QUICK2 ADDED    3/99
!      DIFFUSION REMOVED      5/99
!      Changed for MPI-Parallelization R. Johanni 11/03
!
      USE MO_PARAM1
      USE mo_boundsexch, ONLY : bounds_exch
      USE MO_COMMO1
      USE MO_COMMOAU1
      USE mo_planetary_constants, ONLY: rhoicwa, rhosnwa
      USE MO_PARALLEL

      IMPLICIT NONE

      REAL(wp) :: CURVUPT(IE,JE),CURVLOT(IE,JE),CURVEAT(IE,JE)            &
     &         ,CURVWET(IE,JE),CURVNOT(IE,JE),CURVSOT(IE,JE)
      REAL(wp) :: CURVUPS(IE,JE),CURVLOS(IE,JE),CURVEAS(IE,JE)            &
     &         ,CURVWES(IE,JE),CURVNOS(IE,JE),CURVSOS(IE,JE)
      REAL(wp) :: ZSURF(KE)
      REAL(wp) :: DLXP2(IE,JE),DLYP2(IE,JE)

      REAL(wp), ALLOCATABLE :: dsx(:,:,:), dsy(:,:,:), dtx(:,:,:), dty(:,:,:)
      REAL(wp), ALLOCATABLE :: tinf(:,:,:), sinf(:,:,:)

      REAL(wp) :: gew, sa, sm, ta, tm, turdif, uos, uwe, vno, vsu, wob, wun
      INTEGER :: i, j, k, klo, klolo, kup, kupup

      IF ( iocad == 6 ) THEN
        ALLOCATE ( dsx(ie,je,ke), dsy(ie,je,ke), dtx(ie,je,ke), dty(ie,je,ke) )
      ENDIF

      IF ( iocad == 7 ) THEN
        ALLOCATE ( tinf(ie,je,ke), sinf(ie,je,ke) )
      ENDIF

      DO K=1,KE
       zsurf(k) = 0._wp
      ENDDO
      zsurf(1) = 1._wp
!
!C      WRITE(IO_STDOUT,*)'FS ADVEKTION!!!!'
!---------------------------------------------------------------------
!
!     A.1)
!
!     1. HALF STEP IN TIME
!
!
!---------------------------------------------------------------------
!
!     A.1.1)
!
!----------------------------------------------------------------------
!
!     A.1.2)
!
!     ODD FIELDS
!
!
      DO K=1,KE
!
!
      KUP = MAX(1,K-1)
      KLO = MIN(KE,K+1)
!
!
      DO J=2,JE1
      DO I=2,IE1
!
      TM = THO(I,J,K)
      SM = SAO(I,J,K)
!
      WUN = WO(I,J,K+1)*DLXP(I,J)*DLYP(I,J)                             &
     &     /(DZW(KLO)+DZW(K))
      WOB = WO(I,J,K)*DLXP(I,J)*DLYP(I,J)                               &
     &     /(DZW(KUP)+DZW(K))
      UWE = UKO(I-1,J,K)*DDUO(I-1,J,K)*DLYU(I-1,J)                      &
     &     /(DLXP(I,J)+DLXP(I-1,J))
      UOS = UKO(I,J,K)*DDUO(I,J,K)*DLYU(I,J)                            &
     &     /(DLXP(I,J)+DLXP(I+1,J))
      VSU = VKE(I,J,K)*DDUE(I,J,K)*DLXV(I,J)                            &
     &     /(DLYP(I,J)+DLYP(I,J+1))
      VNO = VKE(I,J-1,K)*DLXV(I,J-1)*DDUE(I,J-1,K)                      &
     &     /(DLYP(I,J)+DLYP(I,J-1))
!
      t1o(i, j, k) = tho(i, j, k) + 0.5_wp * dt * weto(i, j, k) *(      &
     &         UWE * ( DLXP(I,J)*THO(I-1,J,K) +DLXP(I-1,J)*TM )         &
     &     -   UOS * ( DLXP(I+1,J)*TM + DLXP(I,J)*THO(I+1,J,K) )        &
     &     +   VSU * ( DLYP(I,J)*THO(I,J+1,K) + DLYP(I,J+1)*TM )        &
     &     -   VNO * ( DLYP(I,J-1)*TM + DLYP(I,J)*THO(I,J-1,K) )        &
     &     +   WUN * ( DZW(K)*THO(I,J,KLO) + DZW(KLO)*TM )              &
     &     -   WOB * ( DZW(KUP)*TM + DZW(K)*THO(I,J,KUP) )   )          &
!CUWE  INCLUDE ZETA IN UPPERMOST LEVEL THICKNESS 12/99
     &  /((DDPO(I,J,K)+ZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)            &
     &       -RHOSNWA*SICSNO(I,J))+ALMZER)                              &
     &      *DLXP(I,J)*DLYP(I,J))
!
!
      s1o(i, j, k) = sao(i, j, k) + 0.5_wp * dt * weto(i, j, k)*(       &
     &         UWE * ( DLXP(I,J)*SAO(I-1,J,K) +DLXP(I-1,J)*SM )         &
     &     -   UOS * ( DLXP(I+1,J)*SM + DLXP(I,J)*SAO(I+1,J,K) )        &
     &     +   VSU * ( DLYP(I,J)*SAO(I,J+1,K) + DLYP(I,J+1)*SM )        &
     &     -   VNO * ( DLYP(I,J-1)*SM + DLYP(I,J)*SAO(I,J-1,K) )        &
     &     +   WUN * ( DZW(K)*SAO(I,J,KLO) + DZW(KLO)*SM )              &
     &     -   WOB * ( DZW(KUP)*SM + DZW(K)*SAO(I,J,KUP) )   )          &
     &  /((DDPO(I,J,K)+ZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)            &
     &                          -RHOSNWA*SICSNO(I,J))+ALMZER)           &
     &      *DLXP(I,J)*DLYP(I,J))
!
      END DO
      END DO
      END DO
!
      CALL bounds_exch(1,'p',T1O)
      CALL bounds_exch(1,'p',S1O)
!
!----------------------------------------------------------------------
!
!     A.2)
!
!     2. HALF STEP IN TIME
!
!
!----------------------------------------------------------------------
!
      turdif = 1.e-3_wp
!
      IF ( iocad == 7 ) THEN
!
!UWE
!     COMPUTE T AND S OF INFLOWING WATER
!
      DO K=1,KE
       KLO=MIN(K+1,KE)
       KUP=MAX(K-1,1)
       DO J=2,JE-1
        DO I=2,IE-1
         WUN=HALF*(WO(I,J,K+1)+ABS(WO(I,J,K+1)) )*DLXP(I,J)*DLYP(I,J)
         WOB=HALF*(ABS(WO(I,J,K))-WO(I,J,K))*DLXP(I,J)*DLYP(I,J)
         UWE=HALF*(UKO(I-1,J,K)+ABS(UKO(I-1,J,K)))                      &
     &           *DLYU(I-1,J)*DDUO(I-1,J,K)
         UOS=HALF*(ABS(UKO(I,J,K))-UKO(I,J,K))*DLYU(I,J)*DDUO(I,J,K)
         VSU=HALF*(VKE(I,J,K)+ABS(VKE(I,J,K)))*DLYV(I,J)*DDUE(I,J,K)
         VNO=HALF*(ABS(VKE(I,J-1,K))-VKE(I,J-1,K))                      &
     &           *DLYV(I,J-1)*DDUE(I,J-1,K)
         GEW=ALMZER+WUN+WOB+UWE+UOS+VSU+VNO
         TINF(I,J,K)=(ALMZER*T1O(I,J,K)                                 &
     &        +WUN*T1O(I,J,KLO)+WOB*T1O(I,J,KUP)                        &
     &        +UWE*T1O(I-1,J,K)+UOS*T1O(I+1,J,K)                        &
     &        +VSU*T1O(I,J+1,K)+VNO*T1O(I,J-1,K))/GEW
         SINF(I,J,K)=(ALMZER*S1O(I,J,K)                                 &
     &        +WUN*S1O(I,J,KLO)+WOB*S1O(I,J,KUP)                        &
     &        +UWE*S1O(I-1,J,K)+UOS*S1O(I+1,J,K)                        &
     &        +VSU*S1O(I,J+1,K)+VNO*S1O(I,J-1,K))/GEW
        ENDDO
       ENDDO
       DO I=1,IE
        tinf(i, 1, k) = 0._wp
        sinf(i, 1, k) = 0._wp
        tinf(i, je, k) = 0._wp
        sinf(i, je, k) = 0._wp
       ENDDO
      ENDDO
      CALL bounds_exch(1,'p',TINF)
      CALL bounds_exch(1,'p',SINF)
!
      ENDIF ! (iocad==7)
!
!
!     Calculate some auxiliary arrays

      dlxp2(:,:) = 0._wp
      dlyp2(:,:) = 0._wp
      DO J=1,JE-1
      DO I=2,IE-1
        DLXP2(I,J) = DLXP(I,J)+DLXP(I+1,J)
        DLYP2(I,J) = DLYP(I,J)+DLYP(I,J+1)
      ENDDO
      ENDDO

      CALL bounds_exch(1,'u+',DLXP2)
      CALL bounds_exch(1,'v+',DLYP2)

      IF ( iocad == 6 ) THEN
      DO K=1,KE
        DO J=1,JE-1
        DO I=2,IE-1
          DSX(I,J,K) = (S1O(I+1,J,K)-S1O(I,J,K))/(DLXP(I,J)+DLXP(I+1,J))
          DTX(I,J,K) = (T1O(I+1,J,K)-T1O(I,J,K))/(DLXP(I,J)+DLXP(I+1,J))
          DSY(I,J,K) = (S1O(I,J,K)-S1O(I,J+1,K))/(DLYP(I,J)+DLYP(I,J+1))
          DTY(I,J,K) = (T1O(I,J,K)-T1O(I,J+1,K))/(DLYP(I,J)+DLYP(I,J+1))
        ENDDO
        ENDDO
        DO I=2,IE-1
          dsx(i, je, k) = 0._wp
          dtx(i, je, k) = 0._wp
          dsy(i, je, k) = 0._wp
          dty(i, je, k) = 0._wp
        ENDDO
      ENDDO

      CALL bounds_exch(1,'u+',DSX)
      CALL bounds_exch(1,'u+',DTX)
      CALL bounds_exch(1,'v+',DSY)
      CALL bounds_exch(1,'v+',DTY)
      ENDIF ! (iocad==6)

!
      DO K=1,KE
      DO J=1,JE
       DO I=1,IE
        curvnot(i, j) = 0._wp
        curvsot(i, j) = 0._wp
        curvwet(i, j) = 0._wp
        curveat(i, j) = 0._wp
        curvnos(i, j) = 0._wp
        curvsos(i, j) = 0._wp
        curvwes(i, j) = 0._wp
        curveas(i, j) = 0._wp
       ENDDO
      ENDDO
!
      KUP = MAX(1,K-1)
      KLO = MIN(KE,K+1)
!
      KUPUP=MAX(1,K-2)
      KLOLO=MIN(KE,K+2)
      IF(K.EQ.1)THEN
      DO J=1,JE
       DO I=1,IE
         curvupt(i, j) = 0._wp
         curvups(i, j) = 0._wp
         curvlot(i, j) = 0._wp
         curvlos(i, j) = 0._wp
       ENDDO
      ENDDO
      ELSE
      DO J=1,JE
        DO I=1,IE
         CURVUPT(I,J)=CURVLOT(I,J)
         CURVUPS(I,J)=CURVLOS(I,J)
        ENDDO
      ENDDO
      ENDIF

      IF ( iocad == 7 ) THEN
      DO J=1,JE
       DO I=2,IE-1
        curvlot(i, j) = 0._wp
        curvlos(i, j) = 0._wp
        IF (wo(i, j, k+1) .GT. 0._wp)THEN
!
!       UPWARD FLOW
!HH UM SWICH OFF ?
!CWET        IF(K.LE.KE-1.AND.WETO(I,J,KLO).GT.0.5)THEN
!C        IF(K.LE.KE-2.AND.WETO(I,J,KLOLO).GT.0.5)THEN
        CURVLOT(I,J)=((T1O(I,J,K)-T1O(I,J,KLO))/(DZW(K)+DZW(KLO))       &
     &          -(T1O(I,J,KLO)-TINF(I,J,KLO))/(DZW(KLO)+DZW(KLOLO)))    &
     &  *DZW(K)*DZW(KLO)/(DZW(KLOLO)+2._wp*DZW(KLO)+DZW(K))
        CURVLOS(I,J)=((S1O(I,J,K)-S1O(I,J,KLO))/(DZW(K)+DZW(KLO))       &
     &          -(S1O(I,J,KLO)-SINF(I,J,KLO))/(DZW(KLO)+DZW(KLOLO)))    &
     &  *DZW(K)*DZW(KLO)/(DZW(KLOLO)+2._wp*DZW(KLO)+DZW(K))
!CWET        ENDIF
        ELSE
!
!       DOWNWARD FLOW
!
!C        IF(K.GT.1.AND.K.LT.KE)THEN
!C        IF(K.LT.KE)THEN
!C        CURVLOT(I,J)=((T1O(I,J,KUP)-T1O(I,J,K))/(DZW(K)+DZW(KUP))
        CURVLOT(I,J)=((TINF(I,J,K)-T1O(I,J,K))/(DZW(K)+DZW(KUP))        &
     &         -(T1O(I,J,K)-T1O(I,J,KLO))/(DZW(KLO)+DZW(K)))            &
     &  *WETO(I,J,KLO)*DZW(KLO)*DZW(K)/(DZW(KLO)+2._wp*DZW(K)+DZW(KUP))
        CURVLOS(I,J)=((SINF(I,J,K)-S1O(I,J,K))/(DZW(K)+DZW(KUP))        &
     &         -(S1O(I,J,K)-S1O(I,J,KLO))/(DZW(KLO)+DZW(K)))            &
     &  *WETO(I,J,KLO)*DZW(KLO)*DZW(K)/(DZW(KLO)+2._wp*DZW(K)+DZW(KUP))
!C        ENDIF
        ENDIF
       ENDDO
      ENDDO
      ENDIF ! (iocad==7)

      IF ( iocad == 6 ) THEN
      DO J=1,JE
       DO I=2,IE-1
        curvlot(i, j) = 0._wp
        curvlos(i, j) = 0._wp
        IF (wo(i, j, klo) .GE. 0._wp)THEN
!
!       UPWARD FLOW
!
        IF (k .LT. ke-1 .AND. weto(i, j, klolo) .GT. 0.5_wp)THEN
        CURVLOT(I,J)=((T1O(I,J,K)-T1O(I,J,KLO))/(DZW(K)+DZW(KLO))       &
     &          -(T1O(I,J,KLO)-T1O(I,J,KLOLO))/(DZW(KLO)+DZW(KLOLO)))   &
     &  * dzw(k) * dzw(klo) / (dzw(klolo) + 2._wp * dzw(klo) + dzw(k))
        CURVLOS(I,J)=((S1O(I,J,K)-S1O(I,J,KLO))/(DZW(K)+DZW(KLO))       &
     &          -(S1O(I,J,KLO)-S1O(I,J,KLOLO))/(DZW(KLO)+DZW(KLOLO)))   &
     &  * dzw(k) * dzw(klo) / (dzw(klolo) + 2._wp * dzw(klo) + dzw(k))
        ENDIF
        ELSE
!
!       DOWNWARD FLOW
!
!CC        IF(K.GT.1.AND.K.LT.KE)THEN
        CURVLOT(I,J)=((T1O(I,J,KUP)-T1O(I,J,K))/(DZW(K)+DZW(KUP))       &
     &         -(T1O(I,J,K)-T1O(I,J,KLO))/(DZW(KLO)+DZW(K)))            &
     &  * weto(i, j, klo) * dzw(klo) &
     &  * dzw(k) / (dzw(klo) + 2._wp * dzw(k) + dzw(kup))
        CURVLOS(I,J)=((S1O(I,J,KUP)-S1O(I,J,K))/(DZW(K)+DZW(KUP))       &
     &         -(S1O(I,J,K)-S1O(I,J,KLO))/(DZW(KLO)+DZW(K)))            &
     &  * weto(i, j, klo) * dzw(klo) &
     &  * dzw(k) / (dzw(klo) + 2._wp * dzw(k) + dzw(kup))
!CC        ENDIF
        ENDIF
       ENDDO
      ENDDO
      ENDIF ! (iocad==6)
!
      CALL bounds_exch(1,'p',CURVLOT)
      CALL bounds_exch(1,'p',CURVLOS)
!
      IF ( iocad == 7 ) THEN
      DO J=2,JE1
       DO I=2,IE1
!
!       SOUTHERLY POINT
!
        IF (vke(i, j, k) .GT. 0._wp)THEN
!
!        INFLOW FROM SOUTH
!
         CURVSOT(I,J)=                                                  &
     &     ((T1O(I,J,K)-T1O(I,J+1,K))/DLYP2(I,J)                        &
     &      -(T1O(I,J+1,K)-TINF(I,J+1,K))/DLYP2(I,J+1))                 &
     &     *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J)+DLYP2(I,J+1))
!
         CURVSOS(I,J)=                                                  &
     &     ((S1O(I,J,K)-S1O(I,J+1,K))/DLYP2(I,J)                        &
     &      -(S1O(I,J+1,K)-SINF(I,J+1,K))/DLYP2(I,J+1))                 &
     &     *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J)+DLYP2(I,J+1))
!
        ELSE
!
!        OUTFLOW TO SOUTH
!
         CURVSOT(I,J)=                                                  &
     &     ((TINF(I,J,K)-T1O(I,J,K))/DLYP2(I,J-1)                       &
     &      -(T1O(I,J,K)-T1O(I,J+1,K))/DLYP2(I,J))                      &
     &     *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J-1)+DLYP2(I,J))
!
         CURVSOS(I,J)=                                                  &
     &     ((SINF(I,J,K)-S1O(I,J,K))/DLYP2(I,J-1)                       &
     &      -(S1O(I,J,K)-S1O(I,J+1,K))/DLYP2(I,J))                      &
     &     *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J-1)+DLYP2(I,J))
        ENDIF
       ENDDO
      ENDDO

      CALL bounds_exch(1,'v+',CURVSOT)
      call bounds_exch(1,'v+',CURVSOS)

      DO J=2,JE
       DO I=2,IE-1
        CURVNOT(I,J)=CURVSOT(I,J-1)
        CURVNOS(I,J)=CURVSOS(I,J-1)
       ENDDO
      ENDDO

      CALL bounds_exch(1,'v+',CURVNOT)
      call bounds_exch(1,'v+',CURVNOS)

      DO J=2,JE1
       DO I=2,IE1
        IF (uko(i, j, k) .GT. 0._wp) THEN
!
!        INFLOW
!
         CURVEAT(I,J)=                                                  &
     &     ((T1O(I+1,J,K)-T1O(I,J,K))/DLXP2(I,J)                        &
     &      -(T1O(I,J,K)-TINF(I,J,K))/DLXP2(I-1,J))                     &
     &     *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I-1,J)+DLXP2(I,J))
!
         CURVEAS(I,J)=                                                  &
     &     ((S1O(I+1,J,K)-S1O(I,J,K))/DLXP2(I,J)                        &
     &      -(S1O(I,J,K)-SINF(I,J,K))/DLXP2(I-1,J))                     &
     &     *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I-1,J)+DLXP2(I,J))
        ELSE
!
!        OUTFLOW
!
         CURVEAT(I,J)=                                                  &
     &     ((TINF(I+1,J,K)-T1O(I+1,J,K))/DLXP2(I+1,J)                   &
     &      -(T1O(I+1,J,K)-T1O(I,J,K))/DLXP2(I,J))                      &
     &     *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I,J)+DLXP2(I+1,J))
!
         CURVEAS(I,J)=                                                  &
     &     ((SINF(I+1,J,K)-S1O(I+1,J,K))/DLXP2(I+1,J)                   &
     &      -(S1O(I+1,J,K)-S1O(I,J,K))/DLXP2(I,J))                      &
     &     *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I,J)+DLXP2(I+1,J))
        ENDIF
!
       ENDDO
      ENDDO

      CALL bounds_exch(1,'u+',CURVEAT)
      CALL bounds_exch(1,'u+',CURVEAS)

      DO J=1,JE
       DO I=2,IE-1
        CURVWET(I,J)=CURVEAT(I-1,J)
        CURVWES(I,J)=CURVEAS(I-1,J)
       ENDDO
      ENDDO

      CALL bounds_exch(1,'u+',CURVWET)
      CALL bounds_exch(1,'u+',CURVWES)
      ENDIF ! (iocad==7)

      IF ( iocad == 6 ) THEN
      DO J=2,JE1
       DO I=2,IE1
!
!       SOUTHERLY POINT
!
        IF (amsue(i, j, k) .GT. 0.5_wp) THEN
         IF (vke(i, j, k) .GT. 0._wp) THEN
!
!         INFLOW FROM SOUTH
!
          CURVSOT(I,J)=AMSUE(I,J+1,K)*(DTY(I,J,K)-DTY(I,J+1,K))         &
     &                *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J)+DLYP2(I,J+1))
!
          CURVSOS(I,J)=AMSUE(I,J+1,K)*(DSY(I,J,K)-DSY(I,J+1,K))         &
     &                *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J)+DLYP2(I,J+1))
!
         ELSE
!
!         OUTFLOW TO SOUTH
!
          CURVSOT(I,J)=AMSUE(I,J-1,K)*(DTY(I,J-1,K)-DTY(I,J,K))         &
     &                *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J-1)+DLYP2(I,J))
!
          CURVSOS(I,J)=AMSUE(I,J-1,K)*(DSY(I,J-1,K)-DSY(I,J,K))         &
     &                *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J-1)+DLYP2(I,J))
         ENDIF
        ENDIF
       ENDDO
      ENDDO

      CALL bounds_exch(1,'v+',CURVSOT)
      CALL bounds_exch(1,'v+',CURVSOS)

      DO J=2,JE
       DO I=2,IE-1
        CURVNOT(I,J)=CURVSOT(I,J-1)
        CURVNOS(I,J)=CURVSOS(I,J-1)
       ENDDO
      ENDDO

      CALL bounds_exch(1,'v+',CURVNOT)
      CALL bounds_exch(1,'v+',CURVNOS)

      DO J=2,JE1
       DO I=2,IE1
        IF (amsuo(i, j, k) .GT. 0.5_wp)THEN
         IF (uko(i, j, k) .GT. 0._wp) THEN
!
!         INFLOW
!
          CURVEAT(I,J)=AMSUO(I-1,J,K)*(DTX(I,J,K)-DTX(I-1,J,K))         &
     &                *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I-1,J)+DLXP2(I,J))
!
          CURVEAS(I,J)=AMSUO(I-1,J,K)*(DSX(I,J,K)-DSX(I-1,J,K))         &
     &                *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I-1,J)+DLXP2(I,J))
         ELSE
!
!         OUTFLOW
!
          CURVEAT(I,J)=AMSUO(I+1,J,K)*(DTX(I+1,J,K)-DTX(I,J,K))         &
     &                *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I,J)+DLXP2(I+1,J))
!
          CURVEAS(I,J)=AMSUO(I+1,J,K)*(DSX(I+1,J,K)-DSX(I,J,K))         &
     &                *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I,J)+DLXP2(I+1,J))
         ENDIF
        ENDIF
       ENDDO
      ENDDO

      CALL bounds_exch(1,'u+',CURVEAT)
      CALL bounds_exch(1,'u+',CURVEAS)

      DO J=1,JE
       DO I=2,IE-1
        CURVWET(I,J)=CURVEAT(I-1,J)
        CURVWES(I,J)=CURVEAS(I-1,J)
       ENDDO
      ENDDO

      CALL bounds_exch(1,'u+',CURVWET)
      CALL bounds_exch(1,'u+',CURVWES)
      ENDIF ! (iocad==6)
!
!
      DO J=2,JE1
      DO I=2,IE1
!
      TM = T1O(I,J,K)
      SM = S1O(I,J,K)
!
      TA=THO(I,J,K)
      SA=SAO(I,J,K)
      WUN = WO(I,J,K+1)*DLXP(I,J)*DLYP(I,J)                             &
     &     /(DZW(KLO)+DZW(K))
      WOB = WO(I,J,K)*DLXP(I,J)*DLYP(I,J)                               &
     &     /(DZW(KUP)+DZW(K))
      UWE = UKO(I-1,J,K)*DDUO(I-1,J,K)*DLYU(I-1,J)                      &
     &     /(DLXP(I,J)+DLXP(I-1,J))
      UOS = UKO(I,J,K)*DDUO(I,J,K)*DLYU(I,J)                            &
     &     /(DLXP(I,J)+DLXP(I+1,J))
      VSU = VKE(I,J,K)*DDUE(I,J,K)*DLXV(I,J)                            &
     &     /(DLYP(I,J)+DLYP(I,J+1))
      VNO = VKE(I,J-1,K)*DLXV(I,J-1)*DDUE(I,J-1,K)                      &
     &     /(DLYP(I,J)+DLYP(I,J-1))
!
      UK1O(I,J,K) = THO (I,J,K) + DT*WETO(I,J,K)*(                      &
     &         UWE * ( DLXP(I,J)*T1O(I-1,J,K) +DLXP(I-1,J)*TM           &
     &         -(DLXP(I,J)+DLXP(I+1,J))*CURVWET(I,J))                   &
     &     -   UOS * ( DLXP(I+1,J)*TM + DLXP(I,J)*T1O(I+1,J,K)          &
     &         -(DLXP(I,J)+DLXP(I+1,J))*CURVEAT(I,J))                   &
     &     +   VSU * ( DLYP(I,J)*T1O(I,J+1,K) + DLYP(I,J+1)*TM          &
     &         -(DLYP(I,J)+DLYP(I,J-1))*CURVSOT(I,J))                   &
     &     -   VNO * ( DLYP(I,J-1)*TM + DLYP(I,J)*T1O(I,J-1,K)          &
     &         -(DLYP(I,J)+DLYP(I,J-1))*CURVNOT(I,J))                   &
     &     +WUN*(DZW(K)*T1O(I,J,KLO)+DZW(KLO)*TM                        &
     &         -(DZW(KLO)+DZW(K))*CURVLOT(I,J))                         &
     &     -WOB*(DZW(KUP)*TM + DZW(K)*T1O(I,J,KUP)                      &
     &         -(dzw(kup) + dzw(k)) * (1._wp - zsurf(k)) * curvupt(i, j))) &
     &  /((DDPO(I,J,K)+ZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)            &
     &                  -RHOSNWA*SICSNO(I,J))+ALMZER)                   &
     &      *DLXP(I,J)*DLYP(I,J))
!
      VK1E(I,J,K) = SAO (I,J,K) + DT*WETO(I,J,K)*(                      &
     &         UWE * ( DLXP(I,J)*S1O(I-1,J,K) +DLXP(I-1,J)*SM           &
     &         -(DLXP(I,J)+DLXP(I-1,J))*CURVWES(I,J))                   &
     &     -   UOS * ( DLXP(I+1,J)*SM + DLXP(I,J)*S1O(I+1,J,K)          &
     &         -(DLXP(I,J)+DLXP(I+1,J))*CURVEAS(I,J))                   &
     &     +   VSU * ( DLYP(I,J)*S1O(I,J+1,K) + DLYP(I,J+1)*SM          &
     &         -(DLYP(I,J)+DLYP(I,J+1))*CURVSOS(I,J))                   &
     &     -   VNO * ( DLYP(I,J-1)*SM + DLYP(I,J)*S1O(I,J-1,K)          &
     &         -(DLYP(I,J)+DLYP(I,J-1))*CURVNOS(I,J))                   &
     &     +WUN*(DZW(K)*S1O(I,J,KLO)+DZW(KLO)*SM                        &
     &         -(DZW(KLO)+DZW(K))*CURVLOS(I,J))                         &
     &     -WOB*(DZW(KUP)*SM + DZW(K)*S1O(I,J,KUP)                      &
     &         -(DZW(KUP)+DZW(K))*CURVUPS(I,J)))                        &
     & /((DDPO(I,J,K)+ZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)             &
     &       -RHOSNWA*SICSNO(I,J))+ALMZER)                              &
     &       *DLXP(I,J)*DLYP(I,J))
!
      END DO
      END DO
      END DO
!
      CALL bounds_exch(1,'p',UK1O)
      CALL bounds_exch(1,'p',VK1E)
!
       DO K=1,KE
        DO J=1,JE
         DO I=1,IE
          THO(I,J,K)=UK1O(I,J,K)
          SAO(I,J,K)=VK1E(I,J,K)
         ENDDO
        ENDDO
       ENDDO
!
      RETURN
      END
