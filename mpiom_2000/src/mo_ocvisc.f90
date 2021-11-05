MODULE mo_ocvisc

  USE MO_PARAM1
  USE MO_PARALLEL
  USE MO_COMMO1
  USE mo_boundsexch, ONLY: bounds_exch

  IMPLICIT NONE

  PRIVATE

  REAL(wp) :: bofric !<     BOFRIC  FRICTION COEFFICENT AT BOTTOM [M**2/S]
  REAL(wp) :: rayfric !<     COEFFICIENT FOR QUADRATIC BOTTOM FRICTION []

  PUBLIC :: bofric, rayfric, ocvisc

CONTAINS

  SUBROUTINE OCVISC

!-=====================================================================
!
!     SBR OCVISC
!
!     COMPUTES MOMENTUM DIFFUSION
!     VERTICALLY (IMPLICITE)
!     RAYLEIGH BOTTOM FRICTION
!     HORIZONTALLY BIHARMONIC
!
!     LAST MODIFIED UWE MIKOLAJEWICZ 12.12.99
!     RAYLEIGH FRICTION NOW QUADRATIC!!!!
!     DX,DY INCLUDED IN BIHARMONIC FRICTION!!
!
!     MODIFIED 4.2.2000 UWE MIKOLAJEWICZ
!      INCLUDE QUADRATIC BOTTOM FRICTION IN VERTICAL DIFFUSION,
!              REMOVE ERROR
!
!
!     BIHARMONIC DIFFUSION COEFFICIENTS, DX,DY ON PSI-POINTS
!

    INTEGER i,j,k,l
    REAL(wp) AULUX,AULUY,AULVX,AULVY,SPEEDU,SPEEDV
    DIMENSION AULUX(IE,JE),AULUY(IE,JE),AULVX(IE,JE),AULVY(IE,JE)     &
         &        ,SPEEDU(IE,JE),SPEEDV(IE,JE)
    REAL(wp) TRIDSY(ie,ke,3)
    REAL(wp) avup,avlo
!
!OtB_TURBSHEAR ifdef TURBSHEAR
!OtB_TURBSHEAR
!OtB_TURBSHEAR      COEFFICIENT FOR SHEAR DIFFUSION
!OtB_TURBSHEAR
!OtB_TURBSHEAR      DIMENSION TURBE(IE,JE)
!OtB_TURBSHEAR endif
!
!     VERTICAL DIFFUSION OF MOMENTUM
!
!
!
!     SPEED IN BOTTOM LAYER AT OLD TIME STEP
!       ON U AND V POINTS
!

!$OMP PARALLEL PRIVATE(i,j,k,l,avup,avlo,tridsy)

!$OMP DO
      DO J=1,JE
       DO I=1,IE
        speedu(i, j) = 0._wp
        speedv(i, j) = 0._wp
       ENDDO
      ENDDO
!
!$OMP DO
      DO J=2,JE-1
       DO K=1,KE-1
        DO I=2,IE1
         IF(AMSUO(I,J,K)-AMSUO(I,J,K+1).GT.0.5_wp) THEN
          SPEEDU(I,J)=SQRT(UOO(I,J,K)**2+                               &
     &      (0.25_wp * (voe(i, j, k) + voe(i+1, j, k)                   &
     &      +VOE(I,J-1,K)+VOE(I+1,J-1,K)))**2)
          speedu(i, j) = MAX(speedu(i, j), 1.E-3_wp)
         ENDIF
         IF (amsue(i, j, k) - amsue(i, j, k+1).GT.0.5_wp) THEN
          speedv(i, j) = SQRT(voe(i, j, k)**2 +                         &
     &      (0.25_wp * (uoo(i, j, k) + uoo(i-1, j, k)                   &
     &      + uoo(i, j+1, k) + uoo(i-1, j+1, k)))**2)
          speedv(i, j) = MAX(speedv(i, j), 1.E-3_wp)
         ENDIF
        ENDDO
       ENDDO
       K=KE
       DO I=2,IE-1
        IF (amsuo(i, j, k) .GT. 0.5_wp) THEN
          speedu(i, j) = SQRT(uoo(i, j, k)**2 +                         &
     &      (0.25_wp * (voe(i, j, k) + voe(i+1, j, k)                   &
     &      + voe(i, j-1, k) + voe(i+1, j-1, k)))**2)
          speedu(i, j) = MAX(speedu(i, j), 1.E-3_wp)
         ENDIF
         IF (amsue(i, j, k) .GT. 0.5_wp) THEN
          speedv(i, j) = SQRT(voe(i, j, k)**2 +                         &
     &      (0.25_wp * (uoo(i, j, k) + uoo(i-1, j, k)                   &
     &      + uoo(i, j+1, k) + uoo(i-1, j+1, k)))**2)
          speedv(i, j) = MAX(speedv(i, j), 1.E-3_wp)
         ENDIF
       ENDDO
      ENDDO
!
!$OMP SINGLE
      CALL bounds_exch(1,'u',SPEEDU,'ocvisc 1')
      CALL bounds_exch(1,'v',SPEEDV,'ocvisc 2')
!$OMP END SINGLE
!
!     OUTER LOOP
!
!$OMP DO
      DO J=2,JE1
!
      DO K=1,KE-1
        DO I=1,IE
          avup = 0.5_wp * (avo(i, j, k) + avo(i, j+1, k)) * amsue(i, j, k)
          avlo = 0.5_wp * (avo(i, j, k+1) + avo(i, j+1, k+1)) &
               * amsue(i, j, k+1) + (bofric + rayfric * ddue(i, j, k) &
               &                     * speedv(i, j))          &
               * (amsue(i, j, k) - amsue(i, j, k+1))
          !       TRIDSY(I,K,1)= - DT * AVUP   * DWIE(I,J,K) * DI(K)
          !       TRIDSY(I,K,3)= - DT * AVLO * DWIE(I,J,K) * DI(K+1)
          TRIDSY(I,K,1)= - DT * AVUP   * AMSUE(I,J,K) * DI(K)             &
               &          /(ALMZER+DDUE(I,J,K))
          TRIDSY(I,K,3)= - DT * AVLO * AMSUE(I,J,K) * DI(K+1)             &
               &          /(ALMZER+DDUE(I,J,K))
          tridsy(i, k, 2) = 1._wp - tridsy(i, k, 1) - tridsy(i, k, 3)
          TRIDSY(I,K,3)=TRIDSY(I,K,3)*AMSUE(I,J,K+1)
        END DO
      END DO
!
      K=KE
!
      DO I=1,IE
        avup = 0.5_wp * (avo(i, j, k) + avo(i, j+1, k)) * amsue(i, j, k)
        AVLO=(BOFRIC+RAYFRIC*                                             &
     &     DDUE(I,J,K)*SPEEDV(I,J))*AMSUE(I,J,K)
!     TRIDSY(I,K,1)= - DT * AVUP   * DWIE(I,J,K) * DI(K)
!     TRIDSY(I,K,3)= - DT * AVLO * DWIE(I,J,K) * DI(K+1)
        TRIDSY(I,K,1)= - DT * AVUP   * AMSUE(I,J,K) * DI(K)               &
     &        /(ALMZER+DDUE(I,J,K))
        TRIDSY(I,K,3)= - DT * AVLO * AMSUE(I,J,K) * DI(K+1)               &
     &        /(ALMZER+DDUE(I,J,K))

        tridsy(i, k, 2) = 1._wp - tridsy(i, k, 1) - tridsy(i, k, 3)
        tridsy(i, k, 3) = 0._wp
      END DO
!
      DO K=2,KE
        DO I=1,IE
          TRIDSY(I,K-1,1)  = TRIDSY(I,K,1) / TRIDSY(I,K-1,2)
          TRIDSY(I,K,2)    = TRIDSY(I,K,2) -                                &
               &   TRIDSY(I,K-1,3) * TRIDSY(I,K,1) / TRIDSY(I,K-1,2)
        END DO
      END DO
!
      DO K=1,KE
        DO I=1,IE
          VK1E(I,J,K)=VOE(I,J,K)
          UK1O(I,J,K)=UOO(I,J,K)
        END DO
      END DO
!
      DO K=2,KE
        DO I=1,IE
          VK1E(I,J,K)=VK1E(I,J,K)-TRIDSY(I,K-1,1)*VK1E(I,J,K-1)
        END DO
      END DO
!
      K=KE
      DO I=1,IE
        VKE(I,J,K)=AMSUE(I,J,K)*  VK1E(I,J,K) / TRIDSY(I,K,2)
      END DO
!
!
      DO K=1,KE1
        L=KE-K
        DO I=1,IE
          VKE(I,J,L) = ( VK1E(I,J,L) - TRIDSY(I,L,3) * VKE(I,J,L+1) )       &
               &*AMSUE(I,J,L) / TRIDSY(I,L,2)
        END DO
      END DO
!
!     ODD FIELDS
!
      DO K=1,KE-1
        DO I=1,IE1
          avup = 0.5_wp * (avo(i, j, k) + avo(i+1, j, k)) * amsuo(i, j, k)
          avlo = 0.5_wp * (avo(i, j, k+1) + avo(i+1, j, k)) * amsuo(i, j, k+1) &
     &+(BOFRIC+RAYFRIC*DDUO(I,J,K)*SPEEDU(I,J))                         &
     &    *(AMSUO(I,J,K)-AMSUO(I,J,K+1))
!     TRIDSY(I,K,1)= - DT * AVUP   * DWIO(I,J,K) * DI(K)
!     TRIDSY(I,K,3)= - DT * AVLO * DWIO(I,J,K) * DI(K+1)
          TRIDSY(I,K,1)= - DT * AVUP   * AMSUO(I,J,K)* DI(K)                &
               &        /(almzer+dduo(i,j,k))
          TRIDSY(I,K,3)= - DT * AVLO * amsuo(I,J,K) * DI(K+1)               &
               &         /(almzer+dduo(i,j,k))
          tridsy(i, k, 2) = 1._wp - tridsy(i, k, 1) - tridsy(i, k, 3)
          TRIDSY(I,K,3)=TRIDSY(I,K,3)*AMSUO(I,J,K+1)
        END DO
      END DO
!     TRIDSY(IE,K,1)=TRIDSY(2,K,1)
!     TRIDSY(IE,K,2)=TRIDSY(2,K,2)
!     TRIDSY(IE,K,3)=TRIDSY(2,K,3)
!
      K=KE
      DO I=1,IE1
        avup = 0.5_wp * (avo(i, j, k) + avo(i+1, j, k)) * amsuo(i, j, k)
        AVLO=(BOFRIC+RAYFRIC                                              &
     &     *DDUO(I,J,K)*SPEEDU(I,J))*AMSUO(I,J,K)
!     TRIDSY(I,K,1)= - DT * AVUP   * DWIO(I,J,K) * DI(K)
!     TRIDSY(I,K,3)= - DT * AVLO * DWIO(I,J,K) * DI(K+1)
        TRIDSY(I,K,1)= - DT * AVUP   * amsuO(I,J,K) * DI(K)               &
     &         /(almzer+dduo(i,j,k))
        TRIDSY(I,K,3)= - DT * AVLO * amsuO(I,J,K) * DI(K+1)               &
     &         /(almzer+dduo(i,j,k))

        tridsy(i, k, 2) = 1._wp - tridsy(i, k, 1) - tridsy(i, k, 3)
        tridsy(i, k, 3) = 0._wp
      END DO
!
      DO K=2,KE
        DO I=1,IE1
          TRIDSY(I,K-1,1) = TRIDSY(I,K,1) / TRIDSY(I,K-1,2)
          TRIDSY(I,K,2)   = TRIDSY(I,K,2) -                                 &
               &      TRIDSY(I,K-1,3) * TRIDSY(I,K,1) / TRIDSY(I,K-1,2)
        END DO
      END DO
!
      DO K=2,KE
        DO I=1,IE1
          UK1O(I,J,K) = UK1O(I,J,K) - TRIDSY(I,K-1,1) * UK1O(I,J,K-1)
        END DO
      END DO
!
      K=KE
!
      DO I=1,IE1
        UKO(I,J,K) =AMSUO(I,J,K)*   UK1O(I,J,K) / TRIDSY(I,K,2)
      END DO
!
      DO K=1,KE1
        L=KE-K
        DO I=1,IE
          UKO(I,J,L) = ( UK1O(I,J,L) - TRIDSY(I,L,3) * UKO(I,J,L+1) )       &
               &*AMSUO(I,J,L) / TRIDSY(I,L,2)
        END DO
      END DO
!
      END DO
!
!      END VERTICAL DIFFUSION
!
!#ifdef bounds_exch_save
!$OMP SINGLE
      CALL bounds_exch(1,'u',UKO,'ocvisc 3')
      CALL bounds_exch(1,'v',VKE,'ocvisc 4')
!$OMP END SINGLE
!#endif
!
!
!     HORIZONTAL DIFFUSION OF MOMENTUM
!
!OtB_TURBSHEAR ifdef TURBSHEAR
!OtB_TURBSHEAR
!OtB_TURBSHEAR      WRITE(IO_STDOUT,*)'SHEAR-TURBULENCE, NOT MASS CONSERVING!!!!'
!OtB_TURBSHEAR      DO 576 J=1,JE
!OtB_TURBSHEAR       DO 576 I=1,IE
!OtB_TURBSHEAR        TURBE(I,J)=0.
!OtB_TURBSHEAR 576  CONTINUE
!OtB_TURBSHEAR
!OtB_TURBSHEAR endif
!
!     COEFFICIENTS FOR BIHARMONIC HORIZONTAL MOMENTUM DIFFUSION,
!         DEFINED ON P- AND PSI-POINTS
!
!$OMP DO
      DO J=1,JE
       DO I=2,IE1

!OtB_AULREDUV ifdef AULREDUV
!OtB_AULREDUV        AULUX(I,J)=AULAPUV*1.E5*DLXP(I,J)**2*MAX(DLXP(I,J),DLYP(I,J))
!OtB_AULREDUV        AULVY(I,J)=AULAPUV*1.E5*DLYP(I,J)**2*MAX(DLXP(I,J),DLYP(I,J))
!OtB_AULREDUV        AULUY(I,J)=AULAPUV*1.E5*DLYPSI(I,J)**2                          &
!OtB_AULREDUV     &              *MAX(DLXPSI(I,J),DLYPSI(I,J))
!OtB_AULREDUV        AULVX(I,J)=AULAPUV*1.E5*DLXPSI(I,J)**2                          &
!OtB_AULREDUV     &              *MAX(DLXPSI(I,J),DLYPSI(I,J))
!OtB_AULREDUV else
#ifdef AULREDUV
chh        fix for anisotropic grid in the eq refinement
         AULUX(I,J)=AULAPUV*MAX(DLXP(I,J),1.)**3
     x   *min(dlxp(i,j),dlyp(i,j))
        AULVY(I,J)=AULAPUV*MAX(1.,DLYP(I,J))**3
     x   *min(dlxp(i,j),dlyp(i,j))
        AULUY(I,J)=AULAPUV*MAX(1.,DLYPSI(I,J))**3
     x    *min(dlxpsi(i,j),dlypsi(i,j))
        AULVX(I,J)=AULAPUV*MAX(DLXPSI(I,J),1.)**3
     x     *min(dlxpsi(i,j),dlypsi(i,j))
#else
        aulux(i, j) = aulapuv * MAX(dlxp(i, j), 1._wp)**4
        aulvy(i, j) = aulapuv * MAX(1._wp, dlyp(i, j))**4
        auluy(i, j) = aulapuv * MAX(1._wp, dlypsi(i, j))**4
        aulvx(i, j) = aulapuv * MAX(dlxpsi(i, j), 1._wp)**4
#endif /*AULREDUV*/
       ENDDO
      ENDDO
!
!$OMP SINGLE


!#ifdef bounds_exch_save
      CALL bounds_exch(1,'p',AULUX,'ocvisc 7')
      CALL bounds_exch(1,'p',AULUY,'ocvisc 8')
      CALL bounds_exch(1,'s',AULVX,'ocvisc 9')
      CALL bounds_exch(1,'s',AULVY,'ocvisc 10')
!#endif
!$OMP END SINGLE
!
!     MAIN LOOP OVER LEVELS
!
!
!OtB_TURBSHEAR ifdef TURBSHEAR
!OtB_TURBSHEAR
!OtB_TURBSHEAR      SHEAR FRICTION
!OtB_TURBSHEAR
!OtB_TURBSHEAR      TURNED OFF!!!!!
!OtB_TURBSHEAR
!OtB_TURBSHEAR      DO 1501 J=2,JE1
!OtB_TURBSHEAR      DO 1501 I=2,IE1
!OtB_TURBSHEAR      SHEAR=VKE(I+1,J,K)-VKE(I,J,K)                                     &
!OtB_TURBSHEAR     & +(DTDXPE(I,J)/DPYE(I,J))                                         &
!OtB_TURBSHEAR     & *(UKO(I,J,K)-UKO(I,J+1,K))
!OtB_TURBSHEAR      MIT 0.03 GINEG ES EINIGE MONATE GUT, ABER GOLFSTROM MAU
!OtB_TURBSHEAR      TURPOW=SHEAR**4
!OtB_TURBSHEAR UWE   STILLGELEGT
!OtB_TURBSHEAR C      TURBE(I,J)=0.001*TURPOW/(1.+TURPOW)
!OtB_TURBSHEAR      TURBE(I,J)=0.
!OtB_TURBSHEAR 501  CONTINUE
!OtB_TURBSHEAR      IF(ICYCLI.GE.1) CALL PERIO2(881,TURBE)
!OtB_TURBSHEAR      DO 1502 J=2,JE1
!OtB_TURBSHEAR      DO 1502 I=2,IE1
!OtB_TURBSHEAR      S1O(I,J,K)=WETO(I,J,K)*(                                          &
!OtB_TURBSHEAR     &WETO(I+1,J,K)*(TURBE(I,J)+TURBE(I,J-1))*(SAO(I+1,J,K)             &
!OtB_TURBSHEAR     & - SAO(I,J,K)) +WETO(I-1,J,K)*(TURBE(I-1,J)+TURBE(I-1,J-1))       &
!OtB_TURBSHEAR     & *(SAO(I-1,J,K)-SAO(I,J,K))  )
!OtB_TURBSHEAR      T1O(I,J,K)=WETO(I,J,K)*(                                          &
!OtB_TURBSHEAR     &WETO(I+1,J,K)*(TURBE(I,J)+TURBE(I,J-1))*(THO(I+1,J,K)             &
!OtB_TURBSHEAR     & - THO(I,J,K)) +WETO(I-1,J,K)*(TURBE(I-1,J)+TURBE(I-1,J-1))       &
!OtB_TURBSHEAR     & *(THO(I-1,J,K)-THO(I,J,K))  )
!OtB_TURBSHEAR      UK1O(I,J,K)=UKO(I,J,K)                                            &
!OtB_TURBSHEAR     & +TURBE(I,J)*(UKO(I,J+1,K)-UKO(I,J,K))                            &
!OtB_TURBSHEAR     & +TURBE(I,J-1)*(UKO(I,J-1,K)-UKO(I,J,K))
!OtB_TURBSHEAR      4 -1.E-2*(2.*AMSUO(I,J,K)-AMSUO(I-1,J,K)-AMSUO(I+1,J,K))
!OtB_TURBSHEAR      5 *UKO(I,J,K)
!OtB_TURBSHEAR      VK1E(I,J,K)=VKE(I,J,K)                                            &
!OtB_TURBSHEAR     & +TURBE(I,J)*(VKE(I+1,J,K)-VKE(I,J,K))                            &
!OtB_TURBSHEAR     & +TURBE(I-1,J)*(VKE(I-1,J,K)-VKE(I,J,K))
!OtB_TURBSHEAR      4 -1.E-2*(AMSUE(I,J,K)*2.-AMSUE(I,J-1,K)-AMSUE(I,J+1,K))
!OtB_TURBSHEAR      5 *VKE(I,J,K)
!OtB_TURBSHEAR 502  CONTINUE
!OtB_TURBSHEAR      IF(ICYCLI.GE.1) THEN
!OtB_TURBSHEAR      CALL PERIO32(S1O,K)
!OtB_TURBSHEAR      CALL PERIO32(T1O,K)
!OtB_TURBSHEAR      CALL PERIO32(UK1O,K)
!OtB_TURBSHEAR      CALL PERIO32(VK1E,K)
!OtB_TURBSHEAR      ENDIF
!OtB_TURBSHEAR
!OtB_TURBSHEAR      DO 1527 J=1,JE
!OtB_TURBSHEAR      DO 1527 I=1,IE
!OtB_TURBSHEAR      THO(I,J,K)=THO(I,J,K)+T1O(I,J,K)*WETO(I,J,K)
!OtB_TURBSHEAR      SAO(I,J,K)=SAO(I,J,K)+S1O(I,J,K)*WETO(I,J,K)
!OtB_TURBSHEAR      UKO(I,J,K)=UK1O(I,J,K)*AMSUO(I,J,K)
!OtB_TURBSHEAR      VKE(I,J,K)=VK1E(I,J,K)*AMSUE(I,J,K)
!OtB_TURBSHEAR C      UKO(I,J,K)=UK1O(I,J,K)*DDUO(I,J,K)*DLYU(I,J)
!OtB_TURBSHEAR C      VKE(I,J,K)=VK1E(I,J,K)*DDUE(I,J,K)*DLXV(I,J)
!OtB_TURBSHEAR 527  CONTINUE
!OtB_TURBSHEAR 679  CONTINUE
!OtB_TURBSHEAR      IF(ICYCLI.GE.1) THEN
!OtB_TURBSHEAR      CALL PERIO3(115,THO)
!OtB_TURBSHEAR      CALL PERIO3(115,SAO)
!OtB_TURBSHEAR      CALL PERIO3(115,UKO)
!OtB_TURBSHEAR      CALL PERIO3(115,VKE)
!OtB_TURBSHEAR      ENDIF
!OtB_TURBSHEAR endif
!
#ifdef FREESLIP
!     BIHARMONIC HORIZONTAL MOMENTUM DIFFUSION
!
!$OMP DO
      DO K=1,KE
      DO J=2,JE1
      DO I=2,IE1
        UK1O(I,J,K)=AMSUO(I,J,K)*                                       &
     &((AMSUO(I+1,J,K)*(UKO(I+1,J,K)-UKO(I,J,K))/DLXP(I+1,J)            &
     &-AMSUO(I-1,J,K)*(UKO(I,J,K)-UKO(I-1,J,K))/DLXP(I,J))/DLXU(I,J)    &
     &+(AMSUO(I,J-1,K)*(UKO(I,J-1,K)-UKO(I,J,K))/DLYPSI(I,J-1)          &
     &-AMSUO(I,J+1,K)*(UKO(I,J,K)-UKO(I,J+1,K))/DLYPSI(I,J))            &
     &/DLYU(I,J))
!
        VK1E(I,J,K)=AMSUE(I,J,K)*                                       &
     &((AMSUE(I+1,J,K)*(VKE(I+1,J,K)-VKE(I,J,K))/DLXPSI(I,J)            &
     &-AMSUE(I-1,J,K)*(VKE(I,J,K)-VKE(I-1,J,K))/DLXPSI(I-1,J))          &
     &/DLXV(I,J)                                                        &
     &+(AMSUE(I,J-1,K)*(VKE(I,J-1,K)-VKE(I,J,K))/DLYP(I,J)              &
     &-AMSUE(I,J+1,K)*(VKE(I,J,K)-VKE(I,J+1,K))/DLYP(I,J+1))/DLYV(I,J))
!
      ENDDO
      ENDDO
      ENDDO
!$OMP SINGLE
      CALL bounds_exch(1,'u',UK1O,'ocvisc 11')
      CALL bounds_exch(1,'v',VK1E,'ocvisc 12')
!$OMP END SINGLE
!
!$OMP DO
      DO K=1,KE
      DO J=2,JE1
      DO I=2,IE1
        UOO(I,J,K)=(AMSUO(I,J,K)/(DLXU(I,J)*DLYU(I,J)))                 &
     &            *(DLXU(I,J)*DLYU(I,J)*UKO(I,J,K)                      &
     &             -(AULUX(I+1,J)*DLYP(I+1,J)*AMSUO(I+1,J,K)            &
     &                  *(UK1O(I+1,J,K)-UK1O(I,J,K))/DLXP(I+1,J)        &
     &              -AULUX(I,J)*DLYP(I,J)*AMSUO(I-1,J,K)                &
     &                  *(UK1O(I,J,K)-UK1O(I-1,J,K))/DLXP(I,J)          &
     &              +AULUY(I,J-1)*DLXPSI(I,J-1)*AMSUO(I,J-1,K)          &
     &                  *(UK1O(I,J-1,K)-UK1O(I,J,K))/DLYPSI(I,J-1)      &
     &              -AULUY(I,J)*DLXPSI(I,J)*AMSUO(I,J+1,K)              &
     &                  *(UK1O(I,J,K)-UK1O(I,J+1,K))/DLYPSI(I,J)))
!
        VOE(I,J,K)=(AMSUE(I,J,K)/(DLXV(I,J)*DLYV(I,J)))                 &
     &            *(DLXV(I,J)*DLYV(I,J)*VKE(I,J,K)                      &
     &             -(AULVX(I,J)*DLYPSI(I,J)*AMSUE(I+1,J,K)              &
     &                  *(VK1E(I+1,J,K)-VK1E(I,J,K))/DLXPSI(I,J)        &
     &              -AULVX(I-1,J)*DLYPSI(I-1,J)*AMSUE(I-1,J,K)          &
     &                  *(VK1E(I,J,K)-VK1E(I-1,J,K))/DLXPSI(I-1,J)      &
     &              +AULVY(I,J)*DLXP(I,J)*AMSUE(I,J-1,K)                &
     &                  *(VK1E(I,J-1,K)-VK1E(I,J,K))/DLYP(I,J)          &
     &              -AULVY(I,J+1)*DLXP(I,J+1)*AMSUE(I,J+1,K)            &
     &                  *(VK1E(I,J,K)-VK1E(I,J+1,K))/DLYP(I,J+1)))
!
      ENDDO
      ENDDO
      ENDDO
!
#else
!$OMP DO
      DO K=1,KE
      DO J=2,JE1
      DO I=2,IE1
        UK1O(I,J,K)=AMSUO(I,J,K)*                                       &
     &             (((UKO(I+1,J,K)-UKO(I,J,K))/DLXP(I+1,J)              &
     &              -(UKO(I,J,K)-UKO(I-1,J,K))/DLXP(I,J))/DLXU(I,J)     &
     &             +((UKO(I,J-1,K)-UKO(I,J,K))/DLYPSI(I,J-1)            &
     &              -(UKO(I,J,K)-UKO(I,J+1,K))/DLYPSI(I,J))/DLYU(I,J))
!
        VK1E(I,J,K)=AMSUE(I,J,K)*                                       &
     &             (((VKE(I+1,J,K)-VKE(I,J,K))/DLXPSI(I,J)              &
     &              -(VKE(I,J,K)-VKE(I-1,J,K))/DLXPSI(I-1,J))/DLXV(I,J) &
     &             +((VKE(I,J-1,K)-VKE(I,J,K))/DLYP(I,J)                &
     &              -(VKE(I,J,K)-VKE(I,J+1,K))/DLYP(I,J+1))/DLYV(I,J))
!
      ENDDO
      ENDDO
      ENDDO
!#ifdef bounds_exch_save
!$OMP SINGLE
      CALL bounds_exch(1,'u',UK1O,'ocvisc 13')
      CALL bounds_exch(1,'v',VK1E,'ocvisc 14')
!$OMP END SINGLE
!#endif
!
!$OMP DO
      DO K=1,KE
      DO J=2,JE1
      DO I=2,IE1
        UOO(I,J,K)=(AMSUO(I,J,K)/(DLXU(I,J)*DLYU(I,J)))                 &
     &            *(DLXU(I,J)*DLYU(I,J)*UKO(I,J,K)                      &
     &             -(AULUX(I+1,J)*DLYP(I+1,J)                           &
     &                  *(UK1O(I+1,J,K)-UK1O(I,J,K))/DLXP(I+1,J)        &
     &              -AULUX(I,J)*DLYP(I,J)                               &
     &                  *(UK1O(I,J,K)-UK1O(I-1,J,K))/DLXP(I,J)          &
     &              +AULUY(I,J-1)*DLXPSI(I,J-1)                         &
     &                  *(UK1O(I,J-1,K)-UK1O(I,J,K))/DLYPSI(I,J-1)      &
     &              -AULUY(I,J)*DLXPSI(I,J)                             &
     &                  *(UK1O(I,J,K)-UK1O(I,J+1,K))/DLYPSI(I,J)))
!
        VOE(I,J,K)=(AMSUE(I,J,K)/(DLXV(I,J)*DLYV(I,J)))                 &
     &            *(DLXV(I,J)*DLYV(I,J)*VKE(I,J,K)                      &
     &             -(AULVX(I,J)*DLYPSI(I,J)                             &
     &                  *(VK1E(I+1,J,K)-VK1E(I,J,K))/DLXPSI(I,J)        &
     &              -AULVX(I-1,J)*DLYPSI(I-1,J)                         &
     &                  *(VK1E(I,J,K)-VK1E(I-1,J,K))/DLXPSI(I-1,J)      &
     &              +AULVY(I,J)*DLXP(I,J)                               &
     &                  *(VK1E(I,J-1,K)-VK1E(I,J,K))/DLYP(I,J)          &
     &              -AULVY(I,J+1)*DLXP(I,J+1)                           &
     &                  *(VK1E(I,J,K)-VK1E(I,J+1,K))/DLYP(I,J+1)))
!
      ENDDO
      ENDDO
      ENDDO
#endif /*FREESLIP*/
!#ifdef bounds_exch_save
!$OMP SINGLE
      CALL bounds_exch(1,'u',UOO,'ocvisc 15')
      CALL bounds_exch(1,'v',VOE,'ocvisc 16')
!$OMP END SINGLE
!#endif
!
!$OMP DO
      DO K=1,KE
        DO J=1,JE
          DO I=1,IE
            UKO(I,J,K)=UOO(I,J,K)
            VKE(I,J,K)=VOE(I,J,K)
          END DO
        END DO
      END DO
!
!$OMP END PARALLEL
!

    END SUBROUTINE OCVISC
  END MODULE mo_ocvisc
