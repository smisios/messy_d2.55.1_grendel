SUBROUTINE OCVAD(TRF)

  USE mo_kind, ONLY : dp, i4, wp
  USE mo_param1, ONLY : ie, je, ke, ie1, je1, kep
  USE mo_parallel, ONLY : p_joff
  USE mo_boundsexch, ONLY : bounds_exch
  USE mo_commo1, ONLY : dt, iocaduv, dlxv, dlyv, dlypsi, dlxp, &
                        dduo, ddue, amsue, uko, vke, lbounds_exch_tp

  IMPLICIT NONE

  REAL(dp), INTENT(INOUT) :: TRF(IE,JE,KE)
  REAL(dp) :: TRP(IE,JE,KEP),TRM(IE,JE,0:KE)
  REAL(dp) :: TRVOL(IE,JE,KE)
  REAL(dp) :: WTP(IE,JE,KEP),WTM(IE,JE,0:KE)

  REAL(dp) :: abl,zwabl,up,um,scha,sch,schi,rn,uos,uwe,vsu,vno, &
              tf, wf
  INTEGER(i4) :: i, j, k

!$OMP PARALLEL PRIVATE(i,j,k,abl,zwabl,up,um,scha,sch,schi,rn,uos,uwe,vsu,vno)

!$OMP DO
  DO K=1,KE
    DO J=1,JE
      DO I=1,IE
        TRVOL(I,J,K)=DLXV(I,J)*DLYV(I,J)*DDUE(I,J,K)
        trp(i, j, k) = 0._wp
        trm(i, j, k) = 0._wp
        wtp(i, j, k) = 0._wp
        wtm(i, j, k) = 0._wp
      ENDDO
    ENDDO
  ENDDO

! TRANSPORTS IN Z-DIRECTION

  CALL ocvad_z(trf, trvol, wtm, wtp, trm, trp)

! TRANSPORTS IN X-DIRECTION
!$OMP DO
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        ABL=ABS(TRF(I+1,J,K)-TRF(I-1,J,K))
        zwabl = ABS(trf(i+1, j, k) + trf(i-1, j, k) - 2._wp * trf(i, j, k))
        uos = 0.5_wp * (uko(i, j, k) * dduo(i, j, k) &
             &          + uko(i, j+1, k) * dduo(i, j+1, k))
!EMR    UWE=0.5*(UKO(I-1,J,K)*DDUO(I,J,K)+UKO(I,J+1,K)*DDUO(I,J+1,K))
        uwe = 0.5_wp * (uko(i-1, j, k) * dduo(i-1, j, k) &
             &          + uko(i-1, j+1, k) * dduo(i-1, j+1, k))
        up = 0.5_wp * dt * dlypsi(i, j) * (uos + ABS(uos))
        um = 0.5_wp * dt * (ABS(uwe) - uwe) * dlypsi(i-1, j)
        scha = MAX(0._wp, (abl - zwabl) / (abl + 1.E-20_wp))
        !
        sch = MIN(1._wp, scha * trvol(i, j, k) / (up + um + 1.E-20_wp))
        schi = 1._wp - sch
        !
        trp(i, j, k) = up * (schi * trf(i, j, k) &
             + sch * 0.5_wp * (trf(i, j, k) + trf(i+1, j, k)))
        trm(i, j, k) = um * (schi * trf(i, j, k) &
             + sch * 0.5_wp * (trf(i, j, k) + trf(i-1, j, k)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      ENDDO
    ENDDO
  ENDDO
  !
!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch(1,'uu',TRP,'ocvad 1')
  CALL bounds_exch(1,'uu',TRM,'ocvad 2')
  CALL bounds_exch(1,'uu',WTP,'ocvad 3')
  CALL bounds_exch(1,'uu',WTM,'ocvad 4')
!$OMP END SINGLE
!#endif
  !
!$OMP DO
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        IF (amsue(i, j, k) .GT. 0.5_wp) THEN
          !
          RN=TRVOL(I,J,K)+WTP(I-1,J,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I+1,J,K)
          TRF(I,J,K)=(TRF(I,J,K)*TRVOL(I,J,K)+TRP(I-1,J,K)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I+1,J,K))/RN
          TRVOL(I,J,K)=RN
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !
!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch(1,'v+',TRVOL,'ocvad 5')
  CALL bounds_exch(1,'v',TRF,'ocvad 6')
!$OMP END SINGLE
!#endif

! TRANSPORTS IN Z-DIRECTION (for fractional advection with 1/3 of time step)
  IF (iocaduv .EQ. 8)  CALL ocvad_z(trf, trvol, wtm, wtp, trm, trp)

! TRANSPORTS IN Y-DIRECTION

!$OMP DO
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        ABL=ABS(TRF(I,J-1,K)-TRF(I,J+1,K))
        zwabl = ABS(trf(i, j-1, k) + trf(i, j+1, k) - 2._wp * trf(i, j, k))
        vno = 0.5_wp * (vke(i, j, k) * ddue(i, j, k) &
             + vke(i, j-1, k) * ddue(i, j-1, k))
        vsu = 0.5_wp * (vke(i, j, k) * ddue(i, j, k) &
             + vke(i, j+1, k) * ddue(i, j+1, k))
        up = 0.5_wp * dt * dlxp(i, j) * (vno + ABS(vno))
        um = 0.5_wp * dt * dlxp(i, j+1) * (ABS(vsu) - vsu)

        !
        scha = MAX(0._wp, (abl - zwabl) / (abl + 1.E-20_wp))
        !
        sch = MIN(1._wp, scha * trvol(i, j, k) / (up + um + 1.E-20_wp))
        schi = 1._wp - sch
        trp(i, j, k) = up * (schi * trf(i, j, k) &
             + sch * 0.5_wp * (trf(i, j, k) + trf(i, j-1, k)))
        trm(i, j, k) = um * (schi * trf(i, j, k) &
             + sch * 0.5_wp * (trf(i, j, k) + trf(i, j+1, k)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      ENDDO
    ENDDO
  ENDDO
  !
!$OMP SINGLE
!#ifdef bounds_exch_save


  CALL bounds_exch(1,'v+',TRP,'ocvad 7')
  CALL bounds_exch(1,'v+',TRM,'ocvad 8')
  CALL bounds_exch(1,'v+',WTP,'ocvad 9')
  CALL bounds_exch(1,'v+',WTM,'ocvad 10')

!#endif

    if( p_joff.eq.0 .and. lbounds_exch_tp ) then
       do k=1,ke
          do i=2,ie1
             tf=trp(i,1,k)
             wf=wtp(i,1,k)
             !    wtp(i,1,k)=wtm(i,1,k)
             !    trp(i,1,k)=-trm(i,1,k)
             wtm(i,1,k)=wf
             trm(i,1,k)=-tf
          enddo
       enddo
    endif


!$OMP END SINGLE

!$OMP DO
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        IF (amsue(i, j, k) .GT. 0.5_wp) THEN
          !
          RN=TRVOL(I,J,K)+WTP(I,J+1,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J-1,K)
          TRF(I,J,K)=(TRF(I,J,K)*TRVOL(I,J,K)+TRP(I,J+1,K)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I,J-1,K))/RN
          TRVOL(I,J,K)=RN
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !
!$OMP END PARALLEL
!#ifdef bounds_exch_save
  CALL bounds_exch(1,'v+',TRVOL,'ocvad 11')
  CALL bounds_exch(1,'v',TRF,'ocvad 12')
!#endif

! TRANSPORTS IN Z-DIRECTION (for fractional advection with 1/3 of time step)
  IF (iocaduv .EQ. 8)  CALL ocvad_z(trf, trvol, wtm, wtp, trm, trp)

END SUBROUTINE OCVAD

SUBROUTINE ocvad_z(trf, trvol, wtm, wtp, trm, trp)

  USE mo_kind, ONLY : dp, i4, wp
  USE mo_param1, ONLY : ie, je, ke, ie1, je1, kep
  USE mo_boundsexch, ONLY : bounds_exch
  USE mo_commo1, ONLY : dt, iocaduv, dlxv, dlyv, amsue, wo

  IMPLICIT NONE

  REAL(dp), INTENT(INOUT) :: TRF(IE,JE,KE), TRVOL(IE,JE,KE)
  REAL(dp), INTENT(INOUT) :: TRP(IE,JE,KEP),TRM(IE,JE,0:KE)
  REAL(dp), INTENT(INOUT) :: WTP(IE,JE,KEP),WTM(IE,JE,0:KE)

  REAL(dp) :: dtfrac ! time step fraction for vertical transport
  REAL(dp) :: abl,zwabl,wup,wlo,up,um,scha,sch,schi,rn
  INTEGER(i4) :: i, j, k, klo, kup

!$OMP PARALLEL PRIVATE(i,j,k,klo,kup,abl,zwabl,wup,wlo,up,um,scha,sch,schi,rn)

  IF (iocaduv .EQ. 8) THEN
    dtfrac=1._dp/3._dp   ! vertical transport with one third of time step
  ELSE
    dtfrac=1._dp
  ENDIF

!$OMP DO
  DO J=1,JE
    DO I=1,IE
      TRP(I,J,KEP)=0._dp
      TRM(I,J,0)=0._dp
      WTP(I,J,KEP)=0._dp
      WTM(I,J,0)=0._dp
    ENDDO
  ENDDO
  !
!$OMP DO
  DO K=1,KE
    KLO=MIN(K+1,KE)
    KUP=MAX(K-1,1)
    DO J=2,JE1
      DO I=2,IE1
        ABL=ABS(TRF(I,J,KUP)-TRF(I,J,KLO))
        zwabl = ABS(trf(i, j, kup) + trf(i, j, klo) - 2._wp * trf(i, j, k))
        wup = 0.5_wp * (wo(i, j, k) + wo(i, j+1, k))
        wlo = 0.5_wp * (wo(i, j, k+1) + wo(i, j+1, k+1))
        up = 0.5_wp * dt * dlyv(i, j) * dlxv(i, j) * (wup + ABS(wup)) * dtfrac
        um = 0.5_wp * dt * dlyv(i, j) * dlxv(i, j) * (ABS(wlo) - wlo) * dtfrac
        !
        scha = MAX(0._wp, (abl - zwabl) / (abl + 1.E-20_wp))
        !
        sch = MIN(1._wp, scha * trvol(i, j, k) / (up + um + 1.E-20_wp))
        schi = 1._wp - sch
        trp(i, j, k) = up * (schi * trf(i, j, k) &
             &               + sch * 0.5_wp * (trf(i, j, k) + trf(i, j, kup)))
        trm(i, j, k) = um * (schi * trf(i, j, k) &
             &               + sch * 0.5_wp * (trf(i, j, k) + trf(i, j, klo)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      ENDDO
    ENDDO
  ENDDO

  WTP(:,:,1)=0._dp

!$OMP DO
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        IF (amsue(i, j, k) .GT. 0.5_wp) THEN
          !RJ        CONT=TRVOL(I,J,K)*TRF(I,J,K)
          !
          RN=TRVOL(I,J,K)+WTP(I,J,K+1)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J,K-1)
          TRF(I,J,K)=(TRF(I,J,K)*TRVOL(I,J,K)+TRP(I,J,K+1)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I,J,K-1))/RN
          TRVOL(I,J,K)=RN
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !
!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch(1,'v+',TRVOL,'ocvad_z 1')
  CALL bounds_exch(1,'v',TRF,'ocvad_z 2')
!$OMP END SINGLE
!#endif

END SUBROUTINE ocvad_z
