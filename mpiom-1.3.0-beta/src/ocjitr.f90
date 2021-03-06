SUBROUTINE OCJITR
#ifdef GMBOLUS
  USE MO_PARAM1
  USE MO_PARALLEL
  USE MO_COMMO1
  USE MO_COMMOAU1
  USE MO_UNITS
#ifdef PBGC
  USE mo_carbch
  USE mo_sedmnt
  USE mo_biomod
  USE mo_control_bgc
  use mo_param1_bgc 
#endif /*PBGC*/
  !
  !     PARAMETRIZATION OF SUBGRID EDDY EFFECTS ACCORDING GENT AND 
  !     JIM MCWILLIAMS, 
  !      GENT ET AL. 1995, JPO 25,463.
  !UWE
  !     OUTLINE BY E. MAIER-REIMER
  !     IMPLEMENTED BY U. MIKOLAJEWICZ 4/99
  !     MODIFIED 6/2000 INCLUDE SEA LEVEL
  !
  !     BOLX,BOLY  GM COEFFICIENT K
  !     BOUNDARY CONDITION IMPLEMENTED IN DETERMINATION OF RHO-GRADIENTS
  !     FINAL ADVECTION WITH AN UPWIND SCHEME
  !
  implicit none

  real ZSURGM
  integer i,j,k,l
  DIMENSION ZSURGM(KE)
  real tsup_p(ie,je), ssup_p(ie,je), tlow_p(ie,je), slow_p(ie,je)
  real roxo,roxu,royo,royu,rozxo,rozxu,rozyo,rozyu
  real tm,sm,wun,wob,uwe,uos,vsu,vno,dhi
  real bolk,valinf,stabmin
  !
  !UWE VALUE FOR INFINITY
  VALINF=1.E30
  STABMIN=1.E-3
  DO K=1,KE
    ZSURGM(K)=0.
  ENDDO
  ZSURGM(1)=1.
  !
  !
#ifdef GMVISETAL
  !JJ USE VISBECK ET AL. (JPO, 1997) FORMULATION TO CALCULATE
  !   EDDY TRANSFER COEFFICIENTS
  !   BOLX, BOLY ARE CALCULATED IN SBR OCTHER
  !JJ
#else
  !UWE  MAKE BOLUS COEFFICIENT A LINEAR FUNCTION OF DX
  !     BOLX IN M**2/S, BOLK IN M/S
  !HH   SET GM DIFFUSION TO HARMOIC DIFFCOEF FROM NAMELIST
  BOLK=AH00
#ifdef BOLK025
  BOLK=BOLK*0.25
#endif /*BOLK025*/
#ifdef BOLK05
  BOLK=BOLK*0.5
#endif /*BOLK05*/
  DO J=1,JE
    DO I=1,IE
      BOLX(I,J)=MIN(BOLK*DLXU(I,J),2000.)
      BOLY(I,J)=MIN(BOLK*DLYV(I,J),2000.)
    ENDDO
  ENDDO
  !
#endif /*GMVISETAL*/

!$OMP PARALLEL PRIVATE(i,j,k,l,roxo,roxu,royo,royu,rozxo,rozxu,rozyo,rozyu, &
!$OMP    tsup_p,ssup_p,tlow_p,slow_p,tm,sm,wun,wob,uwe,uos,vsu,vno,dhi)

!$OMP DO
  DO K=1,KE
    UPTRT(K)=0.
    DOTRT(K)=0.
    DO J=1,JE
      DO I=1,IE
        UK1O(I,J,K)=0.
        VK1E(I,J,K)=0.
        WGO(I,J,K)=0.
      ENDDO
    ENDDO
  ENDDO
  !
!$OMP DO
  DO K=2,KE-1
    DO J=2,JE1
      DO I=2,IE1
        ROXO=0.5*(RHOO(I+1,J,K)-RHOO(I,J,K)+RHOO(I+1,J,K-1)               &
             -RHOO(I,J,K-1))/DLXP(I,J)
        ROXU=0.5*(RHOO(I+1,J,K)-RHOO(I,J,K)+(RHOO(I+1,J,K+1)              &
             -RHOO(I,J,K+1))*AMSUO(I,J,K+1))/DLXP(I,J)
        ROYO=0.5*(RHOO(I,J,K)-RHOO(I,J+1,K)+RHOO(I,J,K-1)                 &
             -RHOO(I,J+1,K-1))/DLYP(I,J)
        ROYU=0.5*(RHOO(I,J,K)-RHOO(I,J+1,K)+(RHOO(I,J,K+1)                &
             -RHOO(I,J+1,K+1))*AMSUE(I,J,K+1))/DLYP(I,J)
        !
        ROZXO=0.5*(STABIO(I+1,J,K)+STABIO(I,J,K))
        ROZXU=0.5*(STABIO(I+1,J,K+1)+STABIO(I,J,K+1))
        ROZYO=0.5*(STABIO(I,J+1,K)+STABIO(I,J,K))
        ROZYU=0.5*(STABIO(I,J+1,K+1)+STABIO(I,J,K+1))
        !
        ROZXO=MAX(ROZXO,STABMIN)
        ROZXU=MAX(ROZXU,STABMIN)
        ROZYO=MAX(ROZYO,STABMIN)
        ROZYU=MAX(ROZYU,STABMIN)
        IF(AMSUO(I,J,K+1).LT.0.5)ROZXU=VALINF
        IF(AMSUE(I,J,K+1).LT.0.5)ROZYU=VALINF
        !C      UK1O(I,J,K)=DWI(K)*(ROXO/ROZXO-ROXU/ROZXU)*AMSUO(I,J,K)
        !C      VK1E(I,J,K)=DWI(K)*(ROYO/ROZYO-ROYU/ROZYU)*AMSUE(I,J,K)
        !UWE  USE CORRECT LAYER THICKNESS
        !
        UK1O(I,J,K)=-BOLX(I,J)*(ROXO/ROZXO-ROXU/ROZXU)*AMSUO(I,J,K)       &
             /(DDUO(I,J,K)+(1.-AMSUO(I,J,K)))
        VK1E(I,J,K)=-BOLY(I,J)*(ROYO/ROZYO-ROYU/ROZYU)*AMSUE(I,J,K)       &
             /(DDUE(I,J,K)+(1.-AMSUE(I,J,K)))
        !
      ENDDO
    ENDDO
  ENDDO
  !
  !UWE INCLUDE SURFACE AND BOTTOM LAYERS
  !
  !     SURFACE LAYER
  !
  K=1
!$OMP DO
  DO J=2,JE1
    DO I=2,IE1
      ROXO=0.
      ROXU=0.5*(RHOO(I+1,J,K)-RHOO(I,J,K)+(RHOO(I+1,J,K+1)              &
           -RHOO(I,J,K+1))*AMSUO(I,J,K+1))/DLXP(I,J)
      ROYO=0.
      ROYU=0.5*(RHOO(I,J,K)-RHOO(I,J+1,K)+(RHOO(I,J,K+1)                &
           -RHOO(I,J+1,K+1))*AMSUE(I,J,K+1))/DLYP(I,J)
      !
      ROZXO=VALINF
      ROZXU=0.5*(STABIO(I+1,J,K+1)+STABIO(I,J,K+1))
      ROZYO=VALINF
      ROZYU=0.5*(STABIO(I,J+1,K+1)+STABIO(I,J,K+1))
      !
      ROZXU=MAX(ROZXU,STABMIN)
      ROZYU=MAX(ROZYU,STABMIN)
      !
      UK1O(I,J,K)=-BOLX(I,J)*DWI(K)*(ROXO/ROZXO-ROXU/ROZXU)*AMSUO(I,J,K)
      VK1E(I,J,K)=-BOLY(I,J)*DWI(K)*(ROYO/ROZYO-ROYU/ROZYU)*AMSUE(I,J,K)
    ENDDO
  ENDDO
  !
  !     BOTTOM LAYER
  !
  K=KE
!$OMP DO
  DO J=2,JE1
    DO I=2,IE1
      ROXO=0.5*(RHOO(I+1,J,K)-RHOO(I,J,K)+RHOO(I+1,J,K-1)               &
           -RHOO(I,J,K-1))/DLXP(I,J)
      ROXU=0.
      ROYO=0.5*(RHOO(I,J,K)-RHOO(I,J+1,K)+RHOO(I,J,K-1)                 &
           -RHOO(I,J+1,K-1))/DLYP(I,J)
      ROYU=0.
      !
      ROZXO=0.5*(STABIO(I+1,J,K)+STABIO(I,J,K))
      ROZXU=VALINF
      ROZYO=0.5*(STABIO(I,J+1,K)+STABIO(I,J,K))
      ROZYU=VALINF
      !
      ROZXO=MAX(ROZXO,STABMIN)
      ROZYO=MAX(ROZYO,STABMIN)
      !
      UK1O(I,J,K)=-BOLX(I,J)*(ROXO/ROZXO-ROXU/ROZXU)*AMSUO(I,J,K)       &
           /(DDUO(I,J,K)+(1.-AMSUO(I,J,K)))
      VK1E(I,J,K)=-BOLY(I,J)*(ROYO/ROZYO-ROYU/ROZYU)*AMSUE(I,J,K)       &
           /(DDUE(I,J,K)+(1.-AMSUE(I,J,K)))
    ENDDO
  ENDDO
!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('u',UK1O,'ocjitr 1')
  CALL bounds_exch('v',VK1E,'ocjitr 2')
!$OMP END SINGLE
!#endif
  !OtB_TEST_GMBOLUS ifdef TEST_GMBOLUS
  !OtB_TEST_GMBOLUS 
  !OtB_TEST_GMBOLUS UWE TEST
  !OtB_TEST_GMBOLUS 
  !OtB_TEST_GMBOLUS      UMAX=0.
  !OtB_TEST_GMBOLUS      VMAX=0.
  !OtB_TEST_GMBOLUS      DO K=1,KE
  !OtB_TEST_GMBOLUS      DO J=2,JE-1
  !OtB_TEST_GMBOLUS       DO I=2,IE-1
  !OtB_TEST_GMBOLUS        USUM(I,J)=USUM(I,J)+UK1O(I,J,K)*DDUO(I,J,K)
  !OtB_TEST_GMBOLUS        VSUM(I,J)=VSUM(I,J)+VK1E(I,J,K)*DDUE(I,J,K)
  !OtB_TEST_GMBOLUS        IF(ABS(UK1O(I,J,K)).GT.UMAX)THEN
  !OtB_TEST_GMBOLUS         UMAX=ABS(UK1O(I,J,K))
  !OtB_TEST_GMBOLUS         IUM=I
  !OtB_TEST_GMBOLUS         JUM=J
  !OtB_TEST_GMBOLUS         KUM=K
  !OtB_TEST_GMBOLUS        ENDIF
  !OtB_TEST_GMBOLUS        IF(ABS(VK1E(I,J,K)).GT.VMAX)THEN
  !OtB_TEST_GMBOLUS         VMAX=ABS(VK1E(I,J,K))
  !OtB_TEST_GMBOLUS         IVM=I
  !OtB_TEST_GMBOLUS         JVM=J
  !OtB_TEST_GMBOLUS         KVM=K
  !OtB_TEST_GMBOLUS        ENDIF
  !OtB_TEST_GMBOLUS        ENDDO
  !OtB_TEST_GMBOLUS       ENDDO
  !OtB_TEST_GMBOLUS      ENDDO
  !OtB_TEST_GMBOLUS      DO K=1,KE
  !OtB_TEST_GMBOLUS      I=IUM
  !OtB_TEST_GMBOLUS      J=JUM
  !OtB_TEST_GMBOLUS      WRITE(IO_STDOUT,*)'UMAX: ',IUM,JUM,K,UK1O(IUM,JUM,K)              &
  !OtB_TEST_GMBOLUS     &    ,UK1O(IUM,JUM,K)*DDUO(IUM,JUM,K)
  !OtB_TEST_GMBOLUS HH  K=>0 and KE+1 !!!
  !OtB_TEST_GMBOLUS      PRINT*,'TEST_GMBOLUS K=0 and KE+1 !!! CHECK CODE'
  !OtB_TEST_GMBOLUS      stop
  !OtB_TEST_GMBOLUS      ROXO=0.5*(RHOO(I+1,J,K)-RHOO(I,J,K)+RHOO(I+1,J,K-1)               &
  !OtB_TEST_GMBOLUS     & -RHOO(I,J,K-1))/DLXP(I,J)
  !OtB_TEST_GMBOLUS      ROXU=0.5*(RHOO(I+1,J,K)-RHOO(I,J,K)+(RHOO(I+1,J,K+1)              &
  !OtB_TEST_GMBOLUS     & -RHOO(I,J,K+1))*AMSUO(I,J,K+1))/DLXP(I,J)
  !OtB_TEST_GMBOLUS      ROYO=0.5*(RHOO(I,J,K)-RHOO(I,J+1,K)+RHOO(I,J,K-1)                 &
  !OtB_TEST_GMBOLUS     & -RHOO(I,J+1,K-1))/DLYP(I,J)
  !OtB_TEST_GMBOLUS      ROYU=0.5*(RHOO(I,J,K)-RHOO(I,J+1,K)+(RHOO(I,J,K+1)                &
  !OtB_TEST_GMBOLUS     & -RHOO(I,J+1,K+1))*AMSUE(I,J,K+1))/DLYP(I,J)
  !OtB_TEST_GMBOLUS 
  !OtB_TEST_GMBOLUS      ROZXO=0.5*(STABIO(I+1,J,K)+STABIO(I,J,K))
  !OtB_TEST_GMBOLUS      ROZXU=0.5*(STABIO(I+1,J,K+1)+STABIO(I,J,K+1))
  !OtB_TEST_GMBOLUS      ROZYO=0.5*(STABIO(I,J+1,K)+STABIO(I,J,K))
  !OtB_TEST_GMBOLUS      ROZYU=0.5*(STABIO(I,J+1,K+1)+STABIO(I,J,K+1))
  !OtB_TEST_GMBOLUS      WRITE(IO_STDOUT,*)RHOO(I+1,J,K),RHOO(I,J,K)                       &
  !OtB_TEST_GMBOLUS     &  ,STABIO(I+1,J,K),STABIO(I,J,K)
  !OtB_TEST_GMBOLUS C      WRITE(IO_STDOUT,*)ROXO,ROXU,ROZXO,ROZXU
  !OtB_TEST_GMBOLUS      ENDDO
  !OtB_TEST_GMBOLUS      DO K=1,KE
  !OtB_TEST_GMBOLUS      WRITE(IO_STDOUT,*)'VMAX: ',IVM,JVM,K,VK1E(IVM,JVM,K)              &
  !OtB_TEST_GMBOLUS     &                  ,DDUE(IVM,JVM,K)
  !OtB_TEST_GMBOLUS      ENDDO
  !OtB_TEST_GMBOLUS      TUMA=0.
  !OtB_TEST_GMBOLUS      DO J=2,JE-1
  !OtB_TEST_GMBOLUS       DO I=2,IE-1
  !OtB_TEST_GMBOLUS        IF(ABS(USUM(I,J)).GT.TUMA)THEN
  !OtB_TEST_GMBOLUS         II=I
  !OtB_TEST_GMBOLUS         JJ=J
  !OtB_TEST_GMBOLUS         TUMA=ABS(USUM(I,J))
  !OtB_TEST_GMBOLUS        ENDIF
  !OtB_TEST_GMBOLUS        IF(ABS(VSUM(I,J)).GT.TUMA)THEN
  !OtB_TEST_GMBOLUS         II=I
  !OtB_TEST_GMBOLUS         JJ=J
  !OtB_TEST_GMBOLUS         TUMA=ABS(VSUM(I,J))
  !OtB_TEST_GMBOLUS        ENDIF
  !OtB_TEST_GMBOLUS       ENDDO
  !OtB_TEST_GMBOLUS      ENDDO
  !OtB_TEST_GMBOLUS      WRITE(IO_STDOUT,*)'VERTINT: ',II,JJ,USUM(II,JJ)                   &
  !OtB_TEST_GMBOLUS     &                  ,VSUM(II,JJ),DEUTE(II,JJ)
  !OtB_TEST_GMBOLUS      SSU=0.
  !OtB_TEST_GMBOLUS      SSV=0.
  !OtB_TEST_GMBOLUS      DO K=1,KE
  !OtB_TEST_GMBOLUS       SSU=SSU+UK1O(II,JJ,K)*DDUO(II,JJ,K)
  !OtB_TEST_GMBOLUS       SSV=SSV+VK1E(II,JJ,K)*DDUE(II,JJ,K)
  !OtB_TEST_GMBOLUS 
  !OtB_TEST_GMBOLUS       WRITE(IO_STDOUT,*)UK1O(II,JJ,K)                                  &
  !OtB_TEST_GMBOLUS     &   ,DDUO(II,JJ,K),UK1O(II,JJ,K)*DDUO(II,JJ,K)                     &
  !OtB_TEST_GMBOLUS     &,SSU
  !OtB_TEST_GMBOLUS      ENDDO
  !OtB_TEST_GMBOLUS      WRITE(IO_STDOUT,*)'NACH 102'
  !OtB_TEST_GMBOLUS endif /*TEST_GMBOLUS*/
  !
  !     VERTICAL VELOCITY = VERTICAL INTEGRAL OF DIVERGENCE OF
  !                             HORIZONTAL VELOCITY FIELD
  !                         IN OCJITR
  !
!$OMP DO
  DO J=1,JE
    DO I=1,IE
      WGO(I,J,KEP) = ZERO
    ENDDO
  ENDDO
  !
!$OMP DO
  DO J=2,JE1
    DO K=KE,1,-1
      DO I=2,IE1
        WGO(I,J,K) = WGO(I,J,K+1)                                             &
             + DTI*WETO(I,J,K)*(                                              &
             DTDXPO(I,J)   * (   UK1O(I-1,J,K) * DDUO(I-1,J,K)*DLYU(I-1,J)    &
             - UK1O(I,J,K)   * DDUO(I,J,K)*DLYU(I,J)   )/DLYP(I,J)  )         &
             + DTI*WETO(I,J,K)* (                                             &
             + DTDYO(I,J) *                                                   &
             (   VK1E(I,J,K)     * DDUE(I,J,K)*DLXV(I,J)                      &
             - VK1E(I,J-1,K)   * DDUE(I,J-1,K)*DLXV(I,J-1)  )/DLXP(I,J)  )
      ENDDO
    ENDDO
  ENDDO

!$OMP DO
  DO K=1,KE
    UPTRA(K)=0.
    DOTRA(K)=0.
    UPTRT(K)=0.
    DOTRT(K)=0.
    !
    DO J=2,JE1
      DO I=2,IE1
        !
        UPTRT(K)=UPTRT(K)+DLXP(I,J)*DLYP(I,J)*(WGO(I,J,K)+ABS(WGO(I,J,K)))
        DOTRT(K)=DOTRT(K)-DLXP(I,J)*DLYP(I,J)*(WGO(I,J,K)-ABS(WGO(I,J,K)))
        !
      ENDDO
    ENDDO
    !
    UPTRT(K)=UPTRT(K)*1.E-6
    DOTRT(K)=DOTRT(K)*1.E-6
    !
  ENDDO
!$OMP SINGLE
  !
!#ifdef bounds_exch_save
  CALL bounds_exch('p',WGO,'ocjitr 3')
!#endif

  CALL global_sum(UPTRT)
  CALL global_sum(DOTRT)
  !
!HH
!6626 FORMAT(10F8.2)
!  WRITE(IO_STDOUT,6626)UPTRT
!$OMP END SINGLE
  !     
!$OMP DO
  DO K=1,KE
    IF(K.EQ.1) THEN
      DO J=1,JE
        DO I=1,IE
          TSUP_P(I,J)=THO(I,J,1)
          SSUP_P(I,J)=SAO(I,J,1)
        ENDDO
      ENDDO
    ELSE
      DO J=1,JE
        DO I=1,IE
          TSUP_P(I,J)=THO(I,J,K-1)
          SSUP_P(I,J)=SAO(I,J,K-1)
        ENDDO
      ENDDO
    ENDIF
    IF(K.EQ.KE) THEN
      DO J=1,JE
        DO I=1,IE
          TLOW_P(I,J)=THO(I,J,KE)
          SLOW_P(I,J)=SAO(I,J,KE)
        ENDDO
      ENDDO
    ELSE
      DO J=1,JE
        DO I=1,IE
          TLOW_P(I,J)=THO(I,J,K+1)
          SLOW_P(I,J)=SAO(I,J,K+1)
        ENDDO
      ENDDO
    ENDIF
    DO J=2,JE1
      !
      !UWE  MAKE ADVECTION MASS CONSERVING
      !
      DO I=2,IE1
        TM=THO(I,J,K)
        SM=SAO(I,J,K)
        WUN=HALF*(WGO(I,J,K+1)+ABS(WGO(I,J,K+1)))*DLXP(I,J)*DLYP(I,J)
        WOB=HALF*(ABS(WGO(I,J,K))-WGO(I,J,K))*DLXP(I,J)*DLYP(I,J)
        UWE=HALF*(UK1O(I-1,J,K)+ABS(UK1O(I-1,J,K)))                         &
             *DLYU(I-1,J)*DDUO(I-1,J,K)
        UOS=HALF*(ABS(UK1O(I,J,K))-UK1O(I,J,K))                             &
             *DLYU(I,J)*DDUO(I,J,K)
        VSU=HALF*(VK1E(I,J,K)+ABS(VK1E(I,J,K)))                             &
             *DLXV(I,J)*DDUE(I,J,K)
        VNO=HALF*(ABS(VK1E(I,J-1,K))-VK1E(I,J-1,K))                         &
             *DLXV(I,J-1)*DDUE(I,J-1,K)
        T1O(I,J,K)=(TM*DLXP(I,J)*DLYP(I,J)*(DDPO(I,J,K)+ALMZER              &
             +ZSURGM(K)*(ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA))  &
             +DT*                                                           &
             (WOB*(TSUP_P(I,J)-TM)+WUN*(TLOW_P(I,J)-TM)                     &
             +UWE*(THO(I-1,J,K)-TM)+UOS*(THO(I+1,J,K)-TM)                   &
             +VNO*(THO(I,J-1,K)-TM)+VSU*(THO(I,J+1,K)-TM)))                 &
             /(DLXP(I,J)*DLYP(I,J)*(DDPO(I,J,K)+ALMZER                      &
             +ZSURGM(K)*(ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)))
        !
        S1O(I,J,K)=(SM*DLXP(I,J)*DLYP(I,J)*(DDPO(I,J,K)+ALMZER              &
             +ZSURGM(K)*(ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA))  &
             +DT*                                                           &
             (WOB*(SSUP_P(I,J)-SM)+WUN*(SLOW_P(I,J)-SM)                     &
             +UWE*(SAO(I-1,J,K)-SM)+UOS*(SAO(I+1,J,K)-SM)                   &
             +VNO*(SAO(I,J-1,K)-SM)+VSU*(SAO(I,J+1,K)-SM)))                 &
             /(DLXP(I,J)*DLYP(I,J)*(DDPO(I,J,K)+ALMZER                      &
             +ZSURGM(K)*(ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)))
        !
        !
      ENDDO
    ENDDO
  ENDDO
  !
!$OMP DO
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        THO(I,J,K)=T1O(I,J,K)
        SAO(I,J,K)=S1O(I,J,K)
      ENDDO
    ENDDO
  ENDDO
!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('p',THO,'ocjitr 4')
  CALL bounds_exch('p',SAO,'ocjitr 5')
!$OMP END SINGLE
!#endif

#ifdef PBGC
  !
  ! ...for ocean tracer ...
  !
  DO  l=1,nocetra

!$OMP DO
    DO  k=1,ke

      IF(K.EQ.1) THEN
        DO  J=1,JE
          DO  I=1,IE
            tsup_p(i,j)=ocetra(i,j,1,l)
          ENDDO
        ENDDO
      ELSE
        DO  J=1,JE
          DO  I=1,IE
            tsup_p(i,j)=ocetra(i,j,k-1,l)
          ENDDO
        ENDDO
      ENDIF

      IF(K.EQ.KE) THEN
        DO  J=1,JE
          DO  I=1,IE
            tlow_p(i,j)=ocetra(i,j,ke,l)
          ENDDO
        ENDDO
      ELSE
        DO  J=1,JE
          DO  I=1,IE
            tlow_p(i,j)=ocetra(i,j,k+1,l)
          ENDDO
        ENDDO
      ENDIF

      DO  J=2,JE1
        DHI=DTI*DH
        DO  I=2,IE1
          IF ( weto(i,j,k).GT. 0.5 ) THEN
            tm=ocetra(i,j,k,l)
            wun= half*dlxp(i,j)*dlyp(i,j)                            &
                 *(wgo(i,j,k+1)+abs(wgo(i,j,k+1)))                 

            wob= half*dlxp(i,j)*dlyp(i,j)                            &
                 *(abs(wgo(i,j,k))-wgo(i,j,k))
            uwe= half*dlyu(i-1,j)*dduo(i-1,j,k)                      &
                 *(uk1o(i-1,j,k)+abs(uk1o(i-1,j,k)))
            uos= half*dlyu(i,j)*dduo(i,j,k)                          &
                 *(abs(uk1o(i,j,k))-uk1o(i,j,k))
            vsu= half*dlxv(i,j)*ddue(i,j,k)                          &
                 *(vk1e(i,j,k)+abs(vk1e(i,j,k)))
            vno= half*dlxv(i,j-1)*ddue(i,j-1,k)                      &
                 *(abs(vk1e(i,j-1,k))-vk1e(i,j-1,k))
            t1o(i,j,k)=                                                   &
                 (tm*dlxp(i,j)*dlyp(i,j)*(                                &
                 ddpo(i,j,k)+almzer+zsurgm(k)*                            &
                 (zo(i,j)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa) )      &
                 +dt*(wob*(tsup_p(i,j)      -tm)+wun*(tlow_p(i,j)  -tm)   &
                 +  uwe*(ocetra(i-1,j,k,l)-tm)+uos*(ocetra(i+1,j,k,l)-tm) &
                 +vno*(ocetra(i,j-1,k,l)-tm)+vsu*(ocetra(i,j+1,k,l)-tm))) &
                 /(dlxp(i,j)*dlyp(i,j)*(                                  &
                 ddpo(i,j,k)+almzer+zsurgm(k)*                            &
                 (zo(i,j)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa) ) ) 
          ELSE
            t1o(i,j,k)=rmasko
          ENDIF
        ENDDO
      ENDDO
    ENDDO

!$OMP DO
    DO K=1,KE
      DO  J=2,JE1
        DO  I=2,IE1
          ocetra(i,j,k,l)=t1o(i,j,k)
        ENDDO
      ENDDO
    ENDDO

!$OMP SINGLE
    CALL bounds_exch('p',ocetra(:,:,:,l),'ocjitr 6')
!$OMP END SINGLE

  ENDDO
#endif /*PBGC*/

!$OMP END PARALLEL

#endif /*GMBOLUS*/

END SUBROUTINE OCJITR
