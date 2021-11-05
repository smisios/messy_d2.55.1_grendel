SUBROUTINE ocjitr_base
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
  !     Changes by S. J. Lorenz, J.-O. Biesmann, 07/2007
  !     tracer-independent matrices calculated separately


  USE mo_kind, ONLY: dp
  USE mo_param1
  USE mo_parallel
  USE mo_boundsexch, ONLY : bounds_exch
  USE mo_commo1
  USE mo_commoau1
  USE mo_units
#ifdef PBGC
  USE mo_carbch
  USE mo_sedmnt
  USE mo_biomod
  USE mo_control_bgc
  use mo_param1_bgc
#endif /*PBGC*/

  !  #slo#: Calculation of matrices to be used in ocjitr_trf for all tracers
  !         UK1O, VK1E, WGO in mo_commo1
  !
  !
  IMPLICIT NONE
  INTEGER :: i,j,k
  REAL(wp) :: roxo,roxu,royo,royu,rozxo,rozxu,rozyo,rozyu
  REAL(wp) :: bolk
  !
  !UWE VALUE FOR INFINITY
  REAL(wp), PARAMETER :: valinf = 1.E30_wp
  REAL(wp), PARAMETER :: stabmin = 1.E-3_wp
  REAL(wp), ALLOCATABLE :: uptrt(:), dotrt(:)
  ALLOCATE(uptrt(ke), dotrt(ke))
  !
  !
  IF (ibolk .LT. 0) THEN
     !jj use visbeck et al. (jpo, 1997) formulation to calculate
     !   eddy transfer coefficients
     !   bolx, boly are calculated in sbr octher
     !jj
  ELSE
     !uwe  make bolus coefficient a linear function of dx
     !     bolx in m**2/s, bolk in m/s
     !hh   get gm diffusion from  namelist parameter ibolk

     bolk=REAL(ibolk,dp)/4.0e5_dp             ! scaled for 400 km grid space

     if ( icontro .gt.0 ) then
        write(0,*) 'bolk=',bolk
     endif


!$OMP PARALLEL
!$OMP DO
     DO j=1,je
        DO i=1,ie
           bolx(i, j) = MIN(bolk * dlxu(i, j), 2000._wp)
           boly(i, j) = MIN(bolk * dlyv(i, j), 2000._wp)
        ENDDO
     ENDDO
!$OMP END PARALLEL
  ENDIF

!$OMP PARALLEL PRIVATE(i,j,k,roxo,roxu,royo,royu,rozxo,rozxu,rozyo,rozyu)

!$OMP DO
  DO K=1,KE
    uptrt(k) = 0._wp
    dotrt(k) = 0._wp
    DO J=1,JE
      DO I=1,IE
        uk1o(i, j, k) = 0._wp
        vk1e(i, j, k) = 0._wp
        wgo(i, j, k) = 0._wp
      ENDDO
    ENDDO
  ENDDO
  !
!$OMP DO
  DO K=2,KE-1
    DO J=2,JE1
      DO I=2,IE1
        roxo = 0.5_wp * (rhoo(i+1, j, k) - rhoo(i, j, k) + rhoo(i+1, j, k-1) &
             -RHOO(I,J,K-1))/DLXP(I,J)
        roxu = 0.5_wp * (rhoo(i+1, j, k) - rhoo(i, j, k) + (rhoo(i+1, j, k+1) &
             -RHOO(I,J,K+1))*AMSUO(I,J,K+1))/DLXP(I,J)
        royo = 0.5_wp * (rhoo(i, j, k) - rhoo(i, j+1, k) + rhoo(i, j, k-1) &
             -RHOO(I,J+1,K-1))/DLYP(I,J)
        royu = 0.5_wp * (rhoo(i, j, k) - rhoo(i, j+1, k) + (rhoo(i, j, k+1) &
             -RHOO(I,J+1,K+1))*AMSUE(I,J,K+1))/DLYP(I,J)
        !
        rozxo = 0.5_wp * (stabio(i+1, j, k) + stabio(i, j, k))
        rozxu = 0.5_wp * (stabio(i+1, j, k+1) + stabio(i, j, k+1))
        rozyo = 0.5_wp * (stabio(i, j+1, k) + stabio(i, j, k))
        rozyu = 0.5_wp * (stabio(i, j+1, k+1) + stabio(i, j, k+1))
        !
        ROZXO=MAX(ROZXO,STABMIN)
        ROZXU=MAX(ROZXU,STABMIN)
        ROZYO=MAX(ROZYO,STABMIN)
        ROZYU=MAX(ROZYU,STABMIN)
        IF (amsuo(i, j, k+1) .LT. 0.5_wp) rozxu = valinf
        IF (amsue(i, j, k+1) .LT. 0.5_wp) rozyu = valinf
        !C      UK1O(I,J,K)=DWI(K)*(ROXO/ROZXO-ROXU/ROZXU)*AMSUO(I,J,K)
        !C      VK1E(I,J,K)=DWI(K)*(ROYO/ROZYO-ROYU/ROZYU)*AMSUE(I,J,K)
        !UWE  USE CORRECT LAYER THICKNESS
        !
        UK1O(I,J,K)=-BOLX(I,J)*(ROXO/ROZXO-ROXU/ROZXU)*AMSUO(I,J,K)       &
             /(dduo(i, j, k) + (1._wp - amsuo(i, j, k)))
        VK1E(I,J,K)=-BOLY(I,J)*(ROYO/ROZYO-ROYU/ROZYU)*AMSUE(I,J,K)       &
             /(ddue(i, j, k) + (1._wp - amsue(i, j, k)))
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
      roxo = 0._wp
      roxu = 0.5_wp * (rhoo(i+1, j, k) - rhoo(i, j, k) + (rhoo(i+1, j, k+1) &
           - rhoo(i, j, k+1)) * amsuo(i, j, k+1))/dlxp(i, j)
      royo = 0._wp
      royu = 0.5_wp * (rhoo(i, j, k) - rhoo(i, j+1, k) + (rhoo(i, j, k+1)   &
           -RHOO(I,J+1,K+1))*AMSUE(I,J,K+1))/DLYP(I,J)
      !
      ROZXO=VALINF
      rozxu = 0.5_wp * (stabio(i+1, j, k+1) + stabio(i, j, k+1))
      ROZYO=VALINF
      rozyu = 0.5_wp * (stabio(i, j+1, k+1) + stabio(i, j, k+1))
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
      roxo = 0.5_wp * (rhoo(i+1, j, k) - rhoo(i, j, k) + rhoo(i+1, j, k-1) &
           -RHOO(I,J,K-1))/DLXP(I,J)
      roxu = 0._wp
      royo = 0.5_wp * (rhoo(i, j, k) - rhoo(i, j+1, k) + rhoo(i, j, k-1) &
           -RHOO(I,J+1,K-1))/DLYP(I,J)
      royu = 0._wp
      !
      rozxo = 0.5_wp * (stabio(i+1, j, k) + stabio(i, j, k))
      ROZXU=VALINF
      rozyo = 0.5_wp * (stabio(i, j+1, k) + stabio(i, j, k))
      ROZYU=VALINF
      !
      ROZXO=MAX(ROZXO,STABMIN)
      ROZYO=MAX(ROZYO,STABMIN)
      !
      UK1O(I,J,K)=-BOLX(I,J)*(ROXO/ROZXO-ROXU/ROZXU)*AMSUO(I,J,K)       &
           /(dduo(i, j, k) + (1._wp - amsuo(i, j, k)))
      VK1E(I,J,K)=-BOLY(I,J)*(ROYO/ROZYO-ROYU/ROZYU)*AMSUE(I,J,K)       &
           /(ddue(i, j, k) + (1._wp - amsue(i, j, k)))
    ENDDO
  ENDDO
!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch(1,'u',UK1O,'ocjitr_base 1')
  CALL bounds_exch(1,'v',VK1E,'ocjitr_base 2')
!$OMP END SINGLE
!#endif
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
    uptrt(k) = 0._wp
    dotrt(k) = 0._wp
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
    uptrt(k) = uptrt(k) * 1.e-6_wp
    dotrt(k) = dotrt(k) * 1.e-6_wp
    !
  ENDDO

!$OMP SINGLE
!#ifdef bounds_exch_save
  CALL bounds_exch(1,'p',WGO,'ocjitr_base 3')
!#endif
!$OMP END SINGLE

!$OMP END PARALLEL

  CALL global_sum(UPTRT)
  CALL global_sum(DOTRT)

END SUBROUTINE ocjitr_base

SUBROUTINE ocjitr_trf(trf)

  USE mo_param1
  USE mo_boundsexch, ONLY : bounds_exch
  USE mo_commo1
  USE mo_planetary_constants, ONLY: rhoicwa, rhosnwa
  USE mo_units

! #slo#: Calculation of GM diffusion for all tracers, including T and S

  IMPLICIT NONE
  INTEGER :: i,j,k
  REAL(wp) :: trf(ie,je,ke)
  REAL(wp) :: dlxypi(ie,je), dlxyph(ie,je)
  REAL(wp) :: t2o(ie,je,ke)
  REAL(wp) :: tm,wun,wob,uwe,uos,vsu,vno
! REAL(wp) :: rmiss=999.9

! #slo#: geometric calculations (dlxyp) should be moved to a general routine
! #slo#  compare octdiff and ocadpo

#ifndef TRACER_OMP
!$OMP PARALLEL PRIVATE(i,j,k,tm,wun,wob,uwe,uos,vsu,vno)
!$OMP DO
#endif
      DO  J=2,JE1
        DO  I=2,IE1
        dlxypi(i, j) = 1._wp / (dlxp(i, j) * dlyp(i, j))
        dlxyph(i,j) = dlxp(i,j) * dlyp(i,j) * half
        ENDDO
      ENDDO

! for k=1
      k = 1
#ifndef TRACER_OMP
!$OMP DO
#endif
      DO  J=2,JE1
        DO  I=2,IE1
!
! #slo#: optional treatment of dry points: set to rmiss
! #slo#  old ocjitr: THO and SAO - calculate dry and wet points
! #slo#              OCETRA      - set to rmask, only available with PBGC
!         IF ( weto(i,j,k).GT. 0.5 ) THEN
            tm=trf(i,j,k)
            wun= dlxyph(i,j)                                         &
                 *(wgo(i,j,k+1)+abs(wgo(i,j,k+1)))
            wob= dlxyph(i,j)                                         &
                 *(abs(wgo(i,j,k))-wgo(i,j,k))
            uwe= half*dlyu(i-1,j)*dduo(i-1,j,k)                      &
                 *(uk1o(i-1,j,k)+abs(uk1o(i-1,j,k)))
            uos= half*dlyu(i,j)*dduo(i,j,k)                          &
                 *(abs(uk1o(i,j,k))-uk1o(i,j,k))
            vsu= half*dlxv(i,j)*ddue(i,j,k)                          &
                 *(vk1e(i,j,k)+abs(vk1e(i,j,k)))
            vno= half*dlxv(i,j-1)*ddue(i,j-1,k)                      &
                 *(abs(vk1e(i,j-1,k))-vk1e(i,j-1,k))
            t2o(i,j,k)= tm + dt * (                                  &
                   wob*(trf(i,j,k  )-tm) + wun*(trf(i,j,k+1)-tm)     &
                 + uwe*(trf(i-1,j,k)-tm) + uos*(trf(i+1,j,k)-tm)     &
                 + vno*(trf(i,j-1,k)-tm) + vsu*(trf(i,j+1,k)-tm) )   &
                                * dlxypi(i,j)/(ddpo(i,j,k)+almzer+   &
                  (zo(i,j)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa))

!         ELSE
!!          t2o(i,j,k)=tm
!           t2o(i,j,k)=rmiss
!         ENDIF
        ENDDO
      ENDDO

#ifndef TRACER_OMP
!$OMP DO
#endif
    DO  k=2,ke-1
      DO  J=2,JE1
        DO  I=2,IE1
            tm=trf(i,j,k)
            wun= dlxyph(i,j)                                         &
                 *(wgo(i,j,k+1)+abs(wgo(i,j,k+1)))
            wob= dlxyph(i,j)                                         &
                 *(abs(wgo(i,j,k))-wgo(i,j,k))
            uwe= half*dlyu(i-1,j)*dduo(i-1,j,k)                      &
                 *(uk1o(i-1,j,k)+abs(uk1o(i-1,j,k)))
            uos= half*dlyu(i,j)*dduo(i,j,k)                          &
                 *(abs(uk1o(i,j,k))-uk1o(i,j,k))
            vsu= half*dlxv(i,j)*ddue(i,j,k)                          &
                 *(vk1e(i,j,k)+abs(vk1e(i,j,k)))
            vno= half*dlxv(i,j-1)*ddue(i,j-1,k)                      &
                 *(abs(vk1e(i,j-1,k))-vk1e(i,j-1,k))
            t2o(i,j,k)= tm + dt * (                                  &
                   wob*(trf(i,j,k-1)-tm) + wun*(trf(i,j,k+1)-tm)     &
                 + uwe*(trf(i-1,j,k)-tm) + uos*(trf(i+1,j,k)-tm)     &
                 + vno*(trf(i,j-1,k)-tm) + vsu*(trf(i,j+1,k)-tm) )   &
                                * dlxypi(i,j)/(ddpo(i,j,k)+almzer)
        ENDDO
      ENDDO
    ENDDO

! for k=ke
      k=ke
#ifndef TRACER_OMP
!$OMP DO
#endif
      DO  J=2,JE1
        DO  I=2,IE1
            tm=trf(i,j,k)
            wun= dlxyph(i, j)                                        &
                 *(wgo(i,j,k+1)+abs(wgo(i,j,k+1)))
            wob= dlxyph(i, j)                                        &
                 *(abs(wgo(i,j,k))-wgo(i,j,k))
            uwe= half*dlyu(i-1,j)*dduo(i-1,j,k)                      &
                 *(uk1o(i-1,j,k)+abs(uk1o(i-1,j,k)))
            uos= half*dlyu(i,j)*dduo(i,j,k)                          &
                 *(abs(uk1o(i,j,k))-uk1o(i,j,k))
            vsu= half*dlxv(i,j)*ddue(i,j,k)                          &
                 *(vk1e(i,j,k)+abs(vk1e(i,j,k)))
            vno= half*dlxv(i,j-1)*ddue(i,j-1,k)                      &
                 *(abs(vk1e(i,j-1,k))-vk1e(i,j-1,k))
            t2o(i,j,k)= tm + dt * (                                  &
                   wob*(trf(i,j,k-1)-tm) + wun*(trf(i,j,k  )-tm)     &
                 + uwe*(trf(i-1,j,k)-tm) + uos*(trf(i+1,j,k)-tm)     &
                 + vno*(trf(i,j-1,k)-tm) + vsu*(trf(i,j+1,k)-tm) )   &
                                * dlxypi(i,j)/(ddpo(i,j,k)+almzer)
        ENDDO
      ENDDO

#ifndef TRACER_OMP
!$OMP WORKSHARE
#endif
     trf(2:ie1,2:je1,:) = t2o(2:ie1,2:je1,:)
#ifndef TRACER_OMP
!$OMP END WORKSHARE
#endif

!$OMP SINGLE
    CALL bounds_exch(1,'p',trf(:,:,:),'ocjitr_trf 1')
!$OMP END SINGLE

#ifndef TRACER_OMP
!$OMP END PARALLEL
#endif


END SUBROUTINE ocjitr_trf
! END MODULE mo_ocjitr
