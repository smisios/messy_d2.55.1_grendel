SUBROUTINE OCVTOT

  USE MO_PARAM1
  USE MO_MPI
  USE MO_PARALLEL
  USE MO_COMMO1
  USE MO_UNITS
  USE mo_boundsexch, ONLY : bounds_exch

  IMPLICIT NONE

  REAL(wp), ALLOCATABLE :: surdif(:,:), courau(:,:,:), courav(:,:,:), couraw(:,:,:)
  REAL(wp)              :: courau_max, courav_max, couraw_max

  INTEGER :: i, j, k, jb

  REAL(wp) :: vake, uako

  REAL(wp), ALLOCATABLE :: uptrt(:), dotrt(:)

!$OMP PARALLEL PRIVATE(i,j,k,vake,uako)

  IF ( .NOT. lwith_barotropic_stokes_drift) THEN
!$OMP DO
      DO  K=1,KE
        DO  J=1,JE
          DO  I=1,IE
            T1O(I,J,K)=UOO(I,J,K)
            S1O(I,J,K)=VOE(I,J,K)

            VKE(I,J,K)=AMSUE(I,J,K)*(VSE(I,J)+VKE(I,J,K))
            UKO(I,J,K)=AMSUO(I,J,K)*(USO(I,J)+UKO(I,J,K))

            VAKE=VKE(I,J,K)
            UAKO=UKO(I,J,K)

            VKE(I,J,K)=CONN*VKE(I,J,K)+VOE(I,J,K)*CONO
            UKO(I,J,K)=CONN*UKO(I,J,K)+UOO(I,J,K)*CONO

            UOO(I,J,K)=UAKO
            VOE(I,J,K)=VAKE
          ENDDO
        ENDDO
      ENDDO

    ELSE
!$OMP DO
      DO  K=2,KE
        DO  J=1,JE
          DO  I=1,IE

            T1O(I,J,K)=UOO(I,J,K)
            S1O(I,J,K)=VOE(I,J,K)

            VKE(I,J,K)=AMSUE(I,J,K)*(VSE(I,J)+VKE(I,J,K))
            UKO(I,J,K)=AMSUO(I,J,K)*(USO(I,J)+UKO(I,J,K))

            uk1o(i,j,k)=uko(i,j,k)
            vk1e(i,j,k)=vke(i,j,k)

            VKE(I,J,K)=(CONN*VKE(I,J,K)+VOE(I,J,K)*CONO)*ddue(i,j,k)
            UKO(I,J,K)=(CONN*UKO(I,J,K)+UOO(I,J,K)*CONO)*dduo(i,j,k)

          ENDDO
        ENDDO
      ENDDO

      K=1
!$OMP DO
      DO  J=1,JE
        DO  I=1,IE

          T1O(I,J,K)=UOO(I,J,K)
          S1O(I,J,K)=VOE(I,J,K)

          VKE(I,J,K)=AMSUE(I,J,K)*(VSE(I,J)+VKE(I,J,K))
          UKO(I,J,K)=AMSUO(I,J,K)*(USO(I,J)+UKO(I,J,K))

          uk1o(i,j,k)=uko(i,j,k)
          vk1e(i,j,k)=vke(i,j,k)

          VKE(I,J,K)=(surlen(i,j)*CONN*VKE(I,J,K)+VOE(I,J,K)*CONO*ddue(i,j,k))
          UKO(I,J,K)=(surlon(i,j)*CONN*UKO(I,J,K)+UOO(I,J,K)*CONO*dduo(i,j,k))

        ENDDO
      ENDDO

    ENDIF


!$OMP DO
      DO J=1,JE
        DO I=1,IE
          WO(I,J,KEP) = ZERO
        ENDDO
      ENDDO

!
!
!======================================================================
!
!     B)
!
!     VERTICAL VELOCITY = VERTICAL INTEGRAL OF DIVERGENCE OF
!                             HORIZONTAL VELOCITY FIELD
!
!
!RJ       UDIMA=0.
!RJ       K=1

      IF (icontro.NE.0) THEN
        ALLOCATE(uptrt(ke), dotrt(ke))
        uptrt(:) = 0._wp
        dotrt(:) = 0._wp
      END IF

      Jb=2
      IF ( have_g_js .AND. lbounds_exch_tp ) jb=3

      IF ( .NOT. lwith_barotropic_stokes_drift) THEN

!$OMP DO
        DO K=KE,1,-1
          DO J=jb,JE1
            DO I=2,IE1

              wo(i,j,k) = wo(i,j,k+1) + weto(i,j,k) * (                                 &
                   uko(i-1,j,k) * dduo(i-1,j,k)*dlyu(i-1,j)                            &
                   - uko(i,j,k)   * dduo(i,j,k)*dlyu(i,j)                                &
                   + vke(i,j,k)     * ddue(i,j,k)*dlxv(i,j)                              &
                   - vke(i,j-1,k)   * ddue(i,j-1,k)*dlxv(i,j-1)          )*areain(i,j)

              IF (icontro.NE.0) THEN

                uptrt(k)=uptrt(k)+area(i,j)*(wo(i,j,k)+ABS(wo(i,j,k)))
                dotrt(k)=dotrt(k)-area(i,j)*(wo(i,j,k)-ABS(wo(i,j,k)))

              ENDIF

            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL

      ELSE

!$OMP DO
        DO K=KE,1,-1
          DO J=jb,JE1
            DO I=2,IE1

              wo(i,j,k) = wo(i,j,k+1) + weto(i,j,k) * (                                 &
                   uko(i-1,j,k) *dlyu(i-1,j)                            &
                   - uko(i,j,k)   *dlyu(i,j)                                &
                   + vke(i,j,k)   *dlxv(i,j)                              &
                   - vke(i,j-1,k) *dlxv(i,j-1)   )*areain(i,j)

              IF (icontro.NE.0) THEN

                uptrt(k)=uptrt(k)+area(i,j)*(wo(i,j,k)+ABS(wo(i,j,k)))
                dotrt(k)=dotrt(k)-area(i,j)*(wo(i,j,k)-ABS(wo(i,j,k)))

              ENDIF

            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL

        uoo=uk1o
        voe=vk1e

        dduo(:,:,1)=surlon(:,:)
        ddue(:,:,1)=surlen(:,:)
        ! FIXME: why not use almzer here?
        uko = uko / (dduo + 1.e-20_wp)
        vke = vke / (ddue + 1.e-20_wp)

      ENDIF

      CALL bounds_exch(1,'p',WO,'ocvtot 1')

      ! below are diagnostics

      IF (icontro.NE.0) THEN
        ! FIXME: what's the meaning of 1.e-6 here?
        uptrt = uptrt * 1.e-6_wp
        dotrt = dotrt * 1.e-6_wp

        CALL global_sum(uptrt)
        CALL global_sum(dotrt)

        IF (p_pe == p_io) THEN
          WRITE(0,*)'TOTAL UPWELLING [sv]'
          DO k=1,ke
            WRITE(0,'(i2,2f12.2)')k,uptrt(k),dotrt(k)
          ENDDO
        ENDIF

        ALLOCATE (surdif(ie,je))
        surdif(:,:)=wo(:,:,1)*dt-z1o(:,:)
        WRITE(0,*) 'surdif check (values < 1.0e-12): ', p_pe, &
             MINVAL(surdif), MAXVAL(surdif), MAXLOC(surdif)
        DEALLOCATE(surdif)

        ALLOCATE (courau(ie,je,ke), courav(ie,je,ke), couraw(ie,je,ke))
        DO k=1,ke
          courau(:,:,k)=ABS(dt*uko(:,:,k)/dlxu(:,:))
          courav(:,:,k)=ABS(dt*vke(:,:,k)/dlyv(:,:))
          couraw(:,:,k)=ABS(dt*wo(:,:,k)/dz(k))
        ENDDO

        courau_max=MAXVAL(courau)
        courav_max=MAXVAL(courav)
        couraw_max=MAXVAL(couraw)

        IF (courau_max >= 0.95_wp .OR. courav_max >= 0.95_wp &
             .OR. couraw_max >= 0.95_wp) THEN
          WRITE(0,'(a,3i5,3(f7.3,3i5))') 'CFL check: ', p_pe, p_ioff, p_joff, &
               courau_max, MAXLOC(courau), &
               courav_max, MAXLOC(courav), &
               couraw_max, MAXLOC(couraw)
        ELSE
          WRITE(0,*) 'CFL check: ', p_pe, 'OK'
        ENDIF

        DEALLOCATE(courau, courav, couraw)

      ENDIF

    END SUBROUTINE OCVTOT
