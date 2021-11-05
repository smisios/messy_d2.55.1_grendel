SUBROUTINE OCMODMOM

  USE MO_PARAM1
  USE mo_planetary_constants, ONLY: g
  USE MO_PARALLEL
  USE mo_boundsexch, ONLY : bounds_exch
  USE MO_COMMO1
  USE MO_UNITS
  USE mo_grid, ONLY : shaswe, shaswo

#ifdef _PROFILE
  USE mo_profile,      ONLY: trace_start, trace_stop
#endif

  IMPLICIT NONE
  !
  !-=====================================================================
  !
  !    DECOMPOSITION INTO BAROTROPIC AND BAROCLINIC FIELD
  !
  REAL(wp) UCOR(IE,JE),VCOR(IE,JE)
  INTEGER i,j,k,iter
  REAL(wp) :: deuto0, deute0
  !
  !---------------------------------------------------------------------
  !
#ifdef _PROFILE
  CALL trace_start ('ocmodmom', 8)
#endif

  !$OMP PARALLEL PRIVATE(i,j,k,iter)

!$OMP DO
  u1o(:,:) = zero
  u1e(:,:) = zero
  v1o(:,:) = zero
  v1e(:,:) = zero
  uso(:,:) = zero
  vse(:,:) = zero

!$OMP SINGLE
!#ifdef bounds_exch_save
  CALL bounds_exch(1,'u',uko,'ocmodmom 1')
  CALL bounds_exch(1,'v',vke,'ocmodmom 2')
  CALL bounds_exch(1,'u',uoo,'ocmodmom 3')
  CALL bounds_exch(1,'v',voe,'ocmodmom 4')
  CALL bounds_exch(1,'u',uaccel,'ocmodmom 3a')
  CALL bounds_exch(1,'v',vaccel,'ocmodmom 4a')
!#endif
  CALL bounds_exch(1,'p',po,'ocmodmom 5')
  CALL bounds_exch(1,'p',zo,'ocmodmom 6')
!$OMP END SINGLE

  IF ( lwith_barotropic_stokes_drift ) THEN
    uzo(:,:)=zero
    vze(:,:)=zero

    deutio = 0._wp
    deutie = 0._wp

    DO j=2,je1
      DO i=2,ie1
        deuto0=MIN(DEPTO(I,J),DEPTO(I+1,J))
        deute0=MIN(DEPTO(I,J),DEPTO(I,J+1))
        deuto(i, j) = amsuo(i, j, 1) &
             * MAX(0._wp, deuto0 + 0.5_wp * (zo(i, j) + zo(i+1, j)))
        deute(i, j) = amsue(i, j, 1) &
             * MAX(0._wp, deute0 + 0.5_wp * (zo(i, j) + zo(i, j+1)))
        dduo(i, j, 1) = amsuo(i, j, 1) &
             * MAX(0._wp, dzw(1) + 0.5_wp * (zo(i, j) + zo(i+1, j)))
        ddue(i, j, 1) = amsue(i,j,1) &
             * MAX(0._wp, dzw(1) + 0.5_wp * (zo(i, j) + zo(i, j+1)))
      ENDDO
    ENDDO

!$OMP SINGLE
    CALL bounds_exch(1,'u+',deuto,'ocmodmom 5')
    CALL bounds_exch(1,'v+',deute,'ocmodmom 6')
!$OMP END SINGLE

    WHERE (deuto .GT. 0._wp) deutio = 1._wp/deuto
    WHERE (deute .GT. 0._wp) deutie = 1._wp/deute

  ENDIF

!$OMP SINGLE
  CALL bounds_exch(1,'u+',dduo,'ocmodmom 5')
  CALL bounds_exch(1,'v+',ddue,'ocmodmom 6')
!$OMP END SINGLE

!$OMP DO
  DO k=1,ke
    DO j=1,je
      DO i=1,ie
        uko(i,j,k)=uko(i,j,k)*dduo(i,j,k)
        vke(i,j,k)=vke(i,j,k)*ddue(i,j,k)
        uaccel(i,j,k)=uaccel(i,j,k)*dduo(i,j,k)
        vaccel(i,j,k)=vaccel(i,j,k)*ddue(i,j,k)
        uk1o(i,j,k)=uko(i,j,k)
        vk1e(i,j,k)=vke(i,j,k)
        stabio(i, j, k) = 0.001_wp * dz(k) * stabio(i, j, k) * weto(i, j, k)
        po(i,j,k)=po(i,j,k)+g*zo(i,j)
      ENDDO
    ENDDO
  ENDDO
  !
#ifdef _PROFILE
  CALL trace_start ('ocmodmom iteration-loop', 9)
#endif

  iteration_loop: DO ITER=1,12

   DO K=1,KE

    DO J=1,JE
        DO I=1,IE
          UCOR(I,J)=STABN*UK1O(I,J,K)+STABO*UKO(I,J,K)
          VCOR(I,J)=STABN*VK1E(I,J,K)+STABO*VKE(I,J,K)
        ENDDO
      ENDDO


    DO j=2,je1
        DO i=2,ie1

          !cuweneu   include earth curvature
          !c
          uk1o(i,j,k)=uko(i,j,k)+dtdxuo(i,j)*(po(i,j,k)-po(i+1,j,k))   &
               *dduo(i,j,k)  + uaccel(i,j,k)                           &
               + 0.25_wp * dt * (                                      &
               ftwou(i,j)                                              &
               + 0.25_wp * (curvav(i,j)+curvav(i+1,j)+curvav(i,j-1)    &
               +curvav(i+1,j))*ucor(i,j)                             &
               )*(vcor(i,j)+vcor(i+1,j)                            &
               +vcor(i,j-1)+vcor(i+1,j-1))*amsuo(i,j,k)
        ENDDO
      ENDDO

    IF (je1 > 1) THEN ! trick to prevent xlf compiler from fusion of loops
    DO j=2,je1
        DO i=2,ie1

          vk1e(i,j,k)=vke(i,j,k)+dpye(i,j)*(po(i,j+1,k)-po(i,j,k))     &
               *ddue(i,j,k)   + vaccel(i,j,k)                          &
               - 0.25_wp * dt * (ftwov(i,j)                            &
               +curvav(i,j) * 0.25_wp                                  &
               *(ucor(i,j)+ucor(i-1,j)+ucor(i,j+1)               &
               +ucor(i-1,j+1))                                       &
               )*(ucor(i,j)+ucor(i-1,j)                            &
               +ucor(i,j+1)+ucor(i-1,j+1))*amsue(i,j,k)
        ENDDO
      ENDDO
    ENDIF


    ENDDO  ! new k-loop

!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch(1,'u',UK1O,'ocmodmom 9')
    CALL bounds_exch(1,'v',VK1E,'ocmodmom 10')
!$OMP END SINGLE
!#endif

  ENDDO iteration_loop
#ifdef _PROFILE
  CALL trace_stop ('ocmodmom iteration-loop', 9)
#endif


  !     CALCULATION OF BAROTROPIC VELOCITIES ON FIELDS U1 AND V1

!$OMP DO
  DO J=2,JE1
    DO K=1,KE
      DO I=1,IE
        po(i,j,k) = po(i,j,k)-g*zo(i,j)
        v1e(i,j)  = v1e(i,j)+ddue(i,j,k)*voe(i,j,k)
        u1o(i,j)  = u1o(i,j)+dduo(i,j,k)*uoo(i,j,k)
        vse(i,j)  = vse(i,j)+vk1e(i,j,k)*amsue(i,j,k)
        uso(i,j)  = uso(i,j)+uk1o(i,j,k)*amsuo(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!$OMP DO
  DO k=1,ke
    DO j=2,je1
      DO i=1,ie
        uko(i,j,k)=uoo(i,j,k)-amsuo(i,j,k)*deutio(i,j)*u1o(i,j)
        vke(i,j,k)=voe(i,j,k)-amsue(i,j,k)*deutie(i,j)*v1e(i,j)
        uk1o(i,j,k)=uk1o(i,j,k)*(amsuo(i,j,k)/(almzer+dduo(i,j,k))) &
             -amsuo(i,j,k)*deutio(i,j)*uso(i,j)
        vk1e(i,j,k)=vk1e(i,j,k)*(amsue(i,j,k)/(almzer+ddue(i,j,k))) &
             -amsue(i,j,k)*deutie(i,j)*vse(i,j)

        ! lwith_one_layer_shelfs : grid points with only one layer do not have a baroclinic velocity
        uko(i,j,k)=uko(i,j,k)*shaswo(i,j)
        vke(i,j,k)=vke(i,j,k)*shaswe(i,j)
        uk1o(i,j,k)=uk1o(i,j,k)*shaswo(i,j)
        vk1e(i,j,k)=vk1e(i,j,k)*shaswe(i,j)

      ENDDO
    ENDDO
  ENDDO


!$OMP SINGLE
  CALL bounds_exch(1,'u',uk1o,'ocmodmom 9')
  CALL bounds_exch(1,'v',vk1e,'ocmodmom 11')
  CALL bounds_exch(1,'u',uko,'ocmodmom 9')
  CALL bounds_exch(1,'v',vke,'ocmodmom 12')
  CALL bounds_exch(1,'u',u1o,'ocmodmom 9')
  CALL bounds_exch(1,'v',v1e,'ocmodmom 10')
  CALL bounds_exch(1,'p',po,'ocmodmom 10')
  CALL bounds_exch(1,'u',uso,'ocmodmom 9')
  CALL bounds_exch(1,'v',vse,'ocmodmom 10')
!$OMP END SINGLE

  IF ( lwith_barotropic_stokes_drift ) THEN
    uzo=uso
    vze=vse

    uso=uzo*deutio
    vse=vze*deutie

  ENDIF

!$OMP END PARALLEL

#ifdef _PROFILE
  CALL trace_stop ('ocmodmom', 8)
#endif
!
END SUBROUTINE OCMODMOM
