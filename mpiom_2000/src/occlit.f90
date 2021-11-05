
SUBROUTINE occlit

  USE mo_param1
  USE mo_planetary_constants, ONLY: g
  USE mo_parallel
  USE mo_boundsexch, ONLY : bounds_exch
  USE mo_commo1
  USE mo_units
  USE mo_grid, ONLY : shamsuo, shamsue

#ifdef _PROFILE
  USE mo_profile,      ONLY: trace_start, trace_stop
#endif


  IMPLICIT NONE

  INTEGER :: i,j,k,iter,itermax, jb, id(2)
  REAL(wp) :: brems, contra, ove, speed, under, vorw, uko1, vke1
  !
  REAL(wp) :: bruva(ie,je),contrij(ie,je)
  REAL(wp), ALLOCATABLE :: contrij_g(:,:)
  !
  REAL(wp) :: wpo(ie,je,ke) ! changed from kep to ke by LK
  !
  !  ITERATION OF BAROCLINIC SYSTEM DIRECTLY IN BAROCLINIC VELOCITIES


#ifdef _PROFILE
  CALL trace_start ('occlit', 1)
#endif

  IF (p_pe==p_io) THEN
    ALLOCATE(contrij_g(ie_g,je_g))
  ELSE
    ALLOCATE(contrij_g(0,0))
  ENDIF

!$OMP PARALLEL PRIVATE(i,j,k,iter,itermax,under,uko1,vke1,ove,speed,brems,vorw)

#ifdef _PROFILE
  CALL trace_start ('occlit_1st_loop', 29)
#endif
!$OMP DO
  DO k=ke,1,-1

    IF (k.EQ.1) THEN

      bruva(:,:)    = 0._wp
      stabio(:,:,1) = 0._wp
      slow(:,:)     = 0._wp

      IF ( lwith_barotropic_stokes_drift ) THEN

        uko(:,:,1)=uko(:,:,1)*surlon(:,:)
        vke(:,:,1)=vke(:,:,1)*surlen(:,:)
        uk1o(:,:,1)=uk1o(:,:,1)*surlon(:,:)
        vk1e(:,:,1)=vk1e(:,:,1)*surlen(:,:)

      ELSE

        uko(:,:,1)=uko(:,:,1)*dduo(:,:,1)
        vke(:,:,1)=vke(:,:,1)*ddue(:,:,1)
        uk1o(:,:,1)=uk1o(:,:,1)*dduo(:,:,1)
        vk1e(:,:,1)=vk1e(:,:,1)*ddue(:,:,1)

      ENDIF

    ENDIF

    IF (k.GE.2) THEN
      uko(:,:,k)=uko(:,:,k)*dduo(:,:,k)
      vke(:,:,k)=vke(:,:,k)*ddue(:,:,k)
      uk1o(:,:,k)=uk1o(:,:,k)*dduo(:,:,k)
      vk1e(:,:,k)=vk1e(:,:,k)*ddue(:,:,k)
    ENDIF

    DO j=je1,2,-1
      DO i=ie1,2,-1
        wo(i,j,k) = wo(i,j,k+1) + dti * weto(i,j,k) * (                   &
             dtdxpo(i,j) * (   uko(i-1,j,k)*dlyu(i-1,j)                   &
             - uko(i,j,k)*dlyu(i,j)    )/dlyp(i,j)                        &
             + dtdyo(i,j)  * (   vke(i,j,k)*dlxv(i,j)                     &
             - vke(i,j-1,k)*dlxv(i,j-1))/dlxp(i,j)  )
      ENDDO
    ENDDO


#ifdef _PROFILE
  CALL trace_start ('occlit_boundsexchange1', 25)
#endif
!$OMP SINGLE
    CALL bounds_exch(1,'p',wo(:,:,k),'occlit 4')
!$OMP END SINGLE

#ifdef _PROFILE
  CALL trace_stop ('occlit_boundsexchange1', 25)
#endif

  ENDDO

#ifdef _PROFILE
  CALL trace_stop ('occlit_1st_loop', 29)
#endif

#ifdef _PROFILE
  CALL trace_start ('occlit_2nd_loop', 30)
#endif
!$OMP DO
  DO k=1,ke
     DO j=1,je
        DO i=1,ie
           bruva(i,j)=bruva(i,j)+tiestw(k+1)*stabio(i,j,k)
           slow(i,j)=slow(i,j)+stabio(i,j,k)*dt*wo(i,j,k)
           wpo(i,j,k)=slow(i,j)
        ENDDO
     ENDDO
  ENDDO

#ifdef _PROFILE
  CALL trace_stop ('occlit_2nd_loop', 30)
#endif
  !OtB Needs more iterations than 6
  !     itermax=6
  !     itermax=16
  !Uwe says 8 is cheaper
  itermax=8

  IF ( lbounds_exch_tp ) itermax=12
  !
  iteration_loop: DO iter=1,itermax

!$OMP DO
    ucos(:,:)  = 0._wp
    vcos(:,:)  = 0._wp
    pxoin(:,:) = 0._wp
    pyein(:,:) = 0._wp

    IF ( .NOT. lwith_barotropic_stokes_drift ) THEN

#ifdef _PROFILE
  CALL trace_start ('occlit_loop_pxoin', 24)
#endif
!$OMP DO
      DO j=2,je1

        DO k=1,ke
          DO i=2,ie1
            pxoin(i,j)=pxoin(i,j)+dduo(i,j,k)*dtdxuo(i,j)*                   &
                 (wpo(i+1,j,k)-wpo(i,j,k))*amsuo(i,j,k)
            pyein(i,j)=pyein(i,j)+ddue(i,j,k)*dpye(i,j)*                     &
                 (wpo(i,j,k)-wpo(i,j+1,k))*amsue(i,j,k)
          ENDDO
        ENDDO

!$OMP DO
        DO i=2,ie1
          pxoin(i,j)=pxoin(i,j)*deutio(i,j)
          pyein(i,j)=pyein(i,j)*deutie(i,j)
        ENDDO

#ifdef _PROFILE
  CALL trace_start ('occlit_loop_shasw', 20)
#endif
      ! for k=1 use mask shamsu[oe] (merged from amsu[oe] and shasw[oe])
      k=1
        DO i=2,ie1
          uko1       = (uk1o(i,j,k)+stabn*(pxoin(i,j)                    &
               -dtdxuo(i,j)*(wpo(i+1,j,k)-wpo(i,j,k)))*dduo(i,j,k))*shamsuo(i,j)
          vke1       = (vk1e(i,j,k)+stabn*(pyein(i,j)                    &
               -dpye(i,j)*(wpo(i,j,k)-wpo(i,j+1,k)))*ddue(i,j,k))*shamsue(i,j)

          uko(i,j,k) = uko1
          vke(i,j,k) = vke1

          ucos(i,j)  = ucos(i,j) + uko1
          vcos(i,j)  = vcos(i,j) + vke1

        ENDDO

      ! for k>=2 use usual land sea mask on vector points amsu[oe]
      DO k=2,ke
          DO i=2,ie1
            uko1       = (uk1o(i,j,k)+stabn*(pxoin(i,j)                    &
                 -dtdxuo(i,j)*(wpo(i+1,j,k)-wpo(i,j,k)))*dduo(i,j,k))*amsuo(i,j,k)
            vke1       = (vk1e(i,j,k)+stabn*(pyein(i,j)                    &
                 -dpye(i,j)*(wpo(i,j,k)-wpo(i,j+1,k)))*ddue(i,j,k))*amsue(i,j,k)

            uko(i,j,k) = uko1
            vke(i,j,k) = vke1

            ucos(i,j)  = ucos(i,j) + uko1
            vcos(i,j)  = vcos(i,j) + vke1


          ENDDO
      ENDDO

#ifdef _PROFILE
  CALL trace_stop ('occlit_loop_shasw', 20)
#endif

!$OMP DO
#ifdef _PROFILE
  CALL trace_start ('occlit_loop_ucosdeut', 21)
#endif
        DO k=1,ke
          DO i=2,ie1
            uko(i,j,k)=uko(i,j,k)-ucos(i,j)*deutio(i,j)*dduo(i,j,k)
            vke(i,j,k)=vke(i,j,k)-vcos(i,j)*deutie(i,j)*ddue(i,j,k)
          ENDDO
        ENDDO
#ifdef _PROFILE
  CALL trace_stop ('occlit_loop_ucosdeut', 21)
#endif

      ENDDO ! fused j-loop
#ifdef _PROFILE
  CALL trace_stop ('occlit_loop_pxoin', 24)
#endif

    ELSE ! lwith_barotropic_stokes_drift = .true.

      DO j=2,je1

        DO k=2,ke
          DO i=2,ie1
            pxoin(i,j)=pxoin(i,j)+dduo(i,j,k)*dtdxuo(i,j)*                   &
                 (wpo(i+1,j,k)-wpo(i,j,k))*amsuo(i,j,k)
            pyein(i,j)=pyein(i,j)+ddue(i,j,k)*dpye(i,j)*                     &
                 (wpo(i,j,k)-wpo(i,j+1,k))*amsue(i,j,k)
          ENDDO
        ENDDO
        k=1
        DO i=2,ie1
          pxoin(i,j)=pxoin(i,j)+SURLON(i,j)*dtdxuo(i,j)*                   &
               (wpo(i+1,j,k)-wpo(i,j,k))*amsuo(i,j,k)
          pyein(i,j)=pyein(i,j)+SURLEN(i,j)*dpye(i,j)*                     &
               (wpo(i,j,k)-wpo(i,j+1,k))*amsue(i,j,k)
        ENDDO

        DO i=2,ie1
          pxoin(i,j)=pxoin(i,j)*deution(i,j)
          pyein(i,j)=pyein(i,j)*deutien(i,j)
        ENDDO

        ! for k>=2 use usual land sea mask on vector points amsu[oe]
        DO k=2,ke
          DO i=2,ie1
            uko1       = (uk1o(i,j,k)+stabn*(pxoin(i,j)                    &
                 -dtdxuo(i,j)*(wpo(i+1,j,k)-wpo(i,j,k)))*dduo(i,j,k))*amsuo(i,j,k)
            vke1       = (vk1e(i,j,k)+stabn*(pyein(i,j)                    &
                 -dpye(i,j)*(wpo(i,j,k)-wpo(i,j+1,k)))*ddue(i,j,k))*amsue(i,j,k)

            uko(i,j,k) = uko1
            vke(i,j,k) = vke1

            ucos(i,j)  = ucos(i,j) + uko1
            vcos(i,j)  = vcos(i,j) + vke1

          ENDDO
        ENDDO
       ! for k=1 use mask shamsu[oe] (merged from amsu[oe] and shasw[oe])
        k=1
        DO i=2,ie1
          uko1       = (uk1o(i,j,k)+stabn*(pxoin(i,j)                    &
               -dtdxuo(i,j)*(wpo(i+1,j,k)-wpo(i,j,k)))*surlon(i,j))*shamsuo(i,j)
          vke1       = (vk1e(i,j,k)+stabn*(pyein(i,j)                    &
               -dpye(i,j)*(wpo(i,j,k)-wpo(i,j+1,k)))*surlen(i,j))*shamsue(i,j)
          uko(i,j,k) = uko1
          vke(i,j,k) = vke1

          ucos(i,j)  = ucos(i,j) + uko1
          vcos(i,j)  = vcos(i,j) + vke1
        ENDDO

        DO k=2,ke
          DO i=2,ie1
            uko(i,j,k)=uko(i,j,k)-ucos(i,j)*deution(i,j)*dduo(i,j,k)
            vke(i,j,k)=vke(i,j,k)-vcos(i,j)*deutien(i,j)*ddue(i,j,k)
          ENDDO
        ENDDO
        k=1
        DO i=2,ie1
          uko(i,j,k)=uko(i,j,k)-ucos(i,j)*deution(i,j)*surlon(i,j)
          vke(i,j,k)=vke(i,j,k)-vcos(i,j)*deutien(I,J)*surlen(I,J)
        ENDDO

      ENDDO  ! fused j-loop
    ENDIF  ! lwith_barotropic_stokes_drift

!$OMP SINGLE
#ifdef _PROFILE
  CALL trace_start ('occlit_boundsexchange2', 26)
#endif
    CALL bounds_exch(1,'u',uko,'occlit 7')
    CALL bounds_exch(1,'v',vke,'occlit 8')
!$OMP END SINGLE
#ifdef _PROFILE
  CALL trace_stop ('occlit_boundsexchange2', 26)
#endif


    IF (iter == itermax) CYCLE

!$OMP DO
    tlow(:,:)    = 0._wp
    slow(:,:)    = 0._wp
    contrij(:,:) = 0._wp

!$OMP DO
#ifdef _PROFILE
  CALL trace_start ('occlit_loop_t1o', 23)
#endif
      DO k=ke,1,-1
       DO j=je1,2,-1
        DO i=ie1,2,-1

!          t1o(i,j,k)= slow(i,j) + dti * weto(i,j,k) * (                    &
!               dtdxpo(i,j) * ( uko(i-1,j,k)*dlyu(i-1,j)                    &
!                              -uko(i,j,k)  *dlyu(i,j)  )/dlyp(i,j)         &
!             + dtdyo(i,j)  * ( vke(i,j,k)  *dlxv(i,j)                      &
!                              -vke(i,j-1,k)*dlxv(i,j-1))/dlxp(i,j)  )

!emr new formulation
          t1o(i,j,k)= slow(i,j) + weto(i,j,k)* (                    &
       &          uko(i-1,j,k)*dlyu(i-1,j)                    &
                 -uko(i,j,k)  *dlyu(i,j)         &
             +   vke(i,j,k)  *dlxv(i,j)                      &
       &      -vke(i,j-1,k)*dlxv(i,j-1))/(dlyp(i,j)*dlxp(i,j)  )


          slow(i,j)=t1o(i,j,k)
        ENDDO
      ENDDO
    ENDDO
#ifdef _PROFILE
  CALL trace_stop ('occlit_loop_t1o', 23)
#endif

    IF (iter == 1) THEN
      under = 0.0_wp
    ELSE
      under = 1.0_wp
    ENDIF

    jb=2
    IF(p_joff.EQ.0 .AND. lbounds_exch_tp ) jb=3

!$OMP DO
#ifdef _PROFILE
  CALL trace_start ('occlit_loop_contr', 22)
#endif
      DO k=2,ke
       DO j=jb,je1
        DO i=2,ie1
          ove=wpo(i,j,k)
          tlow(i,j)=tlow(i,j)+g*stabio(i,j,k)*dt*(conn*t1o(i,j,k)+wo(i,j,k))
!          speed=2.*g*bruva(i,j)*(dtdxuo(i,j)**2+dtdyo(i,j)**2)
!          brems=under*speed*conn*stabn/(1._wp+speed)
          speed=g*bruva(i,j)*(dtdxuo(i,j)**2+dtdyo(i,j)**2)
          brems = 2._wp * under * speed * conn * stabn/(1._wp + speed)
          vorw = 1._wp - brems
          wpo(i,j,k)=vorw*tlow(i,j)+brems*wpo(i,j,k)
          contrij(i,j)=contrij(i,j)+(ove-wpo(i,j,k))**2
        ENDDO
      ENDDO
    ENDDO
#ifdef _PROFILE
  CALL trace_stop ('occlit_loop_contr', 22)
#endif

!$OMP SINGLE
#ifdef _PROFILE
  CALL trace_start ('occlit_boundsexchange3/gather', 27)
#endif
    CALL bounds_exch(1,'p',wpo,'occlit 10')
    IF (icontro.NE.0) THEN
      CALL gather(contrij,contrij_g,p_io)
    ENDIF
#ifdef _PROFILE
  CALL trace_stop ('occlit_boundsexchange3/gather', 27)
#endif
!$OMP END SINGLE

    IF (icontro.NE.0) THEN
      IF(p_pe==p_io) THEN
        contra = 0._wp
        DO  j=2,je_g-1
          DO  i=2,ie_g-1
            contra=contra+contrij_g(i,j)
          ENDDO
        ENDDO
        WRITE(0,*)'ITER CONTRA: ',iter,contra,MAXLOC(contrij_g)
      ENDIF
    ENDIF

  END DO iteration_loop

!$OMP WORKSHARE

  uko(:,:,2:)=uko(:,:,2:)*amsuo(:,:,2:)/(almzer+dduo(:,:,2:))
  vke(:,:,2:)=vke(:,:,2:)*amsue(:,:,2:)/(almzer+ddue(:,:,2:))
  uk1o(:,:,2:)=uk1o(:,:,2:)*amsuo(:,:,2:)/(almzer+dduo(:,:,2:))
  vk1e(:,:,2:)=vk1e(:,:,2:)*amsue(:,:,2:)/(almzer+ddue(:,:,2:))

  IF ( .NOT. lwith_barotropic_stokes_drift ) THEN

    uko(:,:,1)=uko(:,:,1)*amsuo(:,:,1)/(almzer+dduo(:,:,1))
    vke(:,:,1)=vke(:,:,1)*amsue(:,:,1)/(almzer+ddue(:,:,1))
    uk1o(:,:,1)=uk1o(:,:,1)*amsuo(:,:,1)/(almzer+dduo(:,:,1))
    vk1e(:,:,1)=vk1e(:,:,1)*amsue(:,:,1)/(almzer+ddue(:,:,1))

  ELSE

    uko(:,:,1)=uko(:,:,1)*amsuo(:,:,1)/(almzer+surlon(:,:))
    vke(:,:,1)=vke(:,:,1)*amsue(:,:,1)/(almzer+surlen(:,:))
    uk1o(:,:,1)=uk1o(:,:,1)*amsuo(:,:,1)/(almzer+surlon(:,:))
    vk1e(:,:,1)=vk1e(:,:,1)*amsue(:,:,1)/(almzer+surlen(:,:))

  ENDIF

!$OMP END WORKSHARE



!$OMP END PARALLEL
  !
#ifdef _PROFILE
  CALL trace_start ('occlit_lastgather', 28)
#endif
  CALL gather(contrij,contrij_g,p_io)
#ifdef _PROFILE
  CALL trace_stop ('occlit_lastgather', 28)
#endif

  IF(p_pe==p_io) THEN
    contra = 0._wp
    DO  j=2,je_g-1
      DO  i=2,ie_g-1
        contra=contra+contrij_g(i,j)
      ENDDO
    ENDDO
    !
    WRITE(io_stdout,*)'CONTRA: ',contra
    !
    IF (contra .GT. 1.e-3_wp) THEN

       id=maxloc(contrij_g(2:ie_g-1,2:je_g-1))


      WRITE(io_stdout,*) 'CONTR : ',id(1),id(2)
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) 'CONTRI : ',SUM(contrij_g(2:ie_g-1,2:je_g-1),1)
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) 'CONTRJ : ',SUM(contrij_g(2:ie_g-1,2:je_g-1),2)
      CALL stop_all('Error - No proper solution of the baroclinic subsystem => run aborted ')
    ENDIF
  ENDIF

#ifdef _PROFILE
  CALL trace_stop ('occlit', 1)
#endif



END SUBROUTINE occlit



