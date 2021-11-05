
SUBROUTINE OCCLIT2

  USE MO_PARAM1
  USE mo_planetary_constants, ONLY: g
  USE MO_PARALLEL
  USE mo_boundsexch, ONLY : bounds_exch
  USE MO_COMMO1
  USE MO_UNITS

#ifdef _PROFILE
  USE mo_profile,      ONLY: trace_start, trace_stop
#endif


  IMPLICIT NONE

  INTEGER :: i,j,k,iter,itermax, jb, id(2)
  INTEGER :: JS,JJ,JBLOCK,JSTOP
  REAL(wp) :: brems, contra, ove, speed, under, vorw, uko1, vke1
  !
  REAL(wp) :: bruva(ie,je),contrij(ie,je)
  REAL(wp), POINTER :: contrij_g(:,:)
  !
  REAL(wp) :: wpo(ie,je,ke) ! changed from kep to ke by LK
  !
  !  ITERATION OF BAROCLINIC SYSTEM DIRECTLY IN BAROCLINIC VELOCITIES

#ifdef _PROFILE
  CALL trace_start ('occlit2', 1)
#endif

  jb=2
  if( p_joff.eq.0 .and. lbounds_exch_tp ) jb=3


  if (p_pe==p_io) then
     ALLOCATE(contrij_g(ie_g,je_g))
  else
     contrij_g => NULL()
  endif


!$OMP PARALLEL PRIVATE(i,j,k,iter,itermax,under,uko1,vke1,ove,speed,brems,vorw)

  JBLOCK=MIN(16,JE)

  DO JS=1,JE,JBLOCK

     JJ=JS
     JSTOP=MIN(JS+JBLOCK-1,JE)

     DO K=KE,1,-1
        IF (K.eq.1)THEN
           DO J=JJ,JSTOP
              DO I=1,IE
                 bruva(i, j) = 0._wp
                 stabio(i, j, 1) = 0._wp
                 slow(i, j) = 0._wp
              ENDDO
           ENDDO
        ENDIF
        DO J=JS,JSTOP
           DO I=1,IE
              UKO(I,J,K)=UKO(I,J,K)*DDUO(I,J,K)
              VKE(I,J,K)=VKE(I,J,K)*DDUE(I,J,K)
              UK1O(I,J,K)=UK1O(I,J,K)*DDUO(I,J,K)
              VK1E(I,J,K)=VK1E(I,J,K)*DDUE(I,J,K)
           ENDDO
        ENDDO

        JJ=MAX(2,JS)
        JSTOP=MIN(JS+JBLOCK-1,JE1)

        DO J=JJ,JSTOP
           DO I=2,IE1
!           WO(I,J,K) = WO(I,J,K+1) + DTI * WETO(I,J,K) * (                   &
!                DTDXPO(I,J) * (   UKO(I-1,J,K)*DLYU(I-1,J)                 &
!                - UKO(I,J,K)*DLYU(I,J)    )/DLYP(I,J)      &
!                + DTDYO(I,J)  * (   VKE(I,J,K)*DLXV(I,J)                     &
!                - VKE(I,J-1,K)*DLXV(I,J-1))/DLXP(I,J)  )

! needs to be checked in a longer run

              WO(I,J,K)= WO(I,J,K+1) + weto(i,j,k)* areain(I,J)*       &
                (UKO(I-1,J,K)*DLYU(I-1,J)                          &
                 -UKO(I,J,K)  *DLYU(I,J)                           &
             +   VKE(I,J,K)  *DLXV(I,J)                            &
                -VKE(I,J-1,K)*DLXV(I,J-1))
           ENDDO
        ENDDO
     ENDDO
  ENDDO


!$OMP SINGLE
!     CALL bounds_exch(1,'p',WO(:,:,K),'occlit 4')
     CALL bounds_exch(1,'p',WO,'occlit 4')
!$OMP END SINGLE

     DO K=1,KE
        DO J=1,JE
           DO I=1,IE
              BRUVA(I,J)=BRUVA(I,J)+TIESTW(K+1)*STABIO(I,J,K)
              SLOW(I,J)=SLOW(I,J)+STABIO(I,J,K)*DT*WO(I,J,K)
              WPO(I,J,K)=SLOW(I,J)
           ENDDO
        ENDDO
     ENDDO


  !OtB Needs more iterations than 6
  !     ITERMAX=6
  !     ITERMAX=16
  !Uwe says 8 is cheaper
  ITERMAX=8

  if ( lbounds_exch_tp ) ITERMAX=12

  iteration_loop: DO ITER=1,ITERMAX

    IF (iter == 1) THEN
      under = 0.0_wp
    ELSE
      under = 1.0_wp
    ENDIF

    !
!$OMP DO
    DO J=1,JE
      DO I=1,IE
        ucos(i, j) = 0._wp
        vcos(i, j) = 0._wp
        pxoin(i, j) = 0._wp
        pyein(i, j) = 0._wp
      ENDDO
    ENDDO

!$OMP DO
    DO K=1,KE
       DO J=2,JE1
          DO I=2,IE1
             PXOIN(I,J)=PXOIN(I,J)+DDUO(I,J,K)*DTDXUO(I,J)*                   &
                  (WPO(I+1,J,K)-WPO(I,J,K))*AMSUO(I,J,K)
             PYEIN(I,J)=PYEIN(I,J)+DDUE(I,J,K)*DPYE(I,J)*                     &
                  (WPO(I,J,K)-WPO(I,J+1,K))*AMSUE(I,J,K)
          ENDDO
       ENDDO
    ENDDO

!$OMP DO
    DO J=2,JE1
       DO I=2,IE1
          PXOIN(I,J)=PXOIN(I,J)*DEUTIO(I,J)
          PYEIN(I,J)=PYEIN(I,J)*DEUTIE(I,J)
       ENDDO
    ENDDO

!$OMP DO
    DO K=1,KE
       DO J=2,JE1
          DO I=2,IE1
             uko1       = uk1o(i,j,k)+stabn*(pxoin(i,j)                    &
                  -dtdxuo(i,j)*(wpo(i+1,j,k)-wpo(i,j,k)))*dduo(i,j,k)
             vke1       = vk1e(i,j,k)+stabn*(pyein(i,j)                    &
                  -dpye(i,j)*(wpo(i,j,k)-wpo(i,j+1,k)))*ddue(i,j,k)
             uko(i,j,k) = uko1*amsuo(i,j,k)
             vke(i,j,k) = vke1*amsue(i,j,k)

             ucos(i,j)  = ucos(i,j)+uko(i,j,k)
             vcos(i,j)  = vcos(i,j)+vke(i,j,k)

          ENDDO
       ENDDO
    ENDDO

!$OMP DO
    DO K=1,KE
       DO J=2,JE1
          DO I=2,IE1
          UKO(I,J,K)=UKO(I,J,K)-UCOS(I,J)*DEUTIO(I,J)*DDUO(I,J,K)
          VKE(I,J,K)=VKE(I,J,K)-VCOS(I,J)*DEUTIE(I,J)*DDUE(I,J,K)
        ENDDO
      ENDDO
    ENDDO

!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch(1,'u',UKO,'occlit 7')
    CALL bounds_exch(1,'v',VKE,'occlit 8')
!$OMP END SINGLE
!#endif
    !
    IF (ITER == ITERMAX) CYCLE
    !
!$OMP DO
    DO J=1,JE
      DO I=1,IE
        tlow(i, j) = 0._wp
        slow(i, j) = 0._wp
        contrij(i, j) = 0._wp
      ENDDO
    ENDDO
    !
!$OMP DO
      DO K=KE,1,-1
         DO J=2,JE1
            DO I=2,IE1

!          T1O(I,J,K)= SLOW(I,J) + DTI * WETO(I,J,K) * (                    &
!               DTDXPO(I,J) * ( UKO(I-1,J,K)*DLYU(I-1,J)                    &
!                              -UKO(I,J,K)  *DLYU(I,J)  )/DLYP(I,J)         &
!             + DTDYO(I,J)  * ( VKE(I,J,K)  *DLXV(I,J)                      &
!                              -VKE(I,J-1,K)*DLXV(I,J-1))/DLXP(I,J)  )

! new formulation
          T1O(I,J,K)= SLOW(I,J) + weto(i,j,k)* areain(I,J)*       &
                (UKO(I-1,J,K)*DLYU(I-1,J)                          &
                 -UKO(I,J,K)  *DLYU(I,J)                           &
             +   VKE(I,J,K)  *DLXV(I,J)                            &
                -VKE(I,J-1,K)*DLXV(I,J-1))

!          T1O(I,J,K)= SLOW(I,J) + weto(i,j,k)* (                    &
!       &          UKO(I-1,J,K)*DLYU(I-1,J)                    &
!                 -UKO(I,J,K)  *DLYU(I,J)         &
!             +   VKE(I,J,K)  *DLXV(I,J)                      &
!       &      -VKE(I,J-1,K)*DLXV(I,J-1))/(dlyp(i,j)*DLXP(I,J)  )

          SLOW(I,J)=T1O(I,J,K)
        ENDDO
      ENDDO
    ENDDO

!$OMP DO
      DO K=2,KE
         DO J=jb,JE1
            DO I=2,IE1
          OVE=WPO(I,J,K)
          TLOW(I,J)=TLOW(I,J)+G*STABIO(I,J,K)*DT*(CONN*T1O(I,J,K)+WO(I,J,K))
!          SPEED=2.*G*BRUVA(I,J)*(DTDXUO(I,J)**2+DTDYO(I,J)**2)
!          BREMS=UNDER*SPEED*CONN*STABN/(1.+SPEED)
          SPEED=G*BRUVA(I,J)*(DTDXUO(I,J)**2+DTDYO(I,J)**2)
          brems = 2._wp * under * speed * conn * stabn/(1._wp + speed)
          vorw = 1._wp - brems
          WPO(I,J,K)=VORW*TLOW(I,J)+BREMS*WPO(I,J,K)
          CONTRIJ(I,J)=CONTRIJ(I,J)+(OVE-WPO(I,J,K))**2
        ENDDO
      ENDDO
    ENDDO

!$OMP SINGLE

    CALL bounds_exch(1,'p',WPO,'occlit 10')

    if (icontro.ne.0) then
       CALL gather(contrij,contrij_g,p_io)
    endif
!$OMP END SINGLE

    if (icontro.ne.0) then
       IF(p_pe==p_io) THEN
          contra = 0._wp
          DO  J=2,JE_G-1
             DO  I=2,IE_G-1
                CONTRA=CONTRA+CONTRIJ_G(I,J)
             ENDDO
          ENDDO
          WRITE(0,*)'ITER CONTRA: ',iter,CONTRA,maxloc(CONTRIJ_G)
       ENDIF
    endif

  END DO iteration_loop

!$OMP DO
  DO K=1,KE
    DO J=1,JE
      DO I=1,IE
        uko(i,j,k)=uko(i,j,k)*amsuo(i,j,k)/(almzer+dduo(i,j,k))
        vke(i,j,k)=vke(i,j,k)*amsue(i,j,k)/(almzer+ddue(i,j,k))
        uk1o(i,j,k)=uk1o(i,j,k)*amsuo(i,j,k)/(almzer+dduo(i,j,k))
        vk1e(i,j,k)=vk1e(i,j,k)*amsue(i,j,k)/(almzer+ddue(i,j,k))
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL
  !

!  if (icontro.ne.0) then

  CALL gather(contrij,contrij_g,p_io)

  IF(p_pe==p_io) THEN
    contra = 0._wp
    DO  J=2,JE_G-1
      DO  I=2,IE_G-1
        CONTRA=CONTRA+CONTRIJ_G(I,J)
      ENDDO
    ENDDO
    !
    WRITE(IO_STDOUT,*)'CONTRA: ',CONTRA
    !
    IF (contra .GT. 1.e-3_wp) THEN

       id=maxloc(CONTRIJ_G(2:IE_G-1,2:JE_G-1))


      WRITE(IO_STDOUT,*) 'CONTR : ',id(1),id(2)
      WRITE(IO_STDOUT,*) ' '
      WRITE(IO_STDOUT,*) 'CONTRI : ',SUM(CONTRIJ_G(2:IE_G-1,2:JE_G-1),1)
      WRITE(IO_STDOUT,*) ' '
      WRITE(IO_STDOUT,*) 'CONTRJ : ',SUM(CONTRIJ_G(2:IE_G-1,2:JE_G-1),2)
      CALL STOP_ALL('Error - No proper solution of the baroclinic subsystem => run aborted ')
    ENDIF
  ENDIF

!ENDIF

#ifdef _PROFILE
  CALL trace_stop ('occlit2', 1)
#endif



END SUBROUTINE OCCLIT2



