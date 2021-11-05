      SUBROUTINE BARTIM
!*********************************************************
! UNTERPROGR. ZUR BERECHNUNG DES WASSERSTANDES
! BAROTROPES SYSTEM (IMPLIZITES VERFAHREN,
! BEWEGUNGSGL. UND KONTI. GL.)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      USE mo_kind, ONLY: wp
      USE MO_PARAM1
      USE mo_para2
      USE MO_MPI
      USE MO_PARALLEL
      USE mo_boundsexch, ONLY : bounds_exch
      USE MO_ELICOM
      USE MO_COMMO1
      USE MO_UNITS

      IMPLICIT NONE

      INTEGER I, J, L, LL
      REAL(wp), ALLOCATABLE :: B1O_G(:,:)
      REAL(wp) BBMM, XR
!      integer, save :: ic = 0
!
      IF(p_pe==0) THEN
         ALLOCATE (B1O_G(IE_G,JE_G))
      ELSE
         ALLOCATE (B1O_G(0,0))
      ENDIF

      Z1O(:,:) = 0._wp
      B1O(:,:) = 0._wp

      DO J=2,JE1
        DO I=2,IE1
          B1O(I,J) = WETO(I,J,1) * (                                        &
         & DTDXPO(I,J) * ((CONO*U1O(I-1,J) +  UZO(I-1,J)*CONN)*DLYU(I-1,J)  &
         &       - (U1O(I,J)*CONO  +  UZO(I,J)*CONN)*DLYU(I,J) )/DLYP(I,J)  &
         & + DTDYO(I,J) * (                                                 &
         &  ( V1E(I,J) *CONO  +  CONN*VZE(I,J))*DLXV(I,J)                   &
         &- ( V1E(I,J-1) *CONO+CONN*VZE(I,J-1))*DLXV(I,J-1))/DLXP(I,J)  )

          B1O(I,J) = B1O(I,J)*DLXP(I,J)*DLYP(I,J)
        ENDDO
      ENDDO

      CALL gather(B1O,B1O_G,0)

      IF(p_pe==0) THEN

!      DO 2030 J=1,JE_G
!      DO 2030 I=2,IE_G-1
!      IF(NUM_G(I,J).EQ.0) GO TO 2030
!      L=NUM_G(I,J)
!      B(L)=B1O_G(I,J)
!2030  CONTINUE
!
!if (ic == 1) then
!write (1) SIZE(b)
!write (1) b
!endif

      jloop: DO J=1,JE_G
        iloop: DO I=2,IE_G-1
          IF(NUM_G(I,J).EQ.0) CYCLE iloop
          L=NUM_G(I,J)
          B(L)=B1O_G(I,J) *SKAL(L)
        END DO iloop
      END DO jloop

      DO J=1,MATR-1
        BBMM=B(J)
        DO I=1,KBB
          B(I+J)=B(I+J)-PGL(I,J)*BBMM
        END DO
      END DO

      X(MATR)=B(MATR)/PGL(KM,MATR)

      DO LL=1,MATR-1
        L=MATR-LL
        XR=B(L)
        DO J=1,KBB
          XR=XR-PGL(KM+J,L)*X(L+J)
        END DO
        X(L)=XR/PGL(KM,L)
      END DO

      ENDIF

      CALL p_bcast(X,0)

      DO J=1,JE
        iiloop: DO I=2,IE1
          IF(NUM(I,J).EQ.0) CYCLE iiloop
          L=NUM(I,J)
          Z1O(I,J)=X(L)
        END DO iiloop
      END DO
!
      CALL bounds_exch(1,'p',Z1O,'bartim 1')
!if (ic == 0) then
!open(1,file='rhs-guess.dat',form='unformatted',access='sequential')
!write (0,*) '##################### X',  SIZE(x)
!write (1) SIZE(x)
!write (1) x
!endif
!if (ic == 1) then
!write (1) SIZE(x)
!write (1) x
!close(1)
!stop
!endif

!ic = 1
      RETURN
      END



