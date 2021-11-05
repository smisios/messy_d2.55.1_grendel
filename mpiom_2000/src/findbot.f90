!>
!! Finds deepest wet level on pressure point
!!
!! @author JHJ
!! @date Sept. 2, 1999
!!
SUBROUTINE findbot(kbot,weto,ie,je,ke)

  USE mo_kind, ONLY: wp
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ie, je, ke
  INTEGER, INTENT(INOUT) :: kbot(ie,je)
  REAL(wp), INTENT(IN) :: weto(ie,je,ke)

  INTEGER :: i, j, k

  iloop: DO i=1,ie
     jloop: DO j=1,je
        IF (weto(i,j,1) .LT. 0.5_wp) THEN
           kbot(i,j)=0
           CYCLE jloop
        ENDIF
        kbot(i,j)=ke
        kloop: DO k=2,ke
           IF (weto(i,j,k) .LT. 0.5_wp) THEN
              kbot(i,j)=k-1
              CYCLE jloop
           ENDIF
        END DO kloop
     END DO jloop
  END DO iloop
  RETURN
END SUBROUTINE findbot
