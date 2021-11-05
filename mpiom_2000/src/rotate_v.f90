      SUBROUTINE ROTATE_V(FIELDU,FIELDV,II,JJ)
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1, ONLY: giph, gila, amsue
      USE mo_constants, ONLY: api
      USE MO_UNITS
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ii, jj
      REAL(wp) FIELDU(II,JJ),FIELDV(II,JJ)
      INTEGER :: i, j
      REAL(wp) :: absin, absout, deltxx, deltxy, deltyx, deltyy, proms, &
           uin, uout, vin, vout
      DO J=2,JE1
       DO I=2,IE1
        UIN=FIELDU(I,J)
        VIN=FIELDV(I,J)
        uout = 0.0_wp
        vout = 0.0_wp
        IF (amsue(i, j, 1) .GT. 0.5_wp) THEN
        DELTXX=GILA((2*I),(2*J)+1) - GILA((2*I)-1,(2*J)+1)
        DELTYY=GIPH(2*I,2*J) -GIPH(2*I,(2*J)+1)
        DELTXY=GIPH((2*I),(2*J)+1) - GIPH((2*I)-1,(2*J)+1)
        DELTYX=GILA(2*I,2*J)-GILA(2*I,(2*J)+1)
!::
        IF (deltxx .GT. api) deltxx = deltxx - 2._wp * api
        IF (deltxx .LT. -api) deltxx = deltxx + 2._wp * api
        IF (deltyx .LT. -api) deltyx = deltyx + 2._wp * api
        IF (deltyx .GT. api) deltyx = deltyx - 2._wp * api
!
        PROMS=COS( GIPH(2*I,(2*J)+1) )
        DELTXX=PROMS*DELTXX
        DELTYX=PROMS*DELTYX
!:: PROJECTION BY SCALAR PRODUCT
        UOUT= UIN*DELTXX + VIN*DELTXY
        VOUT= UIN*DELTYX + VIN*DELTYY
!:: TEST ORTHOGONALITY
        IF (deltxx * deltyx + deltyy * deltxy .GT. 0.1_wp) THEN
         WRITE(IO_STDOUT,*)'ORTHOGONAL? ',I,J,                          &
     &       (DELTXX*DELTYX+DELTYY*DELTXY)
        ENDIF
        UOUT=UOUT/SQRT(DELTXX**2+DELTXY**2)
        VOUT=VOUT/SQRT(DELTYY**2+DELTYX**2)
!:: CORRECT SKALE
        ABSIN=SQRT(UIN**2+VIN**2)
        ABSOUT=SQRT(UOUT**2+VOUT**2)
        IF (absout .GT. 1.e-10_wp) THEN
         VOUT=VOUT*ABSIN/ABSOUT
        ELSE
         vout = 0.0_wp
        ENDIF
        ENDIF
!JJ ENDIF AMSUE
        FIELDV(I,J)=VOUT
       ENDDO
      ENDDO
      RETURN
      END
