
      SUBROUTINE ROTATE_U(FIELDU,FIELDV,II,JJ)
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1, ONLY: giph, gila, amsuo
      USE mo_constants, ONLY: api
      USE MO_UNITS, ONLY: io_stdout

      IMPLICIT NONE
      INTEGER, INTENT(in) :: ii, jj
      REAL(wp) FIELDU(II,JJ),FIELDV(II,JJ)

      INTEGER :: i, j
      REAL(wp) :: uout, vout, absin, absout, uin, vin, proms, &
           deltxx, deltxy, deltyx, deltyy
!TEST   WRITE(IO_STDOUT,*) 'IN SBR ROTATE_U'
      DO I=2,IE1
       DO J=2,JE1
        UIN=FIELDU(I,J)
        VIN=FIELDV(I,J)
        uout = 0.0_wp
        vout = 0.0_wp
!TEST   WRITE(IO_STDOUT,*) 'I,J,UIN,VIN',I,J,UIN,VIN
        IF (amsuo(i, j, 1) .GT. 0.5_wp) THEN
        DELTXX=GILA((2*I)+1,2*J) - GILA(2*I,2*J)
        DELTYY=GIPH((2*I)+1,(2*J)-1)-GIPH((2*I)+1,2*J)
        DELTXY=GIPH((2*I)+1,2*J) - GIPH(2*I,2*J)
        DELTYX=GILA((2*I)+1,(2*J)-1)-GILA((2*I)+1,2*J)
        IF (deltxx .GT. api) deltxx = deltxx - 2._wp * api
        IF (deltxx .LT. -api) deltxx = deltxx + 2._wp * api
        IF (deltyx .LT. -api) deltyx = deltyx + 2._wp * api
        IF (deltyx .GT. api) deltyx = deltyx - 2._wp * api
        PROMS=COS(GIPH((2*I)+1,2*J))
        DELTXX=PROMS*DELTXX
        DELTYX=PROMS*DELTYX
!:: PROJECTION BY SCALAR PRODUCT
        UOUT= UIN*DELTXX + VIN*DELTXY
        VOUT= UIN*DELTYX + VIN*DELTYY
!TEST   WRITE(IO_STDOUT,*) 'I,J,UOUT,VOUT',I,J,UOUT,VOUT
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
        IF (absout .GT. 1.E-10_wp) THEN
         UOUT=UOUT*ABSIN/ABSOUT
!         UOUT=UOUT
        ELSE
         uout = 0.0_wp
        ENDIF
        ENDIF
        FIELDU(I,J)=UOUT
       ENDDO
      ENDDO
      RETURN
      END
