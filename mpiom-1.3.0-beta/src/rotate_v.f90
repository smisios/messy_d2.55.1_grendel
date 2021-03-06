      SUBROUTINE ROTATE_V(FIELDU,FIELDV,II,JJ)
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS    
      REAL FIELDU(II,JJ),FIELDV(II,JJ)
      PI=ATAN(1.)*4.0
      DO J=2,JE1
       DO I=2,IE1
        UIN=FIELDU(I,J)
        VIN=FIELDV(I,J)
        UOUT=0.0
        VOUT=0.0
        IF (AMSUE(I,J,1) .GT. 0.5) THEN
        DELTXX=GILA((2*I),(2*J)+1) - GILA((2*I)-1,(2*J)+1)
        DELTYY=GIPH(2*I,2*J) -GIPH(2*I,(2*J)+1)
        DELTXY=GIPH((2*I),(2*J)+1) - GIPH((2*I)-1,(2*J)+1)
        DELTYX=GILA(2*I,2*J)-GILA(2*I,(2*J)+1)
!::
        IF (DELTXX .GT. PI) DELTXX=DELTXX-2.*PI
        IF (DELTXX .LT. -PI) DELTXX=DELTXX+2.*PI
        IF (DELTYX .LT. -PI) DELTYX=DELTYX+2.*PI
        IF (DELTYX .GT. PI) DELTYX=DELTYX-2.*PI
!
        PROMS=COS( GIPH(2*I,(2*J)+1) )
        DELTXX=PROMS*DELTXX
        DELTYX=PROMS*DELTYX
!:: PROJECTION BY SCALAR PRODUCT
        UOUT= UIN*DELTXX + VIN*DELTXY
        VOUT= UIN*DELTYX + VIN*DELTYY
!:: TEST ORTHOGONALITY
        IF ((DELTXX*DELTYX+DELTYY*DELTXY) .GT.0.1) THEN
         WRITE(IO_STDOUT,*)'ORTHOGONAL? ',I,J,                          &
     &       (DELTXX*DELTYX+DELTYY*DELTXY)
        ENDIF
        UOUT=UOUT/SQRT(DELTXX**2+DELTXY**2)
        VOUT=VOUT/SQRT(DELTYY**2+DELTYX**2)
!:: CORRECT SKALE
        ABSIN=SQRT(UIN**2+VIN**2)
        ABSOUT=SQRT(UOUT**2+VOUT**2)
        IF (ABSOUT .GT. 1.E-10) THEN
         VOUT=VOUT*ABSIN/ABSOUT
        ELSE
         VOUT=0.0
        ENDIF
        ENDIF
!JJ ENDIF AMSUE
        FIELDV(I,J)=VOUT
       ENDDO
      ENDDO
      RETURN
      END
