PROGRAM DEWPOINT

  USE mo_kind, ONLY: i4, sp
  !     CALCULATES TDEW FROM SPECIFIC HUMIDITY AND PRESSURE

  INTEGER(i4), PARAMETER :: IE=192,JE=94

  REAL(sp)    :: SHUM(IE,JE),PRESS(IE,JE),TDEW(IE,JE),UWIND(IE,JE), &
                 VWIND(IE,JE),UVWIND(IE,JE),ETEST(IE,JE),QTEST(IE,JE)
  INTEGER(i4) :: ext4_hdr(4)

  OPEN(10,FILE='shum.2m',FORM='UNFORMATTED')
  OPEN(20,FILE='pres.sfc',FORM='UNFORMATTED')
  OPEN(30,FILE='tdew.2m',FORM='UNFORMATTED')
  OPEN(40,FILE='uwnd.10m',FORM='UNFORMATTED')
  OPEN(50,FILE='vwnd.10m',FORM='UNFORMATTED')
  OPEN(60,FILE='uvwnd.10m',FORM='UNFORMATTED')

  READ(*,*)LDAY
  DO L=1,LDAY

    READ(10) ext4_hdr
    READ(10) SHUM
    READ(20) ext4_hdr
    READ(20) PRESS
    READ(40) ext4_hdr
    READ(40) UWIND
    READ(50) ext4_hdr
    READ(50) VWIND

    DO I=1,IE
      DO J=1,JE
        EQP=(SHUM(I,J)*PRESS(I,J))/(0.378*SHUM(I,J)+0.6222)
        ALPHA=(1./7.5)*LOG10(EQP/611.)
        TDEW(I,J)=(273.16-(35.86*alpha))/(1.-alpha)
        UVWIND(I,J)=SQRT(UWIND(I,J)**2+VWIND(I,J)**2)
      ENDDO
    ENDDO

    !       TEST
    DO I=1,IE
      DO J=1,JE

        ETEST(I,J)=611.*10.**(((TDEW(I,J)-273.16) &
             /(TDEW(I,J)-35.86))*7.5)
        QTEST(I,J)=(0.6222*ETEST(I,J)) &
             /(PRESS(I,J)-0.378*ETEST(I,J))
      ENDDO
    ENDDO

    WRITE(70) ext4_hdr(1),51,ext4_hdr(3:4)
    WRITE(70) QTEST
    WRITE(30) ext4_hdr(1),17,ext4_hdr(3:4)
    WRITE(30) TDEW
    WRITE(60) ext4_hdr(1),32,ext4_hdr(3:4)
    WRITE(60) UVWIND

  ENDDO

  CLOSE(10)
  CLOSE(20)
  CLOSE(30)
  CLOSE(40)
  CLOSE(50)
  CLOSE(60)

END PROGRAM DEWPOINT
