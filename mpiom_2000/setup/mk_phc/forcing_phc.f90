PROGRAM FORCING
  USE iso_varying_string
  USE mo_kind, ONLY: dp, sp, i4, i8, wp
  USE mo_constants, ONLY: api, aradtogra
  IMPLICIT NONE
  INTERFACE
    FUNCTION gausslo(i)
      USE mo_kind
      INTEGER, INTENT(in) :: i
      REAL(sp) gausslo
    END FUNCTION gausslo
  END INTERFACE
  INTERFACE
    FUNCTION gaussla(j)
      USE mo_kind
      INTEGER, INTENT(in) :: j
      REAL(sp) gaussla
    END FUNCTION gaussla
  END INTERFACE
  INTEGER, PARAMETER :: IE=360, JE=180
  INTEGER, PARAMETER :: io_in_anta=81, io_in_forcing=20, &
       io_out_diag=54, io_out_forcing=55
  INTEGER :: idepth_lvls, ME, NE, KE

  REAL(sp) GAULAT(IE+2, JE+2), GAULON(IE+2, JE+2)
  REAL(sp) LAND(IE+2, JE+2), DUMMY(IE, JE)
  !  REAL(sp) GAUSSLA(JE), GAUSSLO(IE)
  REAL(sp), ALLOCATABLE :: GILA(:, :), GIPH(:, :), &
       HOPLAT(:, :), HOPLON(:, :)
  REAL(sp) :: S(ie, je), pt(ie, je)
  REAL(sp) :: FIIN(IE+2, JE+2, 1)
  REAL(sp), ALLOCATABLE :: FIOUT(:, :, :)
  REAL(wp), ALLOCATABLE :: TINTER(:, :, :)
  REAL(sp) :: TEMIN(IE+2, JE+2)
  REAL(dp), ALLOCATABLE :: TEMOUT(:, :), TOUTER(:, :)
  REAL(sp), ALLOCATABLE :: ZZOUT(:)

  REAL(sp) :: tiefe, dzl, dzs, p, sumtiefe
  REAL(wp) :: alpha

  INTEGER(i8) :: output_hdr(4)
  INTEGER(i4) :: ext4_hdr(4)
  INTEGER lday, iday, if3up, i, j, k, lev, m, n, verbose
  CHARACTER(len=4) :: code
  LOGICAL :: periodic_bound
  NAMELIST /gridparams/ idepth_lvls, me, ne, ke, code, lday,&
       periodic_bound, verbose
  INTEGER :: ioerrstat
  TYPE(varying_string) :: diag_out, output_fname, &
       input_fname, input1_fname

  idepth_lvls = -1
  me = -1
  ne = -1
  ke = -1
  lday = -1
  code = ''
  periodic_bound = .false.
  verbose = 0
  READ (*, nml=gridparams)
  IF (verbose .GT. 0) PRINT *, 'idepth_lvls=', idepth_lvls, &
       'me=', me, 'ne=', ne, 'ke=', ke, 'lday=', lday, 'code=', code
  IF (idepth_lvls < 1 .OR. me < 1 .OR. ne < 1 .OR. ke < 1 .or. lday < 1) THEN
    PRINT *, &
     'Only positive integral values may be used for idepth_lvls, me, ne and ke!'
    STOP 1
  END IF
  IF (code .EQ. 'temp') THEN
    output_hdr(2) = 2_i8
  ELSE IF (code .EQ. 'salt') THEN
    output_hdr(2) = 5_i8
  ELSE
    PRINT *, 'code must be one of (/''temp'', ''salt''/).'
    STOP 1
  END IF

  ALLOCATE(GILA(2*ME, 2*NE), GIPH(2*ME, 2*NE), HOPLAT(ME, NE), &
       HOPLON(ME, NE), FIOUT(ME,NE,1), TINTER(ME,NE,idepth_lvls), &
       TEMOUT(ME, NE), TOUTER(ME, NE))
  ALLOCATE(zzout(ke))
  sumtiefe = 0._sp
  READ *, zzout
  DO k = 1, ke
    tiefe    = zzout(k)
    zzout(k) = sumtiefe + tiefe / 2._sp
    sumtiefe = sumtiefe + tiefe
    IF (verbose .GT. 0) PRINT *, NINT(sumtiefe), zzout(k)
  ENDDO

  CALL get(diag_out, iostat=ioerrstat)
  IF (ioerrstat .gt. 0) CALL err_exit
  CALL get(output_fname, iostat=ioerrstat)
  IF (ioerrstat .gt. 0) CALL err_exit
  CALL get(input_fname, iostat=ioerrstat)
  IF (ioerrstat .gt. 0) CALL err_exit
  IF (code == 'temp') CALL get(input1_fname, iostat=ioerrstat)
  IF (ioerrstat .gt. 0) CALL err_exit

  !   DO I=1, IE
  !     GAUSSLO(I)=-0.5+I
  !   ENDDO
  !   DO J=1, JE
  !     GAUSSLA(J)=90.5-j
  !   ENDDO

  DO I=1, IE
    DO J=1, JE
      GAULAT(I+1, J+1)=GAUSSLA(J)
      GAULON(I+1, J+1)=GAUSSLO(I)
    END DO
  END DO

  !     ZYKL.RAND OST/WEST
  DO J=1, JE+2
    GAULAT(1, J)=GAULAT(IE+1, J)
    GAULAT(IE+2, J)=GAULAT(2, J)
    GAULON(1, J)=GAULON(IE+1, J)-360._sp
    GAULON(IE+2, J)=GAULON(2, J)+360._sp
  END DO

  !     POLE (SCHUMMELN: OBEN UND UNTEN EINE REIHE HINZU)
  DO I=1, IE+2

    GAULAT(I, 1) = 180._sp - GAULAT(I, 2)
    GAULAT(I, JE+2) = -180._sp - GAULAT(I, JE+1)
    GAULON(I, 1) = GAULON(I, 2)
    GAULON(I, JE+2) = GAULON(I, JE+1)

  END DO

  !     NOCHMAL ZYKL.RAND OST/WEST
  DO J=2, JE+1
    GAULAT(1, J)=GAULAT(IE+1, J)
    GAULAT(IE+2, J)=GAULAT(2, J)
    GAULON(1, J)=GAULON(IE+1, J)-360._sp
    GAULON(IE+2, J)=GAULON(2, J)+360._sp
  END DO

  DO I=1, IE+2
    DO J=1, JE+2
      GAULAT(I, J) = GAULAT(I, J)/REAL(aradtogra, sp)
      GAULON(I, J) = GAULON(I, J)/REAL(aradtogra, sp)
    END DO
  END DO

  OPEN(io_in_anta, FILE='anta.ext4', ACCESS='SEQUENTIAL', &
       FORM='UNFORMATTED', ACTION='read')
  OPEN(io_out_diag, FILE=CHAR(diag_out), ACCESS='SEQUENTIAL', &
       FORM='UNFORMATTED')

  OPEN(io_out_forcing, FILE=CHAR(output_fname), ACCESS='SEQUENTIAL', &
       FORM='UNFORMATTED')

  OPEN(io_in_forcing, FILE=CHAR(input_fname), ACCESS='SEQUENTIAL', &
       FORM='UNFORMATTED', ACTION='read')

  IF (code == 'temp') THEN
    OPEN(21, FILE=CHAR(input1_fname), &
         ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
  END IF

  READ(io_in_anta) ext4_hdr
  READ(io_in_anta) GILA
  READ(io_in_anta) ext4_hdr
  READ(io_in_anta) GIPH

  !H    TAGESWERTE
  DO IDAY=1, LDAY
    IF (verbose .GT. 0) PRINT*, 'DATUM ', IDAY
    !H    READ INPUT FILES
    ext4_hdr(3)=1
    DO LEV=1, idepth_lvls
      IF3UP=MAX(1, ext4_hdr(3))

      READ(io_in_forcing) ext4_hdr
      READ(io_in_forcing)((DUMMY(I, J), I=1, IE), J=JE, 1, -1)
      IF (code == 'temp') THEN
        READ(21) ext4_hdr
        READ(21)((S(I, J), I=1, IE), J=JE, 1, -1)

        p=REAL(ext4_hdr(3), sp)/10._sp

        pt = MERGE(potemf(s, dummy, p), dummy, dummy .LT. 500._sp)
      ENDIF

      DO I=1, IE
        DO J=1, JE
          TEMIN(I+1, J+1) = DUMMY(I, J)
        END DO
      END DO

      !     DATENFELDER VORRBEITEN
      !     ZYKL.RAND OST/WEST

      DO J=1, JE+2
        TEMIN(1, J) = TEMIN(IE+1, J)
        TEMIN(IE+2, J) = TEMIN(2, J)
      ENDDO

      !     POLE (EINE REIHE HINZU)
      DO I=1, IE+2
        TEMIN(I, 1)=TEMIN(I, 2)
        TEMIN(I, JE+2)=TEMIN(I, JE+1)
      ENDDO

      IF (verbose .GT. 0) PRINT*, 'Datum ', iday, ' EINLESEN UEBERLEBT'

      !H    AN DIESER STELLE WIRD JEWEILS FUER U, V UND SKALARE GROESSEN
      !H    INTERPOLIERT, UND ZWAR AM JEWEILS DAZUGEHOERENDEN GITTERPUNKT
      !H    1. FUER SKALARE WIE GEHABT AN DRUCKPUNKT
      IF (verbose .GT. 0) PRINT*, 'SKALARE GROESSE INTERPOLIEREN'

      DO M=1, ME
        DO N=1, NE
          HOPLON(M, N)=GILA(2*M, 2*N)
          HOPLAT(M, N)=GIPH(2*M, 2*N)
          IF (hoplon(m, n) .GT. (2._sp * REAL(api, sp))) &
               hoplon(m, n) = hoplon(m, n) - (2._sp * REAL(api, sp))
          IF (hoplon(m, n).LT. 0._sp) &
               hoplon(m, n)=hoplon(m, n) + (2._sp * REAL(api, sp))
        ENDDO
      ENDDO


      DO I=1, IE+2
        DO J=1, JE+2
          land(i, j) = 0._sp
          IF (TEMIN(I, J) .eq. 999._sp) land(i, j) = 1._sp
          IF (TEMIN(I, J) .eq. 999._sp) TEMIN(I, j) =15._sp
          FIIN(I, J, 1)=TEMIN(I, J)
        ENDDO
      ENDDO
      CALL BLN2HOP(1, IE, JE, ME, NE, GAULAT, GAULON, HOPLAT, HOPLON, &
           FIIN, FIOUT, LAND)

      DO M=1, ME
        DO N=1, NE
          temout(m, n) = REAL(fiout(m, n, 1), wp)
        ENDDO
      ENDDO

      !H    AUSGABE

      IF (periodic_bound) THEN
        !h    PERIODISCHER RAND
        IF (verbose .GT. 1) PRINT*, 'PERIODISCHER RAND'
        DO N=1, NE
          TEMOUT(1, N)=TEMOUT(ME-1, N)
          TEMOUT(ME, N)=TEMOUT(2, N)
        END DO
      ENDIF

      output_hdr(3) = INT(ext4_hdr(3), i8)
      WRITE(io_out_diag) INT(IDAY, i8), output_hdr(2), &
           output_hdr(3), INT(ME*NE, i8)
      WRITE(io_out_diag) TEMOUT

      !HH   INTERPLOATE TO LEVEL
      DO M=1, ME
        DO N=1, NE
          tinter(m, n, lev) = temout(m, n)
        END DO
      END DO

      DO K=1, KE
        IF (NINT(ZZOUT(K)).GT.IF3UP &
             .AND.NINT(ZZOUT(K)).LE.ext4_hdr(3)) THEN
          IF (verbose .GT. 1) PRINT*, k, IF3UP, ZZOUT(K), ext4_hdr(3)
          DZL=REAL(ext4_hdr(3) - IF3UP, sp)
          DZS=REAL(ext4_hdr(3), sp) - ZZOUT(K)
          alpha = REAL(DZS/DZL, wp)
          DO M=1, ME
            DO N=1, NE
              TOUTER(M, N) = (alpha * tinter(m, n, lev-1)) &
                   + ((1._wp - alpha) * tinter(m, n, lev))
            END DO
          END DO
          IF (verbose .GT. 1) PRINT *, lev, alpha, tinter(14, 27, lev-1), &
               TOUTER(14, 27), TINTER(14, 27, lev)
          IF (verbose .GT. 1) PRINT *, 'SCHREIBE: ', k
          WRITE(io_out_forcing) INT(IDAY, i8), output_hdr(2), &
               INT(NINT(ZZOUT(K)), i8), INT(ME*NE, i8)
          WRITE(io_out_forcing) TOUTER
        ELSE
          IF (NINT(ZZOUT(K)).GT.ext4_hdr(3).AND.LEV .EQ. idepth_lvls)THEN
            ALPHA = 0._wp
            IF (verbose .GT. 1) PRINT *, K, IF3UP, ZZOUT(K), ext4_hdr(3)
            DO M=1, ME
              DO N=1, NE
                touter(m, n) = (alpha*tinter(m, n, lev-1)) &
                     + ((1._wp - alpha) * tinter(m, n, lev))
              ENDDO
            ENDDO
            IF (verbose .GT. 1) PRINT *, k, alpha, tinter(14, 27, lev-1), &
                 TOUTER(14, 27), TINTER(14, 27, lev)
            IF (verbose .GT. 1) PRINT*, 'SCHREIBE: ', k
            WRITE(io_out_forcing) INT(IDAY, i8), output_hdr(2), &
                 INT(NINT(ZZOUT(K)), i8), INT(me*ne, i8)
            WRITE(io_out_forcing) TOUTER
          END IF
        END IF
      END DO
      IF3UP=ext4_hdr(3)

    END DO

  END DO

  STOP

CONTAINS
  !> convert from phc data set in-situ temperatures to potential
  !> temperature (needed by mpiom)
  !> from Adrian Gill, Atmosphere-Ocean Dynamics (International
  !>        Geophysics), page 602
  !> @param s salinity
  !> @param t in-situ temperature (in degrees C)
  !> @param p pressure (in bars)
  !> @return potential temperature
  ELEMENTAL FUNCTION potemf(s, t, p)
    USE mo_kind, ONLY: sp
    REAL(sp) :: potemf
    REAL(sp), INTENT(in) :: s, t, p
    REAL(sp) :: a, b, c, d, e
    a = 3.6504e-4_sp + 8.3198e-5_sp * t &
         - 5.4065e-7_sp * t**2 &
         + 4.0274e-9_sp * t**3
    b = (s - 35._sp) * (1.7439e-5_sp - 2.9778e-7_sp * t)
    c = 8.9309e-7_sp  - 3.1628e-8_sp * t &
         + 2.1987e-10_sp * t**2
    d = 4.1057e-9_sp * (s - 35._sp)
    e = -1.6056e-10_sp + 5.0484e-12_sp * t
    potemf = t - a*p - b*p - c*p**2 + d*p**2 - e*p**3
  END FUNCTION potemf

  SUBROUTINE potem(s, t, p, pt)
    USE mo_kind
    IMPLICIT NONE
    INTEGER, PARAMETER :: IE=360, JE=180
    INTEGER :: i, j
    REAL(sp), INTENT(IN) :: s(ie, je), t(ie, je), p
    REAL(sp), INTENT(OUT) :: pt(ie, je)
    DO i=1, ie
      DO j=1, je
        pt(i, j) = potemf(s(i,j), t(i, j), p)
      END DO
    END DO
  END SUBROUTINE potem

END PROGRAM FORCING

!H    ********************************************************

SUBROUTINE BLN2HOP(NUMFI, IE, JE, ME, NE, ALAT2, ALON2, ALAT1, ALON1, &
     SOE, SREG, LAND)
  USE mo_kind
  USE mo_constants, ONLY: api!, aradtogra
  IMPLICIT NONE
  INTEGER, INTENT(in) :: numfi, ie, je, me, ne
  !     BILINEAR INTERPOLATION FROM ONE GRID(IE, JE) TO ANOTHER GRID(ME, NE)

  !      PARAMETER(NB=3, MB=3)

  REAL(sp), INTENT(in) :: ALAT2(IE+2, JE+2), ALON2(IE+2, JE+2), &
       ALAT1(ME, NE), ALON1(ME, NE), &
       LAND(IE+2, JE+2)
  REAL(sp), INTENT(inout) :: SOE(IE+2, JE+2, NUMFI), SREG(ME, NE, NUMFI)
  REAL(sp) :: HX(IE+2, JE+2, NUMFI), HHX(IE+2, JE+2, NUMFI), &
       ALPHA, BETA, HG(IE+2, JE+2), G(IE+2, JE+2), rsumg, &
       WWWALPHA, WWWBETA

  INTEGER :: i, j, l, il, ilu, ir, iter, jlu, jo, ju, m, n, nland
  !     DIFFUSION INTO LAND (NEW 12/99)

  DO L=1, NUMFI
    DO J=1, JE+2
      DO I=1, IE+2
        HHX(I, J, L)=SOE(I, J, L)
        !HH      X                 *((LAND(I, J)-1.))*(-1.)
      ENDDO
    ENDDO
  ENDDO

  DO J=1, JE+2
    DO I=1, IE+2
      HG(I, J) = 1._sp
      IF (LAND(I, J) .ge. 0.5_sp) HG(I, J) = 1.E-12_sp
    ENDDO
  ENDDO

  DO ITER=1, 6000
    DO J=1, JE+2
      DO I=1, IE+2
        G(I, J)=HG(I, J)
        DO L=1, NUMFI
          HX(I, J, l)=HHX(I, J, l)
        ENDDO
      ENDDO
    ENDDO
    DO J=1, JE+2
      JO=MAX(J-1, 1)
      JU=MIN(J+1, JE+2)
      DO I=1, IE+2
        IL=I-1
        IF(IL.LT.1)IL=IE
        IR=I+1
        IF(IR.GT.IE+2)IR=3
        RSUMG = 0._sp
        IF(LAND(I, J) .GE. 0.5_sp)THEN
          RSUMG=(4._sp * G(I, J) &
               + G(IL, J) + G(IR, J) + G(I, JO) + G(I, JU))/8._sp

          HG(I, J)=MIN(RSUMG, 0.125_sp)

          DO L=1, NUMFI
            HHX(I, J, L)=(4._sp*HX(I, J, L)*G(I, J) &
                 +HX(IL, J, L)*G(IL, J) &
                 +HX(IR, J, L)*G(IR, J) &
                 +HX(I, JO, L)*G(I, JO) &
                 +HX(I, JU, L)*G(I, JU))/8._sp

            HHX(I, J, L)=HHX(I, J, L)/RSUMG
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    nland=0
    Do I=1, IE
      DO j=1, je
        IF (HG(i, j) .LE. 2.E-12_sp) nland = nland + 1
      ENDDO
    ENDDO
    !          PRINT*, iter, ' Nland: ', nland
  ENDDO


  DO J=1, JE+2
    DO I=1, IE+2
      G(i, j)=HG(i, j)
      DO L=1, NUMFI
        SOE(I, J, L)=HHX(I, J, L)
        If (ABS(SOE(I, J, L)) .LT. 1.E-12_sp) SOE(I, J, L) = 0._sp
      ENDDO
    ENDDO
  ENDDO

  !
  DO M=1, ME
    DO N=1, NE

      !        PUNKT RECHTS OBEN
      DO J=1, JE+2
        IF (ALAT2(2, J).GE.ALAT1(M, N)) JO=J
      ENDDO
      DO I=1, IE+2
        IF (ALON2(I, 2).LE.ALON1(M, N)) IL=I
      ENDDO

      !       IL=NINT(ALON1(M, N)*DLAMDAI+0.499999999)+1

      !            PRINT*, 'WE ARE AT THE MINIMUM ', M, N, &
      !                     ALAT1(M, N)*REAL(aradtogra, sp)
      !     X, IL, JO, ALAT2(IL, JO)*REAL(aradtogra, sp)

      !       PUNKT RECHTS OBEN --> LINKS UNTEN

      JLU=JO+1
      ILU=IL

      WWWALPHA=ALAT2(ILU, JLU)-ALAT1(M, N)
      WWWBETA=ALON2(ILU, JLU)-ALON1(M, N)
      IF (wwwbeta .GE. REAL(api, sp)) wwwbeta = wwwbeta - 2._sp * REAL(api, sp)
      IF (wwwbeta .LE. REAL(-api, sp)) wwwbeta = wwwbeta + 2._sp * REAL(api, sp)

      ALPHA=(WWWALPHA)/(ALAT2(ILU, JLU)-ALAT2(ILU, JLU-1))
      BETA=(WWWBETA)/(ALON2(ILU, JLU)-ALON2(ILU+1, JLU))
      DO I=1, NUMFI
        SREG(M, N, I)=ALPHA*BETA*SOE(ILU+1, JLU-1, I)*G(ILU+1, JLU-1)   &
             + (1._sp - ALPHA)*(1._sp - BETA)*SOE(ILU, JLU, I)*G(ILU, JLU) &
             + (1._sp - ALPHA)*(BETA)*SOE(ILU+1, JLU, I)*G(ILU+1, JLU) &
             +(ALPHA)*(1._sp - BETA)*SOE(ILU, JLU-1, I)*G(ILU, JLU-1)
        SREG(M, N, I)=SREG(M, N, I)/ &
             (ALPHA * BETA * G(ILU+1, JLU-1) &
             + (1._sp - ALPHA)*(1._sp - BETA)*G(ILU, JLU) &
             + (1._sp - ALPHA)*(BETA) * G(ILU+1, JLU) &
             + (ALPHA) * (1._sp - BETA) * G(ILU, JLU-1))
      ENDDO
      !         *************************************************
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE BLN2HOP

SUBROUTINE err_exit
  PRINT *, 'io error occurred!, aborting'
  STOP 1
END SUBROUTINE err_exit

FUNCTION gausslo(i)
  USE mo_kind
  INTEGER, INTENT(in) :: i
  REAL(sp) gausslo
  gausslo = -0.5_sp + REAL(i, sp)
END FUNCTION gausslo

FUNCTION gaussla(j)
  USE mo_kind
  INTEGER, INTENT(in) :: j
  REAL(sp) gaussla
  GAUSSLA = 90.5_sp - REAL(j, sp)
END FUNCTION gaussla
