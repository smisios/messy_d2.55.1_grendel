MODULE trigonometry
  USE mo_kind, ONLY: i4, sp
  USE mo_planetary_constants, ONLY: rhoref_water
  IMPLICIT NONE
  REAL(sp), PARAMETER :: PI = 3.14159265358979323846_sp, &
       PIDEG = 1.7453292519943E-2_sp, &
       DEGPI = 57.29577951308_sp, &
       GRARAD = 180._sp / PI
END MODULE trigonometry

PROGRAM FORCING
  USE mo_kind, ONLY: i4, sp
  USE trigonometry

  IMPLICIT NONE

  INTEGER(i4), PARAMETER :: IE = 192, JE = 94
  REAL(sp),    PARAMETER :: ZSMALL=1.e-20
  INTEGER(i4)            :: ME, NE, iyear, lday, i, iday, imon, j, m, n, ml, mr
  REAL(sp), ALLOCATABLE  :: GILA(:,:),GIPH(:,:)
  REAL(sp), ALLOCATABLE  :: GIXYZ(:,:,:)
  INTEGER(i4)            :: dummy_i4(4)
  REAL(sp)               :: GAULAT(IE+2,JE+2),GAULON(IE+2,JE+2)
  REAL(sp)               :: LAND(IE+2,JE+2),DUMMY(IE,JE)
  REAL(sp)               :: GAUSSLA(JE),GAUSSLO(IE)

  LOGICAL                :: lbounds_exch_tp

  NAMELIST /forcingparams/ me, ne, lday, iyear, lbounds_exch_tp

  REAL(sp)               :: FIIN(IE+2,JE+2,8)
  REAL(sp), ALLOCATABLE  :: FIOUT(:,:,:)

  REAL(sp) :: UWIN(IE+2,JE+2), &
       VWIN(IE+2,JE+2), &
       TEMIN(IE+2,JE+2), &
       TDEWIN(IE+2,JE+2), &
       SWRADIN(IE+2,JE+2), &
       SLRADIN(IE+2,JE+2), &
       TCDCIN(IE+2,JE+2), &
       PRATEIN(IE+2,JE+2), &
       PRESSIN(IE+2,JE+2), &
       U10IN(IE+2,JE+2), &
       UIN(IE+2,JE+2), &
       VIN(IE+2,JE+2)

  REAL(sp), ALLOCATABLE :: UWOUT(:,:), VWOUT(:,:), TEMOUT(:,:), &
       TDEWOUT(:,:), SWRADOUT(:,:), SLRADOUT(:,:), TCDCOUT(:,:), &
       PRATEOUT(:,:), PRESSOUT(:,:), U10OUT(:,:), UOUT(:,:), &
       VOUT(:,:), UKE(:,:),VKE(:,:),UHOPE(:,:),VHOPE(:,:), &
       UE(:,:),VE(:,:),UH(:,:),VH(:,:), HOPLAT(:,:),HOPLON(:,:)


  REAL(sp) :: SVAL, R_i(3), R_j(3), R_tmp(3), R_poli(2), R_polj(2), R_pol(2), &
       R_iabs, R_jabs, schugeo, schugr, pq

  INTEGER(i4) :: L, IEXT(4)

  DATA SVAL / 0./

  !     NCEP GITTER
  DATA GAUSSLA/88.542, 86.6531, 84.7532, 82.8508, 80.9473, 79.0435, &
       77.1394, 75.2351, 73.3307, 71.4262, 69.5217, 67.6171, &
       65.7125, 63.8079, 61.9033, 59.9986, 58.0939, 56.1893, &
       54.2846, 52.3799, 50.4752, 48.5705, 46.6658, 44.7611, &
       42.8564, 40.9517, 39.047, 37.1422, 35.2375, 33.3328, &
       31.4281, 29.5234, 27.6186, 25.7139, 23.8092, 21.9044, &
       19.9997, 18.095, 16.1902, 14.2855, 12.3808, 10.47604, &
       8.57131, 6.66657, 4.76184, 2.8571, 0.952368, -0.952368, &
       -2.8571, -4.76184, -6.66657, -8.57131, -10.47604, -12.3808, &
       -14.2855, -16.1902, -18.095, -19.9997, -21.9044, -23.8092, &
       -25.7139, -27.6186, -29.5234, -31.4281, -33.3328, -35.2375, &
       -37.1422, -39.047, -40.9517, -42.8564, -44.7611, -46.6658, &
       -48.5705, -50.4752, -52.3799, -54.2846, -56.1893, -58.0939, &
       -59.9986, -61.9033, -63.8079, -65.7125, -67.6171, -69.5217, &
       -71.4262, -73.3307, -75.2351, -77.1394, -79.0435, -80.9473, &
       -82.8508, -84.7532, -86.6531, -88.542 /

  DATA GAUSSLO/0.,1.875,3.75,5.625, 7.5, 9.375, 11.25, 13.125, 15., &
       16.875, 18.75, 20.625, 22.5, 24.375, 26.25, 28.125, 30., 31.875, &
       33.75, 35.625, 37.5, 39.375, 41.25, 43.125, 45., 46.875, 48.75, &
       50.625, 52.5, 54.375, 56.25, 58.125, 60., 61.875, 63.75, 65.625, &
       67.5, 69.375, 71.25, 73.125, 75., 76.875, 78.75, 80.625, 82.5, &
       84.375, 86.25, 88.125, 90., 91.875, 93.75, 95.625, 97.5, 99.375, &
       101.25, 103.125, 105., 106.875, 108.75, 110.625, 112.5, 114.375, &
       116.25, 118.125, 120., 121.875, 123.75, 125.625, 127.5, 129.375, &
       131.25, 133.125, 135., 136.875, 138.75, 140.625, 142.5, 144.375, &
       146.25, 148.125, 150., 151.875, 153.75, 155.625, 157.5, 159.375, &
       161.25, 163.125, 165., 166.875, 168.75, 170.625, 172.5, 174.375, &
       176.25, 178.125, 180., 181.875, 183.75, 185.625, 187.5, 189.375, &
       191.25, 193.125, 195., 196.875, 198.75, 200.625, 202.5, 204.375, &
       206.25, 208.125, 210., 211.875, 213.75, 215.625, 217.5, 219.375, &
       221.25, 223.125, 225., 226.875, 228.75, 230.625, 232.5, 234.375, &
       236.25, 238.125, 240., 241.875, 243.75, 245.625, 247.5, 249.375, &
       251.25, 253.125, 255., 256.875, 258.75, 260.625, 262.5, 264.375, &
       266.25, 268.125, 270., 271.875, 273.75, 275.625, 277.5, 279.375, &
       281.25, 283.125, 285., 286.875, 288.75, 290.625, 292.5, 294.375, &
       296.25, 298.125, 300., 301.875, 303.75, 305.625, 307.5, 309.375, &
       311.25, 313.125, 315., 316.875, 318.75, 320.625, 322.5, 324.375, &
       326.25, 328.125, 330., 331.875, 333.75, 335.625, 337.5, 339.375, &
       341.25, 343.125, 345., 346.875, 348.75, 350.625, 352.5, 354.375, &
       356.25, 358.125 /
  me = -1
  ne = -1
  lday = -1
  READ (*, nml=forcingparams)

  IF (me < 1 .OR. ne < 1 .OR. lday < 1) THEN
    PRINT *, 'Only positive integral values may be used for me, ne and lday!'
    STOP 1
  END IF
  ALLOCATE(GILA(2*ME,2*NE), GIPH(2*ME,2*NE))
  allocate(GIXYZ(3,2*ME,2*NE))
  ALLOCATE(UWOUT(ME,NE), VWOUT(ME,NE), TEMOUT(ME,NE), &
       TDEWOUT(ME,NE), SWRADOUT(ME,NE), SLRADOUT(ME,NE), TCDCOUT(ME,NE), &
       PRATEOUT(ME,NE), PRESSOUT(ME,NE), U10OUT(ME,NE), UOUT(ME,NE), &
       VOUT(ME,NE), UKE(ME,NE),VKE(ME,NE),UHOPE(ME,NE),VHOPE(ME,NE), &
       UE(ME,NE), VE(ME,NE), UH(ME,NE), VH(ME,NE), HOPLAT(ME,NE), HOPLON(ME,NE))
  ALLOCATE(FIOUT(ME,NE,8))

  !     SOME CONSTANTS

  DO I=1,IE
    DO J=1,JE
      GAULAT(I+1,J+1)=GAUSSLA(J)
      GAULON(I+1,J+1)=GAUSSLO(I)
    ENDDO
  ENDDO

  !     ZYKL.RAND OST/WEST
  DO J=1,JE+2
    GAULAT(1,J)=GAULAT(IE+1,J)
    GAULAT(IE+2,J)=GAULAT(2,J)
    GAULON(1,J)=GAULON(IE+1,J)-360.
    GAULON(IE+2,J)=GAULON(2,J)+360.
  ENDDO

  !     POLE (SCHUMMELN: OBEN UND UNTEN EINE REIHE HINZU)
  DO I=1,IE+2

    GAULAT(I,1)=180.-GAULAT(I,2)
    GAULAT(I,JE+2)=-180.-GAULAT(I,JE+1)
    GAULON(I,1)=GAULON(I,2)
    GAULON(I,JE+2)=GAULON(I,JE+1)

  ENDDO

  !     NOCHMAL ZYKL.RAND OST/WEST
  DO J=2,JE+1
    GAULAT(1,J)=GAULAT(IE+1,J)
    GAULAT(IE+2,J)=GAULAT(2,J)
    GAULON(1,J)=GAULON(IE+1,J)-360.
    GAULON(IE+2,J)=GAULON(2,J)+360.
  ENDDO

  DO I=1,IE+2
    DO J=1,JE+2
      GAULAT(I,J)=GAULAT(I,J)/GRARAD
      GAULON(I,J)=GAULON(I,J)/GRARAD
    ENDDO
  ENDDO

  OPEN(81,FILE='anta.ext',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(10,FILE='land.ext',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

  READ(81)(dummy_i4(i),i=1,4)
  READ(81) GILA
  READ(81)(dummy_i4(i),i=1,4)
  READ(81) GIPH
  do i = 1,2*ME
    do j = 1,2*NE
      CALL POLXYZ(1._sp,REAL(GILA(i,j),sp),REAL(GIPH(i,j),sp),R_tmp)
      GIXYZ(:,i,j)=REAL(R_tmp(:),sp)
    enddo
  enddo

  !H    TAGESWERTE
  !      OPEN(1,FILE='FYEAR')
  !      READ(1,*)IYEAR
  !      CLOSE(1)

  READ(10)(dummy_i4(i),i=1,4)
  READ(10)((DUMMY(I,J),I=1,IE),J=1,JE)
  DO I=1,IE
    DO J=1,JE
      LAND(I+1,J+1)=DUMMY(I,J)
      IF (land(I+1,J+1) .GT. 1.e-4) land(I+1,J+1)=1.
    ENDDO
  ENDDO
  !     ZYKL.RAND OST/WEST/NORDPOL/SUEDPOL
  DO I=1,IE+2
    LAND(I,1)=0.
    LAND(I,JE+2)=1.
  ENDDO
  DO J=1,JE+2
    LAND(1,J)=LAND(IE+1,J)
    LAND(IE+2,J)=LAND(2,J)
  ENDDO

  OPEN(2,FILE='air.2m' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(11,FILE='tdew.2m' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(12,FILE='prate.sfc' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(13,FILE='tcdc.eatm' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(14,FILE='uvwnd.10m' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(15,FILE='dswrf.sfc' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(16,FILE='dlwrf.sfc' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(17,FILE='pres.sfc' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(3,FILE='uflx.sfc' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(4,FILE='vflx.sfc' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(18,FILE='uwnd.10m' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(19,FILE='vwnd.10m' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

  OPEN(52,FILE='GIWIX' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(53,FILE='GIWIY' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(54,FILE='GITEM' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(55,FILE='GITDEW' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(56,FILE='GIPREC' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(57,FILE='GICLOUD' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(58,FILE='GIU10' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(59,FILE='GISWRAD' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(60,FILE='GILWRAD' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(61,FILE='GIPRESS' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(62,FILE='GIWX' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(63,FILE='GIWY' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

  DO L=1,LDAY

    !H    READ INPUT FILES

    READ(2)IEXT
    READ(2) DUMMY !((DUMMY(I,J),J=1,JE),I=1,IE)  ! air2m
    DO I=1,IE
      DO J=1,JE
        TEMIN(I+1,J+1)=DUMMY(I,J)
      END DO
    END DO

    READ(11)IEXT
    READ(11) DUMMY !((DUMMY(I,J),J=1,JE),I=1,IE) ! tdew
    DO I=1,IE
      DO J=1,JE
        TDEWIN(I+1,J+1)=DUMMY(I,J)
      END DO
    END DO

    READ(12)IEXT
    READ(12) DUMMY !((DUMMY(I,J),J=1,JE),I=1,IE) ! prate
    DO I=1,IE
      DO J=1,JE
        PRATEIN(I+1,J+1)=DUMMY(I,J)
      END DO
    END DO

    READ(17)IEXT
    READ(17) DUMMY !((DUMMY(I,J),J=1,JE),I=1,IE) ! pres
    DO I=1,IE
      DO J=1,JE
        PRESSIN(I+1,J+1)=DUMMY(I,J)
      END DO
    END DO

    READ(13)IEXT
    READ(13) DUMMY !((DUMMY(I,J),I=1,IE),J=1,JE) ! tcdc
    DO I=1,IE
      DO J=1,JE
        TCDCIN(I+1,J+1)=DUMMY(I,J)
      END DO
    END DO

    READ(14)IEXT
    READ(14) DUMMY !((DUMMY(I,J),J=1,JE),I=1,IE) ! uvwn
    DO I=1,IE
      DO J=1,JE
        U10IN(I+1,J+1)=DUMMY(I,J)
      END DO
    END DO

    READ(15)IEXT
    READ(15) DUMMY !((DUMMY(I,J),J=1,JE),I=1,IE) ! dswr
    DO I=1,IE
      DO J=1,JE
        SWRADIN(I+1,J+1) = MAX(DUMMY(I,J), 0._sp)
      END DO
    END DO

    READ(16)IEXT
    READ(16) DUMMY !((DUMMY(I,J),J=1,JE),I=1,IE) ! dlwr
    DO I=1,IE
      DO J=1,JE
        SLRADIN(I+1,J+1) = MAX(DUMMY(I,J), 0._sp)
      ENDDO
    ENDDO

    READ(3)IEXT
    READ(3) DUMMY !((DUMMY(I,J),J=1,JE),I=1,IE) ! uflx
    DO I=1,IE
      DO J=1,JE
        UWIN(I+1,J+1)= DUMMY(I,J)*(-1.)
      END DO
    END DO

    READ(4)IEXT
    READ(4) DUMMY !((DUMMY(I,J),J=1,JE),I=1,IE) ! vflx
    DO I=1,IE
      DO J=1,JE
        VWIN(I+1,J+1)= DUMMY(I,J)*(-1.)
      END DO
    END DO

    READ(18)IEXT
    READ(18) DUMMY !((DUMMY(I,J),J=1,JE),I=1,IE) ! u
    DO I=1,IE
      DO J=1,JE
        UIN(I+1,J+1)= DUMMY(I,J)*(-1.)
      ENDDO
    ENDDO

    READ(19)IEXT
    READ(19) DUMMY !((DUMMY(I,J),J=1,JE),I=1,IE) ! v
    DO I=1,IE
      DO J=1,JE
        VIN(I+1,J+1)= DUMMY(I,J)*(-1.)
      ENDDO
    ENDDO

    !     DATENFELDER VORRBEITEN
    !     ZYKL.RAND OST/WEST

    DO J=1,JE+2
      UWIN(1,J)=UWIN(IE+1,J)
      UWIN(IE+2,J)=UWIN(2,J)

      VWIN(1,J)=VWIN(IE+1,J)
      VWIN(IE+2,J)=VWIN(2,J)

      UIN(1,J)=UIN(IE+1,J)
      UIN(IE+2,J)=UIN(2,J)

      VIN(1,J)=VIN(IE+1,J)
      VIN(IE+2,J)=VIN(2,J)

      TEMIN(1,J)=TEMIN(IE+1,J)
      TEMIN(IE+2,J)=TEMIN(2,J)

      TDEWIN(1,J)=TDEWIN(IE+1,J)
      TDEWIN(IE+2,J)=TDEWIN(2,J)

      PRATEIN(1,J)=PRATEIN(IE+1,J)
      PRATEIN(IE+2,J)=PRATEIN(2,J)

      TCDCIN(1,J)=TCDCIN(IE+1,J)
      TCDCIN(IE+2,J)=TCDCIN(2,J)

      U10IN(1,J)=U10IN(IE+1,J)
      U10IN(IE+2,J)=U10IN(2,J)

      SWRADIN(1,J)=SWRADIN(IE+1,J)
      SWRADIN(IE+2,J)=SWRADIN(2,J)

      SLRADIN(1,J)=SLRADIN(IE+1,J)
      SLRADIN(IE+2,J)=SLRADIN(2,J)

      PRESSIN(1,J)=PRESSIN(IE+1,J)
      PRESSIN(IE+2,J)=PRESSIN(2,J)

    ENDDO

    !     POLE (EINE REIHE HINZU)
    DO I=1,IE+2

      UWIN(I,1)=UWIN(I,2)
      UWIN(I,JE+2)=UWIN(I,JE+1)

      VWIN(I,1)=VWIN(I,2)
      VWIN(I,JE+2)=VWIN(I,JE+1)

      UIN(I,1)=UIN(I,2)
      UIN(I,JE+2)=UIN(I,JE+1)

      VIN(I,1)=VIN(I,2)
      VIN(I,JE+2)=VIN(I,JE+1)

      TEMIN(I,1)=TEMIN(I,2)
      TEMIN(I,JE+2)=TEMIN(I,JE+1)

      TDEWIN(I,1)=TDEWIN(I,2)
      TDEWIN(I,JE+2)=TDEWIN(I,JE+1)

      PRATEIN(I,1)=PRATEIN(I,2)
      PRATEIN(I,JE+2)=PRATEIN(I,JE+1)

      TCDCIN(I,1)=TCDCIN(I,2)
      TCDCIN(I,JE+2)=TCDCIN(I,JE+1)

      U10IN(I,1)=U10IN(I,2)
      U10IN(I,JE+2)=U10IN(I,JE+1)

      SWRADIN(I,1)=SWRADIN(I,2)
      SWRADIN(I,JE+2)=SWRADIN(I,JE+1)

      SLRADIN(I,1)=SLRADIN(I,2)
      SLRADIN(I,JE+2)=SLRADIN(I,JE+1)

      PRESSIN(I,1)=PRESSIN(I,2)
      PRESSIN(I,JE+2)=PRESSIN(I,JE+1)

    ENDDO

    CALL day_of_year2day_of_month(l, lday, imon, iday)
    WRITE(*,'(a,i4.4,''-'',i2.2,''-'',i2.2,2x,a,i8.8,2x,a)') &
         'Date: ', iyear, imon, iday, 'from file: ', iext(1), 'all fields read'

    !H    AN DIESER STELLE WIRD JEWEILS FUER U,V UND SKALARE GROESSEN
    !H    INTERPOLIERT, UND ZWAR AM JEWEILS DAZUGEHOERENDEN GITTERPUNKT
    !H    1. FUER SKALARE WIE GEHABT AN DRUCKPUNKT
    !      PRINT*,'SKALARE GROESSE INTERPOLIEREN'

    DO M=1,ME
      DO N=1,NE
        HOPLON(M,N)=REAL(GILA(2*M,2*N),sp)
        HOPLAT(M,N)=REAL(GIPH(2*M,2*N),sp)
        IF (HOPLON(M,N) .GT. (2.*PI)) HOPLON(M,N)=HOPLON(M,N)-(2.*PI)
        IF (HOPLON(M,N) .LT. 0.) HOPLON(M,N)=HOPLON(M,N)+(2.*PI)
      ENDDO
    ENDDO

    !      CPUANF=SECOND()
    DO I=1,IE+2
      DO J=1,JE+2
        FIIN(I,J,1)=TEMIN(I,J)
        FIIN(I,J,2)=TDEWIN(I,J)-273.15
        FIIN(I,J,3)=PRATEIN(I,J)
        FIIN(I,J,4)=TCDCIN(I,J)
        FIIN(I,J,5)=U10IN(I,J)
        FIIN(I,J,6)=SWRADIN(I,J)
        FIIN(I,J,7)=SLRADIN(I,J)
        FIIN(I,J,8)=PRESSIN(I,J)
      ENDDO
    ENDDO
    CALL BLN2HOP(8,IE,JE,ME,NE,GAULAT,GAULON,HOPLAT,HOPLON, &
         FIIN,FIOUT,LAND)

    DO M=1,ME
      DO N=1,NE
        TEMOUT(M,N)=FIOUT(M,N,1)
        TDEWOUT(M,N)=FIOUT(M,N,2)+273.15
        PRATEOUT(M,N)=FIOUT(M,N,3)
        TCDCOUT(M,N)=FIOUT(M,N,4)
        U10OUT(M,N)=FIOUT(M,N,5)
        SWRADOUT(M,N)=FIOUT(M,N,6)
        SLRADOUT(M,N)=FIOUT(M,N,7)
        PRESSOUT(M,N)=FIOUT(M,N,8)
      ENDDO
    ENDDO
    !      CPUEND=SECOND()
    ! 1    PRINT*,'CALL BLN2HOP DAUERT ',CPUEND-CPUANF

    !H    AUSGABE
    !h    PERIODISCHER RAND
    DO N=1,NE

      TEMOUT(1,N)=TEMOUT(ME-1,N)
      TEMOUT(ME,N)=TEMOUT(2,N)

      TDEWOUT(1,N)=TDEWOUT(ME-1,N)
      TDEWOUT(ME,N)=TDEWOUT(2,N)

      PRATEOUT(1,N)=PRATEOUT(ME-1,N)
      PRATEOUT(ME,N)=PRATEOUT(2,N)

      TCDCOUT(1,N)=TCDCOUT(ME-1,N)
      TCDCOUT(ME,N)=TCDCOUT(2,N)

      U10OUT(1,N)=U10OUT(ME-1,N)
      U10OUT(ME,N)=U10OUT(2,N)

      SWRADOUT(1,N)=SWRADOUT(ME-1,N)
      SWRADOUT(ME,N)=SWRADOUT(2,N)

      SLRADOUT(1,N)=SLRADOUT(ME-1,N)
      SLRADOUT(ME,N)=SLRADOUT(2,N)

      PRESSOUT(1,N)=PRESSOUT(ME-1,N)
      PRESSOUT(ME,N)=PRESSOUT(2,N)

    ENDDO

    IF (lbounds_exch_tp) THEN
      DO m=2,me-1
        ml=m
        mr=me+1-m
        TEMOUT(ml,2)=TEMOUT(mr,3)         ! synchronise line 2 with line 3
        TEMOUT(ml,1)=TEMOUT(mr,4)         ! synchronise line line 1 with line 4

        TDEWOUT(ml,2)=TDEWOUT(mr,3)
        TDEWOUT(ml,1)=TDEWOUT(mr,4)

        PRATEOUT(ml,2)=PRATEOUT(mr,3)
        PRATEOUT(ml,1)=PRATEOUT(mr,4)

        TCDCOUT(ml,2)=TCDCOUT(mr,3)
        TCDCOUT(ml,1)=TCDCOUT(mr,4)

        U10OUT(ml,2)=U10OUT(mr,3)
        U10OUT(ml,1)=U10OUT(mr,4)

        SWRADOUT(ml,2)=SWRADOUT(mr,3)
        SWRADOUT(ml,1)=SWRADOUT(mr,4)

        SLRADOUT(ml,2)=SLRADOUT(mr,3)
        SLRADOUT(ml,1)=SLRADOUT(mr,4)

        PRESSOUT(ml,2)=PRESSOUT(mr,3)
        PRESSOUT(ml,1)=PRESSOUT(mr,4)
      ENDDO
    ENDIF

    !!!IEXT(1) = IYEAR*10000 + IMON*100 + IDAY
    IEXT(3) = 0
    IEXT(4) = ME*NE
    IEXT(2) = 92
    WRITE(54) IEXT
    WRITE(54) TEMOUT
    IEXT(2)=81
    WRITE(55) IEXT
    WRITE(55) TDEWOUT
    IEXT(2)=260
    WRITE(56) IEXT
    WRITE(56) PRATEOUT
    IEXT(2)=164
    WRITE(57) IEXT
    WRITE(57) TCDCOUT
    IEXT(2)=171
    WRITE(58) IEXT
    WRITE(58) U10OUT
    IEXT(2)=80
    WRITE(59) IEXT
    WRITE(59) SWRADOUT
    IEXT(2)=177
    WRITE(60) IEXT
    WRITE(60) SLRADOUT
    IEXT(2)=151
    WRITE(61) IEXT
    WRITE(61) PRESSOUT

    !H    2. FUER VECTOR U KOMPONENTE AM U-PUNKT
    !      PRINT*,'U-KOMPOMENTE INTERPOLIEREN'

    DO M=2,ME-1
      DO N=1,NE
        HOPLON(M,N)=REAL(GILA((2*M)+1,2*N),sp)
        HOPLAT(M,N)=REAL(GIPH((2*M)+1,2*N),sp)
        IF (HOPLON(M,N) .GT. (2.*PI)) HOPLON(M,N)=HOPLON(M,N)-(2.*PI)
        IF (HOPLON(M,N) .LT. 0.) HOPLON(M,N)=HOPLON(M,N)+(2.*PI)
      ENDDO
    ENDDO

    DO I=1,IE+2
      DO J=1,JE+2
        FIIN(I,J,1)=UWIN(I,J)
        FIIN(I,J,2)=VWIN(I,J)
        FIIN(I,J,3)=UIN(I,J)
        FIIN(I,J,4)=VIN(I,J)
      ENDDO
    ENDDO

    CALL BLN2HOP(4,IE,JE,ME,NE,GAULAT,GAULON,HOPLAT,HOPLON, &
         FIIN,FIOUT,LAND)

    DO M=2,ME-1
      DO N=2,NE-1
        UWOUT(M,N)=FIOUT(M,N,1)
        VWOUT(M,N)=FIOUT(M,N,2)
        UOUT(M,N)=FIOUT(M,N,3)
        VOUT(M,N)=FIOUT(M,N,4)
      ENDDO
    ENDDO

    !H     DREHEN DER LOKALEN VEKTOREN IN X UND Y RICHTUNG
    !H     DIE VEKTOREN (DELTAXX,DELTAXY) UND (DELTAYY,DELTAYX)
    !H     SPANNEN DAS N/S, O/W SYSTEM AUF
    !H     WINDSTRESS (N/kg)/rho

    !       PRINT*,'U-DREHEN'

    DO M=2,ME-1
      DO N=2,NE-1

        IF ((UWOUT(M,N) .NE. SVAL) .AND. (VWOUT(M,N) .NE. SVAL)) THEN

          R_pol(1)=HOPLON(M,N)
          R_pol(2)=HOPLAT(M,N)

          UWOUT(M,N)=UWOUT(M,N)/REAL(rhoref_water, sp)
          VWOUT(M,N)=VWOUT(M,N)/REAL(rhoref_water, sp)

          R_i(:)=REAL(GIXYZ(:,2*M+1,2*N)-GIXYZ(:,2*M,2*N),sp)
          R_iabs=SQRT(SUM(R_i(:)**2))
          R_i(:)=R_i(:)/R_iabs
          CALL dek2geo(R_i,R_tmp,R_pol)
          R_poli=R_tmp(1:2)
          R_j(:)=REAL(GIXYZ(:,2*M+1,2*N-1)-GIXYZ(:,2*M+1,2*N),sp)
          R_jabs=SQRT(SUM(R_j(:)**2))
          R_j=R_j/R_jabs
          IF ( (R_jabs .LE. ZSMALL) .OR. (R_iabs .LE. ZSMALL) ) THEN
            WRITE(*,*) M,N, R_i, R_j
            STOP
          END IF
          CALL dek2geo(R_j,R_tmp,R_pol)
          R_polj=R_tmp(1:2)
          !         write(*,*) R_tmp

          UKE(M,N) = (UWOUT(M,N)*R_poli(1)) + (VWOUT(M,N)*R_poli(2))
          VKE(M,N) = (UWOUT(M,N)*R_polj(1)) + (VWOUT(M,N)*R_polj(2))
          UE(M,N) = (UOUT(M,N)*R_poli(1)) + (VOUT(M,N)*R_poli(2))
          VE(M,N) = (UOUT(M,N)*R_polj(1)) + (VOUT(M,N)*R_polj(2))

          !H       ES WIRD PASSEND SKALIERT

          !H       BETRAG VOR DER DREHUNG
          SCHUGEO = SQRT(UWOUT(M,N)**2 + VWOUT(M,N)**2)

          !H       BETRAG NACH DER DREHUNG
          SCHUGR = SQRT(UKE(M,N)**2+VKE(M,N)**2)

          !H       QUOTIENT ALS KORREKTURFAKTOR

          UHOPE(M,N) = UKE(M,N)*SCHUGEO / SCHUGR
          !
          SCHUGEO = SQRT(UOUT(M,N)**2 + VOUT(M,N)**2)

          !H       BETRAG NACH DER DREHUNG
          SCHUGR = SQRT(UE(M,N)**2+VE(M,N)**2)

          UH(M,N) = UE(M,N)*SCHUGEO / SCHUGR
        ELSE

          UHOPE(M,N)=0.
          UH(M,N)=0.
        ENDIF

      ENDDO
    ENDDO

    !H    AUSGABE
    !h    PERIODISCHER RAND
    DO N=1,NE
      UHOPE(1,N)=UHOPE(ME-1,N)
      UHOPE(ME,N)=UHOPE(2,N)
      UH(1,N)=UH(ME-1,N)
      UH(ME,N)=UH(2,N)
    END DO

    IF (lbounds_exch_tp) THEN
      DO m=2,me-1
        ml=m
        mr=me-m
        UHOPE(ml,2) = -UHOPE(mr,3)    ! syncronise line 2 with line 3
        UHOPE(ml,1) = -UHOPE(mr,4)    ! syncronise line 1 with line 4
      END DO
    ENDIF

    IEXT(2)=52
    WRITE(52) IEXT
    WRITE(52) UHOPE

    IEXT(2)=165
    WRITE(62) IEXT
    WRITE(62) UH


    !H    3. FUER VECTOR V KOMPONENTE AM V-PUNKT
    !      PRINT*,'V-KOMPOMENTE INTERPOLIEREN'

    DO M=1,ME
      DO N=2,NE-1
        HOPLON(M,N)=REAL(GILA(2*M,(2*N)+1),sp)
        HOPLAT(M,N)=REAL(GIPH(2*M,(2*N)+1),sp)
        IF (HOPLON(M,N) .GT. (2.*PI)) HOPLON(M,N)=HOPLON(M,N)-(2.*PI)
        IF (HOPLON(M,N) .LT. 0.) HOPLON(M,N)=HOPLON(M,N)+(2.*PI)
      END DO
    END DO

    DO I=1,IE+2
      DO J=1,JE+2
        FIIN(I,J,1)=UWIN(I,J)
        FIIN(I,J,2)=VWIN(I,J)
        FIIN(I,J,3)=UIN(I,J)
        FIIN(I,J,4)=VIN(I,J)
      END DO
    END DO

    CALL BLN2HOP(4,IE,JE,ME,NE,GAULAT,GAULON,HOPLAT,HOPLON, &
         FIIN,FIOUT,LAND)

    DO M=1,ME
      DO N=1,NE
        UWOUT(M,N)=FIOUT(M,N,1)
        VWOUT(M,N)=FIOUT(M,N,2)
        UOUT(M,N)=FIOUT(M,N,3)
        VOUT(M,N)=FIOUT(M,N,4)
      END DO
    END DO

    !H     DREHEN DER LOKALEN VEKTOREN IN X UND Y RICHTUNG
    !H     DIE VEKTOREN (DELTAXX,DELTAXY) UND (DELTAYY,DELTAYX)
    !H     SPANNEN DAS N/S, O/W SYSTEM AUF

    !       PRINT*,'V-DREHEN'

    DO M=2,ME-1
      DO N=2,NE-1

        IF ((UWOUT(M,N) .NE. SVAL) .AND. (VWOUT(M,N) .NE. SVAL)) THEN

          R_pol(1)=HOPLON(M,N)
          R_pol(2)=HOPLAT(M,N)

          UWOUT(M,N)= UWOUT(M,N)/REAL(rhoref_water, sp)
          VWOUT(M,N)= VWOUT(M,N)/REAL(rhoref_water, sp)

          R_i(:)=REAL(GIXYZ(:,2*M,2*N+1)-GIXYZ(:,2*M-1,2*N+1),sp)
          R_iabs=SQRT(SUM(R_i(:)**2))
          R_i(:)=R_i(:)/R_iabs
          CALL dek2geo(R_i,R_tmp,R_pol)
          R_poli=R_tmp(1:2)
          R_j(:)=REAL(GIXYZ(:,2*M,2*N)-GIXYZ(:,2*M,2*N+1),sp)
          R_jabs=SQRT(SUM(R_j(:)**2))
          R_j=R_j/R_jabs
          IF ( (R_jabs .LE. ZSMALL) .OR. (R_iabs .LE. ZSMALL) ) THEN
            WRITE(*,*) M,N, R_i, R_j
            STOP
          END IF
          CALL dek2geo(R_j,R_tmp,R_pol)
          R_polj=R_tmp(1:2)

          UKE(M,N) = (UWOUT(M,N)*R_poli(1)) + (VWOUT(M,N)*R_poli(2))
          VKE(M,N) = (UWOUT(M,N)*R_polj(1)) + (VWOUT(M,N)*R_polj(2))
          UE(M,N) = (UOUT(M,N)*R_poli(1)) + (VOUT(M,N)*R_poli(2))
          VE(M,N) = (UOUT(M,N)*R_polj(1)) + (VOUT(M,N)*R_polj(2))

          !H       ES WIRD PASSEND SKALIERT

          !H       BETRAG VOR DER DREHUNG
          SCHUGEO = SQRT(UWOUT(M,N)**2 + VWOUT(M,N)**2)

          !H       BETRAG NACH DER DREHUNG
          SCHUGR = SQRT(UKE(M,N)**2+VKE(M,N)**2)

          !H       QUOTIENT ALS KORREKTURFAKTOR

          VHOPE(M,N) = VKE(M,N)*SCHUGEO / SCHUGR
          !
          SCHUGEO = SQRT(UOUT(M,N)**2 + VOUT(M,N)**2)

          !H       BETRAG NACH DER DREHUNG
          SCHUGR = SQRT(UE(M,N)**2+VE(M,N)**2)

          VH(M,N) = VE(M,N)*SCHUGEO / SCHUGR

        ELSE

          VHOPE(M,N)=0.
          VH(M,N)=0.

        ENDIF

      END DO
    END DO


    !H    AUSGABE
    !h    PERIODISCHER RAND
    DO N=1,NE
      VHOPE(1,N)=VHOPE(ME-1,N)
      VHOPE(ME,N)=VHOPE(2,N)
      VH(1,N)=VH(ME-1,N)
      VH(ME,N)=VH(2,N)
    END DO

    IF (lbounds_exch_tp) THEN
      DO m=2,me-1
        ml=m
        mr=me+1-m
        VHOPE(ml,1) = -VHOPE(mr,3)           ! syncronise line 1 with line 3
        pq=0.5*(VHOPE(ml,2)-VHOPE(mr,2))
        vhope(ml,2)= pq                      ! syncronise line 2 with line 2
        vhope(mr,2)=-pq
        ! VHOPE(ml,2) = -VHOPE(mr,2)         ! syncronise line 2 with line 2
      END DO
    END IF


    IEXT(2)=53
    WRITE(53) IEXT
    WRITE(53) VHOPE
    IEXT(2)=166
    WRITE(63) IEXT
    WRITE(63) VH


  END DO

  CLOSE(52)
  CLOSE(53)
  CLOSE(54)
  CLOSE(55)
  CLOSE(56)
  CLOSE(57)
  CLOSE(58)
  CLOSE(59)

END PROGRAM FORCING

!H    ********************************************************

SUBROUTINE BLN2HOP(NUMFI,IE,JE,ME,NE,ALAT2,ALON2,ALAT1,ALON1, &
     SOE,SREG,LAND)

  !     BILINEAR INTERPOLATION FROM ONE GRID(IE,JE) TO ANOTHER GRID(ME,NE)

  !      PARAMETER(NB=3,MB=3)
  USE mo_kind, ONLY: i4, sp
  USE trigonometry
  IMPLICIT NONE
  INTEGER(i4), INTENT(in) :: ie, je, me, ne, numfi
  REAL(sp) ::   SOE(IE+2,JE+2,NUMFI), &
       SREG(ME,NE,NUMFI),G(IE+2,JE+2), &
       ALAT2(IE+2,JE+2),ALON2(IE+2,JE+2), &
       ALAT1(ME,NE),ALON1(ME,NE), &
       ALPHA,BETA, &
       HX(IE+2,JE+2,NUMFI),HHX(IE+2,JE+2,NUMFI), &
       LAND(IE+2,JE+2), &
       HG(IE+2,JE+2)

  INTEGER(i4) :: i, j, l, iter, jo, ju, il, ir, nland, m, n, jlu, ilu
  REAL(sp) :: rsumg, wwwalpha, wwwbeta

  !     DIFFUSION INTO LAND (NEW 12/99)

  DO L=1,NUMFI
    DO J=1,JE+2
      DO I=1,IE+2
        HHX(I,J,L)=SOE(I,J,L)
        !HH      X                 *((LAND(I,J)-1.))*(-1.)
      END DO
    END DO
  END DO

  DO J=1,JE+2
    DO I=1,IE+2
      HG(I,J)=1.
      IF (LAND(I,J) .ge. 0.5) HG(I,J)=1.E-12
    END DO
  END DO

  DO ITER=1,100
    DO J=1,JE+2
      DO I=1,IE+2
        G(I,J)=HG(I,J)
        DO L=1,NUMFI
          HX(I,J,l)=HHX(I,J,l)
        END DO
      END DO
    END DO
    DO J=1,JE+2
      JO=MAX(J-1,1)
      JU=MIN(J+1,JE+2)
      DO I=1,IE+2
        IL=I-1
        IF (IL .LT. 1)IL=IE
        IR=I+1
        IF (IR .GT. IE+2)IR=3
        RSUMG=0.
        IF (LAND(I,J) .GE. 0.5) THEN
          RSUMG=(4.*G(I,J) &
               +G(IL,J)+G(IR,J)+G(I,JO)+G(I,JU))/8.

          HG(I,J) = MIN(RSUMG,0.125_sp)

          DO L=1,NUMFI
            HHX(I,J,L)=(4.*HX(I,J,L)*G(I,J) &
                 +HX(IL,J,L)*G(IL,J) &
                 +HX(IR,J,L)*G(IR,J) &
                 +HX(I,JO,L)*G(I,JO) &
                 +HX(I,JU,L)*G(I,JU))/8.

            HHX(I,J,L)=HHX(I,J,L)/RSUMG
          END DO
        END IF
      END DO
    END DO
    nland=0
    DO I=1,IE
      DO j=1,je
        IF (HG(i,j) .LE. 2.E-12) nland=nland+1
      END DO
    END DO
    !          PRINT*,iter,' Nland: ', nland
  END DO


  DO J=1,JE+2
    DO I=1,IE+2
      G(i,j)=HG(i,j)
      DO L=1,NUMFI
        SOE(I,J,L)=HHX(I,J,L)
        IF (ABS(SOE(I,J,L)) .LT. 1.E-12) SOE(I,J,L)=0.
      END DO
    END DO
  END DO

  !
  DO M=1,ME
    DO N=1,NE

      !        PUNKT RECHTS OBEN
      DO J=1,JE+2
        IF (ALAT2(2,J) .GE. ALAT1(M,N)) JO=J
      END DO
      DO I=1,IE+2
        IF (ALON2(I,2) .LE. ALON1(M,N)) IL=I
      END DO

      !       IL=NINT(ALON1(M,N)*DLAMDAI+0.499999999)+1

      !            PRINT*,'WE ARE AT THE MINIMUM ',M,N,ALAT1(M,N)*GRARAD
      !     X                                    ,IL,JO,ALAT2(IL,JO)*GRARAD

      !       PUNKT RECHTS OBEN --> LINKS UNTEN

      JLU=JO+1
      ILU=IL

      WWWALPHA=ALAT2(ILU,JLU)-ALAT1(M,N)
      WWWBETA=ALON2(ILU,JLU)-ALON1(M,N)
      IF (WWWBETA .GE. PI) WWWBETA=WWWBETA-2.*PI
      IF (WWWBETA .LE. -PI) WWWBETA=WWWBETA+2.*PI

      ALPHA=(WWWALPHA)/(ALAT2(ILU,JLU)-ALAT2(ILU,JLU-1))
      BETA=(WWWBETA)/(ALON2(ILU,JLU)-ALON2(ILU+1,JLU))
      DO I=1,NUMFI
        SREG(M,N,I)=ALPHA*BETA*SOE(ILU+1,JLU-1,I)*G(ILU+1,JLU-1) &
             +(1.-ALPHA)*(1.-BETA)*SOE(ILU,JLU,I)*G(ILU,JLU) &
             +(1.-ALPHA)*(BETA)*SOE(ILU+1,JLU,I)*G(ILU+1,JLU) &
             +(ALPHA)*(1.-BETA)*SOE(ILU,JLU-1,I)*G(ILU,JLU-1)
        SREG(M,N,I)=SREG(M,N,I)/ &
             (ALPHA*BETA*G(ILU+1,JLU-1) &
             +(1.-ALPHA)*(1.-BETA)*G(ILU,JLU) &
             +(1.-ALPHA)*(BETA)*G(ILU+1,JLU) &
             +(ALPHA)*(1.-BETA)*G(ILU,JLU-1))
      END DO
      !         *************************************************
    END DO
  END DO
  RETURN
END SUBROUTINE BLN2HOP


SUBROUTINE XYZPOL(XYZ,R,ALAM,PHI)
  USE mo_kind, ONLY: sp

  IMPLICIT NONE

  REAL(sp) :: R, PHI, ALAM, XYZ(3)

  R = SQRT(XYZ(1)*XYZ(1)+XYZ(2)*XYZ(2)+XYZ(3)*XYZ(3))
  PHI = ASIN(XYZ(3)/R)
  ALAM = ATAN2(XYZ(1),XYZ(2))
  RETURN
END SUBROUTINE XYZPOL


SUBROUTINE POLXYZ(R,ALAM,PHI,XYZ)
  USE mo_kind, ONLY: sp

  IMPLICIT NONE

  REAL(sp) :: R, ALAM, PHI, XYZ(3), ZC

  XYZ(3) = R * SIN(PHI)
  ZC = R * COS(PHI)
  XYZ(1) = ZC * COS(ALAM)
  XYZ(2) = ZC * SIN(ALAM)
  RETURN
END SUBROUTINE POLXYZ
!===========================================================================

SUBROUTINE geo2dek(vc,vg,p)

  ! convert geographical vector components to cartesian

  USE mo_kind, ONLY: sp

  IMPLICIT NONE

  REAL(sp) :: vc(3), vg(3), p(2), lam, phi
  ! vc:  x,y,z vector components
  ! vg:  lam, phi, r vector components
  !  p:  lam, phi (in radias!)

  lam=p(1)
  phi=p(2)
  vc(1)= vg(3) * COS(phi) * COS(lam) - vg(1) * SIN(lam) &
       - vg(2) * SIN(phi) * COS(lam)
  vc(2)= vg(3) * COS(phi) * SIN(lam) + vg(1) * COS(lam) &
       - vg(2) * SIN(phi) * SIN(lam)
  vc(3)= vg(3) * SIN(phi) + vg(2) * COS(phi)
  RETURN
END SUBROUTINE geo2dek
!===========================================================================

SUBROUTINE dek2geo(vc,vg,p)

  ! convert cartesian vector components to geographical

  USE mo_kind, ONLY: sp

  IMPLICIT NONE

  REAL(sp) :: vc(3), vg(3), p(2), lam, phi
  ! vc:  x,y,z vector components
  ! vg:  lam, phi, r vector components
  !  p:  lam, phi (in radians!)

  lam = p(1)
  phi = p(2)
  vg(1) = -vc(1) * SIN(lam) + vc(2) * COS(lam)
  vg(2) = vc(3) * COS(phi) - vc(1) * SIN(phi) * COS(lam) &
       - vc(2) * SIN(phi) * SIN(lam)
  vg(3) = vc(1) * COS(phi) * COS(lam) + vc(2) * COS(phi) * SIN(lam) &
       + vc(3) * SIN(phi)
  !       write(*,*) 'vg=',vg(:)
  RETURN
END SUBROUTINE dek2geo

SUBROUTINE day_of_year2day_of_month(day_of_year, length_of_year_in_days, &
     month, day_of_month)

  USE mo_kind, ONLY: i4
  IMPLICIT NONE
  INTEGER(i4), INTENT(in) :: day_of_year, length_of_year_in_days
  INTEGER(i4), INTENT(out) :: month, day_of_month
  INTEGER(i4), PARAMETER :: days_per_month(12,2) = RESHAPE( &
       source=(/ 31, 28, 31, 30, 31, 30, &
       31, 31, 30, 31, 30, 31, &
       31, 29, 31, 30, 31, 30, &
       31, 31, 30, 31, 30, 31 /), shape=(/ 12, 2 /) )
  INTEGER(i4) :: days_aggregate(1:13)
  INTEGER(i4) :: year_idx, i, mon_count
  IF (length_of_year_in_days == 365) THEN
    year_idx = 1
  ELSE
    year_idx = 2
  END IF
  DO i = 1,13
    days_aggregate(i) = SUM(days_per_month(1:(i - 1), year_idx))
  END DO
  DO mon_count = 1, 12
    IF (day_of_year <= days_aggregate(mon_count + 1)) EXIT
  END DO
  month = mon_count
  day_of_month = day_of_year - days_aggregate(mon_count)
END SUBROUTINE day_of_year2day_of_month
