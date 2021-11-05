PROGRAM FORCING
  USE mo_kind, ONLY: sp, wp
  USE mo_planetary_constants, ONLY: rhoref_water
  IMPLICIT NONE
  INTEGER, PARAMETER :: IE=320, JE=160
  INTEGER :: ME, NE
  REAL(wp) cpuanf, cpuend
  REAL(wp), PARAMETER :: pi = 3.1415927_wp, grarad = 180._wp/pi
  REAL(wp) deltxx, deltxy, deltyx, deltyy
  REAL(wp) GAULAT(IE+2,JE+2),GAULON(IE+2,JE+2)
  REAL(sp) LAND(IE+2,JE+2), DUMMY(IE,JE)
  REAL(wp) GAUSSLA(JE),GAUSSLO(IE)
  INTEGER lday
  NAMELIST /forcingparams/ me, ne, lday
  REAL(sp) :: FIIN(IE+2,JE+2,7)

  REAL(sp) :: UWIN(IE+2,JE+2), &
       VWIN(IE+2,JE+2), &
       TEMIN(IE+2,JE+2), &
       ! SSTIN(IE+2,JE+2),
       PRESSIN(IE+2,JE+2), &
       TDEWIN(IE+2,JE+2), &
       SWRADIN(IE+2,JE+2), &
       DWLWIN(IE+2,JE+2), &
       PRATEIN(IE+2,JE+2), &
       U10IN(IE+2,JE+2)

  REAL(wp), ALLOCATABLE :: UWOUT(:,:), VWOUT(:,:), TEMOUT(:,:), &
       ! SSTOUT(:,:), &
       PRESSOUT(:,:), TDEWOUT(:,:), SWRADOUT(:,:), DWLWOUT(:,:), &
       PRATEOUT(:,:), U10OUT(:,:), UKE(:,:),VKE(:,:), &
       UHOPE(:,:), VHOPE(:,:), HOPLAT(:,:),HOPLON(:,:)

  REAL(sp), ALLOCATABLE :: GILA(:,:),GIPH(:,:)
  REAL(wp), ALLOCATABLE :: FIOUT(:,:,:)

  REAL(wp) SVAL, promeru, promerv, schugeo, schugr
  INTEGER, PARAMETER :: anta_lun=81
  INTEGER :: i, iday, if1, if2, if3, if4, iu1, iu2, iu3, iu4, &
       iv1, iv2, iv3, iv4, j, m, n

  DATA SVAL / 0._wp/

  !     ECMWF GITTER
  DATA GAUSSLA/89.142_wp, 88.029_wp, 86.911_wp, 85.791_wp, 84.670_wp, 83.549_wp, &
       82.428_wp, 81.307_wp, 80.185_wp, 79.064_wp, 77.943_wp, 76.821_wp, &
       75.700_wp, 74.578_wp, 73.457_wp, 72.336_wp, 71.214_wp, 70.093_wp, &
       68.971_wp, 67.850_wp, 66.728_wp, 65.607_wp, 64.485_wp, 63.364_wp, &
       62.242_wp, 61.121_wp, 60.000_wp, 58.878_wp, 57.757_wp, 56.635_wp, &
       55.514_wp, 54.392_wp, 53.271_wp, 52.149_wp, 51.028_wp, 49.906_wp, &
       48.785_wp, 47.663_wp, 46.542_wp, 45.420_wp, 44.299_wp, 43.177_wp, &
       42.056_wp, 40.934_wp, 39.813_wp, 38.691_wp, 37.570_wp, 36.448_wp, &
       35.327_wp, 34.205_wp, 33.084_wp, 31.962_wp, 30.841_wp, 29.719_wp, &
       28.598_wp, 27.476_wp, 26.355_wp, 25.234_wp, 24.112_wp, 22.991_wp, &
       21.869_wp, 20.748_wp, 19.626_wp, 18.505_wp, 17.383_wp, 16.262_wp, &
       15.140_wp, 14.019_wp, 12.897_wp, 11.776_wp, 10.654_wp,  9.533_wp, &
       8.411_wp,  7.290_wp,  6.168_wp,  5.047_wp,  3.925_wp,  2.804_wp, &
       1.682_wp,  0.561_wp, -0.561_wp, -1.682_wp, -2.804_wp, -3.925_wp, &
       -5.047_wp, -6.168_wp, -7.290_wp, -8.411_wp, -9.533_wp,-10.654_wp, &
       -11.776_wp,-12.897_wp,-14.019_wp,-15.140_wp,-16.262_wp,-17.383_wp, &
       -18.505_wp,-19.626_wp,-20.748_wp,-21.869_wp,-22.991_wp,-24.112_wp, &
       -25.234_wp,-26.355_wp,-27.476_wp,-28.598_wp,-29.719_wp,-30.841_wp, &
       -31.962_wp,-33.084_wp,-34.205_wp,-35.327_wp,-36.448_wp,-37.570_wp, &
       -38.691_wp,-39.813_wp,-40.934_wp,-42.056_wp,-43.177_wp,-44.299_wp, &
       -45.420_wp,-46.542_wp,-47.663_wp,-48.785_wp,-49.906_wp,-51.028_wp, &
       -52.149_wp,-53.271_wp,-54.392_wp,-55.514_wp,-56.635_wp,-57.757_wp, &
       -58.878_wp,-60.000_wp,-61.121_wp,-62.242_wp,-63.364_wp,-64.485_wp, &
       -65.607_wp,-66.728_wp,-67.850_wp,-68.971_wp,-70.093_wp,-71.214_wp, &
       -72.336_wp,-73.457_wp,-74.578_wp,-75.700_wp,-76.821_wp,-77.943_wp, &
       -79.064_wp,-80.185_wp,-81.307_wp,-82.428_wp,-83.549_wp,-84.670_wp, &
       -85.791_wp,-86.911_wp,-88.029_wp,-89.142_wp/

  me = -1
  ne = -1
  lday = -1
  READ (*, nml=forcingparams)

  IF (me < 1 .OR. ne < 1 .OR. lday < 1) THEN
    PRINT *, 'Only positive integral values may be used for me, ne and lday!'
    STOP
  END IF
  ALLOCATE(UWOUT(ME,NE), VWOUT(ME,NE), TEMOUT(ME,NE), &
       ! SSTOUT(ME,NE), &
       PRESSOUT(ME,NE), TDEWOUT(ME,NE), SWRADOUT(ME,NE), DWLWOUT(ME,NE), &
       PRATEOUT(ME,NE), U10OUT(ME,NE), UKE(ME,NE),VKE(ME,NE), &
       UHOPE(ME,NE), VHOPE(ME,NE), HOPLAT(ME,NE), HOPLON(ME,NE))
  ALLOCATE(GILA(2*ME,2*NE),GIPH(2*ME,2*NE))
  ALLOCATE(FIOUT(ME, NE, 7))

  GAUSSLO(1)=0._wp
  DO I=2,IE
    GAUSSLO(I)=GAUSSLO(I-1)+1.125_wp
  ENDDO

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
    GAULON(1,J)=GAULON(IE+1,J)-360._wp
    GAULON(IE+2,J)=GAULON(2,J)+360._wp
  ENDDO

  !     POLE (SCHUMMELN: OBEN UND UNTEN EINE REIHE HINZU)
  DO I=1,IE+2

    GAULAT(I,1)=180._wp-GAULAT(I,2)
    GAULAT(I,JE+2)=-180._wp-GAULAT(I,JE+1)
    GAULON(I,1)=GAULON(I,2)
    GAULON(I,JE+2)=GAULON(I,JE+1)

  ENDDO

  !     NOCHMAL ZYKL.RAND OST/WEST
  DO J=2,JE+1
    GAULAT(1,J)=GAULAT(IE+1,J)
    GAULAT(IE+2,J)=GAULAT(2,J)
    GAULON(1,J)=GAULON(IE+1,J)-360._wp
    GAULON(IE+2,J)=GAULON(2,J)+360._wp
  ENDDO

  DO I=1,IE+2
    DO J=1,JE+2
      GAULAT(I,J)=GAULAT(I,J)/GRARAD
      GAULON(I,J)=GAULON(I,J)/GRARAD
    ENDDO
  ENDDO

  OPEN(anta_lun,FILE='anta.ext4',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(52,FILE='GIWIX.ext4' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(53,FILE='GIWIY.ext4' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(54,FILE='GITEM.ext4' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(55,FILE='GITDEW.ext4' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(56,FILE='GIPREC.ext4' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  !      OPEN(57,FILE='GIDWLW',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(57,FILE='GICLOUD.ext4' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(58,FILE='GIU10.ext4' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(59,FILE='GISWRAD.ext4' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  !      OPEN(60,FILE='GISST',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(61,FILE='GIPRESS.ext4' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(10,FILE='land_sea_mask.ECMWF.ext4' &
       ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

  READ(anta_lun)IF1,IF2,IF3,IF4
  READ(anta_lun) GILA
  READ(anta_lun)IF1,IF2,IF3,IF4
  READ(anta_lun) GIPH

  !H    TAGESWERTE
  !cc      READ(*,*)LDAY
  READ(10)IF1,IF2,IF3,IF4
  READ(10)((DUMMY(I,J),I=1,IE),J=1,JE)
  DO I=1,IE
    DO J=1,JE
      LAND(I+1,J+1)=DUMMY(I,J)
      IF (land(I+1,J+1) .GT. 1.e-4_sp) land(i+1, j+1) = 1._sp
    ENDDO
  ENDDO
  !     ZYKL.RAND OST/WEST/NORDPOL/SUEDPOL
  DO I=1,IE+2
    LAND(I,1)=0._sp
    LAND(I,JE+2)=1._sp
  ENDDO
  DO J=1,JE+2
    LAND(1,J)=LAND(IE+1,J)
    LAND(IE+2,J)=LAND(2,J)
  ENDDO

  DO IDAY=1,LDAY
    Print*,'DATUM ',IDAY
    !H    READ INPUT FILES

    !c      OPEN(2,FILE='sfcfullfilt2416704.no_f_norcw.ext4'
    OPEN(2,FILE='2m_temperature.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(2)IF1,IF2,IF3,IF4
    READ(2)((DUMMY(I,J),I=1,IE),J=1,JE)
    DO I=1,IE
      DO J=1,JE
        TEMIN(I+1,J+1)=DUMMY(I,J)-273.15_sp
      ENDDO
    ENDDO

    !c      OPEN(11,FILE='sfcfullfilt2416804.no_f_norcw.ext4'
    OPEN(11,FILE='2m_dewpoint_temperature.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(11)IF1,IF2,IF3,IF4
    READ(11)((DUMMY(I,J),I=1,IE),J=1,JE)
    DO I=1,IE
      DO J=1,JE
        TDEWIN(I+1,J+1)=DUMMY(I,J)
      ENDDO
    ENDDO

    !c      OPEN(12,FILE='sfcfullfilt2442304.no_f_norcw.ext4'
    OPEN(12,FILE='total_precipitation.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(12)IF1,IF2,IF3,IF4
    READ(12)((DUMMY(I,J),I=1,IE),J=1,JE)
    DO I=1,IE
      DO J=1,JE
        PRATEIN(I+1,J+1)=MAX(0._sp, DUMMY(I,J))
      ENDDO
    ENDDO


    !c      OPEN(13,FILE='sfcfullfilt2417904.ext4'
    !c      OPEN(13,FILE='sfcfullfilt2417904.ext4'
    !c     X    ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    !c      READ(13)IF1,IF2,IF3,IF4
    !c      READ(13)((DUMMY(I,J),I=1,IE),J=1,JE)
    !c      DO I=1,IE
    !c        DO J=1,JE
    !c          DWLWIN(I+1,J+1)=DUMMY(I,J)
    !c        ENDDO
    !c      ENDDO

    !c      OPEN(13,FILE='sfcfullfilt2417904.ext4'
    OPEN(13,FILE='total_cloud_cover.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(13)IF1,IF2,IF3,IF4
    READ(13)((DUMMY(I,J),I=1,IE),J=1,JE)
    DO I=1,IE
      DO J=1,JE
        DWLWIN(I+1,J+1)=DUMMY(I,J)
      ENDDO
    ENDDO


    !c      OPEN(14,FILE='sfcfullfilt2466504.no_f_norcw.ext4'
    OPEN(14,FILE='scalar_wind.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(14)IF1,IF2,IF3,IF4
    READ(14)((DUMMY(I,J),I=1,IE),J=1,JE)
    DO I=1,IE
      DO J=1,JE
        U10IN(I+1,J+1)=DUMMY(I,J)
      ENDDO
    ENDDO

    !c      OPEN(15,FILE='sfcfullfilt2427604.no_f_norcw.ext4'
    OPEN(15,FILE='total_solar_radiation.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(15)IF1,IF2,IF3,IF4
    READ(15)((DUMMY(I,J),I=1,IE),J=1,JE)
    DO I=1,IE
      DO J=1,JE
        SWRADIN(I+1,J+1)=DUMMY(I,J)
      ENDDO
    ENDDO

    !c      OPEN(3,FILE='sfcfullfilt2418004.no_f_norcw.ext4'
    OPEN(3,FILE='east_west_stress.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(3)IU1,IU2,IU3,IU4
    READ(3)((DUMMY(I,J),I=1,IE),J=1,JE)
    DO I=1,IE
      DO J=1,JE
        UWIN(I+1,J+1)= DUMMY(I,J)
      ENDDO
    ENDDO


    !c      OPEN(4,FILE='sfcfullfilt2418104.no_f_norcw.ext4'
    OPEN(4,FILE='north_south_stress.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(4)IV1,IV2,IV3,IV4
    READ(4)((DUMMY(I,J),I=1,IE),J=1,JE)
    DO I=1,IE
      DO J=1,JE
        VWIN(I+1,J+1)= DUMMY(I,J)
      ENDDO
    ENDDO

    !c      OPEN(16,FILE='sfcfullfilt2415104.no_f_norcw.ext4'
    OPEN(16,FILE='mean_sea_level_pressure.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(16)IV1,IV2,IV3,IV4
    READ(16)((DUMMY(I,J),I=1,IE),J=1,JE)
    DO I=1,IE
      DO J=1,JE
        PRESSIN(I+1,J+1)= DUMMY(I,J)
      ENDDO
    ENDDO


    !      OPEN(17,FILE='sfcfullfilt2413904.no_f_norcw.ext4'
    !     X    ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    !      READ(17)IV1,IV2,IV3,IV4
    !      READ(17)((DUMMY(I,J),I=1,IE),J=1,JE)
    !      DO I=1,IE
    !        DO J=1,JE
    !          SSTIN(I+1,J+1)= DUMMY(I,J)
    !        ENDDO
    !      ENDDO


    !     DATENFELDER VORRBEITEN
    !     ZYKL.RAND OST/WEST

    DO J=1,JE+2
      UWIN(1,J)=UWIN(IE+1,J)
      UWIN(IE+2,J)=UWIN(2,J)

      VWIN(1,J)=VWIN(IE+1,J)
      VWIN(IE+2,J)=VWIN(2,J)

      TEMIN(1,J)=TEMIN(IE+1,J)
      TEMIN(IE+2,J)=TEMIN(2,J)

      PRESSIN(1,J)=PRESSIN(IE+1,J)
      PRESSIN(IE+2,J)=PRESSIN(2,J)

      !        SSTIN(1,J)=SSTIN(IE+1,J)
      !        SSTIN(IE+2,J)=SSTIN(2,J)

      TDEWIN(1,J)=TDEWIN(IE+1,J)
      TDEWIN(IE+2,J)=TDEWIN(2,J)

      PRATEIN(1,J)=PRATEIN(IE+1,J)
      PRATEIN(IE+2,J)=PRATEIN(2,J)

      DWLWIN(1,J)=DWLWIN(IE+1,J)
      DWLWIN(IE+2,J)=DWLWIN(2,J)

      U10IN(1,J)=U10IN(IE+1,J)
      U10IN(IE+2,J)=U10IN(2,J)

      SWRADIN(1,J)=SWRADIN(IE+1,J)
      SWRADIN(IE+2,J)=SWRADIN(2,J)

    ENDDO

    !     POLE (EINE REIHE HINZU)
    DO I=1,IE+2

      UWIN(I,1)=UWIN(I,2)
      UWIN(I,JE+2)=UWIN(I,JE+1)

      VWIN(I,1)=VWIN(I,2)
      VWIN(I,JE+2)=VWIN(I,JE+1)

      TEMIN(I,1)=TEMIN(I,2)
      TEMIN(I,JE+2)=TEMIN(I,JE+1)

      PRESSIN(I,1)=PRESSIN(I,2)
      PRESSIN(I,JE+2)=PRESSIN(I,JE+1)

      !c          SSTIN(I,1)=SSTIN(I,2)
      !c          SSTIN(I,JE+2)=SSTIN(I,JE+1)

      TDEWIN(I,1)=TDEWIN(I,2)
      TDEWIN(I,JE+2)=TDEWIN(I,JE+1)

      PRATEIN(I,1)=PRATEIN(I,2)
      PRATEIN(I,JE+2)=PRATEIN(I,JE+1)

      DWLWIN(I,1)=DWLWIN(I,2)
      DWLWIN(I,JE+2)=DWLWIN(I,JE+1)

      U10IN(I,1)=U10IN(I,2)
      U10IN(I,JE+2)=U10IN(I,JE+1)

      SWRADIN(I,1)=SWRADIN(I,2)
      SWRADIN(I,JE+2)=SWRADIN(I,JE+1)

    ENDDO

    PRINT*,'Datum ',iday,' EINLESEN UEBERLEBT'

    !H    AN DIESER STELLE WIRD JEWEILS FUER U,V UND SKALARE GROESSEN
    !H    INTERPOLIERT, UND ZWAR AM JEWEILS DAZUGEHOERENDEN GITTERPUNKT
    !H    1. FUER SKALARE WIE GEHABT AN DRUCKPUNKT
    PRINT*,'SKALARE GROESSE INTERPOLIEREN'

    DO M=1,ME
      DO N=1,NE
        HOPLON(M,N) = REAL(GILA(2*M,2*N), wp)
        HOPLAT(M,N) = REAL(GIPH(2*M,2*N), wp)
        IF (HOPLON(M,N) .GT. (2._wp*PI)) HOPLON(M,N)=HOPLON(M,N)-(2._wp*PI)
        IF (HOPLON(M,N) .LT. 0._wp) HOPLON(M,N)=HOPLON(M,N)+(2._wp*PI)
      ENDDO
    ENDDO

    CALL CPU_TIME(cpuanf)
    DO I=1,IE+2
      DO J=1,JE+2
        FIIN(I,J,1)=TEMIN(I,J)
        FIIN(I,J,2)=TDEWIN(I,J)
        FIIN(I,J,3)=PRATEIN(I,J)
        FIIN(I,J,4)=DWLWIN(I,J)
        FIIN(I,J,5)=U10IN(I,J)
        FIIN(I,J,6)=SWRADIN(I,J)
        FIIN(I,J,7)=PRESSIN(I,J)
        !c            FIIN(I,J,8)=SSTIN(I,J)
      ENDDO
    ENDDO
    CALL BLN2HOP(7,IE,JE,ME,NE,GAULAT,GAULON,HOPLAT,HOPLON, &
         FIIN,FIOUT,LAND)

    DO M=1,ME
      DO N=1,NE
        TEMOUT(M,N)=FIOUT(M,N,1)
        TDEWOUT(M,N)=FIOUT(M,N,2)
        PRATEOUT(M,N)=FIOUT(M,N,3)
        DWLWOUT(M,N)=FIOUT(M,N,4)
        U10OUT(M,N)=FIOUT(M,N,5)
        SWRADOUT(M,N)=FIOUT(M,N,6)
        PRESSOUT(M,N)=FIOUT(M,N,7)
        !c            SSTOUT(M,N)=FIOUT(M,N,8)

      ENDDO
    ENDDO
    CALL CPU_TIME(cpuend)
1   PRINT *, 'CALL BLN2HOP DAUERT ', cpuend - cpuanf

    !H    AUSGABE
    !h    PERIODISCHER RAND
    DO N=1,NE

      TEMOUT(1,N)=TEMOUT(ME-1,N)
      TEMOUT(ME,N)=TEMOUT(2,N)

      !c        SSTOUT(1,N)=SStOUT(ME-1,N)
      !c        SSTOUT(ME,N)=SSTOUT(2,N)

      PRESSOUT(1,N)=PRESSOUT(ME-1,N)
      PRESSOUT(ME,N)=PRESSOUT(2,N)

      TDEWOUT(1,N)=TDEWOUT(ME-1,N)
      TDEWOUT(ME,N)=TDEWOUT(2,N)

      PRATEOUT(1,N)=PRATEOUT(ME-1,N)
      PRATEOUT(ME,N)=PRATEOUT(2,N)

      DWLWOUT(1,N)=DWLWOUT(ME-1,N)
      DWLWOUT(ME,N)=DWLWOUT(2,N)

      U10OUT(1,N)=U10OUT(ME-1,N)
      U10OUT(ME,N)=U10OUT(2,N)

      SWRADOUT(1,N)=SWRADOUT(ME-1,N)
      SWRADOUT(ME,N)=SWRADOUT(2,N)

    ENDDO





    IF2=167
    WRITE(54) IDAY,IF2,IF3,(ME*NE)
    WRITE(54) TEMOUT
    IF2=168
    WRITE(55) IDAY,IF2,IF3,(ME*NE)
    WRITE(55) TDEWOUT
    IF2=423
    WRITE(56) IDAY,IF2,IF3,(ME*NE)
    WRITE(56) PRATEOUT
    IF2=177
    WRITE(57) IDAY,IF2,IF3,(ME*NE)
    WRITE(57) DWLWOUT
    IF2=665
    WRITE(58) IDAY,IF2,IF3,(ME*NE)
    WRITE(58) U10OUT
    IF2=276
    WRITE(59) IDAY,IF2,IF3,(ME*NE)
    WRITE(59) SWRADOUT
    IF2=139
    !c      WRITE(60) IDAY,IF2,IF3,(ME*NE)
    !c      WRITE(60) SSTOUT
    IF2=151
    WRITE(61) IDAY,IF2,IF3,(ME*NE)
    WRITE(61) PRESSOUT


    !H    2. FUER VECTOR U KOMPONENTE AM U-PUNKT
    PRINT*,'U-KOMPOMENTE INTERPOLIEREN'

    DO M=2,ME-1
      DO N=1,NE
        HOPLON(M,N)=REAL(GILA((2*M)+1,2*N), wp)
        HOPLAT(M,N)=REAL(GIPH((2*M)+1,2*N), wp)
        IF (HOPLON(M,N) .GT. (2._wp*PI)) HOPLON(M,N)=HOPLON(M,N)-(2._wp*PI)
        IF (HOPLON(M,N) .LT. 0._wp) HOPLON(M,N)=HOPLON(M,N)+(2._wp*PI)
      ENDDO
    ENDDO

    DO I=1,IE+2
      DO J=1,JE+2
        FIIN(I,J,1)=UWIN(I,J)
        FIIN(I,J,2)=VWIN(I,J)
      ENDDO
    ENDDO

    CALL BLN2HOP(2,IE,JE,ME,NE,GAULAT,GAULON,HOPLAT,HOPLON, &
         FIIN,FIOUT,LAND)

    DO M=1,ME
      DO N=1,NE
        UWOUT(M,N)=FIOUT(M,N,1)
        VWOUT(M,N)=FIOUT(M,N,2)
      ENDDO
    ENDDO

    !H     DREHEN DER LOKALEN VEKTOREN IN X UND Y RICHTUNG
    !H     DIE VEKTOREN (DELTAXX,DELTAXY) UND (DELTAYY,DELTAYX)
    !H     SPANNEN DAS N/S, O/W SYSTEM AUF
    !H     WINDSTRESS (N/kg)/rho

    PRINT*,'U-DREHEN'

    DO M=2,ME-1
      DO N=2,NE-1

        IF ((UWOUT(M,N) .NE. SVAL) .AND. (VWOUT(M,N) .NE. SVAL)) THEN

          UWOUT(M,N)=UWOUT(M,N)/rhoref_water
          VWOUT(M,N)=VWOUT(M,N)/rhoref_water

          DELTXX = REAL(GILA((2*M)+1,2*N)-GILA(2*M,2*N), wp)
          DELTYY = REAL(GIPH((2*M)+1,(2*N)-1)-GIPH(((2*M)+1),2*N), wp)
          DELTXY = REAL(GIPH((2*M)+1,2*N)-GIPH(2*M,2*N), wp)
          DELTYX = REAL(GILA((2*M)+1,2*N-1)-GILA((2*M)+1,2*N), wp)

          IF(DELTXX .GT. PI) DELTXX = DELTXX - 2._wp*PI
          IF(DELTXX .LT. -PI) DELTXX = DELTXX + 2._wp*PI
          IF(DELTYX .LT. -PI) DELTYX = DELTYX + 2._wp*PI
          IF(DELTYX .GT. PI) DELTYX = DELTYX - 2._wp*PI

          !H       TEST WINKELTREUE MERKATORPROJECTION

          PROMERU = REAL(COS(GIPH(((2*M)+1),2*N)), wp)
          DELTXX = DELTXX*PROMERU
          DELTYX = DELTYX*PROMERU


          !H       PROJEKTION VON UKE AUS RICHTUNG I,J IN RICHTUNG N/S, O/W
          !H       MITTELS SKALARPRODUKT = A1*B1+A2*B2+...

          UKE(M,N) = (UWOUT(M,N)*DELTXX) + (VWOUT(M,N)*DELTXY)
          VKE(M,N) = (UWOUT(M,N)*DELTYX) + (VWOUT(M,N)*DELTYY)

          !H       ES MUESSTE EIGENTLICH NOCH DURCH DEN BETRAG
          !H       VON DELTA GETEILT WERDEN DA DELTA KEIN EINHEITSVEKTOR IST

          !H       TEST
          IF ((DELTXX*DELTYX+DELTYY*DELTXY) .GT. 0.1_wp) THEN
            PRINT*,'ORTHOGONAL?',M,N,(DELTXX*DELTYX+DELTYY*DELTXY)
          ENDIF

          UKE(M,N) = UKE(M,N) / SQRT(DELTXX**2+DELTXY**2)
          VKE(M,N) = VKE(M,N) / SQRT(DELTYY**2+DELTYX**2)

          !H       ES WIRD PASSEND SKALIERT

          !H       BETRAG VOR DER DREHUNG
          SCHUGEO = SQRT(UWOUT(M,N)**2 + VWOUT(M,N)**2)

          !H       BETRAG NACH DER DREHUNG
          SCHUGR = SQRT(UKE(M,N)**2+VKE(M,N)**2)

          !H       QUOTIENT ALS KORREKTURFAKTOR

          UHOPE(M,N) = UKE(M,N)*SCHUGEO / SCHUGR

        ELSE

          UHOPE(M,N)=0._wp

        ENDIF

      ENDDO
    ENDDO


    !H    AUSGABE
    !h    PERIODISCHER RAND
    DO N=1,NE
      UHOPE(1,N)=UHOPE(ME-1,N)
      UHOPE(ME,N)=UHOPE(2,N)
    ENDDO

    IF2=180
    WRITE(52) IDAY,IF2,IF3,(ME*NE)
    WRITE(52) UHOPE



    !H    3. FUER VECTOR V KOMPONENTE AM V-PUNKT
    PRINT*,'V-KOMPOMENTE INTERPOLIEREN'

    DO M=1,ME
      DO N=2,NE-1
        HOPLON(M,N) = REAL(GILA(2*M,(2*N)+1), wp)
        HOPLAT(M,N) = REAL(GIPH(2*M,(2*N)+1), wp)
        IF (HOPLON(M,N) .GT. (2._wp*PI)) HOPLON(M,N)=HOPLON(M,N)-(2._wp*PI)
        IF (HOPLON(M,N) .LT. 0._wp) HOPLON(M,N)=HOPLON(M,N)+(2._wp*PI)
      ENDDO
    ENDDO

    DO I=1,IE+2
      DO J=1,JE+2
        FIIN(I,J,1)=UWIN(I,J)
        FIIN(I,J,2)=VWIN(I,J)
      ENDDO
    ENDDO

    CALL BLN2HOP(2,IE,JE,ME,NE,GAULAT,GAULON,HOPLAT,HOPLON, &
         FIIN,FIOUT,LAND)

    DO M=1,ME
      DO N=1,NE
        UWOUT(M,N)=FIOUT(M,N,1)
        VWOUT(M,N)=FIOUT(M,N,2)
      ENDDO
    ENDDO

    !H     DREHEN DER LOKALEN VEKTOREN IN X UND Y RICHTUNG
    !H     DIE VEKTOREN (DELTAXX,DELTAXY) UND (DELTAYY,DELTAYX)
    !H     SPANNEN DAS N/S, O/W SYSTEM AUF

    PRINT*,'V-DREHEN'

    DO M=2,ME-1
      DO N=2,NE-1

        IF ((UWOUT(M,N) .NE. SVAL) .AND. (VWOUT(M,N) .NE. SVAL)) THEN

          UWOUT(M,N)= UWOUT(M,N)/rhoref_water
          VWOUT(M,N)= VWOUT(M,N)/rhoref_water

          DELTXX = REAL(GILA(2*M,(2*N)+1)-GILA((2*M)-1,(2*N)+1), wp)
          DELTYY = REAL(GIPH(2*M,2*N)    -GIPH(2*M    ,(2*N)+1), wp)
          DELTXY = REAL(GIPH(2*M,(2*N)+1)-GIPH((2*M)-1,(2*N)+1), wp)
          DELTYX = REAL(GILA(2*M,2*N)    -GILA(2*M    ,(2*N)+1), wp)

          IF(DELTXX .GT. PI) DELTXX = DELTXX - 2._wp*PI
          IF(DELTXX .LT. -PI) DELTXX = DELTXX + 2._wp*PI
          IF(DELTYX .LT. -PI) DELTYX = DELTYX + 2._wp*PI
          IF(DELTYX .GT. PI) DELTYX = DELTYX - 2._wp*PI

          !H       TEST WINKELTREUE MERKATORPROJECTION

          PROMERV = REAL(COS(GIPH(2*M,(2*N)+1)), wp)
          DELTXX = DELTXX*PROMERV
          DELTYX = DELTYX*PROMERV


          !H       PROJEKTION VON VKE AUS RICHTUNG I,J IN RICHTUNG N/S, O/W
          !H       MITTELS SKALARPRODUKT = A1*B1+A2*B2+...

          UKE(M,N) = (UWOUT(M,N)*DELTXX) + (VWOUT(M,N)*DELTXY)
          VKE(M,N) = (UWOUT(M,N)*DELTYX) + (VWOUT(M,N)*DELTYY)

          !H       ES MUESSTE EIGENTLICH NOCH DURCH DEN BETRAG
          !H       VON DELTA GETEILT WERDEN DA DELTA KEIN EINHEITSVEKTOR IST

          !H       TEST
          IF ((DELTXX*DELTYX+DELTYY*DELTXY) .GT. 0.1_wp) THEN
            PRINT*,'ORTHOGONAL?',M,N,(DELTXX*DELTYX+DELTYY*DELTXY)
          ENDIF

          UKE(M,N) = UKE(M,N) / SQRT(DELTXX**2+DELTXY**2)
          VKE(M,N) = VKE(M,N) / SQRT(DELTYY**2+DELTYX**2)

          !H       ES WIRD PASSEND SKALIERT

          !H       BETRAG VOR DER DREHUNG
          SCHUGEO = SQRT(UWOUT(M,N)**2 + VWOUT(M,N)**2)

          !H       BETRAG NACH DER DREHUNG
          SCHUGR = SQRT(UKE(M,N)**2+VKE(M,N)**2)

          !H       QUOTIENT ALS KORREKTURFAKTOR

          VHOPE(M,N) = VKE(M,N)*SCHUGEO / SCHUGR


        ELSE

          VHOPE(M,N)=0._wp

        ENDIF

      ENDDO
    ENDDO


    !H    AUSGABE
    !h    PERIODISCHER RAND
    DO N=1,NE
      VHOPE(1,N)=VHOPE(ME-1,N)
      VHOPE(ME,N)=VHOPE(2,N)
    ENDDO

    if2=181
    WRITE(53) IDAY,IF2,IF3,(ME*NE)
    WRITE(53) VHOPE


  END DO

  STOP

CONTAINS
  SUBROUTINE BLN2HOP(NUMFI,IE,JE,ME,NE,alat2,alon2,alat1,alon1, &
       soe, sreg, land)
    USE mo_kind, ONLY: sp
    IMPLICIT NONE
    INTEGER, INTENT(in) :: numfi, ie, je, me, ne
    REAL(wp), INTENT(in) :: alat1(me, ne), alon1(me, ne), alon2(IE+2,JE+2), &
         alat2(ie+2,je+2)
    REAL(sp), INTENT(inout) :: soe(ie+2, je+2, numfi)
    REAL(wp), INTENT(inout) :: sreg(me,ne,numfi)
    REAL(sp), INTENT(in) :: land(ie+2,je+2)

    !     BILINEAR INTERPOLATION FROM ONE GRID(IE,JE) TO ANOTHER GRID(ME,NE)

    !      PARAMETER(NB=3,MB=3)

    REAL(wp) G(IE+2,JE+2), &
         ALPHA,BETA, &
         HX(IE+2,JE+2,NUMFI),HHX(IE+2,JE+2,NUMFI), &
         HG(IE+2,JE+2)
    REAL(wp) :: rsumg, wwwalpha, wwwbeta
    INTEGER :: i, il, ilu, ir, iter, j, jlu, jo, ju, l, m, n, nland
    !     DIFFUSION INTO LAND (NEW 12/99)

    DO L=1,NUMFI
      DO J=1,JE+2
        DO I=1,IE+2
          HHX(I,J,L) = REAL(SOE(I,J,L), wp)
          !HH      X                 *((LAND(I,J)-1.))*(-1.)
        ENDDO
      ENDDO
    ENDDO

    DO J=1,JE+2
      DO I=1,IE+2
        HG(I,J)=1._wp
        IF (LAND(I,J) .ge. 0.5_sp) HG(I,J)=1.E-12_wp
      ENDDO
    ENDDO

    DO ITER=1,100
      DO J=1,JE+2
        DO I=1,IE+2
          G(I,J)=HG(I,J)
          DO L=1,NUMFI
            HX(I,J,l)=HHX(I,J,l)
          ENDDO
        ENDDO
      ENDDO
      DO J=1,JE+2
        JO=MAX(J-1,1)
        JU=MIN(J+1,JE+2)
        DO I=1,IE+2
          IL=I-1
          IF(IL .LT. 1)IL=IE
          IR=I+1
          IF(IR .GT. IE+2)IR=3
          RSUMG = 0._wp
          IF(LAND(I,J) .GE. 0.5_sp)THEN
            RSUMG=(4._wp*G(I,J) &
                 +G(IL,J)+G(IR,J)+G(I,JO)+G(I,JU))/8._wp

            HG(I,J)=MIN(RSUMG,0.125_wp)

            DO L=1,NUMFI
              HHX(I,J,L)=(4._wp*HX(I,J,L)*G(I,J) &
                   +HX(IL,J,L)*G(IL,J) &
                   +HX(IR,J,L)*G(IR,J) &
                   +HX(I,JO,L)*G(I,JO) &
                   +HX(I,JU,L)*G(I,JU))/8._wp

              HHX(I,J,L)=HHX(I,J,L)/RSUMG
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      nland = 0
      Do I=1,IE
        DO j=1,je
          if (HG(i,j) .le. 2.E-12_wp) nland=nland+1
        ENDDO
      ENDDO
      PRINT*,iter,' Nland: ', nland
    ENDDO


    DO J=1,JE+2
      DO I=1,IE+2
        G(i,j)=HG(i,j)
        DO L=1,NUMFI
          SOE(I,J,L) = REAL(HHX(I,J,L), sp)
          IF (ABS(SOE(I,J,L)) .LT. 1.E-12_sp) SOE(I,J,L)=0._sp
        ENDDO
      ENDDO
    ENDDO

    !
    DO M=1,ME
      DO N=1,NE

        !        PUNKT RECHTS OBEN
        DO J=1,JE+2
          IF (ALAT2(2,J) .GE. ALAT1(M,N)) JO=J
        ENDDO
        DO I=1,IE+2
          IF (ALON2(I,2) .LE. ALON1(M,N)) IL=I
        ENDDO

        !       IL=NINT(ALON1(M,N)*DLAMDAI+0.499999999)+1

        !            PRINT*,'WE ARE AT THE MINIMUM ',M,N,ALAT1(M,N)*GRARAD
        !     X                                    ,IL,JO,ALAT2(IL,JO)*GRARAD

        !       PUNKT RECHTS OBEN --> LINKS UNTEN

        JLU=JO+1
        ILU=IL

        WWWALPHA=ALAT2(ILU,JLU)-ALAT1(M,N)
        WWWBETA=ALON2(ILU,JLU)-ALON1(M,N)
        IF(WWWBETA .GE. PI) WWWBETA=WWWBETA-2._wp*PI
        IF(WWWBETA .LE. -PI) WWWBETA=WWWBETA+2._wp*PI

        ALPHA=(WWWALPHA)/(ALAT2(ILU,JLU)-ALAT2(ILU,JLU-1))
        BETA=(WWWBETA)/(ALON2(ILU,JLU)-ALON2(ILU+1,JLU))
        DO I=1,NUMFI
          SREG(M,N,I)=ALPHA*BETA*REAL(SOE(ILU+1,JLU-1,I), wp)*G(ILU+1,JLU-1) &
               +(1._wp-ALPHA)*(1._wp-BETA)*REAL(SOE(ILU,JLU,I), wp)*G(ILU,JLU) &
               +(1._wp-ALPHA)*(BETA)*REAL(SOE(ILU+1,JLU,I), wp)*G(ILU+1,JLU) &
               +(ALPHA)*(1._wp-BETA)*REAL(SOE(ILU,JLU-1,I), wp)*G(ILU,JLU-1)
          SREG(M,N,I)=SREG(M,N,I)/ &
               (ALPHA*BETA*G(ILU+1,JLU-1) &
               +(1._wp-ALPHA)*(1._wp-BETA)*G(ILU,JLU) &
               +(1._wp-ALPHA)*(BETA)*G(ILU+1,JLU) &
               +(ALPHA)*(1._wp-BETA)*G(ILU,JLU-1))
        ENDDO
        !         *************************************************
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE BLN2HOP

END PROGRAM FORCING
