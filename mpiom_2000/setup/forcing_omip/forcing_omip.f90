PROGRAM FORCING
  USE mo_kind, ONLY: dp, sp, i4, i8, wp
  USE mo_constants, ONLY: aradtogra, api
  IMPLICIT NONE
  INTEGER, PARAMETER :: IE=320, JE=160
  INTEGER :: ME, NE

  !     SOME CONSTANTS
  REAL(sp) :: gaulat(ie+2,je+2), gaulon(IE+2,JE+2), &
       gaussla(je), gausslo(IE)
  REAL(sp) :: land(ie+2, je+2), dummy(ie, je)

  REAL(sp) :: fiin(ie+2, je+2, 7)

  REAL(sp) :: temin(ie+2, je+2), tdewin(ie+2, je+2), pratein(ie+2, je+2), &
       dwlwin(ie+2, je+2), u10in(ie+2, je+2), &
       uwin(ie+2, je+2), vwin(ie+2, je+2), swradin(ie+2, je+2), &
       pressin(ie+2, je+2)
       !c     X         SSTIN(IE+2,JE+2), &

  REAL(dp), ALLOCATABLE :: GILA(:,:),GIPH(:,:)
  REAL(wp), ALLOCATABLE :: HOPLAT(:,:),HOPLON(:,:)
  REAL(wp), ALLOCATABLE :: FIOUT(:,:,:)
  REAL(wp), ALLOCATABLE :: UWOUT(:,:), &
       VWOUT(:,:), &
       temout(:,:), &
       ! SSTOUT(:,:), &
       pressout(:,:), &
       tdewout(:,:), &
       swradout(:,:), &
       dwlwout(:,:), &
       prateout(:,:), &
       u10out(:,:)
  REAL(wp), ALLOCATABLE :: uke(:,:), vke(:,:)
  REAL(wp), ALLOCATABLE :: uhope(:,:), vhope(:,:)
  LOGICAL :: bounds_exch_tp
  REAL(wp) :: deltxx, deltxy, deltyx, deltyy, promeru, promerv, pq, schugeo, &
       schugr
  INTEGER :: i, iday, lday, j, m, ml, mr, n
  INTEGER(i4) :: ext4_hdr(4), ext4_hdr_temp(4)
  INTEGER(i8) :: ext8_hdr(4)
  INTEGER :: verbose
  INTEGER, PARAMETER :: io_in_anta=81, io_out_giwix=52, io_out_giwiy=53, &
       io_in_lsm=10, io_out_gitem=54, io_out_gitdew=55, io_out_giprec=56, &
       io_out_gicloud=57, io_out_giu10=58, io_out_swrad=59, &
       io_out_gipress=61

  NAMELIST /gridparams/ me, ne, bounds_exch_tp, verbose

  REAL(wp) SVAL

  DATA SVAL / 0._wp/

  !     ECMWF GITTER
  DATA GAUSSLA/&
       89.142_sp, 88.029_sp, 86.911_sp, 85.791_sp, 84.670_sp, 83.549_sp, &
       82.428_sp, 81.307_sp, 80.185_sp, 79.064_sp, 77.943_sp, 76.821_sp, &
       75.700_sp, 74.578_sp, 73.457_sp, 72.336_sp, 71.214_sp, 70.093_sp, &
       68.971_sp, 67.850_sp, 66.728_sp, 65.607_sp, 64.485_sp, 63.364_sp, &
       62.242_sp, 61.121_sp, 60.000_sp, 58.878_sp, 57.757_sp, 56.635_sp, &
       55.514_sp, 54.392_sp, 53.271_sp, 52.149_sp, 51.028_sp, 49.906_sp, &
       48.785_sp, 47.663_sp, 46.542_sp, 45.420_sp, 44.299_sp, 43.177_sp, &
       42.056_sp, 40.934_sp, 39.813_sp, 38.691_sp, 37.570_sp, 36.448_sp, &
       35.327_sp, 34.205_sp, 33.084_sp, 31.962_sp, 30.841_sp, 29.719_sp, &
       28.598_sp, 27.476_sp, 26.355_sp, 25.234_sp, 24.112_sp, 22.991_sp, &
       21.869_sp, 20.748_sp, 19.626_sp, 18.505_sp, 17.383_sp, 16.262_sp, &
       15.140_sp, 14.019_sp, 12.897_sp, 11.776_sp, 10.654_sp,  9.533_sp, &
        8.411_sp,  7.290_sp,  6.168_sp,  5.047_sp,  3.925_sp,  2.804_sp, &
        1.682_sp,  0.561_sp, -0.561_sp, -1.682_sp, -2.804_sp, -3.925_sp, &
       -5.047_sp, -6.168_sp, -7.290_sp, -8.411_sp, -9.533_sp,-10.654_sp, &
      -11.776_sp,-12.897_sp,-14.019_sp,-15.140_sp,-16.262_sp,-17.383_sp, &
      -18.505_sp,-19.626_sp,-20.748_sp,-21.869_sp,-22.991_sp,-24.112_sp, &
      -25.234_sp,-26.355_sp,-27.476_sp,-28.598_sp,-29.719_sp,-30.841_sp, &
      -31.962_sp,-33.084_sp,-34.205_sp,-35.327_sp,-36.448_sp,-37.570_sp, &
      -38.691_sp,-39.813_sp,-40.934_sp,-42.056_sp,-43.177_sp,-44.299_sp, &
      -45.420_sp,-46.542_sp,-47.663_sp,-48.785_sp,-49.906_sp,-51.028_sp, &
      -52.149_sp,-53.271_sp,-54.392_sp,-55.514_sp,-56.635_sp,-57.757_sp, &
      -58.878_sp,-60.000_sp,-61.121_sp,-62.242_sp,-63.364_sp,-64.485_sp, &
      -65.607_sp,-66.728_sp,-67.850_sp,-68.971_sp,-70.093_sp,-71.214_sp, &
      -72.336_sp,-73.457_sp,-74.578_sp,-75.700_sp,-76.821_sp,-77.943_sp, &
      -79.064_sp,-80.185_sp,-81.307_sp,-82.428_sp,-83.549_sp,-84.670_sp, &
      -85.791_sp,-86.911_sp,-88.029_sp,-89.142_sp/

  INTERFACE
    SUBROUTINE bln2hop(numfi, ie, je, me, ne, alat2, alon2, alat1, alon1, &
     soe, sreg, land)
      USE mo_kind, ONLY: dp, sp
      INTEGER, INTENT(in) :: numfi, ie, je, me, ne
      REAL(sp), INTENT(inout) :: soe(ie+2, je+2, numfi)
      REAL(dp), INTENT(out) :: sreg(me, ne, numfi)
      REAL(sp), INTENT(in) :: alat2(ie+2, je+2), alon2(ie+2, je+2)
      REAL(dp), intent(in) :: alat1(me, ne), alon1(me, ne)
      REAL(sp), INTENT(in) :: land(ie+2, je+2)
    END SUBROUTINE BLN2HOP
  END INTERFACE

  me = -1
  ne = -1
  bounds_exch_tp = .FALSE.
  verbose = 0
  ! setup basic parameters
  READ (*, nml=gridparams)
  IF (me < 1 .OR. ne < 1) THEN
    PRINT *, 'Only positive integral values may be used for me and ne!'
    STOP
  END IF

  ALLOCATE(GILA(2*ME,2*NE),GIPH(2*ME,2*NE))
  ALLOCATE(HOPLAT(ME,NE),HOPLON(ME,NE))
  ALLOCATE(FIOUT(ME,NE,7))
  ALLOCATE(UWOUT(ME,NE), &
       VWOUT(ME,NE), &
       temout(me, ne), &
       ! SSTOUT(ME,NE), &
       pressout(me, ne), &
       tdewout(me, ne), &
       swradout(me, ne), &
       dwlwout(me, ne), &
       prateout(me, ne), &
       u10out(me, ne))
  ALLOCATE(UKE(ME,NE),VKE(ME,NE),UHOPE(ME,NE),VHOPE(ME,NE))

  gausslo(1) = 0._sp
  DO I=2,IE
    gausslo(i) = gausslo(i-1) + 1.125_sp
  ENDDO

  DO I=1,IE
    DO J=1,JE
      GAULAT(I+1, J+1) = GAUSSLA(J)
      GAULON(I+1, J+1) = GAUSSLO(I)
    ENDDO
  ENDDO

  !     ZYKL.RAND OST/WEST
  DO J=1,JE+2
    GAULAT(1,J)=GAULAT(IE+1,J)
    GAULAT(IE+2,J)=GAULAT(2,J)
    gaulon(1, j) = gaulon(ie+1, j) - 360._sp
    gaulon(ie+2, j) = gaulon(2, j) + 360._sp
  ENDDO

  !     POLE (SCHUMMELN: OBEN UND UNTEN EINE REIHE HINZU)
  DO I=1, IE+2

    gaulat(i, 1) = 180._sp - gaulat(i, 2)
    gaulat(i, je+2) = -180._sp - gaulat(i, je+1)
    gaulon(i, 1) = gaulon(i, 2)
    gaulon(i, je+2) = gaulon(i, je+1)

  ENDDO

  !     NOCHMAL ZYKL.RAND OST/WEST
  DO J=2, JE+1
    GAULAT(1, J) = GAULAT(IE+1, J)
    GAULAT(IE+2, J) = GAULAT(2, J)
    gaulon(1, j) = gaulon(ie+1, j) - 360._sp
    gaulon(ie+2, j) = gaulon(2, j) + 360._sp
  ENDDO

  DO I=1, IE+2
    DO J=1, JE+2
      gaulat(i, j) = gaulat(i, j) / REAL(aradtogra, sp)
      gaulon(i, j) = gaulon(i, j) / REAL(aradtogra, sp)
    ENDDO
  ENDDO

  OPEN(io_in_anta, file='anta.ext', access='SEQUENTIAL', &
       form='UNFORMATTED', action='read')
  OPEN(io_out_giwix, file='GIWIX.ext', &
       access='SEQUENTIAL', form='UNFORMATTED', action='write')
  OPEN(io_out_giwiy, file='GIWIY.ext', access='SEQUENTIAL', &
       form='unformatted', action='write')
  OPEN(io_out_gitem, file='GITEM.ext', access='sequential', &
       form='unformatted', action='write')
  OPEN(io_out_gitdew, file='GITDEW.ext', access='sequential', &
       form='unformatted', action='write')
  OPEN(io_out_giprec, file='GIPREC.ext', access='sequential', &
       form='unformatted', action='write')
  !c      OPEN(57,FILE='GIDWLW',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(io_out_gicloud, file='GICLOUD.ext', access='sequential', &
       form='unformatted', action='write')
  OPEN(io_out_giu10, file='GIU10.ext', access='sequential', &
       form='unformatted', action='write')
  OPEN(io_out_swrad, file='GISWRAD.ext', access='sequential', &
       form='unformatted', action='write')
  !c      OPEN(60,FILE='GISST',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  OPEN(io_out_gipress, file='GIPRESS.ext', access='sequential', &
       form='unformatted', action='write')
  OPEN(io_in_lsm, file='land_sea_mask.ECMWF.ext4', access='sequential', &
       form='unformatted', action='read')

  READ (io_in_anta) ext8_hdr
  READ (io_in_anta) GILA
  READ (io_in_anta) ext8_hdr
  READ (io_in_anta) GIPH

  !H    TAGESWERTE
  !cc      READ(*,*)LDAY
  lday=365
  READ(io_in_lsm) ext4_hdr
  READ(io_in_lsm) ((dummy(i,j), i=1, ie), j=1, je)
  FORALL (I = 1:IE, j=1:JE)
    land(i+1,j+1) = MERGE(1._sp, dummy(i, j), dummy(i, j) .GT. 1.e-4_sp)
  END FORALL
  !     ZYKL.RAND OST/WEST/NORDPOL/SUEDPOL
  DO i=1, ie+2
    land(i, 1) = 0._sp
    land(i, je+2) = 1._sp
  END DO
  DO j=1, je+2
    land(1,j) = land(ie+1,j)
    land(ie+2,j) = land(2,j)
  ENDDO

  DO IDAY=1,LDAY
    IF (verbose .GT. 0) PRINT *, 'DATUM ', IDAY
    !H    READ INPUT FILES

    !c      OPEN(2,FILE='sfcfullfilt2416704.no_f_norcw.ext4'
    OPEN(2,FILE='2m_temperature.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(2) ext4_hdr
    READ(2) ((DUMMY(I,J),I=1,IE),J=1,JE)
    FORALL (I =1:IE, J = 1:JE)
      temin(i+1, j+1) = dummy(i,j) - 273.15_sp
    END FORALL

    !c      OPEN(11,FILE='sfcfullfilt2416804.no_f_norcw.ext4'
    OPEN(11,FILE='2m_dewpoint_temperature.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(11) ext4_hdr
    READ(11)((DUMMY(I,J),I=1,IE),J=1,JE)
    FORALL (i = 1:ie, j = 1:je)
      tdewin(i+1, j+1)=dummy(i, j)
    END FORALL

    !c      OPEN(12,FILE='sfcfullfilt2442304.no_f_norcw.ext4'
    OPEN(12,FILE='total_precipitation.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(12) ext4_hdr
    READ(12)((DUMMY(I,J),I=1,IE),J=1,JE)
    FORALL (i = 1:ie, j = 1:je)
      pratein(i+1,j+1) = MAX(0._sp, dummy(i,j))
    END FORALL


    !c      OPEN(13,FILE='sfcfullfilt2417904.ext4'
    !c      OPEN(13,FILE='sfcfullfilt2417904.ext4'
    !c     X    ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    !c      READ(13) ext4_hdr
    !c      READ(13)((DUMMY(I,J),I=1,IE),J=1,JE)
    !c      DO I=1,IE
    !c        DO J=1,JE
    !c          DWLWIN(I+1,J+1)=DUMMY(I,J)
    !c        ENDDO
    !c      ENDDO

    !c      OPEN(13,FILE='sfcfullfilt2417904.ext4'
    OPEN(13,FILE='total_cloud_cover.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(13) ext4_hdr
    READ(13)((DUMMY(I,J),I=1,IE),J=1,JE)
    FORALL (i = 1:ie, j = 1:je)
      dwlwin(i+1, j+1) = dummy(i, j)
    END FORALL


    !c      OPEN(14,FILE='sfcfullfilt2466504.no_f_norcw.ext4'
    OPEN(14,FILE='scalar_wind.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(14) ext4_hdr
    READ(14)((DUMMY(I,J),I=1,IE),J=1,JE)
    FORALL (i = 1:ie, j = 1:je)
      u10in(i+1, j+1) = dummy(i, j)
    END FORALL

    !c      OPEN(15,FILE='sfcfullfilt2427604.no_f_norcw.ext4'
    OPEN(15,FILE='total_solar_radiation.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(15) ext4_hdr
    READ(15)((DUMMY(I,J),I=1,IE),J=1,JE)
    FORALL (i = 1:ie, j = 1:je)
      swradin(i+1, j+1) = dummy(i, j)
    END FORALL

    !c      OPEN(3,FILE='sfcfullfilt2418004.no_f_norcw.ext4'
    OPEN(3,FILE='east_west_stress.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(3) ext4_hdr_temp
    READ(3)((DUMMY(I,J),I=1,IE),J=1,JE)
    FORALL (i = 1:ie, j = 1:je)
      uwin(i+1, j+1) = dummy(i, j)
    END FORALL


    !c      OPEN(4,FILE='sfcfullfilt2418104.no_f_norcw.ext4'
    OPEN(4,FILE='north_south_stress.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(4) ext4_hdr_temp
    READ(4)((DUMMY(I,J),I=1,IE),J=1,JE)
    FORALL (i = 1:ie, j = 1:je)
      vwin(i+1, j+1) = dummy(i, j)
    END FORALL

    !c      OPEN(16,FILE='sfcfullfilt2415104.no_f_norcw.ext4'
    OPEN(16,FILE='mean_sea_level_pressure.ext' &
         ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    READ(16) ext4_hdr_temp
    READ(16)((DUMMY(I,J),I=1,IE),J=1,JE)
    FORALL (i = 1:ie, j = 1:je)
      pressin(i+1, j+1) = dummy(i, j)
    END FORALL


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
      uwin(1, j) = uwin(ie+1, j)
      uwin(ie+2, j) = uwin(2, j)

      vwin(1, j) = vwin(ie+1, j)
      vwin(ie+2, j) = vwin(2, j)

      temin(1, j) = temin(ie+1, j)
      temin(ie+2, J) = temin(2, j)

      pressin(1, j) = pressin(ie+1, j)
      pressin(ie+2, j) = pressin(2, j)

      !        SSTIN(1,J)=SSTIN(IE+1,J)
      !        SSTIN(IE+2,J)=SSTIN(2,J)

      tdewin(1, j) = tdewin(ie+1, j)
      tdewin(ie+2, j) = tdewin(2, j)

      pratein(1, j) = pratein(ie+1, j)
      pratein(ie+2, j) = pratein(2, j)

      dwlwin(1, j) = dwlwin(ie+1, j)
      dwlwin(ie+2, j) = dwlwin(2, j)

      u10in(1, j) = u10in(ie+1, j)
      u10in(ie+2, j) = u10in(2, j)

      swradin(1, j) = swradin(ie+1, j)
      swradin(ie+2, j) = swradin(2, j)

    ENDDO

    !     POLE (EINE REIHE HINZU)
    DO I=1,IE+2

      uwin(i, 1) = uwin(i, 2)
      uwin(i, je+2) = uwin(i, je+1)

      vwin(i, 1) = vwin(i, 2)
      vwin(i, je+2) = vwin(i, je+1)

      temin(i, 1) = temin(i, 2)
      temin(i, je+2) = temin(i, je+1)

      pressin(i, 1) = pressin(i, 2)
      pressin(i, je+2) = pressin(i, je+1)

      !c          SSTIN(I,1)=SSTIN(I,2)
      !c          SSTIN(I,JE+2)=SSTIN(I,JE+1)

      tdewin(i, 1) = tdewin(i, 2)
      tdewin(i, je+2) = tdewin(i, je+1)

      pratein(i, 1) = pratein(i, 2)
      pratein(i, je+2) = pratein(i, je+1)

      dwlwin(i, 1) = dwlwin(i,2)
      dwlwin(i, je+2) = dwlwin(i, je+1)

      u10in(i, 1) = u10in(i, 2)
      u10in(i, je+2) = u10in(i, je+1)

      swradin(i, 1) = swradin(i, 2)
      swradin(i, je+2) = swradin(i, je+1)

    ENDDO

    IF (verbose .GT. 1) PRINT *, 'Datum ', iday, ' EINLESEN UEBERLEBT'

    !H    AN DIESER STELLE WIRD JEWEILS FUER U,V UND SKALARE GROESSEN
    !H    INTERPOLIERT, UND ZWAR AM JEWEILS DAZUGEHOERENDEN GITTERPUNKT
    !H    1. FUER SKALARE WIE GEHABT AN DRUCKPUNKT
    IF (verbose .GT. 1) PRINT*,'SKALARE GROESSE INTERPOLIEREN'

    DO M=1,ME
      DO N=1,NE
        HOPLON(M,N)=GILA(2*M,2*N)
        HOPLAT(M,N)=GIPH(2*M,2*N)
        IF (HOPLON(M,N).GT.(2._wp*API)) HOPLON(M,N)=HOPLON(M,N)-(2._wp*API)
        IF (HOPLON(M,N).LT.0._wp) HOPLON(M,N)=HOPLON(M,N)+(2._wp*API)
      ENDDO
    ENDDO

    !      CPUANF=SECOND()
    DO I=1,IE+2
      DO J=1,JE+2
        fiin(i, j, 1) = temin(i, j)
        fiin(i, j, 2) = tdewin(i, j)
        fiin(i, j, 3) = pratein(i, j)
        fiin(i, j, 4) = dwlwin(i, j)
        fiin(i, j, 5) = u10in(i, j)
        fiin(i, j, 6) = swradin(i, j)
        FIIN(I,J,7) = pressin(i,j)
        !c            FIIN(I,J,8)=SSTIN(I,J)
      ENDDO
    ENDDO
    CALL BLN2HOP(7,IE,JE,ME,NE,GAULAT,GAULON,HOPLAT,HOPLON, &
         FIIN,FIOUT,LAND)

    DO M=1,ME
      DO N=1,NE
        temout(m, n) = fiout(m, n, 1)
        tdewout(m, n) = fiout(m, n, 2)
        prateout(m, n) = fiout(m, n, 3)
        dwlwout(m, n) = fiout(m, n, 4)
        u10out(m, n) = fiout(m, n, 5)
        swradout(m, n) = fiout(m, n, 6)
        pressout(m,n) = fiout(m,n,7)
        !c            SSTOUT(M,N)=FIOUT(M,N,8)

      ENDDO
    ENDDO
    !      CPUEND=SECOND()
    !      PRINT*,'CALL BLN2HOP DAUERT ',CPUEND-CPUANF

    !H    AUSGABE
    !h    PERIODISCHER RAND
    DO N=1,NE

      temout(1, n) = temout(me-1, n)
      temout(me, n) = temout(2, n)

      !c        SSTOUT(1,N)=SStOUT(ME-1,N)
      !c        SSTOUT(ME,N)=SSTOUT(2,N)

      pressout(1, n) = pressout(me-1, n)
      pressout(me, n) = pressout(2, n)

      tdewout(1, n) = tdewout(me-1, n)
      tdewout(me, n) = tdewout(2, n)

      prateout(1, n) = prateout(me-1, n)
      prateout(me, n) = prateout(2, n)

      dwlwout(1, n) = dwlwout(me-1, n)
      dwlwout(me, n) = dwlwout(2, n)

      u10out(1, n) = u10out(me-1, n)
      u10out(me, n) = u10out(2,n)

      swradout(1, n) = swradout(me-1, n)
      swradout(me, n) = swradout(2, n)

    ENDDO







    IF (bounds_exch_tp) THEN

      !          DO i=2,ie-1
      !             il=i
      !             ir=ie+1-i
      ! syncronise line 1 with line 3
      !             a0(il,2) = yrl3(ir)
      ! syncronise line 1 with line 3
      !             a0(il,1) = yrl4(ir)
      !          END DO


      DO m=2,me-1
        ml=i
        mr=me+1-i
        swradout(ml, 2) = swradout(mr, 3)                            ! syncronise line 1 with line 3
        swradout(ml, 1) = swradout(mr, 4)                            ! syncronise line 2 with line 2

        temout(ml, 2) = temout(mr, 3)                            ! syncronise line 1 with line 3
        temout(ml, 1) = temout(mr, 4)                            ! syncronise line 2 with line 2

        tdewout(ml, 2) = tdewout(mr, 3)                            ! syncronise line 1 with line 3
        tdewout(ml, 1) = tdewout(mr, 4)                            ! syncronise line 2 with line 2

        prateout(ml, 2) = prateout(mr, 3)                            ! syncronise line 1 with line 3
        prateout(ml, 1) = prateout(mr, 4)                            ! syncronise line 2 with line 2

        pressout(ml, 2) = pressout(mr, 3)                            ! syncronise line 1 with line 3
        pressout(ml, 1) = pressout(mr, 4)                            ! syncronise line 2 with line 2

        dwlwout(ml, 2) = dwlwout(mr, 3)                            ! syncronise line 1 with line 3
        dwlwout(ml, 1) = dwlwout(mr, 4)                            ! syncronise line 2 with line 2

        u10out(ml, 2) = u10out(mr, 3)                            ! syncronise line 1 with line 3
        u10out(ml, 1) = u10out(mr, 4)                            ! syncronise line 2 with line 2
      END DO
    ENDIF


    ext8_hdr(1) = INT(iday, i8)
    ext8_hdr(3) = INT(ext4_hdr(3), i8)
    ext8_hdr(4) = INT(me * ne, i8)

    ext8_hdr(2) = 167_i8
    WRITE(io_out_gitem) ext8_hdr
    WRITE(io_out_gitem) temout

    ext8_hdr(2) = 168_i8
    WRITE(io_out_gitdew) ext8_hdr
    WRITE(io_out_gitdew) tdewout

    ext8_hdr(2) = 423_i8
    WRITE(io_out_giprec) ext8_hdr
    WRITE(io_out_giprec) prateout

    ext8_hdr(2) = 177_i8
    WRITE(io_out_gicloud) ext8_hdr
    WRITE(io_out_gicloud) dwlwout

    ext8_hdr(2) = 665_i8
    WRITE(io_out_giu10) ext8_hdr
    WRITE(io_out_giu10) u10out

    ext8_hdr(2) = 276_i8
    WRITE(io_out_swrad) ext8_hdr
    WRITE(io_out_swrad) swradout
    !ext8_hdr(2) = 139_i8
    !c      WRITE(60) ext8_hdr
    !c      WRITE(60) SSTOUT
    ext8_hdr(2) = 151_i8
    WRITE(io_out_gipress) ext8_hdr
    WRITE(io_out_gipress) pressout


    !H    2. FUER VECTOR U KOMPONENTE AM U-PUNKT
    IF (verbose .GT. 1) PRINT*,'U-KOMPOMENTE INTERPOLIEREN'

    DO M=2,ME-1
      DO N=1,NE
        HOPLON(M,N)=GILA((2*M)+1,2*N)
        HOPLAT(M,N)=GIPH((2*M)+1,2*N)
        IF (HOPLON(M,N).GT.(2._wp*API)) HOPLON(M,N)=HOPLON(M,N)-(2._wp*API)
        IF (HOPLON(M,N).LT.0._wp) HOPLON(M,N)=HOPLON(M,N)+(2._wp*API)
      ENDDO
    ENDDO

    DO I=1,IE+2
      DO J=1,JE+2
        fiin(i, j, 1) = uwin(i, j)
        fiin(i, j, 2) = vwin(i, j)
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

    IF (verbose .GT. 1) PRINT*,'U-DREHEN'

    DO M=2,ME-1
      DO N=2,NE-1

        IF ((UWOUT(M,N).NE.SVAL).AND.(VWOUT(M,N).NE.SVAL)) THEN

          UWOUT(M,N)=UWOUT(M,N)/1025._wp
          VWOUT(M,N)=VWOUT(M,N)/1025._wp

          DELTXX=(GILA((2*M)+1,2*N)-GILA(2*M,2*N))
          DELTYY= GIPH((2*M)+1,(2*N)-1)-GIPH(((2*M)+1),2*N)
          DELTXY= GIPH((2*M)+1,2*N)-GIPH(2*M,2*N)
          DELTYX=(GILA((2*M)+1,2*N-1)-GILA((2*M)+1,2*N))

          IF(DELTXX.GT. API) DELTXX = DELTXX - 2._wp*API
          IF(DELTXX.LT.-API) DELTXX = DELTXX + 2._wp*API
          IF(DELTYX.LT.-API) DELTYX = DELTYX + 2._wp*API
          IF(DELTYX.GT. API) DELTYX = DELTYX - 2._wp*API

          !H       TEST WINKELTREUE MERKATORPROJECTION

          PROMERU = COS(GIPH(((2*M)+1),2*N))
          DELTXX=DELTXX*PROMERU
          DELTYX=DELTYX*PROMERU


          !H       PROJEKTION VON UKE AUS RICHTUNG I,J IN RICHTUNG N/S, O/W
          !H       MITTELS SKALARPRODUKT = A1*B1+A2*B2+...

          UKE(M,N) = (UWOUT(M,N)*DELTXX) + (VWOUT(M,N)*DELTXY)
          VKE(M,N) = (UWOUT(M,N)*DELTYX) + (VWOUT(M,N)*DELTYY)

          !H       ES MUESSTE EIGENTLICH NOCH DURCH DEN BETRAG
          !H       VON DELTA GETEILT WERDEN DA DELTA KEIN EINHEITSVEKTOR IST

          !H       TEST
          IF ((DELTXX*DELTYX+DELTYY*DELTXY).GT.0.1_wp) THEN
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


    IF (bounds_exch_tp) THEN

      !          DO i=2,ie-1
      !             il=i
      !             ir=ie-i
      ! syncronise line 2 with line 3
      !             a0(il,2) = -yrl3(ir)
      !          END DO


      DO m=2,me-1
        ml=i
        mr=me-i
        UHOPE(ml,2) = -UHOPE(mr,3)                            ! syncronise line 2 with line 3
        UHOPE(ml,1) = -UHOPE(mr,4)                            ! syncronise line 1 with line 4
      END DO
    ENDIF

    ext8_hdr(1) = INT(iday, i8)
    ext8_hdr(2) = 180_i8
    ext8_hdr(3) = INT(ext4_hdr(3), i8)
    ext8_hdr(4) = INT(me * ne, i8)
    WRITE(io_out_giwix) ext8_hdr
    WRITE(io_out_giwix) uhope



    !H    3. FUER VECTOR V KOMPONENTE AM V-PUNKT
    IF (verbose .GT. 1) PRINT *, 'V-KOMPOMENTE INTERPOLIEREN'

    DO M=1,ME
      DO N=2,NE-1
        HOPLON(M,N)=GILA(2*M,(2*N)+1)
        HOPLAT(M,N)=GIPH(2*M,(2*N)+1)
        IF (HOPLON(M,N).GT.(2._wp*API)) HOPLON(M,N)=HOPLON(M,N)-(2._wp*API)
        IF (HOPLON(M,N).LT.0._wp) HOPLON(M,N)=HOPLON(M,N)+(2._wp*API)
      ENDDO
    ENDDO

    FORALL (i = 1:ie+2, j = 1:je+2)
      fiin(i, j, 1) = uwin(i, j)
      fiin(i, j, 2) = vwin(i, j)
    END FORALL

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

    IF (verbose .GT. 1) PRINT *, 'V-DREHEN'

    DO M=2,ME-1
      DO N=2,NE-1

        IF ((UWOUT(M,N).NE.SVAL).AND.(VWOUT(M,N).NE.SVAL)) THEN

          UWOUT(M,N)= UWOUT(M,N)/1025._wp
          VWOUT(M,N)= VWOUT(M,N)/1025._wp

          DELTXX=(GILA(2*M,(2*N)+1)-GILA((2*M)-1,(2*N)+1))
          DELTYY= GIPH(2*M,2*N)    -GIPH(2*M    ,(2*N)+1)
          DELTXY= GIPH(2*M,(2*N)+1)-GIPH((2*M)-1,(2*N)+1)
          DELTYX=(GILA(2*M,2*N)    -GILA(2*M    ,(2*N)+1))

          IF(DELTXX.GT. API) DELTXX = DELTXX - 2._wp*API
          IF(DELTXX.LT.-API) DELTXX = DELTXX + 2._wp*API
          IF(DELTYX.LT.-API) DELTYX = DELTYX + 2._wp*API
          IF(DELTYX.GT. API) DELTYX = DELTYX - 2._wp*API

          !H       TEST WINKELTREUE MERKATORPROJECTION

          PROMERV = COS(GIPH(2*M,(2*N)+1))
          DELTXX=DELTXX*PROMERV
          DELTYX=DELTYX*PROMERV


          !H       PROJEKTION VON VKE AUS RICHTUNG I,J IN RICHTUNG N/S, O/W
          !H       MITTELS SKALARPRODUKT = A1*B1+A2*B2+...

          UKE(M,N) = (UWOUT(M,N)*DELTXX) + (VWOUT(M,N)*DELTXY)
          VKE(M,N) = (UWOUT(M,N)*DELTYX) + (VWOUT(M,N)*DELTYY)

          !H       ES MUESSTE EIGENTLICH NOCH DURCH DEN BETRAG
          !H       VON DELTA GETEILT WERDEN DA DELTA KEIN EINHEITSVEKTOR IST

          !H       TEST
          IF ((DELTXX*DELTYX+DELTYY*DELTXY).GT.0.1_wp) THEN
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

    IF (bounds_exch_tp) THEN
      !       IF (ttt=='v') THEN
      ! v point with sign change
      !          DO i=2,ie-1
      !             il=i
      !             ir=ie+1-i
      !             pq(1:kk)=0.5_wp*(arr(il,2,1:kk)-yrl2(ir,1:kk))
      !             arr(il,2,1:kk)= pq(1:kk)
      ! syncronise line 2 with line 2
      !             arr(ir,2,1:kk)=-pq(1:kk)
      !          arr(il,1,1:kk)=-yrl3(ir,1:kk)
      !          END DO
      !       ENDIF

      DO m=2,me-1
        ml=i
        mr=me+1-i
        VHOPE(ml,1) = -VHOPE(mr,3)                            ! syncronise line 1 with line 3
        pq=0.5_wp*(VHOPE(ml,2)-VHOPE(mr,2))
        vhope(ml,2)= pq                       ! syncronise line 2 with line 2
        vhope(mr,2)=-pq
        !             VHOPE(ml,2) = -VHOPE(mr,2)                            ! syncronise line 2 with line 2
      END DO
    END IF


    ext8_hdr(1) = INT(iday, i8)
    ext8_hdr(2) = 181_i8
    ext8_hdr(3) = INT(ext4_hdr(3), i8)
    ext8_hdr(4) = INT(me * ne, i8)
    WRITE(io_out_giwiy) ext8_hdr
    WRITE(io_out_giwiy) vhope


  enddo

  stop
end PROGRAM FORCING

!H    ********************************************************

SUBROUTINE BLN2HOP(NUMFI,IE,JE,ME,NE,ALAT2,ALON2,ALAT1,ALON1, &
     SOE,SREG,LAND)

  !     BILINEAR INTERPOLATION FROM ONE GRID(IE,JE) TO ANOTHER GRID(ME,NE)

  !      PARAMETER(NB=3,MB=3)
  USE mo_kind, ONLY: dp, sp, wp
  USE mo_constants, ONLY: api!, aradtogra
  IMPLICIT NONE
  INTEGER, INTENT(in) :: numfi, ie, je, me, ne
  REAL(sp), INTENT(inout) :: soe(ie+2, je+2, numfi)
  REAL(dp), INTENT(out) :: sreg(me, ne, numfi)
  REAL(sp), INTENT(in) :: alat2(ie+2, je+2), alon2(ie+2, je+2)
  REAL(dp), INTENT(in) :: alat1(me, ne), alon1(me, ne)
  REAL(sp), INTENT(in) :: land(ie+2, je+2)
  REAL(sp), PARAMETER :: epsilon=1.e-12_sp
  REAL(wp) G(IE+2,JE+2), &
       alpha, beta, wwwalpha, wwwbeta, &
       HX(IE+2,JE+2,NUMFI),HHX(IE+2,JE+2,NUMFI), &
       HG(IE+2,JE+2), rsumg
  INTEGER :: i, il, ilu, ir, iter, j, jlu, jo, ju, l, m, n, nland

  !     DIFFUSION INTO LAND (NEW 12/99)

  DO L=1,NUMFI
    DO J=1,JE+2
      DO I=1,IE+2
        HHX(I,J,L)=SOE(I,J,L)
        !HH      X                 *((LAND(I,J)-1.))*(-1.)
      ENDDO
    ENDDO
  ENDDO

  DO J=1,JE+2
    DO I=1,IE+2
      hg(i, j) = MERGE(REAL(epsilon, wp), 1._wp, land(i, j) .GE. 0.5_sp)
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
        IF(IL.LT.1)IL=IE
        IR=I+1
        IF(IR.GT.IE+2)IR=3
        RSUMG=0._wp
        IF (LAND(I,J) .GE. 0.5_sp)THEN
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
    nland=0
    Do I=1,IE
      DO j=1,je
        IF (HG(i,j) .LE. 2._wp * REAL(epsilon, dp)) nland=nland+1
      ENDDO
    ENDDO
    !PRINT*,iter,' Nland: ', nland
  ENDDO


  DO J=1,JE+2
    DO I=1,IE+2
      G(i,j)=HG(i,j)
      DO L=1,NUMFI
        SOE(I,J,L)=HHX(I,J,L)
        IF (ABS(SOE(I,J,L)).LT. epsilon) SOE(I,J,L)=0._wp
      ENDDO
    ENDDO
  ENDDO

  !
  DO M=1,ME
    DO N=1,NE

      !        PUNKT RECHTS OBEN
      DO J=1,JE+2
        IF (REAL(alat2(2, j), dp) .GE. alat1(m, n)) jo=j
      ENDDO
      DO I=1,IE+2
        IF (REAL(alon2(i,2), dp) .LE. alon1(m,n)) il=i
      ENDDO

      !       IL=NINT(ALON1(M,N)*DLAMDAI+0.499999999)+1

      !            PRINT*,'WE ARE AT THE MINIMUM ',M,N,ALAT1(M,N)*aradtogra
      !     X                                    ,IL,JO,ALAT2(IL,JO)*aradtogra

      !       PUNKT RECHTS OBEN --> LINKS UNTEN

      JLU=JO+1
      ILU=IL

      WWWALPHA=ALAT2(ILU,JLU)-ALAT1(M,N)
      WWWBETA=ALON2(ILU,JLU)-ALON1(M,N)
      IF(WWWBETA.GE.API) WWWBETA=WWWBETA-2._wp*API
      IF(WWWBETA.LE.-API) WWWBETA=WWWBETA+2._wp*API

      ALPHA=(WWWALPHA)/(ALAT2(ILU,JLU)-ALAT2(ILU,JLU-1))
      BETA=(WWWBETA)/(ALON2(ILU,JLU)-ALON2(ILU+1,JLU))
      DO I=1,NUMFI
        SREG(M,N,I)=ALPHA*BETA*SOE(ILU+1,JLU-1,I)*G(ILU+1,JLU-1)   &
             +(1._wp-ALPHA)*(1._wp-BETA)*SOE(ILU,JLU,I)*G(ILU,JLU) &
             +(1._wp-ALPHA)*(BETA)*SOE(ILU+1,JLU,I)*G(ILU+1,JLU) &
             +(ALPHA)*(1._wp-BETA)*SOE(ILU,JLU-1,I)*G(ILU,JLU-1)
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


! activate once the previous code has been confirmed to work
SUBROUTINE grid_add_periodic_rim_dp(ie, je, a, x_overlap, y_overlap, &
     bounds_exch_tp)
  USE mo_kind, ONLY: dp
  INTEGER, INTENT(in) :: ie, je
  REAL(dp), INTENT(inout) :: a(ie, je)
  INTEGER, INTENT(in) :: x_overlap, y_overlap
  LOGICAL, INTENT(in) :: bounds_exch_tp

  INTEGER :: j, i!, il, ir
  FORALL (j = 1:je, i = 1:x_overlap)
    a(i, j) = a(ie - i, j)
    a(ie - i + 1, j) = a(i + 1, j)
  END FORALL

  IF (bounds_exch_tp) THEN
    FORALL (i = 2:ie - 1,  j = 1:y_overlap)
      ! il = i
      ! ir = ie + 1 - i
      ! syncronise line j with line 2*overlap - j + 1
      a(i, j) = a(ie + 1 - i, 2*y_overlap - j + 1)
    END FORALL
  END IF
END SUBROUTINE grid_add_periodic_rim_dp
