MODULE cgri_shared
  USE mo_kind, ONLY: i4, wp
  INTEGER(kind=i4), PARAMETER :: resolution(2) = (/ 4320, 2160 /)
  REAL(wp), ALLOCATABLE :: feld(:, :)
CONTAINS
  SUBROUTINE cgri_init
    ALLOCATE(feld(resolution(1), resolution(2)))
  END SUBROUTINE cgri_init
END MODULE cgri_shared

PROGRAM cgri
  USE mo_kind
  USE mo_constants, ONLY: api, aradtogra
  USE cgri_shared
  USE iso_varying_string
  IMPLICIT NONE
  REAL(wp), PARAMETER :: reflat=55._wp, reflon=1.45_wp, topfai=resolution(2)/api
  INTEGER, PARAMETER :: trig_tab_res = 629
  INTEGER :: NP, np1, ie, ie1, je, it8, ito, jf, jto, jep, je6
  INTEGER :: i, iant, idim, iet5, ii, ipr1, ipr2, ipr3, &
       irad, irr, j, jant, jd, jet5, jfi, jj, jmid, l, m, nde, nel
  INTEGER :: idr, kitt, ilim!, is
  LOGICAL ifind

  REAL(wp) x, y, xx, yy, a, aa, absch, al, arg, auft, b, bb, be, betra
  REAL(wp) dphih, elp2, elp1, criac, de, deep, dfdr, dphi, dra, dx, dy, &
       elev, f, ga, gril, grill, grilr, griph, phiq, poldi, r, ra, &
       radi, ralt, reff, rr, sipia, tw, zerr, ziel

  REAL(wp), ALLOCATABLE :: dphieq(:), x8(:, :), y8(:, :)
  ! REAL(wp), DIMENSION(trig_tab_res) :: sinat, cosat
  REAL(wp), ALLOCATABLE :: elrad(:)
  REAL(wp), ALLOCATABLE :: giph(:, :), gila(:, :), giltes(:, :), xk(:, :), &
       yk(:, :)
  REAL(wp), ALLOCATABLE :: xh(:, :), yh(:, :)

  INTEGER(kind=i8) :: date, ext8_hdr(4)
  INTEGER :: ioerrstat, verbose
  INTEGER, PARAMETER :: totato_lun=12, anta_lun=13, topo_lun=33
  REAL(wp), ALLOCATABLE :: depto(:, :)
  ! file names
  TYPE(varying_string) :: totato_fname, anta_fname, topo_fname

  NAMELIST /gridparams/ np, jf, date, verbose
  np = -1
  jf = -1
  date = -1_i8
  verbose = 0

  CALL cgri_init

  ! setup basic parameters
  READ (*, nml = gridparams)
  IF (verbose .GT. 0) PRINT *, 'np=', np, 'jf=', jf, 'date=', date
  IF (np < 1 .OR. jf < 1 .OR. date < 1_i8) THEN
    PRINT *, 'Only positive integral values may be used for np, jf and date!'
    STOP 1
  END IF

  CALL get(totato_fname, iostat=ioerrstat)
  IF (ioerrstat .GT. 0 .OR. CHAR(totato_fname) == '') totato_fname='TOTATO'

  CALL get(anta_fname, iostat=ioerrstat)
  IF (ioerrstat .GT. 0 .OR. CHAR(anta_fname) == '') anta_fname='antar8'

  CALL get(topo_fname, iostat=ioerrstat)
  IF (ioerrstat .GT. 0 .OR. CHAR(topo_fname) == '') topo_fname='topo'


  np1 = np + 1
  ie = 4 * np
  ie1 = ie + 2
  je = 2 * np
  it8 = 2 * ie1
  ito = 2 * ie1
  jto = 2 * jf
  jep = je + 2
  je6 = je + 6

  ALLOCATE(dphieq(ito), x8(ito, jto), y8(ito, jto), elrad(jep), &
       giph(ito, jto), gila(ito, jto), giltes(ito, 3),  depto(ie1, jf), &
       xk(ito, jto), yk(ito, jto), xh(ito, jto), yh(ito, jto))

  ext8_hdr(1) = date
  ext8_hdr(2) = INT(np, i8)
  ext8_hdr(3) = 1_i8
  ext8_hdr(4) = INT(ito * jto, i8)

  IF (verbose .GT. 0) PRINT *, ie1, it8, jto, it8 * jto

  OPEN(totato_lun, file=CHAR(totato_fname), form='formatted', action='read')
  OPEN(anta_lun, file=CHAR(anta_fname), form='unformatted', action='write')
  READ(totato_lun, 1200) feld

1200 FORMAT(12f6.0)

  dphih = 0.25_wp * api / REAL(np, wp)
  elp2 = 2._wp / (REAL(np, wp) + 0.5_wp)
  elp1 = 1._wp / (REAL(np, wp) + 0.5_wp)

  IF (verbose .GT. 0) &
       PRINT *, 'latitude of singularity', reflat, reflon * aradtogra
  poldi = 90._wp - reflat
  absch = 2._wp * TAN(poldi * api / 360._wp)
  IF (verbose .GT. 0) PRINT *, 'absch=', absch
  aa = 2._wp - absch
  bb = absch
!   DO is = 1, trig_tab_res
!     arg = 0.01_wp * REAL(is, wp)
!     sinat(is) = SIN(arg)
!     cosat(is) = COS(arg)
!   END DO
  DO irr = 1, jep
    rr = 0.5_wp * REAL(irr - 1, wp)
    criac = 1._wp - 1.e-10_wp
    ra = MIN(elp2 * rr / (1._wp + elp1 * rr), criac)

    elrad(irr) = ra
    a = 2._wp - ra * aa
    b = MAX(2._wp - ra * aa - bb * ra**2, 1.e-9_wp)
    IF (verbose .GT. 1) PRINT *, 'axes', irr, ra, a, b
!     ip = 3
!     DO is = 1, trig_tab_res
!       xx = a * cosat(is)
!       yy = b * sinat(is)
!       CALL plot(xx, yy, ip)
!       ip = 2
!     END DO
  END DO
  IF (verbose .GT. 0) PRINT * , elrad
  auft = api / (4._wp * REAL(np, wp))
  !     construction of orthogonals to ellipses
  DO irad = 1, ito
    !     do irad = 35, 35
    iant = irad
    arg = REAL(irad - 3, wp) * auft
    x = 2._wp * COS(arg)
    y = 2._wp * SIN(arg)
    xx = x
    yy = y
    xk(iant, jep) = xx
    yk(iant, jep) = yy
    x8(iant, jep) = x
    y8(iant, jep) = y
    gila(iant, jep) = arg
    IF (verbose .GT. 1) PRINT * , 'start', irad, arg, x, y
    zerr = 1.05_wp + 0.002_wp * REAL(irad, wp)
    !     call plot(xx, yy, 3)
    ra = 0._wp

    tw = 2._wp
    loop200: DO l = 1, je
      jant = jep - l
      ziel = elrad(l + 1)

      ifind = .FALSE.
      dra = 0.000001_wp
      IF(l .GT. 5) dra = 0.00001_wp
      loop101: DO idr = 1, 100000
        reff = 1._wp
        ralt = ra

        DO kitt = 1, 25

          f =    x**2 / (tw - ra * aa)**2 &
               + y**2 / (tw - ra * aa - bb * ra**2)**2 - 1._wp
          dfdr =   -tw * x**2 / (tw - ra * aa)**3 &
                -   tw * y**2 * (aa + tw * bb * ra) &
                       /  (tw - ra * aa - bb * ra**2)**3
          ra = ra + f / dfdr
        END DO

        IF (ra .GE. ziel) THEN
          ifind = .TRUE.
          reff = (ziel - ralt) / (ra - ralt)
          ra = ziel
        ENDIF
        dx = 2._wp * x / (2._wp - aa * ra)**2
        dy = 2._wp * y / (2._wp - aa * ra - bb * ra**2) ** 2
        !     print *, 'dx', dx, dy
        betra = SQRT(dx**2 + dy**2)
        r = 1._wp
        dx = dra * dx / betra
        dy = dra * dy / betra
        IF(dy * y .GT. 0.9_wp * y**2) r = 0.9_wp * y / dy
        x = x - dx * r * reff
        y = y - dy * r * reff
        xx = x
        yy = y
        IF (ifind) EXIT loop101
      END DO loop101
      IF (l .EQ. je + 1) PRINT *, '?', irad, x, y
      IF (ABS(y) .LT. 1.e-6_wp) y = 0._wp

      xk(iant, jant) = x
      yk(iant, jant) = y
      x8(iant, jant) = x
      y8(iant, jant) = y
      IF (verbose .GT. 1 .AND. ( irad .EQ. 3 .OR. irad .EQ. 403)) &
           PRINT *, 'ax2', irad, x, y

      gila(iant, jant) = ATAN2(y, x)
      radi = SQRT(x**2 + y**2)
      poldi = 2._wp * ATAN(0.5_wp * radi)
      giph(iant, jant) = 0.5_wp * api - poldi

    END DO loop200
  ENDDO
  IF (verbose .GT. 0) PRINT *, 'zum ersten'
  !      do j = 1, jep
  !         print *, 'xk', j, yk(44, j), xk(44, j)**2 + yk(44, j)**2
  !      enddo
  !      do i = 1, 50
  !         print *, i, xk(i, jep), yk(i, jep), xk(i, jep)**2 + yk(i, jep)**2
  !      enddo
  !      print *, 'yk'
  !      print *, (i, yk(i, 42), i = 1, ito)

  DO i = 1, ito
    xk(i, 1) = xk(i, 2)
    yk(i, 1) = 0._wp
    x8(i, 1) = x8(i, 2)
    y8(i, 1) = 0._wp
  ENDDO

  DO j = 1, jto
    DO i = 1, ito
      xh(i, j) = x8(i, j)
      yh(i, j) = y8(i, j)
    ENDDO
  ENDDO
  sipia = 2._wp * SIN(0.15_wp * api)
  ipr1 = 3
  ipr2 = 4 * np + 3
  ipr3 = 8 * np + 3
  idim = INT(0.6_wp * REAL(np, wp))
  DO i = 1, ito
    !     if(abs(yk(i, je)) .lt. sipia) then
    jmid = MIN(ABS(i - ipr1), ABS(i - ipr2), ABS(i - ipr3))
    IF(jmid .LE. idim)THEN
      jfi = MIN(np / 5, idim - jmid)
      jd = 2 * jfi + 1
      IF (verbose .GT. 1) PRINT * , 'jd0:', jfi, np / 5, idim, jmid, idim - jmid
      DO j = 0, jfi - 1
        IF (verbose .GT. 2) PRINT *, 'jd1:', jd - 2 * j - 1
        xk(i, jd - 2 * j - 1) =   0.5_wp * (xh(i, jd - j) &
                                      + xh(i, jd - (j + 1)))
        yk(i, jd - 2 * j - 1) =   0.5_wp * (yh(i, jd - j) &
                                      + yh(i, jd - (j + 1)))
        IF (verbose .GT. 2) PRINT *, 'jd2:', jd - 2 * j - 2
        xk(i, jd - 2 * j - 2) = xh(i, jd - j - 1)
        yk(i, jd - 2 * j - 2) = yh(i, jd - j - 1)
        x8(i, jd - 2 * j - 1) =   0.5_wp * (xh(i, jd - j) &
                                      + xh(i, jd - (j + 1)))
        y8(i, jd - 2 * j - 1) =   0.5_wp * (yh(i, jd - j) &
                                      + yh(i, jd - (j + 1)))
        x8(i, jd - 2 * j - 2) = xh(i, jd - j - 1)
        y8(i, jd - 2 * j - 2) = yh(i, jd - j - 1)
      ENDDO
    ENDIF
  ENDDO

  DO j = jep, 1, - 1
    DO i = 1, ito
      x = xk(i, j)
      y = yk(i, j)
      xk(i, j + 4) = x
      yk(i, j + 4) = y
      gila(i, j + 4) = ATAN2(y, x)
      radi = SQRT(x**2 + y**2)
      poldi = 2._wp * ATAN(0.5_wp * radi)
      giph(i, j + 4) = 0.5_wp * api - poldi
    ENDDO
  ENDDO
  IF (verbose .GT. 0) PRINT *, 'gila'



622 FORMAT(6e10.3)


  DO j = 1, 4
    DO i = 1, ito
      yk(i, 5 - j) = 2._wp * yk(i, 5) - yk(i, 5 + j)
      xk(i, 5 - j) = 2._wp * xk(i, 5) - xk(i, 5 + j)
      !      y8(i, 5 - j) = 2._wp * y8(i, 5) - y8(i, 5 + j)
      !      x8(i, 5 - j) = 2._wp * x8(i, 5) - x8(i, 5 + j)
    ENDDO
  ENDDO

  !!      do j = 1, 5
  DO j = 1, 4
    DO i = 1, ito
      gila(i, j) = ATAN2(yk(i, j), xk(i, j))
      radi = SQRT(xk(i, j)**2 + yk(i, j)**2)
      poldi = 2._wp * ATAN(0.5_wp * radi)
      giph(i, j) = 0.5_wp * api - poldi
    ENDDO
  ENDDO


601 FORMAT('wo', i3, 7e10.3)
  dphi = api / REAL(ie, wp)

  phiq = 0._wp
  DO i = 1, ito
    dphieq(i) = giph(i, je6 - 1) - giph(i, je6)
    phiq = phiq + dphieq(i)
  ENDDO
  phiq = phiq / REAL(ito, wp)
  m = 0
  DO j = je6, jto - 1
    m = m + 1
    DO i = 1, ito
      gila(i, j + 1) = gila(i, j)
      dphi = (40._wp * dphieq(i) + phiq * REAL(m, wp)) / REAL(m + 40, wp)
      giph(i, j + 1) = giph(i, j) - dphi
    ENDDO
  ENDDO

  IF (verbose .GT. 0) PRINT *, 'final rotation'
  DO j = 1, jto
    DO i = 1, ito
      gila(i, j) = MOD(gila(i, j) + reflon + api, 2._wp * api) - api
    ENDDO
  ENDDO

  !HH   Here we need to set the longitudes of the northpole grid point
  !     to avoid problems with remapcon interpolation

  gila(2 * np + 3, 5) = gila(2 * np + 3, 6)
  gila(6 * np + 3, 5) = gila(6 * np + 3, 6)


  IF (verbose .GT. 0) PRINT *, 'np1-1', np1 - 1, gila(np1 - 1, 5) * aradtogra, &
       gila(np1 - 1, 6) * aradtogra
  IF (verbose .GT. 0) PRINT *, 'np1  ', np1, gila(np1, 5) * aradtogra, &
       gila(np1, 6) * aradtogra
  IF (verbose .GT. 0) PRINT *, 'np1+1', np1 + 2, gila(np1 + 2, 5) * aradtogra, &
       GILA(np1 + 2, 6) * aradtogra

  ext8_hdr(2) = 54_i8
  WRITE (anta_lun) ext8_hdr
  WRITE (anta_lun) gila
  ext8_hdr(2) = 55_i8
  WRITE (anta_lun) ext8_hdr
  WRITE (anta_lun) giph

  IF (verbose .GT. 0) PRINT *, 'suchen'

  DO j = 2, jf - 1
    DO i = 2, ie1 - 1

      elev = 0._wp
      deep = 0._wp
      nel = 0
      nde = 0

      DO jj = 1, 10
        DO ii = 1, 10
          al = 0.1_wp * REAL(ii, wp) - 0.05_wp
          ga = 0.1_wp * REAL(jj, wp) - 0.05_wp
          be = 1._wp - al
          de = 1._wp - ga
          gril = topfai &
               * (al * gila(2 * i + 1, 2 * j) &
                  + be * gila(2 * i - 1, 2 * j))
          grilr = gila(2 * i + 1, 2 * j)
          grill = gila(2 * i - 1, 2 * j)
          IF (ABS(grilr - grill) .GT. 2._wp) grill = grilr
          gril = (al * grilr + be * grill) * topfai

          griph = topfai &
               * (ga * giph(2 * i, 2 * j + 1) &
                  + de * giph(2 * i, 2 * j - 1))
          iet5 = INT(MOD(REAL(gril, wp) + 4319.0_wp, REAL(resolution(1), wp)) + 1.0_wp)
          jet5 = MIN(1080 - INT(griph), resolution(2))
          !                  print * , 'jet5', jet5, 1080, griph

          IF (feld(iet5, jet5) .GT. 0._wp) elev = elev + feld(iet5, jet5)
          IF (feld(iet5, jet5) .LE. 0._wp) deep = deep - feld(iet5, jet5)
          IF (feld(iet5, jet5) .GT. 0._wp) nel = nel + 1
          IF (feld(iet5, jet5) .LE. 0._wp) nde = nde + 1

        ENDDO
      ENDDO

      IF(nde .GE. nel .AND. nde .NE. 0) depto(i, j) = deep / REAL(nde, wp)
      IF(nel .GT. nde) depto(i, j) = 0._wp
    ENDDO
  ENDDO

  IF (verbose .GT. 0) PRINT *, 'nach besetzung'
  ext8_hdr(4) = INT(ie1 * jf, i8)
  ext8_hdr(2) = 84_i8
  WRITE (anta_lun) ext8_hdr
  WRITE (anta_lun) depto
  IF (verbose .GT. 0) PRINT *, 'rundplot'
  OPEN(topo_lun, file = CHAR(topo_fname), form='formatted')
  DO ii = 1, ie, 20
    WRITE(topo_lun, *)'streifen', ii
    ilim = MIN(20, ie - ii + 1)
    DO j = 1, jf
      WRITE(topo_lun, 3300)j, (depto(ii + i, j), i = 1, ilim)
    ENDDO
  ENDDO
3300 FORMAT(i5, 20f6.0)

  STOP
END PROGRAM cgri
