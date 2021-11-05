PROGRAM SAM
  USE mo_planetary_constants, ONLY: radius, omega
  USE mo_constants, ONLY: aradtogra!, api
  USE mo_kind, ONLY: i8, wp
  IMPLICIT NONE
  INTEGER :: i, j, ie, je, ito, jto, ie1, je1, mende
  INTEGER(i8) :: ext_hdr(4)
  REAL(WP) :: colaee, colane, colase, colass, colasw, colaum, colaup, colavm, &
       colavp, colazz, cophee, cophne, cophse, cophss, cophsw, cophum, &
       cophup, cophvm, cophvp, cophzz, rmindlxp, rmindlyp, silaee, silane, &
       silase, silass, silasw, silaum, silaup, silavm, silavp, silazz, &
       siphee, siphne, siphse, siphsw, siphss, siphum, siphup, siphvm, &
       siphvp, siphzz, xee, xne, xse, xss, xsw, xum, xup, xvm, xvp, &
       xzz, yee, yne, yse, yss, ysw, yum, yup, yvm, yvp, yzz, zee, &
       zne, zse, zss, zsw, zum, zup, zvm, zvp, zzz
  REAL(WP), ALLOCATABLE :: depto(:, :)
  !REAL(WP), ALLOCATABLE :: MASKE(:, :), MASKO(:, :)
  REAL(WP), ALLOCATABLE :: dlxp(:, :), dlyp(:, :), dlxv(:, :), dlyu(:, :) &
       , dlxu(:, :), dlyv(:, :)
  REAL(WP), ALLOCATABLE :: ftwou(:, :), ftwov(:, :)
  REAL(WP), ALLOCATABLE :: GILA(:, :), GIPH(:, :)
  INTEGER, PARAMETER :: io_out_arcgri=82, io_in_anta=81
  INTEGER :: verbose
  NAMELIST /gridparams/ ie, je, verbose

  ! setup default values
  ie = -1
  je = -1
  verbose = 0
  READ (*, nml=gridparams)
  IF (verbose .GT. 0) PRINT *, 'ie=', ie, 'je=', je
  IF (ie < 1 .OR. je < 1) THEN
    PRINT *, 'Only positive integral values may be used for ie and je!'
    STOP
  END IF

  ito = 2 * ie
  jto = 2 * je
  je1 = je - 1
  ie1 = ie - 1

  !  ALLOCATE(MASKE(IE, JE), MASKO(IE, JE))
  ALLOCATE(dlxp(ie, je), dlyp(ie, je), dlxv(ie, je), dlyu(ie, je), &
       dlxu(ie, je), dlyv(ie, je), ftwou(ie, je), ftwov(ie, je))
  ALLOCATE(GILA(ITO, JTO), GIPH(ITO, JTO))
  ALLOCATE(DEPTO(IE, JE))

  OPEN(io_in_anta, FILE='anta', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
  OPEN(io_out_arcgri, FILE='arcgri', ACCESS='SEQUENTIAL', FORM='UNFORMATTED', &
       action='write')

  IF (verbose .GT. 0) PRINT*, ' vor gila '
  READ(io_in_anta) ext_hdr
  READ(io_in_anta) gila

  IF (verbose .GT. 0) PRINT*, ' vor giph'
  read(io_in_anta) ext_hdr
  read(io_in_anta) giph


  IF (verbose .GT. 0) PRINT*, ' vor depto '
  read(io_in_anta) ext_hdr
  read(io_in_anta) depto


  !      do 7722 j=1, je
  !      do 7722 i = 1, ie
  !      masko(i, j) = 0
  !      if(depto(i, j) .gt. 1.) masko(i, j) = 1
  !7722  continue
  !      do 5867 j = 1, jto
  !      do 5867 i = 2, ito
  !      dellg = gila(i, j) - gila(i - 1, j)
  !      if(abs(dellg) .lt. pi) go to 5867
  !      if(dellg .lt.  - pi) then
  !      do 5868 ii = i, ito
  !       gila(ii, j) = gila(ii, j) + (2.*api)
  !
  !5868   continue
  !      else
  !      do 5869 ii = i, ito
  !      gila(ii, j) = gila(ii, j) - (2.*api)
  !5869  continue
  !      endif
  !5867  continue


  !      mende = 0
  !65543 format(20f6.0)

  DO j = 1, je1
    DO i = 1, ie1



      siphup = sin(giph(2 * i + 1, 2 * j))
      cophup = cos(giph(2 * i + 1, 2 * j))
      siphum = sin(giph(2 * i - 1, 2 * j))
      cophum = cos(giph(2 * i - 1, 2 * j))
      siphvp = sin(giph(2 * i, 2 * j + 1))
      cophvp = cos(giph(2 * i, 2 * j + 1))
      siphvm = sin(giph(2 * i, 2 * j - 1))
      cophvm = cos(giph(2 * i, 2 * j - 1))

      cophzz = cos(giph(2 * i, 2 * j))
      cophsw = cos(giph(2 * i - 1, 2 * j + 1))
      cophss = cos(giph(2 * i, 2 * j + 2))
      cophse = cos(giph(2 * i + 1, 2 * j + 1))
      cophne = cos(giph(2 * i + 1, 2 * j - 1))
      cophee = cos(giph(2 * i + 2, 2 * j))

      siphzz = sin(giph(2 * i, 2 * j))
      siphsw = sin(giph(2 * i - 1, 2 * j + 1))
      siphss = sin(giph(2 * i, 2 * j + 2))
      siphse = sin(giph(2 * i + 1, 2 * j + 1))
      siphne = sin(giph(2 * i + 1, 2 * j - 1))
      siphee = sin(giph(2 * i + 2, 2 * j))

      silaup = sin(gila(2 * i + 1, 2 * j))
      colaup = cos(gila(2 * i + 1, 2 * j))
      silaum = sin(gila(2 * i - 1, 2 * j))
      colaum = cos(gila(2 * i - 1, 2 * j))
      silavp = sin(gila(2 * i, 2 * j + 1))
      colavp = cos(gila(2 * i, 2 * j + 1))
      silavm = sin(gila(2 * i, 2 * j - 1))
      colavm = cos(gila(2 * i, 2 * j - 1))


      colazz = cos(gila(2 * i, 2 * j))
      colasw = cos(gila(2 * i - 1, 2 * j + 1))
      colass = cos(gila(2 * i, 2 * j + 2))
      colase = cos(gila(2 * i + 1, 2 * j + 1))
      colane = cos(gila(2 * i + 1, 2 * j - 1))
      colaee = cos(gila(2 * i + 2, 2 * j))

      silazz = sin(gila(2 * i, 2 * j))
      silasw = sin(gila(2 * i - 1, 2 * j + 1))
      silass = sin(gila(2 * i, 2 * j + 2))
      silase = sin(gila(2 * i + 1, 2 * j + 1))
      silane = sin(gila(2 * i + 1, 2 * j - 1))
      silaee = SIN(gila(2 * i + 2, 2 * j))
      xup = silaup * cophup
      yup = colaup * cophup
      zup = siphup

      xum = silaum * cophum
      yum = colaum * cophum
      zum = siphum

      xvp = silavp * cophvp
      yvp = colavp * cophvp
      zvp = siphvp

      xvm = silavm * cophvm
      yvm = colavm * cophvm
      zvm = siphvm


      xzz = silazz * cophzz
      yzz = colazz * cophzz
      zzz = siphzz


      xee = silaee * cophee
      yee = colaee * cophee
      zee = siphee


      xss = silass * cophss
      yss = colass * cophss
      zss = siphss


      xsw = silasw * cophsw
      ysw = colasw * cophsw
      zsw = siphsw


      xse = silase * cophse
      yse = colase * cophse
      zse = siphse

      xne = silane * cophne
      yne = colane * cophne
      zne = siphne

      IF (verbose .GT. 1) PRINT *, i, j

      !      dlyu(i, j) = radius * acos(xne * xse + yne * yse + zne * zse)
      !      dlxv(i, j) = radius * acos(xsw * xse + ysw * yse + zsw * zse)
      !      dlxu(i, j) = radius * acos(xzz * xee + yzz * yee + zzz * zee)
      !      dlyv(i, j) = radius * acos(xzz * xss + yzz * yss + zzz * zss)
      !      dlxp(i, j) = radius * acos(xup * xum + yup * yum + zup * zum)
      !      dlyp(i, j) = radius * acos(xvp * xvm + yvp * yvm + zvp * zvm)


      dlyu(i, j) = MAX(1._wp, radius * ACOS(MIN((xne * xse + yne * yse + zne * zse), 1._wp)))
      dlxv(i, j) = MAX(1._wp, radius * ACOS(MIN((xsw * xse + ysw * yse + zsw * zse), 1._wp)))
      dlxu(i, j) = MAX(1._wp, radius * ACOS(MIN((xzz * xee + yzz * yee + zzz * zee), 1._wp)))
      dlyv(i, j) = MAX(1._wp, radius * ACOS(MIN((xzz * xss + yzz * yss + zzz * zss), 1._wp)))
      dlxp(i, j) = MAX(1._wp, radius * ACOS(MIN((xup * xum + yup * yum + zup * zum), 1._wp)))
      dlyp(i, j) = MAX(1._wp, radius * ACOS(MIN((xvp * xvm + yvp * yvm + zvp * zvm), 1._wp)))


      IF (verbose .GT. 1) PRINT *, i, j

      IF (dlxv(i, j) .LT. 1._wp) WRITE(6, *)'scheiss dlxv:', i, j, dlxv(i, j), &
           dlyv(i, j), dlxp(i, j), dlxu(i, j), gila(2 * i, 2 * j), giph(2 * i, 2 * j)

      IF (verbose .GT. 0 .AND. i .EQ. 25) THEN
        PRINT *, 'pos', j, gila(2 * i, 2 * j) * aradtogra, &
             giph(2 * i, 2 * j) * aradtogra, &
             dlxp(i, j), dlyp(i, j)
        PRINT *, aradtogra * gila(2 * i, 2 * j - 1), &
             aradtogra * gila(2 * i, 2 * j + 1), &
             aradtogra * giph(2 * i, 2 * j - 1), &
             aradtogra * giph(2 * i, 2 * j + 1)
      END IF
      IF (depto(i, j) .GT. 1._wp) mende = mende + 1
      ftwou(i, j) = 2._wp * omega * sin(giph(2 * i + 1, 2 * j))
      ftwov(i, j) = 2._wp * omega * sin(giph(2 * i, 2 * j + 1))
    END DO
  END DO


  rmindlxp = 10.e9_wp
  rmindlyp = 10.e9_wp
  DO I = 1, IE
    DO J = 1, JE
      IF (depto(i, j) .GT. 0.5_wp) THEN
        rmindlxp = min(dlxp(i, j), rmindlxp)
        rmindlyp = min(dlyp(i, j), rmindlyp)
      END IF
    END DO
  END DO

  IF (verbose .GT. 0) PRINT *, 'min dlxp ', rmindlxp
  IF (verbose .GT. 0) PRINT *, 'min dlyp ', rmindlyp


  IF (verbose .GT. 0) PRINT *, ' sum of depth points ', mende

  WRITE(io_out_arcgri) 0, 507, - 100, ie * je
  WRITE(io_out_arcgri) depto
  WRITE(io_out_arcgri) 0, 85, - 100, ie * je
  WRITE(io_out_arcgri) dlxp
  WRITE(io_out_arcgri) 0, 185, - 100, ie * je
  WRITE(io_out_arcgri) dlxu
  WRITE(io_out_arcgri) 0, 285, - 100, ie * je
  WRITE(io_out_arcgri) dlxv
  WRITE(io_out_arcgri) 0, 86, - 100, ie * je
  WRITE(io_out_arcgri) dlyp
  WRITE(io_out_arcgri) 0, 186, - 100, ie * je
  WRITE(io_out_arcgri) dlyu
  WRITE(io_out_arcgri) 0, 286, - 100, ie * je
  WRITE(io_out_arcgri) dlyv
  WRITE(io_out_arcgri) 0, 175, - 100, ie * je
  WRITE(io_out_arcgri) ftwou
  WRITE(io_out_arcgri) 0, 176, - 100, ie * je
  WRITE(io_out_arcgri) ftwov

END PROGRAM SAM
