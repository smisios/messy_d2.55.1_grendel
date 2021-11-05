PROGRAM runoff_omip
  USE mo_constants, ONLY: aradtogra
  !USE mo_grid_elementals, ONLY: suchij_2d
  USE mo_kind, ONLY: i8, i4, sp, dp, wp
  IMPLICIT NONE
  TYPE ipair
    INTEGER :: a, b
  END TYPE ipair
  INTEGER :: ie, je, ke, ito, jto!, itemp, jtemp
  INTEGER, PARAMETER :: input_dim(2) = (/ 320, 160 /)
  INTEGER :: i, icount, ipos, j, jpos, m, verbose
  INTEGER(i8) :: ext8_hdr(4)
  INTEGER(i4) :: ext4_hdr(4)
  REAL(wp) :: dist, rtest, rtest1
  REAL(wp) ri(input_dim(1), input_dim(2))
  REAL(sp) lsm(input_dim(1), input_dim(2))
  REAL(wp) area(input_dim(1), input_dim(2)), alat(input_dim(1), input_dim(2)), &
       alon(input_dim(1), input_dim(2))
  TYPE(ipair) :: suchij_cache(input_dim(1), input_dim(2))
  REAL(dp), ALLOCATABLE :: GILA(:, :), GIPH(:, :), riv(:, :), alat_out(:,:), &
       alon_out(:, :)
  REAL(dp), ALLOCATABLE :: WETO(:, :, :)
  LOGICAL :: tp_grid
  INTEGER, PARAMETER :: io_in_runoff=11, io_in_lsm=70, io_in_ribek=71, &
       io_in_anta=73, io_out_giriv=12
  NAMELIST /gridparams/ ie, je, ke, verbose, tp_grid
  INTERFACE
    SUBROUTINE suchij(ie, je, ke, alat, alon, k, ipos, jpos, dist, &
         wetref, weto, gila, giph)
      USE mo_kind, ONLY: wp
      REAL(wp), INTENT(in) :: alat, alon, wetref
      INTEGER, INTENT(in) :: k, ie, je, ke
      INTEGER, INTENT(out) :: ipos, jpos
      REAL(wp), INTENT(out) :: dist
      REAL(wp), INTENT(in) :: GILA(2*ie, 2*je), GIPH(2*ie, 2*je)
      REAL(wp), INTENT(in) :: WETO(IE, JE, KE)
    END SUBROUTINE suchij
  END INTERFACE
  ie = -1
  je = -1
  ke = -1
  verbose = 0
  READ (*, nml=gridparams)
  IF (verbose .GT. 0) PRINT *, 'ie=', ie, 'je=', je, 'ke=', ke
  IF (ie < 1 .OR. je < 1 .OR. ke < 1) THEN
    PRINT *, 'Only positive integral values may be used for ie, je and ke!'
    STOP 1
  END IF
  ito=2*ie
  jto=2*je

  ALLOCATE(GILA(ITO, JTO), GIPH(ITO, JTO))

  OPEN(io_in_runoff, file='runoff.ext', &
       form='unformatted', action='read')
  OPEN(io_out_giriv, file='giriv_omip365', form='unformatted', &
       action='write')

  OPEN(io_in_anta, file='anta.ext', form='unformatted', action='read')
  READ(io_in_anta) ext8_hdr
  READ(io_in_anta) GILA
  READ(io_in_anta) ext8_hdr
  READ(io_in_anta) GIPH
  CLOSE(io_in_anta)

  ALLOCATE(alon_out(ie, je), alat_out(ie, je))
  FORALL(i = 1:ie, j = 1:je)
    alat_out(i, j) = giph(2*i, 2*j) * aradtogra
    alon_out(i, j) = gila(2*i, 2*j) * aradtogra
  END FORALL
  !DEALLOCATE(giph, gila)

  ALLOCATE(WETO(IE, JE, KE), riv(ie, je))

  OPEN(io_in_ribek, file='RIBEK.ext', form='unformatted', action='read')
  READ(io_in_ribek) ext8_hdr
  READ(io_in_ribek) ((weto(I, J, 1), I=1, IE), J=1, JE)
  CLOSE(io_in_ribek)

  OPEN(io_in_lsm, file='land_sea_mask.ECMWF.ext4' &
       ,form='unformatted', action='read')
  READ(io_in_lsm) ext4_hdr
  READ(io_in_lsm) lsm
  CLOSE(io_in_lsm)

  !hh     remove lakes
  DO i = 1, input_dim(1)
    DO j = 1, input_dim(2)
      IF (i.GT.39.AND.i.LT.58)THEN
        IF (j.GT.109.AND.j.LT.124)THEN
          lsm(i, j)=1._sp
        END IF
      END IF
      lsm(i, j) = MERGE(1._sp, 0._sp, lsm(i, j) .GT. 0.5_sp)
    END DO
  END DO



  CALL surface_areas(area, alat, alon)

  DO i = 1, input_dim(1)
    IF (verbose .GT. 1) PRINT*, 'suchij:', i,'/',input_dim(1)
    DO j = 1, input_dim(2)
      CALL suchij(ie, je, ke, alat(i, j), alon(i, j), 1, &
           suchij_cache(i, j)%a, suchij_cache(i, j)%b, &
           dist, 1._wp, weto, gila, giph)
    END DO
  END DO

  !DO i = 1, input_dim(1)
  !  DO j = 1, input_dim(2)
  !    CALL suchij_2d(ie, je, alat(i, j), alon(i, j), &
  !         itemp, jtemp, dist, .TRUE., weto, alat_out, alon_out, tp_grid)
  !    IF (itemp /= suchij_cache(i, j)%a &
  !         .OR. jtemp /= suchij_cache(i, j)%b) THEN
  !      PRINT *, 'problem: suchij does not match for (', i, j, '),', &
  !           ' (', itemp, jtemp, ') vs. (', suchij_cache(i, j)%a, &
  !           suchij_cache(i, j)%b, ')'
  !    END IF
  !  END DO
  !END DO

  DO m=1, 365

    IF (verbose .GT. 1) PRINT*, 'day ', m

    rtest=0._wp

    riv = 0._wp

    icount=0

    READ(io_in_runoff) ext8_hdr
    READ(io_in_runoff) ri

    DO i=1, input_dim(1)
      DO j=1, input_dim(2)
        IF (lsm(i, j) .LT. 0.5_sp .AND. ri(i, j) .NE. 0._wp)THEN

          ri(i, j) = ri(i, j) * area(i, j)

          rtest = rtest + ri(i, j)

          icount = icount + 1


          ipos = suchij_cache(i, j)%a
          jpos = suchij_cache(i, j)%b

          riv(ipos, jpos) = riv(ipos, jpos) + ri(i, j)

        END IF
      END DO
    END DO

    rtest1 = SUM(riv(2:ie - 1, 2:je - 1))

    IF (verbose .GT. 1) PRINT *, m, rtest, rtest1, icount
    ext8_hdr(4) = INT(ie * je, i8)
    WRITE(io_out_giriv) ext8_hdr
    WRITE(io_out_giriv) riv
  ENDDO
  CLOSE(io_in_runoff)
  CLOSE(io_out_giriv)
END PROGRAM runoff_omip

!===========================================================
!
! Subroutine to calculate surface areas for each point
! of a T106 Gaussian grid.
!
! Author: S.Legutke DKRZ/MPI
! -------
!
!===========================================================
SUBROUTINE surface_areas(surfatm, ALAT, ALON)
  USE mo_constants, ONLY: agratorad
  USE mo_kind, ONLY: wp
  USE mo_planetary_constants, ONLY: radius
  IMPLICIT NONE
  INTEGER, PARAMETER :: nxatm=320, nyatm=160
  REAL(wp) alonate, alonatw, be, bn, bs, bw, dx
  INTEGER :: i, j, jj, jx, jy
  REAL(wp) :: ALATATM(NYATM), ALONATM(NXATM), SURFATM(NXATM, NYATM), &
       alat(nxatm, nyatm), alon(nxatm, nyatm)

  DATA ALATATM/ &
       -89.14151942647_wp, -88.02942886796_wp, -86.91077081413_wp, -85.79062888364_wp, -84.66992408445_wp, -83.54894691254_wp, &
       -82.42781752401_wp, -81.30659452267_wp, -80.18530987248_wp, -79.06398248141_wp, -77.94262424667_wp, -76.82124302710_wp, &
       -75.69984422201_wp, -74.57843166330_wp, -73.45700814558_wp, -72.33557575491_wp, -71.21413607989_wp, -70.09269035162_wp, &
       -68.97123953894_wp, -67.84978441467_wp, -66.72832560288_wp, -65.60686361301_wp, -64.48539886504_wp, -63.36393170834_wp, &
       -62.24246243589_wp, -61.12099129525_wp, -59.99951849704_wp, -58.87804422158_wp, -57.75656862418_wp, -56.63509183933_wp, &
       -55.51361398408_wp, -54.39213516079_wp, -53.27065545940_wp, -52.14917495922_wp, -51.02769373051_wp, -49.90621183571_wp, &
       -48.78472933054_wp, -47.66324626484_wp, -46.54176268341_wp, -45.42027862655_wp, -44.29879413069_wp, -43.17730922883_wp, &
       -42.05582395094_wp, -40.93433832428_wp, -39.81285237377_wp, -38.69136612220_wp, -37.56987959047_wp, -36.44839279779_wp, &
       -35.32690576187_wp, -34.20541849905_wp, -33.08393102445_wp, -31.96244335209_wp, -30.84095549500_wp, -29.71946746532_wp, &
       -28.59797927436_wp, -27.47649093270_wp, -26.35500245025_wp, -25.23351383632_wp, -24.11202509967_wp, -22.99053624854_wp, &
       -21.86904729073_wp, -20.74755823362_wp, -19.62606908420_wp, -18.50457984914_wp, -17.38309053477_wp, -16.26160114716_wp, &
       -15.14011169211_wp, -14.01862217519_wp, -12.89713260175_wp, -11.77564297696_wp, -10.65415330582_wp, -9.532663593176_wp, &
       -8.411173843743_wp, -7.289684062115_wp, -6.168194252784_wp, -5.046704420157_wp, -3.925214568566_wp, -2.803724702287_wp, &
       -1.682234825547_wp, -0.560744942544_wp, 0.560744942544_wp, 1.682234825547_wp, 2.803724702287_wp, 3.925214568566_wp, &
       5.046704420157_wp, 6.168194252784_wp, 7.289684062115_wp, 8.411173843743_wp, 9.532663593176_wp, 10.65415330582_wp, &
       11.77564297696_wp, 12.89713260175_wp, 14.01862217519_wp, 15.14011169211_wp, 16.26160114716_wp, 17.38309053477_wp, &
       18.50457984914_wp, 19.62606908420_wp, 20.74755823362_wp, 21.86904729073_wp, 22.99053624854_wp, 24.11202509967_wp, &
       25.23351383632_wp, 26.35500245025_wp, 27.47649093270_wp, 28.59797927436_wp, 29.71946746532_wp, 30.84095549500_wp, &
       31.96244335209_wp, 33.08393102445_wp, 34.20541849905_wp, 35.32690576187_wp, 36.44839279779_wp, 37.56987959047_wp, &
       38.69136612220_wp, 39.81285237377_wp, 40.93433832428_wp, 42.05582395094_wp, 43.17730922883_wp, 44.29879413069_wp, &
       45.42027862655_wp, 46.54176268341_wp, 47.66324626484_wp, 48.78472933054_wp, 49.90621183571_wp, 51.02769373051_wp, &
       52.14917495922_wp, 53.27065545940_wp, 54.39213516079_wp, 55.51361398408_wp, 56.63509183933_wp, 57.75656862418_wp, &
       58.87804422158_wp, 59.99951849704_wp, 61.12099129525_wp, 62.24246243589_wp, 63.36393170834_wp, 64.48539886504_wp, &
       65.60686361301_wp, 66.72832560288_wp, 67.84978441467_wp, 68.97123953894_wp, 70.09269035162_wp, 71.21413607989_wp, &
       72.33557575491_wp, 73.45700814558_wp, 74.57843166330_wp, 75.69984422201_wp, 76.82124302710_wp, 77.94262424667_wp, &
       79.06398248141_wp, 80.18530987248_wp, 81.30659452267_wp, 82.42781752401_wp, 83.54894691254_wp, 84.66992408445_wp, &
       85.79062888364_wp, 86.91077081413_wp, 88.02942886796_wp, 89.14151942647_wp  /


  DX  = 360._wp/REAL(NXATM, wp)

  DO J = 1, NYATM
    DO I = 1, NXATM
      jj=161-j
      ALONATM(I) = DX*REAL(I - 1, wp)
      alat(i, j) = ALATATM(jj)
      alon(i, j) = ALONATM(i)
    ENDDO
  ENDDO

  DO JX = 1, NXATM
    DO JY = 1, NYATM

      BN = MERGE(90._wp, (ALATATM(JY + 1) + ALATATM(JY)) * 0.5_wp, JY .EQ. NYATM)
      BS = MERGE(-90._wp, (ALATATM(JY - 1) + ALATATM(JY)) * 0.5_wp, JY .EQ. 1)
      ALONATE = MERGE(ALONATM(1), ALONATM(JX + 1), JX .EQ. NXATM)
      ALONATW = MERGE(ALONATM(NXATM), ALONATM(JX - 1), JX .EQ. 1)

      IF(ALONATE .LE. ALONATM(JX)) ALONATE = ALONATE + 360._wp
      IF(ALONATW .GE. ALONATM(JX)) ALONATW = ALONATW - 360._wp

      BW=MOD((ALONATM(JX)+ALONATW)*0.5_wp, 360._wp)
      BE=MOD((ALONATM(JX)+ALONATE)*0.5_wp, 360._wp)

      SURFATM(JX, JY) = (BE-BW)* radius * agratorad &
           * (BN-BS) * radius * agratorad &
           * COS(ALATATM(JY) * agratorad)
    ENDDO
  ENDDO
END SUBROUTINE surface_areas

SUBROUTINE SUCHIJ(ie, je, ke, alat, alon, k, ipos, jpos, dist, wetref, &
     weto, gila, giph)
  !
  !     subroutine to determine i and j indices of nearest point
  !     from lat and long
  !
  !     input:
  !     alat       latitude (deg.)
  !     alon       longitude (deg.)
  !     k          level (integer)
  !     wetref     if land shall be considered too : 0.
  !                ocean points only                 1.
  !
  !     output:
  !     ipos       i index
  !     jpos       j index
  !     dist       distance from point in m
  !
  USE mo_constants, ONLY: agratorad
  USE mo_kind, ONLY: wp
  IMPLICIT NONE
  REAL(wp), INTENT(in) :: alat, alon, wetref
  INTEGER, INTENT(in) :: k, ie, je, ke
  INTEGER, INTENT(out) :: ipos, jpos
  REAL(wp), INTENT(out) :: dist
  REAL(wp), INTENT(in) :: GILA(2*ie, 2*je), GIPH(2*ie, 2*je)
  REAL(wp), INTENT(in) ::  WETO(IE, JE, KE)
  REAL(wp) :: disti(ie, je), dmini(ie)
  INTEGER :: jmini(ie)
  REAL(wp) :: alam, aphi, ax, ay, az, dmi, xx, yy, zz
  INTEGER :: i, imi, j, jmi
  !       print*, 'debug: in suchij 1'
  aphi = agratorad * alat
  alam = agratorad * alon

  xx=cos(alam)*cos(aphi)
  yy=sin(alam)*cos(aphi)
  zz=sin(aphi)

  dmini=1.e20_wp
  disti=1.e20_wp
  dmi=1.E20_wp
  jmini=1
  imi=1
  jmi=1
  do j=3, je-2
    do i=2, ie-1
      !         print*, 'debug: in suchij 2', i, j, k, weto(i, j, k)
      IF(weto(i, j, k) .GT. (wetref-0.1_wp))THEN
        !           print*, 'debug: in suchij 2', i, j, weto(i, j, k)
        ax=cos(giph(2*i, 2*j))*cos(gila(2*i, 2*j))
        ay=cos(giph(2*i, 2*j))*sin(gila(2*i, 2*j))
        az=sin(giph(2*i, 2*j))
        disti(i, j)=(xx-ax)**2+(yy-ay)**2+(zz-az)**2
      endif
    enddo
  enddo
  do j=3, je-2
    do i=1, ie
      if(disti(i, j).lt.dmini(i))then

        jmini(i)=j
        dmini(i)=disti(i, j)
      endif
    enddo
  enddo

  dmi=1.e20_wp

  do i=1, ie
    if(dmi.gt.dmini(i))then
      dmi=dmini(i)
      imi=i
      jmi=jmini(i)
    endif
  enddo
  dist = SQRT(dmi)*6350000._wp
  ipos=imi
  jpos=jmi
  !       write(6, *)'in suchij ende: ', ipos, jpos, k, weto(ipos, jpos, k)
END SUBROUTINE SUCHIJ


