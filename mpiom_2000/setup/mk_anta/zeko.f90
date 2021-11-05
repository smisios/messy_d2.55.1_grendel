PROGRAM zeko
  USE mo_kind, ONLY: sp, i4, i8, wp
  USE mo_constants, ONLY: api, aradtogra, agratorad
  USE iso_varying_string
  IMPLICIT NONE
  ! .. Parameters ..

  !     dimension of etopo2
  INTEGER, PARAMETER :: me= 10800
  INTEGER, PARAMETER :: ne = 5401
  ! luns
  INTEGER, PARAMETER :: io_out_topo=66, io_in_etopo2=21, io_out_anta=12, &
       io_out_depto=67, io_out_dry=50, io_out_wet=51
  ! .. Local Scalars ..
  COMPLEX(wp) :: ci, cr, czer
  REAL(wp) :: al, al1, al2, ala1, ala2, ala3, ala4, alo1, alo2, alo3, alo4, &
       alr, antopo, betr, dphi, dryd, forz, forzi, geola, &
       geoph, gilama, gilami, giphma, giphmi, &
       onor1, onor2, phi, phi1, phi2, phi_i, pm, &
       redred, rlat1, rlat2, rlon1, rlon2, sud, suw, wetd, x, xi, xku, xp1, &
       xp2, xst, xv1, xv2, xz, y, yi, yku, yp1, yp2, ys, yst, yv1, yv2, yz, &
       z, zku, zp1, zp2, zs, zst, zv1, zv2, zz
  INTEGER :: i, ih, ii, ii1, ii2, imerc, isum, j, jj, k, l, &
       ie, je, ito, jto, iunit, ioerrstat, verbose
  INTEGER(i8) :: ext8_hdr(4)
  INTEGER(i4) :: ext4_hdr(4)
  LOGICAL :: write_ascii_drywet, write_ascii_depto
  ! .. Local Arrays ..

  COMPLEX(wp), ALLOCATABLE :: cgrid(:,:)
  REAL(wp), ALLOCATABLE :: depto(:,:), feld(:,:), gila(:,:), giph(:,:), &
       xsp(:,:), ysp(:,:), zsp(:,:)
  REAL(sp), ALLOCATABLE :: feld4(:, :)
  INTEGER,ALLOCATABLE :: mask(:,:), naf(:,:)
  ! file names
  TYPE(varying_string) :: etopo2_fname, topo_out_fname, tempstr
  ! .. Intrinsic Functions ..
  INTRINSIC acos, aimag, asin, atan2, cos, max, min, &
       mod, nint, REAL, sin, sqrt, tan!, abs
  NAMELIST /gridparams/ ie, je, rlat1, rlon1, rlat2, rlon2, phi, verbose

  ie = -1
  je = -1
  rlat1 = -1.0e10_wp
  rlon1 = -1.0e10_wp
  rlat2 = -1.0e10_wp
  rlon2 = -1.0e10_wp
  phi = -1.0e10_wp
  verbose = 0
  READ (*, nml=gridparams)
  IF (ie == -1 .OR. je == -1 .OR. rlat1 < -1000.0_wp .OR. rlon1 < -1000.0_wp &
       .OR. rlat2 < -1000.0_wp .OR. rlon2 < -1000.0_wp .OR. phi < -1000.0_wp) THEN
    PRINT *, 'illegal/missing value for ie, je, rlat[12], rlon[12] or phi'
    STOP 1
  END IF
  PRINT *, 'etopo2 file name = ?'
  CALL get(etopo2_fname, iostat=ioerrstat)
  IF (ioerrstat .GT. 0 .OR. CHAR(etopo2_fname) == '') &
       etopo2_fname='etopo2.ext'

  PRINT *, 'topo output file name = ?'
  CALL get(topo_out_fname, iostat=ioerrstat)
  IF (ioerrstat .GT. 0 .OR. CHAR(topo_out_fname) == '') &
       topo_out_fname='topo'

  PRINT *, 'write ascii description of dry and wet grid points?'
  CALL get(tempstr, iostat=ioerrstat)
  IF (ioerrstat .GT. 0 .OR. CHAR(tempstr, 1) /= 'y') THEN
    write_ascii_drywet=.FALSE.
  ELSE
    write_ascii_drywet=.TRUE.
  END IF

  PRINT *, 'write ascii representation of depto array?'
  CALL get(tempstr, iostat=ioerrstat)
  IF (ioerrstat .GT. 0 .OR. CHAR(tempstr, 1) /= 'y') THEN
    write_ascii_depto=.FALSE.
  ELSE
    write_ascii_depto=.TRUE.
  END IF

  IF (verbose .GT. 0) THEN
    PRINT *, 'ie = ', ie
    PRINT *, 'je = ', je
    PRINT *, 'rlat1 = ', rlat1
    PRINT *, 'rlon1 = ', rlon1
    PRINT *, 'rlat2 = ', rlat2
    PRINT *, 'rlon2 = ', rlon2
    PRINT *, 'phi   = ', phi
  END IF

  ito=2*ie
  jto=2*je

  !allocate memory
  ALLOCATE (cgrid(ito,0:jto))
  ALLOCATE (mask(ie,je), naf(me,ne))
  ALLOCATE (depto(ie,je), feld(me,ne), gila(ito,jto), giph(ito,jto), &
       xsp(ito,0:jto), ysp(ito,0:jto), zsp(ito,0:jto))

  cr = (1._wp,0._wp)
  ci = (0._wp,1._wp)

  OPEN (io_in_etopo2, file=CHAR(etopo2_fname), form='unformatted', &
        access='sequential', action='read')
  OPEN (io_out_anta, file='anta', form='unformatted', &
        access='sequential', action='write')

  PRINT *, 'read topography .. '
  ! try to read at 4 byte size first
  READ (io_in_etopo2) ext4_hdr
  IF (ext4_hdr(4) == me * ne) THEN
    ALLOCATE(feld4(me, ne))
    READ (io_in_etopo2) feld4
    feld = REAL(feld4, wp)
    DEALLOCATE(feld4)
  ELSE
    REWIND (io_in_etopo2)
    READ (io_in_etopo2) ext8_hdr
    IF (ext8_hdr(4) /= INT(me * ne, i8)) THEN
      PRINT *, 'Unexpected input size: ', ext8_hdr(4)
    END IF
    READ (io_in_etopo2) feld
  END IF
  CLOSE (io_in_etopo2)
  PRINT *, 'read topography .. done'

  !hh  E/W +180deg
  ih = me/2
  DO i = 1, ih
     DO j = 1, ne
       naf(i,j) = INT(feld(i+ih-1,j))
       naf(i+ih-1,j) = INT(feld(i,j))
     END DO
  END DO
  DO i = 1, me
     DO j = 1, ne
       feld(i,j) = REAL(naf(i,j), wp)
     END DO
  END DO

  phi1 = rlat1 * agratorad
  phi2 = rlat2 * agratorad
  al1 = rlon1 * agratorad
  al2 = rlon2 * agratorad

  xp1 = COS(phi1)*COS(al1)
  yp1 = COS(phi1)*SIN(al1)
  zp1 = SIN(phi1)
  xp2 = COS(phi2)*COS(al2)
  yp2 = COS(phi2)*SIN(al2)
  zp2 = SIN(phi2)

  !      tangential point of projection plane
  xz = 0.5_wp*(xp1+xp2)
  yz = 0.5_wp*(yp1+yp2)
  zz = 0.5_wp*(zp1+zp2)

!  PRINT *, ' pol1 ', xp1, yp1, zp1
!  PRINT *, ' pol2 ', xp2, yp2, zp2
  betr = SQRT(xz**2+yz**2+zz**2)
  IF (betr==0._wp) THEN
    PRINT *, 'pure rotation intended '
    PRINT *, ' TOO MUCH ALGEBRA FOR SIMPLE PROBLEM '
    STOP
  END IF
  xz = xz/betr
  yz = yz/betr
  zz = zz/betr
  !      two other vectors for orthonormal base
  xv1 = yp1*zp2 - yp2*zp1
  yv1 = zp1*xp2 - zp2*xp1
  zv1 = xp1*yp2 - xp2*yp1
  betr = SQRT(xv1**2+yv1**2+zv1**2)
  xv1 = xv1/betr
  yv1 = yv1/betr
  zv1 = zv1/betr

  xv2 = yz*zv1 - yv1*zz
  yv2 = zz*xv1 - zv1*xz
  zv2 = xz*yv1 - xv1*yz
  !      normalsystem
!  PRINT *, xz, yz, zz
!  PRINT *, xv1, yv1, zv1
!  PRINT *, xv2, yv2, zv2
  onor1 = xv1**2 + yv1**2 + zv1**2
  onor2 = xv2**2 + yv2**2 + zv2**2
  PRINT *, ' length of normal vectors ', onor1, onor2

  antopo = ACOS(xz*xp1+yz*yp1+zz*zp1)
  redred = 2._wp*TAN(0.5_wp*antopo)

  PRINT *, ' reductionfactor ', redred

  forz = REAL(ie-2, wp)/api

  forzi = 1._wp/forz
  czer = (0._wp,0._wp)

  DO  i = 1, ito
    imerc = 1
    !1 for quadratic grid
    phi_i = phi
    dphi = 0._wp
    DO  j = 0, jto
      IF (imerc==1) THEN
        dphi = COS(phi_i)*forzi
        phi_i = phi_i - dphi
      ELSE
        phi_i = REAL(j-(ie-2)/2, wp) * agratorad
      END IF

      cgrid(i,j) = CMPLX(1.E12_wp, 0._wp, wp)
      al = -forzi*REAL(i, wp)
      z = SIN(phi_i)
      y = COS(phi_i)*SIN(al)
      x = COS(phi_i)*COS(al)
      IF (x/=1._wp) THEN
        ys = y/(1._wp - x)
        zs = z/(1._wp - x)
        cgrid(i,j) = cr*CMPLX(zs, 0._wp, wp) + ci*CMPLX(ys, 0._wp, wp)
      END IF
    ENDDO
  ENDDO

  DO j = 1, jto
    DO i = 1, ito
      cgrid(i,j) = cgrid(i,j)*CMPLX(redred, 0._wp, wp)
      xi = REAL(cgrid(i,j), wp)
      yi = AIMAG(cgrid(i,j))
      xst = xz + xi*xv2 + yi*xv1
      yst = yz + xi*yv2 + yi*yv1
      zst = zz + xi*zv2 + yi*zv1
      alr = 2._wp*(xz**2+yz**2+zz**2+xz*xst+yz*yst+zz*zst)/ &
           ((xz+xst)**2+(yz+yst)**2+(zz+zst)**2)
      xku = -xz + alr*(xst+xz)
      yku = -yz + alr*(yst+yz)
      zku = -zz + alr*(zst+zz)
      xsp(i,j) = xku
      ysp(i,j) = yku
      zsp(i,j) = zku
      betr = xku**2 + yku**2 + zku**2
      geoph = ASIN(zku)
      geola = ATAN2(yku,xku)
      gila(i,j) = geola
      giph(i,j) = geoph
      IF (giph(i,j)>1.6_wp) WRITE (6,*) 'ATTN: giph ', i, j, geoph, zku, yku, xku
    END DO
  END DO

  DO j = 2, je - 1
    DO i = 2, ie - 1

      suw = 0._wp
      sud = 0._wp
      dryd = 0._wp
      wetd = 0._wp
      mask(i,j) = 0
      DO ii = -1, 1
        DO  jj = -1, 1

          geoph = ASIN(zsp(2*i+ii,2*j+jj))
          geola = ATAN2(ysp(2*i+ii,2*j+jj),xsp(2*i+ii,2*j+jj))

          al = geola*REAL(ne, wp)/api + REAL(me, wp)

          pm = (0.5_wp * api -geoph)*REAL(ne, wp)/api

          k = MOD(INT(al),me) + 1
          l = INT(pm)
          IF (i==120) PRINT *, j, ii, jj, geoph, geola, k, l
          IF (REAL(naf(k,l), wp)<0._wp) wetd = wetd - REAL(naf(k,l), wp)
          IF (REAL(naf(k,l), wp)>0._wp) dryd = dryd + REAL(naf(k,l), wp)
          IF (REAL(naf(k,l), wp)<0._wp) suw = suw + 1._wp
          IF (REAL(naf(k,l), wp)>0._wp) sud = sud + 1._wp

        ENDDO
      ENDDO

      IF (wetd>dryd) mask(i,j) = 13
      depto(i,j) = 0._wp
      IF (wetd>dryd) depto(i,j) = wetd/suw

!!$              IF (j==je/2) THEN
!!$                 WRITE (6,*) 'mitte ', i, gila(2*i,2*j)*180./3.14159265, &
!!$                      giph(2*i,2*j)*180./3.14159265, (180./3.14159265)*111.111*SQRT(( &
!!$                      giph(2*i,2*j)-giph(2*i-2,2*j))**2+COS(giph(2*i, &
!!$                      2*j))**2*(gila(2*i,2*j)-gila(2*i-2,2*j))**2), NINT(depto(i,j))
!!$              END IF

    ENDDO
  ENDDO

!!$        j = 2
!!$        DO i = 1, ie
!!$           WRITE (6,*) 'ANFANGSGIPH: ', i, j, ' giph,gila: ', giph(2*i,2*j-1), &
!!$                giph(2*i,2*j), gila(2*i,2*j-1), gila(2*i,2*j)
!!$        END DO

!!$        DO j = 1, 2*je
!!$           DO i = 1, 2*ie
!!$              IF (ABS(giph(i,j))<1.E-5) THEN
!!$                 WRITE (6,*) 'scheiss-giph ', i, j, giph(i,j)
!!$              END IF
!!$           END DO
!!$        END DO
  ext8_hdr(1) = 0_i8
  ext8_hdr(2) = 54_i8
  ext8_hdr(3) = 1_i8

  ext8_hdr(4) = INT(ito * jto, i8)
  WRITE (io_out_anta) ext8_hdr
  WRITE (io_out_anta) ((gila(i,j),i=1,ito),j=1,jto)

  ext8_hdr(2) = 55_i8
  WRITE (io_out_anta) ext8_hdr
  WRITE (io_out_anta) ((giph(i,j),i=1,ito),j=1,jto)

  ext8_hdr(2) = 84_i8
  ext8_hdr(4) = INT(ie * je, i8)
  WRITE (io_out_anta) ext8_hdr
  WRITE (io_out_anta) depto

  DO j = 1, je
    isum = 0
    DO i = 2, ie - 1
      IF (depto(i,j)>1._wp) isum = isum + 1
    END DO
    IF(verbose .GT. 0) WRITE (6,*) 'NO. OF WET POINTS AT J=', j,' IS ',isum
  END DO
  giphma = 0._wp
  giphmi = 999._wp
  gilama = 0._wp
  gilami = 999._wp
  DO j = 1, jto
    DO i = 1, ito + 4
      IF (giph(i,j)>3.1416_wp/2._wp) WRITE (6,*) 'ATTN: GIPH ', i, j, giph(i,j)
      giphma = MAX(giph(i,j),giphma)
      giphmi = MIN(giph(i,j),giphmi)
      gilama = MAX(gila(i,j),gilama)
      gilami = MIN(gila(i,j),gilami)
    END DO
  END DO
  IF (verbose .GT. 0) WRITE (6,*) 'GILA: ', gilami, gilama
  IF (verbose .GT. 0) WRITE (6,*) 'GIPH: ', giphmi, giphma

  DO j = 1, je
    DO i = 1, ie
      depto(i,j) = MAX(0._wp,depto(i,j))

      !              IF (depto(i,j)>0.) THEN
      !                 IF (depto(i,j)<11.) depto(i,j) = 0.
      !              END IF

      !              IF (depto(i,j)>=11.) THEN
      !                 IF (depto(i,j)<22.) depto(i,j) = 22.
      !              END IF

    END DO
  END DO

  DO j = 1, je
    depto(1,j) = depto(ie-1,j)
    depto(ie,j) = depto(2,j)
  END DO

  IF (write_ascii_drywet) THEN
    OPEN(io_out_wet, file='wetpoints.txt', form='formatted', &
         action='write')
    OPEN(io_out_dry, file='drypoints.txt', form='formatted', &
         action='write')
    DO i = 2, ie - 1
      DO j = 2, je - 1

        IF (depto(i,j)<=0.5_wp) THEN
          WRITE (io_out_dry,'(1X,F8.2,1X,F8.2,1X,F6.0)') aradtogra * gila(2*i,2*j),  &
               aradtogra * giph(2*i,2*j), depto(i,j)
        ELSE
          WRITE (io_out_wet,'(1X,F8.2,1X,F8.2,1X,F6.0)') aradtogra * gila(2*i,2*j),  &
               aradtogra * giph(2*i,2*j), depto(i,j)
        END IF
      END DO
    END DO
    CLOSE(io_out_wet)
    CLOSE(io_out_dry)
  END IF

  IF (write_ascii_depto) THEN
    DO I=2,IE-1,2
      iunit=100
      WRITE(iunit,'(1a,1X,F8.2,1X,F8.2,1X,F6.0)')'>', aradtogra * gila(2*i,2)&
           ,aradtogra * giph(2*i,2)                                          &
           ,DEPTO(I,1)
      DO j=1,jE,2
        WRITE(iunit,'(1X,F8.2,1X,F8.2,1X,F6.0)')aradtogra * gila(2*i,2*j)    &
             ,aradtogra * giph(2*i,2*j)                                      &
             ,DEPTO(I,j)
      ENDDO
    ENDDO

    DO j=1,jE,2
      WRITE(iunit,'(1A,1X,F8.2,1X,F8.2,1X,F6.0)')'>',aradtogra * gila(2*i,2) &
           ,aradtogra * giph(2*i,2)                                          &
           ,DEPTO(1,j)
      DO I=1,IE,2
        WRITE(iunit,'(1X,F8.2,1X,F8.2,1X,F6.0)')aradtogra * gila(2*i,2*j)    &
             ,aradtogra * giph(2*i,2*j)                                      &
             ,depto(I,j)
      END DO
    END DO
  END IF


  OPEN (io_out_depto, file='depto', form='unformatted', action='write')
  ext8_hdr(1) = 20021021_i8
  ext8_hdr(2) = 84_i8
  ext8_hdr(3) = 1_i8
  ext8_hdr(4) = INT(ie*je, i8)
  WRITE (io_out_depto) ext8_hdr
  WRITE (io_out_depto) ((depto(i, j), i=1, ie), j=1, je)

  OPEN (io_out_topo, file=char(topo_out_fname), &
        form='formatted', action='write')
  DO ii1 = 2, ie - 1, 20
    ii2 = MIN(ii1+19,ie-1)
    WRITE (io_out_topo,*) 'STREIFEN ', ii1, ii2
    DO j = 1, je
      WRITE (io_out_topo,6388) j, (NINT(depto(i,j)),i=ii1,ii2)
    END DO
6388 FORMAT (I5,20I5)

  END DO


  IF (write_ascii_depto) THEN
    DO i=2,ie-1
      DO j=2,je-1
        IF (depto(i, j) .LE. 0.5_wp) THEN

          ala1 = aradtogra * giph(2*i-1,2*j-1)
          ala2 = aradtogra * giph(2*i+1,2*j-1)
          ala3 = aradtogra * giph(2*i+1,2*j+1)
          ala4 = aradtogra * giph(2*i-1,2*j+1)
          alo1 = aradtogra * gila(2*i-1,2*j-1)
          alo2 = aradtogra * gila(2*i+1,2*j-1)
          alo3 = aradtogra * gila(2*i+1,2*j+1)
          alo4 = aradtogra * gila(2*i-1,2*j+1)

          WRITE(8,'(1A,1X,F8.2,1X,F8.2)') '>',ALO1,ALA1
          WRITE(8,'(1X,F8.2,1X,F8.2)') ALO2,ALA2
          WRITE(8,'(1X,F8.2,1X,F8.2)') ALO3,ALA3
          WRITE(8,'(1X,F8.2,1X,F8.2)') ALO4,ALA4
          WRITE(8,'(1X,F8.2,1X,F8.2)') ALO1,ALA1
        END IF
      END DO
    END DO
  END IF


  STOP
END PROGRAM zeko

