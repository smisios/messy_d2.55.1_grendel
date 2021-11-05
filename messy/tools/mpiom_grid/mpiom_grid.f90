PROGRAM MAIN

    ! ------------------------------------------------------------------------
    !    mpiom_grid / version 0.1 / POZZER ANDREA / MPI-C Mainz, 2009
    ! ------------------------------------------------------------------------
    !
    !   This code is used for preparing new grid (and initial)
    !   conditions for MPIOM.
    !
    !   Be aware: for low resolution grid, some manual adjustment of
    !             the topo file (for straits) may be necessary
    !
    !   Developer:  Andrea Pozzer, MPIC-Mainz / The Cyprus Institute, 2009
    !
    !   Based on original code at the MPI-M Hamburg.
    !   Pease refer to MPIOM author for any reference.
    !
    !
    !   BUGS: probably a lot of problems could arise with the "BIG_ENDIAN"
    !         and "LITTLE ENDIAN" format.
    !
    !   TODO: *) create input files not more machine dependent (netcf format)
    !            so we can get rid off the bugs.....
    !         *) Clean the code (a lot of work)
    !
    ! ----------------------------------------------------------------------

  USE netcdf

  !I/O related
!mz_ap_20141007+
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)
!mz_ap_20141007-

  INTEGER :: status ! status flag
  CHARACTER(LEN=*), PARAMETER :: hline='==========================================================='
  CHARACTER(LEN=*), PARAMETER :: file_nml="grid.nml"
  CHARACTER(LEN=200):: files_path = '' ! CTRL namelsit
  CHARACTER(LEN=200):: phc_files_path = '' ! CTRL namelsit
  CHARACTER(LEN=200):: grid_name = '' ! CTRL namelsit
  INTEGER :: IE,JE, rlat1,rlat2,rlon1,rlon2
  REAL    :: phiread
  INTEGER :: ito,jto
  INTEGER :: levels
  LOGICAl :: l_init_clim
  LOGICAl :: l_tripolar
  INTEGER, PARAMETER :: iou  = 21       ! I/O unit

  ! .. Intrinsic Functions ..
  INTRINSIC abs, acos, aimag, asin, atan, atan2, cos, float, max, min, &
       mod, nint, sin, sqrt, tan, TRIM

!--------------------------------------------------------------------------------
!  READ NAMELIST-FILE
!--------------------------------------------------------------------------------

  CALL read_nml(status, iou, TRIM(file_nml))
  IF (status /= 0) STOP

  ito=2*ie
  jto=2*je

!--------------------------------------------------------------------------------
!  CREATION OF ANTA AND TOPO
!--------------------------------------------------------------------------------

  CALL ANTA

!--------------------------------------------------------------------------------
!  CREATION OF ARCGRI (DLXP & DLYP)
!--------------------------------------------------------------------------------

  CALL ARCGRI

!--------------------------------------------------------------------------------
!  CREATION NETCDF FILES OF GRID (generally not used in MPIOM but useful)
!--------------------------------------------------------------------------------

  CALL ANTA2NC

!--------------------------------------------------------------------------------
!  CREATION SALINITY AND TEMPERATURE INITIAL CONDITIONS
!--------------------------------------------------------------------------------

  CALL INIT

!--------------------------------------------------------------------------------
!  END OF PROGRAM
!--------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE ANTA

  ! .. Parameters ..

  !     dimension of etopo2
  INTEGER,PARAMETER :: me= 10800
  INTEGER,PARAMETER :: ne = 5401
  ! .. GLOBAL Scalars ..
  COMPLEX :: ci, cr, czer
  REAL :: al, al1, al2, alr, antopo, betr, dphi, dryd, forz, forzi, geola, &
       geoph, gilama, gilami, giphma, giphmi, gl1, gl2, gp1, gp2, &
       grarad, onor1, onor2, phi, phi1, phi2, pi, pih, pm, &
       redred, sud, suw, wetd, x, xi, xku, xp1, &
       xp2, xst, xv1, xv2, xz, y, yi, yku, yp1, yp2, ys, yst, yv1, yv2, yz, &
       z, zku, zp1, zp2, zs, zst, zv1, zv2, zz
  INTEGER :: i, idumi1, idumito, idumi3, idumi4, ih, ii, ii1, iito, ii3, &
       ii4, imerc, isum, j, jj, k, l
  ! .. GLOBAL Arrays ..

  COMPLEX, ALLOCATABLE :: cgrid(:,:)
!mz_ap_20141007+
!  REAL,ALLOCATABLE :: depto(:,:), feld(:,:), gila(:,:), giph(:,:), &
!       xsp(:,:), ysp(:,:), zsp(:,:)
  REAL,ALLOCATABLE :: depto(:,:), gila(:,:), giph(:,:), &
       xsp(:,:), ysp(:,:), zsp(:,:)
  REAL(sp),ALLOCATABLE :: feld(:,:)
!mz_ap_20141007-


  INTEGER,ALLOCATABLE :: mask(:,:), naf(:,:)


  !allocate memory
  ALLOCATE (cgrid(ito,0:jto))
  ALLOCATE (mask(ie,je), naf(me,ne))
  ALLOCATE  (depto(ie,je), feld(me,ne), gila(ito,jto), giph(ito,jto), &
       xsp(ito,0:jto), ysp(ito,0:jto), zsp(ito,0:jto))

  cr = (1.,0.)
  ci = (0.,1.)

#ifndef NOENDIANCONVERT
  OPEN (21,file=TRIM(files_path)//'/etopo2.ext',convert='BIG_ENDIAN',form='UNFORMATTED',access='SEQUENTIAL')
  OPEN (12,file='./output/'//TRIM(grid_name)//'_anta',convert='BIG_ENDIAN',form='UNFORMATTED',access='SEQUENTIAL')
#else
  STOP 'ERROR: open(...,convert=BIG_ENDIAN,...) not supported!'
#endif

  PRINT *, 'read topography .. '
  READ (21) idumi1, idumi2, idumi3, idumi4
  READ (21) feld
  PRINT *, 'read topography .. done'

  !hh  E/W +180deg
  ih = me/2
  DO i = 1, ih
     DO j = 1, ne
        naf(i,j) = feld(i+ih-1,j)
        naf(i+ih-1,j) = feld(i,j)
     END DO
  END DO
  DO i = 1, me
     DO j = 1, ne
        feld(i,j) = naf(i,j)
        naf(i,j) = 0.
     END DO
  END DO

  pi = 4.*ATAN(1.)
  pih = 0.5*pi
  gp1 = rlat1
  gp2 = rlat2
  gl1 = rlon1
  gl2 = rlon2
  phi1 = gp1*pi/180.
  phi2 = gp2*pi/180.
  al1 = gl1*pi/180.
  al2 = gl2*pi/180.

  xp1 = COS(phi1)*COS(al1)
  yp1 = COS(phi1)*SIN(al1)
  zp1 = SIN(phi1)
  xp2 = COS(phi2)*COS(al2)
  yp2 = COS(phi2)*SIN(al2)
  zp2 = SIN(phi2)

  !      tangential point of projection plane
  xz = 0.5*(xp1+xp2)
  yz = 0.5*(yp1+yp2)
  zz = 0.5*(zp1+zp2)

!  PRINT *, ' pol1 ', xp1, yp1, zp1
!  PRINT *, ' pol2 ', xp2, yp2, zp2
  betr = SQRT(xz**2+yz**2+zz**2)
  IF (betr==0.) THEN
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
  redred = 2.*TAN(0.5*antopo)

  antopo = ACOS(xz*xp1+yz*yp1+zz*zp1)
  redred = 2.*TAN(0.5*antopo)

  PRINT *, ' reductionfactor ', redred
  phi = 0.2
  PRINT *, ' length of normal vectors ', onor1, onor2

  antopo = ACOS(xz*xp1+yz*yp1+zz*zp1)
  redred = 2.*TAN(0.5*antopo)

  antopo = ACOS(xz*xp1+yz*yp1+zz*zp1)
  redred = 2.*TAN(0.5*antopo)

  PRINT *, ' reductionfactor ', redred
  forz = float(ie-2)/pi

  forzi = 1./forz
  phi = 0.
  czer = (0.,0.)

  DO  i = 1, ito
     imerc = 1
     !1 for quadratic grid
     phi = phiread
     dphi = 0.
     DO  j = 0, jto
        IF (imerc==1) THEN
           dphi = COS(phi)*forzi
           phi = phi - dphi
        ELSE
           phi = (j-(ie-2)/2)*pi/180.
        END IF

        cgrid(i,j) = 1.E12
        al = -forzi*i
        z = SIN(phi)
        y = COS(phi)*SIN(al)
        x = COS(phi)*COS(al)
        IF (x/=1.) THEN
           ys = y/(1-x)
           zs = z/(1-x)
           cgrid(i,j) = cr*zs + ci*ys
        END IF
     ENDDO
  ENDDO

  DO j = 1, jto
     DO i = 1, ito
        cgrid(i,j) = cgrid(i,j)*redred
        xi = REAL(cgrid(i,j))
        yi = AIMAG(cgrid(i,j))
        xst = xz + xi*xv2 + yi*xv1
        yst = yz + xi*yv2 + yi*yv1
        zst = zz + xi*zv2 + yi*zv1
        alr = 2.*(xz**2+yz**2+zz**2+xz*xst+yz*yst+zz*zst)/ &
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
        IF (giph(i,j)>1.6) WRITE (6,*) 'ATTN: giph ', i, j, geoph, zku, yku, xku
 end do
end do

        DO j = 2, je - 1
           DO i = 2, ie - 1

              suw = 0.
              sud = 0.
              dryd = 0.
              wetd = 0.
              mask(i,j) = 0
              DO ii = -1, 1
                 DO  jj = -1, 1

                    geoph = ASIN(zsp(2*i+ii,2*j+jj))
                    geola = ATAN2(ysp(2*i+ii,2*j+jj),xsp(2*i+ii,2*j+jj))

                    al = geola*ne/pi + me

                    pm = (pih-geoph)*ne/pi

                    k = MOD(int(al),me) + 1
                    l = pm
                    ! IF (i==120) PRINT *, j, ii, jj, geoph, geola, k, l
                    IF (feld(k,l)<0.) wetd = wetd - feld(k,l)
                    IF (feld(k,l)>0.) dryd = dryd + feld(k,l)
                    IF (feld(k,l)<0.) suw = suw + 1.
                    IF (feld(k,l)>0.) sud = sud + 1.

                 ENDDO
              ENDDO

              IF (wetd>dryd) mask(i,j) = 13
              depto(i,j) = 0.
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
        WRITE (12) 0, 54, -99, ito*jto
        WRITE (12) ((gila(i,j),i=1,ito),j=1,jto)
        WRITE (12) 0, 55, -99, ito*jto
        WRITE (12) ((giph(i,j),i=1,ito),j=1,jto)
        WRITE (12) 0, 507, -99, ie*je
        WRITE (12) depto

        DO j = 1, je
           isum = 0.
           DO i = 2, ie - 1
              IF (depto(i,j)>1) isum = isum + 1
           END DO
           WRITE (6,*) 'NO. OF WET POINTS AT J=', j,' IS ',isum
        END DO
        giphma = 0.
        giphmi = 999.
        gilama = 0.
        gilami = 999.
        DO j = 1, jto
           DO i = 1, ito
              IF (giph(i,j)>3.1416/2.) WRITE (6,*) 'ATTN: GIPH ', i, j, giph(i,j)
              giphma = MAX(giph(i,j),giphma)
              giphmi = MIN(giph(i,j),giphmi)
              gilama = MAX(gila(i,j),gilama)
              gilami = MIN(gila(i,j),gilami)
           END DO
        END DO
        WRITE (6,*) 'GILA: ', gilami, gilama
        WRITE (6,*) 'GIPH: ', giphmi, giphma

        DO j = 1, je
           DO i = 1, ie
              depto(i,j) = MAX(0.,depto(i,j))

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

        grarad = 180./3.1415927

!        DO i = 2, ie - 1
!           DO j = 2, je - 1
!
!              IF (depto(i,j)<=0.5) THEN
!                 WRITE (50,'(1X,F8.2,1X,F8.2,1X,F6.0)') grarad*gila(2*i,2*j), &
!                      grarad*giph(2*i,2*j), depto(i,j)
!              ELSE
!                 WRITE (51,'(1X,F8.2,1X,F8.2,1X,F6.0)') grarad*gila(2*i,2*j), &
!                      grarad*giph(2*i,2*j), depto(i,j)
!              END IF
!           END DO
!        END DO
!
!       DO I=2,IE-1,2
!       iunit=100
!       WRITE(iunit,'(1a,1X,F8.2,1X,F8.2,1X,F6.0)')'>', grarad*gila(2*i,2)     &
!            ,grarad*giph(2*i,2)                                               &
!                 ,DEPTO(I,1)
!       DO j=1,jE,2
!        WRITE(iunit,'(1X,F8.2,1X,F8.2,1X,F6.0)')grarad*gila(2*i,2*j)    &
!           ,grarad*giph(2*i,2*j)                                       &
!                 ,DEPTO(I,j)
!       enddo
!       enddo
!
!       DO j=1,jE,2
!       WRITE(iunit,'(1A,1X,F8.2,1X,F8.2,1X,F6.0)')'>',grarad*gila(2*i,2) &
!                 ,grarad*giph(2*i,2)                                     &
!                 ,DEPTO(1,j)
!       DO I=1,IE,2
!        WRITE(iunit,'(1X,F8.2,1X,F8.2,1X,F6.0)')grarad*gila(2*i,2*j)    &
!           ,grarad*giph(2*i,2*j)                                       &
!                 ,DEPTO(I,j)
!       enddo
!       enddo



!        OPEN (67,file='./output/depto',form='unformatted')
!        ii1 = 20021021
!        ii2 = 84
!        ii3 = 1
!        ii4 = ie*je
!        WRITE (67) ii1, ii2, ii3, ii4
!        WRITE (67) ((depto(i,j),i=1,ie),j=1,je)
!
        DO ii1 = 2, ie - 1, 20
           ii2 = MIN(ii1+19,ie-1)
#ifdef LITTLE_ENDIAN
#ifndef NOENDIANCONVERT
           OPEN (66,file='./output/'//TRIM(grid_name)//'_topo',convert='BIG_ENDIAN',form='formatted')
#else
           STOP 'ERROR: open(...,convert=BIG_ENDIAN,...) not supported!'
#endif
#else
           OPEN (66,file='./output/'//TRIM(grid_name)//'_topo',form='formatted')
#endif
           WRITE (66,*) 'STREIFEN ', ii1, ii2
           DO j = 1, je
              WRITE (66,6388) j, (NINT(depto(i,j)),i=ii1,ii2)
           END DO
6388       FORMAT (I5,20I5)

        END DO


!        DO I=2,IE-1
!         DO j=2,jE-1
!         tt=DEPTO(I,J)
!            if (tt.le.0.5) then
!            ALA1=GRARAD*GIPH(2*I-1,2*j-1)
!            ALA2=GRARAD*GIPH(2*I+1,2*j-1)
!            ALA3=GRARAD*GIPH(2*I+1,2*j+1)
!            ALA4=GRARAD*GIPH(2*I-1,2*j+1)
!            ALo1=GRARAD*GIla(2*I-1,2*j-1)
!            ALo2=GRARAD*GIla(2*I+1,2*j-1)
!            ALo3=GRARAD*GIla(2*I+1,2*j+1)
!            ALo4=GRARAD*GIla(2*I-1,2*j+1)
!            WRITE(8,'(1A,1X,F8.2,1X,F8.2)') '>',ALO1,ALA1
!            WRITE(8,'(1X,F8.2,1X,F8.2)') ALO2,ALA2
!            WRITE(8,'(1X,F8.2,1X,F8.2)') ALO3,ALA3
!            WRITE(8,'(1X,F8.2,1X,F8.2)') ALO4,ALA4
!            WRITE(8,'(1X,F8.2,1X,F8.2)') ALO1,ALA1
!            endif
!         enddo
!        enddo
! --------------------------------------------------------------------------
! ##########################################################################
! --------------------------------------------------------------------------
  DEALLOCATE (cgrid)
  DEALLOCATE (mask, naf)
  DEALLOCATE (depto, feld, gila, giph, xsp, ysp, zsp)
       ! STOP

  END SUBROUTINE ANTA

  ! ------------------------------------------------------------------------

  SUBROUTINE ARCGRI

      IMPLICIT NONE

      REAL :: MASKO(IE,JE)
      REAL :: dlxp(ie,je),dlyp(ie,je),dlxv(ie,je),dlyu(ie,je) &
              ,dlxu(ie,je),dlyv(ie,je)
      REAL :: ftwou(ie,je),ftwov(ie,je)
      REAL :: GILA(ITO,JTO),GIPH(ITO,JTO)
      REAL :: DEPTO(IE,JE)
      INTEGER :: JE1
      INTEGER :: IE1
      INTEGER :: i,ii,j

! mz_pj_20090730+
!!$      REAL, PARAMETER :: erdrad=6378000.
!!$      REAL, PARAMETER :: PI=4.*ATAN(1.)
!!$      REAL, PARAMETER :: pih=0.5*pi
!!$      REAL, PARAMETER :: zwepi=2.*pi
!!$      REAL, PARAMETER :: pigra=45./atan(1.)
!!$      REAL, PARAMETER :: DPHI=222000.
!!$      REAL, PARAMETER :: omega=2.*pi/86164.
      REAL, PARAMETER :: erdrad=6378000.
      REAL            :: PI
      REAL            :: pih
      REAL            :: zwepi
      REAL            :: pigra
      REAL, PARAMETER :: DPHI=222000.
      REAL            :: omega
! mz_pj_20090730-

      INTEGER :: ia1,ia2,ia3,ia4
      REAL :: siphup
      REAL :: cophup
      REAL :: siphum
      REAL :: cophum
      REAL :: siphvp
      REAL :: cophvp
      REAL :: siphvm
      REAL :: cophvm
      REAL :: cophzz
      REAL :: cophsw
      REAL :: cophss
      REAL :: cophse
      REAL :: cophne
      REAL :: cophee
      REAL :: siphzz
      REAL :: siphsw
      REAL :: siphss
      REAL :: siphse
      REAL :: siphne
      REAL :: siphee
      REAL :: silaup
      REAL :: colaup
      REAL :: silaum
      REAL :: colaum
      REAL :: silavp
      REAL :: colavp
      REAL :: silavm
      REAL :: colavm
      REAL :: colazz
      REAL :: colasw
      REAL :: colass
      REAL :: colase
      REAL :: colane
      REAL :: colaee
      REAL :: silazz
      REAL :: silasw
      REAL :: silass
      REAL :: silase
      REAL :: silane
      REAL :: silaee
      REAL :: xup
      REAL :: yup
      REAL :: zup

      REAL :: xum
      REAL :: yum
      REAL :: zum

      REAL :: xvp
      REAL :: yvp
      REAL :: zvp

      REAL :: xvm
      REAL :: yvm
      REAL :: zvm


      REAL :: xzz
      REAL :: yzz
      REAL :: zzz


      REAL :: xee
      REAL :: yee
      REAL :: zee


      REAL :: xss
      REAL :: yss
      REAL :: zss


      REAL :: xsw
      REAL :: ysw
      REAL :: zsw


      REAL :: xse
      REAL :: yse
      REAL :: zse

      REAL :: xne
      REAL :: yne
      REAL :: zne

      INTEGER :: mende
      REAL :: DELLG
      REAL :: RMINDLXP
      REAL :: RMINDLYP

      ! mz_pj_20090730+
      PI=4.*ATAN(1.)
      pih=0.5*pi
      zwepi=2.*pi
      pigra=45./atan(1.)
      omega=2.*pi/86164.
      ! mz_pj_20090730-

#ifndef NOENDIANCONVERT
      OPEN(81,FILE='./output/'//TRIM(grid_name)//'_anta',convert='BIG_ENDIAN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(82,FILE='./output/'//TRIM(grid_name)//'_arcgri',convert='BIG_ENDIAN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
#else
  STOP 'ERROR: open(...,convert=BIG_ENDIAN,...) not supported!'
#endif

      print*, ' vor gila '
      read(81) ia1,ia2,ia3,ia4
      read(81) gila

      print*, ' vor giph'
      read(81)ia1,ia2,ia3,ia4
      read(81) giph


      print*, ' vor depto '
      read(81)ia1,ia2,ia3,ia4
      read(81) depto

      JE1=JE-1
      IE1=IE-1

      do 7722 j=1,je
      do 7722 i=1,ie
      masko(i,j)=0
      if(depto(i,j).gt.1.) masko(i,j)=1
7722  continue
      do 5867 j=1,jto
      do 5867 i=2,ito
      dellg=gila(i,j)-gila(i-1,j)
      if(abs(dellg).lt.pi) go to 5867
      if(dellg.lt.-pi) then
      do 5868 ii=i,ito
       gila(ii,j)=gila(ii,j)+zwepi

5868   continue
     else
      do 5869 ii=i,ito
      gila(ii,j)=gila(ii,j)-zwepi
5869  continue
      endif
5867  continue


      mende=0
65543 format(20f6.0)

      do 1478 j=1,je1
      do 1478 i=1,ie1


      siphup=sin(giph(2*i+1,2*j))
      cophup=cos(giph(2*i+1,2*j))
      siphum=sin(giph(2*i-1,2*j))
      cophum=cos(giph(2*i-1,2*j))
      siphvp=sin(giph(2*i,2*j+1))
      cophvp=cos(giph(2*i,2*j+1))
      siphvm=sin(giph(2*i,2*j-1))
      cophvm=cos(giph(2*i,2*j-1))

      cophzz=cos(giph(2*i  ,2*j))
      cophsw=cos(giph(2*i-1,2*j+1))
      cophss=cos(giph(2*i  ,2*j+2))
      cophse=cos(giph(2*i+1,2*j+1))
      cophne=cos(giph(2*i+1,2*j-1))
      cophee=cos(giph(2*i+2,2*j))

      siphzz=sin(giph(2*i  ,2*j))
      siphsw=sin(giph(2*i-1,2*j+1))
      siphss=sin(giph(2*i  ,2*j+2))
      siphse=sin(giph(2*i+1,2*j+1))
      siphne=sin(giph(2*i+1,2*j-1))
      siphee=sin(giph(2*i+2,2*j))

      silaup=sin(gila(2*i+1,2*j))
      colaup=cos(gila(2*i+1,2*j))
      silaum=sin(gila(2*i-1,2*j))
      colaum=cos(gila(2*i-1,2*j))
      silavp=sin(gila(2*i,2*j+1))
      colavp=cos(gila(2*i,2*j+1))
      silavm=sin(gila(2*i,2*j-1))
      colavm=cos(gila(2*i,2*j-1))


      colazz=cos(gila(2*i  ,2*j))
      colasw=cos(gila(2*i-1,2*j+1))
      colass=cos(gila(2*i  ,2*j+2))
      colase=cos(gila(2*i+1,2*j+1))
      colane=cos(gila(2*i+1,2*j-1))
      colaee=cos(gila(2*i+2,2*j))

      silazz=sin(gila(2*i  ,2*j))
      silasw=sin(gila(2*i-1,2*j+1))
      silass=sin(gila(2*i  ,2*j+2))
      silase=sin(gila(2*i+1,2*j+1))
      silane=sin(gila(2*i+1,2*j-1))
      silaee=sin(gila(2*i+2,2*j))
      xup=silaup*cophup
      yup=colaup*cophup
      zup=siphup

      xum=silaum*cophum
      yum=colaum*cophum
      zum=siphum

      xvp=silavp*cophvp
      yvp=colavp*cophvp
      zvp=siphvp

      xvm=silavm*cophvm
      yvm=colavm*cophvm
      zvm=siphvm


      xzz=silazz*cophzz
      yzz=colazz*cophzz
      zzz=siphzz


      xee=silaee*cophee
      yee=colaee*cophee
      zee=siphee


      xss=silass*cophss
      yss=colass*cophss
      zss=siphss


      xsw=silasw*cophsw
      ysw=colasw*cophsw
      zsw=siphsw


      xse=silase*cophse
      yse=colase*cophse
      zse=siphse

      xne=silane*cophne
      yne=colane*cophne
      zne=siphne

!      dlyu(i,j)=erdrad*acos(xne*xse+yne*yse+zne*zse)
!      dlxv(i,j)=erdrad*acos(xsw*xse+ysw*yse+zsw*zse)
!      dlxu(i,j)=erdrad*acos(xzz*xee+yzz*yee+zzz*zee)
!      dlyv(i,j)=erdrad*acos(xzz*xss+yzz*yss+zzz*zss)
!      dlxp(i,j)=erdrad*acos(xup*xum+yup*yum+zup*zum)
!      dlyp(i,j)=erdrad*acos(xvp*xvm+yvp*yvm+zvp*zvm)


      dlyu(i,j)=max(1., erdrad*acos(min((xne*xse+yne*yse+zne*zse),1.)))
      dlxv(i,j)=max(1., erdrad*acos(min((xsw*xse+ysw*yse+zsw*zse),1.)))
      dlxu(i,j)=max(1., erdrad*acos(min((xzz*xee+yzz*yee+zzz*zee),1.)))
      dlyv(i,j)=max(1., erdrad*acos(min((xzz*xss+yzz*yss+zzz*zss),1.)))
      dlxp(i,j)=max(1., erdrad*acos(min((xup*xum+yup*yum+zup*zum),1.)))
      dlyp(i,j)=max(1., erdrad*acos(min((xvp*xvm+yvp*yvm+zvp*zvm),1.)))


      !print*,i,j

      !if(dlxv(i,j).lt.1.)write(6,*)'scheiss dlxv:',i,j,dlxv(i,j)         &
      !   ,dlyv(i,j),dlxp(i,j),dlxu(i,j),gila(2*i,2*j),giph(2*i,2*j)

      if(i.eq.25)print*,'pos',j,gila(2*i,2*j)*pigra,giph(2*i,2*j)*pigra, &
                              dlxp(i,j),dlyp(i,j)
      if(i.eq.25)print*,pigra*gila(2*i,2*j-1),pigra*gila(2*i,2*j+1),     &
                        pigra*giph(2*i,2*j-1),pigra*giph(2*i,2*j+1)
      if(depto(i,j).gt.1.) mende=mende+1
      ftwou(i,j)=2.*omega*sin(giph(2*i+1,2*j))
      ftwov(i,j)=2.*omega*sin(giph(2*i,2*j+1))
1478   continue


       rmindlxp=10E9
       rmindlyp=10E9
       DO I=1,IE
        DO J=1,JE
          if (depto(i,j).gt.0.5) then
          rmindlxp=min(dlxp(i,j),rmindlxp)
          rmindlyp=min(dlyp(i,j),rmindlyp)
          endif
        ENDDO
       ENDDO

       PRINT*,'min dlxp ',rmindlxp
       PRINT*,'min dlyp ',rmindlyp


       print*, ' sum of depth points ', mende
       write(82)0,507,-100,ie*je
       write(82)depto
       write(82)0,85,-100,ie*je
       write(82)dlxp
       write(82)0,185,-100,ie*je
       write(82)dlxu
       write(82)0,285,-100,ie*je
       write(82)dlxv
       write(82)0,86,-100,ie*je
       write(82)dlyp
       write(82)0,186,-100,ie*je
       write(82)dlyu
       write(82)0,286,-100,ie*je
       write(82)dlyv
       write(82)0,175,-100,ie*je
       write(82)ftwou
       write(82)0,176,-100,ie*je
       write(82)ftwov

  END SUBROUTINE ARCGRI

  ! ------------------------------------------------------------------------


   SUBROUTINE  INIT

      CHARACTER*256 :: ifilename_anta !(anta related)
      CHARACTER*256 :: ifilename !(anta related)
      CHARACTER*256 :: ofilename !(anta related)
      CHARACTER*10 :: lev_txt
      LOGICAL :: lex   ! file exists
      INTEGER :: IDAY
      INTEGER :: LDAY
      PARAMETER(me=360,ne=180,LE=33)
      INTEGER KE
      REAL :: GILA(2*IE,2*JE),GIPH(2*IE,2*JE)
      REAL :: GAULAT(me+2,ne+2),GAULON(me+2,ne+2)
      REAL(KIND=sp) :: DUMMY(me,ne)
      REAL :: LAND(me+2,ne+2)
      REAL :: GAUSSLA(ne),GAUSSLO(me)
      REAL :: HOPLAT(IE,JE),HOPLON(IE,JE)

      REAL :: S(ie,je),pt(ie,je)

      REAL :: FIIN(me+2,ne+2,1)
      REAL :: FIOUT(IE,JE,1)

      REAL :: TINTER(IE,JE,LE)

      REAL :: TEMIN(me+2,ne+2)

      REAL :: TEMOUT(IE,JE)

      REAL :: TOUTER(IE,JE)

      REAL :: SVAL

      REAL, ALLOCATABLE :: ZZOUT(:)


      ALLOCATE (ZZOUT(levels))

      SVAL = 0.
      KE = levels

      SELECT CASE (levels)

      CASE(3)
         ZZOUT =(/ 20.,20., 5000./)
      CASE(20)
         ZZOUT =(/ 20.,20., 20., 30.,40.,50.,70.       &
                 ,90.,120.,150.,180.,210.,250.,300.    &
             ,400.,500.,600.,700.,900.,1400./)
      CASE(22)
         ZZOUT =(/ 25,75,139,205,287,374,480,593,730,879,1055,   &
                 1249,1478,1731,2027,2357,2741,3171,3668,4229,   &
                 4874,5603 /)
      CASE(40)
         ZZOUT =(/ 12.,10.,10.,10.,10.,10.,13.,15.,20.,25.             &
                  ,30.,35.,40.,45.,50.,55.,60.,70.,80.,90.             &
                  ,100.,110.,120.,130.,140.,150.,170.,180.,190.,200.   &
                  ,220.,250.,270.,300.,350.,400.,450.,500.,500.,600./)

      CASE(80)
         ZZOUT =(/ 12.,10.,10.,10.,10.,10.,10.,11.,11.,12.,          &
                  13.,13.,14.,14.,15.,16.,16.,17.,18.,19.,           &
                  20.,21.,21.,22.,24.,25.,26.,27.,28.,29.,           &
                  31.,32.,34.,35.,37.,39.,40.,42.,44.,46.,           &
                  48.,50.,53.,55.,58.,60.,63.,66.,69.,72.,           &
                  76.,79.,83.,87.,91.,95.,99.,104.,108.,113.,        &
                  119.,124.,130.,136.,142.,149.,155.,163.,170.,178., &
                  186.,195.,204.,213.,223.,233.,244.,255.,267.,279./)

      CASE DEFAULT
         write(*,*) " no vertical levels definitions"
         STOP
      END SELECT


      sumtiefe=0.
      do i=1,levels
      tiefe=zzout(i)
      zzout(i)=sumtiefe+tiefe/2.
      sumtiefe=sumtiefe+tiefe
      print*,nint(sumtiefe),zzout(i)
      enddo


      DO I=1,me
        GAUSSLO(I)=-0.5+I
      ENDDO
      DO J=1,ne
        GAUSSLA(J)=90.5-j
      ENDDO

!     SOIE CONSTANTS
      PI=3.141592656
      GRARAD=180./PI

      DO I=1,me
       DO J=1,ne
        GAULAT(I+1,J+1)=GAUSSLA(J)
        GAULON(I+1,J+1)=GAUSSLO(I)
       ENDDO
      ENDDO

!     ZYKL.RAND OST/WEST
      DO J=1,ne+2
        GAULAT(1,J)=GAULAT(me+1,J)
        GAULAT(me+2,J)=GAULAT(2,J)
        GAULON(1,J)=GAULON(me+1,J)-360.
        GAULON(me+2,J)=GAULON(2,J)+360.
      ENDDO

!     POLE (SCHUMIELN: OBEN UND UNTEN EINE REIHE HINZU)
      DO I=1,me+2

          GAULAT(I,1)=180.-GAULAT(I,2)
          GAULAT(I,ne+2)=-180.-GAULAT(I,ne+1)
          GAULON(I,1)=GAULON(I,2)
          GAULON(I,ne+2)=GAULON(I,ne+1)

      ENDDO

!     NOCHMAL ZYKL.RAND OST/WEST
      DO J=2,ne+1
        GAULAT(1,J)=GAULAT(me+1,J)
        GAULAT(me+2,J)=GAULAT(2,J)
        GAULON(1,J)=GAULON(me+1,J)-360.
        GAULON(me+2,J)=GAULON(2,J)+360
      ENDDO

      DO I=1,me+2
       DO J=1,ne+2
        GAULAT(I,J)=GAULAT(I,J)/GRARAD
        GAULON(I,J)=GAULON(I,J)/GRARAD
       ENDDO
      ENDDO

!-----------------------------------------------------------------------------------
! ANNUAL OR MONTHLY CLIMATOLOGY:
!-----------------------------------------------------------------------------------
     IF (l_init_clim) THEN
        LDAY = 1
        write(*,*) " initial condition from annual climatology"
     ELSE
        LDAY = 12
        write(*,*) " initial condition from monthly climatology"
     ENDIF
!-----------------------------------------------------------------------------------
! TEMPERATURE
!-----------------------------------------------------------------------------------

     ifilename_anta = './output/'//TRIM(grid_name)//'_anta'

     INQUIRE (FILE=TRIM(ifilename_anta), EXIST=lex)
     IF (.NOT.lex) THEN
             WRITE(*,*) ': FILE DOES NOT EXIST (',TRIM(ifilename_anta),')'
        STOP
     END IF

#ifndef NOENDIANCONVERT
     OPEN(81,FILE=TRIM(ifilename_anta), &
           FORM='UNFORMATTED', convert='BIG_ENDIAN')
#else
  STOP 'ERROR: open(...,convert=BIG_ENDIAN,...) not supported!'
#endif
     REWIND(81)

     IF (LDAY == 1) THEN
        ifilename = TRIM(phc_files_path)//'/phc.temp00.98.p.ext4'
     ELSE
        ifilename = TRIM(phc_files_path)//'/phc.tempmon.98.p.ext4'
     ENDIF

     write(lev_txt,'(i2)') levels
     ofilename = './output/'//TRIM(grid_name)//'L'//TRIM(lev_txt)//'_INITEM_PHC'
#ifndef NOENDIANCONVERT
     OPEN(55,FILE=TRIM(ofilename),ACCESS='SEQUENTIAL', &
          FORM='UNFORMATTED',convert='BIG_ENDIAN')
#else
  STOP 'ERROR: open(...,convert=BIG_ENDIAN,...) not supported!'
#endif

     INQUIRE (FILE=TRIM(ifilename), EXIST=lex)
     IF (.NOT.lex) THEN
             WRITE(*,*) ': FILE DOES NOT EXIST (',TRIM(ifilename),')'
        STOP
     END IF
#ifndef NOENDIANCONVERT
     OPEN(20,FILE=TRIM(ifilename),ACCESS='SEQUENTIAL', &
             FORM='UNFORMATTED',convert='BIG_ENDIAN')
#else
  STOP 'ERROR: open(...,convert=BIG_ENDIAN,...) not supported!'
#endif
     REWIND(20)

     ! READ ANTA.....
     READ(81)IF1,IF2,IF3,IF4
     READ(81) GILA
     READ(81)IF1,IF2,IF3,IF4
     READ(81) GIPH


     DO IDAY=1,LDAY
     Print*,'DATUM ',IDAY
!     READ INPUT FILES
     IF3=1
     DO LEV=1,LE
     IF3UP=MAX(1,IF3)

     READ(20) IF1, IF2, IF3, IF4
     READ(20)((DUMMY(I,J),I=1,ME),J=NE,1,-1)

     DO I=1,me
       DO J=1,ne
         TEMIN(I+1,J+1)=DUMMY(I,J)
       ENDDO
     ENDDO

!     DATENFELDER VORRBEITEN
!     ZYKL.RAND OST/WEST

      DO J=1,ne+2
        TEMIN(1,J)=TEMIN(me+1,J)
        TEMIN(me+2,J)=TEMIN(2,J)
      ENDDO

!     POLE (EINE REIHE HINZU)
      DO I=1,me+2
          TEMIN(I,1)=TEMIN(I,2)
          TEMIN(I,ne+2)=TEMIN(I,ne+1)
      ENDDO

      PRINT*,'Datum ',iday,' EINLESEN UEBERLEBT'

!H    AN DIESER STELLE WIRD JEWEILS FUER U,V UND SKALARE GROESSEN
!H    INTERPOLIERT, UND ZWAR AM JEWEILS DAZUGEHOERENDEN GITTERPUNKT
!H    1. FUER SKALARE WIE GEHABT AN DRUCKPUNKT
      PRINT*,'SKALARE GROESSE INTERPOLIEREN'

       DO M=1,IE
         DO N=1,JE
           HOPLON(M,N)=GILA(2*M,2*N)
           HOPLAT(M,N)=GIPH(2*M,2*N)
           IF (HOPLON(M,N).GT.(2.*PI)) HOPLON(M,N)=HOPLON(M,N)-(2.*PI)
           IF (HOPLON(M,N).LT.0.) HOPLON(M,N)=HOPLON(M,N)+(2.*PI)
         ENDDO
       ENDDO


      DO I=1,me+2
        DO J=1,ne+2
            land(i,j)=0.
            IF(TEMIN(I,J).eq.999.)land(I,j)=1.
             IF(TEMIN(I,J).eq.999.) TEMIN(I,j)=15.
            FIIN(I,J,1)=TEMIN(I,J)
        ENDDO
      ENDDO
       CALL BLN2HOP(1,me,ne,IE,JE,GAULAT,GAULON,HOPLAT,HOPLON, &
                    FIIN,FIOUT,LAND,SVAL)

      DO M=1,IE
        DO N=1,JE
            TEMOUT(M,N)=FIOUT(M,N,1)
        ENDDO
      ENDDO

!h    PERIODISCHER RAND
      PRINT*,'PERIODISCHER RAND'
      DO N=1,JE
        TEMOUT(1,N)=TEMOUT(IE-1,N)
        TEMOUT(IE,N)=TEMOUT(2,N)
      ENDDO


!HH   INTERPLOATE TO LEVEL
      DO M=1,IE
        DO N=1,JE
          TINTER(M,N,LEV)=TEMOUT(M,N)
        ENDDO
      ENDDO

      DO KLEV=1,KE
        IF (nint(ZZOUT(KLEV)).GT.IF3UP        &
         .AND.NINT(ZZOUT(KLEV)).LE.IF3) THEN
          PRINT*,klev,IF3UP,ZZOUT(KLEV),IF3
          DZL=FLOAT(IF3-IF3UP)
          DZS=FLOAT(IF3)-ZZOUT(KLEV)
          ALPHA=DZS/DZL
          DO M=1,IE
            DO N=1,JE
              TOUTER(M,N)=(ALPHA*TINTER(M,N,LEV-1))         &
                              +((1.-ALPHA)*TINTER(M,N,LEV))
            ENDDO
          ENDDO
          WRITE(55) IDAY,IF2,NINT(ZZOUT(KLEV)),(ME*NE)
          WRITE(55) TOUTER
         ELSE
           IF (NINT(ZZOUT(KLEV)).GT.IF3.AND.LEV.EQ.LE)then
            ALPHA=0.
            DO M=1,IE
              DO N=1,JE
                TOUTER(M,N)=(ALPHA*TINTER(M,N,LEV-1))       &
                              +((1.-ALPHA)*TINTER(M,N,LEV))
              ENDDO
            ENDDO
            WRITE(55) IDAY,IF2,NINT(ZZOUT(KLEV)),(ME*NE)
            WRITE(55) TOUTER
           ENDIF
         ENDIF
      ENDDO
      IF3UP=IF3

      ENDDO

      ENDDO

      CLOSE(55)
      CLOSE(20)
      CLOSE(81)

!-----------------------------------------------------------------------------------
! SALINITY
!-----------------------------------------------------------------------------------

     ifilename_anta = './output/'//TRIM(grid_name)//'_anta'

     INQUIRE (FILE=TRIM(ifilename_anta), EXIST=lex)
     IF (.NOT.lex) THEN
             WRITE(*,*) ': FILE DOES NOT EXIST (',TRIM(ifilename_anta),')'
        STOP
     END IF

#ifndef NOENDIANCONVERT
     OPEN(81,FILE=TRIM(ifilename_anta), &
           FORM='UNFORMATTED', convert='BIG_ENDIAN')
#else
  STOP 'ERROR: open(...,convert=BIG_ENDIAN,...) not supported!'
#endif
     REWIND(81)

     IF (LDAY == 1) THEN
        ifilename = TRIM(phc_files_path)//'/phc.salt00.98.p.ext4'
     ELSE
        ifilename = TRIM(phc_files_path)//'/phc.saltmon.98.p.ext4'
     ENDIF

     write(lev_txt,'(i2)') levels
     ofilename = './output/'//TRIM(grid_name)//'L'//TRIM(lev_txt)//'_INISAL_PHC'
#ifndef NOENDIANCONVERT
     OPEN(55,FILE=TRIM(ofilename),ACCESS='SEQUENTIAL', &
          FORM='UNFORMATTED',convert='BIG_ENDIAN')
#else
  STOP 'ERROR: open(...,convert=BIG_ENDIAN,...) not supported!'
#endif

     INQUIRE (FILE=TRIM(ifilename), EXIST=lex)
     IF (.NOT.lex) THEN
             WRITE(*,*) ': FILE DOES NOT EXIST (',TRIM(ifilename),')'
        STOP
     END IF
#ifndef NOENDIANCONVERT
     OPEN(20,FILE=TRIM(ifilename),ACCESS='SEQUENTIAL',   &
          FORM='UNFORMATTED',convert='BIG_ENDIAN')
#else
  STOP 'ERROR: open(...,convert=BIG_ENDIAN,...) not supported!'
#endif
     REWIND(20)

     ! READ ANTA.....
     READ(81)IF1,IF2,IF3,IF4
     READ(81) GILA
     READ(81)IF1,IF2,IF3,IF4
     READ(81) GIPH


     DO IDAY=1,LDAY
     Print*,'DATUM ',IDAY
!     READ INPUT FILES
     IF3=1
     DO LEV=1,LE
     IF3UP=MAX(1,IF3)

     READ(20)IF1,IF2,IF3,IF4
     READ(20)((DUMMY(I,J),I=1,me),J=ne,1,-1)

     DO I=1,me
       DO J=1,ne
         TEMIN(I+1,J+1)=DUMMY(I,J)
       ENDDO
     ENDDO

!     DATENFELDER VORRBEITEN
!     ZYKL.RAND OST/WEST

      DO J=1,ne+2
        TEMIN(1,J)=TEMIN(me+1,J)
        TEMIN(me+2,J)=TEMIN(2,J)
      ENDDO

!     POLE (EINE REIHE HINZU)
      DO I=1,me+2
          TEMIN(I,1)=TEMIN(I,2)
          TEMIN(I,ne+2)=TEMIN(I,ne+1)
      ENDDO

      PRINT*,'Datum ',iday,' EINLESEN UEBERLEBT'

!H    AN DIESER STELLE WIRD JEWEILS FUER U,V UND SKALARE GROESSEN
!H    INTERPOLIERT, UND ZWAR AM JEWEILS DAZUGEHOERENDEN GITTERPUNKT
!H    1. FUER SKALARE WIE GEHABT AN DRUCKPUNKT
      PRINT*,'SKALARE GROESSE INTERPOLIEREN'

       DO M=1,IE
         DO N=1,JE
           HOPLON(M,N)=GILA(2*M,2*N)
           HOPLAT(M,N)=GIPH(2*M,2*N)
           IF (HOPLON(M,N).GT.(2.*PI)) HOPLON(M,N)=HOPLON(M,N)-(2.*PI)
           IF (HOPLON(M,N).LT.0.) HOPLON(M,N)=HOPLON(M,N)+(2.*PI)
         ENDDO
       ENDDO


      DO I=1,me+2
        DO J=1,ne+2
            land(i,j)=0.
            IF(TEMIN(I,J).eq.999.)land(I,j)=1.
             IF(TEMIN(I,J).eq.999.) TEMIN(I,j)=15.
            FIIN(I,J,1)=TEMIN(I,J)
        ENDDO
      ENDDO
       CALL BLN2HOP(1,me,ne,IE,JE,GAULAT,GAULON,HOPLAT,HOPLON, &
                    FIIN,FIOUT,LAND,SVAL)

      DO M=1,IE
        DO N=1,JE
            TEMOUT(M,N)=FIOUT(M,N,1)
        ENDDO
      ENDDO

!h    PERIODISCHER RAND
      PRINT*,'PERIODISCHER RAND'
      DO N=1,JE
        TEMOUT(1,N)=TEMOUT(IE-1,N)
        TEMOUT(IE,N)=TEMOUT(2,N)
      ENDDO


!HH   INTERPLOATE TO LEVEL
      DO M=1,IE
        DO N=1,JE
          TINTER(M,N,LEV)=TEMOUT(M,N)
        ENDDO
      ENDDO

      DO KLEV=1,KE
        IF (nint(ZZOUT(KLEV)).GT.IF3UP        &
         .AND.NINT(ZZOUT(KLEV)).LE.IF3) THEN
          PRINT*,klev,IF3UP,ZZOUT(KLEV),IF3
          DZL=FLOAT(IF3-IF3UP)
          DZS=FLOAT(IF3)-ZZOUT(KLEV)
          ALPHA=DZS/DZL
          DO M=1,IE
            DO N=1,JE
              TOUTER(M,N)=(ALPHA*TINTER(M,N,LEV-1))         &
                              +((1.-ALPHA)*TINTER(M,N,LEV))
            ENDDO
          ENDDO
          WRITE(55) IDAY,IF2,NINT(ZZOUT(KLEV)),(ME*NE)
          WRITE(55) TOUTER
         ELSE
           IF (NINT(ZZOUT(KLEV)).GT.IF3.AND.LEV.EQ.LE)then
            ALPHA=0.
            DO M=1,IE
              DO N=1,JE
                TOUTER(M,N)=(ALPHA*TINTER(M,N,LEV-1))       &
                              +((1.-ALPHA)*TINTER(M,N,LEV))
              ENDDO
            ENDDO
            WRITE(55) IDAY,IF2,NINT(ZZOUT(KLEV)),(ME*NE)
            WRITE(55) TOUTER
           ENDIF
         ENDIF
      ENDDO
      IF3UP=IF3

      ENDDO

      ENDDO

      CLOSE(55)
      CLOSE(20)
      CLOSE(81)

      DEALLOCATE(ZZOUT)

  END SUBROUTINE INIT

!H    ********************************************************

      SUBROUTINE BLN2HOP(NUMFI,me,ne,IE,JE,ALAT2,ALON2,ALAT1,ALON1, &
                         SOE,SREG,LAND,SVAL)

!     BILINEAR INTERPOLATION FROM ONE GRID(me,ne) TO ANOTHER GRID(IE,JE)

!      PARAIETER(NB=3,MB=3)

     REAL ::  SOE(me+2,ne+2,NUMFI)
     REAL ::  SREG(IE,JE,NUMFI),G(me+2,ne+2)
     REAL ::  ALAT2(me+2,ne+2),ALON2(me+2,ne+2)
     REAL ::  ALAT1(IE,JE),ALON1(IE,JE)
     REAL ::  SVAL,ALPHA,BETA
     REAL ::  HX(me+2,ne+2,NUMFI),HHX(me+2,ne+2,NUMFI)
     REAL ::  LAND(me+2,ne+2)
     REAL ::  HG(me+2,ne+2)

!     SOIE CONSTANTS
      PI=3.141592656
      GRARAD=180./PI


!     DIFFUSION INTO LAND (NEW 12/99)

      DO L=1,NUMFI
        DO J=1,ne+2
          DO I=1,me+2
            HHX(I,J,L)=SOE(I,J,L)
!HH      X                 *((LAND(I,J)-1.))*(-1.)
          ENDDO
        ENDDO
      ENDDO

      DO J=1,ne+2
        DO I=1,me+2
         HG(I,J)=1.
         IF (LAND(I,J).ge.0.5) HG(I,J)=1.E-12
        ENDDO
      ENDDO

        DO ITER=1,300
          DO J=1,ne+2
            DO I=1,me+2
              G(I,J)=HG(I,J)
              DO L=1,NUMFI
                HX(I,J,l)=HHX(I,J,l)
              ENDDO
            ENDDO
          ENDDO
          DO J=1,ne+2
            JO=MAX(J-1,1)
            JU=MIN(J+1,ne+2)
            DO I=1,me+2
              IL=I-1
              IF(IL.LT.1)IL=me
              IR=I+1
              IF(IR.GT.me+2)IR=3
              RSUMG=0
              IF(LAND(I,J).GE.0.5)THEN
                RSUMG=(4.*G(I,J)                           &
                      +G(IL,J)+G(IR,J)+G(I,JO)+G(I,JU))/8.

                HG(I,J)=MIN(RSUMG,0.125)

                DO L=1,NUMFI
                  HHX(I,J,L)=(4.*HX(I,J,L)*G(I,J)         &
                                +HX(IL,J,L)*G(IL,J)       &
                                +HX(IR,J,L)*G(IR,J)       &
                                +HX(I,JO,L)*G(I,JO)       &
                                +HX(I,JU,L)*G(I,JU))/8.

                  HHX(I,J,L)=HHX(I,J,L)/RSUMG
                ENDDO
              ENDIF
            ENDDO
          ENDDO
          nland=0.
          Do I=1,me
            DO j=1,ne
              if (HG(i,j).le.2.E-12) nland=nland+1
            ENDDO
          ENDDO
!          PRINT*,iter,' Nland: ', nland
        ENDDO


        DO J=1,ne+2
         DO I=1,me+2
          G(i,j)=HG(i,j)
          DO L=1,NUMFI
            SOE(I,J,L)=HHX(I,J,L)
            If (ABS(SOE(I,J,L)).LT.1.E-12) SOE(I,J,L)=0.
          ENDDO
         ENDDO
        ENDDO

        DO M=1,IE
        DO N=1,JE

!        PUNKT RECHTS OBEN
        DO J=1,ne+2
          IF (ALAT2(2,J).GE.ALAT1(M,N)) JO=J
        ENDDO
        DO I=1,me+2
          IF (ALON2(I,2).LE.ALON1(M,N)) IL=I
        ENDDO

!       IL=NINT(ALON1(M,N)*DLAMDAI+0.499999999)+1

!            PRINT*,'WE ARE AT THE MINIMUM ',M,N,ALAT1(M,N)*GRARAD
!     X                                    ,IL,JO,ALAT2(IL,JO)*GRARAD

!       PUNKT RECHTS OBEN --> LINKS UNTEN

        JLU=JO+1
        ILU=IL

        WWWALPHA=ALAT2(ILU,JLU)-ALAT1(M,N)
        WWWBETA=ALON2(ILU,JLU)-ALON1(M,N)
        IF(WWWBETA.GE.PI) WWWBETA=WWWBETA-2.*PI
        IF(WWWBETA.LE.-PI) WWWBETA=WWWBETA+2.*PI

        ALPHA=(WWWALPHA)/(ALAT2(ILU,JLU)-ALAT2(ILU,JLU-1))
        BETA=(WWWBETA)/(ALON2(ILU,JLU)-ALON2(ILU+1,JLU))
        DO I=1,NUMFI
           SREG(M,N,I)=ALPHA*BETA*SOE(ILU+1,JLU-1,I)*G(ILU+1,JLU-1)  &
                    +(1.-ALPHA)*(1.-BETA)*SOE(ILU,JLU,I)*G(ILU,JLU)  &
                    +(1.-ALPHA)*(BETA)*SOE(ILU+1,JLU,I)*G(ILU+1,JLU) &
                    +(ALPHA)*(1.-BETA)*SOE(ILU,JLU-1,I)*G(ILU,JLU-1)
           SREG(M,N,I)=SREG(M,N,I)/                    &
                    (ALPHA*BETA*G(ILU+1,JLU-1)         &
                    +(1.-ALPHA)*(1.-BETA)*G(ILU,JLU)   &
                    +(1.-ALPHA)*(BETA)*G(ILU+1,JLU)    &
                    +(ALPHA)*(1.-BETA)*G(ILU,JLU-1))
         ENDDO
!         *************************************************
      ENDDO
      ENDDO

      END  SUBROUTINE



! FOR TRIPOLAR GRID
!  ! ------------------------------------------------------------------------
!
!  SUBROUTINE cgri
!      parameter(np=${NP},np1=np+1,ie=4*np,ie1=ie+2,je=2*np,it8=2*ie1
!     &  ,ito=2*ie1,jf=${ne},jto=2*jf,np2=2*np,np8=8*np,np88=np8+8
!     &  ,jem=je+1
!     &  ,jep=je+2,je6=je+6)
!
!      real xk,yk,xh,yh,xx,yy,xz,yz,sx,sy
!      dimension dphieq(ito),x8(ito,jto),y8(ito,jto)
!      dimension sinat(629),cosat(629),elrad(jep)
!      dimension xend(160),yend(160),sx(7),sy(7)
!      dimension giph(ito,jto),gila(ito,jto),giltes(ito,3)
!      dimension xk(ito,jto),yk(ito,jto),feld(4320,2160)
!      dimension xh(ito,jto),yh(ito,jto)
!
!      integer*8 ibla(4)
!
!      real depto(ie1,jf)
!      reflat=55.
!      reflon=1.45
!      ibla(1)=${DATE}
!      ibla(2)=${NP}
!      ibla(3)=1
!      ibla(4)=ito*jto
!
!      print*,'np88=',np88,ie1,it8,jto,it8*jto
!
!      open(12,file=TRIM(files_path)//'TOTATO',form='formatted')
!      open(13,file='./output/anta',form='unformatted')
!      read(12,1200) feld
! 1200 format(12f6.0)
!
!
!      pi=4.*atan(1.)
!      pih=0.5*pi
!      topfa=pi/2160.
!      topfai=1./topfa
!      dphih=0.5*pih/np
!      elp2=2./(np+0.5)
!      elp1=1./(np+0.5)
!      zwop=2.*pi
!      bogra=180./pi
!      print*,'latitude of singularity',reflat,reflon*bogra
!      poldi=90.-reflat
!      absch=2.*tan(poldi*pi/360.)
!      print*,'absch=',absch
!      aa=2.-absch
!      bb=absch
!      do i=1,629
!         arg=0.01*i
!         sinat(i)=dsin(arg)
!         cosat(i)=dcos(arg)
!      enddo
!      do irr=1,jep
!         rr=0.5*(irr-1)
!         criac=1.-1.e-10
!         ra=min(elp2*rr/(1.+elp1*rr),criac)
!
!         elrad(irr)=ra
!         a=2.-ra*aa
!         b=max(2.-ra*aa-bb*ra**2,1.e-9)
!         print*,'axes',irr,ra,a,b
!         ip=3
!         do i=1,629
!            xx=a*cosat(i)
!            yy=b*sinat(i)
!!     call plot(xx,yy,ip)
!            ip=2
!         enddo
!      enddo
!      print*,elrad
!      auft=4.*np
!C     construction of orthogonals to ellipses
!      do irad=1,ito
!C     do irad=35,35
!         iant=irad
!         arg=pi*(irad-3)/auft
!         x=2.*dcos(arg)
!         y=2.*dsin(arg)
!         xx=x
!         yy=y
!         xk(iant,jep)=xx
!         yk(iant,jep)=yy
!         x8(iant,jep)=x
!         y8(iant,jep)=y
!!         if(irad.eq.3.or.irad.eq.403)print*,'ax2',irad,x,y
!
!         gila(iant,jep)=arg
!         str=irad
!         print*,'start',irad,arg,x,y
!         zerr=1.05+0.002*irad
!C     call number(zerr*x,zerr*y,0.05,str,0.,-1)
!!     call plot(xx,yy,3)
!         ra=0.
!
!C     do 200 l=1,jem-1
!         tw=2.
!         do 200 l=1,jem-1
!            jant=jep-l
!            ziel=elrad(l+1)
!
!            ifind=0
!            dra=0.000001
!            if(l.gt.5)dra=0.00001
!            do 101 idr=1,100000
!               reff=1.
!               ralt=ra
!
!               do kitt=1,25
!
!                  f=x**2/(tw-ra*aa)**2 + y**2/(tw-ra*aa-bb*ra**2)**2-1.
!                  dfdr=-tw*x**2/(tw-ra*aa)**3
!     1                 -tw*y**2*(aa+tw*bb*ra)/(tw-ra*aa-bb*ra**2)**3
!
!                  ra=ra+f/dfdr
!               enddo
!
!               if(ra.ge.ziel) then
!                  ifind=1
!                  reff=(ziel-ralt)/(ra-ralt)
!                  ra=ziel
!               endif
!               fx=2.*x/(2.-aa*ra)**2
!               fy=2.*y/(2.-aa*ra-bb*ra**2)**2
!               yst=-fx/fy
!               dx=  fx
!               dy=  fy
!C     print*,'dx',dx,dy
!               betra=sqrt(dx**2+dy**2)
!               r=1.
!               dx=dra*dx/betra
!               dy=dra*dy/betra
!               if(dy*y.gt.0.9*y**2) r=0.9*y/dy
!               x=x-dx*r*reff
!               y=y-dy*r*reff
!               xx=x
!               yy=y
!               if(ifind.eq.1) go to 201
! 101        continue
! 201        continue
!            if(l.eq.jem)print*,'?',irad,x,y
!            if(abs(y).lt.1.e-6)y=0.
!
!            xk(iant,jant)=x
!            yk(iant,jant)=y
!            x8(iant,jant)=x
!            y8(iant,jant)=y
!            if(irad.eq.3.or.irad.eq.403)print*,'ax2',irad,x,y
!
!            gila(iant,jant)=atan2(y,x)
!            radi=sqrt(x**2+y**2)
!            poldi=2.*atan(0.5*radi)
!            giph(iant,jant)=pih-poldi
!
! 200     continue
! 99      continue
!      enddo
!      print*,'zum ersten'
!!      do j=1,jep
!!         print*,'xk',j,yk(44,j),xk(44,j)**2+yk(44,j)**2
!!      enddo
!!      do i=1,50
!!         print*,i,xk(i,jep),yk(i,jep),xk(i,jep)**2+yk(i,jep)**2
!!      enddo
!!      print*,'yk'
!!      print*,(i,yk(i,42),i=1,ito)
!
!      do i=1,ito
!         xk(i,1)=xk(i,2)
!         yk(i,1)=0.
!         x8(i,1)=x8(i,2)
!         y8(i,1)=0.
!      enddo
!
!      do j=1,jto
!         do i=1,ito
!            xh(i,j)=x8(i,j)
!            yh(i,j)=y8(i,j)
!         enddo
!      enddo
!      sipia=2.*sin(0.15*pi)
!      ipr1=3
!      ipr2=4*np+3
!      ipr3=8*np+3
!      idim=0.6*np
!      do i=1,ito
!C     if(abs(yk(i,je)).lt.sipia) then
!         jmid=min(abs(i-ipr1),abs(i-ipr2),abs(i-ipr3))
!         if(jmid.le.idim)then
!            jfi=min(np/5,idim-jmid)
!            jd=2*jfi+1
!            print*,'jd0:',jfi,np/5,idim,jmid,idim-jmid
!            do j=0,jfi-1
!               print*,'jd1:',jd-2*j-1
!               xk(i,jd-2*j-1)=0.5*(xh(i,jd-j)+xh(i,jd-(j+1)))
!               yk(i,jd-2*j-1)=0.5*(yh(i,jd-j)+yh(i,jd-(j+1)))
!             print*,'jd2:',jd-2*j-2
!               xk(i,jd-2*j-2)=xh(i,jd-j-1)
!               yk(i,jd-2*j-2)=yh(i,jd-j-1)
!               x8(i,jd-2*j-1)=0.5*(xh(i,jd-j)+xh(i,jd-(j+1)))
!               y8(i,jd-2*j-1)=0.5*(yh(i,jd-j)+yh(i,jd-(j+1)))
!               x8(i,jd-2*j-2)=xh(i,jd-j-1)
!               y8(i,jd-2*j-2)=yh(i,jd-j-1)
!            enddo
!         endif
!      enddo
!
!      do j=jep,1,-1
!         do i=1,ito
!            x=xk(i,j)
!            y=yk(i,j)
!            xk(i,j+4)=x
!            yk(i,j+4)=y
!            gila(i,j+4)=atan2(y,x)
!            radi=sqrt(x**2+y**2)
!            poldi=2.*atan(0.5*radi)
!            giph(i,j+4)=pih-poldi
!         enddo
!      enddo
!      print*,'gila'
!
!
!
! 622  format(6e10.3)
!
!
!      do j=1,4
!      do i=1,ito
!      yk(i,5-j)=2.*yk(i,5)-yk(i,5+j)
!      xk(i,5-j)=2.*xk(i,5)-xk(i,5+j)
!!      y8(i,5-j)=2.*y8(i,5)-y8(i,5+j)
!!      x8(i,5-j)=2.*x8(i,5)-x8(i,5+j)
!      enddo
!      enddo
!
!!!      do j=1,5
!      do j=1,4
!      do i=1,ito
!      gila(i,j)=atan2(yk(i,j),xk(i,j))
!      radi=sqrt(xk(i,j)**2+yk(i,j)**2)
!      poldi=2.*atan(0.5*radi)
!      giph(i,j)=pih-poldi
!      enddo
!      enddo
!
!
! 601  format('wo',i3,7e10.3)
!      dphi=pi/ie
!
!      phiq=0.
!      do i=1,ito
!         dphieq(i)=giph(i,je6-1)-giph(i,je6)
!         phiq=phiq+dphieq(i)
!      enddo
!      phiq=phiq/ito
!      m=0
!      do j=je6,jto-1
!         m=m+1
!         do i=1,ito
!            gila(i,j+1)=gila(i,j)
!            dphi=(40*dphieq(i)+phiq*m)/(m+40)
!            cophi=cos(giph(i,j))
!            giph(i,j+1)=giph(i,j)-dphi
!         enddo
!      enddo
!
!      print*,'final rotation'
!      do j=1,jto
!         do i=1,ito
!            gila(i,j)=mod(gila(i,j)+reflon+pi,zwop)-pi
!         enddo
!      enddo
!
!!HH   Here we need to set the longitudes of the northpole grid point
!!     to avoid problems with remapcon interpolation
!
!      gila(2*np+3,5)=gila(2*np+3,6)
!      gila(6*np+3,5)=gila(6*np+3,6)
!
!
!      print*,'np1-1',np1-1,gila(np1-1,5)*bogra,GILA(np1-1,6)*bogra
!      print*,'np1  ',np1,gila(np1,5)*bogra,GILA(np1,6)*bogra
!      print*,'np1+1',np1+2,gila(np1+2,5)*bogra,GILA(np1+2,6)*bogra
!
!      ibla(2)=54
!      write(13)ibla
!      write(13)gila
!      ibla(2)=55
!      write(13)ibla
!      write(13)giph
!
!      print*,'suchen'
!
!      do j=2,jf-1
!         do i=2,ie1-1
!
!            elev=0.
!            deep=0.
!            nel=0
!            nde=0
!
!            do jj=1,10
!               do ii=1,10
!                  al=0.1*ii-0.05
!                  ga=0.1*jj-0.05
!                  be=1.-al
!                  de=1.-ga
!                  gril=topfai*(al*gila(2*i+1,2*j)+be*gila(2*i-1,2*j))
!                  grilr=gila(2*i+1,2*j)
!                  grill=gila(2*i-1,2*j)
!                  if(abs(grilr-grill).gt.2.)grill=grilr
!                  gril=(al*grilr+be*grill)*topfai
!
!                  griph=topfai*(ga*giph(2*i,2*j+1)+de*giph(2*i,2*j-1))
!                  iet5=int(mod(gril+4319.0,4320.0)+1.0)
!                  jet5=min(1080-int(griph),2160)
!!                  print*,'jet5',jet5,1080,griph
!
!                  if(feld(iet5,jet5).gt.0.)elev=elev+feld(iet5,jet5)
!                  if(feld(iet5,jet5).le.0.)deep=deep-feld(iet5,jet5)
!                  if(feld(iet5,jet5).gt.0.)nel=nel+1
!                  if(feld(iet5,jet5).le.0.)nde=nde+1
!
!               enddo
!            enddo
!
!            if(nde.ge.nel.and.nde.ne.0) depto(i,j)=deep/nde
!            if(nel.gt.nde) depto(i,j)=0.
!         enddo
!      enddo
!
!      print*,'nach besetzung'
!      ibla(4)=ie1*jf
!      ibla(2)=84
!      write(13)ibla
!      write(13)depto
!      print*,'rundplot'
!      open(33,file='./output/topo',form='formatted')
!      do ii=1,ie1-21,20
!         write(33,*)'streifen',ii
!         do j=1,jf
!            write(33,3300)j,(depto(ii+i,j),i=1,20)
!         enddo
!      enddo
! 3300 format(i5,20f6.0)
!
!  END SUBROUTINE cgri

SUBROUTINE anta2nc

  IMPLICIT NONE

  REAL :: DEPTO(IE,JE)

  CHARACTER*256 :: ifilename, ofilename

  LOGICAL :: lex
  INTEGER :: i1, i2, i3, i4
  INTEGER :: ipoints, flen

  INTEGER       ::   i,j, ii, jj         ! looping indicees
  INTEGER       ::   ip1,im1             ! i+1, i-1
  REAL          ::   pi                  ! PI
  REAL          ::   rad2deg             ! 180/PI
  REAL, ALLOCATABLE    ::   lon(:,:)     ! (2*ie,2*je) longitudes of doubled array
  REAL, ALLOCATABLE    ::   lat(:,:)     ! (2*ie,2*je) latitudes of doubled array
  REAL, ALLOCATABLE    ::   deuto(:,:)   ! MASK U
  REAL, ALLOCATABLE    ::   deute(:,:)   ! MASK V
  REAL, ALLOCATABLE    ::   lons(:,:)    ! (ie-2,je) longitudes of scalars
  REAL, ALLOCATABLE    ::   lats(:,:)    ! (ie-2,je) latitudes of scalars
  REAL, ALLOCATABLE    ::   masks(:,:)   ! (ie-2,je) inverse land sea mask
  REAL, ALLOCATABLE    ::   lonu(:,:)    ! (ie-2,je) longitudes of vector u-components
  REAL, ALLOCATABLE    ::   latu(:,:)    ! (ie-2,je) latitudes of vector u-components
  REAL, ALLOCATABLE    ::   masku(:,:)   ! (ie-2,je) inverse land sea mask
  REAL, ALLOCATABLE    ::   lonv(:,:)    ! (ie-2,je) longitudes of vector v-components
  REAL, ALLOCATABLE    ::   latv(:,:)    ! (ie-2,je) latitudes of vector v-components
  REAL, ALLOCATABLE    ::   maskv(:,:)   ! (ie-2,je) inverse land sea mask
  REAL, ALLOCATABLE    ::   lonc(:,:)    ! (ie-2,je) corner grid points of scalar grid
  REAL, ALLOCATABLE    ::   latc(:,:)    ! (ie-2,je) corner grid points of scalar grid
  REAL, ALLOCATABLE    ::   clons(:,:,:) ! (4,ie-2,je) corner longitudes of scalars
  REAL, ALLOCATABLE    ::   clats(:,:,:) ! (4,ie-2,je) corner latitudes of scalars
  REAL, ALLOCATABLE    ::   clonu(:,:,:) ! (4,ie-2,je) corner longitudes of vector u-components
  REAL, ALLOCATABLE    ::   clatu(:,:,:) ! (4,ie-2,je) corner latitudes of vector u-components
  REAL, ALLOCATABLE    ::   clonv(:,:,:) ! (4,ie-2,je) corner longitudes of vector v-components
  REAL, ALLOCATABLE    ::   clatv(:,:,:) ! (4,ie-2,je) corner latitudes of vector v-components

  INTEGER :: fileid, nvid, notsid, ncartid

  ! input filename
  ifilename = './output/'//TRIM(grid_name)//'_anta'
  INQUIRE (FILE=TRIM(ifilename), EXIST=lex)
  IF ( .NOT. lex ) THEN
    WRITE(0,*) 'Could not open file <', ifilename, '>'
    STOP 'Open failed'
  END IF

  ! output filename
  ofilename =  './output/'//TRIM(grid_name)

#ifndef NOENDIANCONVERT
  OPEN(14, FILE=ifilename,convert='BIG_ENDIAN', FORM='UNFORMATTED')
#else
  STOP 'ERROR: open(...,convert=BIG_ENDIAN,...) not supported!'
#endif
  REWIND(14)

  READ(14) i1,i2,i3,i4

  ipoints = i4 / 4

  WRITE(*,*) 'The input file has ', ipoints, ' gridpoints.'

  IF ( ie*je .NE. ipoints ) THEN
    STOP 'Number of longitudes and latitudes not match!'
  END IF

  WRITE(*,*) 'o.k.'

  ALLOCATE(lon(2*ie,2*je))
  ALLOCATE(lat(2*ie,2*je))
  ALLOCATE(deuto(ie,je))
  ALLOCATE(deute(ie,je))
  ALLOCATE(lons(ie,je))
  ALLOCATE(lats(ie,je))
  ALLOCATE(masks(ie,je))
  ALLOCATE(lonu(ie,je))
  ALLOCATE(latu(ie,je))
  ALLOCATE(masku(ie,je))
  ALLOCATE(lonv(ie,je))
  ALLOCATE(latv(ie,je))
  ALLOCATE(maskv(ie,je))
  ALLOCATE(lonc(ie,je))
  ALLOCATE(latc(ie,je))
  ALLOCATE(clons(4,ie,je))
  ALLOCATE(clats(4,ie,je))
  ALLOCATE(clonu(4,ie,je))
  ALLOCATE(clatu(4,ie,je))
  ALLOCATE(clonv(4,ie,je))
  ALLOCATE(clatv(4,ie,je))

  READ(14) lon
  READ(14) i1,i2,i3,i4
  READ(14) lat
  READ(14) i1,i2,i3,i4
  READ(14) depto

  CLOSE(14)

  !
  !-- create arrays of longitudes and latitudes
  !
  !--  follow the OASIS conventions: S -> N
  !
!!  DO j = 1, 2*je
!!    lon(:,j)=gila(:,2*je+1-j)
!!    lat(:,j)=giph(:,2*je+1-j)
!!  ENDDO
  !
  !--  convert from radiant to degree
  !
  pi = ATAN(1.)*4.0
  rad2deg = 180./pi
  lon(:,:) = lon(:,:) * rad2deg
  lat(:,:) = lat(:,:) * rad2deg

  WHERE (lon(:,:) < 0.)
    lon(:,:) = lon(:,:) + 360.
  END WHERE

  !
  !--  extract scalar/vector grid points
  !
  !          1  2  3  4  5  6  ... 2*ie
  !
  !      1   c  v  c  v  c  v
  !      2   u  s  u  s  u  s
  !      3   c  v  c  v  c  v
  !      4   u  s  u  s  u  s
  !      5   c  v  c  v  c  v           c: grid cell corners of scalars
  !      6   u  s  u  s  u  s           v: vector-v
  !      :                              u: vector-u
  !      :                              s: scalar
  !     2*ij
  !

          DO i = 1, ie
             DO j = 1, je
     !--     scalar
                lats(i,j) = lat(i*2,j*2)
                lons(i,j) = lon(i*2,j*2)
     !--     vector - u
                latu(i,j) = lat(i*2-1,j*2)
                lonu(i,j) = lon(i*2-1,j*2)
     !--     vector - v
                latv(i,j) = lat(i*2,j*2-1)
                lonv(i,j) = lon(i*2,j*2-1)
     !--     corners of scalar grid cells
                latc(i,j) = lat(i*2-1,j*2-1)
                lonc(i,j) = lon(i*2-1,j*2-1)
             ENDDO
          ENDDO
     !
     !--  create corner arrays for SCRIP interpolation
     !
          DO i = 1, ie
     !
             im1 = i - 1
             ip1 = i + 1
             IF (im1 == 0) im1 = ie-2
             IF (ip1 == ie+1) ip1 = 3

!-------------------------------------------------
! REVERSED  to  ORIGINAL from mo_coupling.f90
!------------------------------------------------
     !
     !--     scalar
     !
             DO j = 1, je-1
                clons(1,i,j) = lonc(i  ,j)
                clons(2,i,j) = lonc(i  ,j+1)
                clons(3,i,j) = lonc(ip1,j+1)
                clons(4,i,j) = lonc(ip1,j)
                clats(1,i,j) = latc(i  ,j)
                clats(2,i,j) = latc(i  ,j+1)
                clats(3,i,j) = latc(ip1,j+1)
                clats(4,i,j) = latc(ip1,j)
             ENDDO
             clons(1,i,je) = lonc     (i  ,je)
             clons(2,i,je) = lonu(i  ,je)
             clons(3,i,je) = lonu(ip1,je)
             clons(4,i,je) = lonc     (ip1,je)
             clats(1,i,je) = latc     (i  ,je)
             clats(2,i,je) = latu(i  ,je)
             clats(3,i,je) = latu(ip1,je)
             clats(4,i,je) = latc     (ip1,je)

     !
     !--     vector - u
     !
             DO j = 1, je-1
                clonu(1,i,j) = lonv(im1,j  )
                clonu(2,i,j) = lonv(im1,j+1)
                clonu(3,i,j) = lonv(i  ,j+1)
                clonu(4,i,j) = lonv(i  ,j  )
                clatu(1,i,j) = latv(im1,j  )
                clatu(2,i,j) = latv(im1,j+1)
                clatu(3,i,j) = latv(i  ,j+1)
                clatu(4,i,j) = latv(i  ,j  )
             ENDDO
             clonu(1,i,je) = lonv(im1,je)
             clonu(2,i,je) = lons(im1,je)
             clonu(3,i,je) = lons(i  ,je)
             clonu(4,i,je) = lonv(i  ,je)
             clatu(1,i,je) = latv(im1,je)
             clatu(2,i,je) = lats(im1,je)
             clatu(3,i,je) = lats(i  ,je)
             clatu(4,i,je) = latv(i  ,je)
     !
     !--     vector - v
     !
             DO j = 2, je
                clonv(1,i,j) = lonu(i  ,j-1)
                clonv(2,i,j) = lonu(i  ,j  )
                clonv(3,i,j) = lonu(ip1,j  )
                clonv(4,i,j) = lonu(ip1,j-1)
                clatv(1,i,j) = latu(i  ,j-1)
                clatv(2,i,j) = latu(i  ,j  )
                clatv(3,i,j) = latu(ip1,j  )
                clatv(4,i,j) = latu(ip1,j-1)
             ENDDO
             clonv(1,i,1) = lonc     (i  ,1)
             clonv(2,i,1) = lonu(i  ,1)
             clonv(3,i,1) = lonu(ip1,1)
             clonv(4,i,1) = lonc     (ip1,1)
             clatv(1,i,1) = latc     (i  ,1)
             clatv(2,i,1) = latu(i  ,1)
             clatv(3,i,1) = latu(ip1,1)
             clatv(4,i,1) = latc     (ip1,1)
          ENDDO

!  DO i = 1, ie-2
!    ii = (i*2)+2
!    DO j = 1, je
!      jj = j*2
!      !--     scalar
!      lats(i,j)  = lat(ii,jj)
!      lons(i,j)  = lon(ii,jj)
!      !--     vector - u
!      latu(i,j)  = lat(ii+1,jj)
!      lonu(i,j)  = lon(ii+1,jj)
!    ENDDO
!    !--     vector - v
!    DO j = 1, je-1
!      jj = j*2
!      latv(i,j)  = lat(ii,jj+1)
!      lonv(i,j)  = lon(ii,jj+1)
!    ENDDO
!    latv(i,je) = lat(ii,je*2)
!    lonv(i,je) = lon(ii,je*2)
!    !--     corners of scalar grid cells
!    latc(i,1) = lat(ii+1,1)
!    lonc(i,1) = lon(ii+1,1)
!    DO j = 2, je
!      jj = j*2
!      latc(i,j) = lat(ii+1,jj-1)
!      lonc(i,j) = lon(ii+1,jj-1)
!    ENDDO
!  ENDDO
!  !
!  !--  create corner arrays for SCRIP interpolation
!  !
!  DO i = 1, ie-2
!    ii = (i*2)+2
!    im1 = i-1
!    ip1 = i+1
!    IF (im1 == 0) im1 = ie-2
!    IF (ip1 == ie-1) ip1 = 1
!    !
!    !--     scalar
!    !
!    DO j = 1, je-1
!      clons(1,i,j) = lonc(im1,j  )
!      clons(2,i,j) = lonc(im1,j+1)
!      clons(3,i,j) = lonc(i  ,j+1)
!      clons(4,i,j) = lonc(i  ,j  )
!      clats(1,i,j) = latc(im1,j  )
!      clats(2,i,j) = latc(im1,j+1)
!      clats(3,i,j) = latc(i  ,j+1)
!      clats(4,i,j) = latc(i  ,j  )
!    ENDDO
!    clons(1,i,je) = lonc(im1,je)
!    clons(2,i,je) = lonu(im1,je)
!    clons(3,i,je) = lonu(i  ,je)
!    clons(4,i,je) = lonc(i  ,je)
!    clats(1,i,je) = latc(im1,je)
!    clats(2,i,je) = latu(im1,je)
!    clats(3,i,je) = latu(i  ,je)
!    clats(4,i,je) = latc(i  ,je)
!    !
!    !--     vector - u
!    !
!    DO j = 2, je
!      clonu(1,i,j) = lonv(i  ,j-1)
!      clonu(2,i,j) = lonv(i  ,j  )
!      clonu(3,i,j) = lonv(ip1,j  )
!      clonu(4,i,j) = lonv(ip1,j-1)
!      clatu(1,i,j) = latv(i  ,j-1)
!      clatu(2,i,j) = latv(i  ,j  )
!      clatu(3,i,j) = latv(ip1,j  )
!      clatu(4,i,j) = latv(ip1,j-1)
!    ENDDO
!    clonu(1,i,1) = lon(ii  ,1)
!    clonu(2,i,1) = lonv(i  ,1)
!    clonu(3,i,1) = lonv(ip1,1)
!    clonu(4,i,1) = lon(ii+2,1)
!    clatu(1,i,1) = lat(ii  ,1)
!    clatu(2,i,1) = latv(i  ,1)
!    clatu(3,i,1) = latv(ip1,1)
!    clatu(4,i,1) = lat(ii+2,1)
!    !
!    !--     vector - v
!    !
!    DO j = 1, je-1
!      clonv(1,i,j) = lonu(im1,j  )
!      clonv(2,i,j) = lonu(im1,j+1)
!      clonv(3,i,j) = lonu(i  ,j+1)
!      clonv(4,i,j) = lonu(i  ,j  )
!      clatv(1,i,j) = latu(im1,j  )
!      clatv(2,i,j) = latu(im1,j+1)
!      clatv(3,i,j) = latu(i  ,j+1)
!      clatv(4,i,j) = latu(i  ,j  )
!    ENDDO
!    clonv(1,i,je) = lonu(im1,je)
!    clonv(2,i,je) = lonu(im1,je)
!    clonv(3,i,je) = lonu(i  ,je)
!    clonv(4,i,je) = lonu(i  ,je)
!    clatv(1,i,je) = latu(im1,je)
!    clatv(2,i,je) = latu(im1,je)
!    clatv(3,i,je) = latu(i  ,je)
!    clatv(4,i,je) = latu(i  ,je)
!  ENDDO


  masks(:,:) = 1
  masku(:,:) = 1
  maskv(:,:) = 1
!  WHERE (depto(:,:) .gt. 0)
!     masks(:,:) = 1
!  END WHERE
!  WHERE (deuto(:,:) .gt. 0)
!     maskv(:,:) = 1
!  END WHERE
!  WHERE (deute(:,:) .gt. 0)
!     masku(:,:) = 1
!  END WHERE

  ! write scalar grid to ofilenames.nc

  flen = LEN_TRIM(ofilename)

  ofilename(flen+1:flen+5) = 's.nc'
  CALL write_grid(ofilename, je, ie, lats, lons, masks, clats, clons)

  ! write vector-u grid to ofilenameu.nc

  ofilename(flen+1:flen+5) = 'u.nc'
  CALL write_grid(ofilename, je, ie, latu, lonu, masku, clatu, clonu)

  ! write vector-v grid to ofilenamev.nc

  ofilename(flen+1:flen+5) = 'v.nc'
  CALL write_grid(ofilename, je, ie, latv, lonv, maskv, clatv, clonv)


  DEALLOCATE(lon)
  DEALLOCATE(lat)
  DEALLOCATE(deuto)
  DEALLOCATE(deute)
  DEALLOCATE(lons)
  DEALLOCATE(lats)
  DEALLOCATE(masks)
  DEALLOCATE(lonu)
  DEALLOCATE(latu)
  DEALLOCATE(masku)
  DEALLOCATE(lonv)
  DEALLOCATE(latv)
  DEALLOCATE(maskv)
  DEALLOCATE(lonc)
  DEALLOCATE(latc)
  DEALLOCATE(clons)
  DEALLOCATE(clats)
  DEALLOCATE(clonu)
  DEALLOCATE(clatu)
  DEALLOCATE(clonv)
  DEALLOCATE(clatv)

END SUBROUTINE anta2nc

SUBROUTINE nfce(status)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: status

  !INCLUDE 'netcdf.inc'

  IF ( status .NE. NF90_NOERR ) THEN
    WRITE(0,*) NF90_STRERROR(status)
    STOP 'netCDF error'
  END IF

END SUBROUTINE nfce

SUBROUTINE write_grid(ofilename, nlat, nlon, grid_center_lat, grid_center_lon, grid_mask, &
                      grid_corner_lat, grid_corner_lon)

  IMPLICIT NONE

  CHARACTER*256 :: ofilename

  INTEGER :: nlon, nlat

  REAL :: grid_center_lat(nlon, nlat), grid_center_lon(nlon, nlat),  grid_mask(nlon,nlat)
  REAL :: grid_corner_lat(4,nlon, nlat), grid_corner_lon(4,nlon, nlat)


  INTEGER :: grid_size

  INTEGER :: grid_corners, grid_rank, grid_dims(2)
  INTEGER :: fileid, nc_dims_id_3d(3), nc_dims_id_2d(2)
  INTEGER :: nc_gridsize_id, nc_gridcorn_id, nc_gridrank_id, nc_griddims_id
  INTEGER :: nc_gridxsize_id, nc_gridysize_id
  INTEGER :: nc_grdcntrlat_id, nc_grdcntrlon_id, nc_grdmask_id
  INTEGER :: nc_grdcrnrlat_id, nc_grdcrnrlon_id

  INCLUDE 'netcdf.inc'

  grid_size    = nlon*nlat
  grid_dims(1) = nlon
  grid_dims(2) = nlat
  grid_corners = 4
  grid_rank    = 2

  WRITE(*,*) 'write grid to ', ofilename
  !***
  !*** create netCDF dataset for this grid
  !***
  CALL nfce(nf90_create (ofilename, NF90_CLOBBER, fileid))

  CALL nfce(nf90_put_att (fileid, NF90_GLOBAL, 'title', ofilename))

  !***
  !*** define grid size dimension
  !***
  CALL nfce(nf90_def_dim (fileid, 'grid_size', grid_size, nc_gridsize_id))
  CALL nfce(nf90_def_dim (fileid, 'grid_xsize', nlon, nc_gridxsize_id))
  CALL nfce(nf90_def_dim (fileid, 'grid_ysize', nlat, nc_gridysize_id))

  !***
  !*** define grid corner dimension
  !***
  CALL nfce(nf90_def_dim (fileid, 'grid_corners', grid_corners, nc_gridcorn_id))

  !***
  !*** define grid rank dimension
  !***
  CALL nfce(nf90_def_dim (fileid, 'grid_rank', grid_rank, nc_gridrank_id))

  !***
  !*** define grid dimension size array
  !***

  CALL nfce(nf90_def_var (fileid, 'grid_dims', NF90_INT, nc_gridrank_id, nc_griddims_id))

  !***
  !*** define grid center latitude array
  !***
  nc_dims_id_2d(1) = nc_gridxsize_id
  nc_dims_id_2d(2) = nc_gridysize_id

  CALL nfce(nf90_def_var (fileid, 'grid_center_lat', NF90_REAL, nc_dims_id_2d, nc_grdcntrlat_id))

  CALL nfce(nf90_put_att (fileid, nc_grdcntrlat_id, 'units', 'degrees'))
  CALL nfce(nf90_put_att (fileid, nc_grdcntrlat_id, 'bounds', 'grid_corner_lat'))

  !***
  !*** define grid center longitude array
  !***
  CALL nfce(nf90_def_var (fileid, 'grid_center_lon', NF90_REAL, nc_dims_id_2d, nc_grdcntrlon_id))

  CALL nfce(nf90_put_att (fileid, nc_grdcntrlon_id, 'units', 'degrees'))
  CALL nfce(nf90_put_att (fileid, nc_grdcntrlon_id, 'bounds', 'grid_corner_lon'))

  !***
  !*** define grid mask
  !***
  CALL nfce(nf90_def_var (fileid, 'grid_mask', NF90_REAL, nc_dims_id_2d, nc_grdmask_id))

  CALL nfce(nf90_put_att (fileid, nc_grdmask_id, 'units', 'unitless'))
  CALL nfce(nf90_put_att (fileid, nc_grdmask_id, 'coordinates', 'grid_center_lon grid_center_lat'))

  !***
  !*** define grid corner latitude array
  !***
  nc_dims_id_3d(1) = nc_gridcorn_id
  nc_dims_id_3d(2) = nc_gridxsize_id
  nc_dims_id_3d(3) = nc_gridysize_id

  CALL nfce(nf90_def_var (fileid, 'grid_corner_lat', NF90_REAL, nc_dims_id_3d, nc_grdcrnrlat_id))

  CALL nfce(nf90_put_att (fileid, nc_grdcrnrlat_id, 'units','degrees'))

  !***
  !*** define grid corner longitude array
  !***
  CALL nfce(nf90_def_var (fileid, 'grid_corner_lon', NF90_REAL, nc_dims_id_3d, nc_grdcrnrlon_id))

  CALL nfce(nf90_put_att (fileid, nc_grdcrnrlon_id, 'units', 'degrees'))

  !***
  !*** end definition stage
  !***
  CALL nfce(nf90_enddef(fileid))

  !-----------------------------------------------------------------------
  !
  !     write grid data
  !
  !-----------------------------------------------------------------------

  CALL nfce(nf90_put_var(fileid, nc_griddims_id, grid_dims))

  CALL nfce(nf90_put_var(fileid, nc_grdmask_id, grid_mask))

  CALL nfce(nf90_put_var(fileid, nc_grdcntrlat_id, grid_center_lat))

  CALL nfce(nf90_put_var(fileid, nc_grdcntrlon_id, grid_center_lon))

  CALL nfce(nf90_put_var(fileid, nc_grdcrnrlat_id, grid_corner_lat))

  CALL nfce(nf90_put_var(fileid, nc_grdcrnrlon_id, grid_corner_lon))

  CALL nfce(nf90_close(fileid))

END SUBROUTINE write_grid

  ! ------------------------------------------------------------------------

  SUBROUTINE read_nml(status, iou, fname)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    INTEGER,          INTENT(IN)  :: iou
    CHARACTER(LEN=*), INTENT(IN)  :: fname

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_nml'
    LOGICAL :: lex   ! file exists
    INTEGER :: fstat ! file status
    INTEGER :: jt

    NAMELIST /CTRL/ files_path, grid_name, l_tripolar, IE,JE,rlat1,rlon1,rlat2,rlon2, &
                    phiread, levels, l_init_clim, phc_files_path
    status = 1 ! ERROR

    WRITE(*,*) '==========================================================='

    ! CHECK IF FILE EXISTS
    INQUIRE(file=TRIM(fname), exist=lex)
    IF (.NOT.lex) THEN
       WRITE(*,*) substr,': FILE DOES NOT EXIST (',TRIM(fname),')'
       status = 1
       RETURN
    END IF

    ! OPEN FILE
    OPEN(iou,file=TRIM(fname))

    ! READ NEMELIST
    WRITE(*,*) 'READING NAMELIST ''CTRL'''//&
         &' FROM '''//TRIM(fname),''' (unit ',iou,') ...'
    !
    READ(iou, NML=CTRL, IOSTAT=fstat)
    !
    IF (fstat /= 0) THEN
       WRITE(*,*) substr,': READ ERROR IN NAMELIST ''CTRL'' (',TRIM(fname),')'
       status = 3  ! READ ERROR IN NAMELIST
       RETURN
    END IF

    write (*,*) hline
    write (*,*) hline
    write (*,*) "       GRID NAME           "
    write (*,*) "    ", grid_name
    write (*,*) hline
    write (*,*) hline
    write (*,*) "      INPUT FILES PATH:    "
    write (*,*) "    ", files_path
    write (*,*) hline
    write (*,*) "      PHC FILES PATH:      "
    write (*,*) "    ", phc_files_path
    write (*,*) hline
    write (*,*) "         GRID SPECIFIC:    "
    write (*,*) hline
    write (*,*) "                           "
    write (*,*) "         TRIPOLAR             = ", l_tripolar
    write (*,*) "         HORIZONTAL DIMENSION = ", IE
    write (*,*) "         VERTICAL   DIMENSION = ", JE
    write (*,*) "         FIRST POLE LAT       = ", rlat1
    write (*,*) "         FIRST POLE LON       = ", rlon1
    write (*,*) "         SECOND POLE LAT      = ", rlat2
    write (*,*) "         SECOND POLE LON      = ", rlon2
    write (*,*) "         SIZE OF 1ST POLE     = ", phiread
    write (*,*) "         LEVELS               = ", levels
    write (*,*) hline


    ! CLOSE FILE
    CLOSE(iou)
    status = 0

  END SUBROUTINE read_nml
  ! ------------------------------------------------------------------------

      END PROGRAM MAIN
