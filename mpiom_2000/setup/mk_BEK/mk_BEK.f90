PROGRAM sam
  USE mo_constants, ONLY: aradtogra
  USE mo_kind, ONLY: dp, i8!, i4
  IMPLICIT NONE
  INTEGER :: ie, je
  CHARACTER*6 blord
  REAL(dp), ALLOCATABLE :: lat(:,:),lon(:,:)
  REAL(dp), ALLOCATABLE :: depto(:,:)
  REAL(dp), ALLOCATABLE :: ribek(:,:)
  REAL(dp), ALLOCATABLE :: gi(:,:)
  INTEGER :: i, i1, i2, il, ir, ito, itt, &
       j, jjj, jmm, jmmm, jmanf, jmend, jto
  INTEGER(i8) :: ext8_hdr(4)
  INTEGER, PARAMETER :: io_in_topo=10, io_inout_bek=20, &
       io_in_anta=21, io_out_ribek=68
  REAL(dp), PARAMETER :: maxdeg=360._dp, mindeg=0._dp
  CHARACTER(len=262) :: topo_format
  CHARACTER(len=255) :: frmt
  LOGICAL :: bounds_exch_tp
  INTEGER :: verbose
  NAMELIST /gridparams/ ie, je, frmt, bounds_exch_tp, verbose
  ! default values
  frmt = 'F6.0'
  ie = -1
  je = -1
  bounds_exch_tp = .FALSE.
  verbose = 0
  ! setup run-time parameters
  READ (*, nml=gridparams)
  topo_format = '(i5,20' // TRIM(frmt) // ')'
  IF (verbose .GT. 0) PRINT *, 'ie=', ie, 'je=', je, &
       'format=', TRIM(frmt), 'tp=', bounds_exch_tp
  IF (ie < 1 .OR. je < 1) THEN
    PRINT *, 'Only positive integral values may be used for ie and je!'
    STOP 1
  END IF
  ito = 2 * ie
  jto = 2 * je
  ALLOCATE(lat(ie,je), lon(ie,je), &
       depto(ie,je), ribek(ie,je), gi(ito,jto))

  OPEN(io_in_topo, file='topo', &
       access='sequential', form='formatted', action='read')
  DO i1=2,ie-1,20
    i2 = MIN(i1 + 19, ie - 1)
!    WRITE(6,*) 'lesstreifen ', i1, i2
    READ(io_in_topo,*) blord
    DO j=1,je
      READ (io_in_topo, topo_format) jjj, (depto(i,j), i=i1,i2)
    END DO
  END DO
  CLOSE(io_in_topo)

  OPEN(io_inout_bek, file='BEK', form='formatted', action='readwrite')

  DO j=1,je
    DO i=1,ie
      ribek(i,j)=0._dp
      IF(depto(i,j) .GT. 0.5_dp) ribek(i,j)=0._dp
    END DO
  END DO

  jmmm = (je - 1) / 120

  DO jmm=0,jmmm
    jmanf=1+jmm*120
    jmend=MIN((jmm+1)*120,je)
!    DO i=2,ie-1
!      read(io_inout_bek,'(i4,2x,120i1)',err=99,end=99)
!     & ijju,(ribek(i,j),j=jmend,jmanf,-1)
!    END DO
  END DO

  OPEN(io_in_anta, file='anta.ext', form='unformatted', action='read')

  READ(io_in_anta) ext8_hdr
  READ(io_in_anta) gi

  DO j=1,je
    DO i=1,ie
      lon(i,j) = gi(2*i, 2*j) * aradtogra
      IF (lon(i,j) .LT. mindeg) lon(i,j)=lon(i,j)+maxdeg
      IF (lon(i,j) .GT. maxdeg) lon(i,j)=lon(i,j)-maxdeg
    END DO
  END DO

  READ(io_in_anta) ext8_hdr
  READ(io_in_anta) gi
  CLOSE(io_in_anta)

  DO j=1,je
    DO i=1,ie
      lat(i,j) = gi(2*i, 2*j) * aradtogra
    ENDDO
  ENDDO


  DO j=1,je
    ribek(1,j)=ribek(ie-1,j)
    ribek(ie,j)=ribek(2,j)
  END DO

  DO j=2,je-1
    DO i=2,ie-1

      IF (lat(i,j) .LE. -49._dp)ribek(i,j)=5._dp
      IF (lat(i,j) .GT. 63._dp) ribek(i,j)=1._dp
      IF (lat(i,j) .GT. 80._dp) ribek(i,j)=2._dp

      IF(lon(i,j) .GE. 20._dp .AND. lon(i,j) .LE. 290._dp) THEN
        IF(lat(i,j) .LE.  0._dp) ribek(i,j) = 7._dp
        IF(lat(i,j) .GT. 62._dp) ribek(i,j) = 7._dp
        IF(lat(i,j) .GT. 63._dp) ribek(i,j) = 2._dp
      ELSE
        IF(lat(i,j) .LE. -30._dp) ribek(i,j)=5._dp
        IF(ribek(i,j) .NE. 5._dp .AND. lat(i,j) .LE. 30._dp) ribek(i,j)=5._dp
        IF(ribek(i,j) .NE. 5._dp .AND. lat(i,j) .GE. 30._dp) ribek(i,j)=4._dp
        IF(lat(i,j) .GT. 63._dp) ribek(i,j)=1._dp
        IF(lat(i,j) .GT. 80._dp) ribek(i,j)=2._dp
      END IF

!        if(lon(i,j) .le. 30 .or. lon(i,j) .ge. 270.) then
!        if(ribek(i,j) .eq. 7 .and. lat(i,j) .ge. 0.) ribek(i,j)=4.
!        endif
      IF(lon(i,j) .LE. 90._dp .AND. lat(i,j) .GE. 48._dp) THEN
        IF (ribek(i,j) .EQ. 4._dp) ribek(i,j)=1._dp
      ENDIF
      IF(lon(i,j) .GE. 260._dp .AND. lon(i,j) .LE. 316._dp) THEN
        IF (lat(i,j) .GE. 51._dp .AND. lat(i,j) .LE. 80._dp) ribek(i,j)=3._dp
      ENDIF

      IF(lon(i,j) .GE. 0._dp .AND. lon(i,j) .LE. 30._dp) THEN
        IF (lat(i,j) .GE. 47._dp .AND. lat(i,j) .LE. 80._dp) ribek(i,j)=1._dp
      ENDIF

      IF(lon(i,j) .GE. 354._dp .OR. lon(i,j) .LE. 50._dp) THEN
        IF (lat(i,j) .GE. 27._dp .AND. lat(i,j) .LE. 41._dp) THEN
          IF (ribek(i,j) .EQ. 4._dp) ribek(i,j)=9._dp
        ENDIF
      ENDIF
      IF(lon(i,j) .GE. 0._dp .AND. lon(i,j) .LE. 50._dp) THEN
        IF (lat(i,j) .GE. 27._dp .AND. lat(i,j) .LE. 47._dp) THEN
          IF (ribek(i,j) .EQ. 4._dp) ribek(i,j)=9._dp
        ENDIF
      ENDIF

      IF(lat(i,j) .LE. -50._dp) ribek(i,j)=6._dp


    END DO
  END DO


  DO itt=1,5000
    IF (verbose .GT. 1) PRINT *, itt
    DO j=2,je-1
      DO i=2,ie-1
        IF(depto(i,j) .GT. 0.5_dp .AND. ribek(i,j) .LT. 0.5_dp)THEN
          ribek(i,j)=MAX(ribek(i+1,j),ribek(i-1,j) &
               ,ribek(i,j+1),ribek(i,j-1))
        END IF
      END DO
    END DO
  END DO


  DO j=2,je-1
    DO i=2,ie-1
      IF (lat(i,j) .LE. 5._dp .AND. lat(i,j) .GE. -5._dp) THEN
        IF (lon(i,j) .GE. 210._dp .AND. lon(i,j) .LE. 270._dp) THEN
          ribek(i,j)=8._dp
        END IF
      END IF
    END DO
  END DO



  DO j=1,je
    DO i=1,ie
      IF (depto(i,j) .LT. 0.5_dp) ribek(i,j) =  0._dp
      IF (depto(i,j) .LT. 0.5_dp) lat(i,j)   = -9e33_dp
      IF (depto(i,j) .LT. 0.5_dp) lon(i,j)   = -9e33_dp
    END DO
  END DO

  ribek(1,:)=ribek(ie-1,:)
  ribek(ie,:)=ribek(2,:)

  IF (bounds_exch_tp) THEN
    DO i=1,ie
      il=i
      ir=ie+1-i
      ribek(il,2) = ribek(ir,3)  ! syncronise line 2 with line 3
      ribek(il,1) = ribek(ir,4)  ! syncronise line 1 with line 4
    END DO
  END IF

  REWIND io_inout_bek
  DO jmm = 0, jmmm
    jmanf = 1+jmm*120
    jmend = MIN((jmm+1)*120,je)
!    WRITE(6,*)'ribek, jm ',jmm,jmend,jmanf
    DO i = 2, ie-1
      WRITE(io_inout_bek,'(i4,2x,120i1)') i,(INT(ribek(i,j)),j=jmend,jmanf,-1)
    ENDDO
  ENDDO
  CLOSE(io_inout_bek)

  ext8_hdr(1) = 1_i8
  ext8_hdr(2) = 1_i8
  ext8_hdr(3) = 0_i8
  ext8_hdr(4) = INT(ie * je, i8)
  OPEN(io_out_ribek, file='RIBEK.ext', form='unformatted')
  WRITE(io_out_ribek) ext8_hdr
  WRITE(io_out_ribek) ribek
  ext8_hdr(2) = ext8_hdr(2) + 1_i8
  WRITE(io_out_ribek) ext8_hdr
  WRITE(io_out_ribek) lat
  ext8_hdr(2) = ext8_hdr(2) + 1_i8
  WRITE(io_out_ribek) ext8_hdr
  WRITE(io_out_ribek) lon
  CLOSE(io_out_ribek)
END PROGRAM sam
