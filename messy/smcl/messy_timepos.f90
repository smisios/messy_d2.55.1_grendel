! **********************************************************************
MODULE messy_timepos
! **********************************************************************

  ! Time & Position sampling
  !
  ! CORE MODULE (MESSy/SMCL)
  !
  ! Author: Patrick Joeckel, MPICH, March 2006
  !

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: ASSOCIATED, TRIM

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'timepos'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '2.0.1'

  TYPE T_TIMEPOS
     INTEGER           :: yr, mo, dy, ho, mi, se
     REAL(DP)          :: lon, lat, pre
     CHARACTER(LEN=15) :: id 
  END TYPE T_TIMEPOS

  TYPE T_IPOS
     INTEGER,  DIMENSION(4,2)  :: ixy    ! global longitude/latiude indices
     REAL(DP), DIMENSION(4)    :: hw     ! horizontal weight
     LOGICAL                   :: lh = .FALSE. ! horizontal weight present
     INTEGER,  DIMENSION(4,2)  :: ipress ! indices of neighbouring press. lev.
     REAL(DP), DIMENSION(4,2)  :: vw     ! vertical weight
     LOGICAL                   :: lv = .FALSE. ! vertical weight present
  END TYPE T_IPOS

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: fstr = &
       '(i4,i2.2,i2.2,1x,i2.2,i2.2,1x,f9.4,1x,f8.4,1x,i6,1x,a14)'

  PUBLIC :: T_TIMEPOS
  PUBLIC :: T_IPOS
  PUBLIC :: read_timepos
  PUBLIC :: set_hweight
  PUBLIC :: set_vweight
  PUBLIC :: debug_output
  PUBLIC :: interpolate

CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE read_timepos(status, fname, iou, list)

    IMPLICIT NONE

    INTRINSIC :: REAL

    ! I/O
    INTEGER,                       INTENT(OUT) :: status 
    CHARACTER(LEN=*),              INTENT(IN)  :: fname  ! filename
    INTEGER,                       INTENT(IN)  :: iou    ! I/O unit
    TYPE(T_TIMEPOS), DIMENSION(:), POINTER     :: list   ! list of t-pos

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_timepos'
    LOGICAL :: lex   ! file exists ?
    LOGICAL :: lopn  ! file open ?
    INTEGER :: iout  ! unit ?
    INTEGER :: n     ! number of data lines
    INTEGER :: mstat ! memory status
    INTEGER :: fstat ! file status
    INTEGER :: i     ! counter
    INTEGER :: ip    ! integer pressure

    status = 1 ! ERROR

    ! CHECK, IF FILE IS PRESENT
    INQUIRE(file=TRIM(fname), exist=lex, opened=lopn, number=iout)
    IF (.NOT. lex) THEN
       status = -1
       WRITE(*,*) substr,': *** WARNING: File '''//TRIM(fname)//&
            &''' not found !'
       RETURN
    END IF
    IF (lopn) THEN
       WRITE(*,*) substr,': *** ERROR: File  '''//TRIM(fname)//&
            &''' already open on unit ',iout
       RETURN ! ERROR
    END IF

    ! OPEN FILE
    OPEN(unit=iou,file=TRIM(fname))
    ! COUNT LINES
    n = 0
    DO
       READ(iou, *, IOSTAT=fstat)
       IF (fstat < 0) EXIT
       n = n + 1
    END DO
    CLOSE(iou)
    n = n - 1  ! one header line

    ! ALLOCATE DATA
    IF (ASSOCIATED(list)) THEN
       DEALLOCATE(list)
       NULLIFY(list)
    END IF
    ALLOCATE(list(n), STAT=mstat)
    IF (mstat /= 0) THEN
       WRITE(*,*) '*** ERROR: Memory allocation failed !'
       RETURN ! ERROR
    END IF

    ! READ DATA
    OPEN(unit=iou,file=TRIM(fname))

    READ(iou, *) ! one header line
    DO i=1, n

       READ(iou, fstr, IOSTAT=fstat) &
            list(i)%yr,  list(i)%mo,  list(i)%dy,   &
            list(i)%ho,  list(i)%mi,                &
            list(i)%lon, list(i)%lat, ip,           &
            list(i)%id

       list(i)%pre = REAL(ip, dp)
       list(i)%se  = 0

       IF (fstat /= 0) THEN
          WRITE(*,*) substr,': *** READ ERROR in line ',i,' of file '''&
               &//TRIM(fname)//''''
          RETURN ! ERROR
       END IF

    END DO

    ! END
    CLOSE(iou)
    status = 0

  END SUBROUTINE read_timepos
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE set_hweight(status, lon, lat, tpos, ipos)

    ! MESSy
    USE messy_main_tools,       ONLY: nn_index, bilin_weight

    IMPLICIT NONE

    INTRINSIC :: MAX, MIN

    ! I/O
    INTEGER,                    INTENT(OUT)   :: status
    REAL(DP), DIMENSION(:),     INTENT(IN)    :: lon
    REAL(DP), DIMENSION(:),     INTENT(IN)    :: lat
    TYPE(t_timepos),            INTENT(IN)    :: tpos
    TYPE(t_ipos),               INTENT(INOUT) :: ipos

    ! LOCAL
    REAL(DP)                 :: zlon
    INTEGER                  :: jgx1, jgx2, jgy1, jgy2
    INTEGER                  :: jgx, jgy
    REAL(DP), DIMENSION(4,2) :: vn
    REAL(DP), DIMENSION(2)   :: v

    ! CHECKS
    IF ((tpos%lon < -180.0_DP) .OR. (tpos%lon > 360.0_DP)) THEN
       WRITE(*,*) 'ERROR: LONGITUDE OUT OF RANGE'
       status = 1  ! LONGITUDE OUT OF RANGE
       RETURN
    ELSE
       ! CHANGE INTERVAL [-180 ... 180] -> [0 ... 360]
       IF (tpos%lon <= 0.0_DP) THEN
          zlon = tpos%lon + 360.0_DP
       ELSE
          zlon = tpos%lon
       END IF
    END IF
    !
    IF ((tpos%lat < -90.0_DP) .OR. (tpos%lat > 90.0_DP)) THEN
       WRITE(*,*) 'ERROR: LATITUDE OUT OF RANGE'
       status = 2  ! LATITUDE OUT OF RANGE
       RETURN
    END IF

    ! SET INDICES
    CALL nn_index(lon,      zlon, jgx1, jgx2)
    CALL nn_index(lat,  tpos%lat, jgy1, jgy2)

    v(1) = zlon
    v(2) = tpos%lat

    jgx = MIN(jgx1, jgx2)
    jgy = MIN(jgy1, jgy2)
    vn(1,1) = lon(jgx)
    vn(1,2) = lat(jgy)
    ipos%ixy(1,1) = jgx
    ipos%ixy(1,2) = jgy

    jgx = MAX(jgx1, jgx2)
    jgy = MIN(jgy1, jgy2)
    vn(2,1) = lon(jgx)
    vn(2,2) = lat(jgy)
    ipos%ixy(2,1) = jgx
    ipos%ixy(2,2) = jgy

    jgx = MAX(jgx1, jgx2)
    jgy = MAX(jgy1, jgy2)
    vn(3,1) = lon(jgx)
    vn(3,2) = lat(jgy)
    ipos%ixy(3,1) = jgx
    ipos%ixy(3,2) = jgy

    jgx = MIN(jgx1, jgx2)
    jgy = MAX(jgy1, jgy2)
    vn(4,1) = lon(jgx)
    vn(4,2) = lat(jgy)
    ipos%ixy(4,1) = jgx
    ipos%ixy(4,2) = jgy

    CALL bilin_weight(vn, v, ipos%hw)
    
    ipos%lh = .TRUE.

    status = 0    

  END SUBROUTINE set_hweight
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE set_vweight(status, press, tpos, ipos)

    ! MESSy
    USE messy_main_tools,       ONLY: nn_index

    IMPLICIT NONE

    INTRINSIC :: ABS, SIZE

    ! I/O
    INTEGER,                    INTENT(OUT)   :: status
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)    :: press
    TYPE(t_timepos),            INTENT(IN)    :: tpos
    TYPE(t_ipos),               INTENT(INOUT) :: ipos

    ! LOCAL
    INTEGER  :: i
    INTEGER  :: jgx, jgy
    INTEGER  :: k1, k2
    INTEGER  :: nk
    REAL(DP) :: zpre
    REAL(DP) :: wp

    IF (.NOT. ipos%lh) THEN
       WRITE(*,*) 'ERROR: VERTICAL INTERPOALTION'//&
            &' REQUIRES HORIZONTAL INTERPOLATION'
       status = 1
       RETURN
    END IF

    nk = SIZE(press,2)

    DO i=1, 4
       zpre = tpos%pre ! Pa
       jgx = ipos%ixy(i,1)
       jgy = ipos%ixy(i,2)
       !zpre = MAX(zpre, press(jgx,  1, jgy))  ! not above upper layer
       !zpre = MIN(zpre, press(jgx, nk, jgy))  ! not below lowest layer
       IF (zpre < press(jgx,  1, jgy)) THEN
          zpre = press(jgx,  1, jgy)
          k1 = 1
          k2 = 1
          wp = 0.5_dp
       ELSE IF (zpre > press(jgx,  nk, jgy)) THEN
          zpre = press(jgx,  nk, jgy)
          k1 = nk
          k2 = nk
          wp = 0.5_dp
       ELSE
          CALL nn_index(press(jgx,:,jgy), zpre, k1, k2)
          wp = ABS( (zpre - press(jgx,k1,jgy)) / &
               ( press(jgx, k2, jgy) - press(jgx, k1, jgy)) )
       END IF

       ipos%ipress(i,1) = k1
       ipos%ipress(i,2) = k2
       ipos%vw(i,1) = (1._dp - wp)
       ipos%vw(i,2) = wp
    END DO

    ipos%lv = .TRUE.
    status = 0

  END SUBROUTINE set_vweight
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE debug_output(iou, tpos, ipos)

    IMPLICIT NONE

    INTRINSIC :: INT, SUM 

    ! I/O
    INTEGER,          INTENT(IN) :: iou
    TYPE(t_timepos),  INTENT(IN) :: tpos
    TYPE(t_ipos),     INTENT(IN) :: ipos

    ! LOCAL
    INTEGER :: ip
    INTEGER :: j

    ip = INT(tpos%pre)
    WRITE(iou,fstr) &
         tpos%yr,  tpos%mo,  tpos%dy,   &
         tpos%ho,  tpos%mi,             &
         tpos%lon, tpos%lat, ip,        &
         tpos%id

    IF (ipos%lh) THEN
       DO j=1, 4
          WRITE(iou,'(a11,i2,a3,i4,a1,i4,a5,e12.6)') &
               '  (LON,LAT)',j,': (' &
               , ipos%ixy(j,1),',',ipos%ixy(j,2),'); w=',ipos%hw(j)
       END DO
       WRITE(iou,'(14x,a16,f12.6)') 'SUM OF WEIGHTS: ',SUM(ipos%hw(:))
    ELSE
       WRITE(iou,*) '  HORIZONTAL INTERPOLATION UNDEFINED'
    END IF
    
    IF (ipos%lv) THEN
       DO j=1, 4
          WRITE(iou,'(a11,i2,a3,i4,a1,i4,a5,e12.6,1x,e12.6,1x,f12.6)') &
               '  (k1 ,k2 )',j,': (' &
               , ipos%ipress(j,1),',',ipos%ipress(j,2),'); w=',ipos%vw(j,:) &
               , SUM(ipos%vw(j,:))
       END DO
    ELSE
       WRITE(iou,*) '  VERTICAL INTERPOLATION UNDEFINED'
    END IF

  END SUBROUTINE debug_output
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE interpolate(status, zfield, ipos, value)

    IMPLICIT NONE

    ! I/O
    INTEGER,                    INTENT(OUT) :: status
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: zfield
    TYPE(t_ipos),               INTENT(IN)  :: ipos
    REAL(DP),                   INTENT(OUT) :: value

    ! LOCAL
    INTEGER :: i
    INTEGER :: jgx, jgy, k1, k2

    IF ((.NOT. ipos%lh) .OR. (.NOT. ipos%lv)) THEN
       WRITE(*,*) 'ERROR: WEIGHTS NOT SET FOR INTERPOLATION'
       status = 1
       RETURN
    END IF

    value = 0.0
    DO i=1,4
       jgx = ipos%ixy(i,1)
       jgy = ipos%ixy(i,2)
       k1  = ipos%ipress(i,1)
       k2  = ipos%ipress(i,2)
       value = value &
            + ( zfield(jgx,k1,jgy) * ipos%vw(i,1) &
            +   zfield(jgx,k2,jgy) * ipos%vw(i,2) &
              ) * ipos%hw(i)
    END DO

    status = 0

  END SUBROUTINE interpolate
! ----------------------------------------------------------------------

! **********************************************************************
END MODULE messy_timepos
! **********************************************************************
