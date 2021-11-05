! ============================================================================
MODULE messy_main_import_ts
! ============================================================================

  USE messy_main_constants_mem, ONLY: DP, STRLEN_MEDIUM, STRLEN_ULONG &
                                    , STRLEN_VLONG
  USE messy_main_import,        ONLY: modstr
  USE messy_main_grid_netcdf,   ONLY: t_ncatt

  IMPLICIT NONE
  PRIVATE
  SAVE
  INTRINSIC :: NULL, HUGE

  PUBLIC :: DP

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodstr = 'import_ts'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodver = '0.1'

  INTEGER, PARAMETER :: STRLEN_OBJECT  = 2*STRLEN_MEDIUM + 1 + 4

  INTEGER, PARAMETER, PUBLIC :: TS_BD_STOP = 0  ! ERROR, IF OUT OF BOUNDS
  INTEGER, PARAMETER, PUBLIC :: TS_BD_CONT = 1  ! CONTINUE WITH LAST VALUE

  INTEGER, PARAMETER, PUBLIC :: TS_IM_PREV = -1  ! USE PREVIOUS VALUE
  INTEGER, PARAMETER, PUBLIC :: TS_IM_NEXT = 1   ! USE NEXT VALUE
  INTEGER, PARAMETER, PUBLIC :: TS_IM_LINT = 0   ! INTERPOLATE LINEARLY

  TYPE T_TS_IO
     !
     ! NAME OF TIME SERIES
     CHARACTER(LEN=STRLEN_OBJECT) :: name = ''
     ! NAME OF FILE WITH DATA
     CHARACTER(LEN=STRLEN_ULONG)  :: fname = ''
     ! VALID RANGE
     REAL(DP), DIMENSION(2)       :: vr = (/ -HUGE(0.0_dp), HUGE(0.0_dp) /)
     ! WHAT TODO IF BEYOND LOWER/UPPER BOUNDARY
     INTEGER, DIMENSION(2)        :: cnt = (/TS_BD_STOP, TS_BD_STOP/)
     ! INTERPOLATION METHOD
     INTEGER                      :: im = TS_IM_PREV
     !
     ! PICK OUT THIS DATE/TIME (year, month, day, hour, minute, second)
     INTEGER, DIMENSION(6) :: pdt = (/ -1, -1, -1, -1, -1, -1 /)
     !
     ! SHIFT BY THIS NUMBER OF DAYS
     REAL(DP) :: offset = 0.0_dp
     !     
  END TYPE T_TS_IO
  PUBLIC :: T_TS_IO

  TYPE T_TS
     !
     ! IO
     TYPE(T_TS_IO) :: io
     !
     ! NUMBER OF TIME STEPS IN SERIES
     INTEGER :: nt = 0
     !
     ! NUMBER OF PARAMETERS IN SERIES
     INTEGER :: np = 0
     !
     ! TIME AXIS
     REAL(DP), DIMENSION(:), POINTER :: jd => NULL()  ! Julian day + fract.
     !
     ! 'PARAMETER' AXIS
     REAL(DP), DIMENSION(:), POINTER :: par => NULL()
     !
     ! DATA (RANK-1: time, RANK-2: number of parameters)
     REAL(DP), DIMENSION(:,:), POINTER :: data => NULL()
     !
     ! 'CURRENT' VALUE (channel object)
     LOGICAL                         :: lalloc = .FALSE.
     REAL(DP), DIMENSION(:), POINTER :: obj => NULL()
     !
     ! 'FLAG' VALUE (1: OK, 0: OUT OF VALID RANGE)
     REAL(DP), DIMENSION(:), POINTER :: flg => NULL()
     !
     ! NAME OF PARAMETER DIMENSION
     CHARACTER(LEN=STRLEN_MEDIUM) :: dimname  = ''
     !
     ! ATTRIBUTES OF VARIABLE
     TYPE (t_ncatt), DIMENSION(:), POINTER :: varatt => NULL()
     !
     ! ATTRIBUTES OF PARAMETER DIMENSION
     TYPE (t_ncatt), DIMENSION(:), POINTER :: dimvaratt => NULL() 
     !
  END TYPE T_TS
  PUBLIC :: T_TS

  ! WORKSPACE
  INTEGER,       PARAMETER,         PUBLIC :: NMAXTS = 100 ! max. numer
  INTEGER,                          PUBLIC :: NTS = 0      ! actual number
  TYPE(T_TS_IO), DIMENSION(NMAXTS), PUBLIC :: TS
  TYPE(T_TS),    DIMENSION(NMAXTS), PUBLIC :: XTS

  PUBLIC :: its_read_ts
  !PRIVATE :: its_read_ts_ascii
  !PRIVATE :: its_read_ts_netcdf
  PUBLIC :: its_delete_ts
  PUBLIC :: its_set_value_ts
  PUBLIC :: its_copy_io
  PUBLIC :: import_ts_read_nml_ctrl

CONTAINS

  ! --------------------------------------------------------------------------
  SUBROUTINE its_read_ts(status, zts, lalloc)

    USE messy_main_tools,      ONLY: strcrack

    IMPLICIT NONE
    INTRINSIC :: LEN_TRIM, ADJUSTL, TRIM, MINVAL, MAXVAL

    ! I/O
    INTEGER,    INTENT(OUT)   :: status
    TYPE(T_TS), INTENT(INOUT) :: zts
    LOGICAL,    INTENT(IN)    :: lalloc
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'its_read_ts'
    CHARACTER(LEN=STRLEN_ULONG) :: str = ''
    CHARACTER(LEN=STRLEN_ULONG) :: fname = ''
    CHARACTER(LEN=STRLEN_ULONG) :: vname = ''
    LOGICAL                     :: lex = .FALSE.
    CHARACTER(LEN=3)            :: ext = '   '
    INTEGER                     :: n = 0
    ! field of substrings
    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: el => NULL()

    status = 0
    IF (zts%io%name == '') RETURN

    str = ADJUSTL(TRIM(zts%io%fname))

    n = LEN_TRIM(str)
    ext(1:3) = str(n-2:n)

    IF (ext == '.nc') THEN

       CALL strcrack(TRIM(ADJUSTL(zts%io%fname)), '@', el, n)
       IF (n /= 2) THEN
          ! ERROR
          WRITE(*,*) substr &
               , ' ERROR: SYNTAX ERROR IN VARIABLE SELECTION OF TIME SERIES ' &
               , TRIM(zts%io%name)
          RETURN
       END IF
       fname = TRIM(ADJUSTL(el(2)))
       vname = TRIM(ADJUSTL(el(1)))
       DEALLOCATE(el); NULLIFY(el)

    ELSE

       fname = str
       vname = ''

    END IF

    INQUIRE(file=TRIM(fname), exist=lex)
    IF (.NOT. lex) THEN
       WRITE(*,*) substr,' ERROR: file ',TRIM(str),' does not exist!'
       status = 1
       RETURN
    END IF

    IF (ext == '.nc') THEN
       CALL its_read_ts_netcdf(status, zts, fname, vname)
    ELSE
       CALL its_read_ts_ascii(status, zts)
    END IF

    WRITE(*,*) '    ',substr,': RANGE OF JULIAN DAY  : '&
         ,MINVAL(zts%jd),' - ',MAXVAL(zts%jd)
    WRITE(*,*) '    ',substr,': RANGE OF DATA        : '&
         ,MINVAL(zts%data),' - ',MAXVAL(zts%data)
    WRITE(*,*) '    ',substr,': PICK OUT             : ',zts%io%pdt
    WRITE(*,*) '    ',substr,': OFFSET               : ',zts%io%offset,' days'

    ! MEMORY 3
    zts%lalloc = lalloc
    IF (lalloc) THEN
       ALLOCATE(zts%obj(zts%np))
       zts%obj(:) = 0.0_dp
       ALLOCATE(zts%flg(zts%np))
       zts%flg(:) = 0.0_dp
    ENDIF

  END SUBROUTINE its_read_ts
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE its_read_ts_ascii(status, zts)

    USE messy_main_tools, ONLY: find_next_free_unit
    USE messy_main_timer, ONLY: gregor2julian

    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, TRIM

    ! I/O
    INTEGER,    INTENT(OUT)   :: status
    TYPE(T_TS), INTENT(INOUT) :: zts
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'its_read_ts_ascii'
    CHARACTER(LEN=6), DIMENSION(6), PARAMETER :: tstr = (/ &
         'year  ','month ','day   ','hour  ', 'minute', 'second' /)
    INTEGER  :: iou
    INTEGER  :: fstat
    INTEGER  :: flag
    INTEGER  :: yr,mo,dy,hr,mi,se
    INTEGER  :: i
    INTEGER  :: y1, y2

    ! INIT
    status = 1

    WRITE(*,*) '    ',substr,': OPENING FILE ',TRIM(ADJUSTL(zts%io%fname))

    IF (ASSOCIATED(zts%jd)) DEALLOCATE(zts%jd)
    NULLIFY(zts%jd)

    IF (ASSOCIATED(zts%par)) DEALLOCATE(zts%par)
    NULLIFY(zts%par)

    IF (ASSOCIATED(zts%data)) DEALLOCATE(zts%data)
    NULLIFY(zts%data)

    ! OPEN DATASET
    iou = find_next_free_unit(100,200)
    OPEN(unit=iou,file=TRIM(ADJUSTL(zts%io%fname)))

    ! READ DATA
    READ(iou,*) ! header (line 1)
    READ(iou,*) ! header (line 2)
    READ(iou,*) ! header (line 3)
    READ(iou,*) ! header (line 4): note on syntax
    READ(iou,*) flag, y1, y2, zts%np   ! line 5

    IF ((flag < 1) .OR. (flag > 6)) THEN
       ! ERROR
       WRITE(*,*) substr,' ERROR: unknown time resolution flag (<1 or >6)!'  
       RETURN
    ENDIF
    WRITE(*,*) '    ',substr,': TIME RESOLUTION      : ',tstr(flag)
    WRITE(*,*) '    ',substr,': (y1, y2)             : (',y1,y2,')'
    WRITE(*,*) '    ',substr,': NUMBER OF PARAMETERES: ',zts%np

    ! MEMORY 1
    ALLOCATE(zts%par(zts%np))

    READ(iou,*) ! header (line 6)
    READ(iou,*) zts%par(:) ! line 7
    WRITE(*,*) '    ',substr,': AXIS                 : ',zts%par(:)

    READ(iou,*) ! header (line 8)

    ! COUNT NUMBER OF DATA LINES
    zts%nt = 0
    DO
       READ(iou, *, IOSTAT=fstat)
       IF (fstat < 0) EXIT
       zts%nt = zts%nt + 1
    END DO
    WRITE(*,*) '    ',substr,': NUMBER OF TIME STEPS : ',zts%nt
    CLOSE(iou)

    ! MEMORY 2
    ALLOCATE(zts%jd(zts%nt))
    ALLOCATE(zts%data(zts%nt,zts%np))

    ! PRE-FILL
    zts%jd(:)     = 0.0_dp
    zts%data(:,:) = 0.0_dp

    ! READ DATA
    ! ... OPEN FILE
    OPEN(unit=iou,file=TRIM(ADJUSTL(zts%io%fname)))
    ! ... READ HEADER LINES
    DO i=1,8
       READ(iou,*)
    END DO
    ! ... READ DATA LINES
    read_data_lines: DO i=1, zts%nt
       yr = 0
       mo = 1
       dy = 1
       hr = 0
       mi = 0
       se = 0
       SELECT CASE(flag)
       CASE(1)
          READ(iou, *, IOSTAT=fstat) yr, zts%data(i,:)
       CASE(2)
          READ(iou, *, IOSTAT=fstat) yr, mo, zts%data(i,:)
       CASE(3)
          READ(iou, *, IOSTAT=fstat) yr, mo, dy, zts%data(i,:)
       CASE(4)
          READ(iou, *, IOSTAT=fstat) yr, mo, dy, hr, zts%data(i,:)
       CASE(5)
          READ(iou, *, IOSTAT=fstat) yr, mo, dy, hr, mi, zts%data(i,:)
       CASE(6)
          READ(iou, *, IOSTAT=fstat) yr, mo, dy, hr, mi, se, zts%data(i,:)
       END SELECT

       IF (fstat < 0) THEN
          WRITE(*,*) '    ',substr,': END OF FILE REACHED; ',&
               i,' DATA LINES READ'
          EXIT
       ELSEIF (fstat > 0) THEN
          ! add number of header lines
          WRITE(*,*) substr,' ERROR: READ ERROR IN LINE ',i+8
          RETURN ! ERROR
       END IF

       zts%jd(i) = gregor2julian(yr,mo,dy,hr,mi,se)

       ! CHECK MONOTONICITY
       IF (i > 1) THEN
          IF (zts%jd(i)-zts%jd(i-1) <= 0.0_dp) THEN
             ! ERROR
             WRITE(*,*) substr,' ERROR: t(i) <= t(i-1) at LINE ',i+8
             RETURN
          ENDIF
       END IF

    END DO read_data_lines

    ! CLOSE DATASET
    WRITE(*,*) '    ',substr,': CLOSING FILE ',TRIM(ADJUSTL(zts%io%fname))
    CLOSE(iou)

    ! FINAL SETTINGS
    zts%dimname = 'dim_'//TRIM(ADJUSTL(zts%io%name))

    status = 0

  END SUBROUTINE its_read_ts_ascii
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE its_read_ts_netcdf(status, zts, fname, vname)

    USE messy_main_grid_netcdf, ONLY: t_ncvar, import_ncvar, NULL_DIMID       &
                                    , INIT_NCVAR, NULL_VARID, t_ncatt, string &
                                    , DOUBLE_NARRAY,IMPORT_NCATT, COPY_NCATT  &
                                    , element
    USE messy_main_timer,      ONLY: eval_time_str, gregor2julian

    IMPLICIT NONE
    INTRINSIC :: TRIM, ADJUSTL, PRODUCT, NULL

    ! I/O
    INTEGER,    INTENT(OUT)   :: status
    TYPE(T_TS), INTENT(INOUT) :: zts
    CHARACTER(LEN=*), INTENT(IN) :: fname
    CHARACTER(LEN=*), INTENT(IN) :: vname

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'its_read_ts_netcdf'
    !
    TYPE (t_ncvar) :: var, tvar, pvar
    TYPE (t_ncatt) :: att
    INTEGER      :: i, np, nt
    CHARACTER(LEN=STRLEN_VLONG) :: tunit = ''
    INTEGER      :: z_year, z_month, z_day, z_hour, z_min, z_sec
    REAL(DP)     :: z_tuf
    REAL(DP)     :: jd0
    INTEGER, DIMENSION(:), POINTER :: vec => NULL()
    INTEGER      :: idim_time, idim_par

    status = 1 ! ERROR

    WRITE(*,*) '    ',substr,': OPENING FILE ',TRIM(ADJUSTL(fname))
    WRITE(*,*) '    ',substr,': READING VARIABLE ',TRIM(ADJUSTL(vname))

    ! STEP 2: IMPORT VARIABLE ------------------------------------------------
    CALL IMPORT_NCVAR( var, varname=TRIM(ADJUSTL(vname)) &
         , file=TRIM(ADJUSTL(fname)) )

    ! STEP 3: CONSTRUCT TIME AXIS --------------------------------------------
    ! ... check for unlimited ID
    IF (var%uid == NULL_DIMID) THEN
       WRITE(*,*) substr &
            , ' ERROR: NO UNLIMITED DIMENSION FOUND FOR ', TRIM(ADJUSTL(vname))
       CALL INIT_NCVAR(var)
       RETURN
    END IF
    !
    IF (var%ndims > 2) THEN
       WRITE(*,*) substr &
            , ' ERROR: ',TRIM(ADJUSTL(vname)),' HAS MORE THAN 2 DIMENSIONS '&
            , '(INCLUDING TIME AXIS)'
       CALL INIT_NCVAR(var)
       RETURN       
    END IF
    !
    dimension_loop: DO i=1, var%ndims
       unlimited_id: IF (var%dim(i)%fuid) THEN

          idim_time = i
          zts%nt = var%dim(i)%len
          WRITE(*,*) '    ',substr,': NUMBER OF TIME STEPS : ',zts%nt
          ! GET TIME AXIS INFORMATION
          IF (var%dim(i)%varid == NULL_VARID) THEN
             WRITE(*,*) substr &
                  , ' ERROR: NO COORDINATE VARIABLE FOUND FOR DIMENSION ' &
                  , TRIM(var%dim(i)%name)
             CALL INIT_NCVAR(var)
             RETURN
          ENDIF
          ! attribute "units"
          CALL IMPORT_NCATT(att, varid=var%dim(i)%varid, attname='units' &
               , file=TRIM(ADJUSTL(fname)))
          tunit=TRIM(string(att%dat%vc))
          WRITE(*,*) '    ',substr,': TIME UNIT            : ',TRIM(tunit)
          ! convert: string to date/time
!!$          CALL eval_time_str(status, TRIM(ADJUSTL(tunit)) &
          CALL eval_time_str(status, tunit &
               , z_tuf, z_year, z_month, z_day, z_hour, z_min, z_sec)
          IF (status /= 0) THEN
             WRITE(*,*) substr &
                  , ' ERROR: TIME UNIT NOT IN VALID FORMAT FOR VARIABLE ' &
                  , TRIM(var%dim(i)%name)
             CALL INIT_NCVAR(var)
             RETURN
          END IF
          ! convert conversion factor from "to seconds" to "to days" ...
          z_tuf = z_tuf/86400.0_dp
          ! convert: date/time to JulianDay
          jd0 = gregor2julian(z_year, z_month, z_day, z_hour, z_min, z_sec)
          ! ALLOCATE MEMORY FOR TIME AXIS
          ALLOCATE(zts%jd(zts%nt))
          ! READ DATA
          CALL IMPORT_NCVAR( tvar, varid=var%dim(i)%varid &
               , file=TRIM(ADJUSTL(fname)) )
          CALL DOUBLE_NARRAY(tvar%dat)
          zts%jd(:) = tvar%dat%vd(:) * z_tuf + jd0
          CALL INIT_NCVAR(tvar)
          !
          EXIT  ! ONLY ONE UNLIMITED DIMENSION POSSIBLE

       END IF unlimited_id
    END DO dimension_loop

    ! STEP 4: CONSTRUCT PARAMETER AXIS ---------------------------------------
    ! set i to non-unlimited dimension (currently only one allowed !)
    DO i=1, var%ndims
       IF (var%dim(i)%fuid) THEN
          CYCLE
       ELSE
          EXIT
       END IF
    END DO

    IF (var%ndims == 1) THEN
       zts%np = 1  ! SCALAR VARIABLE, IF ONLY 1-D
       idim_par = 0
       zts%dimname = 'd0'
    ELSE
       zts%np = var%dim(i)%len
       idim_par = i
    ENDIF
    WRITE(*,*) '    ',substr,': NUMBER OF PARAMETERES: ',zts%np

    ! ALLOCATE MEMORY FOR PARAMETER AXIS
    ALLOCATE(zts%par(zts%np))
    ! ... pre-fill with index
    DO np=1, zts%np
       zts%par(np) = REAL(np, dp)
    END DO

    ! ... overwrite in case it is present
    IF (var%ndims > 1) THEN
       ! READ DATA
       IF (var%dim(i)%varid /= NULL_VARID) THEN
          CALL IMPORT_NCVAR( pvar, varid=var%dim(i)%varid &
               , file=TRIM(ADJUSTL(fname)) )
          CALL DOUBLE_NARRAY(pvar%dat)
          zts%par(:) = pvar%dat%vd(:)
          ! name of dimension (and dimension variable !)
          zts%dimname = TRIM(ADJUSTL(pvar%name))
          ! attributes of dimension variable
          IF (pvar%natts > 0) THEN
             ALLOCATE(zts%dimvaratt(pvar%natts))
             DO i=1, pvar%natts
                CALL COPY_NCATT(zts%dimvaratt(i), pvar%att(i))
             END DO
          END IF
          ! CLEAN
          CALL INIT_NCVAR(pvar)
       ELSE
          WRITE(*,*) substr &
               , ' WARNING: NO COORDINATE VARIABLE FOUND FOR DIMENSION ' &
               , TRIM(var%dim(i)%name)
       END IF
    END IF
    WRITE(*,*) '    ',substr,': AXIS                 : ',zts%par(:)

    ! STEP 5: READ DATA ------------------------------------------------------
    ! ALLOCATE MEMORY FOR DATA
    ALLOCATE(zts%data(zts%nt,zts%np))
    CALL DOUBLE_NARRAY(var%dat)
    ! attributes of variable
    IF (var%natts > 0) THEN
       ALLOCATE(zts%varatt(var%natts))
       DO i=1, var%natts
          CALL COPY_NCATT(zts%varatt(i), var%att(i))
       END DO
    END IF
    DO i=1, PRODUCT(var%dat%dim(1:var%dat%n))
       CALL element(var%dat%dim, i, vec)
       nt = vec(idim_time)
       IF (idim_par == 0) THEN
          np = 1
       ELSE
          np = vec(idim_par)
       END IF
       zts%data(nt,np) = var%dat%vd(i)
    END DO
    DEALLOCATE(vec); NULLIFY(vec)

    ! STEP 6 : CLEAN ---------------------------------------------------------
    CALL INIT_NCVAR(var)    
    WRITE(*,*) '    ',substr,': DONE WITH FILE ',TRIM(ADJUSTL(zts%io%fname))

    ! OK
    status = 0

  END SUBROUTINE its_read_ts_netcdf
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE its_delete_ts(zts)

    USE messy_main_grid_netcdf, ONLY: INIT_NCATT

    IMPLICIT NONE
    INTRINSIC :: HUGE

    ! I/O
    TYPE(T_TS), INTENT(INOUT) :: zts

    ! LOCAL
    INTEGER :: i

    IF (ASSOCIATED(zts%jd)) DEALLOCATE(zts%jd)
    NULLIFY(zts%jd)

    IF (ASSOCIATED(zts%par)) DEALLOCATE(zts%par)
    NULLIFY(zts%par)

    IF (ASSOCIATED(zts%data)) DEALLOCATE(zts%data)
    NULLIFY(zts%data)

    IF (zts%lalloc) DEALLOCATE(zts%obj)
    NULLIFY(zts%obj)

    IF (zts%lalloc) DEALLOCATE(zts%flg)
    NULLIFY(zts%flg)

    zts%nt = 0
    zts%np = 0

    zts%io%name = ''
    zts%io%fname = ''
    zts%io%vr(1) = -HUGE(0.0_dp)
    zts%io%vr(2) =  HUGE(0.0_dp)
    zts%io%cnt = (/TS_BD_STOP, TS_BD_STOP/)
    zts%io%im = TS_IM_PREV

    zts%dimname = ''
    IF (ASSOCIATED(zts%dimvaratt)) THEN
       DO i=1, SIZE(zts%dimvaratt)
          CALL INIT_NCATT(zts%dimvaratt(i))
       END DO
    ENDIF
    NULLIFY(zts%dimvaratt)

    IF (ASSOCIATED(zts%varatt)) THEN
       DO i=1, SIZE(zts%varatt)
          CALL INIT_NCATT(zts%varatt(i))
       END DO
    ENDIF
    NULLIFY(zts%varatt)

  END SUBROUTINE its_delete_ts
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE its_set_value_ts(status, zts, yr, mo, dy, hr, mi, se)

    USE messy_main_timer, ONLY: gregor2julian

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER,    INTENT(OUT)   :: status
    TYPE(T_TS), INTENT(INOUT) :: zts
    INTEGER,    INTENT(IN)    :: yr, mo, dy, hr, mi, se

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'its_set_value_ts'
    REAL(dp) :: jd
    INTEGER  :: nt, np
    REAL(dp) :: f
    INTEGER  :: i
    INTEGER, DIMENSION(6)  :: zpdt

    status = 1

    zpdt(1) = yr
    zpdt(2) = mo
    zpdt(3) = dy
    zpdt(4) = hr
    zpdt(5) = mi
    zpdt(6) = se

    DO i=1, 6
       IF (zts%io%pdt(i) >= 0) THEN
          zpdt(i) = zts%io%pdt(i)
       ENDIF
    END DO

    jd = gregor2julian(zpdt(1),zpdt(2),zpdt(3),zpdt(4),zpdt(5),zpdt(6)) + &
         zts%io%offset

    ! CHECK BOUNDS
    IF (jd < zts%jd(1)) THEN
       SELECT CASE(zts%io%cnt(1))
       CASE(TS_BD_STOP)
          ! ERROR
          WRITE(*,*) substr,' ERROR: time below lower bound of time series ',&
               TRIM(zts%io%name)
          RETURN
       CASE(TS_BD_CONT)
          zts%obj(:) = zts%data(1,:)
       END SELECT
    ENDIF

    IF (jd > zts%jd(zts%nt)) THEN
       SELECT CASE(zts%io%cnt(2))
       CASE(TS_BD_STOP)
          ! ERROR
          WRITE(*,*) substr,' ERROR: time above upper bound of time series ',&
               TRIM(zts%io%name)
          RETURN
       CASE(TS_BD_CONT)
          zts%obj(:) = zts%data(zts%nt,:)
       END SELECT
    ENDIF

    DO nt=1, zts%nt-1
       IF ( (jd >= zts%jd(nt)) .AND. (jd <= zts%jd(nt+1)) ) EXIT
    END DO

    ! op_pj_20120420+
    ! IF TIME IS OUT OF INTERVAL AND CONTINUATION IS REQUESTED
    ! (Note: in this case nt == zts%nt at end of loop above!),
    ! NO FURTHER PROCESSING (INTERPOLATION) CAN BE PERFORMED,
    ! BUT THE FLAG STILL NEEDS TO BE SET CORRECTLY.

    IF (nt == zts%nt) THEN
       DO np=1, zts%np
          IF ( (zts%obj(np) <= zts%io%vr(1)) .OR. &
               (zts%obj(np) >= zts%io%vr(2)) ) THEN
             ! OUT OF VALID RANGE
             zts%flg(np) = 0.0_dp
          ELSE
             zts%flg(np) = 1.0_dp
          END IF
       END DO
    ELSE
    ! op_pj_20120420-

       SELECT CASE(zts%io%im)
       CASE(TS_IM_PREV)
          !
          zts%obj(:) = zts%data(nt,:)
          DO np=1, zts%np
             IF ( (zts%obj(np) <= zts%io%vr(1)) .OR. &
                  (zts%obj(np) >= zts%io%vr(2)) ) THEN
                ! OUT OF VALID RANGE
                zts%flg(np) = 0.0_dp
             ELSE
                zts%flg(np) = 1.0_dp
             END IF
          END DO
          !
       CASE(TS_IM_NEXT)
          !
          zts%obj(:) = zts%data(nt+1,:)
          DO np=1, zts%np
             IF ( (zts%obj(np) <= zts%io%vr(1)) .OR. &
                  (zts%obj(np) >= zts%io%vr(2)) ) THEN
                ! OUT OF VALID RANGE
                zts%flg(np) = 0.0_dp
             ELSE
                zts%flg(np) = 1.0_dp
             END IF
          END DO
          !
       CASE(TS_IM_LINT)
          !
          f = (jd - zts%jd(nt)) / (zts%jd(nt+1) - zts%jd(nt))
          zts%obj(:) = f * zts%data(nt+1,:) + (1.0_dp -f) * zts%data(nt,:)      
          !
          DO np=1, zts%np
             IF ( (zts%data(nt+1,np) <= zts%io%vr(1)) .OR. &
                  (zts%data(nt+1,np) >= zts%io%vr(2)) .OR. &
                  (zts%data(nt  ,np) <= zts%io%vr(1)) .OR. &
                  (zts%data(nt  ,np) >= zts%io%vr(2)) ) THEN
                ! OUT OF VALID RANGE
                zts%flg(np) = 0.0_dp
             ELSE
                zts%flg(np) = 1.0_dp
             END IF
          END DO
          !
       END SELECT

    END IF ! op_j_20120420

    status = 0

  END SUBROUTINE its_set_value_ts
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE its_copy_io(d, s)

    IMPLICIT NONE

    ! I/O
    TYPE (T_TS_IO), INTENT(OUT) :: d
    TYPE (T_TS_IO), INTENT(IN)  :: s

    d%name   = s%name
    d%fname  = s%fname
    d%vr(:)  = s%vr(:)
    d%cnt(:) = s%cnt(:)
    d%im     = s%im
    d%pdt(:) = s%pdt(:)
    d%offset = s%offset

  END SUBROUTINE its_copy_io
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE import_ts_read_nml_ctrl(status, iou)

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'import_ts_read_nml_ctrl'
    NAMELIST /CTRL_TS/ TS

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    ! NOTE: already at definition

    CALL read_nml_open(lex, substr, iou, 'CTRL_TS', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_TS, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_TS', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0

  END SUBROUTINE import_ts_read_nml_ctrl
  ! --------------------------------------------------------------------------

! ============================================================================
END MODULE messy_main_import_ts
! ============================================================================
