! ============================================================================
MODULE messy_main_import_lt
! ============================================================================

  USE messy_main_constants_mem, ONLY: DP, STRLEN_MEDIUM, STRLEN_ULONG
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY
  USE messy_main_import,        ONLY: modstr

  IMPLICIT NONE
  PRIVATE
  SAVE
  INTRINSIC :: NULL

  PUBLIC :: DP

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodstr = 'import_lt'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodver = '0.2'

  INTEGER, PARAMETER :: STRLEN_OBJECT  = 2*STRLEN_MEDIUM + 1 + 4
  INTEGER, PARAMETER, PUBLIC :: NRANKMAX = 6

  TYPE T_LT_IO
     !
     ! NAME OF LOOKUP TABLE
     CHARACTER(LEN=STRLEN_OBJECT) :: tname = ''
     ! NAME OF FILE WITH DATA, POSSIBLY CONTAINS TEMPLATE(S)
     CHARACTER(LEN=STRLEN_ULONG)  :: fname = ''
     ! NAME OF VARIABLE
     CHARACTER(LEN=STRLEN_MEDIUM) :: vname = ''
     !
  END TYPE T_LT_IO
  PUBLIC :: T_LT_IO

  TYPE T_LT
     !
     ! IO
     TYPE(T_LT_IO) :: io
     !
     LOGICAL :: linit = .FALSE. ! LT was already initialised, to be set after
     !                          ! first read)
     !
     ! NAME OF FILE WITH DATA; TEMPLATES REPLACED
     CHARACTER(LEN=STRLEN_ULONG)  :: fname = ''
     LOGICAL :: ltpl    = .FALSE.  ! file name contains templates and therefore
     !                             ! LT requires updates
     LOGICAL :: lupdate = .TRUE.   ! update required in this time step
     !
     INTEGER :: rank = 0
     !
     ! RANK MUST BE NRANKMAX (see above)
     REAL(DP), DIMENSION(:,:,:,:,:,:), POINTER :: data => NULL()
     !
     INTEGER, DIMENSION(:), POINTER                      :: dimlen  => NULL()
     CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: dimname => NULL()
     CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: dimunit => NULL()
     !
     TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER           :: dimvar => NULL()
     !
  END TYPE T_LT
  PUBLIC :: T_LT

  ! WORKSPACE
  INTEGER,       PARAMETER,         PUBLIC :: NMAXLT = 100 ! max. numer
  INTEGER,                          PUBLIC :: NLT = 0      ! actual number
  TYPE(T_LT_IO), DIMENSION(NMAXLT), PUBLIC :: LT
  TYPE(T_LT),    DIMENSION(NMAXLT), PUBLIC :: XLT

  ! TO BE CALLED FROM SMIL OF SUBMODELS
  PUBLIC :: get_lookup_table

  ! INTERNAL USE
  PUBLIC :: import_lt_read_nml_ctrl
  PUBLIC :: ilt_read_lt_netcdf
  PUBLIC :: ilt_delete_lt
  PUBLIC :: ilt_set_file

CONTAINS

  ! --------------------------------------------------------------------------
  SUBROUTINE import_lt_read_nml_ctrl(status, iou)

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'import_lt_read_nml_ctrl'
    NAMELIST /CTRL_LT/ LT

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    ! NOTE: already at definition

    CALL read_nml_open(lex, substr, iou, 'CTRL_LT', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_LT, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_LT', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0

  END SUBROUTINE import_lt_read_nml_ctrl
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE ilt_read_lt_netcdf(status, id)

    USE messy_main_grid_netcdf, ONLY: t_ncvar, import_ncvar    &
                                    , INIT_NCVAR, NULL_VARID, string &
                                    , DOUBLE_NARRAY, element, VTYPE_CHAR

    IMPLICIT NONE
    INTRINSIC :: TRIM, ADJUSTL, PRODUCT

    ! I/O
    INTEGER,    INTENT(OUT)   :: status
    INTEGER,    INTENT(IN)    :: id

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ilt_read_lt_netcdf'
    !
    TYPE (t_ncvar) :: var, pvar !, tvar, pvar
    INTEGER        :: i, j, k
    INTEGER        :: nval
    INTEGER, DIMENSION(:), POINTER :: vec => NULL()

    status = 1 ! ERROR

    ! STEP 1: INFO -----------------------------------------------------------
    WRITE(*,*) '    ',substr,': TABLE ', &
         TRIM(ADJUSTL(XLT(id)%io%tname)),' ID = ',id
    IF (XLT(id)%linit) THEN
       WRITE(*,*) '    ',substr,': OPENING FILE ', &
            TRIM(ADJUSTL(XLT(id)%fname)), ' FOR UPDATING LT'
    ELSE
       WRITE(*,*) '    ',substr,': OPENING FILE ', &
            TRIM(ADJUSTL(XLT(id)%fname))
    END IF
    WRITE(*,*) '    ',substr,': READING VARIABLE ', &
         TRIM(ADJUSTL(XLT(id)%io%vname))

    ! STEP 2: IMPORT VARIABLE ------------------------------------------------
    CALL IMPORT_NCVAR( var, varname=TRIM(ADJUSTL(XLT(id)%io%vname)) &
         , file=TRIM(ADJUSTL(XLT(id)%fname)) )

    ! STEP 3: CHECK RANK
    IF (var%ndims < 1) THEN
       WRITE(*,*) substr &
            , ' ERROR: ',TRIM(ADJUSTL(XLT(id)%io%vname)) &
            , ' HAS LESS THAN ',1,' DIMENSION'
       CALL INIT_NCVAR(var)
       RETURN
    END IF
    IF (var%ndims > NRANKMAX) THEN
       WRITE(*,*) substr &
            , ' ERROR: ',TRIM(ADJUSTL(XLT(id)%io%vname)) &
            , ' HAS MORE THAN ',NRANKMAX,' DIMENSIONS '
       CALL INIT_NCVAR(var)
       RETURN
    END IF

    ! STEP 4: EXTRACT DIMENSION INFORMATION ---------------------------------
    IF (XLT(id)%linit) THEN
       IF (XLT(id)%rank /= var%ndims) THEN
          WRITE(*,*) '    ',substr,': ERROR: NUMBER OF DIMENSIONS CHANGED '//&
               &'AFTER UPDATE: ',XLT(id)%rank,' -> ',var%ndims
          RETURN
       END IF
    ELSE
       XLT(id)%rank = var%ndims
    END IF
    WRITE(*,*) '    ',substr,': ',XLT(id)%rank,' DIMENSIONS'

    IF (.NOT. XLT(id)%linit) THEN
       ALLOCATE(XLT(id)%dimvar(XLT(id)%rank))
    END IF

    dimension_loop: DO i=1, var%ndims

       ! DIMENSION LENGTH
       IF (XLT(id)%linit) THEN
          IF (XLT(id)%dimlen(i) /= var%dim(i)%len) THEN
             WRITE(*,*) '    ',substr,': ERROR: LENGTH OF DIMENSION ',i&
                  ,' CHANGED '//&
                  &'AFTER UPDATE: ',XLT(id)%dimlen(i),' -> ',var%dim(i)%len
             RETURN
          END IF
       ELSE
          XLT(id)%dimlen(i)  = var%dim(i)%len
       ENDIF

       ! DIMENSION NAME
       IF (XLT(id)%linit) THEN
          IF (TRIM(XLT(id)%dimname(i)) /= TRIM(ADJUSTL(var%dim(i)%name))) THEN
             WRITE(*,*) '    ',substr,': ERROR: NAME OF DIMENSION CHANGED '//&
                  &'AFTER UPDATE: ',TRIM(XLT(id)%dimname(i)),' -> '&
                  ,TRIM(var%dim(i)%name)
             RETURN
          ENDIF
       ELSE
          XLT(id)%dimname(i) = TRIM(ADJUSTL(var%dim(i)%name))
       END IF

       WRITE(*,*) '    ',substr,': DIMENSION ',i,' NAMED '&
            , TRIM(ADJUSTL(var%dim(i)%name)),' HAS LENGTH ',XLT(id)%dimlen(i)

       ! DIMENSION AXES
       IF (.NOT. XLT(id)%linit) THEN
          ALLOCATE(XLT(id)%dimvar(i)%ptr(XLT(id)%dimlen(i)))
       END IF

       IF (var%dim(i)%varid == NULL_VARID) THEN
          ! NO DIMENSION VARIABLE PRESENT: USE INDEX
          WRITE(*,*) '    ',substr,': -> NO DIMENSION VARIABLE PRESENT'
          DO j=1, XLT(id)%dimlen(i)
             XLT(id)%dimvar(i)%ptr(j) = REAL(j,DP)
          END DO
          XLT(id)%dimunit(i) = 'index'
       ELSE
          WRITE(*,*) '    ',substr,': -> READING DIMENSION VARIABLE'
          CALL IMPORT_NCVAR( pvar, varid=var%dim(i)%varid &
               , file=TRIM(ADJUSTL(XLT(id)%fname)) )
          CALL DOUBLE_NARRAY(pvar%dat)
          XLT(id)%dimvar(i)%ptr(:) = pvar%dat%vd(:)
          ! UNIT
          XLT(id)%dimunit(i) = 'unknown'
          DO k=1, pvar%natts
             IF ( (pvar%att(k)%dat%type == VTYPE_CHAR) .AND. &
                  (TRIM(ADJUSTL(pvar%att(k)%name)) == 'units') ) THEN
                XLT(id)%dimunit(i) = TRIM(string(pvar%att(k)%dat%vc))
             END IF
          END DO
          ! CLEAN
          CALL INIT_NCVAR(pvar)
       ENDIF
       WRITE(*,*) '    ',substr,': -> units are ',TRIM(XLT(id)%dimunit(i))
       WRITE(*,*) '    ',substr,': -> values are ',XLT(id)%dimvar(i)%ptr(:)
    END DO dimension_loop

    ! STEP 5: EXTRACT DATA ---------------------------------------------------
    IF (.NOT. XLT(id)%linit) THEN
       ! ALLOCATE MEMORY
       ALLOCATE( XLT(id)%data( &
            XLT(id)%dimlen(1), &
            XLT(id)%dimlen(2), &
            XLT(id)%dimlen(3), &
            XLT(id)%dimlen(4), &
            XLT(id)%dimlen(5), &
            XLT(id)%dimlen(6)  ) )
    END IF

    CALL DOUBLE_NARRAY(var%dat)

    ! should be the same as: PRODUCT(var%dat%dim(1:var%dat%n))
    nval = PRODUCT(XLT(id)%dimlen(:))

    DO i=1, nval
       CALL element(XLT(id)%dimlen(:), i, vec)
       XLT(id)%data( vec(1), vec(2), vec(3), vec(4), vec(5), vec(6) ) = &
            var%dat%vd(i)
    END DO
    DEALLOCATE(vec); NULLIFY(vec)

    ! STEP 6: CLEAN ---------------------------------------------------------
    CALL INIT_NCVAR(var)
    WRITE(*,*) '    ',substr,': DONE WITH TABLE ',&
         TRIM(ADJUSTL(XLT(id)%io%tname))

    ! STEP 7:
    ! LT HAS BEEN INITIALISED
    XLT(id)%linit = .TRUE.

    ! OK
    status = 0

  END SUBROUTINE ilt_read_lt_netcdf
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE ilt_delete_lt(id)

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: id

    ! LOCAL
    INTEGER :: i

    IF (ASSOCIATED(XLT(id)%data)) THEN
       DEALLOCATE(XLT(id)%data) ; NULLIFY(XLT(id)%data)
    END IF

    IF (ASSOCIATED(XLT(id)%dimvar)) THEN
       dimension_loop: DO i=1, XLT(id)%rank
          IF (ASSOCIATED(XLT(id)%dimlen)) THEN
             DEALLOCATE(XLT(id)%dimlen) ; NULLIFY(XLT(id)%dimlen)
          ENDIF
          IF (ASSOCIATED(XLT(id)%dimname)) THEN
             DEALLOCATE(XLT(id)%dimname) ; NULLIFY(XLT(id)%dimname)
          ENDIF
          IF (ASSOCIATED(XLT(id)%dimunit)) THEN
             DEALLOCATE(XLT(id)%dimunit) ; NULLIFY(XLT(id)%dimunit)
          ENDIF
          IF (ASSOCIATED(XLT(id)%dimvar(i)%ptr)) THEN
             DEALLOCATE(XLT(id)%dimvar(i)%ptr); NULLIFY(XLT(id)%dimvar(i)%ptr)
          END IF
       END DO dimension_loop
       DEALLOCATE(XLT(id)%dimvar); NULLIFY(XLT(id)%dimvar)
    END IF

  END SUBROUTINE ilt_delete_lt
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE ilt_set_file(status,id,yr,mo,dy,hr,mi,se,linit)

    USE messy_main_tools, ONLY: int2str

    IMPLICIT NONE
    INTRINSIC :: INDEX, TRIM, ADJUSTL

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: id
    INTEGER, INTENT(IN)  :: yr,mo,dy,hr,mi,se
    LOGICAL, INTENT(IN)  :: linit

    ! LOCAL
    CHARACTER, PARAMETER        :: t = '%'
    CHARACTER(LEN=STRLEN_ULONG) :: fname = ''
    CHARACTER                   :: ch
    CHARACTER(LEN=4)            :: s4
    CHARACTER(LEN=2)            :: s2
    INTEGER                     :: i, j

    IF (linit) THEN
       ! check, if %io%fname contains templates
       XLT(id)%ltpl = (INDEX(TRIM(XLT(id)%io%fname), t) > 0)
    END IF

    IF (XLT(id)%ltpl) THEN
       fname = ''
       i = 0
       j = 0
       DO
          i = i + 1
          IF (i > LEN_TRIM(XLT(id)%io%fname)) EXIT
          IF (XLT(id)%io%fname(i:i) == t) THEN
             i = i + 1
             ch = XLT(id)%io%fname(i:i)
             s4 = ''
             s2 = ''
             SELECT CASE(ch)
             CASE('Y')
                CALL int2str(s4, yr, '0')
             CASE('M')
                CALL int2str(s2, mo, '0')
             CASE('D')
                CALL int2str(s2, dy, '0')
             CASE('h')
                CALL int2str(s2, hr, '0')
             CASE('m')
                CALL int2str(s2, mi, '0')
             CASE('s')
                CALL int2str(s2, se, '0')
             CASE DEFAULT
                ! unknown template
                status = 1
                RETURN
             END SELECT
             j = j + LEN_TRIM(ADJUSTL(s4)) + LEN_TRIM(ADJUSTL(s2))
             fname = TRIM(fname)//TRIM(ADJUSTL(s4))//TRIM(ADJUSTL(s2))
          ELSE
             j = j + 1
             fname(j:j) = XLT(id)%io%fname(i:i)
          END IF
       END DO
       !
       IF (linit) THEN
          XLT(id)%fname = TRIM(fname)
          XLT(id)%lupdate = .FALSE. ! no update directly after initialisation
       ELSE
          ! update required, if file name changes after templates have been
          ! replaced ...
          XLT(id)%lupdate = TRIM(XLT(id)%fname) /= TRIM(fname)
          ! ... with new filename
          IF (XLT(id)%lupdate) XLT(id)%fname = TRIM(fname)
       END IF
       !
    ELSE
       ! no templates: filename is that from namelist without replacements
       XLT(id)%fname = TRIM(XLT(id)%io%fname)
       XLT(id)%lupdate = .FALSE. ! update never required
    END IF

    status = 0

  END SUBROUTINE ilt_set_file
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE get_lookup_table(status, tname, p1, p2, p3, p4, p5, p6, id, &
       rank, dimlen, dimaxe, dimname, dimunit)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT)             :: status
    CHARACTER(LEN=*), INTENT(IN)              :: tname
    REAL(DP), DIMENSION(:),           POINTER, OPTIONAL :: p1  ! INTENT(OUT)
    REAL(DP), DIMENSION(:,:),         POINTER, OPTIONAL :: p2
    REAL(DP), DIMENSION(:,:,:),       POINTER, OPTIONAL :: p3
    REAL(DP), DIMENSION(:,:,:,:),     POINTER, OPTIONAL :: p4
    REAL(DP), DIMENSION(:,:,:,:,:),   POINTER, OPTIONAL :: p5
    REAL(DP), DIMENSION(:,:,:,:,:,:), POINTER, OPTIONAL :: p6
    INTEGER,           INTENT(OUT),            OPTIONAL :: id
    INTEGER,           INTENT(OUT),            OPTIONAL :: rank
    INTEGER,                      DIMENSION(:), POINTER, OPTIONAL :: dimlen
    TYPE(PTR_1D_ARRAY),           DIMENSION(:), POINTER, OPTIONAL :: dimaxe
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER, OPTIONAL :: dimname
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER, OPTIONAL :: dimunit

    ! LOCAL
    INTEGER :: i, zrank

    status = 1

    DO i=1, NLT
       IF ( TRIM(ADJUSTL(tname)) == TRIM(ADJUSTL(XLT(i)%io%tname)) ) THEN
          IF (PRESENT(p1)) p1 => XLT(i)%data(:,1,1,1,1,1)
          IF (PRESENT(p2)) p2 => XLT(i)%data(:,:,1,1,1,1)
          IF (PRESENT(p3)) p3 => XLT(i)%data(:,:,:,1,1,1)
          IF (PRESENT(p4)) p4 => XLT(i)%data(:,:,:,:,1,1)
          IF (PRESENT(p5)) p5 => XLT(i)%data(:,:,:,:,:,1)
          IF (PRESENT(p6)) p6 => XLT(i)%data(:,:,:,:,:,:)
          IF (PRESENT(id)) id = i
          IF (PRESENT(rank)) rank = XLT(i)%rank
          zrank = XLT(i)%rank
          IF (PRESENT(dimlen))  dimlen  => XLT(i)%dimlen(1:zrank)
          IF (PRESENT(dimname)) dimname => XLT(i)%dimname(1:zrank)
          IF (PRESENT(dimunit)) dimunit => XLT(i)%dimunit(1:zrank)
          IF (PRESENT(dimaxe))  dimaxe  => XLT(i)%dimvar
          ! FOUND: RETURN
          status = 0
          RETURN
       ENDIF
    END DO

    ! NOT FOUND: RESET
    IF (PRESENT(p1)) NULLIFY(p1)
    IF (PRESENT(p2)) NULLIFY(p2)
    IF (PRESENT(p3)) NULLIFY(p3)
    IF (PRESENT(p4)) NULLIFY(p4)
    IF (PRESENT(p5)) NULLIFY(p5)
    IF (PRESENT(p6)) NULLIFY(p6)
    IF (PRESENT(id)) id = 0
    IF (PRESENT(rank)) rank = 0
    IF (PRESENT(dimlen)) NULLIFY(dimlen)
    IF (PRESENT(dimname)) NULLIFY(dimname)
    IF (PRESENT(dimunit)) NULLIFY(dimunit)
    IF (PRESENT(dimaxe))  NULLIFY(dimaxe)

  END SUBROUTINE get_lookup_table
  ! --------------------------------------------------------------------------

! ============================================================================
END MODULE messy_main_import_lt
! ============================================================================
